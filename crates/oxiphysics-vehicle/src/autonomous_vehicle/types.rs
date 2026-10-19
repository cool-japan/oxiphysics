//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

// TODO(P4): public API churn
use super::functions::*;
use crate::{VehicleError, VehicleResult};
use rand::{Rng, RngExt};
use std::ops::{Add, Sub};

/// A detected radar target.
#[derive(Debug, Clone, Copy)]
pub struct RadarTarget {
    /// Target position in sensor frame (m).
    pub position: Vec2,
    /// Radial velocity (m/s, positive = receding).
    pub radial_velocity: f64,
    /// Radar cross-section estimate (m²).
    pub rcs: f64,
    /// Signal-to-noise ratio (dB).
    pub snr_db: f64,
}
/// A forward safety corridor computed via forward reachability.
#[derive(Debug, Clone)]
pub struct SafetyEnvelope {
    /// Sequence of safe position intervals (x_min, x_max, y_min, y_max) per time step.
    pub intervals: Vec<(f64, f64, f64, f64)>,
    /// Time horizon (s).
    pub horizon: f64,
    /// Time step (s).
    pub dt: f64,
}
impl SafetyEnvelope {
    /// Compute safety envelope for a bicycle model vehicle.
    ///
    /// Assumes the vehicle can deviate at most `lateral_accel` (m/s²) laterally
    /// and `decel_max` (m/s²) in deceleration.
    pub fn compute(
        initial: BicycleState,
        horizon_s: f64,
        dt: f64,
        lateral_accel: f64,
        decel_max: f64,
    ) -> Self {
        let steps = (horizon_s / dt).ceil() as usize;
        let mut intervals = Vec::with_capacity(steps);
        let mut x_min = initial.x;
        let mut x_max = initial.x;
        let mut y_min = initial.y;
        let mut y_max = initial.y;
        let mut v_min = initial.v;
        let v_max = initial.v;
        for _k in 0..steps {
            let x_fwd = x_max + v_max * dt;
            let x_bwd = x_min + v_min.max(0.0) * dt;
            x_max = x_fwd;
            x_min = x_bwd;
            let lat_spread = 0.5 * lateral_accel * dt * dt;
            y_min -= lat_spread;
            y_max += lat_spread;
            v_min = (v_min - decel_max * dt).max(0.0);
            intervals.push((x_min, x_max, y_min, y_max));
        }
        Self {
            intervals,
            horizon: horizon_s,
            dt,
        }
    }
    /// Check if a point (x, y) is inside the safety envelope at time step k.
    pub fn contains(&self, k: usize, x: f64, y: f64) -> bool {
        if k >= self.intervals.len() {
            return false;
        }
        let (xmin, xmax, ymin, ymax) = self.intervals[k];
        x >= xmin && x <= xmax && y >= ymin && y <= ymax
    }
}
/// A single LiDAR return point.
#[derive(Debug, Clone, Copy)]
pub struct LidarPoint {
    /// 3-D position in sensor frame (m).
    pub position: Vec3,
    /// Intensity (0–1).
    pub intensity: f64,
    /// Channel index.
    pub channel: usize,
}
/// Radar sensor configuration.
#[derive(Debug, Clone)]
pub struct RadarConfig {
    /// Maximum detection range (m).
    pub range_max: f64,
    /// Horizontal field-of-view half-angle (rad).
    pub hfov_half: f64,
    /// Range resolution (m).
    pub range_resolution: f64,
    /// Velocity resolution (m/s).
    pub velocity_resolution: f64,
    /// Range noise std (m).
    pub range_noise_std: f64,
    /// Angle noise std (rad).
    pub angle_noise_std: f64,
    /// Minimum detectable RCS (m²).
    pub min_rcs: f64,
}
impl RadarConfig {
    /// Typical 77 GHz automotive radar.
    pub fn automotive_77ghz() -> Self {
        Self {
            range_max: 200.0,
            hfov_half: 10.0_f64.to_radians(),
            range_resolution: 0.3,
            velocity_resolution: 0.1,
            range_noise_std: 0.1,
            angle_noise_std: 0.5_f64.to_radians(),
            min_rcs: 1.0,
        }
    }
}
/// Simplified linear MPC horizon result.
#[derive(Debug, Clone)]
pub struct MpcResult {
    /// Planned states over the horizon.
    pub states: Vec<BicycleState>,
    /// Optimal control inputs (steer, accel) per step.
    pub controls: Vec<(f64, f64)>,
    /// Total cost.
    pub cost: f64,
}
/// A detected traffic sign.
#[derive(Debug, Clone)]
pub struct TrafficSign {
    /// Estimated position in world frame (m).
    pub position: Vec2,
    /// Sign type.
    pub sign_type: TrafficSignType,
    /// Detection confidence in \[0, 1\].
    pub confidence: f64,
    /// Estimated distance from ego (m).
    pub distance: f64,
}
/// Model Predictive Controller for trajectory tracking.
#[derive(Debug, Clone)]
pub struct Mpc {
    /// Prediction horizon (steps).
    pub horizon: usize,
    /// Time step (s).
    pub dt: f64,
    /// Wheelbase (m).
    pub wheelbase: f64,
    /// Weight on tracking error.
    pub q_xy: f64,
    /// Weight on speed error.
    pub q_v: f64,
    /// Weight on steering input.
    pub r_steer: f64,
    /// Weight on acceleration input.
    pub r_accel: f64,
    /// Maximum steering angle (rad).
    pub max_steer: f64,
    /// Maximum acceleration (m/s²).
    pub max_accel: f64,
}
impl Mpc {
    /// Default MPC parameters.
    pub fn default_params() -> Self {
        Self {
            horizon: 10,
            dt: 0.1,
            wheelbase: 2.8,
            q_xy: 1.0,
            q_v: 0.5,
            r_steer: 0.1,
            r_accel: 0.05,
            max_steer: 0.5,
            max_accel: 3.0,
        }
    }
    /// Solve via a simple shooting optimisation (gradient-free random sampling).
    ///
    /// Samples `n_samples` random control sequences and returns the lowest-cost one.
    pub fn solve(
        &self,
        initial_state: BicycleState,
        reference: &[BicycleState],
        n_samples: usize,
    ) -> MpcResult {
        let mut rng = rand::rng();
        let h = self.horizon.min(reference.len());
        let mut best_cost = f64::INFINITY;
        let mut best_controls: Vec<(f64, f64)> = vec![(0.0, 0.0); h];
        let mut best_states = vec![initial_state; h + 1];
        for _ in 0..n_samples {
            let mut states = Vec::with_capacity(h + 1);
            let mut controls = Vec::with_capacity(h);
            let mut state = initial_state;
            let mut cost = 0.0;
            states.push(state);
            for (_k, ref_k) in reference.iter().enumerate().take(h) {
                let steer: f64 = rng.random_range(-self.max_steer..self.max_steer);
                let accel: f64 = rng.random_range(-self.max_accel..self.max_accel);
                state = state.step(steer, accel, self.wheelbase, self.dt);
                let ref_k = *ref_k;
                let dx = state.x - ref_k.x;
                let dy = state.y - ref_k.y;
                let dv = state.v - ref_k.v;
                cost += self.q_xy * (dx * dx + dy * dy)
                    + self.q_v * dv * dv
                    + self.r_steer * steer * steer
                    + self.r_accel * accel * accel;
                controls.push((steer, accel));
                states.push(state);
            }
            if cost < best_cost {
                best_cost = cost;
                best_controls = controls;
                best_states = states;
            }
        }
        MpcResult {
            states: best_states,
            controls: best_controls,
            cost: best_cost,
        }
    }
}
/// Camera intrinsic parameters.
#[derive(Debug, Clone)]
pub struct CameraIntrinsics {
    /// Focal length x (pixels).
    pub fx: f64,
    /// Focal length y (pixels).
    pub fy: f64,
    /// Principal point x (pixels).
    pub cx: f64,
    /// Principal point y (pixels).
    pub cy: f64,
    /// Image width (pixels).
    pub width: u32,
    /// Image height (pixels).
    pub height: u32,
    /// Near clipping distance (m).
    pub near: f64,
    /// Far clipping distance (m).
    pub far: f64,
}
impl CameraIntrinsics {
    /// Project a 3-D point (camera frame) to pixel coordinates.
    /// Returns `None` if behind camera or outside image bounds.
    pub fn project(&self, p: Vec3) -> Option<(f64, f64)> {
        if p.z < self.near || p.z > self.far {
            return None;
        }
        let u = self.fx * p.x / p.z + self.cx;
        let v = self.fy * p.y / p.z + self.cy;
        if u >= 0.0 && u < self.width as f64 && v >= 0.0 && v < self.height as f64 {
            Some((u, v))
        } else {
            None
        }
    }
    /// Check whether a 3-D point is inside the camera frustum.
    pub fn in_frustum(&self, p: Vec3) -> bool {
        self.project(p).is_some()
    }
}
/// Types of traffic signs the recognition system can detect.
#[derive(Debug, Clone, PartialEq)]
pub enum TrafficSignType {
    /// Speed limit sign with value in km/h.
    SpeedLimit(u32),
    /// Stop sign.
    Stop,
    /// Yield sign.
    Yield,
    /// No-entry sign.
    NoEntry,
    /// Unknown / unrecognised sign.
    Unknown,
}
/// A FastSLAM particle: one pose hypothesis together with its own landmark
/// map (one 2-D Gaussian per landmark).
#[derive(Debug, Clone)]
pub struct SlamParticle {
    /// Hypothesis pose: (x, y, yaw).
    pub pose: (f64, f64, f64),
    /// Weight.
    pub weight: f64,
    /// Landmark map of this particle: a Gaussian position estimate per
    /// observed landmark, built by [`SlamParticle::observe`].
    pub landmarks: Vec<Landmark>,
}
/// Uniform noise sample in `[-half_width, half_width)`; zero for a
/// non-positive (or NaN) half width, for which the range would be empty.
fn symmetric_uniform(rng: &mut impl Rng, half_width: f64) -> f64 {
    if half_width > 0.0 {
        rng.random_range(-half_width..half_width)
    } else {
        0.0
    }
}
impl SlamParticle {
    /// Create a particle at a given pose with uniform weight.
    pub fn new(x: f64, y: f64, yaw: f64, weight: f64) -> Self {
        Self {
            pose: (x, y, yaw),
            weight,
            landmarks: Vec::new(),
        }
    }
    /// Propagate pose with a motion model: delta_s (forward), delta_theta (turn).
    ///
    /// The translation and rotation noise are uniform in
    /// `[-noise_s, noise_s)` and `[-noise_theta, noise_theta)`; a zero noise
    /// level gives deterministic odometry.
    pub fn propagate(
        &mut self,
        delta_s: f64,
        delta_theta: f64,
        noise_s: f64,
        noise_theta: f64,
        rng: &mut impl Rng,
    ) {
        let ns: f64 = symmetric_uniform(rng, noise_s);
        let nt: f64 = symmetric_uniform(rng, noise_theta);
        let (x, y, yaw) = self.pose;
        let yaw_new = yaw + delta_theta + nt;
        let x_new = x + (delta_s + ns) * yaw_new.cos();
        let y_new = y + (delta_s + ns) * yaw_new.sin();
        self.pose = (x_new, y_new, yaw_new);
    }
    /// Update the weight from one scan of range-bearing observations.
    ///
    /// Runs [`SlamParticle::observe`] (landmark map update) and multiplies the
    /// weight by the likelihood of the scan under this particle's pose and
    /// map.  `observations`: `(landmark_id, measured_range, measured_bearing)`.
    pub fn update_weight(
        &mut self,
        observations: &[(usize, f64, f64)],
        range_sigma: f64,
        bearing_sigma: f64,
    ) {
        let log_likelihood = self.observe(observations, range_sigma, bearing_sigma);
        self.weight *= log_likelihood.exp();
    }
    /// FastSLAM 1.0 measurement update of this particle.
    ///
    /// `observations` are `(landmark_id, range, bearing)` triples with the
    /// bearing measured from the particle's heading (rad), and the sensor
    /// noise is `diag(range_sigma², bearing_sigma²)`.  Observations are
    /// processed one after the other:
    ///
    /// * a landmark seen for the first time is initialised at
    ///   `pose + range * (cos, sin)(yaw + bearing)` with covariance
    ///   `G⁻¹ R G⁻ᵀ`, `G` being the Jacobian of the (range, bearing) model
    ///   with respect to the landmark position; it does not change the weight;
    /// * a known landmark is corrected with an EKF update (Joseph form), and
    ///   the log-density of the innovation `N(ν; 0, S)`, `S = G Σ Gᵀ + R`,
    ///   (bearing innovation wrapped to `[-π, π)`) is accumulated.
    ///
    /// Observations with non-finite or negative range, and landmarks exactly at
    /// the sensor position (degenerate bearing), are skipped.  Returns the
    /// total log-likelihood of the scan; the weight itself is not modified.
    pub fn observe(
        &mut self,
        observations: &[(usize, f64, f64)],
        range_sigma: f64,
        bearing_sigma: f64,
    ) -> f64 {
        let noise = [
            range_sigma * range_sigma,
            0.0,
            0.0,
            bearing_sigma * bearing_sigma,
        ];
        let mut log_likelihood = 0.0;
        for &(id, range, bearing) in observations {
            if !range.is_finite() || !bearing.is_finite() || range < 0.0 {
                continue;
            }
            match self.landmarks.iter().position(|l| l.id == id) {
                Some(index) => {
                    log_likelihood += self.correct_landmark(index, range, bearing, &noise);
                }
                None => {
                    let landmark = self.initial_landmark(id, range, bearing, &noise);
                    self.landmarks.push(landmark);
                }
            }
        }
        log_likelihood
    }
    /// Landmark initialised from its first range-bearing sighting.
    fn initial_landmark(&self, id: usize, range: f64, bearing: f64, noise: &[f64; 4]) -> Landmark {
        let (x, y, yaw) = self.pose;
        let (sin_phi, cos_phi) = (yaw + bearing).sin_cos();
        // Inverse of G = [[cos, sin], [-sin/r, cos/r]] at the measured point.
        let g_inv = [cos_phi, -range * sin_phi, sin_phi, range * cos_phi];
        let cov = mat2_mul(&mat2_mul(&g_inv, noise), &mat2_transpose(&g_inv));
        Landmark::new(id, Vec2::new(x + range * cos_phi, y + range * sin_phi), cov)
    }
    /// EKF correction of landmark `index`; returns the innovation
    /// log-likelihood (zero when the observation cannot be linearised).
    fn correct_landmark(
        &mut self,
        index: usize,
        range: f64,
        bearing: f64,
        noise: &[f64; 4],
    ) -> f64 {
        let (x, y, yaw) = self.pose;
        let Some(landmark) = self.landmarks.get_mut(index) else {
            return 0.0;
        };
        let dx = landmark.position.x - x;
        let dy = landmark.position.y - y;
        let q = dx * dx + dy * dy;
        let predicted_range = q.sqrt();
        if predicted_range < 1e-9 {
            return 0.0;
        }
        let g = [dx / predicted_range, dy / predicted_range, -dy / q, dx / q];
        let g_t = mat2_transpose(&g);
        let sigma = landmark.cov;
        let s = mat2_add(&mat2_mul(&mat2_mul(&g, &sigma), &g_t), noise);
        let Some((det_s, s_inv)) = mat2_det_inverse(&s) else {
            return 0.0;
        };
        if det_s <= 0.0 {
            return 0.0;
        }
        let nu = [
            range - predicted_range,
            wrap_angle(bearing - (dy.atan2(dx) - yaw)),
        ];
        // K = Sigma G^T S^-1
        let gain = mat2_mul(&mat2_mul(&sigma, &g_t), &s_inv);
        landmark.position.x += gain[0] * nu[0] + gain[1] * nu[1];
        landmark.position.y += gain[2] * nu[0] + gain[3] * nu[1];
        // Joseph form: Sigma' = (I - K G) Sigma (I - K G)^T + K R K^T.
        let kg = mat2_mul(&gain, &g);
        let i_kg = [1.0 - kg[0], -kg[1], -kg[2], 1.0 - kg[3]];
        let updated = mat2_add(
            &mat2_mul(&mat2_mul(&i_kg, &sigma), &mat2_transpose(&i_kg)),
            &mat2_mul(&mat2_mul(&gain, noise), &mat2_transpose(&gain)),
        );
        let off_diagonal = 0.5 * (updated[1] + updated[2]);
        landmark.cov = [updated[0], off_diagonal, off_diagonal, updated[3]];
        let mahalanobis = nu[0] * (s_inv[0] * nu[0] + s_inv[1] * nu[1])
            + nu[1] * (s_inv[2] * nu[0] + s_inv[3] * nu[1]);
        -std::f64::consts::TAU.ln() - 0.5 * det_s.ln() - 0.5 * mahalanobis
    }
}
/// 2-D point / vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
}
impl Vec2 {
    /// Construct from components.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    /// Zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0)
    }
    /// Euclidean length.
    #[inline]
    pub fn norm(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
    /// Squared length.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.x * self.x + self.y * self.y
    }
    /// Normalise.
    #[inline]
    pub fn normalise(self) -> Self {
        let n = self.norm();
        if n < 1e-300 {
            Self::zero()
        } else {
            Self::new(self.x / n, self.y / n)
        }
    }
    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }
    /// 2-D "cross" (scalar z-component).
    #[inline]
    pub fn cross2(self, rhs: Self) -> f64 {
        self.x * rhs.y - self.y * rhs.x
    }
    /// Scale.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
    /// Distance to another point.
    #[inline]
    pub fn dist(self, rhs: Self) -> f64 {
        (Self::new(self.x - rhs.x, self.y - rhs.y)).norm()
    }
}
impl std::ops::Add for Vec2 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}
impl std::ops::Sub for Vec2 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}
/// Default sensor bearing noise (rad) of a new [`ParticleFilterSlam`].
const DEFAULT_BEARING_SIGMA: f64 = 0.02;
/// Particle-filter SLAM (FastSLAM 1.0).
///
/// A Rao-Blackwellised particle filter: every particle carries a pose
/// hypothesis and, conditioned on that pose, an independent EKF per landmark
/// ([`SlamParticle::landmarks`]).  [`ParticleFilterSlam::predict`] samples the
/// motion model, [`ParticleFilterSlam::update`] corrects the particle maps with
/// range-bearing observations (known landmark identities) and re-weights the
/// particles, and [`ParticleFilterSlam::resample`] concentrates the particle
/// set on the likely pose hypotheses.
#[derive(Debug, Clone)]
pub struct ParticleFilterSlam {
    /// Particle set.
    pub particles: Vec<SlamParticle>,
    /// Motion noise: half width of the uniform translation noise (m).
    pub noise_s: f64,
    /// Motion noise: half width of the uniform rotation noise (rad).
    pub noise_theta: f64,
    /// Sensor range noise (m), standard deviation.
    pub range_sigma: f64,
    /// Sensor bearing noise (rad), standard deviation.
    pub bearing_sigma: f64,
}
impl ParticleFilterSlam {
    /// Create a filter with `n` particles initialised at origin.
    ///
    /// The bearing noise starts at 0.02 rad; change it with
    /// [`ParticleFilterSlam::with_bearing_sigma`].
    pub fn new(n: usize, noise_s: f64, noise_theta: f64, range_sigma: f64) -> Self {
        let particles = (0..n)
            .map(|_| SlamParticle::new(0.0, 0.0, 0.0, 1.0 / n as f64))
            .collect();
        Self {
            particles,
            noise_s,
            noise_theta,
            range_sigma,
            bearing_sigma: DEFAULT_BEARING_SIGMA,
        }
    }
    /// Set the sensor bearing noise (rad, standard deviation).
    pub fn with_bearing_sigma(mut self, bearing_sigma: f64) -> Self {
        self.bearing_sigma = bearing_sigma;
        self
    }
    /// Propagate all particles with odometry.
    pub fn predict(&mut self, delta_s: f64, delta_theta: f64) {
        let mut rng = rand::rng();
        self.predict_with_rng(delta_s, delta_theta, &mut rng);
    }
    /// [`ParticleFilterSlam::predict`] with a caller-supplied random source
    /// (for reproducible runs).
    pub fn predict_with_rng(&mut self, delta_s: f64, delta_theta: f64, rng: &mut impl Rng) {
        for p in &mut self.particles {
            p.propagate(delta_s, delta_theta, self.noise_s, self.noise_theta, rng);
        }
    }
    /// Update the particle maps and weights from range-bearing observations
    /// `(landmark_id, range, bearing)`.
    ///
    /// Every particle runs [`SlamParticle::observe`]; the resulting scan
    /// log-likelihoods are turned into weights relative to the best particle
    /// (log-sum-exp scaling), so that a poor scan cannot underflow all weights
    /// at once, and the weights are normalised.
    pub fn update(&mut self, observations: &[(usize, f64, f64)]) {
        let range_sigma = self.range_sigma;
        let bearing_sigma = self.bearing_sigma;
        let log_likelihoods: Vec<f64> = self
            .particles
            .iter_mut()
            .map(|p| p.observe(observations, range_sigma, bearing_sigma))
            .collect();
        let best = log_likelihoods
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        if best.is_finite() {
            for (p, ll) in self.particles.iter_mut().zip(log_likelihoods.iter()) {
                p.weight *= if ll.is_finite() {
                    (ll - best).exp()
                } else {
                    0.0
                };
            }
        }
        self.normalise_weights();
    }
    /// Normalise particle weights.
    ///
    /// If the total weight is not a finite number above `1e-300` (all
    /// hypotheses have been rejected, or a weight is NaN) the weights are reset
    /// to uniform instead of being divided by a vanishing sum.
    pub fn normalise_weights(&mut self) {
        let n = self.particles.len();
        let sum: f64 = self.particles.iter().map(|p| p.weight).sum();
        if sum.is_finite() && sum > 1e-300 {
            for p in &mut self.particles {
                p.weight /= sum;
            }
        } else {
            let uniform = 1.0 / n.max(1) as f64;
            for p in &mut self.particles {
                p.weight = uniform;
            }
        }
    }
    /// Systematic resampling.
    pub fn resample(&mut self) {
        let mut rng = rand::rng();
        self.resample_with_rng(&mut rng);
    }
    /// [`ParticleFilterSlam::resample`] with a caller-supplied random source.
    ///
    /// The weights are normalised first (a vanishing total resets them to
    /// uniform); afterwards every particle has weight `1 / n`.
    pub fn resample_with_rng(&mut self, rng: &mut impl Rng) {
        let n = self.particles.len();
        if n == 0 {
            return;
        }
        self.normalise_weights();
        let step = 1.0 / n as f64;
        let start: f64 = rng.random_range(0.0_f64..step);
        let mut cumsum = self.particles[0].weight;
        let mut idx = 0;
        let mut new_particles = Vec::with_capacity(n);
        let mut u = start;
        for _ in 0..n {
            while idx + 1 < n && u > cumsum {
                idx += 1;
                cumsum += self.particles[idx].weight;
            }
            new_particles.push(self.particles[idx].clone());
            u += step;
        }
        let w = 1.0 / n as f64;
        for p in &mut new_particles {
            p.weight = w;
        }
        self.particles = new_particles;
    }
    /// Effective sample size `1 / Σ w²` of the normalised weights: `n` for
    /// uniform weights down to `1` when one particle carries all the weight.
    /// Resampling is customarily triggered when it falls below `n / 2`.
    pub fn effective_sample_size(&self) -> f64 {
        let sum: f64 = self.particles.iter().map(|p| p.weight).sum();
        let sum_sq: f64 = self.particles.iter().map(|p| p.weight * p.weight).sum();
        if sum > 0.0 && sum_sq > 0.0 {
            sum * sum / sum_sq
        } else {
            self.particles.len() as f64
        }
    }
    /// The particle with the largest weight (the first one on ties).
    pub fn best_particle(&self) -> Option<&SlamParticle> {
        self.particles
            .iter()
            .fold(None, |best: Option<&SlamParticle>, p| match best {
                Some(b) if b.weight >= p.weight => Some(b),
                _ => Some(p),
            })
    }
    /// Weighted mean estimate of the world position of landmark `id` over all
    /// particles that have observed it, or `None` if no particle has.
    pub fn landmark_estimate(&self, id: usize) -> Option<Vec2> {
        let (mut sum_x, mut sum_y, mut sum_w) = (0.0, 0.0, 0.0);
        let (mut plain_x, mut plain_y, mut count) = (0.0, 0.0, 0.0);
        for p in &self.particles {
            if let Some(l) = p.landmarks.iter().find(|l| l.id == id) {
                sum_x += p.weight * l.position.x;
                sum_y += p.weight * l.position.y;
                sum_w += p.weight;
                plain_x += l.position.x;
                plain_y += l.position.y;
                count += 1.0;
            }
        }
        if sum_w > 0.0 {
            Some(Vec2::new(sum_x / sum_w, sum_y / sum_w))
        } else if count > 0.0 {
            Some(Vec2::new(plain_x / count, plain_y / count))
        } else {
            None
        }
    }
    /// Weighted mean pose estimate.
    pub fn mean_pose(&self) -> (f64, f64, f64) {
        let (mut wx, mut wy, mut ws, mut wc) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for p in &self.particles {
            wx += p.weight * p.pose.0;
            wy += p.weight * p.pose.1;
            ws += p.weight * p.pose.2.sin();
            wc += p.weight * p.pose.2.cos();
        }
        (wx, wy, ws.atan2(wc))
    }
}
/// A pedestrian agent.
#[derive(Debug, Clone)]
pub struct Pedestrian {
    /// Position (m).
    pub position: Vec2,
    /// Velocity (m/s).
    pub velocity: Vec2,
    /// Desired velocity (m/s).
    pub desired_velocity: Vec2,
    /// Mass (kg).
    pub mass: f64,
    /// Radius (m).
    pub radius: f64,
}
impl Pedestrian {
    /// Create a pedestrian.
    pub fn new(position: Vec2, desired_velocity: Vec2, mass: f64, radius: f64) -> Self {
        Self {
            position,
            velocity: Vec2::zero(),
            desired_velocity,
            mass,
            radius,
        }
    }
    /// Integrate social-force model step.
    pub fn step(
        &mut self,
        others: &[Pedestrian],
        obstacles: &[Vec2],
        dt: f64,
        tau: f64,
        a_rep: f64,
        b_rep: f64,
    ) {
        let f_drive = self
            .desired_velocity
            .sub(self.velocity)
            .scale(self.mass / tau);
        let mut f_rep_ped = Vec2::zero();
        for other in others {
            let d = self.position.dist(other.position);
            let sum_r = self.radius + other.radius;
            if d < 1e-12 {
                continue;
            }
            let dir = self.position.sub(other.position).normalise();
            let overlap = sum_r - d;
            let mag =
                a_rep * ((-d / b_rep).exp() + if overlap > 0.0 { 200.0 * overlap } else { 0.0 });
            f_rep_ped = f_rep_ped.add(dir.scale(mag));
        }
        let mut f_rep_obs = Vec2::zero();
        for &obs in obstacles {
            let d = self.position.dist(obs);
            if d < 1e-12 {
                continue;
            }
            let dir = self.position.sub(obs).normalise();
            let mag = a_rep * (-d / b_rep).exp();
            f_rep_obs = f_rep_obs.add(dir.scale(mag));
        }
        let f_total = f_drive.add(f_rep_ped).add(f_rep_obs);
        let accel = f_total.scale(1.0 / self.mass);
        self.velocity = self.velocity.add(accel.scale(dt));
        let speed = self.velocity.norm();
        let max_speed = self.desired_velocity.norm() * 2.0;
        if speed > max_speed && speed > 1e-12 {
            self.velocity = self.velocity.scale(max_speed / speed);
        }
        self.position = self.position.add(self.velocity.scale(dt));
    }
}
/// A polynomial lane model: y = a*x^2 + b*x + c in road-frame.
#[derive(Debug, Clone, Copy)]
pub struct LaneModel {
    /// Quadratic coefficient (1/m, curvature).
    pub a: f64,
    /// Linear coefficient (heading offset, rad).
    pub b: f64,
    /// Constant term (lateral offset, m).
    pub c: f64,
    /// Lane width (m).
    pub width: f64,
}
impl LaneModel {
    /// Evaluate lane centre lateral offset at longitudinal distance x.
    pub fn lateral_offset(&self, x: f64) -> f64 {
        self.a * x * x + self.b * x + self.c
    }
    /// Curvature of lane at x (1/m).
    pub fn curvature(&self, x: f64) -> f64 {
        2.0 * self.a + self.b * 0.0 * x
    }
    /// Fit a polynomial lane model to a set of (x, y) lane points.
    pub fn fit(points: &[(f64, f64)]) -> Option<Self> {
        if points.len() < 3 {
            return None;
        }
        let n = points.len() as f64;
        let (s00, mut s10, mut s20, mut s30, mut s40) = (n, 0.0, 0.0, 0.0, 0.0);
        let (mut r0, mut r1, mut r2) = (0.0, 0.0, 0.0);
        for &(x, y) in points {
            let x2 = x * x;
            let x3 = x2 * x;
            let x4 = x3 * x;
            s10 += x;
            s20 += x2;
            s30 += x3;
            s40 += x4;
            r0 += y;
            r1 += x * y;
            r2 += x2 * y;
        }
        let _ = s00;
        let a_mat = [n, s10, s20, s10, s20, s30, s20, s30, s40];
        let b_vec = [r0, r1, r2];
        solve3x3(&a_mat, &b_vec).map(|[c, b, a]| LaneModel {
            a,
            b,
            c,
            width: 3.5,
        })
    }
}
/// Artificial potential field planner.
#[derive(Debug, Clone)]
pub struct PotentialFieldPlanner {
    /// Attractive gain.
    pub k_att: f64,
    /// Repulsive gain.
    pub k_rep: f64,
    /// Influence radius of obstacles (m).
    pub d_safe: f64,
    /// Maximum step size (m).
    pub step_size: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}
impl PotentialFieldPlanner {
    /// Compute gradient descent path from `start` to `goal` avoiding obstacles.
    pub fn plan(&self, start: Vec2, goal: Vec2, obstacles: &[Vec2]) -> Vec<Vec2> {
        let mut path = vec![start];
        let mut pos = start;
        for _ in 0..self.max_iter {
            if pos.dist(goal) < self.step_size {
                path.push(goal);
                break;
            }
            let to_goal = goal.sub(pos).normalise();
            let f_att = to_goal.scale(self.k_att);
            let mut f_rep = Vec2::zero();
            for &obs in obstacles {
                let d = pos.dist(obs);
                if d < self.d_safe && d > 1e-12 {
                    let away = pos.sub(obs).normalise();
                    let mag = self.k_rep * (1.0 / d - 1.0 / self.d_safe) / (d * d);
                    f_rep = f_rep.add(away.scale(mag));
                }
            }
            let f_total = f_att.add(f_rep);
            let step = f_total.normalise().scale(self.step_size);
            pos = pos.add(step);
            path.push(pos);
        }
        path
    }
}
/// Unscented Kalman filter for a planar unicycle pose \[x, y, yaw\].
///
/// The filter uses the scaled unscented transform with the usual `alpha`,
/// `beta`, `kappa` parameters: `lambda = alpha² (n + kappa) − n` with `n = 3`
/// and `2n + 1 = 7` sigma points built from the columns of the lower Cholesky
/// factor of `(n + lambda) P`.  [`UkfState::predict`] propagates them through
/// the unicycle motion model and [`UkfState::update`] corrects the state with a
/// position measurement.  Yaw residuals are wrapped to `[-π, π)` everywhere.
#[derive(Debug, Clone)]
pub struct UkfState {
    /// Mean state vector \[x, y, yaw\].
    pub mean: [f64; 3],
    /// Covariance (3×3, row-major).
    pub cov: [f64; 9],
    /// UKF alpha parameter.
    pub alpha: f64,
    /// UKF beta parameter.
    pub beta: f64,
    /// UKF kappa parameter.
    pub kappa: f64,
}
/// State dimension of [`UkfState`].
const UKF_DIM: usize = 3;
impl UkfState {
    /// Create a new UKF state at origin with unit covariance.
    pub fn new() -> Self {
        let mut cov = [0.0f64; 9];
        for i in 0..3 {
            cov[i * 3 + i] = 1.0;
        }
        Self {
            mean: [0.0; 3],
            cov,
            alpha: 1e-3,
            beta: 2.0,
            kappa: 0.0,
        }
    }
    /// Spread factor `n + lambda = alpha² (n + kappa)` of the unscented
    /// transform; it must be a positive finite number.
    fn spread(&self) -> VehicleResult<f64> {
        let spread = self.alpha * self.alpha * (UKF_DIM as f64 + self.kappa);
        if spread.is_finite() && spread > 0.0 {
            Ok(spread)
        } else {
            Err(VehicleError::invalid_param(
                "alpha/kappa",
                "the sigma-point spread alpha² (n + kappa) must be positive and finite",
            ))
        }
    }
    /// Mean and covariance weights `(Wm, Wc)` of the 7 sigma points:
    /// `Wm0 = lambda / (n + lambda)`, `Wc0 = Wm0 + 1 − alpha² + beta` and
    /// `Wi = 1 / (2 (n + lambda))` for the other six.
    pub fn weights(&self) -> VehicleResult<([f64; 7], [f64; 7])> {
        let spread = self.spread()?;
        let lambda = spread - UKF_DIM as f64;
        let mut wm = [0.5 / spread; 7];
        let mut wc = [0.5 / spread; 7];
        wm[0] = lambda / spread;
        wc[0] = wm[0] + 1.0 - self.alpha * self.alpha + self.beta;
        Ok((wm, wc))
    }
    /// Generate the 2n+1 sigma points for n=3: the mean, followed by the mean
    /// plus and minus `sqrt(n + lambda)` times each column of the lower
    /// Cholesky factor of the covariance (so the off-diagonal terms of the
    /// covariance are honoured).
    ///
    /// The covariance is read through its lower triangle.  Fails if it is not
    /// positive semi-definite (or not finite) or if `alpha` / `kappa` do not
    /// give a positive spread.
    pub fn sigma_points(&self) -> VehicleResult<[[f64; 3]; 7]> {
        let scale = self.spread()?.sqrt();
        let factor = cholesky3_psd(&self.cov)
            .ok_or_else(|| VehicleError::general("UKF covariance is not positive semi-definite"))?;
        let mut pts = [self.mean; 7];
        for column in 0..UKF_DIM {
            for (row, mean_row) in self.mean.iter().enumerate() {
                let offset = scale * factor[row * UKF_DIM + column];
                pts[1 + column][row] = mean_row + offset;
                pts[1 + UKF_DIM + column][row] = mean_row - offset;
            }
        }
        Ok(pts)
    }
    /// Unscented prediction through the unicycle model with forward speed `v`
    /// (m/s) and yaw rate `omega` (rad/s) over `dt` seconds:
    ///
    /// `x += v cos(yaw) dt`, `y += v sin(yaw) dt`, `yaw += omega dt`.
    ///
    /// The sigma points are propagated individually, the yaw of the mean is
    /// computed from wrapped residuals, the covariance is
    /// `Σ Wc_i (χ_i − mean)(χ_i − mean)ᵀ + diag(q)` and `q` is the additive
    /// per-step process-noise variance of `[x, y, yaw]`.
    pub fn predict(&mut self, dt: f64, v: f64, omega: f64, q: [f64; 3]) -> VehicleResult<()> {
        if !(dt.is_finite() && dt >= 0.0) {
            return Err(VehicleError::invalid_param(
                "dt",
                "must be finite and non-negative",
            ));
        }
        if !(v.is_finite() && omega.is_finite()) {
            return Err(VehicleError::invalid_param("v/omega", "must be finite"));
        }
        if q.iter()
            .any(|variance| !(variance.is_finite() && *variance >= 0.0))
        {
            return Err(VehicleError::invalid_param(
                "q",
                "process-noise variances must be finite and non-negative",
            ));
        }
        let (wm, wc) = self.weights()?;
        let mut points = self.sigma_points()?;
        for point in &mut points {
            let yaw = point[2];
            point[0] += v * yaw.cos() * dt;
            point[1] += v * yaw.sin() * dt;
            point[2] = yaw + omega * dt;
        }
        // Yaw of the mean from wrapped residuals about the propagated centre
        // point: robust against the 2π jump and against the negative centre
        // weight of small-alpha parameterisations.
        let yaw_ref = points[0][2];
        let mut mean = [0.0f64; 3];
        for (point, w) in points.iter().zip(wm.iter()) {
            mean[0] += w * point[0];
            mean[1] += w * point[1];
            mean[2] += w * wrap_angle(point[2] - yaw_ref);
        }
        mean[2] = wrap_angle(yaw_ref + mean[2]);
        let mut cov = [0.0f64; 9];
        for (point, w) in points.iter().zip(wc.iter()) {
            let d = [
                point[0] - mean[0],
                point[1] - mean[1],
                wrap_angle(point[2] - mean[2]),
            ];
            add_scaled_outer3(&mut cov, *w, &d, &d);
        }
        for (i, variance) in q.iter().enumerate() {
            cov[i * UKF_DIM + i] += variance;
        }
        self.mean = mean;
        self.cov = Self::symmetrised(cov);
        Ok(())
    }
    /// Unscented measurement update with a position measurement `z = [x, y]`
    /// (m) of isotropic variance `r` (m²).
    ///
    /// The sigma points are mapped through the measurement function
    /// `h(state) = [x, y]`, giving the predicted measurement, the innovation
    /// covariance `S = Σ Wc_i dz_i dz_iᵀ + r I` and the state-measurement
    /// cross-covariance `Pxz`.  The state is corrected with the Kalman gain
    /// `K = Pxz S⁻¹` and the covariance becomes `P − K S Kᵀ`; the yaw
    /// residuals use wrapped angles.  Fails if `S` is singular.
    pub fn update(&mut self, z: [f64; 2], r: f64) -> VehicleResult<()> {
        if !(z[0].is_finite() && z[1].is_finite()) {
            return Err(VehicleError::invalid_param("z", "must be finite"));
        }
        if !(r.is_finite() && r >= 0.0) {
            return Err(VehicleError::invalid_param(
                "r",
                "measurement variance must be finite and non-negative",
            ));
        }
        let (wm, wc) = self.weights()?;
        let points = self.sigma_points()?;
        let mut z_pred = [0.0f64; 2];
        for (point, w) in points.iter().zip(wm.iter()) {
            z_pred[0] += w * point[0];
            z_pred[1] += w * point[1];
        }
        let mut s = [r, 0.0, 0.0, r];
        let mut p_xz = [0.0f64; 6];
        for (point, w) in points.iter().zip(wc.iter()) {
            let dz = [point[0] - z_pred[0], point[1] - z_pred[1]];
            let dx = [
                point[0] - self.mean[0],
                point[1] - self.mean[1],
                wrap_angle(point[2] - self.mean[2]),
            ];
            s[0] += w * dz[0] * dz[0];
            s[1] += w * dz[0] * dz[1];
            s[2] += w * dz[1] * dz[0];
            s[3] += w * dz[1] * dz[1];
            for (row, dx_i) in p_xz.chunks_exact_mut(2).zip(dx.iter()) {
                row[0] += w * dx_i * dz[0];
                row[1] += w * dx_i * dz[1];
            }
        }
        let (_, s_inv) = mat2_det_inverse(&s)
            .ok_or_else(|| VehicleError::general("UKF innovation covariance is singular"))?;
        // K = Pxz S^-1 (3x2).
        let mut gain = [0.0f64; 6];
        for (k_row, p_row) in gain.chunks_exact_mut(2).zip(p_xz.chunks_exact(2)) {
            k_row[0] = p_row[0] * s_inv[0] + p_row[1] * s_inv[2];
            k_row[1] = p_row[0] * s_inv[1] + p_row[1] * s_inv[3];
        }
        let innovation = [z[0] - z_pred[0], z[1] - z_pred[1]];
        for (m, k_row) in self.mean.iter_mut().zip(gain.chunks_exact(2)) {
            *m += k_row[0] * innovation[0] + k_row[1] * innovation[1];
        }
        self.mean[2] = wrap_angle(self.mean[2]);
        // P -= K S K^T
        for (i, k_i) in gain.chunks_exact(2).enumerate() {
            for (j, k_j) in gain.chunks_exact(2).enumerate() {
                let ks_j = [s[0] * k_j[0] + s[1] * k_j[1], s[2] * k_j[0] + s[3] * k_j[1]];
                self.cov[i * UKF_DIM + j] -= k_i[0] * ks_j[0] + k_i[1] * ks_j[1];
            }
        }
        let mut cov = Self::symmetrised(self.cov);
        for i in 0..UKF_DIM {
            // Roundoff can push a vanishing variance marginally below zero.
            cov[i * UKF_DIM + i] = cov[i * UKF_DIM + i].max(0.0);
        }
        self.cov = cov;
        Ok(())
    }
    /// `(A + Aᵀ) / 2` of a row-major 3x3 matrix.
    fn symmetrised(a: [f64; 9]) -> [f64; 9] {
        let mut out = a;
        for i in 0..UKF_DIM {
            for j in 0..UKF_DIM {
                out[i * UKF_DIM + j] = 0.5 * (a[i * UKF_DIM + j] + a[j * UKF_DIM + i]);
            }
        }
        out
    }
}
/// A 2-D probabilistic occupancy grid using log-odds representation.
#[derive(Debug, Clone)]
pub struct OccupancyGrid {
    /// Log-odds values (row-major, row = y, col = x).
    pub log_odds: Vec<f64>,
    /// Grid width in cells.
    pub width: usize,
    /// Grid height in cells.
    pub height: usize,
    /// Cell size (m).
    pub resolution: f64,
    /// Grid origin (world coords of cell \[0,0\] lower-left corner).
    pub origin: Vec2,
    /// Log-odds clamp maximum.
    pub l_max: f64,
    /// Log-odds clamp minimum.
    pub l_min: f64,
    /// Log-odds for occupied update.
    pub l_occ: f64,
    /// Log-odds for free update.
    pub l_free: f64,
}
impl OccupancyGrid {
    /// Create a new occupancy grid.
    pub fn new(width: usize, height: usize, resolution: f64, origin: Vec2) -> Self {
        Self {
            log_odds: vec![0.0; width * height],
            width,
            height,
            resolution,
            origin,
            l_max: 3.5,
            l_min: -3.5,
            l_occ: 0.85,
            l_free: -0.4,
        }
    }
    /// Convert world position to grid cell indices.  Returns `None` if outside.
    pub fn world_to_cell(&self, p: Vec2) -> Option<(usize, usize)> {
        let rel = p.sub(self.origin);
        let col = (rel.x / self.resolution).floor() as isize;
        let row = (rel.y / self.resolution).floor() as isize;
        if col < 0 || row < 0 || col as usize >= self.width || row as usize >= self.height {
            None
        } else {
            Some((col as usize, row as usize))
        }
    }
    /// Convert cell indices to world position (cell centre).
    pub fn cell_to_world(&self, col: usize, row: usize) -> Vec2 {
        Vec2::new(
            self.origin.x + (col as f64 + 0.5) * self.resolution,
            self.origin.y + (row as f64 + 0.5) * self.resolution,
        )
    }
    /// Linear index from (col, row).
    #[inline]
    fn idx(&self, col: usize, row: usize) -> usize {
        row * self.width + col
    }
    /// Mark a cell as occupied.
    pub fn update_occupied(&mut self, col: usize, row: usize) {
        let i = self.idx(col, row);
        self.log_odds[i] = (self.log_odds[i] + self.l_occ).min(self.l_max);
    }
    /// Mark a cell as free.
    pub fn update_free(&mut self, col: usize, row: usize) {
        let i = self.idx(col, row);
        self.log_odds[i] = (self.log_odds[i] + self.l_free).max(self.l_min);
    }
    /// Get the occupancy probability for a cell.
    pub fn probability(&self, col: usize, row: usize) -> f64 {
        let l = self.log_odds[self.idx(col, row)];
        1.0 / (1.0 + (-l).exp())
    }
    /// Bresenham ray-cast from sensor to hit: mark free cells, then occupied.
    pub fn update_ray(&mut self, sensor: Vec2, hit: Vec2) {
        if let Some((x0, y0)) = self.world_to_cell(sensor)
            && let Some((x1, y1)) = self.world_to_cell(hit)
        {
            let free_cells = bresenham_line(x0 as isize, y0 as isize, x1 as isize, y1 as isize);
            let n = free_cells.len();
            for (i, &(cx, cy)) in free_cells.iter().enumerate() {
                if cx < 0 || cy < 0 {
                    continue;
                }
                let (c, r) = (cx as usize, cy as usize);
                if c < self.width && r < self.height {
                    if i + 1 < n {
                        self.update_free(c, r);
                    } else {
                        self.update_occupied(c, r);
                    }
                }
            }
        }
    }
}
/// A* planner result.
#[derive(Debug, Clone)]
pub struct AStarResult {
    /// Sequence of (col, row) cell indices from start to goal.
    pub path: Vec<(usize, usize)>,
    /// Total path cost.
    pub cost: f64,
}
/// State of the kinematic bicycle model: \[x, y, yaw, speed\].
#[derive(Debug, Clone, Copy, Default)]
pub struct BicycleState {
    /// x position (m).
    pub x: f64,
    /// y position (m).
    pub y: f64,
    /// Heading angle (rad).
    pub yaw: f64,
    /// Speed (m/s).
    pub v: f64,
}
impl BicycleState {
    /// Integrate the bicycle model for one time step.
    ///
    /// Controls: steer_angle (rad), acceleration (m/s²).
    pub fn step(&self, steer: f64, accel: f64, wheelbase: f64, dt: f64) -> Self {
        let v_new = (self.v + accel * dt).max(0.0);
        let beta = steer.atan2(1.0) * 2.0;
        let yaw_rate = self.v * steer.tan() / wheelbase.max(1e-3);
        BicycleState {
            x: self.x + self.v * self.yaw.cos() * dt,
            y: self.y + self.v * self.yaw.sin() * dt,
            yaw: self.yaw + yaw_rate * dt + beta * 0.0,
            v: v_new,
        }
    }
    /// Euclidean distance to another state in position space.
    pub fn pos_dist(&self, other: &BicycleState) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
}
/// LiDAR sensor configuration.
#[derive(Debug, Clone)]
pub struct LidarConfig {
    /// Number of vertical beams (channels).
    pub num_channels: usize,
    /// Number of horizontal rays per rotation.
    pub rays_per_revolution: usize,
    /// Minimum range (m).
    pub range_min: f64,
    /// Maximum range (m).
    pub range_max: f64,
    /// Vertical field of view, half-angle (radians).
    pub vfov_half: f64,
    /// Angular noise standard deviation (radians).
    pub angle_noise_std: f64,
    /// Range noise standard deviation (m).
    pub range_noise_std: f64,
}
impl LidarConfig {
    /// Typical 64-channel spinning LiDAR (Velodyne HDL-64E class).
    pub fn hdl64() -> Self {
        Self {
            num_channels: 64,
            rays_per_revolution: 1800,
            range_min: 0.9,
            range_max: 100.0,
            vfov_half: 26.9_f64.to_radians(),
            angle_noise_std: 0.01_f64.to_radians(),
            range_noise_std: 0.02,
        }
    }
}
/// EKF state for a constant-velocity motion model: \[x, y, yaw, v, yaw_rate\].
#[derive(Debug, Clone)]
pub struct EkfState {
    /// State vector \[x, y, yaw, v, yaw_rate\].
    pub x: [f64; 5],
    /// Covariance matrix (5×5, row-major).
    pub p: [f64; 25],
}
impl EkfState {
    /// Initialise with zero state and small uncertainty.
    pub fn new() -> Self {
        let mut p = [0.0f64; 25];
        for i in 0..5 {
            p[i * 5 + i] = 1.0;
        }
        Self { x: [0.0; 5], p }
    }
    /// Predict step with process noise Q (5×5 row-major).
    pub fn predict(&mut self, dt: f64, q: &[f64; 25]) {
        let x = self.x;
        let v = x[3];
        let yaw = x[2];
        let yr = x[4];
        self.x[0] += v * yaw.cos() * dt;
        self.x[1] += v * yaw.sin() * dt;
        self.x[2] += yr * dt;
        let mut f = [0.0f64; 25];
        for i in 0..5 {
            f[i * 5 + i] = 1.0;
        }
        f[2] = -v * yaw.sin() * dt;
        f[3] = yaw.cos() * dt;
        f[5 + 2] = v * yaw.cos() * dt;
        f[5 + 3] = yaw.sin() * dt;
        f[2 * 5 + 4] = dt;
        let fp = mat5x5_mul(&f, &self.p);
        let ft = mat5x5_transpose(&f);
        let fpft = mat5x5_mul(&fp, &ft);
        for i in 0..25 {
            self.p[i] = fpft[i] + q[i];
        }
    }
    /// Update step with a 2-D position measurement \[z_x, z_y\], noise R (2×2).
    pub fn update_position(&mut self, z: [f64; 2], r: &[f64; 4]) {
        let mut h = [0.0f64; 10];
        h[0] = 1.0;
        h[5 + 1] = 1.0;
        let hpht = mat2x5_mul_5x5_mul_5x2(&h, &self.p);
        let s = [
            hpht[0] + r[0],
            hpht[1] + r[1],
            hpht[2] + r[2],
            hpht[3] + r[3],
        ];
        let ht = mat2x5_transpose(&h);
        let pht = mat5x5_mul_5x2(&self.p, &ht);
        let s_inv = mat2x2_inv(&s);
        let k = mat5x2_mul_2x2(&pht, &s_inv);
        let inn = [z[0] - self.x[0], z[1] - self.x[1]];
        for i in 0..5 {
            self.x[i] += k[i * 2] * inn[0] + k[i * 2 + 1] * inn[1];
        }
        let kh = mat5x2_mul_2x5(&k, &h);
        let mut ikh = [0.0f64; 25];
        for i in 0..5 {
            for j in 0..5 {
                ikh[i * 5 + j] = (if i == j { 1.0 } else { 0.0 }) - kh[i * 5 + j];
            }
        }
        self.p = mat5x5_mul(&ikh, &self.p);
    }
}
/// A landmark in the map: a Gaussian estimate of its world position.
#[derive(Debug, Clone, Copy)]
pub struct Landmark {
    /// Estimated world position (m).
    pub position: Vec2,
    /// Landmark ID.
    pub id: usize,
    /// Position covariance (m²), row-major 2x2 `[xx, xy, yx, yy]`.
    pub cov: [f64; 4],
}
impl Landmark {
    /// Landmark `id` at `position` with position covariance `cov`.
    pub fn new(id: usize, position: Vec2, cov: [f64; 4]) -> Self {
        Self { position, id, cov }
    }
    /// Landmark whose position is known exactly (zero covariance).
    pub fn known(id: usize, position: Vec2) -> Self {
        Self::new(id, position, [0.0; 4])
    }
}
/// A node in the RRT tree.
#[derive(Debug, Clone)]
pub struct RrtNode {
    /// Configuration (x, y) in metres.
    pub q: Vec2,
    /// Parent index in the tree (`usize::MAX` = root).
    pub parent: usize,
}
/// 3-D point / vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}
impl Vec3 {
    /// Construct from components.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    /// Zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
    /// Euclidean length.
    #[inline]
    pub fn norm(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
    /// Normalise.
    #[inline]
    pub fn normalise(self) -> Self {
        let n = self.norm();
        if n < 1e-300 {
            Self::zero()
        } else {
            Self::new(self.x / n, self.y / n, self.z / n)
        }
    }
    /// Scale.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
    /// Project to 2-D (drop z).
    #[inline]
    pub fn xy(self) -> Vec2 {
        Vec2::new(self.x, self.y)
    }
}
impl std::ops::Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl std::ops::Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
/// Rapidly-exploring Random Tree planner.
pub struct RrtPlanner {
    /// Tree nodes.
    pub nodes: Vec<RrtNode>,
    /// Goal tolerance (m).
    pub goal_tolerance: f64,
    /// Maximum extension step (m).
    pub step_size: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Configuration space bounds.
    pub bounds_min: Vec2,
    /// Configuration space bounds max.
    pub bounds_max: Vec2,
}
impl RrtPlanner {
    /// Create a new RRT planner.
    pub fn new(
        start: Vec2,
        step_size: f64,
        goal_tolerance: f64,
        max_iter: usize,
        bounds_min: Vec2,
        bounds_max: Vec2,
    ) -> Self {
        Self {
            nodes: vec![RrtNode {
                q: start,
                parent: usize::MAX,
            }],
            goal_tolerance,
            step_size,
            max_iter,
            bounds_min,
            bounds_max,
        }
    }
    /// Plan a path from the current start (root) to `goal`.
    ///
    /// `is_collision_free(a, b)` must return `true` if the segment a→b is free.
    pub fn plan(
        &mut self,
        goal: Vec2,
        is_collision_free: &dyn Fn(Vec2, Vec2) -> bool,
    ) -> Option<Vec<Vec2>> {
        let mut rng = rand::rng();
        for _iter in 0..self.max_iter {
            let q_rand = if rng.random_range(0.0_f64..1.0_f64) < 0.1 {
                goal
            } else {
                Vec2::new(
                    rng.random_range(self.bounds_min.x..self.bounds_max.x),
                    rng.random_range(self.bounds_min.y..self.bounds_max.y),
                )
            };
            let (nearest_idx, _) = self
                .nodes
                .iter()
                .enumerate()
                .map(|(i, n)| (i, n.q.dist(q_rand)))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .expect("operation should succeed");
            let q_near = self.nodes[nearest_idx].q;
            let dir = q_rand.sub(q_near).normalise();
            let q_new = q_near.add(dir.scale(self.step_size));
            if !is_collision_free(q_near, q_new) {
                continue;
            }
            let new_idx = self.nodes.len();
            self.nodes.push(RrtNode {
                q: q_new,
                parent: nearest_idx,
            });
            if q_new.dist(goal) <= self.goal_tolerance {
                let mut path = vec![goal];
                let mut idx = new_idx;
                while idx != usize::MAX {
                    path.push(self.nodes[idx].q);
                    let p = self.nodes[idx].parent;
                    if p == usize::MAX {
                        break;
                    }
                    idx = p;
                }
                path.reverse();
                return Some(path);
            }
        }
        None
    }
}
