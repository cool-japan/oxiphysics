//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

// TODO(P4): public API churn
use super::functions::*;
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
/// A SLAM particle holding a pose hypothesis and a local landmark map.
#[derive(Debug, Clone)]
pub struct SlamParticle {
    /// Hypothesis pose: (x, y, yaw).
    pub pose: (f64, f64, f64),
    /// Weight.
    pub weight: f64,
    /// Known landmarks observed by this particle.
    pub landmarks: Vec<Landmark>,
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
    pub fn propagate(
        &mut self,
        delta_s: f64,
        delta_theta: f64,
        noise_s: f64,
        noise_theta: f64,
        rng: &mut impl Rng,
    ) {
        let ns: f64 = rng.random_range(-noise_s..noise_s);
        let nt: f64 = rng.random_range(-noise_theta..noise_theta);
        let (x, y, yaw) = self.pose;
        let yaw_new = yaw + delta_theta + nt;
        let x_new = x + (delta_s + ns) * yaw_new.cos();
        let y_new = y + (delta_s + ns) * yaw_new.sin();
        self.pose = (x_new, y_new, yaw_new);
    }
    /// Update weight based on range observations to known landmarks.
    ///
    /// `observations`: (landmark_id, measured_range).
    pub fn update_weight(&mut self, observations: &[(usize, f64)], range_sigma: f64) {
        for &(lid, measured_range) in observations {
            if let Some(lm) = self.landmarks.iter().find(|l| l.id == lid) {
                let (x, y, _) = self.pose;
                let pred_range = ((lm.position.x - x).powi(2) + (lm.position.y - y).powi(2)).sqrt();
                let err = measured_range - pred_range;
                let sigma2 = range_sigma * range_sigma;
                let likelihood = (-0.5 * err * err / sigma2).exp();
                self.weight *= likelihood;
            }
        }
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
/// Particle-filter SLAM.
#[derive(Debug, Clone)]
pub struct ParticleFilterSlam {
    /// Particle set.
    pub particles: Vec<SlamParticle>,
    /// Motion noise (std dev for translation).
    pub noise_s: f64,
    /// Motion noise (std dev for rotation).
    pub noise_theta: f64,
    /// Sensor range noise (m).
    pub range_sigma: f64,
}
impl ParticleFilterSlam {
    /// Create a filter with `n` particles initialised at origin.
    pub fn new(n: usize, noise_s: f64, noise_theta: f64, range_sigma: f64) -> Self {
        let particles = (0..n)
            .map(|_| SlamParticle::new(0.0, 0.0, 0.0, 1.0 / n as f64))
            .collect();
        Self {
            particles,
            noise_s,
            noise_theta,
            range_sigma,
        }
    }
    /// Propagate all particles with odometry.
    pub fn predict(&mut self, delta_s: f64, delta_theta: f64) {
        let mut rng = rand::rng();
        for p in &mut self.particles {
            p.propagate(
                delta_s,
                delta_theta,
                self.noise_s,
                self.noise_theta,
                &mut rng,
            );
        }
    }
    /// Update weights from range observations.
    pub fn update(&mut self, observations: &[(usize, f64)]) {
        for p in &mut self.particles {
            p.update_weight(observations, self.range_sigma);
        }
        self.normalise_weights();
    }
    /// Normalise particle weights.
    pub fn normalise_weights(&mut self) {
        let sum: f64 = self.particles.iter().map(|p| p.weight).sum();
        if sum > 1e-30 {
            for p in &mut self.particles {
                p.weight /= sum;
            }
        }
    }
    /// Systematic resampling.
    pub fn resample(&mut self) {
        let n = self.particles.len();
        if n == 0 {
            return;
        }
        let mut rng = rand::rng();
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
/// Unscented Kalman Filter state (simplified 3-D pose: \[x, y, yaw\]).
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
    /// Generate 2n+1 sigma points for n=3.
    pub fn sigma_points(&self) -> [[f64; 3]; 7] {
        let n = 3usize;
        let lambda = self.alpha * self.alpha * (n as f64 + self.kappa) - n as f64;
        let scale = ((n as f64 + lambda) * self.cov[0]).sqrt();
        let mut pts = [[0.0f64; 3]; 7];
        pts[0] = self.mean;
        for i in 0..n {
            let s = (self.cov[i * 3 + i] * (n as f64 + lambda)).sqrt();
            let mut p = self.mean;
            let mut m = self.mean;
            p[i] += s;
            m[i] -= s;
            pts[1 + i] = p;
            pts[1 + n + i] = m;
        }
        let _ = scale;
        pts
    }
    /// Predict with a constant-velocity model (delta_t).
    pub fn predict(&mut self, delta_t: f64, process_noise: f64) {
        self.mean[0] += 1.0 * delta_t;
        for i in 0..3 {
            self.cov[i * 3 + i] += process_noise;
        }
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
/// A landmark in the map.
#[derive(Debug, Clone, Copy)]
pub struct Landmark {
    /// World position (m).
    pub position: Vec2,
    /// Landmark ID.
    pub id: usize,
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
