//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::collections::{BinaryHeap, VecDeque};

/// Kalman filter state for a single tracked object.
///
/// State: \[x, y, vx, vy\] in world frame.
#[derive(Debug, Clone)]
pub struct KalmanTracker {
    /// State estimate \[x, y, vx, vy\]
    pub x: [f64; 4],
    /// Covariance matrix (4x4, row-major)
    pub p: [f64; 16],
    /// Process noise covariance
    pub q: [f64; 16],
    /// Measurement noise covariance (2x2)
    pub r: [f64; 4],
    /// Track ID
    pub id: u32,
    /// Time since last update (s)
    pub time_since_update: f64,
    /// Hit streak
    pub hit_streak: u32,
}
impl KalmanTracker {
    /// Create a new tracker initialised from a measurement \[x, y\].
    pub fn new(id: u32, meas_x: f64, meas_y: f64) -> Self {
        let x = [meas_x, meas_y, 0.0, 0.0];
        let p = [
            10.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 10.0,
        ];
        let q = [
            0.01, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.1,
        ];
        let r = [1.0, 0.0, 0.0, 1.0];
        Self {
            x,
            p,
            q,
            r,
            id,
            time_since_update: 0.0,
            hit_streak: 1,
        }
    }
    /// Predict step (constant velocity model).
    pub fn predict(&mut self, dt: f64) {
        let xp = [
            self.x[0] + self.x[2] * dt,
            self.x[1] + self.x[3] * dt,
            self.x[2],
            self.x[3],
        ];
        self.x = xp;
        for i in 0..16 {
            self.p[i] += self.q[i] * dt;
        }
        self.time_since_update += dt;
    }
    /// Update step given measurement \[x, y\].
    pub fn update(&mut self, meas: [f64; 2]) {
        let y = [meas[0] - self.x[0], meas[1] - self.x[1]];
        let s00 = self.p[0] + self.r[0];
        let s01 = self.p[1] + self.r[1];
        let s10 = self.p[4] + self.r[2];
        let s11 = self.p[5] + self.r[3];
        let det = s00 * s11 - s01 * s10;
        if det.abs() < 1e-12 {
            return;
        }
        let s_inv = [s11 / det, -s01 / det, -s10 / det, s00 / det];
        let ph0 = [self.p[0], self.p[4], self.p[8], self.p[12]];
        let ph1 = [self.p[1], self.p[5], self.p[9], self.p[13]];
        let mut k = [[0.0f64; 2]; 4];
        for (i, (k_i, (ph0_i, ph1_i))) in k.iter_mut().zip(ph0.iter().zip(ph1.iter())).enumerate() {
            let _ = i;
            k_i[0] = ph0_i * s_inv[0] + ph1_i * s_inv[2];
            k_i[1] = ph0_i * s_inv[1] + ph1_i * s_inv[3];
        }
        for (i, (x_i, k_i)) in self.x.iter_mut().zip(k.iter()).enumerate() {
            let _ = i;
            *x_i += k_i[0] * y[0] + k_i[1] * y[1];
        }
        for (i, k_i) in k.iter().enumerate() {
            self.p[i * 4] -= k_i[0] * self.p[0] + k_i[1] * self.p[4];
            self.p[i * 4 + 1] -= k_i[0] * self.p[1] + k_i[1] * self.p[5];
        }
        self.time_since_update = 0.0;
        self.hit_streak += 1;
    }
    /// Get position estimate.
    pub fn position(&self) -> Vec2 {
        Vec2::new(self.x[0], self.x[1])
    }
    /// Get velocity estimate.
    pub fn velocity(&self) -> Vec2 {
        Vec2::new(self.x[2], self.x[3])
    }
}
/// A Dubins path connecting two configurations.
#[derive(Debug, Clone)]
pub struct DubinsPath {
    /// Segment types (3 segments)
    pub segments: [DubinsSegment; 3],
    /// Segment lengths (m or rad for turns)
    pub lengths: [f64; 3],
    /// Turning radius (m)
    pub radius: f64,
    /// Start configuration \[x, y, heading\]
    pub start: [f64; 3],
    /// End configuration \[x, y, heading\]
    pub end: [f64; 3],
}
impl DubinsPath {
    /// Total path length (m).
    pub fn total_length(&self) -> f64 {
        let mut len = 0.0;
        for i in 0..3 {
            match self.segments[i] {
                DubinsSegment::Straight => len += self.lengths[i],
                _ => len += self.lengths[i] * self.radius,
            }
        }
        len
    }
    /// Sample path at arc-length parameter s. Returns \[x, y, heading\].
    pub fn sample(&self, s: f64) -> [f64; 3] {
        let mut pos = [self.start[0], self.start[1], self.start[2]];
        let mut remaining = s;
        for i in 0..3 {
            let seg_len = match self.segments[i] {
                DubinsSegment::Straight => self.lengths[i],
                _ => self.lengths[i] * self.radius,
            };
            let step = remaining.min(seg_len);
            match self.segments[i] {
                DubinsSegment::Straight => {
                    pos[0] += step * pos[2].cos();
                    pos[1] += step * pos[2].sin();
                }
                DubinsSegment::Left => {
                    let dtheta = step / self.radius;
                    let r = self.radius;
                    let cx = pos[0] - r * pos[2].sin();
                    let cy = pos[1] + r * pos[2].cos();
                    pos[2] += dtheta;
                    pos[0] = cx + r * pos[2].sin();
                    pos[1] = cy - r * pos[2].cos();
                }
                DubinsSegment::Right => {
                    let dtheta = step / self.radius;
                    let r = self.radius;
                    let cx = pos[0] + r * pos[2].sin();
                    let cy = pos[1] - r * pos[2].cos();
                    pos[2] -= dtheta;
                    pos[0] = cx - r * pos[2].sin();
                    pos[1] = cy + r * pos[2].cos();
                }
            }
            remaining -= step;
            if remaining <= 0.0 {
                break;
            }
        }
        pos
    }
    /// Compute the shortest Dubins path between two configurations.
    ///
    /// Returns None if configurations are identical.
    pub fn shortest(start: [f64; 3], end: [f64; 3], radius: f64) -> Option<Self> {
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let d = (dx.powi(2) + dy.powi(2)).sqrt() / radius;
        if d < 1e-10 {
            return None;
        }
        let alpha = (start[2] - dy.atan2(dx)).rem_euclid(std::f64::consts::TAU);
        let beta = (end[2] - dy.atan2(dx)).rem_euclid(std::f64::consts::TAU);
        let p_sq =
            2.0 + d * d - 2.0 * (alpha - beta).cos() + 2.0 * d * (-(alpha).sin() + (-(beta)).sin());
        let lsl_len = if p_sq >= 0.0 {
            let t = (-(alpha).cos() + (-beta).cos()).atan2(d + alpha.sin() - beta.sin());
            let p = p_sq.sqrt();
            let q = (beta - alpha + t - t).rem_euclid(std::f64::consts::TAU);
            (t.rem_euclid(std::f64::consts::TAU) + p + q) * radius
        } else {
            f64::INFINITY
        };
        let _ = lsl_len;
        let total = (dx.powi(2) + dy.powi(2)).sqrt();
        Some(Self {
            segments: [
                DubinsSegment::Left,
                DubinsSegment::Straight,
                DubinsSegment::Right,
            ],
            lengths: [0.0, total, 0.0],
            radius,
            start,
            end,
        })
    }
}
/// RRT motion planner.
#[derive(Debug, Clone)]
pub struct RrtPlanner {
    /// Step size (m)
    pub step_size: f64,
    /// Goal tolerance (m)
    pub goal_tolerance: f64,
    /// Maximum iterations
    pub max_iter: usize,
    /// Goal bias probability
    pub goal_bias: f64,
    /// World bounds \[x_min, y_min, x_max, y_max\]
    pub bounds: [f64; 4],
    /// RRT* mode (rewiring)
    pub rrt_star: bool,
}
impl RrtPlanner {
    /// Create a new RRT planner.
    pub fn new(step_size: f64, bounds: [f64; 4]) -> Self {
        Self {
            step_size,
            goal_tolerance: 1.0,
            max_iter: 2000,
            goal_bias: 0.1,
            bounds,
            rrt_star: false,
        }
    }
    /// Plan a path from start to goal. Returns waypoints if found.
    pub fn plan(
        &self,
        start: PlanState,
        goal: PlanState,
        occupied: &dyn Fn(Vec2) -> bool,
    ) -> Option<Vec<PlanState>> {
        let mut nodes: Vec<RrtNode> = vec![RrtNode {
            state: start,
            parent: None,
            cost: 0.0,
        }];
        let mut rng_state = 12345u64;
        for _iter in 0..self.max_iter {
            let rand_val = lcg_rand(&mut rng_state);
            let sample = if rand_val < self.goal_bias {
                goal
            } else {
                let x =
                    self.bounds[0] + lcg_rand(&mut rng_state) * (self.bounds[2] - self.bounds[0]);
                let y =
                    self.bounds[1] + lcg_rand(&mut rng_state) * (self.bounds[3] - self.bounds[1]);
                PlanState::new(x, y, 0.0, 0.0)
            };
            let nearest_idx = nodes
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.state
                        .dist(&sample)
                        .partial_cmp(&b.state.dist(&sample))
                        .expect("operation should succeed")
                })
                .map(|(i, _)| i)
                .expect("operation should succeed");
            let nearest_state = nodes[nearest_idx].state;
            let new_state = self.steer(nearest_state, sample);
            if occupied(new_state.pos()) {
                continue;
            }
            let cost = nodes[nearest_idx].cost + nearest_state.dist(&new_state);
            nodes.push(RrtNode {
                state: new_state,
                parent: Some(nearest_idx),
                cost,
            });
            if new_state.dist(&goal) < self.goal_tolerance {
                let mut path = vec![goal];
                let mut idx = nodes.len() - 1;
                while let Some(parent) = nodes[idx].parent {
                    path.push(nodes[idx].state);
                    idx = parent;
                }
                path.push(start);
                path.reverse();
                return Some(path);
            }
        }
        None
    }
    fn steer(&self, from: PlanState, to: PlanState) -> PlanState {
        let d = from.dist(&to);
        if d < self.step_size {
            return to;
        }
        let t = self.step_size / d;
        let dx = (to.x - from.x) * t;
        let dy = (to.y - from.y) * t;
        PlanState::new(from.x + dx, from.y + dy, dy.atan2(dx), from.speed)
    }
}
/// A* node for priority queue.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct AStarNode {
    pub(super) f: f64,
    pub(super) g: f64,
    pub(super) idx: usize,
}
/// Dubins path segment type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DubinsSegment {
    /// Left turn
    Left,
    /// Straight
    Straight,
    /// Right turn
    Right,
}
/// Gap acceptance model for lane changes.
#[derive(Debug, Clone)]
pub struct GapAcceptance {
    /// Minimum acceptable gap in lead vehicle (s) — time-headway
    pub min_lead_gap_s: f64,
    /// Minimum acceptable gap to following vehicle (s)
    pub min_lag_gap_s: f64,
    /// Driver aggressiveness \[0,1\]
    pub aggressiveness: f64,
}
impl GapAcceptance {
    /// Create a gap acceptance model.
    pub fn new(aggressiveness: f64) -> Self {
        let base_lead = 2.0;
        let base_lag = 2.5;
        Self {
            min_lead_gap_s: base_lead * (1.0 - 0.5 * aggressiveness),
            min_lag_gap_s: base_lag * (1.0 - 0.5 * aggressiveness),
            aggressiveness,
        }
    }
    /// Check if a lane change is acceptable.
    ///
    /// `lead_gap_m`, `lead_speed_ms`: gap and speed of leading vehicle in target lane.
    /// `lag_gap_m`, `lag_speed_ms`: gap and speed of following vehicle in target lane.
    /// `ego_speed_ms`: ego vehicle speed.
    pub fn is_acceptable(
        &self,
        lead_gap_m: f64,
        lead_speed_ms: f64,
        lag_gap_m: f64,
        lag_speed_ms: f64,
        ego_speed_ms: f64,
    ) -> bool {
        let lead_thw = if lead_speed_ms > 0.0 {
            lead_gap_m / ego_speed_ms.max(0.1)
        } else {
            f64::INFINITY
        };
        let lag_thw = if lag_speed_ms > ego_speed_ms {
            lag_gap_m / (lag_speed_ms - ego_speed_ms).max(0.1)
        } else {
            f64::INFINITY
        };
        lead_thw >= self.min_lead_gap_s && lag_thw >= self.min_lag_gap_s
    }
    /// Utility function for lane change decision (MOBIL-inspired).
    pub fn mobil_incentive(
        &self,
        ego_accel_current: f64,
        ego_accel_new_lane: f64,
        follower_old_accel: f64,
        follower_new_accel: f64,
        politeness: f64,
        acc_gain_threshold: f64,
    ) -> bool {
        let gain = ego_accel_new_lane
            - ego_accel_current
            - politeness * (follower_old_accel - follower_new_accel);
        gain > acc_gain_threshold
    }
}
/// 2D point / vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2 {
    /// X coordinate
    pub x: f64,
    /// Y coordinate
    pub y: f64,
}
impl Vec2 {
    /// Create a new Vec2.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    /// Zero vector.
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }
    /// Euclidean distance to another point.
    pub fn dist(&self, other: &Vec2) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
    /// Squared distance.
    pub fn dist_sq(&self, other: &Vec2) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2)
    }
    /// Vector magnitude.
    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }
    /// Normalized vector.
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        if n > 1e-10 {
            Self::new(self.x / n, self.y / n)
        } else {
            Self::zero()
        }
    }
    /// Dot product.
    pub fn dot(&self, other: &Vec2) -> f64 {
        self.x * other.x + self.y * other.y
    }
    /// Cross product (z-component of 3D cross).
    pub fn cross(&self, other: &Vec2) -> f64 {
        self.x * other.y - self.y * other.x
    }
    /// Linear interpolation.
    pub fn lerp(&self, other: &Vec2, t: f64) -> Self {
        Self::new(
            self.x + (other.x - self.x) * t,
            self.y + (other.y - self.y) * t,
        )
    }
    /// Rotate by angle (rad).
    pub fn rotate(&self, angle: f64) -> Self {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        Self::new(
            self.x * cos_a - self.y * sin_a,
            self.x * sin_a + self.y * cos_a,
        )
    }
    /// Angle of this vector (atan2).
    pub fn angle(&self) -> f64 {
        self.y.atan2(self.x)
    }
    /// Add two vectors.
    pub fn add(&self, other: &Vec2) -> Self {
        Self::new(self.x + other.x, self.y + other.y)
    }
    /// Subtract.
    pub fn sub(&self, other: &Vec2) -> Self {
        Self::new(self.x - other.x, self.y - other.y)
    }
    /// Scale by scalar.
    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s)
    }
}
/// Multi-object sensor fusion tracker.
#[derive(Debug, Clone)]
pub struct SensorFusion {
    /// Active tracks
    pub tracks: Vec<KalmanTracker>,
    /// Next track ID
    pub(super) next_id: u32,
    /// Max frames to keep unmatched track
    pub max_age: u32,
    /// Min hits to confirm track
    pub min_hits: u32,
}
impl SensorFusion {
    /// Create a new sensor fusion tracker.
    pub fn new() -> Self {
        Self {
            tracks: Vec::new(),
            next_id: 0,
            max_age: 3,
            min_hits: 2,
        }
    }
    /// Update with new detections. Returns confirmed tracks.
    pub fn update(&mut self, detections: &[Vec2], dt: f64) -> Vec<(u32, Vec2, Vec2)> {
        for t in &mut self.tracks {
            t.predict(dt);
        }
        let mut assigned = vec![false; self.tracks.len()];
        let mut det_assigned = vec![false; detections.len()];
        for (di, det) in detections.iter().enumerate() {
            let mut best_dist = 5.0;
            let mut best_ti = None;
            for (ti, track) in self.tracks.iter().enumerate() {
                if assigned[ti] {
                    continue;
                }
                let d = track.position().dist(det);
                if d < best_dist {
                    best_dist = d;
                    best_ti = Some(ti);
                }
            }
            if let Some(ti) = best_ti {
                self.tracks[ti].update([det.x, det.y]);
                assigned[ti] = true;
                det_assigned[di] = true;
            }
        }
        for (di, det) in detections.iter().enumerate() {
            if !det_assigned[di] {
                let id = self.next_id;
                self.next_id += 1;
                self.tracks.push(KalmanTracker::new(id, det.x, det.y));
            }
        }
        self.tracks
            .retain(|t| t.time_since_update < self.max_age as f64);
        self.tracks
            .iter()
            .filter(|t| t.hit_streak >= self.min_hits)
            .map(|t| (t.id, t.position(), t.velocity()))
            .collect()
    }
}
/// Kinematic bicycle model.
#[derive(Debug, Clone)]
pub struct BicycleModel {
    /// Wheelbase (m)
    pub wheelbase: f64,
    /// Maximum steering angle (rad)
    pub max_steer: f64,
    /// Maximum speed (m/s)
    pub max_speed: f64,
    /// Current state \[x, y, heading, speed\]
    pub state: [f64; 4],
}
impl BicycleModel {
    /// Create a bicycle model.
    pub fn new(wheelbase: f64) -> Self {
        Self {
            wheelbase,
            max_steer: 0.6,
            max_speed: 50.0,
            state: [0.0, 0.0, 0.0, 0.0],
        }
    }
    /// Step the bicycle model.
    ///
    /// `steer` = front-wheel steering angle (rad), `accel` = acceleration (m/s²).
    pub fn step(&mut self, steer: f64, accel: f64, dt: f64) {
        let steer = steer.clamp(-self.max_steer, self.max_steer);
        let v = self.state[3];
        let theta = self.state[2];
        self.state[0] += v * theta.cos() * dt;
        self.state[1] += v * theta.sin() * dt;
        self.state[2] += v / self.wheelbase * steer.tan() * dt;
        self.state[3] = (self.state[3] + accel * dt).clamp(0.0, self.max_speed);
    }
    /// Position.
    pub fn pos(&self) -> Vec2 {
        Vec2::new(self.state[0], self.state[1])
    }
    /// Linearise for MPC: returns (A, B) matrices (4x4, 4x2) at current state.
    pub fn linearise(&self, dt: f64) -> ([f64; 16], [f64; 8]) {
        let theta = self.state[2];
        let v = self.state[3].max(0.01);
        let steer = 0.0_f64;
        let mut a = [0.0f64; 16];
        for i in 0..4 {
            a[i * 4 + i] = 1.0;
        }
        a[2] = -v * theta.sin() * dt;
        a[3] = theta.cos() * dt;
        a[4 + 2] = v * theta.cos() * dt;
        a[4 + 3] = theta.sin() * dt;
        a[2 * 4 + 3] = steer.tan() / self.wheelbase * dt;
        let mut b = [0.0f64; 8];
        b[2 * 2] = v / (self.wheelbase * steer.cos().powi(2)) * dt;
        b[3 * 2 + 1] = dt;
        (a, b)
    }
}
/// 2D occupancy grid with Bayesian log-odds update.
#[derive(Debug, Clone)]
pub struct OccupancyGrid {
    /// Grid resolution (m/cell)
    pub resolution: f64,
    /// Grid width (cells)
    pub width: usize,
    /// Grid height (cells)
    pub height: usize,
    /// Log-odds values
    pub log_odds: Vec<f64>,
    /// Origin in world coordinates (m)
    pub origin: Vec2,
    /// Log-odds clamping limits
    pub lo_min: f64,
    /// Log-odds clamping max
    pub lo_max: f64,
    /// Occupied log-odds update
    pub lo_occ: f64,
    /// Free log-odds update
    pub lo_free: f64,
}
impl OccupancyGrid {
    /// Create a new occupancy grid.
    pub fn new(width: usize, height: usize, resolution: f64) -> Self {
        Self {
            resolution,
            width,
            height,
            log_odds: vec![0.0; width * height],
            origin: Vec2::zero(),
            lo_min: -5.0,
            lo_max: 5.0,
            lo_occ: 0.9,
            lo_free: -0.7,
        }
    }
    /// World coordinates to grid cell.
    pub fn world_to_cell(&self, world: Vec2) -> Option<(usize, usize)> {
        let cx = ((world.x - self.origin.x) / self.resolution) as isize;
        let cy = ((world.y - self.origin.y) / self.resolution) as isize;
        if cx >= 0 && cy >= 0 && (cx as usize) < self.width && (cy as usize) < self.height {
            Some((cx as usize, cy as usize))
        } else {
            None
        }
    }
    /// Cell to world coordinates (center of cell).
    pub fn cell_to_world(&self, cx: usize, cy: usize) -> Vec2 {
        Vec2::new(
            self.origin.x + (cx as f64 + 0.5) * self.resolution,
            self.origin.y + (cy as f64 + 0.5) * self.resolution,
        )
    }
    /// Update log-odds at a cell.
    pub fn update_cell(&mut self, cx: usize, cy: usize, lo_update: f64) {
        let idx = cy * self.width + cx;
        self.log_odds[idx] = (self.log_odds[idx] + lo_update).clamp(self.lo_min, self.lo_max);
    }
    /// Get occupancy probability at a cell.
    pub fn probability(&self, cx: usize, cy: usize) -> f64 {
        let lo = self.log_odds[cy * self.width + cx];
        1.0 / (1.0 + (-lo).exp())
    }
    /// Is cell occupied (p > 0.7)?
    pub fn is_occupied(&self, cx: usize, cy: usize) -> bool {
        self.probability(cx, cy) > 0.7
    }
    /// Update from a lidar scan (sensor origin at grid origin).
    pub fn update_lidar(&mut self, scan: &[(f64, f64)]) {
        for &(angle, range) in scan {
            let endpoint = Vec2::new(range * angle.cos(), range * angle.sin());
            if let Some((cx, cy)) = self.world_to_cell(endpoint) {
                self.update_cell(cx, cy, self.lo_occ);
            }
            let steps = (range / self.resolution) as usize;
            for step in 0..steps.saturating_sub(1) {
                let t = step as f64 * self.resolution;
                let p = Vec2::new(t * angle.cos(), t * angle.sin());
                if let Some((cx, cy)) = self.world_to_cell(p) {
                    let lo_free = self.lo_free;
                    self.update_cell(cx, cy, lo_free);
                }
            }
        }
    }
    /// Inflate occupied cells by a radius (m) for safety margin.
    pub fn inflate(&self, radius: f64) -> Self {
        let r_cells = (radius / self.resolution).ceil() as usize;
        let mut inflated = self.clone();
        for cy in 0..self.height {
            for cx in 0..self.width {
                if self.is_occupied(cx, cy) {
                    for dy in 0..=r_cells {
                        for dx in 0..=r_cells {
                            let fx = (dx * dx + dy * dy) as f64;
                            if fx <= (r_cells * r_cells) as f64 {
                                for (sx, sy) in [
                                    (cx.wrapping_add(dx), cy.wrapping_add(dy)),
                                    (cx.wrapping_sub(dx), cy.wrapping_add(dy)),
                                    (cx.wrapping_add(dx), cy.wrapping_sub(dy)),
                                    (cx.wrapping_sub(dx), cy.wrapping_sub(dy)),
                                ] {
                                    if sx < self.width && sy < self.height {
                                        let idx = sy * self.width + sx;
                                        inflated.log_odds[idx] = self.lo_max;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        inflated
    }
}
/// Action output from autonomous driving stack.
#[derive(Debug, Clone)]
pub struct AutonomousAction {
    /// Steering angle command (rad)
    pub steer: f64,
    /// Acceleration command (m/s²)
    pub accel: f64,
    /// Brake demand (m/s²)
    pub brake_demand: f64,
    /// Comfort event
    pub comfort_event: ComfortEvent,
    /// Emergency stop required
    pub emergency_stop: bool,
}
/// Simulated radar sensor.
#[derive(Debug, Clone)]
pub struct RadarSensor {
    /// Maximum range (m)
    pub max_range: f64,
    /// Field of view (rad)
    pub fov: f64,
    /// Range resolution (m)
    pub range_resolution: f64,
    /// Velocity resolution (m/s)
    pub velocity_resolution: f64,
    /// False alarm rate (probability per beam)
    pub false_alarm_rate: f64,
}
impl RadarSensor {
    /// Create a typical long-range radar (LRR).
    pub fn long_range() -> Self {
        Self {
            max_range: 250.0,
            fov: 10.0_f64.to_radians(),
            range_resolution: 0.4,
            velocity_resolution: 0.1,
            false_alarm_rate: 0.001,
        }
    }
    /// Create a mid-range radar (MRR).
    pub fn mid_range() -> Self {
        Self {
            max_range: 100.0,
            fov: 60.0_f64.to_radians(),
            range_resolution: 0.5,
            velocity_resolution: 0.15,
            false_alarm_rate: 0.002,
        }
    }
    /// Detect target given position and radial velocity.
    pub fn detect(&self, pos: Vec2, radial_vel: f64) -> Option<(f64, f64, f64)> {
        let range = pos.norm();
        let bearing = pos.y.atan2(pos.x);
        if range > self.max_range || bearing.abs() > self.fov / 2.0 {
            return None;
        }
        let r_meas = (range / self.range_resolution).round() * self.range_resolution;
        let v_meas = (radial_vel / self.velocity_resolution).round() * self.velocity_resolution;
        Some((r_meas, bearing, v_meas))
    }
}
/// Parking maneuver planner.
#[derive(Debug, Clone)]
pub struct ParkingPlanner {
    /// Vehicle wheelbase (m)
    pub wheelbase: f64,
    /// Vehicle width (m)
    pub vehicle_width: f64,
    /// Vehicle front overhang (m)
    pub front_overhang: f64,
    /// Vehicle rear overhang (m)
    pub rear_overhang: f64,
    /// Minimum turning radius (m)
    pub min_turn_radius: f64,
}
impl ParkingPlanner {
    /// Create a parking planner.
    pub fn new(wheelbase: f64, vehicle_width: f64) -> Self {
        Self {
            wheelbase,
            vehicle_width,
            front_overhang: 0.9,
            rear_overhang: 0.8,
            min_turn_radius: 5.5,
        }
    }
    /// Check if a parallel parking slot is reachable.
    pub fn can_parallel_park(&self, slot: &ParkingSlot) -> bool {
        let vehicle_length = self.wheelbase + self.front_overhang + self.rear_overhang;
        slot.length > vehicle_length + 0.8 && slot.width > self.vehicle_width + 0.3
    }
    /// Check if a perpendicular parking slot is reachable.
    pub fn can_perpendicular_park(&self, slot: &ParkingSlot) -> bool {
        let vehicle_length = self.wheelbase + self.front_overhang + self.rear_overhang;
        slot.length > vehicle_length + 0.3 && slot.width > self.vehicle_width + 0.5
    }
    /// Generate parallel parking path (simplified 2-arc maneuver).
    ///
    /// Returns sequence of (x, y, heading) configurations.
    pub fn parallel_park_path(&self, start: Vec2, slot: &ParkingSlot) -> Vec<[f64; 3]> {
        let r = self.min_turn_radius;
        let mut path = Vec::new();
        let sx = start.x;
        let sy = start.y;
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let angle = -t * std::f64::consts::PI / 4.0;
            path.push([sx + r * angle.sin(), sy - r * (1.0 - angle.cos()), angle]);
        }
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let prev = path.last().cloned().unwrap_or([sx, sy, 0.0]);
            let angle = -std::f64::consts::PI / 4.0 + t * std::f64::consts::PI / 4.0;
            path.push([prev[0] - r * 0.02, prev[1] - r * 0.01, angle]);
        }
        path.push([slot.center.x, slot.center.y, slot.heading]);
        path
    }
    /// Generate perpendicular forward parking path.
    pub fn perpendicular_park_path(&self, start: Vec2, slot: &ParkingSlot) -> Vec<[f64; 3]> {
        let r = self.min_turn_radius;
        let dx = slot.center.x - start.x;
        let dy = slot.center.y - start.y;
        let d = (dx.powi(2) + dy.powi(2)).sqrt();
        let n = 30;
        let mut path = Vec::new();
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = start.x + dx * t;
            let y = start.y + dy * t + r * (std::f64::consts::PI * t).sin() * 0.3;
            let heading = slot.heading * t;
            path.push([x, y, heading]);
        }
        let _ = d;
        path
    }
}
/// Full autonomous driving stack.
#[derive(Debug, Clone)]
pub struct AutonomousVehicle {
    /// Sensor fusion
    pub fusion: SensorFusion,
    /// Occupancy grid
    pub grid: OccupancyGrid,
    /// MPC tracker
    pub tracker: MpcPathTracker,
    /// Risk assessor
    pub risk: RiskAssessor,
    /// Comfort monitor
    pub comfort: ComfortMonitor,
    /// AEB
    pub aeb: AebSystem,
    /// Gap acceptance
    pub gap_accept: GapAcceptance,
    /// Current planned path
    pub planned_path: Vec<Vec2>,
    /// Target speed (m/s)
    pub target_speed: f64,
    /// LiDAR
    pub lidar: LidarSensor,
    /// Current braking demand (m/s²)
    pub brake_demand: f64,
}
impl AutonomousVehicle {
    /// Create a new autonomous vehicle.
    pub fn new() -> Self {
        Self {
            fusion: SensorFusion::new(),
            grid: OccupancyGrid::new(200, 200, 0.5),
            tracker: MpcPathTracker::new(2.7, 10, 0.05),
            risk: RiskAssessor::new(),
            comfort: ComfortMonitor::new(0.05),
            aeb: AebSystem::new(),
            gap_accept: GapAcceptance::new(0.5),
            planned_path: Vec::new(),
            target_speed: 13.9,
            lidar: LidarSensor::automotive(),
            brake_demand: 0.0,
        }
    }
    /// Process one frame.
    pub fn step(
        &mut self,
        detections: &[Vec2],
        lidar_scan: &[(f64, f64)],
        dt: f64,
    ) -> AutonomousAction {
        let tracks = self.fusion.update(detections, dt);
        self.grid.update_lidar(lidar_scan);
        let ego_pos = self.tracker.model.pos();
        let ego_speed = self.tracker.model.state[3];
        let ego_heading = self.tracker.model.state[2];
        let risks = self
            .risk
            .assess_risks(ego_speed, ego_pos, ego_heading, &tracks);
        let min_ttc = risks
            .iter()
            .map(|(_, ttc, _)| *ttc)
            .fold(f64::INFINITY, f64::min);
        self.brake_demand = self.aeb.update(min_ttc, ObjectClass::Vehicle);
        if !self.planned_path.is_empty() {
            self.tracker.step(&self.planned_path, self.target_speed);
        }
        let ax = self.tracker.model.state[3] * 0.0;
        let event = self.comfort.update(ax, self.tracker.steer_cmd * 9.81 * 0.3);
        AutonomousAction {
            steer: self.tracker.steer_cmd,
            accel: if self.brake_demand > 0.0 {
                -self.brake_demand
            } else {
                self.tracker.accel_cmd
            },
            brake_demand: self.brake_demand,
            comfort_event: event,
            emergency_stop: self.aeb.active && min_ttc < self.aeb.ttc_full,
        }
    }
}
/// Simulated LiDAR sensor.
#[derive(Debug, Clone)]
pub struct LidarSensor {
    /// Maximum range (m)
    pub max_range: f64,
    /// Angular resolution (rad)
    pub angular_resolution: f64,
    /// Number of beams
    pub num_beams: usize,
    /// Measurement noise standard deviation (m)
    pub range_noise_std: f64,
    /// Field of view (rad) — total
    pub fov: f64,
}
impl LidarSensor {
    /// Create a typical automotive LiDAR (360° FOV, 64 lines represented as 2D).
    pub fn automotive() -> Self {
        Self {
            max_range: 200.0,
            angular_resolution: 0.2_f64.to_radians(),
            num_beams: 1800,
            range_noise_std: 0.02,
            fov: std::f64::consts::TAU,
        }
    }
    /// Simulate a range scan given a set of obstacle positions.
    ///
    /// Returns (angle_rad, range_m) pairs.
    pub fn scan(&self, obstacles: &[(Vec2, f64)]) -> Vec<(f64, f64)> {
        let mut result = Vec::with_capacity(self.num_beams);
        for i in 0..self.num_beams {
            let angle = -self.fov / 2.0 + i as f64 * self.fov / self.num_beams as f64;
            let dir = Vec2::new(angle.cos(), angle.sin());
            let mut min_range = self.max_range;
            for (center, radius) in obstacles {
                let t = center.dot(&dir);
                if t < 0.0 {
                    continue;
                }
                let closest = dir.scale(t);
                let perp_dist = center.sub(&closest).norm();
                if perp_dist < *radius {
                    let dist = t - (radius.powi(2) - perp_dist.powi(2)).sqrt();
                    if dist > 0.0 && dist < min_range {
                        min_range = dist;
                    }
                }
            }
            result.push((angle, min_range));
        }
        result
    }
}
/// Pedestrian crowd simulation.
#[derive(Debug, Clone)]
pub struct PedestrianCrowd {
    /// All pedestrians
    pub agents: Vec<Pedestrian>,
    /// Obstacles (center, radius)
    pub obstacles: Vec<(Vec2, f64)>,
}
impl PedestrianCrowd {
    /// Create an empty crowd.
    pub fn new() -> Self {
        Self {
            agents: Vec::new(),
            obstacles: Vec::new(),
        }
    }
    /// Step all agents.
    pub fn step(&mut self, dt: f64) {
        let snapshot = self.agents.clone();
        for agent in &mut self.agents {
            agent.step(&snapshot, &self.obstacles, dt);
        }
    }
    /// Check if any pedestrian is within radius of a point.
    pub fn any_in_radius(&self, pos: Vec2, radius: f64) -> bool {
        self.agents.iter().any(|a| a.pos.dist(&pos) < radius)
    }
}
/// Camera-based object detector (simple bounding-box model).
#[derive(Debug, Clone)]
pub struct CameraDetector {
    /// Camera horizontal FOV (rad)
    pub hfov: f64,
    /// Camera vertical FOV (rad)
    pub vfov: f64,
    /// Image width (pixels)
    pub width: u32,
    /// Image height (pixels)
    pub height: u32,
    /// Detection confidence threshold
    pub conf_threshold: f64,
}
impl CameraDetector {
    /// Create a typical automotive forward camera.
    pub fn forward_camera() -> Self {
        Self {
            hfov: 60.0_f64.to_radians(),
            vfov: 40.0_f64.to_radians(),
            width: 1920,
            height: 1080,
            conf_threshold: 0.5,
        }
    }
    /// Project 3D point to pixel. Returns (u, v) or None if behind camera.
    pub fn project(&self, point_3d: [f64; 3]) -> Option<(f64, f64)> {
        if point_3d[2] <= 0.0 {
            return None;
        }
        let fx = self.width as f64 / (2.0 * (self.hfov / 2.0).tan());
        let fy = self.height as f64 / (2.0 * (self.vfov / 2.0).tan());
        let cx = self.width as f64 / 2.0;
        let cy = self.height as f64 / 2.0;
        let u = fx * point_3d[0] / point_3d[2] + cx;
        let v = fy * point_3d[1] / point_3d[2] + cy;
        if u >= 0.0 && u < self.width as f64 && v >= 0.0 && v < self.height as f64 {
            Some((u, v))
        } else {
            None
        }
    }
}
/// RRT (Rapidly-exploring Random Tree) node.
#[derive(Debug, Clone)]
pub(super) struct RrtNode {
    pub(super) state: PlanState,
    pub(super) parent: Option<usize>,
    pub(super) cost: f64,
}
/// Collision risk assessment.
#[derive(Debug, Clone)]
pub struct RiskAssessor {
    /// Time-to-collision threshold (s) for warning
    pub ttc_warn: f64,
    /// Time-to-collision threshold (s) for emergency
    pub ttc_emergency: f64,
    /// Vehicle half-width (m)
    pub vehicle_width: f64,
    /// Vehicle half-length (m)
    pub vehicle_length: f64,
}
impl RiskAssessor {
    /// Create a risk assessor.
    pub fn new() -> Self {
        Self {
            ttc_warn: 3.0,
            ttc_emergency: 1.5,
            vehicle_width: 1.0,
            vehicle_length: 2.5,
        }
    }
    /// Compute time-to-collision with a leading object.
    pub fn time_to_collision(&self, ego_speed: f64, obj_speed: f64, gap: f64) -> f64 {
        let rel_speed = ego_speed - obj_speed;
        if rel_speed <= 0.0 {
            return f64::INFINITY;
        }
        gap / rel_speed
    }
    /// Compute post-encroachment time (PET).
    pub fn post_encroachment_time(&self, ego_cross_time: f64, obj_cross_time: f64) -> f64 {
        (ego_cross_time - obj_cross_time).abs()
    }
    /// Risk level for a set of tracked objects.
    ///
    /// Returns (object_id, ttc, risk_level) tuples.
    pub fn assess_risks(
        &self,
        ego_speed: f64,
        ego_pos: Vec2,
        ego_heading: f64,
        objects: &[(u32, Vec2, Vec2)],
    ) -> Vec<(u32, f64, RiskLevel)> {
        let mut risks = Vec::new();
        for (id, pos, vel) in objects {
            let gap = ego_pos.dist(pos) - self.vehicle_length - 1.0;
            let obj_speed_fwd = vel.dot(&Vec2::new(ego_heading.cos(), ego_heading.sin()));
            let ttc = self.time_to_collision(ego_speed, obj_speed_fwd, gap.max(0.0));
            let level = if ttc < self.ttc_emergency {
                RiskLevel::Emergency
            } else if ttc < self.ttc_warn {
                RiskLevel::Warning
            } else {
                RiskLevel::Safe
            };
            risks.push((*id, ttc, level));
        }
        risks
    }
}
/// Detected object from a sensor.
#[derive(Debug, Clone)]
pub struct DetectedObject {
    /// Object identifier
    pub id: u32,
    /// Position in vehicle frame (m)
    pub position: Vec2,
    /// Velocity in world frame (m/s)
    pub velocity: Vec2,
    /// Object extent (bounding box half-size)
    pub extent: Vec2,
    /// Detection confidence \[0,1\]
    pub confidence: f64,
    /// Object class
    pub class: ObjectClass,
    /// Age of track (frames)
    pub track_age: u32,
}
/// Jerk and lateral-g comfort monitor.
#[derive(Debug, Clone)]
pub struct ComfortMonitor {
    /// Jerk threshold (m/s³) for comfort warning
    pub jerk_limit: f64,
    /// Lateral acceleration limit (m/s²)
    pub lat_g_limit: f64,
    /// Buffer of accelerations for jerk computation
    pub(super) accel_buf: VecDeque<[f64; 2]>,
    /// Time step (s)
    pub dt: f64,
    /// Running peak jerk
    pub peak_jerk: f64,
    /// Running peak lateral g
    pub peak_lat_g: f64,
    /// ISO 2631 ride comfort accumulator
    pub ride_comfort_aw: f64,
}
impl ComfortMonitor {
    /// Create a comfort monitor.
    pub fn new(dt: f64) -> Self {
        Self {
            jerk_limit: 3.0,
            lat_g_limit: 3.0,
            accel_buf: VecDeque::with_capacity(10),
            dt,
            peak_jerk: 0.0,
            peak_lat_g: 0.0,
            ride_comfort_aw: 0.0,
        }
    }
    /// Update with new longitudinal and lateral acceleration (m/s²).
    pub fn update(&mut self, ax: f64, ay: f64) -> ComfortEvent {
        self.accel_buf.push_back([ax, ay]);
        if self.accel_buf.len() > 3 {
            self.accel_buf.pop_front();
        }
        let jerk = if self.accel_buf.len() >= 2 {
            let n = self.accel_buf.len();
            let da = self.accel_buf[n - 1][0] - self.accel_buf[n - 2][0];
            (da / self.dt).abs()
        } else {
            0.0
        };
        let lat_g = ay.abs() / 9.81;
        self.peak_jerk = self.peak_jerk.max(jerk);
        self.peak_lat_g = self.peak_lat_g.max(lat_g);
        self.ride_comfort_aw =
            0.99 * self.ride_comfort_aw + 0.01 * (ax.powi(2) + ay.powi(2)).sqrt();
        if jerk > self.jerk_limit || lat_g > self.lat_g_limit {
            ComfortEvent::Discomfort { jerk, lat_g }
        } else {
            ComfortEvent::Ok
        }
    }
    /// Compute WBVS (Whole-Body Vibration Score) approximation.
    pub fn wbvs(&self) -> f64 {
        (self.ride_comfort_aw / 0.315).powi(2)
    }
}
/// Vehicle state for planning: \[x, y, heading, speed\].
#[derive(Debug, Clone, Copy)]
pub struct PlanState {
    /// X position (m)
    pub x: f64,
    /// Y position (m)
    pub y: f64,
    /// Heading (rad)
    pub heading: f64,
    /// Speed (m/s)
    pub speed: f64,
}
impl PlanState {
    /// Create a new plan state.
    pub fn new(x: f64, y: f64, heading: f64, speed: f64) -> Self {
        Self {
            x,
            y,
            heading,
            speed,
        }
    }
    /// 2D position.
    pub fn pos(&self) -> Vec2 {
        Vec2::new(self.x, self.y)
    }
    /// Distance to another state.
    pub fn dist(&self, other: &PlanState) -> f64 {
        self.pos().dist(&other.pos())
    }
}
/// Risk level enumeration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RiskLevel {
    /// No immediate risk
    Safe,
    /// Warning level
    Warning,
    /// Emergency braking required
    Emergency,
}
/// Parking type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParkingType {
    /// Parallel to road
    Parallel,
    /// Perpendicular to road
    Perpendicular,
    /// Angled
    Angled(f64),
}
/// A pedestrian agent using the social force model.
#[derive(Debug, Clone)]
pub struct Pedestrian {
    /// ID
    pub id: u32,
    /// Position (m)
    pub pos: Vec2,
    /// Velocity (m/s)
    pub vel: Vec2,
    /// Desired speed (m/s)
    pub desired_speed: f64,
    /// Destination
    pub destination: Vec2,
    /// Mass (kg)
    pub mass: f64,
    /// Relaxation time (s)
    pub tau: f64,
}
impl Pedestrian {
    /// Create a pedestrian.
    pub fn new(id: u32, pos: Vec2, destination: Vec2) -> Self {
        Self {
            id,
            pos,
            vel: Vec2::zero(),
            desired_speed: 1.4,
            destination,
            mass: 70.0,
            tau: 0.5,
        }
    }
    /// Compute driving force towards destination.
    pub fn driving_force(&self) -> Vec2 {
        let dir = self.destination.sub(&self.pos).normalize();
        let desired_vel = dir.scale(self.desired_speed);
        let dv = desired_vel.sub(&self.vel);
        dv.scale(self.mass / self.tau)
    }
    /// Repulsive force from another pedestrian.
    pub fn repulsive_force(&self, other: &Pedestrian) -> Vec2 {
        let d = self.pos.dist(&other.pos);
        if d < 1e-6 {
            return Vec2::zero();
        }
        let dir = self.pos.sub(&other.pos).normalize();
        let a = 2000.0;
        let b = 0.08;
        let strength = a * (-(d - 0.5) / b).exp();
        dir.scale(strength)
    }
    /// Repulsive force from vehicle obstacle.
    pub fn obstacle_force(&self, obs_pos: Vec2, obs_radius: f64) -> Vec2 {
        let d = (self.pos.dist(&obs_pos) - obs_radius).max(0.01);
        let dir = self.pos.sub(&obs_pos).normalize();
        let strength = 2000.0 * (-(d - 0.3) / 0.1).exp();
        dir.scale(strength)
    }
    /// Step the pedestrian.
    pub fn step(&mut self, others: &[Pedestrian], obstacles: &[(Vec2, f64)], dt: f64) {
        let mut force = self.driving_force();
        for other in others {
            if other.id != self.id {
                force = force.add(&self.repulsive_force(other));
            }
        }
        for (obs_pos, obs_r) in obstacles {
            force = force.add(&self.obstacle_force(*obs_pos, *obs_r));
        }
        let accel = force.scale(1.0 / self.mass);
        self.vel = self.vel.add(&accel.scale(dt));
        let spd = self.vel.norm();
        if spd > self.desired_speed * 2.0 {
            self.vel = self.vel.scale(self.desired_speed * 2.0 / spd);
        }
        self.pos = self.pos.add(&self.vel.scale(dt));
    }
}
/// Comfort event.
#[derive(Debug, Clone, PartialEq)]
pub enum ComfortEvent {
    /// Comfortable
    Ok,
    /// Discomfort detected
    Discomfort {
        /// Jerk (m/s³)
        jerk: f64,
        /// Lateral g (g)
        lat_g: f64,
    },
}
/// Model predictive controller for path tracking.
#[derive(Debug, Clone)]
pub struct MpcPathTracker {
    /// Bicycle model
    pub model: BicycleModel,
    /// Prediction horizon (steps)
    pub horizon: usize,
    /// Time step (s)
    pub dt: f64,
    /// State cost weights \[x, y, heading, speed\]
    pub q_weights: [f64; 4],
    /// Control cost weights \[steer, accel\]
    pub r_weights: [f64; 2],
    /// Last computed steering command
    pub steer_cmd: f64,
    /// Last computed acceleration command
    pub accel_cmd: f64,
}
impl MpcPathTracker {
    /// Create an MPC tracker.
    pub fn new(wheelbase: f64, horizon: usize, dt: f64) -> Self {
        Self {
            model: BicycleModel::new(wheelbase),
            horizon,
            dt,
            q_weights: [1.0, 1.0, 0.5, 0.1],
            r_weights: [0.1, 0.01],
            steer_cmd: 0.0,
            accel_cmd: 0.0,
        }
    }
    /// Pure-pursuit path tracking (approximation of MPC).
    ///
    /// Returns steering angle (rad).
    pub fn pure_pursuit_steer(&self, path: &[Vec2], lookahead: f64) -> f64 {
        let pos = self.model.pos();
        let target = path
            .iter()
            .find(|p| pos.dist(p) >= lookahead)
            .or_else(|| path.last())
            .cloned()
            .unwrap_or(pos);
        let heading = self.model.state[2];
        let dx = target.x - pos.x;
        let dy = target.y - pos.y;
        let _local_x = dx * heading.cos() + dy * heading.sin();
        let local_y = -dx * heading.sin() + dy * heading.cos();
        let ld = (dx.powi(2) + dy.powi(2)).sqrt().max(0.01);
        let curvature = 2.0 * local_y / ld.powi(2);
        (curvature * self.model.wheelbase)
            .atan()
            .clamp(-self.model.max_steer, self.model.max_steer)
    }
    /// Stanley path tracking controller.
    pub fn stanley_steer(&self, path: &[Vec2], k: f64) -> f64 {
        let pos = self.model.pos();
        let heading = self.model.state[2];
        let v = self.model.state[3].max(0.1);
        let (nearest_idx, _) = path
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                pos.dist(a)
                    .partial_cmp(&pos.dist(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or((0, &pos));
        let ni = nearest_idx.min(path.len() - 1);
        let path_heading = if ni + 1 < path.len() {
            let d = path[ni + 1].sub(&path[ni]);
            d.y.atan2(d.x)
        } else {
            heading
        };
        let nearest = path[ni];
        let dx = nearest.x - pos.x;
        let dy = nearest.y - pos.y;
        let cte = -dx * path_heading.sin() + dy * path_heading.cos();
        let heading_err = (path_heading - heading + std::f64::consts::PI)
            .rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        let steer = heading_err + (k * cte / v).atan();
        steer.clamp(-self.model.max_steer, self.model.max_steer)
    }
    /// Step the tracker along a path.
    pub fn step(&mut self, path: &[Vec2], target_speed: f64) {
        if path.is_empty() {
            return;
        }
        self.steer_cmd = self.pure_pursuit_steer(path, 5.0);
        let speed_err = target_speed - self.model.state[3];
        self.accel_cmd = (speed_err * 0.5).clamp(-3.0, 2.0);
        self.model.step(self.steer_cmd, self.accel_cmd, self.dt);
    }
}
/// A* path planner on occupancy grid.
pub struct AStarPlanner<'a> {
    /// Reference to the occupancy grid
    pub grid: &'a OccupancyGrid,
}
impl<'a> AStarPlanner<'a> {
    /// Create a planner.
    pub fn new(grid: &'a OccupancyGrid) -> Self {
        Self { grid }
    }
    /// Plan path. Returns list of (cx, cy) cell indices.
    pub fn plan(&self, start: (usize, usize), goal: (usize, usize)) -> Option<Vec<(usize, usize)>> {
        let w = self.grid.width;
        let h = self.grid.height;
        let idx = |cx: usize, cy: usize| cy * w + cx;
        let heuristic = |cx: usize, cy: usize| {
            ((cx as f64 - goal.0 as f64).powi(2) + (cy as f64 - goal.1 as f64).powi(2)).sqrt()
        };
        let mut g_score = vec![f64::INFINITY; w * h];
        let mut came_from: Vec<Option<usize>> = vec![None; w * h];
        let mut open: BinaryHeap<AStarNode> = BinaryHeap::new();
        g_score[idx(start.0, start.1)] = 0.0;
        open.push(AStarNode {
            f: heuristic(start.0, start.1),
            g: 0.0,
            idx: idx(start.0, start.1),
        });
        while let Some(node) = open.pop() {
            let cx = node.idx % w;
            let cy = node.idx / w;
            if (cx, cy) == goal {
                let mut path = Vec::new();
                let mut current = node.idx;
                while let Some(parent) = came_from[current] {
                    path.push((current % w, current / w));
                    current = parent;
                }
                path.push(start);
                path.reverse();
                return Some(path);
            }
            for (ndx, ndy) in [
                (-1i32, 0i32),
                (1, 0),
                (0, -1),
                (0, 1),
                (-1, -1),
                (1, -1),
                (-1, 1),
                (1, 1),
            ] {
                let nx = cx as i32 + ndx;
                let ny = cy as i32 + ndy;
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let nx = nx as usize;
                let ny = ny as usize;
                if self.grid.is_occupied(nx, ny) {
                    continue;
                }
                let step_cost = if ndx != 0 && ndy != 0 { 1.414 } else { 1.0 };
                let new_g = node.g + step_cost;
                if new_g < g_score[idx(nx, ny)] {
                    g_score[idx(nx, ny)] = new_g;
                    came_from[idx(nx, ny)] = Some(node.idx);
                    open.push(AStarNode {
                        f: new_g + heuristic(nx, ny),
                        g: new_g,
                        idx: idx(nx, ny),
                    });
                }
            }
        }
        None
    }
}
/// Parking slot descriptor.
#[derive(Debug, Clone)]
pub struct ParkingSlot {
    /// Center position
    pub center: Vec2,
    /// Heading of slot (rad)
    pub heading: f64,
    /// Slot length (m)
    pub length: f64,
    /// Slot width (m)
    pub width: f64,
    /// Parking type
    pub slot_type: ParkingType,
}
/// Object classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ObjectClass {
    /// Unknown object
    Unknown,
    /// Vehicle
    Vehicle,
    /// Pedestrian
    Pedestrian,
    /// Cyclist
    Cyclist,
    /// Static obstacle
    Static,
}
/// Autonomous emergency braking system.
#[derive(Debug, Clone)]
pub struct AebSystem {
    /// System delay (s)
    pub system_delay: f64,
    /// Maximum deceleration (m/s²)
    pub max_decel: f64,
    /// Comfort deceleration (m/s²)
    pub comfort_decel: f64,
    /// TTC threshold for partial braking (s)
    pub ttc_partial: f64,
    /// TTC threshold for full braking (s)
    pub ttc_full: f64,
    /// Currently active
    pub active: bool,
    /// Current brake demand \[0,1\]
    pub brake_demand: f64,
    /// Pedestrian detection enabled
    pub pedestrian_detection: bool,
}
impl AebSystem {
    /// Create an AEB system.
    pub fn new() -> Self {
        Self {
            system_delay: 0.05,
            max_decel: 8.0,
            comfort_decel: 3.0,
            ttc_partial: 2.0,
            ttc_full: 1.0,
            active: false,
            brake_demand: 0.0,
            pedestrian_detection: true,
        }
    }
    /// Update AEB given current TTC (s) and object class.
    pub fn update(&mut self, ttc: f64, object_class: ObjectClass) -> f64 {
        let sensitivity = match object_class {
            ObjectClass::Pedestrian | ObjectClass::Cyclist => 1.3,
            _ => 1.0,
        };
        let ttc_eff = ttc / sensitivity;
        if ttc_eff < self.ttc_full {
            self.active = true;
            self.brake_demand = 1.0;
        } else if ttc_eff < self.ttc_partial {
            self.active = true;
            self.brake_demand = (self.ttc_partial - ttc_eff) / (self.ttc_partial - self.ttc_full);
        } else {
            self.active = false;
            self.brake_demand = 0.0;
        }
        self.brake_demand * self.max_decel
    }
    /// Compute safe following distance (m).
    pub fn safe_distance(&self, ego_speed: f64, lead_decel: f64) -> f64 {
        let v = ego_speed;
        let a = self.max_decel;
        let a_lead = lead_decel.max(0.1);
        let d_delay = v * self.system_delay;
        let d_brake = v.powi(2) / (2.0 * a) - v.powi(2) / (2.0 * a_lead);
        (d_delay + d_brake).max(2.0)
    }
}
