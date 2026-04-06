#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Race line optimization and lap simulation.
//!
//! Provides track geometry, vehicle dynamics modeling, lap time simulation,
//! and geometric race line optimization for motorsport applications.

/// A segment of a racing track.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum TrackSegment {
    /// A straight section with a given length in metres.
    Straight {
        /// Length of the straight in metres.
        length: f64,
    },
    /// A corner with a given radius, sweep angle, and direction.
    Corner {
        /// Corner radius in metres.
        radius: f64,
        /// Sweep angle in radians.
        angle_rad: f64,
        /// Direction: 1 = left, -1 = right.
        direction: i8,
    },
}

impl TrackSegment {
    /// Returns the arc length of this segment in metres.
    pub fn length(&self) -> f64 {
        match self {
            TrackSegment::Straight { length } => *length,
            TrackSegment::Corner {
                radius, angle_rad, ..
            } => radius * angle_rad.abs(),
        }
    }
}

/// A racing circuit composed of track segments.
#[derive(Debug, Clone, Default)]
pub struct Track {
    /// Ordered list of track segments.
    pub segments: Vec<TrackSegment>,
}

impl Track {
    /// Creates an empty track.
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a straight section.
    pub fn add_straight(&mut self, length: f64) {
        self.segments.push(TrackSegment::Straight { length });
    }

    /// Appends a corner.
    ///
    /// `angle_deg` is the sweep angle in degrees. `left` selects direction.
    pub fn add_corner(&mut self, radius: f64, angle_deg: f64, left: bool) {
        let angle_rad = angle_deg.to_radians();
        let direction: i8 = if left { 1 } else { -1 };
        self.segments.push(TrackSegment::Corner {
            radius,
            angle_rad,
            direction,
        });
    }

    /// Returns the total arc length of the track in metres.
    pub fn total_length(&self) -> f64 {
        self.segments.iter().map(|s| s.length()).sum()
    }

    /// Returns the minimum corner radius found on the track, or `f64::INFINITY` if none.
    pub fn min_corner_radius(&self) -> f64 {
        self.segments
            .iter()
            .filter_map(|s| match s {
                TrackSegment::Corner { radius, .. } => Some(*radius),
                _ => None,
            })
            .fold(f64::INFINITY, f64::min)
    }
}

/// Point-mass vehicle dynamics parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VehicleDynamics {
    /// Vehicle mass in kg.
    pub mass: f64,
    /// Maximum engine force in N.
    pub max_engine_force: f64,
    /// Maximum brake force in N.
    pub max_brake_force: f64,
    /// Aerodynamic drag coefficient (F_drag = coeff * v²).
    pub aero_drag_coeff: f64,
    /// Maximum lateral acceleration in units of g.
    pub max_lateral_g: f64,
}

const G: f64 = 9.81;

impl VehicleDynamics {
    /// Typical open-wheel formula car parameters.
    pub fn formula_car() -> Self {
        Self {
            mass: 740.0,
            max_engine_force: 8000.0,
            max_brake_force: 16000.0,
            aero_drag_coeff: 0.9,
            max_lateral_g: 4.5,
        }
    }

    /// Typical sports car parameters.
    pub fn sports_car() -> Self {
        Self {
            mass: 1400.0,
            max_engine_force: 5000.0,
            max_brake_force: 10000.0,
            aero_drag_coeff: 0.5,
            max_lateral_g: 1.8,
        }
    }

    /// Maximum speed achievable in a corner of given radius.
    ///
    /// Uses `v = sqrt(mu * g * radius)` where `mu = max_lateral_g`.
    pub fn max_speed_in_corner(&self, radius: f64) -> f64 {
        (self.max_lateral_g * G * radius).sqrt()
    }

    /// Speed after one integration step with engine force minus aero drag.
    ///
    /// `v` is current speed (m/s), `dt` is time step (s).
    pub fn acceleration_limited_speed(&self, v: f64, dt: f64) -> f64 {
        let f_drag = self.aero_drag_coeff * v * v;
        let f_net = self.max_engine_force - f_drag;
        let a = f_net / self.mass;
        (v + a * dt).max(0.0)
    }

    /// Distance required to brake from `v` down to `v_target` (m).
    ///
    /// Uses kinematics with maximum deceleration from `max_brake_force`.
    pub fn braking_distance(&self, v: f64, v_target: f64) -> f64 {
        if v <= v_target {
            return 0.0;
        }
        let a_brake = self.max_brake_force / self.mass;
        (v * v - v_target * v_target) / (2.0 * a_brake)
    }
}

/// Result of a lap simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LapResult {
    /// Total lap time in seconds.
    pub total_time: f64,
    /// Maximum speed reached during the lap (m/s).
    pub max_speed: f64,
    /// Average speed over the lap (m/s).
    pub avg_speed: f64,
    /// Time spent in each segment (s).
    pub sector_times: Vec<f64>,
}

/// Point-mass lap time simulator.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LapSimulator {
    track: Track,
    vehicle: VehicleDynamics,
}

impl LapSimulator {
    /// Creates a new simulator for the given track and vehicle.
    pub fn new(track: Track, vehicle: VehicleDynamics) -> Self {
        Self { track, vehicle }
    }

    /// Computes the maximum entry speed for a segment.
    fn segment_max_speed(&self, seg: &TrackSegment) -> f64 {
        match seg {
            TrackSegment::Straight { .. } => {
                // Terminal velocity: F_engine = F_drag => v = sqrt(F_engine / coeff)
                (self.vehicle.max_engine_force / self.vehicle.aero_drag_coeff).sqrt()
            }
            TrackSegment::Corner { radius, .. } => self.vehicle.max_speed_in_corner(*radius),
        }
    }

    /// Runs the lap simulation and returns a `LapResult`.
    pub fn simulate(&self) -> LapResult {
        let n = self.track.segments.len();
        if n == 0 {
            return LapResult {
                total_time: 0.0,
                max_speed: 0.0,
                avg_speed: 0.0,
                sector_times: vec![],
            };
        }

        // Build speed limit per segment.
        let speed_limits: Vec<f64> = self
            .track
            .segments
            .iter()
            .map(|s| self.segment_max_speed(s))
            .collect();

        // Forward pass – limited by engine/acceleration.
        let mut speeds = vec![0.0f64; n];
        speeds[0] = speed_limits[0].min(10.0); // standing start
        for i in 1..n {
            let len = self.track.segments[i - 1].length().max(1.0);
            let dt = len / speeds[i - 1].max(1.0);
            let accel_v = self.vehicle.acceleration_limited_speed(speeds[i - 1], dt);
            speeds[i] = accel_v.min(speed_limits[i]);
        }

        // Backward pass – braking constraints.
        for i in (0..n - 1).rev() {
            let v_next = speeds[i + 1];
            let seg_len = self.track.segments[i].length();
            let bd = self.vehicle.braking_distance(speeds[i], v_next);
            if bd > seg_len {
                // Must enter slower.
                let a_brake = self.vehicle.max_brake_force / self.vehicle.mass;
                speeds[i] = (v_next * v_next + 2.0 * a_brake * seg_len).sqrt();
            }
            speeds[i] = speeds[i].min(speed_limits[i]);
        }

        // Compute sector times.
        let mut sector_times = Vec::with_capacity(n);
        let mut max_speed: f64 = 0.0;
        let mut total_dist = 0.0f64;

        for (i, seg) in self.track.segments.iter().enumerate() {
            let len = seg.length();
            let v_entry = speeds[i];
            let v_exit = if i + 1 < n { speeds[i + 1] } else { speeds[i] };
            let avg = ((v_entry + v_exit) / 2.0).max(1.0);
            let t = len / avg;
            sector_times.push(t);
            max_speed = max_speed.max(v_entry).max(v_exit);
            total_dist += len;
        }

        let total_time: f64 = sector_times.iter().sum();
        let avg_speed = if total_time > 0.0 {
            total_dist / total_time
        } else {
            0.0
        };

        LapResult {
            total_time,
            max_speed,
            avg_speed,
            sector_times,
        }
    }

    /// Returns the simulated time per segment in seconds.
    pub fn sector_times(&self) -> Vec<f64> {
        self.simulate().sector_times
    }
}

/// Geometric race line optimizer.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RaceLineOptimizer {
    track: Track,
}

impl RaceLineOptimizer {
    /// Creates a new optimizer for the given track.
    pub fn new(track: Track) -> Self {
        Self { track }
    }

    /// Returns a simplified geometric racing line as 2D waypoints.
    ///
    /// For straight segments the waypoint is placed at the centre.
    /// For corners the apex is placed at the inner edge.
    pub fn geometric_racing_line(&self) -> Vec<[f64; 2]> {
        let mut waypoints = Vec::new();
        let mut x = 0.0f64;
        let mut y = 0.0f64;
        let mut heading = 0.0f64; // radians

        for seg in &self.track.segments {
            match seg {
                TrackSegment::Straight { length } => {
                    x += heading.cos() * length;
                    y += heading.sin() * length;
                    waypoints.push([x, y]);
                }
                TrackSegment::Corner {
                    radius,
                    angle_rad,
                    direction,
                } => {
                    // Apex at inner edge: offset inward by `radius`.
                    let lateral = *direction as f64; // +1 = left, -1 = right
                    let perp_heading = heading + lateral * std::f64::consts::FRAC_PI_2;
                    let apex_x = x + perp_heading.cos() * radius;
                    let apex_y = y + perp_heading.sin() * radius;
                    waypoints.push([apex_x, apex_y]);

                    // Advance position along the arc.
                    let arc_len = radius * angle_rad.abs();
                    x += heading.cos() * arc_len;
                    y += heading.sin() * arc_len;
                    heading += lateral * angle_rad.abs();
                }
            }
        }

        waypoints
    }

    /// Applies Chaikin curve subdivision smoothing for `iterations` passes.
    pub fn smoothed_line(&self, waypoints: &[[f64; 2]], iterations: usize) -> Vec<[f64; 2]> {
        let mut pts: Vec<[f64; 2]> = waypoints.to_vec();
        for _ in 0..iterations {
            if pts.len() < 2 {
                break;
            }
            let mut next = Vec::with_capacity(pts.len() * 2);
            for i in 0..pts.len() - 1 {
                let p = pts[i];
                let q = pts[i + 1];
                next.push([0.75 * p[0] + 0.25 * q[0], 0.75 * p[1] + 0.25 * q[1]]);
                next.push([0.25 * p[0] + 0.75 * q[0], 0.25 * p[1] + 0.75 * q[1]]);
            }
            pts = next;
        }
        pts
    }
}

/// Estimates fuel consumption for a given trip.
///
/// Returns kilograms of fuel consumed.
///
/// # Arguments
/// - `distance_m`: trip distance in metres
/// - `avg_speed_mps`: average speed in m/s
/// - `specific_consumption_kg_per_mj`: fuel mass per megajoule of energy (kg/MJ)
/// - `power_kw`: average power output in kilowatts
#[allow(dead_code)]
pub fn fuel_consumption(
    distance_m: f64,
    avg_speed_mps: f64,
    specific_consumption_kg_per_mj: f64,
    power_kw: f64,
) -> f64 {
    if avg_speed_mps <= 0.0 {
        return 0.0;
    }
    let duration_s = distance_m / avg_speed_mps;
    let energy_mj = power_kw * duration_s / 1000.0; // kW * s -> kJ -> /1000 -> MJ
    specific_consumption_kg_per_mj * energy_mj
}

/// Computes tire wear rate per metre.
///
/// Returns a dimensionless wear index.
///
/// # Arguments
/// - `lateral_g`: lateral acceleration in g
/// - `speed`: vehicle speed in m/s
/// - `tire_compound_factor`: compound durability factor (higher = more wear)
#[allow(dead_code)]
pub fn tire_wear_model(lateral_g: f64, speed: f64, tire_compound_factor: f64) -> f64 {
    let lateral_g_sq = lateral_g * lateral_g;
    let speed_factor = speed / 100.0; // normalise to ~100 m/s reference
    tire_compound_factor * lateral_g_sq * (1.0 + speed_factor)
}

// ---------------------------------------------------------------------------
// Track segment classification
// ---------------------------------------------------------------------------

/// Classification of a track segment by its characteristics.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SegmentClass {
    /// Long straight (> 200 m).
    LongStraight,
    /// Short straight (<= 200 m).
    ShortStraight,
    /// High-speed corner (radius > 100 m).
    HighSpeedCorner,
    /// Medium-speed corner (50 m < radius <= 100 m).
    MediumSpeedCorner,
    /// Low-speed corner (radius <= 50 m).
    LowSpeedCorner,
    /// Hairpin (radius <= 20 m and angle >= 120 degrees).
    Hairpin,
}

/// Classify a track segment based on its geometry.
#[allow(dead_code)]
pub fn classify_segment(seg: &TrackSegment) -> SegmentClass {
    match seg {
        TrackSegment::Straight { length } => {
            if *length > 200.0 {
                SegmentClass::LongStraight
            } else {
                SegmentClass::ShortStraight
            }
        }
        TrackSegment::Corner {
            radius, angle_rad, ..
        } => {
            let angle_deg = angle_rad.abs().to_degrees();
            if *radius <= 20.0 && angle_deg >= 120.0 {
                SegmentClass::Hairpin
            } else if *radius <= 50.0 {
                SegmentClass::LowSpeedCorner
            } else if *radius <= 100.0 {
                SegmentClass::MediumSpeedCorner
            } else {
                SegmentClass::HighSpeedCorner
            }
        }
    }
}

/// Classify all segments of a track.
#[allow(dead_code)]
pub fn classify_track(track: &Track) -> Vec<SegmentClass> {
    track.segments.iter().map(classify_segment).collect()
}

// ---------------------------------------------------------------------------
// Apex detection
// ---------------------------------------------------------------------------

/// An apex point detected on the racing line.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ApexPoint {
    /// Segment index in the track.
    pub segment_index: usize,
    /// Position along the racing line.
    pub position: [f64; 2],
    /// Curvature at the apex (1/radius).
    pub curvature: f64,
    /// Whether this is a geometric apex (innermost point).
    pub is_geometric: bool,
}

/// Detect apex points on the racing line from track geometry.
///
/// Returns one apex per corner segment.
#[allow(dead_code)]
pub fn detect_apexes(track: &Track) -> Vec<ApexPoint> {
    let mut apexes = Vec::new();
    let mut x = 0.0_f64;
    let mut y = 0.0_f64;
    let mut heading = 0.0_f64;

    for (i, seg) in track.segments.iter().enumerate() {
        match seg {
            TrackSegment::Straight { length } => {
                x += heading.cos() * length;
                y += heading.sin() * length;
            }
            TrackSegment::Corner {
                radius,
                angle_rad,
                direction,
            } => {
                let lateral = *direction as f64;
                let half_angle = angle_rad.abs() / 2.0;
                // Apex is at the midpoint of the arc, offset inward
                let mid_heading = heading + lateral * half_angle;
                let perp = mid_heading + lateral * std::f64::consts::FRAC_PI_2;
                let apex_x =
                    x + mid_heading.cos() * radius * half_angle + perp.cos() * radius * 0.1;
                let apex_y =
                    y + mid_heading.sin() * radius * half_angle + perp.sin() * radius * 0.1;

                apexes.push(ApexPoint {
                    segment_index: i,
                    position: [apex_x, apex_y],
                    curvature: 1.0 / radius,
                    is_geometric: true,
                });

                let arc_len = radius * angle_rad.abs();
                x += heading.cos() * arc_len;
                y += heading.sin() * arc_len;
                heading += lateral * angle_rad.abs();
            }
        }
    }

    apexes
}

// ---------------------------------------------------------------------------
// Braking point calculation
// ---------------------------------------------------------------------------

/// A braking point on the track.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BrakingPoint {
    /// Segment index before which braking must begin.
    pub corner_index: usize,
    /// Distance before the corner entry where braking starts (m).
    pub braking_distance: f64,
    /// Speed at the start of braking (m/s).
    pub speed_at_braking: f64,
    /// Target corner entry speed (m/s).
    pub target_speed: f64,
}

/// Calculate braking points for all corners on the track.
#[allow(dead_code)]
pub fn compute_braking_points(track: &Track, vehicle: &VehicleDynamics) -> Vec<BrakingPoint> {
    let n = track.segments.len();
    if n == 0 {
        return Vec::new();
    }

    // Forward pass: compute achievable speed at each segment entry.
    let mut entry_speeds = vec![10.0_f64]; // standing start
    for i in 1..n {
        let seg_len = track.segments[i - 1].length().max(1.0);
        let dt = seg_len / entry_speeds[i - 1].max(1.0);
        let accel_v = vehicle.acceleration_limited_speed(entry_speeds[i - 1], dt);
        let limit = match &track.segments[i] {
            TrackSegment::Straight { .. } => {
                (vehicle.max_engine_force / vehicle.aero_drag_coeff).sqrt()
            }
            TrackSegment::Corner { radius, .. } => vehicle.max_speed_in_corner(*radius),
        };
        entry_speeds.push(accel_v.min(limit));
    }

    let mut braking_points = Vec::new();

    for i in 0..n {
        if let TrackSegment::Corner { radius, .. } = &track.segments[i] {
            let target = vehicle.max_speed_in_corner(*radius);
            // Look backward through preceding straights
            let approach_speed = if i > 0 { entry_speeds[i] } else { 10.0 };
            let bd = vehicle.braking_distance(approach_speed, target);

            braking_points.push(BrakingPoint {
                corner_index: i,
                braking_distance: bd,
                speed_at_braking: approach_speed,
                target_speed: target,
            });
        }
    }

    braking_points
}

// ---------------------------------------------------------------------------
// Racing line smoothing (Gaussian)
// ---------------------------------------------------------------------------

/// Apply Gaussian smoothing to a racing line.
///
/// Uses a 1D Gaussian kernel with the given sigma (in number of points).
#[allow(dead_code)]
pub fn gaussian_smooth(waypoints: &[[f64; 2]], sigma: f64) -> Vec<[f64; 2]> {
    if waypoints.len() < 3 || sigma <= 0.0 {
        return waypoints.to_vec();
    }

    let n = waypoints.len();
    let kernel_half = (3.0 * sigma).ceil() as usize;
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut sum_w = 0.0;

        let start = i.saturating_sub(kernel_half);
        let end = (i + kernel_half + 1).min(n);

        for j in start..end {
            let d = (j as f64) - (i as f64);
            let w = (-0.5 * d * d / (sigma * sigma)).exp();
            sum_x += w * waypoints[j][0];
            sum_y += w * waypoints[j][1];
            sum_w += w;
        }

        if sum_w > 1e-15 {
            result.push([sum_x / sum_w, sum_y / sum_w]);
        } else {
            result.push(waypoints[i]);
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Optimal racing line computation (curvature minimization)
// ---------------------------------------------------------------------------

/// Compute the curvature at each point of a polyline.
///
/// Uses three-point curvature estimation. Returns curvature for each interior
/// point (endpoints get zero curvature).
#[allow(dead_code)]
pub fn polyline_curvature(points: &[[f64; 2]]) -> Vec<f64> {
    let n = points.len();
    if n < 3 {
        return vec![0.0; n];
    }

    let mut curvatures = vec![0.0; n];
    for i in 1..n - 1 {
        let p0 = points[i - 1];
        let p1 = points[i];
        let p2 = points[i + 1];

        let ax = p1[0] - p0[0];
        let ay = p1[1] - p0[1];
        let bx = p2[0] - p1[0];
        let by = p2[1] - p1[1];

        let cross = ax * by - ay * bx;
        let la = (ax * ax + ay * ay).sqrt();
        let lb = (bx * bx + by * by).sqrt();
        let lc = ((p2[0] - p0[0]).powi(2) + (p2[1] - p0[1]).powi(2)).sqrt();

        let denom = la * lb * lc;
        if denom > 1e-15 {
            curvatures[i] = (2.0 * cross.abs()) / denom;
        }
    }

    curvatures
}

/// Compute the total curvature (sum of absolute curvatures) of a polyline.
#[allow(dead_code)]
pub fn total_curvature(points: &[[f64; 2]]) -> f64 {
    polyline_curvature(points).iter().sum()
}

/// Compute the total path length of a polyline.
#[allow(dead_code)]
pub fn path_length(points: &[[f64; 2]]) -> f64 {
    if points.len() < 2 {
        return 0.0;
    }
    let mut len = 0.0;
    for i in 1..points.len() {
        let dx = points[i][0] - points[i - 1][0];
        let dy = points[i][1] - points[i - 1][1];
        len += (dx * dx + dy * dy).sqrt();
    }
    len
}

/// Iteratively optimize a racing line by moving points toward the inside
/// of corners to reduce curvature (gradient descent on curvature).
///
/// `alpha` is the learning rate, `iterations` is the number of passes.
#[allow(dead_code)]
pub fn optimize_racing_line(
    waypoints: &[[f64; 2]],
    alpha: f64,
    iterations: usize,
) -> Vec<[f64; 2]> {
    if waypoints.len() < 4 {
        return waypoints.to_vec();
    }

    let mut pts = waypoints.to_vec();
    let n = pts.len();

    for _ in 0..iterations {
        let curvatures = polyline_curvature(&pts);
        let mut new_pts = pts.clone();

        // Don't move first and last points (track boundary conditions)
        for i in 1..n - 1 {
            if curvatures[i] < 1e-10 {
                continue;
            }

            // Move toward the midpoint of neighbors to reduce curvature
            let mid_x = (pts[i - 1][0] + pts[i + 1][0]) * 0.5;
            let mid_y = (pts[i - 1][1] + pts[i + 1][1]) * 0.5;

            let dx = mid_x - pts[i][0];
            let dy = mid_y - pts[i][1];

            new_pts[i][0] += alpha * dx;
            new_pts[i][1] += alpha * dy;
        }

        pts = new_pts;
    }

    pts
}

// ---------------------------------------------------------------------------
// Track analysis utilities
// ---------------------------------------------------------------------------

/// Count the number of corners in a track.
#[allow(dead_code)]
pub fn count_corners(track: &Track) -> usize {
    track
        .segments
        .iter()
        .filter(|s| matches!(s, TrackSegment::Corner { .. }))
        .count()
}

/// Count the number of straights in a track.
#[allow(dead_code)]
pub fn count_straights(track: &Track) -> usize {
    track
        .segments
        .iter()
        .filter(|s| matches!(s, TrackSegment::Straight { .. }))
        .count()
}

/// Maximum corner radius on the track.
#[allow(dead_code)]
pub fn max_corner_radius(track: &Track) -> f64 {
    track
        .segments
        .iter()
        .filter_map(|s| match s {
            TrackSegment::Corner { radius, .. } => Some(*radius),
            _ => None,
        })
        .fold(0.0_f64, f64::max)
}

/// Average corner radius.
#[allow(dead_code)]
pub fn avg_corner_radius(track: &Track) -> f64 {
    let radii: Vec<f64> = track
        .segments
        .iter()
        .filter_map(|s| match s {
            TrackSegment::Corner { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.is_empty() {
        return 0.0;
    }
    radii.iter().sum::<f64>() / radii.len() as f64
}

/// Compute the percentage of track length that is corners.
#[allow(dead_code)]
pub fn corner_percentage(track: &Track) -> f64 {
    let total = track.total_length();
    if total <= 0.0 {
        return 0.0;
    }
    let corner_len: f64 = track
        .segments
        .iter()
        .filter_map(|s| match s {
            TrackSegment::Corner { .. } => Some(s.length()),
            _ => None,
        })
        .sum();
    corner_len / total * 100.0
}

/// Compute the estimated top speed on the track.
#[allow(dead_code)]
pub fn estimated_top_speed(track: &Track, vehicle: &VehicleDynamics) -> f64 {
    let terminal = (vehicle.max_engine_force / vehicle.aero_drag_coeff).sqrt();
    let longest_straight = track
        .segments
        .iter()
        .filter_map(|s| match s {
            TrackSegment::Straight { length } => Some(*length),
            _ => None,
        })
        .fold(0.0_f64, f64::max);

    if longest_straight <= 0.0 {
        return 0.0;
    }

    // Iteratively compute achievable speed over the longest straight
    let dt = 0.1;
    let mut v = 10.0;
    let mut dist = 0.0;
    while dist < longest_straight && v < terminal {
        v = vehicle.acceleration_limited_speed(v, dt);
        dist += v * dt;
    }
    v.min(terminal)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- TrackSegment ---

    #[test]
    fn straight_length() {
        let s = TrackSegment::Straight { length: 300.0 };
        assert!((s.length() - 300.0).abs() < 1e-9);
    }

    #[test]
    fn corner_arc_length() {
        // 90° corner at 50 m radius -> arc = pi/2 * 50
        let s = TrackSegment::Corner {
            radius: 50.0,
            angle_rad: std::f64::consts::FRAC_PI_2,
            direction: 1,
        };
        let expected = 50.0 * std::f64::consts::FRAC_PI_2;
        assert!((s.length() - expected).abs() < 1e-9);
    }

    // --- Track ---

    #[test]
    fn track_new_is_empty() {
        let t = Track::new();
        assert!(t.segments.is_empty());
    }

    #[test]
    fn track_add_straight() {
        let mut t = Track::new();
        t.add_straight(500.0);
        assert_eq!(t.segments.len(), 1);
        assert!((t.total_length() - 500.0).abs() < 1e-9);
    }

    #[test]
    fn track_add_corner() {
        let mut t = Track::new();
        t.add_corner(100.0, 90.0, true);
        assert_eq!(t.segments.len(), 1);
        let expected = 100.0 * std::f64::consts::FRAC_PI_2;
        assert!((t.total_length() - expected).abs() < 1e-9);
    }

    #[test]
    fn track_total_length_mixed() {
        let mut t = Track::new();
        t.add_straight(200.0);
        t.add_corner(50.0, 90.0, false);
        t.add_straight(150.0);
        let arc = 50.0 * std::f64::consts::FRAC_PI_2;
        assert!((t.total_length() - (200.0 + arc + 150.0)).abs() < 1e-6);
    }

    #[test]
    fn track_min_corner_radius_no_corners() {
        let mut t = Track::new();
        t.add_straight(100.0);
        assert_eq!(t.min_corner_radius(), f64::INFINITY);
    }

    #[test]
    fn track_min_corner_radius_multiple() {
        let mut t = Track::new();
        t.add_corner(80.0, 45.0, true);
        t.add_corner(30.0, 90.0, false);
        t.add_corner(120.0, 30.0, true);
        assert!((t.min_corner_radius() - 30.0).abs() < 1e-9);
    }

    // --- VehicleDynamics ---

    #[test]
    fn formula_car_creates() {
        let v = VehicleDynamics::formula_car();
        assert!(v.mass > 0.0);
        assert!(v.max_lateral_g > 0.0);
    }

    #[test]
    fn sports_car_creates() {
        let v = VehicleDynamics::sports_car();
        assert!(v.mass > v.max_engine_force / 10.0);
    }

    #[test]
    fn max_speed_corner_scales_with_radius() {
        let v = VehicleDynamics::formula_car();
        let s1 = v.max_speed_in_corner(50.0);
        let s2 = v.max_speed_in_corner(200.0);
        assert!(s2 > s1);
    }

    #[test]
    fn max_speed_corner_formula() {
        let v = VehicleDynamics::formula_car();
        let expected = (v.max_lateral_g * G * 100.0f64).sqrt();
        let got = v.max_speed_in_corner(100.0);
        assert!((got - expected).abs() < 1e-9);
    }

    #[test]
    fn acceleration_limited_speed_increases() {
        let v = VehicleDynamics::formula_car();
        let v0 = 10.0;
        let v1 = v.acceleration_limited_speed(v0, 0.1);
        assert!(v1 > v0);
    }

    #[test]
    fn braking_distance_zero_when_already_slower() {
        let v = VehicleDynamics::formula_car();
        let d = v.braking_distance(20.0, 30.0);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn braking_distance_positive() {
        let v = VehicleDynamics::formula_car();
        let d = v.braking_distance(60.0, 20.0);
        assert!(d > 0.0);
    }

    #[test]
    fn braking_distance_kinematics() {
        let v = VehicleDynamics::formula_car();
        let v0 = 50.0;
        let vt = 0.0;
        let a = v.max_brake_force / v.mass;
        let expected = (v0 * v0 - vt * vt) / (2.0 * a);
        let got = v.braking_distance(v0, vt);
        assert!((got - expected).abs() < 1e-6);
    }

    // --- LapSimulator ---

    #[test]
    fn lap_sim_empty_track() {
        let t = Track::new();
        let sim = LapSimulator::new(t, VehicleDynamics::formula_car());
        let r = sim.simulate();
        assert_eq!(r.total_time, 0.0);
        assert!(r.sector_times.is_empty());
    }

    #[test]
    fn lap_sim_single_straight() {
        let mut t = Track::new();
        t.add_straight(100.0);
        let sim = LapSimulator::new(t, VehicleDynamics::formula_car());
        let r = sim.simulate();
        assert!(r.total_time > 0.0);
        assert_eq!(r.sector_times.len(), 1);
    }

    #[test]
    fn lap_sim_max_speed_positive() {
        let mut t = Track::new();
        t.add_straight(1000.0);
        t.add_corner(100.0, 90.0, true);
        let sim = LapSimulator::new(t, VehicleDynamics::sports_car());
        let r = sim.simulate();
        assert!(r.max_speed > 0.0);
    }

    #[test]
    fn lap_sim_avg_speed_consistent() {
        let mut t = Track::new();
        t.add_straight(500.0);
        t.add_corner(80.0, 90.0, false);
        let sim = LapSimulator::new(t, VehicleDynamics::formula_car());
        let r = sim.simulate();
        let total_len = sim.track.total_length();
        let expected_avg = total_len / r.total_time;
        assert!((r.avg_speed - expected_avg).abs() < 1e-6);
    }

    #[test]
    fn sector_times_matches_simulate() {
        let mut t = Track::new();
        t.add_straight(300.0);
        t.add_corner(60.0, 45.0, true);
        t.add_straight(200.0);
        let sim = LapSimulator::new(t, VehicleDynamics::formula_car());
        let st = sim.sector_times();
        let r = sim.simulate();
        assert_eq!(st.len(), r.sector_times.len());
        for (a, b) in st.iter().zip(r.sector_times.iter()) {
            assert!((a - b).abs() < 1e-9);
        }
    }

    // --- RaceLineOptimizer ---

    #[test]
    fn geometric_line_empty_track() {
        let opt = RaceLineOptimizer::new(Track::new());
        let pts = opt.geometric_racing_line();
        assert!(pts.is_empty());
    }

    #[test]
    fn geometric_line_single_straight() {
        let mut t = Track::new();
        t.add_straight(100.0);
        let opt = RaceLineOptimizer::new(t);
        let pts = opt.geometric_racing_line();
        assert_eq!(pts.len(), 1);
    }

    #[test]
    fn smoothed_line_increases_points() {
        let mut t = Track::new();
        t.add_straight(100.0);
        t.add_corner(50.0, 90.0, true);
        t.add_straight(100.0);
        let opt = RaceLineOptimizer::new(t);
        let raw = opt.geometric_racing_line();
        let smooth = opt.smoothed_line(&raw, 3);
        assert!(smooth.len() >= raw.len());
    }

    #[test]
    fn smoothed_line_zero_iterations() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let t = Track::new();
        let opt = RaceLineOptimizer::new(t);
        let out = opt.smoothed_line(&pts, 0);
        assert_eq!(out.len(), pts.len());
    }

    // --- fuel_consumption ---

    #[test]
    fn fuel_zero_speed() {
        let f = fuel_consumption(1000.0, 0.0, 0.3, 100.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn fuel_positive() {
        let f = fuel_consumption(5000.0, 50.0, 0.25, 150.0);
        assert!(f > 0.0);
    }

    #[test]
    fn fuel_scales_with_distance() {
        let f1 = fuel_consumption(1000.0, 50.0, 0.25, 100.0);
        let f2 = fuel_consumption(2000.0, 50.0, 0.25, 100.0);
        assert!((f2 - 2.0 * f1).abs() < 1e-9);
    }

    // --- tire_wear_model ---

    #[test]
    fn tire_wear_zero_lateral() {
        let w = tire_wear_model(0.0, 50.0, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn tire_wear_positive() {
        let w = tire_wear_model(2.0, 60.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn tire_wear_scales_with_compound() {
        let w1 = tire_wear_model(1.5, 40.0, 1.0);
        let w2 = tire_wear_model(1.5, 40.0, 2.0);
        assert!((w2 - 2.0 * w1).abs() < 1e-9);
    }

    // --- classify_segment ---

    #[test]
    fn classify_long_straight() {
        let s = TrackSegment::Straight { length: 500.0 };
        assert_eq!(classify_segment(&s), SegmentClass::LongStraight);
    }

    #[test]
    fn classify_short_straight() {
        let s = TrackSegment::Straight { length: 100.0 };
        assert_eq!(classify_segment(&s), SegmentClass::ShortStraight);
    }

    #[test]
    fn classify_high_speed_corner() {
        let s = TrackSegment::Corner {
            radius: 150.0,
            angle_rad: 1.0,
            direction: 1,
        };
        assert_eq!(classify_segment(&s), SegmentClass::HighSpeedCorner);
    }

    #[test]
    fn classify_medium_speed_corner() {
        let s = TrackSegment::Corner {
            radius: 75.0,
            angle_rad: 1.0,
            direction: -1,
        };
        assert_eq!(classify_segment(&s), SegmentClass::MediumSpeedCorner);
    }

    #[test]
    fn classify_low_speed_corner() {
        let s = TrackSegment::Corner {
            radius: 30.0,
            angle_rad: 1.0,
            direction: 1,
        };
        assert_eq!(classify_segment(&s), SegmentClass::LowSpeedCorner);
    }

    #[test]
    fn classify_hairpin() {
        let s = TrackSegment::Corner {
            radius: 15.0,
            angle_rad: std::f64::consts::PI, // 180 degrees
            direction: 1,
        };
        assert_eq!(classify_segment(&s), SegmentClass::Hairpin);
    }

    #[test]
    fn classify_track_mixed() {
        let mut t = Track::new();
        t.add_straight(500.0);
        t.add_corner(30.0, 90.0, true);
        t.add_straight(100.0);
        let classes = classify_track(&t);
        assert_eq!(classes.len(), 3);
        assert_eq!(classes[0], SegmentClass::LongStraight);
        assert_eq!(classes[1], SegmentClass::LowSpeedCorner);
        assert_eq!(classes[2], SegmentClass::ShortStraight);
    }

    // --- apex detection ---

    #[test]
    fn detect_apexes_empty_track() {
        let t = Track::new();
        let apexes = detect_apexes(&t);
        assert!(apexes.is_empty());
    }

    #[test]
    fn detect_apexes_single_corner() {
        let mut t = Track::new();
        t.add_straight(100.0);
        t.add_corner(50.0, 90.0, true);
        t.add_straight(100.0);
        let apexes = detect_apexes(&t);
        assert_eq!(apexes.len(), 1);
        assert!(apexes[0].is_geometric);
        assert!((apexes[0].curvature - 1.0 / 50.0).abs() < 1e-9);
        assert_eq!(apexes[0].segment_index, 1);
    }

    #[test]
    fn detect_apexes_multiple_corners() {
        let mut t = Track::new();
        t.add_corner(80.0, 45.0, true);
        t.add_straight(200.0);
        t.add_corner(30.0, 90.0, false);
        t.add_straight(100.0);
        t.add_corner(120.0, 30.0, true);
        let apexes = detect_apexes(&t);
        assert_eq!(apexes.len(), 3);
    }

    #[test]
    fn detect_apexes_no_corners() {
        let mut t = Track::new();
        t.add_straight(100.0);
        t.add_straight(200.0);
        let apexes = detect_apexes(&t);
        assert!(apexes.is_empty());
    }

    // --- braking points ---

    #[test]
    fn braking_points_empty_track() {
        let t = Track::new();
        let v = VehicleDynamics::formula_car();
        let bp = compute_braking_points(&t, &v);
        assert!(bp.is_empty());
    }

    #[test]
    fn braking_points_straight_only() {
        let mut t = Track::new();
        t.add_straight(1000.0);
        let v = VehicleDynamics::formula_car();
        let bp = compute_braking_points(&t, &v);
        assert!(bp.is_empty());
    }

    #[test]
    fn braking_points_before_corner() {
        let mut t = Track::new();
        t.add_straight(500.0);
        t.add_corner(50.0, 90.0, true);
        let v = VehicleDynamics::formula_car();
        let bp = compute_braking_points(&t, &v);
        assert_eq!(bp.len(), 1);
        assert!(bp[0].braking_distance >= 0.0);
        assert!(bp[0].target_speed > 0.0);
        assert!(bp[0].speed_at_braking >= 0.0);
    }

    #[test]
    fn braking_points_target_speed_matches_corner() {
        let mut t = Track::new();
        t.add_straight(300.0);
        t.add_corner(100.0, 90.0, false);
        let v = VehicleDynamics::formula_car();
        let bp = compute_braking_points(&t, &v);
        let expected = v.max_speed_in_corner(100.0);
        assert!((bp[0].target_speed - expected).abs() < 1e-6);
    }

    // --- gaussian_smooth ---

    #[test]
    fn gaussian_smooth_preserves_length() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0], [4.0, 0.0]];
        let smoothed = gaussian_smooth(&pts, 1.0);
        assert_eq!(smoothed.len(), pts.len());
    }

    #[test]
    fn gaussian_smooth_reduces_variation() {
        let pts = vec![[0.0, 0.0], [1.0, 10.0], [2.0, 0.0], [3.0, 10.0], [4.0, 0.0]];
        let smoothed = gaussian_smooth(&pts, 1.5);
        // After smoothing, the y-range should be reduced
        let y_range_orig: f64 = pts.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max)
            - pts.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        let y_range_smooth: f64 = smoothed
            .iter()
            .map(|p| p[1])
            .fold(f64::NEG_INFINITY, f64::max)
            - smoothed.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        assert!(y_range_smooth < y_range_orig);
    }

    #[test]
    fn gaussian_smooth_too_few_points() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0]];
        let smoothed = gaussian_smooth(&pts, 1.0);
        assert_eq!(smoothed.len(), 2);
    }

    #[test]
    fn gaussian_smooth_zero_sigma() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]];
        let smoothed = gaussian_smooth(&pts, 0.0);
        for (a, b) in pts.iter().zip(smoothed.iter()) {
            assert!((a[0] - b[0]).abs() < 1e-9);
            assert!((a[1] - b[1]).abs() < 1e-9);
        }
    }

    // --- polyline_curvature ---

    #[test]
    fn curvature_straight_line() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]];
        let curv = polyline_curvature(&pts);
        for c in &curv {
            assert!(*c < 1e-10, "straight line curvature should be ~0");
        }
    }

    #[test]
    fn curvature_right_angle() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        let curv = polyline_curvature(&pts);
        assert!(curv[1] > 0.0, "right angle should have positive curvature");
    }

    #[test]
    fn curvature_too_few_points() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0]];
        let curv = polyline_curvature(&pts);
        assert_eq!(curv.len(), 2);
        assert!(curv[0] < 1e-10);
    }

    #[test]
    fn total_curvature_straight() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        assert!(total_curvature(&pts) < 1e-10);
    }

    // --- path_length ---

    #[test]
    fn path_length_straight() {
        let pts = vec![[0.0, 0.0], [3.0, 4.0]];
        assert!((path_length(&pts) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn path_length_empty() {
        let pts: Vec<[f64; 2]> = vec![];
        assert_eq!(path_length(&pts), 0.0);
    }

    #[test]
    fn path_length_multi_segment() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        assert!((path_length(&pts) - 2.0).abs() < 1e-9);
    }

    // --- optimize_racing_line ---

    #[test]
    fn optimize_reduces_curvature() {
        let pts = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.5, 0.5],
            [2.0, 1.0],
            [2.5, 0.5],
            [3.0, 0.0],
            [4.0, 0.0],
        ];
        let orig_curv = total_curvature(&pts);
        let optimized = optimize_racing_line(&pts, 0.3, 50);
        let opt_curv = total_curvature(&optimized);
        assert!(
            opt_curv <= orig_curv + 1e-6,
            "optimized curvature ({opt_curv}) should be <= original ({orig_curv})"
        );
    }

    #[test]
    fn optimize_preserves_endpoints() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0], [4.0, 0.0]];
        let opt = optimize_racing_line(&pts, 0.3, 20);
        assert!((opt[0][0] - pts[0][0]).abs() < 1e-9);
        assert!((opt[0][1] - pts[0][1]).abs() < 1e-9);
        let last = pts.len() - 1;
        assert!((opt[last][0] - pts[last][0]).abs() < 1e-9);
        assert!((opt[last][1] - pts[last][1]).abs() < 1e-9);
    }

    #[test]
    fn optimize_too_few_points() {
        let pts = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]];
        let opt = optimize_racing_line(&pts, 0.3, 10);
        assert_eq!(opt.len(), 3);
    }

    // --- track analysis ---

    #[test]
    fn count_corners_mixed() {
        let mut t = Track::new();
        t.add_straight(100.0);
        t.add_corner(50.0, 90.0, true);
        t.add_straight(200.0);
        t.add_corner(30.0, 45.0, false);
        assert_eq!(count_corners(&t), 2);
        assert_eq!(count_straights(&t), 2);
    }

    #[test]
    fn max_corner_radius_test() {
        let mut t = Track::new();
        t.add_corner(50.0, 90.0, true);
        t.add_corner(120.0, 45.0, false);
        assert!((max_corner_radius(&t) - 120.0).abs() < 1e-9);
    }

    #[test]
    fn avg_corner_radius_test() {
        let mut t = Track::new();
        t.add_corner(50.0, 90.0, true);
        t.add_corner(100.0, 45.0, false);
        assert!((avg_corner_radius(&t) - 75.0).abs() < 1e-9);
    }

    #[test]
    fn avg_corner_radius_no_corners() {
        let mut t = Track::new();
        t.add_straight(100.0);
        assert_eq!(avg_corner_radius(&t), 0.0);
    }

    #[test]
    fn corner_percentage_all_straight() {
        let mut t = Track::new();
        t.add_straight(100.0);
        assert!((corner_percentage(&t)).abs() < 1e-6);
    }

    #[test]
    fn corner_percentage_empty() {
        let t = Track::new();
        assert_eq!(corner_percentage(&t), 0.0);
    }

    #[test]
    fn corner_percentage_mixed() {
        let mut t = Track::new();
        t.add_straight(100.0);
        t.add_corner(50.0, 90.0, true);
        let pct = corner_percentage(&t);
        assert!(pct > 0.0 && pct < 100.0);
    }

    #[test]
    fn estimated_top_speed_positive() {
        let mut t = Track::new();
        t.add_straight(1000.0);
        let v = VehicleDynamics::formula_car();
        let top = estimated_top_speed(&t, &v);
        assert!(top > 10.0);
    }

    #[test]
    fn estimated_top_speed_no_straight() {
        let mut t = Track::new();
        t.add_corner(50.0, 90.0, true);
        let v = VehicleDynamics::formula_car();
        let top = estimated_top_speed(&t, &v);
        assert_eq!(top, 0.0);
    }

    #[test]
    fn estimated_top_speed_longer_straight_is_faster() {
        let v = VehicleDynamics::formula_car();
        let mut t1 = Track::new();
        t1.add_straight(200.0);
        let mut t2 = Track::new();
        t2.add_straight(1000.0);
        assert!(estimated_top_speed(&t2, &v) >= estimated_top_speed(&t1, &v));
    }
}
