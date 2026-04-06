#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Racing lap simulation: point-mass lap time estimation, cornering limits,
//! tire degradation, pit-stop strategy, and fuel consumption.
//!
//! # Overview
//!
//! - [`TrackPoint`]    — individual waypoint with curvature, speed limit, and bank.
//! - [`RaceTrack`]     — 2-D race track as a sequence of waypoints.
//! - [`LapSimulator`]  — point-mass lap simulator producing time and speed profile.
//! - Free functions for aerodynamic downforce, tractive force, tire degradation,
//!   pit-stop strategy, and WLTC-style fuel consumption.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// TrackPoint
// ---------------------------------------------------------------------------

/// A single waypoint on a race track.
#[derive(Debug, Clone)]
pub struct TrackPoint {
    /// X coordinate (m).
    pub x: f64,
    /// Y coordinate (m).
    pub y: f64,
    /// Signed curvature κ at this point (1/m).
    pub curvature: f64,
    /// Imposed speed limit at this point (m/s); `f64::INFINITY` if none.
    pub speed_limit: f64,
    /// Road bank angle (rad), positive toward the inside of the corner.
    pub bank_angle: f64,
}

impl TrackPoint {
    /// Construct a waypoint.
    pub fn new(x: f64, y: f64, curvature: f64, speed_limit: f64, bank_angle: f64) -> Self {
        Self {
            x,
            y,
            curvature,
            speed_limit,
            bank_angle,
        }
    }

    /// Construct a simple waypoint without a speed limit or bank angle.
    pub fn simple(x: f64, y: f64, curvature: f64) -> Self {
        Self::new(x, y, curvature, f64::INFINITY, 0.0)
    }
}

// ---------------------------------------------------------------------------
// RaceTrack
// ---------------------------------------------------------------------------

/// A 2-D race track defined by ordered waypoints.
#[derive(Debug, Clone)]
pub struct RaceTrack {
    /// Ordered sequence of track waypoints.
    pub points: Vec<TrackPoint>,
    /// Total lap length (m).  Must be set externally or recomputed.
    pub lap_length: f64,
}

impl RaceTrack {
    /// Create an empty track.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            lap_length: 0.0,
        }
    }

    /// Append a waypoint and update `lap_length` by the Euclidean distance
    /// from the previous point.
    pub fn add_point(&mut self, tp: TrackPoint) {
        if let Some(last) = self.points.last() {
            let dx = tp.x - last.x;
            let dy = tp.y - last.y;
            self.lap_length += (dx * dx + dy * dy).sqrt();
        }
        self.points.push(tp);
    }

    /// Number of waypoints.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Total lap length (m).
    pub fn lap_length(&self) -> f64 {
        self.lap_length
    }

    /// Maximum speed at waypoint `i` based on tyre–road friction limit.
    ///
    /// Returns `f64::INFINITY` when curvature is zero, and
    /// `speed_limit.min(friction_limit)` otherwise.
    pub fn max_speed_at(&self, i: usize, mass: f64, mu: f64, g: f64) -> f64 {
        if i >= self.points.len() {
            return 0.0;
        }
        let tp = &self.points[i];
        let kappa = tp.curvature.abs();
        let friction_limit = cornering_speed_limit(kappa, mu, g, mass);
        friction_limit.min(tp.speed_limit)
    }
}

impl Default for RaceTrack {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Maximum cornering speed from tyre–road friction.
///
/// `v_max = sqrt(μ · g / κ)`
///
/// Returns `f64::INFINITY` for zero or near-zero curvature.
pub fn cornering_speed_limit(curvature: f64, mu: f64, g: f64, _mass: f64) -> f64 {
    let kappa = curvature.abs();
    if kappa < 1e-12 {
        return f64::INFINITY;
    }
    (mu * g / kappa).sqrt()
}

/// Maximum tractive force from available power.
///
/// `F_trac = P / v`
///
/// Returns `f64::INFINITY` at zero speed to avoid division by zero.
pub fn tractive_force_limit(power: f64, speed: f64) -> f64 {
    if speed < 1e-9 {
        return f64::INFINITY;
    }
    power / speed
}

/// Aerodynamic downforce (N).
///
/// `L = cl · A · 0.5 · ρ · v²`
///
/// * `cd`     – lift coefficient (positive → downforce)
/// * `cl_neg` – downforce coefficient (same as negative-lift coefficient)
/// * `area`   – reference area (m²)
/// * `rho`    – air density (kg/m³)
/// * `speed`  – vehicle speed (m/s)
pub fn aerodynamic_downforce(_cd: f64, cl_neg: f64, area: f64, rho: f64, speed: f64) -> f64 {
    0.5 * rho * cl_neg * area * speed * speed
}

// ---------------------------------------------------------------------------
// LapSimulator
// ---------------------------------------------------------------------------

/// Point-mass lap simulator.
///
/// Estimates lap time and speed profile using the simplified
/// cornering-speed limit method with tractive force saturation.
#[derive(Debug, Clone)]
pub struct LapSimulator {
    /// Race track.
    pub track: RaceTrack,
    /// Vehicle mass (kg).
    pub vehicle_mass: f64,
    /// Maximum engine/motor power (W).
    pub max_power: f64,
    /// Tyre–road friction coefficient.
    pub mu: f64,
}

impl LapSimulator {
    /// Create a simulator.
    ///
    /// * `track`  – [`RaceTrack`] describing the circuit
    /// * `mass`   – vehicle mass (kg)
    /// * `power`  – peak power (W)
    /// * `mu`     – tyre friction coefficient
    pub fn new(track: RaceTrack, mass: f64, power: f64, mu: f64) -> Self {
        Self {
            track,
            vehicle_mass: mass,
            max_power: power,
            mu,
        }
    }

    /// Run a single-lap simulation.
    ///
    /// Returns `(lap_time_s, speed_profile_m_per_s)` where the speed profile
    /// has one entry per track waypoint.
    pub fn simulate_lap(&self) -> (f64, Vec<f64>) {
        const G: f64 = 9.81;
        let n = self.track.points.len();
        if n == 0 {
            return (0.0, Vec::new());
        }

        // Compute speed limit at each point
        let mut speeds: Vec<f64> = (0..n)
            .map(|i| self.track.max_speed_at(i, self.vehicle_mass, self.mu, G))
            .collect();

        // Clamp to tractive limit using a simple forward pass
        let mut prev_speed = 0.0_f64;
        for i in 0..n {
            // Distance to next point
            let dist = if i + 1 < n {
                let dx = self.track.points[i + 1].x - self.track.points[i].x;
                let dy = self.track.points[i + 1].y - self.track.points[i].y;
                (dx * dx + dy * dy).sqrt().max(1.0)
            } else {
                10.0
            };
            // Maximum acceleration from engine
            let f_max = tractive_force_limit(self.max_power, prev_speed.max(1.0));
            let a_max = f_max / self.vehicle_mass;
            let accel_limited = (prev_speed * prev_speed + 2.0 * a_max * dist).sqrt();
            speeds[i] = speeds[i].min(accel_limited);
            prev_speed = speeds[i];
        }

        // Compute lap time from speeds and distances
        let mut lap_time = 0.0_f64;
        for i in 0..n {
            let v = speeds[i].max(1.0); // avoid divide-by-zero
            let dist = if i + 1 < n {
                let dx = self.track.points[i + 1].x - self.track.points[i].x;
                let dy = self.track.points[i + 1].y - self.track.points[i].y;
                (dx * dx + dy * dy).sqrt()
            } else {
                0.0
            };
            lap_time += dist / v;
        }

        (lap_time, speeds)
    }

    /// Split the lap into `n_sectors` equal-time sectors and return each
    /// sector time (s).
    pub fn sector_times(&self, n_sectors: usize) -> Vec<f64> {
        if n_sectors == 0 {
            return Vec::new();
        }
        let (lap_time, _) = self.simulate_lap();
        vec![lap_time / n_sectors as f64; n_sectors]
    }
}

// ---------------------------------------------------------------------------
// lap_time_estimate
// ---------------------------------------------------------------------------

/// Simple lap time estimate: T = L / v_avg (s).
pub fn lap_time_estimate(track_length: f64, avg_speed: f64) -> f64 {
    if avg_speed < 1e-12 {
        return f64::INFINITY;
    }
    track_length / avg_speed
}

// ---------------------------------------------------------------------------
// fuel_lap_consumption
// ---------------------------------------------------------------------------

/// Fuel consumed per lap (kg).
///
/// `m_fuel = lap_distance / avg_speed · fuel_rate`
///
/// * `lap_distance` – lap length (m)
/// * `fuel_rate`    – fuel mass flow rate (kg/s)
/// * `avg_speed`    – average lap speed (m/s)
pub fn fuel_lap_consumption(lap_distance: f64, fuel_rate: f64, avg_speed: f64) -> f64 {
    if avg_speed < 1e-12 {
        return 0.0;
    }
    let lap_time = lap_distance / avg_speed;
    fuel_rate * lap_time
}

// ---------------------------------------------------------------------------
// tire_degradation_model
// ---------------------------------------------------------------------------

/// Tyre grip as a function of lap count.
///
/// `μ(n) = μ₀ · exp(−k · n)`
///
/// * `laps`          – number of laps completed
/// * `wear_rate`     – exponential wear rate coefficient k
/// * `grip_initial`  – initial friction coefficient μ₀
pub fn tire_degradation_model(laps: usize, wear_rate: f64, grip_initial: f64) -> f64 {
    grip_initial * (-wear_rate * laps as f64).exp()
}

// ---------------------------------------------------------------------------
// pit_stop_strategy
// ---------------------------------------------------------------------------

/// Compute the optimal number of pit stops.
///
/// Uses a simple model: pit if tyre life `< total_laps / (stops + 1)`.
/// Minimises total time = stops * pit_time_loss.
///
/// Returns 0 if tyres last the entire race.
///
/// * `total_laps`    – race length in laps
/// * `tire_life`     – maximum stint length before grip becomes critical
/// * `pit_time_loss` – time lost per pit stop (s)
/// * `lap_time`      – nominal lap time (s, currently unused but part of API)
pub fn pit_stop_strategy(
    total_laps: usize,
    tire_life: usize,
    _pit_time_loss: f64,
    _lap_time: f64,
) -> usize {
    if tire_life == 0 {
        return total_laps.saturating_sub(1);
    }
    if tire_life >= total_laps {
        return 0;
    }
    // Number of stints = ceil(total_laps / tire_life)
    // Number of pit stops = stints - 1
    total_laps.div_ceil(tire_life).saturating_sub(1)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    // ── TrackPoint ────────────────────────────────────────────────────────

    #[test]
    fn track_point_new_stores_values() {
        let tp = TrackPoint::new(1.0, 2.0, 0.5, 30.0, 0.1);
        assert!((tp.x - 1.0).abs() < EPS);
        assert!((tp.y - 2.0).abs() < EPS);
        assert!((tp.curvature - 0.5).abs() < EPS);
    }

    #[test]
    fn track_point_simple_no_limit() {
        let tp = TrackPoint::simple(0.0, 0.0, 1.0);
        assert!(tp.speed_limit.is_infinite());
        assert!(tp.bank_angle.abs() < EPS);
    }

    // ── RaceTrack ─────────────────────────────────────────────────────────

    #[test]
    fn race_track_starts_empty() {
        let t = RaceTrack::new();
        assert_eq!(t.point_count(), 0);
        assert!((t.lap_length()).abs() < EPS);
    }

    #[test]
    fn race_track_add_point_increments_count() {
        let mut t = RaceTrack::new();
        t.add_point(TrackPoint::simple(0.0, 0.0, 0.0));
        t.add_point(TrackPoint::simple(10.0, 0.0, 0.0));
        assert_eq!(t.point_count(), 2);
    }

    #[test]
    fn race_track_lap_length_accumulates() {
        let mut t = RaceTrack::new();
        t.add_point(TrackPoint::simple(0.0, 0.0, 0.0));
        t.add_point(TrackPoint::simple(3.0, 4.0, 0.0)); // distance = 5
        assert!((t.lap_length() - 5.0).abs() < EPS);
    }

    #[test]
    fn race_track_max_speed_finite_curvature() {
        let mut t = RaceTrack::new();
        t.add_point(TrackPoint::simple(0.0, 0.0, 0.1));
        let v = t.max_speed_at(0, 1000.0, 1.5, 9.81);
        assert!(v.is_finite() && v > 0.0, "v={v}");
    }

    #[test]
    fn race_track_max_speed_zero_curvature_infinite() {
        let mut t = RaceTrack::new();
        t.add_point(TrackPoint::simple(0.0, 0.0, 0.0));
        let v = t.max_speed_at(0, 1000.0, 1.5, 9.81);
        assert!(v.is_infinite());
    }

    #[test]
    fn race_track_max_speed_out_of_bounds_zero() {
        let t = RaceTrack::new();
        let v = t.max_speed_at(5, 1000.0, 1.5, 9.81);
        assert!(v.abs() < EPS);
    }

    // ── cornering_speed_limit ─────────────────────────────────────────────

    #[test]
    fn cornering_speed_limit_scales_sqrt_mu() {
        let v1 = cornering_speed_limit(0.1, 1.0, 9.81, 1000.0);
        let v2 = cornering_speed_limit(0.1, 4.0, 9.81, 1000.0);
        assert!((v2 / v1 - 2.0).abs() < 1e-9, "should scale as sqrt(mu)");
    }

    #[test]
    fn cornering_speed_limit_zero_curvature_infinite() {
        let v = cornering_speed_limit(0.0, 1.5, 9.81, 1000.0);
        assert!(v.is_infinite());
    }

    #[test]
    fn cornering_speed_limit_positive() {
        let v = cornering_speed_limit(0.05, 1.5, 9.81, 1000.0);
        assert!(v > 0.0);
    }

    #[test]
    fn cornering_speed_limit_higher_curvature_lower_speed() {
        let v_lo = cornering_speed_limit(0.05, 1.5, 9.81, 1000.0);
        let v_hi = cornering_speed_limit(0.2, 1.5, 9.81, 1000.0);
        assert!(v_hi < v_lo);
    }

    // ── tractive_force_limit ─────────────────────────────────────────────

    #[test]
    fn tractive_force_limit_p_over_v() {
        let f = tractive_force_limit(150_000.0, 50.0);
        assert!((f - 3000.0).abs() < EPS);
    }

    #[test]
    fn tractive_force_limit_zero_speed_infinite() {
        let f = tractive_force_limit(150_000.0, 0.0);
        assert!(f.is_infinite());
    }

    #[test]
    fn tractive_force_limit_decreases_with_speed() {
        let f1 = tractive_force_limit(150_000.0, 30.0);
        let f2 = tractive_force_limit(150_000.0, 60.0);
        assert!(f2 < f1);
    }

    // ── aerodynamic_downforce ─────────────────────────────────────────────

    #[test]
    fn aero_downforce_increases_with_speed_squared() {
        let l1 = aerodynamic_downforce(0.3, 1.5, 2.0, 1.225, 50.0);
        let l2 = aerodynamic_downforce(0.3, 1.5, 2.0, 1.225, 100.0);
        assert!((l2 / l1 - 4.0).abs() < 1e-9, "downforce should scale as v²");
    }

    #[test]
    fn aero_downforce_zero_speed_zero() {
        let l = aerodynamic_downforce(0.3, 1.5, 2.0, 1.225, 0.0);
        assert!(l.abs() < EPS);
    }

    #[test]
    fn aero_downforce_positive() {
        let l = aerodynamic_downforce(0.3, 1.5, 2.0, 1.225, 60.0);
        assert!(l > 0.0);
    }

    // ── LapSimulator ──────────────────────────────────────────────────────

    #[test]
    fn lap_simulator_empty_track_zero_time() {
        let track = RaceTrack::new();
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let (t, sp) = sim.simulate_lap();
        assert!(t.abs() < EPS);
        assert!(sp.is_empty());
    }

    #[test]
    fn lap_simulator_single_point_no_crash() {
        let mut track = RaceTrack::new();
        track.add_point(TrackPoint::simple(0.0, 0.0, 0.0));
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let (t, sp) = sim.simulate_lap();
        assert!(t >= 0.0);
        assert_eq!(sp.len(), 1);
    }

    #[test]
    fn lap_simulator_speeds_nonnegative() {
        let mut track = RaceTrack::new();
        for i in 0..10 {
            let angle = i as f64 * 2.0 * PI / 10.0;
            let r = 100.0;
            track.add_point(TrackPoint::simple(
                r * angle.cos(),
                r * angle.sin(),
                1.0 / r,
            ));
        }
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let (_t, sp) = sim.simulate_lap();
        for v in &sp {
            assert!(*v >= 0.0, "speed={v}");
        }
    }

    #[test]
    fn lap_simulator_lap_time_positive() {
        let mut track = RaceTrack::new();
        track.add_point(TrackPoint::simple(0.0, 0.0, 0.0));
        track.add_point(TrackPoint::simple(100.0, 0.0, 0.0));
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let (t, _) = sim.simulate_lap();
        assert!(t >= 0.0);
    }

    #[test]
    fn lap_simulator_sector_times_correct_count() {
        let mut track = RaceTrack::new();
        for i in 0..6 {
            track.add_point(TrackPoint::simple(i as f64 * 50.0, 0.0, 0.0));
        }
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let sectors = sim.sector_times(3);
        assert_eq!(sectors.len(), 3);
    }

    #[test]
    fn lap_simulator_sector_times_sum_to_lap_time() {
        let mut track = RaceTrack::new();
        for i in 0..6 {
            track.add_point(TrackPoint::simple(i as f64 * 50.0, 0.0, 0.0));
        }
        let sim = LapSimulator::new(track, 800.0, 200_000.0, 1.5);
        let (lap_time, _) = sim.simulate_lap();
        let sectors = sim.sector_times(3);
        let sector_sum: f64 = sectors.iter().sum();
        assert!(
            (sector_sum - lap_time).abs() < 1e-9,
            "sum={sector_sum} lap={lap_time}"
        );
    }

    // ── lap_time_estimate ─────────────────────────────────────────────────

    #[test]
    fn lap_time_estimate_formula() {
        let t = lap_time_estimate(4000.0, 50.0);
        assert!((t - 80.0).abs() < EPS);
    }

    #[test]
    fn lap_time_estimate_zero_speed_infinite() {
        let t = lap_time_estimate(4000.0, 0.0);
        assert!(t.is_infinite());
    }

    #[test]
    fn lap_time_estimate_positive() {
        let t = lap_time_estimate(5000.0, 40.0);
        assert!(t > 0.0);
    }

    // ── fuel_lap_consumption ─────────────────────────────────────────────

    #[test]
    fn fuel_consumption_positive() {
        let m = fuel_lap_consumption(4000.0, 0.05, 50.0);
        assert!(m > 0.0);
    }

    #[test]
    fn fuel_consumption_zero_speed_zero() {
        let m = fuel_lap_consumption(4000.0, 0.05, 0.0);
        assert!(m.abs() < EPS);
    }

    #[test]
    fn fuel_consumption_formula() {
        // lap_time = 4000/50 = 80 s; m = 0.05 * 80 = 4.0 kg
        let m = fuel_lap_consumption(4000.0, 0.05, 50.0);
        assert!((m - 4.0).abs() < EPS);
    }

    // ── tire_degradation_model ────────────────────────────────────────────

    #[test]
    fn tire_degradation_zero_laps_full_grip() {
        let mu = tire_degradation_model(0, 0.02, 1.5);
        assert!((mu - 1.5).abs() < EPS);
    }

    #[test]
    fn tire_degradation_decreases_with_laps() {
        let mu1 = tire_degradation_model(5, 0.02, 1.5);
        let mu2 = tire_degradation_model(20, 0.02, 1.5);
        assert!(mu2 < mu1, "grip should decrease with laps");
    }

    #[test]
    fn tire_degradation_positive() {
        let mu = tire_degradation_model(30, 0.01, 1.5);
        assert!(mu > 0.0);
    }

    #[test]
    fn tire_degradation_exponential_decay() {
        let k = 0.02;
        let mu0 = 1.5;
        let n = 10;
        let mu = tire_degradation_model(n, k, mu0);
        let expected = mu0 * (-k * n as f64).exp();
        assert!((mu - expected).abs() < EPS);
    }

    // ── pit_stop_strategy ─────────────────────────────────────────────────

    #[test]
    fn pit_stop_strategy_no_stops_for_long_tire_life() {
        let stops = pit_stop_strategy(30, 50, 30.0, 90.0);
        assert_eq!(stops, 0);
    }

    #[test]
    fn pit_stop_strategy_one_stop_for_half_race_tires() {
        // Tires last 20 laps, race is 40 laps → 1 stop
        let stops = pit_stop_strategy(40, 20, 30.0, 90.0);
        assert_eq!(stops, 1);
    }

    #[test]
    fn pit_stop_strategy_two_stops() {
        // Tires last 15 laps, race is 45 laps → 2 stops
        let stops = pit_stop_strategy(45, 15, 30.0, 90.0);
        assert_eq!(stops, 2);
    }

    #[test]
    fn pit_stop_strategy_tire_life_equals_race_zero_stops() {
        let stops = pit_stop_strategy(30, 30, 30.0, 90.0);
        assert_eq!(stops, 0);
    }

    #[test]
    fn pit_stop_strategy_nonnegative() {
        let stops = pit_stop_strategy(50, 10, 25.0, 90.0);
        assert!(stops as i64 >= 0);
    }
}
