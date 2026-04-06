// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Train and rail vehicle dynamics.
//!
//! Provides models for multi-car train consists, track geometry, wheelset
//! mechanics, Davis resistance, and signaling headway calculations.

#![allow(dead_code)]

// ── TrainCar ──────────────────────────────────────────────────────────────────

/// A single passenger or freight car in a train consist.
#[derive(Debug, Clone)]
pub struct TrainCar {
    /// Mass of the car in kilograms.
    pub mass: f64,
    /// Overall length of the car in metres.
    pub length: f64,
    /// Distance between front and rear bogies (wheel centres) in metres.
    pub wheel_base: f64,
    /// Position of the car centre along the track in metres.
    pub position: f64,
    /// Longitudinal velocity of the car in m/s.
    pub velocity: f64,
    /// Longitudinal acceleration of the car in m/s².
    pub acceleration: f64,
}

impl TrainCar {
    /// Create a new `TrainCar` with the given physical parameters.
    ///
    /// * `mass`       – Mass in kg.
    /// * `length`     – Overall length in m.
    /// * `wheel_base` – Bogie centre distance in m.
    pub fn new(mass: f64, length: f64, wheel_base: f64) -> Self {
        Self {
            mass,
            length,
            wheel_base,
            position: 0.0,
            velocity: 0.0,
            acceleration: 0.0,
        }
    }

    /// Compute the tractive effort (force) produced by the traction motor
    /// at the given normalised throttle setting \[0, 1\].
    ///
    /// A simplified linear model is used: `F = throttle * mass * 0.15 * g`.
    pub fn tractive_effort(&self, throttle: f64) -> f64 {
        let throttle = throttle.clamp(0.0, 1.0);
        throttle * self.mass * 0.15 * 9.81
    }

    /// Compute the braking force produced at the given brake cylinder
    /// pressure (normalised 0–1).
    ///
    /// Uses a simple friction-brake model: `F = p * mass * 0.12 * g`.
    pub fn braking_force(&self, brake_pressure: f64) -> f64 {
        let p = brake_pressure.clamp(0.0, 1.0);
        p * self.mass * 0.12 * 9.81
    }

    /// Compute the total motion resistance of the car (N) at the given
    /// speed (m/s) on a track with the given grade (rise/run, positive
    /// uphill).
    ///
    /// Combines aerodynamic and rolling resistance with grade resistance.
    pub fn resistance(&self, speed: f64, grade: f64) -> f64 {
        let v = speed.abs();
        let a = 6.4;
        let b = 0.09;
        let c = 0.00285;
        let davis = davis_resistance(self.mass, v, a, b, c);
        let grade_force = self.mass * 9.81 * grade;
        davis + grade_force
    }

    /// Integrate the car state forward by `dt` seconds under the given net
    /// force `f_net` (N).
    pub fn integrate(&mut self, dt: f64, f_net: f64) {
        self.acceleration = f_net / self.mass;
        self.velocity += self.acceleration * dt;
        self.position += self.velocity * dt;
    }
}

// ── TrainConsist ──────────────────────────────────────────────────────────────

/// A coupled set of `TrainCar`s forming a complete train consist.
#[derive(Debug, Clone)]
pub struct TrainConsist {
    /// Ordered list of cars from locomotive to rear.
    pub cars: Vec<TrainCar>,
    /// Linear stiffness of the inter-car couplers (N/m).
    pub coupler_stiffness: f64,
    /// Linear damping coefficient of the couplers (N·s/m).
    pub coupler_damping: f64,
}

impl TrainConsist {
    /// Create a new `TrainConsist` from a vector of cars and coupler
    /// parameters.
    pub fn new(cars: Vec<TrainCar>, coupler_stiffness: f64, coupler_damping: f64) -> Self {
        Self {
            cars,
            coupler_stiffness,
            coupler_damping,
        }
    }

    /// Total mass of all cars in the consist (kg).
    pub fn total_mass(&self) -> f64 {
        self.cars.iter().map(|c| c.mass).sum()
    }

    /// Compute the longitudinal force (N) in each coupler between adjacent
    /// cars.  Returns a vector of length `n − 1` where element `i` is the
    /// force between car `i` and car `i+1`.  A positive value means tension;
    /// a negative value means compression.
    pub fn coupler_forces(&self) -> Vec<f64> {
        let n = self.cars.len();
        if n < 2 {
            return Vec::new();
        }
        let mut forces = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let a = &self.cars[i];
            let b = &self.cars[i + 1];
            // Nominal inter-car gap at rest = half-lengths summed.
            let rest_gap = (a.length + b.length) * 0.5;
            let current_gap = b.position - a.position;
            let extension = current_gap - rest_gap;
            let rel_vel = b.velocity - a.velocity;
            let force = self.coupler_stiffness * extension + self.coupler_damping * rel_vel;
            forces.push(force);
        }
        forces
    }

    /// Step the entire consist forward by `dt` seconds.
    ///
    /// The lead car receives `throttle` and every car receives uniform
    /// `brake` pressure.  Coupler forces are distributed accordingly.
    pub fn step(&mut self, dt: f64, throttle: f64, brake: f64) {
        let n = self.cars.len();
        if n == 0 {
            return;
        }
        let coupler = self.coupler_forces();

        // Gather net forces per car.
        let mut net_forces: Vec<f64> = Vec::with_capacity(n);
        for i in 0..n {
            let car = &self.cars[i];
            let speed = car.velocity;
            let grade = 0.0_f64; // flat track for consist-level step
            let res = car.resistance(speed, grade);
            let traction = if i == 0 {
                car.tractive_effort(throttle)
            } else {
                0.0
            };
            let braking = car.braking_force(brake);
            // Sign convention: motion is positive.
            let sign = if speed >= 0.0 { 1.0 } else { -1.0 };
            let mut f_net = traction - braking - sign * res;
            // Apply coupler forces.
            if i > 0 {
                f_net += coupler[i - 1];
            }
            if i < n - 1 {
                f_net -= coupler[i];
            }
            net_forces.push(f_net);
        }

        for (i, car) in self.cars.iter_mut().enumerate() {
            car.integrate(dt, net_forces[i]);
        }
    }

    /// The maximum compressive (bunching) load magnitude across all couplers
    /// (N).
    pub fn bunching_load(&self) -> f64 {
        self.coupler_forces()
            .iter()
            .filter(|&&f| f < 0.0)
            .fold(0.0_f64, |acc, &f| acc.max(-f))
    }

    /// The total run-through brake force across all cars (N) at full brake
    /// pressure.
    pub fn run_through_brake_force(&self) -> f64 {
        self.cars.iter().map(|c| c.braking_force(1.0)).sum()
    }
}

// ── Track ─────────────────────────────────────────────────────────────────────

/// A track profile described by piecewise-constant grade and curvature
/// segments, and a list of station stop positions.
#[derive(Debug, Clone)]
pub struct Track {
    /// Grade profile as `(position, grade)` waypoints.  Grade is rise/run
    /// (positive = uphill).
    pub grade: Vec<(f64, f64)>,
    /// Curvature profile as `(position, curvature)` waypoints (1/radius, m⁻¹).
    pub curvature: Vec<(f64, f64)>,
    /// Positions of station stops along the track (m).
    pub station_positions: Vec<f64>,
}

impl Track {
    /// Construct a new `Track` from the given waypoint vectors and station
    /// positions.
    pub fn new(
        grade: Vec<(f64, f64)>,
        curvature: Vec<(f64, f64)>,
        station_positions: Vec<f64>,
    ) -> Self {
        Self {
            grade,
            curvature,
            station_positions,
        }
    }

    /// Return the grade (rise/run) at the given track position by looking
    /// up the nearest preceding waypoint.  Returns 0.0 for an empty profile.
    pub fn grade_at(&self, pos: f64) -> f64 {
        Self::lookup(&self.grade, pos)
    }

    /// Return the curvature (m⁻¹) at the given track position.
    /// Returns 0.0 for an empty profile.
    pub fn curvature_at(&self, pos: f64) -> f64 {
        Self::lookup(&self.curvature, pos)
    }

    /// Distance from `pos` to the next station ahead (m).  Returns `f64::MAX`
    /// if there are no stations or all are behind `pos`.
    pub fn distance_to_next_station(&self, pos: f64) -> f64 {
        self.station_positions
            .iter()
            .filter(|&&s| s > pos)
            .map(|&s| s - pos)
            .fold(f64::MAX, f64::min)
    }

    fn lookup(table: &[(f64, f64)], pos: f64) -> f64 {
        if table.is_empty() {
            return 0.0;
        }
        let mut last = table[0].1;
        for &(wp, val) in table {
            if pos < wp {
                break;
            }
            last = val;
        }
        last
    }
}

// ── Wheelset ──────────────────────────────────────────────────────────────────

/// A railway wheelset (pair of wheels on a common axle).
#[derive(Debug, Clone)]
pub struct Wheelset {
    /// Static axle load (N).
    pub axle_load: f64,
    /// Wheel rolling radius (m).
    pub wheel_radius: f64,
    /// Kalker creep coefficient (N/m, dimensionless scaling).
    pub creep_coeff: f64,
}

impl Wheelset {
    /// Create a new `Wheelset`.
    ///
    /// * `axle_load`   – Vertical load on the axle (N).
    /// * `wheel_radius` – Nominal rolling radius (m).
    /// * `creep_coeff`  – Kalker dimensionless creep coefficient.
    pub fn new(axle_load: f64, wheel_radius: f64, creep_coeff: f64) -> Self {
        Self {
            axle_load,
            wheel_radius,
            creep_coeff,
        }
    }

    /// Compute the lateral creep force (N) for the given lateral velocity
    /// (m/s) relative to the rail.
    ///
    /// Uses a linear Kalker model: `F = C * Q * (v_lat / r)`.
    pub fn lateral_creep_force(&self, lateral_vel: f64) -> f64 {
        let creep = lateral_vel / self.wheel_radius;
        self.creep_coeff * self.axle_load * creep
    }

    /// Nadal derailment coefficient `Y/Q` (dimensionless).
    ///
    /// Computed from the flange angle (fixed at 60°) and friction coefficient
    /// (0.36) using the Nadal formula.
    pub fn derailment_coefficient(&self) -> f64 {
        let mu = 0.36_f64;
        let delta = 60.0_f64.to_radians();
        (delta.tan() - mu) / (1.0 + mu * delta.tan())
    }

    /// Wear index `T_gamma` (N/mm², proxy for wheel–rail wear rate).
    ///
    /// Simplified: `W = (F_creep^2) / (axle_load * wheel_radius)`.
    pub fn wear_index(&self) -> f64 {
        let f_ref = self.lateral_creep_force(0.001); // 1 mm/s lateral slip
        (f_ref * f_ref) / (self.axle_load * self.wheel_radius + 1e-12)
    }
}

// ── SignalingSystem ───────────────────────────────────────────────────────────

/// A fixed-block signaling system governing safe train separation.
#[derive(Debug, Clone)]
pub struct SignalingSystem {
    /// Length of each signal block (m).
    pub block_length: f64,
    /// Minimum time headway between successive trains (s).
    pub headway: f64,
}

impl SignalingSystem {
    /// Create a new `SignalingSystem`.
    ///
    /// * `block_length` – Fixed block length (m).
    /// * `headway`      – Minimum time headway (s).
    pub fn new(block_length: f64, headway: f64) -> Self {
        Self {
            block_length,
            headway,
        }
    }

    /// Minimum safe braking distance (m) for a train travelling at `speed`
    /// (m/s) with service deceleration `decel` (m/s², positive magnitude).
    ///
    /// `d = v² / (2 * a)` plus one block length as overlap.
    pub fn safe_braking_distance(&self, speed: f64, decel: f64) -> f64 {
        let decel = decel.max(0.01); // prevent division by zero
        let kinematic = (speed * speed) / (2.0 * decel);
        kinematic + self.block_length
    }

    /// Maximum theoretical line capacity (trains per hour) at the given
    /// average speed (m/s).
    pub fn line_capacity(&self, speed: f64) -> f64 {
        let speed = speed.max(0.01);
        let headway_metres = speed * self.headway;
        let effective_spacing = headway_metres.max(self.block_length);
        // Capacity = 3600 / headway_s, or inversely from spacing.
        let _ = effective_spacing; // used for reference; headway-based formula:
        3600.0 / self.headway
    }
}

// ── Davis resistance ─────────────────────────────────────────────────────────

/// Compute the Davis rolling resistance formula (N).
///
/// The classic form is `R = a * W + b * W * v + c * v²` where:
///
/// * `mass`  – Total car mass (kg); weight `W = mass * g`.
/// * `speed` – Speed (m/s).
/// * `a`, `b`, `c` – Empirical Davis coefficients.
pub fn davis_resistance(mass: f64, speed: f64, a: f64, b: f64, c: f64) -> f64 {
    let w = mass * 9.81;
    let v = speed.abs();
    a * w / 1000.0 + b * w * v / 1000.0 + c * v * v
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── TrainCar tests ────────────────────────────────────────────────────

    #[test]
    fn train_car_tractive_effort_zero_throttle() {
        let car = TrainCar::new(50_000.0, 25.0, 14.0);
        assert_eq!(car.tractive_effort(0.0), 0.0);
    }

    #[test]
    fn train_car_tractive_effort_full_throttle() {
        let car = TrainCar::new(50_000.0, 25.0, 14.0);
        let f = car.tractive_effort(1.0);
        // Expected: 1.0 * 50000 * 0.15 * 9.81 = 73575 N
        let expected = 50_000.0 * 0.15 * 9.81;
        assert!((f - expected).abs() < 1.0, "f={f}");
    }

    #[test]
    fn train_car_tractive_effort_clamps_above_one() {
        let car = TrainCar::new(50_000.0, 25.0, 14.0);
        assert_eq!(car.tractive_effort(2.0), car.tractive_effort(1.0));
    }

    #[test]
    fn train_car_braking_force_zero_pressure() {
        let car = TrainCar::new(40_000.0, 20.0, 12.0);
        assert_eq!(car.braking_force(0.0), 0.0);
    }

    #[test]
    fn train_car_braking_force_full_pressure() {
        let car = TrainCar::new(40_000.0, 20.0, 12.0);
        let f = car.braking_force(1.0);
        let expected = 40_000.0 * 0.12 * 9.81;
        assert!((f - expected).abs() < 1.0, "f={f}");
    }

    #[test]
    fn train_car_resistance_uphill_greater_than_flat() {
        let car = TrainCar::new(50_000.0, 25.0, 14.0);
        let r_flat = car.resistance(20.0, 0.0);
        let r_hill = car.resistance(20.0, 0.02); // 2% grade
        assert!(r_hill > r_flat, "uphill resistance should be greater");
    }

    #[test]
    fn train_car_integrate_accelerates_from_rest() {
        let mut car = TrainCar::new(10_000.0, 20.0, 12.0);
        car.integrate(1.0, 5_000.0);
        // a = 5000/10000 = 0.5 m/s²; v = 0.5; x = 0.5
        assert!((car.velocity - 0.5).abs() < 1e-10, "v={}", car.velocity);
        assert!((car.position - 0.5).abs() < 1e-10, "x={}", car.position);
    }

    #[test]
    fn train_car_integrate_decelerates_with_negative_force() {
        let mut car = TrainCar::new(10_000.0, 20.0, 12.0);
        car.velocity = 10.0;
        car.integrate(1.0, -10_000.0);
        // a = -1 m/s²; new v = 9 m/s
        assert!((car.velocity - 9.0).abs() < 1e-10);
    }

    // ── TrainConsist tests ────────────────────────────────────────────────

    #[test]
    fn consist_total_mass_sums_cars() {
        let cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        assert!((consist.total_mass() - 130_000.0).abs() < 1e-6);
    }

    #[test]
    fn consist_coupler_forces_empty_with_one_car() {
        let cars = vec![TrainCar::new(50_000.0, 25.0, 14.0)];
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        assert!(consist.coupler_forces().is_empty());
    }

    #[test]
    fn consist_coupler_forces_tension_when_stretched() {
        let mut cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        // Position car 1 further away than the nominal gap.
        let nominal_gap = (cars[0].length + cars[1].length) * 0.5;
        cars[1].position = nominal_gap + 1.0; // 1 m extra extension
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        let forces = consist.coupler_forces();
        assert_eq!(forces.len(), 1);
        assert!(forces[0] > 0.0, "stretched coupler should be in tension");
    }

    #[test]
    fn consist_coupler_forces_compression_when_bunched() {
        let mut cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        // Position car 1 closer than the nominal gap.
        let nominal_gap = (cars[0].length + cars[1].length) * 0.5;
        cars[1].position = nominal_gap - 1.0; // 1 m compression
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        let forces = consist.coupler_forces();
        assert!(forces[0] < 0.0, "compressed coupler should be negative");
    }

    #[test]
    fn consist_run_through_brake_force_positive() {
        let cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        assert!(consist.run_through_brake_force() > 0.0);
    }

    #[test]
    fn consist_step_advances_lead_car_with_throttle() {
        let cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        let nominal_gap = (25.0 + 20.0) * 0.5;
        let mut consist = TrainConsist::new(cars, 1e6, 5e4);
        consist.cars[1].position = nominal_gap; // place at rest gap
        consist.step(1.0, 0.8, 0.0);
        // Lead car should have moved forward.
        assert!(consist.cars[0].position > 0.0 || consist.cars[0].velocity > 0.0);
    }

    #[test]
    fn consist_bunching_load_zero_when_cars_at_rest_gap() {
        let mut cars = vec![
            TrainCar::new(50_000.0, 25.0, 14.0),
            TrainCar::new(40_000.0, 20.0, 12.0),
        ];
        let nominal_gap = (cars[0].length + cars[1].length) * 0.5;
        cars[1].position = nominal_gap;
        let consist = TrainConsist::new(cars, 1e6, 5e4);
        assert_eq!(consist.bunching_load(), 0.0);
    }

    // ── Track tests ───────────────────────────────────────────────────────

    #[test]
    fn track_grade_at_returns_zero_on_empty_profile() {
        let track = Track::new(vec![], vec![], vec![]);
        assert_eq!(track.grade_at(100.0), 0.0);
    }

    #[test]
    fn track_grade_at_picks_correct_segment() {
        let track = Track::new(
            vec![(0.0, 0.0), (500.0, 0.02), (1000.0, -0.01)],
            vec![],
            vec![],
        );
        assert!((track.grade_at(600.0) - 0.02).abs() < 1e-10);
        assert!((track.grade_at(1200.0) - (-0.01)).abs() < 1e-10);
    }

    #[test]
    fn track_curvature_at_returns_correct_value() {
        let track = Track::new(vec![], vec![(0.0, 0.0), (200.0, 0.005)], vec![]);
        assert!((track.curvature_at(300.0) - 0.005).abs() < 1e-10);
    }

    #[test]
    fn track_distance_to_next_station_forward() {
        let track = Track::new(vec![], vec![], vec![500.0, 1000.0, 2000.0]);
        let d = track.distance_to_next_station(300.0);
        assert!((d - 200.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn track_distance_to_next_station_past_all_stations() {
        let track = Track::new(vec![], vec![], vec![500.0, 1000.0]);
        let d = track.distance_to_next_station(1500.0);
        assert_eq!(d, f64::MAX);
    }

    #[test]
    fn track_distance_to_next_station_no_stations() {
        let track = Track::new(vec![], vec![], vec![]);
        assert_eq!(track.distance_to_next_station(0.0), f64::MAX);
    }

    // ── Wheelset tests ────────────────────────────────────────────────────

    #[test]
    fn wheelset_lateral_creep_force_zero_velocity() {
        let w = Wheelset::new(100_000.0, 0.46, 10.0);
        assert_eq!(w.lateral_creep_force(0.0), 0.0);
    }

    #[test]
    fn wheelset_lateral_creep_force_positive_velocity() {
        let w = Wheelset::new(100_000.0, 0.46, 10.0);
        let f = w.lateral_creep_force(0.01);
        assert!(f > 0.0, "positive lateral vel should give positive force");
    }

    #[test]
    fn wheelset_lateral_creep_force_scales_linearly() {
        let w = Wheelset::new(100_000.0, 0.46, 10.0);
        let f1 = w.lateral_creep_force(0.01);
        let f2 = w.lateral_creep_force(0.02);
        assert!((f2 - 2.0 * f1).abs() < 1e-9, "f1={f1}, f2={f2}");
    }

    #[test]
    fn wheelset_derailment_coefficient_positive() {
        let w = Wheelset::new(100_000.0, 0.46, 10.0);
        let yq = w.derailment_coefficient();
        assert!(yq > 0.0, "derailment coefficient must be positive: {yq}");
    }

    #[test]
    fn wheelset_wear_index_positive() {
        let w = Wheelset::new(100_000.0, 0.46, 10.0);
        assert!(w.wear_index() > 0.0);
    }

    // ── SignalingSystem tests ─────────────────────────────────────────────

    #[test]
    fn signaling_braking_distance_stationary() {
        let sig = SignalingSystem::new(200.0, 90.0);
        let d = sig.safe_braking_distance(0.0, 1.0);
        // kinematic = 0; result = block_length = 200
        assert!((d - 200.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn signaling_braking_distance_moving() {
        let sig = SignalingSystem::new(200.0, 90.0);
        let d = sig.safe_braking_distance(30.0, 1.0);
        // kinematic = 900/2 = 450; total = 650
        assert!((d - 650.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn signaling_line_capacity_positive() {
        let sig = SignalingSystem::new(200.0, 90.0);
        let cap = sig.line_capacity(27.8);
        assert!(cap > 0.0, "capacity={cap}");
        // 3600 / 90 = 40 trains/hr
        assert!(
            (cap - 40.0).abs() < 1e-10,
            "expected 40 trains/hr, got {cap}"
        );
    }

    // ── Davis resistance tests ────────────────────────────────────────────

    #[test]
    fn davis_resistance_zero_speed() {
        // At zero speed the c*v² and b*W*v terms vanish.
        let r = davis_resistance(50_000.0, 0.0, 6.4, 0.09, 0.00285);
        let w = 50_000.0 * 9.81;
        let expected = 6.4 * w / 1000.0;
        assert!((r - expected).abs() < 1e-6, "r={r}");
    }

    #[test]
    fn davis_resistance_increases_with_speed() {
        let r_low = davis_resistance(50_000.0, 10.0, 6.4, 0.09, 0.00285);
        let r_high = davis_resistance(50_000.0, 30.0, 6.4, 0.09, 0.00285);
        assert!(r_high > r_low);
    }

    #[test]
    fn davis_resistance_scales_with_mass() {
        let r1 = davis_resistance(10_000.0, 20.0, 6.4, 0.09, 0.00285);
        let r2 = davis_resistance(20_000.0, 20.0, 6.4, 0.09, 0.00285);
        // Mass-dependent terms double; c*v² term does not depend on mass.
        assert!(r2 > r1);
    }

    #[test]
    fn davis_resistance_non_negative() {
        let r = davis_resistance(50_000.0, 15.0, 6.4, 0.09, 0.00285);
        assert!(r >= 0.0, "resistance must be non-negative");
    }
}
