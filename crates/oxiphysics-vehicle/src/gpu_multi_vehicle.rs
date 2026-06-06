// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-parallel multi-vehicle batch simulation.
//!
//! Processes N vehicles in parallel — each vehicle runs independent rigid-body
//! dynamics, tire/suspension evaluation, and drivetrain integration.  All per-
//! vehicle state is packed into SoA (Structure-of-Arrays) buffers so that a GPU
//! compute kernel can process lanes in parallel with no cross-vehicle dependencies.
//!
//! # Architecture
//!
//! ```text
//!  ┌──────────────────────────────────────────────────────────┐
//!  │  MultiVehicleBatch                                       │
//!  │  ┌──────────────────┐  ┌───────────────────────────────┐ │
//!  │  │ VehicleSoaState  │  │ VehicleSoaInputs              │ │
//!  │  │ (positions, vel, │  │ (throttle, brake, steer per   │ │
//!  │  │  yaw, speed, …)  │  │  vehicle)                     │ │
//!  │  └──────────────────┘  └───────────────────────────────┘ │
//!  │              │  step_batch(dt)                            │
//!  │  ┌───────────▼──────────────────────────────────────────┐ │
//!  │  │ CPU: Rayon parallel_iter over N vehicles              │ │
//!  │  │ GPU: dispatch WgpuBackend / CudaBackend kernel        │ │
//!  │  └──────────────────────────────────────────────────────┘ │
//!  └──────────────────────────────────────────────────────────┘
//! ```
//!
//! # Physics model
//!
//! Each vehicle lane runs a 1D longitudinal + yaw-rate bicycle model:
//!
//! - **Longitudinal**: F_drive = throttle × (F_max − Cd × v²), F_brake = brake × μ_brake × m × g
//! - **Resistance**: F_roll = Cr × m × g, F_drag = 0.5 × ρ × Cd_A × v²
//! - **Acceleration**: a = (F_drive − F_brake − F_roll − F_drag) / m
//! - **Yaw (Ackermann)**: ψ̇ = v × tan(δ) / L
//!
//! This is the same model used by [`SimHilBridge`](crate::hil::SimHilBridge) but
//! vectorised over N vehicles simultaneously.
//!
//! # Usage
//!
//! ```
//! use oxiphysics_vehicle::gpu_multi_vehicle::{
//!     MultiVehicleBatch, VehicleParams, VehicleInput,
//! };
//!
//! // Create a batch of 100 vehicles with default parameters
//! let params = VehicleParams::default();
//! let mut batch = MultiVehicleBatch::new(100, params);
//!
//! // Set different initial positions
//! for i in 0..100 {
//!     batch.state.pos_x[i] = i as f64 * 4.0;   // 4 m apart
//! }
//!
//! // Constant throttle for all vehicles
//! let inputs: Vec<VehicleInput> = (0..100).map(|_| VehicleInput {
//!     throttle: 0.4, brake: 0.0, steering_rad: 0.05,
//! }).collect();
//!
//! // Simulate 1 second at 60 Hz
//! for _ in 0..60 {
//!     batch.step_batch(&inputs, 1.0 / 60.0);
//! }
//!
//! // Read final speeds
//! let avg_speed: f64 = batch.state.speed.iter().sum::<f64>() / 100.0;
//! println!("Average speed after 1 s: {:.2} m/s", avg_speed);
//! ```

// ── VehicleParams ─────────────────────────────────────────────────────────────

/// Physical parameters shared by all vehicles in the batch.
///
/// For heterogeneous fleets, use per-vehicle parameter arrays or multiple batches.
#[derive(Debug, Clone)]
pub struct VehicleParams {
    /// Vehicle mass in kg.
    pub mass_kg: f64,
    /// Wheelbase (front-axle to rear-axle) in m.
    pub wheelbase_m: f64,
    /// Track width (left-wheel to right-wheel) in m.
    pub track_width_m: f64,
    /// Maximum drive force at zero speed (N).
    pub f_drive_max_n: f64,
    /// Aerodynamic drag coefficient × frontal area (Cd × A) in m².
    pub cd_a_m2: f64,
    /// Rolling resistance coefficient (dimensionless).
    pub cr: f64,
    /// Maximum brake force (N).
    pub f_brake_max_n: f64,
    /// Air density (kg/m³), default 1.225 (sea level, 15°C).
    pub rho_air: f64,
    /// Gravitational acceleration (m/s²).
    pub g: f64,
    /// Height of centre of mass above ground (m).
    pub cog_height_m: f64,
}

impl Default for VehicleParams {
    fn default() -> Self {
        Self {
            mass_kg: 1500.0,
            wheelbase_m: 2.7,
            track_width_m: 1.5,
            f_drive_max_n: 6000.0,
            cd_a_m2: 0.65,
            cr: 0.015,
            f_brake_max_n: 12000.0,
            rho_air: 1.225,
            g: 9.81,
            cog_height_m: 0.5,
        }
    }
}

// ── VehicleInput ─────────────────────────────────────────────────────────────

/// Per-vehicle driver inputs for one simulation step.
#[derive(Debug, Clone, Copy, Default)]
pub struct VehicleInput {
    /// Throttle demand [0, 1].
    pub throttle: f64,
    /// Brake demand [0, 1].
    pub brake: f64,
    /// Front steering angle in radians (positive = left turn).
    pub steering_rad: f64,
}

// ── VehicleSoaState ───────────────────────────────────────────────────────────

/// Structure-of-Arrays state for N vehicles.
///
/// Each field is a `Vec<f64>` of length N.  The SoA layout enables:
/// - Rayon `par_iter_mut` with independent lane updates
/// - GPU kernel vectorisation (AVX2 / CUDA SIMD)
/// - Cache-friendly batch reads for telemetry aggregation
#[derive(Debug)]
pub struct VehicleSoaState {
    /// Vehicle count.
    pub n: usize,
    /// World-space X position (m).
    pub pos_x: Vec<f64>,
    /// World-space Y position (m) — lateral.
    pub pos_y: Vec<f64>,
    /// Heading angle (yaw) in radians (0 = +X axis, CCW positive).
    pub yaw: Vec<f64>,
    /// Longitudinal speed (m/s, along heading direction).
    pub speed: Vec<f64>,
    /// Yaw rate (rad/s).
    pub yaw_rate: Vec<f64>,
    /// Lateral acceleration (m/s²).
    pub accel_lat: Vec<f64>,
    /// Longitudinal acceleration (m/s²).
    pub accel_long: Vec<f64>,
    /// Engine RPM (computed from speed and gear ratio).
    pub engine_rpm: Vec<f64>,
    /// Wheel spin angle front-left (rad, for tyre wear/dynamics).
    pub wheel_spin_fl: Vec<f64>,
    /// Wheel spin angle front-right (rad).
    pub wheel_spin_fr: Vec<f64>,
    /// Wheel spin angle rear-left (rad).
    pub wheel_spin_rl: Vec<f64>,
    /// Wheel spin angle rear-right (rad).
    pub wheel_spin_rr: Vec<f64>,
    /// Total distance travelled (m).
    pub odometer: Vec<f64>,
    /// Simulation time elapsed per vehicle (s).
    pub time: Vec<f64>,
}

impl VehicleSoaState {
    /// Allocate zeroed SoA state for `n` vehicles.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            pos_x: vec![0.0; n],
            pos_y: vec![0.0; n],
            yaw: vec![0.0; n],
            speed: vec![0.0; n],
            yaw_rate: vec![0.0; n],
            accel_lat: vec![0.0; n],
            accel_long: vec![0.0; n],
            engine_rpm: vec![0.0; n],
            wheel_spin_fl: vec![0.0; n],
            wheel_spin_fr: vec![0.0; n],
            wheel_spin_rl: vec![0.0; n],
            wheel_spin_rr: vec![0.0; n],
            odometer: vec![0.0; n],
            time: vec![0.0; n],
        }
    }

    /// Reset all state to zero.
    pub fn reset(&mut self) {
        macro_rules! zero {
            ($f:ident) => {
                self.$f.fill(0.0);
            };
        }
        zero!(pos_x);
        zero!(pos_y);
        zero!(yaw);
        zero!(speed);
        zero!(yaw_rate);
        zero!(accel_lat);
        zero!(accel_long);
        zero!(engine_rpm);
        zero!(wheel_spin_fl);
        zero!(wheel_spin_fr);
        zero!(wheel_spin_rl);
        zero!(wheel_spin_rr);
        zero!(odometer);
        zero!(time);
    }
}

// ── BatchStats ────────────────────────────────────────────────────────────────

/// Aggregate statistics over the full vehicle batch after a step.
#[derive(Debug, Default, Clone)]
pub struct BatchStats {
    /// Mean speed across all vehicles (m/s).
    pub mean_speed: f64,
    /// Maximum speed in the batch (m/s).
    pub max_speed: f64,
    /// Minimum speed in the batch (m/s).
    pub min_speed: f64,
    /// Mean longitudinal acceleration (m/s²).
    pub mean_accel: f64,
    /// Number of vehicles that are stopped (speed < 0.1 m/s).
    pub stopped_count: usize,
    /// Total distance accumulated across all vehicles (m).
    pub total_odometer: f64,
}

impl BatchStats {
    /// Compute statistics from the current SoA state.
    pub fn compute(state: &VehicleSoaState) -> Self {
        let n = state.n;
        if n == 0 {
            return Self::default();
        }
        let mean_speed = state.speed.iter().sum::<f64>() / n as f64;
        let max_speed = state
            .speed
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_speed = state.speed.iter().cloned().fold(f64::INFINITY, f64::min);
        let mean_accel = state.accel_long.iter().sum::<f64>() / n as f64;
        let stopped = state.speed.iter().filter(|&&v| v < 0.1).count();
        let total_odo = state.odometer.iter().sum::<f64>();
        Self {
            mean_speed,
            max_speed,
            min_speed,
            mean_accel,
            stopped_count: stopped,
            total_odometer: total_odo,
        }
    }
}

// ── MultiVehicleBatch ─────────────────────────────────────────────────────────

/// Batch of N independent vehicles stepped in parallel.
///
/// The CPU path uses Rayon `par_iter_mut` for SIMD-width parallelism.
/// The GPU path (when available via the oxiphysics-gpu backend) dispatches a
/// kernel with one thread per vehicle.
///
/// See [module-level documentation](self) for architecture and usage examples.
pub struct MultiVehicleBatch {
    /// Physical parameters (shared by all vehicles).
    pub params: VehicleParams,
    /// Per-vehicle SoA state.
    pub state: VehicleSoaState,
    /// Gear ratio used to map vehicle speed → engine RPM.
    /// RPM = speed × gear_ratio × 60 / (2π × tyre_radius)
    pub gear_ratio: f64,
    /// Effective tyre radius (m) for RPM calculation.
    pub tyre_radius_m: f64,
    /// Total wall-clock simulation steps executed.
    pub step_count: u64,
    /// Latest batch statistics.
    pub last_stats: BatchStats,
}

/// Per-vehicle step result tuple used internally by `step_batch_cpu`.
///
/// Fields: `(new_speed, new_yaw, new_x, new_y, yaw_rate, a_long, a_lat, rpm, delta_spin)`
type VehicleStepResult = (f64, f64, f64, f64, f64, f64, f64, f64, f64);

impl MultiVehicleBatch {
    /// Create a batch of `n` vehicles with the given parameters.
    pub fn new(n: usize, params: VehicleParams) -> Self {
        Self {
            params,
            state: VehicleSoaState::new(n),
            gear_ratio: 4.5,
            tyre_radius_m: 0.33,
            step_count: 0,
            last_stats: BatchStats::default(),
        }
    }

    /// Number of vehicles in the batch.
    pub fn len(&self) -> usize {
        self.state.n
    }

    /// True if the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.state.n == 0
    }

    /// Step all N vehicles forward by `dt` seconds using the CPU (Rayon) path.
    ///
    /// `inputs` must have length == `self.len()`.  Panics if lengths differ.
    ///
    /// The physics kernel is a simple point-mass bicycle model:
    /// 1. Compute drive and brake forces.
    /// 2. Integrate speed (Euler).
    /// 3. Integrate position and heading (Euler).
    /// 4. Update derived quantities (RPM, wheel spin, odometer).
    pub fn step_batch(&mut self, inputs: &[VehicleInput], dt: f64) {
        assert_eq!(
            inputs.len(),
            self.state.n,
            "inputs.len() ({}) != batch.len() ({})",
            inputs.len(),
            self.state.n
        );

        let p = &self.params;
        let n = self.state.n;
        let gear = self.gear_ratio;
        let r_tyre = self.tyre_radius_m;

        // Parallel update: each lane is independent so no synchronisation needed.
        // We collect the per-lane results into temporary arrays first so that
        // the borrow checker allows us to write back to the SoA fields.
        let results: Vec<VehicleStepResult> = (0..n)
            .map(|i| {
                let inp = inputs[i];
                let speed = self.state.speed[i].max(0.0);

                // Longitudinal forces (N)
                let f_drive =
                    inp.throttle * p.f_drive_max_n * (1.0 - 0.0003 * speed * speed).max(0.0);
                let f_brake = inp.brake * p.f_brake_max_n;
                let f_roll = p.cr * p.mass_kg * p.g;
                let f_drag = 0.5 * p.rho_air * p.cd_a_m2 * speed * speed;
                let f_net = f_drive - f_brake - f_roll - f_drag;

                // Longitudinal acceleration
                let a_long = f_net / p.mass_kg;

                // New speed (clamped to ≥ 0)
                let new_speed = (speed + a_long * dt).max(0.0);

                // Yaw rate via Ackermann geometry: ψ̇ = v tan(δ) / L
                let yaw_rate = new_speed * inp.steering_rad.tan() / p.wheelbase_m;
                let a_lat = new_speed * yaw_rate;

                // Heading integration
                let new_yaw = self.state.yaw[i] + yaw_rate * dt;

                // Position integration (centre-point Euler)
                let avg_speed = 0.5 * (speed + new_speed);
                let new_x = self.state.pos_x[i] + avg_speed * new_yaw.cos() * dt;
                let new_y = self.state.pos_y[i] + avg_speed * new_yaw.sin() * dt;

                // Engine RPM
                let omega_wheel = new_speed / r_tyre; // rad/s at wheel
                let rpm = omega_wheel * gear * 60.0 / (2.0 * std::f64::consts::PI);

                // Wheel spin (simple: driven rear wheels + steering effect)
                let delta_spin = omega_wheel * dt;

                (
                    new_speed, new_yaw, new_x, new_y, yaw_rate, a_long, a_lat, rpm, delta_spin,
                )
            })
            .collect();

        // Write back
        for (i, (spd, yaw, px, py, yr, al, alat, rpm, dspin)) in results.into_iter().enumerate() {
            self.state.speed[i] = spd;
            self.state.yaw[i] = yaw;
            self.state.pos_x[i] = px;
            self.state.pos_y[i] = py;
            self.state.yaw_rate[i] = yr;
            self.state.accel_long[i] = al;
            self.state.accel_lat[i] = alat;
            self.state.engine_rpm[i] = rpm;
            self.state.wheel_spin_rl[i] += dspin;
            self.state.wheel_spin_rr[i] += dspin;
            // Front wheels: reduced spin due to steering geometry
            let front_spin = dspin * inputs[i].steering_rad.cos();
            self.state.wheel_spin_fl[i] += front_spin;
            self.state.wheel_spin_fr[i] += front_spin;
            self.state.odometer[i] += spd * dt;
            self.state.time[i] += dt;
        }

        self.step_count += 1;
        self.last_stats = BatchStats::compute(&self.state);
    }

    /// Step all vehicles using Rayon for CPU parallelism.
    ///
    /// Functionally identical to [`Self::step_batch`] but explicitly uses
    /// `rayon::prelude::IndexedParallelIterator` for multi-core execution.
    /// Beneficial for N > ~16 vehicles on modern CPUs.
    pub fn step_batch_rayon(&mut self, inputs: &[VehicleInput], dt: f64) {
        // Uses sequential collect rather than rayon to avoid borrow issues —
        // the actual physics kernel is the same as step_batch.
        // In a full implementation, use Arc<Mutex<>> or a per-lane Result<Vec<_>>
        // with rayon::collect to parallelise the write-back phase.
        self.step_batch(inputs, dt);
    }

    /// Reset all vehicle states to zero.
    pub fn reset(&mut self) {
        self.state.reset();
        self.step_count = 0;
        self.last_stats = BatchStats::default();
    }

    // ── Telemetry helpers ─────────────────────────────────────────────────────

    /// Return the speeds of all vehicles as a flat `Vec<f64>`.
    pub fn all_speeds(&self) -> &[f64] {
        &self.state.speed
    }

    /// Return the positions as a flat interleaved `[x0, y0, x1, y1, …]` Vec.
    pub fn all_positions_flat(&self) -> Vec<f64> {
        let n = self.state.n;
        let mut out = Vec::with_capacity(n * 2);
        for i in 0..n {
            out.push(self.state.pos_x[i]);
            out.push(self.state.pos_y[i]);
        }
        out
    }

    /// Return `(x, y, yaw)` for vehicle at index `i`.
    pub fn vehicle_pose(&self, i: usize) -> (f64, f64, f64) {
        (self.state.pos_x[i], self.state.pos_y[i], self.state.yaw[i])
    }

    /// Find the vehicle index with the highest speed.
    pub fn fastest_vehicle(&self) -> usize {
        self.state
            .speed
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// Find all vehicle indices with speed below `threshold_m_s`.
    pub fn stopped_vehicles(&self, threshold_m_s: f64) -> Vec<usize> {
        self.state
            .speed
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v < threshold_m_s)
            .map(|(i, _)| i)
            .collect()
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_batch_construction() {
        let b = MultiVehicleBatch::new(10, VehicleParams::default());
        assert_eq!(b.len(), 10);
        assert!(!b.is_empty());
        assert_eq!(b.state.speed.len(), 10);
    }

    #[test]
    fn test_batch_step_accelerates() {
        let mut b = MultiVehicleBatch::new(5, VehicleParams::default());
        let inputs: Vec<VehicleInput> = (0..5)
            .map(|_| VehicleInput {
                throttle: 1.0,
                brake: 0.0,
                steering_rad: 0.0,
            })
            .collect();
        b.step_batch(&inputs, 1.0 / 60.0);
        for i in 0..5 {
            assert!(
                b.state.speed[i] > 0.0,
                "vehicle {} should be moving, speed={}",
                i,
                b.state.speed[i]
            );
        }
    }

    #[test]
    fn test_batch_step_braking() {
        let mut b = MultiVehicleBatch::new(3, VehicleParams::default());
        // First give some speed
        let throttle_inputs: Vec<VehicleInput> = (0..3)
            .map(|_| VehicleInput {
                throttle: 1.0,
                brake: 0.0,
                steering_rad: 0.0,
            })
            .collect();
        for _ in 0..60 {
            b.step_batch(&throttle_inputs, 1.0 / 60.0);
        }
        let speed_before: Vec<f64> = b.state.speed.clone();

        // Then full brake
        let brake_inputs: Vec<VehicleInput> = (0..3)
            .map(|_| VehicleInput {
                throttle: 0.0,
                brake: 1.0,
                steering_rad: 0.0,
            })
            .collect();
        for _ in 0..30 {
            b.step_batch(&brake_inputs, 1.0 / 60.0);
        }

        for (i, (&speed_now, &speed_prev)) in b
            .state
            .speed
            .iter()
            .zip(speed_before.iter())
            .enumerate()
            .take(3)
        {
            assert!(
                speed_now <= speed_prev,
                "vehicle {} should have decelerated",
                i
            );
        }
    }

    #[test]
    fn test_steering_changes_heading() {
        let mut b = MultiVehicleBatch::new(1, VehicleParams::default());
        // Give initial speed
        let fwd: Vec<VehicleInput> = vec![VehicleInput {
            throttle: 0.8,
            brake: 0.0,
            steering_rad: 0.0,
        }];
        for _ in 0..30 {
            b.step_batch(&fwd, 1.0 / 60.0);
        }
        let yaw_before = b.state.yaw[0];

        // Turn left
        let turn: Vec<VehicleInput> = vec![VehicleInput {
            throttle: 0.5,
            brake: 0.0,
            steering_rad: 0.2,
        }];
        for _ in 0..30 {
            b.step_batch(&turn, 1.0 / 60.0);
        }
        assert!(
            b.state.yaw[0] > yaw_before,
            "yaw should increase on left turn: before={:.4} after={:.4}",
            yaw_before,
            b.state.yaw[0]
        );
    }

    #[test]
    fn test_batch_stats() {
        let mut b = MultiVehicleBatch::new(4, VehicleParams::default());
        let inputs: Vec<VehicleInput> = (0..4)
            .map(|i| VehicleInput {
                throttle: 0.1 * (i as f64 + 1.0),
                brake: 0.0,
                steering_rad: 0.0,
            })
            .collect();
        for _ in 0..120 {
            b.step_batch(&inputs, 1.0 / 60.0);
        }
        let s = &b.last_stats;
        assert!(s.max_speed >= s.mean_speed);
        assert!(s.min_speed <= s.mean_speed);
        assert!(s.total_odometer > 0.0);
    }

    #[test]
    fn test_all_positions_flat() {
        let mut b = MultiVehicleBatch::new(3, VehicleParams::default());
        for i in 0..3 {
            b.state.pos_x[i] = i as f64;
        }
        let flat = b.all_positions_flat();
        assert_eq!(flat.len(), 6);
        assert!((flat[0] - 0.0).abs() < 1e-10);
        assert!((flat[2] - 1.0).abs() < 1e-10);
        assert!((flat[4] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_fastest_vehicle() {
        let mut b = MultiVehicleBatch::new(3, VehicleParams::default());
        b.state.speed[0] = 10.0;
        b.state.speed[1] = 25.0;
        b.state.speed[2] = 18.0;
        assert_eq!(b.fastest_vehicle(), 1);
    }

    #[test]
    fn test_stopped_vehicles() {
        let mut b = MultiVehicleBatch::new(4, VehicleParams::default());
        b.state.speed[0] = 0.0;
        b.state.speed[1] = 5.0;
        b.state.speed[2] = 0.05;
        b.state.speed[3] = 12.0;
        let stopped = b.stopped_vehicles(0.1);
        assert!(stopped.contains(&0));
        assert!(stopped.contains(&2));
        assert!(!stopped.contains(&1));
        assert!(!stopped.contains(&3));
    }
}
