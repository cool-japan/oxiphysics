// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real-time hardware-in-the-loop (HiL) interfaces for vehicle simulation.
//!
//! This module provides the [`HilInterface`] trait and supporting types for
//! connecting the OxiPhysics vehicle simulation to real-time hardware — such as
//! dSpace, NI VeriStand, ETAS LABCAR, or custom embedded controllers.
//!
//! ## Architecture
//!
//! ```text
//!  ┌──────────────────┐    HilChannel I/O    ┌──────────────────────┐
//!  │  Vehicle Solver  │◄────────────────────►│  HilBridge           │
//!  │  (OxiPhysics)    │  dt=fixed step       │  (impl HilInterface) │
//!  └──────────────────┘                      └──────────┬───────────┘
//!                                                       │
//!                                              TCP/UDP / shared-mem
//!                                                       │
//!                                            ┌──────────▼───────────┐
//!                                            │  ECU / Actuator HW   │
//!                                            └──────────────────────┘
//! ```
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_vehicle::hil::{HilChannel, HilConfig, HilInterface, SimHilBridge};
//!
//! let cfg = HilConfig::default();
//! let mut bridge = SimHilBridge::new(cfg);
//!
//! // Inject a steering input
//! bridge.set_input(HilChannel::SteeringAngle, 0.05_f64.to_radians());
//!
//! // Step the HiL bridge one frame
//! bridge.tick(1.0 / 1000.0);
//!
//! // Read back the wheel speed output
//! let omega = bridge.get_output(HilChannel::WheelSpeedFl);
//! assert!(omega.unwrap_or(0.0) >= 0.0);
//! ```

use std::collections::HashMap;

// ── HilChannel ───────────────────────────────────────────────────────────────

/// Named signal channels for the HiL interface.
///
/// Each variant represents a physical signal exchanged between the simulation
/// and external hardware.  Values use SI units (radians, m/s, N, Pa, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HilChannel {
    // ── Actuator inputs (hardware → simulation) ──────────────────────────
    /// Brake pedal position in `[0, 1]` (0 = released, 1 = fully depressed).
    BrakePedal,
    /// Throttle pedal position in `[0, 1]`.
    ThrottlePedal,
    /// Steering angle at the rack (rad, positive = left).
    SteeringAngle,
    /// Clutch pedal position in `[0, 1]`.
    ClutchPedal,
    /// Selected gear index (`-1` = reverse, `0` = neutral, `1..` = forward).
    GearSelector,
    /// Engine torque demand override (Nm).  `NaN` = not overridden.
    EngineTorqueOverride,
    /// Individual wheel brake torque (Nm, 0 = no braking).
    WheelBrakeFl,
    /// Individual wheel brake torque (Nm).
    WheelBrakeFr,
    /// Individual wheel brake torque (Nm).
    WheelBrakeRl,
    /// Individual wheel brake torque (Nm).
    WheelBrakeRr,

    // ── Sensor outputs (simulation → hardware) ───────────────────────────
    /// Front-left wheel speed (rad/s).
    WheelSpeedFl,
    /// Front-right wheel speed (rad/s).
    WheelSpeedFr,
    /// Rear-left wheel speed (rad/s).
    WheelSpeedRl,
    /// Rear-right wheel speed (rad/s).
    WheelSpeedRr,
    /// Vehicle longitudinal acceleration (m/s²).
    AccelLongitudinal,
    /// Vehicle lateral acceleration (m/s²).
    AccelLateral,
    /// Vehicle yaw rate (rad/s).
    YawRate,
    /// Vehicle body roll angle (rad).
    RollAngle,
    /// Vehicle speed at centre of mass (m/s).
    VehicleSpeed,
    /// Engine speed (rpm).
    EngineRpm,
    /// Suspension deflection front-left (m).
    SuspDeflFl,
    /// Suspension deflection front-right (m).
    SuspDeflFr,
    /// Suspension deflection rear-left (m).
    SuspDeflRl,
    /// Suspension deflection rear-right (m).
    SuspDeflRr,
    /// Fuel/energy remaining (J or kg, context-dependent).
    EnergyRemaining,
}

// ── HilConfig ────────────────────────────────────────────────────────────────

/// Configuration for the HiL bridge.
#[derive(Debug, Clone)]
pub struct HilConfig {
    /// Real-time step size (seconds).  Typical values: 0.001 (1 ms) for fast ECU, 0.01 (10 ms).
    pub dt: f64,
    /// Maximum allowed overrun before raising a [`HilError::Overrun`] (in units of `dt`).
    pub overrun_tolerance: u32,
    /// Whether to apply low-pass filtering to input channels.
    pub input_filter_enabled: bool,
    /// Low-pass cut-off frequency for input filtering (Hz).
    pub input_filter_hz: f64,
    /// Whether to log timing statistics.
    pub timing_log: bool,
}

impl Default for HilConfig {
    fn default() -> Self {
        Self {
            dt: 1e-3,
            overrun_tolerance: 5,
            input_filter_enabled: true,
            input_filter_hz: 200.0,
            timing_log: false,
        }
    }
}

// ── HilError ─────────────────────────────────────────────────────────────────

/// Errors that can be returned from the HiL interface.
#[derive(Debug, Clone, PartialEq)]
pub enum HilError {
    /// The simulation step took longer than `overrun_tolerance * dt`.
    Overrun { step: u64, overrun_dt: f64 },
    /// Attempted to read a channel that has no value.
    ChannelNotConnected(HilChannel),
    /// The external hardware connection was lost.
    ConnectionLost(String),
    /// A channel value was outside its valid range.
    OutOfRange { channel: HilChannel, value: f64 },
}

impl std::fmt::Display for HilError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HilError::Overrun { step, overrun_dt } => {
                write!(f, "HiL overrun at step {step}: {overrun_dt:.4} s late")
            }
            HilError::ChannelNotConnected(ch) => {
                write!(f, "HiL channel {ch:?} not connected")
            }
            HilError::ConnectionLost(s) => write!(f, "HiL connection lost: {s}"),
            HilError::OutOfRange { channel, value } => {
                write!(f, "HiL channel {channel:?} value {value} out of range")
            }
        }
    }
}

// ── HilInterface trait ────────────────────────────────────────────────────────

/// Trait that any HiL transport implementation must fulfil.
///
/// Implementors include:
/// * [`SimHilBridge`] — software-in-the-loop (pure simulation, no real HW)
/// * A future `TcpHilBridge` — XCP/CAN-over-TCP to dSpace/ETAS hardware
/// * A future `SharedMemHilBridge` — zero-copy via shared memory on a real-time OS
pub trait HilInterface {
    /// Write an input value to the simulation (hardware → simulation).
    fn set_input(&mut self, channel: HilChannel, value: f64);

    /// Read an output value from the simulation (simulation → hardware).
    ///
    /// Returns `None` if the channel has not yet produced a value.
    fn get_output(&self, channel: HilChannel) -> Option<f64>;

    /// Advance the HiL bridge by exactly one `dt` tick.
    ///
    /// Returns `Ok(())` on success or `Err(HilError::Overrun)` if the
    /// deadline was missed.
    fn tick(&mut self, dt: f64) -> Result<(), HilError>;

    /// Return all channels that currently have connected outputs.
    fn available_outputs(&self) -> Vec<HilChannel>;

    /// Return the cumulative number of completed steps.
    fn step_count(&self) -> u64;

    /// Reset the HiL bridge to its initial state.
    fn reset(&mut self);
}

// ── HilTimingStats ───────────────────────────────────────────────────────────

/// Per-step timing statistics for HiL monitoring.
#[derive(Debug, Clone, Default)]
pub struct HilTimingStats {
    /// Total number of steps executed.
    pub total_steps: u64,
    /// Number of steps that exceeded the deadline.
    pub overruns: u64,
    /// Maximum observed step duration (seconds).
    pub max_step_dt: f64,
    /// Minimum observed step duration (seconds).
    pub min_step_dt: f64,
    /// Mean step duration (seconds).
    pub mean_step_dt: f64,
    /// Accumulated sum for mean computation.
    pub(crate) sum_step_dt: f64,
}

impl HilTimingStats {
    /// Record one step with observed duration `elapsed`.
    pub fn record(&mut self, elapsed: f64) {
        self.total_steps += 1;
        self.sum_step_dt += elapsed;
        self.mean_step_dt = self.sum_step_dt / self.total_steps as f64;
        if elapsed > self.max_step_dt {
            self.max_step_dt = elapsed;
        }
        if self.min_step_dt == 0.0 || elapsed < self.min_step_dt {
            self.min_step_dt = elapsed;
        }
    }
}

// ── SimHilBridge (Software-in-the-Loop) ──────────────────────────────────────

/// Software-in-the-loop (SiL) HiL bridge for testing without real hardware.
///
/// All input/output values are held in memory.  `tick()` applies simple
/// physics models to propagate inputs to outputs (wheel speed from vehicle speed,
/// etc.) so that test harnesses can verify signal flow without a real ECU.
#[derive(Debug)]
pub struct SimHilBridge {
    /// Bridge configuration.
    pub config: HilConfig,
    /// Raw (unfiltered) input values as set by the caller.
    inputs: HashMap<HilChannel, f64>,
    /// Low-pass filtered input values used by the physics model.
    /// Tracks the IIR state so filtering doesn't corrupt the raw inputs.
    filtered_inputs: HashMap<HilChannel, f64>,
    /// Computed output values.
    outputs: HashMap<HilChannel, f64>,
    /// Step counter.
    step_count: u64,
    /// Accumulated simulation time.
    sim_time: f64,
    /// Timing statistics.
    pub timing: HilTimingStats,
    /// Simple vehicle state: longitudinal speed (m/s).
    speed: f64,
    /// Simple yaw rate (rad/s).
    yaw_rate: f64,
    /// Engine RPM.
    engine_rpm: f64,
}

impl SimHilBridge {
    /// Create a new SiL bridge with the given configuration.
    pub fn new(config: HilConfig) -> Self {
        let mut outputs = HashMap::new();
        // Initialise all sensor outputs to zero
        for ch in Self::all_output_channels() {
            outputs.insert(ch, 0.0);
        }
        Self {
            config,
            inputs: HashMap::new(),
            filtered_inputs: HashMap::new(),
            outputs,
            step_count: 0,
            sim_time: 0.0,
            timing: HilTimingStats::default(),
            speed: 0.0,
            yaw_rate: 0.0,
            engine_rpm: 800.0, // idle
        }
    }

    /// All channels that this bridge produces as outputs.
    pub fn all_output_channels() -> Vec<HilChannel> {
        vec![
            HilChannel::WheelSpeedFl,
            HilChannel::WheelSpeedFr,
            HilChannel::WheelSpeedRl,
            HilChannel::WheelSpeedRr,
            HilChannel::AccelLongitudinal,
            HilChannel::AccelLateral,
            HilChannel::YawRate,
            HilChannel::RollAngle,
            HilChannel::VehicleSpeed,
            HilChannel::EngineRpm,
            HilChannel::SuspDeflFl,
            HilChannel::SuspDeflFr,
            HilChannel::SuspDeflRl,
            HilChannel::SuspDeflRr,
            HilChannel::EnergyRemaining,
        ]
    }

    /// Current simulation time (seconds).
    pub fn sim_time(&self) -> f64 {
        self.sim_time
    }

    /// Shortcut to set multiple inputs at once.
    pub fn set_inputs(&mut self, values: &[(HilChannel, f64)]) {
        for &(ch, v) in values {
            self.set_input(ch, v);
        }
    }

    /// Return `true` if a given input channel has been set.
    pub fn has_input(&self, ch: HilChannel) -> bool {
        self.inputs.contains_key(&ch)
    }

    // ── Simple physics model ──────────────────────────────────────────────────

    fn update_physics(&mut self, dt: f64) {
        let throttle = self
            .filtered_inputs
            .get(&HilChannel::ThrottlePedal)
            .copied()
            .unwrap_or(0.0);
        let brake = self
            .filtered_inputs
            .get(&HilChannel::BrakePedal)
            .copied()
            .unwrap_or(0.0);
        let steering = self
            .filtered_inputs
            .get(&HilChannel::SteeringAngle)
            .copied()
            .unwrap_or(0.0);

        // Simplified longitudinal model: F = mass * a
        let mass = 1500.0_f64; // kg
        let drive_force = throttle * 8000.0; // max ~8 kN tractive force
        let brake_force = brake * 15000.0; // max ~15 kN brake force
        let drag = 0.5 * 1.225 * 0.3 * 2.5 * self.speed * self.speed; // ½ρ Cd A v²
        let rolling = 0.015 * mass * 9.81;
        let net_force = drive_force - brake_force - drag - rolling.copysign(self.speed);
        let accel = net_force / mass;

        self.speed += accel * dt;
        self.speed = self.speed.max(0.0); // no reverse (simplified)

        // Yaw rate from Ackermann geometry (simplified)
        let wheelbase = 2.7;
        self.yaw_rate = if self.speed.abs() > 0.1 {
            self.speed * steering.tan() / wheelbase
        } else {
            0.0
        };

        // Engine RPM (simplified: directly coupled to wheel speed with gear ratio 5:1)
        let wheel_radius = 0.33;
        let gear_ratio = 5.0;
        let wheel_omega = self.speed / wheel_radius;
        self.engine_rpm =
            (wheel_omega * gear_ratio * 60.0 / (2.0 * std::f64::consts::PI)).max(800.0); // idle

        // Write outputs
        let wheel_omega_rad = self.speed / wheel_radius;
        *self.outputs.entry(HilChannel::WheelSpeedFl).or_default() = wheel_omega_rad;
        *self.outputs.entry(HilChannel::WheelSpeedFr).or_default() = wheel_omega_rad;
        *self.outputs.entry(HilChannel::WheelSpeedRl).or_default() = wheel_omega_rad;
        *self.outputs.entry(HilChannel::WheelSpeedRr).or_default() = wheel_omega_rad;
        *self
            .outputs
            .entry(HilChannel::AccelLongitudinal)
            .or_default() = accel;
        *self.outputs.entry(HilChannel::AccelLateral).or_default() = self.speed * self.yaw_rate;
        *self.outputs.entry(HilChannel::YawRate).or_default() = self.yaw_rate;
        *self.outputs.entry(HilChannel::VehicleSpeed).or_default() = self.speed;
        *self.outputs.entry(HilChannel::EngineRpm).or_default() = self.engine_rpm;

        // Suspension deflection (simplified spring model)
        let susp_static = mass * 9.81 / 4.0 / 20_000.0; // static deflection with k=20 kN/m
        let susp_lat_transfer = self.outputs[&HilChannel::AccelLateral] * mass / 4.0 / 20_000.0;
        *self.outputs.entry(HilChannel::SuspDeflFl).or_default() = susp_static + susp_lat_transfer;
        *self.outputs.entry(HilChannel::SuspDeflFr).or_default() = susp_static - susp_lat_transfer;
        *self.outputs.entry(HilChannel::SuspDeflRl).or_default() = susp_static + susp_lat_transfer;
        *self.outputs.entry(HilChannel::SuspDeflRr).or_default() = susp_static - susp_lat_transfer;
    }
}

impl HilInterface for SimHilBridge {
    fn set_input(&mut self, channel: HilChannel, value: f64) {
        self.inputs.insert(channel, value);
    }

    fn get_output(&self, channel: HilChannel) -> Option<f64> {
        self.outputs.get(&channel).copied()
    }

    fn tick(&mut self, dt: f64) -> Result<(), HilError> {
        // Apply a first-order IIR low-pass filter to each raw input channel and
        // store the result in `filtered_inputs` for use by the physics model.
        // The filter coefficient α is derived from the bilinear transform:
        //   α = ω_c·dt / (1 + ω_c·dt),  ω_c = 2π·f_c
        // The update rule is:  y_n = α·x_n + (1−α)·y_{n−1}
        // Raw `inputs` are never modified so `set_input` always takes effect
        // on the next tick without being corrupted by repeated filtering.
        if self.config.input_filter_enabled {
            let omega_c = 2.0 * std::f64::consts::PI * self.config.input_filter_hz;
            let alpha = (omega_c * dt) / (1.0 + omega_c * dt);
            let one_minus_alpha = 1.0 - alpha;
            for (&ch, &raw) in &self.inputs {
                let prev = self.filtered_inputs.get(&ch).copied().unwrap_or(raw);
                self.filtered_inputs
                    .insert(ch, alpha * raw + one_minus_alpha * prev);
            }
        } else {
            // No filtering: filtered values mirror raw inputs directly.
            for (&ch, &raw) in &self.inputs {
                self.filtered_inputs.insert(ch, raw);
            }
        }

        self.update_physics(dt);
        self.sim_time += dt;
        self.step_count += 1;
        self.timing.record(dt);

        Ok(())
    }

    fn available_outputs(&self) -> Vec<HilChannel> {
        Self::all_output_channels()
    }

    fn step_count(&self) -> u64 {
        self.step_count
    }

    fn reset(&mut self) {
        self.inputs.clear();
        self.filtered_inputs.clear();
        for v in self.outputs.values_mut() {
            *v = 0.0;
        }
        self.step_count = 0;
        self.sim_time = 0.0;
        self.speed = 0.0;
        self.yaw_rate = 0.0;
        self.engine_rpm = 800.0;
        self.timing = HilTimingStats::default();
    }
}

// ── HilSignalLogger ───────────────────────────────────────────────────────────

/// Records a time-series log of HiL channel values for post-processing.
#[derive(Debug, Default)]
pub struct HilSignalLogger {
    /// Logged (time, channel, value) triples.
    pub log: Vec<(f64, HilChannel, f64)>,
    /// Maximum entries before old entries are dropped (0 = unlimited).
    pub max_entries: usize,
}

impl HilSignalLogger {
    /// Create a logger with optional maximum capacity.
    pub fn new(max_entries: usize) -> Self {
        Self {
            log: Vec::new(),
            max_entries,
        }
    }

    /// Snapshot all current outputs from a bridge into the log at time `t`.
    pub fn snapshot<H: HilInterface>(&mut self, bridge: &H, t: f64) {
        for ch in bridge.available_outputs() {
            if let Some(v) = bridge.get_output(ch) {
                if self.max_entries > 0 && self.log.len() >= self.max_entries {
                    self.log.remove(0);
                }
                self.log.push((t, ch, v));
            }
        }
    }

    /// Return all logged values for a specific channel, sorted by time.
    pub fn channel_history(&self, ch: HilChannel) -> Vec<(f64, f64)> {
        self.log
            .iter()
            .filter(|(_, c, _)| *c == ch)
            .map(|(t, _, v)| (*t, *v))
            .collect()
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sim_hil_bridge_ticks_and_increments_step() {
        let mut bridge = SimHilBridge::new(HilConfig::default());
        bridge.tick(1e-3).unwrap();
        assert_eq!(bridge.step_count(), 1);
        assert!((bridge.sim_time() - 1e-3).abs() < 1e-12);
    }

    #[test]
    fn throttle_increases_speed() {
        let mut bridge = SimHilBridge::new(HilConfig::default());
        bridge.set_input(HilChannel::ThrottlePedal, 1.0);
        for _ in 0..1000 {
            bridge.tick(1e-3).unwrap();
        }
        let speed = bridge.get_output(HilChannel::VehicleSpeed).unwrap();
        assert!(
            speed > 0.1,
            "speed should increase with throttle, got {speed}"
        );
    }

    #[test]
    fn brake_stops_vehicle() {
        let mut bridge = SimHilBridge::new(HilConfig::default());
        // Accelerate first
        bridge.set_input(HilChannel::ThrottlePedal, 1.0);
        for _ in 0..500 {
            bridge.tick(1e-3).unwrap();
        }
        // Then brake
        bridge.set_input(HilChannel::ThrottlePedal, 0.0);
        bridge.set_input(HilChannel::BrakePedal, 1.0);
        for _ in 0..2000 {
            bridge.tick(1e-3).unwrap();
        }
        let speed = bridge.get_output(HilChannel::VehicleSpeed).unwrap();
        assert!(
            speed < 0.01,
            "vehicle should stop after braking, speed = {speed}"
        );
    }

    #[test]
    fn reset_clears_state() {
        let mut bridge = SimHilBridge::new(HilConfig::default());
        bridge.set_input(HilChannel::ThrottlePedal, 1.0);
        for _ in 0..100 {
            bridge.tick(1e-3).unwrap();
        }
        bridge.reset();
        assert_eq!(bridge.step_count(), 0);
        let speed = bridge.get_output(HilChannel::VehicleSpeed).unwrap();
        assert!(
            speed.abs() < 1e-10,
            "speed should be 0 after reset, got {speed}"
        );
    }

    #[test]
    fn signal_logger_records_snapshots() {
        let mut bridge = SimHilBridge::new(HilConfig::default());
        let mut logger = HilSignalLogger::new(0);
        bridge.set_input(HilChannel::ThrottlePedal, 0.5);
        bridge.tick(1e-3).unwrap();
        logger.snapshot(&bridge, bridge.sim_time());
        bridge.tick(1e-3).unwrap();
        logger.snapshot(&bridge, bridge.sim_time());
        let hist = logger.channel_history(HilChannel::VehicleSpeed);
        assert_eq!(hist.len(), 2);
        assert!(
            hist[1].1 >= hist[0].1,
            "speed should not decrease with constant throttle"
        );
    }
}
