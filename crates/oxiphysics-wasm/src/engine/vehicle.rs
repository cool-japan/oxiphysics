// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmVehicleSim` — simple arcade-style vehicle simulation.

use wasm_bindgen::prelude::*;

/// State snapshot of the vehicle.
///
/// Exposed to JavaScript via individual getter methods on `WasmVehicleSim`;
/// the raw `[f64; 3]` arrays are not directly wasm-bindgen compatible.
#[derive(Debug, Clone)]
pub struct VehicleState {
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Heading angle in radians (rotation about Y axis).
    pub heading: f64,
    /// Forward speed (m/s).
    pub speed: f64,
    /// Current gear (1–6).
    pub gear: u32,
    /// Engine RPM.
    pub rpm: f64,
}

/// Simple vehicle simulation suitable for a driving game demo.
///
/// Models: engine torque → wheel force → linear acceleration.
/// Uses a simple Euler integrator. Suspension and tyre slip are abstracted.
#[wasm_bindgen]
#[derive(Debug)]
pub struct WasmVehicleSim {
    position: [f64; 3],
    velocity: [f64; 3],
    heading: f64,
    /// Throttle input 0..1.
    throttle: f64,
    /// Brake input 0..1.
    brake: f64,
    /// Steering input -1..1 (negative = left).
    steering: f64,
    /// Current gear.
    gear: u32,
    /// Engine RPM.
    rpm: f64,
    /// Vehicle mass (kg).
    mass: f64,
    /// Maximum engine torque (N·m).
    max_torque: f64,
    /// Wheel radius (m).
    wheel_radius: f64,
    /// Final drive ratio (engine rpm / wheel rpm).
    final_drive: f64,
    /// Gear ratios for gears 1–6.
    gear_ratios: [f64; 6],
    /// Drag coefficient * frontal area (C_d * A).
    drag_coeff: f64,
    /// Accumulated time.
    time: f64,
}

#[wasm_bindgen]
impl WasmVehicleSim {
    /// Create a car-like vehicle with sensible defaults.
    pub fn new() -> Self {
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            heading: 0.0,
            throttle: 0.0,
            brake: 0.0,
            steering: 0.0,
            gear: 1,
            rpm: 1000.0,
            mass: 1400.0,
            max_torque: 300.0,
            wheel_radius: 0.32,
            final_drive: 3.73,
            gear_ratios: [3.5, 2.0, 1.4, 1.0, 0.75, 0.6],
            drag_coeff: 0.3,
            time: 0.0,
        }
    }

    /// Set throttle (clamped 0..1).
    pub fn set_throttle(&mut self, t: f64) {
        self.throttle = t.clamp(0.0, 1.0);
    }

    /// Set brake (clamped 0..1).
    pub fn set_brake(&mut self, b: f64) {
        self.brake = b.clamp(0.0, 1.0);
    }

    /// Set steering (clamped -1..1).
    pub fn set_steering(&mut self, s: f64) {
        self.steering = s.clamp(-1.0, 1.0);
    }

    /// Set gear (clamped 1..6).
    pub fn set_gear(&mut self, g: u32) {
        self.gear = g.clamp(1, 6);
    }

    /// Set vehicle position.
    pub fn set_position(&mut self, x: f64, y: f64, z: f64) {
        self.position = [x, y, z];
    }

    /// Get heading in radians.
    pub fn heading(&self) -> f64 {
        self.heading
    }

    /// Get current gear.
    pub fn gear(&self) -> u32 {
        self.gear
    }

    /// Get engine RPM.
    pub fn rpm(&self) -> f64 {
        self.rpm
    }

    /// Get forward speed (m/s).
    pub fn speed(&self) -> f64 {
        let vx = self.velocity[0];
        let vz = self.velocity[2];
        (vx * vx + vz * vz).sqrt()
    }

    /// Get position as `Vec<f64>` of `[x, y, z]` (JS-compatible wrapper).
    pub fn get_position_js(&self) -> Vec<f64> {
        self.position.to_vec()
    }

    /// Get velocity as `Vec<f64>` of `[vx, vy, vz]` (JS-compatible wrapper).
    pub fn get_velocity_js(&self) -> Vec<f64> {
        self.velocity.to_vec()
    }

    /// Get full state as a flat `Vec<f64>`:
    /// `[px, py, pz, vx, vy, vz, heading, speed, gear_as_f64, rpm]`.
    pub fn get_state_flat(&self) -> Vec<f64> {
        let speed = self.speed();
        vec![
            self.position[0],
            self.position[1],
            self.position[2],
            self.velocity[0],
            self.velocity[1],
            self.velocity[2],
            self.heading,
            speed,
            f64::from(self.gear),
            self.rpm,
        ]
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }

        let gr = self.gear_ratios[(self.gear - 1) as usize];
        let total_ratio = gr * self.final_drive;

        // Wheel angular velocity from vehicle speed
        let speed =
            (self.velocity[0] * self.velocity[0] + self.velocity[2] * self.velocity[2]).sqrt();
        let wheel_rpm = speed / self.wheel_radius * 60.0 / (2.0 * std::f64::consts::PI);
        self.rpm = (wheel_rpm * total_ratio).clamp(800.0, 8000.0);

        // Engine torque (simple map: peak at mid-RPM).
        // A minimum idle torque is added so the vehicle can pull away from rest.
        let rpm_frac = (self.rpm - 800.0) / 7200.0; // 0..1
        let torque_factor = (4.0 * rpm_frac * (1.0 - rpm_frac)).clamp(0.0, 1.0);
        // Minimum idle factor of 0.2 when throttle is applied, so tractive
        // force exists even at stall RPM.
        let effective_factor = torque_factor.max(if self.throttle > 0.05 { 0.2 } else { 0.0 });
        let engine_torque = self.max_torque * effective_factor * self.throttle;

        // Tractive force
        let tractive_force = engine_torque * total_ratio / self.wheel_radius;

        // Braking force
        let brake_force = self.brake * self.mass * 9.81 * 0.8; // 80% of weight

        // Aero drag
        let drag = self.drag_coeff * 1.225 * speed * speed * 0.5;

        // Net force along heading
        let heading_sin = self.heading.sin();
        let heading_cos = self.heading.cos();
        let net_forward = tractive_force - brake_force - drag;
        let ax = net_forward / self.mass * heading_sin;
        let az = net_forward / self.mass * heading_cos;

        // Integrate velocity
        self.velocity[0] += ax * dt;
        self.velocity[2] += az * dt;

        // Simple friction model (lateral)
        self.velocity[0] *= (1.0 - 0.05 * dt).max(0.0);
        self.velocity[2] *= (1.0 - 0.05 * dt).max(0.0);

        // Steering → heading rate (yaw = v * steer * factor)
        let yaw_rate = speed * self.steering * 0.5;
        self.heading += yaw_rate * dt;

        // Integrate position
        self.position[0] += self.velocity[0] * dt;
        self.position[2] += self.velocity[2] * dt;

        self.time += dt;
    }

    /// Accumulated simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }
}

impl WasmVehicleSim {
    /// Get the current state snapshot (Rust-only; use `get_state_flat` from JS).
    pub fn get_state(&self) -> VehicleState {
        let vx = self.velocity[0];
        let vz = self.velocity[2];
        let speed = (vx * vx + vz * vz).sqrt();
        VehicleState {
            position: self.position,
            velocity: self.velocity,
            heading: self.heading,
            speed,
            gear: self.gear,
            rpm: self.rpm,
        }
    }

    /// Get position as `[x, y, z]` (Rust-only).
    pub fn get_position(&self) -> [f64; 3] {
        self.position
    }

    /// Get velocity as `[vx, vy, vz]` (Rust-only).
    pub fn get_velocity(&self) -> [f64; 3] {
        self.velocity
    }
}

impl Default for WasmVehicleSim {
    fn default() -> Self {
        Self::new()
    }
}
