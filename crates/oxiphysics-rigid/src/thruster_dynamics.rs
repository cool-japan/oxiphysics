// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thruster and propulsion dynamics for underwater/aerial vehicles.
//!
//! Provides first-order thruster models, propeller coefficients, thrust allocation,
//! and a simplified brushless motor electrical model.
//!
//! # Overview
//!
//! - [`Thruster`] — first-order thruster with saturation and response time.
//! - [`ThrusterArray`] — collection of thrusters with force-torque summation.
//! - [`BrushlessMotor`] — electrical brushless-DC motor model.
//! - [`kt_kq_polynomial`] — propeller Kt/Kq as polynomial in advance ratio J.
//! - [`propeller_thrust`] / [`propeller_torque`] — thrust and torque from Kt/Kq.
//! - [`advance_ratio`] — dimensionless advance ratio J.
//! - [`propeller_efficiency`] — open-water propeller efficiency.
//! - [`thrust_allocation_pseudoinverse`] — minimum-norm wrench allocation.
//! - [`required_battery_capacity`] — required battery capacity in Ah.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Thruster
// ─────────────────────────────────────────────────────────────────────────────

/// First-order thruster model with saturation.
#[derive(Debug, Clone)]
pub struct Thruster {
    /// Maximum achievable thrust \[N\].
    pub max_thrust: f64,
    /// Maximum achievable torque \[N·m\].
    pub max_torque: f64,
    /// First-order time constant \[s\].
    pub response_time: f64,
    /// Steady-state efficiency (0..1].
    pub efficiency: f64,
    /// Current thrust output \[N\].
    pub current_thrust: f64,
    /// Current shaft speed \[rev/s\].
    pub current_rpm: f64,
}

impl Thruster {
    /// Create a new thruster.
    ///
    /// * `max_thrust` — saturation thrust \[N\]
    /// * `response_time` — first-order time constant τ \[s\]
    /// * `efficiency` — propulsive efficiency (clamped to (0, 1])
    pub fn new(max_thrust: f64, response_time: f64, efficiency: f64) -> Self {
        Self {
            max_thrust,
            max_torque: max_thrust * 0.1, // heuristic: torque ~ 10 % of thrust * 1 m
            response_time: response_time.max(1e-6),
            efficiency: efficiency.clamp(1e-6, 1.0),
            current_thrust: 0.0,
            current_rpm: 0.0,
        }
    }

    /// Set the desired thrust command \[N\] (clamped to ±max_thrust).
    pub fn set_command(&mut self, cmd: f64) {
        // Stored as desired; step() performs the first-order lag.
        let _clamped = cmd.clamp(-self.max_thrust, self.max_thrust);
        // We store command in current_rpm as a proxy (sign + magnitude).
        self.current_rpm = _clamped;
    }

    /// Advance the thruster state by one time step `dt` \[s\].
    ///
    /// Applies first-order lag: `T_dot = (T_cmd − T) / τ`.
    pub fn step(&mut self, dt: f64) {
        let cmd = self.current_rpm; // desired thrust stored here
        let error = cmd - self.current_thrust;
        self.current_thrust += error * dt / self.response_time;
        // Saturate
        self.current_thrust = self.current_thrust.clamp(-self.max_thrust, self.max_thrust);
    }

    /// Current thrust output \[N\].
    pub fn thrust(&self) -> f64 {
        self.current_thrust
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Propeller free functions
// ─────────────────────────────────────────────────────────────────────────────

/// Propeller thrust and torque coefficients as polynomials in advance ratio J.
///
/// Returns `(Kt, Kq)` using a simple two-term Wageningen B-series approximation:
/// `Kt(J) = 0.3 − 0.25 · J`
/// `Kq(J) = 0.04 − 0.03 · J`
///
/// Negative values are clamped to zero (non-propulsive regime).
pub fn kt_kq_polynomial(j: f64) -> (f64, f64) {
    let kt = (0.3 - 0.25 * j).max(0.0);
    let kq = (0.04 - 0.03 * j).max(1e-9);
    (kt, kq)
}

/// Propeller thrust \[N\].
///
/// `T = Kt · ρ · n² · D⁴`
///
/// * `kt` — thrust coefficient \[-\]
/// * `rho` — fluid density \[kg/m³\]
/// * `n` — shaft speed \[rev/s\]
/// * `d` — propeller diameter \[m\]
pub fn propeller_thrust(kt: f64, rho: f64, n: f64, d: f64) -> f64 {
    kt * rho * n * n * d.powi(4)
}

/// Propeller torque \[N·m\].
///
/// `Q = Kq · ρ · n² · D⁵`
pub fn propeller_torque(kq: f64, rho: f64, n: f64, d: f64) -> f64 {
    kq * rho * n * n * d.powi(5)
}

/// Advance ratio J \[-\].
///
/// `J = Va / (n · D)`
///
/// * `va` — advance velocity (water inflow) \[m/s\]
/// * `n` — shaft speed \[rev/s\]
/// * `d` — propeller diameter \[m\]
pub fn advance_ratio(va: f64, n: f64, d: f64) -> f64 {
    if n.abs() < 1e-12 || d.abs() < 1e-12 {
        return 0.0;
    }
    va / (n * d)
}

/// Open-water propeller efficiency \[-\].
///
/// `η = J · Kt / (2π · Kq)`
///
/// Clamped to \[0, 1\].
pub fn propeller_efficiency(j: f64, kt: f64, kq: f64) -> f64 {
    if kq.abs() < 1e-15 {
        return 0.0;
    }
    (j * kt / (2.0 * PI * kq)).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// ThrusterArray
// ─────────────────────────────────────────────────────────────────────────────

/// Array of thrusters with positions and 6-DOF direction vectors.
#[derive(Debug, Clone)]
pub struct ThrusterArray {
    /// Thrusters with (thruster, position \[m\], 6-DOF wrench column).
    pub thrusters: Vec<(Thruster, [f64; 3], [f64; 6])>,
}

impl ThrusterArray {
    /// Create an empty thruster array.
    pub fn new() -> Self {
        Self {
            thrusters: Vec::new(),
        }
    }

    /// Add a thruster with body-frame position and 6-DOF wrench direction.
    pub fn add_thruster(&mut self, thruster: Thruster, pos: [f64; 3], direction: [f64; 6]) {
        self.thrusters.push((thruster, pos, direction));
    }

    /// Number of thrusters in the array.
    pub fn count(&self) -> usize {
        self.thrusters.len()
    }

    /// Total 6-DOF wrench `[Fx, Fy, Fz, Mx, My, Mz]` produced by all thrusters.
    pub fn total_force_torque(&self) -> [f64; 6] {
        let mut result = [0.0_f64; 6];
        for (thruster, _pos, dir) in &self.thrusters {
            let t = thruster.thrust();
            for i in 0..6 {
                result[i] += t * dir[i];
            }
        }
        result
    }
}

impl Default for ThrusterArray {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Thrust allocation
// ─────────────────────────────────────────────────────────────────────────────

/// Minimum-norm thrust allocation via pseudo-inverse (Moore–Penrose).
///
/// Given the thruster configuration matrix `B` (each row is the 6-DOF wrench
/// column for one thruster) and a desired wrench `[Fx, Fy, Fz, Mx, My, Mz]`,
/// returns the commanded thrust vector `u = B^T (B B^T)^{-1} w`.
///
/// This implementation uses a simple explicit 1-DOF approximation for the
/// common case of a single thruster, and a diagonal approximation otherwise.
/// For full generality a proper linear algebra library is recommended.
pub fn thrust_allocation_pseudoinverse(b_matrix: &[[f64; 6]], wrench: [f64; 6]) -> Vec<f64> {
    let n = b_matrix.len();
    if n == 0 {
        return Vec::new();
    }
    // Compute B * B^T (6×6 matrix is not practical without nalgebra here;
    // instead use the diagonal of B^T B scaled approach)
    // Simplified: u_i = (b_i · w) / (b_i · b_i)
    b_matrix
        .iter()
        .map(|row| {
            let numerator: f64 = row.iter().zip(wrench.iter()).map(|(a, b)| a * b).sum();
            let denominator: f64 = row.iter().map(|a| a * a).sum::<f64>().max(1e-12);
            numerator / denominator
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// BrushlessMotor
// ─────────────────────────────────────────────────────────────────────────────

/// Simplified brushless DC motor (electrical + mechanical model).
#[derive(Debug, Clone)]
pub struct BrushlessMotor {
    /// Motor velocity constant \[rpm/V\] or equivalently \[rad/(s·V)\].
    pub kv: f64,
    /// Phase resistance \[Ω\].
    pub resistance: f64,
    /// Phase inductance \[H\].
    pub inductance: f64,
    /// Rotor inertia \[kg·m²\].
    pub inertia: f64,
    /// Angular velocity \[rad/s\].
    pub omega: f64,
    /// Phase current \[A\].
    pub current: f64,
}

impl BrushlessMotor {
    /// Create a new brushless motor.
    ///
    /// * `kv` — velocity constant \[rpm/V\]; internally converted to \[rad/(s·V)\]
    /// * `resistance` — phase resistance \[Ω\]
    pub fn new(kv: f64, resistance: f64) -> Self {
        Self {
            kv: kv * (2.0 * PI / 60.0), // store as rad/(s·V)
            resistance: resistance.max(1e-9),
            inductance: 1e-4, // 0.1 mH default
            inertia: 1e-5,    // 10 μkg·m² default
            omega: 0.0,
            current: 0.0,
        }
    }

    /// Advance the motor model by one time step `dt` \[s\].
    ///
    /// * `voltage` — applied terminal voltage \[V\]
    /// * `load_torque` — external load torque opposing rotation \[N·m\]
    pub fn step(&mut self, voltage: f64, load_torque: f64, dt: f64) {
        // Back-EMF: e = Kv_rad * omega
        let back_emf = self.kv * self.omega;
        // Electrical: L * di/dt = V - R*i - e
        let di_dt =
            (voltage - self.resistance * self.current - back_emf) / self.inductance.max(1e-12);
        self.current += di_dt * dt;
        self.current = self.current.max(0.0); // no regenerative braking in basic model

        // Motor torque: τ_motor = Kt * i  where Kt = 1 / Kv_rad
        let kt_motor = 1.0 / self.kv.max(1e-12);
        let motor_torque = kt_motor * self.current;
        // Mechanical: J * d(omega)/dt = τ_motor - τ_load
        let d_omega = (motor_torque - load_torque.abs()) / self.inertia.max(1e-15);
        self.omega = (self.omega + d_omega * dt).max(0.0);
    }

    /// Electrical input power \[W\].
    pub fn power_in(&self, voltage: f64) -> f64 {
        voltage * self.current.abs()
    }

    /// Mechanical output power \[W\].
    pub fn power_out(&self) -> f64 {
        let kt_motor = 1.0 / self.kv.max(1e-12);
        kt_motor * self.current * self.omega
    }

    /// Motor efficiency \[-\] (clamped to \[0, 1\]).
    pub fn efficiency(&self, voltage: f64) -> f64 {
        let p_in = self.power_in(voltage);
        if p_in < 1e-12 {
            return 0.0;
        }
        (self.power_out() / p_in).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Battery capacity
// ─────────────────────────────────────────────────────────────────────────────

/// Required battery capacity \[Ah\] for a given mission.
///
/// `C = P_shaft / (V · η) · t_hours`
///
/// where `P_shaft = thrust * v_vehicle` is approximated as `thrust / efficiency`
/// if `v_vehicle` is unknown.
///
/// * `thrust` — required thrust \[N\]
/// * `efficiency` — overall system efficiency \[-\]
/// * `duration_hours` — mission duration \[h\]
/// * `voltage` — battery nominal voltage \[V\]
pub fn required_battery_capacity(
    thrust: f64,
    efficiency: f64,
    duration_hours: f64,
    voltage: f64,
) -> f64 {
    let eta = efficiency.clamp(1e-6, 1.0);
    let power_watts = thrust / eta; // P = F / η (simplified)
    let current_amps = power_watts / voltage.max(1e-6);
    current_amps * duration_hours
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── Thruster ─────────────────────────────────────────────────────────

    #[test]
    fn thruster_initial_thrust_zero() {
        let t = Thruster::new(100.0, 0.1, 0.8);
        assert_eq!(t.thrust(), 0.0);
    }

    #[test]
    fn thruster_saturates_at_max() {
        let mut t = Thruster::new(100.0, 0.01, 0.8);
        t.set_command(9999.0);
        // Run many steps to converge
        for _ in 0..10000 {
            t.step(0.01);
        }
        assert!(t.thrust() <= 100.0 + EPS);
    }

    #[test]
    fn thruster_saturates_negative() {
        let mut t = Thruster::new(100.0, 0.01, 0.8);
        t.set_command(-9999.0);
        for _ in 0..10000 {
            t.step(0.01);
        }
        assert!(t.thrust() >= -100.0 - EPS);
    }

    #[test]
    fn thruster_step_moves_toward_command() {
        let mut t = Thruster::new(100.0, 0.1, 0.8);
        t.set_command(50.0);
        t.step(0.05); // half time-constant
        assert!(t.thrust() > 0.0);
        assert!(t.thrust() < 50.0);
    }

    #[test]
    fn thruster_converges_to_command() {
        let mut t = Thruster::new(100.0, 0.1, 0.8);
        t.set_command(40.0);
        for _ in 0..10000 {
            t.step(0.001);
        }
        assert!((t.thrust() - 40.0).abs() < 0.01);
    }

    #[test]
    fn thruster_zero_command_stays_zero() {
        let mut t = Thruster::new(100.0, 0.1, 0.8);
        t.set_command(0.0);
        for _ in 0..100 {
            t.step(0.01);
        }
        assert!(t.thrust().abs() < EPS);
    }

    // ── advance_ratio ────────────────────────────────────────────────────

    #[test]
    fn advance_ratio_zero_at_zero_speed() {
        assert_eq!(advance_ratio(0.0, 10.0, 0.3), 0.0);
    }

    #[test]
    fn advance_ratio_positive() {
        let j = advance_ratio(2.0, 10.0, 0.3);
        assert!(j > 0.0);
    }

    #[test]
    fn advance_ratio_zero_n() {
        // n = 0 → should return 0 without panic
        assert_eq!(advance_ratio(5.0, 0.0, 0.3), 0.0);
    }

    #[test]
    fn advance_ratio_formula() {
        // J = 3 / (10 * 0.3) = 1.0
        let j = advance_ratio(3.0, 10.0, 0.3);
        assert!((j - 1.0).abs() < EPS);
    }

    // ── kt_kq_polynomial ─────────────────────────────────────────────────

    #[test]
    fn kt_kq_at_j_zero() {
        let (kt, kq) = kt_kq_polynomial(0.0);
        assert!((kt - 0.3).abs() < EPS);
        assert!((kq - 0.04).abs() < EPS);
    }

    #[test]
    fn kt_kq_positive() {
        let (kt, kq) = kt_kq_polynomial(0.5);
        assert!(kt >= 0.0);
        assert!(kq >= 0.0);
    }

    #[test]
    fn kt_kq_decreases_with_j() {
        let (kt0, _) = kt_kq_polynomial(0.0);
        let (kt1, _) = kt_kq_polynomial(0.5);
        assert!(kt1 < kt0);
    }

    // ── propeller_thrust ─────────────────────────────────────────────────

    #[test]
    fn propeller_thrust_positive() {
        let t = propeller_thrust(0.3, 1000.0, 10.0, 0.3);
        assert!(t > 0.0);
    }

    #[test]
    fn propeller_thrust_scales_with_n_squared() {
        let t1 = propeller_thrust(0.3, 1000.0, 10.0, 0.3);
        let t2 = propeller_thrust(0.3, 1000.0, 20.0, 0.3);
        assert!((t2 / t1 - 4.0).abs() < 1e-6);
    }

    // ── propeller_efficiency ─────────────────────────────────────────────

    #[test]
    fn propeller_efficiency_in_zero_one() {
        let (kt, kq) = kt_kq_polynomial(0.5);
        let eta = propeller_efficiency(0.5, kt, kq);
        assert!((0.0..=1.0).contains(&eta));
    }

    #[test]
    fn propeller_efficiency_zero_at_j_zero() {
        let (kt, kq) = kt_kq_polynomial(0.0);
        let eta = propeller_efficiency(0.0, kt, kq);
        assert_eq!(eta, 0.0);
    }

    // ── ThrusterArray ────────────────────────────────────────────────────

    #[test]
    fn thruster_array_count() {
        let mut arr = ThrusterArray::new();
        arr.add_thruster(
            Thruster::new(50.0, 0.1, 0.8),
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        );
        arr.add_thruster(
            Thruster::new(50.0, 0.1, 0.8),
            [0.0, 0.0, -1.0],
            [0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        );
        assert_eq!(arr.count(), 2);
    }

    #[test]
    fn thruster_array_total_zero_at_start() {
        let arr = ThrusterArray::new();
        let w = arr.total_force_torque();
        assert!(w.iter().all(|x| x.abs() < EPS));
    }

    #[test]
    fn thruster_array_total_force_after_step() {
        let mut arr = ThrusterArray::new();
        let mut thr = Thruster::new(50.0, 0.01, 0.8);
        thr.set_command(50.0);
        for _ in 0..1000 {
            thr.step(0.01);
        }
        arr.add_thruster(thr, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let w = arr.total_force_torque();
        assert!(w[0] > 40.0);
    }

    // ── thrust_allocation_pseudoinverse ──────────────────────────────────

    #[test]
    fn thrust_allocation_empty_returns_empty() {
        let result = thrust_allocation_pseudoinverse(&[], [0.0; 6]);
        assert!(result.is_empty());
    }

    #[test]
    fn thrust_allocation_single_thruster() {
        let b: &[[f64; 6]] = &[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
        let wrench = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let u = thrust_allocation_pseudoinverse(b, wrench);
        assert_eq!(u.len(), 1);
        assert!((u[0] - 10.0).abs() < EPS);
    }

    // ── BrushlessMotor ───────────────────────────────────────────────────

    #[test]
    fn motor_initial_state_zero() {
        let m = BrushlessMotor::new(1000.0, 0.1);
        assert_eq!(m.omega, 0.0);
        assert_eq!(m.current, 0.0);
    }

    #[test]
    fn motor_spins_up_under_voltage() {
        let mut m = BrushlessMotor::new(1000.0, 0.1);
        for _ in 0..1000 {
            m.step(12.0, 0.0, 0.001);
        }
        assert!(m.omega > 0.0);
    }

    #[test]
    fn motor_power_in_non_negative() {
        let mut m = BrushlessMotor::new(1000.0, 0.1);
        m.step(12.0, 0.0, 0.001);
        assert!(m.power_in(12.0) >= 0.0);
    }

    #[test]
    fn motor_efficiency_in_zero_one() {
        let mut m = BrushlessMotor::new(1000.0, 0.1);
        for _ in 0..500 {
            m.step(12.0, 0.01, 0.001);
        }
        let eta = m.efficiency(12.0);
        assert!((0.0..=1.0).contains(&eta));
    }

    #[test]
    fn motor_power_in_ge_power_out() {
        let mut m = BrushlessMotor::new(1000.0, 0.1);
        for _ in 0..500 {
            m.step(12.0, 0.01, 0.001);
        }
        assert!(m.power_in(12.0) >= m.power_out() - EPS);
    }

    // ── required_battery_capacity ────────────────────────────────────────

    #[test]
    fn battery_capacity_positive() {
        let cap = required_battery_capacity(100.0, 0.8, 1.0, 14.8);
        assert!(cap > 0.0);
    }

    #[test]
    fn battery_capacity_scales_with_duration() {
        let c1 = required_battery_capacity(100.0, 0.8, 1.0, 14.8);
        let c2 = required_battery_capacity(100.0, 0.8, 2.0, 14.8);
        assert!((c2 - 2.0 * c1).abs() < 1e-6);
    }

    #[test]
    fn battery_capacity_increases_with_thrust() {
        let c1 = required_battery_capacity(100.0, 0.8, 1.0, 14.8);
        let c2 = required_battery_capacity(200.0, 0.8, 1.0, 14.8);
        assert!(c2 > c1);
    }

    #[test]
    fn battery_capacity_decreases_with_efficiency() {
        let c1 = required_battery_capacity(100.0, 0.5, 1.0, 14.8);
        let c2 = required_battery_capacity(100.0, 0.9, 1.0, 14.8);
        assert!(c1 > c2);
    }
}
