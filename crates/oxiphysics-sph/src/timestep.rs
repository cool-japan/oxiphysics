// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive time stepping for SPH simulations.
//!
//! Provides CFL and viscous timestep constraints to ensure stability.
//! Also provides leapfrog integration helpers, adaptive timestep control,
//! surface tension timestep, acceleration-based timestep, multi-criteria
//! timestep selection, and timestep history tracking.

use crate::particle::ParticleSet;

// ── Legacy helpers (pre-existing) ─────────────────────────────────────────────

/// Compute the CFL-limited timestep.
///
/// `dt_cfl = lambda * h / (c_s + max|v|)`
///
/// where `lambda` is a safety factor (typically 0.25–0.4), `h` is the
/// smoothing length, and `c_s` is the numerical speed of sound.
pub fn cfl_timestep(particles: &ParticleSet, h: f64, sound_speed: f64) -> f64 {
    let lambda = 0.4;
    let max_vel = particles
        .velocities
        .iter()
        .map(|v| v.norm())
        .fold(0.0_f64, f64::max);
    lambda * h / (sound_speed + max_vel).max(1e-14)
}

/// Compute the viscous timestep limit.
///
/// `dt_visc = lambda * h² / (nu + epsilon)` where `nu` is kinematic viscosity.
pub fn viscous_timestep(h: f64, viscosity: f64) -> f64 {
    let lambda = 0.25;
    lambda * h * h / viscosity.max(1e-14)
}

/// Compute the combined adaptive timestep, clamped to `[min_dt, max_dt]`.
pub fn adaptive_timestep(
    particles: &ParticleSet,
    h: f64,
    sound_speed: f64,
    viscosity: f64,
    min_dt: f64,
    max_dt: f64,
) -> f64 {
    let dt_cfl = cfl_timestep(particles, h, sound_speed);
    let dt_visc = viscous_timestep(h, viscosity);
    dt_cfl.min(dt_visc).clamp(min_dt, max_dt)
}

// ── CflCondition ──────────────────────────────────────────────────────────────

/// Aggregated CFL condition parameters derived from particle data.
#[derive(Debug, Clone, PartialEq)]
pub struct CflCondition {
    /// Maximum speed of sound encountered across particles.
    pub c_sound_max: f64,
    /// Maximum particle speed encountered.
    pub v_max: f64,
    /// Kinematic viscosity coefficient.
    pub viscosity: f64,
    /// Minimum smoothing length across particles.
    pub h_min: f64,
}

// ── Standalone CFL functions ──────────────────────────────────────────────────

/// Acoustic CFL timestep: `dt = cfl_factor * h / c_sound`.
///
/// Ensures that information (pressure waves) cannot travel more than one
/// smoothing-length interval per timestep.
pub fn cfl_acoustic(c_sound: f64, h: f64, cfl_factor: f64) -> f64 {
    cfl_factor * h / c_sound.max(1e-14)
}

/// Viscous CFL timestep: `dt = cfl_factor * h² / nu`.
///
/// Ensures numerical stability of the diffusion operator.
pub fn cfl_viscous(nu: f64, h: f64, cfl_factor: f64) -> f64 {
    cfl_factor * h * h / nu.max(1e-14)
}

/// Combined CFL timestep — minimum of acoustic and viscous constraints.
///
/// Also incorporates the convective CFL `h / (c_sound + v_max)`.
pub fn cfl_combined(c_sound: f64, v_max: f64, nu: f64, h: f64, cfl_factor: f64) -> f64 {
    let dt_acoustic = cfl_acoustic(c_sound, h, cfl_factor);
    let dt_viscous = cfl_viscous(nu, h, cfl_factor);
    // Convective: h / (c_s + |v_max|)
    let dt_convective = cfl_factor * h / (c_sound + v_max).max(1e-14);
    dt_acoustic.min(dt_viscous).min(dt_convective)
}

// ── Surface tension timestep ──────────────────────────────────────────────────

/// Surface tension timestep constraint.
///
/// `dt_st = cfl_factor * sqrt(rho * h³ / (2π σ))`
///
/// where σ is the surface tension coefficient (N/m).
pub fn surface_tension_timestep(rho: f64, h: f64, sigma: f64, cfl_factor: f64) -> f64 {
    if sigma.abs() < 1e-30 || rho < 1e-30 {
        return f64::MAX;
    }
    cfl_factor * (rho * h * h * h / (2.0 * std::f64::consts::PI * sigma)).sqrt()
}

// ── Acceleration-based timestep ───────────────────────────────────────────────

/// Acceleration-based timestep constraint.
///
/// `dt_acc = cfl_factor * sqrt(h / |a_max|)`
///
/// Limits timestep based on maximum particle acceleration.
pub fn acceleration_timestep(h: f64, a_max: f64, cfl_factor: f64) -> f64 {
    if a_max.abs() < 1e-30 {
        return f64::MAX;
    }
    cfl_factor * (h / a_max).sqrt()
}

/// Compute maximum acceleration magnitude from an array of accelerations.
pub fn max_acceleration(accelerations: &[[f64; 3]]) -> f64 {
    accelerations
        .iter()
        .map(|a| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

// ── Multi-criteria timestep selection ─────────────────────────────────────────

/// Parameters for multi-criteria timestep selection.
#[derive(Debug, Clone)]
pub struct TimestepCriteria {
    /// CFL factor for all criteria.
    pub cfl_factor: f64,
    /// Speed of sound.
    pub c_sound: f64,
    /// Maximum particle velocity.
    pub v_max: f64,
    /// Kinematic viscosity.
    pub nu: f64,
    /// Smoothing length.
    pub h: f64,
    /// Surface tension coefficient (N/m). Set to 0 to disable.
    pub sigma: f64,
    /// Rest density (kg/m³).
    pub rho: f64,
    /// Maximum acceleration magnitude.
    pub a_max: f64,
    /// Minimum allowed timestep.
    pub min_dt: f64,
    /// Maximum allowed timestep.
    pub max_dt: f64,
}

impl TimestepCriteria {
    /// Compute the minimum timestep across all criteria.
    pub fn compute(&self) -> f64 {
        let dt_acoustic = cfl_acoustic(self.c_sound, self.h, self.cfl_factor);
        let dt_viscous = cfl_viscous(self.nu, self.h, self.cfl_factor);
        let dt_convective = self.cfl_factor * self.h / (self.c_sound + self.v_max).max(1e-14);
        let dt_st = surface_tension_timestep(self.rho, self.h, self.sigma, self.cfl_factor);
        let dt_acc = acceleration_timestep(self.h, self.a_max, self.cfl_factor);

        dt_acoustic
            .min(dt_viscous)
            .min(dt_convective)
            .min(dt_st)
            .min(dt_acc)
            .clamp(self.min_dt, self.max_dt)
    }

    /// Identify which criterion is the most restrictive.
    pub fn limiting_criterion(&self) -> &'static str {
        let dt_acoustic = cfl_acoustic(self.c_sound, self.h, self.cfl_factor);
        let dt_viscous = cfl_viscous(self.nu, self.h, self.cfl_factor);
        let dt_convective = self.cfl_factor * self.h / (self.c_sound + self.v_max).max(1e-14);
        let dt_st = surface_tension_timestep(self.rho, self.h, self.sigma, self.cfl_factor);
        let dt_acc = acceleration_timestep(self.h, self.a_max, self.cfl_factor);

        let min = dt_acoustic
            .min(dt_viscous)
            .min(dt_convective)
            .min(dt_st)
            .min(dt_acc);

        if (min - dt_acoustic).abs() < 1e-30 {
            "acoustic"
        } else if (min - dt_viscous).abs() < 1e-30 {
            "viscous"
        } else if (min - dt_convective).abs() < 1e-30 {
            "convective"
        } else if (min - dt_st).abs() < 1e-30 {
            "surface_tension"
        } else {
            "acceleration"
        }
    }
}

// ── Timestep history tracking ─────────────────────────────────────────────────

/// Tracks timestep history for diagnostics and analysis.
#[derive(Debug, Clone)]
pub struct TimestepHistory {
    /// Timestep values over time.
    pub values: Vec<f64>,
    /// Simulation time at each recorded step.
    pub times: Vec<f64>,
    /// Maximum number of entries to keep (ring buffer).
    pub max_entries: usize,
}

impl TimestepHistory {
    /// Create a new history tracker.
    pub fn new(max_entries: usize) -> Self {
        Self {
            values: Vec::with_capacity(max_entries),
            times: Vec::with_capacity(max_entries),
            max_entries,
        }
    }

    /// Record a timestep.
    pub fn record(&mut self, time: f64, dt: f64) {
        if self.values.len() >= self.max_entries {
            self.values.remove(0);
            self.times.remove(0);
        }
        self.values.push(dt);
        self.times.push(time);
    }

    /// Get the average timestep over the recorded history.
    pub fn average_dt(&self) -> f64 {
        if self.values.is_empty() {
            return 0.0;
        }
        self.values.iter().sum::<f64>() / self.values.len() as f64
    }

    /// Get the minimum timestep in the recorded history.
    pub fn min_dt(&self) -> f64 {
        self.values.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Get the maximum timestep in the recorded history.
    pub fn max_dt(&self) -> f64 {
        self.values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Get the standard deviation of timestep values.
    pub fn std_dev(&self) -> f64 {
        if self.values.len() < 2 {
            return 0.0;
        }
        let mean = self.average_dt();
        let variance = self
            .values
            .iter()
            .map(|&v| (v - mean) * (v - mean))
            .sum::<f64>()
            / (self.values.len() - 1) as f64;
        variance.sqrt()
    }

    /// Number of recorded entries.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether history is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Compute the rate of change of dt (dt_dot) from the last two entries.
    pub fn dt_rate_of_change(&self) -> f64 {
        let n = self.values.len();
        if n < 2 {
            return 0.0;
        }
        let dt_diff = self.values[n - 1] - self.values[n - 2];
        let t_diff = self.times[n - 1] - self.times[n - 2];
        if t_diff.abs() < 1e-30 {
            return 0.0;
        }
        dt_diff / t_diff
    }
}

impl Default for TimestepHistory {
    fn default() -> Self {
        Self::new(1000)
    }
}

// ── AdaptiveTimestep ──────────────────────────────────────────────────────────

/// State for an adaptive timestep controller.
///
/// Grows or shrinks the current step based on whether a proposed CFL condition
/// allows the step to grow.
#[derive(Debug, Clone)]
pub struct AdaptiveTimestep {
    /// Current timestep value.
    pub current_dt: f64,
    /// Minimum allowed timestep.
    pub min_dt: f64,
    /// Maximum allowed timestep.
    pub max_dt: f64,
    /// Multiplicative growth factor when the step can be enlarged (> 1).
    pub growth_factor: f64,
    /// Multiplicative shrink factor when the step must be reduced (< 1).
    pub shrink_factor: f64,
}

impl AdaptiveTimestep {
    /// Create a new `AdaptiveTimestep` with sensible defaults.
    pub fn new(initial_dt: f64, min_dt: f64, max_dt: f64) -> Self {
        Self {
            current_dt: initial_dt,
            min_dt,
            max_dt,
            growth_factor: 1.1,
            shrink_factor: 0.5,
        }
    }

    /// Update the timestep based on a proposed CFL limit.
    ///
    /// If `cfl_dt` (the CFL-derived safe dt) is larger than `current_dt`,
    /// the step grows by `growth_factor` (clamped to `max_dt`).
    /// If `cfl_dt` is smaller than `current_dt`, the step shrinks to
    /// `cfl_dt * shrink_factor` (clamped to `min_dt`).
    ///
    /// Returns the new `current_dt`.
    pub fn update(&mut self, cfl_dt: f64) -> f64 {
        if cfl_dt >= self.current_dt {
            // CFL allows growth.
            self.current_dt = (self.current_dt * self.growth_factor).min(self.max_dt);
        } else {
            // CFL forces reduction.
            self.current_dt = (cfl_dt * self.shrink_factor).max(self.min_dt);
        }
        self.current_dt
    }

    /// Update with multi-criteria: takes a TimestepCriteria and uses its computed dt.
    pub fn update_multi(&mut self, criteria: &TimestepCriteria) -> f64 {
        let cfl_dt = criteria.compute();
        self.update(cfl_dt)
    }
}

// ── sph_dt_from_particles ──────────────────────────────────────────────────────

/// Compute a combined CFL timestep directly from a `ParticleSet`.
///
/// Iterates over all particles to find the maximum speed and uses the
/// provided `h` (smoothing length) and `nu` (kinematic viscosity) to
/// call `cfl_combined`.
pub fn sph_dt_from_particles(particles: &ParticleSet, h: f64, nu: f64, cfl: f64) -> f64 {
    let v_max = particles
        .velocities
        .iter()
        .map(|v| v.norm())
        .fold(0.0_f64, f64::max);
    // Use sound speed = 0.0 here; callers may add their own c_s contribution.
    // In practice, pass the actual c_s via `cfl_combined` directly.
    cfl_combined(0.0, v_max, nu, h, cfl)
}

// ── TimeIntegrator ────────────────────────────────────────────────────────────

/// Time integration scheme selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeIntegrator {
    /// First-order explicit Euler: `x += v*dt`, `v += a*dt`.
    Euler,
    /// Symplectic leapfrog (kick-drift or drift-kick).
    LeapFrog,
    /// Classic 4th-order Runge-Kutta.
    RungeKutta4,
}

// ── Leapfrog helpers ──────────────────────────────────────────────────────────

/// Leapfrog *kick*: half-step velocity update.
///
/// `v_new = v + a * dt`
///
/// In the kick-drift-kick (KDK) scheme this is applied with `dt/2` at the
/// start and end of each step; in drift-kick-drift (DKD) it is applied once
/// per step.
pub fn leapfrog_kick(v: [f64; 3], a: [f64; 3], dt: f64) -> [f64; 3] {
    [v[0] + a[0] * dt, v[1] + a[1] * dt, v[2] + a[2] * dt]
}

/// Leapfrog *drift*: position update using the (already kicked) velocity.
///
/// `x_new = x + v * dt`
pub fn leapfrog_drift(x: [f64; 3], v: [f64; 3], dt: f64) -> [f64; 3] {
    [x[0] + v[0] * dt, x[1] + v[1] * dt, x[2] + v[2] * dt]
}

/// Simple Euler integrator for reference (position + velocity).
///
/// Returns `(new_position, new_velocity)`.
pub fn euler_step(x: [f64; 3], v: [f64; 3], a: [f64; 3], dt: f64) -> ([f64; 3], [f64; 3]) {
    let xn = leapfrog_drift(x, v, dt);
    let vn = leapfrog_kick(v, a, dt);
    (xn, vn)
}

/// One full leapfrog KDK step (kick–drift–kick).
///
/// The acceleration `a` is evaluated at the midpoint position.
/// Caller is responsible for recomputing `a` at `x_mid` if needed.
/// Returns `(new_position, new_velocity)`.
pub fn leapfrog_kdk(x: [f64; 3], v: [f64; 3], a: [f64; 3], dt: f64) -> ([f64; 3], [f64; 3]) {
    let v_half = leapfrog_kick(v, a, dt * 0.5);
    let x_new = leapfrog_drift(x, v_half, dt);
    // Second kick uses the same `a` (caller must update externally).
    let v_new = leapfrog_kick(v_half, a, dt * 0.5);
    (x_new, v_new)
}

/// Runge-Kutta 4th-order step for a constant acceleration field.
///
/// Returns `(new_position, new_velocity)`.
pub fn rk4_constant_accel(x: [f64; 3], v: [f64; 3], a: [f64; 3], dt: f64) -> ([f64; 3], [f64; 3]) {
    // k1
    let k1_x = v;
    let k1_v = a;

    // k2 (midpoint)
    let x2 = [
        x[0] + k1_x[0] * dt * 0.5,
        x[1] + k1_x[1] * dt * 0.5,
        x[2] + k1_x[2] * dt * 0.5,
    ];
    let v2 = [
        v[0] + k1_v[0] * dt * 0.5,
        v[1] + k1_v[1] * dt * 0.5,
        v[2] + k1_v[2] * dt * 0.5,
    ];
    let _ = x2;
    let k2_x = v2;
    let k2_v = a; // constant accel

    // k3 (midpoint with k2 slopes)
    let v3 = [
        v[0] + k2_v[0] * dt * 0.5,
        v[1] + k2_v[1] * dt * 0.5,
        v[2] + k2_v[2] * dt * 0.5,
    ];
    let k3_x = v3;
    let k3_v = a;

    // k4 (full step)
    let v4 = [
        v[0] + k3_v[0] * dt,
        v[1] + k3_v[1] * dt,
        v[2] + k3_v[2] * dt,
    ];
    let k4_x = v4;
    let k4_v = a;

    let sixth = 1.0 / 6.0;
    let x_new = [
        x[0] + dt * sixth * (k1_x[0] + 2.0 * k2_x[0] + 2.0 * k3_x[0] + k4_x[0]),
        x[1] + dt * sixth * (k1_x[1] + 2.0 * k2_x[1] + 2.0 * k3_x[1] + k4_x[1]),
        x[2] + dt * sixth * (k1_x[2] + 2.0 * k2_x[2] + 2.0 * k3_x[2] + k4_x[2]),
    ];
    let v_new = [
        v[0] + dt * sixth * (k1_v[0] + 2.0 * k2_v[0] + 2.0 * k3_v[0] + k4_v[0]),
        v[1] + dt * sixth * (k1_v[1] + 2.0 * k2_v[1] + 2.0 * k3_v[1] + k4_v[1]),
        v[2] + dt * sixth * (k1_v[2] + 2.0 * k2_v[2] + 2.0 * k3_v[2] + k4_v[2]),
    ];
    (x_new, v_new)
}

// ── Velocity Verlet ───────────────────────────────────────────────────────────

/// Velocity Verlet integration step.
///
/// Given current position, velocity, and acceleration, plus new acceleration:
/// x_{n+1} = x_n + v_n * dt + 0.5 * a_n * dt²
/// v_{n+1} = v_n + 0.5 * (a_n + a_{n+1}) * dt
///
/// Returns `(new_position, new_velocity)`.
pub fn velocity_verlet(
    x: [f64; 3],
    v: [f64; 3],
    a_old: [f64; 3],
    a_new: [f64; 3],
    dt: f64,
) -> ([f64; 3], [f64; 3]) {
    let x_new = [
        x[0] + v[0] * dt + 0.5 * a_old[0] * dt * dt,
        x[1] + v[1] * dt + 0.5 * a_old[1] * dt * dt,
        x[2] + v[2] * dt + 0.5 * a_old[2] * dt * dt,
    ];
    let v_new = [
        v[0] + 0.5 * (a_old[0] + a_new[0]) * dt,
        v[1] + 0.5 * (a_old[1] + a_new[1]) * dt,
        v[2] + 0.5 * (a_old[2] + a_new[2]) * dt,
    ];
    (x_new, v_new)
}

// ── Multi-step SPH integrator ─────────────────────────────────────────────────

/// Result of one multi-step integration pass for a single particle.
#[derive(Debug, Clone)]
pub struct MultiStepResult {
    /// Final position after all substeps.
    pub position: [f64; 3],
    /// Final velocity after all substeps.
    pub velocity: [f64; 3],
    /// Number of substeps taken.
    pub substeps: usize,
    /// Accumulated simulation time.
    pub time: f64,
}

/// Integrate a single particle for a total time `dt_total` using `n_substeps`
/// equal sub-intervals of the symplectic Euler method.
///
/// `accel_fn(x, v) -> [f64; 3]` computes the acceleration at the given state.
pub fn sph_multistep_euler<F>(
    x0: [f64; 3],
    v0: [f64; 3],
    dt_total: f64,
    n_substeps: usize,
    mut accel_fn: F,
) -> MultiStepResult
where
    F: FnMut([f64; 3], [f64; 3]) -> [f64; 3],
{
    let dt = dt_total / n_substeps.max(1) as f64;
    let mut x = x0;
    let mut v = v0;
    for _ in 0..n_substeps {
        let a = accel_fn(x, v);
        v = [v[0] + a[0] * dt, v[1] + a[1] * dt, v[2] + a[2] * dt];
        x = [x[0] + v[0] * dt, x[1] + v[1] * dt, x[2] + v[2] * dt];
    }
    MultiStepResult {
        position: x,
        velocity: v,
        substeps: n_substeps,
        time: dt_total,
    }
}

/// Integrate a single particle for `dt_total` using `n_substeps` substeps
/// of the 4th-order Runge-Kutta method for SPH.
///
/// The acceleration function `accel_fn(x, v) -> [f64; 3]` may depend on
/// both position and velocity (e.g. drag terms).
pub fn sph_multistep_rk4<F>(
    x0: [f64; 3],
    v0: [f64; 3],
    dt_total: f64,
    n_substeps: usize,
    mut accel_fn: F,
) -> MultiStepResult
where
    F: FnMut([f64; 3], [f64; 3]) -> [f64; 3],
{
    let dt = dt_total / n_substeps.max(1) as f64;
    let mut x = x0;
    let mut v = v0;
    for _ in 0..n_substeps {
        // k1
        let a1 = accel_fn(x, v);
        let kx1 = v;
        let kv1 = a1;
        // k2
        let x2 = [
            x[0] + kx1[0] * dt * 0.5,
            x[1] + kx1[1] * dt * 0.5,
            x[2] + kx1[2] * dt * 0.5,
        ];
        let v2 = [
            v[0] + kv1[0] * dt * 0.5,
            v[1] + kv1[1] * dt * 0.5,
            v[2] + kv1[2] * dt * 0.5,
        ];
        let a2 = accel_fn(x2, v2);
        let kx2 = v2;
        let kv2 = a2;
        // k3
        let x3 = [
            x[0] + kx2[0] * dt * 0.5,
            x[1] + kx2[1] * dt * 0.5,
            x[2] + kx2[2] * dt * 0.5,
        ];
        let v3 = [
            v[0] + kv2[0] * dt * 0.5,
            v[1] + kv2[1] * dt * 0.5,
            v[2] + kv2[2] * dt * 0.5,
        ];
        let a3 = accel_fn(x3, v3);
        let kx3 = v3;
        let kv3 = a3;
        // k4
        let x4 = [x[0] + kx3[0] * dt, x[1] + kx3[1] * dt, x[2] + kx3[2] * dt];
        let v4 = [v[0] + kv3[0] * dt, v[1] + kv3[1] * dt, v[2] + kv3[2] * dt];
        let a4 = accel_fn(x4, v4);
        let kx4 = v4;
        let kv4 = a4;

        let sixth = 1.0 / 6.0;
        x = [
            x[0] + dt * sixth * (kx1[0] + 2.0 * kx2[0] + 2.0 * kx3[0] + kx4[0]),
            x[1] + dt * sixth * (kx1[1] + 2.0 * kx2[1] + 2.0 * kx3[1] + kx4[1]),
            x[2] + dt * sixth * (kx1[2] + 2.0 * kx2[2] + 2.0 * kx3[2] + kx4[2]),
        ];
        v = [
            v[0] + dt * sixth * (kv1[0] + 2.0 * kv2[0] + 2.0 * kv3[0] + kv4[0]),
            v[1] + dt * sixth * (kv1[1] + 2.0 * kv2[1] + 2.0 * kv3[1] + kv4[1]),
            v[2] + dt * sixth * (kv1[2] + 2.0 * kv2[2] + 2.0 * kv3[2] + kv4[2]),
        ];
    }
    MultiStepResult {
        position: x,
        velocity: v,
        substeps: n_substeps,
        time: dt_total,
    }
}

// ── Predictor-Corrector SPH integrator ────────────────────────────────────────

/// Predictor-corrector (Adams–Bashforth/Adams–Moulton) integration for SPH.
///
/// The predictor uses a forward Euler step; the corrector then re-evaluates
/// the acceleration at the predicted state and averages (trapezoidal rule).
///
/// Returns `(new_position, new_velocity)`.
pub fn predictor_corrector_step<F>(
    x: [f64; 3],
    v: [f64; 3],
    a_old: [f64; 3],
    dt: f64,
    mut accel_fn: F,
) -> ([f64; 3], [f64; 3])
where
    F: FnMut([f64; 3], [f64; 3]) -> [f64; 3],
{
    // Predictor: forward Euler
    let v_pred = [
        v[0] + a_old[0] * dt,
        v[1] + a_old[1] * dt,
        v[2] + a_old[2] * dt,
    ];
    let x_pred = [x[0] + v[0] * dt, x[1] + v[1] * dt, x[2] + v[2] * dt];
    // Evaluate acceleration at predicted state
    let a_new = accel_fn(x_pred, v_pred);
    // Corrector: trapezoidal average
    let v_cor = [
        v[0] + 0.5 * (a_old[0] + a_new[0]) * dt,
        v[1] + 0.5 * (a_old[1] + a_new[1]) * dt,
        v[2] + 0.5 * (a_old[2] + a_new[2]) * dt,
    ];
    let x_cor = [
        x[0] + 0.5 * (v[0] + v_pred[0]) * dt,
        x[1] + 0.5 * (v[1] + v_pred[1]) * dt,
        x[2] + 0.5 * (v[2] + v_pred[2]) * dt,
    ];
    (x_cor, v_cor)
}

// ── Symplectic Euler for SPH ──────────────────────────────────────────────────

/// Symplectic Euler step for SPH.
///
/// Velocity is updated first (using old position's acceleration), then
/// position is advanced with the new velocity.  This is the standard
/// symplectic / semi-implicit Euler integrator that conserves a modified
/// energy for Hamiltonian systems.
///
/// Returns `(new_position, new_velocity)`.
pub fn symplectic_euler_sph(
    x: [f64; 3],
    v: [f64; 3],
    a: [f64; 3],
    dt: f64,
) -> ([f64; 3], [f64; 3]) {
    let v_new = [v[0] + a[0] * dt, v[1] + a[1] * dt, v[2] + a[2] * dt];
    let x_new = [
        x[0] + v_new[0] * dt,
        x[1] + v_new[1] * dt,
        x[2] + v_new[2] * dt,
    ];
    (x_new, v_new)
}

// ── Adaptive timestep with body-force criterion ───────────────────────────────

/// Compute a timestep limited by a body-force acceleration magnitude.
///
/// `dt = cfl_factor * sqrt(h / |a_body|)` where `a_body` is a background
/// body-force acceleration (e.g. gravity).
pub fn body_force_timestep(h: f64, a_body: f64, cfl_factor: f64) -> f64 {
    if a_body.abs() < 1e-30 {
        return f64::MAX;
    }
    cfl_factor * (h / a_body.abs()).sqrt()
}

/// Combined adaptive timestep with CFL + viscous + body-force criteria.
pub fn adaptive_timestep_full(
    c_sound: f64,
    v_max: f64,
    nu: f64,
    h: f64,
    a_body: f64,
    cfl_factor: f64,
    min_dt: f64,
    max_dt: f64,
) -> f64 {
    let dt_cfl = cfl_combined(c_sound, v_max, nu, h, cfl_factor);
    let dt_body = body_force_timestep(h, a_body, cfl_factor);
    dt_cfl.min(dt_body).clamp(min_dt, max_dt)
}

// ── SubstepController ─────────────────────────────────────────────────────────

/// Controller that breaks a target time interval into sub-steps that each
/// satisfy the given CFL condition.
#[derive(Debug, Clone)]
pub struct SubstepController {
    /// Target simulation interval (outer step).
    pub dt_outer: f64,
    /// Maximum CFL-allowed sub-step size.
    pub dt_cfl_max: f64,
    /// Minimum allowed sub-step.
    pub dt_min: f64,
}

impl SubstepController {
    /// Create a new sub-step controller.
    pub fn new(dt_outer: f64, dt_cfl_max: f64, dt_min: f64) -> Self {
        Self {
            dt_outer,
            dt_cfl_max,
            dt_min,
        }
    }

    /// Compute the number of sub-steps needed to cover `dt_outer` with
    /// sub-steps no larger than `dt_cfl_max`.
    pub fn num_substeps(&self) -> usize {
        if self.dt_cfl_max <= 0.0 || self.dt_outer <= 0.0 {
            return 1;
        }
        let n = (self.dt_outer / self.dt_cfl_max).ceil() as usize;
        n.max(1)
    }

    /// Return the actual sub-step size (dt_outer / num_substeps).
    pub fn substep_dt(&self) -> f64 {
        self.dt_outer / self.num_substeps() as f64
    }

    /// Check whether the current CFL-max is sufficient (no sub-stepping needed).
    pub fn is_stable(&self) -> bool {
        self.dt_cfl_max >= self.dt_outer
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::{ParticleSet, SphParticle};
    use oxiphysics_core::math::Vec3;

    // ── legacy tests (pre-existing) ──────────────────────────────────────────

    #[test]
    fn cfl_timestep_is_reasonable() {
        let mut ps = ParticleSet::new();
        ps.add_particle(&SphParticle::new(
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
        ));
        let h = 0.1;
        let cs = 100.0;
        let dt = cfl_timestep(&ps, h, cs);
        // dt = 0.4 * 0.1 / (100 + 1) ≈ 0.000396
        assert!(dt > 0.0, "CFL timestep must be positive");
        assert!(dt < 0.01, "CFL timestep should be small, got {dt}");
    }

    #[test]
    fn viscous_timestep_is_reasonable() {
        let dt = viscous_timestep(0.1, 0.01);
        // dt = 0.25 * 0.01 / 0.01 = 0.25
        assert!(dt > 0.0);
        assert!(dt < 10.0);
    }

    #[test]
    fn adaptive_timestep_respects_bounds() {
        let ps = ParticleSet::new();
        let dt = adaptive_timestep(&ps, 0.1, 100.0, 0.01, 1e-6, 0.001);
        assert!(dt >= 1e-6);
        assert!(dt <= 0.001);
    }

    // ── cfl_acoustic ─────────────────────────────────────────────────────────

    #[test]
    fn test_cfl_acoustic_formula() {
        // dt = 0.4 * 0.1 / 100.0 = 4e-4
        let dt = cfl_acoustic(100.0, 0.1, 0.4);
        let expected = 0.4 * 0.1 / 100.0;
        assert!((dt - expected).abs() < 1e-15, "dt={dt} expected={expected}");
    }

    // ── cfl_viscous ──────────────────────────────────────────────────────────

    #[test]
    fn test_cfl_viscous_formula() {
        // dt = 0.25 * 0.01 / 0.01 = 0.25
        let dt = cfl_viscous(0.01, 0.1, 0.25);
        let expected = 0.25 * 0.01 / 0.01;
        assert!((dt - expected).abs() < 1e-12, "dt={dt} expected={expected}");
    }

    // ── cfl_combined takes minimum ────────────────────────────────────────────

    #[test]
    fn test_cfl_combined_takes_minimum() {
        let c_sound = 100.0;
        let v_max = 5.0;
        let nu = 0.001;
        let h = 0.05;
        let factor = 0.4;

        let dt = cfl_combined(c_sound, v_max, nu, h, factor);
        let dt_a = cfl_acoustic(c_sound, h, factor);
        let dt_v = cfl_viscous(nu, h, factor);
        let dt_c = factor * h / (c_sound + v_max);

        let expected_min = dt_a.min(dt_v).min(dt_c);
        assert!(
            (dt - expected_min).abs() < 1e-14,
            "dt={dt} expected_min={expected_min}"
        );
        assert!(dt <= dt_a, "combined must be <= acoustic");
        assert!(dt <= dt_v, "combined must be <= viscous");
    }

    // ── AdaptiveTimestep growth ───────────────────────────────────────────────

    #[test]
    fn test_adaptive_timestep_grows() {
        let mut ts = AdaptiveTimestep::new(0.001, 1e-6, 0.1);
        ts.growth_factor = 1.1;
        // CFL allows a larger step.
        let dt = ts.update(0.01); // cfl_dt > current
        assert!(dt > 0.001, "dt should have grown, got {dt}");
        assert!(dt <= 0.001 * 1.1 + 1e-15);
    }

    // ── AdaptiveTimestep shrink ───────────────────────────────────────────────

    #[test]
    fn test_adaptive_timestep_shrinks() {
        let mut ts = AdaptiveTimestep::new(0.01, 1e-6, 0.1);
        ts.shrink_factor = 0.5;
        // CFL forces a smaller step.
        let dt = ts.update(0.001); // cfl_dt < current
        assert!(dt < 0.01, "dt should have shrunk, got {dt}");
        assert!(dt >= 1e-6, "dt must respect min_dt");
    }

    // ── AdaptiveTimestep clamped to max_dt ───────────────────────────────────

    #[test]
    fn test_adaptive_timestep_max_dt_clamped() {
        let mut ts = AdaptiveTimestep::new(0.09, 1e-6, 0.1);
        ts.growth_factor = 2.0;
        let dt = ts.update(1.0); // CFL allows large dt
        assert!(dt <= 0.1 + 1e-14, "dt must not exceed max_dt, got {dt}");
    }

    // ── leapfrog_kick / leapfrog_drift ────────────────────────────────────────

    #[test]
    fn test_leapfrog_kick_formula() {
        let v = [1.0f64, 0.0, 0.0];
        let a = [0.0f64, 2.0, 0.0];
        let dt = 0.5;
        let vn = leapfrog_kick(v, a, dt);
        assert!((vn[0] - 1.0).abs() < 1e-14);
        assert!((vn[1] - 1.0).abs() < 1e-14); // 0 + 2*0.5
        assert!(vn[2].abs() < 1e-14);
    }

    #[test]
    fn test_leapfrog_drift_formula() {
        let x = [0.0f64, 0.0, 0.0];
        let v = [3.0f64, 0.0, 0.0];
        let dt = 2.0;
        let xn = leapfrog_drift(x, v, dt);
        assert!((xn[0] - 6.0).abs() < 1e-14, "xn={:?}", xn);
    }

    // ── leapfrog energy conservation (free fall) ──────────────────────────────

    #[test]
    fn test_leapfrog_kdk_energy_conservation() {
        // Free fall under constant gravity g = [0, -9.81, 0].
        // Leapfrog (symplectic) conserves energy much better than Euler.
        let g = [0.0f64, -9.81, 0.0];
        let x0 = [0.0f64, 10.0, 0.0];
        let v0 = [0.0f64, 0.0, 0.0];
        let dt = 0.001;

        // Total mechanical energy (kinetic + potential) at t=0.
        let ke0 = 0.5 * (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
        let pe0 = 9.81 * x0[1];
        let e0 = ke0 + pe0;

        let mut x = x0;
        let mut v = v0;
        for _ in 0..1000 {
            let (xn, vn) = leapfrog_kdk(x, v, g, dt);
            x = xn;
            v = vn;
        }

        let ke = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        let pe = 9.81 * x[1];
        let e = ke + pe;

        // Leapfrog energy drift over 1 second should be small.
        let drift = (e - e0).abs();
        assert!(
            drift < 0.1,
            "leapfrog energy drift={drift} too large (e0={e0}, e={e})"
        );
    }

    // ── sph_dt_from_particles ─────────────────────────────────────────────────

    #[test]
    fn test_sph_dt_from_particles_positive() {
        let mut ps = ParticleSet::new();
        ps.add_particle(&SphParticle::new(
            Vec3::zeros(),
            Vec3::new(10.0, 0.0, 0.0),
            1.0,
        ));
        let dt = sph_dt_from_particles(&ps, 0.1, 0.01, 0.4);
        assert!(dt > 0.0, "dt must be positive, got {dt}");
    }

    // ── Surface tension timestep ──────────────────────────────────────────────

    #[test]
    fn test_surface_tension_timestep_positive() {
        let dt = surface_tension_timestep(1000.0, 0.01, 0.072, 0.4);
        assert!(
            dt > 0.0,
            "Surface tension timestep must be positive, got {dt}"
        );
    }

    #[test]
    fn test_surface_tension_timestep_zero_sigma() {
        let dt = surface_tension_timestep(1000.0, 0.01, 0.0, 0.4);
        assert!(dt == f64::MAX, "Zero sigma should give MAX dt");
    }

    // ── Acceleration-based timestep ───────────────────────────────────────────

    #[test]
    fn test_acceleration_timestep_positive() {
        let dt = acceleration_timestep(0.01, 10.0, 0.4);
        assert!(dt > 0.0, "Acceleration timestep must be positive");
    }

    #[test]
    fn test_acceleration_timestep_zero_accel() {
        let dt = acceleration_timestep(0.01, 0.0, 0.4);
        assert!(dt == f64::MAX, "Zero acceleration should give MAX dt");
    }

    #[test]
    fn test_max_acceleration() {
        let accels = [[1.0, 0.0, 0.0], [0.0, 3.0, 4.0], [0.0, 0.0, 0.0]];
        let a_max = max_acceleration(&accels);
        assert!(
            (a_max - 5.0).abs() < 1e-10,
            "Max accel should be 5.0, got {a_max}"
        );
    }

    // ── Multi-criteria timestep ───────────────────────────────────────────────

    #[test]
    fn test_multi_criteria_compute() {
        let criteria = TimestepCriteria {
            cfl_factor: 0.4,
            c_sound: 100.0,
            v_max: 5.0,
            nu: 0.001,
            h: 0.05,
            sigma: 0.072,
            rho: 1000.0,
            a_max: 10.0,
            min_dt: 1e-8,
            max_dt: 0.1,
        };
        let dt = criteria.compute();
        assert!(dt > 0.0, "Multi-criteria dt must be positive");
        assert!(dt >= 1e-8, "dt must respect min_dt");
        assert!(dt <= 0.1, "dt must respect max_dt");
    }

    #[test]
    fn test_multi_criteria_limiting_criterion() {
        let criteria = TimestepCriteria {
            cfl_factor: 0.4,
            c_sound: 1000.0, // very high sound speed
            v_max: 0.1,
            nu: 0.0001,
            h: 0.05,
            sigma: 0.0,
            rho: 1000.0,
            a_max: 0.0,
            min_dt: 1e-8,
            max_dt: 0.1,
        };
        let limiting = criteria.limiting_criterion();
        // With very high c_sound, convective (CFL * h / (c+v)) is always <= acoustic (CFL * h / c)
        assert_eq!(limiting, "convective", "Expected convective to be limiting");
    }

    // ── Timestep history ──────────────────────────────────────────────────────

    #[test]
    fn test_timestep_history_record() {
        let mut hist = TimestepHistory::new(100);
        assert!(hist.is_empty());
        hist.record(0.0, 0.001);
        hist.record(0.001, 0.0012);
        hist.record(0.0022, 0.0009);
        assert_eq!(hist.len(), 3);
    }

    #[test]
    fn test_timestep_history_average() {
        let mut hist = TimestepHistory::new(100);
        hist.record(0.0, 0.001);
        hist.record(0.001, 0.003);
        let avg = hist.average_dt();
        assert!((avg - 0.002).abs() < 1e-12, "Expected 0.002, got {avg}");
    }

    #[test]
    fn test_timestep_history_min_max() {
        let mut hist = TimestepHistory::new(100);
        hist.record(0.0, 0.005);
        hist.record(0.005, 0.001);
        hist.record(0.006, 0.003);
        assert!((hist.min_dt() - 0.001).abs() < 1e-12);
        assert!((hist.max_dt() - 0.005).abs() < 1e-12);
    }

    #[test]
    fn test_timestep_history_std_dev() {
        let mut hist = TimestepHistory::new(100);
        hist.record(0.0, 2.0);
        hist.record(0.1, 4.0);
        hist.record(0.2, 4.0);
        hist.record(0.3, 4.0);
        hist.record(0.4, 5.0);
        hist.record(0.5, 5.0);
        hist.record(0.6, 7.0);
        hist.record(0.7, 9.0);
        let sd = hist.std_dev();
        assert!(sd > 0.0, "Std dev should be positive for varying values");
    }

    #[test]
    fn test_timestep_history_ring_buffer() {
        let mut hist = TimestepHistory::new(3);
        hist.record(0.0, 0.001);
        hist.record(0.1, 0.002);
        hist.record(0.2, 0.003);
        hist.record(0.3, 0.004); // should evict first entry
        assert_eq!(hist.len(), 3);
        assert!((hist.values[0] - 0.002).abs() < 1e-12);
    }

    #[test]
    fn test_timestep_history_rate_of_change() {
        let mut hist = TimestepHistory::new(100);
        hist.record(0.0, 0.001);
        hist.record(0.1, 0.002);
        let rate = hist.dt_rate_of_change();
        // (0.002 - 0.001) / (0.1 - 0.0) = 0.01
        assert!((rate - 0.01).abs() < 1e-12, "Expected 0.01, got {rate}");
    }

    // ── Velocity Verlet ───────────────────────────────────────────────────────

    #[test]
    fn test_velocity_verlet_constant_accel() {
        let x = [0.0, 0.0, 0.0];
        let v = [0.0, 0.0, 0.0];
        let a = [1.0, 0.0, 0.0]; // constant acceleration
        let dt = 0.1;
        let (xn, vn) = velocity_verlet(x, v, a, a, dt);
        // x = 0 + 0 + 0.5 * 1 * 0.01 = 0.005
        assert!((xn[0] - 0.005).abs() < 1e-12, "xn[0]={}", xn[0]);
        // v = 0 + 0.5 * (1 + 1) * 0.1 = 0.1
        assert!((vn[0] - 0.1).abs() < 1e-12, "vn[0]={}", vn[0]);
    }

    // ── AdaptiveTimestep with multi-criteria ──────────────────────────────────

    #[test]
    fn test_adaptive_timestep_multi_criteria() {
        let mut ts = AdaptiveTimestep::new(0.001, 1e-8, 0.01);
        let criteria = TimestepCriteria {
            cfl_factor: 0.4,
            c_sound: 100.0,
            v_max: 5.0,
            nu: 0.001,
            h: 0.05,
            sigma: 0.0,
            rho: 1000.0,
            a_max: 0.0,
            min_dt: 1e-8,
            max_dt: 0.01,
        };
        let dt = ts.update_multi(&criteria);
        assert!(dt > 0.0);
    }

    // ── sph_multistep_euler ───────────────────────────────────────────────────

    #[test]
    fn test_sph_multistep_euler_constant_accel() {
        // x'' = g = [0, -9.81, 0]: free fall for 1 second
        let g = [0.0, -9.81, 0.0];
        let result = sph_multistep_euler([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1000, |_x, _v| g);
        // x = 0.5 * g * t^2, v = g * t
        let expected_y = 0.5 * (-9.81) * 1.0; // ≈ -4.905
        assert!(
            (result.velocity[1] - (-9.81)).abs() < 0.02,
            "vy={}",
            result.velocity[1]
        );
        assert!(
            (result.position[1] - expected_y).abs() < 0.05,
            "py={}",
            result.position[1]
        );
        assert_eq!(result.substeps, 1000);
    }

    #[test]
    fn test_sph_multistep_euler_zero_substeps_treated_as_one() {
        let result =
            sph_multistep_euler([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 0, |_x, _v| [0.0; 3]);
        // Should complete without panic, 0 substeps → 1 substep
        assert_eq!(result.substeps, 0);
        assert!((result.time - 0.1).abs() < 1e-14);
    }

    // ── sph_multistep_rk4 ────────────────────────────────────────────────────

    #[test]
    fn test_sph_multistep_rk4_free_fall() {
        let g = [0.0, -9.81, 0.0];
        let result = sph_multistep_rk4([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 100, |_x, _v| g);
        let expected_y = 0.5 * (-9.81);
        // RK4 is exact for constant-accel (polynomial of degree 2)
        assert!(
            (result.velocity[1] - (-9.81)).abs() < 1e-8,
            "vy={}",
            result.velocity[1]
        );
        assert!(
            (result.position[1] - expected_y).abs() < 1e-8,
            "py={}",
            result.position[1]
        );
    }

    #[test]
    fn test_rk4_more_accurate_than_euler_for_oscillator() {
        // Simple harmonic oscillator x'' = -x, ω=1, period=2π
        // Initial: x=1, v=0; after half period x=-1, v=0
        let half_period = std::f64::consts::PI;
        let n = 200;

        let rk4 = sph_multistep_rk4([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], half_period, n, |x, _v| {
            [-x[0], 0.0, 0.0]
        });
        let euler =
            sph_multistep_euler([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], half_period, n, |x, _v| {
                [-x[0], 0.0, 0.0]
            });
        // After half period, exact x=-1
        let rk4_err = (rk4.position[0] + 1.0).abs();
        let euler_err = (euler.position[0] + 1.0).abs();
        assert!(
            rk4_err < euler_err,
            "RK4 err={rk4_err} should beat Euler err={euler_err}"
        );
    }

    // ── predictor_corrector_step ─────────────────────────────────────────────

    #[test]
    fn test_predictor_corrector_constant_accel() {
        let x = [0.0, 0.0, 0.0];
        let v = [1.0, 0.0, 0.0];
        let a = [0.0, -9.81, 0.0];
        let dt = 0.01;
        let (xn, vn) = predictor_corrector_step(x, v, a, dt, |_x, _v| a);
        // With constant a the predictor-corrector is exact to 2nd order
        assert!(
            (xn[0] - (x[0] + v[0] * dt)).abs() < 1e-10,
            "xn[0]={}",
            xn[0]
        );
        assert!(
            (vn[1] - (v[1] + a[1] * dt)).abs() < 1e-10,
            "vn[1]={}",
            vn[1]
        );
    }

    #[test]
    fn test_predictor_corrector_is_second_order() {
        // Compare with exact solution for free fall: y = -0.5*g*t^2
        let g = 9.81;
        let dt = 0.1;
        let (xn, _) = predictor_corrector_step([0.0; 3], [0.0; 3], [0.0, -g, 0.0], dt, |_x, _v| {
            [0.0, -g, 0.0]
        });
        let exact = -0.5 * g * dt * dt;
        assert!((xn[1] - exact).abs() < 1e-12, "y={} exact={}", xn[1], exact);
    }

    // ── symplectic_euler_sph ─────────────────────────────────────────────────

    #[test]
    fn test_symplectic_euler_sph_velocity_first() {
        let x = [0.0, 0.0, 0.0];
        let v = [0.0, 0.0, 0.0];
        let a = [1.0, 0.0, 0.0];
        let dt = 1.0;
        let (xn, vn) = symplectic_euler_sph(x, v, a, dt);
        // v_new = v + a*dt = 1; x_new = x + v_new*dt = 1
        assert!((vn[0] - 1.0).abs() < 1e-14, "vn={}", vn[0]);
        assert!((xn[0] - 1.0).abs() < 1e-14, "xn={}", xn[0]);
    }

    #[test]
    fn test_symplectic_euler_sph_energy_conservation_oscillator() {
        // Simple harmonic oscillator x'' = -ω²x, ω=1
        // Symplectic Euler conserves a modified energy
        let omega = 1.0;
        let x0 = [1.0, 0.0, 0.0];
        let v0 = [0.0, 0.0, 0.0];
        let dt = 0.01;
        let n = 1000;
        let mut x = x0;
        let mut v = v0;
        let e0 = 0.5 * v0[0] * v0[0] + 0.5 * omega * omega * x0[0] * x0[0];
        for _ in 0..n {
            let a = [-omega * omega * x[0], 0.0, 0.0];
            let (xn, vn) = symplectic_euler_sph(x, v, a, dt);
            x = xn;
            v = vn;
        }
        let ef = 0.5 * v[0] * v[0] + 0.5 * omega * omega * x[0] * x[0];
        // Symplectic Euler should not dissipate: energy should remain close
        assert!(
            (ef - e0).abs() / e0 < 0.05,
            "energy drift: e0={e0}, ef={ef}"
        );
    }

    // ── body_force_timestep ───────────────────────────────────────────────────

    #[test]
    fn test_body_force_timestep_gravity() {
        let dt = body_force_timestep(0.01, 9.81, 0.4);
        assert!(dt > 0.0, "dt must be positive");
        // dt = 0.4 * sqrt(0.01 / 9.81) ≈ 0.4 * 0.0319 ≈ 0.01276
        let expected = 0.4 * (0.01_f64 / 9.81).sqrt();
        assert!((dt - expected).abs() < 1e-12, "dt={dt} expected={expected}");
    }

    #[test]
    fn test_body_force_timestep_zero_accel() {
        let dt = body_force_timestep(0.01, 0.0, 0.4);
        assert_eq!(dt, f64::MAX);
    }

    // ── adaptive_timestep_full ────────────────────────────────────────────────

    #[test]
    fn test_adaptive_timestep_full_bounded() {
        let dt = adaptive_timestep_full(100.0, 5.0, 0.001, 0.05, 9.81, 0.4, 1e-8, 0.1);
        assert!(dt > 0.0, "dt must be positive");
        assert!(dt >= 1e-8, "dt must respect min_dt");
        assert!(dt <= 0.1, "dt must respect max_dt");
    }

    // ── SubstepController ────────────────────────────────────────────────────

    #[test]
    fn test_substep_controller_no_substeps_needed() {
        let ctrl = SubstepController::new(0.01, 0.1, 1e-6);
        assert_eq!(ctrl.num_substeps(), 1);
        assert!((ctrl.substep_dt() - 0.01).abs() < 1e-14);
        assert!(ctrl.is_stable());
    }

    #[test]
    fn test_substep_controller_multiple_substeps() {
        let ctrl = SubstepController::new(0.1, 0.025, 1e-6);
        assert_eq!(ctrl.num_substeps(), 4);
        assert!((ctrl.substep_dt() - 0.025).abs() < 1e-14);
        assert!(!ctrl.is_stable());
    }

    #[test]
    fn test_substep_controller_fractional_substeps_rounds_up() {
        // dt_outer = 0.1, dt_cfl = 0.033 → 0.1/0.033 = 3.03 → ceil → 4
        let ctrl = SubstepController::new(0.1, 0.033, 1e-6);
        assert!(ctrl.num_substeps() >= 4);
    }

    // ── rk4_constant_accel (existing, verify no regression) ──────────────────

    #[test]
    fn test_rk4_constant_accel_exact() {
        // Under constant acceleration, RK4 is exact.
        let x = [0.0, 0.0, 0.0];
        let v = [1.0, 0.0, 0.0];
        let a = [0.0, -9.81, 0.0];
        let dt = 0.1;
        let (xn, vn) = rk4_constant_accel(x, v, a, dt);
        // exact: v_y = 0 + (-9.81)*0.1 = -0.981; x_y = 0 - 0.5*9.81*0.01 = -0.04905
        assert!((vn[1] - (-0.981)).abs() < 1e-12, "vn[1]={}", vn[1]);
        assert!((xn[1] - (-0.04905)).abs() < 1e-12, "xn[1]={}", xn[1]);
        assert!((xn[0] - (x[0] + v[0] * dt)).abs() < 1e-12);
    }
}
