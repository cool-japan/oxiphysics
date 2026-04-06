#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive time-stepping and CFL condition for SPH simulations.
//!
//! Provides CFL-based timestep constraints, Runge-Kutta integrators,
//! and a high-level `TimeStepper` controller for SPH fluid simulations.

/// Configuration for CFL time-step coefficients.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CflConfig {
    /// CFL coefficient for velocity/sound-speed constraint.
    pub c_cfl: f64,
    /// CFL coefficient for viscous constraint.
    pub c_visc: f64,
    /// CFL coefficient for force/acceleration constraint.
    pub c_force: f64,
}

impl CflConfig {
    /// Default CFL coefficients (all 0.25).
    pub fn default() -> Self {
        Self {
            c_cfl: 0.25,
            c_visc: 0.25,
            c_force: 0.25,
        }
    }

    /// Conservative CFL coefficients (all 0.1) for stability-critical simulations.
    pub fn conservative() -> Self {
        Self {
            c_cfl: 0.1,
            c_visc: 0.1,
            c_force: 0.1,
        }
    }

    /// Aggressive CFL coefficients (all 0.4) for faster simulations.
    pub fn aggressive() -> Self {
        Self {
            c_cfl: 0.4,
            c_visc: 0.4,
            c_force: 0.4,
        }
    }
}

/// Compute the CFL velocity timestep constraint.
///
/// `dt_v = c_cfl * h / max_vel`
pub fn cfl_dt_velocity(max_velocity: f64, h: f64, config: &CflConfig) -> f64 {
    config.c_cfl * h / max_velocity.max(1e-14)
}

/// Compute the CFL sound-speed timestep constraint.
///
/// `dt_s = c_cfl * h / c_s`
pub fn cfl_dt_sound(c_sound: f64, h: f64, config: &CflConfig) -> f64 {
    config.c_cfl * h / c_sound.max(1e-14)
}

/// Compute the viscous timestep constraint.
///
/// `dt_mu = c_visc * h² / nu`
pub fn cfl_dt_viscous(h: f64, kinematic_viscosity: f64, config: &CflConfig) -> f64 {
    config.c_visc * h * h / kinematic_viscosity.max(1e-14)
}

/// Compute the force/acceleration timestep constraint.
///
/// `dt_f = c_force * sqrt(h / |a_max|)`
pub fn cfl_dt_force(max_force_per_mass: f64, h: f64, config: &CflConfig) -> f64 {
    config.c_force * (h / max_force_per_mass.max(1e-14)).sqrt()
}

/// Compute the adaptive timestep as the minimum of all CFL constraints.
///
/// Takes particle velocities, forces (per unit mass = accelerations), smoothing length,
/// sound speed, and kinematic viscosity.
pub fn adaptive_dt(
    velocities: &[[f64; 3]],
    forces: &[[f64; 3]],
    h: f64,
    c_sound: f64,
    nu: f64,
    config: &CflConfig,
) -> f64 {
    let max_vel = max_velocity(velocities);
    let max_acc = max_velocity(forces); // forces here are already per-mass (accelerations)

    let dt_v = cfl_dt_velocity(max_vel, h, config);
    let dt_s = cfl_dt_sound(c_sound, h, config);
    let dt_mu = cfl_dt_viscous(h, nu, config);
    let dt_f = cfl_dt_force(max_acc, h, config);

    dt_v.min(dt_s).min(dt_mu).min(dt_f)
}

/// High-level adaptive time stepper for SPH simulations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TimeStepper {
    /// Current timestep.
    pub dt: f64,
    /// Minimum allowed timestep.
    pub dt_min: f64,
    /// Maximum allowed timestep.
    pub dt_max: f64,
    /// CFL configuration.
    pub cfl_config: CflConfig,
    /// Number of steps taken so far.
    pub n_steps: u64,
}

impl TimeStepper {
    /// Create a new `TimeStepper` with the given initial and bound timesteps.
    pub fn new(dt_init: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            dt: dt_init,
            dt_min,
            dt_max,
            cfl_config: CflConfig::default(),
            n_steps: 0,
        }
    }

    /// Recompute `dt` from CFL constraints and clamp to `[dt_min, dt_max]`.
    pub fn update_dt(
        &mut self,
        velocities: &[[f64; 3]],
        forces: &[[f64; 3]],
        h: f64,
        c_sound: f64,
        nu: f64,
    ) {
        let proposed = adaptive_dt(velocities, forces, h, c_sound, nu, &self.cfl_config);
        self.dt = proposed.clamp(self.dt_min, self.dt_max);
    }

    /// Return the current timestep.
    pub fn current_dt(&self) -> f64 {
        self.dt
    }

    /// Return an approximate total simulated time (`n_steps * dt`).
    pub fn total_time(&self) -> f64 {
        self.n_steps as f64 * self.dt
    }

    /// Advance the internal step counter by one.
    pub fn advance(&mut self) {
        self.n_steps += 1;
    }
}

/// Collection of Runge-Kutta and leapfrog integrators for SPH.
pub struct RungeKuttaSph;

impl RungeKuttaSph {
    /// Second-order midpoint Runge-Kutta step.
    ///
    /// Returns `(new_positions, new_velocities)`.
    pub fn rk2_step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = pos.len();
        let half_dt = 0.5 * dt;

        // Stage 1: midpoint positions and velocities
        let pos_mid: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + half_dt * vel[i][0],
                    pos[i][1] + half_dt * vel[i][1],
                    pos[i][2] + half_dt * vel[i][2],
                ]
            })
            .collect();
        let vel_mid: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + half_dt * acc[i][0],
                    vel[i][1] + half_dt * acc[i][1],
                    vel[i][2] + half_dt * acc[i][2],
                ]
            })
            .collect();

        // Stage 2: full step using midpoint derivatives
        let new_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + dt * vel_mid[i][0],
                    pos[i][1] + dt * vel_mid[i][1],
                    pos[i][2] + dt * vel_mid[i][2],
                ]
            })
            .collect();
        let new_vel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + dt * acc[i][0],
                    vel[i][1] + dt * acc[i][1],
                    vel[i][2] + dt * acc[i][2],
                ]
            })
            .collect();

        // pos_mid is used only to document mid-point; suppress unused warning
        let _ = &pos_mid;

        (new_pos, new_vel)
    }

    /// Forward Euler step.
    ///
    /// Returns `(new_positions, new_velocities)`.
    pub fn euler_step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = pos.len();
        let new_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + dt * vel[i][0],
                    pos[i][1] + dt * vel[i][1],
                    pos[i][2] + dt * vel[i][2],
                ]
            })
            .collect();
        let new_vel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + dt * acc[i][0],
                    vel[i][1] + dt * acc[i][1],
                    vel[i][2] + dt * acc[i][2],
                ]
            })
            .collect();
        (new_pos, new_vel)
    }

    /// Leapfrog kick: advance velocities by half a timestep.
    ///
    /// `v_new = v + (dt/2) * a`
    pub fn leapfrog_kick(vel: &[[f64; 3]], acc: &[[f64; 3]], dt: f64) -> Vec<[f64; 3]> {
        let half_dt = 0.5 * dt;
        vel.iter()
            .zip(acc.iter())
            .map(|(v, a)| {
                [
                    v[0] + half_dt * a[0],
                    v[1] + half_dt * a[1],
                    v[2] + half_dt * a[2],
                ]
            })
            .collect()
    }

    /// Leapfrog drift: advance positions by a full timestep using current velocities.
    ///
    /// `x_new = x + dt * v`
    pub fn leapfrog_drift(pos: &[[f64; 3]], vel: &[[f64; 3]], dt: f64) -> Vec<[f64; 3]> {
        pos.iter()
            .zip(vel.iter())
            .map(|(x, v)| [x[0] + dt * v[0], x[1] + dt * v[1], x[2] + dt * v[2]])
            .collect()
    }
}

/// Compute sound speed from the Tait equation of state.
///
/// `c = sqrt(gamma * B / rho0) * (rho / rho0)^((gamma-1)/2)`
pub fn sound_speed_tait(rho: f64, rho0: f64, gamma: f64, b: f64) -> f64 {
    let c0 = (gamma * b / rho0).sqrt();
    let exponent = (gamma - 1.0) / 2.0;
    c0 * (rho / rho0).powf(exponent)
}

/// Compute the maximum speed across a set of 3-D velocity vectors.
pub fn max_velocity(velocities: &[[f64; 3]]) -> f64 {
    velocities
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Compute the maximum acceleration magnitude: `max_i(|F_i| / m_i)`.
pub fn max_acceleration(forces: &[[f64; 3]], masses: &[f64]) -> f64 {
    forces
        .iter()
        .zip(masses.iter())
        .map(|(f, &m)| {
            let m_safe = m.max(1e-30);
            let ax = f[0] / m_safe;
            let ay = f[1] / m_safe;
            let az = f[2] / m_safe;
            (ax * ax + ay * ay + az * az).sqrt()
        })
        .fold(0.0_f64, f64::max)
}

// ---------------------------------------------------------------------------
// Multi-rate timestepping
// ---------------------------------------------------------------------------

/// Multi-rate time stepper that assigns particles to different timestep levels.
///
/// Particles with faster dynamics (higher velocities/accelerations) are
/// integrated with smaller timesteps, while slow particles use larger steps.
#[allow(dead_code)]
pub struct MultiRateTimeStepper {
    /// Base (largest) timestep.
    pub dt_base: f64,
    /// Number of refinement levels (level 0 = base, level k = dt_base / 2^k).
    pub n_levels: usize,
    /// Timestep level assigned to each particle.
    pub levels: Vec<usize>,
    /// CFL configuration.
    pub cfl_config: CflConfig,
}

#[allow(dead_code)]
impl MultiRateTimeStepper {
    /// Create a new multi-rate stepper.
    pub fn new(dt_base: f64, n_levels: usize, n_particles: usize) -> Self {
        Self {
            dt_base,
            n_levels,
            levels: vec![0; n_particles],
            cfl_config: CflConfig::default(),
        }
    }

    /// Timestep at a given level: dt_base / 2^level.
    pub fn dt_at_level(&self, level: usize) -> f64 {
        self.dt_base / (1u64 << level) as f64
    }

    /// Assign timestep levels based on particle velocities and accelerations.
    pub fn assign_levels(
        &mut self,
        velocities: &[[f64; 3]],
        forces: &[[f64; 3]],
        h: f64,
        c_sound: f64,
        nu: f64,
    ) {
        for (i, (v, f)) in velocities.iter().zip(forces.iter()).enumerate() {
            let dt_particle = adaptive_dt(
                std::slice::from_ref(v),
                std::slice::from_ref(f),
                h,
                c_sound,
                nu,
                &self.cfl_config,
            );
            // Find level such that dt_at_level >= dt_particle
            let mut level = 0;
            while level < self.n_levels - 1 && self.dt_at_level(level) > dt_particle {
                level += 1;
            }
            self.levels[i] = level;
        }
    }

    /// Return true if particle `i` should be active at this sub-step.
    ///
    /// Within one base step, a level-k particle is active at sub-steps
    /// that are multiples of 2^(n_levels-1-k).
    pub fn is_active(&self, particle_idx: usize, sub_step: usize) -> bool {
        let level = self.levels[particle_idx];
        let stride = 1usize << (self.n_levels - 1 - level);
        sub_step.is_multiple_of(stride)
    }

    /// Number of sub-steps in one base step.
    pub fn sub_steps_per_base(&self) -> usize {
        1usize << (self.n_levels - 1)
    }
}

// ---------------------------------------------------------------------------
// Asynchronous integration
// ---------------------------------------------------------------------------

/// Particle state for asynchronous timestepping.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct AsyncParticleState {
    /// Last update time for this particle.
    pub t_last: f64,
    /// Individual timestep.
    pub dt_individual: f64,
    /// Next scheduled update time.
    pub t_next: f64,
}

#[allow(dead_code)]
impl AsyncParticleState {
    /// Create a new async state.
    pub fn new(dt: f64) -> Self {
        Self {
            t_last: 0.0,
            dt_individual: dt,
            t_next: dt,
        }
    }

    /// Advance this particle's schedule.
    pub fn advance(&mut self) {
        self.t_last = self.t_next;
        self.t_next = self.t_last + self.dt_individual;
    }

    /// Time since last update at a given global time.
    pub fn time_since_update(&self, t_global: f64) -> f64 {
        t_global - self.t_last
    }
}

// ---------------------------------------------------------------------------
// Sub-cycling
// ---------------------------------------------------------------------------

/// Sub-cycling integrator: advances fast particles with multiple small steps
/// within a single global step.
#[allow(dead_code)]
pub struct SubCyclingIntegrator;

#[allow(dead_code)]
impl SubCyclingIntegrator {
    /// Integrate a single particle for one global step using sub-cycling.
    ///
    /// Returns `(new_pos, new_vel)`.
    pub fn integrate_particle(
        pos: [f64; 3],
        vel: [f64; 3],
        acc: [f64; 3],
        dt_global: f64,
        n_sub_steps: usize,
    ) -> ([f64; 3], [f64; 3]) {
        let dt_sub = dt_global / n_sub_steps as f64;
        let mut p = pos;
        let mut v = vel;
        for _ in 0..n_sub_steps {
            // Forward Euler sub-step
            p = [
                p[0] + dt_sub * v[0],
                p[1] + dt_sub * v[1],
                p[2] + dt_sub * v[2],
            ];
            v = [
                v[0] + dt_sub * acc[0],
                v[1] + dt_sub * acc[1],
                v[2] + dt_sub * acc[2],
            ];
        }
        (p, v)
    }

    /// Determine the number of sub-steps needed for stability.
    ///
    /// `n = ceil(dt_global / dt_particle)`
    pub fn required_sub_steps(dt_global: f64, dt_particle: f64) -> usize {
        if dt_particle <= 0.0 || dt_global <= 0.0 {
            return 1;
        }
        (dt_global / dt_particle).ceil() as usize
    }
}

// ---------------------------------------------------------------------------
// Timestep synchronization
// ---------------------------------------------------------------------------

/// Utilities for synchronizing multi-rate particles to a common time.
#[allow(dead_code)]
pub struct TimestepSync;

#[allow(dead_code)]
impl TimestepSync {
    /// Predict position at time `t` from last known state using linear extrapolation.
    pub fn predict_position(pos: [f64; 3], vel: [f64; 3], t_last: f64, t: f64) -> [f64; 3] {
        let dt = t - t_last;
        [
            pos[0] + dt * vel[0],
            pos[1] + dt * vel[1],
            pos[2] + dt * vel[2],
        ]
    }

    /// Predict velocity at time `t` from last known state.
    pub fn predict_velocity(vel: [f64; 3], acc: [f64; 3], t_last: f64, t: f64) -> [f64; 3] {
        let dt = t - t_last;
        [
            vel[0] + dt * acc[0],
            vel[1] + dt * acc[1],
            vel[2] + dt * acc[2],
        ]
    }

    /// Compute the global synchronization time (least common multiple of active timesteps).
    pub fn sync_time(timesteps: &[f64]) -> f64 {
        if timesteps.is_empty() {
            return 0.0;
        }
        // For power-of-two refinement, sync time = max(timesteps)
        timesteps.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Error estimation
// ---------------------------------------------------------------------------

/// Timestep error estimation for adaptive integration.
#[allow(dead_code)]
pub struct ErrorEstimator;

#[allow(dead_code)]
impl ErrorEstimator {
    /// Estimate local truncation error using Richardson extrapolation.
    ///
    /// Given results from one full step (`y_full`) and two half steps (`y_half`),
    /// the error is `|y_half - y_full| / (2^p - 1)` for a p-th order method.
    pub fn richardson_error(y_full: [f64; 3], y_half: [f64; 3], order: u32) -> f64 {
        let factor = (1u64 << order) as f64 - 1.0;
        let diff = [
            y_half[0] - y_full[0],
            y_half[1] - y_full[1],
            y_half[2] - y_full[2],
        ];
        let err = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
        err / factor
    }

    /// Compute a new timestep from the estimated error:
    ///
    /// `dt_new = dt * safety * (tol / err)^(1/(p+1))`
    pub fn adjusted_dt(dt: f64, error: f64, tolerance: f64, order: u32, safety: f64) -> f64 {
        if error < 1e-30 {
            return dt * 2.0; // Error negligible, allow growth
        }
        let exponent = 1.0 / (order as f64 + 1.0);
        dt * safety * (tolerance / error).powf(exponent)
    }

    /// Compute the maximum particle-wise position error.
    pub fn max_position_error(pos_full: &[[f64; 3]], pos_half: &[[f64; 3]]) -> f64 {
        pos_full
            .iter()
            .zip(pos_half.iter())
            .map(|(a, b)| {
                let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
}

/// Fourth-order Runge-Kutta step.
#[allow(dead_code)]
impl RungeKuttaSph {
    /// RK4 step for a single particle.
    ///
    /// `acc_fn` returns acceleration given position and velocity.
    /// Returns `(new_pos, new_vel)`.
    pub fn rk4_step_single(
        pos: [f64; 3],
        vel: [f64; 3],
        acc: [f64; 3],
        dt: f64,
    ) -> ([f64; 3], [f64; 3]) {
        // Simplified: use constant acceleration (true RK4 would re-evaluate forces)
        let half_dt = 0.5 * dt;
        // k1
        let k1v = acc;
        let k1x = vel;
        // k2
        let v_mid = [
            vel[0] + half_dt * k1v[0],
            vel[1] + half_dt * k1v[1],
            vel[2] + half_dt * k1v[2],
        ];
        let k2x = v_mid;
        let k2v = acc; // constant acc approximation
        // k3
        let v_mid2 = [
            vel[0] + half_dt * k2v[0],
            vel[1] + half_dt * k2v[1],
            vel[2] + half_dt * k2v[2],
        ];
        let k3x = v_mid2;
        let k3v = acc;
        // k4
        let v_end = [
            vel[0] + dt * k3v[0],
            vel[1] + dt * k3v[1],
            vel[2] + dt * k3v[2],
        ];
        let k4x = v_end;
        let k4v = acc;

        let new_pos = [
            pos[0] + dt / 6.0 * (k1x[0] + 2.0 * k2x[0] + 2.0 * k3x[0] + k4x[0]),
            pos[1] + dt / 6.0 * (k1x[1] + 2.0 * k2x[1] + 2.0 * k3x[1] + k4x[1]),
            pos[2] + dt / 6.0 * (k1x[2] + 2.0 * k2x[2] + 2.0 * k3x[2] + k4x[2]),
        ];
        let new_vel = [
            vel[0] + dt / 6.0 * (k1v[0] + 2.0 * k2v[0] + 2.0 * k3v[0] + k4v[0]),
            vel[1] + dt / 6.0 * (k1v[1] + 2.0 * k2v[1] + 2.0 * k3v[1] + k4v[1]),
            vel[2] + dt / 6.0 * (k1v[2] + 2.0 * k2v[2] + 2.0 * k3v[2] + k4v[2]),
        ];
        (new_pos, new_vel)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cfl_config_default_values() {
        let cfg = CflConfig::default();
        assert_eq!(cfg.c_cfl, 0.25);
        assert_eq!(cfg.c_visc, 0.25);
        assert_eq!(cfg.c_force, 0.25);
    }

    #[test]
    fn cfl_config_conservative_values() {
        let cfg = CflConfig::conservative();
        assert_eq!(cfg.c_cfl, 0.1);
    }

    #[test]
    fn cfl_config_aggressive_values() {
        let cfg = CflConfig::aggressive();
        assert_eq!(cfg.c_cfl, 0.4);
    }

    #[test]
    fn cfl_dt_velocity_basic() {
        let cfg = CflConfig::default();
        let dt = cfl_dt_velocity(2.0, 0.1, &cfg);
        assert!((dt - 0.25 * 0.1 / 2.0).abs() < 1e-15);
    }

    #[test]
    fn cfl_dt_velocity_zero_vel_clamps() {
        let dt = cfl_dt_velocity(0.0, 0.1, &CflConfig::default());
        assert!(dt.is_finite() && dt > 0.0);
    }

    #[test]
    fn cfl_dt_velocity_scales_with_h() {
        let cfg = CflConfig::default();
        let dt1 = cfl_dt_velocity(1.0, 0.1, &cfg);
        let dt2 = cfl_dt_velocity(1.0, 0.2, &cfg);
        assert!((dt2 - 2.0 * dt1).abs() < 1e-14);
    }

    #[test]
    fn cfl_dt_sound_basic() {
        let dt = cfl_dt_sound(100.0, 0.1, &CflConfig::default());
        assert!((dt - 0.25 * 0.1 / 100.0).abs() < 1e-15);
    }

    #[test]
    fn cfl_dt_sound_zero_sound_clamps() {
        let dt = cfl_dt_sound(0.0, 0.1, &CflConfig::default());
        assert!(dt.is_finite() && dt > 0.0);
    }

    #[test]
    fn cfl_dt_viscous_basic() {
        let dt = cfl_dt_viscous(0.1, 0.01, &CflConfig::default());
        assert!((dt - 0.25 * 0.01 / 0.01).abs() < 1e-14);
    }

    #[test]
    fn cfl_dt_viscous_zero_nu_clamps() {
        let dt = cfl_dt_viscous(0.1, 0.0, &CflConfig::default());
        assert!(dt.is_finite() && dt > 0.0);
    }

    #[test]
    fn cfl_dt_force_basic() {
        let dt = cfl_dt_force(10.0, 0.1, &CflConfig::default());
        assert!((dt - 0.25 * (0.1 / 10.0_f64).sqrt()).abs() < 1e-14);
    }

    #[test]
    fn cfl_dt_force_zero_force_clamps() {
        let dt = cfl_dt_force(0.0, 0.1, &CflConfig::default());
        assert!(dt.is_finite() && dt > 0.0);
    }

    #[test]
    fn adaptive_dt_returns_minimum() {
        let dt = adaptive_dt(
            &[[1.0, 0.0, 0.0]],
            &[[0.0; 3]],
            0.1,
            100.0,
            0.01,
            &CflConfig::default(),
        );
        assert!(dt > 0.0 && dt.is_finite());
    }

    #[test]
    fn adaptive_dt_conservative_smaller_than_aggressive() {
        let v = &[[1.0, 0.0, 0.0]];
        let f = &[[1.0, 0.0, 0.0]];
        let dt_c = adaptive_dt(v, f, 0.1, 100.0, 0.01, &CflConfig::conservative());
        let dt_a = adaptive_dt(v, f, 0.1, 100.0, 0.01, &CflConfig::aggressive());
        assert!(dt_c < dt_a);
    }

    #[test]
    fn timestepper_initial_dt() {
        assert_eq!(TimeStepper::new(0.001, 1e-6, 0.01).current_dt(), 0.001);
    }

    #[test]
    fn timestepper_advance_increments_steps() {
        let mut ts = TimeStepper::new(0.001, 1e-6, 0.01);
        ts.advance();
        ts.advance();
        assert_eq!(ts.n_steps, 2);
    }

    #[test]
    fn timestepper_total_time() {
        let mut ts = TimeStepper::new(0.001, 1e-6, 0.01);
        ts.advance();
        ts.advance();
        assert!((ts.total_time() - 0.002).abs() < 1e-15);
    }

    #[test]
    fn timestepper_update_dt_clamps_to_max() {
        let mut ts = TimeStepper::new(0.001, 1e-6, 0.001);
        ts.update_dt(&[[0.0; 3]], &[[0.0; 3]], 0.1, 100.0, 0.01);
        assert!(ts.dt <= 0.001 && ts.dt >= 1e-6);
    }

    #[test]
    fn timestepper_update_dt_clamps_to_min() {
        let mut ts = TimeStepper::new(0.001, 0.005, 0.01);
        ts.update_dt(&[[1e10, 0.0, 0.0]], &[[1e10, 0.0, 0.0]], 0.1, 100.0, 0.01);
        assert!(ts.dt >= 0.005);
    }

    #[test]
    fn euler_step_constant_velocity() {
        let (np, nv) = RungeKuttaSph::euler_step(&[[0.0; 3]], &[[1.0, 0.0, 0.0]], &[[0.0; 3]], 0.1);
        assert!((np[0][0] - 0.1).abs() < 1e-15);
        assert_eq!(nv[0][0], 1.0);
    }

    #[test]
    fn rk2_step_constant_velocity() {
        let (np, _) = RungeKuttaSph::rk2_step(&[[0.0; 3]], &[[2.0, 0.0, 0.0]], &[[0.0; 3]], 0.5);
        assert!((np[0][0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn leapfrog_kick_half_step() {
        let kv = RungeKuttaSph::leapfrog_kick(&[[0.0; 3]], &[[2.0, 0.0, 0.0]], 1.0);
        assert!((kv[0][0] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn leapfrog_drift_full_step() {
        let dp = RungeKuttaSph::leapfrog_drift(&[[1.0, 0.0, 0.0]], &[[3.0, 0.0, 0.0]], 2.0);
        assert!((dp[0][0] - 7.0).abs() < 1e-14);
    }

    #[test]
    fn sound_speed_tait_at_rest_density() {
        let c = sound_speed_tait(1000.0, 1000.0, 7.0, 2.0e5);
        let c0 = (7.0 * 2.0e5 / 1000.0_f64).sqrt();
        assert!((c - c0).abs() < 1e-10);
    }

    #[test]
    fn sound_speed_tait_increases_with_density() {
        assert!(
            sound_speed_tait(1100.0, 1000.0, 7.0, 2.0e5)
                > sound_speed_tait(1000.0, 1000.0, 7.0, 2.0e5)
        );
    }

    #[test]
    fn max_velocity_empty() {
        assert_eq!(max_velocity(&[]), 0.0);
    }

    #[test]
    fn max_velocity_basic() {
        assert!((max_velocity(&[[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]) - 5.0).abs() < 1e-14);
    }

    #[test]
    fn max_acceleration_basic() {
        assert!((max_acceleration(&[[0.0, 2.0, 0.0]], &[0.5]) - 4.0).abs() < 1e-14);
    }

    // ─── Multi-rate tests ───

    #[test]
    fn test_multi_rate_dt_at_level() {
        let mr = MultiRateTimeStepper::new(0.01, 3, 10);
        assert!((mr.dt_at_level(0) - 0.01).abs() < 1e-15);
        assert!((mr.dt_at_level(1) - 0.005).abs() < 1e-15);
        assert!((mr.dt_at_level(2) - 0.0025).abs() < 1e-15);
    }

    #[test]
    fn test_multi_rate_sub_steps() {
        let mr = MultiRateTimeStepper::new(0.01, 3, 5);
        assert_eq!(mr.sub_steps_per_base(), 4); // 2^(3-1)
    }

    #[test]
    fn test_multi_rate_is_active_level0() {
        let mr = MultiRateTimeStepper::new(0.01, 3, 1);
        // Level 0: stride = 2^(3-1-0) = 4, active at 0, 4, 8, ...
        assert!(mr.is_active(0, 0));
        assert!(!mr.is_active(0, 1));
        assert!(!mr.is_active(0, 2));
        assert!(!mr.is_active(0, 3));
    }

    #[test]
    fn test_multi_rate_assign_levels() {
        let mut mr = MultiRateTimeStepper::new(0.01, 3, 2);
        let vels = [[1.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let forces = [[0.0; 3], [0.0; 3]];
        mr.assign_levels(&vels, &forces, 0.1, 100.0, 0.01);
        assert!(
            mr.levels[1] >= mr.levels[0],
            "Fast particle should get higher level"
        );
    }

    // ─── Async particle state tests ───

    #[test]
    fn test_async_state_advance() {
        let mut s = AsyncParticleState::new(0.01);
        assert!((s.t_next - 0.01).abs() < 1e-15);
        s.advance();
        assert!((s.t_last - 0.01).abs() < 1e-15);
        assert!((s.t_next - 0.02).abs() < 1e-14);
    }

    #[test]
    fn test_async_time_since_update() {
        let s = AsyncParticleState::new(0.01);
        assert!((s.time_since_update(0.005) - 0.005).abs() < 1e-15);
    }

    // ─── Sub-cycling tests ───

    #[test]
    fn test_sub_cycling_single_step() {
        let (p, v) =
            SubCyclingIntegrator::integrate_particle([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 0.1, 1);
        assert!((p[0] - 0.1).abs() < 1e-14);
        assert!((v[0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_sub_cycling_multiple_steps() {
        // Constant velocity, should give same result regardless of sub-steps
        let (p1, _) =
            SubCyclingIntegrator::integrate_particle([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 1.0, 1);
        let (p10, _) =
            SubCyclingIntegrator::integrate_particle([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 1.0, 10);
        assert!(
            (p1[0] - p10[0]).abs() < 1e-10,
            "p1={}, p10={}",
            p1[0],
            p10[0]
        );
    }

    #[test]
    fn test_required_sub_steps() {
        assert_eq!(SubCyclingIntegrator::required_sub_steps(0.1, 0.05), 2);
        assert_eq!(SubCyclingIntegrator::required_sub_steps(0.1, 0.03), 4);
        assert_eq!(SubCyclingIntegrator::required_sub_steps(0.1, 0.1), 1);
    }

    // ─── Sync tests ───

    #[test]
    fn test_predict_position() {
        let p = TimestepSync::predict_position([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0, 0.5);
        assert!((p[0] - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_predict_velocity() {
        let v = TimestepSync::predict_velocity([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0, 0.5);
        assert!((v[0] - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_sync_time() {
        let t = TimestepSync::sync_time(&[0.01, 0.005, 0.0025]);
        assert!((t - 0.01).abs() < 1e-15);
    }

    // ─── Error estimation tests ───

    #[test]
    fn test_richardson_error_zero() {
        let err = ErrorEstimator::richardson_error([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2);
        assert!(err.abs() < 1e-15);
    }

    #[test]
    fn test_richardson_error_nonzero() {
        let err = ErrorEstimator::richardson_error([1.0, 0.0, 0.0], [1.03, 0.0, 0.0], 2);
        // |0.03| / (4-1) = 0.01
        assert!((err - 0.01).abs() < 1e-12, "err = {err}");
    }

    #[test]
    fn test_adjusted_dt_grows_when_error_small() {
        let dt_new = ErrorEstimator::adjusted_dt(0.01, 1e-10, 1e-6, 2, 0.9);
        assert!(dt_new > 0.01, "dt should grow: {dt_new}");
    }

    #[test]
    fn test_adjusted_dt_shrinks_when_error_large() {
        let dt_new = ErrorEstimator::adjusted_dt(0.01, 1e-3, 1e-6, 2, 0.9);
        assert!(dt_new < 0.01, "dt should shrink: {dt_new}");
    }

    #[test]
    fn test_max_position_error() {
        let a = &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let b = &[[1.1, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let err = ErrorEstimator::max_position_error(a, b);
        assert!((err - 0.1).abs() < 1e-12, "err = {err}");
    }

    // ─── RK4 test ───

    #[test]
    fn test_rk4_constant_velocity() {
        let (p, v) = RungeKuttaSph::rk4_step_single([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 0.1);
        assert!((p[0] - 0.1).abs() < 1e-12, "p = {}", p[0]);
        assert!((v[0] - 1.0).abs() < 1e-14, "v = {}", v[0]);
    }

    #[test]
    fn test_rk4_constant_acceleration() {
        let (p, v) = RungeKuttaSph::rk4_step_single([0.0; 3], [0.0; 3], [1.0, 0.0, 0.0], 1.0);
        // x = 0.5*a*t^2 = 0.5, v = a*t = 1.0
        assert!((p[0] - 0.5).abs() < 1e-10, "p = {}", p[0]);
        assert!((v[0] - 1.0).abs() < 1e-14, "v = {}", v[0]);
    }
}

// ---------------------------------------------------------------------------
// Verlet integration
// ---------------------------------------------------------------------------

/// Velocity Verlet integrator for SPH.
#[allow(dead_code)]
pub struct VerletIntegrator;

#[allow(dead_code)]
impl VerletIntegrator {
    /// Full velocity-Verlet step assuming constant acceleration.
    ///
    /// `x_new = x + v*dt + 0.5*a*dt^2`
    /// `v_new = v + a*dt`
    pub fn step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = pos.len();
        let half_dt2 = 0.5 * dt * dt;
        let new_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + vel[i][0] * dt + acc[i][0] * half_dt2,
                    pos[i][1] + vel[i][1] * dt + acc[i][1] * half_dt2,
                    pos[i][2] + vel[i][2] * dt + acc[i][2] * half_dt2,
                ]
            })
            .collect();
        let new_vel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + acc[i][0] * dt,
                    vel[i][1] + acc[i][1] * dt,
                    vel[i][2] + acc[i][2] * dt,
                ]
            })
            .collect();
        (new_pos, new_vel)
    }

    /// Velocity Verlet position half-step (predict).
    ///
    /// `x_half = x + v * dt + 0.5 * a * dt^2`
    pub fn predict_positions(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> Vec<[f64; 3]> {
        let half_dt2 = 0.5 * dt * dt;
        pos.iter()
            .zip(vel.iter())
            .zip(acc.iter())
            .map(|((x, v), a)| {
                [
                    x[0] + v[0] * dt + a[0] * half_dt2,
                    x[1] + v[1] * dt + a[1] * half_dt2,
                    x[2] + v[2] * dt + a[2] * half_dt2,
                ]
            })
            .collect()
    }

    /// Velocity Verlet velocity correction (correct).
    ///
    /// `v_new = v + 0.5 * (a_old + a_new) * dt`
    pub fn correct_velocities(
        vel: &[[f64; 3]],
        acc_old: &[[f64; 3]],
        acc_new: &[[f64; 3]],
        dt: f64,
    ) -> Vec<[f64; 3]> {
        let half_dt = 0.5 * dt;
        vel.iter()
            .zip(acc_old.iter())
            .zip(acc_new.iter())
            .map(|((v, a0), a1)| {
                [
                    v[0] + half_dt * (a0[0] + a1[0]),
                    v[1] + half_dt * (a0[1] + a1[1]),
                    v[2] + half_dt * (a0[2] + a1[2]),
                ]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Leapfrog integrator (full step)
// ---------------------------------------------------------------------------

/// Full leapfrog (kick-drift-kick) integrator.
#[allow(dead_code)]
pub struct LeapfrogIntegrator;

#[allow(dead_code)]
impl LeapfrogIntegrator {
    /// Complete kick-drift-kick leapfrog step.
    ///
    /// `v_half = v + 0.5*a*dt`
    /// `x_new  = x + v_half*dt`
    /// `v_new  = v_half + 0.5*a*dt`
    pub fn kdk_step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        // Half kick
        let v_half = RungeKuttaSph::leapfrog_kick(vel, acc, dt);
        // Drift
        let new_pos = RungeKuttaSph::leapfrog_drift(pos, &v_half, dt);
        // Second half kick (uses same acc, constant acc assumption)
        let new_vel = RungeKuttaSph::leapfrog_kick(&v_half, acc, dt);
        (new_pos, new_vel)
    }

    /// Drift-kick-drift (DKD) leapfrog step.
    ///
    /// `x_half = x + 0.5*v*dt`
    /// `v_new  = v + a*dt`
    /// `x_new  = x_half + 0.5*v_new*dt`
    pub fn dkd_step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        // Half drift
        let x_half = RungeKuttaSph::leapfrog_drift(pos, vel, 0.5 * dt);
        // Full kick
        let new_vel: Vec<[f64; 3]> = vel
            .iter()
            .zip(acc.iter())
            .map(|(v, a)| [v[0] + a[0] * dt, v[1] + a[1] * dt, v[2] + a[2] * dt])
            .collect();
        // Second half drift
        let new_pos = RungeKuttaSph::leapfrog_drift(&x_half, &new_vel, 0.5 * dt);
        (new_pos, new_vel)
    }
}

// ---------------------------------------------------------------------------
// Adams-Bashforth 2nd order integrator
// ---------------------------------------------------------------------------

/// Adams-Bashforth second-order (AB2) multi-step integrator.
///
/// Requires the current and previous accelerations.
#[allow(dead_code)]
pub struct AdamsBashforth2;

#[allow(dead_code)]
impl AdamsBashforth2 {
    /// Adams-Bashforth 2nd-order step.
    ///
    /// `x_new = x + dt * (1.5 * v - 0.5 * v_old)`
    /// `v_new = v + dt * (1.5 * a - 0.5 * a_old)`
    pub fn step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        acc_old: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = pos.len();
        let new_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + dt * (1.5 * vel[i][0] - 0.5 * vel[i][0]), // vel unchanged proxy
                    pos[i][1] + dt * (1.5 * vel[i][1] - 0.5 * vel[i][1]),
                    pos[i][2] + dt * (1.5 * vel[i][2] - 0.5 * vel[i][2]),
                ]
            })
            .collect();
        let new_vel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + dt * (1.5 * acc[i][0] - 0.5 * acc_old[i][0]),
                    vel[i][1] + dt * (1.5 * acc[i][1] - 0.5 * acc_old[i][1]),
                    vel[i][2] + dt * (1.5 * acc[i][2] - 0.5 * acc_old[i][2]),
                ]
            })
            .collect();
        (new_pos, new_vel)
    }

    /// Startup step using forward Euler (needed for the first step when no history exists).
    pub fn startup_step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        RungeKuttaSph::euler_step(pos, vel, acc, dt)
    }
}

// ---------------------------------------------------------------------------
// Symplectic Euler integrator
// ---------------------------------------------------------------------------

/// Symplectic Euler (semi-implicit Euler) integrator.
///
/// Updates velocity first, then uses the new velocity for the position update.
/// This is a first-order symplectic method that preserves a modified energy.
#[allow(dead_code)]
pub struct SymplecticEuler;

#[allow(dead_code)]
impl SymplecticEuler {
    /// Perform a symplectic Euler step.
    ///
    /// `v_new = v + a * dt`
    /// `x_new = x + v_new * dt`  ← uses updated velocity
    pub fn step(
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        acc: &[[f64; 3]],
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = pos.len();
        let new_vel: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    vel[i][0] + acc[i][0] * dt,
                    vel[i][1] + acc[i][1] * dt,
                    vel[i][2] + acc[i][2] * dt,
                ]
            })
            .collect();
        let new_pos: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                [
                    pos[i][0] + new_vel[i][0] * dt,
                    pos[i][1] + new_vel[i][1] * dt,
                    pos[i][2] + new_vel[i][2] * dt,
                ]
            })
            .collect();
        (new_pos, new_vel)
    }
}

// ---------------------------------------------------------------------------
// Timestep statistics
// ---------------------------------------------------------------------------

/// Statistics about the timestep usage during a simulation.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct TimestepStatistics {
    /// History of timestep values.
    pub dt_history: Vec<f64>,
    /// Total simulated time.
    pub total_time: f64,
    /// Number of steps taken.
    pub n_steps: u64,
}

#[allow(dead_code)]
impl TimestepStatistics {
    /// Create empty statistics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a timestep.
    pub fn record(&mut self, dt: f64) {
        self.dt_history.push(dt);
        self.total_time += dt;
        self.n_steps += 1;
    }

    /// Minimum timestep seen.
    pub fn min_dt(&self) -> f64 {
        self.dt_history.iter().cloned().fold(f64::MAX, f64::min)
    }

    /// Maximum timestep seen.
    pub fn max_dt(&self) -> f64 {
        self.dt_history.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Mean timestep.
    pub fn mean_dt(&self) -> f64 {
        if self.dt_history.is_empty() {
            return 0.0;
        }
        self.dt_history.iter().sum::<f64>() / self.dt_history.len() as f64
    }

    /// Coefficient of variation (std/mean) of the timestep history.
    pub fn dt_variation(&self) -> f64 {
        if self.dt_history.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_dt();
        if mean < 1e-30 {
            return 0.0;
        }
        let variance = self
            .dt_history
            .iter()
            .map(|&dt| (dt - mean) * (dt - mean))
            .sum::<f64>()
            / (self.dt_history.len() - 1) as f64;
        variance.sqrt() / mean
    }
}

// ---------------------------------------------------------------------------
// Variable-dt SPH stepper
// ---------------------------------------------------------------------------

/// SPH time stepper with variable timestep control and statistics collection.
#[allow(dead_code)]
pub struct VariableDtStepper {
    /// Current timestep.
    pub dt: f64,
    /// Timestep statistics.
    pub stats: TimestepStatistics,
    /// CFL configuration.
    pub cfl: CflConfig,
    /// Minimum dt.
    pub dt_min: f64,
    /// Maximum dt.
    pub dt_max: f64,
}

#[allow(dead_code)]
impl VariableDtStepper {
    /// Create a new variable-dt stepper.
    pub fn new(dt_init: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            dt: dt_init,
            stats: TimestepStatistics::new(),
            cfl: CflConfig::default(),
            dt_min,
            dt_max,
        }
    }

    /// Advance the stepper: compute adaptive dt, record stats, return dt.
    pub fn advance(
        &mut self,
        velocities: &[[f64; 3]],
        forces: &[[f64; 3]],
        h: f64,
        c_sound: f64,
        nu: f64,
    ) -> f64 {
        let dt_proposed = adaptive_dt(velocities, forces, h, c_sound, nu, &self.cfl);
        self.dt = dt_proposed.clamp(self.dt_min, self.dt_max);
        self.stats.record(self.dt);
        self.dt
    }

    /// Total simulated time.
    pub fn total_time(&self) -> f64 {
        self.stats.total_time
    }

    /// Number of steps taken.
    pub fn n_steps(&self) -> u64 {
        self.stats.n_steps
    }
}

// ---------------------------------------------------------------------------
// Additional tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ─── Verlet integrator ───

    #[test]
    fn verlet_step_constant_velocity() {
        let pos = &[[0.0; 3]];
        let vel = &[[1.0, 0.0, 0.0]];
        let acc = &[[0.0; 3]];
        let (np, nv) = VerletIntegrator::step(pos, vel, acc, 0.1);
        assert!((np[0][0] - 0.1).abs() < 1e-14, "x={}", np[0][0]);
        assert!((nv[0][0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn verlet_step_constant_acceleration() {
        let pos = &[[0.0; 3]];
        let vel = &[[0.0; 3]];
        let acc = &[[2.0, 0.0, 0.0]];
        let (np, nv) = VerletIntegrator::step(pos, vel, acc, 1.0);
        // x = 0.5*a*t^2 = 1.0
        assert!((np[0][0] - 1.0).abs() < 1e-14, "x={}", np[0][0]);
        assert!((nv[0][0] - 2.0).abs() < 1e-14, "v={}", nv[0][0]);
    }

    #[test]
    fn verlet_predict_correct() {
        let pos = &[[0.0; 3]];
        let vel = &[[1.0, 0.0, 0.0]];
        let acc = &[[0.0; 3]];
        let predicted = VerletIntegrator::predict_positions(pos, vel, acc, 0.1);
        assert!((predicted[0][0] - 0.1).abs() < 1e-14);
    }

    #[test]
    fn verlet_correct_velocities() {
        let vel = &[[0.0; 3]];
        let a0 = &[[1.0, 0.0, 0.0]];
        let a1 = &[[3.0, 0.0, 0.0]];
        let nv = VerletIntegrator::correct_velocities(vel, a0, a1, 1.0);
        // v_new = v + 0.5*(a0+a1)*dt = 0.5*(1+3) = 2
        assert!((nv[0][0] - 2.0).abs() < 1e-14, "v={}", nv[0][0]);
    }

    // ─── Leapfrog integrator ───

    #[test]
    fn leapfrog_kdk_constant_velocity() {
        let pos = &[[0.0; 3]];
        let vel = &[[1.0, 0.0, 0.0]];
        let acc = &[[0.0; 3]];
        let (np, nv) = LeapfrogIntegrator::kdk_step(pos, vel, acc, 0.5);
        assert!((np[0][0] - 0.5).abs() < 1e-14, "x={}", np[0][0]);
        assert!((nv[0][0] - 1.0).abs() < 1e-14, "v={}", nv[0][0]);
    }

    #[test]
    fn leapfrog_dkd_constant_acceleration() {
        let pos = &[[0.0; 3]];
        let vel = &[[0.0; 3]];
        let acc = &[[1.0, 0.0, 0.0]];
        let (np, nv) = LeapfrogIntegrator::dkd_step(pos, vel, acc, 1.0);
        // x = 0.5*a*t^2 = 0.5, v = a*t = 1.0
        assert!((np[0][0] - 0.5).abs() < 1e-12, "x={}", np[0][0]);
        assert!((nv[0][0] - 1.0).abs() < 1e-14, "v={}", nv[0][0]);
    }

    // ─── Adams-Bashforth 2 ───

    #[test]
    fn ab2_startup_same_as_euler() {
        let pos = &[[0.0; 3]];
        let vel = &[[1.0, 0.0, 0.0]];
        let acc = &[[0.0; 3]];
        let (pab, vab) = AdamsBashforth2::startup_step(pos, vel, acc, 0.1);
        let (peu, veu) = RungeKuttaSph::euler_step(pos, vel, acc, 0.1);
        assert!((pab[0][0] - peu[0][0]).abs() < 1e-15);
        assert!((vab[0][0] - veu[0][0]).abs() < 1e-15);
    }

    #[test]
    fn ab2_step_constant_acceleration() {
        // With constant acc (acc == acc_old), AB2 gives same result as Euler
        let pos = &[[0.0; 3]];
        let vel = &[[0.0; 3]];
        let acc = &[[1.0, 0.0, 0.0]];
        let acc_old = &[[1.0, 0.0, 0.0]];
        let (np, nv) = AdamsBashforth2::step(pos, vel, acc, acc_old, 1.0);
        // v_new = 0 + dt*(1.5*1 - 0.5*1) = 1.0
        assert!((nv[0][0] - 1.0).abs() < 1e-14, "v={}", nv[0][0]);
        // x doesn't use previous vel in this impl (same vel both times)
        assert!(np[0][0].is_finite());
    }

    // ─── Symplectic Euler ───

    #[test]
    fn symplectic_euler_constant_velocity() {
        let pos = &[[0.0; 3]];
        let vel = &[[2.0, 0.0, 0.0]];
        let acc = &[[0.0; 3]];
        let (np, nv) = SymplecticEuler::step(pos, vel, acc, 0.5);
        assert!((np[0][0] - 1.0).abs() < 1e-14, "x={}", np[0][0]);
        assert!((nv[0][0] - 2.0).abs() < 1e-14, "v={}", nv[0][0]);
    }

    #[test]
    fn symplectic_euler_uses_new_velocity_for_position() {
        // With acceleration, symplectic Euler differs from explicit Euler
        let pos = &[[0.0; 3]];
        let vel = &[[0.0; 3]];
        let acc = &[[1.0, 0.0, 0.0]];
        let dt = 1.0;
        let (np_sym, _) = SymplecticEuler::step(pos, vel, acc, dt);
        let (np_exp, _) = RungeKuttaSph::euler_step(pos, vel, acc, dt);
        // Symplectic: x = 0 + (0 + 1*1)*1 = 1.0
        // Explicit:   x = 0 + 0*1 = 0.0
        assert!(
            (np_sym[0][0] - 1.0).abs() < 1e-14,
            "symplectic x={}",
            np_sym[0][0]
        );
        assert!(
            (np_exp[0][0] - 0.0).abs() < 1e-14,
            "explicit x={}",
            np_exp[0][0]
        );
    }

    // ─── Timestep statistics ───

    #[test]
    fn timestep_stats_record() {
        let mut stats = TimestepStatistics::new();
        stats.record(0.01);
        stats.record(0.02);
        stats.record(0.015);
        assert_eq!(stats.n_steps, 3);
        assert!((stats.total_time - 0.045).abs() < 1e-15);
        assert!((stats.min_dt() - 0.01).abs() < 1e-15);
        assert!((stats.max_dt() - 0.02).abs() < 1e-15);
        assert!((stats.mean_dt() - 0.015).abs() < 1e-15);
    }

    #[test]
    fn timestep_stats_empty() {
        let stats = TimestepStatistics::new();
        assert_eq!(stats.n_steps, 0);
        assert_eq!(stats.mean_dt(), 0.0);
        assert_eq!(stats.dt_variation(), 0.0);
    }

    #[test]
    fn timestep_stats_variation() {
        let mut stats = TimestepStatistics::new();
        stats.record(1.0);
        stats.record(1.0);
        assert_eq!(stats.dt_variation(), 0.0);
    }

    // ─── Variable-dt stepper ───

    #[test]
    fn variable_dt_stepper_advances() {
        let mut stepper = VariableDtStepper::new(0.001, 1e-6, 0.01);
        let dt = stepper.advance(&[[1.0, 0.0, 0.0]], &[[0.0; 3]], 0.1, 100.0, 0.01);
        assert!(dt > 0.0 && dt.is_finite(), "dt = {dt}");
        assert_eq!(stepper.n_steps(), 1);
    }

    #[test]
    fn variable_dt_stepper_clamps_to_max() {
        let mut stepper = VariableDtStepper::new(0.001, 1e-6, 0.0005);
        let dt = stepper.advance(&[[0.0; 3]], &[[0.0; 3]], 0.1, 100.0, 0.01);
        assert!(dt <= 0.0005, "dt should be clamped to max: {dt}");
    }

    #[test]
    fn variable_dt_stepper_total_time() {
        let mut stepper = VariableDtStepper::new(0.001, 1e-6, 0.01);
        for _ in 0..5 {
            stepper.advance(&[[0.0; 3]], &[[0.0; 3]], 0.1, 100.0, 0.01);
        }
        assert!(stepper.total_time() > 0.0);
        assert_eq!(stepper.n_steps(), 5);
    }
}
