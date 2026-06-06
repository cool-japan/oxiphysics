// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Time integration schemes for molecular dynamics.
//!
//! These integrators advance atom positions and velocities by one time step.

use crate::atom::AtomSet;

/// Trait for MD time integrators.
///
/// The `step` method advances the system by one time step `dt`.
/// `force_fn` is a callback that computes forces on the atoms
/// (it should clear and recompute `atoms.forces`).
pub trait Integrator: Send + Sync {
    /// Advance the system by one time step.
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet));
}

// ---------------------------------------------------------------------------
// Velocity Verlet
// ---------------------------------------------------------------------------

/// Standard velocity Verlet integrator.
///
/// Algorithm:
/// 1. v(t+dt/2) = v(t) + F(t)/(2m) * dt
/// 2. x(t+dt) = x(t) + v(t+dt/2) * dt
/// 3. Compute F(t+dt)
/// 4. v(t+dt) = v(t+dt/2) + F(t+dt)/(2m) * dt
#[derive(Debug, Clone, Default)]
pub struct VelocityVerlet;

impl VelocityVerlet {
    /// Create a new velocity Verlet integrator.
    pub fn new() -> Self {
        Self
    }
}

impl Integrator for VelocityVerlet {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();
        let half_dt = 0.5 * dt;

        // Step 1: half-kick velocities
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
        }

        // Step 2: full-step positions
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * dt;
        }

        // Step 3: recompute forces at new positions
        force_fn(atoms);

        // Step 4: half-kick velocities with new forces
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
        }
    }
}

// ---------------------------------------------------------------------------
// Leapfrog
// ---------------------------------------------------------------------------

/// Leapfrog integrator.
///
/// Velocities are offset by half a time step from positions.
/// On first call, does an initial half-step kick.
///
/// Algorithm each step:
/// 1. v(t+dt/2) = v(t-dt/2) + F(t)/m * dt
/// 2. x(t+dt) = x(t) + v(t+dt/2) * dt
/// 3. Compute F(t+dt)
#[derive(Debug, Clone)]
pub struct LeapfrogIntegrator {
    /// Whether the first step has been taken.
    initialized: bool,
}

impl Default for LeapfrogIntegrator {
    fn default() -> Self {
        Self::new()
    }
}

impl LeapfrogIntegrator {
    /// Create a new leapfrog integrator.
    pub fn new() -> Self {
        Self { initialized: false }
    }
}

impl Integrator for LeapfrogIntegrator {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();

        if !self.initialized {
            // Initial half-step kick: v(-dt/2) -> v(dt/2) using initial forces
            // Actually we just do v(0) + F(0)/(2m)*dt to get v(dt/2)
            let half_dt = 0.5 * dt;
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
            }
            self.initialized = true;
        } else {
            // Full velocity kick
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * dt);
            }
        }

        // Position update
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * dt;
        }

        // Recompute forces
        force_fn(atoms);
    }
}

// ---------------------------------------------------------------------------
// Symplectic Euler
// ---------------------------------------------------------------------------

/// Symplectic Euler (semi-implicit Euler) integrator.
///
/// Algorithm:
/// 1. v(t+dt) = v(t) + F(t)/m * dt
/// 2. x(t+dt) = x(t) + v(t+dt) * dt    ← uses updated velocity
/// 3. Recompute F(t+dt)
///
/// First-order symplectic method. Less accurate than Verlet but simple.
#[derive(Debug, Clone, Default)]
pub struct SymplecticEuler;

impl SymplecticEuler {
    /// Create a new symplectic Euler integrator.
    pub fn new() -> Self {
        Self
    }
}

impl Integrator for SymplecticEuler {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();

        // Update velocities with current forces
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += atoms.forces[i] * (inv_m * dt);
        }

        // Update positions with new velocities
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * dt;
        }

        // Recompute forces
        force_fn(atoms);
    }
}

// ---------------------------------------------------------------------------
// Position Verlet (Störmer-Verlet)
// ---------------------------------------------------------------------------

/// Position Verlet (Störmer-Verlet) integrator.
///
/// This is the position form of the Verlet algorithm:
///   x(t+dt) = 2·x(t) - x(t-dt) + F(t)/m · dt²
///   v(t) = (x(t+dt) - x(t-dt)) / (2·dt)
///
/// Requires storing the previous positions.
#[derive(Debug, Clone)]
pub struct PositionVerlet {
    /// Previous positions (stored from last step).
    prev_positions: Option<Vec<oxiphysics_core::math::Vec3>>,
}

impl Default for PositionVerlet {
    fn default() -> Self {
        Self::new()
    }
}

impl PositionVerlet {
    /// Create a new position Verlet integrator.
    pub fn new() -> Self {
        Self {
            prev_positions: None,
        }
    }
}

impl Integrator for PositionVerlet {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();
        let dt2 = dt * dt;

        if let Some(ref prev) = self.prev_positions {
            // Store current positions before update
            let cur = atoms.positions.clone();

            // Position Verlet: x_new = 2*x - x_old + a*dt²
            for (i, prev_pos) in prev.iter().enumerate() {
                let inv_m = 1.0 / atoms.masses[i];
                let accel = atoms.forces[i] * inv_m;
                let x_new = atoms.positions[i] * 2.0 - *prev_pos + accel * dt2;
                // Estimate velocity: v = (x_new - x_old) / (2*dt)
                atoms.velocities[i] = (x_new - *prev_pos) * (0.5 / dt);
                atoms.positions[i] = x_new;
            }

            self.prev_positions = Some(cur);
        } else {
            // First step: use velocity Verlet as bootstrap
            let prev = atoms.positions.clone();
            let half_dt = 0.5 * dt;

            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
            }
            for i in 0..n {
                atoms.positions[i] += atoms.velocities[i] * dt;
            }

            force_fn(atoms);

            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
            }

            self.prev_positions = Some(prev);
            return;
        }

        // Recompute forces at new positions
        force_fn(atoms);
    }
}

// ---------------------------------------------------------------------------
// Langevin dynamics integrator (BBK / Brünger-Brooks-Karplus)
// ---------------------------------------------------------------------------

/// Langevin dynamics integrator using the BBK (Brünger-Brooks-Karplus) scheme.
///
/// Adds friction and stochastic forces to the velocity Verlet scheme:
///   v(t+dt/2) = v(t) + \[F(t)/m - γ·v(t)\] · dt/2 + σ·√(dt/2)·R
///   x(t+dt)   = x(t) + v(t+dt/2) · dt
///   F(t+dt) computed
///   v(t+dt)   = \[v(t+dt/2) + F(t+dt)/m · dt/2 + σ·√(dt/2)·R'\] / (1 + γ·dt/2)
///
/// where γ is the friction coefficient and σ = √(2·γ·k_B·T/m).
///
/// This implementation uses `rand 0.9.1` API (`rng()`, `random_range`).
#[derive(Debug, Clone)]
pub struct LangevinIntegrator {
    /// Friction coefficient γ (s⁻¹).
    pub gamma: f64,
    /// Target temperature (K).
    pub temperature: f64,
}

/// Boltzmann constant for integrator module (J/K).
const KB_INT: f64 = 1.380_649e-23;

impl LangevinIntegrator {
    /// Create a new Langevin dynamics integrator.
    ///
    /// # Arguments
    /// * `gamma` – friction coefficient (s⁻¹)
    /// * `temperature` – target temperature (K)
    pub fn new(gamma: f64, temperature: f64) -> Self {
        Self { gamma, temperature }
    }

    /// Generate a Gaussian random number using Box-Muller transform.
    fn gaussian(&self) -> f64 {
        use rand::RngExt;
        let mut rng = rand::rng();
        let u1: f64 = rng.random_range(1e-10..1.0_f64);
        let u2: f64 = rng.random_range(0.0..std::f64::consts::TAU);
        (-2.0 * u1.ln()).sqrt() * u2.cos()
    }
}

impl Integrator for LangevinIntegrator {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();
        let half_dt = 0.5 * dt;
        let sqrt_half_dt = half_dt.sqrt();

        // Half-step velocity update with friction and noise
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            let sigma = (2.0 * self.gamma * KB_INT * self.temperature * inv_m).sqrt();

            for d in 0..3 {
                let noise = sigma * sqrt_half_dt * self.gaussian();
                let accel = match d {
                    0 => atoms.forces[i].x * inv_m,
                    1 => atoms.forces[i].y * inv_m,
                    _ => atoms.forces[i].z * inv_m,
                };
                let v_old = match d {
                    0 => atoms.velocities[i].x,
                    1 => atoms.velocities[i].y,
                    _ => atoms.velocities[i].z,
                };
                let v_half = v_old + (accel - self.gamma * v_old) * half_dt + noise;
                match d {
                    0 => atoms.velocities[i].x = v_half,
                    1 => atoms.velocities[i].y = v_half,
                    _ => atoms.velocities[i].z = v_half,
                }
            }
        }

        // Full position update
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * dt;
        }

        // Recompute forces
        force_fn(atoms);

        // Second half-step velocity update
        let denom = 1.0 / (1.0 + self.gamma * half_dt);
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            let sigma = (2.0 * self.gamma * KB_INT * self.temperature * inv_m).sqrt();

            for d in 0..3 {
                let noise = sigma * sqrt_half_dt * self.gaussian();
                let accel = match d {
                    0 => atoms.forces[i].x * inv_m,
                    1 => atoms.forces[i].y * inv_m,
                    _ => atoms.forces[i].z * inv_m,
                };
                let v_half = match d {
                    0 => atoms.velocities[i].x,
                    1 => atoms.velocities[i].y,
                    _ => atoms.velocities[i].z,
                };
                let v_new = (v_half + accel * half_dt + noise) * denom;
                match d {
                    0 => atoms.velocities[i].x = v_new,
                    1 => atoms.velocities[i].y = v_new,
                    _ => atoms.velocities[i].z = v_new,
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// RESPA (multi-timestep)
// ---------------------------------------------------------------------------

/// Reference System Propagator Algorithm (RESPA) multi-timestep integrator.
///
/// RESPA splits forces into fast (inner) and slow (outer) components and
/// integrates them with different time steps to improve efficiency.
///
/// # Reference
/// Tuckerman, Berne & Martyna, J. Chem. Phys. 97(3), 1990 (1992).
///
/// # Current status
/// This is a stub implementation that applies `n_inner` velocity-Verlet inner
/// steps with step size `dt / n_inner` for the fast forces, while the slow
/// forces are applied once per outer step.  A full RESPA implementation would
/// separate the force computation into fast/slow components; here both the
/// single supplied `force_fn` callback covers all forces (equivalent to a
/// plain velocity-Verlet but with a sub-cycling loop).
#[derive(Debug, Clone)]
pub struct RespaIntegrator {
    /// Number of inner (fast) steps per outer step.
    pub n_inner: usize,
}

impl RespaIntegrator {
    /// Create a new RESPA integrator.
    ///
    /// # Arguments
    /// * `n_inner` – Number of inner time steps per outer step (must be >= 1).
    pub fn new(n_inner: usize) -> Self {
        Self {
            n_inner: n_inner.max(1),
        }
    }
}

impl Integrator for RespaIntegrator {
    /// Advance the system by one outer time step using `n_inner` inner steps.
    ///
    /// Each inner step is a standard velocity-Verlet step with
    /// `dt_inner = dt / n_inner`.
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let dt_inner = dt / self.n_inner as f64;
        let half_dt_inner = 0.5 * dt_inner;
        let n = atoms.len();

        for _inner in 0..self.n_inner {
            // Half-kick velocities
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt_inner);
            }

            // Full-step positions
            for i in 0..n {
                atoms.positions[i] += atoms.velocities[i] * dt_inner;
            }

            // Recompute forces
            force_fn(atoms);

            // Second half-kick
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt_inner);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Predictor-Corrector (Adams-Bashforth / Adams-Moulton 2nd order)
// ---------------------------------------------------------------------------

/// Second-order predictor-corrector integrator (Adams-Bashforth predictor +
/// Adams-Moulton corrector).
///
/// **Predictor** (AB2):
///   x*(t+dt) = x(t) + dt · \[3/2 · v(t) − 1/2 · v(t−dt)\]
///   v*(t+dt) = v(t) + dt · \[3/2 · a(t) − 1/2 · a(t−dt)\]
///
/// **Corrector** (AM2):
///   x(t+dt) = x(t) + dt/2 · \[v*(t+dt) + v(t)\]
///   v(t+dt) = v(t) + dt/2 · \[a*(t+dt) + a(t)\]
///
/// Bootstraps via velocity Verlet on the first step.
#[derive(Debug, Clone)]
pub struct PredictorCorrector {
    prev_forces: Option<Vec<oxiphysics_core::math::Vec3>>,
    prev_velocities: Option<Vec<oxiphysics_core::math::Vec3>>,
}

impl Default for PredictorCorrector {
    fn default() -> Self {
        Self::new()
    }
}

impl PredictorCorrector {
    /// Create a new predictor-corrector integrator.
    pub fn new() -> Self {
        Self {
            prev_forces: None,
            prev_velocities: None,
        }
    }
}

impl Integrator for PredictorCorrector {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();

        if self.prev_forces.is_none() {
            // Bootstrap with velocity Verlet
            let half_dt = 0.5 * dt;
            let prev_f = atoms.forces.clone();
            let prev_v = atoms.velocities.clone();

            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
            }
            for i in 0..n {
                atoms.positions[i] += atoms.velocities[i] * dt;
            }
            force_fn(atoms);
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
            }

            self.prev_forces = Some(prev_f);
            self.prev_velocities = Some(prev_v);
            return;
        }

        // Safety: prev_forces/prev_velocities are always Some after the bootstrap path above
        let Some(prev_f) = self.prev_forces.as_ref() else {
            return;
        };
        let Some(prev_v) = self.prev_velocities.as_ref() else {
            return;
        };

        // AB2 predictor
        let mut pos_pred = atoms.positions.clone();
        let mut vel_pred = atoms.velocities.clone();
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            let a_cur = atoms.forces[i] * inv_m;
            let a_prv = prev_f[i] * inv_m;
            vel_pred[i] = atoms.velocities[i] + (a_cur * 1.5 - a_prv * 0.5) * dt;
            pos_pred[i] = atoms.positions[i] + (atoms.velocities[i] * 1.5 - prev_v[i] * 0.5) * dt;
        }

        // Save current state for next step
        let saved_f = atoms.forces.clone();
        let saved_v = atoms.velocities.clone();
        let saved_pos = atoms.positions.clone();

        // Apply predictor positions, compute forces
        atoms.positions[..n].copy_from_slice(&pos_pred[..n]);
        atoms.velocities[..n].copy_from_slice(&vel_pred[..n]);
        force_fn(atoms);

        // AM2 corrector
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            let a_pred = atoms.forces[i] * inv_m;
            let a_cur = saved_f[i] * inv_m;
            atoms.velocities[i] = saved_v[i] + (a_pred + a_cur) * (0.5 * dt);
        }
        // Position corrector: base on saved positions x(t), not the predicted position
        for i in 0..n {
            atoms.positions[i] = saved_pos[i] + (atoms.velocities[i] + saved_v[i]) * (0.5 * dt);
        }
        force_fn(atoms);

        self.prev_forces = Some(saved_f);
        self.prev_velocities = Some(saved_v);
    }
}

// ---------------------------------------------------------------------------
// Operator Splitting (Strang splitting) integrator
// ---------------------------------------------------------------------------

/// Strang-splitting (symmetric operator splitting) integrator.
///
/// Decomposes the time evolution operator as:
///   exp(L·dt) ≈ exp(A·dt/2) · exp(B·dt) · exp(A·dt/2)
///
/// Where A is the "kick" (velocity update) and B is the "drift" (position update).
/// This is algebraically equivalent to velocity Verlet and preserves the
/// symplectic structure.
///
/// This implementation applies an additional split by halving position updates
/// symmetrically, making it explicit for demonstration purposes.
#[derive(Debug, Clone, Default)]
pub struct StrangSplitting;

impl StrangSplitting {
    /// Create a new Strang-splitting integrator.
    pub fn new() -> Self {
        Self
    }
}

impl Integrator for StrangSplitting {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();
        let half_dt = 0.5 * dt;

        // A(dt/2): half-kick
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
        }

        // B(dt/2): half-drift
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * half_dt;
        }

        // Recompute forces at mid-step position (optional sub-step)
        force_fn(atoms);

        // B(dt/2): second half-drift
        for i in 0..n {
            atoms.positions[i] += atoms.velocities[i] * half_dt;
        }

        // Recompute forces at full-step position
        force_fn(atoms);

        // A(dt/2): second half-kick
        for i in 0..n {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt);
        }
    }
}

// ---------------------------------------------------------------------------
// 4th-order Ruth symplectic integrator
// ---------------------------------------------------------------------------

/// 4th-order Ruth symplectic integrator.
///
/// Uses the classical Ruth coefficients:
///   c = \[1/(2(2−2^(1/3))), (1−2^(1/3))/(2(2−2^(1/3))), (1−2^(1/3))/(2(2−2^(1/3))), 1/(2(2−2^(1/3)))\]
///   d = \[1/(2−2^(1/3)), −2^(1/3)/(2−2^(1/3)), 1/(2−2^(1/3)), 0\]
///
/// Reference: Forest & Ruth (1990), Physica D 43, 105.
///
/// This is a composition of 4 velocity Verlet steps with specific sub-step sizes.
/// It achieves 4th-order accuracy for separable Hamiltonians.
#[derive(Debug, Clone)]
pub struct Ruth4Integrator;

impl Default for Ruth4Integrator {
    fn default() -> Self {
        Self::new()
    }
}

impl Ruth4Integrator {
    /// Create a new 4th-order Ruth integrator.
    pub fn new() -> Self {
        Self
    }

    /// Ruth 4th-order coefficients.
    ///
    /// Returns (c_coeffs, d_coeffs) where each has 4 elements.
    fn coefficients() -> ([f64; 4], [f64; 4]) {
        let cbrt2 = 2.0_f64.powf(1.0 / 3.0);
        let w1 = 1.0 / (2.0 - cbrt2);
        let w0 = -cbrt2 * w1;
        // c[i]: position update coefficients
        let c = [w1 / 2.0, (w0 + w1) / 2.0, (w0 + w1) / 2.0, w1 / 2.0];
        // d[i]: velocity (kick) coefficients
        let d = [w1, w0, w1, 0.0];
        (c, d)
    }
}

impl Integrator for Ruth4Integrator {
    fn step(&mut self, atoms: &mut AtomSet, dt: f64, force_fn: &mut dyn FnMut(&mut AtomSet)) {
        let n = atoms.len();
        let (c, d) = Self::coefficients();

        for stage in 0..4 {
            // Drift
            let c_dt = c[stage] * dt;
            for i in 0..n {
                atoms.positions[i] += atoms.velocities[i] * c_dt;
            }

            // Compute forces at new position
            if d[stage].abs() > 1e-14 {
                force_fn(atoms);
                // Kick
                let d_dt = d[stage] * dt;
                for i in 0..n {
                    let inv_m = 1.0 / atoms.masses[i];
                    atoms.velocities[i] += atoms.forces[i] * (inv_m * d_dt);
                }
            }
        }
        // Final force evaluation at new positions
        force_fn(atoms);
    }
}

// ---------------------------------------------------------------------------
// Integration Error Monitor
// ---------------------------------------------------------------------------

/// Monitors integration error by tracking energy drift and step-size
/// recommendations.
///
/// Accumulates statistics about the energy at each step and can suggest
/// whether the time step should be reduced or increased.
#[derive(Debug, Clone)]
pub struct IntegrationErrorMonitor {
    /// Initial total energy (set at first call to `record`).
    pub e0: Option<f64>,
    /// Maximum observed relative energy drift.
    pub max_drift: f64,
    /// Current step count.
    pub step_count: u64,
    /// Running mean of relative energy drift.
    pub mean_drift: f64,
    /// Energy tolerance threshold for step-size advice.
    pub tolerance: f64,
}

impl IntegrationErrorMonitor {
    /// Create a new monitor with the given energy drift tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            e0: None,
            max_drift: 0.0,
            step_count: 0,
            mean_drift: 0.0,
            tolerance,
        }
    }

    /// Record the current total energy.
    ///
    /// Returns the current relative drift |E − E0| / |E0|.
    pub fn record(&mut self, total_energy: f64) -> f64 {
        if self.e0.is_none() {
            self.e0 = Some(total_energy);
            return 0.0;
        }
        // Safety: guaranteed Some by the is_none() check above
        let Some(e0) = self.e0 else {
            return 0.0;
        };
        let drift = if e0.abs() < 1e-30 {
            0.0
        } else {
            (total_energy - e0).abs() / e0.abs()
        };
        if drift > self.max_drift {
            self.max_drift = drift;
        }
        // Online mean update
        self.step_count += 1;
        self.mean_drift += (drift - self.mean_drift) / self.step_count as f64;
        drift
    }

    /// Return `true` if the current max drift exceeds the tolerance.
    pub fn is_drifting(&self) -> bool {
        self.max_drift > self.tolerance
    }

    /// Suggest a scale factor for the time step based on observed drift.
    ///
    /// If drift < tolerance/10, returns 1.1 (increase dt by 10%).
    /// If drift > tolerance, returns 0.5 (halve dt).
    /// Otherwise returns 1.0 (keep dt).
    pub fn suggest_dt_scale(&self) -> f64 {
        if self.max_drift > self.tolerance {
            0.5
        } else if self.max_drift < self.tolerance / 10.0 {
            1.1
        } else {
            1.0
        }
    }

    /// Reset the monitor (e.g. after re-thermalising).
    pub fn reset(&mut self) {
        self.e0 = None;
        self.max_drift = 0.0;
        self.step_count = 0;
        self.mean_drift = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Multi-timestep RESPA with separate fast/slow force functions
// ---------------------------------------------------------------------------

/// Improved RESPA integrator that accepts separate fast and slow force
/// functions.
///
/// The slow forces are applied once at the outer step (split across the
/// two half-kicks), while the fast forces are recomputed at every inner step.
/// This faithfully implements the RESPA scheme from Tuckerman *et al.* (1992).
#[derive(Debug, Clone)]
pub struct RespaImproved {
    /// Number of inner steps per outer step.
    pub n_inner: usize,
}

impl RespaImproved {
    /// Create a new improved RESPA integrator.
    pub fn new(n_inner: usize) -> Self {
        Self {
            n_inner: n_inner.max(1),
        }
    }

    /// Advance the system by one outer step using separate fast/slow force functions.
    ///
    /// # Arguments
    /// * `atoms`        – atom set to advance
    /// * `dt`           – outer time step
    /// * `fast_force_fn`– callback for fast (bonded) forces; cleared and added to `atoms.forces`
    /// * `slow_force_fn`– callback for slow (non-bonded) forces; cleared and added to `atoms.forces`
    pub fn step_split(
        &mut self,
        atoms: &mut AtomSet,
        dt: f64,
        fast_force_fn: &mut dyn FnMut(&mut AtomSet),
        slow_force_fn: &mut dyn FnMut(&mut AtomSet),
    ) {
        let n = atoms.len();
        let dt_inner = dt / self.n_inner as f64;
        let half_dt_inner = 0.5 * dt_inner;
        let half_dt_outer = 0.5 * dt;

        // Compute slow forces (outer)
        slow_force_fn(atoms);
        let slow_forces = atoms.forces.clone();

        // Outer half-kick with slow forces
        for (i, sf) in slow_forces.iter().enumerate() {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += *sf * (inv_m * half_dt_outer);
        }

        // Inner loop: velocity Verlet with fast forces
        for _inner in 0..self.n_inner {
            fast_force_fn(atoms);

            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt_inner);
            }
            for i in 0..n {
                atoms.positions[i] += atoms.velocities[i] * dt_inner;
            }
            fast_force_fn(atoms);
            for i in 0..n {
                let inv_m = 1.0 / atoms.masses[i];
                atoms.velocities[i] += atoms.forces[i] * (inv_m * half_dt_inner);
            }
        }

        // Outer half-kick with slow forces at new position
        slow_force_fn(atoms);
        let slow_forces_new = atoms.forces.clone();
        for (i, sfn) in slow_forces_new.iter().enumerate() {
            let inv_m = 1.0 / atoms.masses[i];
            atoms.velocities[i] += *sfn * (inv_m * half_dt_outer);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::math::Vec3;

    /// NVE energy conservation over 100 steps (required test).
    ///
    /// Two LJ atoms at equilibrium distance with small initial velocities;
    /// total energy must be conserved within 1%.
    #[test]
    fn test_nve_energy_conservation_100_steps() {
        use crate::forcefield::{ForceField, PairForceField};
        use crate::neighbor::PeriodicBox;
        use crate::potential::LennardJones;

        let sigma = 1.0_f64;
        let r_eq = 2.0_f64.powf(1.0 / 6.0) * sigma;
        let mut atoms = AtomSet::new();
        // Place two atoms near equilibrium with small velocities
        atoms.add_atom(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.05, 0.0, 0.0),
            1.0,
            0.0,
            0,
        );
        atoms.add_atom(
            Vec3::new(r_eq, 0.0, 0.0),
            Vec3::new(-0.05, 0.0, 0.0),
            1.0,
            0.0,
            0,
        );

        let pbox = PeriodicBox::cubic(20.0);
        let mut ff = PairForceField::new();
        ff.add_interaction(0, 0, LennardJones::new(1.0, sigma, 5.0));

        // Compute initial forces
        atoms.clear_forces();
        ff.compute_forces(&mut atoms, &pbox);

        let ke0 = atoms.kinetic_energy();
        let pe0 = {
            atoms.clear_forces();

            ff.compute_forces(&mut atoms, &pbox)
        };
        let e0 = ke0 + pe0;

        let mut integrator = VelocityVerlet::new();
        let dt = 0.002;

        for _ in 0..100 {
            let ff_ref = &ff;
            let pbox_ref = &pbox;
            integrator.step(&mut atoms, dt, &mut |a: &mut AtomSet| {
                a.clear_forces();
                ff_ref.compute_forces(a, pbox_ref);
            });
        }

        let ke_f = atoms.kinetic_energy();
        atoms.clear_forces();
        let pe_f = ff.compute_forces(&mut atoms, &pbox);
        let e_f = ke_f + pe_f;

        let drift = (e_f - e0).abs();
        let scale = e0.abs().max(1.0);
        assert!(
            drift / scale < 0.01,
            "NVE energy drift too large: {drift} over scale {scale} (e0={e0}, ef={e_f})"
        );
    }

    /// RESPA with n_inner=1 should behave identically to plain VelocityVerlet.
    #[test]
    fn test_respa_n1_equivalent_to_verlet() {
        let mut atoms_vv = AtomSet::new();
        let mut atoms_rp = AtomSet::new();

        atoms_vv.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
        atoms_rp.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn_vv = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        let mut force_fn_rp = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };

        force_fn_vv(&mut atoms_vv);
        force_fn_rp(&mut atoms_rp);

        let dt = 0.01;
        let mut vv = VelocityVerlet::new();
        let mut rp = RespaIntegrator::new(1);

        for _ in 0..50 {
            vv.step(&mut atoms_vv, dt, &mut force_fn_vv);
            rp.step(&mut atoms_rp, dt, &mut force_fn_rp);
        }

        let diff = (atoms_vv.positions[0] - atoms_rp.positions[0]).norm();
        assert!(
            diff < 1e-10,
            "RESPA(n=1) should match VelocityVerlet; position diff={diff}"
        );
    }

    /// RESPA oscillator: with n_inner=4 should still give reasonable accuracy.
    #[test]
    fn test_respa_harmonic_oscillator_accuracy() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = RespaIntegrator::new(4);
        let dt = 0.001;
        let n_steps = 6283; // ~one period

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        // After ~one period should return close to x=1
        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.05,
            "RESPA harmonic x = {} after one period",
            atoms.positions[0].x
        );
    }

    /// Simple harmonic oscillator: F = -k*x, with k=1, m=1.
    /// x(t) = cos(t), v(t) = -sin(t) for x0=1, v0=0.
    #[test]
    fn test_velocity_verlet_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };

        // Initial force
        force_fn(&mut atoms);

        let mut integrator = VelocityVerlet::new();
        let dt = 0.001;
        let n_steps = 6283; // ~one period

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        // After one period, should return close to (1, 0, 0)
        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.01,
            "x = {}",
            atoms.positions[0].x
        );
    }

    // --- Symplectic Euler tests ---

    #[test]
    fn test_symplectic_euler_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = SymplecticEuler::new();
        let dt = 0.001;
        let n_steps = 6283;

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        // Symplectic Euler is first-order, so less accurate, but should still
        // oscillate around the equilibrium. Check it's within 0.1.
        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.1,
            "Symplectic Euler x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_symplectic_euler_energy_bounded() {
        // Symplectic methods should conserve energy (within bounds)
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        // Initial energy: KE=0, PE=0.5*1*1=0.5
        let e0 = 0.5 * k * atoms.positions[0].x * atoms.positions[0].x;

        let mut integrator = SymplecticEuler::new();
        let dt = 0.001;

        for _ in 0..1000 {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        let ke = 0.5 * atoms.velocities[0].x * atoms.velocities[0].x;
        let pe = 0.5 * k * atoms.positions[0].x * atoms.positions[0].x;
        let e_final = ke + pe;
        let drift = (e_final - e0).abs() / e0;
        assert!(
            drift < 0.05,
            "Symplectic Euler energy drift {drift:.4} should be < 5%"
        );
    }

    // --- Position Verlet tests ---

    #[test]
    fn test_position_verlet_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = PositionVerlet::new();
        let dt = 0.001;
        let n_steps = 6283;

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.02,
            "Position Verlet x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_position_verlet_energy_conservation() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let e0 = 0.5 * k * atoms.positions[0].x.powi(2);
        let mut integrator = PositionVerlet::new();
        let dt = 0.001;

        for _ in 0..1000 {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        let ke = 0.5 * atoms.velocities[0].norm_squared();
        let pe = 0.5 * k * atoms.positions[0].x.powi(2);
        let drift = ((ke + pe) - e0).abs() / e0;
        assert!(drift < 0.02, "Position Verlet energy drift = {drift:.4}");
    }

    // --- Langevin dynamics tests ---

    #[test]
    fn test_langevin_runs_without_crash() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = LangevinIntegrator::new(1.0, 300.0);
        let dt = 0.001;

        for _ in 0..100 {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        // Just check it produced finite positions and velocities
        assert!(
            atoms.positions[0].x.is_finite(),
            "Langevin position should be finite"
        );
        assert!(
            atoms.velocities[0].x.is_finite(),
            "Langevin velocity should be finite"
        );
    }

    #[test]
    fn test_langevin_friction_damps() {
        // High friction should damp the motion
        let mut atoms = AtomSet::new();
        atoms.add_atom(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(10.0, 0.0, 0.0),
            1.0,
            0.0,
            0,
        );

        // No forces, just free particle with friction
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = Vec3::zeros();
        };
        force_fn(&mut atoms);

        // High friction, zero temperature → should damp velocity
        // Use gamma * dt/2 >> 1 so the BBK denominator 1/(1+gamma*dt/2) is small
        let mut integrator = LangevinIntegrator::new(1e4, 0.0);
        let dt = 0.01;

        for _ in 0..500 {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        // Velocity should be heavily damped
        let v = atoms.velocities[0].norm();
        assert!(
            v < 1.0,
            "Langevin with high friction should damp velocity, got |v|={v}"
        );
    }

    // --- RESPA with different inner steps ---

    #[test]
    fn test_respa_n4_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = RespaIntegrator::new(4);
        let dt = 0.004; // outer step, inner = 0.001
        let n_steps = 1571; // ~one period (6283 inner steps)

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.05,
            "RESPA(4) harmonic x = {} after one period",
            atoms.positions[0].x
        );
    }

    // --- PredictorCorrector tests ---

    #[test]
    fn test_predictor_corrector_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = PredictorCorrector::new();
        let dt = 0.001;
        let n_steps = 6283;

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.05,
            "PredictorCorrector x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_predictor_corrector_finite() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.1, 0.0, 0.0),
            1.0,
            0.0,
            0,
        );

        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0];
        };
        force_fn(&mut atoms);

        let mut integrator = PredictorCorrector::new();
        for _ in 0..200 {
            integrator.step(&mut atoms, 0.001, &mut force_fn);
        }
        assert!(atoms.positions[0].x.is_finite());
        assert!(atoms.velocities[0].x.is_finite());
    }

    // --- Strang Splitting tests ---

    #[test]
    fn test_strang_splitting_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = StrangSplitting::new();
        let dt = 0.001;
        let n_steps = 6283;

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.05,
            "StrangSplitting x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_strang_splitting_energy_bounded() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let e0 = 0.5 * atoms.positions[0].x.powi(2);
        let mut integrator = StrangSplitting::new();

        for _ in 0..1000 {
            integrator.step(&mut atoms, 0.001, &mut force_fn);
        }

        let ke = 0.5 * atoms.velocities[0].norm_squared();
        let pe = 0.5 * atoms.positions[0].x.powi(2);
        let drift = ((ke + pe) - e0).abs() / e0;
        assert!(
            drift < 0.05,
            "Strang splitting energy drift {drift:.4} > 5%"
        );
    }

    // --- Ruth 4th-order symplectic tests ---

    #[test]
    fn test_ruth4_harmonic_accuracy() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let mut integrator = Ruth4Integrator::new();
        let dt = 0.01; // larger dt thanks to 4th order
        let n_steps = 629; // ~one period

        for _ in 0..n_steps {
            integrator.step(&mut atoms, dt, &mut force_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.02,
            "Ruth4 harmonic x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_ruth4_energy_conservation() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut force_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * k;
        };
        force_fn(&mut atoms);

        let e0 = 0.5 * k * atoms.positions[0].x.powi(2);
        let mut integrator = Ruth4Integrator::new();

        for _ in 0..500 {
            integrator.step(&mut atoms, 0.01, &mut force_fn);
        }

        let ke = 0.5 * atoms.velocities[0].norm_squared();
        let pe = 0.5 * k * atoms.positions[0].x.powi(2);
        let drift = ((ke + pe) - e0).abs() / e0;
        assert!(drift < 0.01, "Ruth4 energy drift {drift:.6} should be < 1%");
    }

    #[test]
    fn test_ruth4_coefficients_sum_to_one() {
        // Drift coefficients c should sum to 1 (covers a full step)
        let (c, d) = Ruth4Integrator::coefficients();
        let c_sum: f64 = c.iter().sum();
        let d_sum: f64 = d.iter().sum();
        assert!((c_sum - 1.0).abs() < 1e-12, "c sum = {c_sum}");
        assert!((d_sum - 1.0).abs() < 1e-12, "d sum = {d_sum}");
    }

    // --- Integration Error Monitor tests ---

    #[test]
    fn test_error_monitor_zero_drift() {
        let mut monitor = IntegrationErrorMonitor::new(0.01);
        monitor.record(100.0);
        let drift = monitor.record(100.0);
        assert!(drift < 1e-12);
        assert!(!monitor.is_drifting());
    }

    #[test]
    fn test_error_monitor_detects_drift() {
        let mut monitor = IntegrationErrorMonitor::new(0.01);
        monitor.record(100.0);
        monitor.record(102.0); // 2% drift
        assert!(monitor.is_drifting());
        assert_eq!(monitor.suggest_dt_scale(), 0.5);
    }

    #[test]
    fn test_error_monitor_suggests_increase() {
        let mut monitor = IntegrationErrorMonitor::new(0.01);
        monitor.record(100.0);
        monitor.record(100.0 + 1e-5); // tiny drift
        assert_eq!(monitor.suggest_dt_scale(), 1.1);
    }

    #[test]
    fn test_error_monitor_reset() {
        let mut monitor = IntegrationErrorMonitor::new(0.01);
        monitor.record(100.0);
        monitor.record(200.0); // large drift
        monitor.reset();
        assert!(monitor.e0.is_none());
        assert!(monitor.max_drift < 1e-30);
    }

    #[test]
    fn test_error_monitor_step_count() {
        let mut monitor = IntegrationErrorMonitor::new(0.01);
        monitor.record(1.0); // initialises e0
        monitor.record(1.0);
        monitor.record(1.0);
        assert_eq!(monitor.step_count, 2);
    }

    // --- RespaImproved tests ---

    #[test]
    fn test_respa_improved_harmonic() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

        let k = 1.0;
        let mut fast_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * (k * 0.5);
        };
        let mut slow_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * (k * 0.5);
        };

        fast_fn(&mut atoms);

        let mut integrator = RespaImproved::new(4);
        let dt = 0.004;
        let n_steps = 1571;

        for _ in 0..n_steps {
            integrator.step_split(&mut atoms, dt, &mut fast_fn, &mut slow_fn);
        }

        assert!(
            (atoms.positions[0].x - 1.0).abs() < 0.05,
            "RespaImproved harmonic x = {} after one period",
            atoms.positions[0].x
        );
    }

    #[test]
    fn test_respa_improved_finite() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.2, 0.0, 0.0),
            2.0,
            0.0,
            0,
        );

        let mut fast_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * 0.5;
        };
        let mut slow_fn = |atoms: &mut AtomSet| {
            atoms.forces[0] = -atoms.positions[0] * 0.5;
        };

        fast_fn(&mut atoms);

        let mut integrator = RespaImproved::new(2);
        for _ in 0..100 {
            integrator.step_split(&mut atoms, 0.002, &mut fast_fn, &mut slow_fn);
        }

        assert!(atoms.positions[0].x.is_finite());
    }
}
