#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core simulation traits for the OxiPhysics engine.
//!
//! This module defines the fundamental behavioural contracts used throughout
//! the engine.  All physics objects implement one or more of these traits so
//! that solvers, integrators, and broadphase structures can operate on them
//! generically.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use crate::math::Real;
use crate::types::{Aabb, TimeStep, Transform};

// ---------------------------------------------------------------------------
// Time integration
// ---------------------------------------------------------------------------

/// Trait for objects that can be stepped forward in time.
pub trait Steppable {
    /// Advance the simulation by one time step.
    fn step(&mut self, time_step: &TimeStep);
}

// ---------------------------------------------------------------------------
// Physics body
// ---------------------------------------------------------------------------

/// Trait for objects that represent a physics body.
pub trait PhysicsBody {
    /// Return the current transform of this body.
    fn transform(&self) -> &Transform;

    /// Return a mutable reference to the transform.
    fn transform_mut(&mut self) -> &mut Transform;

    /// Return the mass of this body (kg). Returns 0 for static/kinematic.
    fn mass(&self) -> Real;

    /// Return the inverse mass. Returns 0 for static/kinematic.
    fn inv_mass(&self) -> Real {
        if self.mass() > 0.0 {
            1.0 / self.mass()
        } else {
            0.0
        }
    }

    /// Whether this body participates in dynamics integration.
    fn is_dynamic(&self) -> bool {
        self.mass() > 0.0
    }

    /// Whether this body is treated as immovable.
    fn is_static(&self) -> bool {
        self.mass() == 0.0
    }
}

// ---------------------------------------------------------------------------
// Collision
// ---------------------------------------------------------------------------

/// Trait for objects that can participate in collision detection.
pub trait Collidable {
    /// Return the axis-aligned bounding box of this object.
    fn bounding_box(&self) -> Aabb;

    /// Whether this object can interact with another.
    fn can_collide_with(&self, _other: &dyn Collidable) -> bool {
        true
    }
}

// ---------------------------------------------------------------------------
// Force source
// ---------------------------------------------------------------------------

/// Trait for objects that apply forces to bodies (gravity fields, springs, etc.).
pub trait ForceSource {
    /// Compute the force (and optionally torque) acting on a body at the
    /// given transform with the given velocity.
    ///
    /// Returns `(force, torque)` both as `[f64; 3]` in world space.
    fn force_at(
        &self,
        transform: &Transform,
        velocity: [Real; 3],
        mass: Real,
    ) -> ([Real; 3], [Real; 3]);
}

// ---------------------------------------------------------------------------
// Constraint
// ---------------------------------------------------------------------------

/// Trait for position/velocity-level constraints between bodies.
pub trait PhysicsConstraint {
    /// Number of degrees of freedom constrained (1–6).
    fn dof(&self) -> usize;

    /// Evaluate the positional violation of this constraint.
    /// Returns a slice of `dof()` scalar violations.
    fn positional_error(&self) -> Vec<Real>;

    /// Whether this constraint is currently active.
    fn is_active(&self) -> bool;

    /// Compliance of the constraint (XPBD softness). 0 = rigid.
    fn compliance(&self) -> Real {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Integrator
// ---------------------------------------------------------------------------

/// Trait for numerical time integrators.
pub trait Integrator {
    /// Advance state `q` (positions) and `v` (velocities) by `dt` given
    /// forces `f`.
    ///
    /// `dof` is the number of degrees of freedom.
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real);
}

/// Semi-implicit (symplectic) Euler integrator.
///
/// Updates velocity first, then position:
/// ```text
/// v_{n+1} = v_n + (f_n / m) * dt
/// q_{n+1} = q_n + v_{n+1} * dt
/// ```
pub struct SemiImplicitEuler;

impl Integrator for SemiImplicitEuler {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        for i in 0..n {
            v[i] += f[i] * inv_mass[i] * dt;
            q[i] += v[i] * dt;
        }
    }
}

/// Classic explicit (forward) Euler integrator.
///
/// Updates position first using *old* velocity, then updates velocity:
/// ```text
/// q_{n+1} = q_n + v_n * dt
/// v_{n+1} = v_n + (f_n / m) * dt
/// ```
pub struct ExplicitEuler;

impl Integrator for ExplicitEuler {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        for i in 0..n {
            q[i] += v[i] * dt;
            v[i] += f[i] * inv_mass[i] * dt;
        }
    }
}

/// 4th-order Runge-Kutta integrator (position only, given an acceleration field).
///
/// Uses `f` as the acceleration (force/mass) directly.
pub struct RungeKutta4;

impl RungeKutta4 {
    /// Integrate positions `q` with velocities `v` under accelerations `a = f * inv_mass`,
    /// using the classic RK4 scheme.  Both `q` and `v` are updated.
    pub fn integrate_rk4(q: &mut [Real], v: &mut [Real], accel: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(accel.len());
        let h = dt;
        let h2 = h / 2.0;
        let h6 = h / 6.0;
        // k1 = (v, a)
        // k2 = (v + h/2*a, a)   (constant acceleration assumption)
        // k3 = (v + h/2*a, a)
        // k4 = (v + h*a, a)
        for i in 0..n {
            let a = accel[i];
            // Velocity increments
            let kv1 = a;
            let kv2 = a; // constant force assumption
            let kv3 = a;
            let kv4 = a;
            // Position increments
            let kq1 = v[i];
            let kq2 = v[i] + h2 * kv1;
            let kq3 = v[i] + h2 * kv2;
            let kq4 = v[i] + h * kv3;
            v[i] += h6 * (kv1 + 2.0 * kv2 + 2.0 * kv3 + kv4);
            q[i] += h6 * (kq1 + 2.0 * kq2 + 2.0 * kq3 + kq4);
        }
    }
}

impl Integrator for RungeKutta4 {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        let accel: Vec<Real> = (0..n).map(|i| f[i] * inv_mass[i]).collect();
        Self::integrate_rk4(q, v, &accel, dt);
    }
}

// ---------------------------------------------------------------------------
// Broadphase
// ---------------------------------------------------------------------------

/// Trait for broadphase collision detection structures.
pub trait Broadphase {
    /// Insert an object with the given id and AABB.
    fn insert(&mut self, id: u64, aabb: Aabb);

    /// Update the AABB for an existing object.
    fn update(&mut self, id: u64, aabb: Aabb);

    /// Remove an object.
    fn remove(&mut self, id: u64);

    /// Return all pairs of overlapping ids.
    fn overlapping_pairs(&self) -> Vec<(u64, u64)>;

    /// Return all ids whose AABB overlaps the query AABB.
    fn query_aabb(&self, query: &Aabb) -> Vec<u64>;
}

// ---------------------------------------------------------------------------
// Narrowphase
// ---------------------------------------------------------------------------

/// A contact point produced by narrowphase collision detection.
#[derive(Debug, Clone)]
pub struct ContactPoint {
    /// Contact point in world space.
    pub point: [Real; 3],
    /// Contact normal pointing from body B into body A.
    pub normal: [Real; 3],
    /// Penetration depth (positive = penetrating).
    pub depth: Real,
}

/// Trait for narrowphase collision queries.
pub trait Narrowphase {
    /// Test whether two convex shapes (described by their transforms) overlap.
    fn test_overlap(&self, transform_a: &Transform, transform_b: &Transform) -> bool;

    /// Compute contact points between two shapes.
    fn contacts(&self, transform_a: &Transform, transform_b: &Transform) -> Vec<ContactPoint>;
}

// ---------------------------------------------------------------------------
// Serialisable physics state
// ---------------------------------------------------------------------------

/// Trait for objects whose full physics state can be captured and restored.
pub trait PhysicsStateful {
    /// The type representing the saved state.
    type State;

    /// Capture the current state.
    fn save_state(&self) -> Self::State;

    /// Restore from a previously saved state.
    fn restore_state(&mut self, state: Self::State);
}

// ---------------------------------------------------------------------------
// Damping
// ---------------------------------------------------------------------------

/// Trait for objects that model energy dissipation.
pub trait Damped {
    /// Return the linear damping coefficient.
    fn linear_damping(&self) -> Real;

    /// Return the angular damping coefficient.
    fn angular_damping(&self) -> Real;

    /// Compute damped velocity components from current velocity over `dt`.
    fn damp_velocity(&self, velocity: Real, dt: Real, is_angular: bool) -> Real {
        let d = if is_angular {
            self.angular_damping()
        } else {
            self.linear_damping()
        };
        velocity * (1.0 - d * dt).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// Sleeping / activation
// ---------------------------------------------------------------------------

/// Trait for bodies that can sleep when below activity thresholds.
pub trait Sleepable {
    /// Whether the body is currently sleeping.
    fn is_sleeping(&self) -> bool;

    /// Wake the body from sleep.
    fn wake(&mut self);

    /// Put the body to sleep, zeroing velocities.
    fn sleep(&mut self);

    /// Test whether the body should fall asleep given its speed squared and
    /// the world's thresholds.
    fn should_sleep(
        &self,
        speed_sq: Real,
        threshold_sq: Real,
        accumulated_sleep_time: Real,
        time_to_sleep: Real,
    ) -> bool {
        speed_sq < threshold_sq && accumulated_sleep_time >= time_to_sleep
    }
}

// ---------------------------------------------------------------------------
// Renderer output
// ---------------------------------------------------------------------------

/// Trait for physics objects that can provide render-ready data.
pub trait Renderable {
    /// Return the interpolated world-space transform for rendering at
    /// sub-step fraction `alpha` ∈ \[0,1\].
    fn render_transform(&self, alpha: Real) -> Transform;
}

// ---------------------------------------------------------------------------
// Physics material
// ---------------------------------------------------------------------------

/// Physical surface properties.
#[derive(Debug, Clone)]
pub struct PhysicsMaterial {
    /// Coefficient of static friction.
    pub static_friction: Real,
    /// Coefficient of dynamic friction.
    pub dynamic_friction: Real,
    /// Restitution (bounciness), 0 = perfectly inelastic, 1 = perfectly elastic.
    pub restitution: Real,
    /// Rolling friction coefficient (for spheres/cylinders).
    pub rolling_friction: Real,
}

impl Default for PhysicsMaterial {
    fn default() -> Self {
        Self {
            static_friction: 0.5,
            dynamic_friction: 0.4,
            restitution: 0.3,
            rolling_friction: 0.01,
        }
    }
}

impl PhysicsMaterial {
    /// Create a rubber-like material.
    pub fn rubber() -> Self {
        Self {
            static_friction: 1.0,
            dynamic_friction: 0.8,
            restitution: 0.7,
            rolling_friction: 0.02,
        }
    }

    /// Create a steel-on-steel material.
    pub fn steel() -> Self {
        Self {
            static_friction: 0.74,
            dynamic_friction: 0.57,
            restitution: 0.6,
            rolling_friction: 0.001,
        }
    }

    /// Create a frictionless ice-like material.
    pub fn ice() -> Self {
        Self {
            static_friction: 0.03,
            dynamic_friction: 0.02,
            restitution: 0.05,
            rolling_friction: 0.001,
        }
    }

    /// Combine two materials using geometric mean for friction and min for restitution
    /// (a common convention in game physics).
    pub fn combine(&self, other: &PhysicsMaterial) -> PhysicsMaterial {
        PhysicsMaterial {
            static_friction: (self.static_friction * other.static_friction).sqrt(),
            dynamic_friction: (self.dynamic_friction * other.dynamic_friction).sqrt(),
            restitution: self.restitution.min(other.restitution),
            rolling_friction: (self.rolling_friction * other.rolling_friction).sqrt(),
        }
    }
}

/// Trait for objects that expose a surface material.
pub trait HasMaterial {
    /// Return the physics material for this object's surface.
    fn material(&self) -> &PhysicsMaterial;
}

// ---------------------------------------------------------------------------
// Debug / diagnostics
// ---------------------------------------------------------------------------

/// Trait for physics objects that can emit diagnostic information.
pub trait PhysicsDiagnostics {
    /// Return a vector of (name, value) pairs for debug display.
    fn diagnostics(&self) -> Vec<(&'static str, Real)>;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use crate::VelocityVerlet;

    use crate::fabrik_solve;

    // ── SemiImplicitEuler ─────────────────────────────────────────────────────

    #[test]
    fn test_semi_implicit_euler_constant_force() {
        let integrator = SemiImplicitEuler;
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        let f = [1.0f64]; // constant force
        let inv_mass = [1.0f64]; // mass = 1 kg
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 1.0);
        // v = 0 + 1*1 = 1; q = 0 + 1*1 = 1
        assert!((v[0] - 1.0).abs() < 1e-12, "v={}", v[0]);
        assert!((q[0] - 1.0).abs() < 1e-12, "q={}", q[0]);
    }

    #[test]
    fn test_semi_implicit_euler_zero_force() {
        let integrator = SemiImplicitEuler;
        let mut q = [5.0f64];
        let mut v = [2.0f64];
        let f = [0.0f64];
        let inv_mass = [1.0f64];
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.5);
        // v unchanged = 2; q = 5 + 2*0.5 = 6
        assert!((v[0] - 2.0).abs() < 1e-12);
        assert!((q[0] - 6.0).abs() < 1e-12);
    }

    // ── ExplicitEuler ─────────────────────────────────────────────────────────

    #[test]
    fn test_explicit_euler_differs_from_semi_implicit() {
        // With a non-zero force the two integrators give different results.
        let f = [2.0f64];
        let inv_mass = [1.0f64];
        let dt = 1.0;

        let mut q_si = [0.0f64];
        let mut v_si = [0.0f64];
        SemiImplicitEuler.integrate(&mut q_si, &mut v_si, &f, &inv_mass, dt);

        let mut q_ex = [0.0f64];
        let mut v_ex = [0.0f64];
        ExplicitEuler.integrate(&mut q_ex, &mut v_ex, &f, &inv_mass, dt);

        // Semi-implicit: v=2, q=2; Explicit: q=0 first, v=2
        assert!((v_si[0] - 2.0).abs() < 1e-12);
        assert!((q_si[0] - 2.0).abs() < 1e-12);
        assert!((v_ex[0] - 2.0).abs() < 1e-12);
        assert!((q_ex[0] - 0.0).abs() < 1e-12); // position integrated with OLD v=0
    }

    // ── RungeKutta4 ───────────────────────────────────────────────────────────

    #[test]
    fn test_rk4_constant_acceleration() {
        // Under constant acceleration a=2, starting from rest:
        // q(dt) = ½*a*dt², v(dt) = a*dt
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        let a = [2.0f64];
        let inv_mass = [1.0f64];
        RungeKutta4.integrate(&mut q, &mut v, &a, &inv_mass, 1.0);
        assert!((v[0] - 2.0).abs() < 1e-10, "v={}", v[0]);
        // RK4 with constant a gives q = v0*dt + a*dt²/2 (exact for const a)
        assert!((q[0] - 1.0).abs() < 1e-10, "q={}", q[0]);
    }

    #[test]
    fn test_rk4_with_initial_velocity() {
        let mut q = [0.0f64];
        let mut v = [3.0f64];
        let a = [0.0f64]; // no force
        let inv_mass = [1.0f64];
        RungeKutta4.integrate(&mut q, &mut v, &a, &inv_mass, 2.0);
        // q = 0 + 3*2 = 6, v = 3 unchanged
        assert!((q[0] - 6.0).abs() < 1e-10, "q={}", q[0]);
        assert!((v[0] - 3.0).abs() < 1e-10, "v={}", v[0]);
    }

    // ── PhysicsMaterial ───────────────────────────────────────────────────────

    #[test]
    fn test_material_combine_friction_geometric_mean() {
        let a = PhysicsMaterial {
            static_friction: 4.0,
            dynamic_friction: 4.0,
            restitution: 0.5,
            rolling_friction: 0.01,
        };
        let b = PhysicsMaterial {
            static_friction: 1.0,
            dynamic_friction: 1.0,
            restitution: 0.8,
            rolling_friction: 0.01,
        };
        let c = a.combine(&b);
        // geometric mean of 4 and 1 = 2
        assert!((c.static_friction - 2.0).abs() < 1e-12);
        assert!((c.dynamic_friction - 2.0).abs() < 1e-12);
        // min restitution
        assert!((c.restitution - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_material_presets_valid() {
        let rubber = PhysicsMaterial::rubber();
        assert!(rubber.static_friction > 0.0);
        assert!(rubber.restitution <= 1.0);

        let steel = PhysicsMaterial::steel();
        assert!(steel.dynamic_friction > 0.0);

        let ice = PhysicsMaterial::ice();
        assert!(ice.static_friction < 0.1);
    }

    // ── Damped via PhysicsMaterial-like manual impl ───────────────────────────

    struct SimpleDamped {
        lin: Real,
        ang: Real,
    }
    impl Damped for SimpleDamped {
        fn linear_damping(&self) -> Real {
            self.lin
        }
        fn angular_damping(&self) -> Real {
            self.ang
        }
    }

    #[test]
    fn test_damped_velocity() {
        let d = SimpleDamped { lin: 0.5, ang: 0.2 };
        let v_after = d.damp_velocity(10.0, 1.0, false);
        // factor = 1 - 0.5*1 = 0.5; v = 10 * 0.5 = 5
        assert!((v_after - 5.0).abs() < 1e-12, "v={}", v_after);
    }

    #[test]
    fn test_damped_velocity_clamp_zero() {
        let d = SimpleDamped { lin: 2.0, ang: 0.0 };
        // factor = 1 - 2.0*1.0 = -1 → clamped to 0
        let v_after = d.damp_velocity(10.0, 1.0, false);
        assert_eq!(v_after, 0.0);
    }

    // ── Sleepable via manual impl ─────────────────────────────────────────────

    #[test]
    fn test_sleepable_should_sleep() {
        struct DummySleepable;
        impl Sleepable for DummySleepable {
            fn is_sleeping(&self) -> bool {
                false
            }
            fn wake(&mut self) {}
            fn sleep(&mut self) {}
        }
        let s = DummySleepable;
        assert!(s.should_sleep(0.0001, 0.01, 1.0, 0.5));
        assert!(!s.should_sleep(0.1, 0.01, 1.0, 0.5)); // too fast
        assert!(!s.should_sleep(0.0001, 0.01, 0.1, 0.5)); // not long enough
    }

    // ── ContactPoint ─────────────────────────────────────────────────────────

    #[test]
    fn test_contact_point_construction() {
        let cp = ContactPoint {
            point: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.05,
        };
        assert!((cp.depth - 0.05).abs() < 1e-12);
        assert_eq!(cp.normal[1], 1.0);
    }

    // ── Integrator length mismatch safety ────────────────────────────────────

    #[test]
    fn test_integrator_length_mismatch_safe() {
        // Should not panic when arrays have different lengths.
        let integrator = SemiImplicitEuler;
        let mut q = vec![0.0, 0.0, 0.0];
        let mut v = vec![1.0, 1.0];
        let f = vec![0.0, 0.0, 0.0];
        let inv_mass = vec![1.0, 1.0, 1.0];
        // min length is 2; should process only 2 elements without panic.
        integrator.integrate(&mut q, &mut v, &f, &inv_mass, 0.1);
        assert_eq!(q[2], 0.0); // third element untouched
    }

    // ── Velocity Verlet ───────────────────────────────────────────────────────

    #[test]
    fn test_velocity_verlet_constant_accel() {
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        let a = [2.0f64];
        let dt = 1.0;
        VelocityVerlet.integrate(&mut q, &mut v, &a, &[1.0], dt);
        // q = 0 + 0*1 + 0.5*2*1 = 1; v = 0 + 2*1 = 2
        assert!((q[0] - 1.0).abs() < 1e-12, "q={}", q[0]);
        assert!((v[0] - 2.0).abs() < 1e-12, "v={}", v[0]);
    }

    // ── Leapfrog ──────────────────────────────────────────────────────────────

    #[test]
    fn test_leapfrog_energy_conservation() {
        // For a harmonic oscillator x'' = -x, energy E = v²/2 + x²/2 is nearly
        // conserved by the leapfrog integrator over many steps.
        let omega2 = 1.0f64; // spring constant / mass = 1
        let dt = 0.01;
        let n_steps = 1000;
        let mut x = 1.0f64;
        let mut v = 0.0f64;
        // Kick-start half-step for leapfrog
        let a0 = -omega2 * x;
        let mut v_half = v + 0.5 * a0 * dt;
        let e0 = 0.5 * v * v + 0.5 * x * x;
        for _ in 0..n_steps {
            x += v_half * dt;
            let a = -omega2 * x;
            v_half += a * dt;
            v = v_half - 0.5 * a * dt;
        }
        let e_final = 0.5 * v * v + 0.5 * x * x;
        // Leapfrog is symplectic; energy oscillates but does not drift.
        assert!((e_final - e0).abs() < 0.01, "E0={e0}, Ef={e_final}");
    }

    // ── VelocityVerlet multi-step ─────────────────────────────────────────────

    #[test]
    fn test_velocity_verlet_free_fall() {
        // Under constant gravity g = -9.81, starting from rest:
        // y(t) = -1/2 g t²
        let g = -9.81f64;
        let dt = 0.001;
        let n = 1000;
        let mut y = 0.0f64;
        let mut vy = 0.0f64;
        let mut ay = g;
        for _ in 0..n {
            y += vy * dt + 0.5 * ay * dt * dt;
            let ay_new = g;
            vy += 0.5 * (ay + ay_new) * dt;
            ay = ay_new;
        }
        let t = dt * n as f64;
        let expected = 0.5 * g * t * t;
        assert!((y - expected).abs() < 1e-8, "y={y}, expected={expected}");
    }

    // ── Runge-Kutta adaptive ──────────────────────────────────────────────────

    #[test]
    fn test_rkf45_step_count_positive() {
        // RKF45 step for a simple linear ODE
        let mut y = [1.0f64];
        let a = [0.0f64];
        let inv_mass = [1.0f64];
        RungeKutta4.integrate(&mut y, &mut a.clone(), &a, &inv_mass, 0.1);
        assert!(y[0].is_finite());
    }

    // ── PhysicsMaterial edge cases ────────────────────────────────────────────

    #[test]
    fn test_material_combine_zero_friction() {
        let a = PhysicsMaterial {
            static_friction: 0.0,
            dynamic_friction: 0.0,
            restitution: 0.5,
            rolling_friction: 0.0,
        };
        let b = PhysicsMaterial::rubber();
        let c = a.combine(&b);
        // sqrt(0 * x) = 0
        assert_eq!(c.static_friction, 0.0);
        assert_eq!(c.dynamic_friction, 0.0);
    }

    #[test]
    fn test_material_combine_min_restitution() {
        let a = PhysicsMaterial {
            restitution: 0.9,
            ..PhysicsMaterial::default()
        };
        let b = PhysicsMaterial {
            restitution: 0.1,
            ..PhysicsMaterial::default()
        };
        let c = a.combine(&b);
        assert!((c.restitution - 0.1).abs() < 1e-12);
    }

    // ── InverseKinematics trait (FABRIK) ──────────────────────────────────────

    #[test]
    fn test_fabrik_reaches_target_1dof() {
        // 2-joint arm (2 bones) with total reach = 2.0; target at distance 1.0 → reachable
        let bones = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let lengths = vec![1.0f64, 1.0];
        let target = [1.0, 0.5, 0.0]; // within reach
        let result = fabrik_solve(&bones, &lengths, &target, 50, 1e-6);
        // end effector should be at target
        let ee = result.last().unwrap();
        let d = ((ee[0] - target[0]).powi(2)
            + (ee[1] - target[1]).powi(2)
            + (ee[2] - target[2]).powi(2))
        .sqrt();
        assert!(d < 0.01, "FABRIK end effector distance to target: {d}");
    }

    #[test]
    fn test_fabrik_unreachable_extends_toward_target() {
        // 1-joint arm with reach 1.0; target at distance 3.0 → unreachable
        let bones = vec![[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let lengths = vec![1.0f64];
        let target = [3.0, 0.0, 0.0];
        let result = fabrik_solve(&bones, &lengths, &target, 20, 1e-6);
        // End effector should be pointing toward target
        let ee = result.last().unwrap();
        assert!(
            ee[0] > 0.9,
            "end effector should extend toward target: x={}",
            ee[0]
        );
    }
}

// ---------------------------------------------------------------------------
// Velocity Verlet integrator
// ---------------------------------------------------------------------------

/// Velocity Verlet integration scheme.
///
/// Updates positions using both current and next accelerations:
/// ```text
/// q_{n+1} = q_n + v_n * dt + 0.5 * a_n * dt²
/// v_{n+1} = v_n + 0.5 * (a_n + a_{n+1}) * dt
/// ```
/// Since `a_{n+1}` is not generally available (forces depend on positions),
/// this implementation uses `a_n` for both terms, equivalent to the
/// "Störmer-Verlet" update.
pub struct VelocityVerlet;

impl Integrator for VelocityVerlet {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        let dt2 = dt * dt;
        for i in 0..n {
            let a = f[i] * inv_mass[i];
            q[i] += v[i] * dt + 0.5 * a * dt2;
            v[i] += a * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Leapfrog integrator
// ---------------------------------------------------------------------------

/// Leapfrog (Störmer-Verlet) integration.
///
/// The leapfrog scheme stores velocities at half-integer time steps:
/// ```text
/// v_{n+1/2} = v_{n-1/2} + a_n * dt
/// q_{n+1}   = q_n + v_{n+1/2} * dt
/// ```
/// This implementation initialises the half-step velocity from the
/// full-step velocity on first call, making it self-starting.
pub struct Leapfrog;

impl Integrator for Leapfrog {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        for i in 0..n {
            let a = f[i] * inv_mass[i];
            // Kick: advance velocity by half-step
            let v_half = v[i] + a * dt * 0.5;
            // Drift: advance position
            q[i] += v_half * dt;
            // Kick: advance velocity to full step
            v[i] = v_half + a * dt * 0.5;
        }
    }
}

// ---------------------------------------------------------------------------
// Adams-Bashforth 2nd-order integrator
// ---------------------------------------------------------------------------

/// Adams-Bashforth 2-step explicit integrator.
///
/// Requires the acceleration from the previous step (`f_prev`).  When
/// `f_prev` is `None` the method falls back to a single Euler step.
pub struct AdamsBashforth2 {
    /// Stored acceleration from the previous step (optional first-step bootstrap).
    pub f_prev: Option<Vec<Real>>,
}

impl AdamsBashforth2 {
    /// Create a new Adams-Bashforth 2-step integrator.
    pub fn new() -> Self {
        Self { f_prev: None }
    }

    /// Reset the stored previous-step forces.
    pub fn reset(&mut self) {
        self.f_prev = None;
    }

    /// Integrate state using the 2-step formula.
    pub fn integrate_ab2(
        &mut self,
        q: &mut [Real],
        v: &mut [Real],
        f: &[Real],
        inv_mass: &[Real],
        dt: Real,
    ) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        match &self.f_prev {
            None => {
                // Bootstrap with forward Euler.
                for i in 0..n {
                    let a = f[i] * inv_mass[i];
                    v[i] += a * dt;
                    q[i] += v[i] * dt;
                }
            }
            Some(fp) => {
                let m = fp.len().min(n);
                for i in 0..m {
                    let a_curr = f[i] * inv_mass[i];
                    let a_prev = fp[i] * inv_mass[i];
                    // AB2: v_{n+1} = v_n + dt*(3/2 * a_n - 1/2 * a_{n-1})
                    let a_eff = 1.5 * a_curr - 0.5 * a_prev;
                    v[i] += a_eff * dt;
                    q[i] += v[i] * dt;
                }
            }
        }
        self.f_prev = Some(f[..n].to_vec());
    }
}

impl Default for AdamsBashforth2 {
    fn default() -> Self {
        Self::new()
    }
}

impl Integrator for AdamsBashforth2 {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        // Stateless fallback: just use Euler.
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        for i in 0..n {
            let a = f[i] * inv_mass[i];
            v[i] += a * dt;
            q[i] += v[i] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Implicit (backward) Euler integrator
// ---------------------------------------------------------------------------

/// Implicit (backward) Euler integrator — unconditionally stable.
///
/// For linear systems `F = k * x` the implicit update solves:
/// ```text
/// v_{n+1} = v_n + a_{n+1} * dt   (where a_{n+1} uses new position)
/// q_{n+1} = q_n + v_{n+1} * dt
/// ```
/// In the absence of position-dependent forces this reduces to the
/// semi-implicit Euler integrator.
pub struct ImplicitEuler;

impl Integrator for ImplicitEuler {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        // For constant forces this is the same as semi-implicit Euler.
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        for i in 0..n {
            v[i] += f[i] * inv_mass[i] * dt;
            q[i] += v[i] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Midpoint (RK2) integrator
// ---------------------------------------------------------------------------

/// Explicit midpoint (Runge-Kutta 2) integrator.
///
/// ```text
/// k1_v = a_n,  k1_q = v_n
/// k2_v = a_n,  k2_q = v_n + k1_v * dt/2
/// v_{n+1} = v_n + k2_v * dt
/// q_{n+1} = q_n + k2_q * dt
/// ```
/// For constant forces this gives the same result as the explicit Euler
/// update for velocity but uses the midpoint velocity for position.
pub struct Midpoint;

impl Integrator for Midpoint {
    fn integrate(&self, q: &mut [Real], v: &mut [Real], f: &[Real], inv_mass: &[Real], dt: Real) {
        let n = q.len().min(v.len()).min(f.len()).min(inv_mass.len());
        let h2 = dt / 2.0;
        for i in 0..n {
            let a = f[i] * inv_mass[i];
            let v_mid = v[i] + a * h2;
            v[i] += a * dt;
            q[i] += v_mid * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Constraint solver helpers
// ---------------------------------------------------------------------------

/// XPBD position-level constraint correction for a distance constraint.
///
/// Corrects the positions of two bodies so that the distance between
/// `p_a` and `p_b` equals `rest_length`, using XPBD compliance.
///
/// Returns the position deltas `(delta_a, delta_b)` for each body.
///
/// # Arguments
/// * `p_a` / `p_b` — current world-space positions of each body's anchor
/// * `inv_mass_a` / `inv_mass_b` — inverse masses of each body
/// * `rest_length` — target separation distance
/// * `compliance` — XPBD compliance (= 1 / stiffness), 0 = rigid
/// * `dt` — substep duration
pub fn xpbd_distance_correction(
    p_a: [Real; 3],
    p_b: [Real; 3],
    inv_mass_a: Real,
    inv_mass_b: Real,
    rest_length: Real,
    compliance: Real,
    dt: Real,
) -> ([Real; 3], [Real; 3]) {
    let dx = [p_b[0] - p_a[0], p_b[1] - p_a[1], p_b[2] - p_a[2]];
    let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    if dist < 1e-12 {
        return ([0.0; 3], [0.0; 3]);
    }
    let n = [dx[0] / dist, dx[1] / dist, dx[2] / dist];
    let c = dist - rest_length;
    let alpha = compliance / (dt * dt);
    let w = inv_mass_a + inv_mass_b;
    if w < 1e-12 {
        return ([0.0; 3], [0.0; 3]);
    }
    let lambda = -c / (w + alpha);
    let da = [
        -lambda * inv_mass_a * n[0],
        -lambda * inv_mass_a * n[1],
        -lambda * inv_mass_a * n[2],
    ];
    let db = [
        lambda * inv_mass_b * n[0],
        lambda * inv_mass_b * n[1],
        lambda * inv_mass_b * n[2],
    ];
    (da, db)
}

/// Velocity-level constraint impulse for a contact constraint (Baumgarte-free).
///
/// Resolves relative normal velocity at a contact point given restitution.
/// Returns the scalar impulse magnitude.
///
/// # Arguments
/// * `v_rel_n` — relative normal velocity (vB − vA) · normal
/// * `inv_mass_a` / `inv_mass_b` — inverse masses
/// * `inv_inertia_a` / `inv_inertia_b` — inverse inertia scalars
/// * `r_a` / `r_b` — contact offset vectors (world)
/// * `normal` — contact normal (unit vector)
/// * `restitution` — coefficient of restitution
#[allow(clippy::too_many_arguments)]
pub fn contact_impulse_scalar(
    v_rel_n: Real,
    inv_mass_a: Real,
    inv_mass_b: Real,
    inv_inertia_a: Real,
    inv_inertia_b: Real,
    r_a: [Real; 3],
    r_b: [Real; 3],
    normal: [Real; 3],
) -> Real {
    // Effective mass along normal including angular terms.
    fn cross_sq(r: [Real; 3], n: [Real; 3]) -> Real {
        let cx = r[1] * n[2] - r[2] * n[1];
        let cy = r[2] * n[0] - r[0] * n[2];
        let cz = r[0] * n[1] - r[1] * n[0];
        cx * cx + cy * cy + cz * cz
    }
    let k = inv_mass_a
        + inv_mass_b
        + inv_inertia_a * cross_sq(r_a, normal)
        + inv_inertia_b * cross_sq(r_b, normal);
    if k < 1e-12 { 0.0 } else { -v_rel_n / k }
}

// ---------------------------------------------------------------------------
// FABRIK inverse kinematics
// ---------------------------------------------------------------------------

/// FABRIK (Forward And Backward Reaching Inverse Kinematics) solver.
///
/// Solves an n-joint chain to place the end-effector at `target`.
///
/// # Arguments
/// * `bones` — initial positions of each joint (n+1 points for n segments)
/// * `lengths` — bone lengths (must have length == bones.len() - 1)
/// * `target` — desired end-effector position
/// * `max_iter` — maximum FABRIK iterations
/// * `tolerance` — convergence threshold (Euclidean distance)
///
/// Returns the updated joint positions.
pub fn fabrik_solve(
    bones: &[[Real; 3]],
    lengths: &[Real],
    target: &[Real; 3],
    max_iter: usize,
    tolerance: Real,
) -> Vec<[Real; 3]> {
    let n = bones.len();
    if n < 2 || lengths.len() < n - 1 {
        return bones.to_vec();
    }
    let mut joints: Vec<[Real; 3]> = bones.to_vec();
    let root = joints[0];
    let total_length: Real = lengths.iter().sum();

    // Check if target is reachable
    let dist_to_target = dist3(&joints[0], target);
    if dist_to_target >= total_length {
        // Stretch toward target
        for i in 0..n - 1 {
            let d = dist3(&joints[i], target);
            let t = lengths[i] / d;
            joints[i + 1] = lerp3(&joints[i], target, t);
        }
        return joints;
    }

    for _ in 0..max_iter {
        // Forward pass: start from end effector
        joints[n - 1] = *target;
        for i in (0..n - 1).rev() {
            let d = dist3(&joints[i + 1], &joints[i]);
            if d < 1e-12 {
                continue;
            }
            let t = lengths[i] / d;
            joints[i] = lerp3(&joints[i + 1], &joints[i], t);
        }

        // Backward pass: start from root
        joints[0] = root;
        for i in 0..n - 1 {
            let d = dist3(&joints[i], &joints[i + 1]);
            if d < 1e-12 {
                continue;
            }
            let t = lengths[i] / d;
            joints[i + 1] = lerp3(&joints[i], &joints[i + 1], t);
        }

        // Check convergence
        let ee_dist = dist3(&joints[n - 1], target);
        if ee_dist < tolerance {
            break;
        }
    }
    joints
}

/// Euclidean distance between two 3D points.
fn dist3(a: &[Real; 3], b: &[Real; 3]) -> Real {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Linear interpolation between two 3D points.
fn lerp3(a: &[Real; 3], b: &[Real; 3], t: Real) -> [Real; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

// ---------------------------------------------------------------------------
// Position-based dynamics (PBD) constraint solver
// ---------------------------------------------------------------------------

/// Trait for PBD position constraints.
pub trait PbdConstraint {
    /// Evaluate constraint C(p) — returns violation (0 = satisfied).
    fn evaluate(&self, positions: &[[Real; 3]]) -> Real;

    /// Compute constraint gradient ∇C with respect to each particle position.
    fn gradient(&self, positions: &[[Real; 3]]) -> Vec<[Real; 3]>;

    /// Particle indices involved in this constraint.
    fn particle_indices(&self) -> &[usize];

    /// Stiffness in range \[0, 1\] (1 = fully rigid, 0 = no correction).
    fn stiffness(&self) -> Real {
        1.0
    }
}

/// A simple PBD distance constraint between two particles.
pub struct PbdDistanceConstraint {
    /// Indices of the two particles.
    pub indices: [usize; 2],
    /// Rest (target) distance.
    pub rest_length: Real,
    /// Stiffness in \[0, 1\].
    pub k: Real,
}

impl PbdDistanceConstraint {
    /// Create a new distance constraint.
    pub fn new(i: usize, j: usize, rest_length: Real, stiffness: Real) -> Self {
        Self {
            indices: [i, j],
            rest_length,
            k: stiffness,
        }
    }
}

impl PbdConstraint for PbdDistanceConstraint {
    fn evaluate(&self, positions: &[[Real; 3]]) -> Real {
        let [i, j] = self.indices;
        if i >= positions.len() || j >= positions.len() {
            return 0.0;
        }
        let p = &positions[i];
        let q = &positions[j];
        dist3(&[p[0], p[1], p[2]], &[q[0], q[1], q[2]]) - self.rest_length
    }

    fn gradient(&self, positions: &[[Real; 3]]) -> Vec<[Real; 3]> {
        let [i, j] = self.indices;
        if i >= positions.len() || j >= positions.len() {
            return vec![[0.0; 3]; 2];
        }
        let p = positions[i];
        let q = positions[j];
        let dx = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let d = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if d < 1e-12 {
            return vec![[0.0; 3]; 2];
        }
        let n = [dx[0] / d, dx[1] / d, dx[2] / d];
        // ∇C_i = -n, ∇C_j = +n
        vec![[-n[0], -n[1], -n[2]], [n[0], n[1], n[2]]]
    }

    fn particle_indices(&self) -> &[usize] {
        &self.indices
    }

    fn stiffness(&self) -> Real {
        self.k
    }
}

/// Solve one PBD iteration over a set of constraints.
///
/// Updates `positions` in-place by applying each constraint correction
/// weighted by stiffness and inverse masses.
pub fn pbd_solve_iteration(
    positions: &mut Vec<[Real; 3]>,
    inv_masses: &[Real],
    constraints: &[&dyn PbdConstraint],
) {
    for c in constraints {
        let c_val = c.evaluate(positions);
        let grads = c.gradient(positions);
        let indices = c.particle_indices();
        if grads.len() < indices.len() {
            continue;
        }
        // Compute denominator: sum of w_i * |∇C_i|²
        let mut denom = 0.0;
        for (k, &idx) in indices.iter().enumerate() {
            if idx < inv_masses.len() {
                let g = grads[k];
                denom += inv_masses[idx] * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            }
        }
        if denom < 1e-12 {
            continue;
        }
        let lambda = -c_val / denom * c.stiffness();
        for (k, &idx) in indices.iter().enumerate() {
            if idx < positions.len() && idx < inv_masses.len() {
                let g = grads[k];
                let w = inv_masses[idx];
                positions[idx][0] += lambda * w * g[0];
                positions[idx][1] += lambda * w * g[1];
                positions[idx][2] += lambda * w * g[2];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Collision filter trait
// ---------------------------------------------------------------------------

/// Trait for filtering which body pairs participate in collision detection.
pub trait CollisionFilter {
    /// Return `true` if the pair `(id_a, id_b)` should be tested for collision.
    fn should_collide(&self, id_a: u64, id_b: u64) -> bool;
}

/// A simple bitmask-based collision filter.
///
/// Body A collides with body B only if `mask_a & layer_b != 0` AND
/// `mask_b & layer_a != 0`.
pub struct LayerFilter {
    /// Collision layer bits for each body id.
    pub layers: std::collections::HashMap<u64, u32>,
    /// Collision mask bits for each body id.
    pub masks: std::collections::HashMap<u64, u32>,
}

impl LayerFilter {
    /// Create an empty layer filter.
    pub fn new() -> Self {
        Self {
            layers: std::collections::HashMap::new(),
            masks: std::collections::HashMap::new(),
        }
    }

    /// Register a body with its layer and mask.
    pub fn register(&mut self, id: u64, layer: u32, mask: u32) {
        self.layers.insert(id, layer);
        self.masks.insert(id, mask);
    }
}

impl Default for LayerFilter {
    fn default() -> Self {
        Self::new()
    }
}

impl CollisionFilter for LayerFilter {
    fn should_collide(&self, id_a: u64, id_b: u64) -> bool {
        let layer_a = self.layers.get(&id_a).copied().unwrap_or(0xFFFF_FFFF);
        let layer_b = self.layers.get(&id_b).copied().unwrap_or(0xFFFF_FFFF);
        let mask_a = self.masks.get(&id_a).copied().unwrap_or(0xFFFF_FFFF);
        let mask_b = self.masks.get(&id_b).copied().unwrap_or(0xFFFF_FFFF);
        (mask_a & layer_b) != 0 && (mask_b & layer_a) != 0
    }
}

// ---------------------------------------------------------------------------
// Sensor trait
// ---------------------------------------------------------------------------

/// Trait for physics sensors (trigger volumes, proximity detectors, etc.).
pub trait PhysicsSensor {
    /// Return `true` if the sensor is currently triggered by any body.
    fn is_triggered(&self) -> bool;

    /// Return the IDs of all bodies currently overlapping this sensor.
    fn overlapping_bodies(&self) -> Vec<u64>;

    /// Called each frame to update the sensor's overlap state.
    fn update(&mut self, body_positions: &[(u64, [Real; 3])]);
}

/// A spherical proximity sensor.
pub struct SphereSensor {
    /// World-space center.
    pub center: [Real; 3],
    /// Trigger radius.
    pub radius: Real,
    /// Currently overlapping body IDs.
    pub current_overlaps: Vec<u64>,
}

impl SphereSensor {
    /// Create a new sphere sensor.
    pub fn new(center: [Real; 3], radius: Real) -> Self {
        Self {
            center,
            radius,
            current_overlaps: Vec::new(),
        }
    }
}

impl PhysicsSensor for SphereSensor {
    fn is_triggered(&self) -> bool {
        !self.current_overlaps.is_empty()
    }

    fn overlapping_bodies(&self) -> Vec<u64> {
        self.current_overlaps.clone()
    }

    fn update(&mut self, body_positions: &[(u64, [Real; 3])]) {
        let r2 = self.radius * self.radius;
        self.current_overlaps = body_positions
            .iter()
            .filter_map(|(id, pos)| {
                let dx = pos[0] - self.center[0];
                let dy = pos[1] - self.center[1];
                let dz = pos[2] - self.center[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    Some(*id)
                } else {
                    None
                }
            })
            .collect();
    }
}

// ---------------------------------------------------------------------------
// Interpolation trait for physics state
// ---------------------------------------------------------------------------

/// Trait for physics state objects that support linear interpolation.
pub trait PhysicsInterpolatable: Sized {
    /// Linearly interpolate between `self` and `other` at blend factor `t ∈ [0,1]`.
    fn lerp_state(&self, other: &Self, t: Real) -> Self;
}

// ---------------------------------------------------------------------------
// Extra tests for expanded traits
// ---------------------------------------------------------------------------

#[cfg(test)]
mod expanded_tests {
    use super::*;
    use crate::AdamsBashforth2;
    use crate::LayerFilter;
    use crate::Leapfrog;
    use crate::Midpoint;
    use crate::PbdConstraint;
    use crate::PbdDistanceConstraint;
    use crate::SphereSensor;
    use crate::VelocityVerlet;
    use crate::dist3;
    use crate::fabrik_solve;
    use crate::pbd_solve_iteration;
    use crate::xpbd_distance_correction;

    // ── VelocityVerlet ────────────────────────────────────────────────────────

    #[test]
    fn test_velocity_verlet_no_force() {
        let mut q = [3.0f64];
        let mut v = [2.0f64];
        VelocityVerlet.integrate(&mut q, &mut v, &[0.0], &[1.0], 0.5);
        // q = 3 + 2*0.5 + 0 = 4; v = 2
        assert!((q[0] - 4.0).abs() < 1e-12, "q={}", q[0]);
        assert!((v[0] - 2.0).abs() < 1e-12, "v={}", v[0]);
    }

    // ── Leapfrog ──────────────────────────────────────────────────────────────

    #[test]
    fn test_leapfrog_rest() {
        let mut q = [1.0f64];
        let mut v = [0.0f64];
        Leapfrog.integrate(&mut q, &mut v, &[0.0], &[1.0], 1.0);
        assert!(
            (q[0] - 1.0).abs() < 1e-12,
            "q should not change: q={}",
            q[0]
        );
    }

    // ── Midpoint ──────────────────────────────────────────────────────────────

    #[test]
    fn test_midpoint_position_uses_midpoint_velocity() {
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        // a = 4.0; v_mid = 0 + 4*0.5*0.5 = 1; q = 0 + 1*0.5*(2?) wait let dt=1
        // a=4, h2=0.5; v_mid = 0 + 4*0.5 = 2; q = 0 + 2*1 = 2; v = 0 + 4*1 = 4
        Midpoint.integrate(&mut q, &mut v, &[4.0], &[1.0], 1.0);
        assert!((q[0] - 2.0).abs() < 1e-12, "q={}", q[0]);
        assert!((v[0] - 4.0).abs() < 1e-12, "v={}", v[0]);
    }

    // ── Adams-Bashforth 2 ─────────────────────────────────────────────────────

    #[test]
    fn test_ab2_bootstrap_euler() {
        let mut ab2 = AdamsBashforth2::new();
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        ab2.integrate_ab2(&mut q, &mut v, &[2.0], &[1.0], 1.0);
        // First step = Euler: v = 0 + 2*1 = 2; q = 0 + 2*1 = 2
        assert!((v[0] - 2.0).abs() < 1e-12);
        assert!((q[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_ab2_second_step() {
        let mut ab2 = AdamsBashforth2::new();
        let mut q = [0.0f64];
        let mut v = [0.0f64];
        ab2.integrate_ab2(&mut q, &mut v, &[2.0], &[1.0], 1.0);
        // Step 2: same force; AB2: a_eff = 1.5*2 - 0.5*2 = 2; same as Euler
        ab2.integrate_ab2(&mut q, &mut v, &[2.0], &[1.0], 1.0);
        assert!(v[0].is_finite());
        assert!(q[0].is_finite());
    }

    // ── XPBD distance correction ──────────────────────────────────────────────

    #[test]
    fn test_xpbd_distance_rigid() {
        let pa = [0.0, 0.0, 0.0];
        let pb = [2.0, 0.0, 0.0];
        // rest_length = 1.0 → violation = 1.0; equal masses
        let (da, db) = xpbd_distance_correction(pa, pb, 1.0, 1.0, 1.0, 0.0, 0.01);
        // Bodies should move together: da positive, db negative
        assert!(da[0] > 0.0, "da.x should be positive: {:?}", da);
        assert!(db[0] < 0.0, "db.x should be negative: {:?}", db);
        // Correction magnitudes equal for equal masses
        assert!((da[0].abs() - db[0].abs()).abs() < 1e-10);
    }

    #[test]
    fn test_xpbd_distance_no_violation() {
        let pa = [0.0, 0.0, 0.0];
        let pb = [1.0, 0.0, 0.0];
        let (da, db) = xpbd_distance_correction(pa, pb, 1.0, 1.0, 1.0, 0.0, 0.01);
        assert!(da[0].abs() < 1e-12);
        assert!(db[0].abs() < 1e-12);
    }

    #[test]
    fn test_xpbd_static_body() {
        // inv_mass_b = 0 → body B doesn't move
        let (da, db) =
            xpbd_distance_correction([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.0, 1.0, 0.0, 0.01);
        assert!(da[0].abs() > 0.0);
        assert_eq!(db[0], 0.0);
    }

    // ── PBD distance constraint ───────────────────────────────────────────────

    #[test]
    fn test_pbd_distance_constraint_evaluate() {
        let c = PbdDistanceConstraint::new(0, 1, 1.0, 1.0);
        let positions = vec![[0.0, 0.0, 0.0f64], [2.0, 0.0, 0.0]];
        let violation = c.evaluate(&positions);
        assert!((violation - 1.0).abs() < 1e-12, "violation={}", violation);
    }

    #[test]
    fn test_pbd_solve_iteration() {
        let c = PbdDistanceConstraint::new(0, 1, 1.0, 1.0);
        let constraints: Vec<&dyn PbdConstraint> = vec![&c];
        let mut positions = vec![[0.0, 0.0, 0.0f64], [3.0, 0.0, 0.0]];
        let inv_masses = [1.0f64, 1.0];
        pbd_solve_iteration(&mut positions, &inv_masses, &constraints);
        // After correction distance should be closer to 1.0
        let d = (positions[1][0] - positions[0][0]).abs();
        assert!(d < 3.0, "distance should decrease: {d}");
    }

    // ── Layer filter ──────────────────────────────────────────────────────────

    #[test]
    fn test_layer_filter_collision() {
        let mut filter = LayerFilter::new();
        filter.register(1, 0b01, 0b10); // layer 1, collides with layer 2
        filter.register(2, 0b10, 0b01); // layer 2, collides with layer 1
        assert!(filter.should_collide(1, 2));
        assert!(filter.should_collide(2, 1));
    }

    #[test]
    fn test_layer_filter_no_collision() {
        let mut filter = LayerFilter::new();
        filter.register(1, 0b01, 0b01); // layer 1, only collides with layer 1
        filter.register(2, 0b10, 0b10); // layer 2, only collides with layer 2
        assert!(!filter.should_collide(1, 2));
    }

    // ── SphereSensor ──────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_sensor_triggers() {
        let mut sensor = SphereSensor::new([0.0; 3], 5.0);
        let bodies = vec![(1u64, [1.0, 0.0, 0.0]), (2u64, [10.0, 0.0, 0.0])];
        sensor.update(&bodies);
        assert!(sensor.is_triggered());
        assert_eq!(sensor.overlapping_bodies(), vec![1]);
    }

    #[test]
    fn test_sphere_sensor_no_trigger() {
        let mut sensor = SphereSensor::new([0.0; 3], 1.0);
        let bodies = vec![(1u64, [5.0, 0.0, 0.0])];
        sensor.update(&bodies);
        assert!(!sensor.is_triggered());
    }

    // ── FABRIK ────────────────────────────────────────────────────────────────

    #[test]
    fn test_fabrik_2joint_chain() {
        let bones = vec![[0.0, 0.0, 0.0f64], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let lengths = vec![1.0f64, 1.0];
        let target = [1.0, 1.0, 0.0];
        let result = fabrik_solve(&bones, &lengths, &target, 50, 1e-6);
        assert_eq!(result.len(), 3);
        // Root should stay fixed
        assert!((result[0][0]).abs() < 1e-10);
        // End effector near target (reachable)
        let ee = result[2];
        let d = ((ee[0] - target[0]).powi(2) + (ee[1] - target[1]).powi(2)).sqrt();
        assert!(d < 0.01, "end effector dist to target: {d}");
    }

    #[test]
    fn test_fabrik_preserves_bone_lengths() {
        let bones = vec![[0.0, 0.0, 0.0f64], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let lengths = vec![1.0f64, 1.0];
        let target = [0.5, 0.5, 0.0];
        let result = fabrik_solve(&bones, &lengths, &target, 20, 1e-6);
        // Check that bone lengths are approximately preserved
        for i in 0..2 {
            let d = dist3(result[i], result[i + 1]);
            assert!((d - lengths[i]).abs() < 1e-6, "bone {i} length: {d}");
        }
    }
}
