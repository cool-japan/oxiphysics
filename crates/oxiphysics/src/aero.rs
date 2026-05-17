// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aerodynamics — velocity-squared drag and airfoil lift/drag forces.
//!
//! This module computes aerodynamic forces for rigid bodies. It supports:
//!
//! - **Point drag** (`AeroBody`): simple velocity-squared drag opposing
//!   relative motion, parameterised by drag coefficient, cross-sectional
//!   area, and air density.
//! - **Wing surfaces** (`WingSurface`): full lift + drag forces computed
//!   from the angle of attack (AoA) against configurable lift/drag curves.
//!   The centre of pressure is tracked so that off-centre wings produce a
//!   torque about the body centre.
//!
//! ## Types
//!
//! | Type | Role |
//! |---|---|
//! | `LiftCurve` | CL vs AoA relationship (linear, stalled, or lookup table) |
//! | `DragCurve` | CD vs AoA relationship (constant or lookup table) |
//! | `AeroBody` | Bluff-body point drag |
//! | `WingSurface` | Airfoil surface with lift and drag |
//! | `AeroEntry` | One body: position, rotation, velocity, drag body, wings |
//! | `AeroForce` | Net force and torque output |
//! | `AeroSystem` | Collection of entries plus global wind; produces `AeroForce` per entry |
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::aero::{
//!     AeroBody, AeroEntry, AeroSystem,
//!     DragCurve, LiftCurve, WingSurface,
//! };
//!
//! // Identity rotation (body axes aligned with world axes)
//! let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
//!
//! // A flat plate moving at 20 m/s into the wind
//! let entry = AeroEntry {
//!     position: [0.0, 0.0, 0.0],
//!     rotation: identity,
//!     velocity: [20.0, 0.0, 0.0],
//!     drag_body: Some(AeroBody {
//!         drag_coeff: 1.28,
//!         cross_section: 0.05,
//!         air_density: 1.225,
//!     }),
//!     wings: vec![WingSurface {
//!         center_local: [0.0, 0.0, 0.0],
//!         normal_local: [0.0, 1.0, 0.0],
//!         chord: 0.3,
//!         span: 1.2,
//!         lift_curve: LiftCurve::Linear { cl0: 0.0, slope: 5.73 },
//!         drag_curve: DragCurve::Constant(0.02),
//!         air_density: 1.225,
//!     }],
//! };
//!
//! let system = AeroSystem { entries: vec![entry], wind: [0.0; 3] };
//! let forces = system.apply();
//! assert_eq!(forces.len(), 1);
//! assert!(forces[0].force[0] < 0.0, "drag opposes +X motion");
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Internal math helpers (no external deps)
// ---------------------------------------------------------------------------

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Returns `None` if the vector length is below 1e-12.
fn normalize(a: [f64; 3]) -> Option<[f64; 3]> {
    let l = len(a);
    if l < 1e-12 {
        None
    } else {
        Some(scale(a, 1.0 / l))
    }
}

/// Transform a local vector into world space using the 3×3 row-major rotation
/// matrix. `rot[i]` is the i-th row (world-space axis of the body).
///
/// `world_v[j] = Σ_i rot[i][j] · local_v[i]`
fn local_to_world(rot: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        rot[0][0] * v[0] + rot[1][0] * v[1] + rot[2][0] * v[2],
        rot[0][1] * v[0] + rot[1][1] * v[1] + rot[2][1] * v[2],
        rot[0][2] * v[0] + rot[1][2] * v[1] + rot[2][2] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// LiftCurve
// ---------------------------------------------------------------------------

/// Lift coefficient (C_L) as a function of angle of attack (AoA, in radians).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LiftCurve {
    /// Purely linear: `CL = cl0 + slope * alpha`.
    Linear {
        /// Zero-AoA lift coefficient.
        cl0: f64,
        /// Lift-curve slope (1/rad).
        slope: f64,
    },
    /// Linear up to stall, then falls at `after_stall_slope`.
    Stalled {
        /// Stall angle (rad).
        alpha_stall_rad: f64,
        /// Maximum lift coefficient at stall.
        cl_max: f64,
        /// Post-stall slope (usually negative, 1/rad).
        after_stall_slope: f64,
    },
    /// Piecewise-linear lookup table: sorted `(alpha_rad, CL)` pairs.
    Lookup {
        /// Must be sorted by `alpha_rad`.
        samples: Vec<(f64, f64)>,
    },
}

impl LiftCurve {
    /// Evaluate C_L at the given angle of attack `alpha` (radians).
    pub fn evaluate(&self, alpha: f64) -> f64 {
        match self {
            LiftCurve::Linear { cl0, slope } => cl0 + slope * alpha,
            LiftCurve::Stalled {
                alpha_stall_rad,
                cl_max,
                after_stall_slope,
            } => {
                if alpha.abs() < *alpha_stall_rad {
                    // Linear ramp from 0 to cl_max over [0, alpha_stall_rad]
                    // (same for negative alpha, mirrored).
                    let slope = cl_max / alpha_stall_rad;
                    slope * alpha
                } else {
                    // Post-stall: continues from cl_max with after_stall_slope.
                    let sign = alpha.signum();
                    sign * cl_max + after_stall_slope * (alpha - sign * alpha_stall_rad)
                }
            }
            LiftCurve::Lookup { samples } => interpolate_lookup(samples, alpha),
        }
    }
}

// ---------------------------------------------------------------------------
// DragCurve
// ---------------------------------------------------------------------------

/// Drag coefficient (C_D) as a function of angle of attack (AoA, in radians).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DragCurve {
    /// Constant drag coefficient independent of AoA.
    Constant(f64),
    /// Piecewise-linear lookup table: sorted `(alpha_rad, CD)` pairs.
    Lookup {
        /// Must be sorted by `alpha_rad`.
        samples: Vec<(f64, f64)>,
    },
}

impl DragCurve {
    /// Evaluate C_D at the given angle of attack `alpha` (radians).
    pub fn evaluate(&self, alpha: f64) -> f64 {
        match self {
            DragCurve::Constant(c) => *c,
            DragCurve::Lookup { samples } => interpolate_lookup(samples, alpha),
        }
    }
}

// ---------------------------------------------------------------------------
// Lookup-table interpolation (shared by Lift and Drag)
// ---------------------------------------------------------------------------

/// Linear interpolation in a sorted `(x, y)` lookup table.
/// Clamps to endpoint values outside the sample range.
fn interpolate_lookup(samples: &[(f64, f64)], x: f64) -> f64 {
    if samples.is_empty() {
        return 0.0;
    }
    if samples.len() == 1 {
        return samples[0].1;
    }
    // Clamp below range
    if x <= samples[0].0 {
        return samples[0].1;
    }
    // Clamp above range
    if x >= samples[samples.len() - 1].0 {
        return samples[samples.len() - 1].1;
    }
    // Binary search for the bracketing interval
    let idx = samples.partition_point(|&(xi, _)| xi <= x);
    // idx is the first index where xi > x, so the bracket is [idx-1, idx]
    let (x0, y0) = samples[idx - 1];
    let (x1, y1) = samples[idx];
    let t = (x - x0) / (x1 - x0);
    y0 + t * (y1 - y0)
}

// ---------------------------------------------------------------------------
// AeroBody (point drag)
// ---------------------------------------------------------------------------

/// Simple bluff-body with velocity-squared drag and no lift.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AeroBody {
    /// Drag coefficient C_d.
    pub drag_coeff: f64,
    /// Reference cross-sectional area (m²).
    pub cross_section: f64,
    /// Air density ρ (kg/m³). Default sea-level value: 1.225.
    pub air_density: f64,
}

impl Default for AeroBody {
    fn default() -> Self {
        Self {
            drag_coeff: 0.47,
            cross_section: 1.0,
            air_density: 1.225,
        }
    }
}

impl AeroBody {
    /// Compute the drag force vector for the given relative velocity.
    ///
    /// `F_drag = -0.5 * ρ * Cd * A * |v_rel| * v_rel`
    fn drag_force(&self, relative_velocity: [f64; 3]) -> [f64; 3] {
        let speed = len(relative_velocity);
        // The factor `speed * relative_velocity` gives `speed² * direction`.
        scale(
            relative_velocity,
            -0.5 * self.air_density * self.drag_coeff * self.cross_section * speed,
        )
    }
}

// ---------------------------------------------------------------------------
// WingSurface
// ---------------------------------------------------------------------------

/// An airfoil wing surface producing both lift and drag forces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WingSurface {
    /// Centre of pressure in body-local space.
    pub center_local: [f64; 3],
    /// Wing normal direction (lift axis) in body-local space.
    pub normal_local: [f64; 3],
    /// Wing chord length (m).
    pub chord: f64,
    /// Wing span (m). Reference area `S = chord * span`.
    pub span: f64,
    /// Lift coefficient curve.
    pub lift_curve: LiftCurve,
    /// Drag coefficient curve.
    pub drag_curve: DragCurve,
    /// Air density ρ (kg/m³).
    pub air_density: f64,
}

impl WingSurface {
    /// Compute the aerodynamic force and the moment arm `r` (world-space
    /// vector from body centre to the wing's centre of pressure).
    ///
    /// Returns `(force, r)` or `([0;3], [0;3])` when the airspeed is too low.
    fn compute(
        &self,
        body_velocity: [f64; 3],
        wind: [f64; 3],
        rotation: &[[f64; 3]; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let relative_v = sub(body_velocity, wind);
        let speed = len(relative_v);
        if speed < 1e-6 {
            return ([0.0; 3], [0.0; 3]);
        }

        // Wind unit direction
        let wind_dir = scale(relative_v, 1.0 / speed);

        // Transform wing normal from local → world
        let world_normal = local_to_world(rotation, self.normal_local);

        // Angle of attack: sin(α) = dot(-wind_dir, world_normal)
        let sin_alpha = dot(scale(wind_dir, -1.0), world_normal).clamp(-1.0, 1.0);
        let alpha = sin_alpha.asin();

        // Dynamic pressure and reference area
        let q = 0.5 * self.air_density * speed * speed;
        let s = self.chord * self.span;

        let cl = self.lift_curve.evaluate(alpha);
        let cd = self.drag_curve.evaluate(alpha);

        // Lift direction: project world_normal onto the plane ⊥ to wind_dir
        // lift_dir = world_normal - (dot(world_normal, wind_dir)) * wind_dir
        let lift_dir_raw = sub(world_normal, scale(wind_dir, dot(world_normal, wind_dir)));

        let lift_force = if let Some(ld) = normalize(lift_dir_raw) {
            scale(ld, cl * q * s)
        } else {
            [0.0; 3]
        };

        let drag_force = scale(wind_dir, -cd * q * s);

        let wing_force = add(lift_force, drag_force);

        // Moment arm: rotation * center_local (already world-space offset from body centre)
        let r = local_to_world(rotation, self.center_local);

        (wing_force, r)
    }
}

// ---------------------------------------------------------------------------
// AeroEntry
// ---------------------------------------------------------------------------

/// One body participating in the aerodynamics system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AeroEntry {
    /// Body centre position in world space.
    pub position: [f64; 3],
    /// Row-major 3×3 orientation matrix.
    /// `rotation[i]` is the i-th row (world-space body axis i).
    pub rotation: [[f64; 3]; 3],
    /// Body velocity in world space (m/s).
    pub velocity: [f64; 3],
    /// Optional bluff-body drag (applied at centre of mass).
    pub drag_body: Option<AeroBody>,
    /// Wing surfaces attached to this body.
    pub wings: Vec<WingSurface>,
}

// ---------------------------------------------------------------------------
// AeroForce
// ---------------------------------------------------------------------------

/// Aerodynamic force and torque output for one body.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AeroForce {
    /// Net aerodynamic force in world space (N).
    pub force: [f64; 3],
    /// Net aerodynamic torque about the body centre in world space (N·m).
    pub torque: [f64; 3],
}

impl AeroForce {
    fn zero() -> Self {
        Self {
            force: [0.0; 3],
            torque: [0.0; 3],
        }
    }
}

// ---------------------------------------------------------------------------
// AeroSystem
// ---------------------------------------------------------------------------

/// The aerodynamics system: holds all entries and a global wind vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AeroSystem {
    /// All bodies registered with the system.
    pub entries: Vec<AeroEntry>,
    /// Global wind velocity in world space (m/s).
    pub wind: [f64; 3],
}

impl AeroSystem {
    /// Compute one [`AeroForce`] per [`AeroEntry`].
    pub fn apply(&self) -> Vec<AeroForce> {
        self.entries
            .iter()
            .map(|entry| {
                let mut total_force = [0.0; 3];
                let mut total_torque = [0.0; 3];

                // --- Point drag ---
                if let Some(body) = &entry.drag_body {
                    let rel_v = sub(entry.velocity, self.wind);
                    let f = body.drag_force(rel_v);
                    total_force = add(total_force, f);
                    // Point drag acts at CoM → no torque contribution.
                }

                // --- Wing surfaces ---
                for wing in &entry.wings {
                    let (f, r) = wing.compute(entry.velocity, self.wind, &entry.rotation);
                    total_force = add(total_force, f);
                    total_torque = add(total_torque, cross(r, f));
                }

                AeroForce {
                    force: total_force,
                    torque: total_torque,
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: identity rotation (body axes == world axes)
    fn identity_rot() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    // Helper: a simple AeroSystem with a drag body only, no wings.
    fn drag_only_system(velocity: [f64; 3], cd: f64, area: f64) -> AeroSystem {
        AeroSystem {
            entries: vec![AeroEntry {
                position: [0.0; 3],
                rotation: identity_rot(),
                velocity,
                drag_body: Some(AeroBody {
                    drag_coeff: cd,
                    cross_section: area,
                    air_density: 1.225,
                }),
                wings: vec![],
            }],
            wind: [0.0; 3],
        }
    }

    // -----------------------------------------------------------------------
    // 1. Zero velocity → zero force
    // -----------------------------------------------------------------------
    #[test]
    fn test_zero_velocity_zero_force() {
        let sys = drag_only_system([0.0; 3], 1.0, 1.0);
        let forces = sys.apply();
        assert_eq!(forces.len(), 1);
        let f = forces[0].force;
        assert!(
            f[0].abs() < 1e-12 && f[1].abs() < 1e-12 && f[2].abs() < 1e-12,
            "Zero velocity should produce zero drag force, got {f:?}"
        );
    }

    // -----------------------------------------------------------------------
    // 2. Drag opposes velocity (body moving in +X → drag in -X)
    // -----------------------------------------------------------------------
    #[test]
    fn test_drag_opposes_velocity() {
        let sys = drag_only_system([10.0, 0.0, 0.0], 1.0, 1.0);
        let forces = sys.apply();
        assert!(
            forces[0].force[0] < 0.0,
            "Drag should oppose +X motion, got Fx={}",
            forces[0].force[0]
        );
    }

    // -----------------------------------------------------------------------
    // 3. Drag scales as v² (ratio of forces at 10 m/s vs 20 m/s ≈ 4)
    // -----------------------------------------------------------------------
    #[test]
    fn test_drag_scales_v_squared() {
        let sys10 = drag_only_system([10.0, 0.0, 0.0], 1.0, 1.0);
        let sys20 = drag_only_system([20.0, 0.0, 0.0], 1.0, 1.0);

        let f10 = sys10.apply()[0].force[0].abs();
        let f20 = sys20.apply()[0].force[0].abs();

        let ratio = f20 / f10;
        assert!(
            (ratio - 4.0).abs() < 1e-9,
            "Drag force ratio at v=20 vs v=10 should be 4.0, got {ratio}"
        );
    }

    // -----------------------------------------------------------------------
    // 4. Wing CL=0 at α=0 (Linear cl0=0) → lift ≈ 0, only drag
    // -----------------------------------------------------------------------
    #[test]
    fn test_wing_zero_cl_at_zero_alpha() {
        // Normal pointing straight up, velocity in +X → AoA = 0, CL = 0.
        let sys = AeroSystem {
            entries: vec![AeroEntry {
                position: [0.0; 3],
                rotation: identity_rot(),
                velocity: [10.0, 0.0, 0.0],
                drag_body: None,
                wings: vec![WingSurface {
                    center_local: [0.0; 3],
                    normal_local: [0.0, 1.0, 0.0],
                    chord: 1.0,
                    span: 1.0,
                    lift_curve: LiftCurve::Linear {
                        cl0: 0.0,
                        slope: 5.73,
                    },
                    drag_curve: DragCurve::Constant(0.02),
                    air_density: 1.225,
                }],
            }],
            wind: [0.0; 3],
        };
        let forces = sys.apply();
        // Lift (Y component) should be near zero
        assert!(
            forces[0].force[1].abs() < 1e-6,
            "Lift should be ≈0 at α=0 with cl0=0, got Fy={}",
            forces[0].force[1]
        );
        // Drag (X component) should be negative (opposing motion)
        assert!(
            forces[0].force[0] < 0.0,
            "Drag should be negative, got Fx={}",
            forces[0].force[0]
        );
    }

    // -----------------------------------------------------------------------
    // 5. Wing lift direction: normal [0,1,0], body moving slightly downward → positive Fy
    // -----------------------------------------------------------------------
    #[test]
    fn test_wing_lift_direction() {
        // Wing normal up [0,1,0].  The formula computes:
        //   sin_alpha = dot(-wind_dir, world_normal)
        // For alpha > 0 (positive lift) we need:
        //   dot(-wind_dir, [0,1,0]) > 0  →  wind_dir.y < 0
        // i.e. the body moves slightly *downward* (nose-down AoA into the wing).
        let speed = 10.0_f64;
        let angle = 0.2_f64; // 0.2 rad below horizontal
        let sys = AeroSystem {
            entries: vec![AeroEntry {
                position: [0.0; 3],
                rotation: identity_rot(),
                velocity: [speed * angle.cos(), -speed * angle.sin(), 0.0],
                drag_body: None,
                wings: vec![WingSurface {
                    center_local: [0.0; 3],
                    normal_local: [0.0, 1.0, 0.0],
                    chord: 1.0,
                    span: 1.0,
                    lift_curve: LiftCurve::Linear {
                        cl0: 0.0,
                        slope: 5.73,
                    },
                    drag_curve: DragCurve::Constant(0.01),
                    air_density: 1.225,
                }],
            }],
            wind: [0.0; 3],
        };
        let forces = sys.apply();
        // α = asin(sin(angle)) ≈ 0.2 rad → CL > 0 → lift direction ≈ +Y.
        assert!(
            forces[0].force[1] > 0.0,
            "Lift force should have positive Y component for nose-down AoA, got Fy={}",
            forces[0].force[1]
        );
    }

    // -----------------------------------------------------------------------
    // 6. Stall: CL at alpha > alpha_stall < CL at stall
    // -----------------------------------------------------------------------
    #[test]
    fn test_lift_curve_stall() {
        let curve = LiftCurve::Stalled {
            alpha_stall_rad: 0.3,
            cl_max: 1.5,
            after_stall_slope: -2.0,
        };
        let cl_at_stall = curve.evaluate(0.3);
        let cl_post_stall = curve.evaluate(0.5);
        assert!(
            cl_post_stall < cl_at_stall,
            "Post-stall CL ({cl_post_stall}) should be less than CL at stall ({cl_at_stall})"
        );
    }

    // -----------------------------------------------------------------------
    // 7. Lookup interpolation: (0.0, 0.5) and (1.0, 1.5) → evaluate at 0.5 ≈ 1.0
    // -----------------------------------------------------------------------
    #[test]
    fn test_lookup_interpolation() {
        let curve = LiftCurve::Lookup {
            samples: vec![(0.0, 0.5), (1.0, 1.5)],
        };
        let val = curve.evaluate(0.5);
        assert!(
            (val - 1.0).abs() < 1e-10,
            "Linear interpolation at midpoint should give 1.0, got {val}"
        );
    }

    // -----------------------------------------------------------------------
    // 8. Serde round-trip of AeroSystem with a wing entry
    // -----------------------------------------------------------------------
    #[test]
    fn test_serde_round_trip() {
        let original = AeroSystem {
            entries: vec![AeroEntry {
                position: [1.0, 2.0, 3.0],
                rotation: identity_rot(),
                velocity: [5.0, 0.0, 0.0],
                drag_body: Some(AeroBody {
                    drag_coeff: 0.3,
                    cross_section: 0.8,
                    air_density: 1.225,
                }),
                wings: vec![WingSurface {
                    center_local: [0.5, 0.0, 0.0],
                    normal_local: [0.0, 1.0, 0.0],
                    chord: 0.4,
                    span: 2.0,
                    lift_curve: LiftCurve::Lookup {
                        samples: vec![(0.0, 0.0), (0.3, 1.2), (0.5, 0.8)],
                    },
                    drag_curve: DragCurve::Lookup {
                        samples: vec![(0.0, 0.01), (0.5, 0.05)],
                    },
                    air_density: 1.225,
                }],
            }],
            wind: [2.0, 0.0, 0.0],
        };

        let json = serde_json::to_string(&original).expect("serialize AeroSystem");
        let restored: AeroSystem = serde_json::from_str(&json).expect("deserialize AeroSystem");

        assert_eq!(restored.entries.len(), original.entries.len());
        assert_eq!(restored.wind, original.wind);
        let e = &restored.entries[0];
        assert_eq!(e.position, [1.0, 2.0, 3.0]);
        assert_eq!(e.wings.len(), 1);
        assert!((e.wings[0].chord - 0.4).abs() < 1e-15);
    }
}
