// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Archimedes buoyancy and viscous drag for bodies immersed in fluid volumes.
//!
//! ## Physics Model
//!
//! For each body–fluid pair:
//!
//! - **Buoyancy force** (upward): `F_b = ρ_f · V_sub · g`
//!   where `ρ_f` is the fluid density, `V_sub` is the submerged volume, and
//!   `g` is gravitational acceleration.
//! - **Linear drag** (opposing velocity): `F_d = -c_lin · v`
//! - **Angular drag** (opposing rotation): `T_d = -c_ang · ω`
//!
//! Submerged volume is estimated analytically:
//! - **Sphere**: spherical-cap formula `V_cap = π h² (3r − h) / 3`
//!   where `h = clamp(surface_y − (cy − r), 0, 2r)`.
//! - **Box**: axis-aligned fraction `(clamp(surface_y − bottom, 0, 2·hy)) / (2·hy)`.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::buoyancy::{BuoyancyWorld, BuoyantBody, BuoyantShape};
//!
//! let mut world = BuoyancyWorld::new(9.81);
//! // Ocean surface at y = 0, depth down to y = -50
//! world.add_fluid(
//!     [-100.0, -50.0, -100.0],
//!     [ 100.0,   0.0,  100.0],
//!     1025.0, // seawater density  (kg/m³)
//!     0.0,    // surface y
//! );
//!
//! let buoy = BuoyantBody {
//!     id: 42,
//!     position: [0.0, -0.3, 0.0], // partially submerged sphere
//!     velocity:  [0.5, 0.0, 0.0],
//!     angular_velocity: [0.0; 3],
//!     mass: 12.0,
//!     shape: BuoyantShape::Sphere { radius: 0.5 },
//! };
//!
//! let forces = world.apply(&[buoy]);
//! assert!(!forces.is_empty());
//! let (_, bf) = &forces[0];
//! assert!(bf.force[1] > 0.0);           // net upward buoyancy
//! assert!(bf.submerged_fraction > 0.0); // partially submerged
//! ```

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// BuoyantShape
// ---------------------------------------------------------------------------

/// Shape used for displaced-volume estimation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BuoyantShape {
    /// Sphere with the given radius (m).
    Sphere { radius: f64 },
    /// Axis-aligned box with half-extents `[x, y, z]` (m).
    Box { half_extents: [f64; 3] },
}

// ---------------------------------------------------------------------------
// BuoyantBody
// ---------------------------------------------------------------------------

/// A rigid body that can receive buoyancy and drag forces from fluid volumes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuoyantBody {
    /// Unique body identifier.
    pub id: u32,
    /// World-space centroid position (m).
    pub position: [f64; 3],
    /// Linear velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
    /// Body mass (kg).
    pub mass: f64,
    /// Shape used for submerged volume estimation.
    pub shape: BuoyantShape,
}

// ---------------------------------------------------------------------------
// FluidVolume
// ---------------------------------------------------------------------------

/// A rectangular fluid region with uniform physical properties.
///
/// The region is defined by an AABB; the actual fluid surface is given by
/// `surface_y`, which may be at or below `max[1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidVolume {
    /// Unique volume identifier.
    pub id: u32,
    /// AABB minimum corner (m).
    pub min: [f64; 3],
    /// AABB maximum corner (m).
    pub max: [f64; 3],
    /// Fluid density (kg/m³). Water ≈ 1 000, seawater ≈ 1 025, oil ≈ 900,
    /// mercury ≈ 13 600.
    pub density: f64,
    /// World-space Y coordinate of the fluid surface.
    ///
    /// Bodies above this Y level receive no buoyancy. Defaults to `max[1]`.
    pub surface_y: f64,
    /// Linear drag coefficient (N·s/m). Typical range: 0.1 – 2.0.
    pub linear_drag: f64,
    /// Angular drag coefficient (N·m·s/rad). Typical range: 0.05 – 0.5.
    pub angular_drag: f64,
}

impl FluidVolume {
    /// Create a fluid volume with default drag coefficients (0.5 / 0.1).
    pub fn new(min: [f64; 3], max: [f64; 3], density: f64) -> Self {
        let surface_y = max[1];
        Self {
            id: 0,
            min,
            max,
            density,
            surface_y,
            linear_drag: 0.5,
            angular_drag: 0.1,
        }
    }

    /// Return `true` if the world-space point lies inside the AABB (ignoring
    /// `surface_y` — use `BuoyancyForce::submerged_fraction` for surface-aware tests).
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
}

// ---------------------------------------------------------------------------
// BuoyancyForce
// ---------------------------------------------------------------------------

/// The net buoyancy + drag force and torque acting on a body due to one or
/// more fluid volumes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuoyancyForce {
    /// Net force in world space (N): buoyancy (upward) + linear drag.
    pub force: [f64; 3],
    /// Net torque in world space (N·m): angular drag opposing `ω`.
    pub torque: [f64; 3],
    /// Fraction of the body's volume submerged in [0, 1].
    pub submerged_fraction: f64,
}

// ---------------------------------------------------------------------------
// Submerged-volume maths (private helpers)
// ---------------------------------------------------------------------------

/// Volume of a spherical cap of height `h` on a sphere of radius `r`.
///
/// `h` is clamped to `[0, 2r]`; `h = 2r` gives the full sphere.
#[inline]
fn sphere_cap_volume(r: f64, h: f64) -> f64 {
    let h = h.clamp(0.0, 2.0 * r);
    PI * h * h * (3.0 * r - h) / 3.0
}

/// Full sphere volume.
#[inline]
fn sphere_volume(r: f64) -> f64 {
    4.0 / 3.0 * PI * r * r * r
}

/// Submerged volume fraction of a sphere (centre Y = `cy`, radius `r`) given
/// a fluid surface at `surface_y`.
///
/// Returns a value in `[0.0, 1.0]`.
fn sphere_submerged_fraction(cy: f64, r: f64, surface_y: f64) -> f64 {
    let bottom = cy - r;
    let top = cy + r;
    if top <= surface_y {
        // Fully submerged.
        return 1.0;
    }
    if bottom >= surface_y {
        // Fully above the surface.
        return 0.0;
    }
    // Depth of the submerged cap = distance from bottom to surface.
    let h = surface_y - bottom;
    let v_sub = sphere_cap_volume(r, h);
    let v_total = sphere_volume(r);
    (v_sub / v_total).clamp(0.0, 1.0)
}

/// Submerged volume fraction of a box (centre Y = `cy`, half-height `hy`).
fn box_submerged_fraction(cy: f64, hy: f64, surface_y: f64) -> f64 {
    let bottom = cy - hy;
    let top = cy + hy;
    if top <= surface_y {
        return 1.0;
    }
    if bottom >= surface_y {
        return 0.0;
    }
    let submerged_height = (surface_y - bottom).max(0.0);
    (submerged_height / (2.0 * hy)).clamp(0.0, 1.0)
}

/// Compute `(total_volume_m3, submerged_fraction)` for the body in the fluid.
fn displaced_volume(body: &BuoyantBody, fluid: &FluidVolume) -> (f64, f64) {
    match &body.shape {
        BuoyantShape::Sphere { radius } => {
            let r = *radius;
            let total = sphere_volume(r);
            let frac = sphere_submerged_fraction(body.position[1], r, fluid.surface_y);
            (total, frac)
        }
        BuoyantShape::Box { half_extents } => {
            let [hx, hy, hz] = *half_extents;
            let total = 8.0 * hx * hy * hz;
            let frac = box_submerged_fraction(body.position[1], hy, fluid.surface_y);
            (total, frac)
        }
    }
}

// ---------------------------------------------------------------------------
// compute_buoyancy — free function
// ---------------------------------------------------------------------------

/// Compute the buoyancy and drag forces exerted on `body` by `fluid`.
///
/// Returns `None` when the body's centroid is outside the fluid AABB
/// horizontally, or when the computed submerged fraction is zero.
///
/// `gravity` is the gravitational acceleration magnitude (m/s², typically 9.81).
pub fn compute_buoyancy(
    body: &BuoyantBody,
    fluid: &FluidVolume,
    gravity: f64,
) -> Option<BuoyancyForce> {
    // Fast horizontal AABB reject — Y overlap is handled by surface_y.
    let p = body.position;
    if p[0] < fluid.min[0] || p[0] > fluid.max[0] || p[2] < fluid.min[2] || p[2] > fluid.max[2] {
        return None;
    }

    let (total_vol, frac) = displaced_volume(body, fluid);
    if frac <= 0.0 {
        return None;
    }

    let v_sub = total_vol * frac;

    // Archimedes upward force.
    let f_buoy = fluid.density * v_sub * gravity;

    // Linear drag: F_drag = -c_lin · v.
    let [vx, vy, vz] = body.velocity;
    let c = fluid.linear_drag;
    let force = [-c * vx, f_buoy - c * vy, -c * vz];

    // Angular drag: T_drag = -c_ang · ω.
    let [wx, wy, wz] = body.angular_velocity;
    let ca = fluid.angular_drag;
    let torque = [-ca * wx, -ca * wy, -ca * wz];

    Some(BuoyancyForce {
        force,
        torque,
        submerged_fraction: frac,
    })
}

// ---------------------------------------------------------------------------
// BuoyancyWorld
// ---------------------------------------------------------------------------

/// World container managing fluid volumes and computing buoyancy forces.
///
/// ```rust,no_run
/// use oxiphysics::buoyancy::{BuoyancyWorld, BuoyantBody, BuoyantShape};
///
/// let mut world = BuoyancyWorld::new(9.81);
/// world.add_fluid([-10.0, -5.0, -10.0], [10.0, 0.0, 10.0], 1000.0, 0.0);
///
/// let body = BuoyantBody {
///     id: 1,
///     position: [0.0, -0.5, 0.0],
///     velocity: [0.0; 3],
///     angular_velocity: [0.0; 3],
///     mass: 5.0,
///     shape: BuoyantShape::Sphere { radius: 0.5 },
/// };
/// let result = world.apply(&[body]);
/// println!("{result:?}");
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuoyancyWorld {
    fluids: Vec<FluidVolume>,
    next_id: u32,
    /// Gravitational acceleration magnitude (m/s²).  Typically 9.81.
    pub gravity: f64,
}

impl BuoyancyWorld {
    /// Create a world with the given gravitational acceleration.
    pub fn new(gravity: f64) -> Self {
        Self {
            fluids: Vec::new(),
            next_id: 0,
            gravity,
        }
    }

    // -----------------------------------------------------------------------
    // Fluid volume management
    // -----------------------------------------------------------------------

    /// Register a fluid volume and return its ID.
    ///
    /// Default drag coefficients: `linear_drag = 0.5`, `angular_drag = 0.1`.
    pub fn add_fluid(&mut self, min: [f64; 3], max: [f64; 3], density: f64, surface_y: f64) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.fluids.push(FluidVolume {
            id,
            min,
            max,
            density,
            surface_y,
            linear_drag: 0.5,
            angular_drag: 0.1,
        });
        id
    }

    /// Register a fluid volume with explicit drag coefficients.
    pub fn add_fluid_with_drag(
        &mut self,
        min: [f64; 3],
        max: [f64; 3],
        density: f64,
        surface_y: f64,
        linear_drag: f64,
        angular_drag: f64,
    ) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.fluids.push(FluidVolume {
            id,
            min,
            max,
            density,
            surface_y,
            linear_drag,
            angular_drag,
        });
        id
    }

    /// Remove a fluid volume by ID.  No-op if not found.
    pub fn remove_fluid(&mut self, id: u32) {
        self.fluids.retain(|f| f.id != id);
    }

    /// Immutable reference to a fluid volume by ID.
    pub fn fluid(&self, id: u32) -> Option<&FluidVolume> {
        self.fluids.iter().find(|f| f.id == id)
    }

    /// Mutable reference to a fluid volume by ID.
    pub fn fluid_mut(&mut self, id: u32) -> Option<&mut FluidVolume> {
        self.fluids.iter_mut().find(|f| f.id == id)
    }

    /// All fluid volumes whose AABB contains the world-space point `p`.
    pub fn fluids_at_point(&self, p: [f64; 3]) -> Vec<&FluidVolume> {
        self.fluids.iter().filter(|f| f.contains_point(p)).collect()
    }

    /// Number of registered fluid volumes.
    pub fn fluid_count(&self) -> usize {
        self.fluids.len()
    }

    /// Iterate over all fluid volumes.
    pub fn fluids(&self) -> impl Iterator<Item = &FluidVolume> {
        self.fluids.iter()
    }

    // -----------------------------------------------------------------------
    // Force computation
    // -----------------------------------------------------------------------

    /// Compute the combined buoyancy + drag force for every body against every
    /// overlapping fluid volume.
    ///
    /// Returns `(body_id, combined_force)` for bodies that overlap at least one
    /// fluid volume.  Bodies outside all fluids are omitted.
    pub fn apply(&self, bodies: &[BuoyantBody]) -> Vec<(u32, BuoyancyForce)> {
        let g = self.gravity;
        let mut out = Vec::new();
        for body in bodies {
            let mut combined = BuoyancyForce {
                force: [0.0; 3],
                torque: [0.0; 3],
                submerged_fraction: 0.0,
            };
            let mut hit = false;
            for fluid in &self.fluids {
                if let Some(bf) = compute_buoyancy(body, fluid, g) {
                    combined.force[0] += bf.force[0];
                    combined.force[1] += bf.force[1];
                    combined.force[2] += bf.force[2];
                    combined.torque[0] += bf.torque[0];
                    combined.torque[1] += bf.torque[1];
                    combined.torque[2] += bf.torque[2];
                    // Track the worst-case (maximum) submerged fraction across
                    // overlapping volumes.
                    combined.submerged_fraction =
                        combined.submerged_fraction.max(bf.submerged_fraction);
                    hit = true;
                }
            }
            if hit {
                out.push((body.id, combined));
            }
        }
        out
    }

    /// Apply buoyancy to a single body and return the force, or `None` if the
    /// body does not overlap any fluid volume.
    pub fn apply_single(&self, body: &BuoyantBody) -> Option<BuoyancyForce> {
        self.apply(std::slice::from_ref(body))
            .into_iter()
            .next()
            .map(|(_, bf)| bf)
    }

    /// Estimate whether a body is currently floating at the surface.
    ///
    /// A body is considered floating when its submerged fraction is in
    /// `(0.0, threshold)`.  A threshold of `0.6` works well for spheres.
    pub fn is_floating(&self, body: &BuoyantBody, threshold: f64) -> bool {
        if let Some(bf) = self.apply_single(body) {
            bf.submerged_fraction > 0.0 && bf.submerged_fraction < threshold
        } else {
            false
        }
    }
}
