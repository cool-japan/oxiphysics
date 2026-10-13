// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH-rigid coupling: Archimedes buoyancy from a fluid region.
//!
//! [`SphRigidCouplingParams`] describes the region and its fluid;
//! [`PhysicsPipeline`] applies the buoyancy of that fluid to the rigid bodies
//! that dip into it in every step. The force law is documented on
//! [`SphRigidCouplingParams`].

use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_core::{Aabb, BodyHandle};
use oxiphysics_rigid::{BodyState, RigidBodySet};

use super::PhysicsPipeline;
use super::collision::{ColliderEntry, world_aabb};

/// Descriptor for a region where an SPH fluid couples with rigid bodies.
///
/// The region is an axis-aligned box filled with fluid of density
/// [`fluid_density`](Self::fluid_density). While the descriptor is installed
/// in a [`PhysicsPipeline`] (see
/// [`PhysicsPipeline::set_sph_rigid_coupling`]) every step applies the
/// buoyancy of that fluid to the dynamic bodies that dip into the region:
///
/// ```text
/// F = coupling_alpha * fluid_density * |g| * V_displaced,   directed along -g
/// ```
///
/// `V_displaced` is the volume of the collider's shape scaled by the fraction
/// of the collider's world-space bounding box that lies inside the region, so
/// a body that is half out of the fluid is lifted by half the force (exact for
/// boxes aligned with the region, a first-order estimate for round shapes).
///
/// The pipeline does **not** own an SPH simulation: the reverse direction of
/// the coupling, pushing SPH particles away from rigid body surfaces, is not
/// performed, and `fluid_density` is the density you provide (for example an
/// average sampled from your SPH particles) rather than a value the pipeline
/// measures.
#[derive(Debug, Clone)]
pub struct SphRigidCouplingParams {
    /// Average fluid density in the coupling region \[kg/m³\].
    pub fluid_density: Real,
    /// Coupling strength coefficient (0 = no coupling, 1 = full).
    pub coupling_alpha: Real,
    /// Half-extents of the axis-aligned coupling region.
    pub region_half_extents: Vec3,
    /// Centre of the coupling region in world space.
    pub region_centre: Vec3,
}

impl Default for SphRigidCouplingParams {
    fn default() -> Self {
        Self {
            fluid_density: 1000.0,
            coupling_alpha: 1.0,
            region_half_extents: Vec3::new(5.0, 5.0, 5.0),
            region_centre: Vec3::zeros(),
        }
    }
}

impl SphRigidCouplingParams {
    /// Lower and upper corner of the coupling region.
    fn region_bounds(&self) -> (Vec3, Vec3) {
        let half = self.region_half_extents.abs();
        (self.region_centre - half, self.region_centre + half)
    }
}

/// Fraction of the volume of `aabb` that lies inside the box `[low, high]`.
///
/// A degenerate (flat) extent counts as fully inside if it lies within the
/// bounds along that axis.
fn overlap_fraction(aabb: &Aabb, low: &Vec3, high: &Vec3) -> Real {
    let mut fraction = 1.0;
    for ((lo, hi), (region_lo, region_hi)) in aabb
        .min
        .iter()
        .zip(aabb.max.iter())
        .zip(low.iter().zip(high.iter()))
    {
        let extent = hi - lo;
        let inside = if extent > 1e-12 {
            (hi.min(*region_hi) - lo.max(*region_lo)).max(0.0) / extent
        } else if lo >= region_lo && lo <= region_hi {
            1.0
        } else {
            0.0
        };
        fraction *= inside;
        if fraction == 0.0 {
            break;
        }
    }
    fraction
}

/// Apply buoyancy from a uniform fluid region to all dynamic bodies.
///
/// Every dynamic body whose centre lies inside the coupling region receives an
/// upward (+Y) buoyancy impulse that corresponds to displacing a sphere of
/// `effective_radius`, scaled by the fluid density and `coupling_alpha`.
///
/// This is a stand-alone helper that treats every body as the same sphere. The
/// pipeline itself ([`PhysicsPipeline::step`]) does not call it: with
/// [`PhysicsPipeline::set_sph_rigid_coupling`] it applies the buoyancy of each
/// collider's real shape volume along the direction opposite to gravity
/// instead.
///
/// `g_magnitude` – magnitude of gravitational acceleration \[m/s²\].
pub fn apply_buoyancy_impulse(
    bodies: &mut RigidBodySet,
    coupling: &SphRigidCouplingParams,
    effective_radius: Real,
    g_magnitude: Real,
    dt: Real,
) {
    let half = coupling.region_half_extents;
    let centre = coupling.region_centre;
    let rho_f = coupling.fluid_density;
    let alpha = coupling.coupling_alpha;

    // Volume of a sphere: (4/3) π r³
    let submerged_vol = (4.0 / 3.0) * std::f64::consts::PI * effective_radius.powi(3);
    let buoyancy_force_mag = rho_f * g_magnitude * submerged_vol * alpha;

    for (_, body) in bodies.iter_mut() {
        if !body.is_dynamic() {
            continue;
        }
        let p = body.transform.position;
        // Check if inside coupling region.
        if (p.x - centre.x).abs() <= half.x
            && (p.y - centre.y).abs() <= half.y
            && (p.z - centre.z).abs() <= half.z
        {
            // Buoyancy acts upward (+Y).
            let impulse = Vec3::new(0.0, buoyancy_force_mag * dt, 0.0);
            body.velocity += impulse * body.inverse_mass;
        }
    }
}

impl PhysicsPipeline {
    /// Apply the buoyancy of the coupling region to the bodies of `entries`.
    ///
    /// Each non-sensor collider of a dynamic body that overlaps the region
    /// displaces `shape volume * overlap fraction` of fluid and receives the
    /// corresponding buoyancy impulse at its body's centre of mass, directed
    /// against gravity. Without gravity there is no buoyancy.
    ///
    /// A sleeping body is left alone unless the buoyancy exceeds its weight, in
    /// which case it is woken up so it can rise.
    pub(super) fn apply_sph_rigid_coupling(
        &self,
        params: &SphRigidCouplingParams,
        dt: Real,
        bodies: &mut RigidBodySet,
        entries: &[ColliderEntry<'_>],
    ) {
        let gravity = self.config.gravity;
        let g_magnitude = gravity.norm();
        if !(g_magnitude.is_finite() && g_magnitude > 1e-12 && dt.is_finite() && dt > 0.0) {
            return;
        }
        let up = -gravity / g_magnitude;
        let (low, high) = params.region_bounds();
        let density = params.fluid_density.max(0.0);
        let alpha = params.coupling_alpha.max(0.0);

        for entry in entries {
            let collider = entry.collider;
            let Some(handle) = collider.body_handle else {
                continue;
            };
            if collider.is_sensor {
                continue;
            }
            let Some(body) = bodies.get(handle) else {
                continue;
            };
            if !body.is_dynamic() || body.inverse_mass <= 0.0 {
                continue;
            }
            let world = collider.world_transform(&body.transform);
            let fraction = overlap_fraction(&world_aabb(collider, &world), &low, &high);
            if fraction <= 0.0 {
                continue;
            }
            let volume = collider.shape.volume();
            if !volume.is_finite() || volume <= 0.0 {
                continue;
            }

            let force = alpha * density * g_magnitude * volume * fraction;
            let weight = g_magnitude * body.gravity_scale.max(0.0) / body.inverse_mass;
            let sleeping = body.state == BodyState::Sleeping;
            if sleeping && force <= weight {
                continue;
            }
            Self::push_body(bodies, handle, up * (force * dt), sleeping);
        }
    }

    /// Add `impulse` to the linear momentum of body `handle`, waking it if
    /// `wake` is set.
    fn push_body(bodies: &mut RigidBodySet, handle: BodyHandle, impulse: Vec3, wake: bool) {
        if let Some(body) = bodies.get_mut(handle) {
            body.velocity += impulse * body.inverse_mass;
            if wake {
                body.wake_up();
            }
        }
    }
}
