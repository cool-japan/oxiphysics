//! Aerodynamic forces (drag and lift) on cloth/soft-body triangles.
//!
//! The model computes, for each triangle, the relative velocity between the
//! wind and the panel, then splits it into drag (along the panel normal) and
//! lift (perpendicular to both relative velocity and normal).
//!
//! Additional effects:
//! - **Magnus force** -- rotating sphere in a flow.
//! - **Buoyancy** -- Archimedes buoyant force on a submerged volume.
//! - **Panel method** -- basic constant-strength source panel influence.
//! - **Vortex lattice** -- horseshoe vortex lattice element.
//! - **Pressure distribution** -- Bernoulli-based Cp computation.
//! - **Wind load** -- distributed wind load on structural panels.

pub mod functions;
pub(crate) mod helpers;
pub mod types;

// Re-export all public items
pub use functions::*;
pub use types::*;
