//! # ParrinelloRahmanBarostat - Trait Implementations
//!
//! This module contains trait implementations for `ParrinelloRahmanBarostat`.
//!
//! ## Implemented Traits
//!
//! - `Barostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::{Barostat, virial_pressure};
use super::types::ParrinelloRahmanBarostat;

impl Barostat for ParrinelloRahmanBarostat {
    /// One extended-Lagrangian step of the isotropic box edge `L = V^(1/3)`.
    ///
    /// With the instantaneous pressure `P = (2·KE + W) / (3·V)` (pair virial `W`
    /// supplied by the caller) and the target pressure `P₀` (the
    /// `target_pressure` argument if finite, the barostat's own setting
    /// otherwise):
    ///
    /// 1. `L̇ += 3·V·(P − P₀) / W_box · dt` — the box velocity is persistent
    ///    state, so the box has inertia;
    /// 2. `L_new = L + L̇·dt`;
    /// 3. the box and all positions are scaled by `L_new / L`;
    /// 4. all velocities are scaled by `L / L_new` (conservation of the
    ///    canonical momentum of the scaled coordinates).
    ///
    /// The barostat mass `W_box` is `mass`, or the value derived from `tau_p`
    /// when `mass` is not positive (see
    /// [`ParrinelloRahmanBarostat::effective_mass`]).  If no mass can be
    /// determined, or the step would collapse the box, the box is left
    /// unchanged (a collapsing step also zeroes the box velocity).
    fn apply(
        &mut self,
        atoms: &mut AtomSet,
        pbox: &mut PeriodicBox,
        target_pressure: f64,
        virial: f64,
        dt: f64,
    ) {
        let volume = pbox.volume();
        if volume < 1e-30 {
            return;
        }
        let ke = atoms.kinetic_energy();
        let p_current = virial_pressure(volume, ke, virial);
        let p_target = if target_pressure.is_finite() {
            target_pressure
        } else {
            self.target_pressure
        };
        let edge = volume.cbrt();
        let Some(box_mass) = self.effective_mass(ke, edge) else {
            return;
        };
        let new_velocity =
            self.box_velocity + 3.0 * volume * (p_current - p_target) / box_mass * dt;
        let new_edge = edge + new_velocity * dt;
        if !new_edge.is_finite() || new_edge <= 0.0 {
            self.box_velocity = 0.0;
            return;
        }
        self.box_velocity = new_velocity;
        let scale = new_edge / edge;
        pbox.dims *= scale;
        for pos in &mut atoms.positions {
            *pos *= scale;
        }
        let inverse_scale = 1.0 / scale;
        for vel in &mut atoms.velocities {
            *vel *= inverse_scale;
        }
    }
}
