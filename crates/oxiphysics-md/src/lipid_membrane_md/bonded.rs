// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Periodic-box aware bonded interactions of [`CgLipid`].
//!
//! The box of a [`BilayerSimulation`](super::BilayerSimulation) is periodic in
//! x and y.  A bond or angle of a molecule that straddles the boundary has to be
//! evaluated with the minimum image of its bead-to-bead vectors, otherwise the
//! raw coordinate difference is almost a box length and the harmonic bond
//! returns an enormous spurious energy and force.  This module provides
//!
//! - the bonded energies with minimum-image vectors
//!   ([`CgLipid::bond_energy_in_box`], [`CgLipid::angle_energy_in_box`],
//!   [`CgLipid::intramolecular_energy_in_box`]),
//! - the matching analytic forces ([`CgLipid::intramolecular_forces`]), the exact
//!   negative gradient of [`CgLipid::intramolecular_energy_in_box`], and
//! - [`CgLipid::wrap_into_box`], which keeps molecules whole while they diffuse
//!   through the boundary.
//!
//! Only x and y are periodic; z is not wrapped.

use super::{CgLipid, SimBox};

/// Equilibrium length (nm) of the harmonic bonds.
pub(super) const BOND_R0_NM: f64 = 0.47;
/// Force constant (kJ/(mol nm²)) of the harmonic bonds.
pub(super) const BOND_FORCE_CONSTANT: f64 = 3800.0;
/// Force constant (kJ/mol) of the cosine angle potentials.
pub(super) const ANGLE_FORCE_CONSTANT: f64 = 25.0;
/// `cos θ₀` of the angle potentials (straight chain, `θ₀ = 180°`).
pub(super) const ANGLE_COS_THETA0: f64 = -1.0;
/// Shortest bead-to-bead distance (nm) for which an angle is defined.
const MIN_ANGLE_ARM_NM: f64 = 1.0e-12;

/// Multiple of the period `length` that brings the displacement `d` into
/// `[-length/2, length/2]`; zero for a degenerate period.
fn minimum_image_shift(d: f64, length: f64) -> f64 {
    if length > 0.0 && length.is_finite() {
        length * (d / length).round()
    } else {
        0.0
    }
}

/// Multiple of the period `length` that brings the coordinate `x` into
/// `[0, length)`; zero for a degenerate period.
fn wrap_shift(x: f64, length: f64) -> f64 {
    if !(length > 0.0 && length.is_finite()) {
        return 0.0;
    }
    let mut shift = length * (x / length).floor();
    if x - shift >= length {
        shift += length;
    }
    if x - shift < 0.0 {
        shift -= length;
    }
    shift
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `target += scale * direction`.
fn add_scaled(target: &mut [f64; 3], scale: f64, direction: [f64; 3]) {
    for (t, d) in target.iter_mut().zip(direction) {
        *t += scale * d;
    }
}

impl SimBox {
    /// Minimum-image displacement in x and y; z is left untouched because the
    /// box is not periodic along the bilayer normal.
    pub fn minimum_image_xy(&self, displacement: [f64; 3]) -> [f64; 3] {
        [
            displacement[0] - minimum_image_shift(displacement[0], self.lx),
            displacement[1] - minimum_image_shift(displacement[1], self.ly),
            displacement[2],
        ]
    }
}

impl CgLipid {
    /// Minimum-image displacement `r_to - r_from` between two beads, or `None`
    /// if one of the indices does not address a bead.
    fn displacement(&self, from: usize, to: usize, box_: &SimBox) -> Option<[f64; 3]> {
        let (a, b) = (self.beads.get(from)?, self.beads.get(to)?);
        Some(box_.minimum_image_xy([
            b.position[0] - a.position[0],
            b.position[1] - a.position[1],
            b.position[2] - a.position[2],
        ]))
    }

    /// The two arm vectors `(r_i - r_j, r_k - r_j)` of an angle, if both are
    /// longer than [`MIN_ANGLE_ARM_NM`].
    fn angle_arms(
        &self,
        angle: (usize, usize, usize),
        box_: &SimBox,
    ) -> Option<([f64; 3], [f64; 3])> {
        let (i, j, k) = angle;
        let a = self.displacement(j, i, box_)?;
        let b = self.displacement(j, k, box_)?;
        (norm(a) > MIN_ANGLE_ARM_NM && norm(b) > MIN_ANGLE_ARM_NM).then_some((a, b))
    }

    /// Harmonic bond energy (kJ/mol) with minimum-image bond vectors.
    ///
    /// Same potential as [`CgLipid::bond_energy`], but valid for molecules
    /// whose beads sit on different sides of the periodic boundary.
    pub fn bond_energy_in_box(&self, box_: &SimBox) -> f64 {
        self.bonds
            .iter()
            .filter_map(|&(i, j)| self.displacement(i, j, box_))
            .map(|d| 0.5 * BOND_FORCE_CONSTANT * (norm(d) - BOND_R0_NM).powi(2))
            .sum()
    }

    /// Cosine angle energy (kJ/mol) with minimum-image arm vectors.
    ///
    /// Same potential as [`CgLipid::angle_energy`], but valid for molecules
    /// that straddle the periodic boundary.
    pub fn angle_energy_in_box(&self, box_: &SimBox) -> f64 {
        self.angles
            .iter()
            .filter_map(|&angle| self.angle_arms(angle, box_))
            .map(|(a, b)| {
                let cos_theta = dot(a, b) / (norm(a) * norm(b));
                0.5 * ANGLE_FORCE_CONSTANT * (cos_theta - ANGLE_COS_THETA0).powi(2)
            })
            .sum()
    }

    /// Total intramolecular energy (bond + angle, kJ/mol) with minimum-image
    /// vectors.
    pub fn intramolecular_energy_in_box(&self, box_: &SimBox) -> f64 {
        self.bond_energy_in_box(box_) + self.angle_energy_in_box(box_)
    }

    /// Forces (kJ/(mol nm)) on every bead from the bonds and angles of this
    /// lipid: the exact negative gradient of
    /// [`CgLipid::intramolecular_energy_in_box`].
    ///
    /// A bond with `d = r_j - r_i` and `r = |d|` pulls bead `i` along `+d` with
    /// `k (r - r₀)` (and bead `j` along `-d`).  An angle `i-j-k` with arm
    /// vectors `a = r_i - r_j`, `b = r_k - r_j` and `c = cos θ` exerts
    /// `F_i = -k_a (c - c₀) ∂c/∂r_i` with `∂c/∂r_i = b/(|a||b|) - c a/|a|²` (and
    /// the analogous expression for `k`); the central bead takes `-(F_i + F_k)`,
    /// so the net force and torque of the molecule vanish.
    pub fn intramolecular_forces(&self, box_: &SimBox) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0_f64; 3]; self.beads.len()];
        for &(i, j) in &self.bonds {
            let Some(d) = self.displacement(i, j, box_) else {
                continue;
            };
            let r = norm(d).max(1.0e-12);
            let pull = BOND_FORCE_CONSTANT * (r - BOND_R0_NM) / r;
            add_scaled(&mut forces[i], pull, d);
            add_scaled(&mut forces[j], -pull, d);
        }
        for &angle in &self.angles {
            let Some((a, b)) = self.angle_arms(angle, box_) else {
                continue;
            };
            let (i, j, k) = angle;
            let (ra, rb) = (norm(a), norm(b));
            let cos_theta = dot(a, b) / (ra * rb);
            let de_dcos = ANGLE_FORCE_CONSTANT * (cos_theta - ANGLE_COS_THETA0);
            let mut force_i = [0.0_f64; 3];
            let mut force_k = [0.0_f64; 3];
            add_scaled(&mut force_i, -de_dcos / (ra * rb), b);
            add_scaled(&mut force_i, de_dcos * cos_theta / (ra * ra), a);
            add_scaled(&mut force_k, -de_dcos / (ra * rb), a);
            add_scaled(&mut force_k, de_dcos * cos_theta / (rb * rb), b);
            add_scaled(&mut forces[i], 1.0, force_i);
            add_scaled(&mut forces[k], 1.0, force_k);
            add_scaled(&mut forces[j], -1.0, force_i);
            add_scaled(&mut forces[j], -1.0, force_k);
        }
        forces
    }

    /// Keeps the molecule whole and inside the primary cell.
    ///
    /// Every bond is first resolved to the periodic image of the partner bead
    /// closest to the bead it was reached from (walking the bond graph from the
    /// lowest-index bead of each connected fragment), so that raw coordinate
    /// differences inside the molecule are the physical ones again.  Each
    /// fragment is then translated by a whole number of box vectors so that its
    /// first bead lies in `[0, lx) × [0, ly)`.  Beads of a fragment may extend
    /// beyond the box; they are never split off the molecule.
    pub fn wrap_into_box(&mut self, box_: &SimBox) {
        let n = self.beads.len();
        let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(i, j) in &self.bonds {
            if i < n && j < n && i != j {
                neighbours[i].push(j);
                neighbours[j].push(i);
            }
        }
        let mut visited = vec![false; n];
        for start in 0..n {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            let mut fragment = vec![start];
            let mut stack = vec![start];
            while let Some(current) = stack.pop() {
                for &next in &neighbours[current] {
                    if visited[next] {
                        continue;
                    }
                    visited[next] = true;
                    let from = self.beads[current].position;
                    let to = self.beads[next].position;
                    let shift_x = minimum_image_shift(to[0] - from[0], box_.lx);
                    let shift_y = minimum_image_shift(to[1] - from[1], box_.ly);
                    self.beads[next].position[0] -= shift_x;
                    self.beads[next].position[1] -= shift_y;
                    fragment.push(next);
                    stack.push(next);
                }
            }
            let anchor = self.beads[start].position;
            let shift_x = wrap_shift(anchor[0], box_.lx);
            let shift_y = wrap_shift(anchor[1], box_.ly);
            if shift_x != 0.0 || shift_y != 0.0 {
                for &index in &fragment {
                    self.beads[index].position[0] -= shift_x;
                    self.beads[index].position[1] -= shift_y;
                }
            }
        }
    }
}
