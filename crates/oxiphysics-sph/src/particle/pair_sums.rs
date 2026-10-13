//! Neighbour-list based SPH pair sums for [`SphParticleSet`].
//!
//! The weakly-compressible SPH pipeline needs three pair sums per step: the
//! density summation, the symmetric pressure acceleration and the Morris
//! viscosity acceleration.  Every one of them only involves the particles
//! inside the cubic-spline support radius `2h`, so they are evaluated here from
//! neighbour lists (`neighbors[i]` = indices `j ≠ i` with `|xᵢ − xⱼ| ≤ 2h`) in
//! `O(N · k)` for `k` neighbours per particle instead of looping over all `N²`
//! ordered pairs.
//!
//! Build the lists once per step with
//! [`SphParticleSet::support_neighbor_lists`] and hand the same lists to every
//! `*_with_neighbors` method; the plain methods
//! ([`SphParticleSet::compute_sph_densities`],
//! [`SphParticleSet::accumulate_pressure_acceleration`],
//! [`SphParticleSet::accumulate_morris_viscosity`]) build a fresh list on every
//! call and are the convenient entry points when only one sum is needed.
//!
//! Each fluid particle `i` accumulates its own contributions from *all* of its
//! neighbours exactly once (boundary particles act as sources but receive no
//! acceleration), so no pair is ever counted twice.

use super::functions::{cubic_spline_gradient, cubic_spline_kernel};
use super::types::SphParticleSet;
use crate::neighbor::SpatialHash3D;

/// Difference vector `a − b` and its Euclidean length.
#[inline]
fn separation(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], f64) {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    (d, r)
}

impl SphParticleSet {
    /// Neighbour lists for the cubic-spline support radius `2h`.
    ///
    /// Returns one sorted list per particle holding the indices `j ≠ i` with
    /// `|xᵢ − xⱼ| ≤ 2h` (boundary particles included).  The lists are built with
    /// a spatial hash in `O(N)` expected time and are deterministic.  A
    /// non-positive or non-finite `h` yields empty lists.
    pub fn support_neighbor_lists(&self, h: f64) -> Vec<Vec<usize>> {
        let radius = 2.0 * h;
        if !(radius > 0.0 && radius.is_finite()) {
            return vec![Vec::new(); self.len()];
        }
        SpatialHash3D::all_neighbors(&self.positions, radius)
    }

    /// SPH density summation from precomputed neighbour lists.
    ///
    /// `ρᵢ = mᵢ W(0, h) + Σ_{j ∈ N(i)} mⱼ W(|xᵢ − xⱼ|, h)` with the cubic-spline
    /// kernel [`cubic_spline_kernel`].  A list for particle `i` that is missing
    /// from `neighbors`, an index `j` outside the particle set and an `i`
    /// listed among its own neighbours are ignored.
    pub fn compute_sph_densities_with_neighbors(&mut self, h: f64, neighbors: &[Vec<usize>]) {
        let n = self.len();
        let self_weight = cubic_spline_kernel(0.0, h);
        for i in 0..n {
            let pi = self.positions[i];
            let mut density = self.masses[i] * self_weight;
            if let Some(list) = neighbors.get(i) {
                for &j in list {
                    if j == i || j >= n {
                        continue;
                    }
                    let (_, r) = separation(pi, self.positions[j]);
                    density += self.masses[j] * cubic_spline_kernel(r, h);
                }
            }
            self.densities[i] = density;
        }
    }

    /// Symmetric SPH pressure acceleration from precomputed neighbour lists.
    ///
    /// Every fluid particle `i` accumulates
    /// `aᵢ += −Σ_{j ∈ N(i)} mⱼ (Pᵢ/ρᵢ² + Pⱼ/ρⱼ²) ∇ᵢWᵢⱼ`
    /// exactly once for each neighbour `j` (boundary neighbours included);
    /// boundary particles themselves receive no acceleration.  Because the
    /// bracket is symmetric in `i` and `j`, a pair of fluid particles receives
    /// equal and opposite momentum changes `mᵢ aᵢ = −mⱼ aⱼ`.
    ///
    /// Accumulates into `self.accelerations` (does *not* clear first).
    pub fn accumulate_pressure_acceleration_with_neighbors(
        &mut self,
        h: f64,
        neighbors: &[Vec<usize>],
    ) {
        let n = self.len();
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let Some(list) = neighbors.get(i) else {
                continue;
            };
            let pi = self.positions[i];
            let rho_i = self.densities[i].max(1e-20);
            let p_i = self.pressures[i];
            let mut acc = [0.0_f64; 3];
            for &j in list {
                if j == i || j >= n {
                    continue;
                }
                let (d, r) = separation(pi, self.positions[j]);
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let rho_j = self.densities[j].max(1e-20);
                let p_j = self.pressures[j];
                let m_j = self.masses[j];
                let dw = cubic_spline_gradient(r, h);
                let factor = -m_j * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw / r;
                acc[0] += factor * d[0];
                acc[1] += factor * d[1];
                acc[2] += factor * d[2];
            }
            self.accelerations[i][0] += acc[0];
            self.accelerations[i][1] += acc[1];
            self.accelerations[i][2] += acc[2];
        }
    }

    /// Morris et al. (1997) viscous acceleration from precomputed neighbour
    /// lists.
    ///
    /// For every fluid particle `i`
    /// `aᵢ += Σ_{j ∈ N(i)} mⱼ · 2 μ / (ρᵢ ρⱼ) · (dW/dr)/r · (vᵢ − vⱼ)`,
    /// the SPH form of `(μ/ρ) ∇²v` with the *dynamic* viscosity `μ` \[Pa·s\]
    /// (`dW/dr ≤ 0`, so the term pulls `vᵢ` towards the neighbouring
    /// velocities).  For a uniform density `ρ₀` this reduces to the kinematic
    /// viscosity `ν = μ/ρ₀`.
    ///
    /// Accumulates into `self.accelerations`; a vanishing `mu` is a no-op.
    pub fn accumulate_morris_viscosity_with_neighbors(
        &mut self,
        h: f64,
        mu: f64,
        neighbors: &[Vec<usize>],
    ) {
        if mu.abs() < 1e-30 {
            return;
        }
        let n = self.len();
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let Some(list) = neighbors.get(i) else {
                continue;
            };
            let pi = self.positions[i];
            let vi = self.velocities[i];
            let rho_i = self.densities[i].max(1e-20);
            let mut acc = [0.0_f64; 3];
            for &j in list {
                if j == i || j >= n {
                    continue;
                }
                let (_, r) = separation(pi, self.positions[j]);
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let vj = self.velocities[j];
                let rho_j = self.densities[j].max(1e-20);
                let dw = cubic_spline_gradient(r, h);
                let coeff = 2.0 * mu / (rho_i * rho_j) * dw / r;
                let factor = self.masses[j] * coeff;
                acc[0] += factor * (vi[0] - vj[0]);
                acc[1] += factor * (vi[1] - vj[1]);
                acc[2] += factor * (vi[2] - vj[2]);
            }
            self.accelerations[i][0] += acc[0];
            self.accelerations[i][1] += acc[1];
            self.accelerations[i][2] += acc[2];
        }
    }
}
