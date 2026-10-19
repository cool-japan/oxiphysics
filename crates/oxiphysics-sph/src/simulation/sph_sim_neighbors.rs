//! Cached neighbour lists of [`SphSim`](super::SphSim).
//!
//! Every pair sum of the array-based mini simulator (density, pressure,
//! viscosity and surface tension) only involves the particles inside the
//! cubic-spline support radius `2h`.  [`SphSimNeighbors`] finds those
//! neighbours once with a spatial hash in `O(N)` expected time and stores them
//! in the CSR layout of [`FlatNeighborList`] (`offsets` + `indices`), so the sums
//! cost `O(N k)` for `k` neighbours per particle instead of `O(N²)`.
//!
//! The lists remember the positions and kernel radius they were built for, which
//! makes it cheap and safe to decide whether they are still valid: positions are
//! public fields of `SphSim` and may be changed by the caller at any time.

use crate::neighbor::{FlatNeighborList, SpatialHash3D};

/// Neighbour lists of a particle cloud for the cubic-spline support `2h`.
///
/// Row `i` holds the indices `j ≠ i` with `|xᵢ − xⱼ| ≤ 2h` in ascending order.
/// A non-positive or non-finite kernel radius gives empty rows.
#[derive(Debug, Clone)]
pub struct SphSimNeighbors {
    /// Kernel radius `h` the lists were built for.
    kernel_radius: f64,
    /// Particle positions the lists were built for.
    positions: Vec<[f64; 3]>,
    /// The lists in CSR layout.
    list: FlatNeighborList,
}

impl SphSimNeighbors {
    /// Lists for no particles, which are never current for any kernel radius.
    pub fn empty() -> Self {
        Self {
            kernel_radius: f64::NAN,
            positions: Vec::new(),
            list: FlatNeighborList::from_neighbor_lists(&[]),
        }
    }

    /// Build the lists of `positions` for the support radius `2 · kernel_radius`.
    pub fn build(positions: &[[f64; 3]], kernel_radius: f64) -> Self {
        let radius = 2.0 * kernel_radius;
        let rows = if radius > 0.0 && radius.is_finite() {
            SpatialHash3D::all_neighbors(positions, radius)
        } else {
            vec![Vec::new(); positions.len()]
        };
        Self {
            kernel_radius,
            positions: positions.to_vec(),
            list: FlatNeighborList::from_neighbor_lists(&rows),
        }
    }

    /// Whether the lists were built for exactly these `positions` and this
    /// kernel radius, i.e. whether they can be used for them as they are.
    pub fn is_current(&self, positions: &[[f64; 3]], kernel_radius: f64) -> bool {
        self.kernel_radius == kernel_radius && self.positions.as_slice() == positions
    }

    /// The neighbour indices of particle `i`; empty if `i` is out of range.
    pub fn neighbors_of(&self, i: usize) -> &[usize] {
        let end = i
            .checked_add(1)
            .and_then(|next| self.list.offsets.get(next));
        match (self.list.offsets.get(i), end) {
            (Some(&start), Some(&end)) => self.list.flat_neighbors.get(start..end).unwrap_or(&[]),
            _ => &[],
        }
    }

    /// Number of particles the lists were built for.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether the lists were built for no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Total number of stored neighbour entries (each pair is stored twice).
    pub fn entry_count(&self) -> usize {
        self.list.total_pairs()
    }

    /// The kernel radius `h` the lists were built for.
    pub fn kernel_radius(&self) -> f64 {
        self.kernel_radius
    }

    /// The lists in CSR layout (`offsets`, `flat_neighbors`).
    pub fn as_flat(&self) -> &FlatNeighborList {
        &self.list
    }
}

impl Default for SphSimNeighbors {
    fn default() -> Self {
        Self::empty()
    }
}
