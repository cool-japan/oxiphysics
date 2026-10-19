//! Kernel-gradient and neighbour-search helpers for the pressure-solver types.
//!
//! The pair sums of [`IisphMatrixBuilder`](super::IisphMatrixBuilder),
//! [`SphDensityIntegrator`](super::SphDensityIntegrator),
//! [`PcisphCorrector`](super::PcisphCorrector),
//! [`DivergenceFreeProjector`](super::DivergenceFreeProjector) and
//! [`IsphPressure`](super::IsphPressure) all use the exact gradient of the
//! cubic-spline kernel ([`CubicSplineKernel`]) with support radius `2h`.  The
//! gradient is `∇ᵢW(|xᵢ − xⱼ|, h) = W'(r) (xᵢ − xⱼ) / r` with `W'(r) ≤ 0`, so it
//! points from particle `i` towards particle `j`, vanishes for `r ≥ 2h` and is
//! finite (indeed zero) at `r = 0`.
//!
//! Whole-field operations find the neighbours of every particle once, with a
//! spatial hash in `O(N)` expected time, and store them in the CSR layout of
//! [`FlatNeighborList`]; only single-particle queries scan all particles.

use crate::kernel::{CubicSplineKernel, SphKernel};
use crate::neighbor::{FlatNeighborList, SpatialHash3D};

/// Pair distances below this value (m) are treated as coincident particles
/// whose direction is undefined, so the pair is skipped.
pub(super) const MIN_PAIR_DISTANCE: f64 = 1e-14;

/// Separation vector `a − b` and its Euclidean length.
#[inline]
pub(super) fn separation(a: [f64; 3], b: [f64; 3]) -> ([f64; 3], f64) {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    (d, r)
}

/// Cubic-spline kernel gradient `∇ᵢW(r, h) = W'(r) r_ij / r` for the separation
/// `r_ij = xᵢ − xⱼ` of length `r = |r_ij|`.
///
/// The result is the zero vector for coincident particles (`r = 0`), for
/// `r ≥ 2h` and for a non-positive or non-finite `h`.
#[inline]
pub(super) fn kernel_gradient(r_ij: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if r.is_nan() || r <= 0.0 {
        return [0.0; 3];
    }
    let scale = CubicSplineKernel.grad_w(r, h) / r;
    [scale * r_ij[0], scale * r_ij[1], scale * r_ij[2]]
}

/// `(v_j − v_i) · ∇ᵢWᵢⱼ` for the pair of particles at `x_i`, `x_j` with
/// velocities `v_i`, `v_j`.
///
/// Zero for coincident particles and outside the support `r ≥ 2h`.  Weighting
/// it with the neighbour volume `m_j / ρ_j` and summing over `j` gives the SPH
/// velocity divergence `∇·v`; weighting it with `−m_j` gives the continuity
/// rate `Dρ/Dt`.
#[inline]
pub(super) fn relative_velocity_dot_gradient(
    x_i: [f64; 3],
    v_i: [f64; 3],
    x_j: [f64; 3],
    v_j: [f64; 3],
    h: f64,
) -> f64 {
    let (r_ij, r) = separation(x_i, x_j);
    if r < MIN_PAIR_DISTANCE || r >= 2.0 * h {
        return 0.0;
    }
    let grad_w = kernel_gradient(r_ij, r, h);
    (v_j[0] - v_i[0]) * grad_w[0] + (v_j[1] - v_i[1]) * grad_w[1] + (v_j[2] - v_i[2]) * grad_w[2]
}

/// Neighbour lists in CSR layout for the cubic-spline support radius `2h`.
///
/// Row `i` holds the indices `j ≠ i` with `|xᵢ − xⱼ| ≤ 2h` in ascending order.
/// A non-positive or non-finite `h` gives empty rows.
pub(super) fn support_neighbors(positions: &[[f64; 3]], h: f64) -> FlatNeighborList {
    let radius = 2.0 * h;
    if radius > 0.0 && radius.is_finite() {
        FlatNeighborList::from_neighbor_lists(&SpatialHash3D::all_neighbors(positions, radius))
    } else {
        FlatNeighborList::from_neighbor_lists(&vec![Vec::<usize>::new(); positions.len()])
    }
}

/// The neighbour indices of particle `i`; empty if `i` is outside the list.
#[inline]
pub(super) fn neighbor_row(list: &FlatNeighborList, i: usize) -> &[usize] {
    let end = i.checked_add(1).and_then(|next| list.offsets.get(next));
    match (list.offsets.get(i), end) {
        (Some(&start), Some(&end)) => list.flat_neighbors.get(start..end).unwrap_or(&[]),
        _ => &[],
    }
}
