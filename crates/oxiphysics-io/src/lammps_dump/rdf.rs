// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Radial distribution function of a dump frame.

use super::LammpsDumpFrame;
use super::layout::DumpBox;

/// Compute the radial distribution function g(r) from a single frame.
///
/// Returns `(r_bins, g_of_r)` where `r_bins[i]` is the centre of bin `i`
/// and `g_of_r[i]` is g(r) at that bin.
///
/// Centre atoms are those of `type_a` and neighbours those of `type_b`; `0`
/// means "every type" for either argument.  The pair `(i, j)` with `i != j` is
/// counted once for each centre/neighbour assignment, and g(r) is normalised by
/// the number of such ordered pairs per volume, so an ideal gas gives g = 1 for
/// any selection (including the case where `i` itself is also in the
/// neighbour set):
///
/// ```text
/// g(r_i) = hist_i / (rho_pairs * V_shell_i),   rho_pairs = (N_a * N_b - N_ab) / V,
/// V_shell_i = 4 pi / 3 * ((r_i + dr)^3 - r_i^3),   r_i = i * dr
/// ```
///
/// where `N_ab` is the number of atoms that are both centres and neighbours and
/// `r_i` is the inner edge of bin `i`.
///
/// The box is periodic: distances use the minimum image of the (possibly
/// triclinic) cell.  The minimum image is only unique below half the smallest
/// perpendicular cell width, so `r_max` is clamped to that value; the bins
/// then span `[0, r_max]` of the clamped radius and `r_bins` shows the radii
/// actually used.  Atoms with a non-finite coordinate are ignored.
///
/// A frame without a valid box, a non-positive or NaN `r_max`, or an empty
/// selection yields `g = 0` in every bin (an infinite `r_max` simply selects the
/// largest usable radius).  `n_bins == 0` yields empty vectors.
pub fn radial_distribution_function(
    frame: &LammpsDumpFrame,
    type_a: u32,
    type_b: u32,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    if n_bins == 0 {
        return (Vec::new(), Vec::new());
    }
    let dump_box = DumpBox::from_frame(frame);
    let usable_r_max = if dump_box.is_valid() && r_max > 0.0 {
        r_max.min(dump_box.max_periodic_distance())
    } else {
        0.0
    };
    let dr = usable_r_max / n_bins as f64;
    let r_bins: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
    let zeros = vec![0.0; n_bins];
    if usable_r_max <= 0.0 {
        return (r_bins, zeros);
    }

    let in_a = |t: u32| type_a == 0 || t == type_a;
    let in_b = |t: u32| type_b == 0 || t == type_b;
    let atoms: Vec<RdfAtom> = frame
        .types
        .iter()
        .zip(&frame.positions)
        .filter(|(_, p)| p.iter().all(|c| c.is_finite()))
        .map(|(&t, &position)| RdfAtom {
            position,
            center: in_a(t),
            neighbour: in_b(t),
        })
        .collect();
    let n_a = atoms.iter().filter(|a| a.center).count() as f64;
    let n_b = atoms.iter().filter(|a| a.neighbour).count() as f64;
    let n_ab = atoms.iter().filter(|a| a.center && a.neighbour).count() as f64;
    let n_pairs = n_a * n_b - n_ab;
    let volume = dump_box.volume();
    if n_pairs <= 0.0 || !volume.is_finite() || volume <= 0.0 {
        return (r_bins, zeros);
    }

    let hist = pair_histogram(&atoms, &dump_box, usable_r_max, n_bins);
    let pair_density = n_pairs / volume;
    let g_of_r = hist
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let inner = i as f64 * dr;
            let outer = inner + dr;
            let shell_volume = 4.0 * std::f64::consts::PI / 3.0 * (outer.powi(3) - inner.powi(3));
            count as f64 / (pair_density * shell_volume)
        })
        .collect();
    (r_bins, g_of_r)
}

/// One atom as seen by the g(r) calculation.
struct RdfAtom {
    position: [f64; 3],
    /// Belongs to the centre selection.
    center: bool,
    /// Belongs to the neighbour selection.
    neighbour: bool,
}

/// Count ordered centre/neighbour pairs by minimum-image distance.
///
/// Atoms are binned into a periodic grid of cells at least `r_max` thick (in the
/// direction perpendicular to each cell face), so every pair closer than `r_max`
/// lies in the same or an adjacent cell and the work is linear in the number of
/// atoms.
fn pair_histogram(atoms: &[RdfAtom], dump_box: &DumpBox, r_max: f64, n_bins: usize) -> Vec<u64> {
    let grid = CellGrid::build(atoms, dump_box, r_max);
    let r_max_sq = r_max * r_max;
    let bins_per_length = n_bins as f64 / r_max;
    let mut hist = vec![0u64; n_bins];
    for (i, center) in atoms.iter().enumerate() {
        if !center.center {
            continue;
        }
        let home = grid.cell_coords(i);
        let (xs, nx) = axis_neighbours(home[0], grid.dims[0]);
        let (ys, ny) = axis_neighbours(home[1], grid.dims[1]);
        let (zs, nz) = axis_neighbours(home[2], grid.dims[2]);
        for &cx in &xs[..nx] {
            for &cy in &ys[..ny] {
                for &cz in &zs[..nz] {
                    for &j in grid.members(grid.cell_id([cx, cy, cz])) {
                        let j = j as usize;
                        let other = &atoms[j];
                        if j == i || !other.neighbour {
                            continue;
                        }
                        let d = dump_box.minimum_image([
                            other.position[0] - center.position[0],
                            other.position[1] - center.position[1],
                            other.position[2] - center.position[2],
                        ]);
                        let r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if r_sq < r_max_sq {
                            let bin = (r_sq.sqrt() * bins_per_length) as usize;
                            if let Some(slot) = hist.get_mut(bin) {
                                *slot += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    hist
}

/// Distinct cell indices among `c - 1`, `c`, `c + 1` on a periodic axis of `n`
/// cells (fewer than three when the axis has fewer than three cells).
fn axis_neighbours(c: usize, n: usize) -> ([usize; 3], usize) {
    let candidates = [(c + n - 1) % n, c, (c + 1) % n];
    let mut out = [0usize; 3];
    let mut len = 0;
    for v in candidates {
        if !out[..len].contains(&v) {
            out[len] = v;
            len += 1;
        }
    }
    (out, len)
}

/// Periodic cell list over fractional coordinates, stored as compressed rows.
struct CellGrid {
    dims: [usize; 3],
    /// Cell coordinates of every atom.
    coords: Vec<[usize; 3]>,
    /// `start[c]..start[c + 1]` indexes `members` for cell `c`.
    start: Vec<usize>,
    members: Vec<u32>,
}

impl CellGrid {
    fn build(atoms: &[RdfAtom], dump_box: &DumpBox, r_max: f64) -> Self {
        let widths = dump_box.perpendicular_widths();
        // At most a few cells per atom, so a tiny cutoff cannot exhaust memory.
        let max_cells = (4 * atoms.len()).max(64);
        let mut dims = [1usize; 3];
        for (dim, width) in dims.iter_mut().zip(widths) {
            *dim = ((width / r_max).floor() as usize).max(1);
        }
        while dims[0].saturating_mul(dims[1]).saturating_mul(dims[2]) > max_cells {
            let widest = (0..3).max_by_key(|&k| dims[k]).unwrap_or(0);
            dims[widest] = dims[widest].div_ceil(2);
        }
        let origin = dump_box.origin();
        let coords: Vec<[usize; 3]> = atoms
            .iter()
            .map(|a| {
                let s = dump_box.to_fractional_displacement([
                    a.position[0] - origin[0],
                    a.position[1] - origin[1],
                    a.position[2] - origin[2],
                ]);
                let mut cell = [0usize; 3];
                for ((slot, &fraction), &n) in cell.iter_mut().zip(&s).zip(&dims) {
                    let wrapped = fraction - fraction.floor();
                    *slot = ((wrapped * n as f64) as usize).min(n - 1);
                }
                cell
            })
            .collect();
        let n_cells = dims[0] * dims[1] * dims[2];
        let mut start = vec![0usize; n_cells + 1];
        for c in &coords {
            start[(c[0] * dims[1] + c[1]) * dims[2] + c[2] + 1] += 1;
        }
        let mut running = 0usize;
        for slot in start.iter_mut() {
            running += *slot;
            *slot = running;
        }
        let mut fill = start.clone();
        let mut members = vec![0u32; atoms.len()];
        for (i, c) in coords.iter().enumerate() {
            let cell = (c[0] * dims[1] + c[1]) * dims[2] + c[2];
            members[fill[cell]] = i as u32;
            fill[cell] += 1;
        }
        Self {
            dims,
            coords,
            start,
            members,
        }
    }

    fn cell_coords(&self, atom: usize) -> [usize; 3] {
        self.coords[atom]
    }

    fn cell_id(&self, c: [usize; 3]) -> usize {
        (c[0] * self.dims[1] + c[1]) * self.dims[2] + c[2]
    }

    fn members(&self, cell: usize) -> &[u32] {
        &self.members[self.start[cell]..self.start[cell + 1]]
    }
}
