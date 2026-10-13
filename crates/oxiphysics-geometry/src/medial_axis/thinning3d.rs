// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Curve skeleton by topology-preserving 3-D voxel thinning.
//!
//! The solid bounded by the samples is approximated by the union of the
//! samples' inner polar balls (the "power crust" solid of Amenta, Choi &
//! Kolluri 2001), rasterised on a regular grid. The voxel set is then thinned
//! with distance-ordered homotopic thinning (Pudney 1998): voxels are visited
//! in increasing distance to the boundary samples and deleted when they are a
//! *simple point* and not a medial anchor (a voxel holding the centre of an
//! inner polar ball). A voxel `p` is simple
//! (Malandain & Bertrand, "Fast characterization of 3D simple points", 1992)
//! iff
//!
//! * the foreground of `N26*(p)` has exactly one 26-connected component, and
//! * the background of `N18*(p)` has exactly one 6-connected component that is
//!   6-adjacent to `p`,
//!
//! so deleting it changes neither the number of objects, cavities nor tunnels.
//! Deletions are sequential, so every test sees the current voxel set and the
//! scheme is topology preserving. The anchors pin the skeleton to the medial
//! balls (instead of keeping every transient curve end point, which would
//! preserve spurs grown from rasterisation noise), and the grassfire order
//! keeps the connecting curves centred.
//! Its voxels become medial balls whose radius is the distance to the nearest
//! boundary sample, and are skeletonised with 26-adjacency.

use super::construct::{MedialAxisConfig, unresolved};
use super::kdtree::KdTree;
use super::poles::analyze;
use super::skeleton_graph::{SkelNode, skeletonize};
use super::{MedialAxis, dist3};
use crate::error::Result;

/// Dense voxel grid of booleans.
struct Grid {
    dims: [usize; 3],
    origin: [f64; 3],
    h: f64,
    cells: Vec<bool>,
}

impl Grid {
    fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.dims[1] + y) * self.dims[0] + x
    }

    fn center(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + (x as f64 + 0.5) * self.h,
            self.origin[1] + (y as f64 + 0.5) * self.h,
            self.origin[2] + (z as f64 + 0.5) * self.h,
        ]
    }

    /// Index of the voxel containing `p`, if it lies inside the grid.
    fn voxel_of(&self, p: [f64; 3]) -> Option<usize> {
        let mut c = [0usize; 3];
        for k in 0..3 {
            let f = ((p[k] - self.origin[k]) / self.h).floor();
            if !(f >= 0.0 && f < self.dims[k] as f64) {
                return None;
            }
            c[k] = f as usize;
        }
        Some(self.index(c[0], c[1], c[2]))
    }

    fn coords(&self, idx: usize) -> (usize, usize, usize) {
        let x = idx % self.dims[0];
        let y = (idx / self.dims[0]) % self.dims[1];
        let z = idx / (self.dims[0] * self.dims[1]);
        (x, y, z)
    }

    /// 3×3×3 neighbourhood of an interior voxel, indexed
    /// `(dz+1)*9 + (dy+1)*3 + (dx+1)`.
    fn neighbourhood(&self, x: usize, y: usize, z: usize) -> [bool; 27] {
        let mut nb = [false; 27];
        for dz in 0..3 {
            for dy in 0..3 {
                for dx in 0..3 {
                    let idx = self.index(x + dx - 1, y + dy - 1, z + dz - 1);
                    nb[dz * 9 + dy * 3 + dx] = self.cells[idx];
                }
            }
        }
        nb
    }
}

/// Offset triple of neighbourhood slot `k`.
fn offset(k: usize) -> [i32; 3] {
    [
        (k % 3) as i32 - 1,
        ((k / 3) % 3) as i32 - 1,
        (k / 9) as i32 - 1,
    ]
}

const CENTER: usize = 13;
/// Slots of the six face neighbours.
const FACES: [usize; 6] = [4, 10, 12, 14, 16, 22];

/// Number of 26-connected components of the foreground in `N26*(p)`.
fn fg_components_26(nb: &[bool; 27]) -> usize {
    let mut seen = [false; 27];
    let mut count = 0;
    for start in 0..27 {
        if start == CENTER || !nb[start] || seen[start] {
            continue;
        }
        count += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(k) = stack.pop() {
            let a = offset(k);
            for m in 0..27 {
                if m == CENTER || !nb[m] || seen[m] {
                    continue;
                }
                let b = offset(m);
                if (a[0] - b[0]).abs() <= 1 && (a[1] - b[1]).abs() <= 1 && (a[2] - b[2]).abs() <= 1
                {
                    seen[m] = true;
                    stack.push(m);
                }
            }
        }
    }
    count
}

/// Number of 6-connected background components of `N18*(p)` that are
/// 6-adjacent to `p`.
fn bg_components_6(nb: &[bool; 27]) -> usize {
    let in_n18 = |k: usize| {
        let o = offset(k);
        k != CENTER && o.iter().filter(|&&c| c != 0).count() <= 2
    };
    let mut seen = [false; 27];
    let mut count = 0;
    for &start in &FACES {
        if nb[start] || seen[start] {
            continue;
        }
        count += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(k) = stack.pop() {
            let a = offset(k);
            for m in 0..27 {
                if !in_n18(m) || nb[m] || seen[m] {
                    continue;
                }
                let b = offset(m);
                let manhattan = (a[0] - b[0]).abs() + (a[1] - b[1]).abs() + (a[2] - b[2]).abs();
                if manhattan == 1 {
                    seen[m] = true;
                    stack.push(m);
                }
            }
        }
    }
    count
}

/// Whether the centre voxel of `nb` is a 3-D simple point (26/6 topology).
fn is_simple(nb: &[bool; 27]) -> bool {
    fg_components_26(nb) == 1 && bg_components_6(nb) == 1
}

/// Rasterise the union of balls on a grid padded by two background layers.
fn rasterise(boundary: &[[f64; 3]], balls: &[([f64; 3], f64)], resolution: usize) -> Grid {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in boundary {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let longest = (0..3).map(|k| hi[k] - lo[k]).fold(0.0_f64, f64::max);
    let h = (longest / resolution as f64).max(f64::MIN_POSITIVE);
    let pad = 2usize;
    let mut dims = [0usize; 3];
    for k in 0..3 {
        dims[k] = ((hi[k] - lo[k]) / h).ceil() as usize + 1 + 2 * pad;
    }
    let origin = [
        lo[0] - pad as f64 * h,
        lo[1] - pad as f64 * h,
        lo[2] - pad as f64 * h,
    ];
    let mut grid = Grid {
        dims,
        origin,
        h,
        cells: vec![false; dims[0] * dims[1] * dims[2]],
    };
    let to_cell = |v: f64, k: usize| -> usize {
        let c = ((v - origin[k]) / h).floor();
        if c <= 0.0 {
            0
        } else {
            (c as usize).min(dims[k] - 1)
        }
    };
    for &(c, r) in balls {
        let (x0, x1) = (to_cell(c[0] - r, 0), to_cell(c[0] + r, 0));
        let (y0, y1) = (to_cell(c[1] - r, 1), to_cell(c[1] + r, 1));
        let (z0, z1) = (to_cell(c[2] - r, 2), to_cell(c[2] + r, 2));
        let r2 = r * r;
        for z in z0..=z1 {
            for y in y0..=y1 {
                for x in x0..=x1 {
                    let p = grid.center(x, y, z);
                    let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                    if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= r2 {
                        let idx = grid.index(x, y, z);
                        grid.cells[idx] = true;
                    }
                }
            }
        }
    }
    // Keep the two padding layers empty so every foreground voxel has a full
    // 3×3×3 neighbourhood inside the grid.
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let border = x == 0
                    || y == 0
                    || z == 0
                    || x + 1 == dims[0]
                    || y + 1 == dims[1]
                    || z + 1 == dims[2];
                if border {
                    let idx = grid.index(x, y, z);
                    grid.cells[idx] = false;
                }
            }
        }
    }
    grid
}

/// Distance-ordered homotopic thinning with medial anchors (in place).
///
/// Foreground voxels are visited in increasing order of `depth` (distance to
/// the boundary), so the object is peeled like a grassfire and the surviving
/// set is centred (Pudney, "Distance-ordered homotopic thinning", CVIU 1998).
/// A voxel is deleted when it is simple and not an *anchor* (a voxel holding
/// the centre of an inner polar ball, i.e. a medial-ball centre); sweeps
/// repeat until a full sweep deletes nothing. The result is the thinnest
/// topology-preserving voxel set that still contains every medial anchor.
fn thin(grid: &mut Grid, depth: &[f64], anchor: &[bool]) {
    let mut order: Vec<usize> = (0..grid.cells.len()).filter(|&i| grid.cells[i]).collect();
    order.sort_by(|&a, &b| depth[a].total_cmp(&depth[b]).then_with(|| a.cmp(&b)));
    loop {
        let mut changed = false;
        for &idx in &order {
            if !grid.cells[idx] || anchor[idx] {
                continue;
            }
            let (x, y, z) = grid.coords(idx);
            let nb = grid.neighbourhood(x, y, z);
            if !is_simple(&nb) {
                continue;
            }
            grid.cells[idx] = false;
            changed = true;
        }
        order.retain(|&i| grid.cells[i]);
        if !changed {
            break;
        }
    }
}

/// Curve skeleton by voxel thinning of the inner-polar-ball solid.
pub(super) fn thinning_skeleton(
    boundary: &[[f64; 3]],
    cfg: &MedialAxisConfig,
) -> Result<MedialAxis> {
    let pa = analyze(boundary, cfg.min_radius_factor)?;
    let mut pole_tets: Vec<usize> = pa.inner_pole.iter().flatten().copied().collect();
    pole_tets.sort_unstable();
    pole_tets.dedup();
    let balls: Vec<([f64; 3], f64)> = pole_tets
        .iter()
        .filter_map(|&ti| pa.spheres[ti].map(|s| (s.center, s.radius)))
        .collect();
    if balls.is_empty() {
        return Err(unresolved());
    }
    let mut grid = rasterise(boundary, &balls, cfg.thinning_resolution);

    // Depth of every foreground voxel: distance to the nearest sample.
    let used: Vec<usize> = (0..boundary.len())
        .filter(|&i| pa.normals[i].is_some())
        .collect();
    let tree = KdTree::build(boundary, &used);
    let mut depth = vec![0.0_f64; grid.cells.len()];
    for idx in 0..grid.cells.len() {
        if grid.cells[idx] {
            let (x, y, z) = grid.coords(idx);
            let c = grid.center(x, y, z);
            depth[idx] = tree
                .nearest_filtered(c, |_| true)
                .map_or(0.0, |(nearest, _)| dist3(c, boundary[nearest]));
        }
    }
    // Medial anchors: the voxels containing an inner-pole centre.
    let mut anchor = vec![false; grid.cells.len()];
    for &(c, _) in &balls {
        if let Some(idx) = grid.voxel_of(c)
            && grid.cells[idx]
        {
            anchor[idx] = true;
        }
    }
    thin(&mut grid, &depth, &anchor);

    let mut node_of = vec![usize::MAX; grid.cells.len()];
    let mut nodes = Vec::new();
    for idx in 0..grid.cells.len() {
        if !grid.cells[idx] {
            continue;
        }
        let (x, y, z) = grid.coords(idx);
        node_of[idx] = nodes.len();
        nodes.push(SkelNode {
            position: grid.center(x, y, z),
            radius: depth[idx],
            anchor: true,
        });
    }
    if nodes.is_empty() {
        return Err(unresolved());
    }
    let mut edges = Vec::new();
    for idx in 0..grid.cells.len() {
        let a = node_of[idx];
        if a == usize::MAX {
            continue;
        }
        let (x, y, z) = grid.coords(idx);
        // Forward half of the 26-neighbourhood (slots after the centre).
        for k in (CENTER + 1)..27 {
            let o = offset(k);
            let nx = (x as i64 + o[0] as i64) as usize;
            let ny = (y as i64 + o[1] as i64) as usize;
            let nz = (z as i64 + o[2] as i64) as usize;
            let b = node_of[grid.index(nx, ny, nz)];
            if b != usize::MAX {
                edges.push((a, b));
            }
        }
    }
    Ok(skeletonize(
        &nodes,
        &edges,
        cfg.skeleton_params(pa.diagonal),
    ))
}
