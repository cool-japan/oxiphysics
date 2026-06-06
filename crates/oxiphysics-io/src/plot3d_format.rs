// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PLOT3D format for CFD grid and solution files.
//!
//! Supports multi-block structured grids in ASCII PLOT3D `.xyz` and `.q`
//! formats as commonly used in CFD solvers.
//!
//! # File layout
//!
//! ## Grid file (.xyz)
//! ```text
//! `nblocks`
//! `ni1` `nj1` `nk1`  [`ni2` `nj2` `nk2` …]
//! <x values block 1> <y values block 1> <z values block 1>
//! …
//! ```
//!
//! ## Solution file (.q)
//! ```text
//! `nblocks`
//! `ni1` `nj1` `nk1`  [`ni2` `nj2` `nk2` …]
//! `mach` `alpha` `rey` `time`   (per block)
//! `rho` <rho*u> <rho*v> <rho*w> `e`  (per point, per block)
//! ```

use std::fmt::Write as FmtWrite;
use std::str::SplitAsciiWhitespace;

use crate::Error as IoError;

// ── Helper ───────────────────────────────────────────────────────────────────

/// Parse the next f64 token from an iterator.
fn next_f64<'a>(iter: &mut SplitAsciiWhitespace<'a>) -> Result<f64, IoError> {
    iter.next()
        .ok_or_else(|| IoError::Parse("unexpected end of PLOT3D data".into()))?
        .parse::<f64>()
        .map_err(|e| IoError::Parse(e.to_string()))
}

/// Parse the next usize token from an iterator.
fn next_usize<'a>(iter: &mut SplitAsciiWhitespace<'a>) -> Result<usize, IoError> {
    iter.next()
        .ok_or_else(|| IoError::Parse("unexpected end of PLOT3D data".into()))?
        .parse::<usize>()
        .map_err(|e| IoError::Parse(e.to_string()))
}

// ── Grid ─────────────────────────────────────────────────────────────────────

/// A single structured block in a PLOT3D grid.
///
/// Coordinates are stored in row-major order with the fastest index being
/// `i` (x-direction), then `j` (y-direction), then `k` (z-direction).
#[derive(Debug, Clone, PartialEq)]
pub struct Plot3dBlock {
    /// Number of grid points in the i-direction.
    pub ni: usize,
    /// Number of grid points in the j-direction.
    pub nj: usize,
    /// Number of grid points in the k-direction.
    pub nk: usize,
    /// x-coordinates (`ni × nj × nk` values).
    pub x: Vec<f64>,
    /// y-coordinates (`ni × nj × nk` values).
    pub y: Vec<f64>,
    /// z-coordinates (`ni × nj × nk` values).
    pub z: Vec<f64>,
}

impl Plot3dBlock {
    /// Total number of grid points in this block.
    pub fn n_points(&self) -> usize {
        self.ni * self.nj * self.nk
    }
}

/// A multi-block PLOT3D structured grid.
#[derive(Debug, Clone)]
pub struct Plot3dGrid {
    /// Ordered list of structured blocks.
    pub blocks: Vec<Plot3dBlock>,
}

impl Plot3dGrid {
    /// Create a new empty grid.
    pub fn new() -> Self {
        Self { blocks: Vec::new() }
    }

    /// Number of blocks.
    pub fn n_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Dimensions `(ni, nj, nk)` of block `b`.
    ///
    /// Returns `None` if `b` is out of range.
    pub fn block_dimensions(&self, b: usize) -> Option<(usize, usize, usize)> {
        self.blocks.get(b).map(|blk| (blk.ni, blk.nj, blk.nk))
    }
}

impl Default for Plot3dGrid {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse an ASCII PLOT3D `.xyz` file from a string.
///
/// # Errors
/// Returns `Error::Parse` if the input is malformed.
pub fn read_plot3d_grid(src: &str) -> Result<Plot3dGrid, IoError> {
    let mut iter = src.split_ascii_whitespace();
    let nblocks = next_usize(&mut iter)?;
    let mut dims: Vec<(usize, usize, usize)> = Vec::with_capacity(nblocks);
    for _ in 0..nblocks {
        let ni = next_usize(&mut iter)?;
        let nj = next_usize(&mut iter)?;
        let nk = next_usize(&mut iter)?;
        dims.push((ni, nj, nk));
    }
    let mut blocks = Vec::with_capacity(nblocks);
    for (ni, nj, nk) in dims {
        let n = ni * nj * nk;
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        let mut z = Vec::with_capacity(n);
        for _ in 0..n {
            x.push(next_f64(&mut iter)?);
        }
        for _ in 0..n {
            y.push(next_f64(&mut iter)?);
        }
        for _ in 0..n {
            z.push(next_f64(&mut iter)?);
        }
        blocks.push(Plot3dBlock {
            ni,
            nj,
            nk,
            x,
            y,
            z,
        });
    }
    Ok(Plot3dGrid { blocks })
}

/// Write a PLOT3D grid to an ASCII string.
pub fn write_plot3d_grid(grid: &Plot3dGrid) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "{}", grid.n_blocks());
    for blk in &grid.blocks {
        let _ = writeln!(out, "{} {} {}", blk.ni, blk.nj, blk.nk);
    }
    for blk in &grid.blocks {
        let n = blk.n_points();
        for k in 0..n {
            let _ = writeln!(out, "{:.15e}", blk.x[k]);
        }
        for k in 0..n {
            let _ = writeln!(out, "{:.15e}", blk.y[k]);
        }
        for k in 0..n {
            let _ = writeln!(out, "{:.15e}", blk.z[k]);
        }
    }
    out
}

// ── Solution ─────────────────────────────────────────────────────────────────

/// Free-stream reference parameters stored per block in a `.q` file.
#[derive(Debug, Clone, PartialEq)]
pub struct Plot3dFreestream {
    /// Mach number.
    pub mach: f64,
    /// Angle of attack (degrees).
    pub alpha: f64,
    /// Reynolds number.
    pub reynolds: f64,
    /// Non-dimensional time.
    pub time: f64,
}

/// Solution variables at a single grid point.
///
/// The five conservative variables used in the PLOT3D `.q` format are
/// stored: density, three momentum components, and total energy.
#[derive(Debug, Clone, PartialEq)]
pub struct Plot3dQPoint {
    /// Density ρ.
    pub rho: f64,
    /// x-momentum ρ·u.
    pub rho_u: f64,
    /// y-momentum ρ·v.
    pub rho_v: f64,
    /// z-momentum ρ·w.
    pub rho_w: f64,
    /// Total energy per unit volume e.
    pub e: f64,
}

/// A single solution block in a PLOT3D `.q` file.
#[derive(Debug, Clone)]
pub struct Plot3dSolutionBlock {
    /// Grid dimensions matching the corresponding [`Plot3dBlock`].
    pub ni: usize,
    /// Grid dimensions matching the corresponding [`Plot3dBlock`].
    pub nj: usize,
    /// Grid dimensions matching the corresponding [`Plot3dBlock`].
    pub nk: usize,
    /// Free-stream reference parameters.
    pub freestream: Plot3dFreestream,
    /// Solution data (one entry per grid point).
    pub q: Vec<Plot3dQPoint>,
}

impl Plot3dSolutionBlock {
    /// Total number of solution points.
    pub fn n_points(&self) -> usize {
        self.ni * self.nj * self.nk
    }
}

/// A multi-block PLOT3D solution.
#[derive(Debug, Clone)]
pub struct Plot3dSolution {
    /// Ordered list of solution blocks.
    pub blocks: Vec<Plot3dSolutionBlock>,
}

impl Plot3dSolution {
    /// Number of blocks in this solution.
    pub fn n_blocks(&self) -> usize {
        self.blocks.len()
    }
}

/// Parse an ASCII PLOT3D `.q` solution file from a string.
///
/// # Errors
/// Returns `Error::Parse` if the input is malformed.
pub fn read_plot3d_solution(src: &str) -> Result<Plot3dSolution, IoError> {
    let mut iter = src.split_ascii_whitespace();
    let nblocks = next_usize(&mut iter)?;
    let mut dims: Vec<(usize, usize, usize)> = Vec::with_capacity(nblocks);
    for _ in 0..nblocks {
        let ni = next_usize(&mut iter)?;
        let nj = next_usize(&mut iter)?;
        let nk = next_usize(&mut iter)?;
        dims.push((ni, nj, nk));
    }
    let mut blocks = Vec::with_capacity(nblocks);
    for (ni, nj, nk) in dims {
        let mach = next_f64(&mut iter)?;
        let alpha = next_f64(&mut iter)?;
        let reynolds = next_f64(&mut iter)?;
        let time = next_f64(&mut iter)?;
        let n = ni * nj * nk;
        let mut q = Vec::with_capacity(n);
        for _ in 0..n {
            let rho = next_f64(&mut iter)?;
            let rho_u = next_f64(&mut iter)?;
            let rho_v = next_f64(&mut iter)?;
            let rho_w = next_f64(&mut iter)?;
            let e = next_f64(&mut iter)?;
            q.push(Plot3dQPoint {
                rho,
                rho_u,
                rho_v,
                rho_w,
                e,
            });
        }
        blocks.push(Plot3dSolutionBlock {
            ni,
            nj,
            nk,
            freestream: Plot3dFreestream {
                mach,
                alpha,
                reynolds,
                time,
            },
            q,
        });
    }
    Ok(Plot3dSolution { blocks })
}

/// Write a PLOT3D solution to an ASCII string.
pub fn write_plot3d_solution(sol: &Plot3dSolution) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "{}", sol.n_blocks());
    for blk in &sol.blocks {
        let _ = writeln!(out, "{} {} {}", blk.ni, blk.nj, blk.nk);
    }
    for blk in &sol.blocks {
        let fs = &blk.freestream;
        let _ = writeln!(
            out,
            "{:.15e} {:.15e} {:.15e} {:.15e}",
            fs.mach, fs.alpha, fs.reynolds, fs.time
        );
        for pt in &blk.q {
            let _ = writeln!(
                out,
                "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}",
                pt.rho, pt.rho_u, pt.rho_v, pt.rho_w, pt.e
            );
        }
    }
    out
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_single_block_grid() -> Plot3dGrid {
        let blk = Plot3dBlock {
            ni: 2,
            nj: 2,
            nk: 1,
            x: vec![0.0, 1.0, 0.0, 1.0],
            y: vec![0.0, 0.0, 1.0, 1.0],
            z: vec![0.0, 0.0, 0.0, 0.0],
        };
        Plot3dGrid { blocks: vec![blk] }
    }

    fn make_single_block_solution() -> Plot3dSolution {
        let fs = Plot3dFreestream {
            mach: 0.5,
            alpha: 0.0,
            reynolds: 1e6,
            time: 0.0,
        };
        let pts: Vec<Plot3dQPoint> = (0..4)
            .map(|i| Plot3dQPoint {
                rho: 1.0 + i as f64 * 0.1,
                rho_u: 0.5,
                rho_v: 0.0,
                rho_w: 0.0,
                e: 2.5,
            })
            .collect();
        Plot3dSolution {
            blocks: vec![Plot3dSolutionBlock {
                ni: 2,
                nj: 2,
                nk: 1,
                freestream: fs,
                q: pts,
            }],
        }
    }

    #[test]
    fn test_grid_n_blocks() {
        let g = make_single_block_grid();
        assert_eq!(g.n_blocks(), 1);
    }

    #[test]
    fn test_grid_n_points() {
        let g = make_single_block_grid();
        assert_eq!(g.blocks[0].n_points(), 4);
    }

    #[test]
    fn test_block_dimensions() {
        let g = make_single_block_grid();
        assert_eq!(g.block_dimensions(0), Some((2, 2, 1)));
    }

    #[test]
    fn test_block_dimensions_out_of_range() {
        let g = make_single_block_grid();
        assert_eq!(g.block_dimensions(5), None);
    }

    #[test]
    fn test_write_read_grid_roundtrip() {
        let original = make_single_block_grid();
        let s = write_plot3d_grid(&original);
        let parsed = read_plot3d_grid(&s).expect("parse should succeed");
        assert_eq!(parsed.n_blocks(), original.n_blocks());
        let orig_blk = &original.blocks[0];
        let pars_blk = &parsed.blocks[0];
        assert_eq!(pars_blk.ni, orig_blk.ni);
        assert_eq!(pars_blk.nj, orig_blk.nj);
        assert_eq!(pars_blk.nk, orig_blk.nk);
        for i in 0..orig_blk.n_points() {
            assert!((orig_blk.x[i] - pars_blk.x[i]).abs() < 1e-10);
            assert!((orig_blk.y[i] - pars_blk.y[i]).abs() < 1e-10);
            assert!((orig_blk.z[i] - pars_blk.z[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_read_grid_error_on_empty() {
        assert!(read_plot3d_grid("").is_err());
    }

    #[test]
    fn test_read_grid_error_on_truncated() {
        // Only block count with no dimensions
        assert!(read_plot3d_grid("1").is_err());
    }

    #[test]
    fn test_write_grid_contains_nblocks() {
        let g = make_single_block_grid();
        let s = write_plot3d_grid(&g);
        assert!(s.starts_with("1\n"));
    }

    #[test]
    fn test_solution_n_blocks() {
        let sol = make_single_block_solution();
        assert_eq!(sol.n_blocks(), 1);
    }

    #[test]
    fn test_solution_n_points() {
        let sol = make_single_block_solution();
        assert_eq!(sol.blocks[0].n_points(), 4);
    }

    #[test]
    fn test_write_read_solution_roundtrip() {
        let original = make_single_block_solution();
        let s = write_plot3d_solution(&original);
        let parsed = read_plot3d_solution(&s).expect("parse should succeed");
        assert_eq!(parsed.n_blocks(), original.n_blocks());
        let ob = &original.blocks[0];
        let pb = &parsed.blocks[0];
        assert!((ob.freestream.mach - pb.freestream.mach).abs() < 1e-10);
        assert!((ob.freestream.reynolds - pb.freestream.reynolds).abs() < 1.0);
        for i in 0..ob.n_points() {
            assert!((ob.q[i].rho - pb.q[i].rho).abs() < 1e-10);
            assert!((ob.q[i].rho_u - pb.q[i].rho_u).abs() < 1e-10);
            assert!((ob.q[i].e - pb.q[i].e).abs() < 1e-10);
        }
    }

    #[test]
    fn test_read_solution_error_on_empty() {
        assert!(read_plot3d_solution("").is_err());
    }

    #[test]
    fn test_read_solution_error_on_truncated() {
        assert!(read_plot3d_solution("1").is_err());
    }

    #[test]
    fn test_write_solution_contains_mach() {
        let sol = make_single_block_solution();
        let s = write_plot3d_solution(&sol);
        assert!(s.contains("5.000000000000000e-1") || s.contains("5.0000000000000"));
    }

    #[test]
    fn test_multi_block_grid_roundtrip() {
        let blk1 = Plot3dBlock {
            ni: 2,
            nj: 1,
            nk: 1,
            x: vec![0.0, 1.0],
            y: vec![0.0, 0.0],
            z: vec![0.0, 0.0],
        };
        let blk2 = Plot3dBlock {
            ni: 3,
            nj: 1,
            nk: 1,
            x: vec![1.0, 2.0, 3.0],
            y: vec![0.0, 0.0, 0.0],
            z: vec![0.0, 0.0, 0.0],
        };
        let grid = Plot3dGrid {
            blocks: vec![blk1, blk2],
        };
        let s = write_plot3d_grid(&grid);
        let parsed = read_plot3d_grid(&s).unwrap();
        assert_eq!(parsed.n_blocks(), 2);
        assert_eq!(parsed.blocks[1].ni, 3);
    }

    #[test]
    fn test_plot3d_default_grid() {
        let g = Plot3dGrid::default();
        assert_eq!(g.n_blocks(), 0);
    }

    #[test]
    fn test_plot3d_block_clone() {
        let blk = Plot3dBlock {
            ni: 1,
            nj: 1,
            nk: 1,
            x: vec![0.0],
            y: vec![0.0],
            z: vec![0.0],
        };
        let blk2 = blk.clone();
        assert_eq!(blk2.ni, 1);
    }

    #[test]
    fn test_freestream_clone() {
        let fs = Plot3dFreestream {
            mach: 1.0,
            alpha: 2.0,
            reynolds: 3.0,
            time: 4.0,
        };
        let fs2 = fs.clone();
        assert!((fs2.mach - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_qpoint_clone() {
        let pt = Plot3dQPoint {
            rho: 1.0,
            rho_u: 2.0,
            rho_v: 3.0,
            rho_w: 4.0,
            e: 5.0,
        };
        let pt2 = pt.clone();
        assert!((pt2.e - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_grid_debug_format() {
        let g = make_single_block_grid();
        let s = format!("{g:?}");
        assert!(s.contains("Plot3dGrid"));
    }

    #[test]
    fn test_solution_debug_format() {
        let sol = make_single_block_solution();
        let s = format!("{sol:?}");
        assert!(s.contains("Plot3dSolution"));
    }

    #[test]
    fn test_write_grid_has_coordinates() {
        let g = make_single_block_grid();
        let s = write_plot3d_grid(&g);
        // The x=1.0 value must be present as a floating-point token
        assert!(s.contains("1.000000000000000e0") || s.contains("1.0000000000000"));
    }

    #[test]
    fn test_solution_q_length() {
        let sol = make_single_block_solution();
        assert_eq!(sol.blocks[0].q.len(), 4);
    }

    #[test]
    fn test_block_dimensions_3d() {
        let blk = Plot3dBlock {
            ni: 3,
            nj: 4,
            nk: 5,
            x: vec![0.0; 60],
            y: vec![0.0; 60],
            z: vec![0.0; 60],
        };
        let g = Plot3dGrid { blocks: vec![blk] };
        assert_eq!(g.block_dimensions(0), Some((3, 4, 5)));
        assert_eq!(g.blocks[0].n_points(), 60);
    }

    #[test]
    fn test_write_read_3d_block() {
        let n = 8; // 2×2×2
        let blk = Plot3dBlock {
            ni: 2,
            nj: 2,
            nk: 2,
            x: (0..n).map(|i| i as f64 * 0.1).collect(),
            y: (0..n).map(|i| i as f64 * 0.2).collect(),
            z: (0..n).map(|i| i as f64 * 0.3).collect(),
        };
        let grid = Plot3dGrid { blocks: vec![blk] };
        let s = write_plot3d_grid(&grid);
        let parsed = read_plot3d_grid(&s).unwrap();
        for i in 0..n {
            assert!((grid.blocks[0].z[i] - parsed.blocks[0].z[i]).abs() < 1e-10);
        }
    }
}
