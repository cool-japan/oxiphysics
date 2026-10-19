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
//! `rho`(1..n)  <rho*u>(1..n)  <rho*v>(1..n)  <rho*w>(1..n)  `e`(1..n)   (per block)
//! ```
//!
//! The solution of a block is stored *variable-major*, as the standard
//! `((((Q(i,j,k,m), i=1,ni), j=1,nj), k=1,nk), m=1,5)` Fortran loop writes it:
//! all `n = ni*nj*nk` values of density first, then all values of the
//! x-momentum, the y-momentum, the z-momentum and finally the total energy.
//! Within one variable the index `i` runs fastest.  The grid file is likewise
//! stored coordinate-major (all `x`, then all `y`, then all `z`).

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

/// Number of grid points `ni * nj * nk` of a block.
///
/// # Errors
/// Returns `Error::Parse` if the product does not fit in `usize`.
fn block_point_count(ni: usize, nj: usize, nk: usize) -> Result<usize, IoError> {
    ni.checked_mul(nj)
        .and_then(|plane| plane.checked_mul(nk))
        .ok_or_else(|| {
            IoError::Parse(format!(
                "PLOT3D block dimensions {ni} x {nj} x {nk} overflow the point count"
            ))
        })
}

/// Read `nblocks` dimension triples `(ni, nj, nk)`.
///
/// `max_blocks` bounds the up-front allocation: a count that is larger than
/// the input can possibly describe fails on the first missing token instead
/// of reserving memory for it.
fn read_block_dimensions(
    iter: &mut SplitAsciiWhitespace<'_>,
    nblocks: usize,
    max_blocks: usize,
) -> Result<Vec<(usize, usize, usize)>, IoError> {
    let mut dims = Vec::with_capacity(nblocks.min(max_blocks));
    for _ in 0..nblocks {
        let ni = next_usize(iter)?;
        let nj = next_usize(iter)?;
        let nk = next_usize(iter)?;
        dims.push((ni, nj, nk));
    }
    Ok(dims)
}

/// Read `n` consecutive values (one PLOT3D array).
///
/// `max_values` bounds the up-front allocation, see [`read_block_dimensions`].
fn read_values(
    iter: &mut SplitAsciiWhitespace<'_>,
    n: usize,
    max_values: usize,
) -> Result<Vec<f64>, IoError> {
    let mut values = Vec::with_capacity(n.min(max_values));
    for _ in 0..n {
        values.push(next_f64(iter)?);
    }
    Ok(values)
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
    let dims = read_block_dimensions(&mut iter, nblocks, src.len())?;
    let mut blocks = Vec::with_capacity(dims.len());
    for (ni, nj, nk) in dims {
        let n = block_point_count(ni, nj, nk)?;
        let x = read_values(&mut iter, n, src.len())?;
        let y = read_values(&mut iter, n, src.len())?;
        let z = read_values(&mut iter, n, src.len())?;
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
    /// Solution data, one entry per grid point (`ni * nj * nk` entries, the
    /// index `i` runs fastest).  The file stores the variables of all points
    /// one after the other, see [`read_plot3d_solution`].
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
/// Each block starts with its free-stream line (`mach alpha rey time`) and
/// holds the five conservative variables variable-major: `n = ni*nj*nk`
/// densities, then `n` x-momenta, `n` y-momenta, `n` z-momenta and `n` total
/// energies (see the [module documentation](self)).  The values are returned
/// as one [`Plot3dQPoint`] per grid point.
///
/// # Errors
/// Returns `Error::Parse` if the input is malformed, truncated, or declares
/// block dimensions whose point count overflows `usize`.
pub fn read_plot3d_solution(src: &str) -> Result<Plot3dSolution, IoError> {
    let mut iter = src.split_ascii_whitespace();
    let nblocks = next_usize(&mut iter)?;
    let dims = read_block_dimensions(&mut iter, nblocks, src.len())?;
    let mut blocks = Vec::with_capacity(dims.len());
    for (ni, nj, nk) in dims {
        let mach = next_f64(&mut iter)?;
        let alpha = next_f64(&mut iter)?;
        let reynolds = next_f64(&mut iter)?;
        let time = next_f64(&mut iter)?;
        let n = block_point_count(ni, nj, nk)?;
        let rho = read_values(&mut iter, n, src.len())?;
        let rho_u = read_values(&mut iter, n, src.len())?;
        let rho_v = read_values(&mut iter, n, src.len())?;
        let rho_w = read_values(&mut iter, n, src.len())?;
        let e = read_values(&mut iter, n, src.len())?;
        let q = rho
            .into_iter()
            .zip(rho_u)
            .zip(rho_v)
            .zip(rho_w)
            .zip(e)
            .map(|((((rho, rho_u), rho_v), rho_w), e)| Plot3dQPoint {
                rho,
                rho_u,
                rho_v,
                rho_w,
                e,
            })
            .collect();
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

/// Append one value per grid point, taken from `variable`, one per line.
fn write_q_variable(out: &mut String, q: &[Plot3dQPoint], variable: fn(&Plot3dQPoint) -> f64) {
    for pt in q {
        let _ = writeln!(out, "{:.15e}", variable(pt));
    }
}

/// Write a PLOT3D solution to an ASCII string.
///
/// The five conservative variables of a block are written variable-major (all
/// densities, then all x-momenta, ...), the layout [`read_plot3d_solution`]
/// and other PLOT3D tools expect.  Every block's `q` must hold exactly
/// `ni * nj * nk` points; a block with a different number of points produces
/// a file whose arrays do not match its declared dimensions.
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
        write_q_variable(&mut out, &blk.q, |pt| pt.rho);
        write_q_variable(&mut out, &blk.q, |pt| pt.rho_u);
        write_q_variable(&mut out, &blk.q, |pt| pt.rho_v);
        write_q_variable(&mut out, &blk.q, |pt| pt.rho_w);
        write_q_variable(&mut out, &blk.q, |pt| pt.e);
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

    // ── variable-major `.q` layout (F2265) ───────────────────────────────────

    /// A standard formatted PLOT3D function file written by hand: a single
    /// 2x1x1 block with distinct values for every variable.  The Fortran loop
    /// `((((Q(i,j,k,m), i=1,id), j=1,jd), k=1,kd), m=1,5)` puts all densities
    /// first, then all x-momenta, and so on.
    const HAND_WRITTEN_Q: &str = "\
1
2 1 1
0.5 2.0 1.0e6 0.25
1.0 1.1
10.0 11.0
20.0 21.0
30.0 31.0
100.0 110.0
";

    #[test]
    fn read_solution_parses_the_standard_variable_major_layout() {
        let sol = read_plot3d_solution(HAND_WRITTEN_Q).expect("standard file must parse");
        assert_eq!(sol.n_blocks(), 1);
        let blk = &sol.blocks[0];
        assert_eq!((blk.ni, blk.nj, blk.nk), (2, 1, 1));
        assert_eq!(blk.freestream.mach, 0.5);
        assert_eq!(blk.freestream.alpha, 2.0);
        assert_eq!(blk.freestream.reynolds, 1.0e6);
        assert_eq!(blk.freestream.time, 0.25);
        assert_eq!(blk.q.len(), 2);
        assert_eq!(
            blk.q[0],
            Plot3dQPoint {
                rho: 1.0,
                rho_u: 10.0,
                rho_v: 20.0,
                rho_w: 30.0,
                e: 100.0,
            }
        );
        assert_eq!(
            blk.q[1],
            Plot3dQPoint {
                rho: 1.1,
                rho_u: 11.0,
                rho_v: 21.0,
                rho_w: 31.0,
                e: 110.0,
            }
        );
    }

    #[test]
    fn write_solution_emits_all_densities_before_any_momentum() {
        let sol = read_plot3d_solution(HAND_WRITTEN_Q).expect("fixture must parse");
        let text = write_plot3d_solution(&sol);
        let tokens: Vec<&str> = text.split_ascii_whitespace().collect();
        // nblocks, one dimension triple, four free-stream values, 5 x 2 values.
        assert_eq!(tokens.len(), 1 + 3 + 4 + 10);
        assert_eq!(&tokens[..4], ["1", "2", "1", "1"]);
        let values: Vec<f64> = tokens[8..]
            .iter()
            .map(|t| t.parse::<f64>().expect("numeric token"))
            .collect();
        assert_eq!(
            values,
            vec![1.0, 1.1, 10.0, 11.0, 20.0, 21.0, 30.0, 31.0, 100.0, 110.0]
        );
    }

    #[test]
    fn solution_roundtrip_keeps_every_variable_of_a_multi_block_file() {
        let make_block = |ni: usize, nj: usize, nk: usize, base: f64| {
            let n = ni * nj * nk;
            Plot3dSolutionBlock {
                ni,
                nj,
                nk,
                freestream: Plot3dFreestream {
                    mach: 0.8 + base,
                    alpha: 1.5,
                    reynolds: 3.0e6,
                    time: base,
                },
                q: (0..n)
                    .map(|p| {
                        let p = p as f64;
                        Plot3dQPoint {
                            rho: base + p,
                            rho_u: base + 100.0 + p,
                            rho_v: base + 200.0 + p,
                            rho_w: base + 300.0 + p,
                            e: base + 400.0 + p,
                        }
                    })
                    .collect(),
            }
        };
        let original = Plot3dSolution {
            blocks: vec![make_block(3, 2, 2, 1.0), make_block(2, 1, 4, 5000.0)],
        };
        let parsed = read_plot3d_solution(&write_plot3d_solution(&original))
            .expect("written file must parse");
        assert_eq!(parsed.n_blocks(), 2);
        for (ob, pb) in original.blocks.iter().zip(&parsed.blocks) {
            assert_eq!((pb.ni, pb.nj, pb.nk), (ob.ni, ob.nj, ob.nk));
            assert!((pb.freestream.mach - ob.freestream.mach).abs() < 1e-12);
            assert!((pb.freestream.alpha - ob.freestream.alpha).abs() < 1e-12);
            assert!((pb.freestream.reynolds - ob.freestream.reynolds).abs() < 1e-3);
            assert!((pb.freestream.time - ob.freestream.time).abs() < 1e-12);
            assert_eq!(pb.q.len(), ob.q.len());
            for (op, pp) in ob.q.iter().zip(&pb.q) {
                assert!((op.rho - pp.rho).abs() < 1e-9 * op.rho.abs().max(1.0));
                assert!((op.rho_u - pp.rho_u).abs() < 1e-9 * op.rho_u.abs().max(1.0));
                assert!((op.rho_v - pp.rho_v).abs() < 1e-9 * op.rho_v.abs().max(1.0));
                assert!((op.rho_w - pp.rho_w).abs() < 1e-9 * op.rho_w.abs().max(1.0));
                assert!((op.e - pp.e).abs() < 1e-9 * op.e.abs().max(1.0));
            }
        }
    }

    #[test]
    fn solution_missing_a_variable_array_is_an_error() {
        // The last array (energy) is absent.
        let truncated = "1\n2 1 1\n0.5 2.0 1.0e6 0.25\n1.0 1.1\n10.0 11.0\n20.0 21.0\n30.0 31.0\n";
        assert!(matches!(
            read_plot3d_solution(truncated),
            Err(IoError::Parse(_))
        ));
    }

    #[test]
    fn solution_with_overflowing_dimensions_is_a_parse_error() {
        let src = format!("1\n{} 2 2\n0.5 0.0 1.0 0.0\n", usize::MAX);
        match read_plot3d_solution(&src) {
            Err(IoError::Parse(msg)) => assert!(msg.contains("overflow"), "{msg}"),
            other => panic!("expected an overflow parse error, got {other:?}"),
        }
    }

    #[test]
    fn grid_with_overflowing_dimensions_is_a_parse_error() {
        let src = format!("1\n{} {} 3\n", usize::MAX / 2, usize::MAX / 2);
        match read_plot3d_grid(&src) {
            Err(IoError::Parse(msg)) => assert!(msg.contains("overflow"), "{msg}"),
            other => panic!("expected an overflow parse error, got {other:?}"),
        }
    }

    #[test]
    fn huge_declared_dimensions_fail_without_reserving_memory() {
        // 10^12 points do not overflow usize on 64-bit targets, but the input
        // holds a handful of tokens: the reader must report a truncated file
        // instead of trying to allocate terabytes.
        let q_src = "1\n1000000 1000000 1\n0.5 0.0 1.0 0.0\n1.0 2.0\n";
        assert!(matches!(
            read_plot3d_solution(q_src),
            Err(IoError::Parse(_))
        ));
        let g_src = "1\n1000000 1000000 1\n0.0 1.0\n";
        assert!(matches!(read_plot3d_grid(g_src), Err(IoError::Parse(_))));
        // An absurd block count is bounded in the same way.
        assert!(matches!(
            read_plot3d_solution("999999999999\n1 1 1\n"),
            Err(IoError::Parse(_))
        ));
    }
}
