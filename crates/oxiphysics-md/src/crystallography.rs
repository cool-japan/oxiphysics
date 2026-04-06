// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crystallographic structure analysis.
//!
//! Provides tools for:
//! - Unit cell geometry and reciprocal lattice
//! - Miller index d-spacing and diffraction angles
//! - Structure factors and powder diffraction patterns
//! - Crystal system detection
//! - Radial distribution function, coordination number, bond angles
//! - FCC and BCC lattice generation

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// UnitCell
// ---------------------------------------------------------------------------

/// Crystallographic unit cell defined by lattice parameters.
///
/// Angles `alpha`, `beta`, `gamma` are in **radians**.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct UnitCell {
    /// Lattice parameter a (Å).
    pub a: f64,
    /// Lattice parameter b (Å).
    pub b: f64,
    /// Lattice parameter c (Å).
    pub c: f64,
    /// Angle α between b and c (radians).
    pub alpha: f64,
    /// Angle β between a and c (radians).
    pub beta: f64,
    /// Angle γ between a and b (radians).
    pub gamma: f64,
}

impl UnitCell {
    /// Compute the unit cell volume.
    ///
    /// Uses the general triclinic formula:
    /// V = abc√(1 - cos²α - cos²β - cos²γ + 2cosα cosβ cosγ)
    #[allow(dead_code)]
    pub fn volume(&self) -> f64 {
        let ca = self.alpha.cos();
        let cb = self.beta.cos();
        let cc = self.gamma.cos();
        let arg = 1.0 - ca * ca - cb * cb - cc * cc + 2.0 * ca * cb * cc;
        self.a * self.b * self.c * arg.max(0.0).sqrt()
    }

    /// Compute the reciprocal lattice vectors as a 3×3 matrix.
    ///
    /// Returns `[[a*_x, a*_y, a*_z\], [b*_x, ...], [c*_x, ...]]`.
    /// For a cubic cell this simplifies to 2π/a · I.
    #[allow(dead_code)]
    pub fn reciprocal_vectors(&self) -> [[f64; 3]; 3] {
        // Build direct lattice vectors in Cartesian coordinates
        let a_vec = [self.a, 0.0, 0.0];
        let b_vec = [self.b * self.gamma.cos(), self.b * self.gamma.sin(), 0.0];
        let cx = self.c * self.beta.cos();
        let cy = self.c * (self.alpha.cos() - self.beta.cos() * self.gamma.cos())
            / self.gamma.sin().max(1e-30);
        let cz = (self.c * self.c - cx * cx - cy * cy).max(0.0).sqrt();
        let c_vec = [cx, cy, cz];

        let v = self.volume();
        if v < 1e-30 {
            return [[0.0; 3]; 3];
        }

        // a* = 2π (b × c) / V, etc.
        let cross = |u: [f64; 3], w: [f64; 3]| -> [f64; 3] {
            [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ]
        };
        let scale = |v: [f64; 3], s: f64| -> [f64; 3] { [v[0] * s, v[1] * s, v[2] * s] };

        let a_star = scale(cross(b_vec, c_vec), 2.0 * PI / v);
        let b_star = scale(cross(c_vec, a_vec), 2.0 * PI / v);
        let c_star = scale(cross(a_vec, b_vec), 2.0 * PI / v);

        [a_star, b_star, c_star]
    }
}

// ---------------------------------------------------------------------------
// MillerIndex
// ---------------------------------------------------------------------------

/// Miller index notation (h, k, l) for crystallographic planes.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MillerIndex {
    /// h index.
    pub h: i32,
    /// k index.
    pub k: i32,
    /// l index.
    pub l: i32,
}

impl MillerIndex {
    /// Create a new Miller index.
    #[allow(dead_code)]
    pub fn new(h: i32, k: i32, l: i32) -> Self {
        Self { h, k, l }
    }

    /// Compute the d-spacing for a cubic unit cell.
    ///
    /// `d = a / √(h²+k²+l²)` for cubic; general triclinic uses the
    /// full metric tensor calculation.
    #[allow(dead_code)]
    pub fn d_spacing(&self, cell: &UnitCell) -> f64 {
        // General triclinic d-spacing via reciprocal metric tensor
        // d* = |h·a* + k·b* + l·c*|
        let recip = cell.reciprocal_vectors();
        let hf = self.h as f64;
        let kf = self.k as f64;
        let lf = self.l as f64;
        let qx = hf * recip[0][0] + kf * recip[1][0] + lf * recip[2][0];
        let qy = hf * recip[0][1] + kf * recip[1][1] + lf * recip[2][1];
        let qz = hf * recip[0][2] + kf * recip[1][2] + lf * recip[2][2];
        let q_mag = (qx * qx + qy * qy + qz * qz).sqrt();
        if q_mag < 1e-30 {
            f64::INFINITY
        } else {
            2.0 * PI / q_mag
        }
    }

    /// Compute the angle (in radians) between two sets of planes.
    #[allow(dead_code)]
    pub fn angle_between(a: &MillerIndex, b: &MillerIndex, cell: &UnitCell) -> f64 {
        let recip = cell.reciprocal_vectors();
        let q_vec = |m: &MillerIndex| -> [f64; 3] {
            let hf = m.h as f64;
            let kf = m.k as f64;
            let lf = m.l as f64;
            [
                hf * recip[0][0] + kf * recip[1][0] + lf * recip[2][0],
                hf * recip[0][1] + kf * recip[1][1] + lf * recip[2][1],
                hf * recip[0][2] + kf * recip[1][2] + lf * recip[2][2],
            ]
        };
        let qa = q_vec(a);
        let qb = q_vec(b);
        let dot = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2];
        let ma = (qa[0] * qa[0] + qa[1] * qa[1] + qa[2] * qa[2]).sqrt();
        let mb = (qb[0] * qb[0] + qb[1] * qb[1] + qb[2] * qb[2]).sqrt();
        if ma < 1e-30 || mb < 1e-30 {
            return 0.0;
        }
        (dot / (ma * mb)).clamp(-1.0, 1.0).acos()
    }
}

// ---------------------------------------------------------------------------
// BraggPeak
// ---------------------------------------------------------------------------

/// A single Bragg diffraction peak.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BraggPeak {
    /// Momentum transfer |q| = 2π/d (Å⁻¹).
    pub q_value: f64,
    /// Diffracted intensity (structure factor squared).
    pub intensity: f64,
    /// Miller indices of this reflection.
    pub miller: MillerIndex,
}

// ---------------------------------------------------------------------------
// Bragg condition
// ---------------------------------------------------------------------------

/// Compute Bragg scattering angles for a given d-spacing and wavelength.
///
/// Solves `2·d·sin(θ) = n·λ` for n = 1..=5.
/// Returns θ values (in radians) for which a solution exists (sin θ ≤ 1).
#[allow(dead_code)]
pub fn bragg_condition(d_spacing: f64, wavelength: f64) -> Vec<f64> {
    (1..=5)
        .filter_map(|n| {
            let sin_theta = (n as f64 * wavelength) / (2.0 * d_spacing);
            if sin_theta <= 1.0 {
                Some(sin_theta.asin())
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Structure factor
// ---------------------------------------------------------------------------

/// Compute the structure factor F(hkl) = Σ_j exp(2πi (h·x_j + k·y_j + l·z_j)).
///
/// Atoms are given as fractional coordinates `[x, y, z]`.
/// Returns `(Re F, Im F)`.
#[allow(dead_code)]
pub fn structure_factor(atoms: &[[f64; 3]], miller: &MillerIndex, _cell: &UnitCell) -> (f64, f64) {
    let hf = miller.h as f64;
    let kf = miller.k as f64;
    let lf = miller.l as f64;
    let mut re = 0.0f64;
    let mut im = 0.0f64;
    for atom in atoms {
        let phase = 2.0 * PI * (hf * atom[0] + kf * atom[1] + lf * atom[2]);
        re += phase.cos();
        im += phase.sin();
    }
    (re, im)
}

// ---------------------------------------------------------------------------
// Powder diffraction
// ---------------------------------------------------------------------------

/// Generate a simulated powder diffraction pattern.
///
/// Enumerates all (h,k,l) with |h|,|k|,|l| ≤ `hkl_max`,
/// computes d-spacing and Bragg condition, then returns a list of
/// [`BraggPeak`]s with intensity = |F(hkl)|².
#[allow(dead_code)]
pub fn powder_diffraction(
    atoms: &[[f64; 3]],
    cell: &UnitCell,
    hkl_max: i32,
    wavelength: f64,
) -> Vec<BraggPeak> {
    let mut peaks = Vec::new();
    for h in -hkl_max..=hkl_max {
        for k in -hkl_max..=hkl_max {
            for l in -hkl_max..=hkl_max {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let miller = MillerIndex::new(h, k, l);
                let d = miller.d_spacing(cell);
                if !d.is_finite() || d <= wavelength / 2.0 {
                    continue;
                }
                let thetas = bragg_condition(d, wavelength);
                if thetas.is_empty() {
                    continue;
                }
                let (re, im) = structure_factor(atoms, &miller, cell);
                let intensity = re * re + im * im;
                if intensity < 1e-6 {
                    continue; // systematic absence
                }
                let q_val = 2.0 * PI / d;
                peaks.push(BraggPeak {
                    q_value: q_val,
                    intensity,
                    miller,
                });
            }
        }
    }
    peaks.sort_by(|a, b| {
        a.q_value
            .partial_cmp(&b.q_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    peaks
}

// ---------------------------------------------------------------------------
// Crystal system
// ---------------------------------------------------------------------------

/// Crystal system classification.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CrystalSystem {
    /// a=b=c, α=β=γ=90°.
    Cubic,
    /// a=b≠c, α=β=γ=90°.
    Tetragonal,
    /// a≠b≠c, α=β=γ=90°.
    Orthorhombic,
    /// a=b≠c, α=β=90°, γ=120°.
    Hexagonal,
    /// a=b=c, α=β=γ≠90°.
    Trigonal,
    /// a≠b≠c, α=γ=90°, β≠90°.
    Monoclinic,
    /// a≠b≠c, α≠β≠γ≠90°.
    Triclinic,
}

/// Crystal symmetry information.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrystalSymmetry {
    /// The crystal system.
    pub crystal_system: CrystalSystem,
    /// International space group number (1-230).
    pub space_group: u32,
}

/// Detect the crystal system from lattice parameters.
///
/// Uses tolerance `tol` for comparing lengths and angles (in radians).
#[allow(dead_code)]
pub fn detect_crystal_system(cell: &UnitCell, tol: f64) -> CrystalSystem {
    let deg90 = PI / 2.0;
    let deg120 = 2.0 * PI / 3.0;

    let eq = |a: f64, b: f64| (a - b).abs() < tol;

    let a_eq_b = eq(cell.a, cell.b);
    let b_eq_c = eq(cell.b, cell.c);
    let alpha_90 = eq(cell.alpha, deg90);
    let beta_90 = eq(cell.beta, deg90);
    let gamma_90 = eq(cell.gamma, deg90);
    let gamma_120 = eq(cell.gamma, deg120);

    if a_eq_b && b_eq_c && alpha_90 && beta_90 && gamma_90 {
        CrystalSystem::Cubic
    } else if a_eq_b && !b_eq_c && alpha_90 && beta_90 && gamma_90 {
        CrystalSystem::Tetragonal
    } else if a_eq_b && !b_eq_c && alpha_90 && beta_90 && gamma_120 {
        CrystalSystem::Hexagonal
    } else if alpha_90 && beta_90 && gamma_90 {
        CrystalSystem::Orthorhombic
    } else if alpha_90 && gamma_90 && !beta_90 {
        CrystalSystem::Monoclinic
    } else if a_eq_b && b_eq_c && eq(cell.alpha, cell.beta) && eq(cell.beta, cell.gamma) {
        CrystalSystem::Trigonal
    } else {
        CrystalSystem::Triclinic
    }
}

// ---------------------------------------------------------------------------
// Radial distribution function
// ---------------------------------------------------------------------------

/// Compute the radial distribution function g(r) from atomic positions.
///
/// Uses a simple histogram with `n_bins` bins from 0 to `box_size/2`.
/// Returns `(r, g(r))` pairs, normalized so that g→1 at large r.
#[allow(dead_code)]
pub fn radial_distribution_function(
    atoms: &[[f64; 3]],
    box_size: f64,
    n_bins: usize,
) -> Vec<(f64, f64)> {
    let n = atoms.len();
    if n < 2 || n_bins == 0 {
        return vec![(0.0, 0.0); n_bins];
    }
    let r_max = box_size / 2.0;
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];

    for i in 0..n {
        for j in (i + 1)..n {
            let mut dx = atoms[i][0] - atoms[j][0];
            let mut dy = atoms[i][1] - atoms[j][1];
            let mut dz = atoms[i][2] - atoms[j][2];
            // Minimum image convention
            dx -= box_size * (dx / box_size).round();
            dy -= box_size * (dy / box_size).round();
            dz -= box_size * (dz / box_size).round();
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let bin = (r / dr) as usize;
                let bin = bin.min(n_bins - 1);
                hist[bin] += 2; // count both i-j and j-i
            }
        }
    }

    let rho = n as f64 / (box_size * box_size * box_size);
    (0..n_bins)
        .map(|b| {
            let r = (b as f64 + 0.5) * dr;
            let shell_vol = 4.0 * PI * r * r * dr;
            let ideal = rho * shell_vol * n as f64;
            let gr = if ideal > 1e-30 {
                hist[b] as f64 / ideal
            } else {
                0.0
            };
            (r, gr)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Coordination number
// ---------------------------------------------------------------------------

/// Count the number of neighbors within `cutoff` distance for each atom.
///
/// Uses periodic boundary conditions with box of size `box_size` (cubic).
#[allow(dead_code)]
pub fn coordination_number(atoms: &[[f64; 3]], cutoff: f64) -> Vec<usize> {
    let n = atoms.len();
    let mut coord = vec![0usize; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = atoms[i][0] - atoms[j][0];
            let dy = atoms[i][1] - atoms[j][1];
            let dz = atoms[i][2] - atoms[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < cutoff * cutoff {
                coord[i] += 1;
            }
        }
    }
    coord
}

// ---------------------------------------------------------------------------
// Bond angle distribution
// ---------------------------------------------------------------------------

/// Compute the bond angle distribution for all triplets i-j-k within cutoff.
///
/// Returns `(angle_deg, count)` histogram with `n_bins` bins from 0° to 180°.
#[allow(dead_code)]
pub fn bond_angle_distribution(atoms: &[[f64; 3]], cutoff: f64, n_bins: usize) -> Vec<(f64, f64)> {
    let n = atoms.len();
    let mut hist = vec![0u64; n_bins];
    let dtheta = PI / n_bins as f64;

    for j in 0..n {
        // Find neighbors of j
        let mut neighbors = Vec::new();
        for i in 0..n {
            if i == j {
                continue;
            }
            let dx = atoms[i][0] - atoms[j][0];
            let dy = atoms[i][1] - atoms[j][1];
            let dz = atoms[i][2] - atoms[j][2];
            if (dx * dx + dy * dy + dz * dz).sqrt() < cutoff {
                neighbors.push(i);
            }
        }
        // All pairs of neighbors
        for m in 0..neighbors.len() {
            for p in (m + 1)..neighbors.len() {
                let i = neighbors[m];
                let k = neighbors[p];
                let ux = atoms[i][0] - atoms[j][0];
                let uy = atoms[i][1] - atoms[j][1];
                let uz = atoms[i][2] - atoms[j][2];
                let vx = atoms[k][0] - atoms[j][0];
                let vy = atoms[k][1] - atoms[j][1];
                let vz = atoms[k][2] - atoms[j][2];
                let mu = (ux * ux + uy * uy + uz * uz).sqrt();
                let mv = (vx * vx + vy * vy + vz * vz).sqrt();
                if mu < 1e-30 || mv < 1e-30 {
                    continue;
                }
                let cos_t = ((ux * vx + uy * vy + uz * vz) / (mu * mv)).clamp(-1.0, 1.0);
                let theta = cos_t.acos();
                let bin = (theta / dtheta) as usize;
                let bin = bin.min(n_bins - 1);
                hist[bin] += 1;
            }
        }
    }

    let total: u64 = hist.iter().sum();
    (0..n_bins)
        .map(|b| {
            let theta_deg = (b as f64 + 0.5) * dtheta * 180.0 / PI;
            let frac = if total > 0 {
                hist[b] as f64 / total as f64
            } else {
                0.0
            };
            (theta_deg, frac)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// FCC and BCC lattice generators
// ---------------------------------------------------------------------------

/// Generate an FCC lattice with n×n×n conventional unit cells.
///
/// Each conventional FCC cell has 4 atoms at:
/// `(0,0,0)`, `(a/2,a/2,0)`, `(a/2,0,a/2)`, `(0,a/2,a/2)`.
///
/// Returns atom positions in Cartesian coordinates (Å).
#[allow(dead_code)]
pub fn fcc_lattice(a: f64, n: usize) -> Vec<[f64; 3]> {
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut atoms = Vec::new();
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                for b in &basis {
                    atoms.push([
                        (ix as f64 + b[0]) * a,
                        (iy as f64 + b[1]) * a,
                        (iz as f64 + b[2]) * a,
                    ]);
                }
            }
        }
    }
    atoms
}

/// Generate a BCC lattice with n×n×n conventional unit cells.
///
/// Each conventional BCC cell has 2 atoms at:
/// `(0,0,0)` and `(a/2,a/2,a/2)`.
///
/// Returns atom positions in Cartesian coordinates (Å).
#[allow(dead_code)]
pub fn bcc_lattice(a: f64, n: usize) -> Vec<[f64; 3]> {
    let basis = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]];
    let mut atoms = Vec::new();
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                for b in &basis {
                    atoms.push([
                        (ix as f64 + b[0]) * a,
                        (iy as f64 + b[1]) * a,
                        (iz as f64 + b[2]) * a,
                    ]);
                }
            }
        }
    }
    atoms
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic_cell(a: f64) -> UnitCell {
        UnitCell {
            a,
            b: a,
            c: a,
            alpha: PI / 2.0,
            beta: PI / 2.0,
            gamma: PI / 2.0,
        }
    }

    #[test]
    fn test_cubic_volume() {
        let cell = cubic_cell(4.0);
        assert!((cell.volume() - 64.0).abs() < 1e-10);
    }

    #[test]
    fn test_orthorhombic_volume() {
        let cell = UnitCell {
            a: 2.0,
            b: 3.0,
            c: 4.0,
            alpha: PI / 2.0,
            beta: PI / 2.0,
            gamma: PI / 2.0,
        };
        assert!((cell.volume() - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_reciprocal_vectors_cubic() {
        let a = 3.0;
        let cell = cubic_cell(a);
        let recip = cell.reciprocal_vectors();
        // a* = [2π/a, 0, 0]
        assert!((recip[0][0] - 2.0 * PI / a).abs() < 1e-10);
        assert!(recip[0][1].abs() < 1e-10);
        assert!(recip[0][2].abs() < 1e-10);
    }

    #[test]
    fn test_miller_index_new() {
        let m = MillerIndex::new(1, 2, 3);
        assert_eq!(m.h, 1);
        assert_eq!(m.k, 2);
        assert_eq!(m.l, 3);
    }

    #[test]
    fn test_d_spacing_100_cubic() {
        let a = 3.0;
        let cell = cubic_cell(a);
        let m = MillerIndex::new(1, 0, 0);
        let d = m.d_spacing(&cell);
        assert!(
            (d - a).abs() < 1e-8,
            "d(100) should equal a={}, got {}",
            a,
            d
        );
    }

    #[test]
    fn test_d_spacing_110_cubic() {
        let a = 4.0;
        let cell = cubic_cell(a);
        let m = MillerIndex::new(1, 1, 0);
        let d = m.d_spacing(&cell);
        let expected = a / 2.0_f64.sqrt();
        assert!((d - expected).abs() < 1e-8);
    }

    #[test]
    fn test_d_spacing_111_cubic() {
        let a = 3.0;
        let cell = cubic_cell(a);
        let m = MillerIndex::new(1, 1, 1);
        let d = m.d_spacing(&cell);
        let expected = a / 3.0_f64.sqrt();
        assert!((d - expected).abs() < 1e-7);
    }

    #[test]
    fn test_d_spacing_decreases_higher_hkl() {
        let cell = cubic_cell(4.0);
        let d1 = MillerIndex::new(1, 0, 0).d_spacing(&cell);
        let d2 = MillerIndex::new(2, 0, 0).d_spacing(&cell);
        let d3 = MillerIndex::new(3, 0, 0).d_spacing(&cell);
        assert!(d1 > d2 && d2 > d3);
    }

    #[test]
    fn test_angle_between_perpendicular_planes() {
        let cell = cubic_cell(4.0);
        let h = MillerIndex::new(1, 0, 0);
        let k = MillerIndex::new(0, 1, 0);
        let angle = MillerIndex::angle_between(&h, &k, &cell);
        assert!((angle - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_between_same_plane_is_zero() {
        let cell = cubic_cell(4.0);
        let h = MillerIndex::new(1, 1, 0);
        let angle = MillerIndex::angle_between(&h, &h, &cell);
        assert!(angle.abs() < 1e-10);
    }

    #[test]
    fn test_bragg_condition_n1() {
        let d = 2.0;
        let lam = 1.54;
        let thetas = bragg_condition(d, lam);
        assert!(!thetas.is_empty());
        let theta = thetas[0];
        let lhs = 2.0 * d * theta.sin();
        assert!((lhs - lam).abs() < 1e-10);
    }

    #[test]
    fn test_bragg_condition_returns_up_to_5_orders() {
        let d = 10.0;
        let lam = 1.0;
        let thetas = bragg_condition(d, lam);
        assert_eq!(thetas.len(), 5);
    }

    #[test]
    fn test_bragg_condition_no_solution_short_d() {
        // d < lambda/2 -> no solution
        let thetas = bragg_condition(0.3, 1.0);
        assert!(thetas.is_empty());
    }

    #[test]
    fn test_structure_factor_fcc_000() {
        // FCC fractional positions: (0,0,0),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)
        let atoms = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let cell = cubic_cell(4.0);
        // (h+k), (h+l), (k+l) all even -> F = 4
        let miller = MillerIndex::new(2, 0, 0);
        let (re, im) = structure_factor(&atoms, &miller, &cell);
        assert!(
            (re - 4.0).abs() < 1e-10,
            "F(200) for FCC should be 4, got {}",
            re
        );
        assert!(im.abs() < 1e-10);
    }

    #[test]
    fn test_structure_factor_fcc_systematic_absence() {
        let atoms = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let cell = cubic_cell(4.0);
        // (100) is a systematic absence for FCC (mixed h+k,h+l,k+l)
        let miller = MillerIndex::new(1, 0, 0);
        let (re, im) = structure_factor(&atoms, &miller, &cell);
        let intensity = re * re + im * im;
        assert!(
            intensity < 1e-10,
            "F(100) for FCC should be ~0, got {}",
            intensity
        );
    }

    #[test]
    fn test_fcc_atoms_per_unit_cell() {
        let atoms = fcc_lattice(4.0, 1);
        assert_eq!(atoms.len(), 4, "FCC has 4 atoms per unit cell");
    }

    #[test]
    fn test_fcc_lattice_n2() {
        let atoms = fcc_lattice(4.0, 2);
        assert_eq!(atoms.len(), 4 * 8);
    }

    #[test]
    fn test_bcc_atoms_per_unit_cell() {
        let atoms = bcc_lattice(3.0, 1);
        assert_eq!(atoms.len(), 2, "BCC has 2 atoms per unit cell");
    }

    #[test]
    fn test_bcc_lattice_n3() {
        let atoms = bcc_lattice(3.0, 3);
        assert_eq!(atoms.len(), 2 * 27);
    }

    #[test]
    fn test_detect_cubic() {
        let cell = cubic_cell(4.0);
        let sys = detect_crystal_system(&cell, 1e-6);
        assert_eq!(sys, CrystalSystem::Cubic);
    }

    #[test]
    fn test_detect_tetragonal() {
        let cell = UnitCell {
            a: 3.0,
            b: 3.0,
            c: 5.0,
            alpha: PI / 2.0,
            beta: PI / 2.0,
            gamma: PI / 2.0,
        };
        let sys = detect_crystal_system(&cell, 1e-6);
        assert_eq!(sys, CrystalSystem::Tetragonal);
    }

    #[test]
    fn test_detect_orthorhombic() {
        let cell = UnitCell {
            a: 2.0,
            b: 3.0,
            c: 4.0,
            alpha: PI / 2.0,
            beta: PI / 2.0,
            gamma: PI / 2.0,
        };
        let sys = detect_crystal_system(&cell, 1e-6);
        assert_eq!(sys, CrystalSystem::Orthorhombic);
    }

    #[test]
    fn test_detect_monoclinic() {
        let cell = UnitCell {
            a: 2.0,
            b: 3.0,
            c: 4.0,
            alpha: PI / 2.0,
            beta: 1.2,
            gamma: PI / 2.0,
        };
        let sys = detect_crystal_system(&cell, 1e-6);
        assert_eq!(sys, CrystalSystem::Monoclinic);
    }

    #[test]
    fn test_detect_hexagonal() {
        let cell = UnitCell {
            a: 3.0,
            b: 3.0,
            c: 5.0,
            alpha: PI / 2.0,
            beta: PI / 2.0,
            gamma: 2.0 * PI / 3.0,
        };
        let sys = detect_crystal_system(&cell, 1e-6);
        assert_eq!(sys, CrystalSystem::Hexagonal);
    }

    #[test]
    fn test_rdf_fcc_first_peak() {
        // FCC with a=4.0: nearest neighbor at a/sqrt(2) = 2.83
        let atoms = fcc_lattice(4.0, 3);
        let box_size = 4.0 * 3.0;
        let rdf = radial_distribution_function(&atoms, box_size, 100);
        // Find first non-trivial peak
        let max_bin = rdf
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.1.partial_cmp(&y.1.1).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);
        let r_peak = rdf[max_bin].0;
        let expected = 4.0 / 2.0_f64.sqrt();
        assert!(
            (r_peak - expected).abs() < 0.5,
            "RDF first peak at {:.3}, expected {:.3}",
            r_peak,
            expected
        );
    }

    #[test]
    fn test_rdf_length() {
        let atoms = fcc_lattice(4.0, 2);
        let rdf = radial_distribution_function(&atoms, 8.0, 50);
        assert_eq!(rdf.len(), 50);
    }

    #[test]
    fn test_coordination_number_fcc() {
        // FCC with a=4.0: nearest neighbor at a/sqrt(2) = 2.83
        // Use n=4 (4^3=64 unit cells, 256 atoms) so interior atoms dominate
        let a = 4.0_f64;
        let atoms = fcc_lattice(a, 4);
        let cutoff = a / 2.0_f64.sqrt() + 0.1;
        let coord = coordination_number(&atoms, cutoff);
        // Count atoms with CN=12 (interior) vs fewer (surface)
        let max_cn = *coord.iter().max().unwrap_or(&0);
        assert_eq!(
            max_cn, 12,
            "FCC max coordination number should be 12, got {}",
            max_cn
        );
    }

    #[test]
    fn test_coordination_number_bcc() {
        // BCC with a=3.0: nearest neighbor at a*sqrt(3)/2 = 2.60
        let a = 3.0_f64;
        let atoms = bcc_lattice(a, 4);
        let cutoff = a * 3.0_f64.sqrt() / 2.0 + 0.1;
        let coord = coordination_number(&atoms, cutoff);
        let max_cn = *coord.iter().max().unwrap_or(&0);
        assert_eq!(
            max_cn, 8,
            "BCC max coordination number should be 8, got {}",
            max_cn
        );
    }

    #[test]
    fn test_bond_angle_distribution_length() {
        let atoms = fcc_lattice(4.0, 2);
        let bad = bond_angle_distribution(&atoms, 3.5, 36);
        assert_eq!(bad.len(), 36);
    }

    #[test]
    fn test_bond_angle_distribution_normalised() {
        let atoms = fcc_lattice(4.0, 2);
        let bad = bond_angle_distribution(&atoms, 3.5, 36);
        let total: f64 = bad.iter().map(|(_, c)| c).sum();
        if total > 0.0 {
            assert!((total - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_powder_diffraction_returns_peaks() {
        // Simple cubic with one atom at origin
        let atoms = [[0.0, 0.0, 0.0]];
        let cell = cubic_cell(4.0);
        let peaks = powder_diffraction(&atoms, &cell, 2, 1.54);
        assert!(!peaks.is_empty());
    }

    #[test]
    fn test_powder_diffraction_q_sorted() {
        let atoms = [[0.0, 0.0, 0.0]];
        let cell = cubic_cell(4.0);
        let peaks = powder_diffraction(&atoms, &cell, 2, 1.54);
        for w in peaks.windows(2) {
            assert!(w[0].q_value <= w[1].q_value);
        }
    }
}
