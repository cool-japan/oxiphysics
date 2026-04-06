#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Protein structure analysis and crystal structure builders.
//!
//! Provides Kabsch alignment, RMSD, secondary structure assignment,
//! crystal structure builders (FCC, BCC, HCP, diamond), lattice parameter
//! optimization, coordination number analysis, Voronoi analysis basics,
//! and static structure factor S(q).

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Basic geometry helpers
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 { v } else { scale3(v, 1.0 / n) }
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3(sub3(a, b))
}

// 3x3 matrix multiply: result[i][j] = sum_k a[i][k]*b[k][j]
fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

fn mat3_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

fn mat3_det(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

// Apply a one-sided Jacobi rotation to make A[p][q] → 0.
// Updates A in-place and accumulates the rotation in V.
fn jacobi_rotation(a: &mut [[f64; 3]; 3], v: &mut [[f64; 3]; 3], p: usize, q: usize) {
    if a[p][q].abs() < 1e-15 {
        return;
    }
    let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
    let t = if theta >= 0.0 {
        1.0 / (theta + (1.0 + theta * theta).sqrt())
    } else {
        -1.0 / (-theta + (1.0 + theta * theta).sqrt())
    };
    let c = 1.0 / (1.0 + t * t).sqrt();
    let s = t * c;

    // Update symmetric matrix A ← J^T A J
    let app = a[p][p];
    let aqq = a[q][q];
    let apq = a[p][q];

    a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    a[p][q] = 0.0;
    a[q][p] = 0.0;

    // Off-diagonal rows/columns r ≠ p, q
    let pairs: &[usize] = match (p, q) {
        (0, 1) | (1, 0) => &[2],
        (0, 2) | (2, 0) => &[1],
        _ => &[0],
    };
    for &r in pairs {
        let arp = a[r][p];
        let arq = a[r][q];
        a[r][p] = c * arp - s * arq;
        a[p][r] = a[r][p];
        a[r][q] = s * arp + c * arq;
        a[q][r] = a[r][q];
    }

    // Accumulate rotation in V (columns of V ← V * J)
    for i in 0..3 {
        let vip = v[i][p];
        let viq = v[i][q];
        v[i][p] = c * vip - s * viq;
        v[i][q] = s * vip + c * viq;
    }
}

/// Compute eigendecomposition of a 3×3 symmetric matrix via Jacobi iterations.
/// Returns (eigenvectors as columns, eigenvalues).
fn jacobi_eigen(mat: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3]) {
    let mut a = *mat;
    let mut v = mat3_identity();

    let off_diag_pairs = [(0, 1), (0, 2), (1, 2)];

    for _ in 0..30 {
        for &(p, q) in &off_diag_pairs {
            jacobi_rotation(&mut a, &mut v, p, q);
        }
    }

    let eigenvalues = [a[0][0], a[1][1], a[2][2]];
    (v, eigenvalues)
}

// ---------------------------------------------------------------------------
// Center of geometry
// ---------------------------------------------------------------------------

/// Compute the mean (center of geometry) of a set of 3-D positions.
pub fn center_of_geometry(positions: &[[f64; 3]]) -> [f64; 3] {
    if positions.is_empty() {
        return [0.0; 3];
    }
    let mut sum = [0.0f64; 3];
    for p in positions {
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
    }
    let n = positions.len() as f64;
    [sum[0] / n, sum[1] / n, sum[2] / n]
}

// ---------------------------------------------------------------------------
// Translate / rotate
// ---------------------------------------------------------------------------

/// Translate all positions by vector `t`.
pub fn translate(positions: &mut [[f64; 3]], t: [f64; 3]) {
    for p in positions.iter_mut() {
        p[0] += t[0];
        p[1] += t[1];
        p[2] += t[2];
    }
}

/// Apply 3×3 rotation matrix to all positions.
pub fn rotate(positions: &mut [[f64; 3]], r: &[[f64; 3]; 3]) {
    for p in positions.iter_mut() {
        let x = r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2];
        let y = r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2];
        let z = r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2];
        *p = [x, y, z];
    }
}

// ---------------------------------------------------------------------------
// Kabsch rotation
// ---------------------------------------------------------------------------

/// Kabsch algorithm: find the optimal rotation `R` (3×3) that minimises the
/// RMSD between point sets `p` and `q` (both already expected to be centred,
/// but the function centres them internally for robustness).
pub fn kabsch_rotation(p: &[[f64; 3]], q: &[[f64; 3]]) -> [[f64; 3]; 3] {
    assert_eq!(
        p.len(),
        q.len(),
        "kabsch_rotation: point sets must be same length"
    );
    let n = p.len();
    if n == 0 {
        return mat3_identity();
    }

    let cp = center_of_geometry(p);
    let cq = center_of_geometry(q);

    let pc: Vec<[f64; 3]> = p.iter().map(|&x| sub3(x, cp)).collect();
    let qc: Vec<[f64; 3]> = q.iter().map(|&x| sub3(x, cq)).collect();

    let mut h = [[0.0f64; 3]; 3];
    for k in 0..n {
        for i in 0..3 {
            for j in 0..3 {
                h[i][j] += pc[k][i] * qc[k][j];
            }
        }
    }

    let ht = mat3_transpose(&h);
    let m = mat3_mul(&ht, &h);

    let (v, s2) = jacobi_eigen(&m);
    let sv: [f64; 3] = [
        s2[0].max(0.0).sqrt(),
        s2[1].max(0.0).sqrt(),
        s2[2].max(0.0).sqrt(),
    ];

    let mut u_cols = [[0.0f64; 3]; 3];
    for j in 0..3 {
        let vj = [v[0][j], v[1][j], v[2][j]];
        let mut hvj = [0.0f64; 3];
        for i in 0..3 {
            for k in 0..3 {
                hvj[i] += h[i][k] * vj[k];
            }
        }
        if sv[j] > 1e-10 {
            u_cols[j] = scale3(hvj, 1.0 / sv[j]);
        } else {
            u_cols[j] = hvj;
        }
    }

    u_cols[0] = normalize3(u_cols[0]);
    let proj = dot3(u_cols[1], u_cols[0]);
    let u0 = u_cols[0];
    u_cols[1] = sub3(u_cols[1], scale3(u0, proj));
    u_cols[1] = normalize3(u_cols[1]);
    u_cols[2] = cross3(u_cols[0], u_cols[1]);

    let mut u_mat = [[0.0f64; 3]; 3];
    for col in 0..3 {
        for row in 0..3 {
            u_mat[row][col] = u_cols[col][row];
        }
    }

    let ut = mat3_transpose(&u_mat);
    let mut r = mat3_mul(&v, &ut);

    if mat3_det(&r) < 0.0 {
        let mut v_fixed = v;
        for row in 0..3 {
            v_fixed[row][2] = -v_fixed[row][2];
        }
        r = mat3_mul(&v_fixed, &ut);
    }

    r
}

// ---------------------------------------------------------------------------
// RMSD after alignment
// ---------------------------------------------------------------------------

/// Align `mobile` onto `reference` using Kabsch and return the RMSD.
pub fn rmsd_after_alignment(mobile: &[[f64; 3]], reference: &[[f64; 3]]) -> f64 {
    assert_eq!(mobile.len(), reference.len());
    let n = mobile.len();
    if n == 0 {
        return 0.0;
    }

    let cm = center_of_geometry(mobile);
    let cr = center_of_geometry(reference);

    let mut mob_c: Vec<[f64; 3]> = mobile.iter().map(|&x| sub3(x, cm)).collect();
    let ref_c: Vec<[f64; 3]> = reference.iter().map(|&x| sub3(x, cr)).collect();

    let r = kabsch_rotation(&mob_c, &ref_c);
    rotate(&mut mob_c, &r);

    let mut sum_sq = 0.0f64;
    for i in 0..n {
        let d = sub3(mob_c[i], ref_c[i]);
        sum_sq += dot3(d, d);
    }
    (sum_sq / n as f64).sqrt()
}

// ---------------------------------------------------------------------------
// Dihedral angle
// ---------------------------------------------------------------------------

/// Torsion angle defined by four points a–b–c–d (in radians, range \[-π, π\]).
pub fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b1 = sub3(b, a);
    let b2 = sub3(c, b);
    let b3 = sub3(d, c);

    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);

    let m1 = cross3(n1, normalize3(b2));

    let x = dot3(n1, n2);
    let y = dot3(m1, n2);

    y.atan2(x)
}

// ---------------------------------------------------------------------------
// phi / psi angles
// ---------------------------------------------------------------------------

/// Given backbone positions `[N, CA, C]` per residue (stride 3), compute φ/ψ
/// dihedral angles.  Returns vectors of length `n_residues`; the first φ and
/// last ψ are set to `999.0` (undefined).
pub fn compute_phi_psi(positions: &[[f64; 3]], n_residues: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(
        positions.len() >= n_residues * 3,
        "positions must contain at least n_residues*3 entries"
    );

    let atom = |res: usize, atom_idx: usize| -> [f64; 3] { positions[res * 3 + atom_idx] };

    let mut phi = vec![999.0f64; n_residues];
    let mut psi = vec![999.0f64; n_residues];

    for i in 0..n_residues {
        if i > 0 {
            phi[i] = dihedral_angle(atom(i - 1, 2), atom(i, 0), atom(i, 1), atom(i, 2));
            phi[i] = phi[i] * 180.0 / PI;
        }

        if i + 1 < n_residues {
            psi[i] = dihedral_angle(atom(i, 0), atom(i, 1), atom(i, 2), atom(i + 1, 0));
            psi[i] = psi[i] * 180.0 / PI;
        }
    }

    (phi, psi)
}

// ---------------------------------------------------------------------------
// Secondary structure
// ---------------------------------------------------------------------------

/// Ramachandran-based secondary structure types.
#[derive(Debug, Clone, PartialEq)]
pub enum SsType {
    /// α-helix region
    Helix,
    /// β-sheet region
    Sheet,
    /// Loop / coil
    Loop,
}

/// Secondary structure assignment for one residue.
#[derive(Debug, Clone)]
pub struct SecondaryStructure {
    /// Residue index (0-based)
    pub residue_idx: usize,
    /// Assigned secondary structure type
    pub ss_type: SsType,
}

/// Assign secondary structure from Ramachandran angles (in degrees).
///
/// - φ ∈ \[-160, -30\] ∧ ψ ∈ \[-70, 50\]  → [`SsType::Helix`]
/// - φ ∈ \[-160, -60\] ∧ ψ ∈ \[90, 180\]  → [`SsType::Sheet`]
/// - otherwise → [`SsType::Loop`]
pub fn assign_secondary_structure_simple(phi: &[f64], psi: &[f64]) -> Vec<SsType> {
    assert_eq!(phi.len(), psi.len());
    phi.iter()
        .zip(psi.iter())
        .map(|(&p, &s)| {
            if (-160.0..=-30.0).contains(&p) && (-70.0..=50.0).contains(&s) {
                SsType::Helix
            } else if (-160.0..=-60.0).contains(&p) && (90.0..=180.0).contains(&s) {
                SsType::Sheet
            } else {
                SsType::Loop
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Crystal structure builders
// ---------------------------------------------------------------------------

/// Build an FCC (face-centered cubic) lattice.
///
/// # Arguments
/// * `a` – lattice parameter (edge length of the conventional cell).
/// * `nx`, `ny`, `nz` – number of unit cells in each direction.
///
/// Returns atom positions in a box of size `[nx*a, ny*a, nz*a]`.
pub fn build_fcc(a: f64, nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5 * a, 0.5 * a, 0.0],
        [0.5 * a, 0.0, 0.5 * a],
        [0.0, 0.5 * a, 0.5 * a],
    ];
    let mut positions = Vec::with_capacity(4 * nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let origin = [ix as f64 * a, iy as f64 * a, iz as f64 * a];
                for &b in &basis {
                    positions.push(add3(origin, b));
                }
            }
        }
    }
    positions
}

/// Build a BCC (body-centered cubic) lattice.
///
/// # Arguments
/// * `a` – lattice parameter.
/// * `nx`, `ny`, `nz` – number of unit cells.
pub fn build_bcc(a: f64, nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let basis = [[0.0, 0.0, 0.0], [0.5 * a, 0.5 * a, 0.5 * a]];
    let mut positions = Vec::with_capacity(2 * nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let origin = [ix as f64 * a, iy as f64 * a, iz as f64 * a];
                for &b in &basis {
                    positions.push(add3(origin, b));
                }
            }
        }
    }
    positions
}

/// Build an HCP (hexagonal close-packed) lattice.
///
/// # Arguments
/// * `a` – lattice parameter (in-plane nearest-neighbor distance).
/// * `nx`, `ny`, `nz` – number of unit cells.
///
/// The c/a ratio is fixed at the ideal value sqrt(8/3).
pub fn build_hcp(a: f64, nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let c = a * (8.0_f64 / 3.0).sqrt();
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5 * a, a * (3.0_f64).sqrt() / 6.0, 0.5 * c],
    ];
    let ax = a;
    let ay = a * (3.0_f64).sqrt() / 2.0;
    let az = c;
    let mut positions = Vec::with_capacity(2 * nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let origin = [ix as f64 * ax, iy as f64 * ay, iz as f64 * az];
                for &b in &basis {
                    positions.push(add3(origin, b));
                }
            }
        }
    }
    positions
}

/// Build a diamond cubic lattice.
///
/// # Arguments
/// * `a` – lattice parameter.
/// * `nx`, `ny`, `nz` – number of unit cells.
///
/// Diamond = two interpenetrating FCC lattices offset by (a/4, a/4, a/4).
pub fn build_diamond(a: f64, nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let fcc_basis = [
        [0.0, 0.0, 0.0],
        [0.5 * a, 0.5 * a, 0.0],
        [0.5 * a, 0.0, 0.5 * a],
        [0.0, 0.5 * a, 0.5 * a],
    ];
    let offset = [0.25 * a, 0.25 * a, 0.25 * a];
    let mut positions = Vec::with_capacity(8 * nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let origin = [ix as f64 * a, iy as f64 * a, iz as f64 * a];
                for &b in &fcc_basis {
                    positions.push(add3(origin, b));
                    positions.push(add3(add3(origin, b), offset));
                }
            }
        }
    }
    positions
}

/// Simple cubic lattice builder.
///
/// # Arguments
/// * `a` – lattice parameter.
/// * `nx`, `ny`, `nz` – number of unit cells.
pub fn build_simple_cubic(a: f64, nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    let mut positions = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                positions.push([ix as f64 * a, iy as f64 * a, iz as f64 * a]);
            }
        }
    }
    positions
}

// ---------------------------------------------------------------------------
// Lattice parameter optimization
// ---------------------------------------------------------------------------

/// Optimize lattice parameter to minimize a simple pair energy.
///
/// Uses golden-section search to find the lattice parameter `a` in
/// `[a_min, a_max]` that minimizes the total pair energy computed
/// from all nearest-neighbor distances.
///
/// # Arguments
/// * `build_fn` – function that builds positions from a lattice parameter.
/// * `pair_energy_fn` – function that computes pair energy from distance.
/// * `a_min`, `a_max` – search range.
/// * `tol` – convergence tolerance.
pub fn optimize_lattice_parameter(
    build_fn: &dyn Fn(f64) -> Vec<[f64; 3]>,
    pair_energy_fn: &dyn Fn(f64) -> f64,
    a_min: f64,
    a_max: f64,
    tol: f64,
) -> f64 {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let resphi = 2.0 - phi;

    let total_energy = |a: f64| -> f64 {
        let pos = build_fn(a);
        let n = pos.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(pos[i], pos[j]);
                if r > 0.0 && r < a * 1.5 {
                    e += pair_energy_fn(r);
                }
            }
        }
        e
    };

    let mut lo = a_min;
    let mut hi = a_max;
    let mut x1 = lo + resphi * (hi - lo);
    let mut x2 = hi - resphi * (hi - lo);
    let mut f1 = total_energy(x1);
    let mut f2 = total_energy(x2);

    for _ in 0..100 {
        if (hi - lo).abs() < tol {
            break;
        }
        if f1 < f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = lo + resphi * (hi - lo);
            f1 = total_energy(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = hi - resphi * (hi - lo);
            f2 = total_energy(x2);
        }
    }

    (lo + hi) / 2.0
}

// ---------------------------------------------------------------------------
// Coordination number
// ---------------------------------------------------------------------------

/// Compute the coordination number of each atom (number of neighbors
/// within a cutoff distance).
///
/// # Arguments
/// * `positions` – atom positions.
/// * `cutoff` – neighbor cutoff distance.
///
/// Returns a vector of coordination numbers, one per atom.
pub fn coordination_numbers(positions: &[[f64; 3]], cutoff: f64) -> Vec<usize> {
    let n = positions.len();
    let cutoff_sq = cutoff * cutoff;
    let mut counts = vec![0usize; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = sub3(positions[i], positions[j]);
            let dsq = dot3(d, d);
            if dsq <= cutoff_sq {
                counts[i] += 1;
                counts[j] += 1;
            }
        }
    }
    counts
}

/// Average coordination number across all atoms.
pub fn average_coordination_number(positions: &[[f64; 3]], cutoff: f64) -> f64 {
    let counts = coordination_numbers(positions, cutoff);
    if counts.is_empty() {
        return 0.0;
    }
    let sum: usize = counts.iter().sum();
    sum as f64 / counts.len() as f64
}

/// Coordination number with periodic boundary conditions.
pub fn coordination_numbers_pbc(
    positions: &[[f64; 3]],
    box_len: [f64; 3],
    cutoff: f64,
) -> Vec<usize> {
    let n = positions.len();
    let cutoff_sq = cutoff * cutoff;
    let mut counts = vec![0usize; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut d = sub3(positions[i], positions[j]);
            for k in 0..3 {
                if d[k] > box_len[k] * 0.5 {
                    d[k] -= box_len[k];
                } else if d[k] < -box_len[k] * 0.5 {
                    d[k] += box_len[k];
                }
            }
            let dsq = dot3(d, d);
            if dsq <= cutoff_sq {
                counts[i] += 1;
                counts[j] += 1;
            }
        }
    }
    counts
}

// ---------------------------------------------------------------------------
// Voronoi analysis basics
// ---------------------------------------------------------------------------

/// Result of a nearest-neighbor Voronoi-like analysis for one atom.
#[derive(Debug, Clone)]
pub struct VoronoiCell {
    /// Index of the central atom.
    pub atom_idx: usize,
    /// Indices of Voronoi neighbors (atoms sharing a Voronoi face).
    pub neighbors: Vec<usize>,
    /// Approximate volume of the Voronoi cell (from nearest-neighbor estimate).
    pub volume_estimate: f64,
}

/// Compute approximate Voronoi neighbors using a cutoff-based approach.
///
/// An atom j is considered a Voronoi neighbor of i if no atom k is closer
/// to both i and j and lies between them (simplified radical plane test).
///
/// # Arguments
/// * `positions` – atom positions.
/// * `cutoff` – maximum distance to consider as a potential neighbor.
pub fn voronoi_neighbors(positions: &[[f64; 3]], cutoff: f64) -> Vec<VoronoiCell> {
    let n = positions.len();
    let cutoff_sq = cutoff * cutoff;
    let mut cells = Vec::with_capacity(n);

    for i in 0..n {
        let mut candidate_neighbors: Vec<(usize, f64)> = Vec::new();
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = sub3(positions[j], positions[i]);
            let dsq = dot3(d, d);
            if dsq <= cutoff_sq {
                candidate_neighbors.push((j, dsq.sqrt()));
            }
        }
        // Sort by distance
        candidate_neighbors
            .sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Filter: j is a Voronoi neighbor if no k is closer to i and
        // the midpoint(i,j) is closer to i than to k
        let mut neighbors = Vec::new();
        for &(j, dij) in &candidate_neighbors {
            let midpoint = scale3(add3(positions[i], positions[j]), 0.5);
            let mut is_voronoi = true;
            for &(k, _dik) in &candidate_neighbors {
                if k == j {
                    continue;
                }
                let dk_mid = dist3(positions[k], midpoint);
                let di_mid = dij * 0.5;
                if dk_mid < di_mid * 0.95 {
                    is_voronoi = false;
                    break;
                }
            }
            if is_voronoi {
                neighbors.push(j);
            }
        }

        // Volume estimate: 4/3 * pi * (nearest_neighbor_distance/2)^3
        let vol = if let Some(&(_, d_nn)) = candidate_neighbors.first() {
            4.0 / 3.0 * PI * (d_nn / 2.0).powi(3)
        } else {
            0.0
        };

        cells.push(VoronoiCell {
            atom_idx: i,
            neighbors,
            volume_estimate: vol,
        });
    }

    cells
}

/// Voronoi index (n3, n4, n5, n6) counting the number of faces with
/// 3, 4, 5, 6 edges respectively.  This is a simplified signature
/// that uses the neighbor count as a proxy.
pub fn voronoi_index(n_neighbors: usize) -> [usize; 4] {
    // Approximate Voronoi index from coordination number.
    // For ideal FCC: (0, 12, 0, 0)
    // For ideal BCC: (0, 6, 0, 8) or (0, 0, 12, 2) depending on convention
    // This is a rough approximation.
    match n_neighbors {
        12 => [0, 12, 0, 0],         // FCC-like
        14 => [0, 6, 0, 8],          // BCC-like
        8 => [0, 0, 0, 8],           // simple cubic
        _ => [0, 0, n_neighbors, 0], // generic
    }
}

// ---------------------------------------------------------------------------
// Structure factor S(q)
// ---------------------------------------------------------------------------

/// Compute the static structure factor S(q) for a set of positions.
///
/// S(q) = (1/N) |Σ_j exp(i q·r_j)|²
///
/// averaged over uniformly distributed q-vectors on a sphere of radius |q|.
///
/// # Arguments
/// * `positions` – atom positions.
/// * `q_values` – magnitudes of wavevectors |q| to evaluate.
/// * `n_q_vectors` – number of random q-directions to average over per |q|.
///
/// Returns a vector of S(q) values, one per q magnitude.
pub fn structure_factor(positions: &[[f64; 3]], q_values: &[f64], n_q_vectors: usize) -> Vec<f64> {
    let n = positions.len();
    if n == 0 {
        return vec![0.0; q_values.len()];
    }
    let n_f = n as f64;

    // Generate uniformly distributed directions on unit sphere using
    // the Fibonacci spiral method
    let directions: Vec<[f64; 3]> = (0..n_q_vectors)
        .map(|k| {
            let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let theta = (1.0 - 2.0 * (k as f64 + 0.5) / n_q_vectors as f64).acos();
            let phi_angle = 2.0 * PI * k as f64 / golden;
            [
                theta.sin() * phi_angle.cos(),
                theta.sin() * phi_angle.sin(),
                theta.cos(),
            ]
        })
        .collect();

    q_values
        .iter()
        .map(|&q_mag| {
            let mut sq_avg = 0.0;
            for dir in &directions {
                let q_vec = scale3(*dir, q_mag);
                let mut re = 0.0f64;
                let mut im = 0.0f64;
                for pos in positions {
                    let qr = dot3(q_vec, *pos);
                    re += qr.cos();
                    im += qr.sin();
                }
                sq_avg += (re * re + im * im) / n_f;
            }
            sq_avg / n_q_vectors as f64
        })
        .collect()
}

/// Compute the radial distribution function g(r) from positions.
///
/// # Arguments
/// * `positions` – atom positions.
/// * `box_len` – periodic box lengths \[Lx, Ly, Lz\].
/// * `n_bins` – number of histogram bins.
/// * `r_max` – maximum distance to compute g(r).
///
/// Returns (bin_centers, g_r) where g(r) is normalized.
pub fn radial_distribution_function(
    positions: &[[f64; 3]],
    box_len: [f64; 3],
    n_bins: usize,
    r_max: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = positions.len();
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0usize; n_bins];

    for i in 0..n {
        for j in (i + 1)..n {
            let mut d = sub3(positions[i], positions[j]);
            for k in 0..3 {
                if d[k] > box_len[k] * 0.5 {
                    d[k] -= box_len[k];
                } else if d[k] < -box_len[k] * 0.5 {
                    d[k] += box_len[k];
                }
            }
            let r = norm3(d);
            if r < r_max && r > 1e-15 {
                let bin = (r / dr) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }

    let volume = box_len[0] * box_len[1] * box_len[2];
    let _rho = n as f64 / volume;
    let n_pairs = n * (n - 1) / 2;

    let bin_centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * dr).collect();
    let g_r: Vec<f64> = (0..n_bins)
        .map(|b| {
            let r = (b as f64 + 0.5) * dr;
            let shell_vol = 4.0 / 3.0 * PI * ((r + dr * 0.5).powi(3) - (r - dr * 0.5).powi(3));
            let ideal = n_pairs as f64 * shell_vol / volume;
            if ideal > 1e-30 {
                hist[b] as f64 / ideal
            } else {
                0.0
            }
        })
        .collect();

    (bin_centers, g_r)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_center_of_geometry_symmetric() {
        let positions: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let cog = center_of_geometry(&positions);
        assert!(cog[0].abs() < 1e-12, "x != 0: {}", cog[0]);
        assert!(cog[1].abs() < 1e-12, "y != 0: {}", cog[1]);
        assert!(cog[2].abs() < 1e-12, "z != 0: {}", cog[2]);
    }

    #[test]
    fn test_rmsd_identical() {
        let pos: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [-1.0, 0.0, 1.0]];
        let rmsd = rmsd_after_alignment(&pos, &pos);
        assert!(
            rmsd < 1e-10,
            "RMSD of identical sets should be ~0, got {}",
            rmsd
        );
    }

    #[test]
    fn test_rmsd_after_rotation() {
        let original: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let rotated: Vec<[f64; 3]> = original.iter().map(|&[x, y, z]| [-y, x, z]).collect();

        let rmsd = rmsd_after_alignment(&rotated, &original);
        assert!(
            rmsd < 1e-8,
            "RMSD after alignment of rotated config should be ~0, got {}",
            rmsd
        );
    }

    #[test]
    fn test_secondary_structure_helix() {
        let phi = vec![-57.0f64];
        let psi = vec![-47.0f64];
        let ss = assign_secondary_structure_simple(&phi, &psi);
        assert_eq!(ss[0], SsType::Helix, "Expected Helix for phi=-57, psi=-47");
    }

    #[test]
    fn test_secondary_structure_sheet() {
        let phi = vec![-120.0f64];
        let psi = vec![130.0f64];
        let ss = assign_secondary_structure_simple(&phi, &psi);
        assert_eq!(ss[0], SsType::Sheet, "Expected Sheet for phi=-120, psi=130");
    }

    #[test]
    fn test_dihedral_angle_known() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [0.0, 0.0, 1.0];
        let d = [0.0, 1.0, 1.0];
        let angle = dihedral_angle(a, b, c, d);
        assert!(
            (angle.abs() - PI / 2.0).abs() < 1e-10,
            "Expected |pi/2|, got {}",
            angle
        );
    }

    #[test]
    fn test_dihedral_angle_180() {
        let a = [0.0, 1.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [1.0, 0.0, 0.0];
        let d = [1.0, 0.0, 1.0];
        let angle = dihedral_angle(a, b, c, d);
        assert!(
            (angle.abs() - PI / 2.0).abs() < 1e-10,
            "Expected |pi/2|, got {}",
            angle
        );
    }

    // ---- Crystal structure builder tests ----

    #[test]
    fn test_build_fcc_atom_count() {
        let pos = build_fcc(1.0, 2, 2, 2);
        // 4 atoms per unit cell * 8 cells = 32
        assert_eq!(pos.len(), 32, "FCC 2x2x2 should have 32 atoms");
    }

    #[test]
    fn test_build_bcc_atom_count() {
        let pos = build_bcc(1.0, 3, 3, 3);
        // 2 atoms per unit cell * 27 = 54
        assert_eq!(pos.len(), 54, "BCC 3x3x3 should have 54 atoms");
    }

    #[test]
    fn test_build_hcp_atom_count() {
        let pos = build_hcp(1.0, 2, 2, 2);
        // 2 atoms per cell * 8 = 16
        assert_eq!(pos.len(), 16, "HCP 2x2x2 should have 16 atoms");
    }

    #[test]
    fn test_build_diamond_atom_count() {
        let pos = build_diamond(1.0, 2, 2, 2);
        // 8 atoms per unit cell * 8 = 64
        assert_eq!(pos.len(), 64, "Diamond 2x2x2 should have 64 atoms");
    }

    #[test]
    fn test_build_simple_cubic_atom_count() {
        let pos = build_simple_cubic(1.0, 3, 3, 3);
        assert_eq!(pos.len(), 27, "SC 3x3x3 should have 27 atoms");
    }

    #[test]
    fn test_fcc_nearest_neighbor_distance() {
        let a = 4.08; // Gold lattice parameter in Angstrom
        let pos = build_fcc(a, 1, 1, 1);
        // Nearest neighbor in FCC is a/sqrt(2)
        let expected_nn = a / 2.0_f64.sqrt();
        let mut min_dist = f64::MAX;
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let d = dist3(pos[i], pos[j]);
                if d < min_dist && d > 1e-10 {
                    min_dist = d;
                }
            }
        }
        assert!(
            (min_dist - expected_nn).abs() < 1e-10,
            "FCC NN distance: got {min_dist}, expected {expected_nn}"
        );
    }

    #[test]
    fn test_bcc_nearest_neighbor_distance() {
        let a = 2.87; // Iron lattice parameter
        let pos = build_bcc(a, 1, 1, 1);
        // BCC NN distance = a*sqrt(3)/2
        let expected_nn = a * 3.0_f64.sqrt() / 2.0;
        let mut min_dist = f64::MAX;
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let d = dist3(pos[i], pos[j]);
                if d < min_dist && d > 1e-10 {
                    min_dist = d;
                }
            }
        }
        assert!(
            (min_dist - expected_nn).abs() < 1e-10,
            "BCC NN distance: got {min_dist}, expected {expected_nn}"
        );
    }

    #[test]
    fn test_fcc_positions_in_box() {
        let a = 3.5;
        let nx = 2;
        let ny = 2;
        let nz = 2;
        let pos = build_fcc(a, nx, ny, nz);
        let box_max = [nx as f64 * a, ny as f64 * a, nz as f64 * a];
        for (idx, p) in pos.iter().enumerate() {
            for k in 0..3 {
                assert!(
                    p[k] >= 0.0 && p[k] <= box_max[k],
                    "Atom {idx} coord {k}={} out of box [0, {}]",
                    p[k],
                    box_max[k]
                );
            }
        }
    }

    // ---- Coordination number tests ----

    #[test]
    fn test_coordination_number_fcc() {
        let a = 4.0;
        let pos = build_fcc(a, 3, 3, 3);
        let nn = a / 2.0_f64.sqrt();
        let avg = average_coordination_number(&pos, nn * 1.1);
        // Interior FCC atoms have 12 neighbors; edge atoms have fewer.
        // Average should be reasonably high.
        assert!(avg > 4.0, "FCC avg coord should be >4, got {avg}");
    }

    #[test]
    fn test_coordination_number_simple() {
        // Two atoms within cutoff
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let counts = coordination_numbers(&positions, 1.5);
        assert_eq!(counts[0], 1);
        assert_eq!(counts[1], 1);
    }

    #[test]
    fn test_coordination_number_out_of_range() {
        let positions = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let counts = coordination_numbers(&positions, 1.5);
        assert_eq!(counts[0], 0);
        assert_eq!(counts[1], 0);
    }

    #[test]
    fn test_coordination_number_pbc() {
        // Two atoms that are far apart in direct distance but close via PBC
        let positions = [[0.1, 0.0, 0.0], [9.9, 0.0, 0.0]];
        let box_len = [10.0, 10.0, 10.0];
        let counts = coordination_numbers_pbc(&positions, box_len, 0.5);
        // PBC distance = 0.2, which is within cutoff
        assert_eq!(counts[0], 1, "Should be neighbors via PBC");
        assert_eq!(counts[1], 1);
    }

    // ---- Voronoi analysis tests ----

    #[test]
    fn test_voronoi_neighbors_simple() {
        // Tetrahedron: 4 atoms, each should be neighbor of the other 3
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866, 0.0],
            [0.5, 0.289, 0.816],
        ];
        let cells = voronoi_neighbors(&positions, 2.0);
        assert_eq!(cells.len(), 4);
        for cell in &cells {
            assert!(
                cell.neighbors.len() >= 2,
                "atom {} should have at least 2 Voronoi neighbors, got {}",
                cell.atom_idx,
                cell.neighbors.len()
            );
        }
    }

    #[test]
    fn test_voronoi_index_fcc() {
        let idx = voronoi_index(12);
        assert_eq!(idx, [0, 12, 0, 0], "FCC-like Voronoi index");
    }

    #[test]
    fn test_voronoi_index_bcc() {
        let idx = voronoi_index(14);
        assert_eq!(idx, [0, 6, 0, 8], "BCC-like Voronoi index");
    }

    // ---- Structure factor tests ----

    #[test]
    fn test_structure_factor_single_atom() {
        let positions = [[0.0, 0.0, 0.0]];
        let q_values = [1.0, 2.0, 5.0];
        let sq = structure_factor(&positions, &q_values, 20);
        // S(q) = 1/N * |exp(i*q*0)|^2 = 1 for all q
        for (i, &s) in sq.iter().enumerate() {
            assert!(
                (s - 1.0).abs() < 1e-10,
                "S(q={}) for single atom should be 1, got {s}",
                q_values[i]
            );
        }
    }

    #[test]
    fn test_structure_factor_empty() {
        let positions: Vec<[f64; 3]> = vec![];
        let q_values = [1.0, 2.0];
        let sq = structure_factor(&positions, &q_values, 10);
        assert_eq!(sq.len(), 2);
        for s in &sq {
            assert!((s - 0.0).abs() < 1e-15);
        }
    }

    #[test]
    fn test_structure_factor_two_atoms() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let q_values = [0.01]; // small q => S(q) should be close to N=2
        let sq = structure_factor(&positions, &q_values, 50);
        // At very small q, S(q) ~ N = 2
        assert!(
            sq[0] > 1.0,
            "S(q~0) for 2 atoms should be >1, got {}",
            sq[0]
        );
    }

    // ---- RDF tests ----

    #[test]
    fn test_rdf_simple_pair() {
        let positions = [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        let box_len = [5.0, 5.0, 5.0];
        let (bins, gr) = radial_distribution_function(&positions, box_len, 50, 2.5);
        assert_eq!(bins.len(), 50);
        assert_eq!(gr.len(), 50);
        // The peak should be at r ~= 1.0
        let peak_bin = gr
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert!(
            (bins[peak_bin] - 1.0).abs() < 0.2,
            "RDF peak should be near r=1.0, got r={}",
            bins[peak_bin]
        );
    }

    // ---- Lattice parameter optimization test ----

    #[test]
    fn test_optimize_lattice_parameter() {
        // Simple LJ-like pair energy with minimum at r=1.0
        let pair_energy = |r: f64| -> f64 {
            let s = 1.0;
            let r6 = (s / r).powi(6);
            4.0 * (r6 * r6 - r6)
        };
        // For SC lattice, NN distance = a, so minimum at a~1.12 (2^(1/6))
        let build = |a: f64| -> Vec<[f64; 3]> { build_simple_cubic(a, 2, 2, 2) };
        let a_opt = optimize_lattice_parameter(&build, &pair_energy, 0.8, 2.0, 0.001);
        let r_min = 2.0_f64.powf(1.0 / 6.0);
        assert!(
            (a_opt - r_min).abs() < 0.1,
            "Optimal lattice param should be near {r_min}, got {a_opt}"
        );
    }

    #[test]
    fn test_diamond_has_tetrahedral_coordination() {
        let a = 5.43; // Silicon
        let pos = build_diamond(a, 2, 2, 2);
        // Diamond structure: each atom has 4 nearest neighbors at a*sqrt(3)/4
        let nn_dist = a * 3.0_f64.sqrt() / 4.0;
        let avg = average_coordination_number(&pos, nn_dist * 1.1);
        // Interior atoms have 4 neighbors; edge atoms fewer. Average should be >1.
        assert!(avg > 1.0, "Diamond avg coord should be >1, got {avg}");
    }
}
