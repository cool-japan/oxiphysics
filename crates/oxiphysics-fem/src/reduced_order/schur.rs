// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real Schur decomposition and the Bartels–Stewart solution of the continuous
//! Lyapunov equation.
//!
//! Every real square matrix `A` can be written `A = Z T Zᵀ` with an orthogonal
//! `Z` and a **quasi-upper-triangular** `T`: block upper triangular with 1×1
//! blocks (real eigenvalues) and 2×2 blocks (complex conjugate pairs) on the
//! diagonal.  [`RealSchur`] computes it in three stages:
//!
//! 1. Householder reduction to upper Hessenberg form (`O(10/3 n³)`);
//! 2. the implicit double-shift QR iteration of Francis (the algorithm of
//!    EISPACK `hqr2`, Martin–Peters–Wilkinson 1970), which deflates 1×1 and
//!    2×2 blocks and splits 2×2 blocks with real eigenvalues, with exceptional
//!    shifts after 10 and 30 stagnant sweeps;
//! 3. an explicit clean-up of the entries that are zero by construction.
//!
//! With the factorisation, the Lyapunov equation
//!
//! ```text
//! A W + W Aᵀ + Q = 0
//! ```
//!
//! becomes the quasi-triangular equation `T Y + Y Tᵀ = −Zᵀ Q Z` for
//! `Y = Zᵀ W Z`, which the algorithm of Bartels and Stewart (1972) solves block
//! by block from the last block row and column upwards; each block of the
//! solution requires one Sylvester equation of at most 2×2 unknowns.  The total
//! cost is `O(n³)` and the algorithm is backward stable.

use super::{DenseMat, mat_mul, mat_transpose};
use crate::analysis::FemError;

/// Sweeps of the QR iteration allowed per eigenvalue (the sweeps counter is
/// reset whenever a block converges).
const QR_SWEEPS_PER_EIGENVALUE: usize = 60;

/// Real Schur decomposition `A = Z T Zᵀ` of a square matrix.
#[derive(Debug, Clone)]
pub(super) struct RealSchur {
    /// Quasi-upper-triangular factor `T`; everything below the first
    /// sub-diagonal and every sub-diagonal entry not belonging to a 2×2 block
    /// is exactly zero.
    pub(super) t: DenseMat,
    /// Orthogonal factor `Z`.
    pub(super) z: DenseMat,
}

/// Identity matrix of order `n`.
fn identity_matrix(n: usize) -> DenseMat {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}

/// Reduce `h` to upper Hessenberg form by Householder similarity
/// transformations and accumulate them into `z` (`z ← z P₁ P₂ …`), so that
/// `A = z h zᵀ` is preserved.
fn reduce_to_hessenberg(h: &mut [Vec<f64>], z: &mut [Vec<f64>]) {
    let n = h.len();
    for k in 0..n.saturating_sub(2) {
        // The reflector annihilates h[k+2.., k].
        if (k + 2..n).all(|i| h[i][k] == 0.0) {
            continue;
        }
        let scale: f64 = (k + 1..n).map(|i| h[i][k].abs()).sum();
        let mut v: Vec<f64> = (k + 1..n).map(|i| h[i][k] / scale).collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        // Sign chosen against cancellation in v[0] − alpha.
        let alpha = if v[0] > 0.0 { -norm } else { norm };
        v[0] -= alpha;
        let v_norm_sq: f64 = v.iter().map(|x| x * x).sum();
        if v_norm_sq == 0.0 {
            continue;
        }
        let factor = 2.0 / v_norm_sq;
        // h ← P h on the rows k+1.. (columns k..): h -= (2/vᵀv) v (vᵀ h).
        let mut projection = vec![0.0_f64; n - k];
        for (vt, row) in v.iter().zip(h[k + 1..].iter()) {
            for (sum, entry) in projection.iter_mut().zip(row[k..].iter()) {
                *sum += vt * entry;
            }
        }
        for (vt, row) in v.iter().zip(h[k + 1..].iter_mut()) {
            let weight = factor * vt;
            for (entry, sum) in row[k..].iter_mut().zip(projection.iter()) {
                *entry -= weight * sum;
            }
        }
        // h ← h P on the columns k+1.. (all rows), z ← z P likewise.
        for matrix in [&mut *h, &mut *z] {
            for row in matrix.iter_mut() {
                let dot: f64 = v
                    .iter()
                    .enumerate()
                    .map(|(t, vt)| row[k + 1 + t] * vt)
                    .sum();
                let f = factor * dot;
                for (t, vt) in v.iter().enumerate() {
                    row[k + 1 + t] -= f * vt;
                }
            }
        }
        h[k + 1][k] = alpha * scale;
        for row in h.iter_mut().skip(k + 2) {
            row[k] = 0.0;
        }
    }
}

/// Converge the 2×2 block in rows and columns `n-1..=n` of the Hessenberg
/// matrix `h`: a pair of real eigenvalues is rotated to upper triangular form
/// (accumulated into `z`), a complex pair is left as a 2×2 block.  The shifts
/// removed by exceptional-shift sweeps (`shift_total`) are added back to the
/// diagonal.
fn converge_pair(h: &mut [Vec<f64>], z: &mut [Vec<f64>], n: usize, shift_total: f64) {
    let w = h[n][n - 1] * h[n - 1][n];
    let mut p = (h[n - 1][n - 1] - h[n][n]) / 2.0;
    let mut q = p * p + w;
    let mut root = q.abs().sqrt();
    h[n][n] += shift_total;
    h[n - 1][n - 1] += shift_total;
    if q < 0.0 {
        return;
    }
    // Real pair: Givens rotation with the eigenvector of the first eigenvalue.
    root = if p >= 0.0 { p + root } else { p - root };
    let x = h[n][n - 1];
    let s = x.abs() + root.abs();
    p = x / s;
    q = root / s;
    let r = p.hypot(q);
    p /= r;
    q /= r;
    let (rows_above, rows_below) = h.split_at_mut(n);
    for (upper, lower) in rows_above[n - 1][n - 1..]
        .iter_mut()
        .zip(rows_below[0][n - 1..].iter_mut())
    {
        let kept = *upper;
        *upper = q * kept + p * *lower;
        *lower = q * *lower - p * kept;
    }
    for row in h.iter_mut().take(n + 1) {
        let left = row[n - 1];
        row[n - 1] = q * left + p * row[n];
        row[n] = q * row[n] - p * left;
    }
    for row in z.iter_mut() {
        let left = row[n - 1];
        row[n - 1] = q * left + p * row[n];
        row[n] = q * row[n] - p * left;
    }
    h[n][n - 1] = 0.0;
}

/// One implicit double-shift (Francis) QR sweep on the unreduced Hessenberg
/// block `l..=n` (`n ≥ l + 2`) of `h`, accumulated into `z`.
///
/// `sweeps` is the number of sweeps since the last convergence; after 10 and
/// 30 of them an exceptional shift breaks a stagnating iteration (the removed
/// diagonal shift is recorded in `shift_total`).
fn francis_sweep(
    h: &mut [Vec<f64>],
    z: &mut [Vec<f64>],
    l: usize,
    n: usize,
    sweeps: usize,
    shift_total: &mut f64,
) {
    let mut x = h[n][n];
    let mut y = h[n - 1][n - 1];
    let mut w = h[n][n - 1] * h[n - 1][n];
    if sweeps == 10 {
        // Wilkinson's original ad hoc shift.
        *shift_total += x;
        for (i, row) in h.iter_mut().enumerate().take(n + 1) {
            row[i] -= x;
        }
        let s = h[n][n - 1].abs() + h[n - 1][n - 2].abs();
        x = 0.75 * s;
        y = x;
        w = -0.4375 * s * s;
    }
    if sweeps == 30 {
        // The ad hoc shift of MATLAB.
        let half_gap = (y - x) / 2.0;
        let discriminant = half_gap * half_gap + w;
        if discriminant > 0.0 {
            let mut s = discriminant.sqrt();
            if y < x {
                s = -s;
            }
            s = x - w / (half_gap + s);
            for (i, row) in h.iter_mut().enumerate().take(n + 1) {
                row[i] -= s;
            }
            *shift_total += s;
            x = 0.964;
            y = x;
            w = x;
        }
    }

    // Look for two consecutive small sub-diagonal elements.
    let mut m = n - 2;
    let mut p: f64;
    let mut q: f64;
    let mut r: f64;
    loop {
        let diagonal = h[m][m];
        let r_shift = x - diagonal;
        let s_shift = y - diagonal;
        p = (r_shift * s_shift - w) / h[m + 1][m] + h[m][m + 1];
        q = h[m + 1][m + 1] - diagonal - r_shift - s_shift;
        r = h[m + 2][m + 1];
        let scale = p.abs() + q.abs() + r.abs();
        if scale != 0.0 {
            p /= scale;
            q /= scale;
            r /= scale;
        }
        if m == l {
            break;
        }
        let small = h[m][m - 1].abs() * (q.abs() + r.abs())
            < f64::EPSILON
                * (p.abs() * (h[m - 1][m - 1].abs() + diagonal.abs() + h[m + 1][m + 1].abs()));
        if small {
            break;
        }
        m -= 1;
    }
    for i in m + 2..=n {
        h[i][i - 2] = 0.0;
        if i > m + 2 {
            h[i][i - 3] = 0.0;
        }
    }

    // Double QR step on the rows l..=n and the columns m..=n.
    for k in m..n {
        let not_last = k != n - 1;
        if k != m {
            p = h[k][k - 1];
            q = h[k + 1][k - 1];
            r = if not_last { h[k + 2][k - 1] } else { 0.0 };
            x = p.abs() + q.abs() + r.abs();
            if x == 0.0 {
                continue;
            }
            p /= x;
            q /= x;
            r /= x;
        }
        let mut s = (p * p + q * q + r * r).sqrt();
        if p < 0.0 {
            s = -s;
        }
        if s == 0.0 {
            continue;
        }
        if k != m {
            h[k][k - 1] = -s * x;
        } else if l != m {
            h[k][k - 1] = -h[k][k - 1];
        }
        p += s;
        x = p / s;
        y = q / s;
        let third = r / s;
        q /= p;
        r /= p;

        let (rows_above, rows_below) = h.split_at_mut(k + 1);
        let (second, rows_below) = rows_below.split_at_mut(1);
        let first_row = &mut rows_above[k][k..];
        let second_row = &mut second[0][k..];
        if not_last {
            for ((a, b), c) in first_row
                .iter_mut()
                .zip(second_row.iter_mut())
                .zip(rows_below[0][k..].iter_mut())
            {
                let t = *a + q * *b + r * *c;
                *c -= t * third;
                *a -= t * x;
                *b -= t * y;
            }
        } else {
            for (a, b) in first_row.iter_mut().zip(second_row.iter_mut()) {
                let t = *a + q * *b;
                *a -= t * x;
                *b -= t * y;
            }
        }
        for row in h.iter_mut().take(n.min(k + 3) + 1) {
            let mut t = x * row[k] + y * row[k + 1];
            if not_last {
                t += third * row[k + 2];
                row[k + 2] -= t * r;
            }
            row[k] -= t;
            row[k + 1] -= t * q;
        }
        for row in z.iter_mut() {
            let mut t = x * row[k] + y * row[k + 1];
            if not_last {
                t += third * row[k + 2];
                row[k + 2] -= t * r;
            }
            row[k] -= t;
            row[k + 1] -= t * q;
        }
    }
}

/// Francis QR iteration on the Hessenberg matrix `h` until it is quasi upper
/// triangular; the transformations are accumulated into `z`.
///
/// # Errors
/// [`FemError::SolverNotConverged`] if a block does not converge within
/// [`QR_SWEEPS_PER_EIGENVALUE`] sweeps.
fn francis_qr_iteration(h: &mut [Vec<f64>], z: &mut [Vec<f64>]) -> Result<(), FemError> {
    let size = h.len();
    let norm: f64 = h
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row[i.saturating_sub(1)..]
                .iter()
                .map(|v| v.abs())
                .sum::<f64>()
        })
        .sum();
    let mut shift_total = 0.0_f64;
    let mut sweeps = 0_usize;
    // Rows and columns `active..` have converged.
    let mut active = size;
    while active > 0 {
        let n = active - 1;
        // Start of the unreduced block that ends in row n.
        let mut l = n;
        while l > 0 {
            let mut scale = h[l - 1][l - 1].abs() + h[l][l].abs();
            if scale == 0.0 {
                scale = norm;
            }
            if h[l][l - 1].abs() <= f64::EPSILON * scale {
                break;
            }
            l -= 1;
        }
        if l == n {
            // One real eigenvalue.
            h[n][n] += shift_total;
            if n > 0 {
                h[n][n - 1] = 0.0;
            }
            active -= 1;
            sweeps = 0;
        } else if l + 1 == n {
            // A pair of eigenvalues.
            converge_pair(h, z, n, shift_total);
            if n >= 2 {
                h[n - 1][n - 2] = 0.0;
            }
            active -= 2;
            sweeps = 0;
        } else {
            if sweeps >= QR_SWEEPS_PER_EIGENVALUE {
                return Err(FemError::SolverNotConverged);
            }
            francis_sweep(h, z, l, n, sweeps, &mut shift_total);
            sweeps += 1;
        }
    }
    // The bulge chasing leaves stale values below the sub-diagonal.
    for (i, row) in h.iter_mut().enumerate().skip(2) {
        for value in row.iter_mut().take(i - 1) {
            *value = 0.0;
        }
    }
    Ok(())
}

impl RealSchur {
    /// Compute the real Schur decomposition of the square matrix `a`.
    ///
    /// # Errors
    /// * [`FemError::InvalidInput`] if `a` is not square or contains a
    ///   non-finite entry.
    /// * [`FemError::SolverNotConverged`] if the QR iteration fails to converge
    ///   (not observed for finite matrices in practice).
    pub(super) fn new(a: &[Vec<f64>]) -> Result<Self, FemError> {
        let n = a.len();
        if a.iter().any(|row| row.len() != n) {
            return Err(FemError::InvalidInput(
                "the matrix of a Schur decomposition must be square".to_string(),
            ));
        }
        if a.iter().flatten().any(|v| !v.is_finite()) {
            return Err(FemError::InvalidInput(
                "the matrix of a Schur decomposition contains a non-finite entry".to_string(),
            ));
        }
        let mut t: DenseMat = a.to_vec();
        let mut z = identity_matrix(n);
        reduce_to_hessenberg(&mut t, &mut z);
        francis_qr_iteration(&mut t, &mut z)?;
        if t.iter().chain(z.iter()).flatten().any(|v| !v.is_finite()) {
            return Err(FemError::SolverNotConverged);
        }
        Ok(Self { t, z })
    }

    /// Start index and size (1 or 2) of the diagonal blocks of `T`.
    fn blocks(&self) -> Vec<(usize, usize)> {
        let n = self.t.len();
        let mut blocks = Vec::with_capacity(n);
        let mut start = 0;
        while start < n {
            let size = if start + 1 < n && self.t[start + 1][start] != 0.0 {
                2
            } else {
                1
            };
            blocks.push((start, size));
            start += size;
        }
        blocks
    }

    /// Eigenvalues `(real part, imaginary part)` of the decomposed matrix in
    /// the order of the diagonal blocks; a complex pair is listed with the
    /// positive imaginary part first.
    pub(super) fn eigenvalues(&self) -> Vec<(f64, f64)> {
        let mut eigenvalues = Vec::with_capacity(self.t.len());
        for (start, size) in self.blocks() {
            if size == 1 {
                eigenvalues.push((self.t[start][start], 0.0));
                continue;
            }
            let (a, b) = (self.t[start][start], self.t[start][start + 1]);
            let (c, d) = (self.t[start + 1][start], self.t[start + 1][start + 1]);
            let mean = 0.5 * (a + d);
            let half_gap = 0.5 * (a - d);
            let discriminant = half_gap * half_gap + b * c;
            if discriminant < 0.0 {
                let imaginary = (-discriminant).sqrt();
                eigenvalues.push((mean, imaginary));
                eigenvalues.push((mean, -imaginary));
            } else {
                let root = discriminant.sqrt();
                eigenvalues.push((mean + root, 0.0));
                eigenvalues.push((mean - root, 0.0));
            }
        }
        eigenvalues
    }

    /// Check that `A W + W Aᵀ + Q = 0` is uniquely solvable and `A` is
    /// asymptotically stable.
    ///
    /// # Errors
    /// * [`FemError::SingularSystem`] if two eigenvalues (counting an
    ///   eigenvalue with itself) add up to zero within rounding of the size of
    ///   the matrix: the operator `W ↦ A W + W Aᵀ` is singular.
    /// * [`FemError::InvalidInput`] if an eigenvalue has a non-negative real
    ///   part.
    fn check_lyapunov_spectrum(&self) -> Result<(), FemError> {
        let eigenvalues = self.eigenvalues();
        let scale = self.t.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
        let tolerance = 2.0 * eigenvalues.len() as f64 * f64::EPSILON * scale;
        for (i, &(re_i, im_i)) in eigenvalues.iter().enumerate() {
            for &(re_j, im_j) in &eigenvalues[i..] {
                if (re_i + re_j).hypot(im_i + im_j) <= tolerance {
                    return Err(FemError::SingularSystem);
                }
            }
        }
        if eigenvalues.iter().any(|&(re, _)| re >= 0.0) {
            return Err(FemError::InvalidInput(
                "A is not asymptotically stable (an eigenvalue has a non-negative real part)"
                    .to_string(),
            ));
        }
        Ok(())
    }

    /// Solve the Lyapunov equation `A W + W Aᵀ + Q = 0` for `W` with the
    /// Bartels–Stewart algorithm; `A` is the matrix this decomposition was
    /// computed from and `q` is an arbitrary `n × n` matrix.
    ///
    /// For symmetric positive semi-definite `Q` the solution is the symmetric
    /// positive semi-definite matrix `∫₀^∞ e^{At} Q e^{Aᵀt} dt` (up to
    /// rounding).
    ///
    /// # Errors
    /// * [`FemError::InvalidInput`] if `q` is not `n × n` and finite, or `A` is
    ///   not asymptotically stable.
    /// * [`FemError::SingularSystem`] if two eigenvalues of `A` add up to zero.
    pub(super) fn solve_lyapunov(&self, q: &[Vec<f64>]) -> Result<DenseMat, FemError> {
        let n = self.t.len();
        if q.len() != n || q.iter().any(|row| row.len() != n) {
            return Err(FemError::InvalidInput(format!(
                "the right-hand side of the Lyapunov equation must be a {n}x{n} matrix"
            )));
        }
        if q.iter().flatten().any(|v| !v.is_finite()) {
            return Err(FemError::InvalidInput(
                "the right-hand side of the Lyapunov equation contains a non-finite entry"
                    .to_string(),
            ));
        }
        if n == 0 {
            return Ok(Vec::new());
        }
        self.check_lyapunov_spectrum()?;
        let z_transpose = mat_transpose(&self.z);
        // T Y + Y Tᵀ = −Zᵀ Q Z
        let mut y = mat_mul(&mat_mul(&z_transpose, q), &self.z);
        for value in y.iter_mut().flatten() {
            *value = -*value;
        }
        solve_quasi_triangular(&self.t, &self.blocks(), &mut y)?;
        let w = mat_mul(&mat_mul(&self.z, &y), &z_transpose);
        if w.iter().flatten().any(|v| !v.is_finite()) {
            return Err(FemError::SingularSystem);
        }
        Ok(w)
    }
}

/// Solve `T Y + Y Tᵀ = C` for `Y` in place (`y` holds `C` on entry) with
/// `T` quasi upper triangular and `blocks` its diagonal block structure.
///
/// The blocks `Y_IJ` are computed from the last block row and column upwards:
///
/// ```text
/// T_II Y_IJ + Y_IJ T_JJᵀ = C_IJ − Σ_{K>I} T_IK Y_KJ − Σ_{K>J} Y_IK T_JKᵀ.
/// ```
fn solve_quasi_triangular(
    t: &[Vec<f64>],
    blocks: &[(usize, usize)],
    y: &mut [Vec<f64>],
) -> Result<(), FemError> {
    let n = t.len();
    for &(i0, p) in blocks.iter().rev() {
        let i1 = i0 + p;
        for &(j0, q) in blocks.iter().rev() {
            let j1 = j0 + q;
            let mut rhs = [[0.0_f64; 2]; 2];
            for (r, rhs_row) in rhs.iter_mut().enumerate().take(p) {
                for (s, entry) in rhs_row.iter_mut().enumerate().take(q) {
                    let mut value = y[i0 + r][j0 + s];
                    for k in i1..n {
                        value -= t[i0 + r][k] * y[k][j0 + s];
                    }
                    for k in j1..n {
                        value -= y[i0 + r][k] * t[j0 + s][k];
                    }
                    *entry = value;
                }
            }
            let diagonal_i = block_of(t, i0, p);
            let diagonal_j = block_of(t, j0, q);
            let solution = solve_small_sylvester(&diagonal_i, p, &diagonal_j, q, &rhs)?;
            for (r, solution_row) in solution.iter().enumerate().take(p) {
                for (s, &value) in solution_row.iter().enumerate().take(q) {
                    y[i0 + r][j0 + s] = value;
                }
            }
        }
    }
    Ok(())
}

/// The `size × size` diagonal block of `t` starting at `start`, padded with
/// zeros to 2×2.
fn block_of(t: &[Vec<f64>], start: usize, size: usize) -> [[f64; 2]; 2] {
    let mut block = [[0.0_f64; 2]; 2];
    for (r, block_row) in block.iter_mut().enumerate().take(size) {
        for (c, entry) in block_row.iter_mut().enumerate().take(size) {
            *entry = t[start + r][start + c];
        }
    }
    block
}

/// Solve `A X + X Bᵀ = R` for the `p × q` matrix `X` (`p, q ∈ {1, 2}`), where
/// `a` and `b` hold the `p × p` and `q × q` matrices in their leading corners.
///
/// # Errors
/// [`FemError::SingularSystem`] if the `pq × pq` system (the Kronecker sum
/// `I⊗A + B⊗I`) is singular to working precision.
fn solve_small_sylvester(
    a: &[[f64; 2]; 2],
    p: usize,
    b: &[[f64; 2]; 2],
    q: usize,
    rhs: &[[f64; 2]; 2],
) -> Result<[[f64; 2]; 2], FemError> {
    let size = p * q;
    // Augmented matrix [K | r] with the unknowns ordered row-major.
    let mut system = [[0.0_f64; 5]; 4];
    for r in 0..p {
        for s in 0..q {
            let row = r * q + s;
            for l in 0..p {
                system[row][l * q + s] += a[r][l];
            }
            for l in 0..q {
                system[row][r * q + l] += b[s][l];
            }
            system[row][size] = rhs[r][s];
        }
    }
    let scale = system
        .iter()
        .take(size)
        .flat_map(|row| row.iter().take(size))
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    let tolerance = size as f64 * f64::EPSILON * scale;
    for col in 0..size {
        let mut pivot_row = col;
        for row in col + 1..size {
            if system[row][col].abs() > system[pivot_row][col].abs() {
                pivot_row = row;
            }
        }
        if system[pivot_row][col].abs() <= tolerance {
            return Err(FemError::SingularSystem);
        }
        system.swap(col, pivot_row);
        let pivot_values = system[col];
        for target_row in system.iter_mut().take(size).skip(col + 1) {
            let factor = target_row[col] / pivot_values[col];
            if factor != 0.0 {
                for (target, source) in target_row[col..=size]
                    .iter_mut()
                    .zip(pivot_values[col..=size].iter())
                {
                    *target -= factor * source;
                }
            }
        }
    }
    let mut unknowns = [0.0_f64; 4];
    for row in (0..size).rev() {
        let mut value = system[row][size];
        for c in row + 1..size {
            value -= system[row][c] * unknowns[c];
        }
        unknowns[row] = value / system[row][row];
    }
    let mut solution = [[0.0_f64; 2]; 2];
    for r in 0..p {
        for s in 0..q {
            solution[r][s] = unknowns[r * q + s];
        }
    }
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random value in (-1, 1).
    fn lcg_unit(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn random_matrix(rows: usize, cols: usize, seed: u64) -> DenseMat {
        let mut state = seed;
        (0..rows)
            .map(|_| (0..cols).map(|_| lcg_unit(&mut state)).collect())
            .collect()
    }

    fn frobenius(m: &[Vec<f64>]) -> f64 {
        m.iter().flatten().map(|v| v * v).sum::<f64>().sqrt()
    }

    fn trace(m: &[Vec<f64>]) -> f64 {
        m.iter().enumerate().map(|(i, row)| row[i]).sum()
    }

    /// `‖Z T Zᵀ − A‖_F`.
    fn reconstruction_error(a: &[Vec<f64>], schur: &RealSchur) -> f64 {
        let rebuilt = mat_mul(&mat_mul(&schur.z, &schur.t), &mat_transpose(&schur.z));
        let difference: DenseMat = rebuilt
            .iter()
            .zip(a.iter())
            .map(|(x, y)| x.iter().zip(y.iter()).map(|(u, v)| u - v).collect())
            .collect();
        frobenius(&difference)
    }

    /// `‖ZᵀZ − I‖_F`.
    fn orthogonality_error(schur: &RealSchur) -> f64 {
        let product = mat_mul(&mat_transpose(&schur.z), &schur.z);
        let n = product.len();
        let difference: DenseMat = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| product[i][j] - if i == j { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();
        frobenius(&difference)
    }

    /// Assert the structure promised by the decomposition: zeros below the
    /// first sub-diagonal, no two consecutive non-zero sub-diagonal entries and
    /// complex eigenvalues in every 2×2 block.
    fn assert_quasi_triangular(schur: &RealSchur) {
        let t = &schur.t;
        let n = t.len();
        for (i, row) in t.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                if i > j + 1 {
                    assert_eq!(value, 0.0, "T[{i}][{j}] must be exactly zero");
                }
            }
        }
        for i in 0..n.saturating_sub(2) {
            assert!(
                t[i + 1][i] == 0.0 || t[i + 2][i + 1] == 0.0,
                "sub-diagonal entries {i} and {} are both non-zero",
                i + 1
            );
        }
        for (start, size) in schur.blocks() {
            if size == 2 {
                let (a, b) = (t[start][start], t[start][start + 1]);
                let (c, d) = (t[start + 1][start], t[start + 1][start + 1]);
                let gap = 0.5 * (a - d);
                assert!(
                    gap * gap + b * c < 0.0,
                    "2x2 block at {start} has real roots"
                );
            }
        }
    }

    fn assert_valid_decomposition(a: &[Vec<f64>], tolerance: f64) -> RealSchur {
        let schur = RealSchur::new(a).expect("decomposition");
        let n = a.len() as f64;
        let size = frobenius(a).max(1.0);
        let error = reconstruction_error(a, &schur);
        assert!(
            error <= tolerance * n * size,
            "n = {n}: reconstruction error {error}"
        );
        assert!(
            orthogonality_error(&schur) <= tolerance * n,
            "n = {n}: Z is not orthogonal"
        );
        assert_quasi_triangular(&schur);
        // Trace identities: Σλ = tr A and Σλ² = tr A².
        let eigenvalues = schur.eigenvalues();
        assert_eq!(eigenvalues.len(), a.len());
        let sum: f64 = eigenvalues.iter().map(|(re, _)| re).sum();
        assert!((sum - trace(a)).abs() <= tolerance * n * size, "tr A");
        let sum_squares: f64 = eigenvalues.iter().map(|(re, im)| re * re - im * im).sum();
        let a_squared = mat_mul(a, a);
        assert!(
            (sum_squares - trace(&a_squared)).abs() <= tolerance * n * size * size,
            "tr A²"
        );
        schur
    }

    #[test]
    fn random_matrices_of_all_sizes_decompose() {
        for (index, &n) in [1usize, 2, 3, 4, 5, 6, 8, 11, 13, 21, 34, 55, 89]
            .iter()
            .enumerate()
        {
            let a = random_matrix(n, n, 100 + index as u64);
            assert_valid_decomposition(&a, 1e-13);
        }
    }

    #[test]
    fn badly_scaled_matrices_decompose() {
        let n = 12;
        let mut a = random_matrix(n, n, 7);
        for (i, row) in a.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value *= 10f64.powi(i as i32 - j as i32);
            }
        }
        assert_valid_decomposition(&a, 1e-11);
    }

    #[test]
    fn symmetric_matrices_have_a_diagonal_schur_form_with_known_spectrum() {
        // tridiag(1, -2, 1): λ_k = -2 + 2 cos(kπ/(n+1)), k = 1..n.
        let n = 17;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            a[i][i] = -2.0;
            if i + 1 < n {
                a[i][i + 1] = 1.0;
                a[i + 1][i] = 1.0;
            }
        }
        let schur = assert_valid_decomposition(&a, 1e-13);
        let mut computed: Vec<f64> = schur.eigenvalues().iter().map(|(re, _)| *re).collect();
        assert!(schur.eigenvalues().iter().all(|(_, im)| *im == 0.0));
        computed.sort_by(f64::total_cmp);
        let mut exact: Vec<f64> = (1..=n)
            .map(|k| -2.0 + 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos())
            .collect();
        exact.sort_by(f64::total_cmp);
        for (c, e) in computed.iter().zip(exact.iter()) {
            assert!((c - e).abs() < 1e-13, "{c} vs {e}");
        }
        for (i, row) in schur.t.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                if i != j {
                    assert!(value.abs() < 1e-13, "T[{i}][{j}] = {value}");
                }
            }
        }
    }

    #[test]
    fn companion_matrix_recovers_the_roots() {
        // (s+1)(s+2)(s+3)(s+4) = s⁴ + 10 s³ + 35 s² + 50 s + 24.
        let a = vec![
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
            vec![-24.0, -50.0, -35.0, -10.0],
        ];
        let schur = assert_valid_decomposition(&a, 1e-13);
        let mut roots: Vec<f64> = schur.eigenvalues().iter().map(|(re, _)| *re).collect();
        roots.sort_by(f64::total_cmp);
        for (root, exact) in roots.iter().zip([-4.0, -3.0, -2.0, -1.0]) {
            assert!((root - exact).abs() < 1e-10, "{roots:?}");
        }
    }

    #[test]
    fn complex_pairs_are_kept_as_two_by_two_blocks() {
        // Block diagonal: damped oscillators -0.5 ± 3i and -0.1 ± 7i, and a real pole -2.
        let a = vec![
            vec![-0.5, 3.0, 0.0, 0.0, 0.0],
            vec![-3.0, -0.5, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, -0.1, 7.0, 0.0],
            vec![0.0, 0.0, -7.0, -0.1, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, -2.0],
        ];
        let schur = assert_valid_decomposition(&a, 1e-13);
        let mut found: Vec<(f64, f64)> = schur
            .eigenvalues()
            .into_iter()
            .map(|(re, im)| (re, im.abs()))
            .collect();
        found.sort_by(|x, y| x.1.total_cmp(&y.1).then(x.0.total_cmp(&y.0)));
        let expected = [
            (-2.0, 0.0),
            (-0.5, 3.0),
            (-0.5, 3.0),
            (-0.1, 7.0),
            (-0.1, 7.0),
        ];
        for ((re, im), (exact_re, exact_im)) in found.iter().zip(expected) {
            assert!((re - exact_re).abs() < 1e-12 && (im - exact_im).abs() < 1e-12);
        }
        assert_eq!(
            schur.blocks().iter().filter(|(_, size)| *size == 2).count(),
            2
        );
    }

    #[test]
    fn cyclic_permutations_and_rotations_converge() {
        // The cyclic shift has the n-th roots of unity as eigenvalues and defeats an
        // unshifted QR iteration.
        for n in [2usize, 3, 4, 5, 7, 8, 12] {
            let mut a = vec![vec![0.0; n]; n];
            for i in 0..n {
                a[(i + 1) % n][i] = 1.0;
            }
            let schur = assert_valid_decomposition(&a, 1e-13);
            for (re, im) in schur.eigenvalues() {
                assert!((re.hypot(im) - 1.0).abs() < 1e-12, "|λ| = {}", re.hypot(im));
            }
        }
        let rotation = vec![vec![0.0, 1.0], vec![-1.0, 0.0]];
        let schur = assert_valid_decomposition(&rotation, 1e-15);
        let eigenvalues = schur.eigenvalues();
        assert!(
            eigenvalues
                .iter()
                .all(|(re, im)| re.abs() < 1e-15 && (im.abs() - 1.0).abs() < 1e-15)
        );
    }

    #[test]
    fn defective_and_degenerate_matrices_decompose() {
        let zero = vec![vec![0.0; 4]; 4];
        let schur = assert_valid_decomposition(&zero, 1e-15);
        assert!(
            schur
                .eigenvalues()
                .iter()
                .all(|(re, im)| *re == 0.0 && *im == 0.0)
        );
        let identity = identity_matrix(5);
        assert_valid_decomposition(&identity, 1e-15);
        // Jordan block with the triple eigenvalue -1.
        let jordan = vec![
            vec![-1.0, 1.0, 0.0],
            vec![0.0, -1.0, 1.0],
            vec![0.0, 0.0, -1.0],
        ];
        let schur = assert_valid_decomposition(&jordan, 1e-15);
        assert!(
            schur
                .eigenvalues()
                .iter()
                .all(|(re, _)| (re + 1.0).abs() < 1e-6)
        );
        // Lower triangular and rank-one matrices.
        let mut lower = random_matrix(6, 6, 33);
        for (i, row) in lower.iter_mut().enumerate() {
            for value in row.iter_mut().skip(i + 1) {
                *value = 0.0;
            }
        }
        assert_valid_decomposition(&lower, 1e-13);
        let u = random_matrix(7, 1, 34);
        let v = random_matrix(1, 7, 35);
        let rank_one = mat_mul(&u, &v);
        let schur = assert_valid_decomposition(&rank_one, 1e-13);
        let nonzero = schur
            .eigenvalues()
            .iter()
            .filter(|(re, im)| re.hypot(*im) > 1e-12)
            .count();
        assert_eq!(nonzero, 1);
    }

    #[test]
    fn state_matrix_of_a_lightly_damped_structure_has_stable_complex_poles() {
        // x = (q, q'), M = I: A = [[0, I], [-K, -C]], K = tridiag(-1, 2, -1) * 100, C = 1e-3 K.
        let m = 12;
        let n = 2 * m;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..m {
            a[i][m + i] = 1.0;
            for (j, value) in [(i.wrapping_sub(1), -100.0), (i, 200.0), (i + 1, -100.0)] {
                if j < m {
                    a[m + i][j] = -value;
                    a[m + i][m + j] = -1e-3 * value;
                }
            }
        }
        let schur = assert_valid_decomposition(&a, 1e-12);
        let eigenvalues = schur.eigenvalues();
        assert!(eigenvalues.iter().all(|(re, _)| *re < 0.0));
        assert_eq!(
            schur.blocks().iter().filter(|(_, size)| *size == 2).count(),
            m
        );
        // λ² + c_k λ + k_k = 0 with k_k = 400 (1 - cos(kπ/(m+1))): |Im λ| ≈ sqrt(k_k).
        let mut frequencies: Vec<f64> = eigenvalues
            .iter()
            .filter(|(_, im)| *im > 0.0)
            .map(|(_, im)| *im)
            .collect();
        frequencies.sort_by(f64::total_cmp);
        for (k, frequency) in frequencies.iter().enumerate() {
            let stiffness =
                200.0 - 200.0 * ((k + 1) as f64 * std::f64::consts::PI / (m as f64 + 1.0)).cos();
            let exact = (stiffness * (1.0 - 0.25e-6 * stiffness)).sqrt();
            assert!(
                (frequency - exact).abs() < 1e-8 * exact,
                "{frequency} vs {exact}"
            );
        }
    }

    #[test]
    fn empty_and_malformed_matrices_are_handled() {
        let empty = RealSchur::new(&[]).expect("empty");
        assert!(empty.t.is_empty() && empty.z.is_empty());
        assert!(empty.eigenvalues().is_empty());
        assert!(matches!(
            RealSchur::new(&[vec![1.0, 2.0]]),
            Err(FemError::InvalidInput(_))
        ));
        assert!(matches!(
            RealSchur::new(&[vec![1.0, f64::NAN], vec![0.0, 1.0]]),
            Err(FemError::InvalidInput(_))
        ));
        let one = RealSchur::new(&[vec![-3.0]]).expect("scalar");
        assert_eq!(one.t, vec![vec![-3.0]]);
        assert_eq!(one.z, vec![vec![1.0]]);
    }

    // ---- small Sylvester equations ----------------------------------------------

    #[test]
    fn small_sylvester_equations_are_solved_for_every_block_size() {
        let a = [[-1.0, 2.0], [-3.0, -0.5]];
        let b = [[-2.0, 1.0], [-1.5, -0.25]];
        let x_true = [[1.0, -2.0], [0.5, 3.0]];
        for (p, q) in [(1usize, 1usize), (1, 2), (2, 1), (2, 2)] {
            let mut rhs = [[0.0_f64; 2]; 2];
            for r in 0..p {
                for s in 0..q {
                    let ax: f64 = (0..p).map(|l| a[r][l] * x_true[l][s]).sum();
                    let xb: f64 = (0..q).map(|l| x_true[r][l] * b[s][l]).sum();
                    rhs[r][s] = ax + xb;
                }
            }
            let x = solve_small_sylvester(&a, p, &b, q, &rhs).expect("regular");
            for r in 0..p {
                for s in 0..q {
                    assert!((x[r][s] - x_true[r][s]).abs() < 1e-13, "({p}, {q})");
                }
            }
        }
        // λ_a + λ_b = 0 makes the operator singular.
        let singular = solve_small_sylvester(
            &[[1.0, 0.0], [0.0, 0.0]],
            1,
            &[[-1.0, 0.0], [0.0, 0.0]],
            1,
            &[[1.0, 0.0], [0.0, 0.0]],
        );
        assert!(matches!(singular, Err(FemError::SingularSystem)));
    }
}
