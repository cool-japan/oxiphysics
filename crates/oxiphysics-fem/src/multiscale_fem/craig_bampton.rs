// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Craig-Bampton component mode synthesis (fixed-interface substructuring).
//!
//! Every substructure (component) with stiffness `K` and mass `M` is split into
//! boundary DOFs `b` (shared with the rest of the structure) and interior DOFs
//! `i`.  The component displacement is represented by the Craig-Bampton basis
//!
//! ```text
//!     ⎡u_b⎤   ⎡ I  0 ⎤ ⎡u_b⎤
//!     ⎣u_i⎦ = ⎣ Ψ  Φ ⎦ ⎣ η ⎦ = T q,      Ψ = −K_ii⁻¹ K_ib,
//! ```
//!
//! where `Ψ` are the static constraint modes and `Φ` the lowest mass-normalised
//! fixed-interface normal modes (`K_ii φ = ω² M_ii φ`, `φᵀ M_ii φ = 1`).  The
//! reduced matrices are `K_r = Tᵀ K T`, `M_r = Tᵀ M T`.  Components share the
//! same interface DOFs; the assembled reduced system is ordered
//! `[interface DOFs, modal DOFs of component 0, modal DOFs of component 1, …]`.

use crate::eigen::jacobi_eigen_dense;
use crate::error::{Error, Result};
use crate::solvers::CholeskyDense;

/// Relative pivot threshold below which a Cholesky factorisation is treated as
/// numerically singular.
const PIVOT_RTOL: f64 = 1e-12;

/// Relative tolerance for the symmetry check of input matrices.
const SYMMETRY_RTOL: f64 = 1e-9;

/// One substructure for [`SubstructureMethod::from_substructures`]:
/// `(stiffness, mass, boundary_dofs)`.
pub type SubstructureInput<'a> = (&'a [Vec<f64>], &'a [Vec<f64>], &'a [usize]);

/// One Craig-Bampton reduced component.
#[derive(Debug, Clone)]
pub struct CraigBamptonComponent {
    /// Local indices of the boundary (interface) DOFs, in interface order.
    pub boundary_dofs: Vec<usize>,
    /// Local indices of the interior DOFs (all DOFs that are not boundary DOFs).
    pub interior_dofs: Vec<usize>,
    /// Static constraint modes `Ψ = −K_ii⁻¹ K_ib` (`n_interior × n_boundary`).
    pub constraint_modes: Vec<Vec<f64>>,
    /// Mass-normalised fixed-interface normal modes; `normal_modes[k]` is the
    /// k-th mode over the interior DOFs (length `n_interior`).
    pub normal_modes: Vec<Vec<f64>>,
    /// Fixed-interface natural circular frequencies ωₖ \[rad/s\] (ascending).
    pub fixed_interface_freqs: Vec<f64>,
    /// Component reduced stiffness `Tᵀ K T`, ordered `[boundary, modes]`.
    pub reduced_stiffness: Vec<Vec<f64>>,
    /// Component reduced mass `Tᵀ M T`, ordered `[boundary, modes]`.
    pub reduced_mass: Vec<Vec<f64>>,
}

impl CraigBamptonComponent {
    /// Reduce one component given its full stiffness `k`, mass `m`, the local
    /// indices of its boundary DOFs and the number of fixed-interface modes to keep.
    ///
    /// # Errors
    /// * [`Error::InvalidInput`] for non-square / mismatched / asymmetric
    ///   matrices, duplicate boundary DOFs or `n_modes > n_interior`.
    /// * [`Error::DofOutOfRange`] for a boundary DOF outside the matrix.
    /// * [`Error::SingularSystem`] if `K_ii` or `M_ii` is not positive definite
    ///   (the component is not restrained by its boundary DOFs).
    pub fn new(
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        boundary_dofs: &[usize],
        n_modes: usize,
    ) -> Result<Self> {
        let n = k.len();
        validate_square_symmetric(k, n, "stiffness")?;
        validate_square_symmetric(m, n, "mass")?;
        let mut is_boundary = vec![false; n];
        for &b in boundary_dofs {
            if b >= n {
                return Err(Error::DofOutOfRange { dof: b, n });
            }
            if is_boundary[b] {
                return Err(Error::InvalidInput(format!(
                    "boundary DOF {b} is listed more than once"
                )));
            }
            is_boundary[b] = true;
        }
        let interior_dofs: Vec<usize> = (0..n).filter(|&d| !is_boundary[d]).collect();
        let n_b = boundary_dofs.len();
        let n_i = interior_dofs.len();
        if n_modes > n_i {
            return Err(Error::InvalidInput(format!(
                "requested {n_modes} fixed-interface modes but the component has only {n_i} interior DOFs"
            )));
        }

        let k_ii = sub_matrix(k, &interior_dofs, &interior_dofs);
        let m_ii = sub_matrix(m, &interior_dofs, &interior_dofs);

        // Static constraint modes Ψ = −K_ii⁻¹ K_ib (column j = response to u_b = e_j).
        let mut constraint_modes = vec![vec![0.0; n_b]; n_i];
        if n_i > 0 {
            let chol = CholeskyDense::new(n_i);
            let l_kii = cholesky_checked(&chol, &k_ii)?;
            for (j, &b) in boundary_dofs.iter().enumerate() {
                let rhs: Vec<f64> = interior_dofs.iter().map(|&i| -k[i][b]).collect();
                let col = chol.solve(&l_kii, &rhs);
                for (row, &v) in constraint_modes.iter_mut().zip(col.iter()) {
                    row[j] = v;
                }
            }
        }

        // Fixed-interface normal modes.
        let (normal_modes, fixed_interface_freqs) = if n_i > 0 && n_modes > 0 {
            let (lambdas, modes) = generalized_symmetric_eigen(&k_ii, &m_ii)?;
            let freqs: Vec<f64> = lambdas
                .iter()
                .take(n_modes)
                .map(|&l| l.max(0.0).sqrt())
                .collect();
            let kept: Vec<Vec<f64>> = modes.into_iter().take(n_modes).collect();
            (kept, freqs)
        } else {
            (Vec::new(), Vec::new())
        };

        // Transformation T (n × r) in the component's original DOF ordering.
        let r = n_b + n_modes;
        let mut t = vec![vec![0.0; r]; n];
        for (j, &b) in boundary_dofs.iter().enumerate() {
            t[b][j] = 1.0;
        }
        for (row_i, &i) in interior_dofs.iter().enumerate() {
            t[i][..n_b].copy_from_slice(&constraint_modes[row_i]);
            for (mode_k, phi) in normal_modes.iter().enumerate() {
                t[i][n_b + mode_k] = phi[row_i];
            }
        }
        let reduced_stiffness = congruence(k, &t);
        let reduced_mass = congruence(m, &t);

        Ok(Self {
            boundary_dofs: boundary_dofs.to_vec(),
            interior_dofs,
            constraint_modes,
            normal_modes,
            fixed_interface_freqs,
            reduced_stiffness,
            reduced_mass,
        })
    }

    /// Number of boundary DOFs.
    pub fn n_boundary(&self) -> usize {
        self.boundary_dofs.len()
    }

    /// Number of retained fixed-interface modes.
    pub fn n_modes(&self) -> usize {
        self.normal_modes.len()
    }

    /// Recover the full component displacement from boundary displacements
    /// `u_b` and modal amplitudes `eta`: `u_b` on boundary DOFs and
    /// `u_i = Ψ u_b + Φ η` on interior DOFs.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the slice lengths do not match the component.
    pub fn expand(&self, u_b: &[f64], eta: &[f64]) -> Result<Vec<f64>> {
        if u_b.len() != self.n_boundary() || eta.len() != self.n_modes() {
            return Err(Error::InvalidInput(format!(
                "expected {} boundary values and {} modal amplitudes, got {} and {}",
                self.n_boundary(),
                self.n_modes(),
                u_b.len(),
                eta.len()
            )));
        }
        let n = self.boundary_dofs.len() + self.interior_dofs.len();
        let mut u = vec![0.0; n];
        for (&b, &v) in self.boundary_dofs.iter().zip(u_b.iter()) {
            u[b] = v;
        }
        for (row_i, &i) in self.interior_dofs.iter().enumerate() {
            let static_part: f64 = self.constraint_modes[row_i]
                .iter()
                .zip(u_b.iter())
                .map(|(psi, ub)| psi * ub)
                .sum();
            let modal_part: f64 = self
                .normal_modes
                .iter()
                .zip(eta.iter())
                .map(|(phi, e)| phi[row_i] * e)
                .sum();
            u[i] = static_part + modal_part;
        }
        Ok(u)
    }
}

/// Craig-Bampton substructuring for large structural models.
///
/// Reduces each substructure to the shared interface DOFs plus a small set of
/// fixed-interface normal modes, and assembles the coupled reduced system.
#[derive(Debug, Clone)]
pub struct SubstructureMethod {
    /// Number of interface (boundary) DOFs shared by all substructures.
    pub n_interface: usize,
    /// Number of kept interior modes per substructure.
    pub n_modes: usize,
    /// Number of substructures.
    pub n_sub: usize,
    /// Fixed-interface modal circular frequencies (rad/s) for each substructure.
    pub modal_freqs: Vec<Vec<f64>>,
    /// Reduced components (constraint modes, normal modes, reduced matrices).
    pub components: Vec<CraigBamptonComponent>,
    /// Assembled reduced stiffness (`total_reduced_dof²`).
    pub reduced_stiffness_matrix: Vec<Vec<f64>>,
    /// Assembled reduced mass (`total_reduced_dof²`).
    pub reduced_mass_matrix: Vec<Vec<f64>>,
}

impl SubstructureMethod {
    /// Craig-Bampton reduction of a single structure.
    ///
    /// # Arguments
    /// * `k`, `m`         – full symmetric stiffness and mass matrices
    /// * `boundary_dofs`  – DOFs retained physically (interface / attachment DOFs)
    /// * `n_modes`        – number of fixed-interface modes to keep
    ///
    /// # Errors
    /// See [`CraigBamptonComponent::new`].
    pub fn from_matrices(
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        boundary_dofs: &[usize],
        n_modes: usize,
    ) -> Result<Self> {
        let component = CraigBamptonComponent::new(k, m, boundary_dofs, n_modes)?;
        Self::from_components(vec![component])
    }

    /// Reduce several substructures that share the same interface DOFs and
    /// assemble the coupled reduced system.
    ///
    /// Each entry is `(k, m, boundary_dofs)`; `boundary_dofs[j]` is the local DOF
    /// of that substructure attached to global interface DOF `j`, so all
    /// substructures must list the same number of boundary DOFs.
    ///
    /// # Errors
    /// See [`CraigBamptonComponent::new`]; additionally [`Error::InvalidInput`]
    /// if the substructure list is empty or the interface sizes disagree.
    pub fn from_substructures(
        substructures: &[SubstructureInput<'_>],
        n_modes: usize,
    ) -> Result<Self> {
        let components = substructures
            .iter()
            .map(|&(k, m, b)| CraigBamptonComponent::new(k, m, b, n_modes))
            .collect::<Result<Vec<_>>>()?;
        Self::from_components(components)
    }

    /// Assemble already reduced components sharing one interface.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the list is empty, or the components disagree
    /// on the number of interface DOFs or retained modes.
    pub fn from_components(components: Vec<CraigBamptonComponent>) -> Result<Self> {
        let first = components.first().ok_or_else(|| {
            Error::InvalidInput("Craig-Bampton assembly needs at least one substructure".into())
        })?;
        let n_interface = first.n_boundary();
        let n_modes = first.n_modes();
        if components
            .iter()
            .any(|c| c.n_boundary() != n_interface || c.n_modes() != n_modes)
        {
            return Err(Error::InvalidInput(
                "all substructures must share the same interface size and retained mode count"
                    .into(),
            ));
        }
        let n_sub = components.len();
        let n_red = n_interface + n_sub * n_modes;
        let mut k_red = vec![vec![0.0; n_red]; n_red];
        let mut m_red = vec![vec![0.0; n_red]; n_red];
        for (s, comp) in components.iter().enumerate() {
            // Local reduced index → global reduced index.
            let map: Vec<usize> = (0..n_interface)
                .chain((0..n_modes).map(|k| n_interface + s * n_modes + k))
                .collect();
            for (a, &ga) in map.iter().enumerate() {
                for (b, &gb) in map.iter().enumerate() {
                    k_red[ga][gb] += comp.reduced_stiffness[a][b];
                    m_red[ga][gb] += comp.reduced_mass[a][b];
                }
            }
        }
        let modal_freqs = components
            .iter()
            .map(|c| c.fixed_interface_freqs.clone())
            .collect();
        Ok(Self {
            n_interface,
            n_modes,
            n_sub,
            modal_freqs,
            components,
            reduced_stiffness_matrix: k_red,
            reduced_mass_matrix: m_red,
        })
    }

    /// Total reduced DOF count per substructure: n_interface + n_modes.
    pub fn reduced_dof(&self) -> usize {
        self.n_interface + self.n_modes
    }

    /// Total system DOF after Craig-Bampton reduction.
    pub fn total_reduced_dof(&self) -> usize {
        self.n_interface + self.n_sub * self.n_modes
    }

    /// Assembled reduced mass matrix `Σ Tᵀ M T` (`total_reduced_dof²`).
    pub fn reduced_mass(&self) -> &[Vec<f64>] {
        &self.reduced_mass_matrix
    }

    /// Assembled reduced stiffness matrix `Σ Tᵀ K T` (`total_reduced_dof²`).
    pub fn reduced_stiffness(&self) -> &[Vec<f64>] {
        &self.reduced_stiffness_matrix
    }

    /// Natural circular frequencies (rad/s, ascending) of the assembled reduced
    /// system `K_r q = ω² M_r q`.
    ///
    /// # Errors
    /// [`Error::SingularSystem`] if the reduced mass matrix is not positive definite.
    pub fn reduced_frequencies(&self) -> Result<Vec<f64>> {
        let (lambdas, _) =
            generalized_symmetric_eigen(&self.reduced_stiffness_matrix, &self.reduced_mass_matrix)?;
        Ok(lambdas.iter().map(|&l| l.max(0.0).sqrt()).collect())
    }

    /// Modal participation factor of fixed-interface mode `mode` of substructure
    /// `sub` for an interface motion pattern `r` (length `n_interface`):
    ///
    /// Γₖ = φₖᵀ (M_ib + M_ii Ψ) r,
    ///
    /// i.e. the mode/interface coupling block of the reduced mass matrix applied to `r`.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `sub`/`mode` is out of range or `r` has the
    /// wrong length.
    pub fn modal_participation(&self, sub: usize, mode: usize, r: &[f64]) -> Result<f64> {
        let comp = self.components.get(sub).ok_or_else(|| {
            Error::InvalidInput(format!("substructure {sub} out of range ({})", self.n_sub))
        })?;
        if mode >= comp.n_modes() {
            return Err(Error::InvalidInput(format!(
                "mode {mode} out of range ({} retained)",
                comp.n_modes()
            )));
        }
        if r.len() != self.n_interface {
            return Err(Error::InvalidInput(format!(
                "influence vector has length {}, expected {}",
                r.len(),
                self.n_interface
            )));
        }
        let row = &comp.reduced_mass[self.n_interface + mode];
        Ok(row[..self.n_interface]
            .iter()
            .zip(r.iter())
            .map(|(a, b)| a * b)
            .sum())
    }

    /// Modal participation factor Γₖ = φₖᵀ (M_ib + M_ii Ψ) r for a uniform unit
    /// motion of all interface DOFs (`r = 1`), which is the rigid-body
    /// translation influence vector of a scalar-DOF (e.g. axial or spring-chain)
    /// model.  Returns `0.0` for an out-of-range substructure or mode.
    pub fn participation_factor(&self, sub: usize, mode: usize) -> f64 {
        let r = vec![1.0; self.n_interface];
        self.modal_participation(sub, mode, &r).unwrap_or(0.0)
    }

    /// Number of fixed-interface modes retained for substructure `sub`.
    pub fn mode_count(&self, sub: usize) -> usize {
        self.modal_freqs.get(sub).map_or(0, Vec::len)
    }

    /// Recover the full displacement of substructure `sub` from a reduced
    /// global vector `q` (length [`Self::total_reduced_dof`]).
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `sub` is out of range or `q` has the wrong length.
    pub fn expand(&self, sub: usize, q: &[f64]) -> Result<Vec<f64>> {
        let comp = self.components.get(sub).ok_or_else(|| {
            Error::InvalidInput(format!("substructure {sub} out of range ({})", self.n_sub))
        })?;
        if q.len() != self.total_reduced_dof() {
            return Err(Error::InvalidInput(format!(
                "reduced vector has length {}, expected {}",
                q.len(),
                self.total_reduced_dof()
            )));
        }
        let start = self.n_interface + sub * self.n_modes;
        comp.expand(&q[..self.n_interface], &q[start..start + self.n_modes])
    }
}

// ---------------------------------------------------------------------------
// Dense helpers
// ---------------------------------------------------------------------------

fn validate_square_symmetric(a: &[Vec<f64>], n: usize, name: &str) -> Result<()> {
    if a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err(Error::InvalidInput(format!(
            "{name} matrix must be {n}×{n}"
        )));
    }
    let scale = a
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |acc, v| acc.max(v.abs()));
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate().skip(i + 1) {
            if (v - a[j][i]).abs() > SYMMETRY_RTOL * scale.max(f64::MIN_POSITIVE) {
                return Err(Error::InvalidInput(format!(
                    "{name} matrix is not symmetric at ({i}, {j})"
                )));
            }
        }
    }
    Ok(())
}

fn sub_matrix(a: &[Vec<f64>], rows: &[usize], cols: &[usize]) -> Vec<Vec<f64>> {
    rows.iter()
        .map(|&r| cols.iter().map(|&c| a[r][c]).collect())
        .collect()
}

/// Cholesky factor (row-major flat `L`) with a relative pivot check.
fn cholesky_checked(chol: &CholeskyDense, a: &[Vec<f64>]) -> Result<Vec<f64>> {
    let n = chol.n;
    let flat: Vec<f64> = a.iter().flat_map(|row| row.iter().copied()).collect();
    let l = chol.factorize(&flat).map_err(|_| Error::SingularSystem)?;
    let max_diag = (0..n).fold(0.0_f64, |acc, i| acc.max(a[i][i].abs()));
    for i in 0..n {
        let pivot = l[i * n + i];
        if !pivot.is_finite() || pivot * pivot <= PIVOT_RTOL * max_diag {
            return Err(Error::SingularSystem);
        }
    }
    Ok(l)
}

/// `Tᵀ A T` for `A` (n × n) and `T` (n × r), symmetrised.
fn congruence(a: &[Vec<f64>], t: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let r = t.first().map_or(0, Vec::len);
    // AT = A · T  (n × r)
    let mut at = vec![vec![0.0; r]; n];
    for (i, at_row) in at.iter_mut().enumerate() {
        for (k, &a_ik) in a[i].iter().enumerate() {
            if a_ik == 0.0 {
                continue;
            }
            for (dst, &t_kj) in at_row.iter_mut().zip(t[k].iter()) {
                *dst += a_ik * t_kj;
            }
        }
    }
    let mut out = vec![vec![0.0; r]; r];
    for (k, t_row) in t.iter().enumerate() {
        for (p, &t_kp) in t_row.iter().enumerate() {
            if t_kp == 0.0 {
                continue;
            }
            for (dst, &at_kq) in out[p].iter_mut().zip(at[k].iter()) {
                *dst += t_kp * at_kq;
            }
        }
    }
    symmetrised(&out, 1.0)
}

/// `(A + Aᵀ) / (2·scale)` for a square matrix.
fn symmetrised(a: &[Vec<f64>], scale: f64) -> Vec<Vec<f64>> {
    a.iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, &v)| 0.5 * (v + a[j][i]) / scale)
                .collect()
        })
        .collect()
}

/// Solve the symmetric-definite generalized eigenproblem `K φ = λ M φ`.
///
/// With `M = L Lᵀ` the problem becomes `(L⁻¹ K L⁻ᵀ) y = λ y`, solved by cyclic
/// Jacobi; the modes `φ = L⁻ᵀ y` are mass-normalised (`φᵀ M φ = 1`).
/// Returns eigenvalues in ascending order and the matching modes.
pub(crate) fn generalized_symmetric_eigen(
    k: &[Vec<f64>],
    m: &[Vec<f64>],
) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = k.len();
    if m.len() != n {
        return Err(Error::InvalidInput(
            "stiffness and mass matrices differ in size".into(),
        ));
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let chol = CholeskyDense::new(n);
    let l = cholesky_checked(&chol, m)?;
    // Y = L⁻¹ K   (forward substitution column by column)
    let mut y = vec![vec![0.0; n]; n];
    for j in 0..n {
        let col: Vec<f64> = (0..n).map(|i| k[i][j]).collect();
        let sol = forward_substitution(&l, n, &col);
        for (i, v) in sol.into_iter().enumerate() {
            y[i][j] = v;
        }
    }
    // A = L⁻¹ Yᵀ = L⁻¹ K L⁻ᵀ   (K symmetric ⇒ Yᵀ = K L⁻ᵀ)
    let mut a = vec![vec![0.0; n]; n];
    for j in 0..n {
        let col: Vec<f64> = (0..n).map(|i| y[j][i]).collect();
        let sol = forward_substitution(&l, n, &col);
        for (i, v) in sol.into_iter().enumerate() {
            a[i][j] = v;
        }
    }
    let scale = (0..n).fold(0.0_f64, |acc, i| acc.max(a[i][i].abs()));
    let scale = if scale > 0.0 { scale } else { 1.0 };
    let a = symmetrised(&a, scale);
    let (mu, vecs) = jacobi_eigen_dense(&a, n);
    let lambdas: Vec<f64> = mu.iter().map(|&v| v * scale).collect();
    let modes: Vec<Vec<f64>> = vecs
        .iter()
        .map(|yv| backward_substitution_transposed(&l, n, yv))
        .collect();
    Ok((lambdas, modes))
}

/// Solve `L x = b` for row-major lower-triangular `L`.
fn forward_substitution(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * x[k];
        }
        x[i] = s / l[i * n + i];
    }
    x
}

/// Solve `Lᵀ x = b` for row-major lower-triangular `L`.
fn backward_substitution_transposed(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in (i + 1)..n {
            s -= l[k * n + i] * x[k];
        }
        x[i] = s / l[i * n + i];
    }
    x
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Fixed-free spring-mass chain with `n` unit masses and unit springs
    /// (ground spring on DOF 0, free end at DOF n−1).
    fn chain(n: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut k = vec![vec![0.0; n]; n];
        for i in 0..n {
            k[i][i] = if i + 1 < n { 2.0 } else { 1.0 };
            if i + 1 < n {
                k[i][i + 1] = -1.0;
                k[i + 1][i] = -1.0;
            }
        }
        let m: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut row = vec![0.0; n];
                row[i] = 1.0;
                row
            })
            .collect();
        (k, m)
    }

    /// Analytic natural frequencies of the fixed-free unit chain.
    fn chain_freq(n: usize, j: usize) -> f64 {
        let theta = (2 * j + 1) as f64 * std::f64::consts::PI / (2.0 * (2 * n + 1) as f64);
        2.0 * theta.sin()
    }

    #[test]
    fn generalized_eigen_two_by_two() {
        // K = diag(2, 8), M = diag(2, 2) ⇒ λ = 1, 4 with φᵀMφ = 1.
        let k = vec![vec![2.0, 0.0], vec![0.0, 8.0]];
        let m = vec![vec![2.0, 0.0], vec![0.0, 2.0]];
        let (l, v) = generalized_symmetric_eigen(&k, &m).expect("eigen");
        assert!((l[0] - 1.0).abs() < 1e-12 && (l[1] - 4.0).abs() < 1e-12);
        assert!((v[0][0].abs() - 1.0 / 2.0_f64.sqrt()).abs() < 1e-12);
        assert!(v[0][1].abs() < 1e-12);
    }

    #[test]
    fn chain_full_eigen_matches_analytic() {
        let (k, m) = chain(10);
        let (l, _) = generalized_symmetric_eigen(&k, &m).expect("eigen");
        for (j, &lam) in l.iter().enumerate() {
            assert!((lam.sqrt() - chain_freq(10, j)).abs() < 1e-10);
        }
    }

    #[test]
    fn craig_bampton_first_frequency_within_one_percent() {
        // 10-DOF fixed-free chain, tip DOF as interface, 3 fixed-interface modes.
        let (k, m) = chain(10);
        let cb = SubstructureMethod::from_matrices(&k, &m, &[9], 3).expect("reduction");
        assert_eq!(cb.total_reduced_dof(), 4);
        let freqs = cb.reduced_frequencies().expect("frequencies");
        let exact = chain_freq(10, 0);
        assert!(
            ((freqs[0] - exact) / exact).abs() < 0.01,
            "reduced ω₁ = {}, exact {exact}",
            freqs[0]
        );
        // Rayleigh-Ritz: reduced frequencies are upper bounds.
        assert!(freqs[0] >= exact - 1e-12);
    }

    #[test]
    fn craig_bampton_all_modes_is_exact() {
        let (k, m) = chain(10);
        let cb = SubstructureMethod::from_matrices(&k, &m, &[9], 9).expect("reduction");
        let freqs = cb.reduced_frequencies().expect("frequencies");
        assert_eq!(freqs.len(), 10);
        for (j, &w) in freqs.iter().enumerate() {
            assert!((w - chain_freq(10, j)).abs() < 1e-8, "mode {j}: {w}");
        }
    }

    #[test]
    fn craig_bampton_modes_mass_normalised_and_constraint_modes_static() {
        let (k, m) = chain(10);
        let cb = SubstructureMethod::from_matrices(&k, &m, &[9], 4).expect("reduction");
        let comp = &cb.components[0];
        // Reduced modal block: M = I, K = diag(ω²).
        for a in 0..4 {
            for b in 0..4 {
                let m_ab = cb.reduced_mass()[1 + a][1 + b];
                let k_ab = cb.reduced_stiffness()[1 + a][1 + b];
                let delta = if a == b { 1.0 } else { 0.0 };
                assert!((m_ab - delta).abs() < 1e-10);
                let w2 = comp.fixed_interface_freqs[a].powi(2);
                assert!((k_ab - delta * w2).abs() < 1e-10);
            }
            // Static constraint modes decouple stiffness: K_r[mode][boundary] = 0.
            assert!(cb.reduced_stiffness()[1 + a][0].abs() < 1e-10);
        }
        // Fixed-interface chain (fixed at both ends, 9 interior DOFs): Ψ is the
        // linear interpolation i/(9+1)·… of a unit tip displacement.
        for (row, &dof) in comp.constraint_modes.iter().zip(comp.interior_dofs.iter()) {
            let expected = (dof + 1) as f64 / 10.0;
            assert!((row[0] - expected).abs() < 1e-12);
        }
        // Guyan condensed stiffness of a 10-spring series chain = 1/10.
        assert!((cb.reduced_stiffness()[0][0] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn craig_bampton_two_substructures_assembled() {
        // Split the 10-DOF chain at DOF 4 (shared interface node).
        let (k_full, m_full) = chain(10);
        // Substructure A: DOFs 0..=4 with springs ground-0 … 3-4; node 4 gets half mass.
        let mut k_a = vec![vec![0.0; 5]; 5];
        let mut m_a = vec![vec![0.0; 5]; 5];
        for i in 0..5 {
            k_a[i][i] = if i < 4 { 2.0 } else { 1.0 };
            m_a[i][i] = if i < 4 { 1.0 } else { 0.5 };
            if i + 1 < 5 {
                k_a[i][i + 1] = -1.0;
                k_a[i + 1][i] = -1.0;
            }
        }
        // Substructure B: DOFs 4..=9 (local 0..=5); local 0 is the interface.
        let mut k_b = vec![vec![0.0; 6]; 6];
        let mut m_b = vec![vec![0.0; 6]; 6];
        for i in 0..6 {
            k_b[i][i] = if i == 0 || i == 5 { 1.0 } else { 2.0 };
            m_b[i][i] = if i == 0 { 0.5 } else { 1.0 };
            if i + 1 < 6 {
                k_b[i][i + 1] = -1.0;
                k_b[i + 1][i] = -1.0;
            }
        }
        let (exact, _) = generalized_symmetric_eigen(&k_full, &m_full).expect("eigen");
        let exact: Vec<f64> = exact.iter().map(|l| l.sqrt()).collect();
        // A has 4 interior DOFs (all kept), B has 5 (its highest mode is
        // truncated): every reduced frequency is a Rayleigh-Ritz upper bound
        // and the fundamental is captured almost exactly.
        let subs: [SubstructureInput<'_>; 2] = [(&k_a, &m_a, &[4]), (&k_b, &m_b, &[0])];
        let cb = SubstructureMethod::from_substructures(&subs, 4).expect("reduction");
        assert_eq!(cb.total_reduced_dof(), 9);
        let freqs = cb.reduced_frequencies().expect("frequencies");
        for (j, &w) in freqs.iter().enumerate() {
            assert!(w >= exact[j] - 1e-10, "mode {j} must be an upper bound");
        }
        assert!(((freqs[0] - exact[0]) / exact[0]).abs() < 1e-3);

        // Expanding the first reduced mode reproduces a continuous shape at the
        // shared node.
        let (_, modes) =
            generalized_symmetric_eigen(cb.reduced_stiffness(), cb.reduced_mass()).expect("eig");
        let u_a = cb.expand(0, &modes[0]).expect("expand A");
        let u_b = cb.expand(1, &modes[0]).expect("expand B");
        assert!((u_a[4] - u_b[0]).abs() < 1e-14);
    }

    #[test]
    fn craig_bampton_participation_factor_matches_mass_coupling() {
        let (k, m) = chain(10);
        let cb = SubstructureMethod::from_matrices(&k, &m, &[9], 3).expect("reduction");
        let comp = &cb.components[0];
        // Γ₀ = φ₀ᵀ (M_ib + M_ii Ψ) 1 = Σᵢ φ₀[i] Ψ[i]  (M = I, M_ib = 0)
        let expected: f64 = comp.normal_modes[0]
            .iter()
            .zip(comp.constraint_modes.iter())
            .map(|(p, psi)| p * psi[0])
            .sum();
        let gamma = cb.participation_factor(0, 0);
        assert!((gamma - expected).abs() < 1e-12);
        assert!(gamma.abs() > 1e-3);
        assert_eq!(cb.participation_factor(1, 0), 0.0);
        assert!(cb.modal_participation(0, 5, &[1.0]).is_err());
    }

    #[test]
    fn craig_bampton_rejects_bad_input() {
        let (k, m) = chain(4);
        assert!(matches!(
            SubstructureMethod::from_matrices(&k, &m, &[7], 1),
            Err(Error::DofOutOfRange { dof: 7, n: 4 })
        ));
        assert!(SubstructureMethod::from_matrices(&k, &m, &[3, 3], 1).is_err());
        assert!(SubstructureMethod::from_matrices(&k, &m, &[3], 4).is_err());
        // Free-free interior (no ground spring) with no boundary: singular K_ii.
        let k_free = vec![vec![1.0, -1.0], vec![-1.0, 1.0]];
        let m2 = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(matches!(
            SubstructureMethod::from_matrices(&k_free, &m2, &[], 1),
            Err(Error::SingularSystem)
        ));
    }

    #[test]
    fn craig_bampton_mode_count_and_dofs() {
        let (k, m) = chain(8);
        let cb = SubstructureMethod::from_matrices(&k, &m, &[0, 7], 3).expect("reduction");
        assert_eq!(cb.mode_count(0), 3);
        assert_eq!(cb.mode_count(1), 0);
        assert_eq!(cb.reduced_dof(), 5);
        assert_eq!(cb.total_reduced_dof(), 5);
        let freqs = &cb.modal_freqs[0];
        for w in freqs.windows(2) {
            assert!(w[1] > w[0]);
        }
    }
}
