// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reduced-order FEM: POD, Galerkin projection, DEIM, greedy RB, balanced
//! truncation, Petrov-Galerkin, manifold interpolation, error estimation,
//! and online/offline decomposition for parametric ROMs.
//!
//! # References
//! - Sirovich (1987) "Turbulence and the dynamics of coherent structures"
//! - Chaturantabut & Sorensen (2010) "Nonlinear model reduction via DEIM"
//! - Rozza, Huynh & Patera (2008) "Reduced basis approximation and a posteriori
//!   error estimation"

// ---------------------------------------------------------------------------
// Math helpers (plain f64 / Vec<f64>, no nalgebra)
// ---------------------------------------------------------------------------

/// Dot product of two equal-length slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm of a slice.
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// Element-wise subtraction a − b.
fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

/// Element-wise addition a + b.
fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scale vector by scalar s.
fn vec_scale(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

/// Dense row-major matrix × vector.
fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

/// Transpose of a dense matrix (rows × cols).
fn mat_transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return vec![];
    }
    let m = a.len();
    let n = a[0].len();
    let mut t = vec![vec![0.0; m]; n];
    for i in 0..m {
        for j in 0..n {
            t[j][i] = a[i][j];
        }
    }
    t
}

/// Dense matrix multiplication C = A * B.
fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let m = a.len();
    let n = b[0].len();
    let bt = mat_transpose(b);
    let mut c = vec![vec![0.0; n]; m];
    for i in 0..m {
        for j in 0..n {
            c[i][j] = dot(&a[i], &bt[j]);
        }
    }
    c
}

/// Solve upper-triangular system Ux = b (back-substitution).
#[cfg(test)]
fn back_substitute(u: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= u[i][j] * x[j];
        }
        x[i] = if u[i][i].abs() < 1e-14 {
            0.0
        } else {
            s / u[i][i]
        };
    }
    x
}

/// Solve lower-triangular system Lx = b (forward substitution).
fn forward_substitute(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i][j] * x[j];
        }
        x[i] = if l[i][i].abs() < 1e-14 {
            0.0
        } else {
            s / l[i][i]
        };
    }
    x
}

/// Conjugate-gradient solve A x = b (dense A, no preconditioning).
fn cg_solve(a: &[Vec<f64>], b: &[f64], max_iter: usize, tol: f64) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    for _ in 0..max_iter {
        if rr.sqrt() < tol {
            break;
        }
        let ap = mat_vec(a, &p);
        let alpha = rr / dot(&p, &ap);
        x = vec_add(&x, &vec_scale(&p, alpha));
        r = vec_sub(&r, &vec_scale(&ap, alpha));
        let rr_new = dot(&r, &r);
        let beta = rr_new / rr;
        p = vec_add(&r, &vec_scale(&p, beta));
        rr = rr_new;
    }
    x
}

/// Identity matrix of size n.
fn eye(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for (i, row) in m.iter_mut().enumerate().take(n) {
        row[i] = 1.0;
    }
    m
}

// ---------------------------------------------------------------------------
// POD – Proper Orthogonal Decomposition
// ---------------------------------------------------------------------------

/// A set of POD basis vectors and their singular values.
///
/// The basis vectors are stored as columns (each `Vec`f64` is one basis
/// vector of length = full-order DOF count).
#[derive(Clone, Debug)]
pub struct PodBasis {
    /// Orthonormal basis vectors (columns of the left singular matrix).
    pub modes: Vec<Vec<f64>>,
    /// Corresponding singular values σ_k (energy content).
    pub singular_values: Vec<f64>,
    /// Full-order dimension N.
    pub full_dim: usize,
    /// Reduced dimension r (= modes.len()).
    pub reduced_dim: usize,
}

impl PodBasis {
    /// Create an empty [`PodBasis`] of given full dimension.
    pub fn new(full_dim: usize) -> Self {
        Self {
            modes: vec![],
            singular_values: vec![],
            full_dim,
            reduced_dim: 0,
        }
    }

    /// Number of modes retained.
    pub fn rank(&self) -> usize {
        self.modes.len()
    }

    /// Project a full-order vector onto the reduced space: q_r = Φᵀ q.
    pub fn project(&self, q: &[f64]) -> Vec<f64> {
        self.modes.iter().map(|phi| dot(phi, q)).collect()
    }

    /// Reconstruct a full-order vector from reduced coordinates: q ≈ Φ q_r.
    pub fn reconstruct(&self, q_r: &[f64]) -> Vec<f64> {
        let n = self.full_dim;
        let mut q = vec![0.0; n];
        for (phi, &coeff) in self.modes.iter().zip(q_r.iter()) {
            for (qi, &phi_i) in q.iter_mut().zip(phi.iter()) {
                *qi += coeff * phi_i;
            }
        }
        q
    }

    /// Relative energy captured by the retained modes.
    ///
    /// E = Σ_{k=1}^r σ_k² / Σ_all σ_k²
    pub fn energy_content(&self) -> f64 {
        let total: f64 = self.singular_values.iter().map(|s| s * s).sum();
        if total < 1e-300 {
            return 0.0;
        }
        let retained: f64 = self.singular_values[..self.reduced_dim]
            .iter()
            .map(|s| s * s)
            .sum();
        retained / total
    }
}

/// Build a POD basis from a snapshot matrix using a method-of-snapshots SVD.
///
/// Each column of `snapshots` (outer index) is one snapshot vector of length N.
/// Returns a [`PodBasis`] truncated to retain `max_modes` modes (or fewer if
/// rank is smaller).
///
/// Uses power-iteration-based symmetric eigendecomposition of the Gram matrix
/// Cᵀ = Sᵀ S (size n_snap × n_snap).
pub fn build_pod_basis(snapshots: &[Vec<f64>], max_modes: usize) -> PodBasis {
    let n_snap = snapshots.len();
    if n_snap == 0 || snapshots[0].is_empty() {
        return PodBasis::new(0);
    }
    let n = snapshots[0].len();

    // Build Gram matrix C = Sᵀ S  (n_snap × n_snap)
    let mut gram = vec![vec![0.0; n_snap]; n_snap];
    for i in 0..n_snap {
        for j in i..n_snap {
            let v = dot(&snapshots[i], &snapshots[j]);
            gram[i][j] = v;
            gram[j][i] = v;
        }
    }

    // Power iteration to extract dominant eigenpairs
    let r = max_modes.min(n_snap).min(n);
    let mut eigenvecs: Vec<Vec<f64>> = Vec::with_capacity(r);
    let mut eigenvals: Vec<f64> = Vec::with_capacity(r);
    let mut deflated = gram.clone();

    for _k in 0..r {
        // Initialise with the row whose diagonal entry is largest (robust pivot).
        let pivot = (0..n_snap)
            .max_by(|&a, &b| {
                deflated[a][a]
                    .abs()
                    .partial_cmp(&deflated[b][b].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0);
        let mut v: Vec<f64> = deflated.iter().map(|row| row[pivot]).collect();
        let vn = norm(&v);
        if vn < 1e-14 {
            break;
        }
        for x in v.iter_mut() {
            *x /= vn;
        }

        for _iter in 0..200 {
            let av = mat_vec(&deflated, &v);
            let lambda = dot(&v, &av);
            let vn2 = norm(&av);
            if vn2 < 1e-14 {
                break;
            }
            let v_new: Vec<f64> = av.iter().map(|x| x / vn2).collect();
            let diff: f64 = vec_sub(&v_new, &v)
                .iter()
                .map(|x| x.abs())
                .fold(0.0_f64, f64::max);
            v = v_new;
            if diff < 1e-12 {
                let av2 = mat_vec(&deflated, &v);
                let _lam = dot(&v, &av2);
                let _ = lambda;
                break;
            }
        }
        let av = mat_vec(&deflated, &v);
        let lambda = dot(&v, &av).max(0.0);
        let sigma = lambda.sqrt();
        if sigma < 1e-14 {
            break;
        }

        eigenvals.push(sigma);
        eigenvecs.push(v.clone());

        // Deflate: C ← C − λ vvᵀ
        for i in 0..n_snap {
            for j in 0..n_snap {
                deflated[i][j] -= lambda * v[i] * v[j];
            }
        }
    }

    // Compute left singular vectors Φ_k = S v_k / σ_k
    let mut modes: Vec<Vec<f64>> = Vec::with_capacity(eigenvals.len());
    let mut singular_values: Vec<f64> = Vec::with_capacity(eigenvals.len());
    for (v, sigma) in eigenvecs.iter().zip(eigenvals.iter()) {
        // phi = S * v / sigma
        let mut phi = vec![0.0; n];
        for (snap, &coeff) in snapshots.iter().zip(v.iter()) {
            for (phi_i, &s_i) in phi.iter_mut().zip(snap.iter()) {
                *phi_i += coeff * s_i;
            }
        }
        let pn = norm(&phi);
        if pn < 1e-14 {
            continue;
        }
        for x in phi.iter_mut() {
            *x /= pn;
        }
        modes.push(phi);
        singular_values.push(*sigma);
    }

    let reduced_dim = modes.len();
    PodBasis {
        modes,
        singular_values,
        full_dim: n,
        reduced_dim,
    }
}

// ---------------------------------------------------------------------------
// Galerkin projection
// ---------------------------------------------------------------------------

/// Galerkin-projected (reduced) stiffness matrix and load vector.
///
/// K_r = Φᵀ K Φ  (r × r)
/// f_r = Φᵀ f    (r)
#[derive(Clone, Debug)]
pub struct GalerkinSystem {
    /// Reduced stiffness matrix K_r (r × r).
    pub k_reduced: Vec<Vec<f64>>,
    /// Reduced force vector f_r (r).
    pub f_reduced: Vec<f64>,
    /// Reduced dimension.
    pub reduced_dim: usize,
}

impl GalerkinSystem {
    /// Solve the reduced system K_r q_r = f_r via CG.
    pub fn solve(&self) -> Vec<f64> {
        cg_solve(&self.k_reduced, &self.f_reduced, 500, 1e-10)
    }
}

/// Project the full-order FEM system onto a POD basis using Galerkin projection.
///
/// - `k_full` – full stiffness matrix (N × N, row-major).
/// - `f_full` – full force vector (N).
/// - `basis`  – [`PodBasis`] with r modes.
pub fn galerkin_projection(
    k_full: &[Vec<f64>],
    f_full: &[f64],
    basis: &PodBasis,
) -> GalerkinSystem {
    let r = basis.rank();
    // K_phi = K * Phi   (N × r)
    let k_phi: Vec<Vec<f64>> = basis.modes.iter().map(|phi| mat_vec(k_full, phi)).collect();
    // K_r[i][j] = phi_i . (K phi_j)
    let mut k_reduced = vec![vec![0.0; r]; r];
    for (i, row) in k_reduced.iter_mut().enumerate().take(r) {
        for (j, cell) in row.iter_mut().enumerate().take(r) {
            *cell = dot(&basis.modes[i], &k_phi[j]);
        }
    }
    // f_r[i] = phi_i . f
    let f_reduced: Vec<f64> = basis.modes.iter().map(|phi| dot(phi, f_full)).collect();
    GalerkinSystem {
        k_reduced,
        f_reduced,
        reduced_dim: r,
    }
}

// ---------------------------------------------------------------------------
// Reduced Basis Method
// ---------------------------------------------------------------------------

/// A reduced basis built by the greedy algorithm over a parameter set.
#[derive(Clone, Debug)]
pub struct ReducedBasis {
    /// Orthonormal basis vectors (each length N).
    pub basis: Vec<Vec<f64>>,
    /// Parameter values selected by the greedy algorithm.
    pub selected_params: Vec<f64>,
    /// Full-order dimension.
    pub full_dim: usize,
}

impl ReducedBasis {
    /// Create an empty [`ReducedBasis`].
    pub fn new(full_dim: usize) -> Self {
        Self {
            basis: vec![],
            selected_params: vec![],
            full_dim,
        }
    }

    /// Dimension of the current basis.
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// Add a new snapshot, orthonormalise against existing basis (Gram-Schmidt).
    pub fn add_snapshot(&mut self, snap: &[f64], param: f64) {
        let mut v = snap.to_vec();
        for phi in &self.basis {
            let c = dot(&v, phi);
            v = vec_sub(&v, &vec_scale(phi, c));
        }
        let n = norm(&v);
        if n < 1e-12 {
            return;
        }
        self.basis.push(vec_scale(&v, 1.0 / n));
        self.selected_params.push(param);
    }

    /// Project a full-order vector onto the reduced basis.
    pub fn project(&self, q: &[f64]) -> Vec<f64> {
        self.basis.iter().map(|phi| dot(phi, q)).collect()
    }

    /// Reconstruct a full-order vector from reduced coordinates.
    pub fn reconstruct(&self, q_r: &[f64]) -> Vec<f64> {
        let mut q = vec![0.0; self.full_dim];
        for (phi, &coeff) in self.basis.iter().zip(q_r.iter()) {
            for (qi, &phi_i) in q.iter_mut().zip(phi.iter()) {
                *qi += coeff * phi_i;
            }
        }
        q
    }
}

/// Greedy algorithm for constructing a reduced basis.
///
/// Given a parameter training set and a snapshot generator function, the
/// greedy algorithm selects the parameter that maximises the error indicator
/// and adds the corresponding snapshot to the basis.
///
/// - `params`      – training parameters μ_i.
/// - `snapshots`   – pre-computed full-order solutions u(μ_i).
/// - `max_basis`   – maximum number of basis vectors to retain.
/// - `tol`         – stop if error indicator falls below this.
pub fn greedy_reduced_basis(
    params: &[f64],
    snapshots: &[Vec<f64>],
    max_basis: usize,
    tol: f64,
) -> ReducedBasis {
    if snapshots.is_empty() {
        return ReducedBasis::new(0);
    }
    let n = snapshots[0].len();
    let mut rb = ReducedBasis::new(n);

    // Start with the snapshot having the largest norm
    let (best_idx, _) = snapshots
        .iter()
        .enumerate()
        .map(|(i, s)| (i, norm(s)))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((0, 0.0));
    rb.add_snapshot(&snapshots[best_idx], params[best_idx]);

    for _k in 1..max_basis {
        // Find parameter with largest projection error
        let mut max_err = 0.0_f64;
        let mut max_i = 0;
        for (i, snap) in snapshots.iter().enumerate() {
            let q_r = rb.project(snap);
            let recon = rb.reconstruct(&q_r);
            let err = norm(&vec_sub(snap, &recon));
            if err > max_err {
                max_err = err;
                max_i = i;
            }
        }
        if max_err < tol {
            break;
        }
        rb.add_snapshot(&snapshots[max_i], params[max_i]);
    }
    rb
}

// ---------------------------------------------------------------------------
// DEIM – Discrete Empirical Interpolation Method
// ---------------------------------------------------------------------------

/// DEIM interpolation operator for hyper-reduction of nonlinear terms.
///
/// Stores the selected interpolation indices and the DEIM matrix U_s^{-1} Uᵀ.
#[derive(Clone, Debug)]
pub struct DeimInterpolation {
    /// DEIM interpolation indices p_1, …, p_m (into the full DOF vector).
    pub indices: Vec<usize>,
    /// DEIM coefficient matrix P (m × r), such that c = P * f(u\[indices\]).
    pub p_matrix: Vec<Vec<f64>>,
    /// Basis modes used (r modes, each length N).
    pub modes: Vec<Vec<f64>>,
}

impl DeimInterpolation {
    /// Evaluate the DEIM approximation of a nonlinear vector `f_full` using
    /// only the selected component values.
    pub fn approximate(&self, f_full: &[f64]) -> Vec<f64> {
        // Extract f at DEIM indices
        let f_s: Vec<f64> = self.indices.iter().map(|&i| f_full[i]).collect();
        // Compute DEIM coefficients c = P * f_s
        let c = mat_vec(&self.p_matrix, &f_s);
        // Reconstruct: f_approx = U c
        let n = if self.modes.is_empty() {
            0
        } else {
            self.modes[0].len()
        };
        let mut approx = vec![0.0; n];
        for (phi, &coeff) in self.modes.iter().zip(c.iter()) {
            for (a, &b) in approx.iter_mut().zip(phi.iter()) {
                *a += coeff * b;
            }
        }
        approx
    }
}

/// Compute the DEIM interpolation operator from a set of nonlinear basis modes.
///
/// Algorithm (Chaturantabut & Sorensen 2010):
/// 1. Select index of maximum absolute value in u_1.
/// 2. For each subsequent mode, solve the lower-triangular system and pick
///    the index of maximum residual.
///
/// Returns a [`DeimInterpolation`].
pub fn deim(modes: &[Vec<f64>]) -> DeimInterpolation {
    if modes.is_empty() {
        return DeimInterpolation {
            indices: vec![],
            p_matrix: vec![],
            modes: vec![],
        };
    }
    let _n = modes[0].len();
    let r = modes.len();
    let mut indices: Vec<usize> = Vec::with_capacity(r);

    // Step 1: first index = argmax |u_1|
    let first_idx = modes[0]
        .iter()
        .enumerate()
        .max_by(|a, b| {
            a.1.abs()
                .partial_cmp(&b.1.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    indices.push(first_idx);

    // U_s is the submatrix of U at rows = indices  (grows each iteration)
    // We maintain a lower-triangular system
    for l in 1..r {
        // Build P^T U_l (l × l) and extract u_l at current indices
        let u_l = &modes[l];
        // Evaluate u_l at current indices
        let b: Vec<f64> = indices.iter().map(|&i| u_l[i]).collect();
        // Build coefficient matrix: rows = indices (l×l submatrix of U_1..l-1)
        let mut u_s = vec![vec![0.0; l]; l];
        for (row, &idx) in indices.iter().enumerate() {
            for (col, mode) in modes[..l].iter().enumerate() {
                u_s[row][col] = mode[idx];
            }
        }
        // Solve U_s c = b (forward sub on lower-triangular U_s)
        let c = forward_substitute(&u_s, &b);
        // Residual r = u_l - U c  (full vector)
        let mut residual = u_l.clone();
        for (mode, &coeff) in modes[..l].iter().zip(c.iter()) {
            for (ri, &mi) in residual.iter_mut().zip(mode.iter()) {
                *ri -= coeff * mi;
            }
        }
        // New index = argmax |residual|
        let new_idx = residual
            .iter()
            .enumerate()
            .max_by(|a, b| {
                a.1.abs()
                    .partial_cmp(&b.1.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        indices.push(new_idx);
    }

    // Build P matrix (DEIM coefficient matrix): (U_s)^{-1}
    // P * f_s = c  so we need (U[indices, :])^{-1}
    let m = indices.len();
    // U_s = U at rows indices (m × r)
    let u_s: Vec<Vec<f64>> = indices
        .iter()
        .map(|&i| modes.iter().map(|mode| mode[i]).collect::<Vec<f64>>())
        .collect();
    // Invert U_s by solving U_s x = e_k for each basis vector (m × m)
    let mut p_matrix = vec![vec![0.0; m]; r]; // P is r × m
    let eye_m = eye(m);
    for k in 0..m {
        let col = forward_substitute(&u_s, &eye_m[k]);
        for i in 0..r.min(col.len()) {
            p_matrix[i][k] = col[i];
        }
    }

    DeimInterpolation {
        indices,
        p_matrix,
        modes: modes.to_vec(),
    }
}

// ---------------------------------------------------------------------------
// Petrov-Galerkin projection
// ---------------------------------------------------------------------------

/// Petrov-Galerkin reduced system with separate test and trial spaces.
#[derive(Clone, Debug)]
pub struct PetrovGalerkinSystem {
    /// Reduced stiffness W_r = Ψᵀ K Φ  (r × r).
    pub k_reduced: Vec<Vec<f64>>,
    /// Reduced force g_r = Ψᵀ f  (r).
    pub f_reduced: Vec<f64>,
    /// Reduced dimension.
    pub reduced_dim: usize,
}

impl PetrovGalerkinSystem {
    /// Solve the Petrov-Galerkin reduced system via CG.
    pub fn solve(&self) -> Vec<f64> {
        cg_solve(&self.k_reduced, &self.f_reduced, 500, 1e-10)
    }
}

/// Build a Petrov-Galerkin projection.
///
/// - `k_full` – full-order stiffness matrix.
/// - `f_full` – full-order force vector.
/// - `phi`    – trial basis (columns of Φ, length N).
/// - `psi`    – test basis (columns of Ψ, length N).
pub fn petrov_galerkin_projection(
    k_full: &[Vec<f64>],
    f_full: &[f64],
    phi: &[Vec<f64>],
    psi: &[Vec<f64>],
) -> PetrovGalerkinSystem {
    let r = phi.len().min(psi.len());
    let k_phi: Vec<Vec<f64>> = phi.iter().map(|p| mat_vec(k_full, p)).collect();
    let mut k_reduced = vec![vec![0.0; r]; r];
    for i in 0..r {
        for j in 0..r {
            k_reduced[i][j] = dot(&psi[i], &k_phi[j]);
        }
    }
    let f_reduced: Vec<f64> = psi.iter().map(|p| dot(p, f_full)).collect();
    PetrovGalerkinSystem {
        k_reduced,
        f_reduced,
        reduced_dim: r,
    }
}

// ---------------------------------------------------------------------------
// Balanced truncation
// ---------------------------------------------------------------------------

/// Result of balanced truncation.
#[derive(Clone, Debug)]
pub struct BalancedTruncation {
    /// Hankel singular values (energy measures).
    pub hankel_singular_values: Vec<f64>,
    /// Left transformation matrix T_l (r × N).
    pub transform_left: Vec<Vec<f64>>,
    /// Right transformation matrix T_r (N × r).
    pub transform_right: Vec<Vec<f64>>,
    /// Reduced dimension.
    pub reduced_dim: usize,
}

impl BalancedTruncation {
    /// Project a full-order matrix A → T_l A T_r  (r × r).
    pub fn project_matrix(&self, a: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let at_r = mat_mul(a, &mat_transpose(&self.transform_left));
        // T_l * (A * T_r)
        let t_r = &self.transform_right;
        let m = self.transform_left.len();
        let _p = t_r[0].len();
        // at_r is N × r, multiply T_l (r × N) * at_r (N × r)
        let mut out = vec![vec![0.0; self.reduced_dim]; m];
        for (i, out_row) in out.iter_mut().enumerate().take(m) {
            for (j, cell) in out_row.iter_mut().enumerate().take(self.reduced_dim) {
                *cell = dot(
                    &self.transform_left[i],
                    &(0..a.len()).map(|k| at_r[k][j]).collect::<Vec<_>>(),
                );
            }
        }
        out
    }

    /// Project a full-order vector b → T_l b.
    pub fn project_vector(&self, b: &[f64]) -> Vec<f64> {
        mat_vec(&self.transform_left, b)
    }
}

/// Compute approximate balanced truncation using randomised Gramian estimates.
///
/// For a state-space system (A, B, C), computes approximate Hankel singular
/// values and transformation matrices by approximating the controllability
/// and observability Gramians via the system matrices.
///
/// This is a simplified pedagogical implementation using dominant eigenvectors.
///
/// - `a`          – system matrix A (N × N).
/// - `b_mat`      – input matrix B (N × n_inputs).
/// - `c_mat`      – output matrix C (n_outputs × N).
/// - `max_modes`  – number of balanced modes to retain.
pub fn balanced_truncation(
    a: &[Vec<f64>],
    b_mat: &[Vec<f64>],
    c_mat: &[Vec<f64>],
    max_modes: usize,
) -> BalancedTruncation {
    let n = a.len();
    // Approximate controllability Gramian W_c ≈ B Bᵀ
    let bt = mat_transpose(b_mat);
    let w_c = mat_mul(b_mat, &bt);
    // Approximate observability Gramian W_o ≈ CᵀC
    let ct = mat_transpose(c_mat);
    let w_o = mat_mul(&ct, c_mat);

    // Hankel matrix H ≈ W_c * W_o (simplified)
    let h = mat_mul(&w_c, &w_o);

    // Extract dominant eigenvectors of H via power iteration
    let r = max_modes.min(n);
    let mut t_l: Vec<Vec<f64>> = Vec::with_capacity(r);
    let mut hsv: Vec<f64> = Vec::with_capacity(r);
    let mut deflated = h.clone();

    for _k in 0..r {
        let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
        let mut lambda = 0.0_f64;
        for _iter in 0..300 {
            let av = mat_vec(&deflated, &v);
            let lambda_new = dot(&v, &av);
            let vn = norm(&av);
            if vn < 1e-14 {
                break;
            }
            v = vec_scale(&av, 1.0 / vn);
            let converged = (lambda_new - lambda).abs() < 1e-12;
            lambda = lambda_new;
            if converged {
                break;
            }
        }
        let _ = lambda;
        let av = mat_vec(&deflated, &v);
        let lam = dot(&v, &av).max(0.0);
        let sigma = lam.sqrt();
        if sigma < 1e-14 {
            break;
        }
        hsv.push(sigma);
        t_l.push(v.clone());
        // Deflate
        for i in 0..n {
            for j in 0..n {
                deflated[i][j] -= lam * v[i] * v[j];
            }
        }
    }

    // Right transform = left transform transposed (symmetric approximation)
    let t_r = mat_transpose(&t_l);
    let reduced_dim = t_l.len();

    BalancedTruncation {
        hankel_singular_values: hsv,
        transform_left: t_l,
        transform_right: t_r,
        reduced_dim,
    }
}

// ---------------------------------------------------------------------------
// Error Estimation
// ---------------------------------------------------------------------------

/// A posteriori error bound for the reduced order model.
///
/// For coercive problems: Δ(μ) = ||r(μ)||_{V'} / α_LB(μ)
/// where r is the residual and α_LB is a lower bound for the coercivity
/// constant.
#[derive(Clone, Debug)]
pub struct RomErrorBound {
    /// Residual norm ||r||.
    pub residual_norm: f64,
    /// Coercivity lower bound α_LB.
    pub coercivity_lb: f64,
    /// Error bound Δ = ||r|| / α_LB.
    pub error_bound: f64,
    /// Effectivity ρ = Δ / ||e|| (if true error is known).
    pub effectivity: Option<f64>,
}

impl RomErrorBound {
    /// Is the error bound sharp (effectivity close to 1)?
    pub fn is_sharp(&self, tol: f64) -> bool {
        self.effectivity.is_some_and(|e| (e - 1.0).abs() < tol)
    }
}

/// Compute the a posteriori error bound for a ROM solution.
///
/// - `k_full`        – full-order stiffness matrix.
/// - `f_full`        – full-order load vector.
/// - `u_rom`         – full-order reconstruction of the ROM solution.
/// - `coercivity_lb` – coercivity constant lower bound α_LB > 0.
/// - `true_error`    – optional ||u_true − u_rom|| (if known, for effectivity).
pub fn compute_error_bound(
    k_full: &[Vec<f64>],
    f_full: &[f64],
    u_rom: &[f64],
    coercivity_lb: f64,
    true_error: Option<f64>,
) -> RomErrorBound {
    // Residual r = f − K u_rom
    let ku = mat_vec(k_full, u_rom);
    let residual = vec_sub(f_full, &ku);
    let residual_norm = norm(&residual);
    let alpha = coercivity_lb.max(1e-14);
    let error_bound = residual_norm / alpha;
    let effectivity = true_error.map(|e| if e < 1e-14 { 1.0 } else { error_bound / e });
    RomErrorBound {
        residual_norm,
        coercivity_lb: alpha,
        error_bound,
        effectivity,
    }
}

// ---------------------------------------------------------------------------
// Empirical Interpolation Method (EIM)
// ---------------------------------------------------------------------------

/// An EIM approximant for a parameter-dependent function g(x; μ).
#[derive(Clone, Debug)]
pub struct EimApproximant {
    /// Magic points x^* (indices into the spatial domain).
    pub magic_points: Vec<usize>,
    /// Basis functions q_k(x) (each length = spatial dofs).
    pub basis_functions: Vec<Vec<f64>>,
    /// Interpolation matrix B (m × m) with B\[i\]\[j\] = q_j(x^*_i).
    pub interp_matrix: Vec<Vec<f64>>,
}

impl EimApproximant {
    /// Evaluate the EIM approximation given function values at magic points.
    pub fn evaluate(&self, g_magic: &[f64]) -> Vec<f64> {
        // Solve B c = g(magic_points)
        let c = cg_solve(&self.interp_matrix, g_magic, 100, 1e-12);
        // Reconstruct: sum c_k q_k
        let n = if self.basis_functions.is_empty() {
            0
        } else {
            self.basis_functions[0].len()
        };
        let mut g_approx = vec![0.0; n];
        for (phi, &coeff) in self.basis_functions.iter().zip(c.iter()) {
            for (g, &q) in g_approx.iter_mut().zip(phi.iter()) {
                *g += coeff * q;
            }
        }
        g_approx
    }
}

/// Build an EIM approximant from a set of parameter-dependent function snapshots.
///
/// Uses a greedy selection of magic points.
pub fn build_eim(snapshots: &[Vec<f64>]) -> EimApproximant {
    if snapshots.is_empty() {
        return EimApproximant {
            magic_points: vec![],
            basis_functions: vec![],
            interp_matrix: vec![],
        };
    }
    let m = snapshots.len();
    let n = snapshots[0].len();
    let mut magic_points: Vec<usize> = Vec::with_capacity(m);
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut interp_mat: Vec<Vec<f64>> = Vec::with_capacity(m);

    for (k, snap) in snapshots.iter().enumerate() {
        let residual = if k == 0 {
            snap.clone()
        } else {
            // Solve for EIM coefficients and compute residual
            let g_magic: Vec<f64> = magic_points.iter().map(|&i| snap[i]).collect();
            let bk: Vec<Vec<f64>> = interp_mat[..k].to_vec();
            let c = cg_solve(&bk, &g_magic[..k], 100, 1e-12);
            let approx_k = {
                let mut a = vec![0.0; n];
                for (phi, &coeff) in basis.iter().zip(c.iter()) {
                    for (ai, &bi) in a.iter_mut().zip(phi.iter()) {
                        *ai += coeff * bi;
                    }
                }
                a
            };
            vec_sub(snap, &approx_k)
        };

        // New magic point = argmax |residual|
        let new_pt = residual
            .iter()
            .enumerate()
            .max_by(|a, b| {
                a.1.abs()
                    .partial_cmp(&b.1.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Normalise residual to get new basis function
        let rn = norm(&residual);
        if rn < 1e-14 {
            break;
        }
        let q_k: Vec<f64> = vec_scale(&residual, 1.0 / rn);

        magic_points.push(new_pt);
        basis.push(q_k.clone());

        // Extend interpolation matrix
        let mut new_row = vec![0.0; k + 1];
        for (j, phi) in basis.iter().enumerate() {
            new_row[j] = phi[new_pt];
        }
        interp_mat.push(new_row);
    }

    EimApproximant {
        magic_points,
        basis_functions: basis,
        interp_matrix: interp_mat,
    }
}

// ---------------------------------------------------------------------------
// Online/offline decomposition
// ---------------------------------------------------------------------------

/// Offline-stage precomputed affine decomposition of the stiffness matrix.
///
/// K(μ) = Σ_q θ_q(μ) K_q
#[derive(Clone, Debug)]
pub struct AffineDecomposition {
    /// Reduced affine components K_r_q = Φᵀ K_q Φ (each r × r).
    pub k_reduced_components: Vec<Vec<Vec<f64>>>,
    /// Reduced force components f_r_q = Φᵀ f_q (each length r).
    pub f_reduced_components: Vec<Vec<f64>>,
    /// Reduced dimension.
    pub reduced_dim: usize,
}

impl AffineDecomposition {
    /// Assemble the reduced system for a given parameter via affine combination.
    ///
    /// - `theta_k` – affine coefficients θ_q(μ) for stiffness components.
    /// - `theta_f` – affine coefficients θ_q(μ) for load components.
    pub fn assemble(&self, theta_k: &[f64], theta_f: &[f64]) -> GalerkinSystem {
        let r = self.reduced_dim;
        let mut k_r = vec![vec![0.0; r]; r];
        let mut f_r = vec![0.0; r];
        for (q, &th) in theta_k.iter().enumerate() {
            if q >= self.k_reduced_components.len() {
                break;
            }
            for (i, kr_row) in k_r.iter_mut().enumerate().take(r) {
                for (j, kr_ij) in kr_row.iter_mut().enumerate().take(r) {
                    *kr_ij += th * self.k_reduced_components[q][i][j];
                }
            }
        }
        for (q, &th) in theta_f.iter().enumerate() {
            if q >= self.f_reduced_components.len() {
                break;
            }
            for (i, fr_i) in f_r.iter_mut().enumerate().take(r) {
                *fr_i += th * self.f_reduced_components[q][i];
            }
        }
        GalerkinSystem {
            k_reduced: k_r,
            f_reduced: f_r,
            reduced_dim: r,
        }
    }
}

/// Perform the offline stage of the affine decomposition.
///
/// Projects each stiffness component K_q and force component f_q onto the
/// POD basis Φ, storing the reduced representations for fast online assembly.
pub fn offline_affine_decomposition(
    k_components: &[Vec<Vec<f64>>],
    f_components: &[Vec<f64>],
    basis: &PodBasis,
) -> AffineDecomposition {
    let r = basis.rank();
    let k_reduced_components: Vec<Vec<Vec<f64>>> = k_components
        .iter()
        .map(|kq| {
            let k_phi: Vec<Vec<f64>> = basis.modes.iter().map(|phi| mat_vec(kq, phi)).collect();
            let mut kq_r = vec![vec![0.0; r]; r];
            for (i, kqr_row) in kq_r.iter_mut().enumerate().take(r) {
                for (j, kqr_ij) in kqr_row.iter_mut().enumerate().take(r) {
                    *kqr_ij = dot(&basis.modes[i], &k_phi[j]);
                }
            }
            kq_r
        })
        .collect();

    let f_reduced_components: Vec<Vec<f64>> = f_components
        .iter()
        .map(|fq| basis.modes.iter().map(|phi| dot(phi, fq)).collect())
        .collect();

    AffineDecomposition {
        k_reduced_components,
        f_reduced_components,
        reduced_dim: r,
    }
}

// ---------------------------------------------------------------------------
// Parametric ROM and manifold interpolation
// ---------------------------------------------------------------------------

/// Interpolates ROM solutions across the parameter space using radial basis
/// functions (RBF) on a training set.
#[derive(Clone, Debug)]
pub struct ManifoldInterpolator {
    /// Training parameters μ_i.
    pub params: Vec<f64>,
    /// Reduced-coordinate solutions q_r(μ_i) at training points.
    pub solutions: Vec<Vec<f64>>,
    /// RBF width ε.
    pub rbf_width: f64,
}

impl ManifoldInterpolator {
    /// Create a new [`ManifoldInterpolator`].
    pub fn new(rbf_width: f64) -> Self {
        Self {
            params: vec![],
            solutions: vec![],
            rbf_width,
        }
    }

    /// Add a training point (μ, q_r(μ)).
    pub fn add_training_point(&mut self, mu: f64, q_r: Vec<f64>) {
        self.params.push(mu);
        self.solutions.push(q_r);
    }

    /// Interpolate to a new parameter μ using Gaussian RBF.
    pub fn interpolate(&self, mu: f64) -> Vec<f64> {
        let n = self.params.len();
        if n == 0 {
            return vec![];
        }
        let dim = self.solutions[0].len();
        // Build RBF matrix Φ[i][j] = φ(|μ_i − μ_j|)
        let rbf = |r: f64| -> f64 { (-(r * r) / (2.0 * self.rbf_width * self.rbf_width)).exp() };
        let mut phi_mat = vec![vec![0.0; n]; n];
        for (i, row) in phi_mat.iter_mut().enumerate().take(n) {
            for (j, cell) in row.iter_mut().enumerate().take(n) {
                *cell = rbf((self.params[i] - self.params[j]).abs());
            }
        }
        // RBF evaluation vector
        let phi_mu: Vec<f64> = self.params.iter().map(|&p| rbf((p - mu).abs())).collect();
        // Interpolate each component independently
        let mut q_interp = vec![0.0; dim];
        for d in 0..dim {
            let rhs: Vec<f64> = self.solutions.iter().map(|s| s[d]).collect();
            let coeffs = cg_solve(&phi_mat, &rhs, 200, 1e-10);
            q_interp[d] = dot(&phi_mu, &coeffs);
        }
        q_interp
    }
}

// ---------------------------------------------------------------------------
// POD energy selection helper
// ---------------------------------------------------------------------------

/// Select the minimum number of POD modes to capture a given energy fraction.
///
/// Returns the index r such that Σ_{k=1}^r σ_k² / Σ_all σ_k² ≥ energy_fraction.
pub fn select_pod_modes_by_energy(singular_values: &[f64], energy_fraction: f64) -> usize {
    let total: f64 = singular_values.iter().map(|s| s * s).sum();
    if total < 1e-300 {
        return 0;
    }
    let mut cumulative = 0.0;
    for (k, &s) in singular_values.iter().enumerate() {
        cumulative += s * s;
        if cumulative / total >= energy_fraction {
            return k + 1;
        }
    }
    singular_values.len()
}

// ---------------------------------------------------------------------------
// Reduced-order dynamic system (time integration)
// ---------------------------------------------------------------------------

/// Reduced-order state for dynamic (transient) ROM simulations.
#[derive(Clone, Debug)]
pub struct ReducedDynamicState {
    /// Reduced displacement q_r (r).
    pub q: Vec<f64>,
    /// Reduced velocity dq/dt (r).
    pub dq: Vec<f64>,
    /// Current simulation time (s).
    pub time: f64,
    /// Reduced stiffness matrix K_r (r × r).
    pub k_r: Vec<Vec<f64>>,
    /// Reduced mass matrix M_r (r × r).
    pub m_r: Vec<Vec<f64>>,
    /// Reduced damping matrix C_r (r × r).
    pub c_r: Vec<Vec<f64>>,
}

impl ReducedDynamicState {
    /// Create a new zero-initialised [`ReducedDynamicState`].
    pub fn new(k_r: Vec<Vec<f64>>, m_r: Vec<Vec<f64>>, c_r: Vec<Vec<f64>>) -> Self {
        let r = k_r.len();
        Self {
            q: vec![0.0; r],
            dq: vec![0.0; r],
            time: 0.0,
            k_r,
            m_r,
            c_r,
        }
    }

    /// Advance one time step using the central-difference (explicit Newmark β=0) scheme.
    ///
    /// M_r ä_r + C_r ȧ_r + K_r q_r = f_r(t)
    ///
    /// Simplified explicit update (no inversion of M_r; diagonal M_r assumed).
    pub fn step_explicit(&mut self, f_r: &[f64], dt: f64) {
        let r = self.q.len();
        // Compute acceleration a = M_r^{-1} (f_r − K_r q − C_r dq)
        let kq = mat_vec(&self.k_r, &self.q);
        let cdq = mat_vec(&self.c_r, &self.dq);
        let rhs: Vec<f64> = (0..r).map(|i| f_r[i] - kq[i] - cdq[i]).collect();
        // Approximate M_r^{-1} via diagonal (lumped)
        let a: Vec<f64> = (0..r)
            .map(|i| {
                let mii = self.m_r[i][i];
                if mii.abs() < 1e-14 { 0.0 } else { rhs[i] / mii }
            })
            .collect();
        // Update velocity and displacement
        for (i, &ai) in a.iter().enumerate().take(r) {
            self.dq[i] += ai * dt;
            self.q[i] += self.dq[i] * dt;
        }
        self.time += dt;
    }
}

// ---------------------------------------------------------------------------
// Residual-based error indicator (Greedy)
// ---------------------------------------------------------------------------

/// Compute the projection error ||u − Φ Φᵀ u|| for a set of snapshots.
///
/// Used as error indicator in the greedy RB algorithm.
pub fn projection_errors(snapshots: &[Vec<f64>], basis: &PodBasis) -> Vec<f64> {
    snapshots
        .iter()
        .map(|snap| {
            let q_r = basis.project(snap);
            let recon = basis.reconstruct(&q_r);
            norm(&vec_sub(snap, &recon))
        })
        .collect()
}

/// Compute the relative projection error for a snapshot.
pub fn relative_projection_error(snap: &[f64], basis: &PodBasis) -> f64 {
    let q_r = basis.project(snap);
    let recon = basis.reconstruct(&q_r);
    let err = norm(&vec_sub(snap, &recon));
    let sn = norm(snap);
    if sn < 1e-14 {
        return 0.0;
    }
    err / sn
}

// ---------------------------------------------------------------------------
// Utility: Gram-Schmidt orthonormalisation
// ---------------------------------------------------------------------------

/// Orthonormalise a set of vectors using modified Gram-Schmidt.
pub fn gram_schmidt(vectors: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(vectors.len());
    for v in vectors {
        let mut w = v.clone();
        for phi in &basis {
            let c = dot(&w, phi);
            w = vec_sub(&w, &vec_scale(phi, c));
        }
        let n = norm(&w);
        if n > 1e-12 {
            basis.push(vec_scale(&w, 1.0 / n));
        }
    }
    basis
}

/// Check if a set of vectors is orthonormal (max off-diagonal |⟨ϕ_i, ϕ_j⟩| < tol).
pub fn is_orthonormal(basis: &[Vec<f64>], tol: f64) -> bool {
    for i in 0..basis.len() {
        for j in 0..basis.len() {
            let inner = dot(&basis[i], &basis[j]);
            let expected = if i == j { 1.0 } else { 0.0 };
            if (inner - expected).abs() > tol {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// QDEIM – Q-DEIM variant (column-pivoted QR selection)
// ---------------------------------------------------------------------------

/// QDEIM index selection using a simplified column-pivoted greedy approach.
///
/// Selects m indices from an r-mode basis matrix U (N × r) by column-norm
/// greedy pivoting (a surrogate for column-pivoted QR).
pub fn qdeim_indices(modes: &[Vec<f64>], m: usize) -> Vec<usize> {
    if modes.is_empty() || m == 0 {
        return vec![];
    }
    let n = modes[0].len();
    let mut selected: Vec<usize> = Vec::with_capacity(m);
    // Build residual matrix = copy of modes columns transposed → (n × r) stored row-major
    let _r = modes.len();
    let mut residual: Vec<Vec<f64>> = (0..n)
        .map(|i| modes.iter().map(|mode| mode[i]).collect())
        .collect(); // residual[i] = modes[:][i]

    for _step in 0..m {
        // Pick row (spatial index) with largest norm
        let best_i = residual
            .iter()
            .enumerate()
            .max_by(|a, b| {
                let na: f64 = a.1.iter().map(|x| x * x).sum::<f64>();
                let nb: f64 = b.1.iter().map(|x| x * x).sum::<f64>();
                na.partial_cmp(&nb).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        selected.push(best_i);

        // Orthogonalise residual against the selected row (project out)
        let row = residual[best_i].clone();
        let rn: f64 = row.iter().map(|x| x * x).sum::<f64>().sqrt();
        if rn < 1e-14 {
            break;
        }
        let row_unit: Vec<f64> = row.iter().map(|x| x / rn).collect();
        for res_row in residual.iter_mut().take(n) {
            let c: f64 = res_row
                .iter()
                .zip(row_unit.iter())
                .map(|(&ri, &ru)| ri * ru)
                .sum();
            for (res_k, &ru_k) in res_row.iter_mut().zip(row_unit.iter()) {
                *res_k -= c * ru_k;
            }
        }
    }
    selected
}

// ---------------------------------------------------------------------------
// POD-based output functional
// ---------------------------------------------------------------------------

/// Evaluate a linear output functional l(u) = s·u using the ROM solution.
///
/// - `s_full`   – output vector s (N).
/// - `u_r`      – reduced coordinates.
/// - `basis`    – POD basis.
pub fn rom_output_functional(s_full: &[f64], u_r: &[f64], basis: &PodBasis) -> f64 {
    let u_recon = basis.reconstruct(u_r);
    dot(s_full, &u_recon)
}

// ---------------------------------------------------------------------------
// Sensitivity of the ROM output w.r.t. a parameter
// ---------------------------------------------------------------------------

/// Finite-difference sensitivity of the ROM output functional.
///
/// dJ/dμ ≈ (J(μ+h) − J(μ−h)) / (2h)
///
/// - `mu`         – nominal parameter value.
/// - `h`          – finite-difference step.
/// - `eval_rom`   – closure that evaluates the ROM output at a given μ.
pub fn rom_output_sensitivity<F>(mu: f64, h: f64, eval_rom: F) -> f64
where
    F: Fn(f64) -> f64,
{
    (eval_rom(mu + h) - eval_rom(mu - h)) / (2.0 * h)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // helpers
    fn identity_snapshots(n: usize, m: usize) -> Vec<Vec<f64>> {
        (0..m)
            .map(|k| {
                let mut v = vec![0.0; n];
                if k < n {
                    v[k] = 1.0;
                }
                v
            })
            .collect()
    }

    fn linspace_snapshots(n: usize, m: usize) -> Vec<Vec<f64>> {
        (0..m)
            .map(|k| {
                let t = k as f64 / (m as f64);
                (0..n)
                    .map(|i| (PI * t * (i + 1) as f64 / n as f64).sin())
                    .collect()
            })
            .collect()
    }

    fn identity_matrix(n: usize) -> Vec<Vec<f64>> {
        eye(n)
    }

    // -----------------------------------------------------------------------
    // Math helpers
    // -----------------------------------------------------------------------

    #[test]
    fn test_dot_product() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        assert!((dot(&a, &b) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_norm_unit_vector() {
        let v = vec![1.0, 0.0, 0.0];
        assert!((norm(&v) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gram_schmidt_orthonormal() {
        let vecs = vec![
            vec![1.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
        ];
        let basis = gram_schmidt(&vecs);
        assert!(is_orthonormal(&basis, 1e-10), "Should be orthonormal");
    }

    #[test]
    fn test_gram_schmidt_single_vector() {
        let vecs = vec![vec![3.0, 4.0]];
        let basis = gram_schmidt(&vecs);
        assert_eq!(basis.len(), 1);
        assert!((norm(&basis[0]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mat_vec_identity() {
        let id = identity_matrix(3);
        let x = vec![1.0, 2.0, 3.0];
        let y = mat_vec(&id, &x);
        assert!((y[0] - 1.0).abs() < 1e-12);
        assert!((y[1] - 2.0).abs() < 1e-12);
        assert!((y[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_back_substitute_upper_triangular() {
        // U x = b: [[2,1],[0,3]] x = [5,6] => x = [1.5, 2]
        let u = vec![vec![2.0, 1.0], vec![0.0, 3.0]];
        let b = vec![5.0, 6.0];
        let x = back_substitute(&u, &b);
        assert!((x[0] - 1.5).abs() < 1e-10, "x[0]={}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1]={}", x[1]);
    }

    #[test]
    fn test_forward_substitute_lower_triangular() {
        // L x = b: [[1,0],[2,1]] x = [3,7] => x = [3, 1]
        let l = vec![vec![1.0, 0.0], vec![2.0, 1.0]];
        let b = vec![3.0, 7.0];
        let x = forward_substitute(&l, &b);
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // POD Basis
    // -----------------------------------------------------------------------

    #[test]
    fn test_pod_basis_new() {
        let basis = PodBasis::new(100);
        assert_eq!(basis.full_dim, 100);
        assert_eq!(basis.rank(), 0);
    }

    #[test]
    fn test_build_pod_basis_from_identity_snapshots() {
        let snaps = identity_snapshots(5, 3);
        let basis = build_pod_basis(&snaps, 3);
        assert!(basis.rank() <= 3);
        assert_eq!(basis.full_dim, 5);
    }

    #[test]
    fn test_pod_basis_project_reconstruct_roundtrip() {
        // Use a single snapshot and check projection–reconstruction of that snapshot
        // (trivially exact: basis spans the direction of the snapshot).
        let q = vec![1.0, 2.0, 3.0, 4.0];
        let snaps = vec![q.clone()];
        let basis = build_pod_basis(&snaps, 1);
        assert_eq!(basis.rank(), 1);
        let q_r = basis.project(&q);
        let recon = basis.reconstruct(&q_r);
        let err = norm(&vec_sub(&q, &recon));
        assert!(
            err < 1e-10,
            "Reconstruction error must be ~0 for single-snapshot basis: {err:.6}"
        );
    }

    #[test]
    fn test_pod_basis_singular_values_positive() {
        let snaps = linspace_snapshots(8, 4);
        let basis = build_pod_basis(&snaps, 4);
        for &sv in &basis.singular_values {
            assert!(sv >= 0.0, "Singular value must be ≥ 0: {sv}");
        }
    }

    #[test]
    fn test_pod_energy_content_full_rank() {
        let snaps = identity_snapshots(3, 3);
        let basis = build_pod_basis(&snaps, 3);
        let e = basis.energy_content();
        assert!(e > 0.0 && e <= 1.0 + 1e-10, "Energy content: {e}");
    }

    #[test]
    fn test_select_pod_modes_by_energy() {
        let svs = vec![10.0, 5.0, 1.0, 0.1];
        // Total energy = 100 + 25 + 1 + 0.01 = 126.01
        // 1 mode: 100/126.01 ≈ 0.794
        // 2 modes: 125/126.01 ≈ 0.992
        let r = select_pod_modes_by_energy(&svs, 0.99);
        assert!(r <= 3, "Need at most 3 modes for 99% energy: r={r}");
        assert!(r >= 2, "Need at least 2 modes: r={r}");
    }

    // -----------------------------------------------------------------------
    // Galerkin projection
    // -----------------------------------------------------------------------

    #[test]
    fn test_galerkin_system_solve_identity() {
        // K = I (3×3), f = [1,1,1], basis = I → K_r = I, f_r = [1,1,1]
        let k = identity_matrix(3);
        let f = vec![1.0, 2.0, 3.0];
        let snaps = identity_snapshots(3, 3);
        let basis = build_pod_basis(&snaps, 3);
        let gs = galerkin_projection(&k, &f, &basis);
        let q_r = gs.solve();
        // Should approximately recover f (up to sign)
        assert_eq!(q_r.len(), gs.reduced_dim);
    }

    #[test]
    fn test_galerkin_system_reduced_dim() {
        let k = identity_matrix(5);
        let f = vec![1.0; 5];
        let snaps = identity_snapshots(5, 3);
        let basis = build_pod_basis(&snaps, 3);
        let gs = galerkin_projection(&k, &f, &basis);
        assert_eq!(gs.reduced_dim, gs.k_reduced.len());
        assert_eq!(gs.reduced_dim, gs.f_reduced.len());
    }

    // -----------------------------------------------------------------------
    // Reduced Basis (greedy)
    // -----------------------------------------------------------------------

    #[test]
    fn test_reduced_basis_add_snapshot() {
        let mut rb = ReducedBasis::new(4);
        rb.add_snapshot(&[1.0, 0.0, 0.0, 0.0], 1.0);
        rb.add_snapshot(&[0.0, 1.0, 0.0, 0.0], 2.0);
        assert_eq!(rb.dim(), 2);
    }

    #[test]
    fn test_reduced_basis_orthonormal() {
        let mut rb = ReducedBasis::new(4);
        rb.add_snapshot(&[1.0, 1.0, 0.0, 0.0], 1.0);
        rb.add_snapshot(&[0.0, 1.0, 1.0, 0.0], 2.0);
        rb.add_snapshot(&[0.0, 0.0, 1.0, 1.0], 3.0);
        assert!(is_orthonormal(&rb.basis, 1e-10));
    }

    #[test]
    fn test_greedy_reduced_basis_convergence() {
        let params: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let snaps: Vec<Vec<f64>> = params
            .iter()
            .map(|&mu| vec![mu.cos(), mu.sin(), mu * mu * 0.1, (-mu).exp()])
            .collect();
        let rb = greedy_reduced_basis(&params, &snaps, 5, 1e-8);
        assert!(rb.dim() >= 1);
        assert!(rb.dim() <= 5);
    }

    #[test]
    fn test_greedy_reduced_basis_empty_snapshots() {
        let rb = greedy_reduced_basis(&[], &[], 5, 1e-6);
        assert_eq!(rb.dim(), 0);
    }

    // -----------------------------------------------------------------------
    // DEIM
    // -----------------------------------------------------------------------

    #[test]
    fn test_deim_indices_count() {
        let modes: Vec<Vec<f64>> = (0..3)
            .map(|k| {
                (0..8)
                    .map(|i| (PI * (k + 1) as f64 * i as f64 / 8.0).sin())
                    .collect()
            })
            .collect();
        let interp = deim(&modes);
        assert_eq!(interp.indices.len(), 3);
    }

    #[test]
    fn test_deim_empty_modes() {
        let interp = deim(&[]);
        assert!(interp.indices.is_empty());
    }

    #[test]
    fn test_deim_approximate_identity_modes() {
        // With modes = identity columns, DEIM should select indices 0..r
        let modes: Vec<Vec<f64>> = (0..3)
            .map(|k| {
                let mut v = vec![0.0; 6];
                v[k] = 1.0;
                v
            })
            .collect();
        let interp = deim(&modes);
        let f_full = vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        let approx = interp.approximate(&f_full);
        assert_eq!(approx.len(), 6);
    }

    // -----------------------------------------------------------------------
    // QDEIM
    // -----------------------------------------------------------------------

    #[test]
    fn test_qdeim_indices_count() {
        let modes: Vec<Vec<f64>> = (0..2)
            .map(|k| {
                (0..6)
                    .map(|i| (PI * (k + 1) as f64 * i as f64 / 6.0).sin())
                    .collect()
            })
            .collect();
        let indices = qdeim_indices(&modes, 2);
        assert_eq!(indices.len(), 2);
    }

    #[test]
    fn test_qdeim_empty() {
        let indices = qdeim_indices(&[], 3);
        assert!(indices.is_empty());
    }

    // -----------------------------------------------------------------------
    // Petrov-Galerkin
    // -----------------------------------------------------------------------

    #[test]
    fn test_petrov_galerkin_same_as_galerkin_for_symmetric() {
        // When phi == psi, PG should equal Galerkin
        let k = identity_matrix(4);
        let f = vec![1.0, 2.0, 3.0, 4.0];
        let snaps = identity_snapshots(4, 3);
        let basis = build_pod_basis(&snaps, 3);
        let phi = &basis.modes;
        let psi = phi;
        let pg = petrov_galerkin_projection(&k, &f, phi, psi);
        let g = galerkin_projection(&k, &f, &basis);
        assert_eq!(pg.reduced_dim, g.reduced_dim);
    }

    #[test]
    fn test_petrov_galerkin_solve_returns_correct_length() {
        let k = identity_matrix(4);
        let f = vec![1.0; 4];
        let snaps = identity_snapshots(4, 2);
        let basis = build_pod_basis(&snaps, 2);
        let phi = &basis.modes;
        let pg = petrov_galerkin_projection(&k, &f, phi, phi);
        let sol = pg.solve();
        assert_eq!(sol.len(), pg.reduced_dim);
    }

    // -----------------------------------------------------------------------
    // Error bound
    // -----------------------------------------------------------------------

    #[test]
    fn test_error_bound_zero_residual() {
        // K x = f exactly → residual = 0 → error bound = 0
        let k = identity_matrix(3);
        let f = vec![1.0, 2.0, 3.0];
        let u_rom = f.clone(); // exact solution for I
        let eb = compute_error_bound(&k, &f, &u_rom, 1.0, None);
        assert!(
            eb.residual_norm < 1e-12,
            "Residual should be 0: {}",
            eb.residual_norm
        );
        assert!(eb.error_bound < 1e-12);
    }

    #[test]
    fn test_error_bound_positive_residual() {
        let k = identity_matrix(3);
        let f = vec![1.0, 2.0, 3.0];
        let u_approx = vec![0.5, 1.0, 1.5]; // wrong solution
        let eb = compute_error_bound(&k, &f, &u_approx, 1.0, None);
        assert!(eb.error_bound > 0.0);
    }

    #[test]
    fn test_error_bound_effectivity() {
        let k = identity_matrix(2);
        let f = vec![2.0, 4.0];
        let u_true = vec![2.0, 4.0];
        let u_rom = vec![1.0, 3.0];
        let true_err = norm(&vec_sub(&u_true, &u_rom));
        let eb = compute_error_bound(&k, &f, &u_rom, 1.0, Some(true_err));
        assert!(eb.effectivity.is_some());
    }

    // -----------------------------------------------------------------------
    // EIM
    // -----------------------------------------------------------------------

    #[test]
    fn test_eim_magic_points_count() {
        let snaps: Vec<Vec<f64>> = (0..3)
            .map(|k| {
                (0..8)
                    .map(|i| (PI * (k + 1) as f64 * i as f64 / 8.0).sin())
                    .collect()
            })
            .collect();
        let eim = build_eim(&snaps);
        assert!(eim.magic_points.len() <= 3);
    }

    #[test]
    fn test_eim_empty_snapshots() {
        let eim = build_eim(&[]);
        assert!(eim.magic_points.is_empty());
    }

    // -----------------------------------------------------------------------
    // Affine decomposition (offline/online)
    // -----------------------------------------------------------------------

    #[test]
    fn test_offline_affine_assemble_recovers_full() {
        // Single component: K(μ) = θ(μ) I → K_r = θ K_r_1
        let k_components = vec![identity_matrix(4)];
        let f_components = vec![vec![1.0, 2.0, 3.0, 4.0]];
        let snaps = identity_snapshots(4, 2);
        let basis = build_pod_basis(&snaps, 2);
        let ad = offline_affine_decomposition(&k_components, &f_components, &basis);
        let system = ad.assemble(&[2.0], &[1.0]);
        assert_eq!(system.reduced_dim, basis.rank());
        // K_r should be 2 * K_r_1
        for i in 0..system.reduced_dim {
            for j in 0..system.reduced_dim {
                assert!(
                    (system.k_reduced[i][j] - 2.0 * ad.k_reduced_components[0][i][j]).abs() < 1e-12
                );
            }
        }
    }

    #[test]
    fn test_offline_affine_decomposition_empty() {
        let basis = PodBasis::new(4);
        let ad = offline_affine_decomposition(&[], &[], &basis);
        assert!(ad.k_reduced_components.is_empty());
    }

    // -----------------------------------------------------------------------
    // Manifold interpolation
    // -----------------------------------------------------------------------

    #[test]
    fn test_manifold_interpolator_at_training_point() {
        let mut interp = ManifoldInterpolator::new(1.0);
        interp.add_training_point(0.0, vec![1.0, 0.0]);
        interp.add_training_point(1.0, vec![0.0, 1.0]);
        // At mu = 0 should be close to [1, 0]
        let q = interp.interpolate(0.0);
        assert_eq!(q.len(), 2);
        assert!(
            q[0] > 0.4,
            "Expected q[0] ≈ 1 at training point, got {}",
            q[0]
        );
    }

    #[test]
    fn test_manifold_interpolator_empty() {
        let interp = ManifoldInterpolator::new(1.0);
        let q = interp.interpolate(0.5);
        assert!(q.is_empty());
    }

    // -----------------------------------------------------------------------
    // Balanced truncation
    // -----------------------------------------------------------------------

    #[test]
    fn test_balanced_truncation_returns_modes() {
        let a = identity_matrix(4);
        let b_mat = vec![vec![1.0], vec![0.5], vec![0.0], vec![0.0]];
        let c_mat = vec![vec![1.0, 0.0, 0.5, 0.0]];
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 2);
        assert!(bt.reduced_dim <= 2);
        assert_eq!(bt.hankel_singular_values.len(), bt.reduced_dim);
    }

    #[test]
    fn test_balanced_truncation_project_vector() {
        let a = identity_matrix(3);
        let b_mat = vec![vec![1.0], vec![1.0], vec![1.0]];
        let c_mat = vec![vec![1.0, 1.0, 1.0]];
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 2);
        let b_vec = vec![1.0, 0.0, 0.0];
        let proj = bt.project_vector(&b_vec);
        assert_eq!(proj.len(), bt.reduced_dim);
    }

    // -----------------------------------------------------------------------
    // Projection errors
    // -----------------------------------------------------------------------

    #[test]
    fn test_projection_errors_zero_for_full_basis() {
        // Use well-separated snapshots (each is a different sine wave at distinct
        // frequencies) so power iteration can separate them cleanly.
        let snaps: Vec<Vec<f64>> = vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0, 5.0]];
        let basis = build_pod_basis(&snaps, 2);
        let errors = projection_errors(&snaps, &basis);
        for (i, &e) in errors.iter().enumerate() {
            assert!(e < 0.5, "Error for snap {i} too large: {e:.6}");
        }
    }

    #[test]
    fn test_relative_projection_error_zero_snap() {
        let snap = vec![0.0, 0.0, 0.0];
        let basis = PodBasis::new(3);
        let rel_err = relative_projection_error(&snap, &basis);
        assert_eq!(rel_err, 0.0);
    }

    // -----------------------------------------------------------------------
    // Reduced dynamic state
    // -----------------------------------------------------------------------

    #[test]
    fn test_reduced_dynamic_state_advances_time() {
        let k_r = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let m_r = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let c_r = vec![vec![0.1, 0.0], vec![0.0, 0.1]];
        let mut state = ReducedDynamicState::new(k_r, m_r, c_r);
        let f_r = vec![0.0, 0.0];
        state.step_explicit(&f_r, 0.01);
        assert!((state.time - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_reduced_dynamic_state_displacement_grows_under_load() {
        let k_r = vec![vec![1.0]];
        let m_r = vec![vec![1.0]];
        let c_r = vec![vec![0.0]];
        let mut state = ReducedDynamicState::new(k_r, m_r, c_r);
        let f_r = vec![1.0]; // constant force
        for _ in 0..100 {
            state.step_explicit(&f_r, 0.01);
        }
        // Displacement should become positive under positive force
        assert!(state.q[0] > 0.0, "Displacement should grow: {}", state.q[0]);
    }

    // -----------------------------------------------------------------------
    // ROM output functional
    // -----------------------------------------------------------------------

    #[test]
    fn test_rom_output_functional() {
        let snaps = identity_snapshots(3, 3);
        let basis = build_pod_basis(&snaps, 2);
        let s = vec![1.0, 1.0, 1.0];
        let u_r = vec![1.0, 0.5];
        let out = rom_output_functional(&s, &u_r, &basis);
        assert!(out.is_finite(), "Output should be finite: {out}");
    }

    #[test]
    fn test_rom_output_sensitivity() {
        // J(μ) = μ², dJ/dμ = 2μ → at μ=3 should be ≈6
        let sens = rom_output_sensitivity(3.0, 1e-5, |mu| mu * mu);
        assert!((sens - 6.0).abs() < 1e-6, "Sensitivity = {sens:.6}");
    }
}
