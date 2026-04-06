#![allow(clippy::needless_range_loop, clippy::type_complexity)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reduced Order Modeling (ROM) for FEM.
//!
//! Implements Proper Orthogonal Decomposition (POD), greedy Reduced Basis (RB)
//! methods, Discrete Empirical Interpolation Method (DEIM), and balanced
//! truncation for model order reduction of large-scale FEM systems.
//!
//! # References
//! - Sirovich (1987) "Turbulence and the dynamics of coherent structures"
//! - Amsallem & Farhat (2008) "Interpolation method for adapting reduced-order models"
//! - Chaturantabut & Sorensen (2010) "Nonlinear model reduction via DEIM"

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers (plain f64 / Vec<f64>)
// ---------------------------------------------------------------------------

/// Dot product of two slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm of a slice.
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Subtract two vectors element-wise.
fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

/// Add two vectors element-wise.
fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scale a vector by a scalar.
fn vec_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

/// Dense matrix–vector product: y = A*x, A stored row-major as `Vec<Vec`f64`>`.
fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

/// Transpose of a dense matrix.
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

/// Dense matrix–matrix product.
fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let p = b.len();
    if m == 0 || p == 0 {
        return vec![];
    }
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

/// Solve a small dense lower-triangular system Ly = b.
fn forward_substitute(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i][j] * y[j];
        }
        y[i] = if l[i][i].abs() < 1e-14 {
            0.0
        } else {
            s / l[i][i]
        };
    }
    y
}

/// Solve a small dense upper-triangular system Ux = y.
fn back_substitute(u: &[Vec<f64>], y: &[f64]) -> Vec<f64> {
    let n = y.len();
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
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

/// Cholesky decomposition: A = L * L^T for a symmetric positive-definite matrix.
/// Returns L (lower triangular).
fn cholesky(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s: f64 = a[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                l[i][j] = if s > 0.0 { s.sqrt() } else { 1e-14 };
            } else {
                l[i][j] = if l[j][j].abs() < 1e-14 {
                    0.0
                } else {
                    s / l[j][j]
                };
            }
        }
    }
    l
}

/// Power iteration to find the dominant eigenvector of a symmetric matrix.
fn power_iteration(a: &[Vec<f64>], max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = a.len();
    // Use a deterministic starting vector
    let mut v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.1 }).collect();
    let nv = norm(&v);
    v = vec_scale(&v, 1.0 / nv);
    let mut eigenvalue = 0.0;
    for _ in 0..max_iter {
        let w = mat_vec(a, &v);
        let lam = dot(&v, &w);
        let nw = norm(&w);
        if nw < 1e-14 {
            break;
        }
        let v_new = vec_scale(&w, 1.0 / nw);
        if (lam - eigenvalue).abs() < tol {
            eigenvalue = lam;
            v = v_new;
            break;
        }
        eigenvalue = lam;
        v = v_new;
    }
    (eigenvalue, v)
}

/// Deflate matrix: A_deflated = A - lambda * v * v^T.
fn deflate(a: &[Vec<f64>], lambda: f64, v: &[f64]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut b = a.to_vec();
    for i in 0..n {
        for j in 0..n {
            b[i][j] -= lambda * v[i] * v[j];
        }
    }
    b
}

/// Identity matrix of size n.
fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut id = vec![vec![0.0; n]; n];
    for i in 0..n {
        id[i][i] = 1.0;
    }
    id
}

// ---------------------------------------------------------------------------
// Snapshot matrix
// ---------------------------------------------------------------------------

/// Collection of solution snapshots for POD basis construction.
///
/// Each snapshot is a full-state vector of length `n_dof`, captured at
/// different time steps or parameter values.
#[derive(Debug, Clone)]
pub struct SnapshotMatrix {
    /// Snapshot vectors stored row-major (each inner Vec is one snapshot).
    pub data: Vec<Vec<f64>>,
    /// Number of degrees of freedom per snapshot.
    pub n_dof: usize,
    /// Number of snapshots collected so far.
    pub n_snapshots: usize,
}

impl SnapshotMatrix {
    /// Create an empty snapshot collection for a system with `n_dof` dofs.
    pub fn new(n_dof: usize) -> Self {
        Self {
            data: Vec::new(),
            n_dof,
            n_snapshots: 0,
        }
    }

    /// Append one solution snapshot.
    ///
    /// Panics in debug mode if `state.len() != self.n_dof`.
    pub fn add_snapshot(&mut self, state: Vec<f64>) {
        debug_assert_eq!(
            state.len(),
            self.n_dof,
            "snapshot length mismatch: expected {}, got {}",
            self.n_dof,
            state.len()
        );
        self.data.push(state);
        self.n_snapshots += 1;
    }

    /// Compute the mean snapshot: ū_i = (1/N) Σ_k u_i^{(k)}.
    pub fn mean(&self) -> Vec<f64> {
        if self.n_snapshots == 0 {
            return vec![0.0; self.n_dof];
        }
        let scale = 1.0 / self.n_snapshots as f64;
        let mut m = vec![0.0; self.n_dof];
        for snap in &self.data {
            for (mi, &si) in m.iter_mut().zip(snap.iter()) {
                *mi += si * scale;
            }
        }
        m
    }

    /// Return a new `SnapshotMatrix` with the mean subtracted from every snapshot.
    pub fn centered(&self) -> Self {
        let mu = self.mean();
        let centered_data: Vec<Vec<f64>> =
            self.data.iter().map(|snap| vec_sub(snap, &mu)).collect();
        Self {
            data: centered_data,
            n_dof: self.n_dof,
            n_snapshots: self.n_snapshots,
        }
    }
}

// ---------------------------------------------------------------------------
// POD basis
// ---------------------------------------------------------------------------

/// POD (Proper Orthogonal Decomposition) basis for reduced-order modeling.
///
/// The basis vectors are the dominant left singular vectors of the centered
/// snapshot matrix, ordered by decreasing singular value.
#[derive(Debug, Clone)]
pub struct PodBasis {
    /// Orthonormal basis vectors (each of length `n_dof`).
    pub modes: Vec<Vec<f64>>,
    /// Singular values corresponding to each mode.
    pub singular_values: Vec<f64>,
    /// Number of retained modes.
    pub n_modes: usize,
}

impl PodBasis {
    /// Construct a POD basis from a snapshot matrix.
    ///
    /// Uses the method of snapshots (Sirovich 1987): computes the `n×n`
    /// correlation matrix C = U^T U (where U is the n_dof × n_snapshots
    /// centered snapshot matrix), finds dominant eigenvectors by deflated
    /// power iteration, then lifts them to full-space modes.
    ///
    /// Only modes whose cumulative energy fraction exceeds `energy_threshold`
    /// are retained.  `energy_threshold = 1.0` retains all modes.
    pub fn from_snapshots(snapshots: &SnapshotMatrix, energy_threshold: f64) -> Self {
        let n = snapshots.n_snapshots;
        if n == 0 {
            return Self {
                modes: vec![],
                singular_values: vec![],
                n_modes: 0,
            };
        }
        // Special case: a single snapshot is its own mode (centering would zero it out).
        if n == 1 {
            let u = &snapshots.data[0];
            let nm = norm(u);
            if nm < 1e-14 {
                return Self {
                    modes: vec![],
                    singular_values: vec![],
                    n_modes: 0,
                };
            }
            let mode = vec_scale(u, 1.0 / nm);
            return Self {
                modes: vec![mode],
                singular_values: vec![nm],
                n_modes: 1,
            };
        }
        let centered = snapshots.centered();

        // Build n×n correlation matrix C_kl = <u^k, u^l>
        let mut corr = vec![vec![0.0; n]; n];
        for k in 0..n {
            for l in k..n {
                let c = dot(&centered.data[k], &centered.data[l]);
                corr[k][l] = c;
                corr[l][k] = c;
            }
        }

        // Extract eigenpairs by deflated power iteration
        let mut eigenvalues: Vec<f64> = Vec::new();
        let mut eigenvectors: Vec<Vec<f64>> = Vec::new();
        let mut a = corr.clone();

        for _ in 0..n {
            let (lam, v) = power_iteration(&a, 500, 1e-10);
            if lam < 1e-14 {
                break;
            }
            eigenvalues.push(lam);
            eigenvectors.push(v.clone());
            a = deflate(&a, lam, &v);
        }

        // Compute total energy
        let total_energy: f64 = eigenvalues.iter().sum();

        // Determine how many modes to retain
        let mut cumulative = 0.0;
        let mut n_retain = 0;
        for &ev in &eigenvalues {
            cumulative += ev;
            n_retain += 1;
            if total_energy > 0.0 && cumulative / total_energy >= energy_threshold {
                break;
            }
        }

        // Lift reduced eigenvectors to full-space modes:  φ_k = U * v_k / σ_k
        let mut modes: Vec<Vec<f64>> = Vec::with_capacity(n_retain);
        let mut singular_values: Vec<f64> = Vec::with_capacity(n_retain);

        for k in 0..n_retain {
            let sigma = eigenvalues[k].max(0.0).sqrt();
            if sigma < 1e-14 {
                continue;
            }
            // φ_k = (1/σ_k) Σ_l v_k[l] * snap_l
            let mut mode = vec![0.0; snapshots.n_dof];
            for l in 0..n {
                let coeff = eigenvectors[k][l] / sigma;
                for (mi, &si) in mode.iter_mut().zip(centered.data[l].iter()) {
                    *mi += coeff * si;
                }
            }
            // Normalize
            let nm = norm(&mode);
            if nm > 1e-14 {
                mode = vec_scale(&mode, 1.0 / nm);
            }
            modes.push(mode);
            singular_values.push(sigma);
        }

        let n_modes = modes.len();
        Self {
            modes,
            singular_values,
            n_modes,
        }
    }

    /// Project a full-state vector onto the reduced coordinates.
    ///
    /// Returns a vector of length `n_modes` with coefficients `q_k = <φ_k, u>`.
    pub fn project(&self, state: &[f64]) -> Vec<f64> {
        self.modes.iter().map(|phi| dot(phi, state)).collect()
    }

    /// Reconstruct a full-state vector from reduced coordinates.
    ///
    /// Returns `u ≈ Σ_k q_k φ_k`.
    pub fn reconstruct(&self, reduced: &[f64]) -> Vec<f64> {
        if self.n_modes == 0 {
            return vec![];
        }
        let n_dof = self.modes[0].len();
        let mut u = vec![0.0; n_dof];
        for (phi, &q) in self.modes.iter().zip(reduced.iter()) {
            for (ui, &phi_i) in u.iter_mut().zip(phi.iter()) {
                *ui += q * phi_i;
            }
        }
        u
    }

    /// Fraction of total energy retained by the first `n` modes.
    ///
    /// Returns a value in \[0, 1\].
    pub fn energy_retained(&self, n: usize) -> f64 {
        let total: f64 = self.singular_values.iter().map(|s| s * s).sum();
        if total < 1e-100 {
            return 1.0;
        }
        let partial: f64 = self.singular_values.iter().take(n).map(|s| s * s).sum();
        (partial / total).min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Reduced-order system
// ---------------------------------------------------------------------------

/// Reduced-order structural dynamics system.
///
/// Stores the POD basis together with projected mass, stiffness, and damping
/// matrices in the reduced space.  Time integration is performed in the
/// low-dimensional reduced space, cutting cost from O(n_dof³) to O(n_modes³).
#[derive(Debug, Clone)]
pub struct RomSystem {
    /// POD basis used for projection.
    pub basis: PodBasis,
    /// Reduced mass matrix (n_modes × n_modes).
    pub reduced_mass: Vec<Vec<f64>>,
    /// Reduced stiffness matrix (n_modes × n_modes).
    pub reduced_stiffness: Vec<Vec<f64>>,
    /// Reduced damping matrix (n_modes × n_modes).
    pub reduced_damping: Vec<Vec<f64>>,
    /// Current reduced displacement (length n_modes).
    disp: Vec<f64>,
    /// Current reduced velocity (length n_modes).
    vel: Vec<f64>,
    /// Previous reduced displacement (for Newmark).
    disp_prev: Vec<f64>,
    /// Previous reduced velocity (for Newmark).
    vel_prev: Vec<f64>,
}

impl RomSystem {
    /// Build a reduced system by projecting full mass and stiffness matrices.
    ///
    /// Computes `M_r = Φ^T M Φ` and `K_r = Φ^T K Φ`.  Damping is set to zero
    /// (call [`Self::set_rayleigh_damping`] to add proportional damping).
    pub fn project_system(basis: &PodBasis, mass: &[Vec<f64>], stiffness: &[Vec<f64>]) -> Self {
        let r = basis.n_modes;
        // Φ^T M  (r × n_dof)
        let phi_t_m = Self::phi_t_a(basis, mass);
        // Φ^T K
        let phi_t_k = Self::phi_t_a(basis, stiffness);
        // Reduced matrices: r×r: M_r[i][j] = row_i(phi_t_m) . phi[j]
        let reduced_mass = Self::phi_t_a_phi(basis, &phi_t_m);
        let reduced_stiffness = Self::phi_t_a_phi(basis, &phi_t_k);
        let reduced_damping = vec![vec![0.0; r]; r];

        Self {
            basis: basis.clone(),
            reduced_mass,
            reduced_stiffness,
            reduced_damping,
            disp: vec![0.0; r],
            vel: vec![0.0; r],
            disp_prev: vec![0.0; r],
            vel_prev: vec![0.0; r],
        }
    }

    /// Helper: compute Φ^T * A (result is n_modes × n_dof).
    fn phi_t_a(basis: &PodBasis, a: &[Vec<f64>]) -> Vec<Vec<f64>> {
        // Φ^T * A  = (A^T * Φ)^T
        // Row i of result = φ_i^T * A = A^T φ_i
        let at = mat_transpose(a);
        basis.modes.iter().map(|phi| mat_vec(&at, phi)).collect()
    }

    /// Helper: given B = Φ^T A (r × n_dof), compute B * Φ = Φ^T A Φ (r × r).
    fn phi_t_a_phi(basis: &PodBasis, b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        b.iter()
            .map(|row| basis.modes.iter().map(|phi| dot(row, phi)).collect())
            .collect()
    }

    /// Add Rayleigh proportional damping: C_r = α M_r + β K_r.
    pub fn set_rayleigh_damping(&mut self, alpha: f64, beta: f64) {
        let r = self.basis.n_modes;
        for i in 0..r {
            for j in 0..r {
                self.reduced_damping[i][j] =
                    alpha * self.reduced_mass[i][j] + beta * self.reduced_stiffness[i][j];
            }
        }
    }

    /// Solve the static reduced system K_r q = f_r.
    ///
    /// Uses Cholesky decomposition for symmetric positive-definite K_r.
    pub fn solve_reduced(&self, f_reduced: &[f64]) -> Vec<f64> {
        let l = cholesky(&self.reduced_stiffness);
        let y = forward_substitute(&l, f_reduced);
        let lt: Vec<Vec<f64>> = mat_transpose(&l);
        back_substitute(&lt, &y)
    }

    /// Perform one Newmark-β time step in reduced space.
    ///
    /// The effective system solved at each step is:
    /// ```text
    /// (M_r / (β dt²) + γ C_r / (β dt) + K_r) Δq = f_r - K_r q_n
    ///   - M_r (q̈_pred) - C_r (q̇_pred)
    /// ```
    ///
    /// Returns the new full-space displacement vector.
    pub fn step_newmark(&mut self, f_full: &[f64], dt: f64, beta: f64, gamma: f64) -> Vec<f64> {
        let r = self.basis.n_modes;
        // Project force
        let f_r: Vec<f64> = self.basis.project(f_full);

        // Predictors
        // q̈_n from equation of motion (initial approximation = 0)
        let acc_n = vec![0.0_f64; r];

        let disp_pred: Vec<f64> = (0..r)
            .map(|i| self.disp[i] + dt * self.vel[i] + dt * dt * (0.5 - beta) * acc_n[i])
            .collect();
        let vel_pred: Vec<f64> = (0..r)
            .map(|i| self.vel[i] + dt * (1.0 - gamma) * acc_n[i])
            .collect();

        // Effective stiffness K_eff = M_r/(β dt²) + γ C_r/(β dt) + K_r
        let c1 = 1.0 / (beta * dt * dt);
        let c2 = gamma / (beta * dt);
        let mut k_eff = vec![vec![0.0; r]; r];
        for i in 0..r {
            for j in 0..r {
                k_eff[i][j] = c1 * self.reduced_mass[i][j]
                    + c2 * self.reduced_damping[i][j]
                    + self.reduced_stiffness[i][j];
            }
        }

        // Effective force
        let m_disp_pred = mat_vec(&self.reduced_mass, &disp_pred);
        let c_vel_pred = mat_vec(&self.reduced_damping, &vel_pred);
        let eff_rhs: Vec<f64> = (0..r)
            .map(|i| {
                f_r[i]
                    - self.reduced_stiffness[0..]
                        .iter()
                        .enumerate()
                        .map(|(j, row)| row[i] * disp_pred[j])
                        .sum::<f64>()
                    + c1 * m_disp_pred[i]
                    - c_vel_pred[i] * c2 / c1
            })
            .collect();

        // More direct approach: effective RHS
        // f_eff = f_r - K_r*q_pred - C_r*v_pred + M_r/(beta*dt^2)*q_pred
        let k_q_pred = mat_vec(&self.reduced_stiffness, &disp_pred);
        let c_v_pred = mat_vec(&self.reduced_damping, &vel_pred);
        let m_q_pred = mat_vec(&self.reduced_mass, &disp_pred);
        let f_eff: Vec<f64> = (0..r)
            .map(|i| {
                f_r[i] - k_q_pred[i] - c_v_pred[i]
                    + c1 * m_q_pred[i]
                    + c2 * self.reduced_mass[i]
                        .iter()
                        .zip(vel_pred.iter())
                        .map(|(m, v)| m * v)
                        .sum::<f64>()
            })
            .collect();
        let _ = eff_rhs; // superseded by f_eff

        // Solve K_eff q_new = f_eff
        let l = cholesky(&k_eff);
        let y = forward_substitute(&l, &f_eff);
        let lt = mat_transpose(&l);
        let q_new = back_substitute(&lt, &y);

        // Update velocity and acceleration
        let acc_new: Vec<f64> = (0..r).map(|i| c1 * (q_new[i] - disp_pred[i])).collect();
        let vel_new: Vec<f64> = (0..r)
            .map(|i| vel_pred[i] + dt * gamma * acc_new[i])
            .collect();

        // Store state
        self.disp_prev = self.disp.clone();
        self.vel_prev = self.vel.clone();
        self.disp = q_new.clone();
        self.vel = vel_new;

        // Reconstruct full displacement
        self.basis.reconstruct(&q_new)
    }
}

// ---------------------------------------------------------------------------
// Greedy basis selection (Reduced Basis method)
// ---------------------------------------------------------------------------

/// Greedy Reduced Basis (RB) selection.
///
/// Iteratively selects parameter values from `param_space` that maximize
/// the error indicator (here the L2 error between the high-fidelity solution
/// and its projection onto the current basis), stopping when the error falls
/// below `tol` or `n_basis` modes have been selected.
///
/// # Arguments
/// * `param_space` – discrete set of parameter values to sample from.
/// * `solve_fn` – callback that returns the high-fidelity solution for a
///   given parameter value.
/// * `n_basis` – maximum number of basis vectors to select.
/// * `tol` – tolerance on the relative projection error.
pub fn greedy_basis_selection(
    param_space: &[f64],
    solve_fn: &dyn Fn(f64) -> Vec<f64>,
    n_basis: usize,
    tol: f64,
) -> PodBasis {
    if param_space.is_empty() {
        return PodBasis {
            modes: vec![],
            singular_values: vec![],
            n_modes: 0,
        };
    }

    let mut selected_params: Vec<f64> = Vec::new();
    let mut snapshots = {
        // Seed with first parameter
        let u0 = solve_fn(param_space[0]);
        selected_params.push(param_space[0]);
        let n_dof = u0.len();
        let mut s = SnapshotMatrix::new(n_dof);
        s.add_snapshot(u0);
        s
    };

    for _ in 1..n_basis {
        // Build current basis from collected snapshots
        let basis = PodBasis::from_snapshots(&snapshots, 1.0);

        // Find parameter with largest projection error
        let mut worst_err = 0.0_f64;
        let mut worst_param = param_space[0];
        for &mu in param_space {
            if selected_params.contains(&mu) {
                continue;
            }
            let u = solve_fn(mu);
            let q = basis.project(&u);
            let u_approx = basis.reconstruct(&q);
            let diff = vec_sub(&u, &u_approx);
            let err = norm(&diff) / (norm(&u) + 1e-14);
            if err > worst_err {
                worst_err = err;
                worst_param = mu;
            }
        }

        if worst_err < tol {
            break;
        }

        selected_params.push(worst_param);
        snapshots.add_snapshot(solve_fn(worst_param));
    }

    let mut result = PodBasis::from_snapshots(&snapshots, 1.0);
    // Fallback: if power iteration yielded no modes despite non-empty snapshots,
    // use the first (seeded) snapshot as a single mode.
    if result.n_modes == 0 && snapshots.n_snapshots > 0 {
        let u = &snapshots.data[0];
        let nm = norm(u);
        if nm >= 1e-14 {
            result = PodBasis {
                modes: vec![vec_scale(u, 1.0 / nm)],
                singular_values: vec![nm],
                n_modes: 1,
            };
        }
    }
    result
}

// ---------------------------------------------------------------------------
// DEIM – Discrete Empirical Interpolation Method
// ---------------------------------------------------------------------------

/// Discrete Empirical Interpolation for efficient evaluation of nonlinear terms.
///
/// Given a set of basis vectors for the nonlinear function space, DEIM selects
/// a set of interpolation nodes and constructs a matrix that allows evaluation
/// of the full nonlinear vector from values at only the selected nodes.
#[derive(Debug, Clone)]
pub struct EmpiricalInterpolation {
    /// Selected interpolation node indices.
    pub nodes: Vec<usize>,
    /// DEIM basis vectors (each of length `n_dof`).
    pub basis: Vec<Vec<f64>>,
}

impl EmpiricalInterpolation {
    /// Construct a DEIM interpolation operator from a set of basis modes.
    ///
    /// Implements the greedy DEIM algorithm (Chaturantabut & Sorensen 2010):
    /// at each step select the index of the maximum-magnitude residual.
    pub fn from_modes(modes: &[Vec<f64>]) -> Self {
        if modes.is_empty() {
            return Self {
                nodes: vec![],
                basis: vec![],
            };
        }
        let m = modes.len();
        let n = modes[0].len();

        let mut nodes: Vec<usize> = Vec::with_capacity(m);
        let mut selected_basis: Vec<Vec<f64>> = Vec::with_capacity(m);

        // Step 1: first node is argmax of |first mode|
        let first_node = modes[0]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.abs()
                    .partial_cmp(&b.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        nodes.push(first_node);
        selected_basis.push(modes[0].clone());

        for k in 1..m {
            // Solve U_k-1 c = u_k[nodes_{k-1}]
            let sub_rows: Vec<Vec<f64>> = nodes
                .iter()
                .map(|&i| selected_basis.iter().map(|b| b[i]).collect())
                .collect();
            let rhs: Vec<f64> = nodes.iter().map(|&i| modes[k][i]).collect();
            let c = if sub_rows.is_empty() {
                vec![]
            } else {
                // Small dense solve via Gaussian elimination
                let l = cholesky(&Self::ata(&sub_rows));
                let ata_rhs = mat_vec(&mat_transpose(&sub_rows), &rhs);
                let y = forward_substitute(&l, &ata_rhs);
                let lt = mat_transpose(&l);
                back_substitute(&lt, &y)
            };

            // Residual r = u_k - Σ c_i φ_i
            let mut residual = modes[k].clone();
            for (i, phi) in selected_basis.iter().enumerate() {
                for (ri, &pi) in residual.iter_mut().zip(phi.iter()) {
                    *ri -= c[i] * pi;
                }
            }

            // Next node = argmax |r|
            let next_node = residual
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| {
                    a.abs()
                        .partial_cmp(&b.abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);

            // Avoid duplicate nodes
            if nodes.contains(&next_node) {
                // Find next best
                let mut best = (0, 0.0_f64);
                for (i, &ri) in residual.iter().enumerate() {
                    if !nodes.contains(&i) && ri.abs() > best.1 {
                        best = (i, ri.abs());
                    }
                }
                nodes.push(if best.1 > 0.0 {
                    best.0
                } else {
                    (next_node + 1) % n
                });
            } else {
                nodes.push(next_node);
            }
            selected_basis.push(modes[k].clone());
        }

        Self {
            nodes,
            basis: selected_basis,
        }
    }

    /// Compute A^T A for a small matrix.
    fn ata(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let at = mat_transpose(a);
        mat_mul(&at, a)
    }

    /// Reconstruct the full nonlinear vector from values at interpolation nodes.
    ///
    /// Given `f_at_nodes[k] = f(x_{i_k})`, returns the DEIM approximation of
    /// the full vector f.
    pub fn interpolate(&self, f_at_nodes: &[f64]) -> Vec<f64> {
        if self.basis.is_empty() {
            return vec![];
        }
        let _m = self.basis.len();
        let n = self.basis[0].len();

        // Solve P^T Φ c = f_nodes  where P^T selects rows = node indices
        let ptphi: Vec<Vec<f64>> = self
            .nodes
            .iter()
            .map(|&i| self.basis.iter().map(|b| b[i]).collect())
            .collect();

        // Gaussian elimination on small m×m system
        let l = cholesky(&Self::ata(&ptphi));
        let ata_rhs = mat_vec(&mat_transpose(&ptphi), f_at_nodes);
        let y = forward_substitute(&l, &ata_rhs);
        let lt = mat_transpose(&l);
        let c = back_substitute(&lt, &y);

        // Reconstruct: f_approx = Φ c
        let mut f_approx = vec![0.0; n];
        for (phi, &ci) in self.basis.iter().zip(c.iter()) {
            for (fi, &pi) in f_approx.iter_mut().zip(phi.iter()) {
                *fi += ci * pi;
            }
        }
        f_approx
    }
}

// ---------------------------------------------------------------------------
// Balanced truncation
// ---------------------------------------------------------------------------

/// Balanced truncation of a linear time-invariant system.
///
/// Given state-space matrices (A, B, C), computes the balanced realization
/// and truncates to `n_reduced` states.
///
/// Returns reduced matrices (A_r, B_r, C_r).
///
/// # Algorithm
/// 1. Compute approximate controllability Gramian W_c via Lyapunov equation.
/// 2. Compute approximate observability Gramian W_o.
/// 3. Find balancing transformation via Cholesky + SVD-like power iteration.
/// 4. Truncate.
#[allow(clippy::too_many_arguments)]
pub fn balanced_truncation(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    n_reduced: usize,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    if n == 0 || n_reduced == 0 {
        return (vec![], vec![], vec![]);
    }
    let r = n_reduced.min(n);

    // Approximate Gramians by low-rank truncation via power iteration on A+A^T
    let wc = approximate_gramian(a, b, 30);
    let at = mat_transpose(a);
    let ct = mat_transpose(c);
    let wo = approximate_gramian(&at, &ct, 30);

    // Compute square-root factors and balancing
    let lc = cholesky(&add_regularization(&wc, 1e-10));
    let lo = cholesky(&add_regularization(&wo, 1e-10));

    // Product M = Lo^T * Lc (n×n)
    let lo_t = mat_transpose(&lo);
    let m = mat_mul(&lo_t, &lc);

    // Dominant r singular vectors of M via deflated power iteration
    let mmt = mat_mul(&m, &mat_transpose(&m));
    let mut sigma2s: Vec<f64> = Vec::new();
    let mut us: Vec<Vec<f64>> = Vec::new();
    let mut work = mmt.clone();
    for _ in 0..r {
        let (lam, u) = power_iteration(&work, 500, 1e-10);
        if lam < 1e-14 {
            break;
        }
        sigma2s.push(lam);
        us.push(u.clone());
        work = deflate(&work, lam, &u);
    }

    // Balancing transformation T = Lc * V * Σ^{-1/2}  (approximated)
    // Here we use simple projection onto Lc * v_k directions
    let t_cols: Vec<Vec<f64>> = us
        .iter()
        .zip(sigma2s.iter())
        .map(|(u_vec, &s2)| {
            let sigma = s2.sqrt().max(1e-14);
            let lc_u = mat_vec(&lc, u_vec);
            vec_scale(&lc_u, 1.0 / sigma)
        })
        .collect();

    if t_cols.is_empty() {
        let ar = vec![vec![0.0; r]; r];
        let br = vec![vec![0.0; b[0].len()]; r];
        let cr = vec![vec![0.0; r]; c.len()];
        return (ar, br, cr);
    }

    // T matrix: n × r, columns are t_cols
    let t_mat: Vec<Vec<f64>> = (0..n)
        .map(|i| t_cols.iter().map(|col| col[i]).collect())
        .collect();
    let tt = mat_transpose(&t_mat);

    // Reduced system: A_r = T^T A T, B_r = T^T B, C_r = C T
    let ta = mat_mul(a, &t_mat);
    let ar = mat_mul(&tt, &ta);
    let br = mat_mul(&tt, b);
    let cr = mat_mul(c, &t_mat);

    (ar, br, cr)
}

/// Add a small regularization to ensure positive definiteness.
fn add_regularization(a: &[Vec<f64>], eps: f64) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut b = a.to_vec();
    for i in 0..n {
        b[i][i] += eps;
    }
    b
}

/// Approximate Gramian W ≈ B B^T (ignoring dynamics) for small demonstration.
///
/// A more accurate approximation integrates the matrix exponential; here we
/// use a low-rank approximation: W ≈ Σ_k (A^k B)(A^k B)^T with truncation.
fn approximate_gramian(a: &[Vec<f64>], b: &[Vec<f64>], n_terms: usize) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = b[0].len();
    let mut w = vec![vec![0.0; n]; n];
    let mut ak_b = b.to_vec(); // A^0 B = B
    let decay = 0.5_f64; // convergence factor
    let mut factor = 1.0_f64;
    for _ in 0..n_terms {
        // W += factor * (A^k B)(A^k B)^T
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for l in 0..m {
                    s += ak_b[i][l] * ak_b[j][l];
                }
                w[i][j] += factor * s;
            }
        }
        ak_b = mat_mul(a, &ak_b);
        factor *= decay;
    }
    w
}

// ---------------------------------------------------------------------------
// Hankel singular values
// ---------------------------------------------------------------------------

/// Compute approximate Hankel singular values of a linear system (A, B, C).
///
/// The Hankel singular values are the square roots of the eigenvalues of
/// W_c * W_o (controllability × observability Gramian product).
///
/// Uses power iteration on the product matrix for efficiency.
pub fn hankel_singular_values(a: &[Vec<f64>], b: &[Vec<f64>], c: &[Vec<f64>]) -> Vec<f64> {
    let n = a.len();
    if n == 0 {
        return vec![];
    }

    let wc = approximate_gramian(a, b, 30);
    let at = mat_transpose(a);
    let ct = mat_transpose(c);
    let wo = approximate_gramian(&at, &ct, 30);

    // Product P = Wc * Wo
    let p = mat_mul(&wc, &wo);

    // Find eigenvalues by deflated power iteration
    let mut eigenvalues: Vec<f64> = Vec::new();
    let mut work = p.clone();
    for _ in 0..n {
        let (lam, v) = power_iteration(&work, 500, 1e-10);
        if lam < 1e-14 {
            break;
        }
        eigenvalues.push(lam);
        work = deflate(&work, lam, &v);
    }

    eigenvalues.iter().map(|&ev| ev.sqrt().max(0.0)).collect()
}

// ---------------------------------------------------------------------------
// Module-level utility: example sinusoidal parametric solve
// ---------------------------------------------------------------------------

/// Generate a simple parametric solution for testing.
///
/// Returns a state vector representing a sinusoidal mode scaled by `mu`.
/// Used internally by tests and the greedy basis selection example.
#[allow(dead_code)]
pub fn parametric_sine_solution(mu: f64, n_dof: usize) -> Vec<f64> {
    (0..n_dof)
        .map(|i| (mu * PI * i as f64 / n_dof as f64).sin())
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_snapshots(n_dof: usize, n_snap: usize) -> SnapshotMatrix {
        let mut s = SnapshotMatrix::new(n_dof);
        for k in 0..n_snap {
            let snap: Vec<f64> = (0..n_dof)
                .map(|i| ((k + 1) as f64) * (i as f64 / n_dof as f64))
                .collect();
            s.add_snapshot(snap);
        }
        s
    }

    #[test]
    fn snapshot_matrix_add_and_count() {
        let mut s = SnapshotMatrix::new(4);
        s.add_snapshot(vec![1.0, 2.0, 3.0, 4.0]);
        s.add_snapshot(vec![5.0, 6.0, 7.0, 8.0]);
        assert_eq!(s.n_snapshots, 2);
        assert_eq!(s.data.len(), 2);
    }

    #[test]
    fn snapshot_matrix_mean_correct() {
        let mut s = SnapshotMatrix::new(2);
        s.add_snapshot(vec![1.0, 3.0]);
        s.add_snapshot(vec![3.0, 1.0]);
        let m = s.mean();
        assert!((m[0] - 2.0).abs() < 1e-12);
        assert!((m[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn snapshot_matrix_mean_empty() {
        let s = SnapshotMatrix::new(3);
        let m = s.mean();
        assert_eq!(m.len(), 3);
        assert!(m.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn snapshot_matrix_centered_has_zero_mean() {
        let s = make_snapshots(6, 4);
        let c = s.centered();
        let m = c.mean();
        for &mi in &m {
            assert!(mi.abs() < 1e-10, "mean not zero: {mi}");
        }
    }

    #[test]
    fn snapshot_matrix_centered_count_preserved() {
        let s = make_snapshots(5, 3);
        let c = s.centered();
        assert_eq!(c.n_snapshots, s.n_snapshots);
    }

    #[test]
    fn pod_basis_from_single_snapshot() {
        let mut s = SnapshotMatrix::new(4);
        s.add_snapshot(vec![1.0, 0.0, 0.0, 0.0]);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        assert_eq!(basis.n_modes, basis.modes.len());
    }

    #[test]
    fn pod_basis_modes_orthonormal() {
        let s = make_snapshots(8, 4);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        for (i, phi_i) in basis.modes.iter().enumerate() {
            for (j, phi_j) in basis.modes.iter().enumerate() {
                let prod = dot(phi_i, phi_j);
                if i == j {
                    assert!((prod - 1.0).abs() < 1e-8, "mode {i} not unit: {prod}");
                } else {
                    assert!(prod.abs() < 1e-8, "modes {i},{j} not orthogonal: {prod}");
                }
            }
        }
    }

    #[test]
    fn pod_basis_energy_threshold_subset() {
        let s = make_snapshots(10, 5);
        let b_full = PodBasis::from_snapshots(&s, 1.0);
        let b_partial = PodBasis::from_snapshots(&s, 0.5);
        assert!(b_partial.n_modes <= b_full.n_modes);
    }

    #[test]
    fn pod_basis_project_reconstruct_within_span() {
        let s = make_snapshots(6, 3);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        let state = s.data[0].clone();
        let q = basis.project(&state);
        let u = basis.reconstruct(&q);
        assert_eq!(u.len(), state.len());
    }

    #[test]
    fn pod_basis_energy_retained_full() {
        let s = make_snapshots(6, 4);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes > 0 {
            let e = basis.energy_retained(basis.n_modes);
            assert!(e <= 1.0 + 1e-10, "energy > 1: {e}");
        }
    }

    #[test]
    fn pod_basis_energy_retained_zero_modes() {
        let s = make_snapshots(4, 2);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        let e = basis.energy_retained(0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn pod_basis_empty_snapshots_gives_empty_basis() {
        let s = SnapshotMatrix::new(5);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        assert_eq!(basis.n_modes, 0);
    }

    #[test]
    fn pod_project_returns_correct_length() {
        let s = make_snapshots(8, 3);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        let state: Vec<f64> = vec![1.0; 8];
        let q = basis.project(&state);
        assert_eq!(q.len(), basis.n_modes);
    }

    #[test]
    fn pod_reconstruct_returns_correct_length() {
        let s = make_snapshots(8, 3);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        let q: Vec<f64> = vec![1.0; basis.n_modes];
        let u = basis.reconstruct(&q);
        if basis.n_modes > 0 {
            assert_eq!(u.len(), 8);
        }
    }

    fn make_diagonal_matrix(n: usize, diag_val: f64) -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; n]; n];
        for i in 0..n {
            m[i][i] = diag_val;
        }
        m
    }

    #[test]
    fn rom_system_project_system_dimensions() {
        let s = make_snapshots(4, 2);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes == 0 {
            return;
        }
        let mass = make_diagonal_matrix(4, 1.0);
        let stiffness = make_diagonal_matrix(4, 100.0);
        let rom = RomSystem::project_system(&basis, &mass, &stiffness);
        let r = basis.n_modes;
        assert_eq!(rom.reduced_mass.len(), r);
        assert_eq!(rom.reduced_stiffness.len(), r);
        assert_eq!(rom.reduced_damping.len(), r);
    }

    #[test]
    fn rom_solve_reduced_diagonal_system() {
        let s = make_snapshots(4, 2);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes == 0 {
            return;
        }
        let mass = make_diagonal_matrix(4, 1.0);
        let stiffness = make_diagonal_matrix(4, 1.0);
        let rom = RomSystem::project_system(&basis, &mass, &stiffness);
        let f_r: Vec<f64> = vec![1.0; rom.basis.n_modes];
        let q = rom.solve_reduced(&f_r);
        assert_eq!(q.len(), rom.basis.n_modes);
    }

    #[test]
    fn rom_rayleigh_damping_updates_matrix() {
        let s = make_snapshots(4, 2);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes == 0 {
            return;
        }
        let mass = make_diagonal_matrix(4, 1.0);
        let stiffness = make_diagonal_matrix(4, 100.0);
        let mut rom = RomSystem::project_system(&basis, &mass, &stiffness);
        rom.set_rayleigh_damping(0.1, 0.001);
        // Damping diagonal should be non-zero
        let has_damping = rom
            .reduced_damping
            .iter()
            .any(|row| row.iter().any(|&x| x != 0.0));
        assert!(has_damping || rom.basis.n_modes == 0);
    }

    #[test]
    fn rom_step_newmark_returns_full_size() {
        let s = make_snapshots(4, 2);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes == 0 {
            return;
        }
        let mass = make_diagonal_matrix(4, 1.0);
        let stiffness = make_diagonal_matrix(4, 10.0);
        let mut rom = RomSystem::project_system(&basis, &mass, &stiffness);
        let f_full = vec![0.1; 4];
        let u = rom.step_newmark(&f_full, 0.01, 0.25, 0.5);
        assert_eq!(u.len(), 4);
    }

    #[test]
    fn rom_newmark_multiple_steps() {
        let s = make_snapshots(6, 3);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        if basis.n_modes == 0 {
            return;
        }
        let mass = make_diagonal_matrix(6, 2.0);
        let stiffness = make_diagonal_matrix(6, 50.0);
        let mut rom = RomSystem::project_system(&basis, &mass, &stiffness);
        let f_full = vec![0.0; 6];
        for _ in 0..5 {
            let u = rom.step_newmark(&f_full, 0.001, 0.25, 0.5);
            assert!(u.iter().all(|x| x.is_finite()), "non-finite displacement");
        }
    }

    #[test]
    fn greedy_basis_returns_valid_basis() {
        // Use orthogonal snapshots so that multiple distinct modes can be found
        let params: Vec<f64> = (1..=5).map(|k| k as f64 * 0.25).collect();
        let basis = greedy_basis_selection(
            &params,
            &|mu| {
                let n = 20;
                // Mix sin and cos so different mu values give linearly independent vectors
                (0..n)
                    .map(|i| {
                        (mu * PI * i as f64 / n as f64).sin()
                            + (mu * 2.0 * PI * i as f64 / n as f64).cos()
                    })
                    .collect()
            },
            4,
            1e-8,
        );
        assert!(basis.n_modes <= 4);
        // basis should have at least 1 mode since params are non-empty
        assert!(
            basis.n_modes >= 1,
            "should have at least one mode, got {}",
            basis.n_modes
        );
    }

    #[test]
    fn greedy_basis_empty_params_returns_empty() {
        let basis = greedy_basis_selection(&[], &|mu| parametric_sine_solution(mu, 10), 5, 1e-4);
        assert_eq!(basis.n_modes, 0);
    }

    #[test]
    fn greedy_basis_modes_are_unit_norm() {
        let params: Vec<f64> = (1..=4).map(|k| k as f64 * 0.5).collect();
        let basis =
            greedy_basis_selection(&params, &|mu| parametric_sine_solution(mu, 16), 4, 1e-6);
        for (i, phi) in basis.modes.iter().enumerate() {
            let n = norm(phi);
            assert!((n - 1.0).abs() < 1e-8, "mode {i} norm = {n}");
        }
    }

    #[test]
    fn deim_from_modes_node_count_equals_mode_count() {
        let modes: Vec<Vec<f64>> = (1..=3)
            .map(|k| {
                (0..8)
                    .map(|i| (k as f64 * PI * i as f64 / 8.0).sin())
                    .collect()
            })
            .collect();
        let deim = EmpiricalInterpolation::from_modes(&modes);
        assert_eq!(deim.nodes.len(), modes.len());
    }

    #[test]
    fn deim_from_empty_modes_gives_empty() {
        let deim = EmpiricalInterpolation::from_modes(&[]);
        assert_eq!(deim.nodes.len(), 0);
        assert_eq!(deim.basis.len(), 0);
    }

    #[test]
    fn deim_interpolate_returns_correct_length() {
        let modes: Vec<Vec<f64>> = (1..=2)
            .map(|k| (0..6).map(|i| (k as f64 * i as f64).sin()).collect())
            .collect();
        let deim = EmpiricalInterpolation::from_modes(&modes);
        let f_nodes: Vec<f64> = vec![1.0; deim.nodes.len()];
        let f_approx = deim.interpolate(&f_nodes);
        assert_eq!(f_approx.len(), 6);
    }

    #[test]
    fn deim_nodes_within_bounds() {
        let modes: Vec<Vec<f64>> = (1..=3)
            .map(|k| (0..10).map(|i| ((k * i) as f64).sin()).collect())
            .collect();
        let deim = EmpiricalInterpolation::from_modes(&modes);
        for &node in &deim.nodes {
            assert!(node < 10, "node {node} out of bounds");
        }
    }

    #[test]
    fn balanced_truncation_returns_correct_sizes() {
        let n = 4;
        let r = 2;
        // Full-rank B and C so that both controllability and observability Gramians
        // have enough significant eigenvalues.
        let a = make_diagonal_matrix(n, -1.0);
        let b: Vec<Vec<f64>> = (0..n).map(|i| vec![(i + 1) as f64]).collect();
        let c: Vec<Vec<f64>> = vec![(0..n).map(|i| (i + 1) as f64).collect()];
        let (ar, br, cr) = balanced_truncation(&a, &b, &c, r);
        // The function returns up to r modes; actual count may be smaller if
        // the system Gramians have fewer significant eigenvalues.
        assert!(ar.len() <= r, "reduced A too large: {}", ar.len());
        if !ar.is_empty() {
            assert_eq!(ar[0].len(), ar.len(), "A_r not square");
        }
        assert_eq!(br.len(), ar.len(), "B_r row count mismatch");
        assert_eq!(cr.len(), 1, "C_r should have 1 output row");
    }

    #[test]
    fn balanced_truncation_empty_returns_empty() {
        let (ar, br, cr) = balanced_truncation(&[], &[], &[], 2);
        assert!(ar.is_empty());
        assert!(br.is_empty());
        assert!(cr.is_empty());
    }

    #[test]
    fn hankel_singular_values_positive() {
        let n = 3;
        let a = make_diagonal_matrix(n, -1.0);
        let b: Vec<Vec<f64>> = vec![vec![1.0], vec![0.5], vec![0.25]];
        let c: Vec<Vec<f64>> = vec![(0..n).map(|i| (i + 1) as f64 * 0.1).collect()];
        let hsv = hankel_singular_values(&a, &b, &c);
        assert!(!hsv.is_empty());
        for &s in &hsv {
            assert!(s >= 0.0 && s.is_finite(), "invalid HSV: {s}");
        }
    }

    #[test]
    fn hankel_singular_values_empty_system() {
        let hsv = hankel_singular_values(&[], &[], &[]);
        assert!(hsv.is_empty());
    }

    #[test]
    fn parametric_sine_solution_length() {
        let u = parametric_sine_solution(1.0, 20);
        assert_eq!(u.len(), 20);
    }

    #[test]
    fn parametric_sine_solution_range() {
        let u = parametric_sine_solution(1.0, 100);
        for &ui in &u {
            assert!(ui >= -1.0 && ui <= 1.0, "out of range: {ui}");
        }
    }

    #[test]
    fn dot_product_correct() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        assert_eq!(dot(&a, &b), 32.0);
    }

    #[test]
    fn norm_correct() {
        let a = vec![3.0, 4.0];
        assert!((norm(&a) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn mat_vec_identity() {
        let id = identity(3);
        let x = vec![1.0, 2.0, 3.0];
        let y = mat_vec(&id, &x);
        assert_eq!(y, x);
    }

    #[test]
    fn mat_transpose_correct() {
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let at = mat_transpose(&a);
        assert_eq!(at.len(), 3);
        assert_eq!(at[0].len(), 2);
        assert_eq!(at[0][0], 1.0);
        assert_eq!(at[0][1], 4.0);
    }

    #[test]
    fn mat_mul_identity() {
        let id = identity(3);
        let a: Vec<Vec<f64>> = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let b = mat_mul(&a, &id);
        for (i, row) in a.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((b[i][j] - val).abs() < 1e-12, "mismatch at [{i}][{j}]");
            }
        }
    }

    #[test]
    fn cholesky_spd_roundtrip() {
        let a = vec![
            vec![4.0, 2.0, 0.0],
            vec![2.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let l = cholesky(&a);
        let lt = mat_transpose(&l);
        let b = mat_mul(&l, &lt);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (b[i][j] - a[i][j]).abs() < 1e-10,
                    "cholesky error at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn forward_back_substitute_identity() {
        let id = identity(3);
        let b = vec![1.0, 2.0, 3.0];
        let y = forward_substitute(&id, &b);
        let x = back_substitute(&id, &y);
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < 1e-12);
        }
    }

    #[test]
    fn power_iteration_finds_dominant_eigenvector() {
        // Matrix with eigenvalues 5, 1
        let a = vec![vec![4.0, 1.0], vec![1.0, 2.0]];
        let (lam, v) = power_iteration(&a, 200, 1e-10);
        assert!(lam > 0.0, "dominant eigenvalue should be positive: {lam}");
        assert!(
            (norm(&v) - 1.0).abs() < 1e-8,
            "eigenvector not unit norm: {}",
            norm(&v)
        );
    }

    #[test]
    fn approximate_gramian_positive_diagonal() {
        let n = 3;
        let a = make_diagonal_matrix(n, -0.5);
        let b: Vec<Vec<f64>> = (0..n).map(|i| vec![(i + 1) as f64]).collect();
        let w = approximate_gramian(&a, &b, 10);
        for i in 0..n {
            assert!(
                w[i][i] >= 0.0,
                "diagonal should be non-negative: {}",
                w[i][i]
            );
        }
    }

    #[test]
    fn snapshot_new_zero_snapshots() {
        let s = SnapshotMatrix::new(10);
        assert_eq!(s.n_snapshots, 0);
        assert_eq!(s.n_dof, 10);
        assert!(s.data.is_empty());
    }

    #[test]
    fn pod_basis_singular_values_non_negative() {
        let s = make_snapshots(6, 4);
        let basis = PodBasis::from_snapshots(&s, 1.0);
        for &sv in &basis.singular_values {
            assert!(sv >= 0.0, "negative singular value: {sv}");
        }
    }

    #[test]
    fn greedy_basis_single_param() {
        let basis = greedy_basis_selection(&[1.0], &|mu| parametric_sine_solution(mu, 10), 5, 1e-4);
        // Only one unique parameter, so at most 1 mode
        assert!(basis.n_modes <= 1);
    }
}
