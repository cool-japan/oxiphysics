#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Standalone corotational FEM functions: polar decomposition, strain/stress
//! computations, stiffness matrix assembly, solvers, and time-integration helpers.

#![allow(clippy::needless_range_loop)]

use super::math_helpers::{
    det3x3, inv3x3, inv3x3_transpose, isotropic_d_matrix, mul3x3, transpose3x3,
};

// ---------------------------------------------------------------------------
// Polar decomposition utilities (SVD-free, iterative)
// ---------------------------------------------------------------------------

/// Compute the polar decomposition F = R * S using the iterative method of
/// Higham (repeated averaging).  Returns the rotation matrix R such that
/// R^T * R = I and det(R) = +1.
///
/// Convergence is guaranteed for non-singular F.
#[allow(dead_code)]
pub fn polar_decompose_r(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = f;
    for _ in 0..40 {
        let r_inv_t = inv3x3_transpose(r);
        let mut r_next = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r_next[i][j] = 0.5 * (r[i][j] + r_inv_t[i][j]);
            }
        }
        let mut diff = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                diff += (r_next[i][j] - r[i][j]).abs();
            }
        }
        r = r_next;
        if diff < 1e-12 {
            break;
        }
    }
    r
}

/// Compute the symmetric stretch matrix S = R^T * F given F and R.
#[allow(dead_code)]
pub fn polar_decompose_s(r: [[f64; 3]; 3], f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let rt = transpose3x3(r);
    mul3x3(rt, f)
}

// ---------------------------------------------------------------------------
// Corotational strain measure
// ---------------------------------------------------------------------------

/// Compute the corotational (linear) strain tensor:
///
///   eps = sym(R^T * F) - I  =  sym(S) - I
///
/// where S = R^T * F is the symmetric stretch matrix (symmetric by construction
/// after polar decomposition).
#[allow(dead_code)]
pub fn corotational_strain(r: [[f64; 3]; 3], f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let s = polar_decompose_s(r, f);
    // sym(S) = 0.5*(S + S^T); since S should already be symmetric after polar
    // decompose, we symmetrize for numerical safety.
    let mut eps = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            eps[i][j] = 0.5 * (s[i][j] + s[j][i]);
        }
    }
    // Subtract identity
    for i in 0..3 {
        eps[i][i] -= 1.0;
    }
    eps
}

// ---------------------------------------------------------------------------
// Cauchy stress from Lame parameters (linear isotropic)
// ---------------------------------------------------------------------------

/// Compute the Cauchy stress tensor from the corotational strain tensor using
/// linear isotropic elasticity:
///
///   sigma = lambda * tr(eps) * I + 2*mu * eps
///
/// where lambda and mu are the Lame parameters.
#[allow(dead_code)]
pub fn corot_cauchy_stress(eps: [[f64; 3]; 3], lambda: f64, mu: f64) -> [[f64; 3]; 3] {
    let trace = eps[0][0] + eps[1][1] + eps[2][2];
    let mut sigma = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            sigma[i][j] = 2.0 * mu * eps[i][j];
        }
        sigma[i][i] += lambda * trace;
    }
    sigma
}

// ---------------------------------------------------------------------------
// Element stiffness matrix (12x12) for a corotational tet
// ---------------------------------------------------------------------------

/// Compute the 12x12 material stiffness matrix for a linear-elastic
/// tetrahedral element with Lame parameters lambda and mu and volume `vol`.
///
/// The element has 4 nodes x 3 DOF = 12 DOFs.  The returned matrix is stored
/// in row-major order as a flat array of 144 entries.
///
/// Uses the classical linear tetrahedral element B-matrix.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn tet_stiffness_matrix(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    lambda: f64,
    mu: f64,
) -> [f64; 144] {
    // Edge matrix columns: p1-p0, p2-p0, p3-p0
    let dm = [
        [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
        [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
        [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
    ];
    let det = det3x3(dm);
    if det.abs() < 1e-14 {
        return [0.0_f64; 144];
    }
    let vol = det.abs() / 6.0;
    let dm_inv = inv3x3(dm);

    // Shape-function gradients (4x3): each row = [dN/dx, dN/dy, dN/dz]
    // dN1..dN3 come from dm_inv columns; dN0 = -(dN1+dN2+dN3)
    let mut grad = [[0.0_f64; 3]; 4];
    for d in 0..3 {
        grad[1][d] = dm_inv[d][0];
        grad[2][d] = dm_inv[d][1];
        grad[3][d] = dm_inv[d][2];
        grad[0][d] -= dm_inv[d][0] + dm_inv[d][1] + dm_inv[d][2];
    }

    // 6x12 B matrix (Voigt engineering notation: [exx,eyy,ezz,2eyz,2exz,2exy])
    let mut b = [[0.0_f64; 12]; 6];
    for a in 0..4 {
        let col = a * 3;
        let gx = grad[a][0];
        let gy = grad[a][1];
        let gz = grad[a][2];
        b[0][col] = gx;
        b[1][col + 1] = gy;
        b[2][col + 2] = gz;
        b[3][col + 1] = gz;
        b[3][col + 2] = gy;
        b[4][col] = gz;
        b[4][col + 2] = gx;
        b[5][col] = gy;
        b[5][col + 1] = gx;
    }

    // 6x6 constitutive matrix D (isotropic linear elasticity)
    let d = isotropic_d_matrix(lambda, mu);

    // K = vol * B^T * D * B  (12x12)
    // First compute DB = D * B (6x12)
    let mut db = [[0.0_f64; 12]; 6];
    for i in 0..6 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..6 {
                s += d[i][k] * b[k][j];
            }
            db[i][j] = s;
        }
    }
    // Then K = vol * B^T * DB (12x12)
    let mut k_mat = [0.0_f64; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for r in 0..6 {
                s += b[r][i] * db[r][j];
            }
            k_mat[i * 12 + j] = vol * s;
        }
    }
    k_mat
}

// ---------------------------------------------------------------------------
// Corotational internal forces via stiffness matrix
// ---------------------------------------------------------------------------

/// Compute the corotational internal forces for a tetrahedral element.
///
/// Uses the formulation:
///   f_int = R * K_e * (R^T * u_deformed - u_rest)
///
/// where K_e is the linear stiffness, R is the element rotation, and the
/// displacement is measured in the unrotated reference frame.
///
/// Returns 4 nodal force vectors.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn corot_internal_forces(
    rest: [[f64; 3]; 4],
    deformed: [[f64; 3]; 4],
    lambda: f64,
    mu: f64,
) -> [[f64; 3]; 4] {
    let p0 = rest[0];
    let p1 = rest[1];
    let p2 = rest[2];
    let p3 = rest[3];

    // Deformation gradient F = Ds * Dm^{-1}
    let dm = [
        [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
        [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
        [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
    ];
    let det_dm = det3x3(dm);
    if det_dm.abs() < 1e-14 {
        return [[0.0; 3]; 4];
    }
    let dm_inv = inv3x3(dm);

    let q0 = deformed[0];
    let q1 = deformed[1];
    let q2 = deformed[2];
    let q3 = deformed[3];
    let ds = [
        [q1[0] - q0[0], q2[0] - q0[0], q3[0] - q0[0]],
        [q1[1] - q0[1], q2[1] - q0[1], q3[1] - q0[1]],
        [q1[2] - q0[2], q2[2] - q0[2], q3[2] - q0[2]],
    ];
    let f_grad = mul3x3(ds, dm_inv);

    // Polar decompose F -> R
    let r = polar_decompose_r(f_grad);
    let rt = transpose3x3(r);

    // Map deformed positions into unrotated frame: x_bar = R^T * x
    let mut x_bar = [[0.0_f64; 3]; 4];
    for a in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                x_bar[a][i] += rt[i][j] * deformed[a][j];
            }
        }
    }

    // Displacement in unrotated frame relative to rest
    let mut u_bar = [0.0_f64; 12];
    for a in 0..4 {
        for d in 0..3 {
            u_bar[a * 3 + d] = x_bar[a][d] - rest[a][d];
        }
    }

    // K_e (linear stiffness in rest frame)
    let ke = tet_stiffness_matrix(p0, p1, p2, p3, lambda, mu);

    // f_bar = K_e * u_bar  (forces in unrotated frame)
    let mut f_bar = [0.0_f64; 12];
    for i in 0..12 {
        for j in 0..12 {
            f_bar[i] += ke[i * 12 + j] * u_bar[j];
        }
    }

    // Rotate forces back to world frame: f = R * f_bar
    let mut forces = [[0.0_f64; 3]; 4];
    for a in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                forces[a][i] += r[i][j] * f_bar[a * 3 + j];
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// Newmark-beta / HHT-alpha time integration helper
// ---------------------------------------------------------------------------

/// Parameters for HHT-alpha (Hilber-Hughes-Taylor) time integration.
///
/// Special case alpha = 0 gives classic Newmark-beta (beta = 0.25, gamma = 0.5).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HhtAlphaParams {
    /// Algorithmic dissipation parameter in \[-1/3, 0\].
    pub alpha: f64,
    /// Newmark beta parameter.
    pub beta: f64,
    /// Newmark gamma parameter.
    pub gamma: f64,
    /// Time-step size (s).
    pub dt: f64,
}

impl HhtAlphaParams {
    /// Construct HHT-alpha parameters from alpha alone; beta and gamma are derived for
    /// unconditional stability:
    ///
    ///   beta = (1 - alpha)^2 / 4,  gamma = (1 - 2*alpha) / 2
    #[allow(dead_code)]
    pub fn from_alpha(alpha: f64, dt: f64) -> Self {
        let beta = (1.0 - alpha) * (1.0 - alpha) / 4.0;
        let gamma = (1.0 - 2.0 * alpha) / 2.0;
        Self {
            alpha,
            beta,
            gamma,
            dt,
        }
    }

    /// Predict displacement and velocity at t_{n+1} from known quantities.
    ///
    /// Returned as `(u_pred, v_pred)`.
    #[allow(dead_code)]
    pub fn predict(
        &self,
        u: &[[f64; 3]],
        v: &[[f64; 3]],
        a: &[[f64; 3]],
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = u.len();
        let dt = self.dt;
        let mut u_pred = vec![[0.0_f64; 3]; n];
        let mut v_pred = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for d in 0..3 {
                u_pred[i][d] = u[i][d] + dt * v[i][d] + dt * dt * (0.5 - self.beta) * a[i][d];
                v_pred[i][d] = v[i][d] + dt * (1.0 - self.gamma) * a[i][d];
            }
        }
        (u_pred, v_pred)
    }

    /// Correct displacement and velocity given new acceleration a_{n+1}.
    #[allow(dead_code)]
    pub fn correct(
        &self,
        u_pred: &mut Vec<[f64; 3]>,
        v_pred: &mut Vec<[f64; 3]>,
        a_new: &[[f64; 3]],
    ) {
        let dt = self.dt;
        for i in 0..u_pred.len() {
            for d in 0..3 {
                u_pred[i][d] += self.beta * dt * dt * a_new[i][d];
                v_pred[i][d] += self.gamma * dt * a_new[i][d];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Explicit central-difference (leap-frog) integrator
// ---------------------------------------------------------------------------

/// One-step explicit central-difference (leap-frog) integrator.
///
/// Updates positions and half-step velocities in place.
///
/// - `x`:    current positions (updated to t_{n+1}).
/// - `v_half`: half-step velocities (updated from t_{n-1/2} to t_{n+1/2}).
/// - `a`:    accelerations at t_n.
/// - `dt`:   time-step.
/// - `pinned`: nodes to keep fixed.
#[allow(dead_code)]
pub fn central_difference_step(
    x: &mut Vec<[f64; 3]>,
    v_half: &mut Vec<[f64; 3]>,
    a: &[[f64; 3]],
    dt: f64,
    pinned: &[bool],
) {
    for i in 0..x.len() {
        if pinned.get(i).copied().unwrap_or(false) {
            continue;
        }
        for d in 0..3 {
            v_half[i][d] += a[i][d] * dt;
            x[i][d] += v_half[i][d] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Global stiffness matrix assembly
// ---------------------------------------------------------------------------

/// Assemble the global 12-DOF-per-node stiffness matrix into a flat
/// `(3*n_nodes) x (3*n_nodes)` dense matrix (row-major).
///
/// Only non-zero blocks for the 4 nodes of each element are added.
#[allow(dead_code)]
pub fn assemble_global_stiffness(
    positions: &[[f64; 3]],
    elements: &[([usize; 4], f64, f64)], // (indices, lambda, mu)
) -> Vec<f64> {
    let n = positions.len();
    let dof = 3 * n;
    let mut k_global = vec![0.0_f64; dof * dof];
    for (indices, lambda, mu) in elements {
        let [i0, i1, i2, i3] = *indices;
        let p = [positions[i0], positions[i1], positions[i2], positions[i3]];
        let ke = tet_stiffness_matrix(p[0], p[1], p[2], p[3], *lambda, *mu);
        let node_ids = [i0, i1, i2, i3];
        for a in 0..4 {
            for da in 0..3 {
                let row = node_ids[a] * 3 + da;
                for b in 0..4 {
                    for db in 0..3 {
                        let col = node_ids[b] * 3 + db;
                        k_global[row * dof + col] += ke[(a * 3 + da) * 12 + b * 3 + db];
                    }
                }
            }
        }
    }
    k_global
}

// ---------------------------------------------------------------------------
// Conjugate gradient solver (for K*u = f)
// ---------------------------------------------------------------------------

/// Solve the symmetric positive-definite system `A * x = b` using the
/// conjugate-gradient method.
///
/// `a_flat`: row-major dense `n x n` matrix.
/// `b`:      right-hand side vector of length `n`.
/// `tol`:    convergence tolerance on the residual norm.
/// `max_iter`: maximum iterations.
///
/// Returns the solution vector `x`.
#[allow(dead_code)]
pub fn conjugate_gradient(a_flat: &[f64], b: &[f64], tol: f64, max_iter: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rsold: f64 = r.iter().map(|v| v * v).sum();
    for _ in 0..max_iter {
        if rsold.sqrt() < tol {
            break;
        }
        // Ap = A * p
        let mut ap = vec![0.0_f64; n];
        for i in 0..n {
            for j in 0..n {
                ap[i] += a_flat[i * n + j] * p[j];
            }
        }
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, ai)| pi * ai).sum();
        if pap.abs() < 1e-30 {
            break;
        }
        let alpha = rsold / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
        }
        for i in 0..n {
            r[i] -= alpha * ap[i];
        }
        let rsnew: f64 = r.iter().map(|v| v * v).sum();
        let beta = rsnew / rsold;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    x
}

// ---------------------------------------------------------------------------
// Element volume calculation
// ---------------------------------------------------------------------------

/// Compute the signed volume of a tetrahedron given its four vertex positions.
#[allow(dead_code)]
pub fn tet_signed_volume(p: [[f64; 3]; 4]) -> f64 {
    let dm = [
        [p[1][0] - p[0][0], p[2][0] - p[0][0], p[3][0] - p[0][0]],
        [p[1][1] - p[0][1], p[2][1] - p[0][1], p[3][1] - p[0][1]],
        [p[1][2] - p[0][2], p[2][2] - p[0][2], p[3][2] - p[0][2]],
    ];
    det3x3(dm) / 6.0
}

/// Compute the current volume of all elements and return min/max/total.
#[allow(dead_code)]
pub fn element_volume_stats(
    positions: &[[f64; 3]],
    element_indices: &[[usize; 4]],
) -> (f64, f64, f64) {
    let mut min_vol = f64::INFINITY;
    let mut max_vol = f64::NEG_INFINITY;
    let mut total = 0.0;
    for idx in element_indices {
        let p = [
            positions[idx[0]],
            positions[idx[1]],
            positions[idx[2]],
            positions[idx[3]],
        ];
        let v = tet_signed_volume(p).abs();
        min_vol = min_vol.min(v);
        max_vol = max_vol.max(v);
        total += v;
    }
    (min_vol, max_vol, total)
}

// ---------------------------------------------------------------------------
// Gravity force vector
// ---------------------------------------------------------------------------

/// Compute gravitational forces for all nodes given per-node masses.
///
/// `g_vec` is the gravitational acceleration vector (e.g. `[0.0, -9.81, 0.0]`).
#[allow(dead_code)]
pub fn gravity_forces(masses: &[f64], g_vec: [f64; 3]) -> Vec<[f64; 3]> {
    masses
        .iter()
        .map(|&m| [m * g_vec[0], m * g_vec[1], m * g_vec[2]])
        .collect()
}

// ---------------------------------------------------------------------------
// Lumped-mass acceleration from forces
// ---------------------------------------------------------------------------

/// Given lumped masses and applied forces (including internal + external),
/// compute nodal accelerations `a = f / m` for free nodes.
#[allow(dead_code)]
pub fn lumped_mass_accelerations(
    masses: &[f64],
    forces: &[[f64; 3]],
    pinned: &[bool],
) -> Vec<[f64; 3]> {
    let n = masses.len();
    let mut acc = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        if pinned.get(i).copied().unwrap_or(false) {
            continue;
        }
        let m = masses[i].max(1e-30);
        for d in 0..3 {
            acc[i][d] = forces[i][d] / m;
        }
    }
    acc
}
