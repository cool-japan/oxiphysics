//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{CsrMatrix, ImplicitParticle, ImplicitSpring, LineSearch, PcgSolver};

#[inline]
pub(super) fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}
#[cfg(test)]
#[inline]
pub(super) fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
#[inline]
pub(super) fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n > 1e-15 {
        v3_scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}
/// Flat vector dot product.
pub(super) fn vec_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Flat vector norm.
pub(super) fn vec_norm(a: &[f64]) -> f64 {
    vec_dot(a, a).sqrt()
}
pub(super) fn mat3_det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
pub(super) fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}
pub(super) fn mat3_inv_scale(m: [[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}
pub(super) fn mat3_mul_v3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// Solve `H x = b` using matrix-free CG.
///
/// `H` is provided as a function `apply: &[[f64;3\]] -> &[f64] -> Vec`f64`.
///
/// Returns the solution vector.
pub fn matrix_free_cg<F>(
    positions: &[[f64; 3]],
    b: &[f64],
    apply: F,
    max_iter: usize,
    tol: f64,
) -> Vec<f64>
where
    F: Fn(&[[f64; 3]], &[f64]) -> Vec<f64>,
{
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr = vec_dot(&r, &r);
    for _ in 0..max_iter {
        let ap = apply(positions, &p);
        let pap = vec_dot(&p, &ap);
        if pap.abs() < 1e-30 {
            break;
        }
        let alpha = rr / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
        }
        for i in 0..n {
            r[i] -= alpha * ap[i];
        }
        let rr_new = vec_dot(&r, &r);
        if rr_new.sqrt() < tol {
            break;
        }
        let beta = rr_new / rr.max(1e-30);
        rr = rr_new;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
    }
    x
}
/// Compute total elastic potential energy for a spring-mass system.
///
/// `E = sum_springs 0.5 * k * (|x_i - x_j| - L)^2`
pub fn spring_elastic_energy(positions: &[[f64; 3]], springs: &[ImplicitSpring]) -> f64 {
    springs
        .iter()
        .map(|s| {
            let diff = v3_sub(positions[s.i], positions[s.j]);
            let len = v3_norm(diff);
            let stretch = len - s.rest_length;
            0.5 * s.stiffness * stretch * stretch
        })
        .sum()
}
/// Zero out gradient entries corresponding to static (pinned) particles.
pub fn apply_static_mask(grad: &mut [f64], is_static: &[bool]) {
    for (i, &fixed) in is_static.iter().enumerate() {
        if fixed {
            grad[3 * i] = 0.0;
            grad[3 * i + 1] = 0.0;
            grad[3 * i + 2] = 0.0;
        }
    }
}
/// Particles lighter than this are treated as static (pinned).
const STATIC_MASS_THRESHOLD: f64 = 1e-15;
/// Springs shorter than this have no defined direction and exert no force.
const MIN_SPRING_LENGTH: f64 = 1e-15;
/// Add the 3-vector `v` to the degrees of freedom of particle `i` in a flat `3n` vector.
#[inline]
fn add_to_particle(flat: &mut [f64], i: usize, v: [f64; 3]) {
    flat[3 * i] += v[0];
    flat[3 * i + 1] += v[1];
    flat[3 * i + 2] += v[2];
}
/// Coordinate-format (triplet) accumulator for block-assembled sparse matrices.
struct CooBuilder {
    rows: Vec<usize>,
    cols: Vec<usize>,
    vals: Vec<f64>,
}
impl CooBuilder {
    fn with_capacity(capacity: usize) -> Self {
        CooBuilder {
            rows: Vec::with_capacity(capacity),
            cols: Vec::with_capacity(capacity),
            vals: Vec::with_capacity(capacity),
        }
    }
    fn push(&mut self, row: usize, col: usize, value: f64) {
        self.rows.push(row);
        self.cols.push(col);
        self.vals.push(value);
    }
    /// Add `sign * block` with its top-left corner at the degrees of freedom of
    /// the particle pair `(row_particle, col_particle)`.
    fn push_block(
        &mut self,
        row_particle: usize,
        col_particle: usize,
        block: &[[f64; 3]; 3],
        sign: f64,
    ) {
        for (r, block_row) in block.iter().enumerate() {
            for (c, &value) in block_row.iter().enumerate() {
                self.push(3 * row_particle + r, 3 * col_particle + c, sign * value);
            }
        }
    }
    fn into_csr(self, n_dof: usize) -> CsrMatrix {
        CsrMatrix::from_coo(n_dof, n_dof, &self.rows, &self.cols, &self.vals)
    }
}
/// Frozen data of the implicit spring-damper model used by backward Euler.
///
/// The damper of spring `s` acts on the relative velocity along the spring
/// direction.  Backward Euler evaluates that velocity at the end of the step,
/// `v = (x - x^n) / dt`; freezing the spring direction `d^n` at the start of
/// the step keeps the damping force linear in the unknown positions, so it is
/// the exact gradient of the quadratic dissipation potential
/// `c / (2 dt) * (d^n . ((x_i - x_i^n) - (x_j - x_j^n)))^2`.  The step stays a
/// pure minimisation problem (energy line search stays valid) and the damping
/// is unconditionally stable.
struct DampingState {
    /// Particle positions at the start of the step.
    reference: Vec<[f64; 3]>,
    /// Per spring: unit direction at the start of the step, or `None` when the
    /// spring has no (positive) damping or zero length.
    directions: Vec<Option<[f64; 3]>>,
}
impl DampingState {
    fn new(reference: Vec<[f64; 3]>, springs: &[ImplicitSpring]) -> Self {
        let directions = springs
            .iter()
            .map(|s| -> Option<[f64; 3]> {
                let damped = s.damping > 0.0;
                if !damped {
                    return None;
                }
                let diff = v3_sub(*reference.get(s.i)?, *reference.get(s.j)?);
                let len = v3_norm(diff);
                if len < MIN_SPRING_LENGTH {
                    None
                } else {
                    Some(v3_scale(diff, 1.0 / len))
                }
            })
            .collect();
        DampingState {
            reference,
            directions,
        }
    }
    /// Relative displacement of the spring end points since the start of the
    /// step, projected on the frozen spring direction.
    fn relative_displacement(
        &self,
        spring_index: usize,
        spring: &ImplicitSpring,
        x: &[[f64; 3]],
    ) -> Option<(f64, [f64; 3])> {
        let d = (*self.directions.get(spring_index)?)?;
        let di = v3_sub(x[spring.i], self.reference[spring.i]);
        let dj = v3_sub(x[spring.j], self.reference[spring.j]);
        Some((v3_dot(d, v3_sub(di, dj)), d))
    }
}
/// The 3x3 stiffness block of one spring, `K = k [ d d^T + m (I - d d^T) ]`
/// with `m = max(1 - L0 / l, 0)`.
///
/// This is the Hessian of `0.5 k (l - L0)^2` with respect to `x_i`, with the
/// (negative) transverse stiffness of a compressed spring clamped to zero so
/// that the block stays positive semi-definite.
fn spring_stiffness_block(spring: &ImplicitSpring, d: [f64; 3], len: f64) -> [[f64; 3]; 3] {
    let transverse = (1.0 - spring.rest_length / len).max(0.0);
    let mut block = [[0.0f64; 3]; 3];
    for (r, block_row) in block.iter_mut().enumerate() {
        for (c, entry) in block_row.iter_mut().enumerate() {
            let identity = if r == c { transverse } else { 0.0 };
            *entry = spring.stiffness * ((1.0 - transverse) * d[r] * d[c] + identity);
        }
    }
    block
}
/// Assemble `H = M / dt^2 + K(x) [+ D / dt]` as a symmetric CSR matrix.
///
/// Static particles (`mass < 1e-15`) receive an identity row and column and
/// no coupling terms; springs with out-of-range or equal end points are
/// ignored.  With `damping`, every damped spring adds `(c / dt) d^n d^n^T` in
/// the same `+/-` block pattern as its stiffness.
fn assemble_hessian(
    positions: &[[f64; 3]],
    masses: &[f64],
    springs: &[ImplicitSpring],
    dt: f64,
    damping: Option<&DampingState>,
) -> CsrMatrix {
    let n = positions.len();
    let mut coo = CooBuilder::with_capacity(3 * n + 36 * springs.len());
    for (i, &mass) in masses.iter().enumerate().take(n) {
        let diagonal = if mass < STATIC_MASS_THRESHOLD {
            1.0
        } else {
            mass / (dt * dt)
        };
        for axis in 0..3 {
            coo.push(3 * i + axis, 3 * i + axis, diagonal);
        }
    }
    let is_free = |i: usize| masses.get(i).is_some_and(|&m| m >= STATIC_MASS_THRESHOLD);
    for (index, s) in springs.iter().enumerate() {
        let (Some(&xi), Some(&xj)) = (positions.get(s.i), positions.get(s.j)) else {
            continue;
        };
        let (free_i, free_j) = (is_free(s.i), is_free(s.j));
        if s.i == s.j || (!free_i && !free_j) {
            continue;
        }
        let diff = v3_sub(xi, xj);
        let len = v3_norm(diff);
        let mut block = [[0.0f64; 3]; 3];
        let mut contributes = false;
        if len >= MIN_SPRING_LENGTH {
            block = spring_stiffness_block(s, v3_scale(diff, 1.0 / len), len);
            contributes = true;
        }
        if let Some(d) = damping.and_then(|state| state.directions.get(index).copied().flatten()) {
            let scale = s.damping / dt;
            for (r, block_row) in block.iter_mut().enumerate() {
                for (c, entry) in block_row.iter_mut().enumerate() {
                    *entry += scale * d[r] * d[c];
                }
            }
            contributes = true;
        }
        if !contributes {
            continue;
        }
        if free_i {
            coo.push_block(s.i, s.i, &block, 1.0);
        }
        if free_j {
            coo.push_block(s.j, s.j, &block, 1.0);
        }
        if free_i && free_j {
            coo.push_block(s.i, s.j, &block, -1.0);
            coo.push_block(s.j, s.i, &block, -1.0);
        }
    }
    coo.into_csr(3 * n)
}
/// Assemble the sparse system matrix of the implicit spring-mass problem.
///
/// Returns `H = M / dt^2 + K(x)` for the particle positions currently stored
/// in `particles`, where `K` is the Hessian of the spring potential.  A spring
/// `(i, j)` of stiffness `k`, rest length `L0`, current length `l` and unit
/// direction `d = (x_i - x_j) / l` contributes the 3x3 block
///
/// ```text
/// K_s = k [ d d^T + max(1 - L0 / l, 0) (I - d d^T) ]
/// ```
///
/// as `+K_s` at `(i, i)` and `(j, j)` and `-K_s` at `(i, j)` and `(j, i)`.
/// The transverse term is clamped at zero so a compressed spring keeps only
/// its (always positive) axial stiffness: `H` is then symmetric positive
/// definite and can be solved with [`PcgSolver`].
///
/// Static particles (mass below `1e-15`) get an identity row and column (a 1 on
/// the diagonal, no coupling), which keeps `H` symmetric and makes their
/// Newton update zero for a zero right-hand side.  Springs that reference a
/// particle index out of range, join a particle to itself, or currently have
/// zero length are ignored.
pub fn assemble_spring_hessian(
    particles: &[ImplicitParticle],
    springs: &[ImplicitSpring],
    dt: f64,
) -> CsrMatrix {
    let positions: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
    let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    assemble_hessian(&positions, &masses, springs, dt, None)
}
/// Variational description of one implicit time step of a spring-mass system.
///
/// The step is the minimiser of the incremental potential
///
/// ```text
/// Phi(x) = sum_i  m_i / (2 dt^2) |x_i - y_i|^2          (inertia, target y)
///        + sum_s  k_s / 2 (|x_i - x_j| - L0_s)^2          (springs)
///        - sum_i  f_i . (x_i - x0_i)                      (external force)
///        + sum_s  c_s / (2 dt) (d_s^n . dx_ij)^2          (implicit damping, optional)
/// ```
///
/// Energy, gradient and Hessian are mutually consistent (the Hessian is the
/// positive semi-definite projection of the exact one), so a Newton direction
/// from [`ImplicitStepProblem::hessian`] is always a descent direction for
/// [`ImplicitStepProblem::energy`].  The external work is measured from the
/// reference positions `x0` only to keep the magnitudes (and thus the round-off
/// in the line search) small.
pub(super) struct ImplicitStepProblem<'a> {
    springs: &'a [ImplicitSpring],
    masses: Vec<f64>,
    inertia_target: Vec<[f64; 3]>,
    external_force: Vec<[f64; 3]>,
    origin: Vec<[f64; 3]>,
    dt: f64,
    damping: Option<DampingState>,
}
impl<'a> ImplicitStepProblem<'a> {
    /// Build the problem, or return `None` if the inputs are inconsistent: all
    /// per-particle arrays must have `masses.len()` entries, `dt` must be
    /// positive and finite, and every spring must reference valid particles.
    pub(super) fn new(
        masses: Vec<f64>,
        springs: &'a [ImplicitSpring],
        inertia_target: Vec<[f64; 3]>,
        external_force: Vec<[f64; 3]>,
        origin: Vec<[f64; 3]>,
        dt: f64,
        damping_reference: Option<Vec<[f64; 3]>>,
    ) -> Option<Self> {
        let n = masses.len();
        let sizes_match = inertia_target.len() == n
            && external_force.len() == n
            && origin.len() == n
            && damping_reference.as_ref().is_none_or(|r| r.len() == n);
        let dt_valid = dt.is_finite() && dt > 0.0;
        let springs_valid = springs.iter().all(|s| s.i < n && s.j < n);
        if !(sizes_match && dt_valid && springs_valid) {
            return None;
        }
        let damping = damping_reference.map(|reference| DampingState::new(reference, springs));
        Some(ImplicitStepProblem {
            springs,
            masses,
            inertia_target,
            external_force,
            origin,
            dt,
            damping,
        })
    }
    /// `true` if particle `i` is pinned.
    pub(super) fn is_static(&self, i: usize) -> bool {
        self.masses[i] < STATIC_MASS_THRESHOLD
    }
    /// Value of the incremental potential at positions `x`.
    pub(super) fn energy(&self, x: &[[f64; 3]]) -> f64 {
        let dt = self.dt;
        let mut e = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            if self.is_static(i) {
                continue;
            }
            let diff = v3_sub(xi, self.inertia_target[i]);
            e += 0.5 * self.masses[i] / (dt * dt) * v3_dot(diff, diff);
            e -= v3_dot(self.external_force[i], v3_sub(xi, self.origin[i]));
        }
        for (index, s) in self.springs.iter().enumerate() {
            if s.i == s.j || (self.is_static(s.i) && self.is_static(s.j)) {
                continue;
            }
            let stretch = v3_norm(v3_sub(x[s.i], x[s.j])) - s.rest_length;
            e += 0.5 * s.stiffness * stretch * stretch;
            if let Some((rel, _)) = self
                .damping
                .as_ref()
                .and_then(|state| state.relative_displacement(index, s, x))
            {
                e += s.damping / (2.0 * dt) * rel * rel;
            }
        }
        e
    }
    /// Gradient of the incremental potential as a flat `3n` vector (zero on
    /// the degrees of freedom of static particles).
    pub(super) fn gradient(&self, x: &[[f64; 3]]) -> Vec<f64> {
        let dt = self.dt;
        let mut grad = vec![0.0f64; 3 * x.len()];
        for (i, &xi) in x.iter().enumerate() {
            if self.is_static(i) {
                continue;
            }
            let inertial = v3_scale(
                v3_sub(xi, self.inertia_target[i]),
                self.masses[i] / (dt * dt),
            );
            add_to_particle(&mut grad, i, v3_sub(inertial, self.external_force[i]));
        }
        for (index, s) in self.springs.iter().enumerate() {
            let (free_i, free_j) = (!self.is_static(s.i), !self.is_static(s.j));
            if s.i == s.j || (!free_i && !free_j) {
                continue;
            }
            let diff = v3_sub(x[s.i], x[s.j]);
            let len = v3_norm(diff);
            let mut at_i = [0.0f64; 3];
            if len >= MIN_SPRING_LENGTH {
                at_i = v3_scale(diff, s.stiffness * (len - s.rest_length) / len);
            }
            if let Some((rel, d)) = self
                .damping
                .as_ref()
                .and_then(|state| state.relative_displacement(index, s, x))
            {
                at_i = v3_add(at_i, v3_scale(d, s.damping / dt * rel));
            }
            if free_i {
                add_to_particle(&mut grad, s.i, at_i);
            }
            if free_j {
                add_to_particle(&mut grad, s.j, v3_scale(at_i, -1.0));
            }
        }
        grad
    }
    /// Positive semi-definite Hessian of the incremental potential at `x`.
    pub(super) fn hessian(&self, x: &[[f64; 3]]) -> CsrMatrix {
        assemble_hessian(
            x,
            &self.masses,
            self.springs,
            self.dt,
            self.damping.as_ref(),
        )
    }
}
/// Inexact Newton minimisation of an [`ImplicitStepProblem`].
///
/// Every iteration assembles the sparse Hessian, solves `H dx = -g` with the
/// Jacobi-preconditioned conjugate gradient solver `pcg`, and applies the
/// Armijo backtracking `line_search` along `dx`.  A truncated PCG iterate
/// started from zero is still a descent direction for an SPD `H`; should the
/// computed direction ever fail to be one (non-finite values), a Jacobi-scaled
/// steepest-descent direction is used instead.
///
/// Stops when the gradient norm drops below `tolerance`, when no step along
/// the search direction decreases the energy, or after `max_iter` iterations,
/// and returns the number of Newton updates applied to `x`.
pub(super) fn newton_pcg_minimise(
    problem: &ImplicitStepProblem<'_>,
    x: &mut [[f64; 3]],
    max_iter: usize,
    tolerance: f64,
    pcg: &PcgSolver,
    line_search: &LineSearch,
) -> usize {
    let mut iterations = 0usize;
    let mut trial = x.to_vec();
    for _ in 0..max_iter {
        let grad = problem.gradient(x);
        let grad_norm = vec_norm(&grad);
        if !grad_norm.is_finite() || grad_norm < tolerance {
            break;
        }
        let hessian = problem.hessian(x);
        let neg_grad: Vec<f64> = grad.iter().map(|g| -g).collect();
        let mut direction = pcg.solve(&hessian, &neg_grad, None).x;
        let mut slope = vec_dot(&grad, &direction);
        if !(slope < 0.0 && slope.is_finite()) {
            let scaling = hessian.jacobi_preconditioner();
            direction = grad
                .iter()
                .zip(scaling.iter())
                .map(|(g, m)| -g * m)
                .collect();
            slope = vec_dot(&grad, &direction);
            if !(slope < 0.0 && slope.is_finite()) {
                break;
            }
        }
        let f0 = problem.energy(x);
        let slack = 8.0 * f64::EPSILON * (1.0 + f0.abs());
        let mut alpha = 1.0_f64;
        let mut accepted = false;
        for _ in 0..line_search.max_iter {
            for (i, t) in trial.iter_mut().enumerate() {
                *t = if problem.is_static(i) {
                    x[i]
                } else {
                    [
                        x[i][0] + alpha * direction[3 * i],
                        x[i][1] + alpha * direction[3 * i + 1],
                        x[i][2] + alpha * direction[3 * i + 2],
                    ]
                };
            }
            let f_new = problem.energy(&trial);
            if f_new <= f0 + line_search.c1 * alpha * slope + slack {
                accepted = true;
                break;
            }
            alpha *= line_search.rho;
        }
        if !accepted {
            break;
        }
        x.copy_from_slice(&trial);
        iterations += 1;
    }
    iterations
}
pub(super) const _PI: f64 = PI;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerical_softbody::BackwardEulerIntegrator;
    use crate::numerical_softbody::CsrMatrix;
    use crate::numerical_softbody::FastShapeMatchingConstraint;
    use crate::numerical_softbody::FilterLineSearch;
    use crate::numerical_softbody::ImplicitParticle;
    use crate::numerical_softbody::IpcBarrierParams;
    use crate::numerical_softbody::IpcContactSolver;
    use crate::numerical_softbody::LineSearch;
    use crate::numerical_softbody::MatrixFreeHessian;
    use crate::numerical_softbody::NewtonRaphsonSolver;
    use crate::numerical_softbody::PcgParams;
    use crate::numerical_softbody::PcgSolver;
    use crate::numerical_softbody::ProjectiveConstraint;
    use crate::numerical_softbody::ProjectiveDynamicsSolver;
    fn make_two_particle_spring(k: f64, rest: f64) -> (Vec<ImplicitParticle>, Vec<ImplicitSpring>) {
        let particles = vec![
            ImplicitParticle::new([0.0, 0.0, 0.0], 1.0),
            ImplicitParticle::new([2.0, 0.0, 0.0], 1.0),
        ];
        let springs = vec![ImplicitSpring::new(0, 1, rest, k, 0.0)];
        (particles, springs)
    }
    #[test]
    fn test_v3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = v3_add(a, b);
        assert_eq!(c, [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_v3_dot() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_eq!(v3_dot(a, b), 0.0);
        let c = [3.0, 4.0, 0.0];
        assert!((v3_norm(c) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_v3_cross() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = v3_cross(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_csr_matvec_identity() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![1.0, 1.0, 1.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let x = vec![1.0, 2.0, 3.0];
        let y = m.matvec(&x);
        assert!((y[0] - 1.0).abs() < 1e-12);
        assert!((y[1] - 2.0).abs() < 1e-12);
        assert!((y[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_csr_diagonal() {
        let rows = vec![0, 0, 1, 1, 2, 2];
        let cols = vec![0, 1, 0, 1, 1, 2];
        let vals = vec![4.0, -1.0, -1.0, 3.0, -1.0, 5.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let d = m.diagonal();
        assert!((d[0] - 4.0).abs() < 1e-12);
        assert!((d[1] - 3.0).abs() < 1e-12);
        assert!((d[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pcg_2x2() {
        let rows = vec![0, 0, 1, 1];
        let cols = vec![0, 1, 0, 1];
        let vals = vec![4.0, -1.0, -1.0, 3.0];
        let a = CsrMatrix::from_coo(2, 2, &rows, &cols, &vals);
        let b = vec![5.0, -1.0];
        let solver = PcgSolver::new(PcgParams {
            max_iter: 100,
            tolerance: 1e-10,
        });
        let result = solver.solve(&a, &b, None);
        assert!(result.converged, "PCG should converge on a 2x2 SPD system");
        assert!(
            (result.x[0] - 14.0 / 11.0).abs() < 1e-6,
            "x[0] = {}",
            result.x[0]
        );
        assert!(
            (result.x[1] - 1.0 / 11.0).abs() < 1e-6,
            "x[1] = {}",
            result.x[1]
        );
    }
    #[test]
    fn test_pcg_diagonal_3x3() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![2.0, 3.0, 4.0];
        let a = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let b = vec![6.0, 9.0, 12.0];
        let solver = PcgSolver::new(PcgParams::default());
        let result = solver.solve(&a, &b, None);
        assert!(result.converged);
        assert!((result.x[0] - 3.0).abs() < 1e-8);
        assert!((result.x[1] - 3.0).abs() < 1e-8);
        assert!((result.x[2] - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_spring_energy_at_rest() {
        let (particles, springs) = make_two_particle_spring(100.0, 2.0);
        let e = springs[0].potential_energy(&particles);
        assert!(e.abs() < 1e-12, "Energy at rest should be zero, got {e}");
    }
    #[test]
    fn test_spring_force_at_rest() {
        let (particles, springs) = make_two_particle_spring(100.0, 2.0);
        let f = springs[0].force_on_i(&particles);
        assert!(v3_norm(f) < 1e-12, "Force at rest should be zero");
    }
    #[test]
    fn test_spring_force_stretched() {
        let (particles, springs) = make_two_particle_spring(100.0, 1.0);
        let f = springs[0].force_on_i(&particles);
        assert!(
            v3_norm(f) > 0.0,
            "Stretched spring should have nonzero force"
        );
    }
    #[test]
    fn test_backward_euler_gravity() {
        let mut particles = vec![ImplicitParticle::new([0.0, 10.0, 0.0], 1.0)];
        let springs: Vec<ImplicitSpring> = vec![];
        let integrator = BackwardEulerIntegrator::new(10, 1e-6);
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            integrator.step(&mut particles, &springs, gravity, dt);
        }
        assert!(
            particles[0].position[1] < 9.0,
            "Particle should fall under gravity"
        );
    }
    #[test]
    fn test_backward_euler_static_unmoved() {
        let mut particles = vec![ImplicitParticle::new_static([5.0, 5.0, 5.0])];
        let springs: Vec<ImplicitSpring> = vec![];
        let integrator = BackwardEulerIntegrator::new(5, 1e-6);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..10 {
            integrator.step(&mut particles, &springs, gravity, 1.0 / 60.0);
        }
        let pos = particles[0].position;
        assert!((pos[0] - 5.0).abs() < 1e-12);
        assert!((pos[1] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_line_search_decreasing() {
        let ls = LineSearch::default();
        let f0 = 10.0;
        let grad = vec![1.0, 0.0];
        let alpha = ls.backtrack(f0, &grad, |a| f0 - a * 1.0);
        assert!(alpha > 0.0, "Line search should return positive step");
    }
    #[test]
    fn test_ipc_barrier_zero_outside() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        assert_eq!(p.barrier(2e-3), 0.0, "Barrier should be zero for d > d_hat");
        assert_eq!(p.barrier(-1.0), 0.0, "Barrier should be zero for d < 0");
    }
    #[test]
    fn test_ipc_barrier_positive_inside() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let b = p.barrier(5e-4);
        assert!(
            b > 0.0,
            "Barrier should be positive for d in (0, d_hat), got {b}"
        );
    }
    #[test]
    fn test_ipc_barrier_grows_near_zero() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let b1 = p.barrier(8e-4);
        let b2 = p.barrier(1e-4);
        assert!(b2 > b1, "Barrier should grow as d approaches 0");
    }
    #[test]
    fn test_ipc_barrier_gradient_sign() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let g = p.barrier_gradient(5e-4);
        assert!(
            g < 0.0,
            "Barrier gradient should be negative inside active zone, got {g}"
        );
    }
    #[test]
    fn test_ipc_total_energy() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1.0,
        };
        let distances = vec![5e-4, 8e-4, 2e-3];
        let total = p.total_barrier_energy(&distances);
        let expected = p.barrier(5e-4) + p.barrier(8e-4) + p.barrier(2e-3);
        assert!((total - expected).abs() < 1e-15);
    }
    #[test]
    fn test_pd_spring_local_step() {
        let c = ProjectiveConstraint::spring(0, 1, 1.0, 1000.0);
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let proj = c.local_step(&positions);
        let dist = v3_norm(v3_sub(proj[0], proj[1]));
        assert!(
            (dist - 1.0).abs() < 1e-10,
            "Local step should give rest length, got {dist}"
        );
    }
    #[test]
    fn test_pd_anchor_local_step() {
        let target = [1.0, 2.0, 3.0];
        let c = ProjectiveConstraint::anchor(0, target, 1000.0);
        let positions = vec![[0.0, 0.0, 0.0]];
        let proj = c.local_step(&positions);
        assert_eq!(proj[0], target);
    }
    #[test]
    fn test_pd_solver_anchor() {
        let mut solver = ProjectiveDynamicsSolver::new(1, vec![1.0], 1.0 / 60.0, 20);
        solver.add_constraint(ProjectiveConstraint::anchor(0, [0.0, 0.0, 0.0], 1e6));
        let mut positions = vec![[5.0, 0.0, 0.0]];
        let mut velocities = vec![[0.0, 0.0, 0.0]];
        for _ in 0..100 {
            solver.step(&mut positions, &mut velocities, [0.0, 0.0, 0.0]);
        }
        assert!(
            v3_norm(positions[0]) < 1.0,
            "Particle should move close to anchor"
        );
    }
    #[test]
    fn test_pd_inertia_target() {
        let solver = ProjectiveDynamicsSolver::new(1, vec![1.0], 0.1, 5);
        let pos = vec![[0.0, 0.0, 0.0]];
        let vel = vec![[1.0, 0.0, 0.0]];
        let y = solver.inertia_target(&pos, &vel, [0.0, -10.0, 0.0]);
        assert!((y[0][0] - 0.1).abs() < 1e-12);
        assert!((y[0][1] - (-0.1)).abs() < 1e-12);
    }
    #[test]
    fn test_fast_shape_matching_rest() {
        let rest = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let masses = vec![1.0, 1.0, 1.0, 1.0];
        let c = FastShapeMatchingConstraint::new(rest.clone(), masses, 1.0);
        let goals = c.goal_positions(&rest);
        for (i, (g, r)) in goals.iter().zip(rest.iter()).enumerate() {
            let diff = v3_norm(v3_sub(*g, *r));
            assert!(
                diff < 1e-4,
                "Goal {i} should equal rest pos for undeformed config, diff={diff}"
            );
        }
    }
    #[test]
    fn test_center_of_mass() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let com = FastShapeMatchingConstraint::center_of_mass(&positions, &masses);
        assert!(v3_norm(com) < 1e-12, "CoM should be at origin");
    }
    #[test]
    fn test_matrix_free_hessian_diagonal() {
        let masses = vec![2.0, 3.0];
        let springs: Vec<ImplicitSpring> = vec![];
        let hess = MatrixFreeHessian::new(springs, masses, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let v = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let hv = hess.apply(&positions, &v);
        assert!((hv[0] - 2.0).abs() < 1e-12);
        assert!((hv[4] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_spring_energy_levels() {
        let positions_rest = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let positions_stretch = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let springs = vec![ImplicitSpring::new(0, 1, 1.0, 100.0, 0.0)];
        let e_rest = spring_elastic_energy(&positions_rest, &springs);
        let e_stretch = spring_elastic_energy(&positions_stretch, &springs);
        assert!(e_rest < 1e-12, "Energy at rest should be zero");
        assert!(
            (e_stretch - 50.0).abs() < 1e-10,
            "E = 0.5*100*1^2 = 50, got {e_stretch}"
        );
    }
    #[test]
    fn test_ipc_plane_distances() {
        let params = IpcBarrierParams::default();
        let filter = FilterLineSearch::default();
        let solver = IpcContactSolver::new(params, filter);
        let positions = vec![[0.0, 1.0, 0.0], [0.0, -0.5, 0.0]];
        let normal = [0.0, 1.0, 0.0];
        let offset = 0.0;
        let d = solver.plane_distances(&positions, normal, offset);
        assert!((d[0] - 1.0).abs() < 1e-12);
        assert!((d[1] - (-0.5)).abs() < 1e-12);
    }
    #[test]
    fn test_static_mask() {
        let mut grad = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let is_static = vec![true, false];
        apply_static_mask(&mut grad, &is_static);
        assert_eq!(grad[0], 0.0);
        assert_eq!(grad[1], 0.0);
        assert_eq!(grad[2], 0.0);
        assert!((grad[3] - 4.0).abs() < 1e-12);
        assert!((grad[4] - 5.0).abs() < 1e-12);
        assert!((grad[5] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_newton_raphson_spring() {
        let (mut particles, springs) = make_two_particle_spring(100_000.0, 1.0);
        let inertia: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let pos_init: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let e0 = spring_elastic_energy(&pos_init, &springs);
        let solver = NewtonRaphsonSolver::new(50, 1e-8);
        let _iters = solver.minimise(&mut particles, &springs, &inertia, 0.1);
        let pos_final: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let e1 = spring_elastic_energy(&pos_final, &springs);
        assert!(
            e1 < e0,
            "Minimiser should reduce spring energy: e0={e0:.2}, e1={e1:.2}"
        );
    }
    #[test]
    fn test_csr_duplicate_sum() {
        let rows = vec![0, 0, 1];
        let cols = vec![0, 0, 1];
        let vals = vec![1.0, 2.0, 3.0];
        let m = CsrMatrix::from_coo(2, 2, &rows, &cols, &vals);
        let d = m.diagonal();
        assert!(
            (d[0] - 3.0).abs() < 1e-12,
            "Duplicate entries should sum, got {}",
            d[0]
        );
        assert!((d[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_filter_line_search_ccd() {
        let ls = FilterLineSearch::default();
        let positions = vec![[0.0, 1.0, 0.0]];
        let direction = vec![[0.0, -1.0, 0.0]];
        let distances = vec![1.0];
        let grads = vec![-1.0];
        let alpha = ls.ccd_step_fraction(&positions, &direction, &distances, &grads);
        assert!(
            alpha > 0.0 && alpha <= ls.tau_max,
            "CCD step should be in (0, tau_max]"
        );
    }
    #[test]
    fn test_mat3_det_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let d = mat3_det(id);
        assert!((d - 1.0).abs() < 1e-12, "det(I) = 1, got {d}");
    }
    #[test]
    fn test_v3_normalize() {
        let v = [3.0, 4.0, 0.0];
        let n = v3_normalize(v);
        assert!(
            (v3_norm(n) - 1.0).abs() < 1e-12,
            "Normalized vector should have unit length"
        );
    }
    #[test]
    fn test_jacobi_preconditioner() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![4.0, 8.0, 2.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let prec = m.jacobi_preconditioner();
        assert!((prec[0] - 0.25).abs() < 1e-12);
        assert!((prec[1] - 0.125).abs() < 1e-12);
        assert!((prec[2] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_matrix_free_cg_diagonal() {
        let positions: Vec<[f64; 3]> = vec![];
        let b = vec![4.0, 8.0, 16.0];
        let x = matrix_free_cg(
            &positions,
            &b,
            |_pos, v| vec![2.0 * v[0], 4.0 * v[1], 8.0 * v[2]],
            100,
            1e-10,
        );
        assert!((x[0] - 2.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-8, "x[1] = {}", x[1]);
        assert!((x[2] - 2.0).abs() < 1e-8, "x[2] = {}", x[2]);
    }
}
