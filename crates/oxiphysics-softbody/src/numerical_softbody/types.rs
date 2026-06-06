//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// A spring connecting two particles, used in implicit soft body integration.
#[derive(Clone, Debug)]
pub struct ImplicitSpring {
    /// Index of the first particle.
    pub i: usize,
    /// Index of the second particle.
    pub j: usize,
    /// Rest length of the spring (m).
    pub rest_length: f64,
    /// Spring stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}
impl ImplicitSpring {
    /// Create a new spring between particles `i` and `j`.
    pub fn new(i: usize, j: usize, rest_length: f64, stiffness: f64, damping: f64) -> Self {
        ImplicitSpring {
            i,
            j,
            rest_length,
            stiffness,
            damping,
        }
    }
    /// Compute the elastic potential energy of the spring.
    ///
    /// `E = 0.5 * k * (|x_i - x_j| - L)^2`
    pub fn potential_energy(&self, particles: &[ImplicitParticle]) -> f64 {
        let pi = &particles[self.i].position;
        let pj = &particles[self.j].position;
        let diff = v3_sub(*pi, *pj);
        let len = v3_norm(diff);
        let stretch = len - self.rest_length;
        0.5 * self.stiffness * stretch * stretch
    }
    /// Compute the spring force on particle `i` (force on `j` is negated).
    ///
    /// Includes elastic and damping contributions.
    pub fn force_on_i(&self, particles: &[ImplicitParticle]) -> [f64; 3] {
        let pi = &particles[self.i].position;
        let pj = &particles[self.j].position;
        let vi = &particles[self.i].velocity;
        let vj = &particles[self.j].velocity;
        let diff = v3_sub(*pi, *pj);
        let len = v3_norm(diff);
        if len < 1e-15 {
            return [0.0; 3];
        }
        let dir = v3_scale(diff, 1.0 / len);
        let stretch = len - self.rest_length;
        let rel_vel = v3_dot(v3_sub(*vi, *vj), dir);
        let force_mag = -self.stiffness * stretch - self.damping * rel_vel;
        v3_scale(dir, force_mag)
    }
}
/// A single projective constraint.
pub struct ProjectiveConstraint {
    /// The constraint kind and parameters.
    pub kind: ProjectiveConstraintKind,
}
impl ProjectiveConstraint {
    /// Create a spring constraint.
    pub fn spring(i: usize, j: usize, rest_length: f64, weight: f64) -> Self {
        ProjectiveConstraint {
            kind: ProjectiveConstraintKind::Spring {
                i,
                j,
                rest_length,
                weight,
            },
        }
    }
    /// Create an anchor constraint.
    pub fn anchor(i: usize, target: [f64; 3], weight: f64) -> Self {
        ProjectiveConstraint {
            kind: ProjectiveConstraintKind::Anchor { i, target, weight },
        }
    }
    /// **Local step**: project positions and return the auxiliary variable `p`.
    ///
    /// For a spring constraint, `p` is the closest pair of points at rest
    /// distance.  For an anchor, `p` is the target position.
    pub fn local_step(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        match &self.kind {
            ProjectiveConstraintKind::Spring {
                i, j, rest_length, ..
            } => {
                let pi = positions[*i];
                let pj = positions[*j];
                let diff = v3_sub(pi, pj);
                let len = v3_norm(diff);
                if len < 1e-15 {
                    return vec![pi, pj];
                }
                let dir = v3_scale(diff, 1.0 / len);
                let half = v3_scale(dir, 0.5 * rest_length);
                let mid = v3_scale(v3_add(pi, pj), 0.5);
                vec![v3_add(mid, half), v3_sub(mid, half)]
            }
            ProjectiveConstraintKind::Anchor { i: _, target, .. } => vec![*target],
        }
    }
}
/// Armijo–Wolfe line search for soft body energy minimization.
///
/// Used to find a step size that sufficiently reduces the total elastic
/// energy in gradient-based solvers.
pub struct LineSearch {
    /// Armijo (sufficient decrease) parameter, typically 1e-4.
    pub c1: f64,
    /// Maximum number of backtracking iterations.
    pub max_iter: usize,
    /// Step size reduction factor (backtracking).
    pub rho: f64,
}
impl LineSearch {
    /// Backtracking line search.
    ///
    /// Given a function value `f0`, gradient-descent direction `grad` (as a
    /// flat 3n-vector), and a function callback `f`, returns the step size
    /// satisfying the Armijo condition.
    pub fn backtrack<F>(&self, f0: f64, grad: &[f64], f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut alpha = 1.0;
        let derphi0 = -vec_dot(grad, grad);
        for _ in 0..self.max_iter {
            let f_new = f(alpha);
            if f_new <= f0 + self.c1 * alpha * derphi0 {
                return alpha;
            }
            alpha *= self.rho;
        }
        alpha
    }
}
/// Filter line search for IPC (Li et al. 2020).
///
/// Ensures the step does not cause any inter-penetration (distance < 0)
/// by computing the continuous collision detection (CCD) step fraction
/// and using it as an upper bound in the backtracking search.
pub struct FilterLineSearch {
    /// Maximum CCD step fraction (0 < τ ≤ 1).
    pub tau_max: f64,
    /// Armijo sufficient decrease parameter.
    pub c1: f64,
    /// Step reduction factor.
    pub rho: f64,
    /// Maximum number of backtracking steps.
    pub max_iter: usize,
}
impl FilterLineSearch {
    /// Compute the CCD-safe step fraction.
    ///
    /// Given current positions `x`, search direction `d`, and a set of
    /// contact pair distances, returns the largest α ≤ `tau_max` such that
    /// `d(x + α d) > 0` for all contact pairs.
    ///
    /// Uses a conservative linear approximation.
    pub fn ccd_step_fraction(
        &self,
        _x: &[[f64; 3]],
        _d: &[[f64; 3]],
        distances: &[f64],
        distance_gradients: &[f64],
    ) -> f64 {
        let mut alpha = self.tau_max;
        for (&dist, &grad) in distances.iter().zip(distance_gradients.iter()) {
            if grad < 0.0 {
                let t = -0.9 * dist / grad;
                if t < alpha {
                    alpha = t;
                }
            }
        }
        alpha.max(1e-10)
    }
    /// Perform filter line search with energy and CCD constraint.
    ///
    /// Returns the accepted step size.
    pub fn search<F>(&self, alpha_ccd: f64, f0: f64, grad: &[f64], f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut alpha = alpha_ccd.min(1.0);
        let derphi0 = -vec_dot(grad, grad);
        for _ in 0..self.max_iter {
            let f_new = f(alpha);
            if f_new <= f0 + self.c1 * alpha * derphi0 {
                return alpha;
            }
            alpha *= self.rho;
        }
        alpha
    }
}
/// Fast simulation: treats shape matching as a projective constraint.
///
/// This maps the shape-matching local step (polar decomposition) into the
/// projective dynamics framework, enabling fast position-based simulation
/// that is equivalent to the corotational elasticity model in the limit.
///
/// Reference: Müller et al., "Fast Simulation of Inextensible Hair and Fur",
/// adapted to the projective dynamics framework.
pub struct FastShapeMatchingConstraint {
    /// Rest-shape particle positions (relative to center of mass).
    pub rest_positions: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Shape matching stiffness weight.
    pub weight: f64,
}
impl FastShapeMatchingConstraint {
    /// Create a new fast shape matching constraint.
    pub fn new(rest_positions: Vec<[f64; 3]>, masses: Vec<f64>, weight: f64) -> Self {
        FastShapeMatchingConstraint {
            rest_positions,
            masses,
            weight,
        }
    }
    /// Compute center of mass of a set of positions.
    pub fn center_of_mass(positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
        let total_mass: f64 = masses.iter().sum();
        if total_mass < 1e-15 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for (p, &m) in positions.iter().zip(masses.iter()) {
            com = v3_add(com, v3_scale(*p, m));
        }
        v3_scale(com, 1.0 / total_mass)
    }
    /// Compute the 3×3 covariance matrix A_pq.
    ///
    /// `A_pq = sum_i m_i (q_i - q_cm)(r_i - r_cm)^T`
    fn covariance_matrix(
        positions: &[[f64; 3]],
        rest_positions: &[[f64; 3]],
        masses: &[f64],
        com_cur: [f64; 3],
        com_rest: [f64; 3],
    ) -> [[f64; 3]; 3] {
        let mut apq = [[0.0f64; 3]; 3];
        for ((p, r), &m) in positions
            .iter()
            .zip(rest_positions.iter())
            .zip(masses.iter())
        {
            let qi = v3_sub(*p, com_cur);
            let ri = v3_sub(*r, com_rest);
            for row in 0..3 {
                for col in 0..3 {
                    apq[row][col] += m * qi[row] * ri[col];
                }
            }
        }
        apq
    }
    /// Polar decomposition of a 3×3 matrix to extract rotation (iterative).
    ///
    /// Uses the iterative method: R_{k+1} = (R_k + R_k^{-T}) / 2.
    /// If the matrix is near-singular (det ≈ 0) the identity is returned.
    fn polar_decomp_rotation(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let det0 = mat3_det(a);
        if det0.abs() < 1e-10 {
            return id;
        }
        let mut r = a;
        for _ in 0..20 {
            let det = mat3_det(r);
            if det.abs() < 1e-15 {
                return id;
            }
            let inv_t = mat3_transpose(mat3_inv_scale(r, det));
            let mut converged = true;
            for i in 0..3 {
                for j in 0..3 {
                    let new_val = 0.5 * (r[i][j] + inv_t[i][j]);
                    if (new_val - r[i][j]).abs() > 1e-10 {
                        converged = false;
                    }
                    r[i][j] = new_val;
                }
            }
            if converged {
                break;
            }
        }
        r
    }
    /// Local step: compute goal positions using the shape-matched rotation.
    pub fn goal_positions(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let com_cur = Self::center_of_mass(positions, &self.masses);
        let com_rest = Self::center_of_mass(&self.rest_positions, &self.masses);
        let apq = Self::covariance_matrix(
            positions,
            &self.rest_positions,
            &self.masses,
            com_cur,
            com_rest,
        );
        let r = Self::polar_decomp_rotation(apq);
        positions
            .iter()
            .zip(self.rest_positions.iter())
            .map(|(_p, r_rest)| {
                let rel = v3_sub(*r_rest, com_rest);
                let rotated = mat3_mul_v3(r, rel);
                v3_add(com_cur, rotated)
            })
            .collect()
    }
}
/// Result of a PCG solve.
pub struct PcgResult {
    /// Solution vector.
    pub x: Vec<f64>,
    /// Final residual norm.
    pub residual_norm: f64,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Backward Euler time integrator for spring-mass soft bodies.
///
/// Solves the implicit system:
///
/// ```text
/// M (v_{n+1} - v_n) / dt = f(x_{n+1})
/// x_{n+1} = x_n + dt * v_{n+1}
/// ```
///
/// Linearised via Newton-Raphson at each time step.
pub struct BackwardEulerIntegrator {
    /// Maximum Newton-Raphson iterations per step.
    pub max_newton_iter: usize,
    /// Newton convergence tolerance.
    pub newton_tol: f64,
    /// PCG solver for the linear system.
    pub pcg: PcgSolver,
}
impl BackwardEulerIntegrator {
    /// Create a new backward Euler integrator.
    pub fn new(max_newton_iter: usize, newton_tol: f64) -> Self {
        BackwardEulerIntegrator {
            max_newton_iter,
            newton_tol,
            pcg: PcgSolver::new(PcgParams {
                max_iter: 500,
                tolerance: 1e-8,
            }),
        }
    }
    /// Perform one backward Euler step for a spring-mass system.
    ///
    /// Returns the number of Newton iterations performed.
    pub fn step(
        &self,
        particles: &mut [ImplicitParticle],
        springs: &[ImplicitSpring],
        gravity: [f64; 3],
        dt: f64,
    ) -> usize {
        let n = particles.len();
        let dof = 3 * n;
        let x_prev: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let v_prev: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
        let mut x_star: Vec<[f64; 3]> = x_prev
            .iter()
            .zip(v_prev.iter())
            .map(|(x, v)| v3_add(*x, v3_scale(*v, dt)))
            .collect();
        for i in 0..n {
            particles[i].position = x_star[i];
        }
        let mut iter_count = 0;
        for _iter in 0..self.max_newton_iter {
            iter_count += 1;
            let mut residual = vec![0.0f64; dof];
            for i in 0..n {
                if particles[i].is_static() {
                    continue;
                }
                let m = particles[i].mass;
                let v_new = v3_scale(v3_sub(particles[i].position, x_prev[i]), 1.0 / dt);
                let inertia_force = v3_scale(v3_sub(v_new, v_prev[i]), m / dt);
                let fg = v3_scale(gravity, m);
                let mut fs = [0.0f64; 3];
                for s in springs {
                    if s.i == i {
                        let f = s.force_on_i(particles);
                        fs = v3_add(fs, f);
                    } else if s.j == i {
                        let f = s.force_on_i(particles);
                        fs = v3_sub(fs, f);
                    }
                }
                let net = v3_add(fg, fs);
                let g = v3_sub(inertia_force, net);
                residual[3 * i] = g[0];
                residual[3 * i + 1] = g[1];
                residual[3 * i + 2] = g[2];
            }
            let res_norm = vec_norm(&residual);
            if res_norm < self.newton_tol {
                break;
            }
            let mut diag = vec![0.0f64; dof];
            for i in 0..n {
                if particles[i].is_static() {
                    diag[3 * i] = 1.0;
                    diag[3 * i + 1] = 1.0;
                    diag[3 * i + 2] = 1.0;
                } else {
                    let m = particles[i].mass;
                    diag[3 * i] = m / (dt * dt);
                    diag[3 * i + 1] = m / (dt * dt);
                    diag[3 * i + 2] = m / (dt * dt);
                }
            }
            for s in springs {
                let k = s.stiffness;
                for idx in [s.i, s.j] {
                    if !particles[idx].is_static() {
                        diag[3 * idx] += k;
                        diag[3 * idx + 1] += k;
                        diag[3 * idx + 2] += k;
                    }
                }
            }
            let neg_res: Vec<f64> = residual.iter().map(|r| -r).collect();
            let delta_x: Vec<f64> = neg_res
                .iter()
                .zip(diag.iter())
                .map(|(r, d)| r / d.max(1e-15))
                .collect();
            for i in 0..n {
                if particles[i].is_static() {
                    continue;
                }
                x_star[i] = v3_add(
                    x_star[i],
                    [delta_x[3 * i], delta_x[3 * i + 1], delta_x[3 * i + 2]],
                );
                particles[i].position = x_star[i];
            }
        }
        for i in 0..n {
            if !particles[i].is_static() {
                particles[i].velocity =
                    v3_scale(v3_sub(particles[i].position, x_prev[i]), 1.0 / dt);
            }
        }
        iter_count
    }
}
/// Newton-Raphson solver for implicit soft body systems.
///
/// Solves the nonlinear optimality condition arising from backward Euler
/// discretisation of the elastic potential energy.  Uses the PCG solver
/// internally for the linear sub-problem.
pub struct NewtonRaphsonSolver {
    /// Maximum outer Newton iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the gradient norm.
    pub tolerance: f64,
    /// PCG parameters.
    pub pcg_params: PcgParams,
    /// Line search.
    pub line_search: LineSearch,
}
impl NewtonRaphsonSolver {
    /// Create a Newton-Raphson solver with custom parameters.
    pub fn new(max_iter: usize, tolerance: f64) -> Self {
        NewtonRaphsonSolver {
            max_iter,
            tolerance,
            pcg_params: PcgParams::default(),
            line_search: LineSearch::default(),
        }
    }
    /// Perform Newton-Raphson minimisation on a spring-mass system.
    ///
    /// Minimises `E(x) = sum springs + 0.5/dt^2 * ||M(x - y)||^2`
    /// where `y` is the inertia target.
    pub fn minimise(
        &self,
        particles: &mut [ImplicitParticle],
        springs: &[ImplicitSpring],
        inertia_target: &[[f64; 3]],
        dt: f64,
    ) -> usize {
        let n = particles.len();
        let mut iter_count = 0;
        for _iter in 0..self.max_iter {
            iter_count += 1;
            let mut grad = vec![0.0f64; 3 * n];
            for i in 0..n {
                if particles[i].is_static() {
                    continue;
                }
                let m = particles[i].mass;
                let diff = v3_sub(particles[i].position, inertia_target[i]);
                grad[3 * i] += m / (dt * dt) * diff[0];
                grad[3 * i + 1] += m / (dt * dt) * diff[1];
                grad[3 * i + 2] += m / (dt * dt) * diff[2];
            }
            for s in springs {
                if particles[s.i].is_static() && particles[s.j].is_static() {
                    continue;
                }
                let pi = particles[s.i].position;
                let pj = particles[s.j].position;
                let diff = v3_sub(pi, pj);
                let len = v3_norm(diff);
                if len < 1e-15 {
                    continue;
                }
                let dir = v3_scale(diff, 1.0 / len);
                let stretch = len - s.rest_length;
                let fg = v3_scale(dir, s.stiffness * stretch);
                if !particles[s.i].is_static() {
                    grad[3 * s.i] += fg[0];
                    grad[3 * s.i + 1] += fg[1];
                    grad[3 * s.i + 2] += fg[2];
                }
                if !particles[s.j].is_static() {
                    grad[3 * s.j] -= fg[0];
                    grad[3 * s.j + 1] -= fg[1];
                    grad[3 * s.j + 2] -= fg[2];
                }
            }
            let grad_norm = vec_norm(&grad);
            if grad_norm < self.tolerance {
                break;
            }
            let mut diag = vec![1.0f64; 3 * n];
            for i in 0..n {
                if particles[i].is_static() {
                    continue;
                }
                let base = particles[i].mass / (dt * dt);
                diag[3 * i] = base;
                diag[3 * i + 1] = base;
                diag[3 * i + 2] = base;
            }
            for s in springs {
                let k = s.stiffness;
                for idx in [s.i, s.j] {
                    if !particles[idx].is_static() {
                        diag[3 * idx] += k;
                        diag[3 * idx + 1] += k;
                        diag[3 * idx + 2] += k;
                    }
                }
            }
            let direction: Vec<f64> = grad
                .iter()
                .zip(diag.iter())
                .map(|(g, d)| -g / d.max(1e-15))
                .collect();
            let current_energy = self.total_energy(particles, springs, inertia_target, dt);
            let derphi0 = vec_dot(&grad, &direction);
            let c1 = self.line_search.c1;
            let rho = self.line_search.rho;
            let mut alpha = 1.0_f64;
            for _ in 0..self.line_search.max_iter {
                let mut test_particles = particles.to_owned();
                for i in 0..n {
                    if test_particles[i].is_static() {
                        continue;
                    }
                    let dp = [
                        direction[3 * i] * alpha,
                        direction[3 * i + 1] * alpha,
                        direction[3 * i + 2] * alpha,
                    ];
                    test_particles[i].position = v3_add(test_particles[i].position, dp);
                }
                let f_new = self.total_energy(&test_particles, springs, inertia_target, dt);
                if f_new <= current_energy + c1 * alpha * derphi0 {
                    break;
                }
                alpha *= rho;
            }
            for i in 0..n {
                if particles[i].is_static() {
                    continue;
                }
                let dp = [
                    direction[3 * i] * alpha,
                    direction[3 * i + 1] * alpha,
                    direction[3 * i + 2] * alpha,
                ];
                particles[i].position = v3_add(particles[i].position, dp);
            }
        }
        iter_count
    }
    /// Compute the total implicit energy for a particle system.
    fn total_energy(
        &self,
        particles: &[ImplicitParticle],
        springs: &[ImplicitSpring],
        inertia_target: &[[f64; 3]],
        dt: f64,
    ) -> f64 {
        let mut e = 0.0;
        for (i, p) in particles.iter().enumerate() {
            if p.is_static() {
                continue;
            }
            let diff = v3_sub(p.position, inertia_target[i]);
            e += 0.5 * p.mass / (dt * dt) * v3_dot(diff, diff);
        }
        for s in springs {
            e += s.potential_energy(particles);
        }
        e
    }
}
/// Type of projective constraint supported.
#[derive(Clone, Debug)]
pub enum ProjectiveConstraintKind {
    /// Spring / distance constraint (stretch).
    Spring {
        /// Index of first particle.
        i: usize,
        /// Index of second particle.
        j: usize,
        /// Rest length.
        rest_length: f64,
        /// Stiffness weight.
        weight: f64,
    },
    /// Anchor / positional constraint (pinned vertex).
    Anchor {
        /// Index of the particle.
        i: usize,
        /// Target world-space position.
        target: [f64; 3],
        /// Stiffness weight.
        weight: f64,
    },
}
/// Sparse matrix in Compressed Sparse Row (CSR) format.
///
/// Stores only non-zero entries. Used as the system matrix for PCG solvers
/// applied to deformable body dynamics.
pub struct CsrMatrix {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row pointer array of length `nrows + 1`.
    pub row_ptr: Vec<usize>,
    /// Column index array.
    pub col_idx: Vec<usize>,
    /// Value array.
    pub values: Vec<f64>,
}
impl CsrMatrix {
    /// Construct a CSR matrix from COO (coordinate) data.
    ///
    /// Entries with the same `(row, col)` pair are **summed** (assembled).
    pub fn from_coo(
        nrows: usize,
        ncols: usize,
        rows: &[usize],
        cols: &[usize],
        vals: &[f64],
    ) -> Self {
        assert_eq!(rows.len(), cols.len());
        assert_eq!(rows.len(), vals.len());
        let mut row_count = vec![0usize; nrows];
        for &r in rows {
            row_count[r] += 1;
        }
        let mut row_ptr = vec![0usize; nrows + 1];
        for i in 0..nrows {
            row_ptr[i + 1] = row_ptr[i] + row_count[i];
        }
        let nnz = row_ptr[nrows];
        let mut col_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        let mut cursor = row_ptr.clone();
        for k in 0..rows.len() {
            let r = rows[k];
            let pos = cursor[r];
            col_idx[pos] = cols[k];
            values[pos] = vals[k];
            cursor[r] += 1;
        }
        let mut m = CsrMatrix {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            values,
        };
        m.sort_and_sum_duplicates();
        m
    }
    /// Sort column indices within each row and sum duplicate entries.
    fn sort_and_sum_duplicates(&mut self) {
        let nrows = self.nrows;
        let mut new_col: Vec<usize> = Vec::new();
        let mut new_val: Vec<f64> = Vec::new();
        let mut new_ptr = vec![0usize; nrows + 1];
        for r in 0..nrows {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            let mut pairs: Vec<(usize, f64)> = self.col_idx[start..end]
                .iter()
                .zip(self.values[start..end].iter())
                .map(|(&c, &v)| (c, v))
                .collect();
            pairs.sort_by_key(|p| p.0);
            let mut j = 0;
            while j < pairs.len() {
                let c = pairs[j].0;
                let mut v = pairs[j].1;
                let mut k = j + 1;
                while k < pairs.len() && pairs[k].0 == c {
                    v += pairs[k].1;
                    k += 1;
                }
                new_col.push(c);
                new_val.push(v);
                j = k;
            }
            new_ptr[r + 1] = new_col.len();
        }
        self.row_ptr = new_ptr;
        self.col_idx = new_col;
        self.values = new_val;
    }
    /// Sparse matrix–vector product: y = A * x
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols);
        let mut y = vec![0.0f64; self.nrows];
        for (r, y_r) in y.iter_mut().enumerate() {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            for k in start..end {
                *y_r += self.values[k] * x[self.col_idx[k]];
            }
        }
        y
    }
    /// Extract diagonal vector.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0f64; self.nrows];
        for (r, d_r) in d.iter_mut().enumerate() {
            let start = self.row_ptr[r];
            let end = self.row_ptr[r + 1];
            for k in start..end {
                if self.col_idx[k] == r {
                    *d_r = self.values[k];
                }
            }
        }
        d
    }
    /// Construct a diagonal (Jacobi) preconditioner M⁻¹.
    ///
    /// Returns a vector where each entry is `1 / diag[i]`.
    /// Zero diagonal entries are treated as 1 to avoid division by zero.
    pub fn jacobi_preconditioner(&self) -> Vec<f64> {
        self.diagonal()
            .iter()
            .map(|&d| if d.abs() > 1e-15 { 1.0 / d } else { 1.0 })
            .collect()
    }
}
/// State of a soft body particle for implicit time integration.
///
/// Stores position, velocity, and mass for use in backward Euler or other
/// implicit schemes.
#[derive(Clone, Debug)]
pub struct ImplicitParticle {
    /// Current position in world space (m).
    pub position: [f64; 3],
    /// Current velocity (m/s).
    pub velocity: [f64; 3],
    /// Particle mass (kg). Zero means the particle is static (pinned).
    pub mass: f64,
    /// External force accumulator (N).
    pub external_force: [f64; 3],
}
impl ImplicitParticle {
    /// Create a new dynamic particle at the given position.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        ImplicitParticle {
            position,
            velocity: [0.0; 3],
            mass,
            external_force: [0.0; 3],
        }
    }
    /// Create a static (pinned) particle.
    pub fn new_static(position: [f64; 3]) -> Self {
        ImplicitParticle {
            position,
            velocity: [0.0; 3],
            mass: 0.0,
            external_force: [0.0; 3],
        }
    }
    /// Return true if this particle is pinned (zero mass).
    pub fn is_static(&self) -> bool {
        self.mass < 1e-15
    }
}
/// Preconditioned Conjugate Gradient solver for symmetric positive-definite systems.
///
/// Solves `A x = b` using the PCG method with a diagonal (Jacobi) preconditioner.
/// This is the core linear solver used in implicit soft-body time integration.
///
/// # References
/// - Shewchuk, "An Introduction to the Conjugate Gradient Method Without the
///   Agonizing Pain", 1994.
pub struct PcgSolver {
    /// Solver parameters.
    pub params: PcgParams,
}
impl PcgSolver {
    /// Create a new PCG solver with the given parameters.
    pub fn new(params: PcgParams) -> Self {
        PcgSolver { params }
    }
    /// Solve `A x = b` with an optional initial guess `x0`.
    ///
    /// Uses the Jacobi (diagonal) preconditioner extracted from `A`.
    pub fn solve(&self, a: &CsrMatrix, b: &[f64], x0: Option<&[f64]>) -> PcgResult {
        let n = b.len();
        assert_eq!(a.nrows, n);
        assert_eq!(a.ncols, n);
        let m_inv = a.jacobi_preconditioner();
        let mut x: Vec<f64> = match x0 {
            Some(x0) => x0.to_vec(),
            None => vec![0.0f64; n],
        };
        let ax = a.matvec(&x);
        let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, axi)| bi - axi).collect();
        let mut z: Vec<f64> = r.iter().zip(m_inv.iter()).map(|(ri, mi)| ri * mi).collect();
        let mut p = z.clone();
        let mut rz = vec_dot(&r, &z);
        let b_norm = vec_norm(b);
        let tol = if b_norm > 1e-15 {
            self.params.tolerance * b_norm
        } else {
            self.params.tolerance
        };
        for iter in 0..self.params.max_iter {
            let ap = a.matvec(&p);
            let pap = vec_dot(&p, &ap);
            if pap.abs() < 1e-30 {
                break;
            }
            let alpha = rz / pap;
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            for i in 0..n {
                r[i] -= alpha * ap[i];
            }
            let res_norm = vec_norm(&r);
            if res_norm < tol {
                return PcgResult {
                    x,
                    residual_norm: res_norm,
                    iterations: iter + 1,
                    converged: true,
                };
            }
            for i in 0..n {
                z[i] = r[i] * m_inv[i];
            }
            let rz_new = vec_dot(&r, &z);
            let beta = rz_new / rz.max(1e-30);
            rz = rz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
        }
        let res_norm = vec_norm(&r);
        PcgResult {
            x,
            residual_norm: res_norm,
            iterations: self.params.max_iter,
            converged: false,
        }
    }
}
/// Projective Dynamics solver.
///
/// Implements the global/local splitting scheme from Bouaziz et al. 2014.
/// The global step solves a fixed SPD linear system, while the local step
/// independently projects each constraint.
///
/// This enables real-time simulation of elastic solids at low iteration counts.
pub struct ProjectiveDynamicsSolver {
    /// Number of particles (n).
    pub n_particles: usize,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// List of projective constraints.
    pub constraints: Vec<ProjectiveConstraint>,
    /// Time step size (s).
    pub dt: f64,
    /// Number of local/global iterations per time step.
    pub n_iterations: usize,
}
impl ProjectiveDynamicsSolver {
    /// Create a new Projective Dynamics solver.
    pub fn new(n_particles: usize, masses: Vec<f64>, dt: f64, n_iterations: usize) -> Self {
        ProjectiveDynamicsSolver {
            n_particles,
            masses,
            constraints: Vec::new(),
            dt,
            n_iterations,
        }
    }
    /// Add a constraint to the solver.
    pub fn add_constraint(&mut self, c: ProjectiveConstraint) {
        self.constraints.push(c);
    }
    /// Compute the inertia target: `y = q_n + dt * v_n + dt^2 * M^{-1} f_ext`
    pub fn inertia_target(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        gravity: [f64; 3],
    ) -> Vec<[f64; 3]> {
        let dt = self.dt;
        positions
            .iter()
            .zip(velocities.iter())
            .zip(self.masses.iter())
            .map(|((q, v), &m)| {
                let inertia = v3_add(*q, v3_scale(*v, dt));
                if m > 1e-15 {
                    v3_add(inertia, v3_scale(gravity, dt * dt))
                } else {
                    *q
                }
            })
            .collect()
    }
    /// One step of Projective Dynamics.
    ///
    /// Alternates between local (constraint projection) and global
    /// (system solve) steps for `n_iterations` rounds.
    pub fn step(&self, positions: &mut [[f64; 3]], velocities: &mut [[f64; 3]], gravity: [f64; 3]) {
        let n = self.n_particles;
        let dt = self.dt;
        let y = self.inertia_target(positions, velocities, gravity);
        let mut q = y.clone();
        for _iter in 0..self.n_iterations {
            let mut rhs: Vec<[f64; 3]> = (0..n)
                .map(|i| v3_scale(y[i], self.masses[i] / (dt * dt)))
                .collect();
            for c in &self.constraints {
                let proj = c.local_step(&q);
                match &c.kind {
                    ProjectiveConstraintKind::Spring { i, j, weight, .. } => {
                        let w_dt2 = weight / (dt * dt);
                        rhs[*i] = v3_add(rhs[*i], v3_scale(proj[0], w_dt2));
                        rhs[*j] = v3_add(rhs[*j], v3_scale(proj[1], w_dt2));
                    }
                    ProjectiveConstraintKind::Anchor { i, weight, .. } => {
                        let w_dt2 = weight / (dt * dt);
                        rhs[*i] = v3_add(rhs[*i], v3_scale(proj[0], w_dt2));
                    }
                }
            }
            for i in 0..n {
                let mut diag = self.masses[i] / (dt * dt);
                for c in &self.constraints {
                    match &c.kind {
                        ProjectiveConstraintKind::Spring {
                            i: ci,
                            j: cj,
                            weight,
                            ..
                        } => {
                            if *ci == i || *cj == i {
                                diag += weight / (dt * dt);
                            }
                        }
                        ProjectiveConstraintKind::Anchor { i: ci, weight, .. } => {
                            if *ci == i {
                                diag += weight / (dt * dt);
                            }
                        }
                    }
                }
                if diag > 1e-15 {
                    q[i] = v3_scale(rhs[i], 1.0 / diag);
                }
            }
        }
        for i in 0..n {
            if self.masses[i] > 1e-15 {
                velocities[i] = v3_scale(v3_sub(q[i], positions[i]), 1.0 / dt);
                positions[i] = q[i];
            }
        }
    }
}
/// IPC barrier energy parameters.
///
/// The barrier function prevents inter-penetration by adding an energy that
/// grows to infinity as the contact distance approaches zero.
///
/// Reference: Li et al., "Incremental Potential Contact: Intersection- and
/// Inversion-free Large Deformation Dynamics", SIGGRAPH 2020.
pub struct IpcBarrierParams {
    /// Activation distance: barrier activates when gap < d_hat (m).
    pub d_hat: f64,
    /// Barrier stiffness multiplier κ.
    pub kappa: f64,
}
impl IpcBarrierParams {
    /// Evaluate the IPC barrier function b(d, d̂).
    ///
    /// `b(d, d̂) = -(d - d̂)^2 * ln(d / d̂)` for `0 < d < d̂`,
    /// else `0`.
    pub fn barrier(&self, d: f64) -> f64 {
        if d <= 0.0 || d >= self.d_hat {
            return 0.0;
        }
        let x = d / self.d_hat;
        -(d - self.d_hat).powi(2) * x.ln()
    }
    /// Derivative of the barrier function with respect to `d`.
    pub fn barrier_gradient(&self, d: f64) -> f64 {
        if d <= 0.0 || d >= self.d_hat {
            return 0.0;
        }
        let dh = self.d_hat;
        -2.0 * (d - dh) * (d / dh).ln() - (d - dh).powi(2) / d
    }
    /// Second derivative of the barrier function.
    pub fn barrier_hessian(&self, d: f64) -> f64 {
        if d <= 0.0 || d >= self.d_hat {
            return 0.0;
        }
        let dh = self.d_hat;
        -2.0 * (d / dh).ln() - 2.0 * (d - dh) / d - 2.0 * (d - dh) / d + (d - dh).powi(2) / (d * d)
    }
    /// Total barrier energy for a list of contact distances.
    pub fn total_barrier_energy(&self, distances: &[f64]) -> f64 {
        self.kappa * distances.iter().map(|&d| self.barrier(d)).sum::<f64>()
    }
}
/// Parameters for the Preconditioned Conjugate Gradient solver.
pub struct PcgParams {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the residual norm.
    pub tolerance: f64,
}
/// Matrix-free Hessian-vector product for spring-mass systems.
///
/// Computes `H * v` without explicitly forming the Hessian matrix, where
/// `H` is the Hessian of the total implicit energy.  Used in matrix-free PCG.
pub struct MatrixFreeHessian {
    /// Spring definitions.
    pub springs: Vec<ImplicitSpring>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Time step.
    pub dt: f64,
}
impl MatrixFreeHessian {
    /// Create a new matrix-free Hessian.
    pub fn new(springs: Vec<ImplicitSpring>, masses: Vec<f64>, dt: f64) -> Self {
        MatrixFreeHessian {
            springs,
            masses,
            dt,
        }
    }
    /// Compute `H * v` for the implicit system Hessian.
    ///
    /// `H = M/dt^2 + K` where `K` is the stiffness matrix of the springs.
    pub fn apply(&self, positions: &[[f64; 3]], v: &[f64]) -> Vec<f64> {
        let n = positions.len();
        let mut result = vec![0.0f64; 3 * n];
        for i in 0..n {
            let s = self.masses[i] / (self.dt * self.dt);
            result[3 * i] += s * v[3 * i];
            result[3 * i + 1] += s * v[3 * i + 1];
            result[3 * i + 2] += s * v[3 * i + 2];
        }
        for s in &self.springs {
            let pi = positions[s.i];
            let pj = positions[s.j];
            let diff = v3_sub(pi, pj);
            let len = v3_norm(diff);
            if len < 1e-15 {
                continue;
            }
            let d = v3_scale(diff, 1.0 / len);
            let ratio = s.rest_length / len;
            let vi = [v[3 * s.i], v[3 * s.i + 1], v[3 * s.i + 2]];
            let vj = [v[3 * s.j], v[3 * s.j + 1], v[3 * s.j + 2]];
            let dv = v3_sub(vi, vj);
            let ddot = v3_dot(d, dv);
            let ki: [f64; 3] = [
                s.stiffness * (d[0] * ddot + (1.0 - ratio) * (dv[0] - d[0] * ddot)),
                s.stiffness * (d[1] * ddot + (1.0 - ratio) * (dv[1] - d[1] * ddot)),
                s.stiffness * (d[2] * ddot + (1.0 - ratio) * (dv[2] - d[2] * ddot)),
            ];
            result[3 * s.i] += ki[0];
            result[3 * s.i + 1] += ki[1];
            result[3 * s.i + 2] += ki[2];
            result[3 * s.j] -= ki[0];
            result[3 * s.j + 1] -= ki[1];
            result[3 * s.j + 2] -= ki[2];
        }
        result
    }
}
/// Simple IPC contact solver for point-plane contacts.
///
/// Computes barrier energies and gradients for a list of particle-plane
/// contact pairs, and applies the filter line search to avoid penetration.
pub struct IpcContactSolver {
    /// IPC barrier parameters.
    pub params: IpcBarrierParams,
    /// Filter line search.
    pub filter: FilterLineSearch,
}
impl IpcContactSolver {
    /// Create a new IPC contact solver.
    pub fn new(params: IpcBarrierParams, filter: FilterLineSearch) -> Self {
        IpcContactSolver { params, filter }
    }
    /// Compute contact distances for particle-plane contacts.
    ///
    /// Each plane is defined by a normal `n` and offset `d`:
    /// `distance = x · n - d`.
    pub fn plane_distances(
        &self,
        positions: &[[f64; 3]],
        plane_normal: [f64; 3],
        plane_offset: f64,
    ) -> Vec<f64> {
        positions
            .iter()
            .map(|p| v3_dot(*p, plane_normal) - plane_offset)
            .collect()
    }
    /// Compute per-particle barrier energy gradient components.
    ///
    /// Returns a flat 3n vector of the barrier energy gradient.
    pub fn barrier_gradient_flat(
        &self,
        positions: &[[f64; 3]],
        plane_normal: [f64; 3],
        plane_offset: f64,
    ) -> Vec<f64> {
        let n = positions.len();
        let mut grad = vec![0.0f64; 3 * n];
        for (i, p) in positions.iter().enumerate() {
            let d = v3_dot(*p, plane_normal) - plane_offset;
            let db = self.params.barrier_gradient(d);
            let scale = self.params.kappa * db;
            grad[3 * i] += scale * plane_normal[0];
            grad[3 * i + 1] += scale * plane_normal[1];
            grad[3 * i + 2] += scale * plane_normal[2];
        }
        grad
    }
}
