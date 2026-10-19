//! Iterative LQR (iLQR) and Differential Dynamic Programming (DDP) for
//! discrete-time nonlinear systems `x_{k+1} = f(x_k, u_k)` with running cost
//! `l(x, u)` and terminal cost `φ(x)`.
//!
//! Every iteration linearises the dynamics and expands the cost to second
//! order along the current trajectory (central finite differences), runs a
//! backward Riccati-style recursion on the Q-function
//!
//! ```text
//! Qx  = lx  + Fxᵀ Vx                      Qu  = lu  + Fuᵀ Vx
//! Qxx = lxx + Fxᵀ Vxx Fx  [+ Σ Vx_i ∇²ₓₓ f_i]
//! Quu = luu + Fuᵀ Vxx Fu  [+ Σ Vx_i ∇²ᵤᵤ f_i]
//! Qux = lux + Fuᵀ Vxx Fx  [+ Σ Vx_i ∇²ᵤₓ f_i]
//! ```
//!
//! where the bracketed second derivatives of the dynamics are only used by
//! DDP, and finishes with a forward line search on the nonlinear dynamics.
//! The feedback law is `u = ū - α k - K (x - x̄)` with `k = Quu⁻¹ Qu` and
//! `K = Quu⁻¹ Qux` (a Levenberg–Marquardt term `μ I` is added to `Quu`).
//!
//! The regularisation `μ` is adaptive: it grows whenever `Quu + μ I` is not
//! positive definite (Cholesky fails) or no step length produces a sufficient
//! decrease, and it shrinks after every accepted step, so the iteration
//! behaves like a Gauss–Newton / Newton method near a minimum and like a
//! damped gradient method far from it.

use super::functions::{mat_mul_rect, mat_transpose_rect};
use super::types::IlqrResult;
use crate::dense_linalg::{cholesky_factor, cholesky_solve};

/// Relative step of the central finite differences (scaled by `max(1, |z|)`).
const FD_STEP: f64 = 1e-4;
/// Line-search step lengths are `1, 1/2, ..., 2^-(LINE_SEARCH_STEPS - 1)`.
const LINE_SEARCH_STEPS: usize = 11;
/// A step is accepted when the cost decreases by at least this fraction of the
/// decrease predicted by the quadratic model.
const MIN_REDUCTION_RATIO: f64 = 1e-4;
/// Smallest non-zero regularisation reached when growing from zero.
const MU_FLOOR: f64 = 1e-6;
/// Growth / shrink factor of the regularisation schedule.
const MU_FACTOR: f64 = 2.0;
/// The iteration gives up when the regularisation exceeds this value.
const MU_MAX: f64 = 1e10;
/// Relative gradient norm below which the trajectory is stationary.
const STATIONARY_GRADIENT: f64 = 1e-10;

/// Order of the expansion of the dynamics used in the backward pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DynamicsOrder {
    /// Linearised dynamics only: iLQR (Gauss–Newton).
    First,
    /// Linearised dynamics plus the second derivatives contracted with the
    /// value-function gradient: DDP.
    Second,
}

/// Dimensions and tuning of one solve.
#[derive(Debug, Clone, Copy)]
pub(super) struct IlqrSettings {
    pub(super) n: usize,
    pub(super) m: usize,
    pub(super) horizon: usize,
    pub(super) max_iter: usize,
    pub(super) tol: f64,
    pub(super) reg: f64,
    pub(super) order: DynamicsOrder,
}

/// Relative finite-difference step for a coordinate of value `v`.
fn fd_step(v: f64) -> f64 {
    FD_STEP * v.abs().max(1.0)
}

/// Central-difference gradient and Hessian of `f` at `z`.
///
/// The Hessian is assembled from the diagonal second difference
/// `(f(z+h) - 2 f(z) + f(z-h)) / h²` and the four-point mixed formula
/// `(f(++) - f(+-) - f(-+) + f(--)) / (4 h_i h_j)`, so it is symmetric by
/// construction and exact (up to rounding) for quadratic `f`.
pub(super) fn gradient_and_hessian<F>(f: &F, z: &[f64]) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(&[f64]) -> f64,
{
    let d = z.len();
    let h: Vec<f64> = z.iter().map(|&v| fd_step(v)).collect();
    let f0 = f(z);
    let mut work = z.to_vec();
    let mut grad = vec![0.0; d];
    let mut hess = vec![0.0; d * d];
    for i in 0..d {
        work[i] = z[i] + h[i];
        let plus = f(&work);
        work[i] = z[i] - h[i];
        let minus = f(&work);
        work[i] = z[i];
        grad[i] = (plus - minus) / (2.0 * h[i]);
        hess[i * d + i] = (plus - 2.0 * f0 + minus) / (h[i] * h[i]);
    }
    for i in 0..d {
        for j in (i + 1)..d {
            let mut corner = |si: f64, sj: f64| {
                work[i] = z[i] + si * h[i];
                work[j] = z[j] + sj * h[j];
                let value = f(&work);
                work[i] = z[i];
                work[j] = z[j];
                value
            };
            let (pp, pm, mp, mm) = (
                corner(1.0, 1.0),
                corner(1.0, -1.0),
                corner(-1.0, 1.0),
                corner(-1.0, -1.0),
            );
            let value = (pp - pm - mp + mm) / (4.0 * h[i] * h[j]);
            hess[i * d + j] = value;
            hess[j * d + i] = value;
        }
    }
    (grad, hess)
}

/// Central-difference Jacobian `∂f/∂z` (row-major `n x d`) and, when requested,
/// the Hessian of every component `f_r` (`d x d` each) at `z`.
///
/// `f` maps the `d = len(z)` stacked inputs to the `n` outputs.  All Hessians
/// are obtained from the same function evaluations.
fn vector_derivatives<F>(
    f: &F,
    z: &[f64],
    n: usize,
    second_order: bool,
) -> (Vec<f64>, Vec<Vec<f64>>)
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let d = z.len();
    let h: Vec<f64> = z.iter().map(|&v| fd_step(v)).collect();
    let mut work = z.to_vec();
    let mut jac = vec![0.0; n * d];
    let mut plus: Vec<Vec<f64>> = Vec::with_capacity(d);
    let mut minus: Vec<Vec<f64>> = Vec::with_capacity(d);
    for i in 0..d {
        work[i] = z[i] + h[i];
        plus.push(f(&work));
        work[i] = z[i] - h[i];
        minus.push(f(&work));
        work[i] = z[i];
        for r in 0..n {
            jac[r * d + i] = (plus[i][r] - minus[i][r]) / (2.0 * h[i]);
        }
    }
    let mut hessians: Vec<Vec<f64>> = Vec::new();
    if second_order {
        let f0 = f(z);
        hessians = vec![vec![0.0; d * d]; n];
        for i in 0..d {
            for r in 0..n {
                hessians[r][i * d + i] = (plus[i][r] - 2.0 * f0[r] + minus[i][r]) / (h[i] * h[i]);
            }
        }
        for i in 0..d {
            for j in (i + 1)..d {
                let mut corner = |si: f64, sj: f64| {
                    work[i] = z[i] + si * h[i];
                    work[j] = z[j] + sj * h[j];
                    let value = f(&work);
                    work[i] = z[i];
                    work[j] = z[j];
                    value
                };
                let (pp, pm, mp, mm) = (
                    corner(1.0, 1.0),
                    corner(1.0, -1.0),
                    corner(-1.0, 1.0),
                    corner(-1.0, -1.0),
                );
                for r in 0..n {
                    let value = (pp[r] - pm[r] - mp[r] + mm[r]) / (4.0 * h[i] * h[j]);
                    hessians[r][i * d + j] = value;
                    hessians[r][j * d + i] = value;
                }
            }
        }
    }
    (jac, hessians)
}

/// First- and second-order expansion of dynamics and cost at one stage.
pub(super) struct StageDerivatives {
    /// `∂f/∂x` (`n x n`).
    pub(super) fx: Vec<f64>,
    /// `∂f/∂u` (`n x m`).
    pub(super) fu: Vec<f64>,
    /// Hessians of the dynamics components w.r.t. `[x; u]` (`(n+m) x (n+m)`
    /// each); empty unless the expansion is second order.
    pub(super) fzz: Vec<Vec<f64>>,
    pub(super) lx: Vec<f64>,
    pub(super) lu: Vec<f64>,
    pub(super) lxx: Vec<f64>,
    pub(super) luu: Vec<f64>,
    /// `∂²l/∂u∂x` (`m x n`).
    pub(super) lux: Vec<f64>,
}

/// Derivatives along a whole trajectory.
pub(super) struct TrajectoryDerivatives {
    pub(super) stages: Vec<StageDerivatives>,
    /// Gradient of the terminal cost at the final state.
    pub(super) vx_terminal: Vec<f64>,
    /// Hessian of the terminal cost at the final state.
    pub(super) vxx_terminal: Vec<f64>,
}

/// Feedback policy from one backward pass.
pub(super) struct Policy {
    /// Feed-forward terms `k_k = Quu⁻¹ Qu` (length `m` each).
    pub(super) k_ff: Vec<Vec<f64>>,
    /// Feedback gains `K_k = Quu⁻¹ Qux` (`m x n` each).
    pub(super) gains: Vec<Vec<f64>>,
    /// `Σ Quᵀ k`: the predicted decrease is `α dv1 - α²/2 dv2`.
    pub(super) dv1: f64,
    /// `Σ kᵀ Quu k`.
    pub(super) dv2: f64,
    /// Largest `|Qu|` entry, which is the true gradient of the total cost with
    /// respect to the controls along the nominal trajectory.
    pub(super) grad_norm: f64,
}

/// `a + b` for equally long slices.
fn add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Replace the row-major `n x n` matrix by `(A + Aᵀ) / 2`.
fn symmetrize(a: &mut [f64], n: usize) {
    for i in 0..n {
        for j in (i + 1)..n {
            let avg = 0.5 * (a[i * n + j] + a[j * n + i]);
            a[i * n + j] = avg;
            a[j * n + i] = avg;
        }
    }
}

/// The nonlinear optimal-control problem handed to the solver.
struct Ocp<'a, F, G, H> {
    n: usize,
    m: usize,
    horizon: usize,
    x0: &'a [f64],
    dynamics: F,
    running_cost: G,
    terminal_cost: H,
}

impl<F, G, H> Ocp<'_, F, G, H>
where
    F: Fn(&[f64], &[f64]) -> Vec<f64>,
    G: Fn(&[f64], &[f64]) -> f64,
    H: Fn(&[f64]) -> f64,
{
    /// One dynamics step; an output of the wrong length is padded with NaN so
    /// that the cost becomes non-finite and the step is rejected.
    fn step(&self, x: &[f64], u: &[f64]) -> Vec<f64> {
        let mut next = (self.dynamics)(x, u);
        next.resize(self.n, f64::NAN);
        next
    }
    fn rollout(&self, u_seq: &[f64]) -> Vec<f64> {
        let (n, m) = (self.n, self.m);
        let mut traj = vec![0.0; (self.horizon + 1) * n];
        traj[..n].copy_from_slice(self.x0);
        for k in 0..self.horizon {
            let next = self.step(&traj[k * n..(k + 1) * n], &u_seq[k * m..(k + 1) * m]);
            traj[(k + 1) * n..(k + 2) * n].copy_from_slice(&next);
        }
        traj
    }
    fn total_cost(&self, traj: &[f64], u_seq: &[f64]) -> f64 {
        let (n, m) = (self.n, self.m);
        let running: f64 = (0..self.horizon)
            .map(|k| (self.running_cost)(&traj[k * n..(k + 1) * n], &u_seq[k * m..(k + 1) * m]))
            .sum();
        running + (self.terminal_cost)(&traj[self.horizon * n..])
    }
    fn derivatives(
        &self,
        traj: &[f64],
        u_seq: &[f64],
        order: DynamicsOrder,
    ) -> TrajectoryDerivatives {
        let (n, m) = (self.n, self.m);
        let d = n + m;
        let second_order = order == DynamicsOrder::Second;
        let stages = (0..self.horizon)
            .map(|k| {
                let mut z = Vec::with_capacity(d);
                z.extend_from_slice(&traj[k * n..(k + 1) * n]);
                z.extend_from_slice(&u_seq[k * m..(k + 1) * m]);
                let (jac, fzz) = vector_derivatives(
                    &|zz: &[f64]| self.step(&zz[..n], &zz[n..]),
                    &z,
                    n,
                    second_order,
                );
                let (grad, hess) =
                    gradient_and_hessian(&|zz: &[f64]| (self.running_cost)(&zz[..n], &zz[n..]), &z);
                let mut fx = vec![0.0; n * n];
                let mut fu = vec![0.0; n * m];
                for r in 0..n {
                    fx[r * n..(r + 1) * n].copy_from_slice(&jac[r * d..r * d + n]);
                    fu[r * m..(r + 1) * m].copy_from_slice(&jac[r * d + n..(r + 1) * d]);
                }
                let mut lxx = vec![0.0; n * n];
                let mut luu = vec![0.0; m * m];
                let mut lux = vec![0.0; m * n];
                for i in 0..n {
                    lxx[i * n..(i + 1) * n].copy_from_slice(&hess[i * d..i * d + n]);
                }
                for i in 0..m {
                    luu[i * m..(i + 1) * m]
                        .copy_from_slice(&hess[(n + i) * d + n..(n + i + 1) * d]);
                    lux[i * n..(i + 1) * n].copy_from_slice(&hess[(n + i) * d..(n + i) * d + n]);
                }
                StageDerivatives {
                    fx,
                    fu,
                    fzz,
                    lx: grad[..n].to_vec(),
                    lu: grad[n..].to_vec(),
                    lxx,
                    luu,
                    lux,
                }
            })
            .collect();
        let (vx_terminal, vxx_terminal) =
            gradient_and_hessian(&self.terminal_cost, &traj[self.horizon * n..]);
        TrajectoryDerivatives {
            stages,
            vx_terminal,
            vxx_terminal,
        }
    }
    /// Roll the policy out on the nonlinear dynamics with step length `alpha`.
    fn forward_pass(
        &self,
        traj: &[f64],
        u_seq: &[f64],
        policy: &Policy,
        alpha: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let (n, m) = (self.n, self.m);
        let mut x_new = vec![0.0; (self.horizon + 1) * n];
        let mut u_new = vec![0.0; self.horizon * m];
        x_new[..n].copy_from_slice(self.x0);
        for k in 0..self.horizon {
            let dx: Vec<f64> = (0..n).map(|i| x_new[k * n + i] - traj[k * n + i]).collect();
            let feedback = mat_mul_rect(&policy.gains[k], m, n, &dx, 1);
            for j in 0..m {
                u_new[k * m + j] = u_seq[k * m + j] - alpha * policy.k_ff[k][j] - feedback[j];
            }
            let next = self.step(&x_new[k * n..(k + 1) * n], &u_new[k * m..(k + 1) * m]);
            x_new[(k + 1) * n..(k + 2) * n].copy_from_slice(&next);
        }
        (x_new, u_new)
    }
}

/// Backward pass with regularisation `mu`; `None` if `Quu + mu I` is not
/// positive definite at some stage.
pub(super) fn backward_pass(
    s: &IlqrSettings,
    derivs: &TrajectoryDerivatives,
    mu: f64,
) -> Option<Policy> {
    let (n, m) = (s.n, s.m);
    let d = n + m;
    let horizon = derivs.stages.len();
    let mut vx = derivs.vx_terminal.clone();
    let mut vxx = derivs.vxx_terminal.clone();
    symmetrize(&mut vxx, n);
    let mut k_ff = vec![Vec::new(); horizon];
    let mut gains = vec![Vec::new(); horizon];
    let (mut dv1, mut dv2, mut grad_norm) = (0.0, 0.0, 0.0_f64);
    for k in (0..horizon).rev() {
        let st = &derivs.stages[k];
        let fx_t = mat_transpose_rect(&st.fx, n, n);
        let fu_t = mat_transpose_rect(&st.fu, n, m);
        let fxt_vxx = mat_mul_rect(&fx_t, n, n, &vxx, n);
        let fut_vxx = mat_mul_rect(&fu_t, m, n, &vxx, n);
        let mut qxx = add(&st.lxx, &mat_mul_rect(&fxt_vxx, n, n, &st.fx, n));
        let mut quu = add(&st.luu, &mat_mul_rect(&fut_vxx, m, n, &st.fu, m));
        let mut qux = add(&st.lux, &mat_mul_rect(&fut_vxx, m, n, &st.fx, n));
        let fxt_vx = mat_mul_rect(&fx_t, n, n, &vx, 1);
        let fut_vx = mat_mul_rect(&fu_t, m, n, &vx, 1);
        let qx = add(&st.lx, &fxt_vx);
        let qu = add(&st.lu, &fut_vx);
        if s.order == DynamicsOrder::Second {
            for (r, hess) in st.fzz.iter().enumerate() {
                let weight = vx[r];
                if weight == 0.0 {
                    continue;
                }
                for i in 0..n {
                    for j in 0..n {
                        qxx[i * n + j] += weight * hess[i * d + j];
                    }
                }
                for i in 0..m {
                    for j in 0..m {
                        quu[i * m + j] += weight * hess[(n + i) * d + n + j];
                    }
                    for j in 0..n {
                        qux[i * n + j] += weight * hess[(n + i) * d + j];
                    }
                }
            }
        }
        symmetrize(&mut qxx, n);
        symmetrize(&mut quu, m);
        let mut quu_reg = quu.clone();
        for i in 0..m {
            quu_reg[i * m + i] += mu;
        }
        let factor = cholesky_factor(&quu_reg, m)?;
        let k_k = cholesky_solve(&factor, m, &qu, 1);
        let gain = cholesky_solve(&factor, m, &qux, n);
        let gain_t = mat_transpose_rect(&gain, m, n);
        let qux_t = mat_transpose_rect(&qux, m, n);
        // Vx  = Qx + Kᵀ Quu k - Kᵀ Qu - Quxᵀ k
        let quu_k = mat_mul_rect(&quu, m, m, &k_k, 1);
        let kt_term: Vec<f64> = quu_k.iter().zip(qu.iter()).map(|(a, b)| a - b).collect();
        let kt_vec = mat_mul_rect(&gain_t, n, m, &kt_term, 1);
        let qux_k = mat_mul_rect(&qux_t, n, m, &k_k, 1);
        let vx_new: Vec<f64> = (0..n).map(|i| qx[i] + kt_vec[i] - qux_k[i]).collect();
        // Vxx = Qxx + Kᵀ Quu K - Kᵀ Qux - Quxᵀ K
        let quu_gain = mat_mul_rect(&quu, m, m, &gain, n);
        let inner: Vec<f64> = quu_gain
            .iter()
            .zip(qux.iter())
            .map(|(a, b)| a - b)
            .collect();
        let kt_inner = mat_mul_rect(&gain_t, n, m, &inner, n);
        let qux_t_gain = mat_mul_rect(&qux_t, n, m, &gain, n);
        let mut vxx_new: Vec<f64> = (0..n * n)
            .map(|i| qxx[i] + kt_inner[i] - qux_t_gain[i])
            .collect();
        symmetrize(&mut vxx_new, n);
        dv1 += qu.iter().zip(k_k.iter()).map(|(a, b)| a * b).sum::<f64>();
        dv2 += k_k
            .iter()
            .zip(quu_k.iter())
            .map(|(a, b)| a * b)
            .sum::<f64>();
        grad_norm = qu.iter().fold(grad_norm, |acc, v| acc.max(v.abs()));
        vx = vx_new;
        vxx = vxx_new;
        k_ff[k] = k_k;
        gains[k] = gain;
    }
    Some(Policy {
        k_ff,
        gains,
        dv1,
        dv2,
        grad_norm,
    })
}

/// Raise the regularisation after a failed backward pass or line search;
/// `false` once it exceeds [`MU_MAX`].
fn grow_regularisation(mu: &mut f64, scale: &mut f64) -> bool {
    *scale = (*scale * MU_FACTOR).max(MU_FACTOR);
    *mu = (*mu * *scale).max(MU_FLOOR);
    *mu <= MU_MAX
}

/// Relax the regularisation after an accepted step (never below `mu_min`).
fn shrink_regularisation(mu: &mut f64, scale: &mut f64, mu_min: f64) {
    *scale = (*scale / MU_FACTOR).min(1.0 / MU_FACTOR);
    *mu = (*mu * *scale).max(mu_min);
}

/// Solve the optimal-control problem by iLQR or DDP (see the module docs).
///
/// `x0` must have `n` entries; `u_init` holds `horizon * m` controls (missing
/// entries count as zero, surplus entries are ignored).
pub(super) fn solve<F, G, H>(
    s: &IlqrSettings,
    x0: &[f64],
    u_init: &[f64],
    dynamics: F,
    running_cost: G,
    terminal_cost: H,
) -> IlqrResult
where
    F: Fn(&[f64], &[f64]) -> Vec<f64>,
    G: Fn(&[f64], &[f64]) -> f64,
    H: Fn(&[f64]) -> f64,
{
    let (n, m, horizon) = (s.n, s.m, s.horizon);
    let mut start = vec![0.0; n];
    for (dst, src) in start.iter_mut().zip(x0.iter()) {
        *dst = *src;
    }
    let ocp = Ocp {
        n,
        m,
        horizon,
        x0: &start,
        dynamics,
        running_cost,
        terminal_cost,
    };
    let mut u_seq = vec![0.0; horizon * m];
    for (dst, src) in u_seq.iter_mut().zip(u_init.iter()) {
        *dst = *src;
    }
    let mut traj = ocp.rollout(&u_seq);
    let mut cost = ocp.total_cost(&traj, &u_seq);
    let no_gains = || vec![vec![0.0; m * n]; horizon];
    if horizon == 0 || m == 0 || n == 0 {
        return IlqrResult {
            x_traj: traj,
            u_seq,
            cost,
            iterations: 0,
            gains: no_gains(),
            converged: true,
        };
    }
    if !cost.is_finite() {
        return IlqrResult {
            x_traj: traj,
            u_seq,
            cost,
            iterations: 0,
            gains: no_gains(),
            converged: false,
        };
    }

    let mu_min = s.reg.max(0.0);
    let mut mu = mu_min;
    let mut mu_scale = 1.0;
    let mut converged = false;
    let mut iterations = 0;
    let mut derivs = ocp.derivatives(&traj, &u_seq, s.order);
    let mut derivs_current = true;
    'iterations: for iter in 0..s.max_iter {
        iterations = iter + 1;
        if !derivs_current {
            derivs = ocp.derivatives(&traj, &u_seq, s.order);
            derivs_current = true;
        }
        loop {
            let Some(policy) = backward_pass(s, &derivs, mu) else {
                if !grow_regularisation(&mut mu, &mut mu_scale) {
                    break 'iterations;
                }
                continue;
            };
            if policy.grad_norm <= STATIONARY_GRADIENT * (1.0 + cost.abs()) {
                converged = true;
                break 'iterations;
            }
            let mut accepted = false;
            for step in 0..LINE_SEARCH_STEPS {
                let alpha = 0.5_f64.powi(step as i32);
                let (x_new, u_new) = ocp.forward_pass(&traj, &u_seq, &policy, alpha);
                let new_cost = ocp.total_cost(&x_new, &u_new);
                if !new_cost.is_finite() || new_cost >= cost {
                    continue;
                }
                let expected = alpha * policy.dv1 - 0.5 * alpha * alpha * policy.dv2;
                if expected > 0.0 && (cost - new_cost) / expected < MIN_REDUCTION_RATIO {
                    continue;
                }
                let improvement = cost - new_cost;
                traj = x_new;
                u_seq = u_new;
                cost = new_cost;
                derivs_current = false;
                shrink_regularisation(&mut mu, &mut mu_scale, mu_min);
                converged = improvement <= s.tol;
                accepted = true;
                break;
            }
            if accepted {
                if converged {
                    break 'iterations;
                }
                break;
            }
            if policy.dv1 - 0.5 * policy.dv2 <= s.tol {
                // The model predicts no worthwhile decrease: stationary.
                converged = true;
                break 'iterations;
            }
            if !grow_regularisation(&mut mu, &mut mu_scale) {
                break 'iterations;
            }
        }
    }

    // Feedback gains along the returned trajectory.
    if !derivs_current {
        derivs = ocp.derivatives(&traj, &u_seq, s.order);
    }
    let mut gain_mu = mu;
    let mut gain_scale = 1.0;
    let gains = loop {
        match backward_pass(s, &derivs, gain_mu) {
            Some(policy) => break policy.gains,
            None => {
                if !grow_regularisation(&mut gain_mu, &mut gain_scale) {
                    break no_gains();
                }
            }
        }
    };
    IlqrResult {
        x_traj: traj,
        u_seq,
        cost,
        iterations,
        gains,
        converged,
    }
}

/// Backward pass along the trajectory generated by `u_seq` (test support).
#[cfg(test)]
pub(super) fn policy_at<F, G, H>(
    s: &IlqrSettings,
    x0: &[f64],
    u_seq: &[f64],
    dynamics: F,
    running_cost: G,
    terminal_cost: H,
    mu: f64,
) -> Option<Policy>
where
    F: Fn(&[f64], &[f64]) -> Vec<f64>,
    G: Fn(&[f64], &[f64]) -> f64,
    H: Fn(&[f64]) -> f64,
{
    let ocp = Ocp {
        n: s.n,
        m: s.m,
        horizon: s.horizon,
        x0,
        dynamics,
        running_cost,
        terminal_cost,
    };
    let traj = ocp.rollout(u_seq);
    let derivs = ocp.derivatives(&traj, u_seq, s.order);
    backward_pass(s, &derivs, mu)
}

/// Dynamics Jacobian and Hessians of `f` at `z` (test support).
#[cfg(test)]
pub(super) fn dynamics_derivatives<F>(
    f: &F,
    z: &[f64],
    n: usize,
    second_order: bool,
) -> (Vec<f64>, Vec<Vec<f64>>)
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    vector_derivatives(f, z, n, second_order)
}
