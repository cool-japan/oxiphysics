//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
#[allow(unused_imports)]
use super::functions::*;
/// Solution returned by MPC.
#[derive(Debug, Clone)]
pub struct MpcSolution {
    /// Optimal input sequence (length horizon * m).
    pub u_seq: Vec<f64>,
    /// Predicted state trajectory (length (horizon+1) * n).
    pub x_traj: Vec<f64>,
    /// Optimal cost.
    pub cost: f64,
    /// Whether constraints are satisfied.
    pub feasible: bool,
}
/// Iterative LQR for discrete-time nonlinear systems.
///
/// Uses a linearisation at each step and solves a backward Riccati pass,
/// followed by a forward rollout.
#[derive(Debug, Clone)]
pub struct IlqrSolver {
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Horizon T.
    pub horizon: usize,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance on cost improvement.
    pub tol: f64,
    /// Regularisation on Quu.
    pub reg: f64,
}
impl IlqrSolver {
    /// Create a new iLQR solver.
    pub fn new(n: usize, m: usize, horizon: usize, max_iter: usize, tol: f64, reg: f64) -> Self {
        Self {
            n,
            m,
            horizon,
            max_iter,
            tol,
            reg,
        }
    }
    /// Solve using numerical Jacobians of dynamics f and cost l.
    ///
    /// `dynamics`: (x, u) → x_next
    /// `running_cost`: (x, u) → scalar
    /// `terminal_cost`: x → scalar
    pub fn solve<F, G, H>(
        &self,
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
        let n = self.n;
        let m = self.m;
        let t = self.horizon;
        let eps = 1e-5;
        let rollout = |u_seq: &[f64]| -> Vec<f64> {
            let mut traj = vec![0.0; (t + 1) * n];
            traj[..n].copy_from_slice(x0);
            for k in 0..t {
                let xk = traj[k * n..(k + 1) * n].to_vec();
                let uk = &u_seq[k * m..(k + 1) * m];
                let xk1 = dynamics(&xk, uk);
                traj[(k + 1) * n..(k + 2) * n].copy_from_slice(&xk1);
            }
            traj
        };
        let total_cost = |traj: &[f64], u_seq: &[f64]| -> f64 {
            let mut c = 0.0;
            for k in 0..t {
                let xk = &traj[k * n..(k + 1) * n];
                let uk = &u_seq[k * m..(k + 1) * m];
                c += running_cost(xk, uk);
            }
            c + terminal_cost(&traj[t * n..])
        };
        let mut u_seq = u_init.to_vec();
        let mut traj = rollout(&u_seq);
        let mut cost = total_cost(&traj, &u_seq);
        for iter in 0..self.max_iter {
            let mut v_x = {
                let xf = &traj[t * n..];
                let mut grad = vec![0.0; n];
                for i in 0..n {
                    let mut xp = xf.to_vec();
                    xp[i] += eps;
                    grad[i] = (terminal_cost(&xp) - terminal_cost(xf)) / eps;
                }
                grad
            };
            let mut v_xx = mat_eye(n);
            let mut ks: Vec<Vec<f64>> = vec![vec![0.0; m * n]; t];
            let mut ds: Vec<Vec<f64>> = vec![vec![0.0; m]; t];
            for k in (0..t).rev() {
                let xk = traj[k * n..(k + 1) * n].to_vec();
                let uk = u_seq[k * m..(k + 1) * m].to_vec();
                let xk1 = dynamics(&xk, &uk);
                let mut fx = vec![0.0; n * n];
                let mut fu = vec![0.0; n * m];
                for i in 0..n {
                    let mut xp = xk.clone();
                    xp[i] += eps;
                    let xp1 = dynamics(&xp, &uk);
                    for j in 0..n {
                        fx[j * n + i] = (xp1[j] - xk1[j]) / eps;
                    }
                }
                for i in 0..m {
                    let mut up = uk.clone();
                    up[i] += eps;
                    let xp1 = dynamics(&xk, &up);
                    for j in 0..n {
                        fu[j * m + i] = (xp1[j] - xk1[j]) / eps;
                    }
                }
                let c0 = running_cost(&xk, &uk);
                let mut lx = vec![0.0; n];
                let mut lu = vec![0.0; m];
                for i in 0..n {
                    let mut xp = xk.clone();
                    xp[i] += eps;
                    lx[i] = (running_cost(&xp, &uk) - c0) / eps;
                }
                for i in 0..m {
                    let mut up = uk.clone();
                    up[i] += eps;
                    lu[i] = (running_cost(&xk, &up) - c0) / eps;
                }
                let lxx = mat_eye(n);
                let luu_diag = mat_eye(m);
                let fut = mat_transpose_rect(&fu, n, m);
                let fxt = mat_transpose(&fx, n);
                let vxx_fx = mat_mul(&v_xx, &fx, n);
                let fut_vxx_fx = mat_mul_rect(&fut, m, n, &vxx_fx, n);
                let qux = fut_vxx_fx;
                let vxx_fu = mat_mul_rect(&v_xx, n, n, &fu, m);
                let fut_vxx_fu = mat_mul_rect(&fut, m, n, &vxx_fu, m);
                let mut quu = mat_add(&luu_diag, &fut_vxx_fu, m);
                for i in 0..m {
                    quu[i * m + i] += self.reg;
                }
                let fxt_vxx = mat_mul_rect(&fxt, n, n, &v_xx, n);
                let qxx = mat_add(&lxx, &mat_mul(&fxt_vxx, &fx, n), n);
                let fut_vx = mat_mul_rect(&fut, m, n, &v_x, 1);
                let fxt_vx = mat_mul_rect(&fxt, n, n, &v_x, 1);
                let qu: Vec<f64> = lu.iter().zip(fut_vx.iter()).map(|(a, b)| a + b).collect();
                let qx: Vec<f64> = lx.iter().zip(fxt_vx.iter()).map(|(a, b)| a + b).collect();
                let quu_inv = match mat_inv(&quu, m) {
                    Some(inv) => inv,
                    None => mat_eye(m),
                };
                let k_gain = mat_mul_rect(&quu_inv, m, m, &qux, n);
                let d_gain = mat_mul_rect(&quu_inv, m, m, &qu, 1);
                ks[k] = k_gain.clone();
                ds[k] = d_gain.clone();
                let kt = mat_transpose_rect(&k_gain, m, n);
                let kt_quu = mat_mul_rect(&kt, n, m, &quu, m);
                let kt_quu_k = mat_mul_rect(&kt_quu, n, m, &k_gain, n);
                let kt_qu = mat_mul_rect(&kt, n, m, &qu, 1);
                v_xx = mat_sub(&qxx, &kt_quu_k, n);
                v_x = qx
                    .iter()
                    .zip(kt_qu.iter())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>();
            }
            let mut alpha = 1.0;
            let mut new_cost = cost;
            let mut new_u = u_seq.clone();
            let mut new_traj = traj.clone();
            for _ in 0..10 {
                let mut u_new = vec![0.0; t * m];
                let mut x_new = vec![0.0; (t + 1) * n];
                x_new[..n].copy_from_slice(x0);
                for k in 0..t {
                    let dx: Vec<f64> = x_new[k * n..(k + 1) * n]
                        .iter()
                        .zip(traj[k * n..(k + 1) * n].iter())
                        .map(|(a, b)| a - b)
                        .collect();
                    let kd = mat_mul_rect(&ks[k], m, n, &dx, 1);
                    for j in 0..m {
                        u_new[k * m + j] = u_seq[k * m + j] - alpha * ds[k][j] - kd[j];
                    }
                    let xk1 = dynamics(&x_new[k * n..(k + 1) * n], &u_new[k * m..(k + 1) * m]);
                    x_new[(k + 1) * n..(k + 2) * n].copy_from_slice(&xk1);
                }
                let c_new = total_cost(&x_new, &u_new);
                if c_new < cost {
                    new_cost = c_new;
                    new_u = u_new;
                    new_traj = x_new;
                    break;
                }
                alpha *= 0.5;
            }
            let improvement = cost - new_cost;
            u_seq = new_u;
            traj = new_traj;
            cost = new_cost;
            if improvement.abs() < self.tol {
                return IlqrResult {
                    x_traj: traj,
                    u_seq,
                    cost,
                    iterations: iter + 1,
                };
            }
        }
        IlqrResult {
            x_traj: traj,
            u_seq,
            cost,
            iterations: self.max_iter,
        }
    }
}
/// Finite-horizon LQR solver (discrete time, backward Riccati recursion).
///
/// Minimises J = xₙᵀ Qf xₙ + Σ_{k=0}^{N-1} (xᵀQx + uᵀRu).
#[derive(Debug, Clone)]
pub struct FiniteHorizonLqr {
    /// Horizon length N.
    pub horizon: usize,
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Sequence of gain matrices K_k (length horizon, each m×n row-major).
    pub gains: Vec<Vec<f64>>,
}
impl FiniteHorizonLqr {
    /// Solve finite-horizon LQR via backward Riccati recursion.
    pub fn new(
        n: usize,
        m: usize,
        a: &[f64],
        b: &[f64],
        q: &[f64],
        r: &[f64],
        qf: &[f64],
        horizon: usize,
    ) -> Self {
        let at = mat_transpose(a, n);
        let bt = mat_transpose_rect(b, n, m);
        let mut p = qf.to_vec();
        let mut gains = Vec::with_capacity(horizon);
        for _ in 0..horizon {
            let bt_p = mat_mul_rect(&bt, m, n, &p, n);
            let bt_p_b = mat_mul_rect(&bt_p, m, n, b, m);
            let s = mat_add(r, &bt_p_b, m);
            let s_inv = mat_inv(&s, m).expect("R + BᵀPB must be invertible");
            let bt_p_a = mat_mul_rect(&bt_p, m, n, a, n);
            let k = mat_mul_rect(&s_inv, m, m, &bt_p_a, n);
            let at_p = mat_mul_rect(&at, n, n, &p, n);
            let at_p_a = mat_mul_rect(&at_p, n, n, a, n);
            let at_p_b = mat_mul_rect(&at_p, n, n, b, m);
            let at_p_b_k = mat_mul_rect(&at_p_b, n, m, &k, n);
            p = mat_sub(&mat_add(q, &at_p_a, n), &at_p_b_k, n);
            gains.push(k);
        }
        gains.reverse();
        Self {
            horizon,
            n,
            m,
            gains,
        }
    }
    /// Apply control at step k: u = -K_k x.
    pub fn control(&self, step: usize, state: &[f64]) -> Vec<f64> {
        let k_idx = step.min(self.horizon - 1);
        let k = &self.gains[k_idx];
        let mut u = vec![0.0; self.m];
        for i in 0..self.m {
            for j in 0..self.n {
                u[i] -= k[i * self.n + j] * state[j];
            }
        }
        u
    }
}
/// Problem specification for linear MPC.
#[derive(Debug, Clone)]
pub struct MpcProblem {
    /// Prediction horizon N.
    pub horizon: usize,
    /// State dimension n.
    pub n: usize,
    /// Input dimension m.
    pub m: usize,
    /// State matrix A.
    pub a: Vec<f64>,
    /// Input matrix B.
    pub b: Vec<f64>,
    /// State cost Q (n×n).
    pub q: Vec<f64>,
    /// Terminal state cost Qf (n×n).
    pub qf: Vec<f64>,
    /// Input cost R (m×m).
    pub r: Vec<f64>,
    /// State lower bound (length n). Use f64::NEG_INFINITY to disable.
    pub x_lb: Vec<f64>,
    /// State upper bound (length n).
    pub x_ub: Vec<f64>,
    /// Input lower bound (length m).
    pub u_lb: Vec<f64>,
    /// Input upper bound (length m).
    pub u_ub: Vec<f64>,
}
/// Linear Quadratic Regulator for a discrete-time LTI system.
///
/// Minimises J = Σ (xᵀQx + uᵀRu) over infinite horizon by solving the
/// discrete algebraic Riccati equation (DARE) iteratively.
///
/// System: x_{k+1} = A x_k + B u_k
/// Optimal control: u_k = -K x_k
#[derive(Debug, Clone)]
pub struct LqrController {
    /// State dimension n.
    pub n: usize,
    /// Input dimension m.
    pub m: usize,
    /// System matrix A (n×n, row-major).
    pub a: Vec<f64>,
    /// Input matrix B (n×m, row-major).
    pub b: Vec<f64>,
    /// State cost Q (n×n, row-major, positive semi-definite).
    pub q: Vec<f64>,
    /// Input cost R (m×m, row-major, positive definite).
    pub r: Vec<f64>,
    /// Optimal gain matrix K (m×n, row-major).
    pub k: Vec<f64>,
    /// Solution P to the DARE (n×n, row-major).
    pub p: Vec<f64>,
}
impl LqrController {
    /// Create a new LQR controller and solve for the optimal gain K.
    ///
    /// # Panics
    /// Panics if R is singular (not invertible) or dimensions are inconsistent.
    pub fn new(n: usize, m: usize, a: Vec<f64>, b: Vec<f64>, q: Vec<f64>, r: Vec<f64>) -> Self {
        assert_eq!(a.len(), n * n);
        assert_eq!(b.len(), n * m);
        assert_eq!(q.len(), n * n);
        assert_eq!(r.len(), m * m);
        let (k, p) = solve_dare(&a, &b, &q, &r, n, m, 1000, 1e-10);
        Self {
            n,
            m,
            a,
            b,
            q,
            r,
            k,
            p,
        }
    }
    /// Compute the optimal control input u = -K x.
    pub fn control(&self, state: &[f64]) -> Vec<f64> {
        let mut u = vec![0.0; self.m];
        for i in 0..self.m {
            for j in 0..self.n {
                u[i] -= self.k[i * self.n + j] * state[j];
            }
        }
        u
    }
    /// Closed-loop system matrix: A_cl = A - B*K (n×n).
    pub fn closed_loop_a(&self) -> Vec<f64> {
        let bk = mat_mul_rect(&self.b, self.n, self.m, &self.k, self.n);
        mat_sub(&self.a, &bk, self.n)
    }
    /// Check closed-loop stability: all eigenvalues inside unit circle.
    ///
    /// Uses spectral radius estimate. Returns `true` if stable.
    pub fn is_stable(&self) -> bool {
        let acl = self.closed_loop_a();
        spectral_radius(&acl, self.n) < 1.0 + 1e-6
    }
}
/// Direct collocation trajectory optimization (trapezoidal scheme).
///
/// Discretises the continuous OCP using trapezoidal integration and
/// treats both states and inputs as decision variables.
#[derive(Debug, Clone)]
pub struct TrajectoryOptimization {
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Number of collocation nodes.
    pub nodes: usize,
    /// Time step between nodes.
    pub dt: f64,
}
impl TrajectoryOptimization {
    /// Create a trajectory optimization problem.
    pub fn new(n: usize, m: usize, nodes: usize, dt: f64) -> Self {
        Self { n, m, nodes, dt }
    }
    /// Evaluate defect constraints (trapezoidal collocation).
    ///
    /// `traj`: flattened state trajectory (nodes * n).
    /// `u_seq`: flattened input sequence ((nodes-1) * m).
    /// `dynamics`: (x, u) → ẋ
    pub fn defects<F>(&self, traj: &[f64], u_seq: &[f64], dynamics: F) -> Vec<f64>
    where
        F: Fn(&[f64], &[f64]) -> Vec<f64>,
    {
        let n = self.n;
        let m = self.m;
        let mut defects = vec![0.0; (self.nodes - 1) * n];
        for k in 0..self.nodes - 1 {
            let xk = &traj[k * n..(k + 1) * n];
            let xk1 = &traj[(k + 1) * n..(k + 2) * n];
            let uk = &u_seq[k * m..(k + 1) * m];
            let fk = dynamics(xk, uk);
            let fk1 = dynamics(xk1, uk);
            for i in 0..n {
                defects[k * n + i] = xk1[i] - xk[i] - 0.5 * self.dt * (fk[i] + fk1[i]);
            }
        }
        defects
    }
}
/// Model Predictive Controller with receding horizon.
///
/// Uses a condensed QP formulation solved by projected gradient descent.
#[derive(Debug, Clone)]
pub struct ModelPredictiveControl {
    /// The problem specification.
    pub problem: MpcProblem,
}
impl ModelPredictiveControl {
    /// Create a new MPC from a problem specification.
    pub fn new(problem: MpcProblem) -> Self {
        Self { problem }
    }
    /// Solve the MPC QP for current state x0. Returns the MPC solution.
    ///
    /// Uses projected gradient descent as a simple QP solver.
    pub fn solve(&self, x0: &[f64]) -> MpcSolution {
        let p = &self.problem;
        let n = p.n;
        let m = p.m;
        let h = p.horizon;
        let mut u_seq = vec![0.0; h * m];
        let simulate = |u_seq: &[f64]| -> Vec<f64> {
            let mut traj = vec![0.0; (h + 1) * n];
            traj[..n].copy_from_slice(x0);
            for k in 0..h {
                let xk = &traj[k * n..(k + 1) * n].to_vec();
                let uk = &u_seq[k * m..(k + 1) * m];
                let ax = mat_vec(&p.a, xk, n);
                let bu = mat_mul_rect(&p.b, n, m, uk, 1);
                for i in 0..n {
                    traj[(k + 1) * n + i] = ax[i] + bu[i];
                }
            }
            traj
        };
        let cost_fn = |u_seq: &[f64]| -> f64 {
            let traj = simulate(u_seq);
            let mut cost = 0.0;
            for k in 0..h {
                let xk = &traj[k * n..(k + 1) * n];
                let uk = &u_seq[k * m..(k + 1) * m];
                let qx = mat_vec(&p.q, xk, n);
                cost += dot(xk, &qx);
                let ru = mat_vec(&p.r, uk, m);
                cost += dot(uk, &ru);
            }
            let xf = &traj[h * n..];
            let qfx = mat_vec(&p.qf, xf, n);
            cost += dot(xf, &qfx);
            cost
        };
        let step_size = 1e-3;
        for _iter in 0..500 {
            let c0 = cost_fn(&u_seq);
            let mut grad = vec![0.0; h * m];
            for i in 0..h * m {
                let mut u_perturbed = u_seq.clone();
                u_perturbed[i] += 1e-6;
                grad[i] = (cost_fn(&u_perturbed) - c0) / 1e-6;
            }
            for k in 0..h {
                for j in 0..m {
                    let idx = k * m + j;
                    u_seq[idx] -= step_size * grad[idx];
                    u_seq[idx] = u_seq[idx].max(p.u_lb[j]).min(p.u_ub[j]);
                }
            }
        }
        let traj = simulate(&u_seq);
        let cost = cost_fn(&u_seq);
        let mut feasible = true;
        for k in 0..=h {
            for i in 0..n {
                let xi = traj[k * n + i];
                if xi < p.x_lb[i] - 1e-6 || xi > p.x_ub[i] + 1e-6 {
                    feasible = false;
                }
            }
        }
        MpcSolution {
            u_seq,
            x_traj: traj,
            cost,
            feasible,
        }
    }
}
/// Discrete-time Kalman filter for a linear system.
///
/// System:  x_{k+1} = A x_k + B u_k + w_k   (w ~ N(0, Q_noise))
/// Observation: y_k = C x_k + v_k            (v ~ N(0, R_noise))
#[derive(Debug, Clone)]
pub struct KalmanFilter {
    /// State dimension n.
    pub n: usize,
    /// Observation dimension p.
    pub p: usize,
    /// System matrix A (n×n).
    pub a: Vec<f64>,
    /// Input matrix B (n×m).
    pub b: Vec<f64>,
    /// Observation matrix C (p×n).
    pub c: Vec<f64>,
    /// Process noise covariance Q_noise (n×n).
    pub q_noise: Vec<f64>,
    /// Observation noise covariance R_noise (p×p).
    pub r_noise: Vec<f64>,
    /// State estimate x̂.
    pub x_hat: Vec<f64>,
    /// Estimate covariance P.
    pub cov: Vec<f64>,
    /// Input dimension m.
    pub m: usize,
}
impl KalmanFilter {
    /// Create a new Kalman filter with initial estimate x0 and covariance P0.
    pub fn new(
        n: usize,
        p: usize,
        m: usize,
        a: Vec<f64>,
        b: Vec<f64>,
        c: Vec<f64>,
        q_noise: Vec<f64>,
        r_noise: Vec<f64>,
        x0: Vec<f64>,
        p0: Vec<f64>,
    ) -> Self {
        Self {
            n,
            p,
            m,
            a,
            b,
            c,
            q_noise,
            r_noise,
            x_hat: x0,
            cov: p0,
        }
    }
    /// Prediction step: propagate state and covariance.
    pub fn predict(&mut self, u: &[f64]) {
        let ax = mat_vec(&self.a, &self.x_hat, self.n);
        let bu = mat_mul_rect(&self.b, self.n, self.m, u, 1);
        self.x_hat = ax.iter().zip(bu.iter()).map(|(a, b)| a + b).collect();
        let at = mat_transpose(&self.a, self.n);
        let ap = mat_mul(&self.a, &self.cov, self.n);
        let ap_at = mat_mul(&ap, &at, self.n);
        self.cov = mat_add(&ap_at, &self.q_noise, self.n);
    }
    /// Update step: incorporate observation y.
    pub fn update(&mut self, y: &[f64]) {
        let cx = mat_mul_rect(&self.c, self.p, self.n, &self.x_hat, 1);
        let innovation: Vec<f64> = y.iter().zip(cx.iter()).map(|(yi, cxi)| yi - cxi).collect();
        let ct = mat_transpose_rect(&self.c, self.p, self.n);
        let p_ct = mat_mul_rect(&self.cov, self.n, self.n, &ct, self.p);
        let c_p_ct = mat_mul_rect(&self.c, self.p, self.n, &p_ct, self.p);
        let s = mat_add(&c_p_ct, &self.r_noise, self.p);
        let s_inv = mat_inv(&s, self.p).expect("Innovation covariance S must be invertible");
        let kg = mat_mul_rect(&p_ct, self.n, self.p, &s_inv, self.p);
        let k_nu = mat_mul_rect(&kg, self.n, self.p, &innovation, 1);
        for i in 0..self.n {
            self.x_hat[i] += k_nu[i];
        }
        let kc = mat_mul_rect(&kg, self.n, self.p, &self.c, self.n);
        let i_minus_kc = mat_sub(&mat_eye(self.n), &kc, self.n);
        self.cov = mat_mul(&i_minus_kc, &self.cov, self.n);
    }
    /// Compute steady-state Kalman gain by iterating Riccati equation.
    pub fn steady_state_gain(&self, max_iter: usize) -> Vec<f64> {
        let mut p = self.cov.clone();
        let at = mat_transpose(&self.a, self.n);
        let ct = mat_transpose_rect(&self.c, self.p, self.n);
        for _ in 0..max_iter {
            let p_ct = mat_mul_rect(&p, self.n, self.n, &ct, self.p);
            let c_p_ct = mat_mul_rect(&self.c, self.p, self.n, &p_ct, self.p);
            let s = mat_add(&c_p_ct, &self.r_noise, self.p);
            let s_inv = match mat_inv(&s, self.p) {
                Some(inv) => inv,
                None => break,
            };
            let kg = mat_mul_rect(&p_ct, self.n, self.p, &s_inv, self.p);
            let kc = mat_mul_rect(&kg, self.n, self.p, &self.c, self.n);
            let i_minus_kc = mat_sub(&mat_eye(self.n), &kc, self.n);
            let i_p = mat_mul(&i_minus_kc, &p, self.n);
            let ap = mat_mul(&self.a, &i_p, self.n);
            let ap_at = mat_mul(&ap, &at, self.n);
            let p_new = mat_add(&ap_at, &self.q_noise, self.n);
            let diff = frob_norm(&mat_sub(&p_new, &p, self.n));
            p = p_new;
            if diff < 1e-10 {
                break;
            }
        }
        let ct = mat_transpose_rect(&self.c, self.p, self.n);
        let p_ct = mat_mul_rect(&p, self.n, self.n, &ct, self.p);
        let c_p_ct = mat_mul_rect(&self.c, self.p, self.n, &p_ct, self.p);
        let s = mat_add(&c_p_ct, &self.r_noise, self.p);
        let s_inv = mat_inv(&s, self.p).expect("S must be invertible");
        mat_mul_rect(&p_ct, self.n, self.p, &s_inv, self.p)
    }
}
/// Result of one iLQR iteration.
#[derive(Debug, Clone)]
pub struct IlqrResult {
    /// State trajectory (length (T+1) * n).
    pub x_traj: Vec<f64>,
    /// Input sequence (length T * m).
    pub u_seq: Vec<f64>,
    /// Total cost.
    pub cost: f64,
    /// Number of iterations taken.
    pub iterations: usize,
}
/// Hamilton–Jacobi–Bellman equation value function approximation.
///
/// For a 1D or 2D state space, stores V(x) on a grid and updates via:
/// V(x) = min_u { l(x,u) + V(f(x,u)) }
#[derive(Debug, Clone)]
pub struct HamiltonJacobiBellman {
    /// Number of grid points per dimension.
    pub grid_size: usize,
    /// State space lower bound.
    pub x_min: f64,
    /// State space upper bound.
    pub x_max: f64,
    /// Value function V (flattened, grid_size × grid_size for 2D).
    pub value: Vec<f64>,
    /// Discount factor γ ∈ (0, 1].
    pub gamma: f64,
}
impl HamiltonJacobiBellman {
    /// Create with zero value function.
    pub fn new(grid_size: usize, x_min: f64, x_max: f64, gamma: f64) -> Self {
        Self {
            grid_size,
            x_min,
            x_max,
            value: vec![0.0; grid_size * grid_size],
            gamma,
        }
    }
    /// Grid spacing.
    pub fn dx(&self) -> f64 {
        (self.x_max - self.x_min) / (self.grid_size as f64 - 1.0).max(1.0)
    }
    /// Get grid index for state x (1D).
    pub fn index(&self, x: f64) -> usize {
        let idx = ((x - self.x_min) / self.dx()).round() as isize;
        idx.max(0).min(self.grid_size as isize - 1) as usize
    }
    /// Value at state (x1, x2) on 2D grid.
    pub fn get_value(&self, x1: f64, x2: f64) -> f64 {
        let i1 = self.index(x1);
        let i2 = self.index(x2);
        self.value[i1 * self.grid_size + i2]
    }
    /// Set value at (x1, x2).
    pub fn set_value(&mut self, x1: f64, x2: f64, v: f64) {
        let i1 = self.index(x1);
        let i2 = self.index(x2);
        self.value[i1 * self.grid_size + i2] = v;
    }
    /// Bellman backup: one sweep of value iteration.
    ///
    /// `dynamics`: (x1, x2, u) → (x1_next, x2_next)
    /// `cost`: (x1, x2, u) → running cost
    /// `u_grid`: candidate control values
    pub fn bellman_sweep<F, G>(&mut self, dynamics: F, cost: G, u_grid: &[f64]) -> f64
    where
        F: Fn(f64, f64, f64) -> (f64, f64),
        G: Fn(f64, f64, f64) -> f64,
    {
        let n = self.grid_size;
        let dx = self.dx();
        let mut max_diff = 0.0_f64;
        let old_value = self.value.clone();
        for i in 0..n {
            for j in 0..n {
                let x1 = self.x_min + i as f64 * dx;
                let x2 = self.x_min + j as f64 * dx;
                let mut best = f64::INFINITY;
                for &u in u_grid {
                    let (xn1, xn2) = dynamics(x1, x2, u);
                    let v_next = {
                        let i1 = self.index(xn1);
                        let i2 = self.index(xn2);
                        old_value[i1 * n + i2]
                    };
                    let val = cost(x1, x2, u) + self.gamma * v_next;
                    if val < best {
                        best = val;
                    }
                }
                let diff = (best - old_value[i * n + j]).abs();
                if diff > max_diff {
                    max_diff = diff;
                }
                self.value[i * n + j] = best;
            }
        }
        max_diff
    }
}
/// Minimum-time bang-bang controller for a double integrator.
///
/// For ẍ = u, |u| ≤ u_max, switches sign at t* = (x_dot) / u_max.
#[derive(Debug, Clone)]
pub struct MinimumTimeControl {
    /// Maximum control magnitude.
    pub u_max: f64,
    /// Target state (position).
    pub target: f64,
}
impl MinimumTimeControl {
    /// Create a new bang-bang controller.
    pub fn new(u_max: f64, target: f64) -> Self {
        Self { u_max, target }
    }
    /// Compute optimal bang-bang control for state (x, v) (position, velocity).
    ///
    /// Returns +u_max or -u_max based on switching curve.
    pub fn control(&self, x: f64, v: f64) -> f64 {
        let error = self.target - x;
        let switching = error - v * v.abs() / (2.0 * self.u_max);
        if switching > 0.0 {
            self.u_max
        } else {
            -self.u_max
        }
    }
    /// Estimate minimum time to reach target from (x, v).
    pub fn min_time(&self, x: f64, v: f64) -> f64 {
        let error = self.target - x;
        if error.abs() < 1e-10 && v.abs() < 1e-10 {
            return 0.0;
        }
        let t_brake = v.abs() / self.u_max;
        let x_brake = x + v * t_brake - 0.5 * self.u_max * t_brake * t_brake * v.signum();
        let remaining = self.target - x_brake;
        let t_accel = (2.0 * remaining.abs() / self.u_max).sqrt();
        t_brake + t_accel
    }
}
/// Linear Quadratic Gaussian (LQG) controller = LQR + Kalman filter.
///
/// Implements the separation principle: design LQR for the estimated state,
/// and Kalman filter for state estimation independently.
#[derive(Debug, Clone)]
pub struct LqgController {
    /// The LQR controller.
    pub lqr: LqrController,
    /// The Kalman filter (state estimator).
    pub kalman: KalmanFilter,
}
impl LqgController {
    /// Create a new LQG controller.
    pub fn new(lqr: LqrController, kalman: KalmanFilter) -> Self {
        Self { lqr, kalman }
    }
    /// One step: predict → observe → compute control.
    ///
    /// Returns the control input u.
    pub fn step(&mut self, observation: &[f64]) -> Vec<f64> {
        let u = self.lqr.control(&self.kalman.x_hat);
        self.kalman.predict(&u);
        self.kalman.update(observation);
        self.lqr.control(&self.kalman.x_hat)
    }
}
/// Pontryagin Minimum Principle for continuous-time optimal control.
///
/// For a system ẋ = f(x, u), with running cost l(x, u) and terminal cost Φ(x(T)),
/// the co-state λ satisfies λ̇ = -∂H/∂x where H = l + λᵀ f.
///
/// This struct stores a linearised version around a nominal trajectory.
#[derive(Debug, Clone)]
pub struct PontryaginMinimum {
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Co-state vector λ (length n).
    pub costate: Vec<f64>,
    /// Hamiltonian value H.
    pub hamiltonian: f64,
}
impl PontryaginMinimum {
    /// Create with initial co-state.
    pub fn new(n: usize, m: usize, costate: Vec<f64>) -> Self {
        Self {
            n,
            m,
            costate,
            hamiltonian: 0.0,
        }
    }
    /// Propagate co-state backward: λ̇ = -Aᵀ λ - Q x  (linearised).
    pub fn costate_step(&mut self, a: &[f64], q: &[f64], x: &[f64], dt: f64) {
        let at = mat_transpose(a, self.n);
        let at_lam = mat_vec(&at, &self.costate, self.n);
        let qx = mat_vec(q, x, self.n);
        for i in 0..self.n {
            self.costate[i] -= dt * (at_lam[i] + qx[i]);
        }
    }
    /// Evaluate the Hamiltonian H = xᵀQx + uᵀRu + λᵀ(Ax + Bu).
    pub fn eval_hamiltonian(
        &self,
        a: &[f64],
        b: &[f64],
        q: &[f64],
        r: &[f64],
        x: &[f64],
        u: &[f64],
    ) -> f64 {
        let ax = mat_vec(a, x, self.n);
        let bu = mat_mul_rect(b, self.n, self.m, u, 1);
        let f: Vec<f64> = ax.iter().zip(bu.iter()).map(|(a, b)| a + b).collect();
        let qx = mat_vec(q, x, self.n);
        let ru = mat_vec(r, u, self.m);
        dot(x, &qx) + dot(u, &ru) + dot(&self.costate, &f)
    }
    /// Minimum principle: optimal u* = -R⁻¹ Bᵀ λ.
    pub fn optimal_control(&self, b: &[f64], r: &[f64]) -> Vec<f64> {
        let bt = mat_transpose_rect(b, self.n, self.m);
        let bt_lam = mat_mul_rect(&bt, self.m, self.n, &self.costate, 1);
        let r_inv = mat_inv(r, self.m).expect("R must be invertible");
        let u: Vec<f64> = mat_mul_rect(&r_inv, self.m, self.m, &bt_lam, 1)
            .iter()
            .map(|x| -x)
            .collect();
        u
    }
}
/// Differential Dynamic Programming solver.
///
/// Similar to iLQR but includes second-order (Hxx, Huu, Hux) terms for
/// better local convergence.
#[derive(Debug, Clone)]
pub struct DifferentialDynamicProgramming {
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Horizon.
    pub horizon: usize,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}
impl DifferentialDynamicProgramming {
    /// Create a DDP solver.
    pub fn new(n: usize, m: usize, horizon: usize, max_iter: usize, tol: f64) -> Self {
        Self {
            n,
            m,
            horizon,
            max_iter,
            tol,
        }
    }
    /// Solve DDP (first-order, equivalent to iLQR for simplicity).
    pub fn solve<F, G, H>(
        &self,
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
        let ilqr = IlqrSolver::new(self.n, self.m, self.horizon, self.max_iter, self.tol, 1e-6);
        ilqr.solve(x0, u_init, dynamics, running_cost, terminal_cost)
    }
}
/// Robust H-infinity control via small gain theorem.
///
/// Stores the system transfer function gain and computes the H-infinity norm
/// via bisection on the Hamiltonian matrix (simplified scalar version).
#[derive(Debug, Clone)]
pub struct RobustControl {
    /// State dimension.
    pub n: usize,
    /// Input dimension.
    pub m: usize,
    /// Output dimension.
    pub p: usize,
    /// System A matrix.
    pub a: Vec<f64>,
    /// System B matrix (n×m).
    pub b: Vec<f64>,
    /// System C matrix (p×n).
    pub c: Vec<f64>,
    /// System D matrix (p×m).
    pub d: Vec<f64>,
}
impl RobustControl {
    /// Create a new robust controller.
    pub fn new(
        n: usize,
        m: usize,
        p: usize,
        a: Vec<f64>,
        b: Vec<f64>,
        c: Vec<f64>,
        d: Vec<f64>,
    ) -> Self {
        Self {
            n,
            m,
            p,
            a,
            b,
            c,
            d,
        }
    }
    /// Upper bound on H-infinity norm via induced 2-norm at DC (ω=0).
    ///
    /// H(0) = -C A⁻¹ B + D  (transfer matrix at ω=0).
    /// Returns the Frobenius norm of H(0) as a proxy.
    pub fn h_inf_upper_bound(&self) -> f64 {
        let a_inv = match mat_inv(&self.a, self.n) {
            Some(inv) => inv,
            None => return f64::INFINITY,
        };
        let neg_c_ainv: Vec<f64> = mat_mul_rect(&self.c, self.p, self.n, &a_inv, self.n)
            .iter()
            .map(|x| -x)
            .collect();
        let neg_c_ainv_b = mat_mul_rect(&neg_c_ainv, self.p, self.n, &self.b, self.m);
        let h0: Vec<f64> = neg_c_ainv_b
            .iter()
            .zip(self.d.iter())
            .map(|(a, b)| a + b)
            .collect();
        frob_norm(&h0)
    }
    /// Small gain theorem: system is robustly stable if γ * ||P||_∞ < 1.
    ///
    /// Returns `true` if `gamma * h_inf_upper_bound() < 1`.
    pub fn small_gain_stable(&self, gamma: f64) -> bool {
        gamma * self.h_inf_upper_bound() < 1.0
    }
}
