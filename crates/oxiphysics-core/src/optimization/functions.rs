//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{OptResult, Particle};
use rand::RngExt;

/// Gradient descent with fixed learning rate.
///
/// Minimises `f` starting from `x0` using the gradient `grad`.
/// Stops when the L2 norm of the gradient is below `tol`.
pub fn gradient_descent(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lr: f64,
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let mut x = x0;
    let mut converged = false;
    let mut n_iter = 0u32;
    for _ in 0..max_iter {
        n_iter += 1;
        let g = grad(&x);
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        for (xi, gi) in x.iter_mut().zip(g.iter()) {
            *xi -= lr * gi;
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Hyperparameters for the Adam adaptive-moment optimizer.
#[derive(Debug, Clone, Copy)]
pub struct AdamConfig {
    /// Learning rate (step size).
    pub lr: f64,
    /// Exponential decay rate for the first moment estimate.
    pub beta1: f64,
    /// Exponential decay rate for the second moment estimate.
    pub beta2: f64,
    /// Small constant for numerical stability in the denominator.
    pub eps: f64,
    /// Maximum number of iterations.
    pub max_iter: u32,
}

impl Default for AdamConfig {
    fn default() -> Self {
        Self {
            lr: 0.001,
            beta1: 0.9,
            beta2: 0.999,
            eps: 1e-8,
            max_iter: 1000,
        }
    }
}

/// Adam adaptive-moment optimizer.
///
/// Uses the standard Adam update rule with bias correction.
/// Stops when the L2 norm of the gradient falls below the default tolerance
/// `1e-7`, or after `cfg.max_iter` steps.
pub fn adam(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    cfg: AdamConfig,
) -> OptResult {
    let AdamConfig {
        lr,
        beta1,
        beta2,
        eps,
        max_iter,
    } = cfg;
    let n = x0.len();
    let mut x = x0;
    let mut m = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let mut converged = false;
    let mut n_iter = 0u32;
    for t in 1..=max_iter {
        n_iter = t;
        let g = grad(&x);
        let gnorm: f64 = g.iter().map(|gi| gi * gi).sum::<f64>().sqrt();
        if gnorm < 1e-7 {
            converged = true;
            break;
        }
        let b1t = beta1.powi(t as i32);
        let b2t = beta2.powi(t as i32);
        for i in 0..n {
            m[i] = beta1 * m[i] + (1.0 - beta1) * g[i];
            v[i] = beta2 * v[i] + (1.0 - beta2) * g[i] * g[i];
            let m_hat = m[i] / (1.0 - b1t);
            let v_hat = v[i] / (1.0 - b2t);
            x[i] -= lr * m_hat / (v_hat.sqrt() + eps);
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}
pub(super) fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}
/// L-BFGS (Limited-memory Broyden–Fletcher–Goldfarb–Shanno).
///
/// Uses a history of size `m` curvature pairs and Wolfe-condition line search.
pub fn lbfgs(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    m: usize,
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let _n = x0.len();
    let mut x = x0;
    let mut g = grad(&x);
    let mut s_list: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut y_list: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut converged = false;
    let mut n_iter = 0u32;
    for _ in 0..max_iter {
        n_iter += 1;
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        let k = s_list.len();
        let mut q = g.clone();
        let mut alphas = vec![0.0f64; k];
        for i in (0..k).rev() {
            let rho_i = 1.0 / dot(&y_list[i], &s_list[i]);
            let a = rho_i * dot(&s_list[i], &q);
            alphas[i] = a;
            axpy(-a, &y_list[i].clone(), &mut q);
        }
        let mut r = if k > 0 {
            let scale = dot(&s_list[k - 1], &y_list[k - 1])
                / dot(&y_list[k - 1], &y_list[k - 1]).max(1e-15);
            q.iter().map(|qi| scale * qi).collect::<Vec<_>>()
        } else {
            q.clone()
        };
        for i in 0..k {
            let rho_i = 1.0 / dot(&y_list[i], &s_list[i]);
            let beta = rho_i * dot(&y_list[i], &r);
            axpy(alphas[i] - beta, &s_list[i].clone(), &mut r);
        }
        let direction: Vec<f64> = r.iter().map(|ri| -*ri).collect();
        let f0 = f(&x);
        let gd = dot(&g, &direction);
        let alpha =
            backtracking_line_search(&f, &x, &direction, f0, gd, BacktrackingConfig::default());
        let s: Vec<f64> = direction.iter().map(|di| alpha * di).collect();
        let x_new: Vec<f64> = x.iter().zip(s.iter()).map(|(xi, si)| xi + si).collect();
        let g_new = grad(&x_new);
        let y: Vec<f64> = g_new
            .iter()
            .zip(g.iter())
            .map(|(gni, gi)| gni - gi)
            .collect();
        let sy = dot(&s, &y);
        if sy > 1e-15 {
            if s_list.len() == m {
                s_list.remove(0);
                y_list.remove(0);
            }
            s_list.push(s);
            y_list.push(y);
        }
        x = x_new;
        g = g_new;
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Nelder-Mead simplex method (derivative-free).
pub fn nelder_mead(
    f: impl Fn(&[f64]) -> f64,
    x0: Vec<f64>,
    step: f64,
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let n = x0.len();
    let mut simplex: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    simplex.push(x0.clone());
    for i in 0..n {
        let mut v = x0.clone();
        v[i] += step;
        simplex.push(v);
    }
    let mut fvals: Vec<f64> = simplex.iter().map(|v| f(v)).collect();
    let mut n_iter = 0u32;
    let mut converged = false;
    let alpha = 1.0_f64;
    let gamma = 2.0_f64;
    let rho = 0.5_f64;
    let sigma = 0.5_f64;
    for _ in 0..max_iter {
        n_iter += 1;
        let mut order: Vec<usize> = (0..=n).collect();
        order.sort_by(|&a, &b| {
            fvals[a]
                .partial_cmp(&fvals[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        simplex = order.iter().map(|&i| simplex[i].clone()).collect();
        fvals = order.iter().map(|&i| fvals[i]).collect();
        let fmean = fvals.iter().sum::<f64>() / (n + 1) as f64;
        let fstd: f64 =
            (fvals.iter().map(|fi| (fi - fmean).powi(2)).sum::<f64>() / (n + 1) as f64).sqrt();
        if fstd < tol {
            converged = true;
            break;
        }
        let mut centroid = vec![0.0f64; n];
        for v in &simplex[..n] {
            for (c, vi) in centroid.iter_mut().zip(v.iter()) {
                *c += vi;
            }
        }
        for c in centroid.iter_mut() {
            *c /= n as f64;
        }
        let xr: Vec<f64> = centroid
            .iter()
            .zip(simplex[n].iter())
            .map(|(c, w)| c + alpha * (c - w))
            .collect();
        let fr = f(&xr);
        if fr < fvals[0] {
            let xe: Vec<f64> = centroid
                .iter()
                .zip(xr.iter())
                .map(|(c, r)| c + gamma * (r - c))
                .collect();
            let fe = f(&xe);
            if fe < fr {
                simplex[n] = xe;
                fvals[n] = fe;
            } else {
                simplex[n] = xr;
                fvals[n] = fr;
            }
        } else if fr < fvals[n - 1] {
            simplex[n] = xr;
            fvals[n] = fr;
        } else {
            let xc: Vec<f64> = centroid
                .iter()
                .zip(simplex[n].iter())
                .map(|(c, w)| c + rho * (w - c))
                .collect();
            let fc = f(&xc);
            if fc < fvals[n] {
                simplex[n] = xc;
                fvals[n] = fc;
            } else {
                let x0b = simplex[0].clone();
                for i in 1..=n {
                    simplex[i] = x0b
                        .iter()
                        .zip(simplex[i].iter())
                        .map(|(b, v)| b + sigma * (v - b))
                        .collect();
                    fvals[i] = f(&simplex[i]);
                }
            }
        }
    }
    let f_val = fvals[0];
    OptResult {
        x: simplex[0].clone(),
        f_val,
        n_iter,
        converged,
    }
}
/// Golden-section search for a unimodal minimum on `[a, b]`.
///
/// Returns `(x_min, f_min)`.
pub fn golden_section(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> (f64, f64) {
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut lo = a;
    let mut hi = b;
    let mut c = hi - phi * (hi - lo);
    let mut d = lo + phi * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    while (hi - lo).abs() > tol {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - phi * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + phi * (hi - lo);
            fd = f(d);
        }
    }
    let x_min = (lo + hi) / 2.0;
    (x_min, f(x_min))
}
/// Brent's method for 1-D minimization given a bracketing triple `(a, b, c)`
/// where `b` is between `a` and `c` and `f(b) < f(a)`, `f(b) < f(c)`.
///
/// Returns `(x_min, f_min)`.
pub fn brent(f: impl Fn(f64) -> f64, a: f64, b: f64, c: f64, tol: f64) -> (f64, f64) {
    let gold = 0.381_966_011_250_105;
    let tiny = 1e-10;
    let mut x = b;
    let mut w = b;
    let mut v = b;
    let mut fx = f(x);
    let mut fw = fx;
    let mut fv = fx;
    let mut lo = a.min(c);
    let mut hi = a.max(c);
    let mut d = 0.0f64;
    let mut e = 0.0f64;
    for _ in 0..500 {
        let xm = 0.5 * (lo + hi);
        let tol1 = tol * x.abs() + tiny;
        let tol2 = 2.0 * tol1;
        if (x - xm).abs() <= tol2 - 0.5 * (hi - lo) {
            return (x, fx);
        }
        let mut use_gold = true;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let q_val = (x - v) * (fx - fw);
            let p = (x - v) * q_val - (x - w) * r;
            let mut q2 = 2.0 * (q_val - r);
            let p = if q2 > 0.0 { -p } else { p };
            q2 = q2.abs();
            if p.abs() < (0.5 * q2 * e).abs() && p > q2 * (lo - x) && p < q2 * (hi - x) {
                e = d;
                d = p / q2;
                use_gold = false;
                let u = x + d;
                if (u - lo) < tol2 || (hi - u) < tol2 {
                    d = if xm >= x { tol1 } else { -tol1 };
                }
            }
        }
        if use_gold {
            e = if x >= xm { lo - x } else { hi - x };
            d = gold * e;
        }
        let u = x + if d.abs() >= tol1 {
            d
        } else if d >= 0.0 {
            tol1
        } else {
            -tol1
        };
        let fu = f(u);
        if fu <= fx {
            if u < x {
                hi = x;
            } else {
                lo = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                lo = u;
            } else {
                hi = u;
            }
            if fu <= fw || (w - x).abs() < tiny {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || (v - x).abs() < tiny || (v - w).abs() < tiny {
                v = u;
                fv = fu;
            }
        }
    }
    (x, fx)
}
/// Bisection method for root finding on `[a, b]`.
///
/// Requires `f(a)` and `f(b)` to have opposite signs.
/// Returns `Some(root)` when converged, or `None` if the sign condition is not
/// satisfied.
pub fn bisection(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> Option<f64> {
    let mut lo = a;
    let mut hi = b;
    let flo = f(lo);
    let fhi = f(hi);
    if flo * fhi > 0.0 {
        return None;
    }
    let mut flo_cur = flo;
    for _ in 0..1000 {
        if (hi - lo).abs() < tol {
            break;
        }
        let mid = (lo + hi) / 2.0;
        let fmid = f(mid);
        if fmid * flo_cur <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo_cur = fmid;
        }
    }
    Some((lo + hi) / 2.0)
}
/// Central-difference numerical gradient.
pub fn numerical_gradient(f: impl Fn(&[f64]) -> f64, x: &[f64], h: f64) -> Vec<f64> {
    let n = x.len();
    let mut grad = vec![0.0f64; n];
    let mut xp = x.to_vec();
    let mut xm = x.to_vec();
    for i in 0..n {
        xp[i] += h;
        xm[i] -= h;
        grad[i] = (f(&xp) - f(&xm)) / (2.0 * h);
        xp[i] = x[i];
        xm[i] = x[i];
    }
    grad
}
/// Finite-difference Hessian matrix (central differences, second order).
pub fn numerical_hessian(f: impl Fn(&[f64]) -> f64, x: &[f64], h: f64) -> Vec<Vec<f64>> {
    let n = x.len();
    let f0 = f(x);
    let mut hess = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            if i == j {
                let mut xp = x.to_vec();
                let mut xm = x.to_vec();
                xp[i] += h;
                xm[i] -= h;
                hess[i][j] = (f(&xp) - 2.0 * f0 + f(&xm)) / (h * h);
            } else {
                let mut xpp = x.to_vec();
                let mut xpm = x.to_vec();
                let mut xmp = x.to_vec();
                let mut xmm = x.to_vec();
                xpp[i] += h;
                xpp[j] += h;
                xpm[i] += h;
                xpm[j] -= h;
                xmp[i] -= h;
                xmp[j] += h;
                xmm[i] -= h;
                xmm[j] -= h;
                let val = (f(&xpp) - f(&xpm) - f(&xmp) + f(&xmm)) / (4.0 * h * h);
                hess[i][j] = val;
                hess[j][i] = val;
            }
        }
    }
    hess
}
/// Configuration for the Armijo backtracking line search.
#[derive(Debug, Clone, Copy)]
pub struct BacktrackingConfig {
    /// Initial step length.
    pub alpha0: f64,
    /// Shrinkage factor per rejection (`rho ∈ (0, 1)`).
    pub rho: f64,
    /// Armijo sufficient-decrease constant (`c ∈ (0, 1)`).
    pub c: f64,
}

impl Default for BacktrackingConfig {
    fn default() -> Self {
        Self {
            alpha0: 1.0,
            rho: 0.5,
            c: 1e-4,
        }
    }
}

/// Backtracking line search satisfying the Armijo (sufficient-decrease)
/// condition.
///
/// Returns the accepted step length `alpha`.
pub fn backtracking_line_search(
    f: impl Fn(&[f64]) -> f64,
    x: &[f64],
    direction: &[f64],
    f0: f64,
    grad_dot: f64,
    cfg: BacktrackingConfig,
) -> f64 {
    let BacktrackingConfig { alpha0, rho, c } = cfg;
    let mut alpha = alpha0;
    for _ in 0..100 {
        let x_new: Vec<f64> = x
            .iter()
            .zip(direction.iter())
            .map(|(xi, di)| xi + alpha * di)
            .collect();
        if f(&x_new) <= f0 + c * alpha * grad_dot {
            break;
        }
        alpha *= rho;
    }
    alpha
}
/// Box-constrained minimization via projected gradient descent.
///
/// Projects iterates onto `[lb[i\], ub[i]]` after each gradient step.
pub fn constrained_min_box(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lb: &[f64],
    ub: &[f64],
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    for i in 0..n {
        x[i] = x[i].clamp(lb[i], ub[i]);
    }
    let mut converged = false;
    let mut n_iter = 0u32;
    let mut lr = 1.0f64;
    for _ in 0..max_iter {
        n_iter += 1;
        let g = grad(&x);
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        let f0 = f(&x);
        let grad_dot = -g.iter().map(|gi| gi * gi).sum::<f64>();
        let dir: Vec<f64> = g.iter().map(|gi| -*gi).collect();
        let mut accepted = false;
        for _ in 0..30 {
            let mut x_try: Vec<f64> = x
                .iter()
                .zip(dir.iter())
                .map(|(xi, di)| xi + lr * di)
                .collect();
            for i in 0..n {
                x_try[i] = x_try[i].clamp(lb[i], ub[i]);
            }
            if f(&x_try) <= f0 + 1e-4 * lr * grad_dot {
                x = x_try;
                lr *= 1.2;
                lr = lr.min(1.0);
                accepted = true;
                break;
            }
            lr *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Gradient descent with heavy-ball (Polyak) momentum.
///
/// Update rule: `v ← momentum * v - lr * grad(x)`,  `x ← x + v`.
pub fn gradient_descent_momentum(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lr: f64,
    momentum: f64,
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    let mut vel = vec![0.0f64; n];
    let mut converged = false;
    let mut n_iter = 0u32;
    for _ in 0..max_iter {
        n_iter += 1;
        let g = grad(&x);
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        for i in 0..n {
            vel[i] = momentum * vel[i] - lr * g[i];
            x[i] += vel[i];
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Simulated annealing for global minimisation.
///
/// Uses a random walk proposal `x_new = x + N(0, step_size)`.  The Metropolis
/// acceptance criterion allows uphill moves with probability `exp(-ΔE / T)`.
///
/// Temperature follows a geometric schedule: `T ← T * cooling_rate`.
pub fn simulated_annealing(
    f: impl Fn(&[f64]) -> f64,
    x0: Vec<f64>,
    initial_temp: f64,
    cooling_rate: f64,
    step_size: f64,
    max_iter: u32,
) -> OptResult {
    use rand::RngExt;
    let mut x = x0.clone();
    let mut f_cur = f(&x);
    let mut best_x = x.clone();
    let mut best_f = f_cur;
    let mut temp = initial_temp;
    let mut rng = rand::rng();
    for _ in 0..max_iter {
        let x_new: Vec<f64> = x
            .iter()
            .map(|xi| xi + (rng.random::<f64>() - 0.5) * 2.0 * step_size)
            .collect();
        let f_new = f(&x_new);
        let delta = f_new - f_cur;
        if delta < 0.0 || (temp > 1e-300 && rng.random::<f64>() < (-delta / temp).exp()) {
            x = x_new;
            f_cur = f_new;
            if f_cur < best_f {
                best_f = f_cur;
                best_x = x.clone();
            }
        }
        temp *= cooling_rate;
    }
    OptResult {
        x: best_x,
        f_val: best_f,
        n_iter: max_iter,
        converged: temp < 1e-10,
    }
}
/// Configuration for Particle Swarm Optimization.
#[derive(Debug, Clone, Copy)]
pub struct PsoConfig {
    /// Number of particles.
    pub n_particles: usize,
    /// Inertia weight (0.5–0.9 typical).
    pub w: f64,
    /// Cognitive coefficient (≈ 2.0).
    pub c1: f64,
    /// Social coefficient (≈ 2.0).
    pub c2: f64,
    /// Maximum number of iterations.
    pub max_iter: u32,
}

impl Default for PsoConfig {
    fn default() -> Self {
        Self {
            n_particles: 30,
            w: 0.7,
            c1: 2.0,
            c2: 2.0,
            max_iter: 500,
        }
    }
}

/// Particle Swarm Optimization (PSO).
///
/// Minimises `f` using `cfg.n_particles` agents, each initialised
/// uniformly in `[lb[i], ub[i]]`.
pub fn particle_swarm(
    f: impl Fn(&[f64]) -> f64,
    lb: &[f64],
    ub: &[f64],
    cfg: PsoConfig,
) -> OptResult {
    let PsoConfig {
        n_particles,
        w,
        c1,
        c2,
        max_iter,
    } = cfg;
    let n = lb.len();
    let mut rng = rand::rng();
    let mut particles: Vec<Particle> = (0..n_particles)
        .map(|_| {
            let pos: Vec<f64> = (0..n)
                .map(|i| lb[i] + rng.random::<f64>() * (ub[i] - lb[i]))
                .collect();
            let vel: Vec<f64> = (0..n)
                .map(|i| (rng.random::<f64>() - 0.5) * (ub[i] - lb[i]) * 0.1)
                .collect();
            let val = f(&pos);
            Particle {
                best_pos: pos.clone(),
                best_val: val,
                pos,
                vel,
            }
        })
        .collect();
    let init_best = particles
        .iter()
        .min_by(|a, b| {
            a.best_val
                .partial_cmp(&b.best_val)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .expect("particles is non-empty");
    let mut global_best_pos = init_best.best_pos.clone();
    let mut global_best_val = init_best.best_val;
    for _ in 0..max_iter {
        for p in particles.iter_mut() {
            for i in 0..n {
                let r1: f64 = rng.random();
                let r2: f64 = rng.random();
                p.vel[i] = w * p.vel[i]
                    + c1 * r1 * (p.best_pos[i] - p.pos[i])
                    + c2 * r2 * (global_best_pos[i] - p.pos[i]);
                p.pos[i] = (p.pos[i] + p.vel[i]).clamp(lb[i], ub[i]);
            }
            let val = f(&p.pos);
            if val < p.best_val {
                p.best_val = val;
                p.best_pos = p.pos.clone();
                if val < global_best_val {
                    global_best_val = val;
                    global_best_pos = p.pos.clone();
                }
            }
        }
    }
    OptResult {
        x: global_best_pos,
        f_val: global_best_val,
        n_iter: max_iter,
        converged: true,
    }
}
/// Configuration for the genetic algorithm optimizer.
#[derive(Debug, Clone, Copy)]
pub struct GeneticAlgorithmConfig {
    /// Population size.
    pub pop_size: usize,
    /// Probability of crossover per individual.
    pub crossover_rate: f64,
    /// Probability of mutation per gene.
    pub mutation_rate: f64,
    /// Scale of random mutation perturbation.
    pub mutation_scale: f64,
    /// Number of generations.
    pub max_generations: u32,
}

impl Default for GeneticAlgorithmConfig {
    fn default() -> Self {
        Self {
            pop_size: 50,
            crossover_rate: 0.8,
            mutation_rate: 0.1,
            mutation_scale: 0.1,
            max_generations: 200,
        }
    }
}

/// Genetic algorithm for real-valued continuous minimisation.
///
/// Represents each individual as a `Vec<f64>` in `[lb[i], ub[i]]`.
/// Uses tournament selection, arithmetic crossover, and uniform mutation.
pub fn genetic_algorithm(
    f: impl Fn(&[f64]) -> f64,
    lb: &[f64],
    ub: &[f64],
    cfg: GeneticAlgorithmConfig,
) -> OptResult {
    let GeneticAlgorithmConfig {
        pop_size,
        crossover_rate,
        mutation_rate,
        mutation_scale,
        max_generations,
    } = cfg;
    let n = lb.len();
    let mut rng = rand::rng();
    let mut pop: Vec<Vec<f64>> = (0..pop_size)
        .map(|_| {
            (0..n)
                .map(|i| lb[i] + rng.random::<f64>() * (ub[i] - lb[i]))
                .collect()
        })
        .collect();
    let mut fitness: Vec<f64> = pop.iter().map(|x| f(x)).collect();
    for _ in 0..max_generations {
        let mut new_pop: Vec<Vec<f64>> = Vec::with_capacity(pop_size);
        for _ in 0..pop_size {
            let a_idx = {
                let i = (rng.random::<f64>() * pop_size as f64) as usize % pop_size;
                let j = (rng.random::<f64>() * pop_size as f64) as usize % pop_size;
                if fitness[i] < fitness[j] { i } else { j }
            };
            let b_idx = {
                let i = (rng.random::<f64>() * pop_size as f64) as usize % pop_size;
                let j = (rng.random::<f64>() * pop_size as f64) as usize % pop_size;
                if fitness[i] < fitness[j] { i } else { j }
            };
            let mut child = pop[a_idx].clone();
            if rng.random::<f64>() < crossover_rate {
                let alpha: f64 = rng.random();
                child = child
                    .iter()
                    .zip(pop[b_idx].iter())
                    .map(|(&a, &b)| alpha * a + (1.0 - alpha) * b)
                    .collect();
            }
            for i in 0..n {
                if rng.random::<f64>() < mutation_rate {
                    let delta = (rng.random::<f64>() - 0.5) * 2.0 * mutation_scale;
                    child[i] = (child[i] + delta).clamp(lb[i], ub[i]);
                }
            }
            new_pop.push(child);
        }
        pop = new_pop;
        fitness = pop.iter().map(|x| f(x)).collect();
    }
    let best_idx = fitness
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    OptResult {
        x: pop[best_idx].clone(),
        f_val: fitness[best_idx],
        n_iter: max_generations,
        converged: true,
    }
}
/// Configuration for L-BFGS-B (box-constrained L-BFGS).
#[derive(Debug, Clone, Copy)]
pub struct LbfgsbConfig {
    /// Number of correction vectors to store (history size).
    pub m: usize,
    /// Convergence tolerance on the gradient norm.
    pub tol: f64,
    /// Maximum number of iterations.
    pub max_iter: u32,
}

impl Default for LbfgsbConfig {
    fn default() -> Self {
        Self {
            m: 5,
            tol: 1e-8,
            max_iter: 1000,
        }
    }
}

/// L-BFGS-B: box-constrained variant of L-BFGS.
///
/// Enforces `lb[i] ≤ x[i] ≤ ub[i]` by projecting the search direction onto the
/// feasible set after each step (simple projected L-BFGS).
pub fn lbfgsb(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lb: &[f64],
    ub: &[f64],
    cfg: LbfgsbConfig,
) -> OptResult {
    let LbfgsbConfig { m, tol, max_iter } = cfg;
    let mut x: Vec<f64> = x0
        .iter()
        .enumerate()
        .map(|(i, &xi)| xi.clamp(lb[i], ub[i]))
        .collect();
    let mut g = grad(&x);
    let mut s_list: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut y_list: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut converged = false;
    let mut n_iter = 0u32;
    for _ in 0..max_iter {
        n_iter += 1;
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        let k = s_list.len();
        let mut q = g.clone();
        let mut alphas = vec![0.0f64; k];
        for i in (0..k).rev() {
            let sy = dot(&y_list[i], &s_list[i]);
            if sy.abs() < 1e-15 {
                continue;
            }
            let rho_i = 1.0 / sy;
            let a = rho_i * dot(&s_list[i], &q);
            alphas[i] = a;
            axpy(-a, &y_list[i].clone(), &mut q);
        }
        let mut r = if k > 0 {
            let scale = dot(&s_list[k - 1], &y_list[k - 1])
                / dot(&y_list[k - 1], &y_list[k - 1]).max(1e-15);
            q.iter().map(|qi| scale * qi).collect::<Vec<_>>()
        } else {
            q.clone()
        };
        for i in 0..k {
            let sy = dot(&y_list[i], &s_list[i]);
            if sy.abs() < 1e-15 {
                continue;
            }
            let rho_i = 1.0 / sy;
            let beta = rho_i * dot(&y_list[i], &r);
            axpy(alphas[i] - beta, &s_list[i].clone(), &mut r);
        }
        let direction: Vec<f64> = r.iter().map(|ri| -*ri).collect();
        let proj_dir: Vec<f64> = direction
            .iter()
            .enumerate()
            .map(|(i, &di)| {
                if (x[i] <= lb[i] && di < 0.0) || (x[i] >= ub[i] && di > 0.0) {
                    0.0
                } else {
                    di
                }
            })
            .collect();
        let f0 = f(&x);
        let gd = dot(&g, &proj_dir);
        let alpha =
            backtracking_line_search(&f, &x, &proj_dir, f0, gd, BacktrackingConfig::default());
        let s: Vec<f64> = proj_dir.iter().map(|di| alpha * di).collect();
        let x_new: Vec<f64> = x
            .iter()
            .zip(s.iter())
            .enumerate()
            .map(|(i, (xi, si))| (xi + si).clamp(lb[i], ub[i]))
            .collect();
        let g_new = grad(&x_new);
        let y: Vec<f64> = g_new
            .iter()
            .zip(g.iter())
            .map(|(gni, gi)| gni - gi)
            .collect();
        let sy = dot(&s, &y);
        if sy > 1e-15 {
            if s_list.len() == m {
                s_list.remove(0);
                y_list.remove(0);
            }
            s_list.push(s);
            y_list.push(y);
        }
        x = x_new;
        g = g_new;
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Nonlinear conjugate gradient minimization using the Fletcher-Reeves update.
///
/// Uses a backtracking line search at each step. The direction is restarted
/// every `n` iterations.
pub fn conjugate_gradient_minimize<F, G>(
    f: F,
    grad: G,
    x0: Vec<f64>,
    max_iter: usize,
    tol: f64,
) -> Vec<f64>
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let mut x = x0;
    let mut g = grad(&x);
    let mut d: Vec<f64> = g.iter().map(|gi| -*gi).collect();
    for iter in 0..max_iter {
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            break;
        }
        let f0 = f(&x);
        let gd = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum::<f64>();
        if gd >= 0.0 {
            d = g.iter().map(|gi| -*gi).collect();
        }
        let gd = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum::<f64>();
        let alpha = backtracking_line_search(&f, &x, &d, f0, gd, BacktrackingConfig::default());
        let x_new: Vec<f64> = (0..n).map(|i| x[i] + alpha * d[i]).collect();
        let g_new = grad(&x_new);
        let g_new_sq: f64 = g_new.iter().map(|v| v * v).sum();
        let g_sq: f64 = g.iter().map(|v| v * v).sum();
        let beta = if g_sq < 1e-30 || (iter + 1) % n == 0 {
            0.0
        } else {
            g_new_sq / g_sq
        };
        d = (0..n).map(|i| -g_new[i] + beta * d[i]).collect();
        x = x_new;
        g = g_new;
    }
    x
}
/// Strong Wolfe conditions line search.
///
/// Finds a step length `α` satisfying:
/// 1. Sufficient decrease (Armijo): `f(x + α d) ≤ f(x) + c1 * α * ∇f(x)·d`
/// 2. Curvature: `|∇f(x + α d)·d| ≤ c2 * |∇f(x)·d|`
///
/// Returns the accepted step length.
pub fn wolfe_line_search(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x: &[f64],
    direction: &[f64],
    f0: f64,
    g0_dot_d: f64,
    c1: f64,
    c2: f64,
    max_iter: usize,
) -> f64 {
    let n = x.len();
    let mut alpha = 1.0f64;
    let mut alpha_lo = 0.0f64;
    let mut alpha_hi = f64::INFINITY;
    for _ in 0..max_iter {
        let x_new: Vec<f64> = (0..n).map(|i| x[i] + alpha * direction[i]).collect();
        let f_new = f(&x_new);
        if f_new > f0 + c1 * alpha * g0_dot_d {
            alpha_hi = alpha;
            alpha = 0.5 * (alpha_lo + alpha_hi);
            continue;
        }
        let g_new = grad(&x_new);
        let g_new_dot_d: f64 = g_new.iter().zip(direction.iter()).map(|(g, d)| g * d).sum();
        if g_new_dot_d.abs() <= c2 * g0_dot_d.abs() {
            return alpha;
        }
        if g_new_dot_d > 0.0 {
            alpha_hi = alpha;
            alpha = 0.5 * (alpha_lo + alpha_hi);
        } else {
            alpha_lo = alpha;
            if alpha_hi.is_infinite() {
                alpha *= 2.0;
            } else {
                alpha = 0.5 * (alpha_lo + alpha_hi);
            }
        }
    }
    alpha
}
/// Augmented Lagrangian method for equality-constrained minimisation.
///
/// Minimises `f(x)` subject to `c_i(x) = 0` for each constraint in `constraints`.
///
/// Strategy: alternates between unconstrained minimisation of the augmented
/// Lagrangian `L_ρ(x,λ) = f(x) + Σ λ_i c_i(x) + (ρ/2) Σ c_i(x)²` and a
/// dual update `λ_i ← λ_i + ρ * c_i(x)`.
///
/// Each unconstrained sub-problem is solved with a few steps of gradient
/// descent (a full sub-solver can be plugged in via `inner_iter`).
pub fn augmented_lagrangian<F, G, C, DC>(
    f: F,
    grad_f: G,
    constraints: &[C],
    grad_constraints: &[DC],
    x0: Vec<f64>,
    rho0: f64,
    outer_iter: usize,
    inner_iter: usize,
    tol: f64,
) -> OptResult
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
    C: Fn(&[f64]) -> f64,
    DC: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let m = constraints.len();
    let mut x = x0;
    let mut lam = vec![0.0f64; m];
    let mut rho = rho0;
    let mut converged = false;
    let mut total_iter = 0u32;
    for _outer in 0..outer_iter {
        for _inner in 0..inner_iter {
            total_iter += 1;
            let mut g = grad_f(&x);
            for (i, (ci, dci)) in constraints.iter().zip(grad_constraints.iter()).enumerate() {
                let cv = ci(&x);
                let dv = dci(&x);
                let coeff = lam[i] + rho * cv;
                for j in 0..n {
                    g[j] += coeff * dv[j];
                }
            }
            let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if gnorm < tol * 0.1 {
                break;
            }
            let mut lr = 0.1;
            let aug_f = |xp: &[f64]| {
                let mut val = f(xp);
                for (i, ci) in constraints.iter().enumerate() {
                    let cv = ci(xp);
                    val += lam[i] * cv + 0.5 * rho * cv * cv;
                }
                val
            };
            let f0 = aug_f(&x);
            for _ in 0..20 {
                let x_try: Vec<f64> = (0..n).map(|j| x[j] - lr * g[j]).collect();
                if aug_f(&x_try) < f0 {
                    x = x_try;
                    break;
                }
                lr *= 0.5;
            }
        }
        let mut max_viol = 0.0f64;
        for (i, ci) in constraints.iter().enumerate() {
            let cv = ci(&x);
            lam[i] += rho * cv;
            max_viol = max_viol.max(cv.abs());
        }
        if max_viol < tol {
            converged = true;
            break;
        }
        rho *= 2.0;
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter: total_iter,
        converged,
    }
}
/// Cyclic coordinate descent (derivative-free).
///
/// Iterates over each coordinate in turn, performing a golden-section search
/// along that axis.  Useful when gradients are unavailable.
pub fn coordinate_descent(
    f: impl Fn(&[f64]) -> f64,
    x0: Vec<f64>,
    step: f64,
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    let mut converged = false;
    let mut n_iter = 0u32;
    for _iter in 0..max_iter {
        n_iter += 1;
        let x_old = x.clone();
        for d in 0..n {
            let (best_delta, _) = golden_section(
                |delta| {
                    let mut xp = x.clone();
                    xp[d] += delta;
                    f(&xp)
                },
                -step,
                step,
                tol * 1e-2,
            );
            x[d] += best_delta;
        }
        let diff: f64 = x
            .iter()
            .zip(x_old.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        if diff < tol {
            converged = true;
            break;
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Powell's method (derivative-free multidimensional minimisation).
///
/// Maintains a set of conjugate directions and performs line searches
/// along each.  Resets directions periodically to avoid linear dependence.
pub fn powell(f: impl Fn(&[f64]) -> f64, x0: Vec<f64>, tol: f64, max_iter: usize) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    let mut dirs: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut d = vec![0.0; n];
            d[i] = 1.0;
            d
        })
        .collect();
    let mut converged = false;
    let mut n_iter = 0u32;
    for _iter in 0..max_iter {
        n_iter += 1;
        let x_start = x.clone();
        let mut max_decrease = 0.0f64;
        let mut max_dir_idx = 0;
        for (k, dir) in dirs.iter().enumerate() {
            let f_before = f(&x);
            let (alpha, _) = golden_section(
                |a| {
                    let xp: Vec<f64> = x
                        .iter()
                        .zip(dir.iter())
                        .map(|(xi, di)| xi + a * di)
                        .collect();
                    f(&xp)
                },
                -10.0,
                10.0,
                tol * 0.01,
            );
            let xp: Vec<f64> = x
                .iter()
                .zip(dir.iter())
                .map(|(xi, di)| xi + alpha * di)
                .collect();
            let decrease = f_before - f(&xp);
            if decrease > max_decrease {
                max_decrease = decrease;
                max_dir_idx = k;
            }
            x = xp;
        }
        let diff: f64 = x
            .iter()
            .zip(x_start.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        if diff < tol {
            converged = true;
            break;
        }
        let new_dir: Vec<f64> = x.iter().zip(x_start.iter()).map(|(a, b)| a - b).collect();
        let norm: f64 = new_dir.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm > 1e-15 {
            dirs[max_dir_idx] = new_dir.iter().map(|v| v / norm).collect();
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// SGD with cosine annealing learning-rate schedule.
///
/// Useful for noisy / stochastic objective functions.  The gradient `stoch_grad`
/// may return a noisy (mini-batch) estimate.
pub fn sgd_cosine_annealing(
    f: impl Fn(&[f64]) -> f64,
    stoch_grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lr_max: f64,
    lr_min: f64,
    max_iter: u32,
) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    let mut n_iter = 0u32;
    for t in 0..max_iter {
        n_iter = t + 1;
        let cos_factor = 0.5 * (1.0 + (std::f64::consts::PI * t as f64 / max_iter as f64).cos());
        let lr = lr_min + (lr_max - lr_min) * cos_factor;
        let g = stoch_grad(&x);
        for i in 0..n {
            x[i] -= lr * g[i];
        }
    }
    let f_val = f(&x);
    let gnorm: f64 = stoch_grad(&x).iter().map(|v| v * v).sum::<f64>().sqrt();
    OptResult {
        x,
        f_val,
        n_iter,
        converged: gnorm < 1e-4,
    }
}
/// Full-memory BFGS (Broyden-Fletcher-Goldfarb-Shanno).
///
/// Maintains an approximate Hessian inverse `H` (n×n matrix) updated with
/// each iteration.  Falls back to gradient descent when curvature information
/// is insufficient.
pub fn bfgs(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    tol: f64,
    max_iter: usize,
) -> OptResult {
    let n = x0.len();
    let mut x = x0;
    let mut h: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0; n];
            row[i] = 1.0;
            row
        })
        .collect();
    let mut converged = false;
    let mut n_iter = 0u32;
    let mut g = grad(&x);
    for _ in 0..max_iter {
        n_iter += 1;
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        let p: Vec<f64> = (0..n)
            .map(|i| {
                -h[i]
                    .iter()
                    .zip(g.iter())
                    .map(|(hij, gj)| hij * gj)
                    .sum::<f64>()
            })
            .collect();
        let f0 = f(&x);
        let gd: f64 = g.iter().zip(p.iter()).map(|(gi, pi)| gi * pi).sum();
        let alpha = backtracking_line_search(&f, &x, &p, f0, gd, BacktrackingConfig::default());
        let x_new: Vec<f64> = (0..n).map(|i| x[i] + alpha * p[i]).collect();
        let g_new = grad(&x_new);
        let s: Vec<f64> = (0..n).map(|i| x_new[i] - x[i]).collect();
        let y: Vec<f64> = (0..n).map(|i| g_new[i] - g[i]).collect();
        let sy: f64 = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum();
        if sy > 1e-15 {
            let hy: Vec<f64> = (0..n)
                .map(|i| {
                    h[i].iter()
                        .zip(y.iter())
                        .map(|(hij, yj)| hij * yj)
                        .sum::<f64>()
                })
                .collect();
            let ythy: f64 = y.iter().zip(hy.iter()).map(|(yi, hyi)| yi * hyi).sum();
            let rho_bfgs = 1.0 / sy;
            for i in 0..n {
                for j in 0..n {
                    h[i][j] += rho_bfgs * s[i] * s[j] - (hy[i] * s[j] + s[i] * hy[j]) / ythy;
                }
            }
        }
        x = x_new;
        g = g_new;
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Clips the L2 norm of `grad` to `max_norm` (in-place).
///
/// Used in deep learning / physics-ML to stabilise training.
pub fn clip_gradient(grad: &mut [f64], max_norm: f64) {
    let norm: f64 = grad.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm > max_norm && norm > 1e-15 {
        let scale = max_norm / norm;
        for g in grad.iter_mut() {
            *g *= scale;
        }
    }
}
/// Gradient descent with gradient clipping.
pub fn gradient_descent_clipped(
    f: impl Fn(&[f64]) -> f64,
    grad: impl Fn(&[f64]) -> Vec<f64>,
    x0: Vec<f64>,
    lr: f64,
    max_norm: f64,
    tol: f64,
    max_iter: u32,
) -> OptResult {
    let mut x = x0;
    let mut converged = false;
    let mut n_iter = 0u32;
    for _ in 0..max_iter {
        n_iter += 1;
        let mut g = grad(&x);
        clip_gradient(&mut g, max_norm);
        let gnorm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if gnorm < tol {
            converged = true;
            break;
        }
        for (xi, gi) in x.iter_mut().zip(g.iter()) {
            *xi -= lr * gi;
        }
    }
    let f_val = f(&x);
    OptResult {
        x,
        f_val,
        n_iter,
        converged,
    }
}
/// Compute the lower-triangular Cholesky factor L of an n×n matrix stored
/// row-major.  Returns L such that L * L^T ≈ A.
///
/// Falls back to the identity if the matrix is not positive-definite.
pub fn cholesky_lower(a: &[f64], n: usize) -> Vec<f64> {
    let mut l = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let s: f64 = (0..j).map(|k| l[i * n + k] * l[j * n + k]).sum();
            if i == j {
                let d = a[i * n + i] - s;
                l[i * n + j] = if d > 0.0 { d.sqrt() } else { 1e-8 };
            } else {
                let ljj = l[j * n + j];
                l[i * n + j] = if ljj.abs() > f64::EPSILON {
                    (a[i * n + j] - s) / ljj
                } else {
                    0.0
                };
            }
        }
    }
    l
}
/// Forward substitution: solve L * x = b for lower-triangular L (n×n, row-major).
pub fn forward_solve_lower(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut x = vec![0.0_f64; n];
    for i in 0..n {
        let s: f64 = (0..i).map(|j| l[i * n + j] * x[j]).sum();
        let lii = l[i * n + i];
        x[i] = if lii.abs() > f64::EPSILON {
            (b[i] - s) / lii
        } else {
            0.0
        };
    }
    x
}
