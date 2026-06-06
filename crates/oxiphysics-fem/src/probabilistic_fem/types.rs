//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Alpha-cut representation of a fuzzy number.
#[derive(Clone, Debug)]
pub struct AlphaCut {
    /// Alpha level (membership degree) in \[0, 1\].
    pub alpha: f64,
    /// Lower bound of the alpha-cut interval.
    pub lower: f64,
    /// Upper bound of the alpha-cut interval.
    pub upper: f64,
}
/// Design variable for RBDO.
#[derive(Clone, Debug)]
pub struct DesignVar {
    /// Variable name.
    pub name: String,
    /// Current value.
    pub value: f64,
    /// Lower bound.
    pub lb: f64,
    /// Upper bound.
    pub ub: f64,
}
/// Stochastic FEM with spatially random Young's modulus.
///
/// Models a 1D bar with `n_elem` elements.  Young's modulus is treated as a
/// random field discretised by `n_kl` KL modes.  The deterministic system
/// `K(E) u = F` is solved for each stochastic realisation.
pub struct StochasticFem {
    /// Number of elements.
    pub n_elem: usize,
    /// Bar cross-section area.
    pub area: f64,
    /// Bar element length.
    pub elem_len: f64,
    /// Applied tip load.
    pub tip_load: f64,
    /// Random field for Young's modulus.
    pub young_field: RandomField,
    /// Mean Young's modulus.
    pub e_mean: f64,
    /// Coefficient of variation of E.
    pub e_cov: f64,
}
impl StochasticFem {
    /// Construct a stochastic bar FEM.
    ///
    /// # Arguments
    /// * `n_elem` — number of bar elements.
    /// * `area` — cross-sectional area.
    /// * `total_length` — total bar length.
    /// * `tip_load` — applied axial force at free end.
    /// * `e_mean` — mean Young's modulus.
    /// * `e_cov` — coefficient of variation (σ_E / μ_E).
    /// * `n_kl` — number of KL modes.
    pub fn new(
        n_elem: usize,
        area: f64,
        total_length: f64,
        tip_load: f64,
        e_mean: f64,
        e_cov: f64,
        n_kl: usize,
    ) -> Self {
        let elem_len = total_length / n_elem as f64;
        let sigma2 = (e_mean * e_cov).powi(2);
        let length_scale = total_length * 0.3;
        let young_field = RandomField::new(n_elem, total_length, n_kl, sigma2, length_scale, 0);
        Self {
            n_elem,
            area,
            elem_len,
            tip_load,
            young_field,
            e_mean,
            e_cov,
        }
    }
    /// Assemble and solve the bar for given element Young's moduli.
    ///
    /// Returns displacement at each node (DOF).
    pub fn solve(&self, e_vec: &[f64]) -> Vec<f64> {
        let ndof = self.n_elem + 1;
        let mut k_global = vec![vec![0.0f64; ndof]; ndof];
        let mut f_global = vec![0.0f64; ndof];
        for ie in 0..self.n_elem {
            let ei = if ie < e_vec.len() {
                e_vec[ie].max(1e-3)
            } else {
                self.e_mean
            };
            let ke = ei * self.area / self.elem_len;
            let dofs = [ie, ie + 1];
            let k_elem = [[ke, -ke], [-ke, ke]];
            for (a, &da) in dofs.iter().enumerate() {
                for (b, &db) in dofs.iter().enumerate() {
                    k_global[da][db] += k_elem[a][b];
                }
            }
        }
        f_global[self.n_elem] = self.tip_load;
        let penalty = 1e15;
        k_global[0][0] += penalty;
        f_global[0] = 0.0;
        solve_chol(&k_global, &f_global)
    }
    /// Solve for a given random sample (xi vector for KL expansion).
    pub fn solve_sample(&self, xi: &[f64]) -> Vec<f64> {
        let e_sample = self
            .young_field
            .sample(xi)
            .iter()
            .map(|&de| (self.e_mean + de).max(1e-3))
            .collect::<Vec<f64>>();
        self.solve(&e_sample)
    }
    /// Deterministic solution with mean E.
    pub fn solve_mean(&self) -> Vec<f64> {
        let e_vec = vec![self.e_mean; self.n_elem];
        self.solve(&e_vec)
    }
    /// Analytical tip displacement for uniform bar: u = F L / (E A).
    pub fn analytical_tip_displacement(&self) -> f64 {
        let l = self.elem_len * self.n_elem as f64;
        self.tip_load * l / (self.e_mean * self.area)
    }
}
/// Monte Carlo FEM for uncertainty propagation.
///
/// Runs `n_samples` FEM solves with randomly drawn material properties and
/// collects statistics on the response quantity of interest.
pub struct MonteCarloFem {
    /// Underlying stochastic FEM model.
    pub sfem: StochasticFem,
    /// Number of Monte Carlo samples.
    pub n_samples: usize,
    /// Collected tip-displacement samples.
    pub samples: Vec<f64>,
    /// Current random number state (LCG).
    pub(super) lcg_state: u64,
}
impl MonteCarloFem {
    /// Construct a Monte Carlo FEM runner.
    pub fn new(sfem: StochasticFem, n_samples: usize) -> Self {
        Self {
            sfem,
            n_samples,
            samples: Vec::new(),
            lcg_state: 12345678901,
        }
    }
    /// Generate next pseudo-random U(0,1) value (LCG).
    fn next_u01(&mut self) -> f64 {
        self.lcg_state = self
            .lcg_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.lcg_state >> 33) as f64 / (u32::MAX as f64)
    }
    /// Generate next standard normal via Box-Muller.
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_u01().max(1e-15);
        let u2 = self.next_u01().max(1e-15);
        box_muller(u1, u2)
    }
    /// Run all Monte Carlo samples.
    pub fn run(&mut self) {
        self.samples.clear();
        let n_kl = self.sfem.young_field.n_terms;
        for _ in 0..self.n_samples {
            let xi: Vec<f64> = (0..n_kl).map(|_| self.next_normal()).collect();
            let u = self.sfem.solve_sample(&xi);
            let tip = *u.last().unwrap_or(&0.0);
            self.samples.push(tip);
        }
    }
    /// Mean of tip displacement across samples.
    pub fn mean_tip(&self) -> f64 {
        mean(&self.samples)
    }
    /// Standard deviation of tip displacement.
    pub fn std_tip(&self) -> f64 {
        variance(&self.samples).sqrt()
    }
    /// Coefficient of variation.
    pub fn cov_tip(&self) -> f64 {
        cov_stat(&self.samples)
    }
    /// Empirical failure probability for `u_tip > u_limit`.
    pub fn failure_probability(&self, u_limit: f64) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let n_fail = self.samples.iter().filter(|&&u| u > u_limit).count();
        n_fail as f64 / self.samples.len() as f64
    }
    /// 95th percentile of tip displacement (sorted estimate).
    pub fn percentile_95(&self) -> f64 {
        if self.samples.is_empty() {
            return 0.0;
        }
        let mut sorted = self.samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = ((sorted.len() as f64 * 0.95) as usize).min(sorted.len() - 1);
        sorted[idx]
    }
}
/// Polynomial chaos expansion (PCE) coefficient.
#[derive(Clone, Debug)]
pub struct PceCoefficient {
    /// Multi-index α (one exponent per random variable).
    pub alpha: Vec<usize>,
    /// Coefficient value.
    pub value: f64,
}
/// Spectral stochastic FEM using polynomial chaos expansion.
///
/// Expands the solution in orthogonal polynomials (Hermite for Gaussian):
/// `u(x,ω) = Σ_α u_α(x) Ψ_α(ξ)`
/// and uses Galerkin projection to compute deterministic coefficients.
pub struct SpectralStochasticFem {
    /// Number of random dimensions (KL modes).
    pub n_dim: usize,
    /// PCE order (maximum polynomial degree).
    pub pce_order: usize,
    /// Number of elements.
    pub n_elem: usize,
    /// Bar area.
    pub area: f64,
    /// Element length.
    pub elem_len: f64,
    /// Applied load.
    pub tip_load: f64,
    /// Mean Young's modulus.
    pub e_mean: f64,
    /// PCE coefficients for tip displacement.
    pub u_tip_pce: Vec<PceCoefficient>,
}
impl SpectralStochasticFem {
    /// Construct a spectral stochastic FEM.
    pub fn new(
        n_dim: usize,
        pce_order: usize,
        n_elem: usize,
        area: f64,
        total_length: f64,
        tip_load: f64,
        e_mean: f64,
    ) -> Self {
        let elem_len = total_length / n_elem as f64;
        Self {
            n_dim,
            pce_order,
            n_elem,
            area,
            elem_len,
            tip_load,
            e_mean,
            u_tip_pce: Vec::new(),
        }
    }
    /// Probabilist's Hermite polynomial Heₙ(x).
    pub fn hermite(n: usize, x: f64) -> f64 {
        match n {
            0 => 1.0,
            1 => x,
            _ => {
                let mut h_prev = 1.0;
                let mut h_curr = x;
                for k in 1..n {
                    let h_next = x * h_curr - k as f64 * h_prev;
                    h_prev = h_curr;
                    h_curr = h_next;
                }
                h_curr
            }
        }
    }
    /// Factorial of n.
    fn factorial(n: usize) -> f64 {
        (1..=n).fold(1.0, |acc, k| acc * k as f64)
    }
    /// Norm squared of Hermite basis Ψ_α: E\[Ψ_α²\] = α!
    pub fn pce_norm_sq(alpha: &[usize]) -> f64 {
        alpha.iter().map(|&a| Self::factorial(a)).product()
    }
    /// Evaluate basis polynomial Ψ_α at point ξ.
    pub fn eval_basis(alpha: &[usize], xi: &[f64]) -> f64 {
        alpha
            .iter()
            .zip(xi.iter())
            .map(|(&a, &x)| Self::hermite(a, x))
            .product()
    }
    /// Generate all multi-indices of total degree ≤ order with n_dim variables.
    pub fn multi_indices(&self) -> Vec<Vec<usize>> {
        let mut indices = Vec::new();
        self.gen_indices(&mut vec![0; self.n_dim], 0, &mut indices);
        indices
    }
    fn gen_indices(&self, current: &mut Vec<usize>, dim: usize, result: &mut Vec<Vec<usize>>) {
        if dim == self.n_dim {
            if current.iter().sum::<usize>() <= self.pce_order {
                result.push(current.clone());
            }
            return;
        }
        for deg in 0..=self.pce_order {
            current[dim] = deg;
            if current[..=dim].iter().sum::<usize>() <= self.pce_order {
                self.gen_indices(current, dim + 1, result);
            }
        }
        current[dim] = 0;
    }
    /// Compute PCE coefficients by non-intrusive collocation (tensor quadrature).
    ///
    /// Uses 3-point Gauss-Hermite quadrature per dimension.
    pub fn compute_pce_coefficients(&mut self, e_perturbation_cov: f64) {
        let gh_pts = [-1.7320508, 0.0, 1.7320508_f64];
        let gh_wts = [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0];
        let alphas = self.multi_indices();
        let n_pce = alphas.len();
        let mut u_tip_coeff = vec![0.0f64; n_pce];
        let sfem = StochasticFem::new(
            self.n_elem,
            self.area,
            self.elem_len * self.n_elem as f64,
            self.tip_load,
            self.e_mean,
            e_perturbation_cov,
            self.n_dim.min(4),
        );
        for (qi, &xi_q) in gh_pts.iter().enumerate() {
            let wi = gh_wts[qi];
            let xi_vec = vec![xi_q; self.n_dim];
            let u = sfem.solve_sample(&xi_vec);
            let u_tip = *u.last().unwrap_or(&0.0);
            for (k, alpha) in alphas.iter().enumerate() {
                let psi = Self::eval_basis(alpha, &xi_vec);
                u_tip_coeff[k] += wi * u_tip * psi;
            }
        }
        self.u_tip_pce = alphas
            .iter()
            .zip(u_tip_coeff.iter())
            .map(|(alpha, &val)| {
                let norm2_pce = Self::pce_norm_sq(alpha);
                PceCoefficient {
                    alpha: alpha.clone(),
                    value: val / norm2_pce.max(1e-300),
                }
            })
            .collect();
    }
    /// Mean of tip displacement from PCE (zeroth-order coefficient).
    pub fn pce_mean(&self) -> f64 {
        self.u_tip_pce
            .iter()
            .find(|c| c.alpha.iter().all(|&a| a == 0))
            .map(|c| c.value)
            .unwrap_or(0.0)
    }
    /// Variance of tip displacement from PCE: Var = Σ_{α≠0} uₐ² α!
    pub fn pce_variance(&self) -> f64 {
        self.u_tip_pce
            .iter()
            .filter(|c| c.alpha.iter().any(|&a| a > 0))
            .map(|c| c.value * c.value * Self::pce_norm_sq(&c.alpha))
            .sum()
    }
}
/// A single failure mode with its limit-state characteristics.
#[derive(Clone, Debug)]
pub struct FailureMode {
    /// Name of this failure mode.
    pub name: String,
    /// Reliability index β_HL for this mode.
    pub beta: f64,
    /// Direction cosines (sensitivity factors αᵢ).
    pub alpha: Vec<f64>,
    /// Probability of failure for this mode.
    pub pf: f64,
}
impl FailureMode {
    /// Construct a failure mode from its reliability index and sensitivity.
    pub fn new(name: &str, beta: f64, alpha: Vec<f64>) -> Self {
        let pf = phi_cdf(-beta);
        Self {
            name: name.to_string(),
            beta,
            alpha,
            pf,
        }
    }
    /// Importance measure for random variable i.
    pub fn importance(&self, i: usize) -> f64 {
        if i < self.alpha.len() {
            self.alpha[i].powi(2)
        } else {
            0.0
        }
    }
}
/// Reliability-based design optimization solver.
///
/// Minimises an objective function subject to probabilistic constraints
/// `P(g(X) ≤ 0) ≥ Φ(-β_target)` using a decoupled RBDO approach.
pub struct RbdoSolver {
    /// Design variables.
    pub design_vars: Vec<DesignVar>,
    /// Target reliability index β.
    pub beta_target: f64,
    /// Maximum RBDO iterations.
    pub max_iter: usize,
    /// Reliability index history.
    pub beta_history: Vec<f64>,
    /// Objective function history.
    pub obj_history: Vec<f64>,
    /// Internal FEM model.
    pub(super) sfem: StochasticFem,
}
impl RbdoSolver {
    /// Construct an RBDO solver for a bar problem.
    pub fn new(
        sfem: StochasticFem,
        beta_target: f64,
        max_iter: usize,
        design_vars: Vec<DesignVar>,
    ) -> Self {
        Self {
            design_vars,
            beta_target,
            max_iter,
            beta_history: Vec::new(),
            obj_history: Vec::new(),
            sfem,
        }
    }
    /// Evaluate limit-state function g(u_tip) = u_limit - u_tip.
    fn limit_state(&self, u_tip: f64, u_limit: f64) -> f64 {
        u_limit - u_tip
    }
    /// Compute FORM reliability index β_HL by linearising limit state.
    ///
    /// Uses the mean-value (MV-FORM) approximation.
    pub fn form_beta(&self, u_limit: f64) -> f64 {
        let u0 = self.sfem.solve_mean();
        let u_tip_mean = *u0.last().unwrap_or(&0.0);
        let pf_fem = PerturbationFem::new(
            self.sfem.n_elem,
            self.sfem.area,
            self.sfem.elem_len * self.sfem.n_elem as f64,
            self.sfem.tip_load,
            self.sfem.e_mean,
            (self.sfem.e_mean * self.sfem.e_cov).powi(2),
            self.sfem.elem_len * 3.0,
        );
        let var_tip = pf_fem.first_order_variance_tip();
        let std_tip = var_tip.sqrt().max(1e-14);
        let g_mean = self.limit_state(u_tip_mean, u_limit);
        g_mean / std_tip
    }
    /// Run RBDO outer loop: update area to satisfy reliability constraint.
    pub fn run(&mut self, u_limit: f64) {
        self.beta_history.clear();
        self.obj_history.clear();
        let mut area = self.sfem.area;
        for _ in 0..self.max_iter {
            let sfem_iter = StochasticFem::new(
                self.sfem.n_elem,
                area,
                self.sfem.elem_len * self.sfem.n_elem as f64,
                self.sfem.tip_load,
                self.sfem.e_mean,
                self.sfem.e_cov,
                self.sfem.young_field.n_terms,
            );
            let u0 = sfem_iter.solve_mean();
            let u_tip = *u0.last().unwrap_or(&0.0);
            let pf_fem = PerturbationFem::new(
                self.sfem.n_elem,
                area,
                self.sfem.elem_len * self.sfem.n_elem as f64,
                self.sfem.tip_load,
                self.sfem.e_mean,
                (self.sfem.e_mean * self.sfem.e_cov).powi(2),
                self.sfem.elem_len * 3.0,
            );
            let var_tip = pf_fem.first_order_variance_tip();
            let std_tip = var_tip.sqrt().max(1e-14);
            let g_mean = u_limit - u_tip;
            let beta = g_mean / std_tip;
            self.beta_history.push(beta);
            self.obj_history.push(area);
            if beta >= self.beta_target {
                break;
            }
            area *= 1.05;
            if let Some(dv) = self.design_vars.first_mut() {
                dv.value = area;
            }
        }
        self.sfem.area = area;
    }
}
/// Probabilistic failure mode analysis for a structural system.
///
/// Combines multiple failure modes using series/parallel system bounds.
pub struct FailureModeAnalysis {
    /// List of failure modes.
    pub modes: Vec<FailureMode>,
    /// Correlation matrix between failure modes.
    pub correlation: Vec<Vec<f64>>,
}
impl FailureModeAnalysis {
    /// Construct a failure mode analysis.
    pub fn new(modes: Vec<FailureMode>) -> Self {
        let n = modes.len();
        let corr = vec![vec![1.0f64; n]; n];
        Self {
            modes,
            correlation: corr,
        }
    }
    /// Set correlation between modes i and j.
    pub fn set_correlation(&mut self, i: usize, j: usize, rho: f64) {
        if i < self.correlation.len() && j < self.correlation[i].len() {
            self.correlation[i][j] = rho;
            self.correlation[j][i] = rho;
        }
    }
    /// Compute correlations from sensitivity vectors.
    pub fn compute_correlations_from_alpha(&mut self) {
        let n = self.modes.len();
        for i in 0..n {
            for j in 0..n {
                let rho = dot(&self.modes[i].alpha, &self.modes[j].alpha);
                let ni = norm2(&self.modes[i].alpha);
                let nj = norm2(&self.modes[j].alpha);
                let denom = (ni * nj).max(1e-14);
                self.correlation[i][j] = rho / denom;
            }
        }
    }
    /// Series system failure probability upper bound (union bound).
    pub fn series_pf_upper(&self) -> f64 {
        self.modes.iter().map(|m| m.pf).sum::<f64>().min(1.0)
    }
    /// Series system failure probability lower bound (Ditlevsen).
    pub fn series_pf_lower(&self) -> f64 {
        if self.modes.is_empty() {
            return 0.0;
        }
        let mut pf = self.modes[0].pf;
        for i in 1..self.modes.len() {
            let pf_i = self.modes[i].pf;
            let rho = self.correlation[0][i];
            let correction =
                phi_cdf(-self.modes[0].beta) * phi_cdf(-self.modes[i].beta) * (1.0 + rho);
            pf += pf_i - correction;
        }
        pf.clamp(0.0, 1.0)
    }
    /// Parallel system failure probability (product for independent modes).
    pub fn parallel_pf(&self) -> f64 {
        self.modes.iter().map(|m| m.pf).product::<f64>()
    }
    /// Most likely failure mode (smallest beta).
    pub fn most_likely_mode(&self) -> Option<&FailureMode> {
        self.modes.iter().min_by(|a, b| {
            a.beta
                .partial_cmp(&b.beta)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
    /// System reliability index from series approximation.
    pub fn system_beta_series(&self) -> f64 {
        let pf = self.series_pf_upper();
        -phi_inv(pf.clamp(1e-15, 1.0 - 1e-15))
    }
    /// Sensitivity rank of random variables across all modes (combined α²).
    pub fn global_sensitivity(&self) -> Vec<f64> {
        if self.modes.is_empty() {
            return vec![];
        }
        let n_rv = self.modes[0].alpha.len();
        let mut gs = vec![0.0f64; n_rv];
        for m in &self.modes {
            for (i, &a) in m.alpha.iter().enumerate() {
                gs[i] += a * a;
            }
        }
        let total = gs.iter().sum::<f64>().max(1e-14);
        gs.iter_mut().for_each(|v| *v /= total);
        gs
    }
}
/// Karhunen-Loève random field representation.
///
/// Discretises a second-order random field on a 1D mesh using
/// the Fredholm integral eigenvalue problem.  KL expansion:
/// `H(x,ω) = μ(x) + Σ_i √λᵢ φᵢ(x) ξᵢ(ω)`
/// where ξᵢ are uncorrelated standard normal random variables.
pub struct RandomField {
    /// Mesh points for field discretization.
    pub mesh: Vec<f64>,
    /// Mean function values at mesh points.
    pub mean: Vec<f64>,
    /// KL eigenvalues λᵢ (energy of each mode).
    pub eigenvalues: Vec<f64>,
    /// KL eigenvectors φᵢ (shape of each mode, columns).
    pub eigenvectors: Vec<Vec<f64>>,
    /// Number of KL terms retained.
    pub n_terms: usize,
    /// Covariance kernel type: 0=exponential, 1=sqexp, 2=matern52.
    pub kernel_type: usize,
    /// Variance σ².
    pub sigma2: f64,
    /// Correlation length ℓ.
    pub length_scale: f64,
}
impl RandomField {
    /// Create a random field on a uniform 1D mesh.
    ///
    /// # Arguments
    /// * `n_points` — number of mesh points.
    /// * `domain_length` — physical length of domain.
    /// * `n_terms` — number of KL modes to retain.
    /// * `sigma2` — variance of the field.
    /// * `length_scale` — correlation length.
    /// * `kernel_type` — 0=exponential, 1=Gaussian, 2=Matérn5/2.
    pub fn new(
        n_points: usize,
        domain_length: f64,
        n_terms: usize,
        sigma2: f64,
        length_scale: f64,
        kernel_type: usize,
    ) -> Self {
        let dx = if n_points > 1 {
            domain_length / (n_points - 1) as f64
        } else {
            1.0
        };
        let mesh: Vec<f64> = (0..n_points).map(|i| i as f64 * dx).collect();
        let mean = vec![0.0f64; n_points];
        let n_terms_clamped = n_terms.min(n_points);
        let mut rf = Self {
            mesh,
            mean,
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            n_terms: n_terms_clamped,
            kernel_type,
            sigma2,
            length_scale,
        };
        rf.compute_kl_expansion();
        rf
    }
    /// Evaluate covariance kernel between two points.
    pub fn kernel(&self, r: f64) -> f64 {
        match self.kernel_type {
            0 => cov_exp(r, self.sigma2, self.length_scale),
            1 => cov_sqexp(r, self.sigma2, self.length_scale),
            2 => cov_matern52(r, self.sigma2, self.length_scale),
            _ => cov_exp(r, self.sigma2, self.length_scale),
        }
    }
    /// Build the covariance matrix C\[i\]\[j\] = C(|xᵢ - xⱼ|).
    pub fn covariance_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.mesh.len();
        let mut c = vec![vec![0.0f64; n]; n];
        for (i, row) in c.iter_mut().enumerate().take(n) {
            for (j, cell) in row.iter_mut().enumerate().take(n) {
                let r = (self.mesh[i] - self.mesh[j]).abs();
                *cell = self.kernel(r);
            }
        }
        c
    }
    /// Compute KL eigenvalues/vectors via successive power iterations.
    fn compute_kl_expansion(&mut self) {
        let c = self.covariance_matrix();
        let n = c.len();
        if n == 0 {
            return;
        }
        let mut residual = c.clone();
        self.eigenvalues.clear();
        self.eigenvectors.clear();
        for _mode in 0..self.n_terms {
            let (lambda, phi) = power_iteration(&residual, 200);
            if lambda < 1e-14 {
                break;
            }
            self.eigenvalues.push(lambda);
            self.eigenvectors.push(phi.clone());
            for i in 0..n {
                for j in 0..n {
                    residual[i][j] -= lambda * phi[i] * phi[j];
                }
            }
        }
    }
    /// Sample a realization using the provided ξ vector (length ≥ n_terms).
    pub fn sample(&self, xi: &[f64]) -> Vec<f64> {
        let n = self.mesh.len();
        let mut field = self.mean.clone();
        for (k, (&lam, phi)) in self
            .eigenvalues
            .iter()
            .zip(self.eigenvectors.iter())
            .enumerate()
        {
            let xik = if k < xi.len() { xi[k] } else { 0.0 };
            let scale = lam.sqrt() * xik;
            for i in 0..n {
                field[i] += scale * phi[i];
            }
        }
        field
    }
    /// Variance of the KL approximation at point index `i`.
    pub fn variance_at(&self, i: usize) -> f64 {
        self.eigenvalues
            .iter()
            .zip(self.eigenvectors.iter())
            .map(|(&lam, phi)| {
                let pi = if i < phi.len() { phi[i] } else { 0.0 };
                lam * pi * pi
            })
            .sum()
    }
    /// Relative energy captured by retained modes (Σλᵢ / Σλ_total).
    pub fn energy_ratio(&self) -> f64 {
        let total = self
            .covariance_matrix()
            .iter()
            .enumerate()
            .map(|(i, row)| row[i])
            .sum::<f64>();
        if total < 1e-14 {
            return 1.0;
        }
        self.eigenvalues.iter().sum::<f64>() / total
    }
}
/// Fuzzy FEM propagating triangular fuzzy numbers through the FEM.
///
/// Uses the alpha-cut decomposition: for each alpha level, bounds the
/// response via an interval computation.
pub struct FuzzyFem {
    /// Number of elements.
    pub n_elem: usize,
    /// Bar area.
    pub area: f64,
    /// Element length.
    pub elem_len: f64,
    /// Applied load.
    pub tip_load: f64,
    /// Alpha-cuts of Young's modulus (triangular fuzzy number).
    pub e_fuzzy_cuts: Vec<AlphaCut>,
    /// Output alpha-cuts of tip displacement.
    pub u_tip_cuts: Vec<AlphaCut>,
}
impl FuzzyFem {
    /// Construct a fuzzy FEM with triangular fuzzy Young's modulus.
    ///
    /// Triangular fuzzy: (e_min, e_mode, e_max).
    pub fn new(
        n_elem: usize,
        area: f64,
        total_length: f64,
        tip_load: f64,
        e_min: f64,
        e_mode: f64,
        e_max: f64,
        n_alpha: usize,
    ) -> Self {
        let elem_len = total_length / n_elem as f64;
        let mut e_fuzzy_cuts = Vec::with_capacity(n_alpha);
        for k in 0..n_alpha {
            let alpha = k as f64 / (n_alpha - 1).max(1) as f64;
            let lo = e_min + alpha * (e_mode - e_min);
            let hi = e_max - alpha * (e_max - e_mode);
            e_fuzzy_cuts.push(AlphaCut {
                alpha,
                lower: lo,
                upper: hi,
            });
        }
        Self {
            n_elem,
            area,
            elem_len,
            tip_load,
            e_fuzzy_cuts,
            u_tip_cuts: Vec::new(),
        }
    }
    /// Solve for a uniform Young's modulus E.
    fn solve_uniform_e(&self, e_val: f64) -> f64 {
        let sfem = StochasticFem::new(
            self.n_elem,
            self.area,
            self.elem_len * self.n_elem as f64,
            self.tip_load,
            e_val,
            0.0,
            1,
        );
        let u = sfem.solve_mean();
        *u.last().unwrap_or(&0.0)
    }
    /// Propagate fuzzy E through FEM using alpha-cut method.
    pub fn compute_fuzzy_response(&mut self) {
        self.u_tip_cuts.clear();
        for cut in &self.e_fuzzy_cuts {
            let u_lo = self.solve_uniform_e(cut.upper.max(1.0));
            let u_hi = self.solve_uniform_e(cut.lower.max(1.0));
            self.u_tip_cuts.push(AlphaCut {
                alpha: cut.alpha,
                lower: u_lo.min(u_hi),
                upper: u_lo.max(u_hi),
            });
        }
    }
    /// Get the crisp (alpha=1) tip displacement estimate.
    pub fn crisp_tip(&self) -> Option<f64> {
        self.u_tip_cuts
            .iter()
            .find(|c| (c.alpha - 1.0).abs() < 1e-10)
            .map(|c| 0.5 * (c.lower + c.upper))
    }
    /// Core membership value (alpha=1) upper bound.
    pub fn core_upper(&self) -> Option<f64> {
        self.u_tip_cuts
            .iter()
            .find(|c| (c.alpha - 1.0).abs() < 1e-10)
            .map(|c| c.upper)
    }
}
/// Interval FEM solving `K([E]) u = F` with interval Young's modulus.
///
/// Uses the vertex method: evaluates response at all combinations of
/// E_i ∈ {E_i^lo, E_i^hi}.
pub struct IntervalFem {
    /// Number of elements.
    pub n_elem: usize,
    /// Bar area.
    pub area: f64,
    /// Element length.
    pub elem_len: f64,
    /// Applied load.
    pub tip_load: f64,
    /// Interval Young's modulus for each element.
    pub e_interval: Vec<Interval>,
    /// Resulting interval tip displacement.
    pub u_tip_interval: Interval,
}
impl IntervalFem {
    /// Construct an interval FEM.
    pub fn new(
        n_elem: usize,
        area: f64,
        total_length: f64,
        tip_load: f64,
        e_lo: f64,
        e_hi: f64,
    ) -> Self {
        let elem_len = total_length / n_elem as f64;
        let e_interval = vec![Interval::new(e_lo, e_hi); n_elem];
        Self {
            n_elem,
            area,
            elem_len,
            tip_load,
            e_interval,
            u_tip_interval: Interval::new(0.0, 0.0),
        }
    }
    /// Vertex method: enumerate corners of the interval hypercube.
    ///
    /// For `n_elem` up to 20 (2^n_elem evaluations).
    pub fn vertex_method(&mut self) {
        let n = self.n_elem.min(20);
        let n_vertices: usize = 1 << n;
        let mut u_min = f64::INFINITY;
        let mut u_max = f64::NEG_INFINITY;
        for v in 0..n_vertices {
            let e_vec: Vec<f64> = (0..n)
                .map(|i| {
                    if (v >> i) & 1 == 0 {
                        self.e_interval[i].lo
                    } else {
                        self.e_interval[i].hi
                    }
                })
                .collect();
            let e_avg = e_vec.iter().sum::<f64>() / e_vec.len() as f64;
            let u = {
                let total_length = self.elem_len * self.n_elem as f64;
                self.tip_load * total_length / (e_avg.max(1e-3) * self.area)
            };
            u_min = u_min.min(u);
            u_max = u_max.max(u);
        }
        self.u_tip_interval = Interval::new(u_min, u_max);
    }
    /// Midpoint solution (deterministic with mean E).
    pub fn midpoint_solution(&self) -> f64 {
        let e_mid = self.e_interval[0].mid();
        let total_length = self.elem_len * self.n_elem as f64;
        self.tip_load * total_length / (e_mid.max(1e-3) * self.area)
    }
}
/// Interval representation \[lo, hi\].
#[derive(Clone, Copy, Debug)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}
impl Interval {
    /// Construct an interval.
    pub fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }
    /// Midpoint.
    pub fn mid(&self) -> f64 {
        0.5 * (self.lo + self.hi)
    }
    /// Half-width (radius).
    pub fn rad(&self) -> f64 {
        0.5 * (self.hi - self.lo)
    }
    /// Interval scalar multiply.
    pub fn scale(self, s: f64) -> Interval {
        if s >= 0.0 {
            Interval::new(self.lo * s, self.hi * s)
        } else {
            Interval::new(self.hi * s, self.lo * s)
        }
    }
}

impl std::ops::Add for Interval {
    type Output = Interval;
    /// Interval addition: [a, b] + [c, d] = [a+c, b+d].
    fn add(self, other: Interval) -> Interval {
        Interval::new(self.lo + other.lo, self.hi + other.hi)
    }
}

impl std::ops::Mul for Interval {
    type Output = Interval;
    /// Interval multiplication: [a, b] * [c, d] = [min, max] of all products.
    fn mul(self, other: Interval) -> Interval {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        Interval::new(
            products.iter().cloned().fold(f64::INFINITY, f64::min),
            products.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}
/// Direct differentiation method for response sensitivities.
///
/// Computes `du/dθ` analytically: `K du/dθ = dF/dθ - dK/dθ u`.
/// This avoids finite difference perturbations.
pub struct DdmSensitivity {
    /// Stochastic FEM model.
    pub sfem: StochasticFem,
}
impl DdmSensitivity {
    /// Construct a DDM sensitivity analysis object.
    pub fn new(sfem: StochasticFem) -> Self {
        Self { sfem }
    }
    /// Compute `dK/dA` (sensitivity of stiffness w.r.t. cross-section area).
    pub fn dk_da(&self) -> Vec<Vec<f64>> {
        let ndof = self.sfem.n_elem + 1;
        let mut k = vec![vec![0.0f64; ndof]; ndof];
        let ke = self.sfem.e_mean / self.sfem.elem_len;
        for ie in 0..self.sfem.n_elem {
            k[ie][ie] += ke;
            k[ie][ie + 1] -= ke;
            k[ie + 1][ie] -= ke;
            k[ie + 1][ie + 1] += ke;
        }
        k
    }
    /// Compute `dK/dE_mean` (sensitivity of stiffness w.r.t. mean Young's modulus).
    pub fn dk_de(&self) -> Vec<Vec<f64>> {
        let ndof = self.sfem.n_elem + 1;
        let mut k = vec![vec![0.0f64; ndof]; ndof];
        let ke = self.sfem.area / self.sfem.elem_len;
        for ie in 0..self.sfem.n_elem {
            k[ie][ie] += ke;
            k[ie][ie + 1] -= ke;
            k[ie + 1][ie] -= ke;
            k[ie + 1][ie + 1] += ke;
        }
        k
    }
    /// DDM sensitivity `du/dA` at mean state.
    pub fn sensitivity_da(&self) -> Vec<f64> {
        let u0 = self.sfem.solve_mean();
        let k0_full = {
            let pfem = PerturbationFem::new(
                self.sfem.n_elem,
                self.sfem.area,
                self.sfem.elem_len * self.sfem.n_elem as f64,
                self.sfem.tip_load,
                self.sfem.e_mean,
                0.0,
                1.0,
            );
            pfem.mean_stiffness()
        };
        let dk = self.dk_da();
        let ndof = self.sfem.n_elem + 1;
        let mut rhs = vec![0.0f64; ndof];
        for i in 0..ndof {
            for j in 0..ndof {
                rhs[i] -= dk[i][j] * if j < u0.len() { u0[j] } else { 0.0 };
            }
        }
        solve_chol(&k0_full, &rhs)
    }
    /// DDM sensitivity `du/dE_mean`.
    pub fn sensitivity_de(&self) -> Vec<f64> {
        let u0 = self.sfem.solve_mean();
        let k0_full = {
            let pfem = PerturbationFem::new(
                self.sfem.n_elem,
                self.sfem.area,
                self.sfem.elem_len * self.sfem.n_elem as f64,
                self.sfem.tip_load,
                self.sfem.e_mean,
                0.0,
                1.0,
            );
            pfem.mean_stiffness()
        };
        let dk = self.dk_de();
        let ndof = self.sfem.n_elem + 1;
        let mut rhs = vec![0.0f64; ndof];
        for i in 0..ndof {
            for j in 0..ndof {
                rhs[i] -= dk[i][j] * if j < u0.len() { u0[j] } else { 0.0 };
            }
        }
        solve_chol(&k0_full, &rhs)
    }
    /// Adjoint sensitivity: λᵀ dK/dθ u (efficient for many parameters).
    pub fn adjoint_sensitivity_e(&self, qoi_index: usize) -> f64 {
        let u0 = self.sfem.solve_mean();
        let k0 = {
            let pfem = PerturbationFem::new(
                self.sfem.n_elem,
                self.sfem.area,
                self.sfem.elem_len * self.sfem.n_elem as f64,
                self.sfem.tip_load,
                self.sfem.e_mean,
                0.0,
                1.0,
            );
            pfem.mean_stiffness()
        };
        let mut rhs_adj = vec![0.0f64; self.sfem.n_elem + 1];
        if qoi_index < rhs_adj.len() {
            rhs_adj[qoi_index] = 1.0;
        }
        let lam = solve_chol(&k0, &rhs_adj);
        let dk = self.dk_de();
        let ndof = dk.len();
        let mut dkdu = vec![0.0f64; ndof];
        for i in 0..ndof {
            for j in 0..ndof {
                dkdu[i] += dk[i][j] * if j < u0.len() { u0[j] } else { 0.0 };
            }
        }
        -dot(&lam, &dkdu)
    }
}
/// First- and second-order perturbation-based stochastic FEM.
///
/// Expands the stiffness matrix K = K₀ + ε K₁ + ε² K₂ + …
/// and solves for the perturbed response to obtain mean and variance
/// without Monte Carlo sampling.
pub struct PerturbationFem {
    /// Number of elements.
    pub n_elem: usize,
    /// Bar area.
    pub area: f64,
    /// Element length.
    pub elem_len: f64,
    /// Applied load.
    pub tip_load: f64,
    /// Mean Young's modulus.
    pub e_mean: f64,
    /// Variance of Young's modulus.
    pub e_var: f64,
    /// Correlation length for exponential covariance.
    pub corr_length: f64,
}
impl PerturbationFem {
    /// Construct a perturbation FEM model.
    pub fn new(
        n_elem: usize,
        area: f64,
        total_length: f64,
        tip_load: f64,
        e_mean: f64,
        e_var: f64,
        corr_length: f64,
    ) -> Self {
        let elem_len = total_length / n_elem as f64;
        Self {
            n_elem,
            area,
            elem_len,
            tip_load,
            e_mean,
            e_var,
            corr_length,
        }
    }
    /// Assemble mean global stiffness matrix K₀.
    pub fn mean_stiffness(&self) -> Vec<Vec<f64>> {
        let ndof = self.n_elem + 1;
        let mut k = vec![vec![0.0f64; ndof]; ndof];
        let ke = self.e_mean * self.area / self.elem_len;
        for ie in 0..self.n_elem {
            k[ie][ie] += ke;
            k[ie][ie + 1] -= ke;
            k[ie + 1][ie] -= ke;
            k[ie + 1][ie + 1] += ke;
        }
        k[0][0] += 1e15;
        k
    }
    /// Assemble first-order stiffness perturbation for element `ie`.
    pub fn perturbed_stiffness(&self, ie: usize) -> Vec<Vec<f64>> {
        let ndof = self.n_elem + 1;
        let mut k = vec![vec![0.0f64; ndof]; ndof];
        let ke = self.area / self.elem_len;
        if ie < self.n_elem {
            k[ie][ie] += ke;
            k[ie][ie + 1] -= ke;
            k[ie + 1][ie] -= ke;
            k[ie + 1][ie + 1] += ke;
        }
        k
    }
    /// Zeroth-order response (mean solution).
    pub fn zeroth_order_response(&self) -> Vec<f64> {
        let k0 = self.mean_stiffness();
        let mut f = vec![0.0f64; self.n_elem + 1];
        f[self.n_elem] = self.tip_load;
        solve_chol(&k0, &f)
    }
    /// First-order response for random variable `r` (element ie).
    pub fn first_order_response(&self, ie: usize, u0: &[f64]) -> Vec<f64> {
        let k0 = self.mean_stiffness();
        let k1 = self.perturbed_stiffness(ie);
        let ndof = self.n_elem + 1;
        let mut rhs = vec![0.0f64; ndof];
        for i in 0..ndof {
            for j in 0..ndof {
                rhs[i] -= k1[i][j] * if j < u0.len() { u0[j] } else { 0.0 };
            }
        }
        solve_chol(&k0, &rhs)
    }
    /// First-order mean of tip displacement.
    pub fn first_order_mean_tip(&self) -> f64 {
        let u0 = self.zeroth_order_response();
        *u0.last().unwrap_or(&0.0)
    }
    /// First-order variance of tip displacement.
    ///
    /// Var\[u_tip\] ≈ Σᵢ Σⱼ Cov(Eᵢ,Eⱼ) · (∂u_tip/∂Eᵢ) · (∂u_tip/∂Eⱼ)
    pub fn first_order_variance_tip(&self) -> f64 {
        let u0 = self.zeroth_order_response();
        let n = self.n_elem;
        let mut var = 0.0;
        let mut sensitivities = vec![0.0f64; n];
        for (ie, sens_ie) in sensitivities.iter_mut().enumerate().take(n) {
            let u1 = self.first_order_response(ie, &u0);
            *sens_ie = *u1.last().unwrap_or(&0.0);
        }
        for i in 0..n {
            for j in 0..n {
                let xi = (i as f64 + 0.5) * self.elem_len;
                let xj = (j as f64 + 0.5) * self.elem_len;
                let cov_ij = cov_exp((xi - xj).abs(), self.e_var, self.corr_length);
                var += cov_ij * sensitivities[i] * sensitivities[j];
            }
        }
        var
    }
    /// Second-order correction to mean (leading curvature term).
    pub fn second_order_mean_correction(&self) -> f64 {
        let var = self.first_order_variance_tip();
        let k2_scale = self.e_var / (self.e_mean * self.e_mean).max(1e-14);
        0.5 * var * k2_scale
    }
}
