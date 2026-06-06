//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{D2Q9_CX, D2Q9_CY, D2Q9_W, fick_flux};

/// Multi-component LBM with cross-diffusion via Fick's law.
///
/// Each component has its own distribution function array and relaxation time.
/// Cross-diffusion between components is handled via explicit source terms.
#[derive(Debug, Clone)]
pub struct MultiComponentLbm {
    /// Number of grid nodes in x.
    pub nx: usize,
    /// Number of grid nodes in y.
    pub ny: usize,
    /// Number of species components.
    pub ncomp: usize,
    /// Scalar fields φ_α for each species (length `ncomp × nx × ny`).
    pub phi: Vec<Vec<f64>>,
    /// Distribution functions g_α (length `ncomp × nx × ny × 9`).
    pub g: Vec<Vec<f64>>,
    /// Relaxation times τ_α for each species.
    pub tau: Vec<f64>,
    /// Cross-diffusion matrix D_αβ (length `ncomp × ncomp`), stored row-major.
    pub cross_diff: Vec<f64>,
}
impl MultiComponentLbm {
    /// Create a new `MultiComponentLbm` with uniform initial concentrations.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – Grid dimensions.
    /// * `tau`        – Relaxation times for each species (length `ncomp`).
    /// * `phi0`       – Initial concentration for each species (length `ncomp`).
    pub fn new(nx: usize, ny: usize, tau: Vec<f64>, phi0: Vec<f64>) -> Self {
        let ncomp = tau.len();
        assert_eq!(phi0.len(), ncomp, "phi0 length must equal ncomp");
        let n = nx * ny;
        let mut phi = Vec::with_capacity(ncomp);
        let mut g = Vec::with_capacity(ncomp);
        for &phi0_a in phi0.iter() {
            phi.push(vec![phi0_a; n]);
            let mut g_alpha = vec![0.0_f64; n * 9];
            for k in 0..n {
                for (q, &w) in D2Q9_W.iter().enumerate() {
                    g_alpha[k * 9 + q] = w * phi0_a;
                }
            }
            g.push(g_alpha);
        }
        let cross_diff = vec![0.0_f64; ncomp * ncomp];
        Self {
            nx,
            ny,
            ncomp,
            phi,
            g,
            tau,
            cross_diff,
        }
    }
    /// Set the cross-diffusion coefficient between species α and β.
    pub fn set_cross_diff(&mut self, alpha: usize, beta: usize, d_ab: f64) {
        self.cross_diff[alpha * self.ncomp + beta] = d_ab;
    }
    /// Get the cross-diffusion coefficient D_αβ.
    pub fn cross_diff_val(&self, alpha: usize, beta: usize) -> f64 {
        self.cross_diff[alpha * self.ncomp + beta]
    }
    /// Perform one BGK collision for all species without cross-diffusion.
    pub fn collide_uncoupled(&mut self) {
        let n = self.nx * self.ny;
        for alpha in 0..self.ncomp {
            let tau_a = self.tau[alpha];
            for k in 0..n {
                let phi_k = self.phi[alpha][k];
                for (q, &w) in D2Q9_W.iter().enumerate() {
                    let feq = w * phi_k;
                    self.g[alpha][k * 9 + q] -= (self.g[alpha][k * 9 + q] - feq) / tau_a;
                }
                self.phi[alpha][k] = self.g[alpha][k * 9..k * 9 + 9].iter().sum();
            }
        }
    }
    /// Total concentration of species `alpha`.
    pub fn total_concentration(&self, alpha: usize) -> f64 {
        self.phi[alpha].iter().sum()
    }
    /// Sum of all species concentrations at node `k` (total mixture concentration).
    pub fn mixture_concentration(&self, k: usize) -> f64 {
        (0..self.ncomp).map(|a| self.phi[a][k]).sum()
    }
}
/// Turbulent mixing model based on k–ε turbulent diffusivity.
///
/// Computes the turbulent diffusivity `D_t = C_mu * k² / ε` and combines
/// it with a molecular Schmidt number to give the effective species diffusivity.
#[derive(Debug, Clone)]
pub struct TurbulentMixing {
    /// Turbulent Schmidt number Sc_t (dimensionless, typically ~0.7–0.9).
    pub schmidt_number: f64,
    /// Molecular (laminar) diffusivity (m² s⁻¹).
    pub molecular_diffusivity: f64,
    /// C_μ constant (standard value 0.09).
    pub c_mu: f64,
}
impl TurbulentMixing {
    /// Create a new turbulent-mixing model.
    ///
    /// # Arguments
    /// * `schmidt_number`        – Turbulent Schmidt number.
    /// * `molecular_diffusivity` – Laminar diffusivity (m² s⁻¹).
    pub fn new(schmidt_number: f64, molecular_diffusivity: f64) -> Self {
        Self {
            schmidt_number,
            molecular_diffusivity,
            c_mu: 0.09,
        }
    }
    /// Compute turbulent (eddy) diffusivity from k–ε variables.
    ///
    /// `D_t = C_μ k² / ε`
    ///
    /// Returns 0 when ε ≤ 0 or k ≤ 0.
    ///
    /// # Arguments
    /// * `k`       – Turbulent kinetic energy (m² s⁻²).
    /// * `epsilon` – Turbulent dissipation rate (m² s⁻³).
    pub fn turbulent_diffusivity(&self, k: f64, epsilon: f64) -> f64 {
        if k <= 0.0 || epsilon <= 0.0 {
            return 0.0;
        }
        self.c_mu * k * k / epsilon
    }
    /// Effective diffusivity = molecular + turbulent/Sc_t.
    ///
    /// # Arguments
    /// * `k`       – Turbulent kinetic energy (m² s⁻²).
    /// * `epsilon` – Turbulent dissipation rate (m² s⁻³).
    pub fn effective_diffusivity(&self, k: f64, epsilon: f64) -> f64 {
        let d_t = self.turbulent_diffusivity(k, epsilon);
        let sc = self.schmidt_number.max(1e-30);
        self.molecular_diffusivity + d_t / sc
    }
}
/// Advection-diffusion of a passive scalar on a D2Q9 grid.
///
/// The passive scalar `phi` does not feed back to the flow field.
/// A BGK-LBM step evolves the scalar distribution functions `g`.
#[derive(Debug, Clone)]
pub struct PassiveScalarField {
    /// Grid width (nodes in x).
    pub nx: usize,
    /// Grid height (nodes in y).
    pub ny: usize,
    /// Scalar field φ (length `nx * ny`).
    pub phi: Vec<f64>,
    /// Distribution functions for scalar (length `nx * ny * 9`).
    pub g: Vec<f64>,
    /// BGK relaxation time for the scalar.
    pub tau: f64,
    /// Velocity field ux (length `nx * ny`).
    pub ux: Vec<f64>,
    /// Velocity field uy (length `nx * ny`).
    pub uy: Vec<f64>,
}
impl PassiveScalarField {
    /// Create a new `PassiveScalarField` with uniform scalar `phi0`.
    pub fn new(nx: usize, ny: usize, tau: f64, phi0: f64) -> Self {
        let n = nx * ny;
        let mut g = vec![0.0_f64; n * 9];
        for k in 0..n {
            for q in 0..9 {
                g[k * 9 + q] = D2Q9_W[q] * phi0;
            }
        }
        Self {
            nx,
            ny,
            phi: vec![phi0; n],
            g,
            tau,
            ux: vec![0.0_f64; n],
            uy: vec![0.0_f64; n],
        }
    }
    /// Set the velocity field (flat, row-major).
    pub fn set_velocity(&mut self, ux: Vec<f64>, uy: Vec<f64>) {
        self.ux = ux;
        self.uy = uy;
    }
    /// BGK equilibrium for the scalar distribution function.
    #[inline]
    fn geq(&self, phi_k: f64, ux_k: f64, uy_k: f64, q: usize) -> f64 {
        let cu = D2Q9_CX[q] * ux_k + D2Q9_CY[q] * uy_k;
        let uu = ux_k * ux_k + uy_k * uy_k;
        D2Q9_W[q] * phi_k * (1.0_f64 + 3.0_f64 * cu + 4.5_f64 * cu * cu - 1.5_f64 * uu)
    }
    /// Perform one BGK collision step (no streaming) for the scalar.
    pub fn collide(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let phi_k = self.phi[k];
            let ux_k = self.ux[k];
            let uy_k = self.uy[k];
            for q in 0..9 {
                let feq = self.geq(phi_k, ux_k, uy_k, q);
                self.g[k * 9 + q] -= (self.g[k * 9 + q] - feq) / self.tau;
            }
            self.phi[k] = self.g[k * 9..k * 9 + 9].iter().sum();
        }
    }
    /// Total scalar (zeroth moment of g, summed over all nodes).
    pub fn total_scalar(&self) -> f64 {
        self.phi.iter().sum()
    }
    /// Spatial mean of the scalar field.
    pub fn mean_scalar(&self) -> f64 {
        if self.phi.is_empty() {
            return 0.0_f64;
        }
        self.total_scalar() / self.phi.len() as f64
    }
    /// Spatial variance of the scalar field.
    pub fn variance_scalar(&self) -> f64 {
        let mu = self.mean_scalar();
        self.phi.iter().map(|&p| (p - mu) * (p - mu)).sum::<f64>() / self.phi.len() as f64
    }
}
/// Two-component mixture with an inter-diffusion coefficient.
///
/// Models Fick's first law for a binary system A–B.
#[derive(Debug, Clone)]
pub struct BinaryMixture {
    /// Component A.
    pub component_a: MixingComponent,
    /// Component B.
    pub component_b: MixingComponent,
    /// Inter-diffusion coefficient D_AB (m² s⁻¹).
    pub interdiffusion_coeff: f64,
}
impl BinaryMixture {
    /// Create a new binary mixture.
    pub fn new(
        component_a: MixingComponent,
        component_b: MixingComponent,
        interdiffusion_coeff: f64,
    ) -> Self {
        Self {
            component_a,
            component_b,
            interdiffusion_coeff,
        }
    }
    /// Compute the Fickian inter-diffusion flux for a given concentration gradient.
    ///
    /// `J = -D_AB * grad_c`
    ///
    /// # Arguments
    /// * `grad_c` – Concentration gradient of component A (mol m⁻⁴).
    pub fn compute_flux(&self, grad_c: f64) -> f64 {
        fick_flux(self.interdiffusion_coeff, grad_c)
    }
    /// Mole fraction of component A at grid index `idx`.
    ///
    /// Returns 0 when both densities are zero.
    pub fn mole_fraction_a(&self, idx: usize) -> f64 {
        let ca = *self.component_a.density_field.get(idx).unwrap_or(&0.0)
            / self.component_a.molecular_weight.max(1e-30);
        let cb = *self.component_b.density_field.get(idx).unwrap_or(&0.0)
            / self.component_b.molecular_weight.max(1e-30);
        let tot = ca + cb;
        if tot <= 0.0 { 0.0 } else { ca / tot }
    }
}
/// Analysis tools for mixing quality: segregation index, efficiency, Lyapunov.
#[derive(Debug, Clone)]
pub struct MixingAnalysis {
    /// History of variance values (for Lyapunov exponent estimate).
    pub variance_history: Vec<f64>,
    /// History of time stamps.
    pub time_history: Vec<f64>,
}
impl MixingAnalysis {
    /// Create a new `MixingAnalysis` with empty history.
    pub fn new() -> Self {
        Self {
            variance_history: Vec::new(),
            time_history: Vec::new(),
        }
    }
    /// Record the current variance at time `t`.
    pub fn record(&mut self, t: f64, var: f64) {
        self.time_history.push(t);
        self.variance_history.push(var);
    }
    /// Segregation index I_s = σ²(t) / σ²(0).
    ///
    /// Returns 1 at t=0, approaches 0 as mixing proceeds.
    /// Returns `None` if no initial variance recorded or initial variance is zero.
    pub fn segregation_index(&self) -> Option<f64> {
        if self.variance_history.len() < 2 {
            return None;
        }
        let var0 = self.variance_history[0];
        if var0 < f64::EPSILON {
            return None;
        }
        let var_now = self.variance_history[self.variance_history.len() - 1];
        Some(var_now / var0)
    }
    /// Mixing efficiency η = 1 − I_s.
    pub fn mixing_efficiency(&self) -> Option<f64> {
        self.segregation_index().map(|is| 1.0_f64 - is)
    }
    /// Lyapunov exponent estimate from variance decay: λ = -d(ln σ²)/dt.
    ///
    /// Uses linear regression of ln(σ²) vs t. Returns `None` if history
    /// has fewer than 2 points or variance is non-positive.
    pub fn lyapunov_exponent_estimate(&self) -> Option<f64> {
        let n = self.variance_history.len();
        if n < 2 {
            return None;
        }
        let pts: Vec<(f64, f64)> = self
            .time_history
            .iter()
            .zip(self.variance_history.iter())
            .filter(|&(_, &v)| v > 0.0_f64)
            .map(|(&t, &v)| (t, v.ln()))
            .collect();
        if pts.len() < 2 {
            return None;
        }
        let m = pts.len() as f64;
        let sum_t: f64 = pts.iter().map(|(t, _)| t).sum();
        let sum_lv: f64 = pts.iter().map(|(_, lv)| lv).sum();
        let sum_t2: f64 = pts.iter().map(|(t, _)| t * t).sum();
        let sum_tlv: f64 = pts.iter().map(|(t, lv)| t * lv).sum();
        let denom = m * sum_t2 - sum_t * sum_t;
        if denom.abs() < f64::EPSILON {
            return None;
        }
        let slope = (m * sum_tlv - sum_t * sum_lv) / denom;
        Some(-slope)
    }
    /// Return the number of recorded data points.
    pub fn num_records(&self) -> usize {
        self.time_history.len()
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// Taylor dispersion: effective axial dispersion in Poiseuille flow.
///
/// In a tube of radius R with mean velocity U, Taylor (1953) showed that
/// the effective axial diffusivity is:
///
/// `D_eff = D + U² R² / (48 D)`   (for a 2D channel of half-width H)
///
/// In a circular tube: `D_eff = D + U² R² / (48 D)`.
/// In a 2D channel:    `D_eff = D + U² H² / (210 D)`.
#[derive(Debug, Clone)]
pub struct TaylorDispersion {
    /// Molecular diffusivity D (m² s⁻¹).
    pub diffusivity: f64,
    /// Mean velocity U (m s⁻¹).
    pub mean_velocity: f64,
    /// Channel half-width or tube radius R (m).
    pub length_scale: f64,
    /// Geometry: `true` = circular tube, `false` = 2D channel.
    pub is_tube: bool,
}
impl TaylorDispersion {
    /// Create a new `TaylorDispersion` model.
    ///
    /// # Arguments
    /// * `diffusivity`  – Molecular diffusivity (m² s⁻¹).
    /// * `mean_velocity`– Mean axial velocity (m s⁻¹).
    /// * `length_scale` – Tube radius or channel half-width (m).
    /// * `is_tube`      – `true` for a circular tube, `false` for a 2D channel.
    pub fn new(diffusivity: f64, mean_velocity: f64, length_scale: f64, is_tube: bool) -> Self {
        Self {
            diffusivity,
            mean_velocity,
            length_scale,
            is_tube,
        }
    }
    /// Effective axial diffusivity D_eff.
    ///
    /// Circular tube:  `D_eff = D + U² R² / (48 D)`
    /// 2D channel:     `D_eff = D + U² H² / (210 D)`
    ///
    /// Returns `D` when D ≤ 0 to avoid division by zero.
    pub fn effective_diffusivity(&self) -> f64 {
        let d = self.diffusivity;
        if d <= 0.0_f64 {
            return d;
        }
        let u2 = self.mean_velocity * self.mean_velocity;
        let r2 = self.length_scale * self.length_scale;
        if self.is_tube {
            d + u2 * r2 / (48.0_f64 * d)
        } else {
            d + u2 * r2 / (210.0_f64 * d)
        }
    }
    /// Axial dispersion coefficient (effective − molecular).
    pub fn axial_dispersion(&self) -> f64 {
        (self.effective_diffusivity() - self.diffusivity).max(0.0_f64)
    }
    /// Péclet number Pe = U R / D.
    pub fn peclet(&self) -> f64 {
        if self.diffusivity <= 0.0_f64 {
            return f64::INFINITY;
        }
        self.mean_velocity.abs() * self.length_scale / self.diffusivity
    }
    /// Taylor time scale τ_T = R² / D.
    pub fn taylor_time_scale(&self) -> f64 {
        if self.diffusivity <= 0.0_f64 {
            return f64::INFINITY;
        }
        self.length_scale * self.length_scale / self.diffusivity
    }
    /// 1D Gaussian concentration profile at position x and time t.
    ///
    /// `c(x, t) = c0 / sqrt(4π D_eff t) * exp(-(x - U t)² / (4 D_eff t))`
    ///
    /// Returns 0 when t ≤ 0.
    pub fn gaussian_profile(&self, x: f64, t: f64, c0: f64) -> f64 {
        if t <= 0.0_f64 {
            return 0.0_f64;
        }
        let d_eff = self.effective_diffusivity();
        let x_centered = x - self.mean_velocity * t;
        let var = 4.0_f64 * d_eff * t;
        c0 / (std::f64::consts::PI * var).sqrt() * (-x_centered * x_centered / var).exp()
    }
}
/// 1-D advection-diffusion solver for a single species.
///
/// Implements a first-order upwind advection plus central-difference diffusion
/// step, with an optional linear source/sink term.
#[derive(Debug, Clone)]
pub struct SpeciesTransport {
    /// Diffusion coefficient (m² s⁻¹).
    pub diffusivity: f64,
    /// Linear source/sink coefficient (s⁻¹).
    /// The source rate is `source_term * c`.
    pub source_term: f64,
}
impl SpeciesTransport {
    /// Create a new species-transport solver.
    pub fn new(diffusivity: f64, source_term: f64) -> Self {
        Self {
            diffusivity,
            source_term,
        }
    }
    /// Perform one explicit advection-diffusion-source step.
    ///
    /// Uses upwind advection and central-difference diffusion (1-D, Δx = `dx`).
    ///
    /// # Arguments
    /// * `u`  – Advection velocity (m s⁻¹).  Positive = left-to-right.
    /// * `c`  – Current concentration at node (mol m⁻³).
    /// * `dt` – Time step (s).
    /// * `dx` – Grid spacing (m).
    ///
    /// Returns the updated concentration.
    pub fn compute_step(&self, u: f64, c: f64, dt: f64, dx: f64) -> f64 {
        let d = self.diffusivity;
        let diffusion = -2.0 * d * dt / (dx * dx) * c;
        let advection = -u * dt / dx * c;
        let source = self.source_term * c * dt;
        (c + diffusion + advection + source).max(0.0)
    }
}
/// Advection-diffusion LBM solver tracking a scalar concentration field.
///
/// The scalar field `c_scalar` holds local concentrations.  Two distribution
/// function arrays (`f_c` for the scalar, `f_u` for the velocity) are
/// stored; each uses a D2Q9 stencil (9 populations per node).
#[derive(Debug, Clone)]
pub struct MixingLBM {
    /// Number of nodes in x.
    pub nx: usize,
    /// Number of nodes in y.
    pub ny: usize,
    /// Scalar concentration field, length nx*ny.
    pub c_scalar: Vec<f64>,
    /// Scalar distribution functions, length nx*ny*9.
    pub f_c: Vec<f64>,
    /// Velocity distribution functions, length nx*ny*9.
    pub f_u: Vec<f64>,
    /// Péclet number (advection-to-diffusion ratio).
    pub pe: f64,
    /// Reynolds number.
    pub re: f64,
}
impl MixingLBM {
    /// Create a new `MixingLBM` with zero concentration and unit equilibrium.
    pub fn new(nx: usize, ny: usize, pe: f64, re: f64) -> Self {
        let n = nx * ny;
        let w = 1.0 / 9.0;
        let f_u = vec![w; n * 9];
        Self {
            nx,
            ny,
            c_scalar: vec![0.0; n],
            f_c: vec![0.0; n * 9],
            f_u,
            pe,
            re,
        }
    }
    /// Row-major flat index for node (i, j).
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Index into the velocity distribution function for node (i,j), direction q.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }
    /// Index into the scalar distribution function for node (i,j), direction q.
    #[inline]
    pub fn fc(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }
    /// Initialise the concentration field with alternating high/low stripes
    /// along the x direction.
    ///
    /// # Arguments
    /// * `n_stripes` – Number of stripes.  Even stripes get `c = 1`, odd `c = 0`.
    pub fn init_stripe_concentration(&mut self, n_stripes: usize) {
        if n_stripes == 0 {
            return;
        }
        let stripe_width = (self.nx as f64 / n_stripes as f64).max(1.0);
        for j in 0..self.ny {
            for i in 0..self.nx {
                let stripe_idx = (i as f64 / stripe_width) as usize;
                let c = if stripe_idx.is_multiple_of(2) {
                    1.0
                } else {
                    0.0
                };
                let idx = self.index(i, j);
                self.c_scalar[idx] = c;
            }
        }
    }
    /// Initialise the concentration field with a Gaussian blob.
    ///
    /// `c(x,y) = c0 * exp(-((x-cx)^2 + (y-cy)^2) / (2σ^2))`
    ///
    /// # Arguments
    /// * `cx`    – x-coordinate of blob centre.
    /// * `cy`    – y-coordinate of blob centre.
    /// * `sigma` – Width parameter σ.
    /// * `c0`   – Peak concentration.
    pub fn init_gaussian_blob(&mut self, cx: f64, cy: f64, sigma: f64, c0: f64) {
        let s2 = sigma * sigma;
        for j in 0..self.ny {
            for i in 0..self.nx {
                let dx = i as f64 - cx;
                let dy = j as f64 - cy;
                let idx = self.index(i, j);
                self.c_scalar[idx] = c0 * (-(dx * dx + dy * dy) / (2.0 * s2)).exp();
            }
        }
    }
    /// Perform one LBM step for both velocity and scalar fields.
    ///
    /// Applies a simple BGK-relaxation to `f_u` and a passive-scalar
    /// relaxation to `f_c`, both using D2Q9 equilibrium distributions.
    /// The velocity is fixed to zero (pressure-driven diffusion only).
    pub fn step(&mut self) {
        let tau_u = 0.5 + 1.0 / (6.0 * self.re.max(1e-6));
        let d_eff = 1.0 / (6.0 * self.pe.max(1e-6));
        let tau_c = 0.5 + d_eff;
        let n = self.nx * self.ny;
        let w0 = 4.0 / 9.0;
        let ws = 1.0 / 9.0;
        for k in 0..n {
            for q in 0..9 {
                let w_eq = if q == 0 { w0 } else { ws };
                let idx = k * 9 + q;
                self.f_u[idx] += -(self.f_u[idx] - w_eq) / tau_u;
            }
        }
        let old_c = self.c_scalar.clone();
        for (k, &c_k) in old_c.iter().enumerate() {
            let c_eq = c_k / 9.0;
            for q in 0..9 {
                let idx = k * 9 + q;
                self.f_c[idx] += -(self.f_c[idx] - c_eq) / tau_c;
            }
            self.c_scalar[k] = self.f_c[k * 9..k * 9 + 9].iter().sum();
        }
    }
}
/// Double distribution function (DDF): separate `f` for flow, `g` for scalar.
///
/// This struct bundles two D2Q9 distribution function arrays, one for the
/// fluid momentum (`f`) and one for the temperature/species scalar (`g`),
/// as commonly used in thermal and reactive LBM.
#[derive(Debug, Clone)]
pub struct DoubleDistributionFunction {
    /// Number of nodes in x.
    pub nx: usize,
    /// Number of nodes in y.
    pub ny: usize,
    /// Fluid distribution functions (D2Q9, length `nx * ny * 9`).
    pub f: Vec<f64>,
    /// Scalar (temperature/species) distribution functions (length `nx * ny * 9`).
    pub g: Vec<f64>,
    /// BGK relaxation time for the fluid.
    pub tau_f: f64,
    /// BGK relaxation time for the scalar.
    pub tau_g: f64,
    /// Fluid density field (length `nx * ny`).
    pub rho: Vec<f64>,
    /// Scalar field (length `nx * ny`).
    pub phi: Vec<f64>,
}
impl DoubleDistributionFunction {
    /// Create a new `DoubleDistributionFunction` with equilibrium initialisation.
    ///
    /// # Arguments
    /// * `nx`, `ny`   – Grid dimensions.
    /// * `tau_f`      – Relaxation time for the fluid.
    /// * `tau_g`      – Relaxation time for the scalar.
    /// * `rho0`       – Initial uniform fluid density.
    /// * `phi0`       – Initial uniform scalar value.
    pub fn new(nx: usize, ny: usize, tau_f: f64, tau_g: f64, rho0: f64, phi0: f64) -> Self {
        let n = nx * ny;
        let mut f = vec![0.0_f64; n * 9];
        let mut g = vec![0.0_f64; n * 9];
        for k in 0..n {
            for q in 0..9 {
                f[k * 9 + q] = D2Q9_W[q] * rho0;
                g[k * 9 + q] = D2Q9_W[q] * phi0;
            }
        }
        Self {
            nx,
            ny,
            f,
            g,
            tau_f,
            tau_g,
            rho: vec![rho0; n],
            phi: vec![phi0; n],
        }
    }
    /// BGK collision for the fluid `f`.
    pub fn collide_flow(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let rho_k = self.rho[k];
            for (q, &w) in D2Q9_W.iter().enumerate() {
                let feq = w * rho_k;
                self.f[k * 9 + q] -= (self.f[k * 9 + q] - feq) / self.tau_f;
            }
            self.rho[k] = self.f[k * 9..k * 9 + 9].iter().sum();
        }
    }
    /// BGK collision for the scalar `g`.
    pub fn collide_scalar(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let phi_k = self.phi[k];
            for (q, &w) in D2Q9_W.iter().enumerate() {
                let geq = w * phi_k;
                self.g[k * 9 + q] -= (self.g[k * 9 + q] - geq) / self.tau_g;
            }
            self.phi[k] = self.g[k * 9..k * 9 + 9].iter().sum();
        }
    }
    /// Perform one coupled step: collide both `f` and `g`.
    pub fn step(&mut self) {
        self.collide_flow();
        self.collide_scalar();
    }
    /// Total fluid mass (sum of rho field).
    pub fn total_mass(&self) -> f64 {
        self.rho.iter().sum()
    }
    /// Total scalar (sum of phi field).
    pub fn total_scalar(&self) -> f64 {
        self.phi.iter().sum()
    }
}
/// Full 2-D multi-species LBM mixing simulation.
///
/// Advances all component density fields via a simple explicit diffusion step.
#[derive(Debug, Clone)]
pub struct MixingSimulation {
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Chemical species being tracked.
    pub components: Vec<MixingComponent>,
}
impl MixingSimulation {
    /// Create a new mixing simulation.
    pub fn new(nx: usize, ny: usize, components: Vec<MixingComponent>) -> Self {
        Self { nx, ny, components }
    }
    /// Advance the simulation by one explicit diffusion step `dt`.
    ///
    /// Uses a simple 5-point Laplacian (Δx = Δy = 1) with Dirichlet zero
    /// boundary conditions on all edges.
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for comp in &mut self.components {
            let d = comp.diffusivity;
            let old = comp.density_field.clone();
            for j in 1..ny.saturating_sub(1) {
                for i in 1..nx.saturating_sub(1) {
                    let idx = j * nx + i;
                    let lap = old[(j - 1) * nx + i]
                        + old[(j + 1) * nx + i]
                        + old[j * nx + (i - 1)]
                        + old[j * nx + (i + 1)]
                        - 4.0 * old[idx];
                    comp.density_field[idx] = old[idx] + d * dt * lap;
                }
            }
        }
    }
    /// Return the total mass across all components.
    pub fn total_mass(&self) -> f64 {
        self.components.iter().map(|c| c.total_mass()).sum()
    }
}
/// Parameters for LBM-based species mixing / multi-component transport.
#[derive(Debug, Clone)]
pub struct MixingLbmParams {
    /// Molecular diffusivity (m² s⁻¹).
    pub diffusivity: f64,
    /// Schmidt number Sc = ν / D.
    pub schmidt_number: f64,
    /// Number of species components.
    pub ncomp: usize,
    /// Lattice time step (lattice units).
    pub dt_lbm: f64,
    /// Lattice spacing (lattice units, typically 1).
    pub dx_lbm: f64,
}
impl MixingLbmParams {
    /// Create new `MixingLbmParams`.
    ///
    /// # Arguments
    /// * `diffusivity`    – Molecular diffusivity.
    /// * `schmidt_number` – Schmidt number ν/D.
    /// * `ncomp`          – Number of species.
    pub fn new(diffusivity: f64, schmidt_number: f64, ncomp: usize) -> Self {
        Self {
            diffusivity,
            schmidt_number,
            ncomp,
            dt_lbm: 1.0_f64,
            dx_lbm: 1.0_f64,
        }
    }
    /// Relaxation time τ for the scalar distribution function.
    ///
    /// `τ = 0.5 + D / (c_s² Δt)` where c_s² = 1/3.
    pub fn tau_scalar(&self) -> f64 {
        0.5_f64 + 3.0_f64 * self.diffusivity * self.dt_lbm / (self.dx_lbm * self.dx_lbm)
    }
    /// Kinematic viscosity inferred from Sc and D: `ν = Sc · D`.
    pub fn kinematic_viscosity(&self) -> f64 {
        self.schmidt_number * self.diffusivity
    }
    /// Relaxation time for the flow distribution function from viscosity.
    pub fn tau_flow(&self) -> f64 {
        0.5_f64 + 3.0_f64 * self.kinematic_viscosity() * self.dt_lbm / (self.dx_lbm * self.dx_lbm)
    }
}
/// A single chemical species participating in the mixing simulation.
///
/// Stores the spatial density field as a flat row-major `Vec`f64` of length
/// `nx * ny`, along with physical parameters.
#[derive(Debug, Clone)]
pub struct MixingComponent {
    /// Spatial density field, row-major (length `nx * ny`).
    pub density_field: Vec<f64>,
    /// Molecular diffusivity (m² s⁻¹).
    pub diffusivity: f64,
    /// Molecular weight (g mol⁻¹).
    pub molecular_weight: f64,
}
impl MixingComponent {
    /// Create a new component, initialising the density field with
    /// `initial_concentration(x)` evaluated at equi-spaced x ∈ [0, 1].
    ///
    /// # Arguments
    /// * `nx`                   – Number of grid points in x.
    /// * `ny`                   – Number of grid points in y.
    /// * `diffusivity`          – Molecular diffusivity (m² s⁻¹).
    /// * `molecular_weight`     – Molecular weight (g mol⁻¹).
    /// * `initial_concentration`– Function `f(x) -> c` providing the initial
    ///   concentration profile along x ∈ [0, 1].
    pub fn new(
        nx: usize,
        ny: usize,
        diffusivity: f64,
        molecular_weight: f64,
        initial_concentration: impl Fn(f64) -> f64,
    ) -> Self {
        let mut density_field = vec![0.0_f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = if nx > 1 {
                    i as f64 / (nx - 1) as f64
                } else {
                    0.0
                };
                density_field[j * nx + i] = initial_concentration(x);
            }
        }
        Self {
            density_field,
            diffusivity,
            molecular_weight,
        }
    }
    /// Total integrated mass (sum of density field values).
    pub fn total_mass(&self) -> f64 {
        self.density_field.iter().sum()
    }
}
