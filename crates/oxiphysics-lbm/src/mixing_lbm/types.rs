//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    D2Q9_W, d2q9_equilibrium, d2q9_moments, fick_flux, laplacian_d2q9_periodic,
    lattice_diffusivity, stream_populations, sync_populations,
};
use crate::error::{Error, Result};

/// Multi-component LBM with cross-diffusion via Fick's law.
///
/// Every species α lives on its own periodic D2Q9 advection–diffusion lattice
/// at rest, with relaxation time `tau[α]`; its self-diffusivity is
/// `D_αα = c_s² (τ_α − ½)` ([`MultiComponentLbm::species_diffusivity`]).
///
/// Cross-diffusion between components is handled via explicit source terms:
/// Fick's law `J_α = −Σ_β D_αβ ∇φ_β` adds `Σ_{β≠α} D_αβ ∇²φ_β` to the balance
/// of species α, with `∇²` the isotropic D2Q9 Laplacian
/// ([`MultiComponentLbm::add_cross_diffusion_source`]).  The diagonal of
/// `cross_diff` is ignored — the self-diffusion belongs to `tau` — and the
/// total of every species is conserved exactly, with or without
/// cross-diffusion.  The explicit coupling is stable for `D_αβ` well below `1/4`.
///
/// [`MultiComponentLbm::step`] performs one full time step (collision,
/// cross-diffusion source, streaming).  The fields `phi` are the zeroth
/// moments of `g`; they may be edited directly between steps (the next
/// collision moves the difference onto `g`) or set with
/// [`MultiComponentLbm::set_species_field`].
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
    /// Persistent streaming buffer, so that streaming does not allocate.
    scratch: Vec<f64>,
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
            scratch: Vec::new(),
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
    /// Self-diffusivity `D_αα = c_s² (τ_α − ½)` of species `alpha`, or `None`
    /// when there is no such species.
    pub fn species_diffusivity(&self, alpha: usize) -> Option<f64> {
        self.tau.get(alpha).map(|&tau| lattice_diffusivity(tau))
    }
    /// Set the concentration field of species `alpha` (row-major, `nx · ny`
    /// values) and bring its populations to the matching equilibrium.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidSpeciesCount`] when `alpha >= ncomp` and
    /// [`Error::IncompatibleParameters`] when `field.len() != nx * ny`.
    pub fn set_species_field(&mut self, alpha: usize, field: &[f64]) -> Result<()> {
        if alpha >= self.ncomp {
            return Err(Error::InvalidSpeciesCount {
                min_species: alpha + 1,
                actual: self.ncomp,
            });
        }
        let n = self.nx * self.ny;
        if field.len() != n {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "species field needs {n} values for a {} x {} grid, got {}",
                    self.nx,
                    self.ny,
                    field.len()
                ),
            });
        }
        if let Some(phi) = self.phi.get_mut(alpha) {
            phi.clear();
            phi.extend_from_slice(field);
        }
        if let Some(g) = self.g.get_mut(alpha) {
            g.clear();
            g.resize(n * 9, 0.0);
            for (node, &value) in g.chunks_exact_mut(9).zip(field.iter()) {
                for (p, &w) in node.iter_mut().zip(D2Q9_W.iter()) {
                    *p = w * value;
                }
            }
        }
        Ok(())
    }
    /// Perform one BGK collision for all species without cross-diffusion.
    ///
    /// Every species relaxes toward `w_q φ_α` with its own `τ_α`; the
    /// collision conserves each species' mass.
    pub fn collide_uncoupled(&mut self) {
        for ((g, phi), &tau) in self
            .g
            .iter_mut()
            .zip(self.phi.iter_mut())
            .zip(self.tau.iter())
        {
            sync_populations(g, phi);
            for (node, phi_k) in g.chunks_exact_mut(9).zip(phi.iter_mut()) {
                let density: f64 = node.iter().sum();
                for (p, &w) in node.iter_mut().zip(D2Q9_W.iter()) {
                    *p -= (*p - w * density) / tau;
                }
                *phi_k = node.iter().sum();
            }
        }
    }
    /// Whether any off-diagonal entry of `cross_diff` is non-zero, i.e. whether
    /// [`MultiComponentLbm::add_cross_diffusion_source`] changes anything.
    fn has_cross_coupling(&self) -> bool {
        (0..self.ncomp).any(|alpha| {
            (0..self.ncomp).any(|beta| {
                beta != alpha
                    && self
                        .cross_diff
                        .get(alpha * self.ncomp + beta)
                        .is_some_and(|&d_ab| d_ab != 0.0)
            })
        })
    }
    /// Add the cross-diffusion source `Σ_{β≠α} D_αβ ∇²φ_β` to every species α.
    ///
    /// The source is distributed on the populations of species α with the
    /// lattice weights, so it changes the field but not the momentum (there
    /// is none) and the total of each species is unchanged: the periodic
    /// Laplacian sums to zero.  All Laplacians are evaluated from the fields
    /// before any source is added.
    pub fn add_cross_diffusion_source(&mut self) {
        let n = self.nx * self.ny;
        if self.ncomp < 2 || n == 0 || !self.has_cross_coupling() {
            return;
        }
        let laplacians: Vec<Vec<f64>> = self
            .phi
            .iter()
            .map(|phi| laplacian_d2q9_periodic(phi, self.nx, self.ny))
            .collect();
        for (alpha, (g, phi)) in self.g.iter_mut().zip(self.phi.iter_mut()).enumerate() {
            let mut source = vec![0.0_f64; n];
            let mut driven = false;
            for (beta, laplacian) in laplacians.iter().enumerate() {
                if beta == alpha {
                    continue;
                }
                let d_ab = self
                    .cross_diff
                    .get(alpha * self.ncomp + beta)
                    .copied()
                    .unwrap_or(0.0);
                if d_ab == 0.0 {
                    continue;
                }
                driven = true;
                for (s, &lap) in source.iter_mut().zip(laplacian.iter()) {
                    *s += d_ab * lap;
                }
            }
            if !driven {
                continue;
            }
            for ((node, phi_k), &s) in g.chunks_exact_mut(9).zip(phi.iter_mut()).zip(source.iter())
            {
                for (p, &w) in node.iter_mut().zip(D2Q9_W.iter()) {
                    *p += w * s;
                }
                *phi_k = node.iter().sum();
            }
        }
    }
    /// Stream the populations of every species along the D2Q9 links with
    /// periodic wrap-around and refresh the fields `phi`.
    pub fn stream(&mut self) {
        for (g, phi) in self.g.iter_mut().zip(self.phi.iter_mut()) {
            stream_populations(g, &mut self.scratch, self.nx, self.ny);
            for (p, node) in phi.iter_mut().zip(g.chunks_exact(9)) {
                *p = node.iter().sum();
            }
        }
    }
    /// Advance all species by one time step: BGK collision, cross-diffusion
    /// source, periodic streaming.
    pub fn step(&mut self) {
        self.collide_uncoupled();
        self.add_cross_diffusion_source();
        self.stream();
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
/// A BGK-LBM step evolves the scalar distribution functions `g`: collision
/// toward the equilibrium advected with the prescribed velocity field
/// (`ux`, `uy`), then periodic streaming.  The scalar diffuses with
/// `D = c_s² (τ − ½)` ([`PassiveScalarField::diffusivity`]).
///
/// `phi` is the zeroth moment of `g`.  It may be edited directly between
/// steps (the next collision moves the difference onto `g`, keeping the edited
/// mass) or set with [`PassiveScalarField::set_scalar_field`].
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
    /// Persistent streaming buffer, so that a step does not allocate.
    scratch: Vec<f64>,
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
            scratch: Vec::new(),
        }
    }
    /// Set the velocity field (flat, row-major, `nx * ny` values each).
    pub fn set_velocity(&mut self, ux: Vec<f64>, uy: Vec<f64>) {
        self.ux = ux;
        self.uy = uy;
    }
    /// Diffusivity `D = c_s² (τ − ½)` of the scalar.
    pub fn diffusivity(&self) -> f64 {
        lattice_diffusivity(self.tau)
    }
    /// Set the scalar field (row-major, `nx · ny` values) and bring `g` to the
    /// equilibrium advected with the current velocity field.
    ///
    /// # Errors
    ///
    /// [`Error::IncompatibleParameters`] when `phi`, `ux` or `uy` does not have
    /// `nx * ny` values.
    pub fn set_scalar_field(&mut self, phi: &[f64]) -> Result<()> {
        let n = self.nx * self.ny;
        if phi.len() != n || self.ux.len() != n || self.uy.len() != n {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "scalar and velocity fields need {n} values for a {} x {} grid, \
                     got {}, {} and {}",
                    self.nx,
                    self.ny,
                    phi.len(),
                    self.ux.len(),
                    self.uy.len()
                ),
            });
        }
        self.phi.clear();
        self.phi.extend_from_slice(phi);
        self.g.clear();
        self.g.resize(n * 9, 0.0);
        for (k, node) in self.g.chunks_exact_mut(9).enumerate() {
            for (q, p) in node.iter_mut().enumerate() {
                *p = d2q9_equilibrium(q, phi[k], self.ux[k], self.uy[k]);
            }
        }
        Ok(())
    }
    /// Perform one BGK collision step (no streaming) for the scalar.
    ///
    /// Every node relaxes toward the second-order equilibrium advected with
    /// its velocity; the collision conserves the scalar of each node.
    pub fn collide(&mut self) {
        sync_populations(&mut self.g, &self.phi);
        let tau = self.tau;
        for (k, node) in self.g.chunks_exact_mut(9).enumerate() {
            let density: f64 = node.iter().sum();
            let (ux_k, uy_k) = (self.ux[k], self.uy[k]);
            for (q, p) in node.iter_mut().enumerate() {
                *p -= (*p - d2q9_equilibrium(q, density, ux_k, uy_k)) / tau;
            }
            self.phi[k] = node.iter().sum();
        }
    }
    /// Stream the scalar populations along the D2Q9 links with periodic
    /// wrap-around and refresh `phi`.
    pub fn stream(&mut self) {
        stream_populations(&mut self.g, &mut self.scratch, self.nx, self.ny);
        for (p, node) in self.phi.iter_mut().zip(self.g.chunks_exact(9)) {
            *p = node.iter().sum();
        }
    }
    /// Advance the scalar by one time step: BGK collision, then streaming.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
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
/// Two D2Q9 populations evolve on a fully periodic `nx × ny` lattice (lattice
/// units, `Δx = Δt = 1`):
///
/// * `f_u` is the flow: a BGK fluid with relaxation time
///   [`MixingLBM::tau_u`] whose density and velocity are the zeroth and first
///   moment of `f_u`.  [`MixingLBM::new`] starts it at rest with unit density,
///   so the velocity stays zero (pure diffusion) until
///   [`MixingLBM::set_uniform_velocity`] sets a flow.
/// * `f_c` is the passive scalar: a BGK advection–diffusion lattice with
///   relaxation time [`MixingLBM::tau_c`].  Its equilibrium
///   `w_q c (1 + 3 c_q·u + 9/2 (c_q·u)² − 3/2 |u|²)` is advected with the local
///   flow velocity `u`, and the scalar diffuses with `D = c_s² (τ_c − ½)`
///   ([`MixingLBM::scalar_diffusivity`]).
///
/// The relaxation times follow from the dimensionless numbers as
/// `τ_u = ½ + 1/(6 Re)` and `τ_c = ½ + 1/(6 Pe)` (the numbers are clamped
/// below at `1e-6`).
///
/// `c_scalar` is the zeroth moment of `f_c` and is refreshed by every
/// [`MixingLBM::step`].  Set the concentration with
/// [`MixingLBM::init_stripe_concentration`], [`MixingLBM::init_gaussian_blob`]
/// or [`MixingLBM::set_concentration`], which write both; if `c_scalar` is
/// edited directly, the next step moves the difference onto `f_c` so that the
/// edited mass is kept.
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
    /// Persistent streaming buffer, so that a step does not allocate.
    scratch: Vec<f64>,
}
impl MixingLBM {
    /// Create a new `MixingLBM` with zero concentration and a unit-density
    /// fluid at rest.
    pub fn new(nx: usize, ny: usize, pe: f64, re: f64) -> Self {
        let n = nx * ny;
        let mut f_u = vec![0.0_f64; n * 9];
        for node in f_u.chunks_exact_mut(9) {
            node.copy_from_slice(&D2Q9_W);
        }
        Self {
            nx,
            ny,
            c_scalar: vec![0.0; n],
            f_c: vec![0.0; n * 9],
            f_u,
            pe,
            re,
            scratch: Vec::new(),
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
    /// BGK relaxation time of the flow, `τ_u = ½ + 1/(6 Re)`.
    pub fn tau_u(&self) -> f64 {
        0.5 + 1.0 / (6.0 * self.re.max(1e-6))
    }
    /// BGK relaxation time of the scalar, `τ_c = ½ + 1/(6 Pe)`.
    pub fn tau_c(&self) -> f64 {
        0.5 + 1.0 / (6.0 * self.pe.max(1e-6))
    }
    /// Kinematic viscosity of the flow, `c_s² (τ_u − ½)`.
    pub fn flow_viscosity(&self) -> f64 {
        lattice_diffusivity(self.tau_u())
    }
    /// Diffusivity of the scalar, `D = c_s² (τ_c − ½)`.
    pub fn scalar_diffusivity(&self) -> f64 {
        lattice_diffusivity(self.tau_c())
    }
    /// Set the concentration of node `idx` and put `f_c` there at the
    /// equilibrium advected with the local flow velocity.
    fn set_node_concentration(&mut self, idx: usize, c: f64) {
        let (_, ux, uy) = d2q9_moments(&self.f_u[idx * 9..idx * 9 + 9]);
        self.c_scalar[idx] = c;
        for (q, p) in self.f_c[idx * 9..idx * 9 + 9].iter_mut().enumerate() {
            *p = d2q9_equilibrium(q, c, ux, uy);
        }
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
                self.set_node_concentration(idx, c);
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
                self.set_node_concentration(idx, c0 * (-(dx * dx + dy * dy) / (2.0 * s2)).exp());
            }
        }
    }
    /// Set the whole concentration field (row-major, `nx · ny` values) and put
    /// `f_c` at the equilibrium advected with the local flow velocity.
    ///
    /// # Errors
    ///
    /// [`Error::IncompatibleParameters`] when `field.len() != nx * ny`.
    pub fn set_concentration(&mut self, field: &[f64]) -> Result<()> {
        let n = self.nx * self.ny;
        if field.len() != n || self.c_scalar.len() != n || self.f_c.len() != n * 9 {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "concentration field needs {n} values for a {} x {} grid, got {}",
                    self.nx,
                    self.ny,
                    field.len()
                ),
            });
        }
        for (idx, &c) in field.iter().enumerate() {
            self.set_node_concentration(idx, c);
        }
        Ok(())
    }
    /// Put the whole fluid in motion with the uniform velocity `(ux, uy)`.
    ///
    /// The flow populations and the scalar populations of every node are set
    /// to the equilibrium for the current density and concentration at the new
    /// velocity, so the state is a fixed point of the collision.
    pub fn set_uniform_velocity(&mut self, ux: f64, uy: f64) {
        for ((fu, fc), &c) in self
            .f_u
            .chunks_exact_mut(9)
            .zip(self.f_c.chunks_exact_mut(9))
            .zip(self.c_scalar.iter())
        {
            let rho: f64 = fu.iter().sum();
            for (q, (u_q, c_q)) in fu.iter_mut().zip(fc.iter_mut()).enumerate() {
                *u_q = d2q9_equilibrium(q, rho, ux, uy);
                *c_q = d2q9_equilibrium(q, c, ux, uy);
            }
        }
    }
    /// Flow velocity `(u_x, u_y)` at node `(i, j)`, or `None` outside the grid.
    pub fn velocity_at(&self, i: usize, j: usize) -> Option<(f64, f64)> {
        if i >= self.nx || j >= self.ny {
            return None;
        }
        let k = self.index(i, j);
        let node = self.f_u.get(k * 9..k * 9 + 9)?;
        let (_, ux, uy) = d2q9_moments(node);
        Some((ux, uy))
    }
    /// Total scalar, the sum of `c_scalar` over all nodes.
    pub fn total_scalar(&self) -> f64 {
        self.c_scalar.iter().sum()
    }
    /// Advance the flow and the scalar by one time step.
    ///
    /// Both populations relax toward their equilibria — the scalar toward the
    /// one advected with the local velocity of `f_u` — are streamed
    /// periodically, and `c_scalar` is refreshed from `f_c`.  The total scalar
    /// and the total fluid mass and momentum are conserved to rounding error.
    pub fn step(&mut self) {
        let tau_u = self.tau_u();
        let tau_c = self.tau_c();
        sync_populations(&mut self.f_c, &self.c_scalar);
        for (fu, fc) in self
            .f_u
            .chunks_exact_mut(9)
            .zip(self.f_c.chunks_exact_mut(9))
        {
            let (rho, ux, uy) = d2q9_moments(fu);
            let c: f64 = fc.iter().sum();
            for (q, (u_q, c_q)) in fu.iter_mut().zip(fc.iter_mut()).enumerate() {
                *u_q -= (*u_q - d2q9_equilibrium(q, rho, ux, uy)) / tau_u;
                *c_q -= (*c_q - d2q9_equilibrium(q, c, ux, uy)) / tau_c;
            }
        }
        stream_populations(&mut self.f_u, &mut self.scratch, self.nx, self.ny);
        stream_populations(&mut self.f_c, &mut self.scratch, self.nx, self.ny);
        for (c, node) in self.c_scalar.iter_mut().zip(self.f_c.chunks_exact(9)) {
            *c = node.iter().sum();
        }
    }
}
/// Double distribution function (DDF): separate `f` for flow, `g` for scalar.
///
/// This struct bundles two D2Q9 distribution function arrays on a periodic
/// `nx × ny` lattice, one for the fluid momentum (`f`) and one for the
/// temperature/species scalar (`g`), as commonly used in thermal and reactive
/// LBM:
///
/// * `f` is an isothermal BGK fluid with relaxation time `tau_f`.  It relaxes
///   toward the second-order equilibrium of its own density and velocity, so
///   mass and momentum are conserved and the kinematic viscosity is
///   `ν = c_s² (τ_f − ½)` ([`DoubleDistributionFunction::viscosity`]).
/// * `g` is a passive scalar with relaxation time `tau_g`.  It relaxes toward
///   the equilibrium advected with the local velocity of `f` and diffuses with
///   `D = c_s² (τ_g − ½)` ([`DoubleDistributionFunction::scalar_diffusivity`]).
///
/// `rho` and `phi` are the zeroth moments of `f` and `g`.  They may be edited
/// directly between steps (the next collision moves the difference onto the
/// populations) or set per node with [`DoubleDistributionFunction::set_flow_at`]
/// and [`DoubleDistributionFunction::set_scalar_at`].
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
    /// Persistent streaming buffer, so that a step does not allocate.
    scratch: Vec<f64>,
}
impl DoubleDistributionFunction {
    /// Create a new `DoubleDistributionFunction` with equilibrium initialisation.
    ///
    /// The fluid starts at rest.
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
            scratch: Vec::new(),
        }
    }
    /// Kinematic viscosity of the fluid, `ν = c_s² (τ_f − ½)`.
    pub fn viscosity(&self) -> f64 {
        lattice_diffusivity(self.tau_f)
    }
    /// Diffusivity of the scalar, `D = c_s² (τ_g − ½)`.
    pub fn scalar_diffusivity(&self) -> f64 {
        lattice_diffusivity(self.tau_g)
    }
    /// Fluid velocity `(u_x, u_y)` at node `(x, y)`, or `None` outside the grid.
    pub fn velocity_at(&self, x: usize, y: usize) -> Option<(f64, f64)> {
        if x >= self.nx || y >= self.ny {
            return None;
        }
        let k = y * self.nx + x;
        let node = self.f.get(k * 9..k * 9 + 9)?;
        let (_, ux, uy) = d2q9_moments(node);
        Some((ux, uy))
    }
    /// Put node `(x, y)` of the fluid at the equilibrium with density `rho`
    /// and velocity `(ux, uy)`.
    ///
    /// # Errors
    ///
    /// [`Error::IndexOutOfBounds`] when the node is outside the grid.
    pub fn set_flow_at(&mut self, x: usize, y: usize, rho: f64, ux: f64, uy: f64) -> Result<()> {
        Error::check_index(x, y, self.nx, self.ny)?;
        let k = y * self.nx + x;
        for (q, p) in self.f[k * 9..k * 9 + 9].iter_mut().enumerate() {
            *p = d2q9_equilibrium(q, rho, ux, uy);
        }
        self.rho[k] = rho;
        Ok(())
    }
    /// Put node `(x, y)` of the scalar at the equilibrium advected with the
    /// current fluid velocity of that node, with value `phi`.
    ///
    /// # Errors
    ///
    /// [`Error::IndexOutOfBounds`] when the node is outside the grid.
    pub fn set_scalar_at(&mut self, x: usize, y: usize, phi: f64) -> Result<()> {
        Error::check_index(x, y, self.nx, self.ny)?;
        let k = y * self.nx + x;
        let (_, ux, uy) = d2q9_moments(&self.f[k * 9..k * 9 + 9]);
        for (q, p) in self.g[k * 9..k * 9 + 9].iter_mut().enumerate() {
            *p = d2q9_equilibrium(q, phi, ux, uy);
        }
        self.phi[k] = phi;
        Ok(())
    }
    /// BGK collision for the fluid `f`.
    ///
    /// Every node relaxes toward the second-order equilibrium of its own
    /// density and velocity, which conserves its mass and momentum.
    pub fn collide_flow(&mut self) {
        sync_populations(&mut self.f, &self.rho);
        let tau = self.tau_f;
        for (node, rho_k) in self.f.chunks_exact_mut(9).zip(self.rho.iter_mut()) {
            let (rho, ux, uy) = d2q9_moments(node);
            for (q, p) in node.iter_mut().enumerate() {
                *p -= (*p - d2q9_equilibrium(q, rho, ux, uy)) / tau;
            }
            *rho_k = node.iter().sum();
        }
    }
    /// BGK collision for the scalar `g`.
    ///
    /// Every node relaxes toward the equilibrium advected with the local fluid
    /// velocity (the first moment of `f`); the scalar of each node is conserved.
    pub fn collide_scalar(&mut self) {
        sync_populations(&mut self.g, &self.phi);
        let tau = self.tau_g;
        for ((g_node, f_node), phi_k) in self
            .g
            .chunks_exact_mut(9)
            .zip(self.f.chunks_exact(9))
            .zip(self.phi.iter_mut())
        {
            let (_, ux, uy) = d2q9_moments(f_node);
            let density: f64 = g_node.iter().sum();
            for (q, p) in g_node.iter_mut().enumerate() {
                *p -= (*p - d2q9_equilibrium(q, density, ux, uy)) / tau;
            }
            *phi_k = g_node.iter().sum();
        }
    }
    /// Stream `f` and `g` along the D2Q9 links with periodic wrap-around and
    /// refresh `rho` and `phi`.
    pub fn stream(&mut self) {
        stream_populations(&mut self.f, &mut self.scratch, self.nx, self.ny);
        stream_populations(&mut self.g, &mut self.scratch, self.nx, self.ny);
        for (rho, node) in self.rho.iter_mut().zip(self.f.chunks_exact(9)) {
            *rho = node.iter().sum();
        }
        for (phi, node) in self.phi.iter_mut().zip(self.g.chunks_exact(9)) {
            *phi = node.iter().sum();
        }
    }
    /// Perform one coupled step: collide both `f` and `g`, then stream them.
    pub fn step(&mut self) {
        self.collide_flow();
        self.collide_scalar();
        self.stream();
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
