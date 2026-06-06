//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Tracks the SGS energy budget: production, dissipation, and transfer.
#[derive(Debug, Clone)]
pub struct LesEnergyBudget {
    /// Accumulated SGS production Σ P_sgs.
    pub production: f64,
    /// Accumulated SGS dissipation Σ ε_sgs.
    pub dissipation: f64,
    /// Number of particles contributing.
    pub count: usize,
}
impl LesEnergyBudget {
    /// Create a new zeroed energy budget.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add the contribution of a single particle to the budget.
    ///
    /// - `tau` : SGS stress tensor (Pa).
    /// - `s`   : resolved strain-rate tensor (s⁻¹).
    pub fn add_particle(&mut self, tau: &Mat3, s: &Mat3) {
        let p_sgs = -mat3_double_contraction(tau, s);
        self.production += p_sgs;
        self.dissipation += p_sgs.max(0.0);
        self.count += 1;
    }
    /// Mean SGS production per particle.
    pub fn mean_production(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.production / self.count as f64
        }
    }
    /// Reset the budget counters.
    pub fn reset(&mut self) {
        self.production = 0.0;
        self.dissipation = 0.0;
        self.count = 0;
    }
}
/// Van-Driest damped mixing-length model for near-wall turbulence.
#[derive(Debug, Clone)]
pub struct MixingLength {
    /// Dimensional wall distance (m) used for reference (informational).
    pub wall_distance: f64,
    /// Von Kármán constant κ (default 0.41).
    pub von_karman: f64,
    /// Viscous length scale ν/u_τ (m).  Used to convert y → y⁺.
    pub viscous_length: f64,
}
impl MixingLength {
    /// Create a new mixing-length model.
    pub fn new(wall_distance: f64, von_karman: f64, viscous_length: f64) -> Self {
        Self {
            wall_distance,
            von_karman,
            viscous_length,
        }
    }
    /// Compute the van-Driest-damped mixing length at wall-normal distance `y`.
    ///
    /// l_m = κ · y · (1 - exp(-y⁺/26))
    ///
    /// where y⁺ = y / viscous_length.
    pub fn mixing_length(&self, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        let y_plus = y / self.viscous_length;
        self.von_karman * y * (1.0 - (-y_plus / 26.0).exp())
    }
}
/// Sub-Particle Scale (SPS) turbulence model (Gotoh et al. 2001).
///
/// Accounts for sub-particle-scale turbulent stresses via a Smagorinsky-type
/// eddy-viscosity closure.
#[derive(Debug, Clone)]
pub struct SpsModel {
    /// Smagorinsky constant (default 0.12).
    pub cs: f64,
    /// SPS turbulent kinetic-energy constant (default 0.0066).
    pub ci: f64,
    /// Nominal particle spacing Δx (m).
    pub particle_spacing: f64,
}
impl SpsModel {
    /// Create a new SPS model.
    pub fn new(cs: f64, ci: f64, particle_spacing: f64) -> Self {
        Self {
            cs,
            ci,
            particle_spacing,
        }
    }
    /// Compute the symmetric strain-rate tensor S_ij for particle `i`.
    ///
    /// Uses the standard SPH gradient:
    ///   ∂v_α/∂x_β ≈ Σ_j (m_j/ρ_j) (v_j_α - v_i_α) ∂W_ij/∂x_β
    ///
    /// `kernel_grads[k]` is the gradient vector of W(r_ij, h) evaluated at
    /// the position of neighbour `neighbors[k]`, pointing from i → j:
    ///   ∇_i W_ij = (dW/dr) * r_hat_ij  (a \[f64;3\] value).
    pub fn compute_strain_rate_tensor(
        &self,
        i: usize,
        neighbors: &[usize],
        _positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        kernel_grads: &[[f64; 3]],
    ) -> Mat3 {
        let mut dv: Mat3 = mat3_zero();
        let vi = velocities[i];
        let rho_i = densities[i].max(1e-14);
        for (k, &j) in neighbors.iter().enumerate() {
            let rho_j = densities[j].max(1e-14);
            let weight = masses[j] / rho_j;
            let dvij = [
                velocities[j][0] - vi[0],
                velocities[j][1] - vi[1],
                velocities[j][2] - vi[2],
            ];
            let grad = kernel_grads[k];
            for alpha in 0..3 {
                for beta in 0..3 {
                    dv[alpha][beta] += weight * dvij[alpha] * grad[beta];
                }
            }
        }
        let mut s: Mat3 = mat3_zero();
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] = 0.5 * (dv[a][b] + dv[b][a]);
            }
        }
        let _ = rho_i;
        s
    }
    /// Compute the SPS turbulent stress tensor τ^SPS_ij.
    ///
    /// τ_ij = ν_turb * S_ij - (2/3) k_sgs δ_ij
    ///
    /// where
    ///   ν_turb = (Cs · Δx)² · |S|
    ///   k_sgs  = (Ci · Δx)² · |S|²
    pub fn compute_sps_stress(&self, s: &Mat3, _rho: f64) -> Mat3 {
        let dx = self.particle_spacing;
        let s_mag = strain_rate_magnitude(s);
        let cs_dx = self.cs * dx;
        let ci_dx = self.ci * dx;
        let nu_turb = cs_dx * cs_dx * s_mag;
        let k_sgs = ci_dx * ci_dx * s_mag * s_mag;
        let mut tau: Mat3 = mat3_zero();
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                tau[a][b] = nu_turb * s[a][b] - (2.0 / 3.0) * k_sgs * delta;
            }
        }
        tau
    }
    /// Compute SPS turbulent forces for all particles.
    ///
    /// The force on particle `i` from the SPS stress divergence is:
    ///   f_i^SPS = m_i * Σ_j m_j * (τ_i/ρ_i² + τ_j/ρ_j²) · ∇W_ij
    ///
    /// `all_kernel_grads[i][k]` is the kernel gradient for the k-th
    /// neighbour of particle i.
    pub fn apply_sps_forces(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        all_kernel_grads: &[Vec<[f64; 3]>],
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let stresses: Vec<Mat3> = (0..n)
            .map(|i| {
                let s = self.compute_strain_rate_tensor(
                    i,
                    &neighbors[i],
                    positions,
                    velocities,
                    masses,
                    densities,
                    &all_kernel_grads[i],
                );
                self.compute_sps_stress(&s, densities[i])
            })
            .collect();
        for i in 0..n {
            let rho_i = densities[i].max(1e-14);
            let tau_i = &stresses[i];
            for (k, &j) in neighbors[i].iter().enumerate() {
                let rho_j = densities[j].max(1e-14);
                let tau_j = &stresses[j];
                let grad = all_kernel_grads[i][k];
                for beta in 0..3 {
                    let mut contrib = 0.0;
                    for alpha in 0..3 {
                        contrib += (tau_i[beta][alpha] / (rho_i * rho_i)
                            + tau_j[beta][alpha] / (rho_j * rho_j))
                            * grad[alpha];
                    }
                    forces[i][beta] += masses[i] * masses[j] * contrib;
                }
            }
        }
        forces
    }
}
/// Large Eddy Simulation (LES) filter using the Smagorinsky sub-grid model.
///
/// The eddy viscosity is: ν_t = (C_s · h)² · |S|
#[derive(Debug, Clone)]
pub struct LesFilter {
    /// Filter width / smoothing length h (m).
    pub h: f64,
    /// Smagorinsky constant C_s (dimensionless, typically ~0.1–0.2).
    pub c_s: f64,
}
impl LesFilter {
    /// Create a new LES filter.
    pub fn new(h: f64, c_s: f64) -> Self {
        Self { h, c_s }
    }
    /// Smagorinsky eddy viscosity: ν_t = (C_s · h)² · |S|
    pub fn smagorinsky_viscosity(&self, s_mag: f64, _rho: f64) -> f64 {
        let cs_h = self.c_s * self.h;
        cs_h * cs_h * s_mag
    }
    /// Turbulent stress tensor: τ_ij = −2 · ρ · ν_t · S_ij
    pub fn turbulent_stress(&self, s: Mat3, rho: f64) -> Mat3 {
        let s_mag = strain_rate_magnitude(&s);
        let nu_t = self.smagorinsky_viscosity(s_mag, rho);
        let mut tau = mat3_zero();
        for a in 0..3 {
            for b in 0..3 {
                tau[a][b] = -2.0 * rho * nu_t * s[a][b];
            }
        }
        tau
    }
}
/// SPH particle carrying LES turbulent fields.
#[derive(Debug, Clone)]
pub struct SphLesParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Turbulent (eddy) viscosity ν_t (m²/s).
    pub nu_t: f64,
    /// Magnitude of the local strain-rate tensor |S| (1/s).
    pub s_mag: f64,
}
/// Classical Smagorinsky large-eddy simulation (LES) sub-grid-scale (SGS)
/// model.
///
/// The SGS eddy viscosity is:
/// ```text
/// ν_t = (C_s Δ)² |S̄|
/// ```
/// where |S̄| = sqrt(2 S̄_ij S̄_ij) is the resolved strain-rate magnitude and
/// Δ is the filter width (particle spacing).
///
/// Reference: Smagorinsky (1963), "General circulation experiments with the
/// primitive equations", Monthly Weather Review.
#[derive(Debug, Clone)]
pub struct SmagorinskyLes {
    /// Smagorinsky constant C_s (typically 0.1–0.2 for turbulent flows).
    pub cs: f64,
    /// Filter width Δ (particle spacing, m).
    pub delta: f64,
}
impl SmagorinskyLes {
    /// Create a new Smagorinsky LES model.
    pub fn new(cs: f64, delta: f64) -> Self {
        Self { cs, delta }
    }
    /// Compute the SGS eddy viscosity ν_t = (C_s Δ)² |S̄|.
    ///
    /// `strain_rate` is the symmetric resolved strain-rate tensor S̄_ij.
    pub fn eddy_viscosity(&self, strain_rate: &Mat3) -> f64 {
        let s_mag = strain_rate_magnitude(strain_rate);
        let cs_delta_sq = (self.cs * self.delta) * (self.cs * self.delta);
        cs_delta_sq * s_mag
    }
    /// Compute the SGS stress tensor τ_ij = -2 ρ ν_t S̄_ij.
    pub fn sgs_stress(&self, strain_rate: &Mat3, density: f64) -> Mat3 {
        let nu_t = self.eddy_viscosity(strain_rate);
        mat3_scale(-2.0 * density * nu_t, *strain_rate)
    }
    /// Compute the SGS acceleration contribution for particle i:
    ///
    /// ```text
    /// a_sgs_α = (1/ρ_i) Σ_j m_j (τ_i/ρ_i² + τ_j/ρ_j²) · ∇W_ij
    /// ```
    ///
    /// `tau_i` and `tau_j` are the SGS stress tensors; `grad_w` is the
    /// kernel gradient ∇W_ij pointing from i to j.
    pub fn sgs_acceleration(
        &self,
        tau_i: &Mat3,
        rho_i: f64,
        neighbors_tau: &[Mat3],
        neighbors_rho: &[f64],
        neighbor_masses: &[f64],
        kernel_grads: &[[f64; 3]],
    ) -> [f64; 3] {
        let mut acc = [0.0_f64; 3];
        let rho_i2 = rho_i * rho_i;
        for (k, &rho_j) in neighbors_rho.iter().enumerate() {
            let rho_j2 = rho_j * rho_j;
            let m_j = neighbor_masses[k];
            let tau_j = &neighbors_tau[k];
            let grad = kernel_grads[k];
            for alpha in 0..3 {
                let mut sum_beta = 0.0;
                for beta in 0..3 {
                    sum_beta +=
                        (tau_i[alpha][beta] / rho_i2 + tau_j[alpha][beta] / rho_j2) * grad[beta];
                }
                acc[alpha] += m_j * sum_beta;
            }
        }
        acc
    }
}
/// Dynamic Smagorinsky coefficient computed via the Germano–Lilly procedure.
///
/// Uses two filter levels (grid filter Δ and test filter 2Δ) to dynamically
/// adjust the Smagorinsky constant in space and time.
#[derive(Debug, Clone)]
pub struct DynamicSmagorinsky {
    /// Base filter width Δ (particle spacing).
    pub delta: f64,
    /// Ratio of test-filter to grid-filter width (typically 2.0).
    pub test_filter_ratio: f64,
    /// Lower clip for the dynamic coefficient (prevents negative viscosity blow-up).
    pub c_min: f64,
    /// Upper clip for the dynamic coefficient.
    pub c_max: f64,
}
impl DynamicSmagorinsky {
    /// Create a new dynamic Smagorinsky model.
    pub fn new(delta: f64, test_filter_ratio: f64) -> Self {
        Self {
            delta,
            test_filter_ratio,
            c_min: 0.0,
            c_max: 0.25,
        }
    }
    /// Compute the dynamic Smagorinsky constant C_s² via the Germano–Lilly method.
    ///
    /// Arguments:
    /// - `s_grid` : strain-rate tensor at the grid-filter level.
    /// - `s_test` : strain-rate tensor at the test-filter level.
    ///
    /// Returns the clipped dynamic coefficient C_s² (≥ 0).
    pub fn dynamic_coefficient(&self, s_grid: Mat3, s_test: Mat3) -> f64 {
        let delta_hat = self.delta * self.test_filter_ratio;
        let s_bar_mag = strain_rate_magnitude(&s_grid);
        let s_hat_mag = strain_rate_magnitude(&s_test);
        let mut lm = 0.0_f64;
        let mut mm = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let m_ij = 2.0
                    * (delta_hat * delta_hat * s_hat_mag * s_test[i][j]
                        - self.delta * self.delta * s_bar_mag * s_grid[i][j]);
                let l_ij = delta_hat * delta_hat * s_hat_mag * s_test[i][j]
                    - self.delta * self.delta * s_bar_mag * s_grid[i][j];
                lm += l_ij * m_ij;
                mm += m_ij * m_ij;
            }
        }
        if mm < 1e-30 {
            return 0.0;
        }
        (lm / mm).clamp(self.c_min, self.c_max)
    }
    /// Compute the dynamic eddy viscosity ν_t = C_s² · Δ² · |S|.
    pub fn eddy_viscosity(&self, s_grid: Mat3, s_test: Mat3) -> f64 {
        let cs2 = self.dynamic_coefficient(s_grid, s_test);
        let s_mag = strain_rate_magnitude(&s_grid);
        cs2 * self.delta * self.delta * s_mag
    }
    /// Compute the full dynamic SGS stress tensor τ_ij = -2 ρ ν_t S_ij.
    pub fn sgs_stress(&self, s_grid: Mat3, s_test: Mat3, density: f64) -> Mat3 {
        let nu_t = self.eddy_viscosity(s_grid, s_test);
        let mut tau = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                tau[i][j] = -2.0 * density * nu_t * s_grid[i][j];
            }
        }
        tau
    }
}
/// Turbulent diffusion coefficient for the SPH momentum equation.
///
/// Returns the particle contribution to the turbulent stress divergence:
/// (1/ρ) ∂τ_ij/∂x_j ≈ sum_j (m_j / ρ_j) (τ_i^j + τ_i^j) · ∇W_ij
#[derive(Debug, Clone)]
pub struct TurbulentDiffusionSph {
    /// Smagorinsky constant.
    pub cs: f64,
    /// Particle spacing / filter width.
    pub delta: f64,
}
impl TurbulentDiffusionSph {
    /// Create a new turbulent diffusion SPH operator.
    pub fn new(cs: f64, delta: f64) -> Self {
        Self { cs, delta }
    }
    /// Compute turbulent acceleration contribution on particle i from neighbor j.
    ///
    /// a_turb_i = Σ_j (m_j / ρ_j) · (τ_ij_i/ρ_i + τ_ij_j/ρ_j) · ∇W_ij
    ///
    /// Returns the 3-component acceleration vector.
    pub fn turbulent_acceleration(
        &self,
        s_i: Mat3,
        s_j: Mat3,
        rho_i: f64,
        rho_j: f64,
        mass_j: f64,
        grad_w: [f64; 3],
    ) -> [f64; 3] {
        let tau_i = sgs_stress_smagorinsky(s_i, rho_i, self.cs, self.delta);
        let tau_j = sgs_stress_smagorinsky(s_j, rho_j, self.cs, self.delta);
        let mut acc = [0.0_f64; 3];
        for alpha in 0..3 {
            let mut sum = 0.0_f64;
            for beta in 0..3 {
                sum += (tau_i[alpha][beta] / (rho_i * rho_i)
                    + tau_j[alpha][beta] / (rho_j * rho_j))
                    * grad_w[beta];
            }
            acc[alpha] = mass_j * sum;
        }
        acc
    }
}
/// Dynamic Smagorinsky LES model.
///
/// Uses the Germano identity to compute the Smagorinsky constant locally:
///
/// ```text
/// C_s² = < L_ij M_ij > / < M_ij M_ij >
/// ```
///
/// where L_ij is the resolved turbulent stress (Leonard stress) and M_ij is
/// the Germano tensor. A test filter at scale Δ̂ = 2Δ is used.
///
/// Reference: Germano et al. (1991), "A dynamic subgrid-scale eddy viscosity
/// model", Physics of Fluids A.
#[derive(Debug, Clone)]
pub struct DynamicSmagorinskyLes {
    /// Base filter width Δ (particle spacing, m).
    pub delta: f64,
    /// Test filter ratio α = Δ̂/Δ (typically 2).
    pub test_filter_ratio: f64,
    /// Clipping bounds for C_s² to avoid negative viscosity.
    pub cs_sq_min: f64,
    /// Maximum clipping bound for C_s².
    pub cs_sq_max: f64,
}
impl DynamicSmagorinskyLes {
    /// Create a new dynamic Smagorinsky model.
    pub fn new(delta: f64, test_filter_ratio: f64) -> Self {
        Self {
            delta,
            test_filter_ratio,
            cs_sq_min: 0.0,
            cs_sq_max: 0.04,
        }
    }
    /// Compute the dynamic C_s² using averaged Germano identity.
    ///
    /// # Arguments
    /// - `s_bar`     : strain-rate tensor at the grid filter level S̄_ij.
    /// - `s_hat`     : strain-rate tensor at the test-filter level Ŝ_ij.
    /// - `u_bar`     : velocity at grid-filter level (3-vector).
    /// - `u_hat`     : velocity at test-filter level (3-vector).
    ///
    /// # Returns
    /// Clipped value of C_s².
    pub fn compute_cs_sq(&self, s_bar: &Mat3, s_hat: &Mat3) -> f64 {
        let s_bar_mag = strain_rate_magnitude(s_bar);
        let s_hat_mag = strain_rate_magnitude(s_hat);
        let delta_hat = self.test_filter_ratio * self.delta;
        let lhs_scale = delta_hat * delta_hat * s_hat_mag;
        let rhs_scale = self.delta * self.delta * s_bar_mag;
        let mut m_ij = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                m_ij[i][j] = 2.0 * (lhs_scale * s_hat[i][j] - rhs_scale * s_bar[i][j]);
            }
        }
        let l_scale = (delta_hat * delta_hat - self.delta * self.delta) * s_bar_mag;
        let mut l_ij = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                l_ij[i][j] = l_scale * s_bar[i][j];
            }
        }
        let lm = mat3_double_contraction(&l_ij, &m_ij);
        let mm = mat3_double_contraction(&m_ij, &m_ij);
        if mm < 1e-30 {
            return 0.0;
        }
        (lm / mm).clamp(self.cs_sq_min, self.cs_sq_max)
    }
    /// Compute the local eddy viscosity using the dynamic C_s².
    pub fn eddy_viscosity(&self, s_bar: &Mat3, s_hat: &Mat3) -> f64 {
        let cs_sq = self.compute_cs_sq(s_bar, s_hat);
        let s_mag = strain_rate_magnitude(s_bar);
        self.delta * self.delta * cs_sq * s_mag
    }
    /// Return the dynamic SGS stress tensor τ_ij = -2 ρ ν_t S̄_ij.
    pub fn sgs_stress(&self, s_bar: &Mat3, s_hat: &Mat3, density: f64) -> Mat3 {
        let nu_t = self.eddy_viscosity(s_bar, s_hat);
        mat3_scale(-2.0 * density * nu_t, *s_bar)
    }
}
/// Stores a per-particle turbulent viscosity field and provides operations
/// such as SPH-smoothed diffusion and local averaging.
#[derive(Debug, Clone)]
pub struct TurbulentViscosityField {
    /// Per-particle turbulent kinematic viscosity ν_t (m² s⁻¹).
    pub nu_t: Vec<f64>,
    /// Per-particle turbulent kinetic energy k (m² s⁻²).
    pub k: Vec<f64>,
    /// Per-particle specific dissipation rate ω (s⁻¹) or ε (m² s⁻³).
    pub omega: Vec<f64>,
}
impl TurbulentViscosityField {
    /// Allocate a zero-initialised field for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            nu_t: vec![0.0; n],
            k: vec![0.0; n],
            omega: vec![0.0; n],
        }
    }
    /// Update ν_t from the stored k and ω (k-ω closure: ν_t = k/ω).
    pub fn update_from_k_omega(&mut self) {
        for i in 0..self.nu_t.len() {
            self.nu_t[i] = if self.omega[i] > 1e-30 {
                self.k[i].max(0.0) / self.omega[i]
            } else {
                0.0
            };
        }
    }
    /// Kernel-weighted average of ν_t at a point using SPH neighbours.
    ///
    /// ⟨ν_t⟩(x) = Σ_j (m_j/ρ_j) ν_t_j W(|x-x_j|, h)
    pub fn interpolate_nu_t(
        &self,
        pos: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_masses: &[f64],
        sph_densities: &[f64],
        h: f64,
    ) -> f64 {
        let mut num = 0.0_f64;
        for (j, &x_j) in sph_positions.iter().enumerate() {
            let r = {
                let dx = pos[0] - x_j[0];
                let dy = pos[1] - x_j[1];
                let dz = pos[2] - x_j[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            if r > 2.0 * h {
                continue;
            }
            let q = r / h;
            let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
            let w = if q < 1.0 {
                sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
            } else if q < 2.0 {
                let t = 2.0 - q;
                sigma * 0.25 * t * t * t
            } else {
                0.0
            };
            let vol_j = sph_masses[j] / sph_densities[j].max(1e-30);
            num += self.nu_t[j] * w * vol_j;
        }
        num
    }
    /// Clamp all ν_t values to the interval \[0, nu_t_max\].
    pub fn clamp_nu_t(&mut self, nu_t_max: f64) {
        for v in self.nu_t.iter_mut() {
            *v = v.clamp(0.0, nu_t_max);
        }
    }
    /// Compute the average turbulent kinetic energy across all particles.
    pub fn mean_k(&self) -> f64 {
        if self.k.is_empty() {
            return 0.0;
        }
        self.k.iter().sum::<f64>() / self.k.len() as f64
    }
}
/// Sigma turbulence model.
///
/// Uses the singular values of the velocity-gradient tensor to build
/// an eddy viscosity that vanishes at walls and in laminar regions.
/// Nicoud et al. (2011) "Using singular values to build a subgrid-scale
/// model for large eddy simulations."
#[derive(Debug, Clone)]
pub struct SigmaModel {
    /// Sigma model constant C_σ (default 1.35).
    pub c_sigma: f64,
    /// Filter width Δ.
    pub delta: f64,
}
impl SigmaModel {
    /// Create a new Sigma model.
    pub fn new(c_sigma: f64, delta: f64) -> Self {
        Self { c_sigma, delta }
    }
    /// Compute the three invariants of G = gᵀg (where g = velocity-gradient tensor).
    ///
    /// Returns (I1, I2, I3) of the matrix G = gᵀg.
    fn invariants_g(g: &Mat3) -> (f64, f64, f64) {
        let mut big_g = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                for g_k in g.iter() {
                    big_g[i][j] += g_k[i] * g_k[j];
                }
            }
        }
        let i1 = big_g[0][0] + big_g[1][1] + big_g[2][2];
        let i2 = 0.5
            * (i1 * i1
                - (big_g[0][0] * big_g[0][0]
                    + big_g[1][1] * big_g[1][1]
                    + big_g[2][2] * big_g[2][2]
                    + 2.0 * big_g[0][1] * big_g[1][0]
                    + 2.0 * big_g[0][2] * big_g[2][0]
                    + 2.0 * big_g[1][2] * big_g[2][1]));
        let i3 = big_g[0][0] * (big_g[1][1] * big_g[2][2] - big_g[1][2] * big_g[2][1])
            - big_g[0][1] * (big_g[1][0] * big_g[2][2] - big_g[1][2] * big_g[2][0])
            + big_g[0][2] * (big_g[1][0] * big_g[2][1] - big_g[1][1] * big_g[2][0]);
        (i1, i2, i3)
    }
    /// Approximate singular values σ_1 ≥ σ_2 ≥ σ_3 ≥ 0 of g via
    /// the eigenvalues of G = gᵀg (σ_i = sqrt(λ_i)).
    ///
    /// Uses the closed-form cubic formula for the eigenvalues of the
    /// symmetric 3×3 matrix G.
    fn singular_values(g: &Mat3) -> [f64; 3] {
        let (i1, i2, i3) = Self::invariants_g(g);
        let p = i2 - i1 * i1 / 3.0;
        let q = (2.0 * i1 * i1 * i1 / 27.0) - (i1 * i2 / 3.0) + i3;
        let disc = -(4.0 * p * p * p + 27.0 * q * q);
        let mut lambdas = [0.0_f64; 3];
        if disc >= 0.0 {
            let m = 2.0 * (-p / 3.0).sqrt();
            let theta = (3.0 * q / (p * m)).clamp(-1.0, 1.0).acos() / 3.0;
            use std::f64::consts::PI;
            lambdas[0] = m * theta.cos() + i1 / 3.0;
            lambdas[1] = m * (theta - 2.0 * PI / 3.0).cos() + i1 / 3.0;
            lambdas[2] = m * (theta - 4.0 * PI / 3.0).cos() + i1 / 3.0;
        } else {
            let a = (-q / 2.0 + (q * q / 4.0 + p * p * p / 27.0).sqrt()).cbrt();
            let b = if a.abs() > 1e-30 { -p / (3.0 * a) } else { 0.0 };
            lambdas[0] = a + b + i1 / 3.0;
        }
        lambdas.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        [
            lambdas[0].max(0.0).sqrt(),
            lambdas[1].max(0.0).sqrt(),
            lambdas[2].max(0.0).sqrt(),
        ]
    }
    /// Sigma model eddy viscosity:
    /// ν_t = (C_σ Δ)² · D_σ(σ_1, σ_2, σ_3)
    /// where D_σ = σ_3 (σ_1 - σ_2)(σ_2 - σ_3) / σ_1²
    pub fn eddy_viscosity(&self, g: Mat3) -> f64 {
        let [s1, s2, s3] = Self::singular_values(&g);
        if s1 < 1e-30 {
            return 0.0;
        }
        let d_sigma = s3 * (s1 - s2) * (s2 - s3) / (s1 * s1);
        let cw2 = self.c_sigma * self.delta;
        cw2 * cw2 * d_sigma.max(0.0)
    }
}
/// Simplified DES (Detached Eddy Simulation) model.
///
/// Blends between a RANS-like model near walls and LES in the bulk.
/// The blending uses a length-scale comparison: l_RANS = √k / (β* ω),
/// l_LES = C_DES * Δ.  Wherever l_LES < l_RANS the LES branch is active.
#[derive(Debug, Clone)]
pub struct DesModel {
    /// DES constant (default 0.65).
    pub c_des: f64,
    /// Grid spacing / particle spacing Δ.
    pub delta: f64,
    /// Inner k-ω model.
    pub k_omega: KOmegaModel,
}
impl DesModel {
    /// Create a new DES model.
    pub fn new(delta: f64, nu_mol: f64) -> Self {
        Self {
            c_des: 0.65,
            delta,
            k_omega: KOmegaModel::new(nu_mol),
        }
    }
    /// RANS length scale: l_RANS = √k / (β* ω).
    pub fn rans_length(&self, k: f64, omega: f64) -> f64 {
        if omega.abs() < 1e-14 {
            return f64::MAX;
        }
        k.max(0.0).sqrt() / (self.k_omega.beta_star * omega)
    }
    /// LES length scale: l_LES = C_DES * Δ.
    pub fn les_length(&self) -> f64 {
        self.c_des * self.delta
    }
    /// Effective length scale: min(l_RANS, l_LES).
    pub fn effective_length(&self, k: f64, omega: f64) -> f64 {
        self.rans_length(k, omega).min(self.les_length())
    }
    /// Effective eddy viscosity using the DES blending.
    ///
    /// ν_t = min(k/ω, (C_DES Δ)² |S|)
    pub fn eddy_viscosity_des(&self, k: f64, omega: f64, strain_mag: f64) -> f64 {
        let nu_rans = self.k_omega.eddy_viscosity(k, omega);
        let nu_les = (self.c_des * self.delta).powi(2) * strain_mag;
        nu_rans.min(nu_les).max(0.0)
    }
    /// Advance k and ω with DES-modified dissipation.
    pub fn step_des(&self, k: f64, omega: f64, strain_mag: f64, dt: f64) -> (f64, f64) {
        let l_eff = self.effective_length(k, omega);
        let l_rans = self.rans_length(k, omega);
        let f_des = if l_rans.abs() < 1e-14 {
            1.0
        } else {
            (l_rans / l_eff).max(1.0)
        };
        let nu_t = self.k_omega.eddy_viscosity(k, omega);
        let pk = self.k_omega.production_k(nu_t, strain_mag);
        let dk = self.k_omega.beta_star * k * omega * f_des;
        let po = self.k_omega.production_omega(strain_mag);
        let do_ = self.k_omega.dissipation_omega(omega);
        let k_new = (k + dt * (pk - dk)).max(0.0);
        let omega_new = (omega + dt * (po - do_)).max(1e-14);
        (k_new, omega_new)
    }
}
/// DES (Detached Eddy Simulation) detachment sensor.
///
/// Blends between RANS (near walls) and LES (bulk) using the wall distance
/// `d_wall` and the local grid/particle spacing `h`.
/// The blending function is: f_DES = min(1, d_wall / (C_DES · Δ)).
#[derive(Debug, Clone)]
pub struct DesDetachment {
    /// DES constant C_DES (typically 0.65).
    pub c_des: f64,
    /// Maximum LES length scale Δ (particle spacing or filter width).
    pub max_le: f64,
}
impl DesDetachment {
    /// Create a new DES detachment sensor.
    pub fn new(c_des: f64, max_le: f64) -> Self {
        Self { c_des, max_le }
    }
    /// DES blending function.
    ///
    /// Returns a value in `[0, 1]`:
    /// - near 0 → LES region (d_wall is small compared to C_DES · Δ)
    /// - near 1 → RANS region (d_wall >> C_DES · Δ)
    pub fn rans_les_blend(&self, d_wall: f64, h: f64) -> f64 {
        let les_len = self.c_des * h.max(self.max_le);
        if les_len < 1e-14 {
            return 1.0;
        }
        (d_wall / les_len).clamp(0.0, 1.0)
    }
}
/// k-ω RANS turbulence model for SPH simulations.
///
/// Transport equations:
/// ```text
/// Dk/Dt   = P_k - β* k ω + ∇·[(ν + σ_k ν_t) ∇k]
/// Dω/Dt   = α (ω/k) P_k - β ω² + ∇·[(ν + σ_ω ν_t) ∇ω]
/// ν_t     = k / ω
/// P_k     = ν_t |S|²
/// ```
///
/// Reference: Wilcox (1988), "Reassessment of the scale-determining equation
/// for advanced turbulence models", AIAA Journal.
#[derive(Debug, Clone)]
pub struct KOmegaRans {
    /// Closure constant α (≈ 5/9).
    pub alpha: f64,
    /// Closure constant β (≈ 3/40).
    pub beta: f64,
    /// Closure constant β* (≈ 9/100).
    pub beta_star: f64,
    /// Diffusion coefficient σ_k (≈ 0.5).
    pub sigma_k: f64,
    /// Diffusion coefficient σ_ω (≈ 0.5).
    pub sigma_omega: f64,
    /// Kinematic viscosity ν (m² s⁻¹).
    pub nu: f64,
}
impl KOmegaRans {
    /// Create a new k-ω model with custom constants.
    pub fn new(
        alpha: f64,
        beta: f64,
        beta_star: f64,
        sigma_k: f64,
        sigma_omega: f64,
        nu: f64,
    ) -> Self {
        Self {
            alpha,
            beta,
            beta_star,
            sigma_k,
            sigma_omega,
            nu,
        }
    }
    /// Turbulent kinematic viscosity: ν_t = k / ω.
    pub fn turbulent_viscosity(&self, k: f64, omega: f64) -> f64 {
        if omega < 1e-30 {
            return 0.0;
        }
        (k / omega).max(0.0)
    }
    /// Production rate of k: P_k = ν_t |S|².
    pub fn production(&self, k: f64, omega: f64, strain_rate: &Mat3) -> f64 {
        let nu_t = self.turbulent_viscosity(k, omega);
        let s_mag_sq = 2.0 * mat3_double_contraction(strain_rate, strain_rate);
        nu_t * s_mag_sq
    }
    /// Destruction rate of k: D_k = β* k ω.
    pub fn destruction_k(&self, k: f64, omega: f64) -> f64 {
        self.beta_star * k.max(0.0) * omega.max(0.0)
    }
    /// Destruction rate of ω: D_ω = β ω².
    pub fn destruction_omega(&self, omega: f64) -> f64 {
        self.beta * omega.max(0.0) * omega.max(0.0)
    }
    /// Production rate of ω: P_ω = α (ω/k) P_k (if k > 0).
    pub fn production_omega(&self, k: f64, omega: f64, strain_rate: &Mat3) -> f64 {
        if k < 1e-30 {
            return 0.0;
        }
        let p_k = self.production(k, omega, strain_rate);
        self.alpha * (omega / k) * p_k
    }
    /// Source terms for k and ω at a single particle.
    ///
    /// Returns (Dk/Dt_src, Dω/Dt_src) = (P_k - D_k, P_ω - D_ω).
    pub fn source_terms(&self, k: f64, omega: f64, strain_rate: &Mat3) -> (f64, f64) {
        let p_k = self.production(k, omega, strain_rate);
        let d_k = self.destruction_k(k, omega);
        let p_omega = self.production_omega(k, omega, strain_rate);
        let d_omega = self.destruction_omega(omega);
        (p_k - d_k, p_omega - d_omega)
    }
    /// Integrate k and ω forward by dt using explicit Euler with positivity
    /// enforcement.
    ///
    /// Diffusion terms are not included here; they should be added via SPH
    /// gradient operators.
    pub fn integrate_explicit(
        &self,
        k: f64,
        omega: f64,
        strain_rate: &Mat3,
        dt: f64,
    ) -> (f64, f64) {
        let (src_k, src_omega) = self.source_terms(k, omega, strain_rate);
        let k_new = (k + dt * src_k).max(1e-14);
        let omega_new = (omega + dt * src_omega).max(1e-14);
        (k_new, omega_new)
    }
    /// Compute the effective (total) viscosity ν_eff = ν + ν_t.
    pub fn effective_viscosity(&self, k: f64, omega: f64) -> f64 {
        self.nu + self.turbulent_viscosity(k, omega)
    }
    /// Compute the turbulent length scale L = k^(1/2) / ω.
    pub fn turbulent_length_scale(&self, k: f64, omega: f64) -> f64 {
        if omega < 1e-30 {
            return 0.0;
        }
        k.max(0.0).sqrt() / omega
    }
}
/// Algebraic mixing-length turbulence model.
///
/// Eddy viscosity:
/// ```text
/// ν_t = l_m² |S|
/// ```
/// where `l_m` is the mixing length (m) prescribed by the user or a
/// geometry-based formula.
///
/// For boundary layers the Prandtl mixing-length hypothesis gives:
/// ```text
/// l_m = κ y (1 - exp(-y+ / A+))   (van Driest damping)
/// ```
///
/// Reference: Prandtl (1925).
#[derive(Debug, Clone)]
pub struct MixingLengthModel {
    /// Von Kármán constant κ (≈ 0.41).
    pub kappa: f64,
    /// Van Driest damping constant A+ (≈ 26).
    pub a_plus: f64,
    /// Maximum mixing length cap (m).  Set to f64::MAX to disable.
    pub l_max: f64,
}
impl MixingLengthModel {
    /// Create a new mixing-length model.
    pub fn new(kappa: f64, a_plus: f64, l_max: f64) -> Self {
        Self {
            kappa,
            a_plus,
            l_max,
        }
    }
    /// Compute mixing length with van Driest near-wall damping.
    ///
    /// - `y`      : wall-normal distance (m).
    /// - `y_plus` : dimensionless wall distance y+ = y u_τ / ν.
    pub fn mixing_length_van_driest(&self, y: f64, y_plus: f64) -> f64 {
        let damping = 1.0 - (-y_plus / self.a_plus).exp();
        (self.kappa * y * damping).min(self.l_max)
    }
    /// Compute mixing length for free shear flows (no wall damping).
    ///
    /// - `y` : cross-stream coordinate (m).
    /// - `b` : half-width of the shear layer (m).
    pub fn mixing_length_free_shear(&self, b: f64) -> f64 {
        (0.09 * b).min(self.l_max)
    }
    /// Eddy viscosity: ν_t = l_m² |S̄|.
    ///
    /// `l_m` is the mixing length (m) and `strain_rate` is the strain-rate
    /// tensor S̄_ij.
    pub fn eddy_viscosity(&self, l_m: f64, strain_rate: &Mat3) -> f64 {
        let s_mag = strain_rate_magnitude(strain_rate);
        l_m * l_m * s_mag
    }
    /// SGS stress tensor: τ_ij = -2 ρ ν_t S̄_ij.
    pub fn sgs_stress(&self, l_m: f64, strain_rate: &Mat3, density: f64) -> Mat3 {
        let nu_t = self.eddy_viscosity(l_m, strain_rate);
        mat3_scale(-2.0 * density * nu_t, *strain_rate)
    }
}
/// Near-wall turbulence parameters.
#[derive(Debug, Clone)]
pub struct WallTurbulence {
    /// von Kármán constant κ (default 0.41).
    pub kappa: f64,
    /// van Driest damping constant A⁺ (default 26.0).
    pub a_plus: f64,
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
}
impl WallTurbulence {
    /// Create a new near-wall turbulence model.
    pub fn new(kappa: f64, a_plus: f64, nu: f64) -> Self {
        Self { kappa, a_plus, nu }
    }
    /// van Driest damping factor D(y⁺) = 1 - exp(-y⁺/A⁺).
    pub fn van_driest_damping(&self, y_plus: f64) -> f64 {
        1.0 - (-y_plus / self.a_plus).exp()
    }
    /// Mixing-length with van Driest damping: l_m = κ y D(y⁺).
    pub fn mixing_length(&self, y: f64, u_tau: f64) -> f64 {
        let y_plus = y * u_tau / self.nu.max(1e-30);
        self.kappa * y * self.van_driest_damping(y_plus)
    }
    /// Wall eddy viscosity: ν_t = l_m² |S|.
    pub fn wall_eddy_viscosity(&self, y: f64, u_tau: f64, s_mag: f64) -> f64 {
        let lm = self.mixing_length(y, u_tau);
        lm * lm * s_mag
    }
    /// Law-of-the-wall: u⁺ = y⁺ in viscous sublayer (y⁺ < 11.3),
    /// u⁺ = (1/κ) ln(y⁺) + B otherwise (B = 5.2).
    pub fn u_plus(&self, y_plus: f64) -> f64 {
        const VISCOUS_LIMIT: f64 = 11.3;
        const B: f64 = 5.2;
        if y_plus < VISCOUS_LIMIT {
            y_plus
        } else {
            (1.0 / self.kappa) * y_plus.ln() + B
        }
    }
    /// Friction velocity u_τ from wall shear stress τ_w and density ρ.
    pub fn friction_velocity(tau_w: f64, density: f64) -> f64 {
        (tau_w.abs() / density.max(1e-30)).sqrt()
    }
    /// Wall shear stress from law of the wall (Newton iteration to find u_τ).
    ///
    /// Given the fluid velocity u at wall-normal distance y,
    /// iterates to find u_τ such that u/u_τ = u⁺(u_τ·y/ν).
    pub fn wall_shear_stress(&self, u: f64, y: f64, density: f64) -> f64 {
        let mut u_tau = (u * self.nu / y.max(1e-30)).sqrt().max(1e-10);
        for _ in 0..20 {
            let y_plus = u_tau * y / self.nu.max(1e-30);
            let u_plus_calc = self.u_plus(y_plus);
            if u_plus_calc < 1e-14 {
                break;
            }
            u_tau = u / u_plus_calc;
        }
        density * u_tau * u_tau
    }
}
/// WALE (Wall-Adapting Local Eddy-viscosity) sub-grid scale model.
///
/// Nicoud & Ducros (1999).  The WALE model captures the correct near-wall
/// behaviour (ν_t → O(y³)) without requiring explicit wall distance.
#[derive(Debug, Clone)]
pub struct WaleModel {
    /// WALE constant C_w (default 0.325).
    pub c_w: f64,
    /// Filter width Δ (particle spacing).
    pub delta: f64,
}
impl WaleModel {
    /// Create a new WALE model.
    pub fn new(c_w: f64, delta: f64) -> Self {
        Self { c_w, delta }
    }
    /// Compute the traceless symmetric part of the square of the
    /// velocity-gradient tensor g_ij = ∂u_i/∂x_j.
    ///
    /// S^d_ij = ½(g_ik g_kj + g_jk g_ki) - (1/3) δ_ij g_kk² / something
    /// Here we use the standard WALE S^d formula.
    pub fn sd_tensor(&self, g: Mat3) -> Mat3 {
        let mut g2 = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                for (k, &g_ik) in g[i].iter().enumerate() {
                    g2[i][j] += g_ik * g[k][j];
                }
            }
        }
        let tr_g2 = g2[0][0] + g2[1][1] + g2[2][2];
        let mut sd = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                sd[i][j] = 0.5 * (g2[i][j] + g2[j][i]);
                if i == j {
                    sd[i][j] -= tr_g2 / 3.0;
                }
            }
        }
        sd
    }
    /// Compute the WALE eddy viscosity ν_t.
    ///
    /// ν_t = (C_w Δ)² · (S^d_ij S^d_ij)^(3/2) / \[(S_ij S_ij)^(5/2) + (S^d_ij S^d_ij)^(5/4)\]
    pub fn eddy_viscosity(&self, s: Mat3, g: Mat3) -> f64 {
        let sd = self.sd_tensor(g);
        let mut sd_sq = 0.0_f64;
        for row in &sd {
            for &v in row {
                sd_sq += v * v;
            }
        }
        let mut s_sq = 0.0_f64;
        for row in &s {
            for &v in row {
                s_sq += v * v;
            }
        }
        let denom = s_sq.powf(2.5) + sd_sq.powf(1.25);
        if denom < 1e-30 {
            return 0.0;
        }
        let cw2 = self.c_w * self.delta;
        cw2 * cw2 * sd_sq.powf(1.5) / denom
    }
    /// SGS stress tensor: τ_ij = -2 ρ ν_t S_ij.
    pub fn sgs_stress(&self, s: Mat3, g: Mat3, density: f64) -> Mat3 {
        let nu_t = self.eddy_viscosity(s, g);
        let mut tau = mat3_zero();
        for i in 0..3 {
            for j in 0..3 {
                tau[i][j] = -2.0 * density * nu_t * s[i][j];
            }
        }
        tau
    }
}
/// Lightweight k-ω turbulence model state for a single SPH particle.
///
/// Tracks turbulent kinetic energy `k` and specific dissipation `ω`, and
/// exposes the standard Wilcox (1988) model constants.
#[derive(Debug, Clone)]
pub struct KomegaSph {
    /// Turbulent kinetic energy k (m²/s²).
    pub k: f64,
    /// Specific dissipation rate ω (1/s).
    pub omega: f64,
    /// Production coefficient α (Wilcox standard: 5/9).
    pub alpha: f64,
    /// Dissipation coefficient β (Wilcox standard: 3/40).
    pub beta: f64,
    /// Turbulent diffusion coefficient for k.
    pub sigma_k: f64,
    /// Turbulent diffusion coefficient for ω.
    pub sigma_w: f64,
}
impl KomegaSph {
    /// Create a new k-ω state with Wilcox (1988) standard constants.
    pub fn new(k: f64, omega: f64) -> Self {
        Self {
            k,
            omega,
            alpha: 5.0 / 9.0,
            beta: 3.0 / 40.0,
            sigma_k: 0.5,
            sigma_w: 0.5,
        }
    }
    /// Turbulent kinematic viscosity: ν_t = k / ω (clamped to ≥ 0).
    pub fn turbulent_viscosity(&self) -> f64 {
        if self.omega.abs() < 1e-14 {
            return 0.0;
        }
        (self.k / self.omega).max(0.0)
    }
    /// Production term: P = ν_t · |S|²
    pub fn production_term(&self, s_mag_sq: f64) -> f64 {
        self.turbulent_viscosity() * s_mag_sq
    }
}
/// Two-equation k-ε turbulence model adapted for SPH.
///
/// Follows the standard high-Reynolds-number k-ε model
/// (Launder & Spalding 1974) with SPH-adapted diffusion terms.
#[derive(Debug, Clone)]
pub struct KepsilonSph {
    /// Turbulent kinetic energy k (m²/s²).
    pub k: f64,
    /// Turbulent dissipation rate ε (m²/s³).
    pub epsilon: f64,
    /// Model constant C_μ (default 0.09).
    pub c_mu: f64,
    /// Production constant C_1ε (default 1.44).
    pub c1_eps: f64,
    /// Destruction constant C_2ε (default 1.92).
    pub c2_eps: f64,
    /// Turbulent Prandtl number for k (default 1.0).
    pub sigma_k: f64,
    /// Turbulent Prandtl number for ε (default 1.3).
    pub sigma_eps: f64,
}
impl KepsilonSph {
    /// Create a new k-ε SPH model.
    pub fn new(k: f64, epsilon: f64) -> Self {
        Self {
            k,
            epsilon,
            ..Default::default()
        }
    }
    /// Turbulent viscosity ν_t = C_μ k² / ε.
    pub fn turbulent_viscosity(&self) -> f64 {
        if self.epsilon < 1e-30 {
            return 0.0;
        }
        self.c_mu * self.k * self.k / self.epsilon
    }
    /// Turbulent length scale L = C_μ^(3/4) k^(3/2) / ε.
    pub fn turbulent_length_scale(&self) -> f64 {
        if self.epsilon < 1e-30 {
            return 0.0;
        }
        self.c_mu.powf(0.75) * self.k.powf(1.5) / self.epsilon
    }
    /// Turbulent time scale τ = k / ε.
    pub fn turbulent_time_scale(&self) -> f64 {
        if self.epsilon < 1e-30 {
            return 0.0;
        }
        self.k / self.epsilon
    }
    /// Production term P_k = ν_t |S|² (input: strain-rate magnitude squared).
    pub fn production(&self, s_mag_sq: f64) -> f64 {
        self.turbulent_viscosity() * s_mag_sq
    }
    /// Advance k and ε by one explicit time step.
    ///
    /// dkdt   = P_k - ε
    /// dεdt   = (C_1ε P_k - C_2ε ε) / τ
    pub fn step(&mut self, s_mag_sq: f64, dt: f64) {
        let p_k = self.production(s_mag_sq);
        let tau = self.turbulent_time_scale().max(1e-30);
        let dk = (p_k - self.epsilon) * dt;
        let de = (self.c1_eps * p_k - self.c2_eps * self.epsilon) / tau * dt;
        self.k = (self.k + dk).max(1e-30);
        self.epsilon = (self.epsilon + de).max(1e-30);
    }
    /// SGS stress tensor: τ_ij = -2 ρ ν_t S_ij + (2/3) ρ k δ_ij.
    pub fn sgs_stress(&self, s: Mat3, density: f64) -> Mat3 {
        let nu_t = self.turbulent_viscosity();
        let mut tau = mat3_zero();
        let tke_diag = (2.0 / 3.0) * density * self.k;
        for i in 0..3 {
            for j in 0..3 {
                tau[i][j] = -2.0 * density * nu_t * s[i][j];
                if i == j {
                    tau[i][j] += tke_diag;
                }
            }
        }
        tau
    }
}
/// Turbulent Prandtl number model for scalar transport.
///
/// The turbulent Prandtl number Pr_t relates the turbulent momentum and
/// scalar diffusivities:
///
/// ```text
/// α_t = ν_t / Pr_t
/// ```
///
/// where α_t is the turbulent thermal diffusivity.
#[derive(Debug, Clone)]
pub struct TurbulentPrandtlNumber {
    /// Turbulent Prandtl number (typically 0.85–0.9 for temperature).
    pub pr_t: f64,
}
impl TurbulentPrandtlNumber {
    /// Create with a given Pr_t.
    pub fn new(pr_t: f64) -> Self {
        Self { pr_t }
    }
    /// Turbulent thermal diffusivity: α_t = ν_t / Pr_t.
    pub fn thermal_diffusivity(&self, nu_t: f64) -> f64 {
        nu_t / self.pr_t.max(1e-30)
    }
    /// Turbulent heat flux: q_t = -ρ c_p α_t ∇T.
    pub fn turbulent_heat_flux(
        &self,
        nu_t: f64,
        density: f64,
        cp: f64,
        grad_t: [f64; 3],
    ) -> [f64; 3] {
        let alpha_t = self.thermal_diffusivity(nu_t);
        let coeff = -density * cp * alpha_t;
        [coeff * grad_t[0], coeff * grad_t[1], coeff * grad_t[2]]
    }
}
/// Two-equation k-ω turbulence model (Wilcox 1988) adapted for SPH.
///
/// Tracks turbulent kinetic energy (k) and specific dissipation rate (ω)
/// per particle and computes an eddy viscosity ν_t = k / ω.
#[derive(Debug, Clone)]
pub struct KOmegaModel {
    /// Model constant α (production coefficient, default 5/9).
    pub alpha: f64,
    /// Model constant β (dissipation coefficient for ω, default 3/40).
    pub beta: f64,
    /// Model constant β_star (dissipation coefficient for k, default 0.09).
    pub beta_star: f64,
    /// Model constant σ_k (diffusion of k, default 0.5).
    pub sigma_k: f64,
    /// Model constant σ_omega (diffusion of ω, default 0.5).
    pub sigma_omega: f64,
    /// Molecular kinematic viscosity ν (m²/s).
    pub nu_mol: f64,
}
impl KOmegaModel {
    /// Create a new k-ω model with standard constants.
    pub fn new(nu_mol: f64) -> Self {
        Self {
            alpha: 5.0 / 9.0,
            beta: 3.0 / 40.0,
            beta_star: 0.09,
            sigma_k: 0.5,
            sigma_omega: 0.5,
            nu_mol,
        }
    }
    /// Compute eddy viscosity ν_t = k / ω (clamped to non-negative).
    pub fn eddy_viscosity(&self, k: f64, omega: f64) -> f64 {
        if omega.abs() < 1e-14 {
            return 0.0;
        }
        (k / omega).max(0.0)
    }
    /// Compute production term P_k = ν_t * |S|².
    pub fn production_k(&self, nu_t: f64, strain_mag: f64) -> f64 {
        nu_t * strain_mag * strain_mag
    }
    /// Compute dissipation term for k: ε_k = β* k ω.
    pub fn dissipation_k(&self, k: f64, omega: f64) -> f64 {
        self.beta_star * k * omega
    }
    /// Compute production term for ω: P_ω = α * |S|².
    pub fn production_omega(&self, strain_mag: f64) -> f64 {
        self.alpha * strain_mag * strain_mag
    }
    /// Compute dissipation term for ω: ε_ω = β ω².
    pub fn dissipation_omega(&self, omega: f64) -> f64 {
        self.beta * omega * omega
    }
    /// Advance k and ω by one explicit Euler step.
    ///
    /// Returns `(k_new, omega_new)`.
    pub fn step_explicit(&self, k: f64, omega: f64, strain_mag: f64, dt: f64) -> (f64, f64) {
        let nu_t = self.eddy_viscosity(k, omega);
        let pk = self.production_k(nu_t, strain_mag);
        let dk = self.dissipation_k(k, omega);
        let po = self.production_omega(strain_mag);
        let do_ = self.dissipation_omega(omega);
        let k_new = (k + dt * (pk - dk)).max(0.0);
        let omega_new = (omega + dt * (po - do_)).max(1e-14);
        (k_new, omega_new)
    }
}
