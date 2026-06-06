//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Synthetic turbulence inflow generator using the random Fourier modes method.
///
/// Generates a divergence-free turbulent velocity field as a superposition of
/// random Fourier modes consistent with a prescribed energy spectrum (von Kármán).
///
/// Reference: Smirnov, Shi & Celik, ASME J. Fluids Eng. 123, 359 (2001).
#[derive(Debug, Clone)]
pub struct SyntheticTurbulenceInflow {
    /// Number of Fourier modes N.
    pub n_modes: usize,
    /// Turbulence intensity (u'/U).
    pub intensity: f64,
    /// Integral length scale L.
    pub length_scale: f64,
    /// Wavenumber array (N modes).
    pub wavenumbers: Vec<f64>,
    /// Phase angles (N modes, pre-computed).
    pub phases: Vec<f64>,
    /// Amplitude factors (N modes, from spectrum).
    pub amplitudes: Vec<f64>,
}
impl SyntheticTurbulenceInflow {
    /// Create a new synthetic turbulence generator with N modes.
    ///
    /// Uses a simple deterministic initialization based on mode index
    /// to avoid requiring the rand crate.
    pub fn new(n_modes: usize, intensity: f64, length_scale: f64) -> Self {
        let mut wavenumbers = Vec::with_capacity(n_modes);
        let mut phases = Vec::with_capacity(n_modes);
        let mut amplitudes = Vec::with_capacity(n_modes);
        let k_min = 2.0 * std::f64::consts::PI / (10.0 * length_scale);
        let k_max = 2.0 * std::f64::consts::PI / (0.1 * length_scale);
        for n in 0..n_modes {
            let t = (n + 1) as f64 / n_modes as f64;
            let k = k_min * (k_max / k_min).powf(t);
            let phase = 2.0 * std::f64::consts::PI * lcg_pseudo_random(n as u64 + 1);
            let amp = von_karman_energy_spectral_amplitude(k, length_scale, intensity);
            wavenumbers.push(k);
            phases.push(phase);
            amplitudes.push(amp);
        }
        Self {
            n_modes,
            intensity,
            length_scale,
            wavenumbers,
            phases,
            amplitudes,
        }
    }
    /// Generate the fluctuating velocity u' at position (x, y, z) and time t.
    ///
    /// `u'(x,t) = Σ_n 2 a_n cos(k_n · x + ω_n t + φ_n)`
    ///
    /// Returns \[u', v', w'\] velocity fluctuations.
    pub fn velocity_fluctuation(&self, x: f64, y: f64, z: f64, t: f64) -> [f64; 3] {
        let mut u = [0.0f64; 3];
        for n in 0..self.n_modes {
            let k = self.wavenumbers[n];
            let phi = self.phases[n];
            let amp = self.amplitudes[n];
            let omega_n = k / self.length_scale;
            let theta = lcg_pseudo_random(n as u64 + 100) * std::f64::consts::PI;
            let psi = lcg_pseudo_random(n as u64 + 200) * 2.0 * std::f64::consts::PI;
            let kx = k * theta.sin() * psi.cos();
            let ky = k * theta.sin() * psi.sin();
            let kz = k * theta.cos();
            let phase = kx * x + ky * y + kz * z + omega_n * t + phi;
            let cos_ph = phase.cos();
            u[0] += 2.0 * amp * cos_ph * (ky * kz).cos();
            u[1] += 2.0 * amp * cos_ph * (kz * kx).cos();
            u[2] += 2.0 * amp * cos_ph * (kx * ky).cos();
        }
        u
    }
    /// Return the mean turbulence intensity.
    pub fn turbulence_intensity(&self) -> f64 {
        self.intensity
    }
    /// Return the turbulent kinetic energy k = (3/2) (u' rms)².
    pub fn turbulent_ke(&self, u_ref: f64) -> f64 {
        let u_rms = self.intensity * u_ref;
        1.5 * u_rms * u_rms
    }
}
/// Hybrid RANS/LES blending approach.
///
/// Computes an effective eddy viscosity that blends between a RANS model
/// (k-ω) and an LES model (Smagorinsky) based on grid resolution.
///
/// The blending parameter `sigma` = 0 → pure RANS, 1 → pure LES.
#[derive(Debug, Clone)]
pub struct HybridRansLes {
    /// k-ω model state.
    pub k: f64,
    /// Specific dissipation rate.
    pub omega: f64,
    /// Smagorinsky constant.
    pub cs: f64,
    /// k-ω parameters.
    pub komega_params: KOmegaParams,
}
impl HybridRansLes {
    /// Create a new hybrid RANS/LES model.
    pub fn new(k: f64, omega: f64, cs: f64) -> Self {
        Self {
            k,
            omega,
            cs,
            komega_params: KOmegaParams::wilcox_1988(),
        }
    }
    /// Compute the RANS eddy viscosity: ν_t^RANS = k / ω.
    pub fn nu_rans(&self) -> f64 {
        turbulent_viscosity(self.k, self.omega)
    }
    /// Compute the LES (Smagorinsky) eddy viscosity: ν_t^LES = (Cs Δ)² |S|.
    pub fn nu_les(&self, delta: f64, strain_rate: f64) -> f64 {
        (self.cs * delta) * (self.cs * delta) * strain_rate
    }
    /// Compute the blending function sigma based on the Kolmogorov length scale.
    ///
    /// `sigma = exp(-k / (3 nu_t * omega))`
    ///
    /// Approaches 1 (LES) when turbulent length scale < filter width.
    pub fn blending_sigma(&self, delta: f64, strain_rate: f64) -> f64 {
        let nu_t_les = self.nu_les(delta, strain_rate);
        let nu_t_rans = self.nu_rans();
        if nu_t_les >= nu_t_rans { 1.0 } else { 0.0 }
    }
    /// Effective eddy viscosity using linear blending.
    ///
    /// `ν_t = σ * ν_t^LES + (1-σ) * ν_t^RANS`
    pub fn effective_nu(&self, delta: f64, strain_rate: f64, nu_base: f64) -> f64 {
        let sigma = self.blending_sigma(delta, strain_rate);
        let nu_t_les = self.nu_les(delta, strain_rate);
        let nu_t_rans = self.nu_rans();
        nu_base + sigma * nu_t_les + (1.0 - sigma) * nu_t_rans
    }
    /// Update RANS state with explicit Euler.
    pub fn update_rans(&mut self, strain_rate: f64, dt: f64) {
        let prod = k_production(self.nu_rans(), strain_rate * strain_rate);
        let dk = dk_dt(self.k, self.omega, prod, &self.komega_params) * dt;
        let dw = domega_dt(self.k, self.omega, prod, &self.komega_params) * dt;
        self.k = (self.k + dk).max(1e-14);
        self.omega = (self.omega + dw).max(1e-14);
    }
}
/// LES filter width models.
pub struct LesFilterWidth;
impl LesFilterWidth {
    /// Cubic root filter width for a 3D cell:
    ///
    /// `Delta = (dx * dy * dz)^(1/3)`
    pub fn cubic_root(dx: f64, dy: f64, dz: f64) -> f64 {
        (dx * dy * dz).cbrt()
    }
    /// Maximum filter width:
    ///
    /// `Delta = max(dx, dy, dz)`
    pub fn maximum(dx: f64, dy: f64, dz: f64) -> f64 {
        dx.max(dy).max(dz)
    }
    /// Uniform filter width for equal spacing:
    ///
    /// `Delta = dx`
    pub fn uniform(dx: f64) -> f64 {
        dx
    }
    /// Van Driest-damped filter width near walls:
    ///
    /// `Delta_vd = Delta * (1 - exp(-y_plus / A_plus))`
    pub fn van_driest_damped(delta: f64, y_plus: f64, a_plus: f64) -> f64 {
        delta * (1.0 - (-y_plus / a_plus).exp())
    }
}
/// Vreman (2004) sub-grid scale model for LES.
///
/// Reference: Vreman, A.W. (2004). An eddy-viscosity subgrid-scale model for
/// turbulent shear flow: Algebraic theory and applications. *Phys. Fluids*,
/// 16(10):3670–3681.
#[derive(Debug, Clone, Copy)]
pub struct VremanModel {
    /// Vreman constant (typically 2.5 * Cs² ≈ 0.07).
    pub alpha: f64,
}
impl VremanModel {
    /// Create a new Vreman model with the given constant.
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }
    /// Compute the Vreman sub-grid viscosity ν_t.
    ///
    /// # Arguments
    /// * `velocity_gradient` – 3×3 velocity gradient tensor g_{ij} = ∂u_i/∂x_j
    ///   stored as a flat \[9\] array in row-major order.
    /// * `dx` – grid spacing (lattice units)
    ///
    /// # Returns
    /// Sub-grid turbulent viscosity ν_t ≥ 0.
    pub fn vreman_nu_t(&self, velocity_gradient: [f64; 9], _dx: f64) -> f64 {
        let g = velocity_gradient;
        let mut beta = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for m in 0..3 {
                    sum += g[m * 3 + i] * g[m * 3 + j];
                }
                beta[i * 3 + j] = sum;
            }
        }
        let b11 = beta[0];
        let b12 = beta[1];
        let b13 = beta[2];
        let b22 = beta[4];
        let b23 = beta[5];
        let b33 = beta[8];
        let b_beta = b11 * b22 - b12 * b12 + b11 * b33 - b13 * b13 + b22 * b33 - b23 * b23;
        if b_beta <= 0.0 {
            return 0.0;
        }
        let alpha_sq: f64 = g.iter().map(|x| x * x).sum();
        if alpha_sq <= 0.0 {
            return 0.0;
        }
        self.alpha * (b_beta / alpha_sq).sqrt()
    }
}
/// Dynamic Smagorinsky model using the Germano identity.
pub struct DynamicSmagorinsky;
impl DynamicSmagorinsky {
    /// Compute the dynamic constant Cs_dyn from resolved strain-rate
    /// information at two filter levels.
    ///
    /// Uses the Lilly (1992) least-squares formulation:
    ///
    /// `Cs^2 = <L_ij * M_ij> / <M_ij * M_ij>`
    ///
    /// where `l_m` = sum(L_ij * M_ij) and `m_m` = sum(M_ij * M_ij).
    pub fn dynamic_constant(l_m: f64, m_m: f64) -> f64 {
        if m_m.abs() < 1e-30 {
            return 0.0;
        }
        let cs2 = l_m / m_m;
        cs2.max(0.0)
    }
    /// Compute the Leonard stress tensor component L_ij:
    ///
    /// `L_ij = <u_i * u_j> - `u_i` * `u_j`
    ///
    /// where angle brackets denote test filtering.
    pub fn leonard_stress(ui_uj_filtered: f64, ui_filtered: f64, uj_filtered: f64) -> f64 {
        ui_uj_filtered - ui_filtered * uj_filtered
    }
    /// Compute M_ij tensor component:
    ///
    /// `M_ij = 2 * (Delta_test^2 * |S_test| * S_ij_test - Delta^2 * <|S| * S_ij>)`
    pub fn m_tensor_component(
        delta_test_sq: f64,
        s_mag_test: f64,
        s_ij_test: f64,
        delta_sq: f64,
        s_mag_s_ij_filtered: f64,
    ) -> f64 {
        2.0 * (delta_test_sq * s_mag_test * s_ij_test - delta_sq * s_mag_s_ij_filtered)
    }
    /// Effective dynamic viscosity:
    ///
    /// `nu_sgs = Cs_dyn^2 * Delta^2 * |S|`
    pub fn dynamic_viscosity(cs_dyn_sq: f64, delta: f64, strain_rate_mag: f64) -> f64 {
        cs_dyn_sq * delta * delta * strain_rate_mag
    }
}
/// Mixed SGS model combining Smagorinsky and scale-similar parts.
pub struct MixedModel {
    /// Smagorinsky constant.
    pub cs: f64,
    /// Scale-similar coefficient.
    pub c_ss: f64,
}
impl MixedModel {
    /// Create a new mixed model.
    pub fn new(cs: f64, c_ss: f64) -> Self {
        Self { cs, c_ss }
    }
    /// Compute the total SGS viscosity:
    ///
    /// `nu_sgs = (Cs * Delta)^2 * |S|`
    ///
    /// The scale-similar part is added as a stress, not a viscosity.
    pub fn smagorinsky_viscosity(&self, delta: f64, strain_rate: f64) -> f64 {
        (self.cs * delta) * (self.cs * delta) * strain_rate
    }
    /// Total effective viscosity including molecular:
    ///
    /// `nu_eff = nu + nu_sgs`
    pub fn effective_viscosity(&self, nu_base: f64, delta: f64, strain_rate: f64) -> f64 {
        nu_base + self.smagorinsky_viscosity(delta, strain_rate)
    }
}
/// Parameters for the SST k-ω model (Menter 1994).
#[derive(Debug, Clone, Copy)]
pub struct SstParams {
    /// Beta* coefficient (dissipation in k equation).
    pub beta_star: f64,
    /// Sigma_k1 (k equation inner layer).
    pub sigma_k1: f64,
    /// Sigma_k2 (k equation outer layer).
    pub sigma_k2: f64,
    /// Sigma_omega1 (omega equation inner layer).
    pub sigma_omega1: f64,
    /// Sigma_omega2 (omega equation outer layer).
    pub sigma_omega2: f64,
    /// Beta1 (omega destruction inner).
    pub beta1: f64,
    /// Beta2 (omega destruction outer).
    pub beta2: f64,
    /// Alpha1 (omega production inner).
    pub alpha1: f64,
    /// Alpha2 (omega production outer).
    pub alpha2: f64,
}
impl SstParams {
    /// Standard Menter (1994) SST constants.
    pub fn menter_1994() -> Self {
        Self {
            beta_star: 0.09,
            sigma_k1: 0.85,
            sigma_k2: 1.0,
            sigma_omega1: 0.5,
            sigma_omega2: 0.856,
            beta1: 0.075,
            beta2: 0.0828,
            alpha1: 5.0 / 9.0,
            alpha2: 0.44,
        }
    }
    /// Blend SST inner (k-ω) and outer (k-ε) constants using blending function F1.
    ///
    /// `phi_blended = F1 * phi_inner + (1 - F1) * phi_outer`
    pub fn blended_beta(&self, f1: f64) -> f64 {
        f1 * self.beta1 + (1.0 - f1) * self.beta2
    }
    /// Blended alpha for omega production.
    pub fn blended_alpha(&self, f1: f64) -> f64 {
        f1 * self.alpha1 + (1.0 - f1) * self.alpha2
    }
    /// Blended sigma_k.
    pub fn blended_sigma_k(&self, f1: f64) -> f64 {
        f1 * self.sigma_k1 + (1.0 - f1) * self.sigma_k2
    }
}
/// Wall-Adapting Local Eddy-Viscosity (WALE) SGS model.
///
/// The WALE model uses both the strain rate and rotation rate through the
/// traceless symmetric part of the squared velocity gradient g²:
/// `ν_t = (Cw Δ)² (Sd:Sd)^(3/2) / ((S:S)^(5/2) + (Sd:Sd)^(5/4))`
///
/// It naturally yields zero eddy viscosity at the wall (∝ y³) without
/// requiring van Driest damping.
///
/// Reference: Nicoud & Ducros, Flow Turbulence Combust. 62, 183–200 (1999).
#[derive(Debug, Clone, Copy)]
pub struct WaleModel {
    /// WALE model constant Cw (typically 0.5).
    pub cw: f64,
    /// Filter width Δ.
    pub delta: f64,
}
impl WaleModel {
    /// Construct WALE model with default Cw = 0.5.
    pub fn new(delta: f64) -> Self {
        Self { cw: 0.5, delta }
    }
    /// Construct with custom Cw.
    pub fn with_constant(delta: f64, cw: f64) -> Self {
        Self { cw, delta }
    }
    /// Compute the WALE eddy viscosity from the velocity gradient tensor.
    ///
    /// # Arguments
    /// * `g` – velocity gradient tensor as 9-element row-major array
    ///   [∂u/∂x, ∂u/∂y, ∂u/∂z, ∂v/∂x, ∂v/∂y, ∂v/∂z, ∂w/∂x, ∂w/∂y, ∂w/∂z]
    pub fn eddy_viscosity(&self, g: [f64; 9]) -> f64 {
        let sd = self.compute_sd(g);
        let s = strain_rate_tensor(g);
        let sd_sq: f64 = sd.iter().map(|v| v * v).sum();
        let s_sq: f64 = s.iter().map(|v| v * v).sum();
        let numerator = sd_sq.powf(1.5);
        let denominator = s_sq.powf(2.5) + sd_sq.powf(1.25);
        let l = self.cw * self.delta;
        if denominator > 1e-20 {
            l * l * numerator / denominator
        } else {
            0.0
        }
    }
    /// Compute the traceless symmetric part of g² (the Sd tensor).
    pub fn compute_sd(&self, g: [f64; 9]) -> [f64; 9] {
        let mut g2 = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    g2[i * 3 + j] += g[i * 3 + k] * g[k * 3 + j];
                }
            }
        }
        let trace = g2[0] + g2[4] + g2[8];
        let third = trace / 3.0;
        let mut sd = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                let sym = 0.5 * (g2[i * 3 + j] + g2[j * 3 + i]);
                sd[i * 3 + j] = if i == j { sym - third } else { sym };
            }
        }
        sd
    }
    /// Compute the SGS stress magnitude |τ_sgs| = 2 ν_t |S|.
    pub fn sgs_stress_magnitude(&self, g: [f64; 9]) -> f64 {
        let nu_t = self.eddy_viscosity(g);
        let s = strain_rate_tensor(g);
        let s_mag = (s.iter().map(|v| v * v).sum::<f64>() * 2.0).sqrt();
        2.0 * nu_t * s_mag
    }
}
/// Sigma SGS turbulence model.
///
/// Based on the singular values σ₁ ≥ σ₂ ≥ σ₃ ≥ 0 of the velocity gradient tensor:
/// `ν_t = (Cσ Δ)² σ₃(σ₁ − σ₂)(σ₂ − σ₃) / σ₁²`
///
/// The sigma model has the correct near-wall behaviour (ν_t ∝ y³) and
/// vanishes for solid-body rotation and irrotational flow.
///
/// Reference: Nicoud et al., Phys. Fluids 23, 085106 (2011).
#[derive(Debug, Clone, Copy)]
pub struct SigmaModel {
    /// Sigma model constant Cσ (typically 1.35).
    pub c_sigma: f64,
    /// Filter width Δ.
    pub delta: f64,
}
impl SigmaModel {
    /// Construct Sigma model with default Cσ = 1.35.
    pub fn new(delta: f64) -> Self {
        Self {
            c_sigma: 1.35,
            delta,
        }
    }
    /// Compute eddy viscosity from velocity gradient tensor.
    ///
    /// Uses singular values of g_{ij} via G = gᵀg (characteristic polynomial).
    pub fn eddy_viscosity(&self, g: [f64; 9]) -> f64 {
        let (s1, s2, s3) = self.singular_values_g(g);
        if s1 < 1e-20 {
            return 0.0;
        }
        let numer = s3 * (s1 - s2) * (s2 - s3);
        let l = self.c_sigma * self.delta;
        l * l * numer / (s1 * s1)
    }
    /// Compute the three singular values of g (approximated via gᵀg eigenvalues).
    ///
    /// Returns (σ₁, σ₂, σ₃) sorted descending.
    pub fn singular_values_g(&self, g: [f64; 9]) -> (f64, f64, f64) {
        let mut gt_g = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    gt_g[i * 3 + j] += g[k * 3 + i] * g[k * 3 + j];
                }
            }
        }
        let mut ev = eigenvalues_3x3_symmetric(gt_g);
        ev.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        (
            ev[0].max(0.0).sqrt(),
            ev[1].max(0.0).sqrt(),
            ev[2].max(0.0).sqrt(),
        )
    }
}
/// Scale-similar SGS model (Bardina model).
///
/// The scale-similar part of the SGS stress is:
///
/// `tau_ij^ss = C_ss * (u_i_hat * u_j_hat - u_i_hat_hat * u_j_hat_hat)`
///
/// where hat denotes the grid filter and hat-hat the test filter.
pub struct ScaleSimilarModel {
    /// Model coefficient (typically ~1.0).
    pub c_ss: f64,
}
impl ScaleSimilarModel {
    /// Create a new scale-similar model.
    pub fn new(c_ss: f64) -> Self {
        Self { c_ss }
    }
    /// Compute scale-similar stress component.
    pub fn stress_component(&self, ui_hat_uj_hat: f64, ui_hat_hat: f64, uj_hat_hat: f64) -> f64 {
        self.c_ss * (ui_hat_uj_hat - ui_hat_hat * uj_hat_hat)
    }
}
/// Standard k-ε model constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KEpsilonParams {
    /// Eddy-viscosity coefficient c_μ.
    pub c_mu: f64,
    /// C1 coefficient in ε production term.
    pub c1: f64,
    /// C2 coefficient in ε destruction term.
    pub c2: f64,
    /// Diffusion coefficient for k equation.
    pub sigma_k: f64,
    /// Diffusion coefficient for ε equation.
    pub sigma_eps: f64,
}
impl KEpsilonParams {
    /// Standard k-ε constants: c_μ=0.09, C1=1.44, C2=1.92, σ_k=1.0, σ_ε=1.3.
    pub fn standard() -> Self {
        Self {
            c_mu: 0.09,
            c1: 1.44,
            c2: 1.92,
            sigma_k: 1.0,
            sigma_eps: 1.3,
        }
    }
}
/// Dynamic Smagorinsky model using Lilly's least-squares formulation.
///
/// The model constant Cs² is computed dynamically from the Germano identity
/// applied at two filter levels (grid Δ and test Δ̂ = 2Δ).  The least-squares
/// minimization (Lilly 1992) gives:
/// `Cs² = ⟨L:M⟩ / ⟨M:M⟩`
/// where L = T − t (Leonard stress minus SGS stress model at test level)
/// and M = Δ̂² |Ŝ| Ŝ − Δ² |S| Ŝ (moment tensor).
///
/// Reference: Lilly, Phys. Fluids A 4, 633 (1992).
#[derive(Debug, Clone, Copy)]
pub struct DynamicSmagorinskyLilly {
    /// Grid filter width Δ.
    pub delta: f64,
    /// Test filter ratio r = Δ̂/Δ (typically 2).
    pub filter_ratio: f64,
    /// Current dynamic model constant Cs² (updated each time step).
    pub cs2: f64,
    /// Maximum Cs² clamp to prevent instability.
    pub cs2_max: f64,
}
impl DynamicSmagorinskyLilly {
    /// Construct with filter width Δ, filter ratio r, and initial Cs² = 0.
    pub fn new(delta: f64, filter_ratio: f64) -> Self {
        Self {
            delta,
            filter_ratio,
            cs2: 0.0,
            cs2_max: 0.04,
        }
    }
    /// Update Cs² from accumulated L:M and M:M statistics.
    ///
    /// # Arguments
    /// * `lm_sum` – accumulated ⟨L:M⟩ (numerator)
    /// * `mm_sum` – accumulated ⟨M:M⟩ (denominator)
    pub fn update_cs2(&mut self, lm_sum: f64, mm_sum: f64) {
        if mm_sum.abs() > 1e-20 {
            self.cs2 = (lm_sum / mm_sum).clamp(0.0, self.cs2_max);
        } else {
            self.cs2 = 0.0;
        }
    }
    /// Compute the SGS stress from the current Cs².
    ///
    /// τ_sgs = −2 Cs² Δ² |S| S_{ij}
    ///
    /// # Arguments
    /// * `s_mag` – local strain-rate magnitude |S|
    /// * `delta` – optional override filter width (uses self.delta if 0)
    pub fn sgs_stress_magnitude(&self, s_mag: f64) -> f64 {
        2.0 * self.cs2 * self.delta * self.delta * s_mag * s_mag
    }
    /// Compute local eddy viscosity ν_t = Cs² Δ² |S|.
    pub fn eddy_viscosity(&self, s_mag: f64) -> f64 {
        self.cs2 * self.delta * self.delta * s_mag
    }
    /// Compute the Leonard stress L = T̃ − Res(τ_test).
    ///
    /// For a single scalar tensor component.
    pub fn leonard_stress_scalar(&self, uu_hat: f64, u_hat_sq: f64) -> f64 {
        uu_hat - u_hat_sq
    }
    /// Compute the M tensor component: Δ̂²|Ŝ|Ŝ − Δ²|S|Ŝ_hat.
    pub fn m_tensor_scalar(
        &self,
        s_hat_mag: f64,
        s_hat_ij: f64,
        s_mag_hat: f64,
        s_ij_hat: f64,
    ) -> f64 {
        let delta_hat = self.filter_ratio * self.delta;
        delta_hat * delta_hat * s_hat_mag * s_hat_ij
            - self.delta * self.delta * s_mag_hat * s_ij_hat
    }
}
/// Accumulator for turbulence statistics (mean, variance, Reynolds stress).
#[derive(Debug, Clone)]
pub struct TurbulenceStatistics {
    /// Number of samples accumulated.
    pub n_samples: u64,
    /// Running sum of ux.
    pub sum_ux: f64,
    /// Running sum of uy.
    pub sum_uy: f64,
    /// Running sum of ux^2.
    pub sum_ux2: f64,
    /// Running sum of uy^2.
    pub sum_uy2: f64,
    /// Running sum of ux*uy (Reynolds stress).
    pub sum_ux_uy: f64,
    /// Running sum of k.
    pub sum_k: f64,
}
impl TurbulenceStatistics {
    /// Create new empty statistics.
    pub fn new() -> Self {
        Self {
            n_samples: 0,
            sum_ux: 0.0,
            sum_uy: 0.0,
            sum_ux2: 0.0,
            sum_uy2: 0.0,
            sum_ux_uy: 0.0,
            sum_k: 0.0,
        }
    }
    /// Add a sample.
    pub fn add_sample(&mut self, ux: f64, uy: f64, k: f64) {
        self.n_samples += 1;
        self.sum_ux += ux;
        self.sum_uy += uy;
        self.sum_ux2 += ux * ux;
        self.sum_uy2 += uy * uy;
        self.sum_ux_uy += ux * uy;
        self.sum_k += k;
    }
    /// Mean x-velocity.
    pub fn mean_ux(&self) -> f64 {
        if self.n_samples == 0 {
            return 0.0;
        }
        self.sum_ux / self.n_samples as f64
    }
    /// Mean y-velocity.
    pub fn mean_uy(&self) -> f64 {
        if self.n_samples == 0 {
            return 0.0;
        }
        self.sum_uy / self.n_samples as f64
    }
    /// Variance of x-velocity (u'u' Reynolds stress).
    pub fn variance_ux(&self) -> f64 {
        if self.n_samples < 2 {
            return 0.0;
        }
        let n = self.n_samples as f64;
        (self.sum_ux2 / n) - (self.sum_ux / n).powi(2)
    }
    /// Variance of y-velocity.
    pub fn variance_uy(&self) -> f64 {
        if self.n_samples < 2 {
            return 0.0;
        }
        let n = self.n_samples as f64;
        (self.sum_uy2 / n) - (self.sum_uy / n).powi(2)
    }
    /// Reynolds shear stress <u'v'>.
    pub fn reynolds_stress_uv(&self) -> f64 {
        if self.n_samples < 2 {
            return 0.0;
        }
        let n = self.n_samples as f64;
        (self.sum_ux_uy / n) - (self.sum_ux / n) * (self.sum_uy / n)
    }
    /// Turbulent kinetic energy from resolved fluctuations:
    ///
    /// `TKE = 0.5 * (var_ux + var_uy)`
    pub fn resolved_tke(&self) -> f64 {
        0.5 * (self.variance_ux() + self.variance_uy())
    }
    /// Mean modeled TKE.
    pub fn mean_k(&self) -> f64 {
        if self.n_samples == 0 {
            return 0.0;
        }
        self.sum_k / self.n_samples as f64
    }
    /// Turbulence intensity: `I = sqrt(2/3 * TKE) / U_mean`.
    pub fn turbulence_intensity(&self) -> f64 {
        let u_mean = (self.mean_ux().powi(2) + self.mean_uy().powi(2)).sqrt();
        if u_mean < 1e-15 {
            return 0.0;
        }
        let tke = self.resolved_tke();
        (2.0 / 3.0 * tke).sqrt() / u_mean
    }
    /// Reset all statistics.
    pub fn reset(&mut self) {
        *self = Self::new();
    }
}
/// Cell-local turbulence state for the k-ω model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KOmegaCell {
    /// Turbulent kinetic energy k (m²/s²).
    pub k: f64,
    /// Specific dissipation rate ω (1/s).
    pub omega: f64,
}
/// Wilcox k-ω model constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KOmegaParams {
    /// Diffusion coefficient for k equation.
    pub sigma_k: f64,
    /// Diffusion coefficient for ω equation.
    pub sigma_omega: f64,
    /// β* coefficient in k destruction term.
    pub beta_star: f64,
    /// β coefficient in ω destruction term.
    pub beta: f64,
    /// α coefficient in ω production term.
    pub alpha: f64,
}
impl KOmegaParams {
    /// Standard Wilcox (1988) k-ω model constants.
    ///
    /// σ_k = 2, σ_ω = 2, β* = 0.09, β = 0.075, α = 5/9
    pub fn wilcox_1988() -> Self {
        Self {
            sigma_k: 2.0,
            sigma_omega: 2.0,
            beta_star: 0.09,
            beta: 0.075,
            alpha: 5.0 / 9.0,
        }
    }
}
/// LBM-adapted k-ω turbulence model with inline state.
///
/// Holds cell-local turbulent kinetic energy `k` and specific dissipation
/// rate `omega`, along with the Wilcox model constants.
#[derive(Debug, Clone, Copy)]
pub struct KOmegaLBM {
    /// Turbulent kinetic energy k.
    pub k: f64,
    /// Specific dissipation rate omega.
    pub omega: f64,
    /// Model constants.
    pub params: KOmegaParams,
}
impl KOmegaLBM {
    /// Create a new KOmegaLBM with the Wilcox (1988) constants.
    pub fn new(k: f64, omega: f64) -> Self {
        Self {
            k,
            omega,
            params: KOmegaParams::wilcox_1988(),
        }
    }
    /// Compute turbulence production: P_k = (k/ω) · S².
    ///
    /// Returns P_k ≥ 0.  Returns 0 if ω ≤ 0.
    pub fn production(&self, strain_rate: f64) -> f64 {
        if self.omega <= 0.0 {
            return 0.0;
        }
        let nu_t = self.k / self.omega;
        nu_t * strain_rate * strain_rate
    }
    /// Compute destruction term for k: D_k = β* · k · ω.
    pub fn destruction(&self) -> f64 {
        self.params.beta_star * self.k * self.omega
    }
    /// Forward Euler update of k and ω given local strain rate and timestep.
    pub fn update(&mut self, strain_rate: f64, dt: f64) {
        let prod = self.production(strain_rate);
        let dk = (prod - self.destruction()) * dt;
        let domega = if self.k > 0.0 {
            (self.params.alpha * (self.omega / self.k) * prod
                - self.params.beta * self.omega * self.omega)
                * dt
        } else {
            (-self.params.beta * self.omega * self.omega) * dt
        };
        self.k = (self.k + dk).max(1e-14);
        self.omega = (self.omega + domega).max(1e-14);
    }
}
/// 2-D field of k-ω turbulence cells.
pub struct KOmegaField {
    /// Number of cells in x direction.
    pub nx: usize,
    /// Number of cells in y direction.
    pub ny: usize,
    /// Flat array of per-cell turbulence state (row-major: idx = j*nx + i).
    pub cells: Vec<KOmegaCell>,
    /// Model constants.
    pub params: KOmegaParams,
}
impl KOmegaField {
    /// Create a uniform field initialised to (k0, omega0) everywhere.
    pub fn new(nx: usize, ny: usize, k0: f64, omega0: f64) -> Self {
        let cells = vec![
            KOmegaCell {
                k: k0,
                omega: omega0
            };
            nx * ny
        ];
        Self {
            nx,
            ny,
            cells,
            params: KOmegaParams::wilcox_1988(),
        }
    }
    /// Return a reference to the cell at grid position (i, j).
    pub fn get(&self, i: usize, j: usize) -> &KOmegaCell {
        &self.cells[j * self.nx + i]
    }
    /// Forward Euler step for k and ω at each cell.
    ///
    /// `strain_rates_sq` must have length `nx * ny`; element `j*nx + i` is
    /// the local S² value for cell (i, j).
    pub fn step(&mut self, dt: f64, strain_rates_sq: &[f64]) {
        let _nx = self.nx;
        for (cell, &s2) in self.cells.iter_mut().zip(strain_rates_sq.iter()) {
            let k_old = cell.k;
            let w_old = cell.omega;
            let nu_t = turbulent_viscosity(k_old, w_old);
            let prod = k_production(nu_t, s2);
            let dk = dk_dt(k_old, w_old, prod, &self.params) * dt;
            let dw = domega_dt(k_old, w_old, prod, &self.params) * dt;
            cell.k = (k_old + dk).max(1e-14);
            cell.omega = (w_old + dw).max(1e-14);
        }
    }
    /// Effective kinematic viscosity at cell (i, j): ν_eff = ν + ν_t.
    pub fn effective_nu(&self, base_nu: f64, i: usize, j: usize) -> f64 {
        let cell = self.get(i, j);
        let nu_t = turbulent_viscosity(cell.k, cell.omega);
        base_nu + nu_t
    }
    /// Effective relaxation time at cell (i, j): τ = 1 / (3·ν_eff + 0.5).
    ///
    /// This is the BGK formula: τ = 1/ω_lbm where ω_lbm = 1/(3ν + 0.5).
    pub fn effective_tau(&self, base_nu: f64, i: usize, j: usize) -> f64 {
        let nu_eff = self.effective_nu(base_nu, i, j);
        1.0 / (3.0 * nu_eff + 0.5)
    }
}
/// DES-SST model: switches from SST-RANS near walls to LES far from walls.
///
/// The DES length scale modification replaces the RANS length scale
/// `l_RANS = sqrt(k) / (beta* omega)` by `l_DES = min(l_RANS, C_DES * Delta)`.
pub struct DesSstModel {
    /// SST parameters.
    pub params: SstParams,
    /// DES constant (typically 0.61 for k-ω DES).
    pub c_des: f64,
}
impl DesSstModel {
    /// Create a DES-SST model with given DES constant.
    pub fn new(c_des: f64) -> Self {
        Self {
            params: SstParams::menter_1994(),
            c_des,
        }
    }
    /// Compute the DES destruction term for k:
    ///
    /// `D_k^DES = beta* * k * omega * max(1, l_RANS / (C_DES * Delta))`
    pub fn k_destruction_des(&self, k: f64, omega: f64, delta: f64) -> f64 {
        let l_rans = if omega > 0.0 {
            k.sqrt() / (self.params.beta_star * omega)
        } else {
            0.0
        };
        let l_des = self.c_des * delta;
        let psi = if l_des > 1e-30 {
            (l_rans / l_des).max(1.0)
        } else {
            1.0
        };
        self.params.beta_star * k * omega * psi
    }
    /// Effective viscosity using DES-SST blend.
    pub fn effective_nu(&self, k: f64, omega: f64, s_mag: f64, f2: f64, nu: f64) -> f64 {
        let nu_t = sst_eddy_viscosity(k, omega, s_mag, f2);
        nu + nu_t
    }
}
/// Full scale-similar (Bardina) SGS model with mixed approach.
///
/// The model decomposes the SGS stress into a scale-similar part and an
/// eddy-viscosity part:
///
/// `τ_ij = C_ss * τ_ss_ij - 2 * C_v * δ² * |S̃| * S̃_ij`
///
/// where τ_ss is the Bardina scale-similar stress.
#[derive(Debug, Clone)]
pub struct BardinalFullModel {
    /// Scale-similar coefficient (typically 1.0).
    pub c_ss: f64,
    /// Eddy-viscosity correction coefficient.
    pub c_v: f64,
}
impl BardinalFullModel {
    /// Create with default coefficients.
    pub fn new(c_ss: f64, c_v: f64) -> Self {
        Self { c_ss, c_v }
    }
    /// Compute the total SGS stress tensor component tau_ij.
    ///
    /// # Arguments
    /// * `u_hat_i_u_hat_j`    – product of hat-filtered velocities
    /// * `u_hat_hat_i`        – double-filtered velocity component i
    /// * `u_hat_hat_j`        – double-filtered velocity component j
    /// * `delta`              – filter width
    /// * `s_mag`              – |S| magnitude
    /// * `s_ij`               – strain rate component S_ij
    pub fn sgs_stress(
        &self,
        u_hat_i_u_hat_j: f64,
        u_hat_hat_i: f64,
        u_hat_hat_j: f64,
        delta: f64,
        s_mag: f64,
        s_ij: f64,
    ) -> f64 {
        let tau_ss = u_hat_i_u_hat_j - u_hat_hat_i * u_hat_hat_j;
        self.c_ss * tau_ss - 2.0 * self.c_v * delta * delta * s_mag * s_ij
    }
    /// Effective viscosity from eddy-viscosity part only.
    pub fn eddy_viscosity(&self, delta: f64, s_mag: f64) -> f64 {
        self.c_v * delta * delta * s_mag
    }
}
/// Anisotropic Minimum Dissipation (AMD) SGS model.
///
/// The AMD model computes the eddy viscosity by ensuring the resolved kinetic
/// energy dissipation exactly balances the SGS energy production:
/// `ν_t = max(0, −C (δᵢ ∂u / ∂xⱼ)² (∂u / ∂xⱼ)²) / ((∂u / ∂xⱼ)²)²`
///
/// where δᵢ are the mesh spacings and C is a model constant (typically 0.3).
///
/// Reference: Rozema et al., Phys. Fluids 27, 085107 (2015).
#[derive(Debug, Clone, Copy)]
pub struct AmdModel {
    /// AMD model constant C (typically 1/3 for isotropic turbulence).
    pub c: f64,
    /// Grid spacings [Δx, Δy, Δz].
    pub delta: [f64; 3],
}
impl AmdModel {
    /// Construct AMD model with default C = 1/3.
    pub fn new(delta: [f64; 3]) -> Self {
        Self {
            c: 1.0 / 3.0,
            delta,
        }
    }
    /// Construct with custom model constant.
    pub fn with_constant(delta: [f64; 3], c: f64) -> Self {
        Self { c, delta }
    }
    /// Compute the AMD eddy viscosity from the velocity gradient tensor.
    ///
    /// # Arguments
    /// * `g` – velocity gradient ∂u_i/∂x_j as 9-element row-major array
    pub fn eddy_viscosity(&self, g: [f64; 9]) -> f64 {
        let d = &self.delta;
        let mut numer = 0.0f64;
        for i in 0..3 {
            let mut sum_ij = 0.0f64;
            for j in 0..3 {
                sum_ij += g[j * 3 + i] * g[j * 3 + i];
            }
            numer -= d[i] * d[i] * sum_ij;
        }
        numer *= self.c;
        let s = strain_rate_tensor(g);
        let denom = 2.0 * s.iter().map(|v| v * v).sum::<f64>();
        if denom > 1e-20 && numer < 0.0 {
            0.0
        } else if denom > 1e-20 {
            numer / denom
        } else {
            0.0
        }
    }
    /// Compute the anisotropy factor: ratio of largest to smallest grid spacing.
    pub fn anisotropy_factor(&self) -> f64 {
        let max = self.delta.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min = self.delta.iter().cloned().fold(f64::INFINITY, f64::min);
        if min > 1e-20 { max / min } else { 1.0 }
    }
}
/// Large Eddy Simulation model combining Smagorinsky SGS with LBM relaxation.
#[derive(Debug, Clone, Copy)]
pub struct LargeEddySimulation {
    /// Smagorinsky constant Cs (typically 0.1–0.2).
    pub cs: f64,
}
impl LargeEddySimulation {
    /// Create a new LES model with the given Smagorinsky constant.
    pub fn new(cs: f64) -> Self {
        Self { cs }
    }
    /// Compute the effective viscosity: ν_eff = ν + ν_sgs.
    ///
    /// Smagorinsky SGS viscosity:
    /// ```text
    /// ν_sgs = (Cs * dx)² * |S|
    /// ```
    /// where `|S|` is the magnitude of the strain-rate tensor.
    ///
    /// # Arguments
    /// * `strain_rate` – magnitude of the resolved strain-rate tensor |S|
    /// * `nu_base`     – molecular kinematic viscosity
    /// * `dx`          – grid spacing
    pub fn effective_nu(&self, strain_rate: f64, nu_base: f64, dx: f64) -> f64 {
        let nu_sgs = (self.cs * dx) * (self.cs * dx) * strain_rate;
        nu_base + nu_sgs
    }
    /// Compute the effective LBM relaxation frequency from ν_eff.
    ///
    /// Uses the standard BGK relation: ω_eff = 1 / (3 * ν_eff + 0.5).
    pub fn effective_omega(&self, strain_rate: f64, nu_base: f64, dx: f64) -> f64 {
        let nu_eff = self.effective_nu(strain_rate, nu_base, dx);
        1.0 / (3.0 * nu_eff + 0.5)
    }
}
