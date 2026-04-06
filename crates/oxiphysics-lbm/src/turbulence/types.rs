//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Simple constant turbulent Prandtl number for heat/scalar transport.
///
/// Returns the standard LES value `Pr_t = 0.4` when using the WALE model,
/// or 0.9 for the Smagorinsky model.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TurbPrandtlPreset {
    /// Standard LES (WALE) value.
    WaleLes,
    /// Standard LES (Smagorinsky) value.
    SmagorinskyLes,
    /// RANS k-epsilon default.
    KEpsilonRans,
}
#[allow(dead_code)]
impl TurbPrandtlPreset {
    /// Return the preset turbulent Prandtl number.
    pub fn value(self) -> f64 {
        match self {
            TurbPrandtlPreset::WaleLes => 0.4,
            TurbPrandtlPreset::SmagorinskyLes => 0.9,
            TurbPrandtlPreset::KEpsilonRans => 0.85,
        }
    }
}
/// k-ω Shear Stress Transport (SST) model by Menter (1994).
///
/// Combines the k-ω model near walls (using F1 blending) with the k-ε model
/// in free-stream regions. Also enforces a production limiter to prevent
/// excessive build-up of turbulent kinetic energy in stagnation zones.
#[allow(dead_code)]
pub struct KOmegaSst {
    /// Kinematic viscosity ν (m²/s).
    pub nu: f64,
    /// Bradshaw constant a1 (typically 0.31).
    pub a1: f64,
    /// k-ω inner model constant β* (= Cμ = 0.09).
    pub beta_star: f64,
    /// σ_ω2 — outer-model cross-diffusion blending coefficient (0.856).
    pub sigma_w2: f64,
}
impl KOmegaSst {
    /// Create a new k-ω SST model.
    ///
    /// # Arguments
    /// - `nu`  — kinematic viscosity
    /// - `a1`  — Bradshaw limit constant (default 0.31)
    #[allow(dead_code)]
    pub fn new(nu: f64, a1: f64) -> Self {
        Self {
            nu,
            a1,
            beta_star: 0.09,
            sigma_w2: 0.856,
        }
    }
    /// Compute the SST blending function F1 or F2.
    ///
    /// When `use_f1 = true`, computes the F1 function (wall-proximity blending).
    /// When `use_f1 = false`, computes the F2 function (eddy-viscosity limiter).
    ///
    /// # Arguments
    /// - `k`              — turbulent kinetic energy
    /// - `omega`          — specific dissipation rate ω
    /// - `y`              — wall-normal distance
    /// - `dk_domega_dot`  — cross-diffusion term ∇k · ∇ω  (only used for F1)
    /// - `use_f1`         — `true` → F1, `false` → F2
    ///
    /// Returns a value in \[0, 1\].
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn compute_blending_function(
        &self,
        k: f64,
        omega: f64,
        y: f64,
        dk_domega_dot: f64,
        use_f1: bool,
    ) -> f64 {
        if use_f1 {
            sst_blending_f1(k, omega, y, self.nu, dk_domega_dot)
        } else {
            sst_blending_f2(k, omega, y, self.nu)
        }
    }
    /// Apply the SST production limiter.
    ///
    /// Menter's SST model clips production to prevent unbounded growth in
    /// stagnation regions:
    ///
    /// `P_k_lim = min(P_k, C_lim * β* * k * ω)`
    ///
    /// where C_lim = 10 and β* = 0.09.
    ///
    /// Negative production is set to zero (realizability).
    ///
    /// # Arguments
    /// - `p_k`   — raw production term P_k = ν_t * S²
    /// - `k`     — turbulent kinetic energy
    /// - `omega` — specific dissipation rate
    ///
    /// Returns the limited production.
    #[allow(dead_code)]
    pub fn compute_production_limiter(&self, p_k: f64, k: f64, omega: f64) -> f64 {
        const C_LIM: f64 = 10.0;
        let p_k_non_neg = p_k.max(0.0);
        let ceiling = C_LIM * self.beta_star * k.max(0.0) * omega.max(0.0);
        p_k_non_neg.min(ceiling)
    }
}
/// Helper for estimating subgrid-scale (SGS) quantities when bridging
/// Large-Eddy Simulation (LES) and Direct Numerical Simulation (DNS) grids.
///
/// Uses the Kolmogorov scaling law to estimate the kinetic energy carried
/// by subgrid scales between the LES filter width Δ_LES and the DNS
/// resolution Δ_DNS.
#[allow(dead_code)]
pub struct LesToDns {
    /// Kolmogorov constant C_K (typically ~1.5 for 3-D isotropic turbulence).
    pub c_k: f64,
    /// Smagorinsky constant Cs (used to define the effective filter width).
    pub cs: f64,
}
impl LesToDns {
    /// Create a new `LesToDns` bridge object.
    ///
    /// # Arguments
    /// - `c_k` — Kolmogorov constant (default ~1.5)
    /// - `cs`  — Smagorinsky constant (default ~0.1)
    #[allow(dead_code)]
    pub fn new(c_k: f64, cs: f64) -> Self {
        Self { c_k, cs }
    }
    /// Estimate the subgrid-scale kinetic energy `k_sgs`.
    ///
    /// Based on the Kolmogorov inertial-range spectrum integrated between
    /// the DNS wave-number k_DNS = π/Δ_DNS and the LES cut-off k_LES = π/Δ_LES:
    ///
    /// `k_sgs ≈ k_resolved * [ (Δ_LES / Δ_DNS)^(2/3) − 1 ]`
    ///
    /// This is the additional kinetic energy that would be resolved by a
    /// DNS grid but is modelled (filtered out) by the LES grid.
    ///
    /// # Arguments
    /// - `k_resolved`  — resolved (LES) turbulent kinetic energy
    /// - `delta_les`   — LES filter width (grid spacing)
    /// - `delta_dns`   — DNS grid spacing (must be ≤ `delta_les`)
    ///
    /// Returns `k_sgs ≥ 0`.
    #[allow(dead_code)]
    pub fn compute_subgrid_kinetic_energy(
        &self,
        k_resolved: f64,
        delta_les: f64,
        delta_dns: f64,
    ) -> f64 {
        if k_resolved <= 0.0 || delta_dns <= 0.0 || delta_les <= delta_dns {
            return 0.0;
        }
        let ratio = delta_les / delta_dns;
        let scale_factor = ratio.powf(2.0 / 3.0) - 1.0;
        (k_resolved * scale_factor).max(0.0)
    }
}
/// Simple k-epsilon turbulence model state for a single cell.
///
/// Stores the turbulent kinetic energy `k` and its dissipation rate `epsilon`.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct KEpsilonState {
    /// Turbulent kinetic energy.
    pub k: f64,
    /// Turbulent dissipation rate.
    pub epsilon: f64,
}
#[allow(dead_code)]
impl KEpsilonState {
    /// Create a new k-epsilon state.
    pub fn new(k: f64, epsilon: f64) -> Self {
        Self { k, epsilon }
    }
    /// Compute the eddy viscosity: `nu_t = C_mu * k^2 / epsilon`.
    pub fn eddy_viscosity(&self) -> f64 {
        const C_MU: f64 = 0.09;
        if self.epsilon.abs() < 1e-30 {
            return 0.0;
        }
        C_MU * self.k * self.k / self.epsilon
    }
    /// Compute the turbulent time scale: `tau_t = k / epsilon`.
    pub fn time_scale(&self) -> f64 {
        if self.epsilon.abs() < 1e-30 {
            return 0.0;
        }
        self.k / self.epsilon
    }
    /// Compute the turbulent length scale: `l_t = C_mu^{3/4} * k^{3/2} / epsilon`.
    pub fn length_scale(&self) -> f64 {
        const C_MU_34: f64 = 0.09_f64;
        let c_mu_34 = C_MU_34.powf(0.75);
        if self.epsilon.abs() < 1e-30 {
            return 0.0;
        }
        c_mu_34 * self.k.powf(1.5) / self.epsilon
    }
    /// Advance the k-epsilon state by one step using a simplified model.
    ///
    /// `production` = P_k (production of turbulent energy from mean strain),
    /// `dt` = time step.
    pub fn advance(&mut self, production: f64, dt: f64) {
        const C_E1: f64 = 1.44;
        const C_E2: f64 = 1.92;
        const SIGMA_K: f64 = 1.0;
        const SIGMA_E: f64 = 1.3;
        let _ = SIGMA_K;
        let _ = SIGMA_E;
        let dk = (production - self.epsilon) * dt;
        let ts = self.time_scale();
        let de = if ts > 1e-30 {
            (C_E1 * production - C_E2 * self.epsilon) / ts * dt
        } else {
            0.0
        };
        self.k = (self.k + dk).max(1e-20);
        self.epsilon = (self.epsilon + de).max(1e-20);
    }
}
/// Dynamic Smagorinsky model with locally adapted Cs.
///
/// The dynamic procedure computes the optimal `Cs²` at each point by
/// applying a test filter (here: simple box average over nearest neighbors)
/// and invoking the Germano identity.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DynamicSmagorinsky {
    /// Grid spacing (LBM units, typically 1.0).
    pub dx: f64,
    /// Test-filter width relative to grid spacing (typically 2).
    pub filter_ratio: f64,
}
impl DynamicSmagorinsky {
    /// Create a new dynamic Smagorinsky model.
    #[allow(dead_code)]
    pub fn new(dx: f64, filter_ratio: f64) -> Self {
        Self { dx, filter_ratio }
    }
    /// Estimate the local dynamic Cs² from the Germano identity approximation.
    ///
    /// Requires the grid-scale strain rate `s_bar` and the test-scale strain
    /// rate `s_hat` (computed by applying the test filter to the velocity field).
    ///
    /// `Cs² = <L_ij M_ij> / <M_ij M_ij>`
    ///
    /// Here we use a simplified scalar form:
    /// `Cs² ≈ (s_hat - s_bar) / (alpha * s_hat^2)` (clipped to \[0, 0.04\])
    ///
    /// where `alpha = 2 * (filter_ratio² - 1) * dx²`.
    #[allow(dead_code)]
    pub fn compute_cs_sq(&self, s_bar: f64, s_hat: f64) -> f64 {
        let alpha = 2.0 * (self.filter_ratio * self.filter_ratio - 1.0) * self.dx * self.dx;
        if s_hat.abs() < 1e-15 {
            return 0.0;
        }
        let cs_sq = (s_hat - s_bar) / (alpha * s_hat * s_hat);
        cs_sq.clamp(0.0, 0.04)
    }
    /// Compute the effective omega given grid-scale and test-scale strain rates.
    #[allow(dead_code)]
    pub fn effective_omega(&self, omega_base: f64, s_bar: f64, s_hat: f64) -> f64 {
        let cs_sq = self.compute_cs_sq(s_bar, s_hat);
        let cs_dyn = cs_sq.sqrt();
        let tau = 1.0 / omega_base;
        let mixing_length = cs_dyn * self.dx;
        let nu_t = compute_turbulent_viscosity(s_bar, mixing_length);
        let tau_eff = tau + 3.0 * nu_t;
        1.0 / tau_eff.max(0.5 + 1e-10)
    }
}
/// Smagorinsky sub-grid scale (SGS) turbulence model parameters.
///
/// The Smagorinsky model augments the molecular viscosity by a turbulent
/// viscosity `nu_t = (Cs * Delta)^2 * |S|`, where `Delta` is the filter
/// width (grid spacing) and `|S|` is the resolved strain-rate magnitude.
#[derive(Debug, Clone)]
pub struct SmagorinskyModel {
    /// Smagorinsky constant (typically 0.1 – 0.2).
    pub cs_smag: f64,
}
impl SmagorinskyModel {
    /// Create a new Smagorinsky model with the given constant.
    pub fn new(cs_smag: f64) -> Self {
        Self { cs_smag }
    }
    /// Create with the standard default constant `Cs = 0.1`.
    #[allow(dead_code)]
    pub fn default() -> Self {
        Self { cs_smag: 0.1 }
    }
}
