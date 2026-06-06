//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Localization band analysis for strain-softening damage.
///
/// Implements the localization condition (bifurcation) for a damage
/// model: checks whether the acoustic tensor becomes singular.
#[derive(Debug, Clone)]
pub struct DamageBandLocalization {
    /// Current damage variable D.
    pub d: f64,
    /// Softening modulus H_soft (negative for softening).
    pub h_soft: f64,
    /// Young's modulus E.
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
}
impl DamageBandLocalization {
    /// Create a new localization band analysis object.
    pub fn new(d: f64, h_soft: f64, e: f64, nu: f64) -> Self {
        Self { d, h_soft, e, nu }
    }
    /// Localization condition for uniaxial stress.
    ///
    /// Returns `true` if the tangent modulus becomes zero or negative,
    /// indicating onset of localization (strain bifurcation).
    ///
    /// Tangent E_t = (1−D)·E + dD/dε · (−E·ε)
    ///            ≈ (1−D)·E − H_soft
    pub fn localization_condition(&self) -> bool {
        let e_t = (1.0 - self.d) * self.e - self.h_soft.abs();
        e_t <= 0.0
    }
    /// Width of a process zone in a rate-dependent (viscous) regularization.
    ///
    /// w_pz = l_c · π / 2 where l_c = √(E · l²_visc / |H|).
    pub fn process_zone_width(&self, l_visc_sq: f64) -> f64 {
        use std::f64::consts::PI;
        let h_abs = self.h_soft.abs().max(1e-30);
        let lc = ((1.0 - self.d) * self.e * l_visc_sq / h_abs).sqrt();
        lc * PI / 2.0
    }
    /// Critical wave vector for localization:
    ///
    /// k_c = √(−H̃ / (2G))  where H̃ = effective tangent modulus.
    pub fn critical_wave_vector(&self) -> Option<f64> {
        let g = self.e / (2.0 * (1.0 + self.nu));
        let h_eff = (1.0 - self.d) * self.e - self.h_soft.abs();
        if h_eff >= 0.0 {
            return None;
        }
        Some((-h_eff / (2.0 * g)).sqrt())
    }
}
/// Visualization data for damage fields.
#[derive(Debug, Clone)]
pub struct DamageVisualization {
    /// Damage values at each point (0..1).
    pub damage_values: Vec<f64>,
    /// Color mapping: maps damage to RGB \[r, g, b\] in \[0, 1\].
    /// Blue = undamaged, red = fully damaged.
    pub colors: Vec<[f64; 3]>,
}
impl DamageVisualization {
    /// Create visualization data from a set of damage states.
    pub fn from_states(states: &[DamageState]) -> Self {
        let damage_values: Vec<f64> = states.iter().map(|s| s.d).collect();
        let colors: Vec<[f64; 3]> = damage_values
            .iter()
            .map(|&d| Self::damage_to_color(d))
            .collect();
        Self {
            damage_values,
            colors,
        }
    }
    /// Map a damage value \[0, 1\] to an RGB color.
    /// Blue (undamaged) -> yellow (partial) -> red (fully damaged).
    pub fn damage_to_color(d: f64) -> [f64; 3] {
        let d = d.clamp(0.0, 1.0);
        if d < 0.5 {
            let t = d * 2.0;
            [t, t, 1.0 - t]
        } else {
            let t = (d - 0.5) * 2.0;
            [1.0, 1.0 - t, 0.0]
        }
    }
    /// Return the minimum damage value.
    pub fn min_damage(&self) -> f64 {
        self.damage_values
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min)
    }
    /// Return the maximum damage value.
    pub fn max_damage(&self) -> f64 {
        self.damage_values.iter().copied().fold(0.0_f64, f64::max)
    }
    /// Return the average damage value.
    pub fn avg_damage(&self) -> f64 {
        if self.damage_values.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.damage_values.iter().sum();
        sum / self.damage_values.len() as f64
    }
}
/// Coupled creep-damage model (Kachanov-Rabotnov type).
///
/// The Kachanov-Rabotnov model describes creep damage via:
///
/// ```text
/// de_cr/dt = A * sigma^n / (1-D)^m
/// dD/dt = B * sigma^r / (1-D)^k
/// ```
///
/// where A, B, n, m, r, k are material constants.
#[derive(Debug, Clone)]
pub struct DamageCreepCoupling {
    /// Creep rate coefficient A (1/Pa^n/s).
    pub a_creep: f64,
    /// Creep stress exponent n.
    pub n_creep: f64,
    /// Creep damage coupling exponent m.
    pub m_creep: f64,
    /// Damage rate coefficient B (1/Pa^r/s).
    pub b_damage: f64,
    /// Damage stress exponent r.
    pub r_damage: f64,
    /// Damage acceleration exponent k.
    pub k_damage: f64,
}
impl DamageCreepCoupling {
    /// Create a new creep-damage model with Kachanov-Rabotnov parameters.
    pub fn new(
        a_creep: f64,
        n_creep: f64,
        m_creep: f64,
        b_damage: f64,
        r_damage: f64,
        k_damage: f64,
    ) -> Self {
        Self {
            a_creep,
            n_creep,
            m_creep,
            b_damage,
            r_damage,
            k_damage,
        }
    }
    /// Compute the creep strain rate de_cr/dt.
    ///
    /// de_cr/dt = A * |sigma|^n / (1-D)^m
    pub fn creep_strain_rate(&self, sigma: f64, d: f64) -> f64 {
        let omega = 1.0 - d.clamp(0.0, 1.0 - 1e-10);
        self.a_creep * sigma.abs().powf(self.n_creep) / omega.powf(self.m_creep)
    }
    /// Compute the damage rate dD/dt.
    ///
    /// dD/dt = B * |sigma|^r / (1-D)^k
    pub fn damage_rate(&self, sigma: f64, d: f64) -> f64 {
        let omega = 1.0 - d.clamp(0.0, 1.0 - 1e-10);
        self.b_damage * sigma.abs().powf(self.r_damage) / omega.powf(self.k_damage)
    }
    /// Estimate time to rupture using Kachanov's formula:
    ///
    /// t_r = 1 / (B * sigma^r * (k+1))
    pub fn time_to_rupture(&self, sigma: f64) -> f64 {
        let denom = self.b_damage * sigma.abs().powf(self.r_damage) * (self.k_damage + 1.0);
        if denom < 1e-60 {
            f64::INFINITY
        } else {
            1.0 / denom
        }
    }
    /// Perform a single explicit Euler step for creep-damage evolution.
    ///
    /// Returns updated `(epsilon_creep, D)` after time step `dt`.
    pub fn explicit_step(&self, sigma: f64, epsilon_creep: f64, d: f64, dt: f64) -> (f64, f64) {
        let deps_cr = self.creep_strain_rate(sigma, d) * dt;
        let dd = self.damage_rate(sigma, d) * dt;
        (epsilon_creep + deps_cr, (d + dd).clamp(0.0, 1.0))
    }
}
/// Gurson-Tvergaard-Needleman (GTN) model for void growth in ductile metals.
///
/// Yield surface:
/// φ = (σ_eq/σ_0)² + 2 q1 f cosh(1.5 q2 σ_h/σ_0) − (1 + q3 f²)
#[derive(Debug, Clone)]
pub struct GursonModel {
    /// Void volume fraction f ∈ \[0, 1).
    pub f: f64,
    /// GTN parameter q1 (typically ≈ 1.5).
    pub q1: f64,
    /// GTN parameter q2 (typically ≈ 1.0).
    pub q2: f64,
    /// GTN parameter q3 (typically = q1²).
    pub q3: f64,
}
impl GursonModel {
    /// Gurson yield function value.
    ///
    /// * `sigma_eq` – von-Mises equivalent stress
    /// * `sigma_h`  – hydrostatic (mean) stress
    /// * `sigma_0`  – current yield stress of the matrix
    pub fn yield_function(sigma_eq: f64, sigma_h: f64, sigma_0: f64) -> f64 {
        let _ = (sigma_eq, sigma_h, sigma_0);
        0.0
    }
    /// Gurson yield function using this model's parameters.
    pub fn phi(&self, sigma_eq: f64, sigma_h: f64, sigma_0: f64) -> f64 {
        if sigma_0.abs() < 1e-30 {
            return f64::NAN;
        }
        let ratio = sigma_eq / sigma_0;
        let cosh_term = (1.5 * self.q2 * sigma_h / sigma_0).cosh();
        ratio * ratio + 2.0 * self.q1 * self.f * cosh_term - (1.0 + self.q3 * self.f * self.f)
    }
    /// Void growth rate: df/dε_p^vol = (1 - f) * ε_p_vol
    ///
    /// * `eps_p_vol` – volumetric plastic strain rate
    pub fn void_growth_rate(&self, eps_p_vol: f64) -> f64 {
        (1.0 - self.f) * eps_p_vol
    }
}
/// Failure mode at an integration point.
#[derive(Debug, Clone, PartialEq)]
pub enum FailureMode {
    /// Material is undamaged (D approx 0).
    Intact,
    /// Tensile micro-cracking is active.
    DamagingTension {
        /// Current damage level.
        d: f64,
    },
    /// Compressive micro-crushing is active.
    DamagingCompression {
        /// Current damage level.
        d: f64,
    },
    /// Material has lost almost all load-carrying capacity (D -> 1).
    FullyDamaged,
}
impl FailureMode {
    /// Classify the failure mode given scalar damage D and tension fraction alpha_t.
    ///
    /// * D = 0                   -> `Intact`
    /// * D >= 0.99               -> `FullyDamaged`
    /// * alpha_t >= 0.5          -> `DamagingTension`
    /// * otherwise               -> `DamagingCompression`
    pub fn classify(d: f64, alpha_t: f64) -> FailureMode {
        if d <= 0.0 {
            FailureMode::Intact
        } else if d >= 0.99 {
            FailureMode::FullyDamaged
        } else if alpha_t >= 0.5 {
            FailureMode::DamagingTension { d }
        } else {
            FailureMode::DamagingCompression { d }
        }
    }
    /// Return true if the material has any damage.
    pub fn is_damaged(&self) -> bool {
        !matches!(self, FailureMode::Intact)
    }
    /// Return true if the material is fully damaged.
    pub fn is_fully_damaged(&self) -> bool {
        matches!(self, FailureMode::FullyDamaged)
    }
    /// Extract the damage value, or 0.0 for Intact.
    pub fn damage_value(&self) -> f64 {
        match self {
            FailureMode::Intact => 0.0,
            FailureMode::DamagingTension { d } => *d,
            FailureMode::DamagingCompression { d } => *d,
            FailureMode::FullyDamaged => 1.0,
        }
    }
}
/// Effective stiffness tensor for a damaged composite using the
/// Mori-Tanaka homogenization with damage.
///
/// Cracks are modeled as penny-shaped inclusions with zero stiffness,
/// and their effect on the effective moduli is captured via crack density.
#[derive(Debug, Clone)]
pub struct DamageHomogenization {
    /// Matrix Young's modulus E_m.
    pub e_m: f64,
    /// Matrix Poisson's ratio ν_m.
    pub nu_m: f64,
    /// Crack density parameter ρ = N a³ / V (penny-shaped cracks).
    pub crack_density: f64,
}
impl DamageHomogenization {
    /// Create a new damage homogenization model.
    pub fn new(e_m: f64, nu_m: f64, crack_density: f64) -> Self {
        Self {
            e_m,
            nu_m,
            crack_density,
        }
    }
    /// Effective Young's modulus via Budiansky-O'Connell penny-crack formula:
    ///
    /// E_eff/E_m = 1 − ρ · 16(1−ν²) / (9(1−ν/2))
    pub fn effective_youngs_modulus(&self) -> f64 {
        let nu = self.nu_m;
        let factor = 16.0 * (1.0 - nu * nu) / (9.0 * (1.0 - nu / 2.0));
        self.e_m * (1.0 - self.crack_density * factor).max(0.0)
    }
    /// Effective shear modulus via Budiansky-O'Connell:
    ///
    /// G_eff/G_m = 1 − ρ · 32(1−ν)(5−ν) / (45(2−ν))
    pub fn effective_shear_modulus(&self) -> f64 {
        let nu = self.nu_m;
        let g_m = self.e_m / (2.0 * (1.0 + nu));
        let factor = 32.0 * (1.0 - nu) * (5.0 - nu) / (45.0 * (2.0 - nu));
        g_m * (1.0 - self.crack_density * factor).max(0.0)
    }
    /// Effective bulk modulus from effective E and G.
    pub fn effective_bulk_modulus(&self) -> f64 {
        let e_eff = self.effective_youngs_modulus();
        let g_eff = self.effective_shear_modulus();
        if g_eff.abs() < 1e-30 {
            return 0.0;
        }
        e_eff * g_eff / (3.0 * (3.0 * g_eff - e_eff))
    }
    /// Convert crack density to a scalar damage variable:
    ///
    /// D = 1 − E_eff / E_m
    pub fn damage_from_crack_density(&self) -> f64 {
        1.0 - self.effective_youngs_modulus() / self.e_m.max(1e-30)
    }
}
/// Lemaitre damage model (simplified scalar version).
///
/// Damage evolution: dD/dε_p = (σ_eq / S)^s0 * (σ_h / σ_eq)^s0
#[derive(Debug, Clone)]
pub struct LemaitreCDM {
    /// Current damage variable D ∈ \[0, 1\].
    pub d: f64,
    /// Damage energy release rate denominator S.
    pub s: f64,
    /// Damage exponent s0.
    pub s0: f64,
}
impl LemaitreCDM {
    /// Damage evolution rate dD/dε_p.
    ///
    /// * `eps_p`    – equivalent plastic strain increment
    /// * `sigma_eq` – von-Mises equivalent stress
    /// * `sigma_h`  – hydrostatic (mean) stress
    pub fn evolution_rate(&self, eps_p: f64, sigma_eq: f64, sigma_h: f64) -> f64 {
        if self.s.abs() < 1e-30 || sigma_eq.abs() < 1e-30 {
            return 0.0;
        }
        let triax = sigma_h / sigma_eq;
        let y = sigma_eq * sigma_eq / (2.0 * self.s) * (2.0 / 3.0 + 3.0 * triax * triax);
        (y / self.s).powf(self.s0) * eps_p
    }
}
/// Finite element that tracks damage at every Gauss point using the
/// isotropic scalar model.
#[derive(Debug, Clone)]
pub struct DamageMechanicsElement {
    /// Number of Gauss (quadrature) points.
    pub n_gauss_points: usize,
    /// Per-point damage states.
    pub damage_states: Vec<DamageState>,
    /// Material damage model.
    pub material_damage: IsotropicDamage,
    /// Whether this element has been deleted (all points failed).
    pub deleted: bool,
    /// Element ID for tracking.
    pub element_id: usize,
}
impl DamageMechanicsElement {
    /// Construct a new element with `n_gauss` integration points.
    pub fn new(n_gauss: usize, mat: IsotropicDamage) -> Self {
        Self {
            n_gauss_points: n_gauss,
            damage_states: vec![DamageState::default(); n_gauss],
            material_damage: mat,
            deleted: false,
            element_id: 0,
        }
    }
    /// Construct a new element with an ID.
    pub fn with_id(n_gauss: usize, mat: IsotropicDamage, id: usize) -> Self {
        Self {
            n_gauss_points: n_gauss,
            damage_states: vec![DamageState::default(); n_gauss],
            material_damage: mat,
            deleted: false,
            element_id: id,
        }
    }
    /// Update every Gauss point given a slice of strain tensors (one per point).
    ///
    /// # Panics
    /// Panics if `strains.len() != self.n_gauss_points`.
    pub fn update_all(&mut self, strains: &[[f64; 6]]) {
        assert_eq!(
            strains.len(),
            self.n_gauss_points,
            "strains slice length must match n_gauss_points"
        );
        for (state, strain) in self.damage_states.iter_mut().zip(strains.iter()) {
            if state.deleted {
                continue;
            }
            let old_d = state.d;
            state.kappa = IsotropicDamage::update_kappa(state.kappa, strain);
            state.d = self.material_damage.damage_variable(state.kappa);
            state.damage_rate = (state.d - old_d).max(0.0);
        }
    }
    /// Update all Gauss points with rate limiting.
    pub fn update_all_limited(&mut self, strains: &[[f64; 6]], limiter: &DamageRateLimiter) {
        assert_eq!(strains.len(), self.n_gauss_points);
        for (state, strain) in self.damage_states.iter_mut().zip(strains.iter()) {
            if state.deleted {
                continue;
            }
            state.kappa = IsotropicDamage::update_kappa(state.kappa, strain);
            let target_d = self.material_damage.damage_variable(state.kappa);
            limiter.apply(state, target_d);
        }
    }
    /// Maximum damage value across all Gauss points.
    pub fn max_damage(&self) -> f64 {
        self.damage_states
            .iter()
            .map(|s| s.d)
            .fold(0.0_f64, f64::max)
    }
    /// Average damage value across all non-deleted Gauss points.
    pub fn avg_damage(&self) -> f64 {
        let active: Vec<f64> = self
            .damage_states
            .iter()
            .filter(|s| !s.deleted)
            .map(|s| s.d)
            .collect();
        if active.is_empty() {
            return 0.0;
        }
        active.iter().sum::<f64>() / active.len() as f64
    }
    /// Returns `true` when the maximum damage exceeds the supplied threshold.
    pub fn is_failed(&self, threshold: f64) -> bool {
        self.max_damage() >= threshold
    }
    /// Check element deletion: if all Gauss points exceed the threshold,
    /// mark the element as deleted.
    pub fn check_deletion(&mut self, threshold: f64) {
        if self.deleted {
            return;
        }
        let all_failed = self.damage_states.iter().all(|s| s.d >= threshold);
        if all_failed {
            self.deleted = true;
            for state in &mut self.damage_states {
                state.deleted = true;
            }
        }
    }
    /// Delete individual Gauss points that exceed the threshold.
    pub fn delete_failed_points(&mut self, threshold: f64) {
        for state in &mut self.damage_states {
            if state.d >= threshold {
                state.deleted = true;
            }
        }
        if self.damage_states.iter().all(|s| s.deleted) {
            self.deleted = true;
        }
    }
    /// Count the number of deleted Gauss points.
    pub fn deleted_point_count(&self) -> usize {
        self.damage_states.iter().filter(|s| s.deleted).count()
    }
    /// Count the number of active (non-deleted) Gauss points.
    pub fn active_point_count(&self) -> usize {
        self.n_gauss_points - self.deleted_point_count()
    }
    /// Get visualization data for this element.
    pub fn visualization(&self) -> DamageVisualization {
        DamageVisualization::from_states(&self.damage_states)
    }
    /// Get stiffness reduction factors at all Gauss points.
    pub fn stiffness_reduction_factors(&self) -> Vec<f64> {
        self.damage_states
            .iter()
            .map(|s| if s.deleted { 0.0 } else { 1.0 - s.d })
            .collect()
    }
}
/// Parameter bundle for [`CoupledDamagePlasticity::new`].
///
/// Groups the 11 material constants so the constructor stays within the
/// argument-count limit.
#[derive(Debug, Clone, Copy)]
pub struct CoupledDamagePlasticityParams {
    /// Young's modulus E \[Pa\]
    pub e: f64,
    /// Poisson's ratio ν
    pub nu: f64,
    /// Initial yield stress σ_y0 \[Pa\]
    pub sigma_y0: f64,
    /// Isotropic hardening modulus H \[Pa\]
    pub h: f64,
    /// Isotropic hardening saturation stress Q \[Pa\]
    pub q: f64,
    /// Isotropic hardening rate b
    pub b_hard: f64,
    /// Kinematic hardening stiffness C \[Pa\]
    pub c_kin: f64,
    /// Kinematic hardening recall γ
    pub gamma_kin: f64,
    /// Damage energy denominator S \[Pa\]
    pub s_dmg: f64,
    /// Damage exponent s
    pub s_dmg_exp: f64,
    /// Critical damage D_c
    pub d_c: f64,
}

/// Coupled isotropic damage-plasticity model (Lemaitre-Chaboche framework).
///
/// Combines isotropic hardening plasticity with scalar damage using
/// the effective-stress concept: σ̃ = σ / (1 − D).
///
/// # Reference
/// Lemaitre, J. & Chaboche, J.-L. (1990). *Mechanics of Solid Materials*.
/// Cambridge University Press.
#[derive(Debug, Clone)]
pub struct CoupledDamagePlasticity {
    /// Young's modulus E.
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
    /// Initial yield stress σ_y0.
    pub sigma_y0: f64,
    /// Isotropic hardening modulus H.
    pub h: f64,
    /// Isotropic hardening saturation stress Q.
    pub q: f64,
    /// Isotropic hardening rate b.
    pub b_hard: f64,
    /// Kinematic hardening stiffness C.
    pub c_kin: f64,
    /// Kinematic hardening recall γ.
    pub gamma_kin: f64,
    /// Damage energy denominator S.
    pub s_dmg: f64,
    /// Damage exponent s.
    pub s_dmg_exp: f64,
    /// Critical damage D_c.
    pub d_c: f64,
}
impl CoupledDamagePlasticity {
    /// Create a new coupled damage-plasticity model from a
    /// [`CoupledDamagePlasticityParams`] bundle.
    pub fn new(p: CoupledDamagePlasticityParams) -> Self {
        Self {
            e: p.e,
            nu: p.nu,
            sigma_y0: p.sigma_y0,
            h: p.h,
            q: p.q,
            b_hard: p.b_hard,
            c_kin: p.c_kin,
            gamma_kin: p.gamma_kin,
            s_dmg: p.s_dmg,
            s_dmg_exp: p.s_dmg_exp,
            d_c: p.d_c,
        }
    }
    /// Current yield stress including isotropic hardening.
    ///
    /// σ_y(R) = σ_y0 + H·p̄ + Q·(1 − e^{−b·p̄})
    pub fn yield_stress(&self, p_bar: f64) -> f64 {
        self.sigma_y0 + self.h * p_bar + self.q * (1.0 - (-self.b_hard * p_bar).exp())
    }
    /// Von-Mises yield function in effective stress space.
    ///
    /// f = σ̃_eq − σ_y(p̄) ≤ 0 for elastic state.
    ///
    /// `overstress` – overstress relative to von-Mises surface.
    pub fn yield_function(&self, sigma_eff: &[f64; 6], back_stress: &[f64; 6], p_bar: f64) -> f64 {
        let shifted: [f64; 6] = {
            let mut s = [0.0; 6];
            for i in 0..6 {
                s[i] = sigma_eff[i] - back_stress[i];
            }
            s
        };
        von_mises_local(&shifted) - self.yield_stress(p_bar)
    }
    /// Return-mapping (radial return) update for one increment.
    ///
    /// Given trial stress `sigma_trial` and current state, returns
    /// (updated_state, converged).
    pub fn return_mapping(
        &self,
        state: &DamagePlasticityState,
        sigma_trial: &[f64; 6],
    ) -> (DamagePlasticityState, bool) {
        let d = state.d;
        let factor = 1.0 / (1.0 - d).max(1e-12);
        let sigma_eff_trial: [f64; 6] = {
            let mut s = [0.0; 6];
            for i in 0..6 {
                s[i] = sigma_trial[i] * factor;
            }
            s
        };
        let f_trial = self.yield_function(&sigma_eff_trial, &state.back_stress, state.p_bar);
        if f_trial <= 0.0 {
            let mut new_state = state.clone();
            new_state.effective_stress = sigma_eff_trial;
            return (new_state, true);
        }
        let g = self.shear_modulus();
        let n_eff = von_mises_local(&sigma_eff_trial).max(1e-30);
        let n_vec: [f64; 6] = {
            let mut n = [0.0; 6];
            let s: [f64; 6] = deviatoric_stress(&sigma_eff_trial);
            for i in 0..6 {
                n[i] = 1.5 * s[i] / n_eff;
            }
            n
        };
        let mut dp = 0.0f64;
        for _ in 0..50 {
            let p_new = state.p_bar + dp;
            let sy_new = self.yield_stress(p_new);
            let res = n_eff - 3.0 * g * dp - sy_new;
            if res.abs() < 1e-10 * self.sigma_y0 {
                break;
            }
            let dsy = self.h + self.q * self.b_hard * (-self.b_hard * p_new).exp();
            let slope = -3.0 * g - dsy;
            if slope.abs() < 1e-30 {
                break;
            }
            dp -= res / slope;
            dp = dp.max(0.0);
        }
        let mut sigma_eff_new = [0.0f64; 6];
        for i in 0..6 {
            sigma_eff_new[i] = sigma_eff_trial[i] - 2.0 * g * dp * n_vec[i];
        }
        let mut back_new = state.back_stress;
        for i in 0..6 {
            back_new[i] += self.c_kin * dp * n_vec[i] - self.gamma_kin * back_new[i] * dp;
        }
        let sigma_eq = von_mises_local(&sigma_eff_new).max(1e-30);
        let sigma_h = (sigma_eff_new[0] + sigma_eff_new[1] + sigma_eff_new[2]) / 3.0;
        let triax = sigma_h / sigma_eq;
        let g_mod = self.shear_modulus();
        let nu = self.nu;
        let rv = 2.0 / 3.0 * (1.0 + nu) + 3.0 * (1.0 - 2.0 * nu) * triax * triax;
        let _ = g_mod;
        let denom = 2.0 * self.e * (1.0 - d).powi(2).max(1e-30);
        let y = sigma_eq * sigma_eq / denom * rv;
        let dd = (y / self.s_dmg).powf(self.s_dmg_exp) * dp;
        let p_new = state.p_bar + dp;
        let new_d = (d + dd).min(self.d_c).clamp(0.0, 1.0);
        let new_state = DamagePlasticityState {
            d: new_d,
            p_bar: p_new,
            r: state.r + self.h * dp,
            back_stress: back_new,
            effective_stress: sigma_eff_new,
        };
        (new_state, true)
    }
    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.e / (2.0 * (1.0 + self.nu))
    }
}
/// Non-local damage regularization.
///
/// Averages the equivalent strain over a neighborhood to prevent
/// mesh-dependent localization.
#[derive(Debug, Clone)]
pub struct NonLocalDamage {
    /// Characteristic (internal) length l_c (m).
    pub characteristic_length: f64,
}
impl NonLocalDamage {
    /// Create a new non-local damage model.
    pub fn new(characteristic_length: f64) -> Self {
        Self {
            characteristic_length,
        }
    }
    /// Gaussian weight function: w(r) = exp(-r^2 / (2 * l_c^2)).
    pub fn weight(&self, distance: f64) -> f64 {
        let lc = self.characteristic_length;
        (-distance * distance / (2.0 * lc * lc)).exp()
    }
    /// Compute the non-local equivalent strain at a point, given
    /// the local equivalent strains and distances from neighboring points.
    ///
    /// kappa_nl = sum(w_i * kappa_i) / sum(w_i)
    pub fn nonlocal_strain(&self, local_strains: &[f64], distances: &[f64]) -> f64 {
        assert_eq!(local_strains.len(), distances.len());
        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        for (&kappa, &dist) in local_strains.iter().zip(distances.iter()) {
            let w = self.weight(dist);
            weighted_sum += w * kappa;
            weight_sum += w;
        }
        if weight_sum > 1e-30 {
            weighted_sum / weight_sum
        } else {
            0.0
        }
    }
    /// Coupled local/non-local strain:
    /// kappa_coupled = alpha * kappa_local + (1 - alpha) * kappa_nonlocal
    pub fn coupled_strain(&self, kappa_local: f64, kappa_nonlocal: f64, alpha: f64) -> f64 {
        let alpha = alpha.clamp(0.0, 1.0);
        alpha * kappa_local + (1.0 - alpha) * kappa_nonlocal
    }
}
/// Generalized damage evolution law for isotropic damage mechanics.
///
/// Supports exponential, linear, and power-law softening.
#[derive(Debug, Clone)]
pub struct DamageEvolutionLaw {
    /// Damage threshold strain (strain at which damage initiates).
    pub kappa0: f64,
    /// Residual stiffness factor (0 = complete failure, 1 = no softening).
    pub alpha: f64,
    /// Softening exponent or shape parameter (β).
    pub beta: f64,
    /// Softening type.
    pub law_type: DamageLawType,
}
impl DamageEvolutionLaw {
    /// Create an exponential damage law.
    pub fn exponential(kappa0: f64, alpha: f64, beta: f64) -> Self {
        Self {
            kappa0,
            alpha,
            beta,
            law_type: DamageLawType::Exponential,
        }
    }
    /// Create a linear damage law (alpha = ultimate strain).
    pub fn linear(kappa0: f64, kappa_u: f64) -> Self {
        Self {
            kappa0,
            alpha: 1.0,
            beta: kappa_u,
            law_type: DamageLawType::Linear,
        }
    }
    /// Create a power-law damage law.
    pub fn power_law(kappa0: f64, beta: f64) -> Self {
        Self {
            kappa0,
            alpha: 1.0,
            beta,
            law_type: DamageLawType::PowerLaw,
        }
    }
    /// Compute the damage variable D for a given equivalent strain kappa.
    ///
    /// Returns D in \[0, 1\]. If kappa < kappa0, returns 0 (no damage).
    pub fn compute_damage(&self, kappa: f64) -> f64 {
        if kappa <= self.kappa0 {
            return 0.0;
        }
        let d = match self.law_type {
            DamageLawType::Exponential => {
                1.0 - (self.kappa0 / kappa)
                    * (1.0 - self.alpha + self.alpha * (-self.beta * (kappa - self.kappa0)).exp())
            }
            DamageLawType::Linear => {
                let kappa_u = self.beta;
                if kappa >= kappa_u {
                    1.0
                } else {
                    self.alpha * (kappa - self.kappa0) / (kappa_u - self.kappa0)
                }
            }
            DamageLawType::PowerLaw => 1.0 - (self.kappa0 / kappa).powf(self.beta),
        };
        d.clamp(0.0, 1.0)
    }
    /// Compute the damage rate dD/dkappa for a given kappa (for consistent tangent).
    pub fn compute_damage_rate(&self, kappa: f64) -> f64 {
        if kappa <= self.kappa0 {
            return 0.0;
        }
        match self.law_type {
            DamageLawType::Exponential => {
                let ratio = self.kappa0 / kappa;
                let exp_term = (-self.beta * (kappa - self.kappa0)).exp();
                let part1 = ratio / kappa * (1.0 - self.alpha + self.alpha * exp_term);
                let part2 = ratio * self.alpha * self.beta * exp_term;
                part1 + part2
            }
            DamageLawType::Linear => {
                let kappa_u = self.beta;
                if kappa >= kappa_u {
                    0.0
                } else {
                    self.alpha / (kappa_u - self.kappa0)
                }
            }
            DamageLawType::PowerLaw => {
                self.beta * self.kappa0.powf(self.beta) / kappa.powf(self.beta + 1.0)
            }
        }
    }
    /// Update damage state at a Gauss point using history variable.
    ///
    /// Returns updated `DamageState` with new damage and kappa values.
    pub fn update_state(&self, state: &DamageState, strain_equivalent: f64) -> DamageState {
        let kappa_new = state.kappa.max(strain_equivalent);
        let d_new = self.compute_damage(kappa_new).max(state.d);
        DamageState {
            d: d_new,
            kappa: kappa_new,
            plastic_strain: state.plastic_strain,
            accumulated_plastic_strain: state.accumulated_plastic_strain,
            damage_rate: d_new - state.d,
            deleted: d_new >= 0.9999,
        }
    }
}
/// State of a coupled damage-plasticity point.
#[derive(Debug, Clone)]
pub struct DamagePlasticityState {
    /// Damage variable D ∈ \[0, 1\].
    pub d: f64,
    /// Accumulated plastic strain p̄.
    pub p_bar: f64,
    /// Isotropic hardening variable R.
    pub r: f64,
    /// Back-stress (kinematic hardening) in Voigt notation \[αxx, αyy, αzz, αxy, αyz, αxz\].
    pub back_stress: [f64; 6],
    /// Effective (net-section) stress tensor.
    pub effective_stress: [f64; 6],
}
/// Chaboche fatigue damage model.
///
/// Describes high-cycle fatigue damage evolution as a function of the
/// number of cycles and the applied stress amplitude.
///
/// # Reference
/// Chaboche, J.-L. (1988). *Continuum Damage Mechanics: Part I–A Time and Cycle
/// Dependent Damage Model*. J. Appl. Mech. 55(1):59.
#[derive(Debug, Clone)]
pub struct FatigueDamageModel {
    /// Material constant M_0 controlling the endurance limit.
    pub m0: f64,
    /// Exponent beta controlling the damage accumulation rate.
    pub beta: f64,
    /// Endurance limit sigma_u (below which no fatigue occurs).
    pub sigma_u: f64,
    /// Exponent alpha for the stress-life relationship.
    pub alpha: f64,
    /// Mean stress sensitivity parameter b.
    pub b: f64,
}
impl FatigueDamageModel {
    /// Create a new fatigue damage model.
    pub fn new(m0: f64, beta: f64, sigma_u: f64, alpha: f64, b: f64) -> Self {
        Self {
            m0,
            beta,
            sigma_u,
            alpha,
            b,
        }
    }
    /// Fatigue damage increment per cycle dD/dN.
    ///
    /// dD/dN = \[1 - (1 - D)^(beta+1)\]^alpha * (sigma_a / M_0)^beta
    ///
    /// * `sigma_a` – stress amplitude
    /// * `sigma_mean` – mean stress
    /// * `d` – current damage variable D
    pub fn damage_rate_per_cycle(&self, sigma_a: f64, sigma_mean: f64, d: f64) -> f64 {
        if sigma_a <= self.sigma_u {
            return 0.0;
        }
        let m_eff = self.m0 * (1.0 - self.b * sigma_mean / self.sigma_u).max(1e-12);
        let base = (sigma_a / m_eff).powf(self.beta);
        let d_clamped = d.clamp(0.0, 1.0 - 1e-10);
        let factor = (1.0 - (1.0 - d_clamped).powf(self.beta + 1.0)).powf(self.alpha);
        factor * base
    }
    /// Remaining life (number of cycles) from current damage to failure.
    ///
    /// Integrates dD/dN from D to 1 using a simple Euler step estimate.
    /// Returns the estimated number of remaining cycles.
    pub fn remaining_life(
        &self,
        sigma_a: f64,
        sigma_mean: f64,
        d_current: f64,
        n_steps: usize,
    ) -> f64 {
        let mut d = d_current.clamp(0.0, 1.0 - 1e-10);
        let dd = (1.0 - d) / n_steps as f64;
        let mut total_cycles = 0.0f64;
        for _ in 0..n_steps {
            let rate = self.damage_rate_per_cycle(sigma_a, sigma_mean, d);
            if rate < 1e-30 {
                return f64::INFINITY;
            }
            total_cycles += dd / rate;
            d += dd;
            if d >= 1.0 {
                break;
            }
        }
        total_cycles
    }
    /// Number of cycles to failure (Nf) from D = 0.
    pub fn cycles_to_failure(&self, sigma_a: f64, sigma_mean: f64) -> f64 {
        self.remaining_life(sigma_a, sigma_mean, 0.0, 1000)
    }
}
/// Lemaitre-Chaboche ductile damage model.
///
/// # References
/// Lemaitre, J. (1992). *A Course on Damage Mechanics*. Springer.
#[derive(Debug, Clone)]
pub struct LemaitreDamage {
    /// Damage energy release rate denominator S.
    pub s: f64,
    /// Damage exponent r.
    pub r: f64,
    /// Critical damage threshold D_c (typically 0.2 - 0.8).
    pub d_c: f64,
    /// Initial yield stress sigma_Y.
    pub yield_stress: f64,
    /// Isotropic hardening modulus H.
    pub hardening: f64,
}
impl LemaitreDamage {
    /// Triaxiality-modified equivalent (von-Mises) strain.
    ///
    /// e_eq = sigma_eq / (E * (1 - D))   where sigma_eq is the von-Mises stress.
    pub fn equivalent_strain(stress: &[f64; 6], d: f64, e_mod: f64, nu: f64) -> f64 {
        let _ = nu;
        let sigma_eq = von_mises(stress);
        let denom = e_mod * (1.0 - d).max(1e-12);
        sigma_eq / denom
    }
    /// Incremental damage rate  dD/de_p = (Y / S)^r.
    ///
    /// `y`  -- thermodynamic force Y
    /// `dp` -- equivalent plastic strain increment.
    pub fn damage_rate(&self, y: f64, dp: f64) -> f64 {
        (y / self.s).powf(self.r) * dp
    }
    /// Thermodynamic force (energy release rate) conjugate to damage:
    ///
    /// Y = sigma_eq^2 / (2 E (1-D)^2) * Rv
    ///
    /// where Rv = 2/3*(1+nu) + 3*(1-2nu)*(sigma_H/sigma_eq)^2  is the triaxiality function.
    pub fn thermodynamic_force(stress: &[f64; 6], d: f64, e_mod: f64, nu: f64) -> f64 {
        let sigma_eq = von_mises(stress).max(1e-30);
        let sigma_h = hydrostatic(stress);
        let triax = sigma_h / sigma_eq;
        let rv = 2.0 / 3.0 * (1.0 + nu) + 3.0 * (1.0 - 2.0 * nu) * triax * triax;
        let denom = 2.0 * e_mod * (1.0 - d).powi(2).max(1e-30);
        sigma_eq * sigma_eq / denom * rv
    }
    /// Advance damage state at one Gauss point.
    ///
    /// Updates `state.d` using the incremental Lemaitre law.
    /// The state is clamped to \[0, d_c\].
    pub fn update_damage(
        &self,
        state: &mut DamageState,
        stress: &[f64; 6],
        dp: f64,
        e_mod: f64,
        nu: f64,
    ) {
        let y = Self::thermodynamic_force(stress, state.d, e_mod, nu);
        let dd = self.damage_rate(y, dp);
        state.damage_rate = dd;
        state.d = (state.d + dd).min(self.d_c).clamp(0.0, 1.0);
        state.accumulated_plastic_strain += dp;
    }
    /// Advance damage state with rate limiting.
    pub fn update_damage_limited(
        &self,
        state: &mut DamageState,
        stress: &[f64; 6],
        dp: f64,
        e_mod: f64,
        nu: f64,
        limiter: &DamageRateLimiter,
    ) {
        let y = Self::thermodynamic_force(stress, state.d, e_mod, nu);
        let dd = self.damage_rate(y, dp);
        let dd_limited = limiter.limit(dd);
        state.damage_rate = dd_limited;
        state.d = (state.d + dd_limited).min(self.d_c).clamp(0.0, 1.0);
        state.accumulated_plastic_strain += dp;
    }
    /// Effective (net-section) stress:  sigma_tilde = sigma / (1 - D).
    pub fn effective_stress(stress: &[f64; 6], d: f64) -> [f64; 6] {
        let factor = 1.0 / (1.0 - d).max(1e-12);
        let mut out = [0.0_f64; 6];
        for i in 0..6 {
            out[i] = stress[i] * factor;
        }
        out
    }
    /// Damaged elastic stiffness: C_d = (1 - D) * C_intact.
    pub fn damaged_stiffness(d: f64, c_intact: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
        let factor = 1.0 - d;
        let mut out = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                out[i][j] = factor * c_intact[i][j];
            }
        }
        out
    }
}
/// Implicit gradient nonlocal damage model (Peerlings et al., 1996).
///
/// The nonlocal equivalent strain ε̃ satisfies:
///
///   ε̃ − c · ∇²ε̃ = ε̄_local
///
/// where c = l_c² / 2 is the gradient parameter and l_c is the
/// characteristic length.
///
/// # Reference
/// Peerlings, R.H.J., de Borst, R., Brekelmans, W.A.M., & de Vree, J.H.P.
/// (1996). Gradient enhanced damage for quasi-brittle materials.
/// Int. J. Numer. Meth. Engng. 39:3391.
#[derive(Debug, Clone)]
pub struct NonlocalContinuumDamage {
    /// Characteristic length l_c (m).
    pub characteristic_length: f64,
    /// Damage threshold strain ε_0.
    pub epsilon_0: f64,
    /// Residual stress factor (0 = brittle, 1 = no softening).
    pub kappa_d: f64,
    /// Exponential softening slope α.
    pub alpha_soft: f64,
}
impl NonlocalContinuumDamage {
    /// Create a new nonlocal continuum damage model.
    pub fn new(characteristic_length: f64, epsilon_0: f64, kappa_d: f64, alpha_soft: f64) -> Self {
        Self {
            characteristic_length,
            epsilon_0,
            kappa_d,
            alpha_soft,
        }
    }
    /// Gradient parameter c = l_c² / 2.
    pub fn gradient_parameter(&self) -> f64 {
        self.characteristic_length * self.characteristic_length / 2.0
    }
    /// Exponential softening damage function g(κ):
    ///
    /// g(κ) = 1 − ε_0/κ · \[(1 − α) + α · e^{−β(κ − ε_0)}\]
    pub fn damage_loading_function(&self, kappa: f64, beta: f64) -> f64 {
        if kappa <= self.epsilon_0 {
            return 0.0;
        }
        let ratio = self.epsilon_0 / kappa;
        let d = 1.0
            - ratio
                * ((1.0 - self.alpha_soft)
                    + self.alpha_soft * (-(beta * (kappa - self.epsilon_0))).exp());
        d.clamp(0.0, 1.0)
    }
    /// Compute the nonlocal equivalent strain at a set of Gauss points using
    /// a 1D Green's function approximation.
    ///
    /// The Green's function for the implicit gradient equation in 1D is:
    ///
    ///   G(x, x') = 1/(2√c) · e^{−|x−x'|/√c}
    ///
    /// `positions` – 1D positions of Gauss points
    /// `local_strains` – local equivalent strains at those points
    ///
    /// Returns the nonlocal equivalent strains.
    pub fn compute_nonlocal_strains(&self, positions: &[f64], local_strains: &[f64]) -> Vec<f64> {
        assert_eq!(positions.len(), local_strains.len());
        let c = self.gradient_parameter();
        let lc = c.sqrt().max(1e-30);
        let n = positions.len();
        let mut result = vec![0.0f64; n];
        for i in 0..n {
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            for j in 0..n {
                let dist = (positions[i] - positions[j]).abs();
                let w = (-dist / lc).exp() / (2.0 * lc);
                num += w * local_strains[j];
                den += w;
            }
            result[i] = if den > 1e-30 {
                num / den
            } else {
                local_strains[i]
            };
        }
        result
    }
    /// Localization indicator: returns the ratio of the damage zone width to
    /// the characteristic length. A value > 1 indicates mesh-independent
    /// localization zone.
    pub fn localization_indicator(&self, damage_profile: &[f64], dx: f64) -> f64 {
        let threshold = 0.05;
        let n_damaged = damage_profile.iter().filter(|&&d| d > threshold).count();
        let zone_width = n_damaged as f64 * dx;
        zone_width / self.characteristic_length
    }
}
/// Temperature-dependent damage model.
///
/// Material properties degrade with increasing temperature, and
/// thermal gradients can drive additional damage.
#[derive(Debug, Clone)]
pub struct ThermalDamage {
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Melting/critical temperature T_melt (K).
    pub t_melt: f64,
    /// Reference Young's modulus E_0.
    pub e0: f64,
    /// Thermal degradation exponent m.
    pub m_exp: f64,
}
impl ThermalDamage {
    /// Create a new thermal damage model.
    pub fn new(t_ref: f64, t_melt: f64, e0: f64, m_exp: f64) -> Self {
        Self {
            t_ref,
            t_melt,
            e0,
            m_exp,
        }
    }
    /// Temperature-dependent stiffness reduction factor.
    ///
    /// E(T) / E_0 = (1 − (T − T_ref)/(T_melt − T_ref))^m
    pub fn stiffness_factor(&self, temperature: f64) -> f64 {
        let theta = ((temperature - self.t_ref) / (self.t_melt - self.t_ref)).clamp(0.0, 1.0);
        (1.0 - theta).powf(self.m_exp).max(0.0)
    }
    /// Effective Young's modulus at temperature T.
    pub fn effective_modulus(&self, temperature: f64) -> f64 {
        self.e0 * self.stiffness_factor(temperature)
    }
    /// Thermal damage variable D_T = 1 − E(T)/E_0.
    pub fn thermal_damage_variable(&self, temperature: f64) -> f64 {
        1.0 - self.stiffness_factor(temperature)
    }
    /// Combined mechanical-thermal damage:
    ///
    /// D_total = 1 − (1 − D_mech) · (1 − D_T)
    pub fn combined_damage(&self, d_mech: f64, temperature: f64) -> f64 {
        let d_t = self.thermal_damage_variable(temperature);
        (1.0 - (1.0 - d_mech) * (1.0 - d_t)).clamp(0.0, 1.0)
    }
}
/// Crack band model (Bažant–Oh) for quasi-brittle fracture regularization.
#[derive(Debug, Clone)]
pub struct CrackBandModel {
    /// Fracture energy G_f (J/m² or N/m).
    pub gf: f64,
    /// Uniaxial compressive (or tensile peak) strength f_c.
    pub fc: f64,
    /// Crack band width h (element characteristic length).
    pub h: f64,
}
impl CrackBandModel {
    /// Softening slope E_s = −f_c² h / (2 G_f).
    ///
    /// Negative value indicates softening.
    pub fn softening_slope(&self) -> f64 {
        -(self.fc * self.fc * self.h) / (2.0 * self.gf)
    }
}
/// Simple isotropic scalar damage model with exponential softening.
#[derive(Debug, Clone)]
pub struct IsotropicDamage {
    /// Strain threshold at damage initiation kappa_0.
    pub threshold: f64,
    /// Softening parameter A (controls rate of stiffness degradation).
    pub a: f64,
}
impl IsotropicDamage {
    /// Equivalent strain used for initiation criterion:
    ///
    /// e_tilde = sqrt(sum_i `e_i`+^2)
    ///
    /// Only the first three (normal) components contribute.
    pub fn equivalent_strain(strain: &[f64; 6]) -> f64 {
        let pos = |x: f64| x.max(0.0);
        let s = pos(strain[0]).powi(2) + pos(strain[1]).powi(2) + pos(strain[2]).powi(2);
        s.sqrt()
    }
    /// Compute the damage variable D(kappa):
    ///
    /// * D = 0                                          if kappa <= kappa_0
    /// * D = 1 - (kappa_0/kappa)*exp(-A*(kappa - kappa_0))  otherwise
    pub fn damage_variable(&self, kappa: f64) -> f64 {
        if kappa <= self.threshold {
            return 0.0;
        }
        let d = 1.0 - (self.threshold / kappa) * (-(self.a * (kappa - self.threshold))).exp();
        d.clamp(0.0, 1.0)
    }
    /// Advance the history variable kappa: kappa_new = max(kappa_old, e_tilde).
    pub fn update_kappa(current_kappa: f64, current_strain: &[f64; 6]) -> f64 {
        let eps_eq = Self::equivalent_strain(current_strain);
        current_kappa.max(eps_eq)
    }
    /// Damaged (secant) stiffness tensor:  C_d = (1 - D) * C_intact.
    pub fn damaged_stiffness(d: f64, c_intact: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
        let factor = 1.0 - d;
        let mut out = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                out[i][j] = factor * c_intact[i][j];
            }
        }
        out
    }
    /// Compute stiffness reduction factor (1 - D).
    pub fn stiffness_reduction_factor(&self, kappa: f64) -> f64 {
        1.0 - self.damage_variable(kappa)
    }
    /// Tangent stiffness factor (derivative of (1-D)*E with respect to strain).
    /// For the exponential softening model this gives the secant modulus ratio.
    pub fn tangent_factor(&self, kappa: f64) -> f64 {
        if kappa <= self.threshold {
            return 1.0;
        }
        self.stiffness_reduction_factor(kappa)
    }
}
/// Limits the rate at which damage can evolve per step, preventing
/// numerical instabilities from sudden jumps.
#[derive(Debug, Clone)]
pub struct DamageRateLimiter {
    /// Maximum allowed damage increment per step.
    pub max_rate: f64,
}
impl DamageRateLimiter {
    /// Create a new rate limiter.
    pub fn new(max_rate: f64) -> Self {
        Self { max_rate }
    }
    /// Clamp a damage increment to the allowed range.
    pub fn limit(&self, dd: f64) -> f64 {
        dd.clamp(0.0, self.max_rate)
    }
    /// Apply rate limiting to a damage state update.
    pub fn apply(&self, state: &mut DamageState, new_d: f64) {
        let dd = (new_d - state.d).max(0.0);
        let limited_dd = self.limit(dd);
        state.damage_rate = limited_dd;
        state.d = (state.d + limited_dd).clamp(0.0, 1.0);
    }
}
/// Type of damage softening law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DamageLawType {
    /// Exponential softening: D = 1 - (kappa0/kappa) * (1 - alpha + alpha * exp(-beta*(kappa - kappa0)))
    Exponential,
    /// Linear softening: D = alpha * (kappa - kappa0) / (kappa_u - kappa0)
    Linear,
    /// Power-law softening: D = 1 - (kappa0/kappa)^beta
    PowerLaw,
}
/// Mazars concrete damage model with separate tension and compression branches.
///
/// # References
/// Mazars, J. (1984). *Application de la mecanique de l'endommagement au
/// comportement non lineaire et a la rupture du beton de structure*.
/// These de Doctorat d'Etat, Universite Paris VI.
#[derive(Debug, Clone)]
pub struct MazarsDamage {
    /// Compression parameter A_c (approx 1.4).
    pub ac: f64,
    /// Compression exponential slope B_c (approx 1900).
    pub bc: f64,
    /// Tension parameter A_t (approx 0.8).
    pub at: f64,
    /// Tension exponential slope B_t (approx 20 000).
    pub bt: f64,
    /// Initial damage threshold strain kappa_0 (approx 1e-4).
    pub kappa0: f64,
}
impl MazarsDamage {
    /// Mazars equivalent strain:
    ///
    /// e_tilde = sqrt(`e_xx`+^2 + `e_yy`+^2 + `e_zz`+^2)
    ///
    /// where <.>+ = max(., 0).  (Simplified diagonal form.)
    pub fn equivalent_strain_mazars(strain: &[f64; 6]) -> f64 {
        let pos = |x: f64| x.max(0.0);
        let s = pos(strain[0]).powi(2) + pos(strain[1]).powi(2) + pos(strain[2]).powi(2);
        s.sqrt()
    }
    /// Tension damage branch:
    ///
    /// D_t(kappa) = 1 - kappa_0*(1 - A_t)/kappa - A_t*exp(-B_t*(kappa - kappa_0))
    pub fn damage_tension(&self, kappa: f64) -> f64 {
        if kappa <= self.kappa0 {
            return 0.0;
        }
        let d = 1.0
            - self.kappa0 * (1.0 - self.at) / kappa
            - self.at * (-(self.bt * (kappa - self.kappa0))).exp();
        d.clamp(0.0, 1.0)
    }
    /// Compression damage branch:
    ///
    /// D_c(kappa) = 1 - kappa_0*(1 - A_c)/kappa - A_c*exp(-B_c*(kappa - kappa_0))
    pub fn damage_compression(&self, kappa: f64) -> f64 {
        if kappa <= self.kappa0 {
            return 0.0;
        }
        let d = 1.0
            - self.kappa0 * (1.0 - self.ac) / kappa
            - self.ac * (-(self.bc * (kappa - self.kappa0))).exp();
        d.clamp(0.0, 1.0)
    }
    /// Update the damage state for a given strain tensor.
    ///
    /// `alpha_t` in \[0, 1\] is the tension fraction used to blend D_t and D_c:
    ///
    /// D = alpha_t * D_t(kappa) + (1 - alpha_t) * D_c(kappa)
    pub fn update(&self, state: &mut DamageState, strain: &[f64; 6], alpha_t: f64) {
        let eps_eq = Self::equivalent_strain_mazars(strain);
        let kappa = eps_eq.max(state.kappa);
        state.kappa = kappa;
        if kappa <= self.kappa0 {
            state.d = 0.0;
            return;
        }
        let dt = self.damage_tension(kappa);
        let dc = self.damage_compression(kappa);
        let alpha_t = alpha_t.clamp(0.0, 1.0);
        let new_d = (alpha_t * dt + (1.0 - alpha_t) * dc).clamp(0.0, 1.0);
        state.damage_rate = (new_d - state.d).max(0.0);
        state.d = new_d;
    }
    /// Update with rate limiting.
    pub fn update_limited(
        &self,
        state: &mut DamageState,
        strain: &[f64; 6],
        alpha_t: f64,
        limiter: &DamageRateLimiter,
    ) {
        let eps_eq = Self::equivalent_strain_mazars(strain);
        let kappa = eps_eq.max(state.kappa);
        state.kappa = kappa;
        if kappa <= self.kappa0 {
            state.d = 0.0;
            return;
        }
        let dt_val = self.damage_tension(kappa);
        let dc_val = self.damage_compression(kappa);
        let alpha_t = alpha_t.clamp(0.0, 1.0);
        let target_d = (alpha_t * dt_val + (1.0 - alpha_t) * dc_val).clamp(0.0, 1.0);
        limiter.apply(state, target_d);
    }
}
/// State variables stored at every Gauss point.
#[derive(Debug, Clone)]
pub struct DamageState {
    /// Damage variable D in \[0, 1\].  0 = undamaged, 1 = fully damaged.
    pub d: f64,
    /// Maximum equivalent strain ever reached (history / memory variable kappa).
    pub kappa: f64,
    /// Plastic strain tensor in Voigt notation \[e_xx, e_yy, e_zz, g_xy, g_yz, g_xz\].
    pub plastic_strain: [f64; 6],
    /// Accumulated (equivalent) plastic strain p_bar.
    pub accumulated_plastic_strain: f64,
    /// Damage rate from the last update (dD per step).
    pub damage_rate: f64,
    /// Whether this point has been flagged for deletion.
    pub deleted: bool,
}
/// Manages element deletion across a mesh.
#[derive(Debug, Clone)]
pub struct ElementDeletionManager {
    /// Damage threshold for element deletion.
    pub deletion_threshold: f64,
    /// Number of elements deleted so far.
    pub deleted_count: usize,
    /// IDs of deleted elements.
    pub deleted_ids: Vec<usize>,
}
impl ElementDeletionManager {
    /// Create a new deletion manager with the given threshold.
    pub fn new(threshold: f64) -> Self {
        Self {
            deletion_threshold: threshold,
            deleted_count: 0,
            deleted_ids: Vec::new(),
        }
    }
    /// Process a set of elements for deletion.
    pub fn process(&mut self, elements: &mut [DamageMechanicsElement]) {
        for elem in elements.iter_mut() {
            if elem.deleted {
                continue;
            }
            elem.check_deletion(self.deletion_threshold);
            if elem.deleted {
                self.deleted_count += 1;
                self.deleted_ids.push(elem.element_id);
            }
        }
    }
    /// Return the fraction of elements deleted.
    pub fn deletion_ratio(&self, total_elements: usize) -> f64 {
        if total_elements == 0 {
            return 0.0;
        }
        self.deleted_count as f64 / total_elements as f64
    }
}
/// Scalar isotropic damage model storing the current damage variable D,
/// current equivalent strain kappa, and initiation threshold kappa_0.
#[derive(Debug, Clone)]
pub struct ScalarIsotropicDamage {
    /// Current damage variable D ∈ \[0, 1\].
    pub d: f64,
    /// Current (maximum historical) equivalent strain kappa.
    pub kappa: f64,
    /// Damage initiation threshold kappa_0.
    pub kappa_0: f64,
}
impl ScalarIsotropicDamage {
    /// Equivalent strain: Euclidean norm of the strain vector.
    pub fn equivalent_strain(eps: &[f64]) -> f64 {
        eps.iter().map(|&e| e * e).sum::<f64>().sqrt()
    }
    /// Update internal variables given new equivalent strain `eps_eq`.
    ///
    /// Returns the updated damage D.  Uses linear damage loading:
    ///   D = 1 - kappa_0 / kappa  if kappa > kappa_0, else 0.
    pub fn update(&mut self, eps_eq: f64) -> f64 {
        if eps_eq > self.kappa {
            self.kappa = eps_eq;
        }
        if self.kappa > self.kappa_0 {
            let new_d = (1.0 - self.kappa_0 / self.kappa).clamp(0.0, 1.0);
            if new_d > self.d {
                self.d = new_d;
            }
        }
        self.d
    }
}
/// Full Lemaitre-Chaboche continuum damage mechanics model.
///
/// Implements the thermodynamically consistent framework with:
/// - Ductile damage driven by plastic dissipation
/// - Correction factor for triaxiality effects
/// - Damage threshold (no damage below p_D)
/// - Saturated damage at D_c
///
/// # Reference
/// Lemaitre, J. & Desmorat, R. (2005). *Engineering Damage Mechanics*.
/// Springer-Verlag.
#[derive(Debug, Clone)]
pub struct LemaitreChabocheDamage {
    /// Damage strength S.
    pub s_strength: f64,
    /// Damage exponent s.
    pub s: f64,
    /// Plastic strain threshold p_D before damage initiates.
    pub p_d: f64,
    /// Critical damage D_c (rupture criterion).
    pub d_c: f64,
    /// Young's modulus E.
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
}
impl LemaitreChabocheDamage {
    /// Triaxiality function R_v:
    ///
    /// R_v = 2/3(1+ν) + 3(1−2ν)(σ_H/σ_eq)²
    pub fn triaxiality_factor(&self, sigma_eq: f64, sigma_h: f64) -> f64 {
        if sigma_eq.abs() < 1e-30 {
            return 1.0;
        }
        let triax = sigma_h / sigma_eq;
        2.0 / 3.0 * (1.0 + self.nu) + 3.0 * (1.0 - 2.0 * self.nu) * triax * triax
    }
    /// Thermodynamic force Y conjugate to damage:
    ///
    /// Y = σ_eq² · R_v / (2E(1−D)²)
    pub fn damage_energy_release_rate(&self, sigma: &[f64; 6], d: f64) -> f64 {
        let sigma_eq = von_mises_local(sigma).max(1e-30);
        let sigma_h = (sigma[0] + sigma[1] + sigma[2]) / 3.0;
        let rv = self.triaxiality_factor(sigma_eq, sigma_h);
        let denom = 2.0 * self.e * (1.0 - d).powi(2).max(1e-30);
        sigma_eq * sigma_eq * rv / denom
    }
    /// Damage increment dD for a plastic strain increment dp.
    ///
    /// dD/dp = (Y/S)^s if p > p_D, else 0.
    pub fn damage_increment(&self, sigma: &[f64; 6], d: f64, dp: f64, p_bar: f64) -> f64 {
        if p_bar < self.p_d || dp <= 0.0 {
            return 0.0;
        }
        let y = self.damage_energy_release_rate(sigma, d);
        (y / self.s_strength).powf(self.s) * dp
    }
    /// Update a damage state by one increment.
    ///
    /// Returns the updated damage D (clamped to \[0, D_c\]).
    pub fn update(&self, state: &mut DamageState, sigma: &[f64; 6], dp: f64) {
        let dd = self.damage_increment(sigma, state.d, dp, state.accumulated_plastic_strain);
        state.damage_rate = dd;
        state.d = (state.d + dd).min(self.d_c).clamp(0.0, 1.0);
        state.accumulated_plastic_strain += dp;
    }
    /// Effective stress tensor: σ̃ = σ / (1−D).
    pub fn effective_stress(&self, sigma: &[f64; 6], d: f64) -> [f64; 6] {
        let factor = 1.0 / (1.0 - d).max(1e-12);
        let mut out = [0.0f64; 6];
        for i in 0..6 {
            out[i] = sigma[i] * factor;
        }
        out
    }
    /// Check whether the material has reached the rupture criterion D ≥ D_c.
    pub fn is_ruptured(&self, d: f64) -> bool {
        d >= self.d_c
    }
    /// Correction factor for plane-stress fracture mechanics.
    ///
    /// The Lemaitre crack closure parameter h corrects the damage coupling
    /// for compressive triaxiality.
    pub fn closure_corrected_damage(&self, d: f64, sigma: &[f64; 6], h: f64) -> f64 {
        let sigma_h = (sigma[0] + sigma[1] + sigma[2]) / 3.0;
        if sigma_h >= 0.0 { d } else { h * d }
    }
}
