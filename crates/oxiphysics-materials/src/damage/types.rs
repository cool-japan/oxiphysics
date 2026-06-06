//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Bailey-Norton primary creep model.
///
/// ε_cr = B * σ^m * t^p
#[derive(Debug, Clone)]
pub struct BaileyNortonCreep {
    /// Creep coefficient B
    pub b: f64,
    /// Stress exponent m
    pub m: f64,
    /// Time exponent p (< 1 for primary creep)
    pub p: f64,
}
impl BaileyNortonCreep {
    /// Create a new Bailey-Norton creep model.
    pub fn new(b: f64, m: f64, p: f64) -> Self {
        Self { b, m, p }
    }
    /// Creep strain: ε_cr = B * σ^m * t^p
    pub fn creep_strain(&self, sigma: f64, time: f64) -> f64 {
        self.b * sigma.powf(self.m) * time.powf(self.p)
    }
    /// Creep rate: ε̇_cr = B * σ^m * p * t^(p-1)
    pub fn creep_rate(&self, sigma: f64, time: f64) -> f64 {
        self.b * sigma.powf(self.m) * self.p * time.powf(self.p - 1.0)
    }
}
/// General CDM model with selectable damage evolution law.
///
/// Tracks a scalar damage variable D that degrades effective stress
/// and stiffness according to the principle of strain equivalence.
#[derive(Debug, Clone)]
pub struct CdmDamage {
    /// Current damage variable D ∈ \[0, 1\].
    pub d: f64,
    /// Damage initiation strain ε_0.
    pub eps0: f64,
    /// Failure strain ε_f (used in linear and power-law).
    pub eps_f: f64,
    /// Evolution rate parameter (exponential: a, power-law: n).
    pub param: f64,
    /// The damage evolution law to use.
    pub law: DamageEvolutionLaw,
    /// Maximum historical equivalent strain (for irreversibility).
    pub kappa: f64,
}
impl CdmDamage {
    /// Create a new CDM damage model.
    ///
    /// # Arguments
    /// * `eps0` - Damage initiation strain
    /// * `eps_f` - Failure strain (for Linear/PowerLaw)
    /// * `param` - Rate parameter (exponential: a, power-law exponent: n)
    /// * `law` - Which evolution law to use
    pub fn new(eps0: f64, eps_f: f64, param: f64, law: DamageEvolutionLaw) -> Self {
        Self {
            d: 0.0,
            eps0,
            eps_f,
            param,
            law,
            kappa: 0.0,
        }
    }
    /// Update damage based on current equivalent strain.
    ///
    /// Damage only increases (irreversibility via internal variable κ).
    pub fn update(&mut self, equivalent_strain: f64) {
        if equivalent_strain <= self.kappa {
            return;
        }
        self.kappa = equivalent_strain;
        if equivalent_strain <= self.eps0 {
            return;
        }
        let d_new = match self.law {
            DamageEvolutionLaw::Linear => {
                let range = self.eps_f - self.eps0;
                if range.abs() < f64::EPSILON {
                    1.0
                } else {
                    let eps = equivalent_strain;
                    (self.eps_f / eps) * (eps - self.eps0) / range
                }
            }
            DamageEvolutionLaw::Exponential => {
                let delta = equivalent_strain - self.eps0;
                1.0 - (-self.param * delta).exp()
            }
            DamageEvolutionLaw::PowerLaw => {
                let range = self.eps_f - self.eps0;
                if range.abs() < f64::EPSILON {
                    1.0
                } else {
                    let ratio = ((equivalent_strain - self.eps0) / range).min(1.0);
                    ratio.powf(self.param)
                }
            }
        };
        self.d = d_new.clamp(self.d, 1.0);
    }
    /// Effective stiffness: E_eff = (1 - D) * E.
    pub fn effective_stiffness(&self, young_modulus: f64) -> f64 {
        (1.0 - self.d) * young_modulus
    }
    /// Effective stress: σ_eff = σ / (1 - D).
    pub fn effective_stress(&self, nominal_stress: f64) -> f64 {
        let w = 1.0 - self.d;
        if w < 1e-15 {
            f64::INFINITY
        } else {
            nominal_stress / w
        }
    }
    /// Compute the strain localization indicator based on acoustic tensor analysis.
    ///
    /// In continuum damage mechanics, strain localization (shear band initiation)
    /// occurs when the acoustic tensor Q(n) becomes singular. The indicator is:
    ///
    /// L = D * E / ((1-D)² * (1 - 2ν) * (1 + ν) / (3*(1-ν)))
    ///
    /// or more precisely the determinant-like measure derived from the damaged
    /// tangent modulus. Here we use a practical scalar indicator:
    ///
    /// L = D / ((1-D)² + ε_reg)
    ///
    /// scaled by the ratio of shear to bulk modulus terms for 3D isotropy.
    ///
    /// # Arguments
    /// * `young_modulus`  - Undamaged Young's modulus E \[Pa\]
    /// * `poisson_ratio`  - Poisson's ratio ν (0 ≤ ν < 0.5)
    ///
    /// # Returns
    /// Localization indicator L (dimensionless, positive; → ∞ near full damage)
    pub fn compute_localization_indicator(&self, young_modulus: f64, poisson_ratio: f64) -> f64 {
        let nu = poisson_ratio;
        let w = 1.0 - self.d;
        let reg = 1e-15_f64;
        let damage_term = self.d / (w * w + reg);
        let geometric_factor = if (1.0 - nu).abs() > 1e-15 {
            (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu))
        } else {
            1.0
        };
        let _ = young_modulus;
        damage_term * geometric_factor.max(1e-15)
    }
}
/// Non-local integral averaging for damage regularization.
///
/// Computes a weighted average of a local field (e.g., equivalent strain)
/// over a neighborhood defined by a characteristic length l_c.
///
/// ε̄(x) = Σ w(|x - xᵢ|) * ε(xᵢ) / Σ w(|x - xᵢ|)
///
/// where w(r) = exp(-r² / (2*l_c²)) is the Gaussian weight.
#[derive(Debug, Clone)]
pub struct NonLocalRegularization {
    /// Characteristic (internal) length l_c (m).
    pub lc: f64,
}
impl NonLocalRegularization {
    /// Create a new non-local regularization with characteristic length l_c.
    pub fn new(lc: f64) -> Self {
        Self { lc }
    }
    /// Gaussian weight function: w(r) = exp(-r² / (2*l_c²)).
    pub fn weight(&self, distance: f64) -> f64 {
        (-distance * distance / (2.0 * self.lc * self.lc)).exp()
    }
    /// Bell-shaped (polynomial) weight function (compact support at r = l_c).
    ///
    /// w(r) = (1 - (r/l_c)²)² for r < l_c, 0 otherwise
    pub fn weight_bell(&self, distance: f64) -> f64 {
        if distance >= self.lc {
            0.0
        } else {
            let ratio = distance / self.lc;
            let t = 1.0 - ratio * ratio;
            t * t
        }
    }
    /// Compute non-local average of a field at a given point.
    ///
    /// # Arguments
    /// * `point` - Position at which to evaluate the non-local average
    /// * `positions` - Positions of all integration points
    /// * `values` - Local field values at each integration point
    pub fn nonlocal_average(
        &self,
        point: &[f64; 3],
        positions: &[[f64; 3]],
        values: &[f64],
    ) -> f64 {
        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        for (pos, &val) in positions.iter().zip(values.iter()) {
            let dx = point[0] - pos[0];
            let dy = point[1] - pos[1];
            let dz = point[2] - pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let w = self.weight(dist);
            weighted_sum += w * val;
            weight_sum += w;
        }
        if weight_sum.abs() < f64::EPSILON {
            0.0
        } else {
            weighted_sum / weight_sum
        }
    }
    /// Compute non-local averages for all points simultaneously.
    ///
    /// Returns a vector of non-local averages, one per point.
    pub fn nonlocal_average_all(&self, positions: &[[f64; 3]], values: &[f64]) -> Vec<f64> {
        positions
            .iter()
            .map(|pt| self.nonlocal_average(pt, positions, values))
            .collect()
    }
}
/// Isotropic damage model (Lemaitre 1985).
///
/// Damage variable D ∈ \[0, 1\]:
/// - D = 0: undamaged
/// - D = 1: fully damaged (cannot carry load)
///
/// Effective stress: σ_eff = σ / (1 - D)
#[derive(Debug, Clone)]
pub struct IsotropicDamage {
    /// Current damage variable D ∈ \[0, 1\]
    pub d: f64,
    /// Damage initiation threshold Y_0
    pub energy_release_threshold: f64,
    /// Exponential softening rate a
    pub damage_evolution_rate: f64,
}
impl IsotropicDamage {
    /// Create a new isotropic damage model.
    ///
    /// # Arguments
    /// * `y0` - Energy release threshold for damage initiation
    /// * `evolution_rate` - Exponential softening rate `a`
    pub fn new(y0: f64, evolution_rate: f64) -> Self {
        Self {
            d: 0.0,
            energy_release_threshold: y0,
            damage_evolution_rate: evolution_rate,
        }
    }
    /// Strain energy release rate: Y = σ² / (2 * E * (1-D)²)
    pub fn strain_energy_release_rate(&self, sigma: f64, young_modulus: f64) -> f64 {
        let integrity = 1.0 - self.d;
        sigma * sigma / (2.0 * young_modulus * integrity * integrity)
    }
    /// Update damage state based on current stress and stiffness.
    ///
    /// D = 1 - (Y_0/Y) * exp(-a*(Y-Y_0)/Y_0)  if Y > Y_0 and D_new > D
    pub fn update(&mut self, sigma: f64, young_modulus: f64) {
        let y = self.strain_energy_release_rate(sigma, young_modulus);
        let y0 = self.energy_release_threshold;
        if y > y0 {
            let a = self.damage_evolution_rate;
            let d_new = 1.0 - (y0 / y) * (-(a * (y - y0) / y0)).exp();
            let d_new = d_new.clamp(self.d, 1.0);
            self.d = d_new;
        }
    }
    /// Effective (degraded) stiffness: E_eff = E * (1 - D)
    pub fn effective_stiffness(&self, young_modulus: f64) -> f64 {
        young_modulus * (1.0 - self.d)
    }
    /// Effective stress (nominal stress divided by integrity): σ_eff = σ / (1 - D)
    pub fn effective_stress(&self, nominal_stress: f64) -> f64 {
        let integrity = 1.0 - self.d;
        if integrity < 1.0e-15 {
            f64::INFINITY
        } else {
            nominal_stress / integrity
        }
    }
}
/// Progressive failure analysis (PFA) for a composite laminate.
///
/// Each ply is assigned a scalar damage variable d ∈ \[0, 1\].
/// When d = 1 for a ply, it is considered fully failed and its stiffness
/// contribution is set to zero (ply-discount method).
#[derive(Debug, Clone)]
pub struct ProgressiveFailureAnalysis {
    /// Number of plies.
    pub n_plies: usize,
    /// Damage state of each ply (0 = intact, 1 = failed).
    pub ply_damage: Vec<f64>,
    /// Hashin criteria for this laminate.
    pub criteria: HashinFailureCriteria,
}
impl ProgressiveFailureAnalysis {
    /// Create a new PFA model with `n_plies` plies.
    pub fn new(n_plies: usize, x_t: f64, x_c: f64, y_t: f64, y_c: f64, s12: f64) -> Self {
        Self {
            n_plies,
            ply_damage: vec![0.0; n_plies],
            criteria: HashinFailureCriteria::new(x_t, x_c, y_t, y_c, s12),
        }
    }
    /// Degrade ply `i` by setting damage to `d` (irreversible).
    pub fn degrade_ply(&mut self, ply_idx: usize, d: f64) {
        if ply_idx < self.n_plies {
            self.ply_damage[ply_idx] = self.ply_damage[ply_idx].max(d).clamp(0.0, 1.0);
        }
    }
    /// Check if ply `i` is fully failed (d ≥ 1).
    pub fn is_ply_failed(&self, ply_idx: usize) -> bool {
        self.ply_damage.get(ply_idx).copied().unwrap_or(0.0) >= 1.0 - f64::EPSILON
    }
    /// Check if the entire laminate has failed (all plies failed).
    pub fn is_fully_failed(&self) -> bool {
        self.ply_damage.iter().all(|&d| d >= 1.0 - f64::EPSILON)
    }
    /// Mean damage across all plies.
    pub fn mean_damage(&self) -> f64 {
        if self.n_plies == 0 {
            return 0.0;
        }
        self.ply_damage.iter().sum::<f64>() / self.n_plies as f64
    }
    /// Apply ply-level stress state and update damage based on Hashin criteria.
    pub fn apply_stress(
        &mut self,
        ply_idx: usize,
        sigma1: f64,
        sigma2: f64,
        tau12: f64,
        d_increment: f64,
    ) {
        if ply_idx < self.n_plies && self.criteria.is_failed(sigma1, sigma2, tau12) {
            self.degrade_ply(ply_idx, self.ply_damage[ply_idx] + d_increment);
        }
    }
    /// Effective laminate stiffness factor (mean over surviving plies).
    pub fn laminate_stiffness_factor(&self) -> f64 {
        if self.n_plies == 0 {
            return 0.0;
        }
        self.ply_damage.iter().map(|&d| 1.0 - d).sum::<f64>() / self.n_plies as f64
    }
}
/// Damage-induced anisotropy model.
///
/// Maintains separate damage variables for the three principal directions
/// (d11, d22, d33) and the three shear planes (d12, d13, d23).
/// The effective compliance tensor accounts for directional degradation.
#[derive(Debug, Clone)]
pub struct DamageInducedAnisotropy {
    /// Undamaged Young's modulus \[Pa\].
    pub e0: f64,
    /// Undamaged Poisson's ratio.
    pub nu0: f64,
    /// Normal damage in 1-direction.
    pub d11: f64,
    /// Normal damage in 2-direction.
    pub d22: f64,
    /// Normal damage in 3-direction.
    pub d33: f64,
    /// Shear damage in 1-2 plane.
    pub d12: f64,
    /// Shear damage in 1-3 plane.
    pub d13: f64,
    /// Shear damage in 2-3 plane.
    pub d23: f64,
}
impl DamageInducedAnisotropy {
    /// Create a new anisotropic damage model (initially undamaged).
    pub fn new(e0: f64, nu0: f64) -> Self {
        Self {
            e0,
            nu0,
            d11: 0.0,
            d22: 0.0,
            d33: 0.0,
            d12: 0.0,
            d13: 0.0,
            d23: 0.0,
        }
    }
    /// Update shear damage variables (irreversible).
    pub fn update_shear_damage(&mut self, d12: f64, d13: f64, d23: f64) {
        self.d12 = self.d12.max(d12).clamp(0.0, 1.0 - f64::EPSILON);
        self.d13 = self.d13.max(d13).clamp(0.0, 1.0 - f64::EPSILON);
        self.d23 = self.d23.max(d23).clamp(0.0, 1.0 - f64::EPSILON);
    }
    /// Update normal damage variables (irreversible).
    pub fn update_normal_damage(&mut self, d11: f64, d22: f64, d33: f64) {
        self.d11 = self.d11.max(d11).clamp(0.0, 1.0 - f64::EPSILON);
        self.d22 = self.d22.max(d22).clamp(0.0, 1.0 - f64::EPSILON);
        self.d33 = self.d33.max(d33).clamp(0.0, 1.0 - f64::EPSILON);
    }
    /// Effective Young's modulus in direction 1: E1_eff = (1-D11)*E0.
    pub fn e1_eff(&self) -> f64 {
        (1.0 - self.d11) * self.e0
    }
    /// Effective Young's modulus in direction 2: E2_eff = (1-D22)*E0.
    pub fn e2_eff(&self) -> f64 {
        (1.0 - self.d22) * self.e0
    }
    /// Effective shear modulus G12_eff = (1-D12)*G0.
    pub fn effective_shear_modulus_12(&self) -> f64 {
        let g0 = self.e0 / (2.0 * (1.0 + self.nu0));
        (1.0 - self.d12) * g0
    }
    /// Scalar damage intensity (Frobenius norm of damage tensor components).
    pub fn damage_intensity(&self) -> f64 {
        (self.d11.powi(2)
            + self.d22.powi(2)
            + self.d33.powi(2)
            + 2.0 * self.d12.powi(2)
            + 2.0 * self.d13.powi(2)
            + 2.0 * self.d23.powi(2))
        .sqrt()
    }
}
/// Damage evolution law types for general CDM models.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DamageEvolutionLaw {
    /// Linear softening: D grows linearly with equivalent strain.
    Linear,
    /// Exponential softening: D = 1 - exp(-a*(ε-ε0)).
    Exponential,
    /// Power-law softening: D = ((ε-ε0)/(εf-ε0))^n.
    PowerLaw,
}
/// Mazars concrete damage model (biaxial damage).
///
/// Uses positive principal strains to distinguish tension/compression.
#[derive(Debug, Clone)]
pub struct MazarsDamage {
    /// Current damage variable D ∈ \[0, 1\]
    pub d: f64,
    /// Damage threshold (strain) k
    pub k: f64,
    /// Tension parameter A_t
    pub at: f64,
    /// Tension parameter B_t
    pub bt: f64,
    /// Compression parameter A_c
    pub ac: f64,
    /// Compression parameter B_c
    pub bc: f64,
}
impl MazarsDamage {
    /// Create a new Mazars damage model.
    ///
    /// # Arguments
    /// * `k0` - Initial damage threshold (strain)
    /// * `at`, `bt` - Tension damage parameters
    /// * `ac`, `bc` - Compression damage parameters
    pub fn new(k0: f64, at: f64, bt: f64, ac: f64, bc: f64) -> Self {
        Self {
            d: 0.0,
            k: k0,
            at,
            bt,
            ac,
            bc,
        }
    }
    /// Equivalent strain from positive principal strains only.
    ///
    /// ε̃ = sqrt(Σ <ε_i>+²)  where `x`+ = max(x, 0)
    pub fn equivalent_strain(eps_principal: &[f64; 3]) -> f64 {
        let sum_sq: f64 = eps_principal
            .iter()
            .map(|&e| {
                let ep = e.max(0.0);
                ep * ep
            })
            .sum();
        sum_sq.sqrt()
    }
    /// Update damage based on principal strains.
    ///
    /// D = at*(1 - k/ε̃)*exp(-bt*(ε̃-k)) + ac*(1 - k/ε̃)*exp(-bc*(ε̃-k))
    pub fn update(&mut self, eps_principal: &[f64; 3]) {
        let eps_tilde = Self::equivalent_strain(eps_principal);
        if eps_tilde > self.k {
            let k = self.k;
            let factor = 1.0 - k / eps_tilde;
            let delta = eps_tilde - k;
            let d_new = self.at * factor * (-self.bt * delta).exp()
                + self.ac * factor * (-self.bc * delta).exp();
            let d_new = d_new.clamp(self.d, 1.0);
            self.d = d_new;
        }
    }
}
/// Tsai-Wu polynomial failure criterion for composite laminates.
///
/// Combined criterion:
/// F₁σ₁ + F₂σ₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂ = 1  (failure)
///
/// Reference: Tsai & Wu (1971).
#[derive(Debug, Clone)]
pub struct TsaiWuCriterion {
    /// Longitudinal tensile strength X_T \[Pa\].
    pub x_t: f64,
    /// Longitudinal compressive strength X_C \[Pa\].
    pub x_c: f64,
    /// Transverse tensile strength Y_T \[Pa\].
    pub y_t: f64,
    /// Transverse compressive strength Y_C \[Pa\].
    pub y_c: f64,
    /// Shear strength S_12 \[Pa\].
    pub s12: f64,
    /// Interaction coefficient F₁₂* (typically -0.5).
    pub f12_star: f64,
}
impl TsaiWuCriterion {
    /// Create a new Tsai-Wu criterion.
    pub fn new(x_t: f64, x_c: f64, y_t: f64, y_c: f64, s12: f64, f12_star: f64) -> Self {
        Self {
            x_t,
            x_c,
            y_t,
            y_c,
            s12,
            f12_star,
        }
    }
    /// Compute the Tsai-Wu failure index.
    ///
    /// Returns value < 1 → safe, ≥ 1 → failure.
    pub fn failure_index(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let f1 = 1.0 / self.x_t - 1.0 / self.x_c;
        let f2 = 1.0 / self.y_t - 1.0 / self.y_c;
        let f11 = 1.0 / (self.x_t * self.x_c);
        let f22 = 1.0 / (self.y_t * self.y_c);
        let f66 = 1.0 / (self.s12 * self.s12);
        let f12 = self.f12_star / (self.x_t * self.x_c * self.y_t * self.y_c).sqrt();
        f1 * sigma1
            + f2 * sigma2
            + f11 * sigma1 * sigma1
            + f22 * sigma2 * sigma2
            + f66 * tau12 * tau12
            + 2.0 * f12 * sigma1 * sigma2
    }
    /// Strength ratio R: actual load can be multiplied by R to reach failure.
    ///
    /// Solves quadratic: A·R² + B·R - 1 = 0
    pub fn strength_ratio(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let f1 = 1.0 / self.x_t - 1.0 / self.x_c;
        let f2 = 1.0 / self.y_t - 1.0 / self.y_c;
        let f11 = 1.0 / (self.x_t * self.x_c);
        let f22 = 1.0 / (self.y_t * self.y_c);
        let f66 = 1.0 / (self.s12 * self.s12);
        let f12 = self.f12_star / (self.x_t * self.x_c * self.y_t * self.y_c).sqrt();
        let a_coef = f11 * sigma1 * sigma1
            + f22 * sigma2 * sigma2
            + f66 * tau12 * tau12
            + 2.0 * f12 * sigma1 * sigma2;
        let b_coef = f1 * sigma1 + f2 * sigma2;
        if a_coef.abs() < f64::EPSILON {
            if b_coef.abs() < f64::EPSILON {
                return f64::INFINITY;
            }
            return 1.0 / b_coef;
        }
        let discriminant = b_coef * b_coef + 4.0 * a_coef;
        if discriminant < 0.0 {
            return f64::INFINITY;
        }
        (-b_coef + discriminant.sqrt()) / (2.0 * a_coef)
    }
}
/// Repair effectiveness model for damaged composite structures.
///
/// After repair, the material recovers some fraction of its original strength
/// depending on the repair method quality and residual damage level.
///
/// Repair effectiveness η ∈ \[0, 1\]:
/// - η = 1: perfect repair (full strength recovery)
/// - η = 0: no repair (residual damage remains)
#[derive(Debug, Clone)]
pub struct RepairEffectiveness {
    /// Repair quality factor κ ∈ \[0, 1\] (1 = perfect repair quality).
    pub repair_factor: f64,
    /// Residual damage after repair d_res ∈ \[0, 1).
    pub residual_damage: f64,
}
impl RepairEffectiveness {
    /// Create a new repair effectiveness model.
    ///
    /// # Arguments
    /// * `repair_factor` - Quality of repair κ ∈ \[0, 1\]
    /// * `residual_damage` - Remaining damage after repair d_res
    pub fn new(repair_factor: f64, residual_damage: f64) -> Self {
        Self {
            repair_factor: repair_factor.clamp(0.0, 1.0),
            residual_damage: residual_damage.clamp(0.0, 1.0 - f64::EPSILON),
        }
    }
    /// Repair effectiveness η = κ · (1 - d_res).
    pub fn effectiveness(&self) -> f64 {
        self.repair_factor * (1.0 - self.residual_damage)
    }
    /// Recovered strength after repair.
    ///
    /// σ_rep = σ_0 · η
    pub fn strength_recovery(&self, original_strength: f64) -> f64 {
        original_strength * self.effectiveness()
    }
    /// Post-repair stiffness factor relative to undamaged.
    ///
    /// E_rep / E_0 = 1 - d_res · (1 - κ)
    pub fn stiffness_recovery_factor(&self) -> f64 {
        1.0 - self.residual_damage * (1.0 - self.repair_factor)
    }
    /// Check if the repair meets a minimum effectiveness threshold.
    pub fn is_acceptable(&self, min_effectiveness: f64) -> bool {
        self.effectiveness() >= min_effectiveness
    }
}
/// Lemaitre ductile damage model with plastic strain coupling.
///
/// Ḋ = (Y/S)^s * ε̇_p  for ε_p > ε_D
///
/// where Y is the energy release rate, S and s are material parameters,
/// and ε_D is the damage threshold plastic strain.
#[derive(Debug, Clone)]
pub struct LemaitreDuctileDamage {
    /// Current damage variable D ∈ \[0, 1\].
    pub d: f64,
    /// Accumulated plastic strain.
    pub eps_p: f64,
    /// Damage threshold plastic strain ε_D.
    pub eps_d: f64,
    /// Damage strength parameter S (Pa).
    pub s_param: f64,
    /// Damage exponent s.
    pub s_exp: f64,
    /// Critical damage D_c at which rupture occurs.
    pub d_critical: f64,
    /// Poisson's ratio ν (for triaxiality function).
    pub nu: f64,
}
impl LemaitreDuctileDamage {
    /// Create a new Lemaitre ductile damage model.
    pub fn new(eps_d: f64, s_param: f64, s_exp: f64, d_critical: f64, nu: f64) -> Self {
        Self {
            d: 0.0,
            eps_p: 0.0,
            eps_d,
            s_param,
            s_exp,
            d_critical,
            nu,
        }
    }
    /// Stress triaxiality function R_ν.
    ///
    /// R_ν = (2/3)(1 + ν) + 3(1 - 2ν)(σ_H/σ_eq)²
    ///
    /// where σ_H/σ_eq is the triaxiality ratio.
    pub fn triaxiality_function(&self, triaxiality: f64) -> f64 {
        (2.0 / 3.0) * (1.0 + self.nu) + 3.0 * (1.0 - 2.0 * self.nu) * triaxiality * triaxiality
    }
    /// Energy release rate: Y = σ_eq² * R_ν / (2*E*(1-D)²)
    pub fn energy_release_rate(&self, sigma_eq: f64, triaxiality: f64, young_modulus: f64) -> f64 {
        let rv = self.triaxiality_function(triaxiality);
        let w = 1.0 - self.d;
        sigma_eq * sigma_eq * rv / (2.0 * young_modulus * w * w)
    }
    /// Update damage given an increment of plastic strain.
    ///
    /// Ḋ = (Y/S)^s * Δε_p   if ε_p > ε_D and D < D_c
    pub fn update(
        &mut self,
        sigma_eq: f64,
        triaxiality: f64,
        young_modulus: f64,
        delta_eps_p: f64,
    ) {
        self.eps_p += delta_eps_p;
        if self.eps_p <= self.eps_d || self.d >= self.d_critical {
            return;
        }
        let y = self.energy_release_rate(sigma_eq, triaxiality, young_modulus);
        let dd = (y / self.s_param).powf(self.s_exp) * delta_eps_p;
        self.d = (self.d + dd).min(self.d_critical);
    }
    /// Check if rupture has occurred (D >= D_c).
    pub fn is_ruptured(&self) -> bool {
        self.d >= self.d_critical - 1e-15
    }
    /// Compute the stress triaxiality correction factor R_ν.
    ///
    /// The triaxiality correction amplifies the damage driving force under
    /// multiaxial stress states:
    ///
    /// R_ν = (2/3)(1 + ν) + 3(1 - 2ν) η²
    ///
    /// where η = σ_H / σ_eq is the stress triaxiality ratio.
    /// This factor accounts for the accelerated ductile fracture at high
    /// triaxiality (e.g., notch root, crack front).
    ///
    /// # Arguments
    /// * `triaxiality` - η = σ_hydrostatic / σ_von_Mises (dimensionless)
    ///
    /// # Returns
    /// Triaxiality correction factor R_ν ≥ (2/3)(1+ν)
    pub fn compute_triaxiality_correction(&self, triaxiality: f64) -> f64 {
        (2.0 / 3.0) * (1.0 + self.nu) + 3.0 * (1.0 - 2.0 * self.nu) * triaxiality * triaxiality
    }
}
/// Progressive stiffness degradation model for composite materials.
///
/// Tracks separate damage variables D11 (fiber direction) and D22
/// (transverse direction). Effective stiffness:
///
/// E11_eff = (1 - D11) * E11,   E22_eff = (1 - D22) * E22
#[derive(Debug, Clone)]
pub struct ProgressiveDamage {
    /// Undamaged Young's modulus (isotropic assumption) \[Pa\].
    pub e0: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Fiber-direction damage variable D11 ∈ \[0, 1).
    pub d11: f64,
    /// Transverse-direction damage variable D22 ∈ \[0, 1).
    pub d22: f64,
}
impl ProgressiveDamage {
    /// Create a new progressive damage model.
    pub fn new(e0: f64, nu: f64, d11: f64, d22: f64) -> Self {
        Self {
            e0,
            nu,
            d11: d11.clamp(0.0, 1.0 - f64::EPSILON),
            d22: d22.clamp(0.0, 1.0 - f64::EPSILON),
        }
    }
    /// Effective fiber-direction stiffness.
    pub fn e11_eff(&self) -> f64 {
        (1.0 - self.d11) * self.e0
    }
    /// Effective transverse stiffness.
    pub fn e22_eff(&self) -> f64 {
        (1.0 - self.d22) * self.e0
    }
    /// 3×3 effective stiffness tensor (plane stress, Voigt notation: 11, 22, 12).
    ///
    /// Returns `[[C11, C12, 0], [C12, C22, 0], [0, 0, C66]]`.
    pub fn effective_stiffness_tensor(&self) -> [[f64; 3]; 3] {
        let e1 = self.e11_eff();
        let e2 = self.e22_eff();
        let nu12 = self.nu;
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        let c11 = e1 / denom;
        let c22 = e2 / denom;
        let c12 = nu12 * e2 / denom;
        let g0 = self.e0 / (2.0 * (1.0 + self.nu));
        let d_shear = 1.0 - (1.0 - self.d11) * (1.0 - self.d22);
        let c66 = (1.0 - d_shear) * g0;
        [[c11, c12, 0.0], [c12, c22, 0.0], [0.0, 0.0, c66]]
    }
    /// Update damage variables (irreversible – only increases).
    pub fn update_damage(&mut self, d11_new: f64, d22_new: f64) {
        self.d11 = self.d11.max(d11_new).clamp(0.0, 1.0 - f64::EPSILON);
        self.d22 = self.d22.max(d22_new).clamp(0.0, 1.0 - f64::EPSILON);
    }
    /// Residual strength in fiber direction after damage.
    ///
    /// σ_res = σ_0 · (1 - D11)^alpha_f · (1 - D22)^alpha_m
    pub fn residual_strength(&self, sigma0: f64, _alpha_m: f64) -> f64 {
        sigma0 * (1.0 - self.d11)
    }
}
/// Gradient-enhanced damage regularization.
///
/// Implicit gradient formulation:
/// ε̄ - l_c² ∇²ε̄ = ε
///
/// This is a simplified 1D discretization for a single bar element.
#[derive(Debug, Clone)]
pub struct GradientEnhancedDamage {
    /// Characteristic length squared c = l_c².
    pub c: f64,
}
impl GradientEnhancedDamage {
    /// Create gradient-enhanced regularization from characteristic length l_c.
    pub fn new(lc: f64) -> Self {
        Self { c: lc * lc }
    }
    /// Solve the implicit gradient equation in 1D for a uniform bar.
    ///
    /// Discretized with finite differences:
    /// ε̄_i - c*(ε̄_{i-1} - 2*ε̄_i + ε̄_{i+1})/h² = ε_i
    ///
    /// Uses Thomas algorithm (tridiagonal solver).
    pub fn solve_1d(&self, local_strains: &[f64], element_size: f64) -> Vec<f64> {
        let n = local_strains.len();
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![local_strains[0]];
        }
        let h2 = element_size * element_size;
        let alpha = self.c / h2;
        let mut a_coeff = vec![0.0; n];
        let mut b_coeff = vec![0.0; n];
        let mut c_coeff = vec![0.0; n];
        let mut d_vec = vec![0.0; n];
        for i in 0..n {
            b_coeff[i] = 1.0 + 2.0 * alpha;
            d_vec[i] = local_strains[i];
            if i > 0 {
                a_coeff[i] = -alpha;
            }
            if i < n - 1 {
                c_coeff[i] = -alpha;
            }
        }
        b_coeff[0] = 1.0 + alpha;
        b_coeff[n - 1] = 1.0 + alpha;
        let mut cp = vec![0.0; n];
        let mut dp = vec![0.0; n];
        cp[0] = c_coeff[0] / b_coeff[0];
        dp[0] = d_vec[0] / b_coeff[0];
        for i in 1..n {
            let m = b_coeff[i] - a_coeff[i] * cp[i - 1];
            cp[i] = c_coeff[i] / m;
            dp[i] = (d_vec[i] - a_coeff[i] * dp[i - 1]) / m;
        }
        let mut result = vec![0.0; n];
        result[n - 1] = dp[n - 1];
        for i in (0..n - 1).rev() {
            result[i] = dp[i] - cp[i] * result[i + 1];
        }
        result
    }
}
/// Hashin failure criteria for unidirectional fiber composites.
///
/// Four independent failure modes:
/// 1. Fiber tension (σ_1 > 0)
/// 2. Fiber compression (σ_1 < 0)
/// 3. Matrix tension (σ_2 > 0)
/// 4. Matrix compression (σ_2 < 0)
///
/// Reference: Hashin (1980).
#[derive(Debug, Clone)]
pub struct HashinFailureCriteria {
    /// Longitudinal tensile strength X_T \[Pa\].
    pub x_t: f64,
    /// Longitudinal compressive strength X_C \[Pa\].
    pub x_c: f64,
    /// Transverse tensile strength Y_T \[Pa\].
    pub y_t: f64,
    /// Transverse compressive strength Y_C \[Pa\].
    pub y_c: f64,
    /// In-plane shear strength S_12 \[Pa\].
    pub s12: f64,
}
impl HashinFailureCriteria {
    /// Create a new Hashin failure criteria instance.
    pub fn new(x_t: f64, x_c: f64, y_t: f64, y_c: f64, s12: f64) -> Self {
        Self {
            x_t,
            x_c,
            y_t,
            y_c,
            s12,
        }
    }
    /// Fiber tension criterion: f_ft = (σ_1/X_T)^2 + (τ_12/S_12)^2.
    pub fn fiber_tension(&self, sigma1: f64, tau12: f64, _tau13: f64) -> f64 {
        if sigma1 <= 0.0 {
            return 0.0;
        }
        (sigma1 / self.x_t).powi(2) + (tau12 / self.s12).powi(2)
    }
    /// Fiber compression criterion: f_fc = (σ_1/X_C)^2.
    pub fn fiber_compression(&self, sigma1: f64) -> f64 {
        if sigma1 >= 0.0 {
            return 0.0;
        }
        (sigma1 / self.x_c).powi(2)
    }
    /// Matrix tension criterion: f_mt = (σ_2/Y_T)^2 + (τ_12/S_12)^2.
    pub fn matrix_tension(&self, sigma2: f64, tau12: f64) -> f64 {
        if sigma2 <= 0.0 {
            return 0.0;
        }
        (sigma2 / self.y_t).powi(2) + (tau12 / self.s12).powi(2)
    }
    /// Matrix compression criterion (Hashin 1980):
    ///
    /// f_mc = (σ_2 / (2·S_23))^2 + \[(Y_C / (2·S_23))^2 - 1\]·(σ_2/Y_C) + (τ_12/S_12)^2
    pub fn matrix_compression(&self, sigma2: f64, tau12: f64) -> f64 {
        if sigma2 >= 0.0 {
            return 0.0;
        }
        let s23 = self.y_c / 2.0;
        let t1 = (sigma2 / (2.0 * s23)).powi(2);
        let t2 = ((self.y_c / (2.0 * s23)).powi(2) - 1.0) * (sigma2 / self.y_c);
        let t3 = (tau12 / self.s12).powi(2);
        t1 + t2 + t3
    }
    /// Check if any failure mode is triggered (criterion ≥ 1).
    pub fn is_failed(&self, sigma1: f64, sigma2: f64, tau12: f64) -> bool {
        self.fiber_tension(sigma1, tau12, 0.0) >= 1.0
            || self.fiber_compression(sigma1) >= 1.0
            || self.matrix_tension(sigma2, tau12) >= 1.0
            || self.matrix_compression(sigma2, tau12) >= 1.0
    }
}
/// Norton power-law creep model.
///
/// ε̇_cr = A * |σ|^n * sign(σ) * exp(-Q/(R*T))
#[derive(Debug, Clone)]
pub struct NortonCreep {
    /// Pre-exponential factor A
    pub a: f64,
    /// Stress exponent n
    pub n: f64,
    /// Activation energy Q (J/mol)
    pub q: f64,
    /// Gas constant R (J/(mol·K))
    pub r: f64,
}
impl NortonCreep {
    /// Create a new Norton creep model with R = 8.314 J/(mol·K).
    ///
    /// # Arguments
    /// * `a` - Pre-exponential factor
    /// * `n` - Stress exponent
    /// * `q` - Activation energy (J/mol)
    pub fn new(a: f64, n: f64, q: f64) -> Self {
        Self { a, n, q, r: 8.314 }
    }
    /// Creep strain rate: ε̇_cr = A * |σ|^n * sign(σ) * exp(-Q/(R*T))
    pub fn creep_rate(&self, sigma: f64, temperature: f64) -> f64 {
        let sign = sigma.signum();
        let abs_sigma = sigma.abs();
        self.a * abs_sigma.powf(self.n) * sign * (-self.q / (self.r * temperature)).exp()
    }
    /// Creep strain increment in time dt.
    pub fn creep_increment(&self, sigma: f64, temperature: f64, dt: f64) -> f64 {
        self.creep_rate(sigma, temperature) * dt
    }
}
/// Puck failure criteria for unidirectional composites.
///
/// Separates failure into:
/// - Mode A: fiber failure (tension/compression)
/// - Mode B/C: inter-fiber failure (IFF)
///
/// Reference: Puck & Schürmann (1998).
#[derive(Debug, Clone)]
pub struct PuckFailureCriteria {
    /// Longitudinal tensile strength X_T \[Pa\].
    pub x_t: f64,
    /// Longitudinal compressive strength X_C \[Pa\].
    pub x_c: f64,
    /// Transverse tensile strength Y_T \[Pa\].
    pub y_t: f64,
    /// Transverse compressive strength Y_C \[Pa\].
    pub y_c: f64,
    /// In-plane shear strength S_12 \[Pa\].
    pub s12: f64,
}
impl PuckFailureCriteria {
    /// Create a new Puck failure criteria instance.
    pub fn new(x_t: f64, x_c: f64, y_t: f64, y_c: f64, s12: f64) -> Self {
        Self {
            x_t,
            x_c,
            y_t,
            y_c,
            s12,
        }
    }
    /// Fiber failure criterion under longitudinal tension.
    ///
    /// f_FF = (σ_1 + ν_f12·σ_2) / (E_f1 * ε_1f)
    /// Simplified: f_FF = σ_1 / X_T
    pub fn fiber_failure_tension(&self, sigma1: f64, _e_f1: f64, sigma2: f64, _e_f2: f64) -> f64 {
        if sigma1 <= 0.0 {
            return 0.0;
        }
        let correction = 1.0 + sigma2 / self.y_t * 0.01;
        (sigma1 / self.x_t) * correction
    }
    /// Fiber failure criterion under longitudinal compression.
    pub fn fiber_failure_compression(&self, sigma1: f64) -> f64 {
        if sigma1 >= 0.0 {
            return 0.0;
        }
        (-sigma1) / self.x_c
    }
    /// Inter-fiber failure Mode A (transverse tension + shear).
    ///
    /// f_IFF = sqrt((τ_12/S_12)^2 + (1 - p_tt * Y_T / S_12)^2 * (σ_2 / Y_T)^2)
    ///       + p_tt * σ_2 / S_12
    ///
    /// where p_tt is the slope parameter (typically 0.2..0.3).
    pub fn inter_fiber_failure_mode_a(
        &self,
        sigma2: f64,
        tau12: f64,
        p_tt: f64,
        _p_ct: f64,
    ) -> f64 {
        if sigma2 < 0.0 {
            return 0.0;
        }
        let term_tau = (tau12 / self.s12).powi(2);
        let factor = 1.0 - p_tt * self.y_t / self.s12;
        let term_sig = (factor * sigma2 / self.y_t).powi(2);
        (term_tau + term_sig).sqrt() + p_tt * sigma2 / self.s12
    }
    /// Inter-fiber failure Mode B/C (transverse compression + shear).
    pub fn inter_fiber_failure_mode_bc(&self, sigma2: f64, tau12: f64, p_ct: f64) -> f64 {
        if sigma2 >= 0.0 {
            return 0.0;
        }
        let term_tau = (tau12 / self.s12 + p_ct * sigma2 / self.s12).powi(2);
        let term_sig = (sigma2 / self.y_c).powi(2);
        (term_tau + term_sig).sqrt()
    }
    /// Maximum failure index across all modes.
    pub fn max_failure_index(
        &self,
        sigma1: f64,
        sigma2: f64,
        tau12: f64,
        e_f1: f64,
        e_f2: f64,
    ) -> f64 {
        let ff_t = self.fiber_failure_tension(sigma1, e_f1, sigma2, e_f2);
        let ff_c = self.fiber_failure_compression(sigma1);
        let iff_a = self.inter_fiber_failure_mode_a(sigma2, tau12, 0.25, 0.2);
        let iff_bc = self.inter_fiber_failure_mode_bc(sigma2, tau12, 0.2);
        ff_t.max(ff_c).max(iff_a).max(iff_bc)
    }
}
/// Crack band model for mesh-objective softening (Bažant & Oh 1983).
///
/// Regularizes strain-softening by relating the softening curve to fracture
/// energy G_f and the element characteristic length h. This ensures mesh
/// independence of the energy dissipation.
#[derive(Debug, Clone)]
pub struct CrackBandModel {
    /// Fracture energy G_f (J/m² or N/m).
    pub gf: f64,
    /// Tensile strength f_t (Pa).
    pub ft: f64,
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Characteristic element length h (m).
    pub h: f64,
    /// Current damage variable.
    pub d: f64,
    /// Maximum historical strain (for irreversibility).
    pub kappa: f64,
}
impl CrackBandModel {
    /// Create a new crack band model.
    ///
    /// The failure strain is computed as ε_f = 2*G_f / (f_t * h).
    pub fn new(gf: f64, ft: f64, young_modulus: f64, h: f64) -> Self {
        Self {
            gf,
            ft,
            young_modulus,
            h,
            d: 0.0,
            kappa: 0.0,
        }
    }
    /// Damage initiation strain: ε_0 = f_t / E.
    pub fn initiation_strain(&self) -> f64 {
        self.ft / self.young_modulus
    }
    /// Failure strain: ε_f = 2*G_f / (f_t * h).
    ///
    /// This ensures the correct fracture energy is dissipated
    /// regardless of element size.
    pub fn failure_strain(&self) -> f64 {
        2.0 * self.gf / (self.ft * self.h)
    }
    /// Check snap-back condition: ε_f must be > ε_0.
    ///
    /// If h > 2*E*G_f / f_t², snap-back occurs and the element is too large.
    pub fn has_snapback(&self) -> bool {
        self.failure_strain() <= self.initiation_strain()
    }
    /// Maximum allowable element size to avoid snap-back.
    pub fn max_element_size(&self) -> f64 {
        2.0 * self.young_modulus * self.gf / (self.ft * self.ft)
    }
    /// Update damage based on current strain.
    ///
    /// Uses linear softening between ε_0 and ε_f.
    pub fn update(&mut self, strain: f64) {
        if strain <= self.kappa {
            return;
        }
        self.kappa = strain;
        let eps0 = self.initiation_strain();
        let eps_f = self.failure_strain();
        if strain <= eps0 {
            return;
        }
        let d_new = if strain >= eps_f {
            1.0
        } else {
            eps_f / strain * (strain - eps0) / (eps_f - eps0)
        };
        self.d = d_new.clamp(self.d, 1.0);
    }
    /// Stress for a given strain considering damage.
    ///
    /// σ = (1 - D) * E * ε
    pub fn stress(&self, strain: f64) -> f64 {
        (1.0 - self.d) * self.young_modulus * strain
    }
    /// Energy dissipated per unit volume so far.
    ///
    /// For linear softening: g = 0.5 * f_t * (ε_f - ε_0) * D
    /// Total energy dissipated per unit area = g * h ≈ G_f when D=1.
    pub fn dissipated_energy_density(&self) -> f64 {
        let eps0 = self.initiation_strain();
        let eps_f = self.failure_strain();
        0.5 * self.ft * (eps_f - eps0) * self.d
    }
}
/// Gurson-Tvergaard-Needleman (GTN) porous plasticity model.
///
/// Yield function:
/// Φ = (σ_eq/σ_y)² + 2*q₁*f* cosh(3*q₂*σ_H/(2*σ_y)) - 1 - (q₁*f*)² = 0
///
/// where f* is the effective void volume fraction.
#[derive(Debug, Clone)]
pub struct GursonDamage {
    /// Current void volume fraction f ∈ \[0, 1\].
    pub f: f64,
    /// Initial void volume fraction f_0.
    pub f0: f64,
    /// Critical void volume fraction f_c for coalescence.
    pub fc: f64,
    /// Failure void volume fraction f_F.
    pub f_f: f64,
    /// Tvergaard parameter q₁ (typically 1.5).
    pub q1: f64,
    /// Tvergaard parameter q₂ (typically 1.0).
    pub q2: f64,
    /// Tvergaard parameter q₃ (typically q₁²).
    pub q3: f64,
    /// Nucleation amplitude A_n.
    pub a_n: f64,
    /// Mean nucleation strain ε_N.
    pub eps_n: f64,
    /// Standard deviation of nucleation s_N.
    pub s_n: f64,
}
impl GursonDamage {
    /// Create a new GTN model with default Tvergaard parameters.
    pub fn new(f0: f64, fc: f64, f_f: f64, a_n: f64, eps_n: f64, s_n: f64) -> Self {
        let q1 = 1.5;
        let q2 = 1.0;
        Self {
            f: f0,
            f0,
            fc,
            f_f,
            q1,
            q2,
            q3: q1 * q1,
            a_n,
            eps_n,
            s_n,
        }
    }
    /// Effective void volume fraction f* accounting for coalescence.
    ///
    /// f* = f                        if f <= f_c
    /// f* = f_c + (f̄_F - f_c)/(f_F - f_c) * (f - f_c)   if f > f_c
    ///
    /// where f̄_F = (q₁ + sqrt(q₁² - q₃))/q₃ ≈ 1/q₁
    pub fn effective_void_fraction(&self) -> f64 {
        if self.f <= self.fc {
            self.f
        } else {
            let f_bar_f = 1.0 / self.q1;
            let kappa = (f_bar_f - self.fc) / (self.f_f - self.fc);
            self.fc + kappa * (self.f - self.fc)
        }
    }
    /// Evaluate the GTN yield function.
    ///
    /// Φ = (σ_eq/σ_y)² + 2*q₁*f* cosh(3*q₂*σ_H/(2*σ_y)) - 1 - q₃*f*²
    pub fn yield_function(&self, sigma_eq: f64, sigma_h: f64, sigma_y: f64) -> f64 {
        let f_star = self.effective_void_fraction();
        let term1 = (sigma_eq / sigma_y).powi(2);
        let term2 = 2.0 * self.q1 * f_star * (1.5 * self.q2 * sigma_h / sigma_y).cosh();
        let term3 = 1.0 + self.q3 * f_star * f_star;
        term1 + term2 - term3
    }
    /// Void nucleation rate from strain-controlled nucleation (Chu & Needleman).
    ///
    /// ḟ_nucleation = A_n / (s_N * sqrt(2π)) * exp(-0.5*((ε_p - ε_N)/s_N)²) * ε̇_p
    pub fn nucleation_rate(&self, eps_p: f64, eps_p_dot: f64) -> f64 {
        let z = (eps_p - self.eps_n) / self.s_n;
        let coeff = self.a_n / (self.s_n * (2.0 * std::f64::consts::PI).sqrt());
        coeff * (-0.5 * z * z).exp() * eps_p_dot
    }
    /// Void growth rate from plastic dilatation.
    ///
    /// ḟ_growth = (1 - f) * ε̇_kk^p
    pub fn growth_rate(&self, plastic_volumetric_strain_rate: f64) -> f64 {
        (1.0 - self.f) * plastic_volumetric_strain_rate
    }
    /// Update void volume fraction given plastic strain increment.
    pub fn update(
        &mut self,
        eps_p: f64,
        eps_p_dot: f64,
        plastic_volumetric_strain_rate: f64,
        dt: f64,
    ) {
        let f_nucleation = self.nucleation_rate(eps_p, eps_p_dot) * dt;
        let f_growth = self.growth_rate(plastic_volumetric_strain_rate) * dt;
        self.f = (self.f + f_nucleation + f_growth).clamp(0.0, self.f_f);
    }
    /// Check if material has failed (f >= f_F).
    pub fn is_failed(&self) -> bool {
        self.f >= self.f_f - 1e-15
    }
    /// Compute the total void growth rate including nucleation and growth.
    ///
    /// The Rice-Tracey void growth model extended for GTN:
    ///
    /// ḟ_total = ḟ_growth + ḟ_nucleation
    ///
    /// where:
    /// - ḟ_growth = (1 - f) * A * sinh(B * η) * ε̇_p
    /// - ḟ_nucleation = Gaussian nucleation from Chu-Needleman model
    /// - η = σ_H / σ_eq (stress triaxiality)
    /// - A = 0.283 (Rice-Tracey coefficient), B = 1.5 (triaxiality amplifier)
    ///
    /// # Arguments
    /// * `sigma_h`    - Hydrostatic (mean) stress \[Pa\]
    /// * `sigma_eq`   - Von Mises equivalent stress \[Pa\]
    /// * `sigma_y`    - Current yield stress \[Pa\]
    /// * `deps_p`     - Equivalent plastic strain increment (dimensionless)
    ///
    /// # Returns
    /// Total void volume fraction rate df/dε (per unit plastic strain)
    pub fn compute_void_growth_rate(
        &self,
        sigma_h: f64,
        sigma_eq: f64,
        sigma_y: f64,
        deps_p: f64,
    ) -> f64 {
        let triaxiality = if sigma_eq.abs() > 1e-30 {
            sigma_h / sigma_eq
        } else {
            0.0
        };
        let a_rt = 0.283_f64;
        let b_rt = 1.5_f64;
        let growth = (1.0 - self.f) * a_rt * (b_rt * triaxiality).sinh() * deps_p;
        let eps_p_effective = if sigma_eq.abs() > 1e-30 {
            deps_p * sigma_y / sigma_eq
        } else {
            deps_p
        };
        let nucleation = self.nucleation_rate(eps_p_effective, 1.0) * deps_p;
        (growth + nucleation).max(0.0)
    }
}
/// Simple high-cycle fatigue damage accumulation (Miner's rule).
///
/// D = Σ nᵢ / Nᵢ
///
/// where nᵢ is the number of cycles at stress level i and Nᵢ is the
/// fatigue life at that stress level from the S-N curve.
#[derive(Debug, Clone)]
pub struct MinerFatigueDamage {
    /// Accumulated damage D.
    pub d: f64,
    /// S-N curve coefficient C (N = C / S^m).
    pub c: f64,
    /// S-N curve exponent m.
    pub m: f64,
}
impl MinerFatigueDamage {
    /// Create a new Miner fatigue damage model.
    pub fn new(c: f64, m: f64) -> Self {
        Self { d: 0.0, c, m }
    }
    /// Fatigue life N at stress amplitude S from the S-N curve.
    ///
    /// N = C / S^m
    pub fn fatigue_life(&self, stress_amplitude: f64) -> f64 {
        if stress_amplitude.abs() < f64::EPSILON {
            return f64::INFINITY;
        }
        self.c / stress_amplitude.abs().powf(self.m)
    }
    /// Accumulate damage for n cycles at the given stress amplitude.
    pub fn accumulate(&mut self, stress_amplitude: f64, n_cycles: f64) {
        let n_f = self.fatigue_life(stress_amplitude);
        if n_f.is_finite() {
            self.d = (self.d + n_cycles / n_f).min(1.0);
        }
    }
    /// Check if fatigue failure has occurred (D >= 1).
    pub fn is_failed(&self) -> bool {
        self.d >= 1.0 - 1e-15
    }
}
