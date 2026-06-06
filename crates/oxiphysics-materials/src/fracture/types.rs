//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// CTOD (Crack Tip Opening Displacement) fracture criterion.
///
/// Failure is predicted when CTOD reaches a material-dependent critical value.
/// The plastic zone size at the crack tip is computed using the Irwin or Dugdale
/// model.
#[derive(Debug, Clone, Copy)]
pub struct CtodCriterion {
    /// Critical CTOD δ_c (m).
    pub ctod_critical: f64,
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Yield stress σ_ys (Pa).
    pub yield_stress: f64,
}
impl CtodCriterion {
    /// Create a new CTOD criterion model.
    pub fn new(
        ctod_critical: f64,
        young_modulus: f64,
        poisson_ratio: f64,
        yield_stress: f64,
    ) -> Self {
        Self {
            ctod_critical,
            young_modulus,
            poisson_ratio,
            yield_stress,
        }
    }
    /// Structural steel default (δ_c = 0.1 mm, E = 200 GPa, σ_ys = 355 MPa).
    pub fn structural_steel() -> Self {
        Self::new(0.1e-3, 200.0e9, 0.3, 355.0e6)
    }
    /// Compute the Irwin plastic zone radius.
    ///
    /// Plane stress:  r_p = (1/(2π)) * (K_I / σ_ys)²
    /// Plane strain:  r_p = (1/(6π)) * (K_I / σ_ys)²
    pub fn compute_plastic_zone_size(&self, k_i: f64, plane_strain: bool) -> f64 {
        let factor = if plane_strain {
            1.0 / (6.0 * PI)
        } else {
            1.0 / (2.0 * PI)
        };
        factor * (k_i / self.yield_stress).powi(2)
    }
    /// Compute the Dugdale (strip-yield) plastic zone size.
    ///
    /// c = a * (sec(π σ / (2 σ_ys)) - 1)
    pub fn compute_dugdale_plastic_zone(&self, half_crack: f64, applied_stress: f64) -> f64 {
        let arg = PI * applied_stress / (2.0 * self.yield_stress);
        if arg.abs() >= PI / 2.0 - 1e-6 {
            return f64::INFINITY;
        }
        half_crack * (1.0 / arg.cos() - 1.0)
    }
    /// Check whether the CTOD criterion is met (fracture initiates).
    ///
    /// Returns `true` if δ >= δ_c.
    pub fn is_critical(&self, ctod: f64) -> bool {
        ctod >= self.ctod_critical
    }
    /// Compute CTOD from stress intensity factor.
    ///
    /// Plane stress: δ = K_I² / (E σ_ys)
    /// Plane strain: δ = K_I² (1-ν²) / (E σ_ys)
    pub fn ctod_from_sif(&self, k_i: f64, plane_strain: bool) -> f64 {
        let e_prime = if plane_strain {
            self.young_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio)
        } else {
            self.young_modulus
        };
        k_i * k_i / (e_prime * self.yield_stress)
    }
    /// Compute the critical SIF K_Ic from the critical CTOD.
    ///
    /// K_Ic = sqrt(δ_c * E' * σ_ys)
    pub fn critical_sif(&self, plane_strain: bool) -> f64 {
        let e_prime = if plane_strain {
            self.young_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio)
        } else {
            self.young_modulus
        };
        (self.ctod_critical * e_prime * self.yield_stress).sqrt()
    }
}
/// Parameters for the AT2 phase field fracture model.
#[derive(Debug, Clone)]
pub struct PhaseFieldParams {
    /// Critical energy release rate (J/m²)
    pub gc: f64,
    /// Internal length scale (m)
    pub l0: f64,
    /// Young's modulus (Pa)
    pub e: f64,
    /// Poisson's ratio (dimensionless)
    pub nu: f64,
}
impl PhaseFieldParams {
    /// Create new phase field parameters.
    pub fn new(gc: f64, l0: f64, e: f64, nu: f64) -> Self {
        Self { gc, l0, e, nu }
    }
}
/// Parameters for the Gurson-Tvergaard-Needleman ductile fracture model.
#[derive(Debug, Clone)]
pub struct GtnParams {
    /// Tvergaard correction factor q1 (≈1.5)
    pub q1: f64,
    /// Tvergaard correction factor q2 (≈1.0)
    pub q2: f64,
    /// Tvergaard correction factor q3 = q1² (≈2.25)
    pub q3: f64,
    /// Initial void volume fraction
    pub f0: f64,
    /// Critical void fraction at the onset of coalescence
    pub fc: f64,
    /// Void fraction at complete fracture
    pub ff: f64,
    /// Nucleation void volume fraction amplitude
    pub fn_void: f64,
    /// Standard deviation of the nucleation strain distribution
    pub sn: f64,
    /// Mean nucleation strain
    pub en: f64,
    /// Uniaxial yield stress (Pa)
    pub sigma_y: f64,
    /// Young's modulus (Pa)
    pub e: f64,
    /// Poisson's ratio (dimensionless)
    pub nu: f64,
}
impl GtnParams {
    /// Typical GTN parameters for structural steel.
    pub fn steel() -> Self {
        Self {
            q1: 1.5,
            q2: 1.0,
            q3: 2.25,
            f0: 0.001,
            fc: 0.15,
            ff: 0.25,
            fn_void: 0.04,
            sn: 0.1,
            en: 0.3,
            sigma_y: 250.0e6,
            e: 200.0e9,
            nu: 0.3,
        }
    }
}
/// Griffith fracture criterion for brittle materials.
///
/// Fracture occurs when the strain energy release rate G equals or exceeds
/// the critical energy release rate G_c.
///
/// For mode I:   G = K_I² / E'
/// For mode II:  G = K_II² / E'
/// For mode III: G = K_III² / (2 μ)
///
/// where E' = E for plane stress and E' = E / (1-ν²) for plane strain.
#[derive(Debug, Clone, Copy)]
pub struct GriffithCriterion {
    /// Critical energy release rate G_c (J/m²).
    pub g_critical: f64,
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Shear modulus G (Pa) = E / (2(1+ν)).
    pub shear_modulus: f64,
}
impl GriffithCriterion {
    /// Create a Griffith criterion model.
    pub fn new(g_critical: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        let shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
        Self {
            g_critical,
            young_modulus,
            poisson_ratio,
            shear_modulus,
        }
    }
    /// Steel default (G_c = 100 J/m², E = 200 GPa, ν = 0.3).
    pub fn steel() -> Self {
        Self::new(100.0, 200.0e9, 0.3)
    }
    /// Effective modulus E' for plane stress or plane strain.
    pub fn effective_modulus(&self, plane_strain: bool) -> f64 {
        if plane_strain {
            self.young_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio)
        } else {
            self.young_modulus
        }
    }
    /// Compute the energy release rate for mode I loading.
    ///
    /// G_I = K_I² / E'
    pub fn compute_energy_release_rate_mode_i(&self, k_i: f64, plane_strain: bool) -> f64 {
        let e_prime = self.effective_modulus(plane_strain);
        k_i * k_i / e_prime
    }
    /// Compute the energy release rate for mixed mode (I + II + III).
    ///
    /// G = K_I²/E' + K_II²/E' + K_III²/(2μ)
    pub fn compute_energy_release_rate(
        &self,
        k_i: f64,
        k_ii: f64,
        k_iii: f64,
        plane_strain: bool,
    ) -> f64 {
        let e_prime = self.effective_modulus(plane_strain);
        let g_mode1 = k_i * k_i / e_prime;
        let g_mode2 = k_ii * k_ii / e_prime;
        let g_mode3 = k_iii * k_iii / (2.0 * self.shear_modulus);
        g_mode1 + g_mode2 + g_mode3
    }
    /// Check whether the Griffith criterion is satisfied (fracture occurs).
    ///
    /// Returns `true` if G >= G_c.
    pub fn is_critical(&self, k_i: f64, k_ii: f64, k_iii: f64, plane_strain: bool) -> bool {
        self.compute_energy_release_rate(k_i, k_ii, k_iii, plane_strain) >= self.g_critical
    }
    /// Critical SIF K_Ic from G_c (mode I, plane strain).
    ///
    /// K_Ic = sqrt(G_c * E')
    pub fn critical_sif(&self, plane_strain: bool) -> f64 {
        let e_prime = self.effective_modulus(plane_strain);
        (self.g_critical * e_prime).sqrt()
    }
    /// Compute the critical half-crack length a_c for a given applied stress.
    ///
    /// a_c = G_c * E' / (π * σ²)
    pub fn critical_crack_length(&self, stress: f64, plane_strain: bool) -> f64 {
        let e_prime = self.effective_modulus(plane_strain);
        if stress.abs() < f64::EPSILON {
            return f64::INFINITY;
        }
        self.g_critical * e_prime / (PI * stress * stress)
    }
}
/// R-curve model: fracture resistance as a function of crack extension.
///
/// Models the rising resistance to crack growth due to crack-wake processes
/// (bridging, plastic zone development, etc.).
#[derive(Debug, Clone)]
pub struct RCurve {
    /// Initial fracture toughness K_Ic or G_Ic (at initiation).
    pub k_init: f64,
    /// Plateau (steady-state) fracture toughness.
    pub k_plateau: f64,
    /// Characteristic length for the transition (m).
    pub delta_a_char: f64,
}
impl RCurve {
    /// Create a new R-curve model.
    pub fn new(k_init: f64, k_plateau: f64, delta_a_char: f64) -> Self {
        Self {
            k_init,
            k_plateau,
            delta_a_char,
        }
    }
    /// Evaluate the R-curve at a given crack extension Δa.
    ///
    /// K_R(Δa) = K_init + (K_plateau - K_init) * (1 - exp(-Δa / Δa_char))
    pub fn resistance(&self, delta_a: f64) -> f64 {
        if delta_a <= 0.0 {
            return self.k_init;
        }
        self.k_init + (self.k_plateau - self.k_init) * (1.0 - (-delta_a / self.delta_a_char).exp())
    }
    /// Check crack stability: crack grows stably if K_applied ≤ K_R(Δa).
    ///
    /// Returns `true` if the crack is stable (no unstable fracture).
    pub fn is_stable(&self, k_applied: f64, delta_a: f64) -> bool {
        k_applied <= self.resistance(delta_a)
    }
    /// Find the critical crack extension at which instability occurs.
    ///
    /// For a linearly increasing applied K: K = K0 + dK/da * Δa,
    /// find Δa where the driving force tangent equals the R-curve tangent.
    ///
    /// Returns the critical crack extension, or None if always stable.
    pub fn critical_extension(&self, k0: f64, dk_da: f64) -> Option<f64> {
        let delta_k = self.k_plateau - self.k_init;
        if delta_k <= 0.0 || dk_da <= 0.0 {
            if k0 >= self.k_init {
                return Some(0.0);
            }
            return None;
        }
        let slope_init = delta_k / self.delta_a_char;
        if dk_da >= slope_init {
            return Some(0.0);
        }
        let delta_a_crit = -self.delta_a_char * (dk_da / slope_init).ln();
        let k_r = self.resistance(delta_a_crit);
        let k_app = k0 + dk_da * delta_a_crit;
        if k_app >= k_r {
            Some(delta_a_crit)
        } else {
            None
        }
    }
}
impl RCurve {
    /// Simulate stable crack growth under a monotonically increasing applied K.
    ///
    /// Starting from `a0` with applied K = `k0 + dk_da * Δa`, integrate the
    /// stable crack extension Δa until either:
    /// - instability (dK_R/dΔa < dK_applied/dΔa), or
    /// - the maximum crack extension `delta_a_max` is reached.
    ///
    /// Returns a vector of (Δa, K_applied, K_R) tuples.
    pub fn compute_stable_crack_growth(
        &self,
        k0: f64,
        dk_da: f64,
        delta_a_max: f64,
        n_steps: usize,
    ) -> Vec<(f64, f64, f64)> {
        if n_steps == 0 || delta_a_max <= 0.0 {
            return Vec::new();
        }
        let step = delta_a_max / n_steps as f64;
        let mut result = Vec::with_capacity(n_steps + 1);
        for i in 0..=n_steps {
            let da = i as f64 * step;
            let k_app = k0 + dk_da * da;
            let k_r = self.resistance(da);
            result.push((da, k_app, k_r));
            if k_app > k_r + 1e-6 * k_r {
                break;
            }
        }
        result
    }
    /// Minimum applied K0 needed to initiate stable crack growth.
    ///
    /// Crack initiates when K_applied >= K_init (= K_R at Δa = 0).
    pub fn initiation_load(&self) -> f64 {
        self.k_init
    }
    /// Fracture instability load: K at which dK/da = dK_R/dΔa.
    ///
    /// For a linear loading: K = K0 + (dK/da) * Δa with slope `dk_da`,
    /// returns the K level at the tangency point, or None if always stable.
    pub fn instability_load(&self, dk_da: f64) -> Option<f64> {
        let delta_a = self.critical_extension(self.k_init, dk_da)?;
        Some(self.k_init + dk_da * delta_a)
    }
}
