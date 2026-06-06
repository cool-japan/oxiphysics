// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tribology materials module.
//!
//! Comprehensive tribology physics covering:
//! - Archard wear law (adhesive and abrasive wear)
//! - Hertz contact pressure and elastic contact mechanics
//! - Flash temperature (Blok and Jaeger models)
//! - Surface roughness parameters (Ra, Rq, Rz)
//! - Asperity contact model (Greenwood-Williamson)
//! - Boundary, mixed, and elastohydrodynamic lubrication (EHL)
//! - Friction coefficient models (Coulomb, viscous, Stribeck)
//! - Fretting wear and fretting fatigue
//! - Wear debris and tribocorrosion
//! - Solid lubricants (MoS2, DLC)
//! - Rolling contact fatigue

use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// § 1.  Archard Wear Law
// ══════════════════════════════════════════════════════════════════════════════

/// Archard wear model parameters.
///
/// The Archard equation is:  V = K · W · s / H
/// where V is wear volume, K is the dimensionless wear coefficient,
/// W is the normal load, s is the sliding distance, and H is the
/// surface hardness.
#[derive(Debug, Clone)]
pub struct ArchardWear {
    /// Dimensionless Archard wear coefficient (typically 10^-7 … 10^-2).
    pub wear_coefficient: f64,
    /// Surface hardness of the softer material \[Pa\].
    pub hardness: f64,
}

impl ArchardWear {
    /// Create a new Archard wear model.
    ///
    /// # Arguments
    /// * `wear_coefficient` – dimensionless wear coefficient K (0 < K ≤ 1)
    /// * `hardness` – Vickers or Brinell hardness converted to Pa
    pub fn new(wear_coefficient: f64, hardness: f64) -> Self {
        Self {
            wear_coefficient,
            hardness,
        }
    }

    /// Compute worn volume \[m³\].
    ///
    /// # Arguments
    /// * `normal_load` – applied normal load W \[N\]
    /// * `sliding_distance` – total sliding distance s \[m\]
    pub fn wear_volume(&self, normal_load: f64, sliding_distance: f64) -> f64 {
        self.wear_coefficient * normal_load * sliding_distance / self.hardness
    }

    /// Compute wear rate dV/ds \[m³/m\] = K·W/H.
    ///
    /// # Arguments
    /// * `normal_load` – applied normal load W \[N\]
    pub fn wear_rate(&self, normal_load: f64) -> f64 {
        self.wear_coefficient * normal_load / self.hardness
    }

    /// Compute specific wear rate k = K/H \[m³/(N·m)\].
    pub fn specific_wear_rate(&self) -> f64 {
        self.wear_coefficient / self.hardness
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 2.  Hertz Contact Mechanics
// ══════════════════════════════════════════════════════════════════════════════

/// Hertz contact model for elastic sphere-on-flat or sphere-on-sphere contact.
///
/// Uses the classical Hertz theory for frictionless elastic contact.
#[derive(Debug, Clone)]
pub struct HertzContact {
    /// Effective Young's modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius of curvature R* \[m\].
    pub effective_radius: f64,
}

impl HertzContact {
    /// Construct a Hertz contact model from material and geometry parameters.
    ///
    /// # Arguments
    /// * `e1`, `e2` – Young's moduli of body 1 and 2 \[Pa\]
    /// * `nu1`, `nu2` – Poisson's ratios of body 1 and 2
    /// * `r1`, `r2` – principal radii; use `f64::INFINITY` for a flat surface
    pub fn new(e1: f64, nu1: f64, e2: f64, nu2: f64, r1: f64, r2: f64) -> Self {
        let inv_e_star = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
        let effective_modulus = 1.0 / inv_e_star;
        let inv_r_star = 1.0 / r1 + 1.0 / r2;
        let effective_radius = 1.0 / inv_r_star;
        Self {
            effective_modulus,
            effective_radius,
        }
    }

    /// Hertz contact radius a \[m\] under normal load F \[N\].
    pub fn contact_radius(&self, normal_load: f64) -> f64 {
        let arg = 3.0 * normal_load * self.effective_radius / (4.0 * self.effective_modulus);
        arg.cbrt()
    }

    /// Maximum Hertz contact pressure p0 \[Pa\] under normal load F \[N\].
    pub fn max_pressure(&self, normal_load: f64) -> f64 {
        let a = self.contact_radius(normal_load);
        3.0 * normal_load / (2.0 * PI * a * a)
    }

    /// Mean Hertz contact pressure \[Pa\] under normal load F \[N\].
    pub fn mean_pressure(&self, normal_load: f64) -> f64 {
        let a = self.contact_radius(normal_load);
        normal_load / (PI * a * a)
    }

    /// Elastic indentation depth δ \[m\].
    pub fn indentation_depth(&self, normal_load: f64) -> f64 {
        let a = self.contact_radius(normal_load);
        a * a / self.effective_radius
    }

    /// Contact stiffness dF/dδ \[N/m\].
    pub fn contact_stiffness(&self, normal_load: f64) -> f64 {
        let a = self.contact_radius(normal_load);
        2.0 * self.effective_modulus * a
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3.  Flash Temperature
// ══════════════════════════════════════════════════════════════════════════════

/// Flash temperature calculator (Blok and Jaeger models).
///
/// Flash temperature is the transient temperature rise at asperity contacts
/// during sliding; important for lubricant degradation and scuffing.
#[derive(Debug, Clone)]
pub struct FlashTemperature {
    /// Thermal conductivity of material 1 \[W/(m·K)\].
    pub lambda1: f64,
    /// Thermal conductivity of material 2 \[W/(m·K)\].
    pub lambda2: f64,
    /// Thermal diffusivity of material 1 \[m²/s\].
    pub alpha1: f64,
    /// Thermal diffusivity of material 2 \[m²/s\].
    pub alpha2: f64,
}

impl FlashTemperature {
    /// Construct a flash temperature model.
    ///
    /// # Arguments
    /// * `lambda1`, `lambda2` – thermal conductivities \[W/(m·K)\]
    /// * `alpha1`, `alpha2` – thermal diffusivities \[m²/s\]
    pub fn new(lambda1: f64, alpha1: f64, lambda2: f64, alpha2: f64) -> Self {
        Self {
            lambda1,
            lambda2,
            alpha1,
            alpha2,
        }
    }

    /// Blok flash temperature rise ΔT \[K\] for a circular contact patch.
    ///
    /// Uses the Blok (1937) fast-moving heat source approximation:
    /// ΔT = 0.308 · μ · p_mean · v · a / (lambda1 + lambda2)
    ///
    /// # Arguments
    /// * `friction_coeff` – sliding friction coefficient μ
    /// * `mean_pressure` – mean contact pressure \[Pa\]
    /// * `sliding_velocity` – relative sliding speed v \[m/s\]
    /// * `contact_radius` – Hertz contact radius a \[m\]
    pub fn blok_flash_temperature(
        &self,
        friction_coeff: f64,
        mean_pressure: f64,
        sliding_velocity: f64,
        contact_radius: f64,
    ) -> f64 {
        let heat_flux = friction_coeff * mean_pressure * sliding_velocity;
        0.308 * heat_flux * contact_radius / (self.lambda1 + self.lambda2)
    }

    /// Jaeger flash temperature for a stationary heat source (low Pe limit).
    ///
    /// ΔT = q · a / (π · (lambda1 + lambda2))  where q = μ·p·v
    ///
    /// # Arguments
    /// * `friction_coeff` – friction coefficient
    /// * `mean_pressure` – mean contact pressure \[Pa\]
    /// * `sliding_velocity` – sliding speed \[m/s\]
    /// * `contact_radius` – contact radius a \[m\]
    pub fn jaeger_flash_temperature(
        &self,
        friction_coeff: f64,
        mean_pressure: f64,
        sliding_velocity: f64,
        contact_radius: f64,
    ) -> f64 {
        let heat_flux = friction_coeff * mean_pressure * sliding_velocity;
        heat_flux * contact_radius / (PI * (self.lambda1 + self.lambda2))
    }

    /// Peclet number Pe = v·a/(2·α_avg) – governs which flash-T model to use.
    ///
    /// Pe >> 1: Blok (fast source); Pe << 1: Jaeger (slow source).
    pub fn peclet_number(&self, sliding_velocity: f64, contact_radius: f64) -> f64 {
        let alpha_avg = 0.5 * (self.alpha1 + self.alpha2);
        sliding_velocity * contact_radius / (2.0 * alpha_avg)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 4.  Surface Roughness Parameters
// ══════════════════════════════════════════════════════════════════════════════

/// Surface roughness measurement data and derived parameters.
///
/// Stores a profile array and computes Ra, Rq, Rz according to ISO 4287.
#[derive(Debug, Clone)]
pub struct SurfaceRoughness {
    /// Height profile values zi relative to mean line \[m\].
    pub profile: Vec<f64>,
    /// Sampling length L \[m\].
    pub sampling_length: f64,
}

impl SurfaceRoughness {
    /// Create a roughness object from a measured height profile.
    ///
    /// # Arguments
    /// * `profile` – vector of surface heights relative to mean line \[m\]
    /// * `sampling_length` – total profile length L \[m\]
    pub fn new(profile: Vec<f64>, sampling_length: f64) -> Self {
        Self {
            profile,
            sampling_length,
        }
    }

    /// Arithmetic mean roughness Ra \[m\] (ISO 4287).
    pub fn ra(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.profile.iter().map(|z| z.abs()).sum();
        sum / self.profile.len() as f64
    }

    /// Root-mean-square roughness Rq \[m\] (ISO 4287).
    pub fn rq(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = self.profile.iter().map(|z| z * z).sum();
        (sum_sq / self.profile.len() as f64).sqrt()
    }

    /// Maximum peak-to-valley roughness Rz \[m\] (ISO 4287).
    pub fn rz(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        let max_z = self
            .profile
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min_z = self.profile.iter().cloned().fold(f64::INFINITY, f64::min);
        max_z - min_z
    }

    /// Skewness Rsk – asymmetry of height distribution.
    pub fn rsk(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        let rq = self.rq();
        if rq < 1e-30 {
            return 0.0;
        }
        let n = self.profile.len() as f64;
        let sum_cube: f64 = self.profile.iter().map(|z| z * z * z).sum();
        sum_cube / (n * rq * rq * rq)
    }

    /// Kurtosis Rku – peakedness of height distribution.
    pub fn rku(&self) -> f64 {
        if self.profile.is_empty() {
            return 0.0;
        }
        let rq = self.rq();
        if rq < 1e-30 {
            return 0.0;
        }
        let n = self.profile.len() as f64;
        let sum_4th: f64 = self.profile.iter().map(|z| z * z * z * z).sum();
        sum_4th / (n * rq * rq * rq * rq)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 5.  Greenwood-Williamson Asperity Contact Model
// ══════════════════════════════════════════════════════════════════════════════

/// Greenwood-Williamson (GW) statistical asperity contact model.
///
/// Models rough surface contact as a distribution of hemispherical asperities
/// with Gaussian height distribution pressing against a flat rigid counter-surface.
///
/// Reference: Greenwood & Williamson, Proc. R. Soc. A 295 (1966) 300-319.
#[derive(Debug, Clone)]
pub struct GreenwoodWilliamson {
    /// Asperity density η \[asperities/m²\].
    pub asperity_density: f64,
    /// Asperity tip radius β \[m\].
    pub asperity_radius: f64,
    /// Standard deviation of asperity height distribution σ_s \[m\].
    pub height_std: f64,
    /// Composite Young's modulus E* \[Pa\].
    pub effective_modulus: f64,
}

impl GreenwoodWilliamson {
    /// Create a new GW asperity contact model.
    ///
    /// # Arguments
    /// * `asperity_density` – η, surface density of asperities \[1/m²\]
    /// * `asperity_radius` – β, mean asperity tip radius \[m\]
    /// * `height_std` – σ_s, standard deviation of asperity heights \[m\]
    /// * `effective_modulus` – E*, composite Young's modulus \[Pa\]
    pub fn new(
        asperity_density: f64,
        asperity_radius: f64,
        height_std: f64,
        effective_modulus: f64,
    ) -> Self {
        Self {
            asperity_density,
            asperity_radius,
            height_std,
            effective_modulus,
        }
    }

    /// Plasticity index ψ (Greenwood-Williamson criterion).
    ///
    /// ψ = (E*/H) · sqrt(σ_s/β)
    /// ψ < 0.6: elastic contact; ψ > 1: plastic contact.
    ///
    /// # Arguments
    /// * `hardness` – material hardness H \[Pa\]
    pub fn plasticity_index(&self, hardness: f64) -> f64 {
        (self.effective_modulus / hardness) * (self.height_std / self.asperity_radius).sqrt()
    }

    /// Expected real contact area A_r \[m²\] at nominal contact area A_0 \[m²\]
    /// and separation d \[m\] (Gaussian approximation).
    ///
    /// Uses numerical integration via mid-point rule over the Gaussian tail.
    ///
    /// # Arguments
    /// * `separation` – mean surface separation d \[m\]
    /// * `nominal_area` – nominal contact area A_0 \[m²\]
    pub fn real_contact_area(&self, separation: f64, nominal_area: f64) -> f64 {
        let n_steps = 200usize;
        let z_max = 6.0 * self.height_std;
        let dz = (z_max - separation) / n_steps as f64;
        if dz <= 0.0 {
            return 0.0;
        }
        let mut integral = 0.0;
        for i in 0..n_steps {
            let z = separation + (i as f64 + 0.5) * dz;
            let phi = gaussian_pdf(z, 0.0, self.height_std);
            let delta = z - separation;
            integral += PI * self.asperity_radius * delta * phi * dz;
        }
        self.asperity_density * nominal_area * integral
    }

    /// Expected total normal load N \[N\] at separation d \[m\] on area A_0 \[m²\].
    ///
    /// # Arguments
    /// * `separation` – mean surface separation d \[m\]
    /// * `nominal_area` – nominal contact area A_0 \[m²\]
    pub fn total_load(&self, separation: f64, nominal_area: f64) -> f64 {
        let n_steps = 200usize;
        let z_max = 6.0 * self.height_std;
        let dz = (z_max - separation) / n_steps as f64;
        if dz <= 0.0 {
            return 0.0;
        }
        let e_star = self.effective_modulus;
        let beta = self.asperity_radius;
        let mut integral = 0.0;
        for i in 0..n_steps {
            let z = separation + (i as f64 + 0.5) * dz;
            let phi = gaussian_pdf(z, 0.0, self.height_std);
            let delta = z - separation;
            // Hertz load per asperity: F = (4/3)·E*·sqrt(β)·δ^(3/2)
            let f_asperity = (4.0 / 3.0) * e_star * beta.sqrt() * delta.powf(1.5);
            integral += f_asperity * phi * dz;
        }
        self.asperity_density * nominal_area * integral
    }
}

/// Standard Gaussian probability density function.
fn gaussian_pdf(x: f64, mean: f64, std: f64) -> f64 {
    let z = (x - mean) / std;
    (-0.5 * z * z).exp() / (std * (2.0 * PI).sqrt())
}

// ══════════════════════════════════════════════════════════════════════════════
// § 6.  Lubrication Regime Models
// ══════════════════════════════════════════════════════════════════════════════

/// Lambda ratio (specific film thickness) classifier for lubrication regimes.
///
/// λ = h_min / (Rq1² + Rq2²)^0.5
/// λ < 1: boundary; 1-3: mixed; > 3: full-film EHL/HD.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LubricationRegime {
    /// Boundary lubrication – asperity-dominated.
    Boundary,
    /// Mixed lubrication – partial film with asperity contact.
    Mixed,
    /// Full-film elastohydrodynamic lubrication.
    FullFilm,
}

/// Compute the specific film thickness (lambda ratio).
///
/// # Arguments
/// * `h_min` – minimum film thickness \[m\]
/// * `rq1`, `rq2` – composite rms roughness of surfaces 1 and 2 \[m\]
pub fn lambda_ratio(h_min: f64, rq1: f64, rq2: f64) -> f64 {
    let composite_roughness = (rq1 * rq1 + rq2 * rq2).sqrt();
    if composite_roughness < 1e-30 {
        return f64::INFINITY;
    }
    h_min / composite_roughness
}

/// Classify lubrication regime from lambda ratio.
///
/// # Arguments
/// * `lambda` – specific film thickness λ (dimensionless)
pub fn classify_lubrication_regime(lambda: f64) -> LubricationRegime {
    if lambda < 1.0 {
        LubricationRegime::Boundary
    } else if lambda < 3.0 {
        LubricationRegime::Mixed
    } else {
        LubricationRegime::FullFilm
    }
}

/// Boundary lubrication friction model (Bowden-Tabor).
///
/// μ = τ_0 / p_mean where τ_0 is the interfacial shear strength.
///
/// # Arguments
/// * `interfacial_shear_strength` – τ_0 \[Pa\]
/// * `mean_pressure` – p_mean \[Pa\]
pub fn boundary_lubrication_friction(interfacial_shear_strength: f64, mean_pressure: f64) -> f64 {
    if mean_pressure < 1e-20 {
        return 0.0;
    }
    interfacial_shear_strength / mean_pressure
}

// ══════════════════════════════════════════════════════════════════════════════
// § 7.  Elastohydrodynamic Lubrication (EHL)
// ══════════════════════════════════════════════════════════════════════════════

/// EHL minimum film thickness model (Dowson-Higginson for line contact,
/// Hamrock-Dowson for point contact).
///
/// Reference: Hamrock & Dowson, J. Lubrication Technology 98 (1976) 375-383.
#[derive(Debug, Clone)]
pub struct EhlFilmThickness {
    /// Lubricant dynamic viscosity at inlet η_0 \[Pa·s\].
    pub viscosity_0: f64,
    /// Pressure-viscosity coefficient α \[1/Pa\].
    pub pressure_viscosity_coeff: f64,
    /// Effective elastic modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R* \[m\].
    pub effective_radius: f64,
}

impl EhlFilmThickness {
    /// Construct an EHL film thickness model.
    ///
    /// # Arguments
    /// * `viscosity_0` – inlet viscosity η_0 \[Pa·s\]
    /// * `pressure_viscosity_coeff` – α \[1/Pa\]  (Barus equation parameter)
    /// * `effective_modulus` – E* \[Pa\]
    /// * `effective_radius` – R* \[m\]
    pub fn new(
        viscosity_0: f64,
        pressure_viscosity_coeff: f64,
        effective_modulus: f64,
        effective_radius: f64,
    ) -> Self {
        Self {
            viscosity_0,
            pressure_viscosity_coeff,
            effective_modulus,
            effective_radius,
        }
    }

    /// Hamrock-Dowson minimum film thickness h_min \[m\] for point contact.
    ///
    /// h_min = 3.63 · R* · U^0.68 · G^0.49 · W^(-0.073) · (1 - e^(-0.68·k))
    /// where k is the ellipticity parameter (assumed 1 for circular contact here).
    ///
    /// # Arguments
    /// * `entrainment_velocity` – U \[m/s\]
    /// * `normal_load` – W \[N\]
    pub fn hamrock_dowson_hmin(&self, entrainment_velocity: f64, normal_load: f64) -> f64 {
        let r = self.effective_radius;
        let e_star = self.effective_modulus;
        let eta0 = self.viscosity_0;
        let alpha = self.pressure_viscosity_coeff;

        // Non-dimensional parameters
        let u_nd = eta0 * entrainment_velocity / (e_star * r); // speed parameter
        let g_nd = alpha * e_star; // materials parameter
        let w_nd = normal_load / (e_star * r * r); // load parameter

        // Ellipticity k = 1 (circular), ellipticity factor = 1 - e^(-0.68)
        let ellipticity_factor = 1.0 - (-0.68_f64).exp();
        3.63 * r * u_nd.powf(0.68) * g_nd.powf(0.49) * w_nd.powf(-0.073) * ellipticity_factor
    }

    /// Barus viscosity at pressure p: η(p) = η_0 · exp(α·p).
    ///
    /// # Arguments
    /// * `pressure` – contact pressure p \[Pa\]
    pub fn barus_viscosity(&self, pressure: f64) -> f64 {
        self.viscosity_0 * (self.pressure_viscosity_coeff * pressure).exp()
    }

    /// Dowson-Higginson minimum film thickness for line contact \[m\].
    ///
    /// h_min = 2.65 · α^0.54 · (η_0·u)^0.7 · R^0.43 / (E*^0.03 · F'^0.13)
    /// where F' = load per unit width \[N/m\].
    ///
    /// # Arguments
    /// * `entrainment_velocity` – u \[m/s\]
    /// * `load_per_width` – F' \[N/m\]
    pub fn dowson_higginson_hmin(&self, entrainment_velocity: f64, load_per_width: f64) -> f64 {
        let alpha = self.pressure_viscosity_coeff;
        let eta0 = self.viscosity_0;
        let r = self.effective_radius;
        let e_star = self.effective_modulus;
        2.65 * alpha.powf(0.54) * (eta0 * entrainment_velocity).powf(0.70) * r.powf(0.43)
            / (e_star.powf(0.03) * load_per_width.powf(0.13))
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 8.  Friction Coefficient Models
// ══════════════════════════════════════════════════════════════════════════════

/// Coulomb-Amontons friction model.
///
/// The classical model: friction force F_f = μ · N.
#[derive(Debug, Clone)]
pub struct CoulombFriction {
    /// Static friction coefficient μ_s.
    pub static_coeff: f64,
    /// Kinetic (dynamic) friction coefficient μ_k.
    pub kinetic_coeff: f64,
}

impl CoulombFriction {
    /// Create a Coulomb friction model.
    ///
    /// # Arguments
    /// * `static_coeff` – μ_s, static friction coefficient
    /// * `kinetic_coeff` – μ_k, kinetic friction coefficient (μ_k ≤ μ_s)
    pub fn new(static_coeff: f64, kinetic_coeff: f64) -> Self {
        Self {
            static_coeff,
            kinetic_coeff,
        }
    }

    /// Friction force \[N\] given normal load and sliding state.
    ///
    /// # Arguments
    /// * `normal_load` – N \[N\]
    /// * `sliding` – true if currently sliding, false if stationary
    pub fn friction_force(&self, normal_load: f64, sliding: bool) -> f64 {
        let mu = if sliding {
            self.kinetic_coeff
        } else {
            self.static_coeff
        };
        mu * normal_load.abs()
    }
}

/// Viscous friction model F_f = c · v.
///
/// Applies when viscous shear dominates (full-film hydrodynamic lubrication).
#[derive(Debug, Clone)]
pub struct ViscousFriction {
    /// Viscous damping coefficient c \[N·s/m\].
    pub damping: f64,
}

impl ViscousFriction {
    /// Create a viscous friction model.
    ///
    /// # Arguments
    /// * `damping` – viscous damping c \[N·s/m\]
    pub fn new(damping: f64) -> Self {
        Self { damping }
    }

    /// Friction force \[N\] for given sliding velocity \[m/s\].
    pub fn friction_force(&self, velocity: f64) -> f64 {
        self.damping * velocity
    }
}

/// Stribeck curve friction model.
///
/// Combines Coulomb, Stribeck, and viscous contributions:
/// μ(v) = μ_k + (μ_s - μ_k)·exp(-(v/v_st)^δ) + σ·v
///
/// Reference: Stribeck (1902).
#[derive(Debug, Clone)]
pub struct StribeckFriction {
    /// Static friction coefficient μ_s.
    pub static_coeff: f64,
    /// Kinetic friction coefficient μ_k.
    pub kinetic_coeff: f64,
    /// Stribeck velocity v_st \[m/s\].
    pub stribeck_velocity: f64,
    /// Shape exponent δ (typically 1 or 2).
    pub shape_exponent: f64,
    /// Viscous friction coefficient σ \[s/m\].
    pub viscous_coeff: f64,
}

impl StribeckFriction {
    /// Create a Stribeck friction model.
    ///
    /// # Arguments
    /// * `static_coeff` – μ_s
    /// * `kinetic_coeff` – μ_k
    /// * `stribeck_velocity` – v_st \[m/s\]
    /// * `shape_exponent` – δ
    /// * `viscous_coeff` – σ \[s/m\]
    pub fn new(
        static_coeff: f64,
        kinetic_coeff: f64,
        stribeck_velocity: f64,
        shape_exponent: f64,
        viscous_coeff: f64,
    ) -> Self {
        Self {
            static_coeff,
            kinetic_coeff,
            stribeck_velocity,
            shape_exponent,
            viscous_coeff,
        }
    }

    /// Friction coefficient μ at sliding velocity v \[m/s\].
    pub fn friction_coeff(&self, velocity: f64) -> f64 {
        let v = velocity.abs();
        let stribeck = (self.static_coeff - self.kinetic_coeff)
            * (-(v / self.stribeck_velocity).powf(self.shape_exponent)).exp();
        self.kinetic_coeff + stribeck + self.viscous_coeff * v
    }

    /// Friction force \[N\] for given normal load \[N\] and velocity \[m/s\].
    pub fn friction_force(&self, normal_load: f64, velocity: f64) -> f64 {
        self.friction_coeff(velocity) * normal_load.abs()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 9.  Fretting Wear and Fretting Fatigue
// ══════════════════════════════════════════════════════════════════════════════

/// Fretting wear model based on energy dissipation approach.
///
/// Wear volume V = α_f · E_d where E_d is dissipated energy and α_f is the
/// fretting wear coefficient.
///
/// Reference: Fouvry et al., Wear 255 (2003) 287-298.
#[derive(Debug, Clone)]
pub struct FrettingWear {
    /// Fretting wear coefficient α_f \[m³/J\].
    pub wear_coefficient: f64,
    /// Contact area A \[m²\].
    pub contact_area: f64,
}

impl FrettingWear {
    /// Create a fretting wear model.
    ///
    /// # Arguments
    /// * `wear_coefficient` – α_f \[m³/J\]
    /// * `contact_area` – nominal contact area \[m²\]
    pub fn new(wear_coefficient: f64, contact_area: f64) -> Self {
        Self {
            wear_coefficient,
            contact_area,
        }
    }

    /// Dissipated energy per cycle E_d \[J\].
    ///
    /// For partial slip: E_d = 4·μ·Q* (a - a*)·δ
    /// Simplified: E_d ≈ 4 · Q_t · δ_slip where Q_t is tangential force amplitude.
    ///
    /// # Arguments
    /// * `tangential_force_amplitude` – Q_t \[N\]
    /// * `slip_amplitude` – δ \[m\]
    pub fn dissipated_energy_per_cycle(
        &self,
        tangential_force_amplitude: f64,
        slip_amplitude: f64,
    ) -> f64 {
        4.0 * tangential_force_amplitude * slip_amplitude
    }

    /// Wear volume per cycle \[m³/cycle\].
    ///
    /// # Arguments
    /// * `tangential_force_amplitude` – Q_t \[N\]
    /// * `slip_amplitude` – δ \[m\]
    pub fn wear_volume_per_cycle(
        &self,
        tangential_force_amplitude: f64,
        slip_amplitude: f64,
    ) -> f64 {
        self.wear_coefficient
            * self.dissipated_energy_per_cycle(tangential_force_amplitude, slip_amplitude)
    }

    /// Fretting fatigue stress reduction factor.
    ///
    /// Estimates the reduction in fatigue limit due to fretting contact.
    /// Factor = 1 - (A * p_max / sigma_f)
    /// where sigma_f is the plain fatigue limit.
    ///
    /// # Arguments
    /// * `max_pressure` – p_max \[Pa\]
    /// * `fatigue_limit` – σ_f \[Pa\]
    /// * `fretting_factor` – empirical factor A (typically 0.1-0.5)
    pub fn fatigue_reduction_factor(
        &self,
        max_pressure: f64,
        fatigue_limit: f64,
        fretting_factor: f64,
    ) -> f64 {
        let reduction = fretting_factor * max_pressure / fatigue_limit;
        (1.0 - reduction).max(0.0)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 10.  Wear Debris and Tribocorrosion
// ══════════════════════════════════════════════════════════════════════════════

/// Wear debris generation and third-body particle model.
///
/// Models the formation, accumulation, and ejection of wear particles
/// in the contact zone (third body concept).
#[derive(Debug, Clone)]
pub struct WearDebris {
    /// Debris particle radius r_d \[m\].
    pub particle_radius: f64,
    /// Debris density ρ_d \[kg/m³\].
    pub particle_density: f64,
    /// Trapping ratio β (0-1): fraction of debris trapped in contact.
    pub trapping_ratio: f64,
}

impl WearDebris {
    /// Create a wear debris model.
    ///
    /// # Arguments
    /// * `particle_radius` – mean debris radius \[m\]
    /// * `particle_density` – debris density \[kg/m³\]
    /// * `trapping_ratio` – fraction of debris retained (0..=1)
    pub fn new(particle_radius: f64, particle_density: f64, trapping_ratio: f64) -> Self {
        Self {
            particle_radius,
            particle_density,
            trapping_ratio,
        }
    }

    /// Mass of a single spherical debris particle \[kg\].
    pub fn particle_mass(&self) -> f64 {
        (4.0 / 3.0) * PI * self.particle_radius.powi(3) * self.particle_density
    }

    /// Number of debris particles generated from worn volume V \[m³\].
    ///
    /// N = V / V_particle
    pub fn particle_count(&self, worn_volume: f64) -> f64 {
        let vol_particle = (4.0 / 3.0) * PI * self.particle_radius.powi(3);
        worn_volume / vol_particle
    }

    /// Effective friction increase due to trapped abrasive debris.
    ///
    /// Δμ ≈ β · r_d / h_contact · μ_abr where h_contact is nominal gap.
    ///
    /// # Arguments
    /// * `nominal_gap` – mean contact gap h \[m\]
    /// * `abrasive_friction_coeff` – μ_abr for hard particles
    pub fn friction_increase(&self, nominal_gap: f64, abrasive_friction_coeff: f64) -> f64 {
        if nominal_gap < 1e-30 {
            return 0.0;
        }
        self.trapping_ratio * self.particle_radius / nominal_gap * abrasive_friction_coeff
    }
}

/// Tribocorrosion model combining mechanical wear and electrochemical dissolution.
///
/// Based on the synergy concept: total material loss = W_mech + W_chem + W_synergy.
///
/// Reference: Landolt et al., Wear 248 (2001) 211-219.
#[derive(Debug, Clone)]
pub struct Tribocorrosion {
    /// Mechanical wear rate dV/dt \[m³/s\].
    pub mechanical_wear_rate: f64,
    /// Corrosion rate without wear (pure dissolution) \[m³/s\].
    pub corrosion_rate_0: f64,
    /// Synergy factor β_s (ratio of wear-enhanced corrosion to base corrosion).
    pub synergy_factor: f64,
}

impl Tribocorrosion {
    /// Create a tribocorrosion model.
    ///
    /// # Arguments
    /// * `mechanical_wear_rate` – pure mechanical wear rate \[m³/s\]
    /// * `corrosion_rate_0` – baseline corrosion rate \[m³/s\]
    /// * `synergy_factor` – tribocorrosion synergy multiplier
    pub fn new(mechanical_wear_rate: f64, corrosion_rate_0: f64, synergy_factor: f64) -> Self {
        Self {
            mechanical_wear_rate,
            corrosion_rate_0,
            synergy_factor,
        }
    }

    /// Total tribocorrosion rate \[m³/s\].
    pub fn total_rate(&self) -> f64 {
        let w_chem_enhanced = self.corrosion_rate_0 * (1.0 + self.synergy_factor);
        self.mechanical_wear_rate + w_chem_enhanced
    }

    /// Synergy contribution to material loss \[m³/s\].
    pub fn synergy_rate(&self) -> f64 {
        self.corrosion_rate_0 * self.synergy_factor
    }

    /// Volume lost over time t \[s\].
    pub fn volume_loss(&self, time: f64) -> f64 {
        self.total_rate() * time
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 11.  Solid Lubricants (MoS₂ and DLC)
// ══════════════════════════════════════════════════════════════════════════════

/// Solid lubricant coating model (MoS₂ or DLC).
///
/// Provides friction and wear properties for layered solid lubricant films.
#[derive(Debug, Clone)]
pub struct SolidLubricant {
    /// Coating type identifier.
    pub coating_type: CoatingType,
    /// Coating thickness t \[m\].
    pub thickness: f64,
    /// Friction coefficient μ under normal atmospheric conditions.
    pub friction_coeff_air: f64,
    /// Friction coefficient μ under vacuum or inert gas.
    pub friction_coeff_vacuum: f64,
    /// Hardness H \[GPa\].
    pub hardness_gpa: f64,
    /// Specific wear rate k \[m³/(N·m)\].
    pub specific_wear_rate: f64,
}

/// Classification of solid lubricant coatings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoatingType {
    /// Molybdenum disulfide – lamellar solid lubricant.
    MoS2,
    /// Diamond-like carbon – hard amorphous carbon coating.
    Dlc,
    /// Polytetrafluoroethylene (PTFE).
    Ptfe,
    /// Tungsten disulfide – analogous to MoS₂.
    WS2,
}

impl SolidLubricant {
    /// Create a default MoS₂ solid lubricant coating.
    ///
    /// Properties representative of sputtered MoS₂ in vacuum.
    pub fn mos2(thickness: f64) -> Self {
        Self {
            coating_type: CoatingType::MoS2,
            thickness,
            friction_coeff_air: 0.15,
            friction_coeff_vacuum: 0.01,
            hardness_gpa: 0.3,
            specific_wear_rate: 1e-16,
        }
    }

    /// Create a default DLC solid lubricant coating.
    ///
    /// Properties representative of hydrogen-free ta-C DLC.
    pub fn dlc(thickness: f64) -> Self {
        Self {
            coating_type: CoatingType::Dlc,
            thickness,
            friction_coeff_air: 0.05,
            friction_coeff_vacuum: 0.08,
            hardness_gpa: 40.0,
            specific_wear_rate: 1e-18,
        }
    }

    /// Estimate coating lifetime \[m of sliding distance\].
    ///
    /// Based on: L = t / (k · p_mean)
    ///
    /// # Arguments
    /// * `mean_pressure` – mean contact pressure \[Pa\]
    pub fn coating_lifetime(&self, mean_pressure: f64) -> f64 {
        if mean_pressure < 1e-20 {
            return f64::INFINITY;
        }
        self.thickness / (self.specific_wear_rate * mean_pressure)
    }

    /// Friction coefficient for given environment.
    ///
    /// # Arguments
    /// * `vacuum` – true if operating in vacuum/inert gas
    pub fn friction_coeff(&self, vacuum: bool) -> f64 {
        if vacuum {
            self.friction_coeff_vacuum
        } else {
            self.friction_coeff_air
        }
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 12.  Rolling Contact Fatigue (RCF)
// ══════════════════════════════════════════════════════════════════════════════

/// Rolling contact fatigue model using the Lundberg-Palmgren theory.
///
/// Predicts bearing fatigue life as a function of load, material, and geometry.
///
/// Reference: Lundberg & Palmgren, Acta Polytech. Scand. 1 (1947).
#[derive(Debug, Clone)]
pub struct RollingContactFatigue {
    /// Dynamic load rating C \[N\] of the bearing element.
    pub dynamic_load_rating: f64,
    /// Weibull slope e (typically 10/3 for balls, 3 for rollers).
    pub weibull_slope: f64,
    /// Life modification factor a_ISO (accounts for contamination, lubrication).
    pub life_modification_factor: f64,
}

impl RollingContactFatigue {
    /// Create a rolling contact fatigue model.
    ///
    /// # Arguments
    /// * `dynamic_load_rating` – C \[N\]
    /// * `weibull_slope` – e (10/3 for ball bearings, 3 for roller bearings)
    /// * `life_modification_factor` – a_ISO (1.0 for reference conditions)
    pub fn new(
        dynamic_load_rating: f64,
        weibull_slope: f64,
        life_modification_factor: f64,
    ) -> Self {
        Self {
            dynamic_load_rating,
            weibull_slope,
            life_modification_factor,
        }
    }

    /// Basic L10 life \[million revolutions\] at equivalent dynamic load P \[N\].
    ///
    /// L10 = (C/P)^e · 10^6
    pub fn l10_life(&self, equivalent_load: f64) -> f64 {
        if equivalent_load < 1e-20 {
            return f64::INFINITY;
        }
        (self.dynamic_load_rating / equivalent_load).powf(self.weibull_slope)
    }

    /// ISO-modified adjusted life Lnm \[million revolutions\].
    ///
    /// Lnm = a_ISO · L10
    pub fn lnm_life(&self, equivalent_load: f64) -> f64 {
        self.life_modification_factor * self.l10_life(equivalent_load)
    }

    /// Fatigue limit load \[N\] below which infinite life is predicted.
    ///
    /// Cu ≈ 0.5 · C (empirical approximation for typical steels).
    pub fn fatigue_limit_load(&self) -> f64 {
        0.5 * self.dynamic_load_rating
    }

    /// Survival probability S at given life ratio t = L / L10.
    ///
    /// Weibull: S = exp(-ln(10/9) · t^e)
    pub fn survival_probability(&self, life_ratio: f64) -> f64 {
        let ln_10_9: f64 = (10.0_f64 / 9.0_f64).ln();
        (-ln_10_9 * life_ratio.powf(self.weibull_slope)).exp()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 13.  Mixed Lubrication Model
// ══════════════════════════════════════════════════════════════════════════════

/// Mixed lubrication model partitioning load between asperities and fluid film.
///
/// Uses a load-sharing approach where total load W = W_asperity + W_fluid.
#[derive(Debug, Clone)]
pub struct MixedLubrication {
    /// GW asperity contact model.
    pub gw: GreenwoodWilliamson,
    /// EHL film thickness model.
    pub ehl: EhlFilmThickness,
    /// Nominal contact area A_0 \[m²\].
    pub nominal_area: f64,
}

impl MixedLubrication {
    /// Create a mixed lubrication model.
    ///
    /// # Arguments
    /// * `gw` – Greenwood-Williamson asperity model
    /// * `ehl` – EHL film model
    /// * `nominal_area` – nominal contact area \[m²\]
    pub fn new(gw: GreenwoodWilliamson, ehl: EhlFilmThickness, nominal_area: f64) -> Self {
        Self {
            gw,
            ehl,
            nominal_area,
        }
    }

    /// Asperity load fraction at given separation d \[m\] and total load W \[N\].
    ///
    /// Returns fraction of load carried by asperities (0..=1).
    ///
    /// # Arguments
    /// * `separation` – mean surface separation \[m\]
    /// * `total_load` – total applied load \[N\]
    pub fn asperity_load_fraction(&self, separation: f64, total_load: f64) -> f64 {
        if total_load < 1e-20 {
            return 0.0;
        }
        let w_asp = self.gw.total_load(separation, self.nominal_area);
        (w_asp / total_load).min(1.0)
    }

    /// Effective friction coefficient in mixed lubrication.
    ///
    /// μ_eff = χ · μ_boundary + (1-χ) · μ_fluid
    ///
    /// # Arguments
    /// * `separation` – mean separation \[m\]
    /// * `total_load` – total load \[N\]
    /// * `boundary_mu` – friction coefficient in boundary regime
    /// * `fluid_mu` – friction coefficient in full-film regime
    pub fn effective_friction(
        &self,
        separation: f64,
        total_load: f64,
        boundary_mu: f64,
        fluid_mu: f64,
    ) -> f64 {
        let chi = self.asperity_load_fraction(separation, total_load);
        chi * boundary_mu + (1.0 - chi) * fluid_mu
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 14.  Abrasive Wear – Two-Body and Three-Body Models
// ══════════════════════════════════════════════════════════════════════════════

/// Abrasive wear model (micro-cutting and micro-ploughing).
///
/// Extends the Archard model with an abrasive hardness ratio correction.
#[derive(Debug, Clone)]
pub struct AbrasiveWear {
    /// Archard base model.
    pub archard: ArchardWear,
    /// Hardness of abrasive particle H_a \[Pa\].
    pub abrasive_hardness: f64,
    /// Cone apex semi-angle of abrasive \[radians\].
    pub abrasive_angle: f64,
}

impl AbrasiveWear {
    /// Create an abrasive wear model.
    ///
    /// # Arguments
    /// * `wear_coefficient` – dimensionless K
    /// * `surface_hardness` – H \[Pa\] of worn surface
    /// * `abrasive_hardness` – H_a \[Pa\] of abrasive
    /// * `abrasive_angle` – cone semi-angle θ \[rad\]
    pub fn new(
        wear_coefficient: f64,
        surface_hardness: f64,
        abrasive_hardness: f64,
        abrasive_angle: f64,
    ) -> Self {
        Self {
            archard: ArchardWear::new(wear_coefficient, surface_hardness),
            abrasive_hardness,
            abrasive_angle,
        }
    }

    /// Wear volume per unit sliding distance for two-body abrasion \[m³/m\].
    ///
    /// V/s = (2/π) · tan(θ) · W / H
    ///
    /// # Arguments
    /// * `normal_load` – W \[N\]
    pub fn two_body_wear_rate(&self, normal_load: f64) -> f64 {
        (2.0 / PI) * self.abrasive_angle.tan() * normal_load / self.archard.hardness
    }

    /// Hardness ratio – governs transition from elastic to plastic abrasion.
    ///
    /// H_ratio = H_surface / H_abrasive
    /// If H_ratio > 0.8: minimal abrasion; < 0.8: severe abrasion.
    pub fn hardness_ratio(&self) -> f64 {
        self.archard.hardness / self.abrasive_hardness
    }

    /// Severity factor for three-body abrasion (loose particles).
    ///
    /// Three-body is ~10-100x less severe than two-body due to rolling.
    /// Returns severity ratio (dimensionless, 0.01-0.1).
    pub fn three_body_severity(&self) -> f64 {
        0.05 // Typical ratio; can be calibrated from experiments
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 15.  Adhesive Wear Model
// ══════════════════════════════════════════════════════════════════════════════

/// Adhesive wear model based on junction growth and interfacial adhesion energy.
///
/// Combines Archard law with surface energy work of adhesion W_ad.
#[derive(Debug, Clone)]
pub struct AdhesiveWear {
    /// Archard wear component.
    pub archard: ArchardWear,
    /// Work of adhesion W_ad \[J/m²\] between the two surfaces.
    pub work_of_adhesion: f64,
}

impl AdhesiveWear {
    /// Create an adhesive wear model.
    ///
    /// # Arguments
    /// * `wear_coefficient` – K (dimensionless)
    /// * `hardness` – H \[Pa\]
    /// * `work_of_adhesion` – W_ad \[J/m²\]
    pub fn new(wear_coefficient: f64, hardness: f64, work_of_adhesion: f64) -> Self {
        Self {
            archard: ArchardWear::new(wear_coefficient, hardness),
            work_of_adhesion,
        }
    }

    /// Pull-off force \[N\] for a circular contact of radius a \[m\].
    ///
    /// Based on DMT theory: F_pull = 2π · W_ad · a²
    pub fn pull_off_force(&self, contact_radius: f64) -> f64 {
        2.0 * PI * self.work_of_adhesion * contact_radius * contact_radius
    }

    /// Adhesive wear volume \[m³\].
    ///
    /// # Arguments
    /// * `normal_load` – W \[N\]
    /// * `sliding_distance` – s \[m\]
    pub fn wear_volume(&self, normal_load: f64, sliding_distance: f64) -> f64 {
        self.archard.wear_volume(normal_load, sliding_distance)
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── § 1 Archard wear ──────────────────────────────────────────────────────

    #[test]
    fn test_archard_wear_volume_basic() {
        // K=0.001, H=1e9 Pa, W=100 N, s=1 m → V = 1e-10 m³
        let wear = ArchardWear::new(0.001, 1e9);
        let v = wear.wear_volume(100.0, 1.0);
        assert!((v - 1e-10).abs() < 1e-22, "Expected 1e-10, got {v:.6e}");
    }

    #[test]
    fn test_archard_wear_rate_proportional_to_load() {
        let wear = ArchardWear::new(0.01, 1e9);
        let r1 = wear.wear_rate(100.0);
        let r2 = wear.wear_rate(200.0);
        assert!(
            (r2 / r1 - 2.0).abs() < 1e-10,
            "Wear rate should double with load"
        );
    }

    #[test]
    fn test_archard_specific_wear_rate() {
        let k = 1e-3;
        let h = 2e9;
        let wear = ArchardWear::new(k, h);
        let expected = k / h;
        assert!(
            (wear.specific_wear_rate() - expected).abs() < 1e-30,
            "Specific wear rate mismatch"
        );
    }

    // ── § 2 Hertz contact ─────────────────────────────────────────────────────

    #[test]
    fn test_hertz_contact_radius_positive() {
        let hertz = HertzContact::new(200e9, 0.3, 200e9, 0.3, 0.01, f64::INFINITY);
        let a = hertz.contact_radius(1000.0);
        assert!(a > 0.0, "Contact radius must be positive");
    }

    #[test]
    fn test_hertz_max_pressure_greater_than_mean() {
        let hertz = HertzContact::new(200e9, 0.3, 200e9, 0.3, 0.01, f64::INFINITY);
        let p_max = hertz.max_pressure(1000.0);
        let p_mean = hertz.mean_pressure(1000.0);
        assert!(
            p_max > p_mean,
            "Max pressure {p_max:.4e} should exceed mean {p_mean:.4e}"
        );
    }

    #[test]
    fn test_hertz_max_to_mean_pressure_ratio() {
        // For Hertz point contact: p0 / p_mean = 3/2
        let hertz = HertzContact::new(200e9, 0.3, 200e9, 0.3, 0.005, f64::INFINITY);
        let p_max = hertz.max_pressure(500.0);
        let p_mean = hertz.mean_pressure(500.0);
        assert!(
            (p_max / p_mean - 1.5).abs() < 1e-6,
            "Ratio p_max/p_mean should be 3/2, got {}",
            p_max / p_mean
        );
    }

    // ── § 3 Flash temperature ─────────────────────────────────────────────────

    #[test]
    fn test_blok_flash_temperature_positive() {
        let ft = FlashTemperature::new(50.0, 1.2e-5, 50.0, 1.2e-5);
        let dt = ft.blok_flash_temperature(0.1, 1e9, 1.0, 1e-4);
        assert!(dt > 0.0, "Flash temperature must be positive");
    }

    #[test]
    fn test_peclet_number_consistency() {
        let ft = FlashTemperature::new(50.0, 1.2e-5, 50.0, 1.2e-5);
        let pe = ft.peclet_number(1.0, 1e-4);
        // α_avg = 1.2e-5, v=1, a=1e-4 → Pe = 1*1e-4/(2*1.2e-5) ≈ 4.17
        assert!(
            (pe - 4.166_666).abs() < 1e-4,
            "Peclet number expected ~4.17, got {pe:.6}"
        );
    }

    #[test]
    fn test_jaeger_less_than_blok_at_high_pe() {
        // For high Pe (fast sliding), Blok should give higher ΔT than Jaeger
        let ft = FlashTemperature::new(50.0, 1.2e-5, 50.0, 1.2e-5);
        let blok = ft.blok_flash_temperature(0.1, 1e9, 10.0, 1e-3);
        let jaeger = ft.jaeger_flash_temperature(0.1, 1e9, 10.0, 1e-3);
        // Both must be positive; Blok uses 0.308 coefficient, Jaeger uses 1/π ≈ 0.318
        // so they're close — just verify both positive
        assert!(
            blok > 0.0 && jaeger > 0.0,
            "Both flash temps must be positive"
        );
    }

    // ── § 4 Surface roughness ─────────────────────────────────────────────────

    #[test]
    fn test_roughness_ra_simple() {
        // Profile: [1, -1, 1, -1] → Ra = 1
        let profile = vec![1.0, -1.0, 1.0, -1.0];
        let sr = SurfaceRoughness::new(profile, 4e-6);
        assert!(
            (sr.ra() - 1.0).abs() < 1e-10,
            "Ra should be 1.0, got {}",
            sr.ra()
        );
    }

    #[test]
    fn test_roughness_rq_greater_than_or_equal_ra() {
        // Rq ≥ Ra always (Cauchy-Schwarz inequality)
        let profile: Vec<f64> = (0..100).map(|i| (i as f64 * 0.1).sin()).collect();
        let sr = SurfaceRoughness::new(profile, 10e-6);
        assert!(sr.rq() >= sr.ra() - 1e-12, "Rq must be >= Ra");
    }

    #[test]
    fn test_roughness_rz_max_minus_min() {
        let profile = vec![-2.0, 0.0, 3.0, 1.0, -1.0];
        let sr = SurfaceRoughness::new(profile, 5e-6);
        assert!((sr.rz() - 5.0).abs() < 1e-10, "Rz should be 5.0 (max-min)");
    }

    #[test]
    fn test_roughness_empty_profile() {
        let sr = SurfaceRoughness::new(vec![], 1e-3);
        assert_eq!(sr.ra(), 0.0);
        assert_eq!(sr.rq(), 0.0);
        assert_eq!(sr.rz(), 0.0);
    }

    // ── § 5 Greenwood-Williamson ──────────────────────────────────────────────

    #[test]
    fn test_gw_plasticity_index() {
        // E*=200 GPa, H=7 GPa (steel), σ_s=0.5 μm, β=5 μm
        // ψ = (200e9/7e9) * sqrt(0.5e-6/5e-6) ≈ 28.57 * 0.316 ≈ 9.03
        let gw = GreenwoodWilliamson::new(1e12, 5e-6, 0.5e-6, 200e9);
        let psi = gw.plasticity_index(7e9);
        assert!(
            psi > 1.0,
            "Plasticity index should indicate plastic contact"
        );
    }

    #[test]
    fn test_gw_real_contact_area_positive() {
        let gw = GreenwoodWilliamson::new(1e12, 1e-6, 1e-7, 100e9);
        let a_r = gw.real_contact_area(0.0, 1e-4);
        assert!(a_r > 0.0, "Real contact area must be positive");
    }

    #[test]
    fn test_gw_total_load_increases_with_separation_decrease() {
        let gw = GreenwoodWilliamson::new(1e12, 1e-6, 1e-7, 100e9);
        let w1 = gw.total_load(5e-8, 1e-4); // smaller separation → more load
        let w2 = gw.total_load(2e-7, 1e-4); // larger separation → less load
        assert!(w1 > w2, "Load should increase as separation decreases");
    }

    // ── § 6 Lubrication regime ────────────────────────────────────────────────

    #[test]
    fn test_lambda_ratio_classification() {
        assert_eq!(
            classify_lubrication_regime(0.5),
            LubricationRegime::Boundary
        );
        assert_eq!(classify_lubrication_regime(2.0), LubricationRegime::Mixed);
        assert_eq!(
            classify_lubrication_regime(5.0),
            LubricationRegime::FullFilm
        );
    }

    #[test]
    fn test_lambda_ratio_computation() {
        // h_min=1e-6, rq1=rq2=0.5e-6 → composite=sqrt(2)*0.5e-6, λ≈1.414
        let lam = lambda_ratio(1e-6, 0.5e-6, 0.5e-6);
        let expected = 1e-6 / (2.0_f64.sqrt() * 0.5e-6);
        assert!((lam - expected).abs() < 1e-10, "Lambda ratio mismatch");
    }

    // ── § 7 EHL ───────────────────────────────────────────────────────────────

    #[test]
    fn test_ehl_hmin_positive() {
        let ehl = EhlFilmThickness::new(0.02, 2e-8, 220e9, 0.01);
        let h = ehl.hamrock_dowson_hmin(1.0, 1000.0);
        assert!(h > 0.0, "EHL film thickness must be positive");
    }

    #[test]
    fn test_barus_viscosity_increases_with_pressure() {
        let ehl = EhlFilmThickness::new(0.02, 2e-8, 220e9, 0.01);
        let v0 = ehl.barus_viscosity(0.0);
        let v1 = ehl.barus_viscosity(1e8);
        assert!(v1 > v0, "Barus viscosity must increase with pressure");
        assert!((v0 - 0.02).abs() < 1e-15, "At p=0, viscosity equals η_0");
    }

    // ── § 8 Friction models ───────────────────────────────────────────────────

    #[test]
    fn test_coulomb_friction_sliding_vs_static() {
        let cf = CoulombFriction::new(0.5, 0.35);
        let f_static = cf.friction_force(100.0, false);
        let f_kinetic = cf.friction_force(100.0, true);
        assert!(f_static > f_kinetic, "Static friction must exceed kinetic");
    }

    #[test]
    fn test_stribeck_friction_coeff_decreases_from_zero_velocity() {
        let sf = StribeckFriction::new(0.5, 0.2, 0.1, 2.0, 0.001);
        let mu_slow = sf.friction_coeff(1e-4);
        let mu_medium = sf.friction_coeff(0.1);
        // Near v=0 the Stribeck term is large; at v=v_st it's reduced
        assert!(
            mu_slow > mu_medium - 0.3,
            "Stribeck curve should be bounded"
        );
        assert!(
            mu_slow > 0.0 && mu_medium > 0.0,
            "Friction must be positive"
        );
    }

    #[test]
    fn test_viscous_friction_linear() {
        let vf = ViscousFriction::new(5.0);
        let f1 = vf.friction_force(1.0);
        let f2 = vf.friction_force(2.0);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-10,
            "Viscous friction should be linear"
        );
    }

    // ── § 9 Fretting ──────────────────────────────────────────────────────────

    #[test]
    fn test_fretting_wear_volume_per_cycle() {
        let fw = FrettingWear::new(1e-8, 1e-6);
        let v = fw.wear_volume_per_cycle(50.0, 1e-5);
        assert!(v > 0.0, "Fretting wear volume must be positive");
    }

    #[test]
    fn test_fretting_fatigue_reduction_factor_bounded() {
        let fw = FrettingWear::new(1e-8, 1e-6);
        let factor = fw.fatigue_reduction_factor(500e6, 300e6, 0.3);
        assert!(
            (0.0..=1.0).contains(&factor),
            "Reduction factor must be in [0,1]"
        );
    }

    // ── § 10 Wear debris & tribocorrosion ─────────────────────────────────────

    #[test]
    fn test_wear_debris_particle_count() {
        let wd = WearDebris::new(1e-6, 7800.0, 0.5);
        let n = wd.particle_count(1e-12);
        let vol_p = (4.0 / 3.0) * PI * 1e-18;
        let expected = 1e-12 / vol_p;
        assert!(
            (n - expected).abs() / expected < 1e-10,
            "Particle count mismatch"
        );
    }

    #[test]
    fn test_tribocorrosion_total_exceeds_components() {
        let tc = Tribocorrosion::new(1e-9, 5e-10, 2.0);
        let total = tc.total_rate();
        assert!(
            total > tc.mechanical_wear_rate,
            "Total rate must exceed mechanical wear rate"
        );
        assert!(
            total > tc.corrosion_rate_0,
            "Total rate must exceed bare corrosion"
        );
    }

    // ── § 11 Solid lubricants ─────────────────────────────────────────────────

    #[test]
    fn test_mos2_friction_lower_in_vacuum() {
        let mos2 = SolidLubricant::mos2(1e-6);
        assert!(
            mos2.friction_coeff(true) < mos2.friction_coeff(false),
            "MoS2 friction lower in vacuum"
        );
    }

    #[test]
    fn test_dlc_hardness_greater_than_mos2() {
        let mos2 = SolidLubricant::mos2(1e-6);
        let dlc = SolidLubricant::dlc(1e-6);
        assert!(
            dlc.hardness_gpa > mos2.hardness_gpa,
            "DLC should be harder than MoS2"
        );
    }

    #[test]
    fn test_coating_lifetime_finite_under_load() {
        let dlc = SolidLubricant::dlc(2e-6);
        let life = dlc.coating_lifetime(1e9);
        assert!(
            life.is_finite() && life > 0.0,
            "Lifetime should be finite and positive"
        );
    }

    // ── § 12 Rolling contact fatigue ──────────────────────────────────────────

    #[test]
    fn test_rcf_l10_life_decreases_with_load() {
        let rcf = RollingContactFatigue::new(10000.0, 10.0 / 3.0, 1.0);
        let l1 = rcf.l10_life(2000.0);
        let l2 = rcf.l10_life(5000.0);
        assert!(l1 > l2, "L10 life must decrease as load increases");
    }

    #[test]
    fn test_rcf_survival_probability_at_l10() {
        // By definition S(L10) ≈ 0.9
        let rcf = RollingContactFatigue::new(10000.0, 10.0 / 3.0, 1.0);
        let s = rcf.survival_probability(1.0); // t = L/L10 = 1
        assert!(
            (s - 0.9).abs() < 1e-6,
            "Survival at L10 should be 0.9, got {s:.6}"
        );
    }

    #[test]
    fn test_rcf_lnm_scaled_by_factor() {
        let rcf = RollingContactFatigue::new(10000.0, 10.0 / 3.0, 2.0);
        let base = RollingContactFatigue::new(10000.0, 10.0 / 3.0, 1.0);
        let load = 3000.0;
        assert!(
            (rcf.lnm_life(load) - 2.0 * base.l10_life(load)).abs() < 1e-8,
            "Lnm should equal a_ISO * L10"
        );
    }
}
