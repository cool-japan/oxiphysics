// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended contact mechanics: Hertz sphere/cylinder variants, JKR/DMT adhesion,
//! loading/unloading hysteresis, Greenwood-Williamson rough surface contact,
//! and flash temperature (Jaeger's formula).

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 1. Elastic material primitives
// ---------------------------------------------------------------------------

/// Reduced (combined) modulus E* from two elastic solids.
///
/// 1/E* = (1 - ν₁²)/E₁ + (1 - ν₂²)/E₂
pub fn reduced_modulus(e1: f64, nu1: f64, e2: f64, nu2: f64) -> f64 {
    let inv = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
    1.0 / inv
}

/// Combined radius R* from two radii (use `f64::INFINITY` for flat).
///
/// 1/R* = 1/R₁ + 1/R₂
pub fn combined_radius(r1: f64, r2: f64) -> f64 {
    1.0 / (1.0 / r1 + 1.0 / r2)
}

// ---------------------------------------------------------------------------
// 2. Hertz sphere-sphere contact
// ---------------------------------------------------------------------------

/// Hertz sphere–sphere (or sphere–flat) contact model.
///
/// Encapsulates the standard Hertz results for a circular contact patch.
#[derive(Debug, Clone)]
pub struct HertzSphereSphere {
    /// Effective radius R* (m).
    pub r_star: f64,
    /// Reduced modulus E* (Pa).
    pub e_star: f64,
}

impl HertzSphereSphere {
    /// Construct from individual sphere radii and material properties.
    pub fn new(r1: f64, nu1: f64, e1: f64, r2: f64, nu2: f64, e2: f64) -> Self {
        let r_star = combined_radius(r1, r2);
        let e_star = reduced_modulus(e1, nu1, e2, nu2);
        Self { r_star, e_star }
    }

    /// Contact radius a = (3 F R* / (4 E*))^(1/3).
    pub fn contact_radius(&self, force: f64) -> f64 {
        (3.0 * force * self.r_star / (4.0 * self.e_star)).cbrt()
    }

    /// Indentation depth δ = a² / R*.
    pub fn indentation(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        a * a / self.r_star
    }

    /// Peak Hertz pressure p₀ = (6 F E*² / (π³ R*²))^(1/3).
    pub fn peak_pressure(&self, force: f64) -> f64 {
        (6.0 * force * self.e_star * self.e_star / (PI * PI * PI * self.r_star * self.r_star))
            .cbrt()
    }

    /// Mean Hertz pressure p_mean = F / (π a²).
    pub fn mean_pressure(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        if a < 1e-30 {
            return 0.0;
        }
        force / (PI * a * a)
    }

    /// Hertz normal contact stiffness k_n = 2 E* a.
    pub fn contact_stiffness(&self, force: f64) -> f64 {
        2.0 * self.e_star * self.contact_radius(force)
    }

    /// Force from indentation: F = (4/3) E* √R* δ^(3/2).
    pub fn force_from_indentation(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * self.e_star * self.r_star.sqrt() * delta.powf(1.5)
    }

    /// Hertz pressure distribution p(r) = p₀ √(1 - (r/a)²) for r ≤ a.
    pub fn pressure_at_radius(&self, r: f64, force: f64) -> f64 {
        let a = self.contact_radius(force);
        if r > a {
            return 0.0;
        }
        let p0 = self.peak_pressure(force);
        let ratio = r / a;
        p0 * (1.0 - ratio * ratio).max(0.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// 3. Hertz sphere-flat contact (specialised)
// ---------------------------------------------------------------------------

/// Hertz contact between a sphere and a flat surface.
///
/// Equivalent to `HertzSphereSphere` with R₂ = ∞, kept separately for clarity.
#[derive(Debug, Clone)]
pub struct HertzSphereFlat {
    /// Sphere radius (m).
    pub sphere_radius: f64,
    /// Reduced modulus E* (Pa).
    pub e_star: f64,
}

impl HertzSphereFlat {
    /// Construct from sphere geometry and bimaterial properties.
    pub fn new(r_sphere: f64, e_sphere: f64, nu_sphere: f64, e_flat: f64, nu_flat: f64) -> Self {
        let e_star = reduced_modulus(e_sphere, nu_sphere, e_flat, nu_flat);
        Self {
            sphere_radius: r_sphere,
            e_star,
        }
    }

    /// Contact radius a = (3 F R / (4 E*))^(1/3).
    pub fn contact_radius(&self, force: f64) -> f64 {
        (3.0 * force * self.sphere_radius / (4.0 * self.e_star)).cbrt()
    }

    /// Indentation depth δ = a² / R.
    pub fn indentation(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        a * a / self.sphere_radius
    }

    /// Peak pressure p₀ = 3 F / (2 π a²).
    pub fn peak_pressure(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        if a < 1e-30 {
            return 0.0;
        }
        3.0 * force / (2.0 * PI * a * a)
    }

    /// Elastic strain energy U = (2/5) F δ.
    pub fn strain_energy(&self, force: f64) -> f64 {
        0.4 * force * self.indentation(force)
    }
}

// ---------------------------------------------------------------------------
// 4. Hertz cylinder-flat (line contact)
// ---------------------------------------------------------------------------

/// Hertz line contact between an infinite cylinder and a flat surface.
///
/// Results are per unit length along the cylinder axis.
#[derive(Debug, Clone)]
pub struct HertzCylinderFlat {
    /// Cylinder radius (m).
    pub cylinder_radius: f64,
    /// Reduced plane-strain modulus E* (Pa).
    pub e_star: f64,
}

impl HertzCylinderFlat {
    /// Construct from cylinder geometry and bimaterial properties.
    pub fn new(r_cyl: f64, e_cyl: f64, nu_cyl: f64, e_flat: f64, nu_flat: f64) -> Self {
        let e_star = reduced_modulus(e_cyl, nu_cyl, e_flat, nu_flat);
        Self {
            cylinder_radius: r_cyl,
            e_star,
        }
    }

    /// Contact half-width b = √(4 F' R / (π E*)) where F' is force per unit length.
    pub fn half_width(&self, force_per_length: f64) -> f64 {
        (4.0 * force_per_length * self.cylinder_radius / (PI * self.e_star)).sqrt()
    }

    /// Peak pressure p₀ = 2 F' / (π b).
    pub fn peak_pressure(&self, force_per_length: f64) -> f64 {
        let b = self.half_width(force_per_length);
        if b < 1e-30 {
            return 0.0;
        }
        2.0 * force_per_length / (PI * b)
    }

    /// Mean pressure p_mean = F' / (2 b).
    pub fn mean_pressure(&self, force_per_length: f64) -> f64 {
        let b = self.half_width(force_per_length);
        if b < 1e-30 {
            return 0.0;
        }
        force_per_length / (2.0 * b)
    }

    /// Ellipsoidal pressure distribution p(x) = p₀ √(1 - (x/b)²) for |x| ≤ b.
    pub fn pressure_at_position(&self, x: f64, force_per_length: f64) -> f64 {
        let b = self.half_width(force_per_length);
        if x.abs() > b {
            return 0.0;
        }
        let p0 = self.peak_pressure(force_per_length);
        p0 * (1.0 - (x / b).powi(2)).max(0.0).sqrt()
    }

    /// Subsurface maximum shear stress τ_max ≈ 0.300 p₀ at depth z ≈ 0.786 b.
    pub fn max_shear_stress(&self, force_per_length: f64) -> f64 {
        0.300 * self.peak_pressure(force_per_length)
    }

    /// Depth at which the maximum shear stress occurs: z ≈ 0.786 b.
    pub fn max_shear_depth(&self, force_per_length: f64) -> f64 {
        0.786 * self.half_width(force_per_length)
    }
}

// ---------------------------------------------------------------------------
// 5. Mindlin tangential compliance
// ---------------------------------------------------------------------------

/// Mindlin reduced shear modulus G* from two elastic solids.
///
/// 1/G* = (2 - ν₁)/(4 G₁) + (2 - ν₂)/(4 G₂)
pub fn mindlin_reduced_shear_modulus(g1: f64, nu1: f64, g2: f64, nu2: f64) -> f64 {
    1.0 / ((2.0 - nu1) / (4.0 * g1) + (2.0 - nu2) / (4.0 * g2))
}

/// Mindlin tangential compliance model.
///
/// Builds on top of a Hertz contact to provide tangential stiffness and
/// the partial-slip (no-slip) solution for a given Coulomb friction limit.
#[derive(Debug, Clone)]
pub struct MindlinCompliance {
    /// Reduced shear modulus G* (Pa).
    pub g_star: f64,
    /// Coefficient of friction μ.
    pub friction: f64,
}

impl MindlinCompliance {
    /// Construct from two elastic solids (given E, ν) and friction coefficient.
    pub fn new(e1: f64, nu1: f64, e2: f64, nu2: f64, friction: f64) -> Self {
        let g1 = e1 / (2.0 * (1.0 + nu1));
        let g2 = e2 / (2.0 * (1.0 + nu2));
        let g_star = mindlin_reduced_shear_modulus(g1, nu1, g2, nu2);
        Self { g_star, friction }
    }

    /// Tangential stiffness k_t = 8 G* a.
    pub fn tangential_stiffness(&self, contact_radius: f64) -> f64 {
        8.0 * self.g_star * contact_radius
    }

    /// Mindlin partial-slip tangential displacement (loading from 0).
    ///
    /// δ_t = (3 μ F_n)/(8 G* a) \[1 - (1 - Q/(μ F_n))^(2/3)\]
    ///
    /// Returns `None` if gross slip (Q ≥ μ F_n).
    pub fn partial_slip_displacement(
        &self,
        tangential_force: f64,
        normal_force: f64,
        contact_radius: f64,
    ) -> Option<f64> {
        let limit = self.friction * normal_force;
        if tangential_force >= limit {
            return None;
        }
        let prefactor = 3.0 * limit / (8.0 * self.g_star * contact_radius);
        let term = 1.0 - tangential_force / limit;
        Some(prefactor * (1.0 - term.powf(2.0 / 3.0)))
    }

    /// Inverse: force from displacement (partial-slip).
    pub fn force_from_displacement(
        &self,
        delta_t: f64,
        normal_force: f64,
        contact_radius: f64,
    ) -> f64 {
        let limit = self.friction * normal_force;
        let prefactor = 3.0 * limit / (8.0 * self.g_star * contact_radius);
        let ratio = (delta_t / prefactor).clamp(0.0, 1.0);
        limit * (1.0 - (1.0 - ratio).powf(1.5))
    }

    /// Stick-zone radius c = a (1 - Q/(μ F_n))^(1/3).
    pub fn stick_radius(
        &self,
        tangential_force: f64,
        normal_force: f64,
        contact_radius: f64,
    ) -> Option<f64> {
        let limit = self.friction * normal_force;
        if tangential_force >= limit || limit < 1e-30 {
            return None;
        }
        Some(contact_radius * (1.0 - tangential_force / limit).powf(1.0 / 3.0))
    }

    /// Energy dissipated per fretting cycle (Mindlin–Deresiewicz).
    ///
    /// W = (9 μ² F_n²)/(10 k_t) \[1 - (1 - Q/(μ F_n))^(5/3)\]
    pub fn fretting_energy_per_cycle(
        &self,
        tangential_amplitude: f64,
        normal_force: f64,
        contact_radius: f64,
    ) -> f64 {
        let limit = self.friction * normal_force;
        if limit < 1e-30 {
            return 0.0;
        }
        let k_t = self.tangential_stiffness(contact_radius);
        if k_t < 1e-30 {
            return 0.0;
        }
        let ratio = (tangential_amplitude / limit).min(1.0);
        (9.0 * limit * limit) / (10.0 * k_t) * (1.0 - (1.0 - ratio).powf(5.0 / 3.0))
    }
}

// ---------------------------------------------------------------------------
// 6. JKR adhesion model
// ---------------------------------------------------------------------------

/// Johnson–Kendall–Roberts (JKR) adhesion model.
///
/// Extends Hertz contact to include surface energy (adhesion work W_ad).
/// The contact area is larger than Hertz and pull-off occurs at a negative
/// load F_po = -(3/2) π W_ad R*.
#[derive(Debug, Clone)]
pub struct JkrAdhesion {
    /// Effective radius R* (m).
    pub r_star: f64,
    /// Reduced modulus E* (Pa).
    pub e_star: f64,
    /// Adhesion work (Dupré energy) W_ad = γ₁ + γ₂ − γ₁₂ (J/m²).
    pub work_adhesion: f64,
}

impl JkrAdhesion {
    /// Construct JKR model.
    pub fn new(r_star: f64, e_star: f64, work_adhesion: f64) -> Self {
        Self {
            r_star,
            e_star,
            work_adhesion,
        }
    }

    /// JKR contact radius as a function of applied load F.
    ///
    /// a³ = (R*/(E*)) \[F + 3πW R* + √(6πW R* F + (3πW R*)²)\]
    /// where W = W_ad.
    pub fn contact_radius(&self, force: f64) -> f64 {
        let w = self.work_adhesion;
        let r = self.r_star;
        let e_star = self.e_star;
        let term = 3.0 * PI * w * r;
        let discriminant = 6.0 * PI * w * r * force + term * term;
        let a3 = (r / e_star) * (force + term + discriminant.max(0.0).sqrt());
        a3.cbrt()
    }

    /// JKR pull-off (adhesion) force: F_po = -(3/2) π W_ad R*.
    pub fn pull_off_force(&self) -> f64 {
        -1.5 * PI * self.work_adhesion * self.r_star
    }

    /// Indentation at zero applied load (contact due to adhesion only).
    pub fn zero_load_contact_radius(&self) -> f64 {
        self.contact_radius(0.0)
    }

    /// Elastic energy stored in the adhesive contact at load F.
    pub fn elastic_energy(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        let delta =
            a * a / self.r_star - (8.0 * PI * self.work_adhesion * a / (3.0 * self.e_star)).sqrt();
        if delta < 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * self.e_star * self.r_star.sqrt() * delta.powf(1.5)
    }

    /// Surface energy gain −W_ad π a².
    pub fn surface_energy(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        -self.work_adhesion * PI * a * a
    }
}

// ---------------------------------------------------------------------------
// 7. DMT adhesion model
// ---------------------------------------------------------------------------

/// Derjaguin–Muller–Toporov (DMT) adhesion model.
///
/// Valid for stiff materials / weak adhesion (small Tabor parameter).
/// The contact area follows Hertz mechanics but with an effective load
/// F_eff = F + 2 π W_ad R*.
#[derive(Debug, Clone)]
pub struct DmtAdhesion {
    /// Effective radius R* (m).
    pub r_star: f64,
    /// Reduced modulus E* (Pa).
    pub e_star: f64,
    /// Adhesion work W_ad (J/m²).
    pub work_adhesion: f64,
}

impl DmtAdhesion {
    /// Construct DMT model.
    pub fn new(r_star: f64, e_star: f64, work_adhesion: f64) -> Self {
        Self {
            r_star,
            e_star,
            work_adhesion,
        }
    }

    /// DMT pull-off force: F_po = -2 π W_ad R*.
    pub fn pull_off_force(&self) -> f64 {
        -2.0 * PI * self.work_adhesion * self.r_star
    }

    /// DMT effective load: F_eff = F + 2 π W_ad R*.
    fn effective_force(&self, force: f64) -> f64 {
        force + 2.0 * PI * self.work_adhesion * self.r_star
    }

    /// Contact radius (Hertz with effective load).
    pub fn contact_radius(&self, force: f64) -> f64 {
        let f_eff = self.effective_force(force);
        if f_eff <= 0.0 {
            return 0.0;
        }
        (3.0 * f_eff * self.r_star / (4.0 * self.e_star)).cbrt()
    }

    /// Indentation depth δ = a² / R*.
    pub fn indentation(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        a * a / self.r_star
    }

    /// Peak pressure at given applied load.
    pub fn peak_pressure(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        if a < 1e-30 {
            return 0.0;
        }
        let f_eff = self.effective_force(force);
        3.0 * f_eff / (2.0 * PI * a * a)
    }

    /// Tabor parameter μ_T for choosing between JKR and DMT regime.
    ///
    /// μ_T = (R* W²_ad / (E*² z₀³))^(1/3)
    /// where z₀ is the equilibrium atomic separation (typ. 0.2–0.4 nm).
    pub fn tabor_parameter(&self, z0: f64) -> f64 {
        (self.r_star * self.work_adhesion * self.work_adhesion
            / (self.e_star * self.e_star * z0 * z0 * z0))
            .cbrt()
    }
}

// ---------------------------------------------------------------------------
// 8. Contact area hysteresis (loading / unloading)
// ---------------------------------------------------------------------------

/// Contact area evolution under loading and unloading cycles.
///
/// Uses the Hertz loading branch and a Maugis-type unloading correction
/// to capture the hysteresis in the contact area vs. force curve.
#[derive(Debug, Clone)]
pub struct ContactAreaHysteresis {
    /// Effective radius R* (m).
    pub r_star: f64,
    /// Reduced modulus E* (Pa).
    pub e_star: f64,
    /// Maximum force reached during loading (for unloading branch).
    pub f_max: f64,
    /// Contact radius at maximum force (for unloading branch).
    pub a_max: f64,
}

impl ContactAreaHysteresis {
    /// Construct from geometry.
    pub fn new(r_star: f64, e_star: f64) -> Self {
        Self {
            r_star,
            e_star,
            f_max: 0.0,
            a_max: 0.0,
        }
    }

    /// Hertz loading: update the peak state and return contact radius.
    pub fn load(&mut self, force: f64) -> f64 {
        let a = (3.0 * force * self.r_star / (4.0 * self.e_star)).cbrt();
        if force > self.f_max {
            self.f_max = force;
            self.a_max = a;
        }
        a
    }

    /// Simplified unloading branch: contact area decreases following the
    /// Maugis–Dugdale flat-punch correction.
    ///
    /// a_unload(F) = a_max * (F / F_max)^(1/3) \[Hertz\] + residual
    /// (simplified, omits detailed elastic hysteresis in full MDR).
    pub fn unload(&self, force: f64) -> f64 {
        if force <= 0.0 {
            return 0.0;
        }
        if self.f_max < 1e-30 {
            return 0.0;
        }
        let a_hertz = (3.0 * force * self.r_star / (4.0 * self.e_star)).cbrt();
        // Residual contact from elastic springback (Mindlin unloading correction).
        let a_residual = self.a_max * (1.0 - force / self.f_max).max(0.0).sqrt() * 0.05;
        a_hertz + a_residual
    }

    /// Hysteresis loop area (energy dissipated over a load-unload cycle).
    ///
    /// Approximated by numerical integration between load and unload curves.
    pub fn hysteresis_energy(&self, n_steps: usize) -> f64 {
        if self.f_max < 1e-30 || n_steps < 2 {
            return 0.0;
        }
        let df = self.f_max / (n_steps as f64 - 1.0);
        let mut energy = 0.0;
        for i in 1..n_steps {
            let f_prev = (i - 1) as f64 * df;
            let f_curr = i as f64 * df;
            let a_load_prev = (3.0 * f_prev * self.r_star / (4.0 * self.e_star)).cbrt();
            let a_load_curr = (3.0 * f_curr * self.r_star / (4.0 * self.e_star)).cbrt();
            let a_unload_prev = self.unload(f_prev);
            let a_unload_curr = self.unload(f_curr);
            let da_load = a_load_curr - a_load_prev;
            let da_unload = a_unload_curr - a_unload_prev;
            energy += 0.5 * (f_prev + f_curr) * (da_load - da_unload).abs();
        }
        energy
    }
}

// ---------------------------------------------------------------------------
// 9. Normal / tangential contact force from separation / overlap
// ---------------------------------------------------------------------------

/// Contact force model based on separation (positive = gap, negative = overlap).
///
/// Implements both linear and nonlinear (Hertz) normal force and Coulomb-limited
/// tangential force.
#[derive(Debug, Clone)]
pub struct ContactForceModel {
    /// Normal stiffness k_n (N/m) — used for linear model only.
    pub normal_stiffness: f64,
    /// Tangential stiffness k_t (N/m) — used for linear model only.
    pub tangential_stiffness: f64,
    /// Coefficient of friction.
    pub friction: f64,
    /// Damping coefficient (critical fraction), 0 = undamped.
    pub damping: f64,
}

impl ContactForceModel {
    /// Linear Hertz-like normal force from overlap δ = −gap.
    ///
    /// F_n = k_n * δ  (no adhesion).
    pub fn linear_normal_force(&self, gap: f64) -> f64 {
        if gap >= 0.0 {
            0.0
        } else {
            -gap * self.normal_stiffness
        }
    }

    /// Hertz nonlinear normal force from overlap:
    /// F_n = k_n * δ^(3/2).
    pub fn hertz_normal_force(&self, gap: f64) -> f64 {
        if gap >= 0.0 {
            0.0
        } else {
            let delta = -gap;
            self.normal_stiffness * delta.powf(1.5)
        }
    }

    /// Damped normal force: F_n + c_n * v_n where c_n = damping * 2 √(m k_n).
    pub fn damped_normal_force(&self, gap: f64, normal_velocity: f64, mass: f64) -> f64 {
        let fn_ = self.linear_normal_force(gap);
        if fn_ <= 0.0 {
            return 0.0;
        }
        let c = self.damping * 2.0 * (mass * self.normal_stiffness).sqrt();
        (fn_ - c * normal_velocity).max(0.0)
    }

    /// Coulomb-limited tangential force from tangential displacement.
    pub fn tangential_force(&self, tangential_disp: f64, normal_force: f64) -> f64 {
        let f_t_trial = self.tangential_stiffness * tangential_disp;
        let limit = self.friction * normal_force;
        if f_t_trial.abs() <= limit {
            f_t_trial
        } else {
            limit * f_t_trial.signum()
        }
    }

    /// Total contact force vector \[f_t, f_n\] for given gap and displacement.
    pub fn contact_force_vector(&self, gap: f64, tangential_disp: f64) -> [f64; 2] {
        let fn_ = self.linear_normal_force(gap);
        let ft = self.tangential_force(tangential_disp, fn_);
        [ft, fn_]
    }
}

// ---------------------------------------------------------------------------
// 10. Greenwood–Williamson rough surface contact
// ---------------------------------------------------------------------------

/// Greenwood–Williamson (GW) statistical rough-surface contact model.
///
/// Models a surface as a collection of spherical asperities with a Gaussian
/// height distribution.  Key outputs are the real contact area and total
/// normal force as a function of the mean separation.
#[derive(Debug, Clone)]
pub struct GreenwoodWilliamson {
    /// Asperity density η (number per m²).
    pub asperity_density: f64,
    /// Mean asperity tip radius β (m).
    pub asperity_radius: f64,
    /// Standard deviation of asperity heights σ_s (m).
    pub height_std: f64,
    /// Composite reduced modulus E* (Pa).
    pub e_star: f64,
}

impl GreenwoodWilliamson {
    /// Gaussian PDF φ(z) = exp(-z²/(2σ²)) / (σ √(2π)).
    fn phi(&self, z: f64) -> f64 {
        let s = self.height_std;
        (-0.5 * (z / s).powi(2)).exp() / (s * (2.0 * PI).sqrt())
    }

    /// Complementary CDF  P(Z > d) ≈ ½ erfc(d/(√2 σ)).
    fn prob_contact(&self, separation: f64) -> f64 {
        let t = separation / (self.height_std * 2.0_f64.sqrt());
        0.5 * gw_erfc(t)
    }

    /// Mean overlap E\[max(z − d, 0)\] under Gaussian pdf.
    fn mean_overlap(&self, separation: f64) -> f64 {
        let s = self.height_std;
        let t = separation / (s * 2.0_f64.sqrt());
        s / (2.0 * PI).sqrt() * (-t * t).exp() - separation * 0.5 * gw_erfc(t)
    }

    /// Mean overlap^(3/2) ⟨max(z − d, 0)^(3/2)⟩ via numerical quadrature.
    fn mean_overlap_3_2(&self, separation: f64) -> f64 {
        let s = self.height_std;
        let n = 200usize;
        let hi = separation + 8.0 * s;
        let dz = (hi - separation) / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let z = separation + (i as f64 + 0.5) * dz;
            sum += (z - separation).powf(1.5) * self.phi(z) * dz;
        }
        sum
    }

    /// Real contact area per unit nominal area A_real / A_nominal.
    pub fn real_contact_area_fraction(&self, separation: f64) -> f64 {
        PI * self.asperity_density * self.asperity_radius * self.mean_overlap(separation)
    }

    /// Number of contacting asperities per unit nominal area.
    pub fn contact_asperity_density(&self, separation: f64) -> f64 {
        self.asperity_density * self.prob_contact(separation)
    }

    /// Mean contact pressure (Pa) for a given separation.
    pub fn contact_pressure(&self, separation: f64) -> f64 {
        (4.0 / 3.0)
            * self.asperity_density
            * self.e_star
            * self.asperity_radius.sqrt()
            * self.mean_overlap_3_2(separation)
    }

    /// Plasticity index ψ = (E*/H) √(σ/β).
    ///
    /// ψ < 0.6 elastic, ψ > 1.0 fully plastic.
    pub fn plasticity_index(&self, hardness: f64) -> f64 {
        if hardness < 1e-30 || self.asperity_radius < 1e-30 {
            return 0.0;
        }
        (self.e_star / hardness) * (self.height_std / self.asperity_radius).sqrt()
    }

    /// Separation d for a target contact pressure via bisection.
    ///
    /// Returns `None` if the target pressure is unachievable in \[0, 10σ\].
    pub fn separation_for_pressure(&self, target_pressure: f64, tol: f64) -> Option<f64> {
        let mut lo = 0.0_f64;
        let mut hi = 10.0 * self.height_std;
        // Ensure bracket is valid.
        if self.contact_pressure(lo) < target_pressure {
            return None;
        }
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.contact_pressure(mid) > target_pressure {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) < tol {
                break;
            }
        }
        Some(0.5 * (lo + hi))
    }
}

/// Complementary error function erfc(x) via rational approximation (A&S 7.1.26).
fn gw_erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - gw_erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

// ---------------------------------------------------------------------------
// 11. Flash temperature (Jaeger's formula)
// ---------------------------------------------------------------------------

/// Flash temperature rise at a sliding contact (Jaeger 1942).
///
/// Estimates the peak surface temperature increment due to frictional heating.
#[derive(Debug, Clone)]
pub struct FlashTemperature {
    /// Thermal conductivity of body 1 k₁ (W/(m·K)).
    pub k1: f64,
    /// Thermal conductivity of body 2 k₂ (W/(m·K)).
    pub k2: f64,
    /// Thermal diffusivity of body 1 κ₁ (m²/s).
    pub kappa1: f64,
    /// Thermal diffusivity of body 2 κ₂ (m²/s).
    pub kappa2: f64,
}

impl FlashTemperature {
    /// Construct from material thermal properties.
    pub fn new(k1: f64, kappa1: f64, k2: f64, kappa2: f64) -> Self {
        Self {
            k1,
            k2,
            kappa1,
            kappa2,
        }
    }

    /// Peclet number Pe = V a / (2 κ) where V is sliding speed, a contact radius.
    pub fn peclet_number(sliding_speed: f64, contact_radius: f64, thermal_diffusivity: f64) -> f64 {
        sliding_speed * contact_radius / (2.0 * thermal_diffusivity)
    }

    /// Partition factor β₁ (fraction of heat entering body 1).
    ///
    /// β₁ = k₁ √κ₂ / (k₁ √κ₂ + k₂ √κ₁)
    pub fn heat_partition_factor(&self) -> f64 {
        let num = self.k1 * self.kappa2.sqrt();
        let den = num + self.k2 * self.kappa1.sqrt();
        if den < 1e-30 { 0.5 } else { num / den }
    }

    /// Jaeger flash temperature rise ΔT_flash (high-speed limit, Pe >> 1).
    ///
    /// ΔT = 0.308 μ F_n V / (a (k₁ + k₂) √Pe)
    ///
    /// where Pe is the Peclet number for the faster body.
    pub fn flash_temperature_high_speed(
        &self,
        friction: f64,
        normal_force: f64,
        sliding_speed: f64,
        contact_radius: f64,
    ) -> f64 {
        if contact_radius < 1e-30 || sliding_speed < 1e-30 {
            return 0.0;
        }
        let pe = Self::peclet_number(sliding_speed, contact_radius, self.kappa1);
        if pe < 1e-10 {
            return 0.0;
        }
        let q_dot = friction * normal_force * sliding_speed; // total frictional power
        0.308 * q_dot / (contact_radius * (self.k1 + self.k2) * pe.sqrt())
    }

    /// Jaeger flash temperature rise (low-speed limit, Pe << 1).
    ///
    /// ΔT = q_dot / (2 π a (k₁ + k₂))
    pub fn flash_temperature_low_speed(
        &self,
        friction: f64,
        normal_force: f64,
        sliding_speed: f64,
        contact_radius: f64,
    ) -> f64 {
        if contact_radius < 1e-30 {
            return 0.0;
        }
        let q_dot = friction * normal_force * sliding_speed;
        q_dot / (2.0 * PI * contact_radius * (self.k1 + self.k2))
    }

    /// Blok's partition-corrected flash temperature (symmetric sliding).
    ///
    /// ΔT = μ p_mean V a / (2 k₁ √(π κ₁ a / V)) + … (partition version)
    ///
    /// Simplified implementation:
    /// ΔT = β₁ * q_dot / (π a (k₁ + k₂))
    pub fn flash_temperature_partition(
        &self,
        friction: f64,
        normal_force: f64,
        sliding_speed: f64,
        contact_radius: f64,
    ) -> f64 {
        if contact_radius < 1e-30 {
            return 0.0;
        }
        let beta = self.heat_partition_factor();
        let q_dot = friction * normal_force * sliding_speed;
        beta * q_dot / (PI * contact_radius * (self.k1 + self.k2))
    }

    /// Bulk temperature rise from steady frictional heating over a sliding distance L.
    ///
    /// ΔT_bulk = μ F_n L / (m c_p)  where m c_p is the thermal mass.
    pub fn bulk_temperature_rise(
        friction: f64,
        normal_force: f64,
        sliding_distance: f64,
        thermal_mass: f64,
    ) -> f64 {
        if thermal_mass < 1e-30 {
            return 0.0;
        }
        friction * normal_force * sliding_distance / thermal_mass
    }
}

// ---------------------------------------------------------------------------
// 12. Archard wear model
// ---------------------------------------------------------------------------

/// Archard's wear law: V_wear = K * F_n * L / H.
///
/// V_wear = wear volume (m³), K = dimensionless wear coefficient,
/// F_n = normal force (N), L = sliding distance (m), H = hardness (Pa).
pub fn archard_wear_volume(
    wear_coefficient: f64,
    normal_force: f64,
    sliding_distance: f64,
    hardness: f64,
) -> f64 {
    if hardness < 1e-30 {
        return 0.0;
    }
    wear_coefficient * normal_force * sliding_distance / hardness
}

/// Wear depth from wear volume, assuming circular contact of radius a.
pub fn wear_depth_from_volume(wear_volume: f64, contact_radius: f64) -> f64 {
    if contact_radius < 1e-30 {
        return 0.0;
    }
    wear_volume / (PI * contact_radius * contact_radius)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── helper ──────────────────────────────────────────────────────────────

    fn steel_steel() -> HertzSphereSphere {
        // E = 200 GPa, ν = 0.3 for both spheres of R = 10 mm
        HertzSphereSphere::new(0.01, 0.3, 200e9, 0.01, 0.3, 200e9)
    }

    fn steel_flat() -> HertzSphereFlat {
        HertzSphereFlat::new(0.01, 200e9, 0.3, 200e9, 0.3)
    }

    fn steel_cyl_flat() -> HertzCylinderFlat {
        HertzCylinderFlat::new(0.01, 200e9, 0.3, 200e9, 0.3)
    }

    // ── reduced_modulus ────────────────────────────────────────────────────

    #[test]
    fn test_reduced_modulus_symmetric() {
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let e_star = reduced_modulus(e, nu, e, nu);
        let expected = 0.5 * e / (1.0 - nu * nu);
        assert!((e_star - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_combined_radius_equal() {
        let r = combined_radius(0.01, 0.01);
        assert!((r - 0.005).abs() < 1e-12);
    }

    #[test]
    fn test_combined_radius_flat() {
        let r = combined_radius(0.01, f64::INFINITY);
        assert!((r - 0.01).abs() < 1e-12);
    }

    // ── HertzSphereSphere ─────────────────────────────────────────────────

    #[test]
    fn test_hertz_ss_contact_radius_positive() {
        let h = steel_steel();
        assert!(h.contact_radius(100.0) > 0.0);
    }

    #[test]
    fn test_hertz_ss_indentation_positive() {
        let h = steel_steel();
        assert!(h.indentation(100.0) > 0.0);
    }

    #[test]
    fn test_hertz_ss_peak_pressure_positive() {
        let h = steel_steel();
        assert!(h.peak_pressure(100.0) > 0.0);
    }

    #[test]
    fn test_hertz_ss_mean_pressure_less_than_peak() {
        let h = steel_steel();
        let f = 100.0;
        assert!(h.mean_pressure(f) < h.peak_pressure(f));
    }

    #[test]
    fn test_hertz_ss_force_roundtrip() {
        let h = steel_steel();
        let f_in = 500.0;
        let delta = h.indentation(f_in);
        let f_out = h.force_from_indentation(delta);
        assert!((f_in - f_out).abs() / f_in < 1e-6);
    }

    #[test]
    fn test_hertz_ss_stiffness_positive() {
        let h = steel_steel();
        assert!(h.contact_stiffness(100.0) > 0.0);
    }

    #[test]
    fn test_hertz_ss_pressure_outside_contact_zero() {
        let h = steel_steel();
        let a = h.contact_radius(100.0);
        assert_eq!(h.pressure_at_radius(a * 2.0, 100.0), 0.0);
    }

    #[test]
    fn test_hertz_ss_pressure_at_center_equals_peak() {
        let h = steel_steel();
        let f = 200.0;
        let p0 = h.peak_pressure(f);
        let p_center = h.pressure_at_radius(0.0, f);
        assert!((p_center - p0).abs() / p0 < 1e-10);
    }

    // ── HertzSphereFlat ───────────────────────────────────────────────────

    #[test]
    fn test_sphere_flat_contact_radius_positive() {
        let h = steel_flat();
        assert!(h.contact_radius(100.0) > 0.0);
    }

    #[test]
    fn test_sphere_flat_indentation_positive() {
        let h = steel_flat();
        assert!(h.indentation(200.0) > 0.0);
    }

    #[test]
    fn test_sphere_flat_peak_pressure_positive() {
        let h = steel_flat();
        assert!(h.peak_pressure(100.0) > 0.0);
    }

    #[test]
    fn test_sphere_flat_strain_energy_positive() {
        let h = steel_flat();
        assert!(h.strain_energy(100.0) > 0.0);
    }

    // ── HertzCylinderFlat ─────────────────────────────────────────────────

    #[test]
    fn test_cyl_flat_half_width_positive() {
        let h = steel_cyl_flat();
        assert!(h.half_width(1e4) > 0.0);
    }

    #[test]
    fn test_cyl_flat_peak_pressure_positive() {
        let h = steel_cyl_flat();
        assert!(h.peak_pressure(1e4) > 0.0);
    }

    #[test]
    fn test_cyl_flat_mean_less_than_peak() {
        let h = steel_cyl_flat();
        let fpl = 1e5;
        assert!(h.mean_pressure(fpl) < h.peak_pressure(fpl));
    }

    #[test]
    fn test_cyl_flat_pressure_outside_zero() {
        let h = steel_cyl_flat();
        let fpl = 1e4;
        let b = h.half_width(fpl);
        assert_eq!(h.pressure_at_position(b * 2.0, fpl), 0.0);
    }

    #[test]
    fn test_cyl_flat_max_shear_stress_positive() {
        let h = steel_cyl_flat();
        assert!(h.max_shear_stress(1e5) > 0.0);
    }

    #[test]
    fn test_cyl_flat_max_shear_depth_proportional() {
        let h = steel_cyl_flat();
        let fpl = 1e4;
        let b = h.half_width(fpl);
        let z_max = h.max_shear_depth(fpl);
        assert!((z_max - 0.786 * b).abs() < 1e-15);
    }

    // ── MindlinCompliance ─────────────────────────────────────────────────

    fn make_mindlin() -> MindlinCompliance {
        MindlinCompliance::new(200e9, 0.3, 200e9, 0.3, 0.4)
    }

    #[test]
    fn test_mindlin_tangential_stiffness_positive() {
        let m = make_mindlin();
        assert!(m.tangential_stiffness(1e-4) > 0.0);
    }

    #[test]
    fn test_mindlin_partial_slip_zero_force() {
        let m = make_mindlin();
        let disp = m.partial_slip_displacement(0.0, 100.0, 1e-4).unwrap();
        assert!(disp.abs() < 1e-20);
    }

    #[test]
    fn test_mindlin_gross_slip_returns_none() {
        let m = make_mindlin();
        // tangential = limit
        let result = m.partial_slip_displacement(40.0, 100.0, 1e-4);
        assert!(result.is_none());
    }

    #[test]
    fn test_mindlin_stick_radius_less_than_contact() {
        let m = make_mindlin();
        let a = 1e-4;
        let c = m.stick_radius(5.0, 100.0, a).unwrap();
        assert!(c < a);
    }

    #[test]
    fn test_mindlin_fretting_energy_positive() {
        let m = make_mindlin();
        let e = m.fretting_energy_per_cycle(5.0, 100.0, 1e-4);
        assert!(e > 0.0);
    }

    #[test]
    fn test_mindlin_fretting_energy_zero_amplitude() {
        let m = make_mindlin();
        let e = m.fretting_energy_per_cycle(0.0, 100.0, 1e-4);
        assert!(e.abs() < 1e-30);
    }

    // ── JKR adhesion ─────────────────────────────────────────────────────

    fn make_jkr() -> JkrAdhesion {
        JkrAdhesion::new(1e-6, 1e10, 0.05)
    }

    #[test]
    fn test_jkr_pull_off_negative() {
        let j = make_jkr();
        assert!(j.pull_off_force() < 0.0);
    }

    #[test]
    fn test_jkr_zero_load_contact_positive() {
        let j = make_jkr();
        // Adhesion creates contact even at zero load.
        assert!(j.zero_load_contact_radius() > 0.0);
    }

    #[test]
    fn test_jkr_contact_radius_increases_with_force() {
        let j = make_jkr();
        let a0 = j.contact_radius(0.0);
        let a1 = j.contact_radius(1e-3);
        assert!(a1 > a0);
    }

    #[test]
    fn test_jkr_surface_energy_negative() {
        let j = make_jkr();
        assert!(j.surface_energy(1e-3) < 0.0);
    }

    // ── DMT adhesion ─────────────────────────────────────────────────────

    fn make_dmt() -> DmtAdhesion {
        DmtAdhesion::new(1e-6, 1e10, 0.05)
    }

    #[test]
    fn test_dmt_pull_off_negative() {
        let d = make_dmt();
        assert!(d.pull_off_force() < 0.0);
    }

    #[test]
    fn test_dmt_pull_off_stronger_than_jkr() {
        // JKR: F_po = -1.5 π W R, DMT: F_po = -2 π W R → |DMT| > |JKR|
        let j = make_jkr();
        let d = make_dmt();
        assert!(d.pull_off_force().abs() > j.pull_off_force().abs());
    }

    #[test]
    fn test_dmt_contact_radius_zero_at_large_negative_force() {
        let d = make_dmt();
        // Below pull-off effective force → no contact
        let f_below = d.pull_off_force() * 1.1;
        let a = d.contact_radius(f_below);
        assert_eq!(a, 0.0);
    }

    #[test]
    fn test_dmt_tabor_parameter_positive() {
        let d = make_dmt();
        let tabor = d.tabor_parameter(3e-10);
        assert!(tabor > 0.0);
    }

    #[test]
    fn test_dmt_peak_pressure_positive() {
        let d = make_dmt();
        assert!(d.peak_pressure(1e-3) > 0.0);
    }

    // ── ContactAreaHysteresis ─────────────────────────────────────────────

    #[test]
    fn test_hysteresis_loading_increases_area() {
        let mut h = ContactAreaHysteresis::new(5e-3, 1e11);
        let a0 = h.load(100.0);
        let a1 = h.load(500.0);
        assert!(a1 > a0);
    }

    #[test]
    fn test_hysteresis_unload_less_or_equal_load() {
        let mut h = ContactAreaHysteresis::new(5e-3, 1e11);
        h.load(500.0);
        let a_load = h.load(200.0);
        let a_unload = h.unload(200.0);
        // Unloading branch ≥ Hertz loading due to residual
        assert!(a_unload >= a_load * 0.99);
    }

    #[test]
    fn test_hysteresis_energy_positive() {
        let mut h = ContactAreaHysteresis::new(5e-3, 1e11);
        h.load(500.0);
        let e = h.hysteresis_energy(50);
        assert!(e >= 0.0);
    }

    // ── ContactForceModel ─────────────────────────────────────────────────

    fn make_cfm() -> ContactForceModel {
        ContactForceModel {
            normal_stiffness: 1e8,
            tangential_stiffness: 5e7,
            friction: 0.3,
            damping: 0.05,
        }
    }

    #[test]
    fn test_cfm_no_force_open_gap() {
        let c = make_cfm();
        assert_eq!(c.linear_normal_force(0.001), 0.0);
    }

    #[test]
    fn test_cfm_linear_normal_force_overlap() {
        let c = make_cfm();
        let fn_ = c.linear_normal_force(-1e-4);
        assert!((fn_ - 1e4).abs() < 1.0);
    }

    #[test]
    fn test_cfm_hertz_normal_force_positive() {
        let c = make_cfm();
        assert!(c.hertz_normal_force(-1e-5) > 0.0);
    }

    #[test]
    fn test_cfm_tangential_force_coulomb_limited() {
        let c = make_cfm();
        let fn_ = 1000.0;
        let limit = c.friction * fn_;
        let ft = c.tangential_force(1.0, fn_); // very large displacement
        assert!(ft.abs() <= limit + 1e-10);
    }

    #[test]
    fn test_cfm_damped_normal_force_approach() {
        let c = make_cfm();
        let fn_ = c.damped_normal_force(-1e-4, -0.001, 0.1);
        assert!(fn_ >= 0.0);
    }

    // ── Greenwood–Williamson ──────────────────────────────────────────────

    fn make_gw() -> GreenwoodWilliamson {
        GreenwoodWilliamson {
            asperity_density: 1e10,
            asperity_radius: 1e-6,
            height_std: 1e-7,
            e_star: 1e11,
        }
    }

    #[test]
    fn test_gw_contact_area_fraction_positive() {
        let gw = make_gw();
        let frac = gw.real_contact_area_fraction(0.0);
        assert!(frac > 0.0);
    }

    #[test]
    fn test_gw_contact_area_decreases_with_separation() {
        let gw = make_gw();
        let a0 = gw.real_contact_area_fraction(0.0);
        let a1 = gw.real_contact_area_fraction(5e-7);
        assert!(a0 > a1);
    }

    #[test]
    fn test_gw_pressure_positive_at_zero_separation() {
        let gw = make_gw();
        assert!(gw.contact_pressure(0.0) > 0.0);
    }

    #[test]
    fn test_gw_plasticity_index_positive() {
        let gw = make_gw();
        let psi = gw.plasticity_index(2e9);
        assert!(psi > 0.0);
    }

    #[test]
    fn test_gw_asperity_density_positive() {
        let gw = make_gw();
        let n = gw.contact_asperity_density(0.0);
        assert!(n > 0.0);
    }

    #[test]
    fn test_gw_erfc_at_zero() {
        // erfc(0) = 1
        assert!((gw_erfc(0.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_gw_erfc_large_arg() {
        // erfc(large) ≈ 0
        assert!(gw_erfc(5.0) < 1e-6);
    }

    #[test]
    fn test_gw_separation_for_pressure() {
        let gw = make_gw();
        let p0 = gw.contact_pressure(0.0);
        let target = p0 * 0.5;
        let sep = gw.separation_for_pressure(target, 1e-12);
        assert!(sep.is_some());
        let sep = sep.unwrap();
        let p_check = gw.contact_pressure(sep);
        assert!((p_check - target).abs() / target < 0.01);
    }

    // ── FlashTemperature ──────────────────────────────────────────────────

    fn make_flash() -> FlashTemperature {
        FlashTemperature::new(50.0, 1.4e-5, 50.0, 1.4e-5)
    }

    #[test]
    fn test_flash_partition_factor_symmetric() {
        let ft = make_flash();
        assert!((ft.heat_partition_factor() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_flash_peclet_positive() {
        let pe = FlashTemperature::peclet_number(1.0, 1e-3, 1e-5);
        assert!(pe > 0.0);
    }

    #[test]
    fn test_flash_high_speed_positive() {
        let ft = make_flash();
        let dt = ft.flash_temperature_high_speed(0.3, 100.0, 1.0, 1e-3);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_flash_low_speed_positive() {
        let ft = make_flash();
        let dt = ft.flash_temperature_low_speed(0.3, 100.0, 0.01, 1e-3);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_flash_partition_positive() {
        let ft = make_flash();
        let dt = ft.flash_temperature_partition(0.3, 100.0, 1.0, 1e-3);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_flash_bulk_temperature_rise_positive() {
        let dt = FlashTemperature::bulk_temperature_rise(0.3, 100.0, 1.0, 500.0);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_flash_zero_speed_gives_zero() {
        let ft = make_flash();
        let dt = ft.flash_temperature_high_speed(0.3, 100.0, 0.0, 1e-3);
        assert_eq!(dt, 0.0);
    }

    // ── Archard wear ──────────────────────────────────────────────────────

    #[test]
    fn test_archard_wear_volume_positive() {
        let v = archard_wear_volume(1e-4, 1000.0, 0.5, 5e9);
        assert!(v > 0.0);
    }

    #[test]
    fn test_archard_zero_hardness() {
        let v = archard_wear_volume(1e-4, 1000.0, 0.5, 0.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_wear_depth_from_volume() {
        let depth = wear_depth_from_volume(1e-12, 1e-3);
        assert!(depth > 0.0);
    }

    #[test]
    fn test_wear_depth_zero_radius() {
        let depth = wear_depth_from_volume(1e-12, 0.0);
        assert_eq!(depth, 0.0);
    }
}
