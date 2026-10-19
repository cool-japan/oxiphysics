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
/// Johnson's convention (Contact Mechanics, tangential loading of spheres; Mindlin 1949):
///
/// 1/G* = (2 - ν₁)/G₁ + (2 - ν₂)/G₂
///
/// With this definition the initial tangential stiffness of a circular Hertzian
/// contact of radius `a` is k_t = 8 G* a.  For two identical bodies
/// G* = G / (2 (2 - ν)), i.e. k_t = 4 G a / (2 - ν).
pub fn mindlin_reduced_shear_modulus(g1: f64, nu1: f64, g2: f64, nu2: f64) -> f64 {
    1.0 / ((2.0 - nu1) / g1 + (2.0 - nu2) / g2)
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

    /// Tangential displacement at the onset of gross slip,
    /// δ_s = 3 μ F_n / (16 G* a).
    fn saturation_displacement(&self, normal_force: f64, contact_radius: f64) -> f64 {
        3.0 * self.friction * normal_force / (16.0 * self.g_star * contact_radius)
    }

    /// Mindlin partial-slip tangential displacement (loading from 0).
    ///
    /// δ_t = (3 μ F_n)/(16 G* a) \[1 - (1 - Q/(μ F_n))^(2/3)\]
    ///
    /// Its initial slope dδ_t/dQ at Q = 0 is 1/k_t with k_t = 8 G* a.
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
        let prefactor = self.saturation_displacement(normal_force, contact_radius);
        let term = 1.0 - tangential_force / limit;
        Some(prefactor * (1.0 - term.powf(2.0 / 3.0)))
    }

    /// Inverse: force from displacement (partial-slip).
    ///
    /// Q = μ F_n \[1 - (1 - δ_t/δ_s)^(3/2)\] with δ_s = 3 μ F_n/(16 G* a); the
    /// displacement is clamped to \[0, δ_s\] (gross-slip limit Q = μ F_n).
    pub fn force_from_displacement(
        &self,
        delta_t: f64,
        normal_force: f64,
        contact_radius: f64,
    ) -> f64 {
        let limit = self.friction * normal_force;
        let prefactor = self.saturation_displacement(normal_force, contact_radius);
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
    /// For a tangential force oscillating between ±Q* under constant normal
    /// force F_n, with q = Q*/(μ F_n) < 1 (Mindlin & Deresiewicz 1953; Johnson, Contact Mechanics):
    ///
    /// W = (36 μ² F_n²)/(5 k_t) \[1 - (1 - q)^(5/3) - (5q/6) (1 + (1 - q)^(2/3))\]
    ///
    /// which scales as Q*³ for small amplitudes (W → 2 Q*³/(9 μ F_n k_t)).  The
    /// amplitude is capped at the onset of gross slip (q = 1).
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
        let ratio = (tangential_amplitude.abs() / limit).min(1.0);
        mindlin_dissipation_per_cycle(limit, k_t, ratio)
    }
}

/// Bracket function of the Mindlin–Deresiewicz dissipation integral,
/// f(q) = 1 - (1 - q)^(5/3) - (5q/6) (1 + (1 - q)^(2/3)), for q ∈ \[0, 1\].
///
/// f(q) = 5q³/162 + 5q⁴/243 + 7q⁵/486 + … vanishes like q³, so the closed form
/// loses all significant digits to cancellation for small q.  Below `q = 0.2`
/// the (all-positive) binomial series is summed instead, which keeps the
/// relative accuracy at machine precision.
fn mindlin_dissipation_bracket(q: f64) -> f64 {
    if q <= 0.0 {
        return 0.0;
    }
    if q >= 0.2 {
        let rest = 1.0 - q;
        return 1.0 - rest.powf(5.0 / 3.0) - (5.0 * q / 6.0) * (1.0 + rest.powf(2.0 / 3.0));
    }
    let alpha_a = 5.0 / 3.0;
    let alpha_b = 2.0 / 3.0;
    // Binomial coefficients c_j(α) = C(α, j) (-1)^j of (1 - q)^α:
    //   c_a = c_k(5/3) and c_b = c_(k-1)(2/3) when entering iteration k.
    let mut c_a = 5.0 / 9.0; // c_2(5/3)
    let mut c_b = -2.0 / 3.0; // c_1(2/3)
    let mut q_pow = q * q; // q^2
    let mut sum = 0.0;
    for k in 3..80_u32 {
        let kf = f64::from(k);
        c_a *= (kf - 1.0 - alpha_a) / kf;
        c_b *= (kf - 2.0 - alpha_b) / (kf - 1.0);
        q_pow *= q;
        let term = (-c_a - 5.0 / 6.0 * c_b) * q_pow;
        sum += term;
        if term.abs() <= 1e-18 * sum.abs() {
            break;
        }
    }
    sum
}

/// Energy dissipated per full tangential load cycle of a Mindlin contact with
/// Coulomb limit `limit` = μ F_n, initial tangential stiffness `k_t` and
/// normalised force amplitude `ratio` = Q*/(μ F_n) ∈ \[0, 1\].
///
/// W = (36 · limit²)/(5 k_t) \[1 - (1 - q)^(5/3) - (5q/6) (1 + (1 - q)^(2/3))\]
///
/// Derived by integrating the Masing-rule hysteresis loop built from the
/// virgin loading curve δ(Q) = (3 μ F_n/(16 G* a)) \[1 - (1 - Q/(μ F_n))^(2/3)\].
pub fn mindlin_dissipation_per_cycle(limit: f64, k_t: f64, ratio: f64) -> f64 {
    let q = ratio.clamp(0.0, 1.0);
    36.0 * limit * limit / (5.0 * k_t) * mindlin_dissipation_bracket(q)
}

// ---------------------------------------------------------------------------
// 6. JKR adhesion model
// ---------------------------------------------------------------------------

/// JKR contact radius `a` for an applied load `force` (Johnson, Kendall &
/// Roberts 1971).
///
/// a³ = (3 R*/(4 E*)) \[F + 3π W R* + √(6π W R* F + (3π W R*)²)\]
///
/// where `r_star` is the effective radius R*, `e_star` the reduced modulus E*
/// and `work_adhesion` the Dupré work of adhesion W.  For W = 0 this reduces to
/// the Hertz radius (3 F R*/(4 E*))^(1/3); at zero load a³ = 9π W R*²/(2 E*).
///
/// Returns 0 when the surfaces are separated: for loads below the JKR pull-off
/// force −(3/2)π W R* the equation has no real solution, and for a vanishing or
/// negative bracket (Hertz limit under tension) there is no contact patch.
pub fn jkr_contact_radius(force: f64, r_star: f64, e_star: f64, work_adhesion: f64) -> f64 {
    let term = 3.0 * PI * work_adhesion * r_star;
    let discriminant = 6.0 * PI * work_adhesion * r_star * force + term * term;
    // Tolerate round-off exactly at the pull-off load, where the discriminant is 0.
    if discriminant < -1e-12 * term * term {
        return 0.0;
    }
    let a3 = 0.75 * (r_star / e_star) * (force + term + discriminant.max(0.0).sqrt());
    if a3 > 0.0 {
        a3.cbrt()
    } else if a3.is_nan() {
        f64::NAN
    } else {
        0.0
    }
}

/// JKR approach (indentation depth) at contact radius `a`:
/// δ = a²/R* − √(2π W a/E*).
pub fn jkr_indentation(a: f64, r_star: f64, e_star: f64, work_adhesion: f64) -> f64 {
    a * a / r_star - (2.0 * PI * work_adhesion * a / e_star).sqrt()
}

/// Elastic strain energy stored in a JKR contact of radius `a`.
///
/// The JKR pressure distribution is the Hertz pressure of radius `a` minus the
/// pressure under a rigid flat punch of the same radius pressed by
/// δ_p = √(2π W a/E*).  Clapeyron's theorem then gives (Hertz + cross + punch
/// parts)
///
/// U_E = 8 E* a⁵/(15 R*²) − 4 E* a³ δ_p/(3 R*) + 2π W a²
///
/// which is a positive-definite quadratic form in δ_p.  Its derivative with
/// respect to `a` at fixed approach reproduces the JKR equilibrium condition
/// (energy release rate = W).
pub fn jkr_elastic_energy_at_radius(a: f64, r_star: f64, e_star: f64, work_adhesion: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    let delta_punch = (2.0 * PI * work_adhesion * a / e_star).sqrt();
    let hertz = 8.0 * e_star * a.powi(5) / (15.0 * r_star * r_star);
    let cross = 4.0 * e_star * a.powi(3) * delta_punch / (3.0 * r_star);
    let punch = 2.0 * PI * work_adhesion * a * a;
    hertz - cross + punch
}

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
    /// a³ = (3R*/(4E*)) \[F + 3πW R* + √(6πW R* F + (3πW R*)²)\]
    /// where W = W_ad.  Equals the Hertz radius for W = 0 and is 0 below the
    /// pull-off load (see [`jkr_contact_radius`]).
    pub fn contact_radius(&self, force: f64) -> f64 {
        jkr_contact_radius(force, self.r_star, self.e_star, self.work_adhesion)
    }

    /// JKR approach (indentation) δ = a²/R* − √(2π W a/E*) at load F.
    pub fn indentation(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        jkr_indentation(a, self.r_star, self.e_star, self.work_adhesion)
    }

    /// JKR pull-off (adhesion) force: F_po = -(3/2) π W_ad R*.
    pub fn pull_off_force(&self) -> f64 {
        -1.5 * PI * self.work_adhesion * self.r_star
    }

    /// Indentation at zero applied load (contact due to adhesion only).
    pub fn zero_load_contact_radius(&self) -> f64 {
        self.contact_radius(0.0)
    }

    /// Elastic strain energy stored in the adhesive contact at load F.
    ///
    /// U_E = 8E*a⁵/(15R*²) − 4E*a³δ_p/(3R*) + 2πW a² with
    /// δ_p = √(2πW a/E*) (see [`jkr_elastic_energy_at_radius`]).
    pub fn elastic_energy(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        jkr_elastic_energy_at_radius(a, self.r_star, self.e_star, self.work_adhesion)
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

    #[test]
    fn test_mindlin_reduced_shear_modulus_identical_bodies() {
        let g = 80.0e9_f64;
        let nu = 0.3_f64;
        let g_star = mindlin_reduced_shear_modulus(g, nu, g, nu);
        let expected = g / (2.0 * (2.0 - nu));
        assert!(
            (g_star - expected).abs() / expected < 1e-14,
            "G* = {g_star}, expected {expected}"
        );
    }

    #[test]
    fn test_mindlin_tangential_stiffness_matches_johnson_for_identical_steel() {
        // Johnson: k_t = 8 a / ((2-nu1)/G1 + (2-nu2)/G2) = 4 G a / (2 - nu) for equal bodies.
        let (e, nu, a) = (200.0e9_f64, 0.3_f64, 1.0e-4_f64);
        let m = MindlinCompliance::new(e, nu, e, nu, 0.4);
        let g = e / (2.0 * (1.0 + nu));
        let expected = 4.0 * g * a / (2.0 - nu);
        let k_t = m.tangential_stiffness(a);
        assert!(
            (k_t - expected).abs() / expected < 1e-12,
            "k_t = {k_t}, expected 4Ga/(2-nu) = {expected}"
        );
    }

    #[test]
    fn test_mindlin_tangential_stiffness_dissimilar_bodies() {
        let (e1, nu1, e2, nu2, a) = (200.0e9_f64, 0.3_f64, 70.0e9_f64, 0.33_f64, 2.0e-4_f64);
        let m = MindlinCompliance::new(e1, nu1, e2, nu2, 0.3);
        let g1 = e1 / (2.0 * (1.0 + nu1));
        let g2 = e2 / (2.0 * (1.0 + nu2));
        let expected = 8.0 * a / ((2.0 - nu1) / g1 + (2.0 - nu2) / g2);
        let k_t = m.tangential_stiffness(a);
        assert!(
            (k_t - expected).abs() / expected < 1e-12,
            "k_t = {k_t}, expected {expected}"
        );
    }

    #[test]
    fn test_mindlin_initial_slope_equals_inverse_tangential_stiffness() {
        let m = make_mindlin();
        let (normal, a) = (100.0_f64, 1.0e-4_f64);
        let limit = m.friction * normal;
        let dq = 1.0e-6 * limit;
        let delta = m.partial_slip_displacement(dq, normal, a).unwrap();
        let k_fd = dq / delta;
        let k_t = m.tangential_stiffness(a);
        assert!(
            (k_fd - k_t).abs() / k_t < 1e-5,
            "finite-difference slope {k_fd} vs k_t {k_t}"
        );
    }

    #[test]
    fn test_mindlin_saturation_displacement_is_johnson_value() {
        // As Q -> mu F the displacement tends to 3 mu F / (16 G* a).
        let m = make_mindlin();
        let (normal, a) = (100.0_f64, 1.0e-4_f64);
        let limit = m.friction * normal;
        let delta = m
            .partial_slip_displacement(limit * (1.0 - 1.0e-12), normal, a)
            .unwrap();
        let expected = 3.0 * limit / (16.0 * m.g_star * a);
        assert!(
            (delta - expected).abs() / expected < 1e-6,
            "delta = {delta}, expected {expected}"
        );
    }

    #[test]
    fn test_mindlin_force_displacement_roundtrip() {
        let m = make_mindlin();
        let (normal, a) = (250.0_f64, 3.0e-4_f64);
        let limit = m.friction * normal;
        for q in [0.05, 0.3, 0.7, 0.95] {
            let force = q * limit;
            let delta = m.partial_slip_displacement(force, normal, a).unwrap();
            let back = m.force_from_displacement(delta, normal, a);
            assert!(
                (back - force).abs() / force < 1e-10,
                "q = {q}: force {force} came back as {back}"
            );
        }
        // Beyond the saturation displacement the force is capped at mu F_n.
        let delta_big = 10.0 * 3.0 * limit / (16.0 * m.g_star * a);
        assert!((m.force_from_displacement(delta_big, normal, a) - limit).abs() < 1e-12);
        assert_eq!(m.force_from_displacement(-1.0, normal, a), 0.0);
    }

    #[test]
    fn test_mindlin_fretting_energy_small_amplitude_is_cubic() {
        // W -> 2 Q*^3 / (9 mu F k_t) as q -> 0 (not linear in Q*).
        let m = make_mindlin();
        let (normal, a) = (100.0_f64, 1.0e-4_f64);
        let limit = m.friction * normal;
        let k_t = m.tangential_stiffness(a);
        for q in [1.0e-6, 1.0e-4, 1.0e-3] {
            let amplitude = q * limit;
            let w = m.fretting_energy_per_cycle(amplitude, normal, a);
            let asymptote = 2.0 * amplitude.powi(3) / (9.0 * limit * k_t);
            let tol = 2.0 * q + 1.0e-9;
            assert!(
                (w - asymptote).abs() / asymptote < tol,
                "q = {q}: W = {w:e}, asymptote = {asymptote:e}"
            );
        }
        // Doubling a small amplitude multiplies the dissipation by ~8.
        let w1 = m.fretting_energy_per_cycle(1.0e-5 * limit, normal, a);
        let w2 = m.fretting_energy_per_cycle(2.0e-5 * limit, normal, a);
        assert!((w2 / w1 - 8.0).abs() < 1.0e-3, "ratio = {}", w2 / w1);
    }

    #[test]
    fn test_mindlin_fretting_energy_matches_masing_hysteresis_loop() {
        // Independent check: loop area W = 4 Q* d* - 8 * integral_0^Q* d(Q) dQ built from the
        // virgin loading curve (Masing rule), integrated with Simpson's rule.
        let m = make_mindlin();
        let (normal, a) = (100.0_f64, 1.0e-4_f64);
        let limit = m.friction * normal;
        for q_star in [0.2_f64, 0.6, 0.9] {
            let amp = q_star * limit;
            let d_star = m.partial_slip_displacement(amp, normal, a).unwrap();
            let n = 4000_usize;
            let h = amp / n as f64;
            let mut simpson = 0.0;
            for i in 0..=n {
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                let d = m
                    .partial_slip_displacement(i as f64 * h, normal, a)
                    .unwrap();
                simpson += weight * d;
            }
            simpson *= h / 3.0;
            let loop_area = 4.0 * amp * d_star - 8.0 * simpson;
            let w = m.fretting_energy_per_cycle(amp, normal, a);
            assert!(
                (w - loop_area).abs() / loop_area < 1.0e-6,
                "q* = {q_star}: closed form {w:e} vs hysteresis loop {loop_area:e}"
            );
        }
    }

    #[test]
    fn test_mindlin_fretting_energy_gross_slip_limit_and_monotonicity() {
        let m = make_mindlin();
        let (normal, a) = (100.0_f64, 1.0e-4_f64);
        let limit = m.friction * normal;
        let k_t = m.tangential_stiffness(a);
        // q = 1: W = 6/5 (mu F)^2 / k_t; larger amplitudes are capped there.
        let expected = 1.2 * limit * limit / k_t;
        let w1 = m.fretting_energy_per_cycle(limit, normal, a);
        let w2 = m.fretting_energy_per_cycle(3.0 * limit, normal, a);
        assert!((w1 - expected).abs() / expected < 1e-12);
        assert_eq!(w1, w2);
        let mut previous = 0.0;
        for i in 1..=20 {
            let w = m.fretting_energy_per_cycle(i as f64 * 0.05 * limit, normal, a);
            assert!(w > previous, "dissipation must increase with amplitude");
            previous = w;
        }
    }

    #[test]
    fn test_mindlin_dissipation_series_matches_closed_form() {
        // Just below the series/closed-form switch (q = 0.2) both branches must agree.
        for q in [0.05_f64, 0.12, 0.19] {
            let rest = 1.0 - q;
            let closed =
                1.0 - rest.powf(5.0 / 3.0) - (5.0 * q / 6.0) * (1.0 + rest.powf(2.0 / 3.0));
            let series = mindlin_dissipation_bracket(q);
            assert!(
                (series - closed).abs() / closed < 1e-9,
                "q = {q}: series {series:e} vs closed form {closed:e}"
            );
        }
        // Leading behaviour 5 q^3 / 162 for tiny q (closed form is pure round-off here).
        let q = 1.0e-9_f64;
        let series = mindlin_dissipation_bracket(q);
        assert!((series / (5.0 * q.powi(3) / 162.0) - 1.0).abs() < 1e-8);
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

    #[test]
    fn test_jkr_reduces_to_hertz_without_adhesion() {
        let (r_star, e_star) = (5.0e-3_f64, 1.0e11_f64);
        let jkr = JkrAdhesion::new(r_star, e_star, 0.0);
        let hertz = HertzSphereSphere { r_star, e_star };
        for force in [0.1, 1.0, 50.0, 1.0e3] {
            let a_jkr = jkr.contact_radius(force);
            let a_hertz = hertz.contact_radius(force);
            assert!(
                (a_jkr - a_hertz).abs() / a_hertz < 1e-12,
                "F = {force}: JKR {a_jkr} vs Hertz {a_hertz}"
            );
        }
    }

    #[test]
    fn test_jkr_zero_load_radius_cubed() {
        // a0^3 = 9 pi W R^2 / (2 E*)
        let j = make_jkr();
        let a0 = j.contact_radius(0.0);
        let expected = 9.0 * PI * j.work_adhesion * j.r_star * j.r_star / (2.0 * j.e_star);
        assert!(
            (a0.powi(3) - expected).abs() / expected < 1e-12,
            "a0^3 = {:e}, expected {expected:e}",
            a0.powi(3)
        );
    }

    #[test]
    fn test_jkr_radius_satisfies_jkr_load_equation() {
        // P = 4 E* a^3 / (3 R) - sqrt(8 pi W E* a^3) must return the applied load.
        let j = make_jkr();
        let pull_off = j.pull_off_force().abs();
        for force in [-0.99, -0.5, 0.0, 0.7, 4.0, 60.0].map(|f| f * pull_off) {
            let a = j.contact_radius(force);
            let load = 4.0 * j.e_star * a.powi(3) / (3.0 * j.r_star)
                - (8.0 * PI * j.work_adhesion * j.e_star * a.powi(3)).sqrt();
            let scale = pull_off.max(force.abs());
            assert!(
                (load - force).abs() / scale < 1e-9,
                "F = {force:e}: load equation gives {load:e}"
            );
        }
    }

    #[test]
    fn test_jkr_pull_off_radius_and_separation() {
        let j = make_jkr();
        let a0 = j.contact_radius(0.0);
        let a_po = j.contact_radius(j.pull_off_force());
        // a_po^3 = a0^3 / 4
        assert!((a_po.powi(3) / a0.powi(3) - 0.25).abs() < 1e-9);
        // Below the pull-off load there is no contact patch (and never a negative radius).
        assert_eq!(j.contact_radius(1.01 * j.pull_off_force()), 0.0);
        assert_eq!(j.contact_radius(10.0 * j.pull_off_force()), 0.0);
        assert_eq!(j.elastic_energy(10.0 * j.pull_off_force()), 0.0);
    }

    #[test]
    fn test_jkr_indentation_formula() {
        let j = make_jkr();
        let force = 2.0e-3;
        let a = j.contact_radius(force);
        let expected = a * a / j.r_star - (2.0 * PI * j.work_adhesion * a / j.e_star).sqrt();
        assert!((j.indentation(force) - expected).abs() <= 1e-12 * expected.abs());
        // Without adhesion the JKR approach is the Hertz approach a^2/R.
        let hertz_like = JkrAdhesion::new(j.r_star, j.e_star, 0.0);
        let a_h = hertz_like.contact_radius(force);
        assert!((hertz_like.indentation(force) - a_h * a_h / j.r_star).abs() < 1e-18);
    }

    #[test]
    fn test_jkr_elastic_energy_hertz_limit() {
        // W = 0: U_E = (2/5) F delta = 8 E* a^5 / (15 R^2).
        let (r_star, e_star) = (5.0e-3_f64, 1.0e11_f64);
        let jkr = JkrAdhesion::new(r_star, e_star, 0.0);
        let hertz = HertzSphereSphere { r_star, e_star };
        let force = 40.0;
        let expected = 0.4 * force * hertz.indentation(force);
        let u = jkr.elastic_energy(force);
        assert!(
            (u - expected).abs() / expected < 1e-12,
            "U_E = {u:e}, expected {expected:e}"
        );
    }

    #[test]
    fn test_jkr_elastic_energy_zero_load_closed_form() {
        // At F = 0 the Hertz and flat-punch loads cancel: U_E = 4 E* a0^5 / (45 R^2).
        let j = make_jkr();
        let a0 = j.contact_radius(0.0);
        let expected = 4.0 * j.e_star * a0.powi(5) / (45.0 * j.r_star * j.r_star);
        let u = j.elastic_energy(0.0);
        assert!(
            (u - expected).abs() / expected < 1e-10,
            "U_E(0) = {u:e}, expected {expected:e}"
        );
    }

    #[test]
    fn test_jkr_elastic_energy_satisfies_griffith_equilibrium() {
        // General (non-equilibrium) energy of the Hertz-minus-punch superposition at fixed approach:
        //   U(a, d) = 8E*a^5/(15R^2) - 4E*a^3 d2/(3R) + E* a d2^2,  d2 = a^2/R - d.
        // At the JKR approach it must equal elastic_energy(), and dU/da|_d must equal
        // W * d(pi a^2)/da = 2 pi W a (energy release rate = work of adhesion).
        let j = make_jkr();
        let u_general = |a: f64, d: f64| -> f64 {
            let d2 = a * a / j.r_star - d;
            8.0 * j.e_star * a.powi(5) / (15.0 * j.r_star * j.r_star)
                - 4.0 * j.e_star * a.powi(3) * d2 / (3.0 * j.r_star)
                + j.e_star * a * d2 * d2
        };
        let pull_off = j.pull_off_force().abs();
        for force in [-0.9, 0.0, 1.0, 5.0, 20.0].map(|f| f * pull_off) {
            let a = j.contact_radius(force);
            let d = j.indentation(force);
            let u = j.elastic_energy(force);
            assert!(
                (u_general(a, d) - u).abs() / u < 1e-9,
                "F = {force}: general {:e} vs closed form {u:e}",
                u_general(a, d)
            );
            // Central difference; the truncation error is O((h/a)^2) and stays below 1e-6 here.
            let h = 1.0e-4 * a;
            let du_da = (u_general(a + h, d) - u_general(a - h, d)) / (2.0 * h);
            let expected = 2.0 * PI * j.work_adhesion * a;
            assert!(
                (du_da - expected).abs() / expected < 1e-5,
                "F = {force}: dU/da = {du_da:e}, 2 pi W a = {expected:e}"
            );
            assert!(u > 0.0, "strain energy must be positive");
        }
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
