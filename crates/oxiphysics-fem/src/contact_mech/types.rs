//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// Discretized contact patch for numerical pressure distribution.
pub struct ContactPatch {
    /// Radial positions of discretization points (m).
    pub radii: Vec<f64>,
    /// Pressure at each point (Pa).
    pub pressures: Vec<f64>,
    /// Contact radius (m).
    pub contact_radius: f64,
}
impl ContactPatch {
    /// Create a Hertzian pressure patch with n discretization points.
    pub fn hertzian(contact_radius: f64, max_pressure: f64, n_points: usize) -> Self {
        let mut radii = Vec::with_capacity(n_points);
        let mut pressures = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let r = contact_radius * (i as f64) / (n_points as f64 - 1.0).max(1.0);
            let ratio = r / contact_radius;
            let p = if ratio <= 1.0 {
                max_pressure * (1.0 - ratio * ratio).max(0.0).sqrt()
            } else {
                0.0
            };
            radii.push(r);
            pressures.push(p);
        }
        Self {
            radii,
            pressures,
            contact_radius,
        }
    }
    /// Total force by numerical integration of the pressure profile (axisymmetric).
    ///
    /// F = integral 2*pi*r*p(r) dr
    pub fn total_force(&self) -> f64 {
        let n = self.radii.len();
        if n < 2 {
            return 0.0;
        }
        let mut force = 0.0;
        for i in 1..n {
            let dr = self.radii[i] - self.radii[i - 1];
            let r_mid = 0.5 * (self.radii[i] + self.radii[i - 1]);
            let p_mid = 0.5 * (self.pressures[i] + self.pressures[i - 1]);
            force += 2.0 * PI * r_mid * p_mid * dr;
        }
        force
    }
    /// Mean pressure over the contact patch.
    pub fn mean_pressure(&self) -> f64 {
        let area = PI * self.contact_radius * self.contact_radius;
        if area < 1e-30 {
            return 0.0;
        }
        self.total_force() / area
    }
    /// Maximum pressure in the patch.
    pub fn max_pressure(&self) -> f64 {
        self.pressures.iter().cloned().fold(0.0f64, f64::max)
    }
}
/// Linear elastic solid defined by Young's modulus and Poisson's ratio.
pub struct ElasticSolid {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio nu (dimensionless).
    pub poisson_ratio: f64,
}
impl ElasticSolid {
    /// Shear modulus G = E / (2*(1+nu)).
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Bulk modulus K = E / (3*(1-2*nu)).
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
    /// Combined (reduced) modulus E*:
    ///   1/E* = (1 - nu1^2)/E1 + (1 - nu2^2)/E2
    pub fn reduced_modulus(s1: &ElasticSolid, s2: &ElasticSolid) -> f64 {
        let inv = (1.0 - s1.poisson_ratio * s1.poisson_ratio) / s1.youngs_modulus
            + (1.0 - s2.poisson_ratio * s2.poisson_ratio) / s2.youngs_modulus;
        1.0 / inv
    }
    /// Plane-strain modulus E' = E / (1 - nu^2).
    pub fn plane_strain_modulus(&self) -> f64 {
        self.youngs_modulus / (1.0 - self.poisson_ratio * self.poisson_ratio)
    }
}
/// Fretting contact analysis.
pub struct FrettingContact {
    /// Mindlin contact.
    pub mindlin: MindlinContact,
    /// Number of loading cycles.
    pub cycles: u64,
    /// Tangential displacement amplitude (m).
    pub amplitude: f64,
}
impl FrettingContact {
    /// Create a new fretting contact.
    pub fn new(mindlin: MindlinContact, amplitude: f64) -> Self {
        Self {
            mindlin,
            cycles: 0,
            amplitude,
        }
    }
    /// Total energy dissipated over all cycles.
    pub fn total_dissipated_energy(&self, normal_force: f64) -> f64 {
        let per_cycle = self
            .mindlin
            .fretting_energy_per_cycle(self.amplitude, normal_force);
        per_cycle * self.cycles as f64
    }
    /// Advance by one cycle.
    pub fn advance_cycle(&mut self) {
        self.cycles += 1;
    }
    /// Advance by n cycles.
    pub fn advance_cycles(&mut self, n: u64) {
        self.cycles += n;
    }
    /// Fretting wear depth estimate (Archard's law simplified).
    ///
    /// h = K * p_mean * delta_t * N / H
    /// where K is the wear coefficient, H is hardness.
    pub fn wear_depth(
        wear_coefficient: f64,
        mean_pressure: f64,
        sliding_distance_per_cycle: f64,
        num_cycles: u64,
        hardness: f64,
    ) -> f64 {
        if hardness < 1e-30 {
            return 0.0;
        }
        wear_coefficient * mean_pressure * sliding_distance_per_cycle * num_cycles as f64 / hardness
    }
}
/// Maugis-Dugdale adhesion contact — a transition model between DMT and JKR.
///
/// Characterized by the Maugis parameter λ:
///   λ < 0.1 → DMT-like
///   λ > 5.0 → JKR-like
pub struct MaugisDugdaleContact {
    /// Effective radius R* (m).
    pub effective_radius: f64,
    /// Reduced modulus E* (Pa).
    pub reduced_modulus: f64,
    /// Work of adhesion W (J/m²).
    pub work_of_adhesion: f64,
    /// Cohesive stress σ₀ (Pa) in Dugdale zone.
    pub cohesive_stress: f64,
}
impl MaugisDugdaleContact {
    /// Maugis parameter λ = σ₀ * (R* / (π * W * E*²))^(1/3).
    pub fn maugis_parameter(&self) -> f64 {
        let denom = PI * self.work_of_adhesion * self.reduced_modulus * self.reduced_modulus;
        if denom < 1e-60 {
            return 0.0;
        }
        self.cohesive_stress * (self.effective_radius / denom).cbrt()
    }
    /// Approximate JKR pull-off force: -(3/2)πWR*.
    pub fn jkr_pull_off_force(&self) -> f64 {
        -1.5 * PI * self.work_of_adhesion * self.effective_radius
    }
    /// Approximate DMT pull-off force: -2πWR*.
    pub fn dmt_pull_off_force(&self) -> f64 {
        -2.0 * PI * self.work_of_adhesion * self.effective_radius
    }
    /// Interpolated pull-off force using Maugis parameter.
    ///
    /// Uses a heuristic: F_po = F_dmt + (F_jkr - F_dmt) * f(λ)
    /// where f(λ) = 1 - exp(-λ).
    pub fn pull_off_force(&self) -> f64 {
        let lambda = self.maugis_parameter();
        let f_dmt = self.dmt_pull_off_force();
        let f_jkr = self.jkr_pull_off_force();
        f_dmt + (f_jkr - f_dmt) * (1.0 - (-lambda).exp())
    }
    /// Transition contact radius under zero applied force (approximation).
    pub fn zero_force_contact_radius(&self) -> f64 {
        let lambda = self.maugis_parameter();
        let a_dmt = (6.0 * PI * self.work_of_adhesion * self.effective_radius.powi(2)
            / self.reduced_modulus)
            .cbrt();
        let a_jkr = (9.0 * PI * self.work_of_adhesion * self.effective_radius.powi(2)
            / (4.0 * self.reduced_modulus))
            .cbrt();
        let t = (1.0 - (-lambda).exp()).clamp(0.0, 1.0);
        a_dmt + t * (a_jkr - a_dmt)
    }
}
/// Rolling contact mechanics (cylinder or sphere on flat).
pub struct RollingContact {
    /// Hertz contact geometry.
    pub hertz: HertzContact,
    /// Rolling radius (m).
    pub rolling_radius: f64,
    /// Rolling friction coefficient.
    pub rolling_friction: f64,
}
impl RollingContact {
    /// Create a new rolling contact.
    pub fn new(hertz: HertzContact, rolling_radius: f64, rolling_friction: f64) -> Self {
        Self {
            hertz,
            rolling_radius,
            rolling_friction,
        }
    }
    /// Rolling resistance force: F_r = mu_r * F_n.
    pub fn rolling_resistance(&self, normal_force: f64) -> f64 {
        self.rolling_friction * normal_force
    }
    /// Rolling resistance torque: T_r = mu_r * F_n * R.
    pub fn rolling_torque(&self, normal_force: f64) -> f64 {
        self.rolling_friction * normal_force * self.rolling_radius
    }
    /// Creep ratio for a given slip velocity and rolling velocity.
    ///
    /// xi = v_slip / v_rolling
    pub fn creep_ratio(slip_velocity: f64, rolling_velocity: f64) -> f64 {
        if rolling_velocity.abs() < 1e-30 {
            return 0.0;
        }
        slip_velocity / rolling_velocity
    }
    /// Carter's 2-D rolling contact creep force (longitudinal).
    ///
    /// F_x = -mu * F_n * (1 - (1 - xi/xi_s)^2) for xi < xi_s
    /// F_x = -mu * F_n                           for xi >= xi_s
    ///
    /// where xi_s = 3*mu*F_n / (8*G*a*c_11) is the saturation creep.
    pub fn carter_creep_force(&self, creep: f64, normal_force: f64, friction: f64) -> f64 {
        let a = self.hertz.contact_radius(normal_force);
        let g_star = self.hertz.solid_a.shear_modulus();
        let xi_s = if g_star * a > 1e-30 {
            3.0 * friction * normal_force / (8.0 * g_star * a)
        } else {
            1e-6
        };
        let xi = creep.abs();
        let force = if xi < xi_s {
            friction * normal_force * (1.0 - (1.0 - xi / xi_s).powi(2))
        } else {
            friction * normal_force
        };
        -force * creep.signum()
    }
    /// Contact patch half-width for a 2-D (line) contact.
    ///
    /// b = sqrt(4*F*R / (pi*E'*L))
    /// where L is the contact length.
    pub fn line_contact_half_width(force_per_length: f64, radius: f64, e_star: f64) -> f64 {
        (4.0 * force_per_length * radius / (PI * e_star)).sqrt()
    }
}
/// Mindlin tangential contact built on top of a Hertz contact problem.
pub struct MindlinContact {
    /// Underlying Hertz contact geometry/properties.
    pub hertz: HertzContact,
    /// Coefficient of friction mu.
    pub friction: f64,
}
impl MindlinContact {
    /// Reduced shear modulus G*:
    ///   1/G* = (2 - nu1) / (4*G1) + (2 - nu2) / (4*G2)
    fn reduced_shear_modulus(&self) -> f64 {
        let g1 = self.hertz.solid_a.shear_modulus();
        let g2 = self.hertz.solid_b.shear_modulus();
        let nu1 = self.hertz.solid_a.poisson_ratio;
        let nu2 = self.hertz.solid_b.poisson_ratio;
        let inv = (2.0 - nu1) / (4.0 * g1) + (2.0 - nu2) / (4.0 * g2);
        1.0 / inv
    }
    /// Tangential stiffness k_t = 8 * G* * a.
    pub fn tangential_stiffness(&self, normal_force: f64) -> f64 {
        let g_star = self.reduced_shear_modulus();
        let a = self.hertz.contact_radius(normal_force);
        8.0 * g_star * a
    }
    /// Mindlin tangential compliance (partial-slip regime):
    ///   delta_t = (3*mu*F_n / (8*G**a)) * \[1 - (1 - Ft/(mu*Fn))^(2/3)\]
    ///
    /// Returns `None` if `tangential_force >= mu * normal_force` (gross slip).
    pub fn tangential_compliance(&self, tangential_force: f64, normal_force: f64) -> Option<f64> {
        let limit = self.friction * normal_force;
        if tangential_force >= limit {
            return None;
        }
        let g_star = self.reduced_shear_modulus();
        let a = self.hertz.contact_radius(normal_force);
        let prefactor = 3.0 * limit / (8.0 * g_star * a);
        let term = 1.0 - (tangential_force / limit);
        Some(prefactor * (1.0 - term.powf(2.0 / 3.0)))
    }
    /// Inverse Mindlin relation: given tangential displacement delta_t,
    /// returns the tangential force (partial-slip regime).
    pub fn partial_slip_force_disp(&self, delta_t: f64, normal_force: f64) -> f64 {
        let limit = self.friction * normal_force;
        let g_star = self.reduced_shear_modulus();
        let a = self.hertz.contact_radius(normal_force);
        let prefactor = 3.0 * limit / (8.0 * g_star * a);
        let ratio = (delta_t / prefactor).min(1.0);
        limit * (1.0 - (1.0 - ratio).powf(1.5))
    }
    /// Slip annulus inner radius c = a * (1 - Ft/(mu*Fn))^(1/3).
    ///
    /// Returns `None` if in gross slip.
    pub fn stick_radius(&self, tangential_force: f64, normal_force: f64) -> Option<f64> {
        let limit = self.friction * normal_force;
        if tangential_force >= limit || limit < 1e-30 {
            return None;
        }
        let a = self.hertz.contact_radius(normal_force);
        Some(a * (1.0 - tangential_force / limit).powf(1.0 / 3.0))
    }
    /// Energy dissipated per cycle in a Mindlin fretting contact.
    ///
    /// W_d = (9 * mu^2 * F_n^2) / (10 * k_t) * \[1 - (1 - Q/(mu*F_n))^(5/3)\]
    /// where Q is the tangential force amplitude.
    pub fn fretting_energy_per_cycle(&self, tangential_amplitude: f64, normal_force: f64) -> f64 {
        let limit = self.friction * normal_force;
        if limit < 1e-30 {
            return 0.0;
        }
        let k_t = self.tangential_stiffness(normal_force);
        if k_t < 1e-30 {
            return 0.0;
        }
        let ratio = (tangential_amplitude / limit).min(1.0);
        (9.0 * limit * limit) / (10.0 * k_t) * (1.0 - (1.0 - ratio).powf(5.0 / 3.0))
    }
}
/// Analytical subsurface stress field for Hertzian contact.
pub struct HertzianStressField;
impl HertzianStressField {
    /// Maximum shear stress at depth z below the contact centre.
    pub fn subsurface_max_shear_stress(depth: f64, a: f64, p0: f64) -> f64 {
        if a < f64::EPSILON {
            return 0.0;
        }
        let xi = depth / a;
        let atan_inv_xi = if xi < 1e-10 {
            PI / 2.0
        } else {
            (1.0 / xi).atan()
        };
        let denom = (1.0 + xi * xi).sqrt();
        0.5 * p0 * ((1.0 + 2.0 * xi * xi) * atan_inv_xi / denom - 3.0 * xi / 2.0 / denom)
    }
    /// Depth at which the maximum subsurface shear stress occurs:  z_max ≈ 0.48 * a.
    pub fn max_shear_depth(a: f64) -> f64 {
        0.48 * a
    }
    /// Normal stress sigma_z on the axis (r=0) at depth z.
    ///
    /// sigma_z = -p0 / (1 + (z/a)^2)
    pub fn axial_normal_stress(depth: f64, a: f64, p0: f64) -> f64 {
        if a < f64::EPSILON {
            return 0.0;
        }
        let xi = depth / a;
        -p0 / (1.0 + xi * xi)
    }
    /// Von Mises stress at the surface (z = 0) as a function of radial position r.
    pub fn von_mises_at_surface(r: f64, a: f64, p0: f64, nu: f64) -> f64 {
        if a < f64::EPSILON || r > a {
            return 0.0;
        }
        let rho = r / a;
        let rho2 = rho * rho;
        let sq = (1.0 - rho2).sqrt();
        let sigma_z = -p0 * sq;
        let sigma_r;
        let sigma_theta;
        if rho < f64::EPSILON {
            sigma_r = p0 * (-(1.0 - 2.0 * nu) / 2.0 - 1.0);
            sigma_theta = sigma_r;
        } else {
            let factor = (1.0 - 2.0 * nu) / (3.0 * rho2);
            let inner = 1.0 - (1.0 - rho2).powf(1.5);
            sigma_r = p0 * (factor * inner - sq);
            sigma_theta = -p0 * (2.0 * nu * sq + (1.0 - 2.0 * nu) * (1.0 - inner / (rho2)));
            let diff_rz = sigma_r - sigma_z;
            let diff_zt = sigma_z - sigma_theta;
            let diff_tr = sigma_theta - sigma_r;
            return (0.5 * (diff_rz * diff_rz + diff_zt * diff_zt + diff_tr * diff_tr)).sqrt();
        }
        let diff_rz = sigma_r - sigma_z;
        let diff_zt = sigma_z - sigma_theta;
        let diff_tr = sigma_theta - sigma_r;
        (0.5 * (diff_rz * diff_rz + diff_zt * diff_zt + diff_tr * diff_tr)).sqrt()
    }
    /// Maximum subsurface von Mises stress (approximate).
    ///
    /// Occurs at z ≈ 0.48*a with value ≈ 0.31 * p0 (for nu = 0.3).
    pub fn max_subsurface_von_mises(a: f64, p0: f64) -> f64 {
        0.31 * p0 * if a > 0.0 { 1.0 } else { 0.0 }
    }
}
/// Contact loading/unloading curve analysis.
pub struct ContactLoadDisplacement {
    /// Hertz contact.
    pub hertz: HertzContact,
    /// Load steps (N).
    pub loads: Vec<f64>,
}
impl ContactLoadDisplacement {
    /// Compute displacement at each load step.
    pub fn displacements(&self) -> Vec<f64> {
        self.loads
            .iter()
            .map(|&f| self.hertz.hertz_indentation(f))
            .collect()
    }
    /// Compute contact radius at each load step.
    pub fn contact_radii(&self) -> Vec<f64> {
        self.loads
            .iter()
            .map(|&f| self.hertz.contact_radius(f))
            .collect()
    }
    /// Compute stiffness dF/dd at each load step.
    pub fn stiffnesses(&self) -> Vec<f64> {
        self.loads
            .iter()
            .map(|&f| self.hertz.contact_stiffness(f))
            .collect()
    }
    /// Peak load.
    pub fn peak_load(&self) -> f64 {
        self.loads.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
    /// Total elastic work = integral F dd (trapezoidal).
    pub fn elastic_work(&self) -> f64 {
        let disps = self.displacements();
        let n = self.loads.len();
        if n < 2 {
            return 0.0;
        }
        let mut work = 0.0;
        for i in 1..n {
            let dd = disps[i] - disps[i - 1];
            let f_mid = 0.5 * (self.loads[i] + self.loads[i - 1]);
            work += f_mid * dd;
        }
        work
    }
}
/// FEM penalty-method contact constraint.
pub struct PenaltyContactConstraint {
    /// Normal penalty stiffness k_n (N/m).
    pub penalty: f64,
    /// Friction coefficient mu.
    pub friction_coef: f64,
}
impl PenaltyContactConstraint {
    /// Normal contact force from gap:  F_n = max(0, -gap) * penalty_n.
    pub fn normal_force(gap: f64, penalty: f64) -> f64 {
        if gap < 0.0 { -gap * penalty } else { 0.0 }
    }
    /// Contact force vector:
    ///   - normal force from gap,
    ///   - tangential force limited by Coulomb cone.
    pub fn contact_force(
        gap: f64,
        tangential_disp: f64,
        penalty_n: f64,
        penalty_t: f64,
        mu: f64,
    ) -> ([f64; 3], f64) {
        let fn_ = Self::normal_force(gap, penalty_n);
        let ft_trial = penalty_t * tangential_disp;
        let limit = mu * fn_;
        let ft = if ft_trial.abs() <= limit {
            ft_trial
        } else {
            limit * ft_trial.signum()
        };
        ([ft, 0.0, 0.0], fn_)
    }
    /// Augmented Lagrangian normal contact force.
    ///
    /// F_n = max(0, lambda + epsilon * g_n)
    /// where lambda is the Lagrange multiplier and g_n is the gap.
    pub fn augmented_lagrangian_force(lambda: f64, gap: f64, epsilon: f64) -> f64 {
        (lambda - epsilon * gap).max(0.0)
    }
    /// Update Lagrange multiplier for augmented Lagrangian method.
    ///
    /// lambda_new = max(0, lambda + epsilon * g_n)
    pub fn update_lagrange_multiplier(lambda: f64, gap: f64, epsilon: f64) -> f64 {
        (lambda - epsilon * gap).max(0.0)
    }
}
/// Derjaguin-Muller-Toporov (DMT) adhesive contact model.
///
/// Suitable for stiff materials with low adhesion (complements JKR).
pub struct DmtContact {
    /// Underlying Hertz contact.
    pub hertz: HertzContact,
    /// Work of adhesion W (J/m²).
    pub work_of_adhesion: f64,
}
impl DmtContact {
    /// DMT contact radius: same as Hertz but with modified load.
    ///
    /// a = (3(F + 2πWR*) R* / (4E*))^(1/3)
    pub fn contact_radius(&self, normal_force: f64) -> f64 {
        let r_star = self.hertz.effective_radius();
        let e_star = ElasticSolid::reduced_modulus(&self.hertz.solid_a, &self.hertz.solid_b);
        let f_eff = normal_force + 2.0 * PI * self.work_of_adhesion * r_star;
        if f_eff <= 0.0 {
            return 0.0;
        }
        (3.0 * f_eff * r_star / (4.0 * e_star)).cbrt()
    }
    /// DMT pull-off force.
    ///
    /// F_pull_off = -2 * π * W * R*
    pub fn pull_off_force(&self) -> f64 {
        let r_star = self.hertz.effective_radius();
        -2.0 * PI * self.work_of_adhesion * r_star
    }
    /// Contact area at zero applied load.
    pub fn zero_load_contact_area(&self) -> f64 {
        let a = self.contact_radius(0.0);
        PI * a * a
    }
    /// DMT contact pressure distribution (Hertz-like).
    pub fn pressure_at_zero_force(&self) -> f64 {
        let a = self.contact_radius(0.0);
        if a < 1e-30 {
            return 0.0;
        }
        let f_eff = 2.0 * PI * self.work_of_adhesion * self.hertz.effective_radius();
        3.0 * f_eff / (2.0 * PI * a * a)
    }
}
/// Greenwood–Williamson rough surface contact model.
pub struct RoughSurfaceContact {
    /// Asperity surface density N (1/m^2).
    pub asperity_density: f64,
    /// Mean asperity tip radius beta (m).
    pub asperity_radius: f64,
    /// Standard deviation of asperity heights sigma_s (m).
    pub asperity_height_std: f64,
}
impl RoughSurfaceContact {
    /// Gaussian probability density phi(z) of asperity heights.
    fn phi(&self, z: f64) -> f64 {
        let sigma = self.asperity_height_std;
        (-0.5 * (z / sigma).powi(2)).exp() / (sigma * (2.0 * PI).sqrt())
    }
    /// Complementary-CDF  P(Z > d) using erfc approximation.
    fn prob_contact(&self, separation: f64) -> f64 {
        let s = self.asperity_height_std;
        let t = separation / (s * 2.0_f64.sqrt());
        0.5 * erfc(t)
    }
    /// Mean asperity overlap integral  E\[max(z-d, 0)\] under Gaussian pdf.
    fn mean_overlap(&self, separation: f64) -> f64 {
        let s = self.asperity_height_std;
        let t = separation / (s * 2.0_f64.sqrt());
        s / (2.0 * PI).sqrt() * (-t * t).exp() - separation * 0.5 * erfc(t)
    }
    /// Mean asperity overlap^(3/2) integral  E\[max(z-d,0)^(3/2)\] under Gaussian.
    fn mean_overlap_3_2(&self, separation: f64) -> f64 {
        let s = self.asperity_height_std;
        let n = 200usize;
        let lo = separation;
        let hi = separation + 8.0 * s;
        let dz = (hi - lo) / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let z = lo + (i as f64 + 0.5) * dz;
            let overlap = z - separation;
            sum += overlap.powf(1.5) * self.phi(z) * dz;
        }
        sum
    }
    /// Real contact area (m^2) per unit nominal area for a given separation d.
    pub fn contact_area(&self, separation: f64, _solid: &ElasticSolid) -> f64 {
        PI * self.asperity_density * self.asperity_radius * self.mean_overlap(separation)
    }
    /// Normal force (Pa, i.e. contact pressure) for a given separation d.
    pub fn normal_force(&self, separation: f64, solid: &ElasticSolid) -> f64 {
        let e_star = solid.youngs_modulus / (1.0 - solid.poisson_ratio * solid.poisson_ratio);
        (4.0 / 3.0)
            * self.asperity_density
            * e_star
            * self.asperity_radius.sqrt()
            * self.mean_overlap_3_2(separation)
    }
    /// Fraction of nominal area that is in real contact (dimensionless, 0..1).
    pub fn real_contact_area_fraction(&self, separation: f64) -> f64 {
        PI * self.asperity_radius
            * self.prob_contact(separation)
            * self.asperity_density
            * self.asperity_height_std
    }
    /// Plasticity index: ψ = (E*/H) * sqrt(σ_s / β).
    ///
    /// ψ < 0.6 → elastic contact, ψ > 1.0 → plastic contact.
    pub fn plasticity_index(&self, e_star: f64, hardness: f64) -> f64 {
        if hardness < 1e-30 || self.asperity_radius < 1e-30 {
            return 0.0;
        }
        (e_star / hardness) * (self.asperity_height_std / self.asperity_radius).sqrt()
    }
}
/// Extended rough surface contact with plasticity and adhesion.
pub struct ExtendedGwContact {
    /// Base GW model.
    pub gw: RoughSurfaceContact,
    /// Hardness H (Pa) for plasticity criterion.
    pub hardness: f64,
    /// Work of adhesion for adhesive rough contact.
    pub work_of_adhesion: f64,
}
impl ExtendedGwContact {
    /// Adhesive contact force using GW + DMT adhesion per asperity.
    ///
    /// Total adhesive force = N * A * 2 * π * W * β * prob_contact(d)
    pub fn adhesive_force(&self, separation: f64, nominal_area: f64) -> f64 {
        let prob = self.gw.real_contact_area_fraction(separation);
        2.0 * PI
            * self.work_of_adhesion
            * self.gw.asperity_radius
            * self.gw.asperity_density
            * nominal_area
            * prob
    }
    /// Fraction of asperities in plastic deformation.
    ///
    /// Criterion: contact pressure at asperity tip > hardness/3.
    pub fn plastic_fraction(&self, separation: f64, e_star: f64) -> f64 {
        if self.hardness < 1e-30 {
            return 0.0;
        }
        let p_mean = self.gw.normal_force(
            separation,
            &ElasticSolid {
                youngs_modulus: e_star * 0.9,
                poisson_ratio: 0.3,
            },
        );
        (p_mean / (self.hardness / 3.0)).clamp(0.0, 1.0)
    }
    /// Effective friction coefficient (Bowden-Tabor model).
    ///
    /// mu = shear_strength / (H/3)
    pub fn friction_coefficient(&self, shear_strength: f64) -> f64 {
        if self.hardness < 1e-30 {
            return 0.0;
        }
        shear_strength / (self.hardness / 3.0)
    }
}
/// Adhesion model selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdhesionModel {
    /// Hertz (no adhesion)
    Hertz,
    /// JKR adhesion
    Jkr,
    /// DMT adhesion
    Dmt,
}
/// Johnson-Kendall-Roberts (JKR) adhesive contact model.
///
/// Extends Hertz contact to include surface adhesion via surface energy.
pub struct JkrContact {
    /// Underlying Hertz contact.
    pub hertz: HertzContact,
    /// Work of adhesion W (J/m²): W = gamma_1 + gamma_2 - gamma_12.
    pub work_of_adhesion: f64,
}
impl JkrContact {
    /// JKR contact radius under applied load F.
    ///
    /// a³ = (R*/E*) * \[F + 3πWR* + sqrt(6πWR*F + (3πWR*)²)\]
    pub fn contact_radius(&self, normal_force: f64) -> f64 {
        let r_star = self.hertz.effective_radius();
        let e_star = ElasticSolid::reduced_modulus(&self.hertz.solid_a, &self.hertz.solid_b);
        let w = self.work_of_adhesion;
        let term1 = 3.0 * PI * w * r_star;
        let under_sqrt = 6.0 * PI * w * r_star * normal_force + term1 * term1;
        let a3 = (r_star / e_star) * (normal_force + term1 + under_sqrt.sqrt());
        if a3 <= 0.0 {
            return 0.0;
        }
        a3.cbrt()
    }
    /// Pull-off force (critical detachment force) in the JKR model.
    ///
    /// F_pull_off = -(3/2) * π * W * R*
    pub fn pull_off_force(&self) -> f64 {
        let r_star = self.hertz.effective_radius();
        -1.5 * PI * self.work_of_adhesion * r_star
    }
    /// Contact area at zero applied load (spontaneous adhesion).
    pub fn zero_load_contact_area(&self) -> f64 {
        let a = self.contact_radius(0.0);
        PI * a * a
    }
    /// JKR indentation depth.
    ///
    /// delta = a²/R* - sqrt(2πWa/E*)
    pub fn indentation(&self, normal_force: f64) -> f64 {
        let r_star = self.hertz.effective_radius();
        let e_star = ElasticSolid::reduced_modulus(&self.hertz.solid_a, &self.hertz.solid_b);
        let a = self.contact_radius(normal_force);
        if a <= 0.0 {
            return 0.0;
        }
        a * a / r_star - (2.0 * PI * self.work_of_adhesion * a / e_star).sqrt()
    }
    /// Elastic energy stored in JKR contact.
    pub fn stored_elastic_energy(&self, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        let e_star = ElasticSolid::reduced_modulus(&self.hertz.solid_a, &self.hertz.solid_b);
        let r_star = self.hertz.effective_radius();
        let hertz_energy = 8.0 * e_star * a.powi(5) / (15.0 * r_star);
        let adhesion_energy = PI * self.work_of_adhesion * a * a;
        hertz_energy - adhesion_energy
    }
}
/// Node-to-segment contact element (2-D segment in 3-D space).
pub struct NodeToSegmentContact {
    /// Index of the slave (follower) node.
    pub slave_node: usize,
    /// Indices of the two master segment end-nodes.
    pub master_edge: [usize; 2],
    /// Current signed gap (negative = penetration).
    pub gap: f64,
    /// Unit outward normal of the contact surface.
    pub contact_normal: [f64; 3],
    /// Lagrange multiplier – contact pressure (Pa).
    pub lambda: f64,
}
impl NodeToSegmentContact {
    /// Detect contact between a slave node and a master edge segment.
    pub fn detect_contact(
        slave_pos: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
    ) -> (f64, [f64; 3]) {
        let proj = Self::project_onto_segment(slave_pos, master_a, master_b);
        let d = [
            slave_pos[0] - proj[0],
            slave_pos[1] - proj[1],
            slave_pos[2] - proj[2],
        ];
        let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if dist < f64::EPSILON {
            return (0.0, [0.0, 0.0, 1.0]);
        }
        let normal = [d[0] / dist, d[1] / dist, d[2] / dist];
        let seg = [
            master_b[0] - master_a[0],
            master_b[1] - master_a[1],
            master_b[2] - master_a[2],
        ];
        let seg_len = (seg[0] * seg[0] + seg[1] * seg[1] + seg[2] * seg[2]).sqrt();
        let _ = seg_len;
        (dist, normal)
    }
    /// Orthogonal projection of point p onto segment \[a, b\].
    pub fn project_onto_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let ab_sq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        if ab_sq < f64::EPSILON {
            return a;
        }
        let t = ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / ab_sq).clamp(0.0, 1.0);
        [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]]
    }
    /// Closest distance between two line segments.
    pub fn segment_segment_distance(a1: [f64; 3], a2: [f64; 3], b1: [f64; 3], b2: [f64; 3]) -> f64 {
        let d1 = point_segment_distance(a1, b1, b2);
        let d2 = point_segment_distance(a2, b1, b2);
        let d3 = point_segment_distance(b1, a1, a2);
        let d4 = point_segment_distance(b2, a1, a2);
        d1.min(d2).min(d3).min(d4)
    }
}
/// Hertz contact between two elastic bodies (sphere-on-flat or sphere-on-sphere).
pub struct HertzContact {
    /// Radius of body A in metres (use `f64::INFINITY` for a flat surface).
    pub radius_a: f64,
    /// Radius of body B in metres (use `f64::INFINITY` for a flat surface).
    pub radius_b: f64,
    /// Elastic properties of body A.
    pub solid_a: ElasticSolid,
    /// Elastic properties of body B.
    pub solid_b: ElasticSolid,
}
impl HertzContact {
    /// Effective (combined) radius R*:
    ///   1/R* = 1/R_A + 1/R_B
    pub fn effective_radius(&self) -> f64 {
        let inv = 1.0 / self.radius_a + 1.0 / self.radius_b;
        1.0 / inv
    }
    /// Hertz contact radius a = (3*F*R* / (4*E*))^(1/3).
    pub fn contact_radius(&self, normal_force: f64) -> f64 {
        let r_star = self.effective_radius();
        let e_star = ElasticSolid::reduced_modulus(&self.solid_a, &self.solid_b);
        (3.0 * normal_force * r_star / (4.0 * e_star)).cbrt()
    }
    /// Peak (maximum) Hertz pressure p0 = 3*F / (2*pi*a^2).
    pub fn max_pressure(&self, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        3.0 * normal_force / (2.0 * PI * a * a)
    }
    /// Mean Hertz pressure p_mean = F / (pi * a^2).
    pub fn mean_pressure(&self, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        normal_force / (PI * a * a)
    }
    /// Normal contact stiffness k = 2 * E* * a.
    pub fn contact_stiffness(&self, normal_force: f64) -> f64 {
        let e_star = ElasticSolid::reduced_modulus(&self.solid_a, &self.solid_b);
        let a = self.contact_radius(normal_force);
        2.0 * e_star * a
    }
    /// Hertz indentation depth delta = a^2 / R*.
    pub fn hertz_indentation(&self, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        let r_star = self.effective_radius();
        a * a / r_star
    }
    /// Contact area A = pi * a^2.
    pub fn contact_area(&self, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        PI * a * a
    }
    /// Hertz pressure distribution p(r) = p0 * sqrt(1 - (r/a)^2) for r <= a, else 0.
    pub fn pressure_distribution(&self, r: f64, normal_force: f64) -> f64 {
        let a = self.contact_radius(normal_force);
        if r > a {
            return 0.0;
        }
        let p0 = self.max_pressure(normal_force);
        let ratio = r / a;
        p0 * (1.0 - ratio * ratio).sqrt()
    }
    /// Elastic strain energy stored in Hertz contact.
    ///
    /// U = (2/5) * F * delta
    pub fn stored_energy(&self, normal_force: f64) -> f64 {
        let delta = self.hertz_indentation(normal_force);
        0.4 * normal_force * delta
    }
    /// Force-displacement relationship: F = (4/3) * E* * sqrt(R*) * delta^(3/2).
    pub fn force_from_indentation(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let r_star = self.effective_radius();
        let e_star = ElasticSolid::reduced_modulus(&self.solid_a, &self.solid_b);
        (4.0 / 3.0) * e_star * r_star.sqrt() * delta.powf(1.5)
    }
    /// Compute pressure at multiple radial points.
    pub fn pressure_profile(&self, normal_force: f64, n_points: usize) -> Vec<(f64, f64)> {
        let a = self.contact_radius(normal_force);
        let mut profile = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let r = a * (i as f64) / (n_points as f64 - 1.0).max(1.0);
            let p = self.pressure_distribution(r, normal_force);
            profile.push((r, p));
        }
        profile
    }
}
