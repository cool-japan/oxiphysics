//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::erfc;
use std::f64::consts::PI;

/// Archard's wear model: Q = K * W / H * s (volume loss).
#[derive(Debug, Clone, Copy)]
pub struct ArchardWear {
    /// Wear coefficient K (dimensionless, typically 1e-3 to 1e-8).
    pub k: f64,
    /// Hardness of softer material \[Pa\].
    pub hardness: f64,
}
impl ArchardWear {
    /// Create an Archard wear model.
    pub fn new(k: f64, hardness: f64) -> Self {
        Self { k, hardness }
    }
    /// Compute wear rate \[m³/m\] = K * P / H.
    /// `pressure` is normal contact pressure \[Pa\].
    pub fn wear_rate(&self, pressure: f64) -> f64 {
        self.k * pressure / self.hardness
    }
    /// Compute volume loss \[m³\] over sliding distance `s` \[m\] under normal load `w` \[N\].
    pub fn volume_loss(&self, w: f64, s: f64) -> f64 {
        self.k * w * s / self.hardness
    }
    /// Compute linear wear depth \[m\] over area `a` \[m²\] sliding distance `s`.
    pub fn depth_loss(&self, w: f64, s: f64, area: f64) -> f64 {
        if area < 1e-30 {
            return 0.0;
        }
        self.volume_loss(w, s) / area
    }
}
/// Scuffing risk assessment based on contact flash temperature criterion (Blok).
#[derive(Debug, Clone, Copy)]
pub struct ScuffingModel {
    /// Friction coefficient μ.
    pub mu: f64,
    /// Thermal conductivity body 1 k1 \[W/(m·K)\].
    pub k1: f64,
    /// Thermal conductivity body 2 k2 \[W/(m·K)\].
    pub k2: f64,
    /// Density * specific heat product ρcp1 \[J/(m³·K)\].
    pub rho_cp1: f64,
    /// Density * specific heat product ρcp2 \[J/(m³·K)\].
    pub rho_cp2: f64,
    /// Critical scuffing temperature Tc \[K\].
    pub t_critical: f64,
    /// Bulk temperature T_bulk \[K\].
    pub t_bulk: f64,
}
impl ScuffingModel {
    /// Create a typical steel-on-steel scuffing model.
    pub fn steel_steel() -> Self {
        ScuffingModel {
            mu: 0.1,
            k1: 50.0,
            k2: 50.0,
            rho_cp1: 3.75e6,
            rho_cp2: 3.75e6,
            t_critical: 700.0 + 273.15,
            t_bulk: 80.0 + 273.15,
        }
    }
    /// Thermal diffusivity κ \[m²/s\] for body 1.
    pub fn thermal_diffusivity_1(&self) -> f64 {
        self.k1 / self.rho_cp1
    }
    /// Thermal diffusivity κ for body 2.
    pub fn thermal_diffusivity_2(&self) -> f64 {
        self.k2 / self.rho_cp2
    }
    /// Flash temperature rise ΔTf \[K\] using Blok's contact temperature formula.
    /// `q` = heat flux at contact \[W/m²\], `a` = contact radius \[m\], `v` = velocity \[m/s\].
    pub fn flash_temperature_rise(&self, q: f64, a: f64, v: f64) -> f64 {
        if a < 1e-30 || v < 1e-10 {
            return 0.0;
        }
        let beta1 = (self.rho_cp1 * self.k1).sqrt();
        let beta2 = (self.rho_cp2 * self.k2).sqrt();
        1.11 * q * (a / v).sqrt() / (beta1 + beta2)
    }
    /// Contact temperature Tc = T_bulk + ΔTf.
    pub fn contact_temperature(&self, q: f64, a: f64, v: f64) -> f64 {
        self.t_bulk + self.flash_temperature_rise(q, a, v)
    }
    /// Scuffing risk: returns true if contact temperature exceeds critical.
    pub fn scuffing_risk(&self, q: f64, a: f64, v: f64) -> bool {
        self.contact_temperature(q, a, v) > self.t_critical
    }
    /// Safety factor against scuffing (Tc_critical / T_contact).
    pub fn safety_factor(&self, q: f64, a: f64, v: f64) -> f64 {
        let tc = self.contact_temperature(q, a, v);
        if tc < 1.0 {
            return f64::INFINITY;
        }
        self.t_critical / tc
    }
}
/// Rolling element bearing life analysis per ISO 281.
#[derive(Debug, Clone, Copy)]
pub struct BearingLife {
    /// Basic dynamic load rating C \[N\].
    pub c: f64,
    /// Applied equivalent dynamic load P \[N\].
    pub p: f64,
    /// Life exponent p (3 for ball bearings, 10/3 for roller bearings).
    pub life_exponent: f64,
    /// Lubrication factor aISO (1.0 = reference).
    pub a_iso: f64,
    /// Material/contamination factor (0.1–4.0).
    pub a_skf: f64,
}
impl BearingLife {
    /// Create a ball bearing life model.
    pub fn ball_bearing(c: f64, p: f64) -> Self {
        BearingLife {
            c,
            p,
            life_exponent: 3.0,
            a_iso: 1.0,
            a_skf: 1.0,
        }
    }
    /// Create a roller bearing life model.
    pub fn roller_bearing(c: f64, p: f64) -> Self {
        BearingLife {
            c,
            p,
            life_exponent: 10.0 / 3.0,
            a_iso: 1.0,
            a_skf: 1.0,
        }
    }
    /// Basic L10 life in millions of revolutions.
    pub fn l10_mrv(&self) -> f64 {
        if self.p < 1e-10 {
            return f64::INFINITY;
        }
        (self.c / self.p).powf(self.life_exponent)
    }
    /// L10 life in hours given shaft speed n \[rpm\].
    pub fn l10_hours(&self, n_rpm: f64) -> f64 {
        if n_rpm < 1e-10 {
            return f64::INFINITY;
        }
        self.l10_mrv() * 1e6 / (60.0 * n_rpm)
    }
    /// Modified L10m life (SKF): Lnm = a1 * a_ISO * a_SKF * L10.
    /// `a1` = reliability factor (1.0 for 90% reliability).
    pub fn l10m_hours(&self, n_rpm: f64, a1: f64) -> f64 {
        a1 * self.a_iso * self.a_skf * self.l10_hours(n_rpm)
    }
    /// Weibull reliability at life t \[hours\], given L10 (90% reliability).
    /// R(t) = exp\[-(t/L50)^e\] where e ≈ 1.5 for bearings.
    pub fn reliability_at(&self, t_hours: f64, n_rpm: f64) -> f64 {
        let l10 = self.l10_hours(n_rpm);
        if l10 < 1e-10 {
            return 0.0;
        }
        let l50 = 5.0 * l10;
        let weibull_slope = 1.5_f64;
        (-(t_hours / l50).powf(weibull_slope)).exp()
    }
    /// Equivalent dynamic load for combined radial Fr and axial Fa loads.
    /// P = X * Fr + Y * Fa (simplified, X=1, Y=0 for pure radial).
    pub fn equivalent_load(fr: f64, fa: f64, x: f64, y: f64) -> f64 {
        x * fr + y * fa
    }
    /// Speed factor fn = (33.3/n)^(1/3) \[rpm\].
    pub fn speed_factor(&self, n_rpm: f64) -> f64 {
        if n_rpm < 1e-10 {
            return 0.0;
        }
        (33.3 / n_rpm).powf(1.0 / 3.0)
    }
}
/// Results of a Hertz contact calculation.
#[derive(Debug, Clone, Copy)]
pub struct HertzResult {
    /// Contact radius \[m\].
    pub contact_radius: f64,
    /// Maximum contact pressure \[Pa\].
    pub max_pressure: f64,
    /// Contact depth (approach) \[m\].
    pub depth: f64,
    /// Normal force \[N\].
    pub force: f64,
}
/// Thermal wear model including frictional heating and oxidative wear.
#[derive(Debug, Clone, Copy)]
pub struct ThermalWear {
    /// Friction coefficient.
    pub friction_coeff: f64,
    /// Thermal conductivity of material 1 \[W/(m·K)\].
    pub k1: f64,
    /// Thermal conductivity of material 2 \[W/(m·K)\].
    pub k2: f64,
    /// Ambient temperature \[K\].
    pub t_ambient: f64,
    /// Oxidative wear activation energy Q \[J/mol\].
    pub activation_energy: f64,
    /// Pre-exponential wear factor A₀.
    pub wear_factor_a0: f64,
}
impl ThermalWear {
    /// Create a thermal wear model.
    pub fn new(
        friction_coeff: f64,
        k1: f64,
        k2: f64,
        t_ambient: f64,
        activation_energy: f64,
        wear_factor_a0: f64,
    ) -> Self {
        Self {
            friction_coeff,
            k1,
            k2,
            t_ambient,
            activation_energy,
            wear_factor_a0,
        }
    }
    /// Flash temperature rise at contact (Blok's equation), simplified.
    /// `f` = friction force \[N\], `v` = sliding speed \[m/s\], `a` = contact radius \[m\].
    pub fn flash_temperature(&self, f: f64, v: f64, a: f64) -> f64 {
        if a < 1e-30 {
            return self.t_ambient;
        }
        let heat_rate = self.friction_coeff * f * v;
        let delta_t = heat_rate / (PI * a * (self.k1 + self.k2));
        self.t_ambient + delta_t
    }
    /// Oxidative wear rate at temperature T (Arrhenius-type).
    /// Returns wear rate \[m/s\].
    pub fn oxidative_wear_rate(&self, temperature: f64) -> f64 {
        const R_GAS: f64 = 8.314;
        self.wear_factor_a0 * (-(self.activation_energy / (R_GAS * temperature))).exp()
    }
    /// Total wear volume \[m³\] over time `t` \[s\] at contact conditions.
    pub fn total_wear_volume(&self, f: f64, v: f64, a: f64, t: f64) -> f64 {
        let temp = self.flash_temperature(f, v, a);
        let rate = self.oxidative_wear_rate(temp);
        rate * a * a * t
    }
}
/// Two-body and three-body abrasive wear model.
#[derive(Debug, Clone, Copy)]
pub struct AbrasiveWear {
    /// Abrasion wear coefficient Ka (dimensionless).
    pub ka: f64,
    /// Hardness of abraded material \[Pa\].
    pub hardness: f64,
    /// Abrasive particle average size dp \[m\].
    pub particle_size: f64,
    /// Fraction of particles in contact (efficiency factor, 0–1).
    pub efficiency: f64,
}
impl AbrasiveWear {
    /// Create an abrasive wear model.
    pub fn new(ka: f64, hardness: f64, particle_size: f64, efficiency: f64) -> Self {
        AbrasiveWear {
            ka,
            hardness,
            particle_size,
            efficiency,
        }
    }
    /// Two-body abrasive wear volume \[m³\] over distance s \[m\] under load W \[N\].
    pub fn two_body_volume(&self, w: f64, s: f64) -> f64 {
        self.ka * w * s / self.hardness
    }
    /// Three-body abrasive wear volume \[m³\] (factor ~3× lower than two-body).
    pub fn three_body_volume(&self, w: f64, s: f64) -> f64 {
        self.two_body_volume(w, s) / 3.0
    }
    /// Scratch depth per asperity pass \[m\] (Bowden-Tabor model).
    pub fn scratch_depth(&self, w_per_asperity: f64) -> f64 {
        (2.0 * w_per_asperity / (PI * self.hardness)).sqrt()
    }
    /// Volume per scratch \[m³\] for scratch length L \[m\].
    pub fn volume_per_scratch(&self, w_per_asperity: f64, scratch_length: f64) -> f64 {
        let depth = self.scratch_depth(w_per_asperity);
        0.5 * PI * depth.powi(2) * scratch_length
    }
    /// Specific wear rate \[m³/(N·m)\] = Ka/H.
    pub fn specific_wear_rate(&self) -> f64 {
        self.ka / self.hardness
    }
}
/// Elliptical Hertz contact parameters for general curved bodies.
/// Based on Hertz contact theory for two bodies with different principal radii.
#[derive(Debug, Clone, Copy)]
pub struct HertzElliptical {
    /// Effective modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Sum of principal curvatures (1/R1a + 1/R1b + 1/R2a + 1/R2b) \[1/m\].
    pub sum_curvature: f64,
    /// Difference parameter for ellipticity.
    pub curvature_difference: f64,
}
impl HertzElliptical {
    /// Create elliptical contact from two bodies' principal radii.
    /// R1a, R1b = principal radii of body 1; R2a, R2b = principal radii of body 2.
    pub fn new(
        e1: f64,
        nu1: f64,
        e2: f64,
        nu2: f64,
        r1a: f64,
        r1b: f64,
        r2a: f64,
        r2b: f64,
    ) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        let sum_k = 1.0 / r1a + 1.0 / r1b + 1.0 / r2a + 1.0 / r2b;
        let diff_k = (1.0 / r1a - 1.0 / r1b + 1.0 / r2a - 1.0 / r2b)
            .powi(2)
            .sqrt();
        HertzElliptical {
            effective_modulus: e_star,
            sum_curvature: sum_k,
            curvature_difference: diff_k,
        }
    }
    /// Ellipticity parameter k (semi-axis ratio b/a, ≤ 1).
    /// Uses Hamrock approximation: k = (B/A)^(2/3) where A = (Σ-δ)/2, B = (Σ+δ)/2.
    pub fn ellipticity_ratio(&self) -> f64 {
        let cap_b = (self.sum_curvature + self.curvature_difference) / 2.0;
        let cap_a = (self.sum_curvature - self.curvature_difference) / 2.0;
        if cap_a < 1e-30 {
            return 1.0;
        }
        (cap_b / cap_a).powf(2.0 / 3.0).min(1.0)
    }
    /// Contact semi-axis a \[m\] (major) under load F \[N\].
    pub fn semi_axis_a(&self, force: f64) -> f64 {
        let r_eff = 2.0 / self.sum_curvature;
        let h = HertzContact {
            effective_modulus: self.effective_modulus,
            effective_radius: r_eff,
        };
        h.compute(force).contact_radius
    }
    /// Contact semi-axis b \[m\] (minor) = k * a.
    pub fn semi_axis_b(&self, force: f64) -> f64 {
        self.ellipticity_ratio() * self.semi_axis_a(force)
    }
    /// Maximum Hertz pressure p0 \[Pa\] for elliptical contact.
    pub fn max_pressure(&self, force: f64) -> f64 {
        let a = self.semi_axis_a(force);
        let b = self.semi_axis_b(force);
        if a < 1e-30 || b < 1e-30 {
            return 0.0;
        }
        3.0 * force / (2.0 * PI * a * b)
    }
    /// Contact area (ellipse) \[m²\].
    pub fn contact_area(&self, force: f64) -> f64 {
        PI * self.semi_axis_a(force) * self.semi_axis_b(force)
    }
}
/// Gear tooth tribology: pitch line velocity, sliding ratio, contact stress.
#[derive(Debug, Clone)]
pub struct GearTooth {
    /// Module m \[mm\].
    pub module: f64,
    /// Number of teeth on gear 1.
    pub z1: u32,
    /// Number of teeth on gear 2.
    pub z2: u32,
    /// Pressure angle φ \[degrees\].
    pub pressure_angle: f64,
    /// Face width b \[mm\].
    pub face_width: f64,
    /// Transmitted tangential force Wt \[N\].
    pub wt: f64,
    /// Effective elastic coefficient ZE \[(N/mm²)^0.5\].
    pub ze: f64,
    /// Speed of gear 1 n1 \[rpm\].
    pub n1: f64,
}
impl GearTooth {
    /// Create a spur gear pair.
    pub fn spur_pair(
        module: f64,
        z1: u32,
        z2: u32,
        pressure_angle: f64,
        face_width: f64,
        wt: f64,
        n1: f64,
    ) -> Self {
        GearTooth {
            module,
            z1,
            z2,
            pressure_angle,
            face_width,
            wt,
            ze: 191.0,
            n1,
        }
    }
    /// Pitch diameter of gear 1 d1 \[mm\].
    pub fn d1(&self) -> f64 {
        self.module * self.z1 as f64
    }
    /// Pitch diameter of gear 2 d2 \[mm\].
    pub fn d2(&self) -> f64 {
        self.module * self.z2 as f64
    }
    /// Pitch line velocity Vp \[m/s\].
    pub fn pitch_line_velocity(&self) -> f64 {
        PI * self.d1() * self.n1 / (60.0 * 1000.0)
    }
    /// Gear ratio i = z2/z1.
    pub fn gear_ratio(&self) -> f64 {
        self.z2 as f64 / self.z1 as f64
    }
    /// Contact ratio εα (profile contact ratio) for spur gear.
    pub fn contact_ratio(&self) -> f64 {
        let phi = self.pressure_angle.to_radians();
        let r1 = self.d1() / 2.0;
        let r2 = self.d2() / 2.0;
        let ra1 = r1 + self.module;
        let ra2 = r2 + self.module;
        let cp = (self.d1() + self.d2()) / 2.0 * phi.sin();
        let term1 = (ra1.powi(2) - (r1 * phi.cos()).powi(2)).sqrt();
        let term2 = (ra2.powi(2) - (r2 * phi.cos()).powi(2)).sqrt();
        (term1 + term2 - cp) / (PI * self.module * phi.cos())
    }
    /// Hertz contact stress at pitch point σH \[N/mm²\] (AGMA/ISO simplified).
    pub fn hertz_contact_stress(&self) -> f64 {
        let phi = self.pressure_angle.to_radians();
        let d1 = self.d1();
        let i = self.gear_ratio();
        let term = self.wt / (d1 * self.face_width) * (i + 1.0) / i / phi.cos();
        if term < 0.0 {
            return 0.0;
        }
        self.ze * term.sqrt()
    }
    /// Sliding velocity at pitch point offset s \[mm\] from pitch point.
    pub fn sliding_velocity_at(&self, s_mm: f64) -> f64 {
        let _vp = self.pitch_line_velocity() * 1000.0;
        let omega1 = 2.0 * PI * self.n1 / 60.0;
        let i = self.gear_ratio();
        omega1 * s_mm * (1.0 + 1.0 / i) / 1000.0
    }
    /// Specific sliding ratio at arbitrary contact point (simplified AGMA).
    pub fn sliding_ratio(&self) -> f64 {
        let phi = self.pressure_angle.to_radians();
        2.0 * PI * (1.0 / self.z1 as f64 + 1.0 / self.z2 as f64) / (PI * phi.cos())
            * self.contact_ratio()
    }
    /// Power loss from friction \[W\], given friction coefficient μ.
    pub fn friction_power_loss(&self, mu: f64) -> f64 {
        let vs = self.sliding_velocity_at(self.module);
        mu * self.wt * vs.abs()
    }
}
/// Brake pad wear and thermal analysis for disc brakes.
#[derive(Debug, Clone)]
pub struct BrakePad {
    /// Normal force on pad Fn \[N\].
    pub fn_force: f64,
    /// Friction coefficient μ.
    pub mu: f64,
    /// Rotor radius (effective rubbing radius) r \[m\].
    pub rotor_radius: f64,
    /// Pad area A \[m²\].
    pub pad_area: f64,
    /// Thermal conductivity of pad kp \[W/(m·K)\].
    pub k_pad: f64,
    /// Heat partition factor for pad βp (fraction of heat to pad, 0–1).
    pub beta_pad: f64,
    /// Pad wear coefficient Kw \[m²/N\].
    pub kw: f64,
    /// Pad thickness tp \[m\].
    pub tp: f64,
    /// Specific heat of pad cp \[J/(kg·K)\].
    pub cp_pad: f64,
    /// Density of pad \[kg/m³\].
    pub density_pad: f64,
}
impl BrakePad {
    /// Create a typical semi-metallic brake pad.
    pub fn semi_metallic() -> Self {
        BrakePad {
            fn_force: 3000.0,
            mu: 0.40,
            rotor_radius: 0.12,
            pad_area: 50e-4,
            k_pad: 2.0,
            beta_pad: 0.05,
            kw: 1e-14,
            tp: 0.012,
            cp_pad: 900.0,
            density_pad: 2500.0,
        }
    }
    /// Braking torque Tb \[N·m\].
    pub fn braking_torque(&self) -> f64 {
        self.mu * self.fn_force * self.rotor_radius
    }
    /// Frictional power dissipation Q \[W\] at sliding speed v \[m/s\].
    pub fn frictional_power(&self, v: f64) -> f64 {
        self.mu * self.fn_force * v
    }
    /// Average interface temperature rise ΔT \[K\] (steady-state).
    pub fn temperature_rise_steady(&self, v: f64) -> f64 {
        let q = self.frictional_power(v) * self.beta_pad;
        if self.k_pad < 1e-10 || self.pad_area < 1e-10 {
            return 0.0;
        }
        q * self.tp / (self.k_pad * self.pad_area)
    }
    /// Transient temperature rise \[K\] after time t \[s\], simplified lumped capacity.
    pub fn temperature_rise_transient(&self, v: f64, t: f64) -> f64 {
        let q = self.frictional_power(v) * self.beta_pad;
        let mass = self.density_pad * self.pad_area * self.tp;
        if mass < 1e-10 {
            return 0.0;
        }
        q * t / (mass * self.cp_pad)
    }
    /// Pad wear depth \[m\] over sliding distance s \[m\].
    pub fn wear_depth(&self, s: f64) -> f64 {
        self.kw * self.fn_force * s / self.pad_area
    }
    /// Remaining pad life \[m of sliding distance\] from initial thickness tp.
    pub fn remaining_life(&self, current_thickness: f64, min_thickness: f64) -> f64 {
        let wear_allowed = current_thickness - min_thickness;
        if wear_allowed <= 0.0 {
            return 0.0;
        }
        wear_allowed * self.pad_area / (self.kw * self.fn_force)
    }
    /// Contact pressure p \[Pa\] = Fn / A.
    pub fn contact_pressure(&self) -> f64 {
        self.fn_force / self.pad_area
    }
}
/// Lubrication regime based on the Stribeck curve and lambda ratio.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LubricationRegimeKind {
    /// Boundary lubrication (λ < 1).
    Boundary,
    /// Mixed lubrication (1 ≤ λ < 3).
    Mixed,
    /// Elastohydrodynamic lubrication (3 ≤ λ < 10).
    Elastohydrodynamic,
    /// Full-film hydrodynamic lubrication (λ ≥ 10).
    FullFilm,
}
/// Surface roughness parameters from ISO 4287 / ASME B46.1.
#[derive(Debug, Clone)]
pub struct SurfaceRoughness {
    /// Arithmetic mean deviation Ra \[m\].
    pub ra: f64,
    /// Root-mean-square roughness Rq \[m\].
    pub rq: f64,
    /// Maximum height of profile Rz \[m\].
    pub rz: f64,
    /// Skewness Rsk (dimensionless).
    pub rsk: f64,
    /// Kurtosis Rku (dimensionless; 3.0 for Gaussian).
    pub rku: f64,
    /// Mean spacing of profile irregularities RSm \[m\].
    pub rsm: f64,
}
impl SurfaceRoughness {
    /// Create a surface roughness profile.
    pub fn new(ra: f64, rq: f64, rz: f64, rsk: f64, rku: f64, rsm: f64) -> Self {
        SurfaceRoughness {
            ra,
            rq,
            rz,
            rsk,
            rku,
            rsm,
        }
    }
    /// Create a ground steel surface (typical Ra = 0.8 µm).
    pub fn ground_steel() -> Self {
        SurfaceRoughness {
            ra: 0.8e-6,
            rq: 1.0e-6,
            rz: 4.0e-6,
            rsk: 0.0,
            rku: 3.0,
            rsm: 80.0e-6,
        }
    }
    /// Create a lapped/polished surface (Ra ≈ 0.1 µm).
    pub fn lapped() -> Self {
        SurfaceRoughness {
            ra: 0.1e-6,
            rq: 0.12e-6,
            rz: 0.6e-6,
            rsk: -0.2,
            rku: 3.5,
            rsm: 20.0e-6,
        }
    }
    /// Combined roughness of two surfaces (root-sum-square).
    pub fn combined_rq(&self, other: &SurfaceRoughness) -> f64 {
        (self.rq.powi(2) + other.rq.powi(2)).sqrt()
    }
    /// Rq / Ra ratio (should be ~1.25 for Gaussian surface).
    pub fn rq_ra_ratio(&self) -> f64 {
        if self.ra < 1e-30 {
            return 0.0;
        }
        self.rq / self.ra
    }
    /// Surface texture index STI = Rz / Ra (bearing capacity quality indicator).
    pub fn surface_texture_index(&self) -> f64 {
        if self.ra < 1e-30 {
            return 0.0;
        }
        self.rz / self.ra
    }
    /// Abbott-Firestone bearing area curve at normalised depth z (0–1).
    /// Returns fraction of surface in contact (0–1).
    pub fn bearing_area(&self, z_norm: f64) -> f64 {
        let threshold = self.rz * (1.0 - z_norm) - self.ra;
        let z = threshold / (self.rq * 2.0_f64.sqrt());
        0.5 * erfc(z)
    }
    /// Valley depth Rv (estimated as Rz * Rsk factor).
    pub fn valley_depth(&self) -> f64 {
        self.rz * if self.rsk < 0.0 { 0.6 } else { 0.4 }
    }
}
/// Run-in friction evolution: exponential decrease from initial to steady-state friction.
#[derive(Debug, Clone, Copy)]
pub struct FrictionEvolution {
    /// Initial friction coefficient (start of run-in).
    pub mu_initial: f64,
    /// Steady-state friction coefficient.
    pub mu_steady: f64,
    /// Run-in time constant τ \[s or m\].
    pub time_constant: f64,
    /// Current time/distance state.
    pub time: f64,
}
impl FrictionEvolution {
    /// Create a friction evolution model.
    pub fn new(mu_initial: f64, mu_steady: f64, time_constant: f64) -> Self {
        Self {
            mu_initial,
            mu_steady,
            time_constant,
            time: 0.0,
        }
    }
    /// Get friction coefficient at current time `t`.
    pub fn friction_at(&self, t: f64) -> f64 {
        self.mu_steady + (self.mu_initial - self.mu_steady) * (-(t / self.time_constant)).exp()
    }
    /// Step simulation by `dt`.
    pub fn step(&mut self, dt: f64) {
        self.time += dt;
    }
    /// Current friction coefficient.
    pub fn current_friction(&self) -> f64 {
        self.friction_at(self.time)
    }
    /// Time to reach within `epsilon` of steady state.
    pub fn convergence_time(&self, epsilon: f64) -> f64 {
        let range = (self.mu_initial - self.mu_steady).abs();
        if range < 1e-12 {
            return 0.0;
        }
        -self.time_constant * (epsilon / range).ln()
    }
}
/// Lubricant viscosity models including Reynolds and Walther equations.
#[derive(Debug, Clone, Copy)]
pub struct LubricantViscosity {
    /// Dynamic viscosity at reference temperature T0 \[K\]: η0 \[Pa·s\].
    pub eta0: f64,
    /// Reference temperature T0 \[K\].
    pub t0: f64,
    /// Reynolds temperature coefficient β \[1/K\].
    pub beta: f64,
    /// Walther equation constants: log(log(ν+c)) = A - B*log(T).
    pub walther_a: f64,
    /// Walther B constant.
    pub walther_b: f64,
    /// Density at reference temperature ρ0 \[kg/m³\].
    pub density: f64,
}
impl LubricantViscosity {
    /// Create a viscosity model for SAE 10W-40 motor oil (approximate).
    pub fn sae10w40() -> Self {
        LubricantViscosity {
            eta0: 0.1,
            t0: 313.15,
            beta: 0.028,
            walther_a: 10.5,
            walther_b: 3.5,
            density: 870.0,
        }
    }
    /// Create a model for ISO VG 46 hydraulic oil.
    pub fn iso_vg46() -> Self {
        LubricantViscosity {
            eta0: 0.046,
            t0: 313.15,
            beta: 0.025,
            walther_a: 9.8,
            walther_b: 3.3,
            density: 875.0,
        }
    }
    /// Dynamic viscosity at temperature T \[K\] using Reynolds equation.
    pub fn eta_reynolds(&self, t: f64) -> f64 {
        self.eta0 * (-(self.beta * (t - self.t0))).exp()
    }
    /// Kinematic viscosity \[m²/s\] at temperature T \[K\].
    pub fn nu_kinematic(&self, t: f64) -> f64 {
        self.eta_reynolds(t) / self.density
    }
    /// Barus equation for pressure-dependent viscosity η(p) \[Pa·s\].
    pub fn eta_barus(&self, t: f64, alpha_p: f64, pressure: f64) -> f64 {
        self.eta_reynolds(t) * (alpha_p * pressure).exp()
    }
    /// Viscosity index (VI) — higher = less temperature sensitivity.
    /// Simplified: VI = (η40 - η100) / (L - H) where L,H are reference viscosities.
    pub fn viscosity_index(&self) -> f64 {
        let eta40 = self.eta_reynolds(313.15);
        let eta100 = self.eta_reynolds(373.15);
        if eta100 < 1e-10 {
            return 0.0;
        }
        100.0 * (eta40 - eta100) / eta100
    }
    /// Pour point estimate \[K\] below which oil loses flowability.
    pub fn pour_point_estimate(&self) -> f64 {
        let target_eta = 3.0;
        if self.eta0 <= target_eta {
            return self.t0;
        }
        self.t0 - (target_eta / self.eta0).ln() / self.beta
    }
}
/// Fretting fatigue model: stick-slip regimes, fretting loop, fretting fatigue life.
#[derive(Debug, Clone)]
pub struct FrettingFatigue {
    /// Normal load Q \[N\] on contact.
    pub normal_load: f64,
    /// Friction coefficient μ (Coulomb).
    pub mu: f64,
    /// Tangential amplitude δ \[m\].
    pub amplitude: f64,
    /// Contact stiffness kt \[N/m\].
    pub contact_stiffness: f64,
    /// Fatigue strength of material σ_fat \[Pa\].
    pub fatigue_strength: f64,
    /// Paris law exponent m for fretting crack propagation.
    pub paris_m: f64,
    /// Paris law coefficient C \[m/cycle / (MPa√m)^m\].
    pub paris_c: f64,
    /// Initial crack half-length a0 \[m\].
    pub a0: f64,
}
impl FrettingFatigue {
    /// Create a fretting fatigue model.
    pub fn new(
        normal_load: f64,
        mu: f64,
        amplitude: f64,
        contact_stiffness: f64,
        fatigue_strength: f64,
    ) -> Self {
        FrettingFatigue {
            normal_load,
            mu,
            amplitude,
            contact_stiffness,
            fatigue_strength,
            paris_m: 3.0,
            paris_c: 1e-11,
            a0: 10e-6,
        }
    }
    /// Maximum tangential (friction) force in the fretting loop Qmax \[N\].
    pub fn max_tangential_force(&self) -> f64 {
        self.mu * self.normal_load
    }
    /// Stick-slip threshold displacement δs \[m\].
    pub fn stick_slip_threshold(&self) -> f64 {
        self.max_tangential_force() / self.contact_stiffness
    }
    /// Determine if gross slip occurs (δ > δs).
    pub fn is_gross_slip(&self) -> bool {
        self.amplitude > self.stick_slip_threshold()
    }
    /// Fretting loop energy dissipated per cycle W \[J\].
    pub fn energy_per_cycle(&self) -> f64 {
        if self.is_gross_slip() {
            let delta_s = self.stick_slip_threshold();
            4.0 * self.max_tangential_force() * (self.amplitude - delta_s)
        } else {
            let ratio = self.amplitude / self.stick_slip_threshold();
            let w_stick = 4.0 * self.contact_stiffness * self.amplitude.powi(2);
            w_stick * (1.0 - (1.0 - ratio).powi(3)) / 3.0
        }
    }
    /// Stress intensity factor range ΔK \[Pa√m\] for crack size a \[m\].
    pub fn stress_intensity_range(&self, stress_range: f64, a: f64) -> f64 {
        stress_range * (PI * a).sqrt()
    }
    /// Cycles to grow crack from a0 to ac \[m\] (Paris law integration).
    pub fn cycles_to_crack(&self, stress_range: f64, a_critical: f64) -> f64 {
        let m = self.paris_m;
        let c = self.paris_c;
        let f_sig = stress_range * PI.sqrt();
        let exponent = 1.0 - m / 2.0;
        if exponent.abs() < 1e-10 {
            return (a_critical / self.a0).ln() / (c * f_sig.powi(2));
        }
        let denom = c * f_sig.powf(m) * exponent;
        if denom.abs() < 1e-30 {
            return f64::INFINITY;
        }
        (a_critical.powf(exponent) - self.a0.powf(exponent)) / denom
    }
    /// Fretting fatigue factor (stress concentration effect), simplified.
    pub fn fretting_factor(&self) -> f64 {
        let q_max = self.max_tangential_force();
        1.0 + q_max / self.fatigue_strength.max(1.0)
    }
}
/// Hertz contact mechanics for elastic bodies.
#[derive(Debug, Clone, Copy)]
pub struct HertzContact {
    /// Effective modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R* \[m\].
    pub effective_radius: f64,
}
impl HertzContact {
    /// Create from two sphere moduli, Poisson ratios, and radii.
    /// Computes E* = 1 / ((1-ν₁²)/E₁ + (1-ν₂²)/E₂) and R* = R₁*R₂/(R₁+R₂).
    pub fn sphere_sphere(e1: f64, nu1: f64, r1: f64, e2: f64, nu2: f64, r2: f64) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        let r_star = r1 * r2 / (r1 + r2);
        Self {
            effective_modulus: e_star,
            effective_radius: r_star,
        }
    }
    /// Create for sphere on a flat surface (R₂ → ∞).
    pub fn sphere_flat(e1: f64, nu1: f64, r: f64, e2: f64, nu2: f64) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        Self {
            effective_modulus: e_star,
            effective_radius: r,
        }
    }
    /// Create for two parallel cylinders of length L (per unit length analysis).
    /// Uses 2D Hertz: a = sqrt(4*F*R*/(π*E*)).
    pub fn cylinder_cylinder(e1: f64, nu1: f64, r1: f64, e2: f64, nu2: f64, r2: f64) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        let r_star = r1 * r2 / (r1 + r2);
        Self {
            effective_modulus: e_star,
            effective_radius: r_star,
        }
    }
    /// Compute Hertz contact parameters for a given normal force \[N\].
    pub fn compute(&self, force: f64) -> HertzResult {
        let a = ((3.0 * force * self.effective_radius) / (4.0 * self.effective_modulus))
            .powf(1.0 / 3.0);
        let p0 = if a > 1e-30 {
            3.0 * force / (2.0 * PI * a * a)
        } else {
            0.0
        };
        let depth = a * a / self.effective_radius;
        HertzResult {
            contact_radius: a,
            max_pressure: p0,
            depth,
            force,
        }
    }
    /// Compute required force for a given contact radius.
    pub fn force_for_radius(&self, a: f64) -> f64 {
        4.0 * self.effective_modulus * a * a * a / (3.0 * self.effective_radius)
    }
}
/// Greenwood-Williamson multi-asperity contact model.
#[derive(Debug, Clone)]
pub struct AsperityModel {
    /// Asperity density \[1/m²\].
    pub asperity_density: f64,
    /// RMS surface roughness σ \[m\].
    pub rms_roughness: f64,
    /// Asperity tip radius β \[m\].
    pub tip_radius: f64,
    /// Effective elastic modulus \[Pa\].
    pub effective_modulus: f64,
}
impl AsperityModel {
    /// Create a Greenwood-Williamson model.
    pub fn new(
        asperity_density: f64,
        rms_roughness: f64,
        tip_radius: f64,
        effective_modulus: f64,
    ) -> Self {
        Self {
            asperity_density,
            rms_roughness,
            tip_radius,
            effective_modulus,
        }
    }
    /// Plasticity index: Ψ = (E*/H) * sqrt(σ/β).
    pub fn plasticity_index(&self, hardness: f64) -> f64 {
        (self.effective_modulus / hardness) * (self.rms_roughness / self.tip_radius).sqrt()
    }
    /// Expected number of asperities in contact for separation `d` \[m\].
    /// Uses Gaussian distribution of asperity heights.
    pub fn n_contacts(&self, area: f64, separation: f64) -> f64 {
        let z = separation / (2.0_f64.sqrt() * self.rms_roughness);
        let p = 0.5 * erfc(z);
        self.asperity_density * area * p
    }
    /// Bearing area (ratio of real to apparent contact area).
    pub fn bearing_area(&self, separation: f64) -> f64 {
        let z = separation / (2.0_f64.sqrt() * self.rms_roughness);
        0.5 * erfc(z)
    }
    /// Total normal load for given separation d over nominal area A.
    pub fn normal_load(&self, area: f64, separation: f64) -> f64 {
        let n = self.n_contacts(area, separation);
        let avg_deform = self.rms_roughness;
        let f_per =
            (4.0 / 3.0) * self.effective_modulus * self.tip_radius.sqrt() * avg_deform.powf(1.5);
        n * f_per
    }
}
/// Extended Stribeck curve with boundary, mixed, and hydrodynamic regimes.
#[derive(Debug, Clone)]
pub struct StribeckCurve {
    /// Boundary friction coefficient μb.
    pub mu_boundary: f64,
    /// Hydrodynamic (full film) friction coefficient μf.
    pub mu_fullfilm: f64,
    /// Stribeck constant Sc \[Pa·s · m/s / Pa\] = \[m\].
    pub sc: f64,
    /// Exponent controlling sharpness of the Stribeck minimum (typically 0.5–2.0).
    pub exponent: f64,
    /// Minimum friction coefficient (Stribeck minimum), typically between μb and μf.
    pub mu_min: f64,
    /// Viscosity-speed-load number at Stribeck minimum.
    pub sn_min: f64,
}
impl StribeckCurve {
    /// Create a Stribeck curve model.
    pub fn new(mu_boundary: f64, mu_fullfilm: f64, sc: f64) -> Self {
        StribeckCurve {
            mu_boundary,
            mu_fullfilm,
            sc,
            exponent: 1.0,
            mu_min: (mu_boundary + mu_fullfilm) / 2.0,
            sn_min: sc,
        }
    }
    /// Hersey number (Stribeck number) Hs = η * N / P \[m\].
    pub fn hersey_number(&self, eta: f64, n: f64, p: f64) -> f64 {
        if p < 1e-10 {
            return f64::INFINITY;
        }
        eta * n / p
    }
    /// Friction coefficient at Hersey number Hs.
    /// Uses Houpert's analytical Stribeck model.
    pub fn mu_at_hersey(&self, hs: f64) -> f64 {
        let transition_hs = self.sc / 10.0;
        if hs <= transition_hs {
            self.mu_boundary
        } else {
            let x = (hs / self.sc).ln();
            self.mu_min
                + (self.mu_boundary - self.mu_min) * (-x.powi(2) / 2.0).exp()
                + (self.mu_fullfilm - self.mu_min) * (1.0 - (-hs / self.sc).exp())
        }
    }
    /// Regime classification.
    pub fn classify(&self, lambda: f64) -> LubricationRegimeKind {
        if lambda < 1.0 {
            LubricationRegimeKind::Boundary
        } else if lambda < 3.0 {
            LubricationRegimeKind::Mixed
        } else if lambda < 10.0 {
            LubricationRegimeKind::Elastohydrodynamic
        } else {
            LubricationRegimeKind::FullFilm
        }
    }
    /// Viscous shear stress in full-film regime \[Pa\].
    pub fn viscous_shear_stress(&self, eta: f64, v: f64, h: f64) -> f64 {
        if h < 1e-30 {
            return 0.0;
        }
        eta * v / h
    }
}
/// Substrate + coating system with mechanical properties.
#[derive(Debug, Clone)]
pub struct CoatingMaterial {
    /// Substrate Young's modulus \[Pa\].
    pub substrate_modulus: f64,
    /// Coating Young's modulus \[Pa\].
    pub coating_modulus: f64,
    /// Coating thickness \[m\].
    pub coating_thickness: f64,
    /// Vickers hardness of coating \[Pa\].
    pub vickers_hardness: f64,
    /// Adhesion strength (critical load for delamination) \[N\].
    pub adhesion_strength: f64,
    /// Coating density \[kg/m³\].
    pub density: f64,
}
impl CoatingMaterial {
    /// Create a coating system.
    pub fn new(
        substrate_modulus: f64,
        coating_modulus: f64,
        coating_thickness: f64,
        vickers_hardness: f64,
        adhesion_strength: f64,
        density: f64,
    ) -> Self {
        Self {
            substrate_modulus,
            coating_modulus,
            coating_thickness,
            vickers_hardness,
            adhesion_strength,
            density,
        }
    }
    /// Effective composite hardness (depth-dependent, simplified).
    /// Uses Bückle's rule: H_eff decreases from coating toward substrate beyond ~10% thickness.
    pub fn effective_hardness(&self, indent_depth: f64) -> f64 {
        let ratio = indent_depth / self.coating_thickness;
        if ratio < 0.1 {
            self.vickers_hardness
        } else {
            let h_sub = 0.3 * self.vickers_hardness;
            let t = ((ratio - 0.1) / 0.9).min(1.0);
            self.vickers_hardness * (1.0 - t) + h_sub * t
        }
    }
    /// Check if a given scratch load exceeds adhesion limit (delamination risk).
    pub fn will_delaminate(&self, scratch_load: f64) -> bool {
        scratch_load > self.adhesion_strength
    }
    /// Modulus ratio (coating / substrate).
    pub fn modulus_ratio(&self) -> f64 {
        self.coating_modulus / self.substrate_modulus
    }
}
/// Bowden-Tabor adhesion model for friction coefficient.
#[derive(Debug, Clone, Copy)]
pub struct BowdenTaborFriction {
    /// Shear strength of interfacial film s \[Pa\].
    pub shear_strength: f64,
    /// Hardness of softer surface H \[Pa\].
    pub hardness: f64,
    /// Ploughing term contribution (dimensionless, typically 0.01–0.1).
    pub ploughing: f64,
}
impl BowdenTaborFriction {
    /// Create a Bowden-Tabor friction model.
    pub fn new(shear_strength: f64, hardness: f64, ploughing: f64) -> Self {
        BowdenTaborFriction {
            shear_strength,
            hardness,
            ploughing,
        }
    }
    /// Adhesion component of friction coefficient: μ_adh = s / H.
    pub fn mu_adhesion(&self) -> f64 {
        self.shear_strength / self.hardness
    }
    /// Total friction coefficient: μ = μ_adh + μ_plough.
    pub fn mu_total(&self) -> f64 {
        self.mu_adhesion() + self.ploughing
    }
    /// Real contact area fraction: Ar/An = P / H.
    pub fn real_contact_fraction(&self, pressure: f64) -> f64 {
        pressure / self.hardness
    }
    /// Friction force \[N\] for given normal load \[N\].
    pub fn friction_force(&self, normal_load: f64) -> f64 {
        self.mu_total() * normal_load
    }
}
/// Hamrock-Dowson EHL point contact film thickness model.
#[derive(Debug, Clone, Copy)]
pub struct HamrockDowson {
    /// Lubricant dynamic viscosity η0 \[Pa·s\].
    pub eta0: f64,
    /// Pressure-viscosity coefficient α \[1/Pa\].
    pub alpha_visc: f64,
    /// Effective elastic modulus E' \[Pa\].
    pub e_prime: f64,
    /// Effective radius in x-direction Rx \[m\].
    pub rx: f64,
    /// Effective radius in y-direction Ry \[m\].
    pub ry: f64,
}
impl HamrockDowson {
    /// Create a Hamrock-Dowson EHL model.
    pub fn new(eta0: f64, alpha_visc: f64, e_prime: f64, rx: f64, ry: f64) -> Self {
        HamrockDowson {
            eta0,
            alpha_visc,
            e_prime,
            rx,
            ry,
        }
    }
    /// Effective radius Re = sqrt(Rx * Ry).
    pub fn re(&self) -> f64 {
        (self.rx * self.ry).sqrt()
    }
    /// Ellipticity ratio k = Rx/Ry (or Ry/Rx, take >1).
    pub fn ellipticity(&self) -> f64 {
        (self.rx / self.ry).max(self.ry / self.rx)
    }
    /// Dimensionless speed parameter U = η0*u / (E' * Re).
    pub fn speed_param(&self, u: f64) -> f64 {
        self.eta0 * u / (self.e_prime * self.re())
    }
    /// Dimensionless load parameter W = F / (E' * Re²).
    pub fn load_param(&self, force: f64) -> f64 {
        force / (self.e_prime * self.re().powi(2))
    }
    /// Dimensionless material parameter G = α * E'.
    pub fn material_param(&self) -> f64 {
        self.alpha_visc * self.e_prime
    }
    /// Central film thickness hc \[m\] using Hamrock-Dowson (1977).
    /// hc/Re = 2.69 * U^0.67 * G^0.53 * W^(-0.067) * (1 - 0.61*exp(-0.73*k)).
    pub fn central_film_thickness(&self, u: f64, force: f64) -> f64 {
        let uu = self.speed_param(u);
        let ww = self.load_param(force);
        let gg = self.material_param();
        let k = self.ellipticity();
        if ww < 1e-30 {
            return 0.0;
        }
        let hc_re = 2.69
            * uu.powf(0.67)
            * gg.powf(0.53)
            * ww.powf(-0.067)
            * (1.0 - 0.61 * (-0.73 * k).exp());
        hc_re * self.re()
    }
    /// Minimum film thickness hmin \[m\] using Hamrock-Dowson.
    /// hmin/Re = 3.63 * U^0.68 * G^0.49 * W^(-0.073) * (1 - exp(-0.68*k)).
    pub fn minimum_film_thickness(&self, u: f64, force: f64) -> f64 {
        let uu = self.speed_param(u);
        let ww = self.load_param(force);
        let gg = self.material_param();
        let k = self.ellipticity();
        if ww < 1e-30 {
            return 0.0;
        }
        let hmin_re =
            3.63 * uu.powf(0.68) * gg.powf(0.49) * ww.powf(-0.073) * (1.0 - (-0.68 * k).exp());
        hmin_re * self.re()
    }
    /// Lambda ratio Λ = hmin / σ_combined.
    pub fn lambda_ratio(&self, u: f64, force: f64, roughness: f64) -> f64 {
        let hmin = self.minimum_film_thickness(u, force);
        if roughness < 1e-30 {
            return f64::INFINITY;
        }
        hmin / roughness
    }
}
/// Elastohydrodynamic (EHL) film thickness model (Dowson-Higginson).
///
/// Computes the central and minimum film thickness for a lubricated
/// point/line contact using Dowson-Higginson parameters.
#[derive(Debug, Clone, Copy)]
pub struct ElastohydrodynamicFilm {
    /// Ambient dynamic viscosity η₀ \[Pa·s\].
    pub eta0: f64,
    /// Pressure-viscosity coefficient α \[1/Pa\].
    pub alpha: f64,
    /// Entrainment speed u \[m/s\].
    pub u: f64,
    /// Effective radius R* \[m\].
    pub r: f64,
    /// Reduced (plane-strain) modulus E* \[Pa\].
    pub e_star: f64,
}
impl ElastohydrodynamicFilm {
    /// Create an `ElastohydrodynamicFilm`.
    pub fn new(eta0: f64, alpha: f64, u: f64, r: f64, e_star: f64) -> Self {
        Self {
            eta0,
            alpha,
            u,
            r,
            e_star,
        }
    }
    /// Dimensionless speed parameter U = η₀·u / (E*·R).
    fn speed_param(&self) -> f64 {
        self.eta0 * self.u / (self.e_star * self.r)
    }
    /// Dimensionless materials parameter G = α·E*.
    fn material_param(&self) -> f64 {
        self.alpha * self.e_star
    }
    /// Dimensionless load parameter W = F / (E*·R²).
    fn load_param(&self, force: f64) -> f64 {
        force / (self.e_star * self.r * self.r)
    }
    /// Central film thickness \[m\] using Dowson-Higginson formula.
    ///
    /// Hc = 2.69 · U^0.67 · G^0.53 · W^(-0.067) · R
    pub fn central_film_thickness(&self, force: f64) -> f64 {
        let u_dim = self.speed_param();
        let g_dim = self.material_param();
        let w_dim = self.load_param(force).max(1e-30);
        2.69 * u_dim.powf(0.67) * g_dim.powf(0.53) * w_dim.powf(-0.067) * self.r
    }
    /// Minimum film thickness \[m\]; approximately 75 % of central.
    ///
    /// Hmin = 0.75 · Hc  (empirical approximation)
    pub fn minimum_film_thickness(&self, force: f64) -> f64 {
        0.75 * self.central_film_thickness(force)
    }
    /// Specific film thickness (lambda ratio) Λ = h_min / σ_composite.
    ///
    /// `composite_roughness` σ = √(σ₁² + σ₂²) \[m\].
    pub fn lambda_ratio(&self, composite_roughness: f64, force: f64) -> f64 {
        let sigma = composite_roughness.max(1e-30);
        self.minimum_film_thickness(force) / sigma
    }
}
/// DMT adhesive contact model for stiff, hard materials.
#[derive(Debug, Clone, Copy)]
pub struct DerjaguinMullerToporov {
    /// Effective modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R* \[m\].
    pub effective_radius: f64,
    /// Work of adhesion w \[J/m²\].
    pub work_of_adhesion: f64,
}
impl DerjaguinMullerToporov {
    /// Create a DMT model.
    pub fn new(effective_modulus: f64, effective_radius: f64, work_of_adhesion: f64) -> Self {
        Self {
            effective_modulus,
            effective_radius,
            work_of_adhesion,
        }
    }
    /// Pull-off force in DMT model: F_pull = -2 * π * w * R*.
    pub fn pull_off_force(&self) -> f64 {
        -2.0 * PI * self.work_of_adhesion * self.effective_radius
    }
    /// Contact radius: same as Hertz with effective load F + |F_pull|.
    pub fn contact_radius(&self, force: f64) -> f64 {
        let f_pull = self.pull_off_force().abs();
        let f_eff = force + f_pull;
        if f_eff <= 0.0 {
            return 0.0;
        }
        let h = HertzContact {
            effective_modulus: self.effective_modulus,
            effective_radius: self.effective_radius,
        };
        h.compute(f_eff).contact_radius
    }
    /// Tabor parameter (determines JKR vs DMT regime).
    /// μ = (R* * w² / (E*² * z₀³))^(1/3), z₀ = atomic cutoff ≈ 0.3 nm.
    pub fn tabor_parameter(&self, z0: f64) -> f64 {
        let num = self.effective_radius * self.work_of_adhesion * self.work_of_adhesion;
        let den = self.effective_modulus * self.effective_modulus * z0 * z0 * z0;
        (num / den).powf(1.0 / 3.0)
    }
}
/// Archard's wear law model.
///
/// Volume worn Q = k_archard · F · s / H, where F is normal force \[N\],
/// s is sliding distance \[m\], and H is hardness \[Pa\].
#[derive(Debug, Clone, Copy)]
pub struct WearModel {
    /// Dimensionless Archard wear coefficient k \[–\].
    pub k_archard: f64,
    /// Hardness of the softer material \[Pa\].
    pub hardness: f64,
}
impl WearModel {
    /// Create a `WearModel`.
    pub fn new(k_archard: f64, hardness: f64) -> Self {
        Self {
            k_archard,
            hardness,
        }
    }
    /// Wear volume \[m³\] = k · F · s / H.
    pub fn wear_volume(&self, normal_force: f64, sliding_distance: f64) -> f64 {
        self.k_archard * normal_force * sliding_distance / self.hardness
    }
    /// Wear rate \[m³/s\] = k · F · v / H.
    pub fn wear_rate(&self, normal_force: f64, sliding_speed: f64) -> f64 {
        self.k_archard * normal_force * sliding_speed / self.hardness
    }
}
/// Holm's wear model (electrical contact inspired, generalised to tribology).
/// Volume per unit distance = Z * W (Holm equation).
#[derive(Debug, Clone, Copy)]
pub struct HolmWear {
    /// Holm wear coefficient Z \[m²/N\] (= K/H in Archard terms).
    pub z: f64,
}
impl HolmWear {
    /// Create a Holm wear model with given Z coefficient.
    pub fn new(z: f64) -> Self {
        HolmWear { z }
    }
    /// Wear volume rate \[m³/m\] = Z * W.
    pub fn volume_rate(&self, normal_load: f64) -> f64 {
        self.z * normal_load
    }
    /// Volume loss \[m³\] over sliding distance s \[m\].
    pub fn volume_loss(&self, normal_load: f64, s: f64) -> f64 {
        self.z * normal_load * s
    }
    /// Wear depth \[m\] over area A \[m²\] and distance s \[m\].
    pub fn wear_depth(&self, normal_load: f64, s: f64, area: f64) -> f64 {
        if area < 1e-30 {
            return 0.0;
        }
        self.volume_loss(normal_load, s) / area
    }
}
/// 2D Hertz contact for two parallel cylinders (line contact).
#[derive(Debug, Clone, Copy)]
pub struct CylinderContact {
    /// Effective modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R* \[m\].
    pub effective_radius: f64,
    /// Contact half-width a \[m\] per unit length.
    pub contact_half_width: f64,
}
impl CylinderContact {
    /// Create a cylinder-on-cylinder contact (per unit length).
    pub fn new(e1: f64, nu1: f64, r1: f64, e2: f64, nu2: f64, r2: f64) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        let r_star = r1 * r2 / (r1 + r2);
        CylinderContact {
            effective_modulus: e_star,
            effective_radius: r_star,
            contact_half_width: 0.0,
        }
    }
    /// Contact half-width a \[m\] per unit length under load w \[N/m\].
    pub fn half_width(&self, w: f64) -> f64 {
        (4.0 * w * self.effective_radius / (PI * self.effective_modulus)).sqrt()
    }
    /// Maximum pressure p0 \[Pa\] per unit length.
    pub fn max_pressure(&self, w: f64) -> f64 {
        let a = self.half_width(w);
        if a < 1e-30 {
            return 0.0;
        }
        2.0 * w / (PI * a)
    }
    /// Subsurface shear stress τmax \[Pa\] at depth z = 0.786a below surface.
    pub fn max_subsurface_shear(&self, w: f64) -> f64 {
        0.300 * self.max_pressure(w)
    }
    /// Contact pressure distribution at position x from centre \[Pa\].
    pub fn pressure_at(&self, w: f64, x: f64) -> f64 {
        let a = self.half_width(w);
        let p0 = self.max_pressure(w);
        if x.abs() >= a {
            return 0.0;
        }
        p0 * (1.0 - (x / a).powi(2)).sqrt()
    }
    /// Depth of maximum shear stress z_max \[m\].
    pub fn depth_max_shear(&self, w: f64) -> f64 {
        0.786 * self.half_width(w)
    }
}
/// State of a tribological system at a given operating point.
#[derive(Debug, Clone)]
pub struct TribologicalState {
    /// Sliding speed v \[m/s\].
    pub velocity: f64,
    /// Normal load W \[N\].
    pub load: f64,
    /// Contact temperature T \[K\].
    pub temperature: f64,
    /// Film thickness h \[m\].
    pub film_thickness: f64,
    /// Lambda ratio (film parameter) Λ.
    pub lambda: f64,
    /// Current friction coefficient μ.
    pub mu: f64,
    /// Accumulated wear volume Vtotal \[m³\].
    pub wear_volume: f64,
    /// Lubrication regime.
    pub regime: LubricationRegimeKind,
}
impl TribologicalState {
    /// Create an initial tribological state.
    pub fn new(velocity: f64, load: f64) -> Self {
        TribologicalState {
            velocity,
            load,
            temperature: 293.15,
            film_thickness: 0.0,
            lambda: 0.0,
            mu: 0.15,
            wear_volume: 0.0,
            regime: LubricationRegimeKind::Boundary,
        }
    }
    /// Update state with new film thickness and roughness.
    pub fn update_film(&mut self, h: f64, roughness: f64) {
        self.film_thickness = h;
        self.lambda = if roughness > 1e-30 {
            h / roughness
        } else {
            f64::INFINITY
        };
        self.regime = if self.lambda < 1.0 {
            LubricationRegimeKind::Boundary
        } else if self.lambda < 3.0 {
            LubricationRegimeKind::Mixed
        } else if self.lambda < 10.0 {
            LubricationRegimeKind::Elastohydrodynamic
        } else {
            LubricationRegimeKind::FullFilm
        };
    }
    /// Accumulate wear volume over time step dt \[s\] using Archard model.
    pub fn accumulate_wear(&mut self, k: f64, hardness: f64, dt: f64) {
        let rate = k * self.load * self.velocity / hardness;
        self.wear_volume += rate * dt;
    }
    /// Update friction coefficient based on current lambda.
    pub fn update_friction(&mut self, mu_boundary: f64, mu_full: f64) {
        let t = ((self.lambda - 1.0) / 9.0).clamp(0.0, 1.0);
        self.mu = mu_boundary * (1.0 - t) + mu_full * t;
    }
    /// Power dissipation \[W\].
    pub fn power_dissipation(&self) -> f64 {
        self.mu * self.load * self.velocity
    }
}
/// Stick-slip friction oscillator (spring-mass with Coulomb friction on a
/// moving belt).
///
/// The mass is connected to a fixed spring (`stiffness`) and rests on a belt
/// moving at speed `v_drive`.  Friction alternates between static
/// (`mu_s`) and kinetic (`mu_k`) regimes.
#[derive(Debug, Clone)]
pub struct FrictionOscillator {
    /// Mass \[kg\].
    pub mass: f64,
    /// Spring stiffness \[N/m\].
    pub stiffness: f64,
    /// Static friction coefficient.
    pub mu_s: f64,
    /// Kinetic friction coefficient.
    pub mu_k: f64,
}
impl FrictionOscillator {
    /// Create a `FrictionOscillator`.
    pub fn new(mass: f64, stiffness: f64, mu_s: f64, mu_k: f64) -> Self {
        Self {
            mass,
            stiffness,
            mu_s,
            mu_k,
        }
    }
    /// Returns `true` when the mass is sticking to the belt.
    ///
    /// The mass sticks when its velocity equals the belt drive speed.
    pub fn is_sticking(&self, v: f64, v_drive: f64) -> bool {
        (v - v_drive).abs() < 1e-9
    }
    /// Advance one time step using semi-implicit Euler.
    ///
    /// `x` is the current displacement \[m\], `v` is the current velocity
    /// \[m/s\], `v_drive` is the belt speed \[m/s\], `dt` is the time step \[s\].
    ///
    /// Returns `(new_x, new_v)`.
    pub fn step(&mut self, x: f64, v: f64, v_drive: f64, dt: f64) -> (f64, f64) {
        let g = 9.81;
        let f_normal = self.mass * g;
        let f_spring = -self.stiffness * x;
        let relative_v = v - v_drive;
        let f_friction = if relative_v.abs() < 1e-9 {
            let f_static_max = self.mu_s * f_normal;
            (-f_spring).clamp(-f_static_max, f_static_max)
        } else {
            -relative_v.signum() * self.mu_k * f_normal
        };
        let a = (f_spring + f_friction) / self.mass;
        let new_v = v + a * dt;
        let new_x = x + new_v * dt;
        (new_x, new_v)
    }
}
/// JKR adhesive contact model for soft, compliant materials.
#[derive(Debug, Clone, Copy)]
pub struct JohnsonKendallRoberts {
    /// Effective elastic modulus E* \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R* \[m\].
    pub effective_radius: f64,
    /// Work of adhesion w \[J/m²\].
    pub work_of_adhesion: f64,
}
impl JohnsonKendallRoberts {
    /// Create a JKR model.
    pub fn new(effective_modulus: f64, effective_radius: f64, work_of_adhesion: f64) -> Self {
        Self {
            effective_modulus,
            effective_radius,
            work_of_adhesion,
        }
    }
    /// Pull-off force: F_pull = -1.5 * π * w * R*.
    pub fn pull_off_force(&self) -> f64 {
        -1.5 * PI * self.work_of_adhesion * self.effective_radius
    }
    /// Contact radius under JKR theory for applied load F.
    pub fn contact_radius(&self, force: f64) -> f64 {
        let w = self.work_of_adhesion;
        let r = self.effective_radius;
        let e = self.effective_modulus;
        let term1 = 3.0 * PI * w * r;
        let discriminant = 6.0 * PI * w * r * force + term1 * term1;
        if discriminant < 0.0 {
            return 0.0;
        }
        let inner = force + term1 + discriminant.sqrt();
        let a3 = r / e * inner;
        a3.powf(1.0 / 3.0)
    }
    /// Load for a given contact radius (inverse JKR).
    pub fn load_for_radius(&self, a: f64) -> f64 {
        let e = self.effective_modulus;
        let r = self.effective_radius;
        let w = self.work_of_adhesion;
        let hertz_term = 4.0 * e * a * a * a / (3.0 * r);
        let adhesion_term = (8.0 * PI * e * w * a * a * a).sqrt();
        hertz_term - adhesion_term
    }
}
/// Elastohydrodynamic lubrication (EHL) film thickness using Dowson-Higginson formula.
#[derive(Debug, Clone, Copy)]
pub struct EhlFilm {
    /// Ambient viscosity η₀ \[Pa·s\].
    pub viscosity: f64,
    /// Pressure-viscosity coefficient α \[1/Pa\].
    pub pressure_viscosity_coeff: f64,
    /// Effective elastic modulus E' \[Pa\].
    pub effective_modulus: f64,
    /// Effective radius R' \[m\].
    pub effective_radius: f64,
}
impl EhlFilm {
    /// Create an EHL film model.
    pub fn new(
        viscosity: f64,
        pressure_viscosity_coeff: f64,
        effective_modulus: f64,
        effective_radius: f64,
    ) -> Self {
        Self {
            viscosity,
            pressure_viscosity_coeff,
            effective_modulus,
            effective_radius,
        }
    }
    /// Compute central film thickness h_c using Dowson-Higginson (1977).
    /// `u` = entrainment velocity \[m/s\], `w` = load per unit length \[N/m\] for line contact.
    pub fn central_film_thickness_line(&self, u: f64, w: f64) -> f64 {
        let ue = (self.viscosity * u) / (self.effective_modulus * self.effective_radius);
        let we = w / (self.effective_modulus * self.effective_radius);
        let ge = self.pressure_viscosity_coeff * self.effective_modulus;
        if we < 1e-30 {
            return 0.0;
        }
        let h_c_r = 2.65 * ge.powf(0.54) * ue.powf(0.7) / we.powf(0.13);
        h_c_r * self.effective_radius
    }
    /// Compute minimum film thickness h_min for point contact (Hamrock-Dowson).
    /// `u` = entrainment velocity \[m/s\], `w` = normal load \[N\].
    pub fn minimum_film_thickness_point(&self, u: f64, w: f64) -> f64 {
        let ue = (self.viscosity * u) / (self.effective_modulus * self.effective_radius);
        let we = w / (self.effective_modulus * self.effective_radius * self.effective_radius);
        let ge = self.pressure_viscosity_coeff * self.effective_modulus;
        if we < 1e-30 {
            return 0.0;
        }
        let h_min_r = 3.63 * ue.powf(0.68) * ge.powf(0.49) * we.powf(-0.073);
        h_min_r * self.effective_radius
    }
    /// Lambda ratio (film parameter): Λ = h_min / σ.
    pub fn lambda_ratio(&self, u: f64, w: f64, roughness: f64) -> f64 {
        let h_min = self.minimum_film_thickness_point(u, w);
        if roughness < 1e-30 {
            f64::INFINITY
        } else {
            h_min / roughness
        }
    }
}
/// Hertz contact mechanics parameterised by the individual sphere/surface data.
///
/// Stores the raw material and geometry parameters (radii, moduli, Poisson
/// ratios) and derives the reduced quantities on demand.
#[derive(Debug, Clone, Copy)]
pub struct HertzContactParams {
    /// Radius of body 1 \[m\] (use `f64::INFINITY` for a flat).
    pub r1: f64,
    /// Radius of body 2 \[m\] (use `f64::INFINITY` for a flat).
    pub r2: f64,
    /// Young's modulus of body 1 \[Pa\].
    pub e1: f64,
    /// Poisson's ratio of body 1 \[–\].
    pub nu1: f64,
    /// Young's modulus of body 2 \[Pa\].
    pub e2: f64,
    /// Poisson's ratio of body 2 \[–\].
    pub nu2: f64,
}
impl HertzContactParams {
    /// Create a `HertzContactParams` from individual surface parameters.
    pub fn new(r1: f64, r2: f64, e1: f64, nu1: f64, e2: f64, nu2: f64) -> Self {
        Self {
            r1,
            r2,
            e1,
            nu1,
            e2,
            nu2,
        }
    }
    /// Reduced (plane-strain) modulus E* = 1/\[(1−ν₁²)/E₁ + (1−ν₂²)/E₂\].
    pub fn reduced_modulus(&self) -> f64 {
        1.0 / ((1.0 - self.nu1 * self.nu1) / self.e1 + (1.0 - self.nu2 * self.nu2) / self.e2)
    }
    /// Reduced radius R* = R₁·R₂ / (R₁ + R₂).
    ///
    /// When one body is flat (`r = ∞`) this correctly returns the other radius.
    pub fn reduced_radius(&self) -> f64 {
        if self.r1.is_infinite() {
            return self.r2;
        }
        if self.r2.is_infinite() {
            return self.r1;
        }
        self.r1 * self.r2 / (self.r1 + self.r2)
    }
    /// Contact radius a = (3·F·R* / (4·E*))^(1/3) \[m\].
    pub fn contact_radius(&self, force: f64) -> f64 {
        let e_star = self.reduced_modulus();
        let r_star = self.reduced_radius();
        ((3.0 * force * r_star) / (4.0 * e_star)).powf(1.0 / 3.0)
    }
    /// Maximum Hertz pressure p₀ = 3F / (2π·a²) \[Pa\].
    pub fn max_pressure(&self, force: f64) -> f64 {
        let a = self.contact_radius(force);
        if a < 1e-30 {
            return 0.0;
        }
        3.0 * force / (2.0 * PI * a * a)
    }
    /// Contact stiffness k = dF/dδ = 2·E*·a \[N/m\].
    pub fn contact_stiffness(&self, force: f64) -> f64 {
        let e_star = self.reduced_modulus();
        let a = self.contact_radius(force);
        2.0 * e_star * a
    }
    /// Maximum sub-surface shear stress τ_max ≈ 0.31·p₀ at depth z ≈ 0.48·a.
    pub fn max_shear_stress(&self, force: f64) -> f64 {
        0.31 * self.max_pressure(force)
    }
}
/// Stribeck curve model for lubrication regime determination.
#[derive(Debug, Clone)]
pub struct LubricationRegime {
    /// Boundary friction coefficient μ_b.
    pub mu_boundary: f64,
    /// Full-film friction coefficient μ_f.
    pub mu_fullfilm: f64,
    /// Stribeck parameter: characteristic viscosity-speed-load number.
    pub stribeck_constant: f64,
}
impl LubricationRegime {
    /// Create a Stribeck model.
    pub fn new(mu_boundary: f64, mu_fullfilm: f64, stribeck_constant: f64) -> Self {
        Self {
            mu_boundary,
            mu_fullfilm,
            stribeck_constant,
        }
    }
    /// Classify regime from lambda ratio.
    pub fn classify(&self, lambda: f64) -> LubricationRegimeKind {
        if lambda < 1.0 {
            LubricationRegimeKind::Boundary
        } else if lambda < 3.0 {
            LubricationRegimeKind::Mixed
        } else if lambda < 10.0 {
            LubricationRegimeKind::Elastohydrodynamic
        } else {
            LubricationRegimeKind::FullFilm
        }
    }
    /// Compute effective friction coefficient from Stribeck number S = η*v/P.
    /// Uses exponential transition model.
    pub fn friction_coefficient(&self, eta: f64, velocity: f64, pressure: f64) -> f64 {
        if pressure < 1e-10 {
            return self.mu_boundary;
        }
        let stribeck_num = eta * velocity / pressure;
        self.mu_fullfilm
            + (self.mu_boundary - self.mu_fullfilm)
                * (-(stribeck_num / self.stribeck_constant)).exp()
    }
}
