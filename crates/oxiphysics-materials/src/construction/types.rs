//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::terzaghi_bearing_factors;
use std::f64::consts::PI;

/// Concrete cover to reinforcement per ACI 318.
pub struct CoverRequirement;
impl CoverRequirement {
    /// Minimum clear cover (mm) for given exposure class and bar size (mm diameter).
    pub fn minimum_cover(exposure: &ExposureClass, bar_dia: f64) -> f64 {
        match exposure {
            ExposureClass::Interior => {
                if bar_dia <= 16.0 {
                    20.0
                } else {
                    40.0
                }
            }
            ExposureClass::Moderate => 50.0,
            ExposureClass::Severe => 65.0,
            ExposureClass::Submerged => 75.0,
        }
    }
    /// Effective depth given total section depth h (mm) and cover + bar dia/2.
    pub fn effective_depth(h: f64, cover: f64, bar_dia: f64) -> f64 {
        h - cover - bar_dia / 2.0
    }
}
/// Stirrup (shear reinforcement) design for reinforced concrete beams (ACI 318).
#[derive(Debug, Clone)]
pub struct StirrupDesign {
    /// Beam width bw (mm).
    pub bw: f64,
    /// Effective depth d (mm).
    pub d: f64,
    /// Concrete f'c (MPa).
    pub fc: f64,
    /// Steel yield strength fy (MPa).
    pub fy: f64,
    /// Stirrup bar area per leg Av (mm²/leg).
    pub av: f64,
    /// Number of legs per stirrup.
    pub legs: u32,
    /// Stirrup spacing s (mm).
    pub spacing: f64,
}
impl StirrupDesign {
    /// Create a stirrup design.
    pub fn new(bw: f64, d: f64, fc: f64, fy: f64, av: f64, legs: u32, spacing: f64) -> Self {
        StirrupDesign {
            bw,
            d,
            fc,
            fy,
            av,
            legs,
            spacing,
        }
    }
    /// Concrete shear contribution Vc (N) per ACI 318.
    pub fn vc(&self) -> f64 {
        0.17 * self.fc.sqrt() * self.bw * self.d
    }
    /// Steel shear contribution Vs (N) per ACI 318.
    pub fn vs(&self) -> f64 {
        let total_av = self.av * self.legs as f64;
        total_av * self.fy * self.d / self.spacing
    }
    /// Nominal shear capacity Vn = Vc + Vs (N).
    pub fn vn(&self) -> f64 {
        self.vc() + self.vs()
    }
    /// Design shear capacity φVn (φ = 0.75).
    pub fn phi_vn(&self) -> f64 {
        0.75 * self.vn()
    }
    /// Maximum stirrup spacing per ACI 318 (d/2 or 600 mm, whichever smaller).
    pub fn max_spacing(&self) -> f64 {
        (self.d / 2.0).min(600.0)
    }
    /// Minimum stirrup area per ACI 318 (mm²/mm): Av_min/s = max(0.062*sqrt(fc)*bw/fy, 0.35*bw/fy).
    pub fn min_av_per_s(&self) -> f64 {
        let t1 = 0.062 * self.fc.sqrt() * self.bw / self.fy;
        let t2 = 0.35 * self.bw / self.fy;
        t1.max(t2)
    }
    /// Check if stirrup spacing is adequate for Vu.
    pub fn is_adequate(&self, vu: f64) -> bool {
        self.phi_vn() >= vu
    }
}
/// Equal-leg angle section properties.
#[derive(Debug, Clone)]
pub struct AngleSection {
    /// Leg length L (mm).
    pub l: f64,
    /// Leg thickness t (mm).
    pub t: f64,
    /// Yield strength Fy (MPa).
    pub fy: f64,
}
impl AngleSection {
    /// Create an equal-leg angle section.
    pub fn equal_leg(l: f64, t: f64, fy: f64) -> Self {
        AngleSection { l, t, fy }
    }
    /// Cross-sectional area (mm²).
    pub fn area(&self) -> f64 {
        2.0 * self.l * self.t - self.t.powi(2)
    }
    /// Centroid location (mm) from back of angle.
    pub fn centroid(&self) -> f64 {
        let a = self.area();
        if a < 1e-10 {
            return 0.0;
        }
        let a1 = self.l * self.t;
        let x1 = self.t / 2.0;
        let a2 = (self.l - self.t) * self.t;
        let x2 = self.t + (self.l - self.t) / 2.0;
        (a1 * x1 + a2 * x2) / a
    }
    /// Ixx about geometric axis (mm⁴) — approximate.
    pub fn ixx(&self) -> f64 {
        let t = self.t;
        let l = self.l;
        let i_horiz = l * t.powi(3) / 12.0 + l * t * (t / 2.0).powi(2);
        let i_vert = t * (l - t).powi(3) / 12.0 + t * (l - t) * ((l - t) / 2.0 + t).powi(2);
        i_horiz + i_vert
    }
    /// Radius of gyration rz (mm) about minimum axis (approx).
    pub fn rz(&self) -> f64 {
        (self.ixx() / self.area()).sqrt()
    }
    /// Axial capacity Pn (N) for short column.
    pub fn axial_capacity(&self) -> f64 {
        self.area() * self.fy
    }
}
/// Retaining wall earth pressure analysis (Rankine/Coulomb).
pub struct RetainingWall {
    /// Height of wall H (m).
    pub height: f64,
    /// Soil unit weight behind wall γ (kN/m³).
    pub gamma: f64,
    /// Friction angle φ (degrees).
    pub phi_deg: f64,
    /// Wall friction angle δ (degrees).
    pub delta_deg: f64,
    /// Backfill inclination β (degrees).
    pub beta_deg: f64,
}
impl RetainingWall {
    /// Create a retaining wall with given parameters.
    pub fn new(height: f64, gamma: f64, phi_deg: f64) -> Self {
        RetainingWall {
            height,
            gamma,
            phi_deg,
            delta_deg: phi_deg * 2.0 / 3.0,
            beta_deg: 0.0,
        }
    }
    /// Rankine active earth pressure coefficient Ka.
    pub fn rankine_ka(&self) -> f64 {
        let phi = self.phi_deg.to_radians();
        ((PI / 4.0 - phi / 2.0).tan()).powi(2)
    }
    /// Rankine passive earth pressure coefficient Kp.
    pub fn rankine_kp(&self) -> f64 {
        let phi = self.phi_deg.to_radians();
        ((PI / 4.0 + phi / 2.0).tan()).powi(2)
    }
    /// Total active thrust Pa per unit width (kN/m) using Rankine.
    pub fn active_thrust_rankine(&self) -> f64 {
        0.5 * self.gamma * self.height.powi(2) * self.rankine_ka()
    }
    /// Total passive resistance Pp per unit width (kN/m) using Rankine.
    pub fn passive_resistance_rankine(&self) -> f64 {
        0.5 * self.gamma * self.height.powi(2) * self.rankine_kp()
    }
    /// Coulomb active pressure coefficient Ka.
    pub fn coulomb_ka(&self) -> f64 {
        let phi = self.phi_deg.to_radians();
        let delta = self.delta_deg.to_radians();
        let beta = self.beta_deg.to_radians();
        let alpha = PI / 2.0;
        let num = (phi - beta).sin().powi(2);
        let term = 1.0
            + ((phi + delta).sin() * (phi - beta).sin()
                / ((alpha + delta).sin() * (alpha - beta).sin()))
            .sqrt();
        num / ((alpha).sin().powi(2) * (alpha - delta).sin() * term.powi(2))
    }
    /// Check overturning stability factor (should be > 2.0).
    pub fn overturning_factor(&self, footing_width: f64) -> f64 {
        let pa = self.active_thrust_rankine();
        let overturning_moment = pa * self.height / 3.0;
        let restoring_moment = self.gamma * self.height * footing_width * footing_width / 2.0;
        if overturning_moment < 1e-10 {
            return f64::INFINITY;
        }
        restoring_moment / overturning_moment
    }
}
/// Masonry prism test data for determining f'm (compressive strength).
#[derive(Debug, Clone)]
pub struct MasonryPrism {
    /// Unit compressive strength (MPa).
    pub unit_strength: f64,
    /// Mortar type (S = 12.4 MPa, N = 5.2 MPa, M = 17.2 MPa minimum).
    pub mortar_type: String,
    /// Prism height-to-thickness ratio (h/t), typically 2–5.
    pub h_t_ratio: f64,
    /// Measured prism strength fp (MPa).
    pub fp_measured: f64,
}
impl MasonryPrism {
    /// Create a masonry prism.
    pub fn new(unit_strength: f64, mortar_type: &str, h_t_ratio: f64, fp_measured: f64) -> Self {
        MasonryPrism {
            unit_strength,
            mortar_type: mortar_type.to_string(),
            h_t_ratio,
            fp_measured,
        }
    }
    /// Correction factor for h/t ratio (TMS 602): accounts for slenderness.
    pub fn ht_correction_factor(&self) -> f64 {
        if self.h_t_ratio >= 5.0 {
            1.0
        } else if self.h_t_ratio >= 2.0 {
            0.75 + 0.05 * (self.h_t_ratio - 2.0)
        } else {
            0.75
        }
    }
    /// Corrected prism strength f'm (MPa).
    pub fn fm_corrected(&self) -> f64 {
        self.fp_measured * self.ht_correction_factor()
    }
    /// Modulus of elasticity Em per TMS 402: Em = 900 * f'm.
    pub fn em(&self) -> f64 {
        900.0 * self.fm_corrected()
    }
    /// Allowable compressive stress Fa (MPa) for unreinforced masonry (TMS 402).
    pub fn allowable_compressive_stress(&self) -> f64 {
        0.25 * self.fm_corrected()
    }
}
/// Concrete mix design with water-cement ratio, admixtures, and aggregate grading.
#[derive(Debug, Clone)]
pub struct ConcreteMixDesign {
    /// Water-cement ratio (w/c) by mass.
    pub wc_ratio: f64,
    /// Cement content (kg/m³).
    pub cement: f64,
    /// Water content (kg/m³).
    pub water: f64,
    /// Fine aggregate content (kg/m³).
    pub fine_agg: f64,
    /// Coarse aggregate content (kg/m³).
    pub coarse_agg: f64,
    /// Superplasticiser dosage (% bwc).
    pub superplasticiser: f64,
    /// Air entrainment (%).
    pub air_content: f64,
    /// Slump (mm).
    pub slump: f64,
}
impl ConcreteMixDesign {
    /// Create a standard Normal Portland Cement mix.
    pub fn new(
        wc_ratio: f64,
        cement: f64,
        fine_agg: f64,
        coarse_agg: f64,
        air_content: f64,
    ) -> Self {
        let water = wc_ratio * cement;
        ConcreteMixDesign {
            wc_ratio,
            cement,
            water,
            fine_agg,
            coarse_agg,
            superplasticiser: 0.0,
            air_content,
            slump: 75.0,
        }
    }
    /// Create a typical C30 normal-weight mix (approximate ACI 211 proportioning).
    pub fn c30_normal_weight() -> Self {
        ConcreteMixDesign {
            wc_ratio: 0.50,
            cement: 350.0,
            water: 175.0,
            fine_agg: 720.0,
            coarse_agg: 1080.0,
            superplasticiser: 0.0,
            air_content: 2.0,
            slump: 75.0,
        }
    }
    /// Create a high-performance mix with superplasticiser.
    pub fn high_performance() -> Self {
        ConcreteMixDesign {
            wc_ratio: 0.30,
            cement: 500.0,
            water: 150.0,
            fine_agg: 700.0,
            coarse_agg: 1050.0,
            superplasticiser: 1.5,
            air_content: 1.5,
            slump: 180.0,
        }
    }
    /// Fresh concrete unit weight (kg/m³), ignoring air.
    pub fn fresh_unit_weight(&self) -> f64 {
        self.cement + self.water + self.fine_agg + self.coarse_agg
    }
    /// Estimated 28-day compressive strength (MPa) via Abrams' law.
    /// fc = A / B^(w/c), with A=96 MPa, B=8.6 for OPC.
    pub fn abrams_strength(&self) -> f64 {
        let a = 96.0_f64;
        let b = 8.6_f64;
        a / b.powf(self.wc_ratio)
    }
    /// Volume of paste (m³/m³ concrete).
    pub fn paste_volume(&self) -> f64 {
        let rho_c = 3150.0;
        let rho_w = 1000.0;
        self.cement / rho_c + self.water / rho_w
    }
    /// Aggregate-cement ratio by mass.
    pub fn aggregate_cement_ratio(&self) -> f64 {
        (self.fine_agg + self.coarse_agg) / self.cement
    }
    /// Check Duff Abrams water-cement ratio limit for durability (≤ 0.50 for exposure class C2).
    pub fn meets_durability_wc_limit(&self, max_wc: f64) -> bool {
        self.wc_ratio <= max_wc
    }
    /// Fineness modulus of combined aggregate (simplified weighted average).
    /// `fm_fine` = fineness modulus of fine agg, `fm_coarse` = coarse.
    pub fn combined_fineness_modulus(&self, fm_fine: f64, fm_coarse: f64) -> f64 {
        let total = self.fine_agg + self.coarse_agg;
        if total < 1e-6 {
            return 0.0;
        }
        (self.fine_agg * fm_fine + self.coarse_agg * fm_coarse) / total
    }
}
/// Primary consolidation settlement analysis (Terzaghi 1D consolidation).
#[derive(Debug, Clone)]
pub struct TerzaghiSettlement {
    /// Compression index Cc (dimensionless).
    pub cc: f64,
    /// Recompression index Cr (dimensionless).
    pub cr: f64,
    /// Initial void ratio e0.
    pub e0: f64,
    /// Layer thickness H (m).
    pub h: f64,
    /// Initial vertical effective stress σ'v0 (kPa).
    pub sigma_v0: f64,
    /// Preconsolidation pressure σ'p (kPa).
    pub sigma_p: f64,
    /// Coefficient of consolidation cv (m²/year).
    pub cv: f64,
}
impl TerzaghiSettlement {
    /// Create a normally consolidated clay layer.
    pub fn normally_consolidated(cc: f64, e0: f64, h: f64, sigma_v0: f64, cv: f64) -> Self {
        TerzaghiSettlement {
            cc,
            cr: cc / 5.0,
            e0,
            h,
            sigma_v0,
            sigma_p: sigma_v0,
            cv,
        }
    }
    /// Create an over-consolidated clay layer.
    pub fn over_consolidated(
        cc: f64,
        cr: f64,
        e0: f64,
        h: f64,
        sigma_v0: f64,
        sigma_p: f64,
        cv: f64,
    ) -> Self {
        TerzaghiSettlement {
            cc,
            cr,
            e0,
            h,
            sigma_v0,
            sigma_p,
            cv,
        }
    }
    /// Primary settlement Sc (m) for stress increment Δσ (kPa).
    pub fn primary_settlement(&self, delta_sigma: f64) -> f64 {
        let sigma_f = self.sigma_v0 + delta_sigma;
        if self.sigma_v0 >= self.sigma_p {
            self.h * self.cc / (1.0 + self.e0) * (sigma_f / self.sigma_v0).log10()
        } else if sigma_f <= self.sigma_p {
            self.h * self.cr / (1.0 + self.e0) * (sigma_f / self.sigma_v0).log10()
        } else {
            let s1 = self.h * self.cr / (1.0 + self.e0) * (self.sigma_p / self.sigma_v0).log10();
            let s2 = self.h * self.cc / (1.0 + self.e0) * (sigma_f / self.sigma_p).log10();
            s1 + s2
        }
    }
    /// Time factor Tv for a given consolidation degree Uv.
    /// Uses Terzaghi's time factor: Tv ≈ π/4 * Uv² for Uv ≤ 0.6.
    pub fn time_factor(&self, uv: f64) -> f64 {
        if uv <= 0.6 {
            PI / 4.0 * uv.powi(2)
        } else {
            -(1.0 - uv).ln() * 1.781 - 0.933
        }
    }
    /// Time t (years) to reach consolidation degree Uv (0–1) for double drainage.
    pub fn consolidation_time(&self, uv: f64) -> f64 {
        let tv = self.time_factor(uv);
        let hdr = self.h / 2.0;
        tv * hdr.powi(2) / self.cv
    }
    /// Overconsolidation ratio OCR = σ'p / σ'v0.
    pub fn ocr(&self) -> f64 {
        self.sigma_p / self.sigma_v0
    }
}
/// Sieve analysis results for aggregate grading.
#[derive(Debug, Clone)]
pub struct AggregateGrading {
    /// Sieve sizes in mm (ascending order).
    pub sieve_sizes: Vec<f64>,
    /// Percent passing each sieve (0–100).
    pub percent_passing: Vec<f64>,
}
impl AggregateGrading {
    /// Create from parallel sieve / passing arrays.
    pub fn new(sieve_sizes: Vec<f64>, percent_passing: Vec<f64>) -> Self {
        AggregateGrading {
            sieve_sizes,
            percent_passing,
        }
    }
    /// Create a typical ASTM C33 coarse aggregate grading (19 mm max size).
    pub fn astm_c33_coarse() -> Self {
        AggregateGrading {
            sieve_sizes: vec![25.0, 19.0, 12.5, 9.5, 4.75, 2.36],
            percent_passing: vec![100.0, 90.0, 55.0, 35.0, 5.0, 0.0],
        }
    }
    /// Create a typical ASTM C33 fine aggregate grading.
    pub fn astm_c33_fine() -> Self {
        AggregateGrading {
            sieve_sizes: vec![9.5, 4.75, 2.36, 1.18, 0.60, 0.30, 0.15],
            percent_passing: vec![100.0, 95.0, 80.0, 65.0, 45.0, 20.0, 5.0],
        }
    }
    /// Fineness modulus: sum of cumulative % retained on standard sieves / 100.
    pub fn fineness_modulus(&self) -> f64 {
        let sum_retained: f64 = self.percent_passing.iter().map(|p| 100.0 - p).sum();
        sum_retained / 100.0
    }
    /// Maximum aggregate size (mm) — largest sieve with 100% passing.
    pub fn maximum_size(&self) -> f64 {
        for (i, &p) in self.percent_passing.iter().enumerate() {
            if p >= 100.0 {
                return self.sieve_sizes[i];
            }
        }
        *self.sieve_sizes.last().unwrap_or(&0.0)
    }
    /// Nominal maximum aggregate size (mm) — first sieve where passing < 100%.
    pub fn nominal_max_size(&self) -> f64 {
        for (i, &p) in self.percent_passing.iter().enumerate() {
            if p < 100.0 {
                if i > 0 {
                    return self.sieve_sizes[i - 1];
                }
                return self.sieve_sizes[0];
            }
        }
        *self.sieve_sizes.last().unwrap_or(&0.0)
    }
}
/// Structural timber material properties.
pub struct TimberMaterial {
    /// Species grade designation.
    pub species_grade: String,
    /// Modulus of elasticity (MPa).
    pub moe: f64,
    /// Modulus of rupture (MPa).
    pub mor: f64,
    /// Longitudinal Young's modulus (MPa).
    pub e_longitudinal: f64,
    /// Radial Young's modulus (MPa).
    pub e_radial: f64,
    /// Tangential Young's modulus (MPa).
    pub e_tangential: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl TimberMaterial {
    /// Create a Douglas Fir Select Structural timber section.
    pub fn douglas_fir_ss() -> Self {
        TimberMaterial {
            species_grade: "Douglas Fir - Select Structural".to_string(),
            moe: 12_400.0,
            mor: 62.0,
            e_longitudinal: 12_400.0,
            e_radial: 930.0,
            e_tangential: 620.0,
            density: 500.0,
        }
    }
    /// Adjusted bending design value Fb' (MPa) with CD (load duration) and CM (moisture).
    pub fn adjusted_fb(&self, cd: f64, cm: f64) -> f64 {
        self.mor * cd * cm
    }
    /// Adjusted modulus E' (MPa) with CM (moisture) factor.
    pub fn adjusted_e(&self, cm: f64) -> f64 {
        self.moe * cm
    }
    /// Shear modulus G_LR ≈ E_longitudinal / 16 (approximate).
    pub fn shear_modulus_lr(&self) -> f64 {
        self.e_longitudinal / 16.0
    }
}
/// Laminated veneer lumber (LVL) section properties.
#[derive(Debug, Clone)]
pub struct LaminatedVeneerLumber {
    /// Width b (mm).
    pub b: f64,
    /// Depth d (mm).
    pub d: f64,
    /// Reference bending design value Fb (MPa).
    pub fb: f64,
    /// Reference compression parallel Fc (MPa).
    pub fc: f64,
    /// Reference tension parallel Ft (MPa).
    pub ft: f64,
    /// Reference shear Fv (MPa).
    pub fv: f64,
    /// Modulus of elasticity E (MPa).
    pub e: f64,
}
impl LaminatedVeneerLumber {
    /// Create a standard 1.9E Microllam LVL section (Weyerhaeuser product).
    pub fn microllam_1_9e(b: f64, d: f64) -> Self {
        LaminatedVeneerLumber {
            b,
            d,
            fb: 19.3,
            fc: 19.3,
            ft: 11.7,
            fv: 2.07,
            e: 13_100.0,
        }
    }
    /// Cross-sectional area (mm²).
    pub fn area(&self) -> f64 {
        self.b * self.d
    }
    /// Moment of inertia (mm⁴).
    pub fn ix(&self) -> f64 {
        self.b * self.d.powi(3) / 12.0
    }
    /// Section modulus (mm³).
    pub fn sx(&self) -> f64 {
        self.b * self.d.powi(2) / 6.0
    }
    /// Allowable bending moment (N·mm) using full Fb.
    pub fn allowable_moment(&self) -> f64 {
        self.fb * self.sx()
    }
    /// Allowable axial compression (N).
    pub fn allowable_compression(&self) -> f64 {
        self.fc * self.area()
    }
    /// Euler critical buckling load (N) for column, effective length Le (mm).
    pub fn euler_buckling_load(&self, le: f64) -> f64 {
        PI.powi(2) * self.e * self.ix() / le.powi(2)
    }
    /// Depth-to-width ratio (slenderness check, should be ≤ 5 for typical use).
    pub fn d_to_b_ratio(&self) -> f64 {
        self.d / self.b
    }
}
/// Combined geosynthetic reinforcement design for MSE walls.
#[derive(Debug, Clone)]
pub struct GeosyntheticReinforcement {
    /// Geosynthetic type: "geogrid" or "geotextile".
    pub geosyn_type: String,
    /// Allowable tensile force per unit width Ta (kN/m).
    pub ta: f64,
    /// Vertical spacing of layers sv (m).
    pub sv: f64,
    /// Friction angle of fill soil (degrees).
    pub phi_fill: f64,
    /// Unit weight of fill γ (kN/m³).
    pub gamma_fill: f64,
    /// Coverage ratio Rc.
    pub rc: f64,
}
impl GeosyntheticReinforcement {
    /// Create a geosynthetic reinforcement design.
    pub fn new(
        geosyn_type: &str,
        ta: f64,
        sv: f64,
        phi_fill: f64,
        gamma_fill: f64,
        rc: f64,
    ) -> Self {
        GeosyntheticReinforcement {
            geosyn_type: geosyn_type.to_string(),
            ta,
            sv,
            phi_fill,
            gamma_fill,
            rc,
        }
    }
    /// Maximum horizontal stress at depth z (kPa).
    pub fn horizontal_stress(&self, z: f64) -> f64 {
        let ka = ((PI / 4.0 - self.phi_fill.to_radians() / 2.0).tan()).powi(2);
        ka * self.gamma_fill * z
    }
    /// Required strength at depth z per unit width (kN/m).
    pub fn required_strength(&self, z: f64) -> f64 {
        self.horizontal_stress(z) * self.sv
    }
    /// Factor of safety in tension.
    pub fn tension_fos(&self, z: f64) -> f64 {
        let t_req = self.required_strength(z);
        if t_req < 1e-10 {
            return f64::INFINITY;
        }
        self.ta / t_req
    }
    /// Pullout length required (m) per FHWA.
    pub fn pullout_length(&self, z: f64, fos: f64) -> f64 {
        let t_max = self.required_strength(z);
        let sigma_v = self.gamma_fill * z;
        let f_star = 0.67 * self.phi_fill.to_radians().tan();
        let denom = 2.0 * self.rc * f_star * sigma_v;
        if denom < 1e-10 {
            return 1.0;
        }
        (t_max * fos / denom).max(1.0)
    }
}
/// Steel rebar bond and development length (ACI 318).
pub struct SteelRebarBond {
    /// Bar diameter db (mm).
    pub db: f64,
    /// Bar yield strength fy (MPa).
    pub fy: f64,
    /// Concrete compressive strength fc (MPa).
    pub fc: f64,
    /// Cover to center of bar (mm).
    pub cover: f64,
}
impl SteelRebarBond {
    /// Create a new rebar bond object.
    pub fn new(db: f64, fy: f64, fc: f64, cover: f64) -> Self {
        SteelRebarBond { db, fy, fc, cover }
    }
    /// Basic development length ld (mm) per ACI 318-19 (simplified).
    pub fn development_length(&self) -> f64 {
        let psi_t = 1.0;
        let psi_e = 1.0;
        let lambda = 1.0;
        (self.fy * psi_t * psi_e) / (1.1 * lambda * self.fc.sqrt()) * self.db
    }
    /// Standard 90° hook development length ldh (mm).
    pub fn hook_development_length(&self) -> f64 {
        let ldh_basic = 0.24 * self.fy * self.db / (self.fc.sqrt());
        ldh_basic.max(8.0 * self.db).max(150.0)
    }
    /// Bond stress u = fy * As / (π * db * ld) assuming uniform distribution.
    pub fn average_bond_stress(&self) -> f64 {
        let ld = self.development_length();
        self.fy / (PI * ld / self.db)
    }
}
/// Chemical admixture dosage and effect.
#[derive(Debug, Clone)]
pub struct Admixture {
    /// Type of admixture.
    pub admixture_type: AdmixtureType,
    /// Dosage (% by mass of cement).
    pub dosage: f64,
    /// Manufacturer-stated water reduction (%).
    pub water_reduction_pct: f64,
    /// Setting time modification (minutes, positive = retard, negative = accelerate).
    pub setting_time_delta: f64,
}
impl Admixture {
    /// Create a new admixture.
    pub fn new(
        admixture_type: AdmixtureType,
        dosage: f64,
        water_reduction_pct: f64,
        setting_time_delta: f64,
    ) -> Self {
        Admixture {
            admixture_type,
            dosage,
            water_reduction_pct,
            setting_time_delta,
        }
    }
    /// Effective water content after admixture (kg/m³).
    pub fn adjusted_water(&self, base_water: f64) -> f64 {
        base_water * (1.0 - self.water_reduction_pct / 100.0)
    }
    /// Effective w/c ratio after water reduction.
    pub fn adjusted_wc_ratio(&self, base_wc: f64) -> f64 {
        base_wc * (1.0 - self.water_reduction_pct / 100.0)
    }
    /// Is this admixture a set accelerator?
    pub fn is_accelerator(&self) -> bool {
        self.admixture_type == AdmixtureType::Accelerator
    }
}
/// Cross-laminated timber (CLT) panel properties.
#[derive(Debug, Clone)]
pub struct CrossLaminatedTimber {
    /// Panel width (mm).
    pub width: f64,
    /// Panel total thickness (mm).
    pub thickness: f64,
    /// Number of layers (odd number).
    pub n_layers: u32,
    /// Layer thickness (mm), all equal assumed.
    pub layer_thickness: f64,
    /// Parallel-to-grain MOE E0 (MPa).
    pub e0: f64,
    /// Perpendicular-to-grain MOE E90 (MPa).
    pub e90: f64,
    /// Bending strength Fb (MPa) for parallel layers.
    pub fb: f64,
    /// Rolling shear modulus Grt (MPa).
    pub g_rolling: f64,
}
impl CrossLaminatedTimber {
    /// Create a 5-layer CLT panel with standard properties.
    pub fn five_layer(width: f64, total_thickness: f64) -> Self {
        let n_layers = 5u32;
        let layer_t = total_thickness / n_layers as f64;
        CrossLaminatedTimber {
            width,
            thickness: total_thickness,
            n_layers,
            layer_thickness: layer_t,
            e0: 11_000.0,
            e90: 370.0,
            fb: 24.0,
            g_rolling: 65.0,
        }
    }
    /// Effective bending stiffness EIeff (N·mm²) per unit width using gamma method.
    pub fn effective_bending_stiffness(&self) -> f64 {
        let n_par = self.n_layers.div_ceil(2);
        let h = self.layer_thickness;
        let mut ei = 0.0;
        for i in 0..n_par {
            let yi = (self.thickness / 2.0) - h / 2.0 - i as f64 * 2.0 * h;
            ei += self.e0 * self.width * h.powi(3) / 12.0 + self.e0 * self.width * h * yi.powi(2);
        }
        ei
    }
    /// Effective axial stiffness EAeff (N/mm).
    pub fn effective_axial_stiffness(&self) -> f64 {
        let n_par = self.n_layers.div_ceil(2);
        self.e0 * self.width * self.layer_thickness * n_par as f64
    }
    /// Rolling shear stress capacity check τ_max (MPa) per unit width.
    pub fn rolling_shear_capacity(&self) -> f64 {
        self.g_rolling * self.layer_thickness / self.width
    }
    /// Panel area (mm²/m width).
    pub fn area_per_m(&self) -> f64 {
        self.thickness * 1000.0
    }
}
/// Masonry shear wall capacity (TMS 402 / MSJC).
#[derive(Debug, Clone)]
pub struct MasonryShearWall {
    /// Wall length (mm).
    pub length: f64,
    /// Wall thickness (mm).
    pub thickness: f64,
    /// Wall height (mm).
    pub height: f64,
    /// f'm — masonry prism strength (MPa).
    pub fm: f64,
    /// Area of vertical steel (mm²/m length).
    pub av: f64,
    /// Steel yield strength fy (MPa).
    pub fy: f64,
    /// Normal axial stress from vertical loads (MPa).
    pub sigma_v: f64,
}
impl MasonryShearWall {
    /// Create a masonry shear wall.
    pub fn new(
        length: f64,
        thickness: f64,
        height: f64,
        fm: f64,
        av: f64,
        fy: f64,
        sigma_v: f64,
    ) -> Self {
        MasonryShearWall {
            length,
            thickness,
            height,
            fm,
            av,
            fy,
            sigma_v,
        }
    }
    /// Net shear area An (mm²).
    pub fn net_area(&self) -> f64 {
        self.length * self.thickness
    }
    /// In-plane shear capacity Vn (N) per TMS 402 (simplified).
    pub fn in_plane_shear_capacity(&self) -> f64 {
        let m_v_d = (self.height / self.length).min(1.0);
        let an = self.net_area();
        let p = self.sigma_v * an;
        let vnm = (4.0 - 1.75 * m_v_d) * an * self.fm.sqrt() + 0.25 * p;
        let vns = 0.5 * self.av * self.fy * self.length / 1000.0;
        (vnm + vns).min(6.0 * an * self.fm.sqrt())
    }
    /// Out-of-plane flexural capacity Mn (N·mm) per unit height.
    pub fn out_of_plane_moment_capacity(&self) -> f64 {
        let d = self.thickness / 2.0;
        self.av * self.fy * d / 1000.0
    }
    /// Aspect ratio h/l (shear wall slenderness).
    pub fn aspect_ratio(&self) -> f64 {
        self.height / self.length
    }
}
/// Driven pile foundation capacity.
#[derive(Debug, Clone)]
pub struct PileFoundation {
    /// Pile diameter or width (m).
    pub diameter: f64,
    /// Pile length L (m).
    pub length: f64,
    /// Pile type: "concrete", "steel", "timber".
    pub pile_type: String,
    /// Unit skin friction qs (kPa) along pile shaft.
    pub qs: f64,
    /// Unit end bearing qb (kPa) at pile tip.
    pub qb: f64,
    /// Soil cohesion cu (kPa).
    pub cu: f64,
    /// Adhesion factor α.
    pub alpha: f64,
}
impl PileFoundation {
    /// Create a concrete driven pile.
    pub fn concrete_pile(diameter: f64, length: f64, qs: f64, qb: f64) -> Self {
        PileFoundation {
            diameter,
            length,
            pile_type: "concrete".to_string(),
            qs,
            qb,
            cu: qs,
            alpha: 0.5,
        }
    }
    /// Pile perimeter (m).
    pub fn perimeter(&self) -> f64 {
        PI * self.diameter
    }
    /// Pile tip area (m²).
    pub fn tip_area(&self) -> f64 {
        PI * (self.diameter / 2.0).powi(2)
    }
    /// Ultimate skin friction (kN).
    pub fn skin_friction(&self) -> f64 {
        self.qs * self.perimeter() * self.length
    }
    /// Ultimate end bearing (kN).
    pub fn end_bearing(&self) -> f64 {
        self.qb * self.tip_area()
    }
    /// Ultimate pile capacity Qu (kN).
    pub fn ultimate_capacity(&self) -> f64 {
        self.skin_friction() + self.end_bearing()
    }
    /// Allowable pile capacity Qa (kN) with FOS = 2.5.
    pub fn allowable_capacity(&self) -> f64 {
        self.ultimate_capacity() / 2.5
    }
    /// α-method skin friction for cohesive soils (kN).
    pub fn alpha_skin_friction(&self) -> f64 {
        self.alpha * self.cu * self.perimeter() * self.length
    }
    /// Settlement under working load (mm) — elastic compression of pile.
    pub fn elastic_compression(&self, load_kn: f64) -> f64 {
        let e_pile = if self.pile_type == "concrete" {
            25_000.0
        } else {
            200_000.0
        };
        let area_mm2 = self.tip_area() * 1e6;
        let load_n = load_kn * 1000.0;
        let length_mm = self.length * 1000.0;
        load_n * length_mm / (e_pile * area_mm2)
    }
}
/// Concrete material properties per ACI/EN standard.
pub struct ConcreteMaterial {
    /// Characteristic compressive strength f'c (MPa).
    pub fc: f64,
    /// Tensile strength ft (MPa), typically ≈ 0.1 * fc.
    pub ft: f64,
    /// Modulus of elasticity Ec (MPa).
    pub ec: f64,
    /// Poisson's ratio ν.
    pub poisson: f64,
    /// Drying shrinkage strain ε_sh (dimensionless).
    pub shrinkage_strain: f64,
    /// Creep coefficient φ (ratio of creep strain to elastic strain).
    pub creep_coefficient: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl ConcreteMaterial {
    /// Create a normal-weight concrete with characteristic strength `fc` (MPa).
    ///
    /// Ec is estimated by ACI 318: Ec = 4700 * sqrt(fc) (MPa).
    pub fn new(fc: f64) -> Self {
        let ec = 4700.0 * fc.sqrt();
        let ft = 0.1 * fc;
        ConcreteMaterial {
            fc,
            ft,
            ec,
            poisson: 0.2,
            shrinkage_strain: 3e-4,
            creep_coefficient: 2.0,
            density: 2400.0,
        }
    }
    /// Shear modulus G = Ec / (2*(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.ec / (2.0 * (1.0 + self.poisson))
    }
    /// Splitting tensile strength (ACI): fct = 0.56 * sqrt(fc) MPa.
    pub fn splitting_tensile_strength(&self) -> f64 {
        0.56 * self.fc.sqrt()
    }
    /// Ultimate compressive strain (ACI): 0.003.
    pub fn ultimate_compressive_strain(&self) -> f64 {
        0.003
    }
    /// Modulus of rupture (flexural): fr = 0.62 * sqrt(fc) MPa.
    pub fn modulus_of_rupture(&self) -> f64 {
        0.62 * self.fc.sqrt()
    }
}
/// Reinforced concrete section for moment/shear capacity (ACI 318).
pub struct ReinforcedConcrete {
    /// Concrete material.
    pub concrete: ConcreteMaterial,
    /// Width of the compression block (mm).
    pub b: f64,
    /// Effective depth to tension steel (mm).
    pub d: f64,
    /// Area of tension steel (mm²).
    pub as_t: f64,
    /// Yield strength of steel fy (MPa).
    pub fy: f64,
    /// Area of compression steel (mm²), if any.
    pub as_comp: f64,
    /// Depth to compression steel (mm).
    pub d_prime: f64,
}
impl ReinforcedConcrete {
    /// Create a new reinforced concrete section.
    pub fn new(
        concrete: ConcreteMaterial,
        b: f64,
        d: f64,
        as_t: f64,
        fy: f64,
        as_comp: f64,
        d_prime: f64,
    ) -> Self {
        ReinforcedConcrete {
            concrete,
            b,
            d,
            as_t,
            fy,
            as_comp,
            d_prime,
        }
    }
    /// Nominal moment capacity Mn (N·mm) using ACI rectangular stress block.
    pub fn moment_capacity(&self) -> f64 {
        let fc = self.concrete.fc;
        let beta1 = if fc <= 28.0 {
            0.85
        } else {
            (0.85 - 0.05 * (fc - 28.0) / 7.0).max(0.65)
        };
        let a = self.as_t * self.fy / (0.85 * fc * self.b);
        let c = a / beta1;
        let eps_cu = 0.003;
        let eps_prime = eps_cu * (c - self.d_prime) / c;
        let fs_prime = if eps_prime >= self.fy / 200_000.0 {
            self.fy
        } else {
            eps_prime * 200_000.0
        };

        self.as_t * self.fy * (self.d - a / 2.0) + self.as_comp * fs_prime * (self.d - self.d_prime)
    }
    /// Design moment capacity φMn (φ = 0.9 for tension-controlled sections).
    pub fn design_moment_capacity(&self) -> f64 {
        0.9 * self.moment_capacity()
    }
    /// Nominal shear capacity Vn (N) per ACI 318: Vn = Vc + Vs.
    /// Concrete shear: Vc = 0.17 * sqrt(fc) * b * d.
    pub fn shear_capacity(&self) -> f64 {
        0.17 * self.concrete.fc.sqrt() * self.b * self.d
    }
    /// Reinforcement ratio ρ = As / (b * d).
    pub fn reinforcement_ratio(&self) -> f64 {
        self.as_t / (self.b * self.d)
    }
    /// Minimum steel ratio ρ_min per ACI 318.
    pub fn rho_min(&self) -> f64 {
        let term1 = 0.25 * self.concrete.fc.sqrt() / self.fy;
        let term2 = 1.4 / self.fy;
        term1.max(term2)
    }
}
/// Hollow Structural Section (HSS / RHS) properties.
#[derive(Debug, Clone)]
pub struct HssSection {
    /// Outer width B (mm).
    pub b: f64,
    /// Outer height H (mm).
    pub h: f64,
    /// Wall thickness t (mm).
    pub t: f64,
    /// Yield strength Fy (MPa).
    pub fy: f64,
    /// Modulus of elasticity E (MPa).
    pub e: f64,
}
impl HssSection {
    /// Create a rectangular HSS section.
    pub fn rectangular(b: f64, h: f64, t: f64, fy: f64) -> Self {
        HssSection {
            b,
            h,
            t,
            fy,
            e: 200_000.0,
        }
    }
    /// Create a square HSS section.
    pub fn square(b: f64, t: f64, fy: f64) -> Self {
        HssSection {
            b,
            h: b,
            t,
            fy,
            e: 200_000.0,
        }
    }
    /// Cross-sectional area (mm²).
    pub fn area(&self) -> f64 {
        self.b * self.h - (self.b - 2.0 * self.t) * (self.h - 2.0 * self.t)
    }
    /// Moment of inertia about strong axis Ix (mm⁴).
    pub fn ix(&self) -> f64 {
        (self.b * self.h.powi(3) - (self.b - 2.0 * self.t) * (self.h - 2.0 * self.t).powi(3)) / 12.0
    }
    /// Moment of inertia about weak axis Iy (mm⁴).
    pub fn iy(&self) -> f64 {
        (self.h * self.b.powi(3) - (self.h - 2.0 * self.t) * (self.b - 2.0 * self.t).powi(3)) / 12.0
    }
    /// Elastic section modulus Sx (mm³).
    pub fn sx(&self) -> f64 {
        self.ix() / (self.h / 2.0)
    }
    /// Plastic section modulus Zx (mm³) (approximate for rectangular HSS).
    pub fn zx(&self) -> f64 {
        let outer = self.b * self.h.powi(2) / 4.0;
        let inner = (self.b - 2.0 * self.t) * (self.h - 2.0 * self.t).powi(2) / 4.0;
        outer - inner
    }
    /// Nominal moment capacity Mp = Zx * Fy (N·mm).
    pub fn plastic_moment(&self) -> f64 {
        self.zx() * self.fy
    }
    /// Torsional constant J (mm⁴) for closed section.
    pub fn torsional_constant(&self) -> f64 {
        let h_mid = self.h - self.t;
        let b_mid = self.b - self.t;
        2.0 * self.t * b_mid * h_mid * (b_mid * h_mid) / (b_mid + h_mid)
    }
    /// Warping constant Cw ≈ 0 for closed HSS (negligible).
    pub fn warping_constant(&self) -> f64 {
        0.0
    }
    /// Slenderness ratio for local buckling (web): h/t.
    pub fn web_slenderness(&self) -> f64 {
        (self.h - 2.0 * self.t) / self.t
    }
    /// Slenderness ratio for local buckling (flange): b/t.
    pub fn flange_slenderness(&self) -> f64 {
        (self.b - 2.0 * self.t) / self.t
    }
    /// AISC compact limit for HSS flanges: λp = 1.12 * sqrt(E/Fy).
    pub fn compact_flange_limit(&self) -> f64 {
        1.12 * (self.e / self.fy).sqrt()
    }
    /// Check if flange is compact.
    pub fn flange_is_compact(&self) -> bool {
        self.flange_slenderness() <= self.compact_flange_limit()
    }
}
/// Eurocode EN 1990 load combinations (fundamental combination, ULS).
#[derive(Debug, Clone)]
pub struct EurocodeLoad {
    /// Permanent action Gk (kN or kN/m²).
    pub gk: f64,
    /// Variable leading action Qk1 (kN or kN/m²).
    pub qk1: f64,
    /// Variable accompanying action Qk2 (kN or kN/m²).
    pub qk2: f64,
    /// Wind action Wk (kN or kN/m²).
    pub wk: f64,
    /// Snow action Sk (kN or kN/m²).
    pub sk: f64,
    /// Seismic action Ed (kN or kN/m²).
    pub ed: f64,
}
impl EurocodeLoad {
    /// Create a Eurocode load set.
    pub fn new(gk: f64, qk1: f64, qk2: f64, wk: f64, sk: f64, ed: f64) -> Self {
        EurocodeLoad {
            gk,
            qk1,
            qk2,
            wk,
            sk,
            ed,
        }
    }
    /// STR/GEO ULS fundamental combination (Eq. 6.10): γG*Gk + γQ1*Qk1 + Σγ_Qi*ψ0i*Qki.
    /// γG = 1.35, γQ = 1.50, ψ0 = 0.7 for imposed, 0.6 for wind/snow.
    pub fn uls_combo_610(&self) -> f64 {
        1.35 * self.gk + 1.50 * self.qk1 + 1.50 * 0.7 * self.qk2
    }
    /// ULS Eq. 6.10a (γG*Gk + 1.5*ψ0*Qk1): generally less critical.
    pub fn uls_combo_610a(&self) -> f64 {
        1.35 * self.gk + 1.50 * 0.7 * self.qk1
    }
    /// ULS Eq. 6.10b (ξ*γG*Gk + 1.5*Qk1): ξ = 0.85 reduction on permanent action.
    pub fn uls_combo_610b(&self) -> f64 {
        0.85 * 1.35 * self.gk + 1.50 * self.qk1
    }
    /// ULS accidental combination with seismic: Gk + Ed + ψ2*Qk.
    /// ψ2 = 0.3 for residential.
    pub fn uls_seismic(&self) -> f64 {
        self.gk + self.ed + 0.3 * self.qk1
    }
    /// Characteristic SLS combination: Gk + Qk1 + ψ0*Qk2.
    pub fn sls_characteristic(&self) -> f64 {
        self.gk + self.qk1 + 0.7 * self.qk2
    }
    /// Frequent SLS combination: Gk + ψ1*Qk1 + ψ2*Qk2. ψ1=0.5, ψ2=0.3.
    pub fn sls_frequent(&self) -> f64 {
        self.gk + 0.5 * self.qk1 + 0.3 * self.qk2
    }
    /// Quasi-permanent SLS combination: Gk + ψ2*Qk1. ψ2=0.3.
    pub fn sls_quasi_permanent(&self) -> f64 {
        self.gk + 0.3 * self.qk1
    }
    /// Governing ULS combination.
    pub fn governing_uls(&self) -> f64 {
        [
            self.uls_combo_610(),
            self.uls_combo_610a(),
            self.uls_combo_610b(),
        ]
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Add wind to ULS: Gk + Wk combination (uplift check).
    pub fn uls_wind_uplift(&self) -> f64 {
        0.9 * self.gk + 1.50 * self.wk
    }
    /// Snow load ULS: 1.35*Gk + 1.5*Sk.
    pub fn uls_snow(&self) -> f64 {
        1.35 * self.gk + 1.50 * self.sk
    }
}
/// C-channel (American Standard Channel) steel section properties.
#[derive(Debug, Clone)]
pub struct ChannelSection {
    /// Overall depth d (mm).
    pub d: f64,
    /// Flange width bf (mm).
    pub bf: f64,
    /// Flange thickness tf (mm).
    pub tf: f64,
    /// Web thickness tw (mm).
    pub tw: f64,
    /// Yield strength Fy (MPa).
    pub fy: f64,
}
impl ChannelSection {
    /// Create a C-channel section.
    pub fn new(d: f64, bf: f64, tf: f64, tw: f64, fy: f64) -> Self {
        ChannelSection { d, bf, tf, tw, fy }
    }
    /// Cross-sectional area (mm²).
    pub fn area(&self) -> f64 {
        2.0 * self.bf * self.tf + (self.d - 2.0 * self.tf) * self.tw
    }
    /// Shear center location ex (mm) from web face.
    pub fn shear_center_x(&self) -> f64 {
        let bf = self.bf;
        let tf = self.tf;
        let d = self.d;
        let tw = self.tw;
        let hw = d - 2.0 * tf;
        let ix = self.ix();
        let sx = if d > 0.0 { ix / (d / 2.0) } else { 1.0 };
        bf.powi(2) * tf / (2.0 * sx / d) / (1.0 + (hw * tw) / (6.0 * bf * tf))
    }
    /// Moment of inertia Ix (mm⁴).
    pub fn ix(&self) -> f64 {
        let hw = self.d - 2.0 * self.tf;
        self.bf * self.d.powi(3) / 12.0 - (self.bf - self.tw) * hw.powi(3) / 12.0
    }
    /// Elastic section modulus Sx (mm³).
    pub fn sx(&self) -> f64 {
        self.ix() / (self.d / 2.0)
    }
    /// Nominal moment capacity Mn (N·mm).
    pub fn moment_capacity(&self) -> f64 {
        self.sx() * self.fy
    }
    /// Shear capacity (simplified) Vn (N).
    pub fn shear_capacity(&self) -> f64 {
        let hw = self.d - 2.0 * self.tf;
        0.6 * self.fy * hw * self.tw
    }
}
/// Type of chemical admixture.
#[derive(Debug, Clone, PartialEq)]
pub enum AdmixtureType {
    /// Water reducer / plasticiser.
    WaterReducer,
    /// Superplasticiser (high-range water reducer).
    Superplasticiser,
    /// Accelerator (increases early strength).
    Accelerator,
    /// Retarder (extends workability time).
    Retarder,
    /// Air entrainer.
    AirEntrainer,
    /// Shrinkage-reducing admixture.
    ShrinkageReducer,
    /// Corrosion inhibitor.
    CorrosionInhibitor,
}
/// Moisture content effects on timber mechanical properties (NDS Table 4A).
#[derive(Debug, Clone)]
pub struct MoistureCorrectionTimber {
    /// Moisture content MC (%) in service.
    pub mc_service: f64,
    /// Fiber saturation point MC_fsp (%), typically 28–30%.
    pub mc_fsp: f64,
    /// Green MOE (MPa) at MC = MC_fsp.
    pub e_green: f64,
    /// Green Fb (MPa).
    pub fb_green: f64,
}
impl MoistureCorrectionTimber {
    /// Create a moisture correction model for Douglas Fir.
    pub fn douglas_fir() -> Self {
        MoistureCorrectionTimber {
            mc_service: 15.0,
            mc_fsp: 28.0,
            e_green: 11_000.0,
            fb_green: 40.0,
        }
    }
    /// Compute CM factor for modulus of elasticity (NDS Table 4A).
    pub fn cm_factor_e(&self) -> f64 {
        if self.mc_service >= self.mc_fsp {
            return 1.0;
        }
        if self.mc_service <= 19.0 {
            1.0
        } else {
            1.0 - 0.1 * (self.mc_service - 19.0) / (self.mc_fsp - 19.0)
        }
    }
    /// Compute CM factor for bending strength Fb (NDS Table 4A).
    pub fn cm_factor_fb(&self) -> f64 {
        if self.mc_service >= self.mc_fsp {
            return 0.85;
        }
        if self.mc_service <= 19.0 {
            1.0
        } else {
            0.85 + 0.15 * (self.mc_fsp - self.mc_service) / (self.mc_fsp - 19.0)
        }
    }
    /// Adjusted MOE (MPa) at service MC.
    pub fn adjusted_e(&self) -> f64 {
        self.e_green * self.cm_factor_e()
    }
    /// Adjusted Fb (MPa) at service MC.
    pub fn adjusted_fb(&self) -> f64 {
        self.fb_green * self.cm_factor_fb()
    }
    /// Shrinkage coefficient: % shrinkage per % MC change (tangential direction).
    pub fn tangential_shrinkage_per_pct_mc(&self) -> f64 {
        0.29
    }
    /// Dimensional change for MC change from green to service.
    pub fn dimensional_change(&self, dimension_mm: f64) -> f64 {
        let delta_mc = (self.mc_fsp - self.mc_service).max(0.0);
        dimension_mm * self.tangential_shrinkage_per_pct_mc() * delta_mc / 100.0
    }
}
/// Geotextile filter design for drainage and erosion control.
#[derive(Debug, Clone)]
pub struct GeotextileFilter {
    /// Apparent opening size AOS (O95, mm).
    pub aos: f64,
    /// Permittivity ψ (s⁻¹).
    pub permittivity: f64,
    /// Transmissivity θ (m²/s).
    pub transmissivity: f64,
    /// Tensile strength (kN/m).
    pub tensile_strength: f64,
    /// Elongation at break (%).
    pub elongation: f64,
    /// Soil D85 particle size (mm).
    pub soil_d85: f64,
    /// Soil D15 particle size (mm).
    pub soil_d15: f64,
    /// In-plane hydraulic conductivity kp (m/s).
    pub kp: f64,
    /// Normal hydraulic conductivity kn (m/s).
    pub kn: f64,
}
impl GeotextileFilter {
    /// Create a typical woven geotextile for road subbase drainage.
    pub fn woven_road_drainage() -> Self {
        GeotextileFilter {
            aos: 0.212,
            permittivity: 0.5,
            transmissivity: 1e-4,
            tensile_strength: 40.0,
            elongation: 15.0,
            soil_d85: 0.3,
            soil_d15: 0.05,
            kp: 1e-3,
            kn: 1e-4,
        }
    }
    /// Check retention criterion: O95 ≤ B * D85.
    /// B = 1.0–2.0 depending on soil uniformity.
    pub fn retention_criterion_met(&self, b: f64) -> bool {
        self.aos <= b * self.soil_d85
    }
    /// Check permeability criterion: kn ≥ ksoil.
    pub fn permeability_criterion_met(&self, k_soil: f64) -> bool {
        self.kn >= k_soil
    }
    /// Gradient ratio test: index of clogging potential (target ≤ 3).
    pub fn gradient_ratio(&self) -> f64 {
        if self.soil_d15 < 1e-10 {
            return f64::INFINITY;
        }
        self.aos / self.soil_d15
    }
    /// Filter ratio (coarse side): AOS / D85 should be ≤ 2.0 for woven.
    pub fn filter_ratio(&self) -> f64 {
        if self.soil_d85 < 1e-10 {
            return f64::INFINITY;
        }
        self.aos / self.soil_d85
    }
    /// Required width for overlapping in a trench drain (m).
    pub fn overlap_width(&self, trench_depth_m: f64) -> f64 {
        (0.3_f64).max(trench_depth_m * 0.5)
    }
}
/// Geogrid properties for soil reinforcement.
#[derive(Debug, Clone)]
pub struct Geogrid {
    /// Ultimate tensile strength in machine direction (kN/m).
    pub tult_md: f64,
    /// Ultimate tensile strength in cross-machine direction (kN/m).
    pub tult_cmd: f64,
    /// Junction efficiency (%).
    pub junction_efficiency: f64,
    /// Aperture size (mm).
    pub aperture_size: f64,
    /// Long-term design strength (kN/m) after creep reduction.
    pub ltds: f64,
    /// Coverage ratio (ratio of solid area to total area).
    pub coverage_ratio: f64,
}
impl Geogrid {
    /// Create a typical biaxial polypropylene geogrid (BX-1100).
    pub fn bx1100() -> Self {
        Geogrid {
            tult_md: 15.0,
            tult_cmd: 20.0,
            junction_efficiency: 93.0,
            aperture_size: 33.0,
            ltds: 8.0,
            coverage_ratio: 0.35,
        }
    }
    /// Create a uniaxial HDPE geogrid (UX-1500HS).
    pub fn ux1500hs() -> Self {
        Geogrid {
            tult_md: 150.0,
            tult_cmd: 25.0,
            junction_efficiency: 91.0,
            aperture_size: 16.0,
            ltds: 80.0,
            coverage_ratio: 0.70,
        }
    }
    /// Reduction factor for installation damage RF_ID (typically 1.1–1.4).
    /// Long-term allowable strength = LTDS / RF_ID / RF_creep.
    pub fn allowable_strength(&self, rf_id: f64, rf_cr: f64) -> f64 {
        self.ltds / (rf_id * rf_cr)
    }
    /// Interaction coefficient for soil-geogrid friction Ci.
    /// `phi_soil` = soil friction angle (degrees).
    pub fn interaction_coefficient(&self, phi_soil: f64) -> f64 {
        let base_ci = 0.8;
        let phi = phi_soil.to_radians();
        base_ci * phi.tan() / phi.tan()
    }
    /// Passive resistance contribution τ_p (kPa) per geogrid layer.
    pub fn passive_resistance(&self, sigma_v: f64) -> f64 {
        self.ltds * sigma_v / 100.0
    }
}
/// Structural load combinations and load factors (ASCE 7 LRFD).
#[derive(Debug, Clone)]
pub struct StructuralLoad {
    /// Dead load D (kN or kN/m²).
    pub dead: f64,
    /// Live load L (kN or kN/m²).
    pub live: f64,
    /// Wind load W (kN or kN/m²).
    pub wind: f64,
    /// Seismic load E (kN or kN/m²).
    pub seismic: f64,
    /// Snow load S (kN or kN/m²).
    pub snow: f64,
}
impl StructuralLoad {
    /// Create a new structural load set.
    pub fn new(dead: f64, live: f64, wind: f64, seismic: f64, snow: f64) -> Self {
        StructuralLoad {
            dead,
            live,
            wind,
            seismic,
            snow,
        }
    }
    /// LRFD load combination 1: 1.4D.
    pub fn lrfd_combo1(&self) -> f64 {
        1.4 * self.dead
    }
    /// LRFD load combination 2: 1.2D + 1.6L + 0.5S.
    pub fn lrfd_combo2(&self) -> f64 {
        1.2 * self.dead + 1.6 * self.live + 0.5 * self.snow
    }
    /// LRFD load combination 3: 1.2D + 1.0W + 1.0L + 0.5S.
    pub fn lrfd_combo3(&self) -> f64 {
        1.2 * self.dead + 1.0 * self.wind + 1.0 * self.live + 0.5 * self.snow
    }
    /// LRFD load combination 4: 0.9D + 1.0W (overturning check).
    pub fn lrfd_combo4(&self) -> f64 {
        0.9 * self.dead + 1.0 * self.wind
    }
    /// LRFD seismic combination: 1.2D + 1.0E + 1.0L + 0.2S.
    pub fn lrfd_seismic(&self) -> f64 {
        1.2 * self.dead + 1.0 * self.seismic + 1.0 * self.live + 0.2 * self.snow
    }
    /// ASD service load combination: D + L.
    pub fn asd_combo_dl(&self) -> f64 {
        self.dead + self.live
    }
    /// Governing (maximum) LRFD combination.
    pub fn governing_lrfd(&self) -> f64 {
        [
            self.lrfd_combo1(),
            self.lrfd_combo2(),
            self.lrfd_combo3(),
            self.lrfd_combo4(),
            self.lrfd_seismic(),
        ]
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
    }
}
/// Asphalt mixture properties (Marshall mix design).
pub struct AsphaltMixture {
    /// Bitumen content (% by mass of mix).
    pub bitumen_content: f64,
    /// Voids in mineral aggregate VMA (%).
    pub vma: f64,
    /// Voids filled with asphalt VFA (%).
    pub vfa: f64,
    /// Air voids (%).
    pub air_voids: f64,
    /// Dynamic stability (rutting resistance) in passes/mm.
    pub dynamic_stability: f64,
    /// Marshall stability (N).
    pub stability: f64,
    /// Marshall flow (mm).
    pub flow: f64,
}
impl AsphaltMixture {
    /// Create a standard Superpave-12.5 mix design.
    pub fn superpave_12_5() -> Self {
        AsphaltMixture {
            bitumen_content: 5.0,
            vma: 14.0,
            vfa: 72.0,
            air_voids: 4.0,
            dynamic_stability: 1000.0,
            stability: 8000.0,
            flow: 3.0,
        }
    }
    /// Computed bulk density (Gmb) approximation: (100 - AV%) / 100 * Gmm.
    /// `gmm` = theoretical maximum density.
    pub fn bulk_density(&self, gmm: f64) -> f64 {
        (100.0 - self.air_voids) / 100.0 * gmm
    }
    /// Determine if design meets Superpave 4% air voids criterion.
    pub fn meets_air_voids_criterion(&self) -> bool {
        (self.air_voids - 4.0).abs() < 0.5
    }
}
/// Wide-flange (I-beam) or channel steel section properties.
#[derive(Debug, Clone)]
pub struct SteelSection {
    /// Cross-sectional area A (mm²).
    pub area: f64,
    /// Moment of inertia about strong axis Ix (mm⁴).
    pub ix: f64,
    /// Moment of inertia about weak axis Iy (mm⁴).
    pub iy: f64,
    /// Plastic section modulus about strong axis Zx (mm³).
    pub zx: f64,
    /// Plastic section modulus about weak axis Zy (mm³).
    pub zy: f64,
    /// Yield strength Fy (MPa).
    pub fy: f64,
    /// Modulus of elasticity E (MPa).
    pub e: f64,
}
impl SteelSection {
    /// Create an I-section from flange/web dimensions.
    ///
    /// `bf` = flange width, `tf` = flange thickness, `d` = total depth,
    /// `tw` = web thickness, all in mm.
    pub fn i_section(bf: f64, tf: f64, d: f64, tw: f64, fy: f64) -> Self {
        let hw = d - 2.0 * tf;
        let area = 2.0 * bf * tf + hw * tw;
        let ix = bf * d.powi(3) / 12.0 - (bf - tw) * hw.powi(3) / 12.0;
        let iy = 2.0 * tf * bf.powi(3) / 12.0 + hw * tw.powi(3) / 12.0;
        let zx = bf * tf * (d / 2.0 - tf / 2.0) * 2.0 + tw * hw.powi(2) / 4.0;
        let zy = bf.powi(2) * tf / 2.0 + tw.powi(2) * hw / 4.0;
        SteelSection {
            area,
            ix,
            iy,
            zx,
            zy,
            fy,
            e: 200_000.0,
        }
    }
    /// Elastic section modulus Sx = Ix / (d/2).
    pub fn sx(&self, d: f64) -> f64 {
        self.ix / (d / 2.0)
    }
    /// Plastic moment capacity Mp = Zx * Fy (N·mm).
    pub fn plastic_moment(&self) -> f64 {
        self.zx * self.fy
    }
    /// Axial load capacity Pn = A * Fy (N) for short columns (no buckling).
    pub fn axial_capacity(&self) -> f64 {
        self.area * self.fy
    }
    /// Radius of gyration about strong axis rx = sqrt(Ix/A).
    pub fn rx(&self) -> f64 {
        (self.ix / self.area).sqrt()
    }
    /// Radius of gyration about weak axis ry = sqrt(Iy/A).
    pub fn ry(&self) -> f64 {
        (self.iy / self.area).sqrt()
    }
    /// Critical buckling stress (Euler) for slenderness ratio KL/r.
    pub fn elastic_buckling_stress(&self, kl_over_r: f64) -> f64 {
        PI.powi(2) * self.e / kl_over_r.powi(2)
    }
}
/// Foundation soil properties for bearing capacity analysis.
pub struct FoundationSoil {
    /// Undrained shear strength cu (kPa) for cohesive soils.
    pub cu: f64,
    /// Effective internal friction angle φ' (degrees).
    pub phi_deg: f64,
    /// Cohesion c' (kPa) for c-φ soil.
    pub c_prime: f64,
    /// Soil unit weight γ (kN/m³).
    pub gamma: f64,
    /// SPT N-value (blows per 300mm).
    pub spt_n: u32,
    /// Constrained modulus D (MPa) for settlement.
    pub constrained_modulus: f64,
}
impl FoundationSoil {
    /// Create a typical medium stiff clay.
    pub fn medium_clay() -> Self {
        FoundationSoil {
            cu: 50.0,
            phi_deg: 0.0,
            c_prime: 50.0,
            gamma: 18.0,
            spt_n: 10,
            constrained_modulus: 5.0,
        }
    }
    /// Ultimate bearing capacity by Terzaghi (strip footing, general shear).
    ///
    /// `b` = footing width (m), `df` = depth of foundation (m).
    pub fn ultimate_bearing_capacity(&self, b: f64, df: f64) -> f64 {
        let phi = self.phi_deg.to_radians();
        let (nc, nq, ng) = terzaghi_bearing_factors(phi);
        self.c_prime * nc + self.gamma * df * nq + 0.5 * self.gamma * b * ng
    }
    /// Allowable bearing capacity with factor of safety FOS = 3.
    pub fn allowable_bearing_capacity(&self, b: f64, df: f64) -> f64 {
        self.ultimate_bearing_capacity(b, df) / 3.0
    }
    /// Immediate settlement by Boussinesq (elastic, uniform load).
    ///
    /// `q` = net foundation pressure (kPa), `b` = footing width (m),
    /// `es` = Young's modulus of soil (MPa), `nu` = Poisson's ratio.
    pub fn immediate_settlement(&self, q: f64, b: f64, es_mpa: f64, nu: f64) -> f64 {
        let is = 0.82;
        q * b * (1.0 - nu * nu) * is / (es_mpa * 1000.0)
    }
}
/// Pre-stressed concrete section analysis (pre-tensioned and post-tensioned).
#[derive(Debug, Clone)]
pub struct PrestressedConcrete {
    /// Gross cross-section area (mm²).
    pub ag: f64,
    /// Moment of inertia of gross section (mm⁴).
    pub ig: f64,
    /// Distance from centroid to bottom fiber (mm).
    pub yb: f64,
    /// Distance from centroid to top fiber (mm).
    pub yt: f64,
    /// Prestressing steel area (mm²).
    pub aps: f64,
    /// Ultimate strength of prestressing steel fpu (MPa).
    pub fpu: f64,
    /// Initial prestress (MPa) after jacking.
    pub fpi: f64,
    /// Eccentricity of prestress at midspan (mm).
    pub eccentricity: f64,
    /// Concrete f'c (MPa).
    pub fc: f64,
    /// Span length (mm).
    pub span: f64,
}
impl PrestressedConcrete {
    /// Create a new pre-stressed concrete section.
    pub fn new(
        ag: f64,
        ig: f64,
        yb: f64,
        yt: f64,
        aps: f64,
        fpu: f64,
        fpi: f64,
        eccentricity: f64,
        fc: f64,
        span: f64,
    ) -> Self {
        PrestressedConcrete {
            ag,
            ig,
            yb,
            yt,
            aps,
            fpu,
            fpi,
            eccentricity,
            fc,
            span,
        }
    }
    /// Initial prestress force Pi (N).
    pub fn initial_prestress_force(&self) -> f64 {
        self.aps * self.fpi
    }
    /// Effective prestress after losses (using 20% total loss estimate).
    pub fn effective_prestress_force(&self, loss_fraction: f64) -> f64 {
        self.initial_prestress_force() * (1.0 - loss_fraction)
    }
    /// Kern distance (upper and lower kern points) for no-tension design.
    pub fn upper_kern(&self) -> f64 {
        self.ig / (self.ag * self.yb)
    }
    /// Lower kern distance.
    pub fn lower_kern(&self) -> f64 {
        self.ig / (self.ag * self.yt)
    }
    /// Bottom fiber stress at midspan under Pe + M (N/mm²).
    /// `pe` = effective prestress force (N), `m` = applied moment (N·mm).
    pub fn bottom_fiber_stress(&self, pe: f64, m: f64) -> f64 {
        let p_term = -pe / self.ag;
        let e_term = -pe * self.eccentricity * self.yb / self.ig;
        let m_term = m * self.yb / self.ig;
        p_term + e_term + m_term
    }
    /// Top fiber stress at midspan under Pe + M.
    pub fn top_fiber_stress(&self, pe: f64, m: f64) -> f64 {
        let p_term = -pe / self.ag;
        let e_term = pe * self.eccentricity * self.yt / self.ig;
        let m_term = -m * self.yt / self.ig;
        p_term + e_term + m_term
    }
    /// Cracking moment (N·mm): moment at which bottom tensile stress = fr.
    pub fn cracking_moment(&self, pe: f64) -> f64 {
        let fr = 0.62 * self.fc.sqrt();
        let sb = self.ig / self.yb;
        (pe / self.ag + pe * self.eccentricity / sb + fr) * sb
    }
    /// Elastic shortening loss (MPa) for pre-tensioned members (ACI 318 simplified).
    pub fn elastic_shortening_loss(&self) -> f64 {
        let es_steel = 197_000.0;
        let ec = 4700.0 * self.fc.sqrt();
        let n = es_steel / ec;
        let pe = self.initial_prestress_force();
        let fc_cgs = pe / self.ag + pe * self.eccentricity.powi(2) / self.ig;
        n * fc_cgs
    }
    /// Shrinkage loss (MPa) — simplified ACI estimate 70 MPa for pre-tensioned.
    pub fn shrinkage_loss(&self) -> f64 {
        70.0
    }
    /// Creep loss (MPa) — simplified: Ccr * n * fc_cgs.
    pub fn creep_loss(&self) -> f64 {
        let ccr = 2.0;
        let n = 197_000.0 / (4700.0 * self.fc.sqrt());
        let pe = self.effective_prestress_force(0.0);
        let fc_cgs = pe / self.ag + pe * self.eccentricity.powi(2) / self.ig;
        ccr * n * fc_cgs
    }
    /// Total prestress losses (MPa).
    pub fn total_losses(&self) -> f64 {
        self.elastic_shortening_loss() + self.shrinkage_loss() + self.creep_loss()
    }
    /// Nominal flexural strength Mn (N·mm) per ACI 318 (fps from strand stress).
    pub fn nominal_flexural_strength(&self) -> f64 {
        let rho_p = self.aps / (self.ag * 0.8);
        let fps = self.fpu * (1.0 - 0.5 * rho_p * self.fpu / self.fc);
        let dp = self.yb;
        let a = fps * self.aps / (0.85 * self.fc * (self.ag / self.yb));
        fps * self.aps * (dp - a / 2.0)
    }
}
/// Concrete cover requirement based on exposure class.
#[derive(Debug, Clone, PartialEq)]
pub enum ExposureClass {
    /// Interior not exposed to weather.
    Interior,
    /// Exposed to weather — moderate.
    Moderate,
    /// Exposed to deicers or aggressive chemicals.
    Severe,
    /// Submerged or buried.
    Submerged,
}
/// Masonry unit (brick or concrete block) and wall properties.
pub struct MasonryUnit {
    /// Net compressive strength of unit f'm (MPa).
    pub fm: f64,
    /// Mortar type (S, N, or M).
    pub mortar_type: String,
    /// Bond pattern (running bond or stack bond).
    pub bond_pattern: String,
    /// Effective modulus of elasticity Em (MPa) per TMS 402.
    pub em: f64,
    /// Density (kg/m³).
    pub density: f64,
}
impl MasonryUnit {
    /// Create a standard clay brick masonry unit.
    pub fn clay_brick(fm: f64) -> Self {
        let em = 700.0 * fm;
        MasonryUnit {
            fm,
            mortar_type: "S".to_string(),
            bond_pattern: "Running".to_string(),
            em,
            density: 1900.0,
        }
    }
    /// Shear modulus Gv ≈ 0.4 * Em.
    pub fn shear_modulus(&self) -> f64 {
        0.4 * self.em
    }
    /// Modulus of rupture fr ≈ 0.064 * fm (MPa) per TMS 402.
    pub fn modulus_of_rupture(&self) -> f64 {
        0.064 * self.fm
    }
}
/// Glued laminated timber (glulam) section properties per NDS / APA.
#[derive(Debug, Clone)]
pub struct GlulamSection {
    /// Width b (mm).
    pub b: f64,
    /// Total depth d (mm).
    pub d: f64,
    /// Number of laminations.
    pub n_lam: u32,
    /// Thickness of each lamination (mm).
    pub lam_thickness: f64,
    /// Reference bending design value Fb (MPa).
    pub fb: f64,
    /// Reference shear design value Fv (MPa).
    pub fv: f64,
    /// Modulus of elasticity (MPa).
    pub e: f64,
    /// Moisture service condition factor CM.
    pub cm: f64,
}
impl GlulamSection {
    /// Create a new glulam section.
    pub fn new(b: f64, d: f64, n_lam: u32, fb: f64, fv: f64, e: f64) -> Self {
        let lam_thickness = d / n_lam as f64;
        GlulamSection {
            b,
            d,
            n_lam,
            lam_thickness,
            fb,
            fv,
            e,
            cm: 1.0,
        }
    }
    /// Create a 24F-V4 Douglas Fir glulam (275 × 570 mm, 19 lams).
    pub fn df_24f_v4() -> Self {
        GlulamSection {
            b: 275.0,
            d: 570.0,
            n_lam: 19,
            lam_thickness: 30.0,
            fb: 16.5,
            fv: 2.4,
            e: 12_400.0,
            cm: 1.0,
        }
    }
    /// Cross-sectional area (mm²).
    pub fn area(&self) -> f64 {
        self.b * self.d
    }
    /// Moment of inertia Ix (mm⁴).
    pub fn ix(&self) -> f64 {
        self.b * self.d.powi(3) / 12.0
    }
    /// Section modulus Sx (mm³).
    pub fn sx(&self) -> f64 {
        self.b * self.d.powi(2) / 6.0
    }
    /// Volume factor Cv for bending members (NDS).
    /// `l` = beam span (m), reduces for long spans.
    pub fn volume_factor(&self, l_m: f64) -> f64 {
        let kl = 1.09;
        let b_ft = self.b / 25.4 / 12.0;
        let d_ft = self.d / 25.4 / 12.0;
        let l_ft = l_m * 3.28084;
        (kl * 21.0 / l_ft).powf(0.1)
            * (12.0 / (d_ft * 12.0)).powf(0.1)
            * (5.125 / (b_ft * 12.0)).powf(0.1).min(1.0)
    }
    /// Adjusted bending design value Fb' (MPa).
    pub fn adjusted_fb(&self, cd: f64, cv: f64) -> f64 {
        self.fb * cd * self.cm * cv
    }
    /// Allowable moment (N·mm).
    pub fn allowable_moment(&self, cd: f64, cv: f64) -> f64 {
        self.adjusted_fb(cd, cv) * self.sx()
    }
    /// Allowable shear (N).
    pub fn allowable_shear(&self) -> f64 {
        let fv_prime = self.fv * self.cm;
        fv_prime * (2.0 / 3.0) * self.b * self.d
    }
    /// Mid-span deflection under uniform load w (N/mm), simple span.
    pub fn midspan_deflection(&self, w_n_per_mm: f64, span_mm: f64) -> f64 {
        5.0 * w_n_per_mm * span_mm.powi(4) / (384.0 * self.e * self.ix())
    }
}
