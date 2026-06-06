// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrochemical corrosion models.
//!
//! Provides physics-based corrosion models including:
//! - Evans diagram and Butler-Volmer electrode kinetics
//! - Pitting corrosion (repassivation potential, pit growth)
//! - Galvanic corrosion (mixed potential theory)
//! - Corrosion rate in mpy and mm/year
//! - Passivation (Flade potential, passive current density)
//! - Stress corrosion cracking (KIscc, crack velocity)
//! - Crevice corrosion (critical gap geometry)
//! - Dealloying / selective dissolution
//! - Cathodic protection design (impressed-current and sacrificial anode)
//! - Corrosion inhibitor efficiency

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Faraday constant \[C/mol\].
pub const FARADAY: f64 = 96_485.332_12;

/// Universal gas constant \[J/(mol·K)\].
pub const GAS_CONSTANT: f64 = 8.314_462_618;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Clamp a value to the range \[lo, hi\].
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

/// Safe natural logarithm: returns `f64::NEG_INFINITY` for x ≤ 0.
#[inline]
fn safe_ln(x: f64) -> f64 {
    if x <= 0.0 { f64::NEG_INFINITY } else { x.ln() }
}

// ---------------------------------------------------------------------------
// Butler-Volmer electrode kinetics
// ---------------------------------------------------------------------------

/// Butler-Volmer electrode kinetics.
///
/// Models the current–overpotential relationship:
/// `i = i₀ · [exp(αa · F · η / RT) − exp(−αc · F · η / RT)]`
#[derive(Debug, Clone)]
pub struct ButlerVolmer {
    /// Exchange current density i₀ \[A/m²\].
    pub i0: f64,
    /// Anodic charge-transfer coefficient αa (dimensionless, typically 0.5).
    pub alpha_a: f64,
    /// Cathodic charge-transfer coefficient αc (dimensionless, typically 0.5).
    pub alpha_c: f64,
    /// Absolute temperature T \[K\].
    pub temperature: f64,
}

impl ButlerVolmer {
    /// Create a new Butler-Volmer model.
    pub fn new(i0: f64, alpha_a: f64, alpha_c: f64, temperature: f64) -> Self {
        Self {
            i0,
            alpha_a,
            alpha_c,
            temperature,
        }
    }

    /// Compute the dimensionless factor F/(RT) \[V⁻¹\].
    #[inline]
    pub fn f_over_rt(&self) -> f64 {
        FARADAY / (GAS_CONSTANT * self.temperature.max(1e-6))
    }

    /// Current density \[A/m²\] at overpotential `eta` \[V\].
    pub fn current_density(&self, eta: f64) -> f64 {
        let q = self.f_over_rt();
        self.i0 * ((self.alpha_a * q * eta).exp() - (-self.alpha_c * q * eta).exp())
    }

    /// Linearised charge-transfer resistance \[Ω·m²\] near equilibrium (η → 0).
    ///
    /// `Rct = RT / [i₀ · (αa + αc) · F]`
    pub fn charge_transfer_resistance(&self) -> f64 {
        GAS_CONSTANT * self.temperature
            / (self.i0.max(1e-30) * (self.alpha_a + self.alpha_c) * FARADAY)
    }

    /// Anodic Tafel slope \[V/decade\].
    ///
    /// `ba = ln(10) · RT / (αa · F)`
    pub fn tafel_slope_anodic(&self) -> f64 {
        std::f64::consts::LN_10 / (self.alpha_a * self.f_over_rt())
    }

    /// Cathodic Tafel slope \[V/decade\].
    ///
    /// `bc = ln(10) · RT / (αc · F)`
    pub fn tafel_slope_cathodic(&self) -> f64 {
        std::f64::consts::LN_10 / (self.alpha_c * self.f_over_rt())
    }

    /// Overpotential \[V\] for a target current density using the Tafel approximation.
    ///
    /// Valid for |η| >> RT/F (large overpotential limit).
    /// Sign follows the input `i_target`: positive → anodic branch.
    pub fn tafel_overpotential(&self, i_target: f64) -> f64 {
        if i_target > 0.0 {
            self.tafel_slope_anodic() * safe_ln(i_target / self.i0.max(1e-30))
                / std::f64::consts::LN_10
        } else if i_target < 0.0 {
            -self.tafel_slope_cathodic() * safe_ln((-i_target) / self.i0.max(1e-30))
                / std::f64::consts::LN_10
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// Evans diagram (mixed-potential / corrosion potential)
// ---------------------------------------------------------------------------

/// An Evans diagram combining anodic dissolution and cathodic reaction curves.
///
/// The corrosion potential is found at the intersection where
/// `i_anodic(E) + i_cathodic(E) = 0`.
#[derive(Debug, Clone)]
pub struct EvansDiagram {
    /// Anodic half-reaction Butler-Volmer model.
    pub anodic: ButlerVolmer,
    /// Equilibrium potential of the anodic (dissolution) reaction \[V vs SHE\].
    pub e_eq_anodic: f64,
    /// Cathodic half-reaction Butler-Volmer model.
    pub cathodic: ButlerVolmer,
    /// Equilibrium potential of the cathodic (reduction) reaction \[V vs SHE\].
    pub e_eq_cathodic: f64,
}

impl EvansDiagram {
    /// Create a new Evans diagram.
    pub fn new(
        anodic: ButlerVolmer,
        e_eq_anodic: f64,
        cathodic: ButlerVolmer,
        e_eq_cathodic: f64,
    ) -> Self {
        Self {
            anodic,
            e_eq_anodic,
            cathodic,
            e_eq_cathodic,
        }
    }

    /// Net current density \[A/m²\] at electrode potential `e` \[V vs SHE\].
    ///
    /// Computed as the sum of the anodic dissolution current and the cathodic
    /// reduction current (which is negative when E < E_eq_cathodic).
    /// Positive = net anodic (dissolution), negative = net cathodic.
    /// At the corrosion potential the net current equals zero.
    pub fn net_current(&self, e: f64) -> f64 {
        let eta_a = e - self.e_eq_anodic;
        let eta_c = e - self.e_eq_cathodic;
        self.anodic.current_density(eta_a) + self.cathodic.current_density(eta_c)
    }

    /// Find the mixed (corrosion) potential \[V\] by bisection over `[e_lo, e_hi]`.
    ///
    /// Returns `None` if no sign change is found in the interval.
    pub fn corrosion_potential(&self, e_lo: f64, e_hi: f64, tol: f64) -> Option<f64> {
        let fa = self.net_current(e_lo);
        let fb = self.net_current(e_hi);
        if fa * fb > 0.0 {
            return None;
        }
        let mut lo = e_lo;
        let mut hi = e_hi;
        for _ in 0..100 {
            let mid = (lo + hi) * 0.5;
            if hi - lo < tol {
                return Some(mid);
            }
            let fm = self.net_current(mid);
            if fa * fm <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Some((lo + hi) * 0.5)
    }

    /// Corrosion current density \[A/m²\] evaluated at the mixed potential.
    pub fn corrosion_current(&self, e_corr: f64) -> f64 {
        let eta_a = e_corr - self.e_eq_anodic;
        self.anodic.current_density(eta_a).abs()
    }
}

// ---------------------------------------------------------------------------
// Corrosion rate conversions
// ---------------------------------------------------------------------------

/// Convert corrosion current density to mass-loss rate \[g/(m²·s)\].
///
/// # Arguments
/// * `i_corr` – corrosion current density \[A/m²\]
/// * `molar_mass` – metal molar mass \[g/mol\]
/// * `n_electrons` – number of electrons exchanged per metal atom
pub fn corrosion_mass_rate(i_corr: f64, molar_mass: f64, n_electrons: f64) -> f64 {
    i_corr * molar_mass / (n_electrons * FARADAY)
}

/// Convert corrosion current density to penetration rate \[mm/year\].
///
/// # Arguments
/// * `i_corr` – corrosion current density \[A/m²\]
/// * `molar_mass` – metal molar mass \[g/mol\]
/// * `n_electrons` – electrons per metal atom
/// * `density` – metal density \[g/cm³\]
pub fn corrosion_rate_mm_per_year(
    i_corr: f64,
    molar_mass: f64,
    n_electrons: f64,
    density: f64,
) -> f64 {
    // mass rate [g/(m²·s)] / density [g/cm³ = 1e6 g/m³] → thickness rate [m/s]
    // then convert to mm/year
    let mass_rate = corrosion_mass_rate(i_corr, molar_mass, n_electrons);
    let thickness_rate_m_per_s = mass_rate / (density * 1e6); // density in g/m³
    thickness_rate_m_per_s * 1e3 * 3.156e7 // mm/year
}

/// Convert corrosion current density to penetration rate \[mpy\] (mils per year).
///
/// 1 mpy = 0.0254 mm/year.
pub fn corrosion_rate_mpy(i_corr: f64, molar_mass: f64, n_electrons: f64, density: f64) -> f64 {
    corrosion_rate_mm_per_year(i_corr, molar_mass, n_electrons, density) / 0.0254
}

// ---------------------------------------------------------------------------
// Passivation model
// ---------------------------------------------------------------------------

/// Passivation model for active–passive metals.
///
/// Describes the three regimes on the anodic polarisation curve:
/// active dissolution → active-passive transition → passive plateau → transpassive.
#[derive(Debug, Clone)]
pub struct PassivationModel {
    /// Flade potential (active/passive boundary) \[V vs SHE\].
    pub flade_potential: f64,
    /// Passivation current density (peak active) \[A/m²\].
    pub i_passivation: f64,
    /// Passive current density (maintained in passive regime) \[A/m²\].
    pub i_passive: f64,
    /// Transpassive potential (onset of transpassive dissolution) \[V vs SHE\].
    pub e_transpassive: f64,
    /// Butler-Volmer parameters for the active regime.
    pub active_bv: ButlerVolmer,
}

impl PassivationModel {
    /// Create a new passivation model.
    pub fn new(
        flade_potential: f64,
        i_passivation: f64,
        i_passive: f64,
        e_transpassive: f64,
        active_bv: ButlerVolmer,
    ) -> Self {
        Self {
            flade_potential,
            i_passivation,
            i_passive,
            e_transpassive,
            active_bv,
        }
    }

    /// Anodic current density \[A/m²\] at potential `e` \[V vs SHE\].
    ///
    /// Regime selection:
    /// - e < Flade potential → active Butler-Volmer dissolution
    /// - Flade ≤ e < transpassive → passive plateau (`i_passive`)
    /// - e ≥ transpassive → transpassive (linear rise above passive)
    pub fn current_density(&self, e: f64) -> f64 {
        if e < self.flade_potential {
            let eta = e - self.active_bv.f_over_rt() * 0.0; // overpotential relative to E_eq
            self.active_bv.current_density(eta).max(0.0)
        } else if e < self.e_transpassive {
            self.i_passive
        } else {
            // Transpassive: rises from i_passive
            self.i_passive + (e - self.e_transpassive) * self.i_passivation
        }
    }

    /// Check whether the metal is in the passive state at potential `e`.
    pub fn is_passive(&self, e: f64) -> bool {
        e >= self.flade_potential && e < self.e_transpassive
    }

    /// Passive range width \[V\].
    pub fn passive_range(&self) -> f64 {
        (self.e_transpassive - self.flade_potential).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// Pitting corrosion
// ---------------------------------------------------------------------------

/// Pitting corrosion model.
///
/// Characterised by the pitting potential Epit (onset of stable pitting)
/// and the repassivation potential Erp (below which pits repassivate).
#[derive(Debug, Clone)]
pub struct PittingModel {
    /// Pitting potential \[V vs SHE\]: above this potential, stable pits nucleate.
    pub e_pit: f64,
    /// Repassivation potential \[V vs SHE\]: below this, active pits repassivate.
    pub e_rp: f64,
    /// Critical chloride concentration for pitting \[mol/L\].
    pub critical_chloride: f64,
    /// Pit growth rate constant k \[m·s⁻¹·(A/m²)⁻¹\].
    pub pit_growth_k: f64,
}

impl PittingModel {
    /// Create a new pitting model.
    pub fn new(e_pit: f64, e_rp: f64, critical_chloride: f64, pit_growth_k: f64) -> Self {
        Self {
            e_pit,
            e_rp,
            critical_chloride,
            pit_growth_k,
        }
    }

    /// Hysteresis width: Epit − Erp \[V\].
    ///
    /// A larger value indicates greater susceptibility to stable pitting.
    pub fn hysteresis_width(&self) -> f64 {
        (self.e_pit - self.e_rp).max(0.0)
    }

    /// Returns `true` if stable pit growth is expected at the given potential.
    pub fn is_pitting_active(&self, e: f64) -> bool {
        e >= self.e_pit
    }

    /// Returns `true` if an active pit will repassivate at potential `e`.
    pub fn will_repassivate(&self, e: f64) -> bool {
        e < self.e_rp
    }

    /// Hemispherical pit radius \[m\] after time `t` \[s\] at pit current `i_pit` \[A/m²\].
    ///
    /// Uses Faraday's law: V = M·Q/(n·F·ρ) with hemispherical geometry.
    /// `r = (3·M·i_pit·t / (2·π·n·F·ρ))^(1/3)`
    ///
    /// # Arguments
    /// * `i_pit` – pit current density \[A/m²\]
    /// * `t` – time \[s\]
    /// * `molar_mass` – \[g/mol\]
    /// * `n` – electrons per atom
    /// * `density` – \[g/cm³\]
    pub fn pit_radius(&self, i_pit: f64, t: f64, molar_mass: f64, n: f64, density: f64) -> f64 {
        let rho_si = density * 1e6; // g/m³
        let numerator = 3.0 * molar_mass * i_pit * t;
        let denominator = 2.0 * std::f64::consts::PI * n * FARADAY * rho_si;
        if denominator <= 0.0 {
            0.0
        } else {
            (numerator / denominator).cbrt()
        }
    }

    /// Pitting initiation time \[s\] according to a simplified induction model.
    ///
    /// `t_i = A / (c_Cl / c_crit − 1)` where A is a material constant \[s\].
    pub fn induction_time(&self, chloride_conc: f64, material_constant_s: f64) -> Option<f64> {
        let ratio = chloride_conc / self.critical_chloride.max(1e-30);
        if ratio <= 1.0 {
            None // Below threshold — no pitting
        } else {
            Some(material_constant_s / (ratio - 1.0))
        }
    }
}

// ---------------------------------------------------------------------------
// Galvanic corrosion — mixed potential theory
// ---------------------------------------------------------------------------

/// Galvanic pair: two dissimilar metals electrically connected in an electrolyte.
///
/// Computes the galvanic (mixed) potential and the galvanic current density
/// using the mixed-potential (Wagner-Traud) theory.
#[derive(Debug, Clone)]
pub struct GalvanicPair {
    /// Anodic metal (more active, dissolves preferentially).
    pub anode: ButlerVolmer,
    /// Equilibrium potential of the anode \[V vs SHE\].
    pub e_eq_anode: f64,
    /// Cathodic metal (nobler, protected).
    pub cathode: ButlerVolmer,
    /// Equilibrium potential of the cathode \[V vs SHE\].
    pub e_eq_cathode: f64,
    /// Area ratio: cathode area / anode area (dimensionless).
    pub area_ratio: f64,
}

impl GalvanicPair {
    /// Create a new galvanic pair.
    pub fn new(
        anode: ButlerVolmer,
        e_eq_anode: f64,
        cathode: ButlerVolmer,
        e_eq_cathode: f64,
        area_ratio: f64,
    ) -> Self {
        Self {
            anode,
            e_eq_anode,
            cathode,
            e_eq_cathode,
            area_ratio,
        }
    }

    /// Galvanic driving force \[V\]: difference in equilibrium potentials.
    pub fn driving_force(&self) -> f64 {
        (self.e_eq_cathode - self.e_eq_anode).abs()
    }

    /// Net current density at potential `e` \[A/m² referenced to anode area\].
    ///
    /// Sum of anodic dissolution (positive) and scaled cathodic reduction (negative
    /// when E < E_eq_cathode).  Equals zero at the galvanic mixed potential.
    pub fn net_current(&self, e: f64) -> f64 {
        let i_a = self.anode.current_density(e - self.e_eq_anode);
        let i_c = self.cathode.current_density(e - self.e_eq_cathode) * self.area_ratio;
        i_a + i_c
    }

    /// Find the galvanic (mixed) potential \[V\] by bisection.
    pub fn galvanic_potential(&self, e_lo: f64, e_hi: f64, tol: f64) -> Option<f64> {
        let fa = self.net_current(e_lo);
        let fb = self.net_current(e_hi);
        if fa * fb > 0.0 {
            return None;
        }
        let mut lo = e_lo;
        let mut hi = e_hi;
        for _ in 0..100 {
            let mid = (lo + hi) * 0.5;
            if hi - lo < tol {
                return Some(mid);
            }
            let fm = self.net_current(mid);
            if fa * fm <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Some((lo + hi) * 0.5)
    }
}

// ---------------------------------------------------------------------------
// Stress corrosion cracking (SCC)
// ---------------------------------------------------------------------------

/// Stress corrosion cracking model.
///
/// Combines fracture mechanics (KI) with the electrochemical dissolution
/// mechanism (anodic dissolution model).
#[derive(Debug, Clone)]
pub struct SccrModel {
    /// Threshold stress intensity factor KIscc \[MPa·m^0.5\].
    pub k_iscc: f64,
    /// Fracture toughness KIc \[MPa·m^0.5\].
    pub k_ic: f64,
    /// Stage II (plateau) crack velocity \[m/s\].
    pub v_plateau: f64,
    /// Exponent in Stage I power-law velocity: `v = A·(KI − KIscc)^n`.
    pub stage1_exponent: f64,
    /// Coefficient A in Stage I velocity \[m/s · (MPa·m^0.5)^(−n)\].
    pub stage1_coeff: f64,
}

impl SccrModel {
    /// Create a new SCC model.
    pub fn new(
        k_iscc: f64,
        k_ic: f64,
        v_plateau: f64,
        stage1_exponent: f64,
        stage1_coeff: f64,
    ) -> Self {
        Self {
            k_iscc,
            k_ic,
            v_plateau,
            stage1_exponent,
            stage1_coeff,
        }
    }

    /// Crack propagation velocity \[m/s\] at stress intensity `ki` \[MPa·m^0.5\].
    ///
    /// Three-stage model:
    /// - KI < KIscc → 0 (no crack growth)
    /// - KIscc ≤ KI < KIc → Stage I/II (Plateau or power-law)
    /// - KI ≥ KIc → catastrophic (returns `f64::INFINITY`)
    pub fn crack_velocity(&self, ki: f64) -> f64 {
        if ki < self.k_iscc {
            0.0
        } else if ki >= self.k_ic {
            f64::INFINITY
        } else {
            let dk = ki - self.k_iscc;
            let v_stage1 = self.stage1_coeff * dk.powf(self.stage1_exponent);
            v_stage1.min(self.v_plateau)
        }
    }

    /// Time to fracture \[s\] from an initial crack of length `a0` \[m\]
    /// to the critical length `ac` \[m\] under constant stress `sigma` \[MPa\].
    ///
    /// Critical crack length: `ac = (KIc / (Y · σ · √π))²`
    /// Uses simple trapezoidal integration over crack length.
    ///
    /// # Arguments
    /// * `a0` – initial half-crack length \[m\]
    /// * `sigma` – applied stress \[MPa\]
    /// * `geometry_factor_y` – geometry correction factor Y (≈ 1.0 for centre crack)
    /// * `steps` – integration steps
    pub fn time_to_fracture(
        &self,
        a0: f64,
        sigma: f64,
        geometry_factor_y: f64,
        steps: usize,
    ) -> f64 {
        let ac = (self.k_ic / (geometry_factor_y * sigma * std::f64::consts::PI.sqrt())).powi(2);
        if a0 >= ac {
            return 0.0;
        }
        let da = (ac - a0) / steps.max(1) as f64;
        let ki = |a: f64| geometry_factor_y * sigma * (std::f64::consts::PI * a).sqrt();
        let dt = |a: f64| {
            let v = self.crack_velocity(ki(a));
            if v <= 0.0 || v.is_infinite() {
                0.0
            } else {
                da / v
            }
        };
        let mut t = 0.0f64;
        let mut a = a0;
        for _ in 0..steps {
            t += dt(a);
            a += da;
        }
        t
    }

    /// Susceptibility index (0–1): `(KIscc / KIc)`.
    ///
    /// Values close to 1 indicate low susceptibility (wide safe window).
    pub fn susceptibility_index(&self) -> f64 {
        (self.k_iscc / self.k_ic.max(1e-30)).min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Crevice corrosion
// ---------------------------------------------------------------------------

/// Crevice corrosion model based on critical geometry criterion.
///
/// Crevice attack initiates when the IR drop within the crevice
/// is large enough to shift the local potential below the critical crevice
/// potential.
#[derive(Debug, Clone)]
pub struct CreviceModel {
    /// Crevice gap width \[m\].
    pub gap: f64,
    /// Crevice depth \[m\].
    pub depth: f64,
    /// Electrolyte resistivity \[Ω·m\].
    pub resistivity: f64,
    /// Critical crevice potential (below which active dissolution occurs) \[V\].
    pub e_crevice_critical: f64,
    /// Passive current density \[A/m²\].
    pub i_passive: f64,
}

impl CreviceModel {
    /// Create a new crevice model.
    pub fn new(
        gap: f64,
        depth: f64,
        resistivity: f64,
        e_crevice_critical: f64,
        i_passive: f64,
    ) -> Self {
        Self {
            gap,
            depth,
            resistivity,
            e_crevice_critical,
            i_passive,
        }
    }

    /// IR drop along the crevice \[V\] using 1-D ohmic model.
    ///
    /// `ΔV = i_passive · ρ · L² / (2 · g)`
    /// where L = depth, g = gap, ρ = resistivity.
    pub fn ir_drop(&self) -> f64 {
        self.i_passive * self.resistivity * self.depth.powi(2) / (2.0 * self.gap.max(1e-15))
    }

    /// Critical depth \[m\] at which crevice attack initiates for a given external potential.
    ///
    /// `L_crit = sqrt(2 · g · ΔV_crit / (i_passive · ρ))`
    pub fn critical_depth(&self, e_external: f64) -> Option<f64> {
        let dv_crit = e_external - self.e_crevice_critical;
        if dv_crit <= 0.0 {
            return None;
        }
        let arg =
            2.0 * self.gap * dv_crit / (self.i_passive.max(1e-30) * self.resistivity.max(1e-30));
        Some(arg.sqrt())
    }

    /// Geometry factor x = i_passive · ρ · depth² / (2 · gap · ΔV).
    ///
    /// Crevice attack expected when x > 1.
    pub fn geometry_factor(&self, e_external: f64) -> f64 {
        let dv = (e_external - self.e_crevice_critical).max(1e-30);
        self.i_passive * self.resistivity * self.depth.powi(2) / (2.0 * self.gap.max(1e-15) * dv)
    }
}

// ---------------------------------------------------------------------------
// Dealloying / selective dissolution
// ---------------------------------------------------------------------------

/// Selective dissolution (dealloying) model.
///
/// Describes the preferential dissolution of the more active component
/// from a binary alloy (e.g., Cu from brass, Zn from α-brass).
#[derive(Debug, Clone)]
pub struct DealloyingModel {
    /// Initial mole fraction of the active component (0–1).
    pub x0: f64,
    /// Dealloying threshold potential \[V vs SHE\].
    pub e_threshold: f64,
    /// Diffusion coefficient of the active component in the dealloyed layer \[m²/s\].
    pub d_alloy: f64,
    /// Dissolution rate constant \[mol/(m²·s)\] above the threshold.
    pub k_dissolution: f64,
}

impl DealloyingModel {
    /// Create a new dealloying model.
    pub fn new(x0: f64, e_threshold: f64, d_alloy: f64, k_dissolution: f64) -> Self {
        Self {
            x0,
            e_threshold,
            d_alloy,
            k_dissolution,
        }
    }

    /// Returns `true` if dealloying is thermodynamically active at potential `e`.
    pub fn is_active(&self, e: f64) -> bool {
        e >= self.e_threshold
    }

    /// Dealloyed layer thickness \[m\] after time `t` \[s\] (parabolic growth law).
    ///
    /// `δ = sqrt(2 · D · t)` — valid for diffusion-controlled regime.
    pub fn layer_thickness(&self, t: f64) -> f64 {
        (2.0 * self.d_alloy * t).sqrt()
    }

    /// Residual active-component fraction in the dealloyed layer after depth `delta` \[m\].
    ///
    /// Uses an exponential depletion profile: `x(δ) = x0 · exp(−k · δ / D)`.
    pub fn residual_fraction(&self, delta: f64) -> f64 {
        let decay = (self.k_dissolution * delta / self.d_alloy.max(1e-30)).min(700.0);
        self.x0 * (-decay).exp()
    }

    /// Total active-component mass dissolved per unit area \[mol/m²\] up to time `t`.
    pub fn dissolved_moles_per_area(&self, t: f64) -> f64 {
        self.k_dissolution * t
    }
}

// ---------------------------------------------------------------------------
// Cathodic protection
// ---------------------------------------------------------------------------

/// Cathodic protection design by impressed current.
///
/// Determines the required protective current and anode specifications
/// to achieve a target protective potential on a steel structure.
#[derive(Debug, Clone)]
pub struct ImpressedCurrentCp {
    /// Structure surface area to be protected \[m²\].
    pub surface_area: f64,
    /// Protective current density (typically 0.01–0.1 A/m² for bare steel) \[A/m²\].
    pub protective_current_density: f64,
    /// Current utilisation efficiency of the anode (0–1, accounts for stray current).
    pub efficiency: f64,
    /// Design lifetime \[years\].
    pub design_life_years: f64,
}

impl ImpressedCurrentCp {
    /// Create a new impressed-current CP design.
    pub fn new(
        surface_area: f64,
        protective_current_density: f64,
        efficiency: f64,
        design_life_years: f64,
    ) -> Self {
        Self {
            surface_area,
            protective_current_density,
            efficiency,
            design_life_years,
        }
    }

    /// Required total current output \[A\].
    pub fn required_current(&self) -> f64 {
        self.surface_area * self.protective_current_density
    }

    /// Required rectifier output current \[A\] accounting for efficiency.
    pub fn rectifier_current(&self) -> f64 {
        self.required_current() / self.efficiency.max(1e-6)
    }

    /// Total charge \[A·h\] needed over the design lifetime.
    pub fn total_charge_ah(&self) -> f64 {
        self.rectifier_current() * self.design_life_years * 8760.0
    }
}

/// Sacrificial (galvanic) anode cathodic protection design.
///
/// Computes the required anode mass for a given protective current
/// and design life using the electrochemical capacity of the anode material.
#[derive(Debug, Clone)]
pub struct SacrificialAnode {
    /// Electrochemical capacity of the anode material \[A·h/kg\]
    /// (e.g., zinc ≈ 780, aluminium ≈ 2600, magnesium ≈ 1100).
    pub electrochemical_capacity: f64,
    /// Current output per anode \[A\].
    pub current_output: f64,
    /// Design lifetime \[years\].
    pub design_life_years: f64,
    /// Current efficiency (fraction of capacity delivering useful protection, 0–1).
    pub current_efficiency: f64,
}

impl SacrificialAnode {
    /// Create a new sacrificial anode specification.
    pub fn new(
        electrochemical_capacity: f64,
        current_output: f64,
        design_life_years: f64,
        current_efficiency: f64,
    ) -> Self {
        Self {
            electrochemical_capacity,
            current_output,
            design_life_years,
            current_efficiency,
        }
    }

    /// Mass of one anode required \[kg\].
    pub fn anode_mass(&self) -> f64 {
        let ah_required = self.current_output * self.design_life_years * 8760.0;
        ah_required / (self.electrochemical_capacity * self.current_efficiency.max(1e-6))
    }

    /// Number of anodes required to supply `total_current` \[A\].
    pub fn anode_count(&self, total_current: f64) -> u64 {
        let n = (total_current / self.current_output.max(1e-30)).ceil();
        n as u64
    }
}

// ---------------------------------------------------------------------------
// Corrosion inhibitor efficiency
// ---------------------------------------------------------------------------

/// Corrosion inhibitor performance model.
///
/// Inhibitors adsorb onto the metal surface and reduce the exposed area
/// available for electrochemical dissolution.  Langmuir adsorption
/// isotherm is used to model surface coverage.
#[derive(Debug, Clone)]
pub struct CorrosionInhibitor {
    /// Langmuir adsorption equilibrium constant K \[L/mol\].
    pub langmuir_k: f64,
    /// Inhibitor concentration in solution \[mol/L\].
    pub concentration: f64,
    /// Uninhibited corrosion rate (reference) \[mm/year\].
    pub uninhibited_rate: f64,
}

impl CorrosionInhibitor {
    /// Create a new inhibitor model.
    pub fn new(langmuir_k: f64, concentration: f64, uninhibited_rate: f64) -> Self {
        Self {
            langmuir_k,
            concentration,
            uninhibited_rate,
        }
    }

    /// Surface coverage θ via Langmuir isotherm.
    ///
    /// `θ = K·C / (1 + K·C)`
    pub fn surface_coverage(&self) -> f64 {
        let kc = self.langmuir_k * self.concentration;
        kc / (1.0 + kc)
    }

    /// Inhibited corrosion rate \[mm/year\]: `r = r₀ · (1 − θ)`.
    pub fn inhibited_rate(&self) -> f64 {
        self.uninhibited_rate * (1.0 - self.surface_coverage())
    }

    /// Inhibitor efficiency IE \[%\]: `IE = (1 − r/r₀) × 100`.
    pub fn efficiency_percent(&self) -> f64 {
        100.0 * self.surface_coverage()
    }

    /// Minimum inhibitor concentration \[mol/L\] for target efficiency `ie_target` (0–1).
    ///
    /// Inverts the Langmuir isotherm: `C = θ / (K · (1 − θ))`.
    pub fn concentration_for_efficiency(&self, ie_target: f64) -> Option<f64> {
        let theta = clamp(ie_target, 0.0, 0.9999);
        if self.langmuir_k <= 0.0 {
            return None;
        }
        Some(theta / (self.langmuir_k * (1.0 - theta)))
    }
}

// ---------------------------------------------------------------------------
// Nernst equation and equilibrium potential
// ---------------------------------------------------------------------------

/// Nernst equation for equilibrium potential.
///
/// `E = E° + (RT / nF) · ln(a_ox / a_red)`
#[derive(Debug, Clone)]
pub struct NernstEquation {
    /// Standard electrode potential E° \[V vs SHE\].
    pub e_standard: f64,
    /// Number of electrons transferred.
    pub n_electrons: f64,
    /// Absolute temperature \[K\].
    pub temperature: f64,
}

impl NernstEquation {
    /// Create a new Nernst model.
    pub fn new(e_standard: f64, n_electrons: f64, temperature: f64) -> Self {
        Self {
            e_standard,
            n_electrons,
            temperature,
        }
    }

    /// Compute the equilibrium potential \[V\] for given oxidised and reduced activities.
    pub fn equilibrium_potential(&self, a_ox: f64, a_red: f64) -> f64 {
        let rt_nf = GAS_CONSTANT * self.temperature / (self.n_electrons * FARADAY);
        self.e_standard + rt_nf * safe_ln(a_ox / a_red.max(1e-30))
    }

    /// Nernst slope \[V per decade of concentration ratio\].
    pub fn nernst_slope_per_decade(&self) -> f64 {
        GAS_CONSTANT * self.temperature * std::f64::consts::LN_10 / (self.n_electrons * FARADAY)
    }
}

// ---------------------------------------------------------------------------
// Pourbaix diagram helper
// ---------------------------------------------------------------------------

/// A single line (boundary) in a simplified Pourbaix (E–pH) diagram.
///
/// Represents either:
/// - A pH-independent line: `E = const`
/// - A pH-dependent line: `E = a + b·pH`
#[derive(Debug, Clone)]
pub struct PourbaixLine {
    /// Line label (e.g., "Fe/Fe²⁺ dissolution boundary").
    pub label: String,
    /// Constant term a \[V\].
    pub a: f64,
    /// pH slope b \[V/pH unit\].
    pub b: f64,
}

impl PourbaixLine {
    /// Create a new Pourbaix diagram line.
    pub fn new(label: &str, a: f64, b: f64) -> Self {
        Self {
            label: label.to_string(),
            a,
            b,
        }
    }

    /// Potential \[V\] at the given pH.
    pub fn potential(&self, ph: f64) -> f64 {
        self.a + self.b * ph
    }

    /// pH at which this line intersects potential `e`.
    ///
    /// Returns `None` if the line is pH-independent (b ≈ 0).
    pub fn ph_at_potential(&self, e: f64) -> Option<f64> {
        if self.b.abs() < 1e-15 {
            None
        } else {
            Some((e - self.a) / self.b)
        }
    }
}

/// A simplified Pourbaix (E–pH) diagram for a metal.
#[derive(Debug, Clone)]
pub struct PourbaixDiagram {
    /// Collection of boundary lines.
    pub lines: Vec<PourbaixLine>,
    /// Metal name.
    pub metal: String,
}

impl PourbaixDiagram {
    /// Create a new Pourbaix diagram.
    pub fn new(metal: &str) -> Self {
        Self {
            metal: metal.to_string(),
            lines: Vec::new(),
        }
    }

    /// Add a boundary line.
    pub fn add_line(&mut self, line: PourbaixLine) {
        self.lines.push(line);
    }

    /// Build a simplified iron Pourbaix diagram with typical lines.
    pub fn iron() -> Self {
        let mut d = Self::new("Fe");
        // Fe/Fe²⁺: E = -0.44 - 0.0592*pH  (simplified)
        d.add_line(PourbaixLine::new("Fe/Fe2+", -0.44, -0.0592));
        // Fe²⁺/Fe₃O₄: E ≈ -0.059 - 0.0888*pH
        d.add_line(PourbaixLine::new("Fe2+/Fe3O4", -0.059, -0.0888));
        // Fe₃O₄/Fe₂O₃: E ≈ 0.221 - 0.0592*pH
        d.add_line(PourbaixLine::new("Fe3O4/Fe2O3", 0.221, -0.0592));
        // H₂/H₂O: E = 0 - 0.0592*pH
        d.add_line(PourbaixLine::new("H2/H2O", 0.0, -0.0592));
        // O₂/H₂O: E = 1.229 - 0.0592*pH
        d.add_line(PourbaixLine::new("O2/H2O", 1.229, -0.0592));
        d
    }

    /// Evaluate all line potentials at a given pH and return sorted values.
    pub fn potentials_at_ph(&self, ph: f64) -> Vec<(String, f64)> {
        let mut v: Vec<(String, f64)> = self
            .lines
            .iter()
            .map(|l| (l.label.clone(), l.potential(ph)))
            .collect();
        v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        v
    }
}

// ---------------------------------------------------------------------------
// Uniform corrosion under mass transport limitation
// ---------------------------------------------------------------------------

/// Limiting diffusion current density model.
///
/// When the cathodic reaction (e.g., oxygen reduction) is transport-limited,
/// the corrosion rate is controlled by the flux of oxidant to the surface.
/// `i_lim = n · F · D · c_bulk / δ`
#[derive(Debug, Clone)]
pub struct DiffusionLimitedCorrosion {
    /// Number of electrons in the cathodic reaction.
    pub n_electrons: f64,
    /// Diffusion coefficient of the oxidant \[m²/s\].
    pub diffusivity: f64,
    /// Bulk oxidant concentration \[mol/m³\].
    pub c_bulk: f64,
    /// Diffusion layer thickness \[m\].
    pub delta: f64,
}

impl DiffusionLimitedCorrosion {
    /// Create a new diffusion-limited corrosion model.
    pub fn new(n_electrons: f64, diffusivity: f64, c_bulk: f64, delta: f64) -> Self {
        Self {
            n_electrons,
            diffusivity,
            c_bulk,
            delta,
        }
    }

    /// Limiting current density \[A/m²\].
    pub fn limiting_current(&self) -> f64 {
        self.n_electrons * FARADAY * self.diffusivity * self.c_bulk / self.delta.max(1e-15)
    }

    /// Effective corrosion rate \[mm/year\] if the corrosion current equals the limiting current.
    pub fn corrosion_rate_mm_per_year(&self, molar_mass: f64, n_anodic: f64, density: f64) -> f64 {
        corrosion_rate_mm_per_year(self.limiting_current(), molar_mass, n_anodic, density)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn standard_bv() -> ButlerVolmer {
        ButlerVolmer::new(1e-3, 0.5, 0.5, 298.15)
    }

    // ── Butler-Volmer ─────────────────────────────────────────────────────────

    #[test]
    fn test_bv_zero_overpotential() {
        let bv = standard_bv();
        let i = bv.current_density(0.0);
        assert!(i.abs() < 1e-10, "i at η=0 should be 0, got {i}");
    }

    #[test]
    fn test_bv_positive_overpotential_anodic() {
        let bv = standard_bv();
        assert!(bv.current_density(0.1) > 0.0);
    }

    #[test]
    fn test_bv_negative_overpotential_cathodic() {
        let bv = standard_bv();
        assert!(bv.current_density(-0.1) < 0.0);
    }

    #[test]
    fn test_bv_charge_transfer_resistance_positive() {
        let bv = standard_bv();
        assert!(bv.charge_transfer_resistance() > 0.0);
    }

    #[test]
    fn test_bv_tafel_slopes_positive() {
        let bv = standard_bv();
        assert!(bv.tafel_slope_anodic() > 0.0);
        assert!(bv.tafel_slope_cathodic() > 0.0);
    }

    #[test]
    fn test_bv_tafel_slope_symmetric() {
        let bv = standard_bv(); // αa = αc = 0.5
        let diff = (bv.tafel_slope_anodic() - bv.tafel_slope_cathodic()).abs();
        assert!(
            diff < 1e-10,
            "symmetric αa=αc should give equal Tafel slopes, diff={diff}"
        );
    }

    #[test]
    fn test_bv_tafel_overpotential_sign() {
        let bv = standard_bv();
        assert!(bv.tafel_overpotential(1e-2) > 0.0);
        assert!(bv.tafel_overpotential(-1e-2) < 0.0);
        assert_eq!(bv.tafel_overpotential(0.0), 0.0);
    }

    #[test]
    fn test_bv_f_over_rt_positive() {
        let bv = standard_bv();
        assert!(bv.f_over_rt() > 0.0);
    }

    // ── Evans diagram ─────────────────────────────────────────────────────────

    #[test]
    fn test_evans_corrosion_potential_found() {
        // Bracket spans both equilibrium potentials (-0.5 and 0.0)
        // → net current changes sign and corrosion potential exists between them
        let anodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let cathodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let diagram = EvansDiagram::new(anodic, -0.5, cathodic, 0.0);
        let e_corr = diagram.corrosion_potential(-0.5, 0.0, 1e-6);
        assert!(e_corr.is_some());
    }

    #[test]
    fn test_evans_net_current_sign_change() {
        let anodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let cathodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let diagram = EvansDiagram::new(anodic, -0.5, cathodic, 0.0);
        // At e=-0.5 (anode eq): anodic = 0, cathodic(eta=-0.5) < 0 → net < 0
        // At e=0.0 (cathode eq): cathodic = 0, anodic(eta=0.5) > 0 → net > 0
        let ia = diagram.net_current(-0.5);
        let ib = diagram.net_current(0.0);
        assert!(ia < 0.0, "net at e_eq_anode should be < 0, got {ia}");
        assert!(ib > 0.0, "net at e_eq_cathode should be > 0, got {ib}");
    }

    #[test]
    fn test_evans_corrosion_current_positive() {
        let anodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let cathodic = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let diagram = EvansDiagram::new(anodic, -0.5, cathodic, 0.0);
        if let Some(e_corr) = diagram.corrosion_potential(-0.5, 0.0, 1e-6) {
            assert!(diagram.corrosion_current(e_corr) >= 0.0);
        }
    }

    // ── Corrosion rate conversions ────────────────────────────────────────────

    #[test]
    fn test_corrosion_mass_rate_positive() {
        // Iron: M=55.845 g/mol, n=2
        let rate = corrosion_mass_rate(1.0, 55.845, 2.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_corrosion_rate_mm_per_year_positive() {
        let rate = corrosion_rate_mm_per_year(1.0, 55.845, 2.0, 7.87);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_corrosion_rate_mpy_greater_than_mm() {
        let mm = corrosion_rate_mm_per_year(1.0, 55.845, 2.0, 7.87);
        let mpy = corrosion_rate_mpy(1.0, 55.845, 2.0, 7.87);
        assert!(mpy > mm, "mpy={mpy} should be > mm/year={mm}");
    }

    // ── Passivation ───────────────────────────────────────────────────────────

    #[test]
    fn test_passivation_is_passive() {
        let bv = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let model = PassivationModel::new(-0.2, 1e-2, 1e-5, 1.5, bv);
        assert!(model.is_passive(0.5));
        assert!(!model.is_passive(-0.5));
        assert!(!model.is_passive(2.0));
    }

    #[test]
    fn test_passivation_passive_range() {
        let bv = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let model = PassivationModel::new(-0.2, 1e-2, 1e-5, 1.5, bv);
        assert!((model.passive_range() - 1.7).abs() < 1e-10);
    }

    #[test]
    fn test_passivation_current_passive_plateau() {
        let bv = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let model = PassivationModel::new(-0.2, 1e-2, 1e-5, 1.5, bv);
        let i = model.current_density(0.5);
        assert!(
            (i - 1e-5).abs() < 1e-12,
            "expected passive i_passive, got {i}"
        );
    }

    // ── Pitting ───────────────────────────────────────────────────────────────

    #[test]
    fn test_pitting_active() {
        let model = PittingModel::new(0.3, 0.1, 0.01, 1e-8);
        assert!(model.is_pitting_active(0.5));
        assert!(!model.is_pitting_active(0.2));
    }

    #[test]
    fn test_pitting_repassivation() {
        let model = PittingModel::new(0.3, 0.1, 0.01, 1e-8);
        assert!(model.will_repassivate(0.05));
        assert!(!model.will_repassivate(0.15));
    }

    #[test]
    fn test_pitting_pit_radius_grows_with_time() {
        let model = PittingModel::new(0.3, 0.1, 0.01, 1e-8);
        let r1 = model.pit_radius(100.0, 3600.0, 55.845, 2.0, 7.87);
        let r2 = model.pit_radius(100.0, 7200.0, 55.845, 2.0, 7.87);
        assert!(r2 > r1);
    }

    #[test]
    fn test_pitting_induction_time_none_below_threshold() {
        let model = PittingModel::new(0.3, 0.1, 0.01, 1e-8);
        let t = model.induction_time(0.005, 1000.0); // below critical_chloride
        assert!(t.is_none());
    }

    #[test]
    fn test_pitting_induction_time_some_above_threshold() {
        let model = PittingModel::new(0.3, 0.1, 0.01, 1e-8);
        let t = model.induction_time(0.02, 1000.0);
        assert!(t.is_some());
        assert!(t.unwrap() > 0.0);
    }

    #[test]
    fn test_pitting_hysteresis_width() {
        let model = PittingModel::new(0.5, 0.2, 0.01, 1e-8);
        assert!((model.hysteresis_width() - 0.3).abs() < 1e-10);
    }

    // ── Galvanic corrosion ────────────────────────────────────────────────────

    #[test]
    fn test_galvanic_driving_force() {
        let a = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let c = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let pair = GalvanicPair::new(a, -0.44, c, 0.34, 10.0);
        assert!((pair.driving_force() - 0.78).abs() < 1e-10);
    }

    #[test]
    fn test_galvanic_potential_found() {
        // Bracket spans both equilibrium potentials → guaranteed sign change
        let a = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let c = ButlerVolmer::new(1e-4, 0.5, 0.5, 298.15);
        let pair = GalvanicPair::new(a, -0.44, c, 0.34, 1.0);
        // At e=-0.44: anode at eq, cathodic term dominates → net < 0
        // At e=0.34: cathode at eq, anodic term dominates → net > 0
        let e_g = pair.galvanic_potential(-0.44, 0.34, 1e-6);
        assert!(e_g.is_some(), "galvanic potential should be found");
    }

    // ── SCC ───────────────────────────────────────────────────────────────────

    #[test]
    fn test_scc_no_growth_below_kiscc() {
        let model = SccrModel::new(10.0, 50.0, 1e-7, 2.0, 1e-9);
        assert_eq!(model.crack_velocity(5.0), 0.0);
    }

    #[test]
    fn test_scc_infinite_at_kic() {
        let model = SccrModel::new(10.0, 50.0, 1e-7, 2.0, 1e-9);
        assert!(model.crack_velocity(50.0).is_infinite());
    }

    #[test]
    fn test_scc_velocity_finite_midrange() {
        let model = SccrModel::new(10.0, 50.0, 1e-7, 2.0, 1e-9);
        let v = model.crack_velocity(25.0);
        assert!(v.is_finite() && v >= 0.0);
    }

    #[test]
    fn test_scc_susceptibility_index_bounded() {
        let model = SccrModel::new(10.0, 50.0, 1e-7, 2.0, 1e-9);
        let s = model.susceptibility_index();
        assert!((0.0..=1.0).contains(&s));
    }

    #[test]
    fn test_scc_time_to_fracture_positive() {
        let model = SccrModel::new(10.0, 80.0, 1e-7, 2.0, 1e-9);
        let t = model.time_to_fracture(1e-4, 100.0, 1.0, 100);
        assert!(t > 0.0 && t.is_finite());
    }

    // ── Crevice ───────────────────────────────────────────────────────────────

    #[test]
    fn test_crevice_ir_drop_positive() {
        let model = CreviceModel::new(1e-4, 1e-2, 0.1, -0.3, 1e-5);
        assert!(model.ir_drop() > 0.0);
    }

    #[test]
    fn test_crevice_critical_depth_some() {
        let model = CreviceModel::new(1e-4, 1e-2, 0.1, -0.3, 1e-5);
        let d = model.critical_depth(-0.1);
        assert!(d.is_some());
        assert!(d.unwrap() > 0.0);
    }

    #[test]
    fn test_crevice_critical_depth_none_below_critical() {
        let model = CreviceModel::new(1e-4, 1e-2, 0.1, 0.5, 1e-5);
        // e_external = 0.3 < e_crevice_critical = 0.5 → None
        let d = model.critical_depth(0.3);
        assert!(d.is_none());
    }

    // ── Dealloying ────────────────────────────────────────────────────────────

    #[test]
    fn test_dealloying_active() {
        let model = DealloyingModel::new(0.3, 0.1, 1e-15, 1e-8);
        assert!(model.is_active(0.2));
        assert!(!model.is_active(0.05));
    }

    #[test]
    fn test_dealloying_layer_grows_with_time() {
        let model = DealloyingModel::new(0.3, 0.1, 1e-15, 1e-8);
        let d1 = model.layer_thickness(1000.0);
        let d2 = model.layer_thickness(4000.0);
        assert!(d2 > d1);
    }

    #[test]
    fn test_dealloying_residual_fraction_decreases() {
        let model = DealloyingModel::new(0.3, 0.1, 1e-15, 1e-8);
        let x1 = model.residual_fraction(1e-6);
        let x2 = model.residual_fraction(1e-5);
        assert!(x2 <= x1);
    }

    // ── Cathodic protection ───────────────────────────────────────────────────

    #[test]
    fn test_impressed_current_required() {
        let cp = ImpressedCurrentCp::new(500.0, 0.02, 0.9, 20.0);
        assert!((cp.required_current() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_impressed_current_rectifier_greater_than_required() {
        let cp = ImpressedCurrentCp::new(500.0, 0.02, 0.9, 20.0);
        assert!(cp.rectifier_current() >= cp.required_current());
    }

    #[test]
    fn test_sacrificial_anode_mass_positive() {
        let anode = SacrificialAnode::new(780.0, 0.1, 20.0, 0.85);
        assert!(anode.anode_mass() > 0.0);
    }

    #[test]
    fn test_sacrificial_anode_count() {
        let anode = SacrificialAnode::new(780.0, 0.1, 20.0, 0.85);
        let n = anode.anode_count(1.0);
        assert!(n > 0);
    }

    // ── Inhibitor ─────────────────────────────────────────────────────────────

    #[test]
    fn test_inhibitor_coverage_between_0_and_1() {
        let inh = CorrosionInhibitor::new(100.0, 0.01, 5.0);
        let theta = inh.surface_coverage();
        assert!(theta > 0.0 && theta < 1.0);
    }

    #[test]
    fn test_inhibitor_rate_less_than_uninhibited() {
        let inh = CorrosionInhibitor::new(100.0, 0.01, 5.0);
        assert!(inh.inhibited_rate() < inh.uninhibited_rate);
    }

    #[test]
    fn test_inhibitor_efficiency_positive() {
        let inh = CorrosionInhibitor::new(50.0, 0.05, 3.0);
        assert!(inh.efficiency_percent() > 0.0 && inh.efficiency_percent() < 100.0);
    }

    #[test]
    fn test_inhibitor_concentration_for_efficiency() {
        let inh = CorrosionInhibitor::new(100.0, 0.0, 5.0); // concentration = 0 → θ = 0
        let c = inh.concentration_for_efficiency(0.9);
        assert!(c.is_some());
        assert!(c.unwrap() > 0.0);
    }

    // ── Nernst equation ───────────────────────────────────────────────────────

    #[test]
    fn test_nernst_standard_potential_at_unit_activity() {
        let nernst = NernstEquation::new(-0.44, 2.0, 298.15);
        let e = nernst.equilibrium_potential(1.0, 1.0);
        assert!((e - (-0.44)).abs() < 1e-10);
    }

    #[test]
    fn test_nernst_slope_positive() {
        let nernst = NernstEquation::new(0.0, 1.0, 298.15);
        assert!(nernst.nernst_slope_per_decade() > 0.0);
    }

    #[test]
    fn test_nernst_59mv_rule() {
        // At 298.15 K, n=1: slope ≈ 0.05916 V/decade
        let nernst = NernstEquation::new(0.0, 1.0, 298.15);
        let slope = nernst.nernst_slope_per_decade();
        assert!(
            (slope - 0.05916).abs() < 2e-4,
            "Expected ~59.16 mV/decade, got {}",
            slope * 1000.0
        );
    }

    // ── Pourbaix diagram ──────────────────────────────────────────────────────

    #[test]
    fn test_pourbaix_line_potential() {
        let line = PourbaixLine::new("test", 0.0, -0.0592);
        let e = line.potential(7.0);
        assert!((e - (-0.0592 * 7.0)).abs() < 1e-10);
    }

    #[test]
    fn test_pourbaix_line_ph_at_potential() {
        let line = PourbaixLine::new("test", 0.0, -0.0592);
        // E = -0.0592 * pH  →  pH = E / -0.0592
        let e_target = -0.0592 * 7.0; // exact value for pH=7
        let ph = line.ph_at_potential(e_target);
        assert!(ph.is_some());
        assert!(
            (ph.unwrap() - 7.0).abs() < 1e-9,
            "expected pH≈7 but got {}",
            ph.unwrap()
        );
    }

    #[test]
    fn test_pourbaix_iron_lines_count() {
        let d = PourbaixDiagram::iron();
        assert!(d.lines.len() >= 5);
    }

    #[test]
    fn test_pourbaix_potentials_at_ph_sorted() {
        let d = PourbaixDiagram::iron();
        let potentials = d.potentials_at_ph(7.0);
        for i in 1..potentials.len() {
            assert!(potentials[i].1 >= potentials[i - 1].1);
        }
    }

    // ── Diffusion-limited corrosion ───────────────────────────────────────────

    #[test]
    fn test_diffusion_limiting_current_positive() {
        // Oxygen reduction in seawater: D≈2e-9 m²/s, c≈8 mol/m³ (≈8 mg/L), δ≈0.1 mm
        let model = DiffusionLimitedCorrosion::new(4.0, 2e-9, 8.0, 1e-4);
        assert!(model.limiting_current() > 0.0);
    }

    #[test]
    fn test_diffusion_corrosion_rate_positive() {
        let model = DiffusionLimitedCorrosion::new(4.0, 2e-9, 8.0, 1e-4);
        let rate = model.corrosion_rate_mm_per_year(55.845, 2.0, 7.87);
        assert!(rate > 0.0);
    }
}
