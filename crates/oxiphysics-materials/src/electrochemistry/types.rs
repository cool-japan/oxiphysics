//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{FARADAY, GAS_CONSTANT};

/// Evans diagram for mixed potential / corrosion potential analysis.
///
/// In the Evans (mixed potential) model, corrosion occurs at the potential where
/// the anodic current of the metal oxidation equals the cathodic current of the
/// oxidant reduction. This potential is the corrosion (or mixed) potential E_corr.
///
/// The model uses two Butler-Volmer half-reactions:
/// - Anodic (metal dissolution): i_a = i0_a * exp(alpha_a * F * (E - E0_a) / RT)
/// - Cathodic (reduction): i_c = i0_c * exp(-alpha_c * F * (E - E0_c) / RT)
#[derive(Debug, Clone)]
pub struct EvansDiagram {
    /// Exchange current density for anodic reaction \[A/m²\]
    pub i0_anodic: f64,
    /// Standard potential for anodic reaction \[V\]
    pub e0_anodic: f64,
    /// Anodic transfer coefficient αa
    pub alpha_a: f64,
    /// Exchange current density for cathodic reaction \[A/m²\]
    pub i0_cathodic: f64,
    /// Standard potential for cathodic reaction \[V\]
    pub e0_cathodic: f64,
    /// Cathodic transfer coefficient αc
    pub alpha_c: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}
impl EvansDiagram {
    /// Create a new Evans diagram model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        i0_anodic: f64,
        e0_anodic: f64,
        alpha_a: f64,
        i0_cathodic: f64,
        e0_cathodic: f64,
        alpha_c: f64,
        temperature: f64,
    ) -> Self {
        Self {
            i0_anodic,
            e0_anodic,
            alpha_a,
            i0_cathodic,
            e0_cathodic,
            alpha_c,
            temperature,
        }
    }
    /// F/(RT) at current temperature \[1/V\].
    pub fn f_over_rt(&self) -> f64 {
        FARADAY / (GAS_CONSTANT * self.temperature)
    }
    /// Anodic (oxidation) current density \[A/m²\] at potential E \[V\].
    ///
    /// Tafel form: `i_a = i0_a * exp(α_a · F · (E - E0_a) / RT)`
    pub fn anodic_current(&self, e: f64) -> f64 {
        let q = self.f_over_rt();
        self.i0_anodic * (self.alpha_a * q * (e - self.e0_anodic)).exp()
    }
    /// Cathodic (reduction) current density \[A/m²\] at potential E \[V\].
    ///
    /// Tafel form: `i_c = i0_c * exp(-α_c · F · (E - E0_c) / RT)`
    pub fn cathodic_current(&self, e: f64) -> f64 {
        let q = self.f_over_rt();
        self.i0_cathodic * (-self.alpha_c * q * (e - self.e0_cathodic)).exp()
    }
    /// Net current density \[A/m²\] at potential E: `i_net = i_a - i_c`.
    pub fn net_current(&self, e: f64) -> f64 {
        self.anodic_current(e) - self.cathodic_current(e)
    }
    /// Mixed (corrosion) potential \[V\] found by bisection.
    ///
    /// The corrosion potential E_corr satisfies `i_a(E_corr) = i_c(E_corr)`.
    /// Bisects the interval \[e_low, e_high\] to find the zero of `i_net`.
    ///
    /// Returns `None` if sign does not change on the interval.
    pub fn corrosion_potential(&self, e_low: f64, e_high: f64, tol: f64) -> Option<f64> {
        let f_low = self.net_current(e_low);
        let f_high = self.net_current(e_high);
        if f_low * f_high > 0.0 {
            return None;
        }
        let mut lo = e_low;
        let mut hi = e_high;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if (hi - lo) < tol {
                return Some(mid);
            }
            let f_mid = self.net_current(mid);
            if f_mid * self.net_current(lo) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }
    /// Corrosion current density \[A/m²\] at the corrosion potential.
    ///
    /// Returns `None` if corrosion potential cannot be found.
    pub fn corrosion_current(&self, e_low: f64, e_high: f64, tol: f64) -> Option<f64> {
        let e_corr = self.corrosion_potential(e_low, e_high, tol)?;
        Some(self.anodic_current(e_corr))
    }
}
/// Nernst diffusion layer model.
///
/// Models mass-transport limitations across a stagnant diffusion layer.
pub struct DiffusionLayer {
    /// Diffusion layer thickness δ (m)
    pub thickness: f64,
    /// Diffusivity D (m²/s)
    pub diffusivity: f64,
    /// Bulk concentration C_bulk (mol/m³)
    pub bulk_concentration: f64,
}
impl DiffusionLayer {
    /// Create a new `DiffusionLayer`.
    pub fn new(thickness: f64, d: f64, c_bulk: f64) -> Self {
        Self {
            thickness,
            diffusivity: d,
            bulk_concentration: c_bulk,
        }
    }
    /// Limiting current density.
    ///
    /// `j_L = n · F · D · C_bulk / δ`
    ///
    /// # Arguments
    /// * `n_electrons` — number of electrons transferred per mole of reactant
    /// * `f_const`     — Faraday constant (C/mol)
    /// * `area`        — electrode area (m²)
    pub fn limiting_current(&self, n_electrons: u32, f_const: f64, area: f64) -> f64 {
        (n_electrons as f64) * f_const * self.diffusivity * self.bulk_concentration * area
            / self.thickness
    }
    /// Concentration at the electrode surface.
    ///
    /// `C_s = C_bulk − j · δ / (n · F · D)`
    ///
    /// # Arguments
    /// * `current_density` — current density (A/m²)
    /// * `n_electrons`     — electrons per mole
    /// * `f_const`         — Faraday constant (C/mol)
    pub fn concentration_at_surface(
        &self,
        current_density: f64,
        n_electrons: u32,
        f_const: f64,
    ) -> f64 {
        self.bulk_concentration
            - current_density * self.thickness / ((n_electrons as f64) * f_const * self.diffusivity)
    }
    /// Diffusion overpotential.
    ///
    /// `η_d = (RT/F) · ln(C_s / C_bulk)`
    ///
    /// # Arguments
    /// * `c_surface`  — surface concentration (mol/m³)
    /// * `temp`       — temperature (K) — kept for API symmetry; used implicitly via `f_over_rt`
    /// * `f_over_rt`  — F / (R·T) (1/V)
    pub fn diffusion_overpotential(&self, c_surface: f64, _temp: f64, f_over_rt: f64) -> f64 {
        (c_surface / self.bulk_concentration).ln() / f_over_rt
    }
}
/// Peukert's law model for battery discharge capacity.
///
/// Peukert's equation relates the capacity of a battery to the discharge rate:
///
/// `t = C_peukert / I^k`
///
/// where `t` is the discharge time \[h\], `I` is the discharge current \[A\],
/// `k` is the Peukert exponent (1 = ideal, >1 for real batteries, typically 1.1-1.3),
/// and `C_peukert` is the Peukert capacity constant \[Ah·A^(k-1)\].
///
/// The actual available capacity at discharge current I is:
/// `Q_actual = C_nominal * (I_nominal / I)^(k-1)`
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PeukertModel {
    /// Nominal capacity \[Ah\] at the rated discharge rate
    pub capacity_nominal_ah: f64,
    /// Rated discharge current \[A\] at which nominal capacity is measured
    pub i_nominal_a: f64,
    /// Peukert exponent k (dimensionless, typically 1.05–1.4)
    pub peukert_exponent: f64,
}
impl PeukertModel {
    /// Create a new Peukert battery model.
    pub fn new(capacity_nominal_ah: f64, i_nominal_a: f64, peukert_exponent: f64) -> Self {
        Self {
            capacity_nominal_ah,
            i_nominal_a,
            peukert_exponent,
        }
    }
    /// Lead-acid battery typical values: C=100 Ah, I_nom=5 A (C/20 rate), k=1.2
    pub fn lead_acid_100ah() -> Self {
        Self::new(100.0, 5.0, 1.2)
    }
    /// Li-ion battery: close to ideal (k ≈ 1.05)
    pub fn li_ion_10ah() -> Self {
        Self::new(10.0, 0.5, 1.05)
    }
    /// Peukert capacity constant \[Ah · A^(k-1)\].
    ///
    /// `C_p = Q_nom * I_nom^(k-1)`
    pub fn peukert_constant(&self) -> f64 {
        self.capacity_nominal_ah * self.i_nominal_a.powf(self.peukert_exponent - 1.0)
    }
    /// Available capacity \[Ah\] at discharge current `i_a` \[A\].
    ///
    /// `Q = Q_nom * (I_nom / I)^(k-1)`
    pub fn available_capacity_ah(&self, current_a: f64) -> f64 {
        if current_a <= 0.0 {
            return self.capacity_nominal_ah;
        }
        self.capacity_nominal_ah * (self.i_nominal_a / current_a).powf(self.peukert_exponent - 1.0)
    }
    /// Discharge time \[h\] at constant current `i_a` \[A\].
    ///
    /// `t = C_p / I^k`
    pub fn discharge_time_h(&self, current_a: f64) -> f64 {
        if current_a <= 0.0 {
            return f64::INFINITY;
        }
        self.peukert_constant() / current_a.powf(self.peukert_exponent)
    }
    /// C-rate (multiple of nominal capacity drained per hour).
    ///
    /// `C-rate = I / Q_nom`
    pub fn c_rate(&self, current_a: f64) -> f64 {
        current_a / self.capacity_nominal_ah
    }
    /// Capacity fade factor at C-rate compared to rated: Q_actual / Q_nominal.
    pub fn capacity_fade_factor(&self, current_a: f64) -> f64 {
        self.available_capacity_ah(current_a) / self.capacity_nominal_ah
    }
}
/// Butler-Volmer electrode kinetics.
///
/// Models the current–overpotential relationship at an electrode surface via
/// the Butler-Volmer equation:
/// `j = j0 * (exp(αa·F·η/RT) − exp(−αc·F·η/RT))`.
pub struct ButlerVolmerKinetics {
    /// Exchange current density j₀ \[A/m²\]
    pub i0: f64,
    /// Anodic transfer coefficient αa (dimensionless)
    pub alpha_a: f64,
    /// Cathodic transfer coefficient αc (dimensionless)
    pub alpha_c: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}
impl ButlerVolmerKinetics {
    /// Create a new `ButlerVolmerKinetics`.
    pub fn new(i0: f64, alpha_a: f64, alpha_c: f64, temp: f64) -> Self {
        Self {
            i0,
            alpha_a,
            alpha_c,
            temperature: temp,
        }
    }
    /// Compute F/(R·T) \[1/V\] at the stored temperature.
    pub fn f_over_rt(&self) -> f64 {
        FARADAY / (GAS_CONSTANT * self.temperature)
    }
    /// Current density \[A/m²\] at overpotential `eta` \[V\].
    pub fn current_density(&self, eta: f64) -> f64 {
        let q = self.f_over_rt();
        self.i0 * ((self.alpha_a * q * eta).exp() - (-self.alpha_c * q * eta).exp())
    }
    /// Linearised charge-transfer resistance \[Ω·m²\] near equilibrium (η → 0).
    ///
    /// `R_ct = RT / (i0 · (αa + αc) · F)`
    pub fn charge_transfer_resistance(&self) -> f64 {
        GAS_CONSTANT * self.temperature / (self.i0 * (self.alpha_a + self.alpha_c) * FARADAY)
    }
    /// Anodic Tafel slope \[V/decade\].
    ///
    /// `b_a = ln(10) · RT / (αa · F)`
    pub fn tafel_slope_anodic(&self) -> f64 {
        10.0_f64.ln() / (self.alpha_a * self.f_over_rt())
    }
    /// Cathodic Tafel slope \[V/decade\].
    ///
    /// `b_c = ln(10) · RT / (αc · F)`
    pub fn tafel_slope_cathodic(&self) -> f64 {
        10.0_f64.ln() / (self.alpha_c * self.f_over_rt())
    }
}
/// Electrochemical Impedance Spectroscopy (EIS) model.
///
/// Models a simple Randles circuit:
/// Z = R_s + (R_ct ∥ C_dl) + Warburg diffusion element
///
/// where:
/// - R_s = solution (ohmic) resistance
/// - R_ct = charge transfer resistance
/// - C_dl = double-layer capacitance
/// - W = Warburg coefficient
#[derive(Debug, Clone)]
pub struct ImpedanceModel {
    /// Solution resistance R_s \[Ω\]
    pub r_solution: f64,
    /// Charge transfer resistance R_ct \[Ω\]
    pub r_ct: f64,
    /// Double-layer capacitance C_dl \[F\]
    pub c_dl: f64,
    /// Warburg coefficient σ \[Ω·s^{-0.5}\]
    pub warburg_sigma: f64,
}
impl ImpedanceModel {
    /// Create a new Randles circuit impedance model.
    pub fn new(r_solution: f64, r_ct: f64, c_dl: f64, warburg_sigma: f64) -> Self {
        Self {
            r_solution,
            r_ct,
            c_dl,
            warburg_sigma,
        }
    }
    /// Angular frequency \[rad/s\] from frequency f \[Hz\].
    pub fn angular_frequency(f_hz: f64) -> f64 {
        2.0 * std::f64::consts::PI * f_hz
    }
    /// Real part of impedance Z' at angular frequency ω \[Ω\].
    ///
    /// Z' = R_s + (R_ct + σ/√ω) / ((1 + C_dl·ω·σ·√ω)² + (C_dl·ω·R_ct)²)  \[simplified\]
    ///
    /// For the Warburg element: Z_W = σ·(1 - j) / √ω
    pub fn z_real(&self, omega: f64) -> f64 {
        if omega < 1e-15 {
            return self.r_solution + self.r_ct;
        }
        let sqrt_w = omega.sqrt();
        let sigma = self.warburg_sigma;
        let zw_r = sigma / sqrt_w;
        let rw = self.r_ct + zw_r;
        let xc = 1.0 / (omega * self.c_dl);
        let _denom = rw * rw + xc * xc;
        let z_par_r = rw / (1.0 + (rw / xc).powi(2));
        self.r_solution + z_par_r
    }
    /// Imaginary part of impedance Z'' at angular frequency ω \[Ω\].
    ///
    /// Negative for capacitive behaviour.
    pub fn z_imag(&self, omega: f64) -> f64 {
        if omega < 1e-15 {
            return 0.0;
        }
        let sqrt_w = omega.sqrt();
        let sigma = self.warburg_sigma;
        let zw_r = sigma / sqrt_w;
        let zw_i = -sigma / sqrt_w;
        let rw = self.r_ct + zw_r;
        let xw = zw_i;
        let xc = 1.0 / (omega * self.c_dl);
        let denom = 1.0 + (rw / xc).powi(2);

        (xw - rw * rw / xc - xc) / denom
    }
    /// Magnitude |Z| at angular frequency ω \[Ω\].
    pub fn z_magnitude(&self, omega: f64) -> f64 {
        let zr = self.z_real(omega);
        let zi = self.z_imag(omega);
        (zr * zr + zi * zi).sqrt()
    }
    /// Phase angle θ = atan(Z''/Z') \[radians\] at angular frequency ω.
    pub fn phase_angle(&self, omega: f64) -> f64 {
        self.z_imag(omega).atan2(self.z_real(omega))
    }
    /// Nyquist plot point (Z', -Z'') at frequency f \[Hz\].
    pub fn nyquist_point(&self, f_hz: f64) -> (f64, f64) {
        let omega = Self::angular_frequency(f_hz);
        (self.z_real(omega), -self.z_imag(omega))
    }
}
/// Concentration polarisation overpotential model.
///
/// When current flows, the reactant concentration at the electrode surface
/// differs from the bulk. This creates a concentration overpotential:
///
/// `η_conc = (RT/nF) * ln(C_bulk / C_surface)`
///
/// For a simple linear diffusion layer: `C_s = C_bulk * (1 - j/j_L)`.
#[derive(Debug, Clone)]
pub struct ConcentrationPolarisation {
    /// Bulk concentration \[mol/m³\]
    pub c_bulk: f64,
    /// Limiting current density \[A/m²\]
    pub j_limiting: f64,
    /// Number of electrons per reaction
    pub n_electrons: u32,
    /// Temperature \[K\]
    pub temperature: f64,
}
impl ConcentrationPolarisation {
    /// Create a new concentration polarisation model.
    pub fn new(c_bulk: f64, j_limiting: f64, n_electrons: u32, temperature: f64) -> Self {
        Self {
            c_bulk,
            j_limiting,
            n_electrons,
            temperature,
        }
    }
    /// Surface concentration at current density j \[A/m²\].
    ///
    /// `C_s = C_bulk * (1 - j/j_L)`
    pub fn surface_concentration(&self, j: f64) -> f64 {
        (self.c_bulk * (1.0 - j / self.j_limiting)).max(0.0)
    }
    /// Concentration overpotential \[V\] at current density j \[A/m²\].
    ///
    /// `η_conc = -(RT/nF) * ln(1 - j/j_L)`
    pub fn overpotential(&self, j: f64) -> f64 {
        let c_s = self.surface_concentration(j);
        if c_s < 1e-15 {
            return f64::INFINITY;
        }
        let rt_over_nf = GAS_CONSTANT * self.temperature / (self.n_electrons as f64 * FARADAY);
        -rt_over_nf * (c_s / self.c_bulk).ln()
    }
    /// Fraction of limiting current (dimensionless): `j / j_L`.
    pub fn current_fraction(&self, j: f64) -> f64 {
        j / self.j_limiting
    }
}
/// Simplified SEI layer resistance growth model.
///
/// The SEI layer forms during the first charge of a Li-ion cell.
/// Resistance grows with time and temperature via a diffusion-limited mechanism:
///
/// `R_SEI(t) = R_SEI_0 + k_SEI * sqrt(t)`
///
/// where `k_SEI` is the growth rate constant \[Ω·s^{-0.5}\] and depends on temperature
/// via an Arrhenius factor.
#[derive(Debug, Clone)]
pub struct SeiLayer {
    /// Initial SEI resistance \[Ω\]
    pub r_sei_0: f64,
    /// Growth rate constant at T_ref \[Ω·s^{-0.5}\]
    pub k_sei_ref: f64,
    /// Reference temperature \[K\]
    pub t_ref: f64,
    /// Activation energy for SEI growth \[J/mol\]
    pub activation_energy: f64,
}
impl SeiLayer {
    /// Create a new SEI layer model.
    pub fn new(r_sei_0: f64, k_sei_ref: f64, t_ref: f64, activation_energy: f64) -> Self {
        Self {
            r_sei_0,
            k_sei_ref,
            t_ref,
            activation_energy,
        }
    }
    /// SEI growth rate constant at temperature T \[K\].
    pub fn k_sei(&self, temp_k: f64) -> f64 {
        let exponent = -self.activation_energy / GAS_CONSTANT * (1.0 / temp_k - 1.0 / self.t_ref);
        self.k_sei_ref * exponent.exp()
    }
    /// SEI resistance \[Ω\] after time t \[s\] at temperature T \[K\].
    pub fn resistance(&self, t_s: f64, temp_k: f64) -> f64 {
        self.r_sei_0 + self.k_sei(temp_k) * t_s.sqrt()
    }
    /// SEI growth rate dR/dt \[Ω/s\] at time t \[s\] and temperature T \[K\].
    pub fn growth_rate(&self, t_s: f64, temp_k: f64) -> f64 {
        if t_s < 1e-15 {
            return f64::INFINITY;
        }
        0.5 * self.k_sei(temp_k) / t_s.sqrt()
    }
    /// Temperature at which SEI resistance doubles in time `t_s` \[s\].
    ///
    /// Solves: `R_SEI_0 + k_SEI(T) * sqrt(t) = 2 * R_SEI_0`
    pub fn time_to_double_resistance(&self, temp_k: f64) -> f64 {
        if self.r_sei_0 < f64::EPSILON {
            return f64::INFINITY;
        }
        let k = self.k_sei(temp_k);
        (self.r_sei_0 / k).powi(2)
    }
}
/// Faraday's law of electrolysis model for electroplating.
///
/// Mass deposited: `m = (I · t · M) / (n · F)`
/// where M is molar mass \[g/mol\], n is electrons per ion.
#[derive(Debug, Clone)]
pub struct Electroplating {
    /// Molar mass of deposited metal \[g/mol\]
    pub molar_mass: f64,
    /// Number of electrons per ion (valence)
    pub valence: u32,
    /// Current efficiency (0–1, accounts for side reactions)
    pub current_efficiency: f64,
}
impl Electroplating {
    /// Create a new electroplating model.
    pub fn new(molar_mass: f64, valence: u32, current_efficiency: f64) -> Self {
        Self {
            molar_mass,
            valence,
            current_efficiency,
        }
    }
    /// Mass deposited \[g\] for current `i_a` \[A\] over time `t_s` \[s\].
    ///
    /// `m = η_curr · I · t · M / (n · F)`
    pub fn mass_deposited_g(&self, current_a: f64, time_s: f64) -> f64 {
        self.current_efficiency * current_a * time_s * self.molar_mass
            / (self.valence as f64 * FARADAY)
    }
    /// Thickness deposited \[µm\] given anode area `area_m2` \[m²\] and density `rho` \[g/cm³\].
    pub fn thickness_um(
        &self,
        current_a: f64,
        time_s: f64,
        area_m2: f64,
        density_g_cm3: f64,
    ) -> f64 {
        let mass_g = self.mass_deposited_g(current_a, time_s);
        let volume_cm3 = mass_g / density_g_cm3;
        let area_cm2 = area_m2 * 1.0e4;
        let thickness_cm = volume_cm3 / area_cm2;
        thickness_cm * 1.0e4
    }
    /// Required time \[s\] to deposit target mass `m_g` \[g\] at current `i_a` \[A\].
    pub fn time_for_mass(&self, m_g: f64, current_a: f64) -> f64 {
        m_g * self.valence as f64 * FARADAY
            / (self.current_efficiency * current_a * self.molar_mass)
    }
    /// Plating rate \[µm/min\] at current density `j_a_m2` \[A/m²\] and density `rho` \[g/cm³\].
    pub fn plating_rate_um_per_min(&self, j_a_m2: f64, density_g_cm3: f64) -> f64 {
        let j_a_cm2 = j_a_m2 * 1.0e-4;
        let rate_g_cm2_s =
            self.current_efficiency * j_a_cm2 * self.molar_mass / (self.valence as f64 * FARADAY);
        let rate_cm_s = rate_g_cm2_s / density_g_cm3;
        rate_cm_s * 1.0e4 * 60.0
    }
}
/// Simple battery cell model with state-of-charge (SOC) tracking.
pub struct BatteryCell {
    /// Nominal capacity (Ah)
    pub capacity_ah: f64,
    /// Nominal voltage (V)
    pub nominal_voltage: f64,
    /// Internal resistance (Ω)
    pub internal_resistance: f64,
    /// State of charge (0 = empty, 1 = full)
    pub soc: f64,
}
impl BatteryCell {
    /// Create a new fully-charged `BatteryCell`.
    ///
    /// # Arguments
    /// * `capacity`   — capacity (Ah)
    /// * `voltage`    — nominal voltage (V)
    /// * `resistance` — internal resistance (Ω)
    pub fn new(capacity: f64, voltage: f64, resistance: f64) -> Self {
        Self {
            capacity_ah: capacity,
            nominal_voltage: voltage,
            internal_resistance: resistance,
            soc: 1.0,
        }
    }
    /// Open-circuit voltage (linear SOC approximation).
    ///
    /// `V_oc = V_nom · (0.8 + 0.2 · SOC)`
    pub fn open_circuit_voltage(&self) -> f64 {
        self.nominal_voltage * (0.8 + 0.2 * self.soc)
    }
    /// Terminal voltage under load.
    ///
    /// `V = V_oc − I · R_int`  (positive `current_a` = discharge)
    pub fn terminal_voltage(&self, current_a: f64) -> f64 {
        self.open_circuit_voltage() - current_a * self.internal_resistance
    }
    /// Advance the cell by `dt_s` seconds at `current_a` amperes.
    ///
    /// Updates SOC: `SOC -= I · Δt / (capacity · 3600)`
    pub fn discharge(&mut self, current_a: f64, dt_s: f64) {
        self.soc -= current_a * dt_s / (self.capacity_ah * 3600.0);
        self.soc = self.soc.clamp(0.0, 1.0);
    }
    /// Returns `true` when the cell is considered depleted (SOC < 5 %).
    pub fn is_depleted(&self) -> bool {
        self.soc < 0.05
    }
}
/// Galvanic series potential values for common metals in seawater \[V vs SHE\].
///
/// Approximate values; actual values depend on environment and surface condition.
#[allow(dead_code)]
pub struct GalvanicSeriesEntry {
    /// Metal name
    pub name: &'static str,
    /// Approximate potential in seawater \[V vs SHE\]
    pub potential_v: f64,
}
/// Butler-Volmer electrode kinetics model.
///
/// Models the current–overpotential relationship at an electrode surface.
pub struct ElectrodeKinetics {
    /// Exchange current density (A/m²)
    pub exchange_current_density: f64,
    /// Anodic transfer coefficient (dimensionless)
    pub transfer_coefficient_a: f64,
    /// Cathodic transfer coefficient (dimensionless)
    pub transfer_coefficient_c: f64,
    /// Temperature (K)
    pub temperature: f64,
}
impl ElectrodeKinetics {
    /// Create a new `ElectrodeKinetics`.
    ///
    /// # Arguments
    /// * `j0`      — exchange current density (A/m²)
    /// * `alpha_a` — anodic transfer coefficient
    /// * `alpha_c` — cathodic transfer coefficient
    /// * `temp`    — temperature (K)
    pub fn new(j0: f64, alpha_a: f64, alpha_c: f64, temp: f64) -> Self {
        Self {
            exchange_current_density: j0,
            transfer_coefficient_a: alpha_a,
            transfer_coefficient_c: alpha_c,
            temperature: temp,
        }
    }
    /// Butler-Volmer current density.
    ///
    /// `j = j0 * (exp(αa · F · η / RT) − exp(−αc · F · η / RT))`
    ///
    /// # Arguments
    /// * `eta`      — overpotential (V)
    /// * `f_over_rt`— F / (R·T) (1/V)
    pub fn current_density(&self, eta: f64, f_over_rt: f64) -> f64 {
        let j0 = self.exchange_current_density;
        let aa = self.transfer_coefficient_a;
        let ac = self.transfer_coefficient_c;
        j0 * ((aa * f_over_rt * eta).exp() - (-ac * f_over_rt * eta).exp())
    }
    /// Linearised charge-transfer resistance at η = 0.
    ///
    /// `R_ct = RT / (j0 · (αa + αc) · F)` in the symmetric limit.
    /// For the requested simplified form `∂η/∂j|_{η=0} = 1 / (j0 · F/RT)`:
    ///
    /// Returns `1 / (j0 · f_over_rt)` (Ω·m²).
    ///
    /// # Arguments
    /// * `f_over_rt` — F / (R·T) (1/V)
    pub fn linearized_resistance(&self, f_over_rt: f64) -> f64 {
        1.0 / (self.exchange_current_density * f_over_rt)
    }
    /// Anodic Tafel slope.
    ///
    /// `b_a = ln(10) · RT / (αa · F)`
    ///
    /// # Arguments
    /// * `f_over_rt` — F / (R·T) (1/V)
    pub fn tafel_slope_anodic(&self, f_over_rt: f64) -> f64 {
        10.0_f64.ln() / (self.transfer_coefficient_a * f_over_rt)
    }
}
/// Galvanic couple model: two dissimilar metals in electrical contact in an electrolyte.
///
/// The more anodic metal corrodes preferentially. The couple potential and
/// current are determined by the intersection of the polarisation curves.
#[derive(Debug, Clone)]
pub struct GalvanicCouple {
    /// Metal 1 (anode): corrosion potential \[V\]
    pub e_corr_1: f64,
    /// Metal 1 Tafel slope anodic \[V/decade\]
    pub b_a1: f64,
    /// Metal 1 corrosion current \[A\]
    pub i_corr_1: f64,
    /// Metal 2 (cathode): corrosion potential \[V\]
    pub e_corr_2: f64,
    /// Metal 2 Tafel slope cathodic \[V/decade\]
    pub b_c2: f64,
    /// Metal 2 corrosion current \[A\]
    pub i_corr_2: f64,
}
impl GalvanicCouple {
    /// Create a new galvanic couple model.
    pub fn new(
        e_corr_1: f64,
        b_a1: f64,
        i_corr_1: f64,
        e_corr_2: f64,
        b_c2: f64,
        i_corr_2: f64,
    ) -> Self {
        Self {
            e_corr_1,
            b_a1,
            i_corr_1,
            e_corr_2,
            b_c2,
            i_corr_2,
        }
    }
    /// Anodic current of metal 1 at potential E \[A\]: Tafel approximation.
    pub fn anodic_current_1(&self, e: f64) -> f64 {
        self.i_corr_1 * 10.0_f64.powf((e - self.e_corr_1) / self.b_a1)
    }
    /// Cathodic current of metal 2 at potential E \[A\]: Tafel approximation.
    pub fn cathodic_current_2(&self, e: f64) -> f64 {
        self.i_corr_2 * 10.0_f64.powf(-(e - self.e_corr_2) / self.b_c2)
    }
    /// Couple potential \[V\]: found by bisection where `i_a1 = i_c2`.
    pub fn couple_potential(&self, tol: f64) -> Option<f64> {
        let e_lo = self.e_corr_1.min(self.e_corr_2) - 0.1;
        let e_hi = self.e_corr_1.max(self.e_corr_2) + 0.1;
        let f = |e: f64| self.anodic_current_1(e) - self.cathodic_current_2(e);
        let fl = f(e_lo);
        let fh = f(e_hi);
        if fl * fh > 0.0 {
            return None;
        }
        let mut lo = e_lo;
        let mut hi = e_hi;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if (hi - lo) < tol {
                return Some(mid);
            }
            if f(mid) * f(lo) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }
    /// Galvanic current \[A\] at the couple potential.
    pub fn galvanic_current(&self, tol: f64) -> Option<f64> {
        let e_couple = self.couple_potential(tol)?;
        Some(self.anodic_current_1(e_couple))
    }
}
/// Li-ion battery degradation model.
///
/// Models capacity fade and resistance growth over cycling using
/// empirical power-law and SEI-layer growth models.
///
/// Reference: Plett, "Battery Management Systems", 2015.
#[derive(Debug, Clone)]
pub struct LiIonDegradation {
    /// Initial capacity \[Ah\].
    pub capacity_init: f64,
    /// Current capacity \[Ah\].
    pub capacity: f64,
    /// Initial internal resistance \[Ω\].
    pub resistance_init: f64,
    /// Current internal resistance \[Ω\].
    pub resistance: f64,
    /// Total charge throughput \[Ah\].
    pub throughput_ah: f64,
    /// Cycle count.
    pub n_cycles: u64,
    /// Calendar age \[days\].
    pub calendar_days: f64,
}
impl LiIonDegradation {
    /// Create a new, fresh battery cell.
    pub fn new(capacity: f64, resistance: f64) -> Self {
        Self {
            capacity_init: capacity,
            capacity,
            resistance_init: resistance,
            resistance,
            throughput_ah: 0.0,
            n_cycles: 0,
            calendar_days: 0.0,
        }
    }
    /// State of health based on capacity fade (SOH = Q/Q_init).
    pub fn soh_capacity(&self) -> f64 {
        self.capacity / self.capacity_init
    }
    /// State of health based on resistance growth (SOH_R = R_init/R).
    pub fn soh_resistance(&self) -> f64 {
        self.resistance_init / self.resistance
    }
    /// Update capacity fade after `delta_ah` Ah of throughput.
    ///
    /// Simple power-law model: Q(Ah) = Q0 * (1 - k_ah * Ah^0.5)
    /// where k_ah is the fade coefficient per √Ah.
    pub fn update_cycle_fade(&mut self, delta_ah: f64, k_fade: f64) {
        self.throughput_ah += delta_ah;
        self.n_cycles += 1;
        let fade = k_fade * self.throughput_ah.sqrt();
        self.capacity = self.capacity_init * (1.0 - fade).max(0.0);
    }
    /// Update SEI-layer resistance growth.
    ///
    /// R(t) = R0 * (1 + k_sei * sqrt(t))  \[calendar aging\]
    pub fn update_calendar_aging(&mut self, delta_days: f64, k_sei: f64) {
        self.calendar_days += delta_days;
        self.resistance = self.resistance_init * (1.0 + k_sei * self.calendar_days.sqrt());
    }
    /// Check if the battery has reached end of life (SOH < 80%).
    pub fn is_end_of_life(&self) -> bool {
        self.soh_capacity() < 0.8
    }
    /// Remaining useful life estimate \[cycles\] using a linear fade model.
    ///
    /// Extrapolates from current fade rate to SOH = 0.8.
    pub fn estimated_remaining_cycles(&self, k_fade: f64) -> Option<u64> {
        if k_fade < 1e-30 {
            return None;
        }
        let ah_eol = (0.2 / k_fade).powi(2);
        if ah_eol <= self.throughput_ah {
            return Some(0);
        }
        let remaining_ah = ah_eol - self.throughput_ah;
        if self.n_cycles == 0 {
            return None;
        }
        let ah_per_cycle = self.throughput_ah / self.n_cycles as f64;
        Some((remaining_ah / ah_per_cycle) as u64)
    }
}
/// Gouy-Chapman-Stern (GCS) model for the electrical double layer.
///
/// Models the differential capacitance of the diffuse layer as a function
/// of potential difference across the Helmholtz layer.
///
/// The Debye length κ⁻¹ = sqrt(ε·RT / (2·n0·F²))
/// where n0 is bulk ion concentration \[mol/m³\] and ε is permittivity \[F/m\].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DoubleLayerCapacitance {
    /// Bulk ion concentration \[mol/m³\] (symmetric 1:1 electrolyte)
    pub concentration_mol_m3: f64,
    /// Permittivity of solvent \[F/m\] (water ≈ 78.5 × ε_0)
    pub permittivity: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}
impl DoubleLayerCapacitance {
    /// Create a new double-layer model.
    pub fn new(concentration_mol_m3: f64, permittivity: f64, temperature: f64) -> Self {
        Self {
            concentration_mol_m3,
            permittivity,
            temperature,
        }
    }
    /// Aqueous KCl at 0.1 mol/L and 25°C.
    pub fn kcl_01_mol_l() -> Self {
        Self::new(100.0, 78.5 * 8.854e-12, 298.15)
    }
    /// Debye length \[m\]: κ⁻¹ = sqrt(ε·R·T / (2·c·F²))
    pub fn debye_length(&self) -> f64 {
        let numerator = self.permittivity * GAS_CONSTANT * self.temperature;
        let denominator = 2.0 * self.concentration_mol_m3 * FARADAY * FARADAY;
        (numerator / denominator).sqrt()
    }
    /// Diffuse-layer capacitance per unit area \[F/m²\] at the potential of zero charge.
    ///
    /// `C_d = ε / κ⁻¹`
    pub fn capacitance_at_pzc(&self) -> f64 {
        self.permittivity / self.debye_length()
    }
    /// Differential capacitance \[F/m²\] at potential ψ \[V\] (linearised GC model).
    ///
    /// `C(ψ) = ε * κ * cosh(F·ψ / (2·R·T))`
    pub fn differential_capacitance(&self, psi: f64) -> f64 {
        let kappa = 1.0 / self.debye_length();
        let arg = FARADAY * psi / (2.0 * GAS_CONSTANT * self.temperature);
        self.permittivity * kappa * arg.cosh()
    }
}
/// Corrosion kinetics from Tafel polarisation data.
///
/// The corrosion rate is derived from the corrosion current density `i_corr`
/// via: `CR [mm/year] = (i_corr · M_w) / (n · F · ρ) · K`
/// where `K = 3.27 × 10⁻³` (unit conversion constant for SI inputs).
#[derive(Debug, Clone)]
pub struct CorrosionRate {
    /// Corrosion current density i_corr \[A/m²\]
    pub i_corr: f64,
    /// Molar mass of the metal \[g/mol\]
    pub molar_mass: f64,
    /// Number of electrons in the oxidation reaction
    pub n_electrons: u32,
    /// Density of the metal \[g/cm³\]
    pub density_g_cm3: f64,
}
impl CorrosionRate {
    /// Create a new corrosion rate model.
    pub fn new(i_corr: f64, molar_mass: f64, n_electrons: u32, density_g_cm3: f64) -> Self {
        Self {
            i_corr,
            molar_mass,
            n_electrons,
            density_g_cm3,
        }
    }
    /// Corrosion rate \[mm/year\].
    ///
    /// Uses the standard electrochemical formula with unit conversion constant
    /// `K = 3.27 × 10⁻³` mm·g/(µA·cm·year) adapted to SI (A/m²).
    pub fn mm_per_year(&self) -> f64 {
        let i_ua_cm2 = self.i_corr * 0.1;
        0.00327 * i_ua_cm2 * self.molar_mass / ((self.n_electrons as f64) * self.density_g_cm3)
    }
    /// Corrosion current density from Tafel slopes using the Stern-Geary equation.
    ///
    /// `i_corr = B / R_p`  where `B = (b_a · b_c) / (2.303 · (b_a + b_c))`
    ///
    /// # Arguments
    /// * `rp`  — polarisation resistance \[Ω·m²\]
    /// * `b_a` — anodic Tafel slope \[V/decade\]
    /// * `b_c` — cathodic Tafel slope \[V/decade\]
    pub fn from_polarisation_resistance(rp: f64, b_a: f64, b_c: f64) -> f64 {
        let b = (b_a * b_c) / (2.303 * (b_a + b_c));
        b / rp
    }
}
/// Battery equivalent-circuit model with internal resistance and one RC pair.
///
/// Topology: `V_oc − R0 − (R1 ∥ C1) = V_terminal`.
/// The RC pair captures diffusion / charge-transfer relaxation dynamics.
#[derive(Debug, Clone)]
pub struct BatteryModel {
    /// Nominal capacity \[Ah\]
    pub capacity_ah: f64,
    /// State of charge (0 = empty, 1 = full)
    pub soc: f64,
    /// Series (ohmic) resistance R₀ \[Ω\]
    pub r0: f64,
    /// RC branch resistance R₁ \[Ω\]
    pub r1: f64,
    /// RC branch capacitance C₁ \[F\]
    pub c1: f64,
    /// Voltage across the RC branch (state variable) \[V\]
    pub u_rc: f64,
}
impl BatteryModel {
    /// Create a fully charged battery model.
    ///
    /// # Arguments
    /// * `capacity_ah` — capacity \[Ah\]
    /// * `r0`          — series resistance \[Ω\]
    /// * `r1`          — RC branch resistance \[Ω\]
    /// * `c1`          — RC branch capacitance \[F\]
    pub fn new(capacity_ah: f64, r0: f64, r1: f64, c1: f64) -> Self {
        Self {
            capacity_ah,
            soc: 1.0,
            r0,
            r1,
            c1,
            u_rc: 0.0,
        }
    }
    /// Open-circuit voltage from a piecewise-linear SOC-OCV curve.
    ///
    /// Uses a simplified linear model: `V_oc = 3.0 + 1.2 · SOC` (LFP-like).
    pub fn open_circuit_voltage(&self) -> f64 {
        3.0 + 1.2 * self.soc
    }
    /// Terminal voltage under load \[V\].
    ///
    /// `V = V_oc − R0 · I − U_rc`
    pub fn terminal_voltage(&self, current_a: f64) -> f64 {
        self.open_circuit_voltage() - self.r0 * current_a - self.u_rc
    }
    /// Advance simulation by `dt_s` seconds at current `current_a` \[A\].
    ///
    /// Updates SOC and RC branch voltage via first-order Euler integration.
    pub fn step(&mut self, current_a: f64, dt_s: f64) {
        let tau = self.r1 * self.c1;
        self.u_rc += (self.r1 * current_a - self.u_rc) / tau * dt_s;
        self.soc -= current_a * dt_s / (self.capacity_ah * 3600.0);
        self.soc = self.soc.clamp(0.0, 1.0);
    }
    /// Returns `true` when the cell is considered depleted (SOC < 5 %).
    pub fn is_depleted(&self) -> bool {
        self.soc < 0.05
    }
    /// Instantaneous power \[W\] (positive = discharge).
    pub fn power(&self, current_a: f64) -> f64 {
        self.terminal_voltage(current_a) * current_a
    }
}
/// Hydrogen PEM fuel cell polarisation curve model.
///
/// Total overpotential = activation loss + ohmic loss + concentration loss.
#[derive(Debug, Clone)]
pub struct FuelCellStack {
    /// Reversible (Nernst) cell voltage \[V\]
    pub e_rev: f64,
    /// Exchange current density i₀ \[A/m²\]
    pub i0: f64,
    /// Tafel slope for activation loss \[V/decade\] (positive value)
    pub tafel_slope: f64,
    /// Ohmic area-specific resistance \[Ω·m²\]
    pub r_ohmic: f64,
    /// Limiting current density i_L \[A/m²\]
    pub i_lim: f64,
    /// Concentration loss empirical coefficient m \[V\]
    pub m_conc: f64,
    /// Concentration loss empirical coefficient n \[m²/A\]
    pub n_conc: f64,
}
impl FuelCellStack {
    /// Create a new fuel cell stack model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        e_rev: f64,
        i0: f64,
        tafel_slope: f64,
        r_ohmic: f64,
        i_lim: f64,
        m_conc: f64,
        n_conc: f64,
    ) -> Self {
        Self {
            e_rev,
            i0,
            tafel_slope,
            r_ohmic,
            i_lim,
            m_conc,
            n_conc,
        }
    }
    /// Activation overpotential \[V\] at current density `j` \[A/m²\].
    ///
    /// `η_act = tafel_slope · log10(j / i0)`
    pub fn activation_loss(&self, j: f64) -> f64 {
        if j <= 0.0 || j < self.i0 {
            return 0.0;
        }
        self.tafel_slope * (j / self.i0).log10()
    }
    /// Ohmic overpotential \[V\] at current density `j` \[A/m²\].
    ///
    /// `η_ohm = R_ohmic · j`
    pub fn ohmic_loss(&self, j: f64) -> f64 {
        self.r_ohmic * j
    }
    /// Concentration overpotential \[V\] at current density `j` \[A/m²\].
    ///
    /// `η_conc = m · exp(n · j)`
    pub fn concentration_loss(&self, j: f64) -> f64 {
        self.m_conc * (self.n_conc * j).exp()
    }
    /// Cell terminal voltage \[V\] at current density `j` \[A/m²\].
    ///
    /// `V = E_rev − η_act − η_ohm − η_conc`
    ///
    /// Returns 0.0 if j ≥ i_lim (cell has stalled).
    pub fn cell_voltage(&self, j: f64) -> f64 {
        if j >= self.i_lim {
            return 0.0;
        }
        let v =
            self.e_rev - self.activation_loss(j) - self.ohmic_loss(j) - self.concentration_loss(j);
        v.max(0.0)
    }
    /// Power density \[W/m²\] at current density `j` \[A/m²\].
    pub fn power_density(&self, j: f64) -> f64 {
        j * self.cell_voltage(j)
    }
    /// Efficiency relative to reversible voltage.
    ///
    /// `η = V_cell / E_rev`
    pub fn efficiency(&self, j: f64) -> f64 {
        self.cell_voltage(j) / self.e_rev
    }
}
/// Electrolyte ionic conductivity model.
///
/// Provides the Kohlrausch law for strong electrolytes and the
/// Casteel-Amis empirical model for concentrated solutions.
#[derive(Debug, Clone)]
pub struct ElectrolyteConductivity {
    /// Limiting molar conductivity Λ₀ \[S·m²/mol\] (at infinite dilution).
    pub lambda_0: f64,
    /// Kohlrausch coefficient B \[S·m²·mol^{-3/2}\].
    pub kohlrausch_b: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}
impl ElectrolyteConductivity {
    /// Create a new electrolyte conductivity model.
    pub fn new(lambda_0: f64, kohlrausch_b: f64, temperature: f64) -> Self {
        Self {
            lambda_0,
            kohlrausch_b,
            temperature,
        }
    }
    /// Kohlrausch law: Λ(c) = Λ₀ - B·√c \[S·m²/mol\].
    ///
    /// Valid for dilute solutions (c < ~0.01 mol/L).
    pub fn molar_conductivity(&self, conc: f64) -> f64 {
        (self.lambda_0 - self.kohlrausch_b * conc.sqrt()).max(0.0)
    }
    /// Specific conductivity κ = Λ·c \[S/m\].
    pub fn specific_conductivity(&self, conc: f64) -> f64 {
        self.molar_conductivity(conc) * conc
    }
    /// Resistivity ρ = 1/κ \[Ω·m\].
    pub fn resistivity(&self, conc: f64) -> f64 {
        let kappa = self.specific_conductivity(conc);
        if kappa < f64::EPSILON {
            return f64::INFINITY;
        }
        1.0 / kappa
    }
    /// Temperature correction using Arrhenius model.
    ///
    /// κ(T) = κ(T_ref) * exp(-E_a/R * (1/T - 1/T_ref))
    ///
    /// Returns the conductivity at temperature T given conductivity at T_ref.
    pub fn temperature_corrected_conductivity(
        &self,
        kappa_ref: f64,
        t_ref: f64,
        activation_energy: f64,
    ) -> f64 {
        const R: f64 = 8.314_462_618;
        let exponent = -activation_energy / R * (1.0 / self.temperature - 1.0 / t_ref);
        kappa_ref * exponent.exp()
    }
    /// Transference number for the cation (fraction of current carried by cation).
    ///
    /// Uses the ratio of limiting ionic conductivities.
    /// t+ = λ+ / (λ+ + λ-)
    pub fn transference_number_cation(lambda_plus: f64, lambda_minus: f64) -> f64 {
        lambda_plus / (lambda_plus + lambda_minus)
    }
}
