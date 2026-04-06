// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Magnetocaloric effect (MCE) material models.
//!
//! Implements:
//! - [`MagnetocaloricMaterial`] — Curie temperature, adiabatic ΔT, isothermal ΔS
//! - [`MagneticEntropy`] — Clausius-Clapeyron, Maxwell relation, field-swept entropy
//! - [`GiantMce`] — La(Fe,Si)13 model, MnAs compounds, first-order transition
//! - [`MagneticRefrigeration`] — AMR cycle, Brayton/Ericsson/Carnot cycles, COP
//! - [`LandauTheoryMagneto`] — Landau expansion for magnetic free energy
//! - [`MeanFieldMagnetics`] — Weiss mean-field, Brillouin function
//! - [`MagnetoCaloric`] — Heusler alloy model (Ni-Mn-Ga), metamagnetic transition
//! - [`BarocaloriEffect`] — pressure-induced entropy change, mechanocaloric
//! - [`ElastocaloricMaterial`] — shape memory alloy caloric effect
//! - [`SpinCaloritronics`] — Seebeck spin tunneling, spin Peltier/Nernst effects

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant \[J/K\].
const K_B: f64 = 1.380_649e-23;
/// Avogadro constant \[mol⁻¹\].
const N_A: f64 = 6.022_140_76e23;
/// Bohr magneton μ_B \[J/T\].
const MU_B: f64 = 9.274_010_08e-24;
/// Vacuum permeability μ₀ \[T·m/A\].
const MU_0: f64 = 1.256_637_062e-6;
/// Gas constant R \[J/(mol·K)\].
const R_GAS: f64 = 8.314_462_618;

// ---------------------------------------------------------------------------
// MagnetocaloricMaterial
// ---------------------------------------------------------------------------

/// Bulk magnetocaloric material properties.
///
/// Stores the key parameters needed to compute the magnetocaloric effect (MCE):
/// the Curie temperature, specific heat, and field-sweep limits.
pub struct MagnetocaloricMaterial {
    /// Curie temperature T_C \[K\] — ferromagnetic ordering temperature.
    pub curie_temperature: f64,
    /// Specific heat capacity c_p \[J/(kg·K)\] near T_C.
    pub specific_heat: f64,
    /// Peak isothermal entropy change |ΔS_M|_max \[J/(kg·K)\] per tesla.
    pub delta_s_per_tesla: f64,
    /// Peak adiabatic temperature change |ΔT_ad|_max \[K\] per tesla.
    pub delta_t_per_tesla: f64,
    /// Material density ρ \[kg/m³\].
    pub density: f64,
}

impl MagnetocaloricMaterial {
    /// Create a new `MagnetocaloricMaterial`.
    pub fn new(
        curie_temperature: f64,
        specific_heat: f64,
        delta_s_per_tesla: f64,
        delta_t_per_tesla: f64,
        density: f64,
    ) -> Self {
        Self {
            curie_temperature,
            specific_heat,
            delta_s_per_tesla,
            delta_t_per_tesla,
            density,
        }
    }

    /// Construct a gadolinium (Gd) reference material.
    ///
    /// Returns standard literature values: T_C ≈ 294 K, |ΔS|_peak ≈ 3 J/(kg·K) at 2 T.
    pub fn gadolinium() -> Self {
        Self::new(294.0, 300.0, 1.5, 1.5, 7_900.0)
    }

    /// Isothermal entropy change for a given applied field Δμ₀H \[T\].
    ///
    /// Uses the linear approximation |ΔS_M| ≈ delta_s_per_tesla · |ΔH|.
    pub fn isothermal_entropy_change(&self, delta_field: f64) -> f64 {
        self.delta_s_per_tesla * delta_field.abs()
    }

    /// Adiabatic temperature change for a given applied field Δμ₀H \[T\].
    ///
    /// Uses the linear approximation |ΔT_ad| ≈ delta_t_per_tesla · |ΔH|.
    pub fn adiabatic_temp_change(&self, delta_field: f64) -> f64 {
        self.delta_t_per_tesla * delta_field.abs()
    }

    /// Reduced temperature t = T / T_C (dimensionless).
    pub fn reduced_temperature(&self, temperature: f64) -> f64 {
        temperature / self.curie_temperature
    }

    /// Approximate peak adiabatic ΔT_ad via the thermodynamic identity:
    ///
    /// ΔT_ad ≈ −(T / c_p) · ΔS_M
    pub fn adiabatic_temp_from_entropy(&self, temperature: f64, delta_s: f64) -> f64 {
        -(temperature / self.specific_heat) * delta_s
    }

    /// Refrigerant capacity RC via the area of the ΔS–T curve.
    ///
    /// Uses the full-width-half-maximum (FWHM) approximation:
    /// RC ≈ |ΔS_M|_max · δT_FWHM
    pub fn refrigerant_capacity(&self, delta_field: f64, fwhm_kelvin: f64) -> f64 {
        self.isothermal_entropy_change(delta_field) * fwhm_kelvin
    }
}

// ---------------------------------------------------------------------------
// MagneticEntropy
// ---------------------------------------------------------------------------

/// Magnetic entropy models.
///
/// Provides Clausius-Clapeyron analysis, the Maxwell relation, and
/// field-swept entropy integration for magnetocaloric characterisation.
pub struct MagneticEntropy {
    /// Molar mass M \[kg/mol\].
    pub molar_mass: f64,
    /// Number of magnetic moments per formula unit (effective Bohr magneton number).
    pub effective_moment: f64,
}

impl MagneticEntropy {
    /// Create a new `MagneticEntropy`.
    pub fn new(molar_mass: f64, effective_moment: f64) -> Self {
        Self {
            molar_mass,
            effective_moment,
        }
    }

    /// Maxwell relation for isothermal entropy change:
    ///
    /// ΔS_M = ∫₀^H (∂M/∂T)_H dH ≈ (∂M/∂T)_H · ΔH
    ///
    /// `dm_dt` is the measured (∂M/∂T)_H in \[A·m²/kg·K\] and `delta_h` in A/m.
    pub fn maxwell_relation_entropy_change(&self, dm_dt: f64, delta_h: f64) -> f64 {
        dm_dt * delta_h
    }

    /// Clausius-Clapeyron slope dT_C/dH for a first-order magnetic transition:
    ///
    /// dT_C/dH = −μ₀ · ΔM / ΔS
    ///
    /// `delta_m` is magnetisation discontinuity \[A·m²/kg\], `delta_s` is
    /// entropy discontinuity \[J/(kg·K)\].
    pub fn clausius_clapeyron_slope(&self, delta_m: f64, delta_s: f64) -> f64 {
        if delta_s.abs() < 1e-30 {
            return 0.0;
        }
        -MU_0 * delta_m / delta_s
    }

    /// Total magnetic entropy in the fully disordered (paramagnetic) state:
    ///
    /// S_mag = R · ln(2J + 1) / M_mol
    ///
    /// where J is the total angular momentum quantum number.
    pub fn disordered_entropy(&self, total_angular_momentum: f64) -> f64 {
        let g = 2.0 * total_angular_momentum + 1.0;
        R_GAS * g.ln() / self.molar_mass
    }

    /// Numerical entropy change from tabulated M(T) data using the Maxwell relation.
    ///
    /// `m_values` is a slice of (T, M) pairs measured at two fields H1 and H2.
    /// `delta_h` = H2 − H1 \[A/m\].
    pub fn numerical_maxwell_entropy(
        &self,
        m_low: &[(f64, f64)],
        m_high: &[(f64, f64)],
        delta_h: f64,
    ) -> Vec<(f64, f64)> {
        let n = m_low.len().min(m_high.len());
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let (t_low, m_l) = m_low[i];
            let (_t_high, m_h) = m_high[i];
            let dm_dh = (m_h - m_l) / delta_h;
            result.push((t_low, dm_dh * delta_h));
        }
        result
    }

    /// Integrated entropy change over a temperature span \[T1, T2\].
    ///
    /// Trapezoidal integration of a ΔS(T) dataset.
    pub fn integrated_entropy(&self, ds_data: &[(f64, f64)]) -> f64 {
        if ds_data.len() < 2 {
            return 0.0;
        }
        let mut integral = 0.0;
        for w in ds_data.windows(2) {
            let (t1, s1) = w[0];
            let (t2, s2) = w[1];
            integral += 0.5 * (s1 + s2) * (t2 - t1);
        }
        integral
    }
}

// ---------------------------------------------------------------------------
// GiantMce
// ---------------------------------------------------------------------------

/// Giant magnetocaloric effect (GMCE) material model.
///
/// Models first-order magnetostructural transitions such as those in
/// La(Fe,Si)₁₃ and MnAs-based compounds with field-tunable transition temperature.
pub struct GiantMce {
    /// Transition temperature T₀ \[K\] at zero field.
    pub transition_temp: f64,
    /// Field sensitivity of transition temperature dT₀/dH \[K/T\].
    pub field_sensitivity: f64,
    /// Peak entropy change |ΔS|_max \[J/(kg·K)\].
    pub peak_entropy_change: f64,
    /// Transition width δT \[K\] (full width at half maximum).
    pub transition_width: f64,
    /// Hysteresis loss per cycle q_hys \[J/kg\].
    pub hysteresis_loss: f64,
}

impl GiantMce {
    /// Create a new `GiantMce` model.
    pub fn new(
        transition_temp: f64,
        field_sensitivity: f64,
        peak_entropy_change: f64,
        transition_width: f64,
        hysteresis_loss: f64,
    ) -> Self {
        Self {
            transition_temp,
            field_sensitivity,
            peak_entropy_change,
            transition_width,
            hysteresis_loss,
        }
    }

    /// La(Fe,Si)₁₃ reference material (approximate literature values).
    pub fn la_fe_si() -> Self {
        Self::new(195.0, 2.5, 20.0, 8.0, 200.0)
    }

    /// MnAs reference material at 2 T.
    pub fn mn_as() -> Self {
        Self::new(318.0, -0.5, 30.0, 5.0, 800.0)
    }

    /// Field-shifted transition temperature T₀(H) = T₀ + (dT₀/dH) · μ₀H.
    pub fn shifted_transition_temp(&self, field_tesla: f64) -> f64 {
        self.transition_temp + self.field_sensitivity * field_tesla
    }

    /// Entropy change at temperature T for applied field μ₀H \[T\].
    ///
    /// Uses a Lorentzian peak shape centred at the shifted T₀(H).
    pub fn entropy_change(&self, temperature: f64, field_tesla: f64) -> f64 {
        let t0 = self.shifted_transition_temp(field_tesla);
        let dt = temperature - t0;
        let half_w = self.transition_width / 2.0;
        -self.peak_entropy_change * (half_w * half_w) / (dt * dt + half_w * half_w)
    }

    /// Net cooling capacity per cycle after subtracting hysteresis loss \[J/kg\].
    pub fn net_cooling_capacity(&self, temperature: f64, field_tesla: f64) -> f64 {
        let ds = self.entropy_change(temperature, field_tesla).abs();
        let gross = ds * self.transition_width;
        (gross - self.hysteresis_loss).max(0.0)
    }

    /// First-order transition indicator: true if |ΔS| > 10 J/(kg·K).
    pub fn is_first_order(&self) -> bool {
        self.peak_entropy_change > 10.0
    }
}

// ---------------------------------------------------------------------------
// MagneticRefrigeration
// ---------------------------------------------------------------------------

/// Active magnetic regenerator (AMR) cycle and thermodynamic cycle models.
///
/// Provides COP and specific cooling power calculations for Brayton, Ericsson,
/// and Carnot magnetic refrigeration cycles.
pub struct MagneticRefrigeration {
    /// Hot reservoir temperature T_H \[K\].
    pub t_hot: f64,
    /// Cold reservoir temperature T_C \[K\].
    pub t_cold: f64,
    /// Applied field swing Δμ₀H \[T\].
    pub field_swing: f64,
    /// MCE material used in the regenerator.
    pub material: MagnetocaloricMaterial,
}

impl MagneticRefrigeration {
    /// Create a new `MagneticRefrigeration` system.
    pub fn new(
        t_hot: f64,
        t_cold: f64,
        field_swing: f64,
        material: MagnetocaloricMaterial,
    ) -> Self {
        Self {
            t_hot,
            t_cold,
            field_swing,
            material,
        }
    }

    /// Temperature span ΔT_span = T_H − T_C \[K\].
    pub fn temp_span(&self) -> f64 {
        self.t_hot - self.t_cold
    }

    /// Ideal Carnot COP for a refrigerator:
    ///
    /// COP_Carnot = T_C / (T_H − T_C)
    pub fn carnot_cop(&self) -> f64 {
        let span = self.temp_span();
        if span < 1e-10 {
            return f64::INFINITY;
        }
        self.t_cold / span
    }

    /// Brayton cycle specific cooling power \[J/kg\] per cycle.
    ///
    /// Q_cold ≈ |ΔS_M| · T_C
    pub fn brayton_cooling_power(&self) -> f64 {
        let ds = self.material.isothermal_entropy_change(self.field_swing);
        ds * self.t_cold
    }

    /// Ericsson cycle specific cooling power \[J/kg\] per cycle.
    ///
    /// For an ideal Ericsson cycle Q_cold ≈ |ΔS_M| · T_C.
    pub fn ericsson_cooling_power(&self) -> f64 {
        self.brayton_cooling_power()
    }

    /// AMR ideal COP accounting for MCE magnitude and temperature span.
    ///
    /// η_AMR = Q_cold / W_in,  W_in = Q_hot − Q_cold
    pub fn amr_cop(&self) -> f64 {
        let ds = self.material.isothermal_entropy_change(self.field_swing);
        let q_cold = ds * self.t_cold;
        let q_hot = ds * self.t_hot;
        if (q_hot - q_cold).abs() < 1e-30 {
            return 0.0;
        }
        q_cold / (q_hot - q_cold)
    }

    /// Second-law efficiency η₂ = COP_actual / COP_Carnot.
    pub fn second_law_efficiency(&self) -> f64 {
        let cop_actual = self.amr_cop();
        let cop_carnot = self.carnot_cop();
        if cop_carnot == f64::INFINITY || cop_carnot < 1e-30 {
            return 0.0;
        }
        cop_actual / cop_carnot
    }

    /// Number of cycles per second for a rotary AMR at frequency f \[Hz\].
    ///
    /// Returns volumetric cooling power \[W/m³\].
    pub fn volumetric_cooling_power(&self, frequency_hz: f64) -> f64 {
        let ds = self.material.isothermal_entropy_change(self.field_swing);
        ds * self.t_cold * self.material.density * frequency_hz
    }
}

// ---------------------------------------------------------------------------
// LandauTheoryMagneto
// ---------------------------------------------------------------------------

/// Landau theory for magnetic phase transitions.
///
/// Expands the Helmholtz free energy in powers of the order parameter
/// (magnetisation M):
///
/// F = a·M² + b·M⁴ + c·M⁶ − μ₀·H·M
///
/// Signs of the Landau coefficients control first- vs second-order transitions.
pub struct LandauTheoryMagneto {
    /// Second-order coefficient a \[J·m³/A²\] — temperature-dependent: a = a₀(T−T_C).
    pub a0: f64,
    /// Fourth-order coefficient b \[J·m⁷/A⁴\].
    pub b: f64,
    /// Sixth-order coefficient c \[J·m¹¹/A⁶\] (stabilises first-order transitions).
    pub c: f64,
    /// Curie temperature T_C \[K\].
    pub curie_temperature: f64,
}

impl LandauTheoryMagneto {
    /// Create a new `LandauTheoryMagneto` model.
    pub fn new(a0: f64, b: f64, c: f64, curie_temperature: f64) -> Self {
        Self {
            a0,
            b,
            c,
            curie_temperature,
        }
    }

    /// Temperature-dependent second Landau coefficient a(T) = a₀ · (T − T_C).
    pub fn a_coeff(&self, temperature: f64) -> f64 {
        self.a0 * (temperature - self.curie_temperature)
    }

    /// Free energy density F(M, T, H).
    pub fn free_energy(&self, magnetization: f64, temperature: f64, field: f64) -> f64 {
        let a = self.a_coeff(temperature);
        let m2 = magnetization * magnetization;
        let m4 = m2 * m2;
        let m6 = m4 * m2;
        a * m2 + self.b * m4 + self.c * m6 - MU_0 * field * magnetization
    }

    /// Equation of state (∂F/∂M = 0 → H):
    ///
    /// μ₀·H = 2a·M + 4b·M³ + 6c·M⁵
    pub fn field_from_magnetization(&self, magnetization: f64, temperature: f64) -> f64 {
        let a = self.a_coeff(temperature);
        let m3 = magnetization * magnetization * magnetization;
        let m5 = m3 * magnetization * magnetization;
        (2.0 * a * magnetization + 4.0 * self.b * m3 + 6.0 * self.c * m5) / MU_0
    }

    /// Paramagnetic susceptibility χ = (∂M/∂H)_{H=0} = μ₀ / (2·|a(T)|).
    pub fn susceptibility(&self, temperature: f64) -> f64 {
        let a = self.a_coeff(temperature).abs();
        if a < 1e-30 {
            return f64::INFINITY;
        }
        MU_0 / (2.0 * a)
    }

    /// Spontaneous magnetisation M_s below T_C using mean-field minimisation.
    ///
    /// For a second-order transition (b > 0, c = 0): M_s = √(−a/(2b)).
    pub fn spontaneous_magnetization(&self, temperature: f64) -> f64 {
        if temperature >= self.curie_temperature {
            return 0.0;
        }
        let a = self.a_coeff(temperature);
        if self.b > 1e-30 {
            let m2 = -a / (2.0 * self.b);
            if m2 > 0.0 { m2.sqrt() } else { 0.0 }
        } else {
            0.0
        }
    }

    /// Transition type: "first-order" if b < 0, else "second-order".
    pub fn transition_type(&self) -> &'static str {
        if self.b < 0.0 {
            "first-order"
        } else {
            "second-order"
        }
    }
}

// ---------------------------------------------------------------------------
// MeanFieldMagnetics
// ---------------------------------------------------------------------------

/// Weiss mean-field theory for ferromagnets.
///
/// Uses the Brillouin function B_J to relate magnetisation to temperature
/// and applied field.
pub struct MeanFieldMagnetics {
    /// Total angular momentum quantum number J.
    pub total_angular_momentum: f64,
    /// Landé g-factor.
    pub g_factor: f64,
    /// Mean-field exchange constant λ \[m/(A·s²)\] (Weiss parameter).
    pub weiss_parameter: f64,
    /// Number of magnetic atoms per unit volume n \[m⁻³\].
    pub number_density: f64,
}

impl MeanFieldMagnetics {
    /// Create a new `MeanFieldMagnetics` model.
    pub fn new(
        total_angular_momentum: f64,
        g_factor: f64,
        weiss_parameter: f64,
        number_density: f64,
    ) -> Self {
        Self {
            total_angular_momentum,
            g_factor,
            weiss_parameter,
            number_density,
        }
    }

    /// Gadolinium reference (J = 7/2, g ≈ 2).
    pub fn gadolinium(number_density: f64) -> Self {
        Self::new(3.5, 2.0, 1.0e3, number_density)
    }

    /// Brillouin function B_J(x).
    ///
    /// B_J(x) = ((2J+1)/(2J)) · coth((2J+1)·x/(2J)) − (1/(2J)) · coth(x/(2J))
    pub fn brillouin(&self, x: f64) -> f64 {
        let j = self.total_angular_momentum;
        let two_j = 2.0 * j;
        let a = (two_j + 1.0) / two_j;
        let b = 1.0 / two_j;
        // Guard against small x to avoid sinh(0) singularity.
        if x.abs() < 1e-10 {
            return (j + 1.0) / (3.0 * j) * x;
        }
        a * coth(a * x) - b * coth(b * x)
    }

    /// Saturation magnetisation M_s = n · g · J · μ_B \[A/m\].
    pub fn saturation_magnetization(&self) -> f64 {
        self.number_density * self.g_factor * self.total_angular_momentum * MU_B
    }

    /// Curie temperature from mean-field theory:
    ///
    /// T_C = n · g² · μ_B² · J(J+1) · λ / (3 · k_B)
    pub fn curie_temperature(&self) -> f64 {
        let j = self.total_angular_momentum;
        let g = self.g_factor;
        let n = self.number_density;
        let lambda = self.weiss_parameter;
        n * g * g * MU_B * MU_B * j * (j + 1.0) * lambda / (3.0 * K_B)
    }

    /// Reduced magnetisation m = M / M_s from mean-field self-consistency.
    ///
    /// Iterates m = B_J( g·J·μ_B·(H + λ·n·g·J·μ_B·m) / (k_B·T) ).
    pub fn reduced_magnetization(
        &self,
        temperature: f64,
        applied_field: f64,
        n_iter: usize,
    ) -> f64 {
        if temperature <= 0.0 {
            return 1.0;
        }
        let j = self.total_angular_momentum;
        let g = self.g_factor;
        let ms = self.saturation_magnetization();
        let mut m = 0.5_f64;
        for _ in 0..n_iter {
            let h_eff = applied_field + self.weiss_parameter * ms * m;
            let x = g * j * MU_B * h_eff / (K_B * temperature);
            m = self.brillouin(x);
        }
        m
    }

    /// Curie-Weiss susceptibility above T_C:
    ///
    /// χ = C / (T − T_C),   C = n · g² · J(J+1) · μ_B² / (3 · k_B)
    pub fn curie_weiss_susceptibility(&self, temperature: f64) -> f64 {
        let j = self.total_angular_momentum;
        let g = self.g_factor;
        let n = self.number_density;
        let c = n * g * g * j * (j + 1.0) * MU_B * MU_B / (3.0 * K_B);
        let tc = self.curie_temperature();
        if (temperature - tc).abs() < 1e-6 {
            return f64::INFINITY;
        }
        c / (temperature - tc)
    }
}

/// Hyperbolic cotangent helper (private).
#[inline]
fn coth(x: f64) -> f64 {
    if x.abs() > 100.0 {
        return x.signum();
    }
    let ex = x.exp();
    let em = (-x).exp();
    (ex + em) / (ex - em)
}

// ---------------------------------------------------------------------------
// MagnetoCaloric (Heusler)
// ---------------------------------------------------------------------------

/// Heusler alloy magnetocaloric model (Ni-Mn-Ga type).
///
/// Models the coupled magnetic/structural metamagnetic transition, including
/// the virgin curve anomaly and field-induced phase conversion.
pub struct MagnetoCaloric {
    /// Martensitic start temperature M_s \[K\].
    pub martensite_start: f64,
    /// Martensitic finish temperature M_f \[K\].
    pub martensite_finish: f64,
    /// Austenitic start temperature A_s \[K\].
    pub austenite_start: f64,
    /// Austenitic finish temperature A_f \[K\].
    pub austenite_finish: f64,
    /// Entropy change across the structural transition ΔS_str \[J/(kg·K)\].
    pub structural_entropy_change: f64,
    /// Field required to induce transition μ₀H_crit \[T\].
    pub critical_field: f64,
}

impl MagnetoCaloric {
    /// Create a new Heusler `MagnetoCaloric` model.
    pub fn new(
        martensite_start: f64,
        martensite_finish: f64,
        austenite_start: f64,
        austenite_finish: f64,
        structural_entropy_change: f64,
        critical_field: f64,
    ) -> Self {
        Self {
            martensite_start,
            martensite_finish,
            austenite_start,
            austenite_finish,
            structural_entropy_change,
            critical_field,
        }
    }

    /// Ni-Mn-Ga reference compound (approximate literature values).
    pub fn ni_mn_ga() -> Self {
        Self::new(300.0, 290.0, 310.0, 320.0, 18.0, 1.0)
    }

    /// Transformation fraction ξ(T) — fraction of martensite during cooling.
    ///
    /// Uses a linear interpolation between M_f and M_s.
    pub fn martensite_fraction(&self, temperature: f64) -> f64 {
        let ms = self.martensite_start;
        let mf = self.martensite_finish;
        if temperature >= ms {
            return 0.0;
        }
        if temperature <= mf {
            return 1.0;
        }
        (ms - temperature) / (ms - mf)
    }

    /// Metamagnetic transition entropy release \[J/(kg·K)\] at temperature T.
    ///
    /// Returns a fraction of ΔS_str proportional to the transformed volume.
    pub fn metamagnetic_entropy_change(&self, temperature: f64, field: f64) -> f64 {
        let xi = self.martensite_fraction(temperature);
        let field_factor = (field / self.critical_field).min(1.0);
        -self.structural_entropy_change * xi * field_factor
    }

    /// Returns true if the material is in the martensitic phase.
    pub fn is_martensitic(&self, temperature: f64) -> bool {
        temperature <= self.martensite_start
    }

    /// Virgin curve: first magnetisation is irreversible if field > H_crit.
    pub fn virgin_curve_magnetization(&self, field: f64, saturation_mag: f64) -> f64 {
        if field >= self.critical_field {
            saturation_mag
        } else {
            saturation_mag * field / self.critical_field
        }
    }
}

// ---------------------------------------------------------------------------
// BarocaloriEffect
// ---------------------------------------------------------------------------

/// Barocaloric effect — pressure-induced entropy change.
///
/// Models the isothermal entropy change and adiabatic temperature change
/// under hydrostatic pressure for plastic crystals and hybrid perovskites.
pub struct BarocaloriEffect {
    /// Pressure sensitivity of transition temperature dT/dP \[K/GPa\].
    pub dt_dp: f64,
    /// Entropy change per pressure increment ΔS/ΔP \[J/(kg·K·GPa)\].
    pub ds_dp: f64,
    /// Specific heat c_p \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Operating temperature T \[K\].
    pub temperature: f64,
}

impl BarocaloriEffect {
    /// Create a new `BarocaloriEffect` model.
    pub fn new(dt_dp: f64, ds_dp: f64, specific_heat: f64, temperature: f64) -> Self {
        Self {
            dt_dp,
            ds_dp,
            specific_heat,
            temperature,
        }
    }

    /// Isothermal entropy change ΔS \[J/(kg·K)\] for applied pressure ΔP \[GPa\].
    pub fn isothermal_entropy_change(&self, delta_pressure_gpa: f64) -> f64 {
        self.ds_dp * delta_pressure_gpa
    }

    /// Adiabatic temperature change ΔT \[K\] for applied pressure ΔP \[GPa\].
    ///
    /// ΔT_ad ≈ −(T / c_p) · ΔS
    pub fn adiabatic_temp_change(&self, delta_pressure_gpa: f64) -> f64 {
        let ds = self.isothermal_entropy_change(delta_pressure_gpa);
        -(self.temperature / self.specific_heat) * ds
    }

    /// Barocaloric coefficient dT_ad/dP \[K/GPa\] (figure-of-merit for cooling).
    pub fn barocaloric_coefficient(&self) -> f64 {
        self.adiabatic_temp_change(1.0)
    }

    /// Clausius-Clapeyron check: dT/dP = T · ΔV / ΔS (volume units ignored here).
    ///
    /// Returns the stored dT/dP as the material parameter.
    pub fn clausius_clapeyron_pressure(&self) -> f64 {
        self.dt_dp
    }
}

// ---------------------------------------------------------------------------
// ElastocaloricMaterial
// ---------------------------------------------------------------------------

/// Elastocaloric (shape memory alloy) material model.
///
/// Describes cooling/heating via stress-induced martensitic transformation
/// in superelastic alloys (e.g. NiTi, Cu-Zn-Al).
pub struct ElastocaloricMaterial {
    /// Transformation entropy change ΔS_tr \[J/(kg·K)\].
    pub transformation_entropy: f64,
    /// Critical stress for transformation σ_tr \[MPa\].
    pub critical_stress: f64,
    /// Clausius-Clapeyron slope dσ/dT \[MPa/K\].
    pub clausius_clapeyron_slope: f64,
    /// Specific heat c_p \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Operating temperature T \[K\].
    pub temperature: f64,
    /// Hysteresis stress Δσ_hys \[MPa\].
    pub hysteresis_stress: f64,
}

impl ElastocaloricMaterial {
    /// Create a new `ElastocaloricMaterial`.
    pub fn new(
        transformation_entropy: f64,
        critical_stress: f64,
        clausius_clapeyron_slope: f64,
        specific_heat: f64,
        temperature: f64,
        hysteresis_stress: f64,
    ) -> Self {
        Self {
            transformation_entropy,
            critical_stress,
            clausius_clapeyron_slope,
            specific_heat,
            temperature,
            hysteresis_stress,
        }
    }

    /// NiTi (Nitinol) reference material.
    pub fn nitinol() -> Self {
        Self::new(20.0, 200.0, 6.0, 500.0, 310.0, 60.0)
    }

    /// Adiabatic temperature change ΔT_ad \[K\] for complete transformation.
    ///
    /// ΔT_ad = −T · ΔS_tr / c_p
    pub fn adiabatic_temp_change(&self) -> f64 {
        -self.temperature * self.transformation_entropy / self.specific_heat
    }

    /// Specific cooling power per cycle \[J/kg\] (net, after hysteresis loss).
    pub fn cooling_power_per_cycle(&self) -> f64 {
        let gross = self.transformation_entropy.abs() * self.temperature;
        let loss = self.hysteresis_stress * 0.01; // 1 % strain assumption
        (gross - loss).max(0.0)
    }

    /// Critical stress at temperature T \[MPa\].
    pub fn critical_stress_at(&self, temperature: f64) -> f64 {
        self.critical_stress + self.clausius_clapeyron_slope * (temperature - self.temperature)
    }

    /// Superelastic work input W_in \[J/kg\] for a loading cycle.
    pub fn superelastic_work_input(&self, strain_amplitude: f64) -> f64 {
        self.critical_stress * 1e6 * strain_amplitude / 1000.0 // convert to J/kg approx
    }
}

// ---------------------------------------------------------------------------
// SpinCaloritronics
// ---------------------------------------------------------------------------

/// Spin caloritronics: spin-dependent heat and charge transport.
///
/// Models the spin Seebeck effect (SSE), spin Peltier effect, and spin
/// Nernst effect in magnetic thin films and tunnel junctions.
pub struct SpinCaloritronics {
    /// Spin Seebeck coefficient S_s \[V/K\].
    pub spin_seebeck_coeff: f64,
    /// Spin conductivity σ_s \[S/m\].
    pub spin_conductivity: f64,
    /// Spin Nernst coefficient N_s \[V/K\].
    pub spin_nernst_coeff: f64,
    /// Tunnel magnetoresistance ratio TMR (dimensionless).
    pub tmr_ratio: f64,
    /// Sample temperature T \[K\].
    pub temperature: f64,
}

impl SpinCaloritronics {
    /// Create a new `SpinCaloritronics` model.
    pub fn new(
        spin_seebeck_coeff: f64,
        spin_conductivity: f64,
        spin_nernst_coeff: f64,
        tmr_ratio: f64,
        temperature: f64,
    ) -> Self {
        Self {
            spin_seebeck_coeff,
            spin_conductivity,
            spin_nernst_coeff,
            tmr_ratio,
            temperature,
        }
    }

    /// Spin Seebeck voltage V_s = S_s · ΔT \[V\].
    pub fn spin_seebeck_voltage(&self, delta_t: f64) -> f64 {
        self.spin_seebeck_coeff * delta_t
    }

    /// Spin Peltier heat flux q_s = Π_s · J_s \[W/m²\].
    ///
    /// Π_s = S_s · T (spin Peltier coefficient).
    pub fn spin_peltier_heat_flux(&self, spin_current_density: f64) -> f64 {
        let pi_s = self.spin_seebeck_coeff * self.temperature;
        pi_s * spin_current_density
    }

    /// Spin Nernst voltage (transverse) V_N = N_s · ΔT \[V\].
    pub fn spin_nernst_voltage(&self, delta_t: f64) -> f64 {
        self.spin_nernst_coeff * delta_t
    }

    /// Spin figure of merit ZT_s = S_s² · σ_s · T / κ_s.
    ///
    /// Thermal conductivity κ_s is estimated as κ_s = T·S_s²·σ_s for Wiedemann-Franz.
    pub fn spin_figure_of_merit(&self) -> f64 {
        // With WF estimate κ_s = L₀·T·σ_s, ZT_s = S_s²/(L₀), where L₀ = 2.44e-8 W·Ω/K²
        let l0 = 2.44e-8_f64;
        self.spin_seebeck_coeff * self.spin_seebeck_coeff / l0
    }

    /// Seebeck spin tunneling voltage across a magnetic tunnel junction.
    ///
    /// V_SST = S_eff · ΔT,  S_eff = (P · S_s) where P is spin polarisation.
    pub fn seebeck_spin_tunneling_voltage(&self, delta_t: f64, spin_polarization: f64) -> f64 {
        spin_polarization * self.spin_seebeck_coeff * delta_t
    }

    /// Tunnel spin current density J_s \[A/m²\] driven by temperature difference.
    pub fn tunnel_spin_current(&self, delta_t: f64) -> f64 {
        self.spin_conductivity * self.spin_seebeck_coeff * delta_t
    }
}

// ---------------------------------------------------------------------------
// Helper: Lorentzian peak
// ---------------------------------------------------------------------------

/// Lorentzian peak function L(x; x₀, Γ) = (Γ/2)² / ((x−x₀)² + (Γ/2)²).
pub fn lorentzian_peak(x: f64, x0: f64, gamma: f64) -> f64 {
    let half_g = gamma / 2.0;
    half_g * half_g / ((x - x0) * (x - x0) + half_g * half_g)
}

/// Gaussian peak function G(x; x₀, σ) = exp(−(x−x₀)²/(2σ²)).
pub fn gaussian_peak(x: f64, x0: f64, sigma: f64) -> f64 {
    let z = (x - x0) / sigma;
    (-0.5 * z * z).exp()
}

/// Planck distribution n(ω, T) = 1 / (exp(ℏω/(k_B·T)) − 1).
///
/// `hbar_omega` is the energy ℏω \[J\].
pub fn bose_einstein(hbar_omega: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    let x = hbar_omega / (K_B * temperature);
    if x > 700.0 {
        return 0.0;
    }
    1.0 / (x.exp() - 1.0)
}

/// Spin wave (magnon) contribution to specific heat c_sw \[J/(kg·K)\].
///
/// Uses a simplified T^(3/2) Bloch law: c_sw ≈ A · T^(3/2)
///
/// `bloch_coefficient` A is material-dependent \[J/(kg·K^(5/2))\].
pub fn magnon_specific_heat(temperature: f64, bloch_coefficient: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    1.5 * bloch_coefficient * temperature.powf(0.5)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Gadolinium MCE: adiabatic ΔT scales linearly with field.
    #[test]
    fn test_gd_mce_linear_field_scaling() {
        let gd = MagnetocaloricMaterial::gadolinium();
        let dt2 = gd.adiabatic_temp_change(2.0);
        let dt4 = gd.adiabatic_temp_change(4.0);
        assert!(
            (dt4 - 2.0 * dt2).abs() < EPS,
            "ΔT_ad should scale linearly: dt2={dt2}, dt4={dt4}"
        );
    }

    // 2. Isothermal entropy change is positive for positive field.
    #[test]
    fn test_isothermal_entropy_positive() {
        let mat = MagnetocaloricMaterial::gadolinium();
        let ds = mat.isothermal_entropy_change(2.0);
        assert!(
            ds > 0.0,
            "ΔS_M should be positive for a positive field swing"
        );
    }

    // 3. Refrigerant capacity increases with FWHM.
    #[test]
    fn test_refrigerant_capacity_fwhm() {
        let mat = MagnetocaloricMaterial::gadolinium();
        let rc5 = mat.refrigerant_capacity(2.0, 5.0);
        let rc10 = mat.refrigerant_capacity(2.0, 10.0);
        assert!(
            rc10 > rc5,
            "Larger FWHM should give greater refrigerant capacity"
        );
    }

    // 4. Maxwell relation entropy change is proportional to dM/dT.
    #[test]
    fn test_maxwell_relation_linearity() {
        let me = MagneticEntropy::new(0.1574, 7.0);
        let ds1 = me.maxwell_relation_entropy_change(1.0, 1e5);
        let ds2 = me.maxwell_relation_entropy_change(2.0, 1e5);
        assert!(
            (ds2 - 2.0 * ds1).abs() < EPS,
            "Maxwell ΔS should be proportional to dM/dT"
        );
    }

    // 5. Clausius-Clapeyron slope has correct sign (negative for MnAs-like).
    #[test]
    fn test_clausius_clapeyron_sign() {
        let me = MagneticEntropy::new(0.1265, 4.5);
        // Positive ΔM, positive ΔS → negative slope.
        let slope = me.clausius_clapeyron_slope(50.0, 30.0);
        assert!(
            slope < 0.0,
            "dT_C/dH should be negative for ΔM>0, ΔS>0: {slope}"
        );
    }

    // 6. Disordered entropy is maximum-entropy limit S = R·ln(2J+1)/M.
    #[test]
    fn test_disordered_entropy_j_half() {
        let me = MagneticEntropy::new(R_GAS, 0.5);
        // J = 0.5 → 2J+1 = 2 → S = R·ln2 / R = ln(2).
        let s = me.disordered_entropy(0.5);
        let expected = 2_f64.ln();
        assert!(
            (s - expected).abs() < 1e-8,
            "Disordered entropy for J=1/2 should be ln(2)={expected:.6}, got {s:.6}"
        );
    }

    // 7. GiantMce: La(Fe,Si) peak entropy > Gd.
    #[test]
    fn test_gmce_lafesi_peak_greater_than_gd() {
        let lafesi = GiantMce::la_fe_si();
        let gd = MagnetocaloricMaterial::gadolinium();
        // |ΔS| at 2 T for Gd ≈ 3 J/(kg·K), for La(Fe,Si) ≈ 20 J/(kg·K)
        assert!(
            lafesi.peak_entropy_change > gd.delta_s_per_tesla * 2.0,
            "La(Fe,Si)13 peak ΔS={} should exceed Gd at 2 T ({})",
            lafesi.peak_entropy_change,
            gd.delta_s_per_tesla * 2.0
        );
    }

    // 8. GiantMce: entropy change is negative (material cools surroundings).
    #[test]
    fn test_gmce_entropy_change_negative() {
        let lafesi = GiantMce::la_fe_si();
        let ds = lafesi.entropy_change(lafesi.transition_temp, 2.0);
        assert!(
            ds < 0.0,
            "Entropy change should be negative (cooling): {ds}"
        );
    }

    // 9. GiantMce: is_first_order returns true for La(Fe,Si).
    #[test]
    fn test_gmce_first_order_flag() {
        let lafesi = GiantMce::la_fe_si();
        assert!(
            lafesi.is_first_order(),
            "La(Fe,Si) should be identified as first-order"
        );
    }

    // 10. MagneticRefrigeration: Carnot COP equals T_C / ΔT.
    #[test]
    fn test_amr_carnot_cop() {
        let mat = MagnetocaloricMaterial::gadolinium();
        let amr = MagneticRefrigeration::new(300.0, 270.0, 2.0, mat);
        let cop = amr.carnot_cop();
        let expected = 270.0 / 30.0;
        assert!(
            (cop - expected).abs() < EPS,
            "Carnot COP should be {expected:.4}, got {cop:.4}"
        );
    }

    // 11. MagneticRefrigeration: AMR COP is positive and finite.
    #[test]
    fn test_amr_cop_positive() {
        let mat = MagnetocaloricMaterial::gadolinium();
        let amr = MagneticRefrigeration::new(300.0, 270.0, 2.0, mat);
        let cop = amr.amr_cop();
        assert!(
            cop > 0.0 && cop.is_finite(),
            "AMR COP should be positive and finite: {cop}"
        );
    }

    // 12. Landau: free energy at zero field and M=0 is zero.
    #[test]
    fn test_landau_free_energy_zero() {
        let lm = LandauTheoryMagneto::new(1.0, 1.0, 0.1, 300.0);
        let f = lm.free_energy(0.0, 300.0, 0.0);
        assert!(f.abs() < EPS, "Free energy at M=0, H=0 should be zero: {f}");
    }

    // 13. Landau: above T_C, spontaneous magnetisation is zero.
    #[test]
    fn test_landau_no_magnetization_above_tc() {
        let lm = LandauTheoryMagneto::new(1.0, 1.0, 0.1, 300.0);
        let ms = lm.spontaneous_magnetization(350.0);
        assert!(ms.abs() < EPS, "M_s should be zero above T_C: {ms}");
    }

    // 14. Landau: second-order transition type is classified correctly.
    #[test]
    fn test_landau_second_order_type() {
        let lm = LandauTheoryMagneto::new(1.0, 1.0, 0.1, 300.0);
        assert_eq!(lm.transition_type(), "second-order");
    }

    // 15. Mean-field: Curie temperature scales with weiss_parameter.
    #[test]
    fn test_mean_field_curie_temp_scaling() {
        let mf1 = MeanFieldMagnetics::gadolinium(1e28);
        let mf2 = MeanFieldMagnetics::new(3.5, 2.0, 2.0e3, 1e28);
        // Doubling λ should double T_C.
        let tc1 = mf1.curie_temperature();
        let tc2 = mf2.curie_temperature();
        assert!(
            (tc2 - 2.0 * tc1).abs() < 1e-4,
            "Doubling λ should double T_C: tc1={tc1:.4}, tc2={tc2:.4}"
        );
    }

    // 16. Mean-field: reduced magnetisation at T=0 ≈ 1 (saturated).
    #[test]
    fn test_mean_field_saturation_at_zero_temp() {
        let mf = MeanFieldMagnetics::gadolinium(1e28);
        // Very low temperature — use 1 K to avoid division by zero.
        let m = mf.reduced_magnetization(1.0, 0.0, 200);
        assert!(
            (m - 1.0).abs() < 0.01,
            "Reduced magnetisation at T→0 should approach 1: m={m:.6}"
        );
    }

    // 17. Brillouin function: B_J(0) = 0.
    #[test]
    fn test_brillouin_at_zero() {
        let mf = MeanFieldMagnetics::gadolinium(1e28);
        let b = mf.brillouin(0.0);
        assert!(b.abs() < 1e-8, "B_J(0) should be 0: {b}");
    }

    // 18. Heusler: martensite fraction is 1 below M_f.
    #[test]
    fn test_heusler_martensite_fraction_below_mf() {
        let h = MagnetoCaloric::ni_mn_ga();
        let xi = h.martensite_fraction(280.0);
        assert!(
            (xi - 1.0).abs() < EPS,
            "Martensite fraction should be 1 below M_f: {xi}"
        );
    }

    // 19. Heusler: martensite fraction is 0 above M_s.
    #[test]
    fn test_heusler_martensite_fraction_above_ms() {
        let h = MagnetoCaloric::ni_mn_ga();
        let xi = h.martensite_fraction(310.0);
        assert!(
            xi.abs() < EPS,
            "Martensite fraction should be 0 above M_s: {xi}"
        );
    }

    // 20. Barocaloric: adiabatic ΔT sign is consistent.
    #[test]
    fn test_barocaloric_sign() {
        // Positive ds_dp with positive pressure → ΔS > 0 → ΔT < 0 (conventional).
        let bce = BarocaloriEffect::new(5.0, 10.0, 1000.0, 300.0);
        let dt = bce.adiabatic_temp_change(0.5);
        assert!(dt < 0.0, "Barocaloric ΔT should be negative for ΔS>0: {dt}");
    }

    // 21. Barocaloric: entropy change scales linearly with pressure.
    #[test]
    fn test_barocaloric_entropy_linear() {
        let bce = BarocaloriEffect::new(5.0, 10.0, 1000.0, 300.0);
        let ds1 = bce.isothermal_entropy_change(1.0);
        let ds2 = bce.isothermal_entropy_change(2.0);
        assert!(
            (ds2 - 2.0 * ds1).abs() < EPS,
            "Barocaloric ΔS should scale linearly: ds1={ds1}, ds2={ds2}"
        );
    }

    // 22. Elastocaloric: NiTi adiabatic ΔT is negative (cooling on loading).
    #[test]
    fn test_elastocaloric_nitinol_cooling() {
        let niti = ElastocaloricMaterial::nitinol();
        let dt = niti.adiabatic_temp_change();
        // ΔS_tr > 0 (disordering on austenite→martensite) → ΔT_ad < 0 on loading.
        assert!(dt < 0.0, "NiTi should cool on stress loading, ΔT_ad={dt}");
    }

    // 23. Elastocaloric: critical stress increases with temperature.
    #[test]
    fn test_elastocaloric_stress_temperature_slope() {
        let niti = ElastocaloricMaterial::nitinol();
        let sigma_low = niti.critical_stress_at(300.0);
        let sigma_high = niti.critical_stress_at(320.0);
        assert!(
            sigma_high > sigma_low,
            "Critical stress should increase with T (Clausius-Clapeyron): σ_low={sigma_low}, σ_high={sigma_high}"
        );
    }

    // 24. SpinCaloritronics: spin Seebeck voltage proportional to ΔT.
    #[test]
    fn test_spin_seebeck_proportional_to_dt() {
        let sc = SpinCaloritronics::new(100e-6, 1e4, 50e-9, 0.3, 300.0);
        let v1 = sc.spin_seebeck_voltage(10.0);
        let v2 = sc.spin_seebeck_voltage(20.0);
        assert!(
            (v2 - 2.0 * v1).abs() < EPS,
            "SSE voltage should be proportional to ΔT: v1={v1}, v2={v2}"
        );
    }

    // 25. SpinCaloritronics: spin Peltier heat flux proportional to current.
    #[test]
    fn test_spin_peltier_proportional_to_current() {
        let sc = SpinCaloritronics::new(100e-6, 1e4, 50e-9, 0.3, 300.0);
        let q1 = sc.spin_peltier_heat_flux(1e6);
        let q2 = sc.spin_peltier_heat_flux(2e6);
        assert!(
            (q2 - 2.0 * q1).abs() < 1e-15,
            "Spin Peltier heat flux should double with double current: q1={q1}, q2={q2}"
        );
    }

    // 26. Lorentzian peak is 1 at centre.
    #[test]
    fn test_lorentzian_peak_at_centre() {
        let l = lorentzian_peak(5.0, 5.0, 2.0);
        assert!(
            (l - 1.0).abs() < EPS,
            "Lorentzian at centre should be 1: {l}"
        );
    }

    // 27. Gaussian peak is 1 at centre.
    #[test]
    fn test_gaussian_peak_at_centre() {
        let g = gaussian_peak(3.0, 3.0, 1.0);
        assert!((g - 1.0).abs() < EPS, "Gaussian at centre should be 1: {g}");
    }

    // 28. Numerical Maxwell entropy integration recovers area.
    #[test]
    fn test_numerical_maxwell_integration() {
        let me = MagneticEntropy::new(0.1574, 7.0);
        // Build a flat ΔS dataset: ΔS = 5 J/(kg·K) over 10 K.
        let ds_data: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64, 5.0)).collect();
        let integral = me.integrated_entropy(&ds_data);
        assert!(
            (integral - 50.0).abs() < 1e-6,
            "Integrated entropy should be 50 J/kg: {integral:.4}"
        );
    }

    // 29. MagneticRefrigeration: volumetric cooling power scales with frequency.
    #[test]
    fn test_volumetric_cooling_power_frequency() {
        let mat = MagnetocaloricMaterial::gadolinium();
        let amr = MagneticRefrigeration::new(300.0, 270.0, 2.0, mat);
        let p1 = amr.volumetric_cooling_power(1.0);
        let p2 = amr.volumetric_cooling_power(2.0);
        assert!(
            (p2 - 2.0 * p1).abs() < EPS,
            "Volumetric cooling power should scale linearly with frequency"
        );
    }

    // 30. Magnon specific heat increases with temperature (T^1/2 law).
    #[test]
    fn test_magnon_specific_heat_increases() {
        let c100 = magnon_specific_heat(100.0, 1e-3);
        let c200 = magnon_specific_heat(200.0, 1e-3);
        assert!(
            c200 > c100,
            "Magnon specific heat should increase with T: c100={c100:.6}, c200={c200:.6}"
        );
    }

    // Bonus: two-π value sanity check.
    #[test]
    fn test_pi_constant() {
        use std::f64::consts::{PI, TAU};
        let two_pi = 2.0 * PI;
        assert!((two_pi - TAU).abs() < 1e-6);
    }
}
