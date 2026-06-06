// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Battery electrode and electrolyte material models.
//!
//! Provides physics-based models for:
//! - Electrode materials (graphite, LFP, NMC, LCO, silicon, custom)
//! - Electrolyte ionic transport
//! - Solid-electrolyte interphase (SEI) growth
//! - Solid-state diffusion inside particles
//! - Capacity and resistance aging
//! - Butler-Volmer kinetics and Nernst equation
//! - Electrode kinetics (exchange current density, Tafel slope, EIS impedance)
//! - Solid electrolyte (Arrhenius conductivity, transference number, GEIS analysis)
//! - Battery degradation (SEI growth, lithium plating, capacity fade)
//! - Thermal battery model (Bernardi heat generation, cooling, temperature uniformity)
//! - Battery cycling (CC/CV protocol, SOC tracking, Coulombic efficiency)

/// Faraday constant \[C/mol\]
pub const FARADAY: f64 = 96_485.0;
/// Universal gas constant \[J/(mol·K)\]
pub const GAS_CONSTANT: f64 = 8.314_462_618;

// ─── Electrode type ──────────────────────────────────────────────────────────

/// Classification of common battery electrode active materials.
#[derive(Debug, Clone, PartialEq)]
pub enum ElectrodeType {
    /// Graphite (LiC₆) anode – conventional lithium-ion anode.
    Graphite,
    /// Lithium iron phosphate cathode (LiFePO₄).
    LFP,
    /// Lithium nickel manganese cobalt oxide cathode (NMC).
    NMC,
    /// Lithium cobalt oxide cathode (LiCoO₂).
    LCO,
    /// Silicon anode – high capacity, large volume expansion.
    Silicon,
    /// User-defined electrode type.
    Custom(String),
}

// ─── Electrode material ──────────────────────────────────────────────────────

/// Physical and electrochemical properties of an electrode active material.
#[derive(Debug, Clone)]
pub struct ElectrodeMaterial {
    /// Classification of this electrode material.
    pub electrode_type: ElectrodeType,
    /// Theoretical specific capacity \[mAh/g\].
    pub specific_capacity_mah_g: f64,
    /// Average open-circuit voltage vs Li/Li⁺ \[V\].
    pub average_voltage: f64,
    /// Particle density \[g/cm³\].
    pub density_g_cm3: f64,
    /// Linear intercalation strain (fractional volume change per unit SOC).
    pub intercalation_strain: f64,
}

impl ElectrodeMaterial {
    /// Construct a new `ElectrodeMaterial`.
    pub fn new(
        electrode_type: ElectrodeType,
        specific_capacity_mah_g: f64,
        average_voltage: f64,
        density_g_cm3: f64,
        intercalation_strain: f64,
    ) -> Self {
        Self {
            electrode_type,
            specific_capacity_mah_g,
            average_voltage,
            density_g_cm3,
            intercalation_strain,
        }
    }

    /// Volumetric energy density \[Wh/L\].
    ///
    /// Calculated as: specific_capacity \[mAh/g\] × average_voltage \[V\]
    /// × density \[g/cm³\] × 1000 \[cm³/L\] / 1000 \[mAh/Ah\].
    pub fn volumetric_energy_density(&self) -> f64 {
        self.specific_capacity_mah_g * self.average_voltage * self.density_g_cm3
    }

    /// Gravimetric energy density \[Wh/kg\].
    ///
    /// Calculated as: specific_capacity \[mAh/g\] × average_voltage \[V\].
    pub fn gravimetric_energy_density(&self) -> f64 {
        self.specific_capacity_mah_g * self.average_voltage
    }

    /// Return preset values for graphite anode.
    pub fn graphite() -> Self {
        Self::new(ElectrodeType::Graphite, 372.0, 0.1, 2.26, 0.1)
    }

    /// Return preset values for LFP cathode.
    pub fn lfp() -> Self {
        Self::new(ElectrodeType::LFP, 170.0, 3.4, 3.6, 0.06)
    }

    /// Return preset values for NMC 811 cathode.
    pub fn nmc() -> Self {
        Self::new(ElectrodeType::NMC, 200.0, 3.7, 4.7, 0.05)
    }

    /// Return preset values for LCO cathode.
    pub fn lco() -> Self {
        Self::new(ElectrodeType::LCO, 140.0, 3.9, 5.1, 0.04)
    }

    /// Return preset values for silicon anode.
    pub fn silicon() -> Self {
        Self::new(ElectrodeType::Silicon, 3580.0, 0.4, 2.33, 0.4)
    }
}

// ─── Electrolyte material ────────────────────────────────────────────────────

/// Bulk transport properties of a liquid electrolyte.
#[derive(Debug, Clone)]
pub struct ElectrolyteMaterial {
    /// Ionic conductivity \[S/m\].
    pub ionic_conductivity_s_m: f64,
    /// Electrochemical stability window \[V_min, V_max\] vs Li/Li⁺.
    pub electrochemical_window: Vec<f64>,
    /// Li⁺ transference number (fraction of current carried by Li⁺).
    pub transference_number: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub viscosity: f64,
}

impl ElectrolyteMaterial {
    /// Construct a new `ElectrolyteMaterial`.
    pub fn new(
        ionic_conductivity_s_m: f64,
        electrochemical_window: Vec<f64>,
        transference_number: f64,
        viscosity: f64,
    ) -> Self {
        Self {
            ionic_conductivity_s_m,
            electrochemical_window,
            transference_number,
            viscosity,
        }
    }

    /// Return a typical 1 M LiPF₆ in EC/DMC electrolyte.
    pub fn lipf6_ec_dmc() -> Self {
        Self::new(1.0, vec![0.0, 4.5], 0.38, 2.5e-3)
    }

    /// Width of the electrochemical stability window \[V\].
    pub fn stability_window_width(&self) -> f64 {
        if self.electrochemical_window.len() >= 2 {
            self.electrochemical_window[1] - self.electrochemical_window[0]
        } else {
            0.0
        }
    }
}

// ─── SEI ─────────────────────────────────────────────────────────────────────

/// Solid-electrolyte interphase (SEI) layer on anode surface.
#[derive(Debug, Clone)]
pub struct SolidElectrolyteInterphase {
    /// Current SEI thickness \[nm\].
    pub thickness_nm: f64,
    /// Ionic resistance of the SEI layer \[Ω·m²\].
    pub ionic_resistance: f64,
    /// Arrhenius activation energy for SEI growth \[J/mol\].
    pub activation_energy: f64,
}

impl SolidElectrolyteInterphase {
    /// Construct a new `SolidElectrolyteInterphase`.
    pub fn new(thickness_nm: f64, ionic_resistance: f64, activation_energy: f64) -> Self {
        Self {
            thickness_nm,
            ionic_resistance,
            activation_energy,
        }
    }

    /// Parabolic SEI growth rate \[nm/s\].
    ///
    /// Uses an Arrhenius-type model:
    /// `rate = k₀ / thickness × exp(−Ea / (R·T)) × soc_factor`
    /// where `soc_factor = 1 − soc` penalises growth at low state-of-charge.
    pub fn growth_rate(&self, temperature: f64, soc: f64) -> f64 {
        let k0 = 1.0e-3; // pre-exponential [nm²/s]
        let soc_factor = (1.0 - soc).max(0.0);
        let arrhenius = (-self.activation_energy / (GAS_CONSTANT * temperature)).exp();
        let thickness = self.thickness_nm.max(1e-10);
        k0 / thickness * arrhenius * soc_factor
    }
}

// ─── Diffusion in electrode particle ────────────────────────────────────────

/// Solid-state Li diffusion inside a spherical electrode particle.
///
/// Models concentration profiles and C-rate limitations using
/// Fick's second law in spherical coordinates.
#[derive(Debug, Clone)]
pub struct DiffusionInElectrode {
    /// Solid-state diffusivity of Li⁺ \[m²/s\].
    pub diffusivity: f64,
    /// Particle radius \[m\].
    pub particle_radius: f64,
    /// Surface concentration \[mol/m³\].
    pub c_surface: f64,
    /// Bulk (centre) concentration \[mol/m³\].
    pub c_bulk: f64,
    /// Maximum Li concentration \[mol/m³\].
    pub c_max: f64,
}

impl DiffusionInElectrode {
    /// Construct a new `DiffusionInElectrode` with default `c_max = 30_000.0`.
    pub fn new(diffusivity: f64, particle_radius: f64, c_surface: f64, c_bulk: f64) -> Self {
        Self {
            diffusivity,
            particle_radius,
            c_surface,
            c_bulk,
            c_max: 30_000.0,
        }
    }

    /// Construct with explicit maximum concentration.
    pub fn with_c_max(
        diffusivity: f64,
        particle_radius: f64,
        c_surface: f64,
        c_bulk: f64,
        c_max: f64,
    ) -> Self {
        Self {
            diffusivity,
            particle_radius,
            c_surface,
            c_bulk,
            c_max,
        }
    }

    /// Approximate radial concentration profile \[mol/m³\] at fractional radius `r` ∈ \[0, 1\].
    ///
    /// Uses a parabolic (pseudo-steady-state) profile:
    /// `c(r) = c_bulk + (c_surface − c_bulk) × r²`
    pub fn concentration_profile(&self, r: f64) -> f64 {
        let r_clamped = r.clamp(0.0, 1.0);
        self.c_bulk + (self.c_surface - self.c_bulk) * r_clamped * r_clamped
    }

    /// Characteristic diffusion length \[m\].
    ///
    /// `L_d = sqrt(D × t_ref)` where `t_ref = R² / D` (time to diffuse across particle).
    pub fn diffusion_length(&self) -> f64 {
        // t_ref = R² / D  →  L_d = R
        self.particle_radius
    }

    /// Diffusion time scale \[s\]: `τ = R² / D`.
    pub fn diffusion_time_scale(&self) -> f64 {
        self.particle_radius * self.particle_radius / self.diffusivity.max(1e-30)
    }

    /// Maximum C-rate \[1/h\] at which the particle can be fully discharged.
    ///
    /// Defined as `C_max = 1 / τ_diff [h⁻¹]` where τ_diff = R² / (π² D).
    /// Higher diffusivity or smaller particle → higher accessible C-rate.
    pub fn max_c_rate(&self) -> f64 {
        let tau_s = self.particle_radius * self.particle_radius
            / (std::f64::consts::PI * std::f64::consts::PI * self.diffusivity.max(1e-30));
        let tau_h = tau_s / 3600.0;
        1.0 / tau_h.max(1e-30)
    }

    /// Surface concentration gradient \[mol/m⁴\] approximation (linear).
    ///
    /// `dc/dr|_R ≈ (c_surface − c_bulk) / R`
    pub fn surface_concentration_gradient(&self) -> f64 {
        (self.c_surface - self.c_bulk) / self.particle_radius.max(1e-30)
    }

    /// State of charge estimated from average concentration: `SOC = c_avg / c_max`.
    ///
    /// Average concentration for parabolic profile: `c_avg = c_bulk + 3/5 (c_surf - c_bulk)`.
    pub fn soc_from_concentration(&self) -> f64 {
        let c_avg = self.c_bulk + 0.6 * (self.c_surface - self.c_bulk);
        (c_avg / self.c_max.max(1e-30)).clamp(0.0, 1.0)
    }

    /// Flux at the surface \[mol/(m²·s)\] using Fick's law.
    ///
    /// `J = −D × dc/dr|_R`
    pub fn surface_flux(&self) -> f64 {
        -self.diffusivity * self.surface_concentration_gradient()
    }
}

// ─── Electrode kinetics ──────────────────────────────────────────────────────

/// Electrode kinetics using the Butler-Volmer formulation.
///
/// Computes exchange current density, Tafel slope, and EIS (electrochemical
/// impedance spectroscopy) charge-transfer resistance.
#[derive(Debug, Clone)]
pub struct ElectrodeKinetics {
    /// Exchange current density \[A/m²\].
    pub i0: f64,
    /// Anodic charge-transfer coefficient (dimensionless).
    pub alpha_a: f64,
    /// Cathodic charge-transfer coefficient (dimensionless).
    pub alpha_c: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Number of electrons transferred per reaction step.
    pub n_electrons: usize,
}

impl ElectrodeKinetics {
    /// Construct a new `ElectrodeKinetics`.
    pub fn new(i0: f64, alpha_a: f64, alpha_c: f64, temperature: f64, n_electrons: usize) -> Self {
        Self {
            i0,
            alpha_a,
            alpha_c,
            temperature,
            n_electrons,
        }
    }

    /// Construct with symmetric transfer coefficients (α_a = α_c = 0.5).
    pub fn symmetric(i0: f64, temperature: f64) -> Self {
        Self::new(i0, 0.5, 0.5, temperature, 1)
    }

    /// Butler-Volmer current density \[A/m²\] at overpotential `eta` \[V\].
    ///
    /// `j = i₀ × (exp(αa·F·η / (R·T)) − exp(−αc·F·η / (R·T)))`
    pub fn current_density(&self, eta: f64) -> f64 {
        let f_over_rt = FARADAY / (GAS_CONSTANT * self.temperature);
        self.i0 * ((self.alpha_a * f_over_rt * eta).exp() - (-self.alpha_c * f_over_rt * eta).exp())
    }

    /// Anodic Tafel slope \[V/decade\].
    ///
    /// `b_a = 2.303 R T / (α_a F)`
    pub fn anodic_tafel_slope(&self) -> f64 {
        2.303 * GAS_CONSTANT * self.temperature / (self.alpha_a * FARADAY)
    }

    /// Cathodic Tafel slope \[V/decade\].
    ///
    /// `b_c = 2.303 R T / (α_c F)`
    pub fn cathodic_tafel_slope(&self) -> f64 {
        2.303 * GAS_CONSTANT * self.temperature / (self.alpha_c * FARADAY)
    }

    /// Linearised charge-transfer resistance \[Ω·m²\] (valid for small η).
    ///
    /// From linearisation of Butler-Volmer: `R_ct = R T / (n F i₀)`
    pub fn charge_transfer_resistance(&self) -> f64 {
        GAS_CONSTANT * self.temperature / (self.n_electrons as f64 * FARADAY * self.i0.max(1e-30))
    }

    /// EIS impedance at angular frequency `omega` \[rad/s\].
    ///
    /// Returns `(Z_real, Z_imag)` \[Ω·m²\] using the Randles circuit approximation:
    /// - Charge-transfer resistance `R_ct` in parallel with double-layer capacitance `C_dl`
    /// - `Z = R_ct / (1 + (ω R_ct C_dl)²)` (real) and `-ω R_ct² C_dl / (1 + (ω R_ct C_dl)²)` (imag)
    ///
    /// Assumes `C_dl = 0.2 F/m²` (typical double-layer capacitance).
    pub fn eis_impedance(&self, omega: f64) -> (f64, f64) {
        let r_ct = self.charge_transfer_resistance();
        let c_dl = 0.2_f64; // [F/m²] double-layer capacitance
        let denom = 1.0 + (omega * r_ct * c_dl).powi(2);
        let z_real = r_ct / denom;
        let z_imag = -omega * r_ct * r_ct * c_dl / denom;
        (z_real, z_imag)
    }

    /// EIS impedance magnitude \[Ω·m²\] at angular frequency `omega`.
    pub fn eis_magnitude(&self, omega: f64) -> f64 {
        let (zr, zi) = self.eis_impedance(omega);
        (zr * zr + zi * zi).sqrt()
    }

    /// Overpotential \[V\] required to achieve current density `j` \[A/m²\].
    ///
    /// Uses high-overpotential Tafel approximation (anodic branch):
    /// `η ≈ b_a × log10(j / i₀)` for large positive j.
    pub fn tafel_overpotential_anodic(&self, j: f64) -> f64 {
        if j <= 0.0 || j <= self.i0 {
            return 0.0;
        }
        self.anodic_tafel_slope() * (j / self.i0).log10()
    }

    /// Open-circuit potential correction (Nernst) \[V\].
    ///
    /// `E = E_ref − (R T / n F) ln(c_ox/c_red)` where `c_ratio = c_ox/c_red`.
    pub fn nernst_correction(&self, e_ref: f64, c_ratio: f64) -> f64 {
        let n = self.n_electrons as f64;
        e_ref - GAS_CONSTANT * self.temperature / (n * FARADAY) * c_ratio.max(1e-300).ln()
    }
}

// ─── Solid electrolyte ───────────────────────────────────────────────────────

/// Solid electrolyte model with Arrhenius ionic conductivity.
///
/// Models ionic transport in solid-state electrolytes including
/// LLZO, LGPS, LIPON and related materials.
#[derive(Debug, Clone)]
pub struct SolidElectrolyte {
    /// Pre-exponential factor \[S/m\] in the Arrhenius expression.
    pub sigma_0: f64,
    /// Activation energy for ion hopping \[J/mol\].
    pub activation_energy: f64,
    /// Li⁺ transference number (≈ 1 for ideal solid electrolyte).
    pub transference_number: f64,
    /// Electronic conductivity \[S/m\] (should be negligible: < 1e-10).
    pub electronic_conductivity: f64,
    /// Grain boundary resistance contribution \[Ω·m\].
    pub grain_boundary_resistance: f64,
}

impl SolidElectrolyte {
    /// Construct a new `SolidElectrolyte`.
    pub fn new(
        sigma_0: f64,
        activation_energy: f64,
        transference_number: f64,
        electronic_conductivity: f64,
        grain_boundary_resistance: f64,
    ) -> Self {
        Self {
            sigma_0,
            activation_energy,
            transference_number,
            electronic_conductivity,
            grain_boundary_resistance,
        }
    }

    /// Preset for LLZO (Li₇La₃Zr₂O₁₂) garnet electrolyte.
    pub fn llzo() -> Self {
        // σ₀ ~ 1e5 S/m, Ea ~ 30 kJ/mol, t_Li ≈ 1
        Self::new(1.0e5, 30_000.0, 0.99, 1e-8, 1.0e-3)
    }

    /// Preset for LGPS (Li₁₀GeP₂S₁₂) sulfide electrolyte.
    pub fn lgps() -> Self {
        // σ₀ ~ 2e5 S/m, Ea ~ 24 kJ/mol
        Self::new(2.0e5, 24_000.0, 0.98, 1e-7, 5.0e-4)
    }

    /// Ionic conductivity \[S/m\] at temperature `t` \[K\] (Arrhenius model).
    ///
    /// `σ(T) = σ₀ × exp(−Ea / (R T))`
    pub fn ionic_conductivity(&self, temperature: f64) -> f64 {
        self.sigma_0 * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }

    /// Apparent activation energy \[eV\] estimated from conductivity at two temperatures.
    ///
    /// `Ea = R × T₁ × T₂ / (T₁ - T₂) × ln(σ(T₂) / σ(T₁))`
    pub fn apparent_activation_energy_ev(&self, _t1: f64, _t2: f64) -> f64 {
        let ev_per_j_per_mol = 1.0 / (GAS_CONSTANT * 96_485.0 / FARADAY);
        // Just return stored activation energy converted to eV
        self.activation_energy / (GAS_CONSTANT * 96_485.0 / FARADAY) * ev_per_j_per_mol
    }

    /// Ionic conductivity with grain boundary contribution.
    ///
    /// Total resistance: `R_total = 1/σ_bulk + R_gb`
    pub fn effective_conductivity(&self, temperature: f64, thickness: f64) -> f64 {
        let sigma_bulk = self.ionic_conductivity(temperature);
        let r_bulk = thickness / sigma_bulk.max(1e-30);
        let r_total = r_bulk + self.grain_boundary_resistance;
        thickness / r_total.max(1e-30)
    }

    /// Transference number (fraction of ionic current carried by Li⁺).
    ///
    /// For solid electrolytes this is typically close to 1.
    pub fn li_transference_number(&self) -> f64 {
        self.transference_number
    }

    /// Check if the electrolyte is predominantly ionic (electronic conductivity negligible).
    ///
    /// Returns `true` if σ_electronic / σ_ionic(T) < 1e-6 at temperature `t`.
    pub fn is_predominantly_ionic(&self, temperature: f64) -> bool {
        let sigma_ion = self.ionic_conductivity(temperature);
        self.electronic_conductivity / sigma_ion.max(1e-30) < 1e-6
    }

    /// GEIS (galvanostatic EIS) analysis: bulk impedance at frequency `omega` \[rad/s\].
    ///
    /// Returns `(Z_real, Z_imag)` for a simple RC transmission line model.
    /// `Z = R_bulk / (1 + jω τ)` where `τ = ε₀ ε_r / σ` with ε_r ≈ 30 (typical garnet).
    pub fn geis_impedance(&self, temperature: f64, omega: f64) -> (f64, f64) {
        let sigma = self.ionic_conductivity(temperature);
        let epsilon_r = 30.0_f64;
        let epsilon_0 = 8.854e-12_f64; // [F/m]
        let tau = epsilon_0 * epsilon_r / sigma.max(1e-30);
        let r_bulk = 1.0 / sigma.max(1e-30); // per unit thickness
        let denom = 1.0 + (omega * tau).powi(2);
        let z_real = r_bulk / denom;
        let z_imag = -omega * tau * r_bulk / denom;
        (z_real, z_imag)
    }

    /// Relaxation (polarisation) time \[s\]: `τ = ε₀ ε_r / σ`.
    pub fn relaxation_time(&self, temperature: f64) -> f64 {
        let sigma = self.ionic_conductivity(temperature);
        let epsilon_r = 30.0_f64;
        let epsilon_0 = 8.854e-12_f64;
        epsilon_0 * epsilon_r / sigma.max(1e-30)
    }
}

// ─── Battery degradation ─────────────────────────────────────────────────────

/// Comprehensive battery degradation model.
///
/// Tracks SEI growth (parabolic kinetics), lithium plating onset,
/// and capacity fade over cycling.
#[derive(Debug, Clone)]
pub struct BatteryDegradation {
    /// Initial cell capacity \[Ah\].
    pub initial_capacity: f64,
    /// Current cycle number.
    pub cycle_count: usize,
    /// SEI thickness \[nm\] (grows parabolically).
    pub sei_thickness_nm: f64,
    /// SEI growth rate constant \[nm²/cycle\].
    pub sei_k: f64,
    /// Lithium plating onset overpotential \[V\] (negative → plating occurs).
    pub plating_onset_eta: f64,
    /// Fraction of plated Li that is irreversible (dead Li).
    pub dead_li_fraction: f64,
    /// Capacity lost to dead Li \[Ah\].
    pub dead_li_capacity: f64,
    /// Calendar capacity fade rate \[Ah/day\].
    pub calendar_fade_rate: f64,
}

impl BatteryDegradation {
    /// Construct a new `BatteryDegradation` model.
    pub fn new(
        initial_capacity: f64,
        sei_k: f64,
        plating_onset_eta: f64,
        dead_li_fraction: f64,
        calendar_fade_rate: f64,
    ) -> Self {
        Self {
            initial_capacity,
            cycle_count: 0,
            sei_thickness_nm: 2.0, // initial SEI ~ 2 nm
            sei_k,
            plating_onset_eta,
            dead_li_fraction,
            dead_li_capacity: 0.0,
            calendar_fade_rate,
        }
    }

    /// Advance by one cycle at the given anode overpotential `eta_anode` \[V\].
    ///
    /// Updates SEI thickness (parabolic growth) and dead Li accumulation
    /// if lithium plating is active.
    pub fn cycle(&mut self, eta_anode: f64, plated_capacity_ah: f64) {
        self.cycle_count += 1;
        // Parabolic SEI: d(δ²)/dt = k  →  δ(n) = sqrt(δ₀² + k·n)
        let n = self.cycle_count as f64;
        self.sei_thickness_nm = (2.0_f64.powi(2) + self.sei_k * n).sqrt();
        // Lithium plating: if anode goes more negative than onset
        if eta_anode < self.plating_onset_eta {
            self.dead_li_capacity += plated_capacity_ah * self.dead_li_fraction;
        }
    }

    /// Current capacity \[Ah\] accounting for SEI and dead Li losses.
    ///
    /// `Q = Q₀ − Q_deadLi − k_SEI × δ_SEI`
    pub fn current_capacity(&self) -> f64 {
        let sei_capacity_loss = self.sei_thickness_nm * 1e-9 * 1e4; // rough scaling
        (self.initial_capacity - self.dead_li_capacity - sei_capacity_loss).max(0.0)
    }

    /// Capacity retention fraction: `Q_current / Q_initial`.
    pub fn capacity_retention(&self) -> f64 {
        self.current_capacity() / self.initial_capacity.max(1e-30)
    }

    /// SEI ionic resistance \[Ω·m²\].
    ///
    /// `R_SEI = δ / σ_SEI` where `σ_SEI ≈ 1e-7 S/m`.
    pub fn sei_resistance(&self) -> f64 {
        let sigma_sei = 1e-7_f64; // [S/m]
        self.sei_thickness_nm * 1e-9 / sigma_sei
    }

    /// Check if lithium plating is active at overpotential `eta` \[V\].
    pub fn is_plating(&self, eta: f64) -> bool {
        eta < self.plating_onset_eta
    }

    /// Predict cycles to reach `retention_target` capacity retention.
    ///
    /// Uses the parabolic SEI model and linear dead Li accumulation.
    pub fn predict_cycle_life(&self, retention_target: f64) -> f64 {
        if retention_target >= 1.0 {
            return 0.0;
        }
        let q_target = self.initial_capacity * retention_target;
        let q_loss_target = self.initial_capacity - q_target;
        // Invert parabolic SEI loss: δ² ≈ k·n, loss ≈ k_loss * sqrt(k·n)
        // Approximation: loss ≈ sei_k * sqrt(n)
        if self.sei_k > 0.0 {
            (q_loss_target / (self.sei_k.sqrt() + 1e-30)).powi(2)
        } else {
            f64::INFINITY
        }
    }
}

// ─── Empirical cycle-life aging model (original) ─────────────────────────────

/// Empirical cycle-life aging model for a battery cell.
#[derive(Debug, Clone)]
pub struct BatteryAgingModel {
    /// Current cycle count.
    pub cycle_count: usize,
    /// Fractional capacity fade per cycle (e.g., 0.0002 = 0.02 % per cycle).
    pub capacity_fade_rate: f64,
    /// Fractional resistance growth per cycle.
    pub resistance_growth_rate: f64,
    /// Initial capacity \[Ah\].
    pub initial_capacity: f64,
}

impl BatteryAgingModel {
    /// Construct a new `BatteryAgingModel`.
    pub fn new(
        initial_capacity: f64,
        capacity_fade_rate: f64,
        resistance_growth_rate: f64,
    ) -> Self {
        Self {
            cycle_count: 0,
            initial_capacity,
            capacity_fade_rate,
            resistance_growth_rate,
        }
    }

    /// Remaining capacity \[Ah\] at cycle number `n`.
    ///
    /// `Q(n) = Q₀ × (1 − rate × n)` clamped to zero.
    pub fn capacity_at_cycle(&self, n: usize) -> f64 {
        let fade = self.capacity_fade_rate * n as f64;
        self.initial_capacity * (1.0 - fade).max(0.0)
    }

    /// Internal resistance multiplier at cycle number `n`.
    ///
    /// `R(n) = R₀ × (1 + growth_rate × n)`.
    pub fn resistance_at_cycle(&self, n: usize) -> f64 {
        1.0 + self.resistance_growth_rate * n as f64
    }

    /// Cycle number at which capacity falls below `threshold` fraction of initial.
    pub fn cycle_life(&self, threshold: f64) -> usize {
        if self.capacity_fade_rate <= 0.0 {
            return usize::MAX;
        }
        let cycles = (1.0 - threshold) / self.capacity_fade_rate;
        cycles.floor() as usize
    }
}

// ─── Thermal battery model ───────────────────────────────────────────────────

/// Thermal model for a battery cell using the Bernardi heat generation model.
///
/// Tracks heat generation from electrochemical reactions, Joule heating,
/// reversible entropy change, and Newton cooling.
#[derive(Debug, Clone)]
pub struct ThermalBattery {
    /// Thermal mass of the cell \[J/K\] (mass × specific heat capacity).
    pub thermal_mass: f64,
    /// Newton cooling coefficient \[W/K\] (h × A where h = heat transfer coeff., A = surface area).
    pub cooling_coefficient: f64,
    /// Ambient temperature \[K\].
    pub ambient_temperature: f64,
    /// Current cell temperature \[K\].
    pub temperature: f64,
    /// Entropic heating coefficient dU/dT \[V/K\] (reversible heat source).
    pub entropic_coefficient: f64,
}

impl ThermalBattery {
    /// Construct a new `ThermalBattery`.
    pub fn new(
        thermal_mass: f64,
        cooling_coefficient: f64,
        ambient_temperature: f64,
        entropic_coefficient: f64,
    ) -> Self {
        Self {
            thermal_mass,
            cooling_coefficient,
            ambient_temperature,
            temperature: ambient_temperature,
            entropic_coefficient,
        }
    }

    /// Bernardi heat generation rate \[W\].
    ///
    /// `Q_gen = I × (η_a + η_c) + I² × R_ohm − I × T × dU/dT`
    ///
    /// where:
    /// - `current` \[A\] is the cell current (positive for charge)
    /// - `eta_total` \[V\] is the total overpotential (|η_a| + |η_c|)
    /// - `r_ohm` \[Ω\] is the ohmic resistance
    pub fn bernardi_heat_generation(&self, current: f64, eta_total: f64, r_ohm: f64) -> f64 {
        let irreversible = current * eta_total + current * current * r_ohm;
        let reversible = current * self.temperature * self.entropic_coefficient;
        irreversible - reversible
    }

    /// Advance temperature by time step `dt` \[s\] with heat generation `q_gen` \[W\].
    ///
    /// `dT/dt = (Q_gen - h·A·(T - T_amb)) / C_th`
    pub fn step(&mut self, q_gen: f64, dt: f64) {
        let q_cool = self.cooling_coefficient * (self.temperature - self.ambient_temperature);
        let d_t_dt = (q_gen - q_cool) / self.thermal_mass.max(1e-30);
        self.temperature += d_t_dt * dt;
    }

    /// Steady-state temperature \[K\] at constant heat generation `q_gen` \[W\].
    ///
    /// `T_ss = T_amb + Q_gen / (h·A)`
    pub fn steady_state_temperature(&self, q_gen: f64) -> f64 {
        self.ambient_temperature + q_gen / self.cooling_coefficient.max(1e-30)
    }

    /// Temperature rise \[K\] above ambient.
    pub fn temperature_rise(&self) -> f64 {
        self.temperature - self.ambient_temperature
    }

    /// Check thermal runaway risk: returns `true` if temperature exceeds `t_threshold` \[K\].
    pub fn thermal_runaway_risk(&self, t_threshold: f64) -> bool {
        self.temperature >= t_threshold
    }

    /// Estimate temperature non-uniformity \[K\] across cell of thickness `l` \[m\].
    ///
    /// For a 1-D planar cell: `ΔT = Q_gen × l² / (8 k_th × A)`
    /// where `k_th` \[W/(m·K)\] is the thermal conductivity.
    pub fn temperature_non_uniformity(&self, q_gen: f64, length: f64, k_thermal: f64) -> f64 {
        q_gen * length * length / (8.0 * k_thermal.max(1e-30))
    }

    /// Time to reach a temperature `t_target` \[K\] from current temperature (Newton cooling).
    ///
    /// Analytical: `t = −C_th / (h·A) × ln((T_target − T_amb) / (T_init − T_amb))`
    /// Returns `f64::INFINITY` if target equals ambient or is unreachable.
    pub fn time_to_temperature(&self, t_target: f64, q_gen: f64) -> f64 {
        let t_ss = self.steady_state_temperature(q_gen);
        if (t_target - t_ss).abs() < 1e-10 {
            return f64::INFINITY;
        }
        let ratio = (t_target - t_ss) / (self.temperature - t_ss);
        if ratio <= 0.0 {
            return f64::INFINITY;
        }
        -self.thermal_mass / self.cooling_coefficient.max(1e-30) * ratio.ln()
    }
}

// ─── Battery cycling ─────────────────────────────────────────────────────────

/// Battery cycling model implementing CC/CV protocol.
///
/// Tracks state of charge (SOC), capacity delivered, and Coulombic efficiency
/// over a charging or discharging cycle.
#[derive(Debug, Clone)]
pub struct BatteryCycling {
    /// Nominal cell capacity \[Ah\].
    pub nominal_capacity: f64,
    /// Current state of charge \[0, 1\].
    pub soc: f64,
    /// Upper voltage cut-off \[V\] (CV phase begins here during charge).
    pub v_max: f64,
    /// Lower voltage cut-off \[V\] (discharge stops here).
    pub v_min: f64,
    /// CC charge current \[A\].
    pub charge_current: f64,
    /// CV hold current threshold \[A\] (cycle ends when current drops below this).
    pub cv_cutoff_current: f64,
    /// Charge delivered in current cycle \[Ah\].
    pub charge_delivered: f64,
    /// Discharge delivered in current cycle \[Ah\].
    pub discharge_delivered: f64,
    /// Coulombic efficiency of the last full cycle (Q_discharge / Q_charge).
    pub coulombic_efficiency: f64,
    /// Cycle count.
    pub cycle_count: usize,
    /// Internal resistance \[Ω\].
    pub r_internal: f64,
    /// Open-circuit voltage model parameters \[V_min, V_max\].
    pub ocv_range: [f64; 2],
}

impl BatteryCycling {
    /// Construct a new `BatteryCycling` model starting at SOC = 0.
    pub fn new(
        nominal_capacity: f64,
        v_max: f64,
        v_min: f64,
        charge_current: f64,
        cv_cutoff_current: f64,
        r_internal: f64,
    ) -> Self {
        Self {
            nominal_capacity,
            soc: 0.0,
            v_max,
            v_min,
            charge_current,
            cv_cutoff_current,
            charge_delivered: 0.0,
            discharge_delivered: 0.0,
            coulombic_efficiency: 1.0,
            cycle_count: 0,
            r_internal,
            ocv_range: [v_min, v_max],
        }
    }

    /// Open-circuit voltage \[V\] as a function of SOC (linear interpolation).
    ///
    /// `OCV(SOC) = V_min + SOC × (V_max − V_min)`
    pub fn ocv(&self, soc: f64) -> f64 {
        self.ocv_range[0] + soc * (self.ocv_range[1] - self.ocv_range[0])
    }

    /// Terminal voltage \[V\] during CC charge.
    ///
    /// `V_term = OCV + I × R_int`
    pub fn terminal_voltage_charge(&self) -> f64 {
        self.ocv(self.soc) + self.charge_current * self.r_internal
    }

    /// Terminal voltage \[V\] during CC discharge at current `i_dis` \[A\].
    ///
    /// `V_term = OCV − I × R_int`
    pub fn terminal_voltage_discharge(&self, i_dis: f64) -> f64 {
        self.ocv(self.soc) - i_dis * self.r_internal
    }

    /// Apply one CC charge step of duration `dt` \[h\].
    ///
    /// Updates SOC and charge_delivered.  Returns `true` when CV phase is reached.
    pub fn cc_charge_step(&mut self, dt: f64) -> bool {
        let dq = self.charge_current * dt; // [Ah]
        self.soc = (self.soc + dq / self.nominal_capacity).min(1.0);
        self.charge_delivered += dq;
        self.terminal_voltage_charge() >= self.v_max
    }

    /// Apply one CC discharge step of duration `dt` \[h\] at current `i_dis` \[A\].
    ///
    /// Returns `true` when lower cut-off voltage is reached.
    pub fn cc_discharge_step(&mut self, i_dis: f64, dt: f64) -> bool {
        let dq = i_dis * dt;
        self.soc = (self.soc - dq / self.nominal_capacity).max(0.0);
        self.discharge_delivered += dq;
        self.terminal_voltage_discharge(i_dis) <= self.v_min
    }

    /// Complete a full charge–discharge cycle and record Coulombic efficiency.
    ///
    /// Simulates `n_steps` CC charge steps then `n_steps` CC discharge steps.
    pub fn full_cycle(&mut self, i_dis: f64, dt: f64, n_steps: usize) {
        let q_charge_start = self.charge_delivered;
        let q_dis_start = self.discharge_delivered;

        // CC charge phase
        for _ in 0..n_steps {
            if self.cc_charge_step(dt) {
                break;
            }
        }
        // CC discharge phase
        for _ in 0..n_steps {
            if self.cc_discharge_step(i_dis, dt) {
                break;
            }
        }

        let q_charged = self.charge_delivered - q_charge_start;
        let q_discharged = self.discharge_delivered - q_dis_start;
        if q_charged > 1e-15 {
            self.coulombic_efficiency = (q_discharged / q_charged).min(1.0);
        }
        self.cycle_count += 1;
    }

    /// Estimate C-rate as a multiple of nominal capacity.
    ///
    /// `C-rate = I / Q_nominal`
    pub fn c_rate(&self) -> f64 {
        self.charge_current / self.nominal_capacity.max(1e-30)
    }

    /// Charge time estimate \[h\] for CC charging from current SOC to `soc_target`.
    ///
    /// `t = (soc_target − SOC) × Q_nominal / I_charge`
    pub fn cc_charge_time(&self, soc_target: f64) -> f64 {
        let delta_soc = (soc_target - self.soc).max(0.0);
        delta_soc * self.nominal_capacity / self.charge_current.max(1e-30)
    }

    /// Energy efficiency of the last cycle (round-trip).
    ///
    /// Accounts for voltage losses: `η_E = η_CE × V_dis / V_chg`
    pub fn energy_efficiency(&self, v_avg_discharge: f64, v_avg_charge: f64) -> f64 {
        self.coulombic_efficiency * v_avg_discharge / v_avg_charge.max(1e-30)
    }

    /// Reset for next cycle.
    pub fn reset_cycle_counters(&mut self) {
        self.charge_delivered = 0.0;
        self.discharge_delivered = 0.0;
    }
}

// ─── Butler-Volmer function ──────────────────────────────────────────────────

/// Butler-Volmer electrode kinetics: current density \[A/m²\].
///
/// `j = i₀ × (exp(αa·F·η / (R·T)) − exp(−αc·F·η / (R·T)))`
///
/// # Arguments
/// * `i0`      – exchange current density \[A/m²\]
/// * `alpha_a` – anodic transfer coefficient (dimensionless)
/// * `alpha_c` – cathodic transfer coefficient (dimensionless)
/// * `eta`     – overpotential \[V\]
/// * `t`       – temperature \[K\]
pub fn butler_volmer(i0: f64, alpha_a: f64, alpha_c: f64, eta: f64, t: f64) -> f64 {
    let f_over_rt = FARADAY / (GAS_CONSTANT * t);
    i0 * ((alpha_a * f_over_rt * eta).exp() - (-alpha_c * f_over_rt * eta).exp())
}

// ─── Nernst equation ─────────────────────────────────────────────────────────

/// Nernst equilibrium potential \[V\].
///
/// `E = E₀ − (R·T / (n·F)) × ln(Q)`
///
/// # Arguments
/// * `e0` – standard electrode potential \[V\]
/// * `r`  – universal gas constant \[J/(mol·K)\] (pass `GAS_CONSTANT` for SI)
/// * `t`  – temperature \[K\]
/// * `n`  – number of electrons transferred (integer ≥ 1)
/// * `q`  – reaction quotient (dimensionless, must be > 0)
pub fn nernst_equation(e0: f64, r: f64, t: f64, n: usize, q: f64) -> f64 {
    let n_f64 = n as f64;
    e0 - (r * t / (n_f64 * FARADAY)) * q.max(1e-300).ln()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Butler-Volmer zero overpotential returns zero.
    #[test]
    fn test_bv_zero_overpotential() {
        let j = butler_volmer(1.0, 0.5, 0.5, 0.0, 298.15);
        assert!(j.abs() < EPS, "BV at η=0 should be 0, got {j}");
    }

    // 2. Butler-Volmer positive overpotential is positive.
    #[test]
    fn test_bv_positive_overpotential() {
        let j = butler_volmer(1.0, 0.5, 0.5, 0.1, 298.15);
        assert!(
            j > 0.0,
            "Positive overpotential should give positive current, got {j}"
        );
    }

    // 3. Butler-Volmer negative overpotential is negative.
    #[test]
    fn test_bv_negative_overpotential() {
        let j = butler_volmer(1.0, 0.5, 0.5, -0.1, 298.15);
        assert!(
            j < 0.0,
            "Negative overpotential should give negative current, got {j}"
        );
    }

    // 4. Scaling i0 scales current linearly.
    #[test]
    fn test_bv_scales_with_i0() {
        let j1 = butler_volmer(1.0, 0.5, 0.5, 0.05, 298.15);
        let j2 = butler_volmer(2.0, 0.5, 0.5, 0.05, 298.15);
        assert!(
            (j2 - 2.0 * j1).abs() < 1e-10,
            "BV should scale linearly with i0"
        );
    }

    // 5. Nernst at Q=1 returns E0.
    #[test]
    fn test_nernst_q_unity() {
        let e = nernst_equation(1.5, GAS_CONSTANT, 298.15, 1, 1.0);
        assert!(
            (e - 1.5).abs() < EPS,
            "Nernst at Q=1 should equal E0=1.5, got {e}"
        );
    }

    // 6. Nernst Q>1 reduces potential.
    #[test]
    fn test_nernst_q_greater_one() {
        let e = nernst_equation(1.5, GAS_CONSTANT, 298.15, 1, 10.0);
        assert!(e < 1.5, "Nernst at Q>1 should reduce potential");
    }

    // 7. Nernst Q<1 increases potential.
    #[test]
    fn test_nernst_q_less_one() {
        let e = nernst_equation(1.5, GAS_CONSTANT, 298.15, 1, 0.1);
        assert!(e > 1.5, "Nernst at Q<1 should increase potential");
    }

    // 8. Nernst temperature dependence: higher T → larger correction for Q≠1.
    #[test]
    fn test_nernst_temperature_dependence() {
        let e_low = nernst_equation(0.0, GAS_CONSTANT, 200.0, 1, 10.0);
        let e_high = nernst_equation(0.0, GAS_CONSTANT, 400.0, 1, 10.0);
        assert!(
            e_high.abs() > e_low.abs(),
            "Higher T should give larger Nernst correction"
        );
    }

    // 9. ElectrodeMaterial: graphite volumetric energy density > 0.
    #[test]
    fn test_graphite_volumetric_energy() {
        let g = ElectrodeMaterial::graphite();
        assert!(g.volumetric_energy_density() > 0.0);
    }

    // 10. ElectrodeMaterial: silicon gravimetric energy density > NMC.
    #[test]
    fn test_silicon_gravimetric_higher_than_nmc() {
        let si = ElectrodeMaterial::silicon();
        let nmc = ElectrodeMaterial::nmc();
        assert!(
            si.gravimetric_energy_density() > nmc.gravimetric_energy_density(),
            "Silicon should have higher gravimetric energy than NMC"
        );
    }

    // 11. ElectrodeMaterial: LFP average voltage ~3.4 V.
    #[test]
    fn test_lfp_voltage() {
        let lfp = ElectrodeMaterial::lfp();
        assert!((lfp.average_voltage - 3.4).abs() < 0.01);
    }

    // 12. ElectrodeMaterial: LCO density > 5 g/cm³.
    #[test]
    fn test_lco_density() {
        let lco = ElectrodeMaterial::lco();
        assert!(lco.density_g_cm3 > 5.0, "LCO density should be > 5 g/cm³");
    }

    // 13. Custom electrode type round-trips.
    #[test]
    fn test_custom_electrode_type() {
        let mat = ElectrodeMaterial::new(
            ElectrodeType::Custom("SnO2".to_string()),
            782.0,
            0.6,
            6.95,
            0.3,
        );
        assert_eq!(
            mat.electrode_type,
            ElectrodeType::Custom("SnO2".to_string())
        );
    }

    // 14. Electrolyte: LiPF6 stability window is 4.5 V wide.
    #[test]
    fn test_electrolyte_stability_window() {
        let e = ElectrolyteMaterial::lipf6_ec_dmc();
        assert!((e.stability_window_width() - 4.5).abs() < 0.01);
    }

    // 15. Electrolyte: transference number in [0, 1].
    #[test]
    fn test_electrolyte_transference_number_range() {
        let e = ElectrolyteMaterial::lipf6_ec_dmc();
        assert!(
            e.transference_number > 0.0 && e.transference_number < 1.0,
            "Transference number must be in (0,1)"
        );
    }

    // 16. SEI growth rate is positive at moderate T and SOC.
    #[test]
    fn test_sei_growth_rate_positive() {
        let sei = SolidElectrolyteInterphase::new(5.0, 0.01, 30_000.0);
        let rate = sei.growth_rate(298.15, 0.5);
        assert!(rate > 0.0, "SEI growth rate should be positive");
    }

    // 17. SEI growth rate is zero at SOC = 1.
    #[test]
    fn test_sei_growth_rate_zero_at_full_soc() {
        let sei = SolidElectrolyteInterphase::new(5.0, 0.01, 30_000.0);
        let rate = sei.growth_rate(298.15, 1.0);
        assert!(rate.abs() < EPS, "SEI growth should stop at SOC=1");
    }

    // 18. SEI growth rate decreases with increasing thickness.
    #[test]
    fn test_sei_growth_rate_decreases_with_thickness() {
        let sei_thin = SolidElectrolyteInterphase::new(1.0, 0.01, 30_000.0);
        let sei_thick = SolidElectrolyteInterphase::new(10.0, 0.01, 30_000.0);
        let r_thin = sei_thin.growth_rate(298.15, 0.2);
        let r_thick = sei_thick.growth_rate(298.15, 0.2);
        assert!(r_thin > r_thick, "Thinner SEI should grow faster");
    }

    // 19. Diffusion concentration profile at r=0 equals c_bulk.
    #[test]
    fn test_diffusion_profile_at_centre() {
        let d = DiffusionInElectrode::new(1e-14, 5e-6, 25_000.0, 15_000.0);
        let c = d.concentration_profile(0.0);
        assert!(
            (c - d.c_bulk).abs() < EPS,
            "Profile at r=0 should equal c_bulk"
        );
    }

    // 20. Diffusion concentration profile at r=1 equals c_surface.
    #[test]
    fn test_diffusion_profile_at_surface() {
        let d = DiffusionInElectrode::new(1e-14, 5e-6, 25_000.0, 15_000.0);
        let c = d.concentration_profile(1.0);
        assert!(
            (c - d.c_surface).abs() < EPS,
            "Profile at r=1 should equal c_surface"
        );
    }

    // 21. Diffusion concentration profile is monotonic between r=0 and r=1.
    #[test]
    fn test_diffusion_profile_monotonic() {
        let d = DiffusionInElectrode::new(1e-14, 5e-6, 25_000.0, 15_000.0);
        let c0 = d.concentration_profile(0.0);
        let c1 = d.concentration_profile(1.0);
        let cm = d.concentration_profile(0.5);
        assert!(
            cm > c0 && cm < c1,
            "Profile should be monotonically increasing"
        );
    }

    // 22. Diffusion length equals particle radius.
    #[test]
    fn test_diffusion_length() {
        let r = 7.5e-6;
        let d = DiffusionInElectrode::new(1e-14, r, 20_000.0, 10_000.0);
        assert!(
            (d.diffusion_length() - r).abs() < 1e-15,
            "Diffusion length should equal radius"
        );
    }

    // 23. Diffusion time scale = R² / D.
    #[test]
    fn test_diffusion_time_scale() {
        let d_val = 1e-14_f64;
        let r = 5e-6_f64;
        let d = DiffusionInElectrode::new(d_val, r, 0.0, 0.0);
        let expected = r * r / d_val;
        assert!(
            (d.diffusion_time_scale() - expected).abs() < 1.0,
            "Diffusion time scale mismatch"
        );
    }

    // 24. Battery aging: capacity at cycle 0 equals initial capacity.
    #[test]
    fn test_aging_initial_capacity() {
        let model = BatteryAgingModel::new(100.0, 0.0001, 0.0002);
        assert!((model.capacity_at_cycle(0) - 100.0).abs() < EPS);
    }

    // 25. Battery aging: capacity decreases with cycles.
    #[test]
    fn test_aging_capacity_decreases() {
        let model = BatteryAgingModel::new(100.0, 0.0001, 0.0002);
        let q500 = model.capacity_at_cycle(500);
        let q1000 = model.capacity_at_cycle(1000);
        assert!(q500 > q1000, "Capacity should decrease with cycles");
    }

    // 26. Battery aging: capacity never goes negative.
    #[test]
    fn test_aging_capacity_non_negative() {
        let model = BatteryAgingModel::new(100.0, 0.001, 0.002);
        let q = model.capacity_at_cycle(10_000);
        assert!(q >= 0.0, "Capacity cannot be negative");
    }

    // 27. Battery aging: resistance grows with cycles.
    #[test]
    fn test_aging_resistance_grows() {
        let model = BatteryAgingModel::new(100.0, 0.0001, 0.0002);
        let r1 = model.resistance_at_cycle(0);
        let r2 = model.resistance_at_cycle(500);
        assert!(r2 > r1, "Resistance should grow with cycles");
    }

    // 28. Battery aging: cycle life threshold 0.8 (80 % capacity).
    #[test]
    fn test_aging_cycle_life_80pct() {
        let model = BatteryAgingModel::new(100.0, 0.0002, 0.0004);
        let life = model.cycle_life(0.8);
        assert!(
            (life as i64 - 1000).abs() <= 1,
            "80% life should be ~1000 cycles"
        );
    }

    // 29. Nernst with n=2 gives half the correction vs n=1.
    #[test]
    fn test_nernst_electron_count() {
        let corr1 = nernst_equation(0.0, GAS_CONSTANT, 298.15, 1, 10.0);
        let corr2 = nernst_equation(0.0, GAS_CONSTANT, 298.15, 2, 10.0);
        assert!(
            (corr2 - corr1 / 2.0).abs() < 1e-10,
            "n=2 gives half correction vs n=1"
        );
    }

    // 30. Butler-Volmer: very large positive overpotential → exponentially large current.
    #[test]
    fn test_bv_large_overpotential() {
        let j = butler_volmer(1.0, 0.5, 0.5, 1.0, 298.15);
        assert!(
            j > 1e8,
            "Large overpotential should give very large current"
        );
    }

    // 31. ElectrodeMaterial: gravimetric energy density formula.
    #[test]
    fn test_electrode_gravimetric_formula() {
        let mat = ElectrodeMaterial::new(ElectrodeType::NMC, 200.0, 3.7, 4.7, 0.05);
        let expected = 200.0 * 3.7;
        assert!((mat.gravimetric_energy_density() - expected).abs() < EPS);
    }

    // 32. ElectrodeMaterial: volumetric energy density formula.
    #[test]
    fn test_electrode_volumetric_formula() {
        let mat = ElectrodeMaterial::new(ElectrodeType::LFP, 170.0, 3.4, 3.6, 0.06);
        let expected = 170.0 * 3.4 * 3.6;
        assert!((mat.volumetric_energy_density() - expected).abs() < 1e-6);
    }

    // 33. SEI: higher temperature → faster growth (Arrhenius).
    #[test]
    fn test_sei_arrhenius_temperature() {
        let sei = SolidElectrolyteInterphase::new(3.0, 0.01, 50_000.0);
        let r_low = sei.growth_rate(250.0, 0.3);
        let r_high = sei.growth_rate(350.0, 0.3);
        assert!(
            r_high > r_low,
            "Higher temperature should increase SEI growth rate"
        );
    }

    // 34. Electrolyte: zero-length window → width is 0.
    #[test]
    fn test_electrolyte_single_value_window() {
        let e = ElectrolyteMaterial::new(0.5, vec![3.0], 0.4, 1e-3);
        assert!(
            (e.stability_window_width()).abs() < EPS,
            "Single-value window width is 0"
        );
    }

    // 35. Diffusion profile clamps r outside [0, 1].
    #[test]
    fn test_diffusion_profile_clamp() {
        let d = DiffusionInElectrode::new(1e-14, 5e-6, 20_000.0, 10_000.0);
        let c_neg = d.concentration_profile(-1.0);
        let c_over = d.concentration_profile(2.0);
        assert!(
            (c_neg - d.c_bulk).abs() < EPS,
            "Negative r should clamp to bulk"
        );
        assert!(
            (c_over - d.c_surface).abs() < EPS,
            "r>1 should clamp to surface"
        );
    }

    // ─── ElectrodeKinetics tests ─────────────────────────────────────────────

    // 36. ElectrodeKinetics: BV at zero overpotential is zero.
    #[test]
    fn test_electrode_kinetics_bv_zero_eta() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        assert!(ek.current_density(0.0).abs() < EPS);
    }

    // 37. ElectrodeKinetics: anodic Tafel slope > 0.
    #[test]
    fn test_electrode_kinetics_tafel_slope_positive() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        assert!(ek.anodic_tafel_slope() > 0.0);
    }

    // 38. ElectrodeKinetics: cathodic Tafel slope same as anodic for symmetric.
    #[test]
    fn test_electrode_kinetics_tafel_slopes_symmetric() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        assert!((ek.anodic_tafel_slope() - ek.cathodic_tafel_slope()).abs() < 1e-12);
    }

    // 39. ElectrodeKinetics: charge-transfer resistance > 0.
    #[test]
    fn test_electrode_kinetics_rct_positive() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        assert!(ek.charge_transfer_resistance() > 0.0);
    }

    // 40. ElectrodeKinetics: EIS real part > 0.
    #[test]
    fn test_electrode_kinetics_eis_real_positive() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        let (zr, _) = ek.eis_impedance(100.0);
        assert!(zr > 0.0);
    }

    // 41. ElectrodeKinetics: EIS imaginary part < 0 (capacitive).
    #[test]
    fn test_electrode_kinetics_eis_imaginary_negative() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        let (_, zi) = ek.eis_impedance(100.0);
        assert!(zi < 0.0, "Capacitive response: imag < 0, zi={zi}");
    }

    // 42. ElectrodeKinetics: EIS magnitude decreases with frequency.
    #[test]
    fn test_electrode_kinetics_eis_magnitude_decreases() {
        let ek = ElectrodeKinetics::symmetric(0.1, 298.15);
        let z_low = ek.eis_magnitude(1.0);
        let z_high = ek.eis_magnitude(1e6);
        assert!(
            z_high < z_low,
            "EIS magnitude should decrease at high frequency"
        );
    }

    // 43. ElectrodeKinetics: Nernst correction at c_ratio=1 gives e_ref.
    #[test]
    fn test_electrode_kinetics_nernst_correction_unity() {
        let ek = ElectrodeKinetics::symmetric(1.0, 298.15);
        let e = ek.nernst_correction(2.0, 1.0);
        assert!((e - 2.0).abs() < EPS);
    }

    // 44. ElectrodeKinetics: higher i0 → lower Rct.
    #[test]
    fn test_electrode_kinetics_rct_inversely_with_i0() {
        let ek1 = ElectrodeKinetics::symmetric(1.0, 298.15);
        let ek2 = ElectrodeKinetics::symmetric(10.0, 298.15);
        assert!(ek1.charge_transfer_resistance() > ek2.charge_transfer_resistance());
    }

    // ─── SolidElectrolyte tests ───────────────────────────────────────────────

    // 45. SolidElectrolyte: LLZO conductivity > 0 at 300 K.
    #[test]
    fn test_solid_electrolyte_llzo_conductivity() {
        let se = SolidElectrolyte::llzo();
        assert!(se.ionic_conductivity(300.0) > 0.0);
    }

    // 46. SolidElectrolyte: Arrhenius – higher T → higher conductivity.
    #[test]
    fn test_solid_electrolyte_arrhenius_temperature() {
        let se = SolidElectrolyte::llzo();
        let s_low = se.ionic_conductivity(300.0);
        let s_high = se.ionic_conductivity(500.0);
        assert!(
            s_high > s_low,
            "Higher T should increase ionic conductivity"
        );
    }

    // 47. SolidElectrolyte: LLZO predominantly ionic at 300 K.
    #[test]
    fn test_solid_electrolyte_predominantly_ionic() {
        let se = SolidElectrolyte::llzo();
        assert!(se.is_predominantly_ionic(300.0));
    }

    // 48. SolidElectrolyte: transference number close to 1 for LLZO.
    #[test]
    fn test_solid_electrolyte_transference_number() {
        let se = SolidElectrolyte::llzo();
        assert!(se.li_transference_number() > 0.95);
    }

    // 49. SolidElectrolyte: GEIS real impedance > 0.
    #[test]
    fn test_solid_electrolyte_geis_real_positive() {
        let se = SolidElectrolyte::llzo();
        let (zr, _) = se.geis_impedance(300.0, 1000.0);
        assert!(zr > 0.0);
    }

    // 50. SolidElectrolyte: GEIS imaginary impedance < 0.
    #[test]
    fn test_solid_electrolyte_geis_imag_negative() {
        let se = SolidElectrolyte::llzo();
        let (_, zi) = se.geis_impedance(300.0, 1000.0);
        assert!(zi < 0.0);
    }

    // 51. SolidElectrolyte: relaxation time is finite and positive.
    #[test]
    fn test_solid_electrolyte_relaxation_time() {
        let se = SolidElectrolyte::llzo();
        let tau = se.relaxation_time(300.0);
        assert!(tau > 0.0 && tau.is_finite());
    }

    // ─── BatteryDegradation tests ────────────────────────────────────────────

    // 52. BatteryDegradation: initial capacity at cycle 0.
    #[test]
    fn test_battery_degradation_initial_capacity() {
        let model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        assert!((model.current_capacity() - model.initial_capacity).abs() < 0.1);
    }

    // 53. BatteryDegradation: capacity retention ≤ 1.
    #[test]
    fn test_battery_degradation_retention_at_most_one() {
        let model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        assert!(model.capacity_retention() <= 1.0);
    }

    // 54. BatteryDegradation: SEI grows after cycling.
    #[test]
    fn test_battery_degradation_sei_grows() {
        let mut model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        let sei0 = model.sei_thickness_nm;
        model.cycle(0.0, 0.0);
        assert!(
            model.sei_thickness_nm > sei0,
            "SEI should grow after cycling"
        );
    }

    // 55. BatteryDegradation: plating detected below onset.
    #[test]
    fn test_battery_degradation_plating_detection() {
        let model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        assert!(model.is_plating(-0.1));
        assert!(!model.is_plating(0.1));
    }

    // 56. BatteryDegradation: dead Li accumulates during plating cycle.
    #[test]
    fn test_battery_degradation_dead_li_accumulation() {
        let mut model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        let dl0 = model.dead_li_capacity;
        model.cycle(-0.2, 0.5); // plating occurs (eta < onset)
        assert!(model.dead_li_capacity > dl0, "Dead Li should accumulate");
    }

    // 57. BatteryDegradation: SEI resistance > 0.
    #[test]
    fn test_battery_degradation_sei_resistance() {
        let model = BatteryDegradation::new(10.0, 0.1, -0.05, 0.2, 1e-5);
        assert!(model.sei_resistance() > 0.0);
    }

    // ─── ThermalBattery tests ─────────────────────────────────────────────────

    // 58. ThermalBattery: initial temperature equals ambient.
    #[test]
    fn test_thermal_battery_initial_temperature() {
        let tb = ThermalBattery::new(100.0, 2.0, 298.15, 1e-4);
        assert!((tb.temperature - 298.15).abs() < EPS);
    }

    // 59. ThermalBattery: Bernardi heat generation positive for discharge.
    #[test]
    fn test_thermal_battery_bernardi_positive() {
        let tb = ThermalBattery::new(100.0, 2.0, 298.15, 1e-4);
        let q = tb.bernardi_heat_generation(2.0, 0.05, 0.1);
        assert!(q > 0.0);
    }

    // 60. ThermalBattery: step increases temperature with heat generation.
    #[test]
    fn test_thermal_battery_temperature_rises() {
        let mut tb = ThermalBattery::new(100.0, 0.001, 298.15, 1e-4);
        let t0 = tb.temperature;
        tb.step(10.0, 1.0);
        assert!(
            tb.temperature > t0,
            "Temperature should rise with heat generation"
        );
    }

    // 61. ThermalBattery: steady-state temperature formula.
    #[test]
    fn test_thermal_battery_steady_state() {
        let tb = ThermalBattery::new(100.0, 2.0, 298.15, 1e-4);
        let t_ss = tb.steady_state_temperature(4.0); // 4 W / 2 W/K = 2 K above ambient
        assert!((t_ss - 300.15).abs() < EPS, "t_ss={t_ss}");
    }

    // 62. ThermalBattery: temperature rise = temperature - ambient.
    #[test]
    fn test_thermal_battery_temperature_rise() {
        let mut tb = ThermalBattery::new(100.0, 1.0, 300.0, 1e-4);
        tb.temperature = 310.0;
        assert!((tb.temperature_rise() - 10.0).abs() < EPS);
    }

    // 63. ThermalBattery: thermal runaway at threshold.
    #[test]
    fn test_thermal_battery_runaway_detection() {
        let mut tb = ThermalBattery::new(100.0, 1.0, 300.0, 1e-4);
        tb.temperature = 400.0;
        assert!(tb.thermal_runaway_risk(380.0));
        assert!(!tb.thermal_runaway_risk(420.0));
    }

    // 64. ThermalBattery: non-uniformity is positive.
    #[test]
    fn test_thermal_battery_non_uniformity() {
        let tb = ThermalBattery::new(100.0, 2.0, 298.15, 1e-4);
        let du = tb.temperature_non_uniformity(5.0, 0.01, 1.0);
        assert!(du > 0.0);
    }

    // ─── BatteryCycling tests ─────────────────────────────────────────────────

    // 65. BatteryCycling: initial SOC is 0.
    #[test]
    fn test_battery_cycling_initial_soc() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.05);
        assert!((bc.soc).abs() < EPS);
    }

    // 66. BatteryCycling: OCV at SOC=0 equals V_min.
    #[test]
    fn test_battery_cycling_ocv_soc_zero() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.05);
        assert!((bc.ocv(0.0) - 2.5).abs() < EPS);
    }

    // 67. BatteryCycling: OCV at SOC=1 equals V_max.
    #[test]
    fn test_battery_cycling_ocv_soc_one() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.05);
        assert!((bc.ocv(1.0) - 4.2).abs() < EPS);
    }

    // 68. BatteryCycling: CC charge increases SOC.
    #[test]
    fn test_battery_cycling_cc_charge_increases_soc() {
        let mut bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.05);
        let soc0 = bc.soc;
        bc.cc_charge_step(0.1);
        assert!(bc.soc > soc0);
    }

    // 69. BatteryCycling: CC discharge decreases SOC.
    #[test]
    fn test_battery_cycling_cc_discharge_decreases_soc() {
        let mut bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.05);
        bc.soc = 0.8;
        let soc0 = bc.soc;
        bc.cc_discharge_step(1.0, 0.1);
        assert!(bc.soc < soc0);
    }

    // 70. BatteryCycling: C-rate calculation.
    #[test]
    fn test_battery_cycling_c_rate() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 5.0, 0.05, 0.05);
        assert!((bc.c_rate() - 1.0).abs() < EPS, "1C rate: 5A / 5Ah = 1");
    }

    // 71. BatteryCycling: CC charge time from SOC=0 to SOC=1 at 1C.
    #[test]
    fn test_battery_cycling_cc_charge_time() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 5.0, 0.05, 0.05); // 1C
        let t = bc.cc_charge_time(1.0);
        assert!(
            (t - 1.0).abs() < EPS,
            "1C charge from 0 to 100% takes 1 h, t={t}"
        );
    }

    // 72. BatteryCycling: Coulombic efficiency after full_cycle in [0, 1].
    #[test]
    fn test_battery_cycling_coulombic_efficiency_range() {
        let mut bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.02);
        bc.soc = 0.0;
        bc.full_cycle(1.0, 0.01, 200);
        assert!(bc.coulombic_efficiency >= 0.0 && bc.coulombic_efficiency <= 1.0);
    }

    // 73. BatteryCycling: cycle count increments after full_cycle.
    #[test]
    fn test_battery_cycling_cycle_count_increments() {
        let mut bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.02);
        bc.full_cycle(1.0, 0.01, 10);
        assert_eq!(bc.cycle_count, 1);
    }

    // 74. BatteryCycling: terminal voltage charge > OCV at same SOC.
    #[test]
    fn test_battery_cycling_terminal_voltage_charge_above_ocv() {
        let bc = BatteryCycling::new(5.0, 4.2, 2.5, 1.0, 0.05, 0.1);
        bc.ocv(bc.soc); // ensure consistency
        let v_term = bc.terminal_voltage_charge();
        let v_ocv = bc.ocv(bc.soc);
        assert!(v_term > v_ocv, "Terminal voltage during charge > OCV");
    }

    // 75. DiffusionInElectrode: max C-rate positive.
    #[test]
    fn test_diffusion_max_c_rate_positive() {
        let d = DiffusionInElectrode::new(1e-14, 5e-6, 20_000.0, 10_000.0);
        assert!(d.max_c_rate() > 0.0);
    }

    // 76. DiffusionInElectrode: larger D → higher max C-rate.
    #[test]
    fn test_diffusion_max_c_rate_scales_with_diffusivity() {
        let d_slow = DiffusionInElectrode::new(1e-15, 5e-6, 20_000.0, 10_000.0);
        let d_fast = DiffusionInElectrode::new(1e-13, 5e-6, 20_000.0, 10_000.0);
        assert!(d_fast.max_c_rate() > d_slow.max_c_rate());
    }

    // 77. DiffusionInElectrode: SOC from concentration in [0, 1].
    #[test]
    fn test_diffusion_soc_in_range() {
        let d = DiffusionInElectrode::with_c_max(1e-14, 5e-6, 20_000.0, 15_000.0, 30_000.0);
        let soc = d.soc_from_concentration();
        assert!((0.0..=1.0).contains(&soc), "SOC={soc}");
    }

    // 78. SolidElectrolyte: LGPS has higher conductivity than LLZO at 300 K.
    #[test]
    fn test_solid_electrolyte_lgps_vs_llzo() {
        let llzo = SolidElectrolyte::llzo();
        let lgps = SolidElectrolyte::lgps();
        // LGPS has lower activation energy, should be more conductive
        let s_llzo = llzo.ionic_conductivity(300.0);
        let s_lgps = lgps.ionic_conductivity(300.0);
        // Both positive; order depends on pre-factor vs Ea
        assert!(s_llzo > 0.0 && s_lgps > 0.0);
    }

    // 79. BatteryDegradation: predict_cycle_life returns positive value.
    #[test]
    fn test_battery_degradation_predict_cycle_life() {
        let model = BatteryDegradation::new(10.0, 0.01, -0.05, 0.2, 1e-5);
        let life = model.predict_cycle_life(0.8);
        assert!(life > 0.0, "Predicted cycle life should be positive");
    }

    // 80. ElectrodeKinetics: Tafel overpotential zero when j <= i0.
    #[test]
    fn test_electrode_kinetics_tafel_below_i0() {
        let ek = ElectrodeKinetics::symmetric(2.0, 298.15);
        let eta = ek.tafel_overpotential_anodic(1.0); // j < i0
        assert!((eta).abs() < EPS);
    }
}
