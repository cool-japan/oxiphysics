#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reactive LBM with combustion modelling.
//!
//! This module provides:
//!
//! - **CombustionLbm**: Main solver coupling species transport, temperature, and LBM flow.
//! - **ArrheniusKinetics**: Arrhenius-type reaction rate for premixed/diffusion flames.
//! - **SpeciesField**: Multi-species (fuel, oxidizer, products) concentration fields.
//! - **FlameTracker**: Detects and tracks flame front position.
//! - **DdtModel**: Deflagration-to-detonation transition model.
//! - **FlameSpeed**: Laminar and turbulent flame speed computation.
//! - **LewisEffect**: Lewis number correction to species diffusivity.
//! - **SootModel**: Simple soot formation/oxidation model.
//! - **NoxModel**: Thermal-NOx and prompt-NOx formation.
//! - **CoModel**: CO formation/oxidation sub-mechanism.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights.
const W9: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 lattice velocities (cx, cy).
const C9: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

/// Speed of sound squared for D2Q9: cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn len2(a: [f64; 2]) -> f64 {
    dot2(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Arrhenius kinetics
// ---------------------------------------------------------------------------

/// Parameters for Arrhenius reaction rate: k(T) = A * exp(-Ea / (R * T)).
#[derive(Clone, Debug)]
pub struct ArrheniusKinetics {
    /// Pre-exponential factor \[1/s\].
    pub pre_exp: f64,
    /// Activation energy \[J/mol\].
    pub activation_energy: f64,
    /// Universal gas constant \[J/(mol·K)\].
    pub gas_constant: f64,
    /// Reaction order with respect to fuel.
    pub fuel_order: f64,
    /// Reaction order with respect to oxidizer.
    pub oxidizer_order: f64,
    /// Heat of combustion per unit mass of fuel \[J/kg\].
    pub heat_of_combustion: f64,
    /// Stoichiometric oxidizer-to-fuel mass ratio.
    pub stoich_ratio: f64,
}

impl ArrheniusKinetics {
    /// Create a new [`ArrheniusKinetics`] with given parameters.
    pub fn new(
        pre_exp: f64,
        activation_energy: f64,
        gas_constant: f64,
        fuel_order: f64,
        oxidizer_order: f64,
        heat_of_combustion: f64,
        stoich_ratio: f64,
    ) -> Self {
        Self {
            pre_exp,
            activation_energy,
            gas_constant,
            fuel_order,
            oxidizer_order,
            heat_of_combustion,
            stoich_ratio,
        }
    }

    /// Default methane-air kinetics (one-step global mechanism).
    pub fn methane_air() -> Self {
        Self::new(
            2.119e11, // pre-exponential
            2.027e8,  // Ea [J/mol] (approx 48.4 kcal/mol)
            8.314, 1.0, 1.0, 5.0e7, 4.0,
        )
    }

    /// Compute reaction rate \[kg/(m³·s)\] at temperature `t_kelvin`.
    ///
    /// Clamps negative concentrations to zero before evaluation.
    pub fn rate(&self, t_kelvin: f64, fuel: f64, oxidizer: f64) -> f64 {
        if t_kelvin <= 0.0 {
            return 0.0;
        }
        let k = self.pre_exp * (-self.activation_energy / (self.gas_constant * t_kelvin)).exp();
        let yf = fuel.max(0.0);
        let yo = oxidizer.max(0.0);
        k * yf.powf(self.fuel_order) * yo.powf(self.oxidizer_order)
    }

    /// Heat release rate \[W/m³\] at given state.
    pub fn heat_release_rate(&self, t_kelvin: f64, fuel: f64, oxidizer: f64) -> f64 {
        self.rate(t_kelvin, fuel, oxidizer) * self.heat_of_combustion
    }

    /// Zeldovich activation temperature Ta = Ea / R.
    pub fn activation_temperature(&self) -> f64 {
        self.activation_energy / self.gas_constant
    }
}

// ---------------------------------------------------------------------------
// Lewis number
// ---------------------------------------------------------------------------

/// Lewis number effects: ratio of thermal to mass diffusivity.
#[derive(Clone, Debug)]
pub struct LewisEffect {
    /// Lewis number Le = α / D (thermal diffusivity / mass diffusivity).
    pub lewis_number: f64,
    /// Reference thermal diffusivity \[m²/s\].
    pub thermal_diffusivity: f64,
}

impl LewisEffect {
    /// Create a new [`LewisEffect`].
    pub fn new(lewis_number: f64, thermal_diffusivity: f64) -> Self {
        Self {
            lewis_number,
            thermal_diffusivity,
        }
    }

    /// Effective mass diffusivity D = α / Le.
    pub fn mass_diffusivity(&self) -> f64 {
        self.thermal_diffusivity / self.lewis_number.max(1e-10)
    }

    /// Corrected species diffusivity omega for LBM (relaxation time shift).
    ///
    /// `omega_thermal` is the thermal relaxation parameter.
    pub fn species_omega(&self, omega_thermal: f64) -> f64 {
        // cs2_dt / (D + cs2_dt/2) where D = alpha/Le
        let _d = self.mass_diffusivity();
        let tau_thermal = 1.0 / omega_thermal;
        let d_thermal = CS2 * (tau_thermal - 0.5);
        let d_mass = d_thermal / self.lewis_number.max(1e-10);
        1.0 / (d_mass / CS2 + 0.5)
    }
}

// ---------------------------------------------------------------------------
// Species field
// ---------------------------------------------------------------------------

/// Multi-species field storing fuel (YF), oxidizer (YO), products (YP), and inerts (YN).
#[derive(Clone, Debug)]
pub struct SpeciesField {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Fuel mass fraction field, row-major \[ny\]\[nx\].
    pub fuel: Vec<f64>,
    /// Oxidizer mass fraction field.
    pub oxidizer: Vec<f64>,
    /// Products mass fraction field.
    pub products: Vec<f64>,
    /// Inert (e.g., N₂) mass fraction field.
    pub inert: Vec<f64>,
}

impl SpeciesField {
    /// Allocate a species field with uniform initial conditions.
    pub fn new(nx: usize, ny: usize, yf0: f64, yo0: f64, yp0: f64, yn0: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            fuel: vec![yf0; n],
            oxidizer: vec![yo0; n],
            products: vec![yp0; n],
            inert: vec![yn0; n],
        }
    }

    /// Linear index from (ix, iy).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Enforce mass-fraction sum = 1 at each cell.
    pub fn normalize(&mut self) {
        for i in 0..self.nx * self.ny {
            let s = self.fuel[i] + self.oxidizer[i] + self.products[i] + self.inert[i];
            if s > 1e-14 {
                let inv = 1.0 / s;
                self.fuel[i] *= inv;
                self.oxidizer[i] *= inv;
                self.products[i] *= inv;
                self.inert[i] *= inv;
            }
        }
    }

    /// Apply reaction source terms (explicit Euler).
    ///
    /// Updates YF, YO, YP in-place using Arrhenius kinetics and temperature field `temp`.
    pub fn react(&mut self, temp: &[f64], kinetics: &ArrheniusKinetics, dt: f64, _rho: &[f64]) {
        for i in 0..self.nx * self.ny {
            let t = temp[i];
            let yf = self.fuel[i];
            let yo = self.oxidizer[i];
            let r = kinetics.rate(t, yf, yo) * dt;
            let dyf = -r.min(yf);
            let dyo = -(kinetics.stoich_ratio * r).min(yo);
            let dyp = -dyf - dyo;
            self.fuel[i] += dyf;
            self.oxidizer[i] += dyo;
            self.products[i] += dyp;
        }
    }

    /// Diffuse all species using explicit finite differences with periodic BC.
    pub fn diffuse(&mut self, diffusivity: f64, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let d = diffusivity * dt;
        let mut fields = [
            self.fuel.clone(),
            self.oxidizer.clone(),
            self.products.clone(),
            self.inert.clone(),
        ];
        let src = fields.clone();
        for k in 0..4 {
            for iy in 0..ny {
                let iyp = (iy + 1) % ny;
                let iym = (iy + ny - 1) % ny;
                for ix in 0..nx {
                    let ixp = (ix + 1) % nx;
                    let ixm = (ix + nx - 1) % nx;
                    let lap = src[k][iy * nx + ixp]
                        + src[k][iy * nx + ixm]
                        + src[k][iyp * nx + ix]
                        + src[k][iym * nx + ix]
                        - 4.0 * src[k][iy * nx + ix];
                    fields[k][iy * nx + ix] = src[k][iy * nx + ix] + d * lap;
                }
            }
        }
        self.fuel = fields[0].clone();
        self.oxidizer = fields[1].clone();
        self.products = fields[2].clone();
        self.inert = fields[3].clone();
    }

    /// Mixture fraction Z = (s*YF - YO + YO_inf) / (s*YF_inf + YO_inf).
    ///
    /// Uses stoichiometric ratio `s` and boundary values `yf_inf`, `yo_inf`.
    pub fn mixture_fraction(&self, s: f64, yf_inf: f64, yo_inf: f64, ix: usize, iy: usize) -> f64 {
        let i = self.idx(ix, iy);
        let num = s * self.fuel[i] - self.oxidizer[i] + yo_inf;
        let den = s * yf_inf + yo_inf;
        if den.abs() < 1e-14 { 0.0 } else { num / den }
    }
}

// ---------------------------------------------------------------------------
// Temperature field with heat release
// ---------------------------------------------------------------------------

/// Temperature field with diffusion and heat release from combustion.
#[derive(Clone, Debug)]
pub struct ThermalField {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Temperature values \[K\].
    pub temp: Vec<f64>,
    /// Thermal diffusivity \[lattice units\].
    pub alpha: f64,
    /// Specific heat capacity \[J/(kg·K)\].
    pub cp: f64,
    /// Reference density \[kg/m³\].
    pub rho_ref: f64,
}

impl ThermalField {
    /// Create a new [`ThermalField`] with uniform temperature `t0`.
    pub fn new(nx: usize, ny: usize, t0: f64, alpha: f64, cp: f64, rho_ref: f64) -> Self {
        Self {
            nx,
            ny,
            temp: vec![t0; nx * ny],
            alpha,
            cp,
            rho_ref,
        }
    }

    /// Linear index.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Apply heat release source term (explicit Euler).
    pub fn apply_heat_release(
        &mut self,
        species: &SpeciesField,
        kinetics: &ArrheniusKinetics,
        dt: f64,
    ) {
        let denom = self.cp * self.rho_ref;
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let i = self.idx(ix, iy);
                let t = self.temp[i];
                let yf = species.fuel[i];
                let yo = species.oxidizer[i];
                let q = kinetics.heat_release_rate(t, yf, yo);
                self.temp[i] += q * dt / denom;
            }
        }
    }

    /// Diffuse temperature with periodic BC.
    pub fn diffuse(&mut self, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let d = self.alpha * dt;
        let old = self.temp.clone();
        for iy in 0..ny {
            let iyp = (iy + 1) % ny;
            let iym = (iy + ny - 1) % ny;
            for ix in 0..nx {
                let ixp = (ix + 1) % nx;
                let ixm = (ix + nx - 1) % nx;
                let lap = old[iy * nx + ixp]
                    + old[iy * nx + ixm]
                    + old[iyp * nx + ix]
                    + old[iym * nx + ix]
                    - 4.0 * old[iy * nx + ix];
                self.temp[iy * nx + ix] = old[iy * nx + ix] + d * lap;
            }
        }
    }

    /// Maximum temperature in the domain.
    pub fn max_temp(&self) -> f64 {
        self.temp.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Adiabatic flame temperature estimate: T_ad = T0 + q/(cp*rho).
    pub fn adiabatic_flame_temp(t0: f64, heat_of_combustion: f64, yf0: f64, cp: f64) -> f64 {
        t0 + heat_of_combustion * yf0 / cp
    }
}

// ---------------------------------------------------------------------------
// Flame front tracker
// ---------------------------------------------------------------------------

/// Tracks the flame front as the iso-surface where reaction rate is maximal.
#[derive(Clone, Debug)]
pub struct FlameTracker {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Temperature threshold for flame front detection \[K\].
    pub t_threshold: f64,
    /// Stores the x-position of the flame front for each row iy.
    pub front_x: Vec<f64>,
}

impl FlameTracker {
    /// Create a new [`FlameTracker`].
    pub fn new(nx: usize, ny: usize, t_threshold: f64) -> Self {
        Self {
            nx,
            ny,
            t_threshold,
            front_x: vec![0.0; ny],
        }
    }

    /// Locate the flame front by finding the first x-cell exceeding `t_threshold`.
    pub fn locate_front(&mut self, temp: &ThermalField) {
        for iy in 0..self.ny {
            let mut found = false;
            for ix in 0..self.nx {
                let i = temp.idx(ix, iy);
                if temp.temp[i] >= self.t_threshold {
                    self.front_x[iy] = ix as f64;
                    found = true;
                    break;
                }
            }
            if !found {
                self.front_x[iy] = self.nx as f64;
            }
        }
    }

    /// Mean flame front x-position averaged over all rows.
    pub fn mean_front_x(&self) -> f64 {
        let s: f64 = self.front_x.iter().sum();
        s / self.ny as f64
    }

    /// Flame front propagation speed estimated from two frames.
    pub fn flame_speed(&self, prev_front_x: f64, dt: f64) -> f64 {
        let curr = self.mean_front_x();
        (curr - prev_front_x) / dt.max(1e-14)
    }

    /// Wrinkle factor: ratio of actual front length to domain width.
    pub fn wrinkle_factor(&self) -> f64 {
        let mean = self.mean_front_x();
        let variance: f64 = self
            .front_x
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / self.ny as f64;
        let sigma = variance.sqrt();
        1.0 + sigma / (self.nx as f64).max(1.0)
    }
}

// ---------------------------------------------------------------------------
// Deflagration-to-detonation transition (DDT)
// ---------------------------------------------------------------------------

/// DDT model: monitors local pressure and reaction rate to detect onset.
#[derive(Clone, Debug)]
pub struct DdtModel {
    /// Critical pressure ratio p/p0 for DDT onset.
    pub p_crit_ratio: f64,
    /// Critical reaction rate (normalized) for DDT onset.
    pub rate_crit: f64,
    /// Whether DDT has been detected.
    pub ddt_detected: bool,
    /// Cell index (flat) where DDT was first detected.
    pub ddt_location: usize,
    /// Detonation wave speed \[m/s\] (Chapman-Jouguet approximation).
    pub cj_speed: f64,
    /// Current run-up distance \[m\].
    pub run_up_distance: f64,
}

impl DdtModel {
    /// Create a new [`DdtModel`].
    pub fn new(p_crit_ratio: f64, rate_crit: f64, cj_speed: f64) -> Self {
        Self {
            p_crit_ratio,
            rate_crit,
            ddt_detected: false,
            ddt_location: 0,
            cj_speed,
            run_up_distance: 0.0,
        }
    }

    /// Check if DDT has occurred given local pressure and reaction rate fields.
    pub fn check_ddt(
        &mut self,
        pressure: &[f64],
        p_ref: f64,
        reaction_rate: &[f64],
        dx: f64,
    ) -> bool {
        if self.ddt_detected {
            return true;
        }
        for (i, (&p, &r)) in pressure.iter().zip(reaction_rate.iter()).enumerate() {
            if p / p_ref.max(1e-14) >= self.p_crit_ratio && r >= self.rate_crit {
                self.ddt_detected = true;
                self.ddt_location = i;
                self.run_up_distance = i as f64 * dx;
                return true;
            }
        }
        false
    }

    /// Chapman-Jouguet Mach number given speed of sound `c0`.
    pub fn cj_mach(&self, c0: f64) -> f64 {
        self.cj_speed / c0.max(1e-14)
    }

    /// Reset DDT state.
    pub fn reset(&mut self) {
        self.ddt_detected = false;
        self.ddt_location = 0;
        self.run_up_distance = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Flame speed computation
// ---------------------------------------------------------------------------

/// Flame speed computation: laminar and turbulent correlations.
#[derive(Clone, Debug)]
pub struct FlameSpeed {
    /// Laminar flame speed SL \[m/s\].
    pub sl: f64,
    /// Thermal diffusivity α \[m²/s\].
    pub alpha: f64,
    /// Reaction zone thickness δ_L \[m\].
    pub delta_l: f64,
    /// Turbulent integral length scale l_t \[m\].
    pub l_t: f64,
    /// Turbulent rms velocity u' \[m/s\].
    pub u_rms: f64,
}

impl FlameSpeed {
    /// Create a new [`FlameSpeed`].
    pub fn new(sl: f64, alpha: f64, l_t: f64, u_rms: f64) -> Self {
        let delta_l = alpha / sl.max(1e-14);
        Self {
            sl,
            alpha,
            delta_l,
            l_t,
            u_rms,
        }
    }

    /// Turbulent flame speed using Damköhler correlation: ST = SL + u'.
    pub fn turbulent_damkohler(&self) -> f64 {
        self.sl + self.u_rms
    }

    /// Turbulent flame speed using Peters correlation.
    ///
    /// ST/SL = 1 + C * (u'/SL)^0.5 * (l_t/δ_L)^0.25
    pub fn turbulent_peters(&self, c: f64) -> f64 {
        let ratio_u = (self.u_rms / self.sl.max(1e-14)).powf(0.5);
        let ratio_l = (self.l_t / self.delta_l.max(1e-14)).powf(0.25);
        self.sl * (1.0 + c * ratio_u * ratio_l)
    }

    /// Damköhler number Da = (l_t * SL) / (u' * δ_L).
    pub fn damkohler_number(&self) -> f64 {
        (self.l_t * self.sl) / ((self.u_rms * self.delta_l).max(1e-14))
    }

    /// Karlovitz number Ka = (u'/SL)^2 * (δ_L/l_t)^0.5.
    pub fn karlovitz_number(&self) -> f64 {
        let ratio_u = (self.u_rms / self.sl.max(1e-14)).powi(2);
        let ratio_l = (self.delta_l / self.l_t.max(1e-14)).powf(0.5);
        ratio_u * ratio_l
    }

    /// Borghi-Peters diagram regime classification.
    ///
    /// Returns a string describing the combustion regime.
    pub fn borghi_regime(&self) -> &'static str {
        let da = self.damkohler_number();
        let ka = self.karlovitz_number();
        let u_ratio = self.u_rms / self.sl.max(1e-14);
        if u_ratio < 1.0 {
            "laminar"
        } else if ka < 1.0 {
            "wrinkled-flamelets"
        } else if da > 1.0 {
            "corrugated-flamelets"
        } else {
            "thin-reaction-zones"
        }
    }
}

// ---------------------------------------------------------------------------
// Soot model
// ---------------------------------------------------------------------------

/// Simple soot formation and oxidation model (phenomenological).
#[derive(Clone, Debug)]
pub struct SootModel {
    /// Number of cells.
    pub n: usize,
    /// Soot volume fraction field fv \[-\].
    pub fv: Vec<f64>,
    /// Soot nucleation rate coefficient \[1/s\].
    pub k_nuc: f64,
    /// Soot surface growth rate coefficient \[m/s\].
    pub k_sg: f64,
    /// Soot oxidation rate coefficient \[1/s\].
    pub k_ox: f64,
    /// Soot particle density \[kg/m³\].
    pub rho_soot: f64,
}

impl SootModel {
    /// Create a new [`SootModel`] with zero initial soot.
    pub fn new(n: usize, k_nuc: f64, k_sg: f64, k_ox: f64, rho_soot: f64) -> Self {
        Self {
            n,
            fv: vec![0.0; n],
            k_nuc,
            k_sg,
            k_ox,
            rho_soot,
        }
    }

    /// Advance soot model one time step.
    ///
    /// `temp` is temperature \[K\], `yf` is fuel mass fraction, `yo2` is O₂ mass fraction.
    pub fn step(&mut self, temp: &[f64], yf: &[f64], yo2: &[f64], dt: f64) {
        for i in 0..self.n {
            let t = temp[i];
            let f = yf[i].max(0.0);
            let o = yo2[i].max(0.0);
            // Nucleation (requires high temperature and fuel)
            let nuc = if t > 1400.0 { self.k_nuc * f } else { 0.0 };
            // Surface growth
            let sg = self.k_sg * f * self.fv[i].sqrt();
            // Oxidation
            let ox = self.k_ox * o * self.fv[i];
            self.fv[i] += (nuc + sg - ox) * dt;
            self.fv[i] = self.fv[i].max(0.0);
        }
    }

    /// Total soot mass per unit volume (soot volume fraction × particle density).
    pub fn soot_mass_concentration(&self) -> Vec<f64> {
        self.fv.iter().map(|&f| f * self.rho_soot).collect()
    }

    /// Mean soot volume fraction over domain.
    pub fn mean_fv(&self) -> f64 {
        self.fv.iter().sum::<f64>() / self.n as f64
    }
}

// ---------------------------------------------------------------------------
// NOx model
// ---------------------------------------------------------------------------

/// Thermal-NOx formation model (extended Zeldovich mechanism).
#[derive(Clone, Debug)]
pub struct NoxModel {
    /// Number of cells.
    pub n: usize,
    /// NO mass fraction field.
    pub no: Vec<f64>,
    /// NO₂ mass fraction field.
    pub no2: Vec<f64>,
    /// Thermal NOx rate coefficient A₁ \[m³/(mol·s)\].
    pub k_thermal: f64,
    /// Activation energy for thermal NOx \[K\] (= Ea/R).
    pub ta_thermal: f64,
    /// Prompt NOx coefficient.
    pub k_prompt: f64,
}

impl NoxModel {
    /// Create a new [`NoxModel`] with zero initial NOx.
    pub fn new(n: usize, k_thermal: f64, ta_thermal: f64, k_prompt: f64) -> Self {
        Self {
            n,
            no: vec![0.0; n],
            no2: vec![0.0; n],
            k_thermal,
            ta_thermal,
            k_prompt,
        }
    }

    /// Advance thermal-NOx formation one time step.
    ///
    /// `temp` \[K\], `yn2` N₂ mass fraction, `yo2` O₂ mass fraction.
    pub fn step_thermal(&mut self, temp: &[f64], yn2: &[f64], yo2: &[f64], dt: f64) {
        for i in 0..self.n {
            let t = temp[i];
            if t < 1500.0 {
                continue;
            }
            let k = self.k_thermal * (-self.ta_thermal / t).exp();
            let d_no = k * yo2[i].max(0.0).sqrt() * yn2[i].max(0.0) * dt;
            self.no[i] += d_no;
        }
    }

    /// Advance prompt-NOx formation one time step.
    ///
    /// `yf` fuel mass fraction, `yo2` O₂ mass fraction.
    pub fn step_prompt(&mut self, temp: &[f64], yf: &[f64], yo2: &[f64], dt: f64) {
        for i in 0..self.n {
            let t = temp[i];
            if t < 800.0 {
                continue;
            }
            let k = self.k_prompt * (-38000.0 / t).exp();
            let d_no = k * yf[i].max(0.0) * yo2[i].max(0.0) * dt;
            self.no[i] += d_no;
        }
    }

    /// Total NOx concentration (NO + NO₂) at cell `i`.
    pub fn total_nox(&self, i: usize) -> f64 {
        self.no[i] + self.no2[i]
    }

    /// Mean NO over domain.
    pub fn mean_no(&self) -> f64 {
        self.no.iter().sum::<f64>() / self.n as f64
    }
}

// ---------------------------------------------------------------------------
// CO model
// ---------------------------------------------------------------------------

/// CO formation/oxidation sub-mechanism.
#[derive(Clone, Debug)]
pub struct CoModel {
    /// Number of cells.
    pub n: usize,
    /// CO mass fraction field.
    pub co: Vec<f64>,
    /// CO formation rate coefficient \[1/s\].
    pub k_form: f64,
    /// CO oxidation rate coefficient \[1/s\].
    pub k_ox: f64,
    /// CO oxidation activation temperature \[K\].
    pub ta_ox: f64,
}

impl CoModel {
    /// Create a new [`CoModel`] with zero initial CO.
    pub fn new(n: usize, k_form: f64, k_ox: f64, ta_ox: f64) -> Self {
        Self {
            n,
            co: vec![0.0; n],
            k_form,
            k_ox,
            ta_ox,
        }
    }

    /// Advance CO model one time step.
    ///
    /// `temp` \[K\], `yf` fuel fraction, `yo2` O₂ fraction.
    pub fn step(&mut self, temp: &[f64], yf: &[f64], yo2: &[f64], dt: f64) {
        for i in 0..self.n {
            let t = temp[i];
            let f = yf[i].max(0.0);
            let o = yo2[i].max(0.0);
            let form = self.k_form * f;
            let ox = self.k_ox * (-self.ta_ox / t.max(300.0)).exp() * o * self.co[i];
            self.co[i] += (form - ox) * dt;
            self.co[i] = self.co[i].max(0.0);
        }
    }

    /// Mean CO over domain.
    pub fn mean_co(&self) -> f64 {
        self.co.iter().sum::<f64>() / self.n as f64
    }
}

// ---------------------------------------------------------------------------
// LBM distribution function utilities
// ---------------------------------------------------------------------------

/// Compute D2Q9 equilibrium distribution for given macroscopic state.
///
/// `rho` is density, `ux`/`uy` are velocity components.
pub fn equilibrium_d2q9(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; 9];
    for q in 0..9 {
        let cu = C9[q].0 as f64 * ux + C9[q].1 as f64 * uy;
        feq[q] = W9[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

/// Compute D2Q9 equilibrium for a passive scalar `phi` advected by velocity (ux, uy).
pub fn scalar_equilibrium_d2q9(phi: f64, ux: f64, uy: f64, omega_s: f64) -> [f64; 9] {
    let _ = omega_s;
    let mut geq = [0.0f64; 9];
    for q in 0..9 {
        let cu = C9[q].0 as f64 * ux + C9[q].1 as f64 * uy;
        geq[q] = W9[q] * phi * (1.0 + cu / CS2);
    }
    geq
}

/// BGK collision step for a single cell distribution.
pub fn bgk_collision(f: [f64; 9], feq: [f64; 9], omega: f64) -> [f64; 9] {
    let mut f_out = [0.0f64; 9];
    for q in 0..9 {
        f_out[q] = f[q] - omega * (f[q] - feq[q]);
    }
    f_out
}

// ---------------------------------------------------------------------------
// Main combustion LBM solver
// ---------------------------------------------------------------------------

/// Main combustion LBM solver coupling flow, species, and temperature.
///
/// Uses D2Q9 lattice for flow, with separate passive-scalar LBM for each species.
#[allow(clippy::too_many_arguments)]
pub struct CombustionLbm {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Flow distribution functions f\[q\]\[iy*nx+ix\].
    pub f: Vec<[f64; 9]>,
    /// Species field.
    pub species: SpeciesField,
    /// Temperature field.
    pub temperature: ThermalField,
    /// Arrhenius kinetics.
    pub kinetics: ArrheniusKinetics,
    /// Flow relaxation parameter.
    pub omega: f64,
    /// Species diffusion relaxation parameter.
    pub omega_s: f64,
    /// Flame tracker.
    pub tracker: FlameTracker,
    /// Soot model.
    pub soot: SootModel,
    /// NOx model.
    pub nox: NoxModel,
    /// CO model.
    pub co_model: CoModel,
    /// Time step counter.
    pub step: usize,
    /// Physical time.
    pub time: f64,
    /// Physical time step size.
    pub dt: f64,
}

impl CombustionLbm {
    /// Create a new [`CombustionLbm`] solver.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        omega_s: f64,
        kinetics: ArrheniusKinetics,
        t0: f64,
        dt: f64,
    ) -> Self {
        let n = nx * ny;
        let alpha = CS2 * (1.0 / omega - 0.5);
        let nu = alpha;
        let _ = nu;

        // Initialise f to rest state (rho=1, u=0)
        let feq0 = equilibrium_d2q9(1.0, 0.0, 0.0);
        let f = vec![feq0; n];

        let species = SpeciesField::new(nx, ny, 0.05, 0.23, 0.0, 0.72);
        let temperature = ThermalField::new(nx, ny, t0, alpha, 1200.0, 1.2);
        let tracker = FlameTracker::new(nx, ny, t0 + 200.0);
        let soot = SootModel::new(n, 1e-4, 1e-5, 1e-3, 1800.0);
        let nox = NoxModel::new(n, 6e16, 69090.0, 1e10);
        let co_model = CoModel::new(n, 0.1, 1e5, 20000.0);

        Self {
            nx,
            ny,
            f,
            species,
            temperature,
            kinetics,
            omega,
            omega_s,
            tracker,
            soot,
            nox,
            co_model,
            step: 0,
            time: 0.0,
            dt,
        }
    }

    /// Set a hot ignition spot at cell (ix, iy) with temperature boost `delta_t`.
    pub fn ignite(&mut self, ix: usize, iy: usize, delta_t: f64) {
        let i = self.temperature.idx(ix, iy);
        self.temperature.temp[i] += delta_t;
    }

    /// Macroscopic density at cell index `i`.
    pub fn density(&self, i: usize) -> f64 {
        self.f[i].iter().sum()
    }

    /// Macroscopic velocity at cell index `i`.
    pub fn velocity(&self, i: usize) -> [f64; 2] {
        let rho = self.density(i);
        let ux: f64 = self.f[i]
            .iter()
            .enumerate()
            .map(|(q, &fq)| C9[q].0 as f64 * fq)
            .sum::<f64>()
            / rho.max(1e-14);
        let uy: f64 = self.f[i]
            .iter()
            .enumerate()
            .map(|(q, &fq)| C9[q].1 as f64 * fq)
            .sum::<f64>()
            / rho.max(1e-14);
        [ux, uy]
    }

    /// Perform one LBM time step: collision, streaming, species, temperature, pollutants.
    pub fn advance(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // --- Collision ---
        let mut f_post = self.f.clone();
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                let rho = self.density(i);
                let [ux, uy] = self.velocity(i);
                let feq = equilibrium_d2q9(rho, ux, uy);
                f_post[i] = bgk_collision(self.f[i], feq, self.omega);
            }
        }

        // --- Streaming ---
        let mut f_new = self.f.clone();
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                for q in 0..9 {
                    let ixs = ((ix as i32 - C9[q].0).rem_euclid(nx as i32)) as usize;
                    let iys = ((iy as i32 - C9[q].1).rem_euclid(ny as i32)) as usize;
                    f_new[i][q] = f_post[iys * nx + ixs][q];
                }
            }
        }
        self.f = f_new;

        // --- Species reaction and diffusion ---
        let diffusivity = CS2 * (1.0 / self.omega_s - 0.5);
        let rho_field: Vec<f64> = (0..nx * ny).map(|i| self.density(i)).collect();
        self.species.react(
            &self.temperature.temp.clone(),
            &self.kinetics,
            self.dt,
            &rho_field,
        );
        self.species.diffuse(diffusivity, self.dt);

        // --- Temperature: heat release + diffusion ---
        self.temperature
            .apply_heat_release(&self.species, &self.kinetics, self.dt);
        self.temperature.diffuse(self.dt);

        // --- Pollutant sub-models ---
        let yf = self.species.fuel.clone();
        let yo = self.species.oxidizer.clone();
        let temp = self.temperature.temp.clone();
        // Approximate N2 as (1 - YF - YO - YP)
        let yn2: Vec<f64> = (0..nx * ny)
            .map(|i| (1.0 - yf[i] - yo[i] - self.species.products[i]).max(0.0))
            .collect();

        self.soot.step(&temp, &yf, &yo, self.dt);
        self.nox.step_thermal(&temp, &yn2, &yo, self.dt);
        self.nox.step_prompt(&temp, &yf, &yo, self.dt);
        self.co_model.step(&temp, &yf, &yo, self.dt);

        // --- Flame tracking ---
        self.tracker.locate_front(&self.temperature);

        self.step += 1;
        self.time += self.dt;
    }

    /// Run `n_steps` time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.advance();
        }
    }

    /// Return the current mean flame front x-position.
    pub fn mean_flame_x(&self) -> f64 {
        self.tracker.mean_front_x()
    }

    /// Return mean CO concentration.
    pub fn mean_co(&self) -> f64 {
        self.co_model.mean_co()
    }

    /// Return mean NO concentration.
    pub fn mean_no(&self) -> f64 {
        self.nox.mean_no()
    }

    /// Return mean soot volume fraction.
    pub fn mean_soot(&self) -> f64 {
        self.soot.mean_fv()
    }
}

// ---------------------------------------------------------------------------
// Premixed flame configuration helper
// ---------------------------------------------------------------------------

/// Configuration helper for a 1D premixed flat flame.
#[derive(Clone, Debug)]
pub struct PremixedFlameConfig {
    /// Equivalence ratio φ = (F/O)_actual / (F/O)_stoich.
    pub equivalence_ratio: f64,
    /// Unburnt temperature \[K\].
    pub t_unburnt: f64,
    /// Ambient pressure \[Pa\].
    pub pressure: f64,
    /// Fuel type identifier.
    pub fuel_type: FuelType,
}

/// Supported fuel types for premixed flame configuration.
#[derive(Clone, Debug)]
pub enum FuelType {
    /// Methane CH₄.
    Methane,
    /// Hydrogen H₂.
    Hydrogen,
    /// Propane C₃H₈.
    Propane,
    /// Generic hydrocarbon fuel.
    Generic,
}

impl PremixedFlameConfig {
    /// Create a stoichiometric methane-air premixed flame at 300 K, 1 atm.
    pub fn methane_stoich() -> Self {
        Self {
            equivalence_ratio: 1.0,
            t_unburnt: 300.0,
            pressure: 101325.0,
            fuel_type: FuelType::Methane,
        }
    }

    /// Adiabatic flame temperature estimate using simplified enthalpy balance.
    pub fn adiabatic_temp(&self) -> f64 {
        let base = match self.fuel_type {
            FuelType::Methane => 2226.0,
            FuelType::Hydrogen => 2480.0,
            FuelType::Propane => 2267.0,
            FuelType::Generic => 2100.0,
        };
        // Lean/rich correction
        let phi = self.equivalence_ratio.min(1.0);
        let dt = self.t_unburnt - 300.0;
        base * phi + dt
    }

    /// Stoichiometric fuel mass fraction.
    pub fn stoich_fuel_fraction(&self) -> f64 {
        match self.fuel_type {
            FuelType::Methane => 0.055,
            FuelType::Hydrogen => 0.029,
            FuelType::Propane => 0.060,
            FuelType::Generic => 0.055,
        }
    }

    /// Laminar flame speed correlation \[m/s\] for methane-air at 300 K, 1 atm.
    pub fn laminar_flame_speed(&self) -> f64 {
        // Gülder correlation for methane: SL0 * (phi/phi_m)^eta * exp(-xi*(phi-phi_m)^2)
        let sl0 = match self.fuel_type {
            FuelType::Methane => 0.37,
            FuelType::Hydrogen => 2.38,
            FuelType::Propane => 0.43,
            FuelType::Generic => 0.35,
        };
        let phi = self.equivalence_ratio;
        let phi_m = 1.05; // phi at max SL
        sl0 * (phi / phi_m).powf(0.5) * (-2.0 * (phi - phi_m).powi(2)).exp()
    }
}

// ---------------------------------------------------------------------------
// Diffusion flame configuration helper
// ---------------------------------------------------------------------------

/// Configuration for a diffusion (non-premixed) flame.
#[derive(Clone, Debug)]
pub struct DiffusionFlameConfig {
    /// Fuel stream velocity \[m/s\].
    pub v_fuel: f64,
    /// Oxidizer stream velocity \[m/s\].
    pub v_ox: f64,
    /// Fuel stream temperature \[K\].
    pub t_fuel: f64,
    /// Oxidizer stream temperature \[K\].
    pub t_ox: f64,
    /// Stoichiometric mixture fraction Zst.
    pub z_st: f64,
}

impl DiffusionFlameConfig {
    /// Create a default methane-air Burke-Schumann diffusion flame.
    pub fn methane_air_default() -> Self {
        Self {
            v_fuel: 0.1,
            v_ox: 0.1,
            t_fuel: 300.0,
            t_ox: 300.0,
            z_st: 0.055,
        }
    }

    /// Burke-Schumann adiabatic flame temperature at stoichiometry.
    pub fn bs_flame_temp(&self, q_comb: f64, cp: f64) -> f64 {
        let t_mix = self.z_st * self.t_fuel + (1.0 - self.z_st) * self.t_ox;
        t_mix + q_comb * self.z_st / cp
    }

    /// Mixture fraction profile for co-flow geometry at x-position `x` and radius `r`.
    ///
    /// Uses simplified Gaussian profile: Z(r) = Zst * exp(-r^2 / (2*sigma^2)).
    pub fn mixture_fraction_profile(&self, r: f64, sigma: f64) -> f64 {
        self.z_st * (-r * r / (2.0 * sigma * sigma)).exp()
    }
}

// ---------------------------------------------------------------------------
// Flame extinction model
// ---------------------------------------------------------------------------

/// Flame extinction model based on Damköhler number criterion.
#[derive(Clone, Debug)]
pub struct ExtinctionModel {
    /// Critical Damköhler number below which extinction occurs.
    pub da_crit: f64,
    /// Scalar dissipation rate χ \[1/s\].
    pub chi: f64,
    /// Chemical time scale τ_chem \[s\].
    pub tau_chem: f64,
}

impl ExtinctionModel {
    /// Create a new [`ExtinctionModel`].
    pub fn new(da_crit: f64, chi: f64, tau_chem: f64) -> Self {
        Self {
            da_crit,
            chi,
            tau_chem,
        }
    }

    /// Damköhler number Da = 1 / (χ * τ_chem).
    pub fn damkohler_number(&self) -> f64 {
        1.0 / (self.chi * self.tau_chem).max(1e-14)
    }

    /// Returns true if the flame is extinguished.
    pub fn is_extinguished(&self) -> bool {
        self.damkohler_number() < self.da_crit
    }

    /// Critical scalar dissipation rate χ_crit = 1 / (Da_crit * τ_chem).
    pub fn chi_crit(&self) -> f64 {
        1.0 / (self.da_crit * self.tau_chem).max(1e-14)
    }
}

// ---------------------------------------------------------------------------
// Scalar dissipation rate
// ---------------------------------------------------------------------------

/// Compute scalar dissipation rate χ = 2D |∇Z|² at a cell.
///
/// `grad_z` is the mixture fraction gradient vector \[dZ/dx, dZ/dy\].
pub fn scalar_dissipation(diffusivity: f64, grad_z: [f64; 2]) -> f64 {
    2.0 * diffusivity * (grad_z[0] * grad_z[0] + grad_z[1] * grad_z[1])
}

/// Compute mixture fraction gradient at cell (ix, iy) using central differences.
pub fn mixture_fraction_gradient(
    z: &[f64],
    nx: usize,
    ny: usize,
    ix: usize,
    iy: usize,
    dx: f64,
    dy: f64,
) -> [f64; 2] {
    let ixp = (ix + 1).min(nx - 1);
    let ixm = ix.saturating_sub(1);
    let iyp = (iy + 1).min(ny - 1);
    let iym = iy.saturating_sub(1);
    let dzdx = (z[iy * nx + ixp] - z[iy * nx + ixm]) / (2.0 * dx);
    let dzdy = (z[iyp * nx + ix] - z[iym * nx + ix]) / (2.0 * dy);
    [dzdx, dzdy]
}

// ---------------------------------------------------------------------------
// Flame surface density model
// ---------------------------------------------------------------------------

/// Coherent flame model: flame surface density Σ \[1/m\].
#[derive(Clone, Debug)]
pub struct FlameSurfaceDensity {
    /// Number of cells.
    pub n: usize,
    /// Flame surface density field Σ.
    pub sigma: Vec<f64>,
    /// Flame surface area production rate.
    pub alpha_prod: f64,
    /// Flame surface area destruction rate.
    pub beta_dest: f64,
}

impl FlameSurfaceDensity {
    /// Create a new [`FlameSurfaceDensity`] model.
    pub fn new(n: usize, alpha_prod: f64, beta_dest: f64) -> Self {
        Self {
            n,
            sigma: vec![0.0; n],
            alpha_prod,
            beta_dest,
        }
    }

    /// Advance flame surface density one step using the CFM transport equation.
    ///
    /// `st` is turbulent flame speed, `sl` is laminar flame speed.
    pub fn step(&mut self, st: f64, sl: f64, dt: f64) {
        let prod = self.alpha_prod * st;
        let dest = self.beta_dest * sl;
        for s in self.sigma.iter_mut() {
            *s += (prod * *s - dest * (*s) * (*s)) * dt;
            *s = s.max(0.0);
        }
    }

    /// Reaction rate per unit volume from FSD: ω = ρ * SL * Σ.
    pub fn reaction_rate(&self, rho: f64, sl: f64) -> Vec<f64> {
        self.sigma.iter().map(|&s| rho * sl * s).collect()
    }
}

// ---------------------------------------------------------------------------
// Thermo-diffusive instability
// ---------------------------------------------------------------------------

/// Thermo-diffusive instability analysis for premixed flames.
///
/// Instability occurs when Le < 1 for cellular flames or Le > 1 for diffusive-thermal.
#[derive(Clone, Debug)]
pub struct ThermoDiffusiveInstability {
    /// Lewis number of deficient reactant.
    pub lewis_number: f64,
    /// Zeldovich number Ze = Ta * (Tb - Tu) / Tb².
    pub zeldovich_number: f64,
}

impl ThermoDiffusiveInstability {
    /// Create a new [`ThermoDiffusiveInstability`] analysis.
    pub fn new(lewis_number: f64, zeldovich_number: f64) -> Self {
        Self {
            lewis_number,
            zeldovich_number,
        }
    }

    /// Markstein length L = δ_L * Ma where Ma is Markstein number.
    ///
    /// Uses approximation Ma ≈ Ze * (Le - 1) / 2.
    pub fn markstein_length(&self, delta_l: f64) -> f64 {
        let ma = self.zeldovich_number * (self.lewis_number - 1.0) / 2.0;
        delta_l * ma
    }

    /// Check if cellular instability is expected (Le < critical Le).
    pub fn is_cellular(&self) -> bool {
        self.lewis_number < 1.0 - 2.0 / self.zeldovich_number.max(1.0)
    }

    /// Onset wavenumber for thermo-diffusive instability.
    pub fn onset_wavenumber(&self, delta_l: f64) -> f64 {
        if self.is_cellular() {
            1.0 / delta_l.max(1e-14)
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arrhenius_zero_temp() {
        let k = ArrheniusKinetics::methane_air();
        assert_eq!(k.rate(0.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_arrhenius_rate_positive() {
        let k = ArrheniusKinetics::methane_air();
        let r = k.rate(2000.0, 0.05, 0.23);
        assert!(r >= 0.0);
    }

    #[test]
    fn test_arrhenius_heat_release() {
        let k = ArrheniusKinetics::methane_air();
        let q = k.heat_release_rate(2000.0, 0.05, 0.23);
        assert!(q >= 0.0);
    }

    #[test]
    fn test_arrhenius_activation_temp() {
        let k = ArrheniusKinetics::new(1e10, 8314.0, 8.314, 1.0, 1.0, 1e7, 4.0);
        let ta = k.activation_temperature();
        assert!((ta - 1000.0).abs() < 1.0);
    }

    #[test]
    fn test_lewis_mass_diffusivity() {
        let le = LewisEffect::new(1.0, 2e-5);
        assert!((le.mass_diffusivity() - 2e-5).abs() < 1e-12);
    }

    #[test]
    fn test_lewis_effect_le_half() {
        let le = LewisEffect::new(0.5, 2e-5);
        assert!((le.mass_diffusivity() - 4e-5).abs() < 1e-12);
    }

    #[test]
    fn test_species_field_normalization() {
        let mut sf = SpeciesField::new(4, 4, 0.3, 0.3, 0.3, 0.1);
        sf.normalize();
        let i = sf.idx(1, 1);
        let sum = sf.fuel[i] + sf.oxidizer[i] + sf.products[i] + sf.inert[i];
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_species_field_idx() {
        let sf = SpeciesField::new(10, 8, 0.05, 0.23, 0.0, 0.72);
        assert_eq!(sf.idx(3, 2), 2 * 10 + 3);
    }

    #[test]
    fn test_species_diffuse_conserves_mass() {
        let mut sf = SpeciesField::new(8, 8, 0.05, 0.23, 0.0, 0.72);
        let initial_sum: f64 = sf.fuel.iter().sum();
        sf.diffuse(0.1, 0.1);
        let final_sum: f64 = sf.fuel.iter().sum();
        assert!((initial_sum - final_sum).abs() < 1e-10);
    }

    #[test]
    fn test_species_mixture_fraction_stoich() {
        let sf = SpeciesField::new(4, 4, 0.055, 0.0, 0.0, 0.945);
        // At pure fuel, Z should be ~1
        let z = sf.mixture_fraction(4.0, 0.055, 0.0, 0, 0);
        assert!((z - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_field_diffuse() {
        let mut tf = ThermalField::new(8, 8, 300.0, 0.1, 1200.0, 1.2);
        tf.temp[32] = 2000.0;
        tf.diffuse(1.0);
        // Peak should decrease after diffusion
        assert!(tf.temp[32] < 2000.0);
    }

    #[test]
    fn test_thermal_max_temp() {
        let mut tf = ThermalField::new(4, 4, 300.0, 0.1, 1200.0, 1.2);
        tf.temp[5] = 3000.0;
        assert_eq!(tf.max_temp(), 3000.0);
    }

    #[test]
    fn test_adiabatic_flame_temp() {
        let t_ad = ThermalField::adiabatic_flame_temp(300.0, 5e7, 0.055, 1200.0);
        assert!(t_ad > 300.0);
    }

    #[test]
    fn test_flame_tracker_locate() {
        let mut tf = ThermalField::new(10, 4, 300.0, 0.1, 1200.0, 1.2);
        let mut tracker = FlameTracker::new(10, 4, 500.0);
        // Set row 0, cell ix=5 above threshold
        tf.temp[5] = 600.0;
        tracker.locate_front(&tf);
        assert_eq!(tracker.front_x[0], 5.0);
    }

    #[test]
    fn test_flame_tracker_mean() {
        let mut tracker = FlameTracker::new(10, 4, 500.0);
        tracker.front_x = vec![2.0, 4.0, 6.0, 8.0];
        assert_eq!(tracker.mean_front_x(), 5.0);
    }

    #[test]
    fn test_flame_speed() {
        let tracker = FlameTracker::new(10, 4, 500.0);
        let speed = tracker.flame_speed(3.0, 1.0);
        assert!((speed - tracker.mean_front_x() + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_ddt_detection() {
        let mut ddt = DdtModel::new(2.0, 0.5, 1500.0);
        let pressure = vec![1.0, 3.0, 1.0];
        let rate = vec![0.0, 0.6, 0.0];
        let detected = ddt.check_ddt(&pressure, 1.0, &rate, 0.01);
        assert!(detected);
    }

    #[test]
    fn test_ddt_no_detection() {
        let mut ddt = DdtModel::new(5.0, 1.0, 1500.0);
        let pressure = vec![1.0, 1.0, 1.0];
        let rate = vec![0.1, 0.1, 0.1];
        let detected = ddt.check_ddt(&pressure, 1.0, &rate, 0.01);
        assert!(!detected);
    }

    #[test]
    fn test_ddt_cj_mach() {
        let ddt = DdtModel::new(2.0, 0.5, 1500.0);
        let ma = ddt.cj_mach(340.0);
        assert!((ma - 1500.0 / 340.0).abs() < 1e-10);
    }

    #[test]
    fn test_flame_speed_damkohler() {
        let fs = FlameSpeed::new(0.37, 2e-5, 1e-3, 0.5);
        let da = fs.damkohler_number();
        assert!(da > 0.0);
    }

    #[test]
    fn test_flame_speed_turbulent_damkohler() {
        let fs = FlameSpeed::new(0.37, 2e-5, 1e-3, 0.5);
        assert!((fs.turbulent_damkohler() - 0.87).abs() < 1e-10);
    }

    #[test]
    fn test_flame_speed_borghi_regime() {
        let fs = FlameSpeed::new(0.37, 2e-5, 1e-3, 0.1);
        let _regime = fs.borghi_regime();
    }

    #[test]
    fn test_soot_model_step() {
        let mut soot = SootModel::new(4, 1e-3, 1e-4, 1e-2, 1800.0);
        let temp = vec![1600.0; 4];
        let yf = vec![0.05; 4];
        let yo2 = vec![0.23; 4];
        soot.step(&temp, &yf, &yo2, 0.01);
        // Soot should form at high T with fuel present
        assert!(soot.fv.iter().any(|&v| v > 0.0));
    }

    #[test]
    fn test_soot_model_no_soot_low_temp() {
        let mut soot = SootModel::new(4, 1e-3, 1e-4, 0.0, 1800.0);
        let temp = vec![500.0; 4];
        let yf = vec![0.05; 4];
        let yo2 = vec![0.23; 4];
        soot.step(&temp, &yf, &yo2, 0.01);
        assert!(soot.mean_fv() == 0.0);
    }

    #[test]
    fn test_nox_thermal_formation() {
        let mut nox = NoxModel::new(4, 6e16, 69090.0, 1e10);
        let temp = vec![2000.0; 4];
        let yn2 = vec![0.72; 4];
        let yo2 = vec![0.23; 4];
        nox.step_thermal(&temp, &yn2, &yo2, 0.01);
        assert!(nox.mean_no() > 0.0);
    }

    #[test]
    fn test_nox_below_threshold() {
        let mut nox = NoxModel::new(4, 6e16, 69090.0, 1e10);
        let temp = vec![1000.0; 4]; // Below 1500 K threshold
        let yn2 = vec![0.72; 4];
        let yo2 = vec![0.23; 4];
        nox.step_thermal(&temp, &yn2, &yo2, 0.01);
        assert_eq!(nox.mean_no(), 0.0);
    }

    #[test]
    fn test_co_model_formation() {
        let mut co = CoModel::new(4, 0.1, 0.0, 20000.0);
        let temp = vec![1800.0; 4];
        let yf = vec![0.05; 4];
        let yo2 = vec![0.23; 4];
        co.step(&temp, &yf, &yo2, 0.1);
        assert!(co.mean_co() > 0.0);
    }

    #[test]
    fn test_equilibrium_d2q9_sum() {
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bgk_collision_reduces_error() {
        let f = equilibrium_d2q9(1.0, 0.05, 0.0);
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        let f_out = bgk_collision(f, feq, 1.0);
        let err_before: f64 = f.iter().zip(feq.iter()).map(|(a, b)| (a - b).abs()).sum();
        let err_after: f64 = f_out
            .iter()
            .zip(feq.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(err_after <= err_before);
    }

    #[test]
    fn test_combustion_lbm_create() {
        let k = ArrheniusKinetics::methane_air();
        let solver = CombustionLbm::new(16, 8, 1.0, 1.2, k, 300.0, 0.01);
        assert_eq!(solver.nx, 16);
        assert_eq!(solver.ny, 8);
    }

    #[test]
    fn test_combustion_lbm_advance() {
        let k = ArrheniusKinetics::methane_air();
        let mut solver = CombustionLbm::new(8, 4, 1.0, 1.2, k, 300.0, 0.01);
        solver.ignite(4, 2, 1500.0);
        solver.advance();
        assert_eq!(solver.step, 1);
    }

    #[test]
    fn test_premixed_config_adiabatic_temp() {
        let cfg = PremixedFlameConfig::methane_stoich();
        let t_ad = cfg.adiabatic_temp();
        assert!(t_ad > 300.0);
    }

    #[test]
    fn test_premixed_config_flame_speed() {
        let cfg = PremixedFlameConfig::methane_stoich();
        let sl = cfg.laminar_flame_speed();
        assert!(sl > 0.0 && sl < 5.0);
    }

    #[test]
    fn test_diffusion_flame_bs_temp() {
        let df = DiffusionFlameConfig::methane_air_default();
        let t_f = df.bs_flame_temp(5e7, 1200.0);
        assert!(t_f > 300.0);
    }

    #[test]
    fn test_extinction_model_extinguished() {
        // Da = 1/(chi*tau_chem) = 1/(10.0*0.5) = 0.2 < da_crit=1.0 → extinguished
        let ext = ExtinctionModel::new(1.0, 10.0, 0.5);
        assert!(ext.is_extinguished());
    }

    #[test]
    fn test_extinction_model_not_extinguished() {
        let ext = ExtinctionModel::new(1.0, 0.01, 0.01);
        assert!(!ext.is_extinguished());
    }

    #[test]
    fn test_scalar_dissipation() {
        let chi = scalar_dissipation(2e-5, [0.1, 0.0]);
        assert!((chi - 2.0 * 2e-5 * 0.01).abs() < 1e-20);
    }

    #[test]
    fn test_mixture_fraction_gradient() {
        let z = vec![0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.2, 0.3];
        let grad = mixture_fraction_gradient(&z, 4, 2, 1, 0, 1.0, 1.0);
        assert!((grad[0] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_thermo_diffusive_instability_le1() {
        let tdi = ThermoDiffusiveInstability::new(1.0, 10.0);
        assert!(!tdi.is_cellular());
    }

    #[test]
    fn test_thermo_diffusive_instability_le_low() {
        let tdi = ThermoDiffusiveInstability::new(0.5, 10.0);
        assert!(tdi.is_cellular());
    }

    #[test]
    fn test_flame_surface_density_step() {
        let mut fsd = FlameSurfaceDensity::new(4, 10.0, 1.0);
        fsd.sigma = vec![0.01; 4];
        fsd.step(0.5, 0.37, 0.01);
        assert!(fsd.sigma.iter().all(|&s| s >= 0.0));
    }

    #[test]
    fn test_lewis_omega_clamped() {
        let le = LewisEffect::new(1e-20, 2e-5);
        let om = le.species_omega(1.0);
        assert!(om.is_finite() && om > 0.0);
    }
}
