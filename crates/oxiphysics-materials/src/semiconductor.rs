// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Semiconductor physics: band gap, carrier transport, junctions, and device models.
//!
//! This module provides comprehensive semiconductor physics modeling including
//! band structure, carrier concentrations, drift-diffusion transport, p-n junctions,
//! Schottky barriers, recombination mechanisms, MOSFET models, and analysis tools.

use std::f64::consts::PI;

/// Boltzmann constant in eV/K.
pub const K_BOLTZMANN_EV: f64 = 8.617333262e-5;

/// Boltzmann constant in J/K.
pub const K_BOLTZMANN_J: f64 = 1.380649e-23;

/// Elementary charge in C.
pub const Q_ELECTRON: f64 = 1.602176634e-19;

/// Permittivity of free space in F/m.
pub const EPSILON_0: f64 = 8.854187817e-12;

/// Planck constant in J·s.
pub const H_PLANCK: f64 = 6.62607015e-34;

/// Reduced Planck constant in J·s.
pub const HBAR: f64 = 1.054571817e-34;

/// Free electron mass in kg.
pub const M_ELECTRON: f64 = 9.1093837015e-31;

/// Thermal voltage at temperature T: V_T = k_B * T / q.
pub fn thermal_voltage(temperature_k: f64) -> f64 {
    K_BOLTZMANN_EV * temperature_k
}

// ---------------------------------------------------------------------------
// BandStructure
// ---------------------------------------------------------------------------

/// Band structure parameters of a semiconductor.
///
/// Characterizes the energy band structure including band gap,
/// band edges, and effective masses.
#[derive(Debug, Clone, Copy)]
pub struct BandStructure {
    /// Band gap energy in eV.
    pub band_gap_ev: f64,
    /// Conduction band edge energy in eV (relative to vacuum).
    pub ec: f64,
    /// Valence band edge energy in eV (relative to vacuum).
    pub ev: f64,
    /// Electron effective mass ratio (m*/m0).
    pub me_eff: f64,
    /// Hole effective mass ratio (m*/m0).
    pub mh_eff: f64,
    /// Dielectric constant (relative permittivity).
    pub dielectric_constant: f64,
}

impl BandStructure {
    /// Create a new band structure.
    pub fn new(
        band_gap_ev: f64,
        electron_affinity: f64,
        me_eff: f64,
        mh_eff: f64,
        dielectric_constant: f64,
    ) -> Self {
        Self {
            band_gap_ev,
            ec: -electron_affinity,
            ev: -electron_affinity - band_gap_ev,
            me_eff,
            mh_eff,
            dielectric_constant,
        }
    }

    /// Silicon band structure at 300K.
    pub fn silicon() -> Self {
        Self::new(1.12, 4.05, 1.08, 0.56, 11.7)
    }

    /// Germanium band structure at 300K.
    pub fn germanium() -> Self {
        Self::new(0.66, 4.0, 0.55, 0.37, 16.0)
    }

    /// Gallium arsenide band structure at 300K.
    pub fn gallium_arsenide() -> Self {
        Self::new(1.42, 4.07, 0.067, 0.45, 12.9)
    }

    /// Silicon carbide (4H-SiC) band structure.
    pub fn silicon_carbide() -> Self {
        Self::new(3.26, 3.17, 0.39, 1.0, 9.7)
    }

    /// Gallium nitride band structure.
    pub fn gallium_nitride() -> Self {
        Self::new(3.4, 4.1, 0.2, 0.8, 8.9)
    }

    /// Temperature-dependent band gap using Varshni equation.
    ///
    /// Eg(T) = Eg(0) - alpha * T^2 / (T + beta)
    pub fn varshni_band_gap(&self, temperature_k: f64, alpha: f64, beta: f64) -> f64 {
        self.band_gap_ev + alpha * 300.0 * 300.0 / (300.0 + beta)
            - alpha * temperature_k * temperature_k / (temperature_k + beta)
    }

    /// Effective density of states in conduction band (Nc) in cm^-3.
    pub fn nc(&self, temperature_k: f64) -> f64 {
        let factor = 2.0 * PI * self.me_eff * M_ELECTRON * K_BOLTZMANN_J * temperature_k
            / (H_PLANCK * H_PLANCK);
        2.0 * factor.powf(1.5) * 1e-6 // m^-3 to cm^-3
    }

    /// Effective density of states in valence band (Nv) in cm^-3.
    pub fn nv(&self, temperature_k: f64) -> f64 {
        let factor = 2.0 * PI * self.mh_eff * M_ELECTRON * K_BOLTZMANN_J * temperature_k
            / (H_PLANCK * H_PLANCK);
        2.0 * factor.powf(1.5) * 1e-6 // m^-3 to cm^-3
    }

    /// Absolute permittivity in F/m.
    pub fn permittivity(&self) -> f64 {
        self.dielectric_constant * EPSILON_0
    }

    /// Midgap energy (average of Ec and Ev).
    pub fn midgap_energy(&self) -> f64 {
        (self.ec + self.ev) / 2.0
    }
}

// ---------------------------------------------------------------------------
// CarrierConcentration
// ---------------------------------------------------------------------------

/// Doping type of a semiconductor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DopingType {
    /// N-type: donor impurities (excess electrons).
    NType,
    /// P-type: acceptor impurities (excess holes).
    PType,
    /// Intrinsic: undoped semiconductor.
    Intrinsic,
}

/// Carrier concentration calculations for a semiconductor.
///
/// Computes intrinsic carrier concentration, doped carrier densities,
/// and Fermi level positions.
#[derive(Debug, Clone, Copy)]
pub struct CarrierConcentration {
    /// Band structure of the material.
    pub band: BandStructure,
    /// Temperature in Kelvin.
    pub temperature_k: f64,
    /// Donor concentration in cm^-3.
    pub nd: f64,
    /// Acceptor concentration in cm^-3.
    pub na: f64,
}

impl CarrierConcentration {
    /// Create a new carrier concentration model.
    pub fn new(band: BandStructure, temperature_k: f64, nd: f64, na: f64) -> Self {
        Self {
            band,
            temperature_k,
            nd,
            na,
        }
    }

    /// Intrinsic carrier concentration ni in cm^-3.
    ///
    /// ni = sqrt(Nc * Nv) * exp(-Eg / (2 * k_B * T))
    pub fn intrinsic_concentration(&self) -> f64 {
        let nc = self.band.nc(self.temperature_k);
        let nv = self.band.nv(self.temperature_k);
        let vt = thermal_voltage(self.temperature_k);
        (nc * nv).sqrt() * (-self.band.band_gap_ev / (2.0 * vt)).exp()
    }

    /// Fermi-Dirac integral of order 1/2 (approximation).
    ///
    /// F_{1/2}(eta) approximated using the Joyce-Dixon formula.
    pub fn fermi_dirac_half(eta: f64) -> f64 {
        if eta < -5.0 {
            // Non-degenerate: F_{1/2}(eta) ≈ exp(eta) * sqrt(pi)/2
            (PI.sqrt() / 2.0) * eta.exp()
        } else if eta > 5.0 {
            // Highly degenerate: F_{1/2}(eta) ≈ (2/3) * eta^{3/2}
            (2.0 / 3.0) * eta.powf(1.5)
        } else {
            // Intermediate: Padé approximation
            let exp_eta = eta.exp();
            let denom = 1.0 + 0.27 * exp_eta;
            (PI.sqrt() / 2.0) * exp_eta / denom.sqrt()
        }
    }

    /// Electron concentration n0 in cm^-3 (assuming complete ionization).
    pub fn electron_concentration(&self) -> f64 {
        let ni = self.intrinsic_concentration();
        let net_donor = self.nd - self.na;
        if net_donor > 0.0 {
            // N-type
            let discriminant = net_donor * net_donor + 4.0 * ni * ni;
            (net_donor + discriminant.sqrt()) / 2.0
        } else if net_donor < 0.0 {
            // P-type: n = ni^2 / p
            let net_acceptor = -net_donor;
            let discriminant = net_acceptor * net_acceptor + 4.0 * ni * ni;
            let p = (net_acceptor + discriminant.sqrt()) / 2.0;
            ni * ni / p
        } else {
            ni
        }
    }

    /// Hole concentration p0 in cm^-3.
    pub fn hole_concentration(&self) -> f64 {
        let ni = self.intrinsic_concentration();
        let n0 = self.electron_concentration();
        ni * ni / n0
    }

    /// Doping type of the material.
    pub fn doping_type(&self) -> DopingType {
        if (self.nd - self.na).abs() < 1e-6 * (self.nd + self.na + 1.0) {
            DopingType::Intrinsic
        } else if self.nd > self.na {
            DopingType::NType
        } else {
            DopingType::PType
        }
    }

    /// Fermi level position relative to intrinsic Fermi level in eV.
    ///
    /// Ef - Ei = k_B * T * ln(n0 / ni)
    pub fn fermi_level_shift(&self) -> f64 {
        let ni = self.intrinsic_concentration();
        let n0 = self.electron_concentration();
        let vt = thermal_voltage(self.temperature_k);
        vt * (n0 / ni).ln()
    }

    /// Compensation ratio: min(Nd, Na) / max(Nd, Na).
    pub fn compensation_ratio(&self) -> f64 {
        let max_doping = self.nd.max(self.na);
        if max_doping < 1e-10 {
            return 0.0;
        }
        let min_doping = self.nd.min(self.na);
        min_doping / max_doping
    }

    /// Debye length in cm.
    pub fn debye_length(&self) -> f64 {
        let vt = thermal_voltage(self.temperature_k);
        let n0 = self.electron_concentration();
        let p0 = self.hole_concentration();
        let epsilon = self.band.permittivity();
        // L_D = sqrt(epsilon * V_T / (q * (n0 + p0)))
        let total_carriers = (n0 + p0) * 1e6; // cm^-3 to m^-3
        let ld = (epsilon * vt / (Q_ELECTRON * total_carriers)).sqrt();
        ld * 100.0 // m to cm
    }
}

// ---------------------------------------------------------------------------
// DriftDiffusion
// ---------------------------------------------------------------------------

/// Drift-diffusion transport model for semiconductors.
///
/// Models carrier transport including mobility, diffusion, and velocity saturation.
#[derive(Debug, Clone, Copy)]
pub struct DriftDiffusion {
    /// Electron mobility in cm^2/(V·s).
    pub mu_n: f64,
    /// Hole mobility in cm^2/(V·s).
    pub mu_p: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Saturation velocity for electrons in cm/s.
    pub v_sat_n: f64,
    /// Saturation velocity for holes in cm/s.
    pub v_sat_p: f64,
}

impl DriftDiffusion {
    /// Create a new drift-diffusion model.
    pub fn new(mu_n: f64, mu_p: f64, temperature_k: f64, v_sat_n: f64, v_sat_p: f64) -> Self {
        Self {
            mu_n,
            mu_p,
            temperature_k,
            v_sat_n,
            v_sat_p,
        }
    }

    /// Silicon at 300K with typical mobility values.
    pub fn silicon_300k() -> Self {
        Self::new(1400.0, 450.0, 300.0, 1.07e7, 8.37e6)
    }

    /// GaAs at 300K.
    pub fn gaas_300k() -> Self {
        Self::new(8500.0, 400.0, 300.0, 7.7e6, 7.7e6)
    }

    /// Electron diffusion coefficient via Einstein relation: D_n = mu_n * V_T.
    pub fn diffusion_coefficient_n(&self) -> f64 {
        self.mu_n * thermal_voltage(self.temperature_k)
    }

    /// Hole diffusion coefficient via Einstein relation: D_p = mu_p * V_T.
    pub fn diffusion_coefficient_p(&self) -> f64 {
        self.mu_p * thermal_voltage(self.temperature_k)
    }

    /// Verify Einstein relation: D/mu = k_B T / q.
    pub fn einstein_ratio(&self) -> f64 {
        thermal_voltage(self.temperature_k)
    }

    /// Electron drift velocity with saturation effects.
    ///
    /// v_d = mu * E / (1 + mu * |E| / v_sat)
    pub fn electron_drift_velocity(&self, electric_field: f64) -> f64 {
        let e_abs = electric_field.abs();
        let v_low = self.mu_n * e_abs;
        v_low / (1.0 + v_low / self.v_sat_n) * electric_field.signum()
    }

    /// Hole drift velocity with saturation effects.
    pub fn hole_drift_velocity(&self, electric_field: f64) -> f64 {
        let e_abs = electric_field.abs();
        let v_low = self.mu_p * e_abs;
        v_low / (1.0 + v_low / self.v_sat_p) * electric_field.signum()
    }

    /// Electron current density (drift + diffusion) in A/cm^2.
    ///
    /// J_n = q * n * mu_n * E + q * D_n * dn/dx
    pub fn electron_current_density(&self, n: f64, electric_field: f64, dn_dx: f64) -> f64 {
        let drift = Q_ELECTRON * n * self.mu_n * electric_field;
        let diffusion = Q_ELECTRON * self.diffusion_coefficient_n() * dn_dx;
        // Convert from C·cm^-3·cm^2/(V·s)·V/cm = A/cm^2
        drift + diffusion
    }

    /// Hole current density (drift + diffusion) in A/cm^2.
    pub fn hole_current_density(&self, p: f64, electric_field: f64, dp_dx: f64) -> f64 {
        let drift = Q_ELECTRON * p * self.mu_p * electric_field;
        let diffusion = -Q_ELECTRON * self.diffusion_coefficient_p() * dp_dx;
        drift + diffusion
    }

    /// Mobility with doping dependence (Caughey-Thomas model).
    ///
    /// mu = mu_min + (mu_max - mu_min) / (1 + (N / N_ref)^alpha)
    pub fn caughey_thomas_mobility(
        mu_min: f64,
        mu_max: f64,
        n_total: f64,
        n_ref: f64,
        alpha: f64,
    ) -> f64 {
        mu_min + (mu_max - mu_min) / (1.0 + (n_total / n_ref).powf(alpha))
    }

    /// Mean free path estimate from mobility in cm.
    pub fn mean_free_path_n(&self) -> f64 {
        let vth = (3.0 * K_BOLTZMANN_J * self.temperature_k / M_ELECTRON).sqrt();
        let vth_cm = vth * 100.0; // m/s to cm/s
        self.mu_n * M_ELECTRON * vth_cm / (Q_ELECTRON * 100.0)
    }
}

// ---------------------------------------------------------------------------
// PnJunction
// ---------------------------------------------------------------------------

/// P-N junction model.
///
/// Models the electrostatics and I-V characteristics of an abrupt p-n junction.
#[derive(Debug, Clone, Copy)]
pub struct PnJunction {
    /// Band structure.
    pub band: BandStructure,
    /// Donor concentration (n-side) in cm^-3.
    pub nd: f64,
    /// Acceptor concentration (p-side) in cm^-3.
    pub na: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Junction area in cm^2.
    pub area: f64,
}

impl PnJunction {
    /// Create a new p-n junction.
    pub fn new(band: BandStructure, nd: f64, na: f64, temperature_k: f64, area: f64) -> Self {
        Self {
            band,
            nd,
            na,
            temperature_k,
            area,
        }
    }

    /// Intrinsic carrier concentration.
    pub fn ni(&self) -> f64 {
        let cc = CarrierConcentration::new(self.band, self.temperature_k, 0.0, 0.0);
        cc.intrinsic_concentration()
    }

    /// Built-in potential Vbi in V.
    ///
    /// Vbi = (k_B T / q) * ln(Na * Nd / ni^2)
    pub fn built_in_potential(&self) -> f64 {
        let ni = self.ni();
        let vt = thermal_voltage(self.temperature_k);
        vt * (self.na * self.nd / (ni * ni)).ln()
    }

    /// Depletion width W in cm.
    ///
    /// W = sqrt(2 * epsilon * (Vbi - V) * (1/Na + 1/Nd) / q)
    pub fn depletion_width(&self, applied_voltage: f64) -> f64 {
        let vbi = self.built_in_potential();
        let v_eff = vbi - applied_voltage;
        if v_eff <= 0.0 {
            return 0.0;
        }
        let epsilon = self.band.permittivity();
        let factor = 2.0 * epsilon * v_eff / Q_ELECTRON;
        let doping_factor = 1.0 / (self.na * 1e6) + 1.0 / (self.nd * 1e6); // cm^-3 to m^-3
        (factor * doping_factor).sqrt() * 100.0 // m to cm
    }

    /// Depletion width on the n-side in cm.
    pub fn xn(&self, applied_voltage: f64) -> f64 {
        let w = self.depletion_width(applied_voltage);
        w * self.na / (self.na + self.nd)
    }

    /// Depletion width on the p-side in cm.
    pub fn xp(&self, applied_voltage: f64) -> f64 {
        let w = self.depletion_width(applied_voltage);
        w * self.nd / (self.na + self.nd)
    }

    /// Junction capacitance per unit area in F/cm^2.
    ///
    /// C_j = epsilon / W
    pub fn junction_capacitance(&self, applied_voltage: f64) -> f64 {
        let w = self.depletion_width(applied_voltage);
        if w < 1e-20 {
            return 0.0;
        }
        let epsilon = self.band.permittivity();
        epsilon / (w * 1e-2) * 1e-4 // F/m / m -> F/cm^2
    }

    /// Total junction capacitance in F.
    pub fn total_capacitance(&self, applied_voltage: f64) -> f64 {
        self.junction_capacitance(applied_voltage) * self.area
    }

    /// Maximum electric field in V/cm.
    pub fn max_electric_field(&self, applied_voltage: f64) -> f64 {
        let w = self.depletion_width(applied_voltage);
        if w < 1e-20 {
            return 0.0;
        }
        let vbi = self.built_in_potential();
        2.0 * (vbi - applied_voltage) / w
    }

    /// Reverse saturation current I0 in A.
    ///
    /// Uses the Shockley diode equation parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn saturation_current(&self, dn: f64, dp: f64, ln: f64, lp: f64) -> f64 {
        let ni = self.ni();
        let term_n = Q_ELECTRON * dn * ni * ni / (ln * self.nd);
        let term_p = Q_ELECTRON * dp * ni * ni / (lp * self.na);
        self.area * (term_n + term_p)
    }

    /// Current through the junction using Shockley diode equation.
    ///
    /// I = I0 * (exp(V / (n * V_T)) - 1)
    pub fn shockley_current(&self, voltage: f64, i0: f64, ideality_factor: f64) -> f64 {
        let vt = thermal_voltage(self.temperature_k);
        i0 * ((voltage / (ideality_factor * vt)).exp() - 1.0)
    }

    /// Breakdown voltage estimate (empirical, for silicon).
    ///
    /// V_br ≈ 60 * (Eg / 1.1)^{3/2} * (Nb / 10^16)^{-3/4}
    pub fn breakdown_voltage(&self) -> f64 {
        let nb = self.nd.min(self.na);
        60.0 * (self.band.band_gap_ev / 1.1).powf(1.5) * (nb / 1e16).powf(-0.75)
    }
}

// ---------------------------------------------------------------------------
// SchottkyBarrier
// ---------------------------------------------------------------------------

/// Schottky barrier (metal-semiconductor) junction model.
#[derive(Debug, Clone, Copy)]
pub struct SchottkyBarrier {
    /// Metal work function in eV.
    pub metal_work_function: f64,
    /// Semiconductor band structure.
    pub band: BandStructure,
    /// Doping concentration in cm^-3.
    pub doping: f64,
    /// Doping type.
    pub doping_type: DopingType,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Ideality factor (1 for ideal).
    pub ideality_factor: f64,
    /// Junction area in cm^2.
    pub area: f64,
}

impl SchottkyBarrier {
    /// Create a new Schottky barrier.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        metal_work_function: f64,
        band: BandStructure,
        doping: f64,
        doping_type: DopingType,
        temperature_k: f64,
        ideality_factor: f64,
        area: f64,
    ) -> Self {
        Self {
            metal_work_function,
            band,
            doping,
            doping_type,
            temperature_k,
            ideality_factor,
            area,
        }
    }

    /// Barrier height for n-type semiconductor: phi_B = phi_M - chi.
    pub fn barrier_height(&self) -> f64 {
        let electron_affinity = -self.band.ec;
        match self.doping_type {
            DopingType::NType => self.metal_work_function - electron_affinity,
            DopingType::PType => {
                electron_affinity + self.band.band_gap_ev - self.metal_work_function
            }
            DopingType::Intrinsic => self.metal_work_function - electron_affinity,
        }
    }

    /// Built-in potential for the Schottky contact.
    pub fn built_in_potential(&self) -> f64 {
        let vt = thermal_voltage(self.temperature_k);
        let nc = self.band.nc(self.temperature_k);
        self.barrier_height() - vt * (nc / self.doping).ln()
    }

    /// Depletion width in cm.
    pub fn depletion_width(&self, applied_voltage: f64) -> f64 {
        let vbi = self.built_in_potential();
        let v_eff = vbi - applied_voltage;
        if v_eff <= 0.0 {
            return 0.0;
        }
        let epsilon = self.band.permittivity();
        let w = (2.0 * epsilon * v_eff / (Q_ELECTRON * self.doping * 1e6)).sqrt();
        w * 100.0 // m to cm
    }

    /// Richardson constant for the semiconductor in A/(cm^2·K^2).
    pub fn richardson_constant(&self) -> f64 {
        let me = match self.doping_type {
            DopingType::NType | DopingType::Intrinsic => self.band.me_eff,
            DopingType::PType => self.band.mh_eff,
        };
        // A* = 4 * pi * q * m* * k_B^2 / h^3
        // Standard value for free electron: 120 A/(cm^2·K^2)
        120.0 * me
    }

    /// Saturation current density in A/cm^2.
    pub fn saturation_current_density(&self) -> f64 {
        let a_star = self.richardson_constant();
        let phi_b = self.barrier_height();
        let vt = thermal_voltage(self.temperature_k);
        a_star * self.temperature_k * self.temperature_k * (-phi_b / vt).exp()
    }

    /// Current through the Schottky junction in A.
    pub fn current(&self, voltage: f64) -> f64 {
        let js = self.saturation_current_density();
        let vt = thermal_voltage(self.temperature_k);
        let i0 = js * self.area;
        i0 * ((voltage / (self.ideality_factor * vt)).exp() - 1.0)
    }
}

// ---------------------------------------------------------------------------
// RecombinationModel
// ---------------------------------------------------------------------------

/// Recombination mechanisms in semiconductors.
///
/// Models Shockley-Read-Hall (SRH), Auger, and radiative recombination.
#[derive(Debug, Clone, Copy)]
pub struct RecombinationModel {
    /// Electron lifetime in s (for SRH).
    pub tau_n: f64,
    /// Hole lifetime in s (for SRH).
    pub tau_p: f64,
    /// Auger coefficient for electrons in cm^6/s.
    pub cn: f64,
    /// Auger coefficient for holes in cm^6/s.
    pub cp: f64,
    /// Radiative recombination coefficient in cm^3/s.
    pub brad: f64,
    /// Intrinsic carrier concentration in cm^-3.
    pub ni: f64,
}

impl RecombinationModel {
    /// Create a new recombination model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(tau_n: f64, tau_p: f64, cn: f64, cp: f64, brad: f64, ni: f64) -> Self {
        Self {
            tau_n,
            tau_p,
            cn,
            cp,
            brad,
            ni,
        }
    }

    /// Typical silicon recombination parameters at 300K.
    pub fn silicon_300k(ni: f64) -> Self {
        Self::new(
            1e-5,    // tau_n = 10 us
            1e-5,    // tau_p = 10 us
            2.8e-31, // Cn
            9.9e-32, // Cp
            1.1e-14, // Brad (indirect gap, very small)
            ni,
        )
    }

    /// SRH recombination rate in cm^-3/s.
    ///
    /// R_SRH = (n*p - ni^2) / (tau_p*(n + ni) + tau_n*(p + ni))
    pub fn srh_rate(&self, n: f64, p: f64) -> f64 {
        let np = n * p;
        let ni2 = self.ni * self.ni;
        let numerator = np - ni2;
        let denominator = self.tau_p * (n + self.ni) + self.tau_n * (p + self.ni);
        if denominator.abs() < 1e-50 {
            return 0.0;
        }
        numerator / denominator
    }

    /// Auger recombination rate in cm^-3/s.
    ///
    /// R_Auger = (Cn*n + Cp*p) * (n*p - ni^2)
    pub fn auger_rate(&self, n: f64, p: f64) -> f64 {
        let ni2 = self.ni * self.ni;
        (self.cn * n + self.cp * p) * (n * p - ni2)
    }

    /// Radiative recombination rate in cm^-3/s.
    ///
    /// R_rad = Brad * (n*p - ni^2)
    pub fn radiative_rate(&self, n: f64, p: f64) -> f64 {
        let ni2 = self.ni * self.ni;
        self.brad * (n * p - ni2)
    }

    /// Total recombination rate (sum of all mechanisms) in cm^-3/s.
    pub fn total_rate(&self, n: f64, p: f64) -> f64 {
        self.srh_rate(n, p) + self.auger_rate(n, p) + self.radiative_rate(n, p)
    }

    /// Effective minority carrier lifetime in s.
    ///
    /// 1/tau_eff = 1/tau_SRH + 1/tau_Auger + 1/tau_rad
    pub fn effective_lifetime(&self, n: f64, p: f64) -> f64 {
        let r_total = self.total_rate(n, p);
        let ni2 = self.ni * self.ni;
        let delta_n = (n * p - ni2).abs();
        if delta_n < 1e-10 || r_total.abs() < 1e-50 {
            return self.tau_n.min(self.tau_p);
        }
        // Assuming low-level injection: delta_n ≈ delta carrier concentration
        let majority = n.max(p);
        majority / r_total.abs()
    }

    /// Surface recombination velocity model in cm/s.
    ///
    /// R_surface = S * (n_s - n_s0)
    pub fn surface_recombination_rate(
        surface_velocity: f64,
        n_surface: f64,
        n_equilibrium: f64,
    ) -> f64 {
        surface_velocity * (n_surface - n_equilibrium)
    }
}

// ---------------------------------------------------------------------------
// MosfetModel
// ---------------------------------------------------------------------------

/// MOSFET (Metal-Oxide-Semiconductor Field-Effect Transistor) model.
///
/// Implements the square-law MOSFET model with subthreshold behavior.
#[derive(Debug, Clone, Copy)]
pub struct MosfetModel {
    /// Threshold voltage Vth in V.
    pub vth: f64,
    /// Oxide capacitance per unit area Cox in F/cm^2.
    pub cox: f64,
    /// Channel mobility in cm^2/(V·s).
    pub mobility: f64,
    /// Channel width in cm.
    pub width: f64,
    /// Channel length in cm.
    pub length: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Subthreshold swing ideality factor (typically 1.0-1.5).
    pub n_sub: f64,
    /// Channel-length modulation parameter lambda in V^-1.
    pub lambda: f64,
}

impl MosfetModel {
    /// Create a new MOSFET model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        vth: f64,
        cox: f64,
        mobility: f64,
        width: f64,
        length: f64,
        temperature_k: f64,
        n_sub: f64,
        lambda: f64,
    ) -> Self {
        Self {
            vth,
            cox,
            mobility,
            width,
            length,
            temperature_k,
            n_sub,
            lambda,
        }
    }

    /// Typical NMOS transistor parameters.
    pub fn typical_nmos() -> Self {
        Self::new(
            0.5,     // Vth = 0.5V
            3.45e-7, // Cox for 10nm oxide
            400.0,   // mobility
            10.0e-4, // W = 10um
            0.18e-4, // L = 0.18um
            300.0,   // T = 300K
            1.3,     // n_sub
            0.05,    // lambda
        )
    }

    /// Transconductance parameter K = mu * Cox * W / L in A/V^2.
    pub fn transconductance_parameter(&self) -> f64 {
        self.mobility * self.cox * self.width / self.length
    }

    /// Drain current in the linear (triode) region in A.
    ///
    /// Id = K * ((Vgs - Vth) * Vds - Vds^2 / 2) * (1 + lambda * Vds)
    pub fn linear_current(&self, vgs: f64, vds: f64) -> f64 {
        if vgs <= self.vth || vds <= 0.0 {
            return 0.0;
        }
        let k = self.transconductance_parameter();
        let vov = vgs - self.vth;
        let vds_eff = vds.min(vov); // Clamp to saturation boundary
        k * (vov * vds_eff - vds_eff * vds_eff / 2.0) * (1.0 + self.lambda * vds)
    }

    /// Drain current in the saturation region in A.
    ///
    /// Id = (K / 2) * (Vgs - Vth)^2 * (1 + lambda * Vds)
    pub fn saturation_current(&self, vgs: f64, vds: f64) -> f64 {
        if vgs <= self.vth {
            return 0.0;
        }
        let k = self.transconductance_parameter();
        let vov = vgs - self.vth;
        (k / 2.0) * vov * vov * (1.0 + self.lambda * vds)
    }

    /// Drain current combining linear and saturation regions.
    pub fn drain_current(&self, vgs: f64, vds: f64) -> f64 {
        if vgs <= self.vth {
            return self.subthreshold_current(vgs, vds);
        }
        let vov = vgs - self.vth;
        if vds < vov {
            self.linear_current(vgs, vds)
        } else {
            self.saturation_current(vgs, vds)
        }
    }

    /// Subthreshold current in A.
    ///
    /// I_sub = I0 * exp((Vgs - Vth) / (n * V_T)) * (1 - exp(-Vds / V_T))
    pub fn subthreshold_current(&self, vgs: f64, vds: f64) -> f64 {
        let vt = thermal_voltage(self.temperature_k);
        let k = self.transconductance_parameter();
        // Subthreshold reference current
        let i0 = k * (self.n_sub * vt) * (self.n_sub * vt) / 2.0;
        let exp_gate = ((vgs - self.vth) / (self.n_sub * vt)).exp();
        let exp_drain = 1.0 - (-vds / vt).exp();
        i0 * exp_gate * exp_drain
    }

    /// Subthreshold swing in mV/decade.
    pub fn subthreshold_swing(&self) -> f64 {
        let vt = thermal_voltage(self.temperature_k);
        self.n_sub * vt * (10.0_f64).ln() * 1000.0
    }

    /// Transconductance gm = dId/dVgs in the saturation region, in A/V.
    pub fn transconductance(&self, vgs: f64) -> f64 {
        if vgs <= self.vth {
            return 0.0;
        }
        let k = self.transconductance_parameter();
        k * (vgs - self.vth)
    }

    /// Output conductance gds = dId/dVds in saturation, in A/V.
    pub fn output_conductance(&self, vgs: f64) -> f64 {
        if vgs <= self.vth {
            return 0.0;
        }
        let k = self.transconductance_parameter();
        let vov = vgs - self.vth;
        (k / 2.0) * vov * vov * self.lambda
    }

    /// Intrinsic voltage gain Av = gm / gds.
    pub fn intrinsic_gain(&self, vgs: f64) -> f64 {
        let gds = self.output_conductance(vgs);
        if gds.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.transconductance(vgs) / gds
    }

    /// Threshold voltage with body effect.
    ///
    /// Vth = Vth0 + gamma * (sqrt(2*phi_f + Vsb) - sqrt(2*phi_f))
    pub fn body_effect_vth(vth0: f64, gamma: f64, phi_f: f64, vsb: f64) -> f64 {
        vth0 + gamma * ((2.0 * phi_f + vsb).sqrt() - (2.0 * phi_f).sqrt())
    }
}

// ---------------------------------------------------------------------------
// SemiconductorAnalysis
// ---------------------------------------------------------------------------

/// Analysis tools for semiconductor properties.
///
/// Provides calculations for resistivity, sheet resistance, and Hall measurements.
pub struct SemiconductorAnalysis;

impl SemiconductorAnalysis {
    /// Resistivity in Ohm·cm.
    ///
    /// rho = 1 / (q * (n * mu_n + p * mu_p))
    pub fn resistivity(n: f64, mu_n: f64, p: f64, mu_p: f64) -> f64 {
        let sigma = Q_ELECTRON * (n * mu_n + p * mu_p);
        if sigma < 1e-50 {
            return f64::INFINITY;
        }
        1.0 / sigma
    }

    /// Conductivity in (Ohm·cm)^-1.
    pub fn conductivity(n: f64, mu_n: f64, p: f64, mu_p: f64) -> f64 {
        Q_ELECTRON * (n * mu_n + p * mu_p)
    }

    /// Sheet resistance in Ohm/square.
    pub fn sheet_resistance(resistivity: f64, thickness: f64) -> f64 {
        if thickness < 1e-30 {
            return f64::INFINITY;
        }
        resistivity / thickness
    }

    /// Hall coefficient R_H in cm^3/C.
    ///
    /// For a single carrier type (n-type): R_H = -1 / (q * n)
    /// For a single carrier type (p-type): R_H = 1 / (q * p)
    pub fn hall_coefficient_n(n: f64) -> f64 {
        -1.0 / (Q_ELECTRON * n)
    }

    /// Hall coefficient for p-type.
    pub fn hall_coefficient_p(p: f64) -> f64 {
        1.0 / (Q_ELECTRON * p)
    }

    /// Hall coefficient for a mixed carrier system.
    ///
    /// R_H = (p * mu_p^2 - n * mu_n^2) / (q * (p * mu_p + n * mu_n)^2)
    pub fn hall_coefficient_mixed(n: f64, mu_n: f64, p: f64, mu_p: f64) -> f64 {
        let numerator = p * mu_p * mu_p - n * mu_n * mu_n;
        let denominator = Q_ELECTRON * (p * mu_p + n * mu_n).powi(2);
        if denominator.abs() < 1e-50 {
            return 0.0;
        }
        numerator / denominator
    }

    /// Hall mobility mu_H = R_H * sigma in cm^2/(V·s).
    pub fn hall_mobility(hall_coefficient: f64, conductivity: f64) -> f64 {
        (hall_coefficient * conductivity).abs()
    }

    /// Carrier concentration from Hall measurement in cm^-3.
    pub fn carrier_from_hall(hall_coefficient: f64) -> f64 {
        1.0 / (Q_ELECTRON * hall_coefficient.abs())
    }

    /// Four-point probe resistivity measurement.
    ///
    /// For a thin sheet: rho = (pi / ln(2)) * t * V / I * correction_factor
    pub fn four_point_probe(
        voltage: f64,
        current: f64,
        thickness: f64,
        correction_factor: f64,
    ) -> f64 {
        (PI / (2.0_f64).ln()) * thickness * (voltage / current) * correction_factor
    }

    /// Mobility from conductivity and carrier concentration.
    pub fn mobility_from_conductivity(conductivity: f64, carrier_conc: f64) -> f64 {
        conductivity / (Q_ELECTRON * carrier_conc)
    }

    /// Contact resistance for ohmic contacts using TLM model.
    ///
    /// R_total = 2 * R_c + R_sheet * d / W
    pub fn tlm_contact_resistance(
        r_total: f64,
        r_sheet: f64,
        gap_distance: f64,
        width: f64,
    ) -> f64 {
        (r_total - r_sheet * gap_distance / width) / 2.0
    }

    /// Specific contact resistivity from contact resistance.
    pub fn specific_contact_resistivity(r_contact: f64, contact_length: f64, width: f64) -> f64 {
        r_contact * contact_length * width
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn test_thermal_voltage_300k() {
        let vt = thermal_voltage(300.0);
        // ~25.85 mV at 300K
        assert!((vt - 0.02585).abs() < 1e-4);
    }

    #[test]
    fn test_silicon_band_gap() {
        let si = BandStructure::silicon();
        assert!((si.band_gap_ev - 1.12).abs() < EPS);
    }

    #[test]
    fn test_band_edges_consistency() {
        let si = BandStructure::silicon();
        let gap = si.ec - si.ev;
        assert!((gap - si.band_gap_ev).abs() < EPS);
    }

    #[test]
    fn test_nc_nv_positive() {
        let si = BandStructure::silicon();
        assert!(si.nc(300.0) > 0.0);
        assert!(si.nv(300.0) > 0.0);
    }

    #[test]
    fn test_silicon_nc_order_of_magnitude() {
        let si = BandStructure::silicon();
        let nc = si.nc(300.0);
        // Nc for Si ~ 2.8e19 cm^-3
        assert!(nc > 1e18, "Nc too small: {nc}");
        assert!(nc < 1e21, "Nc too large: {nc}");
    }

    #[test]
    fn test_intrinsic_concentration_silicon_300k() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 0.0, 0.0);
        let ni = cc.intrinsic_concentration();
        // ni for Si at 300K ~ 1e10 cm^-3
        assert!(ni > 1e8, "ni too small: {ni}");
        assert!(ni < 1e12, "ni too large: {ni}");
    }

    #[test]
    fn test_intrinsic_np_product() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 0.0, 0.0);
        let n = cc.electron_concentration();
        let p = cc.hole_concentration();
        let ni = cc.intrinsic_concentration();
        // n * p should equal ni^2 for intrinsic
        let ratio = (n * p) / (ni * ni);
        assert!((ratio - 1.0).abs() < 0.01, "np != ni^2, ratio = {ratio}");
    }

    #[test]
    fn test_n_type_doping() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 1e16, 0.0);
        let n = cc.electron_concentration();
        let p = cc.hole_concentration();
        assert!(n > p, "n-type should have n >> p");
        assert!(
            (n - 1e16).abs() / 1e16 < 0.01,
            "n should ≈ Nd for heavy doping"
        );
    }

    #[test]
    fn test_p_type_doping() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 0.0, 1e17);
        let n = cc.electron_concentration();
        let p = cc.hole_concentration();
        assert!(p > n, "p-type should have p >> n");
    }

    #[test]
    fn test_doping_type_classification() {
        let si = BandStructure::silicon();
        let cc_n = CarrierConcentration::new(si, 300.0, 1e16, 0.0);
        assert_eq!(cc_n.doping_type(), DopingType::NType);

        let cc_p = CarrierConcentration::new(si, 300.0, 0.0, 1e16);
        assert_eq!(cc_p.doping_type(), DopingType::PType);

        let cc_i = CarrierConcentration::new(si, 300.0, 0.0, 0.0);
        assert_eq!(cc_i.doping_type(), DopingType::Intrinsic);
    }

    #[test]
    fn test_compensation_ratio() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 1e16, 5e15);
        let kr = cc.compensation_ratio();
        assert!((kr - 0.5).abs() < EPS);
    }

    #[test]
    fn test_fermi_level_shift_n_type() {
        let si = BandStructure::silicon();
        let cc = CarrierConcentration::new(si, 300.0, 1e16, 0.0);
        let shift = cc.fermi_level_shift();
        // For n-type, Ef should be above Ei
        assert!(shift > 0.0, "Fermi level should be above Ei for n-type");
    }

    #[test]
    fn test_einstein_relation() {
        let dd = DriftDiffusion::silicon_300k();
        let dn = dd.diffusion_coefficient_n();
        let ratio = dn / dd.mu_n;
        let vt = thermal_voltage(300.0);
        assert!(
            (ratio - vt).abs() < 1e-5,
            "Einstein relation violated: D/mu = {ratio}, V_T = {vt}"
        );
    }

    #[test]
    fn test_drift_velocity_saturation() {
        let dd = DriftDiffusion::silicon_300k();
        // At very high field, velocity should approach v_sat
        let v = dd.electron_drift_velocity(1e6);
        assert!(
            v < dd.v_sat_n * 1.01,
            "Velocity exceeds saturation: {v} > {}",
            dd.v_sat_n
        );
        assert!(v > dd.v_sat_n * 0.5, "Velocity too low at high field: {v}");
    }

    #[test]
    fn test_drift_velocity_low_field() {
        let dd = DriftDiffusion::silicon_300k();
        // At low field, velocity ≈ mu * E
        let e_field = 100.0; // V/cm, low field
        let v = dd.electron_drift_velocity(e_field);
        let v_expected = dd.mu_n * e_field;
        assert!(
            (v - v_expected).abs() / v_expected < 0.1,
            "Low-field velocity: {v}, expected: {v_expected}"
        );
    }

    #[test]
    fn test_pn_junction_built_in_potential() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e16, 300.0, 1e-4);
        let vbi = pn.built_in_potential();
        // Typical Vbi for Si ≈ 0.6-0.8V
        assert!(vbi > 0.4, "Vbi too low: {vbi}");
        assert!(vbi < 1.0, "Vbi too high: {vbi}");
    }

    #[test]
    fn test_pn_junction_depletion_width() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e16, 300.0, 1e-4);
        let w = pn.depletion_width(0.0);
        // Depletion width should be > 0
        assert!(w > 0.0, "Depletion width should be positive");
        // Typical: sub-micron to few microns
        assert!(w < 1.0, "Depletion width too large: {w} cm");

        // Reverse bias should increase depletion width
        let w_rev = pn.depletion_width(-5.0);
        assert!(w_rev > w, "Reverse bias should widen depletion region");
    }

    #[test]
    fn test_pn_junction_xn_xp_sum() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e18, 300.0, 1e-4);
        let w = pn.depletion_width(0.0);
        let xn = pn.xn(0.0);
        let xp = pn.xp(0.0);
        assert!(((xn + xp) - w).abs() / w < 1e-10, "xn + xp should equal W");
    }

    #[test]
    fn test_shockley_diode_forward_bias() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e16, 300.0, 1e-4);
        let i0 = 1e-12; // 1 pA saturation current
        let i_forward = pn.shockley_current(0.6, i0, 1.0);
        assert!(i_forward > 0.0, "Forward current should be positive");
        assert!(i_forward > i0 * 100.0, "Forward current should be >> I0");
    }

    #[test]
    fn test_shockley_diode_reverse_bias() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e16, 300.0, 1e-4);
        let i0 = 1e-12;
        let i_reverse = pn.shockley_current(-5.0, i0, 1.0);
        // Reverse current should be approximately -I0
        assert!(i_reverse < 0.0, "Reverse current should be negative");
        assert!((i_reverse + i0).abs() < i0 * 0.01, "Reverse current ≈ -I0");
    }

    #[test]
    fn test_junction_capacitance_reverse_bias() {
        let si = BandStructure::silicon();
        let pn = PnJunction::new(si, 1e16, 1e16, 300.0, 1e-4);
        let c0 = pn.junction_capacitance(0.0);
        let c_rev = pn.junction_capacitance(-5.0);
        // Capacitance should decrease with reverse bias
        assert!(c_rev < c0, "Capacitance should decrease with reverse bias");
    }

    #[test]
    fn test_schottky_barrier_height() {
        let si = BandStructure::silicon();
        // Aluminum on n-Si: phi_M = 4.08 eV, chi = 4.05 eV
        let sb = SchottkyBarrier::new(4.08, si, 1e16, DopingType::NType, 300.0, 1.0, 1e-4);
        let phi_b = sb.barrier_height();
        assert!((phi_b - 0.03).abs() < 0.01, "Barrier height: {phi_b}");
    }

    #[test]
    fn test_schottky_current_forward() {
        let si = BandStructure::silicon();
        let sb = SchottkyBarrier::new(4.75, si, 1e16, DopingType::NType, 300.0, 1.05, 1e-4);
        let i_fwd = sb.current(0.3);
        assert!(i_fwd > 0.0, "Forward Schottky current should be positive");
    }

    #[test]
    fn test_srh_recombination_equilibrium() {
        let ni = 1e10;
        let rm = RecombinationModel::silicon_300k(ni);
        // At equilibrium (n=ni, p=ni), recombination rate should be 0
        let rate = rm.srh_rate(ni, ni);
        assert!(
            rate.abs() < 1.0,
            "SRH rate at equilibrium should ≈ 0: {rate}"
        );
    }

    #[test]
    fn test_total_recombination_positive_injection() {
        let ni = 1e10;
        let rm = RecombinationModel::silicon_300k(ni);
        // Excess carriers: n = 1e16, p = 1e16
        let rate = rm.total_rate(1e16, 1e16);
        assert!(
            rate > 0.0,
            "Recombination rate should be positive for injection"
        );
    }

    #[test]
    fn test_mosfet_cutoff() {
        let mos = MosfetModel::typical_nmos();
        let id = mos.drain_current(0.0, 1.0);
        // In subthreshold, current should be very small
        assert!(id < 1e-6, "Cutoff current should be very small: {id}");
    }

    #[test]
    fn test_mosfet_saturation_current() {
        let mos = MosfetModel::typical_nmos();
        let id_sat = mos.saturation_current(1.5, 2.0);
        assert!(id_sat > 0.0, "Saturation current should be positive");
    }

    #[test]
    fn test_mosfet_linear_vs_saturation() {
        let mos = MosfetModel::typical_nmos();
        let vgs = 1.5;
        let vds_linear = 0.2;
        let vds_sat = 2.0;
        let id_lin = mos.drain_current(vgs, vds_linear);
        let id_sat = mos.drain_current(vgs, vds_sat);
        // In saturation, current should be higher but not proportionally
        assert!(id_sat > id_lin, "Saturation current > linear current");
    }

    #[test]
    fn test_subthreshold_swing() {
        let mos = MosfetModel::typical_nmos();
        let ss = mos.subthreshold_swing();
        // Ideal SS ≈ 60 mV/dec at 300K; with n_sub=1.3, ≈ 78 mV/dec
        assert!(ss > 55.0, "SS too low: {ss} mV/dec");
        assert!(ss < 120.0, "SS too high: {ss} mV/dec");
    }

    #[test]
    fn test_mosfet_transconductance() {
        let mos = MosfetModel::typical_nmos();
        let gm = mos.transconductance(1.5);
        assert!(gm > 0.0, "Transconductance should be positive above Vth");
        let gm_cutoff = mos.transconductance(0.3);
        assert_eq!(gm_cutoff, 0.0, "Transconductance should be 0 below Vth");
    }

    #[test]
    fn test_resistivity_calculation() {
        // n-type Si: n=1e16, mu_n=1400, p negligible
        let rho = SemiconductorAnalysis::resistivity(1e16, 1400.0, 0.0, 0.0);
        // rho = 1 / (q * n * mu_n) ≈ 0.446 Ohm·cm
        assert!(rho > 0.1, "Resistivity too low: {rho}");
        assert!(rho < 1.0, "Resistivity too high: {rho}");
    }

    #[test]
    fn test_sheet_resistance() {
        let rho = 0.5; // Ohm·cm
        let thickness = 1e-4; // 1 um
        let rs = SemiconductorAnalysis::sheet_resistance(rho, thickness);
        assert_eq!(rs, 5000.0);
    }

    #[test]
    fn test_hall_coefficient_n_type() {
        let rh = SemiconductorAnalysis::hall_coefficient_n(1e16);
        assert!(rh < 0.0, "Hall coefficient for n-type should be negative");
    }

    #[test]
    fn test_carrier_from_hall() {
        let rh = SemiconductorAnalysis::hall_coefficient_n(1e16);
        let n_measured = SemiconductorAnalysis::carrier_from_hall(rh);
        assert!(
            (n_measured - 1e16).abs() / 1e16 < 0.01,
            "Carrier from Hall: {n_measured}"
        );
    }

    #[test]
    fn test_body_effect_increases_vth() {
        let vth0 = 0.5;
        let gamma = 0.4;
        let phi_f = 0.35;
        let vth_with_body = MosfetModel::body_effect_vth(vth0, gamma, phi_f, 1.0);
        assert!(
            vth_with_body > vth0,
            "Body effect should increase Vth: {} vs {}",
            vth_with_body,
            vth0
        );
    }

    #[test]
    fn test_gaas_higher_mobility() {
        let si = DriftDiffusion::silicon_300k();
        let gaas = DriftDiffusion::gaas_300k();
        assert!(
            gaas.mu_n > si.mu_n,
            "GaAs should have higher electron mobility than Si"
        );
    }

    #[test]
    fn test_varshni_band_gap() {
        let si = BandStructure::silicon();
        // Varshni parameters for Si: alpha=4.73e-4, beta=636
        let eg_0k = si.varshni_band_gap(0.0, 4.73e-4, 636.0);
        let eg_300k = si.varshni_band_gap(300.0, 4.73e-4, 636.0);
        // Eg should increase as T decreases (for most semiconductors)
        assert!(
            eg_0k > eg_300k,
            "Band gap should be larger at 0K: {eg_0k} vs {eg_300k}"
        );
    }
}
