// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Combustion models: species grids, reactor models, flame tracking,
//! Zeldovich NOx, chain branching, ignition delay, laminar flame speed,
//! and multi-species diffusion.

use super::lbm_reactive::arrhenius_rate;
use super::species::{ArrheniusRate, R_GAS};

// ---------------------------------------------------------------------------
// D2Q9 constants (same ordering as entropic.rs)
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights used for passive scalar equilibrium.
const W: [f64; 9] = [
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

/// D2Q9 lattice velocities: (cx, cy) for each direction.
const C: [(i32, i32); 9] = [
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

/// Speed of sound squared for D2Q9: cs^2 = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// SpeciesGrid — multi-species LBM density arrays
// ---------------------------------------------------------------------------

/// Per-cell density arrays for multiple chemical species on a 2-D D2Q9 grid.
///
/// Each species `s` carries:
/// - a full set of 9 distribution functions `g[s][cell][dir]` for LBM transport.
/// - a macroscopic concentration field derived from `sum_i g[s][cell][i]`.
///
/// The equilibrium for each species scalar distribution is:
///
/// ```text
/// g_eq_i(s) = w_i * c_s * (1 + (e_i · u) / cs²)
/// ```
///
/// with BGK relaxation time `tau_s = D_s / cs² + 0.5`.
pub struct SpeciesGrid {
    /// Domain width (cells in x).
    pub nx: usize,
    /// Domain height (cells in y).
    pub ny: usize,
    /// Number of chemical species.
    pub n_species: usize,
    /// Molecular diffusivities, one per species.
    pub diffusivities: Vec<f64>,
    /// Distribution functions: `g[species][cell][direction]`.
    pub g: Vec<Vec<[f64; 9]>>,
}

impl SpeciesGrid {
    /// Create a new `SpeciesGrid` with all concentrations zero.
    ///
    /// `diffusivities` must have length `n_species`.
    pub fn new(nx: usize, ny: usize, diffusivities: Vec<f64>) -> Self {
        let n_species = diffusivities.len();
        let n = nx * ny;
        let g = vec![vec![[0.0_f64; 9]; n]; n_species];
        Self {
            nx,
            ny,
            n_species,
            diffusivities,
            g,
        }
    }

    /// Compute scalar equilibrium distribution for species `s` at concentration `c`
    /// and fluid velocity `[ux, uy]`.
    ///
    /// `g_eq_i = w_i * c * (1 + (e_i · u) / cs²)`
    pub fn equilibrium(c: f64, ux: f64, uy: f64) -> [f64; 9] {
        let mut geq = [0.0_f64; 9];
        for (geq_i, (&c_vel, &w)) in geq.iter_mut().zip(C.iter().zip(W.iter())) {
            let cx = c_vel.0 as f64;
            let cy = c_vel.1 as f64;
            let eu = cx * ux + cy * uy;
            *geq_i = w * c * (1.0 + eu / CS2);
        }
        geq
    }

    /// Return the macroscopic concentration of species `s` at cell index `k`.
    pub fn concentration(&self, s: usize, k: usize) -> f64 {
        self.g[s][k].iter().sum()
    }

    /// Set species `s` at cell `k` to concentration `c` at zero velocity.
    pub fn set_concentration(&mut self, s: usize, k: usize, c: f64) {
        self.g[s][k] = Self::equilibrium(c, 0.0, 0.0);
    }

    /// Perform BGK collision for all species using the supplied per-cell velocities.
    ///
    /// `velocities[k] = [ux, uy]` for each cell `k`.
    pub fn collide(&mut self, velocities: &[[f64; 2]]) {
        for s in 0..self.n_species {
            let omega_s = 1.0 / (self.diffusivities[s] / CS2 + 0.5);
            for (k, g_sk) in self.g[s].iter_mut().enumerate() {
                let c: f64 = g_sk.iter().sum();
                let ux = velocities[k][0];
                let uy = velocities[k][1];
                let geq = Self::equilibrium(c, ux, uy);
                for (i, g_ski) in g_sk.iter_mut().enumerate() {
                    *g_ski -= omega_s * (*g_ski - geq[i]);
                }
            }
        }
    }

    /// Pull-scheme periodic streaming for all species.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for s in 0..self.n_species {
            let g_old = self.g[s].clone();
            for y in 0..ny {
                for x in 0..nx {
                    let k_dst = y * nx + x;
                    for i in 0..9 {
                        let cx = C[i].0;
                        let cy = C[i].1;
                        let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                        let k_src = sy * nx + sx;
                        self.g[s][k_dst][i] = g_old[k_src][i];
                    }
                }
            }
        }
    }

    /// Apply a pointwise reaction source term: concentration of species `s` changes by
    /// `delta_c[k]` at each cell `k`.  Distributed uniformly across all directions.
    pub fn apply_source(&mut self, s: usize, delta_c: &[f64]) {
        let n = self.nx * self.ny;
        assert_eq!(delta_c.len(), n);
        for (k, g_sk) in self.g[s].iter_mut().enumerate() {
            for (i, g_ski) in g_sk.iter_mut().enumerate() {
                *g_ski += W[i] * delta_c[k];
            }
        }
    }

    /// Return total concentration of species `s` across all cells.
    pub fn total_concentration(&self, s: usize) -> f64 {
        let n = self.nx * self.ny;
        (0..n).map(|k| self.concentration(s, k)).sum()
    }
}

// ---------------------------------------------------------------------------
// Single-step combustion: CH4 + 2 O2 → CO2 + 2 H2O
// ---------------------------------------------------------------------------

/// Species indices for the CH4 combustion system.
///
/// Convention used by `CombustionReactor`:
/// 0 = CH4, 1 = O2, 2 = CO2, 3 = H2O
pub mod combustion_indices {
    /// Species index for methane (CH₄).
    pub const CH4: usize = 0;
    /// Species index for oxygen (O₂).
    pub const O2: usize = 1;
    /// Species index for carbon dioxide (CO₂).
    pub const CO2: usize = 2;
    /// Species index for water vapour (H₂O).
    pub const H2O: usize = 3;
    /// Total number of species tracked.
    pub const N_SPECIES: usize = 4;
}

/// Stoichiometric coefficients for CH4 + 2 O2 → CO2 + 2 H2O.
///
/// Negative = consumed, positive = produced.  Units: moles per mole CH4 consumed.
pub const COMBUSTION_STOICH: [f64; 4] = [-1.0, -2.0, 1.0, 2.0];

/// Standard heat of combustion of CH4: ΔH_rxn = –890 kJ/mol (exothermic).
///
/// In lattice units this is scaled to a dimensionless value by the user.
pub const DELTA_H_COMBUSTION: f64 = -890_000.0; // J/mol (physical)

/// Single-step combustion reactor: CH4 + 2 O2 → CO2 + 2 H2O.
///
/// Uses an Arrhenius rate law with temperature feedback:
///
/// ```text
/// r = A * exp(-Ea / (R * T)) * [CH4] * [O2]²
/// ```
///
/// where concentrations are mole fractions or molar concentrations.
/// Temperature is updated each step by the exothermic heat release.
pub struct CombustionReactor {
    /// Number of grid cells.
    pub n_cells: usize,
    /// Species concentrations: `c[species][cell]`, 4 species (CH4, O2, CO2, H2O).
    pub c: [[Vec<f64>; 4]; 1],
    /// Temperature field (K) at each cell.
    pub temperature: Vec<f64>,
    /// Arrhenius pre-exponential factor A (1/s).
    pub pre_exp: f64,
    /// Activation energy Ea (J/mol).
    pub activation_energy: f64,
    /// Heat of reaction (J/mol) — negative = exothermic.
    pub heat_of_reaction: f64,
    /// Heat capacity (J / (mol · K)) — used to convert heat release to ΔT.
    pub heat_capacity: f64,
}

impl CombustionReactor {
    /// Create a combustion reactor for `n_cells` grid cells at initial temperature `t0`.
    ///
    /// Default Arrhenius parameters approximate methane ignition:
    /// A = 2.119 × 10¹¹ s⁻¹, Ea = 202 kJ/mol (simplified one-step).
    pub fn new(n_cells: usize, t0: f64) -> Self {
        Self {
            n_cells,
            c: [std::array::from_fn(|_| vec![0.0_f64; n_cells])],
            temperature: vec![t0; n_cells],
            pre_exp: 2.119e11,
            activation_energy: 202_000.0,
            heat_of_reaction: DELTA_H_COMBUSTION,
            heat_capacity: 1000.0, // J / (mol · K), approximate mixture Cp
        }
    }

    /// Get the concentration of species `s` at cell `k`.
    pub fn get(&self, s: usize, k: usize) -> f64 {
        self.c[0][s][k]
    }

    /// Set the concentration of species `s` at cell `k`.
    pub fn set(&mut self, s: usize, k: usize, v: f64) {
        self.c[0][s][k] = v;
    }

    /// Advance one reaction step using explicit Euler.
    ///
    /// Rate: `r = A * exp(-Ea/(R*T)) * [CH4] * [O2]²`
    ///
    /// Species update: `ΔC_s = stoich_s * r_eff * dt`
    /// Temperature update: `ΔT = -ΔH_rxn * r_eff * dt / Cp`
    ///
    /// The effective rate `r_eff` is limited so that no reactant concentration
    /// goes negative, preserving exact stoichiometry.
    pub fn step(&mut self, dt: f64) {
        for k in 0..self.n_cells {
            let t = self.temperature[k].max(1.0); // guard against T=0
            let k_arr = self.pre_exp * (-self.activation_energy / (R_GAS * t)).exp();
            let c_ch4 = self.c[0][combustion_indices::CH4][k].max(0.0);
            let c_o2 = self.c[0][combustion_indices::O2][k].max(0.0);
            let rate = k_arr * c_ch4 * c_o2 * c_o2;
            // Cap the rate so that no reactant goes negative:
            //   CH4 consumed = 1 * rate * dt <= c_ch4  => rate <= c_ch4 / dt
            //   O2  consumed = 2 * rate * dt <= c_o2   => rate <= c_o2 / (2*dt)
            let max_rate_ch4 = if dt > 1e-30 {
                c_ch4 / dt
            } else {
                f64::INFINITY
            };
            let max_rate_o2 = if dt > 1e-30 {
                c_o2 / (2.0 * dt)
            } else {
                f64::INFINITY
            };
            let r_eff = rate.min(max_rate_ch4).min(max_rate_o2);
            // Update each species using the stoichiometrically consistent rate.
            for (s, stoich) in COMBUSTION_STOICH.iter().enumerate() {
                let delta = stoich * r_eff * dt;
                self.c[0][s][k] = (self.c[0][s][k] + delta).max(0.0);
            }
            // Exothermic heat release: ΔT = -ΔH * r_eff * dt / Cp
            let delta_t = -self.heat_of_reaction * r_eff * dt / self.heat_capacity;
            self.temperature[k] += delta_t;
        }
    }

    /// Total molar concentration of all species at cell `k` (should be conserved).
    pub fn total_moles_at(&self, k: usize) -> f64 {
        (0..4).map(|s| self.c[0][s][k]).sum()
    }
}

// ---------------------------------------------------------------------------
// Temperature-field coupling (source/sink from reactions)
// ---------------------------------------------------------------------------

/// Apply an exothermic heat-release source to a temperature field.
///
/// For each cell `k`:
/// ```text
/// T[k] += delta_h * reaction_rate[k] * dt
/// ```
///
/// `delta_h` is positive for exothermic reactions (releases heat).
pub fn apply_heat_source(temperature: &mut [f64], reaction_rates: &[f64], delta_h: f64, dt: f64) {
    debug_assert_eq!(temperature.len(), reaction_rates.len());
    for (t, &r) in temperature.iter_mut().zip(reaction_rates.iter()) {
        *t += delta_h * r * dt;
    }
}

/// Compute species source/sink terms from a single-step reaction.
///
/// Returns a `Vec<Vec`f64`>` of shape `[n_species][n_cells]` containing
/// `stoich[s] * rate[k]` for each cell.
///
/// `rate_field[k]` is the scalar reaction rate at cell `k`.
/// `stoichiometry` has length `n_species`.
pub fn species_source_terms(stoichiometry: &[f64], rate_field: &[f64]) -> Vec<Vec<f64>> {
    let n_cells = rate_field.len();
    let n_species = stoichiometry.len();
    let mut sources = vec![vec![0.0_f64; n_cells]; n_species];
    for s in 0..n_species {
        for k in 0..n_cells {
            sources[s][k] = stoichiometry[s] * rate_field[k];
        }
    }
    sources
}

// ---------------------------------------------------------------------------
// Flame front tracking via temperature gradient
// ---------------------------------------------------------------------------

/// Track the flame front in a 2-D temperature field.
///
/// The flame front is defined as the location of maximum temperature gradient
/// magnitude `|∇T|` along horizontal slices.
pub struct FlameFrontTracker {
    /// Threshold for `|∇T|` to qualify as "on the flame front".
    pub grad_threshold: f64,
    /// History: `(time, average_x_position)`.
    pub history: Vec<(f64, f64)>,
}

impl FlameFrontTracker {
    /// Create a new flame-front tracker.
    pub fn new(grad_threshold: f64) -> Self {
        Self {
            grad_threshold,
            history: Vec::new(),
        }
    }

    /// Scan the temperature field for the flame front and record its position.
    ///
    /// Uses a central-difference gradient in x; the front position is the
    /// average x of cells where `|dT/dx| >= threshold`.
    pub fn track(&mut self, temperature: &[f64], nx: usize, ny: usize, time: f64) {
        let mut sum_x = 0.0_f64;
        let mut count = 0_usize;
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1).min(nx - 1);
                let im = if i == 0 { 0 } else { i - 1 };
                let grad = (temperature[j * nx + ip] - temperature[j * nx + im]).abs() / 2.0;
                if grad >= self.grad_threshold {
                    sum_x += i as f64;
                    count += 1;
                }
            }
        }
        let avg_x = if count > 0 {
            sum_x / count as f64
        } else {
            f64::NAN
        };
        self.history.push((time, avg_x));
    }

    /// Number of recorded front positions.
    pub fn num_records(&self) -> usize {
        self.history.len()
    }

    /// Estimated flame speed from the last two records.
    pub fn flame_speed(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        let n = self.history.len();
        let (t1, x1) = self.history[n - 2];
        let (t2, x2) = self.history[n - 1];
        let dt = t2 - t1;
        if dt.abs() < 1e-20 {
            return 0.0;
        }
        (x2 - x1) / dt
    }
}

// ---------------------------------------------------------------------------
// Zeldovich reaction mechanism
// ---------------------------------------------------------------------------

/// Zeldovich thermal NO formation mechanism.
///
/// The extended Zeldovich mechanism consists of two elementary reactions:
///   R1: O + N2 ⇌ NO + N   (Zeldovich reaction 1)
///   R2: N + O2 ⇌ NO + O   (Zeldovich reaction 2)
///
/// Together they give the net reaction: N2 + O2 → 2 NO
///
/// Rate constants follow Arrhenius kinetics with temperature-dependent
/// pre-exponentials from the GRI mechanism.
pub struct ZeldovichMechanism {
    /// Arrhenius pre-exponential for R1 forward (cm^3 mol^-1 s^-1).
    pub a1f: f64,
    /// Activation energy for R1 forward (J/mol).
    pub ea1f: f64,
    /// Arrhenius pre-exponential for R1 reverse.
    pub a1r: f64,
    /// Activation energy for R1 reverse (J/mol).
    pub ea1r: f64,
    /// Arrhenius pre-exponential for R2 forward.
    pub a2f: f64,
    /// Activation energy for R2 forward (J/mol).
    pub ea2f: f64,
    /// Arrhenius pre-exponential for R2 reverse.
    pub a2r: f64,
    /// Activation energy for R2 reverse (J/mol).
    pub ea2r: f64,
}

impl ZeldovichMechanism {
    /// Create a Zeldovich mechanism with default rate constants
    /// (simplified air chemistry, Baulch et al. 1994).
    pub fn default_air() -> Self {
        Self {
            a1f: 1.8e14,     // cm^3 mol^-1 s^-1
            ea1f: 318_000.0, // J/mol
            a1r: 3.8e13,
            ea1r: 0.0,
            a2f: 1.8e10,
            ea2f: 20_000.0,
            a2r: 3.8e9,
            ea2r: 167_000.0,
        }
    }

    /// Forward rate constant for reaction 1: O + N2 → NO + N.
    pub fn k1f(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.a1f, self.ea1f, temperature, R_GAS)
    }

    /// Reverse rate constant for reaction 1: NO + N → O + N2.
    pub fn k1r(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.a1r, self.ea1r, temperature, R_GAS)
    }

    /// Forward rate constant for reaction 2: N + O2 → NO + O.
    pub fn k2f(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.a2f, self.ea2f, temperature, R_GAS)
    }

    /// Reverse rate constant for reaction 2: NO + O → N + O2.
    pub fn k2r(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.a2r, self.ea2r, temperature, R_GAS)
    }

    /// Net NO production rate (mol m^-3 s^-1) in quasi-steady-state approximation.
    ///
    /// d\[NO\]/dt ≈ 2 * k1f * \[O\] * \[N2\] * (1 - \[NO\]^2 / (K_eq * \[N2\] * \[O2\]))
    ///             / (1 + k1r * \[NO\] / (k2f * \[O2\] + k2r * \[NO\]))
    ///
    /// Concentrations in mol/m^3.
    pub fn no_production_rate(
        &self,
        temperature: f64,
        c_o: f64,
        c_n2: f64,
        c_o2: f64,
        c_no: f64,
    ) -> f64 {
        let k1f = self.k1f(temperature);
        let k1r = self.k1r(temperature);
        let k2f = self.k2f(temperature);
        let k2r = self.k2r(temperature);
        // Equilibrium constant for overall N2 + O2 → 2NO
        let k_eq_inv = if k1r * k2r > 1e-100 {
            (k1f * k2f) / (k1r * k2r)
        } else {
            1e-30
        };
        let no_sq = c_no * c_no;
        let n2_o2 = c_n2 * c_o2;
        let equilibrium_factor = if k_eq_inv * n2_o2 > 1e-30 {
            1.0 - no_sq / (k_eq_inv * n2_o2 + 1e-100)
        } else {
            1.0
        };
        let denominator = 1.0 + k1r * c_no / (k2f * c_o2 + k2r * c_no + 1e-100);
        2.0 * k1f * c_o * c_n2 * equilibrium_factor / denominator.max(1e-100)
    }

    /// Equilibrium NO mole fraction at a given temperature.
    ///
    /// For N2/O2 mixture: \[NO\]_eq = sqrt(K_eq * \[N2\] * \[O2\]).
    pub fn equilibrium_no(&self, temperature: f64, c_n2: f64, c_o2: f64) -> f64 {
        let k1f = self.k1f(temperature);
        let k1r = self.k1r(temperature).max(1e-100);
        let k2f = self.k2f(temperature);
        let k2r = self.k2r(temperature).max(1e-100);
        let k_eq = (k1f * k2f) / (k1r * k2r);
        (k_eq * c_n2 * c_o2).max(0.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Chain branching mechanism
// ---------------------------------------------------------------------------

/// Chain branching reaction for hydrogen oxidation.
///
/// Simplified 3-step chain mechanism:
///   Initiation:  H2 + O2 → 2 OH   (rate k_i)
///   Branching:   OH + H2 → H2O + H (rate k_b)
///   Termination: H + wall → ½ H2   (rate k_t)
///
/// Net chain-branching factor: phi = k_b - k_t  (explosion if phi > 0).
#[derive(Debug, Clone, Copy)]
pub struct ChainBranching {
    /// Initiation rate constant k_i (mol^-1 m^3 s^-1).
    pub k_init: f64,
    /// Branching rate constant k_b (mol^-1 m^3 s^-1).
    pub k_branch: f64,
    /// Termination rate constant k_term (s^-1, first order wall).
    pub k_term: f64,
    /// Initiation activation energy (J/mol).
    pub ea_init: f64,
    /// Branching activation energy (J/mol).
    pub ea_branch: f64,
}

impl ChainBranching {
    /// Create a chain branching model with given rate constants.
    pub fn new(k_init: f64, ea_init: f64, k_branch: f64, ea_branch: f64, k_term: f64) -> Self {
        Self {
            k_init,
            ea_init,
            k_branch,
            ea_branch,
            k_term,
        }
    }

    /// Temperature-dependent initiation rate (mol^-1 m^3 s^-1).
    pub fn initiation_rate(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.k_init, self.ea_init, temperature, R_GAS)
    }

    /// Temperature-dependent branching rate (mol^-1 m^3 s^-1).
    pub fn branching_rate(&self, temperature: f64) -> f64 {
        arrhenius_rate(self.k_branch, self.ea_branch, temperature, R_GAS)
    }

    /// Net branching factor: phi = k_b * \[H2\] - k_t.
    ///
    /// If phi > 0: chain explosion; phi < 0: chain termination dominates.
    pub fn net_branching_factor(&self, temperature: f64, c_h2: f64) -> f64 {
        self.branching_rate(temperature) * c_h2 - self.k_term
    }

    /// Explosion limit: temperature above which phi > 0 at given \[H2\].
    ///
    /// Finds T_exp by bisection in \[T_low, T_high\].
    pub fn explosion_temperature(&self, c_h2: f64, t_low: f64, t_high: f64) -> f64 {
        let mut lo = t_low;
        let mut hi = t_high;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.net_branching_factor(mid, c_h2) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// OH radical concentration in quasi-steady state.
    ///
    /// \[OH\]_ss = k_i * \[H2\] * \[O2\] / (k_t - k_b * \[H2\])
    /// (only valid below explosion limit)
    pub fn quasi_steady_oh(&self, temperature: f64, c_h2: f64, c_o2: f64) -> Option<f64> {
        let k_i = self.initiation_rate(temperature);
        let denom = self.k_term - self.branching_rate(temperature) * c_h2;
        if denom > 0.0 {
            Some(k_i * c_h2 * c_o2 / denom)
        } else {
            None // explosion regime
        }
    }
}

// ---------------------------------------------------------------------------
// Ignition delay model
// ---------------------------------------------------------------------------

/// Ignition delay time model using Arrhenius correlation.
///
/// tau_ign = A * \[fuel\]^a * \[O2\]^b * exp(Ea / (R * T))
///
/// where tau_ign is the induction time (s).
#[derive(Debug, Clone, Copy)]
pub struct IgnitionDelayModel {
    /// Pre-exponential A (units depend on a, b).
    pub pre_exp: f64,
    /// Fuel concentration exponent a.
    pub fuel_exponent: f64,
    /// Oxygen concentration exponent b.
    pub o2_exponent: f64,
    /// Activation energy Ea (J/mol).
    pub activation_energy: f64,
}

impl IgnitionDelayModel {
    /// Create a new ignition delay model.
    pub fn new(pre_exp: f64, fuel_exponent: f64, o2_exponent: f64, activation_energy: f64) -> Self {
        Self {
            pre_exp,
            fuel_exponent,
            o2_exponent,
            activation_energy,
        }
    }

    /// Compute ignition delay time (s).
    ///
    /// tau = A * \[fuel\]^a * \[O2\]^b * exp(Ea / (R * T))
    pub fn delay_time(&self, temperature: f64, c_fuel: f64, c_o2: f64) -> f64 {
        let exponential = (self.activation_energy / (R_GAS * temperature)).exp();
        self.pre_exp * c_fuel.powf(self.fuel_exponent) * c_o2.powf(self.o2_exponent) * exponential
    }

    /// Effective activation energy from two ignition delay measurements.
    ///
    /// E_a = R * T1 * T2 / (T2 - T1) * ln(tau1 / tau2)
    pub fn effective_activation_energy(&self, t1: f64, t2: f64, tau1: f64, tau2: f64) -> f64 {
        if (t2 - t1).abs() < 1e-10 {
            return self.activation_energy;
        }
        R_GAS * t1 * t2 / (t2 - t1) * (tau1 / tau2).ln()
    }

    /// Dimensionless Damköhler number for ignition:
    ///
    /// Da_ign = tau_flow / tau_ign
    pub fn damkohler_number(&self, tau_flow: f64, temperature: f64, c_fuel: f64, c_o2: f64) -> f64 {
        let tau_ign = self.delay_time(temperature, c_fuel, c_o2);
        if tau_ign < 1e-30 {
            return f64::INFINITY;
        }
        tau_flow / tau_ign
    }
}

// ---------------------------------------------------------------------------
// Laminar flame speed model (Metghalchi-Keck correlation)
// ---------------------------------------------------------------------------

/// Laminar flame speed model.
///
/// Uses the Metghalchi-Keck (1982) correlation:
///
/// S_L = S_L0 * (T/T_ref)^alpha * (P/P_ref)^beta * (1 - 2.1 * f_dil)
///
/// where alpha and beta are fuel-specific exponents.
#[derive(Debug, Clone, Copy)]
pub struct LaminarFlameSpeed {
    /// Reference flame speed S_L0 (m/s) at stoichiometry, T_ref, P_ref.
    pub s_l0: f64,
    /// Temperature exponent alpha.
    pub alpha: f64,
    /// Pressure exponent beta.
    pub beta: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Reference pressure P_ref (Pa).
    pub p_ref: f64,
}

impl LaminarFlameSpeed {
    /// Create a flame speed model.  Default parameters for methane/air.
    pub fn methane_air() -> Self {
        Self {
            s_l0: 0.37, // m/s at stoichiometry
            alpha: 2.2,
            beta: -0.5,
            t_ref: 298.0,    // K
            p_ref: 101325.0, // Pa
        }
    }

    /// Laminar flame speed at given conditions.
    pub fn speed(&self, temperature: f64, pressure: f64, diluent_fraction: f64) -> f64 {
        let t_factor = (temperature / self.t_ref).powf(self.alpha);
        let p_factor = (pressure / self.p_ref).powf(self.beta);
        let dil_factor = 1.0 - 2.1 * diluent_fraction.clamp(0.0, 0.47);
        self.s_l0 * t_factor * p_factor * dil_factor
    }

    /// Flame speed sensitivity to temperature: dS_L/dT.
    pub fn speed_temperature_sensitivity(&self, temperature: f64, pressure: f64) -> f64 {
        let s_l = self.speed(temperature, pressure, 0.0);
        s_l * self.alpha / temperature
    }

    /// Markstein length: L_M = Zeldovich thickness × Markstein number.
    ///
    /// Simplified: L_M = D_th / S_L
    pub fn markstein_length(
        &self,
        thermal_diffusivity: f64,
        temperature: f64,
        pressure: f64,
    ) -> f64 {
        let s_l = self.speed(temperature, pressure, 0.0);
        if s_l < 1e-20 {
            return 0.0;
        }
        thermal_diffusivity / s_l
    }

    /// Flame thickness: delta_f = D_th / S_L.
    pub fn flame_thickness(
        &self,
        thermal_diffusivity: f64,
        temperature: f64,
        pressure: f64,
    ) -> f64 {
        self.markstein_length(thermal_diffusivity, temperature, pressure)
    }
}

// ---------------------------------------------------------------------------
// Elementary reaction rate in LBM (species diffusion with reaction source)
// ---------------------------------------------------------------------------

/// Elementary bimolecular reaction in LBM: A + B → C.
///
/// Uses operator splitting: diffusion (LBM passive scalar) then reaction
/// (explicit Euler source term).
pub struct ElementaryReactionLbm {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Concentration of species A: `c_a[cell]`.
    pub c_a: Vec<f64>,
    /// Concentration of species B: `c_b[cell]`.
    pub c_b: Vec<f64>,
    /// Concentration of product C: `c_c[cell]`.
    pub c_c: Vec<f64>,
    /// Rate constant k (mol^-1 m^3 s^-1 or lattice units).
    pub rate_constant: f64,
    /// Diffusivity for A (lattice units).
    pub d_a: f64,
    /// Diffusivity for B (lattice units).
    pub d_b: f64,
    /// Diffusivity for C (lattice units).
    pub d_c: f64,
}

impl ElementaryReactionLbm {
    /// Create a new elementary reaction solver.
    pub fn new(nx: usize, ny: usize, rate_constant: f64, d_a: f64, d_b: f64, d_c: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            c_a: vec![0.0; n],
            c_b: vec![0.0; n],
            c_c: vec![0.0; n],
            rate_constant,
            d_a,
            d_b,
            d_c,
        }
    }

    /// Diffusion step for a single species (explicit FD, periodic BC).
    fn diffuse(c: &mut [f64], nx: usize, ny: usize, d: f64, dt: f64, dx: f64) {
        let coeff = d * dt / (dx * dx);
        let old = c.to_vec();
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                let lap = old[j * nx + ip] + old[j * nx + im] + old[jp * nx + i] + old[jm * nx + i]
                    - 4.0 * old[j * nx + i];
                c[j * nx + i] = old[j * nx + i] + coeff * lap;
            }
        }
    }

    /// Apply reaction source terms: A + B → C.
    fn react(&mut self, dt: f64) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let delta = self.rate_constant * self.c_a[k] * self.c_b[k] * dt;
            let delta = delta.min(self.c_a[k]).min(self.c_b[k]);
            self.c_a[k] -= delta;
            self.c_b[k] -= delta;
            self.c_c[k] += delta;
        }
    }

    /// One operator-split step: diffuse A, B, C then react.
    pub fn step(&mut self, dt: f64, dx: f64) {
        let (nx, ny) = (self.nx, self.ny);
        let (d_a, d_b, d_c) = (self.d_a, self.d_b, self.d_c);
        Self::diffuse(&mut self.c_a, nx, ny, d_a, dt, dx);
        Self::diffuse(&mut self.c_b, nx, ny, d_b, dt, dx);
        Self::diffuse(&mut self.c_c, nx, ny, d_c, dt, dx);
        self.react(dt);
    }

    /// Total moles of A.
    pub fn total_a(&self) -> f64 {
        self.c_a.iter().sum()
    }
    /// Total moles of B.
    pub fn total_b(&self) -> f64 {
        self.c_b.iter().sum()
    }
    /// Total moles of C.
    pub fn total_c(&self) -> f64 {
        self.c_c.iter().sum()
    }

    /// Atom conservation check: (A + C) = const and (B + C) = const.
    pub fn a_plus_c_total(&self) -> f64 {
        self.c_a
            .iter()
            .zip(self.c_c.iter())
            .map(|(&a, &c)| a + c)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Thermal explosion (Semenov theory)
// ---------------------------------------------------------------------------

/// Semenov thermal explosion model.
///
/// Dimensionless equation:
/// dtheta/dt_tilde = exp(theta / (1 + beta * theta)) - delta * theta
///
/// where:
///   theta = dimensionless temperature excess
///   beta  = T_a / T_ad (temperature ratio, often << 1)
///   delta = heat loss / heat generation parameter
///
/// Critical condition: delta_cr = e^-1 ≈ 0.368 (classical Semenov).
pub struct SemenovExplosion;

impl SemenovExplosion {
    /// Critical Semenov parameter delta_cr for thermal explosion.
    ///
    /// Explosion occurs if delta < delta_cr.
    pub fn critical_delta() -> f64 {
        std::f64::consts::E.recip() // 1/e
    }

    /// Dimensionless induction time to explosion:
    ///
    /// t_ind ≈ ln(delta_cr / (delta_cr - delta)) / delta (near criticality)
    pub fn induction_time(delta: f64) -> Option<f64> {
        let delta_cr = Self::critical_delta();
        if delta >= delta_cr {
            return None;
        } // no explosion
        let arg = delta_cr / (delta_cr - delta);
        Some(arg.ln() / delta.max(1e-20))
    }

    /// Adiabatic temperature rise:
    ///
    /// theta_ad = Q * \[fuel\]_0 / (rho * Cp * T_0)
    pub fn adiabatic_temperature_rise(
        heat_of_reaction: f64,
        c_fuel0: f64,
        density: f64,
        cp: f64,
        t0: f64,
    ) -> f64 {
        if density * cp * t0 < 1e-30 {
            return 0.0;
        }
        heat_of_reaction * c_fuel0 / (density * cp * t0)
    }

    /// Frank-Kamenetskii number (for distributed reaction):
    ///
    /// delta_FK = Q * k(T_wall) * r^2 * E_a / (lambda * R * T_w^2)
    ///
    /// Critical value: 0.878 for sphere, 0.878 for cylinder, 0.878 for slab.
    pub fn frank_kamenetskii(
        heat_of_reaction: f64,
        k_wall: f64,
        radius: f64,
        ea: f64,
        thermal_conductivity: f64,
        t_wall: f64,
    ) -> f64 {
        let numerator = heat_of_reaction * k_wall * radius * radius * ea;
        let denominator = thermal_conductivity * R_GAS * t_wall * t_wall;
        if denominator < 1e-30 {
            return 0.0;
        }
        numerator / denominator
    }
}

// ---------------------------------------------------------------------------
// Species diffusion in LBM (multi-species passive scalar)
// ---------------------------------------------------------------------------

/// Multi-species LBM diffusion solver for reactive mixtures.
///
/// Each species has its own LBM scalar distribution `g_s[cell][dir]`
/// and is coupled through reaction source terms.
pub struct MultiSpeciesDiffusion {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Number of species.
    pub n_species: usize,
    /// Diffusivities for each species.
    pub diffusivities: Vec<f64>,
    /// Distribution functions `g[species][cell][dir]`.
    pub g: Vec<Vec<[f64; 9]>>,
    /// Reaction source terms (updated each step externally).
    pub sources: Vec<Vec<f64>>,
}

impl MultiSpeciesDiffusion {
    /// Create a new multi-species diffusion solver, initialised to zero.
    pub fn new(nx: usize, ny: usize, diffusivities: Vec<f64>) -> Self {
        let n_species = diffusivities.len();
        let n = nx * ny;
        Self {
            nx,
            ny,
            n_species,
            diffusivities,
            g: vec![vec![[0.0_f64; 9]; n]; n_species],
            sources: vec![vec![0.0; n]; n_species],
        }
    }

    /// Equilibrium scalar distribution for concentration c, velocity u.
    fn equilibrium(c: f64, u: [f64; 2]) -> [f64; 9] {
        let mut geq = [0.0_f64; 9];
        for (geq_k, (&c_vel, &w)) in geq.iter_mut().zip(C.iter().zip(W.iter())) {
            let cx = c_vel.0 as f64;
            let cy = c_vel.1 as f64;
            let eu = cx * u[0] + cy * u[1];
            *geq_k = w * c * (1.0 + eu / CS2);
        }
        geq
    }

    /// Concentration of species s at cell idx.
    pub fn concentration(&self, s: usize, idx: usize) -> f64 {
        self.g[s][idx].iter().sum()
    }

    /// Set initial concentration for species s at cell idx.
    pub fn set_concentration(&mut self, s: usize, idx: usize, c: f64) {
        self.g[s][idx] = Self::equilibrium(c, [0.0, 0.0]);
    }

    /// BGK collision for all species with Guo source terms.
    pub fn collide(&mut self, velocities: &[[f64; 2]]) {
        for s in 0..self.n_species {
            let omega_s = 1.0 / (self.diffusivities[s] / CS2 + 0.5);
            for (k, g_sk) in self.g[s].iter_mut().enumerate() {
                let c: f64 = g_sk.iter().sum();
                let u = velocities[k];
                let src_k = self.sources[s][k];
                let geq = Self::equilibrium(c, u);
                for (i, g_ski) in g_sk.iter_mut().enumerate() {
                    *g_ski -= omega_s * (*g_ski - geq[i]);
                    // Add reaction source term (Guo scheme: equally distributed).
                    *g_ski += W[i] * src_k;
                }
            }
        }
    }

    /// Pull-scheme periodic streaming for all species.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for s in 0..self.n_species {
            let g_old = self.g[s].clone();
            for y in 0..ny {
                for x in 0..nx {
                    let k_dst = y * nx + x;
                    for i in 0..9 {
                        let cx = C[i].0;
                        let cy = C[i].1;
                        let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                        let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                        self.g[s][k_dst][i] = g_old[sy * nx + sx][i];
                    }
                }
            }
        }
    }

    /// One step: collide then stream.
    pub fn step(&mut self, velocities: &[[f64; 2]]) {
        self.collide(velocities);
        self.stream();
    }

    /// Total concentration of species s across all cells.
    pub fn total_concentration(&self, s: usize) -> f64 {
        let n = self.nx * self.ny;
        (0..n).map(|k| self.concentration(s, k)).sum()
    }

    /// Set source term for species s at all cells.
    pub fn set_sources(&mut self, s: usize, sources: Vec<f64>) {
        self.sources[s] = sources;
    }
}

// ---------------------------------------------------------------------------
// Ignition delay integration (0-D reactor)
// ---------------------------------------------------------------------------

/// Zero-dimensional adiabatic constant-pressure reactor for ignition delay.
///
/// Solves dT/dt = Q * r(T, c) / Cp and dc/dt = -r(T, c)
/// using explicit Euler until T exceeds T_ign.
pub struct ZeroDReactor {
    /// Current temperature (K).
    pub temperature: f64,
    /// Fuel concentration (mol/m^3).
    pub c_fuel: f64,
    /// Oxidizer concentration (mol/m^3).
    pub c_ox: f64,
    /// Heat of reaction per mole of fuel (J/mol).
    pub heat_of_reaction: f64,
    /// Heat capacity Cp (J/(mol·K)).
    pub cp: f64,
    /// Fuel Arrhenius rate model.
    pub arrhenius: ArrheniusRate,
    /// Elapsed physical time (s).
    pub time: f64,
}

impl ZeroDReactor {
    /// Create a new 0-D reactor.
    pub fn new(
        temperature: f64,
        c_fuel: f64,
        c_ox: f64,
        heat_of_reaction: f64,
        cp: f64,
        arrhenius: ArrheniusRate,
    ) -> Self {
        Self {
            temperature,
            c_fuel,
            c_ox,
            heat_of_reaction,
            cp,
            arrhenius,
            time: 0.0,
        }
    }

    /// Advance one explicit Euler step.
    pub fn step(&mut self, dt: f64) {
        let k = self.arrhenius.rate(self.temperature);
        let rate = k * self.c_fuel * self.c_ox;
        // Limit rate so c_fuel and c_ox stay non-negative
        let max_rate = (self.c_fuel.min(self.c_ox) / dt).max(0.0);
        let r_eff = rate.min(max_rate);
        let delta_c = r_eff * dt;
        self.c_fuel -= delta_c;
        self.c_ox -= delta_c;
        self.temperature += self.heat_of_reaction * r_eff * dt / self.cp;
        self.time += dt;
    }

    /// Integrate until temperature exceeds `t_ign` or `max_time` is reached.
    ///
    /// Returns the ignition delay time (s), or `None` if no ignition within max_time.
    pub fn ignition_delay(&mut self, t_ign: f64, dt: f64, max_time: f64) -> Option<f64> {
        let t0 = self.temperature;
        if t0 >= t_ign {
            return Some(0.0);
        }
        while self.time < max_time {
            self.step(dt);
            if self.temperature >= t_ign {
                return Some(self.time);
            }
        }
        None
    }
}
