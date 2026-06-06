// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Photochemistry molecular dynamics: excited states, nonadiabatic dynamics.
//!
//! This module provides:
//! - [`ExcitedState`]: electronic state representation with excitation energy and oscillator strength
//! - [`FranckCondonFactor`]: overlap integrals for vibronic transitions
//! - [`SurfaceHoppingMd`]: Tully's fewest-switches surface hopping algorithm
//! - [`CasscfDynamics`]: interface to CASSCF (state-averaged) gradient computations
//! - [`PhotodissociationModel`]: bond cleavage upon photoexcitation
//! - [`FluorescenceLifetime`]: radiative and non-radiative decay rates and quantum yield
//! - [`EnergyTransfer`]: Förster (FRET) and Dexter energy transfer
//! - [`ConicalIntersection`]: seam geometry and branching plane vectors
//! - [`AbsorptionSpectrum`]: absorption computed from oscillator strengths with Gaussian broadening
//! - [`EmissionSpectrum`]: fluorescence and phosphorescence spectrum

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// ExcitedState
// ---------------------------------------------------------------------------

/// Represents a single adiabatic electronic excited state.
#[derive(Debug, Clone)]
pub struct ExcitedState {
    /// Zero-based index of the electronic state (0 = ground state).
    pub state_index: usize,
    /// Vertical excitation energy in eV.
    pub excitation_energy_ev: f64,
    /// Oscillator strength (dimensionless), determines absorption intensity.
    pub oscillator_strength: f64,
    /// Adiabatic excitation energy (relaxed excited-state minimum) in eV.
    pub adiabatic_energy_ev: f64,
    /// Spin multiplicity (1=singlet, 3=triplet, …).
    pub multiplicity: u32,
}

impl ExcitedState {
    /// Create a new `ExcitedState`.
    pub fn new(
        state_index: usize,
        excitation_energy_ev: f64,
        oscillator_strength: f64,
        adiabatic_energy_ev: f64,
        multiplicity: u32,
    ) -> Self {
        Self {
            state_index,
            excitation_energy_ev,
            oscillator_strength,
            adiabatic_energy_ev,
            multiplicity,
        }
    }

    /// Stokes shift: difference between vertical and adiabatic excitation energies (eV).
    pub fn stokes_shift(&self) -> f64 {
        self.excitation_energy_ev - self.adiabatic_energy_ev
    }

    /// Returns `true` if this state is a singlet (multiplicity == 1).
    pub fn is_singlet(&self) -> bool {
        self.multiplicity == 1
    }
}

// ---------------------------------------------------------------------------
// FranckCondonFactor
// ---------------------------------------------------------------------------

/// Computes Franck-Condon (FC) overlap integrals for vibronic transitions.
///
/// Uses the displaced harmonic oscillator (Huang-Rhys) model.
#[derive(Debug, Clone)]
pub struct FranckCondonFactor {
    /// Huang-Rhys factor S (dimensionless displacement parameter).
    pub huang_rhys: f64,
    /// Vibrational frequency ω (in eV) — same in ground and excited state.
    pub frequency_ev: f64,
}

impl FranckCondonFactor {
    /// Create a new `FranckCondonFactor` from Huang-Rhys factor and vibrational frequency.
    pub fn new(huang_rhys: f64, frequency_ev: f64) -> Self {
        Self {
            huang_rhys,
            frequency_ev,
        }
    }

    /// Compute |⟨v'|0⟩|² — the FC factor from ground vibrational state to excited level `v_prime`.
    ///
    /// Formula: FC(0→v') = e^{-S} * S^{v'} / v'!
    pub fn fc_factor_from_ground(&self, v_prime: u32) -> f64 {
        let s = self.huang_rhys;
        let vf = v_prime as f64;
        let factorial = (1..=v_prime).map(|k| k as f64).product::<f64>().max(1.0);
        (-s).exp() * s.powf(vf) / factorial
    }

    /// Sum of FC factors over all `n_states` vibrational levels (should approach 1).
    pub fn sum_fc_factors(&self, n_states: u32) -> f64 {
        (0..n_states).map(|v| self.fc_factor_from_ground(v)).sum()
    }

    /// Reorganization energy λ = S * ħω in eV.
    pub fn reorganization_energy(&self) -> f64 {
        self.huang_rhys * self.frequency_ev
    }

    /// 0-0 transition energy given vertical excitation energy in eV.
    pub fn zero_zero_energy(&self, vertical_ev: f64) -> f64 {
        vertical_ev - self.reorganization_energy()
    }
}

// ---------------------------------------------------------------------------
// SurfaceHoppingMd
// ---------------------------------------------------------------------------

/// Tully's fewest-switches surface hopping (FSSH) state.
///
/// Propagates electronic populations and determines hopping probabilities
/// according to the Tully FSSH algorithm.
#[derive(Debug, Clone)]
pub struct SurfaceHoppingMd {
    /// Number of electronic states.
    pub n_states: usize,
    /// Current active (adiabatic) state index.
    pub current_state: usize,
    /// Electronic density matrix elements ρ_{kl} = c_k * c_l* (real part stored as 2D flat vec).
    pub rho_real: Vec<f64>,
    /// Imaginary part of density matrix.
    pub rho_imag: Vec<f64>,
    /// Non-adiabatic coupling vectors d_{kl} · v̇ (scalar projection, n_states × n_states).
    pub nac_velocity: Vec<f64>,
    /// Hamiltonian matrix H_{kl} in eV (diagonal = adiabatic energies, off-diagonal = couplings).
    pub hamiltonian: Vec<f64>,
    /// Time step in femtoseconds.
    pub dt_fs: f64,
}

impl SurfaceHoppingMd {
    /// Create a new `SurfaceHoppingMd` starting on `initial_state`.
    pub fn new(n_states: usize, initial_state: usize, dt_fs: f64) -> Self {
        let n2 = n_states * n_states;
        let mut rho_real = vec![0.0_f64; n2];
        let rho_imag = vec![0.0_f64; n2];
        // Pure state: ρ_{kk} = 1 for current state
        rho_real[initial_state * n_states + initial_state] = 1.0;
        Self {
            n_states,
            current_state: initial_state,
            rho_real,
            rho_imag,
            nac_velocity: vec![0.0_f64; n2],
            hamiltonian: vec![0.0_f64; n2],
            dt_fs,
        }
    }

    /// Set the non-adiabatic coupling velocity projection d_{kl} · v̇.
    pub fn set_nac_velocity(&mut self, k: usize, l: usize, value: f64) {
        self.nac_velocity[k * self.n_states + l] = value;
        self.nac_velocity[l * self.n_states + k] = -value;
    }

    /// Set Hamiltonian element H_{kl}.
    pub fn set_hamiltonian(&mut self, k: usize, l: usize, value: f64) {
        self.hamiltonian[k * self.n_states + l] = value;
    }

    /// Compute hopping probability g_{k→l} from current state `k` to state `l`.
    ///
    /// g_{kl} = max(0, -2 Re(ρ_{kl} d_{kl}·v̇) dt / ρ_{kk})
    pub fn hopping_probability(&self, l: usize) -> f64 {
        let k = self.current_state;
        if k == l {
            return 0.0;
        }
        let n = self.n_states;
        let rho_kk = self.rho_real[k * n + k];
        if rho_kk.abs() < 1e-15 {
            return 0.0;
        }
        let rho_kl_real = self.rho_real[k * n + l];
        let d_kl = self.nac_velocity[k * n + l];
        let numerator = -2.0 * rho_kl_real * d_kl * self.dt_fs;
        (numerator / rho_kk).clamp(0.0, 1.0)
    }

    /// Total hopping probability away from current state (sum over all other states).
    pub fn total_hopping_probability(&self) -> f64 {
        let n = self.n_states;
        (0..n)
            .filter(|&l| l != self.current_state)
            .map(|l| self.hopping_probability(l))
            .sum::<f64>()
            .min(1.0)
    }

    /// Population of state `k` (diagonal density matrix element).
    pub fn population(&self, k: usize) -> f64 {
        self.rho_real[k * self.n_states + k]
    }

    /// Total electronic population (should equal 1 for a normalized wavefunction).
    pub fn total_population(&self) -> f64 {
        (0..self.n_states).map(|k| self.population(k)).sum()
    }

    /// Velocity rescaling factor after a hop from state `from` to state `to`.
    ///
    /// Conserves total energy: kinetic energy is adjusted along the NAC direction.
    /// Returns `None` if the hop is classically forbidden (insufficient kinetic energy).
    pub fn velocity_rescale_factor(
        &self,
        kinetic_energy: f64,
        energy_from: f64,
        energy_to: f64,
    ) -> Option<f64> {
        let delta_e = energy_to - energy_from;
        let new_ke = kinetic_energy - delta_e;
        if new_ke < 0.0 {
            None // frustrated hop
        } else {
            Some((new_ke / kinetic_energy).sqrt())
        }
    }

    /// Perform a stochastic surface hop if random number `r` < hopping probability.
    /// Returns the new state index.
    pub fn attempt_hop(&mut self, r: f64, kinetic_energy: f64, energies: &[f64]) -> usize {
        let n = self.n_states;
        let k = self.current_state;
        let mut cumulative = 0.0;
        for l in 0..n {
            if l == k {
                continue;
            }
            let prob = self.hopping_probability(l);
            cumulative += prob;
            if r < cumulative {
                let e_k = energies[k];
                let e_l = energies[l];
                if self
                    .velocity_rescale_factor(kinetic_energy, e_k, e_l)
                    .is_some()
                {
                    self.current_state = l;
                }
                break;
            }
        }
        self.current_state
    }
}

// ---------------------------------------------------------------------------
// CasscfDynamics
// ---------------------------------------------------------------------------

/// Interface to CASSCF (Complete Active Space SCF, state-averaged) gradient data.
///
/// Stores energies and gradients for each electronic state returned by
/// an external CASSCF electronic-structure code.
#[derive(Debug, Clone)]
pub struct CasscfDynamics {
    /// Number of active orbitals.
    pub n_active: usize,
    /// Number of active electrons.
    pub n_electrons: usize,
    /// State-averaging weights for each state.
    pub state_weights: Vec<f64>,
    /// Energies (eV) for each state from last CASSCF call.
    pub state_energies: Vec<f64>,
    /// Gradients (eV/Å) for each state: flat \[state\]\[atom\]\[xyz\].
    pub state_gradients: Vec<Vec<[f64; 3]>>,
}

impl CasscfDynamics {
    /// Create a new `CasscfDynamics` with uniform state-averaging weights.
    pub fn new(n_active: usize, n_electrons: usize, n_states: usize) -> Self {
        let weight = 1.0 / n_states as f64;
        Self {
            n_active,
            n_electrons,
            state_weights: vec![weight; n_states],
            state_energies: vec![0.0; n_states],
            state_gradients: vec![vec![]; n_states],
        }
    }

    /// Update state energies after an electronic structure call.
    pub fn update_energies(&mut self, energies: &[f64]) {
        self.state_energies = energies.to_vec();
    }

    /// Update gradient for a given state.
    pub fn update_gradient(&mut self, state: usize, grad: Vec<[f64; 3]>) {
        if state < self.state_gradients.len() {
            self.state_gradients[state] = grad;
        }
    }

    /// Energy gap between two states in eV.
    pub fn energy_gap(&self, state_a: usize, state_b: usize) -> f64 {
        (self.state_energies[state_a] - self.state_energies[state_b]).abs()
    }

    /// State-averaged energy (weighted mean over all states).
    pub fn state_averaged_energy(&self) -> f64 {
        self.state_energies
            .iter()
            .zip(&self.state_weights)
            .map(|(e, w)| e * w)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// PhotodissociationModel
// ---------------------------------------------------------------------------

/// Models bond cleavage following photoexcitation.
#[derive(Debug, Clone)]
pub struct PhotodissociationModel {
    /// Bond dissociation energy D_e in eV.
    pub dissociation_energy_ev: f64,
    /// Excitation energy deposited in the system (eV).
    pub excitation_energy_ev: f64,
    /// Branching ratio for each dissociation channel (must sum to ≤ 1).
    pub channel_branching: Vec<f64>,
    /// Asymptotic fragment energies for each channel (eV).
    pub channel_asymptotes: Vec<f64>,
}

impl PhotodissociationModel {
    /// Create a new `PhotodissociationModel`.
    pub fn new(
        dissociation_energy_ev: f64,
        excitation_energy_ev: f64,
        channel_branching: Vec<f64>,
        channel_asymptotes: Vec<f64>,
    ) -> Self {
        Self {
            dissociation_energy_ev,
            excitation_energy_ev,
            channel_branching,
            channel_asymptotes,
        }
    }

    /// Kinetic energy release (KER) for a given channel: E_exc - D_e - E_asymptote.
    pub fn kinetic_energy_release(&self, channel: usize) -> f64 {
        let asymptote = self.channel_asymptotes.get(channel).copied().unwrap_or(0.0);
        (self.excitation_energy_ev - self.dissociation_energy_ev - asymptote).max(0.0)
    }

    /// Total branching ratio (should be ≤ 1; remainder is non-dissociative yield).
    pub fn total_branching(&self) -> f64 {
        self.channel_branching.iter().sum::<f64>().min(1.0)
    }

    /// Whether dissociation is energetically accessible (excitation > dissociation energy).
    pub fn is_dissociative(&self) -> bool {
        self.excitation_energy_ev > self.dissociation_energy_ev
    }

    /// Mean KER averaged over channels weighted by branching ratios.
    pub fn mean_kinetic_energy_release(&self) -> f64 {
        let total_br: f64 = self.channel_branching.iter().sum::<f64>();
        if total_br < 1e-15 {
            return 0.0;
        }
        self.channel_branching
            .iter()
            .enumerate()
            .map(|(i, &br)| br * self.kinetic_energy_release(i))
            .sum::<f64>()
            / total_br
    }
}

// ---------------------------------------------------------------------------
// FluorescenceLifetime
// ---------------------------------------------------------------------------

/// Radiative and non-radiative decay rates and quantum yield.
#[derive(Debug, Clone)]
pub struct FluorescenceLifetime {
    /// Radiative rate k_r (in ns^{-1}).
    pub k_radiative: f64,
    /// Non-radiative rate k_nr (in ns^{-1}).
    pub k_nonradiative: f64,
}

impl FluorescenceLifetime {
    /// Create a new `FluorescenceLifetime`.
    pub fn new(k_radiative: f64, k_nonradiative: f64) -> Self {
        Self {
            k_radiative,
            k_nonradiative,
        }
    }

    /// Total decay rate k_total = k_r + k_nr (ns^{-1}).
    pub fn total_rate(&self) -> f64 {
        self.k_radiative + self.k_nonradiative
    }

    /// Fluorescence lifetime τ = 1 / (k_r + k_nr) in ns.
    pub fn lifetime_ns(&self) -> f64 {
        1.0 / self.total_rate()
    }

    /// Quantum yield Φ = k_r / (k_r + k_nr).
    pub fn quantum_yield(&self) -> f64 {
        self.k_radiative / self.total_rate()
    }

    /// Radiative lifetime τ_rad = 1 / k_r in ns (natural lifetime).
    pub fn radiative_lifetime_ns(&self) -> f64 {
        1.0 / self.k_radiative
    }

    /// Effect of adding a quencher with quenching rate k_q (ns^{-1}) on lifetime.
    pub fn quenched_lifetime_ns(&self, k_q: f64) -> f64 {
        1.0 / (self.total_rate() + k_q)
    }

    /// Effect of quencher on quantum yield.
    pub fn quenched_quantum_yield(&self, k_q: f64) -> f64 {
        self.k_radiative / (self.total_rate() + k_q)
    }
}

// ---------------------------------------------------------------------------
// EnergyTransfer
// ---------------------------------------------------------------------------

/// Förster (FRET) and Dexter energy transfer parameters and rates.
#[derive(Debug, Clone)]
pub struct EnergyTransfer {
    /// Donor fluorescence lifetime in ns (in absence of acceptor).
    pub donor_lifetime_ns: f64,
    /// Förster radius R_0 in nm.
    pub forster_radius_nm: f64,
    /// Donor–acceptor distance R in nm.
    pub distance_nm: f64,
    /// Electronic coupling J for Dexter transfer in eV.
    pub dexter_coupling_ev: f64,
    /// Debye screening length L for Dexter transfer in nm.
    pub dexter_screening_nm: f64,
}

impl EnergyTransfer {
    /// Create a new `EnergyTransfer`.
    pub fn new(
        donor_lifetime_ns: f64,
        forster_radius_nm: f64,
        distance_nm: f64,
        dexter_coupling_ev: f64,
        dexter_screening_nm: f64,
    ) -> Self {
        Self {
            donor_lifetime_ns,
            forster_radius_nm,
            distance_nm,
            dexter_coupling_ev,
            dexter_screening_nm,
        }
    }

    /// Förster energy transfer rate k_FRET = (1/τ_D) * (R_0/R)^6 in ns^{-1}.
    pub fn fret_rate(&self) -> f64 {
        let ratio = self.forster_radius_nm / self.distance_nm;
        (1.0 / self.donor_lifetime_ns) * ratio.powi(6)
    }

    /// FRET efficiency E = k_FRET / (k_FRET + 1/τ_D).
    pub fn fret_efficiency(&self) -> f64 {
        let k_fret = self.fret_rate();
        let k_donor = 1.0 / self.donor_lifetime_ns;
        k_fret / (k_fret + k_donor)
    }

    /// Förster radius R_0 from spectral overlap integral J, κ^2, quantum yield Φ_D, η.
    ///
    /// R_0^6 = (9000 * ln(10) * κ^2 * Φ_D * J) / (128 * π^5 * N_A * η^4)
    /// Returns R_0 in nm for typical molecular parameters.
    pub fn compute_forster_radius(
        quantum_yield_donor: f64,
        spectral_overlap_nm4_per_m: f64,
        kappa_sq: f64,
        refractive_index: f64,
    ) -> f64 {
        // Simplified: R_0 (in nm) ≈ 0.211 * (κ^2 * Φ_D * η^{-4} * J)^{1/6}  [Lakowicz form]
        // with J in units of M^{-1} cm^{-1} nm^4
        let r0_6 = 0.211_f64.powi(6) * kappa_sq * quantum_yield_donor * spectral_overlap_nm4_per_m
            / refractive_index.powi(4);
        r0_6.powf(1.0 / 6.0)
    }

    /// Dexter energy transfer rate: k_Dexter = K * J * exp(-2R/L), returns in arbitrary units.
    pub fn dexter_rate(&self, spectral_overlap_normalized: f64) -> f64 {
        let k_prefactor = self.dexter_coupling_ev.powi(2);
        k_prefactor
            * spectral_overlap_normalized
            * (-2.0 * self.distance_nm / self.dexter_screening_nm).exp()
    }

    /// Distance at which FRET rate equals donor decay rate (= R_0).
    pub fn forster_distance(&self) -> f64 {
        self.forster_radius_nm
    }
}

// ---------------------------------------------------------------------------
// ConicalIntersection
// ---------------------------------------------------------------------------

/// Conical intersection (CI) seam geometry and branching plane vectors.
///
/// The branching plane is spanned by the gradient difference vector (g)
/// and the derivative coupling vector (h).
#[derive(Debug, Clone)]
pub struct ConicalIntersection {
    /// Geometry at the CI point (flat \[x,y,z\] for each atom).
    pub geometry: Vec<[f64; 3]>,
    /// Gradient difference vector g_i = (∂E_1/∂Q_i - ∂E_2/∂Q_i)/2.
    pub gradient_diff: Vec<[f64; 3]>,
    /// Derivative coupling vector h_i = ⟨ψ_1|∂H/∂Q_i|ψ_2⟩.
    pub derivative_coupling: Vec<[f64; 3]>,
    /// Energy at the CI point (eV).
    pub energy_ev: f64,
}

impl ConicalIntersection {
    /// Create a new `ConicalIntersection`.
    pub fn new(
        geometry: Vec<[f64; 3]>,
        gradient_diff: Vec<[f64; 3]>,
        derivative_coupling: Vec<[f64; 3]>,
        energy_ev: f64,
    ) -> Self {
        Self {
            geometry,
            gradient_diff,
            derivative_coupling,
            energy_ev,
        }
    }

    /// Magnitude of the gradient difference vector |g|.
    pub fn g_norm(&self) -> f64 {
        norm_vec3_list(&self.gradient_diff)
    }

    /// Magnitude of the derivative coupling vector |h|.
    pub fn h_norm(&self) -> f64 {
        norm_vec3_list(&self.derivative_coupling)
    }

    /// Pitch of the double cone: Δ = sqrt(g^2 + h^2).
    pub fn cone_pitch(&self) -> f64 {
        (self.g_norm().powi(2) + self.h_norm().powi(2)).sqrt()
    }

    /// Tilt of the cone (asymmetry) = |g| / (|g| + |h|), in \[0, 1\].
    pub fn tilt(&self) -> f64 {
        let g = self.g_norm();
        let h = self.h_norm();
        let denom = g + h;
        if denom < 1e-15 { 0.5 } else { g / denom }
    }
}

/// Compute L2 norm of a list of 3-vectors.
fn norm_vec3_list(vecs: &[[f64; 3]]) -> f64 {
    vecs.iter()
        .map(|v| v[0].powi(2) + v[1].powi(2) + v[2].powi(2))
        .sum::<f64>()
        .sqrt()
}

// ---------------------------------------------------------------------------
// AbsorptionSpectrum
// ---------------------------------------------------------------------------

/// Computes an absorption spectrum from excited state data with Gaussian broadening.
#[derive(Debug, Clone)]
pub struct AbsorptionSpectrum {
    /// Excited states contributing to the absorption.
    pub states: Vec<ExcitedState>,
    /// Gaussian broadening half-width at half-maximum (HWHM) in eV.
    pub broadening_ev: f64,
}

impl AbsorptionSpectrum {
    /// Create a new `AbsorptionSpectrum`.
    pub fn new(states: Vec<ExcitedState>, broadening_ev: f64) -> Self {
        Self {
            states,
            broadening_ev,
        }
    }

    /// Compute absorption intensity at energy `e` (eV).
    ///
    /// σ(E) = Σ_i f_i * Gaussian(E - E_i, σ)
    pub fn intensity_at(&self, energy_ev: f64) -> f64 {
        let sigma = self.broadening_ev / (2.0 * (2.0_f64.ln()).sqrt());
        self.states
            .iter()
            .map(|s| {
                let de = energy_ev - s.excitation_energy_ev;
                s.oscillator_strength * (-0.5 * (de / sigma).powi(2)).exp()
                    / (sigma * (2.0 * PI).sqrt())
            })
            .sum()
    }

    /// Return the energy at which the spectrum is maximum (sampled over a grid).
    ///
    /// Scans `n_points` evenly spaced points in \[e_min, e_max\].
    pub fn absorption_maximum(&self, e_min: f64, e_max: f64, n_points: usize) -> f64 {
        let step = (e_max - e_min) / (n_points - 1) as f64;
        let mut best_e = e_min;
        let mut best_i = f64::NEG_INFINITY;
        for i in 0..n_points {
            let e = e_min + i as f64 * step;
            let intensity = self.intensity_at(e);
            if intensity > best_i {
                best_i = intensity;
                best_e = e;
            }
        }
        best_e
    }

    /// Compute spectrum on a grid; returns (energies, intensities).
    pub fn spectrum_grid(&self, e_min: f64, e_max: f64, n_points: usize) -> (Vec<f64>, Vec<f64>) {
        let step = (e_max - e_min) / (n_points - 1) as f64;
        let energies: Vec<f64> = (0..n_points).map(|i| e_min + i as f64 * step).collect();
        let intensities: Vec<f64> = energies.iter().map(|&e| self.intensity_at(e)).collect();
        (energies, intensities)
    }
}

// ---------------------------------------------------------------------------
// EmissionSpectrum
// ---------------------------------------------------------------------------

/// Fluorescence and phosphorescence emission spectrum.
#[derive(Debug, Clone)]
pub struct EmissionSpectrum {
    /// Excited states that emit (singlets → fluorescence, triplets → phosphorescence).
    pub states: Vec<ExcitedState>,
    /// Gaussian broadening HWHM in eV.
    pub broadening_ev: f64,
    /// Franck-Condon parameters for each state.
    pub fc_factors: Vec<FranckCondonFactor>,
    /// Number of vibrational levels to include in FC progression.
    pub n_vib_levels: u32,
}

impl EmissionSpectrum {
    /// Create a new `EmissionSpectrum`.
    pub fn new(
        states: Vec<ExcitedState>,
        broadening_ev: f64,
        fc_factors: Vec<FranckCondonFactor>,
        n_vib_levels: u32,
    ) -> Self {
        Self {
            states,
            broadening_ev,
            fc_factors,
            n_vib_levels,
        }
    }

    /// Compute emission intensity at `energy_ev` using FC progressions.
    ///
    /// Each vibronic line is at E_00 + v * ω (for v = 0, 1, …) and weighted by FC factor.
    pub fn intensity_at(&self, energy_ev: f64) -> f64 {
        let sigma = self.broadening_ev / (2.0 * (2.0_f64.ln()).sqrt());
        let mut total = 0.0;
        for (idx, state) in self.states.iter().enumerate() {
            if let Some(fc) = self.fc_factors.get(idx) {
                let e_00 = state.adiabatic_energy_ev;
                for v in 0..self.n_vib_levels {
                    let e_line = e_00 - v as f64 * fc.frequency_ev; // emission red-shifts
                    let fc_val = fc.fc_factor_from_ground(v);
                    let de = energy_ev - e_line;
                    total +=
                        state.oscillator_strength * fc_val * (-0.5 * (de / sigma).powi(2)).exp()
                            / (sigma * (2.0 * PI).sqrt());
                }
            }
        }
        total
    }

    /// Energy of maximum emission (0-0 transition of strongest state).
    pub fn emission_maximum_ev(&self) -> f64 {
        self.states
            .iter()
            .max_by(|a, b| {
                a.oscillator_strength
                    .partial_cmp(&b.oscillator_strength)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|s| s.adiabatic_energy_ev)
            .unwrap_or(0.0)
    }

    /// Determine whether a state emits fluorescence (singlet) or phosphorescence (triplet).
    pub fn emission_type(state: &ExcitedState) -> &'static str {
        if state.multiplicity == 1 {
            "fluorescence"
        } else {
            "phosphorescence"
        }
    }
}

// ---------------------------------------------------------------------------
// Helper utilities (private)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- ExcitedState tests ----

    #[test]
    fn test_excited_state_stokes_shift() {
        let s = ExcitedState::new(1, 3.5, 0.8, 3.2, 1);
        let shift = s.stokes_shift();
        assert!((shift - 0.3).abs() < 1e-10, "Stokes shift = {:.6}", shift);
    }

    #[test]
    fn test_excited_state_is_singlet() {
        let s1 = ExcitedState::new(1, 3.5, 0.8, 3.2, 1);
        let t1 = ExcitedState::new(1, 2.5, 0.0, 2.4, 3);
        assert!(s1.is_singlet());
        assert!(!t1.is_singlet());
    }

    #[test]
    fn test_excited_state_ground_zero_osc() {
        let s0 = ExcitedState::new(0, 0.0, 0.0, 0.0, 1);
        assert_eq!(s0.oscillator_strength, 0.0);
    }

    // ---- FranckCondonFactor tests ----

    #[test]
    fn test_fc_v0_is_exp_neg_s() {
        let fc = FranckCondonFactor::new(0.5, 0.2);
        let f0 = fc.fc_factor_from_ground(0);
        assert!((f0 - (-0.5_f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn test_fc_sum_to_one_zero_displacement() {
        // S=0: only v'=0 contributes, FC(0) = 1.0
        let fc = FranckCondonFactor::new(0.0, 0.2);
        let sum = fc.sum_fc_factors(20);
        assert!((sum - 1.0).abs() < 1e-10, "sum = {:.6}", sum);
    }

    #[test]
    fn test_fc_sum_approaches_one() {
        // For S=1.5, sum over enough levels should approach 1
        let fc = FranckCondonFactor::new(1.5, 0.18);
        let sum = fc.sum_fc_factors(30);
        assert!(sum > 0.999, "sum = {:.6}", sum);
        assert!(sum <= 1.0001, "sum = {:.6}", sum);
    }

    #[test]
    fn test_fc_reorganization_energy() {
        let fc = FranckCondonFactor::new(2.0, 0.15);
        assert!((fc.reorganization_energy() - 0.30).abs() < 1e-12);
    }

    #[test]
    fn test_fc_zero_zero_energy() {
        let fc = FranckCondonFactor::new(1.0, 0.2);
        let zze = fc.zero_zero_energy(4.0);
        // 4.0 - 1.0*0.2 = 3.8
        assert!((zze - 3.8).abs() < 1e-12);
    }

    #[test]
    fn test_fc_factor_nonnegative() {
        let fc = FranckCondonFactor::new(2.0, 0.1);
        for v in 0..20u32 {
            assert!(fc.fc_factor_from_ground(v) >= 0.0);
        }
    }

    // ---- SurfaceHoppingMd tests ----

    #[test]
    fn test_shmd_initial_population() {
        let shmd = SurfaceHoppingMd::new(3, 1, 0.5);
        assert!((shmd.population(1) - 1.0).abs() < 1e-12);
        assert!(shmd.population(0).abs() < 1e-12);
        assert!(shmd.population(2).abs() < 1e-12);
    }

    #[test]
    fn test_shmd_total_population_one() {
        let shmd = SurfaceHoppingMd::new(3, 0, 0.5);
        assert!((shmd.total_population() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shmd_hopping_probability_zero_coupling() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5);
        // No coupling set → hopping probability = 0
        let p = shmd.hopping_probability(1);
        assert!(p.abs() < 1e-12, "p = {:.6}", p);
    }

    #[test]
    fn test_shmd_hopping_probability_bounds() {
        let mut shmd = SurfaceHoppingMd::new(2, 0, 0.5);
        shmd.set_nac_velocity(0, 1, 0.1);
        // Set some off-diagonal density matrix element
        shmd.rho_real[1] = 0.3;
        let p = shmd.hopping_probability(1);
        assert!((0.0..=1.0).contains(&p), "p = {:.6}", p);
    }

    #[test]
    fn test_shmd_hop_to_self_is_zero() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5);
        assert_eq!(shmd.hopping_probability(0), 0.0);
    }

    #[test]
    fn test_shmd_velocity_rescale_downward_hop() {
        let shmd = SurfaceHoppingMd::new(2, 1, 0.5);
        // Hop downward: energy_to < energy_from → new_ke > kinetic_energy → Some(>1)
        let factor = shmd.velocity_rescale_factor(5.0, 3.0, 1.0);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f > 1.0, "f = {:.6}", f);
    }

    #[test]
    fn test_shmd_velocity_rescale_frustrated_hop() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5);
        // Hop upward by more than kinetic energy → frustrated hop
        let factor = shmd.velocity_rescale_factor(1.0, 2.0, 10.0);
        assert!(factor.is_none());
    }

    #[test]
    fn test_shmd_velocity_rescale_energy_conservation() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5);
        let ke = 3.0;
        let de = 1.0; // upward hop
        let factor = shmd.velocity_rescale_factor(ke, 0.0, de);
        assert!(factor.is_some());
        let f = factor.unwrap();
        // new_ke = f^2 * ke = ke - de
        let new_ke = f * f * ke;
        assert!((new_ke - (ke - de)).abs() < 1e-10, "new_ke = {:.6}", new_ke);
    }

    #[test]
    fn test_shmd_total_hopping_bounded() {
        let mut shmd = SurfaceHoppingMd::new(3, 0, 0.5);
        shmd.set_nac_velocity(0, 1, 0.2);
        shmd.set_nac_velocity(0, 2, 0.1);
        shmd.rho_real[1] = 0.3;
        shmd.rho_real[2] = 0.2;
        let total = shmd.total_hopping_probability();
        assert!((0.0..=1.0).contains(&total), "total = {:.6}", total);
    }

    // ---- CasscfDynamics tests ----

    #[test]
    fn test_casscf_energy_gap() {
        let mut cas = CasscfDynamics::new(4, 4, 2);
        cas.update_energies(&[0.0, 3.5]);
        assert!((cas.energy_gap(0, 1) - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_casscf_uniform_weights() {
        let cas = CasscfDynamics::new(4, 4, 3);
        let w_sum: f64 = cas.state_weights.iter().sum();
        assert!((w_sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_casscf_state_averaged_energy() {
        let mut cas = CasscfDynamics::new(2, 2, 2);
        cas.update_energies(&[0.0, 2.0]);
        // Uniform weights 0.5 each → avg = 1.0
        let avg = cas.state_averaged_energy();
        assert!((avg - 1.0).abs() < 1e-10, "avg = {:.6}", avg);
    }

    // ---- PhotodissociationModel tests ----

    #[test]
    fn test_photodiss_ker_positive() {
        let pd = PhotodissociationModel::new(4.0, 6.0, vec![0.7, 0.3], vec![0.5, 1.0]);
        let ker0 = pd.kinetic_energy_release(0);
        assert!(ker0 >= 0.0, "KER0 = {:.6}", ker0);
    }

    #[test]
    fn test_photodiss_ker_energy_conservation() {
        // KER = excitation - D_e - asymptote
        let pd = PhotodissociationModel::new(3.0, 5.5, vec![1.0], vec![0.5]);
        let ker = pd.kinetic_energy_release(0);
        // 5.5 - 3.0 - 0.5 = 2.0
        assert!((ker - 2.0).abs() < 1e-10, "KER = {:.6}", ker);
    }

    #[test]
    fn test_photodiss_total_branching_le_one() {
        let pd = PhotodissociationModel::new(4.0, 6.0, vec![0.6, 0.5], vec![0.0, 0.0]);
        assert!(pd.total_branching() <= 1.0);
    }

    #[test]
    fn test_photodiss_dissociative_flag() {
        let pd1 = PhotodissociationModel::new(4.0, 5.0, vec![], vec![]);
        let pd2 = PhotodissociationModel::new(4.0, 3.0, vec![], vec![]);
        assert!(pd1.is_dissociative());
        assert!(!pd2.is_dissociative());
    }

    #[test]
    fn test_photodiss_mean_ker() {
        let pd = PhotodissociationModel::new(2.0, 5.0, vec![0.5, 0.5], vec![0.0, 1.0]);
        let mean_ker = pd.mean_kinetic_energy_release();
        // ch0: 5-2-0=3.0, ch1: 5-2-1=2.0; weighted avg = 2.5
        assert!((mean_ker - 2.5).abs() < 1e-10, "mean_ker = {:.6}", mean_ker);
    }

    // ---- FluorescenceLifetime tests ----

    #[test]
    fn test_fluorescence_quantum_yield_le_one() {
        let fl = FluorescenceLifetime::new(0.1, 0.5);
        let qy = fl.quantum_yield();
        assert!((0.0..=1.0).contains(&qy), "qy = {:.6}", qy);
    }

    #[test]
    fn test_fluorescence_quantum_yield_value() {
        let fl = FluorescenceLifetime::new(0.2, 0.8);
        let qy = fl.quantum_yield();
        // 0.2 / (0.2 + 0.8) = 0.2
        assert!((qy - 0.2).abs() < 1e-10, "qy = {:.6}", qy);
    }

    #[test]
    fn test_fluorescence_lifetime_decreases_with_knr() {
        let fl1 = FluorescenceLifetime::new(0.1, 0.1);
        let fl2 = FluorescenceLifetime::new(0.1, 0.9);
        assert!(fl2.lifetime_ns() < fl1.lifetime_ns());
    }

    #[test]
    fn test_fluorescence_lifetime_formula() {
        let fl = FluorescenceLifetime::new(1.0, 4.0);
        // τ = 1/(1+4) = 0.2
        assert!((fl.lifetime_ns() - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_fluorescence_quenched_yield_lower() {
        let fl = FluorescenceLifetime::new(0.3, 0.2);
        let qy_unquenched = fl.quantum_yield();
        let qy_quenched = fl.quenched_quantum_yield(1.0);
        assert!(qy_quenched < qy_unquenched);
    }

    #[test]
    fn test_fluorescence_radiative_lifetime() {
        let fl = FluorescenceLifetime::new(2.0, 8.0);
        assert!((fl.radiative_lifetime_ns() - 0.5).abs() < 1e-10);
    }

    // ---- EnergyTransfer tests ----

    #[test]
    fn test_fret_rate_r6_dependence() {
        // FRET rate ∝ R^{-6}: doubling R → rate drops by 64
        let et1 = EnergyTransfer::new(10.0, 5.0, 5.0, 0.01, 1.0);
        let et2 = EnergyTransfer::new(10.0, 5.0, 10.0, 0.01, 1.0);
        let ratio = et1.fret_rate() / et2.fret_rate();
        assert!((ratio - 64.0).abs() < 1e-6, "ratio = {:.6}", ratio);
    }

    #[test]
    fn test_fret_efficiency_at_r0_is_half() {
        // At R = R_0, efficiency = 0.5
        let et = EnergyTransfer::new(10.0, 5.0, 5.0, 0.01, 1.0);
        let eff = et.fret_efficiency();
        assert!((eff - 0.5).abs() < 1e-10, "eff = {:.6}", eff);
    }

    #[test]
    fn test_fret_efficiency_range() {
        let et = EnergyTransfer::new(10.0, 5.0, 3.0, 0.01, 1.0);
        let eff = et.fret_efficiency();
        assert!((0.0..=1.0).contains(&eff), "eff = {:.6}", eff);
    }

    #[test]
    fn test_dexter_rate_decreases_with_distance() {
        let et1 = EnergyTransfer::new(10.0, 5.0, 1.0, 0.1, 1.0);
        let et2 = EnergyTransfer::new(10.0, 5.0, 2.0, 0.1, 1.0);
        assert!(et2.dexter_rate(1.0) < et1.dexter_rate(1.0));
    }

    // ---- ConicalIntersection tests ----

    #[test]
    fn test_ci_g_norm_positive() {
        let g = vec![[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let h = vec![[0.0_f64, 1.0, 0.0], [1.0, 0.0, 0.0]];
        let ci = ConicalIntersection::new(vec![], g, h, 3.5);
        assert!(ci.g_norm() > 0.0);
        assert!(ci.h_norm() > 0.0);
    }

    #[test]
    fn test_ci_tilt_range() {
        let g = vec![[1.0_f64, 0.0, 0.0]];
        let h = vec![[0.0_f64, 1.0, 0.0]];
        let ci = ConicalIntersection::new(vec![], g, h, 2.5);
        let tilt = ci.tilt();
        assert!((0.0..=1.0).contains(&tilt), "tilt = {:.6}", tilt);
    }

    #[test]
    fn test_ci_pitch_positive() {
        let g = vec![[0.5_f64, 0.3, 0.0]];
        let h = vec![[0.0_f64, 0.4, 0.2]];
        let ci = ConicalIntersection::new(vec![], g, h, 3.0);
        assert!(ci.cone_pitch() > 0.0);
    }

    // ---- AbsorptionSpectrum tests ----

    #[test]
    fn test_absorption_maximum_near_excitation() {
        let states = vec![ExcitedState::new(1, 4.0, 1.0, 3.7, 1)];
        let spec = AbsorptionSpectrum::new(states, 0.1);
        let max_e = spec.absorption_maximum(3.0, 5.0, 500);
        assert!((max_e - 4.0).abs() < 0.05, "max_e = {:.6}", max_e);
    }

    #[test]
    fn test_absorption_intensity_positive() {
        let states = vec![ExcitedState::new(1, 3.5, 0.5, 3.2, 1)];
        let spec = AbsorptionSpectrum::new(states, 0.1);
        let i = spec.intensity_at(3.5);
        assert!(i > 0.0);
    }

    #[test]
    fn test_absorption_zero_far_from_peak() {
        let states = vec![ExcitedState::new(1, 3.5, 0.5, 3.2, 1)];
        let spec = AbsorptionSpectrum::new(states, 0.05);
        // Very far from peak: ~10 sigma away
        let i = spec.intensity_at(10.0);
        assert!(i < 1e-10, "i = {:.6}", i);
    }

    #[test]
    fn test_absorption_spectrum_grid_length() {
        let states = vec![ExcitedState::new(1, 3.5, 0.5, 3.2, 1)];
        let spec = AbsorptionSpectrum::new(states, 0.1);
        let (energies, intensities) = spec.spectrum_grid(2.0, 6.0, 100);
        assert_eq!(energies.len(), 100);
        assert_eq!(intensities.len(), 100);
    }

    #[test]
    fn test_absorption_proportional_to_oscillator_strength() {
        let states1 = vec![ExcitedState::new(1, 3.5, 1.0, 3.2, 1)];
        let states2 = vec![ExcitedState::new(1, 3.5, 2.0, 3.2, 1)];
        let spec1 = AbsorptionSpectrum::new(states1, 0.1);
        let spec2 = AbsorptionSpectrum::new(states2, 0.1);
        let ratio = spec2.intensity_at(3.5) / spec1.intensity_at(3.5);
        assert!((ratio - 2.0).abs() < 1e-10, "ratio = {:.6}", ratio);
    }

    // ---- EmissionSpectrum tests ----

    #[test]
    fn test_emission_type_singlet() {
        let s = ExcitedState::new(1, 3.5, 0.8, 3.2, 1);
        assert_eq!(EmissionSpectrum::emission_type(&s), "fluorescence");
    }

    #[test]
    fn test_emission_type_triplet() {
        let t = ExcitedState::new(1, 2.5, 0.0, 2.4, 3);
        assert_eq!(EmissionSpectrum::emission_type(&t), "phosphorescence");
    }

    #[test]
    fn test_emission_maximum_ev() {
        let states = vec![
            ExcitedState::new(1, 3.5, 0.1, 3.2, 1),
            ExcitedState::new(2, 4.0, 0.9, 3.8, 1),
        ];
        let fcs = vec![
            FranckCondonFactor::new(0.5, 0.15),
            FranckCondonFactor::new(0.5, 0.15),
        ];
        let em = EmissionSpectrum::new(states, 0.1, fcs, 5);
        // Strongest state (osc=0.9) has adiabatic energy 3.8
        let max_ev = em.emission_maximum_ev();
        assert!((max_ev - 3.8).abs() < 1e-10, "max_ev = {:.6}", max_ev);
    }

    #[test]
    fn test_emission_intensity_positive_at_max() {
        let states = vec![ExcitedState::new(1, 3.5, 0.8, 3.2, 1)];
        let fcs = vec![FranckCondonFactor::new(0.5, 0.15)];
        let em = EmissionSpectrum::new(states, 0.1, fcs, 5);
        let i = em.intensity_at(3.2);
        assert!(i > 0.0);
    }

    #[test]
    fn test_emission_intensity_zero_far_away() {
        let states = vec![ExcitedState::new(1, 3.5, 0.8, 3.2, 1)];
        let fcs = vec![FranckCondonFactor::new(0.5, 0.15)];
        let em = EmissionSpectrum::new(states, 0.05, fcs, 5);
        let i = em.intensity_at(10.0);
        assert!(i < 1e-10);
    }
}
