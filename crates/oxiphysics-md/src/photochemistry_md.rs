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

use crate::error::MdError;
use std::f64::consts::PI;

/// Reduced Planck constant ħ in eV·fs (CODATA 2018: 6.582 119 569 × 10⁻¹⁶ eV·s).
pub const HBAR_EV_FS: f64 = 0.658_211_956_9;

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
/// Carries the complex electronic amplitudes `c_k`, propagates them with
/// [`SurfaceHoppingMd::propagate_electronic`] and determines hopping
/// probabilities according to the Tully FSSH algorithm.
///
/// The amplitudes obey the time-dependent Schrödinger equation in the
/// (adiabatic) basis,
///
/// ```text
/// dc_k/dt = −(i/ħ) Σ_l H_kl c_l − Σ_l σ_kl c_l ,   σ_kl = d_kl · v ,
/// ```
///
/// with the real symmetric Hamiltonian `H` (eV; diagonal = adiabatic
/// energies) and the real antisymmetric scalar non-adiabatic coupling `σ`
/// (fs⁻¹).  The density matrix `ρ_kl = c_k c_l*` is rebuilt from the
/// amplitudes after every propagation.
///
/// Per nuclear step the caller supplies the current energies and couplings
/// ([`set_hamiltonian`](SurfaceHoppingMd::set_hamiltonian),
/// [`set_nac_velocity`](SurfaceHoppingMd::set_nac_velocity)), draws a random
/// number `r ∈ [0, 1)`, calls [`attempt_hop`](SurfaceHoppingMd::attempt_hop)
/// with the hop probabilities of the *current* `ρ`, rescales the nuclear
/// velocities by the returned factor, and then advances the electrons with
/// [`propagate_electronic`](SurfaceHoppingMd::propagate_electronic).
#[derive(Debug, Clone)]
pub struct SurfaceHoppingMd {
    /// Number of electronic states.
    pub n_states: usize,
    /// Current active (adiabatic) state index.
    pub current_state: usize,
    /// Real parts of the electronic amplitudes `c_k`.
    pub amp_real: Vec<f64>,
    /// Imaginary parts of the electronic amplitudes `c_k`.
    pub amp_imag: Vec<f64>,
    /// Electronic density matrix elements ρ_{kl} = c_k * c_l^* (real part, row-major n × n).
    ///
    /// Derived from the amplitudes; rebuilt by `propagate_electronic` and `set_amplitudes`.
    pub rho_real: Vec<f64>,
    /// Imaginary part of the density matrix.
    pub rho_imag: Vec<f64>,
    /// Non-adiabatic coupling vectors d_{kl} · v̇ (scalar projection, n_states × n_states, fs⁻¹).
    ///
    /// Real and antisymmetric (`σ_lk = −σ_kl`); `set_nac_velocity` keeps it so.
    pub nac_velocity: Vec<f64>,
    /// Hamiltonian matrix H_{kl} in eV (diagonal = adiabatic energies, off-diagonal = couplings).
    ///
    /// Real and symmetric; `set_hamiltonian` keeps it so.
    pub hamiltonian: Vec<f64>,
    /// Time step in femtoseconds.
    pub dt_fs: f64,
}

impl SurfaceHoppingMd {
    /// Create a new `SurfaceHoppingMd` starting on `initial_state`.
    ///
    /// The electrons start in the pure adiabatic state `initial_state`
    /// (`c_k = δ_{k,initial_state}`).
    ///
    /// # Errors
    /// Returns an error if `n_states` is zero, if `initial_state >= n_states`,
    /// or if `dt_fs` is negative or not finite.
    pub fn new(n_states: usize, initial_state: usize, dt_fs: f64) -> Result<Self, MdError> {
        if n_states == 0 {
            return Err(MdError::General(
                "surface hopping needs at least one electronic state".to_string(),
            ));
        }
        if initial_state >= n_states {
            return Err(MdError::General(format!(
                "initial state {initial_state} is out of range for {n_states} electronic states"
            )));
        }
        if !dt_fs.is_finite() || dt_fs < 0.0 {
            return Err(MdError::General(format!(
                "surface hopping time step must be finite and non-negative, got {dt_fs} fs"
            )));
        }
        let n2 = n_states * n_states;
        let mut amp_real = vec![0.0_f64; n_states];
        // Pure state: c_k = δ_{k,initial}
        amp_real[initial_state] = 1.0;
        let mut state = Self {
            n_states,
            current_state: initial_state,
            amp_real,
            amp_imag: vec![0.0_f64; n_states],
            rho_real: vec![0.0_f64; n2],
            rho_imag: vec![0.0_f64; n2],
            nac_velocity: vec![0.0_f64; n2],
            hamiltonian: vec![0.0_f64; n2],
            dt_fs,
        };
        state.rebuild_density_matrix();
        Ok(state)
    }

    /// Set the non-adiabatic coupling velocity projection d_{kl} · v̇ (fs⁻¹).
    ///
    /// Stores `σ_kl = value` and `σ_lk = −value` (the coupling is antisymmetric).
    ///
    /// # Panics
    /// Panics if `k` or `l` is not smaller than `n_states`.
    pub fn set_nac_velocity(&mut self, k: usize, l: usize, value: f64) {
        self.nac_velocity[k * self.n_states + l] = value;
        self.nac_velocity[l * self.n_states + k] = -value;
    }

    /// Set Hamiltonian element H_{kl} (eV).
    ///
    /// The electronic Hamiltonian is real symmetric, so both `H_kl` and `H_lk`
    /// are set to `value`.
    ///
    /// # Panics
    /// Panics if `k` or `l` is not smaller than `n_states`.
    pub fn set_hamiltonian(&mut self, k: usize, l: usize, value: f64) {
        self.hamiltonian[k * self.n_states + l] = value;
        self.hamiltonian[l * self.n_states + k] = value;
    }

    /// Replace the electronic amplitudes (for example to start from a coherent superposition).
    ///
    /// `real` and `imag` hold `Re c_k` and `Im c_k`; the density matrix is
    /// rebuilt.  The amplitudes are not renormalised: a physical wavefunction
    /// has `Σ|c_k|² = 1`.
    ///
    /// # Errors
    /// Returns an error if either slice does not have `n_states` entries or
    /// contains a non-finite value.
    pub fn set_amplitudes(&mut self, real: &[f64], imag: &[f64]) -> Result<(), MdError> {
        if real.len() != self.n_states || imag.len() != self.n_states {
            return Err(MdError::General(format!(
                "expected {} amplitudes, got {} real and {} imaginary parts",
                self.n_states,
                real.len(),
                imag.len()
            )));
        }
        if real.iter().chain(imag).any(|x| !x.is_finite()) {
            return Err(MdError::General(
                "electronic amplitudes must be finite".to_string(),
            ));
        }
        self.amp_real.copy_from_slice(real);
        self.amp_imag.copy_from_slice(imag);
        self.rebuild_density_matrix();
        Ok(())
    }

    /// Complex electronic amplitude `c_k` as `(Re, Im)`.
    ///
    /// # Panics
    /// Panics if `k >= n_states`.
    pub fn amplitude(&self, k: usize) -> (f64, f64) {
        (self.amp_real[k], self.amp_imag[k])
    }

    /// Rebuild `ρ_kl = c_k c_l*` from the amplitudes.
    fn rebuild_density_matrix(&mut self) {
        let n = self.n_states;
        for k in 0..n {
            for l in 0..n {
                let (ar, ai) = (self.amp_real[k], self.amp_imag[k]);
                let (br, bi) = (self.amp_real[l], self.amp_imag[l]);
                self.rho_real[k * n + l] = ar * br + ai * bi;
                self.rho_imag[k * n + l] = ai * br - ar * bi;
            }
        }
    }

    /// Time derivative `(da/dt, db/dt)` of the real/imaginary amplitude parts.
    ///
    /// With `c = a + i b` and the reference energy `shift` (eV) removed from the
    /// diagonal of `H`: `da/dt = −σ a + (H'/ħ) b` and `db/dt = −σ b − (H'/ħ) a`.
    fn amplitude_derivative(&self, a: &[f64], b: &[f64], shift: f64) -> (Vec<f64>, Vec<f64>) {
        let n = self.n_states;
        let mut da = vec![0.0_f64; n];
        let mut db = vec![0.0_f64; n];
        for k in 0..n {
            let mut sum_a = 0.0;
            let mut sum_b = 0.0;
            for l in 0..n {
                let sigma = self.nac_velocity[k * n + l];
                let h_kl = self.hamiltonian[k * n + l] - if k == l { shift } else { 0.0 };
                let h_over_hbar = h_kl / HBAR_EV_FS;
                sum_a += -sigma * a[l] + h_over_hbar * b[l];
                sum_b += -sigma * b[l] - h_over_hbar * a[l];
            }
            da[k] = sum_a;
            db[k] = sum_b;
        }
        (da, db)
    }

    /// Advance the electronic amplitudes over one nuclear step `dt_fs`.
    ///
    /// Integrates `dc_k/dt = −(i/ħ) Σ_l H_kl c_l − Σ_l σ_kl c_l` with
    /// `n_sub` classical fourth-order Runge-Kutta sub-steps of length
    /// `dt_fs / n_sub` (`n_sub = 0` is treated as 1), holding `H` and `σ`
    /// fixed at their current values, and rebuilds `ρ = c c†`.
    ///
    /// The population-weighted mean energy `Ē = Σ_k ρ_kk H_kk / Σ_k ρ_kk` is
    /// subtracted from the diagonal of `H` for the integration.  That divides
    /// the global phase `exp(−i Ē t / ħ)` out of every amplitude, which leaves
    /// every `ρ_kl` (and so every population, coherence and hopping
    /// probability) unchanged but keeps the Runge-Kutta sub-steps small
    /// compared with the electronic oscillation periods, so the individual
    /// amplitudes are defined up to a global phase.  Choose `n_sub` so that
    /// `|H_kk − Ē| dt_fs / (ħ n_sub) ≪ 1`; the integrator is not exactly
    /// unitary, so monitor [`total_population`](Self::total_population) for
    /// large energy gaps.
    pub fn propagate_electronic(&mut self, n_sub: usize) {
        let n = self.n_states;
        let n_sub = n_sub.max(1);
        let h = self.dt_fs / n_sub as f64;
        let total: f64 = (0..n).map(|k| self.population(k)).sum();
        let shift = if total > 0.0 {
            (0..n)
                .map(|k| self.population(k) * self.hamiltonian[k * n + k])
                .sum::<f64>()
                / total
        } else {
            0.0
        };
        let mut a = self.amp_real.clone();
        let mut b = self.amp_imag.clone();
        for _ in 0..n_sub {
            let (k1a, k1b) = self.amplitude_derivative(&a, &b, shift);
            let (a2, b2) = Self::shifted(&a, &b, &k1a, &k1b, 0.5 * h);
            let (k2a, k2b) = self.amplitude_derivative(&a2, &b2, shift);
            let (a3, b3) = Self::shifted(&a, &b, &k2a, &k2b, 0.5 * h);
            let (k3a, k3b) = self.amplitude_derivative(&a3, &b3, shift);
            let (a4, b4) = Self::shifted(&a, &b, &k3a, &k3b, h);
            let (k4a, k4b) = self.amplitude_derivative(&a4, &b4, shift);
            for k in 0..n {
                a[k] += h / 6.0 * (k1a[k] + 2.0 * k2a[k] + 2.0 * k3a[k] + k4a[k]);
                b[k] += h / 6.0 * (k1b[k] + 2.0 * k2b[k] + 2.0 * k3b[k] + k4b[k]);
            }
        }
        self.amp_real = a;
        self.amp_imag = b;
        self.rebuild_density_matrix();
    }

    /// `(a + step * da, b + step * db)`.
    fn shifted(a: &[f64], b: &[f64], da: &[f64], db: &[f64], step: f64) -> (Vec<f64>, Vec<f64>) {
        (
            a.iter().zip(da).map(|(x, d)| x + step * d).collect(),
            b.iter().zip(db).map(|(x, d)| x + step * d).collect(),
        )
    }

    /// Compute hopping probability g_{k→l} from current state `k` to state `l`.
    ///
    /// ```text
    /// g_{kl} = max(0, [2 Re(ρ_{kl}) σ_{kl} + (2/ħ) H_{kl} Im(ρ_{kl})] dt / ρ_{kk})
    /// ```
    ///
    /// i.e. `dt` times the rate at which the coupling to state `l` depletes the
    /// population of `k` (`−ρ̇_kk|_l`), normalised by `ρ_kk` and clamped to
    /// `[0, 1]`.  In the adiabatic basis `H_kl = 0` for `k ≠ l` and only the
    /// non-adiabatic coupling term `σ_kl = d_kl · v` remains.
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
        let rho_kl_imag = self.rho_imag[k * n + l];
        let sigma_kl = self.nac_velocity[k * n + l];
        let h_kl = self.hamiltonian[k * n + l];
        let flux = 2.0 * rho_kl_real * sigma_kl + 2.0 * h_kl * rho_kl_imag / HBAR_EV_FS;
        (flux * self.dt_fs / rho_kk).clamp(0.0, 1.0)
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
    /// Conserves total energy: multiplying every nuclear velocity by the
    /// returned factor changes the kinetic energy to `kinetic_energy −
    /// (energy_to − energy_from)`.
    /// Returns `None` if the hop is classically forbidden (insufficient kinetic
    /// energy), or if the kinetic energy is zero and would have to *increase*:
    /// a uniform rescaling cannot create motion from rest.
    pub fn velocity_rescale_factor(
        &self,
        kinetic_energy: f64,
        energy_from: f64,
        energy_to: f64,
    ) -> Option<f64> {
        let delta_e = energy_to - energy_from;
        let new_ke = kinetic_energy - delta_e;
        if new_ke.is_nan() || new_ke < 0.0 {
            return None; // frustrated hop
        }
        if kinetic_energy > 0.0 {
            Some((new_ke / kinetic_energy).sqrt())
        } else if new_ke == 0.0 {
            Some(1.0) // at rest and isoenergetic: nothing to rescale
        } else {
            None
        }
    }

    /// Perform a stochastic surface hop if random number `r` < hopping probability.
    ///
    /// The target state is the one whose cumulative-probability interval
    /// contains `r` (Tully's fewest-switches rule).  A hop is accepted only if
    /// the nuclei can supply the energy difference (`energies` are the
    /// adiabatic energies of all states, `kinetic_energy` the nuclear kinetic
    /// energy, same units).
    ///
    /// Returns `Some((new_state, velocity_scale))` after an accepted hop —
    /// the caller must multiply the nuclear velocities by `velocity_scale` to
    /// conserve the total energy — and `None` if no hop occurred (`r` above
    /// the total probability, a frustrated hop, or fewer than `n_states`
    /// energies supplied); `current_state` is then unchanged.
    pub fn attempt_hop(
        &mut self,
        r: f64,
        kinetic_energy: f64,
        energies: &[f64],
    ) -> Option<(usize, f64)> {
        let n = self.n_states;
        let k = self.current_state;
        let mut cumulative = 0.0;
        for l in 0..n {
            if l == k {
                continue;
            }
            cumulative += self.hopping_probability(l);
            if r < cumulative {
                let energy_from = *energies.get(k)?;
                let energy_to = *energies.get(l)?;
                let scale = self.velocity_rescale_factor(kinetic_energy, energy_from, energy_to)?;
                self.current_state = l;
                return Some((l, scale));
            }
        }
        None
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
        let shmd = SurfaceHoppingMd::new(3, 1, 0.5).unwrap();
        assert!((shmd.population(1) - 1.0).abs() < 1e-12);
        assert!(shmd.population(0).abs() < 1e-12);
        assert!(shmd.population(2).abs() < 1e-12);
    }

    #[test]
    fn test_shmd_total_population_one() {
        let shmd = SurfaceHoppingMd::new(3, 0, 0.5).unwrap();
        assert!((shmd.total_population() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_shmd_hopping_probability_zero_coupling() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        // No coupling set → hopping probability = 0
        let p = shmd.hopping_probability(1);
        assert!(p.abs() < 1e-12, "p = {:.6}", p);
    }

    #[test]
    fn test_shmd_hopping_probability_bounds() {
        let mut shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        shmd.set_nac_velocity(0, 1, 0.1);
        // A coherent superposition gives a non-zero off-diagonal density matrix element.
        shmd.set_amplitudes(&[0.9, 0.3], &[0.0, 0.0]).unwrap();
        let p = shmd.hopping_probability(1);
        assert!((0.0..=1.0).contains(&p), "p = {:.6}", p);
    }

    #[test]
    fn test_shmd_hop_to_self_is_zero() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        assert_eq!(shmd.hopping_probability(0), 0.0);
    }

    #[test]
    fn test_shmd_velocity_rescale_downward_hop() {
        let shmd = SurfaceHoppingMd::new(2, 1, 0.5).unwrap();
        // Hop downward: energy_to < energy_from → new_ke > kinetic_energy → Some(>1)
        let factor = shmd.velocity_rescale_factor(5.0, 3.0, 1.0);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!(f > 1.0, "f = {:.6}", f);
    }

    #[test]
    fn test_shmd_velocity_rescale_frustrated_hop() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        // Hop upward by more than kinetic energy → frustrated hop
        let factor = shmd.velocity_rescale_factor(1.0, 2.0, 10.0);
        assert!(factor.is_none());
    }

    #[test]
    fn test_shmd_velocity_rescale_energy_conservation() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
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
        let mut shmd = SurfaceHoppingMd::new(3, 0, 0.5).unwrap();
        shmd.set_nac_velocity(0, 1, 0.2);
        shmd.set_nac_velocity(0, 2, 0.1);
        shmd.set_amplitudes(&[0.8, 0.3, 0.2], &[0.0, 0.0, 0.0])
            .unwrap();
        let total = shmd.total_hopping_probability();
        assert!((0.0..=1.0).contains(&total), "total = {:.6}", total);
    }

    // ---- FSSH electronic propagation and hopping tests ----

    /// Two-state FSSH system: Hamiltonian elements in eV, NAC `sigma_01` in fs⁻¹.
    fn two_state(h00: f64, h11: f64, h01: f64, sigma_01: f64, dt_fs: f64) -> SurfaceHoppingMd {
        let mut shmd = SurfaceHoppingMd::new(2, 0, dt_fs).unwrap();
        shmd.set_hamiltonian(0, 0, h00);
        shmd.set_hamiltonian(1, 1, h11);
        shmd.set_hamiltonian(0, 1, h01);
        shmd.set_nac_velocity(0, 1, sigma_01);
        shmd
    }

    #[test]
    fn test_shmd_new_rejects_invalid_arguments() {
        assert!(SurfaceHoppingMd::new(0, 0, 0.5).is_err());
        assert!(SurfaceHoppingMd::new(2, 2, 0.5).is_err());
        assert!(SurfaceHoppingMd::new(2, 5, 0.5).is_err());
        assert!(SurfaceHoppingMd::new(2, 0, -0.5).is_err());
        assert!(SurfaceHoppingMd::new(2, 0, f64::NAN).is_err());
        assert!(SurfaceHoppingMd::new(2, 1, 0.5).is_ok());
    }

    #[test]
    fn test_shmd_set_hamiltonian_and_nac_keep_their_symmetry() {
        let mut shmd = SurfaceHoppingMd::new(3, 0, 0.5).unwrap();
        shmd.set_hamiltonian(0, 2, 0.04);
        assert_eq!(shmd.hamiltonian[2], 0.04);
        assert_eq!(shmd.hamiltonian[2 * 3], 0.04);
        shmd.set_nac_velocity(1, 2, 0.07);
        assert_eq!(shmd.nac_velocity[3 + 2], 0.07);
        assert_eq!(shmd.nac_velocity[2 * 3 + 1], -0.07);
    }

    #[test]
    fn test_shmd_set_amplitudes_validates_and_rebuilds_the_density_matrix() {
        let mut shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        assert!(shmd.set_amplitudes(&[1.0], &[0.0, 0.0]).is_err());
        assert!(shmd.set_amplitudes(&[1.0, f64::NAN], &[0.0, 0.0]).is_err());
        // c = (0.6, 0.8 i): ρ_01 = c_0 c_1* = 0.6 * (-0.8 i) = -0.48 i
        shmd.set_amplitudes(&[0.6, 0.0], &[0.0, 0.8]).unwrap();
        assert!((shmd.population(0) - 0.36).abs() < 1e-15);
        assert!((shmd.population(1) - 0.64).abs() < 1e-15);
        assert!(shmd.rho_real[1].abs() < 1e-15);
        assert!((shmd.rho_imag[1] + 0.48).abs() < 1e-15);
        assert!(
            (shmd.rho_imag[2] - 0.48).abs() < 1e-15,
            "ρ must be Hermitian"
        );
        assert_eq!(shmd.amplitude(1), (0.0, 0.8));
    }

    #[test]
    fn test_shmd_rabi_oscillation_matches_analytic_population() {
        // Two degenerate states coupled by H_01 = 0.01 eV and no NAC:
        // starting in state 0, P_1(t) = sin²(H_01 t / ħ).
        let h01 = 0.01;
        let mut shmd = two_state(0.0, 0.0, h01, 0.0, 0.5);
        let omega = h01 / HBAR_EV_FS;
        for step in 1..=200 {
            shmd.propagate_electronic(4);
            let t = step as f64 * 0.5;
            let expected = (omega * t).sin().powi(2);
            assert!(
                (shmd.population(1) - expected).abs() < 1e-8,
                "t = {t} fs: P_1 = {}, expected {expected}",
                shmd.population(1)
            );
            assert!((shmd.total_population() - 1.0).abs() < 1e-9);
        }
        assert!(shmd.population(1) > 0.99, "population has transferred");
    }

    #[test]
    fn test_shmd_detuned_rabi_oscillation_matches_analytic_population() {
        // H = [[0, V], [V, Δ]]: P_1(t) = 4V²/(4V² + Δ²) sin²(½ sqrt(4V² + Δ²) t / ħ).
        let (v, delta) = (0.05, 0.4);
        let mut shmd = two_state(0.0, delta, v, 0.0, 0.5);
        let gen_rabi = (4.0 * v * v + delta * delta).sqrt();
        for step in 1..=200 {
            shmd.propagate_electronic(16);
            let t = step as f64 * 0.5;
            let expected = 4.0 * v * v / (gen_rabi * gen_rabi)
                * (0.5 * gen_rabi * t / HBAR_EV_FS).sin().powi(2);
            assert!(
                (shmd.population(1) - expected).abs() < 1e-7,
                "t = {t} fs: P_1 = {}, expected {expected}",
                shmd.population(1)
            );
        }
    }

    #[test]
    fn test_shmd_nonadiabatic_coupling_rotates_the_amplitudes() {
        // H = 0 and σ_01 = s: c_0 = cos(s t), c_1 = sin(s t), both real.
        let s_nac = 0.03;
        let mut shmd = two_state(0.0, 0.0, 0.0, s_nac, 0.5);
        for _ in 0..100 {
            shmd.propagate_electronic(2);
        }
        let t = 50.0;
        let (c0r, c0i) = shmd.amplitude(0);
        let (c1r, c1i) = shmd.amplitude(1);
        assert!((c0r - (s_nac * t).cos()).abs() < 1e-9, "Re c0 = {c0r}");
        assert!((c1r - (s_nac * t).sin()).abs() < 1e-9, "Re c1 = {c1r}");
        assert!(c0i.abs() < 1e-12 && c1i.abs() < 1e-12);
    }

    /// Three coupled states with generic Hermitian `H` and antisymmetric `σ`.
    fn three_state(dt_fs: f64) -> SurfaceHoppingMd {
        let mut shmd = SurfaceHoppingMd::new(3, 0, dt_fs).unwrap();
        shmd.set_hamiltonian(1, 1, 0.35);
        shmd.set_hamiltonian(2, 2, 0.8);
        shmd.set_hamiltonian(0, 1, 0.05);
        shmd.set_hamiltonian(1, 2, -0.03);
        shmd.set_hamiltonian(0, 2, 0.02);
        shmd.set_nac_velocity(0, 1, 0.04);
        shmd.set_nac_velocity(1, 2, -0.02);
        shmd.set_nac_velocity(0, 2, 0.01);
        shmd.set_amplitudes(&[0.6, 0.0, 0.0], &[0.0, 0.0, 0.8])
            .unwrap();
        shmd
    }

    #[test]
    fn test_shmd_propagation_conserves_population_and_keeps_rho_pure_hermitian() {
        let mut shmd = three_state(0.5);
        for _ in 0..400 {
            shmd.propagate_electronic(16);
        }
        assert!((shmd.total_population() - 1.0).abs() < 1e-6);
        let n = 3;
        for k in 0..n {
            for l in 0..n {
                assert!((shmd.rho_real[k * n + l] - shmd.rho_real[l * n + k]).abs() < 1e-14);
                assert!((shmd.rho_imag[k * n + l] + shmd.rho_imag[l * n + k]).abs() < 1e-14);
                // pure state: |ρ_kl|² = ρ_kk ρ_ll
                let abs2 = shmd.rho_real[k * n + l].powi(2) + shmd.rho_imag[k * n + l].powi(2);
                assert!((abs2 - shmd.population(k) * shmd.population(l)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_shmd_propagation_converges_with_fourth_order() {
        let reference = {
            let mut s = three_state(2.0);
            s.propagate_electronic(512);
            s
        };
        let error = |n_sub: usize| {
            let mut s = three_state(2.0);
            s.propagate_electronic(n_sub);
            (0..3)
                .map(|k| {
                    (s.amp_real[k] - reference.amp_real[k])
                        .abs()
                        .max((s.amp_imag[k] - reference.amp_imag[k]).abs())
                })
                .fold(0.0_f64, f64::max)
        };
        let (e8, e16, e32) = (error(8), error(16), error(32));
        assert!(e8 > 0.0 && e16 < e8 / 10.0, "e8 = {e8}, e16 = {e16}");
        assert!(e32 < e16 / 10.0, "e16 = {e16}, e32 = {e32}");
    }

    #[test]
    fn test_shmd_pure_state_acquires_hopping_probability_after_propagation() {
        // Regression: the amplitudes were never propagated, so ρ_kl stayed 0 for a pure
        // state and every hopping probability was exactly zero.
        let mut shmd = two_state(0.0, 0.0, 0.0, 0.05, 0.5);
        assert_eq!(shmd.total_hopping_probability(), 0.0);
        shmd.propagate_electronic(4);
        assert!(shmd.hopping_probability(1) > 0.0);
        assert!(shmd.total_hopping_probability() <= 1.0);
    }

    #[test]
    fn test_shmd_hopping_probability_equals_population_flux() {
        // c = (cos σt, sin σt): the flux out of state 0 is -dρ_00/dt = 2σ cos(σt) sin(σt),
        // so g_{0→1} = flux dt / ρ_00.
        let (sigma, dt) = (0.03, 0.5);
        let mut shmd = two_state(0.0, 0.0, 0.0, sigma, dt);
        for _ in 0..40 {
            shmd.propagate_electronic(2);
        }
        let t = 20.0;
        let (cos_t, sin_t) = ((sigma * t).cos(), (sigma * t).sin());
        let expected = 2.0 * sigma * cos_t * sin_t * dt / (cos_t * cos_t);
        let g = shmd.hopping_probability(1);
        assert!((g - expected).abs() < 1e-8, "g = {g}, expected {expected}");
        // Population is flowing *into* state 1, so there is no flux back out of it.
        shmd.current_state = 1;
        assert_eq!(shmd.hopping_probability(0), 0.0);
    }

    #[test]
    fn test_shmd_hopping_probability_includes_the_hamiltonian_coupling_term() {
        // H_01 coupling only: c = (cos ωt, -i sin ωt), ρ_01 = i cos sin,
        // g_{0→1} = (2/ħ) H_01 Im ρ_01 dt / ρ_00 = 2 ω tan(ωt) dt.
        let (h01, dt) = (0.02, 0.5);
        let omega = h01 / HBAR_EV_FS;
        let mut shmd = two_state(0.0, 0.0, h01, 0.0, dt);
        for _ in 0..30 {
            shmd.propagate_electronic(4);
        }
        let t = 15.0;
        let expected = 2.0 * omega * (omega * t).tan() * dt;
        let g = shmd.hopping_probability(1);
        assert!((g - expected).abs() < 1e-8, "g = {g}, expected {expected}");
    }

    #[test]
    fn test_shmd_fewest_switches_reproduces_quantum_populations() {
        // An ensemble of trajectories under a constant non-adiabatic coupling: the fraction
        // on state 1 must follow the quantum population sin²(σ t) -- also after the flux
        // reverses at σ t = π/2 (t ≈ 52 fs) and trajectories hop back.
        use rand::{RngExt, SeedableRng};
        let (sigma, dt, n_steps, n_traj) = (0.03, 0.25, 300, 3000);
        let energies = [0.0, 0.0];
        let mut rng = rand::rngs::StdRng::seed_from_u64(1234);
        let mut on_state_1 = 0_usize;
        for _ in 0..n_traj {
            let mut shmd = two_state(0.0, 0.0, 0.0, sigma, dt);
            for _ in 0..n_steps {
                let r: f64 = rng.random();
                let _ = shmd.attempt_hop(r, 10.0, &energies);
                shmd.propagate_electronic(2);
            }
            if shmd.current_state == 1 {
                on_state_1 += 1;
            }
        }
        let fraction = on_state_1 as f64 / n_traj as f64;
        let quantum = (sigma * n_steps as f64 * dt).sin().powi(2);
        assert!(
            (fraction - quantum).abs() < 0.035,
            "surface-hopping fraction {fraction} vs quantum population {quantum}"
        );
    }

    #[test]
    fn test_shmd_attempt_hop_rescales_velocity_and_rejects_frustrated_hops() {
        // c = (0.8, 0.6), σ_01 = 0.2, dt = 1: g_{0→1} = 2 (0.48)(0.2) / 0.64 = 0.3.
        let mut shmd = two_state(0.0, 0.0, 0.0, 0.2, 1.0);
        shmd.set_amplitudes(&[0.8, 0.6], &[0.0, 0.0]).unwrap();
        assert!((shmd.hopping_probability(1) - 0.3).abs() < 1e-12);
        let energies = [0.0, 1.0];
        // r above the hopping probability: no hop.
        assert_eq!(shmd.attempt_hop(0.5, 5.0, &energies), None);
        assert_eq!(shmd.current_state, 0);
        // Frustrated hop: ΔE = 1 eV exceeds the 0.5 eV of kinetic energy.
        assert_eq!(shmd.attempt_hop(0.1, 0.5, &energies), None);
        assert_eq!(shmd.current_state, 0);
        // Allowed hop: the velocities must be scaled by sqrt((4 - 1) / 4).
        let (state, scale) = shmd.attempt_hop(0.1, 4.0, &energies).unwrap();
        assert_eq!(state, 1);
        assert!((scale - 0.75_f64.sqrt()).abs() < 1e-12);
        assert_eq!(shmd.current_state, 1);
        assert!(
            (4.0 * scale * scale + 1.0 - 4.0).abs() < 1e-12,
            "energy conservation"
        );
    }

    #[test]
    fn test_shmd_attempt_hop_without_enough_energies_does_nothing() {
        let mut shmd = two_state(0.0, 0.0, 0.0, 0.2, 1.0);
        shmd.set_amplitudes(&[0.8, 0.6], &[0.0, 0.0]).unwrap();
        assert_eq!(shmd.attempt_hop(0.1, 4.0, &[0.0]), None);
        assert_eq!(shmd.current_state, 0);
    }

    #[test]
    fn test_shmd_velocity_rescale_factor_handles_zero_kinetic_energy() {
        let shmd = SurfaceHoppingMd::new(2, 0, 0.5).unwrap();
        // Nothing to rescale for an isoenergetic hop from rest ...
        assert_eq!(shmd.velocity_rescale_factor(0.0, 1.0, 1.0), Some(1.0));
        // ... but a uniform rescaling cannot create motion from rest, and an uphill hop is forbidden.
        assert_eq!(shmd.velocity_rescale_factor(0.0, 2.0, 1.0), None);
        assert_eq!(shmd.velocity_rescale_factor(0.0, 1.0, 2.0), None);
        assert_eq!(shmd.velocity_rescale_factor(f64::NAN, 1.0, 2.0), None);
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
