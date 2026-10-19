// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermostat kernels of [`MdSimulationPlain`](super::MdSimulationPlain).
//!
//! Two canonical-ensemble thermostats act on the velocities of the
//! self-contained driver, each in the energy unit kJ mol⁻¹ and the time unit
//! ps of the module (see [`config`](super::config)):
//!
//! * [`NoseHooverChainState`] -- a Nosé-Hoover chain (Martyna, Klein &
//!   Tuckerman, J. Chem. Phys. 97, 2635 (1992)) integrated with the
//!   time-reversible Trotter factorisation of Martyna, Tuckerman, Tobias &
//!   Klein (Mol. Phys. 87, 1117 (1996)) and a Suzuki-Yoshida decomposition of
//!   the chain propagator.
//! * [`csvr_kinetic_energy`] -- the exact stochastic velocity-rescaling step of
//!   Bussi, Donadio & Parrinello (J. Chem. Phys. 126, 014101 (2007)).

use super::stochastic::SimRng;

/// Number of thermostats in the Nosé-Hoover chain.
pub const NOSE_HOOVER_CHAIN_LENGTH: usize = 3;

/// Weights `w_j` of the three-term Suzuki-Yoshida decomposition of the chain
/// propagator, `w_1 = w_3 = 1 / (2 - 2^(1/3))` and `w_2 = 1 - 2 w_1`.
///
/// They sum to one, so the three sub-steps `w_j dt` together advance the chain
/// by exactly `dt`, and the error of the factorisation is of fourth order.
fn suzuki_yoshida_weights() -> [f64; 3] {
    let outer = 1.0 / (2.0 - 2.0_f64.cbrt());
    [outer, 1.0 - 2.0 * outer, outer]
}

/// Nosé-Hoover chain: `M = 3` coupled thermostat variables
/// `(η_k, η̇_k)` that exchange energy with the physical kinetic energy.
///
/// The equations of motion are
///
/// ```text
/// v̇_i    = F_i / m_i − η̇_1 v_i
/// Q_1 η̈_1 = 2 K − N_f k_B T − Q_1 η̇_1 η̇_2
/// Q_k η̈_k = Q_{k−1} η̇_{k−1}² − k_B T − Q_k η̇_k η̇_{k+1}      (1 < k < M)
/// Q_M η̈_M = Q_{M−1} η̇_{M−1}² − k_B T
/// ```
///
/// with `K` the kinetic energy of `N_f` degrees of freedom and the masses
/// `Q_1 = N_f k_B T τ²`, `Q_k = k_B T τ²` (`k > 1`), so that `τ` (ps) is the
/// characteristic time of the thermostat.  All energies are in kJ mol⁻¹.
///
/// The chain is propagated by [`NoseHooverChainState::propagate`], which returns
/// the factor by which every particle velocity has to be multiplied.  The
/// extended energy
///
/// ```text
/// H' = K + U + Σ_k ½ Q_k η̇_k² + N_f k_B T η_1 + k_B T Σ_{k>1} η_k
/// ```
///
/// is conserved by the exact flow; [`NoseHooverChainState::bath_energy`]
/// returns the thermostat part of it.
///
/// The chain variables `η` and `η̇` persist between steps; the masses follow
/// from `(N_f, T, τ)` and are re-derived by every call of
/// [`NoseHooverChainState::propagate`].
#[derive(Debug, Clone, PartialEq)]
pub struct NoseHooverChainState {
    /// Chain positions `η_k` (dimensionless).
    pub eta: [f64; NOSE_HOOVER_CHAIN_LENGTH],
    /// Chain velocities `η̇_k` (ps⁻¹).
    pub eta_dot: [f64; NOSE_HOOVER_CHAIN_LENGTH],
    /// Chain masses `Q_k` (kJ mol⁻¹ ps²) of the latest call of `propagate`.
    masses: [f64; NOSE_HOOVER_CHAIN_LENGTH],
    /// Degrees of freedom `N_f` of the latest call of `propagate`.
    n_dof: f64,
    /// Thermal energy `k_B T` (kJ mol⁻¹) of the latest call of `propagate`.
    kt: f64,
}

impl Default for NoseHooverChainState {
    fn default() -> Self {
        Self::new()
    }
}

impl NoseHooverChainState {
    /// Create a chain at rest (`η = η̇ = 0`).
    pub fn new() -> Self {
        Self {
            eta: [0.0; NOSE_HOOVER_CHAIN_LENGTH],
            eta_dot: [0.0; NOSE_HOOVER_CHAIN_LENGTH],
            masses: [0.0; NOSE_HOOVER_CHAIN_LENGTH],
            n_dof: 0.0,
            kt: 0.0,
        }
    }

    /// Chain masses `Q_k` (kJ mol⁻¹ ps²) used by the latest call of
    /// [`NoseHooverChainState::propagate`] (all zero before the first call).
    pub fn masses(&self) -> &[f64; NOSE_HOOVER_CHAIN_LENGTH] {
        &self.masses
    }

    /// Reset the chain variables to rest.
    pub fn reset(&mut self) {
        self.eta = [0.0; NOSE_HOOVER_CHAIN_LENGTH];
        self.eta_dot = [0.0; NOSE_HOOVER_CHAIN_LENGTH];
    }

    /// Set the masses for `n_dof` degrees of freedom at the thermal energy `kt`
    /// (kJ mol⁻¹) with the characteristic time `tau` (ps).
    fn configure(&mut self, n_dof: f64, kt: f64, tau: f64) {
        self.n_dof = n_dof;
        self.kt = kt;
        let base = kt * tau * tau;
        self.masses = [n_dof * base, base, base];
    }

    /// Driving term `G_k` of thermostat `k` (0-based) for the doubled kinetic
    /// energy `twice_kinetic = Σ m v²`.
    fn drive(&self, k: usize, twice_kinetic: f64) -> f64 {
        if k == 0 {
            (twice_kinetic - self.n_dof * self.kt) / self.masses[0]
        } else {
            let previous = self.eta_dot[k - 1];
            (self.masses[k - 1] * previous * previous - self.kt) / self.masses[k]
        }
    }

    /// Friction factor `exp(−η̇_{k+1} t / 4)` that thermostat `k` feels from its
    /// successor (the last thermostat has none).
    fn successor_friction(&self, k: usize, t: f64) -> f64 {
        match self.eta_dot.get(k + 1) {
            Some(&next) => (-0.25 * next * t).exp(),
            None => 1.0,
        }
    }

    /// Propagate the chain for the time `t` (ps) and return the velocity scale
    /// factor.
    ///
    /// * `twice_kinetic` -- `2K = Σ m v²` (kJ mol⁻¹) of the particles *before* the
    ///   call.
    /// * `n_dof` -- number of degrees of freedom `N_f`.
    /// * `kt` -- target thermal energy `k_B T` (kJ mol⁻¹), positive.
    /// * `tau` -- characteristic time of the thermostat (ps), positive.
    ///
    /// The returned factor `s` must multiply every particle velocity; the kinetic
    /// energy afterwards is `s² K`.  A velocity-Verlet step is made reversible by
    /// calling `propagate(dt / 2, …)` before the first and after the second half
    /// kick (with the kinetic energy of that moment).
    ///
    /// The chain time `t` is split into the three Suzuki-Yoshida sub-intervals
    /// `h = w_j t`.  Each sub-interval applies the symmetric (Trotter) sequence of
    /// Martyna *et al.* (1996): half an update of the chain velocities from the
    /// last thermostat down to the first (a friction factor
    /// `exp(−η̇_{k+1} h / 4)` on either side of the kick `G_k h / 2`), the particle
    /// scaling `exp(−η̇_1 h)`, the chain-position drift `η_k += η̇_k h`, and the
    /// second half update of the chain velocities from the first thermostat up to
    /// the last.  The driving terms `G_k` are recomputed from the current state at
    /// every use, including the kinetic energy that has just been scaled.
    pub fn propagate(&mut self, t: f64, twice_kinetic: f64, n_dof: f64, kt: f64, tau: f64) -> f64 {
        self.configure(n_dof, kt, tau);
        let mut twice_kinetic = twice_kinetic;
        let mut scale = 1.0_f64;
        for weight in suzuki_yoshida_weights() {
            let h = weight * t;

            // Chain velocities, last thermostat first.
            for k in (0..NOSE_HOOVER_CHAIN_LENGTH).rev() {
                let friction = self.successor_friction(k, h);
                self.eta_dot[k] *= friction;
                self.eta_dot[k] += 0.5 * h * self.drive(k, twice_kinetic);
                self.eta_dot[k] *= friction;
            }

            // Particle velocities scale with exp(−η̇_1 h).
            let factor = (-self.eta_dot[0] * h).exp();
            scale *= factor;
            twice_kinetic *= factor * factor;

            // Chain positions.
            for k in 0..NOSE_HOOVER_CHAIN_LENGTH {
                self.eta[k] += self.eta_dot[k] * h;
            }

            // Chain velocities, first thermostat first, with the updated state.
            for k in 0..NOSE_HOOVER_CHAIN_LENGTH {
                let friction = self.successor_friction(k, h);
                self.eta_dot[k] *= friction;
                self.eta_dot[k] += 0.5 * h * self.drive(k, twice_kinetic);
                self.eta_dot[k] *= friction;
            }
        }
        scale
    }

    /// Thermostat part of the extended energy (kJ mol⁻¹):
    /// `Σ_k ½ Q_k η̇_k² + N_f k_B T η_1 + k_B T Σ_{k>1} η_k`.
    ///
    /// Add the kinetic and potential energy of the particles to obtain the
    /// conserved quantity `H'`.  The masses are those of the latest call of
    /// [`NoseHooverChainState::propagate`]; before any call the result is zero.
    pub fn bath_energy(&self) -> f64 {
        let mut energy = self.n_dof * self.kt * self.eta[0];
        for k in 0..NOSE_HOOVER_CHAIN_LENGTH {
            energy += 0.5 * self.masses[k] * self.eta_dot[k] * self.eta_dot[k];
            if k > 0 {
                energy += self.kt * self.eta[k];
            }
        }
        energy
    }
}

/// Kinetic energy after one exact step of the canonical-sampling
/// velocity-rescaling (CSVR) thermostat of Bussi, Donadio & Parrinello (2007).
///
/// The kinetic energy `K` of `N_f` degrees of freedom obeys the stochastic
/// differential equation
///
/// ```text
/// dK = (K̄ − K) dt/τ + 2 sqrt(K K̄ / N_f) dW / sqrt(τ),          K̄ = N_f k_B T / 2,
/// ```
///
/// whose stationary law is the canonical gamma distribution with mean `K̄` and
/// variance `2 K̄² / N_f`.  Integrating it exactly over `dt` gives (eq. A7 of
/// the paper)
///
/// ```text
/// K' = c K + (1 − c) (K̄ / N_f) (R_1² + S) + 2 R_1 sqrt(c (1 − c) K K̄ / N_f),
/// ```
///
/// with `c = exp(−dt/τ)`, `R_1 ~ N(0, 1)` and `S ~ χ²(N_f − 1)` independent.
/// This is evaluated in the equivalent, manifestly non-negative form
/// `K' = (sqrt(c K) + R_1 sqrt((1 − c) K̄ / N_f))² + (1 − c) (K̄ / N_f) S`.
///
/// * `kinetic` -- current kinetic energy `K` (kJ mol⁻¹); a system at rest
///   (`K <= 0`) cannot be rescaled and is returned unchanged.
/// * `target_kinetic` -- `K̄ = N_f k_B T / 2` (kJ mol⁻¹).
/// * `n_dof` -- number of degrees of freedom `N_f`, at least one.
/// * `dt`, `tau` -- time step and coupling time (ps).  `tau <= 0` couples
///   instantaneously (`c = 0`, the kinetic energy is redrawn from the
///   canonical distribution); `dt <= 0` leaves `K` unchanged.
///
/// Every velocity is to be multiplied by `sqrt(K' / K)`.
pub fn csvr_kinetic_energy(
    kinetic: f64,
    target_kinetic: f64,
    n_dof: f64,
    dt: f64,
    tau: f64,
    rng: &mut SimRng,
) -> f64 {
    let rescalable = kinetic.is_finite() && kinetic > 0.0 && n_dof >= 1.0 && dt > 0.0;
    if !rescalable {
        return kinetic;
    }
    let decay = if tau > 0.0 { (-dt / tau).exp() } else { 0.0 };
    let noise_energy = (1.0 - decay) * target_kinetic / n_dof;
    let r1 = rng.normal();
    let s = rng.chi_squared(n_dof - 1.0);
    let carried = (decay * kinetic).sqrt() + r1 * noise_energy.sqrt();
    carried * carried + noise_energy * s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suzuki_yoshida_weights_sum_to_one() {
        let w = suzuki_yoshida_weights();
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-15, "sum = {sum}");
        assert!(
            (w[0] - 1.351_207_191_959_657_8).abs() < 1e-14,
            "w1 = {}",
            w[0]
        );
        assert!(
            (w[1] + 1.702_414_383_919_315_6).abs() < 1e-14,
            "w2 = {}",
            w[1]
        );
        assert_eq!(w[0], w[2]);
    }

    #[test]
    fn test_chain_masses_follow_the_definition() {
        let mut chain = NoseHooverChainState::new();
        assert_eq!(chain.masses(), &[0.0; 3]);
        let (n_dof, kt, tau) = (30.0, 2.5, 0.2);
        chain.propagate(0.001, n_dof * kt, n_dof, kt, tau);
        let q = chain.masses();
        assert!((q[0] - n_dof * kt * tau * tau).abs() < 1e-15);
        assert!((q[1] - kt * tau * tau).abs() < 1e-15);
        assert!((q[2] - kt * tau * tau).abs() < 1e-15);
    }

    #[test]
    fn test_chain_at_target_kinetic_energy_is_stationary_to_first_order() {
        // 2K = N_f kT: the first driving term vanishes, the chain stays at rest and
        // does not rescale the particles.
        let mut chain = NoseHooverChainState::new();
        let (n_dof, kt) = (60.0, 1.2);
        let scale = chain.propagate(0.01, n_dof * kt, n_dof, kt, 0.1);
        assert!((scale - 1.0).abs() < 1e-12, "scale = {scale}");
        assert!(chain.eta_dot[0].abs() < 1e-12);
        // The outer thermostats feel the −kT/Q term only through η̇_1² = 0.
        assert!(chain.eta_dot[1] < 0.0, "second thermostat is driven by -kT");
    }

    #[test]
    fn test_hot_system_is_cooled_and_cold_system_is_heated() {
        let (n_dof, kt, tau) = (90.0, 1.0, 0.1);
        let mut hot = NoseHooverChainState::new();
        let mut scale_hot = 1.0;
        let mut twice_k = 2.0 * n_dof * kt; // twice the target
        for _ in 0..50 {
            let s = hot.propagate(0.002, twice_k, n_dof, kt, tau);
            scale_hot *= s;
            twice_k *= s * s;
        }
        assert!(scale_hot < 1.0, "hot system must cool: scale = {scale_hot}");

        let mut cold = NoseHooverChainState::new();
        let mut scale_cold = 1.0;
        let mut twice_k = 0.5 * n_dof * kt; // half the target
        for _ in 0..50 {
            let s = cold.propagate(0.002, twice_k, n_dof, kt, tau);
            scale_cold *= s;
            twice_k *= s * s;
        }
        assert!(
            scale_cold > 1.0,
            "cold system must heat: scale = {scale_cold}"
        );
    }

    #[test]
    fn test_reset_returns_the_chain_to_rest() {
        let mut chain = NoseHooverChainState::new();
        chain.propagate(0.01, 500.0, 30.0, 1.0, 0.1);
        assert!(chain.eta_dot.iter().any(|&v| v != 0.0));
        chain.reset();
        assert_eq!(chain.eta, [0.0; 3]);
        assert_eq!(chain.eta_dot, [0.0; 3]);
    }

    #[test]
    fn test_csvr_returns_input_when_it_cannot_rescale() {
        let mut rng = SimRng::new(1);
        assert_eq!(
            csvr_kinetic_energy(0.0, 10.0, 30.0, 0.01, 0.1, &mut rng),
            0.0
        );
        assert_eq!(
            csvr_kinetic_energy(5.0, 10.0, 30.0, 0.0, 0.1, &mut rng),
            5.0
        );
        assert_eq!(
            csvr_kinetic_energy(5.0, 10.0, 0.0, 0.01, 0.1, &mut rng),
            5.0
        );
        assert!(csvr_kinetic_energy(f64::NAN, 10.0, 30.0, 0.01, 0.1, &mut rng).is_nan());
    }

    #[test]
    fn test_csvr_kinetic_energy_is_never_negative() {
        let mut rng = SimRng::new(11);
        for _ in 0..10_000 {
            let k = csvr_kinetic_energy(1e-6, 50.0, 3.0, 0.01, 0.02, &mut rng);
            assert!(k >= 0.0, "negative kinetic energy {k}");
        }
    }

    #[test]
    fn test_csvr_stationary_distribution_is_the_canonical_gamma() {
        // Iterating the exact step must reproduce mean K̄ and variance 2 K̄²/N_f.
        let n_dof = 12.0_f64;
        let target = 0.5 * n_dof * 1.7; // K̄ for k_B T = 1.7
        let mut rng = SimRng::new(5);
        let mut k = target;
        let (mut sum, mut sum_sq) = (0.0, 0.0);
        let steps = 400_000_usize;
        for _ in 0..steps {
            k = csvr_kinetic_energy(k, target, n_dof, 0.05, 0.1, &mut rng);
            sum += k;
            sum_sq += k * k;
        }
        let mean = sum / steps as f64;
        let var = sum_sq / steps as f64 - mean * mean;
        let expected_var = 2.0 * target * target / n_dof;
        assert!(
            (mean - target).abs() < 0.02 * target,
            "mean K = {mean}, expected {target}"
        );
        assert!(
            (var - expected_var).abs() < 0.08 * expected_var,
            "var K = {var}, expected {expected_var}"
        );
    }

    #[test]
    fn test_csvr_with_instant_coupling_redraws_from_the_gamma_law() {
        // tau = 0: K' = K̄/N_f (R² + χ²(N_f − 1)) = K̄/N_f χ²(N_f), independent of K.
        let n_dof = 20.0_f64;
        let target = 0.5 * n_dof * 2.0;
        let mut rng = SimRng::new(9);
        let samples = 100_000_usize;
        let mean: f64 = (0..samples)
            .map(|_| csvr_kinetic_energy(123.0, target, n_dof, 0.01, 0.0, &mut rng))
            .sum::<f64>()
            / samples as f64;
        assert!(
            (mean - target).abs() < 0.02 * target,
            "mean = {mean}, expected {target}"
        );
    }
}
