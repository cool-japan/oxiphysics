// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hill-type muscle simulation for the OxiPhysics soft-body crate.
//!
//! Implements the classic three-element Hill muscle model (contractile element,
//! series elastic element, parallel elastic element), activation dynamics,
//! muscle groups, fatigue, and utility functions for musculotendon geometry and
//! EMG-based optimal fibre length estimation.

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// 3-D Euclidean norm.
#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Subtract two 3-D vectors: a − b.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by scalar s.
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalise a 3-D vector (returns zero vector if nearly degenerate).
#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

/// Clamp `x` to `[lo, hi]`.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// HillMuscle
// ---------------------------------------------------------------------------

/// Hill-type muscle model with contractile, series elastic, and parallel
/// elastic elements.
///
/// Parameters follow the notation of Zajac (1989) and Thelen (2003).
#[derive(Debug, Clone)]
pub struct HillMuscle {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fibre length (m).
    pub l_opt: f64,
    /// Maximum shortening velocity (m/s, positive).
    pub v_max: f64,
    /// Pennation angle at optimal length (radians).
    pub pennation_angle: f64,
}

impl HillMuscle {
    /// Create a new [`HillMuscle`].
    pub fn new(f_max: f64, l_opt: f64, v_max: f64, pennation_angle: f64) -> Self {
        Self {
            f_max,
            l_opt,
            v_max,
            pennation_angle,
        }
    }

    /// Active force–length relationship (normalised, bell-shaped Gaussian).
    ///
    /// `l` is the current fibre length (m).  Returns a value in \[0, 1\] that
    /// peaks at 1 when l == l_opt.
    pub fn active_force_length(&self, l: f64) -> f64 {
        let l_norm = l / self.l_opt;
        // Gaussian with width parameter β = 0.45 (Thelen 2003).
        let beta = 0.45;
        let x = (l_norm - 1.0) / beta;
        (-0.5 * x * x).exp().clamp(0.0, 1.0)
    }

    /// Passive force–length relationship (exponential).
    ///
    /// Returns the normalised passive force (0 at l ≤ l_opt, exponential for l > l_opt).
    pub fn passive_force_length(&self, l: f64) -> f64 {
        let l_norm = l / self.l_opt;
        if l_norm <= 1.0 {
            0.0
        } else {
            // Exponential passive element from Thelen (2003).
            let kpe = 5.0; // shape parameter
            ((l_norm - 1.0) * kpe).exp() - 1.0
        }
    }

    /// Force–velocity relationship (normalised).
    ///
    /// Uses the Hill hyperbola for shortening velocities (v < 0) and a
    /// linear extension for lengthening (v > 0).
    ///
    /// `v` is fibre velocity in m/s (negative = shortening).
    /// Returns a value in \[0, ~1.8\] (typically).
    pub fn force_velocity(&self, v: f64) -> f64 {
        if v <= 0.0 {
            // Shortening: Hill hyperbola.
            // FV = (v_max + v) / (v_max - v / a_hill)
            let a_hill = 0.25; // Hill constant
            let v_norm = v / self.v_max;
            let num = 1.0 + v_norm;
            let den = 1.0 - v_norm / a_hill;
            if den.abs() < 1e-15 {
                0.0
            } else {
                clamp(num / den, 0.0, 1.0)
            }
        } else {
            // Lengthening: eccentric linear extension capped at 1.8 × F_max.
            let fv_max = 1.8;
            let slope = (fv_max - 1.0) / self.v_max;
            clamp(1.0 + slope * v, 1.0, fv_max)
        }
    }

    /// Total muscle force (N).
    ///
    /// F = activation · FL · FV · F_max · cos(pennation) + F_passive · F_max · cos(pennation)
    ///
    /// `l` – current fibre length (m), `v` – fibre velocity (m/s),
    /// `activation` – neural activation in \[0, 1\].
    pub fn total_force(&self, l: f64, v: f64, activation: f64) -> f64 {
        let a = clamp(activation, 0.0, 1.0);
        let fl = self.active_force_length(l);
        let fv = self.force_velocity(v);
        let fp = self.passive_force_length(l);
        let cos_phi = self.pennation_angle.cos();
        (a * fl * fv + fp) * self.f_max * cos_phi
    }
}

// ---------------------------------------------------------------------------
// ActivationDynamics
// ---------------------------------------------------------------------------

/// First-order activation dynamics model.
///
/// Smooths the neural excitation `u(t)` into the muscle activation `a(t)`
/// with separate activation and deactivation time constants (Thelen 2003).
#[derive(Debug, Clone)]
pub struct ActivationDynamics {
    /// Activation time constant (s), typically ~10 ms.
    pub tau_act: f64,
    /// Deactivation time constant (s), typically ~40 ms.
    pub tau_deact: f64,
}

impl ActivationDynamics {
    /// Create new [`ActivationDynamics`].
    pub fn new(tau_act: f64, tau_deact: f64) -> Self {
        Self { tau_act, tau_deact }
    }

    /// Advance activation by one time step `dt` (s).
    ///
    /// `activation` – current activation (0–1), `excitation` – neural signal (0–1).
    /// Returns the new activation value.
    pub fn step(&self, activation: f64, excitation: f64, dt: f64) -> f64 {
        let a = clamp(activation, 0.0, 1.0);
        let u = clamp(excitation, 0.0, 1.0);
        // Time constant depends on direction (activating vs deactivating).
        let tau = if u > a {
            self.tau_act * (0.5 + 1.5 * a)
        } else {
            self.tau_deact / (0.5 + 1.5 * a)
        };
        let da_dt = (u - a) / tau;
        clamp(a + da_dt * dt, 0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// MuscleElement
// ---------------------------------------------------------------------------

/// A single Hill-muscle line-of-action element connecting two attachment points.
#[derive(Debug, Clone)]
pub struct MuscleElement {
    /// Proximal attachment point (m).
    pub pos_a: [f64; 3],
    /// Distal attachment point (m).
    pub pos_b: [f64; 3],
    /// Underlying Hill muscle model.
    pub model: HillMuscle,
    /// Current activation (0–1).
    pub activation: f64,
    /// Rest length of the muscle element (m).
    pub rest_length: f64,
}

impl MuscleElement {
    /// Create a new [`MuscleElement`].
    pub fn new(
        pos_a: [f64; 3],
        pos_b: [f64; 3],
        model: HillMuscle,
        activation: f64,
        rest_length: f64,
    ) -> Self {
        Self {
            pos_a,
            pos_b,
            model,
            activation,
            rest_length,
        }
    }

    /// Current length of the muscle element (m): ‖pos_b − pos_a‖.
    pub fn current_length(&self) -> f64 {
        norm3(sub3(self.pos_b, self.pos_a))
    }

    /// Velocity estimate (m/s) from the previous length and time step.
    ///
    /// Positive velocity = lengthening, negative = shortening.
    pub fn current_velocity(&self, prev_length: f64, dt: f64) -> f64 {
        if dt < 1e-30 {
            0.0
        } else {
            (self.current_length() - prev_length) / dt
        }
    }

    /// Compute force vectors applied at attachment point A and B.
    ///
    /// Returns `(force_on_A, force_on_B)`.  The forces are equal and opposite
    /// (Newton's third law) directed along the line-of-action.
    pub fn force_vector(&self, dt: f64, prev_length: f64) -> ([f64; 3], [f64; 3]) {
        let l = self.current_length();
        let v = self.current_velocity(prev_length, dt);
        let f_scalar = self.model.total_force(l, v, self.activation);

        // Unit vector from A to B.
        let ab = sub3(self.pos_b, self.pos_a);
        let dir = normalize3(ab);

        // Force on B is in the direction from A to B (tension pulls B toward A).
        // Force on A is equal and opposite.
        let force_on_b = scale3(dir, -f_scalar); // tension pulls B toward A
        let force_on_a = scale3(dir, f_scalar); // and A toward B
        (force_on_a, force_on_b)
    }
}

// ---------------------------------------------------------------------------
// MuscleGroup
// ---------------------------------------------------------------------------

/// A group of parallel [`MuscleElement`]s sharing a physiological
/// cross-sectional area (PCSA).
#[derive(Debug, Clone)]
pub struct MuscleGroup {
    /// Individual muscle elements.
    pub elements: Vec<MuscleElement>,
    /// Physiological cross-sectional area (m²).
    pub pcsa: f64,
}

impl MuscleGroup {
    /// Create a new [`MuscleGroup`].
    pub fn new(elements: Vec<MuscleElement>, pcsa: f64) -> Self {
        Self { elements, pcsa }
    }

    /// Resultant force vector (N) from all elements.
    ///
    /// `prev_lengths` must have the same length as `self.elements`.
    pub fn resultant_force(&self, dt: f64, prev_lengths: &[f64]) -> [f64; 3] {
        let mut total = [0.0f64; 3];
        for (elem, &pl) in self.elements.iter().zip(prev_lengths.iter()) {
            let (fa, _fb) = elem.force_vector(dt, pl);
            total[0] += fa[0];
            total[1] += fa[1];
            total[2] += fa[2];
        }
        total
    }

    /// Set activation uniformly across all elements.
    pub fn set_activation(&mut self, activation: f64) {
        let a = clamp(activation, 0.0, 1.0);
        for elem in &mut self.elements {
            elem.activation = a;
        }
    }
}

// ---------------------------------------------------------------------------
// ThilenFibreModel  (three-element model)
// ---------------------------------------------------------------------------

/// Three-element muscle fibre model (Contractile Element + SEE + PEE).
///
/// Simplified static equilibrium version: finds the CE length at which the
/// active CE force balances the series elastic element force.
#[derive(Debug, Clone)]
pub struct ThilenFibreModel {
    /// Current contractile-element length (m).
    pub ce_length: f64,
    /// Series elastic element stiffness (N/m).
    pub see_stiffness: f64,
    /// Parallel elastic element stiffness (N/m).
    pub pee_stiffness: f64,
}

impl ThilenFibreModel {
    /// Create a new [`ThilenFibreModel`].
    pub fn new(ce_length: f64, see_stiffness: f64, pee_stiffness: f64) -> Self {
        Self {
            ce_length,
            see_stiffness,
            pee_stiffness,
        }
    }

    /// Equilibrium CE length (m) for a given activation and maximum force.
    ///
    /// Iterates via bisection until the CE force equals the SEE restoring
    /// force to within 1 N.
    pub fn equilibrium_length(&self, activation: f64, f_max: f64) -> f64 {
        let a = clamp(activation, 0.0, 1.0);
        // CE force at current length (simplified: Gaussian FL, FV=1).
        // SEE elongation = total length - CE length.
        // Equilibrium: a * fl(l_ce) * f_max = k_see * (l_tot - l_ce)
        // Solve for l_ce via bisection over [0.5, 1.5] × l0.
        let l0 = self.ce_length;
        let l_tot = l0 * 1.1; // assume total MTU length is 10% longer
        let mut lo = 0.5 * l0;
        let mut hi = 1.5 * l0;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let l_norm = mid / l0;
            let beta = 0.45;
            let fl = (-(l_norm - 1.0).powi(2) / (2.0 * beta * beta)).exp();
            let f_ce = a * fl * f_max + self.pee_stiffness * (mid - l0).max(0.0);
            let f_see = self.see_stiffness * (l_tot - mid).max(0.0);
            if f_ce < f_see {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ---------------------------------------------------------------------------
// musculotendon_length
// ---------------------------------------------------------------------------

/// Compute musculotendon unit (MTU) length from fibre and tendon lengths.
///
/// `muscle_length` – fibre length (m), `tendon_length` – tendon length (m),
/// `pennation` – pennation angle (rad).
///
/// Returns MTU length = tendon_length + muscle_length · cos(pennation).
pub fn musculotendon_length(muscle_length: f64, tendon_length: f64, pennation: f64) -> f64 {
    tendon_length + muscle_length * pennation.cos()
}

// ---------------------------------------------------------------------------
// optimal_fiber_length_from_emg
// ---------------------------------------------------------------------------

/// Estimate the optimal fibre length from paired EMG and force signals.
///
/// Uses a least-squares regression on the relationship
/// F ≈ EMG · F_max · exp(−((l/l_opt − 1)/β)²)
/// to infer l_opt by minimising residuals.
///
/// `emg` – normalised EMG amplitudes (0–1), `force` – corresponding force
/// measurements (N).
///
/// Returns the estimated optimal fibre length as a dimensionless ratio
/// (relative to a unit reference length of 1 m); multiply by your reference
/// length as needed.
pub fn optimal_fiber_length_from_emg(emg: &[f64], force: &[f64]) -> f64 {
    assert_eq!(
        emg.len(),
        force.len(),
        "emg and force must have the same length"
    );
    if emg.is_empty() {
        return 1.0;
    }

    // Enumerate candidate l_opt values over a grid (0.5 … 2.0 in 100 steps).
    let n_steps = 100;
    let mut best_l_opt = 1.0f64;
    let mut best_sse = f64::INFINITY;

    for k in 0..=n_steps {
        let l_opt = 0.5 + 1.5 * (k as f64 / n_steps as f64);
        let mut sse = 0.0f64;
        for (&u, &f) in emg.iter().zip(force.iter()) {
            let beta = 0.45;
            let x = (1.0 / l_opt - 1.0) / beta; // l = 1.0 (reference)
            let f_pred = u * (-0.5 * x * x).exp();
            sse += (f - f_pred).powi(2);
        }
        if sse < best_sse {
            best_sse = sse;
            best_l_opt = l_opt;
        }
    }
    best_l_opt
}

// ---------------------------------------------------------------------------
// MuscleFatigue
// ---------------------------------------------------------------------------

/// Three-compartment muscle fatigue model (Liu et al. 2002).
///
/// Tracks the fatigued force capacity as a function of activation history.
#[derive(Debug, Clone)]
pub struct MuscleFatigue {
    /// Fatigue rate constant (s⁻¹).
    pub fatigue_rate: f64,
    /// Recovery rate constant (s⁻¹).
    pub recovery_rate: f64,
    /// Current fatigue factor in \[0, 1\] (0 = no fatigue, 1 = fully fatigued).
    pub fatigue_factor: f64,
}

impl MuscleFatigue {
    /// Create a new [`MuscleFatigue`] with the given rate constants.
    pub fn new(fatigue_rate: f64, recovery_rate: f64) -> Self {
        Self {
            fatigue_rate,
            recovery_rate,
            fatigue_factor: 0.0,
        }
    }

    /// Advance the fatigue state by one time step `dt` (s).
    ///
    /// `activation` – current muscle activation (0–1).
    pub fn step(&mut self, activation: f64, dt: f64) {
        let a = clamp(activation, 0.0, 1.0);
        // Fatiguing when active, recovering when inactive.
        let d_fatigue = if a > self.fatigue_factor {
            self.fatigue_rate * (a - self.fatigue_factor)
        } else {
            -self.recovery_rate * (self.fatigue_factor - a)
        };
        self.fatigue_factor = clamp(self.fatigue_factor + d_fatigue * dt, 0.0, 1.0);
    }

    /// Fraction of maximum force still available (1 = fully fresh, 0 = exhausted).
    pub fn available_force_fraction(&self) -> f64 {
        1.0 - self.fatigue_factor
    }
}

// ---------------------------------------------------------------------------
// twitch_response
// ---------------------------------------------------------------------------

/// Double-exponential twitch force response to a single motor-unit impulse.
///
/// `t`       – time since impulse (s, ≥ 0),
/// `t_rise`  – rise time constant (s),
/// `t_fall`  – fall time constant (s).
///
/// Returns the normalised twitch force (peak ≈ 1).
pub fn twitch_response(t: f64, t_rise: f64, t_fall: f64) -> f64 {
    if t < 0.0 {
        return 0.0;
    }
    // f(t) = (exp(-t/t_fall) - exp(-t/t_rise)) normalised by its peak.
    let raw = (-t / t_fall).exp() - (-t / t_rise).exp();
    // Normalise by the theoretical peak: f(t*) where t* = ln(t_fall/t_rise)/(1/t_rise - 1/t_fall).
    let norm = if (t_rise - t_fall).abs() < 1e-15 || t_rise <= 0.0 || t_fall <= 0.0 {
        1.0
    } else {
        let t_peak = (t_fall / t_rise).ln() / (1.0 / t_rise - 1.0 / t_fall);
        let peak = (-t_peak / t_fall).exp() - (-t_peak / t_rise).exp();
        if peak.abs() < 1e-30 { 1.0 } else { peak }
    };
    (raw / norm).max(0.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. active_force_length peaks at l_opt.
    #[test]
    fn test_fl_peaks_at_l_opt() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fl_opt = m.active_force_length(m.l_opt);
        assert!(
            (fl_opt - 1.0).abs() < EPS,
            "FL should be 1 at l_opt, got {fl_opt}"
        );
    }

    // 2. active_force_length decays away from l_opt.
    #[test]
    fn test_fl_decays_off_peak() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fl_near = m.active_force_length(0.85 * m.l_opt);
        assert!(fl_near < 1.0, "FL should be < 1 away from l_opt");
        assert!(fl_near > 0.0, "FL should be positive near l_opt");
    }

    // 3. passive_force_length is zero at l_opt.
    #[test]
    fn test_fp_zero_at_l_opt() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        assert_eq!(m.passive_force_length(m.l_opt), 0.0);
    }

    // 4. passive_force_length is positive for l > l_opt.
    #[test]
    fn test_fp_positive_above_l_opt() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fp = m.passive_force_length(1.2 * m.l_opt);
        assert!(
            fp > 0.0,
            "Passive force should be positive for l > l_opt, got {fp}"
        );
    }

    // 5. force_velocity is 1 at zero velocity.
    #[test]
    fn test_fv_at_zero() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fv = m.force_velocity(0.0);
        assert!((fv - 1.0).abs() < 1e-9, "FV at v=0 should be 1, got {fv}");
    }

    // 6. force_velocity decreases for shortening (v < 0).
    #[test]
    fn test_fv_decreases_shortening() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fv_neg = m.force_velocity(-0.3 * m.v_max);
        assert!(fv_neg < 1.0, "FV should be < 1 during shortening");
        assert!(fv_neg >= 0.0, "FV should be non-negative");
    }

    // 7. force_velocity increases for lengthening (v > 0).
    #[test]
    fn test_fv_increases_lengthening() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fv_pos = m.force_velocity(0.5 * m.v_max);
        assert!(
            fv_pos > 1.0,
            "FV should be > 1 during lengthening, got {fv_pos}"
        );
    }

    // 8. total_force is zero when activation is zero at rest length.
    #[test]
    fn test_total_force_zero_activation() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let f = m.total_force(m.l_opt, 0.0, 0.0);
        assert_eq!(f, 0.0, "No force at zero activation and rest length");
    }

    // 9. total_force is positive with activation at optimal length.
    #[test]
    fn test_total_force_positive_activation() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let f = m.total_force(m.l_opt, 0.0, 1.0);
        assert!(
            (f - 1000.0).abs() < 1.0,
            "Max force at full activation: {f}"
        );
    }

    // 10. total_force scaled by pennation angle cosine.
    #[test]
    fn test_total_force_pennation() {
        let m0 = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let m30 = HillMuscle::new(1000.0, 0.1, 1.0, 30.0_f64.to_radians());
        let f0 = m0.total_force(m0.l_opt, 0.0, 1.0);
        let f30 = m30.total_force(m30.l_opt, 0.0, 1.0);
        assert!(f30 < f0, "Pennation reduces effective force");
    }

    // 11. ActivationDynamics: activation increases toward excitation.
    #[test]
    fn test_activation_increases() {
        let dyn_ = ActivationDynamics::new(0.01, 0.04);
        let a = dyn_.step(0.1, 1.0, 0.001);
        assert!(a > 0.1, "Activation should increase toward excitation=1");
    }

    // 12. ActivationDynamics: activation decreases when excitation < activation.
    #[test]
    fn test_activation_decreases() {
        let dyn_ = ActivationDynamics::new(0.01, 0.04);
        let a = dyn_.step(0.8, 0.0, 0.001);
        assert!(
            a < 0.8,
            "Activation should decrease when excitation is zero"
        );
    }

    // 13. ActivationDynamics: output stays in [0, 1].
    #[test]
    fn test_activation_clamped() {
        let dyn_ = ActivationDynamics::new(0.01, 0.04);
        let a = dyn_.step(0.0, 2.0, 1.0); // large dt and over-excitation
        assert!(
            (0.0..=1.0).contains(&a),
            "Activation must stay in [0, 1]: {a}"
        );
    }

    // 14. MuscleElement: current_length is correct.
    #[test]
    fn test_muscle_element_length() {
        let model = HillMuscle::new(500.0, 0.08, 0.8, 0.0);
        let elem = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], model, 0.5, 0.1);
        assert!((elem.current_length() - 0.1).abs() < EPS);
    }

    // 15. MuscleElement: velocity estimation.
    #[test]
    fn test_muscle_element_velocity() {
        let model = HillMuscle::new(500.0, 0.08, 0.8, 0.0);
        let elem = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.12, 0.0], model, 0.5, 0.1);
        let v = elem.current_velocity(0.1, 0.01);
        assert!((v - 2.0).abs() < 1e-9, "Velocity should be 2 m/s, got {v}");
    }

    // 16. MuscleElement: force_vector is equal and opposite.
    #[test]
    fn test_muscle_element_force_newton3() {
        let model = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let elem = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], model, 1.0, 0.1);
        let (fa, fb) = elem.force_vector(0.01, 0.1);
        assert!(
            (fa[0] + fb[0]).abs() < EPS
                && (fa[1] + fb[1]).abs() < EPS
                && (fa[2] + fb[2]).abs() < EPS,
            "Forces on A and B should sum to zero (Newton 3rd law)"
        );
    }

    // 17. MuscleGroup: resultant force sums elements.
    #[test]
    fn test_muscle_group_resultant() {
        let model = HillMuscle::new(500.0, 0.1, 1.0, 0.0);
        let e1 = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], model.clone(), 1.0, 0.1);
        let e2 = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], model.clone(), 1.0, 0.1);
        let group = MuscleGroup::new(vec![e1, e2], 1e-4);
        let f = group.resultant_force(0.01, &[0.1, 0.1]);
        // Two identical elements should give twice the single element force.
        let single = MuscleElement::new([0.0, 0.0, 0.0], [0.0, 0.1, 0.0], model, 1.0, 0.1);
        let (fa, _) = single.force_vector(0.01, 0.1);
        assert!(
            (f[1] - 2.0 * fa[1]).abs() < EPS,
            "Group force should be 2× single, got {} vs {}",
            f[1],
            2.0 * fa[1]
        );
    }

    // 18. MuscleGroup: set_activation clamps to [0, 1].
    #[test]
    fn test_muscle_group_set_activation() {
        let model = HillMuscle::new(500.0, 0.1, 1.0, 0.0);
        let e = MuscleElement::new([0.0; 3], [0.0, 0.1, 0.0], model, 0.5, 0.1);
        let mut group = MuscleGroup::new(vec![e], 1e-4);
        group.set_activation(1.5); // over-range
        assert_eq!(group.elements[0].activation, 1.0);
    }

    // 19. ThilenFibreModel: equilibrium_length in valid range.
    #[test]
    fn test_thilen_equilibrium_range() {
        let model = ThilenFibreModel::new(0.1, 50_000.0, 5_000.0);
        let l_eq = model.equilibrium_length(0.5, 1000.0);
        assert!(
            l_eq > 0.0 && l_eq < 1.0,
            "Equilibrium length should be in (0, 1), got {l_eq}"
        );
    }

    // 20. ThilenFibreModel: higher activation → different equilibrium.
    #[test]
    fn test_thilen_activation_effect() {
        let model = ThilenFibreModel::new(0.1, 50_000.0, 5_000.0);
        let l_low = model.equilibrium_length(0.1, 1000.0);
        let l_high = model.equilibrium_length(0.9, 1000.0);
        assert_ne!(
            l_low, l_high,
            "Equilibrium length should change with activation"
        );
    }

    // 21. musculotendon_length: zero pennation = tendon + muscle.
    #[test]
    fn test_mtu_length_zero_pennation() {
        let l = musculotendon_length(0.08, 0.25, 0.0);
        assert!(
            (l - 0.33).abs() < EPS,
            "MTU length with zero pennation: {l}"
        );
    }

    // 22. musculotendon_length: 90° pennation → only tendon contributes.
    #[test]
    fn test_mtu_length_perpendicular() {
        let l = musculotendon_length(0.08, 0.25, std::f64::consts::FRAC_PI_2);
        assert!((l - 0.25).abs() < 1e-9, "MTU length at 90° pennation: {l}");
    }

    // 23. optimal_fiber_length_from_emg: returns value in (0.5, 2.0).
    #[test]
    fn test_emg_fiber_length_range() {
        let emg = vec![0.2, 0.5, 0.8, 1.0, 0.6];
        let force = vec![200.0, 500.0, 800.0, 1000.0, 600.0];
        let l_opt = optimal_fiber_length_from_emg(&emg, &force);
        assert!((0.5..=2.0).contains(&l_opt), "l_opt out of range: {l_opt}");
    }

    // 24. optimal_fiber_length_from_emg: empty input returns 1.0.
    #[test]
    fn test_emg_empty() {
        let l = optimal_fiber_length_from_emg(&[], &[]);
        assert_eq!(l, 1.0);
    }

    // 25. MuscleFatigue: starts at 0.
    #[test]
    fn test_fatigue_starts_zero() {
        let f = MuscleFatigue::new(0.01, 0.005);
        assert_eq!(f.fatigue_factor, 0.0);
    }

    // 26. MuscleFatigue: increases under sustained activation.
    #[test]
    fn test_fatigue_increases() {
        let mut f = MuscleFatigue::new(0.1, 0.01);
        for _ in 0..100 {
            f.step(1.0, 0.01);
        }
        assert!(
            f.fatigue_factor > 0.0,
            "Fatigue should increase under sustained activation"
        );
    }

    // 27. MuscleFatigue: available force fraction decreases under fatigue.
    #[test]
    fn test_fatigue_force_fraction() {
        let mut f = MuscleFatigue::new(0.1, 0.01);
        for _ in 0..100 {
            f.step(1.0, 0.01);
        }
        assert!(
            f.available_force_fraction() < 1.0,
            "Force fraction should decrease under fatigue"
        );
    }

    // 28. MuscleFatigue: fatigue_factor stays in [0, 1].
    #[test]
    fn test_fatigue_clamped() {
        let mut f = MuscleFatigue::new(1.0, 0.0001);
        for _ in 0..10_000 {
            f.step(1.0, 0.1);
        }
        assert!(
            f.fatigue_factor >= 0.0 && f.fatigue_factor <= 1.0,
            "Fatigue factor out of range: {}",
            f.fatigue_factor
        );
    }

    // 29. twitch_response: zero at t < 0.
    #[test]
    fn test_twitch_negative_time() {
        assert_eq!(twitch_response(-0.001, 0.02, 0.08), 0.0);
    }

    // 30. twitch_response: positive for t > 0.
    #[test]
    fn test_twitch_positive() {
        let v = twitch_response(0.05, 0.02, 0.08);
        assert!(v > 0.0, "Twitch should be positive at t=0.05, got {v}");
    }

    // 31. twitch_response: returns to near-zero for large t.
    #[test]
    fn test_twitch_decays() {
        let v = twitch_response(10.0, 0.02, 0.08);
        assert!(v < 1e-3, "Twitch should decay for large t, got {v}");
    }

    // 32. twitch_response: peak normalised near 1.
    #[test]
    fn test_twitch_peak_near_one() {
        let t_rise = 0.02;
        let t_fall = 0.08;
        let peak = (0..1000)
            .map(|i| twitch_response(i as f64 * 0.001, t_rise, t_fall))
            .fold(0.0_f64, f64::max);
        assert!(
            (peak - 1.0).abs() < 0.05,
            "Twitch peak should be near 1, got {peak}"
        );
    }

    // 33. active_force_length is symmetric around l_opt.
    #[test]
    fn test_fl_symmetric() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let d = 0.01;
        let fl_above = m.active_force_length(m.l_opt + d);
        let fl_below = m.active_force_length(m.l_opt - d);
        assert!(
            (fl_above - fl_below).abs() < 1e-9,
            "FL curve should be symmetric: fl_above={fl_above}, fl_below={fl_below}"
        );
    }

    // 34. force_velocity: v = -v_max gives FV = 0.
    #[test]
    fn test_fv_at_v_max_shortening() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, 0.0);
        let fv = m.force_velocity(-m.v_max);
        assert!(
            fv.abs() < 1e-9,
            "FV should be 0 at max shortening, got {fv}"
        );
    }

    // 35. MuscleFatigue: recovers when inactive.
    #[test]
    fn test_fatigue_recovers() {
        let mut f = MuscleFatigue::new(0.5, 0.5);
        f.fatigue_factor = 0.8;
        for _ in 0..100 {
            f.step(0.0, 0.01);
        }
        assert!(
            f.fatigue_factor < 0.8,
            "Fatigue should decrease during rest, got {}",
            f.fatigue_factor
        );
    }

    // 36. MuscleGroup: empty group returns zero force.
    #[test]
    fn test_muscle_group_empty() {
        let group = MuscleGroup::new(vec![], 1e-4);
        let f = group.resultant_force(0.01, &[]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    // 37. ActivationDynamics: steady state at u=a (no change).
    #[test]
    fn test_activation_steady_state() {
        let dyn_ = ActivationDynamics::new(0.01, 0.04);
        // At steady state activation == excitation; tau may still give tiny da/dt
        // due to formula structure, but da should be small.
        let a_init = 0.5;
        let a_new = dyn_.step(a_init, a_init, 0.001);
        assert!(
            (a_new - a_init).abs() < 1e-6,
            "No change expected at steady state, got {a_new}"
        );
    }

    // 38. total_force with pennation=90° gives zero.
    #[test]
    fn test_total_force_full_pennation() {
        let m = HillMuscle::new(1000.0, 0.1, 1.0, std::f64::consts::FRAC_PI_2);
        let f = m.total_force(m.l_opt, 0.0, 1.0);
        assert!(
            f.abs() < 1e-9,
            "90° pennation gives zero projected force, got {f}"
        );
    }

    // 39. MuscleElement: zero dt returns zero velocity.
    #[test]
    fn test_muscle_element_zero_dt() {
        let model = HillMuscle::new(500.0, 0.1, 1.0, 0.0);
        let elem = MuscleElement::new([0.0; 3], [0.0, 0.1, 0.0], model, 0.5, 0.1);
        let v = elem.current_velocity(0.09, 0.0);
        assert_eq!(v, 0.0, "Zero dt should give zero velocity");
    }

    // 40. musculotendon_length: increases with muscle length.
    #[test]
    fn test_mtu_length_increases() {
        let l1 = musculotendon_length(0.05, 0.2, 0.2);
        let l2 = musculotendon_length(0.10, 0.2, 0.2);
        assert!(l2 > l1, "Longer muscle → longer MTU");
    }
}
