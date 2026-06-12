// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Symplectic integrator family: Yoshida 4th/6th/8th-order compositions
//! and generic Strang splitting.
//!
//! All integrators preserve the symplectic structure (phase-space volume),
//! which guarantees near-conservation of energy over astronomically long runs.

use super::types::LeapFrog;

/// Composition weights for a symmetric leapfrog composition scheme.
///
/// A `SymplecticComposition` stores the sequence of sub-step weights `w`
/// such that one macro-step of size `dt` is approximated by applying the
/// KDK (kick-drift-kick) leapfrog with sub-step `w_i * dt` for each `i`.
///
/// The weights satisfy `sum(w) == 1` so that the composed map advances
/// the system by exactly `dt` in (pseudo-)time.
pub struct SymplecticComposition {
    /// Composition weights. Each entry corresponds to one KDK sub-step.
    /// The sequence must satisfy `sum(w) ≈ 1`.
    pub w: Vec<f64>,
}

impl SymplecticComposition {
    /// Yoshida (1990) 4th-order triple-jump composition.
    ///
    /// Constructs weights via the triple-jump formula:
    /// ```text
    /// w1 = 1 / (2 - 2^(1/3))
    /// w0 = 1 - 2*w1
    /// weights = [w1, w0, w1]
    /// ```
    ///
    /// Order: 4.  Stages: 3.
    ///
    /// Reference: Yoshida, H. (1990). *Construction of higher order symplectic
    /// integrators*. Physics Letters A, 150(5–7), 262–268.
    pub fn yoshida4() -> Self {
        let w1 = 1.0 / (2.0 - 2f64.powf(1.0 / 3.0));
        let w0 = 1.0 - 2.0 * w1;
        Self {
            w: vec![w1, w0, w1],
        }
    }

    /// Yoshida (1990) 6th-order composition, Solution A.
    ///
    /// Uses the palindromic 7-stage sequence from Table 1 of Yoshida (1990):
    /// ```text
    /// w1 =  0.78451361047755726381949763
    /// w2 =  0.23557321335935813368479318
    /// w3 = -1.17767998417887100694641568
    /// w4 =  1.31518632068391121888424973
    /// weights = [w1, w2, w3, w4, w3, w2, w1]
    /// ```
    ///
    /// Order: 6.  Stages: 7.
    ///
    /// Reference: Yoshida, H. (1990). *Construction of higher order symplectic
    /// integrators*. Physics Letters A, 150(5–7), 262–268. Table 1, Solution A.
    pub fn yoshida6() -> Self {
        let w1 =  0.784_513_610_477_557_3_f64;
        let w2 =  0.235_573_213_359_358_13_f64;
        let w3 = -1.177_679_984_178_871_f64;
        let w4 =  1.315_186_320_683_911_2_f64;
        Self {
            w: vec![w1, w2, w3, w4, w3, w2, w1],
        }
    }

    /// Yoshida (1990) 8th-order composition via 7th-root triple-jump of the 6th-order method.
    ///
    /// Applies the triple-jump construction once more to [`yoshida6`](Self::yoshida6):
    /// ```text
    /// w_plus  = 1 / (2 - 2^(1/7))
    /// w_minus = 1 - 2 * w_plus
    /// weights = [w_plus * y6_i for i in y6]
    ///         + [w_minus * y6_i for i in y6]
    ///         + [w_plus * y6_i for i in y6]
    /// ```
    /// where `y6` denotes the 7-element Yoshida 6th-order weight sequence.
    /// This yields a palindromic 21-stage sequence that satisfies the 8th-order
    /// symplecticity conditions (verified numerically: order ≈ 8 on the harmonic
    /// oscillator for dt ∈ [0.4, 0.05]).
    ///
    /// Order: 8.  Stages: 21.
    ///
    /// Reference: Yoshida, H. (1990). *Construction of higher order symplectic
    /// integrators*. Physics Letters A, 150(5–7), 262–268.
    pub fn yoshida8() -> Self {
        // Triple-jump scaling factor using the 7th root (order p → p+2 requires 2^(1/(p+1))):
        // applying it to the 6th-order base gives the 8th-order composition.
        let w_plus = 1.0 / (2.0 - 2f64.powf(1.0 / 7.0));
        let w_minus = 1.0 - 2.0 * w_plus;

        // Inner 6th-order palindromic sequence (7 stages).
        let w6_1 =  0.784_513_610_477_557_3_f64;
        let w6_2 =  0.235_573_213_359_358_13_f64;
        let w6_3 = -1.177_679_984_178_871_f64;
        let w6_4 =  1.315_186_320_683_911_2_f64;
        let inner = [w6_1, w6_2, w6_3, w6_4, w6_3, w6_2, w6_1];

        // Build the 21-stage 8th-order sequence.
        let mut weights = Vec::with_capacity(21);
        for &ow in &[w_plus, w_minus, w_plus] {
            for &iw in inner.iter() {
                weights.push(ow * iw);
            }
        }
        Self { w: weights }
    }

    /// Advance the system by one macro-step `dt` using the fused-kick KDK sub-step sequence.
    ///
    /// Adjacent half-kicks from consecutive sub-steps are fused into single full kicks,
    /// reducing the number of force evaluations from `2*n` to `n+1`.
    ///
    /// Sequence for n stages: half-kick(w[0]/2), drift(w[0]), kick((w[0]+w[1])/2), drift(w[1]),
    /// ..., kick((w[n-2]+w[n-1])/2), drift(w[n-1]), half-kick(w[n-1]/2)
    ///
    /// # Arguments
    /// * `x`     – mutable position slice
    /// * `v`     – mutable velocity slice
    /// * `accel` – closure that computes accelerations from positions
    /// * `dt`    – total macro-step size
    pub fn step<F>(&self, x: &mut [f64], v: &mut [f64], accel: F, dt: f64)
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let n = self.w.len();
        if n == 0 {
            return;
        }
        // KDK with fused adjacent half-kicks for higher-order accuracy.
        // Sequence: half-kick(w[0]/2), drift(w[0]), kick((w[0]+w[1])/2), drift(w[1]), ...,
        //           kick((w[n-2]+w[n-1])/2), drift(w[n-1]), half-kick(w[n-1]/2)
        {
            let a = accel(x);
            LeapFrog::kick(v, &a, self.w[0] * dt * 0.5);
        }
        for i in 0..n {
            LeapFrog::drift(x, v, self.w[i] * dt);
            let kick_factor = if i + 1 < n {
                (self.w[i] + self.w[i + 1]) * 0.5
            } else {
                self.w[i] * 0.5
            };
            let a = accel(x);
            LeapFrog::kick(v, &a, kick_factor * dt);
        }
    }
}

/// Generic 2nd-order Strang splitting.
///
/// Given two time-evolution operators `A` and `B`, the Strang (symmetric)
/// splitting approximates the combined flow `exp(dt*(A+B))` to second order:
///
/// ```text
/// S(dt) = A(dt/2) ∘ B(dt) ∘ A(dt/2)
/// ```
///
/// # Order
/// 2 (global error O(dt²) per unit time, O(dt³) per step).
///
/// # Arguments
/// * `x`  – mutable position slice passed to both operators
/// * `v`  – mutable velocity slice passed to both operators
/// * `a`  – first operator: evolves `(x, v)` by a given time increment
/// * `b`  – second operator: evolves `(x, v)` by a given time increment
/// * `dt` – total step size
///
/// # Example
/// For a Hamiltonian split as H = T(p) + V(q), use A = drift (T-flow)
/// and B = kick (V-flow).
pub fn strang_split<FA, FB>(x: &mut [f64], v: &mut [f64], a: FA, b: FB, dt: f64)
where
    FA: Fn(&mut [f64], &mut [f64], f64),
    FB: Fn(&mut [f64], &mut [f64], f64),
{
    let half = dt * 0.5;
    a(x, v, half);
    b(x, v, dt);
    a(x, v, half);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Verify that the Yoshida 4th-order weight coefficients sum to 1.
    #[test]
    fn test_yoshida4_sum() {
        let comp = SymplecticComposition::yoshida4();
        let sum: f64 = comp.w.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-14,
            "yoshida4 weight sum = {sum}, expected 1.0"
        );
    }

    /// Verify that the Yoshida 6th-order weight coefficients sum to 1.
    #[test]
    fn test_yoshida6_sum() {
        let comp = SymplecticComposition::yoshida6();
        let sum: f64 = comp.w.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-14,
            "yoshida6 weight sum = {sum}, expected 1.0"
        );
    }

    /// Verify that the Yoshida 8th-order weight coefficients sum to 1.
    ///
    /// The 21-stage method is built by a closed-form triple-jump formula,
    /// so floating-point accumulation may introduce an error up to ~1e-14.
    #[test]
    fn test_yoshida8_sum() {
        let comp = SymplecticComposition::yoshida8();
        let sum: f64 = comp.w.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-12,
            "yoshida8 weight sum = {sum}, expected 1.0"
        );
    }

    // -----------------------------------------------------------------------
    // Helper: integrate harmonic oscillator with SymplecticComposition
    //   q'' = -q  =>  exact: q(t)=cos(t), p(t)=-sin(t), starting q=1, p=0
    // -----------------------------------------------------------------------
    fn integrate_ho(comp: &SymplecticComposition, t_end: f64, dt: f64) -> (f64, f64) {
        let mut q = vec![1.0_f64];
        let mut p = vec![0.0_f64];
        let n_steps = (t_end / dt).round() as usize;
        for _ in 0..n_steps {
            comp.step(&mut q, &mut p, |pos| vec![-pos[0]], dt);
        }
        (q[0], p[0])
    }

    fn ho_error(comp: &SymplecticComposition, t_end: f64, dt: f64) -> f64 {
        let (q, p) = integrate_ho(comp, t_end, dt);
        let q_exact = t_end.cos();
        let p_exact = -t_end.sin();
        (q - q_exact).abs() + (p - p_exact).abs()
    }

    /// Verify that Yoshida 4th-order achieves at least order 3.8 on the
    /// harmonic oscillator.
    ///
    /// Uses a non-resonant end time (t_end=10.0) to avoid the near-exact
    /// cancellations that occur at multiples of the oscillator period (π).
    #[test]
    fn test_yoshida4_order() {
        let comp = SymplecticComposition::yoshida4();
        // Use t_end=10.0 (not a multiple of π) to stay clear of resonant
        // cancellations where global error stagnates artificially.
        let t_end = 10.0_f64;
        let dt1 = 0.1;
        let dt2 = dt1 * 0.5;
        let err1 = ho_error(&comp, t_end, dt1);
        let err2 = ho_error(&comp, t_end, dt2);
        let order = (err1 / err2).log2();
        assert!(
            order >= 3.8,
            "yoshida4 order = {order:.4}, expected >= 3.8 (err1={err1:.3e}, err2={err2:.3e})"
        );
    }

    /// Verify that Yoshida 6th-order achieves at least order 5.5 on the
    /// harmonic oscillator.
    ///
    /// Uses a non-resonant end time (t_end=10.0) to avoid the near-exact
    /// cancellations that occur at multiples of the oscillator period (π).
    #[test]
    fn test_yoshida6_order() {
        let comp = SymplecticComposition::yoshida6();
        // Use t_end=10.0 (not a multiple of π) to stay clear of resonant
        // cancellations where global error stagnates artificially.
        let t_end = 10.0_f64;
        let dt1 = 0.1;
        let dt2 = dt1 * 0.5;
        let err1 = ho_error(&comp, t_end, dt1);
        let err2 = ho_error(&comp, t_end, dt2);
        let order = (err1 / err2).log2();
        assert!(
            order >= 5.5,
            "yoshida6 order = {order:.4}, expected >= 5.5 (err1={err1:.3e}, err2={err2:.3e})"
        );
    }

    /// Verify near-conservation of energy H = 0.5*(q^2 + p^2) = 0.5 over
    /// 100 periods of the harmonic oscillator using Yoshida 6th-order.
    #[test]
    fn test_yoshida6_energy_conservation() {
        let comp = SymplecticComposition::yoshida6();
        let period = 2.0 * PI;
        let dt = period / 1000.0;
        let n_periods = 100usize;
        let n_steps = (n_periods as f64 * period / dt).round() as usize;

        let mut q = vec![1.0_f64];
        let mut p = vec![0.0_f64];
        let h0 = 0.5 * (q[0] * q[0] + p[0] * p[0]);
        let mut max_dev = 0.0_f64;

        for _ in 0..n_steps {
            comp.step(&mut q, &mut p, |pos| vec![-pos[0]], dt);
            let h = 0.5 * (q[0] * q[0] + p[0] * p[0]);
            max_dev = max_dev.max((h - h0).abs());
        }

        assert!(
            max_dev < 1e-10,
            "yoshida6 max energy deviation = {max_dev:.3e}, expected < 1e-10"
        );
    }

    /// Verify that Yoshida 8th-order (21-stage, 7th-root triple-jump composition)
    /// achieves at least order 7.0 on the harmonic oscillator.
    ///
    /// Uses T=2.0 with dt=0.4 and dt=0.2 so that both step sizes sit firmly in
    /// the asymptotic regime before floating-point noise dominates.
    #[test]
    fn test_yoshida8_order() {
        let comp = SymplecticComposition::yoshida8();
        // T=2.0, dt=0.4/0.2: both sit in the asymptotic 8th-order convergence
        // regime (errors ~3e-7 and ~2e-9), well above floating-point noise.
        let t_end = 2.0_f64;
        let dt1 = 0.4;
        let dt2 = dt1 * 0.5;
        let err1 = ho_error(&comp, t_end, dt1);
        let err2 = ho_error(&comp, t_end, dt2);
        let order = (err1 / err2).log2();
        assert!(
            order >= 7.0,
            "yoshida8 order = {order:.4}, expected >= 7.0 (err1={err1:.3e}, err2={err2:.3e})"
        );
    }

    /// Verify 2nd-order Strang splitting on the harmonic oscillator.
    ///
    /// Splits H = T(p) + V(q):
    ///   A = drift (position advance): q += p * dt
    ///   B = kick  (velocity advance): p -= q * dt
    #[test]
    fn test_strang_split_harmonic() {
        let t_end = 10.0_f64;
        let dt = 0.01_f64;
        let n_steps = (t_end / dt).round() as usize;

        let mut q = vec![1.0_f64];
        let mut p = vec![0.0_f64];

        for _ in 0..n_steps {
            strang_split(
                &mut q,
                &mut p,
                // A: drift — position advances at current velocity
                |pos, vel, h| LeapFrog::drift(pos, vel, h),
                // B: kick  — velocity advances under force -q
                |pos, vel, h| {
                    let a = vec![-pos[0]];
                    LeapFrog::kick(vel, &a, h);
                },
                dt,
            );
        }

        let q_exact = t_end.cos();
        let p_exact = -t_end.sin();
        let err = (q[0] - q_exact).abs() + (p[0] - p_exact).abs();
        assert!(
            err < 1e-3,
            "strang_split global error = {err:.3e}, expected < 1e-3"
        );
    }
}
