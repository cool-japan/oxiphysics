// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Response spectra of a ground-motion record.
//!
//! The relative motion `x(t)` of a single-degree-of-freedom oscillator of natural
//! frequency `ω = 2π/T` and damping ratio `ζ` under the ground acceleration
//! `a_g(t)` obeys
//!
//! ```text
//! ẍ + 2ζω ẋ + ω² x = −a_g(t).
//! ```
//!
//! The record is taken to be piecewise linear between its samples.  For such a
//! load the state `z = [x, ẋ]` at the end of a step of length `h` is an *exact*
//! linear function of the state at its start and the two bracketing load values
//! (Nigam & Jennings, 1969):
//!
//! ```text
//! z_{i+1} = E z_i + f₀ p_i + f₁ p_{i+1},        p = −a_g.
//! ```
//!
//! `E = e^{Ah}` and the load columns `f₀`, `f₁` follow from one matrix
//! exponential of the system extended by the load and its slope (valid for under-,
//! critically and over-damped oscillators alike, `ω = 0` included).  They depend
//! only on `(ω, ζ, h)`, so they are formed once per period and reused while the
//! step does not change.  There is no time-stepping error: the result is exact
//! for the piecewise-linear record, however coarse the steps are.

use rayon::prelude::*;
use std::f64::consts::PI;

/// Ground acceleration record for response spectrum analysis.
#[derive(Debug, Clone)]
pub struct GroundMotion {
    /// Time values (seconds).
    pub time: Vec<f64>,
    /// Ground acceleration values (m/s²).
    pub acceleration: Vec<f64>,
}

/// Response spectrum data: peak response vs. natural period.
///
/// For every period `T` (`ω = 2π/T`) the three spectra are taken from the peak
/// response of the damped oscillator to the record:
///
/// * `spectral_displacement` — peak relative displacement `Sd = max|x|`;
/// * `spectral_velocity` — peak **relative** velocity `Sv = max|ẋ|` (not the
///   pseudo-velocity `ω·Sd`);
/// * `spectral_acceleration` — **pseudo**-spectral acceleration `Sa = ω²·Sd`
///   (not the peak absolute acceleration `max|ẍ + a_g|`; the two coincide for
///   `ζ = 0`).
#[derive(Debug, Clone)]
pub struct ResponseSpectrum {
    /// Natural periods T (seconds).
    pub periods: Vec<f64>,
    /// Spectral acceleration Sa(T) (m/s²), pseudo-acceleration `ω² Sd`.
    pub spectral_acceleration: Vec<f64>,
    /// Spectral velocity Sv(T) (m/s), peak relative velocity.
    pub spectral_velocity: Vec<f64>,
    /// Spectral displacement Sd(T) (m).
    pub spectral_displacement: Vec<f64>,
}

/// Relative tolerance below which a step is considered equal to the cached one.
const STEP_REL_TOL: f64 = 1.0e-9;

type Mat4 = [[f64; 4]; 4];

fn mat4_mul(a: &Mat4, b: &Mat4) -> Mat4 {
    let mut c = [[0.0_f64; 4]; 4];
    for (i, row) in c.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

/// Matrix exponential of a 4 x 4 matrix: Taylor series on the matrix scaled to
/// 1-norm at most `1/2` followed by repeated squaring.
fn mat4_exp(m: &Mat4) -> Mat4 {
    let norm = m
        .iter()
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    if !norm.is_finite() {
        return [[f64::NAN; 4]; 4];
    }
    let mut halvings = 0u32;
    let mut scaled_norm = norm;
    while scaled_norm > 0.5 {
        scaled_norm *= 0.5;
        halvings += 1;
    }
    let factor = 0.5_f64.powi(halvings as i32);
    let mut a = *m;
    for row in a.iter_mut() {
        for v in row.iter_mut() {
            *v *= factor;
        }
    }
    let mut result = [[0.0_f64; 4]; 4];
    let mut term = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        result[i][i] = 1.0;
        term[i][i] = 1.0;
    }
    for k in 1..=18 {
        term = mat4_mul(&term, &a);
        let inv_k = 1.0 / f64::from(k);
        for (result_row, term_row) in result.iter_mut().zip(term.iter_mut()) {
            for (r, t) in result_row.iter_mut().zip(term_row.iter_mut()) {
                *t *= inv_k;
                *r += *t;
            }
        }
    }
    for _ in 0..halvings {
        result = mat4_mul(&result, &result);
    }
    result
}

/// Exact one-step map of the oscillator for a load that is linear within the step.
#[derive(Debug, Clone, Copy)]
struct StepOperator {
    /// Step length `h` the operator was built for.
    step: f64,
    /// State transition `E = e^{A h}` (row-major `[[x←x, x←ẋ], [ẋ←x, ẋ←ẋ]]`).
    transition: [[f64; 2]; 2],
    /// Response to the load value at the start of the step.
    load_start: [f64; 2],
    /// Response to the load value at the end of the step.
    load_end: [f64; 2],
}

impl StepOperator {
    /// Build the operator for the oscillator `(ω, ζ)` and the step `h > 0`.
    ///
    /// The state is extended by the load `p` and its constant slope `q`:
    /// `x' = v`, `v' = −ω²x − 2ζωv + p`, `p' = q`, `q' = 0`, whose exponential
    /// contains `E` and the integrals of the load response.  The variables are
    /// scaled by `s = max(ω, 1/h)` (`v → v/s`, `p → p/s²`, `q → q/s³`) so that the
    /// exponent matrix has entries of order `ωh` and the scaling-and-squaring stays
    /// accurate for very stiff oscillators and for `ω → 0`.
    fn new(omega: f64, zeta: f64, step: f64) -> Self {
        let s = omega.max(1.0 / step);
        let sh = s * step;
        let exponent: Mat4 = [
            [0.0, sh, 0.0, 0.0],
            [
                -omega * omega * step / s,
                -2.0 * zeta * omega * step,
                sh,
                0.0,
            ],
            [0.0, 0.0, 0.0, sh],
            [0.0, 0.0, 0.0, 0.0],
        ];
        let e = mat4_exp(&exponent);
        let transition = [[e[0][0], e[0][1] / s], [e[1][0] * s, e[1][1]]];
        // Response to a load p(τ) = p_start + q τ:  G_p p_start + G_q q.
        let g_p = [e[0][2] / (s * s), e[1][2] / s];
        let g_q = [e[0][3] / (s * s * s), e[1][3] / (s * s)];
        // q = (p_end − p_start)/h
        let load_start = [g_p[0] - g_q[0] / step, g_p[1] - g_q[1] / step];
        let load_end = [g_q[0] / step, g_q[1] / step];
        Self {
            step,
            transition,
            load_start,
            load_end,
        }
    }

    /// Advance `[x, ẋ]` over one step with load `p_start → p_end`.
    fn advance(&self, state: [f64; 2], p_start: f64, p_end: f64) -> [f64; 2] {
        let e = &self.transition;
        [
            e[0][0] * state[0]
                + e[0][1] * state[1]
                + self.load_start[0] * p_start
                + self.load_end[0] * p_end,
            e[1][0] * state[0]
                + e[1][1] * state[1]
                + self.load_start[1] * p_start
                + self.load_end[1] * p_end,
        ]
    }
}

/// A ground-motion record reduced to strictly increasing, finite samples.
struct Record {
    knots: Vec<(f64, f64)>,
    peak_ground_acceleration: f64,
}

impl Record {
    /// Pair `time` with `acceleration` (the shorter of the two decides the length)
    /// and drop samples that are not finite or whose time does not advance.
    fn new(ground_motion: &GroundMotion) -> Self {
        let mut knots: Vec<(f64, f64)> = Vec::with_capacity(
            ground_motion
                .time
                .len()
                .min(ground_motion.acceleration.len()),
        );
        for (&t, &a) in ground_motion
            .time
            .iter()
            .zip(ground_motion.acceleration.iter())
        {
            if !t.is_finite() || !a.is_finite() {
                continue;
            }
            if let Some(&(t_prev, _)) = knots.last()
                && t <= t_prev
            {
                continue;
            }
            knots.push((t, a));
        }
        let peak_ground_acceleration = knots.iter().map(|&(_, a)| a.abs()).fold(0.0_f64, f64::max);
        Self {
            knots,
            peak_ground_acceleration,
        }
    }

    /// Integrate the oscillator `(ω, ζ)` from rest over the record, calling
    /// `visit(sample_index, x, ẋ)` at every sample (the initial rest state is
    /// reported at index 0).
    fn integrate(&self, omega: f64, zeta: f64, mut visit: impl FnMut(usize, f64, f64)) {
        if self.knots.is_empty() {
            return;
        }
        let mut state = [0.0_f64, 0.0_f64];
        visit(0, state[0], state[1]);
        let mut cached: Option<StepOperator> = None;
        for (i, pair) in self.knots.windows(2).enumerate() {
            let (t0, a0) = pair[0];
            let (t1, a1) = pair[1];
            let step = t1 - t0;
            let operator = match cached {
                Some(op) if (step - op.step).abs() <= STEP_REL_TOL * op.step => op,
                _ => {
                    let op = StepOperator::new(omega, zeta, step);
                    cached = Some(op);
                    op
                }
            };
            state = operator.advance(state, -a0, -a1);
            visit(i + 1, state[0], state[1]);
        }
    }

    /// `(Sd, Sv, Sa)` of the oscillator with period `period` and damping ratio `zeta`.
    ///
    /// A period that is not positive (or NaN) is the rigid-oscillator limit
    /// `T → 0`: no relative motion, `Sd = Sv = 0` and `Sa` equal to the peak
    /// ground acceleration.
    fn spectral_point(&self, period: f64, zeta: f64) -> (f64, f64, f64) {
        if period.is_nan() || period <= 0.0 {
            return (0.0, 0.0, self.peak_ground_acceleration);
        }
        let omega = 2.0 * PI / period;
        let mut peak_x = 0.0_f64;
        let mut peak_v = 0.0_f64;
        let mut finite = true;
        self.integrate(omega, zeta, |_, x, v| {
            finite &= x.is_finite() && v.is_finite();
            peak_x = peak_x.max(x.abs());
            peak_v = peak_v.max(v.abs());
        });
        if !finite {
            return (f64::NAN, f64::NAN, f64::NAN);
        }
        (peak_x, peak_v, omega * omega * peak_x)
    }
}

/// Compute a response spectrum from a ground motion record.
///
/// For each natural period T in `periods`:
/// 1. Compute ω = 2π/T
/// 2. Integrate the SDOF equation ẍ + 2ζω ẋ + ω²x = −ẍ_g(t) from rest with the
///    exact piecewise-linear recurrence of Nigam & Jennings (see the module
///    documentation); the record is interpolated linearly between samples
/// 3. Record the peak absolute values of x (Sd), ẋ (Sv) and ω²x (Sa, the
///    pseudo-acceleration) over the record
///
/// `damping_ratio` may be any value `ζ ≥ 0` (under-, critically and over-damped
/// oscillators are all exact); a negative or NaN ratio is treated as `0`.  Samples
/// that are not finite or whose time does not advance are ignored, and a record
/// whose `time` and `acceleration` differ in length is truncated to the shorter.
/// The periods are processed in parallel; the result is deterministic.  A period
/// `T ≤ 0` gives the rigid limit (`Sd = Sv = 0`, `Sa` = peak ground acceleration).
///
/// The peaks are taken at the sample instants, so the record should resolve the
/// shortest period (at least about ten samples per period).
pub fn compute_response_spectrum(
    ground_motion: &GroundMotion,
    periods: &[f64],
    damping_ratio: f64,
) -> ResponseSpectrum {
    let zeta = damping_ratio.max(0.0);
    let record = Record::new(ground_motion);
    let points: Vec<(f64, f64, f64)> = periods
        .par_iter()
        .map(|&period| record.spectral_point(period, zeta))
        .collect();
    let mut sd = Vec::with_capacity(points.len());
    let mut sv = Vec::with_capacity(points.len());
    let mut sa = Vec::with_capacity(points.len());
    for (peak_x, peak_v, pseudo_acceleration) in points {
        sd.push(peak_x);
        sv.push(peak_v);
        sa.push(pseudo_acceleration);
    }
    ResponseSpectrum {
        periods: periods.to_vec(),
        spectral_acceleration: sa,
        spectral_velocity: sv,
        spectral_displacement: sd,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Closed-form relative `(x, ẋ)` of an oscillator started at rest under the
    /// constant ground acceleration `a0`, valid for every damping regime.
    fn step_response(a0: f64, omega: f64, zeta: f64, t: f64) -> (f64, f64) {
        let decay = (-zeta * omega * t).exp();
        let q = 1.0 - zeta * zeta;
        if q > 0.0 {
            let wd = omega * q.sqrt();
            let g = decay * ((wd * t).cos() + zeta * omega / wd * (wd * t).sin());
            (
                -(a0 / (omega * omega)) * (1.0 - g),
                -(a0 / wd) * decay * (wd * t).sin(),
            )
        } else if q == 0.0 {
            let g = decay * (1.0 + omega * t);
            (-(a0 / (omega * omega)) * (1.0 - g), -a0 * t * decay)
        } else {
            let ws = omega * (-q).sqrt();
            let g = decay * ((ws * t).cosh() + zeta * omega / ws * (ws * t).sinh());
            (
                -(a0 / (omega * omega)) * (1.0 - g),
                -(a0 / ws) * decay * (ws * t).sinh(),
            )
        }
    }

    /// Closed-form `(x, ẋ)` under the ramp `a_g = c t` from rest (`ζ < 1`).
    fn ramp_response(c: f64, omega: f64, zeta: f64, t: f64) -> (f64, f64) {
        let wd = omega * (1.0 - zeta * zeta).sqrt();
        let decay = (-zeta * omega * t).exp();
        let (s, co) = ((wd * t).sin(), (wd * t).cos());
        let k = c / (omega * omega);
        let x = -k
            * (t - 2.0 * zeta / omega
                + decay * (2.0 * zeta / omega * co + (2.0 * zeta * zeta - 1.0) / wd * s));
        // d/dt of the expression above
        let dx = -k
            * (1.0
                + decay
                    * (-zeta
                        * omega
                        * (2.0 * zeta / omega * co + (2.0 * zeta * zeta - 1.0) / wd * s)
                        - 2.0 * zeta / omega * wd * s
                        + (2.0 * zeta * zeta - 1.0) * co));
        (x, dx)
    }

    fn history(record: &Record, omega: f64, zeta: f64) -> Vec<(f64, f64)> {
        let mut out = vec![(0.0, 0.0); record.knots.len()];
        record.integrate(omega, zeta, |i, x, v| out[i] = (x, v));
        out
    }

    fn constant_record(a0: f64, dt: f64, n: usize) -> GroundMotion {
        GroundMotion {
            time: (0..n).map(|i| i as f64 * dt).collect(),
            acceleration: vec![a0; n],
        }
    }

    #[test]
    fn constant_acceleration_step_matches_the_closed_form_for_every_damping_regime() {
        let (a0, omega) = (2.5, 2.0 * PI / 0.8);
        let gm = constant_record(a0, 0.01, 401);
        let record = Record::new(&gm);
        for &zeta in &[0.0, 0.05, 0.3, 0.999, 1.0, 1.001, 2.5] {
            let h = history(&record, omega, zeta);
            let scale = a0 / (omega * omega);
            for (i, &(x, v)) in h.iter().enumerate() {
                let t = record.knots[i].0;
                let (xe, ve) = step_response(a0, omega, zeta, t);
                assert!(
                    (x - xe).abs() <= 1e-11 * scale.max(xe.abs()),
                    "zeta {zeta} t {t}: x {x} vs {xe}"
                );
                assert!(
                    (v - ve).abs() <= 1e-11 * (a0 / omega).max(ve.abs()),
                    "zeta {zeta} t {t}: v {v} vs {ve}"
                );
            }
        }
    }

    #[test]
    fn ramp_with_uneven_steps_matches_the_closed_form() {
        // A linear ramp is represented exactly by any set of knots; vary the step so the
        // cached operator has to be rebuilt.
        let (c, omega) = (0.7, 5.0);
        let steps = [0.01, 0.01, 0.02, 0.005, 0.005, 0.04, 0.01, 0.0100001, 0.03];
        let mut time = vec![0.0];
        for k in 0..45 {
            let last = *time.last().unwrap_or(&0.0);
            time.push(last + steps[k % steps.len()]);
        }
        let acceleration: Vec<f64> = time.iter().map(|&t| c * t).collect();
        let record = Record::new(&GroundMotion { time, acceleration });
        for &zeta in &[0.0, 0.05, 0.5] {
            let h = history(&record, omega, zeta);
            for (i, &(x, v)) in h.iter().enumerate() {
                let t = record.knots[i].0;
                let (xe, ve) = ramp_response(c, omega, zeta, t);
                assert!(
                    (x - xe).abs() <= 1e-10 * (1.0 + xe.abs()),
                    "zeta {zeta} t {t}: x {x} vs {xe}"
                );
                assert!(
                    (v - ve).abs() <= 1e-10 * (1.0 + ve.abs()),
                    "zeta {zeta} t {t}: v {v} vs {ve}"
                );
            }
        }
    }

    #[test]
    fn refining_the_record_by_linear_midpoints_does_not_change_the_response() {
        // Inserting the linear midpoint of every interval leaves the piecewise-linear load
        // untouched, so an exact integrator must reproduce the states at the old samples.
        // A time-stepping scheme (Newmark, Euler) changes with the step and fails this.
        let dt = 0.02;
        let n = 301;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time
            .iter()
            .map(|&t| (3.0 * t).sin() + 0.5 * (11.0 * t).cos() - 0.2)
            .collect();
        let coarse = Record::new(&GroundMotion {
            time: time.clone(),
            acceleration: acceleration.clone(),
        });
        let mut fine_time = Vec::new();
        let mut fine_acc = Vec::new();
        for i in 0..n - 1 {
            fine_time.push(time[i]);
            fine_acc.push(acceleration[i]);
            fine_time.push(0.5 * (time[i] + time[i + 1]));
            fine_acc.push(0.5 * (acceleration[i] + acceleration[i + 1]));
        }
        fine_time.push(time[n - 1]);
        fine_acc.push(acceleration[n - 1]);
        let fine = Record::new(&GroundMotion {
            time: fine_time,
            acceleration: fine_acc,
        });
        for &(period, zeta) in &[(0.15, 0.02), (0.6, 0.05), (2.0, 0.0), (0.05, 0.4)] {
            let omega = 2.0 * PI / period;
            let hc = history(&coarse, omega, zeta);
            let hf = history(&fine, omega, zeta);
            let peak = hc.iter().fold(0.0_f64, |m, s| m.max(s.0.abs())).max(1e-300);
            let peak_v = hc.iter().fold(0.0_f64, |m, s| m.max(s.1.abs())).max(1e-300);
            for (i, &(x, v)) in hc.iter().enumerate() {
                let (xf, vf) = hf[2 * i];
                assert!(
                    (x - xf).abs() <= 1e-11 * peak,
                    "T {period} sample {i}: {x} vs {xf}"
                );
                assert!(
                    (v - vf).abs() <= 1e-11 * peak_v,
                    "T {period} sample {i}: {v} vs {vf}"
                );
            }
        }
    }

    #[test]
    fn harmonic_ground_motion_at_resonance_reaches_the_steady_state_amplitude() {
        // ẍ + 2ζω ẋ + ω² x = −sin(ω t) from rest has the closed-form solution
        //   x(t) = A [cos ωt − e^{−ζωt}(cos ω_d t + ζω/ω_d sin ω_d t)],  A = 1/(2 ζ ω²).
        let (zeta, period): (f64, f64) = (0.05, 1.0);
        let omega = 2.0 * PI / period;
        let omega_d = omega * (1.0 - zeta * zeta).sqrt();
        let amplitude = 1.0 / (2.0 * zeta * omega * omega);
        let dt = 0.005; // 200 samples per period
        let n = 8001; // 40 s >> 1/(ζ ω) = 3.2 s: the transient has died out
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time.iter().map(|&t| (omega * t).sin()).collect();
        let gm = GroundMotion { time, acceleration };

        // The whole time history follows the closed form to the (ω dt)²/12 error of
        // interpolating the sine linearly between samples.
        let record = Record::new(&gm);
        let states = history(&record, omega, zeta);
        for (i, &(x, _)) in states.iter().enumerate() {
            let t = record.knots[i].0;
            let exact = amplitude
                * ((omega * t).cos()
                    - (-zeta * omega * t).exp()
                        * ((omega_d * t).cos() + zeta * omega / omega_d * (omega_d * t).sin()));
            assert!(
                (x - exact).abs() < 5.0e-4 * amplitude,
                "t = {t}: x = {x}, closed form {exact}"
            );
        }

        let spectrum = compute_response_spectrum(&gm, &[period], zeta);
        let rel = (spectrum.spectral_displacement[0] - amplitude).abs() / amplitude;
        assert!(
            rel < 5.0e-4,
            "Sd {} vs {amplitude}",
            spectrum.spectral_displacement[0]
        );
        // Pseudo-acceleration amplification at resonance is 1/(2ζ) = 10.
        let rel_a = (spectrum.spectral_acceleration[0] - 1.0 / (2.0 * zeta)).abs() * 2.0 * zeta;
        assert!(rel_a < 5.0e-4, "Sa {}", spectrum.spectral_acceleration[0]);
        // At resonance the peak velocity is ω Sd to within the damping correction.
        let ratio = spectrum.spectral_velocity[0] / (omega * spectrum.spectral_displacement[0]);
        assert!((ratio - 1.0).abs() < 1.0e-3, "Sv/(ω Sd) = {ratio}");
    }

    #[test]
    fn rigid_oscillator_limit_has_zero_displacement_and_peak_ground_acceleration() {
        let dt = 0.005;
        let time: Vec<f64> = (0..2001).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time
            .iter()
            .map(|&t| 3.0 * (2.0 * PI * 1.5 * t).sin() * (-0.3 * t).exp())
            .collect();
        let pga = acceleration.iter().fold(0.0_f64, |m, a| m.max(a.abs()));
        let gm = GroundMotion { time, acceleration };
        let spectrum = compute_response_spectrum(&gm, &[0.0, 1.0e-4, 1.0e-6, -1.0, f64::NAN], 0.05);
        // Explicit T = 0 (and T < 0, NaN): rigid, Sa = PGA.
        for &idx in &[0usize, 3, 4] {
            assert_eq!(spectrum.spectral_displacement[idx], 0.0);
            assert_eq!(spectrum.spectral_velocity[idx], 0.0);
            assert_eq!(spectrum.spectral_acceleration[idx], pga);
        }
        // T → 0: Sd → 0 like 1/ω², Sa → PGA.
        for &idx in &[1usize, 2] {
            let omega = 2.0 * PI / spectrum.periods[idx];
            assert!(spectrum.spectral_displacement[idx] <= 1.0001 * pga / (omega * omega));
            assert!(spectrum.spectral_displacement[idx] < 1.0e-8);
            let rel = (spectrum.spectral_acceleration[idx] - pga).abs() / pga;
            assert!(
                rel < 1.0e-3,
                "T {}: Sa {} vs PGA {pga}",
                spectrum.periods[idx],
                spectrum.spectral_acceleration[idx]
            );
        }
    }

    #[test]
    fn displacement_spectrum_is_continuous_across_critical_damping() {
        // Short transient pulse; Sd(ζ) must vary smoothly through ζ = 1.
        let dt = 0.01;
        let n = 600;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time
            .iter()
            .map(|&t| if t < 0.5 { (PI * t / 0.5).sin() } else { 0.0 })
            .collect();
        let gm = GroundMotion { time, acceleration };
        let periods = [0.3, 1.0, 2.0];
        let sd = |zeta: f64| compute_response_spectrum(&gm, &periods, zeta).spectral_displacement;
        let below = sd(0.999);
        let critical = sd(1.0);
        let above = sd(1.001);
        for i in 0..periods.len() {
            assert!(critical[i] > 0.0);
            assert!(
                (below[i] - critical[i]).abs() < 2e-3 * critical[i],
                "T {}: {} vs {}",
                periods[i],
                below[i],
                critical[i]
            );
            assert!(
                (above[i] - critical[i]).abs() < 2e-3 * critical[i],
                "T {}: {} vs {}",
                periods[i],
                above[i],
                critical[i]
            );
        }
        // More damping, less displacement.
        let heavy = sd(3.0);
        for i in 0..periods.len() {
            assert!(heavy[i] < critical[i]);
        }
    }

    #[test]
    fn zero_frequency_oscillator_integrates_the_ground_motion_twice() {
        // ω = 0, ζ = 0: ẍ = −a_g.  For constant a_g = a: x = −a t²/2, ẋ = −a t.
        let gm = constant_record(1.5, 0.25, 9);
        let record = Record::new(&gm);
        let h = history(&record, 0.0, 0.0);
        for (i, &(x, v)) in h.iter().enumerate() {
            let t = record.knots[i].0;
            assert!((x + 0.75 * t * t).abs() < 1e-13, "x({t}) = {x}");
            assert!((v + 1.5 * t).abs() < 1e-13, "v({t}) = {v}");
        }
        // Infinite period is the same oscillator.
        let spectrum = compute_response_spectrum(&gm, &[f64::INFINITY], 0.1);
        assert!(spectrum.spectral_displacement[0].is_finite());
    }

    #[test]
    fn spectrum_points_are_independent_of_the_parallel_schedule() {
        let dt = 0.01;
        let time: Vec<f64> = (0..400).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time.iter().map(|&t| (7.0 * t).sin() * (1.0 + t)).collect();
        let gm = GroundMotion { time, acceleration };
        let periods: Vec<f64> = (1..=40).map(|i| 0.05 * i as f64).collect();
        let spectrum = compute_response_spectrum(&gm, &periods, 0.05);
        let record = Record::new(&gm);
        for (i, &t) in periods.iter().enumerate() {
            let (sd, sv, sa) = record.spectral_point(t, 0.05);
            assert_eq!(spectrum.spectral_displacement[i], sd);
            assert_eq!(spectrum.spectral_velocity[i], sv);
            assert_eq!(spectrum.spectral_acceleration[i], sa);
        }
        assert_eq!(spectrum.periods, periods);
    }

    #[test]
    fn malformed_records_are_handled_without_panicking() {
        // acceleration shorter than time: truncated, no out-of-bounds access.
        let gm = GroundMotion {
            time: vec![0.0, 0.1, 0.2, 0.3, 0.4],
            acceleration: vec![1.0, 1.0, 1.0],
        };
        let spectrum = compute_response_spectrum(&gm, &[0.5], 0.05);
        let reference = compute_response_spectrum(
            &GroundMotion {
                time: vec![0.0, 0.1, 0.2],
                acceleration: vec![1.0, 1.0, 1.0],
            },
            &[0.5],
            0.05,
        );
        assert_eq!(
            spectrum.spectral_displacement,
            reference.spectral_displacement
        );
        // Duplicate and decreasing times and non-finite samples are skipped.
        let messy = GroundMotion {
            time: vec![0.0, 0.1, 0.1, 0.05, f64::NAN, 0.2],
            acceleration: vec![1.0, 1.0, 7.0, 7.0, 7.0, 1.0],
        };
        let messy_spectrum = compute_response_spectrum(&messy, &[0.5], 0.05);
        assert_eq!(
            messy_spectrum.spectral_displacement,
            reference.spectral_displacement
        );
        // Empty and single-sample records produce zero response.
        let empty = compute_response_spectrum(
            &GroundMotion {
                time: Vec::new(),
                acceleration: Vec::new(),
            },
            &[0.5, 0.0],
            0.05,
        );
        assert_eq!(empty.spectral_displacement, vec![0.0, 0.0]);
        assert_eq!(empty.spectral_acceleration, vec![0.0, 0.0]);
        let single = compute_response_spectrum(
            &GroundMotion {
                time: vec![0.0],
                acceleration: vec![2.0],
            },
            &[0.5, 0.0],
            0.05,
        );
        assert_eq!(single.spectral_displacement, vec![0.0, 0.0]);
        assert_eq!(single.spectral_acceleration, vec![0.0, 2.0]);
        // Negative or NaN damping behaves as zero damping.
        let undamped = compute_response_spectrum(&gm, &[0.5], 0.0);
        let negative = compute_response_spectrum(&gm, &[0.5], -0.3);
        let nan = compute_response_spectrum(&gm, &[0.5], f64::NAN);
        assert_eq!(
            undamped.spectral_displacement,
            negative.spectral_displacement
        );
        assert_eq!(undamped.spectral_displacement, nan.spectral_displacement);
    }

    #[test]
    fn step_operator_reproduces_the_constant_load_response_over_one_step() {
        // Single step from rest under the constant load p: closed form of an undamped
        // oscillator, x = p(1 − cos ωh)/ω², v = p sin(ωh)/ω.
        let (omega, h, p) = (3.0, 0.4, -2.0);
        let op = StepOperator::new(omega, 0.0, h);
        let z = op.advance([0.0, 0.0], p, p);
        assert!((z[0] - p * (1.0 - (omega * h).cos()) / (omega * omega)).abs() < 1e-14);
        assert!((z[1] - p * (omega * h).sin() / omega).abs() < 1e-14);
        // Free vibration: E rotates (x, v/ω) by ωh.
        let free = op.advance([1.0, 0.0], 0.0, 0.0);
        assert!((free[0] - (omega * h).cos()).abs() < 1e-14);
        assert!((free[1] + omega * (omega * h).sin()).abs() < 1e-13);
    }

    #[test]
    fn matrix_exponential_handles_zero_nilpotent_and_rotation_matrices() {
        let zero = [[0.0; 4]; 4];
        let e = mat4_exp(&zero);
        for (i, row) in e.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert_eq!(v, if i == j { 1.0 } else { 0.0 });
            }
        }
        // Nilpotent shift: exp(N) = I + N + N²/2 + N³/6.
        let mut n = [[0.0; 4]; 4];
        n[0][1] = 2.0;
        n[1][2] = 3.0;
        n[2][3] = 4.0;
        let e = mat4_exp(&n);
        assert!((e[0][1] - 2.0).abs() < 1e-14);
        assert!((e[0][2] - 3.0).abs() < 1e-14); // 2·3/2
        assert!((e[0][3] - 4.0).abs() < 1e-13); // 2·3·4/6
        // A large rotation exercises the scaling and squaring.
        let theta = 250.0_f64;
        let mut r = [[0.0; 4]; 4];
        r[0][1] = theta;
        r[1][0] = -theta;
        let e = mat4_exp(&r);
        assert!((e[0][0] - theta.cos()).abs() < 1e-11);
        assert!((e[0][1] - theta.sin()).abs() < 1e-11);
        assert!(mat4_exp(&[[f64::INFINITY; 4]; 4])[0][0].is_nan());
    }

    /// Natural periods (s) of the reference spectra.
    const REFERENCE_PERIODS: [f64; 5] = [0.05, 0.2, 0.77, 1.5, 3.0];
    /// Peak responses `(Sd, Sv, Sa)` of one oscillator.
    type Peaks = (f64, f64, f64);
    /// A damping ratio with the peaks at every reference period.
    type SpectrumRow = (f64, [Peaks; 5]);
    /// Reference spectra of the broadband record of
    /// `broadband_record_matches_an_independent_ode_reference`.  The peaks at the 501
    /// sample instants come from integrating each linear segment of the record with an
    /// adaptive 8th-order Dormand-Prince solver (rtol 1e-13), a computation that
    /// shares nothing with the matrix-exponential recurrence.
    const REFERENCE_SPECTRA: [SpectrumRow; 4] = [
        (
            0.0,
            [
                (1.95225928395e-4, 7.34698233228e-3, 3.08288429134),
                (5.89841250614e-3, 0.127370273481, 5.82149980300),
                (0.264756033888, 2.09419691549, 17.6288569221),
                (0.0932648867596, 0.539784744433, 1.63642228770),
                (0.140826389514, 0.414860931095, 0.617733668325),
            ],
        ),
        (
            0.05,
            [
                (1.66123069538e-4, 5.12531912368e-3, 2.62331036518),
                (4.54847301963e-3, 0.085024942035, 4.48916293327),
                (0.109973091441, 0.934978212364, 7.32259003063),
                (0.0858061809175, 0.503614827403, 1.50555210813),
                (0.133909794358, 0.402887557468, 0.587394086998),
            ],
        ),
        (
            1.0,
            [
                (1.52801563219e-4, 2.10209683266e-3, 2.41294556934),
                (2.12013424687e-3, 0.0210290532151, 2.09248862938),
                (0.0169983973896, 0.120884489358, 1.13184319574),
                (0.0357350343946, 0.208300754575, 0.627005604861),
                (0.0665086585212, 0.259595597594, 0.291739621712),
            ],
        ),
        (
            1.7,
            [
                (1.46577345856e-4, 1.64875511190e-3, 2.31465666841),
                (1.74749467424e-3, 0.0133867055881, 1.72470811277),
                (0.0118200107075, 0.0789294505102, 0.787038824085),
                (0.0246736091423, 0.138653557559, 0.432922242456),
                (0.0472048174947, 0.204184718836, 0.207063499777),
            ],
        ),
    ];

    #[test]
    fn broadband_record_matches_an_independent_ode_reference() {
        let dt = 0.01;
        let time: Vec<f64> = (0..501).map(|i| i as f64 * dt).collect();
        let acceleration: Vec<f64> = time
            .iter()
            .map(|&t| {
                2.0 * (2.0 * PI * 1.3 * t).sin() * (-0.4 * t).exp()
                    + 0.7 * (2.0 * PI * 4.1 * t + 0.3).cos()
            })
            .collect();
        let gm = GroundMotion { time, acceleration };
        for (zeta, rows) in REFERENCE_SPECTRA {
            let spectrum = compute_response_spectrum(&gm, &REFERENCE_PERIODS, zeta);
            for (i, &(sd, sv, sa)) in rows.iter().enumerate() {
                let period = REFERENCE_PERIODS[i];
                for (name, got, want) in [
                    ("Sd", spectrum.spectral_displacement[i], sd),
                    ("Sv", spectrum.spectral_velocity[i], sv),
                    ("Sa", spectrum.spectral_acceleration[i], sa),
                ] {
                    assert!(
                        (got - want).abs() <= 1.0e-9 * want,
                        "zeta {zeta} T {period}: {name} = {got}, reference {want}"
                    );
                }
            }
        }
    }
}
