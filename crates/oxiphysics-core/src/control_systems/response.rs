//! Exact sampled step and impulse responses of SISO continuous-time systems.
//!
//! A proper transfer function `num(s) / den(s)` is first converted to its
//! controllable canonical realisation `(A, B, C, D)`.  The realisation is then
//! discretised with the *exact* zero-order hold,
//!
//! ```text
//! exp([[A, B], [0, 0]] * dt) = [[Ad, Bd], [0, 1]],
//! ```
//!
//! and simulated with `x[k+1] = Ad x[k] + Bd u[k]`, `y[k] = C x[k] + D u[k]`.
//! For a step input the held input *is* the true input, so the samples equal
//! the exact response `y(k * dt)` up to rounding, for every pole
//! configuration: distinct, repeated, complex, unstable and integrating poles
//! are all treated alike (no partial-fraction expansion is needed).
//!
//! The arithmetic is complex so that the same engine serves transfer functions
//! with real coefficients (the usual case) and pole/zero lists that are not
//! closed under complex conjugation; the sampled output is the real part of
//! the response.
//!
//! The canonical form is built from the polynomial coefficients, so the
//! accuracy follows the conditioning of the coefficient representation: it is
//! at rounding level for the low and moderate orders met in control design,
//! and degrades for very high orders with strongly clustered or widely spread
//! poles, exactly as the coefficients themselves do.

use nalgebra::{Complex, DMatrix, DVector};

/// Complex scalar used by the realisation.
pub(crate) type C64 = Complex<f64>;

/// Leading polynomial coefficients whose modulus is below this fraction of the
/// largest coefficient are treated as exact zeros.  This removes the rounding
/// residue that polynomial arithmetic (for example a feedback loop whose
/// highest-order terms cancel) leaves in the leading position.
const LEADING_COEFFICIENT_RTOL: f64 = 1e-12;

/// Lifts real polynomial coefficients (descending powers) into the complex field.
pub(crate) fn real_polynomial(coefficients: &[f64]) -> Vec<C64> {
    coefficients.iter().map(|&c| C64::new(c, 0.0)).collect()
}

/// Expands `leading * prod_i (s - roots[i])` into descending-power coefficients.
pub(crate) fn polynomial_from_roots(leading: C64, roots: &[C64]) -> Vec<C64> {
    let mut coefficients = vec![leading];
    for &root in roots {
        let mut next = vec![C64::new(0.0, 0.0); coefficients.len() + 1];
        for (k, &c) in coefficients.iter().enumerate() {
            next[k] += c;
            next[k + 1] -= c * root;
        }
        coefficients = next;
    }
    coefficients
}

/// Drops leading (highest-power) coefficients that are numerically zero.
fn strip_leading_zeros(coefficients: &[C64]) -> &[C64] {
    let scale = coefficients
        .iter()
        .map(|c| c.norm())
        .fold(0.0_f64, f64::max);
    let threshold = LEADING_COEFFICIENT_RTOL * scale;
    let first = coefficients
        .iter()
        .position(|c| c.norm() > threshold)
        .unwrap_or(coefficients.len());
    &coefficients[first..]
}

/// Number of samples `ceil(t_final / dt)` of a response grid `t_k = k * dt`.
///
/// Returns `None` for a non-positive or non-finite `dt` and for a
/// non-positive or non-finite `t_final`.
pub(crate) fn sample_count(t_final: f64, dt: f64) -> Option<usize> {
    if !(dt.is_finite() && dt > 0.0 && t_final.is_finite() && t_final > 0.0) {
        return None;
    }
    let count = (t_final / dt).ceil();
    if count.is_finite() {
        Some(count as usize)
    } else {
        None
    }
}

/// Controllable canonical realisation of a proper SISO transfer function.
///
/// With the monic denominator `s^n + a_1 s^(n-1) + ... + a_n` the realisation
/// is
///
/// ```text
/// A = [[0, 1, 0, ...], ..., [-a_n, -a_(n-1), ..., -a_1]],  B = e_n,
/// C = [b_n - D a_n, ..., b_1 - D a_1],                      D = b_0,
/// ```
///
/// where `b_0, ..., b_n` are the numerator coefficients (descending powers,
/// zero-padded on the left) divided by the leading denominator coefficient.
#[derive(Debug, Clone)]
pub(crate) struct SisoRealization {
    /// System order `n` (number of states).
    order: usize,
    /// State matrix `A` (`n x n`).
    a: DMatrix<C64>,
    /// Input vector `B` (`n`).
    b: DVector<C64>,
    /// Output row `C`, stored as a column vector (`n`).
    c: DVector<C64>,
    /// Direct feedthrough `D`.
    d: C64,
}

impl SisoRealization {
    /// Builds the realisation of `num(s) / den(s)`; polynomial coefficients
    /// are in descending powers of `s`.
    ///
    /// Returns `None` when a coefficient is not finite, when the denominator
    /// is identically zero, or when the transfer function is improper
    /// (`deg num > deg den`), because the response of an improper system
    /// contains Dirac impulses and cannot be represented by samples.
    pub(crate) fn new(num: &[C64], den: &[C64]) -> Option<Self> {
        if num
            .iter()
            .chain(den.iter())
            .any(|c| !c.re.is_finite() || !c.im.is_finite())
        {
            return None;
        }
        let den = strip_leading_zeros(den);
        let num = strip_leading_zeros(num);
        let lead = *den.first()?;
        if num.len() > den.len() {
            return None;
        }
        let order = den.len() - 1;
        let zero = C64::new(0.0, 0.0);
        let one = C64::new(1.0, 0.0);

        let monic_den: Vec<C64> = den.iter().map(|&c| c / lead).collect();
        let mut monic_num = vec![zero; order + 1];
        let offset = order + 1 - num.len();
        for (i, &c) in num.iter().enumerate() {
            monic_num[offset + i] = c / lead;
        }
        let feedthrough = monic_num[0];

        let a = DMatrix::from_fn(order, order, |i, j| {
            if i + 1 == order {
                -monic_den[order - j]
            } else if j == i + 1 {
                one
            } else {
                zero
            }
        });
        let b = DVector::from_fn(order, |i, _| if i + 1 == order { one } else { zero });
        let c = DVector::from_fn(order, |j, _| {
            monic_num[order - j] - feedthrough * monic_den[order - j]
        });

        Some(Self {
            order,
            a,
            b,
            c,
            d: feedthrough,
        })
    }

    /// Exact zero-order-hold discretisation `(Ad, Bd)` for the sample time `dt`.
    fn discretize(&self, dt: f64) -> Option<(DMatrix<C64>, DVector<C64>)> {
        if !(dt.is_finite() && dt > 0.0) {
            return None;
        }
        let n = self.order;
        if n == 0 {
            return Some((DMatrix::zeros(0, 0), DVector::zeros(0)));
        }
        let mut augmented = DMatrix::from_element(n + 1, n + 1, C64::new(0.0, 0.0));
        for i in 0..n {
            for j in 0..n {
                augmented[(i, j)] = self.a[(i, j)] * dt;
            }
            augmented[(i, n)] = self.b[i] * dt;
        }
        if augmented
            .iter()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return None;
        }
        let exponential = augmented.exp();
        if exponential
            .iter()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return None;
        }
        let ad = DMatrix::from_fn(n, n, |i, j| exponential[(i, j)]);
        let bd = DVector::from_fn(n, |i, _| exponential[(i, n)]);
        Some((ad, bd))
    }

    /// Samples `y(k * dt)`, `k = 0..count`, of the unit-step response.
    ///
    /// The first sample is the direct feedthrough `D = y(0+)`.  `None` is
    /// returned if `dt` is invalid, the discretisation overflows or the
    /// `count` samples cannot be allocated.
    pub(crate) fn step_samples(&self, dt: f64, count: usize) -> Option<Vec<f64>> {
        let (ad, bd) = self.discretize(dt)?;
        let mut state = DVector::from_element(self.order, C64::new(0.0, 0.0));
        let mut samples = Vec::new();
        samples.try_reserve_exact(count).ok()?;
        for _ in 0..count {
            samples.push((self.c.dot(&state) + self.d).re);
            state = &ad * &state + &bd;
        }
        Some(samples)
    }

    /// Samples `y(k * dt)`, `k = 0..count`, of the unit-impulse response.
    ///
    /// The impulse sets the state to `B` at `t = 0+`, after which the system
    /// evolves freely.  A direct-feedthrough term `D * delta(t)` is a Dirac
    /// impulse that cannot be represented by samples and is therefore omitted.
    /// `None` is returned in the same situations as for
    /// [`step_samples`](Self::step_samples).
    pub(crate) fn impulse_samples(&self, dt: f64, count: usize) -> Option<Vec<f64>> {
        let (ad, _) = self.discretize(dt)?;
        let mut state = self.b.clone();
        let mut samples = Vec::new();
        samples.try_reserve_exact(count).ok()?;
        for _ in 0..count {
            samples.push(self.c.dot(&state).re);
            state = &ad * &state;
        }
        Some(samples)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::control_systems::{TransferFunction, impulse_response, step_response};

    /// Largest absolute deviation of `samples[k]` from `reference(k * dt)`.
    fn max_deviation(samples: &[f64], dt: f64, reference: impl Fn(f64) -> f64) -> f64 {
        samples
            .iter()
            .enumerate()
            .map(|(k, &y)| (y - reference(k as f64 * dt)).abs())
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn step_of_two_real_poles_matches_closed_form() {
        // 1 / ((s + 1)(s + 2)):  y(t) = 1/2 - e^-t + 1/2 e^-2t.
        let tf = TransferFunction::new(vec![1.0], vec![1.0, 3.0, 2.0]);
        let dt = 0.01;
        let (t, y) = step_response(&tf, 6.0, dt);
        assert_eq!(t.len(), 600);
        assert_eq!(y.len(), 600);
        let dev = max_deviation(&y, dt, |t| 0.5 - (-t).exp() + 0.5 * (-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
        // The old equal-weight modal sum gave 0.5 - 0.25 e^-t - 0.25 e^-2t.
        assert!((y[100] - (0.5 - (-1.0_f64).exp() + 0.5 * (-2.0_f64).exp())).abs() < 1e-10);
    }

    #[test]
    fn step_of_nontrivial_numerator_uses_the_zero_correctly() {
        // (s + 3) / ((s + 1)(s + 2)):  y(t) = 3/2 - 2 e^-t + 1/2 e^-2t.
        // A reversed C vector in the canonical form would fail this test.
        let tf = TransferFunction::new(vec![1.0, 3.0], vec![1.0, 3.0, 2.0]);
        let dt = 0.02;
        let (_, y) = step_response(&tf, 8.0, dt);
        let dev = max_deviation(&y, dt, |t| 1.5 - 2.0 * (-t).exp() + 0.5 * (-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
        assert!(y[0].abs() < 1e-12, "strictly proper => y(0) = 0");
    }

    #[test]
    fn step_of_non_monic_first_order_system() {
        // 2 / (0.5 s + 1) = 4 / (s + 2):  y(t) = 2 (1 - e^-2t).
        let tf = TransferFunction::first_order(2.0, 0.5);
        let dt = 0.01;
        let (_, y) = step_response(&tf, 3.0, dt);
        let dev = max_deviation(&y, dt, |t| 2.0 * (1.0 - (-2.0 * t).exp()));
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn step_of_underdamped_second_order_system() {
        let (wn, zeta) = (5.0_f64, 0.3_f64);
        let tf = TransferFunction::second_order(wn, zeta);
        let dt = 0.005;
        let (_, y) = step_response(&tf, 4.0, dt);
        let wd = wn * (1.0 - zeta * zeta).sqrt();
        let phi = zeta.acos();
        let analytic = |t: f64| {
            1.0 - (-zeta * wn * t).exp() / (1.0 - zeta * zeta).sqrt() * (wd * t + phi).sin()
        };
        let dev = max_deviation(&y, dt, analytic);
        assert!(dev < 1e-9, "max deviation {dev}");
        // Overshoot of the underdamped step response: exp(-pi zeta / sqrt(1 - zeta^2)).
        let peak = y.iter().cloned().fold(f64::MIN, f64::max);
        let overshoot = (-std::f64::consts::PI * zeta / (1.0 - zeta * zeta).sqrt()).exp();
        assert!(
            (peak - (1.0 + overshoot)).abs() < 2e-3,
            "peak {peak}, expected {}",
            1.0 + overshoot
        );
    }

    #[test]
    fn step_of_critically_damped_system_with_repeated_pole() {
        // wn^2 / (s + wn)^2: y(t) = 1 - (1 + wn t) e^(-wn t).  Partial
        // fractions with distinct poles break down here; the exact ZOH does not.
        let wn = 3.0_f64;
        let tf = TransferFunction::second_order(wn, 1.0);
        let dt = 0.01;
        let (_, y) = step_response(&tf, 5.0, dt);
        let dev = max_deviation(&y, dt, |t| 1.0 - (1.0 + wn * t) * (-wn * t).exp());
        assert!(dev < 1e-9, "max deviation {dev}");
    }

    #[test]
    fn step_of_integrators_and_unstable_pole() {
        let dt = 0.01;
        // 1 / s: y = t.
        let integrator = TransferFunction::new(vec![1.0], vec![1.0, 0.0]);
        let (_, y) = step_response(&integrator, 2.0, dt);
        assert!(max_deviation(&y, dt, |t| t) < 1e-10);
        // 1 / s^2: y = t^2 / 2.
        let double_integrator = TransferFunction::new(vec![1.0], vec![1.0, 0.0, 0.0]);
        let (_, y) = step_response(&double_integrator, 2.0, dt);
        assert!(max_deviation(&y, dt, |t| 0.5 * t * t) < 1e-10);
        // 1 / (s - 1): y = e^t - 1.
        let unstable = TransferFunction::new(vec![1.0], vec![1.0, -1.0]);
        let (_, y) = step_response(&unstable, 2.0, dt);
        assert!(max_deviation(&y, dt, |t| t.exp() - 1.0) < 1e-9);
    }

    #[test]
    fn step_of_biproper_system_starts_at_feedthrough() {
        // (s + 1) / (s + 2) = 1 - 1 / (s + 2):  y(t) = 1/2 + 1/2 e^-2t, y(0+) = 1.
        let tf = TransferFunction::new(vec![1.0, 1.0], vec![1.0, 2.0]);
        let dt = 0.01;
        let (_, y) = step_response(&tf, 3.0, dt);
        assert!((y[0] - 1.0).abs() < 1e-12, "y(0) = {}", y[0]);
        let dev = max_deviation(&y, dt, |t| 0.5 + 0.5 * (-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn impulse_of_first_order_system_is_scaled_exponential() {
        // 3 / (s + 1): g(t) = 3 e^-t (the old code returned e^-t).
        let tf = TransferFunction::new(vec![3.0], vec![1.0, 1.0]);
        let dt = 0.01;
        let (t, y) = impulse_response(&tf, 5.0, dt);
        assert_eq!(t.len(), y.len());
        assert!((y[0] - 3.0).abs() < 1e-12, "g(0+) = {}", y[0]);
        let dev = max_deviation(&y, dt, |t| 3.0 * (-t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn impulse_of_nontrivial_numerator_is_derivative_of_step() {
        // (s + 3) / ((s + 1)(s + 2)): g(t) = 2 e^-t - e^-2t.
        let tf = TransferFunction::new(vec![1.0, 3.0], vec![1.0, 3.0, 2.0]);
        let dt = 0.01;
        let (_, y) = impulse_response(&tf, 6.0, dt);
        let dev = max_deviation(&y, dt, |t| 2.0 * (-t).exp() - (-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn impulse_of_biproper_system_omits_the_dirac_term() {
        // (s + 1) / (s + 2) = 1 - 1 / (s + 2): the regular part is -e^-2t.
        let tf = TransferFunction::new(vec![1.0, 1.0], vec![1.0, 2.0]);
        let dt = 0.01;
        let (_, y) = impulse_response(&tf, 3.0, dt);
        let dev = max_deviation(&y, dt, |t| -(-2.0 * t).exp());
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn impulse_response_of_oscillator_matches_sine() {
        // wn^2 / (s^2 + wn^2): g(t) = wn sin(wn t).
        let wn = 2.0_f64;
        let tf = TransferFunction::new(vec![wn * wn], vec![1.0, 0.0, wn * wn]);
        let dt = 0.01;
        let (_, y) = impulse_response(&tf, 6.0, dt);
        let dev = max_deviation(&y, dt, |t| wn * (wn * t).sin());
        assert!(dev < 1e-9, "max deviation {dev}");
    }

    #[test]
    fn response_time_grid_matches_ceil_of_ratio() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (t, y) = step_response(&tf, 1.0, 0.3);
        assert_eq!(t.len(), 4);
        assert_eq!(y.len(), 4);
        for (k, &tk) in t.iter().enumerate() {
            assert!((tk - k as f64 * 0.3).abs() < 1e-15);
        }
    }

    /// A transfer function with reference step/impulse samples at t = 0.5, 1, 2, 4.
    struct ReferenceCase {
        num: Vec<f64>,
        den: Vec<f64>,
        step: [f64; 4],
        impulse: [f64; 4],
    }

    #[test]
    fn higher_order_systems_match_residue_expansion_references() {
        // Reference values come from an independent partial-fraction
        // evaluation (scipy.signal.residue) of the closed-form responses,
        // cross-checked against scipy.signal.lti; they pin down the full
        // canonical form (A, B, C) for n = 3 and n = 4 with complex poles.
        let dt = 0.01;
        let samples = [50usize, 100, 200, 400]; // t = 0.5, 1, 2, 4
        let cases = [
            ReferenceCase {
                // (s + 1) / (s^4 + 3 s^3 + 9 s^2 + 9 s + 10)
                num: vec![1.0, 1.0],
                den: vec![1.0, 3.0, 9.0, 9.0, 10.0],
                step: [
                    1.52004075548920e-02,
                    7.86611198528086e-02,
                    1.89284227434913e-01,
                    6.94746958322508e-02,
                ],
                impulse: [
                    7.97258636361910e-02,
                    1.58532276009514e-01,
                    1.64423671305674e-02,
                    -1.74559298347799e-02,
                ],
            },
            ReferenceCase {
                // (2 s + 1) / (s^3 + 2 s^2 + 5 s + 3)
                num: vec![2.0, 1.0],
                den: vec![1.0, 2.0, 5.0, 3.0],
                step: [
                    1.81314490965462e-01,
                    4.58245180503746e-01,
                    5.02250745222219e-01,
                    3.26914044351406e-01,
                ],
                impulse: [
                    5.83432835676185e-01,
                    4.34089422282111e-01,
                    -2.44331466798420e-01,
                    7.66669616116195e-02,
                ],
            },
        ];
        for case in cases {
            let tf = TransferFunction::new(case.num.clone(), case.den.clone());
            let (_, step) = step_response(&tf, 4.005, dt);
            let (_, impulse) = impulse_response(&tf, 4.005, dt);
            for (idx, &k) in samples.iter().enumerate() {
                assert!(
                    (step[k] - case.step[idx]).abs() < 1e-10,
                    "den {:?} sample {idx}: step {} vs {}",
                    case.den,
                    step[k],
                    case.step[idx]
                );
                assert!(
                    (impulse[k] - case.impulse[idx]).abs() < 1e-10,
                    "den {:?} sample {idx}: impulse {} vs {}",
                    case.den,
                    impulse[k],
                    case.impulse[idx]
                );
            }
        }
    }

    #[test]
    fn step_matches_integral_of_impulse_response() {
        // Cross-check of the two code paths: for a strictly proper system the
        // step response is the running integral of the impulse response.
        let tf = TransferFunction::new(vec![2.0, 1.0], vec![1.0, 2.0, 5.0, 3.0]);
        let dt = 0.001;
        let (_, step) = step_response(&tf, 4.0, dt);
        let (_, impulse) = impulse_response(&tf, 4.0, dt);
        let mut integral = 0.0;
        for k in 1..step.len() {
            integral += 0.5 * dt * (impulse[k - 1] + impulse[k]);
            assert!(
                (integral - step[k]).abs() < 5e-7,
                "k = {k}: integral {integral} vs step {}",
                step[k]
            );
        }
    }

    #[test]
    fn invalid_requests_give_empty_responses() {
        let tf = TransferFunction::first_order(1.0, 1.0);
        for (t_final, dt) in [
            (1.0, 0.0),
            (1.0, -0.1),
            (1.0, f64::NAN),
            (1.0, f64::INFINITY),
            (0.0, 0.1),
            (-1.0, 0.1),
            (f64::NAN, 0.1),
            (f64::INFINITY, 0.1),
        ] {
            let (t, y) = step_response(&tf, t_final, dt);
            assert!(t.is_empty() && y.is_empty(), "step ({t_final}, {dt})");
            let (t, y) = impulse_response(&tf, t_final, dt);
            assert!(t.is_empty() && y.is_empty(), "impulse ({t_final}, {dt})");
        }
        // Improper (s^2 / (s + 1)): no sampled response exists.
        let improper = TransferFunction::new(vec![1.0, 0.0, 0.0], vec![1.0, 1.0]);
        let (t, y) = step_response(&improper, 1.0, 0.1);
        assert!(t.is_empty() && y.is_empty());
        // Zero denominator.
        let singular = TransferFunction::new(vec![1.0], vec![0.0, 0.0]);
        let (t, y) = step_response(&singular, 1.0, 0.1);
        assert!(t.is_empty() && y.is_empty());
        // Non-finite coefficient.
        let nan_tf = TransferFunction::new(vec![f64::NAN], vec![1.0, 1.0]);
        let (t, y) = impulse_response(&nan_tf, 1.0, 0.1);
        assert!(t.is_empty() && y.is_empty());
    }

    #[test]
    fn unallocatable_sample_counts_give_empty_responses() {
        // ceil(t_final / dt) saturates at usize::MAX: the sample buffer cannot
        // exist, which must be reported as "no response", not as a panic.
        let tf = TransferFunction::first_order(1.0, 1.0);
        let (t, y) = step_response(&tf, 1e300, 1.0);
        assert!(t.is_empty() && y.is_empty());
        let (t, y) = impulse_response(&tf, 1e300, 1.0);
        assert!(t.is_empty() && y.is_empty());
    }

    #[test]
    fn zero_numerator_gives_zero_response() {
        let tf = TransferFunction::new(vec![0.0], vec![1.0, 1.0]);
        let (_, y) = step_response(&tf, 1.0, 0.1);
        assert!(y.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn constant_transfer_function_is_a_pure_gain() {
        let tf = TransferFunction::new(vec![2.5], vec![1.0]);
        let (_, y) = step_response(&tf, 1.0, 0.25);
        assert_eq!(y, vec![2.5; 4]);
        let (_, g) = impulse_response(&tf, 1.0, 0.25);
        assert_eq!(g, vec![0.0; 4]);
    }

    #[test]
    fn leading_rounding_residue_of_denominator_is_ignored() {
        // A feedback loop whose s^2 terms cancel leaves ~1e-17 in front of s^2.
        let tf = TransferFunction::new(vec![1.0], vec![1e-17, 1.0, 1.0]);
        let (_, y) = step_response(&tf, 5.0, 0.01);
        // Effective system 1 / (s + 1): y = 1 - e^-t.
        let dev = max_deviation(&y, 0.01, |t| 1.0 - (-t).exp());
        assert!(dev < 1e-9, "max deviation {dev}");
    }

    #[test]
    fn complex_conjugate_poles_from_roots_give_real_response() {
        // Poles -1 +- 2j, gain 5 (so G(0) = 1): y(t) = 1 - e^-t (cos 2t + sin(2t) / 2).
        let poles = [C64::new(-1.0, 2.0), C64::new(-1.0, -2.0)];
        let den = polynomial_from_roots(C64::new(1.0, 0.0), &poles);
        for c in &den {
            assert!(c.im.abs() < 1e-15, "denominator must be real: {c}");
        }
        let realization = SisoRealization::new(&[C64::new(5.0, 0.0)], &den)
            .expect("proper system must have a realisation");
        let dt = 0.01;
        let y = realization
            .step_samples(dt, 500)
            .expect("valid sample time");
        let dev = max_deviation(&y, dt, |t| {
            1.0 - (-t).exp() * ((2.0 * t).cos() + 0.5 * (2.0 * t).sin())
        });
        assert!(dev < 1e-10, "max deviation {dev}");
    }

    #[test]
    fn polynomial_from_roots_expands_in_descending_powers() {
        let p = polynomial_from_roots(
            C64::new(2.0, 0.0),
            &[C64::new(-1.0, 0.0), C64::new(-2.0, 0.0), C64::new(0.0, 0.0)],
        );
        // 2 (s + 1)(s + 2) s = 2 s^3 + 6 s^2 + 4 s.
        let expected = [2.0, 6.0, 4.0, 0.0];
        assert_eq!(p.len(), 4);
        for (c, e) in p.iter().zip(expected.iter()) {
            assert!((c.re - e).abs() < 1e-14 && c.im.abs() < 1e-14, "{c} vs {e}");
        }
    }

    #[test]
    fn sample_count_rejects_degenerate_grids() {
        assert_eq!(sample_count(1.0, 0.25), Some(4));
        assert_eq!(sample_count(1.0, 0.3), Some(4));
        assert_eq!(sample_count(1.0, 0.0), None);
        assert_eq!(sample_count(0.0, 0.1), None);
        assert_eq!(sample_count(f64::NAN, 0.1), None);
    }
}
