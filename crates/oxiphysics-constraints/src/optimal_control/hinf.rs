//! H-infinity norm of a continuous-time LTI system.
//!
//! For a stable system `G(s) = C (sI - A)⁻¹ B + D`,
//!
//! ```text
//! ‖G‖∞ = sup_ω σ_max(G(jω))
//! ```
//!
//! is computed with the quadratically convergent algorithm of Bruinsma and
//! Steinbuch (1990), which refines the bisection of Boyd, Balakrishnan and
//! Kabamba (1989).  For a level `γ > σ_max(D)` the Hamiltonian matrix
//!
//! ```text
//! H(γ) = [ A + B R⁻¹ Dᵀ C        B R⁻¹ Bᵀ              ]      R = γ² I - Dᵀ D
//!        [ -Cᵀ (I + D R⁻¹ Dᵀ) C  -(A + B R⁻¹ Dᵀ C)ᵀ    ]
//! ```
//!
//! has an eigenvalue `jω` on the imaginary axis exactly when some singular
//! value of `G(jω)` equals `γ`.  Starting from a lower bound `lb` taken from
//! actual frequency-response evaluations, the algorithm tests `H(γ)` at
//! `γ = (1 + 2·tol)·lb`:
//!
//! * no imaginary-axis eigenvalue: `γ` is a certified upper bound and the
//!   norm lies in `[lb, γ]`;
//! * otherwise the norm exceeds `γ` somewhere between the crossing
//!   frequencies; `σ_max` at the crossings and at the midpoints between
//!   consecutive crossings raises `lb` and the test is repeated.
//!
//! The lower bound is always a value of `σ_max(G(jω))` that was really
//! evaluated, so it never exceeds the norm; the upper bound is only reported
//! when the Hamiltonian test certifies it.

use super::functions::{mat_mul_rect, mat_transpose_rect};
use crate::dense_linalg::{
    cholesky_factor, cholesky_solve, lu_solve, real_eigenvalues, singular_values,
};

/// Relative gap between the certified upper bound and the lower bound.
const RELATIVE_TOLERANCE: f64 = 1e-6;
/// An eigenvalue `λ` counts as lying on the imaginary axis when
/// `|Re λ| <= IMAGINARY_AXIS_TOLERANCE · |λ|` (plus a rounding floor).
const IMAGINARY_AXIS_TOLERANCE: f64 = 1e-8;
/// `A` must have all eigenvalues left of `Re λ = -STABILITY_MARGIN · ‖A‖_F`;
/// eigenvalues closer to the axis are indistinguishable from marginal ones at
/// rounding level (computed eigenvalues of an undamped oscillator are off the
/// axis by about `ε ‖A‖`).  The margin is deliberately tiny so that stiff but
/// stable plants (a slow pole next to fast dynamics) are not rejected; a pole
/// that is merely very close to the axis simply yields a huge norm.
const STABILITY_MARGIN: f64 = 1e-12;
/// Upper limit on Bruinsma–Steinbuch iterations (typically 1 to 5 are needed).
const MAX_ITERATIONS: usize = 100;
/// Number of log-spaced frequencies tried when the pole-based start gives zero.
const FALLBACK_SWEEP_POINTS: usize = 61;

/// Bounds on the H-infinity norm of a continuous-time LTI system.
///
/// `lower <= ‖G‖∞ <= upper`.  `lower` is the largest value of
/// `σ_max(G(jω))` that was evaluated, `upper` is certified by the Hamiltonian
/// test and exceeds `lower` by a relative `3e-6` at most.  For a system whose
/// state matrix is not Hurwitz (an eigenvalue on or right of the imaginary axis,
/// up to a rounding-level margin) the norm is infinite: both bounds are
/// `f64::INFINITY`.  `upper` is also infinite if the iteration could not
/// certify a bound (which would indicate a numerical breakdown).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HInfBounds {
    /// Lower bound: the peak of `σ_max(G(jω))` over the evaluated frequencies.
    pub lower: f64,
    /// Certified upper bound on the norm.
    pub upper: f64,
    /// Frequency (rad/s) at which `lower` was attained.  `f64::INFINITY` when
    /// the norm is the high-frequency gain `σ_max(D)`, `f64::NAN` when the
    /// system is not stable.
    pub peak_frequency: f64,
}

impl HInfBounds {
    /// Bounds of a system without a finite H-infinity norm.
    const UNBOUNDED: Self = Self {
        lower: f64::INFINITY,
        upper: f64::INFINITY,
        peak_frequency: f64::NAN,
    };
}

/// Borrowed state-space realisation `(A, B, C, D)` with row-major matrices:
/// `A` is `n x n`, `B` is `n x m`, `C` is `p x n` and `D` is `p x m`.
pub(super) struct StateSpace<'a> {
    pub(super) n: usize,
    pub(super) m: usize,
    pub(super) p: usize,
    pub(super) a: &'a [f64],
    pub(super) b: &'a [f64],
    pub(super) c: &'a [f64],
    pub(super) d: &'a [f64],
}

fn frobenius_norm(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

impl StateSpace<'_> {
    /// All matrices have the declared sizes and only finite entries.
    fn is_well_formed(&self) -> bool {
        self.a.len() == self.n * self.n
            && self.b.len() == self.n * self.m
            && self.c.len() == self.p * self.n
            && self.d.len() == self.p * self.m
            && [self.a, self.b, self.c, self.d]
                .iter()
                .all(|matrix| matrix.iter().all(|v| v.is_finite()))
    }

    /// Largest singular value of `D` (zero for an empty matrix).
    fn sigma_max_feedthrough(&self) -> f64 {
        singular_values(self.d, self.p, self.m)
            .first()
            .copied()
            .unwrap_or(0.0)
    }

    /// `σ_max(G(jω))`.
    ///
    /// `(jωI - A) X = B` is solved as the real system
    /// `[[-A, -ωI], [ωI, -A]] [Xr; Xi] = [B; 0]`, which gives
    /// `G = (C Xr + D) + j C Xi`; the largest singular value of the complex
    /// `p x m` matrix equals that of the real embedding
    /// `[[Re G, -Im G], [Im G, Re G]]`.  Returns `f64::INFINITY` when
    /// `jωI - A` is singular to working precision (a pole on the axis).
    pub(super) fn sigma_max(&self, omega: f64) -> f64 {
        let (n, m, p) = (self.n, self.m, self.p);
        let n2 = 2 * n;
        let mut k = vec![0.0; n2 * n2];
        for i in 0..n {
            for j in 0..n {
                let aij = self.a[i * n + j];
                k[i * n2 + j] = -aij;
                k[(n + i) * n2 + n + j] = -aij;
            }
            k[i * n2 + n + i] = -omega;
            k[(n + i) * n2 + i] = omega;
        }
        let mut rhs = vec![0.0; n2 * m];
        rhs[..n * m].copy_from_slice(self.b);
        let Some(x) = lu_solve(&k, n2, &rhs, m) else {
            return f64::INFINITY;
        };
        let (xr, xi) = x.split_at(n * m);
        let mut g_re = mat_mul_rect(self.c, p, n, xr, m);
        let g_im = mat_mul_rect(self.c, p, n, xi, m);
        for (g, d) in g_re.iter_mut().zip(self.d.iter()) {
            *g += d;
        }
        if p == 1 && m == 1 {
            return g_re[0].hypot(g_im[0]);
        }
        let (rows, cols) = (2 * p, 2 * m);
        let mut embedding = vec![0.0; rows * cols];
        for r in 0..p {
            for c in 0..m {
                let (re, im) = (g_re[r * m + c], g_im[r * m + c]);
                embedding[r * cols + c] = re;
                embedding[r * cols + m + c] = -im;
                embedding[(p + r) * cols + c] = im;
                embedding[(p + r) * cols + m + c] = re;
            }
        }
        singular_values(&embedding, rows, cols)
            .first()
            .copied()
            .unwrap_or(0.0)
    }

    /// The Hamiltonian matrix `H(γ)` (row-major `2n x 2n`); `None` unless
    /// `γ > σ_max(D)` so that `R = γ² I - Dᵀ D` is positive definite.
    pub(super) fn hamiltonian(&self, gamma: f64) -> Option<Vec<f64>> {
        let (n, m, p) = (self.n, self.m, self.p);
        let dt = mat_transpose_rect(self.d, p, m);
        let dtd = mat_mul_rect(&dt, m, p, self.d, m);
        let mut r = vec![0.0; m * m];
        let mut identity_m = vec![0.0; m * m];
        for i in 0..m {
            for j in 0..m {
                r[i * m + j] = -dtd[i * m + j];
            }
            r[i * m + i] += gamma * gamma;
            identity_m[i * m + i] = 1.0;
        }
        let r_factor = cholesky_factor(&r, m)?;
        let r_inv = cholesky_solve(&r_factor, m, &identity_m, m);
        let bt = mat_transpose_rect(self.b, n, m);
        let ct = mat_transpose_rect(self.c, p, n);
        let b_rinv = mat_mul_rect(self.b, n, m, &r_inv, m);
        // F = A + B R⁻¹ Dᵀ C
        let dt_c = mat_mul_rect(&dt, m, p, self.c, n);
        let correction = mat_mul_rect(&b_rinv, n, m, &dt_c, n);
        let f: Vec<f64> = self
            .a
            .iter()
            .zip(correction.iter())
            .map(|(a, c)| a + c)
            .collect();
        // G = B R⁻¹ Bᵀ
        let g = mat_mul_rect(&b_rinv, n, m, &bt, n);
        // Q = Cᵀ (I + D R⁻¹ Dᵀ) C
        let d_rinv = mat_mul_rect(self.d, p, m, &r_inv, m);
        let mut w = mat_mul_rect(&d_rinv, p, m, &dt, p);
        for i in 0..p {
            w[i * p + i] += 1.0;
        }
        let w_c = mat_mul_rect(&w, p, p, self.c, n);
        let q = mat_mul_rect(&ct, n, p, &w_c, n);
        let n2 = 2 * n;
        let mut h = vec![0.0; n2 * n2];
        for i in 0..n {
            for j in 0..n {
                h[i * n2 + j] = f[i * n + j];
                h[i * n2 + n + j] = g[i * n + j];
                h[(n + i) * n2 + j] = -q[i * n + j];
                h[(n + i) * n2 + n + j] = -f[j * n + i];
            }
        }
        Some(h)
    }

    /// Frequencies `ω >= 0` at which `H(γ)` has an eigenvalue `jω` on the
    /// imaginary axis, sorted ascending; `None` if the eigenvalue iteration
    /// failed.
    pub(super) fn crossing_frequencies(&self, gamma: f64) -> Option<Vec<f64>> {
        let h = self.hamiltonian(gamma)?;
        let eigenvalues = real_eigenvalues(&h, 2 * self.n)?;
        let floor = 1e-13 * frobenius_norm(&h);
        let mut frequencies: Vec<f64> = eigenvalues
            .iter()
            .filter(|&&(re, im)| {
                im >= 0.0 && re.abs() <= IMAGINARY_AXIS_TOLERANCE * re.hypot(im) + floor
            })
            .map(|&(_, im)| im)
            .collect();
        frequencies.sort_by(f64::total_cmp);
        frequencies.dedup();
        Some(frequencies)
    }

    /// Bounds on `‖G‖∞`, see the module documentation.
    pub(super) fn h_inf_bounds(&self) -> HInfBounds {
        if !self.is_well_formed() {
            return HInfBounds::UNBOUNDED;
        }
        if self.m == 0 || self.p == 0 {
            return HInfBounds {
                lower: 0.0,
                upper: 0.0,
                peak_frequency: 0.0,
            };
        }
        let sigma_d = self.sigma_max_feedthrough();
        if self.n == 0 {
            return HInfBounds {
                lower: sigma_d,
                upper: sigma_d,
                peak_frequency: f64::INFINITY,
            };
        }
        let Some(poles) = real_eigenvalues(self.a, self.n) else {
            return HInfBounds::UNBOUNDED;
        };
        let a_scale = frobenius_norm(self.a);
        if poles
            .iter()
            .any(|&(re, _)| re >= -STABILITY_MARGIN * a_scale)
        {
            return HInfBounds::UNBOUNDED;
        }

        // Lower bound from DC, the high-frequency gain and the resonance
        // frequencies suggested by the poles: |λ|, |Im λ| and the peak
        // frequency √(|λ|² - 2 Re²) of a second-order mode.
        let mut lower = sigma_d;
        let mut peak = f64::INFINITY;
        let mut candidates = vec![0.0];
        for &(re, im) in &poles {
            let magnitude = re.hypot(im);
            candidates.push(magnitude);
            candidates.push(im.abs());
            candidates.push((magnitude * magnitude - 2.0 * re * re).max(0.0).sqrt());
        }
        for &omega in &candidates {
            let sigma = self.sigma_max(omega);
            if sigma > lower {
                lower = sigma;
                peak = omega;
            }
        }
        if lower == 0.0 {
            // The pole-based start missed every non-zero value: sweep.
            let magnitudes = poles.iter().map(|&(re, im)| re.hypot(im));
            let slowest = magnitudes.clone().fold(f64::INFINITY, f64::min).max(1e-12);
            let fastest = magnitudes.fold(0.0_f64, f64::max).max(slowest);
            let (lo, hi) = ((slowest * 1e-2).ln(), (fastest * 1e2).ln());
            for k in 0..FALLBACK_SWEEP_POINTS {
                let omega = (lo + (hi - lo) * k as f64 / (FALLBACK_SWEEP_POINTS - 1) as f64).exp();
                let sigma = self.sigma_max(omega);
                if sigma > lower {
                    lower = sigma;
                    peak = omega;
                }
            }
            if lower == 0.0 {
                return HInfBounds {
                    lower: 0.0,
                    upper: 0.0,
                    peak_frequency: 0.0,
                };
            }
        }

        for _ in 0..MAX_ITERATIONS {
            let gamma = (1.0 + 2.0 * RELATIVE_TOLERANCE) * lower;
            let Some(crossings) = self.crossing_frequencies(gamma) else {
                return HInfBounds {
                    lower,
                    upper: f64::INFINITY,
                    peak_frequency: peak,
                };
            };
            if crossings.is_empty() {
                return HInfBounds {
                    lower,
                    upper: gamma,
                    peak_frequency: peak,
                };
            }
            // The norm exceeds gamma between the crossings: probe them and the
            // midpoints of consecutive crossings.
            let mut probes = crossings.clone();
            probes.extend(crossings.windows(2).map(|w| 0.5 * (w[0] + w[1])));
            let mut improved = lower;
            let mut improved_at = peak;
            for &omega in &probes {
                let sigma = self.sigma_max(omega);
                if sigma > improved {
                    improved = sigma;
                    improved_at = omega;
                }
            }
            if improved > lower * (1.0 + 0.5 * RELATIVE_TOLERANCE) {
                lower = improved;
                peak = improved_at;
            } else {
                // Only near-double eigenvalues just above a peak that is
                // already within tolerance of `lower`: nothing exceeds gamma.
                return HInfBounds {
                    lower: improved,
                    upper: gamma * (1.0 + RELATIVE_TOLERANCE),
                    peak_frequency: improved_at,
                };
            }
        }
        HInfBounds {
            lower,
            upper: f64::INFINITY,
            peak_frequency: peak,
        }
    }
}
