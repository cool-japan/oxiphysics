// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Complex analysis, special functions, and conformal mappings.
//!
//! # Overview
//!
//! Provides:
//!
//! - [`Complex`] number type with full arithmetic and transcendental functions
//! - Cauchy integral formula and residue computation
//! - Conformal mappings (Joukowski, Möbius transformation)
//! - Contour integration and winding numbers
//! - Gamma function (Lanczos approximation), Beta function
//! - Riemann zeta function (Euler–Maclaurin)
//! - Complete elliptic integrals K(k) and E(k)
//! - Bessel functions J_n(x) and Y_n(x)
//! - Legendre polynomials P_n(x) and Hermite polynomials H_n(x)

use std::f64::consts::{FRAC_2_PI, PI};

// ---------------------------------------------------------------------------
// Complex number type
// ---------------------------------------------------------------------------

/// A complex number z = re + i·im.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// Construct a new complex number.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// The purely real number `re + 0i`.
    pub fn from_real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }

    /// The purely imaginary number `0 + im·i`.
    pub fn from_imag(im: f64) -> Self {
        Self { re: 0.0, im }
    }

    /// The additive identity 0.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    /// The multiplicative identity 1.
    pub fn one() -> Self {
        Self::new(1.0, 0.0)
    }

    /// The imaginary unit i.
    pub fn i() -> Self {
        Self::new(0.0, 1.0)
    }

    /// Absolute value (modulus) |z| = sqrt(re² + im²).
    pub fn abs(&self) -> f64 {
        self.re.hypot(self.im)
    }

    /// Argument (phase angle) arg(z) = atan2(im, re) in radians.
    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Complex conjugate z̄ = re − i·im.
    pub fn conj(&self) -> Self {
        Self::new(self.re, -self.im)
    }

    /// Squared modulus |z|².
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Multiplicative inverse 1/z.
    ///
    /// Returns `NaN`-valued complex if z == 0.
    pub fn inv(&self) -> Self {
        let n = self.norm_sq();
        Self::new(self.re / n, -self.im / n)
    }

    /// Complex exponential e^z = e^re · (cos(im) + i·sin(im)).
    pub fn exp(&self) -> Self {
        let r = self.re.exp();
        Self::new(r * self.im.cos(), r * self.im.sin())
    }

    /// Principal natural logarithm ln(z) = ln|z| + i·arg(z).
    pub fn ln(&self) -> Self {
        Self::new(self.abs().ln(), self.arg())
    }

    /// Principal square root √z.
    pub fn sqrt(&self) -> Self {
        let r = self.abs().sqrt();
        let theta = self.arg() / 2.0;
        Self::new(r * theta.cos(), r * theta.sin())
    }

    /// Raise to a real power: z^n = |z|^n · e^(i·n·arg(z)).
    pub fn pow(&self, n: f64) -> Self {
        if self.re == 0.0 && self.im == 0.0 {
            return Self::zero();
        }
        let r = self.abs().powf(n);
        let theta = self.arg() * n;
        Self::new(r * theta.cos(), r * theta.sin())
    }

    /// Raise to a complex power: z^w = e^(w · ln z).
    pub fn cpow(&self, w: Self) -> Self {
        (w * self.ln()).exp()
    }

    /// Complex sine sin(z) = sin(re)·cosh(im) + i·cos(re)·sinh(im).
    pub fn sin(&self) -> Self {
        Self::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh(),
        )
    }

    /// Complex cosine cos(z) = cos(re)·cosh(im) − i·sin(re)·sinh(im).
    pub fn cos(&self) -> Self {
        Self::new(
            self.re.cos() * self.im.cosh(),
            -self.re.sin() * self.im.sinh(),
        )
    }

    /// Complex hyperbolic sine sinh(z) = sinh(re)·cos(im) + i·cosh(re)·sin(im).
    pub fn sinh(&self) -> Self {
        Self::new(
            self.re.sinh() * self.im.cos(),
            self.re.cosh() * self.im.sin(),
        )
    }

    /// Complex hyperbolic cosine cosh(z) = cosh(re)·cos(im) + i·sinh(re)·sin(im).
    pub fn cosh(&self) -> Self {
        Self::new(
            self.re.cosh() * self.im.cos(),
            self.re.sinh() * self.im.sin(),
        )
    }

    /// Complex tangent tan(z) = sin(z)/cos(z).
    pub fn tan(&self) -> Self {
        self.sin() / self.cos()
    }

    /// Test approximate equality within tolerance.
    pub fn approx_eq(&self, other: Self, tol: f64) -> bool {
        (self.re - other.re).abs() < tol && (self.im - other.im).abs() < tol
    }
}

// Arithmetic operator implementations

impl std::ops::Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

impl std::ops::Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

impl std::ops::Add<f64> for Complex {
    type Output = Self;
    fn add(self, rhs: f64) -> Self {
        Self::new(self.re + rhs, self.im)
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl std::ops::Div<f64> for Complex {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

impl std::fmt::Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.im >= 0.0 {
            write!(f, "{} + {}i", self.re, self.im)
        } else {
            write!(f, "{} - {}i", self.re, -self.im)
        }
    }
}

// ---------------------------------------------------------------------------
// Cauchy integral formula
// ---------------------------------------------------------------------------

/// Compute the Cauchy contour integral ∮ f(z) dz along a polygonal contour.
///
/// The contour is given as a sequence of vertices; the integral is computed
/// using the trapezoidal rule over the straight segments connecting them.
/// The contour is automatically closed (last vertex connected back to first).
///
/// # Arguments
/// * `f`       - analytic function to integrate
/// * `contour` - ordered list of waypoints on the contour (polygonal path)
///
/// # Returns
/// Approximation of the contour integral.
pub fn cauchy_integral(f: impl Fn(Complex) -> Complex, contour: &[Complex]) -> Complex {
    if contour.len() < 2 {
        return Complex::zero();
    }
    let n = contour.len();
    let mut result = Complex::zero();
    for i in 0..n {
        let z0 = contour[i];
        let z1 = contour[(i + 1) % n];
        let dz = z1 - z0;
        // trapezoidal rule: (f(z0) + f(z1)) / 2 * dz
        let mid = (f(z0) + f(z1)) / 2.0;
        result = result + mid * dz;
    }
    result
}

/// Default radius of the circle used by [`residue`].
///
/// Small enough to enclose only the pole of interest for typical functions,
/// large enough that the cancellation in the principal part (terms of size
/// `r^{1-m}` for an order-`m` pole) stays far below the 1e-8 level.
const RESIDUE_DEFAULT_RADIUS: f64 = 1e-3;

/// Minimum number of quadrature nodes on the residue circle.
const RESIDUE_MIN_POINTS: usize = 512;

/// Compute the residue of `f` at a (simple or higher-order) pole.
///
/// The residue is the Laurent coefficient `a_{-1}` of `f` about `pole`, i.e.
///
/// ```text
/// Res_{z=p} f = (1 / 2πi) ∮_{|z-p|=r} f(z) dz
/// ```
///
/// for **any** pole order, as long as the circle encloses no other
/// singularity.  The contour integral is evaluated with the periodic
/// trapezoid rule on the circle parametrisation `z = p + r·e^{iθ}`, which is
/// spectrally accurate for analytic periodic integrands (the quadrature error
/// decays geometrically in the node count `N`, like `(r/R)^N` where `R` is the
/// distance to the nearest other singularity), so the result is accurate to
/// near machine precision for the default radius `r = 10⁻³`.
///
/// Use [`residue_with_radius`] when another singularity lies closer to the
/// pole than the default radius.
///
/// # Arguments
/// * `f`     - meromorphic function
/// * `pole`  - location of the pole
/// * `order` - order of the pole (≥ 1).  It is used to size the quadrature
///   (the rule needs at least `order` nodes to annihilate the aliasing of the
///   principal part); `order == 0` means `f` is regular at `pole`, so the
///   residue is exactly zero.
///
/// # Returns
/// Numerical approximation of the residue.
pub fn residue(f: impl Fn(Complex) -> Complex, pole: Complex, order: usize) -> Complex {
    residue_with_radius(f, pole, order, RESIDUE_DEFAULT_RADIUS)
}

/// Compute the residue of `f` at `pole` using a circle of the given `radius`.
///
/// Identical to [`residue`] but lets the caller choose the integration
/// radius, which must be smaller than the distance from `pole` to any other
/// singularity of `f`.  Very small radii amplify cancellation for high-order
/// poles (the integrand grows like `r^{1-order}`), so prefer the largest
/// radius that still isolates the pole.
///
/// # Arguments
/// * `f`      - meromorphic function
/// * `pole`   - location of the pole
/// * `order`  - order of the pole (see [`residue`])
/// * `radius` - radius of the integration circle; must be finite and positive
///
/// # Returns
/// The residue, or a `NaN`-valued complex number if `radius` is not a finite
/// positive number (consistent with [`Complex::inv`] on invalid input).
pub fn residue_with_radius(
    f: impl Fn(Complex) -> Complex,
    pole: Complex,
    order: usize,
    radius: f64,
) -> Complex {
    if order == 0 {
        return Complex::zero();
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Complex::new(f64::NAN, f64::NAN);
    }
    let n_points = RESIDUE_MIN_POINTS.max(4 * order);
    let mut sum = Complex::zero();
    for k in 0..n_points {
        let theta = 2.0 * PI * (k as f64) / (n_points as f64);
        // Node offset z_k - p = r·e^{iθ_k}; dz = i·(z_k - p)·dθ, so
        // (1/2πi) ∮ f dz ≈ (1/N) Σ f(z_k)·(z_k - p).
        let offset = Complex::new(radius * theta.cos(), radius * theta.sin());
        sum = sum + f(pole + offset) * offset;
    }
    sum / (n_points as f64)
}

// ---------------------------------------------------------------------------
// Conformal mappings
// ---------------------------------------------------------------------------

/// Joukowski conformal mapping: w = z + c²/z.
///
/// Maps circles near the unit circle to airfoil-like shapes.
///
/// # Arguments
/// * `z` - input complex coordinate
/// * `c` - parameter (typically near 1.0)
///
/// # Returns
/// Mapped complex coordinate.
pub fn conformal_map_joukowski(z: Complex, c: f64) -> Complex {
    let c2 = Complex::from_real(c * c);
    z + c2 / z
}

/// Möbius (linear fractional) transformation: w = (a·z + b) / (c·z + d).
///
/// Maps circles/lines to circles/lines and preserves angles (conformal).
///
/// # Arguments
/// * `z` - input complex coordinate
/// * `a`, `b`, `c`, `d` - complex coefficients (ad − bc ≠ 0)
///
/// # Returns
/// Mapped complex coordinate.
pub fn conformal_map_mobius(z: Complex, a: Complex, b: Complex, c: Complex, d: Complex) -> Complex {
    (a * z + b) / (c * z + d)
}

// ---------------------------------------------------------------------------
// Contour integral struct
// ---------------------------------------------------------------------------

/// A polygonal contour in the complex plane.
///
/// The path is represented as an ordered list of vertices; straight-line
/// segments connect consecutive vertices.  The contour is treated as closed
/// when performing integrals or computing winding numbers.
#[derive(Debug, Clone)]
pub struct ContourIntegral {
    /// Ordered list of waypoints defining the polygonal path.
    pub path: Vec<Complex>,
}

impl ContourIntegral {
    /// Construct a new contour from a list of waypoints.
    pub fn new(path: Vec<Complex>) -> Self {
        Self { path }
    }

    /// Compute ∮_Γ f(z) dz by the trapezoidal rule.
    ///
    /// The contour is automatically closed.
    pub fn integrate(&self, f: impl Fn(Complex) -> Complex) -> Complex {
        cauchy_integral(f, &self.path)
    }

    /// Compute the winding number of the contour around point `p`.
    ///
    /// The winding number counts (signed) how many times the closed path
    /// winds around `p` in the counter-clockwise direction.
    ///
    /// # Arguments
    /// * `p` - query point
    ///
    /// # Returns
    /// Winding number (positive = counter-clockwise, negative = clockwise).
    pub fn winding_number(&self, p: Complex) -> i32 {
        if self.path.len() < 2 {
            return 0;
        }
        let n = self.path.len();
        let mut winding = 0.0_f64;
        for i in 0..n {
            let z0 = self.path[i] - p;
            let z1 = self.path[(i + 1) % n] - p;
            // Increment angle: arg(z1/z0)
            let dtheta = (z1 / z0).arg();
            winding += dtheta;
        }
        (winding / (2.0 * PI)).round() as i32
    }
}

// ---------------------------------------------------------------------------
// Gamma function (Lanczos approximation)
// ---------------------------------------------------------------------------

// Lanczos coefficients for g=7, n=9 (Numerical Recipes edition)
const LANCZOS_G: f64 = 7.0;
const LANCZOS_C: [f64; 9] = [
    0.999_999_999_999_809_3,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Compute the Gamma function Γ(z) using the Lanczos approximation.
///
/// Valid for Re(z) > 0.  For Re(z) ≤ 0 the reflection formula
/// Γ(z)·Γ(1−z) = π/sin(πz) is used automatically.
///
/// # Arguments
/// * `z` - complex argument
///
/// # Returns
/// Approximation of Γ(z).
pub fn gamma_lanczos(z: Complex) -> Complex {
    if z.re < 0.5 {
        // Reflection formula: Γ(z) = π / (sin(πz) · Γ(1-z))
        let pi_z = Complex::new(PI * z.re, PI * z.im);
        let sin_piz = pi_z.sin();
        let g1mz = gamma_lanczos(Complex::new(1.0 - z.re, -z.im));
        Complex::from_real(PI) / (sin_piz * g1mz)
    } else {
        let z1 = z - Complex::from_real(1.0);
        let mut x = Complex::from_real(LANCZOS_C[0]);
        for (i, &c) in LANCZOS_C[1..].iter().enumerate() {
            x = x + Complex::from_real(c) / (z1 + Complex::from_real((i + 1) as f64));
        }
        let t = z1 + Complex::from_real(LANCZOS_G + 0.5);
        let sqrt2pi = (2.0 * PI).sqrt();
        Complex::from_real(sqrt2pi) * t.cpow(z1 + Complex::from_real(0.5)) * (-t).exp() * x
    }
}

// ---------------------------------------------------------------------------
// Beta function
// ---------------------------------------------------------------------------

/// Compute the Beta function B(a, b) = Γ(a)·Γ(b) / Γ(a+b).
///
/// Valid for a, b > 0.
///
/// # Arguments
/// * `a` - first parameter (> 0)
/// * `b` - second parameter (> 0)
///
/// # Returns
/// Value of B(a, b).
pub fn beta_function(a: f64, b: f64) -> f64 {
    let ga = gamma_lanczos(Complex::from_real(a));
    let gb = gamma_lanczos(Complex::from_real(b));
    let gab = gamma_lanczos(Complex::from_real(a + b));
    (ga * gb / gab).re
}

// ---------------------------------------------------------------------------
// Riemann zeta function (Euler–Maclaurin)
// ---------------------------------------------------------------------------

/// Approximate the Riemann zeta function ζ(s) using the Euler–Maclaurin formula.
///
/// This is an asymptotic expansion and is most accurate for Re(s) > 1.
/// For values near the critical strip, more terms improve accuracy.
///
/// # Arguments
/// * `s`       - complex argument (Re(s) ≠ 1)
/// * `n_terms` - number of partial-sum terms (20–50 recommended)
///
/// # Returns
/// Approximation of ζ(s).
pub fn zeta_euler_maclaurin(s: Complex, n_terms: usize) -> Complex {
    let n = n_terms.max(2);
    // Partial sum: Σ_{k=1}^{n-1} k^{-s}
    let mut sum = Complex::zero();
    for k in 1..n {
        sum = sum + Complex::from_real(k as f64).cpow(-s);
    }
    // Euler–Maclaurin correction: ½ n^{-s} + n^{1-s}/(s-1)
    let nf = Complex::from_real(n as f64);
    let half_term = nf.cpow(-s) / 2.0;
    let integral_term = nf.cpow(Complex::from_real(1.0) - s) / (s - Complex::from_real(1.0));
    // Bernoulli correction (B_2/2! · s · n^{-s-1})
    let b2_term = Complex::from_real(1.0 / 12.0) * s * nf.cpow(-s - Complex::from_real(1.0));
    sum + half_term + integral_term + b2_term
}

// ---------------------------------------------------------------------------
// Complete elliptic integrals (AGM method)
// ---------------------------------------------------------------------------

/// Complete elliptic integral of the first kind K(k).
///
/// Computed via the arithmetic-geometric mean (AGM) iteration.
///
/// K(k) = ∫_0^{π/2} dθ / √(1 − k²·sin²θ)
///
/// # Arguments
/// * `k` - elliptic modulus, 0 ≤ k < 1
///
/// # Returns
/// Value of K(k).
pub fn elliptic_k(k: f64) -> f64 {
    debug_assert!(k.abs() < 1.0, "modulus |k| must be < 1");
    let mut a = 1.0_f64;
    let mut b = (1.0 - k * k).sqrt();
    for _ in 0..50 {
        let a_new = (a + b) / 2.0;
        let b_new = (a * b).sqrt();
        a = a_new;
        b = b_new;
        if (a - b).abs() < 1e-15 {
            break;
        }
    }
    PI / (2.0 * a)
}

/// Complete elliptic integral of the second kind E(k).
///
/// Computed via a variant of the AGM iteration:
///
/// E(k) = K(k) · (1 − Σ_{n≥0} 2^{n−1} c_n²),
///
/// where `a_0 = 1`, `b_0 = √(1 − k²)`, `c_0 = k` and
/// `c_{n+1} = (a_n − b_n) / 2` accompany the AGM sequence of K(k).
///
/// E(k) = ∫_0^{π/2} √(1 − k²·sin²θ) dθ
///
/// # Arguments
/// * `k` - elliptic modulus, 0 ≤ k < 1
///
/// # Returns
/// Value of E(k).
pub fn elliptic_e(k: f64) -> f64 {
    debug_assert!(k.abs() < 1.0, "modulus |k| must be < 1");
    let mut a = 1.0_f64;
    let mut b = (1.0 - k * k).sqrt();
    let mut c = k;
    // Weight 2^{n-1} of the n-th term; the loop doubles it before use, so the
    // first correction (n = 1) carries weight 2^0 = 1.
    let mut weight = 0.5_f64;
    // n = 0 term: 2^{-1} c_0².
    let mut sum = 1.0 - weight * c * c;
    for _ in 0..50 {
        let a_new = (a + b) / 2.0;
        let b_new = (a * b).sqrt();
        c = (a - b) / 2.0;
        a = a_new;
        b = b_new;
        weight *= 2.0;
        sum -= weight * c * c;
        if c.abs() < 1e-15 {
            break;
        }
    }
    sum * PI / (2.0 * a)
}

// ---------------------------------------------------------------------------
// Bessel functions
// ---------------------------------------------------------------------------

/// Euler–Mascheroni constant γ.
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Largest argument for which the power series of J₀, J₁, Y₀ and Y₁ are used.
///
/// The series terms grow like eˣ before they decay while the sum stays of order
/// one, so the cancellation error grows from about 1e-14 at x = 10 to 1e-9 at
/// x = 20.  Hankel's asymptotic expansion, whose smallest term is about e^(−2x),
/// is the more accurate of the two beyond x ≈ 12; both are good to about 1e-12
/// at the switch.
const BESSEL_SERIES_LIMIT: f64 = 12.0;

/// Terms of the power series; for x ≤ 12 they converge below 1e-17 after about 45.
const BESSEL_SERIES_TERMS: usize = 80;

/// J₀(x) from its power series Σ (−1)^m (x/2)^(2m) / (m!)².
fn bessel_j0_series(x: f64) -> f64 {
    let quarter_x2 = 0.25 * x * x;
    let mut term = 1.0_f64;
    let mut sum = 1.0_f64;
    for m in 1..=BESSEL_SERIES_TERMS {
        let m = m as f64;
        term *= -quarter_x2 / (m * m);
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    sum
}

/// J₁(x) from its power series Σ (−1)^m (x/2)^(2m+1) / (m! (m+1)!).
fn bessel_j1_series(x: f64) -> f64 {
    let quarter_x2 = 0.25 * x * x;
    let mut term = 0.5 * x;
    let mut sum = term;
    for m in 1..=BESSEL_SERIES_TERMS {
        let m = m as f64;
        term *= -quarter_x2 / (m * (m + 1.0));
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    sum
}

/// Y₀(x) for x > 0 from
/// `(2/π) [ (ln(x/2) + γ) J₀(x) + Σ_{m≥1} (−1)^(m+1) H_m (x²/4)^m / (m!)² ]`
/// (Abramowitz & Stegun 9.1.13), `H_m` being the harmonic numbers.
fn bessel_y0_series(x: f64) -> f64 {
    let quarter_x2 = 0.25 * x * x;
    // `term` is (−1)^m (x²/4)^m / (m!)², the signed series term of J₀.
    let mut term = 1.0_f64;
    let mut harmonic = 0.0_f64;
    let mut correction = 0.0_f64;
    for m in 1..=BESSEL_SERIES_TERMS {
        let m = m as f64;
        term *= -quarter_x2 / (m * m);
        harmonic += 1.0 / m;
        let contribution = term * harmonic;
        correction -= contribution;
        if contribution.abs() < 1e-17 * correction.abs() {
            break;
        }
    }
    FRAC_2_PI * (((0.5 * x).ln() + EULER_GAMMA) * bessel_j0_series(x) + correction)
}

/// Y₁(x) for x > 0 from
/// `(2/π) ln(x/2) J₁(x) − 2/(π x) − (1/π) Σ_{k≥0} (−1)^k [ψ(k+1) + ψ(k+2)] (x/2)^(2k+1) / (k! (k+1)!)`
/// (Abramowitz & Stegun 9.1.11) with `ψ(k+1) = −γ + H_k`.
fn bessel_y1_series(x: f64) -> f64 {
    let quarter_x2 = 0.25 * x * x;
    // `term` is the signed series term of J₁: (−1)^k (x/2)^(2k+1) / (k! (k+1)!).
    let mut term = 0.5 * x;
    let mut harmonic_k = 0.0_f64;
    let mut harmonic_k1 = 1.0_f64;
    let mut weighted = (harmonic_k + harmonic_k1) * term;
    for k in 1..=BESSEL_SERIES_TERMS {
        let k = k as f64;
        term *= -quarter_x2 / (k * (k + 1.0));
        harmonic_k = harmonic_k1;
        harmonic_k1 += 1.0 / (k + 1.0);
        let contribution = (harmonic_k + harmonic_k1) * term;
        weighted += contribution;
        if contribution.abs() < 1e-17 * weighted.abs() {
            break;
        }
    }
    // −2γ J₁ from the digamma values is folded into the (ln(x/2) + γ) factor.
    FRAC_2_PI * (((0.5 * x).ln() + EULER_GAMMA) * bessel_j1_series(x) - 1.0 / x) - weighted / PI
}

/// `(J_ν(x), Y_ν(x))` for `ν = order ∈ {0, 1}` and large `x > 0` from Hankel's
/// asymptotic expansion (Abramowitz & Stegun 9.2.5–9.2.10):
///
/// `J_ν = sqrt(2/(πx)) (P cos χ − Q sin χ)`, `Y_ν = sqrt(2/(πx)) (P sin χ + Q cos χ)`,
/// `χ = x − (ν/2 + 1/4) π`, with `P = Σ (−1)^j a_{2j}/x^(2j)`,
/// `Q = Σ (−1)^j a_{2j+1}/x^(2j+1)` and
/// `a_k = (4ν² − 1²)(4ν² − 3²)…(4ν² − (2k−1)²) / (k! 8^k)`.
///
/// The series is asymptotic: summation stops at the smallest term, whose size
/// is about e^(−2x).
fn bessel_hankel(order: u8, x: f64) -> (f64, f64) {
    let mu = 4.0 * f64::from(order) * f64::from(order);
    let mut p = 1.0_f64;
    let mut q = 0.0_f64;
    let mut term = 1.0_f64;
    let mut previous = f64::INFINITY;
    for k in 1..=60_u32 {
        let k_f = f64::from(k);
        let odd = 2.0 * k_f - 1.0;
        term *= (mu - odd * odd) / (8.0 * k_f * x);
        let size = term.abs();
        if size >= previous {
            break;
        }
        previous = size;
        match k % 4 {
            1 => q += term,
            2 => p -= term,
            3 => q -= term,
            _ => p += term,
        }
        if size < 1e-17 {
            break;
        }
    }
    let chi = x - (0.5 * f64::from(order) + 0.25) * PI;
    let amplitude = (FRAC_2_PI / x).sqrt();
    let (sin_chi, cos_chi) = chi.sin_cos();
    (
        amplitude * (p * cos_chi - q * sin_chi),
        amplitude * (p * sin_chi + q * cos_chi),
    )
}

/// J₀(x) for any finite real x.
fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax <= BESSEL_SERIES_LIMIT {
        bessel_j0_series(ax)
    } else {
        bessel_hankel(0, ax).0
    }
}

/// J₁(x) for any finite real x.
fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    let magnitude = if ax <= BESSEL_SERIES_LIMIT {
        bessel_j1_series(ax)
    } else {
        bessel_hankel(1, ax).0
    };
    if x < 0.0 { -magnitude } else { magnitude }
}

/// Y₀(x) for finite x > 0.
fn bessel_y0(x: f64) -> f64 {
    if x <= BESSEL_SERIES_LIMIT {
        bessel_y0_series(x)
    } else {
        bessel_hankel(0, x).1
    }
}

/// Y₁(x) for finite x > 0.
fn bessel_y1(x: f64) -> f64 {
    if x <= BESSEL_SERIES_LIMIT {
        bessel_y1_series(x)
    } else {
        bessel_hankel(1, x).1
    }
}

/// J_n(x) for `n ≥ 2` and `x ≥ n` by upward recurrence from J₀ and J₁.
///
/// `J_{k+1} = (2k/x) J_k − J_{k−1}` is stable while `k < x`, where J_k does
/// not decay.
fn bessel_j_upward(order: u32, x: f64) -> f64 {
    let mut previous = bessel_j0(x);
    let mut current = bessel_j1(x);
    for k in 1..order {
        let next = (2.0 * f64::from(k) / x) * current - previous;
        previous = current;
        current = next;
    }
    current
}

/// J_n(x) for `n ≥ 2` and `1 < x < n` by Miller's backward recurrence.
///
/// Upward recurrence is unstable there (J_k decays while the Y_k admixture of
/// the rounding error grows).  Instead the recurrence runs downward from an even
/// index above `n + sqrt(160 n)` starting from arbitrary values and the result
/// is normalised with `J₀ + 2 Σ_{k≥1} J_{2k} = 1`.  Values are rescaled by 1e-10
/// whenever they exceed 1e10.
fn bessel_j_miller(order: u32, x: f64) -> f64 {
    let n = order as usize;
    let start = 2 * ((n + (160.0 * n as f64).sqrt() as usize) / 2);
    let two_over_x = 2.0 / x;
    let mut next = 0.0_f64;
    let mut current = 1.0_f64;
    // J₀ + J₂ + J₄ + … (start is even, so the seed 1.0 at index `start` counts).
    let mut even_sum = 1.0_f64;
    let mut wanted = 0.0_f64;
    for j in (1..=start).rev() {
        let previous = j as f64 * two_over_x * current - next;
        next = current;
        current = previous;
        if current.abs() > 1e10 {
            current *= 1e-10;
            next *= 1e-10;
            wanted *= 1e-10;
            even_sum *= 1e-10;
        }
        if j == n {
            wanted = next;
        }
        if (j - 1).is_multiple_of(2) {
            even_sum += current;
        }
    }
    // `current` is now J₀, counted once in `even_sum`: J₀ + 2 Σ J_{2k} = 2 even_sum − J₀.
    wanted / (2.0 * even_sum - current)
}

/// J_n(x) for `n ≥ 2` and `0 < x ≤ 1` from its power series
/// `Σ_{m≥0} (−1)^m (x/2)^(2m+n) / (m! (m+n)!)`.
///
/// It converges within a few terms here and, unlike the backward recurrence
/// (whose multiplier `2j/x` overflows for tiny `x`), cannot overflow.  The
/// leading term `(x/2)^n / n!` is built factor by factor, so it underflows
/// gracefully to zero.
fn bessel_j_small_argument(order: u32, x: f64) -> f64 {
    let half_x = 0.5 * x;
    let mut term = 1.0_f64;
    for k in 1..=order {
        term *= half_x / f64::from(k);
        if term == 0.0 {
            return 0.0;
        }
    }
    let quarter_x2 = half_x * half_x;
    let mut sum = term;
    for m in 1..=BESSEL_SERIES_TERMS {
        let m = m as f64;
        term *= -quarter_x2 / (m * (m + f64::from(order)));
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    sum
}

/// `true` when `|J_n(x)| ≤ (e x / 2n)^n` is below the smallest subnormal, so the
/// result is exactly zero in `f64` and no recurrence needs to be run.
fn bessel_j_underflows(order: u32, x: f64) -> bool {
    let n = f64::from(order);
    let ratio = std::f64::consts::E * x / (2.0 * n);
    ratio < 1.0 && n * ratio.ln() < -745.5
}

/// Bessel function of the first kind J_n(x) for integer order `n` and real `x`.
///
/// J₀ and J₁ come from their power series for |x| ≤ 12 and from Hankel's
/// asymptotic expansion beyond.  Higher orders follow by upward recurrence when
/// `|x| ≥ |n|` (stable there), by the power series for `|x| ≤ 1`, and by
/// Miller's backward recurrence, normalised with `J₀ + 2 Σ J_{2k} = 1`, for
/// `1 < |x| < |n|` (where the upward recurrence loses all digits).  Accuracy is
/// about 1e-11 relative to the scale of the result; for very large `|x|` it is limited by the argument reduction of
/// `x − (n/2 + 1/4)π`.
///
/// Negative orders and arguments use `J_{−n}(x) = (−1)^n J_n(x)` and
/// `J_n(−x) = (−1)^n J_n(x)`.
///
/// # Arguments
/// * `n` - integer order
/// * `x` - argument (real)
///
/// # Returns
/// Value of J_n(x); `NaN` for a `NaN` argument and `0` for an infinite one.
pub fn bessel_j(n: i32, x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let order = n.unsigned_abs();
    let flips_sign = order % 2 == 1 && (n < 0) != (x < 0.0);
    let sign = if flips_sign { -1.0 } else { 1.0 };
    if x == 0.0 {
        return if order == 0 { 1.0 } else { 0.0 };
    }
    if x.is_infinite() {
        return 0.0;
    }
    let ax = x.abs();
    let magnitude = match order {
        0 => bessel_j0(ax),
        1 => bessel_j1(ax),
        _ if ax >= f64::from(order) => bessel_j_upward(order, ax),
        _ if bessel_j_underflows(order, ax) => 0.0,
        _ if ax <= 1.0 => bessel_j_small_argument(order, ax),
        _ => bessel_j_miller(order, ax),
    };
    sign * magnitude
}

/// Bessel function of the second kind Y_n(x) (Neumann function) for integer
/// order `n` and `x > 0`.
///
/// Y₀ and Y₁ come from their series (Abramowitz & Stegun 9.1.11 and 9.1.13) for
/// x ≤ 12 and from Hankel's asymptotic expansion beyond; higher orders follow
/// from the upward recurrence `Y_{k+1} = (2k/x) Y_k − Y_{k−1}`, which is stable
/// for Y.  Accuracy is about 1e-11 relative to the scale of the result.
/// `Y_{−n}(x) = (−1)^n Y_n(x)`.
///
/// # Arguments
/// * `n` - integer order
/// * `x` - positive real argument (x > 0)
///
/// # Returns
/// Value of Y_n(x).  Y_n is not real for `x < 0`, so a negative or `NaN`
/// argument gives `NaN`; at `x = 0` the function diverges to `−∞`
/// (`+∞` for a negative odd order) and an infinite argument gives `0`.
pub fn bessel_y(n: i32, x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    let order = n.unsigned_abs();
    let sign = if n < 0 && order % 2 == 1 { -1.0 } else { 1.0 };
    if x == 0.0 {
        return -sign * f64::INFINITY;
    }
    if x.is_infinite() {
        return 0.0;
    }
    let y0 = bessel_y0(x);
    if order == 0 {
        return y0;
    }
    let mut previous = y0;
    let mut current = bessel_y1(x);
    for k in 1..order {
        // Once Y_k has overflowed to −∞ the recurrence would produce NaN.
        if !current.is_finite() {
            break;
        }
        let next = (2.0 * f64::from(k) / x) * current - previous;
        previous = current;
        current = next;
    }
    sign * current
}

// ---------------------------------------------------------------------------
// Orthogonal polynomials
// ---------------------------------------------------------------------------

/// Legendre polynomial P_n(x) evaluated at x ∈ \[−1, 1\].
///
/// Uses the three-term recurrence:
/// (n+1) P_{n+1}(x) = (2n+1) x P_n(x) − n P_{n−1}(x)
///
/// # Arguments
/// * `n` - degree (non-negative)
/// * `x` - argument ∈ \[−1, 1\]
///
/// # Returns
/// Value of P_n(x).
pub fn legendre_p(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p_prev = 1.0_f64;
    let mut p_curr = x;
    for k in 1..n {
        let p_next = ((2 * k + 1) as f64 * x * p_curr - k as f64 * p_prev) / ((k + 1) as f64);
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}

/// Hermite polynomial H_n(x) (physicists' convention).
///
/// Uses the recurrence: H_{n+1}(x) = 2x H_n(x) − 2n H_{n−1}(x)
///
/// # Arguments
/// * `n` - degree (non-negative)
/// * `x` - argument
///
/// # Returns
/// Value of H_n(x).
pub fn hermite_h(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return 2.0 * x;
    }
    let mut h_prev = 1.0_f64;
    let mut h_curr = 2.0 * x;
    for k in 1..n {
        let h_next = 2.0 * x * h_curr - 2.0 * k as f64 * h_prev;
        h_prev = h_curr;
        h_curr = h_next;
    }
    h_curr
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;
    const MED_TOL: f64 = 1e-6;

    // --- Complex arithmetic ---

    #[test]
    fn test_complex_add() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a + b;
        assert!((c.re - 4.0).abs() < TOL);
        assert!((c.im - 6.0).abs() < TOL);
    }

    #[test]
    fn test_complex_mul() {
        // (1+2i)(3+4i) = 3+4i+6i+8i² = 3-8 + 10i = -5+10i
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a * b;
        assert!((c.re - (-5.0)).abs() < TOL);
        assert!((c.im - 10.0).abs() < TOL);
    }

    #[test]
    fn test_complex_div() {
        // (1+2i)/(1+i) = (1+2i)(1-i)/2 = (1+2i-i-2i²)/2 = (3+i)/2
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(1.0, 1.0);
        let c = a / b;
        assert!((c.re - 1.5).abs() < TOL);
        assert!((c.im - 0.5).abs() < TOL);
    }

    #[test]
    fn test_complex_abs_arg() {
        let z = Complex::new(0.0, 1.0);
        assert!((z.abs() - 1.0).abs() < TOL);
        assert!((z.arg() - PI / 2.0).abs() < TOL);
    }

    #[test]
    fn test_complex_exp() {
        // e^(iπ) = -1
        let z = Complex::new(0.0, PI);
        let ez = z.exp();
        assert!((ez.re - (-1.0)).abs() < 1e-14);
        assert!(ez.im.abs() < 1e-14);
    }

    #[test]
    fn test_complex_ln() {
        // ln(e) = 1
        let z = Complex::from_real(std::f64::consts::E);
        let lnz = z.ln();
        assert!((lnz.re - 1.0).abs() < TOL);
        assert!(lnz.im.abs() < TOL);
    }

    #[test]
    fn test_complex_sqrt() {
        // √(-1) = i
        let z = Complex::from_real(-1.0);
        let s = z.sqrt();
        assert!(s.re.abs() < 1e-14);
        assert!((s.im - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_complex_pow() {
        // i² = -1
        let z = Complex::new(0.0, 1.0);
        let z2 = z.pow(2.0);
        assert!((z2.re - (-1.0)).abs() < 1e-14);
        assert!(z2.im.abs() < 1e-14);
    }

    #[test]
    fn test_complex_sin_cos() {
        // sin(0) = 0, cos(0) = 1
        let z = Complex::zero();
        assert!(z.sin().abs() < TOL);
        assert!((z.cos().re - 1.0).abs() < TOL);
    }

    #[test]
    fn test_complex_sinh_cosh() {
        // cosh²(z) - sinh²(z) = 1  (complex identity)
        let z = Complex::new(1.0, 0.5);
        let ch = z.cosh();
        let sh = z.sinh();
        let diff = ch * ch - sh * sh;
        assert!((diff.re - 1.0).abs() < 1e-12, "re diff = {}", diff.re);
        assert!(diff.im.abs() < 1e-12, "im = {}", diff.im);
    }

    #[test]
    fn test_complex_conj() {
        let z = Complex::new(3.0, 4.0);
        let zc = z.conj();
        assert!((zc.re - 3.0).abs() < TOL);
        assert!((zc.im - (-4.0)).abs() < TOL);
    }

    #[test]
    fn test_complex_inv() {
        // 1/(1+i) = (1-i)/2
        let z = Complex::new(1.0, 1.0);
        let inv = z.inv();
        assert!((inv.re - 0.5).abs() < TOL);
        assert!((inv.im - (-0.5)).abs() < TOL);
    }

    // --- Cauchy integral ---

    #[test]
    fn test_cauchy_integral_constant() {
        // ∮ 1 dz around any closed contour = 0
        let n = 64_usize;
        let contour: Vec<Complex> = (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64;
                Complex::new(t.cos(), t.sin())
            })
            .collect();
        let result = cauchy_integral(|_z| Complex::one(), &contour);
        assert!(result.abs() < 1e-10);
    }

    #[test]
    fn test_cauchy_integral_z() {
        // ∮ z dz = 0 (analytic function, closed contour)
        let n = 256_usize;
        let contour: Vec<Complex> = (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64;
                Complex::new(t.cos(), t.sin())
            })
            .collect();
        let result = cauchy_integral(|z| z, &contour);
        assert!(result.abs() < 1e-10);
    }

    #[test]
    fn test_cauchy_integral_one_over_z() {
        // ∮_{|z|=1} 1/z dz = 2πi
        let n = 512_usize;
        let contour: Vec<Complex> = (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64;
                Complex::new(t.cos(), t.sin())
            })
            .collect();
        let result = cauchy_integral(|z| z.inv(), &contour);
        assert!(result.re.abs() < 1e-3);
        assert!((result.im - 2.0 * PI).abs() < 1e-2);
    }

    // --- Residue ---

    #[test]
    fn test_residue_simple_pole() {
        // Res_{z=0} 1/z = 1
        let r = residue(|z| z.inv(), Complex::zero(), 1);
        assert!((r.re - 1.0).abs() < 1e-4);
        assert!(r.im.abs() < 1e-4);
    }

    #[test]
    fn test_residue_simple_pole_is_spectrally_accurate() {
        // The circle trapezoid rule is exact for 1/z up to rounding, unlike a
        // chord (polygon) rule whose error is O(1/N²) ≈ 2.5e-5 at N = 512.
        let r = residue(|z| z.inv(), Complex::zero(), 1);
        assert!((r.re - 1.0).abs() < 1e-10, "re = {}", r.re);
        assert!(r.im.abs() < 1e-10, "im = {}", r.im);
    }

    #[test]
    fn test_residue_double_pole_returns_a_minus_one() {
        // f = 1/z² + 3/z has Res = 3; the leading Laurent coefficient a_{-2}
        // (which the old implementation returned) is 1.
        let r = residue(|z| z.inv() * z.inv() + z.inv() * 3.0, Complex::zero(), 2);
        assert!((r.re - 3.0).abs() < 1e-8, "re = {}", r.re);
        assert!(r.im.abs() < 1e-8, "im = {}", r.im);
    }

    #[test]
    fn test_residue_triple_pole_exponential() {
        // e^z / z³ = 1/z³ + 1/z² + 1/(2z) + ...  ⇒  Res = 1/2.
        let r = residue(|z| z.exp() / (z * z * z), Complex::zero(), 3);
        assert!((r.re - 0.5).abs() < 1e-8, "re = {}", r.re);
        assert!(r.im.abs() < 1e-8, "im = {}", r.im);
    }

    #[test]
    fn test_residue_simple_poles_off_origin() {
        // f = 1/(z(z−1)): Res_{z=1} = +1, Res_{z=0} = −1.
        let f = |z: Complex| (z * (z - Complex::one())).inv();
        let at_one = residue(f, Complex::one(), 1);
        let at_zero = residue(f, Complex::zero(), 1);
        assert!((at_one.re - 1.0).abs() < 1e-8, "re = {}", at_one.re);
        assert!(at_one.im.abs() < 1e-8, "im = {}", at_one.im);
        assert!((at_zero.re + 1.0).abs() < 1e-8, "re = {}", at_zero.re);
        assert!(at_zero.im.abs() < 1e-8, "im = {}", at_zero.im);
    }

    #[test]
    fn test_residue_inverse_sine_squared_vanishes() {
        // 1/sin²(z) = 1/z² + 1/3 + z²/15 + ...: a_{-2} = 1 but Res = 0.
        let r = residue(|z| (z.sin() * z.sin()).inv(), Complex::zero(), 2);
        assert!(r.abs() < 1e-8, "|Res| = {}", r.abs());
    }

    #[test]
    fn test_residue_cotangent_at_pi() {
        // cot(z) has simple poles with residue 1 at every multiple of π.
        let r = residue(|z| z.cos() / z.sin(), Complex::from_real(PI), 1);
        assert!((r.re - 1.0).abs() < 1e-8, "re = {}", r.re);
        assert!(r.im.abs() < 1e-8, "im = {}", r.im);
    }

    #[test]
    fn test_residue_at_complex_pole_location() {
        // f = 2/(z−p)² + (5 − 3i)/(z−p) at p = 1+2i  ⇒  Res = 5 − 3i.
        let p = Complex::new(1.0, 2.0);
        let f = move |z: Complex| {
            let d = z - p;
            (d * d).inv() * 2.0 + d.inv() * Complex::new(5.0, -3.0)
        };
        let r = residue(f, p, 2);
        assert!((r.re - 5.0).abs() < 1e-8, "re = {}", r.re);
        assert!((r.im + 3.0).abs() < 1e-8, "im = {}", r.im);
    }

    #[test]
    fn test_residue_independent_of_declared_order_upper_bound() {
        // Res_{z=0} (1/z² + 3/z) is a property of f, not of the declared order.
        let f = |z: Complex| z.inv() * z.inv() + z.inv() * 3.0;
        let tight = residue(f, Complex::zero(), 2);
        let loose = residue(f, Complex::zero(), 6);
        assert!(tight.approx_eq(loose, 1e-8), "{tight} vs {loose}");
        assert!((loose.re - 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_residue_order_zero_is_zero() {
        // A "pole of order 0" is a regular point: the residue vanishes.
        let r = residue(|z| z.exp(), Complex::zero(), 0);
        assert_eq!(r, Complex::zero());
    }

    #[test]
    fn test_residue_with_radius_isolates_close_poles() {
        // f = 1/(z (z − 4e-4)): Res_{z=0} = −2500, Res_{z=4e-4} = +2500.
        // The default radius (1e-3) encloses both poles, so only their sum (0)
        // is seen; a radius below the pole separation isolates z = 0.
        let f = |z: Complex| (z * (z - Complex::from_real(4e-4))).inv();
        let both = residue(f, Complex::zero(), 1);
        assert!(
            both.abs() < 1e-6,
            "enclosing both poles sums to 0, got {both}"
        );
        let isolated = residue_with_radius(f, Complex::zero(), 1, 1e-4);
        assert!((isolated.re + 2500.0).abs() < 1e-6, "re = {}", isolated.re);
        assert!(isolated.im.abs() < 1e-6, "im = {}", isolated.im);
    }

    #[test]
    fn test_residue_with_radius_rejects_invalid_radius() {
        let f = |z: Complex| z.inv();
        for radius in [0.0, -1e-3, f64::NAN, f64::INFINITY] {
            let r = residue_with_radius(f, Complex::zero(), 1, radius);
            assert!(
                r.re.is_nan() && r.im.is_nan(),
                "radius {radius} must give a NaN residue, got {r}"
            );
        }
    }

    // --- Conformal maps ---

    #[test]
    fn test_joukowski_identity_at_two() {
        // Joukowski of z=2 with c=1: w = 2 + 1/2 = 2.5
        let w = conformal_map_joukowski(Complex::from_real(2.0), 1.0);
        assert!((w.re - 2.5).abs() < TOL);
        assert!(w.im.abs() < TOL);
    }

    #[test]
    fn test_mobius_identity() {
        // Identity: a=1,b=0,c=0,d=1 → w = z
        let z = Complex::new(2.0, 3.0);
        let w = conformal_map_mobius(
            z,
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::one(),
        );
        assert!(z.approx_eq(w, TOL));
    }

    // --- Winding number ---

    #[test]
    fn test_winding_number_inside() {
        let n = 64_usize;
        let circle: Vec<Complex> = (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64;
                Complex::new(t.cos(), t.sin())
            })
            .collect();
        let ci = ContourIntegral::new(circle);
        assert_eq!(ci.winding_number(Complex::zero()), 1);
    }

    #[test]
    fn test_winding_number_outside() {
        let n = 64_usize;
        let circle: Vec<Complex> = (0..n)
            .map(|k| {
                let t = 2.0 * PI * k as f64 / n as f64;
                Complex::new(t.cos(), t.sin())
            })
            .collect();
        let ci = ContourIntegral::new(circle);
        // Point far outside the unit circle
        assert_eq!(ci.winding_number(Complex::new(10.0, 0.0)), 0);
    }

    // --- Gamma function ---

    #[test]
    fn test_gamma_factorial() {
        // Γ(n+1) = n!
        for (n, expected) in [
            (1.0_f64, 1.0),
            (2.0, 1.0),
            (3.0, 2.0),
            (4.0, 6.0),
            (5.0, 24.0),
        ] {
            let g = gamma_lanczos(Complex::from_real(n)).re;
            assert!(
                (g - expected).abs() < MED_TOL,
                "Γ({n}) = {g}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_gamma_half() {
        // Γ(1/2) = √π
        let g = gamma_lanczos(Complex::from_real(0.5)).re;
        assert!((g - PI.sqrt()).abs() < MED_TOL);
    }

    // --- Beta function ---

    #[test]
    fn test_beta_symmetric() {
        // B(a,b) = B(b,a)
        let a = 2.5;
        let b = 3.5;
        assert!((beta_function(a, b) - beta_function(b, a)).abs() < MED_TOL);
    }

    #[test]
    fn test_beta_known_value() {
        // B(1,1) = 1
        assert!((beta_function(1.0, 1.0) - 1.0).abs() < MED_TOL);
        // B(2,2) = 1/6
        assert!((beta_function(2.0, 2.0) - 1.0 / 6.0).abs() < MED_TOL);
    }

    // --- Zeta function ---

    #[test]
    fn test_zeta_known_value() {
        // ζ(2) = π²/6 ≈ 1.6449340668...
        let z = zeta_euler_maclaurin(Complex::from_real(2.0), 50);
        assert!((z.re - PI * PI / 6.0).abs() < 1e-4);
    }

    // --- Elliptic integrals ---

    #[test]
    fn test_elliptic_k_zero() {
        // K(0) = π/2
        assert!((elliptic_k(0.0) - PI / 2.0).abs() < TOL);
    }

    #[test]
    fn test_elliptic_e_zero() {
        // E(0) = π/2
        assert!((elliptic_e(0.0) - PI / 2.0).abs() < TOL);
    }

    #[test]
    fn test_elliptic_ke_relation() {
        // For any k, E(k) ≤ K(k)
        let k = 0.7;
        assert!(elliptic_e(k) <= elliptic_k(k) + 1e-10);
    }

    #[test]
    fn test_elliptic_e_reference_values() {
        // Reference values (mpmath / SciPy `ellipe(k²)`), 1e-12 absolute.
        let cases = [
            (0.5, 1.4674622093394272),
            (std::f64::consts::FRAC_1_SQRT_2, 1.3506438810476755),
            (0.9, 1.1716970527816142),
            (0.99, 1.0284758090288042),
        ];
        for (k, expected) in cases {
            let e = elliptic_e(k);
            assert!(
                (e - expected).abs() < 1e-12,
                "E({k}) = {e}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_elliptic_e_legendre_relation() {
        // Legendre's relation E·K' + E'·K − K·K' = π/2.  At k = k' = 1/√2 it
        // reduces to 2·E·K − K² = π/2, independent of tabulated constants.
        let k = std::f64::consts::FRAC_1_SQRT_2;
        let kk = elliptic_k(k);
        let ee = elliptic_e(k);
        let lhs = 2.0 * ee * kk - kk * kk;
        assert!((lhs - PI / 2.0).abs() < 1e-12, "Legendre residual {lhs}");
    }

    #[test]
    fn test_elliptic_e_legendre_relation_generic_modulus() {
        // Same identity for k ≠ k': E(k)·K(k') + E(k')·K(k) − K(k)·K(k') = π/2.
        for &k in &[0.2_f64, 0.6, 0.85] {
            let kp = (1.0 - k * k).sqrt();
            let lhs = elliptic_e(k) * elliptic_k(kp) + elliptic_e(kp) * elliptic_k(k)
                - elliptic_k(k) * elliptic_k(kp);
            assert!(
                (lhs - PI / 2.0).abs() < 1e-11,
                "Legendre residual {lhs} at k = {k}"
            );
        }
    }

    #[test]
    fn test_elliptic_e_monotone_decreasing_in_modulus() {
        // E(k) decreases from π/2 at k = 0 to 1 as k → 1.
        let mut prev = elliptic_e(0.0);
        for i in 1..=19 {
            let k = f64::from(i) * 0.05;
            let e = elliptic_e(k);
            assert!(e < prev, "E must decrease: E({k}) = {e}, previous {prev}");
            assert!(e > 1.0, "E(k) > 1 for k < 1, got {e} at k = {k}");
            prev = e;
        }
    }

    // --- Bessel functions ---

    #[test]
    fn test_bessel_j0_zero() {
        // J_0(0) = 1
        assert!((bessel_j(0, 0.0) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_bessel_j1_zero() {
        // J_1(0) = 0
        assert!(bessel_j(1, 0.0).abs() < TOL);
    }

    #[test]
    fn test_bessel_j0_known() {
        // J_0(2.4048) ≈ 0  (first zero of J_0 ≈ 2.40483)
        let v = bessel_j(0, 2.40483).abs();
        assert!(v < 1e-4);
    }

    #[test]
    fn test_bessel_recurrence() {
        // J_{n-1}(x) + J_{n+1}(x) = (2n/x) J_n(x)
        let x = 3.0;
        let n = 2;
        let lhs = bessel_j(n - 1, x) + bessel_j(n + 1, x);
        let rhs = (2.0 * n as f64 / x) * bessel_j(n, x);
        assert!((lhs - rhs).abs() < 1e-10);
    }

    // --- Bessel functions: reference values (mpmath, 50 digits) ---

    /// `(n, x, J_n(x), Y_n(x))`, covering the power series (x ≤ 12), Hankel's
    /// expansion (x > 12), upward recurrence (x ≥ n) and Miller's algorithm (x < n).
    const BESSEL_REFERENCE: [(i32, f64, f64, f64); 28] = [
        (0, 0.1, 0.99750156206604, -1.5342386513503667),
        (0, 1.0, 0.7651976865579666, 0.08825696421567696),
        (0, 5.0, -0.1775967713143383, -0.30851762524903376),
        (0, 10.0, -0.24593576445134835, 0.055671167283599395),
        (0, 20.0, 0.16702466434058316, 0.06264059680938383),
        (0, 50.0, 0.055812327669251816, -0.09806499547007708),
        (1, 0.1, 0.049937526036242, -6.4589510947020266),
        (1, 1.0, 0.4400505857449335, -0.7812128213002887),
        (1, 5.0, -0.32757913759146523, 0.14786314339122683),
        (1, 10.0, 0.04347274616886144, 0.24901542420695388),
        (1, 20.0, 0.06683312417585005, -0.1655116143625213),
        (1, 50.0, -0.09751182812517514, -0.05679566856201477),
        (2, 1.0, 0.11490348493190047, -1.6506826068162543),
        (2, 5.0, 0.046565116277752214, 0.36766288260552454),
        (2, 13.0, -0.21774426424195678, 0.04588764784776922),
        (3, 1.0, 0.019563353982668407, -5.821517605964729),
        (3, 0.001, 2.0833332031250032e-11, -5092958815.560503),
        (5, 1.0, 0.00024975773021123444, -260.4058666258122),
        (5, 12.0, -0.07347096310165858, -0.22981794662508243),
        (7, 13.0, -0.24057094958616052, 0.0015698795407253514),
        (10, 1.0, 2.6306151236874534e-10, -121618014.27868919),
        (10, 5.0, 0.0014678026473104741, -25.1291100956101),
        (10, 12.0, 0.3004760352712693, -0.0228763140704997),
        (20, 3.0, 1.2275946737992987e-15, -13113540041757.447),
        (25, 40.0, -0.026360341175918508, 0.14026971952776407),
        (8, 25.0, 0.1530061666573989, 0.05879476686163041),
        (30, 12.0, 2.552259043034417e-10, -45366214.38603198),
        (50, 0.5, 2.590558066078543e-95, -2.4575848224461087e92),
    ];

    /// `|got − want| ≤ 1e-9 |want|`, plus an absolute floor of 1e-13 where the
    /// function oscillates (x ≥ n) and the cancellation error is absolute.
    fn assert_bessel_close(what: &str, n: i32, x: f64, got: f64, want: f64) {
        let floor = if x >= f64::from(n) { 1e-13 } else { 0.0 };
        assert!(
            (got - want).abs() <= 1e-9 * want.abs() + floor,
            "{what}({n}, {x}) = {got:e}, expected {want:e}"
        );
    }

    #[test]
    fn test_bessel_matches_reference_values() {
        for &(n, x, j, y) in &BESSEL_REFERENCE {
            assert_bessel_close("J", n, x, bessel_j(n, x), j);
            assert_bessel_close("Y", n, x, bessel_y(n, x), y);
        }
    }

    #[test]
    fn test_bessel_y0_y1_y2_at_one() {
        // The old series gave Y0(1) = -0.2012 and Y1(1) = -0.5653 (wrong sign of the
        // correction sum and -1/(πx) instead of -2/(πx)).
        assert!((bessel_y(0, 1.0) - 0.0882569642).abs() < 1e-9);
        assert!((bessel_y(1, 1.0) + 0.7812128213).abs() < 1e-9);
        assert!((bessel_y(2, 1.0) + 1.6506826068).abs() < 1e-9);
        assert!((bessel_y(0, 1.0) - 0.088_256_964_215_676_96).abs() < 1e-14);
        assert!((bessel_y(1, 1.0) + 0.781_212_821_300_288_7).abs() < 1e-14);
    }

    #[test]
    fn test_bessel_is_continuous_across_the_series_switch() {
        // mpmath at 12 ∓ 1e-6: the series side and the Hankel side must agree with
        // the exact values to ~1e-11, not merely with each other.
        let below = (
            11.999999,
            0.04768908734969589,
            -0.22344717080042226,
            -0.22523736973346944,
            -0.057098997781814535,
        );
        let above = (
            12.000001,
            0.047689534243904874,
            -0.2234470381806166,
            -0.22523725553503293,
            -0.057099438739903426,
        );
        for (x, j0, j1, y0, y1) in [below, above] {
            assert!((bessel_j(0, x) - j0).abs() < 1e-11, "J0({x})");
            assert!((bessel_j(1, x) - j1).abs() < 1e-11, "J1({x})");
            assert!((bessel_y(0, x) - y0).abs() < 1e-11, "Y0({x})");
            assert!((bessel_y(1, x) - y1).abs() < 1e-11, "Y1({x})");
        }
    }

    #[test]
    fn test_bessel_wronskian() {
        // J_{ν+1}(x) Y_ν(x) − J_ν(x) Y_{ν+1}(x) = 2/(π x) pins J and Y against each
        // other on both sides of the series/Hankel switch.
        for &x in &[0.3, 1.0, 4.0, 8.0, 11.9, 12.1, 20.0, 30.0, 80.0] {
            let expected = 2.0 / (PI * x);
            let tol = 1e-10 * expected.max(1.0);
            for nu in [0, 1] {
                let w =
                    bessel_j(nu + 1, x) * bessel_y(nu, x) - bessel_j(nu, x) * bessel_y(nu + 1, x);
                assert!(
                    (w - expected).abs() < tol,
                    "ν = {nu}, x = {x}: {w} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn test_bessel_recurrence_in_the_miller_regime() {
        // J_{n−1} + J_{n+1} = (2n/x) J_n with every term in the n > x regime, where
        // the upward recurrence lost all digits (J_10(1) ≈ 2.6e-10 against a 1e-8
        // rounding amplification).
        for &(n, x) in &[(5, 2.0), (9, 3.5), (15, 6.0), (30, 11.0)] {
            let lhs = bessel_j(n - 1, x) + bessel_j(n + 1, x);
            let rhs = 2.0 * f64::from(n) / x * bessel_j(n, x);
            assert!(
                (lhs - rhs).abs() <= 1e-10 * rhs.abs(),
                "n = {n}, x = {x}: {lhs:e} vs {rhs:e}"
            );
        }
    }

    #[test]
    fn test_bessel_reflection_symmetries() {
        for &x in &[0.7, 3.0, 15.0] {
            for n in 0..6 {
                let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
                let j = bessel_j(n, x);
                assert_eq!(bessel_j(-n, x), sign * j, "J_-{n}({x})");
                assert_eq!(bessel_j(n, -x), sign * j, "J_{n}(-{x})");
                assert_eq!(bessel_j(-n, -x), j, "J_-{n}(-{x})");
                if n > 0 {
                    assert_eq!(bessel_y(-n, x), sign * bessel_y(n, x), "Y_-{n}({x})");
                }
            }
        }
    }

    #[test]
    fn test_bessel_special_arguments_and_extreme_orders() {
        assert_eq!(bessel_j(0, 0.0), 1.0);
        assert_eq!(bessel_j(5, 0.0), 0.0);
        assert_eq!(bessel_j(0, f64::INFINITY), 0.0);
        assert!(bessel_j(2, f64::NAN).is_nan());
        // Y diverges to -inf at 0 (+inf for a negative odd order), is not real
        // for x < 0 and vanishes at infinity.
        assert_eq!(bessel_y(0, 0.0), f64::NEG_INFINITY);
        assert_eq!(bessel_y(1, 0.0), f64::NEG_INFINITY);
        assert_eq!(bessel_y(-1, 0.0), f64::INFINITY);
        assert!(bessel_y(0, -1.0).is_nan());
        assert!(bessel_y(0, f64::NAN).is_nan());
        assert_eq!(bessel_y(2, f64::INFINITY), 0.0);
        // The recurrence overflows to -inf instead of producing NaN.
        assert_eq!(bessel_y(500, 1e-3), f64::NEG_INFINITY);
        // Orders whose J underflows return 0 immediately, even i32::MIN (whose
        // negation overflows) and i32::MAX.
        assert_eq!(bessel_j(i32::MIN, 5.0), 0.0);
        assert_eq!(bessel_j(i32::MAX, 1.0), 0.0);
        assert_eq!(bessel_j(2000, 3.0), 0.0);
        // Tiny arguments: the backward recurrence would overflow (its multiplier
        // 2j/x is ~1e150) and return NaN; the series gives (x/2)^n / n!.
        let tiny = bessel_j(2, 1e-150);
        assert!(
            (tiny - 1.25e-301).abs() < 1e-9 * 1.25e-301,
            "J_2(1e-150) = {tiny:e}"
        );
        let tinier = bessel_j(3, 1e-100);
        assert!(
            (tinier - 2.0833333333333334e-302).abs() < 1e-9 * 2.0833333333333334e-302,
            "J_3(1e-100) = {tinier:e}"
        );
        assert_eq!(bessel_j(2, 1e-200), 0.0);
        assert!(bessel_j(2, -1e-150) > 0.0, "J_2 is even in x");
        assert!(bessel_j(3, -1e-100) < 0.0, "J_3 is odd in x");
        // Just above the underflow boundary the value is still computed
        // (mpmath: J_150(5) = 8.243842322e-204).
        let j150 = bessel_j(150, 5.0);
        assert!(
            (j150 - 8.243842322e-204).abs() < 1e-9 * 8.243842322e-204,
            "J_150(5) = {j150:e}"
        );
    }

    // --- Legendre polynomials ---

    #[test]
    fn test_legendre_p0_p1() {
        assert!((legendre_p(0, 0.5) - 1.0).abs() < TOL);
        assert!((legendre_p(1, 0.5) - 0.5).abs() < TOL);
    }

    #[test]
    fn test_legendre_p2() {
        // P_2(x) = (3x² - 1)/2
        let x = 0.7;
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((legendre_p(2, x) - expected).abs() < TOL);
    }

    #[test]
    fn test_legendre_orthogonality_approx() {
        // ∫_{-1}^{1} P_2(x) P_0(x) dx = 0 (orthogonality)
        let n = 200_usize;
        let mut sum = 0.0_f64;
        for i in 0..n {
            let x = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
            sum += legendre_p(2, x) * legendre_p(0, x) * (2.0 / n as f64);
        }
        assert!(sum.abs() < 1e-4);
    }

    // --- Hermite polynomials ---

    #[test]
    fn test_hermite_h0_h1() {
        assert!((hermite_h(0, 1.0) - 1.0).abs() < TOL);
        assert!((hermite_h(1, 1.0) - 2.0).abs() < TOL);
    }

    #[test]
    fn test_hermite_h2() {
        // H_2(x) = 4x² - 2
        let x = 1.5;
        let expected = 4.0 * x * x - 2.0;
        assert!((hermite_h(2, x) - expected).abs() < TOL);
    }

    #[test]
    fn test_hermite_h3() {
        // H_3(x) = 8x³ - 12x
        let x = 2.0;
        let expected = 8.0 * x * x * x - 12.0 * x;
        assert!((hermite_h(3, x) - expected).abs() < TOL);
    }
}
