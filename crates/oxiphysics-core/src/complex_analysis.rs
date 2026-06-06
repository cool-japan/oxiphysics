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

use std::f64::consts::PI;

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

/// Compute the residue of `f` at a (simple or higher-order) pole.
///
/// For a pole of order `order` at `pole`, the residue is estimated
/// numerically by computing the Laurent coefficient via a small circle
/// contour of radius ε.  For a simple pole this equals
/// lim_{z→pole} (z − pole) · f(z).
///
/// # Arguments
/// * `f`     - meromorphic function
/// * `pole`  - location of the pole
/// * `order` - order of the pole (≥ 1)
///
/// # Returns
/// Numerical approximation of the residue.
pub fn residue(f: impl Fn(Complex) -> Complex, pole: Complex, order: usize) -> Complex {
    let eps = 1e-5_f64;
    let n_points = 512_usize;
    // ∮_{|z-p|=ε} f(z) / (z-p)^(1-order) dz / (2πi)  ← standard formula
    // For a pole of order m: Res = 1/(m-1)! · d^(m-1)/dz^(m-1) [(z-p)^m f(z)] at z=p
    // We use the numerical contour formula:
    //   Res = 1/(2πi) ∮ f(z) (z-p)^(order-1) dz
    // which equals 1/(2πi) ∮ g(z) dz  where g(z) = f(z)·(z-p)^(order-1)
    let contour: Vec<Complex> = (0..n_points)
        .map(|k| {
            let theta = 2.0 * PI * (k as f64) / (n_points as f64);
            pole + Complex::new(eps * theta.cos(), eps * theta.sin())
        })
        .collect();
    let integral = cauchy_integral(
        |z| {
            let zp = z - pole;
            let factor = zp.pow((order as f64) - 1.0);
            f(z) * factor
        },
        &contour,
    );
    // divide by 2πi
    integral / Complex::new(0.0, 2.0 * PI)
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
/// Computed via a variant of the AGM iteration.
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
    let mut two_n = 1.0_f64;
    let mut sum = 1.0 - 0.5 * c * c;
    for _ in 0..50 {
        let a_new = (a + b) / 2.0;
        let b_new = (a * b).sqrt();
        c = (a - b) / 2.0;
        a = a_new;
        b = b_new;
        two_n *= 2.0;
        sum -= two_n * c * c;
        if c.abs() < 1e-15 {
            break;
        }
    }
    sum * PI / (2.0 * a)
}

// ---------------------------------------------------------------------------
// Bessel functions
// ---------------------------------------------------------------------------

/// Bessel function of the first kind J_n(x).
///
/// Computed via the power-series expansion for small x, and a backward
/// recurrence (Miller's algorithm) for integer order n ≥ 0.
///
/// # Arguments
/// * `n` - integer order
/// * `x` - argument (real)
///
/// # Returns
/// Value of J_n(x).
pub fn bessel_j(n: i32, x: f64) -> f64 {
    if x == 0.0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    if n < 0 {
        // J_{-n}(x) = (-1)^n J_n(x)
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_j(-n, x);
    }
    // For |x| <= 20, use upward recurrence starting from J_0 and J_1
    let j0 = bessel_j0(x);
    let j1 = bessel_j1(x);
    if n == 0 {
        return j0;
    }
    if n == 1 {
        return j1;
    }
    // Upward recurrence: J_{n+1}(x) = (2n/x) J_n(x) - J_{n-1}(x)
    let mut jm1 = j0;
    let mut j_curr = j1;
    for k in 1..n {
        let jp1 = (2.0 * (k as f64) / x) * j_curr - jm1;
        jm1 = j_curr;
        j_curr = jp1;
    }
    j_curr
}

/// Bessel J_0(x) via its own series / Chebyshev-like formula.
fn bessel_j0(x: f64) -> f64 {
    // Series: J_0(x) = Σ_{m=0}^∞ (-1)^m (x/2)^{2m} / (m!)^2
    let x2 = x * x;
    let mut term = 1.0_f64;
    let mut sum = 1.0_f64;
    for m in 1..=40_usize {
        term *= -x2 / (4.0 * (m as f64) * (m as f64));
        sum += term;
        if term.abs() < 1e-16 * sum.abs() {
            break;
        }
    }
    sum
}

/// Bessel J_1(x).
fn bessel_j1(x: f64) -> f64 {
    // Series: J_1(x) = Σ_{m=0}^∞ (-1)^m (x/2)^{2m+1} / (m! (m+1)!)
    let x2 = x * x;
    let mut term = x / 2.0;
    let mut sum = term;
    for m in 1..=40_usize {
        term *= -x2 / (4.0 * (m as f64) * ((m + 1) as f64));
        sum += term;
        if term.abs() < 1e-16 * sum.abs() {
            break;
        }
    }
    sum
}

/// Bessel function of the second kind Y_n(x) (Neumann function).
///
/// Computed via the standard Wronskian-based recurrence from Y_0 and Y_1.
///
/// # Arguments
/// * `n` - non-negative integer order
/// * `x` - positive real argument (x > 0)
///
/// # Returns
/// Value of Y_n(x).
pub fn bessel_y(n: i32, x: f64) -> f64 {
    debug_assert!(x > 0.0, "Y_n(x) requires x > 0");
    if n < 0 {
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_y(-n, x);
    }
    let y0 = bessel_y0(x);
    let y1 = bessel_y1(x);
    if n == 0 {
        return y0;
    }
    if n == 1 {
        return y1;
    }
    let mut ym1 = y0;
    let mut y_curr = y1;
    for k in 1..n {
        let yp1 = (2.0 * (k as f64) / x) * y_curr - ym1;
        ym1 = y_curr;
        y_curr = yp1;
    }
    y_curr
}

/// Bessel Y_0(x).
fn bessel_y0(x: f64) -> f64 {
    // Y_0(x) = (2/π) [ln(x/2) + γ] J_0(x) + correction series
    // Using the standard expansion
    let gamma_euler = 0.577_215_664_901_532_9;
    let j0 = bessel_j0(x);
    let x2 = x * x;
    // Series for the correction: Σ_{m=1}^∞ (-1)^{m+1} H_m (x/2)^{2m} / (m!)^2
    let mut term = 1.0_f64;
    let mut h = 0.0_f64;
    let mut correction = 0.0_f64;
    for m in 1..=40_usize {
        term *= -x2 / (4.0 * (m as f64) * (m as f64));
        h += 1.0 / (m as f64);
        correction += term * h;
        if (term * h).abs() < 1e-16 {
            break;
        }
    }
    (2.0 / PI) * ((x / 2.0).ln() + gamma_euler) * j0 + (2.0 / PI) * correction
}

/// Bessel Y_1(x).
fn bessel_y1(x: f64) -> f64 {
    let gamma_euler = 0.577_215_664_901_532_9;
    let j1 = bessel_j1(x);
    let x2 = x * x;
    // Series correction for Y_1(x)
    let mut t2 = x / 2.0;
    let mut h = 1.0_f64;
    let mut correction = 0.0_f64;
    for m in 1..=40_usize {
        t2 *= -x2 / (4.0 * (m as f64) * ((m + 1) as f64));
        h += 1.0 / ((m + 1) as f64);
        correction += t2 * (h + h - 1.0 / (m as f64));
        if (t2 * h).abs() < 1e-16 * correction.abs().max(1e-30) {
            break;
        }
    }
    (2.0 / PI) * ((x / 2.0).ln() + gamma_euler - 0.5) * j1 - (1.0 / (PI * x))
        + (2.0 / PI) * correction
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
