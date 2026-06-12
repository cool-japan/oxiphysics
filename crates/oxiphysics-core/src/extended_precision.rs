// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Software double-double (`~31` decimal digit) extended precision arithmetic.
//!
//! This module provides a [`Dd`] type representing a number as an unevaluated
//! sum of two non-overlapping IEEE-754 binary64 components, `hi + lo`, where
//! `hi` carries the leading rounded value and `lo` carries the rounding error.
//! Together they deliver roughly 106 bits (about 31 decimal digits) of
//! significand precision while using only ordinary `f64` hardware operations.
//!
//! # Algorithm
//!
//! The construction rests on *error-free transforms* (EFTs): primitives such as
//! [`two_sum`], [`fast_two_sum`], and [`two_prod`] that compute both the rounded
//! result of a floating-point operation *and* its exact rounding error, so that
//! `hi + lo` equals the true mathematical result with no loss. Higher-level
//! arithmetic ([`Dd`] addition, multiplication, division, and square root) is
//! then composed from these EFTs following the canonical quad-double recipes.
//!
//! All transforms are implemented purely in terms of addition, subtraction, and
//! multiplication via Dekker's splitting trick (split constant
//! `2^27 + 1 = 134217729.0`); **no fused multiply-add (FMA) is used anywhere**,
//! guaranteeing bit-for-bit deterministic results across platforms regardless of
//! whether hardware FMA is available.
//!
//! # References
//!
//! * Y. Hida, X. S. Li, and D. H. Bailey, *Library for Double-Double and
//!   Quad-Double Arithmetic* (the QD library), 2007 — the algorithms for
//!   double-double `add`, `mul`, `div`, and Karp's `sqrt`.
//! * T. J. Dekker, *A floating-point technique for extending the available
//!   precision*, Numerische Mathematik 18 (1971), 224–242 — the splitting and
//!   two-product transforms.
//! * D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical
//!   Algorithms*, 3rd ed. — Theorem B (the two-sum transform).

/// Knuth's `TwoSum` error-free transform of a floating-point addition.
///
/// Returns `(hi, lo)` where `hi = fl(a + b)` is the correctly-rounded sum and
/// `lo` is the exact rounding error, so that `hi + lo == a + b` holds exactly in
/// real arithmetic. This costs 6 floating-point operations and, unlike
/// [`fast_two_sum`], imposes no ordering requirement on the magnitudes of the
/// operands.
///
/// # Algorithm
///
/// This is the branch-free transform of Knuth (TAOCP Vol. 2, Theorem B): after
/// forming the rounded sum `s`, the error is recovered by reconstructing both
/// operands from `s` and summing their individual deviations.
///
/// # Examples
///
/// ```
/// use oxiphysics_core::extended_precision::two_sum;
/// let (hi, lo) = two_sum(1e16, 1.0);
/// assert_eq!(hi, 1e16); // the +1.0 is lost from the leading term
/// assert_eq!(lo, 1.0); // ... but recovered exactly in the error term
/// ```
pub fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

/// Dekker's `FastTwoSum` error-free transform of a floating-point addition.
///
/// Returns `(hi, lo)` where `hi = fl(a + b)` and `lo` is the exact rounding
/// error, so that `hi + lo == a + b` holds exactly. This costs only 3
/// floating-point operations, making it the workhorse for renormalizing the
/// components of a [`Dd`].
///
/// # Precondition
///
/// The result is correct **only when `|a| >= |b|`** (equivalently, when the
/// binary exponent of `a` is at least that of `b`). If this ordering does not
/// hold the returned error term is not guaranteed to be exact; use [`two_sum`]
/// instead when the relative magnitudes are unknown.
///
/// # Algorithm
///
/// Because `|a| >= |b|`, the quantity `s - a` is computed without rounding error,
/// so `b - (s - a)` recovers the exact discarded low-order bits (Dekker, 1971).
pub fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let err = b - (s - a);
    (s, err)
}

/// Dekker's splitting of a `f64` into two non-overlapping `~26`-bit halves.
///
/// Returns `(hi, lo)` with `a == hi + lo` exactly, where `hi` holds the leading
/// 26 significand bits and `lo` the trailing 26 bits. Splitting each factor this
/// way lets a product of two `f64`s be evaluated exactly without hardware FMA,
/// which is what [`two_prod`] relies upon.
///
/// # Algorithm
///
/// Multiplying by the magic constant `2^27 + 1 = 134217729.0` and subtracting in
/// the right order isolates the high half via a single controlled rounding
/// (Dekker, 1971). No FMA is involved.
///
/// # Examples
///
/// ```
/// use oxiphysics_core::extended_precision::split;
/// let (hi, lo) = split(1.0 + 1e-10);
/// assert_eq!(hi + lo, 1.0 + 1e-10);
/// ```
pub fn split(a: f64) -> (f64, f64) {
    const SPLITTER: f64 = 134217729.0;
    let c = SPLITTER * a;
    let hi = c - (c - a);
    let lo = a - hi;
    (hi, lo)
}

/// Dekker's `TwoProduct` error-free transform of a floating-point multiplication.
///
/// Returns `(hi, lo)` where `hi = fl(a * b)` is the correctly-rounded product and
/// `lo` is the exact rounding error, so that `hi + lo == a * b` holds exactly in
/// real arithmetic.
///
/// # Algorithm
///
/// This implementation splits each operand with [`split`] and assembles the
/// exact partial products (Dekker, 1971); it deliberately **does not** use a
/// fused multiply-add, preserving cross-platform determinism.
///
/// # Examples
///
/// ```
/// use oxiphysics_core::extended_precision::two_prod;
/// let a = 1.0 + 1e-10;
/// let (hi, lo) = two_prod(a, a);
/// assert_eq!(hi, a * a);
/// assert!(lo != 0.0); // the squaring is not exactly representable
/// ```
pub fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let err = ((ahi * bhi - p) + ahi * blo + alo * bhi) + alo * blo;
    (p, err)
}

/// `TwoDiff` error-free transform of a floating-point subtraction.
///
/// Returns `(hi, lo)` where `hi = fl(a - b)` and `lo` is the exact rounding
/// error, so that `hi + lo == a - b` holds exactly. This is the sign-adjusted
/// counterpart of [`two_sum`] and rounds out the EFT family for completeness.
fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let s = a - b;
    let bb = s - a;
    let err = (a - (s - bb)) - (b + bb);
    (s, err)
}

/// A double-double extended-precision number: an unevaluated sum `hi + lo`.
///
/// The pair `(hi, lo)` represents a single real value with roughly 106 bits of
/// significand precision (about 31 decimal digits). After normalization the two
/// components are non-overlapping, satisfying the invariant
/// `|lo| <= 0.5 * ulp(hi)`.
///
/// Equality and ordering are defined componentwise over `(hi, lo)`; because the
/// components are IEEE-754 floats (admitting NaN and signed zero), this type
/// implements [`PartialEq`] and [`PartialOrd`] but deliberately **not** `Eq` or
/// `Ord`.
#[derive(Clone, Copy, Debug)]
pub struct Dd {
    /// Leading (correctly-rounded) component carrying the most significant bits.
    pub hi: f64,
    /// Trailing (error) component; `|lo| <= 0.5 * ulp(hi)` after normalization.
    pub lo: f64,
}

impl Dd {
    /// The additive identity `0`, with both components zero.
    pub const ZERO: Dd = Dd { hi: 0.0, lo: 0.0 };
    /// The multiplicative identity `1`, with leading component one.
    pub const ONE: Dd = Dd { hi: 1.0, lo: 0.0 };

    /// Constructs a normalized [`Dd`] from a leading and trailing component.
    ///
    /// The arguments are canonicalized with [`two_sum`] so that the resulting
    /// components are non-overlapping and satisfy the `|lo| <= 0.5 * ulp(hi)`
    /// invariant even if the inputs overlap.
    ///
    /// # Examples
    ///
    /// ```
    /// use oxiphysics_core::extended_precision::Dd;
    /// let x = Dd::new(1.0, 1e-20);
    /// assert_eq!(x.hi, 1.0);
    /// ```
    pub fn new(hi: f64, lo: f64) -> Self {
        let (hi, lo) = two_sum(hi, lo);
        Dd { hi, lo }
    }

    /// Lifts an ordinary `f64` into a [`Dd`] with a zero error term.
    ///
    /// The conversion is exact: every `f64` is representable as a `Dd`.
    pub fn from_f64(x: f64) -> Self {
        Dd { hi: x, lo: 0.0 }
    }

    /// Rounds this [`Dd`] to the nearest `f64`, returning the leading component.
    ///
    /// For a normalized value `hi` is already the correctly-rounded `f64`
    /// approximation of `hi + lo`.
    pub fn to_f64(&self) -> f64 {
        self.hi
    }

    /// Returns the absolute value of this [`Dd`].
    ///
    /// Negation is exact, so the result is the exact magnitude of the input.
    pub fn abs(self) -> Self {
        if self.hi < 0.0 { -self } else { self }
    }

    /// Returns `true` if both components are finite (neither infinite nor NaN).
    pub fn is_finite(&self) -> bool {
        self.hi.is_finite() && self.lo.is_finite()
    }

    /// Computes the square root using Karp's single-Newton-step method.
    ///
    /// # Algorithm
    ///
    /// Karp's trick starts from the reciprocal square root `x = 1/sqrt(hi)`
    /// computed in `f64`, forms the approximate root `ax = hi * x`, and performs
    /// one Newton correction in double-double precision to recover the full
    /// `~31` digits (Hida–Li–Bailey QD library). Only a single double-double
    /// multiplication is needed for the correction.
    ///
    /// # Precondition
    ///
    /// For negative `hi` the result is NaN; for `hi == 0.0` the result is
    /// [`Dd::ZERO`]. The function never panics.
    ///
    /// # Examples
    ///
    /// ```
    /// use oxiphysics_core::extended_precision::Dd;
    /// let two = Dd::from_f64(2.0);
    /// let root = two.sqrt();
    /// assert!(((root * root) - two).to_f64().abs() < 1e-30);
    /// ```
    pub fn sqrt(self) -> Dd {
        if self.hi == 0.0 {
            return Dd::ZERO;
        }
        if self.hi < 0.0 {
            return Dd { hi: f64::NAN, lo: f64::NAN };
        }
        let x = 1.0 / self.hi.sqrt();
        let ax = self.hi * x;
        let diff = self - Dd::from_f64(ax) * Dd::from_f64(ax);
        Dd::from_f64(ax) + Dd::from_f64(diff.hi * (x * 0.5))
    }
}

impl std::ops::Neg for Dd {
    type Output = Dd;
    fn neg(self) -> Dd {
        Dd { hi: -self.hi, lo: -self.lo }
    }
}

impl std::ops::Add for Dd {
    type Output = Dd;
    fn add(self, rhs: Dd) -> Dd {
        let (s1, s2) = two_sum(self.hi, rhs.hi);
        let (t1, t2) = two_sum(self.lo, rhs.lo);
        let lo = s2 + t1;
        let (s3, lo2) = fast_two_sum(s1, lo);
        let lo3 = lo2 + t2;
        let (hi, lo) = fast_two_sum(s3, lo3);
        Dd { hi, lo }
    }
}

impl std::ops::Sub for Dd {
    type Output = Dd;
    fn sub(self, rhs: Dd) -> Dd {
        let (s1, s2) = two_diff(self.hi, rhs.hi);
        let (t1, t2) = two_diff(self.lo, rhs.lo);
        let lo = s2 + t1;
        let (s3, lo2) = fast_two_sum(s1, lo);
        let lo3 = lo2 + t2;
        let (hi, lo) = fast_two_sum(s3, lo3);
        Dd { hi, lo }
    }
}

impl std::ops::Mul for Dd {
    type Output = Dd;
    fn mul(self, rhs: Dd) -> Dd {
        let (p1, p2) = two_prod(self.hi, rhs.hi);
        let p2 = p2 + (self.hi * rhs.lo + self.lo * rhs.hi);
        let (hi, lo) = fast_two_sum(p1, p2);
        Dd { hi, lo }
    }
}

impl std::ops::Div for Dd {
    type Output = Dd;
    fn div(self, rhs: Dd) -> Dd {
        let q1 = self.hi / rhs.hi;
        let r = self - rhs * Dd::from_f64(q1);
        let q2 = r.hi / rhs.hi;
        let r2 = r - rhs * Dd::from_f64(q2);
        let q3 = r2.hi / rhs.hi;
        let (hi, lo) = fast_two_sum(q1, q2);
        Dd { hi, lo } + Dd::from_f64(q3)
    }
}

impl std::ops::Add<f64> for Dd {
    type Output = Dd;
    fn add(self, rhs: f64) -> Dd {
        self + Dd::from_f64(rhs)
    }
}

impl std::ops::Sub<f64> for Dd {
    type Output = Dd;
    fn sub(self, rhs: f64) -> Dd {
        self - Dd::from_f64(rhs)
    }
}

impl std::ops::Mul<f64> for Dd {
    type Output = Dd;
    fn mul(self, rhs: f64) -> Dd {
        self * Dd::from_f64(rhs)
    }
}

impl std::ops::Div<f64> for Dd {
    type Output = Dd;
    fn div(self, rhs: f64) -> Dd {
        self / Dd::from_f64(rhs)
    }
}

impl std::ops::Add<Dd> for f64 {
    type Output = Dd;
    fn add(self, rhs: Dd) -> Dd {
        Dd::from_f64(self) + rhs
    }
}

impl std::ops::Sub<Dd> for f64 {
    type Output = Dd;
    fn sub(self, rhs: Dd) -> Dd {
        Dd::from_f64(self) - rhs
    }
}

impl std::ops::Mul<Dd> for f64 {
    type Output = Dd;
    fn mul(self, rhs: Dd) -> Dd {
        Dd::from_f64(self) * rhs
    }
}

impl std::ops::Div<Dd> for f64 {
    type Output = Dd;
    fn div(self, rhs: Dd) -> Dd {
        Dd::from_f64(self) / rhs
    }
}

impl PartialEq for Dd {
    fn eq(&self, other: &Self) -> bool {
        self.hi == other.hi && self.lo == other.lo
    }
}

impl PartialOrd for Dd {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        match self.hi.partial_cmp(&other.hi) {
            Some(core::cmp::Ordering::Equal) => self.lo.partial_cmp(&other.lo),
            ord => ord,
        }
    }
}

impl std::fmt::Display for Dd {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:e} + {:e}", self.hi, self.lo)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_f64() {
        for &x in &[0.0, 1.0, -3.5, 1e16, 1e-16, std::f64::consts::PI] {
            assert_eq!(Dd::from_f64(x).to_f64(), x);
        }
    }

    #[test]
    fn two_sum_recovers_lost_bits() {
        let (h, l) = two_sum(1e16, 1.0);
        assert_eq!(h, 1e16);
        assert_eq!(l, 1.0);
    }

    #[test]
    fn two_diff_is_exact() {
        // 1e16 - 1.0 rounds back to 1e16, so the error term must be -1.0.
        let (h, l) = two_diff(1e16, 1.0);
        assert_eq!(h, 1e16);
        assert_eq!(l, -1.0);
    }

    #[test]
    fn two_prod_non_exact() {
        let a = 1.0 + 1e-10;
        let b = 1.0 + 1e-10;
        let (_, l) = two_prod(a, b);
        assert!(l != 0.0);
        assert_eq!((Dd::from_f64(a) * Dd::from_f64(b)).to_f64(), a * b);
    }

    #[test]
    fn algebraic_identity_hard_case() {
        // Naive f64 ((a + b) - b) would lose the 1.0 because 1.0 is below the
        // ulp of 1e20; the double-double carry recovers it exactly.
        let a = Dd::from_f64(1.0);
        let b = Dd::from_f64(1e20);
        assert!(((a + b) - b - a).to_f64().abs() < 1e-15);
    }

    #[test]
    fn sqrt_of_two() {
        let two = Dd::from_f64(2.0);
        let s = two.sqrt();
        assert!(((s * s) - two).to_f64().abs() < 1e-30);
        assert!((s.to_f64() - 2f64.sqrt()).abs() < 1e-15);
    }

    #[test]
    fn sqrt_edge_cases() {
        assert_eq!(Dd::ZERO.sqrt().to_f64(), 0.0);
        assert!(Dd::from_f64(-1.0).sqrt().to_f64().is_nan());
    }

    #[test]
    fn reciprocal_product_sum() {
        // Hilbert-style framing: each (1/(k+1)) * (k+1) is exactly 1 in
        // double-double, so the accumulated sum lands on n with ~31 digits.
        let n = 12;
        let mut acc = Dd::ZERO;
        for k in 0..n {
            acc = acc + (Dd::ONE / Dd::from_f64((k + 1) as f64)) * Dd::from_f64((k + 1) as f64);
        }
        assert!((acc - Dd::from_f64(n as f64)).to_f64().abs() < 1e-28);
    }

    #[test]
    fn ill_conditioned_cancellation() {
        // Hilbert-matrix-style ill-conditioning: Dd recovers what naive f64
        // loses. Summing [1e30, 1.0, -1e30] must yield exactly 1.0.
        let mut acc = Dd::ZERO;
        for &t in &[1e30, 1.0, -1e30] {
            acc = acc + Dd::from_f64(t);
        }
        assert_eq!(acc.to_f64(), 1.0);

        // The same sum in plain f64 collapses to 0.0 (the 1.0 is annihilated).
        let mut f = 0.0f64;
        for &t in &[1e30, 1.0, -1e30] {
            f += t;
        }
        assert_eq!(f, 0.0);
    }

    #[test]
    fn mixed_f64_operators() {
        assert_eq!((Dd::ONE + 1.0).to_f64(), 2.0);
        assert_eq!((1.0 + Dd::ONE).to_f64(), 2.0);
        assert_eq!((Dd::from_f64(3.0) - 1.0).to_f64(), 2.0);
        assert_eq!((3.0 - Dd::ONE).to_f64(), 2.0);
        assert_eq!((Dd::from_f64(3.0) * 2.0).to_f64(), 6.0);
        assert_eq!((2.0 * Dd::from_f64(3.0)).to_f64(), 6.0);
        assert_eq!((Dd::from_f64(6.0) / 2.0).to_f64(), 3.0);
        assert_eq!((6.0 / Dd::from_f64(2.0)).to_f64(), 3.0);
    }

    #[test]
    fn abs_finite_constants_and_display() {
        assert!(Dd::from_f64(1.0).is_finite());
        assert_eq!(Dd::from_f64(-2.0).abs().to_f64(), 2.0);
        assert_eq!(Dd::ZERO.to_f64(), 0.0);
        assert_eq!(Dd::ONE.to_f64(), 1.0);
        let rendered = format!("{}", Dd::from_f64(1.5));
        assert!(rendered.contains('e'));
    }

    #[test]
    fn ordering_and_equality() {
        assert!(Dd::from_f64(1.0) < Dd::from_f64(2.0));
        assert!(Dd::from_f64(1.0) == Dd::from_f64(1.0));
        let a = Dd::new(1.0, 1e-30);
        let b = Dd::from_f64(1.0);
        assert!(a > b);
    }

    #[test]
    fn new_canonicalizes() {
        let x = Dd::new(1.0, 1e-20);
        assert_eq!(x.hi, 1.0);
        assert_eq!(x.hi + x.lo, 1.0 + 1e-20);
    }

    #[test]
    fn division_recovers_precision() {
        let a = Dd::from_f64(1.0);
        let b = Dd::from_f64(3.0);
        let q = a / b;
        // q * 3 should reconstruct 1.0 to within double-double precision.
        assert!(((q * b) - a).to_f64().abs() < 1e-30);
    }
}
