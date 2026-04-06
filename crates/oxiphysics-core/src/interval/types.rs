//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Linear constraint: `a[0]*x[0] + a[1]*x[1] + ... = rhs` (within tolerance).
pub struct LinearConstraint {
    /// Coefficients (index, coefficient).
    pub coeffs: Vec<(usize, f64)>,
    /// Right-hand side interval.
    pub rhs: Interval,
}
/// A Taylor model represents a polynomial approximation plus an interval
/// remainder: `p(x) + r` where `p` is a polynomial with f64 coefficients
/// and `r` is an interval remainder enclosing the truncation error.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TaylorModel {
    /// Polynomial coefficients in ascending degree order: poly\[i\] = coeff of x^i
    pub poly: Vec<f64>,
    /// Interval remainder bound
    pub remainder: Interval,
}
#[allow(dead_code)]
impl TaylorModel {
    /// Create a new Taylor model with given polynomial coefficients and remainder.
    pub fn new(poly: Vec<f64>, remainder: Interval) -> Self {
        TaylorModel { poly, remainder }
    }
    /// Evaluate the Taylor model at a point interval `x` returning an enclosure.
    /// Uses Horner's method for the polynomial part.
    pub fn evaluate(&self, x: Interval) -> Interval {
        let n = self.poly.len();
        if n == 0 {
            return self.remainder;
        }
        let mut result = Interval::point(self.poly[n - 1]);
        for i in (0..n - 1).rev() {
            result = result * x + Interval::point(self.poly[i]);
        }
        result + self.remainder
    }
    /// Add two Taylor models (polynomial addition + remainder addition).
    pub fn add(&self, other: &TaylorModel) -> TaylorModel {
        let n = self.poly.len().max(other.poly.len());
        let mut poly = vec![0.0; n];
        for (i, &c) in self.poly.iter().enumerate() {
            poly[i] += c;
        }
        for (i, &c) in other.poly.iter().enumerate() {
            poly[i] += c;
        }
        TaylorModel {
            poly,
            remainder: self.remainder + other.remainder,
        }
    }
    /// Scale the Taylor model by a scalar factor.
    pub fn scale(&self, s: f64) -> TaylorModel {
        TaylorModel {
            poly: self.poly.iter().map(|&c| c * s).collect(),
            remainder: self.remainder * Interval::point(s),
        }
    }
}
/// Result of the interval Newton method.
#[derive(Debug, Clone)]
pub enum IntervalNewtonResult {
    /// A root is guaranteed to exist in this interval.
    Root(Interval),
    /// No root exists in the searched interval.
    NoRoot,
    /// Multiple candidate intervals (when bisection was needed).
    Multiple(Vec<Interval>),
}
/// A 3D interval (axis-aligned bounding box) represented as three independent
/// [`Interval`]s along the x, y, and z axes.
#[derive(Debug, Clone, Copy)]
pub struct Interval3 {
    /// Interval along the x axis.
    pub x: Interval,
    /// Interval along the y axis.
    pub y: Interval,
    /// Interval along the z axis.
    pub z: Interval,
}
impl Interval3 {
    /// Create a new `Interval3` from three intervals.
    pub fn new(x: Interval, y: Interval, z: Interval) -> Self {
        Self { x, y, z }
    }
    /// Create from AABB min/max corner arrays.
    pub fn from_aabb(min: [f64; 3], max: [f64; 3]) -> Self {
        Self {
            x: Interval::new(min[0], max[0]),
            y: Interval::new(min[1], max[1]),
            z: Interval::new(min[2], max[2]),
        }
    }
    /// Convert to AABB min/max corner arrays.
    pub fn to_aabb(&self) -> ([f64; 3], [f64; 3]) {
        (
            [self.x.lo, self.y.lo, self.z.lo],
            [self.x.hi, self.y.hi, self.z.hi],
        )
    }
    /// Returns `true` if this `Interval3` overlaps with `other` in all three axes.
    pub fn overlaps_3d(&self, other: &Self) -> bool {
        self.x.overlaps(&other.x) && self.y.overlaps(&other.y) && self.z.overlaps(&other.z)
    }
    /// Hull of two 3D intervals.
    pub fn hull_3d(a: Self, b: Self) -> Self {
        Self {
            x: Interval::hull(a.x, b.x),
            y: Interval::hull(a.y, b.y),
            z: Interval::hull(a.z, b.z),
        }
    }
    /// Intersection of two 3D intervals. Returns `None` if disjoint in any axis.
    pub fn intersection_3d(a: Self, b: Self) -> Option<Self> {
        let x = Interval::intersection(a.x, b.x)?;
        let y = Interval::intersection(a.y, b.y)?;
        let z = Interval::intersection(a.z, b.z)?;
        Some(Self { x, y, z })
    }
    /// Volume of the box.
    pub fn volume(&self) -> f64 {
        self.x.width() * self.y.width() * self.z.width()
    }
    /// Surface area of the box.
    pub fn surface_area(&self) -> f64 {
        let wx = self.x.width();
        let wy = self.y.width();
        let wz = self.z.width();
        2.0 * (wx * wy + wy * wz + wz * wx)
    }
    /// Center point of the box.
    pub fn center(&self) -> [f64; 3] {
        [self.x.midpoint(), self.y.midpoint(), self.z.midpoint()]
    }
    /// Inflate the box by `delta` on each side in each dimension.
    pub fn inflate(&self, delta: f64) -> Self {
        Self {
            x: self.x.inflate(delta),
            y: self.y.inflate(delta),
            z: self.z.inflate(delta),
        }
    }
    /// Check if a point is inside the box.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        self.x.contains(p[0]) && self.y.contains(p[1]) && self.z.contains(p[2])
    }
    /// Compute conservative bounds on `|v|²` for all `v` in this interval box.
    ///
    /// Returns an interval `[lo, hi]` such that `lo ≤ |v|² ≤ hi` for any
    /// point `v` inside `self`.
    pub fn length_sq_bounds(&self) -> Interval {
        let x2 = square(self.x);
        let y2 = square(self.y);
        let z2 = square(self.z);
        Interval {
            lo: x2.lo + y2.lo + z2.lo,
            hi: x2.hi + y2.hi + z2.hi,
        }
    }
}
/// A 2×2 interval matrix `[[a,b],[c,d]]`.
#[derive(Debug, Clone, Copy)]
pub struct IntervalMatrix2 {
    /// Row 0, columns 0 and 1.
    pub row0: [Interval; 2],
    /// Row 1, columns 0 and 1.
    pub row1: [Interval; 2],
}
impl IntervalMatrix2 {
    /// Creates a 2×2 interval matrix from rows.
    pub fn new(row0: [Interval; 2], row1: [Interval; 2]) -> Self {
        Self { row0, row1 }
    }
    /// Point (degenerate) matrix from `f64` values.
    pub fn from_f64(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self {
            row0: [Interval::point(a), Interval::point(b)],
            row1: [Interval::point(c), Interval::point(d)],
        }
    }
    /// Identity matrix `[[1,0],[0,1]]`.
    pub fn identity() -> Self {
        Self::from_f64(1.0, 0.0, 0.0, 1.0)
    }
    /// Matrix–vector product: returns `[row0·v, row1·v]`.
    pub fn mul_vec2(&self, v: [Interval; 2]) -> [Interval; 2] {
        [
            self.row0[0] * v[0] + self.row0[1] * v[1],
            self.row1[0] * v[0] + self.row1[1] * v[1],
        ]
    }
    /// Matrix–matrix product.
    pub fn mul_mat2(&self, other: &Self) -> Self {
        let c00 = self.row0[0] * other.row0[0] + self.row0[1] * other.row1[0];
        let c01 = self.row0[0] * other.row0[1] + self.row0[1] * other.row1[1];
        let c10 = self.row1[0] * other.row0[0] + self.row1[1] * other.row1[0];
        let c11 = self.row1[0] * other.row0[1] + self.row1[1] * other.row1[1];
        Self {
            row0: [c00, c01],
            row1: [c10, c11],
        }
    }
    /// Transpose.
    pub fn transpose2(&self) -> Self {
        Self {
            row0: [self.row0[0], self.row1[0]],
            row1: [self.row0[1], self.row1[1]],
        }
    }
}
/// A closed interval `[lo, hi]` on the real line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}
impl Interval {
    /// Create a new interval `[lo, hi]`. Panics if `lo > hi`.
    pub fn new(lo: f64, hi: f64) -> Self {
        assert!(lo <= hi, "Interval::new: lo ({lo}) must be ≤ hi ({hi})");
        Self { lo, hi }
    }
    /// Create a degenerate (point) interval `[v, v]`.
    pub fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }
    /// The empty interval (used as identity for hull operations).
    pub fn empty() -> Self {
        Self {
            lo: f64::INFINITY,
            hi: f64::NEG_INFINITY,
        }
    }
    /// Returns `true` if this is an empty interval (lo > hi).
    pub fn is_empty(&self) -> bool {
        self.lo > self.hi
    }
    /// Width of the interval: `hi - lo`.
    pub fn width(&self) -> f64 {
        self.hi - self.lo
    }
    /// Midpoint of the interval: `(lo + hi) / 2`.
    pub fn midpoint(&self) -> f64 {
        (self.lo + self.hi) * 0.5
    }
    /// Radius of the interval: `(hi - lo) / 2`.
    pub fn radius(&self) -> f64 {
        (self.hi - self.lo) * 0.5
    }
    /// Returns `true` if `v` is within `[lo, hi]`.
    pub fn contains(&self, v: f64) -> bool {
        v >= self.lo && v <= self.hi
    }
    /// Returns `true` if `other` is entirely contained within `self`.
    pub fn contains_interval(&self, other: &Self) -> bool {
        self.lo <= other.lo && other.hi <= self.hi
    }
    /// Returns `true` if this interval overlaps with `other`.
    ///
    /// Two intervals overlap if they share at least one point.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.lo <= other.hi && other.lo <= self.hi
    }
    /// Smallest interval containing both `a` and `b` (convex hull).
    pub fn hull(a: Self, b: Self) -> Self {
        Self {
            lo: a.lo.min(b.lo),
            hi: a.hi.max(b.hi),
        }
    }
    /// Hull of multiple intervals.
    pub fn hull_many(intervals: &[Self]) -> Self {
        let mut result = Self::empty();
        for &iv in intervals {
            result = Self {
                lo: result.lo.min(iv.lo),
                hi: result.hi.max(iv.hi),
            };
        }
        result
    }
    /// Intersection of two intervals. Returns `None` if disjoint.
    pub fn intersection(a: Self, b: Self) -> Option<Self> {
        let lo = a.lo.max(b.lo);
        let hi = a.hi.min(b.hi);
        if lo <= hi {
            Some(Self { lo, hi })
        } else {
            None
        }
    }
    /// Split this interval at the midpoint, returning (left, right).
    pub fn bisect(&self) -> (Self, Self) {
        let m = self.midpoint();
        (Self::new(self.lo, m), Self::new(m, self.hi))
    }
    /// Split this interval at a given point. Panics if `v` is outside.
    pub fn split_at(&self, v: f64) -> (Self, Self) {
        assert!(
            self.contains(v),
            "split_at: {} not in [{}, {}]",
            v,
            self.lo,
            self.hi
        );
        (Self::new(self.lo, v), Self::new(v, self.hi))
    }
    /// Inflate the interval by `delta` on each side.
    pub fn inflate(&self, delta: f64) -> Self {
        Self {
            lo: self.lo - delta,
            hi: self.hi + delta,
        }
    }
    /// Reciprocal of an interval. Panics if zero is contained.
    pub fn reciprocal(&self) -> Self {
        assert!(
            !self.contains(0.0),
            "reciprocal: interval [{}, {}] contains zero",
            self.lo,
            self.hi
        );
        Self {
            lo: (1.0 / self.hi),
            hi: (1.0 / self.lo),
        }
    }
    /// Integer power of an interval.
    ///
    /// For even n: result lo = 0 if 0 ∈ self, else min(|lo|^n, |hi|^n); hi = max(|lo|^n, |hi|^n).
    /// For odd n: monotone, lo = lo^n, hi = hi^n.
    pub fn powi(&self, n: i32) -> Self {
        if n == 0 {
            return Self::point(1.0);
        }
        if n < 0 {
            return self.reciprocal().powi(-n);
        }
        if n % 2 == 0 {
            let a = self.lo.powi(n).abs();
            let b = self.hi.powi(n).abs();
            let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
                0.0
            } else {
                a.min(b)
            };
            Self { lo, hi: a.max(b) }
        } else {
            Self {
                lo: self.lo.powi(n),
                hi: self.hi.powi(n),
            }
        }
    }
    /// Square root of an interval. Clamps negative lower bound to 0.
    pub fn sqrt(&self) -> Self {
        let lo = self.lo.max(0.0).sqrt();
        let hi = self.hi.max(0.0).sqrt();
        Self { lo, hi }
    }
    /// Absolute value of an interval.
    pub fn abs(&self) -> Self {
        let a = self.lo.abs();
        let b = self.hi.abs();
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            a.min(b)
        };
        Self { lo, hi: a.max(b) }
    }
}
/// An affine form: `x0 + x1*e1 + x2*e2 + ... + xn*en`
///
/// where `ei ∈ [-1, 1]` are noise symbols.
///
/// Affine arithmetic provides tighter bounds than plain interval arithmetic
/// because it tracks correlations between variables through shared noise symbols.
#[derive(Debug, Clone)]
pub struct AffineForm {
    /// Central value.
    pub center: f64,
    /// Coefficients for noise symbols (index = noise symbol id).
    pub terms: Vec<(usize, f64)>,
}
impl AffineForm {
    /// Create from a constant value (no noise).
    pub fn constant(v: f64) -> Self {
        Self {
            center: v,
            terms: Vec::new(),
        }
    }
    /// Create from an interval, introducing a fresh noise symbol.
    pub fn from_interval(iv: Interval) -> Self {
        let center = iv.midpoint();
        let radius = iv.radius();
        if radius < 1e-30 {
            Self::constant(center)
        } else {
            let id = new_noise_id();
            Self {
                center,
                terms: vec![(id, radius)],
            }
        }
    }
    /// Convert back to an interval by summing absolute values of all noise terms.
    pub fn to_interval(&self) -> Interval {
        let radius: f64 = self.terms.iter().map(|(_, c)| c.abs()).sum();
        Interval::new(self.center - radius, self.center + radius)
    }
    /// Affine addition: (a + b).
    pub fn add(&self, other: &Self) -> Self {
        let center = self.center + other.center;
        let mut terms = self.terms.clone();
        for &(id, coeff) in &other.terms {
            if let Some(entry) = terms.iter_mut().find(|(eid, _)| *eid == id) {
                entry.1 += coeff;
            } else {
                terms.push((id, coeff));
            }
        }
        terms.retain(|(_, c)| c.abs() > 1e-30);
        Self { center, terms }
    }
    /// Affine subtraction: (a - b).
    pub fn sub(&self, other: &Self) -> Self {
        let center = self.center - other.center;
        let mut terms = self.terms.clone();
        for &(id, coeff) in &other.terms {
            if let Some(entry) = terms.iter_mut().find(|(eid, _)| *eid == id) {
                entry.1 -= coeff;
            } else {
                terms.push((id, -coeff));
            }
        }
        terms.retain(|(_, c)| c.abs() > 1e-30);
        Self { center, terms }
    }
    /// Affine scalar multiplication.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            center: self.center * s,
            terms: self.terms.iter().map(|&(id, c)| (id, c * s)).collect(),
        }
    }
    /// Affine multiplication (introduces a new noise symbol for the nonlinear part).
    pub fn mul(&self, other: &Self) -> Self {
        let center = self.center * other.center;
        let mut terms = Vec::new();
        for &(id, coeff) in &other.terms {
            terms.push((id, self.center * coeff));
        }
        for &(id, coeff) in &self.terms {
            if let Some(entry) = terms.iter_mut().find(|(eid, _)| *eid == id) {
                entry.1 += other.center * coeff;
            } else {
                terms.push((id, other.center * coeff));
            }
        }
        let ra: f64 = self.terms.iter().map(|(_, c)| c.abs()).sum();
        let rb: f64 = other.terms.iter().map(|(_, c)| c.abs()).sum();
        let nonlinear_error = ra * rb;
        if nonlinear_error > 1e-30 {
            let id = new_noise_id();
            terms.push((id, nonlinear_error));
        }
        terms.retain(|(_, c)| c.abs() > 1e-30);
        Self { center, terms }
    }
}
/// Bound constraint: `x[idx] ∈ [lo, hi]`.
pub struct BoundConstraint {
    /// Variable index.
    pub idx: usize,
    /// Allowed bounds.
    pub bounds: Interval,
}
/// A 3×3 interval matrix stored as three row arrays.
#[derive(Debug, Clone, Copy)]
pub struct IntervalMatrix3 {
    /// Rows 0, 1, 2 (each has 3 columns).
    pub rows: [[Interval; 3]; 3],
}
impl IntervalMatrix3 {
    /// Creates a 3×3 identity interval matrix.
    pub fn identity() -> Self {
        let zero = Interval::point(0.0);
        let one = Interval::point(1.0);
        Self {
            rows: [[one, zero, zero], [zero, one, zero], [zero, zero, one]],
        }
    }
    /// Matrix–vector multiply: `M * v` where `v` is a 3-element interval array.
    pub fn mul_vec3(&self, v: [Interval; 3]) -> [Interval; 3] {
        let dot = |row: &[Interval; 3]| row[0] * v[0] + row[1] * v[1] + row[2] * v[2];
        [dot(&self.rows[0]), dot(&self.rows[1]), dot(&self.rows[2])]
    }
    /// Matrix–matrix multiply.
    pub fn mul_mat3(&self, other: &Self) -> Self {
        let mut res = Self::identity();
        for i in 0..3 {
            for j in 0..3 {
                let mut s = Interval::point(0.0);
                for k in 0..3 {
                    s = s + self.rows[i][k] * other.rows[k][j];
                }
                res.rows[i][j] = s;
            }
        }
        res
    }
    /// Transpose.
    pub fn transpose3(&self) -> Self {
        let mut res = Self::identity();
        for i in 0..3 {
            for j in 0..3 {
                res.rows[i][j] = self.rows[j][i];
            }
        }
        res
    }
}
/// An interval matrix with `rows` rows and `cols` columns.
#[derive(Debug, Clone)]
pub struct IntervalMatrix {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Row-major interval entries.
    pub data: Vec<Interval>,
}
impl IntervalMatrix {
    /// Create a zero interval matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![Interval::point(0.0); rows * cols],
        }
    }
    /// Create an identity interval matrix (must be square).
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.data[i * n + i] = Interval::point(1.0);
        }
        m
    }
    /// Get element at (row, col).
    pub fn get(&self, row: usize, col: usize) -> Interval {
        self.data[row * self.cols + col]
    }
    /// Set element at (row, col).
    pub fn set(&mut self, row: usize, col: usize, val: Interval) {
        self.data[row * self.cols + col] = val;
    }
    /// Matrix-vector multiply: returns an IntervalVec.
    pub fn mul_vec(&self, v: &IntervalVec) -> IntervalVec {
        assert_eq!(
            self.cols,
            v.dim(),
            "IntervalMatrix::mul_vec: dimension mismatch"
        );
        let mut result = Vec::with_capacity(self.rows);
        for i in 0..self.rows {
            let mut sum = Interval::point(0.0);
            for j in 0..self.cols {
                sum = sum + self.get(i, j) * v.data[j];
            }
            result.push(sum);
        }
        IntervalVec { data: result }
    }
    /// Matrix-matrix multiply.
    pub fn mul_mat(&self, other: &Self) -> Self {
        assert_eq!(
            self.cols, other.rows,
            "IntervalMatrix::mul_mat: dimension mismatch"
        );
        let mut result = Self::zeros(self.rows, other.cols);
        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut sum = Interval::point(0.0);
                for k in 0..self.cols {
                    sum = sum + self.get(i, k) * other.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        result
    }
    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut result = Self::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(j, i, self.get(i, j));
            }
        }
        result
    }
}
/// An interval vector of fixed dimension N (stored as a Vec for flexibility).
#[derive(Debug, Clone)]
pub struct IntervalVec {
    /// The components.
    pub data: Vec<Interval>,
}
impl IntervalVec {
    /// Create from a slice of intervals.
    pub fn from_slice(intervals: &[Interval]) -> Self {
        Self {
            data: intervals.to_vec(),
        }
    }
    /// Dimension of the vector.
    pub fn dim(&self) -> usize {
        self.data.len()
    }
    /// Add two interval vectors component-wise.
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(
            self.dim(),
            other.dim(),
            "IntervalVec::add: dimension mismatch"
        );
        Self {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        }
    }
    /// Subtract two interval vectors component-wise.
    pub fn sub(&self, other: &Self) -> Self {
        assert_eq!(
            self.dim(),
            other.dim(),
            "IntervalVec::sub: dimension mismatch"
        );
        Self {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| a - b)
                .collect(),
        }
    }
    /// Interval dot product: sum of component-wise products.
    pub fn dot(&self, other: &Self) -> Interval {
        assert_eq!(
            self.dim(),
            other.dim(),
            "IntervalVec::dot: dimension mismatch"
        );
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a * b)
            .fold(Interval::point(0.0), |acc, x| acc + x)
    }
    /// Scale all components by a scalar interval.
    pub fn scale(&self, s: Interval) -> Self {
        Self {
            data: self.data.iter().map(|&a| a * s).collect(),
        }
    }
    /// Hull of two interval vectors (component-wise hull).
    pub fn hull(&self, other: &Self) -> Self {
        assert_eq!(
            self.dim(),
            other.dim(),
            "IntervalVec::hull: dimension mismatch"
        );
        Self {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| Interval::hull(a, b))
                .collect(),
        }
    }
    /// Maximum width among all components.
    pub fn max_width(&self) -> f64 {
        self.data
            .iter()
            .map(|iv| iv.width())
            .fold(0.0_f64, f64::max)
    }
    /// Midpoint vector.
    pub fn midpoint(&self) -> Vec<f64> {
        self.data.iter().map(|iv| iv.midpoint()).collect()
    }
}
