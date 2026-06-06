//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// A dense N×N complex matrix stored in row-major order.
pub struct ComplexMatN {
    /// Number of rows and columns.
    pub n: usize,
    /// Row-major storage: data\[i*n + j\] is entry (i,j).
    pub data: Vec<Complex>,
}
impl ComplexMatN {
    /// Create a new N×N zero matrix.
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![Complex::zero(); n * n],
        }
    }
    /// Create the N×N identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n);
        for i in 0..n {
            m.data[i * n + i] = Complex::one();
        }
        m
    }
    /// Get entry (i, j).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Complex {
        self.data[i * self.n + j]
    }
    /// Set entry (i, j).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: Complex) {
        self.data[i * self.n + j] = v;
    }
    /// Matrix-matrix multiplication `self * other`.
    pub fn mul(&self, other: &ComplexMatN) -> ComplexMatN {
        assert_eq!(self.n, other.n, "ComplexMatN::mul: dimension mismatch");
        let n = self.n;
        let mut result = ComplexMatN::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let a_ik = self.get(i, k);
                if a_ik.norm_sq() < f64::EPSILON * f64::EPSILON {
                    continue;
                }
                for j in 0..n {
                    let old = result.get(i, j);
                    result.set(i, j, old + a_ik * other.get(k, j));
                }
            }
        }
        result
    }
    /// Conjugate transpose (Hermitian adjoint).
    pub fn conj_transpose(&self) -> ComplexMatN {
        let n = self.n;
        let mut result = ComplexMatN::zeros(n);
        for i in 0..n {
            for j in 0..n {
                result.set(i, j, self.get(j, i).conj());
            }
        }
        result
    }
    /// Frobenius norm: sqrt(Σ |a_ij|²).
    pub fn frobenius_norm(&self) -> f64 {
        self.data.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt()
    }
    /// Trace: Σ a_ii.
    pub fn trace(&self) -> Complex {
        (0..self.n)
            .map(|i| self.get(i, i))
            .fold(Complex::zero(), |acc, z| acc + z)
    }
    /// Matrix-vector product `self * v`.
    pub fn apply(&self, v: &[Complex]) -> Vec<Complex> {
        assert_eq!(v.len(), self.n, "ComplexMatN::apply: dimension mismatch");
        (0..self.n)
            .map(|i| {
                (0..self.n)
                    .map(|j| self.get(i, j) * v[j])
                    .fold(Complex::zero(), |acc, z| acc + z)
            })
            .collect()
    }
    /// Matrix multiplication `self * other` (alias for `mul`).
    pub fn matmul(&self, other: &ComplexMatN) -> ComplexMatN {
        self.mul(other)
    }
    /// Matrix-vector multiplication `self * v` (alias for `apply`).
    pub fn matvec(&self, v: &[Complex]) -> Vec<Complex> {
        self.apply(v)
    }
    /// Conjugate transpose (alias for `conj_transpose`).
    pub fn conjugate_transpose(&self) -> ComplexMatN {
        self.conj_transpose()
    }
    /// Element-wise add two matrices.
    pub fn add(&self, other: &ComplexMatN) -> ComplexMatN {
        assert_eq!(self.n, other.n, "ComplexMatN::add: dimension mismatch");
        ComplexMatN {
            n: self.n,
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        }
    }
    /// Scale all entries by a complex scalar.
    pub fn scale(&self, s: Complex) -> ComplexMatN {
        ComplexMatN {
            n: self.n,
            data: self.data.iter().map(|&z| z * s).collect(),
        }
    }
}
/// A quaternion `w + x·i + y·j + z·k` supporting the full Hamilton algebra.
///
/// This type is intentionally separate from `math::Quat` (which is a nalgebra
/// `UnitQuaternion`) and provides raw (possibly non-unit) quaternion arithmetic
/// together with transcendental operations and SLERP.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuatAlgebra {
    /// Scalar part.
    pub w: f64,
    /// `i` component.
    pub x: f64,
    /// `j` component.
    pub y: f64,
    /// `k` component.
    pub z: f64,
}
impl QuatAlgebra {
    /// Create a quaternion from components.
    #[inline]
    pub fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }
    /// Multiplicative identity: `1 + 0·i + 0·j + 0·k`.
    #[inline]
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0)
    }
    /// Additive identity: `0 + 0·i + 0·j + 0·k`.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }
    /// Construct a unit quaternion from an axis and angle (radians).
    ///
    /// The axis need not be pre-normalised; it is normalised internally.
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if len < f64::EPSILON {
            return Self::identity();
        }
        let half = angle * 0.5;
        let s = half.sin() / len;
        Self::new(half.cos(), axis[0] * s, axis[1] * s, axis[2] * s)
    }
    /// Convert a unit quaternion to a 3×3 rotation matrix.
    pub fn to_rotation_matrix(&self) -> [[f64; 3]; 3] {
        let (w, x, y, z) = (self.w, self.x, self.y, self.z);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }
    /// Rotate a 3-vector by this (assumed unit) quaternion: `q v q*`.
    pub fn rotate_vector(&self, v: [f64; 3]) -> [f64; 3] {
        let r = self.to_rotation_matrix();
        [
            r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
            r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
            r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
        ]
    }
    /// Quaternion conjugate: `w - x·i - y·j - z·k`.
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }
    /// Squared norm: `w² + x² + y² + z²`.
    #[inline]
    pub fn norm_sq(&self) -> f64 {
        self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
    }
    /// Euclidean norm.
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }
    /// Return a normalised copy (unit quaternion).
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        Self::new(self.w / n, self.x / n, self.y / n, self.z / n)
    }
    /// Multiplicative inverse: `q* / |q|²`.
    pub fn inverse(&self) -> Self {
        let n2 = self.norm_sq();
        let c = self.conj();
        Self::new(c.w / n2, c.x / n2, c.y / n2, c.z / n2)
    }
    /// Dot product of two quaternions (treating them as 4-vectors).
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z
    }
    /// Quaternion exponential: `exp(q) = e^w (cos|v| + v̂ sin|v|)`.
    pub fn exp(&self) -> Self {
        let v_norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        let ew = self.w.exp();
        if v_norm < f64::EPSILON {
            return Self::new(ew, 0.0, 0.0, 0.0);
        }
        let s = ew * v_norm.sin() / v_norm;
        Self::new(ew * v_norm.cos(), self.x * s, self.y * s, self.z * s)
    }
    /// Quaternion logarithm: `ln(q) = ln|q| + v̂ arccos(w/|q|)`.
    pub fn ln(&self) -> Self {
        let n = self.norm();
        let v_norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        let ln_n = n.ln();
        if v_norm < f64::EPSILON {
            return Self::new(ln_n, 0.0, 0.0, 0.0);
        }
        let theta = (self.w / n).clamp(-1.0, 1.0).acos();
        let s = theta / v_norm;
        Self::new(ln_n, self.x * s, self.y * s, self.z * s)
    }
    /// Raise to a real power: `q^t = exp(t · ln q)`.
    pub fn pow(&self, t: f64) -> Self {
        let l = self.ln();
        Self::new(l.w * t, l.x * t, l.y * t, l.z * t).exp()
    }
    /// Spherical linear interpolation (SLERP) between two unit quaternions.
    pub fn slerp(&self, other: &Self, t: f64) -> Self {
        let mut dot = self.dot(other).clamp(-1.0, 1.0);
        let other_adj = if dot < 0.0 {
            dot = -dot;
            Self::new(-other.w, -other.x, -other.y, -other.z)
        } else {
            *other
        };
        if dot > 1.0 - f64::EPSILON {
            let q = Self::new(
                self.w + t * (other_adj.w - self.w),
                self.x + t * (other_adj.x - self.x),
                self.y + t * (other_adj.y - self.y),
                self.z + t * (other_adj.z - self.z),
            );
            return q.normalize();
        }
        let theta = dot.acos();
        let sin_theta = theta.sin();
        let scale0 = ((1.0 - t) * theta).sin() / sin_theta;
        let scale1 = (t * theta).sin() / sin_theta;
        Self::new(
            scale0 * self.w + scale1 * other_adj.w,
            scale0 * self.x + scale1 * other_adj.x,
            scale0 * self.y + scale1 * other_adj.y,
            scale0 * self.z + scale1 * other_adj.z,
        )
    }
}
/// A complex number `re + im·i`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}
impl Complex {
    /// Create a new complex number.
    #[inline]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    /// Additive identity: `0 + 0·i`.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0)
    }
    /// Multiplicative identity: `1 + 0·i`.
    #[inline]
    pub fn one() -> Self {
        Self::new(1.0, 0.0)
    }
    /// The imaginary unit: `0 + 1·i`.
    #[inline]
    pub fn i() -> Self {
        Self::new(0.0, 1.0)
    }
    /// Complex conjugate: `re - im·i`.
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.re, -self.im)
    }
    /// Squared modulus: `re² + im²`.
    #[inline]
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    /// Modulus (absolute value): `√(re² + im²)`.
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }
    /// Argument (angle): `atan2(im, re)` in radians.
    #[inline]
    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }
    /// Construct from polar form: `r·(cos θ + i·sin θ)`.
    #[inline]
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }
    /// Complex exponential: `e^(re+im·i) = e^re · (cos(im) + i·sin(im))`.
    pub fn exp(&self) -> Self {
        let factor = self.re.exp();
        Self::new(factor * self.im.cos(), factor * self.im.sin())
    }
    /// Principal complex logarithm: `ln|z| + i·arg(z)`.
    pub fn ln(&self) -> Self {
        Self::new(self.norm().ln(), self.arg())
    }
    /// Principal square root.
    pub fn sqrt(&self) -> Self {
        let r = self.norm();
        let theta = self.arg();
        Self::from_polar(r.sqrt(), theta / 2.0)
    }
    /// Raise to a real power: `z^n = exp(n · ln z)`.
    pub fn pow_f64(&self, n: f64) -> Self {
        (self.ln() * Complex::new(n, 0.0)).exp()
    }
    /// Complex sine: `sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b)`.
    pub fn sin(&self) -> Self {
        Self::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh(),
        )
    }
    /// Complex cosine: `cos(a+bi) = cos(a)cosh(b) - i·sin(a)sinh(b)`.
    pub fn cos(&self) -> Self {
        Self::new(
            self.re.cos() * self.im.cosh(),
            -(self.re.sin() * self.im.sinh()),
        )
    }
}
/// A 2×2 complex matrix stored as `[[a, b\], [c, d]]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComplexMat2 {
    /// Entry (0,0).
    pub a: Complex,
    /// Entry (0,1).
    pub b: Complex,
    /// Entry (1,0).
    pub c: Complex,
    /// Entry (1,1).
    pub d: Complex,
}
impl ComplexMat2 {
    /// Create from four complex entries (row-major: a=top-left, b=top-right,
    /// c=bottom-left, d=bottom-right).
    pub fn new(a: Complex, b: Complex, c: Complex, d: Complex) -> Self {
        Self { a, b, c, d }
    }
    /// 2×2 identity matrix.
    pub fn identity() -> Self {
        Self::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::one(),
        )
    }
    /// Zero matrix.
    pub fn zero() -> Self {
        Self::new(
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
        )
    }
    /// Matrix multiplication `self * other`.
    pub fn mul(&self, other: &ComplexMat2) -> ComplexMat2 {
        ComplexMat2::new(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d,
        )
    }
    /// Matrix-vector multiplication `self * [z1, z2]`.
    pub fn apply(&self, z1: Complex, z2: Complex) -> (Complex, Complex) {
        (self.a * z1 + self.b * z2, self.c * z1 + self.d * z2)
    }
    /// Determinant: `ad - bc`.
    pub fn det(&self) -> Complex {
        self.a * self.d - self.b * self.c
    }
    /// Trace: `a + d`.
    pub fn trace(&self) -> Complex {
        self.a + self.d
    }
    /// Transpose.
    pub fn transpose(&self) -> ComplexMat2 {
        ComplexMat2::new(self.a, self.c, self.b, self.d)
    }
    /// Conjugate transpose (Hermitian adjoint).
    pub fn conjugate_transpose(&self) -> ComplexMat2 {
        ComplexMat2::new(self.a.conj(), self.c.conj(), self.b.conj(), self.d.conj())
    }
    /// Inverse of a 2×2 matrix, or `None` if the matrix is singular.
    pub fn inverse(&self) -> Option<ComplexMat2> {
        let det = self.det();
        if det.norm_sq() < 1e-24 {
            return None;
        }
        let inv_det = Complex::one() / det;
        Some(ComplexMat2::new(
            self.d * inv_det,
            -self.b * inv_det,
            -self.c * inv_det,
            self.a * inv_det,
        ))
    }
    /// Compute the two eigenvalues of the matrix.
    ///
    /// Uses the characteristic equation λ² - tr(A)λ + det(A) = 0.
    pub fn eigenvalues(&self) -> (Complex, Complex) {
        let tr = self.trace();
        let det = self.det();
        let disc = tr * tr - det * Complex::new(4.0, 0.0);
        let sqrt_disc = disc.sqrt();
        let two = Complex::new(2.0, 0.0);
        let l1 = (tr + sqrt_disc) / two;
        let l2 = (tr - sqrt_disc) / two;
        (l1, l2)
    }
    /// Add two matrices.
    pub fn add(&self, other: &ComplexMat2) -> ComplexMat2 {
        ComplexMat2::new(
            self.a + other.a,
            self.b + other.b,
            self.c + other.c,
            self.d + other.d,
        )
    }
    /// Scale by a complex scalar.
    pub fn scale(&self, s: Complex) -> ComplexMat2 {
        ComplexMat2::new(self.a * s, self.b * s, self.c * s, self.d * s)
    }
    /// Frobenius norm: sqrt(sum |a_ij|²).
    pub fn frobenius_norm(&self) -> f64 {
        (self.a.norm_sq() + self.b.norm_sq() + self.c.norm_sq() + self.d.norm_sq()).sqrt()
    }
    /// Check if the matrix is Hermitian (A = A†).
    pub fn is_hermitian(&self, tol: f64) -> bool {
        let ah = self.conjugate_transpose();
        (self.a - ah.a).norm() < tol
            && (self.b - ah.b).norm() < tol
            && (self.c - ah.c).norm() < tol
            && (self.d - ah.d).norm() < tol
    }
    /// Check if the matrix is unitary (A† A = I).
    pub fn is_unitary(&self, tol: f64) -> bool {
        let ah = self.conjugate_transpose();
        let prod = ah.mul(self);
        let id = ComplexMat2::identity();
        (prod.a - id.a).norm() < tol
            && (prod.b - id.b).norm() < tol
            && (prod.c - id.c).norm() < tol
            && (prod.d - id.d).norm() < tol
    }
}
/// A Möbius (fractional linear) transformation `f(z) = (az + b) / (cz + d)`.
///
/// These transformations are the automorphisms of the Riemann sphere and map
/// circles (and lines) to circles (or lines).
#[derive(Debug, Clone, Copy)]
pub struct MobiusTransform {
    /// Numerator coefficient of `z`.
    pub a: Complex,
    /// Numerator constant.
    pub b: Complex,
    /// Denominator coefficient of `z`.
    pub c: Complex,
    /// Denominator constant.
    pub d: Complex,
}
impl MobiusTransform {
    /// Create a new Möbius transform `(az+b)/(cz+d)`.
    pub fn new(a: Complex, b: Complex, c: Complex, d: Complex) -> Self {
        Self { a, b, c, d }
    }
    /// The identity transform `z → z` (a=1, b=0, c=0, d=1).
    pub fn identity() -> Self {
        Self::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::one(),
        )
    }
    /// Translation `z → z + t`.
    pub fn translation(t: Complex) -> Self {
        Self::new(Complex::one(), t, Complex::zero(), Complex::one())
    }
    /// Rotation/dilation by `r` (multiply by r): `z → r*z`.
    pub fn rotation_dilation(r: Complex) -> Self {
        Self::new(r, Complex::zero(), Complex::zero(), Complex::one())
    }
    /// Inversion: `z → 1/z` (a=0, b=1, c=1, d=0).
    pub fn inversion() -> Self {
        Self::new(
            Complex::zero(),
            Complex::one(),
            Complex::one(),
            Complex::zero(),
        )
    }
    /// Apply the transform to `z`.  Returns `None` if `cz + d = 0` (pole).
    pub fn apply(&self, z: Complex) -> Complex {
        let denom = self.c * z + self.d;
        (self.a * z + self.b) / denom
    }
    /// Compose `self` with `other`: `(self ∘ other)(z) = self(other(z))`.
    ///
    /// This corresponds to matrix multiplication of `[[a,b\],[c,d]]`.
    pub fn compose(&self, other: &MobiusTransform) -> MobiusTransform {
        MobiusTransform::new(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d,
        )
    }
    /// Compute the inverse Möbius transform, or `None` if `ad - bc = 0`.
    pub fn inverse(&self) -> Option<MobiusTransform> {
        let det = self.a * self.d - self.b * self.c;
        if det.norm_sq() < 1e-24 {
            return None;
        }
        Some(MobiusTransform::new(self.d, -self.b, -self.c, self.a))
    }
    /// Returns the two fixed points of the transform (solve `f(z) = z`).
    ///
    /// Returns `None` for the identity transform or degenerate cases.
    pub fn fixed_points(&self) -> Option<(Complex, Complex)> {
        let c = self.c;
        if c.norm_sq() < 1e-24 {
            let denom = self.a - self.d;
            if denom.norm_sq() < 1e-24 {
                return None;
            }
            let z0 = self.b / denom;
            return Some((z0, z0));
        }
        let disc = (self.d - self.a) * (self.d - self.a) + Complex::new(4.0, 0.0) * self.b * self.c;
        let sqrt_disc = disc.sqrt();
        let two_c = c * Complex::new(2.0, 0.0);
        let z1 = (self.a - self.d + sqrt_disc) / two_c;
        let z2 = (self.a - self.d - sqrt_disc) / two_c;
        Some((z1, z2))
    }
}
