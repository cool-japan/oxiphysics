//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{Complex, ComplexMat2};
use std::f64::consts::PI;

/// Check whether a function `f: ℂ → ℂ` satisfies the Cauchy-Riemann equations
/// at point `z` using finite differences with step `h`.
///
/// The C-R equations are: ∂u/∂x = ∂v/∂y  and  ∂u/∂y = -∂v/∂x,
/// where f(z) = u + iv.
pub fn cauchy_riemann_check(f: &impl Fn(Complex) -> Complex, z: Complex, h: f64) -> bool {
    let fxp = f(Complex::new(z.re + h, z.im));
    let fxm = f(Complex::new(z.re - h, z.im));
    let fyp = f(Complex::new(z.re, z.im + h));
    let fym = f(Complex::new(z.re, z.im - h));
    let du_dx = (fxp.re - fxm.re) / (2.0 * h);
    let du_dy = (fyp.re - fym.re) / (2.0 * h);
    let dv_dx = (fxp.im - fxm.im) / (2.0 * h);
    let dv_dy = (fyp.im - fym.im) / (2.0 * h);
    let tol = h.sqrt().max(1e-4);
    (du_dx - dv_dy).abs() < tol && (du_dy + dv_dx).abs() < tol
}
/// Joukowski conformal map: `w = z + 1/z`.
///
/// Maps circles to airfoil-like shapes (Joukowski profiles).
pub fn joukowski(z: Complex) -> Complex {
    z + Complex::one() / z
}
/// Schwarz-Christoffel map approximation for a rectangle.
///
/// Maps the upper half-plane to a rectangle-like region.
/// Uses the simplified integral `∫ dz / sqrt((1-z²)(1-k²z²))`.
/// Computed numerically using midpoint quadrature.
pub fn schwarz_christoffel_rectangle(z: Complex, k: f64, n_steps: usize) -> Complex {
    let mut sum = Complex::zero();
    let n = n_steps as f64;
    for i in 0..n_steps {
        let t_mid = z * Complex::new((i as f64 + 0.5) / n, 0.0);
        let dz = z * Complex::new(1.0 / n, 0.0);
        let t2 = t_mid * t_mid;
        let k2 = Complex::new(k * k, 0.0);
        let factor = (Complex::one() - t2) * (Complex::one() - k2 * t2);
        let integrand = Complex::one() / factor.sqrt();
        sum = sum + integrand * dz;
    }
    sum
}
/// Exponential conformal map: `w = exp(z)`.
///
/// Maps horizontal strips of height 2π to annular sectors.
pub fn exponential_conformal(z: Complex) -> Complex {
    z.exp()
}
/// Logarithmic conformal map: `w = ln(z)` (principal branch).
///
/// Inverse of the exponential map; maps annular sectors to strips.
pub fn logarithmic_conformal(z: Complex) -> Complex {
    z.ln()
}
/// Cayley transform: maps the upper half-plane to the unit disk.
///
/// `w = (z - i) / (z + i)`
pub fn cayley_transform(z: Complex) -> Complex {
    let i = Complex::new(0.0, 1.0);
    (z - i) / (z + i)
}
/// Inverse Cayley transform: maps the unit disk to the upper half-plane.
///
/// `w = i (1 + z) / (1 - z)`
pub fn cayley_transform_inverse(z: Complex) -> Complex {
    let i = Complex::new(0.0, 1.0);
    i * (Complex::one() + z) / (Complex::one() - z)
}
/// Cooley-Tukey Radix-2 FFT (in-place, decimation-in-time).
///
/// Input length must be a power of two. If not, the input is zero-padded.
pub fn fft_radix2(xs: &[Complex]) -> Vec<Complex> {
    let n = xs.len().next_power_of_two();
    let mut a: Vec<Complex> = xs.to_vec();
    a.resize(n, Complex::zero());
    fft_inplace(&mut a, false);
    a
}
/// Inverse FFT (Cooley-Tukey Radix-2).
pub fn ifft_radix2(xs: &[Complex]) -> Vec<Complex> {
    let n = xs.len().next_power_of_two();
    let mut a: Vec<Complex> = xs.to_vec();
    a.resize(n, Complex::zero());
    fft_inplace(&mut a, true);
    let inv_n = 1.0 / n as f64;
    a.iter_mut()
        .for_each(|z| *z = Complex::new(z.re * inv_n, z.im * inv_n));
    a
}
pub(super) fn fft_inplace(a: &mut [Complex], inverse: bool) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2usize;
    while len <= n {
        use std::f64::consts::PI;
        let sign = if inverse { 1.0 } else { -1.0 };
        let ang = sign * 2.0 * PI / len as f64;
        let wlen = Complex::from_polar(1.0, ang);
        let mut i = 0;
        while i < n {
            let mut w = Complex::one();
            for k in 0..(len / 2) {
                let u = a[i + k];
                let v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w = w * wlen;
            }
            i += len;
        }
        len <<= 1;
    }
}
/// Compute the DFT of a real-valued signal (returns N/2+1 complex frequencies).
pub fn rfft(xs: &[f64]) -> Vec<Complex> {
    let complex_in: Vec<Complex> = xs.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let out = fft_radix2(&complex_in);
    let half = out.len() / 2 + 1;
    out[..half].to_vec()
}
/// Taylor series approximation of exp(z) = Σ z^k / k! up to `terms` terms.
pub fn complex_exp_taylor(z: Complex, terms: usize) -> Complex {
    let mut result = Complex::one();
    let mut term = Complex::one();
    for k in 1..=terms {
        term = term * z * Complex::new(1.0 / k as f64, 0.0);
        result = result + term;
    }
    result
}
/// Taylor series approximation of sin(z) = Σ (-1)^k z^(2k+1) / (2k+1)!
pub fn complex_sin_taylor(z: Complex, terms: usize) -> Complex {
    let mut result = Complex::zero();
    let mut term = z;
    let z2 = z * z;
    for k in 0..terms {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        result = result + term * Complex::new(sign, 0.0);
        let denom = (2 * k + 2) as f64 * (2 * k + 3) as f64;
        term = term * z2 * Complex::new(1.0 / denom, 0.0);
    }
    result
}
/// Numerically differentiate `f` at `z` using the complex-step method.
///
/// For analytic functions, `f'(z) ≈ Im[f(z + ih)] / h` with high accuracy.
pub fn complex_step_derivative(f: &impl Fn(Complex) -> Complex, z: Complex, h: f64) -> Complex {
    let fz_ih = f(Complex::new(z.re, z.im + h));
    let fz = f(z);
    let dz = Complex::new(0.0, h);
    (fz_ih - fz) / dz
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex::Complex;
    use crate::complex::ComplexMat2;

    use crate::complex::MobiusTransform;
    use crate::complex::QuatAlgebra;

    use crate::complex::fft_radix2;

    use crate::complex::joukowski;

    use std::f64::consts::SQRT_2;
    pub(super) const EPS: f64 = 1e-10;
    pub(super) const EPS_6: f64 = 1e-6;
    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }
    fn approx_quat(a: QuatAlgebra, b: QuatAlgebra, eps: f64) -> bool {
        approx(a.w, b.w, eps)
            && approx(a.x, b.x, eps)
            && approx(a.y, b.y, eps)
            && approx(a.z, b.z, eps)
    }
    #[test]
    fn complex_add() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, -1.0);
        let c = a + b;
        assert!(approx(c.re, 4.0, EPS) && approx(c.im, 1.0, EPS));
    }
    #[test]
    fn complex_sub() {
        let a = Complex::new(5.0, 3.0);
        let b = Complex::new(2.0, 1.0);
        let c = a - b;
        assert!(approx(c.re, 3.0, EPS) && approx(c.im, 2.0, EPS));
    }
    #[test]
    fn complex_mul() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a * b;
        assert!(approx(c.re, -5.0, EPS) && approx(c.im, 10.0, EPS));
    }
    #[test]
    fn complex_div() {
        let a = Complex::new(1.0, 2.0);
        let c = a / a;
        assert!(approx(c.re, 1.0, EPS) && approx(c.im, 0.0, EPS));
    }
    #[test]
    fn complex_neg() {
        let a = Complex::new(3.0, -4.0);
        let c = -a;
        assert!(approx(c.re, -3.0, EPS) && approx(c.im, 4.0, EPS));
    }
    #[test]
    fn complex_conj() {
        let a = Complex::new(3.0, 4.0);
        let c = a.conj();
        assert!(approx(c.re, 3.0, EPS) && approx(c.im, -4.0, EPS));
    }
    #[test]
    fn complex_norm() {
        let a = Complex::new(3.0, 4.0);
        assert!(approx(a.norm(), 5.0, EPS));
    }
    #[test]
    fn complex_norm_sq() {
        let a = Complex::new(3.0, 4.0);
        assert!(approx(a.norm_sq(), 25.0, EPS));
    }
    #[test]
    fn complex_arg() {
        let a = Complex::new(0.0, 1.0);
        assert!(approx(a.arg(), PI / 2.0, EPS));
    }
    #[test]
    fn complex_from_polar() {
        let c = Complex::from_polar(1.0, PI / 2.0);
        assert!(approx(c.re, 0.0, EPS_6) && approx(c.im, 1.0, EPS_6));
    }
    #[test]
    fn complex_exp_euler() {
        let c = Complex::new(0.0, PI).exp();
        assert!(approx(c.re, -1.0, EPS_6) && approx(c.im, 0.0, EPS_6));
    }
    #[test]
    fn complex_ln_then_exp() {
        let z = Complex::new(2.0, 3.0);
        let recovered = z.ln().exp();
        assert!(approx(recovered.re, z.re, EPS_6) && approx(recovered.im, z.im, EPS_6));
    }
    #[test]
    fn complex_sqrt_of_minus_one() {
        let z = Complex::new(-1.0, 0.0);
        let s = z.sqrt();
        assert!(approx(s.re, 0.0, EPS_6) && approx(s.im, 1.0, EPS_6));
    }
    #[test]
    fn complex_sqrt_of_two() {
        let z = Complex::new(2.0, 0.0);
        let s = z.sqrt();
        assert!(approx(s.re, SQRT_2, EPS_6) && approx(s.im, 0.0, EPS_6));
    }
    #[test]
    fn complex_pow_square() {
        let z = Complex::new(1.0, 1.0);
        let p = z.pow_f64(2.0);
        assert!(approx(p.re, 0.0, EPS_6) && approx(p.im, 2.0, EPS_6));
    }
    #[test]
    fn complex_sin_real() {
        let z = Complex::new(PI / 6.0, 0.0);
        let s = z.sin();
        assert!(approx(s.re, 0.5, EPS_6) && approx(s.im, 0.0, EPS_6));
    }
    #[test]
    fn complex_cos_real() {
        let z = Complex::new(0.0, 0.0);
        let c = z.cos();
        assert!(approx(c.re, 1.0, EPS) && approx(c.im, 0.0, EPS));
    }
    #[test]
    fn complex_pythagorean_identity() {
        let z = Complex::new(1.0, 0.5);
        let sin2 = z.sin() * z.sin();
        let cos2 = z.cos() * z.cos();
        let sum = sin2 + cos2;
        assert!(approx(sum.re, 1.0, EPS_6) && approx(sum.im, 0.0, EPS_6));
    }
    #[test]
    fn complex_zero_one_i() {
        assert_eq!(Complex::zero(), Complex::new(0.0, 0.0));
        assert_eq!(Complex::one(), Complex::new(1.0, 0.0));
        assert_eq!(Complex::i(), Complex::new(0.0, 1.0));
    }
    #[test]
    fn complex_mul_associative() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, -1.0);
        let c = Complex::new(-2.0, 5.0);
        let lhs = (a * b) * c;
        let rhs = a * (b * c);
        assert!(approx(lhs.re, rhs.re, EPS_6) && approx(lhs.im, rhs.im, EPS_6));
    }
    #[test]
    fn quat_identity_mul() {
        let q = QuatAlgebra::new(0.5, 0.5, 0.5, 0.5);
        let id = QuatAlgebra::identity();
        assert!(approx_quat(q * id, q, EPS));
        assert!(approx_quat(id * q, q, EPS));
    }
    #[test]
    fn quat_mul_not_commutative() {
        let i = QuatAlgebra::new(0.0, 1.0, 0.0, 0.0);
        let j = QuatAlgebra::new(0.0, 0.0, 1.0, 0.0);
        let k = QuatAlgebra::new(0.0, 0.0, 0.0, 1.0);
        let ij = i * j;
        assert!(approx_quat(ij, k, EPS));
        let ji = j * i;
        assert!(approx_quat(ji, QuatAlgebra::new(0.0, 0.0, 0.0, -1.0), EPS));
    }
    #[test]
    fn quat_mul_associative() {
        let a = QuatAlgebra::new(1.0, 2.0, -1.0, 0.5).normalize();
        let b = QuatAlgebra::new(0.5, -0.5, 1.0, 2.0).normalize();
        let c = QuatAlgebra::new(-1.0, 1.0, 1.0, -1.0).normalize();
        let lhs = (a * b) * c;
        let rhs = a * (b * c);
        assert!(approx_quat(lhs, rhs, EPS_6));
    }
    #[test]
    fn quat_norm_of_unit() {
        let q = QuatAlgebra::from_axis_angle([0.0, 1.0, 0.0], PI / 3.0);
        assert!(approx(q.norm(), 1.0, EPS_6));
    }
    #[test]
    fn quat_inverse() {
        let q = QuatAlgebra::new(1.0, 2.0, 3.0, 4.0);
        let qi = q.inverse();
        let prod = q * qi;
        assert!(approx_quat(prod, QuatAlgebra::identity(), EPS_6));
    }
    #[test]
    fn quat_conj_times_self_is_norm_sq() {
        let q = QuatAlgebra::new(1.0, 2.0, 3.0, 4.0);
        let qc = q.conj();
        let prod = q * qc;
        assert!(approx(prod.w, q.norm_sq(), EPS_6));
        assert!(approx(prod.x, 0.0, EPS_6));
        assert!(approx(prod.y, 0.0, EPS_6));
        assert!(approx(prod.z, 0.0, EPS_6));
    }
    #[test]
    fn quat_rotate_90_around_z() {
        let q = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let v = [1.0_f64, 0.0, 0.0];
        let r = q.rotate_vector(v);
        assert!(approx(r[0], 0.0, EPS_6));
        assert!(approx(r[1], 1.0, EPS_6));
        assert!(approx(r[2], 0.0, EPS_6));
    }
    #[test]
    fn quat_rotation_composition() {
        let q90 = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let q180 = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI);
        let q_composed = q90 * q90;
        let v = [1.0_f64, 0.0, 0.0];
        let r1 = q_composed.rotate_vector(v);
        let r2 = q180.rotate_vector(v);
        assert!(approx(r1[0], r2[0], EPS_6));
        assert!(approx(r1[1], r2[1], EPS_6));
        assert!(approx(r1[2], r2[2], EPS_6));
    }
    #[test]
    fn quat_to_rotation_matrix_orthogonal() {
        let q = QuatAlgebra::from_axis_angle([1.0, 1.0, 0.0], PI / 4.0);
        let r = q.to_rotation_matrix();
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| r[i][k] * r[j][k]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    approx(dot, expected, EPS_6),
                    "R R^T [{i}][{j}] = {dot}, expected {expected}"
                );
            }
        }
    }
    #[test]
    fn quat_exp_ln_roundtrip() {
        let q = QuatAlgebra::from_axis_angle([1.0, 0.0, 0.0], PI / 4.0);
        let recovered = q.ln().exp();
        assert!(approx_quat(recovered, q, EPS_6));
    }
    #[test]
    fn quat_pow_half_is_half_angle() {
        let q = QuatAlgebra::from_axis_angle([0.0, 1.0, 0.0], PI / 2.0);
        let q_half = q.pow(0.5);
        let expected = QuatAlgebra::from_axis_angle([0.0, 1.0, 0.0], PI / 4.0);
        assert!(approx_quat(q_half, expected, EPS_6));
    }
    #[test]
    fn quat_slerp_endpoints() {
        let q0 = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], 0.0);
        let q1 = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let s0 = q0.slerp(&q1, 0.0);
        let s1 = q0.slerp(&q1, 1.0);
        assert!(approx_quat(s0, q0, EPS_6));
        assert!(approx_quat(s1, q1, EPS_6));
    }
    #[test]
    fn quat_slerp_midpoint() {
        let q0 = QuatAlgebra::identity();
        let q1 = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let mid = q0.slerp(&q1, 0.5);
        let expected = QuatAlgebra::from_axis_angle([0.0, 0.0, 1.0], PI / 4.0);
        assert!(approx_quat(mid, expected, EPS_6));
    }
    #[test]
    fn quat_dot_self_is_norm_sq() {
        let q = QuatAlgebra::new(1.0, 2.0, 3.0, 4.0);
        assert!(approx(q.dot(&q), q.norm_sq(), EPS));
    }
    #[test]
    fn quat_normalize_gives_unit() {
        let q = QuatAlgebra::new(1.0, 2.0, 3.0, 4.0);
        assert!(approx(q.normalize().norm(), 1.0, EPS_6));
    }
    #[test]
    fn quat_zero() {
        let z = QuatAlgebra::zero();
        assert!(approx(z.norm(), 0.0, EPS));
    }
    #[test]
    fn test_complex_mat2_identity_mul() {
        let id = ComplexMat2::identity();
        let m = ComplexMat2::new(
            Complex::new(1.0, 2.0),
            Complex::new(3.0, -1.0),
            Complex::new(-2.0, 0.5),
            Complex::new(4.0, 1.0),
        );
        let r = id.mul(&m);
        assert!(approx(r.a.re, m.a.re, EPS_6) && approx(r.d.im, m.d.im, EPS_6));
    }
    #[test]
    fn test_complex_mat2_det_identity() {
        let id = ComplexMat2::identity();
        let det = id.det();
        assert!(approx(det.re, 1.0, EPS) && approx(det.im, 0.0, EPS));
    }
    #[test]
    fn test_complex_mat2_inverse() {
        let m = ComplexMat2::new(
            Complex::new(2.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
        );
        if let Some(inv) = m.inverse() {
            let prod = m.mul(&inv);
            let id = ComplexMat2::identity();
            assert!(approx(prod.a.re, id.a.re, EPS_6));
            assert!(approx(prod.d.re, id.d.re, EPS_6));
        }
    }
    #[test]
    fn test_complex_mat2_eigenvalues() {
        let m = ComplexMat2::new(
            Complex::new(3.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
        );
        let (l1, l2) = m.eigenvalues();
        let trace = 5.0_f64;
        let product = 5.0_f64;
        let disc = (trace * trace - 4.0 * product).sqrt();
        let expected_l1 = (trace + disc) / 2.0;
        let expected_l2 = (trace - disc) / 2.0;
        assert!(approx(l1.re, expected_l1, 1e-6), "l1={}", l1.re);
        assert!(approx(l2.re, expected_l2, 1e-6), "l2={}", l2.re);
    }
    #[test]
    fn test_cauchy_riemann_exp() {
        let z = Complex::new(1.0, 0.5);
        let satisfied = cauchy_riemann_check(&|z| z.exp(), z, 1e-6);
        assert!(satisfied, "exp(z) should satisfy Cauchy-Riemann");
    }
    #[test]
    fn test_cauchy_riemann_z_squared() {
        let z = Complex::new(2.0, 1.0);
        let satisfied = cauchy_riemann_check(&|z| z * z, z, 1e-5);
        assert!(satisfied, "z^2 should satisfy Cauchy-Riemann");
    }
    #[test]
    fn test_cauchy_riemann_conj_fails() {
        let z = Complex::new(1.0, 1.0);
        let satisfied = cauchy_riemann_check(&|z: Complex| z.conj(), z, 1e-5);
        assert!(!satisfied, "conj(z) should NOT satisfy Cauchy-Riemann");
    }
    #[test]
    fn test_mobius_identity() {
        let id = MobiusTransform::identity();
        let z = Complex::new(2.0, 3.0);
        let fz = id.apply(z);
        assert!(approx(fz.re, z.re, EPS_6) && approx(fz.im, z.im, EPS_6));
    }
    #[test]
    fn test_mobius_inversion_z() {
        let f = MobiusTransform::new(
            Complex::zero(),
            Complex::one(),
            Complex::one(),
            Complex::zero(),
        );
        let z = Complex::new(2.0, 0.0);
        let fz = f.apply(z);
        assert!(approx(fz.re, 0.5, EPS_6), "1/2 = 0.5, got {}", fz.re);
    }
    #[test]
    fn test_mobius_compose_with_inverse() {
        let f = MobiusTransform::new(
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        );
        if let Some(inv) = f.inverse() {
            let composed = f.compose(&inv);
            let z = Complex::new(1.5, 0.7);
            let fz = composed.apply(z);
            let id_z = MobiusTransform::identity().apply(z);
            assert!(approx(fz.re, id_z.re, 1e-8) && approx(fz.im, id_z.im, 1e-8));
        }
    }
    #[test]
    fn test_joukowski_symmetry() {
        let z = Complex::new(2.0, 1.0);
        let w1 = joukowski(z);
        let w2 = joukowski(z.conj()).conj();
        assert!(approx(w1.re, w2.re, EPS_6) && approx(w1.im, w2.im, EPS_6));
    }
    #[test]
    fn test_joukowski_real_axis_maps_real() {
        let z = Complex::new(3.0, 0.0);
        let w = joukowski(z);
        assert!(w.im.abs() < EPS_6, "im={}", w.im);
    }
    #[test]
    fn test_fft_size_1() {
        let xs = vec![Complex::new(5.0, 3.0)];
        let out = fft_radix2(&xs);
        assert_eq!(out.len(), 1);
        assert!(approx(out[0].re, 5.0, EPS_6) && approx(out[0].im, 3.0, EPS_6));
    }
    #[test]
    fn test_fft_size_4_dc() {
        let val = Complex::new(1.0, 0.0);
        let xs = vec![val; 4];
        let out = fft_radix2(&xs);
        assert!(approx(out[0].re, 4.0, EPS_6), "DC={}", out[0].re);
        assert!(approx(out[1].re, 0.0, EPS_6) && approx(out[1].im, 0.0, EPS_6));
    }
    #[test]
    fn test_fft_parseval() {
        let xs = vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 1.0),
            Complex::new(-1.0, 0.0),
            Complex::new(0.0, -1.0),
        ];
        let out = fft_radix2(&xs);
        let energy_time: f64 = xs.iter().map(|z| z.norm_sq()).sum();
        let energy_freq: f64 = out.iter().map(|z| z.norm_sq()).sum::<f64>() / xs.len() as f64;
        assert!(
            approx(energy_time, energy_freq, 1e-9),
            "Parseval: {} vs {}",
            energy_time,
            energy_freq
        );
    }
    #[test]
    fn test_complex_power_series_exp() {
        let z = Complex::new(0.5, 0.3);
        let approx_exp = complex_exp_taylor(z, 15);
        let exact = z.exp();
        assert!(
            approx(approx_exp.re, exact.re, 1e-8),
            "re diff={}",
            (approx_exp.re - exact.re).abs()
        );
        assert!(
            approx(approx_exp.im, exact.im, 1e-8),
            "im diff={}",
            (approx_exp.im - exact.im).abs()
        );
    }
}
/// Compute the N×N DFT matrix W where W\[j\]\[k\] = exp(-2πi·j·k/N).
///
/// The matrix form of the DFT: X = W·x.
#[allow(dead_code)]
pub fn dft_matrix(n: usize) -> Vec<Vec<Complex>> {
    let mut mat = vec![vec![Complex::zero(); n]; n];
    for j in 0..n {
        for k in 0..n {
            let angle = -2.0 * PI * (j * k) as f64 / n as f64;
            mat[j][k] = Complex::from_polar(1.0, angle);
        }
    }
    mat
}
/// Apply the DFT matrix to a vector (O(N²) naive DFT).
///
/// Returns X\[k\] = Σ x\[n\] · exp(-2πi·k·n/N).
#[allow(dead_code)]
pub fn naive_dft(xs: &[Complex]) -> Vec<Complex> {
    let n = xs.len();
    (0..n)
        .map(|k| {
            let mut sum = Complex::zero();
            for (j, &x) in xs.iter().enumerate() {
                let angle = -2.0 * PI * (k * j) as f64 / n as f64;
                sum = sum + x * Complex::from_polar(1.0, angle);
            }
            sum
        })
        .collect()
}
/// Returns the N-th roots of unity: exp(2πi·k/N) for k=0,…,N-1.
#[allow(dead_code)]
pub fn roots_of_unity(n: usize) -> Vec<Complex> {
    (0..n)
        .map(|k| Complex::from_polar(1.0, 2.0 * PI * k as f64 / n as f64))
        .collect()
}
/// Returns the primitive N-th root of unity: exp(2πi/N).
#[allow(dead_code)]
pub fn primitive_root_of_unity(n: usize) -> Complex {
    Complex::from_polar(1.0, 2.0 * PI / n as f64)
}
/// Numerically integrate `f(z)` along a circular contour of radius `r`
/// centered at `center`, using `n_points` quadrature points.
///
/// Returns ∮_C f(z) dz ≈ Σ f(z_k) · Δz_k.
#[allow(dead_code)]
pub fn contour_integral_circle(
    f: &impl Fn(Complex) -> Complex,
    center: Complex,
    r: f64,
    n_points: usize,
) -> Complex {
    let n = n_points.max(4);
    let mut sum = Complex::zero();
    for k in 0..n {
        let theta = 2.0 * PI * k as f64 / n as f64;
        let theta_next = 2.0 * PI * (k + 1) as f64 / n as f64;
        let z_mid = center + Complex::from_polar(r, (theta + theta_next) / 2.0);
        let dz = Complex::from_polar(r, theta_next) - Complex::from_polar(r, theta);
        sum = sum + f(z_mid) * dz;
    }
    sum
}
/// Numerically compute the residue of `f` at `pole` using the contour integral.
///
/// For a simple pole: Res(f, z0) = (1/(2πi)) ∮_{|z-z0|=r} f(z) dz.
#[allow(dead_code)]
pub fn residue_by_contour(
    f: &impl Fn(Complex) -> Complex,
    pole: Complex,
    r: f64,
    n_points: usize,
) -> Complex {
    let integral = contour_integral_circle(f, pole, r, n_points);
    let two_pi_i = Complex::new(0.0, 2.0 * PI);
    integral / two_pi_i
}
/// Find all roots of a monic polynomial using the Durand-Kerner (Weierstrass) method.
///
/// The polynomial is given as `coeffs[0] + coeffs[1]*z + ... + coeffs[n]*z^n`
/// with the leading coefficient being monic (implicit 1).
///
/// `coeffs` should have length `n` (for a degree-n polynomial with leading coeff 1).
/// Returns `n` approximate roots.
#[allow(dead_code)]
pub fn durand_kerner_roots(coeffs: &[Complex], max_iter: usize) -> Vec<Complex> {
    let n = coeffs.len();
    if n == 0 {
        return vec![];
    }
    let omega = Complex::from_polar(0.4, 2.0 * PI / n as f64);
    let mut roots: Vec<Complex> = (0..n)
        .map(|k| omega.pow_f64(k as f64) * Complex::new(0.4, 0.0))
        .collect();
    let r0 = (coeffs[0].norm() + 1.0).powf(1.0 / n as f64).max(0.5);
    for k in 0..n {
        let angle = 2.0 * PI * k as f64 / n as f64 + 0.1;
        roots[k] = Complex::from_polar(r0, angle);
    }
    let eval_poly = |z: Complex| -> Complex {
        let mut result = Complex::one();
        let mut p = Complex::zero();
        for c in coeffs.iter() {
            p = p + *c * result;
            result = result * z;
        }
        p + result
    };
    for _ in 0..max_iter {
        let old_roots = roots.clone();
        for i in 0..n {
            let pz = eval_poly(old_roots[i]);
            let mut denom = Complex::one();
            for j in 0..n {
                if i != j {
                    denom = denom * (old_roots[i] - old_roots[j]);
                }
            }
            if denom.norm_sq() > 1e-30 {
                roots[i] = old_roots[i] - pz / denom;
            }
        }
    }
    roots
}
/// Solve the complex linear system `A·x = b` using Gaussian elimination
/// with partial pivoting.
///
/// `a` is an N×N matrix in row-major format (Vec of rows).
/// Returns the solution vector `x`, or `None` if the system is singular.
#[allow(dead_code)]
pub fn complex_gaussian_elimination(a: &[Vec<Complex>], b: &[Complex]) -> Option<Vec<Complex>> {
    let n = b.len();
    assert_eq!(a.len(), n, "complex_gaussian_elimination: A must be n×n");
    for row in a {
        assert_eq!(row.len(), n, "complex_gaussian_elimination: A must be n×n");
    }
    let mut mat: Vec<Vec<Complex>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, &bi)| {
            let mut r = row.clone();
            r.push(bi);
            r
        })
        .collect();
    for col in 0..n {
        let (pivot_row, _) = (col..n)
            .map(|r| (r, mat[r][col].norm_sq()))
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;
        if mat[pivot_row][col].norm_sq() < 1e-30 {
            return None;
        }
        mat.swap(col, pivot_row);
        let pivot = mat[col][col];
        for j in col..=n {
            let v = mat[col][j];
            mat[col][j] = v / pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = mat[row][col];
                for j in col..=n {
                    let v = mat[col][j] * factor;
                    let old = mat[row][j];
                    mat[row][j] = old - v;
                }
            }
        }
    }
    Some((0..n).map(|i| mat[i][n]).collect())
}
#[cfg(test)]
mod tests_new_complex {

    use crate::complex::Complex;

    use crate::complex::ComplexMatN;

    use crate::complex::cayley_transform;
    use crate::complex::cayley_transform_inverse;
    use crate::complex::complex_gaussian_elimination;
    use crate::complex::complex_sin_taylor;
    use crate::complex::complex_step_derivative;
    use crate::complex::contour_integral_circle;

    use crate::complex::dft_matrix;
    use crate::complex::durand_kerner_roots;

    use crate::complex::fft_radix2;
    use crate::complex::functions::PI;

    use crate::complex::ifft_radix2;
    use crate::complex::joukowski;

    use crate::complex::naive_dft;

    use crate::complex::primitive_root_of_unity;

    use crate::complex::residue_by_contour;
    use crate::complex::rfft;
    use crate::complex::roots_of_unity;

    use crate::complex::schwarz_christoffel_rectangle;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }
    fn capprox(a: Complex, b: Complex, eps: f64) -> bool {
        approx(a.re, b.re, eps) && approx(a.im, b.im, eps)
    }
    #[test]
    fn test_dft_matrix_size() {
        let m = dft_matrix(4);
        assert_eq!(m.len(), 4);
        assert_eq!(m[0].len(), 4);
    }
    #[test]
    fn test_dft_matrix_unit_entries() {
        let m = dft_matrix(4);
        for row in &m {
            for entry in row {
                assert!(
                    approx(entry.norm(), 1.0, 1e-10),
                    "DFT matrix entries should be on unit circle"
                );
            }
        }
    }
    #[test]
    fn test_naive_dft_matches_fft() {
        let xs = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];
        let fft_out = fft_radix2(&xs);
        let dft_out = naive_dft(&xs);
        for (f, d) in fft_out.iter().zip(dft_out.iter()) {
            assert!(
                approx(f.re, d.re, 1e-8),
                "FFT re={} vs DFT re={}",
                f.re,
                d.re
            );
            assert!(
                approx(f.im, d.im, 1e-8),
                "FFT im={} vs DFT im={}",
                f.im,
                d.im
            );
        }
    }
    #[test]
    fn test_naive_dft_dc_component() {
        let xs = vec![Complex::new(1.0, 0.0); 4];
        let out = naive_dft(&xs);
        assert!(
            approx(out[0].re, 4.0, 1e-10),
            "DC component should be N*value"
        );
        assert!(
            approx(out[1].norm(), 0.0, 1e-8),
            "non-DC should be zero for constant input"
        );
    }
    #[test]
    fn test_complex_mat_n_identity_mul() {
        let id = ComplexMatN::identity(3);
        let mut m = ComplexMatN::zeros(3);
        for i in 0..3 {
            for j in 0..3 {
                m.set(i, j, Complex::new((i * 3 + j) as f64, 0.0));
            }
        }
        let result = id.mul(&m);
        for i in 0..3 {
            for j in 0..3 {
                assert!(approx(result.get(i, j).re, m.get(i, j).re, 1e-10));
            }
        }
    }
    #[test]
    fn test_complex_mat_n_trace_identity() {
        let id = ComplexMatN::identity(4);
        let tr = id.trace();
        assert!(approx(tr.re, 4.0, 1e-10) && approx(tr.im, 0.0, 1e-10));
    }
    #[test]
    fn test_complex_mat_n_frobenius_identity() {
        let id = ComplexMatN::identity(3);
        assert!(approx(id.frobenius_norm(), 3.0_f64.sqrt(), 1e-10));
    }
    #[test]
    fn test_complex_mat_n_conj_transpose() {
        let mut m = ComplexMatN::zeros(2);
        m.set(0, 0, Complex::new(1.0, 2.0));
        m.set(0, 1, Complex::new(3.0, -1.0));
        m.set(1, 0, Complex::new(0.0, 4.0));
        m.set(1, 1, Complex::new(5.0, 0.0));
        let mh = m.conj_transpose();
        assert!(capprox(mh.get(0, 1), m.get(1, 0).conj(), 1e-10));
    }
    #[test]
    fn test_complex_mat_n_apply_identity() {
        let id = ComplexMatN::identity(3);
        let v = vec![
            Complex::new(1.0, 2.0),
            Complex::new(-1.0, 0.5),
            Complex::new(0.0, 3.0),
        ];
        let result = id.apply(&v);
        for (a, b) in result.iter().zip(v.iter()) {
            assert!(capprox(*a, *b, 1e-10));
        }
    }
    #[test]
    fn test_roots_of_unity_count() {
        let roots = roots_of_unity(6);
        assert_eq!(roots.len(), 6);
    }
    #[test]
    fn test_roots_of_unity_on_unit_circle() {
        let roots = roots_of_unity(8);
        for r in &roots {
            assert!(
                approx(r.norm(), 1.0, 1e-10),
                "root should be on unit circle, norm={}",
                r.norm()
            );
        }
    }
    #[test]
    fn test_roots_of_unity_sum_to_zero() {
        let roots = roots_of_unity(5);
        let sum = roots.iter().fold(Complex::zero(), |acc, &r| acc + r);
        assert!(
            approx(sum.norm(), 0.0, 1e-10),
            "sum of roots of unity should be 0, norm={}",
            sum.norm()
        );
    }
    #[test]
    fn test_roots_of_unity_first_is_one() {
        let roots = roots_of_unity(4);
        assert!(approx(roots[0].re, 1.0, 1e-10) && approx(roots[0].im, 0.0, 1e-10));
    }
    #[test]
    fn test_primitive_root_of_unity_order() {
        let n = 6usize;
        let omega = primitive_root_of_unity(n);
        let mut w = omega;
        for _ in 1..n {
            w = w * omega;
        }
        assert!(
            approx(w.re, 1.0, 1e-10) && approx(w.im, 0.0, 1e-10),
            "omega^n should be 1, got ({}, {})",
            w.re,
            w.im
        );
    }
    #[test]
    fn test_contour_integral_analytic_function_zero() {
        let center = Complex::new(0.0, 0.0);
        let integral = contour_integral_circle(&|z: Complex| z * z, center, 1.0, 1000);
        assert!(
            approx(integral.norm(), 0.0, 1e-6),
            "integral of analytic fn should be 0, got {}",
            integral.norm()
        );
    }
    #[test]
    fn test_contour_integral_1_over_z() {
        let center = Complex::new(0.0, 0.0);
        let integral =
            contour_integral_circle(&|z: Complex| Complex::one() / z, center, 1.0, 10000);
        assert!(
            approx(integral.re, 0.0, 1e-3),
            "real part should be 0, got {}",
            integral.re
        );
        assert!(
            approx(integral.im, 2.0 * PI, 1e-3),
            "imag part should be 2π, got {}",
            integral.im
        );
    }
    #[test]
    fn test_residue_1_over_z() {
        let center = Complex::new(0.0, 0.0);
        let res = residue_by_contour(&|z: Complex| Complex::one() / z, center, 1.0, 10000);
        assert!(
            approx(res.re, 1.0, 1e-3),
            "residue should be 1, got {}",
            res.re
        );
    }
    #[test]
    fn test_complex_gauss_identity_system() {
        let n = 3;
        let a: Vec<Vec<Complex>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            Complex::one()
                        } else {
                            Complex::zero()
                        }
                    })
                    .collect()
            })
            .collect();
        let b = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
        ];
        let x = complex_gaussian_elimination(&a, &b).unwrap();
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!(capprox(*xi, *bi, 1e-10), "I·x=b should give x=b");
        }
        let _ = a.len();
    }
    #[test]
    fn test_complex_gauss_2x2_real() {
        let a = vec![
            vec![Complex::new(2.0, 0.0), Complex::new(1.0, 0.0)],
            vec![Complex::new(1.0, 0.0), Complex::new(3.0, 0.0)],
        ];
        let b = vec![Complex::new(5.0, 0.0), Complex::new(10.0, 0.0)];
        let x = complex_gaussian_elimination(&a, &b).unwrap();
        assert!(
            approx(x[0].re, 1.0, 1e-8) && approx(x[0].im, 0.0, 1e-8),
            "x[0]={}",
            x[0].re
        );
        assert!(
            approx(x[1].re, 3.0, 1e-8) && approx(x[1].im, 0.0, 1e-8),
            "x[1]={}",
            x[1].re
        );
    }
    #[test]
    fn test_complex_gauss_singular_returns_none() {
        let a = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)],
            vec![Complex::new(2.0, 0.0), Complex::new(4.0, 0.0)],
        ];
        let b = vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)];
        assert!(
            complex_gaussian_elimination(&a, &b).is_none(),
            "singular system should return None"
        );
    }
    #[test]
    fn test_complex_gauss_complex_rhs() {
        let a = vec![vec![Complex::new(1.0, 1.0)]];
        let b = vec![Complex::new(2.0, 2.0)];
        let x = complex_gaussian_elimination(&a, &b).unwrap();
        assert!(
            approx(x[0].re, 2.0, 1e-8) && approx(x[0].im, 0.0, 1e-8),
            "x = ({}, {})",
            x[0].re,
            x[0].im
        );
    }
    #[test]
    fn test_fft_then_ifft_roundtrip() {
        let xs = vec![
            Complex::new(1.0, 0.5),
            Complex::new(-1.0, 2.0),
            Complex::new(0.5, -0.5),
            Complex::new(2.0, 1.0),
        ];
        let freq = fft_radix2(&xs);
        let recovered = ifft_radix2(&freq);
        for (orig, rec) in xs.iter().zip(recovered.iter()) {
            assert!(
                approx(orig.re, rec.re, 1e-8),
                "re: {} vs {}",
                orig.re,
                rec.re
            );
            assert!(
                approx(orig.im, rec.im, 1e-8),
                "im: {} vs {}",
                orig.im,
                rec.im
            );
        }
    }
    #[test]
    fn test_rfft_length() {
        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let out = rfft(&xs);
        assert_eq!(
            out.len(),
            xs.len() / 2 + 1,
            "rfft should return N/2+1 components"
        );
    }
    #[test]
    fn test_complex_sin_taylor_matches_exact() {
        let z = Complex::new(0.3, 0.2);
        let approx_sin = complex_sin_taylor(z, 12);
        let exact = z.sin();
        assert!(approx(approx_sin.re, exact.re, 1e-8));
        assert!(approx(approx_sin.im, exact.im, 1e-8));
    }
    #[test]
    fn test_complex_step_derivative_for_polynomial() {
        let z = Complex::new(2.0, 0.0);
        let deriv = complex_step_derivative(&|z| z * z * z, z, 1e-8);
        let expected = Complex::new(12.0, 0.0);
        assert!(
            approx(deriv.re, expected.re, 1e-4),
            "derivative re={} expected {}",
            deriv.re,
            expected.re
        );
    }
    #[test]
    fn test_cayley_transform_unit_disk() {
        let z = Complex::new(0.0, 1.0);
        let w = cayley_transform(z);
        assert!(
            approx(w.norm(), 0.0, 1e-10),
            "cayley(i) should be 0, got |w|={}",
            w.norm()
        );
    }
    #[test]
    fn test_cayley_inverse_roundtrip() {
        let z = Complex::new(1.5, 2.0);
        let w = cayley_transform(z);
        let z_back = cayley_transform_inverse(w);
        assert!(
            approx(z_back.re, z.re, 1e-8) && approx(z_back.im, z.im, 1e-8),
            "cayley inverse roundtrip failed: ({}, {})",
            z_back.re,
            z_back.im
        );
    }
    #[test]
    fn test_joukowski_unit_circle_chord() {
        let z = Complex::from_polar(1.0, PI / 3.0);
        let w = joukowski(z);
        let expected_re = 2.0 * (PI / 3.0).cos();
        assert!(approx(w.re, expected_re, 1e-10));
        assert!(
            approx(w.im, 0.0, 1e-10),
            "unit circle joukowski should be real, got im={}",
            w.im
        );
    }
    #[test]
    fn test_schwarz_christoffel_at_origin_is_zero() {
        let result = schwarz_christoffel_rectangle(Complex::zero(), 0.5, 100);
        assert!(approx(result.norm(), 0.0, 1e-10));
    }
    #[test]
    fn test_durand_kerner_quadratic() {
        let coeffs = vec![Complex::new(-1.0, 0.0), Complex::new(0.0, 0.0)];
        let roots = durand_kerner_roots(&coeffs, 200);
        assert_eq!(roots.len(), 2);
        let norms: Vec<f64> = roots.iter().map(|r| (r.norm() - 1.0).abs()).collect();
        for n in &norms {
            assert!(*n < 0.2, "root magnitude should be ~1, error={n}");
        }
    }
    #[test]
    fn test_dft_matrix_row0_all_ones() {
        let m = dft_matrix(4);
        for entry in &m[0] {
            assert!(approx(entry.re, 1.0, 1e-10) && approx(entry.im, 0.0, 1e-10));
        }
    }
}
/// Compute the Short-Time Fourier Transform of a real signal.
///
/// # Arguments
/// * `signal` – real-valued input signal.
/// * `window_size` – number of samples per window (power of two recommended).
/// * `hop_size` – step between successive windows (≤ window_size).
/// * `window_fn` – window function applied to each frame (e.g., Hann window).
///
/// Returns a `Vec<Vec`Complex`>` where each inner Vec is the FFT of one frame.
/// The outer dimension is time (frames), inner is frequency.
#[allow(dead_code)]
pub fn stft(
    signal: &[f64],
    window_size: usize,
    hop_size: usize,
    window_fn: &[f64],
) -> Vec<Vec<Complex>> {
    if signal.is_empty() || window_size == 0 || hop_size == 0 {
        return Vec::new();
    }
    let hop = hop_size.max(1);
    let wlen = window_fn.len().min(window_size);
    let mut frames = Vec::new();
    let mut start = 0;
    while start < signal.len() {
        let mut frame: Vec<Complex> = (0..window_size)
            .map(|i| {
                let sig_i = start + i;
                let sig_val = if sig_i < signal.len() {
                    signal[sig_i]
                } else {
                    0.0
                };
                let win_val = if i < wlen { window_fn[i] } else { 1.0 };
                Complex::new(sig_val * win_val, 0.0)
            })
            .collect();
        let n = frame.len().next_power_of_two();
        frame.resize(n, Complex::zero());
        fft_inplace(&mut frame, false);
        frames.push(frame);
        start += hop;
    }
    frames
}
/// Generate a Hann window of length `n`.
///
/// `w[k] = sin²(π k / (n-1))` for k = 0..n-1.
#[allow(dead_code)]
pub fn hann_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|k| {
            let x = std::f64::consts::PI * k as f64 / (n - 1) as f64;
            x.sin().powi(2)
        })
        .collect()
}
/// Generate a Hamming window of length `n`.
///
/// `w[k] = 0.54 - 0.46 cos(2π k / (n-1))`
#[allow(dead_code)]
pub fn hamming_window(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }
    (0..n)
        .map(|k| 0.54 - 0.46 * (2.0 * std::f64::consts::PI * k as f64 / (n - 1) as f64).cos())
        .collect()
}
/// Generate a rectangular (boxcar) window of length `n`.
#[allow(dead_code)]
pub fn rectangular_window(n: usize) -> Vec<f64> {
    vec![1.0; n]
}
/// Compute the inverse STFT (overlap-add reconstruction).
///
/// # Arguments
/// * `frames` – STFT frames from [`stft`] (each of length window_size).
/// * `hop_size` – hop size used when computing the STFT.
/// * `window_fn` – the same window used for analysis (for proper reconstruction).
///
/// Returns the reconstructed real signal.
#[allow(dead_code)]
pub fn istft(frames: &[Vec<Complex>], hop_size: usize, window_fn: &[f64]) -> Vec<f64> {
    if frames.is_empty() {
        return Vec::new();
    }
    let window_size = frames[0].len();
    let hop = hop_size.max(1);
    let n_out = hop * frames.len() + window_size;
    let mut signal = vec![0.0f64; n_out];
    let mut weight = vec![0.0f64; n_out];
    let wlen = window_fn.len().min(window_size);
    for (fi, frame) in frames.iter().enumerate() {
        let mut f = frame.clone();
        fft_inplace(&mut f, true);
        let inv_n = 1.0 / f.len() as f64;
        let start = fi * hop;
        for k in 0..window_size.min(f.len()) {
            let win = if k < wlen { window_fn[k] } else { 1.0 };
            if start + k < signal.len() {
                signal[start + k] += f[k].re * inv_n * win;
                weight[start + k] += win * win;
            }
        }
    }
    for i in 0..signal.len() {
        if weight[i] > 1e-12 {
            signal[i] /= weight[i];
        }
    }
    signal
}
/// Compute the DCT-II of a real signal (the "standard" DCT used in JPEG/MP3).
///
/// `X[k] = 2 Σ_{n=0}^{N-1} x[n] cos(π(2n+1)k / (2N))`
///
/// Uses a naive O(N²) implementation.
#[allow(dead_code)]
pub fn dct_ii(xs: &[f64]) -> Vec<f64> {
    let n = xs.len();
    if n == 0 {
        return Vec::new();
    }
    let pi = std::f64::consts::PI;
    (0..n)
        .map(|k| {
            let sum: f64 = xs
                .iter()
                .enumerate()
                .map(|(m, &x)| x * (pi * (2 * m + 1) as f64 * k as f64 / (2.0 * n as f64)).cos())
                .sum();
            2.0 * sum
        })
        .collect()
}
/// Compute the inverse DCT-II (DCT-III).
///
/// `x[n] = (1/N) (X[0]/2 + Σ_{k=1}^{N-1} X[k] cos(π(2n+1)k / (2N)))`
#[allow(dead_code)]
pub fn idct_ii(xs: &[f64]) -> Vec<f64> {
    let n = xs.len();
    if n == 0 {
        return Vec::new();
    }
    let pi = std::f64::consts::PI;
    (0..n)
        .map(|m| {
            let dc = xs[0] / 2.0;
            let sum: f64 = (1..n)
                .map(|k| xs[k] * (pi * (2 * m + 1) as f64 * k as f64 / (2.0 * n as f64)).cos())
                .sum();
            (dc + sum) / n as f64
        })
        .collect()
}
/// Compute the DCT-I of a real signal.
///
/// `X[k] = x[0] + (-1)^k x[N-1] + 2 Σ_{n=1}^{N-2} x[n] cos(π n k / (N-1))`
///
/// The DCT-I is its own inverse (when normalized).
#[allow(dead_code)]
pub fn dct_i(xs: &[f64]) -> Vec<f64> {
    let n = xs.len();
    if n < 2 {
        return xs.to_vec();
    }
    let pi = std::f64::consts::PI;
    let nm1 = (n - 1) as f64;
    (0..n)
        .map(|k| {
            let mut sum = xs[0] + if k % 2 == 0 { xs[n - 1] } else { -xs[n - 1] };
            for m in 1..(n - 1) {
                sum += 2.0 * xs[m] * (pi * m as f64 * k as f64 / nm1).cos();
            }
            sum
        })
        .collect()
}
/// Schur-like decomposition of a 2×2 complex matrix.
///
/// Returns `(T, U)` where:
/// - `T` is upper-triangular (quasi-Schur form)
/// - `U` is unitary such that `A = U T U†`
///
/// For 2×2 matrices the decomposition is computed directly from eigenvalues.
/// Returns `None` if the matrix is degenerate.
#[allow(dead_code)]
pub fn schur_2x2(m: &ComplexMat2) -> Option<(ComplexMat2, ComplexMat2)> {
    let (l1, l2) = m.eigenvalues();
    let al1 = ComplexMat2::new(m.a - l1, m.b, m.c, m.d - l1);
    let (u0, u1) = if m.b.norm() > 1e-12 {
        (m.b, l1 - m.a)
    } else if m.c.norm() > 1e-12 {
        (l1 - m.d, m.c)
    } else {
        let u = ComplexMat2::identity();
        let t = ComplexMat2::new(l1, m.b, Complex::zero(), l2);
        return Some((t, u));
    };
    let _ = al1;
    let norm = (u0.norm_sq() + u1.norm_sq()).sqrt();
    if norm < 1e-12 {
        return None;
    }
    let v0 = Complex::new(u0.re / norm, u0.im / norm);
    let v1 = Complex::new(u1.re / norm, u1.im / norm);
    let u = ComplexMat2::new(v0, -v1.conj(), v1, v0.conj());
    let uh = u.conjugate_transpose();
    let t = uh.mul(m).mul(&u);
    Some((t, u))
}
/// Compute the power spectral density (PSD) of a real signal.
///
/// Returns `|X[k]|²` for k = 0..N/2+1.
#[allow(dead_code)]
pub fn power_spectral_density(xs: &[f64]) -> Vec<f64> {
    let freq = rfft(xs);
    freq.iter().map(|z| z.norm_sq()).collect()
}
/// Compute the magnitude spectrum of a real signal.
///
/// Returns `|X[k]|` for k = 0..N/2+1.
#[allow(dead_code)]
pub fn magnitude_spectrum(xs: &[f64]) -> Vec<f64> {
    let freq = rfft(xs);
    freq.iter().map(|z| z.norm()).collect()
}
/// Compute the phase spectrum of a real signal.
///
/// Returns `arg(X[k])` for k = 0..N/2+1.
#[allow(dead_code)]
pub fn phase_spectrum(xs: &[f64]) -> Vec<f64> {
    let freq = rfft(xs);
    freq.iter().map(|z| z.arg()).collect()
}
/// Convolve two real signals using FFT.
///
/// Returns the linear convolution `a * b` of length `a.len() + b.len() - 1`.
#[allow(dead_code)]
pub fn fft_convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let n = out_len.next_power_of_two();
    let mut fa: Vec<Complex> = a.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fa.resize(n, Complex::zero());
    let mut fb: Vec<Complex> = b.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fb.resize(n, Complex::zero());
    fft_inplace(&mut fa, false);
    fft_inplace(&mut fb, false);
    let mut fc: Vec<Complex> = fa.iter().zip(fb.iter()).map(|(&ai, &bi)| ai * bi).collect();
    fft_inplace(&mut fc, true);
    let inv_n = 1.0 / n as f64;
    fc[..out_len].iter().map(|z| z.re * inv_n).collect()
}
/// Compute Parseval's theorem check: sum of |x\[n\]|² == (1/N) sum |X\[k\]|².
///
/// Returns `(time_energy, freq_energy)`.  These should be approximately equal.
#[allow(dead_code)]
pub fn parseval_check(xs: &[f64]) -> (f64, f64) {
    let time_energy: f64 = xs.iter().map(|&x| x * x).sum();
    let complex_in: Vec<Complex> = xs.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let freq = fft_radix2(&complex_in);
    let n = freq.len() as f64;
    let freq_energy: f64 = freq.iter().map(|z| z.norm_sq()).sum::<f64>() / n;
    (time_energy, freq_energy)
}
/// Evaluate a polynomial with complex coefficients at complex point `z`.
///
/// Coefficients are in descending order: `coeffs[0] z^n + ... + coeffs[n]`.
#[allow(dead_code)]
pub fn poly_eval(coeffs: &[Complex], z: Complex) -> Complex {
    if coeffs.is_empty() {
        return Complex::zero();
    }
    let mut result = coeffs[0];
    for &c in &coeffs[1..] {
        result = result * z + c;
    }
    result
}
/// Compute the derivative of a polynomial with complex coefficients.
///
/// Returns the derivative coefficients (one fewer term).
#[allow(dead_code)]
pub fn poly_derivative(coeffs: &[Complex]) -> Vec<Complex> {
    if coeffs.len() <= 1 {
        return Vec::new();
    }
    let n = coeffs.len() - 1;
    (0..n)
        .map(|i| coeffs[i] * Complex::new((n - i) as f64, 0.0))
        .collect()
}
/// Multiply two polynomials with complex coefficients.
///
/// Coefficients in descending order.
#[allow(dead_code)]
pub fn poly_multiply(a: &[Complex], b: &[Complex]) -> Vec<Complex> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let na = a.len();
    let nb = b.len();
    let mut result = vec![Complex::zero(); na + nb - 1];
    for i in 0..na {
        for j in 0..nb {
            result[i + j] = result[i + j] + a[i] * b[j];
        }
    }
    result
}
/// Evaluate the Z-transform of a finite-length sequence at complex frequency `z`.
///
/// `X(z) = Σ_{n=0}^{N-1} x[n] z^{-n}`
#[allow(dead_code)]
pub fn z_transform_eval(xs: &[f64], z: Complex) -> Complex {
    let mut result = Complex::zero();
    let mut zn = Complex::one();
    let inv_z = if z.norm_sq() > 1e-24 {
        Complex::one() / z
    } else {
        Complex::zero()
    };
    for &x in xs {
        result = result + Complex::new(x, 0.0) * zn;
        zn = zn * inv_z;
    }
    result
}
/// Evaluate the Laplace transform of an exponentially decaying signal
/// `f(t) = x * exp(-alpha * t)` at complex frequency `s`.
///
/// `F(s) = x / (s + alpha)` for `Re(s + alpha) > 0`.
#[allow(dead_code)]
pub fn laplace_exp_decay(amplitude: f64, alpha: f64, s: Complex) -> Complex {
    let denom = s + Complex::new(alpha, 0.0);
    if denom.norm_sq() < 1e-24 {
        return Complex::zero();
    }
    Complex::new(amplitude, 0.0) / denom
}
#[cfg(test)]
mod tests_new_complex2 {

    use crate::complex::Complex;
    use crate::complex::ComplexMat2;
    use crate::complex::ComplexMatN;

    use crate::complex::dct_i;
    use crate::complex::dct_ii;

    use crate::complex::fft_convolve;

    use crate::complex::functions::PI;
    use crate::complex::hamming_window;
    use crate::complex::hann_window;
    use crate::complex::idct_ii;

    use crate::complex::laplace_exp_decay;
    use crate::complex::magnitude_spectrum;

    use crate::complex::parseval_check;
    use crate::complex::phase_spectrum;
    use crate::complex::poly_derivative;
    use crate::complex::poly_eval;
    use crate::complex::poly_multiply;
    use crate::complex::power_spectral_density;

    use crate::complex::rectangular_window;

    use crate::complex::schur_2x2;

    use crate::complex::stft;
    use crate::complex::z_transform_eval;
    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }
    #[test]
    fn test_stft_produces_frames() {
        let signal: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();
        let window = hann_window(16);
        let frames = stft(&signal, 16, 8, &window);
        assert!(!frames.is_empty(), "STFT should produce at least one frame");
    }
    #[test]
    fn test_stft_frame_count() {
        let signal: Vec<f64> = vec![0.0; 64];
        let window = rectangular_window(16);
        let frames = stft(&signal, 16, 16, &window);
        assert!(
            frames.len() >= 4,
            "should have at least 4 frames, got {}",
            frames.len()
        );
    }
    #[test]
    fn test_stft_dc_signal() {
        let signal: Vec<f64> = vec![1.0; 16];
        let window = rectangular_window(16);
        let frames = stft(&signal, 16, 16, &window);
        assert!(!frames.is_empty());
        let frame = &frames[0];
        let dc = frame[0].norm();
        let max_ac = frame[1..frame.len() / 2]
            .iter()
            .map(|z| z.norm())
            .fold(0.0f64, f64::max);
        assert!(dc > max_ac, "DC should dominate: dc={dc}, max_ac={max_ac}");
    }
    #[test]
    fn test_hann_window_zero_endpoints() {
        let w = hann_window(8);
        assert!((w[0]).abs() < 1e-10, "Hann window should start at 0");
        assert!(!w.is_empty());
    }
    #[test]
    fn test_hann_window_length() {
        let w = hann_window(16);
        assert_eq!(w.len(), 16);
    }
    #[test]
    fn test_hamming_window_nonzero_endpoints() {
        let w = hamming_window(8);
        assert!(w[0] > 0.0 && w[0] < 0.2, "Hamming endpoint = {}", w[0]);
    }
    #[test]
    fn test_rectangular_window_all_ones() {
        let w = rectangular_window(8);
        for &wi in &w {
            assert!((wi - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_dct_ii_dc_component() {
        let n = 8;
        let xs = vec![1.0f64; n];
        let out = dct_ii(&xs);
        assert!(
            approx(out[0], 2.0 * n as f64, 1e-8),
            "DCT-II[0] of constant 1 should be 2N={}, got {}",
            2 * n,
            out[0]
        );
    }
    #[test]
    fn test_dct_ii_idct_roundtrip() {
        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let dct = dct_ii(&xs);
        let recovered = idct_ii(&dct);
        assert_eq!(recovered.len(), xs.len());
        for (a, b) in xs.iter().zip(recovered.iter()) {
            assert!(approx(*a, *b, 1e-8), "DCT-II roundtrip: {} vs {}", a, b);
        }
    }
    #[test]
    fn test_dct_ii_orthogonality() {
        let n = 8;
        let xs1: Vec<f64> = (0..n)
            .map(|m| (PI * (2 * m + 1) as f64 * 1.0 / (2.0 * n as f64)).cos())
            .collect();
        let xs2: Vec<f64> = (0..n)
            .map(|m| (PI * (2 * m + 1) as f64 * 2.0 / (2.0 * n as f64)).cos())
            .collect();
        let dot: f64 = xs1.iter().zip(xs2.iter()).map(|(a, b)| a * b).sum();
        assert!(
            dot.abs() < 1e-8,
            "DCT basis vectors should be orthogonal, dot={dot}"
        );
    }
    #[test]
    fn test_dct_ii_single_element() {
        let xs = vec![3.0];
        let out = dct_ii(&xs);
        assert_eq!(out.len(), 1);
        assert!(
            approx(out[0], 6.0, 1e-10),
            "DCT-II of [3] should be [6], got {}",
            out[0]
        );
    }
    #[test]
    fn test_dct_i_self_inverse() {
        let xs = vec![1.0, 2.0, 3.0, 2.0, 1.0];
        let dct1 = dct_i(&xs);
        let dct2 = dct_i(&dct1);
        let scale = 2.0 * (xs.len() - 1) as f64;
        for (a, b) in xs.iter().zip(dct2.iter()) {
            assert!(
                approx(*a * scale, *b, 1e-6),
                "DCT-I self-inverse: {} vs {}",
                a * scale,
                b
            );
        }
    }
    #[test]
    fn test_schur_2x2_upper_triangular() {
        let m = ComplexMat2::new(
            Complex::new(3.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
        );
        if let Some((t, _u)) = schur_2x2(&m) {
            assert!(
                t.c.norm() < 0.1,
                "Schur form should have near-zero lower-left, got {}",
                t.c.norm()
            );
        }
    }
    #[test]
    fn test_schur_2x2_diagonal_matrix() {
        let m = ComplexMat2::new(
            Complex::new(2.0, 0.0),
            Complex::zero(),
            Complex::zero(),
            Complex::new(3.0, 0.0),
        );
        let result = schur_2x2(&m);
        assert!(result.is_some(), "diagonal matrix Schur should succeed");
    }
    #[test]
    fn test_schur_2x2_eigenvalue_trace() {
        let m = ComplexMat2::new(
            Complex::new(4.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
        );
        if let Some((t, _u)) = schur_2x2(&m) {
            let tr_m = m.trace();
            let tr_t = t.trace();
            assert!(
                approx(tr_m.re, tr_t.re, 1e-6) && approx(tr_m.im, tr_t.im, 1e-6),
                "traces should match: ({},{}) vs ({},{})",
                tr_m.re,
                tr_m.im,
                tr_t.re,
                tr_t.im
            );
        }
    }
    #[test]
    fn test_complex_matn_identity_mul() {
        let id = ComplexMatN::identity(3);
        let m = ComplexMatN::identity(3);
        let result = id.matmul(&m);
        assert_eq!(result.n, 3);
        for i in 0..3 {
            assert!(approx(result.get(i, i).re, 1.0, 1e-12));
        }
    }
    #[test]
    fn test_complex_matn_trace() {
        let mut m = ComplexMatN::zeros(3);
        for i in 0..3 {
            m[(i, i)] = Complex::new(i as f64 + 1.0, 0.0);
        }
        let tr = m.trace();
        assert!(
            approx(tr.re, 6.0, 1e-12),
            "trace should be 6, got {}",
            tr.re
        );
    }
    #[test]
    fn test_complex_matn_frobenius_identity() {
        let id = ComplexMatN::identity(4);
        let frob = id.frobenius_norm();
        assert!(
            approx(frob, 2.0, 1e-10),
            "Frobenius norm of 4x4 identity should be 2, got {frob}"
        );
    }
    #[test]
    fn test_complex_matn_matvec() {
        let id = ComplexMatN::identity(3);
        let v = vec![
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
        ];
        let result = id.matvec(&v);
        for (a, b) in v.iter().zip(result.iter()) {
            assert!(approx(a.re, b.re, 1e-12));
        }
    }
    #[test]
    fn test_complex_matn_conjugate_transpose() {
        let mut m = ComplexMatN::zeros(2);
        m[(0, 1)] = Complex::new(1.0, 2.0);
        let mh = m.conjugate_transpose();
        assert!(approx(mh.get(1, 0).re, 1.0, 1e-12));
        assert!(approx(mh.get(1, 0).im, -2.0, 1e-12));
    }
    #[test]
    fn test_parseval_theorem() {
        let xs: Vec<f64> = (0..8).map(|i| (i as f64).sin()).collect();
        let (te, fe) = parseval_check(&xs);
        assert!(
            approx(te, fe, 1e-8),
            "Parseval: time_energy={te}, freq_energy={fe}"
        );
    }
    #[test]
    fn test_power_spectral_density_dc() {
        let xs = vec![1.0f64; 8];
        let psd = power_spectral_density(&xs);
        assert!(!psd.is_empty());
        let max_val = *psd.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        assert!((max_val - psd[0]).abs() < 1e-8, "DC should dominate PSD");
    }
    #[test]
    fn test_magnitude_spectrum_length() {
        let xs = vec![0.0f64; 16];
        let mag = magnitude_spectrum(&xs);
        assert_eq!(mag.len(), 9, "rfft of 16 samples should have 16/2+1=9 bins");
    }
    #[test]
    fn test_fft_convolve_unit_impulse() {
        let xs = vec![1.0, 2.0, 3.0];
        let impulse = vec![1.0];
        let result = fft_convolve(&xs, &impulse);
        assert_eq!(result.len(), xs.len());
        for (a, b) in xs.iter().zip(result.iter()) {
            assert!(
                approx(*a, *b, 1e-8),
                "convolve with impulse: {} vs {}",
                a,
                b
            );
        }
    }
    #[test]
    fn test_fft_convolve_commutativity() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![0.5, 1.0];
        let ab = fft_convolve(&a, &b);
        let ba = fft_convolve(&b, &a);
        assert_eq!(ab.len(), ba.len());
        for (x, y) in ab.iter().zip(ba.iter()) {
            assert!(
                approx(*x, *y, 1e-8),
                "convolution should be commutative: {} vs {}",
                x,
                y
            );
        }
    }
    #[test]
    fn test_poly_eval_constant() {
        let coeffs = vec![Complex::new(5.0, 0.0)];
        let v = poly_eval(&coeffs, Complex::new(3.0, 0.0));
        assert!(approx(v.re, 5.0, 1e-12));
    }
    #[test]
    fn test_poly_eval_linear() {
        let coeffs = vec![Complex::new(2.0, 0.0), Complex::new(3.0, 0.0)];
        let v = poly_eval(&coeffs, Complex::new(1.0, 0.0));
        assert!(approx(v.re, 5.0, 1e-12));
    }
    #[test]
    fn test_poly_eval_quadratic() {
        let coeffs = vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(-1.0, 0.0),
        ];
        let v1 = poly_eval(&coeffs, Complex::new(1.0, 0.0));
        let v2 = poly_eval(&coeffs, Complex::new(-1.0, 0.0));
        assert!(approx(v1.re, 0.0, 1e-12));
        assert!(approx(v2.re, 0.0, 1e-12));
    }
    #[test]
    fn test_poly_derivative_constant() {
        let coeffs = vec![Complex::new(5.0, 0.0)];
        let deriv = poly_derivative(&coeffs);
        assert!(deriv.is_empty(), "derivative of constant should be empty");
    }
    #[test]
    fn test_poly_derivative_linear() {
        let coeffs = vec![Complex::new(3.0, 0.0), Complex::new(2.0, 0.0)];
        let deriv = poly_derivative(&coeffs);
        assert_eq!(deriv.len(), 1);
        assert!(approx(deriv[0].re, 3.0, 1e-12));
    }
    #[test]
    fn test_poly_multiply_monomial() {
        let a = vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)];
        let b = vec![Complex::new(1.0, 0.0), Complex::new(-1.0, 0.0)];
        let product = poly_multiply(&a, &b);
        assert_eq!(product.len(), 3);
        assert!(approx(product[0].re, 1.0, 1e-12));
        assert!(approx(product[1].re, 0.0, 1e-12));
        assert!(approx(product[2].re, -1.0, 1e-12));
    }
    #[test]
    fn test_z_transform_impulse() {
        let xs = vec![1.0, 0.0, 0.0, 0.0];
        let result = z_transform_eval(&xs, Complex::new(2.0, 0.0));
        assert!(approx(result.re, 1.0, 1e-12), "Z{{delta}} should be 1");
    }
    #[test]
    fn test_z_transform_delay() {
        let xs = vec![0.0, 1.0, 0.0];
        let result = z_transform_eval(&xs, Complex::new(2.0, 0.0));
        assert!(
            approx(result.re, 0.5, 1e-12),
            "Z{{delta[n-1]}} at z=2 should be 0.5, got {}",
            result.re
        );
    }
    #[test]
    fn test_laplace_exp_decay_at_zero() {
        let alpha = 2.0;
        let result = laplace_exp_decay(1.0, alpha, Complex::zero());
        let expected = 1.0 / alpha;
        assert!(
            approx(result.re, expected, 1e-12),
            "Laplace F(0) should be 1/alpha={expected}"
        );
    }
    #[test]
    fn test_laplace_exp_decay_pole_returns_zero() {
        let alpha = 1.0;
        let result = laplace_exp_decay(1.0, alpha, Complex::new(-alpha, 0.0));
        assert!(result.norm() == 0.0, "at pole should return zero");
    }
    #[test]
    fn test_phase_spectrum_length() {
        let xs: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let phase = phase_spectrum(&xs);
        assert_eq!(phase.len(), 5, "phase of 8 samples should have 5 bins");
    }
}
