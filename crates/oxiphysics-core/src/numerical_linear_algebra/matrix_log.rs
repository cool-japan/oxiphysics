//! Matrix logarithm by inverse scaling and squaring with a diagonal Padé
//! approximant evaluated in partial-fraction form.
//!
//! For `‖X‖ < 1` the logarithm is `log(I + X) = ∫₀¹ X (I + tX)⁻¹ dt`.  Replacing
//! the integral by the `m`-point Gauss-Legendre rule with nodes `tⱼ` and weights
//! `wⱼ` on `[0, 1]` gives exactly the `[m/m]` Padé approximant of `log(1 + x)`,
//!
//! ```text
//! r_m(X) = Σⱼ wⱼ X (I + tⱼ X)⁻¹,
//! ```
//!
//! which for `m = 7` is accurate to about the unit round-off for `‖X‖₁ ≤ 0.25`
//! (Al-Mohy & Higham, 2012).  Larger arguments are first reduced with
//! `A ← A^(1/2)` (`k` times) and the result is multiplied by `2^k`.

use std::f64::consts::PI;

use super::types::MatrixFunctions;

/// Order `m` of the Padé approximant (number of Gauss-Legendre nodes).
const PADE_ORDER: usize = 7;
/// Largest `‖A^(1/2^k) − I‖₁` for which the order-7 approximant is accurate to
/// double precision.
const PADE_NORM_LIMIT: f64 = 0.25;
/// Safety cap on the number of square roots.
const MAX_SQUARE_ROOTS: u32 = 64;

/// Matrix logarithm of the row-major `n × n` matrix `a`.
pub(super) fn logm_dense(a: &[f64], n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    let identity = MatrixFunctions::eye(n);
    let mut current = a.to_vec();
    let mut square_roots = 0u32;
    while square_roots < MAX_SQUARE_ROOTS {
        let norm = one_norm(&subtract(&current, &identity), n);
        if norm <= PADE_NORM_LIMIT || !norm.is_finite() {
            break;
        }
        current = MatrixFunctions::sqrtm(&current, n);
        square_roots += 1;
    }
    let x = subtract(&current, &identity);
    let scale = 2.0_f64.powi(square_roots as i32);
    pade_log_one_plus(&x, n)
        .into_iter()
        .map(|v| v * scale)
        .collect()
}

/// `log(I + X)` for a matrix `X` of small norm via the partial-fraction form of
/// the diagonal Padé approximant.
fn pade_log_one_plus(x: &[f64], n: usize) -> Vec<f64> {
    let mut result = vec![0.0_f64; n * n];
    for (node, weight) in gauss_legendre_unit_interval(PADE_ORDER) {
        let mut lhs = MatrixFunctions::eye(n);
        for (entry, xv) in lhs.iter_mut().zip(x.iter()) {
            *entry += node * xv;
        }
        // (I + t X) Z = X, i.e. Z = X (I + t X)⁻¹ because the two commute.
        let z = solve_dense(&lhs, x, n);
        for (r, zv) in result.iter_mut().zip(z.iter()) {
            *r += weight * zv;
        }
    }
    result
}

/// `a - b` for equally long slices.
fn subtract(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

/// Matrix 1-norm (largest absolute column sum) of a row-major `n × n` matrix.
fn one_norm(a: &[f64], n: usize) -> f64 {
    (0..n)
        .map(|j| (0..n).map(|i| a[i * n + j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

/// Solves `A Z = B` for row-major `n × n` matrices by Gaussian elimination with
/// partial pivoting (no explicit inverse).
fn solve_dense(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut lu = a.to_vec();
    let mut z = b.to_vec();
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_abs = lu[col * n + col].abs();
        for row in (col + 1)..n {
            let candidate = lu[row * n + col].abs();
            if candidate > pivot_abs {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if pivot_row != col {
            for j in 0..n {
                lu.swap(col * n + j, pivot_row * n + j);
                z.swap(col * n + j, pivot_row * n + j);
            }
        }
        let pivot = if lu[col * n + col] == 0.0 {
            f64::MIN_POSITIVE
        } else {
            lu[col * n + col]
        };
        for row in (col + 1)..n {
            let factor = lu[row * n + col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for j in col..n {
                let upper = lu[col * n + j];
                lu[row * n + j] -= factor * upper;
            }
            for j in 0..n {
                let upper = z[col * n + j];
                z[row * n + j] -= factor * upper;
            }
        }
    }
    for j in 0..n {
        for i in (0..n).rev() {
            let mut sum = z[i * n + j];
            for k in (i + 1)..n {
                sum -= lu[i * n + k] * z[k * n + j];
            }
            let diagonal = if lu[i * n + i] == 0.0 {
                f64::MIN_POSITIVE
            } else {
                lu[i * n + i]
            };
            z[i * n + j] = sum / diagonal;
        }
    }
    z
}

/// `order`-point Gauss-Legendre rule on `[0, 1]` as `(node, weight)` pairs, with
/// the nodes found by Newton's method on the Legendre polynomial.
fn gauss_legendre_unit_interval(order: usize) -> Vec<(f64, f64)> {
    let m = order as f64;
    let mut rule = Vec::with_capacity(order);
    for i in 0..order {
        let mut z = (PI * (i as f64 + 0.75) / (m + 0.5)).cos();
        let mut derivative = 1.0_f64;
        for _ in 0..100 {
            let mut previous = 1.0_f64;
            let mut current = z;
            for degree in 2..=order {
                let d = degree as f64;
                let next = ((2.0 * d - 1.0) * z * current - (d - 1.0) * previous) / d;
                previous = current;
                current = next;
            }
            derivative = m * (z * current - previous) / (z * z - 1.0);
            let step = current / derivative;
            z -= step;
            if step.abs() < 1.0e-16 {
                break;
            }
        }
        let weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
        rule.push((0.5 * (z + 1.0), 0.5 * weight));
    }
    rule
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_legendre_rule_integrates_polynomials_exactly() {
        for order in [1_usize, 2, 3, 5, 7, 10] {
            let rule = gauss_legendre_unit_interval(order);
            assert_eq!(rule.len(), order);
            let weight_sum: f64 = rule.iter().map(|(_, w)| w).sum();
            assert!((weight_sum - 1.0).abs() < 1e-14, "order {order}");
            // Exact for every monomial of degree <= 2*order - 1: integral t^p = 1/(p+1).
            for power in 0..(2 * order) {
                let quadrature: f64 = rule.iter().map(|(t, w)| w * t.powi(power as i32)).sum();
                let exact = 1.0 / (power as f64 + 1.0);
                assert!(
                    (quadrature - exact).abs() < 1e-14,
                    "order {order}, t^{power}: {quadrature} vs {exact}"
                );
            }
            assert!(rule.iter().all(|(t, w)| *t > 0.0 && *t < 1.0 && *w > 0.0));
        }
    }

    #[test]
    fn seven_point_pade_matches_the_scalar_logarithm() {
        // Scalar case n = 1: r_7(x) must agree with ln(1 + x) to rounding error on
        // the whole admissible range |x| <= 0.25, including negative arguments.
        for i in -25_i32..=25 {
            let x = f64::from(i) / 100.0;
            let approx = pade_log_one_plus(&[x], 1)[0];
            let exact = x.ln_1p();
            assert!(
                (approx - exact).abs() <= 5e-16 * exact.abs() + 2e-17,
                "x = {x}: {approx} vs {exact}"
            );
        }
    }

    #[test]
    fn solve_dense_matches_a_hand_computed_system() {
        // [[2, 1], [1, 3]] Z = [[1, 0], [0, 1]]  =>  Z = 1/5 [[3, -1], [-1, 2]].
        let z = solve_dense(&[2.0, 1.0, 1.0, 3.0], &[1.0, 0.0, 0.0, 1.0], 2);
        let expected = [0.6, -0.2, -0.2, 0.4];
        for (got, want) in z.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-14, "{got} vs {want}");
        }
        // A system that needs a row swap (zero leading pivot).
        let z = solve_dense(&[0.0, 1.0, 1.0, 0.0], &[2.0, 0.0, 0.0, 3.0], 2);
        let expected = [0.0, 3.0, 2.0, 0.0];
        for (got, want) in z.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-14, "{got} vs {want}");
        }
    }
}
