//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Interval, IntervalNewtonResult};

/// Interval sqrt: `sqrt([lo, hi]) = [sqrt(lo), sqrt(hi)]`.
///
/// Panics if `lo < 0`.
pub fn sqrt_interval(a: Interval) -> Interval {
    assert!(a.lo >= 0.0, "sqrt_interval: lo must be ≥ 0, got {}", a.lo);
    Interval {
        lo: a.lo.sqrt(),
        hi: a.hi.sqrt(),
    }
}
/// Interval absolute value.
///
/// - If `lo >= 0`: `[lo, hi]`
/// - If `hi <= 0`: `[-hi, -lo]`
/// - Otherwise (0 ∈ interval): `[0, max(|lo|, |hi|)]`
pub fn abs_interval(a: Interval) -> Interval {
    if a.lo >= 0.0 {
        a
    } else if a.hi <= 0.0 {
        Interval {
            lo: -a.hi,
            hi: -a.lo,
        }
    } else {
        Interval {
            lo: 0.0,
            hi: (-a.lo).max(a.hi),
        }
    }
}
/// Tighter squaring of an interval (tighter than `a * a`).
///
/// - If `0 ∈ a`: `[0, max(lo², hi²)]`
/// - Otherwise: `[min(lo², hi²), max(lo², hi²)]`
pub fn square(a: Interval) -> Interval {
    let lo_sq = a.lo * a.lo;
    let hi_sq = a.hi * a.hi;
    if a.lo <= 0.0 && a.hi >= 0.0 {
        Interval {
            lo: 0.0,
            hi: lo_sq.max(hi_sq),
        }
    } else {
        Interval {
            lo: lo_sq.min(hi_sq),
            hi: lo_sq.max(hi_sq),
        }
    }
}
/// Interval power: compute `a^n` for non-negative integer `n`.
pub fn pow_interval(a: Interval, n: u32) -> Interval {
    match n {
        0 => Interval::point(1.0),
        1 => a,
        2 => square(a),
        _ => {
            if n.is_multiple_of(2) {
                let abs_a = abs_interval(a);
                Interval {
                    lo: abs_a.lo.powi(n as i32),
                    hi: abs_a.hi.powi(n as i32),
                }
            } else {
                Interval {
                    lo: a.lo.powi(n as i32),
                    hi: a.hi.powi(n as i32),
                }
            }
        }
    }
}
/// Interval sin: conservative enclosure of sin over `[lo, hi]`.
pub fn sin_interval(a: Interval) -> Interval {
    let sin_lo = a.lo.sin();
    let sin_hi = a.hi.sin();
    let mut lo = sin_lo.min(sin_hi);
    let mut hi = sin_lo.max(sin_hi);
    let two_pi = 2.0 * std::f64::consts::PI;
    let half_pi = std::f64::consts::FRAC_PI_2;
    let k_start = ((a.lo - half_pi) / two_pi).ceil() as i64;
    let k_end = ((a.hi - half_pi) / two_pi).floor() as i64;
    if k_start <= k_end {
        hi = 1.0;
    }
    let k_start2 = ((a.lo + half_pi) / two_pi).ceil() as i64;
    let k_end2 = ((a.hi + half_pi) / two_pi).floor() as i64;
    if k_start2 <= k_end2 {
        lo = -1.0;
    }
    Interval { lo, hi }
}
/// Interval cos: conservative enclosure of cos over `[lo, hi]`.
pub fn cos_interval(a: Interval) -> Interval {
    let shifted = Interval::new(
        a.lo + std::f64::consts::FRAC_PI_2,
        a.hi + std::f64::consts::FRAC_PI_2,
    );
    sin_interval(shifted)
}
/// Interval exp: `exp([lo, hi]) = [exp(lo), exp(hi)]`.
pub fn exp_interval(a: Interval) -> Interval {
    Interval {
        lo: a.lo.exp(),
        hi: a.hi.exp(),
    }
}
/// Interval ln: `ln([lo, hi]) = [ln(lo), ln(hi)]`. Panics if `lo <= 0`.
pub fn ln_interval(a: Interval) -> Interval {
    assert!(a.lo > 0.0, "ln_interval: lo must be > 0, got {}", a.lo);
    Interval {
        lo: a.lo.ln(),
        hi: a.hi.ln(),
    }
}
/// Perform interval Newton iteration on a function `f` with derivative `f'`.
///
/// Given `x_k`, the next iterate is:
///   N(x_k) = mid(x_k) - f(mid(x_k)) / F'(x_k)
/// intersected with x_k.
///
/// Uses an explicit work-stack to avoid recursive monomorphization issues.
///
/// - `f`: function evaluated at a point
/// - `f_interval`: function evaluated over an interval (inclusion function)
/// - `df_interval`: derivative evaluated over an interval
/// - `initial`: starting interval
/// - `tol`: width tolerance for convergence
/// - `max_iter`: maximum iterations per sub-interval
pub fn interval_newton(
    f: &dyn Fn(f64) -> f64,
    _f_interval: &dyn Fn(Interval) -> Interval,
    df_interval: &dyn Fn(Interval) -> Interval,
    initial: Interval,
    tol: f64,
    max_iter: usize,
) -> IntervalNewtonResult {
    let mut stack: Vec<(Interval, usize)> = vec![(initial, max_iter)];
    let mut results: Vec<Interval> = Vec::new();
    while let Some((mut x, iters_left)) = stack.pop() {
        let mut converged = false;
        for _ in 0..iters_left {
            if x.width() < tol {
                results.push(x);
                converged = true;
                break;
            }
            let m = x.midpoint();
            let fm = f(m);
            let df = df_interval(x);
            if df.contains(0.0) {
                let (left, right) = x.bisect();
                let sub_iter = iters_left / 2;
                if sub_iter > 0 {
                    stack.push((left, sub_iter));
                    stack.push((right, sub_iter));
                }
                converged = true;
                break;
            }
            let correction = Interval::point(fm) / df;
            let newton_step = Interval::point(m) - correction;
            match Interval::intersection(x, newton_step) {
                Some(new_x) => {
                    if new_x.width() >= x.width() * 0.99 {
                        let (left, right) = x.bisect();
                        let sub_iter = iters_left / 2;
                        if sub_iter > 0 {
                            stack.push((left, sub_iter));
                            stack.push((right, sub_iter));
                        }
                        converged = true;
                        break;
                    }
                    x = new_x;
                }
                None => {
                    converged = true;
                    break;
                }
            }
        }
        if !converged {
            results.push(x);
        }
    }
    if results.is_empty() {
        IntervalNewtonResult::NoRoot
    } else if results.len() == 1 {
        IntervalNewtonResult::Root(
            results
                .into_iter()
                .next()
                .expect("results has exactly one element"),
        )
    } else {
        IntervalNewtonResult::Multiple(results)
    }
}
/// Find all root-containing sub-intervals using bisection with an inclusion function.
///
/// Subdivides `initial` into pieces of width ≤ `tol`, keeping only those
/// sub-intervals where `f_interval` contains zero.
pub fn interval_bisection<F>(
    f_interval: F,
    initial: Interval,
    tol: f64,
    max_depth: usize,
) -> Vec<Interval>
where
    F: Fn(Interval) -> Interval,
{
    let mut stack = vec![(initial, 0usize)];
    let mut results = Vec::new();
    while let Some((iv, depth)) = stack.pop() {
        let fiv = f_interval(iv);
        if !fiv.contains(0.0) {
            continue;
        }
        if iv.width() <= tol || depth >= max_depth {
            results.push(iv);
            continue;
        }
        let (left, right) = iv.bisect();
        stack.push((left, depth + 1));
        stack.push((right, depth + 1));
    }
    results
}
/// Global counter for noise symbol IDs.
pub(super) static NOISE_COUNTER: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
/// Allocate a new noise symbol ID.
pub(super) fn new_noise_id() -> usize {
    NOISE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}
/// Simple constraint propagation over interval boxes.
///
/// An `IntervalConstraint` narrows an interval box by propagating
/// known relationships between variables.
pub trait IntervalConstraint {
    /// Try to narrow the given interval box. Returns `true` if narrowing occurred.
    fn propagate(&self, box_: &mut [Interval]) -> bool;
}
/// Run constraint propagation to a fixed point (or max iterations).
pub fn propagate_constraints(
    constraints: &[&dyn IntervalConstraint],
    box_: &mut [Interval],
    max_iter: usize,
) -> usize {
    let mut total_iters = 0;
    for _ in 0..max_iter {
        let mut any_changed = false;
        for c in constraints {
            if c.propagate(box_) {
                any_changed = true;
            }
        }
        total_iters += 1;
        if !any_changed {
            break;
        }
    }
    total_iters
}
/// Conservative interval bound for `tan(a)`.
///
/// Returns `[-∞, +∞]` if the interval contains any odd multiple of `π/2`
/// (where tan is discontinuous), otherwise evaluates endpoint bounds.
pub fn tan_interval(a: Interval) -> Interval {
    use std::f64::consts::FRAC_PI_2;
    let lo_k = (a.lo / std::f64::consts::PI - 0.5).ceil() as i64;
    let hi_k = (a.hi / std::f64::consts::PI - 0.5).floor() as i64;
    if lo_k <= hi_k {
        return Interval {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        };
    }
    let t_lo = a.lo.tan();
    let t_hi = a.hi.tan();
    let _ = FRAC_PI_2;
    Interval {
        lo: t_lo.min(t_hi),
        hi: t_lo.max(t_hi),
    }
}
/// Conservative interval for `atan(a)`.
///
/// `atan` is monotone on all of ℝ so `[atan(lo), atan(hi)]`.
pub fn atan_interval(a: Interval) -> Interval {
    Interval {
        lo: a.lo.atan(),
        hi: a.hi.atan(),
    }
}
/// Interval min: `min([a.lo,a.hi], [b.lo,b.hi]) = [min(a.lo,b.lo), min(a.hi,b.hi)]`.
pub fn interval_min(a: Interval, b: Interval) -> Interval {
    Interval {
        lo: a.lo.min(b.lo),
        hi: a.hi.min(b.hi),
    }
}
/// Interval max: `max([a.lo,a.hi], [b.lo,b.hi]) = [max(a.lo,b.lo), max(a.hi,b.hi)]`.
pub fn interval_max(a: Interval, b: Interval) -> Interval {
    Interval {
        lo: a.lo.max(b.lo),
        hi: a.hi.max(b.hi),
    }
}
/// Returns a tighter interval bound for `f` on `[a.lo, a.hi]` by evaluating
/// at `n` equally-spaced sample points and taking their hull.
///
/// This is a purely rigorous-*sampling* approach — it gives an inner bound
/// (the true enclosure is at least as wide), but is useful for comparison or
/// for functions where monotonicity is hard to determine analytically.
pub fn sample_enclosure(f: impl Fn(f64) -> f64, a: Interval, n: usize) -> Interval {
    assert!(n >= 2, "need at least 2 sample points");
    let step = a.width() / (n - 1) as f64;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for i in 0..n {
        let x = a.lo + i as f64 * step;
        let y = f(x);
        if y < lo {
            lo = y;
        }
        if y > hi {
            hi = y;
        }
    }
    Interval { lo, hi }
}
/// Recursively splits `iv` until all sub-intervals have width ≤ `max_width`.
///
/// Returns the list of sub-intervals.
pub fn subdivide(iv: Interval, max_width: f64) -> Vec<Interval> {
    if iv.width() <= max_width {
        return vec![iv];
    }
    let (left, right) = iv.bisect();
    let mut result = subdivide(left, max_width);
    result.extend(subdivide(right, max_width));
    result
}
/// Evaluates `f` on each sub-interval and returns the hull of all results.
///
/// Achieves tighter enclosure than evaluating `f` on the original interval
/// when `f` suffers from the dependency problem.
pub fn range_by_subdivision(
    f: impl Fn(Interval) -> Interval,
    iv: Interval,
    max_width: f64,
) -> Interval {
    let parts = subdivide(iv, max_width);
    parts.iter().fold(Interval::empty(), |acc, &sub| {
        let fv = f(sub);
        if acc.is_empty() {
            fv
        } else {
            Interval::hull(acc, fv)
        }
    })
}
/// First-order Taylor enclosure of `f` around the midpoint of `iv`.
///
/// `f_mid` = `f(midpoint(iv))` (a point), `df_iv` = bound on `f'` over `iv`.
///
/// Returns `f_mid + df_iv * (iv - midpoint(iv))`.
pub fn taylor_enclosure(f_mid: f64, df_iv: Interval, iv: Interval) -> Interval {
    let m = iv.midpoint();
    let delta = iv - Interval::point(m);
    Interval::point(f_mid) + df_iv * delta
}
/// Krawczyk operator for root verification.
///
/// Given `f : ℝ → ℝ`, approximate Jacobian inverse `y ≈ 1/f'(m)`, and
/// interval `iv`, computes `K = m - y*f(m) + (I - y*f'(iv))*(iv - m)`.
///
/// If `K ⊆ iv` then `iv` contains exactly one root.  If `K ∩ iv = ∅`
/// then `iv` contains no root.
pub fn krawczyk_operator(f_mid: f64, y: f64, df_iv: Interval, iv: Interval) -> Interval {
    let m = iv.midpoint();
    let radius = iv - Interval::point(m);
    let coeff = Interval::point(1.0) - Interval::point(y) * df_iv;
    Interval::point(m - y * f_mid) + coeff * radius
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::interval::AffineForm;
    use crate::interval::BoundConstraint;
    use crate::interval::Interval3;
    use crate::interval::IntervalMatrix;
    use crate::interval::IntervalMatrix2;
    use crate::interval::IntervalMatrix3;
    use crate::interval::IntervalVec;
    use crate::interval::LinearConstraint;
    use crate::interval::TaylorModel;
    #[test]
    fn test_add() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        assert_eq!(a + b, Interval::new(4.0, 6.0));
    }
    #[test]
    fn test_mul_positive() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        assert_eq!(a * b, Interval::new(3.0, 8.0));
    }
    #[test]
    fn test_mul_mixed_signs() {
        let a = Interval::new(-2.0, 1.0);
        let b = Interval::new(-2.0, 1.0);
        let result = a * b;
        assert_eq!(result, Interval::new(-2.0, 4.0));
    }
    #[test]
    fn test_sqrt_interval() {
        let a = Interval::new(4.0, 9.0);
        assert_eq!(sqrt_interval(a), Interval::new(2.0, 3.0));
    }
    #[test]
    fn test_overlaps_true() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 4.0);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_overlaps_false() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_interval3_overlap() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Interval3::from_aabb([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.overlaps_3d(&b));
        let c = Interval3::from_aabb([3.0, 0.0, 0.0], [5.0, 2.0, 2.0]);
        assert!(!a.overlaps_3d(&c));
    }
    #[test]
    fn test_sub() {
        let a = Interval::new(3.0, 5.0);
        let b = Interval::new(1.0, 2.0);
        assert_eq!(a - b, Interval::new(1.0, 4.0));
    }
    #[test]
    fn test_neg() {
        let a = Interval::new(1.0, 3.0);
        assert_eq!(-a, Interval::new(-3.0, -1.0));
    }
    #[test]
    fn test_square() {
        let a = Interval::new(-2.0, 1.0);
        assert_eq!(square(a), Interval::new(0.0, 4.0));
        let b = Interval::new(1.0, 2.0);
        assert_eq!(square(b), Interval::new(1.0, 4.0));
        let c = Interval::new(-2.0, -1.0);
        assert_eq!(square(c), Interval::new(1.0, 4.0));
    }
    #[test]
    fn test_abs_interval() {
        assert_eq!(
            abs_interval(Interval::new(1.0, 3.0)),
            Interval::new(1.0, 3.0)
        );
        assert_eq!(
            abs_interval(Interval::new(-3.0, -1.0)),
            Interval::new(1.0, 3.0)
        );
        assert_eq!(
            abs_interval(Interval::new(-2.0, 3.0)),
            Interval::new(0.0, 3.0)
        );
    }
    #[test]
    fn test_hull() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 5.0);
        assert_eq!(Interval::hull(a, b), Interval::new(1.0, 5.0));
    }
    #[test]
    fn test_width_midpoint_contains() {
        let a = Interval::new(2.0, 6.0);
        assert_eq!(a.width(), 4.0);
        assert_eq!(a.midpoint(), 4.0);
        assert!(a.contains(3.0));
        assert!(!a.contains(7.0));
    }
    #[test]
    fn test_length_sq_bounds() {
        let iv = Interval3::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let bounds = iv.length_sq_bounds();
        assert_eq!(bounds.lo, 0.0);
        assert_eq!(bounds.hi, 3.0);
    }
    #[test]
    fn test_from_to_aabb() {
        let min = [1.0, 2.0, 3.0];
        let max = [4.0, 5.0, 6.0];
        let iv = Interval3::from_aabb(min, max);
        let (out_min, out_max) = iv.to_aabb();
        assert_eq!(out_min, min);
        assert_eq!(out_max, max);
    }
    #[test]
    fn test_div() {
        let a = Interval::new(2.0, 6.0);
        let b = Interval::new(1.0, 2.0);
        let r = a / b;
        assert_eq!(r, Interval::new(1.0, 6.0));
    }
    #[test]
    fn test_scalar_mul() {
        let a = Interval::new(1.0, 3.0);
        assert_eq!(a * 2.0, Interval::new(2.0, 6.0));
        assert_eq!(a * (-1.0), Interval::new(-3.0, -1.0));
        assert_eq!(a * 0.0, Interval::new(0.0, 0.0));
    }
    #[test]
    fn test_reciprocal() {
        let a = Interval::new(2.0, 4.0);
        let r = a.reciprocal();
        assert!((r.lo - 0.25).abs() < 1e-15);
        assert!((r.hi - 0.5).abs() < 1e-15);
    }
    #[test]
    fn test_intersection() {
        let a = Interval::new(1.0, 4.0);
        let b = Interval::new(2.0, 5.0);
        let c = Interval::intersection(a, b).unwrap();
        assert_eq!(c, Interval::new(2.0, 4.0));
        let d = Interval::new(5.0, 6.0);
        assert!(Interval::intersection(a, d).is_none());
    }
    #[test]
    fn test_bisect() {
        let a = Interval::new(0.0, 4.0);
        let (left, right) = a.bisect();
        assert_eq!(left, Interval::new(0.0, 2.0));
        assert_eq!(right, Interval::new(2.0, 4.0));
    }
    #[test]
    fn test_split_at() {
        let a = Interval::new(0.0, 10.0);
        let (l, r) = a.split_at(3.0);
        assert_eq!(l, Interval::new(0.0, 3.0));
        assert_eq!(r, Interval::new(3.0, 10.0));
    }
    #[test]
    fn test_inflate() {
        let a = Interval::new(1.0, 3.0);
        let b = a.inflate(0.5);
        assert_eq!(b, Interval::new(0.5, 3.5));
    }
    #[test]
    fn test_radius() {
        let a = Interval::new(1.0, 5.0);
        assert!((a.radius() - 2.0).abs() < 1e-15);
    }
    #[test]
    fn test_contains_interval() {
        let a = Interval::new(0.0, 10.0);
        let b = Interval::new(2.0, 8.0);
        assert!(a.contains_interval(&b));
        assert!(!b.contains_interval(&a));
    }
    #[test]
    fn test_hull_many() {
        let intervals = vec![
            Interval::new(1.0, 3.0),
            Interval::new(5.0, 7.0),
            Interval::new(-1.0, 2.0),
        ];
        let h = Interval::hull_many(&intervals);
        assert_eq!(h, Interval::new(-1.0, 7.0));
    }
    #[test]
    fn test_empty_interval() {
        let e = Interval::empty();
        assert!(e.is_empty());
        let a = Interval::new(1.0, 2.0);
        assert!(!a.is_empty());
    }
    #[test]
    fn test_pow_interval() {
        let a = Interval::new(-2.0, 3.0);
        assert_eq!(pow_interval(a, 0), Interval::point(1.0));
        assert_eq!(pow_interval(a, 1), a);
        let sq = pow_interval(a, 2);
        assert!((sq.lo - 0.0).abs() < 1e-15);
        assert!((sq.hi - 9.0).abs() < 1e-15);
    }
    #[test]
    fn test_exp_interval() {
        let a = Interval::new(0.0, 1.0);
        let e = exp_interval(a);
        assert!((e.lo - 1.0).abs() < 1e-10);
        assert!((e.hi - std::f64::consts::E).abs() < 1e-10);
    }
    #[test]
    fn test_ln_interval() {
        let a = Interval::new(1.0, std::f64::consts::E);
        let l = ln_interval(a);
        assert!((l.lo - 0.0).abs() < 1e-10);
        assert!((l.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sin_interval_full_period() {
        let a = Interval::new(0.0, 2.0 * std::f64::consts::PI);
        let s = sin_interval(a);
        assert!((s.lo - (-1.0)).abs() < 1e-10);
        assert!((s.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sin_interval_small() {
        let a = Interval::new(0.0, std::f64::consts::FRAC_PI_6);
        let s = sin_interval(a);
        assert!(s.lo <= 0.0 + 1e-10);
        assert!((s.hi - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_cos_interval() {
        let a = Interval::new(0.0, std::f64::consts::PI);
        let c = cos_interval(a);
        assert!((c.lo - (-1.0)).abs() < 1e-10);
        assert!((c.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_interval_bisection_finds_root() {
        let f_iv = |iv: Interval| square(iv) - Interval::point(2.0);
        let roots = interval_bisection(f_iv, Interval::new(0.0, 2.0), 0.01, 20);
        assert!(!roots.is_empty());
        let sqrt2 = 2.0_f64.sqrt();
        let found = roots.iter().any(|r| r.contains(sqrt2));
        assert!(found, "sqrt(2) not found in any root interval: {:?}", roots);
    }
    #[test]
    fn test_interval_bisection_no_root() {
        let f_iv = |iv: Interval| square(iv) + Interval::point(1.0);
        let roots = interval_bisection(f_iv, Interval::new(-5.0, 5.0), 0.01, 20);
        assert!(roots.is_empty(), "should find no roots, got {:?}", roots);
    }
    #[test]
    fn test_interval_newton_basic() {
        let f = |x: f64| x * x - 4.0;
        let f_iv = |iv: Interval| square(iv) - Interval::point(4.0);
        let df_iv = |iv: Interval| iv * 2.0;
        let result = interval_newton(&f, &f_iv, &df_iv, Interval::new(1.0, 3.0), 1e-10, 50);
        match result {
            IntervalNewtonResult::Root(r) => {
                assert!(r.contains(2.0), "root interval {:?} should contain 2.0", r);
                assert!(
                    r.width() < 1e-8,
                    "root should be tight, width={}",
                    r.width()
                );
            }
            other => panic!("expected Root, got {:?}", other),
        }
    }
    #[test]
    fn test_interval_vec_operations() {
        let a = IntervalVec::from_slice(&[Interval::new(1.0, 2.0), Interval::new(3.0, 4.0)]);
        let b = IntervalVec::from_slice(&[Interval::new(0.5, 1.5), Interval::new(2.0, 3.0)]);
        let sum = a.add(&b);
        assert_eq!(sum.data[0], Interval::new(1.5, 3.5));
        assert_eq!(sum.data[1], Interval::new(5.0, 7.0));
        let diff = a.sub(&b);
        assert_eq!(diff.data[0], Interval::new(-0.5, 1.5));
    }
    #[test]
    fn test_interval_vec_dot() {
        let a = IntervalVec::from_slice(&[Interval::point(1.0), Interval::point(2.0)]);
        let b = IntervalVec::from_slice(&[Interval::point(3.0), Interval::point(4.0)]);
        let d = a.dot(&b);
        assert!((d.lo - 11.0).abs() < 1e-15);
        assert!((d.hi - 11.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_vec_hull() {
        let a = IntervalVec::from_slice(&[Interval::new(1.0, 3.0)]);
        let b = IntervalVec::from_slice(&[Interval::new(2.0, 5.0)]);
        let h = a.hull(&b);
        assert_eq!(h.data[0], Interval::new(1.0, 5.0));
    }
    #[test]
    fn test_interval_vec_max_width() {
        let v = IntervalVec::from_slice(&[
            Interval::new(0.0, 1.0),
            Interval::new(0.0, 3.0),
            Interval::new(0.0, 2.0),
        ]);
        assert!((v.max_width() - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_vec_midpoint() {
        let v = IntervalVec::from_slice(&[Interval::new(0.0, 4.0), Interval::new(2.0, 6.0)]);
        let mid = v.midpoint();
        assert!((mid[0] - 2.0).abs() < 1e-15);
        assert!((mid[1] - 4.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix_identity() {
        let id = IntervalMatrix::identity(3);
        assert_eq!(id.get(0, 0), Interval::point(1.0));
        assert_eq!(id.get(0, 1), Interval::point(0.0));
        assert_eq!(id.get(1, 1), Interval::point(1.0));
    }
    #[test]
    fn test_interval_matrix_mul_vec() {
        let mut m = IntervalMatrix::zeros(2, 2);
        m.set(0, 0, Interval::point(1.0));
        m.set(0, 1, Interval::point(2.0));
        m.set(1, 0, Interval::point(3.0));
        m.set(1, 1, Interval::point(4.0));
        let v = IntervalVec::from_slice(&[Interval::point(1.0), Interval::point(1.0)]);
        let r = m.mul_vec(&v);
        assert!((r.data[0].lo - 3.0).abs() < 1e-15);
        assert!((r.data[1].lo - 7.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix_mul_mat() {
        let id = IntervalMatrix::identity(2);
        let mut m = IntervalMatrix::zeros(2, 2);
        m.set(0, 0, Interval::point(5.0));
        m.set(0, 1, Interval::point(6.0));
        m.set(1, 0, Interval::point(7.0));
        m.set(1, 1, Interval::point(8.0));
        let r = id.mul_mat(&m);
        assert_eq!(r.get(0, 0), Interval::point(5.0));
        assert_eq!(r.get(1, 1), Interval::point(8.0));
    }
    #[test]
    fn test_interval_matrix_transpose() {
        let mut m = IntervalMatrix::zeros(2, 3);
        m.set(0, 2, Interval::point(42.0));
        let t = m.transpose();
        assert_eq!(t.rows, 3);
        assert_eq!(t.cols, 2);
        assert_eq!(t.get(2, 0), Interval::point(42.0));
    }
    #[test]
    fn test_affine_form_constant() {
        let a = AffineForm::constant(3.0);
        let iv = a.to_interval();
        assert!((iv.lo - 3.0).abs() < 1e-15);
        assert!((iv.hi - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_affine_from_interval() {
        let iv = Interval::new(1.0, 5.0);
        let a = AffineForm::from_interval(iv);
        let back = a.to_interval();
        assert!((back.lo - 1.0).abs() < 1e-10);
        assert!((back.hi - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_affine_add_sub() {
        let iv1 = Interval::new(1.0, 3.0);
        let iv2 = Interval::new(2.0, 4.0);
        let a = AffineForm::from_interval(iv1);
        let b = AffineForm::from_interval(iv2);
        let sum = a.add(&b);
        let sum_iv = sum.to_interval();
        assert!(sum_iv.lo <= 3.0 + 1e-10);
        assert!(sum_iv.hi >= 7.0 - 1e-10);
        let diff = a.sub(&b);
        let diff_iv = diff.to_interval();
        assert!(diff_iv.lo <= -3.0 + 1e-10);
        assert!(diff_iv.hi >= 1.0 - 1e-10);
    }
    #[test]
    fn test_affine_correlation_tighter() {
        let iv = Interval::new(1.0, 3.0);
        let a = AffineForm::from_interval(iv);
        let diff = a.sub(&a);
        let diff_iv = diff.to_interval();
        assert!(
            diff_iv.width() < 1e-10,
            "x - x should be tight, got {:?}",
            diff_iv
        );
    }
    #[test]
    fn test_affine_scale() {
        let a = AffineForm::from_interval(Interval::new(1.0, 3.0));
        let scaled = a.scale(2.0);
        let iv = scaled.to_interval();
        assert!((iv.lo - 2.0).abs() < 1e-10);
        assert!((iv.hi - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_affine_mul() {
        let a = AffineForm::from_interval(Interval::new(1.0, 3.0));
        let b = AffineForm::from_interval(Interval::new(2.0, 4.0));
        let prod = a.mul(&b);
        let iv = prod.to_interval();
        assert!(iv.lo <= 2.0 + 1e-10, "lo={} should be <= 2", iv.lo);
        assert!(iv.hi >= 12.0 - 1e-10, "hi={} should be >= 12", iv.hi);
    }
    #[test]
    fn test_constraint_propagation_bound() {
        let mut box_ = vec![Interval::new(-10.0, 10.0), Interval::new(-10.0, 10.0)];
        let c = BoundConstraint {
            idx: 0,
            bounds: Interval::new(0.0, 5.0),
        };
        let changed = c.propagate(&mut box_);
        assert!(changed);
        assert_eq!(box_[0], Interval::new(0.0, 5.0));
        assert_eq!(box_[1], Interval::new(-10.0, 10.0));
    }
    #[test]
    fn test_constraint_propagation_linear() {
        let mut box_ = vec![Interval::new(0.0, 10.0), Interval::new(0.0, 10.0)];
        let c = LinearConstraint {
            coeffs: vec![(0, 1.0), (1, 1.0)],
            rhs: Interval::new(3.0, 3.0),
        };
        c.propagate(&mut box_);
        assert!(box_[0].lo >= -0.01);
        assert!(box_[0].hi <= 3.01);
    }
    #[test]
    fn test_propagate_constraints_fixed_point() {
        let mut box_ = vec![Interval::new(-10.0, 10.0)];
        let c1 = BoundConstraint {
            idx: 0,
            bounds: Interval::new(0.0, 5.0),
        };
        let c2 = BoundConstraint {
            idx: 0,
            bounds: Interval::new(2.0, 8.0),
        };
        let constraints: Vec<&dyn IntervalConstraint> = vec![&c1, &c2];
        let iters = propagate_constraints(&constraints, &mut box_, 100);
        assert!(iters <= 3);
        assert!((box_[0].lo - 2.0).abs() < 1e-10);
        assert!((box_[0].hi - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_interval3_hull() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Interval3::from_aabb([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        let h = Interval3::hull_3d(a, b);
        let (min, max) = h.to_aabb();
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [3.0, 3.0, 3.0]);
    }
    #[test]
    fn test_interval3_intersection() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Interval3::from_aabb([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let c = Interval3::intersection_3d(a, b).unwrap();
        let (min, max) = c.to_aabb();
        assert_eq!(min, [1.0, 1.0, 1.0]);
        assert_eq!(max, [2.0, 2.0, 2.0]);
        let d = Interval3::from_aabb([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        assert!(Interval3::intersection_3d(a, d).is_none());
    }
    #[test]
    fn test_interval3_volume_surface_area() {
        let iv = Interval3::from_aabb([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        assert!((iv.volume() - 24.0).abs() < 1e-15);
        assert!((iv.surface_area() - 52.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval3_center() {
        let iv = Interval3::from_aabb([1.0, 2.0, 3.0], [5.0, 6.0, 7.0]);
        let c = iv.center();
        assert_eq!(c, [3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_interval3_inflate() {
        let iv = Interval3::from_aabb([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let inf = iv.inflate(0.5);
        let (min, max) = inf.to_aabb();
        assert_eq!(min, [0.5, 0.5, 0.5]);
        assert_eq!(max, [2.5, 2.5, 2.5]);
    }
    #[test]
    fn test_interval3_contains_point() {
        let iv = Interval3::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(iv.contains_point([0.5, 0.5, 0.5]));
        assert!(!iv.contains_point([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_tan_interval_monotone() {
        let a = Interval::new(0.0, std::f64::consts::FRAC_PI_4);
        let t = tan_interval(a);
        assert!(t.lo <= 0.0 + 1e-10);
        assert!((t.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_tan_interval_with_pole() {
        let a = Interval::new(0.0, std::f64::consts::PI);
        let t = tan_interval(a);
        assert_eq!(t.lo, f64::NEG_INFINITY);
        assert_eq!(t.hi, f64::INFINITY);
    }
    #[test]
    fn test_atan_interval_monotone() {
        let a = Interval::new(-1.0, 1.0);
        let t = atan_interval(a);
        let expected_lo = (-1.0_f64).atan();
        let expected_hi = 1.0_f64.atan();
        assert!((t.lo - expected_lo).abs() < 1e-12);
        assert!((t.hi - expected_hi).abs() < 1e-12);
    }
    #[test]
    fn test_interval_min() {
        let a = Interval::new(1.0, 5.0);
        let b = Interval::new(3.0, 7.0);
        let m = interval_min(a, b);
        assert_eq!(m, Interval::new(1.0, 5.0));
    }
    #[test]
    fn test_interval_max() {
        let a = Interval::new(1.0, 5.0);
        let b = Interval::new(3.0, 7.0);
        let m = interval_max(a, b);
        assert_eq!(m, Interval::new(3.0, 7.0));
    }
    #[test]
    fn test_interval_matrix2_identity_mul_vec() {
        let id = IntervalMatrix2::identity();
        let v = [Interval::point(3.0), Interval::point(5.0)];
        let r = id.mul_vec2(v);
        assert!((r[0].lo - 3.0).abs() < 1e-15);
        assert!((r[1].lo - 5.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix2_mul_mat_identity() {
        let id = IntervalMatrix2::identity();
        let m = IntervalMatrix2::from_f64(1.0, 2.0, 3.0, 4.0);
        let r = id.mul_mat2(&m);
        assert!((r.row0[0].lo - 1.0).abs() < 1e-15);
        assert!((r.row0[1].lo - 2.0).abs() < 1e-15);
        assert!((r.row1[0].lo - 3.0).abs() < 1e-15);
        assert!((r.row1[1].lo - 4.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix2_mul_vec_general() {
        let m = IntervalMatrix2::from_f64(2.0, 0.0, 0.0, 3.0);
        let v = [Interval::point(1.0), Interval::point(2.0)];
        let r = m.mul_vec2(v);
        assert!((r[0].lo - 2.0).abs() < 1e-15);
        assert!((r[1].lo - 6.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix2_transpose() {
        let m = IntervalMatrix2::from_f64(1.0, 2.0, 3.0, 4.0);
        let t = m.transpose2();
        assert!((t.row0[1].lo - 3.0).abs() < 1e-15);
        assert!((t.row1[0].lo - 2.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix3_identity_mul_vec() {
        let id = IntervalMatrix3::identity();
        let v = [
            Interval::point(1.0),
            Interval::point(2.0),
            Interval::point(3.0),
        ];
        let r = id.mul_vec3(v);
        assert!((r[0].lo - 1.0).abs() < 1e-15);
        assert!((r[1].lo - 2.0).abs() < 1e-15);
        assert!((r[2].lo - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_interval_matrix3_mul_mat_identity() {
        let id = IntervalMatrix3::identity();
        let r = id.mul_mat3(&id);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (r.rows[i][j].lo - expected).abs() < 1e-12,
                    "r[{i}][{j}].lo={}",
                    r.rows[i][j].lo
                );
            }
        }
    }
    #[test]
    fn test_interval_matrix3_transpose() {
        let id = IntervalMatrix3::identity();
        let t = id.transpose3();
        for i in 0..3 {
            for j in 0..3 {
                assert!((t.rows[i][j].lo - id.rows[i][j].lo).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_subdivide_width_limit() {
        let iv = Interval::new(0.0, 1.0);
        let parts = subdivide(iv, 0.25);
        for p in &parts {
            assert!(p.width() <= 0.25 + 1e-12, "width={}", p.width());
        }
        let hull = parts.iter().fold(Interval::empty(), |acc, &p| {
            if acc.is_empty() {
                p
            } else {
                Interval::hull(acc, p)
            }
        });
        assert!((hull.lo - 0.0).abs() < 1e-12);
        assert!((hull.hi - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_range_by_subdivision_tighter_than_naive() {
        let iv = Interval::new(0.0, 1.0);
        let f_iv = |x: Interval| x * (Interval::point(1.0) - x);
        let naive = f_iv(iv);
        let tight = range_by_subdivision(f_iv, iv, 0.1);
        assert!(
            tight.hi <= naive.hi + 1e-10,
            "tight.hi={} naive.hi={}",
            tight.hi,
            naive.hi
        );
    }
    #[test]
    fn test_taylor_enclosure_linear_exact() {
        let iv = Interval::new(1.0, 3.0);
        let enc = taylor_enclosure(4.0, Interval::point(2.0), iv);
        assert!((enc.lo - 2.0).abs() < 1e-12);
        assert!((enc.hi - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_enclosure_contains_true_values() {
        let iv = Interval::new(1.0, 3.0);
        let df_iv = Interval::new(2.0, 6.0);
        let enc = taylor_enclosure(4.0, df_iv, iv);
        assert!(enc.lo <= 1.0, "enc.lo={} should be ≤ 1", enc.lo);
        assert!(enc.hi >= 9.0, "enc.hi={} should be ≥ 9", enc.hi);
    }
    #[test]
    fn test_krawczyk_contains_root() {
        let iv = Interval::new(1.5, 2.5);
        let m = iv.midpoint();
        let f_mid = m * m - 4.0;
        let y = 1.0 / (2.0 * m);
        let df_iv = Interval::new(2.0 * 1.5, 2.0 * 2.5);
        let k = krawczyk_operator(f_mid, y, df_iv, iv);
        assert!(
            k.lo <= 2.0 + 1e-10 && k.hi >= 2.0 - 1e-10,
            "Krawczyk interval {:?} should contain root 2",
            k
        );
    }
    #[test]
    fn test_krawczyk_subset_of_iv() {
        let iv = Interval::new(1.5, 2.5);
        let m = iv.midpoint();
        let f_mid = m * m - 4.0;
        let y = 1.0 / (2.0 * m);
        let df_iv = Interval::new(2.0 * 1.5, 2.0 * 2.5);
        let k = krawczyk_operator(f_mid, y, df_iv, iv);
        assert!(
            k.lo >= iv.lo - 1e-10 && k.hi <= iv.hi + 1e-10,
            "K={:?} should be ⊆ iv={:?}",
            k,
            iv
        );
    }
    #[test]
    fn test_sample_enclosure_monotone() {
        let iv = Interval::new(0.0, 1.0);
        let enc = sample_enclosure(|x| x, iv, 10);
        assert!((enc.lo - 0.0).abs() < 1e-10);
        assert!((enc.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sample_enclosure_parabola() {
        let iv = Interval::new(0.0, 3.0);
        let enc = sample_enclosure(|x| x * x, iv, 100);
        assert!(enc.lo <= 0.01);
        assert!((enc.hi - 9.0).abs() < 0.1);
    }
    #[test]
    fn test_interval_add_basic() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        let c = a + b;
        assert!((c.lo - 4.0).abs() < 1e-12);
        assert!((c.hi - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_sub_basic() {
        let a = Interval::new(5.0, 7.0);
        let b = Interval::new(1.0, 2.0);
        let c = a - b;
        assert!((c.lo - 3.0).abs() < 1e-12);
        assert!((c.hi - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_mul_positive() {
        let a = Interval::new(2.0, 3.0);
        let b = Interval::new(4.0, 5.0);
        let c = a * b;
        assert!((c.lo - 8.0).abs() < 1e-12);
        assert!((c.hi - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_mul_mixed_sign() {
        let a = Interval::new(-1.0, 2.0);
        let b = Interval::new(-1.0, 3.0);
        let c = a * b;
        assert!(
            c.lo <= -2.0 - 1e-12 || (c.lo + 3.0).abs() < 1e-12,
            "lo should cover -3, got {}",
            c.lo
        );
        assert!(c.hi >= 6.0 - 1e-12, "hi should be >= 6, got {}", c.hi);
    }
    #[test]
    fn test_interval_div_basic() {
        let a = Interval::new(6.0, 8.0);
        let b = Interval::new(2.0, 4.0);
        let c = a / b;
        assert!(c.lo <= 1.5 + 1e-12);
        assert!(c.hi >= 4.0 - 1e-12);
    }
    #[test]
    fn test_interval_neg() {
        let a = Interval::new(-1.0, 3.0);
        let b = -a;
        assert!((b.lo - (-3.0)).abs() < 1e-12);
        assert!((b.hi - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_contains_point() {
        let a = Interval::new(1.0, 5.0);
        assert!(a.contains(3.0));
        assert!(!a.contains(6.0));
        assert!(a.contains(1.0));
        assert!(a.contains(5.0));
    }
    #[test]
    fn test_interval_overlaps() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 5.0);
        let c = Interval::new(4.0, 6.0);
        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c));
    }
    #[test]
    fn test_interval_hull() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(5.0, 7.0);
        let h = Interval::hull(a, b);
        assert!((h.lo - 1.0).abs() < 1e-12);
        assert!((h.hi - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_intersection_overlapping() {
        let a = Interval::new(1.0, 5.0);
        let b = Interval::new(3.0, 7.0);
        let c = Interval::intersection(a, b).expect("should intersect");
        assert!((c.lo - 3.0).abs() < 1e-12);
        assert!((c.hi - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_intersection_disjoint() {
        let a = Interval::new(1.0, 2.0);
        let b = Interval::new(3.0, 4.0);
        assert!(Interval::intersection(a, b).is_none());
    }
    #[test]
    fn test_interval_powi_even() {
        let a = Interval::new(-2.0, 3.0);
        let p = a.powi(2);
        assert!(p.lo <= 1e-10, "lo={}", p.lo);
        assert!((p.hi - 9.0).abs() < 1e-10, "hi={}", p.hi);
    }
    #[test]
    fn test_interval_powi_odd() {
        let a = Interval::new(-1.0, 2.0);
        let p = a.powi(3);
        assert!((p.lo - (-1.0)).abs() < 1e-10);
        assert!((p.hi - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_interval_sqrt() {
        let a = Interval::new(4.0, 9.0);
        let s = a.sqrt();
        assert!((s.lo - 2.0).abs() < 1e-10);
        assert!((s.hi - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_interval_abs() {
        let a = Interval::new(-3.0, 2.0);
        let b = a.abs();
        assert!(b.lo >= 0.0, "abs lo should be >= 0");
        assert!(b.hi >= 2.0, "abs hi should be >= 2.0, got {}", b.hi);
    }
    #[test]
    fn test_interval_midpoint() {
        let a = Interval::new(2.0, 6.0);
        assert!((a.midpoint() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_radius() {
        let a = Interval::new(2.0, 6.0);
        assert!((a.radius() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_bisect() {
        let a = Interval::new(0.0, 4.0);
        let (left, right) = a.bisect();
        assert!((left.lo - 0.0).abs() < 1e-12);
        assert!((left.hi - 2.0).abs() < 1e-12);
        assert!((right.lo - 2.0).abs() < 1e-12);
        assert!((right.hi - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_inflate() {
        let a = Interval::new(1.0, 3.0);
        let b = a.inflate(0.5);
        assert!((b.lo - 0.5).abs() < 1e-12);
        assert!((b.hi - 3.5).abs() < 1e-12);
    }
    #[test]
    fn test_cos_interval_at_zero() {
        let a = Interval::new(0.0, 0.0);
        let c = cos_interval(a);
        assert!((c.lo - 1.0).abs() < 1e-10);
        assert!((c.hi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_exp_interval_basic() {
        let a = Interval::new(0.0, 1.0);
        let e = exp_interval(a);
        assert!((e.lo - 1.0).abs() < 1e-10);
        assert!((e.hi - std::f64::consts::E).abs() < 1e-8);
    }
    #[test]
    fn test_ln_interval_basic() {
        let a = Interval::new(1.0, std::f64::consts::E);
        let l = ln_interval(a);
        assert!((l.lo - 0.0).abs() < 1e-10);
        assert!((l.hi - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_interval_newton_x_squared_minus_2() {
        let result = interval_newton(
            &|x| x * x - 2.0,
            &|iv: Interval| iv * iv - Interval::point(2.0),
            &|iv: Interval| Interval::point(2.0) * iv,
            Interval::new(1.0, 2.0),
            1e-8,
            50,
        );
        match result {
            IntervalNewtonResult::Root(r) => {
                assert!(
                    r.contains(2.0_f64.sqrt()),
                    "root interval should contain sqrt(2)"
                );
            }
            IntervalNewtonResult::Multiple(v) => {
                let any_contains = v.iter().any(|iv| iv.contains(2.0_f64.sqrt()));
                assert!(any_contains, "one interval should contain sqrt(2)");
            }
            IntervalNewtonResult::NoRoot => panic!("expected a root"),
        }
    }
    #[test]
    fn test_interval3_overlaps() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Interval3::from_aabb([0.5, 0.5, 0.5], [2.0, 2.0, 2.0]);
        assert!(a.overlaps_3d(&b));
    }
    #[test]
    fn test_interval3_no_overlap() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Interval3::from_aabb([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps_3d(&b));
    }
    #[test]
    fn test_interval3_volume_exact() {
        let a = Interval3::from_aabb([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        assert!((a.volume() - 24.0).abs() < 1e-12, "volume={}", a.volume());
    }
    #[test]
    fn test_taylor_model_evaluate() {
        let tm = TaylorModel::new(vec![1.0, 1.0], Interval::new(-0.01, 0.01));
        let ev = tm.evaluate(Interval::point(0.0));
        assert!(ev.contains(1.0), "Taylor model at 0 should contain 1");
    }
    #[test]
    fn test_taylor_model_add() {
        let a = TaylorModel::new(vec![1.0, 2.0], Interval::new(-0.01, 0.01));
        let b = TaylorModel::new(vec![3.0, 1.0], Interval::new(-0.01, 0.01));
        let c = a.add(&b);
        assert!((c.poly[0] - 4.0).abs() < 1e-12, "constant term 1+3=4");
        assert!((c.poly[1] - 3.0).abs() < 1e-12, "linear term 2+1=3");
    }
    #[test]
    fn test_taylor_model_scale() {
        let tm = TaylorModel::new(vec![2.0, 4.0], Interval::new(-0.1, 0.1));
        let scaled = tm.scale(3.0);
        assert!((scaled.poly[0] - 6.0).abs() < 1e-12);
        assert!((scaled.poly[1] - 12.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_bisection_x_squared_minus_1() {
        let roots = interval_bisection(
            |iv: Interval| iv * iv - Interval::point(1.0),
            Interval::new(-2.0, 0.0),
            0.01,
            20,
        );
        assert!(!roots.is_empty(), "should find root of x^2-1 near -1");
        let any_contains = roots.iter().any(|iv| iv.contains(-1.0));
        assert!(any_contains, "a root interval should contain -1");
    }
    #[test]
    fn test_affine_form_from_interval_roundtrip() {
        let iv = Interval::new(1.0, 3.0);
        let af = AffineForm::from_interval(iv);
        let back = af.to_interval();
        assert!(back.lo <= 1.0 + 1e-10, "lo={}", back.lo);
        assert!(back.hi >= 3.0 - 1e-10, "hi={}", back.hi);
    }
    #[test]
    fn test_affine_form_add_constants() {
        let a = AffineForm::constant(2.0);
        let b = AffineForm::constant(3.0);
        let c = a.add(&b);
        let iv = c.to_interval();
        assert!((iv.lo - 5.0).abs() < 1e-12);
        assert!((iv.hi - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_affine_form_scale_constant() {
        let a = AffineForm::constant(4.0);
        let b = a.scale(2.5);
        let iv = b.to_interval();
        assert!((iv.lo - 10.0).abs() < 1e-12);
        assert!((iv.hi - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_matrix2_diagonal_mul_mat() {
        let a = IntervalMatrix2::from_f64(2.0, 0.0, 0.0, 3.0);
        let b = IntervalMatrix2::from_f64(1.0, 0.0, 0.0, 1.0);
        let c = a.mul_mat2(&b);
        assert!((c.row0[0].lo - 2.0).abs() < 1e-12);
        assert!((c.row1[1].lo - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_interval_matrix3_identity_vec3_passthrough() {
        let m = IntervalMatrix3::identity();
        let v = [
            Interval::new(1.0, 2.0),
            Interval::new(3.0, 4.0),
            Interval::new(5.0, 6.0),
        ];
        let result = m.mul_vec3(v);
        for i in 0..3 {
            assert!((result[i].lo - v[i].lo).abs() < 1e-12);
            assert!((result[i].hi - v[i].hi).abs() < 1e-12);
        }
    }
}
