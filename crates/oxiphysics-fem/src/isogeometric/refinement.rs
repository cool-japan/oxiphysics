//! Geometry-preserving degree elevation of NURBS curves.
//!
//! The control points of the elevated curve are **not** samples of the original
//! curve: they are computed exactly (to rounding) from the original control net
//! with Piegl & Tiller, *The NURBS Book* (2nd ed.), Algorithm A5.9
//! (`DegreeElevateCurve`).  The algorithm works on homogeneous control points
//! `(w·x, w·y, w·z, w)`, so rational curves such as exact conic arcs are
//! elevated without any change of shape.
//!
//! A5.9 decomposes the curve into Bezier segments by (virtual) knot insertion,
//! elevates every Bezier segment with the closed-form Bernstein coefficients,
//! and removes the inserted knots again so that the elevated curve keeps the
//! original continuity `C^(p - m)` at every knot of multiplicity `m`.

use super::types::NurbsCurve;
use crate::{Error, Result};

/// Homogeneous control point `(w·x, w·y, w·z, w)`.
type HomogeneousPoint = [f64; 4];

/// Relative tolerance used to decide whether two knot values are equal.
const KNOT_TOLERANCE: f64 = 1e-14;

/// `true` when two knot values coincide within [`KNOT_TOLERANCE`].
fn knots_equal(a: f64, b: f64) -> bool {
    (a - b).abs() <= KNOT_TOLERANCE * a.abs().max(b.abs()).max(1.0)
}

/// Binomial coefficient `C(n, k)` as a floating point number (`0` for `k > n`).
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut value = 1.0;
    for i in 0..k {
        value = value * (n - i) as f64 / (i + 1) as f64;
    }
    value
}

/// Affine combination `a·p + b·q` of two homogeneous points.
fn combine(a: f64, p: HomogeneousPoint, b: f64, q: HomogeneousPoint) -> HomogeneousPoint {
    [
        a * p[0] + b * q[0],
        a * p[1] + b * q[1],
        a * p[2] + b * q[2],
        a * p[3] + b * q[3],
    ]
}

/// Bezier degree elevation coefficients `bezalfs[i][j]` of A5.9.
///
/// Row `i` (`0..=p + t`) gives the weights with which the `p + 1` Bezier points
/// of the original segment contribute to elevated point `i`:
/// `C(p, j)·C(t, i - j) / C(p + t, i)`.
fn bezier_elevation_coefficients(p: usize, t: usize) -> Vec<Vec<f64>> {
    let ph = p + t;
    let ph2 = ph / 2;
    let mut bezalfs = vec![vec![0.0_f64; p + 1]; ph + 1];
    bezalfs[0][0] = 1.0;
    bezalfs[ph][p] = 1.0;
    for (i, row) in bezalfs.iter_mut().enumerate().take(ph2 + 1).skip(1) {
        let inverse = 1.0 / binomial(ph, i);
        for (j, value) in row
            .iter_mut()
            .enumerate()
            .take(p.min(i) + 1)
            .skip(i.saturating_sub(t))
        {
            *value = inverse * binomial(p, j) * binomial(t, i - j);
        }
    }
    // The coefficients are symmetric: row `i` mirrors row `ph - i`.
    for i in (ph2 + 1)..ph {
        let mirrored = bezalfs[ph - i].clone();
        for (j, value) in bezalfs[i]
            .iter_mut()
            .enumerate()
            .take(p.min(i) + 1)
            .skip(i.saturating_sub(t))
        {
            *value = mirrored[p - j];
        }
    }
    bezalfs
}

/// Elevate a clamped homogeneous B-spline curve of degree `p` by `t` degrees.
///
/// `knots` must be a clamped, non-decreasing knot vector of length
/// `ctrl.len() + p + 1` whose interior knots have multiplicity at most `p + 1`
/// (a knot of multiplicity `p + 1` is a jump of the curve and is handled too).
/// Returns the elevated knot vector (every distinct knot gains `t` copies) and
/// the elevated homogeneous control points.
///
/// This is Algorithm A5.9 of Piegl & Tiller with two details made explicit:
/// the first Bezier segment is told apart from a segment that follows a jump
/// (`prev_r == None` versus `prev_r < 0`), so that a jump contributes both of
/// its end points.
fn elevate_homogeneous(
    p: usize,
    knots: &[f64],
    ctrl: &[HomogeneousPoint],
    t: usize,
) -> (Vec<f64>, Vec<HomogeneousPoint>) {
    let n = ctrl.len() - 1;
    let m = n + p + 1;
    let ph = p + t;
    let bezalfs = bezier_elevation_coefficients(p, t);

    // Sizes of the output: every distinct knot gains `t` copies.
    let distinct = 1 + knots
        .windows(2)
        .filter(|w| !knots_equal(w[0], w[1]))
        .count();
    let knots_out = knots.len() + t * distinct;
    let mut uh = vec![0.0_f64; knots_out];
    let mut qw = vec![[0.0_f64; 4]; knots_out - ph - 1];

    let mut kind = ph + 1;
    let mut prev_r: Option<isize> = None;
    let mut a = p;
    let mut b = p + 1;
    let mut cind = 1;
    let mut ua = knots[0];
    qw[0] = ctrl[0];
    uh[..=ph].fill(ua);

    // Bezier points of the current segment, of the next one, and elevated ones.
    let mut bpts: Vec<HomogeneousPoint> = ctrl[..=p].to_vec();
    let mut next_bpts = vec![[0.0_f64; 4]; p + 1];
    let mut ebpts = vec![[0.0_f64; 4]; ph + 1];
    let mut alfs = vec![0.0_f64; p + 1];

    while b < m {
        let group_start = b;
        while b < m && knots_equal(knots[b], knots[b + 1]) {
            b += 1;
        }
        let mul = b - group_start + 1;
        let ub = knots[b];
        let oldr = prev_r.unwrap_or(-1);
        let r = p as isize - mul as isize;

        // First and last elevated Bezier point that must be stored.
        let lbz: usize = match prev_r {
            Some(previous) if previous > 0 => ((previous + 2) / 2) as usize,
            Some(previous) if previous < 0 => 0,
            _ => 1,
        };
        let rbz: usize = if r > 0 {
            ph - ((r + 1) / 2) as usize
        } else {
            ph
        };

        if r > 0 {
            // Insert the knot `ub` r times to obtain the Bezier segment.
            let ru = r as usize;
            let numer = ub - ua;
            for k in ((mul + 1)..=p).rev() {
                alfs[k - mul - 1] = numer / (knots[a + k] - ua);
            }
            for j in 1..=ru {
                let save = ru - j;
                let s = mul + j;
                for k in (s..=p).rev() {
                    let alpha = alfs[k - s];
                    bpts[k] = combine(alpha, bpts[k], 1.0 - alpha, bpts[k - 1]);
                }
                next_bpts[save] = bpts[p];
            }
        }

        // Degree elevate the Bezier segment (only points lbz..=ph are used).
        for (i, ebpt) in ebpts.iter_mut().enumerate().skip(lbz) {
            let mut acc = [0.0_f64; 4];
            for (j, bpt) in bpts
                .iter()
                .enumerate()
                .take(p.min(i) + 1)
                .skip(i.saturating_sub(t))
            {
                let c = bezalfs[i][j];
                for (acc_q, bpt_q) in acc.iter_mut().zip(bpt.iter()) {
                    *acc_q += c * bpt_q;
                }
            }
            *ebpt = acc;
        }

        if oldr > 1 {
            // The knot ua has to be removed `oldr` times again; removal step `tr`
            // sweeps the control points `first ..= last` towards each other.
            let kind_i = kind as isize;
            let den = ub - ua;
            let bet = (ub - uh[kind - 1]) / den;
            for tr in 1..oldr {
                let mut i = kind_i - 1 - tr;
                let mut j = kind_i + tr - 1;
                let mut kj = j - kind_i + 1;
                while j - i > tr {
                    if i < cind as isize {
                        let iu = i as usize;
                        let alf = (ub - uh[iu]) / (ua - uh[iu]);
                        qw[iu] = combine(alf, qw[iu], 1.0 - alf, qw[iu - 1]);
                    }
                    if j >= lbz as isize {
                        let kju = kj as usize;
                        if j - tr <= kind_i - ph as isize + oldr {
                            let gam = (ub - uh[(j - tr) as usize]) / den;
                            ebpts[kju] = combine(gam, ebpts[kju], 1.0 - gam, ebpts[kju + 1]);
                        } else {
                            ebpts[kju] = combine(bet, ebpts[kju], 1.0 - bet, ebpts[kju + 1]);
                        }
                    }
                    i += 1;
                    j -= 1;
                    kj -= 1;
                }
            }
        }

        // Load the knot ua.
        if a != p {
            for _ in 0..(ph as isize - oldr) {
                uh[kind] = ua;
                kind += 1;
            }
        }
        // Load the control points of this segment into the output.
        for &ebpt in &ebpts[lbz..=rbz] {
            qw[cind] = ebpt;
            cind += 1;
        }

        if b < m {
            // Set up for the next pass through the loop.
            let inserted = r.max(0) as usize;
            bpts[..inserted].copy_from_slice(&next_bpts[..inserted]);
            for (j, bpt) in bpts.iter_mut().enumerate().skip(inserted) {
                *bpt = ctrl[b - p + j];
            }
            a = b;
            b += 1;
            ua = ub;
        } else {
            // End knot.
            uh[kind..=kind + ph].fill(ub);
        }
        prev_r = Some(r);
    }
    (uh, qw)
}

/// Check that `curve` is a well-formed clamped NURBS curve that
/// [`elevate_homogeneous`] can process.
fn validate_curve(curve: &NurbsCurve) -> Result<()> {
    let p = curve.basis.degree;
    let knots = &curve.basis.knot;
    let n_ctrl = curve.control_points.len();
    if p == 0 {
        return Err(Error::InvalidInput(
            "degree elevation needs a curve of degree >= 1".to_string(),
        ));
    }
    if n_ctrl < p + 1 {
        return Err(Error::InvalidInput(format!(
            "a degree {p} curve needs at least {} control points, got {n_ctrl}",
            p + 1
        )));
    }
    if curve.weights.len() != n_ctrl {
        return Err(Error::InvalidInput(format!(
            "{} weights for {n_ctrl} control points",
            curve.weights.len()
        )));
    }
    if knots.len() != n_ctrl + p + 1 {
        return Err(Error::InvalidInput(format!(
            "knot vector has {} entries, expected {} (control points + degree + 1)",
            knots.len(),
            n_ctrl + p + 1
        )));
    }
    if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[1] < w[0]) {
        return Err(Error::InvalidInput(
            "knot vector must be finite and non-decreasing".to_string(),
        ));
    }
    let last = knots.len() - 1;
    if !knots[..=p].iter().all(|&k| knots_equal(k, knots[0]))
        || !knots[last - p..]
            .iter()
            .all(|&k| knots_equal(k, knots[last]))
    {
        return Err(Error::InvalidInput(
            "knot vector must be clamped (end knots of multiplicity degree + 1)".to_string(),
        ));
    }
    // Exactly `p + 1` copies at each end (this also rules out an empty domain).
    if knots_equal(knots[p], knots[p + 1]) || knots_equal(knots[last - p - 1], knots[last - p]) {
        return Err(Error::InvalidInput(
            "end knots must have multiplicity exactly degree + 1".to_string(),
        ));
    }
    let mut i = p + 1;
    while i < last - p {
        let mut multiplicity = 1;
        while i + multiplicity < knots.len() && knots_equal(knots[i + multiplicity], knots[i]) {
            multiplicity += 1;
        }
        if multiplicity > p + 1 {
            return Err(Error::InvalidInput(format!(
                "interior knot {} has multiplicity {multiplicity} > degree + 1 = {}",
                knots[i],
                p + 1
            )));
        }
        i += multiplicity;
    }
    if curve.weights.iter().any(|&w| !w.is_finite() || w <= 0.0) {
        return Err(Error::InvalidInput(
            "NURBS weights must be finite and positive".to_string(),
        ));
    }
    if curve
        .control_points
        .iter()
        .any(|c| c.iter().any(|x| !x.is_finite()))
    {
        return Err(Error::InvalidInput(
            "control points must be finite".to_string(),
        ));
    }
    Ok(())
}

/// Degree elevation for NURBS curves.
///
/// Elevates the polynomial degree of a NURBS curve by `times` **without
/// changing the geometry**, using Piegl & Tiller's Algorithm A5.9 on homogeneous
/// control points (see the [module documentation](self)).  The curve must have a
/// clamped knot vector; this is what every constructor in this crate produces.
#[derive(Debug, Clone)]
pub struct DegreeElevation;

impl DegreeElevation {
    /// Elevate the degree of `curve` by `times`, validating the input.
    ///
    /// The result represents exactly the same curve (same shape *and* same
    /// parametrisation): its degree is `p + times`, every distinct knot gains
    /// `times` copies (so the continuity at each knot is unchanged), and the
    /// control points and weights are the exact elevated ones.  Polynomial
    /// curves (all weights equal) keep exactly their weights.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when the curve is not a well-formed clamped NURBS
    /// curve: degree 0, mismatching numbers of control points / weights / knots,
    /// a decreasing or non-finite knot vector, end knots of multiplicity
    /// smaller than `p + 1`, an interior knot of multiplicity above `p + 1`, or
    /// non-positive / non-finite weights.
    pub fn try_elevate(curve: &NurbsCurve, times: usize) -> Result<NurbsCurve> {
        validate_curve(curve)?;
        if times == 0 {
            return Ok(curve.clone());
        }
        let p = curve.basis.degree;
        let uniform_weight = curve.weights.iter().all(|&w| w == curve.weights[0]);
        let homogeneous: Vec<HomogeneousPoint> = curve
            .control_points
            .iter()
            .zip(curve.weights.iter())
            .map(|(c, &w)| {
                // Equal weights cancel: elevate the polynomial curve itself.
                let w = if uniform_weight { 1.0 } else { w };
                [w * c[0], w * c[1], w * c[2], w]
            })
            .collect();
        let (knots, elevated) = elevate_homogeneous(p, &curve.basis.knot, &homogeneous, times);

        let mut points = Vec::with_capacity(elevated.len());
        let mut weights = Vec::with_capacity(elevated.len());
        for q in &elevated {
            if uniform_weight {
                points.push([q[0], q[1], q[2]]);
                weights.push(curve.weights[0]);
            } else {
                let w = q[3];
                if !w.is_finite() || w <= 0.0 {
                    return Err(Error::General(
                        "degree elevation produced a non-positive weight".to_string(),
                    ));
                }
                points.push([q[0] / w, q[1] / w, q[2] / w]);
                weights.push(w);
            }
        }
        Ok(NurbsCurve::new(points, weights, p + times, knots))
    }

    /// Elevate the degree of `curve` by one, validating the input.
    ///
    /// Shorthand for [`try_elevate`](Self::try_elevate) with `times = 1`.
    ///
    /// # Errors
    ///
    /// See [`try_elevate`](Self::try_elevate).
    pub fn try_elevate_once(curve: &NurbsCurve) -> Result<NurbsCurve> {
        Self::try_elevate(curve, 1)
    }

    /// Elevate degree of a NURBS curve by 1.
    ///
    /// Returns a new `NurbsCurve` with degree p+1 representing exactly the same
    /// geometry (Piegl & Tiller A5.9, rational curves included).
    ///
    /// This is the infallible form of [`try_elevate_once`](Self::try_elevate_once).
    /// For a curve that is not well formed (see [`try_elevate`](Self::try_elevate))
    /// every control point and weight of the returned curve is `NaN`, so a failed
    /// elevation can never be mistaken for a valid curve.
    pub fn elevate_once(curve: &NurbsCurve) -> NurbsCurve {
        Self::elevate(curve, 1)
    }

    /// Elevate degree by `times` degrees in a single pass of A5.9.
    ///
    /// Infallible form of [`try_elevate`](Self::try_elevate); an invalid curve
    /// yields a curve whose control points and weights are all `NaN`.
    pub fn elevate(curve: &NurbsCurve, times: usize) -> NurbsCurve {
        match Self::try_elevate(curve, times) {
            Ok(elevated) => elevated,
            Err(_) => {
                let n = curve.control_points.len();
                NurbsCurve {
                    control_points: vec![[f64::NAN; 3]; n],
                    weights: vec![f64::NAN; n],
                    basis: curve.basis.clone(),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isogeometric::KnotRefinement;

    /// Quadratic rational Bezier representing the unit quarter circle.
    fn quarter_circle() -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve::new(
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            vec![1.0, w, 1.0],
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        )
    }

    /// Deterministic pseudo-random number in `[0, 1)`.
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn next_usize(&mut self, lo: usize, hi: usize) -> usize {
            lo + ((self.next_f64() * (hi - lo + 1) as f64) as usize).min(hi - lo)
        }
    }

    fn max_deviation(a: &NurbsCurve, b: &NurbsCurve, samples: &[f64]) -> f64 {
        samples
            .iter()
            .map(|&t| {
                let p = a.point_at(t);
                let q = b.point_at(t);
                ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
            })
            .fold(0.0, f64::max)
    }

    fn uniform_samples(lo: f64, hi: f64, count: usize) -> Vec<f64> {
        (0..count)
            .map(|i| lo + (hi - lo) * i as f64 / (count - 1) as f64)
            .collect()
    }

    #[test]
    fn quarter_circle_elevated_to_cubic_is_exact() {
        let circle = quarter_circle();
        let elevated = DegreeElevation::try_elevate_once(&circle).expect("valid curve");
        assert_eq!(elevated.basis.degree, 3);
        assert_eq!(elevated.control_points.len(), 4);
        assert_eq!(
            elevated.basis.knot,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        );

        let samples = uniform_samples(0.0, 1.0, 101);
        let deviation = max_deviation(&circle, &elevated, &samples);
        assert!(deviation < 1e-12, "geometry changed by {deviation}");
        for &t in &samples {
            let q = elevated.point_at(t);
            let radius = (q[0] * q[0] + q[1] * q[1]).sqrt();
            assert!((radius - 1.0).abs() < 1e-12, "radius {radius} at t = {t}");
            assert_eq!(q[2], 0.0);
        }
        // Parametrisation is preserved as well, not only the shape.
        let mid = elevated.point_at(0.5);
        let reference = circle.point_at(0.5);
        for d in 0..3 {
            assert!((mid[d] - reference[d]).abs() < 1e-14);
        }
    }

    #[test]
    fn elevated_control_points_are_not_curve_samples() {
        // The old implementation placed samples of the curve at uniform
        // parameters in the control net; the exact elevated net of a conic lies
        // *outside* the circle (it is the convex-hull-like control polygon).
        let elevated = DegreeElevation::elevate_once(&quarter_circle());
        let outside = elevated
            .control_points
            .iter()
            .filter(|c| (c[0] * c[0] + c[1] * c[1]).sqrt() > 1.0 + 1e-6)
            .count();
        assert!(outside >= 2, "control polygon must bulge out of the circle");
        let w = elevated.weights.clone();
        // end weights stay 1 and the interior ones are the elevated blend
        assert!((w[0] - 1.0).abs() < 1e-14 && (w[3] - 1.0).abs() < 1e-14);
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let expected_inner = (1.0 / 3.0) * 1.0 + (2.0 / 3.0) * s;
        assert!(
            (w[1] - expected_inner).abs() < 1e-14,
            "{} vs {expected_inner}",
            w[1]
        );
        assert!((w[2] - expected_inner).abs() < 1e-14);
    }

    #[test]
    fn quarter_circle_elevated_repeatedly_stays_a_circle() {
        let circle = quarter_circle();
        let samples = uniform_samples(0.0, 1.0, 101);
        let single_pass = DegreeElevation::elevate(&circle, 3);
        assert_eq!(single_pass.basis.degree, 5);
        let mut step_by_step = circle.clone();
        for _ in 0..3 {
            step_by_step = DegreeElevation::elevate_once(&step_by_step);
        }
        assert!(max_deviation(&circle, &single_pass, &samples) < 1e-12);
        assert!(max_deviation(&circle, &step_by_step, &samples) < 1e-12);
        assert_eq!(
            single_pass.control_points.len(),
            step_by_step.control_points.len()
        );
        for (a, b) in single_pass
            .control_points
            .iter()
            .zip(step_by_step.control_points.iter())
        {
            for d in 0..3 {
                assert!((a[d] - b[d]).abs() < 1e-12, "{a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn zero_elevations_return_the_same_curve() {
        let circle = quarter_circle();
        let same = DegreeElevation::try_elevate(&circle, 0).expect("valid");
        assert_eq!(same.basis.degree, 2);
        assert_eq!(same.control_points, circle.control_points);
        assert_eq!(same.weights, circle.weights);
    }

    #[test]
    fn polynomial_curve_keeps_unit_weights_exactly() {
        let curve = NurbsCurve::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [3.0, 2.0, 1.0],
                [4.0, 0.0, 1.0],
            ],
            vec![1.0; 4],
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        );
        let elevated = DegreeElevation::elevate_once(&curve);
        assert_eq!(elevated.weights, vec![1.0; 5]);
        // Cubic Bezier elevation: Q_i = (i/4) P_{i-1} + (1 - i/4) P_i.
        let expected_q1 = [0.25 * 0.0 + 0.75 * 1.0, 0.25 * 0.0 + 0.75 * 2.0, 0.0];
        let q1 = elevated.control_points[1];
        for d in 0..3 {
            assert!((q1[d] - expected_q1[d]).abs() < 1e-15, "{q1:?}");
        }
        assert_eq!(elevated.control_points[0], curve.control_points[0]);
        assert_eq!(elevated.control_points[4], curve.control_points[3]);
    }

    #[test]
    fn straight_line_stays_straight_and_uniformly_parametrised() {
        let line = NurbsCurve::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0, 1.0],
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        );
        let elevated = DegreeElevation::elevate_once(&line);
        for &t in &uniform_samples(0.0, 1.0, 41) {
            let q = elevated.point_at(t);
            assert!((q[0] - 2.0 * t).abs() < 1e-14, "x({t}) = {}", q[0]);
            assert_eq!(q[1], 0.0);
            assert_eq!(q[2], 0.0);
        }
    }

    #[test]
    fn interior_knots_of_various_multiplicity_are_preserved() {
        // cubic, 9 control points, interior knots of multiplicity 1 and 2,
        // non-trivial weights, 3-D.
        let points: Vec<[f64; 3]> = (0..9)
            .map(|i| {
                let x = i as f64;
                [
                    x * 0.7 + (x * 1.3).sin(),
                    (x * 0.9).cos() * 2.0,
                    (x * 0.37).sin() * x,
                ]
            })
            .collect();
        let weights: Vec<f64> = (0..9).map(|i| 0.6 + 0.35 * ((i * 7 % 5) as f64)).collect();
        let knots = vec![
            0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.55, 0.7, 0.9, 1.0, 1.0, 1.0, 1.0,
        ];
        let curve = NurbsCurve::new(points, weights, 3, knots.clone());
        let samples = uniform_samples(0.0, 1.0, 301);
        for times in 1..=3 {
            let elevated = DegreeElevation::try_elevate(&curve, times).expect("valid curve");
            assert_eq!(elevated.basis.degree, 3 + times);
            let deviation = max_deviation(&curve, &elevated, &samples);
            assert!(deviation < 1e-12, "times {times}: deviation {deviation}");
            // Each distinct knot gained `times` copies: multiplicities unchanged
            // relative to degree (continuity preserved).
            let distinct = {
                let mut d = knots.clone();
                d.dedup();
                d.len()
            };
            assert_eq!(
                elevated.basis.knot.len(),
                knots.len() + times * distinct,
                "knot count"
            );
            assert_eq!(
                elevated.control_points.len(),
                elevated.basis.knot.len() - (3 + times) - 1
            );
        }
    }

    #[test]
    fn elevated_knot_vector_matches_p_refine_knot() {
        let curve = NurbsCurve::new(
            (0..6).map(|i| [i as f64, (i * i) as f64, 0.0]).collect(),
            vec![1.0, 2.0, 1.5, 1.0, 0.8, 1.0],
            2,
            vec![0.0, 0.0, 0.0, 0.3, 0.3, 0.8, 1.0, 1.0, 1.0],
        );
        let elevated = DegreeElevation::elevate_once(&curve);
        let expected = KnotRefinement::p_refine_knot(&curve.basis.knot, 2);
        assert_eq!(elevated.basis.knot, expected);
    }

    #[test]
    fn continuity_at_knots_is_unchanged() {
        // C^1 cubic with a single interior knot of multiplicity 2: after the
        // elevation the first derivative must still be continuous there (and a
        // multiplicity-3 knot, a C^0 corner, must stay a corner).
        let smooth = NurbsCurve::new(
            (0..6)
                .map(|i| {
                    let x = i as f64;
                    [x, (x * 0.8).sin(), 0.0]
                })
                .collect(),
            vec![1.0; 6],
            3,
            vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
        );
        let elevated = DegreeElevation::elevate_once(&smooth);
        let eps = 1e-7;
        let left = elevated.derivative_at(0.5 - eps);
        let right = elevated.derivative_at(0.5 + eps);
        for d in 0..3 {
            assert!(
                (left[d] - right[d]).abs() < 1e-5,
                "tangent jump {left:?} {right:?}"
            );
        }
        let reference = smooth.derivative_at(0.5 - eps);
        for d in 0..3 {
            assert!((left[d] - reference[d]).abs() < 1e-5);
        }

        let corner = NurbsCurve::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 1.0, 0.0],
                [4.0, 0.0, 0.0],
            ],
            vec![1.0; 5],
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        );
        let corner_elevated = DegreeElevation::elevate_once(&corner);
        let l = corner_elevated.derivative_at(0.5 - eps);
        let r = corner_elevated.derivative_at(0.5 + eps);
        let jump: f64 = (0..3).map(|d| (l[d] - r[d]).abs()).sum();
        assert!(jump > 1e-3, "C^0 corner must not be smoothed away");
        assert!(max_deviation(&corner, &corner_elevated, &uniform_samples(0.0, 1.0, 201)) < 1e-12);
    }

    #[test]
    fn curve_with_a_jump_is_elevated_piecewise() {
        // Interior knot of multiplicity p + 1 = 3: two disconnected Beziers.
        let curve = NurbsCurve::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 1.0],
                [4.0, 2.0, 1.0],
                [5.0, 0.0, 1.0],
            ],
            vec![1.0, 2.0, 1.0, 1.0, 0.5, 1.0],
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
        );
        let elevated = DegreeElevation::try_elevate_once(&curve).expect("valid curve");
        assert_eq!(elevated.control_points.len(), 8);
        assert_eq!(
            elevated.basis.knot,
            vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0]
        );
        let samples: Vec<f64> = uniform_samples(0.0, 1.0, 201)
            .into_iter()
            .filter(|t| (t - 0.5).abs() > 1e-9)
            .collect();
        assert!(max_deviation(&curve, &elevated, &samples) < 1e-12);
    }

    #[test]
    fn random_curves_keep_their_shape() {
        let mut rng = Lcg(0x5eed_1234_abcd_ef01);
        for case in 0..400 {
            let p = rng.next_usize(1, 5);
            let interior_groups = rng.next_usize(0, 4);
            let allow_jump = rng.next_f64() < 0.3;
            let mut knots = vec![-1.5; p + 1];
            let mut value = -1.5;
            for _ in 0..interior_groups {
                value += 0.1 + rng.next_f64();
                let max_multiplicity = if allow_jump { p + 1 } else { p };
                let multiplicity = rng.next_usize(1, max_multiplicity);
                knots.extend(std::iter::repeat_n(value, multiplicity));
            }
            value += 0.1 + rng.next_f64();
            knots.extend(std::iter::repeat_n(value, p + 1));
            let n_ctrl = knots.len() - p - 1;
            let points: Vec<[f64; 3]> = (0..n_ctrl)
                .map(|_| {
                    [
                        rng.next_f64() * 4.0 - 2.0,
                        rng.next_f64() * 4.0 - 2.0,
                        rng.next_f64() * 4.0 - 2.0,
                    ]
                })
                .collect();
            let weights: Vec<f64> = (0..n_ctrl).map(|_| 0.3 + 2.0 * rng.next_f64()).collect();
            let lo = knots[0];
            let hi = knots[knots.len() - 1];
            let curve = NurbsCurve::new(points, weights, p, knots.clone());
            let times = rng.next_usize(1, 3);
            let elevated = DegreeElevation::try_elevate(&curve, times).expect("valid curve");
            assert_eq!(elevated.basis.degree, p + times);
            assert_eq!(
                elevated.basis.knot.len(),
                elevated.control_points.len() + p + times + 1
            );
            let samples: Vec<f64> = (0..80)
                .map(|i| lo + (hi - lo) * (i as f64 + 0.37) / 80.0)
                .filter(|t| knots.iter().all(|k| (k - t).abs() > 1e-7))
                .collect();
            let deviation = max_deviation(&curve, &elevated, &samples);
            assert!(
                deviation < 1e-10,
                "case {case}: p {p}, times {times}: {deviation}"
            );
        }
    }

    #[test]
    fn invalid_curves_are_rejected() {
        let valid = quarter_circle();

        let mut degree_zero = valid.clone();
        degree_zero.basis.degree = 0;
        assert!(DegreeElevation::try_elevate_once(&degree_zero).is_err());

        let mut bad_weights = valid.clone();
        bad_weights.weights.pop();
        assert!(DegreeElevation::try_elevate_once(&bad_weights).is_err());

        let mut negative_weight = valid.clone();
        negative_weight.weights[1] = -1.0;
        assert!(DegreeElevation::try_elevate_once(&negative_weight).is_err());

        let mut decreasing = valid.clone();
        decreasing.basis.knot = vec![0.0, 0.0, 0.0, 1.0, 0.5, 1.0];
        assert!(DegreeElevation::try_elevate_once(&decreasing).is_err());

        let mut unclamped = valid.clone();
        unclamped.basis.knot = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        assert!(DegreeElevation::try_elevate_once(&unclamped).is_err());

        let mut wrong_length = valid.clone();
        wrong_length.basis.knot = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        assert!(DegreeElevation::try_elevate_once(&wrong_length).is_err());

        let over_multiple = NurbsCurve::new(
            vec![[0.0; 3]; 7],
            vec![1.0; 7],
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
        );
        assert!(DegreeElevation::try_elevate_once(&over_multiple).is_err());
    }

    #[test]
    fn infallible_elevation_poisons_invalid_curves_with_nan() {
        let mut bad = quarter_circle();
        bad.weights[0] = 0.0;
        let result = DegreeElevation::elevate_once(&bad);
        assert!(
            result
                .control_points
                .iter()
                .all(|c| c.iter().all(|x| x.is_nan()))
        );
        assert!(result.weights.iter().all(|w| w.is_nan()));
        assert!(result.point_at(0.5)[0].is_nan());
    }
}
