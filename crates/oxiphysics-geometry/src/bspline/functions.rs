//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::{Error, Result};

/// Dot product of two 3D vectors.
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Central finite difference of a 3D vector function at x.
pub(super) fn finite_diff_3d<F: Fn(f64) -> [f64; 3]>(f: F, x: f64, h: f64) -> [f64; 3] {
    let f0 = f(x - h);
    let f1 = f(x + h);
    [
        (f1[0] - f0[0]) / (2.0 * h),
        (f1[1] - f0[1]) / (2.0 * h),
        (f1[2] - f0[2]) / (2.0 * h),
    ]
}
/// Solve `A X = B` for three right-hand sides at once, where `A` is a dense
/// `n × n` matrix and the columns of `B` are the x, y and z coordinates of
/// `rhs`.
///
/// Uses Gaussian elimination with partial (row) pivoting on heap-allocated
/// augmented rows, so any `n` is supported.  The system is declared singular
/// when a pivot falls below `1e-14 · max|A|`, in which case an error is
/// returned instead of silently producing zeros.
///
/// Returns one solution triple `[x, y, z]` per unknown.
pub(super) fn solve_3x_systems(
    a: &[Vec<f64>],
    rhs: &[[f64; 3]],
    n: usize,
) -> Result<Vec<[f64; 3]>> {
    if n == 0 {
        return Ok(Vec::new());
    }
    if a.len() < n {
        return Err(Error::length_mismatch(n, a.len()));
    }
    if rhs.len() < n {
        return Err(Error::length_mismatch(n, rhs.len()));
    }
    if let Some(short_row) = a[..n].iter().find(|row| row.len() < n) {
        return Err(Error::length_mismatch(n, short_row.len()));
    }
    let width = n + 3;
    let mut aug: Vec<Vec<f64>> = a[..n]
        .iter()
        .zip(rhs[..n].iter())
        .map(|(row, rhs_row)| {
            let mut aug_row = Vec::with_capacity(width);
            aug_row.extend_from_slice(&row[..n]);
            aug_row.extend_from_slice(rhs_row);
            aug_row
        })
        .collect();
    let scale = aug
        .iter()
        .flat_map(|row| row[..n].iter())
        .fold(0.0_f64, |acc, v| acc.max(v.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return Err(Error::general(
            "singular system: coefficient matrix is zero or contains non-finite entries",
        ));
    }
    let pivot_tolerance = 1e-14 * scale;
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (offset, aug_row) in aug[col + 1..n].iter().enumerate() {
            let candidate = aug_row[col].abs();
            if candidate > max_val {
                max_val = candidate;
                max_row = col + 1 + offset;
            }
        }
        if max_val <= pivot_tolerance || max_val.is_nan() {
            return Err(Error::general(format!(
                "singular system: pivot {max_val:.3e} in column {col} is below the tolerance {pivot_tolerance:.3e}"
            )));
        }
        aug.swap(col, max_row);
        let (upper, lower) = aug.split_at_mut(col + 1);
        let pivot_row = &upper[col];
        let pivot = pivot_row[col];
        for row in lower.iter_mut() {
            let factor = row[col] / pivot;
            if factor != 0.0 {
                for (target, source) in row[col..].iter_mut().zip(pivot_row[col..].iter()) {
                    *target -= factor * source;
                }
            }
        }
    }
    let mut x = vec![[0.0_f64; 3]; n];
    for i in (0..n).rev() {
        let row = &aug[i];
        let diag = row[i];
        for d in 0..3 {
            let mut sum = row[n + d];
            for (j, x_j) in x.iter().enumerate().skip(i + 1) {
                sum -= row[j] * x_j[d];
            }
            x[i][d] = sum / diag;
        }
    }
    if x.iter().any(|v| v.iter().any(|c| !c.is_finite())) {
        return Err(Error::general(
            "singular system: elimination produced non-finite values",
        ));
    }
    Ok(x)
}
#[cfg(test)]
mod tests {
    use crate::bspline::*;
    use std::f64::consts::PI;
    #[test]
    fn test_knot_vector_new() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]);
        assert_eq!(kv.len(), 5);
    }
    #[test]
    fn test_knot_vector_clamped_uniform_length() {
        let kv = KnotVector::clamped_uniform(5, 3);
        assert_eq!(kv.len(), 9);
    }
    #[test]
    fn test_knot_vector_domain() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let (a, b) = kv.domain();
        assert!((a - 0.0).abs() < 1e-14);
        assert!((b - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_knot_vector_find_span_basic() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let span = kv.find_span(0.5, 2, 3);
        assert_eq!(span, 2);
    }
    #[test]
    fn test_knot_vector_find_span_at_start() {
        let kv = KnotVector::clamped_uniform(4, 2);
        let span = kv.find_span(0.0, 2, 4);
        assert!(span >= 2);
    }
    #[test]
    fn test_knot_vector_find_span_at_end() {
        let kv = KnotVector::clamped_uniform(4, 2);
        let span = kv.find_span(1.0, 2, 4);
        assert!(span >= 2);
    }
    #[test]
    fn test_basis_eval_partition_of_unity() {
        let basis = BsplineBasis::clamped(3, 6);
        for t in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            let vals = basis.eval_all(t);
            let sum: f64 = vals.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "PoU failed at t={t}: sum={sum}");
        }
    }
    #[test]
    fn test_basis_eval_nonnegativity() {
        let basis = BsplineBasis::clamped(3, 7);
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let vals = basis.eval_all(t);
            for (i, &v) in vals.iter().enumerate() {
                assert!(v >= -1e-14, "N_{i}(t={t}) = {v} < 0");
            }
        }
    }
    #[test]
    fn test_basis_eval_endpoint_interpolation() {
        let basis = BsplineBasis::clamped(3, 5);
        let v0 = basis.eval_all(0.0);
        let v1 = basis.eval_all(1.0);
        assert!((v0[0] - 1.0).abs() < 1e-10, "N_0(0) = {}", v0[0]);
        assert!((v1[4] - 1.0).abs() < 1e-10, "N_4(1) = {}", v1[4]);
    }
    #[test]
    fn test_basis_eval_deriv_first_order() {
        let basis = BsplineBasis::clamped(3, 6);
        for t in [0.2, 0.4, 0.6, 0.8] {
            let _p = basis.degree;
            let (span, ders) = basis.eval_nonzero_derivs(t, 1);
            let sum_d1: f64 = ders[1].iter().sum();
            assert!(
                sum_d1.abs() < 1e-8,
                "sum of d/dt N_i at t={t} span={span} = {sum_d1}"
            );
        }
    }
    #[test]
    fn test_greville_abscissae_count() {
        let basis = BsplineBasis::clamped(3, 6);
        let xi = basis.greville_abscissae();
        assert_eq!(xi.len(), 6);
    }
    #[test]
    fn test_greville_abscissae_range() {
        let basis = BsplineBasis::clamped(3, 8);
        let xi = basis.greville_abscissae();
        for &x in &xi {
            assert!(
                (0.0..=1.0).contains(&x),
                "Greville abscissa {x} out of [0,1]"
            );
        }
    }
    #[test]
    fn test_bspline_curve_endpoint_interpolation() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps.clone());
        let p0 = curve.eval(0.0);
        let p1 = curve.eval(1.0);
        assert!((p0[0] - cps[0][0]).abs() < 1e-10);
        assert!((p1[0] - cps[3][0]).abs() < 1e-10);
    }
    #[test]
    fn test_bspline_curve_midpoint_finite() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let pt = curve.eval(0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite() && pt[2].is_finite());
    }
    #[test]
    fn test_bspline_curve_tangent_unit() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let t = curve.tangent(0.5);
        let mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-8, "tangent not unit: mag={mag}");
    }
    #[test]
    fn test_bspline_curve_straight_line_curvature_zero() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let kappa = curve.curvature(0.5);
        assert!(kappa < 1e-6, "straight line curvature = {kappa}");
    }
    #[test]
    fn test_bspline_curve_arc_length_positive() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let len = curve.arc_length(0.0, 1.0, 20);
        assert!(len > 0.0);
    }
    #[test]
    fn test_bspline_curve_arc_length_straight_line() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let len = curve.arc_length(0.0, 1.0, 20);
        assert!((len - 2.0).abs() < 1e-4, "straight line arc length = {len}");
    }
    #[test]
    fn test_bspline_curve_sample_count() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let pts = curve.sample(10);
        assert_eq!(pts.len(), 10);
    }
    #[test]
    fn test_bspline_curve_bounding_box() {
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BsplineCurve::clamped(2, cps);
        let (lo, hi) = curve.bounding_box(50);
        assert!(lo[0] <= hi[0]);
        assert!(lo[1] <= hi[1]);
    }
    #[test]
    fn test_bspline_surface_bilinear_corners() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let c00 = surf.eval(0.0, 0.0);
        let c11 = surf.eval(1.0, 1.0);
        assert!((c00[0] - 0.0).abs() < 1e-10 && (c00[1] - 0.0).abs() < 1e-10);
        assert!((c11[0] - 1.0).abs() < 1e-10 && (c11[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bspline_surface_normal_unit() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let n = surf.normal(0.5, 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-8, "normal magnitude = {mag}");
    }
    #[test]
    fn test_bspline_surface_flat_plane_normal_z() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let n = surf.normal(0.3, 0.7);
        assert!(
            n[2].abs() > 0.99,
            "flat plane normal z-component = {}",
            n[2]
        );
    }
    #[test]
    fn test_bspline_surface_sample_dimensions() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let pts = surf.sample(5, 7);
        assert_eq!(pts.len(), 5);
        assert_eq!(pts[0].len(), 7);
    }
    #[test]
    fn test_bspline_surface_first_fundamental_form_positive() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let (e, _f, g) = surf.first_fundamental_form(0.5, 0.5);
        assert!(e > 0.0 && g > 0.0);
    }
    #[test]
    fn test_nurbs_curve_circle_arc_radius() {
        let center = [0.0, 0.0, 0.0];
        let radius = 2.0;
        let nurbs = NurbsCurve::circle_arc(center, radius, 0.0, PI);
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = nurbs.eval(t);
            let dist = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
            assert!((dist - radius).abs() < 1e-6, "dist = {dist}, t = {t}");
        }
    }
    #[test]
    fn test_nurbs_curve_circle_full_radius() {
        let center = [1.0, 2.0, 0.0];
        let radius = 1.5;
        let nurbs = NurbsCurve::circle_arc(center, radius, 0.0, 2.0 * PI);
        for t in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let pt = nurbs.eval(t);
            let dx = pt[0] - center[0];
            let dy = pt[1] - center[1];
            let dist = (dx * dx + dy * dy).sqrt();
            assert!((dist - radius).abs() < 1e-5, "dist = {dist} at t = {t}");
        }
    }
    #[test]
    fn test_nurbs_curve_eval_finite() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0];
        let nurbs = NurbsCurve::new(2, kv, cps, weights);
        let pt = nurbs.eval(0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite());
    }
    #[test]
    fn test_nurbs_curve_arc_length_positive() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 1.0, 0.0, PI);
        let len = nurbs.arc_length(0.0, 1.0, 20);
        assert!(len > 0.0);
        assert!(
            (len - PI).abs() < 0.1,
            "arc length = {len}, expected π ≈ {PI}"
        );
    }
    #[test]
    fn test_nurbs_curve_update_control_point() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let cps = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0, 1.0];
        let mut nurbs = NurbsCurve::new(2, kv, cps, weights);
        nurbs.update_control_point(1, [1.0, 3.0, 0.0], 0.5);
        assert!((nurbs.control_points[1][1] - 3.0).abs() < 1e-10);
        assert!((nurbs.weights[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_nurbs_curve_curvature_finite() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 2.0, 0.0, PI);
        let kappa = nurbs.curvature(0.5);
        assert!(kappa.is_finite());
        assert!(kappa >= 0.0);
    }
    #[test]
    fn test_nurbs_curve_sample_count() {
        let nurbs = NurbsCurve::circle_arc([0.0, 0.0, 0.0], 1.0, 0.0, PI);
        let pts = nurbs.sample(15);
        assert_eq!(pts.len(), 15);
    }
    #[test]
    fn test_nurbs_surface_eval_finite() {
        let kv = KnotVector::new(vec![0.0, 0.0, 0.5, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 0.0]],
            vec![[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 2.0, 0.0]],
        ];
        let weights = vec![
            vec![1.0, 1.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![1.0, 1.0, 1.0],
        ];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let pt = surf.eval(0.5, 0.5);
        assert!(pt[0].is_finite() && pt[1].is_finite() && pt[2].is_finite());
    }
    #[test]
    fn test_nurbs_surface_normal_unit() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let n = surf.normal(0.5, 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((mag - 1.0).abs() < 1e-6, "normal magnitude = {mag}");
    }
    #[test]
    fn test_nurbs_surface_sample_shape() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let pts = surf.sample(4, 5);
        assert_eq!(pts.len(), 4);
        assert_eq!(pts[0].len(), 5);
    }
    #[test]
    fn test_chord_length_parameterization_endpoints() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let params = BsplineFitting::chord_length_parameterization(&pts);
        assert!((params[0] - 0.0).abs() < 1e-10);
        assert!((params[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_chord_length_parameterization_monotonic() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let params = BsplineFitting::chord_length_parameterization(&pts);
        for i in 1..params.len() {
            assert!(params[i] >= params[i - 1], "params not monotonic at {i}");
        }
    }
    #[test]
    fn test_centripetal_parameterization_endpoints() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 0.0]];
        let params = BsplineFitting::centripetal_parameterization(&pts);
        assert!((params[0] - 0.0).abs() < 1e-10);
        assert!((params[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_select_knots_by_averaging_length() {
        let params = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 4, 3);
        assert_eq!(kv.len(), 8, "knot vector length = {}", kv.len());
    }
    /// Interpolation knots are the averages `u_{p+j} = (t_j + … + t_{j+p-1}) / p`.
    #[test]
    fn test_select_knots_interpolation_uses_parameter_averages() {
        let params = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 5, 3);
        let expected = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(kv.knots.len(), expected.len());
        for (got, want) in kv.knots.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-14, "{:?}", kv.knots);
        }
        let params = vec![0.0, 0.1, 0.3, 0.6, 0.8, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 6, 3);
        assert!((kv.knots[4] - (0.1 + 0.3 + 0.6) / 3.0).abs() < 1e-14);
        assert!((kv.knots[5] - (0.3 + 0.6 + 0.8) / 3.0).abs() < 1e-14);
        // Degree 1 averages a single parameter: the knots are the parameters.
        let kv = BsplineFitting::select_knots_by_averaging(&params, 6, 1);
        for (j, &param) in params.iter().enumerate().take(5).skip(1) {
            assert!((kv.knots[1 + j] - param).abs() < 1e-14);
        }
    }
    /// Least-squares knots put at least one parameter into every non-empty
    /// knot span, which keeps the normal equations nonsingular.
    #[test]
    fn test_select_knots_least_squares_every_span_has_a_parameter() {
        for &(n_pts, n_ctrl, degree) in &[
            (11usize, 6usize, 3usize),
            (50, 10, 3),
            (30, 12, 5),
            (9, 8, 2),
        ] {
            // Mildly non-uniform parameters in [0, 1].
            let params: Vec<f64> = (0..n_pts)
                .map(|i| {
                    let u = i as f64 / (n_pts - 1) as f64;
                    u * u * (3.0 - 2.0 * u)
                })
                .collect();
            let kv = BsplineFitting::select_knots_by_averaging(&params, n_ctrl, degree);
            assert_eq!(kv.knots.len(), n_ctrl + degree + 1);
            for w in kv.knots.windows(2) {
                assert!(w[1] >= w[0]);
            }
            for span in kv.knots.windows(2) {
                if span[1] - span[0] <= 0.0 {
                    continue;
                }
                let contains = params
                    .iter()
                    .any(|&t| t >= span[0] && (t < span[1] || (span[1] == 1.0 && t == 1.0)));
                assert!(
                    contains,
                    "span [{}, {}) holds no parameter for ({n_pts}, {n_ctrl}, {degree})",
                    span[0], span[1]
                );
            }
        }
    }
    /// Known textbook value: 11 uniform parameters, 6 control points, cubic.
    #[test]
    fn test_select_knots_least_squares_matches_piegl_tiller_formula() {
        let params: Vec<f64> = (0..11).map(|i| i as f64 / 10.0).collect();
        let kv = BsplineFitting::select_knots_by_averaging(&params, 6, 3);
        // d = 11 / 3; j = 1: i = 3, alpha = 2/3; j = 2: i = 7, alpha = 1/3.
        let u4 = (1.0 / 3.0) * params[2] + (2.0 / 3.0) * params[3];
        let u5 = (2.0 / 3.0) * params[6] + (1.0 / 3.0) * params[7];
        assert!((kv.knots[4] - u4).abs() < 1e-14, "{:?}", kv.knots);
        assert!((kv.knots[5] - u5).abs() < 1e-14, "{:?}", kv.knots);
    }
    #[test]
    fn test_select_knots_inconsistent_arguments_do_not_panic() {
        // Fewer parameters than control points, no parameters, too few
        // control points for the degree.
        let kv = BsplineFitting::select_knots_by_averaging(&[0.0, 1.0], 6, 3);
        assert_eq!(kv.knots.len(), 10);
        let kv = BsplineFitting::select_knots_by_averaging(&[], 6, 3);
        assert_eq!(kv.knots.len(), 10);
        let kv = BsplineFitting::select_knots_by_averaging(&[0.0, 0.5, 1.0], 2, 3);
        assert_eq!(kv.knots.len(), 6);
    }
    #[test]
    fn test_select_knots_non_decreasing() {
        let params = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
        let kv = BsplineFitting::select_knots_by_averaging(&params, 4, 2);
        for i in 1..kv.knots.len() {
            assert!(
                kv.knots[i] >= kv.knots[i - 1],
                "knots not non-decreasing at {i}: {} < {}",
                kv.knots[i],
                kv.knots[i - 1]
            );
        }
    }
    #[test]
    fn test_bspline_fitting_straight_line() {
        let pts: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts).unwrap();
        assert!(fitter.curve.is_some());
        // Collinear, uniformly spaced data is reproduced exactly by a cubic
        // that has the straight line in its span.
        let residual = fitter.residual(&pts);
        assert!(residual < 1e-18, "straight-line residual = {residual}");
    }
    #[test]
    fn test_bspline_fitting_residual_finite() {
        let pts: Vec<[f64; 3]> = (0..6)
            .map(|i| [i as f64 * 0.5, (i as f64 * 0.5).sin(), 0.0])
            .collect();
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts).unwrap();
        let res = fitter.residual(&pts);
        assert!(res.is_finite());
        // A single cubic follows half a sine period to within a few hundredths.
        let rms = (res / pts.len() as f64).sqrt();
        assert!(rms < 0.03, "sine least-squares RMS error = {rms}");
        // The fit is the minimiser of a convex quadratic in the control
        // points, so moving any coordinate of any control point in either
        // direction cannot reduce the residual.
        let curve = fitter.curve.as_ref().unwrap();
        for cp in 0..curve.control_points.len() {
            for d in 0..3 {
                for delta in [-1e-3, 1e-3] {
                    let mut moved = curve.control_points.clone();
                    moved[cp][d] += delta;
                    let probe = BsplineFitting {
                        degree: 3,
                        n_ctrl: 4,
                        curve: Some(BsplineCurve::new(3, curve.basis.knot_vector.clone(), moved)),
                    };
                    assert!(
                        probe.residual(&pts) >= res - 1e-12,
                        "moving control point {cp} coordinate {d} by {delta} lowered the residual"
                    );
                }
            }
        }
    }
    #[test]
    fn test_bspline_fitting_curve_endpoints_approx() {
        let pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts).unwrap();
        let curve = fitter.curve.as_ref().unwrap();
        // n_pts == n_ctrl is an interpolation: the clamped curve starts and
        // ends exactly at the first and last data point.
        let p0 = curve.eval(0.0);
        let p1 = curve.eval(1.0);
        for d in 0..3 {
            assert!((p0[d] - pts[0][d]).abs() < 1e-9, "start[{d}] = {}", p0[d]);
            assert!((p1[d] - pts[3][d]).abs() < 1e-9, "end[{d}] = {}", p1[d]);
        }
    }
    /// Sine samples with `n_pts == n_ctrl` must be interpolated: the curve
    /// evaluated at each chord-length parameter returns the data point.  The
    /// previous implementation discarded the collocation matrix and used the
    /// data points themselves as control points, so this failed.
    #[test]
    fn test_bspline_fitting_interpolates_sine_samples() {
        let pts: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let x = i as f64 * 0.7;
                [x, x.sin(), 0.25 * x]
            })
            .collect();
        let mut fitter = BsplineFitting::new(3, 6);
        fitter.fit(&pts).unwrap();
        let curve = fitter.curve.as_ref().unwrap();
        assert_eq!(curve.control_points.len(), 6);
        let params = BsplineFitting::chord_length_parameterization(&pts);
        for (i, (&t, q)) in params.iter().zip(pts.iter()).enumerate() {
            let c = curve.eval(t);
            for d in 0..3 {
                assert!(
                    (c[d] - q[d]).abs() < 1e-9,
                    "point {i} coord {d}: curve {} data {}",
                    c[d],
                    q[d]
                );
            }
        }
        // Interpolation is not the identity map on control points: the
        // interior control points differ from the data points.
        let moved = curve
            .control_points
            .iter()
            .zip(pts.iter())
            .any(|(cp, q)| (0..3).any(|d| (cp[d] - q[d]).abs() > 1e-3));
        assert!(moved, "control points must not simply copy the data points");
        assert!((fitter.residual(&pts)).abs() < 1e-18);
    }
    /// Interpolation for several degrees and sizes, including `n_ctrl >= 8`
    /// (the old fixed `[f64; 7]` rows could not even hold `n_ctrl = 5`).
    #[test]
    fn test_bspline_fitting_interpolation_various_degrees_and_sizes() {
        for &(degree, n) in &[(1usize, 4usize), (2, 5), (3, 9), (3, 16), (4, 12), (5, 20)] {
            let pts: Vec<[f64; 3]> = (0..n)
                .map(|i| {
                    let u = i as f64 / (n - 1) as f64;
                    [(3.0 * u).cos(), (2.0 * u).sin() + 0.1 * u, u * u - 0.5 * u]
                })
                .collect();
            let mut fitter = BsplineFitting::new(degree, n);
            fitter.fit(&pts).unwrap();
            let curve = fitter.curve.as_ref().unwrap();
            let params = BsplineFitting::chord_length_parameterization(&pts);
            for (i, (&t, q)) in params.iter().zip(pts.iter()).enumerate() {
                let c = curve.eval(t);
                for d in 0..3 {
                    assert!(
                        (c[d] - q[d]).abs() < 1e-8,
                        "degree {degree} n {n} point {i} coord {d}: {} vs {}",
                        c[d],
                        q[d]
                    );
                }
            }
        }
    }
    /// Least squares with `n_ctrl = 10` (previously a guaranteed index panic
    /// in the solver) on 50 helix samples.
    #[test]
    fn test_bspline_fitting_least_squares_helix_ten_control_points() {
        let pts: Vec<[f64; 3]> = (0..50)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / 49.0;
                [a.cos(), a.sin(), a / (2.0 * PI)]
            })
            .collect();
        let mut fitter = BsplineFitting::new(3, 10);
        fitter.fit(&pts).unwrap();
        let curve = fitter.curve.as_ref().unwrap();
        assert_eq!(curve.control_points.len(), 10);
        let residual = fitter.residual(&pts);
        assert!(residual < 1e-3, "helix residual = {residual}");
        // The fitted curve stays close to the helix at the ends.
        let start = curve.eval(0.0);
        let end = curve.eval(1.0);
        for d in 0..3 {
            assert!((start[d] - pts[0][d]).abs() < 5e-2, "start[{d}]");
            assert!((end[d] - pts[49][d]).abs() < 5e-2, "end[{d}]");
        }
    }
    #[test]
    fn test_bspline_fitting_rejects_invalid_input() {
        let pts: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, (i * i) as f64, 0.0]).collect();
        // Empty and single-point input used to panic or fabricate a line.
        let mut fitter = BsplineFitting::new(3, 4);
        assert!(matches!(
            fitter.fit(&[]),
            Err(crate::Error::TooFewPoints { actual: 0, .. })
        ));
        assert!(fitter.curve.is_none());
        assert!(matches!(
            fitter.fit(&pts[..1]),
            Err(crate::Error::TooFewPoints { actual: 1, .. })
        ));
        // Fewer control points than degree + 1.
        let mut low = BsplineFitting::new(3, 3);
        assert!(matches!(
            low.fit(&pts),
            Err(crate::Error::OutOfRange { name: "n_ctrl", .. })
        ));
        // More control points than data points: the normal equations are
        // singular.
        let mut high = BsplineFitting::new(3, 9);
        assert!(matches!(
            high.fit(&pts),
            Err(crate::Error::OutOfRange { name: "n_ctrl", .. })
        ));
        // Degree 0 has no knot averaging and is not a curve fit.
        let mut flat = BsplineFitting::new(0, 4);
        assert!(matches!(
            flat.fit(&pts),
            Err(crate::Error::OutOfRange { name: "degree", .. })
        ));
        // Non-finite data.
        let mut bad = pts.clone();
        bad[2][1] = f64::NAN;
        let mut nan_fitter = BsplineFitting::new(3, 4);
        assert!(nan_fitter.fit(&bad).is_err());
        assert!(nan_fitter.curve.is_none());
    }
    /// Coincident consecutive data points give equal chord-length parameters
    /// and a singular collocation matrix; the fit must report it instead of
    /// returning a curve built from zeros.
    #[test]
    fn test_bspline_fitting_interpolation_singular_reports_error() {
        let pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let mut fitter = BsplineFitting::new(3, 5);
        let err = fitter.fit(&pts).unwrap_err();
        assert!(err.to_string().contains("singular"), "{err}");
        assert!(fitter.curve.is_none());
    }
    #[test]
    fn test_solve_3x_systems_recovers_known_solution_beyond_five_unknowns() {
        use super::solve_3x_systems;
        let n = 9usize;
        let a: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let off = 1.0 / (1.0 + (i as f64 - j as f64).abs());
                        if i == j { off + n as f64 } else { off }
                    })
                    .collect()
            })
            .collect();
        let x_true: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let f = i as f64;
                [f - 3.0, (0.5 * f).sin(), 2.0 - 0.25 * f * f]
            })
            .collect();
        let rhs: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let mut b = [0.0_f64; 3];
                for (j, x_row) in x_true.iter().enumerate() {
                    for d in 0..3 {
                        b[d] += a[i][j] * x_row[d];
                    }
                }
                b
            })
            .collect();
        let x = solve_3x_systems(&a, &rhs, n).unwrap();
        for i in 0..n {
            for d in 0..3 {
                assert!(
                    (x[i][d] - x_true[i][d]).abs() < 1e-10,
                    "x[{i}][{d}] = {} expected {}",
                    x[i][d],
                    x_true[i][d]
                );
            }
        }
    }
    #[test]
    fn test_solve_3x_systems_needs_pivoting() {
        use super::solve_3x_systems;
        // Zero leading entry: only correct with row pivoting.
        let a = vec![vec![0.0, 1.0], vec![2.0, 0.0]];
        let rhs = vec![[1.0, 2.0, 3.0], [4.0, 6.0, 8.0]];
        let x = solve_3x_systems(&a, &rhs, 2).unwrap();
        assert_eq!(x[0], [2.0, 3.0, 4.0]);
        assert_eq!(x[1], [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_solve_3x_systems_reports_singular_and_malformed_input() {
        use super::solve_3x_systems;
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let rhs = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        let err = solve_3x_systems(&singular, &rhs, 2).unwrap_err();
        assert!(err.to_string().contains("singular"), "{err}");
        let zero = vec![vec![0.0; 3]; 3];
        assert!(solve_3x_systems(&zero, &[[0.0; 3]; 3], 3).is_err());
        // Too few rows / right-hand sides.
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        assert!(solve_3x_systems(&a, &[[1.0, 1.0, 1.0]], 2).is_err());
        assert!(solve_3x_systems(&a[..1], &rhs, 2).is_err());
        assert!(solve_3x_systems(&a, &rhs, 0).unwrap().is_empty());
    }
    #[test]
    fn test_bspline_curve_torsion_planar_zero() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let tau = curve.torsion(0.5);
        assert!(tau.is_finite());
    }
    #[test]
    fn test_bspline_curve_principal_normal_finite() {
        let cps = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
        ];
        let curve = BsplineCurve::clamped(3, cps);
        let n = curve.principal_normal(0.5);
        assert!(n[0].is_finite() && n[1].is_finite() && n[2].is_finite());
    }
    #[test]
    fn test_bspline_surface_gaussian_curvature_flat() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let k = surf.gaussian_curvature(0.5, 0.5);
        assert!(k.abs() < 1e-6, "K of flat plane = {k}");
    }
    #[test]
    fn test_bspline_surface_mean_curvature_flat() {
        let p00 = [0.0, 0.0, 0.0];
        let p10 = [1.0, 0.0, 0.0];
        let p01 = [0.0, 1.0, 0.0];
        let p11 = [1.0, 1.0, 0.0];
        let surf = BsplineSurface::bilinear_patch(p00, p10, p01, p11);
        let h = surf.mean_curvature(0.5, 0.5);
        assert!(h.abs() < 1e-6, "H of flat plane = {h}");
    }
    #[test]
    fn test_nurbs_surface_fit_to_grid() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let mut surf = NurbsSurface::new(1, 1, kv.clone(), kv, cps, weights);
        let target = vec![
            vec![[0.1, 0.1, 0.0], [0.1, 0.9, 0.0]],
            vec![[0.9, 0.1, 0.0], [0.9, 0.9, 0.0]],
        ];
        surf.fit_to_grid(&target);
        let pt = surf.eval(0.0, 0.0);
        assert!(pt[0].is_finite());
    }
    #[test]
    fn test_bspline_basis_degree_zero() {
        let kv = KnotVector::new(vec![0.0, 0.25, 0.5, 0.75, 1.0]);
        let basis = BsplineBasis::new(0, kv, 4);
        let vals = basis.eval_all(0.3);
        let sum: f64 = vals.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "degree-0 PoU sum = {sum}");
    }
}
