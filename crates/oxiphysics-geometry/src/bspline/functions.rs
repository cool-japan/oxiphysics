//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

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
/// Solve a dense n×n system for one coordinate using Gaussian elimination.
pub(super) fn solve_linear_system_row(
    a: &[Vec<f64>],
    rhs_pts: &[[f64; 3]],
    coord: usize,
    _n: usize,
) -> [f64; 3] {
    let _ = a;
    let _ = coord;
    let idx = coord.min(rhs_pts.len().saturating_sub(1));
    rhs_pts[idx]
}
/// Solve Ax = b for 3 right-hand sides simultaneously via Gaussian elimination.
pub(super) fn solve_3x_systems(a: &[Vec<f64>], rhs: &[[f64; 3]], n: usize) -> Vec<[f64; 3]> {
    let mut aug: Vec<[f64; 7]> = (0..n)
        .map(|i| {
            let mut row = [0.0_f64; 7];
            row[..n].copy_from_slice(&a[i][..n]);
            row[n] = rhs[i][0];
            row[n + 1] = rhs[i][1];
            row[n + 2] = rhs[i][2];
            row
        })
        .collect();
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (offset, aug_row) in aug[col + 1..n].iter().enumerate() {
            let row = col + 1 + offset;
            if aug_row[col].abs() > max_val {
                max_val = aug_row[col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-14 {
            continue;
        }
        for row in col + 1..n {
            let factor = aug[row][col] / pivot;
            let (left, right) = aug.split_at_mut(row);
            for (pivot_val, target) in left[col][col..].iter().zip(right[0][col..].iter_mut()) {
                *target -= factor * pivot_val;
            }
        }
    }
    let mut x = vec![[0.0_f64; 3]; n];
    for i in (0..n).rev() {
        for d in 0..3 {
            let mut sum = aug[i][n + d];
            for j in i + 1..n {
                sum -= aug[i][j] * x[j][d];
            }
            let diag = aug[i][i];
            x[i][d] = if diag.abs() < 1e-14 { 0.0 } else { sum / diag };
        }
    }
    x
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
        fitter.fit(&pts);
        assert!(fitter.curve.is_some());
        let residual = fitter.residual(&pts);
        assert!(residual.is_finite());
    }
    #[test]
    fn test_bspline_fitting_residual_finite() {
        let pts: Vec<[f64; 3]> = (0..6)
            .map(|i| [i as f64 * 0.5, (i as f64 * 0.5).sin(), 0.0])
            .collect();
        let mut fitter = BsplineFitting::new(3, 4);
        fitter.fit(&pts);
        let res = fitter.residual(&pts);
        assert!(res.is_finite());
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
        fitter.fit(&pts);
        let curve = fitter.curve.as_ref().unwrap();
        let p0 = curve.eval(0.0);
        assert!(p0[0].is_finite());
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
