//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::BezierCurve;

#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    add3(scale3(a, 1.0 - t), scale3(b, t))
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let len = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(a, 1.0 / len)
    }
}
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Returns a unit vector perpendicular to `v`.
pub(super) fn arbitrary_perp(v: [f64; 3]) -> [f64; 3] {
    let candidates: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for &c in &candidates {
        let proj = dot3(v, c);
        let perp = sub3(c, scale3(v, proj));
        let len = dot3(perp, perp).sqrt();
        if len > 1e-8 {
            return scale3(perp, 1.0 / len);
        }
    }
    [1.0, 0.0, 0.0]
}
/// Approximate a circular arc as a sequence of quadratic Bezier segments.
///
/// Returns the control points for `n_segments` quadratic Bezier arcs,
/// laid out as `[P0, P1, P2, P1', P2', ...]` (each segment shares its last point).
pub fn bezier_arc_approximation(
    center: [f64; 3],
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    n_segments: usize,
) -> Vec<[f64; 3]> {
    if n_segments == 0 {
        return Vec::new();
    }
    let total_angle = end_angle - start_angle;
    let seg_angle = total_angle / n_segments as f64;
    let mut pts = Vec::with_capacity(n_segments * 2 + 1);
    for k in 0..n_segments {
        let a0 = start_angle + k as f64 * seg_angle;
        let a1 = a0 + seg_angle;
        let a_mid = (a0 + a1) / 2.0;
        let p0 = [
            center[0] + radius * a0.cos(),
            center[1] + radius * a0.sin(),
            center[2],
        ];
        let r_mid = radius / a_mid.cos().max(1e-12) * (seg_angle / 2.0).cos();
        let _ = r_mid;
        let t = (seg_angle / 2.0).tan();
        let p1 = [
            center[0] + radius * a_mid.cos() - radius * t * a_mid.sin(),
            center[1] + radius * a_mid.sin() + radius * t * a_mid.cos(),
            center[2],
        ];
        if k == 0 {
            pts.push(p0);
        }
        pts.push(p1);
        let p2 = [
            center[0] + radius * a1.cos(),
            center[1] + radius * a1.sin(),
            center[2],
        ];
        pts.push(p2);
    }
    pts
}
/// Sample a curve uniformly in parameter space, returning `n` points.
///
/// `eval` maps parameter `t ∈ [0, 1]` to a 3-D point.
/// Returns `n` points at parameters `0, 1/(n-1), ..., 1` (arc-length
/// parameterization is approximated by uniform parameter sampling).
pub fn sample_curve_uniform(eval: &dyn Fn(f64) -> [f64; 3], n: usize) -> Vec<[f64; 3]> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![eval(0.0)];
    }
    (0..n).map(|i| eval(i as f64 / (n - 1) as f64)).collect()
}
/// Approximate the curvature κ of a `BezierCurve` at parameter `t`.
///
/// Uses the formula κ = |C'(t) × C''(t)| / |C'(t)|³.
pub fn bezier_curvature(curve: &BezierCurve, t: f64) -> f64 {
    let eps = 1e-5_f64;
    let d1 = curve.derivative(t);
    let t0 = (t - eps).max(0.0);
    let t1 = (t + eps).min(1.0);
    let d1_0 = curve.derivative(t0);
    let d1_1 = curve.derivative(t1);
    let d2 = scale3(sub3(d1_1, d1_0), 1.0 / (t1 - t0));
    let cross_mag = {
        let c = cross3(d1, d2);
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
    };
    let d1_mag = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
    if d1_mag < 1e-12 {
        return 0.0;
    }
    cross_mag / (d1_mag * d1_mag * d1_mag)
}
/// Approximate the torsion τ of a `BezierCurve` at parameter `t`.
///
/// Uses the scalar triple product of C', C'', C''' divided by |C' × C''|².
pub fn bezier_torsion(curve: &BezierCurve, t: f64) -> f64 {
    let eps = 1e-5_f64;
    let d1 = curve.derivative(t);
    let t0 = (t - eps).max(0.0);
    let t1 = (t + eps).min(1.0);
    let d1_0 = curve.derivative(t0);
    let d1_1 = curve.derivative(t1);
    let d2 = scale3(sub3(d1_1, d1_0), 1.0 / (t1 - t0));
    let t00 = (t - 2.0 * eps).max(0.0);
    let t11 = (t + 2.0 * eps).min(1.0);
    let d1_00 = curve.derivative(t00);
    let d1_11 = curve.derivative(t11);
    let d2_0 = scale3(sub3(d1, d1_00), 1.0 / eps);
    let d2_1 = scale3(sub3(d1_11, d1), 1.0 / eps);
    let d3 = scale3(sub3(d2_1, d2_0), 1.0 / (2.0 * eps));
    let cross = cross3(d1, d2);
    let cross_mag2 = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
    if cross_mag2 < 1e-20 {
        return 0.0;
    }
    let dot_val = cross[0] * d3[0] + cross[1] * d3[1] + cross[2] * d3[2];
    dot_val / cross_mag2
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bspline::BsplineSurface;
    use crate::bspline::KnotVector;
    use crate::bspline::NurbsCurve;
    use crate::bspline::NurbsSurface;
    use crate::parametric::Arc;
    use crate::parametric::BSpline;
    use crate::parametric::BezierSurface;
    use crate::parametric::CatmullRomSpline;
    use crate::parametric::FrenetFrame;
    use crate::parametric::HermiteCurve;
    use crate::parametric::LoftSurface;
    use crate::parametric::RevolutionSurface;
    use crate::parametric::SweptSurface;
    use crate::parametric::TangentFrame;
    use crate::parametric::TubeGeometry;
    #[test]
    fn bezier_degree1_endpoints() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [4.0, 2.0, 0.0];
        let curve = BezierCurve::new(vec![p0, p1]);
        let at0 = curve.evaluate(0.0);
        let at1 = curve.evaluate(1.0);
        let at_half = curve.evaluate(0.5);
        for k in 0..3 {
            assert!((at0[k] - p0[k]).abs() < 1e-10);
            assert!((at1[k] - p1[k]).abs() < 1e-10);
            assert!((at_half[k] - (p0[k] + p1[k]) / 2.0).abs() < 1e-10);
        }
    }
    #[test]
    fn bezier_degree2_midpoint() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 2.0, 0.0];
        let p2 = [2.0, 0.0, 0.0];
        let curve = BezierCurve::new(vec![p0, p1, p2]);
        let expected = [1.0, 1.0, 0.0];
        let result = curve.evaluate(0.5);
        for k in 0..3 {
            assert!(
                (result[k] - expected[k]).abs() < 1e-10,
                "component {k}: got {}",
                result[k]
            );
        }
    }
    #[test]
    fn bezier_arc_length_line() {
        let curve = BezierCurve::new(vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let len = curve.arc_length(1000);
        assert!((len - 3.0).abs() < 1e-6, "expected ≈3.0, got {len}");
    }
    #[test]
    fn bezier_surface_flat_normal() {
        let grid: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|i| {
                (0..4)
                    .map(|j| [i as f64 / 3.0, j as f64 / 3.0, 0.0])
                    .collect()
            })
            .collect();
        let surf = BezierSurface::new_bicubic(grid);
        let n = surf.normal(0.5, 0.5);
        assert!(n[0].abs() < 1e-4, "nx={}", n[0]);
        assert!(n[1].abs() < 1e-4, "ny={}", n[1]);
        assert!((n[2].abs() - 1.0).abs() < 1e-4, "nz={}", n[2]);
    }
    #[test]
    fn bspline_clamped_endpoints() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 0.0, 1.0],
            [4.0, 1.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(3, pts.clone());
        let at0 = spline.eval(0.0);
        let at1 = spline.eval(1.0);
        for k in 0..3 {
            assert!(
                (at0[k] - pts[0][k]).abs() < 1e-10,
                "start: axis {k}: {}",
                at0[k]
            );
            assert!(
                (at1[k] - pts[pts.len() - 1][k]).abs() < 1e-10,
                "end: axis {k}: {}",
                at1[k]
            );
        }
    }
    #[test]
    fn bspline_knot_insert_preserves_curve() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(2, pts);
        let new_spline = spline.insert_knot(0.5);
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let p_orig = spline.eval(t);
            let p_new = new_spline.eval(t);
            for k in 0..3 {
                assert!(
                    (p_orig[k] - p_new[k]).abs() < 1e-9,
                    "t={t} axis {k}: orig={} new={}",
                    p_orig[k],
                    p_new[k]
                );
            }
        }
    }
    #[test]
    fn bspline_degree1_is_polyline() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(1, pts.clone());
        let p0 = spline.eval(0.0);
        let p1 = spline.eval(1.0);
        for k in 0..3 {
            assert!((p0[k] - pts[0][k]).abs() < 1e-9);
            assert!((p1[k] - pts[pts.len() - 1][k]).abs() < 1e-9);
        }
    }
    #[test]
    fn hermite_interpolates_endpoints() {
        let pts = vec![[0.0, 0.0, 0.0], [2.0, 3.0, 0.0], [4.0, 0.0, 0.0]];
        let tgts = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]];
        let curve = HermiteCurve::new(pts.clone(), tgts);
        for (i, &pt) in pts.iter().enumerate() {
            let evaled = curve.eval(i as f64);
            for k in 0..3 {
                assert!(
                    (evaled[k] - pt[k]).abs() < 1e-10,
                    "seg {i} axis {k}: {}",
                    evaled[k]
                );
            }
        }
    }
    #[test]
    fn hermite_c0_continuity() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let tgts = vec![
            [0.5, 0.5, 0.0],
            [0.5, -0.5, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, -0.5, 0.0],
        ];
        let curve = HermiteCurve::new(pts, tgts);
        for knot in [1.0, 2.0] {
            let eps = 1e-8;
            let left = curve.eval(knot - eps);
            let right = curve.eval(knot + eps);
            for k in 0..3 {
                assert!((left[k] - right[k]).abs() < 1e-5, "C0 at {knot} axis {k}");
            }
        }
    }
    #[test]
    fn hermite_derivative_at_knot() {
        let pts = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tgts = vec![[3.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let curve = HermiteCurve::new(pts, tgts.clone());
        let deriv = curve.derivative(0.0);
        for k in 0..3 {
            assert!(
                (deriv[k] - tgts[0][k]).abs() < 1e-10,
                "axis {k}: {}",
                deriv[k]
            );
        }
    }
    #[test]
    fn catmull_rom_interpolates_interior() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let spline = CatmullRomSpline::new(pts.clone());
        for (i, &pt) in pts.iter().enumerate() {
            let t = spline.knots[i];
            let evaled = spline.eval(t);
            for k in 0..3 {
                assert!(
                    (evaled[k] - pt[k]).abs() < 1e-6,
                    "point {i} axis {k}: got {} expected {}",
                    evaled[k],
                    pt[k]
                );
            }
        }
    }
    #[test]
    fn catmull_rom_range_clamp() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let spline = CatmullRomSpline::new(pts);
        let p_neg = spline.eval(-0.1);
        let p_over = spline.eval(1.5);
        for k in 0..3 {
            assert!(p_neg[k].is_finite());
            assert!(p_over[k].is_finite());
        }
    }
    #[test]
    fn frenet_tangent_is_unit() {
        let circle_pts: Vec<[f64; 3]> = (0..=20)
            .map(|i| {
                let t = i as f64 / 20.0 * std::f64::consts::TAU;
                [t.cos(), t.sin(), 0.0]
            })
            .collect();
        let curve = BezierCurve::new(circle_pts.clone());
        let frame = FrenetFrame::from_bezier(&curve, 0.5);
        let len =
            (frame.tangent[0].powi(2) + frame.tangent[1].powi(2) + frame.tangent[2].powi(2)).sqrt();
        assert!((len - 1.0).abs() < 1e-6, "tangent not unit: len={len}");
    }
    #[test]
    fn frenet_triad_orthonormal() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BezierCurve::new(pts);
        let frame = FrenetFrame::from_bezier(&curve, 0.5);
        let t = frame.tangent;
        let n = frame.normal;
        let b = frame.binormal;
        assert!(dot3(t, n).abs() < 1e-6, "T·N = {}", dot3(t, n));
        assert!(dot3(t, b).abs() < 1e-6, "T·B = {}", dot3(t, b));
        assert!(dot3(n, b).abs() < 1e-6, "N·B = {}", dot3(n, b));
        for (name, v) in [("T", t), ("N", n), ("B", b)] {
            let len = (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
            assert!((len - 1.0).abs() < 1e-6, "{name} not unit: {len}");
        }
    }
    #[test]
    fn tube_vertex_count() {
        let curve = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let tube = TubeGeometry::from_bezier(&curve, 0.1, 8, 12);
        assert_eq!(tube.vertices.len(), 8 * 12, "vertex count");
        assert_eq!(tube.triangles.len(), (8 - 1) * 12 * 2, "triangle count");
    }
    #[test]
    fn tube_vertices_near_spine() {
        let pts = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let curve = BezierCurve::new(pts);
        let radius = 0.5;
        let tube = TubeGeometry::from_bezier(&curve, radius, 10, 8);
        for v in &tube.vertices {
            let r = (v[1].powi(2) + v[2].powi(2)).sqrt();
            assert!(
                (r - radius).abs() < 1e-8,
                "vertex [{},{},{}] radius {r} != {radius}",
                v[0],
                v[1],
                v[2]
            );
        }
    }
    #[test]
    fn bezier_derivative_line() {
        let curve = BezierCurve::new(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        for &t in &[0.0, 0.5, 1.0] {
            let d = curve.derivative(t);
            assert!((d[0] - 2.0).abs() < 1e-10, "dx at t={t}: {}", d[0]);
            assert!(d[1].abs() < 1e-10);
            assert!(d[2].abs() < 1e-10);
        }
    }
    #[test]
    fn bezier_eval_at_t0_returns_first_control_point() {
        let p0 = [1.0, 2.0, 3.0];
        let p1 = [4.0, 5.0, 6.0];
        let p2 = [7.0, 8.0, 9.0];
        let curve = BezierCurve::new(vec![p0, p1, p2]);
        let result = curve.eval(0.0);
        for k in 0..3 {
            assert!(
                (result[k] - p0[k]).abs() < 1e-10,
                "axis {k}: {} != {}",
                result[k],
                p0[k]
            );
        }
    }
    #[test]
    fn bezier_eval_at_t1_returns_last_control_point() {
        let p0 = [1.0, 2.0, 3.0];
        let p1 = [4.0, 5.0, 6.0];
        let p2 = [7.0, 8.0, 9.0];
        let curve = BezierCurve::new(vec![p0, p1, p2]);
        let result = curve.eval(1.0);
        for k in 0..3 {
            assert!(
                (result[k] - p2[k]).abs() < 1e-10,
                "axis {k}: {} != {}",
                result[k],
                p2[k]
            );
        }
    }
    #[test]
    fn bspline_order_clamped() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 0.0, 1.0],
            [4.0, 1.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(3, pts.clone());
        let start = spline.eval(0.0);
        let end = spline.eval(1.0);
        for k in 0..3 {
            assert!(
                (start[k] - pts[0][k]).abs() < 1e-9,
                "start axis {k}: got {} expected {}",
                start[k],
                pts[0][k]
            );
            assert!(
                (end[k] - pts[pts.len() - 1][k]).abs() < 1e-9,
                "end axis {k}: got {} expected {}",
                end[k],
                pts[pts.len() - 1][k]
            );
        }
    }
    #[test]
    fn nurbs_surface_eval_non_nan() {
        let mut control_points: Vec<Vec<[f64; 3]>> = Vec::new();
        let mut weights: Vec<Vec<f64>> = Vec::new();
        for i in 0..3 {
            let mut cp_row = Vec::new();
            let mut w_row = Vec::new();
            for j in 0..3 {
                cp_row.push([i as f64, j as f64, 0.0]);
                w_row.push(1.0);
            }
            control_points.push(cp_row);
            weights.push(w_row);
        }
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let surf = NurbsSurface::new(
            2,
            2,
            KnotVector::new(knots.clone()),
            KnotVector::new(knots),
            control_points,
            weights,
        );
        for &u in &[0.0, 0.5, 1.0] {
            for &v in &[0.0, 0.5, 1.0] {
                let p = surf.eval(u, v);
                for &val in &p {
                    assert!(val.is_finite(), "NurbsSurface eval({u},{v}) returned NaN");
                }
            }
        }
    }
    #[test]
    fn bezier_arc_approximation_points_on_circle() {
        let center = [0.0, 0.0, 0.0];
        let radius = 2.0;
        let pts = bezier_arc_approximation(center, radius, 0.0, std::f64::consts::PI, 4);
        assert!(!pts.is_empty(), "arc approximation should return points");
        let first = pts[0];
        let last = *pts.last().unwrap();
        let r0 = (first[0].powi(2) + first[1].powi(2)).sqrt();
        let r1 = (last[0].powi(2) + last[1].powi(2)).sqrt();
        assert!(
            (r0 - radius).abs() < 1e-9,
            "first point not on circle: r={r0}"
        );
        assert!(
            (r1 - radius).abs() < 1e-9,
            "last point not on circle: r={r1}"
        );
    }
    #[test]
    fn sample_curve_uniform_count() {
        let n = 10usize;
        let pts = sample_curve_uniform(&|t| [t, t * 2.0, 0.0], n);
        assert_eq!(pts.len(), n);
        assert!((pts[0][0] - 0.0).abs() < 1e-10);
        assert!((pts[n - 1][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn nurbs_basis_partition_of_unity() {
        let degree = 3_usize;
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        let n_ctrl = 5_usize;
        for &t in &[0.0, 0.25, 0.5, 0.75] {
            let sum: f64 = (0..n_ctrl)
                .map(|i| NurbsCurve::b_spline_basis(i, degree, t, &knots))
                .sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "partition of unity failed at t={t}: sum={sum}"
            );
        }
        let sum_at_1: f64 = (0..n_ctrl)
            .map(|i| NurbsCurve::b_spline_basis(i, degree, 1.0, &knots))
            .sum();
        assert!(
            (sum_at_1 - 1.0).abs() < 1e-10,
            "partition of unity failed at t=1.0: sum={sum_at_1}"
        );
    }
    #[test]
    fn arc_endpoints_on_circle() {
        let center = [1.0, 2.0, 0.0];
        let radius = 3.0;
        let arc = Arc::in_xy_plane(center, radius, 0.0, std::f64::consts::PI);
        let p0 = arc.evaluate(0.0);
        let p1 = arc.evaluate(1.0);
        let r0 = dist3(p0, center);
        let r1 = dist3(p1, center);
        assert!((r0 - radius).abs() < 1e-9, "start not on circle: r={r0}");
        assert!((r1 - radius).abs() < 1e-9, "end not on circle: r={r1}");
    }
    #[test]
    fn arc_length_half_circle() {
        let radius = 2.0;
        let arc = Arc::in_xy_plane([0.0; 3], radius, 0.0, std::f64::consts::PI);
        let expected = std::f64::consts::PI * radius;
        let len = arc.arc_length();
        assert!(
            (len - expected).abs() < 1e-9,
            "expected {expected}, got {len}"
        );
    }
    #[test]
    fn arc_sample_count() {
        let arc = Arc::in_xy_plane([0.0; 3], 1.0, 0.0, std::f64::consts::TAU);
        let pts = arc.sample(10);
        assert_eq!(pts.len(), 10);
    }
    #[test]
    fn arc_tangent_perpendicular_to_radius() {
        let center = [0.0; 3];
        let arc = Arc::in_xy_plane(center, 1.0, 0.0, std::f64::consts::TAU);
        for i in 0..4 {
            let t = i as f64 / 4.0;
            let p = arc.evaluate(t);
            let tan = arc.tangent(t);
            let radial = sub3(p, center);
            let d = dot3(radial, tan);
            assert!(
                d.abs() < 1e-6,
                "tangent not perpendicular to radius at t={t}: dot={d}"
            );
        }
    }
    #[test]
    fn loft_surface_at_v0_is_curve_a() {
        let ca = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let cb = BezierCurve::new(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]);
        let loft = LoftSurface::new(ca, cb);
        for &u in &[0.0, 0.5, 1.0] {
            let p = loft.evaluate(u, 0.0);
            let expected = loft.curve_a.evaluate(u);
            for k in 0..3 {
                assert!(
                    (p[k] - expected[k]).abs() < 1e-9,
                    "v=0 loft should equal curve_a"
                );
            }
        }
    }
    #[test]
    fn loft_surface_at_v1_is_curve_b() {
        let ca = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let cb = BezierCurve::new(vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]);
        let loft = LoftSurface::new(ca, cb);
        for &u in &[0.0, 0.5, 1.0] {
            let p = loft.evaluate(u, 1.0);
            let expected = loft.curve_b.evaluate(u);
            for k in 0..3 {
                assert!(
                    (p[k] - expected[k]).abs() < 1e-9,
                    "v=1 loft should equal curve_b"
                );
            }
        }
    }
    #[test]
    fn loft_surface_normal_finite() {
        let ca = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let cb = BezierCurve::new(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]);
        let loft = LoftSurface::new(ca, cb);
        let n = loft.normal(0.5, 0.5);
        for &v in &n {
            assert!(v.is_finite(), "loft normal must be finite: {v}");
        }
    }
    #[test]
    fn loft_sample_grid_dimensions() {
        let ca = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let cb = BezierCurve::new(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]);
        let loft = LoftSurface::new(ca, cb);
        let grid = loft.sample_grid(4, 5);
        assert_eq!(grid.len(), 4);
        assert_eq!(grid[0].len(), 5);
    }
    #[test]
    fn swept_surface_at_u0_is_at_spine_start() {
        let spine = BezierCurve::new(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let profile = BezierCurve::new(vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let swept = SweptSurface::new(spine, profile);
        let p00 = swept.evaluate(0.0, 0.0);
        for (k, &p00k) in p00.iter().enumerate() {
            assert!(p00k.abs() < 1e-9, "p00[{k}]={}", p00k);
        }
    }
    #[test]
    fn swept_surface_normal_is_unit() {
        let spine = BezierCurve::new(vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        let profile = BezierCurve::new(vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let swept = SweptSurface::new(spine, profile);
        let n = swept.normal(0.5, 0.5);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-5 || len < 1e-9,
            "normal not unit: len={len}"
        );
    }
    #[test]
    fn revolution_surface_full_circle_radius() {
        let profile = BezierCurve::new(vec![[2.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let rev = RevolutionSurface::new(profile);
        for i in 0..8 {
            let u = i as f64 / 8.0;
            let p = rev.evaluate(u, 0.5);
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!((r - 2.0).abs() < 1e-9, "radius at u={u}: {r}");
        }
    }
    #[test]
    fn revolution_surface_sample_grid_dimensions() {
        let profile = BezierCurve::new(vec![[1.0, 0.0, 0.0], [1.0, 0.0, 2.0]]);
        let rev = RevolutionSurface::new(profile);
        let grid = rev.sample_grid(6, 4);
        assert_eq!(grid.len(), 6);
        assert_eq!(grid[0].len(), 4);
    }
    #[test]
    fn tangent_frame_bezier_surface_normal_unit() {
        let grid: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|i| {
                (0..4)
                    .map(|j| [i as f64 / 3.0, j as f64 / 3.0, 0.0])
                    .collect()
            })
            .collect();
        let surf = BezierSurface::new_bicubic(grid);
        let frame = TangentFrame::from_bezier_surface(&surf, 0.5, 0.5);
        let n_len =
            (frame.normal[0].powi(2) + frame.normal[1].powi(2) + frame.normal[2].powi(2)).sqrt();
        assert!((n_len - 1.0).abs() < 1e-5, "normal not unit: len={n_len}");
    }
    #[test]
    fn tangent_frame_loft_orthogonality() {
        let ca = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let cb = BezierCurve::new(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]);
        let loft = LoftSurface::new(ca, cb);
        let frame = TangentFrame::from_loft_surface(&loft, 0.5, 0.5);
        let d = dot3(frame.tangent_u, frame.tangent_v);
        assert!(d.abs() < 1e-5, "tangent_u . tangent_v = {d} (should be ~0)");
    }
    #[test]
    fn bspline_derivative_finite() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(3, pts);
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let d = spline.derivative(t);
            for (k, &dk) in d.iter().enumerate() {
                assert!(
                    dk.is_finite(),
                    "derivative[{k}] not finite at t={t}: {}",
                    dk
                );
            }
        }
    }
    #[test]
    fn bspline_arc_length_line_approx() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let spline = BSpline::clamped_uniform(1, pts);
        let len = spline.arc_length(100);
        assert!(
            (len - 4.0).abs() < 1e-3,
            "arc length should be ~4.0, got {len}"
        );
    }
    #[test]
    fn bezier_curvature_straight_line_is_zero() {
        let curve = BezierCurve::new(vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]);
        let kappa = bezier_curvature(&curve, 0.5);
        assert!(
            kappa < 1e-6,
            "straight line curvature should be 0, got {kappa}"
        );
    }
    #[test]
    fn bezier_curvature_circle_approx_positive() {
        let circle_pts: Vec<[f64; 3]> = (0..=8)
            .map(|i| {
                let t = i as f64 / 8.0 * std::f64::consts::TAU;
                [t.cos(), t.sin(), 0.0]
            })
            .collect();
        let curve = BezierCurve::new(circle_pts);
        let kappa = bezier_curvature(&curve, 0.5);
        assert!(
            kappa.is_finite() && kappa >= 0.0,
            "curvature should be non-negative: {kappa}"
        );
    }
    #[test]
    fn bezier_torsion_planar_curve_near_zero() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BezierCurve::new(pts);
        let tau = bezier_torsion(&curve, 0.5);
        assert!(tau.is_finite(), "torsion must be finite: {tau}");
        assert!(
            tau.abs() < 1.0,
            "planar curve torsion should be small: {tau}"
        );
    }
    fn flat_bspline_surface() -> BsplineSurface {
        let net: Vec<Vec<[f64; 3]>> = (0..3)
            .map(|i| (0..3).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        BsplineSurface::new(
            2,
            2,
            KnotVector::new(knots.clone()),
            KnotVector::new(knots),
            net,
        )
    }
    #[test]
    fn bspline_surface_eval_finite() {
        let surf = flat_bspline_surface();
        for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
                let p = surf.eval(u, v);
                for (k, &pk) in p.iter().enumerate() {
                    assert!(pk.is_finite(), "eval({u},{v})[{k}] not finite: {}", pk);
                }
            }
        }
    }
    #[test]
    fn bspline_surface_flat_gaussian_curvature_near_zero() {
        let surf = flat_bspline_surface();
        let (k_gauss, _k_mean) = surf.compute_curvature(0.5, 0.5);
        assert!(
            k_gauss.abs() < 1e-3,
            "flat surface Gaussian curvature should be ~0, got {k_gauss}"
        );
    }
    #[test]
    fn bspline_surface_flat_mean_curvature_near_zero() {
        let surf = flat_bspline_surface();
        let (_k_gauss, k_mean) = surf.compute_curvature(0.5, 0.5);
        assert!(
            k_mean.abs() < 1e-3,
            "flat surface mean curvature should be ~0, got {k_mean}"
        );
    }
    #[test]
    fn bspline_surface_curvature_returns_finite() {
        let surf = flat_bspline_surface();
        let (k_g, k_h) = surf.compute_curvature(0.3, 0.7);
        assert!(
            k_g.is_finite(),
            "Gaussian curvature should be finite: {k_g}"
        );
        assert!(k_h.is_finite(), "mean curvature should be finite: {k_h}");
    }
    #[test]
    fn bspline_surface_refine_knots_increases_control_points() {
        let surf = flat_bspline_surface();
        let refined = surf.refine_knots(&[0.5], &[0.5]);
        let n_u_before = surf.control_points.len();
        let n_u_after = refined.control_points.len();
        assert!(
            n_u_after >= n_u_before,
            "refine_knots should not reduce control point count in u: {n_u_after} < {n_u_before}"
        );
    }
    #[test]
    fn bspline_surface_refine_knots_preserves_geometry() {
        let surf = flat_bspline_surface();
        let refined = surf.refine_knots(&[0.5], &[0.5]);
        for &(u, v) in &[(0.25, 0.25), (0.5, 0.5), (0.75, 0.75)] {
            let p_orig = surf.eval(u, v);
            let p_refined = refined.eval(u, v);
            for k in 0..3 {
                assert!(
                    (p_orig[k] - p_refined[k]).abs() < 1e-6,
                    "eval changed after refinement at ({u},{v})[{k}]: {} vs {}",
                    p_orig[k],
                    p_refined[k]
                );
            }
        }
    }
    fn simple_nurbs_line() -> NurbsCurve {
        NurbsCurve::new(
            1,
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]),
            vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        )
    }
    #[test]
    fn nurbs_arc_length_line() {
        let curve = simple_nurbs_line();
        let len = curve.arc_length(0.0, 1.0, 200);
        assert!(
            (len - 4.0).abs() < 1e-3,
            "NURBS line arc length should be ~4.0, got {len}"
        );
    }
    #[test]
    fn nurbs_arc_length_reparametrize_count() {
        let curve = simple_nurbs_line();
        let pts = curve.compute_arc_length_reparametrize(200, 10);
        assert_eq!(
            pts.len(),
            10,
            "should return exactly 10 reparametrized points"
        );
    }
    #[test]
    fn nurbs_arc_length_reparametrize_uniform_spacing() {
        let curve = simple_nurbs_line();
        let pts = curve.compute_arc_length_reparametrize(500, 5);
        assert_eq!(pts.len(), 5);
        for i in 1..pts.len() {
            let d = dist3(pts[i - 1], pts[i]);
            assert!(
                (d - 1.0).abs() < 5e-2,
                "inter-point distance should be ~1.0, got {d} between pts[{}] and pts[{}]",
                i - 1,
                i
            );
        }
    }
    #[test]
    fn nurbs_arc_length_table_monotone() {
        let curve = simple_nurbs_line();
        let (arc_lengths, _params) = curve.arc_length_table(100);
        for w in arc_lengths.windows(2) {
            assert!(
                w[1] >= w[0] - 1e-12,
                "arc length table must be non-decreasing"
            );
        }
    }
}
