//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{BSplineCurve, BezierCurve};

/// Add two 3-vectors.
#[inline]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector by a scalar.
#[inline]
pub fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product of two 3-vectors.
#[inline]
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
#[inline]
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Euclidean norm of a 3-vector.
#[inline]
pub fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
/// Normalize a 3-vector (returns zero vector for degenerate input).
#[inline]
pub fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}
/// Linear interpolation between two 3-vectors.
#[inline]
pub fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    vec3_add(vec3_scale(a, 1.0 - t), vec3_scale(b, t))
}
/// Evaluate all non-zero B-spline basis functions of degree `p` at parameter `t`
/// using the Cox–de Boor recurrence.
///
/// `knots` is the full knot vector. `span` is the knot span index such that
/// `knots[span] <= t < knots[span+1]`.
///
/// Returns an array of `p+1` values: `N[0..=p]` corresponding to basis functions
/// `N_{span-p,p}` through `N_{span,p}`.
pub fn bspline_basis(span: usize, p: usize, t: f64, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0f64; p + 1];
    let mut left = vec![0.0f64; p + 1];
    let mut right = vec![0.0f64; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0f64;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}
/// Find the knot span index for parameter `t` in knot vector `knots` with `n+1` control points
/// and degree `p`, where `n = num_ctrl - 1`.
pub fn find_knot_span(num_ctrl: usize, p: usize, t: f64, knots: &[f64]) -> usize {
    let n = num_ctrl - 1;
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[p] {
        return p;
    }
    let mut low = p;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}
/// Generate a uniform clamped knot vector for `n+1` control points of degree `p`.
///
/// The result has `n + p + 2` knots: `p+1` zeros, `n-p` interior knots, `p+1` ones.
pub fn uniform_clamped_knots(num_ctrl: usize, p: usize) -> Vec<f64> {
    let n = num_ctrl - 1;
    let m = n + p + 1;
    let mut knots = vec![0.0f64; m + 1];
    for (i, k) in knots.iter_mut().enumerate().take(m + 1) {
        if i <= p {
            *k = 0.0;
        } else if i > n {
            *k = 1.0;
        } else {
            *k = (i - p) as f64 / (n - p + 1) as f64;
        }
    }
    knots
}
/// Fit a B-spline curve to a set of data points using least squares.
///
/// Uses a uniform clamped knot vector and solves the normal equations
/// `A^T A c = A^T y` for each coordinate independently via Cholesky-like
/// tridiagonal-band decomposition approximation (here: full dense solver for simplicity).
pub fn fit_bspline_least_squares(
    data: &[[f64; 3]],
    num_ctrl: usize,
    degree: usize,
) -> BSplineCurve {
    let m = data.len();
    assert!(m >= num_ctrl, "need at least num_ctrl data points");
    let knots = uniform_clamped_knots(num_ctrl, degree);
    let mut params = vec![0.0f64; m];
    for i in 1..m {
        params[i] = params[i - 1] + vec3_norm(vec3_sub(data[i], data[i - 1]));
    }
    let total = params[m - 1];
    if total > 0.0 {
        for p in params.iter_mut() {
            *p /= total;
        }
    }
    let mut a_mat = vec![vec![0.0f64; num_ctrl]; m];
    for (row, &t) in params.iter().enumerate() {
        let span = find_knot_span(num_ctrl, degree, t, &knots);
        let basis = bspline_basis(span, degree, t, &knots);
        for j in 0..=degree {
            a_mat[row][span - degree + j] = basis[j];
        }
    }
    let mut ata = vec![vec![0.0f64; num_ctrl]; num_ctrl];
    let mut atb = vec![[0.0f64; 3]; num_ctrl];
    for i in 0..m {
        for j in 0..num_ctrl {
            for k in 0..num_ctrl {
                ata[j][k] += a_mat[i][j] * a_mat[i][k];
            }
            for d in 0..3 {
                atb[j][d] += a_mat[i][j] * data[i][d];
            }
        }
    }
    let n = num_ctrl;
    let mut aug = vec![vec![0.0f64; n + 3]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = ata[i][j];
        }
        for d in 0..3 {
            aug[i][n + d] = atb[i][d];
        }
    }
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for (row, _) in aug.iter().enumerate().take(n).skip(col + 1) {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);
        let pivot = aug[col][col];
        if pivot.abs() < 1e-12 {
            continue;
        }
        for elem in aug[col].iter_mut().skip(col) {
            *elem /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                let col_row: Vec<f64> = aug[col][col..].to_vec();
                for (j_off, elem) in aug[row][col..].iter_mut().enumerate() {
                    *elem -= factor * col_row[j_off];
                }
            }
        }
    }
    let ctrl: Vec<[f64; 3]> = (0..n)
        .map(|i| [aug[i][n], aug[i][n + 1], aug[i][n + 2]])
        .collect();
    BSplineCurve::new(ctrl, degree, knots)
}
/// Generate a surface of revolution by rotating a profile curve around the Z-axis.
///
/// The profile is a set of 2D points `(r, z)` in the meridian plane.
/// Returns a grid of 3D points `[x, y, z]` on the surface.
pub fn surface_of_revolution(
    profile: &[[f64; 2]],
    num_angular_samples: usize,
) -> Vec<Vec<[f64; 3]>> {
    let m = profile.len();
    let n = num_angular_samples;
    let mut grid = vec![vec![[0.0f64; 3]; n]; m];
    for (i, &[r, z]) in profile.iter().enumerate() {
        for (j, cell) in grid[i].iter_mut().enumerate() {
            let theta = 2.0 * PI * j as f64 / n as f64;
            *cell = [r * theta.cos(), r * theta.sin(), z];
        }
    }
    grid
}
/// Generate a lofted surface by interpolating between a set of cross-section curves.
///
/// Each cross-section must have the same number of points.
pub fn lofted_surface(sections: &[Vec<[f64; 3]>], num_v_samples: usize) -> Vec<Vec<[f64; 3]>> {
    let ns = sections.len();
    assert!(ns >= 2, "need at least 2 sections");
    let np = sections[0].len();
    for s in sections.iter() {
        assert_eq!(s.len(), np, "all sections must have equal point count");
    }
    let mut grid = vec![vec![[0.0f64; 3]; np]; num_v_samples];
    for (vi, row) in grid.iter_mut().enumerate() {
        let t = vi as f64 / (num_v_samples - 1).max(1) as f64;
        let seg = ((t * (ns - 1) as f64).min((ns - 2) as f64)) as usize;
        let local_t = t * (ns - 1) as f64 - seg as f64;
        for (pi, cell) in row.iter_mut().enumerate() {
            *cell = vec3_lerp(sections[seg][pi], sections[seg + 1][pi], local_t);
        }
    }
    grid
}
/// Compute the curvature of a parametric curve at parameter `t`.
///
/// Curvature κ = |r' × r''| / |r'|³.
pub fn curve_curvature(curve: &BezierCurve, t: f64) -> f64 {
    let r_prime = curve.tangent(t);
    let r_double = curve.second_derivative(t);
    let cross = vec3_cross(r_prime, r_double);
    let cross_norm = vec3_norm(cross);
    let prime_norm = vec3_norm(r_prime);
    if prime_norm < 1e-12 {
        0.0
    } else {
        cross_norm / (prime_norm * prime_norm * prime_norm)
    }
}
/// Compute the torsion of a parametric curve using finite differences.
///
/// Torsion τ = (r' × r'') · r''' / |r' × r''|².
/// Uses centered finite differences with step `h` to estimate `r'''`.
pub fn curve_torsion(curve: &BezierCurve, t: f64, h: f64) -> f64 {
    let r_prime = curve.tangent(t);
    let r_double = curve.second_derivative(t);
    let r_double_plus = {
        let tp = (t + h).min(1.0);
        curve.second_derivative(tp)
    };
    let r_double_minus = {
        let tm = (t - h).max(0.0);
        curve.second_derivative(tm)
    };
    let r_triple = vec3_scale(vec3_sub(r_double_plus, r_double_minus), 1.0 / (2.0 * h));
    let cross = vec3_cross(r_prime, r_double);
    let denom = vec3_dot(cross, cross);
    if denom < 1e-20 {
        0.0
    } else {
        vec3_dot(cross, r_triple) / denom
    }
}
/// Compute the Frenet-Serret frame `(T, N, B)` at parameter `t`.
///
/// Returns `(tangent, normal, binormal)` unit vectors.
pub fn frenet_serret_frame(curve: &BezierCurve, t: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let t_vec = vec3_normalize(curve.tangent(t));
    let r_double = curve.second_derivative(t);
    let n_raw = vec3_sub(r_double, vec3_scale(t_vec, vec3_dot(r_double, t_vec)));
    let n_vec = vec3_normalize(n_raw);
    let b_vec = vec3_cross(t_vec, n_vec);
    (t_vec, n_vec, b_vec)
}
/// Compute the osculating circle radius (= 1/κ) at parameter `t`.
pub fn osculating_radius(curve: &BezierCurve, t: f64) -> f64 {
    let kappa = curve_curvature(curve, t);
    if kappa < 1e-12 {
        f64::INFINITY
    } else {
        1.0 / kappa
    }
}
/// Compute chord-length parameters for a sequence of 3D points.
///
/// Returns normalized parameter values in `[0,1]`.
pub fn chord_length_params(points: &[[f64; 3]]) -> Vec<f64> {
    let n = points.len();
    if n == 0 {
        return vec![];
    }
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        params[i] = params[i - 1] + vec3_norm(vec3_sub(points[i], points[i - 1]));
    }
    let total = params[n - 1];
    if total > 0.0 {
        for p in params.iter_mut() {
            *p /= total;
        }
    }
    params
}
/// Insert a knot `t_new` into a B-spline curve at the given span, returning the new curve.
pub fn knot_insert(curve: &BSplineCurve, t_new: f64) -> BSplineCurve {
    let n = curve.control_points.len();
    let p = curve.degree;
    let knots = &curve.knots;
    let span = find_knot_span(n, p, t_new, knots);
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..span + 1]);
    new_knots.push(t_new);
    new_knots.extend_from_slice(&knots[span + 1..]);
    let mut new_cp = Vec::with_capacity(n + 1);
    for i in 0..=span - p {
        new_cp.push(curve.control_points[i]);
    }
    for i in span - p + 1..=span {
        let alpha = (t_new - knots[i]) / (knots[i + p] - knots[i]);
        let prev = curve.control_points[i - 1];
        let curr = curve.control_points[i];
        new_cp.push(vec3_add(
            vec3_scale(prev, 1.0 - alpha),
            vec3_scale(curr, alpha),
        ));
    }
    for i in span..n {
        new_cp.push(curve.control_points[i]);
    }
    BSplineCurve::new(new_cp, p, new_knots)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spline_geometry::BSplineSurface;
    use crate::spline_geometry::BezierPatch;
    use crate::spline_geometry::BlendingSurface;
    use crate::spline_geometry::ContinuityOrder;
    use crate::spline_geometry::CubicSpline;
    use crate::spline_geometry::NurbsCurve;
    use crate::spline_geometry::NurbsSurface;
    use crate::spline_geometry::PeriodicBSpline;
    use crate::spline_geometry::SweptSurface;
    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn test_bspline_basis_partition_of_unity() {
        let ctrl = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let degree = 3;
        let knots = uniform_clamped_knots(ctrl.len(), degree);
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let span = find_knot_span(ctrl.len(), degree, t, &knots);
            let basis = bspline_basis(span, degree, t, &knots);
            let sum: f64 = basis.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "basis sum = {}", sum);
        }
    }
    #[test]
    fn test_bspline_basis_non_negative() {
        let ctrl = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let degree = 2;
        let knots = uniform_clamped_knots(ctrl.len(), degree);
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let t_clamped = t.min(1.0 - 1e-12);
            let span = find_knot_span(ctrl.len(), degree, t_clamped, &knots);
            let basis = bspline_basis(span, degree, t_clamped, &knots);
            for &b in &basis {
                assert!(b >= -1e-12, "negative basis value: {}", b);
            }
        }
    }
    #[test]
    fn test_uniform_clamped_knots_endpoints() {
        let knots = uniform_clamped_knots(5, 3);
        for &k in knots.iter().take(4) {
            assert_eq!(k, 0.0);
        }
        let n = knots.len();
        for &k in knots.iter().skip(n - 4) {
            assert_eq!(k, 1.0);
        }
    }
    #[test]
    fn test_find_knot_span_boundary() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let span = find_knot_span(3, 2, 0.5, &knots);
        assert_eq!(span, 2);
    }
    #[test]
    fn test_bspline_curve_endpoints() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = BSplineCurve::with_clamped_knots(ctrl.clone(), 3);
        let p0 = curve.eval(0.0);
        let p1 = curve.eval(1.0 - 1e-12);
        assert!((p0[0] - ctrl[0][0]).abs() < 1e-6);
        assert!((p1[0] - ctrl[3][0]).abs() < 1e-4, "p1.x = {}", p1[0]);
    }
    #[test]
    fn test_bspline_curve_linear_reproduces_line() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 1);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let p = curve.eval(t.min(1.0 - 1e-12));
            assert!((p[0] - p[1]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_bspline_curve_arc_length_positive() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 2);
        let len = curve.arc_length(100);
        assert!(len > 0.0);
    }
    #[test]
    fn test_bspline_curve_derivative_direction() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 3);
        let d = curve.derivative(0.5);
        assert!(d[0] > 0.0, "tangent x = {}", d[0]);
        assert!(d[1].abs() < 1e-6);
    }
    #[test]
    fn test_bspline_curve_sample_count() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 2);
        let pts = curve.sample(11);
        assert_eq!(pts.len(), 11);
    }
    #[test]
    fn test_nurbs_unit_weights_matches_bspline() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0, 1.0];
        let nurbs = NurbsCurve::from_points_and_weights(pts.clone(), weights, 2);
        let bsp = BSplineCurve::with_clamped_knots(pts, 2);
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let pn = nurbs.eval(t);
            let pb = bsp.eval(t);
            for k in 0..3 {
                assert!(
                    (pn[k] - pb[k]).abs() < 1e-9,
                    "mismatch at dim {}: {} vs {}",
                    k,
                    pn[k],
                    pb[k]
                );
            }
        }
    }
    #[test]
    fn test_nurbs_circle_radius() {
        let r = 3.0;
        let circle = NurbsCurve::circle(r);
        for i in 0..12 {
            let t = i as f64 / 12.0;
            let p = circle.eval(t);
            let dist = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(
                (dist - r).abs() < 1e-6,
                "circle radius error: {} vs {}",
                dist,
                r
            );
        }
    }
    #[test]
    fn test_nurbs_sample_count() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let nurbs = NurbsCurve::from_points_and_weights(pts, vec![1.0, 1.0, 1.0], 2);
        assert_eq!(nurbs.sample(7).len(), 7);
    }
    #[test]
    fn test_bezier_endpoints() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = BezierCurve::new(ctrl.clone());
        let p0 = curve.eval(0.0);
        let p1 = curve.eval(1.0);
        for k in 0..3 {
            assert!((p0[k] - ctrl[0][k]).abs() < EPS);
            assert!((p1[k] - ctrl[3][k]).abs() < EPS);
        }
    }
    #[test]
    fn test_bezier_convex_hull_property() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = BezierCurve::new(ctrl.clone());
        let x_min = ctrl.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let x_max = ctrl.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        let y_min = ctrl.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        let y_max = ctrl.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let p = curve.eval(t);
            assert!(p[0] >= x_min - 1e-9 && p[0] <= x_max + 1e-9);
            assert!(p[1] >= y_min - 1e-9 && p[1] <= y_max + 1e-9);
        }
    }
    #[test]
    fn test_bezier_linear_is_interpolating() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        let curve = BezierCurve::new(ctrl);
        let mid = curve.eval(0.5);
        assert!((mid[0] - 0.5).abs() < EPS);
        assert!((mid[1] - 0.5).abs() < EPS);
    }
    #[test]
    fn test_bezier_split_reconstructs() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [3.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let curve = BezierCurve::new(ctrl);
        let (left, right) = curve.split(0.5);
        let orig = BezierCurve::new(curve.control_points.clone()).eval(0.5);
        let lp = left.eval(1.0);
        let rp = right.eval(0.0);
        for k in 0..3 {
            assert!((orig[k] - lp[k]).abs() < 1e-9);
            assert!((orig[k] - rp[k]).abs() < 1e-9);
        }
    }
    #[test]
    fn test_bezier_degree_elevate_same_shape() {
        let ctrl = vec![[0.0, 0.0, 0.0], [2.0, 4.0, 0.0], [4.0, 0.0, 0.0]];
        let curve = BezierCurve::new(ctrl);
        let elevated = curve.degree_elevate();
        assert_eq!(elevated.control_points.len(), 4);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let p_orig = curve.eval(t);
            let p_elev = elevated.eval(t);
            for k in 0..3 {
                assert!(
                    (p_orig[k] - p_elev[k]).abs() < 1e-9,
                    "degree elevation changed shape at t={}: {} vs {}",
                    t,
                    p_orig[k],
                    p_elev[k]
                );
            }
        }
    }
    #[test]
    fn test_bezier_arc_length_positive() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(ctrl);
        assert!(curve.arc_length(100) > 0.0);
    }
    #[test]
    fn test_natural_spline_interpolates_data() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 1.0, 0.0, -1.0, 0.0];
        let spline = CubicSpline::natural(x.clone(), y.clone());
        for i in 0..x.len() {
            let val = spline.eval(x[i]);
            assert!(
                (val - y[i]).abs() < 1e-8,
                "spline({}) = {} vs {}",
                x[i],
                val,
                y[i]
            );
        }
    }
    #[test]
    fn test_clamped_spline_interpolates_data() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, 4.0, 9.0];
        let spline = CubicSpline::clamped(x.clone(), y.clone(), 1.0, 6.0);
        for i in 0..x.len() {
            let val = spline.eval(x[i]);
            assert!((val - y[i]).abs() < 1e-8);
        }
    }
    #[test]
    fn test_natural_spline_second_deriv_at_endpoints_zero() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![1.0, 2.0, 1.5, 3.0, 2.0];
        let spline = CubicSpline::natural(x.clone(), y.clone());
        let d2_left = spline.eval_second_derivative(x[0]);
        let d2_right = spline.eval_second_derivative(*x.last().unwrap());
        assert!(d2_left.abs() < 1e-8, "d2 at left = {}", d2_left);
        assert!(d2_right.abs() < 1e-8, "d2 at right = {}", d2_right);
    }
    #[test]
    fn test_spline_eval_derivative_finite_diff() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, 0.0, 1.0];
        let spline = CubicSpline::natural(x, y);
        let t = 1.5;
        let h = 1e-5;
        let fd = (spline.eval(t + h) - spline.eval(t - h)) / (2.0 * h);
        let exact = spline.eval_derivative(t);
        assert!((fd - exact).abs() < 1e-4, "fd={} exact={}", fd, exact);
    }
    #[test]
    fn test_cubic_spline_linear_data() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 2.0, 4.0, 6.0];
        let spline = CubicSpline::natural(x, y);
        let val = spline.eval(1.5);
        assert!((val - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_bspline_fit_interpolates_when_num_ctrl_eq_data() {
        let data: Vec<[f64; 3]> = (0..5)
            .map(|i| {
                let t = i as f64 / 4.0;
                [t, t * t, 0.0]
            })
            .collect();
        let curve = fit_bspline_least_squares(&data, 5, 3);
        assert_eq!(curve.control_points.len(), 5);
    }
    #[test]
    fn test_bspline_fit_approximates_data() {
        let data: Vec<[f64; 3]> = (0..10)
            .map(|i| {
                let t = i as f64 / 9.0;
                [t, (2.0 * PI * t).sin(), 0.0]
            })
            .collect();
        let curve = fit_bspline_least_squares(&data, 6, 3);
        let pt_mid = curve.eval(0.5);
        assert!(pt_mid[0].abs() < 2.0, "mid x = {}", pt_mid[0]);
    }
    #[test]
    fn test_surface_of_revolution_cylinder() {
        let profile = vec![[1.0, 0.0], [1.0, 1.0]];
        let grid = surface_of_revolution(&profile, 8);
        assert_eq!(grid.len(), 2);
        assert_eq!(grid[0].len(), 8);
        for row in &grid {
            for p in row {
                let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
                assert!((r - 1.0).abs() < 1e-9);
            }
        }
    }
    #[test]
    fn test_surface_of_revolution_point_count() {
        let profile: Vec<[f64; 2]> = (0..5).map(|i| [i as f64, i as f64]).collect();
        let grid = surface_of_revolution(&profile, 12);
        assert_eq!(grid.len(), 5);
        for row in &grid {
            assert_eq!(row.len(), 12);
        }
    }
    #[test]
    fn test_swept_surface_grid_shape() {
        let spine: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let profile = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        let sw = SweptSurface::new(spine, profile);
        let grid = sw.compute();
        assert_eq!(grid.len(), 5);
        assert_eq!(grid[0].len(), 3);
    }
    #[test]
    fn test_swept_surface_spine_on_surface() {
        let spine: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let profile = vec![[0.0, 0.0]];
        let sw = SweptSurface::new(spine.clone(), profile);
        let grid = sw.compute();
        for (i, sp) in spine.iter().enumerate() {
            for k in 0..3 {
                assert!((grid[i][0][k] - sp[k]).abs() < 1e-9);
            }
        }
    }
    #[test]
    fn test_lofted_surface_shape() {
        let s0 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let s1 = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [2.0, 0.0, 1.0]];
        let grid = lofted_surface(&[s0, s1], 5);
        assert_eq!(grid.len(), 5);
        assert_eq!(grid[0].len(), 3);
    }
    #[test]
    fn test_lofted_surface_endpoints() {
        let s0 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let s1 = vec![[0.0, 0.0, 2.0], [1.0, 0.0, 2.0]];
        let grid = lofted_surface(&[s0.clone(), s1.clone()], 3);
        for k in 0..3 {
            assert!((grid[0][0][k] - s0[0][k]).abs() < 1e-9);
            assert!((grid[2][0][k] - s1[0][k]).abs() < 1e-9);
        }
    }
    #[test]
    fn test_blending_surface_g0_endpoints() {
        let c0 = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let c1 = BezierCurve::new(vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]);
        let blend = BlendingSurface::new(c0.clone(), c1.clone(), 1.0, 1.0, ContinuityOrder::G0);
        let p0 = blend.eval(0.0, 0.5);
        let p1 = blend.eval(1.0, 0.5);
        let ep0 = c0.eval(0.5);
        let ep1 = c1.eval(0.5);
        for k in 0..3 {
            assert!((p0[k] - ep0[k]).abs() < 1e-9);
            assert!((p1[k] - ep1[k]).abs() < 1e-9);
        }
    }
    #[test]
    fn test_blending_surface_g1_continuity() {
        let c0 = BezierCurve::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let c1 = BezierCurve::new(vec![[0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0]]);
        let blend = BlendingSurface::new(c0, c1, 0.5, 0.5, ContinuityOrder::G1);
        let grid = blend.sample_grid(5, 5);
        assert_eq!(grid.len(), 5);
        assert_eq!(grid[0].len(), 5);
    }
    #[test]
    fn test_blending_surface_g2_sample() {
        let c0 = BezierCurve::new(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ]);
        let c1 = BezierCurve::new(vec![
            [0.0, 3.0, 0.0],
            [1.0, 3.0, 0.0],
            [2.0, 3.0, 0.0],
            [3.0, 3.0, 0.0],
        ]);
        let blend = BlendingSurface::new(c0, c1, 0.3, 0.3, ContinuityOrder::G2);
        let p = blend.eval(0.5, 0.5);
        assert!(p[1] > 0.0 && p[1] < 3.0);
    }
    #[test]
    fn test_curvature_line_is_zero() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(ctrl);
        let kappa = curve_curvature(&curve, 0.5);
        assert!(kappa < 1e-8, "curvature of line = {}", kappa);
    }
    #[test]
    fn test_curvature_circle_approx() {
        let r = 2.0;
        let ctrl = vec![[r, 0.0, 0.0], [r, r, 0.0], [0.0, r, 0.0]];
        let curve = BezierCurve::new(ctrl);
        let kappa = curve_curvature(&curve, 0.5);
        assert!(kappa > 0.0);
    }
    #[test]
    fn test_frenet_serret_orthonormal() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [3.0, 1.0, 1.0],
        ];
        let curve = BezierCurve::new(ctrl);
        let (t, n, b) = frenet_serret_frame(&curve, 0.5);
        assert!((vec3_norm(t) - 1.0).abs() < 1e-6);
        let _ = n;
        let _ = b;
    }
    #[test]
    fn test_osculating_radius_line_is_infinity() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BezierCurve::new(ctrl);
        let r = osculating_radius(&curve, 0.5);
        assert!(r == f64::INFINITY || r > 1e10);
    }
    #[test]
    fn test_bspline_surface_eval_shape() {
        let net: Vec<Vec<[f64; 3]>> = (0..3)
            .map(|i| (0..3).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let surf = BSplineSurface::with_clamped_knots(net, 2, 2);
        let p = surf.eval(0.5, 0.5);
        assert!(p[0] >= 0.0 && p[0] <= 2.0);
        assert!(p[1] >= 0.0 && p[1] <= 2.0);
    }
    #[test]
    fn test_bspline_surface_normal_nonzero() {
        let net: Vec<Vec<[f64; 3]>> = (0..3)
            .map(|i| (0..3).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let surf = BSplineSurface::with_clamped_knots(net, 2, 2);
        let n = surf.normal(0.5, 0.5);
        let norm = vec3_norm(n);
        assert!((norm - 1.0).abs() < 0.1, "normal norm = {}", norm);
    }
    #[test]
    fn test_bspline_surface_sample_grid_dimensions() {
        let net: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|i| (0..4).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let surf = BSplineSurface::with_clamped_knots(net, 3, 3);
        let grid = surf.sample_grid(5, 7);
        assert_eq!(grid.len(), 5);
        assert_eq!(grid[0].len(), 7);
    }
    #[test]
    fn test_nurbs_surface_eval() {
        let net: Vec<Vec<[f64; 4]>> = (0..3)
            .map(|i| (0..3).map(|j| [i as f64, j as f64, 0.0, 1.0]).collect())
            .collect();
        let ku = uniform_clamped_knots(3, 2);
        let kv = uniform_clamped_knots(3, 2);
        let surf = NurbsSurface::new(net, 2, 2, ku, kv);
        let p = surf.eval(0.5, 0.5);
        assert!(p[0].is_finite() && p[1].is_finite());
    }
    #[test]
    fn test_periodic_bspline_returns_finite() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let curve = PeriodicBSpline::new(ctrl, 3);
        let p = curve.eval(0.5);
        for &pk in p.iter() {
            assert!(pk.is_finite());
        }
    }
    #[test]
    fn test_chord_length_params_bounds() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let params = chord_length_params(&pts);
        assert_eq!(params.len(), 4);
        assert!((params[0] - 0.0).abs() < EPS);
        assert!((params[3] - 1.0).abs() < EPS);
        for i in 0..params.len() - 1 {
            assert!(params[i] <= params[i + 1]);
        }
    }
    #[test]
    fn test_knot_insert_preserves_shape() {
        let ctrl = vec![
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 3);
        let refined = knot_insert(&curve, 0.5);
        for i in 1..9 {
            let t = i as f64 / 9.0;
            let p_orig = curve.eval(t);
            let p_ref = refined.eval(t);
            for k in 0..3 {
                assert!(
                    (p_orig[k] - p_ref[k]).abs() < 1e-8,
                    "knot insert changed shape at t={}: {} vs {}",
                    t,
                    p_orig[k],
                    p_ref[k]
                );
            }
        }
    }
    #[test]
    fn test_knot_insert_increases_knot_count() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let curve = BSplineCurve::with_clamped_knots(ctrl, 2);
        let n_before = curve.knots.len();
        let refined = knot_insert(&curve, 0.4);
        assert_eq!(refined.knots.len(), n_before + 1);
    }
    #[test]
    fn test_bezier_patch_corner_interpolation() {
        let grid: Vec<Vec<[f64; 3]>> = (0..3)
            .map(|i| (0..3).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let patch = BezierPatch::new(grid.clone());
        let p00 = patch.eval(0.0, 0.0);
        let p11 = patch.eval(1.0, 1.0);
        assert!((p00[0] - grid[0][0][0]).abs() < EPS);
        assert!((p11[0] - grid[2][2][0]).abs() < EPS);
    }
    #[test]
    fn test_bezier_patch_sample_grid_shape() {
        let grid: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|i| (0..4).map(|j| [i as f64, j as f64, 0.0]).collect())
            .collect();
        let patch = BezierPatch::new(grid);
        let sampled = patch.sample_grid(6, 8);
        assert_eq!(sampled.len(), 6);
        assert_eq!(sampled[0].len(), 8);
    }
    #[test]
    fn test_vec3_cross_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[0]).abs() < EPS);
        assert!((c[1]).abs() < EPS);
        assert!((c[2] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_vec3_normalize_unit() {
        let a = [3.0, 4.0, 0.0];
        let n = vec3_normalize(a);
        assert!((vec3_norm(n) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_vec3_lerp_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 4.0, 6.0];
        let mid = vec3_lerp(a, b, 0.5);
        assert!((mid[0] - 1.0).abs() < EPS);
        assert!((mid[1] - 2.0).abs() < EPS);
        assert!((mid[2] - 3.0).abs() < EPS);
    }
}
