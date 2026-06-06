//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Compute cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Compute squared norm of a 3-vector.
#[inline]
pub fn norm_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}
/// Compute norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    norm_sq3(a).sqrt()
}
/// Normalize a 3-vector; returns zero vector if near-zero.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}
/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Find the knot span index i such that `knot[i] <= t < knot[i+1]`.
///
/// Uses binary search. For t == knot\[n+1\] (end of parameter domain),
/// returns n (the last internal span index) per NURBS convention.
///
/// # Arguments
/// * `n` – number of basis functions minus 1 (n = num_control_points - 1)
/// * `p` – polynomial degree
/// * `t` – parameter value in \[knot\[p\\], knot\[n+1\]]
/// * `knot` – knot vector of length n+p+2
pub fn find_knot_span(n: usize, p: usize, t: f64, knot: &[f64]) -> usize {
    if (t - knot[n + 1]).abs() < 1e-14 {
        return n;
    }
    let mut low = p;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while t < knot[mid] || t >= knot[mid + 1] {
        if t < knot[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
        if mid >= knot.len() - 1 {
            break;
        }
    }
    mid
}
/// Compute B-spline basis functions N_{i,p}(t) using Cox-de Boor recursion.
///
/// Returns a vector of p+1 values: N_{i-p,p}(t), ..., N_{i,p}(t)
/// where `i` is the knot span index.
///
/// # Arguments
/// * `i` – knot span index (result of `find_knot_span`)
/// * `p` – polynomial degree
/// * `t` – parameter value
/// * `knot` – knot vector
pub fn basis_functions(i: usize, p: usize, t: f64, knot: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0f64; p + 1];
    let mut left = vec![0.0f64; p + 1];
    let mut right = vec![0.0f64; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = t - knot[i + 1 - j];
        right[j] = knot[i + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}
/// Compute B-spline basis function derivatives up to order `d`.
///
/// Returns a 2D array `ders[k][j]` = d^k N_{i-p+j,p}(t) / dt^k
/// for k=0..d and j=0..p.
///
/// # Arguments
/// * `i` – knot span index
/// * `p` – polynomial degree
/// * `t` – parameter value
/// * `knot` – knot vector
/// * `d` – maximum derivative order
pub fn basis_derivatives(i: usize, p: usize, t: f64, knot: &[f64], d: usize) -> Vec<Vec<f64>> {
    let d_ord = d.min(p);
    let mut ders = vec![vec![0.0f64; p + 1]; d_ord + 1];
    let mut ndu = vec![vec![0.0f64; p + 1]; p + 1];
    let mut left = vec![0.0f64; p + 1];
    let mut right = vec![0.0f64; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knot[i + 1 - j];
        right[j] = knot[i + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = vec![vec![0.0f64; p + 1]; 2];
    for r in 0..=p {
        let mut s1 = 0usize;
        let mut s2 = 1usize;
        a[0][0] = 1.0;
        for k in 1..=d_ord {
            let mut d_val = 0.0;
            let rk = r as isize - k as isize;
            let pk = p as isize - k as isize;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[(pk + 1) as usize][rk as usize];
                d_val = a[s2][0] * ndu[rk as usize][(pk) as usize];
            }
            let j1 = if rk >= -1 {
                1usize
            } else {
                (-(rk + 1)) as usize
            };
            let j2 = if (r as isize - 1) <= pk {
                k - 1
            } else {
                (p as isize - r as isize) as usize
            };
            for j in j1..=j2 {
                let idx_rk = (rk + j as isize) as usize;
                let idx_pk = (pk + 1) as usize;
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[idx_pk][idx_rk];
                d_val += a[s2][j] * ndu[idx_rk][pk as usize];
            }
            if (r as isize) <= pk {
                a[s2][k] = -a[s1][k - 1] / ndu[(pk + 1) as usize][r];
                d_val += a[s2][k] * ndu[r][pk as usize];
            }
            ders[k][r] = d_val;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut rp = p as f64;
    for (k, ders_k) in ders.iter_mut().enumerate().take(d_ord + 1).skip(1) {
        for val in ders_k.iter_mut().take(p + 1) {
            *val *= rp;
        }
        rp *= (p - k) as f64;
    }
    ders
}
/// Evaluate a NURBS curve point at parameter t.
///
/// Given control points, weights, knot vector, degree and parameter,
/// returns the point in 3D space using homogeneous coordinates.
///
/// # Arguments
/// * `ctrl` – control points as slice of \[f64; 3\]
/// * `weights` – corresponding weights
/// * `knot` – knot vector
/// * `p` – degree
/// * `t` – parameter
pub fn nurbs_point(ctrl: &[[f64; 3]], weights: &[f64], knot: &[f64], p: usize, t: f64) -> [f64; 3] {
    let n = ctrl.len() - 1;
    let span = find_knot_span(n, p, t, knot);
    let basis = basis_functions(span, p, t, knot);
    let mut w_sum = 0.0;
    let mut pw = [0.0f64; 3];
    for (j, &basis_j) in basis.iter().enumerate().take(p + 1) {
        let idx = span - p + j;
        let w = weights[idx] * basis_j;
        pw[0] += w * ctrl[idx][0];
        pw[1] += w * ctrl[idx][1];
        pw[2] += w * ctrl[idx][2];
        w_sum += w;
    }
    if w_sum.abs() < 1e-15 {
        return pw;
    }
    [pw[0] / w_sum, pw[1] / w_sum, pw[2] / w_sum]
}
/// Return Gauss-Legendre quadrature points and weights for n_pts points.
///
/// Points are in \[-1, 1\].
pub fn gauss_legendre_1d(n_pts: usize) -> (Vec<f64>, Vec<f64>) {
    match n_pts {
        1 => (vec![0.0], vec![2.0]),
        2 => (
            vec![-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()],
            vec![1.0, 1.0],
        ),
        3 => (
            vec![-0.7745966692414834, 0.0, 0.7745966692414834],
            vec![
                0.5555555555555556,
                0.888_888_888_888_889,
                0.5555555555555556,
            ],
        ),
        4 => (
            vec![
                -0.8611363115940526,
                -0.3399810435848563,
                0.3399810435848563,
                0.8611363115940526,
            ],
            vec![
                0.3478548451374538,
                0.6521451548625461,
                0.6521451548625461,
                0.3478548451374538,
            ],
        ),
        _ => {
            let mut pts = Vec::with_capacity(n_pts);
            let mut wts = Vec::with_capacity(n_pts);
            let h = 2.0 / n_pts as f64;
            for i in 0..n_pts {
                pts.push(-1.0 + (i as f64 + 0.5) * h);
                wts.push(h);
            }
            (pts, wts)
        }
    }
}
/// Create an open uniform (clamped) B-spline knot vector.
///
/// The knot vector has p+1 repetitions at 0 and 1, with uniform interior knots.
///
/// # Arguments
/// * `n` – number of control points minus 1
/// * `p` – polynomial degree
pub fn open_uniform_knot(n: usize, p: usize) -> Vec<f64> {
    let m = n + p + 1;
    let mut knot = vec![0.0f64; m + 1];
    for knot_i in &mut knot[(m - p)..=m] {
        *knot_i = 1.0;
    }
    let n_interior = m - 2 * p;
    if n_interior > 0 {
        for i in 1..=n_interior {
            knot[p + i] = i as f64 / (n_interior + 1) as f64;
        }
    }
    knot
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::isogeometric::*;
    fn make_line_curve() -> NurbsCurve {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0, 1.0];
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        NurbsCurve::new(ctrl, weights, 2, knot)
    }
    fn make_flat_surface() -> NurbsSurface {
        let ctrl = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        NurbsSurface::new(ctrl, weights, 1, 1, knot.clone(), knot)
    }
    #[test]
    fn test_find_knot_span_basic() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        let n = 4;
        let p = 2;
        let span = find_knot_span(n, p, 0.5, &knot);
        assert_eq!(span, 2);
        let span2 = find_knot_span(n, p, 1.5, &knot);
        assert_eq!(span2, 3);
    }
    #[test]
    fn test_find_knot_span_end() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let n = 2;
        let p = 2;
        let span = find_knot_span(n, p, 1.0, &knot);
        assert_eq!(span, n);
    }
    #[test]
    fn test_basis_functions_partition_of_unity() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let p = 2;
        let n = 2;
        for t in [0.0, 0.25, 0.5, 0.75, 1.0 - 1e-10] {
            let span = find_knot_span(n, p, t, &knot);
            let basis = basis_functions(span, p, t, &knot);
            let sum: f64 = basis.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-12,
                "partition of unity failed at t={}, sum={}",
                t,
                sum
            );
        }
    }
    #[test]
    fn test_basis_derivatives_degree_one() {
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let p = 1;
        let n = 1;
        let t = 0.5;
        let span = find_knot_span(n, p, t, &knot);
        let ders = basis_derivatives(span, p, t, &knot, 1);
        let d_sum: f64 = ders[1].iter().sum();
        assert!(
            d_sum.abs() < 1e-12,
            "derivative sum should be 0, got {}",
            d_sum
        );
    }
    #[test]
    fn test_nurbs_point_line() {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let p = make_nurbs_point_at(0.5, &ctrl, &weights, &knot, 1);
        assert!((p[0] - 0.5).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
        assert!(p[2].abs() < 1e-10);
    }
    fn make_nurbs_point_at(
        t: f64,
        ctrl: &[[f64; 3]],
        weights: &[f64],
        knot: &[f64],
        p: usize,
    ) -> [f64; 3] {
        nurbs_point(ctrl, weights, knot, p, t)
    }
    #[test]
    fn test_bspline_basis_new() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let basis = BsplineBasis::new(2, knot);
        assert_eq!(basis.degree, 2);
        assert_eq!(basis.num_basis, 3);
    }
    #[test]
    fn test_bspline_basis_domain() {
        let knot = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let basis = BsplineBasis::new(1, knot);
        let (t0, t1) = basis.domain();
        assert!((t0 - 0.0).abs() < 1e-14);
        assert!((t1 - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_nurbs_curve_point_at_endpoints() {
        let curve = make_line_curve();
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0 - 1e-10);
        assert!((p0[0] - 0.0).abs() < 1e-6, "start point x should be 0");
        assert!(p1[0] > 1.5, "end point x should be near 2");
    }
    #[test]
    fn test_nurbs_curve_derivative_nonzero() {
        let curve = make_line_curve();
        let d = curve.derivative_at(0.5);
        assert!(
            d[0].abs() > 0.1,
            "derivative x should be nonzero, got {:?}",
            d
        );
    }
    #[test]
    fn test_nurbs_curve_insert_knot() {
        let curve = make_line_curve();
        let new_curve = curve.insert_knot(0.5);
        assert_eq!(
            new_curve.control_points.len(),
            curve.control_points.len() + 1
        );
        let p_orig = curve.point_at(0.3);
        let p_new = new_curve.point_at(0.3);
        for (&orig_i, &new_i) in p_orig.iter().zip(p_new.iter()) {
            assert!(
                (orig_i - new_i).abs() < 1e-10,
                "geometry changed after knot insert"
            );
        }
    }
    #[test]
    fn test_nurbs_surface_point_at_corners() {
        let surf = make_flat_surface();
        let p00 = surf.point_at(0.0, 0.0);
        let p11 = surf.point_at(1.0 - 1e-10, 1.0 - 1e-10);
        assert!((p00[0] - 0.0).abs() < 1e-6);
        assert!((p00[1] - 0.0).abs() < 1e-6);
        assert!(p11[0] > 0.5);
        assert!(p11[1] > 0.5);
    }
    #[test]
    fn test_nurbs_surface_normal_direction() {
        let surf = make_flat_surface();
        let (_, _, n) = surf.normal(0.5, 0.5);
        assert!(
            n[2].abs() > 0.5,
            "normal should be in z direction, got {:?}",
            n
        );
    }
    #[test]
    fn test_nurbs_volume_creation() {
        let ctrl: Vec<Vec<Vec<[f64; 3]>>> = vec![
            vec![
                vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                vec![[0.0, 1.0, 0.0], [0.0, 1.0, 1.0]],
            ],
            vec![
                vec![[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
                vec![[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
            ],
        ];
        let weights = vec![
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        ];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let vol = NurbsVolume::new(ctrl, weights, 1, 1, 1, knot.clone(), knot.clone(), knot);
        let p = vol.point_at(0.5, 0.5, 0.5);
        assert!((p[0] - 0.5).abs() < 0.1);
        assert!((p[1] - 0.5).abs() < 0.1);
        assert!((p[2] - 0.5).abs() < 0.1);
    }
    #[test]
    fn test_nurbs_volume_jacobian() {
        let ctrl: Vec<Vec<Vec<[f64; 3]>>> = vec![
            vec![
                vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                vec![[0.0, 1.0, 0.0], [0.0, 1.0, 1.0]],
            ],
            vec![
                vec![[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
                vec![[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
            ],
        ];
        let weights = vec![
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        ];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let vol = NurbsVolume::new(ctrl, weights, 1, 1, 1, knot.clone(), knot.clone(), knot);
        let j_det = vol.jacobian_det(0.5, 0.5, 0.5);
        assert!(
            j_det > 0.0,
            "Jacobian determinant should be positive for unit cube"
        );
    }
    #[test]
    fn test_open_uniform_knot() {
        let knot = open_uniform_knot(3, 2);
        assert_eq!(knot.len(), 7);
        assert_eq!(knot[0], 0.0);
        assert_eq!(knot[1], 0.0);
        assert_eq!(knot[2], 0.0);
        assert_eq!(knot[5], 1.0);
        assert_eq!(knot[6], 1.0);
        for &k in &knot {
            assert!(
                (0.0..=1.0).contains(&k),
                "knot value should be in [0,1]: {}",
                k
            );
        }
        for i in 1..knot.len() {
            assert!(knot[i] >= knot[i - 1], "knot should be non-decreasing");
        }
    }
    #[test]
    fn test_iga_element_shape_functions_sum() {
        let surf = make_flat_surface();
        let connectivity = vec![0, 1, 2, 3];
        let elem = IgaElement::new(0, connectivity, surf, 3, 3);
        let (r, _, _) = elem.shape_functions(0.5, 0.5);
        let sum: f64 = r.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "shape functions should sum to 1, got {}",
            sum
        );
    }
    #[test]
    fn test_iga_assembler_element_stiffness_size() {
        let surf = make_flat_surface();
        let connectivity = vec![0, 1, 2, 3];
        let elem = IgaElement::new(0, connectivity, surf, 3, 3);
        let assembler = IgaAssembler::new(2e11, 0.3, 7800.0, 0.01, 3);
        let k = assembler.element_stiffness(&elem, (0.0, 1.0), (0.0, 1.0));
        assert_eq!(k.len(), 64, "stiffness matrix should be 8x8");
    }
    #[test]
    fn test_iga_assembler_element_mass_positive() {
        let surf = make_flat_surface();
        let connectivity = vec![0, 1, 2, 3];
        let elem = IgaElement::new(0, connectivity, surf, 3, 3);
        let assembler = IgaAssembler::new(2e11, 0.3, 7800.0, 0.01, 3);
        let m = assembler.element_mass(&elem, (0.0, 1.0), (0.0, 1.0));
        let sum: f64 = m.iter().sum();
        assert!(sum > 0.0, "mass matrix sum should be positive");
    }
    #[test]
    fn test_patch_connectivity_new() {
        let conn = PatchConnectivity::new(
            (0, 1),
            CouplingType::Penalty(1e6),
            vec![(0, 3), (1, 2)],
            0,
            100,
        );
        assert_eq!(conn.num_pairs(), 2);
        assert!(conn.is_interface_a(0));
        assert!(!conn.is_interface_a(5));
    }
    #[test]
    fn test_patch_connectivity_penalty() {
        let conn = PatchConnectivity::new((0, 1), CouplingType::Penalty(1e6), vec![(0, 1)], 0, 50);
        assert!((conn.penalty_coefficient(0) - 1e6).abs() < 1.0);
    }
    #[test]
    fn test_knot_refinement_h_refine() {
        let curve = make_line_curve();
        let refined = KnotRefinement::h_refine_curve(&curve, &[0.25, 0.5, 0.75]);
        assert!(refined.control_points.len() > curve.control_points.len());
        let p_orig = curve.point_at(0.4);
        let p_ref = refined.point_at(0.4);
        for i in 0..3 {
            assert!(
                (p_orig[i] - p_ref[i]).abs() < 1e-8,
                "h-refinement changed geometry"
            );
        }
    }
    #[test]
    fn test_knot_refinement_p_refine() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let new_knot = KnotRefinement::p_refine_knot(&knot, 2);
        assert!(!new_knot.is_empty());
    }
    #[test]
    fn test_knot_refinement_k_refine() {
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let new_knot = KnotRefinement::k_refine(&knot, 1, 2);
        assert!(new_knot.len() > knot.len());
    }
    #[test]
    fn test_iga_contact_gap_frictionless() {
        let s1 = make_flat_surface();
        let ctrl2 = vec![
            vec![[0.0, 0.0, -0.1], [1.0, 0.0, -0.1]],
            vec![[0.0, 1.0, -0.1], [1.0, 1.0, -0.1]],
        ];
        let w2 = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let s2 = NurbsSurface::new(ctrl2, w2, 1, 1, knot.clone(), knot);
        let contact = IgaContactFormulation::new(1e5, 0.0, s1, s2, 2);
        let (gap, n) = contact.compute_gap(0.5, 0.5);
        let _ = gap;
        let n_len = norm3(n);
        assert!(n_len > 0.5, "normal should be non-zero");
    }
    #[test]
    fn test_iga_adaptive_creation() {
        let surf = make_flat_surface();
        let adaptive = IgaAdaptive::new(surf);
        assert_eq!(adaptive.num_elements(), 1);
        assert_eq!(adaptive.max_level(), 0);
    }
    #[test]
    fn test_iga_adaptive_refinement() {
        let surf = make_flat_surface();
        let mut adaptive = IgaAdaptive::new(surf);
        adaptive.refine_element(0);
        assert_eq!(
            adaptive.num_elements(),
            4,
            "one element should split into 4"
        );
        assert_eq!(adaptive.max_level(), 1);
    }
    #[test]
    fn test_iga_adaptive_error_indicator() {
        let surf = make_flat_surface();
        let adaptive = IgaAdaptive::new(surf);
        let err = adaptive.compute_error_indicator(0);
        assert!(err >= 0.0, "error indicator should be non-negative");
    }
    #[test]
    fn test_gauss_legendre_weights_sum() {
        for n in 1..=4 {
            let (_, w) = gauss_legendre_1d(n);
            let sum: f64 = w.iter().sum();
            assert!(
                (sum - 2.0).abs() < 1e-10,
                "Gauss weights sum should be 2 for n={}",
                n
            );
        }
    }
    #[test]
    fn test_bspline_basis_evaluate() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let basis = BsplineBasis::new(2, knot);
        let (span, vals) = basis.evaluate(0.5);
        let sum: f64 = vals.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-12,
            "basis should sum to 1 at t=0.5, span={}",
            span
        );
    }
    #[test]
    fn test_nurbs_curve_arc_length() {
        let curve = make_line_curve();
        let len = curve.arc_length(100);
        assert!(
            (len - 2.0).abs() < 0.05,
            "arc length should be ~2, got {}",
            len
        );
    }
    #[test]
    fn test_basis_functions_at_zero() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let p = 2;
        let n = 2;
        let span = find_knot_span(n, p, 0.0, &knot);
        let basis = basis_functions(span, p, 0.0, &knot);
        assert!((basis[0] - 1.0).abs() < 1e-12, "N_0,2(0) should be 1");
    }
}
/// Evaluate a cubic B-spline basis function from 5 local knot values.
///
/// Implements Cox-de Boor for a single cubic B-spline from knot sequence.
pub fn cubic_bspline_basis(t: f64, knots: &[f64; 5]) -> f64 {
    bspline_n(t, knots, 0, 3)
}
/// Recursive B-spline basis evaluation N_{i,p}(t) over a local knot array.
pub fn bspline_n(t: f64, knots: &[f64; 5], i: usize, p: usize) -> f64 {
    if p == 0 {
        if t >= knots[i] && t < knots[i + 1] {
            1.0
        } else {
            0.0
        }
    } else {
        let denom1 = knots[i + p] - knots[i];
        let denom2 = knots[i + p + 1] - knots[i + 1];
        let left = if denom1.abs() > 1e-14 {
            (t - knots[i]) / denom1 * bspline_n(t, knots, i, p - 1)
        } else {
            0.0
        };
        let right = if denom2.abs() > 1e-14 {
            (knots[i + p + 1] - t) / denom2 * bspline_n(t, knots, i + 1, p - 1)
        } else {
            0.0
        };
        left + right
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::isogeometric::*;
    fn make_line_curve() -> NurbsCurve {
        let ctrl = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let weights = vec![1.0, 1.0, 1.0];
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        NurbsCurve::new(ctrl, weights, 2, knot)
    }
    fn make_flat_surface() -> NurbsSurface {
        let ctrl = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        NurbsSurface::new(ctrl, weights, 1, 1, knot.clone(), knot)
    }
    fn make_unit_volume() -> NurbsVolume {
        let ctrl: Vec<Vec<Vec<[f64; 3]>>> = vec![
            vec![
                vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                vec![[0.0, 1.0, 0.0], [0.0, 1.0, 1.0]],
            ],
            vec![
                vec![[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
                vec![[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
            ],
        ];
        let weights = vec![
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        ];
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        NurbsVolume::new(ctrl, weights, 1, 1, 1, knot.clone(), knot.clone(), knot)
    }
    #[test]
    fn bezier_extraction_basic_creation() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let be = BezierExtraction::compute(&knot, 2);
        assert_eq!(be.degree, 2);
        assert!(be.num_elements >= 1);
    }
    #[test]
    fn bezier_extraction_operator_size() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let be = BezierExtraction::compute(&knot, 2);
        let op = be.operator(0);
        assert_eq!(op.len(), 9);
    }
    #[test]
    fn bezier_extraction_ctrl_extraction() {
        let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let be = BezierExtraction::compute(&knot, 2);
        let ctrl_pts = vec![1.0, 2.0, 3.0];
        let bezier_ctrl = be.extract_bezier_ctrl(0, &ctrl_pts, 0);
        assert_eq!(bezier_ctrl.len(), 3);
    }
    #[test]
    fn bezier_extraction_multi_element() {
        let knot = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        let be = BezierExtraction::compute(&knot, 2);
        assert!(be.num_elements >= 1);
    }
    #[test]
    fn brep_add_face() {
        let mut brep = NurbsBRep::new();
        let surf = make_flat_surface();
        let idx = brep.add_face(surf);
        assert_eq!(idx, 0);
        assert_eq!(brep.num_faces(), 1);
    }
    #[test]
    fn brep_add_edge_and_vertex() {
        let mut brep = NurbsBRep::new();
        let v0 = brep.add_vertex([0.0, 0.0, 0.0]);
        let v1 = brep.add_vertex([1.0, 0.0, 0.0]);
        let curve = make_line_curve();
        let e0 = brep.add_edge(curve, v0, v1);
        assert_eq!(e0, 0);
        assert_eq!(brep.num_vertices(), 2);
        assert_eq!(brep.num_edges(), 1);
    }
    #[test]
    fn brep_validity_empty() {
        let brep = NurbsBRep::new();
        assert!(brep.is_valid());
    }
    #[test]
    fn brep_total_area_positive() {
        let mut brep = NurbsBRep::new();
        let surf = make_flat_surface();
        brep.add_face(surf);
        let area = brep.total_surface_area(3);
        assert!(area > 0.0, "area should be positive, got {area}");
    }
    #[test]
    fn brep_attach_edge_to_face() {
        let mut brep = NurbsBRep::new();
        let surf = make_flat_surface();
        let face_idx = brep.add_face(surf);
        let v0 = brep.add_vertex([0.0, 0.0, 0.0]);
        let v1 = brep.add_vertex([1.0, 0.0, 0.0]);
        let curve = make_line_curve();
        let edge_idx = brep.add_edge(curve, v0, v1);
        brep.attach_edge_to_face(face_idx, edge_idx);
        assert_eq!(brep.face_edges[face_idx].len(), 1);
        assert!(brep.is_valid());
    }
    #[test]
    fn param_quality_flat_surface_valid() {
        let surf = make_flat_surface();
        let q = ParameterizationQuality::analyze_surface(&surf, 4);
        assert!(
            q.is_valid(),
            "flat surface parameterization should be valid"
        );
    }
    #[test]
    fn param_quality_min_leq_max() {
        let surf = make_flat_surface();
        let q = ParameterizationQuality::analyze_surface(&surf, 4);
        assert!(q.min_jacobian <= q.max_jacobian);
    }
    #[test]
    fn param_quality_score_in_range() {
        let surf = make_flat_surface();
        let q = ParameterizationQuality::analyze_surface(&surf, 4);
        let score = q.quality_score();
        assert!((0.0..=1.0).contains(&score), "score={score}");
    }
    #[test]
    fn param_quality_num_samples() {
        let surf = make_flat_surface();
        let q = ParameterizationQuality::analyze_surface(&surf, 5);
        assert_eq!(q.num_samples, 25);
    }
    #[test]
    fn tspline_control_point_creation() {
        let s_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let t_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let cp = TSplineControlPoint::new([1.0, 0.0, 0.0], 1.0, s_knots, t_knots);
        assert!((cp.weight - 1.0).abs() < 1e-14);
    }
    #[test]
    fn tspline_surface_point_at() {
        let s_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let t_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let cp = TSplineControlPoint::new([1.0, 1.0, 0.0], 1.0, s_knots, t_knots);
        let surf = TSplineSurface::new(vec![cp], [0.0, 1.0, 0.0, 1.0]);
        let _pt = surf.point_at(0.5, 0.5);
        assert_eq!(surf.num_control_points(), 1);
    }
    #[test]
    fn tspline_no_t_junctions_uniform() {
        let s_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let t_knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        let cps = vec![
            TSplineControlPoint::new([0.0, 0.0, 0.0], 1.0, s_knots, t_knots),
            TSplineControlPoint::new([1.0, 0.0, 0.0], 1.0, s_knots, t_knots),
        ];
        let surf = TSplineSurface::new(cps, [0.0, 1.0, 0.0, 1.0]);
        assert!(
            !surf.has_t_junctions(),
            "uniform knots should have no T-junctions"
        );
    }
    #[test]
    fn tspline_has_t_junctions_different_knots() {
        let sk1 = [0.0, 0.0, 0.5, 1.0, 1.0];
        let tk1 = [0.0, 0.0, 0.5, 1.0, 1.0];
        let sk2 = [0.0, 0.0, 0.3, 1.0, 1.0];
        let tk2 = [0.0, 0.0, 0.5, 1.0, 1.0];
        let cps = vec![
            TSplineControlPoint::new([0.0, 0.0, 0.0], 1.0, sk1, tk1),
            TSplineControlPoint::new([1.0, 0.0, 0.0], 1.0, sk2, tk2),
        ];
        let surf = TSplineSurface::new(cps, [0.0, 1.0, 0.0, 1.0]);
        assert!(surf.has_t_junctions());
    }
    #[test]
    fn degree_elevation_increases_degree() {
        let curve = make_line_curve();
        let elevated = DegreeElevation::elevate_once(&curve);
        assert_eq!(elevated.basis.degree, 3);
    }
    #[test]
    fn degree_elevation_twice() {
        let curve = make_line_curve();
        let elevated = DegreeElevation::elevate(&curve, 2);
        assert_eq!(elevated.basis.degree, 4);
    }
    #[test]
    fn degree_elevation_geometry_preserved() {
        let curve = make_line_curve();
        let elevated = DegreeElevation::elevate_once(&curve);
        let p_orig = curve.point_at(0.5);
        let p_elev = elevated.point_at(0.5);
        assert!(
            (p_orig[0] - p_elev[0]).abs() < 0.5,
            "x should be close: orig={}, elev={}",
            p_orig[0],
            p_elev[0]
        );
    }
    #[test]
    fn volume_mass_matrix_size() {
        let vol = make_unit_volume();
        let assembler = IgaVolumeMassMatrix::new(7800.0, 2);
        let m = assembler.assemble(&vol);
        assert_eq!(
            m.len(),
            24 * 24,
            "mass matrix should be 24x24 for unit volume"
        );
    }
    #[test]
    fn volume_mass_matrix_positive_diagonal() {
        let vol = make_unit_volume();
        let assembler = IgaVolumeMassMatrix::new(1.0, 3);
        let m = assembler.assemble(&vol);
        let n_dof = 24;
        for i in 0..n_dof {
            assert!(
                m[i * n_dof + i] >= 0.0,
                "diagonal entry {i} should be non-negative"
            );
        }
    }
    #[test]
    fn k_refinement_strategy_apply() {
        let knot = vec![0.0, 0.0, 1.0, 1.0];
        let strategy = KRefinementStrategy::new(2, 1, true, 0.01);
        let new_knot = strategy.apply(&knot, 1);
        assert!(
            new_knot.len() > knot.len(),
            "k-refinement should grow knot vector"
        );
    }
    #[test]
    fn k_refinement_dof_estimate() {
        let strategy = KRefinementStrategy::new(3, 2, true, 0.01);
        let dofs = strategy.estimate_dofs_1d(5);
        assert!(dofs >= 5, "estimated DOFs should be >= original");
    }
    #[test]
    fn k_refinement_possible_check() {
        let strategy = KRefinementStrategy::new(2, 1, true, 0.1);
        assert!(strategy.is_refinement_possible(0.5));
        assert!(!strategy.is_refinement_possible(0.05));
    }
    #[test]
    fn mortar_coupling_matrix_size() {
        let s1 = make_flat_surface();
        let s2 = make_flat_surface();
        let coupling = MortarPatchCoupling::new(s1, s2, 3, 1.0);
        let d = coupling.coupling_matrix();
        assert_eq!(d.len(), 4);
    }
    #[test]
    fn mortar_interface_gap_zero_for_identical() {
        let s1 = make_flat_surface();
        let s2 = make_flat_surface();
        let coupling = MortarPatchCoupling::new(s1, s2, 5, 1.0);
        let gaps = coupling.interface_gap(5);
        assert_eq!(gaps.len(), 5);
        for (_t, gap) in &gaps {
            let gap_len = (gap[0] * gap[0] + gap[1] * gap[1] + gap[2] * gap[2]).sqrt();
            assert!(
                gap_len < 1e-10,
                "identical surfaces should have zero gap, got {gap_len}"
            );
        }
    }
}
#[cfg(test)]
mod tests_nurbs_basis_and_assembly {

    use crate::isogeometric::*;
    fn linear_basis() -> NurbsBasis {
        NurbsBasis::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![1.0, 1.0])
    }
    fn quadratic_basis() -> NurbsBasis {
        NurbsBasis::new(2, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0])
    }
    #[test]
    fn nurbs_basis_num_basis() {
        let b = quadratic_basis();
        assert_eq!(b.num_basis(), 3);
    }
    #[test]
    fn nurbs_basis_degree() {
        let b = quadratic_basis();
        assert_eq!(b.degree(), 2);
    }
    #[test]
    fn nurbs_basis_domain() {
        let b = quadratic_basis();
        let (lo, hi) = b.domain();
        assert!((lo - 0.0).abs() < 1e-12);
        assert!((hi - 1.0).abs() < 1e-12);
    }
    #[test]
    fn nurbs_basis_partition_of_unity() {
        let b = quadratic_basis();
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let (_span, r) = b.evaluate(t.min(1.0 - 1e-12));
            let sum: f64 = r.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "partition of unity at t={t}: sum={sum}"
            );
        }
    }
    #[test]
    fn nurbs_basis_non_negative() {
        let b = quadratic_basis();
        for i in 0..=10 {
            let t = (i as f64 / 10.0).min(1.0 - 1e-12);
            let (_span, r) = b.evaluate(t);
            for &v in &r {
                assert!(v >= -1e-12, "basis value should be non-negative, got {v}");
            }
        }
    }
    #[test]
    fn nurbs_basis_linear_partition_of_unity() {
        let b = linear_basis();
        for i in 0..=5 {
            let t = (i as f64 / 5.0).min(1.0 - 1e-12);
            let (_span, r) = b.evaluate(t);
            let sum: f64 = r.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn nurbs_basis_evaluate_rational_returns_span_r_dr() {
        let b = quadratic_basis();
        let (span, r, dr) = b.evaluate_rational(0.5);
        assert!(span < b.bspline.knot.len());
        assert_eq!(r.len(), b.degree() + 1);
        assert_eq!(dr.len(), b.degree() + 1);
    }
    #[test]
    fn nurbs_basis_derivative_sum_zero() {
        let b = quadratic_basis();
        let (_, _, dr) = b.evaluate_rational(0.5);
        let sum_dr: f64 = dr.iter().sum();
        assert!(
            sum_dr.abs() < 1e-9,
            "sum of derivatives should be 0, got {sum_dr}"
        );
    }
    #[test]
    fn nurbs_basis_curve_point_endpoints() {
        let b = quadratic_basis();
        let cps = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let pt0 = b.curve_point(&cps, 1e-12);
        let pt1 = b.curve_point(&cps, 1.0 - 1e-12);
        assert!(pt0[0] < 0.1, "start near first CP");
        assert!(pt1[0] > 1.9, "end near last CP");
    }
    #[test]
    fn nurbs_basis_de_boor_midpoint() {
        let b = quadratic_basis();
        let cvs = [0.0, 1.0, 2.0];
        let val = b.de_boor_scalar(&cvs, 0.5);
        assert!(
            (val - 1.0).abs() < 0.5,
            "midpoint value should be near 1.0, got {val}"
        );
    }
    #[test]
    fn nurbs_basis_insert_knot_increases_num_basis() {
        let b = quadratic_basis();
        let b2 = b.insert_knot(0.5);
        assert!(
            b2.num_basis() > b.num_basis(),
            "knot insertion increases num_basis"
        );
    }
    #[test]
    fn nurbs_basis_uniform_open_knot_length() {
        let knot = NurbsBasis::uniform_open_knot(2, 4);
        assert_eq!(knot.len(), 4 + 2 + 1, "knot length = n_basis + degree + 1");
    }
    #[test]
    fn nurbs_basis_uniform_open_knot_clamped() {
        let p = 3;
        let n = 5;
        let knot = NurbsBasis::uniform_open_knot(p, n);
        for &k in knot.iter().take(p) {
            assert_eq!(k, 0.0, "start of knot clamped to 0");
        }
        let m = knot.len();
        for &k in knot.iter().skip(m - p) {
            assert_eq!(k, 1.0, "end of knot clamped to 1");
        }
    }
    #[test]
    fn iga_assembly_new_empty() {
        let asm = IgaAssembly::new(4);
        assert_eq!(asm.n_dof, 4);
        assert!(asm.k_coo.is_empty());
        assert!(asm.m_coo.is_empty());
    }
    #[test]
    fn iga_assembly_global_stiffness_size() {
        let asm = IgaAssembly::new(3);
        let k = asm.global_stiffness();
        assert_eq!(k.len(), 9);
    }
    #[test]
    fn iga_assembly_global_mass_size() {
        let asm = IgaAssembly::new(3);
        let m = asm.global_mass();
        assert_eq!(m.len(), 9);
    }
    #[test]
    fn iga_assembly_add_stiffness_accumulates() {
        let mut asm = IgaAssembly::new(3);
        let dof_map = vec![0, 1];
        let k_loc = vec![1.0, 0.5, 0.5, 2.0];
        asm.add_stiffness(&dof_map, &k_loc);
        let k = asm.global_stiffness();
        assert!((k[0] - 1.0).abs() < 1e-12, "K[0,0] = 1.0");
        assert!((k[1] - 0.5).abs() < 1e-12, "K[0,1] = 0.5");
    }
    #[test]
    fn iga_assembly_add_mass_accumulates() {
        let mut asm = IgaAssembly::new(3);
        let dof_map = vec![1, 2];
        let m_loc = vec![2.0, 0.0, 0.0, 3.0];
        asm.add_mass(&dof_map, &m_loc);
        let m = asm.global_mass();
        assert!((m[4] - 2.0).abs() < 1e-12, "M[1,1] = 2.0");
        assert!((m[8] - 3.0).abs() < 1e-12, "M[2,2] = 3.0");
    }
    #[test]
    fn iga_assembly_clear_resets() {
        let mut asm = IgaAssembly::new(2);
        asm.add_stiffness(&[0, 1], &[1.0, 0.0, 0.0, 1.0]);
        asm.clear();
        assert!(asm.k_coo.is_empty());
        assert!(asm.m_coo.is_empty());
    }
    #[test]
    fn iga_assembly_global_stiffness_symmetric() {
        let mut asm = IgaAssembly::new(2);
        let k_loc = vec![4.0, 2.0, 2.0, 6.0];
        asm.add_stiffness(&[0, 1], &k_loc);
        let k = asm.global_stiffness();
        assert!((k[1] - k[2]).abs() < 1e-12, "K should be symmetric");
    }
    #[test]
    fn iga_assembly_add_stiffness_multiple_elements() {
        let mut asm = IgaAssembly::new(3);
        asm.add_stiffness(&[0, 1], &[1.0, 0.0, 0.0, 1.0]);
        asm.add_stiffness(&[1, 2], &[2.0, 0.0, 0.0, 2.0]);
        let k = asm.global_stiffness();
        assert!((k[4] - 3.0).abs() < 1e-12, "K[1,1] = 3.0 from two elements");
    }
}
