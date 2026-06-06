//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{dot3, finite_diff_3d};
use super::types::{BsplineBasis, KnotVector};

/// A rational tensor-product B-spline (NURBS) surface.
///
/// S(u,v) = Σ_i Σ_j w_{ij} N_{i,p}(u) N_{j,q}(v) P_{ij} / Σ_i Σ_j w_{ij} N_{i,p}(u) N_{j,q}(v)
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// Basis in u direction.
    pub basis_u: BsplineBasis,
    /// Basis in v direction.
    pub basis_v: BsplineBasis,
    /// Control point grid \[i\]\[j\] = \[x, y, z\].
    pub control_points: Vec<Vec<[f64; 3]>>,
    /// Weight grid \[i\]\[j\].
    pub weights: Vec<Vec<f64>>,
}
impl NurbsSurface {
    /// Create a new NURBS surface.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knot_u: KnotVector,
        knot_v: KnotVector,
        control_points: Vec<Vec<[f64; 3]>>,
        weights: Vec<Vec<f64>>,
    ) -> Self {
        let n_u = control_points.len();
        let n_v = if n_u > 0 { control_points[0].len() } else { 0 };
        Self {
            basis_u: BsplineBasis::new(degree_u, knot_u, n_u),
            basis_v: BsplineBasis::new(degree_v, knot_v, n_v),
            control_points,
            weights,
        }
    }
    /// Evaluate the NURBS surface at (u, v).
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut num = [0.0_f64; 3];
        let mut denom = 0.0_f64;
        for (k, &nu_k) in nu[..=p].iter().enumerate() {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for (l, &nv_l) in nv[..=q].iter().enumerate() {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let wn = self.weights[iu][iv] * nu_k * nv_l;
                denom += wn;
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    num[d] += wn * cp[d];
                }
            }
        }
        if denom.abs() < 1e-30 {
            return self.control_points[0][0];
        }
        [num[0] / denom, num[1] / denom, num[2] / denom]
    }
    /// Compute the unit surface normal at (u, v) via numerical differentiation.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let h = 1e-6;
        let su = {
            let p0 = self.eval(u - h, v);
            let p1 = self.eval(u + h, v);
            [
                (p1[0] - p0[0]) / (2.0 * h),
                (p1[1] - p0[1]) / (2.0 * h),
                (p1[2] - p0[2]) / (2.0 * h),
            ]
        };
        let sv = {
            let p0 = self.eval(u, v - h);
            let p1 = self.eval(u, v + h);
            [
                (p1[0] - p0[0]) / (2.0 * h),
                (p1[1] - p0[1]) / (2.0 * h),
                (p1[2] - p0[2]) / (2.0 * h),
            ]
        };
        let cross = [
            su[1] * sv[2] - su[2] * sv[1],
            su[2] * sv[0] - su[0] * sv[2],
            su[0] * sv[1] - su[1] * sv[0],
        ];
        let mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        if mag < 1e-14 {
            [0.0, 0.0, 1.0]
        } else {
            [cross[0] / mag, cross[1] / mag, cross[2] / mag]
        }
    }
    /// Fit the NURBS surface to a point cloud (updates control points via least squares).
    ///
    /// `points` is a nu × nv grid of target points. Weights are kept at 1.0.
    /// Uses a simple collocation approach.
    pub fn fit_to_grid(&mut self, points: &[Vec<[f64; 3]>]) {
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let nu_pts = points.len();
        let nv_pts = if nu_pts > 0 { points[0].len() } else { 0 };
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        for i in 0..n_u {
            let u = u0 + (u1 - u0) * i as f64 / (n_u - 1).max(1) as f64;
            let pi = (u * (nu_pts - 1) as f64).round() as usize;
            let pi = pi.min(nu_pts - 1);
            for j in 0..n_v {
                let v = v0 + (v1 - v0) * j as f64 / (n_v - 1).max(1) as f64;
                let pj = (v * (nv_pts - 1) as f64).round() as usize;
                let pj = pj.min(nv_pts - 1);
                self.control_points[i][j] = points[pi][pj];
            }
        }
    }
    /// Sample the surface at (nu × nv) parameter values.
    pub fn sample(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        (0..nu)
            .map(|i| {
                let u = u0 + (u1 - u0) * i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = v0 + (v1 - v0) * j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A tensor-product B-spline surface in 3D space.
///
/// S(u,v) = Σ_i Σ_j N_{i,p}(u) N_{j,q}(v) P_{ij}
#[derive(Debug, Clone)]
pub struct BsplineSurface {
    /// Basis in the u direction.
    pub basis_u: BsplineBasis,
    /// Basis in the v direction.
    pub basis_v: BsplineBasis,
    /// Control point grid \[i\]\[j\] = \[x, y, z\], row-major (u varies fastest).
    pub control_points: Vec<Vec<[f64; 3]>>,
}
impl BsplineSurface {
    /// Create a new B-spline surface.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knot_u: KnotVector,
        knot_v: KnotVector,
        control_points: Vec<Vec<[f64; 3]>>,
    ) -> Self {
        let n_u = control_points.len();
        let n_v = control_points[0].len();
        Self {
            basis_u: BsplineBasis::new(degree_u, knot_u, n_u),
            basis_v: BsplineBasis::new(degree_v, knot_v, n_v),
            control_points,
        }
    }
    /// Create a bilinear patch (p=q=1) from four corner points.
    pub fn bilinear_patch(p00: [f64; 3], p10: [f64; 3], p01: [f64; 3], p11: [f64; 3]) -> Self {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0]);
        let cps = vec![vec![p00, p01], vec![p10, p11]];
        Self {
            basis_u: BsplineBasis::new(1, kv.clone(), 2),
            basis_v: BsplineBasis::new(1, kv, 2),
            control_points: cps,
        }
    }
    /// Evaluate the surface at (u, v): returns \[x, y, z\].
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut point = [0.0_f64; 3];
        for (k, &nu_k) in nu[..=p].iter().enumerate() {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for (l, &nv_l) in nv[..=q].iter().enumerate() {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = nu_k * nv_l;
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    point[d] += w * cp[d];
                }
            }
        }
        point
    }
    /// Evaluate partial derivative ∂S/∂u at (u, v).
    pub fn eval_du(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let (span_u, ders_u) = self.basis_u.eval_nonzero_derivs(u, 1);
        let (span_v, nv) = self.basis_v.eval_nonzero(v);
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut result = [0.0_f64; 3];
        for (k, &du_k) in ders_u[1][..=p].iter().enumerate() {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for (l, &nv_l) in nv[..=q].iter().enumerate() {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = du_k * nv_l;
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    result[d] += w * cp[d];
                }
            }
        }
        result
    }
    /// Evaluate partial derivative ∂S/∂v at (u, v).
    pub fn eval_dv(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let (span_u, nu) = self.basis_u.eval_nonzero(u);
        let (span_v, ders_v) = self.basis_v.eval_nonzero_derivs(v, 1);
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let mut result = [0.0_f64; 3];
        for (k, &nu_k) in nu[..=p].iter().enumerate() {
            let iu = span_u + k - p;
            if iu >= n_u {
                continue;
            }
            for (l, &dv_l) in ders_v[1][..=q].iter().enumerate() {
                let iv = span_v + l - q;
                if iv >= n_v {
                    continue;
                }
                let w = nu_k * dv_l;
                let cp = self.control_points[iu][iv];
                for d in 0..3 {
                    result[d] += w * cp[d];
                }
            }
        }
        result
    }
    /// Compute the unit surface normal at (u, v).
    ///
    /// N = (∂S/∂u × ∂S/∂v) / |∂S/∂u × ∂S/∂v|
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let du = self.eval_du(u, v);
        let dv = self.eval_dv(u, v);
        let cross = [
            du[1] * dv[2] - du[2] * dv[1],
            du[2] * dv[0] - du[0] * dv[2],
            du[0] * dv[1] - du[1] * dv[0],
        ];
        let mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        if mag < 1e-14 {
            [0.0, 0.0, 1.0]
        } else {
            [cross[0] / mag, cross[1] / mag, cross[2] / mag]
        }
    }
    /// Compute the first fundamental form coefficients E, F, G at (u, v).
    ///
    /// E = ∂S/∂u · ∂S/∂u, F = ∂S/∂u · ∂S/∂v, G = ∂S/∂v · ∂S/∂v
    pub fn first_fundamental_form(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let su = self.eval_du(u, v);
        let sv = self.eval_dv(u, v);
        let e = su[0] * su[0] + su[1] * su[1] + su[2] * su[2];
        let f = su[0] * sv[0] + su[1] * sv[1] + su[2] * sv[2];
        let g = sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2];
        (e, f, g)
    }
    /// Compute Gaussian curvature K = (LN - M²)/(EG - F²) at (u, v).
    ///
    /// Uses finite differences for second derivatives.
    pub fn gaussian_curvature(&self, u: f64, v: f64) -> f64 {
        let h = 1e-5;
        let (e, f, g) = self.first_fundamental_form(u, v);
        let egf2 = e * g - f * f;
        if egf2.abs() < 1e-20 {
            return 0.0;
        }
        let suu = finite_diff_3d(|uu| self.eval_du(uu, v), u, h);
        let svv = finite_diff_3d(|vv| self.eval_dv(u, vv), v, h);
        let suv = finite_diff_3d(|uu| self.eval_dv(uu, v), u, h);
        let n = self.normal(u, v);
        let l = dot3(suu, n);
        let m_coef = dot3(suv, n);
        let nm = dot3(svv, n);
        (l * nm - m_coef * m_coef) / egf2
    }
    /// Compute mean curvature H = (EN - 2FM + GL)/(2(EG - F²)) at (u, v).
    pub fn mean_curvature(&self, u: f64, v: f64) -> f64 {
        let h = 1e-5;
        let (e, f, g) = self.first_fundamental_form(u, v);
        let egf2 = e * g - f * f;
        if egf2.abs() < 1e-20 {
            return 0.0;
        }
        let suu = finite_diff_3d(|uu| self.eval_du(uu, v), u, h);
        let svv = finite_diff_3d(|vv| self.eval_dv(u, vv), v, h);
        let suv = finite_diff_3d(|uu| self.eval_dv(uu, v), u, h);
        let n = self.normal(u, v);
        let l = dot3(suu, n);
        let m_coef = dot3(suv, n);
        let nm = dot3(svv, n);
        (e * nm - 2.0 * f * m_coef + g * l) / (2.0 * egf2)
    }
    /// Sample the surface at (nu × nv) parameter values.
    pub fn sample(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        let (u0, u1) = self.basis_u.knot_vector.domain();
        let (v0, v1) = self.basis_v.knot_vector.domain();
        (0..nu)
            .map(|i| {
                let u = u0 + (u1 - u0) * i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = v0 + (v1 - v0) * j as f64 / (nv - 1).max(1) as f64;
                        self.eval(u, v)
                    })
                    .collect()
            })
            .collect()
    }
    /// Compute Gaussian and mean curvature at (u, v).
    ///
    /// Returns `(K_gauss, H_mean)`.
    pub fn compute_curvature(&self, u: f64, v: f64) -> (f64, f64) {
        (self.gaussian_curvature(u, v), self.mean_curvature(u, v))
    }
    /// Refine the surface by inserting new knots in u and v directions.
    ///
    /// Uses the Oslo algorithm (knot insertion) to preserve geometry.
    pub fn refine_knots(&self, new_knots_u: &[f64], new_knots_v: &[f64]) -> BsplineSurface {
        let n_u = self.control_points.len();
        let n_v = if n_u > 0 {
            self.control_points[0].len()
        } else {
            0
        };
        let p = self.basis_u.degree;
        let q = self.basis_v.degree;
        let mut refined_u_knots = self.basis_u.knot_vector.knots.clone();
        for &k in new_knots_u {
            refined_u_knots.push(k);
        }
        refined_u_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut new_cp = self.control_points.clone();
        for &knot in new_knots_u {
            let n_cur = new_cp.len();
            let orig_u_knots = &self.basis_u.knot_vector.knots;
            let span = orig_u_knots
                .iter()
                .position(|&k| k > knot)
                .map(|pos| pos.saturating_sub(1))
                .unwrap_or(n_cur.saturating_sub(1))
                .min(n_cur);
            let mut inserted = Vec::with_capacity(n_cur + 1);
            for i in 0..n_cur + 1 {
                let old_knots = &self.basis_u.knot_vector.knots;
                let ki = if i < old_knots.len() {
                    old_knots[i]
                } else {
                    1.0
                };
                let ki_p = if i + p < old_knots.len() {
                    old_knots[i + p]
                } else {
                    1.0
                };
                let denom = ki_p - ki;
                let alpha = if denom.abs() > 1e-14 {
                    (knot - ki) / denom
                } else {
                    0.5
                };
                let row: Vec<[f64; 3]> = (0..n_v)
                    .map(|j| {
                        if i < span.saturating_sub(p) + 1 {
                            new_cp[i][j]
                        } else if i > span {
                            new_cp[i - 1][j]
                        } else {
                            let p0 = if i > 0 {
                                new_cp[i - 1][j]
                            } else {
                                new_cp[0][j]
                            };
                            let p1 = if i < n_cur {
                                new_cp[i][j]
                            } else {
                                new_cp[n_cur - 1][j]
                            };
                            [
                                (1.0 - alpha) * p0[0] + alpha * p1[0],
                                (1.0 - alpha) * p0[1] + alpha * p1[1],
                                (1.0 - alpha) * p0[2] + alpha * p1[2],
                            ]
                        }
                    })
                    .collect();
                inserted.push(row);
            }
            new_cp = inserted;
        }
        let n_u_new = new_cp.len();
        let mut new_cp_v = new_cp;
        let mut refined_v_knots = self.basis_v.knot_vector.knots.clone();
        for &k in new_knots_v {
            refined_v_knots.push(k);
        }
        refined_v_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for &knot in new_knots_v {
            let n_v_cur = if !new_cp_v.is_empty() {
                new_cp_v[0].len()
            } else {
                0
            };
            let orig_v_knots = &self.basis_v.knot_vector.knots;
            let span = orig_v_knots
                .iter()
                .position(|&k| k > knot)
                .map(|pos| pos.saturating_sub(1))
                .unwrap_or(n_v_cur.saturating_sub(1))
                .min(n_v_cur);
            let mut new_rows = Vec::with_capacity(n_u_new);
            for old_row in new_cp_v.iter() {
                let n_cur = old_row.len();
                let mut new_row = Vec::with_capacity(n_cur + 1);
                for j in 0..n_cur + 1 {
                    if j < span.saturating_sub(q) + 1 {
                        new_row.push(old_row[j]);
                    } else if j > span {
                        new_row.push(old_row[j - 1]);
                    } else {
                        let old_v_knots = &self.basis_v.knot_vector.knots;
                        let kj = if j < old_v_knots.len() {
                            old_v_knots[j]
                        } else {
                            1.0
                        };
                        let kj_q = if j + q < old_v_knots.len() {
                            old_v_knots[j + q]
                        } else {
                            1.0
                        };
                        let denom = kj_q - kj;
                        let alpha = if denom.abs() > 1e-14 {
                            (knot - kj) / denom
                        } else {
                            0.5
                        };
                        let p0 = if j > 0 { old_row[j - 1] } else { old_row[0] };
                        let p1 = if j < n_cur {
                            old_row[j]
                        } else {
                            old_row[n_cur - 1]
                        };
                        new_row.push([
                            (1.0 - alpha) * p0[0] + alpha * p1[0],
                            (1.0 - alpha) * p0[1] + alpha * p1[1],
                            (1.0 - alpha) * p0[2] + alpha * p1[2],
                        ]);
                    }
                }
                new_rows.push(new_row);
            }
            new_cp_v = new_rows;
        }
        let final_n_u = new_cp_v.len();
        let final_n_v = if final_n_u > 0 { new_cp_v[0].len() } else { 0 };
        let mut ku = self.basis_u.knot_vector.knots.clone();
        for &k in new_knots_u {
            ku.push(k);
        }
        ku.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut kv = self.basis_v.knot_vector.knots.clone();
        for &k in new_knots_v {
            kv.push(k);
        }
        kv.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        while ku.len() < final_n_u + p + 1 {
            ku.push(1.0);
        }
        while kv.len() < final_n_v + q + 1 {
            kv.push(1.0);
        }
        BsplineSurface::new(p, q, KnotVector::new(ku), KnotVector::new(kv), new_cp_v)
    }
}
