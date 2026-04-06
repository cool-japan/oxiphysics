//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    add3, arbitrary_perp, cross3, dist3, dot3, lerp3, normalize3, scale3, sub3,
};

/// A tube mesh built by extruding a circular cross-section along a spine curve.
///
/// Uses the Frenet frame at each spine sample to orient the cross-section
/// circles, then connects them into a triangle strip.
pub struct TubeGeometry {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (each triple is one triangle, CCW from outside).
    pub triangles: Vec<[usize; 3]>,
}
impl TubeGeometry {
    /// Build a tube along a `BezierCurve`.
    ///
    /// # Parameters
    /// - `curve`:          The spine curve.
    /// - `radius`:         Radius of the circular cross-section.
    /// - `n_spine`:        Number of spine samples (at least 2).
    /// - `n_radial`:       Number of vertices around each ring (at least 3).
    pub fn from_bezier(curve: &BezierCurve, radius: f64, n_spine: usize, n_radial: usize) -> Self {
        assert!(n_spine >= 2, "n_spine must be >= 2");
        assert!(n_radial >= 3, "n_radial must be >= 3");
        use std::f64::consts::TAU;
        let mut vertices: Vec<[f64; 3]> = Vec::with_capacity(n_spine * n_radial);
        let mut triangles: Vec<[usize; 3]> = Vec::new();
        for spine_i in 0..n_spine {
            let t = spine_i as f64 / (n_spine - 1) as f64;
            let center = curve.evaluate(t);
            let frame = FrenetFrame::from_bezier(curve, t);
            for rad_j in 0..n_radial {
                let angle = TAU * rad_j as f64 / n_radial as f64;
                let (sin_a, cos_a) = angle.sin_cos();
                let offset = add3(
                    scale3(frame.normal, radius * cos_a),
                    scale3(frame.binormal, radius * sin_a),
                );
                vertices.push(add3(center, offset));
            }
        }
        for spine_i in 0..n_spine - 1 {
            let ring0 = spine_i * n_radial;
            let ring1 = (spine_i + 1) * n_radial;
            for rad_j in 0..n_radial {
                let next_j = (rad_j + 1) % n_radial;
                let a = ring0 + rad_j;
                let b = ring0 + next_j;
                let c = ring1 + rad_j;
                let d = ring1 + next_j;
                triangles.push([a, b, d]);
                triangles.push([a, d, c]);
            }
        }
        TubeGeometry {
            vertices,
            triangles,
        }
    }
}
/// An orthonormal tangent frame at a surface point.
#[derive(Debug, Clone)]
pub struct TangentFrame {
    /// Surface normal (unit vector).
    pub normal: [f64; 3],
    /// Tangent along u direction (unit vector).
    pub tangent_u: [f64; 3],
    /// Tangent along v direction (unit vector, ≈ normal × tangent_u).
    pub tangent_v: [f64; 3],
}
impl TangentFrame {
    /// Compute the tangent frame of a `BezierSurface` at `(u, v)`.
    pub fn from_bezier_surface(surf: &BezierSurface, u: f64, v: f64) -> Self {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(0.0);
        let u1 = (u + eps).min(1.0);
        let v0 = (v - eps).max(0.0);
        let v1 = (v + eps).min(1.0);
        let du_raw = scale3(
            sub3(surf.evaluate(u1, v), surf.evaluate(u0, v)),
            1.0 / (u1 - u0),
        );
        let dv_raw = scale3(
            sub3(surf.evaluate(u, v1), surf.evaluate(u, v0)),
            1.0 / (v1 - v0),
        );
        let normal = normalize3(cross3(du_raw, dv_raw));
        let tangent_u = normalize3(du_raw);
        let tangent_v = normalize3(cross3(normal, tangent_u));
        TangentFrame {
            normal,
            tangent_u,
            tangent_v,
        }
    }
    /// Compute the tangent frame of a `LoftSurface` at `(u, v)`.
    pub fn from_loft_surface(surf: &LoftSurface, u: f64, v: f64) -> Self {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(0.0);
        let u1 = (u + eps).min(1.0);
        let v0 = (v - eps).max(0.0);
        let v1 = (v + eps).min(1.0);
        let du_raw = scale3(
            sub3(surf.evaluate(u1, v), surf.evaluate(u0, v)),
            1.0 / (u1 - u0),
        );
        let dv_raw = scale3(
            sub3(surf.evaluate(u, v1), surf.evaluate(u, v0)),
            1.0 / (v1 - v0),
        );
        let normal = normalize3(cross3(du_raw, dv_raw));
        let tangent_u = normalize3(du_raw);
        let tangent_v = normalize3(cross3(normal, tangent_u));
        TangentFrame {
            normal,
            tangent_u,
            tangent_v,
        }
    }
}
/// A NURBS surface defined by a control net of weighted 3-D points.
///
/// Each element `control_net[i][j]` is `[x, y, z, w]` (homogeneous coordinates).
pub struct NurbsSurface {
    /// Control net: `control_net[i][j] = [x, y, z, w]`.
    pub control_net: Vec<Vec<[f64; 4]>>,
    /// Knot vector in the u direction.
    pub u_knots: Vec<f64>,
    /// Knot vector in the v direction.
    pub v_knots: Vec<f64>,
    /// Degree in the u direction.
    pub degree_u: usize,
    /// Degree in the v direction.
    pub degree_v: usize,
}
impl NurbsSurface {
    /// Construct a new NURBS surface.
    pub fn new(
        control_net: Vec<Vec<[f64; 4]>>,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        degree_u: usize,
        degree_v: usize,
    ) -> Self {
        NurbsSurface {
            control_net,
            u_knots,
            v_knots,
            degree_u,
            degree_v,
        }
    }
    /// Evaluate the NURBS surface at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let n_u = self.control_net.len();
        let n_v = if n_u == 0 {
            0
        } else {
            self.control_net[0].len()
        };
        let mut num = [0.0f64; 3];
        let mut den = 0.0f64;
        for i in 0..n_u {
            let bu = NurbsCurve::b_spline_basis(i, self.degree_u, u, &self.u_knots);
            for j in 0..n_v {
                let bv = NurbsCurve::b_spline_basis(j, self.degree_v, v, &self.v_knots);
                let pt = self.control_net[i][j];
                let w = pt[3];
                let wb = w * bu * bv;
                num[0] += wb * pt[0];
                num[1] += wb * pt[1];
                num[2] += wb * pt[2];
                den += wb;
            }
        }
        if den.abs() < 1e-12 {
            [0.0, 0.0, 0.0]
        } else {
            [num[0] / den, num[1] / den, num[2] / den]
        }
    }
}
/// A tensor-product B-spline surface defined by a rectangular control net.
///
/// `control_net[i][j]` is the `(i, j)` control point.  Knot vectors in
/// each parameter direction are stored separately.
pub struct BsplineSurface {
    /// Control net — outer index u, inner index v.
    pub control_net: Vec<Vec<[f64; 3]>>,
    /// Knot vector in the u direction.
    pub u_knots: Vec<f64>,
    /// Knot vector in the v direction.
    pub v_knots: Vec<f64>,
    /// Polynomial degree in u.
    pub degree_u: usize,
    /// Polynomial degree in v.
    pub degree_v: usize,
}
impl BsplineSurface {
    /// Construct a `BsplineSurface`.
    pub fn new(
        control_net: Vec<Vec<[f64; 3]>>,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        degree_u: usize,
        degree_v: usize,
    ) -> Self {
        BsplineSurface {
            control_net,
            u_knots,
            v_knots,
            degree_u,
            degree_v,
        }
    }
    /// Evaluate the surface at parameter `(u, v)` using tensor-product de Boor.
    pub fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let n_u = self.control_net.len();
        let n_v = if n_u == 0 {
            0
        } else {
            self.control_net[0].len()
        };
        let mut result = [0.0_f64; 3];
        let mut denom = 0.0_f64;
        for i in 0..n_u {
            let bu = NurbsCurve::b_spline_basis(i, self.degree_u, u, &self.u_knots);
            for j in 0..n_v {
                let bv = NurbsCurve::b_spline_basis(j, self.degree_v, v, &self.v_knots);
                let b = bu * bv;
                let pt = self.control_net[i][j];
                result[0] += b * pt[0];
                result[1] += b * pt[1];
                result[2] += b * pt[2];
                denom += b;
            }
        }
        if denom.abs() < 1e-14 {
            [0.0; 3]
        } else {
            scale3(result, 1.0 / denom)
        }
    }
    /// Partial derivative with respect to u (finite differences).
    fn deriv_u(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(self.u_knots[0]);
        let u1 = (u + eps).min(*self.u_knots.last().unwrap_or(&1.0));
        scale3(sub3(self.eval(u1, v), self.eval(u0, v)), 1.0 / (u1 - u0))
    }
    /// Partial derivative with respect to v (finite differences).
    fn deriv_v(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let v0 = (v - eps).max(self.v_knots[0]);
        let v1 = (v + eps).min(*self.v_knots.last().unwrap_or(&1.0));
        scale3(sub3(self.eval(u, v1), self.eval(u, v0)), 1.0 / (v1 - v0))
    }
    /// Second partial derivative with respect to u (finite differences).
    fn deriv_uu(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(self.u_knots[0]);
        let u1 = (u + eps).min(*self.u_knots.last().unwrap_or(&1.0));
        let du0 = self.deriv_u(u0, v);
        let du1 = self.deriv_u(u1, v);
        scale3(sub3(du1, du0), 1.0 / (u1 - u0))
    }
    /// Second partial derivative with respect to v (finite differences).
    fn deriv_vv(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let v0 = (v - eps).max(self.v_knots[0]);
        let v1 = (v + eps).min(*self.v_knots.last().unwrap_or(&1.0));
        let dv0 = self.deriv_v(u, v0);
        let dv1 = self.deriv_v(u, v1);
        scale3(sub3(dv1, dv0), 1.0 / (v1 - v0))
    }
    /// Mixed partial derivative ∂²S/∂u∂v (finite differences).
    fn deriv_uv(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let v0 = (v - eps).max(self.v_knots[0]);
        let v1 = (v + eps).min(*self.v_knots.last().unwrap_or(&1.0));
        let du_v0 = self.deriv_u(u, v0);
        let du_v1 = self.deriv_u(u, v1);
        scale3(sub3(du_v1, du_v0), 1.0 / (v1 - v0))
    }
    /// Compute Gaussian curvature K and mean curvature H at `(u, v)`.
    ///
    /// Uses the first and second fundamental forms computed via finite
    /// differences of the surface position.
    ///
    /// Returns `(K, H)`.  Both values may be `NaN` at degenerate points.
    pub fn compute_curvature(&self, u: f64, v: f64) -> (f64, f64) {
        let su = self.deriv_u(u, v);
        let sv = self.deriv_v(u, v);
        let suu = self.deriv_uu(u, v);
        let svv = self.deriv_vv(u, v);
        let suv = self.deriv_uv(u, v);
        let normal_raw = cross3(su, sv);
        let normal_len =
            (normal_raw[0].powi(2) + normal_raw[1].powi(2) + normal_raw[2].powi(2)).sqrt();
        if normal_len < 1e-14 {
            return (0.0, 0.0);
        }
        let n = scale3(normal_raw, 1.0 / normal_len);
        let big_e = dot3(su, su);
        let big_f = dot3(su, sv);
        let big_g = dot3(sv, sv);
        let w = big_e * big_g - big_f * big_f;
        if w.abs() < 1e-20 {
            return (0.0, 0.0);
        }
        let big_l = dot3(suu, n);
        let big_m = dot3(suv, n);
        let big_n = dot3(svv, n);
        let gaussian = (big_l * big_n - big_m * big_m) / w;
        let mean = (big_e * big_n - 2.0 * big_f * big_m + big_g * big_l) / (2.0 * w);
        (gaussian, mean)
    }
    /// H-refinement via knot insertion in both parameter directions.
    ///
    /// `new_u_knots` and `new_v_knots` are sorted lists of knot values to
    /// insert (repeated insertions are handled by the Boehm algorithm).
    ///
    /// Returns a new `BsplineSurface` with the refined control net.
    pub fn refine_knots(&self, new_u_knots: &[f64], new_v_knots: &[f64]) -> Self {
        let refined_u = self.refine_u(new_u_knots);
        refined_u.refine_v(new_v_knots)
    }
    /// Insert a single knot `t` in the u direction using the Boehm algorithm
    /// applied to each row of the control net independently.
    fn insert_u_knot(&self, t: f64) -> Self {
        let n_u = self.control_net.len();
        let n_v = if n_u == 0 {
            0
        } else {
            self.control_net[0].len()
        };
        let p = self.degree_u;
        let knots = &self.u_knots;
        let mut k = p;
        for idx in p..(knots.len() - p - 1) {
            if knots[idx] <= t && t < knots[idx + 1] {
                k = idx;
                break;
            }
        }
        let mut new_knots = knots.clone();
        new_knots.insert(k + 1, t);
        let mut new_net: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n_u + 1);
        for j in 0..n_v {
            let row: Vec<[f64; 3]> = (0..n_u).map(|i| self.control_net[i][j]).collect();
            let mut new_row: Vec<[f64; 3]> = Vec::with_capacity(n_u + 1);
            for i in 0..=(n_u) {
                if i <= k - p {
                    new_row.push(row[i]);
                } else if i <= k {
                    let alpha = if (knots[i + p] - knots[i]).abs() < 1e-14 {
                        0.0
                    } else {
                        (t - knots[i]) / (knots[i + p] - knots[i])
                    };
                    let prev = if i == 0 { [0.0; 3] } else { row[i - 1] };
                    let curr = row[i.min(row.len() - 1)];
                    new_row.push(add3(scale3(prev, 1.0 - alpha), scale3(curr, alpha)));
                } else {
                    new_row.push(row[(i - 1).min(row.len() - 1)]);
                }
            }
            for i in 0..=(n_u) {
                while new_net.len() <= i {
                    new_net.push(Vec::with_capacity(n_v));
                }
                new_net[i].push(if i < new_row.len() {
                    new_row[i]
                } else {
                    [0.0; 3]
                });
            }
            let _ = row;
        }
        BsplineSurface {
            control_net: new_net,
            u_knots: new_knots,
            v_knots: self.v_knots.clone(),
            degree_u: self.degree_u,
            degree_v: self.degree_v,
        }
    }
    /// Refine all u-direction knots in `new_knots` sequentially.
    fn refine_u(&self, new_knots: &[f64]) -> Self {
        let mut surf = BsplineSurface {
            control_net: self.control_net.clone(),
            u_knots: self.u_knots.clone(),
            v_knots: self.v_knots.clone(),
            degree_u: self.degree_u,
            degree_v: self.degree_v,
        };
        for &t in new_knots {
            surf = surf.insert_u_knot(t);
        }
        surf
    }
    /// Insert a single knot `t` in the v direction using the Boehm algorithm
    /// applied to each column of the control net independently.
    fn insert_v_knot(&self, t: f64) -> Self {
        let n_u = self.control_net.len();
        let p = self.degree_v;
        let knots = &self.v_knots;
        let mut k = p;
        for idx in p..(knots.len().saturating_sub(p + 1)) {
            if knots[idx] <= t && t < knots[idx + 1] {
                k = idx;
                break;
            }
        }
        let mut new_knots = knots.clone();
        new_knots.insert(k + 1, t);
        let mut new_net: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n_u);
        for i in 0..n_u {
            let row = &self.control_net[i];
            let n_v = row.len();
            let mut new_row: Vec<[f64; 3]> = Vec::with_capacity(n_v + 1);
            for j in 0..=(n_v) {
                if j <= k.saturating_sub(p) {
                    new_row.push(row[j.min(n_v.saturating_sub(1))]);
                } else if j <= k {
                    let alpha = if (knots[j + p] - knots[j]).abs() < 1e-14 {
                        0.0
                    } else {
                        (t - knots[j]) / (knots[j + p] - knots[j])
                    };
                    let prev = if j == 0 {
                        [0.0; 3]
                    } else {
                        row[(j - 1).min(n_v - 1)]
                    };
                    let curr = row[j.min(n_v - 1)];
                    new_row.push(add3(scale3(prev, 1.0 - alpha), scale3(curr, alpha)));
                } else {
                    new_row.push(row[(j - 1).min(n_v - 1)]);
                }
            }
            new_net.push(new_row);
        }
        BsplineSurface {
            control_net: new_net,
            u_knots: self.u_knots.clone(),
            v_knots: new_knots,
            degree_u: self.degree_u,
            degree_v: self.degree_v,
        }
    }
    /// Refine all v-direction knots in `new_knots` sequentially.
    fn refine_v(&self, new_knots: &[f64]) -> Self {
        let mut surf = BsplineSurface {
            control_net: self.control_net.clone(),
            u_knots: self.u_knots.clone(),
            v_knots: self.v_knots.clone(),
            degree_u: self.degree_u,
            degree_v: self.degree_v,
        };
        for &t in new_knots {
            surf = surf.insert_v_knot(t);
        }
        surf
    }
}
/// A piecewise cubic Hermite spline interpolating through control points
/// with user-supplied tangent vectors at each knot.
///
/// Guarantees C0 continuity at joints; C1 continuity if tangents are chosen
/// consistently (e.g., Catmull–Rom style).
pub struct HermiteCurve {
    /// Interpolation knots (positions).
    pub points: Vec<[f64; 3]>,
    /// Tangent vectors at each knot, matching `points` in length.
    pub tangents: Vec<[f64; 3]>,
}
impl HermiteCurve {
    /// Construct a Hermite spline from points and tangents.
    ///
    /// # Panics
    /// Panics if `points.len() != tangents.len()` or fewer than 2 points.
    pub fn new(points: Vec<[f64; 3]>, tangents: Vec<[f64; 3]>) -> Self {
        assert!(points.len() >= 2, "HermiteCurve needs at least 2 points");
        assert_eq!(
            points.len(),
            tangents.len(),
            "points and tangents must have the same length"
        );
        HermiteCurve { points, tangents }
    }
    /// Evaluate the spline at a global parameter `t ∈ [0, n_segments]`.
    ///
    /// The integer part selects the segment; the fractional part is the
    /// local Hermite parameter.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let n = self.points.len() - 1;
        let t = t.clamp(0.0, n as f64);
        let seg = (t.floor() as usize).min(n - 1);
        let u = t - seg as f64;
        let p0 = self.points[seg];
        let p1 = self.points[seg + 1];
        let m0 = self.tangents[seg];
        let m1 = self.tangents[seg + 1];
        let h00 = 2.0 * u * u * u - 3.0 * u * u + 1.0;
        let h10 = u * u * u - 2.0 * u * u + u;
        let h01 = -2.0 * u * u * u + 3.0 * u * u;
        let h11 = u * u * u - u * u;
        let mut result = [0.0f64; 3];
        for i in 0..3 {
            result[i] = h00 * p0[i] + h10 * m0[i] + h01 * p1[i] + h11 * m1[i];
        }
        result
    }
    /// Derivative of the spline at global parameter `t`.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let n = self.points.len() - 1;
        let t = t.clamp(0.0, n as f64);
        let seg = (t.floor() as usize).min(n - 1);
        let u = t - seg as f64;
        let p0 = self.points[seg];
        let p1 = self.points[seg + 1];
        let m0 = self.tangents[seg];
        let m1 = self.tangents[seg + 1];
        let dh00 = 6.0 * u * u - 6.0 * u;
        let dh10 = 3.0 * u * u - 4.0 * u + 1.0;
        let dh01 = -6.0 * u * u + 6.0 * u;
        let dh11 = 3.0 * u * u - 2.0 * u;
        let mut result = [0.0f64; 3];
        for i in 0..3 {
            result[i] = dh00 * p0[i] + dh10 * m0[i] + dh01 * p1[i] + dh11 * m1[i];
        }
        result
    }
}
/// A swept surface: a 2-D profile curve (in a local frame) swept along a
/// spine `BezierCurve`. The profile is evaluated at parameter `v` and
/// translated along the spine at parameter `u`.
///
/// The swept point is:
/// `spine(u)  +  profile(v)`
/// (no rotation — this is a translational sweep).
pub struct SweptSurface {
    /// Spine curve parameterised by `u`.
    pub spine: BezierCurve,
    /// Cross-section profile parameterised by `v`.
    pub profile: BezierCurve,
}
impl SweptSurface {
    /// Create a new swept surface.
    pub fn new(spine: BezierCurve, profile: BezierCurve) -> Self {
        SweptSurface { spine, profile }
    }
    /// Evaluate the surface at `(u, v) ∈ [0, 1]²`.
    ///
    /// The result is `spine(u) + (profile(v) - profile(0))`.
    /// This centres the profile at the spine point.
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let sp = self.spine.evaluate(u);
        let prof = self.profile.evaluate(v);
        let prof0 = self.profile.evaluate(0.0);
        add3(sp, sub3(prof, prof0))
    }
    /// Surface normal at `(u, v)` via finite differences.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(0.0);
        let u1 = (u + eps).min(1.0);
        let v0 = (v - eps).max(0.0);
        let v1 = (v + eps).min(1.0);
        let du = scale3(
            sub3(self.evaluate(u1, v), self.evaluate(u0, v)),
            1.0 / (u1 - u0),
        );
        let dv = scale3(
            sub3(self.evaluate(u, v1), self.evaluate(u, v0)),
            1.0 / (v1 - v0),
        );
        normalize3(cross3(du, dv))
    }
}
/// A parametric circular arc in 3-D space.
///
/// The arc lies in a plane defined by two orthogonal unit vectors `u_axis` and
/// `v_axis`. The arc sweeps from `start_angle` to `end_angle` (in radians)
/// around `center` at the given `radius`.
pub struct Arc {
    /// Centre point of the arc.
    pub center: [f64; 3],
    /// Radius of the arc.
    pub radius: f64,
    /// Start angle in radians.
    pub start_angle: f64,
    /// End angle in radians.
    pub end_angle: f64,
    /// First axis of the arc plane (should be a unit vector).
    pub u_axis: [f64; 3],
    /// Second axis of the arc plane (should be a unit vector perpendicular to `u_axis`).
    pub v_axis: [f64; 3],
}
impl Arc {
    /// Construct an arc in the XY plane.
    pub fn in_xy_plane(center: [f64; 3], radius: f64, start_angle: f64, end_angle: f64) -> Self {
        Arc {
            center,
            radius,
            start_angle,
            end_angle,
            u_axis: [1.0, 0.0, 0.0],
            v_axis: [0.0, 1.0, 0.0],
        }
    }
    /// Evaluate the arc at parameter `t ∈ [0, 1]`.
    ///
    /// `t = 0` gives `start_angle`, `t = 1` gives `end_angle`.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let angle = self.start_angle + t * (self.end_angle - self.start_angle);
        let (sin_a, cos_a) = angle.sin_cos();
        add3(
            self.center,
            add3(
                scale3(self.u_axis, self.radius * cos_a),
                scale3(self.v_axis, self.radius * sin_a),
            ),
        )
    }
    /// Tangent vector at parameter `t` (not normalised).
    pub fn tangent(&self, t: f64) -> [f64; 3] {
        let angle = self.start_angle + t * (self.end_angle - self.start_angle);
        let da = self.end_angle - self.start_angle;
        let (sin_a, cos_a) = angle.sin_cos();
        add3(
            scale3(self.u_axis, -self.radius * sin_a * da),
            scale3(self.v_axis, self.radius * cos_a * da),
        )
    }
    /// Arc length (exact: radius × |end_angle − start_angle|).
    pub fn arc_length(&self) -> f64 {
        self.radius * (self.end_angle - self.start_angle).abs()
    }
    /// Sample the arc at `n` uniformly-spaced parameter values `t ∈ [0, 1]`.
    pub fn sample(&self, n: usize) -> Vec<[f64; 3]> {
        if n == 0 {
            return Vec::new();
        }
        (0..n)
            .map(|i| self.evaluate(i as f64 / (n - 1).max(1) as f64))
            .collect()
    }
}
/// A lofted surface interpolated between two Bezier curves.
///
/// For parameter `(u, v)`, the surface point is:
/// `lerp( curve_a.evaluate(u),  curve_b.evaluate(u),  v )`.
pub struct LoftSurface {
    /// The "bottom" profile curve (at v = 0).
    pub curve_a: BezierCurve,
    /// The "top" profile curve (at v = 1).
    pub curve_b: BezierCurve,
}
impl LoftSurface {
    /// Create a new loft surface from two Bezier curves.
    pub fn new(curve_a: BezierCurve, curve_b: BezierCurve) -> Self {
        LoftSurface { curve_a, curve_b }
    }
    /// Evaluate the lofted surface at `(u, v) ∈ [0, 1]²`.
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let pa = self.curve_a.evaluate(u);
        let pb = self.curve_b.evaluate(u);
        lerp3(pa, pb, v)
    }
    /// Surface normal at `(u, v)` via finite differences: ∂P/∂u × ∂P/∂v.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(0.0);
        let u1 = (u + eps).min(1.0);
        let v0 = (v - eps).max(0.0);
        let v1 = (v + eps).min(1.0);
        let du = scale3(
            sub3(self.evaluate(u1, v), self.evaluate(u0, v)),
            1.0 / (u1 - u0),
        );
        let dv = scale3(
            sub3(self.evaluate(u, v1), self.evaluate(u, v0)),
            1.0 / (v1 - v0),
        );
        normalize3(cross3(du, dv))
    }
    /// Sample the surface as a `nu × nv` grid, returned row-major.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.evaluate(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A tensor-product Bezier surface (typically bicubic: rows = cols = 4).
pub struct BezierSurface {
    /// Control grid, indexed as `control_grid[row][col]`.
    pub control_grid: Vec<Vec<[f64; 3]>>,
    /// Number of rows in the control grid.
    pub rows: usize,
    /// Number of columns in the control grid.
    pub cols: usize,
}
impl BezierSurface {
    /// Create a bicubic (4×4) Bezier surface from a control grid.
    pub fn new_bicubic(control_grid: Vec<Vec<[f64; 3]>>) -> Self {
        let rows = control_grid.len();
        assert!(rows > 0, "control grid must have at least one row");
        let cols = control_grid[0].len();
        assert!(cols > 0, "control grid must have at least one column");
        BezierSurface {
            control_grid,
            rows,
            cols,
        }
    }
    /// Evaluate the surface at `(u, v) ∈ [0, 1]²`.
    ///
    /// Apply de Casteljau row-wise (over u) to get `rows` intermediate points,
    /// then apply de Casteljau over those points along v.
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let row_pts: Vec<[f64; 3]> = self
            .control_grid
            .iter()
            .map(|row| {
                let curve = BezierCurve::new(row.clone());
                curve.evaluate(u)
            })
            .collect();
        let col_curve = BezierCurve::new(row_pts);
        col_curve.evaluate(v)
    }
    /// Surface normal at `(u, v)`: ∂P/∂u × ∂P/∂v, normalized.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let du = {
            let u0 = (u - eps).max(0.0);
            let u1 = (u + eps).min(1.0);
            let p0 = self.evaluate(u0, v);
            let p1 = self.evaluate(u1, v);
            scale3(sub3(p1, p0), 1.0 / (u1 - u0))
        };
        let dv = {
            let v0 = (v - eps).max(0.0);
            let v1 = (v + eps).min(1.0);
            let p0 = self.evaluate(u, v0);
            let p1 = self.evaluate(u, v1);
            scale3(sub3(p1, p0), 1.0 / (v1 - v0))
        };
        normalize3(cross3(du, dv))
    }
}
/// A non-uniform B-spline curve defined by a knot vector and control points.
///
/// Uses the Cox–de Boor recursion for evaluation and supports knot insertion.
pub struct BSpline {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
    /// Knot vector (must have length = `control_points.len() + degree + 1`).
    pub knot_vector: Vec<f64>,
    /// Polynomial degree.
    pub degree: usize,
}
impl BSpline {
    /// Construct a new B-spline curve.
    ///
    /// # Panics
    /// Panics if `knot_vector.len() != control_points.len() + degree + 1`.
    pub fn new(degree: usize, control_points: Vec<[f64; 3]>, knot_vector: Vec<f64>) -> Self {
        assert_eq!(
            knot_vector.len(),
            control_points.len() + degree + 1,
            "knot_vector length must equal n_ctrl + degree + 1"
        );
        BSpline {
            control_points,
            knot_vector,
            degree,
        }
    }
    /// Create a clamped uniform B-spline (interpolates first and last control points).
    pub fn clamped_uniform(degree: usize, control_points: Vec<[f64; 3]>) -> Self {
        let n = control_points.len();
        let m = n + degree + 1;
        let interior = m - 2 * (degree + 1);
        let mut knots = Vec::with_capacity(m);
        knots.extend(std::iter::repeat_n(0.0_f64, degree + 1));
        for i in 1..=(interior) {
            knots.push(i as f64 / (interior + 1) as f64);
        }
        knots.extend(std::iter::repeat_n(1.0_f64, degree + 1));
        Self::new(degree, control_points, knots)
    }
    /// Evaluate the B-spline at parameter `t` using the de Boor algorithm.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let p = self.degree;
        let knots = &self.knot_vector;
        let t = t.clamp(
            *knots.first().unwrap_or(&0.0),
            *knots.last().unwrap_or(&1.0),
        );
        let k = self.find_span(t);
        let mut d: Vec<[f64; 3]> = (0..=(p)).map(|j| self.control_points[k - p + j]).collect();
        for r in 1..=(p) {
            for j in (r..=(p)).rev() {
                let i = k - p + j;
                let denom = knots[i + p - r + 1] - knots[i];
                let alpha = if denom.abs() < 1e-12 {
                    0.0
                } else {
                    (t - knots[i]) / denom
                };
                d[j] = lerp3(d[j - 1], d[j], alpha);
            }
        }
        d[p]
    }
    /// Find the knot span index for `t` (index k such that knots\[k\] <= t < knots\[k+1\]).
    fn find_span(&self, t: f64) -> usize {
        let n = self.control_points.len();
        let p = self.degree;
        let knots = &self.knot_vector;
        if (t - knots[n]).abs() < 1e-12 {
            let mut k = n - 1;
            while k > p && (knots[k] - knots[n]).abs() < 1e-12 {
                k -= 1;
            }
            return k;
        }
        let mut lo = p;
        let mut hi = n;
        let mut mid = (lo + hi) / 2;
        while t < knots[mid] || t >= knots[mid + 1] {
            if t < knots[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
            mid = (lo + hi) / 2;
            if lo + 1 >= hi {
                break;
            }
        }
        mid
    }
    /// Insert knot `t_new` into the knot vector (Boehm's algorithm).
    ///
    /// Returns a new `BSpline` with one additional control point and knot.
    pub fn insert_knot(&self, t_new: f64) -> Self {
        let p = self.degree;
        let knots = &self.knot_vector;
        let pts = &self.control_points;
        let n = pts.len();
        let k = self.find_span(t_new);
        let mut new_knots = Vec::with_capacity(knots.len() + 1);
        new_knots.extend_from_slice(&knots[..=k]);
        new_knots.push(t_new);
        new_knots.extend_from_slice(&knots[k + 1..]);
        let mut new_pts: Vec<[f64; 3]> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            if i <= k - p {
                new_pts.push(pts[i]);
            } else if i <= k {
                let alpha_denom = knots[i + p] - knots[i];
                let alpha = if alpha_denom.abs() < 1e-12 {
                    1.0
                } else {
                    (t_new - knots[i]) / alpha_denom
                };
                let prev = if i == 0 { [0.0; 3] } else { pts[i - 1] };
                let curr = pts[i];
                new_pts.push(lerp3(prev, curr, alpha));
            } else {
                new_pts.push(pts[i - 1]);
            }
        }
        BSpline {
            control_points: new_pts,
            knot_vector: new_knots,
            degree: p,
        }
    }
}
impl BSpline {
    /// Insert knot `t` in place, modifying this B-spline (Boehm's algorithm).
    pub fn insert_knot_inplace(&mut self, t: f64) {
        let new_spline = self.insert_knot(t);
        self.control_points = new_spline.control_points;
        self.knot_vector = new_spline.knot_vector;
    }
}
impl BSpline {
    /// Derivative of the B-spline at parameter `t` using finite differences.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let knots = &self.knot_vector;
        let t_min = *knots.first().unwrap_or(&0.0);
        let t_max = *knots.last().unwrap_or(&1.0);
        let t0 = (t - eps).max(t_min);
        let t1 = (t + eps).min(t_max);
        let p0 = self.eval(t0);
        let p1 = self.eval(t1);
        scale3(sub3(p1, p0), 1.0 / (t1 - t0))
    }
    /// Arc length approximation by sampling `n_samples` parameter values.
    pub fn arc_length(&self, n_samples: usize) -> f64 {
        assert!(n_samples >= 2, "arc_length needs at least 2 samples");
        let knots = &self.knot_vector;
        let t_min = *knots.first().unwrap_or(&0.0);
        let t_max = *knots.last().unwrap_or(&1.0);
        let mut len = 0.0_f64;
        let mut prev = self.eval(t_min);
        for k in 1..=n_samples {
            let t = t_min + (t_max - t_min) * k as f64 / n_samples as f64;
            let curr = self.eval(t);
            len += dist3(prev, curr);
            prev = curr;
        }
        len
    }
}
/// A Bezier curve of degree `n = control_points.len() - 1`.
pub struct BezierCurve {
    /// Control points of the curve.
    pub control_points: Vec<[f64; 3]>,
}
impl BezierCurve {
    /// Create a new Bezier curve from a list of control points.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        assert!(
            !points.is_empty(),
            "BezierCurve needs at least one control point"
        );
        BezierCurve {
            control_points: points,
        }
    }
    /// Degree of the curve (= number of control points − 1).
    pub fn degree(&self) -> usize {
        self.control_points.len() - 1
    }
    /// Evaluate the curve at parameter `t ∈ [0, 1]` using de Casteljau's algorithm.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let mut pts = self.control_points.clone();
        let n = pts.len();
        for r in 1..n {
            for i in 0..(n - r) {
                pts[i] = lerp3(pts[i], pts[i + 1], t);
            }
        }
        pts[0]
    }
    /// Tangent vector at `t` — the derivative Bezier of degree `n-1`.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let n = self.degree();
        if n == 0 {
            return [0.0, 0.0, 0.0];
        }
        let deriv_pts: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                scale3(
                    sub3(self.control_points[i + 1], self.control_points[i]),
                    n as f64,
                )
            })
            .collect();
        let deriv_curve = BezierCurve::new(deriv_pts);
        deriv_curve.evaluate(t)
    }
    /// Approximate arc length by sampling `n_samples` equidistant parameter values.
    pub fn arc_length(&self, n_samples: usize) -> f64 {
        assert!(n_samples >= 2, "arc_length needs at least 2 samples");
        let mut length = 0.0;
        let mut prev = self.evaluate(0.0);
        for k in 1..=n_samples {
            let t = k as f64 / n_samples as f64;
            let curr = self.evaluate(t);
            length += dist3(prev, curr);
            prev = curr;
        }
        length
    }
}
impl BezierCurve {
    /// Alias for `evaluate(t)` — de Casteljau evaluation at `t ∈ [0, 1]`.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        self.evaluate(t)
    }
}
/// The Frenet–Serret frame at a point on a curve.
#[derive(Debug, Clone)]
pub struct FrenetFrame {
    /// Unit tangent vector T = C'(t) / |C'(t)|.
    pub tangent: [f64; 3],
    /// Unit normal vector N = T'(t) / |T'(t)|.
    pub normal: [f64; 3],
    /// Binormal vector B = T × N.
    pub binormal: [f64; 3],
}
impl FrenetFrame {
    /// Compute the Frenet frame at parameter `t` for a curve whose tangent and
    /// second derivative are provided as closures.
    ///
    /// `tangent_fn(t)` should return `C'(t)` and `dtangent_fn(t)` should
    /// return `C''(t)`.  Both are then normalized internally.
    pub fn compute(
        tangent_fn: impl Fn(f64) -> [f64; 3],
        dtangent_fn: impl Fn(f64) -> [f64; 3],
        t: f64,
    ) -> Self {
        let tan_raw = tangent_fn(t);
        let tangent = normalize3(tan_raw);
        let dtan = dtangent_fn(t);
        let proj = dot3(dtan, tangent);
        let normal_raw = sub3(dtan, scale3(tangent, proj));
        let normal = if dot3(normal_raw, normal_raw).sqrt() < 1e-12 {
            arbitrary_perp(tangent)
        } else {
            normalize3(normal_raw)
        };
        let binormal = normalize3(cross3(tangent, normal));
        FrenetFrame {
            tangent,
            normal,
            binormal,
        }
    }
    /// Compute the Frenet frame for a `BezierCurve` at `t` using finite
    /// differences for the second derivative.
    pub fn from_bezier(curve: &BezierCurve, t: f64) -> Self {
        let eps = 1e-5_f64;
        let tangent_fn = |u: f64| curve.derivative(u);
        let dtangent_fn = |u: f64| {
            let t0 = (u - eps).max(0.0);
            let t1 = (u + eps).min(1.0);
            let d = sub3(curve.derivative(t1), curve.derivative(t0));
            scale3(d, 1.0 / (t1 - t0))
        };
        Self::compute(tangent_fn, dtangent_fn, t)
    }
}
/// A surface of revolution: a profile curve in the XZ plane rotated around
/// the Z axis.
///
/// `u ∈ [0, 1]` sweeps the full 2π revolution (or a fraction thereof),
/// `v ∈ [0, 1]` traverses the profile curve.
pub struct RevolutionSurface {
    /// Profile curve (in XZ plane: profile\[.\]\[0\] = X offset, profile\[.\]\[2\] = Z height).
    pub profile: BezierCurve,
    /// Total sweep angle in radians (default 2π for a full revolution).
    pub sweep_angle: f64,
}
impl RevolutionSurface {
    /// Create a surface of revolution with a full 2π sweep.
    pub fn new(profile: BezierCurve) -> Self {
        RevolutionSurface {
            profile,
            sweep_angle: std::f64::consts::TAU,
        }
    }
    /// Create a surface of revolution with an explicit sweep angle.
    pub fn with_angle(profile: BezierCurve, sweep_angle: f64) -> Self {
        RevolutionSurface {
            profile,
            sweep_angle,
        }
    }
    /// Evaluate at `(u, v)`:
    /// - `v` selects a point on the profile `(r, 0, z)`,
    /// - `u` rotates it by `u * sweep_angle` around the Z axis.
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let p = self.profile.evaluate(v);
        let r = p[0];
        let z = p[2];
        let angle = u * self.sweep_angle;
        let (sin_a, cos_a) = angle.sin_cos();
        [r * cos_a, r * sin_a, z]
    }
    /// Surface normal at `(u, v)` via finite differences.
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let eps = 1e-5_f64;
        let u0 = (u - eps).max(0.0);
        let u1 = (u + eps).min(1.0);
        let v0 = (v - eps).max(0.0);
        let v1 = (v + eps).min(1.0);
        let du = scale3(
            sub3(self.evaluate(u1, v), self.evaluate(u0, v)),
            1.0 / (u1 - u0),
        );
        let dv = scale3(
            sub3(self.evaluate(u, v1), self.evaluate(u, v0)),
            1.0 / (v1 - v0),
        );
        normalize3(cross3(du, dv))
    }
    /// Sample the surface as a `nu × nv` grid.
    pub fn sample_grid(&self, nu: usize, nv: usize) -> Vec<Vec<[f64; 3]>> {
        (0..nu)
            .map(|i| {
                let u = i as f64 / (nu - 1).max(1) as f64;
                (0..nv)
                    .map(|j| {
                        let v = j as f64 / (nv - 1).max(1) as f64;
                        self.evaluate(u, v)
                    })
                    .collect()
            })
            .collect()
    }
}
/// A Non-Uniform Rational B-Spline (NURBS) curve.
pub struct NurbsCurve {
    /// Control points.
    pub control_points: Vec<[f64; 3]>,
    /// Rational weights, one per control point.
    pub weights: Vec<f64>,
    /// Knot vector.
    pub knot_vector: Vec<f64>,
    /// Degree of the B-spline basis.
    pub degree: usize,
}
impl NurbsCurve {
    /// Construct a NURBS curve.
    pub fn new(
        degree: usize,
        control_points: Vec<[f64; 3]>,
        weights: Vec<f64>,
        knot_vector: Vec<f64>,
    ) -> Self {
        assert_eq!(
            control_points.len(),
            weights.len(),
            "control_points and weights must have the same length"
        );
        NurbsCurve {
            control_points,
            weights,
            knot_vector,
            degree,
        }
    }
    /// Cox–de Boor B-spline basis function N_{i,p}(t).
    pub fn b_spline_basis(i: usize, p: usize, t: f64, knots: &[f64]) -> f64 {
        if p == 0 {
            let last = *knots.last().unwrap_or(&0.0);
            if knots[i] <= t && t < knots[i + 1] {
                return 1.0;
            }
            if (t - last).abs() < 1e-12 && (knots[i + 1] - last).abs() < 1e-12 {
                return 1.0;
            }
            return 0.0;
        }
        let denom1 = knots[i + p] - knots[i];
        let denom2 = knots[i + p + 1] - knots[i + 1];
        let term1 = if denom1.abs() < 1e-12 {
            0.0
        } else {
            (t - knots[i]) / denom1 * Self::b_spline_basis(i, p - 1, t, knots)
        };
        let term2 = if denom2.abs() < 1e-12 {
            0.0
        } else {
            (knots[i + p + 1] - t) / denom2 * Self::b_spline_basis(i + 1, p - 1, t, knots)
        };
        term1 + term2
    }
    /// Find the knot span index for parameter `t` using binary search.
    ///
    /// Returns index `k` such that `knots[k] <= t < knots[k+1]`.
    pub fn find_span(n: usize, p: usize, t: f64, knots: &[f64]) -> usize {
        if (t - knots[n + 1]).abs() < 1e-12 {
            return n;
        }
        let mut lo = p;
        let mut hi = n + 1;
        let mut mid = (lo + hi) / 2;
        while t < knots[mid] || t >= knots[mid + 1] {
            if t < knots[mid] {
                hi = mid;
            } else {
                lo = mid;
            }
            mid = (lo + hi) / 2;
        }
        mid
    }
    /// Evaluate the NURBS curve at parameter `t`.
    ///
    /// Uses the rational formula: `Σ(w_i · N_{i,p}(t) · P_i) / Σ(w_i · N_{i,p}(t))`.
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let p = self.degree;
        let knots = &self.knot_vector;
        let mut numerator = [0.0_f64; 3];
        let mut denominator = 0.0_f64;
        for i in 0..n {
            let basis = Self::b_spline_basis(i, p, t, knots);
            let w = self.weights[i];
            let wn = w * basis;
            numerator[0] += wn * self.control_points[i][0];
            numerator[1] += wn * self.control_points[i][1];
            numerator[2] += wn * self.control_points[i][2];
            denominator += wn;
        }
        if denominator.abs() < 1e-12 {
            [0.0, 0.0, 0.0]
        } else {
            scale3(numerator, 1.0 / denominator)
        }
    }
}
impl NurbsCurve {
    /// Compute the arc length of the NURBS curve over `[t_start, t_end]`
    /// using `n_samples` trapezoid-rule segments.
    pub fn arc_length(&self, t_start: f64, t_end: f64, n_samples: usize) -> f64 {
        let n = n_samples.max(2);
        let mut total = 0.0_f64;
        let mut prev = self.evaluate(t_start);
        for k in 1..=n {
            let t = t_start + (t_end - t_start) * k as f64 / n as f64;
            let curr = self.evaluate(t);
            total += dist3(prev, curr);
            prev = curr;
        }
        total
    }
    /// Build an arc-length reparametrization table with `n_samples` samples.
    ///
    /// Returns a pair of `Vec`f64` — `(arc_lengths, params)` — where
    /// `arc_lengths\[i\]` is the arc length at parameter `params\[i\]`.  The
    /// parameter range is taken from the knot vector.
    pub fn arc_length_table(&self, n_samples: usize) -> (Vec<f64>, Vec<f64>) {
        let n = n_samples.max(2);
        let t_min = *self.knot_vector.first().unwrap_or(&0.0);
        let t_max = *self.knot_vector.last().unwrap_or(&1.0);
        let mut params = Vec::with_capacity(n + 1);
        let mut arc_lengths = Vec::with_capacity(n + 1);
        let mut total = 0.0_f64;
        let mut prev = self.evaluate(t_min);
        params.push(t_min);
        arc_lengths.push(0.0);
        for k in 1..=n {
            let t = t_min + (t_max - t_min) * k as f64 / n as f64;
            let curr = self.evaluate(t);
            total += dist3(prev, curr);
            prev = curr;
            params.push(t);
            arc_lengths.push(total);
        }
        (arc_lengths, params)
    }
    /// Invert the arc-length table: find the parameter `t` corresponding to
    /// a given arc length `s` by linear interpolation.
    fn invert_arc_length_table(s: f64, arc_lengths: &[f64], params: &[f64]) -> f64 {
        if s <= arc_lengths[0] {
            return params[0];
        }
        let last = *arc_lengths.last().expect("collection should not be empty");
        if s >= last {
            return *params.last().expect("collection should not be empty");
        }
        let mut lo = 0_usize;
        let mut hi = arc_lengths.len() - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if arc_lengths[mid] <= s {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t_lo = params[lo];
        let t_hi = params[hi];
        let s_lo = arc_lengths[lo];
        let s_hi = arc_lengths[hi];
        let frac = if (s_hi - s_lo).abs() < 1e-14 {
            0.0
        } else {
            (s - s_lo) / (s_hi - s_lo)
        };
        t_lo + frac * (t_hi - t_lo)
    }
    /// Reparametrize by arc length: sample the curve at `n_out` equally-spaced
    /// arc-length values and return the corresponding 3-D positions.
    ///
    /// Internally builds an arc-length table with `n_table` entries.
    pub fn compute_arc_length_reparametrize(&self, n_table: usize, n_out: usize) -> Vec<[f64; 3]> {
        let (arc_lengths, params) = self.arc_length_table(n_table);
        let total_len = *arc_lengths.last().unwrap_or(&0.0);
        if total_len < 1e-14 || n_out == 0 {
            return Vec::new();
        }
        (0..n_out)
            .map(|k| {
                let s = total_len * k as f64 / (n_out - 1).max(1) as f64;
                let t = Self::invert_arc_length_table(s, &arc_lengths, &params);
                self.evaluate(t)
            })
            .collect()
    }
}
/// A Catmull–Rom spline with centripetal parameterization (α = 0.5).
///
/// Interpolates through all control points.  Tangents are computed
/// automatically using the Catmull–Rom formula so the spline is C1.
pub struct CatmullRomSpline {
    /// Interpolation points.
    pub points: Vec<[f64; 3]>,
    /// Knot sequence for centripetal parameterization (α = 0.5).
    pub knots: Vec<f64>,
}
impl CatmullRomSpline {
    /// Construct a centripetal Catmull–Rom spline from the given points.
    ///
    /// # Panics
    /// Panics if fewer than 2 points are supplied.
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        assert!(
            points.len() >= 2,
            "CatmullRomSpline needs at least 2 points"
        );
        let knots = Self::centripetal_knots(&points);
        CatmullRomSpline { points, knots }
    }
    /// Centripetal knot sequence: t_0 = 0, t_{i+1} = t_i + |P_{i+1} - P_i|^0.5.
    fn centripetal_knots(pts: &[[f64; 3]]) -> Vec<f64> {
        let mut t = vec![0.0f64];
        for i in 1..pts.len() {
            let d = dist3(pts[i], pts[i - 1]);
            t.push(t[i - 1] + d.sqrt().max(1e-12));
        }
        let last = *t.last().expect("collection should not be empty");
        if last > 1e-12 {
            t.iter_mut().for_each(|v| *v /= last);
        }
        t
    }
    /// Evaluate the spline at global parameter `t ∈ \[0, 1\]`.
    pub fn eval(&self, t: f64) -> [f64; 3] {
        let pts = &self.points;
        let knots = &self.knots;
        let n = pts.len();
        let t = t.clamp(0.0, 1.0);
        let mut seg = 0usize;
        for i in 0..n - 1 {
            if t <= knots[i + 1] + 1e-12 {
                seg = i;
                break;
            }
            seg = n - 2;
        }
        let i0 = if seg == 0 { 0 } else { seg - 1 };
        let i1 = seg;
        let i2 = (seg + 1).min(n - 1);
        let i3 = (seg + 2).min(n - 1);
        let t0 = knots[i0];
        let t1 = knots[i1];
        let t2 = knots[i2];
        let t3 = knots[i3];
        Self::barry_goldman(pts[i0], pts[i1], pts[i2], pts[i3], t0, t1, t2, t3, t)
    }
    /// Barry–Goldman non-uniform Catmull–Rom evaluation.
    fn barry_goldman(
        p0: [f64; 3],
        p1: [f64; 3],
        p2: [f64; 3],
        p3: [f64; 3],
        t0: f64,
        t1: f64,
        t2: f64,
        t3: f64,
        t: f64,
    ) -> [f64; 3] {
        fn interp(a: [f64; 3], b: [f64; 3], ta: f64, tb: f64, t: f64) -> [f64; 3] {
            let span = tb - ta;
            if span.abs() < 1e-12 {
                return a;
            }
            lerp3(a, b, (t - ta) / span)
        }
        let a1 = interp(p0, p1, t0, t1, t);
        let a2 = interp(p1, p2, t1, t2, t);
        let a3 = interp(p2, p3, t2, t3, t);
        let b1 = interp(a1, a2, t0, t2, t);
        let b2 = interp(a2, a3, t1, t3, t);
        interp(b1, b2, t1, t2, t)
    }
}
