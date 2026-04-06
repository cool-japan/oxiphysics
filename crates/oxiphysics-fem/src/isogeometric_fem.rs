#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Isogeometric Analysis (IGA) with NURBS-based FEM.
//!
//! This module provides a complete isogeometric analysis framework:
//! - **B-spline basis** via Cox-de Boor recursion with derivative evaluation
//! - **NURBS curves** with weighted control points and knot vectors
//! - **NURBS surfaces** (bivariate tensor-product NURBS)
//! - **IGA elements** with NURBS shape functions, Jacobian, Gauss integration
//! - **IGA assembler** for global stiffness matrix K and load vector f
//! - **Boundary conditions** for IGA (Dirichlet and Neumann)
//! - **Trimmed NURBS** with trimming curve support

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

#[allow(unused_imports)]
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Dot product of two slices.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm of a slice.
#[inline]
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// 2D cross product (scalar result).
#[inline]
fn cross2d(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// 3D cross product.
#[inline]
fn cross3d(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Distance between two 3D points.
#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ---------------------------------------------------------------------------
// Gauss quadrature
// ---------------------------------------------------------------------------

/// Gauss-Legendre quadrature points and weights on \[-1, 1\].
///
/// Returns `(points, weights)` for the given number of points (1..=5).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    match n {
        1 => (vec![0.0], vec![2.0]),
        2 => {
            let x = 1.0 / 3.0_f64.sqrt();
            (vec![-x, x], vec![1.0, 1.0])
        }
        3 => {
            let x = (3.0 / 5.0_f64).sqrt();
            (vec![-x, 0.0, x], vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
        }
        4 => {
            let a = (3.0 / 7.0 - 2.0 / 7.0 * (6.0 / 5.0_f64).sqrt()).sqrt();
            let b = (3.0 / 7.0 + 2.0 / 7.0 * (6.0 / 5.0_f64).sqrt()).sqrt();
            let wa = (18.0 + 30.0_f64.sqrt()) / 36.0;
            let wb = (18.0 - 30.0_f64.sqrt()) / 36.0;
            (vec![-b, -a, a, b], vec![wb, wa, wa, wb])
        }
        _ => {
            // 5-point
            let x1 = 0.0;
            let x2 = (5.0 - 2.0 * (10.0 / 7.0_f64).sqrt()).sqrt() / 3.0;
            let x3 = (5.0 + 2.0 * (10.0 / 7.0_f64).sqrt()).sqrt() / 3.0;
            let w1 = 128.0 / 225.0;
            let w2 = (322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0;
            let w3 = (322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0;
            (vec![-x3, -x2, x1, x2, x3], vec![w3, w2, w1, w2, w3])
        }
    }
}

/// Map Gauss points from \[-1,1\] to \[a,b\].
pub fn gauss_map(points: &[f64], weights: &[f64], a: f64, b: f64) -> (Vec<f64>, Vec<f64>) {
    let half_len = (b - a) / 2.0;
    let mid = (a + b) / 2.0;
    let pts: Vec<f64> = points.iter().map(|xi| mid + half_len * xi).collect();
    let wts: Vec<f64> = weights.iter().map(|w| w * half_len).collect();
    (pts, wts)
}

// ---------------------------------------------------------------------------
// B-spline Basis
// ---------------------------------------------------------------------------

/// B-spline basis function evaluator using Cox-de Boor recursion.
///
/// Supports arbitrary degree and open/periodic knot vectors.
#[derive(Debug, Clone)]
pub struct BSplineBasis {
    /// Knot vector (non-decreasing sequence).
    pub knots: Vec<f64>,
    /// Polynomial degree.
    pub degree: usize,
}

impl BSplineBasis {
    /// Create a new B-spline basis with given knot vector and degree.
    pub fn new(knots: Vec<f64>, degree: usize) -> Self {
        Self { knots, degree }
    }

    /// Create a uniform open (clamped) knot vector for `n` basis functions of degree `p`.
    #[allow(clippy::same_item_push)]
    pub fn uniform_open(n: usize, p: usize) -> Self {
        let m = n + p + 1;
        let mut knots = Vec::with_capacity(m);
        for _ in 0..=p {
            knots.push(0.0);
        }
        let interior = n - p;
        for i in 1..interior {
            knots.push(i as f64 / interior as f64);
        }
        for _ in 0..=p {
            knots.push(1.0);
        }
        Self { knots, degree: p }
    }

    /// Number of basis functions.
    pub fn num_basis(&self) -> usize {
        self.knots.len() - self.degree - 1
    }

    /// Evaluate basis function N_{i,p}(t) using Cox-de Boor recursion.
    pub fn evaluate(&self, i: usize, t: f64) -> f64 {
        self.cox_de_boor(i, self.degree, t)
    }

    /// Cox-de Boor recursion for N_{i,p}(t).
    fn cox_de_boor(&self, i: usize, p: usize, t: f64) -> f64 {
        if p == 0 {
            // N_{i,0}(t) = 1 if knot[i] <= t < knot[i+1], else 0
            // Special case for the last interval: include the right endpoint
            if i + 1 < self.knots.len() {
                let left = self.knots[i];
                let right = self.knots[i + 1];
                if (left <= t && t < right)
                    || (t == right && (right - self.knots[self.knots.len() - 1]).abs() < 1e-14)
                {
                    return 1.0;
                }
            }
            return 0.0;
        }

        let mut result = 0.0;
        let denom1 = self.knots[i + p] - self.knots[i];
        if denom1.abs() > 1e-14 {
            result += (t - self.knots[i]) / denom1 * self.cox_de_boor(i, p - 1, t);
        }
        let denom2 = self.knots[i + p + 1] - self.knots[i + 1];
        if denom2.abs() > 1e-14 {
            result += (self.knots[i + p + 1] - t) / denom2 * self.cox_de_boor(i + 1, p - 1, t);
        }
        result
    }

    /// Evaluate all non-zero basis functions at parameter t.
    ///
    /// Returns vector of (index, value) pairs.
    pub fn evaluate_all(&self, t: f64) -> Vec<(usize, f64)> {
        let n = self.num_basis();
        let mut result = Vec::new();
        for i in 0..n {
            let val = self.evaluate(i, t);
            if val.abs() > 1e-15 {
                result.push((i, val));
            }
        }
        result
    }

    /// Evaluate all basis functions at parameter t as a dense vector.
    pub fn evaluate_all_dense(&self, t: f64) -> Vec<f64> {
        let n = self.num_basis();
        (0..n).map(|i| self.evaluate(i, t)).collect()
    }

    /// Derivative of basis function N_{i,p}(t).
    ///
    /// Uses the recursive derivative formula:
    /// N'_{i,p}(t) = p / (knot\[i+p\] - knot\[i\]) * N_{i,p-1}(t)
    ///             - p / (knot\[i+p+1\] - knot\[i+1\]) * N_{i+1,p-1}(t)
    pub fn derivative(&self, i: usize, t: f64) -> f64 {
        if self.degree == 0 {
            return 0.0;
        }
        let p = self.degree;
        let mut result = 0.0;

        let denom1 = self.knots[i + p] - self.knots[i];
        if denom1.abs() > 1e-14 {
            result += p as f64 / denom1 * self.cox_de_boor(i, p - 1, t);
        }
        let denom2 = self.knots[i + p + 1] - self.knots[i + 1];
        if denom2.abs() > 1e-14 {
            result -= p as f64 / denom2 * self.cox_de_boor(i + 1, p - 1, t);
        }
        result
    }

    /// Evaluate all basis function derivatives at parameter t.
    pub fn derivative_all(&self, t: f64) -> Vec<f64> {
        let n = self.num_basis();
        (0..n).map(|i| self.derivative(i, t)).collect()
    }

    /// Second derivative of basis function N_{i,p}(t).
    pub fn second_derivative(&self, i: usize, t: f64) -> f64 {
        if self.degree < 2 {
            return 0.0;
        }
        let p = self.degree;

        // Build a lower-degree basis for computing derivatives of (p-1) degree functions.
        let lower = BSplineBasis {
            knots: self.knots.clone(),
            degree: p - 1,
        };

        let mut result = 0.0;
        let denom1 = self.knots[i + p] - self.knots[i];
        if denom1.abs() > 1e-14 {
            result += p as f64 / denom1 * lower.derivative(i, t);
        }
        let denom2 = self.knots[i + p + 1] - self.knots[i + 1];
        if denom2.abs() > 1e-14 {
            result -= p as f64 / denom2 * lower.derivative(i + 1, t);
        }
        result
    }

    /// Find the knot span index for parameter t.
    ///
    /// Returns the index `i` such that `knots[i] <= t < knots[i+1]`.
    pub fn find_span(&self, t: f64) -> usize {
        let n = self.num_basis();
        if t >= self.knots[n] {
            return n - 1;
        }
        if t <= self.knots[self.degree] {
            return self.degree;
        }
        // Binary search
        let mut low = self.degree;
        let mut high = n;
        let mut mid = (low + high) / 2;
        while t < self.knots[mid] || t >= self.knots[mid + 1] {
            if t < self.knots[mid] {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low + high) / 2;
        }
        mid
    }

    /// Greville abscissae (average knot positions for each basis function).
    pub fn greville_abscissae(&self) -> Vec<f64> {
        let n = self.num_basis();
        let p = self.degree;
        (0..n)
            .map(|i| {
                let sum: f64 = (1..=p).map(|j| self.knots[i + j]).sum();
                sum / p as f64
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// NURBS Curve
// ---------------------------------------------------------------------------

/// A single NURBS control point with weight.
#[derive(Debug, Clone, Copy)]
pub struct WeightedPoint {
    /// Position in 3D.
    pub pos: [f64; 3],
    /// NURBS weight (must be positive).
    pub weight: f64,
}

impl WeightedPoint {
    /// Create a new weighted control point.
    pub fn new(pos: [f64; 3], weight: f64) -> Self {
        Self { pos, weight }
    }

    /// Create a unit-weight control point.
    pub fn unweighted(pos: [f64; 3]) -> Self {
        Self { pos, weight: 1.0 }
    }
}

/// A NURBS (Non-Uniform Rational B-Spline) curve.
///
/// Defined by control points with weights and a B-spline basis.
#[derive(Debug, Clone)]
pub struct NurbsCurve {
    /// B-spline basis.
    pub basis: BSplineBasis,
    /// Control points with weights.
    pub control_points: Vec<WeightedPoint>,
}

impl NurbsCurve {
    /// Create a new NURBS curve.
    pub fn new(basis: BSplineBasis, control_points: Vec<WeightedPoint>) -> Self {
        assert_eq!(
            basis.num_basis(),
            control_points.len(),
            "Number of basis functions must equal number of control points"
        );
        Self {
            basis,
            control_points,
        }
    }

    /// Evaluate the curve at parameter t.
    ///
    /// Uses the rational form: C(t) = sum(N_i * w_i * P_i) / sum(N_i * w_i)
    pub fn evaluate(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let mut numerator = [0.0; 3];
        let mut denominator = 0.0;

        for i in 0..n {
            let basis_val = self.basis.evaluate(i, t);
            let w = self.control_points[i].weight;
            let nw = basis_val * w;
            denominator += nw;
            for d in 0..3 {
                numerator[d] += nw * self.control_points[i].pos[d];
            }
        }

        if denominator.abs() < 1e-14 {
            return [0.0; 3];
        }
        [
            numerator[0] / denominator,
            numerator[1] / denominator,
            numerator[2] / denominator,
        ]
    }

    /// Evaluate the NURBS basis (rational) functions at parameter t.
    ///
    /// R_i(t) = N_i(t) * w_i / sum(N_j(t) * w_j)
    pub fn rational_basis(&self, t: f64) -> Vec<f64> {
        let n = self.control_points.len();
        let mut nw: Vec<f64> = Vec::with_capacity(n);
        let mut sum_nw = 0.0;

        for i in 0..n {
            let val = self.basis.evaluate(i, t) * self.control_points[i].weight;
            nw.push(val);
            sum_nw += val;
        }

        if sum_nw.abs() < 1e-14 {
            return vec![0.0; n];
        }
        nw.iter().map(|v| v / sum_nw).collect()
    }

    /// Evaluate the derivative of the curve at parameter t.
    ///
    /// Uses the quotient rule for the rational form.
    pub fn derivative(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len();
        let mut a = [0.0; 3]; // numerator: sum(N_i * w_i * P_i)
        let mut a_prime = [0.0; 3]; // derivative of numerator
        let mut w_sum = 0.0; // denominator: sum(N_i * w_i)
        let mut w_prime = 0.0; // derivative of denominator

        for i in 0..n {
            let ni = self.basis.evaluate(i, t);
            let ni_prime = self.basis.derivative(i, t);
            let wi = self.control_points[i].weight;
            let pi = self.control_points[i].pos;

            w_sum += ni * wi;
            w_prime += ni_prime * wi;
            for d in 0..3 {
                a[d] += ni * wi * pi[d];
                a_prime[d] += ni_prime * wi * pi[d];
            }
        }

        if w_sum.abs() < 1e-14 {
            return [0.0; 3];
        }

        // Quotient rule: C'(t) = (A'(t) * W(t) - A(t) * W'(t)) / W(t)^2
        let w2 = w_sum * w_sum;
        [
            (a_prime[0] * w_sum - a[0] * w_prime) / w2,
            (a_prime[1] * w_sum - a[1] * w_prime) / w2,
            (a_prime[2] * w_sum - a[2] * w_prime) / w2,
        ]
    }

    /// Compute the arc length of the curve via Gauss quadrature.
    pub fn arc_length(&self, num_gauss: usize) -> f64 {
        let (pts, wts) = gauss_legendre(num_gauss);
        let t_min = self.basis.knots[self.basis.degree];
        let t_max = self.basis.knots[self.basis.knots.len() - self.basis.degree - 1];
        let (mapped_pts, mapped_wts) = gauss_map(&pts, &wts, t_min, t_max);

        let mut length = 0.0;
        for (pt, wt) in mapped_pts.iter().zip(mapped_wts.iter()) {
            let d = self.derivative(*pt);
            let speed = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            length += speed * wt;
        }
        length
    }

    /// Create a NURBS circle of given radius in the xy-plane.
    pub fn circle(radius: f64) -> Self {
        let w = 1.0 / 2.0_f64.sqrt();
        let r = radius;
        let control_points = vec![
            WeightedPoint::new([r, 0.0, 0.0], 1.0),
            WeightedPoint::new([r, r, 0.0], w),
            WeightedPoint::new([0.0, r, 0.0], 1.0),
            WeightedPoint::new([-r, r, 0.0], w),
            WeightedPoint::new([-r, 0.0, 0.0], 1.0),
            WeightedPoint::new([-r, -r, 0.0], w),
            WeightedPoint::new([0.0, -r, 0.0], 1.0),
            WeightedPoint::new([r, -r, 0.0], w),
            WeightedPoint::new([r, 0.0, 0.0], 1.0),
        ];
        let knots = vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ];
        let basis = BSplineBasis::new(knots, 2);
        Self::new(basis, control_points)
    }
}

// ---------------------------------------------------------------------------
// NURBS Surface
// ---------------------------------------------------------------------------

/// A bivariate NURBS surface (tensor-product).
///
/// Defined by a 2D grid of weighted control points and two B-spline bases.
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// B-spline basis in the u-direction.
    pub basis_u: BSplineBasis,
    /// B-spline basis in the v-direction.
    pub basis_v: BSplineBasis,
    /// Control points stored in row-major order (u varies fastest).
    /// Total count = num_u * num_v.
    pub control_points: Vec<WeightedPoint>,
    /// Number of control points in u-direction.
    pub num_u: usize,
    /// Number of control points in v-direction.
    pub num_v: usize,
}

impl NurbsSurface {
    /// Create a new NURBS surface.
    pub fn new(
        basis_u: BSplineBasis,
        basis_v: BSplineBasis,
        control_points: Vec<WeightedPoint>,
        num_u: usize,
        num_v: usize,
    ) -> Self {
        assert_eq!(basis_u.num_basis(), num_u);
        assert_eq!(basis_v.num_basis(), num_v);
        assert_eq!(control_points.len(), num_u * num_v);
        Self {
            basis_u,
            basis_v,
            control_points,
            num_u,
            num_v,
        }
    }

    /// Get control point at (i, j).
    pub fn control_point(&self, i: usize, j: usize) -> &WeightedPoint {
        &self.control_points[j * self.num_u + i]
    }

    /// Evaluate the surface at parameters (u, v).
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let mut numerator = [0.0; 3];
        let mut denominator = 0.0;

        for j in 0..self.num_v {
            let nv = self.basis_v.evaluate(j, v);
            for i in 0..self.num_u {
                let nu = self.basis_u.evaluate(i, u);
                let cp = self.control_point(i, j);
                let nw = nu * nv * cp.weight;
                denominator += nw;
                for d in 0..3 {
                    numerator[d] += nw * cp.pos[d];
                }
            }
        }

        if denominator.abs() < 1e-14 {
            return [0.0; 3];
        }
        [
            numerator[0] / denominator,
            numerator[1] / denominator,
            numerator[2] / denominator,
        ]
    }

    /// Partial derivative with respect to u at (u, v).
    pub fn deriv_u(&self, u: f64, v: f64) -> [f64; 3] {
        let mut a = [0.0; 3];
        let mut a_u = [0.0; 3];
        let mut w_sum = 0.0;
        let mut w_u = 0.0;

        for j in 0..self.num_v {
            let nv = self.basis_v.evaluate(j, v);
            for i in 0..self.num_u {
                let nu = self.basis_u.evaluate(i, u);
                let dnu = self.basis_u.derivative(i, u);
                let cp = self.control_point(i, j);
                let wi = cp.weight;

                w_sum += nu * nv * wi;
                w_u += dnu * nv * wi;
                for d in 0..3 {
                    a[d] += nu * nv * wi * cp.pos[d];
                    a_u[d] += dnu * nv * wi * cp.pos[d];
                }
            }
        }

        if w_sum.abs() < 1e-14 {
            return [0.0; 3];
        }
        let w2 = w_sum * w_sum;
        [
            (a_u[0] * w_sum - a[0] * w_u) / w2,
            (a_u[1] * w_sum - a[1] * w_u) / w2,
            (a_u[2] * w_sum - a[2] * w_u) / w2,
        ]
    }

    /// Partial derivative with respect to v at (u, v).
    pub fn deriv_v(&self, u: f64, v: f64) -> [f64; 3] {
        let mut a = [0.0; 3];
        let mut a_v = [0.0; 3];
        let mut w_sum = 0.0;
        let mut w_v = 0.0;

        for j in 0..self.num_v {
            let nv = self.basis_v.evaluate(j, v);
            let dnv = self.basis_v.derivative(j, v);
            for i in 0..self.num_u {
                let nu = self.basis_u.evaluate(i, u);
                let cp = self.control_point(i, j);
                let wi = cp.weight;

                w_sum += nu * nv * wi;
                w_v += nu * dnv * wi;
                for d in 0..3 {
                    a[d] += nu * nv * wi * cp.pos[d];
                    a_v[d] += nu * dnv * wi * cp.pos[d];
                }
            }
        }

        if w_sum.abs() < 1e-14 {
            return [0.0; 3];
        }
        let w2 = w_sum * w_sum;
        [
            (a_v[0] * w_sum - a[0] * w_v) / w2,
            (a_v[1] * w_sum - a[1] * w_v) / w2,
            (a_v[2] * w_sum - a[2] * w_v) / w2,
        ]
    }

    /// Compute the surface normal at (u, v).
    pub fn normal(&self, u: f64, v: f64) -> [f64; 3] {
        let du = self.deriv_u(u, v);
        let dv = self.deriv_v(u, v);
        let n = cross3d(du, dv);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-14 {
            return [0.0, 0.0, 1.0];
        }
        [n[0] / len, n[1] / len, n[2] / len]
    }

    /// Create a flat rectangular NURBS surface in the xy-plane.
    pub fn flat_rectangle(width: f64, height: f64, nu: usize, nv: usize, degree: usize) -> Self {
        let basis_u = BSplineBasis::uniform_open(nu, degree);
        let basis_v = BSplineBasis::uniform_open(nv, degree);
        let mut cps = Vec::with_capacity(nu * nv);
        for j in 0..nv {
            let y = height * j as f64 / (nv - 1) as f64;
            for i in 0..nu {
                let x = width * i as f64 / (nu - 1) as f64;
                cps.push(WeightedPoint::unweighted([x, y, 0.0]));
            }
        }
        Self::new(basis_u, basis_v, cps, nu, nv)
    }
}

// ---------------------------------------------------------------------------
// IGA Element
// ---------------------------------------------------------------------------

/// Material properties for IGA analysis.
#[derive(Debug, Clone, Copy)]
pub struct IgaMaterial {
    /// Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Thickness (for 2D plane stress/strain).
    pub thickness: f64,
}

impl IgaMaterial {
    /// Create a new material.
    pub fn new(young: f64, poisson: f64, thickness: f64) -> Self {
        Self {
            young,
            poisson,
            thickness,
        }
    }

    /// Plane stress constitutive matrix (3x3 Voigt).
    pub fn plane_stress_matrix(&self) -> [[f64; 3]; 3] {
        let e = self.young;
        let nu = self.poisson;
        let factor = e / (1.0 - nu * nu);
        [
            [factor, factor * nu, 0.0],
            [factor * nu, factor, 0.0],
            [0.0, 0.0, factor * (1.0 - nu) / 2.0],
        ]
    }

    /// Plane strain constitutive matrix (3x3 Voigt).
    pub fn plane_strain_matrix(&self) -> [[f64; 3]; 3] {
        let e = self.young;
        let nu = self.poisson;
        let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        [
            [factor * (1.0 - nu), factor * nu, 0.0],
            [factor * nu, factor * (1.0 - nu), 0.0],
            [0.0, 0.0, factor * (1.0 - 2.0 * nu) / 2.0],
        ]
    }
}

/// An isogeometric element defined over a knot span.
///
/// Uses NURBS shape functions and Gauss quadrature for integration.
#[derive(Debug, Clone)]
pub struct IgaElement {
    /// Indices of active basis functions for this element.
    pub basis_indices: Vec<usize>,
    /// Parametric domain \[u_min, u_max, v_min, v_max\].
    pub param_domain: [f64; 4],
    /// Number of Gauss points per direction.
    pub num_gauss: usize,
}

impl IgaElement {
    /// Create a new IGA element.
    pub fn new(basis_indices: Vec<usize>, param_domain: [f64; 4], num_gauss: usize) -> Self {
        Self {
            basis_indices,
            param_domain,
            num_gauss,
        }
    }

    /// Compute NURBS shape functions at parametric point (u, v) for a surface.
    pub fn shape_functions(&self, surface: &NurbsSurface, u: f64, v: f64) -> Vec<f64> {
        let nu_vals: Vec<f64> = (0..surface.num_u)
            .map(|i| surface.basis_u.evaluate(i, u))
            .collect();
        let nv_vals: Vec<f64> = (0..surface.num_v)
            .map(|j| surface.basis_v.evaluate(j, v))
            .collect();

        let mut weighted = Vec::new();
        let mut w_sum = 0.0;
        for j in 0..surface.num_v {
            for i in 0..surface.num_u {
                let idx = j * surface.num_u + i;
                let val = nu_vals[i] * nv_vals[j] * surface.control_points[idx].weight;
                weighted.push(val);
                w_sum += val;
            }
        }

        if w_sum.abs() < 1e-14 {
            return vec![0.0; weighted.len()];
        }
        weighted.iter().map(|v| v / w_sum).collect()
    }

    /// Compute the Jacobian matrix at parametric point (u, v).
    ///
    /// Returns 2x2 Jacobian for a 2D surface embedded in 2D.
    pub fn jacobian_2d(&self, surface: &NurbsSurface, u: f64, v: f64) -> [[f64; 2]; 2] {
        let du = surface.deriv_u(u, v);
        let dv = surface.deriv_v(u, v);
        [
            [du[0], dv[0]], // dx/du, dx/dv
            [du[1], dv[1]], // dy/du, dy/dv
        ]
    }

    /// Determinant of the 2x2 Jacobian.
    pub fn jacobian_det(&self, surface: &NurbsSurface, u: f64, v: f64) -> f64 {
        let j = self.jacobian_2d(surface, u, v);
        j[0][0] * j[1][1] - j[0][1] * j[1][0]
    }

    /// Compute the B-matrix (strain-displacement) at (u, v).
    ///
    /// Returns B as a flat vector of \[3 x (2*n_basis)\] entries for 2D plane stress.
    pub fn b_matrix(&self, surface: &NurbsSurface, u: f64, v: f64) -> (Vec<f64>, f64) {
        let jac = self.jacobian_2d(surface, u, v);
        let det_j = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        if det_j.abs() < 1e-14 {
            let n_total = surface.num_u * surface.num_v;
            return (vec![0.0; 3 * 2 * n_total], det_j);
        }

        let inv_j = [
            [jac[1][1] / det_j, -jac[0][1] / det_j],
            [-jac[1][0] / det_j, jac[0][0] / det_j],
        ];

        // Compute dN/du and dN/dv for all basis functions
        let nu_vals: Vec<f64> = (0..surface.num_u)
            .map(|i| surface.basis_u.evaluate(i, u))
            .collect();
        let nv_vals: Vec<f64> = (0..surface.num_v)
            .map(|j| surface.basis_v.evaluate(j, v))
            .collect();
        let dnu_vals: Vec<f64> = (0..surface.num_u)
            .map(|i| surface.basis_u.derivative(i, u))
            .collect();
        let dnv_vals: Vec<f64> = (0..surface.num_v)
            .map(|j| surface.basis_v.derivative(j, v))
            .collect();

        // Compute weighted derivatives
        let n_total = surface.num_u * surface.num_v;
        let mut nw = vec![0.0; n_total];
        let mut dnw_du = vec![0.0; n_total];
        let mut dnw_dv = vec![0.0; n_total];
        let mut w_sum = 0.0;
        let mut dw_du = 0.0;
        let mut dw_dv = 0.0;

        for j in 0..surface.num_v {
            for i in 0..surface.num_u {
                let idx = j * surface.num_u + i;
                let wi = surface.control_points[idx].weight;
                nw[idx] = nu_vals[i] * nv_vals[j] * wi;
                dnw_du[idx] = dnu_vals[i] * nv_vals[j] * wi;
                dnw_dv[idx] = nu_vals[i] * dnv_vals[j] * wi;
                w_sum += nw[idx];
                dw_du += dnw_du[idx];
                dw_dv += dnw_dv[idx];
            }
        }

        // Rational derivatives: dR/du = (dNw/du * W - Nw * dW/du) / W^2
        let w2 = w_sum * w_sum;
        let mut dr_du = vec![0.0; n_total];
        let mut dr_dv = vec![0.0; n_total];
        for idx in 0..n_total {
            dr_du[idx] = (dnw_du[idx] * w_sum - nw[idx] * dw_du) / w2;
            dr_dv[idx] = (dnw_dv[idx] * w_sum - nw[idx] * dw_dv) / w2;
        }

        // Physical derivatives via inverse Jacobian
        let mut dn_dx = vec![0.0; n_total];
        let mut dn_dy = vec![0.0; n_total];
        for idx in 0..n_total {
            dn_dx[idx] = inv_j[0][0] * dr_du[idx] + inv_j[0][1] * dr_dv[idx];
            dn_dy[idx] = inv_j[1][0] * dr_du[idx] + inv_j[1][1] * dr_dv[idx];
        }

        // B-matrix: 3 rows x (2*n_total) columns
        // B = [dN1/dx  0      dN2/dx  0     ...]
        //     [0       dN1/dy 0       dN2/dy...]
        //     [dN1/dy  dN1/dx dN2/dy  dN2/dx...]
        let mut b = vec![0.0; 3 * 2 * n_total];
        for idx in 0..n_total {
            // Row 0
            b[2 * idx] = dn_dx[idx];
            // Row 1
            b[(2 * n_total) + 2 * idx + 1] = dn_dy[idx];
            // Row 2
            b[2 * (2 * n_total) + 2 * idx] = dn_dy[idx];
            b[2 * (2 * n_total) + 2 * idx + 1] = dn_dx[idx];
        }

        (b, det_j)
    }

    /// Compute element stiffness matrix using Gauss quadrature.
    ///
    /// Returns a dense matrix of size (2*n_dof) x (2*n_dof) flattened.
    pub fn stiffness_matrix(&self, surface: &NurbsSurface, material: &IgaMaterial) -> Vec<f64> {
        let n_total = surface.num_u * surface.num_v;
        let ndof = 2 * n_total;
        let mut k = vec![0.0; ndof * ndof];

        let d_mat = material.plane_stress_matrix();

        let (gp_ref, gw_ref) = gauss_legendre(self.num_gauss);
        let (gp_u, gw_u) = gauss_map(&gp_ref, &gw_ref, self.param_domain[0], self.param_domain[1]);
        let (gp_v, gw_v) = gauss_map(&gp_ref, &gw_ref, self.param_domain[2], self.param_domain[3]);

        for (qu, wu) in gp_u.iter().zip(gw_u.iter()) {
            for (qv, wv) in gp_v.iter().zip(gw_v.iter()) {
                let (b_flat, det_j) = self.b_matrix(surface, *qu, *qv);
                if det_j.abs() < 1e-14 {
                    continue;
                }
                let weight = wu * wv * det_j.abs() * material.thickness;

                // K += B^T * D * B * weight
                // B is 3 x ndof, D is 3 x 3
                for row in 0..ndof {
                    for col in row..ndof {
                        let mut val = 0.0;
                        for s in 0..3 {
                            let b_s_row = b_flat[s * ndof + row];
                            for r in 0..3 {
                                let b_r_col = b_flat[r * ndof + col];
                                val += b_s_row * d_mat[s][r] * b_r_col;
                            }
                        }
                        val *= weight;
                        k[row * ndof + col] += val;
                        if row != col {
                            k[col * ndof + row] += val;
                        }
                    }
                }
            }
        }

        k
    }
}

// ---------------------------------------------------------------------------
// IGA Assembler
// ---------------------------------------------------------------------------

/// Assembled IGA system: stiffness matrix and load vector.
#[derive(Debug, Clone)]
pub struct IgaSystem {
    /// Global stiffness matrix (dense, n x n).
    pub stiffness: Vec<f64>,
    /// Global load vector (n).
    pub load: Vec<f64>,
    /// Total DOFs.
    pub ndof: usize,
}

/// Assembles the global IGA system from elements.
#[derive(Debug, Clone)]
pub struct IgaAssembler {
    /// Total number of control points (DOFs / 2 for 2D).
    pub num_control_points: usize,
    /// Material properties.
    pub material: IgaMaterial,
    /// Number of Gauss points per direction.
    pub num_gauss: usize,
}

impl IgaAssembler {
    /// Create a new assembler.
    pub fn new(num_control_points: usize, material: IgaMaterial, num_gauss: usize) -> Self {
        Self {
            num_control_points,
            material,
            num_gauss,
        }
    }

    /// Assemble the global stiffness matrix and load vector.
    ///
    /// For a single-patch surface, the element stiffness is the global stiffness.
    pub fn assemble(&self, surface: &NurbsSurface) -> IgaSystem {
        let ndof = 2 * self.num_control_points;
        let mut stiffness = vec![0.0; ndof * ndof];
        let load = vec![0.0; ndof];

        // Create elements from non-zero knot spans
        let elements = self.create_elements(surface);

        for elem in &elements {
            let k_elem = elem.stiffness_matrix(surface, &self.material);
            // For single-patch: element DOFs map directly to global DOFs
            let n_elem_dof = 2 * surface.num_u * surface.num_v;
            for i in 0..n_elem_dof {
                for j in 0..n_elem_dof {
                    stiffness[i * ndof + j] += k_elem[i * n_elem_dof + j];
                }
            }
        }

        IgaSystem {
            stiffness,
            load,
            ndof,
        }
    }

    /// Create IGA elements from non-zero knot spans.
    pub fn create_elements(&self, surface: &NurbsSurface) -> Vec<IgaElement> {
        let mut elements = Vec::new();
        let knots_u = &surface.basis_u.knots;
        let knots_v = &surface.basis_v.knots;

        for i in 0..knots_u.len() - 1 {
            if (knots_u[i + 1] - knots_u[i]).abs() < 1e-14 {
                continue;
            }
            for j in 0..knots_v.len() - 1 {
                if (knots_v[j + 1] - knots_v[j]).abs() < 1e-14 {
                    continue;
                }
                let n_total = surface.num_u * surface.num_v;
                let indices: Vec<usize> = (0..n_total).collect();
                let domain = [knots_u[i], knots_u[i + 1], knots_v[j], knots_v[j + 1]];
                elements.push(IgaElement::new(indices, domain, self.num_gauss));
            }
        }
        elements
    }

    /// Apply a uniform body force (e.g., gravity) to the load vector.
    pub fn apply_body_force(
        &self,
        surface: &NurbsSurface,
        system: &mut IgaSystem,
        force: [f64; 2],
    ) {
        let elements = self.create_elements(surface);
        let (gp_ref, gw_ref) = gauss_legendre(self.num_gauss);

        for elem in &elements {
            let (gp_u, gw_u) =
                gauss_map(&gp_ref, &gw_ref, elem.param_domain[0], elem.param_domain[1]);
            let (gp_v, gw_v) =
                gauss_map(&gp_ref, &gw_ref, elem.param_domain[2], elem.param_domain[3]);

            for (qu, wu) in gp_u.iter().zip(gw_u.iter()) {
                for (qv, wv) in gp_v.iter().zip(gw_v.iter()) {
                    let shapes = elem.shape_functions(surface, *qu, *qv);
                    let det_j = elem.jacobian_det(surface, *qu, *qv);
                    let weight = wu * wv * det_j.abs() * self.material.thickness;

                    for (idx, n_val) in shapes.iter().enumerate() {
                        system.load[2 * idx] += n_val * force[0] * weight;
                        system.load[2 * idx + 1] += n_val * force[1] * weight;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// IGA Boundary Conditions
// ---------------------------------------------------------------------------

/// Type of boundary condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IgaBcType {
    /// Fixed displacement (Dirichlet).
    Dirichlet,
    /// Applied traction (Neumann).
    Neumann,
}

/// A boundary condition for IGA.
#[derive(Debug, Clone)]
pub struct IgaBoundaryCondition {
    /// DOF indices affected.
    pub dof_indices: Vec<usize>,
    /// Prescribed values.
    pub values: Vec<f64>,
    /// Type of BC.
    pub bc_type: IgaBcType,
}

impl IgaBoundaryCondition {
    /// Create a new boundary condition.
    pub fn new(dof_indices: Vec<usize>, values: Vec<f64>, bc_type: IgaBcType) -> Self {
        assert_eq!(dof_indices.len(), values.len());
        Self {
            dof_indices,
            values,
            bc_type,
        }
    }

    /// Apply Dirichlet BC to the system using the penalty method.
    pub fn apply_dirichlet(&self, system: &mut IgaSystem, penalty: f64) {
        if self.bc_type != IgaBcType::Dirichlet {
            return;
        }
        for (dof, val) in self.dof_indices.iter().zip(self.values.iter()) {
            let d = *dof;
            system.stiffness[d * system.ndof + d] += penalty;
            system.load[d] += penalty * val;
        }
    }

    /// Apply Neumann BC (add to load vector).
    pub fn apply_neumann(&self, system: &mut IgaSystem) {
        if self.bc_type != IgaBcType::Neumann {
            return;
        }
        for (dof, val) in self.dof_indices.iter().zip(self.values.iter()) {
            system.load[*dof] += val;
        }
    }

    /// Apply this BC to the system (dispatches based on type).
    pub fn apply(&self, system: &mut IgaSystem, penalty: f64) {
        match self.bc_type {
            IgaBcType::Dirichlet => self.apply_dirichlet(system, penalty),
            IgaBcType::Neumann => self.apply_neumann(system),
        }
    }

    /// Create a fixed boundary condition (zero displacement) for given DOFs.
    pub fn fixed(dof_indices: Vec<usize>) -> Self {
        let values = vec![0.0; dof_indices.len()];
        Self::new(dof_indices, values, IgaBcType::Dirichlet)
    }

    /// Create a traction BC for given DOFs with specified forces.
    pub fn traction(dof_indices: Vec<usize>, forces: Vec<f64>) -> Self {
        Self::new(dof_indices, forces, IgaBcType::Neumann)
    }
}

// ---------------------------------------------------------------------------
// Trimmed NURBS
// ---------------------------------------------------------------------------

/// A trimming curve defined in the parametric space of a NURBS surface.
#[derive(Debug, Clone)]
pub struct TrimmingCurve {
    /// Control points (u, v) in parametric space.
    pub control_points: Vec<[f64; 2]>,
    /// Whether this is an outer boundary (true) or a hole (false).
    pub is_outer: bool,
}

impl TrimmingCurve {
    /// Create a new trimming curve.
    pub fn new(control_points: Vec<[f64; 2]>, is_outer: bool) -> Self {
        Self {
            control_points,
            is_outer,
        }
    }

    /// Check if a parametric point (u, v) is inside this trimming curve.
    ///
    /// Uses the winding number algorithm.
    pub fn contains(&self, u: f64, v: f64) -> bool {
        let n = self.control_points.len();
        if n < 3 {
            return false;
        }
        let mut winding = 0i32;
        for i in 0..n {
            let p1 = self.control_points[i];
            let p2 = self.control_points[(i + 1) % n];
            if p1[1] <= v {
                if p2[1] > v {
                    let vt = (v - p1[1]) / (p2[1] - p1[1]);
                    let x_intersect = p1[0] + vt * (p2[0] - p1[0]);
                    if u < x_intersect {
                        winding += 1;
                    }
                }
            } else if p2[1] <= v {
                let vt = (v - p1[1]) / (p2[1] - p1[1]);
                let x_intersect = p1[0] + vt * (p2[0] - p1[0]);
                if u < x_intersect {
                    winding -= 1;
                }
            }
        }
        winding != 0
    }

    /// Compute approximate perimeter of the trimming curve.
    pub fn perimeter(&self) -> f64 {
        let n = self.control_points.len();
        let mut length = 0.0;
        for i in 0..n {
            let p1 = self.control_points[i];
            let p2 = self.control_points[(i + 1) % n];
            let dx = p2[0] - p1[0];
            let dy = p2[1] - p1[1];
            length += (dx * dx + dy * dy).sqrt();
        }
        length
    }

    /// Compute approximate area inside the trimming curve via shoelace formula.
    pub fn area(&self) -> f64 {
        let n = self.control_points.len();
        let mut area = 0.0;
        for i in 0..n {
            let p1 = self.control_points[i];
            let p2 = self.control_points[(i + 1) % n];
            area += p1[0] * p2[1] - p2[0] * p1[1];
        }
        (area / 2.0).abs()
    }
}

/// A trimmed NURBS surface.
///
/// Consists of a base NURBS surface and a set of trimming curves
/// that define the active region.
#[derive(Debug, Clone)]
pub struct TrimmedNurbs {
    /// The base NURBS surface.
    pub surface: NurbsSurface,
    /// Trimming curves.
    pub trimming_curves: Vec<TrimmingCurve>,
}

impl TrimmedNurbs {
    /// Create a new trimmed NURBS surface.
    pub fn new(surface: NurbsSurface, trimming_curves: Vec<TrimmingCurve>) -> Self {
        Self {
            surface,
            trimming_curves,
        }
    }

    /// Check if a parametric point (u, v) is in the active (untrimmed) region.
    pub fn is_active(&self, u: f64, v: f64) -> bool {
        let mut inside_outer = false;
        let mut in_hole = false;

        for curve in &self.trimming_curves {
            if curve.is_outer {
                if curve.contains(u, v) {
                    inside_outer = true;
                }
            } else if curve.contains(u, v) {
                in_hole = true;
            }
        }

        // If no outer boundary, assume entire domain is active
        let has_outer = self.trimming_curves.iter().any(|c| c.is_outer);
        if !has_outer {
            inside_outer = true;
        }

        inside_outer && !in_hole
    }

    /// Evaluate the surface at (u, v) if the point is in the active region.
    pub fn evaluate(&self, u: f64, v: f64) -> Option<[f64; 3]> {
        if self.is_active(u, v) {
            Some(self.surface.evaluate(u, v))
        } else {
            None
        }
    }

    /// Compute the approximate trimmed area using Monte Carlo sampling.
    pub fn trimmed_area_mc(&self, n_samples: usize) -> f64 {
        let u_min = self.surface.basis_u.knots[self.surface.basis_u.degree];
        let u_max = self.surface.basis_u.knots
            [self.surface.basis_u.knots.len() - self.surface.basis_u.degree - 1];
        let v_min = self.surface.basis_v.knots[self.surface.basis_v.degree];
        let v_max = self.surface.basis_v.knots
            [self.surface.basis_v.knots.len() - self.surface.basis_v.degree - 1];

        let total_area = (u_max - u_min) * (v_max - v_min);
        let mut rng = rand::rng();
        let mut count = 0;
        for _ in 0..n_samples {
            let u: f64 = rand::RngExt::random_range(&mut rng, u_min..u_max);
            let v: f64 = rand::RngExt::random_range(&mut rng, v_min..v_max);
            if self.is_active(u, v) {
                count += 1;
            }
        }
        total_area * count as f64 / n_samples as f64
    }

    /// Compute integration points for the trimmed surface.
    ///
    /// Returns (u, v, weight) triples for active Gauss points.
    pub fn integration_points(&self, num_gauss: usize) -> Vec<(f64, f64, f64)> {
        let (gp_ref, gw_ref) = gauss_legendre(num_gauss);

        let knots_u = &self.surface.basis_u.knots;
        let knots_v = &self.surface.basis_v.knots;
        let mut points = Vec::new();

        for i in 0..knots_u.len() - 1 {
            if (knots_u[i + 1] - knots_u[i]).abs() < 1e-14 {
                continue;
            }
            let (gp_u, gw_u) = gauss_map(&gp_ref, &gw_ref, knots_u[i], knots_u[i + 1]);
            for j in 0..knots_v.len() - 1 {
                if (knots_v[j + 1] - knots_v[j]).abs() < 1e-14 {
                    continue;
                }
                let (gp_v, gw_v) = gauss_map(&gp_ref, &gw_ref, knots_v[j], knots_v[j + 1]);
                for (qu, wu) in gp_u.iter().zip(gw_u.iter()) {
                    for (qv, wv) in gp_v.iter().zip(gw_v.iter()) {
                        if self.is_active(*qu, *qv) {
                            points.push((*qu, *qv, wu * wv));
                        }
                    }
                }
            }
        }
        points
    }
}

// ---------------------------------------------------------------------------
// Knot Refinement
// ---------------------------------------------------------------------------

/// Knot insertion (h-refinement) for B-spline basis.
///
/// Inserts a new knot into the knot vector and updates control points.
pub fn knot_insertion_curve(curve: &NurbsCurve, t_new: f64) -> NurbsCurve {
    let p = curve.basis.degree;
    let knots = &curve.basis.knots;
    let n = curve.control_points.len();

    // Find span
    let mut k = 0;
    for i in 0..knots.len() - 1 {
        if knots[i] <= t_new && t_new < knots[i + 1] {
            k = i;
            break;
        }
    }

    // New knot vector
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(t_new);
    new_knots.extend_from_slice(&knots[k + 1..]);

    // New control points (in homogeneous coords)
    let mut new_cps = Vec::with_capacity(n + 1);
    for i in 0..=n {
        if i <= k.saturating_sub(p) {
            new_cps.push(curve.control_points[i]);
        } else if i > k {
            if i - 1 < curve.control_points.len() {
                new_cps.push(curve.control_points[i - 1]);
            }
        } else {
            let idx = i;
            if idx < n && idx > 0 {
                let alpha = if (knots[idx + p] - knots[idx]).abs() > 1e-14 {
                    (t_new - knots[idx]) / (knots[idx + p] - knots[idx])
                } else {
                    0.5
                };
                let cp1 = &curve.control_points[idx - 1];
                let cp2 = &curve.control_points[idx];
                let new_pos = [
                    (1.0 - alpha) * cp1.pos[0] * cp1.weight + alpha * cp2.pos[0] * cp2.weight,
                    (1.0 - alpha) * cp1.pos[1] * cp1.weight + alpha * cp2.pos[1] * cp2.weight,
                    (1.0 - alpha) * cp1.pos[2] * cp1.weight + alpha * cp2.pos[2] * cp2.weight,
                ];
                let new_w = (1.0 - alpha) * cp1.weight + alpha * cp2.weight;
                new_cps.push(WeightedPoint::new(
                    [new_pos[0] / new_w, new_pos[1] / new_w, new_pos[2] / new_w],
                    new_w,
                ));
            } else if idx < n {
                new_cps.push(curve.control_points[idx]);
            }
        }
    }

    let new_basis = BSplineBasis::new(new_knots, p);
    NurbsCurve::new(new_basis, new_cps)
}

/// Degree elevation (p-refinement): raises degree by 1.
///
/// Returns a new basis with degree p+1 and appropriately modified knot vector.
pub fn degree_elevation(basis: &BSplineBasis) -> BSplineBasis {
    let p = basis.degree;
    let knots = &basis.knots;

    // Insert each unique interior knot once more and increase degree
    let mut unique_knots = Vec::new();
    let mut prev = f64::NEG_INFINITY;
    for &k in knots {
        if (k - prev).abs() > 1e-14 {
            unique_knots.push(k);
            prev = k;
        }
    }

    // New knot vector: each unique knot gets one additional copy
    let mut new_knots = Vec::new();
    for &k in knots {
        new_knots.push(k);
    }
    // Add extra copies at the endpoints for the degree increase
    new_knots.insert(0, knots[0]);
    new_knots.push(*knots.last().expect("knots is non-empty"));
    // Add one more copy of each unique interior knot
    for &uk in &unique_knots[1..unique_knots.len() - 1] {
        new_knots.push(uk);
    }
    new_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    BSplineBasis::new(new_knots, p + 1)
}

// ---------------------------------------------------------------------------
// Simple IGA solver (CG for Ku = f)
// ---------------------------------------------------------------------------

/// Solve a dense linear system Ax = b using conjugate gradient.
///
/// Returns the solution vector.
pub fn solve_cg(a: &[f64], b: &[f64], n: usize, max_iter: usize, tol: f64) -> Vec<f64> {
    let mut x = vec![0.0; n];
    let mut r: Vec<f64> = b.to_vec();
    let mut p = r.clone();
    let mut rs_old: f64 = r.iter().map(|v| v * v).sum();

    for _iter in 0..max_iter {
        if rs_old.sqrt() < tol {
            break;
        }

        // Compute Ap
        let mut ap = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                ap[i] += a[i * n + j] * p[j];
            }
        }

        let p_ap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
        if p_ap.abs() < 1e-30 {
            break;
        }
        let alpha = rs_old / p_ap;

        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        let rs_new: f64 = r.iter().map(|v| v * v).sum();
        if rs_new.sqrt() < tol {
            break;
        }

        let beta = rs_new / rs_old;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }
    x
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bspline_partition_of_unity() {
        let basis = BSplineBasis::uniform_open(5, 2);
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let vals = basis.evaluate_all_dense(t);
            let sum: f64 = vals.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "Partition of unity failed at t={t}: sum={sum}"
            );
        }
    }

    #[test]
    fn test_bspline_nonnegativity() {
        let basis = BSplineBasis::uniform_open(6, 3);
        for i in 0..20 {
            let t = i as f64 / 20.0;
            let vals = basis.evaluate_all_dense(t);
            for (j, v) in vals.iter().enumerate() {
                assert!(*v >= -1e-14, "Basis function {j} negative at t={t}: {v}");
            }
        }
    }

    #[test]
    fn test_bspline_degree_zero() {
        let knots = vec![0.0, 1.0, 2.0, 3.0];
        let basis = BSplineBasis::new(knots, 0);
        assert_eq!(basis.num_basis(), 3);
        assert!((basis.evaluate(0, 0.5) - 1.0).abs() < 1e-14);
        assert!(basis.evaluate(1, 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_bspline_derivative_sum_zero() {
        // Sum of derivatives of a partition of unity should be zero
        let basis = BSplineBasis::uniform_open(5, 2);
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let derivs = basis.derivative_all(t);
            let sum: f64 = derivs.iter().sum();
            assert!(
                sum.abs() < 1e-8,
                "Derivative sum should be 0 at t={t}: sum={sum}"
            );
        }
    }

    #[test]
    fn test_bspline_find_span() {
        let basis = BSplineBasis::uniform_open(5, 2);
        let span = basis.find_span(0.5);
        assert!(span >= basis.degree);
        assert!(span < basis.num_basis());
    }

    #[test]
    fn test_bspline_greville() {
        let basis = BSplineBasis::uniform_open(4, 2);
        let grev = basis.greville_abscissae();
        assert_eq!(grev.len(), 4);
        // Greville abscissae should be in [0, 1]
        for g in &grev {
            assert!(*g >= 0.0 && *g <= 1.0);
        }
    }

    #[test]
    fn test_nurbs_curve_line() {
        // A straight line from (0,0,0) to (1,0,0) with uniform weights
        let basis = BSplineBasis::uniform_open(2, 1);
        let cps = vec![
            WeightedPoint::unweighted([0.0, 0.0, 0.0]),
            WeightedPoint::unweighted([1.0, 0.0, 0.0]),
        ];
        let curve = NurbsCurve::new(basis, cps);
        let mid = curve.evaluate(0.5);
        assert!((mid[0] - 0.5).abs() < 1e-10);
        assert!(mid[1].abs() < 1e-10);
    }

    #[test]
    fn test_nurbs_circle_radius() {
        let r = 2.0;
        let circle = NurbsCurve::circle(r);
        // Sample points and check radius
        for i in 1..8 {
            let t = i as f64 / 8.0;
            let pt = circle.evaluate(t);
            let dist = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
            assert!(
                (dist - r).abs() < 1e-6,
                "Circle point at t={t} has radius {dist}, expected {r}"
            );
        }
    }

    #[test]
    fn test_nurbs_circle_closed() {
        let circle = NurbsCurve::circle(1.0);
        let start = circle.evaluate(0.0);
        let end = circle.evaluate(1.0);
        assert!((start[0] - end[0]).abs() < 1e-10);
        assert!((start[1] - end[1]).abs() < 1e-10);
    }

    #[test]
    fn test_nurbs_curve_derivative() {
        let basis = BSplineBasis::uniform_open(3, 2);
        let cps = vec![
            WeightedPoint::unweighted([0.0, 0.0, 0.0]),
            WeightedPoint::unweighted([0.5, 1.0, 0.0]),
            WeightedPoint::unweighted([1.0, 0.0, 0.0]),
        ];
        let curve = NurbsCurve::new(basis, cps);
        let d = curve.derivative(0.5);
        // Should be tangent at midpoint
        assert!(
            d[0].abs() > 1e-6 || d[1].abs() > 1e-6,
            "Derivative should be non-zero"
        );
    }

    #[test]
    fn test_nurbs_rational_basis_partition() {
        let circle = NurbsCurve::circle(1.0);
        for i in 1..8 {
            let t = i as f64 / 8.0;
            let r_basis = circle.rational_basis(t);
            let sum: f64 = r_basis.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-10,
                "Rational basis partition of unity failed at t={t}"
            );
        }
    }

    #[test]
    fn test_surface_flat_eval() {
        let surf = NurbsSurface::flat_rectangle(2.0, 1.0, 3, 3, 1);
        let pt = surf.evaluate(0.5, 0.5);
        assert!((pt[0] - 1.0).abs() < 1e-10);
        assert!((pt[1] - 0.5).abs() < 1e-10);
        assert!(pt[2].abs() < 1e-10);
    }

    #[test]
    fn test_surface_corners() {
        let surf = NurbsSurface::flat_rectangle(3.0, 2.0, 4, 4, 2);
        let p00 = surf.evaluate(0.0, 0.0);
        let p10 = surf.evaluate(1.0, 0.0);
        let p01 = surf.evaluate(0.0, 1.0);
        let p11 = surf.evaluate(1.0, 1.0);
        assert!((p00[0]).abs() < 1e-10);
        assert!((p10[0] - 3.0).abs() < 1e-10);
        assert!((p01[1] - 2.0).abs() < 1e-10);
        assert!((p11[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_surface_normal_flat() {
        let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 3, 3, 1);
        let n = surf.normal(0.5, 0.5);
        // Normal should be (0, 0, +/-1) for a flat xy surface
        assert!(n[2].abs() > 0.99, "Normal z-component should be ~1");
    }

    #[test]
    fn test_stiffness_symmetry() {
        let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 3, 3, 1);
        let mat = IgaMaterial::new(1e6, 0.3, 0.01);
        let elem = IgaElement::new(vec![], [0.0, 0.5, 0.0, 0.5], 2);
        let k = elem.stiffness_matrix(&surf, &mat);
        let ndof = 2 * 9;
        for i in 0..ndof {
            for j in i + 1..ndof {
                let diff = (k[i * ndof + j] - k[j * ndof + i]).abs();
                let max_val = k[i * ndof + j].abs().max(k[j * ndof + i].abs()).max(1e-30);
                assert!(
                    diff / max_val < 1e-8 || diff < 1e-20,
                    "K not symmetric at ({i},{j}): {} vs {}",
                    k[i * ndof + j],
                    k[j * ndof + i]
                );
            }
        }
    }

    #[test]
    fn test_stiffness_positive_diagonal() {
        let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 3, 3, 1);
        let mat = IgaMaterial::new(1e6, 0.3, 0.01);
        let elem = IgaElement::new(vec![], [0.0, 0.5, 0.0, 0.5], 2);
        let k = elem.stiffness_matrix(&surf, &mat);
        let ndof = 2 * 9;
        let mut has_positive = false;
        for i in 0..ndof {
            assert!(k[i * ndof + i] >= -1e-14, "Negative diagonal at {i}");
            if k[i * ndof + i] > 1e-10 {
                has_positive = true;
            }
        }
        assert!(has_positive, "Should have some positive diagonal entries");
    }

    #[test]
    fn test_material_plane_stress() {
        let mat = IgaMaterial::new(200e9, 0.3, 1.0);
        let d = mat.plane_stress_matrix();
        // D should be symmetric
        assert!((d[0][1] - d[1][0]).abs() < 1e-3);
        // D11 = D22
        assert!((d[0][0] - d[1][1]).abs() < 1e-3);
    }

    #[test]
    fn test_material_plane_strain() {
        let mat = IgaMaterial::new(200e9, 0.3, 1.0);
        let d = mat.plane_strain_matrix();
        assert!((d[0][1] - d[1][0]).abs() < 1e-3);
    }

    #[test]
    fn test_bc_dirichlet_penalty() {
        let ndof = 4;
        let mut system = IgaSystem {
            stiffness: vec![0.0; ndof * ndof],
            load: vec![0.0; ndof],
            ndof,
        };
        let bc = IgaBoundaryCondition::fixed(vec![0, 1]);
        bc.apply(&mut system, 1e10);
        assert!(system.stiffness[0] > 1e9);
        assert!(system.stiffness[ndof + 1] > 1e9);
    }

    #[test]
    fn test_bc_neumann() {
        let ndof = 4;
        let mut system = IgaSystem {
            stiffness: vec![0.0; ndof * ndof],
            load: vec![0.0; ndof],
            ndof,
        };
        let bc = IgaBoundaryCondition::traction(vec![2, 3], vec![100.0, -50.0]);
        bc.apply(&mut system, 1e10);
        assert!((system.load[2] - 100.0).abs() < 1e-10);
        assert!((system.load[3] - (-50.0)).abs() < 1e-10);
    }

    #[test]
    fn test_trimming_curve_circle() {
        // Create a circular trimming curve
        let n = 32;
        let cps: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                [0.5 + 0.3 * theta.cos(), 0.5 + 0.3 * theta.sin()]
            })
            .collect();
        let curve = TrimmingCurve::new(cps, true);
        assert!(curve.contains(0.5, 0.5)); // Center inside
        assert!(!curve.contains(0.0, 0.0)); // Outside
    }

    #[test]
    fn test_trimming_curve_area() {
        // Square trimming curve
        let cps = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let curve = TrimmingCurve::new(cps, true);
        let area = curve.area();
        assert!((area - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_trimmed_nurbs_active_region() {
        let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 3, 3, 1);
        let outer = TrimmingCurve::new(vec![[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]], true);
        let trimmed = TrimmedNurbs::new(surf, vec![outer]);
        assert!(trimmed.is_active(0.5, 0.5));
        assert!(!trimmed.is_active(0.05, 0.05));
    }

    #[test]
    fn test_gauss_legendre_exactness() {
        // 2-point rule should integrate x^3 exactly? No, up to degree 2*n-1=3
        let (pts, wts) = gauss_legendre(2);
        // Integrate x^2 on [-1,1]: exact = 2/3
        let integral: f64 = pts.iter().zip(wts.iter()).map(|(x, w)| x * x * w).sum();
        assert!((integral - 2.0 / 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_gauss_map() {
        let (pts, wts) = gauss_legendre(3);
        let (mapped_pts, mapped_wts) = gauss_map(&pts, &wts, 0.0, 2.0);
        // Integrate 1 on [0, 2]: should be 2
        let integral: f64 = mapped_wts.iter().sum();
        assert!((integral - 2.0).abs() < 1e-14);
        // All points should be in [0, 2]
        for p in &mapped_pts {
            assert!(*p >= 0.0 && *p <= 2.0);
        }
    }

    #[test]
    fn test_assembler_creates_elements() {
        let surf = NurbsSurface::flat_rectangle(1.0, 1.0, 4, 4, 2);
        let mat = IgaMaterial::new(1e6, 0.3, 0.01);
        let assembler = IgaAssembler::new(16, mat, 2);
        let elements = assembler.create_elements(&surf);
        assert!(!elements.is_empty());
    }

    #[test]
    fn test_solve_cg_identity() {
        // Solve I * x = b
        let n = 3;
        let a = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let b = vec![1.0, 2.0, 3.0];
        let x = solve_cg(&a, &b, n, 100, 1e-12);
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((x[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_knot_insertion_preserves_curve() {
        let basis = BSplineBasis::uniform_open(4, 2);
        let cps = vec![
            WeightedPoint::unweighted([0.0, 0.0, 0.0]),
            WeightedPoint::unweighted([1.0, 1.0, 0.0]),
            WeightedPoint::unweighted([2.0, 0.0, 0.0]),
            WeightedPoint::unweighted([3.0, 1.0, 0.0]),
        ];
        let curve = NurbsCurve::new(basis, cps);
        let refined = knot_insertion_curve(&curve, 0.5);

        // Evaluate at several points - should be the same
        for i in 0..5 {
            let t = i as f64 / 5.0 + 0.1;
            let p1 = curve.evaluate(t);
            let p2 = refined.evaluate(t);
            assert!(
                (p1[0] - p2[0]).abs() < 1e-6 && (p1[1] - p2[1]).abs() < 1e-6,
                "Knot insertion changed curve at t={t}"
            );
        }
    }

    #[test]
    fn test_degree_elevation() {
        let basis = BSplineBasis::uniform_open(4, 2);
        let elevated = degree_elevation(&basis);
        assert_eq!(elevated.degree, 3);
        assert!(elevated.knots.len() > basis.knots.len());
    }

    #[test]
    fn test_second_derivative() {
        let basis = BSplineBasis::uniform_open(5, 3);
        // Second derivative of cubic spline should exist
        let d2 = basis.second_derivative(2, 0.5);
        // Just check it returns a finite value
        assert!(d2.is_finite());
    }

    #[test]
    fn test_arc_length_line() {
        // Straight line from (0,0,0) to (1,0,0)
        let basis = BSplineBasis::uniform_open(2, 1);
        let cps = vec![
            WeightedPoint::unweighted([0.0, 0.0, 0.0]),
            WeightedPoint::unweighted([1.0, 0.0, 0.0]),
        ];
        let curve = NurbsCurve::new(basis, cps);
        let length = curve.arc_length(5);
        assert!(
            (length - 1.0).abs() < 1e-6,
            "Line length should be 1.0, got {length}"
        );
    }

    #[test]
    fn test_surface_deriv_flat() {
        let surf = NurbsSurface::flat_rectangle(2.0, 3.0, 3, 3, 1);
        let du = surf.deriv_u(0.5, 0.5);
        let dv = surf.deriv_v(0.5, 0.5);
        // du should be roughly (2, 0, 0) and dv roughly (0, 3, 0)
        assert!((du[0] - 2.0).abs() < 1e-6);
        assert!(du[1].abs() < 1e-6);
        assert!((dv[1] - 3.0).abs() < 1e-6);
        assert!(dv[0].abs() < 1e-6);
    }

    #[test]
    fn test_helper_cross3d() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3d(a, b);
        assert!((c[2] - 1.0).abs() < 1e-14);
    }
}
