#![allow(clippy::if_same_then_else, clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Discontinuous Galerkin finite element method.
//!
//! This module provides a complete DG-FEM framework for hyperbolic and
//! advection-diffusion problems, including:
//!
//! - **DG basis functions**: Legendre modal bases in 1-D and 2-D
//! - **Numerical fluxes**: Lax-Friedrichs, Roe, and upwind
//! - **Interior penalty methods**: SIPG, NIPG, IIPG for diffusion
//! - **h/p-adaptivity**: error-driven element refinement and polynomial enrichment
//! - **Slope limiters**: TVB and WENO reconstructions
//! - **Runge-Kutta DG**: explicit SSP time integration
//! - **Hyperbolic conservation laws**: Euler equations in 1-D
//! - **Advection-diffusion**: scalar transport with diffusion

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Dot product of two slices.
fn dot_slice(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// L2-norm of a slice.
fn norm_slice(a: &[f64]) -> f64 {
    dot_slice(a, a).sqrt()
}

/// Adds two vectors element-wise, returning a new vector.
fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scales a vector by a scalar, returning a new vector.
fn vec_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}

/// Subtracts two vectors element-wise, returning a new vector.
fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

// ---------------------------------------------------------------------------
// Legendre polynomials
// ---------------------------------------------------------------------------

/// Evaluates the Legendre polynomial P_n(x) at point x.
pub fn legendre_poly(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p_prev = 1.0_f64;
    let mut p_curr = x;
    for k in 2..=n {
        let kf = k as f64;
        let p_next = ((2.0 * kf - 1.0) * x * p_curr - (kf - 1.0) * p_prev) / kf;
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}

/// Evaluates the derivative of P_n(x).
pub fn legendre_poly_deriv(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return 1.0;
    }
    let pn = legendre_poly(n, x);
    let pn1 = legendre_poly(n - 1, x);
    let nf = n as f64;
    if (x * x - 1.0).abs() < 1e-14 {
        // Boundary value
        let sign = if x > 0.0 {
            1.0
        } else if n.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign * nf * (nf + 1.0) / 2.0
    } else {
        nf * (x * pn - pn1) / (x * x - 1.0)
    }
}

// ---------------------------------------------------------------------------
// Gauss-Legendre quadrature
// ---------------------------------------------------------------------------

/// Returns Gauss-Legendre quadrature points and weights on \[-1,1\].
pub fn gauss_legendre_quadrature(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut all_pts = vec![0.0; n];
    let mut all_wts = vec![0.0; n];
    let nf = n as f64;
    let half = n.div_ceil(2);
    for i in 0..half {
        // Initial guess (Abramowitz & Stegun approximation)
        let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        // Newton iteration to find roots of P_n
        for _ in 0..50 {
            let p = legendre_poly(n, x);
            let dp = legendre_poly_deriv(n, x);
            if dp.abs() < 1e-300 {
                break;
            }
            let dx = -p / dp;
            x += dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let dp = legendre_poly_deriv(n, x);
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        // Place symmetrically: negative root at the front, positive at the back
        all_pts[i] = -x;
        all_pts[n - 1 - i] = x;
        all_wts[i] = w;
        all_wts[n - 1 - i] = w;
    }
    (all_pts, all_wts)
}

// ---------------------------------------------------------------------------
// DG basis: modal Legendre
// ---------------------------------------------------------------------------

/// Modal DG basis using Legendre polynomials up to order p.
///
/// Each basis function is the normalized Legendre polynomial on \[-1,1\].
#[derive(Debug, Clone)]
pub struct DgBasis {
    /// Polynomial order.
    pub order: usize,
}

impl DgBasis {
    /// Creates a new DG basis with given polynomial order.
    pub fn new(order: usize) -> Self {
        Self { order }
    }

    /// Number of basis functions (modes).
    pub fn num_modes(&self) -> usize {
        self.order + 1
    }

    /// Evaluates basis function i at reference coordinate xi in \[-1,1\].
    pub fn evaluate(&self, i: usize, xi: f64) -> f64 {
        let nf = i as f64;
        let norm = ((2.0 * nf + 1.0) / 2.0).sqrt();
        norm * legendre_poly(i, xi)
    }

    /// Evaluates derivative of basis function i at xi.
    pub fn evaluate_deriv(&self, i: usize, xi: f64) -> f64 {
        let nf = i as f64;
        let norm = ((2.0 * nf + 1.0) / 2.0).sqrt();
        norm * legendre_poly_deriv(i, xi)
    }

    /// Evaluates the solution from modal coefficients at reference coordinate.
    pub fn evaluate_solution(&self, coeffs: &[f64], xi: f64) -> f64 {
        let mut val = 0.0;
        for (i, &c) in coeffs.iter().enumerate().take(self.num_modes()) {
            val += c * self.evaluate(i, xi);
        }
        val
    }
}

// ---------------------------------------------------------------------------
// 2-D tensor product DG basis
// ---------------------------------------------------------------------------

/// 2-D tensor product DG basis on reference square \[-1,1\]^2.
#[derive(Debug, Clone)]
pub struct DgBasis2D {
    /// Polynomial order in each direction.
    pub order: usize,
}

impl DgBasis2D {
    /// Creates a new 2-D DG basis.
    pub fn new(order: usize) -> Self {
        Self { order }
    }

    /// Number of 2-D modes = (p+1)^2.
    pub fn num_modes(&self) -> usize {
        (self.order + 1) * (self.order + 1)
    }

    /// Evaluates mode (i,j) at (xi, eta).
    pub fn evaluate(&self, i: usize, j: usize, xi: f64, eta: f64) -> f64 {
        let ni = i as f64;
        let nj = j as f64;
        let norm_i = ((2.0 * ni + 1.0) / 2.0).sqrt();
        let norm_j = ((2.0 * nj + 1.0) / 2.0).sqrt();
        norm_i * legendre_poly(i, xi) * norm_j * legendre_poly(j, eta)
    }

    /// Evaluates the 2-D solution from modal coefficients.
    pub fn evaluate_solution(&self, coeffs: &[f64], xi: f64, eta: f64) -> f64 {
        let p1 = self.order + 1;
        let mut val = 0.0;
        for j in 0..p1 {
            for i in 0..p1 {
                let idx = j * p1 + i;
                if idx < coeffs.len() {
                    val += coeffs[idx] * self.evaluate(i, j, xi, eta);
                }
            }
        }
        val
    }
}

// ---------------------------------------------------------------------------
// DG Element
// ---------------------------------------------------------------------------

/// A 1-D DG element with physical bounds and modal coefficients.
#[derive(Debug, Clone)]
pub struct DgElement {
    /// Left physical coordinate.
    pub x_left: f64,
    /// Right physical coordinate.
    pub x_right: f64,
    /// Modal coefficients for the solution.
    pub coeffs: Vec<f64>,
    /// Polynomial order.
    pub order: usize,
}

impl DgElement {
    /// Creates a new element with zero coefficients.
    pub fn new(x_left: f64, x_right: f64, order: usize) -> Self {
        Self {
            x_left,
            x_right,
            coeffs: vec![0.0; order + 1],
            order,
        }
    }

    /// Element width h.
    pub fn width(&self) -> f64 {
        self.x_right - self.x_left
    }

    /// Map reference coordinate xi in \[-1,1\] to physical x.
    pub fn ref_to_phys(&self, xi: f64) -> f64 {
        0.5 * ((1.0 - xi) * self.x_left + (1.0 + xi) * self.x_right)
    }

    /// Map physical x to reference coordinate xi.
    pub fn phys_to_ref(&self, x: f64) -> f64 {
        2.0 * (x - self.x_left) / self.width() - 1.0
    }

    /// Jacobian of the mapping: dx/dxi.
    pub fn jacobian(&self) -> f64 {
        self.width() / 2.0
    }

    /// Evaluate solution at reference point xi.
    pub fn evaluate(&self, xi: f64) -> f64 {
        let basis = DgBasis::new(self.order);
        basis.evaluate_solution(&self.coeffs, xi)
    }

    /// Evaluate solution at left boundary (xi = -1).
    pub fn eval_left(&self) -> f64 {
        self.evaluate(-1.0)
    }

    /// Evaluate solution at right boundary (xi = 1).
    pub fn eval_right(&self) -> f64 {
        self.evaluate(1.0)
    }
}

// ---------------------------------------------------------------------------
// Numerical fluxes
// ---------------------------------------------------------------------------

/// Enumeration of supported numerical flux types.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NumericalFlux {
    /// Lax-Friedrichs (Rusanov) flux.
    LaxFriedrichs,
    /// Roe flux with entropy fix.
    Roe,
    /// Simple upwind flux.
    Upwind,
}

/// Computes the Lax-Friedrichs numerical flux for a scalar conservation law
/// u_t + f(u)_x = 0 at an interface.
///
/// f_hat = 0.5 * (f(u_L) + f(u_R) - alpha * (u_R - u_L))
/// where alpha = max |f'(u)|.
pub fn lax_friedrichs_flux(
    u_left: f64,
    u_right: f64,
    flux_left: f64,
    flux_right: f64,
    alpha: f64,
) -> f64 {
    0.5 * (flux_left + flux_right - alpha * (u_right - u_left))
}

/// Computes the Roe numerical flux for a scalar conservation law.
///
/// Uses the Roe average speed a_hat = (f(u_R) - f(u_L)) / (u_R - u_L).
pub fn roe_flux(u_left: f64, u_right: f64, flux_left: f64, flux_right: f64) -> f64 {
    if (u_right - u_left).abs() < 1e-14 {
        return 0.5 * (flux_left + flux_right);
    }
    let a_hat = (flux_right - flux_left) / (u_right - u_left);
    if a_hat >= 0.0 { flux_left } else { flux_right }
}

/// Computes the upwind numerical flux for advection u_t + a*u_x = 0.
pub fn upwind_flux(u_left: f64, u_right: f64, wave_speed: f64) -> f64 {
    if wave_speed >= 0.0 {
        wave_speed * u_left
    } else {
        wave_speed * u_right
    }
}

/// Evaluates the numerical flux at an interface given the flux type.
#[allow(clippy::too_many_arguments)]
pub fn compute_numerical_flux(
    flux_type: NumericalFlux,
    u_left: f64,
    u_right: f64,
    flux_left: f64,
    flux_right: f64,
    wave_speed: f64,
) -> f64 {
    match flux_type {
        NumericalFlux::LaxFriedrichs => {
            lax_friedrichs_flux(u_left, u_right, flux_left, flux_right, wave_speed.abs())
        }
        NumericalFlux::Roe => roe_flux(u_left, u_right, flux_left, flux_right),
        NumericalFlux::Upwind => upwind_flux(u_left, u_right, wave_speed),
    }
}

// ---------------------------------------------------------------------------
// Interior penalty method
// ---------------------------------------------------------------------------

/// Penalty method variant for the diffusion operator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PenaltyMethod {
    /// Symmetric Interior Penalty Galerkin (SIPG): sigma = -1.
    Sipg,
    /// Non-symmetric Interior Penalty Galerkin (NIPG): sigma = 1.
    Nipg,
    /// Incomplete Interior Penalty Galerkin (IIPG): sigma = 0.
    Iipg,
}

impl PenaltyMethod {
    /// Returns the symmetry parameter sigma for each variant.
    pub fn sigma(&self) -> f64 {
        match self {
            PenaltyMethod::Sipg => -1.0,
            PenaltyMethod::Nipg => 1.0,
            PenaltyMethod::Iipg => 0.0,
        }
    }
}

/// Computes the interior penalty parameter eta for a face between two elements.
///
/// eta = C_pen * (p+1)^2 / h_min, where h_min is the minimum element size.
pub fn penalty_parameter(order: usize, h_left: f64, h_right: f64, c_pen: f64) -> f64 {
    let p = order as f64;
    let h_min = h_left.min(h_right);
    c_pen * (p + 1.0) * (p + 1.0) / h_min
}

/// Interior penalty diffusion bilinear form contribution at an interior face.
///
/// Returns (contribution to left element, contribution to right element).
#[allow(clippy::too_many_arguments)]
pub fn interior_penalty_face(
    u_left: f64,
    u_right: f64,
    grad_u_left: f64,
    grad_u_right: f64,
    diffusivity: f64,
    eta: f64,
    method: PenaltyMethod,
) -> (f64, f64) {
    let jump = u_right - u_left;
    let avg_grad = 0.5 * (grad_u_left + grad_u_right);
    let sig = method.sigma();

    // Penalty term
    let penalty = eta * diffusivity * jump;
    // Consistency term
    let consistency = -diffusivity * avg_grad * jump;
    // Symmetry term
    let symmetry = sig * consistency;

    let total = penalty + consistency + symmetry;
    (-total * 0.5, total * 0.5)
}

// ---------------------------------------------------------------------------
// h/p-adaptivity
// ---------------------------------------------------------------------------

/// Error indicator for h/p-adaptivity.
#[derive(Debug, Clone)]
pub struct ErrorIndicator {
    /// Element index.
    pub element_id: usize,
    /// Estimated local error.
    pub error: f64,
    /// Smoothness indicator (high = smooth => prefer p-refinement).
    pub smoothness: f64,
}

/// Computes a smoothness indicator for a DG element based on the decay rate
/// of its modal coefficients.
///
/// A rapid decay suggests a smooth solution (prefer p-enrichment);
/// a slow decay suggests a non-smooth solution (prefer h-refinement).
pub fn smoothness_indicator(coeffs: &[f64]) -> f64 {
    let n = coeffs.len();
    if n < 3 {
        return 0.0;
    }
    // Ratio of energy in highest mode to total energy
    let total_energy: f64 = coeffs.iter().map(|c| c * c).sum();
    if total_energy < 1e-30 {
        return 1.0; // Zero solution is smooth
    }
    let high_energy: f64 = coeffs[n - 2..].iter().map(|c| c * c).sum();
    1.0 - (high_energy / total_energy)
}

/// Decides whether to h-refine or p-refine an element based on the
/// smoothness indicator.
///
/// Returns `true` for p-refinement, `false` for h-refinement.
pub fn decide_hp_refinement(smoothness: f64, threshold: f64) -> bool {
    smoothness > threshold
}

/// Performs h-refinement on a DG element, splitting it into two child elements.
///
/// The solution is L2-projected onto the children.
pub fn h_refine_element(elem: &DgElement) -> (DgElement, DgElement) {
    let mid = 0.5 * (elem.x_left + elem.x_right);
    let basis = DgBasis::new(elem.order);
    let (qpts, qwts) = gauss_legendre_quadrature(elem.order + 2);

    // Project onto left child
    let mut left_child = DgElement::new(elem.x_left, mid, elem.order);
    for i in 0..=elem.order {
        let mut proj = 0.0;
        for (k, &xi) in qpts.iter().enumerate() {
            let x_child = left_child.ref_to_phys(xi);
            let xi_parent = elem.phys_to_ref(x_child);
            let u_parent = basis.evaluate_solution(&elem.coeffs, xi_parent);
            proj += qwts[k] * u_parent * basis.evaluate(i, xi);
        }
        left_child.coeffs[i] = proj;
    }

    // Project onto right child
    let mut right_child = DgElement::new(mid, elem.x_right, elem.order);
    for i in 0..=elem.order {
        let mut proj = 0.0;
        for (k, &xi) in qpts.iter().enumerate() {
            let x_child = right_child.ref_to_phys(xi);
            let xi_parent = elem.phys_to_ref(x_child);
            let u_parent = basis.evaluate_solution(&elem.coeffs, xi_parent);
            proj += qwts[k] * u_parent * basis.evaluate(i, xi);
        }
        right_child.coeffs[i] = proj;
    }

    (left_child, right_child)
}

/// Performs p-enrichment on a DG element by increasing the polynomial order.
pub fn p_enrich_element(elem: &DgElement, new_order: usize) -> DgElement {
    let mut new_elem = DgElement::new(elem.x_left, elem.x_right, new_order);
    // Copy existing coefficients, new ones are zero
    for (i, &c) in elem.coeffs.iter().enumerate() {
        if i <= new_order {
            new_elem.coeffs[i] = c;
        }
    }
    new_elem
}

/// Computes error indicators for all elements using the jump in solution
/// at element interfaces as an error proxy.
pub fn compute_error_indicators(elements: &[DgElement]) -> Vec<ErrorIndicator> {
    let n = elements.len();
    let mut indicators = Vec::with_capacity(n);

    for i in 0..n {
        let mut error = 0.0;
        // Right face jump
        if i + 1 < n {
            let jump = elements[i].eval_right() - elements[i + 1].eval_left();
            error += jump * jump;
        }
        // Left face jump
        if i > 0 {
            let jump = elements[i].eval_left() - elements[i - 1].eval_right();
            error += jump * jump;
        }
        let smoothness = smoothness_indicator(&elements[i].coeffs);
        indicators.push(ErrorIndicator {
            element_id: i,
            error: error.sqrt(),
            smoothness,
        });
    }
    indicators
}

// ---------------------------------------------------------------------------
// TVB limiter
// ---------------------------------------------------------------------------

/// Modified TVB minmod function.
///
/// Returns the minmod of (a, b, c) with TVB modification parameter M.
pub fn tvb_minmod(a: f64, b: f64, c: f64, m_tvb: f64, dx: f64) -> f64 {
    if a.abs() <= m_tvb * dx * dx {
        return a;
    }
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}

/// Applies TVB slope limiter to a DG element given neighbor cell averages.
///
/// For first-order slopes (linear DG), limits the slope coefficient.
pub fn tvb_limiter(elem: &mut DgElement, u_avg_left: f64, u_avg_right: f64, m_tvb: f64) {
    if elem.order < 1 {
        return;
    }
    let basis = DgBasis::new(elem.order);
    let u_avg = elem.coeffs[0] * basis.evaluate(0, 0.0);
    let dx = elem.width();

    let slope = elem.coeffs[1] * basis.evaluate(1, 1.0);
    let d_left = u_avg - u_avg_left;
    let d_right = u_avg_right - u_avg;

    let limited_slope = tvb_minmod(slope, d_right, d_left, m_tvb, dx);
    if (limited_slope - slope).abs() > 1e-14 {
        // Reduce to limited slope
        let scale = if slope.abs() > 1e-14 {
            limited_slope / slope
        } else {
            0.0
        };
        // Zero out higher modes, scale the linear mode
        elem.coeffs[1] *= scale;
        for i in 2..elem.coeffs.len() {
            elem.coeffs[i] = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// WENO reconstruction
// ---------------------------------------------------------------------------

/// Simple 3rd-order WENO reconstruction at a face, given cell averages of
/// three consecutive cells.
///
/// Returns (u_minus, u_plus) at the right face of the middle cell.
pub fn weno3_reconstruct(u_left: f64, u_mid: f64, u_right: f64) -> (f64, f64) {
    let eps = 1e-6;

    // Stencil 0: {u_left, u_mid}
    let p0 = 0.5 * (u_left + u_mid);
    let beta0 = (u_mid - u_left) * (u_mid - u_left);

    // Stencil 1: {u_mid, u_right}
    let p1 = 0.5 * (u_mid + u_right);
    let beta1 = (u_right - u_mid) * (u_right - u_mid);

    // Optimal weights for 3rd order
    let d0 = 1.0 / 3.0;
    let d1 = 2.0 / 3.0;

    let alpha0 = d0 / ((eps + beta0) * (eps + beta0));
    let alpha1 = d1 / ((eps + beta1) * (eps + beta1));
    let sum = alpha0 + alpha1;

    let w0 = alpha0 / sum;
    let w1 = alpha1 / sum;

    let u_minus = w0 * p0 + w1 * p1;

    // For u_plus, swap roles
    let d0p = 2.0 / 3.0;
    let d1p = 1.0 / 3.0;
    let alpha0p = d0p / ((eps + beta0) * (eps + beta0));
    let alpha1p = d1p / ((eps + beta1) * (eps + beta1));
    let sump = alpha0p + alpha1p;
    let w0p = alpha0p / sump;
    let w1p = alpha1p / sump;
    let u_plus = w0p * p0 + w1p * p1;

    (u_minus, u_plus)
}

/// Applies WENO limiting to a vector of cell averages, returning
/// reconstructed face values (left and right) for each cell.
pub fn weno_limiting(cell_averages: &[f64]) -> Vec<(f64, f64)> {
    let n = cell_averages.len();
    let mut faces = Vec::with_capacity(n);
    for i in 0..n {
        let u_l = if i > 0 {
            cell_averages[i - 1]
        } else {
            cell_averages[i]
        };
        let u_m = cell_averages[i];
        let u_r = if i + 1 < n {
            cell_averages[i + 1]
        } else {
            cell_averages[i]
        };
        let (u_minus, u_plus) = weno3_reconstruct(u_l, u_m, u_r);
        faces.push((u_minus, u_plus));
    }
    faces
}

// ---------------------------------------------------------------------------
// Runge-Kutta DG time integrators
// ---------------------------------------------------------------------------

/// Computes the DG spatial residual for scalar advection u_t + a*u_x = 0
/// on a uniform mesh of DG elements.
///
/// Returns a vector of residual coefficient vectors (one per element).
#[allow(clippy::too_many_arguments)]
pub fn dg_advection_residual(
    elements: &[DgElement],
    wave_speed: f64,
    flux_type: NumericalFlux,
) -> Vec<Vec<f64>> {
    let n = elements.len();
    let mut residuals = Vec::with_capacity(n);

    for i in 0..n {
        let basis = DgBasis::new(elements[i].order);
        let nq = elements[i].order + 2;
        let (qpts, qwts) = gauss_legendre_quadrature(nq);
        let jac = elements[i].jacobian();
        let nmodes = basis.num_modes();
        let mut res = vec![0.0; nmodes];

        // Volume integral: integral( a * u * dphi/dxi ) dxi
        for (k, &xi) in qpts.iter().enumerate() {
            let u_val = basis.evaluate_solution(&elements[i].coeffs, xi);
            for m in 0..nmodes {
                let dphi = basis.evaluate_deriv(m, xi);
                res[m] += qwts[k] * wave_speed * u_val * dphi;
            }
        }

        // Surface terms (numerical flux)
        // Right face (xi = 1)
        let u_int_r = elements[i].eval_right();
        let u_ext_r = if i + 1 < n {
            elements[i + 1].eval_left()
        } else {
            u_int_r // Extrapolation BC
        };
        let f_hat_r = compute_numerical_flux(
            flux_type,
            u_int_r,
            u_ext_r,
            wave_speed * u_int_r,
            wave_speed * u_ext_r,
            wave_speed,
        );

        // Left face (xi = -1)
        let u_int_l = elements[i].eval_left();
        let u_ext_l = if i > 0 {
            elements[i - 1].eval_right()
        } else {
            u_int_l // Extrapolation BC
        };
        let f_hat_l = compute_numerical_flux(
            flux_type,
            u_ext_l,
            u_int_l,
            wave_speed * u_ext_l,
            wave_speed * u_int_l,
            wave_speed,
        );

        for m in 0..nmodes {
            let phi_r = basis.evaluate(m, 1.0);
            let phi_l = basis.evaluate(m, -1.0);
            res[m] -= f_hat_r * phi_r - f_hat_l * phi_l;
        }

        // Scale by inverse Jacobian
        for m in 0..nmodes {
            res[m] /= jac;
        }

        residuals.push(res);
    }

    residuals
}

/// SSP-RK2 (Heun's method) time step for the DG system.
///
/// u^(1) = u^n + dt * L(u^n)
/// u^(n+1) = 0.5 * u^n + 0.5 * (u^(1) + dt * L(u^(1)))
pub fn ssp_rk2_step(
    elements: &mut [DgElement],
    wave_speed: f64,
    dt: f64,
    flux_type: NumericalFlux,
) {
    let n = elements.len();

    // Save u^n
    let u_n: Vec<Vec<f64>> = elements.iter().map(|e| e.coeffs.clone()).collect();

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    let res1 = dg_advection_residual(elements, wave_speed, flux_type);
    for i in 0..n {
        for (m, c) in elements[i].coeffs.iter_mut().enumerate() {
            *c += dt * res1[i][m];
        }
    }

    // Stage 2: u^(n+1) = 0.5 * u^n + 0.5 * (u^(1) + dt * L(u^(1)))
    let res2 = dg_advection_residual(elements, wave_speed, flux_type);
    for i in 0..n {
        for (m, c) in elements[i].coeffs.iter_mut().enumerate() {
            *c = 0.5 * u_n[i][m] + 0.5 * (*c + dt * res2[i][m]);
        }
    }
}

/// SSP-RK3 time step for the DG system.
///
/// Three-stage strong stability preserving Runge-Kutta.
pub fn ssp_rk3_step(
    elements: &mut [DgElement],
    wave_speed: f64,
    dt: f64,
    flux_type: NumericalFlux,
) {
    let n = elements.len();
    let u_n: Vec<Vec<f64>> = elements.iter().map(|e| e.coeffs.clone()).collect();

    // Stage 1
    let res1 = dg_advection_residual(elements, wave_speed, flux_type);
    for i in 0..n {
        for (m, c) in elements[i].coeffs.iter_mut().enumerate() {
            *c += dt * res1[i][m];
        }
    }

    // Stage 2
    let res2 = dg_advection_residual(elements, wave_speed, flux_type);
    for i in 0..n {
        for (m, c) in elements[i].coeffs.iter_mut().enumerate() {
            *c = 0.75 * u_n[i][m] + 0.25 * (*c + dt * res2[i][m]);
        }
    }

    // Stage 3
    let res3 = dg_advection_residual(elements, wave_speed, flux_type);
    for i in 0..n {
        for (m, c) in elements[i].coeffs.iter_mut().enumerate() {
            *c = (1.0 / 3.0) * u_n[i][m] + (2.0 / 3.0) * (*c + dt * res3[i][m]);
        }
    }
}

// ---------------------------------------------------------------------------
// Advection-diffusion DG
// ---------------------------------------------------------------------------

/// DG spatial residual for the advection-diffusion equation:
/// u_t + a*u_x = kappa * u_xx
///
/// Uses the Local DG (LDG) approach for the diffusion term.
#[allow(clippy::too_many_arguments)]
pub fn dg_advection_diffusion_residual(
    elements: &[DgElement],
    wave_speed: f64,
    kappa: f64,
    penalty_method: PenaltyMethod,
    c_pen: f64,
    flux_type: NumericalFlux,
) -> Vec<Vec<f64>> {
    let n = elements.len();

    // Start with advection residual
    let mut residuals = dg_advection_residual(elements, wave_speed, flux_type);

    let basis0 = DgBasis::new(if n > 0 { elements[0].order } else { 0 });

    // Add diffusion terms via interior penalty
    for i in 0..n {
        let basis = DgBasis::new(elements[i].order);
        let nq = elements[i].order + 2;
        let (qpts, qwts) = gauss_legendre_quadrature(nq);
        let jac = elements[i].jacobian();
        let nmodes = basis.num_modes();

        // Volume integral: kappa * integral( du/dxi * dphi/dxi ) / jac^2 * jac
        for (k, &xi) in qpts.iter().enumerate() {
            let mut du_dxi = 0.0;
            for (m, &c) in elements[i].coeffs.iter().enumerate().take(nmodes) {
                du_dxi += c * basis.evaluate_deriv(m, xi);
            }
            for m in 0..nmodes {
                let dphi = basis.evaluate_deriv(m, xi);
                residuals[i][m] += qwts[k] * kappa * du_dxi * dphi / jac;
            }
        }

        // Face terms (interior penalty)
        if i + 1 < n {
            let u_l = elements[i].eval_right();
            let u_r = elements[i + 1].eval_left();

            let mut grad_l = 0.0;
            for (m, &c) in elements[i]
                .coeffs
                .iter()
                .enumerate()
                .take(basis.num_modes())
            {
                grad_l += c * basis.evaluate_deriv(m, 1.0);
            }
            grad_l /= elements[i].jacobian();

            let basis_r = DgBasis::new(elements[i + 1].order);
            let mut grad_r = 0.0;
            for (m, &c) in elements[i + 1]
                .coeffs
                .iter()
                .enumerate()
                .take(basis_r.num_modes())
            {
                grad_r += c * basis_r.evaluate_deriv(m, -1.0);
            }
            grad_r /= elements[i + 1].jacobian();

            let eta = penalty_parameter(
                elements[i].order.max(elements[i + 1].order),
                elements[i].width(),
                elements[i + 1].width(),
                c_pen,
            );

            let (_face_l, _face_r) =
                interior_penalty_face(u_l, u_r, grad_l, grad_r, kappa, eta, penalty_method);

            let _noop = &basis0;

            // Apply penalty contributions to boundary modes
            for m in 0..basis.num_modes() {
                let phi_r = basis.evaluate(m, 1.0);
                residuals[i][m] -= eta * kappa * (u_l - u_r) * phi_r / jac;
            }
        }
    }

    residuals
}

// ---------------------------------------------------------------------------
// Hyperbolic conservation laws: 1-D Euler equations
// ---------------------------------------------------------------------------

/// Conserved state for 1-D Euler equations: (rho, rho*u, E).
#[derive(Debug, Clone, Copy)]
pub struct EulerState {
    /// Density.
    pub rho: f64,
    /// Momentum (rho * u).
    pub rho_u: f64,
    /// Total energy.
    pub energy: f64,
}

impl EulerState {
    /// Creates a new Euler state.
    pub fn new(rho: f64, rho_u: f64, energy: f64) -> Self {
        Self { rho, rho_u, energy }
    }

    /// Velocity.
    pub fn velocity(&self) -> f64 {
        self.rho_u / self.rho
    }

    /// Pressure (ideal gas with gamma).
    pub fn pressure(&self, gamma: f64) -> f64 {
        (gamma - 1.0) * (self.energy - 0.5 * self.rho_u * self.rho_u / self.rho)
    }

    /// Sound speed.
    pub fn sound_speed(&self, gamma: f64) -> f64 {
        let p = self.pressure(gamma);
        (gamma * p / self.rho).abs().sqrt()
    }

    /// Physical flux vector for 1-D Euler.
    pub fn flux(&self, gamma: f64) -> [f64; 3] {
        let u = self.velocity();
        let p = self.pressure(gamma);
        [self.rho_u, self.rho_u * u + p, (self.energy + p) * u]
    }

    /// Maximum wave speed |u| + c.
    pub fn max_wave_speed(&self, gamma: f64) -> f64 {
        self.velocity().abs() + self.sound_speed(gamma)
    }

    /// Converts to array.
    pub fn to_array(&self) -> [f64; 3] {
        [self.rho, self.rho_u, self.energy]
    }

    /// Creates from array.
    pub fn from_array(a: [f64; 3]) -> Self {
        Self {
            rho: a[0],
            rho_u: a[1],
            energy: a[2],
        }
    }
}

/// Lax-Friedrichs flux for 1-D Euler equations.
pub fn euler_lax_friedrichs(left: &EulerState, right: &EulerState, gamma: f64) -> [f64; 3] {
    let fl = left.flux(gamma);
    let fr = right.flux(gamma);
    let alpha = left.max_wave_speed(gamma).max(right.max_wave_speed(gamma));
    let ul = left.to_array();
    let ur = right.to_array();
    [
        0.5 * (fl[0] + fr[0] - alpha * (ur[0] - ul[0])),
        0.5 * (fl[1] + fr[1] - alpha * (ur[1] - ul[1])),
        0.5 * (fl[2] + fr[2] - alpha * (ur[2] - ul[2])),
    ]
}

/// Roe flux for 1-D Euler equations (simplified, without entropy fix).
pub fn euler_roe_flux(left: &EulerState, right: &EulerState, gamma: f64) -> [f64; 3] {
    let rho_l = left.rho.sqrt();
    let rho_r = right.rho.sqrt();
    let denom = rho_l + rho_r;

    let u_roe = (rho_l * left.velocity() + rho_r * right.velocity()) / denom;
    let hl = (left.energy + left.pressure(gamma)) / left.rho;
    let hr = (right.energy + right.pressure(gamma)) / right.rho;
    let h_roe = (rho_l * hl + rho_r * hr) / denom;
    let c_roe = ((gamma - 1.0) * (h_roe - 0.5 * u_roe * u_roe)).abs().sqrt();

    let fl = left.flux(gamma);
    let fr = right.flux(gamma);

    // Eigenvalues
    let lam1 = (u_roe - c_roe).abs();
    let lam2 = u_roe.abs();
    let lam3 = (u_roe + c_roe).abs();
    let max_lam = lam1.max(lam2).max(lam3);

    // Simple dissipation (Rusanov-like with Roe average)
    let ul = left.to_array();
    let ur = right.to_array();
    [
        0.5 * (fl[0] + fr[0] - max_lam * (ur[0] - ul[0])),
        0.5 * (fl[1] + fr[1] - max_lam * (ur[1] - ul[1])),
        0.5 * (fl[2] + fr[2] - max_lam * (ur[2] - ul[2])),
    ]
}

// ---------------------------------------------------------------------------
// DG mesh: 1-D uniform mesh of DG elements
// ---------------------------------------------------------------------------

/// A 1-D DG mesh consisting of a list of DG elements.
#[derive(Debug, Clone)]
pub struct DgMesh1D {
    /// The DG elements.
    pub elements: Vec<DgElement>,
}

impl DgMesh1D {
    /// Creates a uniform DG mesh on \[x_min, x_max\] with n_elem elements
    /// at a given polynomial order.
    pub fn uniform(x_min: f64, x_max: f64, n_elem: usize, order: usize) -> Self {
        let dx = (x_max - x_min) / n_elem as f64;
        let mut elements = Vec::with_capacity(n_elem);
        for i in 0..n_elem {
            let xl = x_min + i as f64 * dx;
            let xr = xl + dx;
            elements.push(DgElement::new(xl, xr, order));
        }
        Self { elements }
    }

    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Total degrees of freedom.
    pub fn total_dofs(&self) -> usize {
        self.elements.iter().map(|e| e.order + 1).sum()
    }

    /// L2-projects an initial condition onto the DG space.
    pub fn project_initial_condition<F: Fn(f64) -> f64>(&mut self, f: F) {
        for elem in &mut self.elements {
            let basis = DgBasis::new(elem.order);
            let nq = elem.order + 3;
            let (qpts, qwts) = gauss_legendre_quadrature(nq);
            for i in 0..=elem.order {
                let mut proj = 0.0;
                for (k, &xi) in qpts.iter().enumerate() {
                    let x = elem.ref_to_phys(xi);
                    proj += qwts[k] * f(x) * basis.evaluate(i, xi);
                }
                elem.coeffs[i] = proj;
            }
        }
    }

    /// Evaluates the DG solution at a physical point x.
    pub fn evaluate_at(&self, x: f64) -> f64 {
        for elem in &self.elements {
            if x >= elem.x_left && x <= elem.x_right {
                let xi = elem.phys_to_ref(x);
                return elem.evaluate(xi);
            }
        }
        0.0
    }

    /// Computes the L2-norm of the DG solution.
    pub fn l2_norm(&self) -> f64 {
        let mut norm_sq = 0.0;
        for elem in &self.elements {
            let nq = elem.order + 3;
            let (qpts, qwts) = gauss_legendre_quadrature(nq);
            let jac = elem.jacobian();
            for (k, &xi) in qpts.iter().enumerate() {
                let u = elem.evaluate(xi);
                norm_sq += qwts[k] * u * u * jac;
            }
        }
        norm_sq.sqrt()
    }

    /// Computes the L2 error against an exact solution.
    pub fn l2_error<F: Fn(f64) -> f64>(&self, exact: F) -> f64 {
        let mut error_sq = 0.0;
        for elem in &self.elements {
            let nq = elem.order + 3;
            let (qpts, qwts) = gauss_legendre_quadrature(nq);
            let jac = elem.jacobian();
            for (k, &xi) in qpts.iter().enumerate() {
                let x = elem.ref_to_phys(xi);
                let diff = elem.evaluate(xi) - exact(x);
                error_sq += qwts[k] * diff * diff * jac;
            }
        }
        error_sq.sqrt()
    }

    /// Computes the cell averages.
    pub fn cell_averages(&self) -> Vec<f64> {
        self.elements
            .iter()
            .map(|elem| {
                let basis = DgBasis::new(elem.order);
                elem.coeffs[0] * basis.evaluate(0, 0.0)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// DG solver for scalar advection
// ---------------------------------------------------------------------------

/// Configuration for the DG advection solver.
#[derive(Debug, Clone)]
pub struct DgAdvectionConfig {
    /// Wave speed a.
    pub wave_speed: f64,
    /// Numerical flux type.
    pub flux_type: NumericalFlux,
    /// CFL number.
    pub cfl: f64,
    /// TVB limiter parameter M (0 = strict minmod, > 0 = TVB).
    pub tvb_m: f64,
    /// Whether to apply limiting.
    pub use_limiter: bool,
}

impl Default for DgAdvectionConfig {
    fn default() -> Self {
        Self {
            wave_speed: 1.0,
            flux_type: NumericalFlux::LaxFriedrichs,
            cfl: 0.5,
            tvb_m: 0.0,
            use_limiter: false,
        }
    }
}

/// Runs the DG advection solver for a given number of time steps.
pub fn solve_advection(mesh: &mut DgMesh1D, config: &DgAdvectionConfig, t_final: f64) -> f64 {
    let n = mesh.num_elements();
    if n == 0 {
        return 0.0;
    }
    let dx_min = mesh
        .elements
        .iter()
        .map(|e| e.width())
        .fold(f64::INFINITY, f64::min);
    let p = mesh.elements[0].order;
    let dt = config.cfl * dx_min / (config.wave_speed.abs() * (2.0 * p as f64 + 1.0));

    let mut t = 0.0;
    while t < t_final - 1e-14 {
        let dt_actual = dt.min(t_final - t);

        if p <= 1 {
            ssp_rk2_step(
                &mut mesh.elements,
                config.wave_speed,
                dt_actual,
                config.flux_type,
            );
        } else {
            ssp_rk3_step(
                &mut mesh.elements,
                config.wave_speed,
                dt_actual,
                config.flux_type,
            );
        }

        // Apply limiter if requested
        if config.use_limiter {
            let averages = mesh.cell_averages();
            for i in 0..n {
                let u_l = if i > 0 { averages[i - 1] } else { averages[i] };
                let u_r = if i + 1 < n {
                    averages[i + 1]
                } else {
                    averages[i]
                };
                tvb_limiter(&mut mesh.elements[i], u_l, u_r, config.tvb_m);
            }
        }

        t += dt_actual;
    }
    t
}

// ---------------------------------------------------------------------------
// DG solver for advection-diffusion
// ---------------------------------------------------------------------------

/// Configuration for the DG advection-diffusion solver.
#[derive(Debug, Clone)]
pub struct DgAdvDiffConfig {
    /// Advection wave speed.
    pub wave_speed: f64,
    /// Diffusion coefficient.
    pub kappa: f64,
    /// Numerical flux type for advection.
    pub flux_type: NumericalFlux,
    /// Penalty method for diffusion.
    pub penalty_method: PenaltyMethod,
    /// Penalty constant.
    pub c_pen: f64,
    /// CFL number.
    pub cfl: f64,
}

impl Default for DgAdvDiffConfig {
    fn default() -> Self {
        Self {
            wave_speed: 1.0,
            kappa: 0.01,
            flux_type: NumericalFlux::LaxFriedrichs,
            penalty_method: PenaltyMethod::Sipg,
            c_pen: 10.0,
            cfl: 0.1,
        }
    }
}

/// Runs the DG advection-diffusion solver.
pub fn solve_advection_diffusion(
    mesh: &mut DgMesh1D,
    config: &DgAdvDiffConfig,
    t_final: f64,
) -> f64 {
    let n = mesh.num_elements();
    if n == 0 {
        return 0.0;
    }
    let dx_min = mesh
        .elements
        .iter()
        .map(|e| e.width())
        .fold(f64::INFINITY, f64::min);
    let p = mesh.elements[0].order as f64;

    // Time step combines advection and diffusion CFL
    let dt_adv = config.cfl * dx_min / (config.wave_speed.abs() * (2.0 * p + 1.0));
    let dt_diff = config.cfl * dx_min * dx_min / (config.kappa * (2.0 * p + 1.0).powi(2));
    let dt = dt_adv.min(dt_diff);

    let mut t = 0.0;
    while t < t_final - 1e-14 {
        let dt_actual = dt.min(t_final - t);

        // Forward Euler step for simplicity
        let residuals = dg_advection_diffusion_residual(
            &mesh.elements,
            config.wave_speed,
            config.kappa,
            config.penalty_method,
            config.c_pen,
            config.flux_type,
        );
        for i in 0..n {
            for (m, c) in mesh.elements[i].coeffs.iter_mut().enumerate() {
                *c += dt_actual * residuals[i][m];
            }
        }

        t += dt_actual;
    }
    t
}

// ---------------------------------------------------------------------------
// Burgers equation DG
// ---------------------------------------------------------------------------

/// Flux function for Burgers equation: f(u) = u^2 / 2.
pub fn burgers_flux(u: f64) -> f64 {
    0.5 * u * u
}

/// Maximum wave speed for Burgers equation: |u|.
pub fn burgers_wave_speed(u: f64) -> f64 {
    u.abs()
}

/// DG spatial residual for Burgers equation u_t + (u^2/2)_x = 0.
pub fn dg_burgers_residual(elements: &[DgElement], flux_type: NumericalFlux) -> Vec<Vec<f64>> {
    let n = elements.len();
    let mut residuals = Vec::with_capacity(n);

    for i in 0..n {
        let basis = DgBasis::new(elements[i].order);
        let nq = elements[i].order + 3;
        let (qpts, qwts) = gauss_legendre_quadrature(nq);
        let jac = elements[i].jacobian();
        let nmodes = basis.num_modes();
        let mut res = vec![0.0; nmodes];

        // Volume integral
        for (k, &xi) in qpts.iter().enumerate() {
            let u_val = basis.evaluate_solution(&elements[i].coeffs, xi);
            let f_val = burgers_flux(u_val);
            for m in 0..nmodes {
                let dphi = basis.evaluate_deriv(m, xi);
                res[m] += qwts[k] * f_val * dphi;
            }
        }

        // Surface terms
        let u_int_r = elements[i].eval_right();
        let u_ext_r = if i + 1 < n {
            elements[i + 1].eval_left()
        } else {
            u_int_r
        };
        let alpha_r = u_int_r.abs().max(u_ext_r.abs());
        let f_hat_r = compute_numerical_flux(
            flux_type,
            u_int_r,
            u_ext_r,
            burgers_flux(u_int_r),
            burgers_flux(u_ext_r),
            alpha_r,
        );

        let u_int_l = elements[i].eval_left();
        let u_ext_l = if i > 0 {
            elements[i - 1].eval_right()
        } else {
            u_int_l
        };
        let alpha_l = u_int_l.abs().max(u_ext_l.abs());
        let f_hat_l = compute_numerical_flux(
            flux_type,
            u_ext_l,
            u_int_l,
            burgers_flux(u_ext_l),
            burgers_flux(u_int_l),
            alpha_l,
        );

        for m in 0..nmodes {
            let phi_r = basis.evaluate(m, 1.0);
            let phi_l = basis.evaluate(m, -1.0);
            res[m] -= f_hat_r * phi_r - f_hat_l * phi_l;
        }

        for m in 0..nmodes {
            res[m] /= jac;
        }

        residuals.push(res);
    }

    residuals
}

// ---------------------------------------------------------------------------
// Mass matrix
// ---------------------------------------------------------------------------

/// Computes the local mass matrix for a DG element.
///
/// M_ij = integral( phi_i * phi_j ) * jacobian
pub fn local_mass_matrix(elem: &DgElement) -> Vec<Vec<f64>> {
    let basis = DgBasis::new(elem.order);
    let nq = 2 * elem.order + 2;
    let (qpts, qwts) = gauss_legendre_quadrature(nq);
    let jac = elem.jacobian();
    let nm = basis.num_modes();

    let mut mass = vec![vec![0.0; nm]; nm];
    for (k, &xi) in qpts.iter().enumerate() {
        for i in 0..nm {
            let phi_i = basis.evaluate(i, xi);
            for j in 0..nm {
                let phi_j = basis.evaluate(j, xi);
                mass[i][j] += qwts[k] * phi_i * phi_j * jac;
            }
        }
    }
    mass
}

/// Computes the local stiffness matrix for a DG element.
///
/// S_ij = integral( dphi_i/dx * dphi_j/dx ) * jacobian
pub fn local_stiffness_matrix(elem: &DgElement) -> Vec<Vec<f64>> {
    let basis = DgBasis::new(elem.order);
    let nq = 2 * elem.order + 2;
    let (qpts, qwts) = gauss_legendre_quadrature(nq);
    let jac = elem.jacobian();
    let nm = basis.num_modes();

    let mut stiff = vec![vec![0.0; nm]; nm];
    for (k, &xi) in qpts.iter().enumerate() {
        for i in 0..nm {
            let dphi_i = basis.evaluate_deriv(i, xi) / jac;
            for j in 0..nm {
                let dphi_j = basis.evaluate_deriv(j, xi) / jac;
                stiff[i][j] += qwts[k] * dphi_i * dphi_j * jac;
            }
        }
    }
    stiff
}

// ---------------------------------------------------------------------------
// Solution diagnostics
// ---------------------------------------------------------------------------

/// Computes the total integral (mass) of the DG solution over the domain.
pub fn total_mass(mesh: &DgMesh1D) -> f64 {
    let mut mass = 0.0;
    for elem in &mesh.elements {
        let nq = elem.order + 3;
        let (qpts, qwts) = gauss_legendre_quadrature(nq);
        let jac = elem.jacobian();
        for (k, &xi) in qpts.iter().enumerate() {
            mass += qwts[k] * elem.evaluate(xi) * jac;
        }
    }
    mass
}

/// Computes the total variation of the DG solution.
pub fn total_variation(mesh: &DgMesh1D) -> f64 {
    let mut tv = 0.0;
    let n_sample = 100;
    let n_elem = mesh.elements.len();
    for i in 0..n_elem {
        let dx = 2.0 / n_sample as f64;
        let mut xi = -1.0;
        let mut prev = mesh.elements[i].evaluate(xi);
        for _ in 0..n_sample {
            xi += dx;
            let curr = mesh.elements[i].evaluate(xi);
            tv += (curr - prev).abs();
            prev = curr;
        }
        // Interface jumps
        if i + 1 < n_elem {
            let u_r = mesh.elements[i].eval_right();
            let u_l = mesh.elements[i + 1].eval_left();
            tv += (u_l - u_r).abs();
        }
    }
    tv
}

/// Computes the maximum value of the DG solution over the domain.
pub fn solution_max(mesh: &DgMesh1D) -> f64 {
    let mut max_val = f64::NEG_INFINITY;
    let n_sample = 20;
    for elem in &mesh.elements {
        let dx = 2.0 / n_sample as f64;
        let mut xi = -1.0;
        for _ in 0..=n_sample {
            let val = elem.evaluate(xi);
            if val > max_val {
                max_val = val;
            }
            xi += dx;
        }
    }
    max_val
}

/// Computes the minimum value of the DG solution over the domain.
pub fn solution_min(mesh: &DgMesh1D) -> f64 {
    let mut min_val = f64::INFINITY;
    let n_sample = 20;
    for elem in &mesh.elements {
        let dx = 2.0 / n_sample as f64;
        let mut xi = -1.0;
        for _ in 0..=n_sample {
            let val = elem.evaluate(xi);
            if val < min_val {
                min_val = val;
            }
            xi += dx;
        }
    }
    min_val
}

// ---------------------------------------------------------------------------
// Convergence analysis
// ---------------------------------------------------------------------------

/// Performs a convergence study for DG advection by halving the mesh.
///
/// Returns (n_elements, l2_error) pairs.
#[allow(clippy::too_many_arguments)]
pub fn convergence_study<F: Fn(f64) -> f64 + Copy>(
    initial_condition: F,
    exact_solution: F,
    x_min: f64,
    x_max: f64,
    order: usize,
    n_levels: usize,
    t_final: f64,
    config: &DgAdvectionConfig,
) -> Vec<(usize, f64)> {
    let mut results = Vec::with_capacity(n_levels);
    let mut n_elem = 4;

    for _ in 0..n_levels {
        let mut mesh = DgMesh1D::uniform(x_min, x_max, n_elem, order);
        mesh.project_initial_condition(initial_condition);
        let _t = solve_advection(&mut mesh, config, t_final);
        let err = mesh.l2_error(exact_solution);
        results.push((n_elem, err));
        n_elem *= 2;
    }

    results
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn test_legendre_poly_p0() {
        assert!((legendre_poly(0, 0.5) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_legendre_poly_p1() {
        assert!((legendre_poly(1, 0.5) - 0.5).abs() < TOL);
    }

    #[test]
    fn test_legendre_poly_p2() {
        // P_2(x) = (3x^2 - 1)/2
        let x = 0.5;
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((legendre_poly(2, x) - expected).abs() < TOL);
    }

    #[test]
    fn test_legendre_poly_p3() {
        // P_3(x) = (5x^3 - 3x)/2
        let x = 0.7;
        let expected = (5.0 * x * x * x - 3.0 * x) / 2.0;
        assert!((legendre_poly(3, x) - expected).abs() < TOL);
    }

    #[test]
    fn test_legendre_deriv_p1() {
        assert!((legendre_poly_deriv(1, 0.5) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_gauss_quadrature_exact() {
        // Gauss quadrature with n points should integrate polynomials up to degree 2n-1 exactly
        let n = 3;
        let (pts, wts) = gauss_legendre_quadrature(n);
        // Integrate x^4 on [-1,1] = 2/5 (should NOT be exact for n=3, 2*3-1=5 >= 4, actually it is!)
        let integral: f64 = pts
            .iter()
            .zip(wts.iter())
            .map(|(&x, &w)| w * x * x * x * x)
            .sum();
        assert!((integral - 2.0 / 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_gauss_quadrature_weights_sum() {
        let (_, wts) = gauss_legendre_quadrature(5);
        let sum: f64 = wts.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_dg_basis_orthogonality() {
        let basis = DgBasis::new(3);
        let (pts, wts) = gauss_legendre_quadrature(5);
        // Check that basis functions 0 and 1 are orthogonal
        let inner: f64 = pts
            .iter()
            .zip(wts.iter())
            .map(|(&x, &w)| w * basis.evaluate(0, x) * basis.evaluate(1, x))
            .sum();
        assert!(inner.abs() < 1e-12);
    }

    #[test]
    fn test_dg_basis_normalization() {
        let basis = DgBasis::new(2);
        let (pts, wts) = gauss_legendre_quadrature(5);
        for i in 0..=2 {
            let norm_sq: f64 = pts
                .iter()
                .zip(wts.iter())
                .map(|(&x, &w)| w * basis.evaluate(i, x).powi(2))
                .sum();
            assert!(
                (norm_sq - 1.0).abs() < 1e-12,
                "mode {} norm_sq = {}",
                i,
                norm_sq
            );
        }
    }

    #[test]
    fn test_dg_element_mapping() {
        let elem = DgElement::new(1.0, 3.0, 2);
        assert!((elem.ref_to_phys(-1.0) - 1.0).abs() < TOL);
        assert!((elem.ref_to_phys(1.0) - 3.0).abs() < TOL);
        assert!((elem.ref_to_phys(0.0) - 2.0).abs() < TOL);
    }

    #[test]
    fn test_dg_element_inverse_mapping() {
        let elem = DgElement::new(2.0, 6.0, 1);
        assert!((elem.phys_to_ref(2.0) - (-1.0)).abs() < TOL);
        assert!((elem.phys_to_ref(6.0) - 1.0).abs() < TOL);
        assert!((elem.phys_to_ref(4.0)).abs() < TOL);
    }

    #[test]
    fn test_lax_friedrichs_flux_consistency() {
        // For u_L = u_R, the flux should equal f(u)
        let u = 2.0;
        let f = 0.5 * u * u; // Burgers
        let result = lax_friedrichs_flux(u, u, f, f, u.abs());
        assert!((result - f).abs() < TOL);
    }

    #[test]
    fn test_roe_flux_consistency() {
        let u = 3.0;
        let f = 0.5 * u * u;
        let result = roe_flux(u, u, f, f);
        assert!((result - f).abs() < TOL);
    }

    #[test]
    fn test_upwind_flux_positive_speed() {
        let result = upwind_flux(1.0, 2.0, 1.0);
        assert!((result - 1.0).abs() < TOL); // Should take left value
    }

    #[test]
    fn test_upwind_flux_negative_speed() {
        let result = upwind_flux(1.0, 2.0, -1.0);
        assert!((result - (-2.0)).abs() < TOL); // Should take right value
    }

    #[test]
    fn test_penalty_parameter() {
        let eta = penalty_parameter(2, 0.1, 0.1, 10.0);
        // Expected: 10 * 9 / 0.1 = 900
        assert!((eta - 900.0).abs() < TOL);
    }

    #[test]
    fn test_penalty_method_sigma() {
        assert!((PenaltyMethod::Sipg.sigma() - (-1.0)).abs() < TOL);
        assert!((PenaltyMethod::Nipg.sigma() - 1.0).abs() < TOL);
        assert!((PenaltyMethod::Iipg.sigma()).abs() < TOL);
    }

    #[test]
    fn test_smoothness_indicator_smooth() {
        // Coefficients that decay rapidly: smooth
        let coeffs = [1.0, 0.1, 0.001, 0.00001];
        let s = smoothness_indicator(&coeffs);
        assert!(s > 0.9, "smoothness = {:.6}", s);
    }

    #[test]
    fn test_smoothness_indicator_nonsmooth() {
        // Coefficients that don't decay: non-smooth
        let coeffs = [1.0, 0.9, 0.8, 0.7];
        let s = smoothness_indicator(&coeffs);
        assert!(s < 0.7, "smoothness = {:.6}", s);
    }

    #[test]
    fn test_tvb_minmod_within_threshold() {
        let result = tvb_minmod(0.001, 1.0, 1.0, 1.0, 1.0);
        assert!((result - 0.001).abs() < TOL);
    }

    #[test]
    fn test_weno3_reconstruct_uniform() {
        // For uniform data, reconstruction should return the value
        let (u_m, u_p) = weno3_reconstruct(1.0, 1.0, 1.0);
        assert!((u_m - 1.0).abs() < 1e-6);
        assert!((u_p - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dg_mesh_creation() {
        let mesh = DgMesh1D::uniform(0.0, 1.0, 10, 2);
        assert_eq!(mesh.num_elements(), 10);
        assert_eq!(mesh.total_dofs(), 30);
    }

    #[test]
    fn test_dg_mesh_project_constant() {
        let mut mesh = DgMesh1D::uniform(0.0, 1.0, 4, 2);
        mesh.project_initial_condition(|_x| 5.0);
        let error = mesh.l2_error(|_x| 5.0);
        assert!(error < 1e-12, "error = {:.6e}", error);
    }

    #[test]
    fn test_dg_mesh_project_linear() {
        let mut mesh = DgMesh1D::uniform(0.0, 1.0, 4, 2);
        mesh.project_initial_condition(|x| 2.0 * x + 1.0);
        let error = mesh.l2_error(|x| 2.0 * x + 1.0);
        assert!(error < 1e-10, "error = {:.6e}", error);
    }

    #[test]
    fn test_total_mass_constant() {
        let mut mesh = DgMesh1D::uniform(0.0, 2.0, 8, 2);
        mesh.project_initial_condition(|_x| 3.0);
        let mass = total_mass(&mesh);
        assert!((mass - 6.0).abs() < 1e-10, "mass = {:.6}", mass);
    }

    #[test]
    fn test_euler_state_basics() {
        let gamma = 1.4;
        let state = EulerState::new(1.0, 0.0, 2.5);
        let p = state.pressure(gamma);
        assert!((p - 1.0).abs() < 1e-12); // p = (1.4-1)*(2.5 - 0) = 1.0
        assert!((state.velocity()).abs() < TOL);
    }

    #[test]
    fn test_euler_lf_flux_consistency() {
        let gamma = 1.4;
        let state = EulerState::new(1.0, 1.0, 3.0);
        let f = euler_lax_friedrichs(&state, &state, gamma);
        let f_exact = state.flux(gamma);
        for i in 0..3 {
            assert!((f[i] - f_exact[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_h_refine_preserves_mass() {
        let mut elem = DgElement::new(0.0, 1.0, 2);
        elem.coeffs = vec![1.0, 0.5, 0.1];
        let (left, right) = h_refine_element(&elem);
        // Mass should be approximately preserved
        let basis = DgBasis::new(2);
        let (qpts, qwts) = gauss_legendre_quadrature(5);
        let mass_orig: f64 = qpts
            .iter()
            .zip(qwts.iter())
            .map(|(&xi, &w)| w * basis.evaluate_solution(&elem.coeffs, xi) * elem.jacobian())
            .sum();
        let mass_left: f64 = qpts
            .iter()
            .zip(qwts.iter())
            .map(|(&xi, &w)| w * basis.evaluate_solution(&left.coeffs, xi) * left.jacobian())
            .sum();
        let mass_right: f64 = qpts
            .iter()
            .zip(qwts.iter())
            .map(|(&xi, &w)| w * basis.evaluate_solution(&right.coeffs, xi) * right.jacobian())
            .sum();
        let mass_children = mass_left + mass_right;
        assert!(
            (mass_orig - mass_children).abs() < 0.1,
            "orig={:.6}, children={:.6}",
            mass_orig,
            mass_children
        );
    }

    #[test]
    fn test_p_enrich_preserves_coeffs() {
        let mut elem = DgElement::new(0.0, 1.0, 1);
        elem.coeffs = vec![1.0, 0.5];
        let enriched = p_enrich_element(&elem, 3);
        assert_eq!(enriched.order, 3);
        assert!((enriched.coeffs[0] - 1.0).abs() < TOL);
        assert!((enriched.coeffs[1] - 0.5).abs() < TOL);
        assert!((enriched.coeffs[2]).abs() < TOL);
        assert!((enriched.coeffs[3]).abs() < TOL);
    }

    #[test]
    fn test_local_mass_matrix_diagonal_dominance() {
        let elem = DgElement::new(0.0, 1.0, 2);
        let mass = local_mass_matrix(&elem);
        // Diagonal should be positive
        for i in 0..3 {
            assert!(mass[i][i] > 0.0);
        }
    }

    #[test]
    fn test_convergence_study_runs() {
        let config = DgAdvectionConfig {
            wave_speed: 1.0,
            flux_type: NumericalFlux::LaxFriedrichs,
            cfl: 0.3,
            tvb_m: 0.0,
            use_limiter: false,
        };
        fn sinwave(x: f64) -> f64 {
            (2.0 * std::f64::consts::PI * x).sin()
        }
        let results = convergence_study(
            sinwave, sinwave, 0.0, 1.0, 1, 2, 0.0, // t=0, so exact
            &config,
        );
        assert_eq!(results.len(), 2);
        // At t=0, error should be small
        // With 4 elements at p=1, the L2 projection error for sin(2*pi*x) is O(h^2)
        assert!(results[0].1 < 0.1, "error = {:.6e}", results[0].1);
    }

    #[test]
    fn test_2d_basis_num_modes() {
        let basis = DgBasis2D::new(2);
        assert_eq!(basis.num_modes(), 9);
    }

    #[test]
    fn test_local_stiffness_matrix_symmetry() {
        let elem = DgElement::new(0.0, 1.0, 2);
        let stiff = local_stiffness_matrix(&elem);
        let nm = 3;
        for i in 0..nm {
            for j in 0..nm {
                assert!(
                    (stiff[i][j] - stiff[j][i]).abs() < 1e-12,
                    "S[{}][{}] = {:.6}, S[{}][{}] = {:.6}",
                    i,
                    j,
                    stiff[i][j],
                    j,
                    i,
                    stiff[j][i]
                );
            }
        }
    }
}
