// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Local, sparse isogeometric assembly of plane-stress elasticity.
//!
//! A bivariate NURBS surface of degrees `(p, q)` has exactly `(p + 1)(q + 1)`
//! non-zero basis functions on every knot-span rectangle (an IGA *element*).
//! Each element therefore carries only that local connectivity, builds a local
//! stiffness matrix of size `2(p + 1)(q + 1)` and scatters it into a sparse
//! [`CsrMatrix`](crate::sparse::CsrMatrix).  The memory and work of an assembly
//! are linear in the number of elements, and the resulting system is solved with
//! Jacobi-preconditioned conjugate gradients
//! ([`pcg_solve`](crate::solvers::pcg_solve)).
//!
//! The local rational basis functions and their parametric derivatives are
//! evaluated with the non-recursive Algorithms A2.2 / A2.3 of Piegl & Tiller
//! ([`basis_derivatives`]) and the physical gradients follow from the inverse
//! transpose of the Jacobian, `∇ₓR = J⁻ᵀ ∇_ξR`.

use super::{BSplineBasis, IgaElement, IgaMaterial, NurbsSurface, gauss_legendre, gauss_map};
use crate::isogeometric::functions::{basis_derivatives, find_knot_span};
use crate::solvers::{jacobi_preconditioner, pcg_solve};
use crate::sparse::CsrMatrix;
use crate::{Error, Result};

/// Relative threshold below which a Jacobian determinant is treated as singular.
const SINGULAR_JACOBIAN: f64 = 1e-14;

/// Rational basis functions of one IGA element and their parametric derivatives.
///
/// All three vectors have one entry per active control point, in the order of
/// [`IgaElement::basis_indices`]: the `u` index runs fastest.
#[derive(Debug, Clone)]
pub struct ElementBasis {
    /// Rational basis functions `R_a` (they sum to one).
    pub values: Vec<f64>,
    /// Derivatives `∂R_a/∂u`.
    pub du: Vec<f64>,
    /// Derivatives `∂R_a/∂v`.
    pub dv: Vec<f64>,
}

/// Evaluate the non-zero univariate B-splines of the span that contains `span_point`.
///
/// Returns `(span, N, dN/dt)` at parameter `t` (`t` may lie anywhere in the span;
/// the polynomial piece of the span is evaluated).
fn span_basis(
    basis: &BSplineBasis,
    num_basis: usize,
    span_point: f64,
    t: f64,
) -> Result<(usize, Vec<f64>, Vec<f64>)> {
    let degree = basis.degree;
    if num_basis < degree + 1 || basis.knots.len() != num_basis + degree + 1 {
        return Err(Error::InvalidInput(format!(
            "inconsistent B-spline basis: {num_basis} functions, degree {degree}, {} knots",
            basis.knots.len()
        )));
    }
    let (lo, hi) = (basis.knots[degree], basis.knots[num_basis]);
    let tol = 1e-12 * (hi - lo).abs().max(1.0);
    if !(span_point >= lo - tol && span_point <= hi + tol) {
        return Err(Error::InvalidInput(format!(
            "parameter {span_point} lies outside the knot domain [{lo}, {hi}]"
        )));
    }
    let span = find_knot_span(num_basis - 1, degree, span_point, &basis.knots);
    let ders = basis_derivatives(span, degree, t, &basis.knots, 1);
    let values = ders[0].clone();
    let derivatives = ders
        .get(1)
        .cloned()
        .unwrap_or_else(|| vec![0.0; degree + 1]);
    Ok((span, values, derivatives))
}

impl IgaElement {
    /// Rational basis functions and parametric derivatives of this element at
    /// the parametric point `(u, v)`.
    ///
    /// Only the `(p + 1)(q + 1)` basis functions that are non-zero on the
    /// element's knot span are evaluated (Piegl & Tiller A2.2 / A2.3); the
    /// vectors follow [`IgaElement::basis_indices`].  The span is taken from the
    /// element's parametric domain, so `(u, v)` may lie on the element boundary.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when `basis_indices` is not the connectivity of
    /// the knot span `param_domain` of `surface` (wrong length or wrong control
    /// points), when the surface bases are inconsistent, or when the NURBS
    /// weight function vanishes at `(u, v)`.
    pub fn local_basis(&self, surface: &NurbsSurface, u: f64, v: f64) -> Result<ElementBasis> {
        let p = surface.basis_u.degree;
        let q = surface.basis_v.degree;
        let n_loc = (p + 1) * (q + 1);
        if self.basis_indices.len() != n_loc {
            return Err(Error::InvalidInput(format!(
                "element has {} active control points, a ({p}, {q}) surface needs {n_loc}",
                self.basis_indices.len()
            )));
        }
        if surface.control_points.len() != surface.num_u * surface.num_v {
            return Err(Error::InvalidInput(format!(
                "surface stores {} control points for a {} x {} net",
                surface.control_points.len(),
                surface.num_u,
                surface.num_v
            )));
        }
        let u_mid = 0.5 * (self.param_domain[0] + self.param_domain[1]);
        let v_mid = 0.5 * (self.param_domain[2] + self.param_domain[3]);
        let (span_u, nu, dnu) = span_basis(&surface.basis_u, surface.num_u, u_mid, u)?;
        let (span_v, nv, dnv) = span_basis(&surface.basis_v, surface.num_v, v_mid, v)?;

        // The connectivity must be the one of this knot span.
        for b in 0..=q {
            for a in 0..=p {
                let expected = (span_v - q + b) * surface.num_u + (span_u - p + a);
                if self.basis_indices[b * (p + 1) + a] != expected {
                    return Err(Error::InvalidInput(
                        "element connectivity does not match its parametric domain".to_string(),
                    ));
                }
            }
        }

        let mut nw = vec![0.0; n_loc];
        let mut dnw_du = vec![0.0; n_loc];
        let mut dnw_dv = vec![0.0; n_loc];
        let mut w_sum = 0.0;
        let mut dw_du = 0.0;
        let mut dw_dv = 0.0;
        for b in 0..=q {
            for a in 0..=p {
                let loc = b * (p + 1) + a;
                let weight = surface.control_points[self.basis_indices[loc]].weight;
                nw[loc] = nu[a] * nv[b] * weight;
                dnw_du[loc] = dnu[a] * nv[b] * weight;
                dnw_dv[loc] = nu[a] * dnv[b] * weight;
                w_sum += nw[loc];
                dw_du += dnw_du[loc];
                dw_dv += dnw_dv[loc];
            }
        }
        if w_sum.abs() < 1e-14 || !w_sum.is_finite() {
            return Err(Error::InvalidInput(format!(
                "the NURBS weight function vanishes at ({u}, {v})"
            )));
        }

        // Rational derivatives: dR/du = (dNw/du * W - Nw * dW/du) / W^2
        let inv_w = 1.0 / w_sum;
        let inv_w2 = inv_w * inv_w;
        let values: Vec<f64> = nw.iter().map(|n| n * inv_w).collect();
        let du: Vec<f64> = dnw_du
            .iter()
            .zip(nw.iter())
            .map(|(d, n)| (d * w_sum - n * dw_du) * inv_w2)
            .collect();
        let dv: Vec<f64> = dnw_dv
            .iter()
            .zip(nw.iter())
            .map(|(d, n)| (d * w_sum - n * dw_dv) * inv_w2)
            .collect();
        Ok(ElementBasis { values, du, dv })
    }

    /// Jacobian `∂(x, y)/∂(u, v)` from an already evaluated local basis.
    ///
    /// Returns `[[∂x/∂u, ∂x/∂v], [∂y/∂u, ∂y/∂v]]` built from the active control
    /// points only (`O((p + 1)(q + 1))`).
    pub fn local_jacobian(&self, surface: &NurbsSurface, basis: &ElementBasis) -> [[f64; 2]; 2] {
        let mut jac = [[0.0; 2]; 2];
        for (loc, &global) in self.basis_indices.iter().enumerate() {
            let pos = surface.control_points[global].pos;
            jac[0][0] += basis.du[loc] * pos[0];
            jac[0][1] += basis.dv[loc] * pos[0];
            jac[1][0] += basis.du[loc] * pos[1];
            jac[1][1] += basis.dv[loc] * pos[1];
        }
        jac
    }

    /// NURBS shape functions of this element at the parametric point `(u, v)`.
    ///
    /// One value per active control point, in the order of
    /// [`IgaElement::basis_indices`]; the values sum to one.
    ///
    /// # Errors
    ///
    /// See [`local_basis`](Self::local_basis).
    pub fn shape_functions(&self, surface: &NurbsSurface, u: f64, v: f64) -> Result<Vec<f64>> {
        Ok(self.local_basis(surface, u, v)?.values)
    }

    /// Strain-displacement matrix and Jacobian determinant at `(u, v)`.
    ///
    /// The matrix is `None` when the Jacobian is singular.
    fn b_matrix_checked(
        &self,
        surface: &NurbsSurface,
        u: f64,
        v: f64,
    ) -> Result<(Option<Vec<f64>>, f64)> {
        let basis = self.local_basis(surface, u, v)?;
        let jac = self.local_jacobian(surface, &basis);
        let det_j = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        if is_singular(&jac, det_j) {
            return Ok((None, det_j));
        }
        let n_loc = self.basis_indices.len();
        let ndof = 2 * n_loc;
        let inv_j = [
            [jac[1][1] / det_j, -jac[0][1] / det_j],
            [-jac[1][0] / det_j, jac[0][0] / det_j],
        ];
        let mut b = vec![0.0; 3 * ndof];
        for a in 0..n_loc {
            // ∇ₓR = J⁻ᵀ ∇_ξR
            let dn_dx = inv_j[0][0] * basis.du[a] + inv_j[1][0] * basis.dv[a];
            let dn_dy = inv_j[0][1] * basis.du[a] + inv_j[1][1] * basis.dv[a];
            // B = [dN/dx  0     ]
            //     [0      dN/dy ]
            //     [dN/dy  dN/dx ]
            b[2 * a] = dn_dx;
            b[ndof + 2 * a + 1] = dn_dy;
            b[2 * ndof + 2 * a] = dn_dy;
            b[2 * ndof + 2 * a + 1] = dn_dx;
        }
        Ok((Some(b), det_j))
    }

    /// Compute the B-matrix (strain-displacement) at `(u, v)`.
    ///
    /// Returns `(B, det J)` where `B` is a flat row-major `3 × 2n` matrix
    /// (`n` = number of active control points, columns `2a` / `2a + 1` for the x /
    /// y displacement of control point `a`) for 2-D plane problems, and `det J`
    /// is the Jacobian determinant of the map `(u, v) → (x, y)`.  The physical
    /// gradients are `∇ₓR = J⁻ᵀ ∇_ξR`, i.e.
    ///
    /// `∂R/∂x = (J⁻¹)₀₀ ∂R/∂u + (J⁻¹)₁₀ ∂R/∂v`,
    /// `∂R/∂y = (J⁻¹)₀₁ ∂R/∂u + (J⁻¹)₁₁ ∂R/∂v`.
    ///
    /// For a singular Jacobian the returned `B` is zero.
    ///
    /// # Errors
    ///
    /// See [`local_basis`](Self::local_basis).
    pub fn b_matrix(&self, surface: &NurbsSurface, u: f64, v: f64) -> Result<(Vec<f64>, f64)> {
        let (b, det_j) = self.b_matrix_checked(surface, u, v)?;
        let b = b.unwrap_or_else(|| vec![0.0; 3 * 2 * self.basis_indices.len()]);
        Ok((b, det_j))
    }

    /// Compute the element stiffness matrix with Gauss quadrature.
    ///
    /// Returns the *local* matrix of size `2n × 2n` (flattened row-major),
    /// `n = (p + 1)(q + 1)` being the number of active control points; row /
    /// column `2a + c` belongs to component `c` of local control point `a`
    /// (global control point `basis_indices[a]`).  Plane stress is assumed.
    /// Quadrature points with a singular Jacobian do not contribute.
    ///
    /// # Errors
    ///
    /// See [`local_basis`](Self::local_basis).
    pub fn stiffness_matrix(
        &self,
        surface: &NurbsSurface,
        material: &IgaMaterial,
    ) -> Result<Vec<f64>> {
        let ndof = 2 * self.basis_indices.len();
        let mut k = vec![0.0; ndof * ndof];
        let d_mat = material.plane_stress_matrix();

        let (gp_ref, gw_ref) = gauss_legendre(self.num_gauss);
        let (gp_u, gw_u) = gauss_map(&gp_ref, &gw_ref, self.param_domain[0], self.param_domain[1]);
        let (gp_v, gw_v) = gauss_map(&gp_ref, &gw_ref, self.param_domain[2], self.param_domain[3]);

        let mut db = vec![0.0; 3 * ndof];
        for (qu, wu) in gp_u.iter().zip(gw_u.iter()) {
            for (qv, wv) in gp_v.iter().zip(gw_v.iter()) {
                let (b, det_j) = self.b_matrix_checked(surface, *qu, *qv)?;
                let Some(b) = b else {
                    continue;
                };
                let weight = wu * wv * det_j.abs() * material.thickness;

                // DB = D * B  (3 x ndof)
                for r in 0..3 {
                    for c in 0..ndof {
                        db[r * ndof + c] = (0..3).map(|s| d_mat[r][s] * b[s * ndof + c]).sum();
                    }
                }
                // K += weight * B^T * (D B); exploit symmetry
                for row in 0..ndof {
                    for col in row..ndof {
                        let val: f64 = (0..3)
                            .map(|r| b[r * ndof + row] * db[r * ndof + col])
                            .sum::<f64>()
                            * weight;
                        k[row * ndof + col] += val;
                        if row != col {
                            k[col * ndof + row] += val;
                        }
                    }
                }
            }
        }
        Ok(k)
    }
}

/// `true` when a 2×2 Jacobian is (numerically) singular or not finite.
fn is_singular(jac: &[[f64; 2]; 2], det_j: f64) -> bool {
    let scale = (jac[0][0] * jac[1][1]).abs() + (jac[0][1] * jac[1][0]).abs();
    !det_j.is_finite() || det_j.abs() <= SINGULAR_JACOBIAN * scale
}

// ---------------------------------------------------------------------------
// IGA Assembler
// ---------------------------------------------------------------------------

/// Assembled IGA system: sparse stiffness matrix and load vector.
///
/// The unknowns are ordered `2 * control_point + component`, the control points
/// in the row-major order of the surface (`u` fastest).
#[derive(Debug, Clone)]
pub struct IgaSystem {
    /// Global stiffness matrix (sparse CSR, `ndof × ndof`).
    pub stiffness: CsrMatrix,
    /// Global load vector (n).
    pub load: Vec<f64>,
    /// Total DOFs.
    pub ndof: usize,
}

impl IgaSystem {
    /// Solve `K u = f` with Jacobi-preconditioned conjugate gradients.
    ///
    /// `tolerance` is the relative residual `‖f - K u‖ / ‖f‖`.  Apply the
    /// boundary conditions first: the free-floating system is singular.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidInput`] when the matrix and the load vector disagree
    ///   with `ndof`.
    /// * [`Error::NotConverged`] when CG does not reach `tolerance` within
    ///   `max_iter` iterations (typically an unconstrained rigid-body mode).
    /// * [`Error::SingularSystem`] when the solution contains a non-finite value.
    pub fn solve(&self, max_iter: usize, tolerance: f64) -> Result<Vec<f64>> {
        if self.stiffness.nrows != self.ndof
            || self.stiffness.ncols != self.ndof
            || self.load.len() != self.ndof
        {
            return Err(Error::InvalidInput(format!(
                "system of {} DOFs holds a {} x {} matrix and a load of length {}",
                self.ndof,
                self.stiffness.nrows,
                self.stiffness.ncols,
                self.load.len()
            )));
        }
        let precond = jacobi_preconditioner(&self.stiffness);
        let (solution, stats) =
            pcg_solve(&self.stiffness, &self.load, &precond, max_iter, tolerance);
        if !stats.converged {
            let load_norm = self.load.iter().map(|x| x * x).sum::<f64>().sqrt();
            return Err(Error::NotConverged {
                solver: "PCG".to_string(),
                iterations: stats.iterations,
                residual: stats.residual / load_norm.max(f64::MIN_POSITIVE),
            });
        }
        if solution.iter().any(|x| !x.is_finite()) {
            return Err(Error::SingularSystem);
        }
        Ok(solution)
    }

    /// The stiffness matrix as a dense row-major `ndof × ndof` array.
    ///
    /// Intended for small systems (tests, debugging): it allocates `ndof²`
    /// entries.
    pub fn dense_stiffness(&self) -> Vec<f64> {
        let mut dense = vec![0.0; self.ndof * self.ndof];
        for row in 0..self.stiffness.nrows {
            for idx in self.stiffness.row_ptr[row]..self.stiffness.row_ptr[row + 1] {
                dense[row * self.ndof + self.stiffness.col_indices[idx]] =
                    self.stiffness.values[idx];
            }
        }
        dense
    }
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

    /// Check that `surface` is consistent and matches `num_control_points`.
    fn validate_surface(&self, surface: &NurbsSurface) -> Result<()> {
        let net = surface.num_u * surface.num_v;
        if surface.control_points.len() != net {
            return Err(Error::InvalidInput(format!(
                "surface stores {} control points for a {} x {} net",
                surface.control_points.len(),
                surface.num_u,
                surface.num_v
            )));
        }
        if self.num_control_points != net {
            return Err(Error::InvalidInput(format!(
                "assembler is configured for {} control points but the surface has {} x {} = {net}",
                self.num_control_points, surface.num_u, surface.num_v
            )));
        }
        if self.num_gauss == 0 {
            return Err(Error::InvalidInput(
                "at least one Gauss point per direction is required".to_string(),
            ));
        }
        Ok(())
    }

    /// Assemble the global stiffness matrix and an (empty) load vector.
    ///
    /// Every knot span of the surface is an element with its own
    /// `(p + 1)(q + 1)` active control points; the local stiffness matrices are
    /// scattered into a sparse CSR matrix, so the cost is linear in the number
    /// of elements.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when `num_control_points` differs from
    /// `num_u * num_v` of the surface, when the surface is inconsistent or when
    /// an element cannot be evaluated (see [`IgaElement::local_basis`]).
    pub fn assemble(&self, surface: &NurbsSurface) -> Result<IgaSystem> {
        self.validate_surface(surface)?;
        let ndof = 2 * self.num_control_points;
        let elements = self.create_elements(surface);

        let n_loc = (surface.basis_u.degree + 1) * (surface.basis_v.degree + 1);
        let mut triplets = Vec::with_capacity(elements.len() * (2 * n_loc) * (2 * n_loc));
        for elem in &elements {
            let k_elem = elem.stiffness_matrix(surface, &self.material)?;
            let n_elem_dof = 2 * elem.basis_indices.len();
            for (i, &gi) in elem.basis_indices.iter().enumerate() {
                for (j, &gj) in elem.basis_indices.iter().enumerate() {
                    for ci in 0..2 {
                        for cj in 0..2 {
                            let value = k_elem[(2 * i + ci) * n_elem_dof + 2 * j + cj];
                            if value != 0.0 {
                                triplets.push((2 * gi + ci, 2 * gj + cj, value));
                            }
                        }
                    }
                }
            }
        }
        Ok(IgaSystem {
            stiffness: CsrMatrix::from_triplets(ndof, ndof, &triplets),
            load: vec![0.0; ndof],
            ndof,
        })
    }

    /// Create IGA elements from non-zero knot spans.
    ///
    /// Each element knows the `(p + 1)(q + 1)` control points whose basis
    /// functions are non-zero on its knot-span rectangle (global index
    /// `j * num_u + i`, `u` fastest).
    pub fn create_elements(&self, surface: &NurbsSurface) -> Vec<IgaElement> {
        let mut elements = Vec::new();
        let p = surface.basis_u.degree;
        let q = surface.basis_v.degree;
        let knots_u = &surface.basis_u.knots;
        let knots_v = &surface.basis_v.knots;
        if surface.num_u < p + 1
            || surface.num_v < q + 1
            || knots_u.len() != surface.num_u + p + 1
            || knots_v.len() != surface.num_v + q + 1
        {
            return elements;
        }

        // Spans of the valid parameter domain: p <= i <= num_u - 1.
        for i in p..surface.num_u {
            if (knots_u[i + 1] - knots_u[i]).abs() < 1e-14 {
                continue;
            }
            for j in q..surface.num_v {
                if (knots_v[j + 1] - knots_v[j]).abs() < 1e-14 {
                    continue;
                }
                let mut indices = Vec::with_capacity((p + 1) * (q + 1));
                for b in 0..=q {
                    for a in 0..=p {
                        indices.push((j - q + b) * surface.num_u + (i - p + a));
                    }
                }
                let domain = [knots_u[i], knots_u[i + 1], knots_v[j], knots_v[j + 1]];
                elements.push(IgaElement::new(indices, domain, self.num_gauss));
            }
        }
        elements
    }

    /// Apply a uniform body force (e.g., gravity) to the load vector.
    ///
    /// Adds `∫ R_a f t dA` to the two load entries of every control point; the
    /// resultant of the added load is `f · area · thickness`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] when the surface does not match the assembler, when
    /// `system.load` has the wrong length, or when an element cannot be
    /// evaluated.
    pub fn apply_body_force(
        &self,
        surface: &NurbsSurface,
        system: &mut IgaSystem,
        force: [f64; 2],
    ) -> Result<()> {
        self.validate_surface(surface)?;
        if system.load.len() != 2 * self.num_control_points {
            return Err(Error::InvalidInput(format!(
                "load vector has {} entries, expected {}",
                system.load.len(),
                2 * self.num_control_points
            )));
        }
        let elements = self.create_elements(surface);
        let (gp_ref, gw_ref) = gauss_legendre(self.num_gauss);

        for elem in &elements {
            let (gp_u, gw_u) =
                gauss_map(&gp_ref, &gw_ref, elem.param_domain[0], elem.param_domain[1]);
            let (gp_v, gw_v) =
                gauss_map(&gp_ref, &gw_ref, elem.param_domain[2], elem.param_domain[3]);

            for (qu, wu) in gp_u.iter().zip(gw_u.iter()) {
                for (qv, wv) in gp_v.iter().zip(gw_v.iter()) {
                    let basis = elem.local_basis(surface, *qu, *qv)?;
                    let jac = elem.local_jacobian(surface, &basis);
                    let det_j = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
                    let weight = wu * wv * det_j.abs() * self.material.thickness;

                    for (loc, &global) in elem.basis_indices.iter().enumerate() {
                        system.load[2 * global] += basis.values[loc] * force[0] * weight;
                        system.load[2 * global + 1] += basis.values[loc] * force[1] * weight;
                    }
                }
            }
        }
        Ok(())
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

    /// Check that `system` is consistent and that every DOF index addresses it.
    fn check_range(&self, system: &IgaSystem) -> Result<()> {
        if system.stiffness.nrows != system.ndof || system.load.len() != system.ndof {
            return Err(Error::InvalidInput(format!(
                "system of {} DOFs holds a matrix with {} rows and a load of length {}",
                system.ndof,
                system.stiffness.nrows,
                system.load.len()
            )));
        }
        match self.dof_indices.iter().find(|&&dof| dof >= system.ndof) {
            Some(&dof) => Err(Error::DofOutOfRange {
                dof,
                n: system.ndof,
            }),
            None => Ok(()),
        }
    }

    /// Apply Dirichlet BC to the system using the penalty method.
    ///
    /// Adds `penalty` to the diagonal entry of every constrained DOF and
    /// `penalty * value` to its load entry.  Nothing is modified for a
    /// Neumann condition.
    ///
    /// # Errors
    ///
    /// [`Error::DofOutOfRange`] when a DOF index is not smaller than
    /// `system.ndof`; the system is left untouched in that case.
    pub fn apply_dirichlet(&self, system: &mut IgaSystem, penalty: f64) -> Result<()> {
        if self.bc_type != IgaBcType::Dirichlet {
            return Ok(());
        }
        self.check_range(system)?;
        for (dof, val) in self.dof_indices.iter().zip(self.values.iter()) {
            let d = *dof;
            system.stiffness.add_to(d, d, penalty);
            system.load[d] += penalty * val;
        }
        Ok(())
    }

    /// Apply Neumann BC (add to load vector).
    ///
    /// # Errors
    ///
    /// [`Error::DofOutOfRange`] when a DOF index is not smaller than
    /// `system.ndof`; the system is left untouched in that case.
    pub fn apply_neumann(&self, system: &mut IgaSystem) -> Result<()> {
        if self.bc_type != IgaBcType::Neumann {
            return Ok(());
        }
        self.check_range(system)?;
        for (dof, val) in self.dof_indices.iter().zip(self.values.iter()) {
            system.load[*dof] += val;
        }
        Ok(())
    }

    /// Apply this BC to the system (dispatches based on type).
    ///
    /// # Errors
    ///
    /// See [`apply_dirichlet`](Self::apply_dirichlet) and
    /// [`apply_neumann`](Self::apply_neumann).
    pub fn apply(&self, system: &mut IgaSystem, penalty: f64) -> Result<()> {
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
