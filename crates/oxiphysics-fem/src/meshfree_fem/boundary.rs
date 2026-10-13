// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Essential (Dirichlet) boundary conditions for meshfree MLS approximations.
//!
//! The MLS shape functions of an EFG discretisation do not have the Kronecker
//! delta property, `φ_I(x_J) ≠ δ_IJ`.  The solution vector therefore holds
//! *parameters* `û_I` and the field is `u_h(x) = Σ_I φ_I(x) û_I`; fixing a
//! parameter to a value `g` does **not** make the field take that value at the
//! node.  All methods in this module consequently act on the field itself, at
//! boundary points `x_b` (see [`EssentialBc`]):
//!
//! * [`apply_penalty_bc`] adds the boundary-integral penalty
//!   `α ∫ (u_h − g)² dΓ`:  `K_IJ += α w φ_I(x_b) φ_J(x_b)`.
//! * [`apply_lagrange_bc`] appends the point-collocated constraint rows
//!   `Σ_J φ_J(x_b) û_J = g`.
//! * [`apply_nitsche_bc`] adds the symmetric Nitsche terms (consistency,
//!   adjoint consistency and stabilisation) built from the shape functions *and*
//!   their gradients at the boundary points; it is not a penalty method.
//!
//! [`essential_bcs_on_segment`] creates the boundary points of a straight edge.

use super::functions::{Point2, dist2, mls_shape_2d, plane_stress_d};
use super::types::{EfgMaterial, EssentialBc, KernelType, MlsShape, SupportDomain};
use crate::{Error, Result};

/// Number of nearest nodes averaged by [`local_node_spacing`].
const SPACING_NEIGHBOURS: usize = 4;

/// Validate the dense system arrays against the cloud and return `ndof`.
fn check_system(stiffness_len: usize, force_len: usize, domain: &SupportDomain) -> Result<usize> {
    let ndof = 2 * domain.node_count();
    if stiffness_len != ndof * ndof || force_len != ndof {
        return Err(Error::InvalidInput(format!(
            "a cloud of {} nodes needs a {ndof} x {ndof} stiffness matrix and a force vector of \
             length {ndof}, got {stiffness_len} and {force_len} entries",
            domain.node_count()
        )));
    }
    Ok(ndof)
}

/// Validate one boundary condition.
fn check_bc(bc: &EssentialBc, index: usize) -> Result<()> {
    if bc.dof > 1 {
        return Err(Error::InvalidInput(format!(
            "essential BC {index}: DOF {} is not 0 (x) or 1 (y)",
            bc.dof
        )));
    }
    if !bc.pos.iter().all(|c| c.is_finite()) || !bc.value.is_finite() {
        return Err(Error::InvalidInput(format!(
            "essential BC {index}: position and value must be finite"
        )));
    }
    if !bc.weight.is_finite() || bc.weight <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "essential BC {index}: the boundary weight must be positive, got {}",
            bc.weight
        )));
    }
    Ok(())
}

/// MLS shape functions at the position of a boundary condition.
fn shape_at(
    domain: &SupportDomain,
    ktype: KernelType,
    bc: &EssentialBc,
    index: usize,
) -> Result<MlsShape> {
    check_bc(bc, index)?;
    mls_shape_2d(bc.pos, domain, ktype).ok_or_else(|| {
        Error::InvalidInput(format!(
            "essential BC {index}: the MLS shape functions are undefined at ({}, {}) \
             (fewer than 3 nodes in range or a singular moment matrix)",
            bc.pos[0], bc.pos[1]
        ))
    })
}

/// Characteristic node spacing `h` around `pt`: the mean distance from `pt` to its
/// nearest nodes (at most four, a node coinciding with `pt` is skipped).
///
/// Used as the mesh size in Nitsche's stabilisation `β E / h`.  Returns `0.0`
/// when the cloud has no node distinct from `pt`.
pub fn local_node_spacing(domain: &SupportDomain, pt: Point2) -> f64 {
    let mut distances: Vec<f64> = domain.nodes.iter().map(|n| dist2(pt, n.pos)).collect();
    distances.sort_by(|a, b| a.total_cmp(b));
    let largest = distances.last().copied().unwrap_or(0.0);
    let cutoff = 1e-9 * largest;
    let nearest: Vec<f64> = distances
        .into_iter()
        .filter(|&d| d > cutoff)
        .take(SPACING_NEIGHBOURS)
        .collect();
    if nearest.is_empty() {
        0.0
    } else {
        nearest.iter().sum::<f64>() / nearest.len() as f64
    }
}

/// Apply essential BCs via the penalty method to a dense stiffness matrix and
/// force vector.
///
/// For every condition `(x_b, w, dof, g)` the penalty functional
/// `½ α w (u_h(x_b)[dof] − g)²` is added, i.e.
///
/// `K[(I, dof), (J, dof)] += α w φ_I(x_b) φ_J(x_b)`,
/// `f[(I, dof)]           += α w φ_I(x_b) g`
///
/// for all nodes `I`, `J` whose support contains `x_b`.  Since MLS shape
/// functions are not interpolatory, the update couples every node in range and
/// not just the diagonal of one node.  `penalty` has the dimension of a stiffness
/// per unit boundary length squared (`E / L` in 2-D); take it large compared with
/// the stiffness diagonal (`10⁴ … 10⁸ ×`).
///
/// # Errors
///
/// [`Error::InvalidInput`] when the arrays do not match the cloud, when
/// `penalty` is not positive and finite, when a condition is malformed
/// (DOF > 1, non-positive weight, non-finite data) or when the MLS shape
/// functions cannot be evaluated at a boundary point.  The system is left
/// unchanged on error.
pub fn apply_penalty_bc(
    stiffness: &mut [f64],
    force: &mut [f64],
    domain: &SupportDomain,
    ktype: KernelType,
    bcs: &[EssentialBc],
    penalty: f64,
) -> Result<()> {
    let ndof = check_system(stiffness.len(), force.len(), domain)?;
    if !penalty.is_finite() || penalty <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "the penalty parameter must be positive and finite, got {penalty}"
        )));
    }
    let shapes = bcs
        .iter()
        .enumerate()
        .map(|(i, bc)| shape_at(domain, ktype, bc, i))
        .collect::<Result<Vec<MlsShape>>>()?;
    for (bc, shape) in bcs.iter().zip(shapes.iter()) {
        let scale = penalty * bc.weight;
        for (a, &node_a) in shape.neighbours.iter().enumerate() {
            let row = 2 * node_a + bc.dof;
            force[row] += scale * shape.phi[a] * bc.value;
            for (b, &node_b) in shape.neighbours.iter().enumerate() {
                stiffness[row * ndof + 2 * node_b + bc.dof] += scale * shape.phi[a] * shape.phi[b];
            }
        }
    }
    Ok(())
}

/// Apply essential BCs via the Lagrange multiplier method.
///
/// Augments the system from `(K, f)` of size `n` to size `n + m`, `m` being the
/// number of conditions.  Constraint `c` is the point-collocated condition
/// `Σ_J φ_J(x_c) û_{J,dof} = g_c`, so its row holds the *shape function values*
/// `φ_J(x_c)` of all nodes in range (not a single unit entry).  The boundary
/// weights of the conditions are not used.  Returns the augmented system
/// `[[K, Bᵀ], [B, 0]]`, `[f; g]`.
///
/// The constraints must be independent: at most as many points per component as
/// the boundary trace of the approximation can resolve (e.g. the node positions
/// of a boundary), otherwise the augmented matrix is singular.
///
/// # Errors
///
/// [`Error::InvalidInput`] when the arrays do not match the cloud, when a
/// condition is malformed or when the MLS shape functions cannot be evaluated at
/// a boundary point.
pub fn apply_lagrange_bc(
    stiffness: &[f64],
    force: &[f64],
    domain: &SupportDomain,
    ktype: KernelType,
    bcs: &[EssentialBc],
) -> Result<(Vec<f64>, Vec<f64>)> {
    let ndof = check_system(stiffness.len(), force.len(), domain)?;
    let m = bcs.len();
    let n_aug = ndof + m;
    let mut k_aug = vec![0.0; n_aug * n_aug];
    let mut f_aug = vec![0.0; n_aug];
    for i in 0..ndof {
        k_aug[i * n_aug..i * n_aug + ndof].copy_from_slice(&stiffness[i * ndof..(i + 1) * ndof]);
    }
    f_aug[..ndof].copy_from_slice(force);

    for (c, bc) in bcs.iter().enumerate() {
        let shape = shape_at(domain, ktype, bc, c)?;
        for (a, &node) in shape.neighbours.iter().enumerate() {
            let column = 2 * node + bc.dof;
            k_aug[(ndof + c) * n_aug + column] = shape.phi[a];
            k_aug[column * n_aug + ndof + c] = shape.phi[a];
        }
        f_aug[ndof + c] = bc.value;
    }
    Ok((k_aug, f_aug))
}

/// Traction operator `T_I = n·D·B_I` of node `I` at a boundary point.
///
/// `T_I[c][j]` is the traction component `c` that a unit displacement of node
/// `I` in direction `j` produces on a surface with unit normal `normal`
/// (`σ·n` with `σ = D ε`, `ε = B_I e_j`).
fn traction_operator(
    dphi_dx: f64,
    dphi_dy: f64,
    normal: Point2,
    d_mat: &[[f64; 3]; 3],
) -> [[f64; 2]; 2] {
    // B_I (3 x 2): engineering strains (ε_xx, ε_yy, γ_xy) of a unit displacement
    let b = [[dphi_dx, 0.0], [0.0, dphi_dy], [dphi_dy, dphi_dx]];
    // n·σ with σ = (σ_xx, σ_yy, τ_xy): t_x = n_x σ_xx + n_y τ_xy, t_y = n_x τ_xy + n_y σ_yy
    let n_op = [[normal[0], 0.0, normal[1]], [0.0, normal[1], normal[0]]];
    let mut t = [[0.0; 2]; 2];
    for (c, t_row) in t.iter_mut().enumerate() {
        for (j, t_cj) in t_row.iter_mut().enumerate() {
            let mut value = 0.0;
            for (k, n_ck) in n_op[c].iter().enumerate() {
                for (l, b_lj) in b.iter().map(|row| row[j]).enumerate() {
                    value += n_ck * d_mat[k][l] * b_lj;
                }
            }
            *t_cj = value;
        }
    }
    t
}

/// Apply essential BCs via Nitsche's method.
///
/// Modifies the dense stiffness matrix and force vector in place with the
/// symmetric Nitsche terms for linear elasticity.  For a condition on component
/// `c` at the boundary point `x_b` (normal `n`, weight `w`, value `g`) the weak
/// form gains
///
/// ```text
/// −∫ (σ(u_h)·n)_c v_c  −  ∫ (σ(v_h)·n)_c (u_h,c − g)  +  γ ∫ (u_h,c − g) v_c
/// ```
///
/// i.e. consistency, adjoint-consistency and stabilisation terms, which in
/// matrix form read (thickness `t`, `T_I = n·D·B_I` the traction operator of
/// node `I`):
///
/// ```text
/// K_IJ += t w [ γ φ_I φ_J  −  φ_I T_J  −  T_Iᵀ φ_J ]      (restricted to component c)
/// f_I  += t w [ γ φ_I g    −  T_Iᵀ g ]
/// ```
///
/// The stabilisation is `γ = β D₁₁ / h` with `D₁₁` the plane-stress modulus and
/// `h` the local node spacing ([`local_node_spacing`]).  `beta` is
/// dimensionless; values from about 10 to 100 give a stable, accurate method (it
/// must be large enough for coercivity but, unlike a penalty, need not be huge).
/// The resulting contribution is symmetric and, unlike [`apply_penalty_bc`], it
/// contains the consistency terms built from the shape-function gradients, so a
/// field that satisfies the boundary condition and the equilibrium equations is
/// reproduced exactly.
///
/// # Errors
///
/// [`Error::InvalidInput`] when the arrays do not match the cloud, when `beta` is
/// not positive and finite, when a condition is malformed (a non-unit normal,
/// DOF > 1, non-positive weight, non-finite data), when the node spacing at a
/// boundary point cannot be determined, or when the MLS shape functions cannot be
/// evaluated at a boundary point.  The system is left unchanged on error.
pub fn apply_nitsche_bc(
    stiffness: &mut [f64],
    force: &mut [f64],
    domain: &SupportDomain,
    ktype: KernelType,
    mat: &EfgMaterial,
    bcs: &[EssentialBc],
    beta: f64,
) -> Result<()> {
    let ndof = check_system(stiffness.len(), force.len(), domain)?;
    if !beta.is_finite() || beta <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "the Nitsche parameter must be positive and finite, got {beta}"
        )));
    }
    let d_mat = plane_stress_d(mat.young, mat.poisson);
    let modulus = d_mat[0][0];

    // Evaluate everything first so that an error leaves the system untouched.
    struct PointData {
        shape: MlsShape,
        traction: Vec<[[f64; 2]; 2]>,
        gamma: f64,
    }
    let mut data = Vec::with_capacity(bcs.len());
    for (index, bc) in bcs.iter().enumerate() {
        let shape = shape_at(domain, ktype, bc, index)?;
        let normal_length = (bc.normal[0] * bc.normal[0] + bc.normal[1] * bc.normal[1]).sqrt();
        if !normal_length.is_finite() || (normal_length - 1.0).abs() > 1e-6 {
            return Err(Error::InvalidInput(format!(
                "essential BC {index}: Nitsche's method needs a unit outward normal, got ({}, {})",
                bc.normal[0], bc.normal[1]
            )));
        }
        let spacing = local_node_spacing(domain, bc.pos);
        if spacing <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "essential BC {index}: cannot determine the node spacing at ({}, {})",
                bc.pos[0], bc.pos[1]
            )));
        }
        let traction = shape
            .dphi_dx
            .iter()
            .zip(shape.dphi_dy.iter())
            .map(|(&dx, &dy)| traction_operator(dx, dy, bc.normal, &d_mat))
            .collect();
        data.push(PointData {
            shape,
            traction,
            gamma: beta * modulus / spacing,
        });
    }

    for (bc, point) in bcs.iter().zip(data.iter()) {
        let c = bc.dof;
        let scale = bc.weight * mat.thickness;
        let phi = &point.shape.phi;
        for (a, &node_a) in point.shape.neighbours.iter().enumerate() {
            let t_a = &point.traction[a];
            for i in 0..2 {
                // adjoint-consistency right-hand side: − T_Iᵀ g
                force[2 * node_a + i] -= scale * t_a[c][i] * bc.value;
            }
            // stabilisation right-hand side: γ φ_I g
            force[2 * node_a + c] += scale * point.gamma * phi[a] * bc.value;

            for (b, &node_b) in point.shape.neighbours.iter().enumerate() {
                let t_b = &point.traction[b];
                // consistency: −φ_I T_J, rows restricted to the constrained component
                for (j, t_bj) in t_b[c].iter().enumerate() {
                    stiffness[(2 * node_a + c) * ndof + 2 * node_b + j] -= scale * phi[a] * t_bj;
                }
                // adjoint consistency: −T_Iᵀ φ_J, columns restricted to component c
                for (i, t_ai) in t_a[c].iter().enumerate() {
                    stiffness[(2 * node_a + i) * ndof + 2 * node_b + c] -= scale * t_ai * phi[b];
                }
                // stabilisation: γ φ_I φ_J on the constrained component
                stiffness[(2 * node_a + c) * ndof + 2 * node_b + c] +=
                    scale * point.gamma * phi[a] * phi[b];
            }
        }
    }
    Ok(())
}

/// Gauss-Legendre nodes and weights on `[-1, 1]` for `n >= 1` points.
fn gauss_legendre_nodes(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    let nf = n as f64;
    for i in 0..n {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        let mut derivative = 1.0;
        for _ in 0..100 {
            let (mut p1, mut p2) = (1.0, 0.0);
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
            }
            derivative = nf * (z * p1 - p2) / (z * z - 1.0);
            let step = p1 / derivative;
            z -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        nodes.push(z);
        weights.push(2.0 / ((1.0 - z * z) * derivative * derivative));
    }
    (nodes, weights)
}

/// Essential boundary conditions along the straight edge `start → end`.
///
/// The edge is cut into `n_pieces` equal panels with `n_gauss` Gauss-Legendre
/// points each; every point becomes one [`EssentialBc`] per prescribed
/// component, carrying the quadrature weight (its share of the edge length)
/// and the **outward normal**.  The domain must lie to the *left* of
/// `start → end` (traverse the boundary counter-clockwise), so the outward
/// normal is the right-hand normal `(dy, −dx) / |d|`.
///
/// `displacement(x)` returns, for a boundary point `x`, the prescribed
/// `[u_x, u_y]` (`None` leaves a component free, e.g. for a roller).  Use at
/// least one panel per node spacing so that the boundary integrals are resolved.
/// An empty vector is returned for a degenerate edge or zero panels / points.
pub fn essential_bcs_on_segment(
    start: Point2,
    end: Point2,
    n_pieces: usize,
    n_gauss: usize,
    displacement: impl Fn(Point2) -> [Option<f64>; 2],
) -> Vec<EssentialBc> {
    let direction = [end[0] - start[0], end[1] - start[1]];
    let length = (direction[0] * direction[0] + direction[1] * direction[1]).sqrt();
    if !length.is_finite() || length <= 0.0 || n_pieces == 0 || n_gauss == 0 {
        return Vec::new();
    }
    let normal = [direction[1] / length, -direction[0] / length];
    let (xi, wi) = gauss_legendre_nodes(n_gauss);
    let panel = 1.0 / n_pieces as f64;
    let mut bcs = Vec::with_capacity(n_pieces * n_gauss * 2);
    for piece in 0..n_pieces {
        let t0 = piece as f64 * panel;
        for (x, w) in xi.iter().zip(wi.iter()) {
            let t = t0 + 0.5 * panel * (x + 1.0);
            let pos = [start[0] + t * direction[0], start[1] + t * direction[1]];
            let weight = 0.5 * panel * w * length;
            let prescribed = displacement(pos);
            for (dof, value) in prescribed.iter().enumerate() {
                if let Some(value) = value {
                    bcs.push(EssentialBc::new(pos, normal, weight, dof, *value));
                }
            }
        }
    }
    bcs
}
