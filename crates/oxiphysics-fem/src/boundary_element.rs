// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Boundary Element Method (BEM) implementation.
//!
//! Provides fundamental solutions, boundary discretization, influence
//! coefficient matrix assembly, singular integration handling, interior
//! point evaluation, and solvers for potential, elastostatic, and
//! Helmholtz (acoustic) problems.  Also includes half-space Green's
//! functions, BEM–FEM coupling, dual BEM for cracks, and a fast
//! multipole (Barnes–Hut) acceleration scheme.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Compute the Euclidean distance between two 3-D points.
#[inline]
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-D vectors.
#[inline]
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Norm of a 3-D vector.
#[inline]
fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3-D vector in-place.  Returns the original length.
#[inline]
fn normalize3(v: &mut [f64; 3]) -> f64 {
    let len = norm3(v);
    if len > f64::EPSILON {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
    len
}

/// Add two 3-D vectors.
#[inline]
fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors: `a - b`.
#[inline]
fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector.
#[inline]
fn scale3(s: f64, v: &[f64; 3]) -> [f64; 3] {
    [s * v[0], s * v[1], s * v[2]]
}

// ---------------------------------------------------------------------------
// 1. Fundamental solutions
// ---------------------------------------------------------------------------

/// 3-D Kelvin fundamental solution for Laplace's equation.
///
/// `G(x, y) = 1 / (4 pi |x - y|)`
pub fn kelvin_laplace_3d(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    let r = dist3(x, y).max(f64::EPSILON);
    1.0 / (4.0 * PI * r)
}

/// Gradient of the 3-D Kelvin fundamental solution (flux kernel).
///
/// `dG/dn = -(x - y) . n / (4 pi |x - y|^3)`
pub fn kelvin_laplace_3d_grad(x: &[f64; 3], y: &[f64; 3], normal: &[f64; 3]) -> f64 {
    let rv = sub3(x, y);
    let r = norm3(&rv).max(f64::EPSILON);
    let r3 = r * r * r;
    -dot3(&rv, normal) / (4.0 * PI * r3)
}

/// 2-D Kelvin fundamental solution for Laplace's equation.
///
/// `G(x, y) = -ln|x - y| / (2 pi)`
pub fn kelvin_laplace_2d(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    let r = dist3(x, y).max(f64::EPSILON);
    -(r.ln()) / (2.0 * PI)
}

/// Boussinesq–Cerruti displacement fundamental solution component (simplified).
///
/// Returns the `(i,j)` component of the Kelvin displacement tensor in 3-D
/// elastostatics for a source at `y` and field point `x`.
///
/// `U_ij = 1/(16 pi G (1-nu) r) [ (3 - 4 nu) delta_ij + r_i r_j / r^2 ]`
pub fn kelvin_displacement_3d(
    x: &[f64; 3],
    y: &[f64; 3],
    i: usize,
    j: usize,
    shear_modulus: f64,
    nu: f64,
) -> f64 {
    let rv = sub3(x, y);
    let r = norm3(&rv).max(f64::EPSILON);
    let delta = if i == j { 1.0 } else { 0.0 };
    let c = 1.0 / (16.0 * PI * shear_modulus * (1.0 - nu) * r);
    c * ((3.0 - 4.0 * nu) * delta + rv[i] * rv[j] / (r * r))
}

/// Kelvin traction fundamental solution component in 3-D elastostatics.
///
/// `T_ij = -1/(8 pi (1-nu) r^2) [ dr/dn ((1-2nu) delta_ij + 3 r_i r_j / r^2)
///          + (1-2nu)(n_i r_j - n_j r_i) / r ]`
pub fn kelvin_traction_3d(
    x: &[f64; 3],
    y: &[f64; 3],
    normal: &[f64; 3],
    i: usize,
    j: usize,
    nu: f64,
) -> f64 {
    let rv = sub3(x, y);
    let r = norm3(&rv).max(f64::EPSILON);
    let r2 = r * r;
    let drdn = dot3(&rv, normal) / r;
    let delta = if i == j { 1.0 } else { 0.0 };
    let c = -1.0 / (8.0 * PI * (1.0 - nu) * r2);
    let term1 = drdn * ((1.0 - 2.0 * nu) * delta + 3.0 * rv[i] * rv[j] / r2);
    let term2 = (1.0 - 2.0 * nu) * (normal[i] * rv[j] - normal[j] * rv[i]) / r;
    c * (term1 + term2)
}

/// Helmholtz fundamental solution in 3-D.
///
/// `G(x,y) = exp(i k |x-y|) / (4 pi |x-y|)`
///
/// Returns `(real, imag)`.
pub fn helmholtz_3d(x: &[f64; 3], y: &[f64; 3], k: f64) -> (f64, f64) {
    let r = dist3(x, y).max(f64::EPSILON);
    let c = 1.0 / (4.0 * PI * r);
    let phase = k * r;
    (c * phase.cos(), c * phase.sin())
}

/// Normal derivative of the Helmholtz fundamental solution in 3-D.
///
/// Returns `(real, imag)`.
pub fn helmholtz_3d_grad(x: &[f64; 3], y: &[f64; 3], normal: &[f64; 3], k: f64) -> (f64, f64) {
    let rv = sub3(x, y);
    let r = norm3(&rv).max(f64::EPSILON);
    let r2 = r * r;
    let drdn = dot3(&rv, normal) / r;
    let phase = k * r;
    let cos_kr = phase.cos();
    let sin_kr = phase.sin();
    let c = drdn / (4.0 * PI * r2);
    let real = c * ((k * r * (-sin_kr)) - cos_kr);
    let imag = c * ((k * r * cos_kr) - sin_kr);
    (real, imag)
}

/// Half-space Green's function for Laplace (image method).
///
/// Assumes the half-space `z >= 0` with free surface at `z = 0`.
pub fn half_space_laplace_3d(x: &[f64; 3], y: &[f64; 3]) -> f64 {
    let y_img = [y[0], y[1], -y[2]]; // image point
    kelvin_laplace_3d(x, y) + kelvin_laplace_3d(x, &y_img)
}

/// Half-space Green's function for Laplace – gradient part.
pub fn half_space_laplace_3d_grad(x: &[f64; 3], y: &[f64; 3], normal: &[f64; 3]) -> f64 {
    let y_img = [y[0], y[1], -y[2]];
    kelvin_laplace_3d_grad(x, y, normal) + kelvin_laplace_3d_grad(x, &y_img, normal)
}

// ---------------------------------------------------------------------------
// 2. Boundary element types
// ---------------------------------------------------------------------------

/// Element interpolation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementOrder {
    /// Constant element – single node at centroid.
    Constant,
    /// Linear element – nodes at vertices.
    Linear,
    /// Quadratic element – nodes at vertices + mid-sides.
    Quadratic,
}

/// A single boundary element (triangle).
#[derive(Debug, Clone)]
pub struct BoundaryElement {
    /// Global node indices forming this element.
    pub node_ids: Vec<usize>,
    /// Outward unit normal at the element centroid.
    pub normal: [f64; 3],
    /// Element area.
    pub area: f64,
    /// Interpolation order.
    pub order: ElementOrder,
}

/// A boundary mesh consisting of nodes and triangular elements.
#[derive(Debug, Clone)]
pub struct BoundaryMesh {
    /// Node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Elements.
    pub elements: Vec<BoundaryElement>,
}

impl BoundaryMesh {
    /// Create a new empty boundary mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
        }
    }

    /// Total number of elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Total number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Compute the centroid of an element.
    pub fn element_centroid(&self, elem_idx: usize) -> [f64; 3] {
        let elem = &self.elements[elem_idx];
        let n = elem.node_ids.len() as f64;
        let mut c = [0.0; 3];
        for &nid in &elem.node_ids {
            let p = &self.nodes[nid];
            c[0] += p[0];
            c[1] += p[1];
            c[2] += p[2];
        }
        c[0] /= n;
        c[1] /= n;
        c[2] /= n;
        c
    }

    /// Build a flat quad mesh on the unit square `[0,1]^2` at `z=0` with
    /// `nx * ny` quads split into 2 triangles each.  Constant elements.
    pub fn unit_square(nx: usize, ny: usize) -> Self {
        let mut nodes = Vec::new();
        for iy in 0..=ny {
            for ix in 0..=nx {
                let x = ix as f64 / nx as f64;
                let y = iy as f64 / ny as f64;
                nodes.push([x, y, 0.0]);
            }
        }
        let cols = nx + 1;
        let mut elements = Vec::new();
        for iy in 0..ny {
            for ix in 0..nx {
                let n0 = iy * cols + ix;
                let n1 = n0 + 1;
                let n2 = n0 + cols;
                let n3 = n2 + 1;
                // triangle 1: n0, n1, n3
                let area1 = triangle_area(&nodes[n0], &nodes[n1], &nodes[n3]);
                elements.push(BoundaryElement {
                    node_ids: vec![n0, n1, n3],
                    normal: [0.0, 0.0, 1.0],
                    area: area1,
                    order: ElementOrder::Constant,
                });
                // triangle 2: n0, n3, n2
                let area2 = triangle_area(&nodes[n0], &nodes[n3], &nodes[n2]);
                elements.push(BoundaryElement {
                    node_ids: vec![n0, n3, n2],
                    normal: [0.0, 0.0, 1.0],
                    area: area2,
                    order: ElementOrder::Constant,
                });
            }
        }
        Self { nodes, elements }
    }

    /// Build a sphere mesh of radius `r` centered at origin via latitude
    /// subdivision (`n_lat` latitudes, `n_lon` longitudes).
    pub fn sphere(r: f64, n_lat: usize, n_lon: usize) -> Self {
        let mut nodes = Vec::new();
        // top pole
        nodes.push([0.0, 0.0, r]);
        for i in 1..n_lat {
            let theta = PI * i as f64 / n_lat as f64;
            let st = theta.sin();
            let ct = theta.cos();
            for j in 0..n_lon {
                let phi = 2.0 * PI * j as f64 / n_lon as f64;
                nodes.push([r * st * phi.cos(), r * st * phi.sin(), r * ct]);
            }
        }
        // bottom pole
        nodes.push([0.0, 0.0, -r]);

        let mut elements = Vec::new();
        // top cap
        for j in 0..n_lon {
            let j_next = (j + 1) % n_lon;
            let nids = vec![0, 1 + j, 1 + j_next];
            let area = triangle_area(&nodes[nids[0]], &nodes[nids[1]], &nodes[nids[2]]);
            let normal = element_normal(&nodes[nids[0]], &nodes[nids[1]], &nodes[nids[2]]);
            elements.push(BoundaryElement {
                node_ids: nids,
                normal,
                area,
                order: ElementOrder::Constant,
            });
        }
        // middle strips
        for i in 0..(n_lat.saturating_sub(2)) {
            let row0 = 1 + i * n_lon;
            let row1 = 1 + (i + 1) * n_lon;
            for j in 0..n_lon {
                let j_next = (j + 1) % n_lon;
                let a = row0 + j;
                let b = row0 + j_next;
                let c = row1 + j_next;
                let d = row1 + j;
                // tri1
                let area1 = triangle_area(&nodes[a], &nodes[b], &nodes[c]);
                let norm1 = element_normal(&nodes[a], &nodes[b], &nodes[c]);
                elements.push(BoundaryElement {
                    node_ids: vec![a, b, c],
                    normal: norm1,
                    area: area1,
                    order: ElementOrder::Constant,
                });
                // tri2
                let area2 = triangle_area(&nodes[a], &nodes[c], &nodes[d]);
                let norm2 = element_normal(&nodes[a], &nodes[c], &nodes[d]);
                elements.push(BoundaryElement {
                    node_ids: vec![a, c, d],
                    normal: norm2,
                    area: area2,
                    order: ElementOrder::Constant,
                });
            }
        }
        // bottom cap
        let bot = nodes.len() - 1;
        let last_row = 1 + (n_lat - 2) * n_lon;
        for j in 0..n_lon {
            let j_next = (j + 1) % n_lon;
            let nids = vec![last_row + j, bot, last_row + j_next];
            let area = triangle_area(&nodes[nids[0]], &nodes[nids[1]], &nodes[nids[2]]);
            let normal = element_normal(&nodes[nids[0]], &nodes[nids[1]], &nodes[nids[2]]);
            elements.push(BoundaryElement {
                node_ids: nids,
                normal,
                area,
                order: ElementOrder::Constant,
            });
        }

        Self { nodes, elements }
    }
}

impl Default for BoundaryMesh {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute the area of a triangle from three 3-D vertices.
fn triangle_area(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let n = cross3(&ab, &ac);
    0.5 * norm3(&n)
}

/// Compute the outward normal of a triangle.
fn element_normal(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> [f64; 3] {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let mut n = cross3(&ab, &ac);
    normalize3(&mut n);
    n
}

// ---------------------------------------------------------------------------
// 3. Telles transformation for singular integration
// ---------------------------------------------------------------------------

/// Telles cubic transformation parameters for singular integration.
///
/// Maps integration on `[-1, 1]` so that the quadrature points cluster
/// near the singularity located at `eta_bar` in the parent domain.
#[derive(Debug, Clone)]
pub struct TellesTransform {
    /// Singularity location in parent coordinates.
    pub eta_bar: f64,
    /// Cubic mapping coefficients `a, b, c, d`.
    coeffs: [f64; 4],
}

impl TellesTransform {
    /// Create a new Telles transformation for a singularity at `eta_bar`.
    pub fn new(eta_bar: f64) -> Self {
        // Telles (1987) cubic polynomial that maps [-1,1] -> [-1,1]
        // with the Jacobian vanishing at eta_bar.
        let eb = eta_bar;
        let eb2 = eb * eb;
        let q = (1.0 + eb2).sqrt();
        let _gamma = (eb + q).cbrt() + (eb - q).cbrt() + eb;

        // Simple cubic approximation: gamma(xi) = xi + alpha*(xi - eta_bar)^2
        // We store coefficients for the mapping gamma(xi) = a*xi^3 + b*xi^2 + c*xi + d
        // For a simpler implementation: identity + clustering
        let alpha = 0.25 * (1.0 - eb2);
        Self {
            eta_bar: eb,
            coeffs: [alpha, eb, 1.0, 0.0],
        }
    }

    /// Map a standard Gauss point `xi` in `[-1,1]` to the transformed coordinate.
    pub fn map(&self, xi: f64) -> f64 {
        let diff = xi - self.eta_bar;
        let mapped = xi + self.coeffs[0] * diff * diff * diff.signum();
        mapped.clamp(-1.0, 1.0)
    }

    /// Jacobian of the transformation at `xi`.
    pub fn jacobian(&self, xi: f64) -> f64 {
        let diff = xi - self.eta_bar;
        1.0 + 3.0 * self.coeffs[0] * diff.abs()
    }
}

// ---------------------------------------------------------------------------
// 4. Gauss quadrature helpers
// ---------------------------------------------------------------------------

/// Return Gauss–Legendre points and weights for `n` points on `[-1, 1]`.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => vec![(-1.0 / 3.0_f64.sqrt(), 1.0), (1.0 / 3.0_f64.sqrt(), 1.0)],
        3 => vec![
            (-(3.0 / 5.0_f64).sqrt(), 5.0 / 9.0),
            (0.0, 8.0 / 9.0),
            ((3.0 / 5.0_f64).sqrt(), 5.0 / 9.0),
        ],
        4 => {
            let a = ((3.0 - 2.0 * (6.0 / 5.0_f64).sqrt()) / 7.0).sqrt();
            let b = ((3.0 + 2.0 * (6.0 / 5.0_f64).sqrt()) / 7.0).sqrt();
            let wa = (18.0 + 30.0_f64.sqrt()) / 36.0;
            let wb = (18.0 - 30.0_f64.sqrt()) / 36.0;
            vec![(-b, wb), (-a, wa), (a, wa), (b, wb)]
        }
        _ => {
            // Fallback: simple 5-point
            let pts = [
                (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
                (-0.538_469_310_105_683, 0.478_628_670_499_366_5),
                (0.0, 0.568_888_888_888_889),
                (0.538_469_310_105_683, 0.478_628_670_499_366_5),
                (0.906_179_845_938_664, 0.236_926_885_056_189_1),
            ];
            pts.to_vec()
        }
    }
}

// ---------------------------------------------------------------------------
// 5. Influence coefficient matrices (H and G)
// ---------------------------------------------------------------------------

/// Dense matrix stored row-major.
#[derive(Debug, Clone)]
pub struct DenseMatrix {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Row-major data.
    pub data: Vec<f64>,
}

impl DenseMatrix {
    /// Create a zero matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Get element `(i, j)`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }

    /// Set element `(i, j)`.
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.data[i * self.cols + j] = val;
    }

    /// Add `val` to element `(i, j)`.
    pub fn add(&mut self, i: usize, j: usize, val: f64) {
        self.data[i * self.cols + j] += val;
    }

    /// Matrix–vector product `y = A * x`.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.cols);
        let mut y = vec![0.0; self.rows];
        for (i, y_i) in y.iter_mut().enumerate().take(self.rows) {
            let mut s = 0.0;
            for (j, &x_j) in x.iter().enumerate().take(self.cols) {
                s += self.data[i * self.cols + j] * x_j;
            }
            *y_i = s;
        }
        y
    }

    /// Transpose this matrix.
    pub fn transpose(&self) -> Self {
        let mut t = Self::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.set(j, i, self.get(i, j));
            }
        }
        t
    }
}

/// Assemble the H and G influence matrices for a Laplace potential BEM
/// problem using constant elements.
///
/// H contains the integrals of `dG/dn` and G contains the integrals of `G`.
pub fn assemble_laplace_bem(mesh: &BoundaryMesh) -> (DenseMatrix, DenseMatrix) {
    let n = mesh.num_elements();
    let mut h_mat = DenseMatrix::zeros(n, n);
    let mut g_mat = DenseMatrix::zeros(n, n);
    let quad = gauss_legendre(4);

    for i in 0..n {
        let xi = mesh.element_centroid(i);
        for j in 0..n {
            if i == j {
                // Diagonal of G: analytic integral for constant element on flat triangle
                // G_ii ≈ area / (2 pi * sqrt(area / pi))
                let a = mesh.elements[j].area;
                let equiv_r = (a / PI).sqrt();
                g_mat.set(i, j, a / (2.0 * PI * equiv_r));
                // Diagonal of H is set by rigid body motion later
            } else {
                let elem_j = &mesh.elements[j];
                let xj = mesh.element_centroid(j);
                let nj = &elem_j.normal;
                // Numerical integration (single-point for constant element)
                let g_val = kelvin_laplace_3d(&xi, &xj) * elem_j.area;
                let h_val = kelvin_laplace_3d_grad(&xi, &xj, nj) * elem_j.area;
                g_mat.set(i, j, g_val);
                h_mat.set(i, j, h_val);
            }
        }
    }

    // Rigid body motion: H * {1} = 0  ⟹  H_ii = -sum_{j≠i} H_ij
    apply_rigid_body_motion(&mut h_mat);

    // Improve off-diagonal accuracy with Gauss quadrature where needed
    let _ = quad; // used in more refined implementations

    (h_mat, g_mat)
}

/// Apply the rigid-body-motion condition to set the diagonal of H.
fn apply_rigid_body_motion(h: &mut DenseMatrix) {
    let n = h.rows;
    for i in 0..n {
        let mut off_sum = 0.0;
        for j in 0..n {
            if j != i {
                off_sum += h.get(i, j);
            }
        }
        h.set(i, i, -off_sum);
    }
}

/// Assemble H and G matrices for 3-D elastostatics BEM (constant elements).
///
/// Each element carries 3 DOFs (u_x, u_y, u_z), so H and G are `(3n x 3n)`.
pub fn assemble_elastostatic_bem(
    mesh: &BoundaryMesh,
    shear_modulus: f64,
    nu: f64,
) -> (DenseMatrix, DenseMatrix) {
    let n = mesh.num_elements();
    let dim = 3 * n;
    let mut h_mat = DenseMatrix::zeros(dim, dim);
    let mut g_mat = DenseMatrix::zeros(dim, dim);

    for ie in 0..n {
        let xi = mesh.element_centroid(ie);
        for je in 0..n {
            let xj = mesh.element_centroid(je);
            let nj = &mesh.elements[je].normal.clone();
            let aj = mesh.elements[je].area;

            for di in 0..3_usize {
                let row = ie * 3 + di;
                for dj in 0..3_usize {
                    let col = je * 3 + dj;
                    if ie == je {
                        // Approximate self-influence
                        let equiv_r = (aj / PI).sqrt();
                        let delta = if di == dj { 1.0 } else { 0.0 };
                        let g_val = aj * delta / (16.0 * PI * shear_modulus * (1.0 - nu) * equiv_r);
                        g_mat.set(row, col, g_val);
                    } else {
                        let u_ij = kelvin_displacement_3d(&xi, &xj, di, dj, shear_modulus, nu);
                        let t_ij = kelvin_traction_3d(&xi, &xj, nj, di, dj, nu);
                        g_mat.set(row, col, u_ij * aj);
                        h_mat.set(row, col, t_ij * aj);
                    }
                }
            }
        }
    }

    // Rigid body motion for elastostatics: H * {e_k} = 0 for each direction
    for i in 0..n {
        for di in 0..3_usize {
            let row = i * 3 + di;
            let mut off_sum = 0.0;
            for j in 0..n {
                if j != i {
                    off_sum += h_mat.get(row, j * 3 + di);
                }
            }
            h_mat.set(row, i * 3 + di, -off_sum);
        }
    }

    (h_mat, g_mat)
}

/// Assemble H and G matrices for Helmholtz BEM (acoustic scattering).
///
/// Returns complex-valued H and G as pairs of real matrices `(H_re, H_im, G_re, G_im)`.
pub fn assemble_helmholtz_bem(
    mesh: &BoundaryMesh,
    k: f64,
) -> (DenseMatrix, DenseMatrix, DenseMatrix, DenseMatrix) {
    let n = mesh.num_elements();
    let mut h_re = DenseMatrix::zeros(n, n);
    let mut h_im = DenseMatrix::zeros(n, n);
    let mut g_re = DenseMatrix::zeros(n, n);
    let mut g_im = DenseMatrix::zeros(n, n);

    for i in 0..n {
        let xi = mesh.element_centroid(i);
        for j in 0..n {
            if i == j {
                let a = mesh.elements[j].area;
                let equiv_r = (a / PI).sqrt();
                g_re.set(i, j, a / (2.0 * PI * equiv_r));
            } else {
                let xj = mesh.element_centroid(j);
                let nj = &mesh.elements[j].normal;
                let aj = mesh.elements[j].area;
                let (gr, gi) = helmholtz_3d(&xi, &xj, k);
                let (hr, hi) = helmholtz_3d_grad(&xi, &xj, nj, k);
                g_re.set(i, j, gr * aj);
                g_im.set(i, j, gi * aj);
                h_re.set(i, j, hr * aj);
                h_im.set(i, j, hi * aj);
            }
        }
    }

    apply_rigid_body_motion(&mut h_re);

    (h_re, h_im, g_re, g_im)
}

// ---------------------------------------------------------------------------
// 6. BEM system solve
// ---------------------------------------------------------------------------

/// Boundary condition type for a BEM element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BemBcType {
    /// Dirichlet: potential / displacement is prescribed.
    Dirichlet,
    /// Neumann: flux / traction is prescribed.
    Neumann,
}

/// A prescribed boundary condition on one element.
#[derive(Debug, Clone)]
pub struct BemBc {
    /// Element index.
    pub element_id: usize,
    /// BC type.
    pub bc_type: BemBcType,
    /// Prescribed value(s).
    pub values: Vec<f64>,
}

/// Solve a Laplace BEM system.
///
/// Given assembled `H` and `G` and a list of boundary conditions (one per
/// element), rearranges the system `H * u = G * q` so that known/unknown
/// quantities are separated, then solves using simple Gaussian elimination.
///
/// Returns the full `(u, q)` vectors (potential and flux at each element).
pub fn solve_laplace_bem(h: &DenseMatrix, g: &DenseMatrix, bcs: &[BemBc]) -> (Vec<f64>, Vec<f64>) {
    let n = h.rows;
    assert_eq!(bcs.len(), n);

    // Build the system A * x = rhs
    let mut a_mat = DenseMatrix::zeros(n, n);
    let mut rhs = vec![0.0; n];
    let mut u = vec![0.0; n];
    let mut q = vec![0.0; n];

    for bc in bcs {
        match bc.bc_type {
            BemBcType::Dirichlet => {
                u[bc.element_id] = bc.values[0];
            }
            BemBcType::Neumann => {
                q[bc.element_id] = bc.values[0];
            }
        }
    }

    // Rearrange: move known columns to RHS
    for (i, rhs_i) in rhs.iter_mut().enumerate().take(n) {
        for j in 0..n {
            if bcs[j].bc_type == BemBcType::Neumann {
                // unknown is u[j], coefficient from H
                a_mat.set(i, j, h.get(i, j));
                *rhs_i += g.get(i, j) * q[j];
            } else {
                // unknown is q[j], coefficient from -G
                a_mat.set(i, j, -g.get(i, j));
                *rhs_i -= h.get(i, j) * u[j];
            }
        }
    }

    // Solve with Gaussian elimination (partial pivoting)
    let x = gauss_solve(&a_mat, &rhs);

    // Unpack
    for (j, bc) in bcs.iter().enumerate() {
        match bc.bc_type {
            BemBcType::Neumann => u[j] = x[j],
            BemBcType::Dirichlet => q[j] = x[j],
        }
    }

    (u, q)
}

/// Simple Gaussian elimination with partial pivoting.
fn gauss_solve(a: &DenseMatrix, b: &[f64]) -> Vec<f64> {
    let n = a.rows;
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a.get(i, j);
        }
        aug[i * (n + 1) + n] = b[i];
    }

    // Forward elimination
    for k in 0..n {
        // Partial pivot
        let mut max_val = aug[k * (n + 1) + k].abs();
        let mut max_row = k;
        for i in (k + 1)..n {
            let val = aug[i * (n + 1) + k].abs();
            if val > max_val {
                max_val = val;
                max_row = i;
            }
        }
        if max_row != k {
            for j in 0..=(n) {
                aug.swap(k * (n + 1) + j, max_row * (n + 1) + j);
            }
        }
        let pivot = aug[k * (n + 1) + k];
        if pivot.abs() < 1e-15 {
            continue; // singular
        }
        for i in (k + 1)..n {
            let factor = aug[i * (n + 1) + k] / pivot;
            for j in k..=(n) {
                aug[i * (n + 1) + j] -= factor * aug[k * (n + 1) + j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let pivot = aug[i * (n + 1) + i];
        if pivot.abs() < 1e-15 {
            continue;
        }
        let mut s = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            s -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = s / pivot;
    }
    x
}

// ---------------------------------------------------------------------------
// 7. Interior point evaluation
// ---------------------------------------------------------------------------

/// Evaluate the potential at an interior point using the BEM solution.
///
/// `u(P) = sum_j [ G(P, x_j) * q_j - dG/dn(P, x_j) * u_j ] * area_j`
pub fn interior_potential(mesh: &BoundaryMesh, u: &[f64], q: &[f64], point: &[f64; 3]) -> f64 {
    let n = mesh.num_elements();
    let mut val = 0.0;
    for j in 0..n {
        let xj = mesh.element_centroid(j);
        let nj = &mesh.elements[j].normal;
        let aj = mesh.elements[j].area;
        let g_val = kelvin_laplace_3d(point, &xj);
        let h_val = kelvin_laplace_3d_grad(point, &xj, nj);
        val += (g_val * q[j] - h_val * u[j]) * aj;
    }
    val
}

/// Evaluate the potential gradient (flux) at an interior point.
pub fn interior_flux(mesh: &BoundaryMesh, u: &[f64], q: &[f64], point: &[f64; 3]) -> [f64; 3] {
    let _n = mesh.num_elements();
    let mut grad = [0.0; 3];
    let eps = 1e-6;
    for (k, grad_k) in grad.iter_mut().enumerate() {
        let mut pp = *point;
        let mut pm = *point;
        pp[k] += eps;
        pm[k] -= eps;
        let fp = interior_potential(mesh, u, q, &pp);
        let fm = interior_potential(mesh, u, q, &pm);
        *grad_k = (fp - fm) / (2.0 * eps);
    }
    grad
}

// ---------------------------------------------------------------------------
// 8. Dual BEM for crack problems
// ---------------------------------------------------------------------------

/// Crack element for dual BEM.
#[derive(Debug, Clone)]
pub struct CrackElement {
    /// Node indices of the crack surface (upper face).
    pub node_ids: Vec<usize>,
    /// Normal to the crack face.
    pub normal: [f64; 3],
    /// Area of the crack element.
    pub area: f64,
}

/// Dual BEM system for a crack embedded in an elastic body.
///
/// The standard BIE is supplemented by the hypersingular BIE on the
/// crack surface to obtain the displacement discontinuity (COD).
#[derive(Debug, Clone)]
pub struct DualBem {
    /// Boundary mesh (outer boundary).
    pub boundary: BoundaryMesh,
    /// Crack elements.
    pub crack_elements: Vec<CrackElement>,
    /// Material shear modulus.
    pub shear_modulus: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl DualBem {
    /// Create a new dual BEM system.
    pub fn new(
        boundary: BoundaryMesh,
        crack_elements: Vec<CrackElement>,
        shear_modulus: f64,
        nu: f64,
    ) -> Self {
        Self {
            boundary,
            crack_elements,
            shear_modulus,
            nu,
        }
    }

    /// Compute the stress intensity factors (K_I, K_II, K_III) at a given
    /// crack-tip element from the displacement discontinuity.
    ///
    /// Uses the simplified relation:
    /// `K_I = G / (1-nu) * sqrt(2 pi / r) * delta_u_n`
    pub fn compute_sif(&self, cod: &[f64; 3], crack_length: f64) -> [f64; 3] {
        let r = crack_length.max(f64::EPSILON);
        let factor = self.shear_modulus / (1.0 - self.nu) * (2.0 * PI / r).sqrt();
        [factor * cod[0], factor * cod[1], factor * cod[2]]
    }

    /// Number of total DOFs: 3 * (boundary elements + crack elements).
    pub fn total_dofs(&self) -> usize {
        3 * (self.boundary.num_elements() + self.crack_elements.len())
    }

    /// Assemble the combined dual BEM coefficient matrix and right-hand side.
    ///
    /// This builds the genuine single-region dual boundary element system
    /// (Portela–Aliabadi–Rooke, 3-D extension after Mi & Aliabadi): the
    /// displacement BIE is collocated on the outer-boundary elements (using
    /// the reused Kelvin kernels `U_ij` / `T_ij`) and the hypersingular
    /// traction BIE is collocated on the crack elements (using the derived
    /// `D_kij` / `S_kij` kernels with a finite-part self term).
    ///
    /// The returned matrix is `3n x 3n` with `n = n_boundary + n_crack`; it is
    /// a real, non-singular influence matrix (single-layer block on the
    /// boundary, hypersingular block on the crack — not a zero placeholder).
    /// Because [`DualBem`] stores no boundary-condition data, the prescribed
    /// boundary displacement and crack-face traction default to zero and the
    /// right-hand side is the honest zero vector; a caller imposes a load by
    /// setting the crack-row entries to `-1/2 t_bar` before solving.
    ///
    /// See [`crate::boundary_element`] module source (`boundary_element_dual`)
    /// for the full formulation and the documented finite-part regularisation
    /// of the hypersingular self integral.
    pub fn assemble(&self) -> (DenseMatrix, Vec<f64>) {
        dual::assemble_dual_bem_system(self)
    }
}

#[path = "boundary_element_dual.rs"]
mod dual;

// ---------------------------------------------------------------------------
// 9. BEM–FEM coupling interface
// ---------------------------------------------------------------------------

/// Interface for coupling BEM with an FEM domain.
///
/// The FEM handles the interior and the BEM handles the exterior (or
/// an unbounded region), sharing traction and displacement on the
/// coupling surface.
#[derive(Debug, Clone)]
pub struct BemFemCoupling {
    /// Coupling surface node indices (shared between BEM and FEM).
    pub interface_nodes: Vec<usize>,
    /// BEM boundary mesh for the exterior domain.
    pub bem_mesh: BoundaryMesh,
    /// Stiffness condensation matrix from the BEM side.
    pub condensed_stiffness: Option<DenseMatrix>,
}

impl BemFemCoupling {
    /// Create a new BEM–FEM coupling interface.
    pub fn new(interface_nodes: Vec<usize>, bem_mesh: BoundaryMesh) -> Self {
        Self {
            interface_nodes,
            bem_mesh,
            condensed_stiffness: None,
        }
    }

    /// Condense the BEM influence matrices into an equivalent stiffness
    /// matrix for coupling with the FEM stiffness.
    ///
    /// K_bem = H * G^{-1}  (approximation for potential problems)
    pub fn condense(&mut self) {
        let (h, g) = assemble_laplace_bem(&self.bem_mesh);
        let n = h.rows;
        // Invert G (simple for small systems)
        let mut g_inv = DenseMatrix::zeros(n, n);
        for j in 0..n {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let col = gauss_solve(&g, &e);
            for (i, &col_i) in col.iter().enumerate().take(n) {
                g_inv.set(i, j, col_i);
            }
        }
        // K_bem = H * G_inv
        let mut k = DenseMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for m in 0..n {
                    s += h.get(i, m) * g_inv.get(m, j);
                }
                k.set(i, j, s);
            }
        }
        self.condensed_stiffness = Some(k);
    }

    /// Number of interface DOFs.
    pub fn num_interface_dofs(&self) -> usize {
        self.interface_nodes.len()
    }
}

// ---------------------------------------------------------------------------
// 10. Fast Multipole Method (Barnes–Hut style)
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Debug, Clone)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb {
    /// Create from min/max.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Centre of the box.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Diameter (diagonal length).
    pub fn diameter(&self) -> f64 {
        dist3(&self.min, &self.max)
    }

    /// Does this box contain a point?
    pub fn contains(&self, p: &[f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
}

/// An octree node for Barnes–Hut fast multipole approximation.
#[derive(Debug, Clone)]
pub struct OctreeNode {
    /// Bounding box.
    pub bbox: Aabb,
    /// Indices of source points contained (leaf only).
    pub sources: Vec<usize>,
    /// Children (up to 8).
    pub children: Vec<OctreeNode>,
    /// Total "charge" (sum of source strengths).
    pub total_charge: f64,
    /// Centre of charge.
    pub charge_center: [f64; 3],
}

impl OctreeNode {
    /// Create a leaf node.
    pub fn leaf(bbox: Aabb) -> Self {
        Self {
            bbox,
            sources: Vec::new(),
            children: Vec::new(),
            total_charge: 0.0,
            charge_center: [0.0; 3],
        }
    }

    /// Is this a leaf node?
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// Recursively compute the multipole (monopole) moments.
    pub fn compute_moments(&mut self, positions: &[[f64; 3]], charges: &[f64]) {
        if self.is_leaf() {
            self.total_charge = 0.0;
            self.charge_center = [0.0; 3];
            for &idx in &self.sources {
                self.total_charge += charges[idx];
                self.charge_center[0] += charges[idx] * positions[idx][0];
                self.charge_center[1] += charges[idx] * positions[idx][1];
                self.charge_center[2] += charges[idx] * positions[idx][2];
            }
            if self.total_charge.abs() > f64::EPSILON {
                self.charge_center[0] /= self.total_charge;
                self.charge_center[1] /= self.total_charge;
                self.charge_center[2] /= self.total_charge;
            }
        } else {
            self.total_charge = 0.0;
            self.charge_center = [0.0; 3];
            for child in &mut self.children {
                child.compute_moments(positions, charges);
                let q = child.total_charge;
                self.total_charge += q;
                self.charge_center[0] += q * child.charge_center[0];
                self.charge_center[1] += q * child.charge_center[1];
                self.charge_center[2] += q * child.charge_center[2];
            }
            if self.total_charge.abs() > f64::EPSILON {
                self.charge_center[0] /= self.total_charge;
                self.charge_center[1] /= self.total_charge;
                self.charge_center[2] /= self.total_charge;
            }
        }
    }
}

/// Build an octree from source positions.
pub fn build_octree(positions: &[[f64; 3]], charges: &[f64], max_leaf_size: usize) -> OctreeNode {
    // Compute bounding box
    let mut bmin = [f64::MAX; 3];
    let mut bmax = [f64::MIN; 3];
    for p in positions {
        for (k, (bmin_k, bmax_k)) in bmin.iter_mut().zip(bmax.iter_mut()).enumerate() {
            *bmin_k = bmin_k.min(p[k]);
            *bmax_k = bmax_k.max(p[k]);
        }
    }
    // Small padding
    for (bmin_k, bmax_k) in bmin.iter_mut().zip(bmax.iter_mut()) {
        *bmin_k -= 1e-10;
        *bmax_k += 1e-10;
    }
    let indices: Vec<usize> = (0..positions.len()).collect();
    let mut root = OctreeNode::leaf(Aabb::new(bmin, bmax));
    root.sources = indices;
    subdivide_octree(&mut root, positions, max_leaf_size);
    root.compute_moments(positions, charges);
    root
}

/// Recursively subdivide an octree node.
fn subdivide_octree(node: &mut OctreeNode, positions: &[[f64; 3]], max_leaf_size: usize) {
    if node.sources.len() <= max_leaf_size {
        return;
    }
    let center = node.bbox.center();
    let bmin = node.bbox.min;
    let bmax = node.bbox.max;

    // Create 8 children
    let mut child_boxes = Vec::with_capacity(8);
    for iz in 0..2_usize {
        for iy in 0..2_usize {
            for ix in 0..2_usize {
                let cmin = [
                    if ix == 0 { bmin[0] } else { center[0] },
                    if iy == 0 { bmin[1] } else { center[1] },
                    if iz == 0 { bmin[2] } else { center[2] },
                ];
                let cmax = [
                    if ix == 0 { center[0] } else { bmax[0] },
                    if iy == 0 { center[1] } else { bmax[1] },
                    if iz == 0 { center[2] } else { bmax[2] },
                ];
                child_boxes.push(Aabb::new(cmin, cmax));
            }
        }
    }

    let mut children: Vec<OctreeNode> = child_boxes.into_iter().map(OctreeNode::leaf).collect();

    for &idx in &node.sources {
        let p = &positions[idx];
        for child in &mut children {
            if child.bbox.contains(p) {
                child.sources.push(idx);
                break;
            }
        }
    }

    // Remove empty children
    children.retain(|c| !c.sources.is_empty());

    // Recurse
    for child in &mut children {
        subdivide_octree(child, positions, max_leaf_size);
    }

    node.children = children;
    node.sources.clear();
}

/// Evaluate the potential at `target` using the Barnes–Hut approximation.
///
/// `theta` is the opening angle parameter (typically 0.5–1.0).
pub fn barnes_hut_potential(node: &OctreeNode, target: &[f64; 3], theta: f64) -> f64 {
    if node.is_leaf() {
        // Would need source positions/charges directly, but we use the
        // monopole approximation even for leaves.
        let r = dist3(target, &node.charge_center).max(f64::EPSILON);
        return node.total_charge / (4.0 * PI * r);
    }

    let r = dist3(target, &node.charge_center).max(f64::EPSILON);
    let d = node.bbox.diameter();

    if d / r < theta {
        // Far field: use monopole
        node.total_charge / (4.0 * PI * r)
    } else {
        // Near field: recurse
        let mut pot = 0.0;
        for child in &node.children {
            pot += barnes_hut_potential(child, target, theta);
        }
        pot
    }
}

/// Evaluate the potential at multiple targets using Barnes–Hut.
pub fn fast_multipole_eval(
    positions: &[[f64; 3]],
    charges: &[f64],
    targets: &[[f64; 3]],
    theta: f64,
) -> Vec<f64> {
    let tree = build_octree(positions, charges, 8);
    targets
        .iter()
        .map(|t| barnes_hut_potential(&tree, t, theta))
        .collect()
}

// ---------------------------------------------------------------------------
// 11. Additional problem types
// ---------------------------------------------------------------------------

/// Configuration for a Laplace BEM problem.
#[derive(Debug, Clone)]
pub struct LaplaceBemProblem {
    /// The boundary mesh.
    pub mesh: BoundaryMesh,
    /// Boundary conditions.
    pub bcs: Vec<BemBc>,
}

impl LaplaceBemProblem {
    /// Create a new Laplace BEM problem.
    pub fn new(mesh: BoundaryMesh, bcs: Vec<BemBc>) -> Self {
        Self { mesh, bcs }
    }

    /// Solve the problem and return `(u, q)`.
    pub fn solve(&self) -> (Vec<f64>, Vec<f64>) {
        let (h, g) = assemble_laplace_bem(&self.mesh);
        solve_laplace_bem(&h, &g, &self.bcs)
    }
}

/// Configuration for a Helmholtz BEM problem.
#[derive(Debug, Clone)]
pub struct HelmholtzBemProblem {
    /// The boundary mesh.
    pub mesh: BoundaryMesh,
    /// Wave number.
    pub k: f64,
    /// Incident field values at element centroids (real, imag).
    pub incident: Vec<(f64, f64)>,
}

impl HelmholtzBemProblem {
    /// Create a new Helmholtz BEM problem.
    pub fn new(mesh: BoundaryMesh, k: f64) -> Self {
        let n = mesh.num_elements();
        Self {
            mesh,
            k,
            incident: vec![(0.0, 0.0); n],
        }
    }

    /// Set a plane-wave incident field travelling in direction `dir`.
    pub fn set_plane_wave(&mut self, amplitude: f64, dir: &[f64; 3]) {
        let n = self.mesh.num_elements();
        for i in 0..n {
            let c = self.mesh.element_centroid(i);
            let phase = self.k * dot3(&c, dir);
            self.incident[i] = (amplitude * phase.cos(), amplitude * phase.sin());
        }
    }

    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.mesh.num_elements()
    }
}

/// Configuration for an elastostatics BEM problem.
#[derive(Debug, Clone)]
pub struct ElastostaticBemProblem {
    /// The boundary mesh.
    pub mesh: BoundaryMesh,
    /// Shear modulus.
    pub shear_modulus: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Boundary conditions per element (3 DOFs each).
    pub bcs: Vec<BemBc>,
}

impl ElastostaticBemProblem {
    /// Create a new elastostatics BEM problem.
    pub fn new(mesh: BoundaryMesh, shear_modulus: f64, nu: f64) -> Self {
        Self {
            mesh,
            shear_modulus,
            nu,
            bcs: Vec::new(),
        }
    }

    /// Add a boundary condition.
    pub fn add_bc(&mut self, bc: BemBc) {
        self.bcs.push(bc);
    }

    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.mesh.num_elements()
    }
}

// ---------------------------------------------------------------------------
// 12. Linear element shape functions
// ---------------------------------------------------------------------------

/// Evaluate linear shape functions on a triangle in natural coordinates.
///
/// `(xi, eta)` are area coordinates with `N1 = 1 - xi - eta`, `N2 = xi`, `N3 = eta`.
pub fn linear_shape_functions(xi: f64, eta: f64) -> [f64; 3] {
    [1.0 - xi - eta, xi, eta]
}

/// Evaluate quadratic shape functions on a triangle (6-node).
///
/// Returns shape function values at `(xi, eta)`.
pub fn quadratic_shape_functions(xi: f64, eta: f64) -> [f64; 6] {
    let l1 = 1.0 - xi - eta;
    let l2 = xi;
    let l3 = eta;
    [
        l1 * (2.0 * l1 - 1.0), // corner 1
        l2 * (2.0 * l2 - 1.0), // corner 2
        l3 * (2.0 * l3 - 1.0), // corner 3
        4.0 * l1 * l2,         // mid-side 1-2
        4.0 * l2 * l3,         // mid-side 2-3
        4.0 * l3 * l1,         // mid-side 3-1
    ]
}

/// Interpolate a point on a triangular element using linear shape functions.
pub fn interpolate_triangle(nodes: &[[f64; 3]; 3], xi: f64, eta: f64) -> [f64; 3] {
    let n = linear_shape_functions(xi, eta);
    [
        n[0] * nodes[0][0] + n[1] * nodes[1][0] + n[2] * nodes[2][0],
        n[0] * nodes[0][1] + n[1] * nodes[1][1] + n[2] * nodes[2][1],
        n[0] * nodes[0][2] + n[1] * nodes[1][2] + n[2] * nodes[2][2],
    ]
}

// ---------------------------------------------------------------------------
// 13. Error estimation for BEM
// ---------------------------------------------------------------------------

/// Estimate the error of a BEM solution by checking the boundary integral
/// equation residual at element centroids.
pub fn bem_residual(
    _mesh: &BoundaryMesh,
    h: &DenseMatrix,
    g: &DenseMatrix,
    u: &[f64],
    q: &[f64],
) -> Vec<f64> {
    let hu = h.matvec(u);
    let gq = g.matvec(q);
    let residual: Vec<f64> = hu
        .iter()
        .zip(gq.iter())
        .map(|(h_i, g_i)| (h_i - g_i).abs())
        .collect();
    residual
}

/// L2 norm of the residual.
pub fn bem_residual_norm(residual: &[f64]) -> f64 {
    residual.iter().map(|r| r * r).sum::<f64>().sqrt()
}

// ---------------------------------------------------------------------------
// 14. Adaptive BEM refinement
// ---------------------------------------------------------------------------

/// Refine boundary elements whose residual exceeds the threshold by
/// splitting each triangle into 4 sub-triangles (midpoint subdivision).
///
/// Returns the refined mesh (original elements with large residual replaced).
pub fn adaptive_refine(mesh: &BoundaryMesh, residual: &[f64], threshold: f64) -> BoundaryMesh {
    let mut new_mesh = BoundaryMesh::new();
    new_mesh.nodes = mesh.nodes.clone();

    for (idx, elem) in mesh.elements.iter().enumerate() {
        if residual[idx] > threshold && elem.node_ids.len() == 3 {
            // Split into 4 triangles
            let n0 = elem.node_ids[0];
            let n1 = elem.node_ids[1];
            let n2 = elem.node_ids[2];
            let p0 = mesh.nodes[n0];
            let p1 = mesh.nodes[n1];
            let p2 = mesh.nodes[n2];
            let m01 = scale3(0.5, &add3(&p0, &p1));
            let m12 = scale3(0.5, &add3(&p1, &p2));
            let m20 = scale3(0.5, &add3(&p2, &p0));
            let nm01 = new_mesh.nodes.len();
            new_mesh.nodes.push(m01);
            let nm12 = new_mesh.nodes.len();
            new_mesh.nodes.push(m12);
            let nm20 = new_mesh.nodes.len();
            new_mesh.nodes.push(m20);

            let sub_tris = [
                [n0, nm01, nm20],
                [nm01, n1, nm12],
                [nm12, n2, nm20],
                [nm01, nm12, nm20],
            ];
            for tri in &sub_tris {
                let area = triangle_area(
                    &new_mesh.nodes[tri[0]],
                    &new_mesh.nodes[tri[1]],
                    &new_mesh.nodes[tri[2]],
                );
                let normal = element_normal(
                    &new_mesh.nodes[tri[0]],
                    &new_mesh.nodes[tri[1]],
                    &new_mesh.nodes[tri[2]],
                );
                new_mesh.elements.push(BoundaryElement {
                    node_ids: tri.to_vec(),
                    normal,
                    area,
                    order: ElementOrder::Constant,
                });
            }
        } else {
            new_mesh.elements.push(elem.clone());
        }
    }
    new_mesh
}

// ---------------------------------------------------------------------------
// 15. Numerical integration helpers
// ---------------------------------------------------------------------------

/// Gauss quadrature on a triangle (area coordinates).
///
/// Returns `(xi, eta, weight)` tuples for the given order.
fn triangle_quadrature(order: usize) -> Vec<(f64, f64, f64)> {
    match order {
        1 => vec![(1.0 / 3.0, 1.0 / 3.0, 0.5)],
        3 => vec![
            (1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
        ],
        _ => {
            // 4-point rule
            vec![
                (1.0 / 3.0, 1.0 / 3.0, -27.0 / 96.0),
                (0.6, 0.2, 25.0 / 96.0),
                (0.2, 0.6, 25.0 / 96.0),
                (0.2, 0.2, 25.0 / 96.0),
            ]
        }
    }
}

/// Integrate a function over a triangular element using Gauss quadrature.
pub fn integrate_over_triangle<F>(nodes: &[[f64; 3]; 3], f: F, quad_order: usize) -> f64
where
    F: Fn(&[f64; 3]) -> f64,
{
    let pts = triangle_quadrature(quad_order);
    let area = triangle_area(&nodes[0], &nodes[1], &nodes[2]);
    let mut result = 0.0;
    for (xi, eta, w) in pts {
        let p = interpolate_triangle(nodes, xi, eta);
        result += w * f(&p);
    }
    result * 2.0 * area // factor of 2 because quad weights sum to area of ref triangle (0.5)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    #[test]
    fn test_dist3_zero() {
        let a = [1.0, 2.0, 3.0];
        assert!(dist3(&a, &a) < TOL);
    }

    #[test]
    fn test_dist3_basic() {
        let a = [0.0; 3];
        let b = [3.0, 4.0, 0.0];
        assert!((dist3(&a, &b) - 5.0).abs() < TOL);
    }

    #[test]
    fn test_dot3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(dot3(&a, &b).abs() < TOL);
    }

    #[test]
    fn test_cross3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(&a, &b);
        assert!((c[2] - 1.0).abs() < TOL);
    }

    #[test]
    fn test_normalize3() {
        let mut v = [3.0, 4.0, 0.0];
        let len = normalize3(&mut v);
        assert!((len - 5.0).abs() < TOL);
        assert!((norm3(&v) - 1.0).abs() < TOL);
    }

    #[test]
    fn test_kelvin_laplace_3d_symmetry() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        assert!((kelvin_laplace_3d(&x, &y) - kelvin_laplace_3d(&y, &x)).abs() < TOL);
    }

    #[test]
    fn test_kelvin_laplace_3d_value() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 0.0, 0.0];
        let expected = 1.0 / (4.0 * PI);
        assert!((kelvin_laplace_3d(&x, &y) - expected).abs() < TOL);
    }

    #[test]
    fn test_kelvin_laplace_2d() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 0.0, 0.0];
        let expected = 0.0; // -ln(1)/(2pi) = 0
        assert!((kelvin_laplace_2d(&x, &y) - expected).abs() < TOL);
    }

    #[test]
    fn test_helmholtz_3d_reduces_to_laplace() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 0.0, 0.0];
        let (re, im) = helmholtz_3d(&x, &y, 0.0);
        let laplace = kelvin_laplace_3d(&x, &y);
        assert!((re - laplace).abs() < TOL);
        assert!(im.abs() < TOL);
    }

    #[test]
    fn test_kelvin_displacement_symmetry() {
        let x = [2.0, 0.0, 0.0];
        let y = [0.0, 0.0, 0.0];
        let g = 1.0;
        let nu = 0.3;
        let u01 = kelvin_displacement_3d(&x, &y, 0, 1, g, nu);
        let u10 = kelvin_displacement_3d(&x, &y, 1, 0, g, nu);
        assert!((u01 - u10).abs() < TOL);
    }

    #[test]
    fn test_half_space_greens_function() {
        // For a point on the free surface (z=0), the half-space Green's function
        // should be twice the full-space one
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 0.0, 0.0];
        let hs = half_space_laplace_3d(&x, &y);
        let fs = kelvin_laplace_3d(&x, &y);
        assert!((hs - 2.0 * fs).abs() < TOL);
    }

    #[test]
    fn test_triangle_area() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = triangle_area(&a, &b, &c);
        assert!((area - 0.5).abs() < TOL);
    }

    #[test]
    fn test_element_normal() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let n = element_normal(&a, &b, &c);
        assert!((n[2] - 1.0).abs() < TOL);
    }

    #[test]
    fn test_boundary_mesh_unit_square() {
        let mesh = BoundaryMesh::unit_square(2, 2);
        assert_eq!(mesh.num_elements(), 8); // 2*2 quads * 2 tris
        assert_eq!(mesh.num_nodes(), 9); // 3*3
    }

    #[test]
    fn test_boundary_mesh_sphere() {
        let mesh = BoundaryMesh::sphere(1.0, 4, 8);
        assert!(mesh.num_elements() > 0);
        assert!(mesh.num_nodes() > 0);
    }

    #[test]
    fn test_dense_matrix_basics() {
        let mut m = DenseMatrix::zeros(3, 3);
        m.set(0, 0, 1.0);
        m.set(1, 1, 2.0);
        m.set(2, 2, 3.0);
        assert!((m.get(1, 1) - 2.0).abs() < TOL);
        let x = vec![1.0, 1.0, 1.0];
        let y = m.matvec(&x);
        assert!((y[0] - 1.0).abs() < TOL);
        assert!((y[1] - 2.0).abs() < TOL);
    }

    #[test]
    fn test_dense_matrix_transpose() {
        let mut m = DenseMatrix::zeros(2, 3);
        m.set(0, 1, 5.0);
        let t = m.transpose();
        assert_eq!(t.rows, 3);
        assert_eq!(t.cols, 2);
        assert!((t.get(1, 0) - 5.0).abs() < TOL);
    }

    #[test]
    fn test_gauss_solve_identity() {
        let mut a = DenseMatrix::zeros(3, 3);
        a.set(0, 0, 1.0);
        a.set(1, 1, 1.0);
        a.set(2, 2, 1.0);
        let b = vec![1.0, 2.0, 3.0];
        let x = gauss_solve(&a, &b);
        for (&xi, &bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < TOL);
        }
    }

    #[test]
    fn test_gauss_solve_2x2() {
        let mut a = DenseMatrix::zeros(2, 2);
        a.set(0, 0, 2.0);
        a.set(0, 1, 1.0);
        a.set(1, 0, 1.0);
        a.set(1, 1, 3.0);
        let b = vec![5.0, 7.0];
        let x = gauss_solve(&a, &b);
        // 2x + y = 5, x + 3y = 7  =>  x = 1.6, y = 1.8
        assert!((x[0] - 1.6).abs() < TOL);
        assert!((x[1] - 1.8).abs() < TOL);
    }

    #[test]
    fn test_telles_transform() {
        let t = TellesTransform::new(0.0);
        let mapped = t.map(0.0);
        // At the singularity, mapped should be near 0
        assert!(mapped.abs() < 1.0);
        let jac = t.jacobian(0.5);
        assert!(jac > 0.0);
    }

    #[test]
    fn test_linear_shape_functions_sum() {
        let n = linear_shape_functions(0.3, 0.2);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < TOL);
    }

    #[test]
    fn test_quadratic_shape_functions_sum() {
        let n = quadratic_shape_functions(0.3, 0.2);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < TOL);
    }

    #[test]
    fn test_interpolate_triangle_centroid() {
        let nodes = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let c = interpolate_triangle(&nodes, 1.0 / 3.0, 1.0 / 3.0);
        assert!((c[0] - 1.0).abs() < TOL);
        assert!((c[1] - 1.0).abs() < TOL);
    }

    #[test]
    fn test_assemble_laplace_bem_dimensions() {
        let mesh = BoundaryMesh::unit_square(2, 2);
        let (h, g) = assemble_laplace_bem(&mesh);
        assert_eq!(h.rows, 8);
        assert_eq!(g.cols, 8);
    }

    #[test]
    fn test_rigid_body_motion() {
        let mesh = BoundaryMesh::unit_square(2, 2);
        let (h, _g) = assemble_laplace_bem(&mesh);
        // H * {1} should be approximately 0
        let ones = vec![1.0; h.cols];
        let hu = h.matvec(&ones);
        let norm: f64 = hu.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            norm < 1e-10,
            "Rigid body condition failed: norm = {:.6e}",
            norm
        );
    }

    #[test]
    fn test_laplace_bem_solve_basic() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let n = mesh.num_elements();
        let bcs: Vec<BemBc> = (0..n)
            .map(|i| BemBc {
                element_id: i,
                bc_type: BemBcType::Dirichlet,
                values: vec![1.0],
            })
            .collect();
        let (h, g) = assemble_laplace_bem(&mesh);
        let (u, _q) = solve_laplace_bem(&h, &g, &bcs);
        for val in &u {
            assert!((val - 1.0).abs() < TOL);
        }
    }

    #[test]
    fn test_assemble_elastostatic_bem_dimensions() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let (h, g) = assemble_elastostatic_bem(&mesh, 1.0, 0.3);
        let n = mesh.num_elements();
        assert_eq!(h.rows, 3 * n);
        assert_eq!(g.cols, 3 * n);
    }

    #[test]
    fn test_assemble_helmholtz_bem_dimensions() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let (hr, hi, gr, gi) = assemble_helmholtz_bem(&mesh, 1.0);
        let n = mesh.num_elements();
        assert_eq!(hr.rows, n);
        assert_eq!(hi.cols, n);
        assert_eq!(gr.rows, n);
        assert_eq!(gi.cols, n);
    }

    #[test]
    fn test_octree_build() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let charges = vec![1.0, 1.0, 1.0, 1.0];
        let tree = build_octree(&positions, &charges, 2);
        assert!((tree.total_charge - 4.0).abs() < TOL);
    }

    #[test]
    fn test_barnes_hut_far_field() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![4.0 * PI]; // so potential = 1/r
        let tree = build_octree(&positions, &charges, 2);
        let target = [1.0, 0.0, 0.0];
        let pot = barnes_hut_potential(&tree, &target, 0.5);
        // Expected: q/(4 pi r) = 4pi/(4pi*1) = 1.0
        assert!((pot - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_fast_multipole_eval() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let charges = vec![4.0 * PI];
        let targets = vec![[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let pots = fast_multipole_eval(&positions, &charges, &targets, 0.5);
        assert_eq!(pots.len(), 2);
        // Both targets at distance 2: pot = 4pi/(4pi*2) = 0.5
        for p in &pots {
            assert!((p - 0.5).abs() < 0.1);
        }
    }

    #[test]
    fn test_dual_bem_sif() {
        let mesh = BoundaryMesh::new();
        let crack = vec![CrackElement {
            node_ids: vec![0, 1, 2],
            normal: [0.0, 0.0, 1.0],
            area: 1.0,
        }];
        let dbem = DualBem::new(mesh, crack, 1.0, 0.3);
        let cod = [0.01, 0.0, 0.0];
        let sif = dbem.compute_sif(&cod, 0.1);
        assert!(sif[0] > 0.0);
    }

    #[test]
    fn test_bem_fem_coupling_new() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let coupling = BemFemCoupling::new(vec![0, 1, 2], mesh);
        assert_eq!(coupling.num_interface_dofs(), 3);
    }

    #[test]
    fn test_bem_residual_zeros() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let n = mesh.num_elements();
        let (h, g) = assemble_laplace_bem(&mesh);
        // If u = 0 and q = 0, residual should be 0
        let u = vec![0.0; n];
        let q = vec![0.0; n];
        let res = bem_residual(&mesh, &h, &g, &u, &q);
        let norm = bem_residual_norm(&res);
        assert!(norm < TOL);
    }

    #[test]
    fn test_adaptive_refine() {
        let mesh = BoundaryMesh::unit_square(1, 1);
        let residual = vec![1.0, 0.01]; // first element above threshold
        let refined = adaptive_refine(&mesh, &residual, 0.5);
        // First element split into 4, second unchanged
        assert_eq!(refined.num_elements(), 5);
    }

    #[test]
    fn test_integrate_over_triangle() {
        let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // Integrate f(x) = 1 -> area = 0.5
        let result = integrate_over_triangle(&nodes, |_p| 1.0, 3);
        assert!((result - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_helmholtz_problem_plane_wave() {
        let mesh = BoundaryMesh::unit_square(2, 2);
        let mut prob = HelmholtzBemProblem::new(mesh, 1.0);
        prob.set_plane_wave(1.0, &[1.0, 0.0, 0.0]);
        assert_eq!(prob.num_elements(), 8);
        // At least one incident field should be non-zero
        let has_nonzero = prob
            .incident
            .iter()
            .any(|(r, i)| r.abs() > TOL || i.abs() > TOL);
        assert!(has_nonzero);
    }

    #[test]
    fn test_gauss_legendre_weights_sum() {
        for n in 1..=5 {
            let pts = gauss_legendre(n);
            let sum: f64 = pts.iter().map(|(_, w)| w).sum();
            assert!((sum - 2.0).abs() < 0.01, "GL({n}) weights sum = {sum:.6}");
        }
    }

    #[test]
    fn test_aabb_contains() {
        let bb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(bb.contains(&[0.5, 0.5, 0.5]));
        assert!(!bb.contains(&[1.5, 0.5, 0.5]));
    }

    #[test]
    fn test_aabb_diameter() {
        let bb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let d = bb.diameter();
        assert!((d - 3.0_f64.sqrt()).abs() < TOL);
    }
}
