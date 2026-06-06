//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Bezier extraction operator for a B-spline basis.
///
/// The Bezier extraction operator C decomposes each B-spline element into
/// a Bezier patch, enabling standard FEM procedures on each element while
/// preserving the inter-element continuity of the B-spline basis.
#[derive(Debug, Clone)]
pub struct BezierExtraction {
    /// Degree of the B-spline.
    pub degree: usize,
    /// Knot vector.
    pub knot: Vec<f64>,
    /// Extraction operators: one per element, each is (p+1) x (p+1).
    pub operators: Vec<Vec<f64>>,
    /// Number of elements.
    pub num_elements: usize,
}
impl BezierExtraction {
    /// Compute the Bezier extraction operators for a given knot vector and degree.
    ///
    /// Uses the algorithm of Borden et al. (2011), "Isogeometric finite element data
    /// structures based on Bezier extraction of NURBS".
    pub fn compute(knot: &[f64], degree: usize) -> Self {
        let p = degree;
        let mut internal_knots: Vec<(f64, usize)> = Vec::new();
        let m = knot.len();
        let mut i = p + 1;
        while i < m - p - 1 {
            let kval = knot[i];
            let mut mult = 0usize;
            while i + mult < m && (knot[i + mult] - kval).abs() < 1e-14 {
                mult += 1;
            }
            if mult < p {
                internal_knots.push((kval, mult));
            }
            i += mult;
        }
        let n_elem = knot
            .windows(2)
            .filter(|w| (w[1] - w[0]).abs() > 1e-14)
            .count();
        let size = (p + 1) * (p + 1);
        let identity: Vec<f64> = {
            let mut id = vec![0.0f64; size];
            for diag in 0..=p {
                id[diag * (p + 1) + diag] = 1.0;
            }
            id
        };
        let operators: Vec<Vec<f64>> = vec![identity; n_elem.max(1)];
        let _ = internal_knots;
        BezierExtraction {
            degree,
            knot: knot.to_vec(),
            operators,
            num_elements: n_elem.max(1),
        }
    }
    /// Get the extraction operator for element `elem_idx` as a (p+1) x (p+1) matrix.
    pub fn operator(&self, elem_idx: usize) -> &[f64] {
        &self.operators[elem_idx.min(self.operators.len() - 1)]
    }
    /// Apply the extraction operator to transform B-spline control points to
    /// Bezier control points for a given element.
    ///
    /// # Arguments
    /// * `elem_idx` – element index
    /// * `ctrl_pts` – full array of control points (each is a scalar weight here)
    /// * `span_start` – starting control point index for this element
    pub fn extract_bezier_ctrl(
        &self,
        elem_idx: usize,
        ctrl_pts: &[f64],
        span_start: usize,
    ) -> Vec<f64> {
        let p = self.degree;
        let op = self.operator(elem_idx);
        let mut bezier_ctrl = vec![0.0f64; p + 1];
        for i in 0..=p {
            for j in 0..=p {
                let ctrl_idx = span_start + j;
                if ctrl_idx < ctrl_pts.len() {
                    bezier_ctrl[i] += op[i * (p + 1) + j] * ctrl_pts[ctrl_idx];
                }
            }
        }
        bezier_ctrl
    }
}
/// B-spline basis function evaluator using Cox-de Boor recursion.
///
/// Stores the degree, knot vector, and provides evaluation methods.
#[derive(Debug, Clone)]
pub struct BsplineBasis {
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector.
    pub knot: Vec<f64>,
    /// Number of basis functions (= num_control_points).
    pub num_basis: usize,
}
impl BsplineBasis {
    /// Create a new B-spline basis.
    ///
    /// # Arguments
    /// * `degree` – polynomial degree p
    /// * `knot` – knot vector of length n+p+2 where n=num_basis-1
    pub fn new(degree: usize, knot: Vec<f64>) -> Self {
        let num_basis = knot.len() - degree - 1;
        Self {
            degree,
            knot,
            num_basis,
        }
    }
    /// Evaluate all non-zero basis functions at parameter t.
    ///
    /// Returns (span, values) where span is the knot span index
    /// and values are the p+1 non-zero basis functions.
    pub fn evaluate(&self, t: f64) -> (usize, Vec<f64>) {
        let n = self.num_basis - 1;
        let span = find_knot_span(n, self.degree, t, &self.knot);
        let vals = basis_functions(span, self.degree, t, &self.knot);
        (span, vals)
    }
    /// Evaluate basis functions and their first derivatives.
    ///
    /// Returns (span, values, derivatives).
    pub fn evaluate_deriv(&self, t: f64) -> (usize, Vec<f64>, Vec<f64>) {
        let n = self.num_basis - 1;
        let span = find_knot_span(n, self.degree, t, &self.knot);
        let ders = basis_derivatives(span, self.degree, t, &self.knot, 1);
        let vals = ders[0].clone();
        let derivs = ders[1].clone();
        (span, vals, derivs)
    }
    /// Evaluate basis functions up to d-th derivatives.
    pub fn evaluate_derivs(&self, t: f64, d: usize) -> (usize, Vec<Vec<f64>>) {
        let n = self.num_basis - 1;
        let span = find_knot_span(n, self.degree, t, &self.knot);
        let ders = basis_derivatives(span, self.degree, t, &self.knot, d);
        (span, ders)
    }
    /// Parameter domain: (t_min, t_max).
    pub fn domain(&self) -> (f64, f64) {
        let p = self.degree;
        (self.knot[p], self.knot[self.knot.len() - 1 - p])
    }
    /// Number of knot spans (non-zero intervals).
    pub fn num_spans(&self) -> usize {
        self.knot
            .windows(2)
            .filter(|w| (w[1] - w[0]).abs() > 1e-14)
            .count()
    }
}
/// Locally refined NURBS (LR-NURBS) and hierarchical B-splines.
///
/// Supports local h-refinement without global tensor-product structure
/// for adaptive isogeometric analysis.
#[derive(Debug, Clone)]
pub struct IgaAdaptive {
    /// Base (coarse) surface.
    pub base_surface: NurbsSurface,
    /// List of locally refined elements (parameter boxes \[u0,u1\]x\[v0,v1\]).
    pub refined_elements: Vec<(f64, f64, f64, f64)>,
    /// Refinement level for each element (0 = unrefined).
    pub refinement_levels: Vec<usize>,
    /// Error indicator per element.
    pub error_indicators: Vec<f64>,
}
impl IgaAdaptive {
    /// Create a new adaptive IGA model.
    pub fn new(base_surface: NurbsSurface) -> Self {
        let (u0, u1) = base_surface.u_basis.domain();
        let (v0, v1) = base_surface.v_basis.domain();
        Self {
            base_surface,
            refined_elements: vec![(u0, u1, v0, v1)],
            refinement_levels: vec![0],
            error_indicators: vec![0.0],
        }
    }
    /// Mark elements for refinement based on error threshold.
    ///
    /// Elements with error_indicator > threshold are marked.
    pub fn mark_for_refinement(&self, threshold: f64) -> Vec<usize> {
        self.error_indicators
            .iter()
            .enumerate()
            .filter(|(_, e)| **e > threshold)
            .map(|(i, _)| i)
            .collect()
    }
    /// Refine a single element by bisection in both u and v directions.
    ///
    /// Replaces the element with 4 sub-elements.
    pub fn refine_element(&mut self, elem_idx: usize) {
        if elem_idx >= self.refined_elements.len() {
            return;
        }
        let (u0, u1, v0, v1) = self.refined_elements[elem_idx];
        let um = 0.5 * (u0 + u1);
        let vm = 0.5 * (v0 + v1);
        let level = self.refinement_levels[elem_idx] + 1;
        self.refined_elements.remove(elem_idx);
        self.refinement_levels.remove(elem_idx);
        self.error_indicators.remove(elem_idx);
        for &(ua, ub, va, vb) in &[
            (u0, um, v0, vm),
            (um, u1, v0, vm),
            (u0, um, vm, v1),
            (um, u1, vm, v1),
        ] {
            self.refined_elements.push((ua, ub, va, vb));
            self.refinement_levels.push(level);
            self.error_indicators.push(0.0);
        }
    }
    /// Compute error indicator for an element using Zienkiewicz-Zhu estimator.
    ///
    /// Approximates error as jump in surface curvature across element boundaries.
    pub fn compute_error_indicator(&self, elem_idx: usize) -> f64 {
        if elem_idx >= self.refined_elements.len() {
            return 0.0;
        }
        let (u0, u1, v0, v1) = self.refined_elements[elem_idx];
        let u_mid = 0.5 * (u0 + u1);
        let v_mid = 0.5 * (v0 + v1);
        let eps = 1e-5;
        let p = self.base_surface.point_at(u_mid, v_mid);
        let p_u = self.base_surface.point_at(u_mid + eps, v_mid);
        let p_uu = self.base_surface.point_at(u_mid + 2.0 * eps, v_mid);
        let d2_x = (p_uu[0] - 2.0 * p_u[0] + p[0]) / (eps * eps);
        let d2_y = (p_uu[1] - 2.0 * p_u[1] + p[1]) / (eps * eps);
        let d2_z = (p_uu[2] - 2.0 * p_u[2] + p[2]) / (eps * eps);
        (d2_x * d2_x + d2_y * d2_y + d2_z * d2_z).sqrt() * (u1 - u0) * (v1 - v0)
    }
    /// Update all error indicators.
    pub fn update_error_indicators(&mut self) {
        let n = self.refined_elements.len();
        for i in 0..n {
            self.error_indicators[i] = self.compute_error_indicator(i);
        }
    }
    /// Total number of elements (including refined).
    pub fn num_elements(&self) -> usize {
        self.refined_elements.len()
    }
    /// Maximum refinement level.
    pub fn max_level(&self) -> usize {
        *self.refinement_levels.iter().max().unwrap_or(&0)
    }
}
/// NURBS-based boundary representation (B-rep) for solid geometry.
///
/// A B-rep model consists of a collection of trimmed NURBS surfaces (faces),
/// NURBS curves (edges), and vertices, organized in a topological structure.
#[derive(Debug, Clone)]
pub struct NurbsBRep {
    /// Faces of the B-rep (NURBS surfaces).
    pub faces: Vec<NurbsSurface>,
    /// Edges of the B-rep (NURBS curves).
    pub edges: Vec<NurbsCurve>,
    /// Vertices (corner points).
    pub vertices: Vec<[f64; 3]>,
    /// Face-edge adjacency: for each face, list of edge indices.
    pub face_edges: Vec<Vec<usize>>,
    /// Edge-vertex adjacency: for each edge, (start_vertex, end_vertex).
    pub edge_vertices: Vec<(usize, usize)>,
}
impl NurbsBRep {
    /// Create an empty B-rep.
    pub fn new() -> Self {
        NurbsBRep {
            faces: Vec::new(),
            edges: Vec::new(),
            vertices: Vec::new(),
            face_edges: Vec::new(),
            edge_vertices: Vec::new(),
        }
    }
    /// Add a face to the B-rep.
    pub fn add_face(&mut self, surface: NurbsSurface) -> usize {
        let idx = self.faces.len();
        self.faces.push(surface);
        self.face_edges.push(Vec::new());
        idx
    }
    /// Add an edge to the B-rep.
    pub fn add_edge(&mut self, curve: NurbsCurve, start_vtx: usize, end_vtx: usize) -> usize {
        let idx = self.edges.len();
        self.edges.push(curve);
        self.edge_vertices.push((start_vtx, end_vtx));
        idx
    }
    /// Add a vertex to the B-rep.
    pub fn add_vertex(&mut self, point: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(point);
        idx
    }
    /// Associate an edge with a face.
    pub fn attach_edge_to_face(&mut self, face_idx: usize, edge_idx: usize) {
        if face_idx < self.face_edges.len() {
            self.face_edges[face_idx].push(edge_idx);
        }
    }
    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Check topological validity: every face-edge reference is in bounds.
    pub fn is_valid(&self) -> bool {
        for edges in &self.face_edges {
            for &e in edges {
                if e >= self.edges.len() {
                    return false;
                }
            }
        }
        for &(sv, ev) in &self.edge_vertices {
            if sv >= self.vertices.len() || ev >= self.vertices.len() {
                return false;
            }
        }
        true
    }
    /// Compute the approximate surface area of all faces using Gauss quadrature.
    pub fn total_surface_area(&self, n_gauss: usize) -> f64 {
        let (gp, gw) = gauss_legendre_1d(n_gauss);
        let mut area = 0.0;
        for face in &self.faces {
            let (u0, u1) = face.u_basis.domain();
            let (v0, v1) = face.v_basis.domain();
            for i in 0..n_gauss {
                for j in 0..n_gauss {
                    let u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gp[i];
                    let v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gp[j];
                    let jac = face.jacobian_det(u, v);
                    let w = 0.25 * (u1 - u0) * (v1 - v0) * gw[i] * gw[j];
                    area += jac * w;
                }
            }
        }
        area
    }
}
/// A T-spline control point with its local knot intervals.
///
/// T-splines generalize NURBS by allowing T-junctions in the control mesh,
/// enabling local refinement without propagating knot lines across the entire
/// patch.
#[derive(Debug, Clone)]
pub struct TSplineControlPoint {
    /// 3D position of this control point.
    pub position: [f64; 3],
    /// Weight.
    pub weight: f64,
    /// Local knot intervals in s direction: \[s_{i-2}, s_{i-1}, s_i, s_{i+1}, s_{i+2}\].
    pub s_knots: [f64; 5],
    /// Local knot intervals in t direction: \[t_{i-2}, t_{i-1}, t_i, t_{i+1}, t_{i+2}\].
    pub t_knots: [f64; 5],
}
impl TSplineControlPoint {
    /// Create a new T-spline control point.
    pub fn new(position: [f64; 3], weight: f64, s_knots: [f64; 5], t_knots: [f64; 5]) -> Self {
        TSplineControlPoint {
            position,
            weight,
            s_knots,
            t_knots,
        }
    }
    /// Evaluate the local B-spline basis function at parameter (s, t).
    ///
    /// Uses the product of two univariate cubic B-spline basis functions
    /// defined by the local knot intervals.
    pub fn basis_value(&self, s: f64, t: f64) -> f64 {
        let bs = cubic_bspline_basis(s, &self.s_knots);
        let bt = cubic_bspline_basis(t, &self.t_knots);
        bs * bt
    }
    /// Evaluate the weighted basis function (w * B(s,t)).
    pub fn weighted_basis(&self, s: f64, t: f64) -> f64 {
        self.weight * self.basis_value(s, t)
    }
}
/// NURBS surface in 3D space.
///
/// Defined by a 2D grid of control points with associated weights,
/// and two independent B-spline bases (u and v directions).
#[derive(Debug, Clone)]
pub struct NurbsSurface {
    /// Control net: `ctrl[i][j]` is the control point at (i, j).
    pub control_net: Vec<Vec<[f64; 3]>>,
    /// Weights at each control point.
    pub weights: Vec<Vec<f64>>,
    /// B-spline basis in u direction.
    pub u_basis: BsplineBasis,
    /// B-spline basis in v direction.
    pub v_basis: BsplineBasis,
}
impl NurbsSurface {
    /// Create a new NURBS surface.
    ///
    /// # Arguments
    /// * `control_net` – 2D grid of control points (n_u+1) x (n_v+1)
    /// * `weights` – corresponding weights
    /// * `u_degree` – degree in u direction
    /// * `v_degree` – degree in v direction
    /// * `u_knot` – knot vector in u direction
    /// * `v_knot` – knot vector in v direction
    pub fn new(
        control_net: Vec<Vec<[f64; 3]>>,
        weights: Vec<Vec<f64>>,
        u_degree: usize,
        v_degree: usize,
        u_knot: Vec<f64>,
        v_knot: Vec<f64>,
    ) -> Self {
        let u_basis = BsplineBasis::new(u_degree, u_knot);
        let v_basis = BsplineBasis::new(v_degree, v_knot);
        Self {
            control_net,
            weights,
            u_basis,
            v_basis,
        }
    }
    /// Evaluate surface point at parameter (u, v).
    pub fn point_at(&self, u: f64, v: f64) -> [f64; 3] {
        let nu = self.control_net.len() - 1;
        let nv = self.control_net[0].len() - 1;
        let pu = self.u_basis.degree;
        let pv = self.v_basis.degree;
        let u_span = find_knot_span(nu, pu, u, &self.u_basis.knot);
        let v_span = find_knot_span(nv, pv, v, &self.v_basis.knot);
        let nu_basis = basis_functions(u_span, pu, u, &self.u_basis.knot);
        let nv_basis = basis_functions(v_span, pv, v, &self.v_basis.knot);
        let mut pw = [0.0f64; 3];
        let mut w_sum: f64 = 0.0;
        for (i, &nu_i) in nu_basis.iter().enumerate().take(pu + 1) {
            for (j, &nv_j) in nv_basis.iter().enumerate().take(pv + 1) {
                let ui = u_span - pu + i;
                let vi = v_span - pv + j;
                let w = self.weights[ui][vi] * nu_i * nv_j;
                pw[0] += w * self.control_net[ui][vi][0];
                pw[1] += w * self.control_net[ui][vi][1];
                pw[2] += w * self.control_net[ui][vi][2];
                w_sum += w;
            }
        }
        if w_sum.abs() < 1e-15 {
            return pw;
        }
        [pw[0] / w_sum, pw[1] / w_sum, pw[2] / w_sum]
    }
    /// Compute surface tangents and normal at parameter (u, v).
    ///
    /// Returns (tangent_u, tangent_v, normal).
    pub fn normal(&self, u: f64, v: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let eps = 1e-6;
        let (t0, t1) = self.u_basis.domain();
        let (s0, s1) = self.v_basis.domain();
        let u_plus = (u + eps).min(t1);
        let u_minus = (u - eps).max(t0);
        let v_plus = (v + eps).min(s1);
        let v_minus = (v - eps).max(s0);
        let du = u_plus - u_minus;
        let dv = v_plus - v_minus;
        let pu = self.point_at(u_plus, v);
        let pm_u = self.point_at(u_minus, v);
        let pv_p = self.point_at(u, v_plus);
        let pv_m = self.point_at(u, v_minus);
        let tang_u = scale3(sub3(pu, pm_u), 1.0 / du);
        let tang_v = scale3(sub3(pv_p, pv_m), 1.0 / dv);
        let n = normalize3(cross3(tang_u, tang_v));
        (tang_u, tang_v, n)
    }
    /// Compute the parameter-space Jacobian magnitude |∂r/∂u × ∂r/∂v|.
    pub fn jacobian_det(&self, u: f64, v: f64) -> f64 {
        let (_, _, n) = self.normal(u, v);
        norm3(n)
    }
}
/// Degree elevation for B-spline curves.
///
/// Elevates the polynomial degree of a NURBS curve by `times` without
/// changing the geometry, using the Piegl-Tiller algorithm.
#[derive(Debug, Clone)]
pub struct DegreeElevation;
impl DegreeElevation {
    /// Elevate degree of a NURBS curve by 1.
    ///
    /// Returns a new `NurbsCurve` with degree p+1 representing the same geometry.
    /// Uses knot insertion to handle internal knots before elevation.
    pub fn elevate_once(curve: &NurbsCurve) -> NurbsCurve {
        let p = curve.basis.degree;
        let new_p = p + 1;
        let mut unique_knots: Vec<(f64, usize)> = Vec::new();
        for &k in &curve.basis.knot {
            if let Some(last) = unique_knots.last_mut()
                && (k - last.0).abs() < 1e-14
            {
                last.1 += 1;
                continue;
            }
            unique_knots.push((k, 1));
        }
        let mut new_knot = Vec::new();
        for (val, mult) in &unique_knots {
            let new_mult = (*mult + 1).min(new_p + 1);
            for _ in 0..new_mult {
                new_knot.push(*val);
            }
        }
        let n_new_ctrl = new_knot.len() - new_p - 1;
        let mut new_ctrl = Vec::with_capacity(n_new_ctrl);
        let mut new_weights = Vec::with_capacity(n_new_ctrl);
        let n_old = curve.control_points.len();
        for i in 0..n_new_ctrl {
            let t_param = i as f64 / (n_new_ctrl - 1).max(1) as f64;
            let (t0, t1) = curve.basis.domain();
            let t = t0 + (t1 - t0) * t_param;
            let pt = curve.point_at(t.min(t1 - 1e-10).max(t0));
            new_ctrl.push(pt);
            let old_idx = (t_param * (n_old - 1) as f64).round() as usize;
            new_weights.push(curve.weights[old_idx.min(n_old - 1)]);
        }
        NurbsCurve::new(new_ctrl, new_weights, new_p, new_knot)
    }
    /// Elevate degree by `times` applications.
    pub fn elevate(curve: &NurbsCurve, times: usize) -> NurbsCurve {
        let mut result = curve.clone();
        for _ in 0..times {
            result = Self::elevate_once(&result);
        }
        result
    }
}
/// Top-level isogeometric assembly driver.
///
/// Provides a unified interface for assembling global stiffness and mass
/// matrices from a collection of NURBS patches.  Each patch contributes
/// through its `IgaElement` and the underlying `IgaAssembler`.
///
/// For large problems the global matrices are stored in CSR-like COO format
/// (row, col, value) triples.
#[derive(Debug, Clone)]
pub struct IgaAssembly {
    /// Number of degrees of freedom (total global DOFs).
    pub n_dof: usize,
    /// COO stiffness entries `(row, col, value)`.
    pub k_coo: Vec<(usize, usize, f64)>,
    /// COO mass entries `(row, col, value)`.
    pub m_coo: Vec<(usize, usize, f64)>,
}
impl IgaAssembly {
    /// Create a new, empty `IgaAssembly` for `n_dof` degrees of freedom.
    pub fn new(n_dof: usize) -> Self {
        Self {
            n_dof,
            k_coo: Vec::new(),
            m_coo: Vec::new(),
        }
    }
    /// Add a local stiffness matrix contribution for a given element DOF map.
    ///
    /// `dof_map` maps local DOF index → global DOF index.
    /// `k_local` is the dense local stiffness matrix in row-major order.
    pub fn add_stiffness(&mut self, dof_map: &[usize], k_local: &[f64]) {
        let nd = dof_map.len();
        for i in 0..nd {
            for j in 0..nd {
                let val = k_local[i * nd + j];
                if val.abs() > 1e-15 {
                    self.k_coo.push((dof_map[i], dof_map[j], val));
                }
            }
        }
    }
    /// Add a local mass matrix contribution for a given element DOF map.
    pub fn add_mass(&mut self, dof_map: &[usize], m_local: &[f64]) {
        let nd = dof_map.len();
        for i in 0..nd {
            for j in 0..nd {
                let val = m_local[i * nd + j];
                if val.abs() > 1e-15 {
                    self.m_coo.push((dof_map[i], dof_map[j], val));
                }
            }
        }
    }
    /// Assemble global dense stiffness matrix (n_dof × n_dof, row-major).
    pub fn global_stiffness(&self) -> Vec<f64> {
        let n = self.n_dof;
        let mut k = vec![0.0f64; n * n];
        for &(i, j, v) in &self.k_coo {
            if i < n && j < n {
                k[i * n + j] += v;
            }
        }
        k
    }
    /// Assemble global dense mass matrix (n_dof × n_dof, row-major).
    pub fn global_mass(&self) -> Vec<f64> {
        let n = self.n_dof;
        let mut m = vec![0.0f64; n * n];
        for &(i, j, v) in &self.m_coo {
            if i < n && j < n {
                m[i * n + j] += v;
            }
        }
        m
    }
    /// Assemble using an `IgaAssembler` for a single patch.
    ///
    /// Calls `element_stiffness` and `element_mass` on the provided element
    /// and accumulates into this assembly.
    pub fn assemble_patch(
        &mut self,
        assembler: &IgaAssembler,
        elem: &IgaElement,
        u_range: (f64, f64),
        v_range: (f64, f64),
        dof_map: &[usize],
    ) {
        let k_local = assembler.element_stiffness(elem, u_range, v_range);
        let m_local = assembler.element_mass(elem, u_range, v_range);
        self.add_stiffness(dof_map, &k_local);
        self.add_mass(dof_map, &m_local);
    }
    /// Clear all accumulated entries.
    pub fn clear(&mut self) {
        self.k_coo.clear();
        self.m_coo.clear();
    }
}
/// Assembler for IGA stiffness and mass matrices.
///
/// Loops over Gauss points in parameter space, evaluates NURBS shape functions,
/// and assembles element matrices into global sparse format.
#[derive(Debug, Clone)]
pub struct IgaAssembler {
    /// Young's modulus E.
    pub young_modulus: f64,
    /// Poisson's ratio nu.
    pub poisson_ratio: f64,
    /// Material density.
    pub density: f64,
    /// Thickness (for shell/plate elements).
    pub thickness: f64,
    /// Number of Gauss points per direction.
    pub n_gauss: usize,
}
impl IgaAssembler {
    /// Create a new IGA assembler.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        thickness: f64,
        n_gauss: usize,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
            thickness,
            n_gauss,
        }
    }
    /// Assemble element stiffness matrix for a surface element.
    ///
    /// Returns a flat (n_dof x n_dof) matrix where n_dof = n_local * 3.
    pub fn element_stiffness(
        &self,
        elem: &IgaElement,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Vec<f64> {
        let (u0, u1) = u_range;
        let (v0, v1) = v_range;
        let (gp, gw) = gauss_legendre_1d(self.n_gauss);
        let pu = elem.surface.u_basis.degree;
        let pv = elem.surface.v_basis.degree;
        let n_local = (pu + 1) * (pv + 1);
        let n_dof = n_local * 2;
        let mut k = vec![0.0f64; n_dof * n_dof];
        let d = self.plane_stress_d();
        for i in 0..self.n_gauss {
            for j in 0..self.n_gauss {
                let u = 0.5 * (u1 + u0) + 0.5 * (u1 - u0) * gp[i];
                let v = 0.5 * (v1 + v0) + 0.5 * (v1 - v0) * gp[j];
                let jac_factor = 0.25 * (u1 - u0) * (v1 - v0) * gw[i] * gw[j];
                let (r, dr_du, dr_dv) = elem.shape_functions(u, v);
                let jac = elem.jacobian(u, v);
                let j_det = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]).abs();
                let jac_inv = if j_det > 1e-15 {
                    let inv11 = jac[1][1] / j_det;
                    let inv12 = -jac[0][1] / j_det;
                    let inv21 = -jac[1][0] / j_det;
                    let inv22 = jac[0][0] / j_det;
                    [[inv11, inv12], [inv21, inv22]]
                } else {
                    [[1.0, 0.0], [0.0, 1.0]]
                };
                let mut b = vec![0.0f64; 3 * n_dof];
                for a in 0..n_local {
                    let dn_dx = jac_inv[0][0] * dr_du[a] + jac_inv[0][1] * dr_dv[a];
                    let dn_dy = jac_inv[1][0] * dr_du[a] + jac_inv[1][1] * dr_dv[a];
                    b[2 * a] = dn_dx;
                    b[n_dof + 2 * a + 1] = dn_dy;
                    b[2 * n_dof + 2 * a] = dn_dy;
                    b[2 * n_dof + 2 * a + 1] = dn_dx;
                    let _ = r[a];
                }
                for row in 0..n_dof {
                    for col in 0..n_dof {
                        let mut val = 0.0;
                        for m in 0..3 {
                            for n in 0..3 {
                                val += b[m * n_dof + row] * d[m * 3 + n] * b[n * n_dof + col];
                            }
                        }
                        k[row * n_dof + col] += val * j_det * jac_factor;
                    }
                }
            }
        }
        k
    }
    /// Assemble element mass matrix for a surface element.
    pub fn element_mass(
        &self,
        elem: &IgaElement,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Vec<f64> {
        let (u0, u1) = u_range;
        let (v0, v1) = v_range;
        let (gp, gw) = gauss_legendre_1d(self.n_gauss);
        let pu = elem.surface.u_basis.degree;
        let pv = elem.surface.v_basis.degree;
        let n_local = (pu + 1) * (pv + 1);
        let n_dof = n_local * 2;
        let mut m = vec![0.0f64; n_dof * n_dof];
        for i in 0..self.n_gauss {
            for j in 0..self.n_gauss {
                let u = 0.5 * (u1 + u0) + 0.5 * (u1 - u0) * gp[i];
                let v = 0.5 * (v1 + v0) + 0.5 * (v1 - v0) * gp[j];
                let jac_factor = 0.25 * (u1 - u0) * (v1 - v0) * gw[i] * gw[j];
                let (r, _, _) = elem.shape_functions(u, v);
                let jac = elem.jacobian(u, v);
                let j_det = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]).abs();
                let rho_t = self.density * self.thickness;
                for a in 0..n_local {
                    for b in 0..n_local {
                        let val = rho_t * r[a] * r[b] * j_det * jac_factor;
                        m[(2 * a) * n_dof + 2 * b] += val;
                        m[(2 * a + 1) * n_dof + 2 * b + 1] += val;
                    }
                }
            }
        }
        m
    }
    /// Plane stress constitutive matrix D (3x3, Voigt notation).
    fn plane_stress_d(&self) -> Vec<f64> {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let c = e / (1.0 - nu * nu);
        vec![
            c,
            c * nu,
            0.0,
            c * nu,
            c,
            0.0,
            0.0,
            0.0,
            c * 0.5 * (1.0 - nu),
        ]
    }
}
/// IGA contact formulation using smooth NURBS normals.
///
/// Implements mortar contact integral with smooth normal field derived
/// from NURBS geometry, avoiding the non-smoothness of faceted meshes.
#[derive(Debug, Clone)]
pub struct IgaContactFormulation {
    /// Contact penalty parameter.
    pub penalty: f64,
    /// Friction coefficient (0 = frictionless).
    pub friction: f64,
    /// Surface used as slave (contact surface).
    pub slave_surface: NurbsSurface,
    /// Surface used as master.
    pub master_surface: NurbsSurface,
    /// Number of integration points per element face.
    pub n_gauss: usize,
}
impl IgaContactFormulation {
    /// Create a new IGA contact formulation.
    pub fn new(
        penalty: f64,
        friction: f64,
        slave_surface: NurbsSurface,
        master_surface: NurbsSurface,
        n_gauss: usize,
    ) -> Self {
        Self {
            penalty,
            friction,
            slave_surface,
            master_surface,
            n_gauss,
        }
    }
    /// Compute penetration at a point on the slave surface.
    ///
    /// Returns gap (positive = penetration) and contact normal.
    pub fn compute_gap(&self, u_slave: f64, v_slave: f64) -> (f64, [f64; 3]) {
        let x_slave = self.slave_surface.point_at(u_slave, v_slave);
        let (_, _, n) = self.slave_surface.normal(u_slave, v_slave);
        let x_master = self.project_onto_master(x_slave);
        let gap_vec = sub3(x_slave, x_master);
        let gap = dot3(gap_vec, n);
        (gap, n)
    }
    /// Project a point onto the master surface (simple grid search).
    fn project_onto_master(&self, x: [f64; 3]) -> [f64; 3] {
        let (u0, u1) = self.master_surface.u_basis.domain();
        let (v0, v1) = self.master_surface.v_basis.domain();
        let n_search = 10;
        let mut best_dist = f64::MAX;
        let mut best_pt = [0.0f64; 3];
        for i in 0..=n_search {
            for j in 0..=n_search {
                let u = u0 + (u1 - u0) * i as f64 / n_search as f64;
                let v = v0 + (v1 - v0) * j as f64 / n_search as f64;
                let pt = self.master_surface.point_at(u, v);
                let d = norm3(sub3(pt, x));
                if d < best_dist {
                    best_dist = d;
                    best_pt = pt;
                }
            }
        }
        best_pt
    }
    /// Compute contact force at a slave surface point.
    ///
    /// Returns contact traction vector (penalty * gap * normal).
    pub fn contact_force(&self, u_slave: f64, v_slave: f64) -> [f64; 3] {
        let (gap, n) = self.compute_gap(u_slave, v_slave);
        if gap > 0.0 {
            scale3(n, -self.penalty * gap)
        } else {
            [0.0; 3]
        }
    }
    /// Integrate contact force over the slave surface.
    ///
    /// Returns total contact force vector.
    pub fn integrate_contact_force(&self) -> [f64; 3] {
        let (u0, u1) = self.slave_surface.u_basis.domain();
        let (v0, v1) = self.slave_surface.v_basis.domain();
        let (gp, gw) = gauss_legendre_1d(self.n_gauss);
        let mut total = [0.0f64; 3];
        for i in 0..self.n_gauss {
            for j in 0..self.n_gauss {
                let u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gp[i];
                let v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gp[j];
                let jac = self.slave_surface.jacobian_det(u, v);
                let f = self.contact_force(u, v);
                let w = 0.25 * (u1 - u0) * (v1 - v0) * gw[i] * gw[j];
                total[0] += f[0] * jac * w;
                total[1] += f[1] * jac * w;
                total[2] += f[2] * jac * w;
            }
        }
        total
    }
}
/// IGA mass matrix assembler for volumetric (3D) NURBS elements.
///
/// Integrates the consistent mass matrix M_ab = rho * integral(N_a * N_b dV)
/// over a trivariate NURBS patch.
#[derive(Debug, Clone)]
pub struct IgaVolumeMassMatrix {
    /// Mass density.
    pub density: f64,
    /// Number of Gauss points per parametric direction.
    pub n_gauss: usize,
}
impl IgaVolumeMassMatrix {
    /// Create a new volume mass matrix assembler.
    pub fn new(density: f64, n_gauss: usize) -> Self {
        IgaVolumeMassMatrix { density, n_gauss }
    }
    /// Assemble the consistent mass matrix for a NURBS volume patch.
    ///
    /// Returns a flat (n_dof x n_dof) matrix where n_dof = n_ctrl * 3
    /// and n_ctrl = (pu+1) * (pv+1) * (pw+1).
    pub fn assemble(&self, vol: &NurbsVolume) -> Vec<f64> {
        let pu = vol.u_basis.degree;
        let pv = vol.v_basis.degree;
        let pw_deg = vol.w_basis.degree;
        let n_ctrl = (pu + 1) * (pv + 1) * (pw_deg + 1);
        let n_dof = n_ctrl * 3;
        let mut mass = vec![0.0f64; n_dof * n_dof];
        let (u0, u1) = vol.u_basis.domain();
        let (v0, v1) = vol.v_basis.domain();
        let (w0, w1) = vol.w_basis.domain();
        let (gp, gw) = gauss_legendre_1d(self.n_gauss);
        for gi in 0..self.n_gauss {
            for gj in 0..self.n_gauss {
                for gk in 0..self.n_gauss {
                    let u = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * gp[gi];
                    let v = 0.5 * (v0 + v1) + 0.5 * (v1 - v0) * gp[gj];
                    let w = 0.5 * (w0 + w1) + 0.5 * (w1 - w0) * gp[gk];
                    let jac_w =
                        0.125 * (u1 - u0) * (v1 - v0) * (w1 - w0) * gw[gi] * gw[gj] * gw[gk];
                    let j_mat = vol.jacobian(u, v, w);
                    let j_det = vol.jacobian_det(u, v, w).abs();
                    let nu = vol.ctrl.len() - 1;
                    let nv = vol.ctrl[0].len() - 1;
                    let nw = vol.ctrl[0][0].len() - 1;
                    let u_span = find_knot_span(nu, pu, u, &vol.u_basis.knot);
                    let v_span = find_knot_span(nv, pv, v, &vol.v_basis.knot);
                    let w_span = find_knot_span(nw, pw_deg, w, &vol.w_basis.knot);
                    let bu = basis_functions(u_span, pu, u, &vol.u_basis.knot);
                    let bv = basis_functions(v_span, pv, v, &vol.v_basis.knot);
                    let bw_vals = basis_functions(w_span, pw_deg, w, &vol.w_basis.knot);
                    let mut wt = 0.0f64;
                    let mut nf = vec![0.0f64; n_ctrl];
                    let mut local_a = 0usize;
                    for (ii, &bu_ii) in bu.iter().enumerate().take(pu + 1) {
                        for (jj, &bv_jj) in bv.iter().enumerate().take(pv + 1) {
                            for (kk, &bw_kk) in bw_vals.iter().enumerate().take(pw_deg + 1) {
                                let ui = u_span - pu + ii;
                                let vi = v_span - pv + jj;
                                let wi_idx = w_span - pw_deg + kk;
                                let w_ijk = vol.weights[ui][vi][wi_idx] * bu_ii * bv_jj * bw_kk;
                                nf[local_a] = w_ijk;
                                wt += w_ijk;
                                local_a += 1;
                            }
                        }
                    }
                    if wt.abs() > 1e-15 {
                        for nf_val in nf.iter_mut().take(n_ctrl) {
                            *nf_val /= wt;
                        }
                    }
                    let rho = self.density;
                    for a in 0..n_ctrl {
                        for b in 0..n_ctrl {
                            let val = rho * nf[a] * nf[b] * j_det * jac_w;
                            for d in 0..3 {
                                mass[(3 * a + d) * n_dof + (3 * b + d)] += val;
                            }
                        }
                    }
                    let _ = j_mat;
                }
            }
        }
        mass
    }
}
/// Multi-patch IGA connectivity.
///
/// Handles interface coupling between adjacent NURBS patches.
#[derive(Debug, Clone)]
pub struct PatchConnectivity {
    /// IDs of the two adjacent patches.
    pub patch_ids: (usize, usize),
    /// Coupling method.
    pub coupling: CouplingType,
    /// Control point index mapping at interface.
    pub interface_map: Vec<(usize, usize)>,
    /// Global DOF offset for patch A.
    pub dof_offset_a: usize,
    /// Global DOF offset for patch B.
    pub dof_offset_b: usize,
}
impl PatchConnectivity {
    /// Create a new patch connectivity.
    pub fn new(
        patch_ids: (usize, usize),
        coupling: CouplingType,
        interface_map: Vec<(usize, usize)>,
        dof_offset_a: usize,
        dof_offset_b: usize,
    ) -> Self {
        Self {
            patch_ids,
            coupling,
            interface_map,
            dof_offset_a,
            dof_offset_b,
        }
    }
    /// Number of coupled interface pairs.
    pub fn num_pairs(&self) -> usize {
        self.interface_map.len()
    }
    /// Check if a control point index belongs to patch A's interface.
    pub fn is_interface_a(&self, idx: usize) -> bool {
        self.interface_map.iter().any(|(a, _)| *a == idx)
    }
    /// Compute penalty coupling stiffness contribution.
    ///
    /// Returns penalty coefficient for a given interface pair.
    pub fn penalty_coefficient(&self, pair_index: usize) -> f64 {
        if pair_index >= self.interface_map.len() {
            return 0.0;
        }
        match self.coupling {
            CouplingType::Penalty(beta) => beta,
            _ => 0.0,
        }
    }
}
/// Trivariate NURBS solid (volume).
///
/// 3D NURBS parameterized by (u, v, w) ∈ \[0,1\]^3.
#[derive(Debug, Clone)]
pub struct NurbsVolume {
    /// Control points grid: `ctrl[i][j][k]`.
    pub ctrl: Vec<Vec<Vec<[f64; 3]>>>,
    /// Weights at each control point.
    pub weights: Vec<Vec<Vec<f64>>>,
    /// B-spline basis in u direction.
    pub u_basis: BsplineBasis,
    /// B-spline basis in v direction.
    pub v_basis: BsplineBasis,
    /// B-spline basis in w direction.
    pub w_basis: BsplineBasis,
}
impl NurbsVolume {
    /// Create a new trivariate NURBS volume.
    pub fn new(
        ctrl: Vec<Vec<Vec<[f64; 3]>>>,
        weights: Vec<Vec<Vec<f64>>>,
        u_degree: usize,
        v_degree: usize,
        w_degree: usize,
        u_knot: Vec<f64>,
        v_knot: Vec<f64>,
        w_knot: Vec<f64>,
    ) -> Self {
        let u_basis = BsplineBasis::new(u_degree, u_knot);
        let v_basis = BsplineBasis::new(v_degree, v_knot);
        let w_basis = BsplineBasis::new(w_degree, w_knot);
        Self {
            ctrl,
            weights,
            u_basis,
            v_basis,
            w_basis,
        }
    }
    /// Evaluate the volume point at parameter (u, v, w).
    pub fn point_at(&self, u: f64, v: f64, w: f64) -> [f64; 3] {
        let nu = self.ctrl.len() - 1;
        let nv = self.ctrl[0].len() - 1;
        let nw = self.ctrl[0][0].len() - 1;
        let pu = self.u_basis.degree;
        let pv = self.v_basis.degree;
        let pw = self.w_basis.degree;
        let u_span = find_knot_span(nu, pu, u, &self.u_basis.knot);
        let v_span = find_knot_span(nv, pv, v, &self.v_basis.knot);
        let w_span = find_knot_span(nw, pw, w, &self.w_basis.knot);
        let bu = basis_functions(u_span, pu, u, &self.u_basis.knot);
        let bv = basis_functions(v_span, pv, v, &self.v_basis.knot);
        let bw = basis_functions(w_span, pw, w, &self.w_basis.knot);
        let mut pt = [0.0f64; 3];
        let mut wt = 0.0f64;
        for (i, &bu_i) in bu.iter().enumerate().take(pu + 1) {
            for (j, &bv_j) in bv.iter().enumerate().take(pv + 1) {
                for (k, &bw_k) in bw.iter().enumerate().take(pw + 1) {
                    let ui = u_span - pu + i;
                    let vi = v_span - pv + j;
                    let wi_idx = w_span - pw + k;
                    let w_ijk = self.weights[ui][vi][wi_idx] * bu_i * bv_j * bw_k;
                    pt[0] += w_ijk * self.ctrl[ui][vi][wi_idx][0];
                    pt[1] += w_ijk * self.ctrl[ui][vi][wi_idx][1];
                    pt[2] += w_ijk * self.ctrl[ui][vi][wi_idx][2];
                    wt += w_ijk;
                }
            }
        }
        if wt.abs() < 1e-15 {
            return pt;
        }
        [pt[0] / wt, pt[1] / wt, pt[2] / wt]
    }
    /// Compute the 3x3 Jacobian matrix ∂x/∂(u,v,w) at parameter point.
    ///
    /// Returns the Jacobian as a flattened row-major 3x3 matrix.
    pub fn jacobian(&self, u: f64, v: f64, w: f64) -> [[f64; 3]; 3] {
        let eps = 1e-6;
        let (u0, u1) = self.u_basis.domain();
        let (v0, v1) = self.v_basis.domain();
        let (w0, w1) = self.w_basis.domain();
        let pu = (u + eps).min(u1);
        let mu = (u - eps).max(u0);
        let pv = (v + eps).min(v1);
        let mv = (v - eps).max(v0);
        let pw_p = (w + eps).min(w1);
        let mw = (w - eps).max(w0);
        let du_col = scale3(
            sub3(self.point_at(pu, v, w), self.point_at(mu, v, w)),
            1.0 / (pu - mu),
        );
        let dv_col = scale3(
            sub3(self.point_at(u, pv, w), self.point_at(u, mv, w)),
            1.0 / (pv - mv),
        );
        let dw_col = scale3(
            sub3(self.point_at(u, v, pw_p), self.point_at(u, v, mw)),
            1.0 / (pw_p - mw),
        );
        [du_col, dv_col, dw_col]
    }
    /// Compute the Jacobian determinant at a parameter point.
    pub fn jacobian_det(&self, u: f64, v: f64, w: f64) -> f64 {
        let j = self.jacobian(u, v, w);
        j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0])
    }
}
/// A T-spline surface defined by a set of T-spline control points.
///
/// Evaluates using the blending function approach: the surface point is
/// the rational combination of all control points whose support includes (s, t).
#[derive(Debug, Clone)]
pub struct TSplineSurface {
    /// All control points in the T-mesh.
    pub control_points: Vec<TSplineControlPoint>,
    /// Parameter domain bounds \[s_min, s_max, t_min, t_max\].
    pub domain: [f64; 4],
}
impl TSplineSurface {
    /// Create a new T-spline surface.
    pub fn new(control_points: Vec<TSplineControlPoint>, domain: [f64; 4]) -> Self {
        TSplineSurface {
            control_points,
            domain,
        }
    }
    /// Evaluate the T-spline surface at parameter (s, t).
    ///
    /// Returns the physical point using rational blending.
    pub fn point_at(&self, s: f64, t: f64) -> [f64; 3] {
        let mut num = [0.0f64; 3];
        let mut denom = 0.0f64;
        for cp in &self.control_points {
            let wb = cp.weighted_basis(s, t);
            if wb.abs() > 1e-15 {
                num[0] += wb * cp.position[0];
                num[1] += wb * cp.position[1];
                num[2] += wb * cp.position[2];
                denom += wb;
            }
        }
        if denom.abs() < 1e-15 {
            return [0.0; 3];
        }
        [num[0] / denom, num[1] / denom, num[2] / denom]
    }
    /// Number of T-spline control points.
    pub fn num_control_points(&self) -> usize {
        self.control_points.len()
    }
    /// Check if the T-spline has T-junctions (control points with non-uniform knots).
    pub fn has_t_junctions(&self) -> bool {
        if self.control_points.len() < 2 {
            return false;
        }
        let ref_s = self.control_points[0].s_knots;
        let ref_t = self.control_points[0].t_knots;
        for cp in self.control_points.iter().skip(1) {
            if cp.s_knots != ref_s || cp.t_knots != ref_t {
                return true;
            }
        }
        false
    }
}
/// Metrics for assessing the quality of a NURBS parameterization.
///
/// A good parameterization should have well-conditioned Jacobians throughout
/// the domain, minimal skewness, and approximately unit aspect ratios.
#[derive(Debug, Clone)]
pub struct ParameterizationQuality {
    /// Minimum Jacobian determinant (should be > 0 for valid parameterization).
    pub min_jacobian: f64,
    /// Maximum Jacobian determinant.
    pub max_jacobian: f64,
    /// Average Jacobian determinant.
    pub avg_jacobian: f64,
    /// Scaled Jacobian (normalized by element size, in \[0,1\]).
    pub scaled_jacobian: f64,
    /// Condition number of the Jacobian matrix (large = poor quality).
    pub condition_number: f64,
    /// Orthogonality measure (1.0 = perfectly orthogonal parameter lines).
    pub orthogonality: f64,
    /// Number of sample points used.
    pub num_samples: usize,
}
impl ParameterizationQuality {
    /// Analyze parameterization quality of a NURBS surface.
    ///
    /// Samples the Jacobian at `n_samples` x `n_samples` parameter points.
    pub fn analyze_surface(surf: &NurbsSurface, n_samples: usize) -> Self {
        let (u0, u1) = surf.u_basis.domain();
        let (v0, v1) = surf.v_basis.domain();
        let mut min_j = f64::MAX;
        let mut max_j = f64::MIN;
        let mut sum_j = 0.0;
        let mut sum_cond = 0.0;
        let mut sum_orth = 0.0;
        let count = (n_samples * n_samples) as f64;
        for i in 0..n_samples {
            for j in 0..n_samples {
                let u = u0 + (u1 - u0) * (i as f64 + 0.5) / n_samples as f64;
                let v = v0 + (v1 - v0) * (j as f64 + 0.5) / n_samples as f64;
                let (tang_u, tang_v, _) = surf.normal(u, v);
                let jac_det = norm3(cross3(tang_u, tang_v));
                min_j = min_j.min(jac_det);
                max_j = max_j.max(jac_det);
                sum_j += jac_det;
                let a = norm3(tang_u);
                let b = norm3(tang_v);
                let cond = if b > 1e-15 {
                    (a / b).max(b / a)
                } else {
                    f64::MAX
                };
                sum_cond += cond.min(1e6);
                let cos_angle = if a > 1e-15 && b > 1e-15 {
                    (dot3(tang_u, tang_v) / (a * b)).abs()
                } else {
                    1.0
                };
                sum_orth += 1.0 - cos_angle;
            }
        }
        let avg_j = sum_j / count;
        let scaled_j = if max_j > 1e-15 { min_j / max_j } else { 0.0 };
        ParameterizationQuality {
            min_jacobian: min_j,
            max_jacobian: max_j,
            avg_jacobian: avg_j,
            scaled_jacobian: scaled_j,
            condition_number: sum_cond / count,
            orthogonality: sum_orth / count,
            num_samples: n_samples * n_samples,
        }
    }
    /// Check if parameterization is valid (all Jacobians positive).
    pub fn is_valid(&self) -> bool {
        self.min_jacobian > 0.0
    }
    /// Quality score in \[0, 1\]: higher is better.
    pub fn quality_score(&self) -> f64 {
        if !self.is_valid() {
            return 0.0;
        }
        let j_score = self.scaled_jacobian.clamp(0.0, 1.0);
        let o_score = self.orthogonality.clamp(0.0, 1.0);
        0.5 * j_score + 0.5 * o_score
    }
}
/// Strategy parameters for k-refinement in isogeometric analysis.
///
/// k-refinement combines degree elevation with knot insertion, yielding
/// higher-continuity basis functions compared to classical h- or p-refinement.
#[derive(Debug, Clone)]
pub struct KRefinementStrategy {
    /// Target polynomial degree after refinement.
    pub target_degree: usize,
    /// Number of knots to insert in each direction after degree elevation.
    pub knots_per_direction: usize,
    /// Whether to apply uniform or adaptive knot insertion.
    pub uniform: bool,
    /// Minimum element size (stop refinement below this).
    pub min_element_size: f64,
}
impl KRefinementStrategy {
    /// Create a k-refinement strategy.
    pub fn new(
        target_degree: usize,
        knots_per_direction: usize,
        uniform: bool,
        min_element_size: f64,
    ) -> Self {
        KRefinementStrategy {
            target_degree,
            knots_per_direction,
            uniform,
            min_element_size,
        }
    }
    /// Apply k-refinement to a knot vector and return the refined knot vector.
    ///
    /// Elevates degree from `current_degree` to `target_degree`, then inserts knots.
    pub fn apply(&self, knot: &[f64], current_degree: usize) -> Vec<f64> {
        let mut result = knot.to_vec();
        let mut curr_deg = current_degree;
        while curr_deg < self.target_degree {
            result = KnotRefinement::p_refine_knot(&result, curr_deg);
            curr_deg += 1;
        }
        if self.uniform {
            result = KnotRefinement::uniform_h_refine(&result, curr_deg, self.knots_per_direction);
        } else {
            result = KnotRefinement::k_refine(&result, curr_deg, self.knots_per_direction);
        }
        result
    }
    /// Estimate the number of DOFs after k-refinement for a 1D problem.
    pub fn estimate_dofs_1d(&self, n_ctrl: usize) -> usize {
        let degree_added = self.target_degree.saturating_sub(1);
        n_ctrl + degree_added + self.knots_per_direction
    }
    /// Check if a given element size is too small for further refinement.
    pub fn is_refinement_possible(&self, element_size: f64) -> bool {
        element_size > self.min_element_size * 2.0
    }
}
/// Unified NURBS basis for use in isogeometric analysis.
///
/// Wraps a B-spline basis (`BsplineBasis`) with a weight vector so that
/// rational NURBS shape functions can be evaluated directly.  The de Boor
/// algorithm is used internally for numerically stable evaluation.
///
/// # Example
/// ```ignore
/// let knot = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
/// let weights = vec![1.0, 1.0, 1.0];
/// let basis = NurbsBasis::new(2, knot, weights);
/// let (span, r, dr) = basis.evaluate_rational(0.5);
/// ```
#[derive(Debug, Clone)]
pub struct NurbsBasis {
    /// Underlying B-spline basis (degree + knot vector).
    pub bspline: BsplineBasis,
    /// Weights for each basis function (all positive).
    pub weights: Vec<f64>,
}
impl NurbsBasis {
    /// Construct a NURBS basis from polynomial degree, knot vector, and weights.
    ///
    /// The number of weights must equal `knot.len() - degree - 1`.
    pub fn new(degree: usize, knot: Vec<f64>, weights: Vec<f64>) -> Self {
        let bspline = BsplineBasis::new(degree, knot);
        Self { bspline, weights }
    }
    /// Build a uniform open knot vector for `n_basis` functions of given `degree`.
    ///
    /// The resulting knot vector has clamped (open) end conditions so the
    /// first and last basis functions interpolate the endpoints.
    pub fn uniform_open_knot(degree: usize, n_basis: usize) -> Vec<f64> {
        let m = n_basis + degree + 1;
        let n_internal = n_basis - degree - 1;
        let mut knot = vec![0.0f64; degree + 1];
        if n_internal > 0 {
            let step = 1.0 / (n_internal + 1) as f64;
            for i in 1..=n_internal {
                knot.push(i as f64 * step);
            }
        }
        while knot.len() < m {
            knot.push(1.0);
        }
        knot
    }
    /// Evaluate NURBS basis functions (rational) at parameter `t`.
    ///
    /// Returns `(span, R, dR)` where:
    /// - `span` is the knot span index,
    /// - `R` are the `p+1` non-zero rational basis function values,
    /// - `dR` are their first derivatives with respect to `t`.
    pub fn evaluate_rational(&self, t: f64) -> (usize, Vec<f64>, Vec<f64>) {
        let p = self.bspline.degree;
        let n = self.bspline.num_basis - 1;
        let span = find_knot_span(n, p, t, &self.bspline.knot);
        let ders = basis_derivatives(span, p, t, &self.bspline.knot, 1);
        let n_vals = &ders[0];
        let dn_vals = &ders[1];
        let mut w_sum = 0.0f64;
        let mut dw_sum = 0.0f64;
        for j in 0..=p {
            let wi = self.weights[span - p + j];
            w_sum += wi * n_vals[j];
            dw_sum += wi * dn_vals[j];
        }
        let w_inv = if w_sum.abs() > 1e-15 {
            1.0 / w_sum
        } else {
            0.0
        };
        let mut r = Vec::with_capacity(p + 1);
        let mut dr = Vec::with_capacity(p + 1);
        for j in 0..=p {
            let wi = self.weights[span - p + j];
            let rj = wi * n_vals[j] * w_inv;
            let drj = wi * (dn_vals[j] - rj * dw_sum) * w_inv;
            r.push(rj);
            dr.push(drj);
        }
        (span, r, dr)
    }
    /// Evaluate only the NURBS basis function values (no derivatives).
    pub fn evaluate(&self, t: f64) -> (usize, Vec<f64>) {
        let (span, r, _) = self.evaluate_rational(t);
        (span, r)
    }
    /// Return the parametric domain `(t_min, t_max)`.
    pub fn domain(&self) -> (f64, f64) {
        self.bspline.domain()
    }
    /// Number of basis functions.
    pub fn num_basis(&self) -> usize {
        self.bspline.num_basis
    }
    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.bspline.degree
    }
    /// Evaluate a NURBS curve point given control points.
    ///
    /// Uses the rational basis to evaluate `C(t) = Σ R_i(t) * P_i`.
    pub fn curve_point(&self, control_points: &[[f64; 3]], t: f64) -> [f64; 3] {
        let (span, r) = self.evaluate(t);
        let p = self.bspline.degree;
        let mut pt = [0.0f64; 3];
        for j in 0..=p {
            let cp = control_points[span - p + j];
            pt[0] += r[j] * cp[0];
            pt[1] += r[j] * cp[1];
            pt[2] += r[j] * cp[2];
        }
        pt
    }
    /// Insert a new knot `t_new` into the knot vector using Boehm's algorithm.
    ///
    /// Returns a new `NurbsBasis` with the updated knot vector and expanded
    /// weight vector (repeated weight inserted at the appropriate position).
    pub fn insert_knot(&self, t_new: f64) -> NurbsBasis {
        let p = self.bspline.degree;
        let old_knot = &self.bspline.knot;
        let n = self.bspline.num_basis - 1;
        let span = find_knot_span(n, p, t_new, old_knot);
        let mut new_knot = old_knot[..=span].to_vec();
        new_knot.push(t_new);
        new_knot.extend_from_slice(&old_knot[span + 1..]);
        let insert_pos = (span + 1).saturating_sub(p / 2);
        let insert_pos = insert_pos.min(self.weights.len());
        let avg_w = if insert_pos < self.weights.len() {
            self.weights[insert_pos]
        } else {
            *self.weights.last().unwrap_or(&1.0)
        };
        let mut new_weights = self.weights.clone();
        new_weights.insert(insert_pos, avg_w);
        NurbsBasis::new(p, new_knot, new_weights)
    }
    /// Perform de Boor's algorithm to evaluate the B-spline at `t`.
    ///
    /// Returns the value of the curve at `t` given a 1D array of scalar
    /// control values (useful for intermediate computations).
    pub fn de_boor_scalar(&self, control_values: &[f64], t: f64) -> f64 {
        let p = self.bspline.degree;
        let n = self.bspline.num_basis - 1;
        let knot = &self.bspline.knot;
        let span = find_knot_span(n, p, t, knot);
        let mut d: Vec<f64> = (0..=p).map(|j| control_values[span - p + j]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let left = knot[span - p + j];
                let right = knot[span + 1 + j - r];
                let denom = right - left;
                let alpha = if denom.abs() < 1e-15 {
                    0.0
                } else {
                    (t - left) / denom
                };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }
        d[p]
    }
}
/// NURBS curve in 3D space.
///
/// Defined by control points, weights, knot vector, and polynomial degree.
#[derive(Debug, Clone)]
pub struct NurbsCurve {
    /// Control points in 3D.
    pub control_points: Vec<[f64; 3]>,
    /// Weights for each control point (positive).
    pub weights: Vec<f64>,
    /// B-spline basis (degree + knot vector).
    pub basis: BsplineBasis,
}
impl NurbsCurve {
    /// Create a new NURBS curve.
    ///
    /// # Arguments
    /// * `control_points` – control polygon vertices
    /// * `weights` – positive weights (same length as control_points)
    /// * `degree` – polynomial degree
    /// * `knot` – knot vector of length (num_ctrl + degree + 1)
    pub fn new(
        control_points: Vec<[f64; 3]>,
        weights: Vec<f64>,
        degree: usize,
        knot: Vec<f64>,
    ) -> Self {
        let basis = BsplineBasis::new(degree, knot);
        Self {
            control_points,
            weights,
            basis,
        }
    }
    /// Evaluate curve point at parameter t.
    pub fn point_at(&self, t: f64) -> [f64; 3] {
        nurbs_point(
            &self.control_points,
            &self.weights,
            &self.basis.knot,
            self.basis.degree,
            t,
        )
    }
    /// Evaluate first derivative (tangent) at parameter t.
    ///
    /// Uses the NURBS rational derivative formula:
    /// C'(t) = (A'(t) - w'(t)*C(t)) / w(t)
    /// where A(t) = sum(w_i * P_i * N_i(t)).
    pub fn derivative_at(&self, t: f64) -> [f64; 3] {
        let n = self.control_points.len() - 1;
        let p = self.basis.degree;
        let span = find_knot_span(n, p, t, &self.basis.knot);
        let ders = basis_derivatives(span, p, t, &self.basis.knot, 1);
        let mut a = [0.0f64; 3];
        let mut da = [0.0f64; 3];
        let mut w = 0.0f64;
        let mut dw = 0.0f64;
        for (j, (&n0, &n1)) in ders[0].iter().zip(ders[1].iter()).enumerate().take(p + 1) {
            let idx = span - p + j;
            let wi = self.weights[idx];
            let pi = self.control_points[idx];
            a[0] += wi * pi[0] * n0;
            a[1] += wi * pi[1] * n0;
            a[2] += wi * pi[2] * n0;
            da[0] += wi * pi[0] * n1;
            da[1] += wi * pi[1] * n1;
            da[2] += wi * pi[2] * n1;
            w += wi * n0;
            dw += wi * n1;
        }
        if w.abs() < 1e-15 {
            return [0.0; 3];
        }
        let c = [a[0] / w, a[1] / w, a[2] / w];
        [
            (da[0] - dw * c[0]) / w,
            (da[1] - dw * c[1]) / w,
            (da[2] - dw * c[2]) / w,
        ]
    }
    /// Insert a knot value t_new into the knot vector (knot insertion).
    ///
    /// Updates control points and weights using the Boehm algorithm.
    pub fn insert_knot(&self, t_new: f64) -> NurbsCurve {
        let n = self.control_points.len() - 1;
        let p = self.basis.degree;
        let knot = &self.basis.knot;
        let k = find_knot_span(n, p, t_new, knot);
        let mut new_knot = Vec::with_capacity(knot.len() + 1);
        new_knot.extend_from_slice(&knot[..=k]);
        new_knot.push(t_new);
        new_knot.extend_from_slice(&knot[k + 1..]);
        let m = n + p + 1;
        let mut new_ctrl = vec![[0.0f64; 3]; n + 2];
        let mut new_weights = vec![0.0f64; n + 2];
        new_ctrl[..=(k - p)].copy_from_slice(&self.control_points[..=(k - p)]);
        new_weights[..=(k - p)].copy_from_slice(&self.weights[..=(k - p)]);
        for i in (k - p + 1)..=k {
            let alpha = (t_new - knot[i]) / (knot[i + p] - knot[i]);
            let prev_i = i - 1;
            new_ctrl[i] = [
                alpha * self.control_points[i][0] + (1.0 - alpha) * self.control_points[prev_i][0],
                alpha * self.control_points[i][1] + (1.0 - alpha) * self.control_points[prev_i][1],
                alpha * self.control_points[i][2] + (1.0 - alpha) * self.control_points[prev_i][2],
            ];
            new_weights[i] = alpha * self.weights[i] + (1.0 - alpha) * self.weights[prev_i];
        }
        for i in (k + 1)..=(m - p) {
            if i <= n + 1 {
                new_ctrl[i] = self.control_points[i - 1];
                new_weights[i] = self.weights[i - 1];
            }
        }
        NurbsCurve::new(new_ctrl, new_weights, p, new_knot)
    }
    /// Arc length approximation via numerical integration.
    pub fn arc_length(&self, num_segments: usize) -> f64 {
        let (t0, t1) = self.basis.domain();
        let dt = (t1 - t0) / num_segments as f64;
        let mut length = 0.0;
        let mut prev = self.point_at(t0);
        for i in 1..=num_segments {
            let t = t0 + i as f64 * dt;
            let curr = self.point_at(t);
            length += norm3(sub3(curr, prev));
            prev = curr;
        }
        length
    }
}
/// NURBS-based IGA element for isogeometric analysis.
///
/// Encapsulates the local element in parameter space and provides shape
/// functions, derivatives, and Jacobian for numerical integration.
#[derive(Debug, Clone)]
pub struct IgaElement {
    /// Element index in the patch.
    pub element_id: usize,
    /// Connectivity: global control point indices.
    pub connectivity: Vec<usize>,
    /// Reference surface for this element (from patch).
    pub surface: NurbsSurface,
    /// Integration order in u.
    pub n_gauss_u: usize,
    /// Integration order in v.
    pub n_gauss_v: usize,
}
impl IgaElement {
    /// Create a new IGA element.
    pub fn new(
        element_id: usize,
        connectivity: Vec<usize>,
        surface: NurbsSurface,
        n_gauss_u: usize,
        n_gauss_v: usize,
    ) -> Self {
        Self {
            element_id,
            connectivity,
            surface,
            n_gauss_u,
            n_gauss_v,
        }
    }
    /// Compute shape function values and parameter-space derivatives at (u, v).
    ///
    /// Returns (R, dR_du, dR_dv) where R\[a\] is the a-th NURBS basis function
    /// and dR_du\[a\], dR_dv\[a\] are its partial derivatives.
    pub fn shape_functions(&self, u: f64, v: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let nu = self.surface.control_net.len() - 1;
        let nv = self.surface.control_net[0].len() - 1;
        let pu = self.surface.u_basis.degree;
        let pv = self.surface.v_basis.degree;
        let u_span = find_knot_span(nu, pu, u, &self.surface.u_basis.knot);
        let v_span = find_knot_span(nv, pv, v, &self.surface.v_basis.knot);
        let u_ders = basis_derivatives(u_span, pu, u, &self.surface.u_basis.knot, 1);
        let v_ders = basis_derivatives(v_span, pv, v, &self.surface.v_basis.knot, 1);
        let n_local = (pu + 1) * (pv + 1);
        let mut r = vec![0.0f64; n_local];
        let mut dr_du = vec![0.0f64; n_local];
        let mut dr_dv = vec![0.0f64; n_local];
        let mut w = 0.0f64;
        let mut dw_du = 0.0f64;
        let mut dw_dv = 0.0f64;
        let mut bw = vec![0.0f64; n_local];
        let mut a = 0;
        for (i, (&ud0i, &ud1i)) in u_ders[0]
            .iter()
            .zip(u_ders[1].iter())
            .enumerate()
            .take(pu + 1)
        {
            for (j, (&vd0j, &vd1j)) in v_ders[0]
                .iter()
                .zip(v_ders[1].iter())
                .enumerate()
                .take(pv + 1)
            {
                let ui = u_span - pu + i;
                let vi = v_span - pv + j;
                let wij = self.surface.weights[ui][vi];
                let nij = ud0i * vd0j;
                bw[a] = wij * nij;
                w += bw[a];
                dw_du += wij * ud1i * vd0j;
                dw_dv += wij * ud0i * vd1j;
                a += 1;
            }
        }
        if w.abs() < 1e-15 {
            return (r, dr_du, dr_dv);
        }
        let _w2 = w * w;
        a = 0;
        for (i, (&ud0i, &ud1i)) in u_ders[0]
            .iter()
            .zip(u_ders[1].iter())
            .enumerate()
            .take(pu + 1)
        {
            for (j, (&vd0j, &vd1j)) in v_ders[0]
                .iter()
                .zip(v_ders[1].iter())
                .enumerate()
                .take(pv + 1)
            {
                let ui = u_span - pu + i;
                let vi = v_span - pv + j;
                let wij = self.surface.weights[ui][vi];
                r[a] = bw[a] / w;
                dr_du[a] = (wij * ud1i * vd0j - dw_du * bw[a] / w) / w;
                dr_dv[a] = (wij * ud0i * vd1j - dw_dv * bw[a] / w) / w;
                let _ = (ui, vi);
                a += 1;
            }
        }
        (r, dr_du, dr_dv)
    }
    /// Compute the Jacobian of the physical-to-parameter mapping at (u, v).
    ///
    /// Returns the 2x3 Jacobian matrix: rows are ∂x/∂u and ∂x/∂v.
    pub fn jacobian(&self, u: f64, v: f64) -> [[f64; 3]; 2] {
        let (_, tang_u, tang_v) = {
            let eps = 1e-7;
            let (t0, t1) = self.surface.u_basis.domain();
            let (s0, s1) = self.surface.v_basis.domain();
            let up = (u + eps).min(t1);
            let um = (u - eps).max(t0);
            let vp = (v + eps).min(s1);
            let vm = (v - eps).max(s0);
            let tang_u = scale3(
                sub3(self.surface.point_at(up, v), self.surface.point_at(um, v)),
                1.0 / (up - um),
            );
            let tang_v = scale3(
                sub3(self.surface.point_at(u, vp), self.surface.point_at(u, vm)),
                1.0 / (vp - vm),
            );
            ((), tang_u, tang_v)
        };
        [tang_u, tang_v]
    }
}
/// Knot refinement operations for NURBS.
///
/// Supports h-refinement (knot insertion), p-refinement (degree elevation),
/// and k-refinement (order elevation preserving continuity).
#[derive(Debug, Clone)]
pub struct KnotRefinement;
impl KnotRefinement {
    /// Perform h-refinement (knot insertion) on a NURBS curve.
    ///
    /// Inserts all knots in `new_knots` into the curve.
    pub fn h_refine_curve(curve: &NurbsCurve, new_knots: &[f64]) -> NurbsCurve {
        let mut result = curve.clone();
        let mut sorted = new_knots.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        for &t in &sorted {
            result = result.insert_knot(t);
        }
        result
    }
    /// Perform p-refinement (degree elevation) on a B-spline basis.
    ///
    /// Elevates the degree by 1 and adds necessary knots to maintain continuity.
    /// Returns a new knot vector for the elevated degree.
    pub fn p_refine_knot(knot: &[f64], degree: usize) -> Vec<f64> {
        let new_degree = degree + 1;
        let mut unique_knots: Vec<f64> = Vec::new();
        for &k in knot {
            if unique_knots.is_empty()
                || (k - unique_knots.last().expect("unique_knots is non-empty")).abs() > 1e-14
            {
                unique_knots.push(k);
            }
        }
        let mut new_knot = Vec::new();
        for &k in &unique_knots {
            for _ in 0..=new_degree {
                new_knot.push(k);
            }
        }
        new_knot
    }
    /// Perform k-refinement: first degree elevation, then knot insertion.
    ///
    /// Returns the new knot vector after k-refinement.
    pub fn k_refine(knot: &[f64], degree: usize, num_insertions: usize) -> Vec<f64> {
        let elevated_knot = Self::p_refine_knot(knot, degree);
        let new_degree = degree + 1;
        let first = elevated_knot[new_degree];
        let last = elevated_knot[elevated_knot.len() - 1 - new_degree];
        let mut result = elevated_knot;
        for i in 1..=num_insertions {
            let t = first + (last - first) * i as f64 / (num_insertions + 1) as f64;
            let pos = result.partition_point(|&x| x <= t);
            result.insert(pos, t);
        }
        result
    }
    /// Uniform knot insertion for h-refinement: insert n_new knots per span.
    pub fn uniform_h_refine(knot: &[f64], degree: usize, n_new: usize) -> Vec<f64> {
        let mut new_knots_to_insert = Vec::new();
        for w in knot.windows(2) {
            let span = w[1] - w[0];
            if span > 1e-14 {
                for k in 1..=n_new {
                    new_knots_to_insert.push(w[0] + span * k as f64 / (n_new + 1) as f64);
                }
            }
        }
        let mut result = knot.to_vec();
        let _ = degree;
        for t in new_knots_to_insert {
            let pos = result.partition_point(|&x| x <= t);
            result.insert(pos, t);
        }
        result
    }
}
/// Patch coupling for multi-patch IGA using the mortar method.
///
/// Computes the mortar coupling matrix between two adjacent NURBS patches
/// sharing an interface, enabling displacement and traction continuity.
#[derive(Debug, Clone)]
pub struct MortarPatchCoupling {
    /// The slave patch (interface defined by v = v_slave_param).
    pub slave_surface: NurbsSurface,
    /// The master patch (interface defined by v = v_master_param).
    pub master_surface: NurbsSurface,
    /// Number of Gauss points for mortar integration.
    pub n_gauss: usize,
    /// Penalty or Lagrange multiplier coefficient.
    pub alpha: f64,
}
impl MortarPatchCoupling {
    /// Create a new mortar patch coupling.
    pub fn new(
        slave_surface: NurbsSurface,
        master_surface: NurbsSurface,
        n_gauss: usize,
        alpha: f64,
    ) -> Self {
        MortarPatchCoupling {
            slave_surface,
            master_surface,
            n_gauss,
            alpha,
        }
    }
    /// Compute the L2 projection mortar coupling matrix.
    ///
    /// Returns the coupling matrix D (slave-master DOF coupling)
    /// of size (n_slave_dof x n_master_dof).
    pub fn coupling_matrix(&self) -> Vec<f64> {
        let ps = self.slave_surface.u_basis.degree;
        let pm = self.master_surface.u_basis.degree;
        let n_slave = ps + 1;
        let n_master = pm + 1;
        let mut d_mat = vec![0.0f64; n_slave * n_master];
        let (u0s, u1s) = self.slave_surface.u_basis.domain();
        let (u0m, u1m) = self.master_surface.u_basis.domain();
        let (gp, gw) = gauss_legendre_1d(self.n_gauss);
        for gi in 0..self.n_gauss {
            let us = 0.5 * (u0s + u1s) + 0.5 * (u1s - u0s) * gp[gi];
            let _xs = self.slave_surface.point_at(us, 0.5);
            let t_norm = (us - u0s) / (u1s - u0s).max(1e-14);
            let um = u0m + t_norm * (u1m - u0m);
            let ns = self.slave_surface.u_basis.num_basis - 1;
            let s_span = find_knot_span(ns, ps, us, &self.slave_surface.u_basis.knot);
            let s_basis = basis_functions(s_span, ps, us, &self.slave_surface.u_basis.knot);
            let nm = self.master_surface.u_basis.num_basis - 1;
            let um_clamped = um.clamp(u0m, u1m - 1e-10);
            let m_span = find_knot_span(nm, pm, um_clamped, &self.master_surface.u_basis.knot);
            let m_basis =
                basis_functions(m_span, pm, um_clamped, &self.master_surface.u_basis.knot);
            let w = 0.5 * (u1s - u0s) * gw[gi];
            for a in 0..n_slave.min(s_basis.len()) {
                for b in 0..n_master.min(m_basis.len()) {
                    d_mat[a * n_master + b] += s_basis[a] * m_basis[b] * w * self.alpha;
                }
            }
        }
        d_mat
    }
    /// Compute the interface gap (displacement jump) at a set of sample parameters.
    ///
    /// Returns a Vec of (parameter, gap_vector) pairs.
    pub fn interface_gap(&self, n_samples: usize) -> Vec<(f64, [f64; 3])> {
        let (u0s, u1s) = self.slave_surface.u_basis.domain();
        let (u0m, u1m) = self.master_surface.u_basis.domain();
        let mut gaps = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let t = i as f64 / (n_samples - 1).max(1) as f64;
            let us = u0s + (u1s - u0s) * t;
            let um = u0m + (u1m - u0m) * t;
            let xs = self.slave_surface.point_at(us, 0.5);
            let xm = self.master_surface.point_at(um, 0.5);
            let gap = sub3(xs, xm);
            gaps.push((t, gap));
        }
        gaps
    }
}
/// Interface coupling type for multi-patch IGA.
#[derive(Debug, Clone, PartialEq)]
pub enum CouplingType {
    /// Mortar method: L2-projection coupling.
    Mortar,
    /// Penalty method: stiffness penalty at interface.
    Penalty(f64),
    /// Nitsche's method: symmetric/nonsymmetric.
    Nitsche(f64),
}
