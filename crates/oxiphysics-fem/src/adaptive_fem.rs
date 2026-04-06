#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive FEM: a posteriori error estimation and mesh refinement.
//!
//! This module provides:
//!
//! - **ZzErrorEstimator**: Zienkiewicz-Zhu superconvergent patch recovery (ZZ/SPR).
//! - **ResidualEstimator**: Element residual-based error estimator.
//! - **HpRefinementStrategy**: hp-adaptive strategy combining h- and p-refinement.
//! - **MeshQuality**: aspect ratio, skewness, and Jacobian-based quality metrics.
//! - **AdaptiveQuadrature**: adaptive Gaussian quadrature with error control.
//! - **ConvergenceEstimator**: convergence rate estimation from successive meshes.
//! - **RefinementIndicator**: gradient-based and energy-based indicators.
//! - **CoarseningStrategy**: safe coarsening of over-refined regions.
//! - **ParallelAdaptiveRefinement**: work-stealing parallel h-refinement skeleton.
//! - **AdaptiveFemDriver**: top-level adaptive loop controller.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Math helpers (plain arrays, no nalgebra)
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Triangular element
// ---------------------------------------------------------------------------

/// A triangular (2D) or surface element for adaptive FEM.
#[derive(Clone, Debug)]
pub struct TriElement {
    /// Global node indices of the three vertices.
    pub nodes: [usize; 3],
    /// Positions of the three vertices \[\[x,y,0\\]; 3].
    pub coords: [[f64; 3]; 3],
    /// Polynomial degree p used in this element.
    pub degree: usize,
    /// Computed element diameter h_e.
    pub h_e: f64,
    /// Element area.
    pub area: f64,
}

impl TriElement {
    /// Create a new [`TriElement`] and compute h_e and area.
    pub fn new(nodes: [usize; 3], coords: [[f64; 3]; 3], degree: usize) -> Self {
        let e01 = sub3(coords[1], coords[0]);
        let e12 = sub3(coords[2], coords[1]);
        let e20 = sub3(coords[0], coords[2]);
        let h_e = len3(e01).max(len3(e12)).max(len3(e20));
        let cr = cross3(e01, sub3(coords[2], coords[0]));
        let area = 0.5 * len3(cr).abs();
        Self {
            nodes,
            coords,
            degree,
            h_e,
            area,
        }
    }

    /// Centroid of the element.
    pub fn centroid(&self) -> [f64; 3] {
        let c = self.coords;
        [
            (c[0][0] + c[1][0] + c[2][0]) / 3.0,
            (c[0][1] + c[1][1] + c[2][1]) / 3.0,
            (c[0][2] + c[1][2] + c[2][2]) / 3.0,
        ]
    }

    /// Midpoint of edge (i, j) where i,j in {0,1,2}.
    pub fn edge_midpoint(&self, i: usize, j: usize) -> [f64; 3] {
        let a = self.coords[i];
        let b = self.coords[j];
        scale3(add3(a, b), 0.5)
    }

    /// Longest edge index pair (i,j).
    pub fn longest_edge(&self) -> (usize, usize) {
        let l01 = len3(sub3(self.coords[1], self.coords[0]));
        let l12 = len3(sub3(self.coords[2], self.coords[1]));
        let l20 = len3(sub3(self.coords[0], self.coords[2]));
        if l01 >= l12 && l01 >= l20 {
            (0, 1)
        } else if l12 >= l20 {
            (1, 2)
        } else {
            (2, 0)
        }
    }
}

// ---------------------------------------------------------------------------
// Mesh quality metrics
// ---------------------------------------------------------------------------

/// Mesh quality metrics: aspect ratio, skewness, and minimum angle.
#[derive(Clone, Debug)]
pub struct MeshQuality {
    /// Aspect ratio: longest edge / shortest altitude.
    pub aspect_ratio: f64,
    /// Skewness: (optimal angle - actual minimum angle) / optimal angle.
    pub skewness: f64,
    /// Minimum interior angle \[radians\].
    pub min_angle: f64,
    /// Scaled Jacobian (1.0 = perfect, < 0 = inverted).
    pub scaled_jacobian: f64,
}

impl MeshQuality {
    /// Compute quality metrics for a triangle with vertices a, b, c.
    pub fn compute_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        let ab = sub3(b, a);
        let bc = sub3(c, b);
        let ca = sub3(a, c);
        let lab = len3(ab);
        let lbc = len3(bc);
        let lca = len3(ca);
        let longest = lab.max(lbc).max(lca);
        let _shortest = lab.min(lbc).min(lca);

        // Area
        let cr = cross3(ab, sub3(c, a));
        let area = 0.5 * len3(cr).abs();

        // Altitude from each vertex
        let alt_c = if lab > 1e-14 { 2.0 * area / lab } else { 0.0 };
        let alt_a = if lbc > 1e-14 { 2.0 * area / lbc } else { 0.0 };
        let alt_b = if lca > 1e-14 { 2.0 * area / lca } else { 0.0 };
        let min_alt = alt_c.min(alt_a).min(alt_b);
        let aspect_ratio = if min_alt > 1e-14 {
            longest / min_alt
        } else {
            f64::MAX
        };

        // Angles using dot product
        let cos_a = dot3(ab, scale3(ca, -1.0)) / (lab * lca).max(1e-14);
        let cos_b = dot3(scale3(ab, -1.0), bc) / (lab * lbc).max(1e-14);
        let cos_c = dot3(scale3(bc, -1.0), ca) / (lbc * lca).max(1e-14);
        let angle_a = cos_a.clamp(-1.0, 1.0).acos();
        let angle_b = cos_b.clamp(-1.0, 1.0).acos();
        let angle_c = cos_c.clamp(-1.0, 1.0).acos();
        let min_angle = angle_a.min(angle_b).min(angle_c);

        let optimal_angle = std::f64::consts::PI / 3.0; // 60° for equilateral
        let skewness = (optimal_angle - min_angle) / optimal_angle;

        // Scaled Jacobian ≈ 2*area / (lab*lbc) (simplified)
        let scaled_jacobian = if lab * lbc > 1e-14 {
            2.0 * area / (lab * lbc)
        } else {
            0.0
        };

        Self {
            aspect_ratio,
            skewness,
            min_angle,
            scaled_jacobian,
        }
    }

    /// Returns true if this element is acceptable (aspect_ratio < threshold).
    pub fn is_acceptable(&self, max_aspect_ratio: f64) -> bool {
        self.aspect_ratio < max_aspect_ratio && self.scaled_jacobian > 0.0
    }

    /// Quality score in \[0,1\]: 1 = perfect equilateral.
    pub fn quality_score(&self) -> f64 {
        (1.0 - self.skewness.clamp(0.0, 1.0)) * self.scaled_jacobian.clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Zienkiewicz-Zhu superconvergent patch recovery
// ---------------------------------------------------------------------------

/// Zienkiewicz-Zhu (ZZ-SPR) error estimator.
///
/// Computes recovered stress/gradient at nodes via least-squares patch fitting,
/// then estimates element errors from the difference between raw and recovered values.
#[derive(Clone, Debug)]
pub struct ZzErrorEstimator {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Number of elements.
    pub n_elements: usize,
    /// Recovered gradient at each node (x,y components).
    pub recovered_grad: Vec<[f64; 2]>,
    /// Element error indicators η_e.
    pub eta: Vec<f64>,
    /// Global error estimate (L2 norm of element errors).
    pub global_error: f64,
}

impl ZzErrorEstimator {
    /// Create a new [`ZzErrorEstimator`].
    pub fn new(n_nodes: usize, n_elements: usize) -> Self {
        Self {
            n_nodes,
            n_elements,
            recovered_grad: vec![[0.0; 2]; n_nodes],
            eta: vec![0.0; n_elements],
            global_error: 0.0,
        }
    }

    /// Perform SPR recovery: for each node, average gradients of surrounding elements.
    ///
    /// `elem_grads[e]` is the raw gradient \[dU/dx, dU/dy\] in element e.
    /// `connectivity[e]` lists the three node indices of element e.
    pub fn recover_gradients(&mut self, elem_grads: &[[f64; 2]], connectivity: &[[usize; 3]]) {
        let mut count = vec![0usize; self.n_nodes];
        for g in self.recovered_grad.iter_mut() {
            *g = [0.0; 2];
        }
        for (e, conn) in connectivity.iter().enumerate() {
            if e >= elem_grads.len() {
                break;
            }
            for &n in conn.iter() {
                if n < self.n_nodes {
                    self.recovered_grad[n][0] += elem_grads[e][0];
                    self.recovered_grad[n][1] += elem_grads[e][1];
                    count[n] += 1;
                }
            }
        }
        for n in 0..self.n_nodes {
            let c = count[n].max(1) as f64;
            self.recovered_grad[n][0] /= c;
            self.recovered_grad[n][1] /= c;
        }
    }

    /// Compute element error indicators.
    ///
    /// For each element, interpolates recovered gradient at centroid and
    /// computes ||σ* - σ_h||_e * sqrt(area).
    pub fn compute_element_errors(
        &mut self,
        elem_grads: &[[f64; 2]],
        connectivity: &[[usize; 3]],
        areas: &[f64],
    ) {
        let mut global_sq = 0.0;
        for e in 0..self.n_elements.min(elem_grads.len()) {
            let conn = connectivity[e];
            // Interpolate recovered gradient at centroid (average of node values)
            let rec = [
                (self.recovered_grad[conn[0].min(self.n_nodes - 1)][0]
                    + self.recovered_grad[conn[1].min(self.n_nodes - 1)][0]
                    + self.recovered_grad[conn[2].min(self.n_nodes - 1)][0])
                    / 3.0,
                (self.recovered_grad[conn[0].min(self.n_nodes - 1)][1]
                    + self.recovered_grad[conn[1].min(self.n_nodes - 1)][1]
                    + self.recovered_grad[conn[2].min(self.n_nodes - 1)][1])
                    / 3.0,
            ];
            let raw = elem_grads[e];
            let area = if e < areas.len() { areas[e] } else { 1.0 };
            let diff_sq = (rec[0] - raw[0]).powi(2) + (rec[1] - raw[1]).powi(2);
            self.eta[e] = (diff_sq * area).sqrt();
            global_sq += self.eta[e].powi(2);
        }
        self.global_error = global_sq.sqrt();
    }

    /// Effectivity index: ratio of estimated to exact error (approximate).
    ///
    /// `exact_error` should be supplied from an analytical solution if available.
    pub fn effectivity_index(&self, exact_error: f64) -> f64 {
        self.global_error / exact_error.max(1e-14)
    }

    /// Return indices of elements where η_e > θ * η_max (Dörfler marking).
    pub fn dorfer_mark(&self, theta: f64) -> Vec<usize> {
        let eta_max = self.eta.iter().cloned().fold(0.0f64, f64::max);
        self.eta
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, e)| *e > theta * eta_max)
            .map(|(i, _)| i)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Residual-based error estimator
// ---------------------------------------------------------------------------

/// Residual-based a posteriori error estimator (element residual + jump terms).
#[derive(Clone, Debug)]
pub struct ResidualEstimator {
    /// Number of elements.
    pub n_elements: usize,
    /// Interior element residual contribution R_K for each element.
    pub r_element: Vec<f64>,
    /// Edge jump contribution J_e for each element.
    pub r_jump: Vec<f64>,
    /// Combined error indicator η_e = h_e * R_K + h_e^0.5 * J_e.
    pub eta: Vec<f64>,
    /// Global error estimate.
    pub global_error: f64,
}

impl ResidualEstimator {
    /// Create a new [`ResidualEstimator`].
    pub fn new(n_elements: usize) -> Self {
        Self {
            n_elements,
            r_element: vec![0.0; n_elements],
            r_jump: vec![0.0; n_elements],
            eta: vec![0.0; n_elements],
            global_error: 0.0,
        }
    }

    /// Compute error indicators from residuals and jumps.
    ///
    /// `h_e[e]` is the element size, `r_e[e]` the volumetric residual magnitude,
    /// `j_e[e]` the edge flux jump magnitude.
    pub fn compute(&mut self, h_e: &[f64], r_e: &[f64], j_e: &[f64]) {
        let mut global_sq = 0.0;
        for k in 0..self.n_elements {
            let h = if k < h_e.len() { h_e[k] } else { 1.0 };
            let r = if k < r_e.len() { r_e[k] } else { 0.0 };
            let j = if k < j_e.len() { j_e[k] } else { 0.0 };
            self.r_element[k] = h * r;
            self.r_jump[k] = h.sqrt() * j;
            self.eta[k] = self.r_element[k] + self.r_jump[k];
            global_sq += self.eta[k].powi(2);
        }
        self.global_error = global_sq.sqrt();
    }

    /// Maximum element indicator.
    pub fn max_eta(&self) -> f64 {
        self.eta.iter().cloned().fold(0.0f64, f64::max)
    }

    /// Fraction of elements exceeding threshold `tol`.
    pub fn fraction_exceeding(&self, tol: f64) -> f64 {
        let cnt = self.eta.iter().filter(|&&e| e > tol).count();
        cnt as f64 / self.n_elements as f64
    }
}

// ---------------------------------------------------------------------------
// H-refinement: longest-edge bisection
// ---------------------------------------------------------------------------

/// Result of h-refinement: new nodes and connectivity.
#[derive(Clone, Debug)]
pub struct HRefinementResult {
    /// New node coordinates (existing + midpoints).
    pub nodes: Vec<[f64; 3]>,
    /// New triangle connectivity.
    pub elements: Vec<[usize; 3]>,
    /// Polynomial degrees for new elements.
    pub degrees: Vec<usize>,
    /// Number of refinements performed.
    pub n_refined: usize,
}

/// Perform longest-edge bisection h-refinement on marked elements.
///
/// `nodes` are the current nodal coordinates, `elements` are triangles,
/// `degrees` is per-element polynomial degree, `marked` are element indices to refine.
#[allow(clippy::too_many_arguments)]
pub fn h_refine_longest_edge(
    nodes: &[[f64; 3]],
    elements: &[[usize; 3]],
    degrees: &[usize],
    marked: &[usize],
) -> HRefinementResult {
    let mut new_nodes: Vec<[f64; 3]> = nodes.to_vec();
    let mut new_elements: Vec<[usize; 3]> = Vec::new();
    let mut new_degrees: Vec<usize> = Vec::new();
    let mut n_refined = 0;

    // Track midpoint nodes to avoid duplicates: edge (min,max) -> node index
    let mut midpoints: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();

    let marked_set: std::collections::HashSet<usize> = marked.iter().cloned().collect();

    for (e, &conn) in elements.iter().enumerate() {
        let deg = if e < degrees.len() { degrees[e] } else { 1 };
        if !marked_set.contains(&e) {
            new_elements.push(conn);
            new_degrees.push(deg);
            continue;
        }
        // Find longest edge
        let elem = TriElement::new(conn, [nodes[conn[0]], nodes[conn[1]], nodes[conn[2]]], deg);
        let (i, j) = elem.longest_edge();
        let ni = conn[i];
        let nj = conn[j];
        let key = (ni.min(nj), ni.max(nj));

        let mid_idx = if let Some(&m) = midpoints.get(&key) {
            m
        } else {
            let mid = scale3(add3(nodes[ni], nodes[nj]), 0.5);
            let m = new_nodes.len();
            new_nodes.push(mid);
            midpoints.insert(key, m);
            m
        };

        // k = third vertex (not i or j)
        let k = 3 - i - j; // works for {0,1,2} if i+j = 0+1=1 → k=2, etc.
        let nk = conn[k];

        new_elements.push([ni, mid_idx, nk]);
        new_elements.push([mid_idx, nj, nk]);
        new_degrees.push(deg);
        new_degrees.push(deg);
        n_refined += 1;
    }

    HRefinementResult {
        nodes: new_nodes,
        elements: new_elements,
        degrees: new_degrees,
        n_refined,
    }
}

// ---------------------------------------------------------------------------
// P-refinement
// ---------------------------------------------------------------------------

/// P-refinement: increase polynomial degree on marked elements.
#[derive(Clone, Debug)]
pub struct PRefinement {
    /// Maximum allowed polynomial degree.
    pub p_max: usize,
    /// Current polynomial degrees per element.
    pub degrees: Vec<usize>,
}

impl PRefinement {
    /// Create a new [`PRefinement`] with uniform degree `p0`.
    pub fn new(n_elements: usize, p0: usize, p_max: usize) -> Self {
        Self {
            p_max,
            degrees: vec![p0; n_elements],
        }
    }

    /// Increase degree by 1 on marked elements (capped at p_max).
    pub fn refine_marked(&mut self, marked: &[usize]) {
        for &e in marked {
            if e < self.degrees.len() {
                self.degrees[e] = (self.degrees[e] + 1).min(self.p_max);
            }
        }
    }

    /// Number of DOF per element for degree p (triangular: (p+1)(p+2)/2).
    pub fn dof_per_element(p: usize) -> usize {
        (p + 1) * (p + 2) / 2
    }

    /// Total DOFs (approximate; counts element-interior DOFs only).
    pub fn total_interior_dof(&self) -> usize {
        self.degrees.iter().map(|&p| Self::dof_per_element(p)).sum()
    }
}

// ---------------------------------------------------------------------------
// HP-refinement strategy
// ---------------------------------------------------------------------------

/// Smoothness indicator threshold type.
#[derive(Clone, Debug)]
pub enum SmoothnessIndicator {
    /// Solution is smooth in this element: prefer p-refinement.
    Smooth,
    /// Solution is non-smooth (e.g., has a singularity): prefer h-refinement.
    NonSmooth,
}

/// HP-adaptive strategy combining h- and p-refinement.
#[derive(Clone, Debug)]
pub struct HpRefinementStrategy {
    /// Error tolerances per element.
    pub eta: Vec<f64>,
    /// Per-element polynomial degrees.
    pub degrees: Vec<usize>,
    /// Smoothness classification per element.
    pub smoothness: Vec<SmoothnessIndicator>,
    /// Error threshold θ for marking.
    pub theta: f64,
    /// Maximum polynomial degree.
    pub p_max: usize,
}

impl HpRefinementStrategy {
    /// Create a new [`HpRefinementStrategy`].
    pub fn new(n_elements: usize, theta: f64, p_max: usize) -> Self {
        Self {
            eta: vec![0.0; n_elements],
            degrees: vec![1; n_elements],
            smoothness: (0..n_elements)
                .map(|_| SmoothnessIndicator::Smooth)
                .collect(),
            theta,
            p_max,
        }
    }

    /// Classify smoothness using a decay rate of Legendre coefficients.
    ///
    /// `coeff_decay[e]` is the estimated coefficient decay rate for element e.
    /// Positive decay → smooth; negative or slow → non-smooth.
    pub fn classify_smoothness(&mut self, coeff_decay: &[f64]) {
        for (e, s) in self.smoothness.iter_mut().enumerate() {
            let decay = if e < coeff_decay.len() {
                coeff_decay[e]
            } else {
                1.0
            };
            *s = if decay > 0.5 {
                SmoothnessIndicator::Smooth
            } else {
                SmoothnessIndicator::NonSmooth
            };
        }
    }

    /// Mark elements for refinement (Dörfler strategy).
    pub fn mark_elements(&self) -> Vec<usize> {
        let eta_max = self.eta.iter().cloned().fold(0.0f64, f64::max);
        self.eta
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, e)| *e > self.theta * eta_max)
            .map(|(i, _)| i)
            .collect()
    }

    /// Partition marked elements into h-refine and p-refine sets.
    pub fn partition_refinement(&self, marked: &[usize]) -> (Vec<usize>, Vec<usize>) {
        let mut h_set = Vec::new();
        let mut p_set = Vec::new();
        for &e in marked {
            match &self.smoothness[e] {
                SmoothnessIndicator::Smooth if self.degrees[e] < self.p_max => {
                    p_set.push(e);
                }
                _ => {
                    h_set.push(e);
                }
            }
        }
        (h_set, p_set)
    }

    /// Expected DOF reduction factor when using hp vs h-only.
    pub fn hp_efficiency_gain(&self, h_only_dof: usize, hp_dof: usize) -> f64 {
        if hp_dof == 0 {
            0.0
        } else {
            h_only_dof as f64 / hp_dof as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Adaptive quadrature
// ---------------------------------------------------------------------------

/// Gauss-Legendre quadrature points and weights on \[-1, 1\] (n=5).
const GL_POINTS_5: [f64; 5] = [
    -0.906_179_845_938_664,
    -0.538_469_310_105_683,
    0.0,
    0.538_469_310_105_683,
    0.906_179_845_938_664,
];
const GL_WEIGHTS_5: [f64; 5] = [
    0.236_926_885_056_189,
    0.478_628_670_499_367,
    0.568_888_888_888_889,
    0.478_628_670_499_367,
    0.236_926_885_056_189,
];

/// Adaptive quadrature using recursive interval bisection.
#[derive(Clone, Debug)]
pub struct AdaptiveQuadrature {
    /// Absolute tolerance for local error.
    pub tol: f64,
    /// Maximum recursion depth.
    pub max_depth: usize,
    /// Total function evaluations performed.
    pub n_evals: usize,
}

impl AdaptiveQuadrature {
    /// Create a new [`AdaptiveQuadrature`].
    pub fn new(tol: f64, max_depth: usize) -> Self {
        Self {
            tol,
            max_depth,
            n_evals: 0,
        }
    }

    /// Integrate `f` over \[a, b\] using GL5 with adaptive bisection.
    pub fn integrate<F: Fn(f64) -> f64>(&mut self, f: &F, a: f64, b: f64) -> f64 {
        self.n_evals = 0;
        self.integrate_recursive(f, a, b, 0)
    }

    fn integrate_recursive<F: Fn(f64) -> f64>(
        &mut self,
        f: &F,
        a: f64,
        b: f64,
        depth: usize,
    ) -> f64 {
        let coarse = self.gl5(f, a, b);
        if depth >= self.max_depth {
            return coarse;
        }
        let mid = 0.5 * (a + b);
        let fine = self.gl5(f, a, mid) + self.gl5(f, mid, b);
        let err = (fine - coarse).abs();
        if err < self.tol * (1.0 + fine.abs()) {
            fine
        } else {
            self.integrate_recursive(f, a, mid, depth + 1)
                + self.integrate_recursive(f, mid, b, depth + 1)
        }
    }

    fn gl5<F: Fn(f64) -> f64>(&mut self, f: &F, a: f64, b: f64) -> f64 {
        let half = 0.5 * (b - a);
        let center = 0.5 * (a + b);
        let mut sum = 0.0;
        for i in 0..5 {
            let x = center + half * GL_POINTS_5[i];
            sum += GL_WEIGHTS_5[i] * f(x);
            self.n_evals += 1;
        }
        sum * half
    }

    /// 2D integration over triangle (0,0)-(1,0)-(0,1) using Dunavant quadrature (3-point).
    pub fn integrate_triangle<F: Fn(f64, f64) -> f64>(&mut self, f: &F, area: f64) -> f64 {
        // 3-point Dunavant rule on reference triangle
        let pts = [
            (1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0),
        ];
        let w = 1.0 / 3.0;
        let sum: f64 = pts.iter().map(|&(u, v)| f(u, v)).sum();
        self.n_evals += 3;
        sum * w * area
    }
}

// ---------------------------------------------------------------------------
// Convergence rate estimator
// ---------------------------------------------------------------------------

/// Estimates algebraic convergence rate from successive mesh refinements.
#[derive(Clone, Debug)]
pub struct ConvergenceEstimator {
    /// History of (n_dof, error) pairs.
    pub history: Vec<(usize, f64)>,
}

impl ConvergenceEstimator {
    /// Create a new [`ConvergenceEstimator`].
    pub fn new() -> Self {
        Self {
            history: Vec::new(),
        }
    }

    /// Record a new data point (n_dof, error).
    pub fn record(&mut self, n_dof: usize, error: f64) {
        self.history.push((n_dof, error));
    }

    /// Estimate convergence rate p from the last two data points.
    ///
    /// Assumes error ~ C * h^p ~ C * N^{-p/d} (d = dimension).
    /// Returns (rate, C) or None if insufficient data.
    pub fn convergence_rate(&self, dimension: usize) -> Option<(f64, f64)> {
        let n = self.history.len();
        if n < 2 {
            return None;
        }
        let (n1, e1) = self.history[n - 2];
        let (n2, e2) = self.history[n - 1];
        if n1 == 0 || n2 == 0 || e1 <= 0.0 || e2 <= 0.0 {
            return None;
        }
        let d = dimension as f64;
        let rate = -(e2 / e1).ln() / ((n2 as f64 / n1 as f64).ln() / d);
        let h1 = (n1 as f64).powf(-1.0 / d);
        let c = e1 / h1.powf(rate);
        Some((rate, c))
    }

    /// Predict error for target DOF count.
    pub fn predict_error(&self, n_target: usize, dimension: usize) -> Option<f64> {
        let (rate, c) = self.convergence_rate(dimension)?;
        let h = (n_target as f64).powf(-1.0 / dimension as f64);
        Some(c * h.powf(rate))
    }

    /// Determine if convergence goal `tol` is reached.
    pub fn is_converged(&self, tol: f64) -> bool {
        self.history.last().map(|(_, e)| *e < tol).unwrap_or(false)
    }
}

impl Default for ConvergenceEstimator {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Refinement indicators
// ---------------------------------------------------------------------------

/// Gradient-based refinement indicator for stress-driven or energy-driven adaptivity.
#[derive(Clone, Debug)]
pub struct RefinementIndicator {
    /// Number of elements.
    pub n_elements: usize,
    /// Per-element stress gradient magnitude ||∇σ||.
    pub stress_grad: Vec<f64>,
    /// Per-element energy error density.
    pub energy_error: Vec<f64>,
    /// Combined indicator.
    pub indicator: Vec<f64>,
    /// Weighting between stress gradient and energy error (0=energy only, 1=stress only).
    pub alpha: f64,
}

impl RefinementIndicator {
    /// Create a new [`RefinementIndicator`].
    pub fn new(n_elements: usize, alpha: f64) -> Self {
        Self {
            n_elements,
            stress_grad: vec![0.0; n_elements],
            energy_error: vec![0.0; n_elements],
            indicator: vec![0.0; n_elements],
            alpha,
        }
    }

    /// Compute the combined indicator from stress gradients and energy errors.
    pub fn compute(&mut self) {
        let sg_max = self
            .stress_grad
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
            .max(1e-14);
        let ee_max = self
            .energy_error
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
            .max(1e-14);
        for e in 0..self.n_elements {
            let sg_norm = self.stress_grad[e] / sg_max;
            let ee_norm = self.energy_error[e] / ee_max;
            self.indicator[e] = self.alpha * sg_norm + (1.0 - self.alpha) * ee_norm;
        }
    }

    /// Mark elements where indicator > theta.
    pub fn mark(&self, theta: f64) -> Vec<usize> {
        self.indicator
            .iter()
            .cloned()
            .enumerate()
            .filter(|(_, v)| *v > theta)
            .map(|(i, _)| i)
            .collect()
    }

    /// Compute energy error density e_e = (sigma_h - sigma_rec)^T * D^{-1} * (sigma_h - sigma_rec) * area.
    ///
    /// For simplicity uses scalar energy: e_e = (err_x^2 + err_y^2) * area.
    pub fn compute_energy_error(
        &mut self,
        elem_stress: &[[f64; 2]],
        rec_stress: &[[f64; 2]],
        areas: &[f64],
    ) {
        for e in 0..self.n_elements.min(elem_stress.len()) {
            let area = if e < areas.len() { areas[e] } else { 1.0 };
            let dx = elem_stress[e][0] - rec_stress[e][0];
            let dy = elem_stress[e][1] - rec_stress[e][1];
            self.energy_error[e] = (dx * dx + dy * dy) * area;
        }
    }
}

// ---------------------------------------------------------------------------
// Coarsening strategy
// ---------------------------------------------------------------------------

/// Coarsening strategy: merge over-refined elements to reduce DOF count.
#[derive(Clone, Debug)]
pub struct CoarseningStrategy {
    /// Coarsening threshold: elements with η_e < tol_coarsen are candidates.
    pub tol_coarsen: f64,
    /// Minimum polynomial degree.
    pub p_min: usize,
    /// Minimum element size h_min.
    pub h_min: f64,
}

impl CoarseningStrategy {
    /// Create a new [`CoarseningStrategy`].
    pub fn new(tol_coarsen: f64, p_min: usize, h_min: f64) -> Self {
        Self {
            tol_coarsen,
            p_min,
            h_min,
        }
    }

    /// Mark elements eligible for coarsening (η_e < tol and h_e > h_min).
    pub fn mark_coarsen(&self, eta: &[f64], h_e: &[f64]) -> Vec<usize> {
        eta.iter()
            .enumerate()
            .filter(|(e, err)| **err < self.tol_coarsen && *e < h_e.len() && h_e[*e] > self.h_min)
            .map(|(i, _)| i)
            .collect()
    }

    /// Reduce polynomial degree by 1 on marked elements (floor at p_min).
    pub fn p_coarsen(&self, degrees: &mut [usize], marked: &[usize]) {
        for &e in marked {
            if e < degrees.len() && degrees[e] > self.p_min {
                degrees[e] -= 1;
            }
        }
    }

    /// Compute DOF savings from coarsening (approximate).
    pub fn dof_savings(&self, degrees: &[usize], marked: &[usize]) -> usize {
        marked
            .iter()
            .filter(|&&e| e < degrees.len() && degrees[e] > self.p_min)
            .map(|&e| {
                let p = degrees[e];
                PRefinement::dof_per_element(p) - PRefinement::dof_per_element(p - 1)
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Parallel adaptive refinement skeleton
// ---------------------------------------------------------------------------

/// Work item for parallel h-refinement.
#[derive(Clone, Debug)]
pub struct RefinementTask {
    /// Element index to be refined.
    pub element_idx: usize,
    /// Priority (higher = refine first).
    pub priority: f64,
}

/// Parallel adaptive refinement coordinator (work-stealing skeleton).
///
/// In a real implementation this would use a thread pool; here it processes
/// tasks sequentially for correctness and portability.
#[derive(Clone, Debug)]
pub struct ParallelAdaptiveRefinement {
    /// Pending refinement tasks.
    pub tasks: Vec<RefinementTask>,
    /// Number of threads (stored for future parallelism).
    pub n_threads: usize,
    /// Number of completed refinements.
    pub n_completed: usize,
}

impl ParallelAdaptiveRefinement {
    /// Create a new [`ParallelAdaptiveRefinement`].
    pub fn new(n_threads: usize) -> Self {
        Self {
            tasks: Vec::new(),
            n_threads,
            n_completed: 0,
        }
    }

    /// Queue elements for refinement.
    pub fn queue(&mut self, marked: &[usize], eta: &[f64]) {
        for &e in marked {
            let priority = if e < eta.len() { eta[e] } else { 0.0 };
            self.tasks.push(RefinementTask {
                element_idx: e,
                priority,
            });
        }
        // Sort by descending priority
        self.tasks.sort_by(|a, b| {
            b.priority
                .partial_cmp(&a.priority)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Process all queued tasks (sequential execution).
    ///
    /// Returns element indices in processing order.
    pub fn process_all(&mut self) -> Vec<usize> {
        let result: Vec<usize> = self.tasks.iter().map(|t| t.element_idx).collect();
        self.n_completed += self.tasks.len();
        self.tasks.clear();
        result
    }

    /// Estimated parallel speedup (Amdahl's law, assuming 90% parallelizable).
    pub fn amdahl_speedup(&self) -> f64 {
        let p = 0.9;
        let n = self.n_threads as f64;
        1.0 / ((1.0 - p) + p / n)
    }
}

// ---------------------------------------------------------------------------
// Error norm computation
// ---------------------------------------------------------------------------

/// Compute L2 norm of error given exact and approximate solutions at quadrature points.
pub fn l2_error_norm(exact: &[f64], approx: &[f64], weights: &[f64]) -> f64 {
    let sum: f64 = exact
        .iter()
        .zip(approx.iter())
        .zip(weights.iter())
        .map(|((e, a), w)| w * (e - a).powi(2))
        .sum();
    sum.sqrt()
}

/// Compute H1 seminorm of error.
///
/// `grad_exact` and `grad_approx` are \[N x 2\] gradient arrays (flattened as \[grad_x, grad_y\] pairs).
pub fn h1_seminorm_error(
    grad_exact: &[[f64; 2]],
    grad_approx: &[[f64; 2]],
    weights: &[f64],
) -> f64 {
    let sum: f64 = grad_exact
        .iter()
        .zip(grad_approx.iter())
        .zip(weights.iter())
        .map(|((ge, ga), w)| {
            let dx = ge[0] - ga[0];
            let dy = ge[1] - ga[1];
            w * (dx * dx + dy * dy)
        })
        .sum();
    sum.sqrt()
}

// ---------------------------------------------------------------------------
// Optimal refinement ratio predictor
// ---------------------------------------------------------------------------

/// Given current error `e0` and convergence rate `p`, compute the element size
/// ratio needed to achieve target error `e_target` in one refinement step.
///
/// Returns h_new / h_old.
pub fn optimal_refinement_ratio(e0: f64, e_target: f64, convergence_rate: f64) -> f64 {
    if e0 <= 0.0 || e_target <= 0.0 || convergence_rate <= 0.0 {
        return 0.5;
    }
    (e_target / e0).powf(1.0 / convergence_rate)
}

// ---------------------------------------------------------------------------
// Adaptive FEM driver
// ---------------------------------------------------------------------------

/// Top-level adaptive FEM loop controller.
///
/// Orchestrates error estimation, marking, refinement, and convergence checks.
#[derive(Clone, Debug)]
pub struct AdaptiveFemDriver {
    /// Error tolerance for stopping criterion.
    pub tol: f64,
    /// Maximum number of adaptive iterations.
    pub max_iter: usize,
    /// Dörfler marking fraction θ.
    pub theta_mark: f64,
    /// Coarsening fraction.
    pub theta_coarsen: f64,
    /// Current iteration counter.
    pub iteration: usize,
    /// Convergence history.
    pub convergence: ConvergenceEstimator,
    /// Use hp-refinement (true) or h-only (false).
    pub use_hp: bool,
}

impl AdaptiveFemDriver {
    /// Create a new [`AdaptiveFemDriver`].
    pub fn new(
        tol: f64,
        max_iter: usize,
        theta_mark: f64,
        theta_coarsen: f64,
        use_hp: bool,
    ) -> Self {
        Self {
            tol,
            max_iter,
            theta_mark,
            theta_coarsen,
            iteration: 0,
            convergence: ConvergenceEstimator::new(),
            use_hp,
        }
    }

    /// Check stopping criterion.
    pub fn should_stop(&self, global_error: f64) -> bool {
        self.iteration >= self.max_iter || global_error < self.tol
    }

    /// Record an adaptive iteration result.
    pub fn record_iteration(&mut self, n_dof: usize, global_error: f64) {
        self.convergence.record(n_dof, global_error);
        self.iteration += 1;
    }

    /// Compute suggested new error tolerance for next refinement level.
    ///
    /// Uses predictor-corrector: reduce by factor sqrt(theta_mark).
    pub fn next_tolerance(&self) -> f64 {
        self.tol * self.theta_mark.sqrt()
    }

    /// Return estimated number of remaining iterations to convergence.
    pub fn estimated_remaining_iters(&self, global_error: f64) -> usize {
        if global_error <= self.tol {
            return 0;
        }
        let ratio = global_error / self.tol;
        // Assume 50% error reduction per iteration
        (ratio.log2() / 1.0_f64.log2()).ceil() as usize + 1
    }
}

// ---------------------------------------------------------------------------
// Recovery-based stress smoothing
// ---------------------------------------------------------------------------

/// Recovered (smoothed) stress field using node-averaging SPR.
#[derive(Clone, Debug)]
pub struct RecoveredStress {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Recovered stress components \[σ_xx, σ_yy, σ_xy\] at each node.
    pub sigma: Vec<[f64; 3]>,
}

impl RecoveredStress {
    /// Create a new [`RecoveredStress`].
    pub fn new(n_nodes: usize) -> Self {
        Self {
            n_nodes,
            sigma: vec![[0.0; 3]; n_nodes],
        }
    }

    /// Recover stress by averaging element stresses at shared nodes.
    ///
    /// `elem_stress[e]` is \[σ_xx, σ_yy, σ_xy\] for element e.
    pub fn recover(&mut self, elem_stress: &[[f64; 3]], connectivity: &[[usize; 3]]) {
        let mut count = vec![0usize; self.n_nodes];
        for s in self.sigma.iter_mut() {
            *s = [0.0; 3];
        }

        for (e, conn) in connectivity.iter().enumerate() {
            if e >= elem_stress.len() {
                break;
            }
            for &n in conn.iter() {
                if n < self.n_nodes {
                    self.sigma[n][0] += elem_stress[e][0];
                    self.sigma[n][1] += elem_stress[e][1];
                    self.sigma[n][2] += elem_stress[e][2];
                    count[n] += 1;
                }
            }
        }
        for n in 0..self.n_nodes {
            let c = count[n].max(1) as f64;
            self.sigma[n][0] /= c;
            self.sigma[n][1] /= c;
            self.sigma[n][2] /= c;
        }
    }

    /// Von Mises equivalent stress at node n.
    pub fn von_mises(&self, n: usize) -> f64 {
        if n >= self.n_nodes {
            return 0.0;
        }
        let s = self.sigma[n];
        let vm_sq = s[0] * s[0] - s[0] * s[1] + s[1] * s[1] + 3.0 * s[2] * s[2];
        vm_sq.max(0.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Anisotropic refinement metric
// ---------------------------------------------------------------------------

/// Anisotropic mesh metric tensor M for directional refinement.
///
/// M = R * diag(1/h1^2, 1/h2^2) * R^T where R is the principal direction matrix.
#[derive(Clone, Debug)]
pub struct AnisotropicMetric {
    /// Principal stretch in direction 1 \[m\].
    pub h1: f64,
    /// Principal stretch in direction 2 \[m\].
    pub h2: f64,
    /// Principal direction angle \[rad\].
    pub angle: f64,
}

impl AnisotropicMetric {
    /// Create a new [`AnisotropicMetric`].
    pub fn new(h1: f64, h2: f64, angle: f64) -> Self {
        Self { h1, h2, angle }
    }

    /// Evaluate the metric at a point: returns 2x2 matrix as \[m00, m01, m10, m11\].
    pub fn metric_tensor(&self) -> [f64; 4] {
        let ca = self.angle.cos();
        let sa = self.angle.sin();
        let l1 = 1.0 / self.h1.powi(2).max(1e-20);
        let l2 = 1.0 / self.h2.powi(2).max(1e-20);
        let m00 = l1 * ca * ca + l2 * sa * sa;
        let m01 = (l1 - l2) * ca * sa;
        let m11 = l1 * sa * sa + l2 * ca * ca;
        [m00, m01, m01, m11]
    }

    /// Edge length in the metric: sqrt(e^T M e) where e is edge vector.
    pub fn edge_length_metric(&self, e: [f64; 2]) -> f64 {
        let m = self.metric_tensor();
        let me0 = m[0] * e[0] + m[1] * e[1];
        let me1 = m[2] * e[0] + m[3] * e[1];
        (e[0] * me0 + e[1] * me1).max(0.0).sqrt()
    }

    /// Anisotropy ratio h1/h2.
    pub fn anisotropy_ratio(&self) -> f64 {
        self.h1 / self.h2.max(1e-14)
    }
}

// ---------------------------------------------------------------------------
// Interpolation error estimate
// ---------------------------------------------------------------------------

/// Interpolation error estimate: ||u - u_h||_K ≤ C * h_K^{p+1} * |u|_{p+1,K}.
///
/// Returns the estimated interpolation error for element K.
pub fn interpolation_error_estimate(
    h_e: f64,
    degree: usize,
    seminorm_p1: f64,
    constant_c: f64,
) -> f64 {
    constant_c * h_e.powi((degree + 1) as i32) * seminorm_p1
}

/// Compute element size from area: h_e = sqrt(2 * area) (equilateral triangle scaling).
pub fn element_size_from_area(area: f64) -> f64 {
    (2.0 * area).sqrt()
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tri_element_area() {
        let coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let e = TriElement::new([0, 1, 2], coords, 1);
        assert!((e.area - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_tri_element_centroid() {
        let coords = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let e = TriElement::new([0, 1, 2], coords, 1);
        let c = e.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tri_element_h_e() {
        let coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let e = TriElement::new([0, 1, 2], coords, 1);
        assert!(e.h_e > 1.0 && e.h_e < 2.0); // hypotenuse sqrt(2)
    }

    #[test]
    fn test_mesh_quality_equilateral() {
        let s = 3.0_f64.sqrt() / 2.0;
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, s, 0.0];
        let q = MeshQuality::compute_triangle(a, b, c);
        assert!(q.aspect_ratio < 1.2);
        assert!(q.quality_score() > 0.8);
    }

    #[test]
    fn test_mesh_quality_is_acceptable() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.5, 0.0];
        let q = MeshQuality::compute_triangle(a, b, c);
        assert!(q.is_acceptable(10.0));
    }

    #[test]
    fn test_mesh_quality_degenerate() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0]; // collinear
        let q = MeshQuality::compute_triangle(a, b, c);
        assert!(!q.is_acceptable(10.0));
    }

    #[test]
    fn test_zz_estimator_recover_uniform() {
        let mut zz = ZzErrorEstimator::new(3, 2);
        let grads = [[1.0f64, 0.0], [1.0, 0.0]];
        let conn = [[0usize, 1, 2], [0, 1, 2]];
        zz.recover_gradients(&grads, &conn);
        assert!((zz.recovered_grad[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zz_estimator_element_errors_zero_for_exact() {
        let mut zz = ZzErrorEstimator::new(3, 2);
        let grads = [[1.0f64, 0.0], [1.0, 0.0]];
        let conn = [[0usize, 1, 2], [0, 1, 2]];
        let areas = [0.5, 0.5];
        zz.recover_gradients(&grads, &conn);
        zz.compute_element_errors(&grads, &conn, &areas);
        assert!(zz.global_error < 1e-10);
    }

    #[test]
    fn test_zz_dorfer_mark() {
        let mut zz = ZzErrorEstimator::new(4, 4);
        zz.eta = vec![0.1, 0.5, 0.2, 0.9];
        let marked = zz.dorfer_mark(0.5);
        assert!(marked.contains(&3));
        assert!(marked.contains(&1));
    }

    #[test]
    fn test_residual_estimator_compute() {
        let mut est = ResidualEstimator::new(3);
        let h = [0.1, 0.1, 0.1];
        let r = [1.0, 2.0, 3.0];
        let j = [0.5, 0.5, 0.5];
        est.compute(&h, &r, &j);
        assert!(est.global_error > 0.0);
        assert!(est.max_eta() == est.eta[2]);
    }

    #[test]
    fn test_h_refine_longest_edge() {
        let nodes: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let elements: Vec<[usize; 3]> = vec![[0, 1, 2]];
        let degrees = vec![1usize];
        let result = h_refine_longest_edge(&nodes, &elements, &degrees, &[0]);
        assert_eq!(result.n_refined, 1);
        assert!(result.nodes.len() > 3);
        assert_eq!(result.elements.len(), 2);
    }

    #[test]
    fn test_h_refine_unmarked_unchanged() {
        let nodes: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let elements: Vec<[usize; 3]> = vec![[0, 1, 2]];
        let degrees = vec![1usize];
        let result = h_refine_longest_edge(&nodes, &elements, &degrees, &[]);
        assert_eq!(result.n_refined, 0);
        assert_eq!(result.elements.len(), 1);
    }

    #[test]
    fn test_p_refinement_increase() {
        let mut pr = PRefinement::new(4, 1, 5);
        pr.refine_marked(&[0, 2]);
        assert_eq!(pr.degrees[0], 2);
        assert_eq!(pr.degrees[1], 1);
        assert_eq!(pr.degrees[2], 2);
    }

    #[test]
    fn test_p_refinement_cap() {
        let mut pr = PRefinement::new(2, 5, 5);
        pr.refine_marked(&[0]);
        assert_eq!(pr.degrees[0], 5); // Capped
    }

    #[test]
    fn test_p_refinement_dof() {
        // p=1: 3 nodes; p=2: 6 nodes
        assert_eq!(PRefinement::dof_per_element(1), 3);
        assert_eq!(PRefinement::dof_per_element(2), 6);
    }

    #[test]
    fn test_hp_strategy_partition() {
        let mut hp = HpRefinementStrategy::new(4, 0.5, 5);
        hp.eta = vec![0.1, 0.9, 0.3, 0.8];
        hp.smoothness[1] = SmoothnessIndicator::Smooth;
        hp.smoothness[3] = SmoothnessIndicator::NonSmooth;
        let marked = hp.mark_elements();
        let (h_set, p_set) = hp.partition_refinement(&marked);
        assert!(!h_set.is_empty() || !p_set.is_empty());
    }

    #[test]
    fn test_adaptive_quadrature_x2() {
        let mut aq = AdaptiveQuadrature::new(1e-10, 20);
        let result = aq.integrate(&|x: f64| x * x, 0.0, 1.0);
        assert!((result - 1.0 / 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_adaptive_quadrature_triangle() {
        let mut aq = AdaptiveQuadrature::new(1e-10, 20);
        // Integral of 1 over reference triangle (area = 0.5)
        let result = aq.integrate_triangle(&|_u, _v| 1.0, 0.5);
        assert!((result - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_convergence_estimator_rate() {
        let mut ce = ConvergenceEstimator::new();
        ce.record(100, 0.1);
        ce.record(400, 0.05);
        let (rate, _) = ce.convergence_rate(2).unwrap();
        assert!(rate > 0.0);
    }

    #[test]
    fn test_convergence_estimator_is_converged() {
        let mut ce = ConvergenceEstimator::new();
        ce.record(100, 1e-10);
        assert!(ce.is_converged(1e-8));
    }

    #[test]
    fn test_convergence_estimator_insufficient() {
        let ce = ConvergenceEstimator::new();
        assert!(ce.convergence_rate(2).is_none());
    }

    #[test]
    fn test_refinement_indicator_mark() {
        let mut ri = RefinementIndicator::new(4, 0.5);
        ri.stress_grad = vec![0.1, 0.9, 0.2, 0.8];
        ri.energy_error = vec![0.1, 0.9, 0.2, 0.8];
        ri.compute();
        let marked = ri.mark(0.5);
        assert!(!marked.is_empty());
    }

    #[test]
    fn test_coarsening_strategy_mark() {
        let cs = CoarseningStrategy::new(0.01, 1, 0.01);
        let eta = vec![0.001, 0.1, 0.005, 0.2];
        let h = vec![0.05, 0.05, 0.05, 0.05];
        let marked = cs.mark_coarsen(&eta, &h);
        assert!(marked.contains(&0));
        assert!(marked.contains(&2));
        assert!(!marked.contains(&1));
    }

    #[test]
    fn test_coarsening_strategy_p_coarsen() {
        let cs = CoarseningStrategy::new(0.01, 1, 0.01);
        let mut degrees = vec![3, 2, 3, 1];
        cs.p_coarsen(&mut degrees, &[0, 2]);
        assert_eq!(degrees[0], 2);
        assert_eq!(degrees[3], 1); // At p_min, unchanged
    }

    #[test]
    fn test_parallel_adaptive_queue_process() {
        let mut par = ParallelAdaptiveRefinement::new(4);
        par.queue(&[0, 2, 3], &[0.1, 0.5, 0.3, 0.9]);
        let processed = par.process_all();
        assert_eq!(processed.len(), 3);
        assert_eq!(par.n_completed, 3);
    }

    #[test]
    fn test_parallel_adaptive_speedup() {
        let par = ParallelAdaptiveRefinement::new(4);
        let s = par.amdahl_speedup();
        assert!(s > 1.0 && s < 5.0);
    }

    #[test]
    fn test_l2_error_norm() {
        let exact = [1.0, 2.0, 3.0];
        let approx = [1.0, 2.0, 3.0];
        let weights = [1.0, 1.0, 1.0];
        assert!(l2_error_norm(&exact, &approx, &weights) < 1e-14);
    }

    #[test]
    fn test_l2_error_norm_nonzero() {
        let exact = [2.0, 2.0];
        let approx = [1.0, 1.0];
        let weights = [1.0, 1.0];
        let err = l2_error_norm(&exact, &approx, &weights);
        assert!((err - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_h1_seminorm() {
        let ge = [[1.0f64, 0.0], [1.0, 0.0]];
        let ga = [[0.9, 0.0], [0.9, 0.0]];
        let w = [1.0, 1.0];
        let err = h1_seminorm_error(&ge, &ga, &w);
        assert!((err - (2.0 * 0.01_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_optimal_refinement_ratio() {
        let ratio = optimal_refinement_ratio(0.1, 0.05, 2.0);
        assert!((ratio - 0.5_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_adaptive_fem_driver_stop() {
        let driver = AdaptiveFemDriver::new(1e-6, 10, 0.5, 0.1, true);
        assert!(driver.should_stop(1e-10));
        assert!(!driver.should_stop(1e-3));
    }

    #[test]
    fn test_adaptive_fem_driver_record() {
        let mut driver = AdaptiveFemDriver::new(1e-6, 10, 0.5, 0.1, false);
        driver.record_iteration(100, 0.01);
        assert_eq!(driver.iteration, 1);
        assert_eq!(driver.convergence.history.len(), 1);
    }

    #[test]
    fn test_recovered_stress_von_mises() {
        let mut rs = RecoveredStress::new(2);
        rs.sigma[0] = [100.0, 50.0, 30.0];
        let vm = rs.von_mises(0);
        assert!(vm > 0.0);
    }

    #[test]
    fn test_anisotropic_metric_isotropic() {
        let m = AnisotropicMetric::new(0.1, 0.1, 0.0);
        let mt = m.metric_tensor();
        // Isotropic: diagonal, m00 = m11 = 1/h^2
        assert!((mt[0] - mt[3]).abs() < 1e-10);
        assert!(mt[1].abs() < 1e-10);
    }

    #[test]
    fn test_anisotropic_metric_edge_length() {
        let m = AnisotropicMetric::new(0.1, 0.1, 0.0);
        let len = m.edge_length_metric([0.1, 0.0]);
        assert!((len - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_interpolation_error_estimate() {
        let err = interpolation_error_estimate(0.1, 1, 10.0, 1.0);
        assert!((err - 0.01 * 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_element_size_from_area() {
        // Unit right triangle area = 0.5
        let h = element_size_from_area(0.5);
        assert!((h - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zz_effectivity_index() {
        let mut zz = ZzErrorEstimator::new(3, 2);
        zz.global_error = 0.1;
        let eff = zz.effectivity_index(0.1);
        assert!((eff - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_residual_fraction_exceeding() {
        let mut est = ResidualEstimator::new(4);
        est.eta = vec![0.1, 0.5, 0.2, 0.9];
        let frac = est.fraction_exceeding(0.4);
        assert!((frac - 0.5).abs() < 1e-10);
    }
}
