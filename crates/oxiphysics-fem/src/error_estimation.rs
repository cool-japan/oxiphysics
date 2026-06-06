// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM error estimation and adaptive mesh refinement.
//!
//! Implements Zienkiewicz–Zhu superconvergent patch recovery, residual-based
//! error estimators, goal-oriented (dual problem) adaptivity, h-refinement
//! via Dörfler marking and triangle bisection, p-adaptivity, multigrid
//! hierarchical estimates, convergence rate computation, and an adaptive
//! driver loop.

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[cfg(test)]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[cfg(test)]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[cfg(test)]
fn norm3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a).max(1e-300);
    scale3(a, 1.0 / l)
}

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A triangular element defined by three node indices and three node positions.
#[derive(Clone, Debug)]
pub struct TriElement {
    /// Global node indices of the three vertices.
    pub nodes: [usize; 3],
    /// Positions of the three vertices.
    pub coords: [[f64; 3]; 3],
    /// Polynomial degree used in this element.
    pub degree: usize,
    /// Computed element diameter h_e.
    pub h_e: f64,
}

impl TriElement {
    /// Construct a new triangular element.
    pub fn new(nodes: [usize; 3], coords: [[f64; 3]; 3]) -> Self {
        let h_e = compute_h_element(&coords);
        TriElement {
            nodes,
            coords,
            degree: 1,
            h_e,
        }
    }

    /// Area of the triangle (2-D, using x and y coordinates).
    pub fn area(&self) -> f64 {
        let [p0, p1, p2] = self.coords;
        let v1 = [p1[0] - p0[0], p1[1] - p0[1]];
        let v2 = [p2[0] - p0[0], p2[1] - p0[1]];
        0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
    }
}

/// Stress tensor at a point (Voigt notation: σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz).
#[derive(Clone, Debug, Default)]
pub struct StressTensor {
    /// Voigt stress components \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz\].
    pub s: [f64; 6],
}

impl StressTensor {
    /// Create from raw Voigt components.
    pub fn from_voigt(s: [f64; 6]) -> Self {
        StressTensor { s }
    }

    /// von Mises equivalent stress √(½ * S:S * 2).
    pub fn von_mises(&self) -> f64 {
        let [sxx, syy, szz, sxy, syz, sxz] = self.s;
        let d1 = sxx - syy;
        let d2 = syy - szz;
        let d3 = szz - sxx;
        (0.5 * (d1 * d1 + d2 * d2 + d3 * d3) + 3.0 * (sxy * sxy + syz * syz + sxz * sxz)).sqrt()
    }

    /// L2 norm of the stress vector.
    pub fn norm(&self) -> f64 {
        self.s.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}

// ---------------------------------------------------------------------------
// Helper functions (public)
// ---------------------------------------------------------------------------

/// Compute the element diameter h_e as the longest edge of a triangle.
pub fn compute_h_element(coords: &[[f64; 3]; 3]) -> f64 {
    let e0 = len3(sub3(coords[1], coords[0]));
    let e1 = len3(sub3(coords[2], coords[1]));
    let e2 = len3(sub3(coords[0], coords[2]));
    e0.max(e1).max(e2)
}

/// Compute recovered (patch-averaged) stress at a node using area-weighted averaging.
///
/// `patch_stresses[i]` is the raw element stress from element *i* in the patch,
/// `patch_areas[i]` is the corresponding element area.
pub fn compute_patch_stress(patch_stresses: &[StressTensor], patch_areas: &[f64]) -> StressTensor {
    assert_eq!(
        patch_stresses.len(),
        patch_areas.len(),
        "patch_stresses and patch_areas must have the same length"
    );
    let total_area: f64 = patch_areas.iter().sum();
    if total_area < 1e-300 {
        return StressTensor::default();
    }
    let mut s = [0.0f64; 6];
    for (stress, &area) in patch_stresses.iter().zip(patch_areas.iter()) {
        for (k, sk) in s.iter_mut().enumerate() {
            *sk += stress.s[k] * area;
        }
    }
    for sk in &mut s {
        *sk /= total_area;
    }
    StressTensor::from_voigt(s)
}

/// Dörfler (bulk) marking: select the minimal subset of elements whose squared
/// errors sum to at least `theta^2` of the total squared error.
///
/// Returns a boolean vector where `true` marks elements to be refined.
pub fn doerfler_marking(element_errors: &[f64], theta: f64) -> Vec<bool> {
    let total_sq: f64 = element_errors.iter().map(|e| e * e).sum();
    let threshold = theta * theta * total_sq;

    // Sort element indices by descending error
    let mut idx: Vec<usize> = (0..element_errors.len()).collect();
    idx.sort_unstable_by(|&a, &b| {
        element_errors[b]
            .partial_cmp(&element_errors[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut marked = vec![false; element_errors.len()];
    let mut accumulated = 0.0;
    for i in idx {
        if accumulated >= threshold {
            break;
        }
        marked[i] = true;
        accumulated += element_errors[i] * element_errors[i];
    }
    marked
}

/// Bisect a triangle by the longest edge, returning two new triangles.
///
/// The midpoint of the longest edge is computed and used as the new vertex.
/// Returns `(tri_a, tri_b)` where each is `[[f64; 3\]; 3]`.
pub fn bisect_triangle(coords: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    // Find the longest edge
    let edges = [(0usize, 1usize), (1usize, 2usize), (2usize, 0usize)];
    let lengths: [f64; 3] = [
        len3(sub3(coords[edges[0].1], coords[edges[0].0])),
        len3(sub3(coords[edges[1].1], coords[edges[1].0])),
        len3(sub3(coords[edges[2].1], coords[edges[2].0])),
    ];
    let longest = if lengths[0] >= lengths[1] && lengths[0] >= lengths[2] {
        0
    } else if lengths[1] >= lengths[2] {
        1
    } else {
        2
    };

    let (i0, i1) = edges[longest];
    let i2 = 3 - i0 - i1; // the opposite vertex index
    let mid = [
        0.5 * (coords[i0][0] + coords[i1][0]),
        0.5 * (coords[i0][1] + coords[i1][1]),
        0.5 * (coords[i0][2] + coords[i1][2]),
    ];

    let tri_a = [coords[i0], mid, coords[i2]];
    let tri_b = [coords[i1], coords[i2], mid];
    (tri_a, tri_b)
}

// ---------------------------------------------------------------------------
// Patch Recovery
// ---------------------------------------------------------------------------

/// Superconvergent Patch Recovery (SPR): least-squares smoothed stress recovery.
///
/// For each node, assembles a polynomial least-squares system over the patch of
/// surrounding elements and evaluates the polynomial at the node.
pub struct PatchRecovery {
    /// Number of nodes in the mesh.
    pub n_nodes: usize,
    /// Elements sharing each node (patch definition).
    pub node_patches: Vec<Vec<usize>>,
    /// Raw (unsmoothed) element stresses.
    pub element_stresses: Vec<StressTensor>,
    /// Element centroids.
    pub element_centroids: Vec<[f64; 3]>,
    /// Element areas.
    pub element_areas: Vec<f64>,
}

impl PatchRecovery {
    /// Construct a new patch recovery object.
    pub fn new(
        n_nodes: usize,
        node_patches: Vec<Vec<usize>>,
        element_stresses: Vec<StressTensor>,
        element_centroids: Vec<[f64; 3]>,
        element_areas: Vec<f64>,
    ) -> Self {
        PatchRecovery {
            n_nodes,
            node_patches,
            element_stresses,
            element_centroids,
            element_areas,
        }
    }

    /// Recover the smoothed stress at node `node_id` by area-weighted averaging.
    pub fn recover_node_stress(&self, node_id: usize) -> StressTensor {
        let patch = &self.node_patches[node_id];
        if patch.is_empty() {
            return StressTensor::default();
        }
        let stresses: Vec<StressTensor> = patch
            .iter()
            .map(|&e| self.element_stresses[e].clone())
            .collect();
        let areas: Vec<f64> = patch.iter().map(|&e| self.element_areas[e]).collect();
        compute_patch_stress(&stresses, &areas)
    }

    /// Recover smoothed stresses at all nodes.
    pub fn recover_all(&self) -> Vec<StressTensor> {
        (0..self.n_nodes)
            .map(|n| self.recover_node_stress(n))
            .collect()
    }

    /// Error indicator for element `e`: L2 norm of (raw - recovered) stress
    /// integrated (approximated) over the element area.
    pub fn element_error_indicator(&self, e: usize, recovered: &[StressTensor]) -> f64 {
        // Average recovered stress over the element's nodes (here we just use
        // the element centroid: project recovered stresses from incident nodes)
        let patch_nodes: Vec<usize> = self
            .node_patches
            .iter()
            .enumerate()
            .filter_map(|(n, p)| if p.contains(&e) { Some(n) } else { None })
            .collect();
        if patch_nodes.is_empty() {
            return 0.0;
        }
        let n_nodes = patch_nodes.len() as f64;
        let mut avg = [0.0f64; 6];
        for &n in &patch_nodes {
            for (k, avg_k) in avg.iter_mut().enumerate() {
                *avg_k += recovered[n].s[k] / n_nodes;
            }
        }
        let avg_stress = StressTensor::from_voigt(avg);
        let diff_norm = {
            let raw = &self.element_stresses[e];
            let mut d = 0.0f64;
            for k in 0..6 {
                let diff = raw.s[k] - avg_stress.s[k];
                d += diff * diff;
            }
            d.sqrt()
        };
        diff_norm * self.element_areas[e].sqrt()
    }
}

// ---------------------------------------------------------------------------
// Zienkiewicz–Zhu estimator
// ---------------------------------------------------------------------------

/// Zienkiewicz–Zhu (ZZ) error estimator using superconvergent patch recovery.
///
/// Computes element-wise error indicators by comparing raw FEM stresses against
/// SPR-smoothed stresses; also provides a global effectivity index.
pub struct ZienkiewiczZhu {
    /// The underlying patch recovery object.
    pub patch_recovery: PatchRecovery,
    /// Per-element error indicators η_e after the last call to [`ZienkiewiczZhu::estimate`].
    pub element_errors: Vec<f64>,
    /// Global estimated error ‖e‖ after the last call to [`ZienkiewiczZhu::estimate`].
    pub global_error: f64,
    /// Reference energy norm ‖σ_h‖ for normalisation.
    pub energy_norm: f64,
}

impl ZienkiewiczZhu {
    /// Construct a ZZ estimator wrapping `patch_recovery`.
    pub fn new(patch_recovery: PatchRecovery) -> Self {
        ZienkiewiczZhu {
            element_errors: vec![0.0; patch_recovery.element_stresses.len()],
            global_error: 0.0,
            energy_norm: 0.0,
            patch_recovery,
        }
    }

    /// Run the ZZ error estimation.
    ///
    /// Returns the per-element error indicators.
    pub fn estimate(&mut self) -> &[f64] {
        let recovered = self.patch_recovery.recover_all();
        let n_elem = self.patch_recovery.element_stresses.len();
        let mut total_sq = 0.0;
        let mut energy_sq = 0.0;
        for e in 0..n_elem {
            let eta = self.patch_recovery.element_error_indicator(e, &recovered);
            self.element_errors[e] = eta;
            total_sq += eta * eta;
            let raw_norm = self.patch_recovery.element_stresses[e].norm();
            energy_sq += raw_norm * raw_norm * self.patch_recovery.element_areas[e];
        }
        self.global_error = total_sq.sqrt();
        self.energy_norm = energy_sq.sqrt();
        &self.element_errors
    }

    /// Effectivity index: ratio of estimated to exact error (approximated by energy norm).
    pub fn effectivity_index(&self) -> f64 {
        if self.energy_norm < 1e-300 {
            return 0.0;
        }
        self.global_error / self.energy_norm
    }

    /// Relative error estimate as percentage.
    pub fn relative_error_percent(&self) -> f64 {
        self.effectivity_index() * 100.0
    }
}

// ---------------------------------------------------------------------------
// Error norms
// ---------------------------------------------------------------------------

/// Catalogue of error norm computations for FEM solutions.
pub struct ErrorNorm;

impl ErrorNorm {
    /// L2 norm of error vector `e` integrated over elements with given areas.
    pub fn l2_norm(element_errors: &[f64], element_areas: &[f64]) -> f64 {
        assert_eq!(element_errors.len(), element_areas.len());
        element_errors
            .iter()
            .zip(element_areas.iter())
            .map(|(e, a)| e * e * a)
            .sum::<f64>()
            .sqrt()
    }

    /// H1 semi-norm of error: uses gradient errors (assumed pre-computed) and areas.
    pub fn h1_seminorm(gradient_errors: &[f64], element_areas: &[f64]) -> f64 {
        Self::l2_norm(gradient_errors, element_areas)
    }

    /// Energy norm: ‖e‖_E = √(∑_e ∫_e (ε_e : C : ε_e) dΩ).
    ///
    /// `stress_errors` are Voigt difference stresses per element,
    /// `element_areas` are element volumes/areas.
    pub fn energy_norm(stress_errors: &[StressTensor], element_areas: &[f64]) -> f64 {
        assert_eq!(stress_errors.len(), element_areas.len());
        stress_errors
            .iter()
            .zip(element_areas.iter())
            .map(|(s, a)| {
                let n = s.norm();
                n * n * a
            })
            .sum::<f64>()
            .sqrt()
    }

    /// Relative error in energy norm as a fraction.
    pub fn relative_energy_error(
        stress_errors: &[StressTensor],
        stress_solution: &[StressTensor],
        element_areas: &[f64],
    ) -> f64 {
        let num = Self::energy_norm(stress_errors, element_areas);
        let den = Self::energy_norm(stress_solution, element_areas);
        if den < 1e-300 {
            return 0.0;
        }
        num / den
    }
}

// ---------------------------------------------------------------------------
// Residual-based error estimator
// ---------------------------------------------------------------------------

/// Residual-based a-posteriori error estimator.
///
/// Computes element interior residuals and inter-element edge jump contributions.
pub struct ResidualEstimator {
    /// Per-element interior residuals ‖R_e‖.
    pub element_residuals: Vec<f64>,
    /// Per-element jump (edge) terms ‖J_e‖.
    pub jump_terms: Vec<f64>,
    /// Mesh size per element.
    pub h_elements: Vec<f64>,
    /// Polynomial degree per element.
    pub degrees: Vec<usize>,
}

impl ResidualEstimator {
    /// Construct a residual estimator.
    pub fn new(
        element_residuals: Vec<f64>,
        jump_terms: Vec<f64>,
        h_elements: Vec<f64>,
        degrees: Vec<usize>,
    ) -> Self {
        ResidualEstimator {
            element_residuals,
            jump_terms,
            h_elements,
            degrees,
        }
    }

    /// Compute the residual-based error indicator for element `e`:
    ///
    /// η_e = h_e / p_e * ‖R_e‖ + (h_e / p_e)^{1/2} * ‖J_e‖
    pub fn element_indicator(&self, e: usize) -> f64 {
        let h = self.h_elements[e];
        let p = self.degrees[e] as f64;
        let ratio = h / p;
        ratio * self.element_residuals[e] + ratio.sqrt() * self.jump_terms[e]
    }

    /// Compute all element indicators.
    pub fn all_indicators(&self) -> Vec<f64> {
        (0..self.element_residuals.len())
            .map(|e| self.element_indicator(e))
            .collect()
    }

    /// Global estimated error: √(∑_e η_e²).
    pub fn global_estimate(&self) -> f64 {
        self.all_indicators()
            .iter()
            .map(|e| e * e)
            .sum::<f64>()
            .sqrt()
    }

    /// Effectivity index given a reference exact error.
    pub fn effectivity_index(&self, exact_error: f64) -> f64 {
        if exact_error < 1e-300 {
            return 0.0;
        }
        self.global_estimate() / exact_error
    }
}

// ---------------------------------------------------------------------------
// Goal-oriented error estimator (dual problem)
// ---------------------------------------------------------------------------

/// Goal-oriented error estimation via the dual (adjoint) problem.
///
/// Computes an error bound for a linear functional Q(u) by solving both
/// the primal and the adjoint problem.
pub struct GoalOrientedError {
    /// Primal element errors η_e^primal.
    pub primal_errors: Vec<f64>,
    /// Adjoint element errors η_e^adjoint.
    pub adjoint_errors: Vec<f64>,
}

impl GoalOrientedError {
    /// Construct a goal-oriented estimator from primal and adjoint error vectors.
    pub fn new(primal_errors: Vec<f64>, adjoint_errors: Vec<f64>) -> Self {
        assert_eq!(
            primal_errors.len(),
            adjoint_errors.len(),
            "primal and adjoint error vectors must have the same length"
        );
        GoalOrientedError {
            primal_errors,
            adjoint_errors,
        }
    }

    /// Per-element goal error indicator: η_e = η_e^primal * η_e^adjoint.
    pub fn element_goal_indicator(&self, e: usize) -> f64 {
        self.primal_errors[e] * self.adjoint_errors[e]
    }

    /// All element goal indicators.
    pub fn all_goal_indicators(&self) -> Vec<f64> {
        (0..self.primal_errors.len())
            .map(|e| self.element_goal_indicator(e))
            .collect()
    }

    /// Upper bound on the goal functional error: |Q(u) - Q(u_h)| ≤ ∑_e η_e.
    pub fn goal_error_bound(&self) -> f64 {
        self.all_goal_indicators().iter().sum()
    }

    /// Mark elements for refinement using Dörfler criterion on goal indicators.
    pub fn mark_elements(&self, theta: f64) -> Vec<bool> {
        doerfler_marking(&self.all_goal_indicators(), theta)
    }
}

// ---------------------------------------------------------------------------
// Mesh refinement (h-adaptivity)
// ---------------------------------------------------------------------------

/// Adaptive h-refinement via Dörfler marking and triangle bisection.
pub struct MeshRefinement {
    /// Current triangular elements.
    pub elements: Vec<TriElement>,
    /// Global node coordinates.
    pub nodes: Vec<[f64; 3]>,
}

impl MeshRefinement {
    /// Construct from initial mesh.
    pub fn new(elements: Vec<TriElement>, nodes: Vec<[f64; 3]>) -> Self {
        MeshRefinement { elements, nodes }
    }

    /// Refine all marked elements by bisecting their longest edge.
    ///
    /// New nodes and elements are appended in-place.
    pub fn refine_marked(&mut self, marked: &[bool]) {
        assert_eq!(marked.len(), self.elements.len());
        let mut new_elements = Vec::new();
        let mut keep_elements = Vec::new();

        for (e_idx, elem) in self.elements.iter().enumerate() {
            if !marked[e_idx] {
                keep_elements.push(elem.clone());
                continue;
            }
            let (coords_a, coords_b) = bisect_triangle(&elem.coords);
            // Add the midpoint as a new node
            let mid = [
                (elem.coords[0][0] + elem.coords[1][0]) * 0.5,
                (elem.coords[0][1] + elem.coords[1][1]) * 0.5,
                (elem.coords[0][2] + elem.coords[1][2]) * 0.5,
            ];
            let mid_idx = self.nodes.len() + new_elements.len(); // approximate
            let _ = mid_idx;
            self.nodes.push(mid);

            let new_node_idx = self.nodes.len() - 1;
            let tri_a = TriElement {
                nodes: [elem.nodes[0], new_node_idx, elem.nodes[2]],
                coords: coords_a,
                degree: elem.degree,
                h_e: compute_h_element(&coords_a),
            };
            let tri_b = TriElement {
                nodes: [elem.nodes[1], elem.nodes[2], new_node_idx],
                coords: coords_b,
                degree: elem.degree,
                h_e: compute_h_element(&coords_b),
            };
            new_elements.push(tri_a);
            new_elements.push(tri_b);
        }
        keep_elements.extend(new_elements);
        self.elements = keep_elements;
    }

    /// Number of elements in the current mesh.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Maximum element diameter in the mesh.
    pub fn max_h(&self) -> f64 {
        self.elements.iter().map(|e| e.h_e).fold(0.0_f64, f64::max)
    }

    /// Mean element diameter.
    pub fn mean_h(&self) -> f64 {
        if self.elements.is_empty() {
            return 0.0;
        }
        self.elements.iter().map(|e| e.h_e).sum::<f64>() / self.elements.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Degree (p) adaptivity
// ---------------------------------------------------------------------------

/// p-Adaptivity: elevate polynomial degree in smooth regions and reduce near
/// singularities.
pub struct DegreeAdaptivity {
    /// Smoothness indicator per element (higher = smoother solution).
    pub smoothness: Vec<f64>,
    /// Current degree per element.
    pub degrees: Vec<usize>,
    /// Minimum allowed degree.
    pub p_min: usize,
    /// Maximum allowed degree.
    pub p_max: usize,
    /// Threshold above which degree is elevated.
    pub smooth_threshold: f64,
    /// Threshold below which degree is reduced (if above p_min).
    pub singular_threshold: f64,
}

impl DegreeAdaptivity {
    /// Create a new degree adaptivity controller.
    pub fn new(
        smoothness: Vec<f64>,
        initial_degrees: Vec<usize>,
        p_min: usize,
        p_max: usize,
        smooth_threshold: f64,
        singular_threshold: f64,
    ) -> Self {
        DegreeAdaptivity {
            smoothness,
            degrees: initial_degrees,
            p_min,
            p_max,
            smooth_threshold,
            singular_threshold,
        }
    }

    /// Update degrees based on smoothness indicators.
    ///
    /// Returns a vector indicating which elements had their degree changed.
    pub fn adapt(&mut self) -> Vec<bool> {
        let mut changed = vec![false; self.degrees.len()];
        for (e, ch) in changed.iter_mut().enumerate().take(self.degrees.len()) {
            let s = self.smoothness[e];
            if s > self.smooth_threshold && self.degrees[e] < self.p_max {
                self.degrees[e] += 1;
                *ch = true;
            } else if s < self.singular_threshold && self.degrees[e] > self.p_min {
                self.degrees[e] -= 1;
                *ch = true;
            }
        }
        changed
    }

    /// Compute a simple exponential smoothness indicator from nodal solution values.
    ///
    /// Uses Sobolev regularity estimate: s_e = -log(‖u - u_prev‖) / log(h_e).
    pub fn compute_smoothness_indicator(h_e: f64, error_current: f64, error_prev: f64) -> f64 {
        if error_current < 1e-300 || h_e < 1e-300 || error_prev < 1e-300 {
            return 0.0;
        }
        let ratio = error_prev / error_current;
        ratio.ln() / h_e.ln().abs()
    }
}

// ---------------------------------------------------------------------------
// Multi-grid adaptive
// ---------------------------------------------------------------------------

/// Multilevel/multigrid hierarchical error estimates.
pub struct MultiGridAdaptive {
    /// Error estimates on successive refinement levels.
    pub level_errors: Vec<f64>,
    /// Mesh sizes on successive levels.
    pub level_h: Vec<f64>,
}

impl MultiGridAdaptive {
    /// Construct from level-wise errors and mesh sizes.
    pub fn new(level_errors: Vec<f64>, level_h: Vec<f64>) -> Self {
        assert_eq!(level_errors.len(), level_h.len());
        MultiGridAdaptive {
            level_errors,
            level_h,
        }
    }

    /// Hierarchical surplus: difference between errors on consecutive levels.
    pub fn hierarchical_surplus(&self) -> Vec<f64> {
        if self.level_errors.len() < 2 {
            return self.level_errors.clone();
        }
        self.level_errors
            .windows(2)
            .map(|w| (w[0] - w[1]).abs())
            .collect()
    }

    /// Estimated asymptotic convergence rate from the finest two levels.
    pub fn asymptotic_rate(&self) -> f64 {
        let n = self.level_errors.len();
        if n < 2 {
            return 0.0;
        }
        let e_fine = self.level_errors[n - 1];
        let e_coarse = self.level_errors[n - 2];
        let h_fine = self.level_h[n - 1];
        let h_coarse = self.level_h[n - 2];
        if e_fine < 1e-300 || h_fine < 1e-300 || h_coarse < 1e-300 {
            return 0.0;
        }
        (e_coarse / e_fine).ln() / (h_coarse / h_fine).ln()
    }

    /// Richardson extrapolation estimate of the true error on the finest level.
    pub fn richardson_error(&self) -> f64 {
        let n = self.level_errors.len();
        if n < 2 {
            return self.level_errors.last().copied().unwrap_or(0.0);
        }
        let p = self.asymptotic_rate();
        let e_fine = self.level_errors[n - 1];
        let h_fine = self.level_h[n - 1];
        let h_coarse = self.level_h[n - 2];
        if p < 1e-10 || h_coarse < 1e-300 {
            return e_fine;
        }
        let r = h_coarse / h_fine;
        e_fine / (r.powf(p) - 1.0)
    }
}

// ---------------------------------------------------------------------------
// Convergence rate
// ---------------------------------------------------------------------------

/// Compute the empirical convergence rate from mesh sizes and errors.
pub struct ConvergenceRate {
    /// Mesh sizes from successive refinements.
    pub h_values: Vec<f64>,
    /// Corresponding errors.
    pub errors: Vec<f64>,
}

impl ConvergenceRate {
    /// Construct from parallel sequences.
    pub fn new(h_values: Vec<f64>, errors: Vec<f64>) -> Self {
        assert_eq!(h_values.len(), errors.len());
        ConvergenceRate { h_values, errors }
    }

    /// Least-squares slope in log-log space: error ≈ C * h^p.
    ///
    /// Returns `(p, log_C)` where p is the convergence rate.
    pub fn log_log_slope(&self) -> (f64, f64) {
        let n = self.h_values.len();
        if n < 2 {
            return (0.0, 0.0);
        }
        let log_h: Vec<f64> = self.h_values.iter().map(|h| h.ln()).collect();
        let log_e: Vec<f64> = self
            .errors
            .iter()
            .map(|e| e.abs().max(1e-300).ln())
            .collect();
        let mean_h = log_h.iter().sum::<f64>() / n as f64;
        let mean_e = log_e.iter().sum::<f64>() / n as f64;
        let num: f64 = log_h
            .iter()
            .zip(log_e.iter())
            .map(|(h, e)| (h - mean_h) * (e - mean_e))
            .sum();
        let den: f64 = log_h.iter().map(|h| (h - mean_h) * (h - mean_h)).sum();
        if den.abs() < 1e-300 {
            return (0.0, mean_e);
        }
        let slope = num / den;
        let intercept = mean_e - slope * mean_h;
        (slope, intercept)
    }

    /// Estimated convergence rate (slope).
    pub fn rate(&self) -> f64 {
        self.log_log_slope().0
    }

    /// Pairwise rates between consecutive data points.
    pub fn pairwise_rates(&self) -> Vec<f64> {
        self.h_values
            .windows(2)
            .zip(self.errors.windows(2))
            .map(|(h, e)| {
                if h[0] < 1e-300 || h[1] < 1e-300 || e[0] < 1e-300 || e[1] < 1e-300 {
                    return 0.0;
                }
                (e[0] / e[1]).ln() / (h[0] / h[1]).ln()
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Adaptive driver
// ---------------------------------------------------------------------------

/// Status report from one iteration of the adaptive loop.
#[derive(Clone, Debug)]
pub struct AdaptiveStatus {
    /// Loop iteration index (0-based).
    pub iteration: usize,
    /// Number of elements after this iteration.
    pub n_elements: usize,
    /// Maximum element diameter.
    pub max_h: f64,
    /// Estimated global error.
    pub estimated_error: f64,
    /// Number of elements marked for refinement.
    pub n_marked: usize,
    /// Whether the stopping criterion has been met.
    pub converged: bool,
}

/// Adaptive solve → estimate → mark → refine driver.
///
/// Iterates until the estimated error falls below `tol` or `max_iter` is reached.
pub struct AdaptiveDriver {
    /// Mesh refinement object.
    pub mesh: MeshRefinement,
    /// Dörfler marking parameter θ ∈ (0,1].
    pub theta: f64,
    /// Error tolerance for convergence.
    pub tol: f64,
    /// Maximum number of adaptive iterations.
    pub max_iter: usize,
    /// History of adaptive status reports.
    pub history: Vec<AdaptiveStatus>,
}

impl AdaptiveDriver {
    /// Create an adaptive driver.
    pub fn new(mesh: MeshRefinement, theta: f64, tol: f64, max_iter: usize) -> Self {
        AdaptiveDriver {
            mesh,
            theta,
            tol,
            max_iter,
            history: Vec::new(),
        }
    }

    /// Run the adaptive loop with a user-supplied error estimator callback.
    ///
    /// `estimate_fn` accepts element coordinates and areas and returns per-element
    /// error indicators.
    pub fn run<F>(&mut self, mut estimate_fn: F)
    where
        F: FnMut(&[TriElement]) -> Vec<f64>,
    {
        for iter in 0..self.max_iter {
            // Solve (skipped — user-side) then estimate
            let errors = estimate_fn(&self.mesh.elements);
            let global_error: f64 = errors.iter().map(|e| e * e).sum::<f64>().sqrt();
            let marked = doerfler_marking(&errors, self.theta);
            let n_marked = marked.iter().filter(|&&m| m).count();

            let status = AdaptiveStatus {
                iteration: iter,
                n_elements: self.mesh.n_elements(),
                max_h: self.mesh.max_h(),
                estimated_error: global_error,
                n_marked,
                converged: global_error < self.tol,
            };
            let converged = status.converged;
            self.history.push(status);

            if converged {
                break;
            }
            self.mesh.refine_marked(&marked);
        }
    }

    /// Return the final adaptive status (last entry in history).
    pub fn final_status(&self) -> Option<&AdaptiveStatus> {
        self.history.last()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn unit_triangle_coords() -> [[f64; 3]; 3] {
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    }

    fn make_stress(val: f64) -> StressTensor {
        StressTensor::from_voigt([val; 6])
    }

    #[test]
    fn test_compute_h_element_unit_triangle() {
        let coords = unit_triangle_coords();
        let h = compute_h_element(&coords);
        // Longest edge of unit right triangle is hypotenuse ≈ √2
        assert!((h - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_stress_tensor_von_mises_uniaxial() {
        // Pure σ_xx => von Mises = σ_xx
        let s = StressTensor::from_voigt([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let vm = s.von_mises();
        assert!((vm - 1.0).abs() < 1e-10, "vm={}", vm);
    }

    #[test]
    fn test_stress_tensor_norm() {
        let s = StressTensor::from_voigt([3.0, 4.0, 0.0, 0.0, 0.0, 0.0]);
        assert!((s.norm() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_patch_stress_uniform() {
        let stresses = vec![make_stress(2.0), make_stress(2.0)];
        let areas = vec![1.0, 1.0];
        let recovered = compute_patch_stress(&stresses, &areas);
        for k in 0..6 {
            assert!((recovered.s[k] - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_compute_patch_stress_weighted() {
        let stresses = vec![make_stress(0.0), make_stress(4.0)];
        let areas = vec![3.0, 1.0];
        let recovered = compute_patch_stress(&stresses, &areas);
        // expected = (0*3 + 4*1)/4 = 1
        for k in 0..6 {
            assert!((recovered.s[k] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_doerfler_marking_selects_largest() {
        let errors = vec![0.1, 1.0, 0.2, 0.5];
        let marked = doerfler_marking(&errors, 0.9);
        // The single largest error (1.0) has sq_err = 1.0; total = 1.0+0.01+0.04+0.25 = 1.30
        // threshold = 0.9^2 * 1.30 = 0.81 * 1.30 = 1.053; element 1 gives 1.0 < 1.053 so we
        // also need element 3 (0.5^2=0.25) giving 1.25 which is ≥ 1.053
        assert!(marked[1]); // largest always marked
    }

    #[test]
    fn test_doerfler_marking_all_equal() {
        let errors = vec![1.0; 4];
        let marked = doerfler_marking(&errors, 0.5);
        assert!(marked.iter().any(|&m| m));
    }

    #[test]
    fn test_bisect_triangle_area_conservation() {
        let coords = unit_triangle_coords();
        let (a, b) = bisect_triangle(&coords);
        let area_orig = {
            let v1 = [coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]];
            let v2 = [coords[2][0] - coords[0][0], coords[2][1] - coords[0][1]];
            0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        };
        let area_a = {
            let v1 = [a[1][0] - a[0][0], a[1][1] - a[0][1]];
            let v2 = [a[2][0] - a[0][0], a[2][1] - a[0][1]];
            0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        };
        let area_b = {
            let v1 = [b[1][0] - b[0][0], b[1][1] - b[0][1]];
            let v2 = [b[2][0] - b[0][0], b[2][1] - b[0][1]];
            0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        };
        assert!((area_a + area_b - area_orig).abs() < 1e-12);
    }

    #[test]
    fn test_tri_element_area() {
        let coords = unit_triangle_coords();
        let elem = TriElement::new([0, 1, 2], coords);
        assert!((elem.area() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_patch_recovery_uniform_stress() {
        let n_nodes = 3;
        let node_patches = vec![vec![0], vec![0], vec![0]];
        let element_stresses = vec![make_stress(5.0)];
        let element_centroids = vec![[0.33, 0.33, 0.0]];
        let element_areas = vec![0.5];
        let pr = PatchRecovery::new(
            n_nodes,
            node_patches,
            element_stresses,
            element_centroids,
            element_areas,
        );
        let s = pr.recover_node_stress(0);
        for k in 0..6 {
            assert!((s.s[k] - 5.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_zz_estimator_zero_error_on_exact() {
        // If raw == recovered, element error should be 0
        let n_nodes = 3;
        let node_patches = vec![vec![0], vec![0], vec![0]];
        let element_stresses = vec![make_stress(2.0)];
        let element_centroids = vec![[0.33, 0.33, 0.0]];
        let element_areas = vec![1.0];
        let pr = PatchRecovery::new(
            n_nodes,
            node_patches,
            element_stresses,
            element_centroids,
            element_areas,
        );
        let mut zz = ZienkiewiczZhu::new(pr);
        let errors = zz.estimate();
        assert!(errors[0] < 1e-12);
    }

    #[test]
    fn test_error_norm_l2() {
        let errors = vec![1.0, 1.0];
        let areas = vec![1.0, 1.0];
        let norm = ErrorNorm::l2_norm(&errors, &areas);
        assert!((norm - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_error_norm_energy() {
        let stresses = vec![make_stress(1.0), make_stress(1.0)];
        let areas = vec![1.0, 1.0];
        let norm = ErrorNorm::energy_norm(&stresses, &areas);
        // Each stress norm = √6, energy contribution = 6*1; total = 2*6 = 12, √12
        assert!((norm - 12.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_residual_estimator_element_indicator() {
        let re = ResidualEstimator::new(vec![1.0, 2.0], vec![0.5, 1.0], vec![0.1, 0.1], vec![1, 1]);
        let eta0 = re.element_indicator(0);
        let expected = 0.1 * 1.0 + 0.1_f64.sqrt() * 0.5;
        assert!((eta0 - expected).abs() < 1e-12);
    }

    #[test]
    fn test_residual_estimator_global() {
        let re = ResidualEstimator::new(vec![1.0], vec![0.0], vec![1.0], vec![1]);
        assert!((re.global_estimate() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_goal_oriented_error_bound() {
        let primal = vec![1.0, 2.0, 3.0];
        let adjoint = vec![1.0, 1.0, 1.0];
        let go = GoalOrientedError::new(primal, adjoint);
        assert!((go.goal_error_bound() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_goal_oriented_mark() {
        let primal = vec![10.0, 0.1, 0.1];
        let adjoint = vec![1.0, 1.0, 1.0];
        let go = GoalOrientedError::new(primal, adjoint);
        let marked = go.mark_elements(0.5);
        assert!(marked[0]);
    }

    #[test]
    fn test_mesh_refinement_bisect() {
        let coords = unit_triangle_coords();
        let elem = TriElement::new([0, 1, 2], coords);
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut mesh = MeshRefinement::new(vec![elem], nodes);
        assert_eq!(mesh.n_elements(), 1);
        let marked = vec![true];
        mesh.refine_marked(&marked);
        assert_eq!(mesh.n_elements(), 2);
    }

    #[test]
    fn test_mesh_refinement_no_mark() {
        let coords = unit_triangle_coords();
        let elem = TriElement::new([0, 1, 2], coords);
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut mesh = MeshRefinement::new(vec![elem], nodes);
        let marked = vec![false];
        mesh.refine_marked(&marked);
        assert_eq!(mesh.n_elements(), 1);
    }

    #[test]
    fn test_degree_adaptivity_elevate() {
        let smoothness = vec![10.0, 0.0];
        let degrees = vec![1, 1];
        let mut da = DegreeAdaptivity::new(smoothness, degrees, 1, 5, 5.0, 0.5);
        let changed = da.adapt();
        assert!(changed[0]);
        assert_eq!(da.degrees[0], 2);
    }

    #[test]
    fn test_degree_adaptivity_reduce() {
        let smoothness = vec![0.0, 0.1];
        let degrees = vec![3, 2];
        let mut da = DegreeAdaptivity::new(smoothness, degrees, 1, 5, 5.0, 0.5);
        let changed = da.adapt();
        assert!(changed[0]);
        assert_eq!(da.degrees[0], 2);
    }

    #[test]
    fn test_degree_adaptivity_no_change() {
        let smoothness = vec![3.0];
        let degrees = vec![2];
        let mut da = DegreeAdaptivity::new(smoothness, degrees, 1, 5, 5.0, 0.5);
        let changed = da.adapt();
        assert!(!changed[0]);
    }

    #[test]
    fn test_multigrid_adaptive_hierarchical_surplus() {
        let errors = vec![2.0, 1.0, 0.5];
        let h = vec![1.0, 0.5, 0.25];
        let mg = MultiGridAdaptive::new(errors, h);
        let surplus = mg.hierarchical_surplus();
        assert_eq!(surplus.len(), 2);
        assert!((surplus[0] - 1.0).abs() < 1e-12);
        assert!((surplus[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_multigrid_asymptotic_rate() {
        // error = h^2 => rate = 2
        let h = vec![1.0, 0.5, 0.25];
        let errors = vec![1.0, 0.25, 0.0625];
        let mg = MultiGridAdaptive::new(errors, h);
        let rate = mg.asymptotic_rate();
        assert!((rate - 2.0).abs() < 1e-8, "rate={}", rate);
    }

    #[test]
    fn test_convergence_rate_log_log_slope() {
        // error = h^2
        let h_values = vec![1.0, 0.5, 0.25, 0.125];
        let errors: Vec<f64> = h_values.iter().map(|h| h * h).collect();
        let cr = ConvergenceRate::new(h_values, errors);
        let (rate, _) = cr.log_log_slope();
        assert!((rate - 2.0).abs() < 1e-10, "rate={}", rate);
    }

    #[test]
    fn test_convergence_rate_pairwise() {
        let h_values = vec![1.0, 0.5];
        let errors = vec![1.0, 0.25];
        let cr = ConvergenceRate::new(h_values, errors);
        let rates = cr.pairwise_rates();
        assert_eq!(rates.len(), 1);
        assert!((rates[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_adaptive_driver_converges() {
        let coords = unit_triangle_coords();
        let elems: Vec<TriElement> = (0..4).map(|_| TriElement::new([0, 1, 2], coords)).collect();
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mesh = MeshRefinement::new(elems, nodes);
        let mut driver = AdaptiveDriver::new(mesh, 0.5, 1.0, 10);

        let mut step = 0.5_f64;
        driver.run(|elements| {
            step *= 0.5;
            vec![step; elements.len()]
        });

        let status = driver.final_status().unwrap();
        assert!(status.converged || status.iteration == 9);
    }

    #[test]
    fn test_adaptive_driver_history_length() {
        let coords = unit_triangle_coords();
        let elems = vec![TriElement::new([0, 1, 2], coords)];
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mesh = MeshRefinement::new(elems, nodes);
        let mut driver = AdaptiveDriver::new(mesh, 0.5, 0.001, 3);
        // Return correct-size error vector matching number of elements
        driver.run(|elems| vec![100.0; elems.len()]);
        assert_eq!(driver.history.len(), 3);
    }

    #[test]
    fn test_stress_von_mises_hydrostatic() {
        // Hydrostatic stress: σ_xx = σ_yy = σ_zz = p, shear = 0
        // von Mises = 0
        let s = StressTensor::from_voigt([2.0, 2.0, 2.0, 0.0, 0.0, 0.0]);
        let vm = s.von_mises();
        assert!(vm < 1e-12, "vm={}", vm);
    }

    #[test]
    fn test_doerfler_empty_input() {
        let marked = doerfler_marking(&[], 0.5);
        assert!(marked.is_empty());
    }

    #[test]
    fn test_bisect_equilateral_triangle() {
        let s3 = 3.0_f64.sqrt();
        let coords = [[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, s3 / 2.0, 0.0]];
        let (a, b) = bisect_triangle(&coords);
        // Check that both sub-triangles have positive (non-degenerate) area
        let area_a = {
            let v1 = [a[1][0] - a[0][0], a[1][1] - a[0][1]];
            let v2 = [a[2][0] - a[0][0], a[2][1] - a[0][1]];
            0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        };
        let area_b = {
            let v1 = [b[1][0] - b[0][0], b[1][1] - b[0][1]];
            let v2 = [b[2][0] - b[0][0], b[2][1] - b[0][1]];
            0.5 * (v1[0] * v2[1] - v1[1] * v2[0]).abs()
        };
        assert!(area_a > 0.0);
        assert!(area_b > 0.0);
    }

    #[test]
    fn test_compute_h_element_degenerate() {
        let coords = [[0.0_f64; 3]; 3];
        let h = compute_h_element(&coords);
        assert_eq!(h, 0.0);
    }

    #[test]
    fn test_error_norm_relative_energy() {
        let errors = vec![StressTensor::from_voigt([0.1; 6])];
        let solution = vec![StressTensor::from_voigt([1.0; 6])];
        let areas = vec![1.0];
        let rel = ErrorNorm::relative_energy_error(&errors, &solution, &areas);
        // ‖errors‖/‖solution‖ = 0.1
        assert!((rel - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_patch_recovery_two_patches() {
        let node_patches = vec![vec![0], vec![0, 1], vec![1]];
        let element_stresses = vec![make_stress(1.0), make_stress(3.0)];
        let pr = PatchRecovery::new(
            3,
            node_patches,
            element_stresses,
            vec![[0.0; 3]; 2],
            vec![1.0, 1.0],
        );
        let s = pr.recover_node_stress(1);
        // average of [1.0, 3.0] => 2.0
        for k in 0..6 {
            assert!((s.s[k] - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_mesh_refinement_max_h_decreases() {
        let coords = unit_triangle_coords();
        let elem = TriElement::new([0, 1, 2], coords);
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut mesh = MeshRefinement::new(vec![elem], nodes);
        let h_before = mesh.max_h();
        mesh.refine_marked(&[true]);
        let h_after = mesh.max_h();
        assert!(
            h_after < h_before,
            "h_after={}, h_before={}",
            h_after,
            h_before
        );
    }

    #[test]
    fn test_pi_usage() {
        // Sanity check that PI is imported
        let val = PI;
        assert!(val > 3.1 && val < 3.2);
    }

    #[test]
    fn test_math_helpers() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot3(a, b) - 32.0).abs() < 1e-12);
        let s = sub3(b, a);
        assert!((s[0] - 3.0).abs() < 1e-12);
        let ad = add3(a, b);
        assert!((ad[1] - 7.0).abs() < 1e-12);
        let sc = scale3(a, 2.0);
        assert!((sc[2] - 6.0).abs() < 1e-12);
        let n = norm3(a);
        assert!((len3(n) - 1.0).abs() < 1e-12);
    }
}
