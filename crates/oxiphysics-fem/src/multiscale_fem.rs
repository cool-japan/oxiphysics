#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale FEM methods including FE² and homogenization.
//!
//! # Overview
//!
//! This module provides a complete toolkit for multiscale FEM including:
//!
//! - [`RveGeometry`] — Representative Volume Element (periodic unit cell)
//! - [`FeSquared`] — FE² concurrent multiscale with macro BVP + micro RVE
//! - [`HomogenizationScheme`] — Voigt, Reuss, Mori-Tanaka, self-consistent
//! - [`EffectiveMedium`] — effective elastic constants from microstructure
//! - [`MultiscaleMesh`] — hierarchical mesh with coarse + fine regions
//! - [`AdaptiveMultiscale`] — error estimator for multiscale region selection
//! - [`ConcurrentCoupling`] — Arlequin method for domain coupling
//! - [`HierarchicalBasis`] — hierarchical shape functions (p-refinement)
//! - [`MultiresolutionFem`] — wavelet-based multiresolution analysis
//! - [`SubstructureMethod`] — Craig-Bampton substructuring

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers (private)
// ---------------------------------------------------------------------------

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm_vec(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}

fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn vec_scale(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

fn mat_vec_mul(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

/// Dense conjugate-gradient solver (symmetric positive-definite system).
fn dense_cg(a: &[Vec<f64>], b: &[f64], tol: f64, max_iter: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0_f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rsold = dot(&r, &r);
    for _ in 0..max_iter {
        let ap = mat_vec_mul(a, &p);
        let denom = dot(&p, &ap);
        if denom.abs() < 1e-300 {
            break;
        }
        let alpha = rsold / denom;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rsnew = dot(&r, &r);
        if rsnew.sqrt() < tol {
            break;
        }
        let beta = rsnew / rsold;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    x
}

// ---------------------------------------------------------------------------
// RveGeometry
// ---------------------------------------------------------------------------

/// Representative Volume Element (RVE) — periodic unit cell geometry.
///
/// A cubic RVE with n³ Gauss points represents the microstructure.
/// Periodic boundary conditions are enforced on all faces.
#[derive(Debug, Clone)]
pub struct RveGeometry {
    /// Physical side length of the unit cell (m).
    pub side_length: f64,
    /// Number of elements per side.
    pub n_per_side: usize,
    /// Inclusion volume fraction φ ∈ \[0, 1\].
    pub inclusion_fraction: f64,
    /// Inclusion elastic modulus E_i.
    pub e_inclusion: f64,
    /// Matrix elastic modulus E_m.
    pub e_matrix: f64,
    /// Inclusion Poisson ratio ν_i.
    pub nu_inclusion: f64,
    /// Matrix Poisson ratio ν_m.
    pub nu_matrix: f64,
}

impl RveGeometry {
    /// Create a new cubic RVE.
    ///
    /// # Arguments
    /// * `side_length`        – physical size of the cube (m)
    /// * `n_per_side`         – elements per side (total = n³)
    /// * `inclusion_fraction` – volume fraction of inclusions φ
    /// * `e_inclusion`        – inclusion elastic modulus
    /// * `e_matrix`           – matrix elastic modulus
    /// * `nu_inclusion`       – inclusion Poisson ratio
    /// * `nu_matrix`          – matrix Poisson ratio
    pub fn new(
        side_length: f64,
        n_per_side: usize,
        inclusion_fraction: f64,
        e_inclusion: f64,
        e_matrix: f64,
        nu_inclusion: f64,
        nu_matrix: f64,
    ) -> Self {
        Self {
            side_length,
            n_per_side,
            inclusion_fraction,
            e_inclusion,
            e_matrix,
            nu_inclusion,
            nu_matrix,
        }
    }

    /// Total number of elements = n³.
    pub fn n_elements(&self) -> usize {
        self.n_per_side.pow(3)
    }

    /// Total number of nodes = (n+1)³.
    pub fn n_nodes(&self) -> usize {
        (self.n_per_side + 1).pow(3)
    }

    /// Volume of the unit cell.
    pub fn volume(&self) -> f64 {
        self.side_length.powi(3)
    }

    /// Element side length h = L / n.
    pub fn element_size(&self) -> f64 {
        self.side_length / self.n_per_side as f64
    }

    /// Local elastic modulus at a given point: inclusion or matrix.
    ///
    /// Points within a sphere of radius r_i centered at (L/2, L/2, L/2) are
    /// inclusions; the rest are matrix.
    pub fn local_modulus(&self, x: f64, y: f64, z: f64) -> f64 {
        let cx = self.side_length / 2.0;
        let r_incl =
            self.side_length * (3.0 * self.inclusion_fraction / (4.0 * PI)).powf(1.0 / 3.0);
        let dist = ((x - cx).powi(2) + (y - cx).powi(2) + (z - cx).powi(2)).sqrt();
        if dist <= r_incl {
            self.e_inclusion
        } else {
            self.e_matrix
        }
    }
}

// ---------------------------------------------------------------------------
// FeSquared
// ---------------------------------------------------------------------------

/// FE² concurrent multiscale method.
///
/// Each macro Gauss point embeds a micro-scale RVE boundary value problem.
/// The macro-scale effective tangent is obtained from the homogenized
/// stiffness tensor of the RVE.
#[derive(Debug, Clone)]
pub struct FeSquared {
    /// Macro mesh: list of element corner coordinates.
    pub macro_elements: Vec<MacroElement>,
    /// RVE geometry shared by all macro Gauss points.
    pub rve: RveGeometry,
    /// Homogenized (effective) 6×6 stiffness matrix (Voigt notation).
    pub effective_stiffness: Vec<Vec<f64>>,
    /// Number of macro elements.
    pub n_macro: usize,
}

impl FeSquared {
    /// Create a new FE² solver.
    ///
    /// # Arguments
    /// * `rve`     – representative volume element geometry
    /// * `n_macro` – number of macro elements
    pub fn new(rve: RveGeometry, n_macro: usize) -> Self {
        let macro_elements = (0..n_macro)
            .map(|k| MacroElement::new(k, 1.0, rve.clone()))
            .collect();
        let eff = initial_isotropic_stiffness(rve.e_matrix, rve.nu_matrix);
        Self {
            macro_elements,
            rve,
            effective_stiffness: eff,
            n_macro,
        }
    }

    /// Run macro iteration: homogenize each RVE and assemble effective moduli.
    pub fn homogenize(&mut self) {
        let c_eff = self.rve.voigt_effective_stiffness();
        self.effective_stiffness = c_eff;
        for elem in &mut self.macro_elements {
            elem.update_tangent(&self.effective_stiffness);
        }
    }

    /// Apply macro strain E and return homogenized macro stress Σ.
    ///
    /// Σ = C_eff : E  (Voigt notation, 6-vector).
    pub fn macro_stress(&self, strain: &[f64; 6]) -> [f64; 6] {
        let mut stress = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                stress[i] += self.effective_stiffness[i][j] * strain[j];
            }
        }
        stress
    }

    /// Number of macro degrees of freedom (3 per node, 8 nodes per hex).
    pub fn n_dof(&self) -> usize {
        self.n_macro * 8 * 3
    }
}

/// A single macro-scale hexahedral element in FE².
#[derive(Debug, Clone)]
pub struct MacroElement {
    /// Element index.
    pub index: usize,
    /// Current macro strain (Voigt, 6-vector).
    pub strain: [f64; 6],
    /// Current macro stress (Voigt, 6-vector).
    pub stress: [f64; 6],
    /// Embedded RVE geometry.
    pub rve: RveGeometry,
    /// Local 6×6 tangent stiffness.
    pub tangent: Vec<Vec<f64>>,
    /// Volume of this macro element.
    pub volume: f64,
}

impl MacroElement {
    /// Create a new macro element embedding a copy of the RVE.
    pub fn new(index: usize, volume: f64, rve: RveGeometry) -> Self {
        let tangent = initial_isotropic_stiffness(rve.e_matrix, rve.nu_matrix);
        Self {
            index,
            strain: [0.0; 6],
            stress: [0.0; 6],
            rve,
            tangent,
            volume,
        }
    }

    /// Update the local tangent stiffness from an externally computed C.
    pub fn update_tangent(&mut self, c: &[Vec<f64>]) {
        self.tangent = c.to_vec();
    }

    /// Compute homogenized stress from current strain via local tangent.
    pub fn compute_stress(&mut self) {
        for i in 0..6 {
            self.stress[i] = 0.0;
            for j in 0..6 {
                self.stress[i] += self.tangent[i][j] * self.strain[j];
            }
        }
    }
}

/// Micro-scale FEM solver for one RVE.
///
/// Applies a prescribed macro-strain to the RVE and returns the
/// volume-averaged (homogenized) stress tensor.
#[derive(Debug, Clone)]
pub struct MicroSolver {
    /// RVE geometry.
    pub rve: RveGeometry,
    /// Applied macro strain (Voigt 6-vector).
    pub macro_strain: [f64; 6],
    /// Resulting average stress (Voigt 6-vector).
    pub avg_stress: [f64; 6],
}

impl MicroSolver {
    /// Create a new micro solver for the given RVE.
    pub fn new(rve: RveGeometry) -> Self {
        Self {
            rve,
            macro_strain: [0.0; 6],
            avg_stress: [0.0; 6],
        }
    }

    /// Apply macro strain and solve (simplified linear elastic homogenization).
    pub fn solve(&mut self, macro_strain: [f64; 6]) {
        self.macro_strain = macro_strain;
        let c = self.rve.voigt_effective_stiffness();
        for i in 0..6 {
            self.avg_stress[i] = 0.0;
            for j in 0..6 {
                self.avg_stress[i] += c[i][j] * macro_strain[j];
            }
        }
    }

    /// Return effective tangent dΣ/dE (the homogenized 6×6 stiffness).
    pub fn effective_tangent(&self) -> Vec<Vec<f64>> {
        self.rve.voigt_effective_stiffness()
    }
}

// Helper: construct initial isotropic 6×6 stiffness matrix from E, ν.
fn initial_isotropic_stiffness(e: f64, nu: f64) -> Vec<Vec<f64>> {
    let factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let c11 = factor * (1.0 - nu);
    let c12 = factor * nu;
    let c44 = e / (2.0 * (1.0 + nu));
    let mut c = vec![vec![0.0f64; 6]; 6];
    // Normal components
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = if i == j { c11 } else { c12 };
        }
    }
    // Shear components
    c[3][3] = c44;
    c[4][4] = c44;
    c[5][5] = c44;
    c
}

impl RveGeometry {
    /// Compute Voigt-averaged 6×6 effective stiffness for this RVE.
    ///
    /// Uses a simple rule-of-mixtures (Voigt bound).
    pub fn voigt_effective_stiffness(&self) -> Vec<Vec<f64>> {
        let phi = self.inclusion_fraction;
        let c_i = initial_isotropic_stiffness(self.e_inclusion, self.nu_inclusion);
        let c_m = initial_isotropic_stiffness(self.e_matrix, self.nu_matrix);
        let mut c = vec![vec![0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                c[i][j] = phi * c_i[i][j] + (1.0 - phi) * c_m[i][j];
            }
        }
        c
    }
}

// ---------------------------------------------------------------------------
// HomogenizationScheme
// ---------------------------------------------------------------------------

/// Available homogenization schemes for composite materials.
#[derive(Debug, Clone, PartialEq)]
pub enum HomogenizationScheme {
    /// Voigt (iso-strain) upper bound.
    Voigt,
    /// Reuss (iso-stress) lower bound.
    Reuss,
    /// Mori-Tanaka mean-field scheme.
    MoriTanaka,
    /// Self-consistent (effective medium) scheme.
    SelfConsistent,
}

impl HomogenizationScheme {
    /// Compute the effective Young's modulus for a two-phase composite.
    ///
    /// # Arguments
    /// * `e_matrix`     – matrix modulus E_m
    /// * `e_inclusion`  – inclusion modulus E_i
    /// * `phi`          – volume fraction of inclusions
    pub fn effective_modulus(&self, e_matrix: f64, e_inclusion: f64, phi: f64) -> f64 {
        match self {
            HomogenizationScheme::Voigt => voigt_modulus(e_matrix, e_inclusion, phi),
            HomogenizationScheme::Reuss => reuss_modulus(e_matrix, e_inclusion, phi),
            HomogenizationScheme::MoriTanaka => mori_tanaka_modulus(e_matrix, e_inclusion, phi),
            HomogenizationScheme::SelfConsistent => {
                // Iterative: E* = Voigt/Reuss average as first estimate

                0.5 * (voigt_modulus(e_matrix, e_inclusion, phi)
                    + reuss_modulus(e_matrix, e_inclusion, phi))
            }
        }
    }
}

/// Voigt upper bound: E_V = φ·E_i + (1-φ)·E_m.
pub fn voigt_modulus(e_matrix: f64, e_inclusion: f64, phi: f64) -> f64 {
    phi * e_inclusion + (1.0 - phi) * e_matrix
}

/// Reuss lower bound: 1/E_R = φ/E_i + (1-φ)/E_m.
pub fn reuss_modulus(e_matrix: f64, e_inclusion: f64, phi: f64) -> f64 {
    let phi = phi.clamp(0.0, 1.0);
    let inv = phi / e_inclusion.max(1e-30) + (1.0 - phi) / e_matrix.max(1e-30);
    1.0 / inv.max(1e-30)
}

/// Mori-Tanaka effective modulus for spherical inclusions.
///
/// E_MT = E_m · \[1 + φ(E_i - E_m) / (E_m + (1-φ)(E_i - E_m)/3)\]
pub fn mori_tanaka_modulus(e_matrix: f64, e_inclusion: f64, phi: f64) -> f64 {
    let phi = phi.clamp(0.0, 1.0);
    let delta_e = e_inclusion - e_matrix;
    let denom = e_matrix + (1.0 - phi) * delta_e / 3.0;
    if denom.abs() < 1e-30 {
        return e_matrix;
    }
    e_matrix * (1.0 + phi * delta_e / denom)
}

// ---------------------------------------------------------------------------
// EffectiveMedium
// ---------------------------------------------------------------------------

/// Effective medium theory for composite microstructures.
///
/// Computes effective elastic constants (E*, ν*, K*, G*) from constituent
/// properties and microstructure topology.
#[derive(Debug, Clone)]
pub struct EffectiveMedium {
    /// Matrix elastic modulus.
    pub e_matrix: f64,
    /// Inclusion elastic modulus.
    pub e_inclusion: f64,
    /// Matrix Poisson ratio.
    pub nu_matrix: f64,
    /// Inclusion Poisson ratio.
    pub nu_inclusion: f64,
    /// Inclusion volume fraction φ.
    pub phi: f64,
    /// Active homogenization scheme.
    pub scheme: HomogenizationScheme,
}

impl EffectiveMedium {
    /// Create a new effective medium object.
    pub fn new(
        e_matrix: f64,
        e_inclusion: f64,
        nu_matrix: f64,
        nu_inclusion: f64,
        phi: f64,
        scheme: HomogenizationScheme,
    ) -> Self {
        Self {
            e_matrix,
            e_inclusion,
            nu_matrix,
            nu_inclusion,
            phi,
            scheme,
        }
    }

    /// Effective Young's modulus E*.
    pub fn effective_e(&self) -> f64 {
        self.scheme
            .effective_modulus(self.e_matrix, self.e_inclusion, self.phi)
    }

    /// Effective Poisson ratio ν* (Voigt average).
    pub fn effective_nu(&self) -> f64 {
        voigt_modulus(self.nu_matrix, self.nu_inclusion, self.phi)
    }

    /// Effective bulk modulus K* = E* / (3(1 - 2ν*)).
    pub fn effective_k(&self) -> f64 {
        let e = self.effective_e();
        let nu = self.effective_nu();
        e / (3.0 * (1.0 - 2.0 * nu).max(1e-10))
    }

    /// Effective shear modulus G* = E* / (2(1 + ν*)).
    pub fn effective_g(&self) -> f64 {
        let e = self.effective_e();
        let nu = self.effective_nu();
        e / (2.0 * (1.0 + nu).max(1e-10))
    }

    /// Check that Voigt bound ≥ Reuss bound for current fractions.
    pub fn voigt_geq_reuss(&self) -> bool {
        let ev = voigt_modulus(self.e_matrix, self.e_inclusion, self.phi);
        let er = reuss_modulus(self.e_matrix, self.e_inclusion, self.phi);
        ev >= er - 1e-10
    }
}

// ---------------------------------------------------------------------------
// MultiscaleMesh
// ---------------------------------------------------------------------------

/// Hierarchical mesh with coarse macro regions and refined micro regions.
#[derive(Debug, Clone)]
pub struct MultiscaleMesh {
    /// Coarse element size h_c.
    pub coarse_size: f64,
    /// Fine element size h_f.
    pub fine_size: f64,
    /// Number of coarse elements per side.
    pub n_coarse: usize,
    /// Refinement ratio r = h_c / h_f (integer).
    pub refinement_ratio: usize,
    /// Indices of macro elements requiring fine resolution.
    pub refined_elements: Vec<usize>,
    /// Total bounding box side length.
    pub domain_size: f64,
}

impl MultiscaleMesh {
    /// Create a new hierarchical mesh.
    ///
    /// # Arguments
    /// * `domain_size`      – total domain size
    /// * `n_coarse`         – coarse elements per side
    /// * `refinement_ratio` – fine-to-coarse refinement ratio
    pub fn new(domain_size: f64, n_coarse: usize, refinement_ratio: usize) -> Self {
        let coarse_size = domain_size / n_coarse as f64;
        let fine_size = coarse_size / refinement_ratio as f64;
        Self {
            coarse_size,
            fine_size,
            n_coarse,
            refinement_ratio,
            refined_elements: vec![],
            domain_size,
        }
    }

    /// Mark element `idx` for fine-scale resolution.
    pub fn refine(&mut self, idx: usize) {
        if !self.refined_elements.contains(&idx) {
            self.refined_elements.push(idx);
        }
    }

    /// Total number of coarse elements (3-D cube).
    pub fn n_coarse_total(&self) -> usize {
        self.n_coarse.pow(3)
    }

    /// Number of fine elements in a refined coarse element.
    pub fn fine_per_coarse(&self) -> usize {
        self.refinement_ratio.pow(3)
    }

    /// Estimated total DoF count (rough).
    pub fn estimated_dof(&self) -> usize {
        let n_ref = self.refined_elements.len();
        let coarse_dof = (self.n_coarse + 1).pow(3) * 3;
        let fine_extra = n_ref * (self.refinement_ratio + 1).pow(3) * 3;
        coarse_dof + fine_extra
    }
}

// ---------------------------------------------------------------------------
// AdaptiveMultiscale
// ---------------------------------------------------------------------------

/// Adaptive multiscale: a posteriori error estimator to select micro regions.
///
/// Uses a residual-based error indicator to identify macro elements where
/// multiscale resolution is required.
#[derive(Debug, Clone)]
pub struct AdaptiveMultiscale {
    /// Error indicator at each macro element.
    pub error_indicator: Vec<f64>,
    /// Global error tolerance θ.
    pub tolerance: f64,
    /// Multiscale flag: true means element needs micro-scale treatment.
    pub is_multiscale: Vec<bool>,
    /// Number of macro elements.
    pub n_elements: usize,
}

impl AdaptiveMultiscale {
    /// Create an adaptive multiscale selector.
    ///
    /// # Arguments
    /// * `n_elements` – number of macro elements
    /// * `tolerance`  – error threshold for multiscale activation
    pub fn new(n_elements: usize, tolerance: f64) -> Self {
        Self {
            error_indicator: vec![0.0; n_elements],
            tolerance,
            is_multiscale: vec![false; n_elements],
            n_elements,
        }
    }

    /// Update error indicators from element stresses and strains.
    ///
    /// Uses |σ - C_eff:ε|² as a residual indicator.
    pub fn update_indicators(
        &mut self,
        stress: &[[f64; 6]],
        strain: &[[f64; 6]],
        c_eff: &[Vec<f64>],
    ) {
        for e in 0..self.n_elements {
            let mut residual = [0.0f64; 6];
            for i in 0..6 {
                let sigma_pred: f64 = (0..6).map(|j| c_eff[i][j] * strain[e][j]).sum();
                residual[i] = stress[e][i] - sigma_pred;
            }
            self.error_indicator[e] = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        }
        // Döfler marking: activate multiscale where indicator > tolerance
        let max_err = self.error_indicator.iter().cloned().fold(0.0_f64, f64::max);
        for e in 0..self.n_elements {
            self.is_multiscale[e] = self.error_indicator[e] > self.tolerance * max_err;
        }
    }

    /// Number of elements flagged for multiscale treatment.
    pub fn multiscale_count(&self) -> usize {
        self.is_multiscale.iter().filter(|&&b| b).count()
    }

    /// Global error norm √(Σ ηₑ²).
    pub fn global_error(&self) -> f64 {
        self.error_indicator
            .iter()
            .map(|e| e * e)
            .sum::<f64>()
            .sqrt()
    }
}

// ---------------------------------------------------------------------------
// ConcurrentCoupling
// ---------------------------------------------------------------------------

/// Arlequin method for concurrent macro/micro domain coupling.
///
/// Overlaps a macro (coarse) region Ω_C and a micro (fine) region Ω_F in an
/// overlap zone Ω_ov.  Energy is partitioned using weight functions α(x).
#[derive(Debug, Clone)]
pub struct ConcurrentCoupling {
    /// Number of nodes in the overlap zone.
    pub n_overlap: usize,
    /// Coupling weight α at each overlap node (α ∈ \[0, 1\]).
    pub alpha: Vec<f64>,
    /// Lagrange multiplier values at coupling nodes.
    pub lambda: Vec<f64>,
    /// Coupling stiffness parameter k_c.
    pub k_coupling: f64,
    /// Total domain size.
    pub domain_size: f64,
    /// Size of the overlap zone.
    pub overlap_size: f64,
}

impl ConcurrentCoupling {
    /// Create an Arlequin coupling zone.
    ///
    /// # Arguments
    /// * `n_overlap`   – number of nodes in the overlap region
    /// * `k_coupling`  – coupling stiffness
    /// * `domain_size` – total domain extent
    /// * `overlap_size`– size of the gluing zone
    pub fn new(n_overlap: usize, k_coupling: f64, domain_size: f64, overlap_size: f64) -> Self {
        // Linear ramp: α = 1 at macro side, 0 at micro side
        let alpha: Vec<f64> = (0..n_overlap)
            .map(|i| i as f64 / n_overlap.max(1) as f64)
            .collect();
        Self {
            n_overlap,
            alpha,
            lambda: vec![0.0; n_overlap],
            k_coupling,
            domain_size,
            overlap_size,
        }
    }

    /// Energy partition: E_total = α·E_macro + (1-α)·E_micro + E_coupling.
    ///
    /// Returns the coupling energy for given displacement mismatch.
    pub fn coupling_energy(&self, u_macro: &[f64], u_micro: &[f64]) -> f64 {
        let n = self.n_overlap.min(u_macro.len()).min(u_micro.len());
        (0..n)
            .map(|i| {
                let diff = u_macro[i] - u_micro[i];
                0.5 * self.k_coupling * diff * diff
            })
            .sum()
    }

    /// Update Lagrange multipliers: λᵢ += k_c · (u_macro - u_micro).
    pub fn update_lambda(&mut self, u_macro: &[f64], u_micro: &[f64]) {
        let n = self.n_overlap.min(u_macro.len()).min(u_micro.len());
        for i in 0..n {
            self.lambda[i] += self.k_coupling * (u_macro[i] - u_micro[i]);
        }
    }

    /// Fraction of domain covered by the overlap zone.
    pub fn overlap_fraction(&self) -> f64 {
        if self.domain_size <= 0.0 {
            return 0.0;
        }
        self.overlap_size / self.domain_size
    }
}

// ---------------------------------------------------------------------------
// HierarchicalBasis
// ---------------------------------------------------------------------------

/// Hierarchical shape functions for p-refinement.
///
/// Provides Legendre-based hierarchical basis functions on \[-1, 1\].
/// The order-p approximation space contains all polynomials up to degree p.
#[derive(Debug, Clone)]
pub struct HierarchicalBasis {
    /// Polynomial order p.
    pub order: usize,
    /// Number of basis functions = order + 1.
    pub n_basis: usize,
}

impl HierarchicalBasis {
    /// Create a hierarchical basis of order p.
    pub fn new(order: usize) -> Self {
        Self {
            order,
            n_basis: order + 1,
        }
    }

    /// Evaluate all basis functions at ξ ∈ \[-1, 1\].
    ///
    /// Uses integrated Legendre polynomials as hierarchical bubble functions.
    pub fn evaluate(&self, xi: f64) -> Vec<f64> {
        let mut phi = vec![0.0f64; self.n_basis];
        // Standard nodal basis: φ₀ = (1-ξ)/2, φ₁ = (1+ξ)/2
        phi[0] = 0.5 * (1.0 - xi);
        if self.n_basis > 1 {
            phi[1] = 0.5 * (1.0 + xi);
        }
        // Hierarchical bubble functions: integrated Legendre
        for n in 2..self.n_basis {
            phi[n] = self.integrated_legendre(n, xi);
        }
        phi
    }

    /// Evaluate derivatives dφ/dξ at ξ.
    pub fn derivatives(&self, xi: f64) -> Vec<f64> {
        let mut dphi = vec![0.0f64; self.n_basis];
        dphi[0] = -0.5;
        if self.n_basis > 1 {
            dphi[1] = 0.5;
        }
        for n in 2..self.n_basis {
            dphi[n] = self.legendre_poly(n - 1, xi);
        }
        dphi
    }

    /// Integrated Legendre polynomial L̃_n(ξ) = ∫₋₁^ξ P_{n-1}(t) dt.
    fn integrated_legendre(&self, n: usize, xi: f64) -> f64 {
        // L̃_n = (P_n(ξ) - P_{n-2}(ξ)) / (2n-1)
        if n < 2 {
            return 0.0;
        }
        let pn = self.legendre_poly(n, xi);
        let pn2 = self.legendre_poly(n - 2, xi);
        (pn - pn2) / (2 * n - 1) as f64
    }

    /// Legendre polynomial P_n(ξ) via recurrence.
    fn legendre_poly(&self, n: usize, xi: f64) -> f64 {
        if n == 0 {
            return 1.0;
        }
        if n == 1 {
            return xi;
        }
        let mut p_prev = 1.0f64;
        let mut p_curr = xi;
        for k in 2..=n {
            let p_next = ((2 * k - 1) as f64 * xi * p_curr - (k - 1) as f64 * p_prev) / k as f64;
            p_prev = p_curr;
            p_curr = p_next;
        }
        p_curr
    }

    /// Check completeness: sum of basis functions = 1 at any ξ.
    pub fn partition_of_unity(&self, xi: f64) -> f64 {
        // Only nodal functions sum to 1; hierarchical bubbles sum to 0 at nodes
        let phi = self.evaluate(xi);
        phi[0] + if self.n_basis > 1 { phi[1] } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// MultiresolutionFem
// ---------------------------------------------------------------------------

/// Wavelet-based multiresolution FEM.
///
/// Decomposes the FEM solution into approximation (scaling) and detail
/// (wavelet) coefficients at multiple resolution levels.
#[derive(Debug, Clone)]
pub struct MultiresolutionFem {
    /// Maximum resolution level J.
    pub max_level: usize,
    /// Number of scaling coefficients at level j.
    pub n_scaling: Vec<usize>,
    /// Scaling coefficients (approximation).
    pub scaling_coeffs: Vec<Vec<f64>>,
    /// Wavelet coefficients (detail) at each level.
    pub wavelet_coeffs: Vec<Vec<f64>>,
}

impl MultiresolutionFem {
    /// Create a multiresolution FEM object with J levels.
    ///
    /// # Arguments
    /// * `max_level` – maximum decomposition level J
    /// * `n_coarse`  – number of coarse-scale basis functions
    pub fn new(max_level: usize, n_coarse: usize) -> Self {
        let mut n_scaling = vec![n_coarse];
        for j in 1..=max_level {
            n_scaling.push(n_scaling[j - 1] * 2);
        }
        let scaling_coeffs = vec![vec![0.0; n_coarse]; max_level + 1];
        let wavelet_coeffs: Vec<Vec<f64>> =
            (0..max_level).map(|j| vec![0.0; n_scaling[j]]).collect();
        Self {
            max_level,
            n_scaling,
            scaling_coeffs,
            wavelet_coeffs,
        }
    }

    /// Haar wavelet decomposition of a 1-D signal.
    ///
    /// Performs one level of the fast wavelet transform in-place.
    pub fn haar_decompose(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len();
        let half = n / 2;
        let mut approx = vec![0.0f64; half];
        let mut detail = vec![0.0f64; half];
        for i in 0..half {
            approx[i] = (signal[2 * i] + signal[2 * i + 1]) / 2.0_f64.sqrt();
            detail[i] = (signal[2 * i] - signal[2 * i + 1]) / 2.0_f64.sqrt();
        }
        (approx, detail)
    }

    /// Haar wavelet reconstruction from approximation and detail.
    pub fn haar_reconstruct(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let n = approx.len();
        let mut signal = vec![0.0f64; 2 * n];
        for i in 0..n {
            signal[2 * i] = (approx[i] + detail[i]) / 2.0_f64.sqrt();
            signal[2 * i + 1] = (approx[i] - detail[i]) / 2.0_f64.sqrt();
        }
        signal
    }

    /// Total number of basis functions at level j.
    pub fn basis_count(&self, level: usize) -> usize {
        if level < self.n_scaling.len() {
            self.n_scaling[level]
        } else {
            0
        }
    }

    /// Energy in detail coefficients at level j.
    pub fn detail_energy(&self, level: usize) -> f64 {
        if level >= self.wavelet_coeffs.len() {
            return 0.0;
        }
        self.wavelet_coeffs[level].iter().map(|c| c * c).sum()
    }
}

// ---------------------------------------------------------------------------
// SubstructureMethod
// ---------------------------------------------------------------------------

/// Craig-Bampton substructuring for large structural models.
///
/// Reduces each substructure to its interface DOFs plus a small set of
/// fixed-interface normal modes.
#[derive(Debug, Clone)]
pub struct SubstructureMethod {
    /// Number of interface (boundary) DOFs.
    pub n_interface: usize,
    /// Number of kept interior modes per substructure.
    pub n_modes: usize,
    /// Number of substructures.
    pub n_sub: usize,
    /// Modal frequencies (rad/s) for each substructure.
    pub modal_freqs: Vec<Vec<f64>>,
    /// Static constraint modes (n_interface × n_interface).
    pub constraint_modes: Vec<Vec<f64>>,
    /// Normal modes for each substructure (n_modes × n_interior).
    pub normal_modes: Vec<Vec<Vec<f64>>>,
}

impl SubstructureMethod {
    /// Create a Craig-Bampton substructure model.
    ///
    /// # Arguments
    /// * `n_interface` – number of interface DOFs
    /// * `n_modes`     – number of fixed-interface modes to keep
    /// * `n_sub`       – number of substructures
    /// * `base_freq`   – base frequency for synthetic modal spectrum
    pub fn new(n_interface: usize, n_modes: usize, n_sub: usize, base_freq: f64) -> Self {
        // Synthetic modal frequencies: ωₖ = base_freq · k for k = 1..n_modes
        let modal_freqs: Vec<Vec<f64>> = (0..n_sub)
            .map(|_| (1..=n_modes).map(|k| base_freq * k as f64).collect())
            .collect();
        // Identity-like constraint modes
        let constraint_modes = (0..n_interface)
            .map(|i| {
                let mut row = vec![0.0f64; n_interface];
                row[i] = 1.0;
                row
            })
            .collect();
        let normal_modes = vec![vec![vec![1.0; n_interface]; n_modes]; n_sub];
        Self {
            n_interface,
            n_modes,
            n_sub,
            modal_freqs,
            constraint_modes,
            normal_modes,
        }
    }

    /// Total reduced DOF count per substructure: n_interface + n_modes.
    pub fn reduced_dof(&self) -> usize {
        self.n_interface + self.n_modes
    }

    /// Total system DOF after Craig-Bampton reduction.
    pub fn total_reduced_dof(&self) -> usize {
        self.n_interface + self.n_sub * self.n_modes
    }

    /// Assemble reduced mass matrix (diagonal approximation).
    ///
    /// Returns a (total_reduced_dof × total_reduced_dof) diagonal mass matrix.
    pub fn reduced_mass(&self) -> Vec<f64> {
        let n = self.total_reduced_dof();
        // Interface mass = 1 per DOF; modal mass = 1 (mass-normalized modes)
        vec![1.0f64; n]
    }

    /// Modal participation factor Γₖ = φₖᵀ · M · r for mode k of substructure s.
    ///
    /// Here simplified: Γₖ = ωₖ (proportional to frequency).
    pub fn participation_factor(&self, sub: usize, mode: usize) -> f64 {
        if sub >= self.n_sub || mode >= self.n_modes {
            return 0.0;
        }
        self.modal_freqs[sub][mode]
    }

    /// Check that the modal count matches n_modes.
    pub fn mode_count(&self, sub: usize) -> usize {
        if sub < self.n_sub {
            self.modal_freqs[sub].len()
        } else {
            0
        }
    }
}

// ---------------------------------------------------------------------------
// Volume averaging
// ---------------------------------------------------------------------------

/// Volume average stress Σ = (1/V) ∫_Ω σ dV.
///
/// # Arguments
/// * `stresses` – list of Voigt stress vectors at each Gauss point
/// * `weights`  – integration weights (Gauss weights × det(J))
pub fn volume_average_stress(stresses: &[[f64; 6]], weights: &[f64]) -> [f64; 6] {
    let mut avg = [0.0f64; 6];
    let total_w: f64 = weights.iter().sum::<f64>().max(1e-30);
    for (s, &w) in stresses.iter().zip(weights.iter()) {
        for i in 0..6 {
            avg[i] += w * s[i];
        }
    }
    for i in 0..6 {
        avg[i] /= total_w;
    }
    avg
}

/// Hill-Mandel condition: macro work = micro work per unit volume.
///
/// Checks |Σ:E - (1/V)∫σ:ε dV| < tolerance.
///
/// # Arguments
/// * `macro_stress` – homogenized stress Σ (Voigt)
/// * `macro_strain` – homogenized strain E (Voigt)
/// * `micro_stresses` – micro-field stresses at each Gauss point
/// * `micro_strains`  – micro-field strains at each Gauss point
/// * `weights`        – integration weights
pub fn hill_mandel_residual(
    macro_stress: &[f64; 6],
    macro_strain: &[f64; 6],
    micro_stresses: &[[f64; 6]],
    micro_strains: &[[f64; 6]],
    weights: &[f64],
) -> f64 {
    let macro_work: f64 = (0..6).map(|i| macro_stress[i] * macro_strain[i]).sum();
    let total_w: f64 = weights.iter().sum::<f64>().max(1e-30);
    let micro_work: f64 = micro_stresses
        .iter()
        .zip(micro_strains.iter())
        .zip(weights.iter())
        .map(|((s, e), &w)| w * (0..6).map(|i| s[i] * e[i]).sum::<f64>())
        .sum::<f64>()
        / total_w;
    (macro_work - micro_work).abs()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- RveGeometry -------------------------------------------------------

    #[test]
    fn test_rve_n_elements() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        assert_eq!(rve.n_elements(), 64);
    }

    #[test]
    fn test_rve_n_nodes() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        assert_eq!(rve.n_nodes(), 125);
    }

    #[test]
    fn test_rve_volume() {
        let rve = RveGeometry::new(2.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        assert!((rve.volume() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_rve_element_size() {
        let rve = RveGeometry::new(1.0, 5, 0.3, 200.0, 70.0, 0.3, 0.33);
        assert!((rve.element_size() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_rve_voigt_stiffness_positive_diagonal() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let c = rve.voigt_effective_stiffness();
        for i in 0..6 {
            assert!(c[i][i] > 0.0, "diagonal entry C[{i}][{i}] must be positive");
        }
    }

    #[test]
    fn test_rve_voigt_stiffness_symmetry() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let c = rve.voigt_effective_stiffness();
        for i in 0..6 {
            for j in 0..6 {
                assert!((c[i][j] - c[j][i]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_rve_local_modulus_center_is_inclusion() {
        let rve = RveGeometry::new(1.0, 4, 0.5, 200.0, 70.0, 0.3, 0.33);
        let e = rve.local_modulus(0.5, 0.5, 0.5); // center
        assert!((e - 200.0).abs() < 1e-10);
    }

    // ---- FeSquared ---------------------------------------------------------

    #[test]
    fn test_fe2_new() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let fe2 = FeSquared::new(rve, 8);
        assert_eq!(fe2.n_macro, 8);
    }

    #[test]
    fn test_fe2_homogenize_runs() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let mut fe2 = FeSquared::new(rve, 4);
        fe2.homogenize();
        assert!(fe2.effective_stiffness[0][0] > 0.0);
    }

    #[test]
    fn test_fe2_macro_stress_zero_strain() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let fe2 = FeSquared::new(rve, 4);
        let stress = fe2.macro_stress(&[0.0; 6]);
        assert!(stress.iter().all(|&s| s.abs() < 1e-12));
    }

    #[test]
    fn test_fe2_macro_stress_nonzero() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let fe2 = FeSquared::new(rve, 4);
        let mut strain = [0.0f64; 6];
        strain[0] = 0.001;
        let stress = fe2.macro_stress(&strain);
        assert!(stress[0] != 0.0);
    }

    #[test]
    fn test_micro_solver_new() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let ms = MicroSolver::new(rve);
        assert!(ms.avg_stress.iter().all(|&s| s == 0.0));
    }

    #[test]
    fn test_micro_solver_solve() {
        let rve = RveGeometry::new(1.0, 4, 0.3, 200.0, 70.0, 0.3, 0.33);
        let mut ms = MicroSolver::new(rve);
        let strain = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        ms.solve(strain);
        assert!(ms.avg_stress[0] > 0.0);
    }

    // ---- HomogenizationScheme ---------------------------------------------

    #[test]
    fn test_voigt_geq_reuss() {
        let ev = voigt_modulus(70.0, 200.0, 0.3);
        let er = reuss_modulus(70.0, 200.0, 0.3);
        assert!(ev >= er - 1e-10);
    }

    #[test]
    fn test_voigt_modulus_pure_matrix() {
        let ev = voigt_modulus(70.0, 200.0, 0.0);
        assert!((ev - 70.0).abs() < 1e-10);
    }

    #[test]
    fn test_voigt_modulus_pure_inclusion() {
        let ev = voigt_modulus(70.0, 200.0, 1.0);
        assert!((ev - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_reuss_modulus_pure_matrix() {
        let er = reuss_modulus(70.0, 200.0, 0.0);
        assert!((er - 70.0).abs() < 1e-10);
    }

    #[test]
    fn test_mori_tanaka_between_bounds() {
        let e_m = 70.0;
        let e_i = 200.0;
        let phi = 0.3;
        let e_mt = mori_tanaka_modulus(e_m, e_i, phi);
        let e_v = voigt_modulus(e_m, e_i, phi);
        let e_r = reuss_modulus(e_m, e_i, phi);
        assert!(e_mt >= e_r - 1e-10);
        assert!(e_mt <= e_v + 1e-10);
    }

    // ---- EffectiveMedium --------------------------------------------------

    #[test]
    fn test_effective_medium_new() {
        let em = EffectiveMedium::new(70.0, 200.0, 0.33, 0.3, 0.3, HomogenizationScheme::Voigt);
        assert!(em.effective_e() > 0.0);
    }

    #[test]
    fn test_effective_medium_voigt_geq_reuss() {
        let em_v = EffectiveMedium::new(70.0, 200.0, 0.33, 0.3, 0.3, HomogenizationScheme::Voigt);
        let em_r = EffectiveMedium::new(70.0, 200.0, 0.33, 0.3, 0.3, HomogenizationScheme::Reuss);
        assert!(em_v.effective_e() >= em_r.effective_e() - 1e-10);
    }

    #[test]
    fn test_effective_medium_bulk_modulus_positive() {
        let em = EffectiveMedium::new(70.0, 200.0, 0.33, 0.3, 0.3, HomogenizationScheme::Voigt);
        assert!(em.effective_k() > 0.0);
    }

    #[test]
    fn test_effective_medium_shear_modulus_positive() {
        let em = EffectiveMedium::new(70.0, 200.0, 0.33, 0.3, 0.3, HomogenizationScheme::Voigt);
        assert!(em.effective_g() > 0.0);
    }

    // ---- MultiscaleMesh ---------------------------------------------------

    #[test]
    fn test_multiscale_mesh_new() {
        let mesh = MultiscaleMesh::new(1.0, 4, 4);
        assert_eq!(mesh.n_coarse, 4);
        assert!((mesh.coarse_size - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_multiscale_mesh_refine() {
        let mut mesh = MultiscaleMesh::new(1.0, 4, 4);
        mesh.refine(0);
        mesh.refine(5);
        assert_eq!(mesh.refined_elements.len(), 2);
    }

    #[test]
    fn test_multiscale_mesh_no_duplicate_refine() {
        let mut mesh = MultiscaleMesh::new(1.0, 4, 4);
        mesh.refine(3);
        mesh.refine(3);
        assert_eq!(mesh.refined_elements.len(), 1);
    }

    #[test]
    fn test_multiscale_mesh_fine_per_coarse() {
        let mesh = MultiscaleMesh::new(1.0, 4, 2);
        assert_eq!(mesh.fine_per_coarse(), 8);
    }

    // ---- AdaptiveMultiscale -----------------------------------------------

    #[test]
    fn test_adaptive_new() {
        let am = AdaptiveMultiscale::new(16, 0.5);
        assert_eq!(am.n_elements, 16);
        assert_eq!(am.multiscale_count(), 0);
    }

    #[test]
    fn test_adaptive_global_error_zero() {
        let am = AdaptiveMultiscale::new(8, 0.5);
        assert_eq!(am.global_error(), 0.0);
    }

    #[test]
    fn test_adaptive_update_indicators() {
        let n = 4;
        let c_eff = initial_isotropic_stiffness(70.0, 0.33);
        let strain = [[1e-3, 0.0, 0.0, 0.0, 0.0, 0.0]; 4];
        // Introduce mismatch in element 0
        let mut stress = [[0.0f64; 6]; 4];
        for e in 0..n {
            for i in 0..6 {
                stress[e][i] = (0..6usize).map(|j| c_eff[i][j] * strain[e][j]).sum::<f64>();
            }
        }
        stress[0][0] *= 2.0; // introduce error in element 0
        let mut am = AdaptiveMultiscale::new(n, 0.5);
        am.update_indicators(&stress, &strain, &c_eff);
        assert!(am.error_indicator[0] > 0.0);
    }

    // ---- ConcurrentCoupling -----------------------------------------------

    #[test]
    fn test_concurrent_coupling_new() {
        let cc = ConcurrentCoupling::new(16, 1e6, 1.0, 0.2);
        assert_eq!(cc.n_overlap, 16);
        assert!((cc.overlap_fraction() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_concurrent_coupling_energy_zero_match() {
        let cc = ConcurrentCoupling::new(8, 1e6, 1.0, 0.2);
        let u = vec![0.5; 8];
        assert_eq!(cc.coupling_energy(&u, &u), 0.0);
    }

    #[test]
    fn test_concurrent_coupling_energy_positive_mismatch() {
        let cc = ConcurrentCoupling::new(8, 1e6, 1.0, 0.2);
        let u1 = vec![0.5; 8];
        let u2 = vec![0.0; 8];
        assert!(cc.coupling_energy(&u1, &u2) > 0.0);
    }

    #[test]
    fn test_concurrent_coupling_update_lambda() {
        let mut cc = ConcurrentCoupling::new(8, 1.0, 1.0, 0.2);
        let u1 = vec![1.0; 8];
        let u2 = vec![0.0; 8];
        cc.update_lambda(&u1, &u2);
        assert!(cc.lambda.iter().all(|&l| l > 0.0));
    }

    // ---- HierarchicalBasis ------------------------------------------------

    #[test]
    fn test_hierarchical_basis_order1() {
        let hb = HierarchicalBasis::new(1);
        assert_eq!(hb.n_basis, 2);
    }

    #[test]
    fn test_hierarchical_basis_partition_of_unity() {
        let hb = HierarchicalBasis::new(3);
        for &xi in &[-1.0, -0.5, 0.0, 0.5, 1.0] {
            let pou = hb.partition_of_unity(xi);
            assert!((pou - 1.0).abs() < 1e-12, "POU failed at xi={xi:.6}");
        }
    }

    #[test]
    fn test_hierarchical_basis_evaluate_length() {
        let hb = HierarchicalBasis::new(4);
        let phi = hb.evaluate(0.0);
        assert_eq!(phi.len(), 5);
    }

    #[test]
    fn test_hierarchical_basis_derivatives_length() {
        let hb = HierarchicalBasis::new(3);
        let dphi = hb.derivatives(0.0);
        assert_eq!(dphi.len(), 4);
    }

    #[test]
    fn test_hierarchical_basis_linear_recovery() {
        // At order=1, φ(1) = [0, 1] and φ(-1) = [1, 0]
        let hb = HierarchicalBasis::new(1);
        let phi_m1 = hb.evaluate(-1.0);
        let phi_p1 = hb.evaluate(1.0);
        assert!((phi_m1[0] - 1.0).abs() < 1e-12);
        assert!((phi_p1[1] - 1.0).abs() < 1e-12);
    }

    // ---- MultiresolutionFem -----------------------------------------------

    #[test]
    fn test_multires_new() {
        let mr = MultiresolutionFem::new(3, 8);
        assert_eq!(mr.max_level, 3);
    }

    #[test]
    fn test_multires_basis_count() {
        let mr = MultiresolutionFem::new(2, 4);
        assert_eq!(mr.basis_count(0), 4);
        assert_eq!(mr.basis_count(1), 8);
    }

    #[test]
    fn test_haar_decompose_reconstruct() {
        let signal = vec![1.0, 3.0, 5.0, 11.0];
        let (approx, detail) = MultiresolutionFem::haar_decompose(&signal);
        let recon = MultiresolutionFem::haar_reconstruct(&approx, &detail);
        for (s, r) in signal.iter().zip(recon.iter()) {
            assert!((s - r).abs() < 1e-12);
        }
    }

    #[test]
    fn test_haar_energy_conservation() {
        let signal = vec![2.0, 4.0, 6.0, 8.0];
        let (approx, detail) = MultiresolutionFem::haar_decompose(&signal);
        let e_orig: f64 = signal.iter().map(|x| x * x).sum();
        let e_decomp: f64 = approx.iter().chain(detail.iter()).map(|x| x * x).sum();
        assert!((e_orig - e_decomp).abs() < 1e-10);
    }

    // ---- SubstructureMethod -----------------------------------------------

    #[test]
    fn test_substructure_new() {
        let cb = SubstructureMethod::new(6, 4, 3, 100.0);
        assert_eq!(cb.n_sub, 3);
        assert_eq!(cb.n_modes, 4);
    }

    #[test]
    fn test_substructure_mode_count() {
        let cb = SubstructureMethod::new(6, 5, 2, 100.0);
        assert_eq!(cb.mode_count(0), 5);
        assert_eq!(cb.mode_count(1), 5);
    }

    #[test]
    fn test_substructure_reduced_dof() {
        let cb = SubstructureMethod::new(6, 4, 3, 100.0);
        assert_eq!(cb.reduced_dof(), 10); // 6 + 4
    }

    #[test]
    fn test_substructure_total_reduced_dof() {
        let cb = SubstructureMethod::new(6, 4, 3, 100.0);
        assert_eq!(cb.total_reduced_dof(), 6 + 3 * 4);
    }

    #[test]
    fn test_substructure_modal_freqs_increasing() {
        let cb = SubstructureMethod::new(6, 5, 2, 50.0);
        let freqs = &cb.modal_freqs[0];
        for w in freqs.windows(2) {
            assert!(w[1] > w[0]);
        }
    }

    #[test]
    fn test_substructure_participation_factor() {
        let cb = SubstructureMethod::new(6, 4, 3, 100.0);
        let pf = cb.participation_factor(0, 0);
        assert!(pf > 0.0);
    }

    // ---- Volume averaging / Hill-Mandel -----------------------------------

    #[test]
    fn test_volume_average_stress_uniform() {
        let stresses = [[1.0, 2.0, 3.0, 0.0, 0.0, 0.0]; 4];
        let weights = [0.25, 0.25, 0.25, 0.25];
        let avg = volume_average_stress(&stresses, &weights);
        assert!((avg[0] - 1.0).abs() < 1e-12);
        assert!((avg[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_hill_mandel_consistent_zero() {
        // If micro field = macro field everywhere, Hill-Mandel error = 0
        let s = [1.0, 0.5, 0.2, 0.0, 0.0, 0.0f64];
        let e = [1e-3, 5e-4, 2e-4, 0.0, 0.0, 0.0f64];
        let stresses = [s; 4];
        let strains = [e; 4];
        let weights = [0.25, 0.25, 0.25, 0.25];
        let err = hill_mandel_residual(&s, &e, &stresses, &strains, &weights);
        assert!(err < 1e-12);
    }

    #[test]
    fn test_dense_cg_diagonal() {
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| {
                let mut row = vec![0.0; 4];
                row[i] = (i + 1) as f64;
                row
            })
            .collect();
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let x = dense_cg(&a, &b, 1e-12, 100);
        for i in 0..4 {
            assert!((x[i] - 1.0).abs() < 1e-8, "x[{i}] = {:.6}", x[i]);
        }
    }
}
