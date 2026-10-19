// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Inverse-homogenisation topology optimisation of 2-D periodic metamaterial
//! unit cells with a prescribed effective Poisson's ratio.
//!
//! The implementation follows the energy-based homogenisation approach of
//! Xia & Breitkopf, *"Design of materials using topology optimization and
//! energy-based homogenization approach in Matlab"*, Struct. Multidisc.
//! Optim. 52 (2015) 1229–1241:
//!
//! * the unit cell `Y = [0,1]²` is discretised by an `n × n` periodic mesh of
//!   bilinear (Q4) plane-stress elements;
//! * the element Young's modulus is interpolated with SIMP,
//!   `E(ρ) = E_min + ρᵖ (E₀ − E_min)`;
//! * the homogenised stiffness `Cᴴ` (Voigt, engineering shear) is obtained from
//!   three periodic unit-strain problems `K χⁱ = −fⁱ` and the mutual element
//!   energies `Cᴴ_ij = |Y|⁻¹ Σ_e (u_eⁱ)ᵀ K_e (u_eʲ)`, `u_eⁱ = u₀ⁱ + χ_eⁱ`;
//! * by self-adjointness of the energy formulation the sensitivities are
//!   `∂Cᴴ_ij/∂ρ_e = |Y|⁻¹ p ρ_eᵖ⁻¹ (E₀ − E_min) (u_eⁱ)ᵀ K_e⁰ (u_eʲ)`;
//! * a periodic linear density filter (radius `r_min`, default 1.5 elements),
//!   a smoothed Heaviside projection with β-continuation and an optional
//!   robust (eroded / intermediate / dilated) formulation regularise the
//!   design, which starts from a uniform field with a weakened central disc;
//! * the design is updated with the Method of Moving Asymptotes (Svanberg
//!   1987) under a move limit and a volume-fraction constraint, and the final
//!   layout is obtained by thresholding the physical densities at 0.5.
//!
//! The objective is
//!
//! ```text
//! f = (νᴴ − ν*)² + w · A₀ / A,    νᴴ = Cᴴ₁₂ / A,   A = ½ (Cᴴ₁₁ + Cᴴ₂₂)
//! ```
//!
//! i.e. the squared mismatch between the homogenised Poisson's ratio and the
//! target plus a small stiffness barrier.  The literal form
//! `(Cᴴ₁₂ − ν* Cᴴ₁₁)²` is minimised by the trivial all-void (disconnected)
//! cell, and the pure ratio mismatch still admits disconnected optima; the
//! barrier keeps the load paths connected.  Because the sensitivities of this
//! objective have both signs, the optimality-criteria fixed point used for
//! compliance problems is not a descent method here, so MMA is used instead.
//!
//! The global periodic stiffness matrix is factorised with a banded Cholesky
//! decomposition: rows and columns of the periodic node grid are numbered in
//! "folded" order (0, n−1, 1, n−2, …) so that the periodic wrap-around couples
//! only nodes at most two positions apart, giving a half-bandwidth of `O(4 n)`.

use crate::error::{Error, Result};

/// Number of degrees of freedom of a bilinear plane element.
const ELEMENT_DOFS: usize = 8;

/// Gauss–Legendre abscissa for 2-point quadrature.
const GAUSS_POINT: f64 = 0.577_350_269_189_625_8;

/// Spatial symmetry imposed on a periodic unit-cell design.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellSymmetry {
    /// No symmetry constraint.
    None,
    /// Mirror symmetry about both mid-lines of the cell (orthotropic design).
    Orthotropic,
    /// Mirror symmetry about both mid-lines and both diagonals (square/D₄
    /// symmetry).  Guarantees `Cᴴ₁₁ = Cᴴ₂₂` and `Cᴴ₁₃ = Cᴴ₂₃ = 0`.
    Square,
}

/// SIMP material model for a two-phase (solid/void) plane-stress design.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimpMaterial {
    /// Young's modulus of the solid phase `E₀` (Pa).
    pub solid_modulus: f64,
    /// Young's modulus assigned to void elements `E_min` (Pa), `0 < E_min < E₀`.
    pub void_modulus: f64,
    /// Poisson's ratio of the solid phase `ν₀`.
    pub solid_poisson: f64,
    /// SIMP penalisation exponent `p ≥ 1`.
    pub penalty: f64,
}

impl Default for SimpMaterial {
    fn default() -> Self {
        Self {
            solid_modulus: 1.0,
            void_modulus: 1e-9,
            solid_poisson: 0.3,
            penalty: 3.0,
        }
    }
}

impl SimpMaterial {
    /// Validate the material parameters.
    pub fn validate(&self) -> Result<()> {
        check_range(
            "solid_modulus",
            self.solid_modulus,
            f64::MIN_POSITIVE,
            f64::MAX,
        )?;
        check_range(
            "void_modulus",
            self.void_modulus,
            f64::MIN_POSITIVE,
            f64::MAX,
        )?;
        if self.void_modulus >= self.solid_modulus {
            return Err(Error::ParameterOutOfRange {
                param: "void_modulus",
                value: self.void_modulus,
                min: f64::MIN_POSITIVE,
                max: self.solid_modulus,
            });
        }
        // Plane-stress constitutive matrix is positive definite for |ν| < 1;
        // physically admissible isotropic solids satisfy −1 < ν ≤ 0.5.
        check_open_closed("solid_poisson", self.solid_poisson, -1.0, 0.5)?;
        check_range("penalty", self.penalty, 1.0, 16.0)?;
        Ok(())
    }

    /// SIMP-interpolated Young's modulus for density `rho ∈ [0, 1]`.
    pub fn modulus(&self, rho: f64) -> f64 {
        self.void_modulus + rho.powf(self.penalty) * (self.solid_modulus - self.void_modulus)
    }

    /// Derivative `dE/dρ` of the SIMP interpolation.
    pub fn modulus_derivative(&self, rho: f64) -> f64 {
        if rho <= 0.0 {
            if self.penalty <= 1.0 {
                self.solid_modulus - self.void_modulus
            } else {
                0.0
            }
        } else {
            self.penalty * rho.powf(self.penalty - 1.0) * (self.solid_modulus - self.void_modulus)
        }
    }
}

/// Homogenised in-plane stiffness of a periodic unit cell.
///
/// Voigt ordering `(σ₁₁, σ₂₂, σ₁₂) = C · (ε₁₁, ε₂₂, γ₁₂)` with engineering
/// shear strain `γ₁₂ = 2 ε₁₂`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomogenizedStiffness {
    /// Homogenised 3 × 3 stiffness matrix `Cᴴ` (Pa).
    pub c: [[f64; 3]; 3],
}

impl HomogenizedStiffness {
    /// Compliance matrix `Sᴴ = (Cᴴ)⁻¹`.
    pub fn compliance(&self) -> Result<[[f64; 3]; 3]> {
        invert_3x3(&self.c)
    }

    /// Poisson's ratio `ν₁₂ = −S₂₁ / S₁₁` (transverse contraction along 2 for
    /// uniaxial stress along 1).
    pub fn poisson_ratio_12(&self) -> Result<f64> {
        let s = self.compliance()?;
        Ok(-s[1][0] / s[0][0])
    }

    /// Poisson's ratio `ν₂₁ = −S₁₂ / S₂₂` (transverse contraction along 1 for
    /// uniaxial stress along 2).
    pub fn poisson_ratio_21(&self) -> Result<f64> {
        let s = self.compliance()?;
        Ok(-s[0][1] / s[1][1])
    }

    /// Effective Young's modulus along direction 1, `E₁ = 1 / S₁₁`.
    pub fn youngs_modulus_1(&self) -> Result<f64> {
        let s = self.compliance()?;
        Ok(1.0 / s[0][0])
    }

    /// Effective Young's modulus along direction 2, `E₂ = 1 / S₂₂`.
    pub fn youngs_modulus_2(&self) -> Result<f64> {
        let s = self.compliance()?;
        Ok(1.0 / s[1][1])
    }

    /// Effective in-plane shear modulus `G₁₂ = 1 / S₃₃`.
    pub fn shear_modulus(&self) -> Result<f64> {
        let s = self.compliance()?;
        Ok(1.0 / s[2][2])
    }

    /// Stiffness-ratio Poisson measure `νᴴ = C₁₂ / ½(C₁₁ + C₂₂)` used as the
    /// optimisation objective.  Equals `ν₁₂ = ν₂₁` for square-symmetric cells.
    pub fn stiffness_poisson_ratio(&self) -> f64 {
        let avg = 0.5 * (self.c[0][0] + self.c[1][1]);
        if avg.abs() < f64::MIN_POSITIVE {
            return f64::NAN;
        }
        self.c[0][1] / avg
    }
}

/// Configuration of the inverse-homogenisation Poisson's-ratio optimisation.
#[derive(Debug, Clone, PartialEq)]
pub struct AuxeticTopologyConfig {
    /// Number of elements along each side of the square unit cell (`n ≥ 4`).
    pub resolution: usize,
    /// Target effective Poisson's ratio `ν*`, `−1 < ν* < 1`.
    pub target_poisson: f64,
    /// Prescribed solid volume fraction, `0 < V < 1`.
    pub volume_fraction: f64,
    /// SIMP material description.
    pub material: SimpMaterial,
    /// Density-filter radius in element widths (`r_min > 0`; `≤ 1` disables
    /// smoothing).
    pub filter_radius: f64,
    /// MMA move limit per iteration, `0 < m ≤ 1`.
    pub move_limit: f64,
    /// Maximum number of optimisation iterations (0 evaluates the initial
    /// design only).
    pub max_iterations: usize,
    /// Convergence tolerance on the maximum design-variable change.
    pub tolerance: f64,
    /// Final Heaviside projection sharpness `β_max` reached by continuation
    /// (β starts at 1 and doubles whenever the design converges).  `None`
    /// disables projection, so the filtered densities are used directly.
    pub projection_beta_max: Option<f64>,
    /// Number of iterations after which β is doubled even if the design has
    /// not yet settled (`≥ 1`).
    pub projection_interval: usize,
    /// Robust formulation (Wang, Lazarov & Sigmund 2011): with `Some(δ)` the
    /// objective is the sum over the dilated (η = 0.5 − δ), intermediate
    /// (η = 0.5) and eroded (η = 0.5 + δ) projections, which removes one-node
    /// hinges and makes the thresholded layout perform like the optimised
    /// design.  Requires `projection_beta_max`; `0 < δ < 0.5`.
    pub robust_threshold_offset: Option<f64>,
    /// Weight `w ≥ 0` of the stiffness barrier `w · A₀ / A` added to the
    /// Poisson-ratio mismatch (`A = ½(Cᴴ₁₁ + Cᴴ₂₂)`, `A₀` its initial value).
    /// The barrier keeps the load paths connected; `0` gives the pure ratio
    /// objective, which admits a degenerate disconnected (all-void) optimum.
    pub stiffness_weight: f64,
    /// Symmetry imposed on the design.
    pub symmetry: CellSymmetry,
}

impl Default for AuxeticTopologyConfig {
    fn default() -> Self {
        Self {
            resolution: 32,
            target_poisson: -0.3,
            volume_fraction: 0.4,
            material: SimpMaterial::default(),
            filter_radius: 1.5,
            move_limit: 0.2,
            max_iterations: 300,
            tolerance: 0.01,
            projection_beta_max: Some(32.0),
            projection_interval: 40,
            robust_threshold_offset: Some(0.2),
            stiffness_weight: 0.01,
            symmetry: CellSymmetry::Square,
        }
    }
}

impl AuxeticTopologyConfig {
    /// Validate every parameter of the configuration.
    pub fn validate(&self) -> Result<()> {
        if self.resolution < 4 || self.resolution > 512 {
            return Err(Error::ParameterOutOfRange {
                param: "resolution",
                value: self.resolution as f64,
                min: 4.0,
                max: 512.0,
            });
        }
        check_open("target_poisson", self.target_poisson, -1.0, 1.0)?;
        check_open("volume_fraction", self.volume_fraction, 0.0, 1.0)?;
        self.material.validate()?;
        check_range(
            "filter_radius",
            self.filter_radius,
            f64::MIN_POSITIVE,
            self.resolution as f64,
        )?;
        check_open_closed("move_limit", self.move_limit, 0.0, 1.0)?;
        check_range("tolerance", self.tolerance, f64::MIN_POSITIVE, 1.0)?;
        if let Some(beta) = self.projection_beta_max {
            check_range("projection_beta_max", beta, 1.0, 1024.0)?;
        }
        if self.projection_interval == 0 {
            return Err(Error::ParameterOutOfRange {
                param: "projection_interval",
                value: 0.0,
                min: 1.0,
                max: f64::INFINITY,
            });
        }
        check_range("stiffness_weight", self.stiffness_weight, 0.0, f64::MAX)?;
        if let Some(delta) = self.robust_threshold_offset {
            check_open("robust_threshold_offset", delta, 0.0, 0.5)?;
            if self.projection_beta_max.is_none() {
                return Err(Error::MissingProperty("projection_beta_max"));
            }
        }
        Ok(())
    }
}

/// Outcome of [`optimize_auxetic_topology`].
#[derive(Debug, Clone)]
pub struct AuxeticTopologyResult {
    /// Physical (filtered and projected) element densities, `densities[row][col]`.
    pub densities: Vec<Vec<f64>>,
    /// Binary layout obtained by thresholding the physical densities at 0.5
    /// (`true` = solid), `layout[row][col]`.
    pub layout: Vec<Vec<bool>>,
    /// Homogenised stiffness of the grey-scale optimised design.
    pub stiffness: HomogenizedStiffness,
    /// Effective Poisson's ratio `ν₁₂` of the grey-scale design.
    pub effective_poisson: f64,
    /// Homogenised stiffness of the thresholded binary layout (void = `E_min`).
    pub layout_stiffness: HomogenizedStiffness,
    /// Effective Poisson's ratio `ν₁₂` of the binary layout.
    pub layout_poisson: f64,
    /// Solid volume fraction of the binary layout.
    pub layout_volume_fraction: f64,
    /// Number of optimisation iterations performed.
    pub iterations: usize,
    /// `true` when the design change dropped below the tolerance at the final
    /// projection sharpness before `max_iterations` was reached.
    pub converged: bool,
    /// Objective value `(νᴴ − ν*)² + w A₀ / A` at every iteration.
    pub objective_history: Vec<f64>,
}

/// Run energy-based inverse-homogenisation topology optimisation of a periodic
/// square unit cell for a target effective Poisson's ratio.
///
/// The initial design is the prescribed volume fraction with a weakened
/// central circular region (radius `n/3`), which breaks the uniform-density
/// stationary point exactly as in Xia & Breitkopf (2015).
pub fn optimize_auxetic_topology(config: &AuxeticTopologyConfig) -> Result<AuxeticTopologyResult> {
    config.validate()?;
    let n = config.resolution;
    let n_elem = n * n;
    let mesh = PeriodicQ4Mesh::new(n, config.material.solid_poisson);
    let filter = PeriodicDensityFilter::new(n, config.filter_radius);
    let symmetry = SymmetryMap::new(n, config.symmetry);
    let target = config.target_poisson;
    let vol_frac = config.volume_fraction;

    let mut x = initial_design(n, vol_frac);
    symmetry.symmetrize(&mut x);
    let beta_max = config.projection_beta_max;
    let mut beta = beta_max.map(|_| 1.0);
    let mut history = Vec::with_capacity(config.max_iterations);
    let mut converged = false;
    let mut iterations = 0;
    let mut since_beta_change = 0;
    let mut reference_stiffness: Option<f64> = None;
    let mut objective_scale: Option<f64> = None;
    let mut mma = MmaState::new(&x);

    // Projection thresholds of the designs entering the objective: the
    // intermediate (blueprint) design alone, or the dilated / intermediate /
    // eroded triple of the robust formulation.
    let thresholds: Vec<f64> = match config.robust_threshold_offset {
        Some(delta) => vec![0.5 - delta, 0.5, 0.5 + delta],
        None => vec![0.5],
    };
    let inv_n = 1.0 / n_elem as f64;

    while iterations < config.max_iterations {
        let rho_tilde = filter.apply(&x);
        let mut objective = 0.0;
        let mut d_obj_tilde = vec![0.0; n_elem];
        let mut volume = 0.0;
        let mut dv_tilde = vec![0.0; n_elem];
        for &eta in &thresholds {
            let (rho_phys, dproj) = project(&rho_tilde, beta, eta);
            let analysis = mesh.homogenize(&rho_phys, &config.material, true)?;
            if (eta - 0.5).abs() < f64::EPSILON {
                let mean_normal = 0.5 * (analysis.stiffness.c[0][0] + analysis.stiffness.c[1][1]);
                reference_stiffness.get_or_insert(mean_normal);
                volume = rho_phys.iter().sum::<f64>() * inv_n;
                for (dv, p) in dv_tilde.iter_mut().zip(&dproj) {
                    *dv = p * inv_n;
                }
            }
            let a_ref = reference_stiffness.unwrap_or(1.0);
            let (f_k, d_obj) =
                poisson_objective(&analysis, target, config.stiffness_weight, a_ref)?;
            objective += f_k;
            // Chain rule through the projection.
            for ((acc, d), p) in d_obj_tilde.iter_mut().zip(&d_obj).zip(&dproj) {
                *acc += d * p;
            }
        }
        history.push(objective);

        // Chain rule through the (self-adjoint) density filter.
        let mut dc = filter.apply_transpose(&d_obj_tilde);
        let mut dv = filter.apply_transpose(&dv_tilde);
        symmetry.symmetrize(&mut dc);
        symmetry.symmetrize(&mut dv);

        // MMA on scaled functions: objective × 10 / f₀(x⁰), volume constraint
        // g = V(x)/V* − 1 ≤ 0.
        let scale = 10.0 / *objective_scale.get_or_insert(objective.abs().max(1e-12));
        let df0: Vec<f64> = dc.iter().map(|v| v * scale).collect();
        let g = volume / vol_frac - 1.0;
        let dg: Vec<f64> = dv.iter().map(|v| v / vol_frac).collect();
        let mut x_new = mma.update(&x, &df0, g, &dg, config.move_limit)?;
        symmetry.symmetrize(&mut x_new);
        let change = x
            .iter()
            .zip(&x_new)
            .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
        x = x_new;
        iterations += 1;
        since_beta_change += 1;

        // Converge only on a settled, volume-feasible design that has had a few
        // iterations to adapt to the current projection sharpness.
        let settled = change < config.tolerance && g <= 1e-3 && since_beta_change >= 5;
        match (beta, beta_max) {
            (Some(b), Some(b_max))
                if b < b_max && (settled || since_beta_change >= config.projection_interval) =>
            {
                beta = Some((2.0 * b).min(b_max));
                since_beta_change = 0;
            }
            _ if settled => {
                converged = true;
                break;
            }
            _ => {}
        }
    }

    let rho_tilde = filter.apply(&x);
    let (rho_phys, _) = project(&rho_tilde, beta, 0.5);
    let final_analysis = mesh.homogenize(&rho_phys, &config.material, false)?;
    let stiffness = final_analysis.stiffness;
    let effective_poisson = stiffness.poisson_ratio_12()?;

    let layout_flat: Vec<bool> = rho_phys.iter().map(|&r| r > 0.5).collect();
    let layout_rho: Vec<f64> = layout_flat
        .iter()
        .map(|&s| if s { 1.0 } else { 0.0 })
        .collect();
    let layout_analysis = mesh.homogenize(&layout_rho, &config.material, false)?;
    let layout_stiffness = layout_analysis.stiffness;
    let layout_poisson = layout_stiffness.poisson_ratio_12()?;
    let layout_volume_fraction = layout_rho.iter().sum::<f64>() / n_elem as f64;

    Ok(AuxeticTopologyResult {
        densities: to_grid(&rho_phys, n),
        layout: to_grid(&layout_flat, n),
        stiffness,
        effective_poisson,
        layout_stiffness,
        layout_poisson,
        layout_volume_fraction,
        iterations,
        converged,
        objective_history: history,
    })
}

/// Compute the homogenised plane-stress stiffness of a periodic unit cell.
///
/// `densities[row][col]` holds element densities in `[0, 1]` on a square
/// `n × n` grid (`n ≥ 2`); the element modulus follows the SIMP law of
/// `material`.  Uses the same energy-based periodic finite-element
/// homogenisation as [`optimize_auxetic_topology`].
pub fn homogenize_periodic_plane_stress(
    densities: &[Vec<f64>],
    material: &SimpMaterial,
) -> Result<HomogenizedStiffness> {
    material.validate()?;
    let n = densities.len();
    if n < 2 {
        return Err(Error::ParameterOutOfRange {
            param: "resolution",
            value: n as f64,
            min: 2.0,
            max: f64::INFINITY,
        });
    }
    let mut flat = Vec::with_capacity(n * n);
    for row in densities {
        if row.len() != n {
            return Err(Error::DimensionMismatch {
                expected: format!("{n} x {n} density grid"),
                actual: format!("row of length {}", row.len()),
            });
        }
        for &rho in row {
            if !rho.is_finite() || !(0.0..=1.0).contains(&rho) {
                return Err(Error::ParameterOutOfRange {
                    param: "density",
                    value: rho,
                    min: 0.0,
                    max: 1.0,
                });
            }
            flat.push(rho);
        }
    }
    let mesh = PeriodicQ4Mesh::new(n, material.solid_poisson);
    Ok(mesh.homogenize(&flat, material, false)?.stiffness)
}

/// Generate a binary periodic unit cell whose homogenised Poisson's ratio
/// approaches `target_poisson`.
///
/// Convenience wrapper around [`optimize_auxetic_topology`] using the default
/// configuration (32 × 32 cell, volume fraction 0.4, SIMP p = 3, filter radius
/// 1.5, move limit 0.2, square symmetry) with at most `iterations`
/// optimisation iterations.  `iterations = 0` returns the 0.5-thresholded
/// initial design; because that grey field (density 0.4 with a weakened
/// central disc) lies below the threshold everywhere, the layout is entirely
/// void (`false`).
///
/// Returns the thresholded layout `layout[row][col]` (`true` = solid).
pub fn topology_optimized_auxetic(
    target_poisson: f64,
    iterations: usize,
) -> Result<Vec<Vec<bool>>> {
    let config = AuxeticTopologyConfig {
        target_poisson,
        max_iterations: iterations,
        ..AuxeticTopologyConfig::default()
    };
    Ok(optimize_auxetic_topology(&config)?.layout)
}

// ─────────────────────────────────────────────────────────────────────────────
// Parameter validation helpers
// ─────────────────────────────────────────────────────────────────────────────

fn check_range(param: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value.is_finite() && value >= min && value <= max {
        Ok(())
    } else {
        Err(Error::ParameterOutOfRange {
            param,
            value,
            min,
            max,
        })
    }
}

fn check_open(param: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value.is_finite() && value > min && value < max {
        Ok(())
    } else {
        Err(Error::ParameterOutOfRange {
            param,
            value,
            min,
            max,
        })
    }
}

fn check_open_closed(param: &'static str, value: f64, min: f64, max: f64) -> Result<()> {
    if value.is_finite() && value > min && value <= max {
        Ok(())
    } else {
        Err(Error::ParameterOutOfRange {
            param,
            value,
            min,
            max,
        })
    }
}

fn invert_3x3(m: &[[f64; 3]; 3]) -> Result<[[f64; 3]; 3]> {
    let cof = [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
        ],
        [
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
        ],
        [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ];
    let det = m[0][0] * cof[0][0] + m[0][1] * cof[0][1] + m[0][2] * cof[0][2];
    let scale = m.iter().flatten().fold(0.0_f64, |acc, v| acc.max(v.abs()));
    if !det.is_finite() || det.abs() <= 1e-14 * scale.powi(3) || scale == 0.0 {
        return Err(Error::Singularity {
            context: "homogenized stiffness inverse",
            detail: format!("determinant {det:.3e} is not invertible"),
        });
    }
    let mut inv = [[0.0; 3]; 3];
    for (i, row) in inv.iter_mut().enumerate() {
        for (j, val) in row.iter_mut().enumerate() {
            // Inverse = adjugate / det, adjugate = cofactorᵀ.
            *val = cof[j][i] / det;
        }
    }
    Ok(inv)
}

fn to_grid<T: Copy>(flat: &[T], n: usize) -> Vec<Vec<T>> {
    flat.chunks(n).map(<[T]>::to_vec).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Design initialisation, projection, objective and MMA update
// ─────────────────────────────────────────────────────────────────────────────

/// Uniform `vol_frac` with a central disc (radius n/3) at half that density.
fn initial_design(n: usize, vol_frac: f64) -> Vec<f64> {
    let centre = 0.5 * n as f64;
    let radius = n as f64 / 3.0;
    let mut x = vec![vol_frac; n * n];
    for (idx, val) in x.iter_mut().enumerate() {
        let row = (idx / n) as f64 + 0.5;
        let col = (idx % n) as f64 + 0.5;
        let dist = ((row - centre).powi(2) + (col - centre).powi(2)).sqrt();
        if dist < radius {
            *val = 0.5 * vol_frac;
        }
    }
    x
}

/// Smoothed Heaviside projection `ρ̄ = H_{β,η}(ρ̃)` with threshold `eta`.
///
/// Returns the projected densities and `dρ̄/dρ̃`.  With `beta = None` the
/// projection is the identity.
fn project(rho_tilde: &[f64], beta: Option<f64>, eta: f64) -> (Vec<f64>, Vec<f64>) {
    match beta {
        None => (rho_tilde.to_vec(), vec![1.0; rho_tilde.len()]),
        Some(b) => {
            let denom = (b * eta).tanh() + (b * (1.0 - eta)).tanh();
            let base = (b * eta).tanh();
            let mut rho = Vec::with_capacity(rho_tilde.len());
            let mut deriv = Vec::with_capacity(rho_tilde.len());
            for &r in rho_tilde {
                let t = (b * (r - eta)).tanh();
                rho.push(((base + t) / denom).clamp(0.0, 1.0));
                deriv.push(b * (1.0 - t * t) / denom);
            }
            (rho, deriv)
        }
    }
}

/// Objective `f = (νᴴ − ν*)² + w A₀ / A` with `νᴴ = C₁₂ / A`,
/// `A = ½(C₁₁ + C₂₂)`, and its derivative with respect to the physical
/// densities:
///
/// `df = 2 (νᴴ − ν*) / A · (dC₁₂ − νᴴ dA) − w A₀ / A² · dA`.
fn poisson_objective(
    analysis: &HomogenizationAnalysis,
    target: f64,
    stiffness_weight: f64,
    reference_stiffness: f64,
) -> Result<(f64, Vec<f64>)> {
    let c = &analysis.stiffness.c;
    let avg = 0.5 * (c[0][0] + c[1][1]);
    if !(avg.is_finite() && avg > 0.0) {
        return Err(Error::Singularity {
            context: "Poisson-ratio objective",
            detail: format!("non-positive mean normal stiffness {avg:.3e}"),
        });
    }
    let nu = c[0][1] / avg;
    let mismatch = nu - target;
    let barrier = stiffness_weight * reference_stiffness / avg;
    let objective = mismatch * mismatch + barrier;
    let factor = 2.0 * mismatch / avg;
    let barrier_factor = barrier / avg;
    let grad = analysis
        .sensitivities
        .iter()
        .map(|d| {
            // d = [dC11, dC22, dC33, dC12, dC13, dC23]
            let d_avg = 0.5 * (d[0] + d[1]);
            factor * (d[3] - nu * d_avg) - barrier_factor * d_avg
        })
        .collect();
    Ok((objective, grad))
}

/// Method of Moving Asymptotes (Svanberg 1987) for one objective and a single
/// inequality constraint `g(x) ≤ 0` with box bounds `0 ≤ x ≤ 1`.
///
/// The objective of the inverse-homogenisation problem is non-monotone in the
/// densities (its sensitivities have both signs), for which the classical
/// optimality-criteria fixed point is not a descent method; MMA builds a
/// convex separable approximation from first-order information and moving
/// asymptotes, and its single-constraint subproblem is solved exactly through
/// its one-dimensional dual.
struct MmaState {
    low: Vec<f64>,
    upp: Vec<f64>,
    x_old1: Vec<f64>,
    x_old2: Vec<f64>,
    iteration: usize,
}

impl MmaState {
    const ASY_INIT: f64 = 0.5;
    const ASY_INCR: f64 = 1.2;
    const ASY_DECR: f64 = 0.7;
    const ALBEFA: f64 = 0.1;
    const RAA0: f64 = 1e-5;

    fn new(x: &[f64]) -> Self {
        Self {
            low: vec![0.0; x.len()],
            upp: vec![1.0; x.len()],
            x_old1: x.to_vec(),
            x_old2: x.to_vec(),
            iteration: 0,
        }
    }

    /// One MMA step.  `df0` is the objective gradient, `g` the constraint
    /// value and `dg` its gradient at `x`.
    fn update(
        &mut self,
        x: &[f64],
        df0: &[f64],
        g: f64,
        dg: &[f64],
        move_limit: f64,
    ) -> Result<Vec<f64>> {
        if df0.iter().chain(dg).any(|v| !v.is_finite()) || !g.is_finite() {
            return Err(Error::NumericalDivergence {
                solver: "auxetic topology MMA update",
                iterations: self.iteration,
                residual: g,
            });
        }
        self.iteration += 1;
        // Asymptote update (bounds xmin = 0, xmax = 1).
        for (j, &xj) in x.iter().enumerate() {
            if self.iteration <= 2 {
                self.low[j] = xj - Self::ASY_INIT;
                self.upp[j] = xj + Self::ASY_INIT;
            } else {
                let trend = (xj - self.x_old1[j]) * (self.x_old1[j] - self.x_old2[j]);
                let factor = if trend > 0.0 {
                    Self::ASY_INCR
                } else if trend < 0.0 {
                    Self::ASY_DECR
                } else {
                    1.0
                };
                self.low[j] =
                    (xj - factor * (self.x_old1[j] - self.low[j])).clamp(xj - 10.0, xj - 0.01);
                self.upp[j] =
                    (xj + factor * (self.upp[j] - self.x_old1[j])).clamp(xj + 0.01, xj + 10.0);
            }
        }
        // Move limits and convex approximation coefficients.
        let n = x.len();
        let mut alpha = vec![0.0; n];
        let mut beta = vec![0.0; n];
        let mut p0 = vec![0.0; n];
        let mut q0 = vec![0.0; n];
        let mut p1 = vec![0.0; n];
        let mut q1 = vec![0.0; n];
        let mut r1 = g;
        for j in 0..n {
            let (xj, lj, uj) = (x[j], self.low[j], self.upp[j]);
            alpha[j] = (lj + Self::ALBEFA * (xj - lj))
                .max(xj - move_limit)
                .max(0.0);
            beta[j] = (uj - Self::ALBEFA * (uj - xj))
                .min(xj + move_limit)
                .min(1.0);
            let ux2 = (uj - xj) * (uj - xj);
            let xl2 = (xj - lj) * (xj - lj);
            let (gp, gm) = (df0[j].max(0.0), (-df0[j]).max(0.0));
            p0[j] = ux2 * (1.001 * gp + 0.001 * gm + Self::RAA0);
            q0[j] = xl2 * (0.001 * gp + 1.001 * gm + Self::RAA0);
            let (hp, hm) = (dg[j].max(0.0), (-dg[j]).max(0.0));
            p1[j] = ux2 * (1.001 * hp + 0.001 * hm + Self::RAA0);
            q1[j] = xl2 * (0.001 * hp + 1.001 * hm + Self::RAA0);
            r1 -= p1[j] / (uj - xj) + q1[j] / (xj - lj);
        }
        let low = &self.low;
        let upp = &self.upp;
        let primal = |lambda: f64| -> Vec<f64> {
            (0..n)
                .map(|j| {
                    let pj = (p0[j] + lambda * p1[j]).sqrt();
                    let qj = (q0[j] + lambda * q1[j]).sqrt();
                    let xj = (pj * low[j] + qj * upp[j]) / (pj + qj);
                    xj.clamp(alpha[j], beta[j])
                })
                .collect()
        };
        let approx_constraint = |xs: &[f64]| -> f64 {
            r1 + xs
                .iter()
                .enumerate()
                .map(|(j, &xj)| p1[j] / (upp[j] - xj) + q1[j] / (xj - low[j]))
                .sum::<f64>()
        };
        // Dual: find λ ≥ 0 with complementary slackness on the approximated
        // constraint (it is non-increasing in λ).
        let mut x_new = primal(0.0);
        if approx_constraint(&x_new) > 0.0 {
            let mut lam_lo = 0.0;
            let mut lam_hi = 1.0;
            let mut expansions = 0;
            while approx_constraint(&primal(lam_hi)) > 0.0 && expansions < 200 {
                lam_lo = lam_hi;
                lam_hi *= 4.0;
                expansions += 1;
            }
            for _ in 0..200 {
                let mid = 0.5 * (lam_lo + lam_hi);
                if approx_constraint(&primal(mid)) > 0.0 {
                    lam_lo = mid;
                } else {
                    lam_hi = mid;
                }
                if lam_hi - lam_lo <= 1e-12 * lam_hi.max(1.0) {
                    break;
                }
            }
            x_new = primal(lam_hi);
        }
        self.x_old2 = std::mem::replace(&mut self.x_old1, x.to_vec());
        Ok(x_new)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Symmetry
// ─────────────────────────────────────────────────────────────────────────────

/// Element orbits under the imposed symmetry group.
struct SymmetryMap {
    orbits: Vec<Vec<usize>>,
}

impl SymmetryMap {
    fn new(n: usize, symmetry: CellSymmetry) -> Self {
        let mut orbit_of = vec![usize::MAX; n * n];
        let mut orbits = Vec::new();
        if symmetry == CellSymmetry::None {
            return Self { orbits };
        }
        for row in 0..n {
            for col in 0..n {
                let idx = row * n + col;
                if orbit_of[idx] != usize::MAX {
                    continue;
                }
                let mirror_r = n - 1 - row;
                let mirror_c = n - 1 - col;
                let mut images = vec![
                    (row, col),
                    (mirror_r, col),
                    (row, mirror_c),
                    (mirror_r, mirror_c),
                ];
                if symmetry == CellSymmetry::Square {
                    images.extend([
                        (col, row),
                        (mirror_c, row),
                        (col, mirror_r),
                        (mirror_c, mirror_r),
                    ]);
                }
                let mut members: Vec<usize> = images.iter().map(|&(r, c)| r * n + c).collect();
                members.sort_unstable();
                members.dedup();
                let orbit_id = orbits.len();
                for &m in &members {
                    orbit_of[m] = orbit_id;
                }
                orbits.push(members);
            }
        }
        Self { orbits }
    }

    /// Replace each value by the mean over its symmetry orbit.
    fn symmetrize(&self, values: &mut [f64]) {
        for orbit in &self.orbits {
            let mean = orbit.iter().map(|&i| values[i]).sum::<f64>() / orbit.len() as f64;
            for &i in orbit {
                values[i] = mean;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic density filter
// ─────────────────────────────────────────────────────────────────────────────

/// Linear "hat" density filter `w = max(0, r_min − d)` on the periodic grid.
struct PeriodicDensityFilter {
    n: usize,
    /// `(row offset, col offset, normalised weight)`.
    kernel: Vec<(isize, isize, f64)>,
}

impl PeriodicDensityFilter {
    fn new(n: usize, radius: f64) -> Self {
        let reach = radius.ceil() as isize;
        let mut kernel = Vec::new();
        for di in -reach..=reach {
            for dj in -reach..=reach {
                let dist = ((di * di + dj * dj) as f64).sqrt();
                let w = radius - dist;
                if w > 0.0 {
                    kernel.push((di, dj, w));
                }
            }
        }
        if kernel.is_empty() {
            kernel.push((0, 0, 1.0));
        }
        let total: f64 = kernel.iter().map(|k| k.2).sum();
        for k in &mut kernel {
            k.2 /= total;
        }
        Self { n, kernel }
    }

    fn wrap(&self, base: usize, offset: isize) -> usize {
        let n = self.n as isize;
        (base as isize + offset).rem_euclid(n) as usize
    }

    /// `ρ̃ = H x` (weights are uniform-normalised, so `H` is a symmetric
    /// circulant operator and `Hᵀ = H`).
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0; n * n];
        for (idx, val) in out.iter_mut().enumerate() {
            let row = idx / n;
            let col = idx % n;
            *val = self
                .kernel
                .iter()
                .map(|&(di, dj, w)| w * x[self.wrap(row, di) * n + self.wrap(col, dj)])
                .sum();
        }
        out
    }

    /// `Hᵀ g` — equal to `H g` because the kernel is point-symmetric.
    fn apply_transpose(&self, g: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0; n * n];
        for (idx, val) in out.iter_mut().enumerate() {
            let row = idx / n;
            let col = idx % n;
            *val = self
                .kernel
                .iter()
                .map(|&(di, dj, w)| w * g[self.wrap(row, -di) * n + self.wrap(col, -dj)])
                .sum();
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic Q4 finite-element homogenisation
// ─────────────────────────────────────────────────────────────────────────────

/// Result of one homogenisation analysis.
struct HomogenizationAnalysis {
    stiffness: HomogenizedStiffness,
    /// Per-element `[dC11, dC22, dC33, dC12, dC13, dC23] / dρ_e`.
    sensitivities: Vec<[f64; 6]>,
}

/// Periodic `n × n` mesh of unit-modulus square Q4 plane-stress elements.
struct PeriodicQ4Mesh {
    n: usize,
    /// Unit-modulus element stiffness matrix (element side length 1).
    ke: [[f64; ELEMENT_DOFS]; ELEMENT_DOFS],
    /// Global DOF indices of every element (folded node ordering).
    elem_dofs: Vec<[usize; ELEMENT_DOFS]>,
    /// Nodal displacements of the three macroscopic unit-strain fields.
    u0: [[f64; ELEMENT_DOFS]; 3],
    /// `K_e⁰ u₀ⁱ` for the three unit-strain fields.
    ke_u0: [[f64; ELEMENT_DOFS]; 3],
    n_dofs: usize,
    half_bandwidth: usize,
}

impl PeriodicQ4Mesh {
    fn new(n: usize, nu: f64) -> Self {
        let ke = q4_plane_stress_stiffness(nu);
        // Folded ordering (0, n−1, 1, n−2, …) of rows and columns keeps
        // periodically adjacent rows/columns at most two positions apart.
        let row_pos = |r: usize| -> usize {
            if 2 * r < n {
                2 * r
            } else {
                2 * (n - 1 - r) + 1
            }
        };
        let node_index = |r: usize, c: usize| -> usize { row_pos(r % n) * n + row_pos(c % n) };
        let mut elem_dofs = Vec::with_capacity(n * n);
        let mut half_bandwidth = 0;
        for i in 0..n {
            for j in 0..n {
                // Local node order: (x, y) = (0,0), (1,0), (1,1), (0,1);
                // row index i runs along y, column index j along x.
                let nodes = [
                    node_index(i, j),
                    node_index(i, j + 1),
                    node_index(i + 1, j + 1),
                    node_index(i + 1, j),
                ];
                let mut dofs = [0usize; ELEMENT_DOFS];
                for (k, &node) in nodes.iter().enumerate() {
                    dofs[2 * k] = 2 * node;
                    dofs[2 * k + 1] = 2 * node + 1;
                }
                let lo = dofs.iter().copied().min().unwrap_or(0);
                let hi = dofs.iter().copied().max().unwrap_or(0);
                half_bandwidth = half_bandwidth.max(hi - lo);
                elem_dofs.push(dofs);
            }
        }
        // Macroscopic unit strains ε = (1,0,0), (0,1,0), γ = 1 on local nodes.
        let coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let mut u0 = [[0.0; ELEMENT_DOFS]; 3];
        for (k, &(x, y)) in coords.iter().enumerate() {
            u0[0][2 * k] = x;
            u0[1][2 * k + 1] = y;
            u0[2][2 * k] = 0.5 * y;
            u0[2][2 * k + 1] = 0.5 * x;
        }
        let mut ke_u0 = [[0.0; ELEMENT_DOFS]; 3];
        for (case, out) in ke_u0.iter_mut().enumerate() {
            for (row, val) in out.iter_mut().enumerate() {
                *val = ke[row].iter().zip(&u0[case]).map(|(k, u)| k * u).sum();
            }
        }
        Self {
            n,
            ke,
            elem_dofs,
            u0,
            ke_u0,
            n_dofs: 2 * n * n,
            half_bandwidth,
        }
    }

    /// Homogenise the cell for physical densities `rho` (row-major).
    fn homogenize(
        &self,
        rho: &[f64],
        material: &SimpMaterial,
        with_sensitivities: bool,
    ) -> Result<HomogenizationAnalysis> {
        let moduli: Vec<f64> = rho.iter().map(|&r| material.modulus(r)).collect();
        // Node 0 (position 0, DOFs 0 and 1) is pinned to remove the rigid
        // translations of the periodic fluctuation field.
        let pinned = |dof: usize| dof < 2;
        let mut k_glob = BandedSpdMatrix::zeros(self.n_dofs, self.half_bandwidth);
        let mut rhs = vec![vec![0.0; self.n_dofs]; 3];
        for (dofs, &e_mod) in self.elem_dofs.iter().zip(&moduli) {
            for (a, &ga) in dofs.iter().enumerate() {
                if pinned(ga) {
                    continue;
                }
                for (b, &gb) in dofs.iter().enumerate() {
                    if pinned(gb) || gb > ga {
                        continue;
                    }
                    k_glob.add(ga, gb, e_mod * self.ke[a][b]);
                }
                for (case_rhs, case_ke_u0) in rhs.iter_mut().zip(&self.ke_u0) {
                    case_rhs[ga] -= e_mod * case_ke_u0[a];
                }
            }
        }
        k_glob.add(0, 0, 1.0);
        k_glob.add(1, 1, 1.0);
        k_glob.factorize()?;
        let mut chi = Vec::with_capacity(3);
        for case_rhs in rhs {
            chi.push(k_glob.solve(case_rhs));
        }

        let area = (self.n * self.n) as f64;
        let mut c = [[0.0; 3]; 3];
        let mut sensitivities = if with_sensitivities {
            Vec::with_capacity(self.elem_dofs.len())
        } else {
            Vec::new()
        };
        // (i, j) pairs in the order [11, 22, 33, 12, 13, 23].
        const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for ((dofs, &e_mod), &r) in self.elem_dofs.iter().zip(&moduli).zip(rho) {
            let mut u_e = [[0.0; ELEMENT_DOFS]; 3];
            for (case, ue) in u_e.iter_mut().enumerate() {
                for (k, (val, &g)) in ue.iter_mut().zip(dofs).enumerate() {
                    *val = self.u0[case][k] + chi[case][g];
                }
            }
            let mut ke_ue = [[0.0; ELEMENT_DOFS]; 3];
            for (out, ue) in ke_ue.iter_mut().zip(&u_e) {
                for (val, ke_row) in out.iter_mut().zip(&self.ke) {
                    *val = ke_row.iter().zip(ue).map(|(k, u)| k * u).sum();
                }
            }
            let mut energies = [0.0; 6];
            for (energy, &(i, j)) in energies.iter_mut().zip(&PAIRS) {
                *energy = u_e[i]
                    .iter()
                    .zip(&ke_ue[j])
                    .map(|(u, k)| u * k)
                    .sum::<f64>();
            }
            for (energy, &(i, j)) in energies.iter().zip(&PAIRS) {
                let contrib = e_mod * energy / area;
                c[i][j] += contrib;
                if i != j {
                    c[j][i] += contrib;
                }
            }
            if with_sensitivities {
                let de = material.modulus_derivative(r) / area;
                let mut sens = [0.0; 6];
                for (s, energy) in sens.iter_mut().zip(&energies) {
                    *s = de * energy;
                }
                sensitivities.push(sens);
            }
        }
        if c.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::NumericalDivergence {
                solver: "periodic homogenization",
                iterations: 0,
                residual: f64::NAN,
            });
        }
        Ok(HomogenizationAnalysis {
            stiffness: HomogenizedStiffness { c },
            sensitivities,
        })
    }
}

/// Unit-modulus stiffness matrix of a unit square bilinear plane-stress
/// element, integrated with 2 × 2 Gauss quadrature (exact for this element).
///
/// Local node order (x, y): (0,0), (1,0), (1,1), (0,1); DOFs `(uₓ, u_y)` per node.
fn q4_plane_stress_stiffness(nu: f64) -> [[f64; ELEMENT_DOFS]; ELEMENT_DOFS] {
    let factor = 1.0 / (1.0 - nu * nu);
    let d = [
        [factor, factor * nu, 0.0],
        [factor * nu, factor, 0.0],
        [0.0, 0.0, factor * 0.5 * (1.0 - nu)],
    ];
    let xi_nodes = [-1.0, 1.0, 1.0, -1.0];
    let eta_nodes = [-1.0, -1.0, 1.0, 1.0];
    let mut ke = [[0.0; ELEMENT_DOFS]; ELEMENT_DOFS];
    for &xi in &[-GAUSS_POINT, GAUSS_POINT] {
        for &eta in &[-GAUSS_POINT, GAUSS_POINT] {
            // Mapping [-1,1]² → [0,1]²: dx/dξ = dy/dη = ½, det J = ¼.
            let mut b = [[0.0; ELEMENT_DOFS]; 3];
            for (k, (&xk, &ek)) in xi_nodes.iter().zip(&eta_nodes).enumerate() {
                let dn_dxi = 0.25 * xk * (1.0 + ek * eta);
                let dn_deta = 0.25 * ek * (1.0 + xk * xi);
                let dn_dx = 2.0 * dn_dxi;
                let dn_dy = 2.0 * dn_deta;
                b[0][2 * k] = dn_dx;
                b[1][2 * k + 1] = dn_dy;
                b[2][2 * k] = dn_dy;
                b[2][2 * k + 1] = dn_dx;
            }
            let det_j = 0.25;
            // ke += Bᵀ D B det J
            let mut db = [[0.0; ELEMENT_DOFS]; 3];
            for (db_row, d_row) in db.iter_mut().zip(&d) {
                for (col, val) in db_row.iter_mut().enumerate() {
                    *val = d_row
                        .iter()
                        .zip(&b)
                        .map(|(dv, b_row)| dv * b_row[col])
                        .sum();
                }
            }
            for (a, ke_row) in ke.iter_mut().enumerate() {
                for (c, val) in ke_row.iter_mut().enumerate() {
                    let sum: f64 = b
                        .iter()
                        .zip(&db)
                        .map(|(b_row, db_row)| b_row[a] * db_row[c])
                        .sum();
                    *val += sum * det_j;
                }
            }
        }
    }
    ke
}

// ─────────────────────────────────────────────────────────────────────────────
// Banded symmetric positive-definite matrix with in-place Cholesky
// ─────────────────────────────────────────────────────────────────────────────

/// Lower band of a symmetric matrix; row `i` stores columns `i − bw ..= i`.
struct BandedSpdMatrix {
    size: usize,
    bw: usize,
    data: Vec<f64>,
}

impl BandedSpdMatrix {
    fn zeros(size: usize, bw: usize) -> Self {
        Self {
            size,
            bw,
            data: vec![0.0; size * (bw + 1)],
        }
    }

    /// Storage index of entry `(i, j)` with `j ≤ i ≤ j + bw`.
    fn index(&self, i: usize, j: usize) -> usize {
        i * (self.bw + 1) + (j + self.bw - i)
    }

    /// Add `value` to entry `(i, j)` of the lower triangle (`j ≤ i`).
    fn add(&mut self, i: usize, j: usize, value: f64) {
        let idx = self.index(i, j);
        self.data[idx] += value;
    }

    /// In-place Cholesky factorisation `A = L Lᵀ`.
    fn factorize(&mut self) -> Result<()> {
        let width = self.bw + 1;
        for i in 0..self.size {
            let lo = i.saturating_sub(self.bw);
            for j in lo..=i {
                let row_i = i * width + (lo + self.bw - i);
                let row_j = j * width + (lo + self.bw - j);
                let len = j - lo;
                let dot: f64 = self.data[row_i..row_i + len]
                    .iter()
                    .zip(&self.data[row_j..row_j + len])
                    .map(|(a, b)| a * b)
                    .sum();
                let idx = self.index(i, j);
                let value = self.data[idx] - dot;
                if i == j {
                    if !(value.is_finite() && value > 0.0) {
                        return Err(Error::Singularity {
                            context: "periodic homogenization stiffness",
                            detail: format!("non-positive Cholesky pivot {value:.3e} at row {i}"),
                        });
                    }
                    self.data[idx] = value.sqrt();
                } else {
                    let diag = self.data[self.index(j, j)];
                    self.data[idx] = value / diag;
                }
            }
        }
        Ok(())
    }

    /// Solve `L Lᵀ x = b` using the factor computed by [`Self::factorize`].
    fn solve(&self, mut b: Vec<f64>) -> Vec<f64> {
        let width = self.bw + 1;
        // Forward substitution L y = b.
        for i in 0..self.size {
            let lo = i.saturating_sub(self.bw);
            let row = i * width + (lo + self.bw - i);
            let dot: f64 = self.data[row..row + (i - lo)]
                .iter()
                .zip(&b[lo..i])
                .map(|(l, y)| l * y)
                .sum();
            b[i] = (b[i] - dot) / self.data[self.index(i, i)];
        }
        // Backward substitution Lᵀ x = y (column-oriented).
        for i in (0..self.size).rev() {
            let xi = b[i] / self.data[self.index(i, i)];
            b[i] = xi;
            let lo = i.saturating_sub(self.bw);
            let row = i * width + (lo + self.bw - i);
            for (y, l) in b[lo..i].iter_mut().zip(&self.data[row..row + (i - lo)]) {
                *y -= l * xi;
            }
        }
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane_stress_d(e: f64, nu: f64) -> [[f64; 3]; 3] {
        let f = e / (1.0 - nu * nu);
        [
            [f, f * nu, 0.0],
            [f * nu, f, 0.0],
            [0.0, 0.0, f * 0.5 * (1.0 - nu)],
        ]
    }

    fn central_hole(n: usize, radius: f64) -> Vec<Vec<f64>> {
        let c = 0.5 * n as f64;
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let d =
                            ((i as f64 + 0.5 - c).powi(2) + (j as f64 + 0.5 - c).powi(2)).sqrt();
                        if d < radius { 0.0 } else { 1.0 }
                    })
                    .collect()
            })
            .collect()
    }

    fn fast_config(n: usize, target: f64) -> AuxeticTopologyConfig {
        AuxeticTopologyConfig {
            resolution: n,
            target_poisson: target,
            ..AuxeticTopologyConfig::default()
        }
    }

    #[test]
    fn element_stiffness_is_symmetric_with_rigid_body_null_space() {
        let ke = q4_plane_stress_stiffness(0.3);
        for (a, row) in ke.iter().enumerate() {
            for (b, &kab) in row.iter().enumerate() {
                assert!((kab - ke[b][a]).abs() < 1e-14);
            }
        }
        // Rigid translations and infinitesimal rotation produce zero forces.
        let coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let mut modes = [[0.0; ELEMENT_DOFS]; 3];
        for (k, &(x, y)) in coords.iter().enumerate() {
            modes[0][2 * k] = 1.0;
            modes[1][2 * k + 1] = 1.0;
            modes[2][2 * k] = -y;
            modes[2][2 * k + 1] = x;
        }
        for mode in &modes {
            for row in &ke {
                let f: f64 = row.iter().zip(mode).map(|(k, u)| k * u).sum();
                assert!(f.abs() < 1e-13, "rigid mode force {f}");
            }
        }
    }

    #[test]
    fn banded_cholesky_solves_spd_system() {
        // Symmetric diagonally dominant band matrix, half-bandwidth 3.
        let size = 12;
        let bw = 3;
        let mut dense = vec![vec![0.0; size]; size];
        let mut band = BandedSpdMatrix::zeros(size, bw);
        let entry = |i: usize, j: usize| -> f64 {
            let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
            if hi - lo > bw {
                0.0
            } else if hi == lo {
                10.0 + hi as f64
            } else {
                1.0 / (1.0 + (hi + 2 * lo) as f64)
            }
        };
        for (i, row) in dense.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = entry(i, j);
                if j <= i && i - j <= bw {
                    band.add(i, j, *val);
                }
            }
        }
        let x_true: Vec<f64> = (0..size).map(|i| (i as f64 * 0.37).sin() + 0.5).collect();
        let b: Vec<f64> = dense
            .iter()
            .map(|row| row.iter().zip(&x_true).map(|(a, x)| a * x).sum())
            .collect();
        band.factorize().expect("SPD factorisation");
        let x = band.solve(b);
        for (xi, ti) in x.iter().zip(&x_true) {
            assert!((xi - ti).abs() < 1e-12, "{xi} vs {ti}");
        }
    }

    #[test]
    fn banded_cholesky_rejects_indefinite_matrix() {
        let mut band = BandedSpdMatrix::zeros(2, 1);
        band.add(0, 0, 1.0);
        band.add(1, 0, 2.0);
        band.add(1, 1, 1.0);
        assert!(band.factorize().is_err());
    }

    #[test]
    fn solid_cell_reproduces_plane_stress_stiffness() {
        let material = SimpMaterial {
            solid_modulus: 2.5,
            ..SimpMaterial::default()
        };
        for n in [3, 6, 9] {
            let rho = vec![vec![1.0; n]; n];
            let ch = homogenize_periodic_plane_stress(&rho, &material).expect("homogenize");
            let d = plane_stress_d(2.5, 0.3);
            for (i, (c_row, d_row)) in ch.c.iter().zip(&d).enumerate() {
                for (j, (&cij, &dij)) in c_row.iter().zip(d_row).enumerate() {
                    assert!(
                        (cij - dij).abs() < 1e-10,
                        "n={n} C[{i}][{j}] = {cij} vs {dij}"
                    );
                }
            }
            let nu = ch.poisson_ratio_12().expect("nu");
            assert!((nu - 0.3).abs() < 1e-10);
            let e1 = ch.youngs_modulus_1().expect("E1");
            assert!((e1 - 2.5).abs() < 1e-9);
            let g = ch.shear_modulus().expect("G");
            assert!((g - 2.5 / 2.6).abs() < 1e-9);
        }
    }

    #[test]
    fn layered_cell_matches_voigt_and_reuss_bounds() {
        // Horizontal layers (rows) alternating solid / soft: stiffness along
        // the layers (direction 1) follows the Voigt average for the
        // in-plane-stretch term, normal to the layers the Reuss average.
        let n = 8;
        let material = SimpMaterial {
            solid_modulus: 1.0,
            void_modulus: 0.1,
            solid_poisson: 0.0,
            penalty: 1.0,
        };
        let rho: Vec<Vec<f64>> = (0..n)
            .map(|i| vec![if i < n / 2 { 1.0 } else { 0.0 }; n])
            .collect();
        let ch = homogenize_periodic_plane_stress(&rho, &material).expect("homogenize");
        // ν₀ = 0 ⇒ the phases decouple: C11 = mean(E), C22 = 1/mean(1/E).
        let voigt = 0.5 * (1.0 + 0.1);
        let reuss = 1.0 / (0.5 * (1.0 + 10.0));
        assert!((ch.c[0][0] - voigt).abs() < 1e-10, "C11 = {}", ch.c[0][0]);
        assert!((ch.c[1][1] - reuss).abs() < 1e-10, "C22 = {}", ch.c[1][1]);
        assert!(ch.c[0][1].abs() < 1e-12);
    }

    #[test]
    fn porous_cell_is_symmetric_square_and_softer() {
        let n = 12;
        let rho = central_hole(n, 3.5);
        let material = SimpMaterial::default();
        let ch = homogenize_periodic_plane_stress(&rho, &material).expect("homogenize");
        for i in 0..3 {
            for j in 0..3 {
                assert!((ch.c[i][j] - ch.c[j][i]).abs() < 1e-12);
            }
        }
        assert!((ch.c[0][0] - ch.c[1][1]).abs() < 1e-10 * ch.c[0][0]);
        assert!(ch.c[0][2].abs() < 1e-10 && ch.c[1][2].abs() < 1e-10);
        let d = plane_stress_d(1.0, 0.3);
        assert!(ch.c[0][0] < d[0][0] && ch.c[0][0] > 0.2 * d[0][0]);
        let nu = ch.poisson_ratio_12().expect("nu");
        assert!(nu > 0.0 && nu < 0.3, "porous ν = {nu}");
    }

    #[test]
    fn stiffness_sensitivities_match_finite_differences() {
        let n = 6;
        let material = SimpMaterial::default();
        let mesh = PeriodicQ4Mesh::new(n, material.solid_poisson);
        let rho: Vec<f64> = (0..n * n)
            .map(|k| 0.3 + 0.6 * ((k as f64 * 0.731).sin() * 0.5 + 0.5))
            .collect();
        let base = mesh.homogenize(&rho, &material, true).expect("analysis");
        const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        let h = 1e-6;
        for &e in &[0usize, 7, 20, 35] {
            let mut plus = rho.clone();
            plus[e] += h;
            let mut minus = rho.clone();
            minus[e] -= h;
            let cp = mesh
                .homogenize(&plus, &material, false)
                .expect("plus")
                .stiffness
                .c;
            let cm = mesh
                .homogenize(&minus, &material, false)
                .expect("minus")
                .stiffness
                .c;
            for (k, &(i, j)) in PAIRS.iter().enumerate() {
                let fd = (cp[i][j] - cm[i][j]) / (2.0 * h);
                let an = base.sensitivities[e][k];
                assert!(
                    (fd - an).abs() < 1e-6 * an.abs().max(1e-3),
                    "element {e} pair ({i},{j}): fd={fd} analytic={an}"
                );
            }
        }
    }

    #[test]
    fn objective_gradient_matches_finite_differences() {
        let n = 6;
        let material = SimpMaterial::default();
        let mesh = PeriodicQ4Mesh::new(n, material.solid_poisson);
        let rho: Vec<f64> = (0..n * n)
            .map(|k| 0.25 + 0.7 * ((k as f64 * 1.37).cos() * 0.5 + 0.5))
            .collect();
        let (target, w, a_ref) = (-0.3, 0.01, 0.2);
        let base = mesh.homogenize(&rho, &material, true).expect("analysis");
        let (_, grad) = poisson_objective(&base, target, w, a_ref).expect("objective");
        let h = 1e-6;
        for &e in &[1usize, 13, 30] {
            let mut plus = rho.clone();
            plus[e] += h;
            let mut minus = rho.clone();
            minus[e] -= h;
            let ap = mesh.homogenize(&plus, &material, true).expect("plus");
            let am = mesh.homogenize(&minus, &material, true).expect("minus");
            let fp = poisson_objective(&ap, target, w, a_ref).expect("fp").0;
            let fm = poisson_objective(&am, target, w, a_ref).expect("fm").0;
            let fd = (fp - fm) / (2.0 * h);
            assert!(
                (fd - grad[e]).abs() < 1e-5 * grad[e].abs().max(1e-4),
                "element {e}: fd={fd} analytic={}",
                grad[e]
            );
        }
    }

    #[test]
    fn density_filter_preserves_mean_and_is_self_adjoint() {
        let n = 7;
        let filter = PeriodicDensityFilter::new(n, 1.5);
        let x: Vec<f64> = (0..n * n).map(|k| ((k * 37) % 11) as f64 / 10.0).collect();
        let y: Vec<f64> = (0..n * n).map(|k| ((k * 13) % 7) as f64 / 6.0).collect();
        let hx = filter.apply(&x);
        let mean_x: f64 = x.iter().sum::<f64>() / x.len() as f64;
        let mean_hx: f64 = hx.iter().sum::<f64>() / hx.len() as f64;
        assert!((mean_x - mean_hx).abs() < 1e-12);
        let hty = filter.apply_transpose(&y);
        let lhs: f64 = hx.iter().zip(&y).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(&hty).map(|(a, b)| a * b).sum();
        assert!((lhs - rhs).abs() < 1e-12);
    }

    #[test]
    fn mma_solves_separable_constrained_quadratic() {
        // min Σ (x_j − a_j)²  s.t. mean(x) ≤ 0.3, 0 ≤ x ≤ 1.
        // KKT: x_j = clamp(a_j − μ, 0, 1) with mean = 0.3.
        let a = [0.9, 0.8, 0.5, 0.1];
        let mut x = vec![0.25; 4];
        let mut mma = MmaState::new(&x);
        for _ in 0..200 {
            let df0: Vec<f64> = x.iter().zip(&a).map(|(xi, ai)| 2.0 * (xi - ai)).collect();
            let g = x.iter().sum::<f64>() / 4.0 / 0.3 - 1.0;
            let dg = vec![1.0 / (4.0 * 0.3); 4];
            x = mma.update(&x, &df0, g, &dg, 0.2).expect("mma");
        }
        // μ solves Σ clamp(a − μ) = 1.2 → (0.9−μ)+(0.8−μ)+(0.5−μ) = 1.2 → μ = 1/3.
        let mu = 1.0 / 3.0;
        let expected = [0.9 - mu, 0.8 - mu, 0.5 - mu, 0.0];
        for (xi, ei) in x.iter().zip(&expected) {
            assert!((xi - ei).abs() < 2e-3, "x = {x:?}");
        }
    }

    #[test]
    fn symmetry_map_produces_d4_invariant_fields() {
        let n = 6;
        let map = SymmetryMap::new(n, CellSymmetry::Square);
        let mut v: Vec<f64> = (0..n * n).map(|k| k as f64).collect();
        map.symmetrize(&mut v);
        for i in 0..n {
            for j in 0..n {
                let val = v[i * n + j];
                assert_eq!(val, v[j * n + i]);
                assert_eq!(val, v[(n - 1 - i) * n + j]);
                assert_eq!(val, v[i * n + (n - 1 - j)]);
            }
        }
    }

    #[test]
    fn invalid_configurations_are_rejected() {
        let bad_target = fast_config(16, -1.2);
        assert!(optimize_auxetic_topology(&bad_target).is_err());
        let bad_volume = AuxeticTopologyConfig {
            volume_fraction: 0.0,
            ..fast_config(16, -0.3)
        };
        assert!(optimize_auxetic_topology(&bad_volume).is_err());
        let robust_without_projection = AuxeticTopologyConfig {
            projection_beta_max: None,
            ..fast_config(16, -0.3)
        };
        assert!(optimize_auxetic_topology(&robust_without_projection).is_err());
        let tiny = fast_config(2, -0.3);
        assert!(optimize_auxetic_topology(&tiny).is_err());
        assert!(topology_optimized_auxetic(f64::NAN, 10).is_err());
        let ragged = vec![vec![1.0; 4], vec![1.0; 3], vec![1.0; 4], vec![1.0; 4]];
        assert!(homogenize_periodic_plane_stress(&ragged, &SimpMaterial::default()).is_err());
    }

    #[test]
    fn zero_iterations_returns_initial_design() {
        let config = AuxeticTopologyConfig {
            max_iterations: 0,
            ..fast_config(12, -0.3)
        };
        let res = optimize_auxetic_topology(&config).expect("evaluation");
        assert_eq!(res.iterations, 0);
        assert!(res.objective_history.is_empty());
        assert_eq!(res.layout.len(), 12);
        // Initial grey field (ρ = 0.4 with a 0.2 central disc) lies below the
        // 0.5 threshold everywhere.
        assert!(res.layout.iter().flatten().all(|&s| !s));
        let grid = topology_optimized_auxetic(-0.3, 0).expect("wrapper");
        assert_eq!(grid.len(), 32);
        assert!(grid.iter().all(|row| row.len() == 32));
    }

    #[test]
    fn negative_target_yields_auxetic_binary_layout() {
        let res = optimize_auxetic_topology(&fast_config(24, -0.3)).expect("optimisation");
        assert!(res.converged, "iterations = {}", res.iterations);
        assert!(
            res.effective_poisson < -0.2,
            "grey ν = {}",
            res.effective_poisson
        );
        // The thresholded binary layout (void = E_min) is re-homogenised.
        assert!(
            res.layout_poisson < -0.15,
            "layout ν = {}",
            res.layout_poisson
        );
        assert!(res.layout_stiffness.c[0][1] < 0.0);
        // Connected load paths: stiffness far above the void modulus.
        assert!(res.layout_stiffness.c[0][0] > 1e-3);
        assert!((res.layout_volume_fraction - 0.4).abs() < 0.05);
        let solid = res.layout.iter().flatten().filter(|&&s| s).count();
        assert!(solid > 0 && solid < 24 * 24);
        // Objective decreased from the initial design.
        let first = res.objective_history.first().copied().unwrap_or(f64::NAN);
        let last = res.objective_history.last().copied().unwrap_or(f64::NAN);
        assert!(last < first, "objective {first} -> {last}");
    }

    #[test]
    fn positive_target_differs_from_auxetic_design() {
        let pos = optimize_auxetic_topology(&fast_config(20, 0.3)).expect("positive");
        assert!(
            pos.layout_poisson > 0.2,
            "layout ν = {}",
            pos.layout_poisson
        );
        assert!((pos.effective_poisson - 0.3).abs() < 0.05);
        let neg = optimize_auxetic_topology(&fast_config(20, -0.1)).expect("negative");
        assert!(
            neg.layout_poisson < 0.0,
            "layout ν = {}",
            neg.layout_poisson
        );
        assert_ne!(pos.layout, neg.layout);
    }
}
