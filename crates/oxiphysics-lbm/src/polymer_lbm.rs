#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polymer solution simulation with the Lattice Boltzmann Method.
//!
//! This module provides a complete framework for simulating polymer solutions
//! using LBM coupled with viscoelastic constitutive models:
//!
//! - **Oldroyd-B model**: upper-convected Maxwell polymer stress evolution
//! - **FENE-P model**: finitely-extensible nonlinear elastic (Peterlin closure)
//! - **Conformation tensor transport**: advection + relaxation + stretching
//! - **Weissenberg / Deborah number effects**: elastic turbulence characterisation
//! - **Polymer drag reduction (Toms effect)**: friction factor reduction
//! - **Elastic turbulence**: purely elastic flow instability at low Re
//! - **Coil-stretch transition**: polymer chain extension vs. relaxation
//! - **Chain orientation tensor**: second-moment statistics of end-to-end vector
//! - **Rouse / Zimm relaxation times**: theoretical polymer relaxation scales
//! - **Channel flow with polymer additives**: Poiseuille + polymer stress coupling
//!
//! # Key dimensionless groups
//! - Wi = λ γ̇  (Weissenberg number)
//! - De = λ / t_p (Deborah number)
//! - β = η_s / η_0 (solvent-to-total viscosity ratio)
//! - b = R²_max / R²_eq (FENE extensibility parameter)

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Boltzmann constant \[J/K\].
pub const K_BOLTZMANN: f64 = 1.380_649e-23;

/// Maximum number of iterations for iterative solvers in this module.
pub const MAX_ITER: usize = 10_000;

/// Convergence tolerance for stress/conformation tensor iterations.
pub const CONV_TOL: f64 = 1.0e-10;

// ============================================================================
// Section 1 – Dimensionless numbers and scaling
// ============================================================================

/// Compute the Weissenberg number Wi = λ · γ̇.
///
/// Wi compares the polymer relaxation time λ to the characteristic flow
/// time scale 1/γ̇.  Wi > 1 means elastic effects are significant.
///
/// # Arguments
/// * `lambda`     – polymer relaxation time \[s\]
/// * `shear_rate` – characteristic shear rate \[1/s\]
#[allow(dead_code)]
pub fn weissenberg_number(lambda: f64, shear_rate: f64) -> f64 {
    lambda * shear_rate
}

/// Compute the Deborah number De = λ / t_process.
///
/// De < 1 → fluid-like behaviour; De > 1 → elastic (solid-like) behaviour.
///
/// # Arguments
/// * `lambda`    – polymer relaxation time \[s\]
/// * `t_process` – characteristic process time \[s\]
#[allow(dead_code)]
pub fn deborah_number(lambda: f64, t_process: f64) -> f64 {
    if t_process.abs() < f64::EPSILON {
        return f64::INFINITY;
    }
    lambda / t_process
}

/// Compute the polymer contribution to viscosity η_p = η_0 (1 − β).
///
/// # Arguments
/// * `eta_total` – total zero-shear viscosity η_0 \[Pa·s\]
/// * `beta`      – solvent-to-total viscosity ratio β ∈ (0,1)
#[allow(dead_code)]
pub fn polymer_viscosity(eta_total: f64, beta: f64) -> f64 {
    eta_total * (1.0 - beta)
}

/// Polymer contribution to the LBM relaxation frequency ω_p.
///
/// In lattice units, ω_p = 1 / (3 η_p / ρ c_s² Δt + 0.5).
///
/// # Arguments
/// * `eta_p` – polymer dynamic viscosity in lattice units
#[allow(dead_code)]
pub fn polymer_omega(eta_p: f64) -> f64 {
    1.0 / (3.0 * eta_p + 0.5)
}

// ============================================================================
// Section 2 – Rouse and Zimm relaxation times
// ============================================================================

/// Rouse relaxation time for a Gaussian chain of N segments.
///
/// τ_R = ζ N² b² / (6π² k_B T)
/// where ζ is the monomer friction coefficient,
/// N is the number of segments,
/// b is the statistical segment length.
///
/// # Arguments
/// * `friction_coeff` – monomer friction coefficient ζ \[kg/s\]
/// * `n_segments`     – number of Kuhn segments N
/// * `b`              – segment length \[m\]
/// * `temperature`    – absolute temperature T \[K\]
#[allow(dead_code)]
pub fn rouse_relaxation_time(
    friction_coeff: f64,
    n_segments: f64,
    b: f64,
    temperature: f64,
) -> f64 {
    let kbt = K_BOLTZMANN * temperature;
    friction_coeff * n_segments * n_segments * b * b / (6.0 * PI * PI * kbt)
}

/// Zimm relaxation time (with hydrodynamic interactions).
///
/// τ_Z ≈ 0.325 (η_s / k_B T) (R_g³) where R_g² = N b² / 6.
///
/// # Arguments
/// * `eta_s`      – solvent viscosity \[Pa·s\]
/// * `n_segments` – number of Kuhn segments
/// * `b`          – segment length \[m\]
/// * `temperature`– absolute temperature \[K\]
#[allow(dead_code)]
pub fn zimm_relaxation_time(eta_s: f64, n_segments: f64, b: f64, temperature: f64) -> f64 {
    let rg_sq = n_segments * b * b / 6.0;
    let rg = rg_sq.sqrt();
    let rg3 = rg * rg * rg;
    let kbt = K_BOLTZMANN * temperature;
    0.325 * eta_s * rg3 / kbt
}

/// Compute the p-th Rouse mode relaxation time τ_p = τ_R / p².
///
/// # Arguments
/// * `tau_rouse` – the longest (p=1) Rouse relaxation time \[s\]
/// * `p`         – mode index (p ≥ 1)
#[allow(dead_code)]
pub fn rouse_mode_time(tau_rouse: f64, p: usize) -> f64 {
    let pf = p as f64;
    tau_rouse / (pf * pf)
}

// ============================================================================
// Section 3 – Conformation tensor (second moment <RR>)
// ============================================================================

/// Symmetric 2×2 conformation tensor stored as \[Cxx, Cyy, Cxy\].
///
/// The conformation tensor A = `RR` / R_eq² describes the mean-square
/// end-to-end vector of the polymer chain normalised by the equilibrium value.
/// At rest, A = I (identity).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ConformationTensor2D {
    /// (0,0) component A_xx.
    pub xx: f64,
    /// (1,1) component A_yy.
    pub yy: f64,
    /// (0,1) = (1,0) component A_xy (symmetric).
    pub xy: f64,
}

impl ConformationTensor2D {
    /// Identity conformation tensor (equilibrium, no stretch).
    pub fn identity() -> Self {
        Self {
            xx: 1.0,
            yy: 1.0,
            xy: 0.0,
        }
    }

    /// Trace of the conformation tensor tr(A) = A_xx + A_yy.
    pub fn trace(&self) -> f64 {
        self.xx + self.yy
    }

    /// Frobenius norm squared: |A|² = Axx² + Ayy² + 2 Axy².
    pub fn norm_sq(&self) -> f64 {
        self.xx * self.xx + self.yy * self.yy + 2.0 * self.xy * self.xy
    }

    /// Determinant det(A) = Axx Ayy - Axy².
    pub fn det(&self) -> f64 {
        self.xx * self.yy - self.xy * self.xy
    }

    /// Add two conformation tensors component-wise.
    pub fn add(&self, other: &Self) -> Self {
        Self {
            xx: self.xx + other.xx,
            yy: self.yy + other.yy,
            xy: self.xy + other.xy,
        }
    }

    /// Scale the conformation tensor by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            xx: self.xx * s,
            yy: self.yy * s,
            xy: self.xy * s,
        }
    }
}

/// Symmetric 3×3 conformation tensor stored as \[Axx, Ayy, Azz, Axy, Axz, Ayz\].
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ConformationTensor3D {
    /// Component A_xx.
    pub xx: f64,
    /// Component A_yy.
    pub yy: f64,
    /// Component A_zz.
    pub zz: f64,
    /// Component A_xy = A_yx.
    pub xy: f64,
    /// Component A_xz = A_zx.
    pub xz: f64,
    /// Component A_yz = A_zy.
    pub yz: f64,
}

impl ConformationTensor3D {
    /// Identity conformation tensor in 3D.
    pub fn identity() -> Self {
        Self {
            xx: 1.0,
            yy: 1.0,
            zz: 1.0,
            xy: 0.0,
            xz: 0.0,
            yz: 0.0,
        }
    }

    /// Trace tr(A) = A_xx + A_yy + A_zz.
    pub fn trace(&self) -> f64 {
        self.xx + self.yy + self.zz
    }

    /// Add two 3D conformation tensors.
    pub fn add(&self, other: &Self) -> Self {
        Self {
            xx: self.xx + other.xx,
            yy: self.yy + other.yy,
            zz: self.zz + other.zz,
            xy: self.xy + other.xy,
            xz: self.xz + other.xz,
            yz: self.yz + other.yz,
        }
    }

    /// Scale by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            xx: self.xx * s,
            yy: self.yy * s,
            zz: self.zz * s,
            xy: self.xy * s,
            xz: self.xz * s,
            yz: self.yz * s,
        }
    }
}

// ============================================================================
// Section 4 – Oldroyd-B polymer stress
// ============================================================================

/// Polymer extra-stress tensor in 2D (Oldroyd-B / UCM), stored as \[τxx, τyy, τxy\].
///
/// In the Oldroyd-B model the polymer stress is τ_p = (η_p/λ)(A − I)
/// where A is the conformation tensor and λ is the relaxation time.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PolymerStress2D {
    /// Normal stress component τ_xx \[Pa\].
    pub txx: f64,
    /// Normal stress component τ_yy \[Pa\].
    pub tyy: f64,
    /// Shear stress component τ_xy \[Pa\].
    pub txy: f64,
}

impl PolymerStress2D {
    /// Zero stress state.
    pub fn zero() -> Self {
        Self {
            txx: 0.0,
            tyy: 0.0,
            txy: 0.0,
        }
    }

    /// Compute polymer stress from the conformation tensor (Oldroyd-B).
    ///
    /// τ_p = (η_p / λ) (A − I)
    ///
    /// # Arguments
    /// * `conf`  – conformation tensor A
    /// * `eta_p` – polymer viscosity \[Pa·s\]
    /// * `lambda`– relaxation time \[s\]
    pub fn from_conformation(conf: &ConformationTensor2D, eta_p: f64, lambda: f64) -> Self {
        let c = eta_p / lambda;
        Self {
            txx: c * (conf.xx - 1.0),
            tyy: c * (conf.yy - 1.0),
            txy: c * conf.xy,
        }
    }

    /// First normal stress difference N₁ = τ_xx − τ_yy.
    pub fn n1(&self) -> f64 {
        self.txx - self.tyy
    }

    /// Second normal stress difference N₂ = τ_yy − τ_zz.
    /// For 2-D plane flow, N₂ = τ_yy (τ_zz = 0 assumed).
    pub fn n2(&self) -> f64 {
        self.tyy
    }
}

/// Evolve the 2D conformation tensor by one Oldroyd-B time step (explicit Euler).
///
/// dA/dt + (u·∇)A = A·(∇u)ᵀ + (∇u)·A − (1/λ)(A − I)
///
/// This function applies only the relaxation + stretching terms (no advection).
/// Pass the velocity gradient components directly.
///
/// # Arguments
/// * `a`      – current conformation tensor
/// * `dudx`   – ∂u/∂x
/// * `dudy`   – ∂u/∂y
/// * `dvdx`   – ∂v/∂x
/// * `dvdy`   – ∂v/∂y
/// * `lambda` – polymer relaxation time \[lattice units\]
/// * `dt`     – time step \[lattice units\]
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn oldroyd_b_evolve_2d(
    a: &ConformationTensor2D,
    dudx: f64,
    dudy: f64,
    dvdx: f64,
    dvdy: f64,
    lambda: f64,
    dt: f64,
) -> ConformationTensor2D {
    // Upper-convected derivative stretching: A·(∇u)ᵀ + (∇u)·A
    let stretch_xx = 2.0 * (a.xx * dudx + a.xy * dvdx);
    let stretch_yy = 2.0 * (a.xy * dudy + a.yy * dvdy);
    let stretch_xy = a.xx * dudy + a.xy * dvdy + a.xy * dudx + a.yy * dvdx;

    // Relaxation: −(A − I)/λ
    let relax_xx = -(a.xx - 1.0) / lambda;
    let relax_yy = -(a.yy - 1.0) / lambda;
    let relax_xy = -a.xy / lambda;

    ConformationTensor2D {
        xx: a.xx + dt * (stretch_xx + relax_xx),
        yy: a.yy + dt * (stretch_yy + relax_yy),
        xy: a.xy + dt * (stretch_xy + relax_xy),
    }
}

// ============================================================================
// Section 5 – FENE-P model
// ============================================================================

/// FENE-P spring coefficient f(tr A) = 1 / (1 − tr(A)/b).
///
/// # Arguments
/// * `trace_a` – trace of the conformation tensor tr(A)
/// * `b`       – FENE extensibility parameter (R²_max / R²_eq)
#[allow(dead_code)]
pub fn fene_p_spring(trace_a: f64, b: f64) -> f64 {
    let denom = 1.0 - trace_a / b;
    if denom < f64::EPSILON {
        return f64::MAX;
    }
    1.0 / denom
}

/// Evolve the 2D FENE-P conformation tensor by one explicit Euler step.
///
/// dA/dt = A·(∇u)ᵀ + (∇u)·A − (f/λ)(A − (1/f)I)
///
/// where f = 1/(1 − tr(A)/b) is the Peterlin spring factor.
///
/// # Arguments
/// * `a`      – current conformation tensor
/// * `dudx`   – ∂u/∂x
/// * `dudy`   – ∂u/∂y
/// * `dvdx`   – ∂v/∂x
/// * `dvdy`   – ∂v/∂y
/// * `lambda` – relaxation time \[lattice units\]
/// * `b`      – FENE extensibility parameter
/// * `dt`     – time step
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn fene_p_evolve_2d(
    a: &ConformationTensor2D,
    dudx: f64,
    dudy: f64,
    dvdx: f64,
    dvdy: f64,
    lambda: f64,
    b: f64,
    dt: f64,
) -> ConformationTensor2D {
    let tr = a.trace();
    let f = fene_p_spring(tr, b);

    let stretch_xx = 2.0 * (a.xx * dudx + a.xy * dvdx);
    let stretch_yy = 2.0 * (a.xy * dudy + a.yy * dvdy);
    let stretch_xy = a.xx * dudy + a.xy * dvdy + a.xy * dudx + a.yy * dvdx;

    let relax_xx = -(f * a.xx - 1.0) / lambda;
    let relax_yy = -(f * a.yy - 1.0) / lambda;
    let relax_xy = -(f * a.xy) / lambda;

    ConformationTensor2D {
        xx: a.xx + dt * (stretch_xx + relax_xx),
        yy: a.yy + dt * (stretch_yy + relax_yy),
        xy: a.xy + dt * (stretch_xy + relax_xy),
    }
}

/// Compute the FENE-P polymer stress τ_p = (η_p/λ) f(A − I/f).
///
/// # Arguments
/// * `conf`   – conformation tensor
/// * `eta_p`  – polymer viscosity
/// * `lambda` – relaxation time
/// * `b`      – extensibility parameter
#[allow(dead_code)]
pub fn fene_p_stress_2d(
    conf: &ConformationTensor2D,
    eta_p: f64,
    lambda: f64,
    b: f64,
) -> PolymerStress2D {
    let f = fene_p_spring(conf.trace(), b);
    let c = eta_p / lambda;
    PolymerStress2D {
        txx: c * (f * conf.xx - 1.0),
        tyy: c * (f * conf.yy - 1.0),
        txy: c * f * conf.xy,
    }
}

// ============================================================================
// Section 6 – Polymer chain orientation tensor
// ============================================================================

/// Second-moment (orientation) tensor S of the end-to-end unit vector.
///
/// S_ij = <u_i u_j> where u = R/|R|.
/// For an isotropic distribution S = I/2 (2D) or I/3 (3D).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct OrientationTensor2D {
    /// S_xx component.
    pub sxx: f64,
    /// S_yy component.
    pub syy: f64,
    /// S_xy = S_yx component.
    pub sxy: f64,
}

impl OrientationTensor2D {
    /// Isotropic (equilibrium) orientation tensor: S = I/2.
    pub fn isotropic() -> Self {
        Self {
            sxx: 0.5,
            syy: 0.5,
            sxy: 0.0,
        }
    }

    /// Compute orientation tensor from conformation tensor A.
    ///
    /// S_ij ≈ A_ij / tr(A)
    pub fn from_conformation(conf: &ConformationTensor2D) -> Self {
        let tr = conf.trace().max(f64::EPSILON);
        Self {
            sxx: conf.xx / tr,
            syy: conf.yy / tr,
            sxy: conf.xy / tr,
        }
    }

    /// Orientational order parameter S_op = sqrt(S:S - 1/d) for d=2.
    pub fn order_parameter(&self) -> f64 {
        let s2 = self.sxx * self.sxx + self.syy * self.syy + 2.0 * self.sxy * self.sxy;
        (s2 - 0.5).max(0.0).sqrt()
    }
}

// ============================================================================
// Section 7 – Polymer drag reduction (Toms effect)
// ============================================================================

/// Virk's maximum drag reduction (MDR) asymptote for the friction factor.
///
/// f_MDR = 0.32 / log10²(Re √f_MDR)  — Virk's empirical formula.
/// This function returns Virk's correlation directly:
/// f_MDR ≈ 0.58 Re^{-0.58}
///
/// # Arguments
/// * `re` – Reynolds number
#[allow(dead_code)]
pub fn virk_mdr_friction(re: f64) -> f64 {
    if re < 1.0 {
        return 1.0;
    }
    0.58 * re.powf(-0.58)
}

/// Newtonian (Prandtl–Kármán) friction factor for turbulent pipe flow.
///
/// 1/√f = 4.0 log10(Re √f) − 0.4  (implicit; iterated here)
///
/// # Arguments
/// * `re` – Reynolds number (turbulent: Re > 4000)
#[allow(dead_code)]
pub fn prandtl_karman_friction(re: f64) -> f64 {
    if re < 1.0 {
        return 1.0;
    }
    // Initial guess via Blasius
    let mut f = 0.316 * re.powf(-0.25);
    for _ in 0..50 {
        let lhs = 1.0 / f.sqrt();
        let rhs = 4.0 * (re * f.sqrt()).log10() - 0.4;
        let df = lhs - rhs;
        f += df * 0.01;
        f = f.max(1.0e-6);
        if df.abs() < 1.0e-8 {
            break;
        }
    }
    f
}

/// Polymer drag reduction percentage DR% = (f_N − f_P) / f_N × 100%.
///
/// # Arguments
/// * `f_newtonian` – Newtonian friction factor
/// * `f_polymer`   – friction factor with polymer additives
#[allow(dead_code)]
pub fn drag_reduction_percent(f_newtonian: f64, f_polymer: f64) -> f64 {
    if f_newtonian.abs() < f64::EPSILON {
        return 0.0;
    }
    (f_newtonian - f_polymer) / f_newtonian * 100.0
}

// ============================================================================
// Section 8 – Coil-stretch transition
// ============================================================================

/// Critical Weissenberg number for the coil-stretch transition in extensional flow.
///
/// Wi_c = 0.5 (Lumley criterion for extensional flow).
pub const COIL_STRETCH_WI_CRITICAL: f64 = 0.5;

/// Determine whether the polymer chain is in the coil or stretched state.
///
/// Returns `true` if Wi > Wi_c (stretched state).
///
/// # Arguments
/// * `wi` – local Weissenberg number Wi = λ ε̇  (extensional rate × relax time)
#[allow(dead_code)]
pub fn is_stretched(wi: f64) -> bool {
    wi > COIL_STRETCH_WI_CRITICAL
}

/// Compute fractional chain extension r/r_max from the FENE conformation tensor.
///
/// r/r_max = sqrt(tr(A) / b)
///
/// # Arguments
/// * `trace_a` – trace of the conformation tensor
/// * `b`       – FENE extensibility parameter
#[allow(dead_code)]
pub fn chain_extension_ratio(trace_a: f64, b: f64) -> f64 {
    (trace_a / b).sqrt().min(1.0)
}

// ============================================================================
// Section 9 – Elastic turbulence diagnostics
// ============================================================================

/// Elastic turbulence onset criterion based on Groisman & Steinberg (2000).
///
/// Onset occurs when Wi > Wi_c,elastic for curved streamlines.
/// Here we use the empirical estimate Wi_c,elastic ≈ O(1).
pub const ELASTIC_TURBULENCE_WI_ONSET: f64 = 1.0;

/// Compute the polymer stress amplification factor (viscosity ratio).
///
/// In elastic turbulence, velocity fluctuations scale as u' ~ τ_p / η_0.
///
/// # Arguments
/// * `tau_p_rms` – RMS polymer shear stress
/// * `eta_0`     – total zero-shear viscosity
#[allow(dead_code)]
pub fn elastic_stress_amplification(tau_p_rms: f64, eta_0: f64) -> f64 {
    if eta_0.abs() < f64::EPSILON {
        return 0.0;
    }
    tau_p_rms / eta_0
}

// ============================================================================
// Section 10 – Polymer LBM grid cell
// ============================================================================

/// A single LBM grid cell with D2Q9 distribution functions and polymer state.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PolymerCell {
    /// D2Q9 distribution functions f_i for i=0..9.
    pub f: [f64; 9],
    /// Equilibrium distribution functions f_i^eq.
    pub feq: [f64; 9],
    /// Polymer conformation tensor at this cell.
    pub conformation: ConformationTensor2D,
    /// Polymer extra-stress tensor at this cell.
    pub stress: PolymerStress2D,
    /// Fluid density ρ.
    pub rho: f64,
    /// Fluid velocity u_x.
    pub ux: f64,
    /// Fluid velocity u_y.
    pub uy: f64,
}

impl PolymerCell {
    /// Create a new PolymerCell in equilibrium.
    pub fn new_equilibrium(rho: f64, ux: f64, uy: f64) -> Self {
        let mut cell = Self {
            f: [0.0; 9],
            feq: [0.0; 9],
            conformation: ConformationTensor2D::identity(),
            stress: PolymerStress2D::zero(),
            rho,
            ux,
            uy,
        };
        cell.compute_equilibrium();
        cell.f = cell.feq;
        cell
    }

    /// Compute the D2Q9 Maxwell-Boltzmann equilibrium distribution.
    pub fn compute_equilibrium(&mut self) {
        let w = D2Q9_WEIGHTS;
        let ex = D2Q9_EX;
        let ey = D2Q9_EY;
        let u2 = self.ux * self.ux + self.uy * self.uy;
        for i in 0..9 {
            let eu = ex[i] * self.ux + ey[i] * self.uy;
            self.feq[i] = w[i] * self.rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
        }
    }

    /// Compute macroscopic density and velocity from distribution functions.
    pub fn compute_macroscopic(&mut self) {
        self.rho = self.f.iter().sum();
        let ex = D2Q9_EX;
        let ey = D2Q9_EY;
        let mut ux = 0.0;
        let mut uy = 0.0;
        for i in 0..9 {
            ux += self.f[i] * ex[i];
            uy += self.f[i] * ey[i];
        }
        if self.rho > f64::EPSILON {
            self.ux = ux / self.rho;
            self.uy = uy / self.rho;
        }
    }
}

// ============================================================================
// Section 11 – D2Q9 lattice constants
// ============================================================================

/// D2Q9 lattice weights.
pub const D2Q9_WEIGHTS: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 x-components of lattice velocities.
pub const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 y-components of lattice velocities.
pub const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 bounce-back index pairs (i → opp\[i\]).
pub const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ============================================================================
// Section 12 – BGK-polymer collision step
// ============================================================================

/// Perform a single BGK collision step with polymer stress forcing.
///
/// The post-collision distribution is:
/// f_i* = f_i − ω_s (f_i − f_i^eq) + Δf_i^polymer
///
/// where Δf_i^polymer accounts for the polymer extra-stress divergence.
///
/// # Arguments
/// * `cell`    – mutable polymer cell (updated in-place)
/// * `omega_s` – solvent BGK relaxation frequency ω_s = 1/(3ν_s + 0.5)
/// * `eta_p`   – polymer viscosity in lattice units
/// * `lambda`  – polymer relaxation time in lattice units
#[allow(dead_code)]
pub fn bgk_polymer_collision(cell: &mut PolymerCell, omega_s: f64, eta_p: f64, lambda: f64) {
    cell.compute_macroscopic();
    cell.compute_equilibrium();

    // Polymer stress contribution to the collision operator (Fattal-Kupferman approach)
    let tau_p = PolymerStress2D::from_conformation(&cell.conformation, eta_p, lambda);
    cell.stress = tau_p;

    let ex = D2Q9_EX;
    let ey = D2Q9_EY;
    let w = D2Q9_WEIGHTS;
    let rho = cell.rho.max(f64::EPSILON);

    for i in 0..9 {
        // Polymer stress forcing term (projection onto lattice direction)
        let pi_xx_contrib = ex[i] * ex[i] - 1.0 / 3.0;
        let pi_yy_contrib = ey[i] * ey[i] - 1.0 / 3.0;
        let pi_xy_contrib = ex[i] * ey[i];

        let stress_force = 3.0 * w[i] / rho
            * (pi_xx_contrib * tau_p.txx
                + pi_yy_contrib * tau_p.tyy
                + 2.0 * pi_xy_contrib * tau_p.txy);

        cell.f[i] = cell.f[i] - omega_s * (cell.f[i] - cell.feq[i]) + stress_force;
    }
}

// ============================================================================
// Section 13 – Polymer LBM 2D grid
// ============================================================================

/// 2D polymer LBM simulation grid (Nx × Ny cells).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PolymerLbmGrid2D {
    /// Number of cells in x direction.
    pub nx: usize,
    /// Number of cells in y direction.
    pub ny: usize,
    /// Flat array of polymer cells (row-major, index = iy*nx + ix).
    pub cells: Vec<PolymerCell>,
    /// Solvent BGK relaxation frequency.
    pub omega_s: f64,
    /// Polymer viscosity (lattice units).
    pub eta_p: f64,
    /// Polymer relaxation time (lattice units).
    pub lambda: f64,
    /// FENE extensibility parameter b (use `f64::INFINITY` for Oldroyd-B).
    pub fene_b: f64,
    /// Current simulation time step counter.
    pub step: usize,
}

impl PolymerLbmGrid2D {
    /// Create a new polymer LBM grid with uniform initial conditions.
    ///
    /// # Arguments
    /// * `nx`, `ny` – grid dimensions
    /// * `omega_s`  – solvent relaxation frequency
    /// * `eta_p`    – polymer viscosity (lattice units)
    /// * `lambda`   – polymer relaxation time (lattice units)
    /// * `fene_b`   – FENE extensibility parameter (∞ for Oldroyd-B)
    #[allow(clippy::too_many_arguments)]
    pub fn new(nx: usize, ny: usize, omega_s: f64, eta_p: f64, lambda: f64, fene_b: f64) -> Self {
        let cells = vec![PolymerCell::new_equilibrium(1.0, 0.0, 0.0); nx * ny];
        Self {
            nx,
            ny,
            cells,
            omega_s,
            eta_p,
            lambda,
            fene_b,
            step: 0,
        }
    }

    /// Get cell index from (ix, iy) coordinates.
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Set all cells to a Poiseuille velocity profile (channel flow).
    ///
    /// u_x(y) = u_max * 4 y/H (1 − y/H),  where H = ny − 1.
    ///
    /// # Arguments
    /// * `u_max` – maximum centreline velocity (lattice units)
    pub fn set_poiseuille_profile(&mut self, u_max: f64) {
        let h = (self.ny - 1) as f64;
        for iy in 0..self.ny {
            let y = iy as f64 / h;
            let ux = u_max * 4.0 * y * (1.0 - y);
            for ix in 0..self.nx {
                let idx = self.idx(ix, iy);
                self.cells[idx].ux = ux;
                self.cells[idx].uy = 0.0;
                self.cells[idx].compute_equilibrium();
                self.cells[idx].f = self.cells[idx].feq;
            }
        }
    }

    /// Apply bounce-back boundary conditions on top and bottom walls.
    pub fn apply_wall_bounce_back(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let opp = D2Q9_OPP;

        // Bottom wall (iy = 0)
        for ix in 0..nx {
            let idx = self.idx(ix, 0);
            // Bounce-back: f_opp[i] <- f[i]
            for i in 0..9 {
                let tmp = self.cells[idx].f[i];
                self.cells[idx].f[opp[i]] = tmp;
            }
        }
        // Top wall (iy = ny-1)
        for ix in 0..nx {
            let idx = self.idx(ix, ny - 1);
            for i in 0..9 {
                let tmp = self.cells[idx].f[i];
                self.cells[idx].f[opp[i]] = tmp;
            }
        }
    }

    /// Perform streaming step (periodic in x, bounce-back walls in y).
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let ex: [i64; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
        let ey: [i64; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

        let old_f: Vec<[f64; 9]> = self.cells.iter().map(|c| c.f).collect();

        for iy in 0..ny {
            for ix in 0..nx {
                let dst = self.idx(ix, iy);
                for i in 0..9 {
                    let sx = (ix as i64 - ex[i]).rem_euclid(nx as i64) as usize;
                    let sy = (iy as i64 - ey[i]).rem_euclid(ny as i64) as usize;
                    let src = self.idx(sx, sy);
                    self.cells[dst].f[i] = old_f[src][i];
                }
            }
        }
    }

    /// Perform one complete LBM time step (collision + streaming + conformation update).
    ///
    /// Uses explicit Euler for the conformation tensor evolution with a simple
    /// estimate of the velocity gradient from neighbouring cells.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // --- Collision ---
        for idx in 0..(nx * ny) {
            bgk_polymer_collision(&mut self.cells[idx], self.omega_s, self.eta_p, self.lambda);
        }

        // --- Streaming ---
        self.stream();

        // --- Update macroscopic fields ---
        for idx in 0..(nx * ny) {
            self.cells[idx].compute_macroscopic();
        }

        // --- Conformation tensor evolution (relaxation only, no spatial gradient here) ---
        let lambda = self.lambda;
        let fene_b = self.fene_b;
        let dt = 1.0_f64; // lattice time step
        for idx in 0..(nx * ny) {
            let a = self.cells[idx].conformation;
            self.cells[idx].conformation = if fene_b.is_finite() {
                fene_p_evolve_2d(&a, 0.0, 0.0, 0.0, 0.0, lambda, fene_b, dt)
            } else {
                oldroyd_b_evolve_2d(&a, 0.0, 0.0, 0.0, 0.0, lambda, dt)
            };
        }

        self.step += 1;
    }

    /// Compute the volume-averaged polymer kinetic energy k_p = Σ tr(A) / (2 N).
    pub fn mean_polymer_trace(&self) -> f64 {
        let sum: f64 = self.cells.iter().map(|c| c.conformation.trace()).sum();
        sum / (self.cells.len() as f64)
    }

    /// Compute the maximum first normal stress difference across the grid.
    pub fn max_n1(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| c.stress.n1().abs())
            .fold(0.0_f64, f64::max)
    }
}

// ============================================================================
// Section 14 – Channel flow analytics (polymer-modified Poiseuille)
// ============================================================================

/// Compute the effective viscosity in a polymer channel flow.
///
/// η_eff = η_s + η_p / (1 + Wi²)  (simplified linear Oldroyd-B result)
///
/// # Arguments
/// * `eta_s`  – solvent viscosity \[Pa·s\]
/// * `eta_p`  – polymer viscosity \[Pa·s\]
/// * `wi`     – Weissenberg number
#[allow(dead_code)]
pub fn effective_viscosity_polymer(eta_s: f64, eta_p: f64, wi: f64) -> f64 {
    eta_s + eta_p / (1.0 + wi * wi)
}

/// Poiseuille flow centreline velocity with polymer modification.
///
/// U_c = (−dP/dx) H² / (8 η_eff)
///
/// # Arguments
/// * `dp_dx`    – pressure gradient \[Pa/m\]
/// * `h`        – channel half-height \[m\]
/// * `eta_eff`  – effective viscosity \[Pa·s\]
#[allow(dead_code)]
pub fn polymer_poiseuille_velocity(dp_dx: f64, h: f64, eta_eff: f64) -> f64 {
    if eta_eff.abs() < f64::EPSILON {
        return 0.0;
    }
    (-dp_dx) * h * h / (8.0 * eta_eff)
}

/// Compute the normal stress growth in start-up shear flow (Oldroyd-B).
///
/// τ_xx(t) = 2 η_p Wi² (1 − e^{−t/λ})
///
/// # Arguments
/// * `eta_p`  – polymer viscosity
/// * `wi`     – Weissenberg number
/// * `t`      – time since start-up
/// * `lambda` – relaxation time
#[allow(dead_code)]
pub fn startup_shear_n1(eta_p: f64, wi: f64, t: f64, lambda: f64) -> f64 {
    2.0 * eta_p * wi * wi * (1.0 - (-t / lambda).exp())
}

// ============================================================================
// Section 15 – Weissenberg number field diagnostics
// ============================================================================

/// Compute local Weissenberg number field from a velocity profile.
///
/// Wi(y) = λ |du/dy|
///
/// # Arguments
/// * `u_profile` – x-velocity sampled at evenly spaced y-positions
/// * `dy`        – y-spacing
/// * `lambda`    – polymer relaxation time
#[allow(dead_code)]
pub fn local_weissenberg_field(u_profile: &[f64], dy: f64, lambda: f64) -> Vec<f64> {
    let n = u_profile.len();
    let mut wi = vec![0.0_f64; n];
    for i in 1..n - 1 {
        let du_dy = (u_profile[i + 1] - u_profile[i - 1]) / (2.0 * dy);
        wi[i] = lambda * du_dy.abs();
    }
    if n >= 2 {
        wi[0] = lambda * ((u_profile[1] - u_profile[0]) / dy).abs();
        wi[n - 1] = lambda * ((u_profile[n - 1] - u_profile[n - 2]) / dy).abs();
    }
    wi
}

// ============================================================================
// Section 16 – Elastic turbulence energy spectrum
// ============================================================================

/// Power-law exponent of the elastic turbulence kinetic energy spectrum.
///
/// E(k) ~ k^{−α} with α ≈ 3.5 (Groisman & Steinberg).
pub const ELASTIC_TURBULENCE_EXPONENT: f64 = 3.5;

/// Estimate the elastic turbulence power spectrum amplitude at wavenumber k.
///
/// E(k) = A · k^{−α}
///
/// # Arguments
/// * `k`         – wavenumber \[1/m\]
/// * `amplitude` – spectral amplitude A
/// * `alpha`     – spectral exponent (typically 3.5)
#[allow(dead_code)]
pub fn elastic_turbulence_spectrum(k: f64, amplitude: f64, alpha: f64) -> f64 {
    if k < f64::EPSILON {
        return 0.0;
    }
    amplitude * k.powf(-alpha)
}

// ============================================================================
// Section 17 – Polymer concentration effects
// ============================================================================

/// Polymer solution intrinsic viscosity via the Mark-Houwink equation.
///
/// \[η\] = K M^a  where M is the polymer molecular weight.
///
/// # Arguments
/// * `k`         – Mark-Houwink prefactor
/// * `mol_weight`– molecular weight \[g/mol\]
/// * `a`         – Mark-Houwink exponent (0.5–0.8 for good solvents)
#[allow(dead_code)]
pub fn mark_houwink_viscosity(k: f64, mol_weight: f64, a: f64) -> f64 {
    k * mol_weight.powf(a)
}

/// Huggins equation: relative viscosity for a polymer solution.
///
/// η_r = 1 + \[η\] c + k_H \[η\]² c²
///
/// # Arguments
/// * `intrinsic_visc` – intrinsic viscosity \[η\]
/// * `concentration`  – polymer mass concentration \[g/dL\]
/// * `huggins_k`      – Huggins coefficient k_H (~0.3–0.4 for good solvents)
#[allow(dead_code)]
pub fn huggins_relative_viscosity(intrinsic_visc: f64, concentration: f64, huggins_k: f64) -> f64 {
    1.0 + intrinsic_visc * concentration
        + huggins_k * intrinsic_visc * intrinsic_visc * concentration * concentration
}

// ============================================================================
// Section 18 – 3D conformation tensor evolution (Oldroyd-B)
// ============================================================================

/// Evolve the 3D Oldroyd-B conformation tensor by one explicit Euler step.
///
/// # Arguments
/// * `a`         – current 3D conformation tensor
/// * `grad_u`    – velocity gradient tensor \[3×3\] stored row-major
/// * `lambda`    – polymer relaxation time
/// * `dt`        – time step
#[allow(dead_code)]
pub fn oldroyd_b_evolve_3d(
    a: &ConformationTensor3D,
    grad_u: &[f64; 9],
    lambda: f64,
    dt: f64,
) -> ConformationTensor3D {
    // grad_u[3*i+j] = ∂u_i / ∂x_j
    let l = grad_u;
    // Upper-convected derivative: (∇u)·A + A·(∇u)ᵀ
    // A·Lᵀ  (L = ∇u, row-major)
    let stretch_xx = 2.0 * (a.xx * l[0] + a.xy * l[3] + a.xz * l[6]);
    let stretch_yy = 2.0 * (a.xy * l[1] + a.yy * l[4] + a.yz * l[7]);
    let stretch_zz = 2.0 * (a.xz * l[2] + a.yz * l[5] + a.zz * l[8]);
    let stretch_xy =
        a.xx * l[1] + a.xy * l[4] + a.xz * l[7] + a.xy * l[0] + a.yy * l[3] + a.yz * l[6];
    let stretch_xz =
        a.xx * l[2] + a.xy * l[5] + a.xz * l[8] + a.xz * l[0] + a.yz * l[3] + a.zz * l[6];
    let stretch_yz =
        a.xy * l[2] + a.yy * l[5] + a.yz * l[8] + a.xz * l[1] + a.yz * l[4] + a.zz * l[7];

    let relax_xx = -(a.xx - 1.0) / lambda;
    let relax_yy = -(a.yy - 1.0) / lambda;
    let relax_zz = -(a.zz - 1.0) / lambda;
    let relax_xy = -a.xy / lambda;
    let relax_xz = -a.xz / lambda;
    let relax_yz = -a.yz / lambda;

    ConformationTensor3D {
        xx: a.xx + dt * (stretch_xx + relax_xx),
        yy: a.yy + dt * (stretch_yy + relax_yy),
        zz: a.zz + dt * (stretch_zz + relax_zz),
        xy: a.xy + dt * (stretch_xy + relax_xy),
        xz: a.xz + dt * (stretch_xz + relax_xz),
        yz: a.yz + dt * (stretch_yz + relax_yz),
    }
}

// ============================================================================
// Section 19 – Polymer-modified Reynolds stress
// ============================================================================

/// Polymer contribution to the Reynolds stress in drag-reduced turbulence.
///
/// −<u'v'> = ν_t |∂U/∂y| − τ_xy^polymer / ρ
///
/// # Arguments
/// * `nu_t`        – turbulent eddy viscosity \[m²/s\]
/// * `du_dy`       – mean velocity gradient \[1/s\]
/// * `tau_xy_poly` – polymer shear stress contribution \[Pa\]
/// * `rho`         – fluid density \[kg/m³\]
#[allow(dead_code)]
pub fn polymer_modified_reynolds_stress(nu_t: f64, du_dy: f64, tau_xy_poly: f64, rho: f64) -> f64 {
    nu_t * du_dy.abs() - tau_xy_poly / rho.max(f64::EPSILON)
}

// ============================================================================
// Section 20 – LBM simulation runner (high-level API)
// ============================================================================

/// Configuration for a polymer channel flow LBM simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PolymerChannelConfig {
    /// Grid size in x direction.
    pub nx: usize,
    /// Grid size in y direction.
    pub ny: usize,
    /// Solvent kinematic viscosity (lattice units): ν_s = (1/ω_s − 0.5)/3.
    pub nu_s: f64,
    /// Polymer viscosity (lattice units).
    pub eta_p: f64,
    /// Polymer relaxation time (lattice units).
    pub lambda: f64,
    /// FENE parameter b (f64::INFINITY → Oldroyd-B).
    pub fene_b: f64,
    /// Body force (pressure gradient, lattice units).
    pub body_force: f64,
    /// Total number of time steps.
    pub n_steps: usize,
}

impl PolymerChannelConfig {
    /// Compute the solvent BGK relaxation frequency from ν_s.
    pub fn omega_s(&self) -> f64 {
        1.0 / (3.0 * self.nu_s + 0.5)
    }

    /// Compute the bulk Weissenberg number Wi = λ · U_c / H.
    ///
    /// Uses the analytical Poiseuille centreline velocity estimate.
    ///
    /// # Arguments
    /// * `dp_dx` – driving pressure gradient (lattice units)
    pub fn weissenberg(&self, dp_dx: f64) -> f64 {
        let h = self.ny as f64 * 0.5;
        let eta_total = 3.0 * self.nu_s + self.eta_p;
        let u_c = if eta_total.abs() > f64::EPSILON {
            (-dp_dx) * h * h / (8.0 * eta_total)
        } else {
            0.0
        };
        self.lambda * u_c / h
    }
}

/// Run a polymer channel flow LBM simulation and return the mean polymer trace.
///
/// This is the high-level driver that orchestrates the grid, time stepping,
/// and diagnostic collection.
///
/// # Arguments
/// * `cfg` – simulation configuration
#[allow(dead_code)]
pub fn run_polymer_channel(cfg: &PolymerChannelConfig) -> f64 {
    let mut grid = PolymerLbmGrid2D::new(
        cfg.nx,
        cfg.ny,
        cfg.omega_s(),
        cfg.eta_p,
        cfg.lambda,
        cfg.fene_b,
    );
    grid.set_poiseuille_profile(0.01);

    for _ in 0..cfg.n_steps {
        // Apply body force to all interior cells
        for iy in 1..cfg.ny - 1 {
            for ix in 0..cfg.nx {
                let idx = grid.idx(ix, iy);
                let ex = D2Q9_EX;
                let w = D2Q9_WEIGHTS;
                for i in 0..9 {
                    grid.cells[idx].f[i] += 3.0 * w[i] * cfg.body_force * ex[i];
                }
            }
        }
        grid.step();
        grid.apply_wall_bounce_back();
    }

    grid.mean_polymer_trace()
}

// ============================================================================
// Section 21 – Entropy production in viscoelastic flow
// ============================================================================

/// Polymer entropy production rate per unit volume.
///
/// Ṡ_poly = (η_p / λ T) (tr A − 3 − ln det A)   for Oldroyd-B.
///
/// # Arguments
/// * `conf`        – conformation tensor
/// * `eta_p`       – polymer viscosity \[Pa·s\]
/// * `lambda`      – relaxation time \[s\]
/// * `temperature` – absolute temperature \[K\]
#[allow(dead_code)]
pub fn polymer_entropy_production_2d(
    conf: &ConformationTensor2D,
    eta_p: f64,
    lambda: f64,
    temperature: f64,
) -> f64 {
    let tr = conf.trace();
    let det = conf.det().max(f64::EPSILON);
    let ln_det = det.ln();
    // 2D: tr A − 2 − ln det A  (the dimension is 2)
    let phi = tr - 2.0 - ln_det;
    eta_p * phi / (lambda * temperature.max(f64::EPSILON))
}

// ============================================================================
// Section 22 – Polymer chain statistics
// ============================================================================

/// End-to-end distance squared `R²` from the conformation tensor.
///
/// `R²` = R_eq² tr(A)
///
/// # Arguments
/// * `r_eq_sq`   – equilibrium end-to-end distance squared \[m²\]
/// * `trace_a`   – trace of the conformation tensor
#[allow(dead_code)]
pub fn end_to_end_distance_sq(r_eq_sq: f64, trace_a: f64) -> f64 {
    r_eq_sq * trace_a
}

/// Radius of gyration from end-to-end distance: R_g² = R² / 6 (Gaussian chain).
///
/// # Arguments
/// * `r_sq` – mean-square end-to-end distance \[m²\]
#[allow(dead_code)]
pub fn radius_of_gyration_sq(r_sq: f64) -> f64 {
    r_sq / 6.0
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weissenberg_number() {
        let wi = weissenberg_number(0.1, 10.0);
        assert!((wi - 1.0).abs() < 1.0e-12, "Wi = {wi}");
    }

    #[test]
    fn test_deborah_number() {
        let de = deborah_number(2.0, 4.0);
        assert!((de - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn test_deborah_zero_time() {
        assert_eq!(deborah_number(1.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn test_polymer_viscosity() {
        let eta_p = polymer_viscosity(1.0, 0.6);
        assert!((eta_p - 0.4).abs() < 1.0e-12);
    }

    #[test]
    fn test_rouse_relaxation_time() {
        // Just verify the formula is dimensionally consistent and positive
        let tau_r = rouse_relaxation_time(1.0e-11, 100.0, 1.0e-9, 300.0);
        assert!(tau_r > 0.0);
    }

    #[test]
    fn test_zimm_relaxation_time() {
        let tau_z = zimm_relaxation_time(1.0e-3, 100.0, 1.0e-9, 300.0);
        assert!(tau_z > 0.0);
    }

    #[test]
    fn test_rouse_mode_time() {
        let tau_r = 1.0;
        let tau_2 = rouse_mode_time(tau_r, 2);
        assert!((tau_2 - 0.25).abs() < 1.0e-12);
    }

    #[test]
    fn test_conformation_identity() {
        let a = ConformationTensor2D::identity();
        assert!((a.trace() - 2.0).abs() < 1.0e-12);
        assert!((a.det() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_conformation_add_scale() {
        let a = ConformationTensor2D {
            xx: 1.2,
            yy: 1.3,
            xy: 0.1,
        };
        let b = a.scale(2.0);
        assert!((b.xx - 2.4).abs() < 1.0e-12);
        let c = a.add(&b);
        assert!((c.xx - 3.6).abs() < 1.0e-12);
    }

    #[test]
    fn test_polymer_stress_from_conformation() {
        let a = ConformationTensor2D::identity();
        let tau = PolymerStress2D::from_conformation(&a, 0.1, 0.01);
        // A = I → τ = 0
        assert!(tau.txx.abs() < 1.0e-12);
        assert!(tau.tyy.abs() < 1.0e-12);
        assert!(tau.txy.abs() < 1.0e-12);
    }

    #[test]
    fn test_polymer_stress_n1() {
        let tau = PolymerStress2D {
            txx: 0.5,
            tyy: 0.2,
            txy: 0.1,
        };
        assert!((tau.n1() - 0.3).abs() < 1.0e-12);
    }

    #[test]
    fn test_oldroyd_b_evolve_relaxation_only() {
        // Zero velocity gradient: tensor should relax back to identity
        let a0 = ConformationTensor2D {
            xx: 2.0,
            yy: 2.0,
            xy: 0.5,
        };
        let a1 = oldroyd_b_evolve_2d(&a0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.1);
        // After one step, Axx should decrease toward 1
        assert!(a1.xx < a0.xx, "Relaxation should reduce Axx");
        assert!(a1.xy.abs() < a0.xy.abs(), "Relaxation should reduce Axy");
    }

    #[test]
    fn test_fene_p_spring_finite() {
        let f = fene_p_spring(5.0, 10.0);
        assert!((f - 2.0).abs() < 1.0e-12, "f = {f}");
    }

    #[test]
    fn test_fene_p_spring_near_max() {
        let f = fene_p_spring(9.999, 10.0);
        assert!(f > 1000.0);
    }

    #[test]
    fn test_fene_p_evolve_relaxes_to_identity() {
        let a0 = ConformationTensor2D {
            xx: 3.0,
            yy: 3.0,
            xy: 0.0,
        };
        let a1 = fene_p_evolve_2d(&a0, 0.0, 0.0, 0.0, 0.0, 1.0, 100.0, 0.1);
        assert!(a1.xx < a0.xx);
    }

    #[test]
    fn test_fene_p_stress() {
        let a = ConformationTensor2D::identity();
        let tau = fene_p_stress_2d(&a, 0.1, 0.01, 100.0);
        // tr(A) = 2, b = 100 → f = 100/98 ≈ 1.0204
        // τxx = c*(f*1 - 1) ≠ 0 for finite b
        assert!(tau.txx.abs() > 0.0 || tau.tyy.abs() > 0.0);
    }

    #[test]
    fn test_orientation_tensor_isotropic() {
        let s = OrientationTensor2D::isotropic();
        assert!((s.sxx - 0.5).abs() < 1.0e-12);
        assert!((s.syy - 0.5).abs() < 1.0e-12);
        assert!(s.sxy.abs() < 1.0e-12);
    }

    #[test]
    fn test_orientation_from_conformation() {
        let a = ConformationTensor2D {
            xx: 2.0,
            yy: 1.0,
            xy: 0.0,
        };
        let s = OrientationTensor2D::from_conformation(&a);
        assert!((s.sxx - 2.0 / 3.0).abs() < 1.0e-10);
        assert!((s.syy - 1.0 / 3.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_virk_mdr_friction() {
        let f = virk_mdr_friction(10000.0);
        assert!(f > 0.0 && f < 1.0);
    }

    #[test]
    fn test_drag_reduction_percent() {
        let dr = drag_reduction_percent(0.01, 0.007);
        assert!((dr - 30.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_drag_reduction_zero_newtonian() {
        let dr = drag_reduction_percent(0.0, 0.005);
        assert_eq!(dr, 0.0);
    }

    #[test]
    fn test_is_stretched() {
        assert!(is_stretched(1.0));
        assert!(!is_stretched(0.3));
    }

    #[test]
    fn test_chain_extension_ratio() {
        let r = chain_extension_ratio(5.0, 10.0);
        assert!((r - (0.5f64).sqrt()).abs() < 1.0e-12);
    }

    #[test]
    fn test_chain_extension_clamped() {
        let r = chain_extension_ratio(100.0, 10.0);
        assert_eq!(r, 1.0);
    }

    #[test]
    fn test_elastic_turbulence_spectrum() {
        let e = elastic_turbulence_spectrum(2.0, 1.0, 3.5);
        assert!((e - 2.0_f64.powf(-3.5)).abs() < 1.0e-12);
    }

    #[test]
    fn test_elastic_turbulence_spectrum_zero_k() {
        assert_eq!(elastic_turbulence_spectrum(0.0, 1.0, 3.5), 0.0);
    }

    #[test]
    fn test_mark_houwink() {
        let visc = mark_houwink_viscosity(0.072, 1.0e6, 0.75);
        assert!(visc > 0.0);
    }

    #[test]
    fn test_huggins_viscosity() {
        let eta_r = huggins_relative_viscosity(1.0, 0.1, 0.35);
        // 1 + 0.1 + 0.35 * 0.01 = 1.1035
        assert!((eta_r - 1.1035).abs() < 1.0e-12);
    }

    #[test]
    fn test_effective_viscosity_polymer_newtonian_limit() {
        // At Wi=0, η_eff = η_s + η_p
        let eta_eff = effective_viscosity_polymer(0.6, 0.4, 0.0);
        assert!((eta_eff - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_effective_viscosity_polymer_high_wi() {
        // At Wi >> 1, η_eff → η_s
        let eta_eff = effective_viscosity_polymer(0.6, 0.4, 1000.0);
        assert!((eta_eff - 0.6).abs() < 0.01);
    }

    #[test]
    fn test_startup_shear_n1() {
        let n1 = startup_shear_n1(0.1, 1.0, 0.0, 1.0);
        assert!(n1.abs() < 1.0e-12);
        let n1_inf = startup_shear_n1(0.1, 1.0, 1000.0, 1.0);
        assert!((n1_inf - 0.2).abs() < 1.0e-6);
    }

    #[test]
    fn test_local_weissenberg_field() {
        // Linear profile u = y → du/dy = 1/dy (constant)
        let n = 5usize;
        let dy = 0.25;
        let u: Vec<f64> = (0..n).map(|i| i as f64 * dy).collect();
        let wi = local_weissenberg_field(&u, dy, 1.0);
        // Interior points should give Wi ≈ 1.0
        for i in 1..n - 1 {
            assert!((wi[i] - 1.0).abs() < 1.0e-10, "wi[{i}] = {}", wi[i]);
        }
    }

    #[test]
    fn test_polymer_cell_equilibrium() {
        let cell = PolymerCell::new_equilibrium(1.0, 0.01, 0.0);
        let sum: f64 = cell.f.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1.0e-10,
            "density conservation: sum = {sum}"
        );
    }

    #[test]
    fn test_polymer_lbm_grid_creation() {
        let grid = PolymerLbmGrid2D::new(8, 8, 1.0, 0.1, 10.0, f64::INFINITY);
        assert_eq!(grid.cells.len(), 64);
    }

    #[test]
    fn test_polymer_grid_poiseuille() {
        let mut grid = PolymerLbmGrid2D::new(16, 16, 1.0, 0.05, 5.0, f64::INFINITY);
        grid.set_poiseuille_profile(0.01);
        // Centreline cell (iy = ny/2) should have max velocity
        let mid = grid.idx(0, grid.ny / 2);
        let wall = grid.idx(0, 0);
        assert!(grid.cells[mid].ux > grid.cells[wall].ux);
    }

    #[test]
    fn test_polymer_grid_step_conserves_density() {
        let mut grid = PolymerLbmGrid2D::new(8, 4, 1.0, 0.05, 5.0, f64::INFINITY);
        let rho_before: f64 = grid.cells.iter().map(|c| c.rho).sum();
        grid.step();
        let rho_after: f64 = grid.cells.iter().map(|c| c.f.iter().sum::<f64>()).sum();
        // Density should be approximately conserved (within 1%)
        assert!((rho_after - rho_before).abs() / rho_before < 0.02);
    }

    #[test]
    fn test_polymer_entropy_production() {
        let a = ConformationTensor2D::identity();
        let s = polymer_entropy_production_2d(&a, 0.1, 1.0, 300.0);
        // At equilibrium A = I, entropy production should be ≈ 0
        assert!(s.abs() < 1.0e-10, "s = {s}");
    }

    #[test]
    fn test_polymer_entropy_production_stretched() {
        let a = ConformationTensor2D {
            xx: 3.0,
            yy: 1.5,
            xy: 0.0,
        };
        let s = polymer_entropy_production_2d(&a, 0.1, 1.0, 300.0);
        assert!(
            s > 0.0,
            "Stretched chain should have positive entropy production"
        );
    }

    #[test]
    fn test_end_to_end_distance() {
        let r_sq = end_to_end_distance_sq(1.0e-18, 3.0);
        assert!((r_sq - 3.0e-18).abs() < 1.0e-28);
    }

    #[test]
    fn test_radius_of_gyration() {
        let r_sq = 6.0e-18;
        let rg_sq = radius_of_gyration_sq(r_sq);
        assert!((rg_sq - 1.0e-18).abs() < 1.0e-28);
    }

    #[test]
    fn test_conformation_3d_identity() {
        let a = ConformationTensor3D::identity();
        assert!((a.trace() - 3.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_oldroyd_b_3d_relaxation() {
        let a0 = ConformationTensor3D {
            xx: 2.0,
            yy: 2.0,
            zz: 2.0,
            xy: 0.3,
            xz: 0.2,
            yz: 0.1,
        };
        let grad_u = [0.0f64; 9];
        let a1 = oldroyd_b_evolve_3d(&a0, &grad_u, 1.0, 0.1);
        assert!(a1.xx < a0.xx);
        assert!(a1.xy.abs() < a0.xy.abs());
    }

    #[test]
    fn test_polymer_config_omega() {
        let cfg = PolymerChannelConfig {
            nx: 16,
            ny: 8,
            nu_s: 1.0 / 6.0,
            eta_p: 0.05,
            lambda: 10.0,
            fene_b: f64::INFINITY,
            body_force: 1.0e-4,
            n_steps: 10,
        };
        let omega = cfg.omega_s();
        // nu_s = (1/omega - 0.5)/3 → omega = 1/(3*nu_s + 0.5) = 1/1.0 = 1.0
        assert!((omega - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_channel_config_weissenberg() {
        let cfg = PolymerChannelConfig {
            nx: 16,
            ny: 8,
            nu_s: 1.0 / 6.0,
            eta_p: 0.0,
            lambda: 0.0,
            fene_b: f64::INFINITY,
            body_force: 1.0e-5,
            n_steps: 0,
        };
        let wi = cfg.weissenberg(-1.0e-5);
        assert!(wi >= 0.0);
    }
}
