// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SIMD-accelerated bulk material evaluation using Structure-of-Arrays (SoA) layout.
//!
//! This module provides batch (N-point) evaluation kernels for the most
//! compute-intensive constitutive relations in OxiPhysics.  Data is arranged in
//! SoA format so that LLVM's auto-vectorizer can emit wide SIMD instructions
//! (SSE2 / AVX2 / NEON) on the inner loops.
//!
//! ## Constitutive kernels
//!
//! | Function | Model | Description |
//! |---|---|---|
//! | [`elastic_stress_batch`] | Isotropic linear elastic | σ = λ tr(ε) I + 2μ ε |
//! | [`von_mises_yield_batch`] | Von Mises plasticity | f = √(3/2 s:s) − σ_y |
//! | [`neo_hookean_stress_batch`] | Neo-Hookean hyperelastic | First Piola–Kirchhoff via F |
//! | [`viscoplastic_rate_batch`] | Perzyna viscoplasticity | γ̇ = ⟨f/σ_y⟩ⁿ / η |
//! | [`miner_damage_batch`] | Miner's rule fatigue | D += n/N_f per cycle block |
//! | [`thermal_expansion_stress_batch`] | Thermo-elastic | σ_th = −3K α ΔT I |
//!
//! ## SoA storage
//!
//! [`SoaMaterialPoints`] stores the full state of N integration points with each
//! field in its own contiguous `Vec<f64>` so the batch kernels read and write
//! aligned slices that the compiler can vectorize.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_materials::simd_paths::{SoaMaterialPoints, elastic_stress_batch};
//!
//! let n = 8;
//! let mut pts = SoaMaterialPoints::new(n);
//! // Set strains: ε_xx = 0.001 for all points
//! for i in 0..n { pts.strain_xx[i] = 0.001; }
//!
//! let e  = 210e9_f64;
//! let nu = 0.3;
//! elastic_stress_batch(&mut pts, e, nu);
//!
//! // σ_xx should equal E/(1-ν²) * ε_xx for plane stress
//! assert!(pts.stress_xx[0] > 0.0);
//! ```

// ── SoaMaterialPoints ─────────────────────────────────────────────────────────

/// Structure-of-Arrays storage for `n` material integration points.
///
/// All arrays have the same length `n` so batch kernels can iterate over them
/// with simple index loops that LLVM can auto-vectorize.
#[derive(Debug, Clone, Default)]
pub struct SoaMaterialPoints {
    /// Number of material points.
    pub n: usize,

    // ── Strain tensor components (Voigt notation, symmetric) ──────────────
    /// Normal strain εxx.
    pub strain_xx: Vec<f64>,
    /// Normal strain εyy.
    pub strain_yy: Vec<f64>,
    /// Normal strain εzz.
    pub strain_zz: Vec<f64>,
    /// Shear strain εxy (engineering = 2 * tensorial).
    pub strain_xy: Vec<f64>,
    /// Shear strain εyz.
    pub strain_yz: Vec<f64>,
    /// Shear strain εxz.
    pub strain_xz: Vec<f64>,

    // ── Stress tensor components ──────────────────────────────────────────
    /// Normal stress σxx (output from constitutive kernels).
    pub stress_xx: Vec<f64>,
    /// Normal stress σyy.
    pub stress_yy: Vec<f64>,
    /// Normal stress σzz.
    pub stress_zz: Vec<f64>,
    /// Shear stress σxy.
    pub stress_xy: Vec<f64>,
    /// Shear stress σyz.
    pub stress_yz: Vec<f64>,
    /// Shear stress σxz.
    pub stress_xz: Vec<f64>,

    // ── Scalar state fields ───────────────────────────────────────────────
    /// Accumulated equivalent plastic strain p̄.
    pub plastic_strain: Vec<f64>,
    /// Current yield stress σ_y (may harden).
    pub yield_stress: Vec<f64>,
    /// Scalar damage variable D ∈ [0, 1].
    pub damage: Vec<f64>,
    /// Temperature (K).
    pub temperature: Vec<f64>,
    /// Accumulated fatigue damage (Miner's rule) ∈ [0, 1].
    pub fatigue_damage: Vec<f64>,

    // ── Deformation gradient (row-major 3×3) ─────────────────────────────
    /// F₁₁ components (one per point).
    pub f11: Vec<f64>,
    /// F₁₂.
    pub f12: Vec<f64>,
    /// F₁₃.
    pub f13: Vec<f64>,
    /// F₂₁.
    pub f21: Vec<f64>,
    /// F₂₂.
    pub f22: Vec<f64>,
    /// F₂₃.
    pub f23: Vec<f64>,
    /// F₃₁.
    pub f31: Vec<f64>,
    /// F₃₂.
    pub f32: Vec<f64>,
    /// F₃₃.
    pub f33: Vec<f64>,
}

impl SoaMaterialPoints {
    /// Allocate `n` material points, all zero-initialised.
    pub fn new(n: usize) -> Self {
        // Identity deformation gradients
        let mut pts = Self {
            n,
            strain_xx: vec![0.0; n],
            strain_yy: vec![0.0; n],
            strain_zz: vec![0.0; n],
            strain_xy: vec![0.0; n],
            strain_yz: vec![0.0; n],
            strain_xz: vec![0.0; n],
            stress_xx: vec![0.0; n],
            stress_yy: vec![0.0; n],
            stress_zz: vec![0.0; n],
            stress_xy: vec![0.0; n],
            stress_yz: vec![0.0; n],
            stress_xz: vec![0.0; n],
            plastic_strain: vec![0.0; n],
            yield_stress: vec![250e6; n], // default 250 MPa
            damage: vec![0.0; n],
            temperature: vec![293.15; n], // room temperature K
            fatigue_damage: vec![0.0; n],
            f11: vec![1.0; n],
            f12: vec![0.0; n],
            f13: vec![0.0; n],
            f21: vec![0.0; n],
            f22: vec![1.0; n],
            f23: vec![0.0; n],
            f31: vec![0.0; n],
            f32: vec![0.0; n],
            f33: vec![1.0; n],
        };
        // yield_stress default overrides
        pts.yield_stress.fill(250e6);
        pts
    }

    /// Reset all stresses to zero.
    pub fn zero_stress(&mut self) {
        self.stress_xx.fill(0.0);
        self.stress_yy.fill(0.0);
        self.stress_zz.fill(0.0);
        self.stress_xy.fill(0.0);
        self.stress_yz.fill(0.0);
        self.stress_xz.fill(0.0);
    }

    /// Compute Von Mises equivalent stress for each point.
    ///
    /// σ_vm = √( ½ [(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²] + 3(σxy² + σyz² + σxz²) )
    pub fn von_mises_stress(&self) -> Vec<f64> {
        let n = self.n;
        let mut vm = Vec::with_capacity(n);
        for i in 0..n {
            let dxy = self.stress_xx[i] - self.stress_yy[i];
            let dyz = self.stress_yy[i] - self.stress_zz[i];
            let dzx = self.stress_zz[i] - self.stress_xx[i];
            let shear = self.stress_xy[i] * self.stress_xy[i]
                + self.stress_yz[i] * self.stress_yz[i]
                + self.stress_xz[i] * self.stress_xz[i];
            vm.push((0.5 * (dxy * dxy + dyz * dyz + dzx * dzx) + 3.0 * shear).sqrt());
        }
        vm
    }

    /// Mean (hydrostatic) stress for each point.
    pub fn hydrostatic_stress(&self) -> Vec<f64> {
        (0..self.n)
            .map(|i| (self.stress_xx[i] + self.stress_yy[i] + self.stress_zz[i]) / 3.0)
            .collect()
    }
}

// ── Lamé parameters ───────────────────────────────────────────────────────────

/// Compute Lamé parameters (λ, μ) from Young's modulus E and Poisson ratio ν.
#[inline]
pub fn lame(e: f64, nu: f64) -> (f64, f64) {
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    (lambda, mu)
}

// ── elastic_stress_batch ──────────────────────────────────────────────────────

/// Compute isotropic linear elastic stress for all `n` material points.
///
/// σ = λ tr(ε) I + 2μ ε
///
/// Updates `pts.stress_{xx,yy,zz,xy,yz,xz}` in place.  The loop is written in
/// a form that LLVM can auto-vectorize with 4-wide f64 SIMD lanes.
pub fn elastic_stress_batch(pts: &mut SoaMaterialPoints, e: f64, nu: f64) {
    let (lambda, mu) = lame(e, nu);
    let n = pts.n;
    let two_mu = 2.0 * mu;

    // Batch loop: auto-vectorized by LLVM into AVX2 (4×f64) or SSE2 (2×f64)
    for i in 0..n {
        let tr = pts.strain_xx[i] + pts.strain_yy[i] + pts.strain_zz[i];
        let ltr = lambda * tr;
        pts.stress_xx[i] = ltr + two_mu * pts.strain_xx[i];
        pts.stress_yy[i] = ltr + two_mu * pts.strain_yy[i];
        pts.stress_zz[i] = ltr + two_mu * pts.strain_zz[i];
        pts.stress_xy[i] = mu * pts.strain_xy[i]; // engineering shear
        pts.stress_yz[i] = mu * pts.strain_yz[i];
        pts.stress_xz[i] = mu * pts.strain_xz[i];
    }
}

/// Same as [`elastic_stress_batch`] but processes 4 points at once with
/// manually unrolled 4-wide scalar arithmetic (auto-vectorizer helper).
///
/// Use for very large point clouds where the batch size is always a multiple of 4.
#[inline(always)]
pub fn elastic_stress_batch_x4(pts: &mut SoaMaterialPoints, e: f64, nu: f64) {
    let (lambda, mu) = lame(e, nu);
    let two_mu = 2.0 * mu;
    let n = pts.n;
    let chunks = n / 4;

    for c in 0..chunks {
        let base = c * 4;
        // Process 4 points per iteration — LLVM fuses into a single AVX2 kernel
        let tr0 = pts.strain_xx[base] + pts.strain_yy[base] + pts.strain_zz[base];
        let tr1 = pts.strain_xx[base + 1] + pts.strain_yy[base + 1] + pts.strain_zz[base + 1];
        let tr2 = pts.strain_xx[base + 2] + pts.strain_yy[base + 2] + pts.strain_zz[base + 2];
        let tr3 = pts.strain_xx[base + 3] + pts.strain_yy[base + 3] + pts.strain_zz[base + 3];

        let ltr0 = lambda * tr0;
        let ltr1 = lambda * tr1;
        let ltr2 = lambda * tr2;
        let ltr3 = lambda * tr3;

        pts.stress_xx[base] = ltr0 + two_mu * pts.strain_xx[base];
        pts.stress_xx[base + 1] = ltr1 + two_mu * pts.strain_xx[base + 1];
        pts.stress_xx[base + 2] = ltr2 + two_mu * pts.strain_xx[base + 2];
        pts.stress_xx[base + 3] = ltr3 + two_mu * pts.strain_xx[base + 3];

        pts.stress_yy[base] = ltr0 + two_mu * pts.strain_yy[base];
        pts.stress_yy[base + 1] = ltr1 + two_mu * pts.strain_yy[base + 1];
        pts.stress_yy[base + 2] = ltr2 + two_mu * pts.strain_yy[base + 2];
        pts.stress_yy[base + 3] = ltr3 + two_mu * pts.strain_yy[base + 3];

        pts.stress_zz[base] = ltr0 + two_mu * pts.strain_zz[base];
        pts.stress_zz[base + 1] = ltr1 + two_mu * pts.strain_zz[base + 1];
        pts.stress_zz[base + 2] = ltr2 + two_mu * pts.strain_zz[base + 2];
        pts.stress_zz[base + 3] = ltr3 + two_mu * pts.strain_zz[base + 3];

        pts.stress_xy[base] = mu * pts.strain_xy[base];
        pts.stress_xy[base + 1] = mu * pts.strain_xy[base + 1];
        pts.stress_xy[base + 2] = mu * pts.strain_xy[base + 2];
        pts.stress_xy[base + 3] = mu * pts.strain_xy[base + 3];

        pts.stress_yz[base] = mu * pts.strain_yz[base];
        pts.stress_yz[base + 1] = mu * pts.strain_yz[base + 1];
        pts.stress_yz[base + 2] = mu * pts.strain_yz[base + 2];
        pts.stress_yz[base + 3] = mu * pts.strain_yz[base + 3];

        pts.stress_xz[base] = mu * pts.strain_xz[base];
        pts.stress_xz[base + 1] = mu * pts.strain_xz[base + 1];
        pts.stress_xz[base + 2] = mu * pts.strain_xz[base + 2];
        pts.stress_xz[base + 3] = mu * pts.strain_xz[base + 3];
    }

    // Scalar tail
    for i in (chunks * 4)..n {
        let tr = pts.strain_xx[i] + pts.strain_yy[i] + pts.strain_zz[i];
        let ltr = lambda * tr;
        pts.stress_xx[i] = ltr + two_mu * pts.strain_xx[i];
        pts.stress_yy[i] = ltr + two_mu * pts.strain_yy[i];
        pts.stress_zz[i] = ltr + two_mu * pts.strain_zz[i];
        pts.stress_xy[i] = mu * pts.strain_xy[i];
        pts.stress_yz[i] = mu * pts.strain_yz[i];
        pts.stress_xz[i] = mu * pts.strain_xz[i];
    }
}

// ── von_mises_yield_batch ─────────────────────────────────────────────────────

/// Evaluate the Von Mises yield function for each material point.
///
/// `f_i = σ_vm_i − σ_y_i`
///
/// * Positive values indicate plastic flow.
/// * `pts.stress_*` must be populated before calling (e.g. via [`elastic_stress_batch`]).
///
/// Returns a `Vec<f64>` of yield function values.
pub fn von_mises_yield_batch(pts: &SoaMaterialPoints) -> Vec<f64> {
    let vm = pts.von_mises_stress();
    (0..pts.n).map(|i| vm[i] - pts.yield_stress[i]).collect()
}

// ── neo_hookean_stress_batch ──────────────────────────────────────────────────

/// Compute the Cauchy stress for a Neo-Hookean hyperelastic material at each point.
///
/// The Neo-Hookean model:
/// ```text
/// Ψ = μ/2 (tr(b) − 3) − μ ln(J) + λ/2 (ln J)²
/// σ = (1/J) [ μ (b − I) + λ ln(J) I ]
/// ```
/// where `b = F Fᵀ` (left Cauchy-Green tensor) and `J = det(F)`.
///
/// Updates `pts.stress_*` from the stored deformation gradients `pts.f*`.
pub fn neo_hookean_stress_batch(pts: &mut SoaMaterialPoints, e: f64, nu: f64) {
    let (lambda, mu) = lame(e, nu);
    let n = pts.n;

    for i in 0..n {
        // Deformation gradient F (row-major)
        let f = [
            [pts.f11[i], pts.f12[i], pts.f13[i]],
            [pts.f21[i], pts.f22[i], pts.f23[i]],
            [pts.f31[i], pts.f32[i], pts.f33[i]],
        ];

        // J = det(F)
        let j = det3(f);
        if j.abs() < 1e-15 {
            continue; // degenerate element, skip
        }
        let ln_j = j.ln();
        let j_inv = 1.0 / j;

        // b = F Fᵀ (left Cauchy-Green)
        let b = mat_mat_t(f, f);

        // σ = (1/J) [ μ (b − I) + λ ln(J) I ]
        pts.stress_xx[i] = j_inv * (mu * (b[0][0] - 1.0) + lambda * ln_j);
        pts.stress_yy[i] = j_inv * (mu * (b[1][1] - 1.0) + lambda * ln_j);
        pts.stress_zz[i] = j_inv * (mu * (b[2][2] - 1.0) + lambda * ln_j);
        pts.stress_xy[i] = j_inv * mu * b[0][1];
        pts.stress_yz[i] = j_inv * mu * b[1][2];
        pts.stress_xz[i] = j_inv * mu * b[0][2];
    }
}

// ── viscoplastic_rate_batch ───────────────────────────────────────────────────

/// Compute the viscoplastic strain rate for each point (Perzyna model).
///
/// `γ̇_i = ⟨f_i / σ_y_i⟩ⁿ / η`
///
/// where `⟨x⟩ = max(0, x)` (Macaulay bracket).
///
/// * `yield_values` — yield function values (from [`von_mises_yield_batch`])
/// * `n_exp`        — rate sensitivity exponent (typical 1–10)
/// * `eta`          — viscosity parameter (typical 1e-3 – 1e-1)
///
/// Returns `Vec<f64>` of plastic strain rates.
pub fn viscoplastic_rate_batch(
    pts: &SoaMaterialPoints,
    yield_values: &[f64],
    n_exp: f64,
    eta: f64,
) -> Vec<f64> {
    assert_eq!(yield_values.len(), pts.n);
    let eta_inv = 1.0 / eta.max(1e-300);
    (0..pts.n)
        .map(|i| {
            let f = yield_values[i].max(0.0);
            let overstress = f / pts.yield_stress[i].max(1e-15);
            overstress.powf(n_exp) * eta_inv
        })
        .collect()
}

// ── miner_damage_batch ────────────────────────────────────────────────────────

/// Accumulate Miner's rule fatigue damage for a block of `n_applied` loading cycles.
///
/// `D_i += n_applied / N_f_i`
///
/// where `N_f_i` is the number of cycles to failure for the current stress amplitude
/// computed via the Basquin relation:
///
/// `N_f = (σ_f' / σ_a)^(1/b)`
///
/// * `sigma_f_prime` — fatigue strength coefficient
/// * `b`             — Basquin exponent (typical −0.05 to −0.12)
/// * `n_applied`     — number of cycles applied this block
///
/// Clips `D_i` to `[0, 1]`.  Points where `D_i == 1.0` have failed.
pub fn miner_damage_batch(pts: &mut SoaMaterialPoints, sigma_f_prime: f64, b: f64, n_applied: f64) {
    let vm = pts.von_mises_stress();
    for (i, &sigma_v) in vm.iter().enumerate() {
        let sigma_a = sigma_v.max(1e-15);
        // N_f from Basquin: σ_a = σ_f' * (2 N_f)^b  →  N_f = ½ (σ_a / σ_f')^(1/b)
        // b is negative (typical -0.05 to -0.12), so 1/b is also negative, which means
        // higher stress amplitude → smaller N_f → more damage per cycle (physically correct).
        let n_f = 0.5 * (sigma_a / sigma_f_prime.max(1e-15)).powf(1.0 / b);
        if n_f > 1e-15 {
            pts.fatigue_damage[i] = (pts.fatigue_damage[i] + n_applied / n_f).min(1.0);
        }
    }
}

// ── thermal_expansion_stress_batch ────────────────────────────────────────────

/// Add the thermal expansion stress increment to existing stress.
///
/// `Δσ_th = −3 K α ΔT I`
///
/// where `K = E / (3(1 − 2ν))` is the bulk modulus, `α` is the coefficient of
/// thermal expansion, and `ΔT = T − T_ref`.
///
/// * `alpha` — coefficient of thermal expansion (1/K), e.g. 12e-6 for steel
/// * `t_ref` — stress-free reference temperature (K)
pub fn thermal_expansion_stress_batch(
    pts: &mut SoaMaterialPoints,
    e: f64,
    nu: f64,
    alpha: f64,
    t_ref: f64,
) {
    let k_bulk = e / (3.0 * (1.0 - 2.0 * nu));
    let n = pts.n;
    for i in 0..n {
        let dt = pts.temperature[i] - t_ref;
        let th_stress = -3.0 * k_bulk * alpha * dt;
        pts.stress_xx[i] += th_stress;
        pts.stress_yy[i] += th_stress;
        pts.stress_zz[i] += th_stress;
    }
}

// ── return_mapping_batch ──────────────────────────────────────────────────────

/// Radial return mapping for isotropic J2 plasticity (batch, linear isotropic hardening).
///
/// Updates both the stress and accumulated plastic strain in place.
///
/// Algorithm (per point):
/// 1. Elastic trial stress (already in `pts.stress_*`)
/// 2. Compute trial Von Mises: `σ_vm_trial`
/// 3. If trial < yield: elastic step, do nothing
/// 4. Else: plastic correction — scale deviatoric back to yield surface
///
/// * `e`, `nu`  — elastic constants
/// * `h`        — isotropic hardening modulus (H = dσ_y / dp̄)
pub fn return_mapping_batch(pts: &mut SoaMaterialPoints, e: f64, nu: f64, h: f64) {
    let mu = e / (2.0 * (1.0 + nu));
    let n = pts.n;

    for i in 0..n {
        let vm_trial = {
            let dxy = pts.stress_xx[i] - pts.stress_yy[i];
            let dyz = pts.stress_yy[i] - pts.stress_zz[i];
            let dzx = pts.stress_zz[i] - pts.stress_xx[i];
            (0.5 * (dxy * dxy + dyz * dyz + dzx * dzx)
                + 3.0
                    * (pts.stress_xy[i] * pts.stress_xy[i]
                        + pts.stress_yz[i] * pts.stress_yz[i]
                        + pts.stress_xz[i] * pts.stress_xz[i]))
                .sqrt()
        };

        let sigma_y = pts.yield_stress[i];
        if vm_trial <= sigma_y + 1e-12 {
            continue; // elastic
        }

        // Plastic multiplier Δγ (linear hardening)
        let dg = (vm_trial - sigma_y) / (3.0 * mu + h);

        // Hydrostatic stress (unchanged in J2)
        let p = (pts.stress_xx[i] + pts.stress_yy[i] + pts.stress_zz[i]) / 3.0;

        // Deviatoric stress components
        let sx = pts.stress_xx[i] - p;
        let sy = pts.stress_yy[i] - p;
        let sz = pts.stress_zz[i] - p;
        let txy = pts.stress_xy[i];
        let tyz = pts.stress_yz[i];
        let txz = pts.stress_xz[i];

        // Scale deviatoric back: s_corrected = s_trial * (1 - 3μΔγ/σ_vm_trial)
        let scale = (1.0 - 3.0 * mu * dg / vm_trial.max(1e-15)).max(0.0);
        pts.stress_xx[i] = p + sx * scale;
        pts.stress_yy[i] = p + sy * scale;
        pts.stress_zz[i] = p + sz * scale;
        pts.stress_xy[i] = txy * scale;
        pts.stress_yz[i] = tyz * scale;
        pts.stress_xz[i] = txz * scale;

        // Update internal variables
        pts.plastic_strain[i] += dg;
        pts.yield_stress[i] = sigma_y + h * dg; // linear isotropic hardening
    }
}

// ── Matrix helpers ─────────────────────────────────────────────────────────────

/// Determinant of a 3×3 matrix (row-major).
#[inline(always)]
fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Matrix product `A * Bᵀ` for 3×3 matrices.
#[inline(always)]
fn mat_mat_t(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[j][k]; // A * B^T
            }
        }
    }
    c
}

// ── MaterialBatchStats ────────────────────────────────────────────────────────

/// Summary statistics for a batch of material points.
#[derive(Debug, Clone)]
pub struct MaterialBatchStats {
    /// Number of points.
    pub n: usize,
    /// Number of plastically active points (yield function > 0).
    pub plastic_count: usize,
    /// Number of failed points (damage >= 1.0 or fatigue >= 1.0).
    pub failed_count: usize,
    /// Maximum Von Mises stress.
    pub max_von_mises: f64,
    /// Mean Von Mises stress.
    pub mean_von_mises: f64,
    /// Maximum accumulated plastic strain.
    pub max_plastic_strain: f64,
    /// Maximum fatigue damage.
    pub max_fatigue_damage: f64,
}

impl MaterialBatchStats {
    /// Compute statistics from the current state of `pts`.
    pub fn compute(pts: &SoaMaterialPoints, yield_values: &[f64]) -> Self {
        let vm = pts.von_mises_stress();
        let n = pts.n;
        let max_vm = vm.iter().cloned().fold(0.0_f64, f64::max);
        let mean_vm = vm.iter().sum::<f64>() / n.max(1) as f64;
        let plastic_count = yield_values.iter().filter(|&&f| f > 0.0).count();
        let failed_count = (0..n)
            .filter(|&i| pts.damage[i] >= 1.0 || pts.fatigue_damage[i] >= 1.0)
            .count();
        let max_ps = pts.plastic_strain.iter().cloned().fold(0.0_f64, f64::max);
        let max_fd = pts.fatigue_damage.iter().cloned().fold(0.0_f64, f64::max);
        Self {
            n,
            plastic_count,
            failed_count,
            max_von_mises: max_vm,
            mean_von_mises: mean_vm,
            max_plastic_strain: max_ps,
            max_fatigue_damage: max_fd,
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const E_STEEL: f64 = 210e9;
    const NU_STEEL: f64 = 0.3;

    #[test]
    fn elastic_uniaxial_stress() {
        let n = 4;
        let mut pts = SoaMaterialPoints::new(n);
        let eps = 0.001;
        for i in 0..n {
            pts.strain_xx[i] = eps;
        }
        elastic_stress_batch(&mut pts, E_STEEL, NU_STEEL);
        // σxx = E ε / (1 - ν²) is NOT the 3D formula; in 3D: σxx = λε + 2με = (λ + 2μ)ε
        let (lambda, mu) = lame(E_STEEL, NU_STEEL);
        let expected_xx = (lambda + 2.0 * mu) * eps;
        for i in 0..n {
            assert!(
                (pts.stress_xx[i] - expected_xx).abs() < 1.0, // within 1 Pa
                "σxx[{i}] = {:.3e}, expected {:.3e}",
                pts.stress_xx[i],
                expected_xx
            );
        }
    }

    #[test]
    fn x4_matches_scalar() {
        let n = 8;
        let mut pts1 = SoaMaterialPoints::new(n);
        let mut pts2 = SoaMaterialPoints::new(n);
        for i in 0..n {
            pts1.strain_xx[i] = 0.001 * (i + 1) as f64;
            pts2.strain_xx[i] = pts1.strain_xx[i];
        }
        elastic_stress_batch(&mut pts1, E_STEEL, NU_STEEL);
        elastic_stress_batch_x4(&mut pts2, E_STEEL, NU_STEEL);
        for i in 0..n {
            assert!(
                (pts1.stress_xx[i] - pts2.stress_xx[i]).abs() < 1e-6,
                "x4 mismatch at {i}: {} vs {}",
                pts1.stress_xx[i],
                pts2.stress_xx[i]
            );
        }
    }

    #[test]
    fn neo_hookean_identity_gives_zero_stress() {
        let n = 2;
        let mut pts = SoaMaterialPoints::new(n); // F = I by default
        neo_hookean_stress_batch(&mut pts, E_STEEL, NU_STEEL);
        for i in 0..n {
            assert!(
                pts.stress_xx[i].abs() < 1e-3,
                "σxx should be ~0 at identity F"
            );
        }
    }

    #[test]
    fn return_mapping_elastic_stays_elastic() {
        let mut pts = SoaMaterialPoints::new(1);
        pts.strain_xx[0] = 1e-6; // tiny strain, well below yield
        elastic_stress_batch(&mut pts, E_STEEL, NU_STEEL);
        let before = pts.plastic_strain[0];
        return_mapping_batch(&mut pts, E_STEEL, NU_STEEL, 2e9);
        assert_eq!(pts.plastic_strain[0], before, "should remain elastic");
    }

    #[test]
    fn return_mapping_plastic_updates_state() {
        let mut pts = SoaMaterialPoints::new(1);
        // Force high stress: set stress directly above yield
        pts.stress_xx[0] = 500e6; // 500 MPa >> 250 MPa yield
        pts.yield_stress[0] = 250e6;
        return_mapping_batch(&mut pts, E_STEEL, NU_STEEL, 0.0); // H=0 (perfect plasticity)
        // After return mapping, Von Mises should be ≤ yield + small tolerance
        let vm = pts.von_mises_stress()[0];
        assert!(
            vm <= pts.yield_stress[0] + 1e6,
            "stress not returned to yield surface: σ_vm={vm:.3e}"
        );
        assert!(
            pts.plastic_strain[0] > 0.0,
            "plastic strain should accumulate"
        );
    }

    #[test]
    fn miner_damage_accumulates() {
        let mut pts = SoaMaterialPoints::new(2);
        pts.stress_xx[0] = 300e6; // high amplitude
        pts.stress_xx[1] = 100e6; // low amplitude
        miner_damage_batch(&mut pts, 1000e6, -0.1, 100.0);
        assert!(
            pts.fatigue_damage[0] > pts.fatigue_damage[1],
            "high-stress point should accumulate more damage"
        );
    }

    #[test]
    fn thermal_expansion_compressive_on_heating() {
        let mut pts = SoaMaterialPoints::new(1);
        // Start with zero stress, heat the point
        pts.temperature[0] = 393.15; // 120 °C
        thermal_expansion_stress_batch(&mut pts, E_STEEL, NU_STEEL, 12e-6, 293.15);
        // Constrained heating → compressive hydrostatic stress
        assert!(
            pts.stress_xx[0] < 0.0,
            "constrained thermal expansion should be compressive"
        );
    }

    #[test]
    fn von_mises_zero_for_hydrostatic() {
        let mut pts = SoaMaterialPoints::new(1);
        // Pure hydrostatic: σxx = σyy = σzz = 100 MPa
        pts.stress_xx[0] = 100e6;
        pts.stress_yy[0] = 100e6;
        pts.stress_zz[0] = 100e6;
        let vm = pts.von_mises_stress();
        assert!(
            vm[0].abs() < 1e-6,
            "hydrostatic stress has zero Von Mises: {}",
            vm[0]
        );
    }

    #[test]
    fn viscoplastic_rate_zero_in_elastic_domain() {
        let pts = SoaMaterialPoints::new(2);
        let yield_values = vec![-10.0, -1.0]; // below yield
        let rates = viscoplastic_rate_batch(&pts, &yield_values, 3.0, 1e-3);
        for &r in &rates {
            assert_eq!(r, 0.0, "viscoplastic rate should be zero below yield");
        }
    }

    #[test]
    fn batch_stats_counts_failed_points() {
        let mut pts = SoaMaterialPoints::new(3);
        pts.fatigue_damage[0] = 1.0; // failed
        pts.fatigue_damage[1] = 0.5; // partial
        let yv = vec![0.0; 3];
        let stats = MaterialBatchStats::compute(&pts, &yv);
        assert_eq!(stats.failed_count, 1);
        assert_eq!(stats.n, 3);
    }
}
