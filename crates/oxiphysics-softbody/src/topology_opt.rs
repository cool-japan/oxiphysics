// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topology optimisation for soft-body / deformable structures.
//!
//! Implements the Solid Isotropic Material with Penalisation (SIMP) method
//! together with the optimality-criteria (OC) update rule, sensitivity
//! filtering, Heaviside projection, and basic manufacturing constraints.
//!
//! The implementation is intentionally self-contained (no nalgebra, only
//! plain `f64` arrays and `Vec`f64`).

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Clamp a value to [lo, hi].
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// Euclidean distance between two 2-D grid points.
#[inline]
fn dist2d(xi: usize, yi: usize, xj: usize, yj: usize) -> f64 {
    let dx = xi as f64 - xj as f64;
    let dy = yi as f64 - yj as f64;
    (dx * dx + dy * dy).sqrt()
}

// ---------------------------------------------------------------------------
// DensityField
// ---------------------------------------------------------------------------

/// Design variable field ρ ∈ [0, 1] defined on a 2-D Cartesian grid.
///
/// Each element e has a density `rho\[e\]`.  The SIMP penalisation gives the
/// effective Young's modulus  E(ρ) = E₀ · ρᵖ.
#[derive(Debug, Clone)]
pub struct DensityField {
    /// Number of elements in the x-direction.
    pub nx: usize,
    /// Number of elements in the y-direction.
    pub ny: usize,
    /// Density values, row-major (index = iy * nx + ix).
    pub rho: Vec<f64>,
    /// SIMP penalisation exponent p (typically 3).
    pub penalty: f64,
    /// Baseline Young's modulus E₀ (Pa).
    pub e0: f64,
    /// Minimum modulus for void elements (prevents singularity).
    pub e_min: f64,
}

impl DensityField {
    /// Create a uniform density field with initial density `rho0`.
    pub fn new(nx: usize, ny: usize, rho0: f64, penalty: f64, e0: f64, e_min: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            rho: vec![rho0.clamp(0.0, 1.0); n],
            penalty,
            e0,
            e_min,
        }
    }

    /// Total number of elements.
    pub fn num_elements(&self) -> usize {
        self.nx * self.ny
    }

    /// SIMP effective modulus for element `e`.
    pub fn effective_modulus(&self, e: usize) -> f64 {
        self.e_min + (self.e0 - self.e_min) * self.rho[e].powf(self.penalty)
    }

    /// Derivative of SIMP modulus w.r.t. density ∂E/∂ρ for element `e`.
    pub fn dmodulus_drho(&self, e: usize) -> f64 {
        (self.e0 - self.e_min) * self.penalty * self.rho[e].powf(self.penalty - 1.0)
    }

    /// Current volume fraction V = (Σ ρ) / N.
    pub fn volume_fraction(&self) -> f64 {
        let sum: f64 = self.rho.iter().sum();
        sum / self.rho.len() as f64
    }

    /// Clamp all densities to [rho_min, rho_max].
    pub fn clamp_densities(&mut self, rho_min: f64, rho_max: f64) {
        for r in &mut self.rho {
            *r = clamp(*r, rho_min, rho_max);
        }
    }

    /// Element index from (ix, iy) grid coordinates.
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }
}

// ---------------------------------------------------------------------------
// ComplianceObjective
// ---------------------------------------------------------------------------

/// Structural compliance objective C = F^T u.
///
/// For a given density field and element stiffness sensitivities this
/// structure stores the total compliance and the per-element sensitivities.
#[derive(Debug, Clone)]
pub struct ComplianceObjective {
    /// Total structural compliance C = Σ_e u_e^T K_e u_e.
    pub compliance: f64,
    /// Per-element sensitivity ∂C/∂ρ_e (always ≤ 0 for maximising stiffness).
    pub sensitivity: Vec<f64>,
}

impl ComplianceObjective {
    /// Construct with a known compliance and sensitivity vector.
    pub fn new(compliance: f64, sensitivity: Vec<f64>) -> Self {
        Self {
            compliance,
            sensitivity,
        }
    }

    /// Compute sensitivity from element strain energy density `ue_ke_ue` (u_e^T K_e u_e
    /// for each element, assumed pre-computed by the FEA back-end) and the
    /// SIMP density field.
    ///
    /// ∂C/∂ρ_e = -p · E₀ · ρ_e^{p−1} · u_e^T k_e u_e / E_e
    ///          = -(p / ρ_e) · u_e^T K_e u_e
    pub fn from_strain_energy(density: &DensityField, ue_ke_ue: &[f64]) -> Self {
        assert_eq!(density.num_elements(), ue_ke_ue.len());
        let n = density.num_elements();
        let mut compliance = 0.0_f64;
        let mut sensitivity = vec![0.0_f64; n];
        for e in 0..n {
            let rho = density.rho[e].max(1e-12);
            compliance += density.effective_modulus(e) / density.e0 * ue_ke_ue[e];
            sensitivity[e] = -(density.penalty / rho) * ue_ke_ue[e];
        }
        Self {
            compliance,
            sensitivity,
        }
    }
}

// ---------------------------------------------------------------------------
// VolumeConstraint
// ---------------------------------------------------------------------------

/// Volume fraction constraint  V* = V_target.
///
/// Tracks whether the current design satisfies the target volume fraction and
/// provides Lagrange-multiplier estimates for the OC update.
#[derive(Debug, Clone, Copy)]
pub struct VolumeConstraint {
    /// Target volume fraction V* ∈ (0, 1].
    pub v_target: f64,
}

impl VolumeConstraint {
    /// Create a volume constraint with target fraction `v_target`.
    pub fn new(v_target: f64) -> Self {
        Self {
            v_target: v_target.clamp(1e-4, 1.0),
        }
    }

    /// Constraint violation: g = V_current − V_target (> 0 means over-filled).
    pub fn violation(&self, density: &DensityField) -> f64 {
        density.volume_fraction() - self.v_target
    }

    /// Returns `true` when the constraint is satisfied within tolerance.
    pub fn is_satisfied(&self, density: &DensityField, tol: f64) -> bool {
        self.violation(density).abs() <= tol
    }

    /// Bisection search for the Lagrange multiplier λ such that the OC update
    /// produces exactly the target volume fraction.
    ///
    /// `sensitivity` is ∂C/∂ρ per element (negative values expected).
    pub fn bisect_lambda(
        &self,
        density: &DensityField,
        sensitivity: &[f64],
        move_limit: f64,
    ) -> f64 {
        let mut lo = 1e-10_f64;
        let mut hi = 1e10_f64;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let vf = trial_volume(density, sensitivity, mid, move_limit);
            if vf > self.v_target {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / (lo + hi + 1e-30) < 1e-12 {
                break;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Helper: volume fraction produced by a trial Lagrange multiplier λ.
fn trial_volume(density: &DensityField, sensitivity: &[f64], lambda: f64, move_limit: f64) -> f64 {
    let n = density.num_elements();
    let sum: f64 = density
        .rho
        .iter()
        .zip(sensitivity.iter())
        .map(|(&rho, &sens)| {
            let s = (-sens / lambda).sqrt();
            clamp(s * rho, rho - move_limit, rho + move_limit).clamp(1e-3, 1.0)
        })
        .sum();
    sum / n as f64
}

// ---------------------------------------------------------------------------
// OcMethod
// ---------------------------------------------------------------------------

/// Optimality-criteria density update.
///
/// ρ_new = clamp( ρ · √(−∂C/∂ρ / λ), [ρ−m, ρ+m] ∩ [ρ_min, 1] )
///
/// where λ is the Lagrange multiplier from the volume constraint bisection
/// and m is the move limit.
#[derive(Debug, Clone, Copy)]
pub struct OcMethod {
    /// Move limit m (maximum allowed density change per iteration).
    pub move_limit: f64,
    /// Minimum density ρ_min (prevents void singularity).
    pub rho_min: f64,
}

impl OcMethod {
    /// Construct with given move limit and minimum density.
    pub fn new(move_limit: f64, rho_min: f64) -> Self {
        Self {
            move_limit,
            rho_min,
        }
    }

    /// Apply the OC update to `density` in-place, returning the new densities.
    pub fn update(&self, density: &mut DensityField, sensitivity: &[f64], lambda: f64) -> Vec<f64> {
        let n = density.num_elements();
        let mut new_rho = vec![0.0_f64; n];
        for e in 0..n {
            let rho = density.rho[e];
            let b = (-sensitivity[e] / lambda).max(0.0).sqrt();
            let rho_new = clamp(b * rho, rho - self.move_limit, rho + self.move_limit)
                .clamp(self.rho_min, 1.0);
            new_rho[e] = rho_new;
            density.rho[e] = rho_new;
        }
        new_rho
    }

    /// Maximum density change between two consecutive iterates.
    pub fn max_change(old_rho: &[f64], new_rho: &[f64]) -> f64 {
        old_rho
            .iter()
            .zip(new_rho.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// FilterSensitivity
// ---------------------------------------------------------------------------

/// Linear density filter for sensitivities.
///
/// Smooths the raw sensitivity field with a linear weight kernel based on
/// distance:  w(e, i) = max(0, r_min − dist(e, i)).
///
/// The filtered sensitivity is:
///   ∂C̃/∂ρ_e = (1/ρ_e) · Σ_i w(e,i) ρ_i (∂C/∂ρ_i) / Σ_i w(e,i)
#[derive(Debug, Clone, Copy)]
pub struct FilterSensitivity {
    /// Filter radius in element units.
    pub r_min: f64,
}

impl FilterSensitivity {
    /// Construct with filter radius `r_min` (in grid cells).
    pub fn new(r_min: f64) -> Self {
        Self { r_min }
    }

    /// Apply the density filter to the sensitivity field.
    ///
    /// Returns a new sensitivity vector of the same length.
    pub fn apply(&self, density: &DensityField, sensitivity: &[f64]) -> Vec<f64> {
        let nx = density.nx;
        let ny = density.ny;
        let n = nx * ny;
        let mut filtered = vec![0.0_f64; n];
        let r_ceil = self.r_min.ceil() as usize;

        for iy in 0..ny {
            for ix in 0..nx {
                let e = density.idx(ix, iy);
                let mut num = 0.0_f64;
                let mut denom = 0.0_f64;

                let iy_lo = iy.saturating_sub(r_ceil);
                let iy_hi = (iy + r_ceil + 1).min(ny);
                let ix_lo = ix.saturating_sub(r_ceil);
                let ix_hi = (ix + r_ceil + 1).min(nx);

                for jy in iy_lo..iy_hi {
                    for jx in ix_lo..ix_hi {
                        let d = dist2d(ix, iy, jx, jy);
                        if d < self.r_min {
                            let w = self.r_min - d;
                            let i = density.idx(jx, jy);
                            num += w * density.rho[i] * sensitivity[i];
                            denom += w * density.rho[i];
                        }
                    }
                }
                filtered[e] = if denom.abs() < 1e-30 {
                    sensitivity[e]
                } else {
                    num / (density.rho[e].max(1e-12) * denom) * density.rho[e]
                };
            }
        }
        filtered
    }
}

// ---------------------------------------------------------------------------
// TopOptHistory
// ---------------------------------------------------------------------------

/// History of compliance and volume fraction values per iteration.
#[derive(Debug, Clone, Default)]
pub struct TopOptHistory {
    /// Compliance value at each iteration.
    pub compliance: Vec<f64>,
    /// Volume fraction at each iteration.
    pub volume: Vec<f64>,
    /// Maximum density change at each iteration.
    pub max_change: Vec<f64>,
}

impl TopOptHistory {
    /// Construct an empty history.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an iteration record.
    pub fn push(&mut self, compliance: f64, volume: f64, max_change: f64) {
        self.compliance.push(compliance);
        self.volume.push(volume);
        self.max_change.push(max_change);
    }

    /// Number of recorded iterations.
    pub fn len(&self) -> usize {
        self.compliance.len()
    }

    /// Returns `true` when no iterations have been recorded.
    pub fn is_empty(&self) -> bool {
        self.compliance.is_empty()
    }

    /// Check convergence: max_change < tol for the last iteration.
    pub fn is_converged(&self, tol: f64) -> bool {
        self.max_change.last().is_some_and(|&c| c < tol)
    }

    /// Compliance improvement over the last `window` iterations.
    pub fn compliance_improvement(&self, window: usize) -> f64 {
        let n = self.compliance.len();
        if n < 2 {
            return f64::MAX;
        }
        let w = window.min(n - 1);
        let old = self.compliance[n - 1 - w];
        let new = self.compliance[n - 1];
        (old - new).abs()
    }
}

// ---------------------------------------------------------------------------
// HeavisideProjection
// ---------------------------------------------------------------------------

/// Smooth Heaviside projection for crisp 0/1 densities.
///
/// H_β(ρ) = (tanh(β·η) + tanh(β·(ρ−η))) / (tanh(β·η) + tanh(β·(1−η)))
///
/// As β → ∞ the projection becomes a step function at threshold η.
#[derive(Debug, Clone, Copy)]
pub struct HeavisideProjection {
    /// Sharpness parameter β (start low, ramp up during optimisation).
    pub beta: f64,
    /// Threshold η ∈ (0, 1); densities below η → 0, above → 1.
    pub eta: f64,
}

impl HeavisideProjection {
    /// Construct a Heaviside projection.
    pub fn new(beta: f64, eta: f64) -> Self {
        Self {
            beta,
            eta: eta.clamp(1e-4, 1.0 - 1e-4),
        }
    }

    /// Evaluate the smooth Heaviside at `rho`.
    pub fn project(&self, rho: f64) -> f64 {
        let b = self.beta;
        let e = self.eta;
        let denom = (b * e).tanh() + (b * (1.0 - e)).tanh();
        if denom.abs() < 1e-30 {
            return rho;
        }
        ((b * e).tanh() + (b * (rho - e)).tanh()) / denom
    }

    /// Derivative dH/dρ.
    pub fn derivative(&self, rho: f64) -> f64 {
        let b = self.beta;
        let e = self.eta;
        let denom = (b * e).tanh() + (b * (1.0 - e)).tanh();
        if denom.abs() < 1e-30 {
            return 1.0;
        }
        let sech2 = 1.0 - (b * (rho - e)).tanh().powi(2);
        b * sech2 / denom
    }

    /// Apply projection to every element in a density field (in-place).
    pub fn apply_field(&self, density: &mut DensityField) {
        for r in &mut density.rho {
            *r = self.project(*r);
        }
    }

    /// Increase β by multiplier (continuation strategy).
    pub fn increase_beta(&mut self, factor: f64) {
        self.beta *= factor;
    }
}

// ---------------------------------------------------------------------------
// ManufacturingConstraint
// ---------------------------------------------------------------------------

/// Manufacturing constraints: minimum length scale and overhang angle.
///
/// After convergence, elements below the minimum length scale are snapped to
/// void; elements violating the overhang constraint are penalised.
#[derive(Debug, Clone, Copy)]
pub struct ManufacturingConstraint {
    /// Minimum solid length scale (in grid cells); features narrower than
    /// this are removed.
    pub min_length_scale: f64,
    /// Maximum allowable overhang angle from vertical (degrees).
    pub max_overhang_deg: f64,
}

impl ManufacturingConstraint {
    /// Construct manufacturing constraints.
    pub fn new(min_length_scale: f64, max_overhang_deg: f64) -> Self {
        Self {
            min_length_scale,
            max_overhang_deg: max_overhang_deg.clamp(0.0, 90.0),
        }
    }

    /// Remove features smaller than `min_length_scale` by eroding then
    /// dilating the binary (thresholded at 0.5) density field.
    ///
    /// Returns the number of elements modified.
    pub fn enforce_min_length(&self, density: &mut DensityField, threshold: f64) -> usize {
        let nx = density.nx;
        let ny = density.ny;
        let n = nx * ny;
        let r = self.min_length_scale as usize + 1;
        let binary: Vec<f64> = density
            .rho
            .iter()
            .map(|&v| if v >= threshold { 1.0 } else { 0.0 })
            .collect();

        // Erosion: element is solid only if all neighbours within r are solid.
        let mut eroded = vec![0.0_f64; n];
        for iy in 0..ny {
            for ix in 0..nx {
                let mut all_solid = true;
                'outer: for jy in iy.saturating_sub(r)..((iy + r + 1).min(ny)) {
                    for jx in ix.saturating_sub(r)..((ix + r + 1).min(nx)) {
                        if dist2d(ix, iy, jx, jy) <= self.min_length_scale
                            && binary[jy * nx + jx] < 0.5
                        {
                            all_solid = false;
                            break 'outer;
                        }
                    }
                }
                eroded[iy * nx + ix] = if all_solid { 1.0 } else { 0.0 };
            }
        }

        // Dilation: element is solid if any neighbour within r of eroded is solid.
        let mut dilated = vec![0.0_f64; n];
        for iy in 0..ny {
            for ix in 0..nx {
                let mut any_solid = false;
                'outer2: for jy in iy.saturating_sub(r)..((iy + r + 1).min(ny)) {
                    for jx in ix.saturating_sub(r)..((ix + r + 1).min(nx)) {
                        if dist2d(ix, iy, jx, jy) <= self.min_length_scale
                            && eroded[jy * nx + jx] > 0.5
                        {
                            any_solid = true;
                            break 'outer2;
                        }
                    }
                }
                dilated[iy * nx + ix] = if any_solid { 1.0 } else { 0.0 };
            }
        }

        // Count changes and apply.
        let count = dilated
            .iter()
            .zip(density.rho.iter())
            .filter(|&(&d, &r)| (d - r).abs() > 0.5)
            .count();
        for (rho, &d) in density.rho.iter_mut().zip(dilated.iter()) {
            *rho = d;
        }
        count
    }

    /// Penalise overhanging elements: reduce density of elements whose
    /// support angle exceeds `max_overhang_deg`.
    ///
    /// "Build direction" is assumed to be +Y (upward).
    /// Returns the number of penalised elements.
    pub fn enforce_overhang(&self, density: &mut DensityField, threshold: f64) -> usize {
        let nx = density.nx;
        let ny = density.ny;
        let max_tan = self.max_overhang_deg.to_radians().tan();
        let mut count = 0;

        for iy in 1..ny {
            for ix in 0..nx {
                let e = density.idx(ix, iy);
                if density.rho[e] < threshold {
                    continue;
                }
                // Check support from row below.
                let support = if ix > 0 {
                    density.rho[density.idx(ix - 1, iy - 1)]
                        .max(density.rho[density.idx(ix, iy - 1)])
                } else {
                    density.rho[density.idx(ix, iy - 1)]
                };
                // Overhang check: if no support below and angle exceeded.
                if support < threshold {
                    // Estimate angle from horizontal span needed.
                    let span_x = 1.0_f64; // one cell horizontal
                    let span_y = 1.0_f64; // one cell vertical
                    if span_x / span_y > max_tan {
                        density.rho[e] *= 0.5; // penalise
                        count += 1;
                    }
                }
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// TopOptSolver
// ---------------------------------------------------------------------------

/// Full topology optimisation solver.
///
/// Orchestrates the FEA sensitivity computation → filtering → OC update →
/// convergence check loop.  The FEA back-end is provided as a closure so
/// that the solver remains physics-agnostic.
#[derive(Debug, Clone)]
pub struct TopOptSolver {
    /// OC update method.
    pub oc: OcMethod,
    /// Sensitivity filter.
    pub filter: FilterSensitivity,
    /// Volume constraint.
    pub volume: VolumeConstraint,
    /// Convergence tolerance on max density change.
    pub tol: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
}

impl TopOptSolver {
    /// Construct a solver with default SIMP parameters.
    pub fn new(v_target: f64, r_min: f64, move_limit: f64, tol: f64, max_iter: usize) -> Self {
        Self {
            oc: OcMethod::new(move_limit, 1e-3),
            filter: FilterSensitivity::new(r_min),
            volume: VolumeConstraint::new(v_target),
            tol,
            max_iter,
        }
    }

    /// Run the optimisation loop.
    ///
    /// `fea_fn(density) -> (compliance, ue_ke_ue)` is a closure that performs
    /// one forward FEA and returns the total compliance and per-element strain
    /// energy densities.
    pub fn optimise<F>(&self, density: &mut DensityField, mut fea_fn: F) -> TopOptHistory
    where
        F: FnMut(&DensityField) -> (f64, Vec<f64>),
    {
        let mut history = TopOptHistory::new();
        let mut prev_rho = density.rho.clone();

        for _iter in 0..self.max_iter {
            // FEA step.
            let (compliance, ue_ke_ue) = fea_fn(density);

            // Sensitivity analysis.
            let obj = ComplianceObjective::from_strain_energy(density, &ue_ke_ue);

            // Sensitivity filtering.
            let filtered_sens = self.filter.apply(density, &obj.sensitivity);

            // Lagrange multiplier from volume constraint bisection.
            let lambda = self
                .volume
                .bisect_lambda(density, &filtered_sens, self.oc.move_limit);

            // OC density update.
            let new_rho = self.oc.update(density, &filtered_sens, lambda);

            // Convergence metrics.
            let max_chg = OcMethod::max_change(&prev_rho, &new_rho);
            let vf = density.volume_fraction();
            history.push(compliance, vf, max_chg);

            prev_rho = new_rho;

            if max_chg < self.tol {
                break;
            }
        }
        history
    }
}

// ---------------------------------------------------------------------------
// TopOptVisualizer
// ---------------------------------------------------------------------------

/// Export a density field as a grayscale ASCII grid.
///
/// White (255) = solid (ρ = 1), black (0) = void (ρ = 0).
pub struct TopOptVisualizer;

impl TopOptVisualizer {
    /// Convert density field to a flat vector of u8 pixel values (row-major,
    /// top row first).  Each element maps to a single pixel.
    pub fn to_grayscale(density: &DensityField) -> Vec<u8> {
        let mut pixels = Vec::with_capacity(density.rho.len());
        // Rows in the density field are stored iy=0 at bottom; flip for image.
        for iy in (0..density.ny).rev() {
            for ix in 0..density.nx {
                let e = density.idx(ix, iy);
                let grey = (density.rho[e].clamp(0.0, 1.0) * 255.0).round() as u8;
                pixels.push(grey);
            }
        }
        pixels
    }

    /// Render density field to an ASCII art string.
    ///
    /// Uses `' '` for void and `'#'` for solid (threshold 0.5).
    pub fn to_ascii(density: &DensityField, threshold: f64) -> String {
        let mut out = String::with_capacity(density.nx * density.ny + density.ny);
        for iy in (0..density.ny).rev() {
            for ix in 0..density.nx {
                let e = density.idx(ix, iy);
                out.push(if density.rho[e] >= threshold {
                    '#'
                } else {
                    ' '
                });
            }
            out.push('\n');
        }
        out
    }

    /// Export density field as a simple PGM header + pixel bytes string.
    pub fn to_pgm_bytes(density: &DensityField) -> Vec<u8> {
        let header = format!("P5\n{} {}\n255\n", density.nx, density.ny);
        let mut out: Vec<u8> = header.into_bytes();
        out.extend(Self::to_grayscale(density));
        out
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // DensityField
    // -----------------------------------------------------------------------

    // 1. Uniform initialisation.
    #[test]
    fn test_density_field_uniform_init() {
        let df = DensityField::new(10, 10, 0.5, 3.0, 1.0, 1e-3);
        assert_eq!(df.num_elements(), 100);
        for &r in &df.rho {
            assert!((r - 0.5).abs() < 1e-12);
        }
    }

    // 2. Effective modulus at ρ=1 equals E0.
    #[test]
    fn test_effective_modulus_solid() {
        let mut df = DensityField::new(2, 2, 1.0, 3.0, 200e9, 1e3);
        df.rho[0] = 1.0;
        assert!((df.effective_modulus(0) - 200e9).abs() < 1.0);
    }

    // 3. Effective modulus at ρ=0 is approximately e_min.
    #[test]
    fn test_effective_modulus_void() {
        let mut df = DensityField::new(2, 2, 0.0, 3.0, 200e9, 1e3);
        df.rho[0] = 0.0;
        assert!((df.effective_modulus(0) - 1e3).abs() < 1.0);
    }

    // 4. SIMP monotonicity: higher density → higher modulus.
    #[test]
    fn test_simp_monotone() {
        let mut df = DensityField::new(1, 1, 0.0, 3.0, 1.0, 0.0);
        df.rho[0] = 0.0;
        let mut prev = df.effective_modulus(0);
        for i in 1..=10 {
            df.rho[0] = i as f64 / 10.0;
            let curr = df.effective_modulus(0);
            assert!(
                curr >= prev,
                "Modulus not monotone at rho={}: curr={curr} prev={prev}",
                df.rho[0]
            );
            prev = curr;
        }
    }

    // 5. Volume fraction of uniform field = initial rho.
    #[test]
    fn test_volume_fraction_uniform() {
        let df = DensityField::new(4, 5, 0.4, 3.0, 1.0, 0.0);
        assert!((df.volume_fraction() - 0.4).abs() < 1e-12);
    }

    // 6. idx returns correct flat index.
    #[test]
    fn test_idx() {
        let df = DensityField::new(5, 5, 0.5, 3.0, 1.0, 0.0);
        assert_eq!(df.idx(2, 3), 3 * 5 + 2);
    }

    // 7. clamp_densities enforces bounds.
    #[test]
    fn test_clamp_densities() {
        let mut df = DensityField::new(3, 3, 0.5, 3.0, 1.0, 0.0);
        df.rho[0] = -0.2;
        df.rho[1] = 1.8;
        df.clamp_densities(0.0, 1.0);
        assert!((df.rho[0] - 0.0).abs() < 1e-12);
        assert!((df.rho[1] - 1.0).abs() < 1e-12);
    }

    // 8. dmodulus_drho is positive for p > 1.
    #[test]
    fn test_dmodulus_drho_positive() {
        let mut df = DensityField::new(1, 1, 0.5, 3.0, 1.0, 0.0);
        df.rho[0] = 0.5;
        assert!(df.dmodulus_drho(0) > 0.0);
    }

    // -----------------------------------------------------------------------
    // ComplianceObjective
    // -----------------------------------------------------------------------

    // 9. Sensitivity from strain energy has correct sign.
    #[test]
    fn test_sensitivity_sign() {
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 1e-6);
        let ue_ke_ue = vec![1.0; 16];
        let obj = ComplianceObjective::from_strain_energy(&df, &ue_ke_ue);
        for &s in &obj.sensitivity {
            assert!(s < 0.0, "Sensitivity should be negative: {s}");
        }
    }

    // 10. Compliance is non-negative.
    #[test]
    fn test_compliance_nonneg() {
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 1e-6);
        let ue_ke_ue = vec![2.0; 16];
        let obj = ComplianceObjective::from_strain_energy(&df, &ue_ke_ue);
        assert!(obj.compliance >= 0.0);
    }

    // 11. Zero strain energy → zero compliance.
    #[test]
    fn test_zero_strain_energy() {
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let ue_ke_ue = vec![0.0; 16];
        let obj = ComplianceObjective::from_strain_energy(&df, &ue_ke_ue);
        assert!(obj.compliance.abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // VolumeConstraint
    // -----------------------------------------------------------------------

    // 12. Violation is zero for matching volume fraction.
    #[test]
    fn test_volume_violation_zero() {
        let vc = VolumeConstraint::new(0.5);
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        assert!(vc.violation(&df).abs() < 1e-10);
    }

    // 13. is_satisfied returns false when over-volume.
    #[test]
    fn test_volume_not_satisfied() {
        let vc = VolumeConstraint::new(0.3);
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        assert!(!vc.is_satisfied(&df, 0.01));
    }

    // 14. bisect_lambda returns positive value.
    #[test]
    fn test_bisect_lambda_positive() {
        let vc = VolumeConstraint::new(0.4);
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let sensitivity = vec![-1.0; 16];
        let lambda = vc.bisect_lambda(&df, &sensitivity, 0.2);
        assert!(lambda > 0.0);
    }

    // -----------------------------------------------------------------------
    // OcMethod
    // -----------------------------------------------------------------------

    // 15. OC update respects move limit.
    #[test]
    fn test_oc_move_limit() {
        let oc = OcMethod::new(0.2, 1e-3);
        let mut df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let old_rho = df.rho.clone();
        let sensitivity = vec![-1.0; 16];
        let lambda = 1.0;
        let new_rho = oc.update(&mut df, &sensitivity, lambda);
        for (o, n) in old_rho.iter().zip(new_rho.iter()) {
            assert!(
                (o - n).abs() <= 0.2 + 1e-12,
                "Move limit violated: |{o} - {n}| > 0.2"
            );
        }
    }

    // 16. OC update keeps densities in [rho_min, 1].
    #[test]
    fn test_oc_density_bounds() {
        let oc = OcMethod::new(0.5, 1e-3);
        let mut df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let sensitivity = vec![-0.001; 16];
        let lambda = 1e-8; // very small → B is huge → rho pushed to 1
        let new_rho = oc.update(&mut df, &sensitivity, lambda);
        for &r in &new_rho {
            assert!(r <= 1.0 + 1e-12, "Density exceeded 1: {r}");
            assert!(r >= 1e-3 - 1e-12, "Density below rho_min: {r}");
        }
    }

    // 17. max_change returns correct value.
    #[test]
    fn test_max_change() {
        let old = vec![0.5, 0.3, 0.7];
        let new = vec![0.6, 0.3, 0.5];
        let mc = OcMethod::max_change(&old, &new);
        assert!((mc - 0.2).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // FilterSensitivity
    // -----------------------------------------------------------------------

    // 18. Filtered sensitivity has same length as input.
    #[test]
    fn test_filter_length() {
        let flt = FilterSensitivity::new(1.5);
        let df = DensityField::new(5, 5, 0.5, 3.0, 1.0, 0.0);
        let sens = vec![-1.0; 25];
        let filtered = flt.apply(&df, &sens);
        assert_eq!(filtered.len(), 25);
    }

    // 19. Uniform sensitivity field stays uniform after filtering.
    #[test]
    fn test_filter_uniform() {
        let flt = FilterSensitivity::new(2.0);
        let df = DensityField::new(8, 8, 0.5, 3.0, 1.0, 0.0);
        let sens = vec![-2.0; 64];
        let filtered = flt.apply(&df, &sens);
        for &f in &filtered {
            assert!(
                (f - (-2.0)).abs() < 1e-8,
                "Uniform field should stay uniform, got {f}"
            );
        }
    }

    // 20. Filter with r_min < 1 acts as identity.
    #[test]
    fn test_filter_r_min_small() {
        let flt = FilterSensitivity::new(0.4); // less than 1 cell
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let sens: Vec<f64> = (0..16).map(|i| -(i as f64)).collect();
        let filtered = flt.apply(&df, &sens);
        // Each element is only influenced by itself since r_min < 1.
        for (i, (&s, &f)) in sens.iter().zip(filtered.iter()).enumerate() {
            assert!((f - s).abs() < 1e-8, "Element {i}: expected {s}, got {f}");
        }
    }

    // -----------------------------------------------------------------------
    // TopOptHistory
    // -----------------------------------------------------------------------

    // 21. Empty history is not converged.
    #[test]
    fn test_history_empty_not_converged() {
        let h = TopOptHistory::new();
        assert!(!h.is_converged(1e-3));
        assert!(h.is_empty());
    }

    // 22. Converged when last change is below tol.
    #[test]
    fn test_history_converged() {
        let mut h = TopOptHistory::new();
        h.push(1.0, 0.5, 1e-4);
        assert!(h.is_converged(1e-3));
    }

    // 23. Not converged when change is above tol.
    #[test]
    fn test_history_not_converged() {
        let mut h = TopOptHistory::new();
        h.push(1.0, 0.5, 0.1);
        assert!(!h.is_converged(1e-3));
    }

    // 24. compliance_improvement returns MAX for single entry.
    #[test]
    fn test_compliance_improvement_single() {
        let mut h = TopOptHistory::new();
        h.push(5.0, 0.4, 0.1);
        assert_eq!(h.compliance_improvement(3), f64::MAX);
    }

    // 25. compliance_improvement correct for multiple entries.
    #[test]
    fn test_compliance_improvement_multi() {
        let mut h = TopOptHistory::new();
        h.push(10.0, 0.5, 0.1);
        h.push(8.0, 0.5, 0.05);
        h.push(7.0, 0.5, 0.02);
        assert!((h.compliance_improvement(1) - 1.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // HeavisideProjection
    // -----------------------------------------------------------------------

    // 26. H(0) ≈ 0 for large beta.
    #[test]
    fn test_heaviside_zero() {
        let h = HeavisideProjection::new(50.0, 0.5);
        let val = h.project(0.0);
        assert!(val < 0.05, "H(0) should be near 0, got {val}");
    }

    // 27. H(1) ≈ 1 for large beta.
    #[test]
    fn test_heaviside_one() {
        let h = HeavisideProjection::new(50.0, 0.5);
        let val = h.project(1.0);
        assert!(val > 0.95, "H(1) should be near 1, got {val}");
    }

    // 28. H(eta) ≈ 0.5 (by construction).
    #[test]
    fn test_heaviside_midpoint() {
        let eta = 0.5;
        let h = HeavisideProjection::new(10.0, eta);
        let val = h.project(eta);
        assert!((val - 0.5).abs() < 1e-6, "H(eta) should be 0.5, got {val}");
    }

    // 29. Derivative is non-negative.
    #[test]
    fn test_heaviside_derivative_nonneg() {
        let h = HeavisideProjection::new(8.0, 0.5);
        for i in 0..=10 {
            let rho = i as f64 / 10.0;
            assert!(h.derivative(rho) >= 0.0);
        }
    }

    // 30. increase_beta multiplies beta.
    #[test]
    fn test_increase_beta() {
        let mut h = HeavisideProjection::new(2.0, 0.5);
        h.increase_beta(3.0);
        assert!((h.beta - 6.0).abs() < 1e-12);
    }

    // 31. apply_field projects all densities to [0, 1].
    #[test]
    fn test_heaviside_apply_field_bounded() {
        let h = HeavisideProjection::new(10.0, 0.5);
        let mut df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        h.apply_field(&mut df);
        for &r in &df.rho {
            assert!(
                (0.0..=1.0).contains(&r),
                "Projected density out of bounds: {r}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // ManufacturingConstraint
    // -----------------------------------------------------------------------

    // 32. enforce_min_length runs without panic.
    #[test]
    fn test_manufacturing_no_panic() {
        let mc = ManufacturingConstraint::new(1.5, 45.0);
        let mut df = DensityField::new(6, 6, 0.7, 3.0, 1.0, 0.0);
        let _changes = mc.enforce_min_length(&mut df, 0.5);
        // Just check no panic and densities are binary.
        for &r in &df.rho {
            assert!(r == 0.0 || r == 1.0, "Result should be binary: {r}");
        }
    }

    // 33. enforce_overhang runs without panic.
    #[test]
    fn test_overhang_no_panic() {
        let mc = ManufacturingConstraint::new(1.0, 45.0);
        let mut df = DensityField::new(5, 5, 0.8, 3.0, 1.0, 0.0);
        let _count = mc.enforce_overhang(&mut df, 0.5);
    }

    // 34. max_overhang_deg is clamped to [0, 90].
    #[test]
    fn test_overhang_angle_clamped() {
        let mc = ManufacturingConstraint::new(1.0, 120.0);
        assert_eq!(mc.max_overhang_deg, 90.0);
    }

    // -----------------------------------------------------------------------
    // TopOptSolver
    // -----------------------------------------------------------------------

    // 35. Solver with trivial FEA runs without panic and produces history.
    #[test]
    fn test_solver_trivial_fea() {
        let mut df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 1e-6);
        let solver = TopOptSolver::new(0.4, 1.5, 0.2, 1e-3, 20);
        let history = solver.optimise(&mut df, |_density| {
            let n = 16;
            (1.0, vec![1.0; n])
        });
        assert!(
            !history.is_empty(),
            "History should have at least one entry"
        );
    }

    // 36. Volume fraction approaches target after optimisation.
    #[test]
    fn test_solver_volume_target() {
        let target = 0.4;
        let mut df = DensityField::new(6, 6, 0.6, 3.0, 1.0, 1e-6);
        let solver = TopOptSolver::new(target, 1.5, 0.2, 1e-4, 50);
        let _history = solver.optimise(&mut df, |d| {
            let n = d.num_elements();
            let ue_ke_ue: Vec<f64> = (0..n).map(|e| d.effective_modulus(e)).collect();
            let compliance: f64 = ue_ke_ue.iter().sum();
            (compliance, ue_ke_ue)
        });
        let vf = df.volume_fraction();
        assert!(
            (vf - target).abs() < 0.1,
            "Volume fraction {vf} should be near target {target}"
        );
    }

    // -----------------------------------------------------------------------
    // TopOptVisualizer
    // -----------------------------------------------------------------------

    // 37. Grayscale output has correct length.
    #[test]
    fn test_visualizer_grayscale_length() {
        let df = DensityField::new(8, 6, 0.5, 3.0, 1.0, 0.0);
        let pixels = TopOptVisualizer::to_grayscale(&df);
        assert_eq!(pixels.len(), 48);
    }

    // 38. Solid field (ρ=1) → all pixels 255.
    #[test]
    fn test_visualizer_solid_white() {
        let df = DensityField::new(4, 4, 1.0, 3.0, 1.0, 0.0);
        let pixels = TopOptVisualizer::to_grayscale(&df);
        assert!(pixels.iter().all(|&p| p == 255));
    }

    // 39. Void field (ρ=0) → all pixels 0.
    #[test]
    fn test_visualizer_void_black() {
        let df = DensityField::new(4, 4, 0.0, 3.0, 1.0, 0.0);
        let pixels = TopOptVisualizer::to_grayscale(&df);
        assert!(pixels.iter().all(|&p| p == 0));
    }

    // 40. ASCII output has correct number of lines.
    #[test]
    fn test_visualizer_ascii_line_count() {
        let df = DensityField::new(5, 3, 0.5, 3.0, 1.0, 0.0);
        let ascii = TopOptVisualizer::to_ascii(&df, 0.5);
        assert_eq!(ascii.lines().count(), 3);
    }

    // 41. ASCII: solid elements produce '#', void produce ' '.
    #[test]
    fn test_visualizer_ascii_chars() {
        let mut df = DensityField::new(2, 1, 0.0, 3.0, 1.0, 0.0);
        df.rho[0] = 1.0;
        df.rho[1] = 0.0;
        let ascii = TopOptVisualizer::to_ascii(&df, 0.5);
        assert!(ascii.contains('#'));
        assert!(ascii.contains(' '));
    }

    // 42. PGM output starts with "P5".
    #[test]
    fn test_visualizer_pgm_header() {
        let df = DensityField::new(4, 4, 0.5, 3.0, 1.0, 0.0);
        let pgm = TopOptVisualizer::to_pgm_bytes(&df);
        assert!(pgm.starts_with(b"P5"), "PGM must start with P5 header");
    }

    // -----------------------------------------------------------------------
    // Integration-style tests
    // -----------------------------------------------------------------------

    // 43. OC update + volume constraint converges volume fraction.
    #[test]
    fn test_oc_volume_converge() {
        let target = 0.35;
        let mut df = DensityField::new(10, 10, 0.5, 3.0, 1.0, 1e-6);
        let vc = VolumeConstraint::new(target);
        let oc = OcMethod::new(0.2, 1e-3);
        let sensitivity = vec![-1.0; 100];

        for _ in 0..100 {
            let lambda = vc.bisect_lambda(&df, &sensitivity, 0.2);
            oc.update(&mut df, &sensitivity, lambda);
        }
        let vf = df.volume_fraction();
        assert!(
            (vf - target).abs() < 0.02,
            "Volume fraction {vf} should approach target {target}"
        );
    }

    // 44. DensityField with all ones → volume fraction = 1.
    #[test]
    fn test_volume_fraction_all_ones() {
        let df = DensityField::new(5, 5, 1.0, 3.0, 1.0, 0.0);
        assert!((df.volume_fraction() - 1.0).abs() < 1e-12);
    }

    // 45. DensityField with all zeros → volume fraction = 0.
    #[test]
    fn test_volume_fraction_all_zeros() {
        let df = DensityField::new(5, 5, 0.0, 3.0, 1.0, 0.0);
        assert!(df.volume_fraction().abs() < 1e-12);
    }

    // 46. Filter preserves sign of sensitivity.
    #[test]
    fn test_filter_preserves_sign() {
        let flt = FilterSensitivity::new(1.5);
        let df = DensityField::new(5, 5, 0.5, 3.0, 1.0, 0.0);
        let sens = vec![-1.0; 25];
        let filtered = flt.apply(&df, &sens);
        for &f in &filtered {
            assert!(f <= 0.0, "Filtered sensitivity should be non-positive: {f}");
        }
    }

    // 47. Heaviside projection is monotone.
    #[test]
    fn test_heaviside_monotone() {
        let h = HeavisideProjection::new(8.0, 0.5);
        let mut prev = h.project(0.0);
        for i in 1..=20 {
            let rho = i as f64 / 20.0;
            let curr = h.project(rho);
            assert!(curr >= prev - 1e-12, "Heaviside must be monotone");
            prev = curr;
        }
    }

    // 48. SIMP derivative for p=1 equals (E0 - E_min).
    #[test]
    fn test_simp_derivative_p1() {
        let e0 = 100.0;
        let e_min = 0.1;
        let mut df = DensityField::new(1, 1, 0.5, 1.0, e0, e_min);
        df.rho[0] = 0.7; // ρ^(p-1) = ρ^0 = 1
        let expected = e0 - e_min;
        assert!((df.dmodulus_drho(0) - expected).abs() < 1e-8);
    }

    // 49. TopOptHistory len matches number of push calls.
    #[test]
    fn test_history_len() {
        let mut h = TopOptHistory::new();
        for i in 0..7 {
            h.push(i as f64, 0.5, 0.01 * i as f64);
        }
        assert_eq!(h.len(), 7);
    }

    // 50. Solver with fully-solid design and zero move limit stays at rho=1.
    #[test]
    fn test_solver_no_update_at_max_density() {
        let mut df = DensityField::new(4, 4, 1.0, 3.0, 1.0, 1e-6);
        let solver = TopOptSolver::new(1.0, 1.5, 0.0, 1e-8, 5);
        let _h = solver.optimise(&mut df, |_d| (1.0, vec![1.0; 16]));
        for &r in &df.rho {
            assert!(
                (r - 1.0).abs() < 0.01,
                "Dense field with zero move limit should stay near 1, got {r}"
            );
        }
    }
}
