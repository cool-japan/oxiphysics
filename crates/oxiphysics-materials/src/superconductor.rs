// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Superconductor material models.
//!
//! Implements Ginzburg–Landau phenomenological theory, BCS gap equations,
//! London electrodynamics, flux-vortex lattice mechanics, and Type I / II
//! classification for superconducting materials.
//!
//! All quantities use SI units unless explicitly noted.

use crate::{Error, Result, require_finite, require_positive};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Flux quantum Φ₀ = h/(2e) (Wb).
pub const FLUX_QUANTUM: f64 = 2.067833848e-15;
/// Boltzmann constant k_B (J/K).
const K_B: f64 = 1.380649e-23;
/// Threshold: kappa > 1/√2 → Type II.
const KAPPA_THRESHOLD: f64 = std::f64::consts::FRAC_1_SQRT_2;

// ---------------------------------------------------------------------------
// Superconductor type
// ---------------------------------------------------------------------------

/// Classification of superconducting materials.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SuperconductorType {
    /// Type I: complete Meissner effect, single critical field H_c.
    TypeI,
    /// Type II: mixed Abrikosov state between H_c1 and H_c2.
    TypeII,
    /// High-temperature cuprate / unconventional superconductor (extreme Type II).
    HighTc,
}

// ---------------------------------------------------------------------------
// SuperconductorProps
// ---------------------------------------------------------------------------

/// Intrinsic physical properties of a superconducting material.
#[derive(Clone, Debug)]
pub struct SuperconductorProps {
    /// Critical temperature T_c (K) below which superconductivity exists.
    pub critical_temperature: f64,
    /// Lower critical field H_c1 (A/m) – vortex entry field (Type II).
    pub critical_field_hc1: f64,
    /// Upper critical field H_c2 (A/m) – superconductivity destroyed above this.
    pub critical_field_hc2: f64,
    /// London penetration depth λ (m): characteristic magnetic field screening length.
    pub london_penetration_depth: f64,
    /// Ginzburg–Landau coherence length ξ (m): spatial scale of order parameter.
    pub coherence_length: f64,
    /// Ginzburg–Landau parameter κ = λ/ξ (dimensionless).
    pub kappa: f64,
}

impl SuperconductorProps {
    /// Construct a `SuperconductorProps` with all fields specified.
    pub fn new(
        critical_temperature: f64,
        critical_field_hc1: f64,
        critical_field_hc2: f64,
        london_penetration_depth: f64,
        coherence_length: f64,
    ) -> Self {
        let kappa = london_penetration_depth / coherence_length.max(1e-300);
        Self {
            critical_temperature,
            critical_field_hc1,
            critical_field_hc2,
            london_penetration_depth,
            coherence_length,
            kappa,
        }
    }

    /// Example Nb (niobium) properties at T → 0.
    pub fn niobium() -> Self {
        Self::new(9.25, 1.58e5, 2.44e5, 39e-9, 38e-9)
    }

    /// Example Pb (lead) – Type I superconductor.
    pub fn lead() -> Self {
        // Lead is Type I: κ ≈ 0.48 < 1/√2.
        Self::new(7.19, 6.4e4, 6.4e4, 37e-9, 83e-9)
    }

    /// Example YBCO (yttrium barium copper oxide) – HiTc superconductor.
    pub fn ybco() -> Self {
        Self::new(93.0, 2.0e7, 1.2e8, 150e-9, 1.5e-9)
    }

    /// Return the classification of this material.
    pub fn superconductor_type(&self) -> SuperconductorType {
        if self.kappa > 10.0 {
            SuperconductorType::HighTc
        } else if is_type_ii(self.kappa) {
            SuperconductorType::TypeII
        } else {
            SuperconductorType::TypeI
        }
    }
}

// ---------------------------------------------------------------------------
// Ginzburg–Landau model
// ---------------------------------------------------------------------------

/// Largest Newton update `max_i |δψ_i|` at which [`GinzburgLandauModel::solve_gl`]
/// accepts the discrete Ginzburg–Landau equation as solved.
const GL_NEWTON_TOLERANCE: f64 = 1.0e-12;
/// Smallest fraction of a Newton step the backtracking line search may fall back to.
const GL_MIN_STEP_FRACTION: f64 = 1.0e-10;

/// Ginzburg–Landau phenomenological model for an order parameter field.
///
/// The order parameter |ψ|² represents the local superfluid density
/// (|ψ|² = 1 in the fully superconducting state, 0 in the normal state).
/// The profile is defined on a uniform grid of `n` points with spacing `dx`
/// (metres), `x_i = i·dx`, so the coherence length ξ and the grid spacing
/// enter the discrete equation only through the ratio `ξ/dx`.
#[derive(Clone, Debug)]
pub struct GinzburgLandauModel {
    /// Spatial grid of order parameter magnitudes |ψ| (dimensionless, 0–1).
    pub order_parameter: Vec<f64>,
    /// Corresponding magnetic field values B (T) on the same grid.
    pub magnetic_field: Vec<f64>,
    /// Grid spacing Δx (m) between consecutive points of the profile.
    pub dx: f64,
    /// Ginzburg–Landau coherence length ξ (m).
    pub coherence_length: f64,
    /// London penetration depth λ (m).
    pub penetration_depth: f64,
    /// GL parameter κ = λ/ξ.
    pub kappa: f64,
}

impl GinzburgLandauModel {
    /// Create a new GL model with `n` grid points of spacing `dx` (m).
    pub fn new(n: usize, dx: f64, coherence_length: f64, penetration_depth: f64) -> Self {
        let kappa = penetration_depth / coherence_length.max(1e-300);
        Self {
            order_parameter: vec![1.0; n],
            magnetic_field: vec![0.0; n],
            dx,
            coherence_length,
            penetration_depth,
            kappa,
        }
    }

    /// Solve the stationary, zero-field Ginzburg–Landau equation
    /// `ξ² ψ'' + ψ (1 − ψ²) = 0` with a Newton iteration and return the
    /// order-parameter profile `|ψ(x)|`.
    ///
    /// The boundary conditions are `|ψ(0)| = 0` (surface) and
    /// `|ψ(L)| = 1` (bulk) with `L = (n − 1)·dx`.  Second-order central
    /// differences with `r = (ξ/dx)²` give the nonlinear system
    /// `R_i = r (ψ_{i−1} − 2ψ_i + ψ_{i+1}) + ψ_i − ψ_i³ = 0`, whose
    /// tridiagonal Jacobian (off-diagonals `r`, diagonal `−2r + 1 − 3ψ_i²`) is
    /// solved with the Thomas algorithm at every iteration; a backtracking line
    /// search on `max|R|` keeps the iteration globally convergent.  The solution
    /// approaches the analytic kink `ψ(x) = tanh(x / (√2 ξ))` with an error
    /// of order `(dx/ξ)²`.
    ///
    /// The current `order_parameter` (default: all ones) is the starting guess;
    /// on success it is replaced by the converged profile.
    ///
    /// `n_iter` is the maximum number of Newton iterations; the iteration stops
    /// as soon as the largest update falls below `1e-12`.
    ///
    /// # Errors
    /// * [`Error::ParameterOutOfRange`] / [`Error::General`] for a non-positive
    ///   or non-finite `dx` or `coherence_length`.
    /// * [`Error::Singularity`] if the Newton system has a zero pivot.
    /// * [`Error::NumericalDivergence`] if the iteration does not converge
    ///   within `n_iter` iterations or the line search stalls.  The stored
    ///   profile is left unchanged in every error case.
    pub fn solve_gl(&mut self, n_iter: usize) -> Result<Vec<f64>> {
        require_finite("dx", self.dx)?;
        require_positive("dx", self.dx)?;
        require_finite("coherence_length", self.coherence_length)?;
        require_positive("coherence_length", self.coherence_length)?;
        let n = self.order_parameter.len();
        if n < 2 {
            return Ok(self.order_parameter.clone());
        }
        let ratio = self.coherence_length / self.dx;
        let r = ratio * ratio;
        require_finite("(coherence_length / dx)^2", r)?;

        // Starting guess with the Dirichlet conditions imposed.
        let mut psi = self.order_parameter.clone();
        for value in psi.iter_mut() {
            if !value.is_finite() {
                *value = 1.0;
            }
        }
        psi[0] = 0.0;
        psi[n - 1] = 1.0;
        let m = n - 2;
        if m == 0 {
            self.order_parameter = psi.clone();
            return Ok(psi);
        }

        let mut residual = vec![0.0_f64; m];
        let mut trial_residual = vec![0.0_f64; m];
        let mut diagonal = vec![0.0_f64; m];
        let mut delta = vec![0.0_f64; m];
        let mut scratch = vec![0.0_f64; m];
        let mut trial = psi.clone();
        let mut last_norm = f64::INFINITY;

        for _ in 0..n_iter {
            let norm = gl_residual(&psi, r, &mut residual);
            last_norm = norm;
            for (d, &p) in diagonal.iter_mut().zip(psi[1..n - 1].iter()) {
                *d = -2.0 * r + 1.0 - 3.0 * p * p;
            }
            for (d, &res) in delta.iter_mut().zip(residual.iter()) {
                *d = -res;
            }
            solve_constant_offdiagonal_tridiagonal(r, &diagonal, &mut delta, &mut scratch)?;
            let step = delta.iter().fold(0.0_f64, |acc, d| acc.max(d.abs()));
            if !step.is_finite() {
                return Err(Error::diverged("ginzburg_landau_newton", 0, norm));
            }
            if step < GL_NEWTON_TOLERANCE {
                for (p, d) in psi[1..n - 1].iter_mut().zip(delta.iter()) {
                    *p += d;
                }
                self.order_parameter = psi.clone();
                return Ok(psi);
            }
            // Backtracking line search on the max-norm of the residual.
            let mut fraction = 1.0_f64;
            loop {
                for ((t, &p), &d) in trial[1..n - 1]
                    .iter_mut()
                    .zip(psi[1..n - 1].iter())
                    .zip(delta.iter())
                {
                    *t = p + fraction * d;
                }
                let trial_norm = gl_residual(&trial, r, &mut trial_residual);
                if trial_norm <= (1.0 - 1.0e-4 * fraction) * norm {
                    break;
                }
                fraction *= 0.5;
                if fraction < GL_MIN_STEP_FRACTION {
                    return Err(Error::diverged("ginzburg_landau_newton", 0, norm));
                }
            }
            psi.copy_from_slice(&trial);
        }
        Err(Error::diverged("ginzburg_landau_newton", n_iter, last_norm))
    }
}

/// Residual `R_i = r (ψ_{i−1} − 2ψ_i + ψ_{i+1}) + ψ_i − ψ_i³` of the discrete
/// Ginzburg–Landau equation at the interior nodes of `psi`; returns `max|R_i|`.
fn gl_residual(psi: &[f64], r: f64, out: &mut [f64]) -> f64 {
    let mut norm = 0.0_f64;
    for (res, w) in out.iter_mut().zip(psi.windows(3)) {
        let p = w[1];
        *res = r * (w[0] - 2.0 * p + w[2]) + p - p * p * p;
        norm = norm.max(res.abs());
    }
    norm
}

/// Thomas algorithm for the tridiagonal system with main diagonal `diagonal`
/// and both off-diagonals equal to `off`.  The right-hand side `rhs` is
/// overwritten with the solution; `scratch` must have the same length.
fn solve_constant_offdiagonal_tridiagonal(
    off: f64,
    diagonal: &[f64],
    rhs: &mut [f64],
    scratch: &mut [f64],
) -> Result<()> {
    let m = diagonal.len();
    if m == 0 {
        return Ok(());
    }
    let zero_pivot = || {
        Error::singular(
            "ginzburg_landau_newton",
            "zero or non-finite pivot in the tridiagonal Newton system",
        )
    };
    let mut pivot = diagonal[0];
    if pivot == 0.0 || !pivot.is_finite() {
        return Err(zero_pivot());
    }
    scratch[0] = off / pivot;
    rhs[0] /= pivot;
    for i in 1..m {
        pivot = diagonal[i] - off * scratch[i - 1];
        if pivot == 0.0 || !pivot.is_finite() {
            return Err(zero_pivot());
        }
        scratch[i] = off / pivot;
        rhs[i] = (rhs[i] - off * rhs[i - 1]) / pivot;
    }
    for i in (0..m - 1).rev() {
        rhs[i] -= scratch[i] * rhs[i + 1];
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// BCS gap equation
// ---------------------------------------------------------------------------

/// BCS (Bardeen–Cooper–Schrieffer) energy gap model.
///
/// The BCS gap Δ(T) determines the binding energy of Cooper pairs.
#[derive(Clone, Debug)]
pub struct BcsGap {
    /// Zero-temperature gap energy Δ(0) (J).
    pub gap_energy: f64,
    /// Electronic density of states at the Fermi level N(0) (J⁻¹ m⁻³).
    pub density_of_states: f64,
    /// Critical temperature T_c (K).
    pub critical_temperature: f64,
}

impl BcsGap {
    /// Construct from zero-temperature gap, density of states, and T_c.
    pub fn new(gap_energy: f64, density_of_states: f64, critical_temperature: f64) -> Self {
        Self {
            gap_energy,
            density_of_states,
            critical_temperature,
        }
    }

    /// Construct from T_c using the BCS relation Δ(0) = 1.764 k_B T_c.
    pub fn from_tc(critical_temperature: f64, density_of_states: f64) -> Self {
        let gap = 1.764 * K_B * critical_temperature;
        Self::new(gap, density_of_states, critical_temperature)
    }

    /// Compute the temperature-dependent gap Δ(T) using the BCS interpolation:
    ///
    /// Δ(T) ≈ Δ(0) · tanh(1.74 √((T_c/T) - 1))  for T < T_c
    ///
    /// Returns 0.0 for T ≥ T_c.
    pub fn compute_gap(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return self.gap_energy;
        }
        if temperature >= self.critical_temperature {
            return 0.0;
        }
        let ratio = self.critical_temperature / temperature;
        let arg = 1.74 * (ratio - 1.0).sqrt();
        self.gap_energy * arg.tanh()
    }

    /// Normalised gap Δ(T)/Δ(0).
    pub fn normalised_gap(&self, temperature: f64) -> f64 {
        if self.gap_energy < 1e-300 {
            return 0.0;
        }
        self.compute_gap(temperature) / self.gap_energy
    }
}

// ---------------------------------------------------------------------------
// Flux vortex
// ---------------------------------------------------------------------------

/// A single Abrikosov flux vortex carrying one flux quantum Φ₀.
#[derive(Clone, Debug)]
pub struct FluxVortex {
    /// 2-D position (x, y) of the vortex core in metres.
    pub position: [f64; 2],
    /// Quantised magnetic flux carried by this vortex (Wb).
    ///
    /// For a singly quantised vortex this equals [`FLUX_QUANTUM`].
    pub flux_quantum: f64,
}

impl FluxVortex {
    /// Create a singly-quantised vortex at `position`.
    pub fn new(position: [f64; 2]) -> Self {
        Self {
            position,
            flux_quantum: FLUX_QUANTUM,
        }
    }

    /// Euclidean distance between this vortex and another.
    pub fn distance_to(&self, other: &FluxVortex) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        (dx * dx + dy * dy).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Vortex lattice
// ---------------------------------------------------------------------------

/// An Abrikosov vortex lattice: a collection of flux vortices in a Type II
/// superconductor at fields H_c1 < H < H_c2.
#[derive(Clone, Debug)]
pub struct VortexLattice {
    /// All vortices in this lattice.
    pub vortices: Vec<FluxVortex>,
    /// London penetration depth λ (m).
    pub penetration_depth: f64,
}

impl VortexLattice {
    /// Create an empty lattice.
    pub fn new(penetration_depth: f64) -> Self {
        Self {
            vortices: Vec::new(),
            penetration_depth,
        }
    }

    /// Add a vortex to the lattice.
    pub fn add_vortex(&mut self, vortex: FluxVortex) {
        self.vortices.push(vortex);
    }

    /// Build a triangular (Abrikosov) vortex lattice with spacing `a` over an
    /// `nx × ny` grid.
    pub fn triangular_lattice(nx: usize, ny: usize, spacing: f64, penetration_depth: f64) -> Self {
        let mut lattice = Self::new(penetration_depth);
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * spacing + if j % 2 == 1 { 0.5 * spacing } else { 0.0 };
                let y = j as f64 * spacing * (3.0_f64.sqrt() / 2.0);
                lattice.add_vortex(FluxVortex::new([x, y]));
            }
        }
        lattice
    }

    /// Compute the pairwise inter-vortex repulsive force on each vortex.
    ///
    /// The London interaction force between two vortices at separation r is:
    ///
    /// F(r) = (Φ₀²)/(2π μ₀ λ²) · K₁(r/λ) / r
    ///
    /// Here we use the approximation K₁(r/λ) ≈ (λ/r) · exp(-r/λ) for r > 0.
    ///
    /// Returns a `Vec` of 2-D force vectors `[fx, fy]` (N/m) per vortex.
    pub fn compute_inter_vortex_force(&self) -> Vec<[f64; 2]> {
        let n = self.vortices.len();
        let mut forces = vec![[0.0_f64; 2]; n];
        let mu0 = 4.0 * PI * 1e-7_f64;
        let lambda = self.penetration_depth;
        // Prefactor C = Φ₀²/(2π μ₀ λ²).
        let prefactor = FLUX_QUANTUM.powi(2) / (2.0 * PI * mu0 * lambda.powi(2));
        for (i, force_i) in forces.iter_mut().enumerate() {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.vortices[i].position[0] - self.vortices[j].position[0];
                let dy = self.vortices[i].position[1] - self.vortices[j].position[1];
                let r = (dx * dx + dy * dy).sqrt().max(1e-15);
                let xi = r / lambda;
                // K₁ approximation.
                let k1_approx = (1.0 / xi) * (-xi).exp();
                let mag = prefactor * k1_approx / r;
                force_i[0] += mag * dx / r;
                force_i[1] += mag * dy / r;
            }
        }
        forces
    }

    /// Total applied flux density B = n_v · Φ₀ (T·m²) where n_v is the
    /// vortex areal density.  Assumes vortices are distributed over area
    /// given by `cell_area` (m²).
    pub fn flux_density(&self, cell_area: f64) -> f64 {
        (self.vortices.len() as f64) * FLUX_QUANTUM / cell_area.max(1e-300)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// London equation screening factor for a uniform applied field.
///
/// The London equation ∇²B = B/λ² gives exponential screening:
///
/// B(x) = `field` · exp(−x/`lambda`)
///
/// This function returns the factor exp(−|field|/`lambda`) as a
/// dimensionless attenuation coefficient for normalised depth units.
pub fn london_equation_factor(lambda: f64, field: f64) -> f64 {
    if lambda <= 0.0 {
        return 0.0;
    }
    (-field.abs() / lambda).exp()
}

/// Returns `true` if the Ginzburg–Landau parameter κ = λ/ξ exceeds 1/√2,
/// indicating a Type II superconductor with an Abrikosov mixed state.
pub fn is_type_ii(kappa: f64) -> bool {
    kappa > KAPPA_THRESHOLD
}

/// Upper critical field H_c2 = Φ₀ / (2π μ₀ ξ²) (SI, A/m).
pub fn upper_critical_field(coherence_length: f64) -> f64 {
    let mu0 = 4.0 * PI * 1e-7_f64;
    FLUX_QUANTUM / (2.0 * PI * mu0 * coherence_length.powi(2).max(1e-300))
}

/// Thermodynamic critical field H_c = Φ₀ / (2π√2 μ₀ λ ξ) (A/m).
pub fn thermodynamic_critical_field(lambda: f64, coherence_length: f64) -> f64 {
    let mu0 = 4.0 * PI * 1e-7_f64;
    FLUX_QUANTUM / (2.0 * PI * 2.0_f64.sqrt() * mu0 * lambda * coherence_length.max(1e-300))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // 1. FLUX_QUANTUM has the accepted value.
    #[test]
    fn test_flux_quantum_value() {
        assert!((FLUX_QUANTUM - 2.067833848e-15).abs() < 1e-22);
    }

    // 2. is_type_ii: kappa < threshold → false.
    #[test]
    fn test_is_type_ii_false() {
        assert!(!is_type_ii(0.3)); // deep Type I
    }

    // 3. is_type_ii: kappa > threshold → true.
    #[test]
    fn test_is_type_ii_true() {
        assert!(is_type_ii(1.0)); // clear Type II
    }

    // 4. Boundary kappa ≈ 1/√2.
    #[test]
    fn test_is_type_ii_boundary() {
        let kappa = 1.0 / 2.0_f64.sqrt();
        assert!(!is_type_ii(kappa)); // exactly at threshold → not Type II
        assert!(is_type_ii(kappa + 1e-9));
    }

    // 5. SuperconductorProps::kappa = lambda / xi.
    #[test]
    fn test_props_kappa_computation() {
        let props = SuperconductorProps::new(9.25, 1.58e5, 2.44e5, 40e-9, 20e-9);
        assert!((props.kappa - 2.0).abs() < 1e-9);
    }

    // 6. Niobium is Type II.
    #[test]
    fn test_niobium_type_ii() {
        let nb = SuperconductorProps::niobium();
        assert_eq!(nb.superconductor_type(), SuperconductorType::TypeII);
    }

    // 7. Lead is Type I.
    #[test]
    fn test_lead_type_i() {
        let pb = SuperconductorProps::lead();
        assert_eq!(pb.superconductor_type(), SuperconductorType::TypeI);
    }

    // 8. YBCO is HighTc.
    #[test]
    fn test_ybco_high_tc() {
        let y = SuperconductorProps::ybco();
        assert_eq!(y.superconductor_type(), SuperconductorType::HighTc);
    }

    // 9. GL model initialises to |ψ| = 1 everywhere.
    #[test]
    fn test_gl_model_initial_order_parameter() {
        let gl = GinzburgLandauModel::new(10, 1e-9, 10e-9, 40e-9);
        assert!(gl.order_parameter.iter().all(|&v| (v - 1.0).abs() < EPS));
        assert_eq!(gl.dx, 1e-9);
        assert!((gl.kappa - 4.0).abs() < 1e-12);
    }

    // 10. solve_gl boundary conditions: |ψ(0)| = 0, |ψ(n-1)| = 1.
    #[test]
    fn test_gl_solve_boundary_conditions() {
        let mut gl = GinzburgLandauModel::new(20, 0.05, 0.1, 0.4);
        let psi = gl.solve_gl(200).unwrap();
        assert!(psi[0].abs() < EPS, "Surface must be 0, got {}", psi[0]);
        assert!(
            (psi[19] - 1.0).abs() < EPS,
            "Bulk must be 1, got {}",
            psi[19]
        );
    }

    // 11. solve_gl output values are in [0, 1].
    #[test]
    fn test_gl_solve_range() {
        let mut gl = GinzburgLandauModel::new(30, 0.05, 0.1, 0.4);
        let psi = gl.solve_gl(100).unwrap();
        for &v in &psi {
            assert!(
                (0.0..=1.0).contains(&v),
                "Order parameter out of range: {v}"
            );
        }
    }

    // 12. solve_gl is monotonically non-decreasing from surface to bulk.
    #[test]
    fn test_gl_solve_monotone() {
        let mut gl = GinzburgLandauModel::new(30, 0.05, 0.1, 0.4);
        let psi = gl.solve_gl(500).unwrap();
        for w in psi.windows(2) {
            assert!(w[1] >= w[0] - 1e-8, "Not monotone: {} > {}", w[0], w[1]);
        }
    }

    // 12a. The profile follows the analytic kink tanh(x / (sqrt(2) xi)) when the
    // grid spacing (in metres) resolves the coherence length (in metres).
    #[test]
    fn test_gl_solve_reproduces_tanh_kink_in_physical_units() {
        let xi = 10e-9;
        let dx = xi / 10.0;
        let mut gl = GinzburgLandauModel::new(401, dx, xi, 40e-9);
        let psi = gl.solve_gl(100).unwrap();
        let mut worst = 0.0_f64;
        for (i, &v) in psi.iter().enumerate() {
            let exact = (i as f64 * dx / (std::f64::consts::SQRT_2 * xi)).tanh();
            worst = worst.max((v - exact).abs());
        }
        assert!(worst < 1.0e-3, "max deviation from tanh profile = {worst}");
        // The first interior node must sit on the ramp, not on the bulk plateau.
        let expected_first = (dx / (std::f64::consts::SQRT_2 * xi)).tanh();
        assert!(
            (psi[1] - expected_first).abs() < 1.0e-3 && psi[1] < 0.1,
            "psi[1] = {} (kink value {expected_first})",
            psi[1]
        );
        // The stored profile is the returned one.
        assert_eq!(gl.order_parameter, psi);
    }

    // 12b. Only the ratio xi/dx matters: the same profile for any length unit.
    #[test]
    fn test_gl_solve_is_invariant_under_length_rescaling() {
        let mut reference = GinzburgLandauModel::new(101, 0.5e-9, 5.0e-9, 20.0e-9);
        let mut rescaled = GinzburgLandauModel::new(101, 0.05, 0.5, 2.0);
        let a = reference.solve_gl(100).unwrap();
        let b = rescaled.solve_gl(100).unwrap();
        for (x, y) in a.iter().zip(b.iter()) {
            assert!((x - y).abs() < 1.0e-12, "{x} vs {y}");
        }
    }

    // 12c. The converged profile satisfies the discrete GL equation to round-off.
    #[test]
    fn test_gl_solve_satisfies_discrete_equation() {
        let xi = 3.0e-9;
        let dx = 0.75e-9;
        let r = (xi / dx) * (xi / dx);
        let mut gl = GinzburgLandauModel::new(161, dx, xi, 12e-9);
        let psi = gl.solve_gl(50).unwrap();
        for w in psi.windows(3) {
            let res = r * (w[0] - 2.0 * w[1] + w[2]) + w[1] - w[1].powi(3);
            assert!(res.abs() < 1.0e-9, "residual {res} at psi = {}", w[1]);
        }
    }

    // 12d. Halving dx quarters the error: the scheme is second order.
    #[test]
    fn test_gl_solve_is_second_order_accurate() {
        let xi = 1.0_f64;
        let error = |points: usize, cells_per_xi: f64| {
            let dx = xi / cells_per_xi;
            let mut gl = GinzburgLandauModel::new(points, dx, xi, 4.0);
            let psi = gl.solve_gl(100).unwrap();
            psi.iter()
                .enumerate()
                .map(|(i, &v)| (v - (i as f64 * dx / (std::f64::consts::SQRT_2 * xi)).tanh()).abs())
                .fold(0.0_f64, f64::max)
        };
        // Both grids cover L = 40 xi.
        let coarse = error(201, 5.0);
        let fine = error(401, 10.0);
        let ratio = coarse / fine;
        assert!(
            (3.5..4.5).contains(&ratio),
            "error ratio {ratio} (coarse {coarse}, fine {fine}) should be ~4"
        );
    }

    // 12e. Unresolved coherence length (xi << dx): the profile jumps to the bulk value.
    #[test]
    fn test_gl_solve_unresolved_kink_is_a_step() {
        let mut gl = GinzburgLandauModel::new(12, 1.0, 0.05, 0.2);
        let psi = gl.solve_gl(100).unwrap();
        assert_eq!(psi[0], 0.0);
        for &v in &psi[2..] {
            assert!((v - 1.0).abs() < 1.0e-3, "bulk value expected, got {v}");
        }
    }

    // 12f. A converged profile is a fixed point: one Newton step confirms it.
    #[test]
    fn test_gl_solve_warm_start_is_a_fixed_point() {
        let mut gl = GinzburgLandauModel::new(81, 0.1, 1.0, 4.0);
        let first = gl.solve_gl(50).unwrap();
        let again = gl.solve_gl(1).unwrap();
        for (x, y) in first.iter().zip(again.iter()) {
            assert!((x - y).abs() < 1.0e-11, "{x} vs {y}");
        }
    }

    // 12g. Non-convergence within the iteration cap is an error and leaves the state alone.
    #[test]
    fn test_gl_solve_reports_non_convergence() {
        use crate::Error;
        let mut gl = GinzburgLandauModel::new(401, 0.1, 1.0, 4.0);
        let before = gl.order_parameter.clone();
        assert!(matches!(
            gl.solve_gl(1),
            Err(Error::NumericalDivergence { .. })
        ));
        assert!(matches!(
            gl.solve_gl(0),
            Err(Error::NumericalDivergence { .. })
        ));
        assert_eq!(gl.order_parameter, before);
        assert!(gl.solve_gl(100).is_ok());
    }

    // 12h. Invalid grid parameters are rejected.
    #[test]
    fn test_gl_solve_rejects_invalid_grid() {
        for (dx, xi) in [
            (0.0, 1.0),
            (-1.0, 1.0),
            (f64::NAN, 1.0),
            (1.0, 0.0),
            (1.0, -2.0),
            (1.0, f64::INFINITY),
            (1.0e-200, 1.0e200),
        ] {
            let mut gl = GinzburgLandauModel::new(10, dx, xi, 1.0);
            assert!(
                gl.solve_gl(10).is_err(),
                "dx={dx}, xi={xi} must be rejected"
            );
        }
    }

    // 12i. Degenerate grids.
    #[test]
    fn test_gl_solve_degenerate_grids() {
        let mut empty = GinzburgLandauModel::new(0, 1.0, 1.0, 2.0);
        assert!(empty.solve_gl(5).unwrap().is_empty());
        let mut single = GinzburgLandauModel::new(1, 1.0, 1.0, 2.0);
        assert_eq!(single.solve_gl(5).unwrap(), vec![1.0]);
        let mut pair = GinzburgLandauModel::new(2, 1.0, 1.0, 2.0);
        assert_eq!(pair.solve_gl(5).unwrap(), vec![0.0, 1.0]);
        let mut triple = GinzburgLandauModel::new(3, 1.0, 1.0, 2.0);
        let psi = triple.solve_gl(50).unwrap();
        assert_eq!((psi[0], psi[2]), (0.0, 1.0));
        assert!(psi[1] > 0.0 && psi[1] < 1.0);
    }

    // 13. BcsGap::from_tc gives Δ(0) ≈ 1.764 k_B T_c.
    #[test]
    fn test_bcs_gap_from_tc() {
        let tc = 9.25;
        let bcs = BcsGap::from_tc(tc, 1.0);
        let expected = 1.764 * K_B * tc;
        assert!((bcs.gap_energy - expected).abs() < 1e-28);
    }

    // 14. BCS gap at T = 0 equals gap_energy.
    #[test]
    fn test_bcs_gap_at_zero_temp() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        assert!((bcs.compute_gap(0.0) - bcs.gap_energy).abs() < EPS);
    }

    // 15. BCS gap vanishes at T = T_c.
    #[test]
    fn test_bcs_gap_at_tc() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        assert!(bcs.compute_gap(9.25).abs() < EPS);
    }

    // 16. BCS gap vanishes above T_c.
    #[test]
    fn test_bcs_gap_above_tc() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        assert!(bcs.compute_gap(100.0).abs() < EPS);
    }

    // 17. BCS gap is monotonically decreasing with temperature.
    #[test]
    fn test_bcs_gap_monotone() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        let temps: Vec<f64> = (0..=9).map(|i| i as f64).collect();
        let gaps: Vec<f64> = temps.iter().map(|&t| bcs.compute_gap(t)).collect();
        for w in gaps.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-15,
                "Gap not monotone: {} -> {}",
                w[0],
                w[1]
            );
        }
    }

    // 18. Normalised gap at T = 0 is 1.0.
    #[test]
    fn test_bcs_normalised_gap_at_zero() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        assert!((bcs.normalised_gap(0.0) - 1.0).abs() < EPS);
    }

    // 19. Normalised gap at T_c is 0.
    #[test]
    fn test_bcs_normalised_gap_at_tc() {
        let bcs = BcsGap::from_tc(9.25, 1.0);
        assert!(bcs.normalised_gap(9.25).abs() < EPS);
    }

    // 20. FluxVortex carries one flux quantum.
    #[test]
    fn test_flux_vortex_carries_one_quantum() {
        let v = FluxVortex::new([0.0, 0.0]);
        assert!((v.flux_quantum - FLUX_QUANTUM).abs() < 1e-25);
    }

    // 21. FluxVortex distance computation.
    #[test]
    fn test_flux_vortex_distance() {
        let v1 = FluxVortex::new([0.0, 0.0]);
        let v2 = FluxVortex::new([3.0, 4.0]);
        assert!((v1.distance_to(&v2) - 5.0).abs() < EPS);
    }

    // 22. VortexLattice starts empty.
    #[test]
    fn test_vortex_lattice_empty() {
        let lat = VortexLattice::new(40e-9);
        assert!(lat.vortices.is_empty());
    }

    // 23. add_vortex increases count.
    #[test]
    fn test_vortex_lattice_add() {
        let mut lat = VortexLattice::new(40e-9);
        lat.add_vortex(FluxVortex::new([0.0, 0.0]));
        lat.add_vortex(FluxVortex::new([1e-6, 0.0]));
        assert_eq!(lat.vortices.len(), 2);
    }

    // 24. Triangular lattice has nx * ny vortices.
    #[test]
    fn test_triangular_lattice_count() {
        let lat = VortexLattice::triangular_lattice(4, 3, 100e-9, 40e-9);
        assert_eq!(lat.vortices.len(), 12);
    }

    // 25. Inter-vortex force length matches vortex count.
    #[test]
    fn test_inter_vortex_force_length() {
        let lat = VortexLattice::triangular_lattice(3, 3, 200e-9, 40e-9);
        let forces = lat.compute_inter_vortex_force();
        assert_eq!(forces.len(), 9);
    }

    // 26. Net force on a symmetric pair has equal and opposite x-components.
    //
    // Vortex 0 is at origin, vortex 1 is at +x.  Repulsion pushes vortex 0
    // in the −x direction and vortex 1 in the +x direction.
    #[test]
    fn test_inter_vortex_force_newton_third() {
        let mut lat = VortexLattice::new(40e-9);
        lat.add_vortex(FluxVortex::new([0.0, 0.0]));
        lat.add_vortex(FluxVortex::new([200e-9, 0.0]));
        let f = lat.compute_inter_vortex_force();
        // Force on 0 should be −x (pushed away from vortex 1 in +x).
        assert!(
            f[0][0] < 0.0,
            "Vortex 0 should be pushed -x, got {}",
            f[0][0]
        );
        // Force on 1 should be +x (pushed away from vortex 0 in -x direction).
        assert!(
            f[1][0] > 0.0,
            "Vortex 1 should be pushed +x, got {}",
            f[1][0]
        );
        // Newton's 3rd law: equal and opposite.
        assert!(
            (f[0][0] + f[1][0]).abs() < 1e-3 * f[0][0].abs(),
            "Newton's 3rd law violated: {} vs {}",
            f[0][0],
            f[1][0]
        );
    }

    // 27. london_equation_factor returns 1.0 at zero field.
    #[test]
    fn test_london_factor_zero_field() {
        let f = london_equation_factor(40e-9, 0.0);
        assert!((f - 1.0).abs() < EPS);
    }

    // 28. london_equation_factor decays with field depth.
    #[test]
    fn test_london_factor_decay() {
        let lambda = 40e-9;
        let f1 = london_equation_factor(lambda, lambda);
        let f2 = london_equation_factor(lambda, 2.0 * lambda);
        assert!(f2 < f1, "Factor should decay: f1={f1}, f2={f2}");
        // f1 ≈ e^-1, f2 ≈ e^-2.
        assert!((f1 - (-1.0_f64).exp()).abs() < 1e-12);
    }

    // 29. london_equation_factor is zero for non-positive lambda.
    #[test]
    fn test_london_factor_zero_lambda() {
        assert!((london_equation_factor(0.0, 1.0)).abs() < EPS);
        assert!((london_equation_factor(-1.0, 1.0)).abs() < EPS);
    }

    // 30. upper_critical_field scales as 1/ξ².
    #[test]
    fn test_upper_critical_field_scaling() {
        let hc2_a = upper_critical_field(10e-9);
        let hc2_b = upper_critical_field(20e-9); // double ξ
        // H_c2 ∝ 1/ξ² → halving ξ quadruples H_c2.
        assert!(
            (hc2_a / hc2_b - 4.0).abs() < 1e-8,
            "H_c2(10nm)/H_c2(20nm) should be 4, got {}",
            hc2_a / hc2_b
        );
    }

    // 31. thermodynamic_critical_field is positive.
    #[test]
    fn test_thermodynamic_critical_field_positive() {
        let hc = thermodynamic_critical_field(40e-9, 20e-9);
        assert!(hc > 0.0);
    }

    // 32. SuperconductorType enum variant equality.
    #[test]
    fn test_superconductor_type_eq() {
        assert_eq!(SuperconductorType::TypeI, SuperconductorType::TypeI);
        assert_ne!(SuperconductorType::TypeI, SuperconductorType::TypeII);
    }

    // 33. VortexLattice flux density scales with vortex count.
    #[test]
    fn test_flux_density_scaling() {
        let mut lat1 = VortexLattice::new(40e-9);
        let mut lat2 = VortexLattice::new(40e-9);
        lat1.add_vortex(FluxVortex::new([0.0, 0.0]));
        lat2.add_vortex(FluxVortex::new([0.0, 0.0]));
        lat2.add_vortex(FluxVortex::new([1e-6, 0.0]));
        let area = 1e-12;
        assert!((lat2.flux_density(area) - 2.0 * lat1.flux_density(area)).abs() < 1e-25);
    }

    // 34. BCS gap at mid-temperature is between 0 and gap_energy.
    #[test]
    fn test_bcs_gap_mid_temperature() {
        let tc = 9.25;
        let bcs = BcsGap::from_tc(tc, 1.0);
        let gap_mid = bcs.compute_gap(tc / 2.0);
        assert!(
            gap_mid > 0.0 && gap_mid < bcs.gap_energy,
            "Mid-T gap should be in (0, Δ₀): {gap_mid}"
        );
    }
}
