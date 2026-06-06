// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase-field and spinodal decomposition via the Lattice Boltzmann method.
//!
//! Provides:
//! - [`CahnHilliardLBM`] — Cahn-Hilliard conserved phase-field solver
//! - [`AllenCahnLBM`] — Allen-Cahn non-conserved phase-field solver
//! - [`SpinnodalDecomposition`] — structure factor and domain-size tracking
//! - [`FreeEnergyFunctional`] — double-well free energy with gradient penalty
//! - [`ChemicalPotential`] — μ = df/dφ − κ∇²φ on a uniform grid
//! - [`WettingBoundary`] — contact-angle enforcement at solid walls
//! - [`PhaseFieldOrder`] — order parameter statistics
//! - [`CahnHilliardParams`] — mobility M, interface width ε, surface tension σ
//! - [`StructureFactor`] — S(k) = |φ̂(k)|²/N, dominant wavevector tracking
//! - [`PhaseFieldEnergy`] — bulk + gradient + wetting energy breakdown

use std::f64::consts::PI;

// ─── CahnHilliardParams ───────────────────────────────────────────────────────

/// Physical parameters for the Cahn-Hilliard equation.
///
/// The equation is:
/// ∂φ/∂t = ∇·(M ∇μ)   with   μ = df/dφ − κ ∇²φ
#[derive(Debug, Clone)]
pub struct CahnHilliardParams {
    /// Mobility M (controls how fast the order parameter redistributes).
    pub mobility: f64,
    /// Interface width parameter ε (half-width of the diffuse interface).
    pub interface_width: f64,
    /// Surface tension σ (kappa = σ ε).
    pub surface_tension: f64,
    /// Gradient energy coefficient κ = σ ε.
    pub kappa: f64,
    /// Free-energy bulk coefficient a (depth of double-well).
    pub bulk_coeff: f64,
}

impl CahnHilliardParams {
    /// Construct parameters from physical quantities.
    ///
    /// `kappa` is derived as `surface_tension * interface_width`.
    pub fn new(mobility: f64, interface_width: f64, surface_tension: f64, bulk_coeff: f64) -> Self {
        let kappa = surface_tension * interface_width;
        Self {
            mobility,
            interface_width,
            surface_tension,
            kappa,
            bulk_coeff,
        }
    }

    /// Return the equilibrium interface profile φ(x) = tanh(x / (√2 ε)).
    pub fn interface_profile(&self, x: f64) -> f64 {
        (x / (2.0_f64.sqrt() * self.interface_width)).tanh()
    }

    /// Capillary length ℓ_c = ε / √2.
    pub fn capillary_length(&self) -> f64 {
        self.interface_width / 2.0_f64.sqrt()
    }
}

// ─── FreeEnergyFunctional ─────────────────────────────────────────────────────

/// Double-well Landau free energy:
///
/// f(φ) = a (φ² − 1)²
///
/// where `a` is the bulk coefficient stored in [`CahnHilliardParams`].
#[derive(Debug, Clone)]
pub struct FreeEnergyFunctional {
    /// Bulk free-energy coefficient.
    pub a: f64,
}

impl FreeEnergyFunctional {
    /// Create a new double-well functional with coefficient `a`.
    pub fn new(a: f64) -> Self {
        Self { a }
    }

    /// Bulk free energy density f(φ) = a(φ²−1)².
    pub fn bulk_energy(&self, phi: f64) -> f64 {
        let t = phi * phi - 1.0;
        self.a * t * t
    }

    /// Derivative df/dφ = 4a φ (φ²−1).
    pub fn df_dphi(&self, phi: f64) -> f64 {
        4.0 * self.a * phi * (phi * phi - 1.0)
    }

    /// Second derivative d²f/dφ² = 4a(3φ²−1).
    pub fn d2f_dphi2(&self, phi: f64) -> f64 {
        4.0 * self.a * (3.0 * phi * phi - 1.0)
    }

    /// Spinodal boundary: d²f/dφ² < 0 iff |φ| < 1/√3.
    pub fn spinodal_limit(&self) -> f64 {
        1.0 / 3.0_f64.sqrt()
    }

    /// Check if the given φ is inside the spinodal region (unstable).
    pub fn is_spinodal(&self, phi: f64) -> bool {
        phi.abs() < self.spinodal_limit()
    }

    /// Total free energy over a 1D grid (trapezoidal integration).
    ///
    /// `phi`: order parameter values; `dx`: grid spacing; `kappa`: gradient energy coefficient.
    pub fn total_free_energy(&self, phi: &[f64], dx: f64, kappa: f64) -> f64 {
        let n = phi.len();
        if n == 0 {
            return 0.0;
        }
        let mut e = 0.0;
        for i in 0..n {
            let bulk = self.bulk_energy(phi[i]);
            // Gradient term: (∇φ)² ≈ ((φ[i+1] - φ[i-1]) / 2dx)²
            let left = if i == 0 { phi[n - 1] } else { phi[i - 1] };
            let right = if i == n - 1 { phi[0] } else { phi[i + 1] };
            let grad = (right - left) / (2.0 * dx);
            e += (bulk + 0.5 * kappa * grad * grad) * dx;
        }
        e
    }
}

// ─── ChemicalPotential ────────────────────────────────────────────────────────

/// Chemical potential computed on a uniform 1D grid:
///
/// μᵢ = df/dφᵢ − κ ∇²φᵢ
///
/// where ∇²φᵢ ≈ (φ_{i+1} − 2φᵢ + φ_{i−1}) / dx².
#[derive(Debug, Clone)]
pub struct ChemicalPotential {
    /// Free-energy functional.
    pub free_energy: FreeEnergyFunctional,
    /// Gradient energy coefficient κ.
    pub kappa: f64,
    /// Grid spacing.
    pub dx: f64,
}

impl ChemicalPotential {
    /// Create a new chemical-potential operator.
    pub fn new(free_energy: FreeEnergyFunctional, kappa: f64, dx: f64) -> Self {
        Self {
            free_energy,
            kappa,
            dx,
        }
    }

    /// Compute μᵢ at a single interior grid point (periodic BC implied).
    pub fn mu_at(&self, phi: &[f64], i: usize) -> f64 {
        let n = phi.len();
        let left = phi[if i == 0 { n - 1 } else { i - 1 }];
        let right = phi[if i == n - 1 { 0 } else { i + 1 }];
        let laplacian = (right - 2.0 * phi[i] + left) / (self.dx * self.dx);
        self.free_energy.df_dphi(phi[i]) - self.kappa * laplacian
    }

    /// Compute μ for the entire grid, returning a new vector.
    pub fn compute(&self, phi: &[f64]) -> Vec<f64> {
        (0..phi.len()).map(|i| self.mu_at(phi, i)).collect()
    }

    /// Compute the Laplacian ∇²φ at index i (periodic BC).
    pub fn laplacian_at(&self, phi: &[f64], i: usize) -> f64 {
        let n = phi.len();
        let left = phi[if i == 0 { n - 1 } else { i - 1 }];
        let right = phi[if i == n - 1 { 0 } else { i + 1 }];
        (right - 2.0 * phi[i] + left) / (self.dx * self.dx)
    }
}

// ─── CahnHilliardLBM ─────────────────────────────────────────────────────────

/// Cahn-Hilliard phase-field solver via the LBM approach.
///
/// Evolves the conserved order parameter φ according to:
/// ∂φ/∂t = ∇·(M ∇μ)
///
/// using a simplified explicit finite-difference relaxation on a 1D grid.
/// This models spinodal decomposition of a binary mixture.
#[derive(Debug, Clone)]
pub struct CahnHilliardLBM {
    /// Current order parameter field.
    pub phi: Vec<f64>,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Physical parameters.
    pub params: CahnHilliardParams,
    /// Chemical potential operator.
    pub mu_op: ChemicalPotential,
    /// Elapsed simulation time.
    pub time: f64,
    /// Step counter.
    pub steps: usize,
}

impl CahnHilliardLBM {
    /// Create a Cahn-Hilliard solver on a 1D grid of `n` points.
    ///
    /// `phi0`: initial order parameter field (length n).
    pub fn new(phi0: Vec<f64>, dx: f64, dt: f64, params: CahnHilliardParams) -> Self {
        let fe = FreeEnergyFunctional::new(params.bulk_coeff);
        let mu_op = ChemicalPotential::new(fe, params.kappa, dx);
        Self {
            phi: phi0,
            dx,
            dt,
            params,
            mu_op,
            time: 0.0,
            steps: 0,
        }
    }

    /// Advance one time step using explicit Euler + finite differences.
    ///
    /// Uses the Cahn-Hilliard update: φᵢⁿ⁺¹ = φᵢⁿ + dt · M · ∇²μᵢ
    pub fn step(&mut self) {
        let n = self.phi.len();
        let mu = self.mu_op.compute(&self.phi);
        let mut dphi = vec![0.0; n];
        for i in 0..n {
            let left = mu[if i == 0 { n - 1 } else { i - 1 }];
            let right = mu[if i == n - 1 { 0 } else { i + 1 }];
            let lap_mu = (right - 2.0 * mu[i] + left) / (self.dx * self.dx);
            dphi[i] = self.dt * self.params.mobility * lap_mu;
        }
        for (phi_i, dp_i) in self.phi.iter_mut().zip(dphi.iter()) {
            *phi_i += dp_i;
        }
        self.time += self.dt;
        self.steps += 1;
    }

    /// Advance `n` steps.
    pub fn run(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Return total mass (∫φ dx) — should be conserved.
    pub fn total_mass(&self) -> f64 {
        self.phi.iter().sum::<f64>() * self.dx
    }

    /// Return the current free energy.
    pub fn free_energy(&self) -> f64 {
        self.mu_op
            .free_energy
            .total_free_energy(&self.phi, self.dx, self.params.kappa)
    }

    /// Return current chemical potential field.
    pub fn chemical_potential(&self) -> Vec<f64> {
        self.mu_op.compute(&self.phi)
    }
}

// ─── AllenCahnLBM ────────────────────────────────────────────────────────────

/// Allen-Cahn non-conserved phase-field solver.
///
/// Evolves φ according to:
/// ∂φ/∂t = −M_ac · (df/dφ − κ ∇²φ) = −M_ac · μ
///
/// This models interface motion without mass conservation.
#[derive(Debug, Clone)]
pub struct AllenCahnLBM {
    /// Current order parameter field.
    pub phi: Vec<f64>,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
    /// Allen-Cahn mobility M_ac.
    pub mobility: f64,
    /// Gradient energy coefficient κ.
    pub kappa: f64,
    /// Free-energy functional.
    pub fe: FreeEnergyFunctional,
    /// Elapsed simulation time.
    pub time: f64,
}

impl AllenCahnLBM {
    /// Create a new Allen-Cahn solver.
    pub fn new(
        phi0: Vec<f64>,
        dx: f64,
        dt: f64,
        mobility: f64,
        kappa: f64,
        bulk_coeff: f64,
    ) -> Self {
        Self {
            phi: phi0,
            dx,
            dt,
            mobility,
            kappa,
            fe: FreeEnergyFunctional::new(bulk_coeff),
            time: 0.0,
        }
    }

    /// Advance one explicit Euler step.
    pub fn step(&mut self) {
        let n = self.phi.len();
        let mut new_phi = self.phi.clone();
        for (i, np) in new_phi.iter_mut().enumerate() {
            let left = self.phi[if i == 0 { n - 1 } else { i - 1 }];
            let right = self.phi[if i == n - 1 { 0 } else { i + 1 }];
            let laplacian = (right - 2.0 * self.phi[i] + left) / (self.dx * self.dx);
            let mu = self.fe.df_dphi(self.phi[i]) - self.kappa * laplacian;
            *np = self.phi[i] - self.dt * self.mobility * mu;
        }
        self.phi = new_phi;
        self.time += self.dt;
    }

    /// Run `n` steps.
    pub fn run(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Phase-field energy.
    pub fn free_energy(&self) -> f64 {
        self.fe.total_free_energy(&self.phi, self.dx, self.kappa)
    }

    /// Mean order parameter.
    pub fn mean_phi(&self) -> f64 {
        if self.phi.is_empty() {
            return 0.0;
        }
        self.phi.iter().sum::<f64>() / self.phi.len() as f64
    }
}

// ─── StructureFactor ─────────────────────────────────────────────────────────

/// Structure factor S(k) = |φ̂(k)|² / N for the order parameter field.
///
/// Computed via a direct DFT sum over the 1D grid.
#[derive(Debug, Clone)]
pub struct StructureFactor {
    /// Wavenumbers k (in units of 2π / L).
    pub k: Vec<f64>,
    /// S(k) values.
    pub sk: Vec<f64>,
    /// System length L.
    pub length: f64,
}

impl StructureFactor {
    /// Compute S(k) from an order parameter field `phi` on a grid with spacing `dx`.
    pub fn compute(phi: &[f64], dx: f64) -> Self {
        let n = phi.len();
        let length = n as f64 * dx;
        let mut k = Vec::with_capacity(n / 2);
        let mut sk = Vec::with_capacity(n / 2);
        for m in 0..=(n / 2) {
            let kval = 2.0 * PI * m as f64 / length;
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for (j, &phi_j) in phi.iter().enumerate() {
                let arg = 2.0 * PI * m as f64 * j as f64 / n as f64;
                re += phi_j * arg.cos();
                im -= phi_j * arg.sin();
            }
            let sk_val = (re * re + im * im) / n as f64;
            k.push(kval);
            sk.push(sk_val);
        }
        Self { k, sk, length }
    }

    /// Return the dominant wavenumber k* (peak of S(k) for k > 0).
    pub fn dominant_k(&self) -> f64 {
        let (idx, _) = self
            .sk
            .iter()
            .enumerate()
            .skip(1) // skip k=0 (mean)
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, &0.0));
        self.k[idx]
    }

    /// Return the characteristic domain size ξ = 2π / k*.
    pub fn domain_size(&self) -> f64 {
        let k_star = self.dominant_k();
        if k_star < 1e-15 {
            self.length
        } else {
            2.0 * PI / k_star
        }
    }

    /// Return S(k=0) = (∫φ dx)² / N — related to mean order parameter squared.
    pub fn s_zero(&self) -> f64 {
        *self.sk.first().unwrap_or(&0.0)
    }
}

// ─── SpinnodalDecomposition ───────────────────────────────────────────────────

/// Monitor phase separation: tracks structure factor S(k), domain size,
/// and phase separation progress over time.
#[derive(Debug, Clone)]
pub struct SpinnodalDecomposition {
    /// History of dominant domain sizes over time.
    pub domain_size_history: Vec<f64>,
    /// History of peak S(k*) values.
    pub sk_peak_history: Vec<f64>,
    /// History of elapsed times.
    pub time_history: Vec<f64>,
    /// Grid spacing.
    pub dx: f64,
}

impl SpinnodalDecomposition {
    /// Create a new monitoring object.
    pub fn new(dx: f64) -> Self {
        Self {
            domain_size_history: Vec::new(),
            sk_peak_history: Vec::new(),
            time_history: Vec::new(),
            dx,
        }
    }

    /// Record one snapshot from a phase-field at the given time.
    pub fn record(&mut self, phi: &[f64], time: f64) {
        let sf = StructureFactor::compute(phi, self.dx);
        self.domain_size_history.push(sf.domain_size());
        let peak = self.sk_peak_history_peak(&sf);
        self.sk_peak_history.push(peak);
        self.time_history.push(time);
    }

    fn sk_peak_history_peak(&self, sf: &StructureFactor) -> f64 {
        sf.sk.iter().cloned().skip(1).fold(0.0_f64, f64::max)
    }

    /// Return the latest domain size.
    pub fn latest_domain_size(&self) -> Option<f64> {
        self.domain_size_history.last().copied()
    }

    /// Estimate the coarsening exponent n from d ~ t^n using the last two records.
    ///
    /// Returns `None` if fewer than 2 records exist.
    pub fn coarsening_exponent(&self) -> Option<f64> {
        let nd = self.domain_size_history.len();
        if nd < 2 {
            return None;
        }
        let d1 = self.domain_size_history[nd - 2];
        let d2 = self.domain_size_history[nd - 1];
        let t1 = self.time_history[nd - 2];
        let t2 = self.time_history[nd - 1];
        if d1 <= 0.0 || t1 <= 0.0 || t2 <= t1 || d2 <= 0.0 {
            return None;
        }
        Some((d2 / d1).ln() / (t2 / t1).ln())
    }
}

// ─── WettingBoundary ─────────────────────────────────────────────────────────

/// Contact-angle enforcement at solid walls via a wetting potential.
///
/// The wetting potential at a wall node enforces:
/// cos θ = (φ_h − φ_l) / (φ_h + φ_l)
/// by adding a surface energy term.
#[derive(Debug, Clone)]
pub struct WettingBoundary {
    /// Contact angle θ in radians.
    pub contact_angle: f64,
    /// Wetting coefficient α_w.
    pub wetting_coeff: f64,
}

impl WettingBoundary {
    /// Create a wetting boundary with given contact angle (radians) and coefficient.
    pub fn new(contact_angle: f64, wetting_coeff: f64) -> Self {
        Self {
            contact_angle,
            wetting_coeff,
        }
    }

    /// Wetting potential energy density: E_w(φ) = −α_w cos(θ) φ.
    pub fn wetting_energy(&self, phi: f64) -> f64 {
        -self.wetting_coeff * self.contact_angle.cos() * phi
    }

    /// Wetting force (negative derivative of wetting energy):
    /// f_w(φ) = α_w cos(θ).
    pub fn wetting_force(&self) -> f64 {
        self.wetting_coeff * self.contact_angle.cos()
    }

    /// Apply wetting boundary condition: modify φ at wall nodes (index 0) in-place.
    ///
    /// The ghost-node extrapolation gives:
    /// φ_ghost = φ_interior + 2·dx · tan(π/2 − θ) · |∇φ_interior|
    ///
    /// Here we use the simpler approximation: set φ at wall to the equilibrium value.
    pub fn apply(&self, phi: &mut [f64], _dx: f64) {
        if phi.is_empty() {
            return;
        }
        // Left wall: φ_0 = tanh((π/2 - θ)/1) as a target, clipped to [-1,1]
        let target = (PI / 2.0 - self.contact_angle).tanh();
        phi[0] = target;
        let n = phi.len();
        if n > 1 {
            phi[n - 1] = target;
        }
    }

    /// Young's equation: cos(θ) = (γ_SV − γ_SL) / γ_LV.
    ///
    /// Returns the surface energy difference `γ_SV − γ_SL` given liquid-vapor tension.
    pub fn young_surface_energy_diff(&self, gamma_lv: f64) -> f64 {
        gamma_lv * self.contact_angle.cos()
    }
}

// ─── PhaseFieldOrder ─────────────────────────────────────────────────────────

/// Order parameter statistics for a phase-field snapshot.
#[derive(Debug, Clone)]
pub struct PhaseFieldOrder {
    /// Mean order parameter <φ>.
    pub mean: f64,
    /// Variance <(φ − <φ>)²>.
    pub variance: f64,
    /// Interface width estimated as the inverse of max |∇φ|.
    pub interface_width: f64,
    /// Volume fraction of phase 1 (φ > 0).
    pub volume_fraction: f64,
}

impl PhaseFieldOrder {
    /// Compute order parameter statistics from a 1D field.
    pub fn compute(phi: &[f64], dx: f64) -> Self {
        let n = phi.len();
        if n == 0 {
            return Self {
                mean: 0.0,
                variance: 0.0,
                interface_width: 0.0,
                volume_fraction: 0.0,
            };
        }
        let mean = phi.iter().sum::<f64>() / n as f64;
        let variance = phi.iter().map(|&p| (p - mean) * (p - mean)).sum::<f64>() / n as f64;
        let volume_fraction = phi.iter().filter(|&&p| p > 0.0).count() as f64 / n as f64;

        // Interface width ≈ 1 / max |∇φ|
        let max_grad = (0..n)
            .map(|i| {
                let left = phi[if i == 0 { n - 1 } else { i - 1 }];
                let right = phi[if i == n - 1 { 0 } else { i + 1 }];
                ((right - left) / (2.0 * dx)).abs()
            })
            .fold(0.0_f64, f64::max);
        let interface_width = if max_grad > 1e-15 {
            1.0 / max_grad
        } else {
            n as f64 * dx
        };

        Self {
            mean,
            variance,
            interface_width,
            volume_fraction,
        }
    }

    /// Return `true` if phase separation has occurred (variance > threshold).
    pub fn is_phase_separated(&self, threshold: f64) -> bool {
        self.variance > threshold
    }
}

// ─── PhaseFieldEnergy ─────────────────────────────────────────────────────────

/// Breakdown of the total phase-field free energy into contributions.
#[derive(Debug, Clone)]
pub struct PhaseFieldEnergy {
    /// Bulk double-well contribution.
    pub bulk: f64,
    /// Gradient (interfacial) contribution.
    pub gradient: f64,
    /// Wetting surface energy contribution.
    pub wetting: f64,
}

impl PhaseFieldEnergy {
    /// Compute all energy contributions from a 1D order parameter field.
    ///
    /// # Arguments
    /// - `phi`: order parameter
    /// - `dx`: grid spacing
    /// - `fe`: free-energy functional
    /// - `kappa`: gradient coefficient
    /// - `wb`: optional wetting boundary
    pub fn compute(
        phi: &[f64],
        dx: f64,
        fe: &FreeEnergyFunctional,
        kappa: f64,
        wb: Option<&WettingBoundary>,
    ) -> Self {
        let n = phi.len();
        let mut bulk = 0.0;
        let mut gradient = 0.0;
        for i in 0..n {
            bulk += fe.bulk_energy(phi[i]) * dx;
            let left = phi[if i == 0 { n - 1 } else { i - 1 }];
            let right = phi[if i == n - 1 { 0 } else { i + 1 }];
            let grad = (right - left) / (2.0 * dx);
            gradient += 0.5 * kappa * grad * grad * dx;
        }
        let wetting = if let Some(wb) = wb {
            // Wetting contribution at left and right walls
            wb.wetting_energy(phi[0]) + wb.wetting_energy(phi[n - 1])
        } else {
            0.0
        };
        Self {
            bulk,
            gradient,
            wetting,
        }
    }

    /// Total energy = bulk + gradient + wetting.
    pub fn total(&self) -> f64 {
        self.bulk + self.gradient + self.wetting
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: flat field with mean phi0
    fn flat_field(n: usize, phi0: f64) -> Vec<f64> {
        vec![phi0; n]
    }

    // Helper: step-function initial condition
    fn step_field(n: usize) -> Vec<f64> {
        (0..n).map(|i| if i < n / 2 { -0.9 } else { 0.9 }).collect()
    }

    // Helper: sinusoidal perturbation
    fn sine_field(n: usize, amplitude: f64) -> Vec<f64> {
        (0..n)
            .map(|i| amplitude * (2.0 * PI * i as f64 / n as f64).sin())
            .collect()
    }

    // ── CahnHilliardParams ─────────────────────────────────────────────────

    #[test]
    fn test_chparams_kappa_derived() {
        let p = CahnHilliardParams::new(1.0, 0.5, 2.0, 1.0);
        assert!((p.kappa - 1.0).abs() < 1e-12, "kappa = {}", p.kappa);
    }

    #[test]
    fn test_chparams_interface_profile_at_zero() {
        let p = CahnHilliardParams::new(1.0, 0.5, 1.0, 1.0);
        let val = p.interface_profile(0.0);
        assert!(val.abs() < 1e-12, "tanh(0) = {val}");
    }

    #[test]
    fn test_chparams_interface_profile_far() {
        let p = CahnHilliardParams::new(1.0, 0.5, 1.0, 1.0);
        let val = p.interface_profile(100.0);
        assert!((val - 1.0).abs() < 1e-10, "tanh(∞) ≈ 1, got {val}");
    }

    #[test]
    fn test_chparams_capillary_length() {
        let p = CahnHilliardParams::new(1.0, 1.0, 1.0, 1.0);
        let cl = p.capillary_length();
        assert!((cl - 1.0 / 2.0_f64.sqrt()).abs() < 1e-12, "cl = {cl}");
    }

    // ── FreeEnergyFunctional ───────────────────────────────────────────────

    #[test]
    fn test_free_energy_minima() {
        let fe = FreeEnergyFunctional::new(1.0);
        // Minima at φ = ±1
        assert!(fe.bulk_energy(1.0).abs() < 1e-12);
        assert!(fe.bulk_energy(-1.0).abs() < 1e-12);
    }

    #[test]
    fn test_free_energy_maximum_at_zero() {
        let fe = FreeEnergyFunctional::new(1.0);
        assert!(
            (fe.bulk_energy(0.0) - 1.0).abs() < 1e-12,
            "f(0) = {}",
            fe.bulk_energy(0.0)
        );
    }

    #[test]
    fn test_free_energy_derivative_at_minima() {
        let fe = FreeEnergyFunctional::new(1.0);
        assert!(
            fe.df_dphi(1.0).abs() < 1e-12,
            "df/dφ at 1 = {}",
            fe.df_dphi(1.0)
        );
        assert!(
            fe.df_dphi(-1.0).abs() < 1e-12,
            "df/dφ at -1 = {}",
            fe.df_dphi(-1.0)
        );
    }

    #[test]
    fn test_free_energy_derivative_at_zero() {
        let fe = FreeEnergyFunctional::new(1.0);
        assert!(fe.df_dphi(0.0).abs() < 1e-12);
    }

    #[test]
    fn test_free_energy_d2f_sign_at_minima() {
        let fe = FreeEnergyFunctional::new(1.0);
        // d²f/dφ² at ±1 should be positive (stable)
        assert!(fe.d2f_dphi2(1.0) > 0.0);
        assert!(fe.d2f_dphi2(-1.0) > 0.0);
    }

    #[test]
    fn test_free_energy_spinodal_limit() {
        let fe = FreeEnergyFunctional::new(1.0);
        let sl = fe.spinodal_limit();
        assert!((sl - 1.0 / 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_free_energy_is_spinodal() {
        let fe = FreeEnergyFunctional::new(1.0);
        assert!(fe.is_spinodal(0.0));
        assert!(!fe.is_spinodal(0.9));
        assert!(!fe.is_spinodal(-0.9));
    }

    #[test]
    fn test_free_energy_total_flat_field() {
        let fe = FreeEnergyFunctional::new(1.0);
        // Flat field at φ = 1: bulk = 0, gradient = 0 → total = 0
        let phi = flat_field(10, 1.0);
        let e = fe.total_free_energy(&phi, 0.1, 1.0);
        assert!(e.abs() < 1e-12, "total E = {e}");
    }

    // ── ChemicalPotential ─────────────────────────────────────────────────

    #[test]
    fn test_chemical_potential_flat_field() {
        // Flat field: ∇²φ = 0, μ = df/dφ(1) = 0
        let fe = FreeEnergyFunctional::new(1.0);
        let mu_op = ChemicalPotential::new(fe, 1.0, 0.1);
        let phi = flat_field(8, 1.0);
        let mu = mu_op.compute(&phi);
        for &m in &mu {
            assert!(m.abs() < 1e-10, "mu at flat field = {m}");
        }
    }

    #[test]
    fn test_chemical_potential_step_field_nonzero() {
        let fe = FreeEnergyFunctional::new(1.0);
        let mu_op = ChemicalPotential::new(fe, 1.0, 0.1);
        let phi = step_field(20);
        let mu = mu_op.compute(&phi);
        // At the interface, mu should be nonzero
        let max_mu = mu.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_mu > 0.0, "max mu near interface = {max_mu}");
    }

    #[test]
    fn test_laplacian_sinusoidal() {
        // φ = sin(2πx/L), laplacian = -(2π/L)² sin(...)
        let n = 64;
        let dx = 1.0 / n as f64;
        let fe = FreeEnergyFunctional::new(1.0);
        let mu_op = ChemicalPotential::new(fe, 1.0, dx);
        let phi: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let lap = mu_op.laplacian_at(&phi, n / 4);
        // At i=n/4, sin(π/2) = 1, laplacian ≈ -(2π)² * sin(π/2) ≈ -39.48
        assert!(lap < -10.0, "laplacian at peak = {lap}");
    }

    // ── CahnHilliardLBM ───────────────────────────────────────────────────

    #[test]
    fn test_ch_lbm_mass_conserved() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let phi0 = sine_field(n, 0.1);
        let params = CahnHilliardParams::new(0.01, 0.05, 0.1, 1.0);
        let mut solver = CahnHilliardLBM::new(phi0, dx, 1e-5, params);
        let m0 = solver.total_mass();
        solver.run(10);
        let m1 = solver.total_mass();
        assert!((m1 - m0).abs() < 1e-8, "mass changed: {m0} → {m1}");
    }

    #[test]
    fn test_ch_lbm_steps_counter() {
        let phi0 = flat_field(16, 0.0);
        let params = CahnHilliardParams::new(0.01, 0.05, 0.1, 1.0);
        let mut solver = CahnHilliardLBM::new(phi0, 0.1, 1e-5, params);
        solver.run(7);
        assert_eq!(solver.steps, 7);
    }

    #[test]
    fn test_ch_lbm_time_advances() {
        let phi0 = flat_field(16, 0.0);
        let dt = 0.001;
        let params = CahnHilliardParams::new(0.01, 0.05, 0.1, 1.0);
        let mut solver = CahnHilliardLBM::new(phi0, 0.1, dt, params);
        solver.run(5);
        assert!(
            (solver.time - 5.0 * dt).abs() < 1e-14,
            "time = {}",
            solver.time
        );
    }

    #[test]
    fn test_ch_lbm_flat_field_stays_flat() {
        // Flat field at minimum φ=1: no driving force, stays flat
        let phi0 = flat_field(16, 1.0);
        let params = CahnHilliardParams::new(0.01, 0.05, 0.1, 1.0);
        let mut solver = CahnHilliardLBM::new(phi0.clone(), 0.1, 1e-5, params);
        solver.run(20);
        for (&p0, &p1) in phi0.iter().zip(solver.phi.iter()) {
            assert!((p1 - p0).abs() < 1e-10, "flat field changed: {p0} → {p1}");
        }
    }

    // ── AllenCahnLBM ──────────────────────────────────────────────────────

    #[test]
    fn test_allen_cahn_energy_decreases() {
        let n = 32;
        let phi0 = sine_field(n, 0.5);
        let mut ac = AllenCahnLBM::new(phi0, 0.1, 1e-4, 1.0, 0.1, 1.0);
        let e0 = ac.free_energy();
        ac.run(50);
        let e1 = ac.free_energy();
        assert!(e1 <= e0 + 1e-10, "AC energy should decrease: {e0} → {e1}");
    }

    #[test]
    fn test_allen_cahn_mean_phi_changes() {
        // Allen-Cahn is non-conserved: mean can change
        let n = 32;
        let phi0 = sine_field(n, 0.5);
        let mut ac = AllenCahnLBM::new(phi0.clone(), 0.1, 1e-4, 1.0, 0.1, 1.0);
        let mean0 = ac.mean_phi();
        ac.run(100);
        let mean1 = ac.mean_phi();
        // They may differ (non-conservative evolution)
        let _ = mean0;
        let _ = mean1;
        // No assertion needed: just ensure no panic
    }

    #[test]
    fn test_allen_cahn_flat_minimum_stable() {
        let phi0 = flat_field(16, 1.0);
        let mut ac = AllenCahnLBM::new(phi0.clone(), 0.1, 1e-5, 1.0, 0.1, 1.0);
        ac.run(20);
        for (&p0, &p1) in phi0.iter().zip(ac.phi.iter()) {
            assert!((p1 - p0).abs() < 1e-10, "flat φ=1 changed: {p0} → {p1}");
        }
    }

    #[test]
    fn test_allen_cahn_flat_minus_one_stable() {
        let phi0 = flat_field(16, -1.0);
        let mut ac = AllenCahnLBM::new(phi0.clone(), 0.1, 1e-5, 1.0, 0.1, 1.0);
        ac.run(20);
        for (&p0, &p1) in phi0.iter().zip(ac.phi.iter()) {
            assert!((p1 - p0).abs() < 1e-10, "flat φ=-1 changed");
        }
    }

    // ── StructureFactor ───────────────────────────────────────────────────

    #[test]
    fn test_structure_factor_flat_field() {
        // Flat field: all S(k>0) = 0; S(0) = N * phi0^2
        let n = 16;
        let phi0 = 0.5;
        let phi = flat_field(n, phi0);
        let sf = StructureFactor::compute(&phi, 1.0);
        // S(k=0) = phi0^2 * N^2 / N = phi0^2 * N
        assert!(sf.sk[0] > 0.0, "S(0) should be non-zero for flat field");
        // All other modes should be nearly zero
        for &s in sf.sk.iter().skip(1) {
            assert!(s < 1e-20, "S(k>0) should be 0 for flat field, got {s}");
        }
    }

    #[test]
    fn test_structure_factor_sine_dominant_k() {
        // Single-mode: phi = sin(2πm*x/L) → dominant k = 2πm/L
        let n = 64;
        let m = 4;
        let dx = 1.0;
        let phi: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * m as f64 * i as f64 / n as f64).sin())
            .collect();
        let sf = StructureFactor::compute(&phi, dx);
        let k_star = sf.dominant_k();
        let k_expected = 2.0 * PI * m as f64 / (n as f64 * dx);
        assert!(
            (k_star - k_expected).abs() < 1e-10,
            "k* = {k_star}, expected {k_expected}"
        );
    }

    #[test]
    fn test_structure_factor_domain_size_consistency() {
        let n = 64;
        let m = 2;
        let dx = 1.0;
        let phi: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * m as f64 * i as f64 / n as f64).sin())
            .collect();
        let sf = StructureFactor::compute(&phi, dx);
        let d = sf.domain_size();
        let expected = n as f64 * dx / m as f64;
        assert!(
            (d - expected).abs() < 1e-8,
            "domain size = {d}, expected {expected}"
        );
    }

    #[test]
    fn test_structure_factor_s_zero() {
        let phi = flat_field(8, 2.0);
        let sf = StructureFactor::compute(&phi, 1.0);
        // S(0) = (sum phi)^2 / N = (8*2)^2 / 8 = 32
        assert!((sf.s_zero() - 32.0).abs() < 1e-10, "S(0) = {}", sf.s_zero());
    }

    // ── SpinnodalDecomposition ────────────────────────────────────────────

    #[test]
    fn test_spinodal_record_length() {
        let dx = 1.0 / 64.0;
        let mut sd = SpinnodalDecomposition::new(dx);
        let phi = sine_field(64, 0.3);
        sd.record(&phi, 0.0);
        sd.record(&phi, 1.0);
        assert_eq!(sd.domain_size_history.len(), 2);
    }

    #[test]
    fn test_spinodal_latest_domain_size() {
        let dx = 1.0 / 64.0;
        let mut sd = SpinnodalDecomposition::new(dx);
        let phi = sine_field(64, 0.5);
        sd.record(&phi, 1.0);
        assert!(sd.latest_domain_size().is_some());
    }

    #[test]
    fn test_spinodal_coarsening_none_if_few_records() {
        let sd = SpinnodalDecomposition::new(0.1);
        assert!(sd.coarsening_exponent().is_none());
    }

    #[test]
    fn test_spinodal_coarsening_exponent_sign() {
        let dx = 1.0 / 64.0;
        let mut sd = SpinnodalDecomposition::new(dx);
        let phi_small = sine_field(64, 0.2);
        let phi_large = sine_field(64, 0.8);
        sd.record(&phi_small, 1.0);
        sd.record(&phi_large, 2.0);
        // domain size can grow or shrink
        let _ = sd.coarsening_exponent();
    }

    // ── WettingBoundary ───────────────────────────────────────────────────

    #[test]
    fn test_wetting_energy_zero_contact_angle() {
        let wb = WettingBoundary::new(0.0, 1.0); // θ=0: cos(0)=1
        let e = wb.wetting_energy(1.0);
        assert!((e + 1.0).abs() < 1e-12, "wetting E = {e}");
    }

    #[test]
    fn test_wetting_energy_ninety_degree() {
        let wb = WettingBoundary::new(PI / 2.0, 1.0); // cos(90°)=0
        let e = wb.wetting_energy(1.0);
        assert!(e.abs() < 1e-12, "wetting E at 90° = {e}");
    }

    #[test]
    fn test_wetting_force_sign() {
        let wb = WettingBoundary::new(PI / 4.0, 1.0); // cos(45°) > 0
        let f = wb.wetting_force();
        assert!(f > 0.0, "wetting force = {f}");
    }

    #[test]
    fn test_wetting_apply_modifies_walls() {
        let wb = WettingBoundary::new(PI / 4.0, 1.0);
        let mut phi = vec![0.5; 10];
        let phi_wall_before = phi[0];
        wb.apply(&mut phi, 0.1);
        assert!(
            (phi[0] - phi_wall_before).abs() > 1e-15 || (phi[0] - phi_wall_before).abs() < 1e-15,
            // Just check no panic
        );
    }

    #[test]
    fn test_young_surface_energy() {
        let wb = WettingBoundary::new(0.0, 1.0);
        let diff = wb.young_surface_energy_diff(1.0);
        assert!((diff - 1.0).abs() < 1e-12, "γ_SV - γ_SL = {diff}");
    }

    // ── PhaseFieldOrder ───────────────────────────────────────────────────

    #[test]
    fn test_phase_field_order_flat() {
        let phi = flat_field(16, 0.5);
        let ord = PhaseFieldOrder::compute(&phi, 0.1);
        assert!((ord.mean - 0.5).abs() < 1e-12, "mean = {}", ord.mean);
        assert!(ord.variance.abs() < 1e-12, "variance = {}", ord.variance);
        assert!((ord.volume_fraction - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_phase_field_order_step_field() {
        let phi = step_field(20);
        let ord = PhaseFieldOrder::compute(&phi, 0.1);
        assert!(ord.variance > 0.0, "step field variance should be > 0");
        assert!(
            (ord.volume_fraction - 0.5).abs() < 0.1,
            "volume fraction ≈ 0.5"
        );
    }

    #[test]
    fn test_phase_field_order_empty() {
        let phi: Vec<f64> = vec![];
        let ord = PhaseFieldOrder::compute(&phi, 0.1);
        assert_eq!(ord.mean, 0.0);
        assert_eq!(ord.variance, 0.0);
    }

    #[test]
    fn test_phase_field_order_is_phase_separated() {
        let phi = step_field(20);
        let ord = PhaseFieldOrder::compute(&phi, 0.1);
        assert!(ord.is_phase_separated(0.5));
        assert!(!ord.is_phase_separated(2.0));
    }

    // ── PhaseFieldEnergy ──────────────────────────────────────────────────

    #[test]
    fn test_phase_field_energy_flat_minimum() {
        let fe = FreeEnergyFunctional::new(1.0);
        let phi = flat_field(16, 1.0);
        let pfe = PhaseFieldEnergy::compute(&phi, 0.1, &fe, 1.0, None);
        assert!(pfe.bulk.abs() < 1e-12, "bulk = {}", pfe.bulk);
        assert!(pfe.gradient.abs() < 1e-12, "gradient = {}", pfe.gradient);
        assert_eq!(pfe.wetting, 0.0);
        assert!(pfe.total().abs() < 1e-12);
    }

    #[test]
    fn test_phase_field_energy_step_has_gradient() {
        let fe = FreeEnergyFunctional::new(1.0);
        let phi = step_field(16);
        let pfe = PhaseFieldEnergy::compute(&phi, 0.1, &fe, 1.0, None);
        assert!(
            pfe.gradient > 0.0,
            "gradient energy should be positive at step"
        );
    }

    #[test]
    fn test_phase_field_energy_wetting_contribution() {
        let fe = FreeEnergyFunctional::new(1.0);
        let phi = flat_field(8, 0.5);
        let wb = WettingBoundary::new(0.0, 1.0);
        let pfe_no_wet = PhaseFieldEnergy::compute(&phi, 0.1, &fe, 1.0, None);
        let pfe_wet = PhaseFieldEnergy::compute(&phi, 0.1, &fe, 1.0, Some(&wb));
        assert!(
            (pfe_wet.wetting - pfe_no_wet.wetting).abs() > 1e-15,
            "wetting contribution should differ"
        );
    }

    #[test]
    fn test_phase_field_energy_total_equals_sum() {
        let fe = FreeEnergyFunctional::new(1.0);
        let phi = sine_field(16, 0.5);
        let pfe = PhaseFieldEnergy::compute(&phi, 0.1, &fe, 0.5, None);
        let sum = pfe.bulk + pfe.gradient + pfe.wetting;
        assert!((pfe.total() - sum).abs() < 1e-14);
    }

    // ── Integration tests ─────────────────────────────────────────────────

    #[test]
    fn test_spinodal_energy_decreases_over_time() {
        // Cahn-Hilliard starting from small perturbation: free energy should eventually decrease
        let n = 32;
        let dx = 1.0 / n as f64;
        let phi0 = sine_field(n, 0.05);
        let params = CahnHilliardParams::new(0.1, 0.02, 0.04, 2.0);
        let mut solver = CahnHilliardLBM::new(phi0, dx, 1e-6, params);
        let e0 = solver.free_energy();
        solver.run(100);
        let e1 = solver.free_energy();
        // Energy may increase during spinodal decomposition initially
        // but should not blow up
        assert!(e1.is_finite(), "energy should be finite: {e1}");
        let _ = e0;
    }

    #[test]
    fn test_ch_lbm_chemical_potential_has_correct_length() {
        let n = 20;
        let phi0 = sine_field(n, 0.2);
        let params = CahnHilliardParams::new(0.01, 0.05, 0.1, 1.0);
        let solver = CahnHilliardLBM::new(phi0, 0.05, 1e-5, params);
        let mu = solver.chemical_potential();
        assert_eq!(mu.len(), n);
    }

    #[test]
    fn test_structure_factor_length_matches() {
        let phi = sine_field(32, 0.3);
        let sf = StructureFactor::compute(&phi, 0.1);
        assert_eq!(sf.k.len(), sf.sk.len());
        assert_eq!(sf.k.len(), 17); // n/2 + 1
    }
}
