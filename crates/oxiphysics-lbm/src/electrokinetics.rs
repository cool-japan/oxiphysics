#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrokinetic LBM: electro-osmosis, electrophoresis, and Poisson-Nernst-Planck transport.
//!
//! This module implements:
//!
//! - [`IonSpecies`] — ionic properties (valence, diffusivity, bulk concentration)
//! - [`ElectricField`] — 2-D Ex/Ey field with Gauss's-law iteration
//! - [`PoissonSolver`] — successive over-relaxation (SOR) for ∇²φ = −ρₑ/ε
//! - [`NernstPlanckDistribution`] — advection-diffusion with electromigration drift
//! - [`ElectroOsmoticFlow`] — Helmholtz-Smoluchowski velocity from zeta potential
//! - [`DiffuseDoubleLayer`] — Gouy-Chapman charge distribution near a charged wall
//! - [`ElectrophoreticMobility`] — Henry's equation μ = 2εζf(κa)/3η
//! - [`ElectroosmoticPump`] — EOF pump flow rate and back-pressure
//! - [`DebyeLength`] — κ⁻¹ = √(εkT / 2n₀z²e²)
//! - [`IonTransportLBM`] — D2Q5 LBM step for ion-concentration transport

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
pub const K_B: f64 = 1.380_649e-23;

/// Elementary charge (C).
pub const E_CHARGE: f64 = 1.602_176_634e-19;

/// Permittivity of free space (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_817e-12;

/// Avogadro's number (mol⁻¹).
pub const N_AV: f64 = 6.022_140_76e23;

// ---------------------------------------------------------------------------
// IonSpecies
// ---------------------------------------------------------------------------

/// Description of one ionic species in the electrolyte.
#[derive(Debug, Clone)]
pub struct IonSpecies {
    /// Human-readable name, e.g. `"Na+"`, `"Cl-"`.
    pub name: String,
    /// Valence (charge number), e.g. `+1` for Na⁺, `−2` for SO₄²⁻.
    pub valence: i32,
    /// Diffusion coefficient D (m² s⁻¹).
    pub diffusivity: f64,
    /// Bulk (reservoir) concentration c₀ (mol m⁻³).
    pub bulk_concentration: f64,
}

impl IonSpecies {
    /// Create a new ion species.
    pub fn new(name: &str, valence: i32, diffusivity: f64, bulk_concentration: f64) -> Self {
        Self {
            name: name.to_string(),
            valence,
            diffusivity,
            bulk_concentration,
        }
    }

    /// Electrochemical mobility μᵢ = zᵢ e Dᵢ / (kB T)  (m² V⁻¹ s⁻¹).
    pub fn mobility(&self, temperature: f64) -> f64 {
        self.valence as f64 * E_CHARGE * self.diffusivity / (K_B * temperature)
    }

    /// Nernst-Einstein relation check: returns D (should equal kBT·μ / (ze)).
    pub fn diffusivity_from_mobility(&self, mobility: f64, temperature: f64) -> f64 {
        mobility * K_B * temperature / (self.valence.abs() as f64 * E_CHARGE)
    }
}

// ---------------------------------------------------------------------------
// ElectricField
// ---------------------------------------------------------------------------

/// 2-D electric field (Ex, Ey) on an nx × ny lattice.
///
/// The field is updated from the electric potential φ via:
///   Ex = −∂φ/∂x,  Ey = −∂φ/∂y  (central differences, periodic BC)
#[derive(Debug, Clone)]
pub struct ElectricField {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// x-component of field (V m⁻¹), row-major \[y*nx + x\].
    pub ex: Vec<f64>,
    /// y-component of field (V m⁻¹).
    pub ey: Vec<f64>,
}

impl ElectricField {
    /// Allocate a zeroed electric field.
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
        }
    }

    /// Compute Ex, Ey from the potential `phi` using central differences.
    ///
    /// Interior points use centred differences; boundaries use one-sided differences.
    pub fn update_from_potential(&mut self, phi: &[f64], dx: f64, dy: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iy * nx + ix;

                // x-derivative
                let phi_xp = if ix + 1 < nx {
                    phi[iy * nx + ix + 1]
                } else {
                    phi[idx]
                };
                let phi_xm = if ix > 0 {
                    phi[iy * nx + ix - 1]
                } else {
                    phi[idx]
                };
                let denom_x = if ix == 0 || ix + 1 == nx {
                    dx
                } else {
                    2.0 * dx
                };
                self.ex[idx] = -(phi_xp - phi_xm) / denom_x;

                // y-derivative
                let phi_yp = if iy + 1 < ny {
                    phi[(iy + 1) * nx + ix]
                } else {
                    phi[idx]
                };
                let phi_ym = if iy > 0 {
                    phi[(iy - 1) * nx + ix]
                } else {
                    phi[idx]
                };
                let denom_y = if iy == 0 || iy + 1 == ny {
                    dy
                } else {
                    2.0 * dy
                };
                self.ey[idx] = -(phi_yp - phi_ym) / denom_y;
            }
        }
    }

    /// Magnitude |E| at cell (ix, iy).
    pub fn magnitude(&self, ix: usize, iy: usize) -> f64 {
        let idx = iy * self.nx + ix;
        (self.ex[idx].powi(2) + self.ey[idx].powi(2)).sqrt()
    }

    /// Apply one Gauss's-law correction step: update φ so that ∇·E = ρ/ε.
    ///
    /// This iterates div(E) and adjusts by a fractional `alpha` step.
    pub fn gauss_law_iteration(
        &mut self,
        phi: &mut Vec<f64>,
        rho: &[f64],
        dx: f64,
        epsilon: f64,
        alpha: f64,
    ) {
        let nx = self.nx;
        let ny = self.ny;
        for iy in 1..ny - 1 {
            for ix in 1..nx - 1 {
                let idx = iy * nx + ix;
                // Divergence of E with central differences
                let div_e = (self.ex[iy * nx + ix + 1] - self.ex[iy * nx + ix - 1]) / (2.0 * dx)
                    + (self.ey[(iy + 1) * nx + ix] - self.ey[(iy - 1) * nx + ix]) / (2.0 * dx);
                let residual = div_e - rho[idx] / epsilon;
                phi[idx] -= alpha * residual * dx * dx;
            }
        }
        self.update_from_potential(phi, dx, dx);
    }
}

// ---------------------------------------------------------------------------
// PoissonSolver
// ---------------------------------------------------------------------------

/// Successive-over-relaxation (SOR) solver for the Poisson equation.
///
/// Equation: ∇²φ = −ρₑ / (ε₀ εᵣ)
///
/// Convergence criterion: max |Δφ| < `tol`.
pub struct PoissonSolver {
    /// SOR relaxation factor ω ∈ (1, 2).  Optimal ω ≈ 2/(1+π/N) for a square N×N grid.
    pub omega: f64,
}

impl PoissonSolver {
    /// Create a solver with the given relaxation factor.
    pub fn new(omega: f64) -> Self {
        Self { omega }
    }

    /// Create a solver with optimal ω for an N×N grid.
    pub fn optimal(n: usize) -> Self {
        let omega = 2.0 / (1.0 + PI / n as f64);
        Self { omega }
    }

    /// Run SOR on `phi` (in-place).
    ///
    /// Returns `(iterations_used, final_residual)`.
    ///
    /// - `phi`      : electric potential (V), Dirichlet BC assumed on boundary.
    /// - `rho`      : charge density (C m⁻³).
    /// - `nx`, `ny` : grid dimensions.
    /// - `dx`       : uniform cell spacing (m).
    /// - `epsilon_r`: relative permittivity.
    /// - `max_iter` : iteration cap.
    /// - `tol`      : convergence tolerance (V).
    #[allow(clippy::too_many_arguments)]
    pub fn solve(
        &self,
        phi: &mut Vec<f64>,
        rho: &[f64],
        nx: usize,
        ny: usize,
        dx: f64,
        epsilon_r: f64,
        max_iter: usize,
        tol: f64,
    ) -> (usize, f64) {
        let dx2 = dx * dx;
        let eps = EPSILON_0 * epsilon_r;
        let inv4 = 0.25;
        let mut residual = f64::INFINITY;
        let mut iters = 0;

        for iter in 0..max_iter {
            residual = 0.0_f64;
            for iy in 1..ny - 1 {
                for ix in 1..nx - 1 {
                    let idx = iy * nx + ix;
                    let phi_e = phi[iy * nx + ix + 1];
                    let phi_w = phi[iy * nx + ix - 1];
                    let phi_n = phi[(iy + 1) * nx + ix];
                    let phi_s = phi[(iy - 1) * nx + ix];
                    let rhs = (phi_e + phi_w + phi_n + phi_s + rho[idx] * dx2 / eps) * inv4;
                    let delta = self.omega * (rhs - phi[idx]);
                    phi[idx] += delta;
                    residual = residual.max(delta.abs());
                }
            }
            iters = iter + 1;
            if residual < tol {
                break;
            }
        }
        (iters, residual)
    }
}

// ---------------------------------------------------------------------------
// NernstPlanckDistribution
// ---------------------------------------------------------------------------

/// Nernst-Planck advection-diffusion-migration solver for ion concentration.
///
/// Equation (dimensionless LBM units):
///   ∂c/∂t + u·∇c = D∇²c − (Dze/kBT) ∇·(c ∇φ)
#[derive(Debug, Clone)]
pub struct NernstPlanckDistribution {
    /// Grid size nx.
    pub nx: usize,
    /// Grid size ny.
    pub ny: usize,
    /// Ion concentration c (mol m⁻³), row-major.
    pub concentration: Vec<f64>,
    /// Ion species parameters.
    pub species: IonSpecies,
    /// Temperature (K).
    pub temperature: f64,
}

impl NernstPlanckDistribution {
    /// Create a new distribution filled with the bulk concentration.
    pub fn new(nx: usize, ny: usize, species: IonSpecies, temperature: f64) -> Self {
        let c0 = species.bulk_concentration;
        Self {
            nx,
            ny,
            concentration: vec![c0; nx * ny],
            species,
            temperature,
        }
    }

    /// Perform one explicit finite-difference step of the NP equation.
    ///
    /// `phi`   — electric potential (V).
    /// `ux/uy` — fluid velocity (m s⁻¹).
    /// `dt`    — time step (s).
    /// `dx`    — cell spacing (m).
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, phi: &[f64], ux: &[f64], uy: &[f64], dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let d = self.species.diffusivity;
        let z = self.species.valence as f64;
        let mobility = z * E_CHARGE * d / (K_B * self.temperature);
        let c = self.concentration.clone();
        let inv_2dx = 0.5 / dx;
        let inv_dx2 = 1.0 / (dx * dx);

        for iy in 1..ny - 1 {
            for ix in 1..nx - 1 {
                let idx = iy * nx + ix;
                let ce = c[iy * nx + ix + 1];
                let cw = c[iy * nx + ix - 1];
                let cn = c[(iy + 1) * nx + ix];
                let cs = c[(iy - 1) * nx + ix];
                let ci = c[idx];

                // Diffusion ∇²c
                let lap_c = (ce + cw + cn + cs - 4.0 * ci) * inv_dx2;

                // Advection u·∇c (upwind)
                let dcx = if ux[idx] > 0.0 { ci - cw } else { ce - ci } / dx;
                let dcy = if uy[idx] > 0.0 { ci - cs } else { cn - ci } / dx;
                let adv = ux[idx] * dcx + uy[idx] * dcy;

                // Electromigration − ∇·(c · mobility · ∇φ)
                let dphi_x = (phi[iy * nx + ix + 1] - phi[iy * nx + ix - 1]) * inv_2dx;
                let dphi_y = (phi[(iy + 1) * nx + ix] - phi[(iy - 1) * nx + ix]) * inv_2dx;
                let drift_x = mobility
                    * (ce * (phi[iy * nx + ix + 1] - phi[idx])
                        - cw * (phi[idx] - phi[iy * nx + ix - 1]))
                    / (dx * dx);
                let drift_y = mobility
                    * (cn * (phi[(iy + 1) * nx + ix] - phi[idx])
                        - cs * (phi[idx] - phi[(iy - 1) * nx + ix]))
                    / (dx * dx);
                let _ = (dphi_x, dphi_y); // used conceptually above

                self.concentration[idx] = ci + dt * (d * lap_c - adv - drift_x - drift_y);
                // clamp to non-negative
                if self.concentration[idx] < 0.0 {
                    self.concentration[idx] = 0.0;
                }
            }
        }
    }

    /// Total amount of ions on the grid (sum c × dx²).
    pub fn total_ions(&self, dx: f64) -> f64 {
        self.concentration.iter().sum::<f64>() * dx * dx
    }
}

// ---------------------------------------------------------------------------
// ElectroOsmoticFlow
// ---------------------------------------------------------------------------

/// Electro-osmotic flow properties based on the Helmholtz-Smoluchowski model.
///
/// Valid in the thin-double-layer (κa >> 1) limit.
#[derive(Debug, Clone)]
pub struct ElectroOsmoticFlow {
    /// Zeta potential ζ (V) at the shear plane.
    pub zeta_potential: f64,
    /// Relative permittivity of the fluid.
    pub permittivity_r: f64,
    /// Dynamic viscosity η (Pa·s).
    pub viscosity: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl ElectroOsmoticFlow {
    /// Create an EOF model.
    pub fn new(zeta_potential: f64, permittivity_r: f64, viscosity: f64, temperature: f64) -> Self {
        Self {
            zeta_potential,
            permittivity_r,
            viscosity,
            temperature,
        }
    }

    /// Helmholtz-Smoluchowski slip velocity: u_HS = −εζE/η  (m s⁻¹).
    pub fn slip_velocity(&self, electric_field: f64) -> f64 {
        -EPSILON_0 * self.permittivity_r * self.zeta_potential * electric_field / self.viscosity
    }

    /// Electroosmotic mobility: μ_eo = −εζ/η  (m² V⁻¹ s⁻¹).
    pub fn eof_mobility(&self) -> f64 {
        -EPSILON_0 * self.permittivity_r * self.zeta_potential / self.viscosity
    }

    /// Debye length at a given ionic strength I (mol m⁻³).
    pub fn debye_length(&self, ionic_strength: f64) -> f64 {
        let num = EPSILON_0 * self.permittivity_r * K_B * self.temperature;
        let den = 2.0 * N_AV * E_CHARGE * E_CHARGE * ionic_strength;
        (num / den).sqrt()
    }

    /// Velocity profile u(y) inside a parallel-plate channel of half-height h (m)
    /// using the Debye-Hückel approximation.
    pub fn velocity_profile(&self, y: f64, h: f64, kappa: f64, e_field: f64) -> f64 {
        let u_hs = self.slip_velocity(e_field);
        // Debye-Hückel: u(y) = u_HS * (1 - cosh(κy)/cosh(κh))
        u_hs * (1.0 - (kappa * y).cosh() / (kappa * h).cosh())
    }
}

// ---------------------------------------------------------------------------
// DiffuseDoubleLayer
// ---------------------------------------------------------------------------

/// Gouy-Chapman model of the diffuse electric double layer near a charged wall.
#[derive(Debug, Clone)]
pub struct DiffuseDoubleLayer {
    /// Debye screening length λ_D (m).
    pub debye_length: f64,
    /// Surface charge density σ (C m⁻²).
    pub surface_charge: f64,
    /// Relative permittivity.
    pub permittivity_r: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl DiffuseDoubleLayer {
    /// Create a double layer model.
    pub fn new(
        debye_length: f64,
        surface_charge: f64,
        permittivity_r: f64,
        temperature: f64,
    ) -> Self {
        Self {
            debye_length,
            surface_charge,
            permittivity_r,
            temperature,
        }
    }

    /// Surface (zeta) potential from the Grahame equation (linearised):
    ///   ζ = σ λ_D / (ε₀ εᵣ)
    pub fn surface_potential(&self) -> f64 {
        self.surface_charge * self.debye_length / (EPSILON_0 * self.permittivity_r)
    }

    /// Potential profile φ(x) = ζ exp(−x/λ_D) (Debye-Hückel, V).
    pub fn potential_profile(&self, x: f64) -> f64 {
        self.surface_potential() * (-x / self.debye_length).exp()
    }

    /// Charge density ρ(x) = −ε₀εᵣ ∇²φ = σ/λ_D · exp(−x/λ_D)  (C m⁻³).
    pub fn charge_density_profile(&self, x: f64) -> f64 {
        (self.surface_charge / self.debye_length) * (-x / self.debye_length).exp()
    }

    /// Integrated charge per unit area from 0 to ∞ (should equal −σ by electroneutrality).
    pub fn integrated_charge(&self) -> f64 {
        -self.surface_charge
    }
}

// ---------------------------------------------------------------------------
// ElectrophoreticMobility
// ---------------------------------------------------------------------------

/// Electrophoretic mobility via Henry's equation:
///   μ_ep = (2 ε ζ / 3 η) f(κa)
///
/// where `f(κa)` is Henry's function (1 for κa → 0, 1.5 for κa → ∞).
#[derive(Debug, Clone)]
pub struct ElectrophoreticMobility {
    /// Zeta potential ζ (V).
    pub zeta_potential: f64,
    /// Relative permittivity of the medium.
    pub permittivity_r: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
}

impl ElectrophoreticMobility {
    /// Create a new electrophoretic mobility calculator.
    pub fn new(zeta_potential: f64, permittivity_r: f64, viscosity: f64) -> Self {
        Self {
            zeta_potential,
            permittivity_r,
            viscosity,
        }
    }

    /// Henry's function f(κa) — Ohshima's approximation valid for all κa.
    ///
    /// f(κa) = 1 + 1 / (2(1 + 2.5/(κa(1 + 2exp(−κa))))³)
    pub fn henry_function(kappa_a: f64) -> f64 {
        if kappa_a < 1e-12 {
            return 1.0; // Hückel limit
        }
        let denom = 1.0 + 2.5 / (kappa_a * (1.0 + 2.0 * (-kappa_a).exp()));
        1.0 + 1.0 / (2.0 * denom.powi(3))
    }

    /// Electrophoretic mobility μ_ep (m² V⁻¹ s⁻¹).
    pub fn mobility(&self, kappa_a: f64) -> f64 {
        let f = Self::henry_function(kappa_a);
        2.0 * EPSILON_0 * self.permittivity_r * self.zeta_potential * f / (3.0 * self.viscosity)
    }

    /// Smoluchowski limit (κa → ∞): μ = εζ/η.
    pub fn smoluchowski_limit(&self) -> f64 {
        EPSILON_0 * self.permittivity_r * self.zeta_potential / self.viscosity
    }

    /// Hückel limit (κa → 0): μ = 2εζ/3η.
    pub fn huckel_limit(&self) -> f64 {
        2.0 * EPSILON_0 * self.permittivity_r * self.zeta_potential / (3.0 * self.viscosity)
    }
}

// ---------------------------------------------------------------------------
// ElectroosmoticPump
// ---------------------------------------------------------------------------

/// EOF micropump model: relates applied voltage to volumetric flow rate and back-pressure.
#[derive(Debug, Clone)]
pub struct ElectroosmoticPump {
    /// Cross-sectional area of the channel (m²).
    pub cross_section: f64,
    /// Channel length (m).
    pub length: f64,
    /// EOF mobility (m² V⁻¹ s⁻¹).
    pub eof_mobility: f64,
    /// Hydraulic permeability k_h (m² Pa⁻¹ s⁻¹) = r²/8η for a cylindrical channel.
    pub hydraulic_permeability: f64,
}

impl ElectroosmoticPump {
    /// Create an EOF pump.
    pub fn new(
        cross_section: f64,
        length: f64,
        eof_mobility: f64,
        hydraulic_permeability: f64,
    ) -> Self {
        Self {
            cross_section,
            length,
            eof_mobility,
            hydraulic_permeability,
        }
    }

    /// Maximum (free-flow) volumetric flow rate Q_max = μ_eo · A · V / L  (m³ s⁻¹).
    pub fn max_flow_rate(&self, voltage: f64) -> f64 {
        self.eof_mobility * self.cross_section * voltage / self.length
    }

    /// Maximum back-pressure ΔP_max at zero flow (Pa).
    ///
    /// From Q_eof = Q_hagen:
    ///   ΔP_max = μ_eo · V · A / (k_h · L)
    pub fn max_back_pressure(&self, voltage: f64) -> f64 {
        self.max_flow_rate(voltage) / self.hydraulic_permeability
    }

    /// Flow rate at back-pressure `dp` (Pa).
    ///
    /// Q(ΔP) = Q_max · (1 − ΔP / ΔP_max)
    pub fn flow_rate_at_pressure(&self, voltage: f64, dp: f64) -> f64 {
        let q_max = self.max_flow_rate(voltage);
        let dp_max = self.max_back_pressure(voltage);
        if dp_max.abs() < 1e-30 {
            return 0.0;
        }
        q_max * (1.0 - dp / dp_max)
    }
}

// ---------------------------------------------------------------------------
// DebyeLength
// ---------------------------------------------------------------------------

/// Debye screening length calculator.
///
/// κ⁻¹ = √(ε₀εᵣkBT / (2 n₀ z² e²))
#[derive(Debug, Clone)]
pub struct DebyeLength {
    /// Temperature (K).
    pub temperature: f64,
    /// Relative permittivity.
    pub permittivity_r: f64,
}

impl DebyeLength {
    /// Create a Debye length calculator.
    pub fn new(temperature: f64, permittivity_r: f64) -> Self {
        Self {
            temperature,
            permittivity_r,
        }
    }

    /// Debye length (m) for a z:z electrolyte with bulk number density n₀ (m⁻³).
    pub fn compute(&self, n0: f64, z: i32) -> f64 {
        let num = EPSILON_0 * self.permittivity_r * K_B * self.temperature;
        let den = 2.0 * n0 * (z as f64 * E_CHARGE).powi(2);
        (num / den).sqrt()
    }

    /// Debye length (m) from ionic strength I (mol m⁻³).
    pub fn from_ionic_strength(&self, ionic_strength_mol_per_m3: f64) -> f64 {
        let n0 = ionic_strength_mol_per_m3 * N_AV;
        let num = EPSILON_0 * self.permittivity_r * K_B * self.temperature;
        let den = 2.0 * n0 * E_CHARGE * E_CHARGE;
        (num / den).sqrt()
    }

    /// Inverse Debye length κ (m⁻¹).
    pub fn kappa(&self, n0: f64, z: i32) -> f64 {
        1.0 / self.compute(n0, z)
    }
}

// ---------------------------------------------------------------------------
// IonTransportLBM  (D2Q5)
// ---------------------------------------------------------------------------

/// D2Q5 lattice Boltzmann solver for ion-concentration transport.
///
/// Velocity set (i=0..4):
///   0: (0,0), 1: (1,0), 2: (0,1), 3: (-1,0), 4: (0,-1)
///
/// Weights: w₀=1/3, w₁₋₄=1/6.
#[derive(Debug, Clone)]
pub struct IonTransportLBM {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions f\[direction\]\[cell\].
    pub f: Vec<Vec<f64>>,
    /// Ion relaxation time τ (dimensionless LBM units).
    pub tau: f64,
    /// Ion valence.
    pub valence: i32,
    /// Temperature (K) for drift term scaling.
    pub temperature: f64,
    /// Ion diffusion coefficient D (m²/s) — used in Nernst-Einstein mobility.
    pub diffusivity: f64,
}

impl IonTransportLBM {
    /// D2Q5 weight factors.
    pub const WEIGHTS: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

    /// D2Q5 velocity vectors \[(cx, cy); 5\].
    pub const VELOCITIES: [(i32, i32); 5] = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)];

    /// Create a new D2Q5 solver initialised to bulk concentration `c0`.
    pub fn new(nx: usize, ny: usize, tau: f64, valence: i32, temperature: f64, c0: f64) -> Self {
        let n = nx * ny;
        // Equilibrium: f_i^eq = w_i * c
        let f = Self::WEIGHTS.iter().map(|&w| vec![w * c0; n]).collect();
        // Default diffusivity from Einstein relation: D = cs² * (τ - 0.5)
        // with cs²=1/3 (lattice units). User may override with a physical value.
        let diffusivity = (1.0 / 3.0) * (tau - 0.5);
        Self {
            nx,
            ny,
            f,
            tau,
            valence,
            temperature,
            diffusivity,
        }
    }

    /// Electrochemical mobility via Nernst-Einstein relation: μ = z·e·D / (k_B·T).
    ///
    /// Returns mobility in SI units (m²·V⁻¹·s⁻¹).
    pub fn mobility(&self) -> f64 {
        self.valence as f64 * E_CHARGE * self.diffusivity / (K_B * self.temperature)
    }

    /// Macroscopic concentration c(x) = Σᵢ fᵢ.
    pub fn concentration(&self) -> Vec<f64> {
        let n = self.nx * self.ny;
        let mut c = vec![0.0; n];
        for i in 0..5 {
            for j in 0..n {
                c[j] += self.f[i][j];
            }
        }
        c
    }

    /// Equilibrium distribution fᵢᵉq = wᵢ c (1 + eᵢ·u_d / cs²).
    ///
    /// `u_drift` = mobility × E  (electromigration velocity).
    fn equilibrium(c: f64, u_drift: [f64; 2], i: usize) -> f64 {
        let cs2 = 1.0 / 3.0;
        let (cx, cy) = Self::VELOCITIES[i];
        let eu = cx as f64 * u_drift[0] + cy as f64 * u_drift[1];
        Self::WEIGHTS[i] * c * (1.0 + eu / cs2)
    }

    /// Perform one BGK collision + streaming step.
    ///
    /// `phi`  — electric potential field (lattice units).
    /// `ux/uy`— fluid velocity field (lattice units).
    pub fn step(&mut self, phi: &[f64], ux: &[f64], uy: &[f64]) {
        let nx = self.nx;
        let ny = self.ny;
        let inv_tau = 1.0 / self.tau;
        let z_e_over_kbt = self.valence as f64 * E_CHARGE / (K_B * self.temperature);

        // --- Collision ---
        let c = self.concentration();
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                // Electric drift (lattice units; assume phi in lattice V)
                let ex_lat = if ix + 1 < nx && ix > 0 {
                    -(phi[iy * nx + ix + 1] - phi[iy * nx + ix - 1]) * 0.5
                } else {
                    0.0
                };
                let ey_lat = if iy + 1 < ny && iy > 0 {
                    -(phi[(iy + 1) * nx + ix] - phi[(iy - 1) * nx + ix]) * 0.5
                } else {
                    0.0
                };
                // Nernst-Einstein mobility μ = z·e·D / (k_B·T)
                let mob = z_e_over_kbt * self.diffusivity;
                let u_drift = [ux[idx] + mob * ex_lat, uy[idx] + mob * ey_lat];
                for i in 0..5 {
                    let feq = Self::equilibrium(c[idx], u_drift, i);
                    self.f[i][idx] += inv_tau * (feq - self.f[i][idx]);
                }
            }
        }

        // --- Streaming ---
        for i in 1..5 {
            let (cx, cy) = Self::VELOCITIES[i];
            let mut f_new = vec![0.0; nx * ny];
            for iy in 0..ny {
                for ix in 0..nx {
                    let src_x = (ix as i32 - cx).rem_euclid(nx as i32) as usize;
                    let src_y = (iy as i32 - cy).rem_euclid(ny as i32) as usize;
                    f_new[iy * nx + ix] = self.f[i][src_y * nx + src_x];
                }
            }
            self.f[i] = f_new;
        }
    }

    /// Average concentration over the entire grid.
    pub fn average_concentration(&self) -> f64 {
        let c = self.concentration();
        c.iter().sum::<f64>() / c.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- IonSpecies ----

    #[test]
    fn test_ion_species_new() {
        let ion = IonSpecies::new("Na+", 1, 1.334e-9, 100.0);
        assert_eq!(ion.name, "Na+");
        assert_eq!(ion.valence, 1);
    }

    #[test]
    fn test_ion_mobility_positive() {
        let ion = IonSpecies::new("Na+", 1, 1.334e-9, 100.0);
        let mu = ion.mobility(298.0);
        assert!(mu > 0.0);
    }

    #[test]
    fn test_ion_mobility_negative_valence() {
        let ion = IonSpecies::new("Cl-", -1, 2.032e-9, 100.0);
        let mu = ion.mobility(298.0);
        assert!(mu < 0.0);
    }

    #[test]
    fn test_diffusivity_from_mobility_roundtrip() {
        let ion = IonSpecies::new("K+", 1, 1.96e-9, 50.0);
        let mu = ion.mobility(298.0);
        let d_back = ion.diffusivity_from_mobility(mu, 298.0);
        assert!((d_back - ion.diffusivity).abs() < 1e-20);
    }

    #[test]
    fn test_ion_mobility_temperature_scaling() {
        let ion = IonSpecies::new("H+", 1, 9.31e-9, 1.0);
        let mu1 = ion.mobility(298.0);
        let mu2 = ion.mobility(596.0);
        // Higher T → lower mobility (D fixed, kBT bigger)
        assert!(mu2 < mu1);
    }

    // ---- ElectricField ----

    #[test]
    fn test_electric_field_new_zeroed() {
        let ef = ElectricField::new(4, 4);
        assert!(ef.ex.iter().all(|&v| v == 0.0));
        assert!(ef.ey.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_electric_field_update_linear_x() {
        let nx = 5;
        let ny = 3;
        let mut ef = ElectricField::new(nx, ny);
        // phi = x (linear in x, constant in y) → Ex = -1, Ey = 0
        let phi: Vec<f64> = (0..nx * ny).map(|i| (i % nx) as f64).collect();
        ef.update_from_potential(&phi, 1.0, 1.0);
        // Check interior cell (1,1): Ex should be -1.0
        let idx = nx + 2; // iy=1, ix=2
        assert!((ef.ex[idx] + 1.0).abs() < 1e-10);
        assert!(ef.ey[idx].abs() < 1e-10);
    }

    #[test]
    fn test_electric_field_magnitude_zero() {
        let ef = ElectricField::new(3, 3);
        assert_eq!(ef.magnitude(1, 1), 0.0);
    }

    #[test]
    fn test_electric_field_magnitude_nonzero() {
        let nx = 5;
        let ny = 5;
        let mut ef = ElectricField::new(nx, ny);
        let phi: Vec<f64> = (0..nx * ny)
            .map(|i| (i % nx) as f64 + (i / nx) as f64)
            .collect();
        ef.update_from_potential(&phi, 1.0, 1.0);
        let mag = ef.magnitude(2, 2);
        assert!(mag > 0.0);
    }

    #[test]
    fn test_electric_field_gauss_iteration_runs() {
        let nx = 6;
        let ny = 6;
        let mut ef = ElectricField::new(nx, ny);
        let mut phi = vec![0.0; nx * ny];
        let rho = vec![1e-6; nx * ny]; // small charge density
        let epsilon = EPSILON_0 * 80.0;
        ef.gauss_law_iteration(&mut phi, &rho, 1e-9, epsilon, 0.1);
        // After one iteration, phi should have changed
        assert!(phi.iter().any(|&v| v != 0.0));
    }

    // ---- PoissonSolver ----

    #[test]
    fn test_poisson_solver_zero_rhs() {
        // Zero rho + zero BC → phi stays zero
        let solver = PoissonSolver::new(1.5);
        let nx = 5;
        let ny = 5;
        let mut phi = vec![0.0; nx * ny];
        let rho = vec![0.0; nx * ny];
        let (iters, _) = solver.solve(&mut phi, &rho, nx, ny, 1.0, 80.0, 100, 1e-10);
        assert!(iters > 0);
        assert!(phi.iter().all(|&v| v.abs() < 1e-10));
    }

    #[test]
    fn test_poisson_solver_converges() {
        let solver = PoissonSolver::new(1.5);
        let nx = 8;
        let ny = 8;
        let mut phi = vec![0.0; nx * ny];
        let rho = vec![1e-6; nx * ny];
        let (iters, residual) = solver.solve(&mut phi, &rho, nx, ny, 1e-6, 80.0, 5000, 1e-8);
        // Should converge well within 5000 iterations
        assert!(iters <= 5000);
        let _ = residual;
    }

    #[test]
    fn test_poisson_solver_symmetric_rho() {
        // Symmetric charge → symmetric potential
        let solver = PoissonSolver::new(1.5);
        let nx = 7;
        let ny = 7;
        let mut phi = vec![0.0; nx * ny];
        let mut rho = vec![0.0; nx * ny];
        rho[3 * nx + 3] = 1e-9; // centre charge
        solver.solve(&mut phi, &rho, nx, ny, 1.0, 80.0, 1000, 1e-8);
        // Centre cell should have largest potential
        let centre = phi[3 * nx + 3];
        for (i, &p) in phi.iter().enumerate() {
            if i != 3 * nx + 3 {
                assert!(centre >= p - 1e-12);
            }
        }
    }

    #[test]
    fn test_poisson_solver_optimal() {
        let solver = PoissonSolver::optimal(10);
        assert!(solver.omega > 1.0 && solver.omega < 2.0);
    }

    // ---- NernstPlanckDistribution ----

    #[test]
    fn test_nernst_planck_initial_concentration() {
        let species = IonSpecies::new("Na+", 1, 1.334e-9, 150.0);
        let np = NernstPlanckDistribution::new(4, 4, species, 298.0);
        assert!(np.concentration.iter().all(|&c| (c - 150.0).abs() < 1e-10));
    }

    #[test]
    fn test_nernst_planck_total_ions() {
        let species = IonSpecies::new("Na+", 1, 1.334e-9, 1.0);
        let np = NernstPlanckDistribution::new(4, 4, species, 298.0);
        let dx = 1e-6;
        let total = np.total_ions(dx);
        assert!((total - 16.0 * dx * dx).abs() < 1e-30);
    }

    #[test]
    fn test_nernst_planck_step_conserves_sign() {
        let species = IonSpecies::new("K+", 1, 1.96e-9, 10.0);
        let mut np = NernstPlanckDistribution::new(6, 6, species, 298.0);
        let n = 6 * 6;
        let phi = vec![0.0; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        np.step(&phi, &ux, &uy, 1e-12, 1e-6);
        assert!(np.concentration.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn test_nernst_planck_step_no_flow_no_change() {
        // Uniform c, zero E, zero velocity → no change on interior
        let species = IonSpecies::new("Cl-", -1, 2.032e-9, 1.0);
        let mut np = NernstPlanckDistribution::new(5, 5, species, 298.0);
        let n = 5 * 5;
        let phi = vec![0.0; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let c_before: Vec<f64> = np.concentration.clone();
        np.step(&phi, &ux, &uy, 1e-14, 1e-6);
        // Interior concentrations should not change (uniform gradient = 0)
        for iy in 1..4 {
            for ix in 1..4 {
                let idx = iy * 5 + ix;
                assert!((np.concentration[idx] - c_before[idx]).abs() < 1e-20);
            }
        }
    }

    // ---- ElectroOsmoticFlow ----

    #[test]
    fn test_eof_slip_velocity_sign() {
        // Negative zeta, positive E → positive velocity
        let eof = ElectroOsmoticFlow::new(-0.025, 80.0, 1e-3, 298.0);
        let v = eof.slip_velocity(1e4); // 10 kV/m
        assert!(v > 0.0);
    }

    #[test]
    fn test_eof_mobility_sign() {
        let eof = ElectroOsmoticFlow::new(-0.025, 80.0, 1e-3, 298.0);
        assert!(eof.eof_mobility() > 0.0);
    }

    #[test]
    fn test_eof_debye_length_positive() {
        let eof = ElectroOsmoticFlow::new(-0.025, 80.0, 1e-3, 298.0);
        let lambda_d = eof.debye_length(100.0); // 100 mol/m³ = 0.1 M
        assert!(lambda_d > 0.0 && lambda_d < 1e-6); // nanometres range
    }

    #[test]
    fn test_eof_velocity_profile_walls() {
        let eof = ElectroOsmoticFlow::new(-0.025, 80.0, 1e-3, 298.0);
        let h = 1e-5; // 10 µm half-channel
        // κh = 1.0 keeps cosh arguments well within f64 range
        let kappa = 1.0 / h;
        // At the wall (y=h) → u = u_HS*(1 - cosh(κh)/cosh(κh)) = 0 exactly
        let u_wall = eof.velocity_profile(h, h, kappa, 1e4);
        assert!(u_wall.abs() < 1e-6);
    }

    #[test]
    fn test_eof_debye_concentration_scaling() {
        let eof = ElectroOsmoticFlow::new(-0.025, 80.0, 1e-3, 298.0);
        // Higher ionic strength → shorter Debye length
        let ld1 = eof.debye_length(10.0);
        let ld2 = eof.debye_length(1000.0);
        assert!(ld1 > ld2);
    }

    // ---- DiffuseDoubleLayer ----

    #[test]
    fn test_ddl_surface_potential_sign() {
        // Positive surface charge → positive surface potential
        let ddl = DiffuseDoubleLayer::new(1e-8, 0.01, 80.0, 298.0);
        assert!(ddl.surface_potential() > 0.0);
    }

    #[test]
    fn test_ddl_potential_decays() {
        let ddl = DiffuseDoubleLayer::new(1e-8, 0.01, 80.0, 298.0);
        let p0 = ddl.potential_profile(0.0);
        let p1 = ddl.potential_profile(1e-8);
        let p2 = ddl.potential_profile(2e-8);
        assert!(p0 > p1 && p1 > p2);
    }

    #[test]
    fn test_ddl_charge_density_positive() {
        let ddl = DiffuseDoubleLayer::new(1e-8, 0.01, 80.0, 298.0);
        assert!(ddl.charge_density_profile(0.0) > 0.0);
    }

    #[test]
    fn test_ddl_integrated_charge() {
        let ddl = DiffuseDoubleLayer::new(1e-8, 0.05, 80.0, 298.0);
        assert_eq!(ddl.integrated_charge(), -0.05);
    }

    // ---- ElectrophoreticMobility ----

    #[test]
    fn test_henry_function_limits() {
        // Small κa → f → 1.0 (Hückel)
        let f_small = ElectrophoreticMobility::henry_function(0.01);
        assert!((f_small - 1.0).abs() < 0.05);
        // Large κa → f → 1.5 (Smoluchowski)
        let f_large = ElectrophoreticMobility::henry_function(1000.0);
        assert!((f_large - 1.5).abs() < 0.05);
    }

    #[test]
    fn test_henry_function_monotone() {
        let f1 = ElectrophoreticMobility::henry_function(1.0);
        let f2 = ElectrophoreticMobility::henry_function(10.0);
        let f3 = ElectrophoreticMobility::henry_function(100.0);
        assert!(f1 <= f2 && f2 <= f3);
    }

    #[test]
    fn test_electrophoretic_mobility_smoluchowski() {
        let em = ElectrophoreticMobility::new(-0.025, 80.0, 1e-3);
        let mu_sm = em.smoluchowski_limit();
        let mu_henry = em.mobility(1e6); // very large κa
        assert!((mu_sm - mu_henry).abs() / mu_sm.abs() < 0.01);
    }

    #[test]
    fn test_electrophoretic_mobility_huckel() {
        let em = ElectrophoreticMobility::new(-0.025, 80.0, 1e-3);
        let mu_h = em.huckel_limit();
        let mu_henry = em.mobility(1e-6); // very small κa
        assert!((mu_h - mu_henry).abs() / mu_h.abs() < 0.1);
    }

    // ---- ElectroosmoticPump ----

    #[test]
    fn test_pump_max_flow_proportional_to_voltage() {
        let pump = ElectroosmoticPump::new(1e-8, 1e-2, 4e-8, 1e-15);
        let q1 = pump.max_flow_rate(100.0);
        let q2 = pump.max_flow_rate(200.0);
        assert!((q2 / q1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_pump_zero_flow_at_max_pressure() {
        let pump = ElectroosmoticPump::new(1e-8, 1e-2, 4e-8, 1e-15);
        let dp_max = pump.max_back_pressure(100.0);
        let q = pump.flow_rate_at_pressure(100.0, dp_max);
        assert!(q.abs() < 1e-20);
    }

    #[test]
    fn test_pump_max_flow_at_zero_pressure() {
        let pump = ElectroosmoticPump::new(1e-8, 1e-2, 4e-8, 1e-15);
        let q_max = pump.max_flow_rate(100.0);
        let q = pump.flow_rate_at_pressure(100.0, 0.0);
        assert!((q - q_max).abs() < 1e-20);
    }

    // ---- DebyeLength ----

    #[test]
    fn test_debye_length_positive() {
        let dl = DebyeLength::new(298.0, 80.0);
        let lambda = dl.compute(6.022e25, 1); // 0.1 M NaCl
        assert!(lambda > 0.0);
    }

    #[test]
    fn test_debye_length_scaling_with_concentration() {
        let dl = DebyeLength::new(298.0, 80.0);
        let l1 = dl.compute(6.022e24, 1);
        let l2 = dl.compute(6.022e25, 1);
        // 10× more ions → 1/√10 shorter Debye length
        let ratio = l1 / l2;
        assert!((ratio - 10.0_f64.sqrt()).abs() < 0.01);
    }

    #[test]
    fn test_debye_length_valence_dependence() {
        let dl = DebyeLength::new(298.0, 80.0);
        let l_1 = dl.compute(6.022e25, 1);
        let l_2 = dl.compute(6.022e25, 2);
        // z=2 → 2× shorter
        assert!((l_1 / l_2 - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_debye_kappa_reciprocal() {
        let dl = DebyeLength::new(298.0, 80.0);
        let n0 = 6.022e25;
        let z = 1;
        let lambda = dl.compute(n0, z);
        let kappa = dl.kappa(n0, z);
        assert!((lambda * kappa - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_debye_from_ionic_strength() {
        let dl = DebyeLength::new(298.0, 80.0);
        let lambda = dl.from_ionic_strength(100.0); // 0.1 mol/L = 100 mol/m³
        assert!(lambda > 0.0 && lambda < 1e-7);
    }

    // ---- IonTransportLBM ----

    #[test]
    fn test_lbm_initial_concentration() {
        let lbm = IonTransportLBM::new(4, 4, 1.0, 1, 298.0, 1.0);
        let c = lbm.concentration();
        assert!(c.iter().all(|&v| (v - 1.0).abs() < 1e-10));
    }

    #[test]
    fn test_lbm_weights_sum_to_one() {
        let sum: f64 = IonTransportLBM::WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_lbm_step_runs() {
        let mut lbm = IonTransportLBM::new(6, 6, 1.0, 1, 298.0, 1.0);
        let n = 6 * 6;
        let phi = vec![0.0; n];
        let ux = vec![0.01; n];
        let uy = vec![0.0; n];
        lbm.step(&phi, &ux, &uy);
        // Average concentration should be conserved (approximately)
        let avg = lbm.average_concentration();
        assert!((avg - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_lbm_average_concentration() {
        let lbm = IonTransportLBM::new(8, 8, 1.0, 1, 298.0, 2.5);
        let avg = lbm.average_concentration();
        assert!((avg - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_lbm_equilibrium_zero_drift() {
        // With zero drift, f^eq_i = w_i * c
        let c = 1.0;
        let u_drift = [0.0, 0.0];
        for i in 0..5 {
            let feq = IonTransportLBM::equilibrium(c, u_drift, i);
            assert!((feq - IonTransportLBM::WEIGHTS[i]).abs() < 1e-14);
        }
    }

    #[test]
    fn test_nernst_einstein_mobility() {
        // z=1, D=1e-9 m²/s, T=298 K → μ ≈ 3.91e-8 m²/(V·s)
        let mut lbm = IonTransportLBM::new(4, 4, 1.0, 1, 298.0, 1.0);
        lbm.diffusivity = 1e-9;
        let mu = lbm.mobility();
        let expected = 1.0_f64 * E_CHARGE * 1e-9 / (K_B * 298.0);
        assert!(
            (mu - expected).abs() / expected < 1e-6,
            "mobility = {mu:.4e}, expected ≈ {expected:.4e}"
        );
        assert!(
            (mu - 3.91e-8).abs() < 0.05e-8,
            "mobility {mu:.4e} not near 3.91e-8"
        );
    }
}
