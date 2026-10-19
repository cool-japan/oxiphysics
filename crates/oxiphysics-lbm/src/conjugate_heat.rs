// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conjugate heat transfer using the Lattice Boltzmann Method.
//!
//! This module implements coupled fluid-solid heat transfer using:
//! - Fluid domain: D2Q9 LBM for velocity + D2Q5 for temperature
//! - Solid domain: finite-difference heat diffusion
//! - Interface: temperature and flux continuity conditions
//!
//! The interface between a fluid node `F` and a solid node `S` lies half-way between the two
//! node centres. Flux continuity `k_f (T_F − T_i) = k_s (T_i − T_S)` fixes the interface
//! temperature `T_i = (k_f T_F + k_s T_S) / (k_f + k_s)`, which acts as a Dirichlet value for both
//! sides: an anti-bounce-back link for the fluid temperature populations and a ghost value
//! `2 T_i − T_S` in the finite-difference stencil of the solid node. The node temperatures
//! themselves are never overwritten, so the interface flux
//! `q = 2 k_f k_s (T_S − T_F) / ((k_f + k_s) Δx)` is non-zero whenever the two sides differ and
//! a steady state satisfies `k_f ∇T_f = k_s ∇T_s`. Fluid–solid links without a registered
//! interface are adiabatic; the fluid flow bounces back from every solid link.
//!
//! Key references:
//! - Wang, J., et al. (2007). A lattice Boltzmann algorithm for fluid-solid
//!   conjugate heat transfer. *Int. J. Therm. Sci.* 46, 228–234.
//! - Karani, H., & Huber, C. (2015). Lattice Boltzmann formulation for
//!   conjugate heat transfer in heterogeneous media. *Phys. Rev. E* 91, 023304.

use crate::error::{Error, Result};
use crate::lattice::{D2Q9_OPPOSITES, D2Q9_VELOCITIES, D2Q9_WEIGHTS};

/// Speed of sound squared (cs² = 1/3).
const CS2: f64 = 1.0 / 3.0;

/// D2Q5 weights for temperature distribution function.
pub const D2Q5_WEIGHTS: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

/// D2Q5 velocity vectors `[cx, cy]`.
///
/// Order: rest, E, N, W, S.
pub const D2Q5_VELOCITIES: [[i32; 2]; 5] = [
    [0, 0],  // 0: rest
    [1, 0],  // 1: E
    [0, 1],  // 2: N
    [-1, 0], // 3: W
    [0, -1], // 4: S
];

/// Opposite direction indices for D2Q5.
pub const D2Q5_OPPOSITES: [usize; 5] = [0, 3, 4, 1, 2];

// ---------------------------------------------------------------------------
// SolidRegion
// ---------------------------------------------------------------------------

/// Material properties of a solid region for conjugate heat transfer.
///
/// Stores the thermophysical properties required for the solid heat diffusion
/// equation: ρ c_p ∂T/∂t = ∇·(k ∇T).
#[derive(Debug, Clone)]
pub struct SolidRegion {
    /// Thermal conductivity k (W/m·K).
    pub thermal_conductivity: f64,
    /// Specific heat capacity c_p (J/kg·K).
    pub specific_heat: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Region label / identifier.
    pub label: String,
}

impl SolidRegion {
    /// Create a new solid region with given thermophysical properties.
    pub fn new(thermal_conductivity: f64, specific_heat: f64, density: f64, label: &str) -> Self {
        Self {
            thermal_conductivity,
            specific_heat,
            density,
            label: label.to_string(),
        }
    }

    /// Thermal diffusivity α = k / (ρ c_p) (m²/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Thermal effusivity e = sqrt(k ρ c_p) (W·s^0.5/m²·K).
    pub fn thermal_effusivity(&self) -> f64 {
        (self.thermal_conductivity * self.density * self.specific_heat).sqrt()
    }
}

// ---------------------------------------------------------------------------
// FluidSolidInterface
// ---------------------------------------------------------------------------

/// Interface node connecting a fluid LBM node to a solid diffusion node.
///
/// At the interface the temperature continuity condition T_f = T_s and the
/// heat-flux continuity condition k_f (∂T/∂n)_f = k_s (∂T/∂n)_s are enforced.
/// The interface is located half-way between the fluid node and the solid node, which must be
/// axis-adjacent lattice neighbours (the link direction is derived from the two node indices).
#[derive(Debug, Clone)]
pub struct FluidSolidInterface {
    /// Flat index of the fluid-side node.
    pub fluid_index: usize,
    /// Flat index of the solid-side node.
    pub solid_index: usize,
    /// Outward unit normal from fluid to solid `[nx, ny]` (informational; the solver derives the
    /// link direction from the node indices).
    pub normal: [f64; 2],
    /// Thermal conductivity of the fluid side (W/m·K).
    pub k_fluid: f64,
    /// Thermal conductivity of the solid side (W/m·K).
    pub k_solid: f64,
    /// Interface temperature (updated each step).
    pub interface_temperature: f64,
}

impl FluidSolidInterface {
    /// Construct a new interface node.
    pub fn new(
        fluid_index: usize,
        solid_index: usize,
        normal: [f64; 2],
        k_fluid: f64,
        k_solid: f64,
    ) -> Self {
        Self {
            fluid_index,
            solid_index,
            normal,
            k_fluid,
            k_solid,
            interface_temperature: 0.0,
        }
    }

    /// Compute the interface temperature from the node temperatures with the
    /// conductivity-weighted mean that balances the fluxes of an interface half-way between
    /// the nodes, `k_f (T_f − T_i) = k_s (T_i − T_s)`:
    ///
    /// ```text
    /// T_int = (k_f * T_f + k_s * T_s) / (k_f + k_s)
    /// ```
    pub fn compute_interface_temperature(&mut self, t_fluid: f64, t_solid: f64) {
        self.interface_temperature =
            (self.k_fluid * t_fluid + self.k_solid * t_solid) / (self.k_fluid + self.k_solid);
    }

    /// Heat flux continuity residual (should be zero at convergence).
    pub fn flux_residual(&self, grad_t_fluid: f64, grad_t_solid: f64) -> f64 {
        self.k_fluid * grad_t_fluid - self.k_solid * grad_t_solid
    }

    /// Heat flux through the interface (W/m²) for the node temperatures `t_fluid` and
    /// `t_solid` at node spacing `dx`, positive when heat flows from the solid into the fluid.
    ///
    /// Both half-cells carry the same flux,
    /// `k_f (T_i − T_f) / (Δx/2) = k_s (T_s − T_i) / (Δx/2)`, which gives
    ///
    /// ```text
    /// q = 2 k_f k_s (T_s − T_f) / ((k_f + k_s) Δx)
    /// ```
    pub fn heat_flux(&self, t_fluid: f64, t_solid: f64, dx: f64) -> f64 {
        2.0 * self.k_fluid * self.k_solid * (t_solid - t_fluid)
            / ((self.k_fluid + self.k_solid) * dx)
    }
}

// ---------------------------------------------------------------------------
// ThermalLbmD2Q5
// ---------------------------------------------------------------------------

/// D2Q5 temperature distribution function LBM solver.
///
/// The advection-diffusion equation for temperature is solved using a 5-speed
/// lattice with populations `h_i`.  The equilibrium is:
///
/// ```text
/// h_eq_i = w_i * T * (1 + (e_i · u) / cs²)
/// ```
///
/// and the BGK collision is:
///
/// ```text
/// h_post_i = h_i - (h_i - h_eq_i) / tau_t
/// ```
#[derive(Debug, Clone)]
pub struct ThermalLbmD2Q5 {
    /// Temperature populations `h[node][direction]`.
    pub h: Vec<[f64; 5]>,
    /// Macroscopic temperature at each node.
    pub temperature: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Thermal relaxation time τ_t.
    pub tau_t: f64,
}

impl ThermalLbmD2Q5 {
    /// Create a new D2Q5 temperature solver initialised to `t_init` at rest.
    pub fn new(nx: usize, ny: usize, tau_t: f64, t_init: f64) -> Self {
        let n = nx * ny;
        let h_eq: [f64; 5] = {
            let mut arr = [0.0; 5];
            for (i, w) in D2Q5_WEIGHTS.iter().enumerate() {
                arr[i] = w * t_init;
            }
            arr
        };
        Self {
            h: vec![h_eq; n],
            temperature: vec![t_init; n],
            nx,
            ny,
            tau_t,
        }
    }

    /// Flat node index for position `(x, y)`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Thermal diffusivity α = (τ_t − 0.5) / 3.
    pub fn alpha(&self) -> f64 {
        (self.tau_t - 0.5) / 3.0
    }

    /// D2Q5 equilibrium distribution for direction `i`.
    pub fn equilibrium(&self, temp: f64, u: [f64; 2], i: usize) -> f64 {
        let w = D2Q5_WEIGHTS[i];
        let c = D2Q5_VELOCITIES[i];
        let e_dot_u = c[0] as f64 * u[0] + c[1] as f64 * u[1];
        w * temp * (1.0 + e_dot_u / CS2)
    }

    /// BGK collision on all nodes given fluid velocities.
    pub fn collide(&mut self, velocities: &[[f64; 2]]) {
        let n = self.nx * self.ny;
        for (k, vel) in velocities.iter().enumerate().take(n) {
            let t = self.temperature[k];
            let u = *vel;
            let tau = self.tau_t;
            // Precompute equilibrium values to avoid conflicting borrows
            let heq: [f64; 5] = std::array::from_fn(|i| self.equilibrium(t, u, i));
            for (h_ki, heq_i) in self.h[k].iter_mut().zip(heq.iter()) {
                *h_ki -= (*h_ki - heq_i) / tau;
            }
        }
    }

    /// Streaming step with periodic boundaries.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut h_new = vec![[0.0_f64; 5]; n];

        for y in 0..ny {
            for x in 0..nx {
                let k_src = self.idx(x, y);
                for i in 0..5 {
                    let c = D2Q5_VELOCITIES[i];
                    let xd = ((x as isize + c[0] as isize).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + c[1] as isize).rem_euclid(ny as isize)) as usize;
                    let k_dst = self.idx(xd, yd);
                    h_new[k_dst][i] = self.h[k_src][i];
                }
            }
        }
        self.h = h_new;
    }

    /// Compute macroscopic temperature from zeroth moment: T = Σ_i h_i.
    pub fn compute_temperature(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            self.temperature[k] = self.h[k].iter().sum();
        }
    }

    /// Apply Dirichlet temperature boundary at node `idx`.
    pub fn apply_dirichlet(&mut self, idx: usize, t_wall: f64, u: [f64; 2]) {
        for i in 0..5 {
            self.h[idx][i] = self.equilibrium(t_wall, u, i);
        }
        self.temperature[idx] = t_wall;
    }

    /// Apply adiabatic (zero-flux) boundary via anti-bounce-back at node `idx`.
    pub fn apply_adiabatic(&mut self, idx: usize, opposite_idx: usize) {
        // Zero-flux: copy from neighbour so no gradient forms.
        self.h[idx] = self.h[opposite_idx];
        self.temperature[idx] = self.temperature[opposite_idx];
    }
}

// ---------------------------------------------------------------------------
// DoubleDistribution
// ---------------------------------------------------------------------------

/// Dual-population LBM state combining D2Q9 fluid and D2Q5 temperature.
///
/// The velocity field is computed from the standard D2Q9 distributions `f_i`,
/// while the temperature field uses the separate D2Q5 distributions `h_i`.
#[derive(Debug, Clone)]
pub struct DoubleDistribution {
    /// D2Q9 fluid distributions `f[node][9]`.
    pub f: Vec<[f64; 9]>,
    /// D2Q5 temperature distributions `h[node][5]`.
    pub h: Vec<[f64; 5]>,
    /// Macroscopic density at each node.
    pub rho: Vec<f64>,
    /// Macroscopic velocity at each node `[ux, uy]`.
    pub u: Vec<[f64; 2]>,
    /// Macroscopic temperature at each node.
    pub temperature: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Fluid relaxation time τ.
    pub tau: f64,
    /// Thermal relaxation time τ_t.
    pub tau_t: f64,
}

impl DoubleDistribution {
    /// Create a new double-distribution state with uniform initial conditions.
    pub fn new(nx: usize, ny: usize, tau: f64, tau_t: f64, rho0: f64, t0: f64) -> Self {
        let n = nx * ny;
        let f_eq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * rho0);
        let h_eq: [f64; 5] = D2Q5_WEIGHTS.map(|w| w * t0);
        Self {
            f: vec![f_eq; n],
            h: vec![h_eq; n],
            rho: vec![rho0; n],
            u: vec![[0.0; 2]; n],
            temperature: vec![t0; n],
            nx,
            ny,
            tau,
            tau_t,
        }
    }

    /// Flat node index.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// D2Q9 equilibrium distribution.
    pub fn f_equilibrium(&self, rho: f64, u: [f64; 2], i: usize) -> f64 {
        let w = D2Q9_WEIGHTS[i];
        let c = D2Q9_VELOCITIES[i];
        let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
        let u2 = u[0] * u[0] + u[1] * u[1];
        w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }

    /// D2Q5 equilibrium distribution.
    pub fn h_equilibrium(&self, temp: f64, u: [f64; 2], i: usize) -> f64 {
        let w = D2Q5_WEIGHTS[i];
        let c = D2Q5_VELOCITIES[i];
        let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
        w * temp * (1.0 + eu / CS2)
    }

    /// BGK collision for both fluid and temperature populations on every node.
    pub fn collide(&mut self) {
        self.collide_nodes(|_| true);
    }

    /// BGK collision on the fluid nodes (`is_fluid[k] == true`) only; the populations of solid
    /// nodes are left untouched. Nodes beyond the end of `is_fluid` count as solid.
    pub fn collide_masked(&mut self, is_fluid: &[bool]) {
        self.collide_nodes(|k| is_fluid.get(k).copied().unwrap_or(false));
    }

    fn collide_nodes<P: Fn(usize) -> bool>(&mut self, is_active: P) {
        let n = self.nx * self.ny;
        for k in 0..n {
            if !is_active(k) {
                continue;
            }
            let rho = self.rho[k];
            let u = self.u[k];
            let t = self.temperature[k];
            for i in 0..9 {
                let feq = self.f_equilibrium(rho, u, i);
                self.f[k][i] -= (self.f[k][i] - feq) / self.tau;
            }
            for i in 0..5 {
                let heq = self.h_equilibrium(t, u, i);
                self.h[k][i] -= (self.h[k][i] - heq) / self.tau_t;
            }
        }
    }

    /// Streaming step with periodic boundaries for both populations on every node.
    pub fn stream(&mut self) {
        self.stream_nodes(|_| true, &[]);
    }

    /// Streaming step on the fluid nodes with solid links (`is_fluid[k] == false` neighbours).
    ///
    /// A population that would stream from a fluid node into a solid node is reflected back
    /// into the fluid node with the reversed direction (half-way bounce-back), so no mass or
    /// heat is exchanged with the solid by the lattice itself and the flow sees a no-slip wall.
    /// For the temperature populations `wall_temperature[k][i] = Some(T_w)` turns the link of
    /// fluid node `k` in D2Q5 direction `i` into a Dirichlet wall of temperature `T_w`
    /// (anti-bounce-back `h_ī = −h*_i + 2 w_i T_w`); `None` keeps the plain bounce-back, i.e. an
    /// adiabatic wall. The populations of solid nodes are left untouched. Nodes beyond the end
    /// of `is_fluid` count as solid.
    pub fn stream_masked(&mut self, is_fluid: &[bool], wall_temperature: &[[Option<f64>; 5]]) {
        self.stream_nodes(
            |k| is_fluid.get(k).copied().unwrap_or(false),
            wall_temperature,
        );
    }

    fn stream_nodes<P: Fn(usize) -> bool>(
        &mut self,
        is_active: P,
        wall_temperature: &[[Option<f64>; 5]],
    ) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = self.f.clone();
        let mut h_new = self.h.clone();

        for y in 0..ny {
            for x in 0..nx {
                let k_src = self.idx(x, y);
                if !is_active(k_src) {
                    continue;
                }
                for i in 0..9 {
                    let c = D2Q9_VELOCITIES[i];
                    let xd = ((x as isize + c[0] as isize).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + c[1] as isize).rem_euclid(ny as isize)) as usize;
                    let k_dst = self.idx(xd, yd);
                    if is_active(k_dst) {
                        f_new[k_dst][i] = self.f[k_src][i];
                    } else {
                        f_new[k_src][D2Q9_OPPOSITES[i]] = self.f[k_src][i];
                    }
                }
                for i in 0..5 {
                    let c = D2Q5_VELOCITIES[i];
                    let xd = ((x as isize + c[0] as isize).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + c[1] as isize).rem_euclid(ny as isize)) as usize;
                    let k_dst = self.idx(xd, yd);
                    if is_active(k_dst) {
                        h_new[k_dst][i] = self.h[k_src][i];
                    } else {
                        let wall = wall_temperature.get(k_src).and_then(|links| links[i]);
                        h_new[k_src][D2Q5_OPPOSITES[i]] = match wall {
                            Some(t_wall) => -self.h[k_src][i] + 2.0 * D2Q5_WEIGHTS[i] * t_wall,
                            None => self.h[k_src][i],
                        };
                    }
                }
            }
        }
        self.f = f_new;
        self.h = h_new;
    }

    /// Compute macroscopic fields from moments on every node.
    pub fn compute_macroscopic(&mut self) {
        self.macroscopic_nodes(|_| true);
    }

    /// Compute the macroscopic fields of the fluid nodes only; solid nodes keep their values.
    pub fn compute_macroscopic_masked(&mut self, is_fluid: &[bool]) {
        self.macroscopic_nodes(|k| is_fluid.get(k).copied().unwrap_or(false));
    }

    fn macroscopic_nodes<P: Fn(usize) -> bool>(&mut self, is_active: P) {
        let n = self.nx * self.ny;
        for k in 0..n {
            if !is_active(k) {
                continue;
            }
            let rho: f64 = self.f[k].iter().sum();
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for (f_ki, c) in self.f[k].iter().zip(D2Q9_VELOCITIES.iter()) {
                ux += c[0] as f64 * f_ki;
                uy += c[1] as f64 * f_ki;
            }
            self.rho[k] = rho;
            self.u[k] = [ux / rho, uy / rho];
            self.temperature[k] = self.h[k].iter().sum();
        }
    }
}

// ---------------------------------------------------------------------------
// ThermalBoundary
// ---------------------------------------------------------------------------

/// Types of thermal boundary conditions.
#[derive(Debug, Clone)]
pub enum ThermalBoundary {
    /// Dirichlet: fixed temperature T_wall.
    ConstantTemperature(f64),
    /// Neumann: fixed heat flux q_wall (W/m²).
    ConstantFlux(f64),
    /// Adiabatic: zero heat flux (q = 0).
    Adiabatic,
    /// Robin: mixed h*(T - T_inf) type condition; h_conv and T_inf.
    Robin {
        /// Convective heat transfer coefficient (W/m²·K).
        h_conv: f64,
        /// Far-field temperature (K).
        t_inf: f64,
    },
}

impl ThermalBoundary {
    /// Apply this boundary condition to a D2Q5 node.
    ///
    /// Returns the effective wall temperature for distribution initialisation.
    pub fn effective_temperature(&self, t_node: f64, alpha: f64) -> f64 {
        match self {
            Self::ConstantTemperature(t_w) => *t_w,
            Self::ConstantFlux(q) => t_node + q * alpha,
            Self::Adiabatic => t_node,
            Self::Robin { h_conv, t_inf } => t_node + h_conv * (t_inf - t_node) * alpha,
        }
    }
}

// ---------------------------------------------------------------------------
// PcmMaterial
// ---------------------------------------------------------------------------

/// Phase-change material (PCM) properties for latent heat treatment.
///
/// The enthalpy method is used:
/// ```text
/// H = c_p * T                         (T < T_melt)
/// H = c_p * T_melt + L * f_l           (T = T_melt, 0 ≤ f_l ≤ 1)
/// H = c_p * T_melt + L + c_p*(T-T_melt) (T > T_melt)
/// ```
/// where f_l is the liquid fraction.
#[derive(Debug, Clone)]
pub struct PcmMaterial {
    /// Latent heat of fusion L (J/kg).
    pub latent_heat: f64,
    /// Melting temperature T_melt (K).
    pub melting_temperature: f64,
    /// Mushy zone half-width ΔT (K); melting occurs over \[T_melt - ΔT, T_melt + ΔT\].
    pub mushy_zone_width: f64,
    /// Solid-phase specific heat c_p,s (J/kg·K).
    pub cp_solid: f64,
    /// Liquid-phase specific heat c_p,l (J/kg·K).
    pub cp_liquid: f64,
    /// Solid density ρ_s (kg/m³).
    pub density_solid: f64,
    /// Liquid density ρ_l (kg/m³).
    pub density_liquid: f64,
    /// Solid thermal conductivity k_s (W/m·K).
    pub k_solid: f64,
    /// Liquid thermal conductivity k_l (W/m·K).
    pub k_liquid: f64,
}

impl PcmMaterial {
    /// Create a new PCM material specification.
    pub fn new(
        latent_heat: f64,
        melting_temperature: f64,
        mushy_zone_width: f64,
        cp_solid: f64,
        cp_liquid: f64,
        density_solid: f64,
        density_liquid: f64,
        k_solid: f64,
        k_liquid: f64,
    ) -> Self {
        Self {
            latent_heat,
            melting_temperature,
            mushy_zone_width,
            cp_solid,
            cp_liquid,
            density_solid,
            density_liquid,
            k_solid,
            k_liquid,
        }
    }

    /// Liquid fraction f_l ∈ \[0, 1\] for a given temperature.
    pub fn liquid_fraction(&self, temperature: f64) -> f64 {
        let t_low = self.melting_temperature - self.mushy_zone_width;
        let t_high = self.melting_temperature + self.mushy_zone_width;
        if temperature <= t_low {
            0.0
        } else if temperature >= t_high {
            1.0
        } else {
            (temperature - t_low) / (2.0 * self.mushy_zone_width)
        }
    }

    /// Effective specific heat (includes latent heat contribution in mushy zone).
    pub fn effective_cp(&self, temperature: f64) -> f64 {
        let fl = self.liquid_fraction(temperature);
        if fl <= 0.0 {
            self.cp_solid
        } else if fl >= 1.0 {
            self.cp_liquid
        } else {
            // Mushy zone: blend + latent heat smeared over 2*ΔT
            let cp_blend = (1.0 - fl) * self.cp_solid + fl * self.cp_liquid;
            cp_blend + self.latent_heat / (2.0 * self.mushy_zone_width)
        }
    }

    /// Effective thermal conductivity (blend between solid and liquid).
    pub fn effective_conductivity(&self, temperature: f64) -> f64 {
        let fl = self.liquid_fraction(temperature);
        (1.0 - fl) * self.k_solid + fl * self.k_liquid
    }

    /// Enthalpy H (J/kg) at a given temperature.
    pub fn enthalpy(&self, temperature: f64) -> f64 {
        let fl = self.liquid_fraction(temperature);
        self.cp_solid * temperature.min(self.melting_temperature)
            + fl * self.latent_heat
            + if temperature > self.melting_temperature {
                self.cp_liquid * (temperature - self.melting_temperature)
            } else {
                0.0
            }
    }
}

// ---------------------------------------------------------------------------
// NaturalConvection
// ---------------------------------------------------------------------------

/// Natural convection parameters for Boussinesq-approximated LBM.
///
/// The buoyancy force is applied as a body force to the fluid momentum:
/// ```text
/// F_y = -ρ g β (T - T_ref)
/// ```
#[derive(Debug, Clone)]
pub struct NaturalConvection {
    /// Gravitational acceleration magnitude (m/s² in lattice units).
    pub g: f64,
    /// Thermal expansion coefficient β (1/K).
    pub beta: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Hot wall temperature T_h (K).
    pub t_hot: f64,
    /// Cold wall temperature T_c (K).
    pub t_cold: f64,
    /// Rayleigh number Ra = g β ΔT L³ / (ν α).
    pub rayleigh: f64,
    /// Prandtl number Pr = ν / α.
    pub prandtl: f64,
}

impl NaturalConvection {
    /// Create a Rayleigh-Bénard natural convection configuration.
    ///
    /// `length` is the characteristic length L (cavity height in lattice units).
    /// `nu` is kinematic viscosity, `alpha` is thermal diffusivity.
    pub fn new(
        g: f64,
        beta: f64,
        t_hot: f64,
        t_cold: f64,
        length: f64,
        nu: f64,
        alpha: f64,
    ) -> Self {
        let dt = t_hot - t_cold;
        let t_ref = (t_hot + t_cold) * 0.5;
        let rayleigh = g * beta * dt * length.powi(3) / (nu * alpha);
        let prandtl = nu / alpha;
        Self {
            g,
            beta,
            t_ref,
            t_hot,
            t_cold,
            rayleigh,
            prandtl,
        }
    }

    /// Compute Boussinesq body force in y-direction for given temperature.
    pub fn buoyancy_force(&self, temperature: f64) -> f64 {
        boussinesq_force(self.g, self.beta, temperature, self.t_ref)
    }
}

// ---------------------------------------------------------------------------
// ForcedConvection
// ---------------------------------------------------------------------------

/// Forced convection parameters for channel flow with heated wall.
#[derive(Debug, Clone)]
pub struct ForcedConvection {
    /// Bulk velocity U_bulk (lattice units).
    pub u_bulk: f64,
    /// Channel half-height H (lattice units).
    pub half_height: f64,
    /// Wall temperature T_wall (K).
    pub t_wall: f64,
    /// Inlet temperature T_in (K).
    pub t_in: f64,
    /// Hydraulic diameter D_h = 2H for a channel.
    pub hydraulic_diameter: f64,
    /// Kinematic viscosity ν.
    pub nu: f64,
    /// Thermal diffusivity α.
    pub alpha: f64,
}

impl ForcedConvection {
    /// Create a new forced convection configuration.
    pub fn new(u_bulk: f64, half_height: f64, t_wall: f64, t_in: f64, nu: f64, alpha: f64) -> Self {
        let hydraulic_diameter = 4.0 * half_height;
        Self {
            u_bulk,
            half_height,
            t_wall,
            t_in,
            hydraulic_diameter,
            nu,
            alpha,
        }
    }

    /// Reynolds number Re = U_bulk * D_h / ν.
    pub fn reynolds(&self) -> f64 {
        self.u_bulk * self.hydraulic_diameter / self.nu
    }

    /// Prandtl number Pr = ν / α.
    pub fn prandtl(&self) -> f64 {
        prandtl_number(self.nu, self.alpha)
    }

    /// Bulk mean Nusselt number using Dittus-Boelter correlation (turbulent):
    /// Nu = 0.023 Re^0.8 Pr^0.4.
    pub fn nusselt_dittus_boelter(&self) -> f64 {
        let re = self.reynolds();
        let pr = self.prandtl();
        0.023 * re.powf(0.8) * pr.powf(0.4)
    }

    /// Laminar Nusselt number for thermally developed channel flow: Nu = 7.54.
    pub fn nusselt_laminar_channel(&self) -> f64 {
        7.54
    }

    /// Convective heat transfer coefficient h = Nu * k / D_h.
    pub fn convection_coefficient(&self, k_fluid: f64) -> f64 {
        let nu = self.nusselt_dittus_boelter();
        nu * k_fluid / self.hydraulic_diameter
    }
}

// ---------------------------------------------------------------------------
// ConjugateResult
// ---------------------------------------------------------------------------

/// Results from a conjugate heat transfer simulation step.
#[derive(Debug, Clone)]
pub struct ConjugateResult {
    /// Temperature field at all nodes (flat, fluid+solid).
    pub temperature: Vec<f64>,
    /// Heat flux at interface nodes (W/m²).
    pub interface_heat_flux: Vec<f64>,
    /// Average Nusselt number over heated surface.
    pub nusselt_number: f64,
    /// Average convective heat transfer coefficient h (W/m²·K).
    pub convective_coefficient: f64,
    /// Maximum temperature in domain.
    pub t_max: f64,
    /// Minimum temperature in domain.
    pub t_min: f64,
}

impl ConjugateResult {
    /// Construct a result from a temperature field and interface data.
    pub fn from_fields(
        temperature: Vec<f64>,
        interface_heat_flux: Vec<f64>,
        nusselt_number: f64,
        convective_coefficient: f64,
    ) -> Self {
        let t_max = temperature
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let t_min = temperature.iter().cloned().fold(f64::INFINITY, f64::min);
        Self {
            temperature,
            interface_heat_flux,
            nusselt_number,
            convective_coefficient,
            t_max,
            t_min,
        }
    }

    /// Mean temperature over entire domain.
    pub fn mean_temperature(&self) -> f64 {
        if self.temperature.is_empty() {
            return 0.0;
        }
        self.temperature.iter().sum::<f64>() / self.temperature.len() as f64
    }
}

// ---------------------------------------------------------------------------
// ConjugateHeatSolver
// ---------------------------------------------------------------------------

/// Full conjugate heat transfer solver using a dual grid.
///
/// - Fluid domain: D2Q9 for velocity + D2Q5 for temperature (LBM) on the fluid nodes only;
///   every fluid–solid link bounces the flow back (no slip) and reflects the temperature
///   populations (adiabatic, or Dirichlet at the interface temperature where an interface is
///   registered).
/// - Solid domain: explicit finite-difference diffusion on the same grid. The outermost ring of
///   grid nodes is never updated and acts as a fixed-temperature (Dirichlet) boundary at its
///   initial value.
/// - Interface nodes track temperature and enforce flux continuity (see the module
///   documentation); the node temperatures are not overwritten.
///
/// The solid temperature field uses the physical `dx`/`dt` of the solver, the fluid lattice
/// works in lattice units; for a consistent setup express the solid properties in lattice
/// units too (`dx = dt = 1`) and use `k_fluid = ρ c_p α_f` with `α_f = (τ_t − ½)/3` for the
/// fluid conductivity of the interfaces.
#[derive(Debug, Clone)]
pub struct ConjugateHeatSolver {
    /// Double-distribution LBM state.
    pub lbm: DoubleDistribution,
    /// Temperature in solid nodes (flat array, same size as lbm grid).
    pub solid_temp: Vec<f64>,
    /// Boolean mask: true = fluid node, false = solid node.
    pub is_fluid: Vec<bool>,
    /// Interface node pairs (fluid_idx, solid_idx, normal).
    pub interfaces: Vec<FluidSolidInterface>,
    /// Solid material properties.
    pub solid: SolidRegion,
    /// Lattice spacing Δx (for physical heat flux computation).
    pub dx: f64,
    /// Time step Δt (lattice time units).
    pub dt: f64,
    /// Current simulation time.
    pub time: f64,
    /// Step counter.
    pub step: u64,
}

/// Largest explicit-diffusion number `α Δt / Δx²` of the 5-point solid stencil in 2-D.
const MAX_SOLID_DIFFUSION_NUMBER: f64 = 0.25;

/// D2Q5 direction from node `from` to its axis-adjacent (periodic) neighbour `to`, if any.
fn link_direction(nx: usize, ny: usize, from: usize, to: usize) -> Option<usize> {
    if nx == 0 || ny == 0 || from >= nx * ny || to >= nx * ny {
        return None;
    }
    let (x, y) = (from % nx, from / nx);
    for (direction, c) in D2Q5_VELOCITIES.iter().enumerate().skip(1) {
        let xn = ((x as isize + c[0] as isize).rem_euclid(nx as isize)) as usize;
        let yn = ((y as isize + c[1] as isize).rem_euclid(ny as isize)) as usize;
        if yn * nx + xn == to {
            return Some(direction);
        }
    }
    None
}

impl ConjugateHeatSolver {
    /// Create a new conjugate heat solver.
    ///
    /// `is_fluid[k] = true` marks fluid nodes; solid nodes use finite-difference diffusion.
    ///
    /// # Errors
    ///
    /// Returns an error if the grid is empty, `is_fluid` does not hold `nx * ny` flags, `tau`
    /// or `tau_t` is not above ½, `dx` or `dt` is not positive and finite, or the explicit solid
    /// diffusion number `α_s Δt / Δx²` exceeds ¼ (the stability limit of the 5-point stencil).
    pub fn new(
        nx: usize,
        ny: usize,
        tau: f64,
        tau_t: f64,
        rho0: f64,
        t0: f64,
        is_fluid: Vec<bool>,
        solid: SolidRegion,
        dx: f64,
        dt: f64,
    ) -> Result<Self> {
        Error::check_grid(nx, ny)?;
        Error::check_tau(tau)?;
        Error::check_tau(tau_t)?;
        if is_fluid.len() != nx * ny {
            return Err(Error::ComponentCountMismatch {
                expected: nx * ny,
                found: is_fluid.len(),
                field: "is_fluid".to_string(),
            });
        }
        Error::check_range(dx, f64::MIN_POSITIVE, f64::MAX, "dx")?;
        Error::check_range(dt, f64::MIN_POSITIVE, f64::MAX, "dt")?;
        let diffusion_number = solid.thermal_diffusivity() * dt / (dx * dx);
        Error::check_range(
            diffusion_number,
            0.0,
            MAX_SOLID_DIFFUSION_NUMBER,
            "solid diffusion number alpha_s*dt/dx^2",
        )?;
        let lbm = DoubleDistribution::new(nx, ny, tau, tau_t, rho0, t0);
        let solid_temp = vec![t0; nx * ny];
        Ok(Self {
            lbm,
            solid_temp,
            is_fluid,
            interfaces: Vec::new(),
            solid,
            dx,
            dt,
            time: 0.0,
            step: 0,
        })
    }

    /// Explicit solid diffusion number `r = α_s Δt / Δx²`.
    pub fn solid_diffusion_number(&self) -> f64 {
        self.solid.thermal_diffusivity() * self.dt / (self.dx * self.dx)
    }

    /// Add an interface between a fluid node and an axis-adjacent solid node.
    ///
    /// # Errors
    ///
    /// Returns an error if a node index is outside the grid, `fluid_index` is not a fluid node,
    /// `solid_index` is not a solid node, the two nodes are not axis-adjacent lattice
    /// neighbours, a conductivity is not positive and finite, or the interface makes the
    /// explicit solid update unstable: with the ghost value `2 T_i − T_s` the solid node keeps
    /// the weight `1 − r (4 − β)` of its own temperature, `β = (k_s − k_f)/(k_s + k_f)`, which
    /// must not be negative.
    pub fn add_interface(&mut self, iface: FluidSolidInterface) -> Result<()> {
        let n = self.lbm.nx * self.lbm.ny;
        if iface.fluid_index >= n || iface.solid_index >= n {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "interface node index out of range: fluid {}, solid {} on a grid of {n} nodes",
                    iface.fluid_index, iface.solid_index
                ),
            });
        }
        if !self
            .is_fluid
            .get(iface.fluid_index)
            .copied()
            .unwrap_or(false)
        {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "interface fluid node {} is not a fluid node",
                    iface.fluid_index
                ),
            });
        }
        if self
            .is_fluid
            .get(iface.solid_index)
            .copied()
            .unwrap_or(true)
        {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "interface solid node {} is not a solid node",
                    iface.solid_index
                ),
            });
        }
        if link_direction(
            self.lbm.nx,
            self.lbm.ny,
            iface.fluid_index,
            iface.solid_index,
        )
        .is_none()
        {
            return Err(Error::IncompatibleParameters {
                message: format!(
                    "interface nodes {} and {} are not axis-adjacent lattice neighbours",
                    iface.fluid_index, iface.solid_index
                ),
            });
        }
        Error::check_range(
            iface.k_fluid,
            f64::MIN_POSITIVE,
            f64::MAX,
            "interface k_fluid",
        )?;
        Error::check_range(
            iface.k_solid,
            f64::MIN_POSITIVE,
            f64::MAX,
            "interface k_solid",
        )?;
        let beta = (iface.k_solid - iface.k_fluid) / (iface.k_solid + iface.k_fluid);
        Error::check_range(
            self.solid_diffusion_number(),
            0.0,
            1.0 / (4.0 - beta.min(0.0)),
            "solid diffusion number alpha_s*dt/dx^2 at the interface",
        )?;
        self.interfaces.push(iface);
        Ok(())
    }

    /// Perform one full conjugate heat transfer step.
    ///
    /// 1. Interface temperatures from the current node temperatures
    /// 2. LBM collide + stream on the fluid nodes (bounce-back / anti-bounce-back at the solid links)
    /// 3. Solid finite-difference diffusion with the interface temperatures as ghost values
    /// 4. Refresh the interface temperatures from the updated node temperatures
    pub fn step(&mut self) {
        // 1. Interface temperatures at the start of the step.
        self.update_interface_temperatures();
        let walls = self.wall_links();

        // 2. LBM fluid step on the fluid nodes only.
        self.lbm.collide_masked(&self.is_fluid);
        self.lbm.stream_masked(&self.is_fluid, &walls);
        self.lbm.compute_macroscopic_masked(&self.is_fluid);

        // 3. Solid diffusion (explicit FD, 5-point stencil).
        self.diffuse_solid(&walls);

        // 4. Interface temperatures consistent with the new node temperatures.
        self.update_interface_temperatures();

        self.time += self.dt;
        self.step += 1;
    }

    /// Interface temperatures of the registered interfaces from the current node temperatures.
    fn update_interface_temperatures(&mut self) {
        for iface in &mut self.interfaces {
            if let (Some(&t_fluid), Some(&t_solid)) = (
                self.lbm.temperature.get(iface.fluid_index),
                self.solid_temp.get(iface.solid_index),
            ) {
                iface.compute_interface_temperature(t_fluid, t_solid);
            }
        }
    }

    /// Dirichlet wall temperature of every fluid link that carries a registered interface:
    /// `wall_links()[fluid_node][d2q5_direction] = Some(interface temperature)`.
    fn wall_links(&self) -> Vec<[Option<f64>; 5]> {
        let (nx, ny) = (self.lbm.nx, self.lbm.ny);
        let mut walls = vec![[None; 5]; nx * ny];
        for iface in &self.interfaces {
            let fluid_ok = self
                .is_fluid
                .get(iface.fluid_index)
                .copied()
                .unwrap_or(false);
            let solid_ok = !self
                .is_fluid
                .get(iface.solid_index)
                .copied()
                .unwrap_or(true);
            if !(fluid_ok && solid_ok) {
                continue;
            }
            if let Some(direction) = link_direction(nx, ny, iface.fluid_index, iface.solid_index) {
                walls[iface.fluid_index][direction] = Some(iface.interface_temperature);
            }
        }
        walls
    }

    /// Explicit finite-difference heat diffusion in solid nodes.
    ///
    /// A fluid neighbour contributes the ghost value `2 T_i − T_s` of its registered interface
    /// (flux `2 k_s (T_i − T_s)/Δx` through the half-cell) or, without an interface, the node's
    /// own temperature (adiabatic link).
    fn diffuse_solid(&mut self, walls: &[[Option<f64>; 5]]) {
        let nx = self.lbm.nx;
        let ny = self.lbm.ny;
        let r = self.solid_diffusion_number();

        let old = self.solid_temp.clone();
        for y in 1..ny.saturating_sub(1) {
            for x in 1..nx.saturating_sub(1) {
                let k = y * nx + x;
                if self.is_fluid[k] {
                    continue;
                }
                let t_c = old[k];
                // Neighbour value seen from the solid node through direction `d` (D2Q5 index
                // pointing from the solid node to the neighbour node `idx`).
                let neighbour = |idx: usize, d: usize| -> f64 {
                    if self.is_fluid[idx] {
                        match walls.get(idx).and_then(|links| links[D2Q5_OPPOSITES[d]]) {
                            Some(t_interface) => 2.0 * t_interface - t_c,
                            None => t_c,
                        }
                    } else {
                        old[idx]
                    }
                };
                let t_e = neighbour(k + 1, 1);
                let t_n = neighbour(k + nx, 2);
                let t_w = neighbour(k - 1, 3);
                let t_s = neighbour(k - nx, 4);
                self.solid_temp[k] = t_c + r * (t_e + t_w + t_n + t_s - 4.0 * t_c);
            }
        }
    }

    /// Temperature of every node: the LBM temperature on fluid nodes, the finite-difference
    /// temperature on solid nodes.
    pub fn temperature_field(&self) -> Vec<f64> {
        (0..self.lbm.nx * self.lbm.ny)
            .map(|k| {
                if self.is_fluid.get(k).copied().unwrap_or(false) {
                    self.lbm.temperature[k]
                } else {
                    self.solid_temp[k]
                }
            })
            .collect()
    }

    /// Compute and return simulation results.
    ///
    /// The interface heat fluxes are positive when heat flows from the solid into the fluid
    /// (see [`FluidSolidInterface::heat_flux`]). The convection coefficient is
    /// `h = |⟨q⟩| / |t_wall − t_ref|` with the mean `⟨q⟩` of the signed fluxes of the registered
    /// interfaces (register only the heated surface, or interfaces that transfer heat in
    /// opposite directions offset each other), and the Nusselt number is `h L / k_fluid` with
    /// the domain width `L`.
    pub fn results(&self, t_ref: f64, t_wall: f64, k_fluid: f64) -> ConjugateResult {
        let temperature = self.temperature_field();
        let interface_heat_flux: Vec<f64> = self
            .interfaces
            .iter()
            .map(|iface| {
                match (
                    self.lbm.temperature.get(iface.fluid_index),
                    self.solid_temp.get(iface.solid_index),
                ) {
                    (Some(&t_fluid), Some(&t_solid)) => iface.heat_flux(t_fluid, t_solid, self.dx),
                    _ => 0.0,
                }
            })
            .collect();

        let mean_flux = if interface_heat_flux.is_empty() {
            0.0
        } else {
            interface_heat_flux.iter().sum::<f64>() / interface_heat_flux.len() as f64
        };
        let dt_drive = (t_wall - t_ref).abs().max(1e-12);
        let h_conv = mean_flux.abs() / dt_drive;
        let nu = nusselt_number_forced(h_conv, self.lbm.nx as f64 * self.dx, k_fluid);

        ConjugateResult::from_fields(temperature, interface_heat_flux, nu, h_conv)
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Boussinesq buoyancy body force: F_y = −g β (T − T_ref).
///
/// Returns the force per unit volume in the gravity direction (y positive up).
pub fn boussinesq_force(g: f64, beta: f64, temperature: f64, t_ref: f64) -> f64 {
    -g * beta * (temperature - t_ref)
}

/// Forced convection Nusselt number: Nu = h L / k.
///
/// `h` is convective coefficient, `L` characteristic length, `k` fluid conductivity.
pub fn nusselt_number_forced(h: f64, length: f64, k_fluid: f64) -> f64 {
    h * length / k_fluid
}

/// Natural convection Nusselt number from Churchill-Chu correlation for a vertical plate:
///
/// ```text
/// Nu = (0.825 + 0.387 Ra^(1/6) / (1 + (0.492/Pr)^(9/16))^(8/27))^2
/// ```
pub fn nusselt_number_natural(rayleigh: f64, prandtl: f64) -> f64 {
    let bracket = 1.0 + (0.492 / prandtl).powf(9.0 / 16.0);
    let denom = bracket.powf(8.0 / 27.0);
    let inner = 0.825 + 0.387 * rayleigh.powf(1.0 / 6.0) / denom;
    inner * inner
}

/// Prandtl number: Pr = ν / α.
pub fn prandtl_number(nu: f64, alpha: f64) -> f64 {
    nu / alpha
}

/// Stanton number: St = Nu / (Re Pr).
pub fn stanton_number(nu: f64, re: f64, pr: f64) -> f64 {
    nu / (re * pr)
}

/// Biot number: Bi = h L / k_solid.
pub fn biot_number(h_conv: f64, length: f64, k_solid: f64) -> f64 {
    h_conv * length / k_solid
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod interface_tests;

#[cfg(test)]
mod tests {
    use super::*;

    // --- SolidRegion ---

    #[test]
    fn test_solid_region_diffusivity() {
        let s = SolidRegion::new(1.0, 500.0, 8000.0, "steel");
        let alpha = s.thermal_diffusivity();
        assert!((alpha - 1.0 / (500.0 * 8000.0)).abs() < 1e-15);
    }

    #[test]
    fn test_solid_region_effusivity() {
        let s = SolidRegion::new(1.0, 500.0, 8000.0, "steel");
        let e = s.thermal_effusivity();
        assert!((e - (1.0 * 500.0 * 8000.0_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_solid_region_label() {
        let s = SolidRegion::new(0.6, 4182.0, 997.0, "water");
        assert_eq!(s.label, "water");
    }

    // --- FluidSolidInterface ---

    #[test]
    fn test_interface_temperature() {
        let mut iface = FluidSolidInterface::new(0, 1, [1.0, 0.0], 0.6, 16.0);
        iface.compute_interface_temperature(300.0, 400.0);
        let expected = (0.6 * 300.0 + 16.0 * 400.0) / (0.6 + 16.0);
        assert!((iface.interface_temperature - expected).abs() < 1e-10);
    }

    #[test]
    fn test_flux_residual_zero() {
        let iface = FluidSolidInterface::new(0, 1, [1.0, 0.0], 1.0, 1.0);
        let res = iface.flux_residual(5.0, 5.0);
        assert!((res).abs() < 1e-15);
    }

    #[test]
    fn test_flux_residual_nonzero() {
        let iface = FluidSolidInterface::new(0, 1, [0.0, 1.0], 2.0, 4.0);
        let res = iface.flux_residual(3.0, 1.0);
        // 2*3 - 4*1 = 6 - 4 = 2
        assert!((res - 2.0).abs() < 1e-10);
    }

    // --- ThermalLbmD2Q5 ---

    #[test]
    fn test_d2q5_new_temperature_uniform() {
        let solver = ThermalLbmD2Q5::new(4, 4, 1.0, 1.5);
        for &t in &solver.temperature {
            assert!((t - 1.5).abs() < 1e-10);
        }
    }

    #[test]
    fn test_d2q5_equilibrium_rest() {
        let solver = ThermalLbmD2Q5::new(4, 4, 1.0, 1.0);
        // At rest, h_eq_0 = w_0 * T = 1/3 * 1.0
        let heq = solver.equilibrium(1.0, [0.0, 0.0], 0);
        assert!((heq - 1.0 / 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_d2q5_equilibrium_sum_equals_temperature() {
        let solver = ThermalLbmD2Q5::new(4, 4, 1.0, 2.5);
        let t = 2.5;
        let u = [0.05, 0.02];
        let sum: f64 = (0..5).map(|i| solver.equilibrium(t, u, i)).sum();
        assert!((sum - t).abs() < 1e-12);
    }

    #[test]
    fn test_d2q5_alpha() {
        let solver = ThermalLbmD2Q5::new(4, 4, 1.0, 1.0);
        let alpha = solver.alpha();
        assert!((alpha - (1.0 - 0.5) / 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_d2q5_collide_stream_temperature_conserved() {
        let mut solver = ThermalLbmD2Q5::new(8, 8, 1.0, 1.0);
        let vels = vec![[0.0_f64; 2]; 64];
        solver.collide(&vels);
        solver.stream();
        solver.compute_temperature();
        let sum: f64 = solver.temperature.iter().sum();
        // Total temperature should be conserved (= 64 nodes * 1.0)
        assert!((sum - 64.0).abs() < 1e-8);
    }

    #[test]
    fn test_d2q5_dirichlet() {
        let mut solver = ThermalLbmD2Q5::new(4, 4, 1.0, 1.0);
        solver.apply_dirichlet(0, 2.0, [0.0, 0.0]);
        assert!((solver.temperature[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_d2q5_idx() {
        let solver = ThermalLbmD2Q5::new(4, 5, 1.0, 1.0);
        assert_eq!(solver.idx(0, 0), 0);
        assert_eq!(solver.idx(3, 4), 4 * 4 + 3);
    }

    // --- DoubleDistribution ---

    #[test]
    fn test_double_distribution_init() {
        let dd = DoubleDistribution::new(4, 4, 1.0, 1.0, 1.0, 1.0);
        assert_eq!(dd.rho.len(), 16);
        assert_eq!(dd.temperature.len(), 16);
        for &r in &dd.rho {
            assert!((r - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_double_distribution_collide_stream() {
        let mut dd = DoubleDistribution::new(8, 8, 1.0, 1.0, 1.0, 1.0);
        dd.collide();
        dd.stream();
        dd.compute_macroscopic();
        let sum_rho: f64 = dd.rho.iter().sum();
        assert!((sum_rho - 64.0).abs() < 1e-8);
    }

    #[test]
    fn test_double_distribution_h_equilibrium_sum() {
        let dd = DoubleDistribution::new(4, 4, 1.0, 1.0, 1.0, 2.0);
        let t = 2.0;
        let u = [0.0, 0.0];
        let sum: f64 = (0..5).map(|i| dd.h_equilibrium(t, u, i)).sum();
        assert!((sum - t).abs() < 1e-12);
    }

    // --- ThermalBoundary ---

    #[test]
    fn test_thermal_boundary_dirichlet() {
        let bc = ThermalBoundary::ConstantTemperature(300.0);
        let t_eff = bc.effective_temperature(250.0, 0.1);
        assert!((t_eff - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_boundary_adiabatic() {
        let bc = ThermalBoundary::Adiabatic;
        let t_eff = bc.effective_temperature(300.0, 0.1);
        assert!((t_eff - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_boundary_flux() {
        let bc = ThermalBoundary::ConstantFlux(100.0);
        let t_eff = bc.effective_temperature(300.0, 0.01);
        // 300 + 100 * 0.01 = 301
        assert!((t_eff - 301.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_boundary_robin() {
        let bc = ThermalBoundary::Robin {
            h_conv: 10.0,
            t_inf: 400.0,
        };
        // T_eff = 300 + 10*(400-300)*0.1 = 300 + 100 = 400... wait
        // = 300 + 10*(400-300)*0.1 = 300 + 10*100*0.1 = 300 + 100 = 400
        let t_eff = bc.effective_temperature(300.0, 0.1);
        assert!((t_eff - 400.0).abs() < 1e-10);
    }

    // --- PcmMaterial ---

    #[test]
    fn test_pcm_liquid_fraction_below_melt() {
        let pcm = PcmMaterial::new(
            200_000.0, 300.0, 5.0, 2000.0, 2100.0, 900.0, 850.0, 0.2, 0.15,
        );
        assert!((pcm.liquid_fraction(290.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_pcm_liquid_fraction_above_melt() {
        let pcm = PcmMaterial::new(
            200_000.0, 300.0, 5.0, 2000.0, 2100.0, 900.0, 850.0, 0.2, 0.15,
        );
        assert!((pcm.liquid_fraction(310.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pcm_liquid_fraction_at_melt() {
        let pcm = PcmMaterial::new(
            200_000.0, 300.0, 5.0, 2000.0, 2100.0, 900.0, 850.0, 0.2, 0.15,
        );
        // At T_melt, mushy zone midpoint -> f_l = 0.5
        assert!((pcm.liquid_fraction(300.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_pcm_enthalpy_increases_with_temperature() {
        let pcm = PcmMaterial::new(
            200_000.0, 300.0, 5.0, 2000.0, 2100.0, 900.0, 850.0, 0.2, 0.15,
        );
        let h1 = pcm.enthalpy(280.0);
        let h2 = pcm.enthalpy(320.0);
        assert!(h2 > h1);
    }

    #[test]
    fn test_pcm_effective_conductivity_blend() {
        let pcm = PcmMaterial::new(
            200_000.0, 300.0, 5.0, 2000.0, 2100.0, 900.0, 850.0, 0.2, 0.15,
        );
        let k = pcm.effective_conductivity(300.0);
        // At midpoint: 0.5*0.2 + 0.5*0.15 = 0.175
        assert!((k - 0.175).abs() < 1e-12);
    }

    // --- NaturalConvection ---

    #[test]
    fn test_natural_convection_rayleigh() {
        let nc = NaturalConvection::new(9.81e-4, 1e-3, 310.0, 290.0, 10.0, 1e-3, 1e-4);
        assert!(nc.rayleigh > 0.0);
    }

    #[test]
    fn test_natural_convection_prandtl() {
        let nc = NaturalConvection::new(9.81e-4, 1e-3, 310.0, 290.0, 10.0, 6e-4, 1e-4);
        assert!((nc.prandtl - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_boussinesq_force_at_ref() {
        let force = boussinesq_force(9.81e-4, 1e-3, 300.0, 300.0);
        assert!((force).abs() < 1e-20);
    }

    #[test]
    fn test_boussinesq_force_hot() {
        let force = boussinesq_force(9.81, 1e-3, 310.0, 300.0);
        // -9.81 * 1e-3 * 10 = -0.0981
        assert!((force + 0.0981).abs() < 1e-10);
    }

    // --- ForcedConvection ---

    #[test]
    fn test_forced_convection_reynolds() {
        let fc = ForcedConvection::new(1.0, 10.0, 400.0, 300.0, 1e-3, 1e-4);
        // D_h = 4*H = 40; Re = 1.0 * 40 / 1e-3 = 40000
        let re = fc.reynolds();
        assert!((re - 40000.0).abs() < 1e-6);
    }

    #[test]
    fn test_forced_convection_prandtl() {
        let fc = ForcedConvection::new(1.0, 10.0, 400.0, 300.0, 7e-4, 1e-4);
        assert!((fc.prandtl() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_forced_convection_nusselt_laminar() {
        let fc = ForcedConvection::new(0.1, 5.0, 400.0, 300.0, 1e-3, 1e-4);
        assert!((fc.nusselt_laminar_channel() - 7.54).abs() < 1e-12);
    }

    // --- Helper functions ---

    #[test]
    fn test_nusselt_number_forced() {
        let nu = nusselt_number_forced(100.0, 0.1, 0.5);
        assert!((nu - 100.0 * 0.1 / 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_nusselt_number_natural_positive() {
        let nu = nusselt_number_natural(1e6, 0.71);
        assert!(nu > 0.0);
    }

    #[test]
    fn test_prandtl_number() {
        let pr = prandtl_number(1e-6, 1e-7);
        assert!((pr - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_stanton_number() {
        let st = stanton_number(50.0, 1000.0, 0.7);
        assert!((st - 50.0 / 700.0).abs() < 1e-12);
    }

    #[test]
    fn test_biot_number() {
        let bi = biot_number(100.0, 0.01, 50.0);
        assert!((bi - 100.0 * 0.01 / 50.0).abs() < 1e-12);
    }

    // --- ConjugateHeatSolver ---

    #[test]
    fn test_conjugate_solver_step() {
        let nx = 8;
        let ny = 8;
        let n = nx * ny;
        let is_fluid: Vec<bool> = (0..n)
            .map(|k| {
                let x = k % nx;
                x < 6 // fluid for x < 6, solid for x >= 6
            })
            .collect();
        let solid = SolidRegion::new(16.0, 500.0, 8000.0, "steel");
        let mut solver =
            ConjugateHeatSolver::new(nx, ny, 1.0, 1.0, 1.0, 1.0, is_fluid, solid, 1.0, 1.0)
                .unwrap();
        solver.step();
        assert_eq!(solver.step, 1);
        assert!((solver.time - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_conjugate_solver_results() {
        let nx = 4;
        let ny = 4;
        let n = nx * ny;
        let is_fluid = vec![true; n];
        let solid = SolidRegion::new(1.0, 1000.0, 1000.0, "test");
        let solver =
            ConjugateHeatSolver::new(nx, ny, 1.0, 1.0, 1.0, 1.0, is_fluid, solid, 1.0, 1.0)
                .unwrap();
        let res = solver.results(1.0, 2.0, 0.6);
        assert_eq!(res.temperature.len(), n);
        assert!(res.t_max >= res.t_min);
    }

    #[test]
    fn test_conjugate_result_mean_temperature() {
        let temps = vec![1.0, 2.0, 3.0, 4.0];
        let res = ConjugateResult::from_fields(temps, vec![], 1.0, 1.0);
        assert!((res.mean_temperature() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_d2q5_weights_sum_to_one() {
        let sum: f64 = D2Q5_WEIGHTS.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14);
    }
}
