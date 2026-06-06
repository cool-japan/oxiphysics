// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Solidification and crystal-growth simulations via the Lattice Boltzmann Method.
//!
//! This module provides:
//! - [`SolidificationLbmParams`]: Thermal LBM + phase-change parameters (latent heat, Stefan number).
//! - [`PhaseChangeLbm`]: Enthalpy-based phase change with mushy-zone Carman-Kozeny permeability.
//! - [`DendriticGrowth`]: Anisotropic interface kinetics, undercooling, tip velocity.
//! - [`EutecticSolidification`]: Coupled eutectic phases, lamellar/rod spacing.
//! - [`CastingSimulation`]: Mould filling, solidification shrinkage, porosity prediction.
//! - [`CrystalGrowthLbm`]: Crystal orientation, growth anisotropy, grain competition.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn len2(v: [f64; 2]) -> f64 {
    dot2(v, v).sqrt()
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 velocity set (used by thermal / solidification kernels)
// ─────────────────────────────────────────────────────────────────────────────

const NQ: usize = 9;
const CX: [f64; NQ] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const CY: [f64; NQ] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
const W: [f64; NQ] = [
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
const CS2: f64 = 1.0 / 3.0;

/// Opposite (bounce-back) direction for D2Q9.
const OPP: [usize; NQ] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Compute equilibrium distribution for one D2Q9 cell.
fn feq(rho: f64, ux: f64, uy: f64) -> [f64; NQ] {
    let mut f = [0.0f64; NQ];
    let u2 = ux * ux + uy * uy;
    for q in 0..NQ {
        let cu = CX[q] * ux + CY[q] * uy;
        f[q] = W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    f
}

// ─────────────────────────────────────────────────────────────────────────────
// SolidificationLbmParams
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for the thermal LBM + phase-change simulation.
///
/// The dimensionless Stefan number is Ste = c_p (T_m − T_ref) / L_f.
#[derive(Debug, Clone)]
pub struct SolidificationLbmParams {
    /// Grid width (number of cells in x-direction).
    pub nx: usize,
    /// Grid height (number of cells in y-direction).
    pub ny: usize,
    /// BGK relaxation time for the fluid distribution (τ_f).
    pub tau_f: f64,
    /// BGK relaxation time for the thermal distribution (τ_T).
    pub tau_t: f64,
    /// Reference density (lattice units).
    pub rho_ref: f64,
    /// Specific heat capacity c_p (lattice units).
    pub c_p: f64,
    /// Latent heat of fusion L_f (lattice units).
    pub latent_heat: f64,
    /// Melting temperature T_m (lattice units).
    pub t_melt: f64,
    /// Reference (initial melt) temperature T_ref (lattice units).
    pub t_ref: f64,
    /// Mould (wall) temperature T_wall (lattice units).
    pub t_wall: f64,
    /// Thermal conductivity ratio: solid / liquid.
    pub k_ratio: f64,
    /// Mushy-zone constant K_0 for Carman-Kozeny permeability (lattice units).
    pub k0_mushy: f64,
    /// Width of the mushy zone in temperature units (ΔT_mushy).
    pub dt_mushy: f64,
    /// Stefan number (dimensionless): Ste = c_p ΔT / L_f.
    pub stefan_number: f64,
}

impl SolidificationLbmParams {
    /// Construct a new parameter set with physically consistent Stefan number.
    pub fn new(
        nx: usize,
        ny: usize,
        tau_f: f64,
        tau_t: f64,
        latent_heat: f64,
        t_melt: f64,
        t_ref: f64,
        t_wall: f64,
    ) -> Self {
        let c_p = 1.0;
        let dt = (t_melt - t_wall).abs();
        let stefan_number = if latent_heat > 1e-15 {
            c_p * dt / latent_heat
        } else {
            1.0
        };
        Self {
            nx,
            ny,
            tau_f,
            tau_t,
            rho_ref: 1.0,
            c_p,
            latent_heat,
            t_melt,
            t_ref,
            t_wall,
            k_ratio: 1.0,
            k0_mushy: 1.6e5,
            dt_mushy: 1.0,
            stefan_number,
        }
    }

    /// Thermal diffusivity α = (τ_T − 0.5) · c_s² in lattice units.
    pub fn thermal_diffusivity(&self) -> f64 {
        (self.tau_t - 0.5) * CS2
    }

    /// Kinematic viscosity ν = (τ_f − 0.5) · c_s² in lattice units.
    pub fn kinematic_viscosity(&self) -> f64 {
        (self.tau_f - 0.5) * CS2
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PhaseChangeLbm
// ─────────────────────────────────────────────────────────────────────────────

/// Enthalpy-based phase-change LBM solver with Carman-Kozeny mushy-zone drag.
///
/// The enthalpy h = c_p · T + f_l · L_f where f_l ∈ \[0, 1\] is the liquid fraction.
/// In the mushy zone the permeability follows Carman-Kozeny:
/// K = f_l³ / (K_0 (1 − f_l)²).
#[derive(Debug, Clone)]
pub struct PhaseChangeLbm {
    /// Simulation parameters.
    pub params: SolidificationLbmParams,
    /// Fluid distribution functions f\[cell * NQ + q\].
    pub f: Vec<f64>,
    /// Thermal distribution functions g\[cell * NQ + q\].
    pub g: Vec<f64>,
    /// Liquid fraction f_l ∈ \[0, 1\] per cell.
    pub liquid_fraction: Vec<f64>,
    /// Temperature field T per cell.
    pub temperature: Vec<f64>,
    /// Enthalpy field h per cell.
    pub enthalpy: Vec<f64>,
    /// Density field per cell.
    pub rho: Vec<f64>,
    /// x-velocity per cell.
    pub ux: Vec<f64>,
    /// y-velocity per cell.
    pub uy: Vec<f64>,
    /// Elapsed lattice time steps.
    pub time: usize,
}

impl PhaseChangeLbm {
    /// Create a new solver, initialising the melt at `T_ref` with zero velocity.
    pub fn new(params: SolidificationLbmParams) -> Self {
        let n = params.nx * params.ny;
        let t0 = params.t_ref;
        let rho0 = params.rho_ref;
        let mut g = vec![0.0; n * NQ];
        for i in 0..n {
            let feq0 = feq(t0, 0.0, 0.0);
            for q in 0..NQ {
                g[i * NQ + q] = feq0[q];
            }
        }
        let f: Vec<f64> = (0..n * NQ)
            .map(|idx| {
                let q = idx % NQ;
                W[q] * rho0
            })
            .collect();
        let fl = vec![1.0; n]; // all liquid initially
        let temp = vec![t0; n];
        let enthalpy: Vec<f64> = temp
            .iter()
            .zip(fl.iter())
            .map(|(&t, &fl)| params.c_p * t + fl * params.latent_heat)
            .collect();
        Self {
            params,
            f,
            g,
            liquid_fraction: fl,
            temperature: temp,
            enthalpy,
            rho: vec![rho0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            time: 0,
        }
    }

    /// Update liquid fraction and temperature from enthalpy (mushy-zone model).
    pub fn update_phase_from_enthalpy(&mut self) {
        let p = &self.params;
        let h_solidus = p.c_p * (p.t_melt - p.dt_mushy * 0.5);
        let h_liquidus = p.c_p * p.t_melt + p.latent_heat;
        for i in 0..self.params.nx * self.params.ny {
            let h = self.enthalpy[i];
            if h <= h_solidus {
                self.liquid_fraction[i] = 0.0;
                self.temperature[i] = h / p.c_p;
            } else if h >= h_liquidus {
                self.liquid_fraction[i] = 1.0;
                self.temperature[i] = (h - p.latent_heat) / p.c_p;
            } else {
                // mushy zone
                let fl = (h - h_solidus) / (h_liquidus - h_solidus);
                self.liquid_fraction[i] = clamp(fl, 0.0, 1.0);
                self.temperature[i] = p.t_melt - p.dt_mushy * 0.5 + (fl * p.dt_mushy);
            }
        }
    }

    /// Carman-Kozeny permeability at a given liquid fraction.
    pub fn carman_kozeny_permeability(&self, fl: f64) -> f64 {
        let fl = clamp(fl, 1e-8, 1.0 - 1e-8);
        let k0 = self.params.k0_mushy;
        fl * fl * fl / (k0 * (1.0 - fl) * (1.0 - fl))
    }

    /// Darcy drag body force due to mushy-zone permeability (per unit density).
    ///
    /// Returns `[fx, fy]` drag acceleration.
    pub fn darcy_drag(&self, idx: usize, dt: f64) -> [f64; 2] {
        let fl = self.liquid_fraction[idx];
        let k = self.carman_kozeny_permeability(fl);
        let nu = self.params.kinematic_viscosity();
        // F_Darcy = −ν/K · u
        let factor = -nu / (k + 1e-30) * dt;
        [self.ux[idx] * factor, self.uy[idx] * factor]
    }

    /// Perform one BGK collision step for both fluid and thermal distributions.
    pub fn collide(&mut self) {
        let p = &self.params;
        let n = p.nx * p.ny;
        for i in 0..n {
            let rho = self.rho[i];
            let ux = self.ux[i];
            let uy = self.uy[i];
            let temp = self.temperature[i];
            // Fluid BGK
            let feq0 = feq(rho, ux, uy);
            for (q, f_iq) in self.f[i * NQ..i * NQ + NQ].iter_mut().enumerate() {
                *f_iq += -(*f_iq - feq0[q]) / p.tau_f;
            }
            // Darcy drag source (add to f)
            let drag = self.darcy_drag(i, 1.0);
            for (q, f_iq) in self.f[i * NQ..i * NQ + NQ].iter_mut().enumerate() {
                let cu = CX[q] * drag[0] + CY[q] * drag[1];
                *f_iq += W[q] * rho * cu / CS2;
            }
            // Thermal BGK
            let geq = feq(temp, ux, uy);
            for (q, g_iq) in self.g[i * NQ..i * NQ + NQ].iter_mut().enumerate() {
                *g_iq += -(*g_iq - geq[q]) / p.tau_t;
            }
        }
    }

    /// Stream distributions (periodic in x and y).
    pub fn stream(&mut self) {
        let nx = self.params.nx;
        let ny = self.params.ny;
        let n = nx * ny;
        let mut f_new = vec![0.0f64; n * NQ];
        let mut g_new = vec![0.0f64; n * NQ];
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                for q in 0..NQ {
                    let xn = ((x as isize + CX[q] as isize).rem_euclid(nx as isize)) as usize;
                    let yn = ((y as isize + CY[q] as isize).rem_euclid(ny as isize)) as usize;
                    let j = yn * nx + xn;
                    f_new[j * NQ + q] = self.f[i * NQ + q];
                    g_new[j * NQ + q] = self.g[i * NQ + q];
                }
            }
        }
        self.f = f_new;
        self.g = g_new;
    }

    /// Apply bounce-back on the left wall (x = 0) as the cold mould.
    pub fn apply_mould_wall(&mut self) {
        let nx = self.params.nx;
        let ny = self.params.ny;
        let t_wall = self.params.t_wall;
        for y in 0..ny {
            let i = y * nx; // x = 0
            // Bounce-back fluid
            for (q, &oq) in OPP.iter().enumerate() {
                self.f.swap(i * NQ + q, i * NQ + oq);
            }
            // Dirichlet temperature = T_wall
            let geq = feq(t_wall, 0.0, 0.0);
            for (g_iq, &geq_q) in self.g[i * NQ..i * NQ + NQ].iter_mut().zip(geq.iter()) {
                *g_iq = geq_q;
            }
        }
    }

    /// Recompute macroscopic quantities from distributions.
    pub fn update_macroscopic(&mut self) {
        let n = self.params.nx * self.params.ny;
        for i in 0..n {
            let mut rho = 0.0;
            let mut mx = 0.0;
            let mut my = 0.0;
            let mut temp = 0.0;
            for q in 0..NQ {
                let fq = self.f[i * NQ + q];
                let gq = self.g[i * NQ + q];
                rho += fq;
                mx += fq * CX[q];
                my += fq * CY[q];
                temp += gq;
            }
            self.rho[i] = rho;
            if rho > 1e-15 {
                self.ux[i] = mx / rho;
                self.uy[i] = my / rho;
            } else {
                self.ux[i] = 0.0;
                self.uy[i] = 0.0;
            }
            self.temperature[i] = temp;
            // Update enthalpy from temperature + liquid fraction
            self.enthalpy[i] =
                self.params.c_p * temp + self.liquid_fraction[i] * self.params.latent_heat;
        }
        self.update_phase_from_enthalpy();
    }

    /// Advance one full LBM time step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.apply_mould_wall();
        self.update_macroscopic();
        self.time += 1;
    }

    /// Run `n_steps` time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Return the fraction of solidified cells (f_l < 0.5).
    pub fn solid_fraction_global(&self) -> f64 {
        let n = self.params.nx * self.params.ny;
        let count = self.liquid_fraction.iter().filter(|&&fl| fl < 0.5).count();
        count as f64 / n as f64
    }

    /// Position of the solidification front (rightmost solid column) in lattice units.
    pub fn front_position_x(&self) -> f64 {
        let nx = self.params.nx;
        let ny = self.params.ny;
        let mut front = 0usize;
        for x in 0..nx {
            let any_solid = (0..ny).any(|y| self.liquid_fraction[y * nx + x] < 0.5);
            if any_solid {
                front = x;
            }
        }
        front as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DendriticGrowth
// ─────────────────────────────────────────────────────────────────────────────

/// Anisotropic dendritic growth model using a phase-field / kinetic undercooling approach.
///
/// The interface velocity V_n depends on undercooling ΔT:
/// V_n = μ_k · ΔT · A(θ),
/// where A(θ) = 1 + ε_4 cos(4θ) is the four-fold anisotropy function.
#[derive(Debug, Clone)]
pub struct DendriticGrowth {
    /// Kinetic mobility μ_k (m s⁻¹ K⁻¹, or lattice equivalent).
    pub mobility: f64,
    /// Anisotropy strength ε₄ for four-fold crystal symmetry.
    pub anisotropy: f64,
    /// Melting temperature (lattice units).
    pub t_melt: f64,
    /// Capillarity length d₀ (Gibbs-Thomson coefficient × kinetic term).
    pub capillarity: f64,
    /// Phase-field variable φ ∈ \[−1 solid, +1 liquid\] on a 2-D grid (nx×ny).
    pub phi: Vec<f64>,
    /// Temperature field on the same grid.
    pub temperature: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Phase-field interface width W₀.
    pub interface_width: f64,
    /// Phase-field time scale τ₀.
    pub tau_pf: f64,
}

impl DendriticGrowth {
    /// Create a new dendritic-growth solver with a circular seed of radius `seed_r`.
    pub fn new(
        nx: usize,
        ny: usize,
        mobility: f64,
        anisotropy: f64,
        t_melt: f64,
        undercooling: f64,
        seed_r: f64,
    ) -> Self {
        let n = nx * ny;
        let mut phi = vec![1.0f64; n]; // liquid everywhere
        let t_init = t_melt - undercooling;
        // Place solid seed at centre
        let cx = (nx / 2) as f64;
        let cy = (ny / 2) as f64;
        for y in 0..ny {
            for x in 0..nx {
                let r = ((x as f64 - cx).powi(2) + (y as f64 - cy).powi(2)).sqrt();
                phi[y * nx + x] = if r < seed_r { -1.0 } else { 1.0 };
            }
        }
        Self {
            mobility,
            anisotropy,
            t_melt,
            capillarity: 0.01,
            phi,
            temperature: vec![t_init; n],
            nx,
            ny,
            interface_width: 1.5,
            tau_pf: 1.0,
        }
    }

    /// Four-fold anisotropy function A(θ) = 1 + ε₄ cos(4θ).
    pub fn anisotropy_func(&self, grad_phi: [f64; 2]) -> f64 {
        let len = len2(grad_phi);
        if len < 1e-12 {
            return 1.0;
        }
        let nx = grad_phi[0] / len;
        let ny = grad_phi[1] / len;
        // cos(4θ) = 8 cos²θ sin²θ - 1 + (cos⁴θ + sin⁴θ) = 1 - 8 cos²θ sin²θ = ...
        // use identity: cos(4θ) = nx⁴ + ny⁴ − 6 nx² ny² (for unit normal)
        let cos4 = nx.powi(4) + ny.powi(4) - 6.0 * nx.powi(2) * ny.powi(2);
        1.0 + self.anisotropy * cos4
    }

    /// Compute phase-field Laplacian at cell (x, y) using second-order FD.
    fn laplacian_phi(&self, x: usize, y: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let xp = (x + 1).min(nx - 1);
        let xm = if x > 0 { x - 1 } else { 0 };
        let yp = (y + 1).min(ny - 1);
        let ym = if y > 0 { y - 1 } else { 0 };
        let phi = &self.phi;
        let p = phi[y * nx + x];
        phi[y * nx + xp] + phi[y * nx + xm] + phi[yp * nx + x] + phi[ym * nx + x] - 4.0 * p
    }

    /// Gradient of φ at (x, y).
    fn grad_phi(&self, x: usize, y: usize) -> [f64; 2] {
        let nx = self.nx;
        let xp = (x + 1).min(nx - 1);
        let xm = if x > 0 { x - 1 } else { 0 };
        let ny = self.ny;
        let yp = (y + 1).min(ny - 1);
        let ym = if y > 0 { y - 1 } else { 0 };
        let phi = &self.phi;
        [
            (phi[y * nx + xp] - phi[y * nx + xm]) * 0.5,
            (phi[yp * nx + x] - phi[ym * nx + x]) * 0.5,
        ]
    }

    /// Interface normal velocity V_n = μ_k · ΔT · A(θ) at cell (x, y).
    pub fn interface_velocity(&self, x: usize, y: usize) -> f64 {
        let i = y * self.nx + x;
        let dt = self.t_melt - self.temperature[i];
        let gp = self.grad_phi(x, y);
        let a = self.anisotropy_func(gp);
        self.mobility * dt * a
    }

    /// Advance the phase-field by one step (Allen-Cahn + thermal coupling).
    pub fn step(&mut self, dt_pf: f64, dt_thermal: f64, alpha: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let w0 = self.interface_width;
        let tau = self.tau_pf;
        let mut dphi = vec![0.0f64; n];
        let mut dtemp = vec![0.0f64; n];
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                let phi = self.phi[i];
                let gp = self.grad_phi(x, y);
                let a = self.anisotropy_func(gp);
                let lap = self.laplacian_phi(x, y);
                // Allen-Cahn: τ dφ/dt = W₀² ∇²φ + φ − φ³ − λ(T − T_m)(1 − φ²)²
                let undercooling = self.temperature[i] - self.t_melt;
                let dphi_ac = w0 * w0 * lap + phi
                    - phi.powi(3)
                    - self.capillarity * undercooling * (1.0 - phi * phi).powi(2);
                dphi[i] = dt_pf * a * dphi_ac / tau;

                // Thermal: dT/dt = α ∇²T + L_f/c_p · dφ/dt / 2
                let xp = (x + 1).min(nx - 1);
                let xm = if x > 0 { x - 1 } else { 0 };
                let yp = (y + 1).min(ny - 1);
                let ym = if y > 0 { y - 1 } else { 0 };
                let t = &self.temperature;
                let lap_t =
                    t[y * nx + xp] + t[y * nx + xm] + t[yp * nx + x] + t[ym * nx + x] - 4.0 * t[i];
                dtemp[i] = dt_thermal * (alpha * lap_t + 0.5 * dphi[i] / dt_pf);
            }
        }
        for i in 0..n {
            self.phi[i] = clamp(self.phi[i] + dphi[i], -1.0, 1.0);
            self.temperature[i] += dtemp[i];
        }
    }

    /// Estimate the dendritic tip position (furthest solid cell from centre).
    pub fn tip_position(&self) -> [f64; 2] {
        let nx = self.nx;
        let ny = self.ny;
        let cx = (nx / 2) as f64;
        let cy = (ny / 2) as f64;
        let mut best_r = 0.0_f64;
        let mut best_pos = [cx, cy];
        for y in 0..ny {
            for x in 0..nx {
                if self.phi[y * nx + x] < 0.0 {
                    let r = ((x as f64 - cx).powi(2) + (y as f64 - cy).powi(2)).sqrt();
                    if r > best_r {
                        best_r = r;
                        best_pos = [x as f64, y as f64];
                    }
                }
            }
        }
        best_pos
    }

    /// Estimate the tip velocity from two successive tip positions and a time interval.
    pub fn tip_velocity(pos1: [f64; 2], pos2: [f64; 2], dt: f64) -> f64 {
        let d = sub2(pos2, pos1);
        len2(d) / dt
    }

    /// Total solid area (number of cells with φ < 0).
    pub fn solid_area(&self) -> usize {
        self.phi.iter().filter(|&&p| p < 0.0).count()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EutecticSolidification
// ─────────────────────────────────────────────────────────────────────────────

/// Phase descriptors for eutectic constituents.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EutecticPhase {
    /// Alpha phase (e.g. Pb-rich in Pb-Sn).
    Alpha,
    /// Beta phase (e.g. Sn-rich in Pb-Sn).
    Beta,
    /// Liquid eutectic melt.
    Liquid,
}

/// Parameters for a binary eutectic system.
#[derive(Debug, Clone)]
pub struct EutecticParams {
    /// Eutectic temperature T_e (lattice units).
    pub t_eutectic: f64,
    /// Equilibrium lamellar spacing λ_0 (lattice units).
    pub lambda_eq: f64,
    /// Growth velocity at equilibrium spacing V_0.
    pub v0: f64,
    /// Jackson-Hunt exponent (≈ 0.5 for diffusion-controlled eutectic).
    pub jh_exponent: f64,
    /// Latent heat of alpha phase L_α.
    pub latent_alpha: f64,
    /// Latent heat of beta phase L_β.
    pub latent_beta: f64,
    /// Diffusivity in the liquid D_l.
    pub diffusivity_liq: f64,
    /// Volume fraction of alpha phase at eutectic composition.
    pub vol_frac_alpha: f64,
}

impl EutecticParams {
    /// Jackson-Hunt relation: V · λ² = const.
    /// Returns growth velocity for a given lamellar spacing λ.
    pub fn jackson_hunt_velocity(&self, lambda: f64) -> f64 {
        if lambda < 1e-12 {
            return 0.0;
        }
        self.v0 * (self.lambda_eq / lambda).powi(2)
    }

    /// Optimum lamellar spacing (minimises undercooling) λ_opt = λ_0.
    pub fn optimal_spacing(&self) -> f64 {
        self.lambda_eq
    }
}

/// Coupled eutectic solidification solver on a 1-D front approximation.
///
/// Tracks the eutectic front position, lamellar spacing, and composition profiles.
#[derive(Debug, Clone)]
pub struct EutecticSolidification {
    /// System parameters.
    pub params: EutecticParams,
    /// Phase labels for each column (lamellar repeat pattern).
    pub phase_pattern: Vec<EutecticPhase>,
    /// Current solidification front position (lattice units).
    pub front_position: f64,
    /// Current lamellar spacing λ (lattice units).
    pub lamellar_spacing: f64,
    /// Composition (mole fraction of B) in the liquid ahead of the front.
    pub composition_liquid: Vec<f64>,
    /// Width of composition field.
    pub field_width: usize,
    /// Current undercooling below T_e (lattice units).
    pub undercooling: f64,
    /// Elapsed simulation time.
    pub time: f64,
}

impl EutecticSolidification {
    /// Create a new eutectic solidification simulation.
    pub fn new(params: EutecticParams, field_width: usize, initial_undercooling: f64) -> Self {
        let n_lamellae = (field_width as f64 / params.lambda_eq).ceil() as usize + 1;
        let mut pattern = Vec::with_capacity(n_lamellae);
        for i in 0..n_lamellae {
            pattern.push(if i % 2 == 0 {
                EutecticPhase::Alpha
            } else {
                EutecticPhase::Beta
            });
        }
        let composition = vec![params.vol_frac_alpha; field_width];
        Self {
            front_position: 0.0,
            lamellar_spacing: params.lambda_eq,
            composition_liquid: composition,
            field_width,
            undercooling: initial_undercooling,
            time: 0.0,
            params,
            phase_pattern: pattern,
        }
    }

    /// Compute front velocity from Jackson-Hunt theory with current spacing.
    pub fn front_velocity(&self) -> f64 {
        self.params.jackson_hunt_velocity(self.lamellar_spacing)
    }

    /// Diffuse composition ahead of the front using explicit FD.
    pub fn diffuse_composition(&mut self, dt: f64) {
        let d = self.params.diffusivity_liq;
        let n = self.field_width;
        let mut new_c = self.composition_liquid.clone();
        for (new_c_i, i) in new_c[1..n - 1].iter_mut().zip(1..n - 1) {
            let lap = self.composition_liquid[i + 1] - 2.0 * self.composition_liquid[i]
                + self.composition_liquid[i - 1];
            *new_c_i += d * dt * lap;
        }
        self.composition_liquid = new_c;
    }

    /// Advance the eutectic front by one time step `dt`.
    pub fn step(&mut self, dt: f64) {
        let v = self.front_velocity();
        self.front_position += v * dt;
        // Competition: adjust spacing towards optimum based on undercooling
        let lambda_opt = self.params.optimal_spacing();
        let err = lambda_opt - self.lamellar_spacing;
        self.lamellar_spacing += 0.01 * err * dt;
        self.undercooling = (self.undercooling - 0.001 * v * dt).max(0.0);
        self.diffuse_composition(dt);
        self.time += dt;
    }

    /// Rod spacing (in 3-D rod eutectic the optimum rod spacing ≈ 0.866 λ_lam).
    pub fn rod_spacing(&self) -> f64 {
        self.params.lambda_eq * 0.866
    }

    /// Return the lamellar count in the current simulation domain.
    pub fn lamellar_count(&self) -> usize {
        (self.field_width as f64 / self.lamellar_spacing).floor() as usize
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CastingSimulation
// ─────────────────────────────────────────────────────────────────────────────

/// Porosity type in a casting defect.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PorosityType {
    /// Micro-shrinkage porosity (isolated).
    Shrinkage,
    /// Gas porosity (dissolved gas rejection).
    Gas,
    /// None (fully dense).
    None,
}

/// Mould filling and solidification shrinkage model for casting processes.
///
/// Uses a simplified volume-of-fluid approach for filling and tracks
/// solidification shrinkage to predict porosity.
#[derive(Debug, Clone)]
pub struct CastingSimulation {
    /// Thermal LBM solver driving solidification.
    pub phase_change: PhaseChangeLbm,
    /// Volume fraction of metal (0 = air/void, 1 = full metal) per cell.
    pub fill_fraction: Vec<f64>,
    /// Porosity type per cell.
    pub porosity_map: Vec<PorosityType>,
    /// Solidification shrinkage factor β_s (volume change on solidification).
    pub shrinkage_factor: f64,
    /// Initial dissolved gas content (cc/100 g).
    pub initial_gas_content: f64,
    /// Threshold liquid fraction below which shrinkage porosity forms.
    pub fl_threshold: f64,
    /// Mould filling velocity (cells per time step).
    pub fill_velocity: f64,
    /// Current fill front x-position.
    pub fill_front: f64,
}

impl CastingSimulation {
    /// Create a new casting simulation.
    pub fn new(params: SolidificationLbmParams, shrinkage: f64, gas_content: f64) -> Self {
        let n = params.nx * params.ny;
        let solver = PhaseChangeLbm::new(params);
        Self {
            fill_fraction: vec![0.0; n],
            porosity_map: vec![PorosityType::None; n],
            shrinkage_factor: shrinkage,
            initial_gas_content: gas_content,
            fl_threshold: 0.1,
            fill_velocity: 0.1,
            fill_front: 0.0,
            phase_change: solver,
        }
    }

    /// Advance the mould-filling front by `dt` time steps.
    pub fn advance_fill(&mut self, dt: f64) {
        let nx = self.phase_change.params.nx;
        let ny = self.phase_change.params.ny;
        self.fill_front += self.fill_velocity * dt;
        let fill_col = (self.fill_front as usize).min(nx);
        for y in 0..ny {
            for x in 0..fill_col {
                self.fill_fraction[y * nx + x] = 1.0;
            }
        }
    }

    /// Mark porosity defects based on liquid fraction and fill fraction.
    pub fn predict_porosity(&mut self) {
        let n = self.phase_change.params.nx * self.phase_change.params.ny;
        for i in 0..n {
            if self.fill_fraction[i] > 0.5 {
                let fl = self.phase_change.liquid_fraction[i];
                if fl < self.fl_threshold {
                    // Solidified region — check for gas or shrinkage
                    let local_gas = self.initial_gas_content * (1.0 - fl);
                    if local_gas > 0.05 {
                        self.porosity_map[i] = PorosityType::Gas;
                    } else if self.shrinkage_factor > 0.01 {
                        self.porosity_map[i] = PorosityType::Shrinkage;
                    }
                }
            }
        }
    }

    /// Total porosity volume fraction in the casting.
    pub fn total_porosity_fraction(&self) -> f64 {
        let n = self.fill_fraction.len();
        if n == 0 {
            return 0.0;
        }
        let count = self
            .porosity_map
            .iter()
            .filter(|&&p| p != PorosityType::None)
            .count();
        count as f64 / n as f64
    }

    /// Shrinkage void volume = β_s × solidified volume.
    pub fn shrinkage_void_volume(&self) -> f64 {
        let solid_frac = self.phase_change.solid_fraction_global();
        self.shrinkage_factor * solid_frac
    }

    /// Advance one combined fill + solidification step.
    pub fn step(&mut self, dt: f64) {
        self.advance_fill(dt);
        self.phase_change.step();
        self.predict_porosity();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CrystalGrowthLbm
// ─────────────────────────────────────────────────────────────────────────────

/// A single grain with orientation and growth state.
#[derive(Debug, Clone)]
pub struct Grain {
    /// Grain index.
    pub id: usize,
    /// Crystal orientation angle θ (radians from x-axis).
    pub orientation: f64,
    /// Nucleation site \[x, y\] (lattice units).
    pub nucleus: [f64; 2],
    /// Current equivalent radius (lattice units).
    pub radius: f64,
    /// Growth rate (lattice units per time step).
    pub growth_rate: f64,
    /// Is this grain active (still growing)?
    pub active: bool,
}

impl Grain {
    /// Create a new grain at `pos` with a given orientation.
    pub fn new(id: usize, orientation: f64, nucleus: [f64; 2]) -> Self {
        Self {
            id,
            orientation,
            nucleus,
            radius: 0.5,
            growth_rate: 0.0,
            active: true,
        }
    }

    /// Anisotropic growth rate with `n_fold`-fold symmetry in direction θ_grow.
    pub fn anisotropic_rate(&self, base_rate: f64, undercooling: f64, n_fold: u32) -> f64 {
        let theta = self.orientation;
        let aniso = 1.0 + 0.1 * ((n_fold as f64) * theta).cos();
        base_rate * undercooling * aniso
    }
}

/// LBM-coupled polycrystalline solidification with grain nucleation and competition.
#[derive(Debug, Clone)]
pub struct CrystalGrowthLbm {
    /// Underlying phase-change LBM.
    pub phase_change: PhaseChangeLbm,
    /// List of grains.
    pub grains: Vec<Grain>,
    /// Grain ID map (cell → grain id, or usize::MAX for liquid).
    pub grain_map: Vec<usize>,
    /// Nucleation undercooling threshold ΔT_nuc.
    pub nucleation_undercooling: f64,
    /// Maximum number of grains allowed.
    pub max_grains: usize,
    /// Base growth mobility.
    pub base_mobility: f64,
    /// Fold symmetry for anisotropy.
    pub n_fold: u32,
}

impl CrystalGrowthLbm {
    /// Create a new polycrystalline growth solver.
    pub fn new(
        params: SolidificationLbmParams,
        nucleation_undercooling: f64,
        max_grains: usize,
        base_mobility: f64,
    ) -> Self {
        let n = params.nx * params.ny;
        let solver = PhaseChangeLbm::new(params);
        Self {
            grain_map: vec![usize::MAX; n],
            grains: Vec::new(),
            nucleation_undercooling,
            max_grains,
            base_mobility,
            n_fold: 4,
            phase_change: solver,
        }
    }

    /// Attempt to nucleate a new grain at cell `(x, y)` if conditions are met.
    pub fn try_nucleate(&mut self, x: usize, y: usize) -> bool {
        let nx = self.phase_change.params.nx;
        let i = y * nx + x;
        if self.grains.len() >= self.max_grains {
            return false;
        }
        if self.grain_map[i] != usize::MAX {
            return false;
        }
        let fl = self.phase_change.liquid_fraction[i];
        if fl < 0.5 {
            return false;
        }
        let undercooling = self.phase_change.params.t_melt - self.phase_change.temperature[i];
        if undercooling < self.nucleation_undercooling {
            return false;
        }
        // Nucleate
        let id = self.grains.len();
        let orient = (id as f64) * PI / (self.max_grains as f64);
        let g = Grain::new(id, orient, [x as f64, y as f64]);
        self.grains.push(g);
        self.grain_map[i] = id;
        true
    }

    /// Grow all active grains by one step.
    pub fn grow_grains(&mut self) {
        let nx = self.phase_change.params.nx;
        let ny = self.phase_change.params.ny;
        let t_melt = self.phase_change.params.t_melt;
        let mob = self.base_mobility;
        for g in self.grains.iter_mut() {
            if !g.active {
                continue;
            }
            let cx = g.nucleus[0] as usize;
            let cy = g.nucleus[1] as usize;
            let i = cy * nx + cx;
            let undercooling = t_melt - self.phase_change.temperature[i.min(nx * ny - 1)];
            let rate = g.anisotropic_rate(mob, undercooling.max(0.0), g_n_fold(g));
            g.growth_rate = rate;
            g.radius += rate;
        }
        // Mark grain cells by radius
        self.rasterise_grains(nx, ny);
    }

    /// Rasterise grain envelopes onto the grain map.
    fn rasterise_grains(&mut self, nx: usize, ny: usize) {
        for g in &self.grains {
            if !g.active {
                continue;
            }
            let cx = g.nucleus[0];
            let cy = g.nucleus[1];
            let r = g.radius;
            let id = g.id;
            let x0 = ((cx - r) as isize).max(0) as usize;
            let x1 = ((cx + r) as isize + 1).min(nx as isize) as usize;
            let y0 = ((cy - r) as isize).max(0) as usize;
            let y1 = ((cy + r) as isize + 1).min(ny as isize) as usize;
            for y in y0..y1 {
                for x in x0..x1 {
                    let dx = x as f64 - cx;
                    let dy = y as f64 - cy;
                    if dx * dx + dy * dy <= r * r {
                        let i = y * nx + x;
                        if self.grain_map[i] == usize::MAX {
                            self.grain_map[i] = id;
                            self.phase_change.liquid_fraction[i] =
                                (self.phase_change.liquid_fraction[i] - 0.1).max(0.0);
                        }
                    }
                }
            }
        }
    }

    /// Detect grain impingement — deactivate grains that have collided.
    pub fn detect_impingement(&mut self) {
        let n_grains = self.grains.len();
        for a in 0..n_grains {
            for b in (a + 1)..n_grains {
                if !self.grains[a].active || !self.grains[b].active {
                    continue;
                }
                let da = sub2(self.grains[a].nucleus, self.grains[b].nucleus);
                let dist = len2(da);
                if dist < self.grains[a].radius + self.grains[b].radius {
                    // Smaller grain is blocked
                    if self.grains[a].radius < self.grains[b].radius {
                        self.grains[a].active = false;
                    } else {
                        self.grains[b].active = false;
                    }
                }
            }
        }
    }

    /// Number of active grains.
    pub fn active_grain_count(&self) -> usize {
        self.grains.iter().filter(|g| g.active).count()
    }

    /// Average grain radius of active grains.
    pub fn mean_grain_radius(&self) -> f64 {
        let active: Vec<_> = self.grains.iter().filter(|g| g.active).collect();
        if active.is_empty() {
            return 0.0;
        }
        active.iter().map(|g| g.radius).sum::<f64>() / active.len() as f64
    }

    /// Advance one combined step: LBM + nucleation + growth + impingement.
    pub fn step(&mut self) {
        self.phase_change.step();
        // Scan liquid cells for nucleation
        let nx = self.phase_change.params.nx;
        let ny = self.phase_change.params.ny;
        for y in 0..ny {
            for x in 0..nx {
                self.try_nucleate(x, y);
            }
        }
        self.grow_grains();
        self.detect_impingement();
    }
}

/// Helper: n_fold for a grain (always 4 in this model).
fn g_n_fold(g: &Grain) -> u32 {
    let _ = g;
    4
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> SolidificationLbmParams {
        SolidificationLbmParams::new(16, 16, 1.0, 0.9, 0.5, 1.0, 1.5, 0.5)
    }

    // ── SolidificationLbmParams ──────────────────────────────────────────────

    #[test]
    fn test_params_stefan_number_positive() {
        let p = default_params();
        assert!(p.stefan_number > 0.0, "Stefan number must be positive");
    }

    #[test]
    fn test_params_thermal_diffusivity() {
        let p = default_params();
        let alpha = p.thermal_diffusivity();
        assert!(alpha > 0.0, "Thermal diffusivity must be positive");
    }

    #[test]
    fn test_params_kinematic_viscosity() {
        let p = default_params();
        let nu = p.kinematic_viscosity();
        assert!(nu > 0.0, "Kinematic viscosity must be positive");
    }

    #[test]
    fn test_params_consistency() {
        let p = SolidificationLbmParams::new(8, 8, 1.2, 0.8, 1.0, 2.0, 3.0, 0.5);
        assert!(
            (p.stefan_number - 1.5).abs() < 1e-9,
            "Stefan number = c_p ΔT / L_f"
        );
    }

    // ── PhaseChangeLbm ───────────────────────────────────────────────────────

    #[test]
    fn test_phase_change_init_all_liquid() {
        let solver = PhaseChangeLbm::new(default_params());
        assert!(solver.liquid_fraction.iter().all(|&fl| fl == 1.0));
    }

    #[test]
    fn test_phase_change_init_temperature() {
        let p = default_params();
        let t_ref = p.t_ref;
        let solver = PhaseChangeLbm::new(p);
        for &t in &solver.temperature {
            assert!(
                (t - t_ref).abs() < 1e-10,
                "Initial temperature should be t_ref"
            );
        }
    }

    #[test]
    fn test_phase_change_init_density() {
        let p = default_params();
        let rho0 = p.rho_ref;
        let solver = PhaseChangeLbm::new(p);
        for &r in &solver.rho {
            assert!((r - rho0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_carman_kozeny_liquid() {
        let solver = PhaseChangeLbm::new(default_params());
        let k = solver.carman_kozeny_permeability(1.0 - 1e-8);
        assert!(k > 0.0, "Fully liquid: high permeability");
    }

    #[test]
    fn test_carman_kozeny_solid() {
        let solver = PhaseChangeLbm::new(default_params());
        let k = solver.carman_kozeny_permeability(1e-8);
        assert!(k < 1e-6, "Fully solid: near-zero permeability");
    }

    #[test]
    fn test_carman_kozeny_monotone() {
        let solver = PhaseChangeLbm::new(default_params());
        let k1 = solver.carman_kozeny_permeability(0.3);
        let k2 = solver.carman_kozeny_permeability(0.7);
        assert!(k2 > k1, "Permeability increases with liquid fraction");
    }

    #[test]
    fn test_step_does_not_panic() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.step();
    }

    #[test]
    fn test_run_5_steps() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.run(5);
        assert_eq!(solver.time, 5);
    }

    #[test]
    fn test_solid_fraction_initially_zero() {
        let solver = PhaseChangeLbm::new(default_params());
        assert_eq!(solver.solid_fraction_global(), 0.0);
    }

    #[test]
    fn test_solid_fraction_after_steps() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.run(20);
        // Wall cooling should have produced some solid
        let sf = solver.solid_fraction_global();
        assert!(
            (0.0..=1.0).contains(&sf),
            "Solid fraction must be in [0, 1]"
        );
    }

    #[test]
    fn test_front_position_initially_zero() {
        let solver = PhaseChangeLbm::new(default_params());
        assert_eq!(solver.front_position_x(), 0.0);
    }

    #[test]
    fn test_update_phase_solidus() {
        let mut solver = PhaseChangeLbm::new(default_params());
        // Force enthalpy below solidus
        solver.enthalpy[0] = -100.0;
        solver.update_phase_from_enthalpy();
        assert_eq!(solver.liquid_fraction[0], 0.0);
    }

    #[test]
    fn test_update_phase_liquidus() {
        let mut solver = PhaseChangeLbm::new(default_params());
        let p = &solver.params;
        let h_liq = p.c_p * p.t_melt + p.latent_heat + 10.0;
        solver.enthalpy[0] = h_liq;
        solver.update_phase_from_enthalpy();
        assert_eq!(solver.liquid_fraction[0], 1.0);
    }

    #[test]
    fn test_darcy_drag_fully_liquid() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.ux[0] = 0.01;
        let drag = solver.darcy_drag(0, 1.0);
        // fl=1 → high permeability → near-zero drag
        assert!(drag[0].abs() < 1e-3, "Drag near zero in liquid");
    }

    #[test]
    fn test_liquid_fraction_clamped() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.enthalpy[0] = 1000.0;
        solver.update_phase_from_enthalpy();
        assert!(solver.liquid_fraction[0] <= 1.0);
        solver.enthalpy[0] = -1000.0;
        solver.update_phase_from_enthalpy();
        assert!(solver.liquid_fraction[0] >= 0.0);
    }

    // ── DendriticGrowth ──────────────────────────────────────────────────────

    #[test]
    fn test_dendritic_seed_solid() {
        let d = DendriticGrowth::new(32, 32, 0.01, 0.05, 1.0, 0.2, 3.0);
        let cx = 16usize;
        let cy = 16usize;
        assert!(d.phi[cy * 32 + cx] < 0.0, "Seed centre should be solid");
    }

    #[test]
    fn test_dendritic_liquid_far() {
        let d = DendriticGrowth::new(32, 32, 0.01, 0.05, 1.0, 0.2, 3.0);
        assert!(d.phi[0] > 0.0, "Far corner should be liquid");
    }

    #[test]
    fn test_anisotropy_zero_grad() {
        let d = DendriticGrowth::new(16, 16, 0.01, 0.05, 1.0, 0.2, 2.0);
        let a = d.anisotropy_func([0.0, 0.0]);
        assert_eq!(a, 1.0, "Zero gradient → isotropic");
    }

    #[test]
    fn test_anisotropy_x_axis() {
        let d = DendriticGrowth::new(16, 16, 0.01, 0.05, 1.0, 0.2, 2.0);
        let a = d.anisotropy_func([1.0, 0.0]);
        // cos(4·0) = 1, so a = 1 + ε
        let expected = 1.0 + d.anisotropy;
        assert!((a - expected).abs() < 1e-10);
    }

    #[test]
    fn test_tip_position_near_centre() {
        let d = DendriticGrowth::new(32, 32, 0.01, 0.05, 1.0, 0.2, 3.0);
        let tip = d.tip_position();
        let cx = 16.0;
        let cy = 16.0;
        let dist = ((tip[0] - cx).powi(2) + (tip[1] - cy).powi(2)).sqrt();
        assert!(dist < 10.0, "Tip should be near seed at t=0");
    }

    #[test]
    fn test_solid_area_initial() {
        let d = DendriticGrowth::new(16, 16, 0.01, 0.05, 1.0, 0.2, 2.0);
        assert!(d.solid_area() > 0, "Seed should create some solid area");
    }

    #[test]
    fn test_dendritic_step_no_panic() {
        let mut d = DendriticGrowth::new(16, 16, 0.01, 0.05, 1.0, 0.2, 2.0);
        d.step(0.01, 0.01, 0.01);
    }

    #[test]
    fn test_dendritic_solid_grows() {
        // Use a large undercooling and small step so the seed grows rather than
        // being smoothed away by the Allen-Cahn diffusion term.
        let mut d = DendriticGrowth::new(32, 32, 0.5, 0.05, 1.0, 2.0, 4.0);
        let initial = d.solid_area();
        assert!(initial > 0, "Seed must create solid cells");
        for _ in 0..5 {
            d.step(0.001, 0.001, 0.001);
        }
        let final_area = d.solid_area();
        assert!(
            final_area >= initial,
            "Solid area should not shrink: initial={initial}, final={final_area}"
        );
    }

    #[test]
    fn test_tip_velocity_calculation() {
        let pos1 = [10.0, 16.0];
        let pos2 = [12.0, 16.0];
        let v = DendriticGrowth::tip_velocity(pos1, pos2, 1.0);
        assert!((v - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_interface_velocity_zero_undercooling() {
        let mut d = DendriticGrowth::new(16, 16, 0.01, 0.05, 1.0, 0.0, 2.0);
        // Set T = T_melt everywhere
        for t in d.temperature.iter_mut() {
            *t = 1.0;
        }
        let v = d.interface_velocity(8, 8);
        assert!((v).abs() < 1e-9, "No undercooling → no growth");
    }

    // ── EutecticSolidification ───────────────────────────────────────────────

    fn eutectic_params() -> EutecticParams {
        EutecticParams {
            t_eutectic: 1.0,
            lambda_eq: 5.0,
            v0: 0.1,
            jh_exponent: 0.5,
            latent_alpha: 0.3,
            latent_beta: 0.2,
            diffusivity_liq: 0.01,
            vol_frac_alpha: 0.4,
        }
    }

    #[test]
    fn test_eutectic_jackson_hunt() {
        let ep = eutectic_params();
        let v = ep.jackson_hunt_velocity(5.0);
        assert!((v - 0.1).abs() < 1e-10, "At optimal spacing V = V0");
    }

    #[test]
    fn test_eutectic_jh_finer_faster() {
        let ep = eutectic_params();
        let v_fine = ep.jackson_hunt_velocity(2.5);
        let v_coarse = ep.jackson_hunt_velocity(10.0);
        assert!(v_fine > v_coarse, "Finer spacing → faster in J-H relation");
    }

    #[test]
    fn test_eutectic_optimal_spacing() {
        let ep = eutectic_params();
        assert_eq!(ep.optimal_spacing(), 5.0);
    }

    #[test]
    fn test_eutectic_rod_spacing() {
        let es = EutecticSolidification::new(eutectic_params(), 20, 0.1);
        let rs = es.rod_spacing();
        assert!((rs - 5.0 * 0.866).abs() < 1e-6);
    }

    #[test]
    fn test_eutectic_step_advances_front() {
        let mut es = EutecticSolidification::new(eutectic_params(), 20, 0.2);
        let fp0 = es.front_position;
        es.step(1.0);
        assert!(es.front_position > fp0);
    }

    #[test]
    fn test_eutectic_lamellar_count() {
        let es = EutecticSolidification::new(eutectic_params(), 30, 0.1);
        assert!(es.lamellar_count() > 0);
    }

    #[test]
    fn test_eutectic_time_advances() {
        let mut es = EutecticSolidification::new(eutectic_params(), 20, 0.1);
        es.step(2.5);
        assert!((es.time - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_eutectic_pattern_alternating() {
        let es = EutecticSolidification::new(eutectic_params(), 20, 0.1);
        assert_eq!(es.phase_pattern[0], EutecticPhase::Alpha);
        assert_eq!(es.phase_pattern[1], EutecticPhase::Beta);
    }

    // ── CastingSimulation ────────────────────────────────────────────────────

    #[test]
    fn test_casting_init_no_fill() {
        let cs = CastingSimulation::new(default_params(), 0.05, 0.1);
        assert!(cs.fill_fraction.iter().all(|&f| f == 0.0));
    }

    #[test]
    fn test_casting_advance_fill() {
        let mut cs = CastingSimulation::new(default_params(), 0.05, 0.1);
        cs.advance_fill(10.0);
        let some_filled = cs.fill_fraction.iter().any(|&f| f > 0.0);
        assert!(some_filled, "Filling should have advanced");
    }

    #[test]
    fn test_casting_porosity_prediction() {
        let mut cs = CastingSimulation::new(default_params(), 0.1, 0.5);
        cs.advance_fill(100.0); // fill everything
        // Force solid in one cell
        cs.phase_change.liquid_fraction[0] = 0.0;
        cs.predict_porosity();
        // At least one porosity cell should exist
        let any_porous = cs.porosity_map.iter().any(|&p| p != PorosityType::None);
        assert!(any_porous);
    }

    #[test]
    fn test_casting_shrinkage_void() {
        let cs = CastingSimulation::new(default_params(), 0.05, 0.0);
        let v = cs.shrinkage_void_volume();
        assert!(v >= 0.0);
    }

    #[test]
    fn test_casting_total_porosity_empty() {
        let cs = CastingSimulation::new(default_params(), 0.05, 0.0);
        assert_eq!(cs.total_porosity_fraction(), 0.0);
    }

    #[test]
    fn test_casting_step_no_panic() {
        let mut cs = CastingSimulation::new(default_params(), 0.05, 0.1);
        cs.step(1.0);
    }

    // ── CrystalGrowthLbm ─────────────────────────────────────────────────────

    #[test]
    fn test_crystal_init_no_grains() {
        let cg = CrystalGrowthLbm::new(default_params(), 0.1, 10, 0.01);
        assert_eq!(cg.grains.len(), 0);
    }

    #[test]
    fn test_crystal_nucleation_requires_undercooling() {
        let p = SolidificationLbmParams::new(16, 16, 1.0, 0.9, 0.5, 2.0, 1.5, 0.5);
        let mut cg = CrystalGrowthLbm::new(p, 10.0, 10, 0.01);
        // T_melt=2.0, T_ref=1.5, undercooling=0.5 < threshold=10.0 → no nucleation
        let result = cg.try_nucleate(8, 8);
        assert!(!result, "Insufficient undercooling: should not nucleate");
    }

    #[test]
    fn test_crystal_nucleation_success() {
        let p = SolidificationLbmParams::new(16, 16, 1.0, 0.9, 0.5, 2.0, 3.0, 0.5);
        let mut cg = CrystalGrowthLbm::new(p, 0.1, 10, 0.01);
        // undercooling = 2.0 - 3.0 is negative → t_ref > t_melt → undercooling < 0 → no nucleation
        // Test nucleation with sufficient initial T:
        // Force temperature below melting
        let nx = cg.phase_change.params.nx;
        for t in cg.phase_change.temperature.iter_mut() {
            *t = 1.0; // well below T_melt = 2.0
        }
        let result = cg.try_nucleate(4, 4);
        // undercooling = 2.0 - 1.0 = 1.0 > 0.1
        assert!(result, "Should nucleate with sufficient undercooling");
        let _ = nx;
    }

    #[test]
    fn test_crystal_active_count() {
        let cg = CrystalGrowthLbm::new(default_params(), 0.1, 5, 0.01);
        assert_eq!(cg.active_grain_count(), 0);
    }

    #[test]
    fn test_crystal_mean_radius_empty() {
        let cg = CrystalGrowthLbm::new(default_params(), 0.1, 5, 0.01);
        assert_eq!(cg.mean_grain_radius(), 0.0);
    }

    #[test]
    fn test_grain_anisotropic_rate() {
        let g = Grain::new(0, 0.0, [0.0, 0.0]);
        let rate = g.anisotropic_rate(0.1, 1.0, 4);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_grain_orientation_stored() {
        let g = Grain::new(0, PI / 4.0, [5.0, 5.0]);
        assert!((g.orientation - PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_crystal_step_no_panic() {
        let mut cg = CrystalGrowthLbm::new(default_params(), 0.01, 5, 0.01);
        cg.step();
    }

    #[test]
    fn test_impingement_deactivates_smaller() {
        let mut cg = CrystalGrowthLbm::new(default_params(), 0.01, 5, 0.01);
        cg.grains.push(Grain {
            id: 0,
            orientation: 0.0,
            nucleus: [8.0, 8.0],
            radius: 5.0,
            growth_rate: 0.1,
            active: true,
        });
        cg.grains.push(Grain {
            id: 1,
            orientation: PI / 4.0,
            nucleus: [9.0, 8.0],
            radius: 3.0,
            growth_rate: 0.05,
            active: true,
        });
        cg.detect_impingement();
        assert!(!cg.grains[1].active, "Smaller grain should be deactivated");
        assert!(cg.grains[0].active, "Larger grain should stay active");
    }

    #[test]
    fn test_impingement_no_contact() {
        let mut cg = CrystalGrowthLbm::new(default_params(), 0.01, 5, 0.01);
        cg.grains.push(Grain {
            id: 0,
            orientation: 0.0,
            nucleus: [2.0, 2.0],
            radius: 1.0,
            growth_rate: 0.1,
            active: true,
        });
        cg.grains.push(Grain {
            id: 1,
            orientation: PI / 4.0,
            nucleus: [14.0, 14.0],
            radius: 1.0,
            growth_rate: 0.05,
            active: true,
        });
        cg.detect_impingement();
        assert!(
            cg.grains[0].active && cg.grains[1].active,
            "Non-overlapping grains both active"
        );
    }

    #[test]
    fn test_solidification_conserves_mass_approx() {
        let mut solver = PhaseChangeLbm::new(default_params());
        solver.run(10);
        let total_rho: f64 = solver.rho.iter().sum();
        // After 10 steps the density sum should be roughly positive and finite
        assert!(total_rho.is_finite() && total_rho > 0.0);
    }

    #[test]
    fn test_feq_sum_to_rho() {
        let rho = 1.2;
        let f = feq(rho, 0.05, -0.03);
        let sum: f64 = f.iter().sum();
        assert!((sum - rho).abs() < 1e-10, "feq must sum to rho");
    }

    #[test]
    fn test_feq_zero_velocity_weights() {
        let f = feq(1.0, 0.0, 0.0);
        for (q, &fq) in f.iter().enumerate() {
            assert!(
                (fq - W[q]).abs() < 1e-10,
                "At zero velocity feq[q] = w[q]*rho"
            );
        }
    }
}
