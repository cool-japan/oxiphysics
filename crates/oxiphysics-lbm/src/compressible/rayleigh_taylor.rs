// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rayleigh-Taylor instability with a two-population lattice Boltzmann
//! method.
//!
//! * **Flow** — D2Q9 BGK at `cs² = 1/3` with Guo forcing.  The buoyancy force
//!   is `F = −(ρ_mix − ρ_ref) g ŷ` (Boussinesq: the flow population carries a
//!   reference density ≈ 1 and the mixture density enters only through the
//!   body force), with `ρ_mix = ρ_light + φ(ρ_heavy − ρ_light)` and `ρ_ref` the
//!   domain mean of `ρ_mix`.  With `ρ_heavy = 1 + At`, `ρ_light = 1 − At` the
//!   inviscid linear growth rate is exactly `σ = sqrt(At g k)`.
//! * **Phase** — conservative Allen-Cahn equation (Chiu & Lin 2011)
//!   `∂φ/∂t + ∇·(φu) = ∇·[M(∇φ − (4/W) φ(1−φ) n̂)]`, solved with a D2Q9
//!   population whose equilibrium carries the sharpening flux
//!   `h_eq_i = w_i [φ + c_i·(φu + M (4/W) φ(1−φ) n̂)/cs²]`,
//!   `M = cs²(1/ω_φ − 1/2)`.
//!
//! Boundaries: periodic in x, half-way bounce-back (no-slip, no-flux) walls
//! below `y = 0` and above `y = ny − 1`.  Gravity points along −y; the heavy
//! fluid (φ = 1) initially sits on top of the light fluid.  All quantities are
//! in lattice units (`Δx = Δt = 1`).

use super::{EX9, EY9, LATTICE_CS2, W9};
use crate::error::{Error, Result};

/// Opposite D2Q9 directions.
const OPP9: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Numerical and initial-condition parameters of [`RayleighTaylorLbm`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayleighTaylorConfig {
    /// Relaxation frequency of the flow population (`ν = cs²(1/ω − 1/2)`).
    pub omega_flow: f64,
    /// Phase-field mobility `M` (sets the phase relaxation frequency).
    pub mobility: f64,
    /// Interface width `W` in lattice cells.
    pub interface_width: f64,
    /// Initial interface displacement amplitude in lattice cells.
    pub perturbation_amplitude: f64,
    /// Number of wavelengths of the initial perturbation across `nx`.
    pub perturbation_mode: usize,
}

impl Default for RayleighTaylorConfig {
    fn default() -> Self {
        Self {
            omega_flow: 1.9,
            mobility: 0.02,
            interface_width: 4.0,
            perturbation_amplitude: 1.0,
            perturbation_mode: 1,
        }
    }
}

/// Rayleigh-Taylor instability (two-population LBM, Boussinesq buoyancy).
///
/// See the module documentation for the model.
#[derive(Debug, Clone)]
pub struct RayleighTaylorLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Mixture density `ρ_light + φ(ρ_heavy − ρ_light)`.
    pub rho: Vec<Vec<f64>>,
    /// Phase field (0 = light, 1 = heavy).
    pub phase: Vec<Vec<f64>>,
    /// Velocity field (Guo-corrected).
    pub u: Vec<Vec<[f64; 2]>>,
    /// Atwood number At = (rho2-rho1)/(rho2+rho1).
    pub atwood: f64,
    /// Gravitational acceleration (lattice units), acting along −y.
    pub gravity: f64,
    /// Elapsed time in lattice steps.
    pub time: f64,
    /// Numerical parameters.
    pub config: RayleighTaylorConfig,
    /// Flow populations f\[y\]\[x\]\[q\].
    pub f: Vec<Vec<Vec<f64>>>,
    /// Phase populations h\[y\]\[x\]\[q\].
    pub h: Vec<Vec<Vec<f64>>>,
    /// Lattice (pressure) density of the flow population, `p = ρ_lbm cs²`.
    pub lattice_density: Vec<Vec<f64>>,
    /// Reusable streaming buffer (avoids per-step allocation).
    scratch: Vec<Vec<Vec<f64>>>,
}

impl RayleighTaylorLbm {
    /// Create a new Rayleigh-Taylor simulation with the default
    /// [`RayleighTaylorConfig`].
    pub fn new(nx: usize, ny: usize, atwood: f64, gravity: f64) -> Result<Self> {
        Self::with_config(nx, ny, atwood, gravity, RayleighTaylorConfig::default())
    }

    /// Create a new Rayleigh-Taylor simulation.
    ///
    /// Errors if `nx < 1` or `ny < 3`, `atwood ∉ [0, 1)`, `gravity` is
    /// negative or non-finite, `omega_flow ∉ (0, 2)`, or the mobility or
    /// interface width is not positive.
    pub fn with_config(
        nx: usize,
        ny: usize,
        atwood: f64,
        gravity: f64,
        config: RayleighTaylorConfig,
    ) -> Result<Self> {
        if nx < 1 || ny < 3 {
            return Err(Error::InvalidGridDimensions { nx, ny });
        }
        let check = |name: &str, value: f64, min: f64, max: f64, ok: bool| -> Result<()> {
            if ok && value.is_finite() {
                Ok(())
            } else {
                Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min,
                    max,
                })
            }
        };
        check("atwood", atwood, 0.0, 1.0, (0.0..1.0).contains(&atwood))?;
        check("gravity", gravity, 0.0, f64::MAX, gravity >= 0.0)?;
        check(
            "omega_flow",
            config.omega_flow,
            0.0,
            2.0,
            config.omega_flow > 0.0 && config.omega_flow < 2.0,
        )?;
        check(
            "mobility",
            config.mobility,
            0.0,
            f64::MAX,
            config.mobility > 0.0,
        )?;
        check(
            "interface_width",
            config.interface_width,
            0.0,
            f64::MAX,
            config.interface_width > 0.0,
        )?;
        check(
            "perturbation_amplitude",
            config.perturbation_amplitude,
            f64::MIN,
            f64::MAX,
            true,
        )?;

        let mut sim = Self {
            nx,
            ny,
            rho: vec![vec![0.0; nx]; ny],
            phase: vec![vec![0.0; nx]; ny],
            u: vec![vec![[0.0; 2]; nx]; ny],
            atwood,
            gravity,
            time: 0.0,
            config,
            f: vec![vec![vec![0.0; 9]; nx]; ny],
            h: vec![vec![vec![0.0; 9]; nx]; ny],
            lattice_density: vec![vec![1.0; nx]; ny],
            scratch: vec![vec![vec![0.0; 9]; nx]; ny],
        };
        sim.initialize();
        Ok(sim)
    }

    /// Heavy-fluid density `1 + At`.
    pub fn rho_heavy(&self) -> f64 {
        1.0 + self.atwood
    }

    /// Light-fluid density `1 − At`.
    pub fn rho_light(&self) -> f64 {
        1.0 - self.atwood
    }

    /// Wavenumber of the initial perturbation, `2π m / nx`.
    pub fn perturbation_wavenumber(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.config.perturbation_mode as f64 / self.nx as f64
    }

    /// Kinematic viscosity of the flow population.
    pub fn viscosity(&self) -> f64 {
        LATTICE_CS2 * (1.0 / self.config.omega_flow - 0.5)
    }

    /// Relaxation frequency of the phase population, from `M = cs²(1/ω − 1/2)`.
    pub fn omega_phase(&self) -> f64 {
        1.0 / (self.config.mobility / LATTICE_CS2 + 0.5)
    }

    /// Inviscid linear growth rate `sqrt(At g k)`.
    pub fn growth_rate(&self, k: f64) -> f64 {
        (self.atwood * self.gravity * k).sqrt()
    }

    /// Linear growth rate including viscous damping,
    /// `sqrt(At g k + ν²k⁴) − νk²`.
    pub fn viscous_growth_rate(&self, k: f64) -> f64 {
        let nu = self.viscosity();
        (self.atwood * self.gravity * k + nu * nu * k.powi(4)).sqrt() - nu * k * k
    }

    /// Compute bubble rise velocity (scaled).
    pub fn bubble_velocity(&self) -> f64 {
        // V_b = sqrt(At*g*L/(1+At)) where L = domain height
        let l = self.ny as f64;
        (self.atwood * self.gravity * l / (1.0 + self.atwood)).sqrt()
    }

    fn mixture_density(&self, phi: f64) -> f64 {
        self.rho_light() + phi * (self.rho_heavy() - self.rho_light())
    }

    /// Mean mixture density (reference density of the Boussinesq force).
    fn reference_density(&self) -> f64 {
        let total: f64 = self.rho.iter().flatten().sum();
        total / (self.nx * self.ny) as f64
    }

    /// Buoyancy force field `F_y = −(ρ_mix − ρ_ref) g`.
    fn buoyancy(&self) -> Vec<Vec<f64>> {
        let rho_ref = self.reference_density();
        self.rho
            .iter()
            .map(|row| row.iter().map(|&r| -(r - rho_ref) * self.gravity).collect())
            .collect()
    }

    /// Set the initial interface, hydrostatic lattice density and
    /// equilibrium populations.
    fn initialize(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let k = self.perturbation_wavenumber();
        let y_mid = 0.5 * (ny as f64 - 1.0);
        let width = self.config.interface_width;
        let amp = self.config.perturbation_amplitude;
        for y in 0..ny {
            for x in 0..nx {
                let y_int = y_mid + amp * (k * x as f64).cos();
                let phi = 0.5 * (1.0 + (2.0 * (y as f64 - y_int) / width).tanh());
                self.phase[y][x] = phi;
                self.rho[y][x] = self.mixture_density(phi);
            }
        }
        // Hydrostatic lattice density: cs² ∂ρ_lbm/∂y = F_y (trapezoidal rule).
        let force = self.buoyancy();
        for x in 0..nx {
            self.lattice_density[0][x] = 1.0;
        }
        for (y, pair) in force.windows(2).enumerate() {
            let below = self.lattice_density[y].clone();
            for (x, ((f_lo, f_hi), ld)) in pair[0]
                .iter()
                .zip(pair[1].iter())
                .zip(self.lattice_density[y + 1].iter_mut())
                .enumerate()
            {
                *ld = below[x] + 0.5 * (f_lo + f_hi) / LATTICE_CS2;
            }
        }
        let sum: f64 = self.lattice_density.iter().flatten().sum();
        let shift = 1.0 - sum / (nx * ny) as f64;
        for row in self.lattice_density.iter_mut() {
            for r in row.iter_mut() {
                *r += shift;
            }
        }
        let normals = self.interface_normals();
        let omega_phi = self.omega_phase();
        for y in 0..ny {
            for x in 0..nx {
                let rho = self.lattice_density[y][x];
                let shifted = [0.0, -0.5 * force[y][x] / rho];
                self.f[y][x] = flow_equilibrium(rho, shifted).to_vec();
                self.u[y][x] = [0.0, 0.0];
                self.h[y][x] = self
                    .phase_equilibrium(self.phase[y][x], [0.0, 0.0], normals[y][x], omega_phi)
                    .to_vec();
            }
        }
    }

    /// Phase value with mirror (zero-gradient) treatment beyond the walls
    /// and periodic wrap in x.
    fn phase_at(&self, x: i64, y: i64) -> f64 {
        let xi = x.rem_euclid(self.nx as i64) as usize;
        let yi = y.clamp(0, self.ny as i64 - 1) as usize;
        self.phase[yi][xi]
    }

    /// Unit interface normals `∇φ/|∇φ|` from the isotropic D2Q9 gradient
    /// stencil (zero where the gradient vanishes).
    fn interface_normals(&self) -> Vec<Vec<[f64; 2]>> {
        let mut normals = vec![vec![[0.0f64; 2]; self.nx]; self.ny];
        for (y, row) in normals.iter_mut().enumerate() {
            for (x, n) in row.iter_mut().enumerate() {
                let mut gx = 0.0;
                let mut gy = 0.0;
                for q in 1..9 {
                    let phi =
                        self.phase_at(x as i64 + i64::from(EX9[q]), y as i64 + i64::from(EY9[q]));
                    gx += W9[q] * f64::from(EX9[q]) * phi;
                    gy += W9[q] * f64::from(EY9[q]) * phi;
                }
                gx /= LATTICE_CS2;
                gy /= LATTICE_CS2;
                let mag = (gx * gx + gy * gy).sqrt();
                if mag > 1e-12 {
                    *n = [gx / mag, gy / mag];
                }
            }
        }
        normals
    }

    /// Conservative Allen-Cahn equilibrium.
    fn phase_equilibrium(
        &self,
        phi: f64,
        u: [f64; 2],
        normal: [f64; 2],
        omega_phi: f64,
    ) -> [f64; 9] {
        let mobility = LATTICE_CS2 * (1.0 / omega_phi - 0.5);
        let sharpen = mobility * 4.0 / self.config.interface_width * phi * (1.0 - phi);
        let jx = phi * u[0] + sharpen * normal[0];
        let jy = phi * u[1] + sharpen * normal[1];
        let mut heq = [0.0f64; 9];
        for (q, hq) in heq.iter_mut().enumerate() {
            let cj = f64::from(EX9[q]) * jx + f64::from(EY9[q]) * jy;
            *hq = W9[q] * (phi + cj / LATTICE_CS2);
        }
        heq
    }

    /// Advance one lattice time step (collide, stream, update macroscopic
    /// fields).
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega = self.config.omega_flow;
        let omega_phi = self.omega_phase();
        let force = self.buoyancy();
        let normals = self.interface_normals();

        // Collision.
        for y in 0..ny {
            for x in 0..nx {
                let rho = self.lattice_density[y][x];
                let u = self.u[y][x];
                let feq = flow_equilibrium(rho, u);
                let fy = force[y][x];
                for (q, fq) in self.f[y][x].iter_mut().enumerate() {
                    let cx = f64::from(EX9[q]);
                    let cy = f64::from(EY9[q]);
                    let cu = cx * u[0] + cy * u[1];
                    let src = (1.0 - 0.5 * omega)
                        * W9[q]
                        * ((cy - u[1]) / LATTICE_CS2 + cu * cy / (LATTICE_CS2 * LATTICE_CS2))
                        * fy;
                    *fq = *fq - omega * (*fq - feq[q]) + src;
                }
                let heq = self.phase_equilibrium(self.phase[y][x], u, normals[y][x], omega_phi);
                for (q, hq) in self.h[y][x].iter_mut().enumerate() {
                    *hq -= omega_phi * (*hq - heq[q]);
                }
            }
        }

        // Streaming with periodic x and half-way bounce-back walls in y.
        stream_with_walls(&self.f, &mut self.scratch, nx, ny);
        std::mem::swap(&mut self.f, &mut self.scratch);
        stream_with_walls(&self.h, &mut self.scratch, nx, ny);
        std::mem::swap(&mut self.h, &mut self.scratch);

        // Macroscopic fields.
        for y in 0..ny {
            for x in 0..nx {
                let phi: f64 = self.h[y][x].iter().sum();
                self.phase[y][x] = phi;
                self.rho[y][x] = self.mixture_density(phi);
            }
        }
        let force = self.buoyancy();
        for (y, force_row) in force.iter().enumerate() {
            for (x, &fy) in force_row.iter().enumerate() {
                let cell = &self.f[y][x];
                let rho: f64 = cell.iter().sum();
                let mut mx = 0.0;
                let mut my = 0.0;
                for (q, &fq) in cell.iter().enumerate() {
                    mx += f64::from(EX9[q]) * fq;
                    my += f64::from(EY9[q]) * fq;
                }
                self.lattice_density[y][x] = rho;
                if rho > 0.0 {
                    self.u[y][x] = [mx / rho, (my + 0.5 * fy) / rho];
                }
            }
        }
        self.time += 1.0;
    }

    /// Total amount of heavy phase `Σ φ` (conserved by the scheme).
    pub fn total_phase(&self) -> f64 {
        self.phase.iter().flatten().sum()
    }

    /// Interface height (φ = 1/2 contour) in each column, linearly
    /// interpolated; `None` for a column without a light→heavy crossing.
    pub fn interface_heights(&self) -> Vec<Option<f64>> {
        (0..self.nx)
            .map(|x| {
                (0..self.ny - 1).find_map(|y| {
                    let lo = self.phase[y][x];
                    let hi = self.phase[y + 1][x];
                    if lo < 0.5 && hi >= 0.5 {
                        Some(y as f64 + (0.5 - lo) / (hi - lo))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    /// Fourier amplitude of the interface height for mode `m`,
    /// `(2/nx) Σ h(x) cos(2π m x/nx)`.  `None` if some column has no
    /// interface crossing.
    pub fn interface_amplitude(&self, mode: usize) -> Option<f64> {
        let k = 2.0 * std::f64::consts::PI * mode as f64 / self.nx as f64;
        let heights = self.interface_heights();
        let mut acc = 0.0;
        for (x, h) in heights.iter().enumerate() {
            acc += (*h)? * (k * x as f64).cos();
        }
        Some(2.0 * acc / self.nx as f64)
    }
}

/// Second-order D2Q9 equilibrium at `cs² = 1/3`.
fn flow_equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
    let usq = u[0] * u[0] + u[1] * u[1];
    let mut feq = [0.0f64; 9];
    for (q, fq) in feq.iter_mut().enumerate() {
        let cu = f64::from(EX9[q]) * u[0] + f64::from(EY9[q]) * u[1];
        *fq = W9[q] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * usq);
    }
    feq
}

/// Stream populations into `out` with periodic x and half-way bounce-back
/// at the bottom and top walls.  Every slot of `out` is overwritten.
fn stream_with_walls(pops: &[Vec<Vec<f64>>], out: &mut [Vec<Vec<f64>>], nx: usize, ny: usize) {
    for (y, row) in pops.iter().enumerate() {
        for (x, cell) in row.iter().enumerate() {
            for (q, &value) in cell.iter().enumerate() {
                let yn = y as i64 + i64::from(EY9[q]);
                if yn < 0 || yn >= ny as i64 {
                    out[y][x][OPP9[q]] = value;
                } else {
                    let xn = (x as i64 + i64::from(EX9[q])).rem_euclid(nx as i64) as usize;
                    out[yn as usize][xn][q] = value;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_degenerate_grids() {
        for (nx, ny) in [(16, 0), (16, 1), (16, 2), (0, 32)] {
            let r = RayleighTaylorLbm::new(nx, ny, 0.1, 1e-4);
            assert!(
                matches!(r, Err(Error::InvalidGridDimensions { .. })),
                "nx={nx} ny={ny}"
            );
        }
        assert!(RayleighTaylorLbm::new(16, 3, 0.1, 1e-4).is_ok());
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(RayleighTaylorLbm::new(16, 32, 1.0, 1e-4).is_err());
        assert!(RayleighTaylorLbm::new(16, 32, -0.1, 1e-4).is_err());
        assert!(RayleighTaylorLbm::new(16, 32, 0.1, -1.0).is_err());
        assert!(RayleighTaylorLbm::new(16, 32, 0.1, f64::NAN).is_err());
        let cfg = RayleighTaylorConfig {
            omega_flow: 2.0,
            ..RayleighTaylorConfig::default()
        };
        assert!(RayleighTaylorLbm::with_config(16, 32, 0.1, 1e-4, cfg).is_err());
    }

    #[test]
    fn heavy_fluid_starts_on_top() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.2, 1e-4).expect("valid");
        assert!(rt.phase[31][0] > 0.99);
        assert!(rt.phase[0][0] < 0.01);
        assert!((rt.rho[31][0] - 1.2).abs() < 1e-2);
        assert!((rt.rho[0][0] - 0.8).abs() < 1e-2);
        let amp = rt.interface_amplitude(1).expect("interface present");
        assert!((amp - 1.0).abs() < 0.05, "initial amplitude {amp}");
    }

    #[test]
    fn hydrostatic_flat_interface_stays_at_rest() {
        let cfg = RayleighTaylorConfig {
            perturbation_amplitude: 0.0,
            ..RayleighTaylorConfig::default()
        };
        let mut rt = RayleighTaylorLbm::with_config(8, 32, 0.3, 1e-4, cfg).expect("valid");
        for _ in 0..200 {
            rt.step();
        }
        let vmax =
            rt.u.iter()
                .flatten()
                .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
                .fold(0.0f64, f64::max);
        assert!(vmax < 1e-4, "spurious velocity {vmax}");
    }

    #[test]
    fn phase_is_conserved_and_transported() {
        let mut rt = RayleighTaylorLbm::new(32, 64, 0.5, 2e-4).expect("valid");
        let total0 = rt.total_phase();
        let phase0 = rt.phase.clone();
        for _ in 0..300 {
            rt.step();
        }
        let total1 = rt.total_phase();
        assert!(
            (total1 - total0).abs() < 1e-9 * total0,
            "phase drift {}",
            total1 - total0
        );
        let changed = rt
            .phase
            .iter()
            .flatten()
            .zip(phase0.iter().flatten())
            .filter(|(a, b)| (*a - *b).abs() > 1e-3)
            .count();
        assert!(changed > 10, "phase never evolved");
        // Density follows the phase field.
        let (y, x) = (40, 5);
        let expected = rt.rho_light() + rt.phase[y][x] * (rt.rho_heavy() - rt.rho_light());
        assert!((rt.rho[y][x] - expected).abs() < 1e-12);
    }

    /// Least-squares fit of `A(t) = a cosh(σt) + b sinh(σt)` (the general
    /// linear solution, both growing and decaying modes) by a scan over σ.
    fn fit_growth_rate(times: &[f64], amps: &[f64], sigma_guess: f64) -> f64 {
        let mut best = (f64::INFINITY, 0.0);
        for i in 1..=4000 {
            let sigma = sigma_guess * f64::from(i) / 2000.0;
            let c: Vec<f64> = times.iter().map(|t| (sigma * t).cosh()).collect();
            let s: Vec<f64> = times.iter().map(|t| (sigma * t).sinh()).collect();
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
            let (cc, ss, cs) = (dot(&c, &c), dot(&s, &s), dot(&c, &s));
            let (ca, sa) = (dot(&c, amps), dot(&s, amps));
            let det = cc * ss - cs * cs;
            let a = (ca * ss - sa * cs) / det;
            let b = (sa * cc - ca * cs) / det;
            let resid: f64 = (0..times.len())
                .map(|n| (a * c[n] + b * s[n] - amps[n]).powi(2))
                .sum();
            if resid < best.0 {
                best = (resid, sigma);
            }
        }
        best.1
    }

    #[test]
    fn linear_growth_rate_matches_theory() {
        // λ = 128 cells keeps k·W small (diffuse-interface correction) and
        // ω = 1.95 keeps viscous damping at ~0.5% of σ.
        let cfg = RayleighTaylorConfig {
            omega_flow: 1.95,
            mobility: 0.01,
            interface_width: 4.0,
            perturbation_amplitude: 1.0,
            perturbation_mode: 1,
        };
        let mut rt = RayleighTaylorLbm::with_config(128, 128, 0.5, 2e-4, cfg).expect("valid");
        let k = rt.perturbation_wavenumber();
        let theory = rt.growth_rate(k);
        let mut times = Vec::new();
        let mut amps = Vec::new();
        for t in 0..=1000usize {
            if t % 20 == 0 {
                times.push(t as f64);
                amps.push(rt.interface_amplitude(1).expect("interface"));
            }
            rt.step();
        }
        // Still in the linear regime (k·A < 0.3) and clearly growing.
        let last = amps.last().copied().unwrap_or(0.0);
        assert!(last > 2.5 && k * last < 0.3, "final amplitude {last}");
        let sigma = fit_growth_rate(&times, &amps, theory);
        let rel = (sigma - theory).abs() / theory;
        assert!(rel < 0.15, "σ = {sigma}, theory = {theory}, rel = {rel}");
    }
}
