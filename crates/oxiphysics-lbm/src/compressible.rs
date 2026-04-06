#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compressible flow LBM.
//!
//! Provides high-Mach LBM with extended equilibrium, Euler equations via LBM,
//! shock capturing, compressible Navier-Stokes, and detonation wave simulation.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

/// High-Mach LBM simulation with extended equilibrium.
///
/// Uses higher-order equilibrium to capture compressible effects.
#[derive(Debug, Clone)]
pub struct CompressibleLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation parameter omega.
    pub omega: f64,
    /// Distribution functions f\[y\]\[x\]\[q\].
    pub f: Vec<Vec<Vec<f64>>>,
    /// Density field rho\[y\]\[x\].
    pub rho: Vec<Vec<f64>>,
    /// Velocity field u\[y\]\[x\]\[2\].
    pub u: Vec<Vec<[f64; 2]>>,
    /// Temperature field T\[y\]\[x\].
    pub temp: Vec<Vec<f64>>,
    /// Speed of sound cs = sqrt(gamma*R*T).
    pub cs: f64,
    /// Adiabatic exponent gamma.
    pub gamma: f64,
    /// Mach number correction order.
    pub ma_order: usize,
}

/// D2Q9 weights.
const W9: [f64; 9] = [
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
/// D2Q9 velocity ex.
const EX9: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 velocity ey.
const EY9: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Compute standard D2Q9 equilibrium with compressible correction.
pub fn equilibrium_compressible(rho: f64, ux: f64, uy: f64, cs2: f64) -> [f64; 9] {
    let mut feq = [0.0f64; 9];
    let usq = ux * ux + uy * uy;
    let inv_cs2 = 1.0 / cs2;
    let inv_cs4 = inv_cs2 * inv_cs2;
    for q in 0..9 {
        let eu = EX9[q] as f64 * ux + EY9[q] as f64 * uy;
        // Extended equilibrium (second-order Ma)
        feq[q] = W9[q]
            * rho
            * (1.0
                + eu * inv_cs2
                + 0.5 * (eu * eu * inv_cs4 - usq * inv_cs2)
                + eu * eu * eu * inv_cs4 / 6.0
                - eu * usq * inv_cs4 / 2.0);
    }
    feq
}

/// Compute local Mach number.
///
/// Returns Ma = |u| / cs.
pub fn local_mach(ux: f64, uy: f64, cs: f64) -> f64 {
    let u_mag = (ux * ux + uy * uy).sqrt();
    u_mag / cs
}

/// Compute entropy generation rate.
///
/// Returns sigma = rho * T * Ds/Dt estimate.
pub fn entropy_generation(rho: f64, temp: f64, du_dx: f64, dv_dy: f64, mu: f64) -> f64 {
    let phi_visc = 2.0 * mu * (du_dx * du_dx + dv_dy * dv_dy);
    phi_visc / (rho * temp + 1e-30)
}

impl CompressibleLbm {
    /// Create a new compressible LBM simulation.
    pub fn new(nx: usize, ny: usize, omega: f64, gamma: f64, ma_order: usize) -> Self {
        let cs = (gamma / 3.0).sqrt();
        let f = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let rho = vec![vec![1.0f64; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let temp = vec![vec![1.0f64; nx]; ny];
        let mut sim = Self {
            nx,
            ny,
            omega,
            f,
            rho,
            u,
            temp,
            cs,
            gamma,
            ma_order,
        };
        sim.initialize();
        sim
    }

    /// Initialize equilibrium distributions.
    pub fn initialize(&mut self) {
        let cs2 = self.cs * self.cs;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let feq = equilibrium_compressible(1.0, 0.0, 0.0, cs2);
                self.f[y][x] = feq.to_vec();
            }
        }
    }

    /// BGK collision step.
    pub fn collide(&mut self) {
        let cs2 = self.cs * self.cs;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.rho[y][x];
                let [ux, uy] = self.u[y][x];
                let feq = equilibrium_compressible(rho, ux, uy, cs2);
                for q in 0..9 {
                    self.f[y][x][q] = self.f[y][x][q] * (1.0 - self.omega) + self.omega * feq[q];
                }
            }
        }
    }

    /// Streaming step.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let xn = (x as i32 + EX9[q]).rem_euclid(nx as i32) as usize;
                    let yn = (y as i32 + EY9[q]).rem_euclid(ny as i32) as usize;
                    f_new[yn][xn][q] = self.f[y][x][q];
                }
            }
        }
        self.f = f_new;
    }

    /// Update macroscopic variables.
    pub fn update_macros(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for q in 0..9 {
                    rho += self.f[y][x][q];
                    ux += EX9[q] as f64 * self.f[y][x][q];
                    uy += EY9[q] as f64 * self.f[y][x][q];
                }
                self.rho[y][x] = rho.max(1e-10);
                if rho > 1e-10 {
                    self.u[y][x] = [ux / rho, uy / rho];
                }
                // Temperature from equation of state: p = rho*cs2
                self.temp[y][x] = rho * self.cs * self.cs / (self.gamma - 1.0 + 1e-30);
            }
        }
    }

    /// Full time step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }

    /// Local Mach number at (x, y).
    pub fn mach_at(&self, x: usize, y: usize) -> f64 {
        let [ux, uy] = self.u[y][x];
        local_mach(ux, uy, self.cs)
    }

    /// Maximum Mach number in the domain.
    pub fn max_mach(&self) -> f64 {
        let mut ma_max = 0.0f64;
        for y in 0..self.ny {
            for x in 0..self.nx {
                ma_max = ma_max.max(self.mach_at(x, y));
            }
        }
        ma_max
    }
}

/// Euler equations via LBM (inviscid compressible flow).
///
/// Uses zero viscosity limit of LBM (omega=2) for Euler equations.
#[derive(Debug, Clone)]
pub struct EulerLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Energy distribution functions g.
    pub g: Vec<Vec<Vec<f64>>>,
    /// Density.
    pub rho: Vec<Vec<f64>>,
    /// Velocity.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Energy density.
    pub energy: Vec<Vec<f64>>,
    /// Speed of sound.
    pub cs: f64,
}

impl EulerLbm {
    /// Create a new Euler LBM.
    pub fn new(nx: usize, ny: usize, cs: f64) -> Self {
        let f = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let g = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let rho = vec![vec![1.0f64; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let energy = vec![vec![1.5f64; nx]; ny];
        let mut sim = Self {
            nx,
            ny,
            f,
            g,
            rho,
            u,
            energy,
            cs,
        };
        sim.initialize();
        sim
    }

    /// Initialize distributions.
    pub fn initialize(&mut self) {
        let cs2 = self.cs * self.cs;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let feq = equilibrium_compressible(1.0, 0.0, 0.0, cs2);
                self.f[y][x] = feq.to_vec();
                self.g[y][x] = feq.to_vec(); // simple energy init
            }
        }
    }

    /// Collision (omega=2 for inviscid Euler).
    pub fn collide(&mut self) {
        let cs2 = self.cs * self.cs;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.rho[y][x];
                let [ux, uy] = self.u[y][x];
                let feq = equilibrium_compressible(rho, ux, uy, cs2);
                self.f[y][x].copy_from_slice(&feq); // omega = 2: full relaxation
            }
        }
    }

    /// Streaming step.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let xn = (x as i32 + EX9[q]).rem_euclid(nx as i32) as usize;
                    let yn = (y as i32 + EY9[q]).rem_euclid(ny as i32) as usize;
                    f_new[yn][xn][q] = self.f[y][x][q];
                }
            }
        }
        self.f = f_new;
    }

    /// Update macroscopic variables.
    pub fn update_macros(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for q in 0..9 {
                    rho += self.f[y][x][q];
                    ux += EX9[q] as f64 * self.f[y][x][q];
                    uy += EY9[q] as f64 * self.f[y][x][q];
                }
                self.rho[y][x] = rho.max(1e-10);
                if rho > 1e-10 {
                    self.u[y][x] = [ux / rho, uy / rho];
                    let ke = 0.5 * rho * (ux * ux + uy * uy) / (rho * rho);
                    self.energy[y][x] = rho * self.cs * self.cs + ke;
                }
            }
        }
    }

    /// Full Euler step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }

    /// Compute Roe flux at interface between left and right states.
    pub fn roe_flux_x(
        &self,
        rho_l: f64,
        u_l: [f64; 2],
        p_l: f64,
        rho_r: f64,
        u_r: [f64; 2],
        p_r: f64,
    ) -> [f64; 4] {
        euler_flux_roe(rho_l, u_l, p_l, rho_r, u_r, p_r)
    }
}

/// Shock-capturing LBM with artificial viscosity.
///
/// Adds bulk viscosity term and artificial dissipation near shocks.
#[derive(Debug, Clone)]
pub struct ShockCaptureLbm {
    /// Underlying compressible LBM.
    pub lbm: CompressibleLbm,
    /// Artificial viscosity coefficient.
    pub nu_art: f64,
    /// Shock sensor threshold.
    pub shock_threshold: f64,
    /// Shock indicator field.
    pub shock_flag: Vec<Vec<bool>>,
}

impl ShockCaptureLbm {
    /// Create a new shock-capturing LBM.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        gamma: f64,
        nu_art: f64,
        shock_threshold: f64,
    ) -> Self {
        let lbm = CompressibleLbm::new(nx, ny, omega, gamma, 2);
        let shock_flag = vec![vec![false; nx]; ny];
        Self {
            lbm,
            nu_art,
            shock_threshold,
            shock_flag,
        }
    }

    /// Detect shocks using pressure gradient sensor.
    pub fn detect_shocks(&mut self) {
        let nx = self.lbm.nx;
        let ny = self.lbm.ny;
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let rho_c = self.lbm.rho[y][x];
                let rho_r = self.lbm.rho[y][x + 1];
                let rho_l = self.lbm.rho[y][x - 1];
                let sensor = ((rho_r - 2.0 * rho_c + rho_l).abs()) / (rho_c + 1e-10);
                self.shock_flag[y][x] = sensor > self.shock_threshold;
            }
        }
    }

    /// Apply artificial viscosity at shock cells.
    pub fn apply_artificial_viscosity(&mut self) {
        let nx = self.lbm.nx;
        let ny = self.lbm.ny;
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                if self.shock_flag[y][x] {
                    // Add diffusion: smear distribution functions
                    for q in 0..9 {
                        let lap = self.lbm.f[y][x + 1][q] - 2.0 * self.lbm.f[y][x][q]
                            + self.lbm.f[y][x - 1][q];
                        self.lbm.f[y][x][q] += self.nu_art * lap;
                    }
                }
            }
        }
    }

    /// Full time step with shock capturing.
    pub fn step(&mut self) {
        self.lbm.step();
        self.detect_shocks();
        self.apply_artificial_viscosity();
    }

    /// Count number of shock cells.
    pub fn shock_cell_count(&self) -> usize {
        self.shock_flag
            .iter()
            .flat_map(|row| row.iter())
            .filter(|&&b| b)
            .count()
    }
}

/// Compressible Navier-Stokes via double distribution (density + energy).
///
/// Uses separate f (density/momentum) and g (energy) distributions.
#[derive(Debug, Clone)]
pub struct NsCompressibleLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Density/momentum distribution f.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Energy distribution g.
    pub g: Vec<Vec<Vec<f64>>>,
    /// Relaxation for f.
    pub omega_f: f64,
    /// Relaxation for g.
    pub omega_g: f64,
    /// Density field.
    pub rho: Vec<Vec<f64>>,
    /// Velocity field.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Temperature field.
    pub temp: Vec<Vec<f64>>,
    /// Prandtl number.
    pub pr: f64,
    /// Gamma.
    pub gamma: f64,
}

impl NsCompressibleLbm {
    /// Create a new NS compressible LBM.
    pub fn new(nx: usize, ny: usize, omega_f: f64, omega_g: f64, pr: f64, gamma: f64) -> Self {
        let f = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let g = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let rho = vec![vec![1.0f64; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let temp = vec![vec![1.0f64; nx]; ny];
        let mut sim = Self {
            nx,
            ny,
            f,
            g,
            omega_f,
            omega_g,
            rho,
            u,
            temp,
            pr,
            gamma,
        };
        sim.initialize();
        sim
    }

    /// Initialize distributions.
    pub fn initialize(&mut self) {
        let cs2 = 1.0 / 3.0;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let feq = equilibrium_compressible(1.0, 0.0, 0.0, cs2);
                self.f[y][x] = feq.to_vec();
                // Energy distribution: g_eq = (cs2 + 0.5*u^2) * feq + internal energy correction
                let geq: Vec<f64> = feq.iter().map(|&fq| fq * 1.5 * cs2).collect();
                self.g[y][x] = geq;
            }
        }
    }

    /// Collision step for both f and g.
    pub fn collide(&mut self) {
        let cs2 = 1.0 / 3.0;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.rho[y][x];
                let [ux, uy] = self.u[y][x];
                let t = self.temp[y][x];
                let feq = equilibrium_compressible(rho, ux, uy, cs2);
                let e_tot = 0.5 * (ux * ux + uy * uy) + t * cs2 / (self.gamma - 1.0);
                let geq: Vec<f64> = feq
                    .iter()
                    .enumerate()
                    .map(|(q, &fq)| {
                        let eu = EX9[q] as f64 * ux + EY9[q] as f64 * uy;
                        fq * (e_tot + eu)
                    })
                    .collect();
                for q in 0..9 {
                    self.f[y][x][q] =
                        self.f[y][x][q] * (1.0 - self.omega_f) + self.omega_f * feq[q];
                    self.g[y][x][q] =
                        self.g[y][x][q] * (1.0 - self.omega_g) + self.omega_g * geq[q];
                }
            }
        }
    }

    /// Streaming step.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let mut g_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let xn = (x as i32 + EX9[q]).rem_euclid(nx as i32) as usize;
                    let yn = (y as i32 + EY9[q]).rem_euclid(ny as i32) as usize;
                    f_new[yn][xn][q] = self.f[y][x][q];
                    g_new[yn][xn][q] = self.g[y][x][q];
                }
            }
        }
        self.f = f_new;
        self.g = g_new;
    }

    /// Update macroscopic variables.
    pub fn update_macros(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                let mut energy = 0.0;
                for q in 0..9 {
                    rho += self.f[y][x][q];
                    ux += EX9[q] as f64 * self.f[y][x][q];
                    uy += EY9[q] as f64 * self.f[y][x][q];
                    energy += self.g[y][x][q];
                }
                self.rho[y][x] = rho.max(1e-10);
                if rho > 1e-10 {
                    self.u[y][x] = [ux / rho, uy / rho];
                    let ke = 0.5 * (ux * ux + uy * uy) / rho;
                    let internal = (energy - ke).max(0.0);
                    self.temp[y][x] = internal * (self.gamma - 1.0) / rho;
                }
            }
        }
    }

    /// Full NS compressible step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }
}

/// Mach number correction for extended equilibrium.
///
/// Provides higher-order terms for compressible equilibrium.
#[derive(Debug, Clone)]
pub struct MaCorrection {
    /// Maximum Mach number order (2 or 3).
    pub order: usize,
    /// Reference speed of sound.
    pub cs: f64,
}

impl MaCorrection {
    /// Create a new Mach correction calculator.
    pub fn new(order: usize, cs: f64) -> Self {
        Self { order, cs }
    }

    /// Compute corrected equilibrium at given Ma.
    pub fn corrected_feq(&self, rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let cs2 = self.cs * self.cs;
        let mut feq = [0.0f64; 9];
        let usq = ux * ux + uy * uy;
        for q in 0..9 {
            let eu = EX9[q] as f64 * ux + EY9[q] as f64 * uy;
            let base = 1.0 + eu / cs2 + 0.5 * (eu * eu / (cs2 * cs2) - usq / cs2);
            let correction = if self.order >= 3 {
                eu * eu * eu / (6.0 * cs2 * cs2 * cs2) - eu * usq / (2.0 * cs2 * cs2)
            } else {
                0.0
            };
            feq[q] = W9[q] * rho * (base + correction);
        }
        feq
    }

    /// Estimate Mach correction error at given Ma.
    pub fn truncation_error(&self, ma: f64) -> f64 {
        ma.powi(self.order as i32 + 1)
    }
}

/// 1D Riemann problem (Sod shock tube test).
///
/// Implements contact discontinuity and shock wave evolution.
#[derive(Debug, Clone)]
pub struct ShockWave1D {
    /// Number of cells.
    pub n: usize,
    /// Density profile.
    pub rho: Vec<f64>,
    /// Velocity profile.
    pub vel: Vec<f64>,
    /// Pressure profile.
    pub pressure: Vec<f64>,
    /// Energy profile.
    pub energy: Vec<f64>,
    /// Grid spacing dx.
    pub dx: f64,
    /// Time step dt.
    pub dt: f64,
    /// Adiabatic index gamma.
    pub gamma: f64,
    /// Current time.
    pub time: f64,
}

impl ShockWave1D {
    /// Create a new Sod shock tube problem.
    pub fn new_sod(n: usize, dx: f64, dt: f64, gamma: f64) -> Self {
        let mut rho = vec![0.0f64; n];
        let mut pressure = vec![0.0f64; n];
        let vel = vec![0.0f64; n];
        let mid = n / 2;
        for i in 0..n {
            if i < mid {
                rho[i] = 1.0;
                pressure[i] = 1.0;
            } else {
                rho[i] = 0.125;
                pressure[i] = 0.1;
            }
        }
        let energy: Vec<f64> = (0..n)
            .map(|i| pressure[i] / (gamma - 1.0) + 0.5 * rho[i] * vel[i] * vel[i])
            .collect();
        Self {
            n,
            rho,
            vel,
            pressure,
            energy,
            dx,
            dt,
            gamma,
            time: 0.0,
        }
    }

    /// Compute local speed of sound.
    pub fn sound_speed(&self, i: usize) -> f64 {
        let p = self.pressure[i].max(0.0);
        let rho = self.rho[i].max(1e-10);
        (self.gamma * p / rho).sqrt()
    }

    /// Compute CFL number.
    pub fn cfl(&self) -> f64 {
        let max_speed: f64 = (0..self.n)
            .map(|i| self.vel[i].abs() + self.sound_speed(i))
            .fold(0.0f64, f64::max);
        max_speed * self.dt / self.dx
    }

    /// First-order Godunov step (upwind).
    pub fn step_godunov(&mut self) {
        let n = self.n;
        let dt = self.dt;
        let dx = self.dx;
        let g = self.gamma;
        let mut rho_new = self.rho.clone();
        let mut vel_new = self.vel.clone();
        let mut e_new = self.energy.clone();

        for i in 1..n - 1 {
            // Left and right states
            let rho_l = self.rho[i - 1];
            let rho_r = self.rho[i];
            let u_l = self.vel[i - 1];
            let u_r = self.vel[i];
            let p_l = (g - 1.0) * (self.energy[i - 1] - 0.5 * rho_l * u_l * u_l);
            let p_r = (g - 1.0) * (self.energy[i] - 0.5 * rho_r * u_r * u_r);

            // Simple upwind flux
            let flux_rho_l = if u_l > 0.0 { rho_l * u_l } else { rho_r * u_r };
            let flux_mom_l = if u_l > 0.0 {
                rho_l * u_l * u_l + p_l
            } else {
                rho_r * u_r * u_r + p_r
            };
            let flux_e_l = if u_l > 0.0 {
                (self.energy[i - 1] + p_l) * u_l
            } else {
                (self.energy[i] + p_r) * u_r
            };

            let rho_rr = self.rho[(i + 1).min(n - 1)];
            let u_rr = self.vel[(i + 1).min(n - 1)];
            let p_rr = (g - 1.0) * (self.energy[(i + 1).min(n - 1)] - 0.5 * rho_rr * u_rr * u_rr);
            let flux_rho_r = if u_r > 0.0 {
                rho_r * u_r
            } else {
                rho_rr * u_rr
            };
            let flux_mom_r = if u_r > 0.0 {
                rho_r * u_r * u_r + p_r
            } else {
                rho_rr * u_rr * u_rr + p_rr
            };
            let flux_e_r = if u_r > 0.0 {
                (self.energy[i] + p_r) * u_r
            } else {
                (self.energy[(i + 1).min(n - 1)] + p_rr) * u_rr
            };

            rho_new[i] = self.rho[i] - dt / dx * (flux_rho_r - flux_rho_l);
            let mom = self.rho[i] * self.vel[i] - dt / dx * (flux_mom_r - flux_mom_l);
            e_new[i] = self.energy[i] - dt / dx * (flux_e_r - flux_e_l);
            rho_new[i] = rho_new[i].max(1e-10);
            vel_new[i] = mom / rho_new[i];
            e_new[i] = e_new[i].max(1e-10);
        }
        self.rho = rho_new;
        self.vel = vel_new;
        self.energy = e_new;
        self.pressure = (0..n)
            .map(|i| {
                let p =
                    (g - 1.0) * (self.energy[i] - 0.5 * self.rho[i] * self.vel[i] * self.vel[i]);
                p.max(0.0)
            })
            .collect();
        self.time += dt;
    }

    /// Run for n steps.
    pub fn run(&mut self, nsteps: usize) {
        for _ in 0..nsteps {
            self.step_godunov();
        }
    }
}

/// Rayleigh-Taylor instability in compressible regime.
///
/// Simulates gravitational instability between dense and light fluids.
#[derive(Debug, Clone)]
pub struct RayleighTaylorLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Density field (two-fluid mixture).
    pub rho: Vec<Vec<f64>>,
    /// Phase field (0 = light, 1 = heavy).
    pub phase: Vec<Vec<f64>>,
    /// Velocity field.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Atwood number At = (rho2-rho1)/(rho2+rho1).
    pub atwood: f64,
    /// Gravity.
    pub gravity: f64,
    /// Time.
    pub time: f64,
    /// Time step.
    pub dt: f64,
}

impl RayleighTaylorLbm {
    /// Create a new Rayleigh-Taylor simulation.
    pub fn new(nx: usize, ny: usize, atwood: f64, gravity: f64, dt: f64) -> Self {
        let rho_heavy = 1.0 + atwood;
        let rho_light = 1.0 - atwood;
        let mut rho = vec![vec![0.0f64; nx]; ny];
        let mut phase = vec![vec![0.0f64; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                let perturbation = 0.05 * (2.0 * std::f64::consts::PI * x as f64 / nx as f64).cos();
                let interface = ny as f64 / 2.0 + perturbation * ny as f64 * 0.1;
                if (y as f64) < interface {
                    rho[y][x] = rho_heavy;
                    phase[y][x] = 1.0;
                } else {
                    rho[y][x] = rho_light;
                    phase[y][x] = 0.0;
                }
            }
        }
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        Self {
            nx,
            ny,
            rho,
            phase,
            u,
            atwood,
            gravity,
            time: 0.0,
            dt,
        }
    }

    /// Compute growth rate of RT instability.
    pub fn growth_rate(&self, k: f64) -> f64 {
        // gamma = sqrt(At * g * k)
        (self.atwood * self.gravity * k).sqrt()
    }

    /// Compute bubble rise velocity (scaled).
    pub fn bubble_velocity(&self) -> f64 {
        // V_b = sqrt(At*g*L/(1+At)) where L = domain height
        let l = self.ny as f64;
        (self.atwood * self.gravity * l / (1.0 + self.atwood)).sqrt()
    }

    /// Advance one step (simple density diffusion model).
    pub fn step(&mut self) {
        let dt = self.dt;
        let g = self.gravity;
        for y in 1..self.ny - 1 {
            for x in 0..self.nx {
                let rho_below = self.rho[y - 1][x];
                let rho_above = self.rho[y + 1][x];
                // Buoyancy: dense fluid sinks
                let f_buoy = (rho_below - rho_above) * g * 0.5;
                self.u[y][x][1] += f_buoy * dt / self.rho[y][x];
            }
        }
        self.time += dt;
    }
}

/// Supersonic inlet boundary condition.
///
/// Specifies inflow Mach number, angle of attack, density, and pressure.
#[derive(Debug, Clone)]
pub struct SupersonicInlet {
    /// Inflow Mach number.
    pub mach: f64,
    /// Angle of attack in radians.
    pub angle: f64,
    /// Inflow density.
    pub rho_in: f64,
    /// Inflow pressure.
    pub p_in: f64,
    /// Adiabatic index.
    pub gamma: f64,
    /// Speed of sound.
    pub cs: f64,
}

impl SupersonicInlet {
    /// Create a new supersonic inlet.
    pub fn new(mach: f64, angle: f64, rho_in: f64, p_in: f64, gamma: f64) -> Self {
        let cs = (gamma * p_in / rho_in).sqrt();
        Self {
            mach,
            angle,
            rho_in,
            p_in,
            gamma,
            cs,
        }
    }

    /// Inflow velocity vector \[ux, uy\].
    pub fn velocity(&self) -> [f64; 2] {
        let u_mag = self.mach * self.cs;
        [u_mag * self.angle.cos(), u_mag * self.angle.sin()]
    }

    /// Apply inlet BC to distribution functions at column x=0.
    pub fn apply(&self, f: &mut Vec<Vec<Vec<f64>>>, ny: usize) {
        let cs2 = self.cs * self.cs;
        let [ux, uy] = self.velocity();
        let feq = equilibrium_compressible(self.rho_in, ux, uy, cs2);
        for y in 0..ny {
            f[y][0] = feq.to_vec();
        }
    }

    /// Total pressure from isentropic relations.
    pub fn total_pressure(&self) -> f64 {
        self.p_in
            * (1.0 + (self.gamma - 1.0) / 2.0 * self.mach * self.mach)
                .powf(self.gamma / (self.gamma - 1.0))
    }

    /// Total temperature ratio T0/T.
    pub fn total_temperature_ratio(&self) -> f64 {
        1.0 + (self.gamma - 1.0) / 2.0 * self.mach * self.mach
    }
}

/// Compressible flow statistics.
///
/// Computes Mach number, pressure ratio, density ratio, and entropy production.
#[derive(Debug, Clone)]
pub struct CompressibleStats {
    /// Density field snapshot.
    pub rho: Vec<Vec<f64>>,
    /// Velocity field snapshot.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Pressure field snapshot.
    pub p: Vec<Vec<f64>>,
    /// Speed of sound field.
    pub cs_field: Vec<Vec<f64>>,
    /// Gamma.
    pub gamma: f64,
    /// Reference density.
    pub rho_ref: f64,
    /// Reference pressure.
    pub p_ref: f64,
}

impl CompressibleStats {
    /// Create new stats from snapshot.
    pub fn new(rho: Vec<Vec<f64>>, u: Vec<Vec<[f64; 2]>>, p: Vec<Vec<f64>>, gamma: f64) -> Self {
        let ny = rho.len();
        let nx = if ny > 0 { rho[0].len() } else { 0 };
        let rho_ref = 1.0;
        let p_ref = 1.0;
        let cs_field = (0..ny)
            .map(|y| {
                (0..nx)
                    .map(|x| (gamma * p[y][x] / rho[y][x].max(1e-10)).sqrt())
                    .collect()
            })
            .collect();
        Self {
            rho,
            u,
            p,
            cs_field,
            gamma,
            rho_ref,
            p_ref,
        }
    }

    /// Local Mach number at (x, y).
    pub fn mach_at(&self, x: usize, y: usize) -> f64 {
        let [ux, uy] = self.u[y][x];
        let u_mag = (ux * ux + uy * uy).sqrt();
        let cs = self.cs_field[y][x];
        u_mag / cs.max(1e-30)
    }

    /// Max Mach number in domain.
    pub fn max_mach(&self) -> f64 {
        let ny = self.rho.len();
        let nx = if ny > 0 { self.rho[0].len() } else { 0 };
        let mut ma_max = 0.0f64;
        for y in 0..ny {
            for x in 0..nx {
                ma_max = ma_max.max(self.mach_at(x, y));
            }
        }
        ma_max
    }

    /// Pressure ratio p/p_ref.
    pub fn pressure_ratio(&self, x: usize, y: usize) -> f64 {
        self.p[y][x] / self.p_ref
    }

    /// Density ratio rho/rho_ref.
    pub fn density_ratio(&self, x: usize, y: usize) -> f64 {
        self.rho[y][x] / self.rho_ref
    }

    /// Average entropy production (simplified).
    pub fn avg_entropy_production(&self) -> f64 {
        let ny = self.rho.len();
        let nx = if ny > 0 { self.rho[0].len() } else { 0 };
        let mut sum = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let s =
                    (self.p[y][x] / self.p_ref) / (self.rho[y][x] / self.rho_ref).powf(self.gamma);
                sum += s.ln().abs();
            }
        }
        sum / (nx * ny).max(1) as f64
    }
}

/// Detonation wave LBM simulation.
///
/// Models detonation wave with Chapman-Jouguet conditions and ZND structure.
#[derive(Debug, Clone)]
pub struct DetonationLbm {
    /// Number of cells.
    pub n: usize,
    /// Density profile.
    pub rho: Vec<f64>,
    /// Velocity profile.
    pub vel: Vec<f64>,
    /// Pressure profile.
    pub pressure: Vec<f64>,
    /// Reaction progress variable (0=unburnt, 1=burnt).
    pub progress: Vec<f64>,
    /// Heat release Q.
    pub heat_release: f64,
    /// Activation energy Ea.
    pub activation_energy: f64,
    /// Pre-exponential factor A.
    pub reaction_rate_a: f64,
    /// Gamma.
    pub gamma: f64,
    /// CJ detonation velocity.
    pub dcj: f64,
    /// Grid spacing.
    pub dx: f64,
    /// Time step.
    pub dt: f64,
}

impl DetonationLbm {
    /// Create a new detonation simulation.
    pub fn new(
        n: usize,
        dx: f64,
        dt: f64,
        gamma: f64,
        heat_release: f64,
        activation_energy: f64,
        reaction_rate_a: f64,
    ) -> Self {
        let rho = vec![1.0f64; n];
        let vel = vec![0.0f64; n];
        let pressure = vec![1.0f64; n];
        let mut progress = vec![0.0f64; n];
        // Initiate detonation at left boundary
        progress[0] = 1.0;
        progress[1] = 0.8;
        let dcj = Self::cj_velocity(gamma, heat_release);
        Self {
            n,
            rho,
            vel,
            pressure,
            progress,
            heat_release,
            activation_energy,
            reaction_rate_a,
            gamma,
            dcj,
            dx,
            dt,
        }
    }

    /// Chapman-Jouguet detonation velocity.
    pub fn cj_velocity(gamma: f64, q: f64) -> f64 {
        // DCJ = sqrt(2*(gamma^2-1)*Q) approximately
        (2.0 * (gamma * gamma - 1.0) * q).sqrt()
    }

    /// CJ pressure ratio.
    pub fn cj_pressure_ratio(gamma: f64) -> f64 {
        // p2/p1 = (2*gamma*M_cj^2 - (gamma-1)) / (gamma+1)
        (2.0 * gamma - (gamma - 1.0)) / (gamma + 1.0)
    }

    /// Arrhenius reaction rate.
    pub fn reaction_rate(&self, temp: f64, progress: f64) -> f64 {
        if progress >= 1.0 {
            return 0.0;
        }
        self.reaction_rate_a * (-self.activation_energy / temp.max(0.01)).exp() * (1.0 - progress)
    }

    /// Advance detonation one time step.
    pub fn step(&mut self) {
        let dt = self.dt;
        let dx = self.dx;
        let g = self.gamma;
        let q = self.heat_release;

        // Update reaction progress
        for i in 0..self.n {
            let temp = self.pressure[i] / (self.rho[i] + 1e-10) * g;
            let rate = self.reaction_rate(temp, self.progress[i]);
            self.progress[i] = (self.progress[i] + rate * dt).min(1.0);
        }

        // Simple advection of detonation front
        let vel_det = self.dcj;
        for i in 1..self.n {
            let rho_l = self.rho[i - 1];
            let p_l = self.pressure[i - 1];
            let flux = rho_l * vel_det;
            let heat = q * (self.progress[i - 1] - self.progress[i]).max(0.0);
            self.rho[i] = (self.rho[i] - dt / dx * (flux - rho_l * vel_det)).max(0.1);
            self.pressure[i] = (self.pressure[i] + heat * (g - 1.0) * dt).max(0.0);
            let _ = p_l;
        }
    }

    /// Induction zone length estimate.
    pub fn induction_length(&self) -> f64 {
        let dcj = self.dcj;
        let ea = self.activation_energy;
        let a = self.reaction_rate_a;
        dcj / (a * (-ea).exp())
    }
}

/// Roe flux splitting for Euler equations.
///
/// Returns approximate Riemann flux \[F_rho, F_rhou, F_rhov, F_E\].
pub fn euler_flux_roe(
    rho_l: f64,
    u_l: [f64; 2],
    p_l: f64,
    rho_r: f64,
    u_r: [f64; 2],
    p_r: f64,
) -> [f64; 4] {
    let gamma = 1.4;
    let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * (u_l[0] * u_l[0] + u_l[1] * u_l[1]);
    let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * (u_r[0] * u_r[0] + u_r[1] * u_r[1]);
    let h_l = (e_l + p_l) / rho_l;
    let h_r = (e_r + p_r) / rho_r;

    // Roe averages
    let sqrt_rho_l = rho_l.sqrt();
    let sqrt_rho_r = rho_r.sqrt();
    let denom = sqrt_rho_l + sqrt_rho_r;
    let u_roe = (sqrt_rho_l * u_l[0] + sqrt_rho_r * u_r[0]) / denom;
    let v_roe = (sqrt_rho_l * u_l[1] + sqrt_rho_r * u_r[1]) / denom;
    let h_roe = (sqrt_rho_l * h_l + sqrt_rho_r * h_r) / denom;
    let cs_roe = ((gamma - 1.0) * (h_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe)))
        .max(0.0)
        .sqrt();

    // Average flux
    let f_l = [
        rho_l * u_l[0],
        rho_l * u_l[0] * u_l[0] + p_l,
        rho_l * u_l[0] * u_l[1],
        (e_l + p_l) * u_l[0],
    ];
    let f_r = [
        rho_r * u_r[0],
        rho_r * u_r[0] * u_r[0] + p_r,
        rho_r * u_r[0] * u_r[1],
        (e_r + p_r) * u_r[0],
    ];

    // Lax-Friedrichs upwinding (simplified Roe)
    let lambda = u_roe.abs() + cs_roe;
    [
        0.5 * (f_l[0] + f_r[0]) - 0.5 * lambda * (rho_r - rho_l),
        0.5 * (f_l[1] + f_r[1]) - 0.5 * lambda * (rho_r * u_r[0] - rho_l * u_l[0]),
        0.5 * (f_l[2] + f_r[2]) - 0.5 * lambda * (rho_r * u_r[1] - rho_l * u_l[1]),
        0.5 * (f_l[3] + f_r[3]) - 0.5 * lambda * (e_r - e_l),
    ]
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compressible_lbm_creation() {
        let sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        assert_eq!(sim.nx, 8);
        assert_eq!(sim.ny, 8);
    }

    #[test]
    fn test_compressible_lbm_step() {
        let mut sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
                assert!(sim.rho[y][x].is_finite());
            }
        }
    }

    #[test]
    fn test_equilibrium_compressible_sum_rho() {
        let feq = equilibrium_compressible(1.0, 0.1, 0.05, 1.0 / 3.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_local_mach_zero_velocity() {
        let ma = local_mach(0.0, 0.0, 1.0);
        assert!(ma.abs() < 1e-15);
    }

    #[test]
    fn test_local_mach_unit() {
        let ma = local_mach(1.0, 0.0, 1.0);
        assert!((ma - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_entropy_generation_zero_shear() {
        let s = entropy_generation(1.0, 1.0, 0.0, 0.0, 1e-3);
        assert!(s.abs() < 1e-15);
    }

    #[test]
    fn test_euler_lbm_creation() {
        let sim = EulerLbm::new(8, 8, 1.0 / 3.0_f64.sqrt());
        assert_eq!(sim.nx, 8);
    }

    #[test]
    fn test_euler_lbm_step() {
        let mut sim = EulerLbm::new(8, 8, 1.0 / 3.0_f64.sqrt());
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
            }
        }
    }

    #[test]
    fn test_shock_capture_creation() {
        let sim = ShockCaptureLbm::new(8, 8, 1.7, 1.4, 0.01, 0.1);
        assert_eq!(sim.lbm.nx, 8);
    }

    #[test]
    fn test_shock_capture_step() {
        let mut sim = ShockCaptureLbm::new(8, 8, 1.7, 1.4, 0.01, 0.1);
        sim.step();
        assert!(sim.shock_cell_count() <= sim.lbm.nx * sim.lbm.ny);
    }

    #[test]
    fn test_ns_compressible_step() {
        let mut sim = NsCompressibleLbm::new(8, 8, 1.7, 1.7, 0.71, 1.4);
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
            }
        }
    }

    #[test]
    fn test_ma_correction_creation() {
        let mc = MaCorrection::new(2, 1.0 / 3.0_f64.sqrt());
        assert_eq!(mc.order, 2);
    }

    #[test]
    fn test_ma_correction_feq_sum() {
        let mc = MaCorrection::new(2, 1.0 / 3.0_f64.sqrt());
        let feq = mc.corrected_feq(1.0, 0.05, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_shock_wave_1d_creation() {
        let sw = ShockWave1D::new_sod(100, 0.01, 1e-4, 1.4);
        assert_eq!(sw.n, 100);
    }

    #[test]
    fn test_shock_wave_1d_step() {
        let mut sw = ShockWave1D::new_sod(100, 0.01, 1e-5, 1.4);
        sw.run(10);
        assert!(sw.time > 0.0);
    }

    #[test]
    fn test_shock_wave_cfl() {
        let sw = ShockWave1D::new_sod(100, 0.01, 1e-5, 1.4);
        let cfl = sw.cfl();
        assert!(cfl > 0.0);
    }

    #[test]
    fn test_rayleigh_taylor_creation() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 9.81, 0.001);
        assert_eq!(rt.nx, 16);
    }

    #[test]
    fn test_rayleigh_taylor_growth_rate() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 9.81, 0.001);
        let gr = rt.growth_rate(1.0);
        assert!(gr > 0.0);
    }

    #[test]
    fn test_supersonic_inlet_velocity() {
        let inlet = SupersonicInlet::new(2.0, 0.0, 1.0, 1.0 / 1.4, 1.4);
        let v = inlet.velocity();
        assert!(v[0] > 0.0);
        assert!(v[1].abs() < 1e-10);
    }

    #[test]
    fn test_supersonic_inlet_total_pressure() {
        let inlet = SupersonicInlet::new(2.0, 0.0, 1.0, 1.0, 1.4);
        let pt = inlet.total_pressure();
        assert!(pt > inlet.p_in); // total > static
    }

    #[test]
    fn test_compressible_stats_mach() {
        let n = 4;
        let rho = vec![vec![1.0f64; n]; n];
        let u = vec![vec![[0.5f64, 0.0f64]; n]; n];
        let p = vec![vec![1.0 / 3.0_f64; n]; n];
        let stats = CompressibleStats::new(rho, u, p, 1.4);
        let ma = stats.mach_at(0, 0);
        assert!(ma > 0.0);
    }

    #[test]
    fn test_detonation_cj_velocity() {
        let dcj = DetonationLbm::cj_velocity(1.4, 10.0);
        assert!(dcj > 0.0);
    }

    #[test]
    fn test_detonation_step() {
        let mut det = DetonationLbm::new(50, 0.01, 1e-5, 1.4, 10.0, 5.0, 1e6);
        det.step();
        for &rho in &det.rho {
            assert!(rho > 0.0);
        }
    }

    #[test]
    fn test_euler_flux_roe_symmetry() {
        let rho = 1.0;
        let u = [0.1f64, 0.0f64];
        let p = 1.0 / 3.0;
        let flux = euler_flux_roe(rho, u, p, rho, u, p);
        // Same state on both sides: flux should be exact physical flux
        let expected_rho_flux = rho * u[0];
        assert!((flux[0] - expected_rho_flux).abs() < 1e-10);
    }

    #[test]
    fn test_max_mach_compressible() {
        let mut sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        sim.u[0][0] = [0.1, 0.0];
        let ma = sim.max_mach();
        assert!(ma >= 0.0);
    }

    #[test]
    fn test_rt_bubble_velocity() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 9.81, 0.001);
        let vb = rt.bubble_velocity();
        assert!(vb > 0.0);
    }
}
