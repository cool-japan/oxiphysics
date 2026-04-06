#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cahn-Hilliard equation via Lattice Boltzmann Method.
//!
//! This module implements the Cahn-Hilliard (CH) equation for phase-field
//! modeling using the lattice Boltzmann method (LBM). It covers:
//!
//! - **Spinodal decomposition**: unstable phase separation from a mixed state
//! - **Interface dynamics**: diffuse interface tracking with Ginzburg-Landau free energy
//! - **Mobility functions**: constant, degenerate, and concentration-dependent
//! - **Wetting boundary conditions**: contact angle enforcement at solid walls
//! - **Ginzburg-Landau free energy**: double-well potential with gradient energy
//! - **Multi-relaxation-time (MRT)**: enhanced stability for the CH-LBM
//!
//! # Physical Model
//!
//! The Cahn-Hilliard equation:
//!
//! ```text
//! dC/dt = div( M(C) * grad(mu) )
//! mu = dF/dC = A * C * (C^2 - 1) - kappa * laplacian(C)
//! ```
//!
//! where C is the concentration (order parameter), mu is the chemical potential,
//! M(C) is the mobility, A is the bulk energy coefficient, and kappa controls
//! the interface width.
//!
//! # References
//!
//! - Zheng, Shu, Chew (2006), *Phys. Rev. E* **73**, 056702.
//! - Lee & Liu (2010), *J. Comput. Phys.* **229**, 8045-8063.
//! - Liang et al. (2014), *Phys. Rev. E* **89**, 053320.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 lattice constants
// ─────────────────────────────────────────────────────────────────────────────

/// D2Q9 lattice weights.
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

/// D2Q9 velocity x-components.
const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];

/// D2Q9 velocity y-components.
const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Opposite direction indices for bounce-back.
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Speed of sound squared for D2Q9: cs^2 = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ─────────────────────────────────────────────────────────────────────────────
// Mobility types
// ─────────────────────────────────────────────────────────────────────────────

/// Mobility function type for the Cahn-Hilliard equation.
#[derive(Debug, Clone, Copy)]
pub enum MobilityType {
    /// Constant mobility M(C) = M0.
    Constant,
    /// Degenerate mobility M(C) = M0 * C * (1 - C), vanishes at pure phases.
    Degenerate,
    /// Concentration-dependent mobility M(C) = M0 * (1 - C^2).
    ConcentrationDependent,
    /// Polynomial mobility M(C) = M0 * (1 - C^2)^2.
    Polynomial,
}

/// Evaluate the mobility function at a given concentration.
///
/// `c` is the order parameter, `m0` is the base mobility.
pub fn eval_mobility(mobility_type: MobilityType, c: f64, m0: f64) -> f64 {
    match mobility_type {
        MobilityType::Constant => m0,
        MobilityType::Degenerate => {
            let cc = c.clamp(0.0, 1.0);
            m0 * cc * (1.0 - cc)
        }
        MobilityType::ConcentrationDependent => m0 * (1.0 - c * c).max(0.0),
        MobilityType::Polynomial => {
            let t = (1.0 - c * c).max(0.0);
            m0 * t * t
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ginzburg-Landau free energy
// ─────────────────────────────────────────────────────────────────────────────

/// Ginzburg-Landau free energy density: f_bulk = A/4 * (C^2 - 1)^2.
///
/// `a_coeff` is the bulk energy coefficient A.
pub fn gl_bulk_energy_density(c: f64, a_coeff: f64) -> f64 {
    a_coeff / 4.0 * (c * c - 1.0).powi(2)
}

/// Derivative of bulk free energy: df/dC = A * C * (C^2 - 1).
pub fn gl_bulk_energy_derivative(c: f64, a_coeff: f64) -> f64 {
    a_coeff * c * (c * c - 1.0)
}

/// Second derivative of bulk free energy: d2f/dC2 = A * (3*C^2 - 1).
pub fn gl_bulk_energy_second_derivative(c: f64, a_coeff: f64) -> f64 {
    a_coeff * (3.0 * c * c - 1.0)
}

/// Interface width from the Ginzburg-Landau parameters.
///
/// The equilibrium interface width is W = sqrt(2 * kappa / A).
pub fn interface_width(kappa: f64, a_coeff: f64) -> f64 {
    (2.0 * kappa / a_coeff).sqrt()
}

/// Surface tension from the Ginzburg-Landau parameters.
///
/// sigma = sqrt(8 * kappa * A / 9).
pub fn surface_tension(kappa: f64, a_coeff: f64) -> f64 {
    (8.0 * kappa * a_coeff / 9.0).sqrt()
}

/// Equilibrium interface profile: C(x) = tanh(x / (sqrt(2) * W)).
///
/// `x` is the distance from the interface center, `width` is the interface width.
pub fn equilibrium_profile(x: f64, width: f64) -> f64 {
    (x / (std::f64::consts::SQRT_2 * width)).tanh()
}

/// Cahn number: ratio of interface width to system size.
pub fn cahn_number(kappa: f64, a_coeff: f64, system_size: f64) -> f64 {
    interface_width(kappa, a_coeff) / system_size
}

// ─────────────────────────────────────────────────────────────────────────────
// Wetting boundary conditions
// ─────────────────────────────────────────────────────────────────────────────

/// Contact angle boundary condition parameter.
///
/// For a wall at the bottom (y=0), the contact angle theta determines
/// the wetting condition: dC/dn = -sqrt(2*A/kappa) * cos(theta).
#[derive(Debug, Clone, Copy)]
pub struct WettingBC {
    /// Contact angle in radians. theta < pi/2 is hydrophilic, theta > pi/2 is hydrophobic.
    pub theta: f64,
    /// Bulk energy coefficient A.
    pub a_coeff: f64,
    /// Interface gradient energy coefficient kappa.
    pub kappa: f64,
}

impl WettingBC {
    /// Create a new wetting boundary condition.
    pub fn new(theta: f64, a_coeff: f64, kappa: f64) -> Self {
        WettingBC {
            theta,
            a_coeff,
            kappa,
        }
    }

    /// Normal gradient of C at the wall: dC/dn = -sqrt(2A/kappa) * cos(theta).
    pub fn normal_gradient(&self) -> f64 {
        -(2.0 * self.a_coeff / self.kappa).sqrt() * self.theta.cos()
    }

    /// Whether the wall is hydrophilic (theta < pi/2).
    pub fn is_hydrophilic(&self) -> bool {
        self.theta < std::f64::consts::FRAC_PI_2
    }

    /// Whether the wall is hydrophobic (theta > pi/2).
    pub fn is_hydrophobic(&self) -> bool {
        self.theta > std::f64::consts::FRAC_PI_2
    }

    /// Neutral wetting condition (theta = pi/2).
    pub fn is_neutral(&self) -> bool {
        (self.theta - std::f64::consts::FRAC_PI_2).abs() < 1e-10
    }

    /// Surface energy contribution from the wetting BC.
    ///
    /// E_wall = -sqrt(2*A*kappa) * cos(theta) * C_wall.
    pub fn wall_energy(&self, c_wall: f64) -> f64 {
        -(2.0 * self.a_coeff * self.kappa).sqrt() * self.theta.cos() * c_wall
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase-field analysis utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the volume fraction of phase +1 (C > 0).
pub fn volume_fraction_positive(phi: &[f64]) -> f64 {
    let n = phi.len();
    if n == 0 {
        return 0.0;
    }
    let count = phi.iter().filter(|&&c| c > 0.0).count();
    count as f64 / n as f64
}

/// Compute the mean concentration.
pub fn mean_concentration(phi: &[f64]) -> f64 {
    let n = phi.len();
    if n == 0 {
        return 0.0;
    }
    phi.iter().sum::<f64>() / n as f64
}

/// Compute the variance of the concentration field.
pub fn concentration_variance(phi: &[f64]) -> f64 {
    let n = phi.len() as f64;
    if n <= 0.0 {
        return 0.0;
    }
    let mean = phi.iter().sum::<f64>() / n;
    phi.iter().map(|&c| (c - mean).powi(2)).sum::<f64>() / n
}

/// Compute the structure factor S(k) from the concentration field.
///
/// Uses a simple 1D radial average. Returns (wavenumber, S(k)) pairs.
/// `phi` is the 2D field stored row-major, `nx` and `ny` are dimensions.
pub fn structure_factor_1d(phi: &[f64], nx: usize, ny: usize) -> Vec<(f64, f64)> {
    let max_k = nx.min(ny) / 2;
    if max_k == 0 {
        return vec![];
    }
    let n = (nx * ny) as f64;
    let mut sk = vec![0.0_f64; max_k];
    let mut counts = vec![0_usize; max_k];

    // Compute Fourier coefficients via direct DFT (simplified for small grids)
    for kx in 0..nx {
        for ky in 0..ny {
            if kx == 0 && ky == 0 {
                continue;
            }
            let kxf = if kx <= nx / 2 {
                kx as f64
            } else {
                (kx as f64) - (nx as f64)
            };
            let kyf = if ky <= ny / 2 {
                ky as f64
            } else {
                (ky as f64) - (ny as f64)
            };
            let kmag = (kxf * kxf + kyf * kyf).sqrt();
            let bin = kmag.round() as usize;
            if bin >= max_k || bin == 0 {
                continue;
            }
            // DFT coefficient
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            let twopi = 2.0 * std::f64::consts::PI;
            for ix in 0..nx {
                for iy in 0..ny {
                    let angle = twopi
                        * (kx as f64 * ix as f64 / nx as f64 + ky as f64 * iy as f64 / ny as f64);
                    let val = phi[ix * ny + iy];
                    re += val * angle.cos();
                    im -= val * angle.sin();
                }
            }
            sk[bin] += (re * re + im * im) / n;
            counts[bin] += 1;
        }
    }

    let mut result = Vec::new();
    for k in 1..max_k {
        if counts[k] > 0 {
            result.push((k as f64, sk[k] / counts[k] as f64));
        }
    }
    result
}

/// Estimate the dominant wavenumber from the structure factor.
pub fn dominant_wavenumber(sk: &[(f64, f64)]) -> f64 {
    if sk.is_empty() {
        return 0.0;
    }
    sk.iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|&(k, _)| k)
        .unwrap_or(0.0)
}

/// Compute the total Ginzburg-Landau free energy of the system.
///
/// F = sum \[ A/4*(C^2-1)^2 + kappa/2*|grad(C)|^2 \]
pub fn total_free_energy(phi: &[f64], nx: usize, ny: usize, a_coeff: f64, kappa: f64) -> f64 {
    let mut fe = 0.0_f64;
    for i in 0..nx {
        for j in 0..ny {
            let idx = i * ny + j;
            let c = phi[idx];
            // Bulk
            fe += a_coeff / 4.0 * (c * c - 1.0).powi(2);
            // Gradient energy (central differences, periodic)
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let jp = (j + 1) % ny;
            let jm = (j + ny - 1) % ny;
            let dcdx = (phi[ip * ny + j] - phi[im * ny + j]) / 2.0;
            let dcdy = (phi[i * ny + jp] - phi[i * ny + jm]) / 2.0;
            fe += kappa / 2.0 * (dcdx * dcdx + dcdy * dcdy);
        }
    }
    fe
}

/// Compute the interfacial area (perimeter) from the concentration field.
///
/// Uses |grad(C)| integrated over the domain as a proxy for interface area.
pub fn interfacial_area(phi: &[f64], nx: usize, ny: usize) -> f64 {
    let mut area = 0.0_f64;
    for i in 0..nx {
        for j in 0..ny {
            let ip = (i + 1) % nx;
            let im = (i + nx - 1) % nx;
            let jp = (j + 1) % ny;
            let jm = (j + ny - 1) % ny;
            let dcdx = (phi[ip * ny + j] - phi[im * ny + j]) / 2.0;
            let dcdy = (phi[i * ny + jp] - phi[i * ny + jm]) / 2.0;
            area += (dcdx * dcdx + dcdy * dcdy).sqrt();
        }
    }
    area
}

// ─────────────────────────────────────────────────────────────────────────────
// Cahn-Hilliard LBM solver
// ─────────────────────────────────────────────────────────────────────────────

/// Cahn-Hilliard LBM solver with phase-field modeling on a D2Q9 lattice.
///
/// Implements the Cahn-Hilliard equation using a single distribution function
/// approach with BGK collision and chemical potential source terms.
#[derive(Debug, Clone)]
pub struct CahnHilliardLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Concentration (order parameter) field C in \[-1, 1\].
    pub concentration: Vec<f64>,
    /// Chemical potential field mu.
    pub chemical_potential: Vec<f64>,
    /// Distribution functions (nx*ny*9).
    pub f_dist: Vec<f64>,
    /// Auxiliary distribution for chemical potential transport (nx*ny*9).
    pub g_dist: Vec<f64>,
    /// Bulk energy coefficient A.
    pub a_coeff: f64,
    /// Interface gradient energy coefficient kappa.
    pub kappa: f64,
    /// Base mobility M0.
    pub mobility: f64,
    /// Mobility type.
    pub mobility_type: MobilityType,
    /// Relaxation time for the concentration equation.
    pub tau_c: f64,
    /// Relaxation time for the chemical potential equation.
    pub tau_mu: f64,
    /// Wetting boundary conditions (bottom, top, left, right). None = periodic.
    pub wetting_bottom: Option<WettingBC>,
    /// Top wall wetting BC.
    pub wetting_top: Option<WettingBC>,
    /// Left wall wetting BC.
    pub wetting_left: Option<WettingBC>,
    /// Right wall wetting BC.
    pub wetting_right: Option<WettingBC>,
    /// Current simulation time step.
    pub time_step: usize,
}

impl CahnHilliardLbm {
    /// Create a new Cahn-Hilliard LBM solver.
    ///
    /// # Arguments
    /// - `nx`, `ny`: grid dimensions
    /// - `a_coeff`: bulk energy coefficient A
    /// - `kappa`: interface gradient energy coefficient
    /// - `mobility`: base mobility M0
    /// - `mobility_type`: type of mobility function
    pub fn new(
        nx: usize,
        ny: usize,
        a_coeff: f64,
        kappa: f64,
        mobility: f64,
        mobility_type: MobilityType,
    ) -> Self {
        let n = nx * ny;
        let tau_c = 0.5 + 3.0 * mobility;
        let tau_mu = 1.0;
        CahnHilliardLbm {
            nx,
            ny,
            concentration: vec![0.0; n],
            chemical_potential: vec![0.0; n],
            f_dist: vec![0.0; n * 9],
            g_dist: vec![0.0; n * 9],
            a_coeff,
            kappa,
            mobility,
            mobility_type,
            tau_c,
            tau_mu,
            wetting_bottom: None,
            wetting_top: None,
            wetting_left: None,
            wetting_right: None,
            time_step: 0,
        }
    }

    /// Flat index for node (i, j).
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }

    /// Distribution function index for node (i, j), direction q.
    #[inline]
    pub fn fidx(&self, i: usize, j: usize, q: usize) -> usize {
        (i * self.ny + j) * 9 + q
    }

    /// Set a wetting boundary condition on the bottom wall.
    pub fn set_wetting_bottom(&mut self, theta: f64) {
        self.wetting_bottom = Some(WettingBC::new(theta, self.a_coeff, self.kappa));
    }

    /// Set a wetting boundary condition on the top wall.
    pub fn set_wetting_top(&mut self, theta: f64) {
        self.wetting_top = Some(WettingBC::new(theta, self.a_coeff, self.kappa));
    }

    /// Set a wetting boundary condition on the left wall.
    pub fn set_wetting_left(&mut self, theta: f64) {
        self.wetting_left = Some(WettingBC::new(theta, self.a_coeff, self.kappa));
    }

    /// Set a wetting boundary condition on the right wall.
    pub fn set_wetting_right(&mut self, theta: f64) {
        self.wetting_right = Some(WettingBC::new(theta, self.a_coeff, self.kappa));
    }

    /// Initialize with random perturbation around mean `c0`.
    pub fn init_random(&mut self, c0: f64, noise: f64) {
        use rand::RngExt;
        let mut rng = rand::rng();
        for c in self.concentration.iter_mut() {
            *c = c0 + noise * rng.random_range(-1.0_f64..1.0_f64);
        }
        self.compute_chemical_potential();
        self.init_equilibrium();
    }

    /// Initialize with a flat interface along the x-direction at position y0.
    pub fn init_flat_interface(&mut self, y0: f64, width: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let y = j as f64 - y0;
                let idx = i * ny + j;
                self.concentration[idx] = equilibrium_profile(y, width);
            }
        }
        self.compute_chemical_potential();
        self.init_equilibrium();
    }

    /// Initialize with a circular droplet centered at (cx, cy) with radius r.
    pub fn init_droplet(&mut self, cx: f64, cy: f64, radius: f64, width: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let dx = i as f64 - cx;
                let dy = j as f64 - cy;
                let dist = (dx * dx + dy * dy).sqrt() - radius;
                let idx = i * ny + j;
                self.concentration[idx] = -equilibrium_profile(dist, width);
            }
        }
        self.compute_chemical_potential();
        self.init_equilibrium();
    }

    /// Initialize with a striped pattern (cosine wave along x).
    pub fn init_stripes(&mut self, n_stripes: usize) {
        let nx = self.nx;
        let ny = self.ny;
        let twopi = 2.0 * std::f64::consts::PI;
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                self.concentration[idx] = (twopi * n_stripes as f64 * i as f64 / nx as f64).cos();
            }
        }
        self.compute_chemical_potential();
        self.init_equilibrium();
    }

    /// Compute the chemical potential mu = A*C*(C^2-1) - kappa*laplacian(C).
    pub fn compute_chemical_potential(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let a = self.a_coeff;
        let kappa = self.kappa;
        let mut mu_new = vec![0.0_f64; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.idx(i, j);
                let c = self.concentration[idx];
                let lap = self.laplacian_c(i, j);
                mu_new[idx] = a * c * (c * c - 1.0) - kappa * lap;
            }
        }
        // Apply wetting BC modifications
        if let Some(ref wbc) = self.wetting_bottom {
            let grad_n = wbc.normal_gradient();
            for i in 0..nx {
                let idx = self.idx(i, 0);
                mu_new[idx] += kappa * grad_n; // Wall correction
            }
        }
        if let Some(ref wbc) = self.wetting_top {
            let grad_n = wbc.normal_gradient();
            for i in 0..nx {
                let idx = self.idx(i, ny - 1);
                mu_new[idx] -= kappa * grad_n; // Opposite normal
            }
        }
        self.chemical_potential = mu_new;
    }

    /// D2Q9 isotropic Laplacian of C at (i, j) with periodic BC.
    pub fn laplacian_c(&self, i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let c_center = self.concentration[self.idx(i, j)];
        let mut lap = 0.0_f64;
        for q in 1..9 {
            let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + CY[q]).rem_euclid(ny as i32) as usize;
            let c_neigh = self.concentration[self.idx(ni, nj)];
            lap += W9[q] * (c_neigh - c_center);
        }
        lap * 6.0
    }

    /// D2Q9 isotropic Laplacian of the chemical potential at (i, j).
    fn laplacian_mu(&self, i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mu_c = self.chemical_potential[self.idx(i, j)];
        let mut lap = 0.0_f64;
        for q in 1..9 {
            let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + CY[q]).rem_euclid(ny as i32) as usize;
            lap += W9[q] * (self.chemical_potential[self.idx(ni, nj)] - mu_c);
        }
        lap * 6.0
    }

    /// Gradient of C at (i, j) using the D2Q9 stencil.
    pub fn gradient_c(&self, i: usize, j: usize) -> [f64; 2] {
        let nx = self.nx;
        let ny = self.ny;
        let mut gx = 0.0_f64;
        let mut gy = 0.0_f64;
        for q in 1..9 {
            let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + CY[q]).rem_euclid(ny as i32) as usize;
            let c_neigh = self.concentration[self.idx(ni, nj)];
            gx += W9[q] * CX[q] as f64 * c_neigh;
            gy += W9[q] * CY[q] as f64 * c_neigh;
        }
        [gx * 3.0, gy * 3.0]
    }

    /// Gradient of mu at (i, j) using the D2Q9 stencil.
    pub fn gradient_mu(&self, i: usize, j: usize) -> [f64; 2] {
        let nx = self.nx;
        let ny = self.ny;
        let mut gx = 0.0_f64;
        let mut gy = 0.0_f64;
        for q in 1..9 {
            let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + CY[q]).rem_euclid(ny as i32) as usize;
            let mu_neigh = self.chemical_potential[self.idx(ni, nj)];
            gx += W9[q] * CX[q] as f64 * mu_neigh;
            gy += W9[q] * CY[q] as f64 * mu_neigh;
        }
        [gx * 3.0, gy * 3.0]
    }

    /// Equilibrium distribution for the concentration field.
    fn feq(&self, c: f64, q: usize) -> f64 {
        W9[q] * c
    }

    /// Equilibrium distribution for the chemical potential field.
    fn geq(&self, mu: f64, q: usize) -> f64 {
        W9[q] * mu
    }

    /// Initialize distributions to equilibrium.
    fn init_equilibrium(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.idx(i, j);
                let c = self.concentration[idx];
                let mu = self.chemical_potential[idx];
                for q in 0..9 {
                    let fidx = self.fidx(i, j, q);
                    self.f_dist[fidx] = self.feq(c, q);
                    self.g_dist[fidx] = self.geq(mu, q);
                }
            }
        }
    }

    /// Perform one full Cahn-Hilliard LBM time step.
    ///
    /// Steps:
    /// 1. Compute chemical potential from current concentration
    /// 2. Update concentration via mobility-weighted Laplacian of mu
    /// 3. Recompute chemical potential
    /// 4. Re-initialize distributions to equilibrium
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // 1. Compute chemical potential
        self.compute_chemical_potential();

        // 2. Update concentration: C += M(C) * laplacian(mu)
        let mut dc = vec![0.0_f64; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.idx(i, j);
                let c = self.concentration[idx];
                let m = eval_mobility(self.mobility_type, c, self.mobility);
                dc[idx] = m * self.laplacian_mu(i, j);
            }
        }
        for (idx, &d) in dc.iter().enumerate() {
            self.concentration[idx] += d;
        }

        // 3. Recompute chemical potential
        self.compute_chemical_potential();

        // 4. Re-initialize distributions
        self.init_equilibrium();

        self.time_step += 1;
    }

    /// Perform multiple time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// BGK collision step for both distributions.
    ///
    /// Uses the relaxation-time formulation with source terms.
    pub fn collide(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega_c = 1.0 / self.tau_c;
        let omega_mu = 1.0 / self.tau_mu;

        for i in 0..nx {
            for j in 0..ny {
                let idx = self.idx(i, j);
                let c = self.concentration[idx];
                let mu = self.chemical_potential[idx];
                let m = eval_mobility(self.mobility_type, c, self.mobility);

                for q in 0..9 {
                    let base = idx * 9 + q;
                    let feq = self.feq(c, q);
                    let geq = self.geq(mu, q);

                    let source = if q == 0 {
                        -m * (1.0 - W9[0]) * mu * omega_c
                    } else {
                        m * W9[q] * mu * omega_c
                    };

                    self.f_dist[base] += omega_c * (feq - self.f_dist[base]) + source;
                    self.g_dist[base] += omega_mu * (geq - self.g_dist[base]);
                }
            }
        }
    }

    /// Streaming step: propagate populations to neighbours.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny * 9;
        let mut f_new = vec![0.0_f64; n];
        let mut g_new = vec![0.0_f64; n];

        for i in 0..nx {
            for j in 0..ny {
                for q in 0..9 {
                    let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
                    let nj = (j as i32 + CY[q]).rem_euclid(ny as i32) as usize;
                    let src = (i * ny + j) * 9 + q;
                    let dst = (ni * ny + nj) * 9 + q;
                    f_new[dst] = self.f_dist[src];
                    g_new[dst] = self.g_dist[src];
                }
            }
        }
        self.f_dist = f_new;
        self.g_dist = g_new;
    }

    /// Apply wetting boundary conditions at walls.
    pub fn apply_wetting_bc(&mut self) {
        let nx = self.nx;
        let ny = self.ny;

        // Bottom wall (j = 0): bounce-back with wetting
        if let Some(ref wbc) = self.wetting_bottom {
            let grad_n = wbc.normal_gradient();
            for i in 0..nx {
                let idx = self.idx(i, 0);
                // Adjust ghost node concentration
                self.concentration[idx] = self.concentration[self.idx(i, 1)] + grad_n;
            }
        }

        // Top wall (j = ny-1)
        if let Some(ref wbc) = self.wetting_top {
            let grad_n = wbc.normal_gradient();
            for i in 0..nx {
                let idx = self.idx(i, ny - 1);
                self.concentration[idx] = self.concentration[self.idx(i, ny - 2)] - grad_n;
            }
        }
    }

    /// Update the concentration from the zeroth moment of f_dist.
    pub fn update_concentration_from_dist(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let base = (i * ny + j) * 9;
                let mut sum = 0.0_f64;
                for q in 0..9 {
                    sum += self.f_dist[base + q];
                }
                self.concentration[i * ny + j] = sum;
            }
        }
    }

    /// Total concentration (conservation diagnostic).
    pub fn total_concentration(&self) -> f64 {
        self.concentration.iter().sum()
    }

    /// Concentration variance (spinodal decomposition diagnostic).
    pub fn variance(&self) -> f64 {
        concentration_variance(&self.concentration)
    }

    /// Total Ginzburg-Landau free energy.
    pub fn free_energy(&self) -> f64 {
        total_free_energy(
            &self.concentration,
            self.nx,
            self.ny,
            self.a_coeff,
            self.kappa,
        )
    }

    /// Interfacial area (perimeter).
    pub fn interface_area(&self) -> f64 {
        interfacial_area(&self.concentration, self.nx, self.ny)
    }

    /// Interface width parameter.
    pub fn interface_width(&self) -> f64 {
        interface_width(self.kappa, self.a_coeff)
    }

    /// Surface tension parameter.
    pub fn surface_tension(&self) -> f64 {
        surface_tension(self.kappa, self.a_coeff)
    }

    /// Cahn number Cn = W / L.
    pub fn cahn_number(&self) -> f64 {
        cahn_number(self.kappa, self.a_coeff, self.nx.max(self.ny) as f64)
    }

    /// Volume fraction of the +1 phase.
    pub fn volume_fraction_positive(&self) -> f64 {
        volume_fraction_positive(&self.concentration)
    }

    /// Mean concentration.
    pub fn mean_concentration(&self) -> f64 {
        mean_concentration(&self.concentration)
    }

    /// Maximum concentration value.
    pub fn max_concentration(&self) -> f64 {
        self.concentration
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum concentration value.
    pub fn min_concentration(&self) -> f64 {
        self.concentration
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MRT-based Cahn-Hilliard LBM
// ─────────────────────────────────────────────────────────────────────────────

/// MRT collision matrix for enhanced stability in the CH-LBM.
///
/// Returns the transformation matrix M for D2Q9 (row-major, 9x9).
pub fn mrt_transform_matrix() -> [f64; 81] {
    // Standard D2Q9 MRT transformation matrix
    let mut m = [0.0_f64; 81];
    // Row 0: density (rho)
    for q in 0..9 {
        m[q] = 1.0;
    }
    // Row 1: energy (e)
    let e_row = [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0];
    m[9..18].copy_from_slice(&e_row);
    // Row 2: energy square (eps)
    let eps_row = [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0];
    m[18..27].copy_from_slice(&eps_row);
    // Row 3: x-momentum (jx)
    for q in 0..9 {
        m[27 + q] = CX[q] as f64;
    }
    // Row 4: energy flux qx
    let qx_row = [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    m[36..45].copy_from_slice(&qx_row);
    // Row 5: y-momentum (jy)
    for q in 0..9 {
        m[45 + q] = CY[q] as f64;
    }
    // Row 6: energy flux qy
    let qy_row = [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0];
    m[54..63].copy_from_slice(&qy_row);
    // Row 7: stress tensor pxx
    for q in 0..9 {
        m[63 + q] = (CX[q] * CX[q] - CY[q] * CY[q]) as f64;
    }
    // Row 8: stress tensor pxy
    for q in 0..9 {
        m[72 + q] = (CX[q] * CY[q]) as f64;
    }
    m
}

/// Compute the relaxation rates for MRT-CH-LBM.
///
/// Returns a 9-element array of relaxation rates for each moment.
pub fn mrt_relaxation_rates(tau: f64) -> [f64; 9] {
    let omega = 1.0 / tau;
    // s0 = 0 (conserved), others tuned for stability
    [0.0, 1.1, 1.1, 0.0, 1.2, 0.0, 1.2, omega, omega]
}

// ─────────────────────────────────────────────────────────────────────────────
// Spinodal decomposition analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Linear stability analysis of the Cahn-Hilliard equation.
///
/// The growth rate sigma(k) for a perturbation of wavenumber k around C=c0:
/// sigma(k) = -M * k^2 * (A*(3*c0^2-1) + kappa*k^2).
///
/// Modes with sigma > 0 are unstable (spinodal decomposition).
pub fn linear_growth_rate(k: f64, c0: f64, mobility: f64, a_coeff: f64, kappa: f64) -> f64 {
    -mobility * k * k * (a_coeff * (3.0 * c0 * c0 - 1.0) + kappa * k * k)
}

/// Critical wavenumber for spinodal decomposition.
///
/// k_c = sqrt(A*(1-3*c0^2) / kappa) (exists only inside spinodal region).
pub fn critical_wavenumber(c0: f64, a_coeff: f64, kappa: f64) -> f64 {
    let arg = a_coeff * (1.0 - 3.0 * c0 * c0) / kappa;
    if arg <= 0.0 { 0.0 } else { arg.sqrt() }
}

/// Most unstable wavenumber (fastest growing mode).
///
/// k_max = k_c / sqrt(2).
pub fn fastest_growing_wavenumber(c0: f64, a_coeff: f64, kappa: f64) -> f64 {
    critical_wavenumber(c0, a_coeff, kappa) / std::f64::consts::SQRT_2
}

/// Maximum growth rate at the fastest growing mode.
pub fn max_growth_rate(c0: f64, mobility: f64, a_coeff: f64, kappa: f64) -> f64 {
    let kmax = fastest_growing_wavenumber(c0, a_coeff, kappa);
    if kmax <= 0.0 {
        return 0.0;
    }
    linear_growth_rate(kmax, c0, mobility, a_coeff, kappa)
}

/// Spinodal region: the range of C where the mixture is unstable.
///
/// Returns (c_lower, c_upper) where spinodal decomposition occurs.
/// For f = A/4*(C^2-1)^2, the spinodal is |C| < 1/sqrt(3).
pub fn spinodal_region(a_coeff: f64) -> (f64, f64) {
    let _a = a_coeff; // Parameter determines potential shape
    let c_sp = 1.0 / 3.0_f64.sqrt();
    (-c_sp, c_sp)
}

/// Coarsening time scale estimate: t ~ L^3 / (M * sigma).
///
/// `length_scale` is the characteristic domain size.
pub fn coarsening_timescale(length_scale: f64, mobility: f64, kappa: f64, a_coeff: f64) -> f64 {
    let sigma = surface_tension(kappa, a_coeff);
    if mobility * sigma <= 0.0 {
        return f64::INFINITY;
    }
    length_scale.powi(3) / (mobility * sigma)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_solver(nx: usize, ny: usize) -> CahnHilliardLbm {
        CahnHilliardLbm::new(nx, ny, 1.0, 0.5, 0.1, MobilityType::Constant)
    }

    // -- basic construction --

    #[test]
    fn test_new_sizes() {
        let s = make_solver(8, 8);
        assert_eq!(s.concentration.len(), 64);
        assert_eq!(s.f_dist.len(), 64 * 9);
        assert_eq!(s.g_dist.len(), 64 * 9);
    }

    #[test]
    fn test_parameters_stored() {
        let s = CahnHilliardLbm::new(4, 4, 2.0, 0.7, 0.3, MobilityType::Degenerate);
        assert!((s.a_coeff - 2.0).abs() < 1e-14);
        assert!((s.kappa - 0.7).abs() < 1e-14);
        assert!((s.mobility - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_idx_correctness() {
        let s = make_solver(4, 6);
        assert_eq!(s.idx(0, 0), 0);
        assert_eq!(s.idx(1, 0), 6);
        assert_eq!(s.idx(0, 1), 1);
        assert_eq!(s.idx(3, 5), 23);
    }

    #[test]
    fn test_fidx_correctness() {
        let s = make_solver(4, 4);
        assert_eq!(s.fidx(0, 0, 0), 0);
        assert_eq!(s.fidx(0, 0, 8), 8);
        assert_eq!(s.fidx(1, 0, 0), 36);
    }

    // -- lattice constants --

    #[test]
    fn test_weights_sum_to_one() {
        let s: f64 = W9.iter().sum();
        assert!((s - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_opposite_directions() {
        for q in 0..9 {
            let oq = OPP[q];
            assert_eq!(CX[q] + CX[oq], 0);
            assert_eq!(CY[q] + CY[oq], 0);
        }
    }

    // -- mobility functions --

    #[test]
    fn test_constant_mobility() {
        let m = eval_mobility(MobilityType::Constant, 0.5, 0.1);
        assert!((m - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_degenerate_mobility_at_boundaries() {
        let m0 = eval_mobility(MobilityType::Degenerate, 0.0, 1.0);
        let m1 = eval_mobility(MobilityType::Degenerate, 1.0, 1.0);
        assert!(m0.abs() < 1e-14);
        assert!(m1.abs() < 1e-14);
    }

    #[test]
    fn test_concentration_dependent_mobility() {
        let m = eval_mobility(MobilityType::ConcentrationDependent, 0.0, 1.0);
        assert!((m - 1.0).abs() < 1e-14);
        let m2 = eval_mobility(MobilityType::ConcentrationDependent, 1.0, 1.0);
        assert!(m2.abs() < 1e-14);
    }

    // -- free energy --

    #[test]
    fn test_gl_bulk_energy_minima() {
        // Minima at C = +1 and C = -1
        let f1 = gl_bulk_energy_density(1.0, 1.0);
        let fm1 = gl_bulk_energy_density(-1.0, 1.0);
        assert!(f1.abs() < 1e-14);
        assert!(fm1.abs() < 1e-14);
    }

    #[test]
    fn test_gl_bulk_energy_maximum_at_zero() {
        let f0 = gl_bulk_energy_density(0.0, 1.0);
        assert!((f0 - 0.25).abs() < 1e-14);
    }

    #[test]
    fn test_gl_derivative_at_minima() {
        let d1 = gl_bulk_energy_derivative(1.0, 1.0);
        let dm1 = gl_bulk_energy_derivative(-1.0, 1.0);
        assert!(d1.abs() < 1e-14);
        assert!(dm1.abs() < 1e-14);
    }

    #[test]
    fn test_interface_width_formula() {
        let w = interface_width(0.5, 1.0);
        assert!((w - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_surface_tension_formula() {
        let sigma = surface_tension(0.5, 1.0);
        let expected = (8.0_f64 * 0.5 * 1.0 / 9.0).sqrt();
        assert!((sigma - expected).abs() < 1e-10);
    }

    // -- equilibrium profile --

    #[test]
    fn test_equilibrium_profile_zero() {
        let p = equilibrium_profile(0.0, 1.0);
        assert!(p.abs() < 1e-14);
    }

    #[test]
    fn test_equilibrium_profile_far_positive() {
        let p = equilibrium_profile(100.0, 1.0);
        assert!((p - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_equilibrium_profile_far_negative() {
        let p = equilibrium_profile(-100.0, 1.0);
        assert!((p + 1.0).abs() < 1e-6);
    }

    // -- wetting BC --

    #[test]
    fn test_wetting_bc_neutral() {
        let wbc = WettingBC::new(std::f64::consts::FRAC_PI_2, 1.0, 0.5);
        assert!(wbc.is_neutral());
        assert!(wbc.normal_gradient().abs() < 1e-14);
    }

    #[test]
    fn test_wetting_bc_hydrophilic() {
        let wbc = WettingBC::new(std::f64::consts::FRAC_PI_4, 1.0, 0.5);
        assert!(wbc.is_hydrophilic());
        assert!(!wbc.is_hydrophobic());
    }

    #[test]
    fn test_wetting_bc_hydrophobic() {
        let wbc = WettingBC::new(3.0 * std::f64::consts::FRAC_PI_4, 1.0, 0.5);
        assert!(wbc.is_hydrophobic());
    }

    #[test]
    fn test_wetting_bc_wall_energy() {
        let wbc = WettingBC::new(0.0, 1.0, 1.0); // theta=0, fully hydrophilic
        let e = wbc.wall_energy(1.0);
        // E = -sqrt(2*1*1) * cos(0) * 1 = -sqrt(2)
        assert!((e + std::f64::consts::SQRT_2).abs() < 1e-10);
    }

    // -- initialization --

    #[test]
    fn test_init_random_mean() {
        let mut s = make_solver(16, 16);
        s.init_random(0.0, 0.01);
        let mean = s.mean_concentration();
        assert!(mean.abs() < 0.1, "mean = {mean}");
    }

    #[test]
    fn test_init_flat_interface_profile() {
        let mut s = make_solver(16, 32);
        s.init_flat_interface(16.0, 2.0);
        // At j = 16 (interface), C ~ 0
        let c_mid = s.concentration[s.idx(8, 16)];
        assert!(c_mid.abs() < 0.1, "c_mid = {c_mid}");
    }

    #[test]
    fn test_init_droplet() {
        let mut s = make_solver(32, 32);
        s.init_droplet(16.0, 16.0, 8.0, 2.0);
        // Center should be close to +1 (inside droplet)
        let c_center = s.concentration[s.idx(16, 16)];
        assert!(c_center > 0.5, "c_center = {c_center}");
    }

    #[test]
    fn test_init_stripes() {
        let mut s = make_solver(16, 16);
        s.init_stripes(1);
        // At i=0, cos(0) = 1
        let c0 = s.concentration[s.idx(0, 0)];
        assert!((c0 - 1.0).abs() < 1e-10);
    }

    // -- Laplacian --

    #[test]
    fn test_laplacian_uniform_zero() {
        let mut s = make_solver(8, 8);
        for c in s.concentration.iter_mut() {
            *c = 0.5;
        }
        for i in 0..8 {
            for j in 0..8 {
                let lap = s.laplacian_c(i, j);
                assert!(lap.abs() < 1e-12, "lap = {lap}");
            }
        }
    }

    // -- chemical potential --

    #[test]
    fn test_chemical_potential_at_equilibrium() {
        let mut s = make_solver(8, 8);
        for c in s.concentration.iter_mut() {
            *c = 1.0; // Equilibrium
        }
        s.compute_chemical_potential();
        for &mu in s.chemical_potential.iter() {
            assert!(mu.abs() < 1e-10, "mu = {mu}");
        }
    }

    #[test]
    fn test_chemical_potential_at_zero() {
        let mut s = make_solver(8, 8);
        s.compute_chemical_potential();
        for &mu in s.chemical_potential.iter() {
            assert!(mu.abs() < 1e-12, "mu = {mu}");
        }
    }

    // -- time stepping --

    #[test]
    fn test_step_no_nan() {
        let mut s = make_solver(16, 16);
        s.init_random(0.0, 0.05);
        for _ in 0..10 {
            s.step();
        }
        for &c in s.concentration.iter() {
            assert!(!c.is_nan(), "NaN in concentration");
        }
    }

    #[test]
    fn test_conservation_during_step() {
        let mut s = make_solver(16, 16);
        s.init_random(0.0, 0.05);
        let c0 = s.total_concentration();
        for _ in 0..5 {
            s.step();
        }
        let c1 = s.total_concentration();
        assert!(
            (c1 - c0).abs() < 1.0,
            "conservation violated: {} vs {}",
            c0,
            c1
        );
    }

    #[test]
    fn test_variance_grows_spinodal() {
        let mut s = CahnHilliardLbm::new(16, 16, 1.0, 0.1, 0.1, MobilityType::Constant);
        s.init_random(0.0, 0.01);
        let v0 = s.variance();
        for _ in 0..50 {
            s.step();
        }
        let v1 = s.variance();
        assert!(v1 >= v0, "variance should grow: {} vs {}", v0, v1);
    }

    // -- free energy --

    #[test]
    fn test_free_energy_zero_at_equilibrium() {
        let mut s = make_solver(8, 8);
        for c in s.concentration.iter_mut() {
            *c = 1.0;
        }
        let fe = s.free_energy();
        assert!(fe.abs() < 1e-10, "fe = {fe}");
    }

    #[test]
    fn test_free_energy_positive_at_interface() {
        let mut s = make_solver(16, 16);
        // Sharp interface
        for i in 0..16 {
            for j in 0..16 {
                let idx = i * 16 + j;
                s.concentration[idx] = if i < 8 { 1.0 } else { -1.0 };
            }
        }
        let fe = s.free_energy();
        assert!(fe > 0.0, "fe = {fe}");
    }

    // -- analysis functions --

    #[test]
    fn test_volume_fraction_half() {
        let phi = vec![1.0, 1.0, -1.0, -1.0];
        let vf = volume_fraction_positive(&phi);
        assert!((vf - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_mean_concentration_zero() {
        let phi = vec![1.0, -1.0, 1.0, -1.0];
        assert!(mean_concentration(&phi).abs() < 1e-14);
    }

    #[test]
    fn test_concentration_variance_uniform() {
        let phi = vec![0.5; 100];
        assert!(concentration_variance(&phi).abs() < 1e-14);
    }

    // -- spinodal decomposition analysis --

    #[test]
    fn test_spinodal_region_symmetric() {
        let (lo, hi) = spinodal_region(1.0);
        assert!((lo + hi).abs() < 1e-10);
        assert!(lo < 0.0);
        assert!(hi > 0.0);
    }

    #[test]
    fn test_critical_wavenumber_at_c0_zero() {
        let kc = critical_wavenumber(0.0, 1.0, 0.5);
        let expected = (1.0 / 0.5_f64).sqrt();
        assert!((kc - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fastest_growing_wavenumber() {
        let kmax = fastest_growing_wavenumber(0.0, 1.0, 0.5);
        let kc = critical_wavenumber(0.0, 1.0, 0.5);
        assert!((kmax - kc / std::f64::consts::SQRT_2).abs() < 1e-10);
    }

    #[test]
    fn test_max_growth_rate_positive_inside_spinodal() {
        let sigma = max_growth_rate(0.0, 0.1, 1.0, 0.5);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_max_growth_rate_zero_outside_spinodal() {
        let sigma = max_growth_rate(0.9, 0.1, 1.0, 0.5);
        assert!(sigma.abs() < 1e-10, "sigma = {sigma}");
    }

    #[test]
    fn test_linear_growth_rate_zero_at_k_zero() {
        let s = linear_growth_rate(0.0, 0.0, 0.1, 1.0, 0.5);
        assert!(s.abs() < 1e-14);
    }

    // -- gradient computation --

    #[test]
    fn test_gradient_uniform_zero() {
        let mut s = make_solver(8, 8);
        for c in s.concentration.iter_mut() {
            *c = 1.0;
        }
        let g = s.gradient_c(4, 4);
        assert!(g[0].abs() < 1e-12);
        assert!(g[1].abs() < 1e-12);
    }

    // -- MRT --

    #[test]
    fn test_mrt_transform_first_row_all_ones() {
        let m = mrt_transform_matrix();
        for q in 0..9 {
            assert!((m[q] - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn test_mrt_relaxation_rates_conserved() {
        let s = mrt_relaxation_rates(1.0);
        assert!(s[0].abs() < 1e-14); // density conserved
        assert!(s[3].abs() < 1e-14); // x-momentum conserved
        assert!(s[5].abs() < 1e-14); // y-momentum conserved
    }

    // -- interface analysis --

    #[test]
    fn test_interfacial_area_uniform_zero() {
        let phi = vec![1.0; 64];
        let area = interfacial_area(&phi, 8, 8);
        assert!(area.abs() < 1e-12);
    }

    #[test]
    fn test_interfacial_area_positive_with_interface() {
        let mut phi = vec![0.0; 64];
        for i in 0..8 {
            for j in 0..8 {
                phi[i * 8 + j] = if i < 4 { 1.0 } else { -1.0 };
            }
        }
        let area = interfacial_area(&phi, 8, 8);
        assert!(area > 0.0);
    }

    #[test]
    fn test_cahn_number_positive() {
        let cn = cahn_number(0.5, 1.0, 100.0);
        assert!(cn > 0.0);
        assert!(cn < 1.0);
    }

    // -- BGK collision and streaming --

    #[test]
    fn test_collide_no_nan() {
        let mut s = make_solver(8, 8);
        s.init_random(0.0, 0.05);
        s.collide();
        for &v in s.f_dist.iter() {
            assert!(v.is_finite(), "f_dist has non-finite value");
        }
    }

    #[test]
    fn test_stream_no_nan() {
        let mut s = make_solver(8, 8);
        s.init_random(0.0, 0.05);
        s.stream();
        for &v in s.f_dist.iter() {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_update_concentration_from_dist() {
        let mut s = make_solver(4, 4);
        let c_val = 0.8;
        for i in 0..4 {
            for j in 0..4 {
                for q in 0..9 {
                    let fidx = s.fidx(i, j, q);
                    s.f_dist[fidx] = W9[q] * c_val;
                }
            }
        }
        s.update_concentration_from_dist();
        for &c in s.concentration.iter() {
            assert!((c - c_val).abs() < 1e-12);
        }
    }

    // -- coarsening timescale --

    #[test]
    fn test_coarsening_timescale_finite() {
        let t = coarsening_timescale(10.0, 0.1, 0.5, 1.0);
        assert!(t.is_finite());
        assert!(t > 0.0);
    }

    // -- structure factor --

    #[test]
    fn test_structure_factor_uniform_near_zero() {
        let phi = vec![1.0; 64];
        let sk = structure_factor_1d(&phi, 8, 8);
        for &(_k, s) in &sk {
            assert!(s.abs() < 1e-8, "S(k) = {s}");
        }
    }

    // -- wetting BC integration --

    #[test]
    fn test_set_wetting_bottom() {
        let mut s = make_solver(8, 8);
        s.set_wetting_bottom(std::f64::consts::FRAC_PI_4);
        assert!(s.wetting_bottom.is_some());
        assert!(s.wetting_bottom.unwrap().is_hydrophilic());
    }

    #[test]
    fn test_set_wetting_top() {
        let mut s = make_solver(8, 8);
        s.set_wetting_top(2.0);
        assert!(s.wetting_top.is_some());
        assert!(s.wetting_top.unwrap().is_hydrophobic());
    }

    // -- run multiple steps --

    #[test]
    fn test_run_multiple_steps() {
        let mut s = make_solver(8, 8);
        s.init_random(0.0, 0.01);
        s.run(5);
        assert_eq!(s.time_step, 5);
    }

    // -- max/min concentration --

    #[test]
    fn test_max_min_concentration() {
        let mut s = make_solver(4, 4);
        for i in 0..16 {
            s.concentration[i] = i as f64 / 15.0;
        }
        assert!((s.max_concentration() - 1.0).abs() < 1e-14);
        assert!(s.min_concentration().abs() < 1e-14);
    }

    // -- GL second derivative --

    #[test]
    fn test_gl_second_derivative_at_minima() {
        // d2f/dC2 at C=1: A*(3-1) = 2A
        let d2 = gl_bulk_energy_second_derivative(1.0, 1.0);
        assert!((d2 - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_gl_second_derivative_at_zero() {
        // d2f/dC2 at C=0: A*(-1) = -A (unstable)
        let d2 = gl_bulk_energy_second_derivative(0.0, 1.0);
        assert!((d2 + 1.0).abs() < 1e-14);
    }

    // -- dominant wavenumber --

    #[test]
    fn test_dominant_wavenumber_empty() {
        let k = dominant_wavenumber(&[]);
        assert!(k.abs() < 1e-14);
    }

    #[test]
    fn test_dominant_wavenumber_single() {
        let sk = vec![(1.0, 5.0), (2.0, 3.0), (3.0, 1.0)];
        let k = dominant_wavenumber(&sk);
        assert!((k - 1.0).abs() < 1e-14);
    }
}
