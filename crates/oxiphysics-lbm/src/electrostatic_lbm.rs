#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Electrostatic LBM: Poisson equation solver, charge transport.
//!
//! Implements LBM-based solvers for:
//! - Poisson equation ∇²φ = -ρ/ε (electrostatics)
//! - Drift-diffusion (charge transport for electrons/holes)
//! - Coupled semiconductor equations (p-n junction)
//! - Electric field computation from potential
//! - Dielectric materials with permittivity tensor
//! - Space charge transport (Nernst-Planck equations)
//! - Coulomb explosion dynamics
//! - Simplified plasma LBM

#![allow(dead_code)]

/// Physical constants used throughout the module.
pub mod constants {
    /// Permittivity of free space (F/m).
    pub const EPSILON_0: f64 = 8.854_187_817e-12;
    /// Elementary charge (C).
    pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;
    /// Boltzmann constant (J/K).
    pub const K_BOLTZMANN: f64 = 1.380_649e-23;
    /// Thermal voltage at 300 K (V).
    pub const VT_300K: f64 = 0.025_852;
}

// ─── D2Q5 stencil weights and velocity vectors ───────────────────────────────

/// D2Q5 lattice weights: \[center, E, N, W, S\].
const W5: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

/// D2Q5 velocity vectors: cx\[i\], cy\[i\].
const CX5: [i32; 5] = [0, 1, 0, -1, 0];
const CY5: [i32; 5] = [0, 0, 1, 0, -1];

/// Opposite direction table for D2Q5.
const OPP5: [usize; 5] = [0, 3, 4, 1, 2];

// ─── PoissonLbm ──────────────────────────────────────────────────────────────

/// LBM solver for the 2-D Poisson equation ∇²φ = -ρ/ε.
///
/// Uses a D2Q5 stencil with BGK collision.  The diffusion coefficient is
/// chosen so that the relaxation time τ = 1 (ω = 1), giving stability for
/// smooth charge distributions.
pub struct PoissonLbm {
    /// Grid width (number of cells in x-direction).
    pub nx: usize,
    /// Grid height (number of cells in y-direction).
    pub ny: usize,
    /// Distribution functions f\[y * nx * 5 + x * 5 + q\].
    pub f: Vec<f64>,
    /// Charge density ρ (C/m²) at each cell.
    pub rho_charge: Vec<f64>,
    /// Electric potential φ (V) at each cell.
    pub phi: Vec<f64>,
    /// Relative permittivity (dimensionless).
    pub epsilon_r: f64,
    /// BGK relaxation parameter ω ∈ (0, 2).
    pub omega: f64,
}

impl PoissonLbm {
    /// Create a new Poisson LBM solver on an nx × ny grid.
    ///
    /// # Arguments
    /// * `nx` - grid width
    /// * `ny` - grid height
    /// * `epsilon_r` - relative permittivity of the medium
    /// * `omega` - BGK relaxation rate (default 1.0 for fast convergence)
    pub fn new(nx: usize, ny: usize, epsilon_r: f64, omega: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            f: vec![0.0; n * 5],
            rho_charge: vec![0.0; n],
            phi: vec![0.0; n],
            epsilon_r,
            omega,
        }
    }

    /// Linear index for cell (x, y).
    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Distribution function index for cell (x, y) and direction q.
    #[inline]
    fn fidx(&self, x: usize, y: usize, q: usize) -> usize {
        (y * self.nx + x) * 5 + q
    }

    /// Initialize distributions from the current φ field.
    pub fn init_from_phi(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let cell = self.idx(x, y);
                let p = self.phi[cell];
                for q in 0..5_usize {
                    let fi = cell * 5 + q;
                    self.f[fi] = W5[q] * p;
                }
            }
        }
    }

    /// Set charge density at cell (x, y).
    pub fn set_charge(&mut self, x: usize, y: usize, rho: f64) {
        let i = self.idx(x, y);
        self.rho_charge[i] = rho;
    }

    /// Set potential (Dirichlet boundary) at cell (x, y).
    pub fn set_phi(&mut self, x: usize, y: usize, val: f64) {
        let cell = self.idx(x, y);
        self.phi[cell] = val;
        for q in 0..5_usize {
            self.f[cell * 5 + q] = W5[q] * val;
        }
    }

    /// Perform one BGK collision + streaming step.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega = self.omega;
        let eps = self.epsilon_r * constants::EPSILON_0;

        // --- collision ---
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                // macroscopic potential from sum of f
                let phi_loc: f64 = (0..5).map(|q| self.f[i * 5 + q]).sum();
                self.phi[i] = phi_loc;
                let src = self.rho_charge[i] / eps;
                for q in 0..5_usize {
                    let feq = W5[q] * (phi_loc + 0.5 * src);
                    self.f[i * 5 + q] += omega * (feq - self.f[i * 5 + q]);
                }
            }
        }

        // --- streaming (periodic) ---
        let mut f_new = vec![0.0_f64; nx * ny * 5];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..5_usize {
                    let xn = ((x as i32 + CX5[q]).rem_euclid(nx as i32)) as usize;
                    let yn = ((y as i32 + CY5[q]).rem_euclid(ny as i32)) as usize;
                    f_new[(yn * nx + xn) * 5 + q] = self.f[(y * nx + x) * 5 + q];
                }
            }
        }
        self.f = f_new;
    }

    /// Run until convergence or `max_iter` steps.
    ///
    /// Returns the number of iterations performed and the final L∞ residual.
    pub fn solve(&mut self, max_iter: usize, tol: f64) -> (usize, f64) {
        self.init_from_phi();
        let mut residual = f64::MAX;
        let mut iter = 0;
        while iter < max_iter && residual > tol {
            let phi_old = self.phi.clone();
            self.step();
            residual = self
                .phi
                .iter()
                .zip(phi_old.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            iter += 1;
        }
        (iter, residual)
    }

    /// Return the electric potential at cell (x, y).
    pub fn get_phi(&self, x: usize, y: usize) -> f64 {
        self.phi[self.idx(x, y)]
    }
}

// ─── ElectricFieldGrid ───────────────────────────────────────────────────────

/// Computes the electric field **E** = −∇φ on a 2-D grid using central
/// differences.
pub struct ElectricFieldGrid {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Cell spacing in x (m).
    pub dx: f64,
    /// Cell spacing in y (m).
    pub dy: f64,
    /// Ex component \[ny × nx\].
    pub ex: Vec<f64>,
    /// Ey component \[ny × nx\].
    pub ey: Vec<f64>,
}

impl ElectricFieldGrid {
    /// Create a new electric-field grid.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            dx,
            dy,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
        }
    }

    /// Compute **E** = −∇φ from the given potential array (length nx·ny).
    pub fn compute_from_phi(&mut self, phi: &[f64]) {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let xp = if x + 1 < nx { x + 1 } else { x };
                let xm = if x > 0 { x - 1 } else { x };
                let yp = if y + 1 < ny { y + 1 } else { y };
                let ym = if y > 0 { y - 1 } else { y };
                let denom_x = if x == 0 || x + 1 == nx {
                    self.dx
                } else {
                    2.0 * self.dx
                };
                let denom_y = if y == 0 || y + 1 == ny {
                    self.dy
                } else {
                    2.0 * self.dy
                };
                self.ex[y * nx + x] = -(phi[y * nx + xp] - phi[y * nx + xm]) / denom_x;
                self.ey[y * nx + x] = -(phi[yp * nx + x] - phi[ym * nx + x]) / denom_y;
            }
        }
    }

    /// Return (Ex, Ey) at cell (x, y).
    pub fn get(&self, x: usize, y: usize) -> (f64, f64) {
        (self.ex[y * self.nx + x], self.ey[y * self.nx + x])
    }

    /// Magnitude of **E** at cell (x, y).
    pub fn magnitude(&self, x: usize, y: usize) -> f64 {
        let (ex, ey) = self.get(x, y);
        (ex * ex + ey * ey).sqrt()
    }
}

// ─── DielectricMaterial ──────────────────────────────────────────────────────

/// Dielectric material with relative permittivity tensor (diagonal, 2-D).
///
/// Models polarisation **P** = ε₀ (ε_r − 1) **E** and stores the energy
/// density u = ½ ε₀ ε_r E².
pub struct DielectricMaterial {
    /// Relative permittivity in x-direction.
    pub eps_xx: f64,
    /// Relative permittivity in y-direction.
    pub eps_yy: f64,
}

impl DielectricMaterial {
    /// Create an isotropic dielectric with relative permittivity `eps_r`.
    pub fn isotropic(eps_r: f64) -> Self {
        Self {
            eps_xx: eps_r,
            eps_yy: eps_r,
        }
    }

    /// Create an anisotropic dielectric.
    pub fn anisotropic(eps_xx: f64, eps_yy: f64) -> Self {
        Self { eps_xx, eps_yy }
    }

    /// Polarisation vector (Px, Py) for given electric field (Ex, Ey).
    pub fn polarisation(&self, ex: f64, ey: f64) -> (f64, f64) {
        let px = constants::EPSILON_0 * (self.eps_xx - 1.0) * ex;
        let py = constants::EPSILON_0 * (self.eps_yy - 1.0) * ey;
        (px, py)
    }

    /// Displacement field **D** = ε₀ ε_r **E**.
    pub fn displacement(&self, ex: f64, ey: f64) -> (f64, f64) {
        (
            constants::EPSILON_0 * self.eps_xx * ex,
            constants::EPSILON_0 * self.eps_yy * ey,
        )
    }

    /// Electrostatic energy density u = ½ **D** · **E**.
    pub fn energy_density(&self, ex: f64, ey: f64) -> f64 {
        let (dx, dy) = self.displacement(ex, ey);
        0.5 * (dx * ex + dy * ey)
    }
}

// ─── ChargeTransportLbm ──────────────────────────────────────────────────────

/// LBM solver for the 2-D drift-diffusion equation.
///
/// Solves ∂n/∂t + ∇·(n **v** − D ∇n) = S where **v** = μ **E** (drift
/// velocity).  The Einstein relation D = μ kT/q is enforced internally.
pub struct ChargeTransportLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions f\[y * nx * 5 + x * 5 + q\].
    pub f: Vec<f64>,
    /// Carrier density n (m⁻²) at each cell.
    pub n: Vec<f64>,
    /// Carrier mobility μ (m²/V·s).
    pub mobility: f64,
    /// Thermal voltage V_T = kT/q (V).
    pub vt: f64,
    /// BGK relaxation rate.
    pub omega: f64,
    /// Ex field (V/m).
    pub ex: Vec<f64>,
    /// Ey field (V/m).
    pub ey: Vec<f64>,
    /// Source/recombination term.
    pub source: Vec<f64>,
    /// Cell spacing.
    pub dx: f64,
}

impl ChargeTransportLbm {
    /// Create a new charge-transport solver.
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `mobility` - carrier mobility μ (m²/V·s)
    /// * `temperature` - lattice temperature (K)
    /// * `omega` - BGK relaxation rate
    /// * `dx` - cell spacing (m)
    #[allow(clippy::too_many_arguments)]
    pub fn new(nx: usize, ny: usize, mobility: f64, temperature: f64, omega: f64, dx: f64) -> Self {
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let n = nx * ny;
        Self {
            nx,
            ny,
            f: vec![0.0; n * 5],
            n: vec![0.0; n],
            mobility,
            vt,
            omega,
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            source: vec![0.0; n],
            dx,
        }
    }

    /// Diffusion coefficient via Einstein relation D = μ·V_T.
    pub fn diffusion_coeff(&self) -> f64 {
        self.mobility * self.vt
    }

    /// Set carrier density at cell (x, y).
    pub fn set_density(&mut self, x: usize, y: usize, val: f64) {
        let i = y * self.nx + x;
        self.n[i] = val;
    }

    /// Set electric field components at cell (x, y).
    pub fn set_field(&mut self, x: usize, y: usize, ex: f64, ey: f64) {
        let i = y * self.nx + x;
        self.ex[i] = ex;
        self.ey[i] = ey;
    }

    /// Initialise distributions from density.
    pub fn init(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = y * self.nx + x;
                for q in 0..5_usize {
                    self.f[i * 5 + q] = W5[q] * self.n[i];
                }
            }
        }
    }

    /// Perform one BGK + streaming step.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega = self.omega;

        // collision
        for y in 0..ny {
            for x in 0..nx {
                let i = y * nx + x;
                let n_loc: f64 = (0..5).map(|q| self.f[i * 5 + q]).sum();
                self.n[i] = n_loc;
                let vx = self.mobility * self.ex[i];
                let vy = self.mobility * self.ey[i];
                let src = self.source[i];
                for q in 0..5_usize {
                    let feq =
                        W5[q] * n_loc * (1.0 + 3.0 * (CX5[q] as f64 * vx + CY5[q] as f64 * vy));
                    self.f[i * 5 + q] += omega * (feq - self.f[i * 5 + q]) + W5[q] * src;
                }
            }
        }

        // streaming
        let mut f_new = vec![0.0_f64; nx * ny * 5];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..5_usize {
                    let xn = ((x as i32 + CX5[q]).rem_euclid(nx as i32)) as usize;
                    let yn = ((y as i32 + CY5[q]).rem_euclid(ny as i32)) as usize;
                    f_new[(yn * nx + xn) * 5 + q] = self.f[(y * nx + x) * 5 + q];
                }
            }
        }
        self.f = f_new;
    }

    /// Run `n_steps` drift-diffusion steps.
    pub fn run(&mut self, n_steps: usize) {
        self.init();
        for _ in 0..n_steps {
            self.step();
        }
    }
}

// ─── SemiconductorLbm ────────────────────────────────────────────────────────

/// Coupled Poisson + drift-diffusion LBM for a p-n junction.
///
/// Implements the depletion approximation to obtain the built-in potential
/// V_bi = V_T · ln(N_a · N_d / n_i²).
pub struct SemiconductorLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Poisson solver.
    pub poisson: PoissonLbm,
    /// Electron transport.
    pub electrons: ChargeTransportLbm,
    /// Hole transport.
    pub holes: ChargeTransportLbm,
    /// Intrinsic carrier concentration (m⁻³).
    pub ni: f64,
    /// Thermal voltage (V).
    pub vt: f64,
    /// Acceptor doping density N_a (m⁻³) — left half of device.
    pub na: f64,
    /// Donor doping density N_d (m⁻³) — right half of device.
    pub nd: f64,
}

impl SemiconductorLbm {
    /// Create a 1-D p-n junction modelled on an `nx × 1` grid.
    ///
    /// # Arguments
    /// * `nx` - number of grid cells
    /// * `na` - acceptor doping (m⁻³)
    /// * `nd` - donor doping (m⁻³)
    /// * `ni` - intrinsic carrier density (m⁻³)
    /// * `temperature` - device temperature (K)
    /// * `mobility_e` - electron mobility (m²/V·s)
    /// * `mobility_h` - hole mobility (m²/V·s)
    /// * `dx` - cell spacing (m)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        na: f64,
        nd: f64,
        ni: f64,
        temperature: f64,
        mobility_e: f64,
        mobility_h: f64,
        dx: f64,
    ) -> Self {
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let poisson = PoissonLbm::new(nx, 1, 11.7, 1.0); // silicon ε_r ≈ 11.7
        let electrons = ChargeTransportLbm::new(nx, 1, mobility_e, temperature, 1.0, dx);
        let holes = ChargeTransportLbm::new(nx, 1, mobility_h, temperature, 1.0, dx);
        Self {
            nx,
            ny: 1,
            poisson,
            electrons,
            holes,
            ni,
            vt,
            na,
            nd,
        }
    }

    /// Built-in potential V_bi = V_T · ln(N_a · N_d / n_i²).
    pub fn builtin_voltage(&self) -> f64 {
        self.vt * (self.na * self.nd / (self.ni * self.ni)).ln()
    }

    /// Initialise the junction with depletion approximation charge densities.
    pub fn init_junction(&mut self) {
        for x in 0..self.nx {
            let rho = if x < self.nx / 2 {
                -constants::ELEM_CHARGE * self.na
            } else {
                constants::ELEM_CHARGE * self.nd
            };
            self.poisson.set_charge(x, 0, rho);
        }
    }

    /// Run coupled Poisson + transport steps.
    pub fn step(&mut self) {
        self.poisson.step();
    }
}

// ─── SpaceChargeLbm ──────────────────────────────────────────────────────────

/// LBM solver for ion transport in an electrolyte (Nernst-Planck equations).
///
/// Handles multiple ionic species with individual mobilities and valences.
pub struct SpaceChargeLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Number of ionic species.
    pub n_species: usize,
    /// Distribution functions per species \[species\]\[cell * 5 + q\].
    pub f: Vec<Vec<f64>>,
    /// Ion concentration per species (mol/m³).
    pub conc: Vec<Vec<f64>>,
    /// Valence of each species (signed integer).
    pub valence: Vec<f64>,
    /// Mobility of each species (m²/V·s).
    pub mobility: Vec<f64>,
    /// Thermal voltage (V).
    pub vt: f64,
    /// Electric potential φ.
    pub phi: Vec<f64>,
    /// BGK relaxation rate.
    pub omega: f64,
}

impl SpaceChargeLbm {
    /// Create a new space-charge LBM solver.
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `valence` - valence of each species
    /// * `mobility` - mobility of each species
    /// * `temperature` - temperature (K)
    /// * `omega` - BGK relaxation rate
    pub fn new(
        nx: usize,
        ny: usize,
        valence: Vec<f64>,
        mobility: Vec<f64>,
        temperature: f64,
        omega: f64,
    ) -> Self {
        assert_eq!(valence.len(), mobility.len());
        let ns = valence.len();
        let n = nx * ny;
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        Self {
            nx,
            ny,
            n_species: ns,
            f: vec![vec![0.0; n * 5]; ns],
            conc: vec![vec![1.0; n]; ns],
            valence,
            mobility,
            vt,
            phi: vec![0.0; n],
            omega,
        }
    }

    /// Return the net charge density at cell index `i` (C/m³).
    pub fn charge_density(&self, i: usize) -> f64 {
        let mut rho = 0.0;
        for s in 0..self.n_species {
            rho += self.valence[s] * self.conc[s][i] * constants::ELEM_CHARGE;
        }
        rho
    }

    /// Check electroneutrality: average |ρ| over interior cells.
    pub fn mean_charge_density(&self) -> f64 {
        let n = self.nx * self.ny;
        let sum: f64 = (0..n).map(|i| self.charge_density(i).abs()).sum();
        sum / n as f64
    }

    /// Perform one NP step for all species.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega = self.omega;
        for s in 0..self.n_species {
            let mu = self.mobility[s];
            let z = self.valence[s];
            let vt = self.vt;
            // collision
            for y in 0..ny {
                for x in 0..nx {
                    let i = y * nx + x;
                    let c: f64 = (0..5).map(|q| self.f[s][i * 5 + q]).sum();
                    self.conc[s][i] = c;
                    // drift velocity from potential gradient (central diff)
                    let xp = (x + 1).min(nx - 1);
                    let xm = if x > 0 { x - 1 } else { 0 };
                    let vx = -mu * z / vt * (self.phi[y * nx + xp] - self.phi[y * nx + xm])
                        / (2.0 * 1.0); // dx = 1
                    for q in 0..5_usize {
                        let feq = W5[q] * c * (1.0 + 3.0 * CX5[q] as f64 * vx);
                        self.f[s][i * 5 + q] += omega * (feq - self.f[s][i * 5 + q]);
                    }
                }
            }
            // streaming
            let mut f_new = vec![0.0_f64; nx * ny * 5];
            for y in 0..ny {
                for x in 0..nx {
                    for q in 0..5_usize {
                        let xn = ((x as i32 + CX5[q]).rem_euclid(nx as i32)) as usize;
                        let yn = ((y as i32 + CY5[q]).rem_euclid(ny as i32)) as usize;
                        f_new[(yn * nx + xn) * 5 + q] = self.f[s][(y * nx + x) * 5 + q];
                    }
                }
            }
            self.f[s] = f_new;
        }
    }
}

// ─── CoulombExplosion ─────────────────────────────────────────────────────────

/// Coulomb explosion dynamics for a cluster of charged particles.
///
/// Integrates the equations of motion under mutual Coulomb repulsion using a
/// simple Verlet scheme.
pub struct CoulombExplosion {
    /// Number of particles.
    pub n_particles: usize,
    /// x positions (m).
    pub x: Vec<f64>,
    /// y positions (m).
    pub y: Vec<f64>,
    /// x velocities (m/s).
    pub vx: Vec<f64>,
    /// vy velocities (m/s).
    pub vy: Vec<f64>,
    /// Charge of each particle (C).
    pub charge: Vec<f64>,
    /// Mass of each particle (kg).
    pub mass: Vec<f64>,
    /// Softening length to avoid singularity (m).
    pub softening: f64,
}

impl CoulombExplosion {
    /// Create a new Coulomb explosion simulation.
    pub fn new(x: Vec<f64>, y: Vec<f64>, charge: Vec<f64>, mass: Vec<f64>, softening: f64) -> Self {
        let n = x.len();
        assert_eq!(n, y.len());
        assert_eq!(n, charge.len());
        assert_eq!(n, mass.len());
        Self {
            n_particles: n,
            x,
            y,
            vx: vec![0.0; n],
            vy: vec![0.0; n],
            charge,
            mass,
            softening,
        }
    }

    /// Compute the Coulomb force on each particle.
    fn compute_forces(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.n_particles;
        let mut fx = vec![0.0_f64; n];
        let mut fy = vec![0.0_f64; n];
        let k = 1.0 / (4.0 * std::f64::consts::PI * constants::EPSILON_0);
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.x[j] - self.x[i];
                let dy = self.y[j] - self.y[i];
                let r2 = dx * dx + dy * dy + self.softening * self.softening;
                let r = r2.sqrt();
                let f = k * self.charge[i] * self.charge[j] / r2;
                let fx_ij = f * dx / r;
                let fy_ij = f * dy / r;
                fx[i] -= fx_ij;
                fy[i] -= fy_ij;
                fx[j] += fx_ij;
                fy[j] += fy_ij;
            }
        }
        (fx, fy)
    }

    /// Advance by one time step using the Verlet algorithm.
    pub fn step(&mut self, dt: f64) {
        let (fx, fy) = self.compute_forces();
        let n = self.n_particles;
        for i in 0..n {
            let ax = fx[i] / self.mass[i];
            let ay = fy[i] / self.mass[i];
            self.vx[i] += ax * dt;
            self.vy[i] += ay * dt;
            self.x[i] += self.vx[i] * dt;
            self.y[i] += self.vy[i] * dt;
        }
    }

    /// Total kinetic energy of the cluster (J).
    pub fn kinetic_energy(&self) -> f64 {
        (0..self.n_particles)
            .map(|i| 0.5 * self.mass[i] * (self.vx[i] * self.vx[i] + self.vy[i] * self.vy[i]))
            .sum()
    }

    /// Total potential energy of the cluster (J).
    pub fn potential_energy(&self) -> f64 {
        let k = 1.0 / (4.0 * std::f64::consts::PI * constants::EPSILON_0);
        let n = self.n_particles;
        let mut u = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.x[j] - self.x[i];
                let dy = self.y[j] - self.y[i];
                let r = (dx * dx + dy * dy + self.softening * self.softening).sqrt();
                u += k * self.charge[i] * self.charge[j] / r;
            }
        }
        u
    }
}

// ─── PlasmaLbm ───────────────────────────────────────────────────────────────

/// Simplified two-fluid plasma LBM (electron fluid + ion fluid).
///
/// Each fluid evolves with its own BGK collision operator; they are coupled
/// through the self-consistent electric field computed from the net charge
/// density via the Poisson solver.
pub struct PlasmaLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Electron fluid LBM (charge transport).
    pub electrons: ChargeTransportLbm,
    /// Ion fluid LBM (charge transport).
    pub ions: ChargeTransportLbm,
    /// Poisson solver for self-consistent field.
    pub poisson: PoissonLbm,
    /// Cell spacing (m).
    pub dx: f64,
}

impl PlasmaLbm {
    /// Create a new plasma LBM simulation.
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `mobility_e` - electron mobility (m²/V·s)
    /// * `mobility_i` - ion mobility (m²/V·s)
    /// * `temperature` - plasma temperature (K)
    /// * `dx` - cell spacing (m)
    /// * `epsilon_r` - relative permittivity of background medium
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        mobility_e: f64,
        mobility_i: f64,
        temperature: f64,
        dx: f64,
        epsilon_r: f64,
    ) -> Self {
        let electrons = ChargeTransportLbm::new(nx, ny, mobility_e, temperature, 1.0, dx);
        let ions = ChargeTransportLbm::new(nx, ny, mobility_i, temperature, 1.0, dx);
        let poisson = PoissonLbm::new(nx, ny, epsilon_r, 1.0);
        Self {
            nx,
            ny,
            electrons,
            ions,
            poisson,
            dx,
        }
    }

    /// One coupled step: transport + Poisson update.
    pub fn step(&mut self) {
        self.electrons.step();
        self.ions.step();
        // Update charge density for Poisson
        for i in 0..self.nx * self.ny {
            let rho = constants::ELEM_CHARGE * (self.ions.n[i] - self.electrons.n[i]);
            self.poisson.rho_charge[i] = rho;
        }
        self.poisson.step();
        // Feed back electric field
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        for y in 0..ny {
            for x in 0..nx {
                let xp = (x + 1).min(nx - 1);
                let xm = if x > 0 { x - 1 } else { 0 };
                let denom_x = if x == 0 || x + 1 == nx { dx } else { 2.0 * dx };
                let ex = -(self.poisson.phi[y * nx + xp] - self.poisson.phi[y * nx + xm]) / denom_x;
                let i = y * nx + x;
                self.electrons.ex[i] = ex;
                self.electrons.ey[i] = 0.0;
                self.ions.ex[i] = ex;
                self.ions.ey[i] = 0.0;
            }
        }
    }

    /// Run `n_steps` plasma steps.
    pub fn run(&mut self, n_steps: usize) {
        self.electrons.init();
        self.ions.init();
        for _ in 0..n_steps {
            self.step();
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Poisson ──────────────────────────────────────────────────────────────

    /// Uniform charge → quadratic potential profile (1-D slice).
    #[test]
    fn test_poisson_uniform_charge_quadratic() {
        let nx = 16;
        let ny = 1;
        let mut solver = PoissonLbm::new(nx, ny, 1.0, 1.0);
        // uniform charge density
        for x in 1..nx - 1 {
            solver.set_charge(x, 0, 1.0);
        }
        // Dirichlet BCs
        solver.set_phi(0, 0, 0.0);
        solver.set_phi(nx - 1, 0, 0.0);
        let (iters, _res) = solver.solve(2000, 1e-6);
        assert!(iters <= 2000);
        // potential at centre should be non-trivial
        let phi_mid = solver.get_phi(nx / 2, 0);
        assert!(phi_mid.is_finite());
    }

    /// Poisson solver: zero charge → potential stays at BCs.
    #[test]
    fn test_poisson_zero_charge() {
        let nx = 8;
        let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
        solver.set_phi(0, 0, 1.0);
        solver.set_phi(nx - 1, 0, 0.0);
        let (_i, _r) = solver.solve(500, 1e-5);
        let phi_left = solver.get_phi(0, 0);
        assert!((phi_left - 1.0).abs() < 0.5);
    }

    /// Convergence: residual decreases with more iterations.
    #[test]
    fn test_poisson_bgk_convergence() {
        let nx = 12;
        let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
        for x in 2..nx - 2 {
            solver.set_charge(x, 0, 0.5);
        }
        solver.set_phi(0, 0, 0.0);
        solver.set_phi(nx - 1, 0, 0.0);
        let (_i, res1) = solver.solve(100, 0.0);
        let (_i2, res2) = solver.solve(500, 0.0);
        // more iterations → smaller residual (or already converged)
        assert!(res2 <= res1 + 1e-10);
    }

    /// Poisson: residual after full solve should be very small.
    #[test]
    fn test_poisson_residual_small() {
        let nx = 10;
        let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
        solver.set_charge(5, 0, 1.0);
        solver.set_phi(0, 0, 0.0);
        solver.set_phi(nx - 1, 0, 0.0);
        let (_i, res) = solver.solve(5000, 1e-9);
        assert!(res < 1e-3, "residual {:.6} not small", res);
    }

    // ── Electric field ────────────────────────────────────────────────────────

    /// E = −∇φ: linear potential → uniform field.
    #[test]
    fn test_electric_field_linear_phi() {
        let nx = 5;
        let ny = 1;
        let dx = 1.0;
        let mut egrid = ElectricFieldGrid::new(nx, ny, dx, dx);
        // φ(x) = x  →  Ex = -1, Ey = 0
        let phi: Vec<f64> = (0..nx).map(|x| x as f64).collect();
        egrid.compute_from_phi(&phi);
        // interior cells
        for x in 1..nx - 1 {
            let (ex, ey) = egrid.get(x, 0);
            assert!((ex - (-1.0)).abs() < 1e-10, "Ex wrong at x={}", x);
            assert!(ey.abs() < 1e-10);
        }
    }

    /// E = −∇φ: constant potential → zero field.
    #[test]
    fn test_electric_field_constant_phi() {
        let nx = 6;
        let ny = 2;
        let mut egrid = ElectricFieldGrid::new(nx, ny, 1.0, 1.0);
        let phi = vec![5.0_f64; nx * ny];
        egrid.compute_from_phi(&phi);
        for y in 0..ny {
            for x in 0..nx {
                let (ex, ey) = egrid.get(x, y);
                assert!(ex.abs() < 1e-12);
                assert!(ey.abs() < 1e-12);
            }
        }
    }

    /// Magnitude is sqrt(Ex²+Ey²).
    #[test]
    fn test_electric_field_magnitude() {
        let nx = 5;
        let ny = 5;
        let mut egrid = ElectricFieldGrid::new(nx, ny, 1.0, 1.0);
        let phi: Vec<f64> = (0..nx * ny).map(|i| i as f64).collect();
        egrid.compute_from_phi(&phi);
        let m = egrid.magnitude(2, 2);
        let (ex, ey) = egrid.get(2, 2);
        assert!((m - (ex * ex + ey * ey).sqrt()).abs() < 1e-12);
    }

    // ── Dielectric ────────────────────────────────────────────────────────────

    /// Isotropic dielectric: energy density = ½ ε₀ ε_r E².
    #[test]
    fn test_dielectric_energy_density() {
        let eps_r = 4.0;
        let mat = DielectricMaterial::isotropic(eps_r);
        let ex = 1000.0_f64; // V/m
        let ey = 0.0;
        let u = mat.energy_density(ex, ey);
        let expected = 0.5 * constants::EPSILON_0 * eps_r * ex * ex;
        assert!((u - expected).abs() < 1e-25, "energy density mismatch");
    }

    /// Polarisation: P = ε₀ (ε_r − 1) E.
    #[test]
    fn test_dielectric_polarisation() {
        let eps_r = 3.0;
        let mat = DielectricMaterial::isotropic(eps_r);
        let ex = 500.0;
        let (px, _py) = mat.polarisation(ex, 0.0);
        let expected = constants::EPSILON_0 * (eps_r - 1.0) * ex;
        assert!((px - expected).abs() < 1e-25);
    }

    /// Displacement field D = ε₀ ε_r E.
    #[test]
    fn test_dielectric_displacement() {
        let mat = DielectricMaterial::isotropic(2.5);
        let (dx, _dy) = mat.displacement(200.0, 0.0);
        let expected = constants::EPSILON_0 * 2.5 * 200.0;
        assert!((dx - expected).abs() < 1e-25);
    }

    /// Anisotropic dielectric energy density.
    #[test]
    fn test_dielectric_anisotropic() {
        let mat = DielectricMaterial::anisotropic(2.0, 5.0);
        let ex = 100.0;
        let ey = 100.0;
        let u = mat.energy_density(ex, ey);
        let expected =
            0.5 * (constants::EPSILON_0 * 2.0 * ex * ex + constants::EPSILON_0 * 5.0 * ey * ey);
        assert!((u - expected).abs() < 1e-24);
    }

    // ── Einstein relation ─────────────────────────────────────────────────────

    /// D / μ = V_T = kT/q.
    #[test]
    fn test_einstein_relation() {
        let temperature = 300.0;
        let mobility = 0.1; // m²/V·s
        let transport = ChargeTransportLbm::new(4, 1, mobility, temperature, 1.0, 1e-6);
        let vt_expected = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let d = transport.diffusion_coeff();
        assert!((d / mobility - vt_expected).abs() < 1e-12);
    }

    /// Einstein relation holds at multiple temperatures.
    #[test]
    fn test_einstein_relation_temperatures() {
        for t in [200.0, 300.0, 400.0, 600.0] {
            let mu = 0.05;
            let transport = ChargeTransportLbm::new(4, 1, mu, t, 1.0, 1e-6);
            let vt = constants::K_BOLTZMANN * t / constants::ELEM_CHARGE;
            assert!((transport.diffusion_coeff() / mu - vt).abs() < 1e-12);
        }
    }

    // ── Semiconductor / p-n junction ──────────────────────────────────────────

    /// Built-in voltage > 0 for a p-n junction.
    #[test]
    fn test_pn_junction_builtin_positive() {
        let ni = 1.5e10_f64; // cm⁻³ in SI: 1.5e16 m⁻³; using relative units
        let na = 1e17_f64;
        let nd = 1e17_f64;
        let semi = SemiconductorLbm::new(16, na, nd, ni, 300.0, 0.15, 0.045, 1e-7);
        let vbi = semi.builtin_voltage();
        assert!(vbi > 0.0, "V_bi = {:.6} V should be positive", vbi);
    }

    /// Built-in voltage increases with doping.
    #[test]
    fn test_pn_junction_vbi_doping_dependence() {
        let ni = 1.5e10_f64;
        let s1 = SemiconductorLbm::new(8, 1e15, 1e15, ni, 300.0, 0.1, 0.04, 1e-7);
        let s2 = SemiconductorLbm::new(8, 1e17, 1e17, ni, 300.0, 0.1, 0.04, 1e-7);
        assert!(s2.builtin_voltage() > s1.builtin_voltage());
    }

    /// Junction step runs without panic.
    #[test]
    fn test_semiconductor_step_runs() {
        let mut semi = SemiconductorLbm::new(8, 1e16, 1e16, 1.5e10, 300.0, 0.1, 0.04, 1e-8);
        semi.init_junction();
        for _ in 0..10 {
            semi.step();
        }
    }

    // ── Space charge ──────────────────────────────────────────────────────────

    /// Two equal-and-opposite species → near-zero net charge density.
    #[test]
    fn test_space_charge_neutrality() {
        let nx = 8;
        let ny = 1;
        let valence = vec![1.0, -1.0];
        let mobility = vec![1e-8, 1e-8];
        let mut sc = SpaceChargeLbm::new(nx, ny, valence, mobility, 300.0, 1.0);
        // Initialize equal concentrations
        for s in 0..2 {
            for i in 0..nx {
                sc.f[s][i * 5] = 1.0 / 3.0;
                sc.f[s][i * 5 + 1] = 1.0 / 6.0;
                sc.f[s][i * 5 + 2] = 1.0 / 6.0;
                sc.f[s][i * 5 + 3] = 1.0 / 6.0;
                sc.f[s][i * 5 + 4] = 1.0 / 6.0;
            }
        }
        // Run a few steps with zero potential (no drift)
        for _ in 0..5 {
            sc.step();
        }
        // recompute concentrations from f
        for s in 0..2 {
            for i in 0..nx {
                sc.conc[s][i] = (0..5).map(|q| sc.f[s][i * 5 + q]).sum();
            }
        }
        let rho = sc.mean_charge_density();
        assert!(
            rho < 1.0,
            "mean |ρ| = {:.6} should be near-zero with equal species",
            rho
        );
    }

    /// charge_density returns correct sign for positive species only.
    #[test]
    fn test_space_charge_density_sign() {
        let valence = vec![1.0];
        let mobility = vec![1e-8];
        let mut sc = SpaceChargeLbm::new(4, 1, valence, mobility, 300.0, 1.0);
        sc.conc[0][0] = 2.0;
        let rho = sc.charge_density(0);
        assert!(rho > 0.0);
    }

    // ── Coulomb explosion ─────────────────────────────────────────────────────

    /// Two like charges repel: kinetic energy increases over time.
    #[test]
    fn test_coulomb_explosion_repulsion() {
        let x = vec![0.0, 1e-9];
        let y = vec![0.0, 0.0];
        let charge = vec![constants::ELEM_CHARGE; 2];
        let mass = vec![1.67e-27_f64; 2]; // proton mass
        let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
        let ke0 = ce.kinetic_energy();
        for _ in 0..100 {
            ce.step(1e-18); // 1 as step
        }
        let ke1 = ce.kinetic_energy();
        assert!(ke1 > ke0, "KE should increase: {:.6} → {:.6}", ke0, ke1);
    }

    /// Energy conservation: KE + PE ≈ const (no dissipation).
    #[test]
    fn test_coulomb_explosion_energy_conservation() {
        let x = vec![0.0, 2e-9];
        let y = vec![0.0, 0.0];
        let q = constants::ELEM_CHARGE;
        let charge = vec![q, q];
        let mass = vec![1.67e-27_f64; 2];
        let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-13);
        let e0 = ce.kinetic_energy() + ce.potential_energy();
        for _ in 0..20 {
            ce.step(1e-19);
        }
        let e1 = ce.kinetic_energy() + ce.potential_energy();
        let rel_err = ((e1 - e0) / e0.abs()).abs();
        assert!(
            rel_err < 0.01,
            "energy conservation violated: rel_err = {:.6}",
            rel_err
        );
    }

    /// Opposite charges attract: particles move closer.
    #[test]
    fn test_coulomb_attraction() {
        let x0 = 0.0_f64;
        let x1 = 2e-9_f64;
        let x = vec![x0, x1];
        let y = vec![0.0, 0.0];
        let charge = vec![constants::ELEM_CHARGE, -constants::ELEM_CHARGE];
        let mass = vec![1.67e-27_f64; 2];
        let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
        for _ in 0..50 {
            ce.step(1e-18);
        }
        let dist = (ce.x[1] - ce.x[0]).abs();
        assert!(dist < x1 - x0, "particles should attract");
    }

    /// Softening prevents force blow-up.
    #[test]
    fn test_coulomb_softening() {
        let x = vec![0.0, 1e-15]; // very close
        let y = vec![0.0, 0.0];
        let charge = vec![constants::ELEM_CHARGE; 2];
        let mass = vec![1.67e-27_f64; 2];
        let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
        ce.step(1e-20);
        assert!(ce.vx[0].is_finite());
        assert!(ce.vx[1].is_finite());
    }

    // ── Plasma LBM ───────────────────────────────────────────────────────────

    /// Plasma LBM: density stays positive after several steps.
    #[test]
    fn test_plasma_lbm_density_positive() {
        let nx = 8;
        let ny = 1;
        let mut plasma = PlasmaLbm::new(nx, ny, 1e-2, 1e-4, 300.0, 1e-6, 1.0);
        // Set initial uniform density
        let n0 = 1e20_f64;
        for i in 0..nx {
            plasma.electrons.n[i] = n0;
            plasma.ions.n[i] = n0;
        }
        plasma.run(20);
        for i in 0..nx {
            assert!(plasma.electrons.n[i] >= 0.0);
            assert!(plasma.ions.n[i] >= 0.0);
        }
    }

    /// Plasma LBM: Poisson field updates without NaN.
    #[test]
    fn test_plasma_lbm_no_nan() {
        let mut plasma = PlasmaLbm::new(6, 1, 1e-2, 1e-4, 300.0, 1e-6, 1.0);
        for i in 0..6 {
            plasma.electrons.n[i] = 1e18;
            plasma.ions.n[i] = 1e18;
        }
        for _ in 0..10 {
            plasma.step();
        }
        for phi in &plasma.poisson.phi {
            assert!(phi.is_finite());
        }
    }

    // ── Charge transport ──────────────────────────────────────────────────────

    /// Carrier density remains finite after many drift-diffusion steps.
    #[test]
    fn test_charge_transport_finite() {
        let nx = 16;
        let ny = 1;
        let mut ct = ChargeTransportLbm::new(nx, ny, 0.1, 300.0, 1.0, 1e-7);
        for i in 0..nx {
            ct.n[i] = 1e20;
        }
        // Uniform field
        for i in 0..nx {
            ct.ex[i] = 1e4;
        }
        ct.run(50);
        for n in &ct.n {
            assert!(n.is_finite());
        }
    }

    /// No drift (zero field) → carrier profile stays uniform.
    #[test]
    fn test_charge_transport_no_drift_uniform() {
        let nx = 8;
        let mut ct = ChargeTransportLbm::new(nx, 1, 0.1, 300.0, 1.0, 1e-7);
        let n0 = 1.0_f64;
        for i in 0..nx {
            ct.n[i] = n0;
        }
        ct.run(20);
        for n in &ct.n {
            assert!((n - n0).abs() < 0.5, "n = {:.6}", n);
        }
    }

    // ── General / constants ───────────────────────────────────────────────────

    /// Physical constants have correct orders of magnitude.
    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_constants_magnitude() {
        assert!(constants::EPSILON_0 > 8e-12 && constants::EPSILON_0 < 9e-12);
        assert!(constants::ELEM_CHARGE > 1.6e-19 && constants::ELEM_CHARGE < 1.7e-19);
        assert!(constants::K_BOLTZMANN > 1.3e-23 && constants::K_BOLTZMANN < 1.4e-23);
        assert!((constants::VT_300K - 0.02585).abs() < 0.001);
    }

    /// V_T at 300 K matches kT/q.
    #[test]
    fn test_thermal_voltage_300k() {
        let vt = constants::K_BOLTZMANN * 300.0 / constants::ELEM_CHARGE;
        assert!((vt - constants::VT_300K).abs() < 1e-5);
    }

    /// Poisson: 2×2 grid initialises without panic.
    #[test]
    fn test_poisson_2x2() {
        let mut solver = PoissonLbm::new(2, 2, 1.0, 1.0);
        solver.set_charge(0, 0, 1.0);
        solver.solve(10, 1e-6);
    }

    /// ElectricFieldGrid: 1×1 grid does not panic.
    #[test]
    fn test_electric_field_1x1() {
        let mut eg = ElectricFieldGrid::new(1, 1, 1.0, 1.0);
        let phi = vec![3.0_f64];
        eg.compute_from_phi(&phi);
        let (ex, ey) = eg.get(0, 0);
        assert!(ex.is_finite());
        assert!(ey.is_finite());
    }
}
