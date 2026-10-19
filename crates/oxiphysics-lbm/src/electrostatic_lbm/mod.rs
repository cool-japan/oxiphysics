// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Electrostatic LBM: Poisson equation solver, charge transport.
//!
//! Implements LBM-based solvers for:
//! - Poisson equation ∇²φ = -ρ/ε (electrostatics) with Dirichlet cells
//! - Drift-diffusion (charge transport for electrons/holes) with the
//!   Einstein relation D = μ·V_T fixing the physical time step
//! - Coupled Poisson + drift-diffusion semiconductor model (p-n junction)
//! - Electric field computation from potential
//! - Dielectric materials with permittivity tensor
//! - Space charge transport (Nernst-Planck equations)
//! - Coulomb explosion dynamics
//! - Two-fluid electrostatic plasma LBM
//!
//! # Lattice ↔ physical units
//!
//! All transport solvers use a D2Q5 BGK lattice whose lattice diffusivity is
//! `D_lu = (1/ω − ½)/3` (lattice speed of sound² = 1/3).  A physical
//! diffusivity `D` on a grid of spacing `Δx` is matched by the time step
//! `Δt = D_lu·Δx²/D`; a physical drift velocity `v` then becomes the lattice
//! velocity `u = v·Δt/Δx`.  For charge carriers `D = μ·V_T`, so the lattice
//! drift `u = s·E·Δx·D_lu/V_T` (s = sign of the carrier charge) depends only on
//! the potential drop per cell in units of the thermal voltage.
//!
//! # Exponentially fitted drift equilibrium
//!
//! Carrier and ion transport use moving equilibrium populations built from
//! `e^{±a/2}` (the rest population takes the remainder), fitted so that the
//! zero-current steady state of the BGK chain reproduces the Boltzmann factor
//! of the potential profile exactly, where `a` is the potential-energy drop
//! per cell in units of kT.  This is the lattice analogue of
//! Scharfetter–Gummel exponential fitting: the plain linear equilibrium
//! `w·n·(1 + 3c·u)` accumulates an `O(a³)` error per cell, which ruins
//! minority-carrier densities across a depletion region.  For weak fields
//! the fitted equilibrium reduces to the linear one.

mod plasma;
mod poisson;
mod semiconductor;
mod space_charge;
mod transport;

#[cfg(test)]
mod regression_tests;
#[cfg(test)]
mod tests;

pub use plasma::PlasmaLbm;
pub use poisson::PoissonLbm;
pub use semiconductor::SemiconductorLbm;
pub use space_charge::SpaceChargeLbm;
pub use transport::ChargeTransportLbm;

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
    /// Avogadro constant (1/mol).
    pub const AVOGADRO: f64 = 6.022_140_76e23;
    /// Faraday constant F = e·N_A (C/mol).
    pub const FARADAY: f64 = ELEM_CHARGE * AVOGADRO;
    /// Relative permittivity of crystalline silicon.
    pub const SILICON_EPSILON_R: f64 = 11.7;
}

// ─── D2Q5 stencil weights and velocity vectors ───────────────────────────────

/// D2Q5 lattice weights: \[center, E, N, W, S\].
const W5: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

/// D2Q5 velocity vectors: cx\[i\], cy\[i\].
const CX5: [i32; 5] = [0, 1, 0, -1, 0];
const CY5: [i32; 5] = [0, 0, 1, 0, -1];

/// Lattice diffusivity `D_lu = (1/ω − ½)/3` of a D2Q5 BGK scheme.
#[inline]
fn lattice_diffusivity(omega: f64) -> f64 {
    (1.0 / omega - 0.5) / 3.0
}

/// Periodic D2Q5 streaming of `f` (layout `cell * 5 + q`) into `out`.
fn stream_periodic(f: &[f64], nx: usize, ny: usize, out: &mut [f64]) {
    let nxi = nx as i64;
    let nyi = ny as i64;
    for y in 0..ny {
        for x in 0..nx {
            let src = (y * nx + x) * 5;
            for q in 0..5_usize {
                let xn = (x as i64 + i64::from(CX5[q])).rem_euclid(nxi) as usize;
                let yn = (y as i64 + i64::from(CY5[q])).rem_euclid(nyi) as usize;
                out[(yn * nx + xn) * 5 + q] = f[src + q];
            }
        }
    }
}

/// Standard D2Q5 weight of each moving population (rest weight 1/3).
const STANDARD_MOVING_WEIGHT: f64 = 1.0 / 6.0;

/// Largest admissible weight of each moving population (rest weight 0).
const MAX_MOVING_WEIGHT: f64 = 0.25;

/// Lattice diffusivity `D_lu = 2w·(1/ω − ½)` of a D2Q5 BGK scheme whose four
/// moving populations carry the weight `w` (`c_s² = 2w`).
#[inline]
fn weighted_lattice_diffusivity(omega: f64, weight: f64) -> f64 {
    2.0 * weight * (1.0 / omega - 0.5)
}

/// Scale (hx, hy) down uniformly so that `cosh hx + cosh hy ≤ c_max`, which
/// keeps the rest population non-negative.
fn limit_half_drifts(hx: f64, hy: f64, c_max: f64) -> (f64, f64) {
    if hx.cosh() + hy.cosh() <= c_max {
        return (hx, hy);
    }
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if (mid * hx).cosh() + (mid * hy).cosh() <= c_max {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo * hx, lo * hy)
}

/// Exponentially fitted D2Q5 drift equilibrium.
///
/// `n` is the density, (ax, ay) the drift strength per cell in thermal units
/// (for a carrier of charge sign s: `a = s·E·Δx/V_T`, the potential-energy
/// drop per cell in units of kT), ω the BGK rate and `w` the weight of each
/// moving population.  With `h = a/2` and `m = 1 − ω` the moving populations
/// along each axis are
///
/// ```text
/// f₊ = n·w·(e^{h} − m·e^{−h})/ω,   f₋ = n·w·(e^{−h} − m·e^{h})/ω
/// ```
///
/// and the rest population takes the remainder, so the zeroth moment is
/// exactly `n`.  For a BGK chain the zero-flux steady state links neighbours
/// through `P = f₊ + m·f₋ = n·w·e^{h}` and `M = f₋ + m·f₊ = n·w·e^{−h}`, hence
/// `n_{x+1}/n_x = exp((a_x + a_{x+1})/2)`: the accumulated density ratio
/// reproduces the Boltzmann factor of the potential profile exactly and
/// pointwise, for every ω and w.  The first moment `2w·(2 − ω)·n·sinh(h)/(2ω)`
/// reduces to `n·D_lu·a` (`D_lu = 2w(1/ω − ½)`) for weak fields, i.e. the
/// drift velocity `u = v·Δt/Δx`.  Drops that would make a population
/// negative (beyond `−V_T·ln(1 − ω)` for ω < 1, or a rest population below
/// zero) are limited.
fn drift_equilibrium(n: f64, ax: f64, ay: f64, omega: f64, weight: f64) -> [f64; 5] {
    if weight <= 0.0 {
        return [n, 0.0, 0.0, 0.0, 0.0];
    }
    let m = 1.0 - omega;
    // rest ≥ 0  ⇔  cosh hx + cosh hy ≤ 1/(2w)
    let c_max = 0.5 / weight;
    let axis_cap = (c_max - 1.0).max(1.0).acosh();
    let h_cap = if m > 0.0 {
        (-0.5 * m.ln()).clamp(0.0, axis_cap)
    } else {
        axis_cap
    };
    let (hx, hy) = limit_half_drifts(
        (0.5 * ax).max(-h_cap).min(h_cap),
        (0.5 * ay).max(-h_cap).min(h_cap),
        c_max,
    );
    let scale = n * weight / omega;
    let pair = |h: f64| {
        let (ep, em) = (h.exp(), (-h).exp());
        (scale * (ep - m * em), scale * (em - m * ep))
    };
    let (east, west) = pair(hx);
    let (north, south) = pair(hy);
    [n - (east + west + north + south), east, north, west, south]
}

/// Neighbour indices and stencil span along one axis.
///
/// Periodic axes use wrapped central differences; non-periodic axes use
/// central differences in the interior and one-sided differences at the two
/// edges.  The returned span is the number of cells between the two
/// neighbours (never zero, so a single-cell axis yields a zero gradient).
#[inline]
fn axis_stencil(len: usize, k: usize, periodic: bool) -> (usize, usize, f64) {
    if len < 2 {
        return (k, k, 1.0);
    }
    if periodic && len >= 3 {
        let kp = (k + 1) % len;
        let km = (k + len - 1) % len;
        return (kp, km, 2.0);
    }
    let kp = if k + 1 < len { k + 1 } else { k };
    let km = if k > 0 { k - 1 } else { k };
    (kp, km, (kp - km) as f64)
}

/// Electric field **E** = −∇φ at cell (x, y) of an `nx × ny` potential array
/// with uniform spacing `dx`.
fn neg_gradient(
    phi: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    x: usize,
    y: usize,
    periodic: bool,
) -> (f64, f64) {
    let (xp, xm, sx) = axis_stencil(nx, x, periodic);
    let (yp, ym, sy) = axis_stencil(ny, y, periodic);
    let ex = -(phi[y * nx + xp] - phi[y * nx + xm]) / (sx * dx);
    let ey = -(phi[yp * nx + x] - phi[ym * nx + x]) / (sy * dx);
    (ex, ey)
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
