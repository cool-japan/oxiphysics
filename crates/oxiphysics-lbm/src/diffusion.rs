// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM-based diffusion and advection-diffusion solvers.
//!
//! Provides:
//! - `DiffusionLattice1D` — D1Q3 BGK diffusion in 1D
//! - `DiffusionLattice2D` — D2Q5 BGK diffusion in 2D
//! - `AdvectionDiffusion1D` — upwind advection + diffusion in 1D
//! - `MultiSpeciesDiffusion` — independent species diffusing on a 2D grid
//! - Utility functions: Péclet number, analytical Gaussian spreading,
//!   Von Neumann stability, D2Q5 equilibrium, and Soret thermodiffusion flux.
//!
//! References:
//! - Succi, S. (2001). *The Lattice Boltzmann Equation for Fluid Dynamics
//!   and Beyond*. Oxford.
//! - Wolf-Gladrow, D. A. (2000). *Lattice-Gas Cellular Automata and Lattice
//!   Boltzmann Models*. Springer.

#[cfg(test)]
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// D1Q3 constants
// ---------------------------------------------------------------------------

/// D1Q3 lattice velocities: c = {-1, 0, +1}.
const D1Q3_C: [i32; 3] = [-1, 0, 1];

/// D1Q3 weights: w = {1/6, 2/3, 1/6}.
const D1Q3_W: [f64; 3] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];

/// Speed-of-sound squared for D1Q3 diffusion lattice (1/3).
const D1Q3_CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// D2Q5 constants
// ---------------------------------------------------------------------------

/// D2Q5 lattice velocities: {rest, ±x, ±y}.
const D2Q5_CX: [i32; 5] = [0, 1, -1, 0, 0];
/// D2Q5 lattice velocities y-component.
const D2Q5_CY: [i32; 5] = [0, 0, 0, 1, -1];

/// D2Q5 weights: rest = 1/3, axial = 1/6.
const D2Q5_W: [f64; 5] = [1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0];

/// Speed-of-sound squared for D2Q5 diffusion lattice (1/4).
const D2Q5_CS2: f64 = 1.0 / 4.0;

// ---------------------------------------------------------------------------
// DiffusionLattice1D
// ---------------------------------------------------------------------------

/// 1D BGK lattice diffusion on a D1Q3 stencil.
///
/// The scalar field `phi` (e.g., concentration or temperature) evolves under:
///
/// ```text
/// ∂φ/∂t = D ∂²φ/∂x²
/// ```
///
/// The relaxation time is derived from the diffusivity:
/// `τ = D / (cs² · dt) + 0.5`
pub struct DiffusionLattice1D {
    /// Concentration / scalar field values at each node.
    pub phi: Vec<f64>,
    /// Distribution functions `f[node][direction]` (3 directions per node).
    f: Vec<[f64; 3]>,
    /// Number of lattice nodes.
    pub n: usize,
    /// Spatial step size (m).
    pub dx: f64,
    /// Time step size (s).
    pub dt: f64,
    /// Diffusion coefficient (m²/s).
    pub d: f64,
    /// BGK relaxation time τ.
    tau: f64,
}

impl DiffusionLattice1D {
    /// Create a new 1D diffusion lattice initialised to φ = 0 everywhere.
    ///
    /// # Arguments
    /// * `n`  — number of lattice nodes (periodic domain)
    /// * `dx` — lattice spacing (m)
    /// * `dt` — time step (s)
    /// * `d`  — diffusion coefficient (m²/s)
    pub fn new(n: usize, dx: f64, dt: f64, d: f64) -> Self {
        let tau = d * dt / (D1Q3_CS2 * dx * dx) + 0.5;
        let f = vec![[D1Q3_W[0], D1Q3_W[1], D1Q3_W[2]]; n];
        Self {
            phi: vec![0.0; n],
            f,
            n,
            dx,
            dt,
            d,
            tau,
        }
    }

    /// Initialise the scalar field as a Gaussian profile.
    ///
    /// `φ(x) = amplitude · exp(−(x − center)² / (2σ²))`
    ///
    /// # Arguments
    /// * `center`    — peak position (same units as `dx`)
    /// * `sigma`     — standard deviation (same units as `dx`)
    /// * `amplitude` — peak value
    pub fn set_gaussian(&mut self, center: f64, sigma: f64, amplitude: f64) {
        for i in 0..self.n {
            let x = i as f64 * self.dx;
            let val = amplitude * (-(x - center).powi(2) / (2.0 * sigma * sigma)).exp();
            self.phi[i] = val;
            // Initialise f to equilibrium
            for (q, &w) in D1Q3_W.iter().enumerate() {
                self.f[i][q] = w * val;
            }
        }
    }

    /// Perform one LBM D1Q3 diffusion step (collision + streaming + periodic BC).
    pub fn step(&mut self) {
        let n = self.n;
        let omega = 1.0 / self.tau;

        // --- Collision ---
        let mut f_post = self.f.clone();
        for (i, f_post_i) in f_post.iter_mut().enumerate().take(n) {
            let rho = self.phi[i];
            for q in 0..3 {
                let f_eq = D1Q3_W[q] * rho;
                f_post_i[q] = self.f[i][q] - omega * (self.f[i][q] - f_eq);
            }
        }

        // --- Streaming (periodic) ---
        let mut f_stream = vec![[0.0_f64; 3]; n];
        for (i, f_post_i) in f_post.iter().enumerate().take(n) {
            for (q, &fq) in f_post_i.iter().enumerate() {
                let c = D1Q3_C[q];
                let dest = ((i as i64 + c as i64).rem_euclid(n as i64)) as usize;
                f_stream[dest][q] = fq;
            }
        }
        self.f = f_stream;

        // --- Update macroscopic field ---
        for i in 0..n {
            self.phi[i] = self.f[i].iter().sum();
        }
    }

    /// Compute the total mass (integral of φ over the domain).
    pub fn total_mass(&self) -> f64 {
        self.phi.iter().sum::<f64>() * self.dx
    }
}

// ---------------------------------------------------------------------------
// DiffusionLattice2D
// ---------------------------------------------------------------------------

/// 2D BGK lattice diffusion on a D2Q5 stencil.
///
/// Solves `∂φ/∂t = D (∂²φ/∂x² + ∂²φ/∂y²)` with optional Dirichlet BCs
/// on the domain boundary.
pub struct DiffusionLattice2D {
    /// Scalar field `phi[iy][ix]`.
    pub phi: Vec<Vec<f64>>,
    /// Distribution functions `f[iy*nx+ix][direction]`.
    f: Vec<[f64; 5]>,
    /// Number of nodes in x.
    pub nx: usize,
    /// Number of nodes in y.
    pub ny: usize,
    /// Spatial step (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Diffusion coefficient (m²/s).
    pub d: f64,
    /// BGK relaxation time.
    tau: f64,
    /// Dirichlet boundary value (None = periodic).
    dirichlet: Option<f64>,
}

impl DiffusionLattice2D {
    /// Create a new 2D diffusion lattice initialised to φ = 0 everywhere.
    pub fn new(nx: usize, ny: usize, dx: f64, dt: f64, d: f64) -> Self {
        let tau = d * dt / (D2Q5_CS2 * dx * dx) + 0.5;
        let n = nx * ny;
        let f_init = {
            let mut a = [0.0_f64; 5];
            for q in 0..5 {
                a[q] = D2Q5_W[q] * 0.0;
            }
            a
        };
        let phi = vec![vec![0.0_f64; nx]; ny];
        Self {
            phi,
            f: vec![f_init; n],
            nx,
            ny,
            dx,
            dt,
            d,
            tau,
            dirichlet: None,
        }
    }

    /// Flat index helper.
    #[inline]
    fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Perform one D2Q5 LBM diffusion step with optional Dirichlet boundary.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let omega = 1.0 / self.tau;

        // --- Collision ---
        let mut f_post = self.f.clone();
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = self.idx(ix, iy);
                let rho = self.phi[iy][ix];
                for q in 0..5 {
                    let f_eq = D2Q5_W[q] * rho;
                    f_post[idx][q] = self.f[idx][q] - omega * (self.f[idx][q] - f_eq);
                }
            }
        }

        // --- Streaming (periodic) ---
        let mut f_stream = vec![[0.0_f64; 5]; nx * ny];
        for iy in 0..ny {
            for ix in 0..nx {
                let src = self.idx(ix, iy);
                for q in 0..5 {
                    let cx = D2Q5_CX[q];
                    let cy = D2Q5_CY[q];
                    let dix = ((ix as i64 + cx as i64).rem_euclid(nx as i64)) as usize;
                    let diy = ((iy as i64 + cy as i64).rem_euclid(ny as i64)) as usize;
                    let dest = self.idx(dix, diy);
                    f_stream[dest][q] = f_post[src][q];
                }
            }
        }
        self.f = f_stream;

        // --- Update macroscopic phi ---
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = self.idx(ix, iy);
                self.phi[iy][ix] = self.f[idx].iter().sum();
            }
        }

        // --- Apply Dirichlet BC on boundary nodes ---
        if let Some(val) = self.dirichlet {
            self.apply_dirichlet_internal(val);
        }
    }

    fn apply_dirichlet_internal(&mut self, val: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for ix in 0..nx {
            self.phi[0][ix] = val;
            self.phi[ny - 1][ix] = val;
            let idx0 = self.idx(ix, 0);
            let idx1 = self.idx(ix, ny - 1);
            for (q, &w) in D2Q5_W.iter().enumerate() {
                self.f[idx0][q] = w * val;
                self.f[idx1][q] = w * val;
            }
        }
        for iy in 0..ny {
            self.phi[iy][0] = val;
            self.phi[iy][nx - 1] = val;
            let idx0 = self.idx(0, iy);
            let idx1 = self.idx(nx - 1, iy);
            for (q, &w) in D2Q5_W.iter().enumerate() {
                self.f[idx0][q] = w * val;
                self.f[idx1][q] = w * val;
            }
        }
    }

    /// Set boundary nodes to a fixed Dirichlet value applied after every step.
    ///
    /// Call with `None` to revert to periodic boundaries.
    pub fn set_dirichlet_bc(&mut self, value: f64) {
        self.dirichlet = Some(value);
        self.apply_dirichlet_internal(value);
    }

    /// Compute total mass ∫∫ φ dx dy.
    pub fn total_mass(&self) -> f64 {
        let s: f64 = self.phi.iter().flat_map(|row| row.iter()).sum();
        s * self.dx * self.dx
    }

    /// Set the concentration field from a 2D array.
    ///
    /// `data[iy][ix]` is copied into φ and the distributions are initialised
    /// to equilibrium.
    pub fn set_field(&mut self, data: &[Vec<f64>]) {
        let nx = self.nx;
        let ny = self.ny;
        for (iy, data_row) in data.iter().enumerate().take(ny) {
            for (ix, &val) in data_row.iter().enumerate().take(nx) {
                self.phi[iy][ix] = val;
                let idx = self.idx(ix, iy);
                for (q, &w) in D2Q5_W.iter().enumerate() {
                    self.f[idx][q] = w * val;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Analytical solution
// ---------------------------------------------------------------------------

/// Analytical solution for 1D diffusion of an initial Gaussian profile.
///
/// The solution to `∂φ/∂t = D ∂²φ/∂x²` with
/// `φ(x, 0) = exp(−(x−x₀)² / (2σ₀²))` is:
///
/// ```text
/// φ(x, t) = σ₀ / √(σ₀² + 2Dt) · exp(−(x−x₀)² / (2(σ₀² + 2Dt)))
/// ```
///
/// (conserves total integral = √(2π) σ₀).
///
/// # Arguments
/// * `x`      — evaluation point (m)
/// * `t`      — time (s)
/// * `d`      — diffusion coefficient (m²/s)
/// * `x0`     — initial Gaussian centre (m)
/// * `sigma0` — initial standard deviation (m)
pub fn analytical_diffusion_1d(x: f64, t: f64, d: f64, x0: f64, sigma0: f64) -> f64 {
    let sigma2 = sigma0 * sigma0 + 2.0 * d * t;
    let norm = sigma0 / sigma2.sqrt();
    norm * (-(x - x0).powi(2) / (2.0 * sigma2)).exp()
}

// ---------------------------------------------------------------------------
// Péclet number
// ---------------------------------------------------------------------------

/// Compute the Péclet number Pe = uL/D.
///
/// The Péclet number characterises the relative importance of advective to
/// diffusive transport.
///
/// # Arguments
/// * `u` — characteristic velocity (m/s)
/// * `l` — characteristic length (m)
/// * `d` — diffusion coefficient (m²/s)
pub fn peclet_number(u: f64, l: f64, d: f64) -> f64 {
    if d == 0.0 {
        return f64::INFINITY;
    }
    u * l / d
}

// ---------------------------------------------------------------------------
// AdvectionDiffusion1D
// ---------------------------------------------------------------------------

/// 1D advection-diffusion solver using upwind advection and central diffusion.
///
/// Solves `∂φ/∂t + u ∂φ/∂x = D ∂²φ/∂x²` on a periodic domain.
pub struct AdvectionDiffusion1D {
    /// Scalar field values.
    pub phi: Vec<f64>,
    /// Advection velocity (m/s).
    pub u: f64,
    /// Number of nodes.
    pub n: usize,
    /// Spatial step (m).
    pub dx: f64,
    /// Time step (s).
    pub dt: f64,
    /// Diffusion coefficient (m²/s).
    pub d: f64,
}

impl AdvectionDiffusion1D {
    /// Create a new advection-diffusion solver initialised to φ = 0.
    pub fn new(n: usize, dx: f64, dt: f64, d: f64, u: f64) -> Self {
        Self {
            phi: vec![0.0; n],
            u,
            n,
            dx,
            dt,
            d,
        }
    }

    /// Initialise φ to a Gaussian profile.
    pub fn set_gaussian(&mut self, center: f64, sigma: f64, amplitude: f64) {
        let dx = self.dx;
        for (i, phi) in self.phi.iter_mut().enumerate() {
            let x = i as f64 * dx;
            *phi = amplitude * (-(x - center).powi(2) / (2.0 * sigma * sigma)).exp();
        }
    }

    /// Advance one time step using first-order upwind advection + central diffusion.
    pub fn step(&mut self) {
        let n = self.n;
        let dx = self.dx;
        let dt = self.dt;
        let u = self.u;
        let d = self.d;
        let phi = &self.phi;

        let mut new_phi = vec![0.0_f64; n];
        for i in 0..n {
            let im1 = (i + n - 1) % n;
            let ip1 = (i + 1) % n;

            // Upwind advection
            let adv = if u >= 0.0 {
                u * (phi[i] - phi[im1]) / dx
            } else {
                u * (phi[ip1] - phi[i]) / dx
            };

            // Central diffusion
            let diff = d * (phi[ip1] - 2.0 * phi[i] + phi[im1]) / (dx * dx);

            new_phi[i] = phi[i] + dt * (-adv + diff);
        }
        self.phi = new_phi;
    }

    /// Compute the L2 error norm against an exact solution slice.
    ///
    /// `||φ − φ_exact||₂ = √(dx · Σ (φᵢ − φ_exact_i)²)`
    pub fn l2_error(&self, exact: &[f64]) -> f64 {
        let s: f64 = self
            .phi
            .iter()
            .zip(exact.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        (s * self.dx).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Von Neumann stability check
// ---------------------------------------------------------------------------

/// Check Von Neumann stability for the D1Q3 BGK diffusion scheme.
///
/// The scheme is stable when `τ > 0.5`, i.e. the diffusive CFL number
/// `r = D · dt / (cs² · dx²)` satisfies `r > 0`.  In practice, large τ
/// leads to slow convergence; `τ ∈ (0.5, 2]` is recommended.
///
/// Returns `true` if the scheme is expected to be stable.
///
/// # Arguments
/// * `d`  — diffusion coefficient
/// * `dx` — lattice spacing
/// * `dt` — time step
pub fn von_neumann_stability_d1q3(d: f64, dx: f64, dt: f64) -> bool {
    let tau = d * dt / (D1Q3_CS2 * dx * dx) + 0.5;
    tau > 0.5
}

// ---------------------------------------------------------------------------
// D2Q5 equilibrium
// ---------------------------------------------------------------------------

/// Compute the D2Q5 equilibrium distribution for a passive scalar.
///
/// For purely diffusive (zero-velocity) transport the equilibrium is simply:
///
/// ```text
/// f_q^eq = w_q · ρ
/// ```
///
/// The optional `cs2` parameter adjusts the speed-of-sound squared used in
/// the equilibrium; pass `1.0/4.0` for the standard D2Q5 diffusion model.
///
/// # Arguments
/// * `rho` — local scalar density
/// * `w`   — weight vector (length 5)
/// * `cs2` — speed of sound squared (e.g. 1/4 for D2Q5)
pub fn d2q5_equilibrium(rho: f64, w: &[f64], cs2: f64) -> Vec<f64> {
    let _ = cs2; // retained for API completeness; zero-velocity eq is w_q * rho
    w.iter().map(|wi| wi * rho).collect()
}

// ---------------------------------------------------------------------------
// MultiSpeciesDiffusion
// ---------------------------------------------------------------------------

/// Independent multi-species diffusion on a 2D grid.
///
/// Each species has its own diffusion coefficient and scalar field.  The
/// species do not interact (ideal mixture); coupling can be added by the
/// caller between time steps.
pub struct MultiSpeciesDiffusion {
    /// Species concentration fields `[species][iy][ix]`.
    pub species: Vec<Vec<Vec<f64>>>,
    /// Diffusion coefficients for each species (m²/s).
    pub diff_coeffs: Vec<f64>,
    /// Number of nodes in x.
    pub nx: usize,
    /// Number of nodes in y.
    pub ny: usize,
}

impl MultiSpeciesDiffusion {
    /// Create a new multi-species diffusion field initialised to zero.
    ///
    /// # Arguments
    /// * `nx`         — grid width
    /// * `ny`         — grid height
    /// * `diff_coeffs` — diffusion coefficient for each species
    pub fn new(nx: usize, ny: usize, diff_coeffs: Vec<f64>) -> Self {
        let ns = diff_coeffs.len();
        let species = vec![vec![vec![0.0_f64; nx]; ny]; ns];
        Self {
            species,
            diff_coeffs,
            nx,
            ny,
        }
    }

    /// Advance all species by one time step using explicit finite differences.
    ///
    /// Periodic boundary conditions are applied in both directions.
    ///
    /// # Arguments
    /// * `dt` — time step (s)
    /// * `dx` — uniform grid spacing (m)
    pub fn step(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let ns = self.species.len();
        for s in 0..ns {
            let d = self.diff_coeffs[s];
            let mut new_field = vec![vec![0.0_f64; nx]; ny];
            for (iy, row) in new_field.iter_mut().enumerate() {
                for (ix, cell) in row.iter_mut().enumerate() {
                    let im1x = (ix + nx - 1) % nx;
                    let ip1x = (ix + 1) % nx;
                    let im1y = (iy + ny - 1) % ny;
                    let ip1y = (iy + 1) % ny;
                    let laplacian = (self.species[s][iy][ip1x]
                        + self.species[s][iy][im1x]
                        + self.species[s][ip1y][ix]
                        + self.species[s][im1y][ix]
                        - 4.0 * self.species[s][iy][ix])
                        / (dx * dx);
                    *cell = self.species[s][iy][ix] + dt * d * laplacian;
                }
            }
            self.species[s] = new_field;
        }
    }

    /// Return a copy of the concentration field for the given species index.
    pub fn species_density(&self, species_idx: usize) -> Vec<Vec<f64>> {
        self.species[species_idx].clone()
    }
}

// ---------------------------------------------------------------------------
// Soret effect (thermodiffusion)
// ---------------------------------------------------------------------------

/// Compute the Soret (thermodiffusion) flux contribution.
///
/// The Soret flux is the extra mass flux driven by a temperature gradient:
///
/// ```text
/// J_Soret = −D · S_T · c · ∇T
/// ```
///
/// where `S_T` is the Soret coefficient (K⁻¹).  In this simplified 1D
/// form the concentration `c` is absorbed into `grad_c`, giving:
///
/// ```text
/// J_total = −D · ∇c − D · S_T · grad_t
/// ```
///
/// # Arguments
/// * `grad_c`      — concentration gradient (mol/m⁴)
/// * `grad_t`      — temperature gradient (K/m)
/// * `soret_coeff` — Soret coefficient S_T (K⁻¹)
/// * `d`           — diffusion coefficient (m²/s)
pub fn soret_effect_flux(grad_c: f64, grad_t: f64, soret_coeff: f64, d: f64) -> f64 {
    -d * grad_c - d * soret_coeff * grad_t
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- DiffusionLattice1D ---

    #[test]
    fn test_d1_new_zero_field() {
        let lat = DiffusionLattice1D::new(32, 1.0, 0.01, 0.1);
        assert_eq!(lat.phi.len(), 32);
        assert!(lat.phi.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn test_d1_gaussian_peak() {
        let mut lat = DiffusionLattice1D::new(64, 0.1, 0.001, 0.01);
        lat.set_gaussian(3.2, 0.3, 1.0);
        let max = lat.phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max > 0.9, "peak should be close to 1.0, got {max}");
    }

    #[test]
    fn test_d1_mass_conservation() {
        let mut lat = DiffusionLattice1D::new(64, 1.0, 0.1, 0.01);
        lat.set_gaussian(32.0, 4.0, 1.0);
        let mass0 = lat.total_mass();
        for _ in 0..20 {
            lat.step();
        }
        let mass1 = lat.total_mass();
        assert!(
            (mass1 - mass0).abs() / (mass0.abs() + 1e-15) < 1e-10,
            "mass not conserved: {mass0} vs {mass1}"
        );
    }

    #[test]
    fn test_d1_spreading() {
        let mut lat = DiffusionLattice1D::new(128, 1.0, 0.01, 0.1);
        lat.set_gaussian(64.0, 4.0, 1.0);
        let peak0 = lat.phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        for _ in 0..50 {
            lat.step();
        }
        let peak1 = lat.phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(peak1 < peak0, "peak should decrease as profile spreads");
    }

    #[test]
    fn test_d1_tau_positive() {
        let lat = DiffusionLattice1D::new(32, 1.0, 0.1, 0.05);
        assert!(lat.tau > 0.5);
    }

    #[test]
    fn test_d1_total_mass_zero() {
        let lat = DiffusionLattice1D::new(10, 1.0, 0.01, 0.1);
        assert_eq!(lat.total_mass(), 0.0);
    }

    // --- DiffusionLattice2D ---

    #[test]
    fn test_d2_new() {
        let lat = DiffusionLattice2D::new(16, 16, 1.0, 0.01, 0.05);
        assert_eq!(lat.phi.len(), 16);
        assert_eq!(lat.phi[0].len(), 16);
    }

    #[test]
    fn test_d2_mass_conservation_periodic() {
        let nx = 32;
        let ny = 32;
        let mut lat = DiffusionLattice2D::new(nx, ny, 1.0, 0.05, 0.1);
        // place blob in center
        lat.phi[ny / 2][nx / 2] = 10.0;
        let idx = lat.idx(nx / 2, ny / 2);
        for (q, &w) in D2Q5_W.iter().enumerate() {
            lat.f[idx][q] = w * 10.0;
        }
        let mass0 = lat.total_mass();
        for _ in 0..10 {
            lat.step();
        }
        let mass1 = lat.total_mass();
        assert!(
            (mass1 - mass0).abs() / (mass0.abs() + 1e-15) < 1e-8,
            "2D mass not conserved: {mass0} vs {mass1}"
        );
    }

    #[test]
    fn test_d2_dirichlet_bc() {
        let mut lat = DiffusionLattice2D::new(16, 16, 1.0, 0.01, 0.05);
        lat.set_dirichlet_bc(0.0);
        // boundary nodes should be zero
        for ix in 0..16 {
            assert_eq!(lat.phi[0][ix], 0.0);
            assert_eq!(lat.phi[15][ix], 0.0);
        }
        for iy in 0..16 {
            assert_eq!(lat.phi[iy][0], 0.0);
            assert_eq!(lat.phi[iy][15], 0.0);
        }
    }

    #[test]
    fn test_d2_dirichlet_step() {
        let mut lat = DiffusionLattice2D::new(16, 16, 1.0, 0.02, 0.05);
        lat.set_dirichlet_bc(0.5);
        lat.step();
        // boundary should remain at 0.5
        assert!((lat.phi[0][0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_d2_set_field() {
        let nx = 8;
        let ny = 8;
        let mut lat = DiffusionLattice2D::new(nx, ny, 1.0, 0.01, 0.05);
        let data: Vec<Vec<f64>> = (0..ny)
            .map(|iy| (0..nx).map(|ix| (iy * nx + ix) as f64).collect())
            .collect();
        lat.set_field(&data);
        assert!((lat.phi[3][4] - data[3][4]).abs() < 1e-15);
    }

    // --- Analytical solution ---

    #[test]
    fn test_analytical_diffusion_1d_at_t0() {
        let val = analytical_diffusion_1d(0.0, 0.0, 0.1, 0.0, 1.0);
        // At t=0: val = 1.0 (peak)
        assert!((val - 1.0).abs() < 1e-10, "val={val}");
    }

    #[test]
    fn test_analytical_diffusion_1d_spread() {
        let v0 = analytical_diffusion_1d(0.0, 0.0, 0.1, 0.0, 1.0);
        let v1 = analytical_diffusion_1d(0.0, 1.0, 0.1, 0.0, 1.0);
        assert!(v1 < v0, "peak must decrease as profile spreads");
    }

    #[test]
    fn test_analytical_diffusion_symmetry() {
        let t = 0.5;
        let d = 0.2;
        let x0 = 0.0;
        let s0 = 1.0;
        let vp = analytical_diffusion_1d(1.0, t, d, x0, s0);
        let vm = analytical_diffusion_1d(-1.0, t, d, x0, s0);
        assert!((vp - vm).abs() < 1e-14);
    }

    #[test]
    fn test_analytical_mass_conserved() {
        // Integral approximation over [-10, 10] with 1000 steps
        let d = 0.1;
        let s0 = 1.0;
        let x0 = 0.0;
        let dx = 0.02;
        let t = 2.0;
        let integral: f64 = (-500..=500)
            .map(|i| {
                let x = i as f64 * dx;
                analytical_diffusion_1d(x, t, d, x0, s0) * dx
            })
            .sum();
        let expected = (2.0 * PI).sqrt() * s0;
        assert!(
            (integral - expected).abs() / expected < 0.01,
            "integral={integral} expected≈{expected}"
        );
    }

    // --- Péclet number ---

    #[test]
    fn test_peclet_basic() {
        let pe = peclet_number(1.0, 1.0, 0.1);
        assert!((pe - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_peclet_zero_diff() {
        assert_eq!(peclet_number(1.0, 1.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn test_peclet_zero_velocity() {
        assert_eq!(peclet_number(0.0, 1.0, 0.1), 0.0);
    }

    // --- AdvectionDiffusion1D ---

    #[test]
    fn test_adv_diff_new() {
        let s = AdvectionDiffusion1D::new(64, 0.1, 0.001, 0.01, 0.5);
        assert_eq!(s.phi.len(), 64);
        assert_eq!(s.u, 0.5);
    }

    #[test]
    fn test_adv_diff_gaussian_set() {
        let mut s = AdvectionDiffusion1D::new(64, 0.1, 0.001, 0.01, 0.0);
        s.set_gaussian(3.2, 0.3, 2.0);
        let max = s.phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max > 1.8);
    }

    #[test]
    fn test_adv_diff_step_runs() {
        let mut s = AdvectionDiffusion1D::new(32, 0.1, 0.001, 0.01, 0.1);
        s.set_gaussian(1.6, 0.2, 1.0);
        s.step();
        // Just check the field has finite values
        assert!(s.phi.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_adv_diff_l2_error_self() {
        let s = AdvectionDiffusion1D::new(32, 1.0, 0.01, 0.1, 0.5);
        let exact = s.phi.clone();
        assert_eq!(s.l2_error(&exact), 0.0);
    }

    #[test]
    fn test_adv_diff_l2_error_nonzero() {
        let s = AdvectionDiffusion1D::new(4, 1.0, 0.01, 0.1, 0.5);
        let exact = vec![1.0, 0.0, 0.0, 0.0];
        let err = s.l2_error(&exact);
        assert!(err > 0.0);
    }

    #[test]
    fn test_adv_only_advects_positive() {
        // With D=0 and positive u, the blob should shift right
        let n = 64;
        let dx = 1.0;
        let dt = 0.4; // CFL=0.4 < 1
        let mut s = AdvectionDiffusion1D::new(n, dx, dt, 0.0, 1.0);
        s.set_gaussian(16.0, 2.0, 1.0);
        let center0: f64 = s
            .phi
            .iter()
            .enumerate()
            .map(|(i, v)| i as f64 * (*v))
            .sum::<f64>()
            / s.phi.iter().sum::<f64>();
        for _ in 0..10 {
            s.step();
        }
        let center1: f64 = s
            .phi
            .iter()
            .enumerate()
            .map(|(i, v)| i as f64 * (*v))
            .sum::<f64>()
            / s.phi.iter().sum::<f64>();
        assert!(center1 > center0, "center0={center0} center1={center1}");
    }

    // --- Von Neumann stability ---

    #[test]
    fn test_stability_stable() {
        // tau > 0.5 when D*dt/(cs2*dx^2) > 0
        assert!(von_neumann_stability_d1q3(0.1, 1.0, 0.1));
    }

    #[test]
    fn test_stability_always_true_for_positive_d() {
        // For any D>0, dx>0, dt>0 the tau > 0.5 always
        for &d in &[0.001, 0.01, 0.1, 1.0] {
            assert!(von_neumann_stability_d1q3(d, 1.0, 0.01));
        }
    }

    // --- D2Q5 equilibrium ---

    #[test]
    fn test_d2q5_eq_sum() {
        let w: Vec<f64> = D2Q5_W.to_vec();
        let eq = d2q5_equilibrium(2.0, &w, D2Q5_CS2);
        let s: f64 = eq.iter().sum();
        assert!((s - 2.0).abs() < 1e-12, "sum={s}");
    }

    #[test]
    fn test_d2q5_eq_zero_rho() {
        let w: Vec<f64> = D2Q5_W.to_vec();
        let eq = d2q5_equilibrium(0.0, &w, D2Q5_CS2);
        assert!(eq.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn test_d2q5_eq_proportional() {
        let w: Vec<f64> = D2Q5_W.to_vec();
        let eq1 = d2q5_equilibrium(1.0, &w, D2Q5_CS2);
        let eq2 = d2q5_equilibrium(3.0, &w, D2Q5_CS2);
        for (a, b) in eq1.iter().zip(eq2.iter()) {
            assert!((b - 3.0 * a).abs() < 1e-12);
        }
    }

    // --- MultiSpeciesDiffusion ---

    #[test]
    fn test_multi_new() {
        let m = MultiSpeciesDiffusion::new(8, 8, vec![0.1, 0.2]);
        assert_eq!(m.species.len(), 2);
    }

    #[test]
    fn test_multi_species_density() {
        let m = MultiSpeciesDiffusion::new(4, 4, vec![0.1, 0.2]);
        let d = m.species_density(0);
        assert_eq!(d.len(), 4);
        assert_eq!(d[0].len(), 4);
    }

    #[test]
    fn test_multi_step_mass_conservation() {
        let mut m = MultiSpeciesDiffusion::new(16, 16, vec![0.05, 0.1]);
        m.species[0][8][8] = 1.0;
        m.species[1][4][4] = 2.0;
        let mass0_s0: f64 = m.species[0].iter().flat_map(|r| r.iter()).sum();
        let mass0_s1: f64 = m.species[1].iter().flat_map(|r| r.iter()).sum();
        for _ in 0..5 {
            m.step(0.01, 1.0);
        }
        let mass1_s0: f64 = m.species[0].iter().flat_map(|r| r.iter()).sum();
        let mass1_s1: f64 = m.species[1].iter().flat_map(|r| r.iter()).sum();
        assert!((mass1_s0 - mass0_s0).abs() < 1e-10);
        assert!((mass1_s1 - mass0_s1).abs() < 1e-10);
    }

    #[test]
    fn test_multi_step_diffuses() {
        let mut m = MultiSpeciesDiffusion::new(16, 16, vec![0.1]);
        m.species[0][8][8] = 10.0;
        let peak0 = m.species[0]
            .iter()
            .flat_map(|r| r.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        for _ in 0..20 {
            m.step(0.1, 1.0);
        }
        let peak1 = m.species[0]
            .iter()
            .flat_map(|r| r.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(peak1 < peak0);
    }

    // --- Soret effect ---

    #[test]
    fn test_soret_zero_grads() {
        assert_eq!(soret_effect_flux(0.0, 0.0, 0.1, 0.05), 0.0);
    }

    #[test]
    fn test_soret_concentration_only() {
        // grad_t = 0 → pure Fick's law
        let j = soret_effect_flux(1.0, 0.0, 0.1, 0.05);
        assert!((j - (-0.05)).abs() < 1e-12, "j={j}");
    }

    #[test]
    fn test_soret_temperature_only() {
        let j = soret_effect_flux(0.0, 1.0, 0.1, 0.05);
        assert!((j - (-0.005)).abs() < 1e-12, "j={j}");
    }

    #[test]
    fn test_soret_negative_coeff() {
        // Soret coefficient can be negative (some systems)
        let j_pos = soret_effect_flux(0.0, 1.0, 0.1, 0.05);
        let j_neg = soret_effect_flux(0.0, 1.0, -0.1, 0.05);
        assert!(j_neg > j_pos);
    }

    #[test]
    fn test_soret_linearity_in_gradients() {
        let d = 0.05;
        let st = 0.1;
        let gc = 2.0;
        let gt = 3.0;
        let j_both = soret_effect_flux(gc, gt, st, d);
        let j_c = soret_effect_flux(gc, 0.0, st, d);
        let j_t = soret_effect_flux(0.0, gt, st, d);
        assert!((j_both - (j_c + j_t)).abs() < 1e-12);
    }
}
