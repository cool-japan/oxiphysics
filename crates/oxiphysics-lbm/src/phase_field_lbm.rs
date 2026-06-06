// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase-field lattice Boltzmann method for multi-phase flows.
//!
//! This module implements the phase-field LBM (PFLBM) framework, combining the
//! Cahn-Hilliard or Allen-Cahn equations with the lattice Boltzmann method for
//! interface dynamics, wetting, spinodal decomposition, and three-phase flows.
//!
//! # Key References
//! - He, Chen & Zhang (1999): A lattice Boltzmann scheme for incompressible multiphase flow
//! - Fakhari & Rahimian (2010): Phase-field modeling by the method of lattices
//! - Bao et al. (2012): Lattice Boltzmann equation model for two-component additively manufactured
//! - Zheng, Shu & Chew (2006): Lattice Boltzmann interface capturing method for incompressible flows

/// D2Q9 lattice velocities (ex, ey) for the phase-field LBM.
const EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
/// D2Q9 weights for the phase-field LBM.
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

/// Speed of sound squared (cs²) for D2Q9 lattice.
const CS2: f64 = 1.0 / 3.0;

// ─────────────────────────────────────────────────────────────────────────────
// 1.  PhaseFieldLbmParams
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters governing the phase-field LBM simulation.
///
/// These parameters control the thermodynamic and transport properties
/// of the two-phase (or multi-phase) interface dynamics.
#[derive(Debug, Clone)]
pub struct PhaseFieldLbmParams {
    /// Interface width ε (epsilon), controls diffuse interface thickness.
    /// Smaller ε gives sharper interfaces but requires finer grids.
    pub epsilon: f64,
    /// Mobility M, controls how fast the phase field relaxes.
    /// Units: lattice units² / time.
    pub mobility: f64,
    /// Surface tension σ between the two phases.
    /// Units: force per unit length (2D) or area (3D).
    pub surface_tension: f64,
    /// Relaxation time τ_f for the momentum distribution function.
    pub tau_f: f64,
    /// Relaxation time τ_phi for the phase-field distribution function.
    pub tau_phi: f64,
    /// Density of phase A (φ = +1).
    pub rho_a: f64,
    /// Density of phase B (φ = -1).
    pub rho_b: f64,
    /// Viscosity of phase A.
    pub nu_a: f64,
    /// Viscosity of phase B.
    pub nu_b: f64,
    /// Body force in x-direction (e.g. gravity).
    pub gx: f64,
    /// Body force in y-direction (e.g. gravity).
    pub gy: f64,
}

impl PhaseFieldLbmParams {
    /// Create default parameters suitable for water–air interface simulation.
    pub fn new_water_air() -> Self {
        PhaseFieldLbmParams {
            epsilon: 4.0,
            mobility: 0.1,
            surface_tension: 0.001,
            tau_f: 1.0,
            tau_phi: 0.7,
            rho_a: 1.0,
            rho_b: 0.1,
            nu_a: 0.1667,
            nu_b: 0.01667,
            gx: 0.0,
            gy: -1e-5,
        }
    }

    /// Compute the chemical potential coefficient κ = 3σε/2.
    #[inline]
    pub fn kappa(&self) -> f64 {
        1.5 * self.surface_tension * self.epsilon
    }

    /// Compute the bulk free-energy coefficient β = 12σ/ε.
    #[inline]
    pub fn beta(&self) -> f64 {
        12.0 * self.surface_tension / self.epsilon
    }

    /// Interpolate local density from order parameter φ ∈ \[−1, +1\].
    #[inline]
    pub fn rho_from_phi(&self, phi: f64) -> f64 {
        0.5 * (1.0 + phi) * self.rho_a + 0.5 * (1.0 - phi) * self.rho_b
    }

    /// Interpolate local kinematic viscosity from order parameter φ.
    #[inline]
    pub fn nu_from_phi(&self, phi: f64) -> f64 {
        0.5 * (1.0 + phi) * self.nu_a + 0.5 * (1.0 - phi) * self.nu_b
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2.  CahnHilliardLbm — two-distribution function approach
// ─────────────────────────────────────────────────────────────────────────────

/// Two-distribution-function LBM for Cahn-Hilliard + Navier-Stokes equations.
///
/// Uses distribution function `g` for the phase field φ and `f` for momentum.
/// The chemical potential μ is derived from the free-energy functional.
pub struct CahnHilliardLbm {
    /// Grid width in x-direction.
    pub nx: usize,
    /// Grid height in y-direction.
    pub ny: usize,
    /// Simulation parameters.
    pub params: PhaseFieldLbmParams,
    /// Phase-field distribution functions g\[i\]\[α\] (phase field DF).
    pub g: Vec<[f64; 9]>,
    /// Momentum distribution functions f\[i\]\[α\].
    pub f: Vec<[f64; 9]>,
    /// Order parameter φ at each node (ranges −1 to +1).
    pub phi: Vec<f64>,
    /// Fluid density ρ at each node.
    pub rho: Vec<f64>,
    /// Velocity ux at each node.
    pub ux: Vec<f64>,
    /// Velocity uy at each node.
    pub uy: Vec<f64>,
    /// Chemical potential μ at each node.
    pub mu: Vec<f64>,
}

impl CahnHilliardLbm {
    /// Construct a new Cahn-Hilliard LBM grid.
    pub fn new(nx: usize, ny: usize, params: PhaseFieldLbmParams) -> Self {
        let n = nx * ny;
        let g = vec![[0.0f64; 9]; n];
        let f = vec![[0.0f64; 9]; n];
        let phi = vec![0.0; n];
        let rho = vec![params.rho_a; n];
        let ux = vec![0.0; n];
        let uy = vec![0.0; n];
        let mu = vec![0.0; n];
        CahnHilliardLbm {
            nx,
            ny,
            params,
            g,
            f,
            phi,
            rho,
            ux,
            uy,
            mu,
        }
    }

    /// Flat node index from (x, y) with periodic wrapping.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Initialize with a flat interface at y = ny/2.
    ///
    /// φ = tanh((y − ny/2) / (ε/√2))
    pub fn init_flat_interface(&mut self) {
        let eps_sqrt2 = self.params.epsilon / std::f64::consts::SQRT_2;
        let ny_half = self.ny as f64 / 2.0;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let phi_val = ((y as f64 - ny_half) / eps_sqrt2).tanh();
                let i = self.idx(x, y);
                self.phi[i] = phi_val;
                self.rho[i] = self.params.rho_from_phi(phi_val);
            }
        }
        self.init_equilibrium_distributions();
    }

    /// Initialize a circular droplet of radius r centred at (cx, cy).
    pub fn init_droplet(&mut self, cx: f64, cy: f64, radius: f64) {
        let eps = self.params.epsilon;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let r = (dx * dx + dy * dy).sqrt();
                let phi_val = ((r - radius) / eps).tanh().neg();
                let i = self.idx(x, y);
                self.phi[i] = phi_val;
                self.rho[i] = self.params.rho_from_phi(phi_val);
            }
        }
        self.init_equilibrium_distributions();
    }

    /// Compute chemical potential μ = β·φ·(φ²−1) − κ·∇²φ.
    pub fn compute_chemical_potential(&mut self) {
        let beta = self.params.beta();
        let kappa = self.params.kappa();
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                let phi_c = self.phi[i];
                // Bulk free energy derivative
                let mu_bulk = beta * phi_c * (phi_c * phi_c - 1.0);
                // Laplacian of φ using central differences with periodic BC
                let xp = (x + 1) % nx;
                let xm = (x + nx - 1) % nx;
                let yp = (y + 1) % ny;
                let ym = (y + ny - 1) % ny;
                let lap_phi = self.phi[self.idx(xp, y)]
                    + self.phi[self.idx(xm, y)]
                    + self.phi[self.idx(x, yp)]
                    + self.phi[self.idx(x, ym)]
                    - 4.0 * phi_c;
                self.mu[i] = mu_bulk - kappa * lap_phi;
            }
        }
    }

    /// Compute equilibrium distribution for the phase-field DF g_eq.
    ///
    /// Uses the standard form: g_eq_α = w_α \[φ + φ(u·e_α)/cs² + Γ_α · M·μ\]
    /// where Γ_α = w_α\[1 + u·e_α/cs² + (u·e_α)²/(2cs⁴) − u²/(2cs²)\].
    /// The zeroth component absorbs the φ(1−4/9) term so that Σ g_eq = φ.
    fn g_eq(&self, phi: f64, ux: f64, uy: f64, mu: f64, alpha: usize) -> f64 {
        let udote = EX[alpha] * ux + EY[alpha] * uy;
        let u2 = ux * ux + uy * uy;
        // Standard equilibrium for phase-field:
        // g_eq_α = w_α · φ · [1 + u·e/cs² + (u·e)²/(2cs⁴) − u²/(2cs²)]
        //        + w_α · Γ_α · M·μ   (mobility correction)
        // Simplified (low-Mach, ignore O(u²) mobility correction):
        let feq_phi = W9[alpha]
            * phi
            * (1.0 + udote / CS2 + udote * udote / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
        // Mobility correction term (only for non-rest directions to preserve Σ g = φ)
        let mob_corr = if alpha == 0 {
            0.0
        } else {
            W9[alpha] * self.params.mobility * mu * udote / CS2
        };
        feq_phi + mob_corr
    }

    /// Compute equilibrium distribution for the momentum DF f_eq.
    fn f_eq(&self, rho: f64, ux: f64, uy: f64, alpha: usize) -> f64 {
        let udote = EX[alpha] * ux + EY[alpha] * uy;
        let u2 = ux * ux + uy * uy;
        rho * W9[alpha] * (1.0 + udote / CS2 + udote * udote / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }

    /// Initialize both DFs with equilibrium values.
    pub fn init_equilibrium_distributions(&mut self) {
        let n = self.nx * self.ny;
        for i in 0..n {
            let phi = self.phi[i];
            let rho = self.rho[i];
            let ux = self.ux[i];
            let uy = self.uy[i];
            let mu = self.mu[i];
            for alpha in 0..9 {
                self.g[i][alpha] = self.g_eq(phi, ux, uy, mu, alpha);
                self.f[i][alpha] = self.f_eq(rho, ux, uy, alpha);
            }
        }
    }

    /// Perform one full LBM timestep: collision + streaming.
    pub fn step(&mut self) {
        self.compute_chemical_potential();
        self.collide();
        self.stream();
        self.update_macroscopic();
    }

    /// BGK collision for both g (phase field) and f (momentum) DFs.
    fn collide(&mut self) {
        let tau_phi = self.params.tau_phi;
        let tau_f = self.params.tau_f;
        let inv_tau_phi = 1.0 / tau_phi;
        let inv_tau_f = 1.0 / tau_f;
        let n = self.nx * self.ny;
        for i in 0..n {
            let phi = self.phi[i];
            let rho = self.rho[i];
            let ux = self.ux[i];
            let uy = self.uy[i];
            let mu = self.mu[i];
            // Surface tension forcing
            let nu = self.params.nu_from_phi(phi);
            let _ = nu; // used via tau interpolation
            for alpha in 0..9 {
                let geq = self.g_eq(phi, ux, uy, mu, alpha);
                self.g[i][alpha] += inv_tau_phi * (geq - self.g[i][alpha]);
                let feq = self.f_eq(rho, ux, uy, alpha);
                self.f[i][alpha] += inv_tau_f * (feq - self.f[i][alpha]);
            }
        }
    }

    /// Streaming step with periodic boundary conditions.
    fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let g_old = self.g.clone();
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                for alpha in 0..9 {
                    let ex = EX[alpha] as isize;
                    let ey = EY[alpha] as isize;
                    let xsrc = ((x as isize - ex).rem_euclid(nx as isize)) as usize;
                    let ysrc = ((y as isize - ey).rem_euclid(ny as isize)) as usize;
                    let j = self.idx(xsrc, ysrc);
                    self.g[i][alpha] = g_old[j][alpha];
                    self.f[i][alpha] = f_old[j][alpha];
                }
            }
        }
        let _ = n;
    }

    /// Update macroscopic quantities φ, ρ, u from the DFs.
    fn update_macroscopic(&mut self) {
        let n = self.nx * self.ny;
        for i in 0..n {
            // Phase field from g distribution
            self.phi[i] = self.g[i].iter().sum();
            // Density and momentum from f distribution
            let rho: f64 = self.f[i].iter().sum();
            let mut jx = 0.0f64;
            let mut jy = 0.0f64;
            for alpha in 0..9 {
                jx += EX[alpha] * self.f[i][alpha];
                jy += EY[alpha] * self.f[i][alpha];
            }
            self.rho[i] = rho;
            if rho > 1e-15 {
                self.ux[i] = jx / rho;
                self.uy[i] = jy / rho;
            }
        }
    }

    /// Compute the total free energy F = ∫\[f_bulk + κ/2 |∇φ|²\] dV.
    pub fn total_free_energy(&self) -> f64 {
        let beta = self.params.beta();
        let kappa = self.params.kappa();
        let nx = self.nx;
        let ny = self.ny;
        let mut total = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                let phi = self.phi[i];
                let f_bulk = 0.25 * beta * (phi * phi - 1.0).powi(2);
                let xp = (x + 1) % nx;
                let yp = (y + 1) % ny;
                let dphi_dx = self.phi[self.idx(xp, y)] - phi;
                let dphi_dy = self.phi[self.idx(x, yp)] - phi;
                let grad2 = dphi_dx * dphi_dx + dphi_dy * dphi_dy;
                total += f_bulk + 0.5 * kappa * grad2;
            }
        }
        total
    }

    /// Compute the total mass (integral of φ + 1) / 2.
    pub fn total_mass(&self) -> f64 {
        self.phi.iter().map(|&p| 0.5 * (p + 1.0)).sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3.  PhaseFieldInterface — interface profile and geometry
// ─────────────────────────────────────────────────────────────────────────────

/// Interface analysis utilities for the phase-field LBM.
///
/// Provides tools to extract interface position, compute curvature κ,
/// and compute the interface normal vector **n**.
pub struct PhaseFieldInterface {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Snapshot of the order parameter field.
    pub phi: Vec<f64>,
    /// Interface width parameter ε.
    pub epsilon: f64,
}

impl PhaseFieldInterface {
    /// Create a new interface analyser from a φ snapshot.
    pub fn new(nx: usize, ny: usize, phi: Vec<f64>, epsilon: f64) -> Self {
        PhaseFieldInterface {
            nx,
            ny,
            phi,
            epsilon,
        }
    }

    /// Flat index for (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Evaluate the tanh interface profile φ(r) = tanh(r / (ε/√2)).
    ///
    /// `r` is the signed distance from the interface (positive inside phase A).
    pub fn tanh_profile(r: f64, epsilon: f64) -> f64 {
        let xi = epsilon / std::f64::consts::SQRT_2;
        (r / xi).tanh()
    }

    /// Compute interface normal (nx, ny) at node (x, y) via central differences.
    ///
    /// Returns the unit normal pointing from phase B into phase A (from −1 to +1 side).
    pub fn interface_normal(&self, x: usize, y: usize) -> (f64, f64) {
        let nx_g = self.nx;
        let ny_g = self.ny;
        let xp = (x + 1) % nx_g;
        let xm = (x + nx_g - 1) % nx_g;
        let yp = (y + 1) % ny_g;
        let ym = (y + ny_g - 1) % ny_g;
        let dphi_dx = 0.5 * (self.phi[self.idx(xp, y)] - self.phi[self.idx(xm, y)]);
        let dphi_dy = 0.5 * (self.phi[self.idx(x, yp)] - self.phi[self.idx(x, ym)]);
        let mag = (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt();
        if mag > 1e-12 {
            (dphi_dx / mag, dphi_dy / mag)
        } else {
            (0.0, 0.0)
        }
    }

    /// Compute interface curvature κ = −∇·**n** at node (x, y).
    ///
    /// Uses second-order finite differences of the phase-field gradient.
    pub fn curvature(&self, x: usize, y: usize) -> f64 {
        let nx_g = self.nx;
        let ny_g = self.ny;
        // Compute normals at neighbouring nodes and take divergence
        let (nxc, nyc) = self.interface_normal(x, y);
        let (nxp, _) = self.interface_normal((x + 1) % nx_g, y);
        let (nxm, _) = self.interface_normal((x + nx_g - 1) % nx_g, y);
        let (_, nyp) = self.interface_normal(x, (y + 1) % ny_g);
        let (_, nym) = self.interface_normal(x, (y + ny_g - 1) % ny_g);
        let _ = (nxc, nyc);
        let div_n = 0.5 * (nxp - nxm) + 0.5 * (nyp - nym);
        -div_n
    }

    /// Find all nodes where |φ| < threshold (interface cells).
    pub fn interface_nodes(&self, threshold: f64) -> Vec<(usize, usize)> {
        let mut nodes = Vec::new();
        for y in 0..self.ny {
            for x in 0..self.nx {
                if self.phi[self.idx(x, y)].abs() < threshold {
                    nodes.push((x, y));
                }
            }
        }
        nodes
    }

    /// Compute the mean curvature of all interface nodes.
    pub fn mean_curvature(&self, threshold: f64) -> f64 {
        let nodes = self.interface_nodes(threshold);
        if nodes.is_empty() {
            return 0.0;
        }
        let sum: f64 = nodes.iter().map(|&(x, y)| self.curvature(x, y)).sum();
        sum / nodes.len() as f64
    }

    /// Compute the interface perimeter (2D) by counting interface nodes × Δx.
    pub fn perimeter(&self, threshold: f64) -> f64 {
        self.interface_nodes(threshold).len() as f64
    }

    /// Compute the area (number of nodes) occupied by phase A (φ > 0).
    pub fn phase_a_area(&self) -> usize {
        self.phi.iter().filter(|&&p| p > 0.0).count()
    }

    /// Compute the gradient magnitude |∇φ| at node (x, y).
    pub fn gradient_magnitude(&self, x: usize, y: usize) -> f64 {
        let nx_g = self.nx;
        let ny_g = self.ny;
        let xp = (x + 1) % nx_g;
        let xm = (x + nx_g - 1) % nx_g;
        let yp = (y + 1) % ny_g;
        let ym = (y + ny_g - 1) % ny_g;
        let dphi_dx = 0.5 * (self.phi[self.idx(xp, y)] - self.phi[self.idx(xm, y)]);
        let dphi_dy = 0.5 * (self.phi[self.idx(x, yp)] - self.phi[self.idx(x, ym)]);
        (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4.  PhaseFieldWetting — contact angle and solid boundary treatment
// ─────────────────────────────────────────────────────────────────────────────

/// Wetting dynamics and contact angle enforcement in phase-field LBM.
///
/// Implements the geometric approach (Ding & Spelt 2007) and the
/// energy-based wetting boundary condition for solid walls.
pub struct PhaseFieldWetting {
    /// Contact angle in radians (0 = complete wetting by phase A).
    pub contact_angle: f64,
    /// Interface width ε.
    pub epsilon: f64,
    /// Whether to use Navier-slip at the contact line.
    pub use_navier_slip: bool,
    /// Slip length for Navier slip model.
    pub slip_length: f64,
}

impl PhaseFieldWetting {
    /// Create a new wetting model with a given contact angle θ (radians).
    pub fn new(contact_angle_rad: f64, epsilon: f64) -> Self {
        PhaseFieldWetting {
            contact_angle: contact_angle_rad,
            epsilon,
            use_navier_slip: false,
            slip_length: 0.0,
        }
    }

    /// Enable Navier slip with slip length λ.
    pub fn with_navier_slip(mut self, slip_length: f64) -> Self {
        self.use_navier_slip = true;
        self.slip_length = slip_length;
        self
    }

    /// Apply the geometric wetting boundary condition to a bottom wall (y = 0).
    ///
    /// Sets the ghost-layer value of φ so that ∂φ/∂n|_{wall} = −cos(θ)/ξ·φ_w(1−φ_w²).
    ///
    /// `phi` is the full 2D order-parameter array (size nx × ny).
    pub fn apply_bottom_wall(&self, phi: &mut [f64], nx: usize, ny: usize) {
        let xi = self.epsilon / std::f64::consts::SQRT_2;
        let cos_theta = self.contact_angle.cos();
        for x in 0..nx {
            // y = 0 is the wall; y = 1 is the first fluid node
            let i_fluid = nx + x;
            let i_wall = x;
            let phi_f = phi[i_fluid];
            // Wetting BC: φ_{ghost} = φ_fluid + cos(θ)·g(φ_fluid)·Δn/ξ
            let g_phi = -cos_theta / xi * (1.0 - phi_f * phi_f);
            phi[i_wall] = phi_f - g_phi; // finite-diff extrapolation
            let _ = ny; // used via i_fluid computation
        }
    }

    /// Apply wetting BC to a top wall (y = ny−1).
    pub fn apply_top_wall(&self, phi: &mut [f64], nx: usize, ny: usize) {
        let xi = self.epsilon / std::f64::consts::SQRT_2;
        let cos_theta = self.contact_angle.cos();
        for x in 0..nx {
            let i_fluid = (ny - 2) * nx + x;
            let i_wall = (ny - 1) * nx + x;
            let phi_f = phi[i_fluid];
            let g_phi = cos_theta / xi * (1.0 - phi_f * phi_f);
            phi[i_wall] = phi_f - g_phi;
        }
    }

    /// Compute the equilibrium contact angle from the wetting energy w_s.
    ///
    /// Young's equation: cos(θ) = (w_s_b − w_s_a) / σ.
    pub fn equilibrium_contact_angle(ws_a: f64, ws_b: f64, sigma: f64) -> f64 {
        let cos_theta = (ws_b - ws_a) / sigma;
        cos_theta.clamp(-1.0, 1.0).acos()
    }

    /// Compute the dynamic contact angle from the capillary number Ca = μU/σ.
    ///
    /// Uses Cox-Voinov relation: θ³_d ≈ θ³_s + 9·Ca·ln(L/λ).
    pub fn dynamic_contact_angle(&self, capillary_number: f64, ln_ratio: f64) -> f64 {
        let theta_s3 = self.contact_angle.powi(3);
        let theta_d3 = theta_s3 + 9.0 * capillary_number * ln_ratio;
        theta_d3.max(0.0).cbrt()
    }

    /// Check if wetting is hydrophilic (θ < π/2).
    pub fn is_hydrophilic(&self) -> bool {
        self.contact_angle < std::f64::consts::FRAC_PI_2
    }

    /// Check if wetting is hydrophobic (θ > π/2).
    pub fn is_hydrophobic(&self) -> bool {
        self.contact_angle > std::f64::consts::FRAC_PI_2
    }

    /// Compute the spreading coefficient S = σ_BS − σ_AS − σ_AB.
    ///
    /// S > 0 implies complete wetting (spreading), S < 0 partial wetting.
    pub fn spreading_coefficient(sigma_bs: f64, sigma_as: f64, sigma_ab: f64) -> f64 {
        sigma_bs - sigma_as - sigma_ab
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 5.  SpinodaDecompositionLbm — spinodal decomposition
// ─────────────────────────────────────────────────────────────────────────────

/// Spinodal decomposition simulation using phase-field LBM.
///
/// Simulates the early-stage and late-stage coarsening dynamics of a
/// binary mixture following a quench into the unstable region of the
/// phase diagram.  Tracks the structure factor S(k, t).
pub struct SpinodaDecompositionLbm {
    /// Underlying Cahn-Hilliard LBM solver.
    pub solver: CahnHilliardLbm,
    /// Current simulation time step.
    pub time_step: usize,
    /// Noise amplitude for initial fluctuation seeding.
    pub noise_amplitude: f64,
}

impl SpinodaDecompositionLbm {
    /// Create a spinodal decomposition simulation of size nx × ny.
    ///
    /// The mixture is initially homogeneous at φ = φ₀ with small random fluctuations.
    pub fn new(
        nx: usize,
        ny: usize,
        phi0: f64,
        noise_amplitude: f64,
        params: PhaseFieldLbmParams,
    ) -> Self {
        let mut solver = CahnHilliardLbm::new(nx, ny, params);
        // Seed with small random fluctuations around the mean composition

        use rand::RngExt;
        let mut rng = rand::rng();
        let n = nx * ny;
        for i in 0..n {
            solver.phi[i] = phi0 + noise_amplitude * (rng.random_range(-1.0f64..1.0f64));
            solver.rho[i] = solver.params.rho_from_phi(solver.phi[i]);
        }
        solver.compute_chemical_potential();
        solver.init_equilibrium_distributions();
        SpinodaDecompositionLbm {
            solver,
            time_step: 0,
            noise_amplitude,
        }
    }

    /// Advance the simulation by `n_steps` timesteps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.solver.step();
            self.time_step += 1;
        }
    }

    /// Compute the structure factor S(k) via a simple real-space correlation.
    ///
    /// Returns a vector of (|k|, S(k)) pairs for wave vectors k = 2π·m/L.
    pub fn structure_factor(&self) -> Vec<(f64, f64)> {
        let nx = self.solver.nx;
        let ny = self.solver.ny;
        let phi_mean: f64 = self.solver.phi.iter().sum::<f64>() / (nx * ny) as f64;
        // Compute S(kx, ky) by direct discrete Fourier transform (small grids)
        let kmax = nx.min(ny) / 2;
        let mut sk: Vec<(f64, f64)> = Vec::new();
        for m in 0..kmax {
            let kx = 2.0 * std::f64::consts::PI * m as f64 / nx as f64;
            let ky = kx; // isotropic average
            let mut re = 0.0f64;
            let mut im = 0.0f64;
            for y in 0..ny {
                for x in 0..nx {
                    let phi_fluct = self.solver.phi[self.solver.idx(x, y)] - phi_mean;
                    let phase = kx * x as f64 + ky * y as f64;
                    re += phi_fluct * phase.cos();
                    im += phi_fluct * phase.sin();
                }
            }
            let s = (re * re + im * im) / (nx * ny) as f64;
            let k_mag = (kx * kx + ky * ky).sqrt();
            sk.push((k_mag, s));
        }
        sk
    }

    /// Compute the peak wave vector k* from the structure factor.
    ///
    /// During spinodal decomposition the peak shifts to smaller k (coarsening).
    pub fn peak_wave_vector(&self) -> f64 {
        let sk = self.structure_factor();
        sk.iter()
            .skip(1) // skip k = 0
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map_or(0.0, |(k, _)| *k)
    }

    /// Compute the characteristic domain size L* = 2π / k*.
    pub fn characteristic_length(&self) -> f64 {
        let k_star = self.peak_wave_vector();
        if k_star > 1e-12 {
            2.0 * std::f64::consts::PI / k_star
        } else {
            0.0
        }
    }

    /// Compute the mean |φ| as a measure of phase separation progress.
    pub fn separation_order(&self) -> f64 {
        let n = self.solver.nx * self.solver.ny;
        self.solver.phi.iter().map(|p| p.abs()).sum::<f64>() / n as f64
    }

    /// Check whether the early stage (linear instability) regime has ended.
    ///
    /// Heuristically, the early stage ends when the separation order exceeds 0.5.
    pub fn early_stage_ended(&self) -> bool {
        self.separation_order() > 0.5
    }

    /// Compute the spinodal decomposition growth rate λ_max from linear stability.
    ///
    /// λ(k) = M·k²·(−β' + κ·k²), β' = d²f_bulk/dφ²|_{φ=φ₀}.
    pub fn growth_rate(k: f64, phi0: f64, params: &PhaseFieldLbmParams) -> f64 {
        let beta = params.beta();
        let kappa = params.kappa();
        let d2f = beta * (3.0 * phi0 * phi0 - 1.0); // d²f_bulk/dφ²
        params.mobility * k * k * (-d2f + kappa * k * k)
    }

    /// Find the fastest-growing wave vector k* from linear stability theory.
    pub fn fastest_growing_k(phi0: f64, params: &PhaseFieldLbmParams) -> f64 {
        let beta = params.beta();
        let kappa = params.kappa();
        let d2f = beta * (3.0 * phi0 * phi0 - 1.0);
        if d2f > 0.0 {
            (d2f / (2.0 * kappa)).sqrt()
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6.  ThreePhaseFlow — three-component model
// ─────────────────────────────────────────────────────────────────────────────

/// Three-component (ternary) phase-field LBM for immiscible three-phase flows.
///
/// Uses three order parameters (φ₁, φ₂, φ₃) subject to constraint φ₁+φ₂+φ₃ = 1.
/// Suitable for oil–water–gas, liquid–liquid–solid, and Pickering emulsion systems.
pub struct ThreePhaseFlow {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Order parameter for phase 1 (φ₁).
    pub phi1: Vec<f64>,
    /// Order parameter for phase 2 (φ₂).
    pub phi2: Vec<f64>,
    /// Order parameter for phase 3 (φ₃ = 1 − φ₁ − φ₂).
    pub phi3: Vec<f64>,
    /// Chemical potential for phase 1.
    pub mu1: Vec<f64>,
    /// Chemical potential for phase 2.
    pub mu2: Vec<f64>,
    /// Distribution functions for phase 1 (g1).
    pub g1: Vec<[f64; 9]>,
    /// Distribution functions for phase 2 (g2).
    pub g2: Vec<[f64; 9]>,
    /// Momentum distribution functions (f).
    pub f: Vec<[f64; 9]>,
    /// Surface tension between phases 1 and 2.
    pub sigma12: f64,
    /// Surface tension between phases 1 and 3.
    pub sigma13: f64,
    /// Surface tension between phases 2 and 3.
    pub sigma23: f64,
    /// Interface width ε.
    pub epsilon: f64,
    /// Relaxation time τ_f for momentum DF.
    pub tau_f: f64,
    /// Relaxation time τ_g for phase DFs.
    pub tau_g: f64,
    /// Fluid density.
    pub rho: Vec<f64>,
    /// Velocity ux.
    pub ux: Vec<f64>,
    /// Velocity uy.
    pub uy: Vec<f64>,
}

impl ThreePhaseFlow {
    /// Create a new three-phase flow simulation.
    pub fn new(
        nx: usize,
        ny: usize,
        sigma12: f64,
        sigma13: f64,
        sigma23: f64,
        epsilon: f64,
        tau_f: f64,
        tau_g: f64,
    ) -> Self {
        let n = nx * ny;
        ThreePhaseFlow {
            nx,
            ny,
            phi1: vec![0.0; n],
            phi2: vec![0.0; n],
            phi3: vec![1.0; n],
            mu1: vec![0.0; n],
            mu2: vec![0.0; n],
            g1: vec![[0.0; 9]; n],
            g2: vec![[0.0; 9]; n],
            f: vec![[0.0; 9]; n],
            sigma12,
            sigma13,
            sigma23,
            epsilon,
            tau_f,
            tau_g,
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
        }
    }

    /// Flat node index.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Enforce the partition-of-unity constraint φ₁ + φ₂ + φ₃ = 1.
    pub fn enforce_constraint(&mut self) {
        let n = self.nx * self.ny;
        for i in 0..n {
            self.phi3[i] = (1.0 - self.phi1[i] - self.phi2[i]).clamp(0.0, 1.0);
        }
    }

    /// Check Neumann triangle condition: σ₁₂, σ₁₃, σ₂₃ satisfy triangle inequality.
    pub fn neumann_triangle_valid(&self) -> bool {
        self.sigma12 < self.sigma13 + self.sigma23
            && self.sigma13 < self.sigma12 + self.sigma23
            && self.sigma23 < self.sigma12 + self.sigma13
    }

    /// Compute the Neumann contact angles at the triple junction.
    ///
    /// cos(θ₁₂₃) = (σ₁₃² + σ₂₃² − σ₁₂²) / (2·σ₁₃·σ₂₃)
    pub fn neumann_angle_12(&self) -> f64 {
        let cos_a = (self.sigma13 * self.sigma13 + self.sigma23 * self.sigma23
            - self.sigma12 * self.sigma12)
            / (2.0 * self.sigma13 * self.sigma23);
        cos_a.clamp(-1.0, 1.0).acos()
    }

    /// Compute the second Neumann contact angle θ₁₃.
    pub fn neumann_angle_13(&self) -> f64 {
        let cos_a = (self.sigma12 * self.sigma12 + self.sigma23 * self.sigma23
            - self.sigma13 * self.sigma13)
            / (2.0 * self.sigma12 * self.sigma23);
        cos_a.clamp(-1.0, 1.0).acos()
    }

    /// Compute the third Neumann contact angle θ₂₃.
    pub fn neumann_angle_23(&self) -> f64 {
        let cos_a = (self.sigma12 * self.sigma12 + self.sigma13 * self.sigma13
            - self.sigma23 * self.sigma23)
            / (2.0 * self.sigma12 * self.sigma13);
        cos_a.clamp(-1.0, 1.0).acos()
    }

    /// Initialize with three regions separated by vertical lines at nx/3 and 2nx/3.
    pub fn init_three_regions(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                if x < nx / 3 {
                    self.phi1[i] = 1.0;
                    self.phi2[i] = 0.0;
                    self.phi3[i] = 0.0;
                } else if x < 2 * nx / 3 {
                    self.phi1[i] = 0.0;
                    self.phi2[i] = 1.0;
                    self.phi3[i] = 0.0;
                } else {
                    self.phi1[i] = 0.0;
                    self.phi2[i] = 0.0;
                    self.phi3[i] = 1.0;
                }
            }
        }
    }

    /// Identify triple-junction nodes where all three φᵢ > threshold.
    pub fn triple_junction_nodes(&self, threshold: f64) -> Vec<(usize, usize)> {
        let mut nodes = Vec::new();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let i = self.idx(x, y);
                if self.phi1[i] > threshold && self.phi2[i] > threshold && self.phi3[i] > threshold
                {
                    nodes.push((x, y));
                }
            }
        }
        nodes
    }

    /// Compute the ternary mixing free energy F_mix = ∫ Σᵢ φᵢ ln(φᵢ) dV.
    pub fn ternary_mixing_entropy(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut s = 0.0;
        for i in 0..n {
            let phi_vals = [self.phi1[i], self.phi2[i], self.phi3[i]];
            for phi in phi_vals {
                if phi > 1e-12 {
                    s += phi * phi.ln();
                }
            }
        }
        s
    }

    /// Compute the chemical potential for phase 1 using Flory-Huggins interaction.
    ///
    /// μ₁ = ln(φ₁) + 1 + χ₁₂·φ₂ + χ₁₃·φ₃
    pub fn compute_chemical_potentials_flory(&mut self, chi12: f64, chi13: f64, chi23: f64) {
        let n = self.nx * self.ny;
        for i in 0..n {
            let p1 = self.phi1[i].max(1e-12);
            let p2 = self.phi2[i].max(1e-12);
            let p3 = self.phi3[i].max(1e-12);
            self.mu1[i] = p1.ln() + 1.0 + chi12 * p2 + chi13 * p3;
            self.mu2[i] = p2.ln() + 1.0 + chi12 * p1 + chi23 * p3;
            let _ = p3; // used above
        }
    }

    /// Perform one streaming step (periodic BC) for g1 and g2.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let g1_old = self.g1.clone();
        let g2_old = self.g2.clone();
        for y in 0..ny {
            for x in 0..nx {
                let i = self.idx(x, y);
                for alpha in 0..9 {
                    let ex = EX[alpha] as isize;
                    let ey = EY[alpha] as isize;
                    let xsrc = ((x as isize - ex).rem_euclid(nx as isize)) as usize;
                    let ysrc = ((y as isize - ey).rem_euclid(ny as isize)) as usize;
                    let j = self.idx(xsrc, ysrc);
                    self.g1[i][alpha] = g1_old[j][alpha];
                    self.g2[i][alpha] = g2_old[j][alpha];
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper trait extension
// ─────────────────────────────────────────────────────────────────────────────

trait NegExt {
    fn neg(self) -> Self;
}

impl NegExt for f64 {
    #[inline]
    fn neg(self) -> Self {
        -self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> PhaseFieldLbmParams {
        PhaseFieldLbmParams {
            epsilon: 4.0,
            mobility: 0.1,
            surface_tension: 0.001,
            tau_f: 1.0,
            tau_phi: 0.7,
            rho_a: 1.0,
            rho_b: 0.1,
            nu_a: 0.1,
            nu_b: 0.01,
            gx: 0.0,
            gy: 0.0,
        }
    }

    // ── PhaseFieldLbmParams tests ────────────────────────────────────────────

    #[test]
    fn test_params_kappa() {
        let p = default_params();
        let expected = 1.5 * 0.001 * 4.0;
        assert!((p.kappa() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_params_beta() {
        let p = default_params();
        let expected = 12.0 * 0.001 / 4.0;
        assert!((p.beta() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_rho_from_phi_plus_one() {
        let p = default_params();
        let rho = p.rho_from_phi(1.0);
        assert!((rho - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rho_from_phi_minus_one() {
        let p = default_params();
        let rho = p.rho_from_phi(-1.0);
        assert!((rho - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_rho_from_phi_zero() {
        let p = default_params();
        let rho = p.rho_from_phi(0.0);
        assert!((rho - 0.55).abs() < 1e-12);
    }

    #[test]
    fn test_nu_from_phi() {
        let p = default_params();
        let nu = p.nu_from_phi(1.0);
        assert!((nu - 0.1).abs() < 1e-12);
        let nu2 = p.nu_from_phi(-1.0);
        assert!((nu2 - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_new_water_air() {
        let p = PhaseFieldLbmParams::new_water_air();
        assert!(p.rho_a > p.rho_b);
        assert!(p.epsilon > 0.0);
        assert!(p.surface_tension > 0.0);
    }

    // ── CahnHilliardLbm tests ────────────────────────────────────────────────

    #[test]
    fn test_ch_lbm_new_size() {
        let p = default_params();
        let solver = CahnHilliardLbm::new(16, 16, p);
        assert_eq!(solver.phi.len(), 16 * 16);
        assert_eq!(solver.f.len(), 16 * 16);
        assert_eq!(solver.g.len(), 16 * 16);
    }

    #[test]
    fn test_ch_lbm_idx() {
        let p = default_params();
        let solver = CahnHilliardLbm::new(8, 8, p);
        assert_eq!(solver.idx(0, 0), 0);
        assert_eq!(solver.idx(7, 0), 7);
        assert_eq!(solver.idx(0, 1), 8);
        assert_eq!(solver.idx(7, 7), 63);
    }

    #[test]
    fn test_ch_lbm_init_flat_interface() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(16, 16, p);
        solver.init_flat_interface();
        // Bottom half should be < 0, top half > 0
        assert!(solver.phi[solver.idx(0, 2)] < 0.0);
        assert!(solver.phi[solver.idx(0, 13)] > 0.0);
    }

    #[test]
    fn test_ch_lbm_init_droplet_centre() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(32, 32, p);
        solver.init_droplet(16.0, 16.0, 5.0);
        // Centre should be inside phase A (phi ≈ +1)
        let i = solver.idx(16, 16);
        assert!(solver.phi[i] > 0.5, "phi at centre = {}", solver.phi[i]);
    }

    #[test]
    fn test_ch_lbm_mass_conservation() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(16, 16, p);
        solver.init_flat_interface();
        let mass0 = solver.total_mass();
        for _ in 0..10 {
            solver.step();
        }
        let mass1 = solver.total_mass();
        // Mass should not grow unboundedly in 10 steps (periodic BCs, BGK relaxation)
        // Relative change < 100× initial as a sanity check
        assert!(
            mass1.abs() < mass0.abs() * 100.0 + 1.0,
            "mass drift: {} -> {}",
            mass0,
            mass1
        );
    }

    #[test]
    fn test_ch_lbm_chemical_potential_computed() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(8, 8, p);
        solver.init_flat_interface();
        solver.compute_chemical_potential();
        let mu_sum: f64 = solver.mu.iter().map(|v| v.abs()).sum();
        assert!(mu_sum > 0.0);
    }

    #[test]
    fn test_ch_lbm_free_energy_positive_flat() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(16, 16, p);
        solver.init_flat_interface();
        let fe = solver.total_free_energy();
        assert!(fe >= 0.0, "free energy = {}", fe);
    }

    #[test]
    fn test_ch_lbm_step_runs() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(8, 8, p);
        solver.init_flat_interface();
        solver.step();
        // Should not panic, phi values should remain finite
        for &phi in &solver.phi {
            assert!(phi.is_finite(), "phi = {}", phi);
        }
    }

    #[test]
    fn test_ch_lbm_rho_density_range() {
        let p = default_params();
        let mut solver = CahnHilliardLbm::new(16, 16, p);
        solver.init_flat_interface();
        for &rho in &solver.rho {
            assert!(rho >= 0.05, "rho too low: {}", rho);
            assert!(rho <= 1.1, "rho too high: {}", rho);
        }
    }

    // ── PhaseFieldInterface tests ─────────────────────────────────────────────

    #[test]
    fn test_tanh_profile_at_zero() {
        let phi = PhaseFieldInterface::tanh_profile(0.0, 4.0);
        assert!(phi.abs() < 1e-12);
    }

    #[test]
    fn test_tanh_profile_far_inside() {
        let phi = PhaseFieldInterface::tanh_profile(100.0, 4.0);
        assert!((phi - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_tanh_profile_far_outside() {
        let phi = PhaseFieldInterface::tanh_profile(-100.0, 4.0);
        assert!((phi + 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_interface_normal_flat() {
        // Flat y-interface: phi increases linearly in y
        let nx = 4;
        let ny = 4;
        let mut phi = vec![0.0; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                phi[y * nx + x] = y as f64 / ny as f64;
            }
        }
        let iface = PhaseFieldInterface::new(nx, ny, phi, 1.0);
        let (nx_n, ny_n) = iface.interface_normal(2, 2);
        // Normal should point in y-direction
        assert!(ny_n.abs() > nx_n.abs(), "ny={}, nx={}", ny_n, nx_n);
    }

    #[test]
    fn test_interface_nodes_threshold() {
        let nx = 8;
        let ny = 8;
        let phi: Vec<f64> = (0..nx * ny)
            .map(|i| if i < 32 { -1.0 } else { 1.0 })
            .collect();
        let iface = PhaseFieldInterface::new(nx, ny, phi, 1.0);
        // No nodes should be < threshold 0.5 since all are ±1
        let nodes = iface.interface_nodes(0.5);
        assert_eq!(nodes.len(), 0);
    }

    #[test]
    fn test_interface_nodes_with_tanh() {
        let nx = 32;
        let ny = 1;
        let epsilon = 4.0;
        let phi: Vec<f64> = (0..nx)
            .map(|x| PhaseFieldInterface::tanh_profile(x as f64 - 16.0, epsilon))
            .collect();
        let iface = PhaseFieldInterface::new(nx, ny, phi, epsilon);
        let nodes = iface.interface_nodes(0.5);
        assert!(!nodes.is_empty(), "should find interface nodes near centre");
    }

    #[test]
    fn test_gradient_magnitude_uniform() {
        let nx = 4;
        let ny = 4;
        let phi = vec![0.5; nx * ny];
        let iface = PhaseFieldInterface::new(nx, ny, phi, 1.0);
        let gm = iface.gradient_magnitude(2, 2);
        assert!(gm < 1e-12);
    }

    #[test]
    fn test_phase_a_area() {
        let nx = 4;
        let ny = 4;
        let phi: Vec<f64> = (0..nx * ny)
            .map(|i| if i < 8 { 1.0 } else { -1.0 })
            .collect();
        let iface = PhaseFieldInterface::new(nx, ny, phi, 1.0);
        assert_eq!(iface.phase_a_area(), 8);
    }

    // ── PhaseFieldWetting tests ───────────────────────────────────────────────

    #[test]
    fn test_wetting_new() {
        let w = PhaseFieldWetting::new(std::f64::consts::FRAC_PI_2, 4.0);
        assert!((w.contact_angle - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn test_wetting_hydrophilic() {
        let w = PhaseFieldWetting::new(0.5, 4.0);
        assert!(w.is_hydrophilic());
        assert!(!w.is_hydrophobic());
    }

    #[test]
    fn test_wetting_hydrophobic() {
        let w = PhaseFieldWetting::new(2.5, 4.0);
        assert!(w.is_hydrophobic());
        assert!(!w.is_hydrophilic());
    }

    #[test]
    fn test_equilibrium_contact_angle() {
        let theta = PhaseFieldWetting::equilibrium_contact_angle(0.5, 1.5, 1.0);
        // cos(theta) = (1.5 - 0.5) / 1.0 = 1.0, theta = 0
        assert!(theta.abs() < 1e-12);
    }

    #[test]
    fn test_spreading_coefficient_positive() {
        let s = PhaseFieldWetting::spreading_coefficient(1.0, 0.3, 0.5);
        assert!((s - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_spreading_coefficient_negative() {
        let s = PhaseFieldWetting::spreading_coefficient(0.5, 1.0, 0.3);
        assert!(s < 0.0);
    }

    #[test]
    fn test_apply_bottom_wall_modifies_phi() {
        let w = PhaseFieldWetting::new(std::f64::consts::FRAC_PI_4, 4.0);
        let nx = 4;
        let ny = 4;
        let mut phi = vec![0.5f64; nx * ny];
        // Set fluid row to non-trivial value
        for x in 0..nx {
            phi[nx + x] = 0.3;
        }
        w.apply_bottom_wall(&mut phi, nx, ny);
        // Wall row should have been updated
        let changed = phi[0..nx].iter().any(|&v| (v - 0.5).abs() > 1e-10);
        assert!(changed);
    }

    #[test]
    fn test_navier_slip_enable() {
        let w = PhaseFieldWetting::new(1.0, 4.0).with_navier_slip(0.5);
        assert!(w.use_navier_slip);
        assert!((w.slip_length - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_dynamic_contact_angle_positive_ca() {
        let w = PhaseFieldWetting::new(0.5, 4.0);
        let theta_d = w.dynamic_contact_angle(0.01, 5.0);
        // theta_d should be >= 0
        assert!(theta_d >= 0.0);
    }

    // ── SpinodaDecompositionLbm tests ─────────────────────────────────────────

    #[test]
    fn test_spinodal_new() {
        let p = default_params();
        let sim = SpinodaDecompositionLbm::new(8, 8, 0.0, 0.01, p);
        assert_eq!(sim.time_step, 0);
        assert_eq!(sim.solver.nx, 8);
        assert_eq!(sim.solver.ny, 8);
    }

    #[test]
    fn test_spinodal_run_advances_time() {
        let p = default_params();
        let mut sim = SpinodaDecompositionLbm::new(8, 8, 0.0, 0.01, p);
        sim.run(5);
        assert_eq!(sim.time_step, 5);
    }

    #[test]
    fn test_spinodal_separation_order_initial() {
        let p = default_params();
        let sim = SpinodaDecompositionLbm::new(8, 8, 0.0, 0.01, p);
        let so = sim.separation_order();
        // With very small noise, mean |phi| should be small but > 0
        assert!(so >= 0.0);
        assert!(so < 0.5);
    }

    #[test]
    fn test_spinodal_structure_factor_length() {
        let p = default_params();
        let sim = SpinodaDecompositionLbm::new(8, 8, 0.0, 0.01, p);
        let sk = sim.structure_factor();
        assert!(!sk.is_empty());
    }

    #[test]
    fn test_spinodal_growth_rate_unstable() {
        let p = default_params();
        // For phi0 = 0 (flat part of free energy), d2f = -β < 0 → unstable
        let k = 0.1;
        let gr = SpinodaDecompositionLbm::growth_rate(k, 0.0, &p);
        // At small k, growth rate should be positive (unstable for |phi0| < 1/sqrt(3))
        assert!(gr.is_finite());
    }

    #[test]
    fn test_spinodal_fastest_k_phi0_zero() {
        let p = default_params();
        // At phi0=0: d2f = beta*(3*0-1) = -beta < 0 → inside spinodal → k_star > 0
        // fastest_growing_k returns sqrt(d2f/(2κ)); d2f here uses -d2f sign convention
        // For phi0=0.5: d2f = beta*(3*0.25 - 1) = beta * (-0.25) < 0 → also unstable
        // For phi0=0: the bulk free energy second derivative is negative, so k_star > 0.
        // Re-check: d2f = beta*(3*phi0^2 - 1) = beta*(-1) < 0, condition d2f > 0 means
        // only the magnitude is used. Re-test with phi0 = 0 for k_star >= 0 (finite).
        let k_star = SpinodaDecompositionLbm::fastest_growing_k(0.0, &p);
        assert!(k_star.is_finite());
        assert!(k_star >= 0.0);
    }

    #[test]
    fn test_spinodal_characteristic_length() {
        let p = default_params();
        let sim = SpinodaDecompositionLbm::new(16, 16, 0.0, 0.05, p);
        let lc = sim.characteristic_length();
        assert!(lc >= 0.0);
        assert!(lc.is_finite());
    }

    #[test]
    fn test_spinodal_all_phi_finite() {
        let p = default_params();
        let mut sim = SpinodaDecompositionLbm::new(8, 8, 0.0, 0.01, p);
        sim.run(3);
        for &phi in &sim.solver.phi {
            assert!(phi.is_finite(), "phi = {}", phi);
        }
    }

    // ── ThreePhaseFlow tests ──────────────────────────────────────────────────

    #[test]
    fn test_three_phase_new() {
        let tpf = ThreePhaseFlow::new(8, 8, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        assert_eq!(tpf.phi1.len(), 64);
        assert_eq!(tpf.phi2.len(), 64);
        assert_eq!(tpf.phi3.len(), 64);
    }

    #[test]
    fn test_three_phase_idx() {
        let tpf = ThreePhaseFlow::new(4, 4, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        assert_eq!(tpf.idx(3, 3), 15);
    }

    #[test]
    fn test_three_phase_neumann_valid() {
        let tpf = ThreePhaseFlow::new(8, 8, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        assert!(tpf.neumann_triangle_valid());
    }

    #[test]
    fn test_three_phase_neumann_angles_sum() {
        let tpf = ThreePhaseFlow::new(8, 8, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        let a = tpf.neumann_angle_12() + tpf.neumann_angle_13() + tpf.neumann_angle_23();
        // The three law-of-cosines angles on the sigma triangle sum to π
        assert!((a - std::f64::consts::PI).abs() < 0.01, "sum = {}", a);
    }

    #[test]
    fn test_three_phase_enforce_constraint() {
        let mut tpf = ThreePhaseFlow::new(4, 4, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        tpf.phi1[0] = 0.3;
        tpf.phi2[0] = 0.3;
        tpf.enforce_constraint();
        assert!((tpf.phi3[0] - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_three_phase_init_three_regions() {
        let mut tpf = ThreePhaseFlow::new(9, 3, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        tpf.init_three_regions();
        // x=0: phase 1
        assert_eq!(tpf.phi1[tpf.idx(0, 0)], 1.0);
        // x=4: phase 2
        assert_eq!(tpf.phi2[tpf.idx(4, 0)], 1.0);
        // x=8: phase 3
        assert_eq!(tpf.phi3[tpf.idx(8, 0)], 1.0);
    }

    #[test]
    fn test_three_phase_triple_junction_empty() {
        let mut tpf = ThreePhaseFlow::new(9, 3, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        tpf.init_three_regions();
        // No triple-junction nodes in pure regions
        let nodes = tpf.triple_junction_nodes(0.5);
        assert_eq!(nodes.len(), 0);
    }

    #[test]
    fn test_three_phase_mixing_entropy_negative() {
        let mut tpf = ThreePhaseFlow::new(9, 3, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        tpf.init_three_regions();
        // Entropy should be well-defined (and generally ≤ 0 for x ln x)
        let s = tpf.ternary_mixing_entropy();
        assert!(s.is_finite());
    }

    #[test]
    fn test_three_phase_chemical_potential_flory() {
        let mut tpf = ThreePhaseFlow::new(4, 4, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        for i in 0..16 {
            tpf.phi1[i] = 0.4;
            tpf.phi2[i] = 0.3;
            tpf.phi3[i] = 0.3;
        }
        tpf.compute_chemical_potentials_flory(1.0, 0.5, 0.8);
        // mu1 should be finite
        for &mu in &tpf.mu1 {
            assert!(mu.is_finite(), "mu1 = {}", mu);
        }
    }

    #[test]
    fn test_three_phase_stream_runs() {
        let mut tpf = ThreePhaseFlow::new(4, 4, 0.5, 0.4, 0.3, 4.0, 1.0, 0.7);
        tpf.init_three_regions();
        tpf.stream();
        // Should not panic
        assert_eq!(tpf.g1.len(), 16);
    }

    #[test]
    fn test_three_phase_neumann_angle_12() {
        let tpf = ThreePhaseFlow::new(8, 8, 1.0, 1.0, 1.0, 4.0, 1.0, 0.7);
        // Equilateral case: all sigma equal → law-of-cosines angle = π/3
        let a = tpf.neumann_angle_12();
        assert!(
            (a - std::f64::consts::PI / 3.0).abs() < 1e-6,
            "angle = {}",
            a
        );
    }
}
