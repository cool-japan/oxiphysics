#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// Parameters for the Cahn-Hilliard phase field model.
pub struct PhaseFieldParams {
    /// Cahn-Hilliard mobility M
    pub mobility: f64,
    /// Interface width W (epsilon)
    pub interface_width: f64,
    /// Surface tension sigma
    pub surface_tension: f64,
    /// Gradient energy coefficient kappa = sigma * W
    pub kappa: f64,
    /// Bulk free energy coefficient alpha
    pub alpha: f64,
    /// Number of components
    pub n_components: usize,
}
/// 2D Cahn-Hilliard phase field on a uniform grid.
///
/// The order parameter `phi` ranges from 0 (bulk phase B) to 1 (bulk phase A).
pub struct PhaseField {
    /// Grid size in x direction
    pub nx: usize,
    /// Grid size in y direction
    pub ny: usize,
    /// Phase field order parameter (flattened nx*ny, 0 to 1)
    pub phi: Vec<f64>,
    /// Chemical potential
    pub mu: Vec<f64>,
    /// Model parameters
    pub params: PhaseFieldParams,
}
impl PhaseField {
    /// Create a new PhaseField with all phi = 0.
    pub fn new(nx: usize, ny: usize, params: PhaseFieldParams) -> Self {
        let n = nx * ny;
        PhaseField {
            nx,
            ny,
            phi: vec![0.0; n],
            mu: vec![0.0; n],
            params,
        }
    }
    /// Flat index for grid cell (x, y).
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Initialize a circular droplet using the tanh profile.
    ///
    /// phi = 0.5 * (1 - tanh(2*(r - radius)/W))
    pub fn initialize_circle(&mut self, cx: f64, cy: f64, radius: f64) {
        let w = self.params.interface_width;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let r = (dx * dx + dy * dy).sqrt();
                let phi_val = 0.5 * (1.0 - ((2.0 * (r - radius)) / w).tanh());
                let k = self.idx(x, y);
                self.phi[k] = phi_val;
            }
        }
    }
    /// Compute the chemical potential using finite differences.
    ///
    /// mu = alpha * phi*(1-phi)*(1-2*phi) - kappa * laplacian(phi)
    ///
    /// Uses zero-flux (Neumann) boundary conditions.
    pub fn compute_chemical_potential(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let alpha = self.params.alpha;
        let kappa = self.params.kappa;
        let phi = self.phi.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let p = phi[k];
                let bulk = alpha * p * (1.0 - p) * (1.0 - 2.0 * p);
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let lap = phi[self.idx(xm, y)]
                    + phi[self.idx(xp, y)]
                    + phi[self.idx(x, ym)]
                    + phi[self.idx(x, yp)]
                    - 4.0 * p;
                self.mu[k] = bulk - kappa * lap;
            }
        }
    }
    /// Advance one Cahn-Hilliard time step (forward Euler).
    ///
    /// d(phi)/dt = M * laplacian(mu)
    pub fn step(&mut self, dt: f64) {
        self.compute_chemical_potential();
        let nx = self.nx;
        let ny = self.ny;
        let mobility = self.params.mobility;
        let mu = self.mu.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let lap_mu = mu[self.idx(xm, y)]
                    + mu[self.idx(xp, y)]
                    + mu[self.idx(x, ym)]
                    + mu[self.idx(x, yp)]
                    - 4.0 * mu[k];
                self.phi[k] += dt * mobility * lap_mu;
            }
        }
    }
    /// Compute the volume fraction (spatial average of phi).
    pub fn volume_fraction(&self) -> f64 {
        let n = self.phi.len();
        if n == 0 {
            return 0.0;
        }
        self.phi.iter().sum::<f64>() / n as f64
    }
    /// Compute the interface energy.
    ///
    /// E = 0.5 * kappa * sum(|grad phi|^2) * dx^2
    ///
    /// Uses dx = 1 (lattice units).
    pub fn interface_energy(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let kappa = self.params.kappa;
        let dx = 1.0_f64;
        let mut energy = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xp = if x < nx - 1 { x + 1 } else { x };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let dphi_dx = self.phi[self.idx(xp, y)] - self.phi[k];
                let dphi_dy = self.phi[self.idx(x, yp)] - self.phi[k];
                energy += dphi_dx * dphi_dx + dphi_dy * dphi_dy;
            }
        }
        0.5 * kappa * energy * dx * dx
    }
    /// Compute the total free energy (bulk + interface).
    ///
    /// ```text
    /// F = sum_cells [ alpha/4 * phi^2 * (1-phi)^2 + kappa/2 * |grad phi|^2 ]
    /// ```
    pub fn total_free_energy(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let alpha = self.params.alpha;
        let kappa = self.params.kappa;
        let mut energy = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let p = self.phi[k];
                let bulk = 0.25 * alpha * p * p * (1.0 - p) * (1.0 - p);
                let xp = if x < nx - 1 { x + 1 } else { x };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let dphi_dx = self.phi[self.idx(xp, y)] - p;
                let dphi_dy = self.phi[self.idx(x, yp)] - p;
                let grad_sq = dphi_dx * dphi_dx + dphi_dy * dphi_dy;
                energy += bulk + 0.5 * kappa * grad_sq;
            }
        }
        energy
    }
    /// Compute the gradient of phi at each cell using central differences.
    ///
    /// Returns a Vec of \[dphi_dx, dphi_dy\] for each cell.
    pub fn compute_gradient(&self) -> Vec<[f64; 2]> {
        let nx = self.nx;
        let ny = self.ny;
        let mut grad = vec![[0.0_f64; 2]; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let dx_span = (xp - xm).max(1) as f64;
                let dy_span = (yp - ym).max(1) as f64;
                grad[k][0] = (self.phi[self.idx(xp, y)] - self.phi[self.idx(xm, y)]) / dx_span;
                grad[k][1] = (self.phi[self.idx(x, yp)] - self.phi[self.idx(x, ym)]) / dy_span;
            }
        }
        grad
    }
    /// Compute the Laplacian of phi at each cell.
    pub fn compute_laplacian(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let mut lap = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                lap[k] = self.phi[self.idx(xm, y)]
                    + self.phi[self.idx(xp, y)]
                    + self.phi[self.idx(x, ym)]
                    + self.phi[self.idx(x, yp)]
                    - 4.0 * self.phi[k];
            }
        }
        lap
    }
    /// Compute the interface normal direction at each cell.
    ///
    /// n = grad(phi) / |grad(phi)| (set to zero where |grad(phi)| is too small).
    pub fn compute_interface_normal(&self) -> Vec<[f64; 2]> {
        let grad = self.compute_gradient();
        let mut normals = vec![[0.0_f64; 2]; self.nx * self.ny];
        for k in 0..grad.len() {
            let mag = (grad[k][0] * grad[k][0] + grad[k][1] * grad[k][1]).sqrt();
            if mag > 1e-12 {
                normals[k][0] = grad[k][0] / mag;
                normals[k][1] = grad[k][1] / mag;
            }
        }
        normals
    }
    /// Compute the interface curvature at each cell.
    ///
    /// kappa = -div(n) where n = grad(phi)/|grad(phi)|
    ///
    /// Returns the curvature at each cell.
    pub fn compute_curvature(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let normals = self.compute_interface_normal();
        let mut curvature = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let dx_span = (xp - xm).max(1) as f64;
                let dy_span = (yp - ym).max(1) as f64;
                let dn_x_dx = (normals[self.idx(xp, y)][0] - normals[self.idx(xm, y)][0]) / dx_span;
                let dn_y_dy = (normals[self.idx(x, yp)][1] - normals[self.idx(x, ym)][1]) / dy_span;
                curvature[k] = -(dn_x_dx + dn_y_dy);
            }
        }
        curvature
    }
    /// Compute surface tension force from the phase field (CSF approach).
    ///
    /// F_st = sigma * kappa * grad(phi) * delta(phi)
    ///
    /// where delta(phi) is approximated by |grad(phi)|.
    /// Returns force per cell as \[Fx, Fy\].
    pub fn surface_tension_force(&self) -> Vec<[f64; 2]> {
        let sigma = self.params.surface_tension;
        let curvature = self.compute_curvature();
        let grad = self.compute_gradient();
        let n = self.nx * self.ny;
        let mut force = vec![[0.0_f64; 2]; n];
        for k in 0..n {
            let mag = (grad[k][0] * grad[k][0] + grad[k][1] * grad[k][1]).sqrt();
            let scale = sigma * curvature[k] * mag;
            if mag > 1e-12 {
                force[k][0] = scale * grad[k][0] / mag;
                force[k][1] = scale * grad[k][1] / mag;
            }
        }
        force
    }
}
/// Parameters for the Swift-Orlandini-Yeomans (SOY) free-energy LBM model.
///
/// The SOY model (Swift, Orlandini, Osborn, Yeomans 1996) provides a
/// thermodynamically consistent multiphase LBM by incorporating a modified
/// equilibrium distribution that enforces the correct pressure tensor.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SoyModel {
    /// Landau parameter A (< 0 for phase separation).
    pub a: f64,
    /// Landau parameter B (> 0 for stability).
    pub b: f64,
    /// Gradient energy coefficient kappa.
    pub kappa: f64,
    /// Lattice sound speed squared cs2 = 1/3 in LBM units.
    pub cs2: f64,
}
impl SoyModel {
    /// Create a new SOY model.
    #[allow(dead_code)]
    pub fn new(a: f64, b: f64, kappa: f64) -> Self {
        Self {
            a,
            b,
            kappa,
            cs2: 1.0 / 3.0,
        }
    }
    /// Bulk chemical potential: `mu_bulk = A * phi + B * phi^3`.
    #[allow(dead_code)]
    pub fn mu_bulk(&self, phi: f64) -> f64 {
        self.a * phi + self.b * phi * phi * phi
    }
    /// Full chemical potential including gradient term:
    /// `mu = A * phi + B * phi^3 - kappa * lap(phi)`.
    #[allow(dead_code)]
    pub fn chemical_potential(&self, phi: f64, laplacian_phi: f64) -> f64 {
        self.mu_bulk(phi) - self.kappa * laplacian_phi
    }
    /// Isotropic bulk pressure: `p_bulk = cs^2 * rho + A/2 * phi^2 + 3B/4 * phi^4`.
    #[allow(dead_code)]
    pub fn bulk_pressure(&self, rho: f64, phi: f64) -> f64 {
        self.cs2 * rho + 0.5 * self.a * phi * phi + 0.75 * self.b * phi.powi(4)
    }
    /// Interface pressure correction:
    /// `P_interface = kappa * phi * lap(phi) + kappa/2 * |grad phi|^2`.
    #[allow(dead_code)]
    pub fn interface_pressure(&self, phi: f64, laplacian_phi: f64, grad_phi_sq: f64) -> f64 {
        self.kappa * phi * laplacian_phi + 0.5 * self.kappa * grad_phi_sq
    }
    /// Total thermodynamic pressure (bulk + interface corrections).
    #[allow(dead_code)]
    pub fn total_pressure(&self, rho: f64, phi: f64, laplacian_phi: f64, grad_phi_sq: f64) -> f64 {
        self.bulk_pressure(rho, phi) + self.interface_pressure(phi, laplacian_phi, grad_phi_sq)
    }
    /// SOY equilibrium distribution for direction `i` in D2Q9.
    ///
    /// `f_i^eq = w_i * [rho + rho * (c_i . u) / cs^2 + rho*(c_i.u)^2/(2cs^4)
    ///                   - rho * u^2/(2cs^2) + (P - rho*cs^2)/cs^2 * (c_i c_i - cs^2 I) ]`
    ///
    /// Here we expose the pressure-tensor correction term that distinguishes SOY
    /// from standard BGK.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn soy_equilibrium_correction(
        &self,
        w: f64,
        _rho: f64,
        phi: f64,
        laplacian_phi: f64,
        grad_phi_sq: f64,
        cx: f64,
        cy: f64,
    ) -> f64 {
        let p_ex = self.interface_pressure(phi, laplacian_phi, grad_phi_sq)
            + 0.5 * self.a * phi * phi
            + 0.75 * self.b * phi.powi(4);
        let c2 = cx * cx + cy * cy;
        let aniso = c2 - self.cs2;
        w * p_ex * aniso / (self.cs2 * self.cs2)
    }
    /// Equilibrium order-parameter distribution `g_i^eq` for the SOY phi equation.
    ///
    /// `g_i^eq = w_i * [phi + phi * (c_i . u) / cs^2 + Gamma * mu * (c_i c_i - cs^2 I)]`
    ///
    /// where `Gamma` is the mobility coefficient.
    #[allow(dead_code)]
    pub fn order_param_equilibrium(
        &self,
        w: f64,
        phi: f64,
        ux: f64,
        uy: f64,
        cx: f64,
        cy: f64,
        mu: f64,
        gamma: f64,
    ) -> f64 {
        let cu = cx * ux + cy * uy;
        let c2 = cx * cx + cy * cy;
        w * (phi + phi * cu / self.cs2 + gamma * mu * (c2 - self.cs2) / self.cs2)
    }
    /// Equilibrium droplet radius from the Young-Laplace balance.
    ///
    /// `R_eq = 2 * sigma / ΔP` where `sigma = sqrt(-8*kappa*A^3/(9*B^2))`.
    #[allow(dead_code)]
    pub fn equilibrium_radius(&self, delta_p: f64) -> f64 {
        if self.a >= 0.0 || self.b <= 0.0 || self.kappa <= 0.0 || delta_p < 1e-30 {
            return f64::INFINITY;
        }
        let sigma = (-8.0 * self.kappa * self.a.powi(3) / (9.0 * self.b * self.b)).sqrt();
        2.0 * sigma / delta_p
    }
    /// Coexistence order parameter values (Maxwell construction).
    ///
    /// Returns `(phi_minus, phi_plus)` with `phi_minus = -sqrt(-A/B)`.
    #[allow(dead_code)]
    pub fn coexistence_values(&self) -> (f64, f64) {
        if self.a < 0.0 && self.b > 0.0 {
            let phi_eq = (-self.a / self.b).sqrt();
            (-phi_eq, phi_eq)
        } else {
            (0.0, 0.0)
        }
    }
    /// Interface width parameter xi = sqrt(-kappa / A) (valid when A < 0).
    #[allow(dead_code)]
    pub fn interface_width_xi(&self) -> f64 {
        if self.a < 0.0 && self.kappa > 0.0 {
            (-self.kappa / self.a).sqrt()
        } else {
            0.0
        }
    }
    /// Surface tension: `sigma = sqrt(-8 * kappa * A^3 / (9 * B^2))`.
    #[allow(dead_code)]
    pub fn surface_tension(&self) -> f64 {
        if self.a < 0.0 && self.b > 0.0 && self.kappa > 0.0 {
            (-8.0 * self.kappa * self.a.powi(3) / (9.0 * self.b * self.b)).sqrt()
        } else {
            0.0
        }
    }
}
/// Coupled phase-field Lattice Boltzmann solver utilities.
///
/// Provides methods for computing interface geometric quantities and
/// applying contact-angle boundary conditions in diffuse-interface
/// (Cahn-Hilliard / Allen-Cahn) LBM simulations.
#[allow(dead_code)]
pub struct PhaseFieldLbm {
    /// Grid width (cells).
    pub nx: usize,
    /// Grid height (cells).
    pub ny: usize,
    /// Interface width parameter ε (lattice units).
    pub epsilon: f64,
    /// Surface tension σ (force per unit length).
    pub sigma: f64,
}
impl PhaseFieldLbm {
    /// Create a new `PhaseFieldLbm` helper.
    ///
    /// # Arguments
    /// - `nx`, `ny`  — grid dimensions
    /// - `epsilon`   — diffuse interface half-width ε (lattice units)
    /// - `sigma`     — surface tension coefficient σ
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, epsilon: f64, sigma: f64) -> Self {
        Self {
            nx,
            ny,
            epsilon,
            sigma,
        }
    }
    /// Compute the effective diffuse interface width from the phase-field profile.
    ///
    /// In the Cahn-Hilliard model the equilibrium profile is
    /// `φ(x) = 0.5 * (1 − tanh(x / (2ε)))`, so the transition from φ = 0.1
    /// to φ = 0.9 spans a distance `W = 4ε * atanh(0.8)`.
    ///
    /// For a given column `x_col`, this method estimates the interface width
    /// by scanning downward through the y-direction and measuring the distance
    /// over which φ changes from `phi_lo` to `phi_hi`.
    ///
    /// Returns the measured width in lattice units, or 0.0 if no interface
    /// is found.
    ///
    /// # Arguments
    /// - `phi`   — phase-field order parameter (flattened, row-major)
    /// - `x_col` — column index at which to measure
    /// - `phi_lo`, `phi_hi` — lower and upper threshold values (e.g. 0.1, 0.9)
    #[allow(dead_code)]
    pub fn compute_interface_width(
        &self,
        phi: &[f64],
        x_col: usize,
        phi_lo: f64,
        phi_hi: f64,
    ) -> f64 {
        let idx = |x: usize, y: usize| y * self.nx + x;
        let mut y_lo: Option<usize> = None;
        let mut y_hi: Option<usize> = None;
        for y in 0..self.ny {
            let p = phi[idx(x_col, y)];
            if y_lo.is_none() && p >= phi_lo {
                y_lo = Some(y);
            }
            if y_hi.is_none() && p >= phi_hi {
                y_hi = Some(y);
                break;
            }
        }
        match (y_lo, y_hi) {
            (Some(lo), Some(hi)) if hi > lo => (hi - lo) as f64,
            _ => 0.0,
        }
    }
    /// Apply a static contact-angle boundary condition at a solid wall.
    ///
    /// The wetting boundary condition (Jacqmin 2000) modifies the phase-field
    /// gradient normal to the wall to enforce a contact angle θ_c:
    ///
    /// `∂φ/∂n|_wall = −(σ_wA − σ_wB) / σ * f_wall(φ_wall)`
    ///
    /// where `f_wall(φ) = cos(θ_c)` for a simplified linear wetting model.
    ///
    /// This method sets the ghost-layer values for row `y = 0` (bottom wall)
    /// using a first-order approximation:
    ///
    /// `φ[x, -1] ← φ[x, 1] − 2 * cos(θ_contact) * |∇_eq|`
    ///
    /// where `|∇_eq| = (2/W) * φ_wall * (1 − φ_wall)` follows the
    /// equilibrium tanh profile.
    ///
    /// # Arguments
    /// - `phi`         — mutable phase-field array (nx * ny, row-major)
    /// - `theta_deg`   — contact angle in degrees \[0°, 180°\]
    /// - `wall_row`    — row index of the solid wall (usually 0)
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn apply_contact_angle(&self, phi: &mut [f64], theta_deg: f64, wall_row: usize) {
        let theta = theta_deg.to_radians();
        let cos_theta = theta.cos();
        let w = 4.0 * self.epsilon;
        let nx = self.nx;
        let ny = self.ny;
        for x in 0..nx {
            let k_wall = wall_row * nx + x;
            let phi_w = phi[k_wall].clamp(0.0, 1.0);
            let grad_eq = (2.0 / w) * phi_w * (1.0 - phi_w);
            let y_above = (wall_row + 1).min(ny - 1);
            let k_above = y_above * nx + x;
            let delta = cos_theta * grad_eq;
            phi[k_wall] = (phi[k_above] - delta).clamp(0.0, 1.0);
        }
    }
    /// Compute the capillary (Laplace) pressure jump across a curved interface.
    ///
    /// For a 2-D interface the Young-Laplace equation gives:
    ///
    /// `ΔP = σ * κ`
    ///
    /// where κ is the local curvature (positive for concave-upward droplet).
    /// The curvature is estimated from the phase-field via the divergence of
    /// the normalised gradient:
    ///
    /// `κ = −∇ · (∇φ / |∇φ|)`
    ///
    /// This method returns a `Vec`f64` of `σ * κ` at every grid cell.
    /// Cells with |∇φ| < `grad_thresh` are set to zero (bulk regions).
    ///
    /// # Arguments
    /// - `phi`        — phase-field order parameter (nx * ny)
    /// - `grad_thresh`— minimum gradient magnitude to compute curvature
    ///   (avoids division by zero in bulk; typically 1e-6)
    #[allow(dead_code)]
    pub fn compute_capillary_pressure(&self, phi: &[f64], grad_thresh: f64) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let idx = |x: usize, y: usize| y * nx + x;
        let mut dp = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let xp = (x + 1) % nx;
                let xm = (x + nx - 1) % nx;
                let yp = (y + 1) % ny;
                let ym = (y + ny - 1) % ny;
                let dphi_dx = 0.5 * (phi[idx(xp, y)] - phi[idx(xm, y)]);
                let dphi_dy = 0.5 * (phi[idx(x, yp)] - phi[idx(x, ym)]);
                let grad_mag = (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt();
                if grad_mag < grad_thresh {
                    continue;
                }
                let nx_n = dphi_dx / grad_mag;
                let ny_n = dphi_dy / grad_mag;
                let compute_normal_x = |xi: usize, yi: usize| -> f64 {
                    let xip = (xi + 1) % nx;
                    let xim = (xi + nx - 1) % nx;
                    let yip = (yi + 1) % ny;
                    let yim = (yi + ny - 1) % ny;
                    let gx = 0.5 * (phi[idx(xip, yi)] - phi[idx(xim, yi)]);
                    let gy = 0.5 * (phi[idx(xi, yip)] - phi[idx(xi, yim)]);
                    let gm = (gx * gx + gy * gy).sqrt();
                    if gm < grad_thresh { nx_n } else { gx / gm }
                };
                let compute_normal_y = |xi: usize, yi: usize| -> f64 {
                    let xip = (xi + 1) % nx;
                    let xim = (xi + nx - 1) % nx;
                    let yip = (yi + 1) % ny;
                    let yim = (yi + ny - 1) % ny;
                    let gx = 0.5 * (phi[idx(xip, yi)] - phi[idx(xim, yi)]);
                    let gy = 0.5 * (phi[idx(xi, yip)] - phi[idx(xi, yim)]);
                    let gm = (gx * gx + gy * gy).sqrt();
                    if gm < grad_thresh { ny_n } else { gy / gm }
                };
                let dn_x_dx = 0.5 * (compute_normal_x(xp, y) - compute_normal_x(xm, y));
                let dn_y_dy = 0.5 * (compute_normal_y(x, yp) - compute_normal_y(x, ym));
                let curvature = -(dn_x_dx + dn_y_dy);
                dp[idx(x, y)] = self.sigma * curvature;
            }
        }
        dp
    }
}
/// Allen-Cahn interface tracking model (simpler phase field).
pub struct AllenCahn {
    /// Grid size in x direction
    pub nx: usize,
    /// Grid size in y direction
    pub ny: usize,
    /// Phase field order parameter (flattened nx*ny)
    pub phi: Vec<f64>,
    /// Interface width parameter epsilon
    pub epsilon: f64,
    /// Mobility coefficient
    pub mobility: f64,
}
impl AllenCahn {
    /// Create a new AllenCahn field with all phi = 0.
    pub fn new(nx: usize, ny: usize, epsilon: f64, mobility: f64) -> Self {
        AllenCahn {
            nx,
            ny,
            phi: vec![0.0; nx * ny],
            epsilon,
            mobility,
        }
    }
    /// Flat index for grid cell (x, y).
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Initialize with a flat tanh interface at y = y_interface.
    ///
    /// phi(y) = 0.5 * (1 + tanh((y - y_interface) / (2*epsilon)))
    pub fn initialize_tanh_interface(&mut self, y_interface: f64, dx: f64) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let y_phys = y as f64 * dx;
                let phi_val = 0.5 * (1.0 + ((y_phys - y_interface) / (2.0 * self.epsilon)).tanh());
                let k = self.idx(x, y);
                self.phi[k] = phi_val;
            }
        }
    }
    /// Initialize a circular droplet using tanh profile.
    pub fn initialize_circle(&mut self, cx: f64, cy: f64, radius: f64) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let r = (dx * dx + dy * dy).sqrt();
                let phi_val = 0.5 * (1.0 - ((r - radius) / (2.0 * self.epsilon)).tanh());
                let k = self.idx(x, y);
                self.phi[k] = phi_val;
            }
        }
    }
    /// Advance one Allen-Cahn time step (forward Euler).
    ///
    /// d(phi)/dt = M * (epsilon^2 * laplacian(phi) - f'(phi))
    ///
    /// where f'(phi) = 2*phi*(1-phi)*(1-2*phi)  (double-well derivative)
    ///
    /// Uses zero-flux (Neumann) boundary conditions.
    pub fn step(&mut self, dt: f64, _dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let eps2 = self.epsilon * self.epsilon;
        let mobility = self.mobility;
        let phi_old = self.phi.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let p = phi_old[k];
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let lap = phi_old[self.idx(xm, y)]
                    + phi_old[self.idx(xp, y)]
                    + phi_old[self.idx(x, ym)]
                    + phi_old[self.idx(x, yp)]
                    - 4.0 * p;
                let df = 2.0 * p * (1.0 - p) * (1.0 - 2.0 * p);
                self.phi[k] += dt * mobility * (eps2 * lap - df);
            }
        }
    }
    /// Allen-Cahn step with advection by velocity field.
    ///
    /// d(phi)/dt + u . grad(phi) = M * (epsilon^2 * laplacian(phi) - f'(phi))
    pub fn step_with_advection(&mut self, dt: f64, ux: &[f64], uy: &[f64]) {
        let nx = self.nx;
        let ny = self.ny;
        let eps2 = self.epsilon * self.epsilon;
        let mobility = self.mobility;
        let phi_old = self.phi.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let p = phi_old[k];
                let xm = if x > 0 { x - 1 } else { x };
                let xp = if x < nx - 1 { x + 1 } else { x };
                let ym = if y > 0 { y - 1 } else { y };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let lap = phi_old[self.idx(xm, y)]
                    + phi_old[self.idx(xp, y)]
                    + phi_old[self.idx(x, ym)]
                    + phi_old[self.idx(x, yp)]
                    - 4.0 * p;
                let df = 2.0 * p * (1.0 - p) * (1.0 - 2.0 * p);
                let dphi_dx = if ux[k] > 0.0 {
                    p - phi_old[self.idx(xm, y)]
                } else {
                    phi_old[self.idx(xp, y)] - p
                };
                let dphi_dy = if uy[k] > 0.0 {
                    p - phi_old[self.idx(x, ym)]
                } else {
                    phi_old[self.idx(x, yp)] - p
                };
                let advection = ux[k] * dphi_dx + uy[k] * dphi_dy;
                self.phi[k] += dt * (mobility * (eps2 * lap - df) - advection);
            }
        }
    }
    /// Compute the sum of all phi values.
    pub fn total_phi(&self) -> f64 {
        self.phi.iter().sum()
    }
    /// Compute the interface area (length in 2D) by summing |grad phi|.
    pub fn interface_length(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut length = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                let xp = if x < nx - 1 { x + 1 } else { x };
                let yp = if y < ny - 1 { y + 1 } else { y };
                let dphi_dx = self.phi[self.idx(xp, y)] - self.phi[k];
                let dphi_dy = self.phi[self.idx(x, yp)] - self.phi[k];
                length += (dphi_dx * dphi_dx + dphi_dy * dphi_dy).sqrt();
            }
        }
        length
    }
}
