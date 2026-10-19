//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::grid::LbmGrid2D;
use crate::lattice::{CS2, D2Q9_VELOCITIES, D2Q9_WEIGHTS};

/// Pseudo-potential function type for the Shan-Chen model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PsiType {
    /// Linear: `psi(rho) = rho`.
    Linear,
    /// Exponential (Shan-Chen original): `psi(rho) = rho_0 * (1 - exp(-rho/rho_0))`.
    Exponential,
    /// Sukop-Thorne: `psi(rho) = psi_0 * exp(-rho_0 / rho)`.
    SukopThorne,
    /// Yuan-Schaefer: `psi(rho) = sqrt(2 * (p_eos - rho * cs^2) / (G * cs^2))`.
    YuanSchaefer,
}
/// Free-energy multiphase LBM model (Swift et al.).
///
/// Uses a Landau-Ginzburg free energy functional of the form:
/// `F = integral [ A/2 * phi^2 + B/4 * phi^4 - kappa/2 * (grad phi)^2 ] dx`
///
/// The pressure tensor includes a thermodynamic correction that enforces
/// mechanical equilibrium across the interface.
#[derive(Debug, Clone)]
pub struct FreeEnergyModel {
    /// Surface tension parameter kappa (controls interface width).
    pub kappa: f64,
    /// Landau parameter A (negative for phase separation).
    pub a: f64,
    /// Landau parameter B (positive for stability).
    pub b: f64,
}
impl FreeEnergyModel {
    /// Create a new free-energy model.
    pub fn new(kappa: f64, a: f64, b: f64) -> Self {
        Self { kappa, a, b }
    }
    /// Compute the thermodynamic pressure including gradient corrections.
    ///
    /// `p = cs² * rho + A/2 * phi² + 3B/4 * phi⁴ - kappa * phi * laplacian(phi) - kappa/2 * (grad phi)²`
    ///
    /// In the bulk approximation (no gradients), this reduces to:
    /// `p_bulk = cs² * rho + A/2 * rho² + B/4 * rho⁴`
    ///
    /// The `grad_rho_sq` argument is `|grad rho|²`.
    pub fn compute_free_energy_pressure(&self, rho: f64, grad_rho_sq: f64) -> f64 {
        use crate::lattice::CS2;
        let p_bulk = CS2 * rho + 0.5 * self.a * rho * rho + 0.25 * self.b * rho * rho * rho * rho;
        let p_grad = -0.5 * self.kappa * grad_rho_sq;
        p_bulk + p_grad
    }
    /// Equilibrium density for the liquid phase (rho > 0, phi > 0).
    ///
    /// At equilibrium: `A * phi + B * phi^3 = 0` → `phi = sqrt(-A/B)`.
    pub fn rho_liquid(&self) -> f64 {
        if self.a < 0.0 && self.b > 0.0 {
            (-self.a / self.b).sqrt()
        } else {
            1.0
        }
    }
    /// Equilibrium density for the vapour phase (symmetric: `rho_vapor = -rho_liquid`).
    pub fn rho_vapor(&self) -> f64 {
        -self.rho_liquid()
    }
    /// Interface width parameter: `xi = sqrt(-2 kappa / A)` (valid when A < 0).
    pub fn interface_width(&self) -> f64 {
        if self.a < 0.0 && self.kappa > 0.0 {
            (-2.0 * self.kappa / self.a).sqrt()
        } else {
            0.0
        }
    }
    /// Surface tension: `sigma = sqrt(-8 kappa A^3 / (9 B^2))`.
    pub fn surface_tension(&self) -> f64 {
        if self.a < 0.0 && self.b > 0.0 && self.kappa > 0.0 {
            (-8.0 * self.kappa * self.a.powi(3) / (9.0 * self.b * self.b)).sqrt()
        } else {
            0.0
        }
    }
    /// Compute the bulk free energy density: `f(phi) = A/2 * phi^2 + B/4 * phi^4`.
    pub fn bulk_free_energy_density(&self, phi: f64) -> f64 {
        0.5 * self.a * phi * phi + 0.25 * self.b * phi * phi * phi * phi
    }
    /// Compute the chemical potential: `mu = A * phi + B * phi^3 - kappa * lap(phi)`.
    pub fn chemical_potential(&self, phi: f64, laplacian_phi: f64) -> f64 {
        self.a * phi + self.b * phi * phi * phi - self.kappa * laplacian_phi
    }
    /// Compute the pressure tensor diagonal (isotropic part) for the free-energy model.
    ///
    /// `P_iso = p_bulk + kappa * phi * lap(phi) + kappa/2 * |grad phi|^2`
    pub fn pressure_tensor_isotropic(
        &self,
        rho: f64,
        phi: f64,
        laplacian_phi: f64,
        grad_phi_sq: f64,
    ) -> f64 {
        let p_bulk = CS2 * rho + 0.5 * self.a * phi * phi + 0.75 * self.b * phi.powi(4);
        p_bulk + self.kappa * phi * laplacian_phi + 0.5 * self.kappa * grad_phi_sq
    }
}
/// Shan-Chen multiphase model parameters (2D).
#[derive(Debug, Clone)]
pub struct ShanChenModel {
    /// Coupling constant `G`. Negative values produce attractive
    /// interactions (phase separation).
    pub g_coupling: f64,
    /// Reference density for the pseudo-potential `rho_0`.
    pub rho_0: f64,
}
impl ShanChenModel {
    /// Create a new Shan-Chen model.
    pub fn new(g_coupling: f64, rho_0: f64) -> Self {
        Self { g_coupling, rho_0 }
    }
    /// Pseudo-potential function: `psi(rho) = rho_0 * (1 - exp(-rho / rho_0))`.
    pub fn psi(&self, rho: f64) -> f64 {
        self.rho_0 * (1.0 - (-rho / self.rho_0).exp())
    }
    /// Compute the Shan-Chen interaction force at cell (x, y).
    ///
    /// `F = -G * psi(x,y) * sum_i w_i * psi(x+c_i, y+c_i) * c_i`
    ///
    /// Returns `(Fx, Fy)`.
    pub fn interaction_force(&self, grid: &LbmGrid2D, x: usize, y: usize) -> (f64, f64) {
        let nx = grid.nx;
        let ny = grid.ny;
        let q = grid.lattice.q();
        let k = grid.idx(x, y);
        let psi_here = self.psi(grid.rho[k]);
        let mut fx = 0.0;
        let mut fy = 0.0;
        for i in 1..q {
            let c = grid.lattice.velocity_2d(i);
            let cx = c[0];
            let cy = c[1];
            let nx_i = ((x as i64 + cx as i64).rem_euclid(nx as i64)) as usize;
            let ny_i = ((y as i64 + cy as i64).rem_euclid(ny as i64)) as usize;
            let k_nb = grid.idx(nx_i, ny_i);
            let psi_nb = self.psi(grid.rho[k_nb]);
            let w = grid.lattice.weight(i);
            fx += w * psi_nb * cx as f64;
            fy += w * psi_nb * cy as f64;
        }
        fx *= -self.g_coupling * psi_here;
        fy *= -self.g_coupling * psi_here;
        (fx, fy)
    }
    /// Apply the Shan-Chen forcing to the velocity field.
    ///
    /// The effective velocity used in the equilibrium is shifted:
    /// `u_eq = u + tau * F / rho`
    ///
    /// This modifies `grid.ux` and `grid.uy` in place (call before collision).
    pub fn apply_force(&self, grid: &mut LbmGrid2D, tau: f64) {
        let nx = grid.nx;
        let ny = grid.ny;
        let mut forces = Vec::with_capacity(nx * ny);
        for y in 0..ny {
            for x in 0..nx {
                forces.push(self.interaction_force(grid, x, y));
            }
        }
        for y in 0..ny {
            for x in 0..nx {
                let k = grid.idx(x, y);
                let rho = grid.rho[k];
                if rho.abs() > 1e-15 {
                    let (fx, fy) = forces[k];
                    grid.ux[k] += tau * fx / rho;
                    grid.uy[k] += tau * fy / rho;
                }
            }
        }
    }
}
/// Parameters for a spinodal decomposition simulation.
#[derive(Debug, Clone)]
pub struct SpinodалDecompositionParams {
    /// Cahn-Hilliard mobility M.
    pub mobility: f64,
    /// Landau parameter A (< 0 for unstable homogeneous state).
    pub a: f64,
    /// Landau parameter B (> 0 for stability).
    pub b: f64,
    /// Interface parameter kappa (> 0).
    pub kappa: f64,
    /// Random noise amplitude for initial perturbation.
    pub noise_amplitude: f64,
}
impl SpinodалDecompositionParams {
    /// Create default parameters for spinodal decomposition.
    pub fn new(a: f64, b: f64, kappa: f64, mobility: f64) -> Self {
        Self {
            mobility,
            a,
            b,
            kappa,
            noise_amplitude: 0.01,
        }
    }
    /// Dominant wavenumber at onset: k* = sqrt(-A / (2*kappa)).
    pub fn dominant_wavenumber(&self) -> f64 {
        if self.a < 0.0 && self.kappa > 0.0 {
            (-self.a / (2.0 * self.kappa)).sqrt()
        } else {
            0.0
        }
    }
    /// Dominant wavelength: lambda* = 2*pi / k*.
    pub fn dominant_wavelength(&self) -> f64 {
        let k = self.dominant_wavenumber();
        if k > 1e-30 {
            std::f64::consts::TAU / k
        } else {
            f64::INFINITY
        }
    }
    /// Maximum growth rate sigma_max = M * A^2 / (4 * kappa).
    pub fn max_growth_rate(&self) -> f64 {
        if self.a < 0.0 && self.kappa > 0.0 {
            self.mobility * self.a * self.a / (4.0 * self.kappa)
        } else {
            0.0
        }
    }
    /// Coarseming time scale: t_c = 1 / sigma_max.
    pub fn coarsening_time(&self) -> f64 {
        let sigma = self.max_growth_rate();
        if sigma > 1e-30 {
            1.0 / sigma
        } else {
            f64::INFINITY
        }
    }
}
/// State of a spherically symmetric bubble.
#[derive(Debug, Clone)]
pub struct BubbleState {
    /// Current bubble radius (m).
    pub radius: f64,
    /// Current wall velocity dR/dt (m/s).
    pub radius_rate: f64,
    /// Equilibrium radius (m).
    pub r_eq: f64,
    /// Liquid pressure far from bubble.
    pub p_inf: f64,
    /// Liquid density.
    pub rho_liquid: f64,
    /// Surface tension.
    pub sigma: f64,
    /// Liquid dynamic viscosity.
    pub mu_liquid: f64,
}
impl BubbleState {
    /// Create a new bubble near equilibrium.
    pub fn new(r_eq: f64, p_inf: f64, rho_liquid: f64, sigma: f64, mu_liquid: f64) -> Self {
        Self {
            radius: r_eq,
            radius_rate: 0.0,
            r_eq,
            p_inf,
            rho_liquid,
            sigma,
            mu_liquid,
        }
    }
    /// Gas pressure inside the bubble at an arbitrary radius `radius`
    /// (isothermal gas, polytropic exponent 1).
    fn gas_pressure_at(&self, radius: f64) -> f64 {
        let p0 = self.p_inf + 2.0 * self.sigma / self.r_eq;
        p0 * (self.r_eq / radius.max(1e-30)).powi(3)
    }
    /// Gas pressure inside bubble using isothermal assumption.
    ///
    /// `P_gas = P_inf + 2*sigma/R_eq` at equilibrium, then `P_gas * R_eq^3 = P_gas0 * R^3`.
    pub fn gas_pressure(&self) -> f64 {
        self.gas_pressure_at(self.radius)
    }
    /// Rayleigh-Plesset wall acceleration for an arbitrary state `(radius, radius_rate)`.
    ///
    /// Solves the Rayleigh-Plesset equation for `R_tt`:
    /// `rho * (R * R_tt + (3/2) * R_t^2) = p_gas(R) - p_inf - 2*sigma/R - 4*mu*R_t/R`
    /// i.e. `R_tt = ((p_gas - p_inf - 2*sigma/R - 4*mu*R_t/R) / rho - (3/2) * R_t^2) / R`.
    fn acceleration_at(&self, radius: f64, radius_rate: f64) -> f64 {
        let r = radius.max(1e-30);
        let p_wall = self.gas_pressure_at(r)
            - self.p_inf
            - 2.0 * self.sigma / r
            - 4.0 * self.mu_liquid * radius_rate / r;
        (p_wall / self.rho_liquid.max(1e-30) - 1.5 * radius_rate * radius_rate) / r
    }
    /// Rayleigh-Plesset wall acceleration `d²R/dt²` of the current state.
    ///
    /// The full Rayleigh-Plesset equation including surface tension and the
    /// Newtonian viscous stress of the liquid is used:
    ///
    /// `rho * (R * R_tt + (3/2) * R_t^2) = P_gas - P_inf - 2*sigma/R - 4*mu*R_t/R`
    ///
    /// which is solved for `R_tt`:
    /// `R_tt = ((P_gas - P_inf - 2*sigma/R - 4*mu*R_t/R) / rho - (3/2) * R_t^2) / R`.
    /// The inertia term `(3/2) * R_t^2` is *not* divided by the density.
    pub fn acceleration(&self) -> f64 {
        self.acceleration_at(self.radius, self.radius_rate)
    }
    /// Advance the state `(R, R_t)` by one time step `dt` with the classical
    /// fourth-order Runge-Kutta scheme.
    ///
    /// The step is fourth-order accurate and, for an inviscid bubble
    /// (`mu_liquid = 0`), conserves the Rayleigh-Plesset energy
    /// `2*pi*rho*R^3*R_t^2 - 4*pi*P_gas0*R_eq^3*ln(R) + 4*pi/3*P_inf*R^3 + 4*pi*sigma*R^2`
    /// up to its `O(dt^4)` truncation error.
    pub fn step(&mut self, dt: f64) {
        let r0 = self.radius;
        let v0 = self.radius_rate;
        let k1_r = v0;
        let k1_v = self.acceleration_at(r0, v0);
        let k2_r = v0 + 0.5 * dt * k1_v;
        let k2_v = self.acceleration_at(r0 + 0.5 * dt * k1_r, k2_r);
        let k3_r = v0 + 0.5 * dt * k2_v;
        let k3_v = self.acceleration_at(r0 + 0.5 * dt * k2_r, k3_r);
        let k4_r = v0 + dt * k3_v;
        let k4_v = self.acceleration_at(r0 + dt * k3_r, k4_r);
        self.radius = (r0 + dt / 6.0 * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r)).max(1e-30);
        self.radius_rate = v0 + dt / 6.0 * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);
    }
    /// Collapse time estimation (Rayleigh formula for pressure-driven collapse).
    ///
    /// `t_c = 0.915 * R_0 * sqrt(rho / (P_inf - P_gas))`
    pub fn collapse_time_estimate(&self) -> f64 {
        let p_drive = self.p_inf - self.gas_pressure();
        if p_drive < 1e-30 || self.rho_liquid < 1e-30 {
            return f64::INFINITY;
        }
        0.915 * self.r_eq * (self.rho_liquid / p_drive).sqrt()
    }
}
/// Parameters for a binary Shan-Chen multi-component simulation.
///
/// Component A (index 0) and component B (index 1) interact via a
/// cross-coupling constant `g_cc` (the A-B coupling).  Like-component coupling
/// constants `g_aa` (A-A) and `g_bb` (B-B) are used for self-interactions
/// (density ratio tuning).  The force on component `s` is
/// `F_s(x) = -psi_s(x) * sum_s' g_ss' * sum_i w_i psi_s'(x + c_i) c_i`.
#[derive(Debug, Clone)]
pub struct MultiComponentSC {
    /// Like-component coupling for fluid A (g_aa; typically 0 or small).
    pub g_aa: f64,
    /// Like-component coupling for fluid B (g_bb; typically 0 or small).
    pub g_bb: f64,
    /// Cross-component coupling (g_ab = g_cc).  Negative → immiscibility.
    pub g_cc: f64,
    /// Reference density for the pseudo-potential of each component.
    pub rho_0: f64,
}
impl MultiComponentSC {
    /// Create a new multi-component Shan-Chen model.
    pub fn new(g_aa: f64, g_bb: f64, g_cc: f64, rho_0: f64) -> Self {
        Self {
            g_aa,
            g_bb,
            g_cc,
            rho_0,
        }
    }
    /// Pseudo-potential: `psi(rho) = rho_0 * (1 - exp(-rho/rho_0))`.
    pub fn psi(&self, rho: f64) -> f64 {
        self.rho_0 * (1.0 - (-rho / self.rho_0).exp())
    }
    /// Weighted neighbour sum `S(x) = sum_i w_i * psi(x + c_i) * c_i` of a
    /// pseudo-potential field on a periodic `nx x ny` D2Q9 grid.
    ///
    /// Returns `(S_x, S_y)`.
    fn neighbour_sum(psi: &[f64], x: usize, y: usize, nx: usize, ny: usize) -> (f64, f64) {
        let mut sum_x = 0.0_f64;
        let mut sum_y = 0.0_f64;
        for i in 1..9 {
            let c = D2Q9_VELOCITIES[i];
            let xi = ((x as i64 + c[0] as i64).rem_euclid(nx as i64)) as usize;
            let yi = ((y as i64 + c[1] as i64).rem_euclid(ny as i64)) as usize;
            let w = D2Q9_WEIGHTS[i];
            let psi_nb = psi[yi * nx + xi];
            sum_x += w * psi_nb * c[0] as f64;
            sum_y += w * psi_nb * c[1] as f64;
        }
        (sum_x, sum_y)
    }
    /// Compute the total Shan-Chen interaction force on component A at (x, y).
    ///
    /// `F_A(x) = -psi_A(x) * [ g_aa * sum_i w_i psi_A(x+c_i) c_i + g_cc * sum_i w_i psi_B(x+c_i) c_i ]`
    ///
    /// The like-component (A-A) term is weighted by the pseudo-potential of
    /// component A at the node, and the cross (A-B) term by the same factor,
    /// so that the force on A never depends on `psi_B(x)` itself.
    /// Returns `(Fx_A, Fy_A)`.
    pub fn interaction_force(
        &self,
        psi_a: &[f64],
        psi_b: &[f64],
        x: usize,
        y: usize,
        nx: usize,
        ny: usize,
    ) -> (f64, f64) {
        let psi_a_here = psi_a[y * nx + x];
        let (sa_x, sa_y) = Self::neighbour_sum(psi_a, x, y, nx, ny);
        let (sb_x, sb_y) = Self::neighbour_sum(psi_b, x, y, nx, ny);
        (
            -psi_a_here * (self.g_aa * sa_x + self.g_cc * sb_x),
            -psi_a_here * (self.g_aa * sa_y + self.g_cc * sb_y),
        )
    }
    /// Compute the total Shan-Chen interaction force on component B at (x, y).
    ///
    /// `F_B(x) = -psi_B(x) * [ g_bb * sum_i w_i psi_B(x+c_i) c_i + g_cc * sum_i w_i psi_A(x+c_i) c_i ]`
    ///
    /// This is the mirror image of [`MultiComponentSC::interaction_force`]: it
    /// uses the like-component coupling `g_bb` and the same cross coupling
    /// `g_cc`.  Returns `(Fx_B, Fy_B)`.
    pub fn interaction_force_b(
        &self,
        psi_a: &[f64],
        psi_b: &[f64],
        x: usize,
        y: usize,
        nx: usize,
        ny: usize,
    ) -> (f64, f64) {
        let psi_b_here = psi_b[y * nx + x];
        let (sa_x, sa_y) = Self::neighbour_sum(psi_a, x, y, nx, ny);
        let (sb_x, sb_y) = Self::neighbour_sum(psi_b, x, y, nx, ny);
        (
            -psi_b_here * (self.g_bb * sb_x + self.g_cc * sa_x),
            -psi_b_here * (self.g_bb * sb_y + self.g_cc * sa_y),
        )
    }
    /// Compute the interaction forces on both components at (x, y).
    ///
    /// Returns `((Fx_A, Fy_A), (Fx_B, Fy_B))`.  On a periodic grid the sum of
    /// both forces over all nodes vanishes (momentum conservation).
    pub fn interaction_forces(
        &self,
        psi_a: &[f64],
        psi_b: &[f64],
        x: usize,
        y: usize,
        nx: usize,
        ny: usize,
    ) -> ((f64, f64), (f64, f64)) {
        (
            self.interaction_force(psi_a, psi_b, x, y, nx, ny),
            self.interaction_force_b(psi_a, psi_b, x, y, nx, ny),
        )
    }
    /// Estimate the equilibrium density ratio `rho_A / rho_B` at saturation.
    ///
    /// For the symmetric binary model, the density ratio is controlled by
    /// adjusting `g_cc`.  This returns a heuristic estimate:
    ///
    /// `rho_ratio ~ exp(2 * |g_cc| * rho_0 / 3)`
    pub fn density_ratio_estimate(&self) -> f64 {
        (2.0 * self.g_cc.abs() * self.rho_0 / 3.0).exp()
    }
}
