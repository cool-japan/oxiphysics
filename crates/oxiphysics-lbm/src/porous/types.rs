//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{CX, CY, OPP, W};

/// Analytical Darcy flow model.
///
/// Provides closed-form expressions for Darcy velocity, pressure drop, and the
/// Darcy-regime Reynolds number.
pub struct DarcyFlow {
    /// Intrinsic permeability K \[m^2\].
    pub permeability: f64,
    /// Dynamic viscosity mu \[Pa*s\].
    pub viscosity: f64,
}
impl DarcyFlow {
    /// Create a new `DarcyFlow` with the given permeability and viscosity.
    pub fn new(k: f64, mu: f64) -> Self {
        Self {
            permeability: k,
            viscosity: mu,
        }
    }
    /// Darcy velocity q = -K/mu * grad(P).
    ///
    /// Returns the magnitude of the Darcy flux for a given `pressure_gradient`.
    pub fn velocity(&self, pressure_gradient: f64) -> f64 {
        -self.permeability / self.viscosity * pressure_gradient
    }
    /// Pressure drop dP = mu * L * v / K for a channel of length `length`
    /// carrying Darcy velocity `velocity`.
    pub fn pressure_drop(&self, length: f64, velocity: f64) -> f64 {
        self.viscosity * length * velocity / self.permeability
    }
    /// Pore-scale Reynolds number Re = rho * v * d / mu.
    ///
    /// For Darcy's law to be valid, Re should be << 1.
    pub fn reynolds_number(&self, velocity: f64, pore_size: f64, density: f64) -> f64 {
        density * velocity * pore_size / self.viscosity
    }
}
/// Extended Darcy model including Forchheimer (inertial) correction.
///
/// -grad(P) = (mu/K) * u + (rho * beta / sqrt(K)) * |u| * u
///
/// The Forchheimer coefficient beta is also known as the
/// non-Darcy coefficient or inertial resistance coefficient.
pub struct ForchheimerFlow {
    /// Intrinsic permeability K (m^2).
    pub permeability: f64,
    /// Forchheimer inertial coefficient beta (dimensionless).
    pub beta: f64,
    /// Fluid dynamic viscosity mu (Pa·s).
    pub viscosity: f64,
    /// Fluid density rho (kg/m^3).
    pub density: f64,
}
impl ForchheimerFlow {
    /// Create a new Forchheimer flow model.
    pub fn new(permeability: f64, beta: f64, viscosity: f64, density: f64) -> Self {
        Self {
            permeability,
            beta,
            viscosity,
            density,
        }
    }
    /// Forchheimer number: ratio of inertial to viscous drag.
    ///
    /// Fo = rho * beta * sqrt(K) * |u| / mu
    ///
    /// Darcy regime: Fo << 1; Forchheimer regime: Fo ~ 1.
    pub fn forchheimer_number(&self, velocity: f64) -> f64 {
        self.density * self.beta * self.permeability.sqrt() * velocity.abs() / self.viscosity
    }
    /// Total pressure gradient (Pa/m) for 1-D flow at velocity `u`.
    ///
    /// -dP/dx = (mu/K) * u + (rho * beta / sqrt(K)) * |u| * u
    pub fn pressure_gradient(&self, velocity: f64) -> f64 {
        let darcy_term = self.viscosity / self.permeability * velocity;
        let inertial_term =
            self.density * self.beta / self.permeability.sqrt() * velocity * velocity.abs();
        darcy_term + inertial_term
    }
    /// Ergun-type pressure drop across a packed bed of length L.
    ///
    /// dP = mu * L * u / K + rho * beta / sqrt(K) * L * u^2
    pub fn pressure_drop(&self, velocity: f64, length: f64) -> f64 {
        self.pressure_gradient(velocity) * length
    }
    /// Inertial correction factor C_F (Ergun equation).
    ///
    /// `C_F = 0.143 / phi^0.5`  (Kozeny-Carman based correlation, Ward 1964)
    pub fn inertial_coefficient_ward(phi: f64) -> f64 {
        0.143 / phi.max(1e-12).sqrt()
    }
    /// Effective permeability (m^2) in Ergun correlation:
    ///
    /// K = phi^3 * d_p^2 / (150 * (1 - phi)^2)
    pub fn ergun_permeability(phi: f64, dp: f64) -> f64 {
        let one_minus = (1.0 - phi).max(1e-12);
        phi * phi * phi * dp * dp / (150.0 * one_minus * one_minus)
    }
}
/// D2Q9 LBM solver with Darcy-Brinkman-Forchheimer body force.
///
/// Extends `PorousMediumD2Q9` with the Guo forcing scheme so that
/// inertial (Forchheimer) resistance is included for high-Re porous flow.
pub struct DbfLbmD2Q9 {
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Distribution functions `f[cell][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Node type for each lattice site.
    pub node_types: Vec<NodeType>,
    /// Kinematic viscosity nu.
    pub nu: f64,
    /// Forchheimer body force model (applied at porous nodes).
    pub body_force: Option<ForchheimerBodyForce>,
}
impl DbfLbmD2Q9 {
    /// Create a new DBF-LBM grid of size `nx x ny`, all fluid nodes.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let feq = PorousMediumD2Q9::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![feq; n],
            node_types: (0..n).map(|_| NodeType::Fluid).collect(),
            nu,
            body_force: None,
        }
    }
    /// Set a Forchheimer body force model for porous nodes.
    pub fn with_body_force(mut self, bf: ForchheimerBodyForce) -> Self {
        self.body_force = Some(bf);
        self
    }
    /// Linear index.
    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Mark a node as solid.
    pub fn set_solid(&mut self, i: usize, j: usize) {
        let idx = self.idx(i, j);
        self.node_types[idx] = NodeType::Solid;
    }
    /// Mark a node as porous with porosity `phi`.
    pub fn set_porous(&mut self, i: usize, j: usize, phi: f64) {
        let idx = self.idx(i, j);
        self.node_types[idx] = NodeType::PorousMedia(phi);
    }
    /// BGK relaxation parameter.
    #[inline]
    fn omega(&self) -> f64 {
        1.0 / (3.0 * self.nu + 0.5)
    }
    /// Compute macroscopic density and velocity at a cell.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let fi = &self.f[idx];
        let rho: f64 = fi.iter().sum();
        let ux = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CX[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        let uy = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CY[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        (rho, [ux, uy])
    }
    /// Collision step with optional Guo forcing at porous nodes.
    pub fn collide(&mut self) {
        let omega = self.omega();
        let n = self.nx * self.ny;
        for idx in 0..n {
            match self.node_types[idx] {
                NodeType::Solid => {
                    let (rho, _) = self.macros(idx);
                    self.f[idx] = PorousMediumD2Q9::equilibrium(rho, [0.0, 0.0]);
                }
                NodeType::Fluid => {
                    let (rho, u) = self.macros(idx);
                    let feq = PorousMediumD2Q9::equilibrium(rho, u);
                    for k in 0..9 {
                        self.f[idx][k] += omega * (feq[k] - self.f[idx][k]);
                    }
                }
                NodeType::PorousMedia(phi) => {
                    let (rho, u) = self.macros(idx);
                    let u_eff = [u[0] * phi, u[1] * phi];
                    let feq = PorousMediumD2Q9::equilibrium(rho, u_eff);
                    for k in 0..9 {
                        self.f[idx][k] += omega * (feq[k] - self.f[idx][k]);
                    }
                    if let Some(ref bf) = self.body_force {
                        let fi_force = bf.guo_forcing_d2q9(u[0], u[1], omega);
                        for k in 0..9 {
                            self.f[idx][k] += fi_force[k];
                        }
                    }
                }
            }
        }
    }
    /// Streaming step (pull, periodic).
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for j in 0..ny {
            for i in 0..nx {
                let dst = j * nx + i;
                for k in 0..9 {
                    let si = (i as isize - CX[k] as isize).rem_euclid(nx as isize) as usize;
                    let sj = (j as isize - CY[k] as isize).rem_euclid(ny as isize) as usize;
                    let src = sj * nx + si;
                    if matches!(self.node_types[src], NodeType::Solid) {
                        self.f[dst][k] = f_old[dst][OPP[k]];
                    } else {
                        self.f[dst][k] = f_old[src][k];
                    }
                }
            }
        }
    }
    /// One full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }
    /// Set uniform initial state.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = PorousMediumD2Q9::equilibrium(rho, u);
        for fi in self.f.iter_mut() {
            *fi = feq;
        }
    }
    /// Total mass.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
    /// Average velocity magnitude over non-solid nodes.
    pub fn average_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut sum_u = 0.0;
        let mut count = 0;
        for idx in 0..n {
            if !matches!(self.node_types[idx], NodeType::Solid) {
                let (_, u) = self.macros(idx);
                sum_u += (u[0] * u[0] + u[1] * u[1]).sqrt();
                count += 1;
            }
        }
        if count > 0 { sum_u / count as f64 } else { 0.0 }
    }
}
/// Tortuosity models for porous media.
///
/// Tortuosity `tau` relates the actual diffusion path length to the straight
/// distance.  Effective diffusivity: D_eff = D * phi / tau.
pub struct TortuosityModel;
impl TortuosityModel {
    /// Bruggeman correlation: tau = phi^{-0.5} (for spheres).
    ///
    /// D_eff = D * phi^{1.5}
    pub fn bruggeman(phi: f64) -> f64 {
        phi.max(1e-12).powf(-0.5)
    }
    /// Millington-Quirk (1961) model for gas diffusion:
    ///
    /// tau = phi^{-1/3}
    pub fn millington_quirk(phi: f64) -> f64 {
        phi.max(1e-12).powf(-1.0 / 3.0)
    }
    /// Weissberg (1963) model for sphere packings:
    ///
    /// tau = 1 - 0.5 * ln(phi)
    pub fn weissberg(phi: f64) -> f64 {
        1.0 - 0.5 * phi.max(1e-12).ln()
    }
    /// Effective diffusivity D_eff = D * phi / tau for a given tortuosity model.
    pub fn effective_diffusivity(d: f64, phi: f64, tau: f64) -> f64 {
        d * phi / tau.max(1e-12)
    }
    /// Bruggeman effective diffusivity: D_eff = D * phi^{1.5}.
    pub fn bruggeman_diffusivity(d: f64, phi: f64) -> f64 {
        d * phi.max(0.0).powf(1.5)
    }
}
/// Brinkman effective viscosity model for porous media.
///
/// In Brinkman's extension of Darcy's law, the effective viscosity mu_eff
/// replaces the fluid viscosity mu in viscous stress terms.  Several
/// correlations exist; this struct provides the most common ones.
pub struct BrinkmanEffectiveViscosity;
impl BrinkmanEffectiveViscosity {
    /// Simple inverse-porosity model: mu_eff = mu / phi.
    ///
    /// Used in basic Brinkman–Darcy coupling.
    pub fn inverse_porosity(mu: f64, phi: f64) -> f64 {
        mu / phi.max(1e-30)
    }
    /// Lundgren's (1972) formula for dilute suspensions:
    ///
    /// `mu_eff = mu * (1 + 5/2 * (1 - phi))`
    pub fn lundgren(mu: f64, phi: f64) -> f64 {
        mu * (1.0 + 2.5 * (1.0 - phi))
    }
    /// Brinkman–Hadamard correction (sphere packing):
    ///
    /// `mu_eff = mu * phi / (1 - 1.5*sqrt(1-phi))`
    pub fn brinkman_hadamard(mu: f64, phi: f64) -> f64 {
        let denom = 1.0 - 1.5 * (1.0 - phi).sqrt();
        if denom.abs() < 1e-12 {
            return mu;
        }
        mu * phi / denom
    }
    /// Ratio mu_eff / mu for a given phi.
    pub fn viscosity_ratio(phi: f64) -> f64 {
        1.0 / phi.max(1e-30)
    }
}
/// Analytical results for pressure-driven flow in a porous channel
/// (Brinkman equation in 1D).
///
/// Brinkman equation: mu_eff * d^2u/dy^2 - mu/K * u = dP/dx
///
/// Solution: u(y) = (dP/dx * K / mu) * \[ 1 - cosh(y/sqrt(K)) / cosh(H/sqrt(K)) \]
pub struct PorousChannelFlow {
    /// Dynamic viscosity mu (Pa·s).
    pub mu: f64,
    /// Effective viscosity mu_eff (Pa·s).
    pub mu_eff: f64,
    /// Intrinsic permeability K (m^2).
    pub k: f64,
    /// Channel half-width H (m).
    pub half_width: f64,
}
impl PorousChannelFlow {
    /// Create a new porous channel flow model.
    pub fn new(mu: f64, mu_eff: f64, k: f64, half_width: f64) -> Self {
        Self {
            mu,
            mu_eff,
            k,
            half_width,
        }
    }
    /// Brinkman penetration depth: alpha = sqrt(K / mu_eff).
    pub fn penetration_depth(&self) -> f64 {
        (self.k / self.mu_eff).sqrt()
    }
    /// Velocity profile u(y) for pressure gradient dp_dx.
    ///
    /// y is measured from the channel centre.
    pub fn velocity_profile(&self, y: f64, dp_dx: f64) -> f64 {
        let alpha = self.penetration_depth();
        let u_darcy = -dp_dx * self.k / self.mu;
        let cosh_h = (self.half_width / alpha).cosh();
        let cosh_y = (y / alpha).cosh();
        u_darcy * (1.0 - cosh_y / cosh_h)
    }
    /// Average velocity for pressure gradient dp_dx.
    ///
    /// `u` = u_darcy * \[ 1 - alpha/H * tanh(H/alpha) \]
    pub fn average_velocity(&self, dp_dx: f64) -> f64 {
        let alpha = self.penetration_depth();
        let u_darcy = -dp_dx * self.k / self.mu;
        u_darcy * (1.0 - alpha / self.half_width * (self.half_width / alpha).tanh())
    }
    /// Wall shear stress: tau_w = mu_eff * du/dy|_{y=H}.
    pub fn wall_shear_stress(&self, dp_dx: f64) -> f64 {
        let alpha = self.penetration_depth();
        let u_darcy = -dp_dx * self.k / self.mu;
        self.mu_eff * (-u_darcy / alpha * (self.half_width / alpha).tanh())
    }
    /// Effective Darcy permeability of the channel:
    /// K_eff = `u` * mu / (-dp_dx)  (should recover K for deep penetration)
    pub fn effective_permeability(&self, dp_dx: f64) -> f64 {
        if dp_dx.abs() < 1e-30 {
            return self.k;
        }
        let u_avg = self.average_velocity(dp_dx);
        -u_avg * self.mu / dp_dx
    }
}
/// Classification of each lattice node.
pub enum NodeType {
    /// Regular fluid node -- standard BGK collision.
    Fluid,
    /// Solid obstacle -- full bounce-back, zero velocity.
    Solid,
    /// Porous medium node.  The wrapped value is the porosity phi in (0, 1].
    PorousMedia(f64),
}
/// D2Q9 LBM solver with support for porous media.
pub struct PorousMediumD2Q9 {
    /// Number of lattice nodes in the x direction.
    pub nx: usize,
    /// Number of lattice nodes in the y direction.
    pub ny: usize,
    /// Distribution functions, indexed as `f[idx][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Node type for each lattice site.
    pub node_types: Vec<NodeType>,
    /// Kinematic viscosity nu.
    pub nu: f64,
}
impl PorousMediumD2Q9 {
    /// Create a new grid of size `nx x ny`, all nodes set to `Fluid`.
    /// Distributions are initialised to the rest-state equilibrium (rho = 1, u = 0).
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let f_eq = Self::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![f_eq; n],
            node_types: (0..n).map(|_| NodeType::Fluid).collect(),
            nu,
        }
    }
    /// Linear index for node `(i, j)`.
    #[inline]
    pub(crate) fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Mark node `(i, j)` as a solid obstacle.
    pub fn set_solid(&mut self, i: usize, j: usize) {
        let idx = self.idx(i, j);
        self.node_types[idx] = NodeType::Solid;
    }
    /// Mark node `(i, j)` as porous medium with porosity `porosity` in (0, 1].
    pub fn set_porous(&mut self, i: usize, j: usize, porosity: f64) {
        let idx = self.idx(i, j);
        self.node_types[idx] = NodeType::PorousMedia(porosity);
    }
    /// Compute macroscopic density rho and velocity u at node `idx`.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let fi = &self.f[idx];
        let rho: f64 = fi.iter().sum();
        let ux = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CX[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        let uy = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CY[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        (rho, [ux, uy])
    }
    /// D2Q9 Maxwell-Boltzmann equilibrium distribution.
    pub fn equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
        let cs2 = 1.0 / 3.0;
        let ux = u[0];
        let uy = u[1];
        let u2 = ux * ux + uy * uy;
        let mut feq = [0.0_f64; 9];
        for k in 0..9 {
            let cu = CX[k] as f64 * ux + CY[k] as f64 * uy;
            feq[k] = W[k] * rho * (1.0 + cu / cs2 + cu * cu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
        }
        feq
    }
    /// BGK relaxation parameter omega from kinematic viscosity nu.
    #[inline]
    fn omega(&self) -> f64 {
        1.0 / (3.0 * self.nu + 0.5)
    }
    /// Collision step.
    ///
    /// - `Fluid`: standard BGK.
    /// - `Solid`: distributions set to equilibrium at zero velocity.
    /// - `PorousMedia(phi)`: BGK with effective velocity `u_eff = u * phi` (Brinkman drag).
    pub fn collide(&mut self) {
        let omega = self.omega();
        let n = self.nx * self.ny;
        for idx in 0..n {
            match self.node_types[idx] {
                NodeType::Solid => {
                    let (rho, _) = self.macros(idx);
                    self.f[idx] = Self::equilibrium(rho, [0.0, 0.0]);
                }
                NodeType::Fluid => {
                    let (rho, u) = self.macros(idx);
                    let feq = Self::equilibrium(rho, u);
                    for k in 0..9 {
                        self.f[idx][k] += omega * (feq[k] - self.f[idx][k]);
                    }
                }
                NodeType::PorousMedia(phi) => {
                    let (rho, u) = self.macros(idx);
                    let u_eff = [u[0] * phi, u[1] * phi];
                    let feq = Self::equilibrium(rho, u_eff);
                    for k in 0..9 {
                        self.f[idx][k] += omega * (feq[k] - self.f[idx][k]);
                    }
                }
            }
        }
    }
    /// Streaming step (pull scheme) with periodic boundaries.
    /// Solid neighbours trigger bounce-back: the population is reflected back.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for j in 0..ny {
            for i in 0..nx {
                let dst = j * nx + i;
                for k in 0..9 {
                    let si = (i as isize - CX[k] as isize).rem_euclid(nx as isize) as usize;
                    let sj = (j as isize - CY[k] as isize).rem_euclid(ny as isize) as usize;
                    let src = sj * nx + si;
                    if matches!(self.node_types[src], NodeType::Solid) {
                        self.f[dst][k] = f_old[dst][OPP[k]];
                    } else {
                        self.f[dst][k] = f_old[src][k];
                    }
                }
            }
        }
    }
    /// One full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }
    /// Set all nodes to a uniform density `rho` and velocity `u`.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = Self::equilibrium(rho, u);
        for fi in self.f.iter_mut() {
            *fi = feq;
        }
    }
    /// Estimate permeability using the Kozeny-Carman relation.
    ///
    /// Counts fluid/porous nodes to determine the bulk fluid fraction phi, then
    /// applies K ~ phi^3 / (5*(1-phi)^2) * d^2 with particle diameter d = 1.
    pub fn permeability_estimate(&self) -> f64 {
        let n = self.nx * self.ny;
        let fluid_count = self
            .node_types
            .iter()
            .filter(|nt| matches!(nt, NodeType::Fluid | NodeType::PorousMedia(_)))
            .count();
        let phi = fluid_count as f64 / n as f64;
        let one_minus_phi = 1.0 - phi;
        if one_minus_phi < 1e-12 {
            return phi * phi * phi / (5.0 * 1e-24);
        }
        phi * phi * phi / (5.0 * one_minus_phi * one_minus_phi)
    }
    /// Compute average velocity magnitude over all fluid nodes.
    pub fn average_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut sum_u = 0.0;
        let mut count = 0;
        for idx in 0..n {
            if !matches!(self.node_types[idx], NodeType::Solid) {
                let (_, u) = self.macros(idx);
                sum_u += (u[0] * u[0] + u[1] * u[1]).sqrt();
                count += 1;
            }
        }
        if count > 0 { sum_u / count as f64 } else { 0.0 }
    }
    /// Total mass (sum of all densities).
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
}
/// Representative Elementary Volume analysis for porous media.
///
/// The REV is the smallest volume over which a measurement represents the
/// bulk behavior of the porous medium.  Below the REV size, statistics are
/// not statistically representative; above the REV size, macroscopic
/// homogeneity is assumed.
pub struct RepresentativeElementaryVolume {
    /// Side length of the cubic REV (m).
    pub side_length: f64,
    /// Porosity measured within this REV.
    pub porosity: f64,
    /// Number of pores counted in the REV.
    pub pore_count: usize,
    /// Mean pore radius (m).
    pub mean_pore_radius: f64,
}
impl RepresentativeElementaryVolume {
    /// Create a new REV descriptor.
    pub fn new(side_length: f64, porosity: f64, pore_count: usize, mean_pore_radius: f64) -> Self {
        Self {
            side_length,
            porosity,
            pore_count,
            mean_pore_radius,
        }
    }
    /// Volume of the REV (m^3).
    pub fn volume(&self) -> f64 {
        self.side_length.powi(3)
    }
    /// Fluid (pore) volume within the REV (m^3).
    pub fn pore_volume(&self) -> f64 {
        self.porosity * self.volume()
    }
    /// Solid volume within the REV (m^3).
    pub fn solid_volume(&self) -> f64 {
        (1.0 - self.porosity) * self.volume()
    }
    /// Specific surface area S_v = total pore surface / REV volume (m^-1).
    ///
    /// Approximates pores as spheres: S_v ≈ 3 * n_pores * r^2 / V_rev.
    pub fn specific_surface_area(&self) -> f64 {
        let v = self.volume();
        if v < 1e-30 {
            return 0.0;
        }
        3.0 * self.pore_count as f64 * self.mean_pore_radius * self.mean_pore_radius / v
    }
    /// Hydraulic diameter D_h = 4 * phi / S_v.
    pub fn hydraulic_diameter(&self) -> f64 {
        let sv = self.specific_surface_area();
        if sv < 1e-30 {
            return 0.0;
        }
        4.0 * self.porosity / sv
    }
    /// Kozeny constant from tortuosity `tau` and specific surface area.
    ///
    /// k_oz = phi^3 / (k_kz * (1-phi)^2) where k_kz = 180 standard.
    pub fn kozeny_constant_estimate(&self) -> f64 {
        let sv = self.specific_surface_area();
        if sv < 1e-30 || self.porosity >= 1.0 {
            return 0.0;
        }
        let phi = self.porosity;
        let one_minus_phi = 1.0 - phi;
        phi * phi * phi / (one_minus_phi * one_minus_phi * sv * sv)
    }
    /// Check whether the REV is large enough relative to the pore size.
    ///
    /// Returns `true` if `side_length / mean_pore_radius >= threshold` (default 10).
    pub fn is_sufficient(&self, threshold: f64) -> bool {
        if self.mean_pore_radius < 1e-30 {
            return false;
        }
        self.side_length / self.mean_pore_radius >= threshold
    }
}
/// Apply Darcy-Brinkman-Forchheimer resistance as a body force in LBM.
///
/// The resistance force per unit volume is:
/// F = -( mu/K + rho*C_F*|u|/sqrt(K) ) * u
///
/// This is added as a source term in the LBM collision (e.g., Guo forcing).
pub struct ForchheimerBodyForce {
    /// Intrinsic permeability K (lattice units).
    pub permeability: f64,
    /// Forchheimer coefficient C_F (dimensionless).
    pub forchheimer_coeff: f64,
    /// Fluid kinematic viscosity nu (lattice units).
    pub nu: f64,
    /// Fluid density rho (lattice units).
    pub rho: f64,
}
impl ForchheimerBodyForce {
    /// Create a new Forchheimer body force model.
    pub fn new(permeability: f64, forchheimer_coeff: f64, nu: f64, rho: f64) -> Self {
        Self {
            permeability,
            forchheimer_coeff,
            nu,
            rho,
        }
    }
    /// Compute force vector \[Fx, Fy\] per unit volume at a given velocity \[ux, uy\].
    #[allow(clippy::too_many_arguments)]
    pub fn force(&self, ux: f64, uy: f64) -> [f64; 2] {
        let u_mag = (ux * ux + uy * uy).sqrt();
        let mu = self.rho * self.nu;
        let darcy = mu / self.permeability;
        let forchheimer = self.rho * self.forchheimer_coeff / self.permeability.sqrt() * u_mag;
        let coeff = -(darcy + forchheimer);
        [coeff * ux, coeff * uy]
    }
    /// Compute the D2Q9 Guo forcing distribution for the porous resistance.
    ///
    /// fi_force_i = w_i * (1 - omega/2) * \[ (c_i - u) / cs^2 + (c_i . u) * c_i / cs^4 \] . F
    ///
    /// Returns the 9 forcing terms to be added to the post-collision distributions.
    pub fn guo_forcing_d2q9(&self, ux: f64, uy: f64, omega: f64) -> [f64; 9] {
        let [fx, fy] = self.force(ux, uy);
        let cs2 = 1.0 / 3.0;
        let cs4 = cs2 * cs2;
        let factor = 1.0 - 0.5 * omega;
        let mut fi = [0.0_f64; 9];
        for k in 0..9 {
            let cx = CX[k] as f64;
            let cy = CY[k] as f64;
            let cu = cx * ux + cy * uy;
            let term1 = (cx - ux) * fx + (cy - uy) * fy;
            let term2 = cu * (cx * fx + cy * fy);
            fi[k] = W[k] * factor * (term1 / cs2 + term2 / cs4);
        }
        fi
    }
}
/// Gray lattice Boltzmann method for porous media.
///
/// In the gray LB method, each node has a "grayness" parameter `n_s` in \[0, 1\]:
///   - n_s = 0: fully fluid (transparent)
///   - n_s = 1: fully solid (opaque)
///   - 0 < n_s < 1: partially solid (gray)
///
/// The collision operator is modified:
///   f_i(x, t+1) = (1 - n_s) * \[f_i - omega*(f_i - feq_i)\]  +  n_s * f_{opp}
///
/// This interpolates between standard BGK and bounce-back.
pub struct GrayLatticeBoltzmann {
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Distribution functions.
    pub f: Vec<[f64; 9]>,
    /// Grayness parameter for each node (0 = fluid, 1 = solid).
    pub ns: Vec<f64>,
    /// Kinematic viscosity nu.
    pub nu: f64,
}
impl GrayLatticeBoltzmann {
    /// Create a new gray LB grid, initialised to equilibrium at rest.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let feq = PorousMediumD2Q9::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![feq; n],
            ns: vec![0.0; n],
            nu,
        }
    }
    /// Linear index.
    #[inline]
    pub(crate) fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Set the grayness at node (i, j).
    pub fn set_grayness(&mut self, i: usize, j: usize, ns: f64) {
        let idx = self.idx(i, j);
        self.ns[idx] = ns.clamp(0.0, 1.0);
    }
    /// BGK relaxation parameter.
    #[inline]
    fn omega(&self) -> f64 {
        1.0 / (3.0 * self.nu + 0.5)
    }
    /// Compute macroscopic density and velocity at node idx.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let fi = &self.f[idx];
        let rho: f64 = fi.iter().sum();
        let ux = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CX[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        let uy = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * CY[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        (rho, [ux, uy])
    }
    /// Gray collision step.
    ///
    /// f_i <- (1-ns) * (f_i - omega*(f_i - feq_i)) + ns * f_opp
    pub fn collide(&mut self) {
        let omega = self.omega();
        let n = self.nx * self.ny;
        for idx in 0..n {
            let ns = self.ns[idx];
            let (rho, u) = self.macros(idx);
            let feq = PorousMediumD2Q9::equilibrium(rho, u);
            let mut f_new = [0.0_f64; 9];
            for k in 0..9 {
                let f_bgk = self.f[idx][k] - omega * (self.f[idx][k] - feq[k]);
                let f_bb = self.f[idx][OPP[k]];
                f_new[k] = (1.0 - ns) * f_bgk + ns * f_bb;
            }
            self.f[idx] = f_new;
        }
    }
    /// Streaming step (pull scheme, periodic).
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for j in 0..ny {
            for i in 0..nx {
                let dst = j * nx + i;
                for k in 0..9 {
                    let si = (i as isize - CX[k] as isize).rem_euclid(nx as isize) as usize;
                    let sj = (j as isize - CY[k] as isize).rem_euclid(ny as isize) as usize;
                    let src = sj * nx + si;
                    self.f[dst][k] = f_old[src][k];
                }
            }
        }
    }
    /// One full step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }
    /// Total mass.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
    /// Average velocity magnitude over fluid-like nodes (ns < 0.5).
    pub fn average_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut sum_u = 0.0;
        let mut count = 0;
        for idx in 0..n {
            if self.ns[idx] < 0.5 {
                let (_, u) = self.macros(idx);
                sum_u += (u[0] * u[0] + u[1] * u[1]).sqrt();
                count += 1;
            }
        }
        if count > 0 { sum_u / count as f64 } else { 0.0 }
    }
    /// Set uniform initial state.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = PorousMediumD2Q9::equilibrium(rho, u);
        for fi in self.f.iter_mut() {
            *fi = feq;
        }
    }
}
/// Partial bounce-back model for porous media.
///
/// In the partial bounce-back method, a fraction `sigma` of each distribution
/// is reflected (bounce-back) and the complement `(1 - sigma)` passes through:
///
///   f_i(x, t+1) = (1-sigma) * f_i_streamed + sigma * f_opp_pre_stream
///
/// The parameter `sigma` is related to the permeability:
///   sigma = (tau - 0.5) / (tau + 0.5 * K^{-1})
///
/// Higher sigma means lower permeability (more reflection).
pub struct PartialBounceBack {
    /// Reflection coefficient sigma in \[0, 1\].
    pub sigma: f64,
}
impl PartialBounceBack {
    /// Create from a reflection coefficient directly.
    pub fn new(sigma: f64) -> Self {
        Self {
            sigma: sigma.clamp(0.0, 1.0),
        }
    }
    /// Create from permeability and relaxation time.
    ///
    /// sigma = (tau - 0.5) / (tau - 0.5 + K)  (simplified relation)
    pub fn from_permeability(permeability: f64, tau: f64) -> Self {
        let num = tau - 0.5;
        let den = num + permeability;
        let sigma = if den.abs() > 1e-15 { num / den } else { 1.0 };
        Self::new(sigma)
    }
    /// Apply partial bounce-back to a single cell's distributions.
    ///
    /// `f_pre` is the pre-streaming distribution, `f_streamed` is the
    /// post-streaming distribution.  The result is stored in `f_out`.
    pub fn apply(&self, f_pre: &[f64; 9], f_streamed: &[f64; 9], f_out: &mut [f64; 9]) {
        let s = self.sigma;
        for k in 0..9 {
            f_out[k] = (1.0 - s) * f_streamed[k] + s * f_pre[OPP[k]];
        }
    }
    /// Effective permeability implied by this sigma value and tau.
    pub fn effective_permeability(&self, tau: f64) -> f64 {
        if self.sigma.abs() < 1e-15 {
            return f64::INFINITY;
        }
        (tau - 0.5) * (1.0 - self.sigma) / self.sigma
    }
}
/// Unified porous media property calculator.
///
/// Provides methods for computing permeability (Kozeny-Carman),
/// Forchheimer inertia corrections, and effective diffusivity using
/// tortuosity models, all in a single convenient struct.
#[allow(dead_code)]
pub struct PorousMedia {
    /// Porosity φ ∈ (0, 1].
    pub porosity: f64,
    /// Effective particle/grain diameter d_p (m).
    pub particle_diameter: f64,
    /// Free (molecular) diffusivity D₀ (m²/s).
    pub free_diffusivity: f64,
}
impl PorousMedia {
    /// Create a new `PorousMedia` instance.
    ///
    /// # Arguments
    /// - `porosity`          — void fraction φ ∈ (0, 1]
    /// - `particle_diameter` — characteristic grain size d_p (m)
    /// - `free_diffusivity`  — molecular diffusivity D₀ (m²/s)
    #[allow(dead_code)]
    pub fn new(porosity: f64, particle_diameter: f64, free_diffusivity: f64) -> Self {
        Self {
            porosity,
            particle_diameter,
            free_diffusivity,
        }
    }
    /// Compute the intrinsic permeability using the Kozeny-Carman relation.
    ///
    /// The classical Kozeny-Carman equation:
    ///
    /// `K = (d_p² * φ³) / (180 * (1 − φ)²)`
    ///
    /// where the constant 180 = 36 × k_KC with Kozeny constant k_KC = 5.
    ///
    /// Returns permeability K (m²). Returns 0.0 for φ ≥ 1.
    ///
    /// # Arguments
    /// - `kozeny_const` — Kozeny-Blake constant (default 180; pass 0.0 to use
    ///   the default value of 180)
    #[allow(dead_code)]
    pub fn compute_kozeny_carman(&self, kozeny_const: f64) -> f64 {
        let phi = self.porosity;
        if phi <= 0.0 || phi >= 1.0 {
            return 0.0;
        }
        let c_kc = if kozeny_const > 0.0 {
            kozeny_const
        } else {
            180.0
        };
        let dp = self.particle_diameter;
        let one_minus_phi = 1.0 - phi;
        (dp * dp * phi * phi * phi) / (c_kc * one_minus_phi * one_minus_phi)
    }
    /// Apply the Forchheimer inertia correction to Darcy's law.
    ///
    /// At high Reynolds numbers the pressure gradient–velocity relationship
    /// deviates from linear (Darcy) behaviour. The Forchheimer extension is:
    ///
    /// `−∇P = (μ/K) * u + β_F * ρ * u²`
    ///
    /// where `β_F = C_F / sqrt(K)` is the Forchheimer coefficient and
    /// `C_F ≈ 0.55` is the Ergun inertial constant.
    ///
    /// This method returns the **total** pressure gradient magnitude (Pa/m)
    /// for a given velocity magnitude.
    ///
    /// # Arguments
    /// - `velocity`   — superficial (Darcy) velocity magnitude |u| (m/s)
    /// - `viscosity`  — dynamic viscosity μ (Pa·s)
    /// - `density`    — fluid density ρ (kg/m³)
    /// - `c_f`        — Forchheimer coefficient C_F (dimensionless; default ~0.55
    ///   if ≤ 0)
    #[allow(dead_code)]
    pub fn compute_forchheimer_correction(
        &self,
        velocity: f64,
        viscosity: f64,
        density: f64,
        c_f: f64,
    ) -> f64 {
        let k = self.compute_kozeny_carman(0.0);
        if k <= 0.0 {
            return 0.0;
        }
        let cf = if c_f > 0.0 { c_f } else { 0.55 };
        let darcy_term = viscosity / k * velocity;
        let beta_f = cf / k.sqrt();
        let forchheimer_term = beta_f * density * velocity * velocity;
        darcy_term + forchheimer_term
    }
    /// Compute the effective diffusivity using the Bruggeman tortuosity model.
    ///
    /// The effective diffusivity in a porous medium accounts for two effects:
    /// 1. Reduced cross-sectional area (porosity factor φ).
    /// 2. Increased path length (tortuosity τ > 1).
    ///
    /// `D_eff = D₀ * φ / τ`
    ///
    /// Supported tortuosity models:
    /// - `"bruggeman"` (default): τ = φ^(−0.5)  →  D_eff = D₀ * φ^1.5
    /// - `"millington_quirk"`:    τ = φ^(−1/3)  →  D_eff = D₀ * φ^(7/3)
    /// - `"weissberg"`:           τ = 1 − 0.5*ln(φ)
    /// - `"constant"`:            τ from `tortuosity_value` argument
    ///
    /// # Arguments
    /// - `model`           — tortuosity model name (see above)
    /// - `tortuosity_value`— used only for `model = "constant"`
    ///
    /// Returns D_eff (m²/s).
    #[allow(dead_code)]
    pub fn compute_effective_diffusivity(&self, model: &str, tortuosity_value: f64) -> f64 {
        let phi = self.porosity.max(1e-30);
        let tau = match model {
            "millington_quirk" => TortuosityModel::millington_quirk(phi),
            "weissberg" => TortuosityModel::weissberg(phi),
            "constant" => tortuosity_value.max(1.0),
            _ => TortuosityModel::bruggeman(phi),
        };
        TortuosityModel::effective_diffusivity(self.free_diffusivity, phi, tau)
    }
}
/// 2D permeability tensor (symmetric 2x2) for anisotropic porous media.
///
/// Stored as \[K_xx, K_xy, K_yy\] (symmetric upper triangle).
#[derive(Debug, Clone, Copy)]
pub struct PermeabilityTensor2D {
    /// K_xx component.
    pub kxx: f64,
    /// K_xy = K_yx component.
    pub kxy: f64,
    /// K_yy component.
    pub kyy: f64,
}
impl PermeabilityTensor2D {
    /// Create a new permeability tensor.
    pub fn new(kxx: f64, kxy: f64, kyy: f64) -> Self {
        Self { kxx, kxy, kyy }
    }
    /// Isotropic tensor: K_xx = K_yy = k, K_xy = 0.
    pub fn isotropic(k: f64) -> Self {
        Self {
            kxx: k,
            kxy: 0.0,
            kyy: k,
        }
    }
    /// Apply tensor to a pressure gradient vector \[dpx, dpy\].
    ///
    /// Darcy velocity: u = -(K/mu) * grad(P)
    pub fn darcy_velocity(&self, dpx: f64, dpy: f64, mu: f64) -> [f64; 2] {
        let ux = -(self.kxx * dpx + self.kxy * dpy) / mu;
        let uy = -(self.kxy * dpx + self.kyy * dpy) / mu;
        [ux, uy]
    }
    /// Determinant of the permeability tensor.
    pub fn determinant(&self) -> f64 {
        self.kxx * self.kyy - self.kxy * self.kxy
    }
    /// Geometric mean permeability: sqrt(det(K)).
    pub fn geometric_mean(&self) -> f64 {
        self.determinant().max(0.0).sqrt()
    }
    /// Effective scalar permeability (harmonic mean of principal values).
    ///
    /// lambda_1, lambda_2 = eigenvalues of K.
    /// K_eff = 2 * lambda_1 * lambda_2 / (lambda_1 + lambda_2)
    pub fn harmonic_mean(&self) -> f64 {
        let trace = self.kxx + self.kyy;
        let disc = ((self.kxx - self.kyy).powi(2) + 4.0 * self.kxy * self.kxy).sqrt();
        let lam1 = 0.5 * (trace + disc);
        let lam2 = 0.5 * (trace - disc);
        let sum = lam1 + lam2;
        if sum.abs() < 1e-30 {
            return 0.0;
        }
        2.0 * lam1 * lam2 / sum
    }
}
/// Darcy-Brinkman-Forchheimer (DBF) porous media model.
///
/// Extends the Darcy model with viscous (Brinkman) and inertial
/// (Forchheimer) corrections:
///
///   -grad(p) = (mu/K) * u + (rho * C_F / sqrt(K)) * |u| * u - mu_eff * laplacian(u)
///
/// where:
///   - K is the intrinsic permeability
///   - C_F is the Forchheimer coefficient (dimensionless)
///   - mu_eff is the effective (Brinkman) viscosity
pub struct DarcyBrinkmanForchheimer {
    /// Intrinsic permeability K.
    pub permeability: f64,
    /// Forchheimer coefficient C_F (dimensionless, typically 0.05-0.55).
    pub forchheimer_coeff: f64,
    /// Fluid dynamic viscosity mu.
    pub viscosity: f64,
    /// Fluid density rho.
    pub density: f64,
    /// Porosity phi.
    pub porosity: f64,
}
impl DarcyBrinkmanForchheimer {
    /// Create a new DBF model.
    pub fn new(
        permeability: f64,
        forchheimer_coeff: f64,
        viscosity: f64,
        density: f64,
        porosity: f64,
    ) -> Self {
        Self {
            permeability,
            forchheimer_coeff,
            viscosity,
            density,
            porosity,
        }
    }
    /// Compute the Darcy drag force per unit volume: (mu/K) * u.
    pub fn darcy_drag(&self, velocity: f64) -> f64 {
        self.viscosity / self.permeability * velocity
    }
    /// Compute the Forchheimer (inertial) drag per unit volume:
    ///   (rho * C_F / sqrt(K)) * |u| * u
    pub fn forchheimer_drag(&self, velocity: f64) -> f64 {
        self.density * self.forchheimer_coeff / self.permeability.sqrt() * velocity * velocity.abs()
    }
    /// Total resistance force per unit volume (Darcy + Forchheimer).
    pub fn total_resistance(&self, velocity: f64) -> f64 {
        self.darcy_drag(velocity) + self.forchheimer_drag(velocity)
    }
    /// Compute the Forchheimer number: Fo = C_F * sqrt(K) * rho * |u| / mu.
    ///
    /// When Fo >> 1, inertial effects dominate; Fo << 1 means Darcy regime.
    pub fn forchheimer_number(&self, velocity: f64) -> f64 {
        self.forchheimer_coeff * self.permeability.sqrt() * self.density * velocity.abs()
            / self.viscosity
    }
    /// Effective viscosity ratio (Brinkman correction): mu_eff / mu.
    ///
    /// Common approximation: mu_eff = mu / porosity.
    pub fn effective_viscosity_ratio(&self) -> f64 {
        1.0 / self.porosity
    }
    /// Compute the LBM body force for a single cell to represent the
    /// porous resistance at the given velocity.
    ///
    /// Returns the force per unit volume (to be applied as a source term).
    pub fn lbm_body_force(&self, ux: f64, uy: f64) -> [f64; 2] {
        let u_mag = (ux * ux + uy * uy).sqrt();
        let darcy = self.viscosity / self.permeability;
        let forch = self.density * self.forchheimer_coeff / self.permeability.sqrt() * u_mag;
        let total = darcy + forch;
        [-total * ux, -total * uy]
    }
}
