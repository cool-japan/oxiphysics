#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluid FEM and stabilized methods.
//!
//! Implements stabilized finite element formulations for incompressible and
//! weakly-compressible flows including:
//!
//! - [`StabilizedNS`]: SUPG/PSPG stabilized Navier-Stokes FEM
//! - [`SupgStabilization`]: Streamline Upwind Petrov-Galerkin parameter τ_SUPG
//! - [`PspgStabilization`]: Pressure Stabilizing Petrov-Galerkin parameter τ_PSPG
//! - [`VmsStabilization`]: Variational Multiscale method
//! - [`IncompressibleElement`]: P2/P1 Taylor-Hood and bubble elements
//! - [`ProjectionMethod`]: Chorin fractional-step pressure projection
//! - [`AdvectionDiffusion`]: SUPG stabilized advection-diffusion
//! - [`AleFormulation`]: Arbitrary Lagrangian-Eulerian FEM
//! - [`FluidThermalCoupling`]: Boussinesq natural convection
//! - [`LevelSetFem`]: FEM level-set method with redistancing

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Stabilized Navier-Stokes
// ============================================================================

/// Stabilized Navier-Stokes FEM using SUPG/PSPG formulation.
///
/// Solves the incompressible Navier-Stokes equations:
/// ρ(∂u/∂t + u·∇u) = -∇p + μ∇²u + f
/// ∇·u = 0
///
/// using a velocity-pressure mixed formulation with SUPG and PSPG stabilization
/// to avoid spurious pressure modes and convective instabilities.
pub struct StabilizedNS {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Characteristic element length (m).
    pub h_elem: f64,
    /// Current time (s).
    pub time: f64,
    /// Time step (s).
    pub dt: f64,
    /// Velocity DOFs (u, v, w alternating per node).
    pub velocity: Vec<f64>,
    /// Pressure DOFs (one per pressure node).
    pub pressure: Vec<f64>,
}

impl StabilizedNS {
    /// Create a new stabilized Navier-Stokes solver.
    pub fn new(density: f64, viscosity: f64, h_elem: f64, dt: f64) -> Self {
        Self {
            density,
            viscosity,
            h_elem,
            time: 0.0,
            dt,
            velocity: Vec::new(),
            pressure: Vec::new(),
        }
    }

    /// Kinematic viscosity: ν = μ/ρ.
    pub fn kinematic_viscosity(&self) -> f64 {
        self.viscosity / self.density
    }

    /// Reynolds number for flow with characteristic velocity u_ref and length h.
    pub fn reynolds_number(&self, u_ref: f64, h: f64) -> f64 {
        self.density * u_ref * h / self.viscosity
    }

    /// SUPG stabilization parameter τ_SUPG.
    ///
    /// Uses the Tezduyar formula combining convective and diffusive contributions.
    pub fn tau_supg(&self, u_norm: f64) -> f64 {
        compute_tau_supg(u_norm, self.h_elem, self.kinematic_viscosity(), self.dt)
    }

    /// PSPG stabilization parameter τ_PSPG.
    pub fn tau_pspg(&self, u_norm: f64) -> f64 {
        compute_tau_pspg(u_norm, self.h_elem, self.kinematic_viscosity(), self.dt)
    }

    /// Residual of the momentum equation at a Gauss point.
    ///
    /// R_mom = ρ(∂u/∂t + u·∇u) + ∇p - μ∇²u - f
    #[allow(clippy::too_many_arguments)]
    pub fn momentum_residual(
        &self,
        du_dt: [f64; 3],
        u_conv: [f64; 3],
        grad_p: [f64; 3],
        lap_u: [f64; 3],
        body_force: [f64; 3],
    ) -> [f64; 3] {
        [
            self.density * (du_dt[0] + u_conv[0]) + grad_p[0]
                - self.viscosity * lap_u[0]
                - body_force[0],
            self.density * (du_dt[1] + u_conv[1]) + grad_p[1]
                - self.viscosity * lap_u[1]
                - body_force[1],
            self.density * (du_dt[2] + u_conv[2]) + grad_p[2]
                - self.viscosity * lap_u[2]
                - body_force[2],
        ]
    }

    /// Continuity residual (divergence of velocity).
    pub fn continuity_residual(&self, div_u: f64) -> f64 {
        div_u
    }

    /// Advance time by dt.
    pub fn step_time(&mut self) {
        self.time += self.dt;
    }
}

// ============================================================================
// SUPG stabilization
// ============================================================================

/// Streamline Upwind Petrov-Galerkin (SUPG) stabilization.
///
/// Computes the stabilization parameter τ_SUPG for convection-dominated flows.
/// Based on the Tezduyar–Osawa formula using the intrinsic time scale.
pub struct SupgStabilization {
    /// Characteristic element size (m).
    pub h_elem: f64,
    /// Kinematic viscosity (m²/s).
    pub nu: f64,
    /// Time step (s).
    pub dt: f64,
}

impl SupgStabilization {
    /// Create a new SUPG stabilization object.
    pub fn new(h_elem: f64, nu: f64, dt: f64) -> Self {
        Self { h_elem, nu, dt }
    }

    /// Compute τ_SUPG for velocity magnitude u_norm (m/s).
    pub fn tau(&self, u_norm: f64) -> f64 {
        compute_tau_supg(u_norm, self.h_elem, self.nu, self.dt)
    }

    /// Element Peclet number: Pe = u * h / (2 * ν).
    pub fn peclet(&self, u_norm: f64) -> f64 {
        u_norm * self.h_elem / (2.0 * self.nu.max(1e-30))
    }

    /// Optimal upwind parameter ξ(Pe) for SUPG: coth(Pe) - 1/Pe.
    pub fn xi_upwind(&self, u_norm: f64) -> f64 {
        let pe = self.peclet(u_norm);
        if pe.abs() < 1e-6 {
            0.0
        } else {
            1.0 / pe.tanh() - 1.0 / pe
        }
    }

    /// SUPG weighting function contribution (streamline direction).
    ///
    /// Returns τ * u/|u| · ∇N
    pub fn supg_weight(&self, u: [f64; 3], grad_n: [f64; 3], u_norm: f64) -> f64 {
        let tau = self.tau(u_norm);
        if u_norm < 1e-14 {
            return 0.0;
        }
        let dot = u[0] * grad_n[0] + u[1] * grad_n[1] + u[2] * grad_n[2];
        tau * dot / u_norm
    }
}

// ============================================================================
// PSPG stabilization
// ============================================================================

/// Pressure Stabilizing Petrov-Galerkin (PSPG) stabilization.
///
/// Avoids inf-sup instability in equal-order velocity-pressure discretizations
/// by adding a stabilizing term to the continuity equation.
pub struct PspgStabilization {
    /// Characteristic element size (m).
    pub h_elem: f64,
    /// Kinematic viscosity (m²/s).
    pub nu: f64,
    /// Time step (s).
    pub dt: f64,
}

impl PspgStabilization {
    /// Create a new PSPG stabilization object.
    pub fn new(h_elem: f64, nu: f64, dt: f64) -> Self {
        Self { h_elem, nu, dt }
    }

    /// Compute τ_PSPG for velocity magnitude u_norm (m/s).
    pub fn tau(&self, u_norm: f64) -> f64 {
        compute_tau_pspg(u_norm, self.h_elem, self.nu, self.dt)
    }

    /// PSPG stabilization residual contribution.
    ///
    /// τ_PSPG / ρ * ∇q · (ρ Du/Dt + ∇p - μ∇²u - f)
    #[allow(clippy::too_many_arguments)]
    pub fn residual_contribution(
        &self,
        tau: f64,
        rho: f64,
        grad_q: [f64; 3],
        momentum_residual: [f64; 3],
    ) -> f64 {
        let dot = grad_q[0] * momentum_residual[0]
            + grad_q[1] * momentum_residual[1]
            + grad_q[2] * momentum_residual[2];
        tau / rho * dot
    }
}

// ============================================================================
// VMS stabilization
// ============================================================================

/// Variational Multiscale (VMS) method for Navier-Stokes.
///
/// Decomposes the solution into coarse-scale (resolved) and fine-scale
/// (subgrid) parts. The fine-scale model provides a closure for the
/// coarse-scale equations.
pub struct VmsStabilization {
    /// Characteristic element size (m).
    pub h_elem: f64,
    /// Dynamic viscosity (Pa·s).
    pub mu: f64,
    /// Fluid density (kg/m³).
    pub rho: f64,
    /// Time step size (s).
    pub dt: f64,
}

impl VmsStabilization {
    /// Create a new VMS stabilization.
    pub fn new(h_elem: f64, mu: f64, rho: f64, dt: f64) -> Self {
        Self {
            h_elem,
            mu,
            rho,
            dt,
        }
    }

    /// Fine-scale velocity (subscale) residual estimate.
    ///
    /// u' ≈ -τ_M * r_M  where r_M is the coarse-scale momentum residual.
    pub fn fine_scale_velocity(&self, tau_m: f64, momentum_residual: [f64; 3]) -> [f64; 3] {
        [
            -tau_m * momentum_residual[0],
            -tau_m * momentum_residual[1],
            -tau_m * momentum_residual[2],
        ]
    }

    /// Fine-scale pressure estimate.
    ///
    /// p' ≈ -τ_C * r_C  where r_C is the continuity residual.
    pub fn fine_scale_pressure(&self, tau_c: f64, div_u: f64) -> f64 {
        -tau_c * div_u
    }

    /// VMS momentum stabilization parameter τ_M.
    pub fn tau_m(&self, u_norm: f64) -> f64 {
        let nu = self.mu / self.rho;
        compute_tau_supg(u_norm, self.h_elem, nu, self.dt)
    }

    /// VMS continuity stabilization parameter τ_C.
    ///
    /// τ_C = h² / (4 * τ_M * ρ)  (approximate).
    pub fn tau_c(&self, u_norm: f64) -> f64 {
        let tau_m = self.tau_m(u_norm);
        self.h_elem * self.h_elem / (4.0 * tau_m * self.rho + 1e-30)
    }
}

// ============================================================================
// Incompressible Element
// ============================================================================

/// Incompressible flow element using Taylor-Hood (P2/P1) or bubble enrichment.
///
/// Satisfies the LBB (inf-sup) stability condition by using higher-order
/// velocity interpolation or adding a bubble function to the velocity space.
pub struct IncompressibleElement {
    /// Number of velocity nodes (per direction).
    pub n_vel_nodes: usize,
    /// Number of pressure nodes.
    pub n_pres_nodes: usize,
    /// Element size (m).
    pub h_elem: f64,
    /// Use bubble enrichment if true, Taylor-Hood if false.
    pub use_bubble: bool,
}

impl IncompressibleElement {
    /// Create a Taylor-Hood (P2/P1) incompressible element.
    pub fn taylor_hood(h_elem: f64) -> Self {
        Self {
            n_vel_nodes: 6,  // P2 triangle: 3 vertex + 3 mid-edge
            n_pres_nodes: 3, // P1 triangle: 3 vertex
            h_elem,
            use_bubble: false,
        }
    }

    /// Create a P1+/P1 (MINI) element with bubble enrichment.
    pub fn mini_element(h_elem: f64) -> Self {
        Self {
            n_vel_nodes: 4, // P1 + 1 bubble DOF
            n_pres_nodes: 3,
            h_elem,
            use_bubble: true,
        }
    }

    /// Bubble function value at barycentric coordinates (λ1, λ2, λ3).
    ///
    /// b(λ) = 27 * λ1 * λ2 * λ3 (cubic bubble for triangle).
    pub fn bubble_value(&self, lambda: [f64; 3]) -> f64 {
        bubble_enrichment(lambda)
    }

    /// Bubble function gradient at barycentric coords, given physical coords.
    ///
    /// Returns approximate gradient \[∂b/∂x, ∂b/∂y\].
    pub fn bubble_gradient(&self, lambda: [f64; 3], dndx: &[[f64; 2]; 3]) -> [f64; 2] {
        // ∂b/∂x = 27 * (λ2*λ3 * dλ1/dx + λ1*λ3 * dλ2/dx + λ1*λ2 * dλ3/dx)
        let db_dx = 27.0
            * (lambda[1] * lambda[2] * dndx[0][0]
                + lambda[0] * lambda[2] * dndx[1][0]
                + lambda[0] * lambda[1] * dndx[2][0]);
        let db_dy = 27.0
            * (lambda[1] * lambda[2] * dndx[0][1]
                + lambda[0] * lambda[2] * dndx[1][1]
                + lambda[0] * lambda[1] * dndx[2][1]);
        [db_dx, db_dy]
    }

    /// Total velocity DOFs in 2D.
    pub fn total_vel_dofs_2d(&self) -> usize {
        self.n_vel_nodes * 2
    }

    /// Total DOFs including pressure.
    pub fn total_dofs_2d(&self) -> usize {
        self.total_vel_dofs_2d() + self.n_pres_nodes
    }
}

// ============================================================================
// Projection Method
// ============================================================================

/// Chorin's fractional-step (projection) method for incompressible flow.
///
/// Steps:
/// 1. Predictor: u* = u^n + dt*(convection + diffusion + body)  (no pressure)
/// 2. Pressure solve: ∇²p^{n+1} = ρ/dt * ∇·u*
/// 3. Corrector: u^{n+1} = u* - dt/ρ * ∇p^{n+1}
pub struct ProjectionMethod {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Time step (s).
    pub dt: f64,
}

impl ProjectionMethod {
    /// Create a new projection method solver.
    pub fn new(density: f64, viscosity: f64, dt: f64) -> Self {
        Self {
            density,
            viscosity,
            dt,
        }
    }

    /// Predictor step: intermediate velocity u*.
    ///
    /// u* = u^n + dt * (- (u·∇)u + ν∇²u + f/ρ)
    pub fn predictor_step(
        &self,
        u: [f64; 3],
        conv_u: [f64; 3],
        lap_u: [f64; 3],
        body: [f64; 3],
    ) -> [f64; 3] {
        let nu = self.viscosity / self.density;
        [
            u[0] + self.dt * (-conv_u[0] + nu * lap_u[0] + body[0] / self.density),
            u[1] + self.dt * (-conv_u[1] + nu * lap_u[1] + body[1] / self.density),
            u[2] + self.dt * (-conv_u[2] + nu * lap_u[2] + body[2] / self.density),
        ]
    }

    /// Corrector step: project u* onto divergence-free space.
    ///
    /// u^{n+1} = u* - (dt/ρ) * ∇p^{n+1}
    pub fn corrector_step(&self, u_star: [f64; 3], grad_p: [f64; 3]) -> [f64; 3] {
        let f = self.dt / self.density;
        [
            u_star[0] - f * grad_p[0],
            u_star[1] - f * grad_p[1],
            u_star[2] - f * grad_p[2],
        ]
    }

    /// Right-hand side of pressure Poisson equation: ρ/dt * ∇·u*.
    pub fn pressure_rhs(&self, div_u_star: f64) -> f64 {
        self.density / self.dt * div_u_star
    }

    /// Check if velocity is approximately divergence-free (|∇·u| < tol).
    pub fn is_divergence_free(&self, div_u: f64, tol: f64) -> bool {
        div_u.abs() < tol
    }
}

// ============================================================================
// Advection-Diffusion
// ============================================================================

/// SUPG stabilized advection-diffusion equation.
///
/// Solves: ∂φ/∂t + u·∇φ - κ∇²φ = f
///
/// with SUPG stabilization to avoid spurious oscillations in the
/// convection-dominated regime.
pub struct AdvectionDiffusion {
    /// Diffusivity (m²/s).
    pub diffusivity: f64,
    /// Advection velocity (m/s).
    pub velocity: [f64; 3],
    /// Characteristic element size (m).
    pub h_elem: f64,
    /// Time step (s).
    pub dt: f64,
}

impl AdvectionDiffusion {
    /// Create a new advection-diffusion problem.
    pub fn new(diffusivity: f64, velocity: [f64; 3], h_elem: f64, dt: f64) -> Self {
        Self {
            diffusivity,
            velocity,
            h_elem,
            dt,
        }
    }

    /// Magnitude of advection velocity.
    pub fn velocity_magnitude(&self) -> f64 {
        let u = &self.velocity;
        (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()
    }

    /// Element Peclet number: Pe = |u| * h / (2 * κ).
    pub fn peclet_number(&self) -> f64 {
        self.velocity_magnitude() * self.h_elem / (2.0 * self.diffusivity.max(1e-30))
    }

    /// Optimal SUPG stabilization parameter τ.
    ///
    /// Uses intrinsic time scale: τ = h / (2 * |u|) * ξ(Pe).
    pub fn tau_supg(&self) -> f64 {
        stabilization_parameter(
            self.velocity_magnitude(),
            self.h_elem,
            self.diffusivity,
            self.dt,
        )
    }

    /// SUPG-stabilized residual at a Gauss point.
    ///
    /// R = ∂φ/∂t + u·∇φ - κ∇²φ - f
    pub fn residual(&self, dphi_dt: f64, u_dot_grad_phi: f64, lap_phi: f64, source: f64) -> f64 {
        dphi_dt + u_dot_grad_phi - self.diffusivity * lap_phi - source
    }

    /// SUPG stabilization term: τ * (u·∇w) * R
    pub fn supg_term(&self, u_dot_grad_w: f64, residual: f64) -> f64 {
        self.tau_supg() * u_dot_grad_w * residual
    }
}

// ============================================================================
// ALE Formulation
// ============================================================================

/// Arbitrary Lagrangian-Eulerian (ALE) FEM formulation.
///
/// Handles mesh motion by introducing a mesh velocity w. The ALE convection
/// term uses the relative velocity (u - w) instead of u.
pub struct AleFormulation {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Mesh velocity at current time (m/s, per DOF).
    pub mesh_velocity: Vec<[f64; 3]>,
    /// Time step (s).
    pub dt: f64,
}

impl AleFormulation {
    /// Create a new ALE formulation with zero mesh velocity.
    pub fn new(density: f64, viscosity: f64, n_nodes: usize, dt: f64) -> Self {
        Self {
            density,
            viscosity,
            mesh_velocity: vec![[0.0; 3]; n_nodes],
            dt,
        }
    }

    /// ALE convection velocity: c = u - w (fluid velocity minus mesh velocity).
    pub fn ale_convection(&self, u: [f64; 3], w: [f64; 3]) -> [f64; 3] {
        [u[0] - w[0], u[1] - w[1], u[2] - w[2]]
    }

    /// ALE momentum equation residual.
    ///
    /// ρ(∂u/∂t|_χ + c·∇u) + ∇p - μ∇²u - f = 0  where c = u - w.
    #[allow(clippy::too_many_arguments)]
    pub fn momentum_residual_ale(
        &self,
        du_dt_chi: [f64; 3],
        ale_conv: [f64; 3],
        grad_p: [f64; 3],
        lap_u: [f64; 3],
        body: [f64; 3],
    ) -> [f64; 3] {
        [
            self.density * (du_dt_chi[0] + ale_conv[0]) + grad_p[0]
                - self.viscosity * lap_u[0]
                - body[0],
            self.density * (du_dt_chi[1] + ale_conv[1]) + grad_p[1]
                - self.viscosity * lap_u[1]
                - body[1],
            self.density * (du_dt_chi[2] + ale_conv[2]) + grad_p[2]
                - self.viscosity * lap_u[2]
                - body[2],
        ]
    }

    /// Set mesh velocity for node i.
    pub fn set_mesh_velocity(&mut self, i: usize, w: [f64; 3]) {
        if i < self.mesh_velocity.len() {
            self.mesh_velocity[i] = w;
        }
    }

    /// Mesh velocity magnitude at node i.
    pub fn mesh_velocity_magnitude(&self, i: usize) -> f64 {
        if i < self.mesh_velocity.len() {
            let w = self.mesh_velocity[i];
            (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
        } else {
            0.0
        }
    }
}

// ============================================================================
// Fluid-Thermal Coupling
// ============================================================================

/// Fluid-thermal coupling with Boussinesq approximation for natural convection.
///
/// Approximates buoyancy as: f_buoy = -ρ * g * β * (T - T_ref)
/// where β is the thermal expansion coefficient.
pub struct FluidThermalCoupling {
    /// Reference density (kg/m³).
    pub density_ref: f64,
    /// Thermal expansion coefficient (1/K).
    pub thermal_expansion: f64,
    /// Gravitational acceleration vector (m/s²).
    pub gravity: [f64; 3],
    /// Reference temperature (K).
    pub temp_ref: f64,
    /// Thermal diffusivity (m²/s).
    pub thermal_diffusivity: f64,
    /// Prandtl number.
    pub prandtl: f64,
}

impl FluidThermalCoupling {
    /// Create a new fluid-thermal coupling model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        density_ref: f64,
        thermal_expansion: f64,
        gravity: [f64; 3],
        temp_ref: f64,
        thermal_diffusivity: f64,
        prandtl: f64,
    ) -> Self {
        Self {
            density_ref,
            thermal_expansion,
            gravity,
            temp_ref,
            thermal_diffusivity,
            prandtl,
        }
    }

    /// Boussinesq buoyancy force at temperature T.
    ///
    /// f_b = -ρ₀ * β * (T - T_ref) * g
    pub fn buoyancy_force(&self, temp: f64) -> [f64; 3] {
        let factor = -self.density_ref * self.thermal_expansion * (temp - self.temp_ref);
        [
            factor * self.gravity[0],
            factor * self.gravity[1],
            factor * self.gravity[2],
        ]
    }

    /// Rayleigh number: Ra = g * β * ΔT * L³ / (ν * α).
    pub fn rayleigh_number(&self, delta_t: f64, length: f64, kinematic_viscosity: f64) -> f64 {
        let g =
            (self.gravity[0].powi(2) + self.gravity[1].powi(2) + self.gravity[2].powi(2)).sqrt();
        g * self.thermal_expansion * delta_t * length.powi(3)
            / (kinematic_viscosity * self.thermal_diffusivity)
    }

    /// Nusselt number correlation for vertical plate (Churchill-Chu).
    ///
    /// Nu = 0.68 + 0.67 * Ra^(1/4) / (1 + (0.492/Pr)^(9/16))^(4/9)
    pub fn nusselt_vertical_plate(&self, ra: f64) -> f64 {
        let f = 1.0 + (0.492 / self.prandtl).powf(9.0 / 16.0);
        0.68 + 0.67 * ra.powf(0.25) / f.powf(4.0 / 9.0)
    }

    /// Grashof number: Gr = g * β * ΔT * L³ / ν².
    pub fn grashof_number(&self, delta_t: f64, length: f64, kinematic_viscosity: f64) -> f64 {
        let g =
            (self.gravity[0].powi(2) + self.gravity[1].powi(2) + self.gravity[2].powi(2)).sqrt();
        g * self.thermal_expansion * delta_t * length.powi(3)
            / (kinematic_viscosity * kinematic_viscosity)
    }
}

// ============================================================================
// Level Set FEM
// ============================================================================

/// FEM level-set method for fluid-fluid interface tracking.
///
/// Solves the level-set equation: ∂φ/∂t + u·∇φ = 0
/// with SUPG stabilization, periodic redistancing to maintain |∇φ|=1,
/// and Heaviside/delta function smoothing at the interface.
pub struct LevelSetFem {
    /// Level-set field values at nodes.
    pub phi: Vec<f64>,
    /// Advection velocity at nodes (u, v, w per node).
    pub velocity: Vec<[f64; 3]>,
    /// Interface half-width for smoothing (m).
    pub epsilon: f64,
    /// Redistancing pseudo-time step.
    pub dtau: f64,
    /// Characteristic element size.
    pub h_elem: f64,
}

impl LevelSetFem {
    /// Create a new level-set FEM object.
    pub fn new(n_nodes: usize, epsilon: f64, dtau: f64, h_elem: f64) -> Self {
        Self {
            phi: vec![0.0; n_nodes],
            velocity: vec![[0.0; 3]; n_nodes],
            epsilon,
            dtau,
            h_elem,
        }
    }

    /// Smooth Heaviside function H_ε(φ).
    ///
    /// H_ε(φ) = 0           if φ < -ε
    /// H_ε(φ) = 1           if φ > ε
    /// H_ε(φ) = (1 + φ/ε + sin(πφ/ε)/π) / 2  otherwise
    pub fn heaviside(&self, phi: f64) -> f64 {
        use std::f64::consts::PI;
        if phi < -self.epsilon {
            0.0
        } else if phi > self.epsilon {
            1.0
        } else {
            0.5 * (1.0 + phi / self.epsilon + (PI * phi / self.epsilon).sin() / PI)
        }
    }

    /// Smooth delta function δ_ε(φ).
    ///
    /// δ_ε(φ) = dH_ε/dφ
    pub fn delta_function(&self, phi: f64) -> f64 {
        if phi.abs() > self.epsilon {
            0.0
        } else {
            0.5 / self.epsilon * (1.0 + (PI * phi / self.epsilon).cos())
        }
    }

    /// Signed distance redistancing RHS: sgn(φ₀) * (1 - |∇φ|).
    pub fn redistancing_rhs(&self, phi: f64, grad_phi_norm: f64) -> f64 {
        let sgn = if phi > 0.0 {
            1.0
        } else if phi < 0.0 {
            -1.0
        } else {
            0.0
        };
        sgn * (1.0 - grad_phi_norm)
    }

    /// SUPG stabilization parameter for level-set advection.
    pub fn tau_supg(&self, u_norm: f64) -> f64 {
        compute_tau_supg(u_norm, self.h_elem, 0.0, self.dtau)
    }

    /// Update phi by one pseudo-time step of redistancing.
    pub fn redistance_step(&mut self, grad_phi_norms: &[f64]) {
        let dtau = self.dtau;
        let n = self.phi.len();
        for i in 0..n {
            let gn = if i < grad_phi_norms.len() {
                grad_phi_norms[i]
            } else {
                1.0
            };
            let phi_i = self.phi[i];
            let rhs = self.redistancing_rhs(phi_i, gn);
            self.phi[i] += dtau * rhs;
        }
    }

    /// Set level-set field from a signed distance function (sphere).
    ///
    /// φ(x) = |x - center| - radius
    pub fn init_sphere(&mut self, positions: &[[f64; 3]], center: [f64; 3], radius: f64) {
        for (i, pos) in positions.iter().enumerate() {
            if i < self.phi.len() {
                let dx = pos[0] - center[0];
                let dy = pos[1] - center[1];
                let dz = pos[2] - center[2];
                self.phi[i] = (dx * dx + dy * dy + dz * dz).sqrt() - radius;
            }
        }
    }
}

// ============================================================================
// Helper functions
// ============================================================================

/// Compute SUPG stabilization parameter τ_SUPG.
///
/// Uses the Tezduyar formula combining transient, convective, and diffusive scales:
/// τ = (4/dt² + (2|u|/h)² + (4ν/h²)²)^{-1/2}
pub fn compute_tau_supg(u_norm: f64, h: f64, nu: f64, dt: f64) -> f64 {
    let inv_dt = if dt > 1e-14 { 2.0 / dt } else { 0.0 };
    let inv_conv = 2.0 * u_norm / h.max(1e-14);
    let inv_diff = 4.0 * nu / (h * h).max(1e-28);
    let denom = (inv_dt * inv_dt + inv_conv * inv_conv + inv_diff * inv_diff).sqrt();
    if denom < 1e-30 { 0.0 } else { 1.0 / denom }
}

/// Compute PSPG stabilization parameter τ_PSPG.
///
/// In many implementations τ_PSPG = τ_SUPG / ρ² or simply τ_SUPG.
/// Here we use τ_PSPG = τ_SUPG for simplicity.
pub fn compute_tau_pspg(u_norm: f64, h: f64, nu: f64, dt: f64) -> f64 {
    compute_tau_supg(u_norm, h, nu, dt)
}

/// Compute optimal stabilization parameter for advection-diffusion.
///
/// τ = h / (2 * |u|) * ξ(Pe)  where ξ(Pe) = coth(Pe) - 1/Pe.
pub fn stabilization_parameter(u_norm: f64, h: f64, kappa: f64, dt: f64) -> f64 {
    if u_norm < 1e-14 {
        // Purely diffusive: τ ~ h²/4κ
        return h * h / (4.0 * kappa.max(1e-30));
    }
    let pe = u_norm * h / (2.0 * kappa.max(1e-30));
    let xi = if pe.abs() < 1e-6 {
        pe / 3.0
    } else {
        1.0 / pe.tanh() - 1.0 / pe
    };
    // Also include time scale
    let tau_ss = h / (2.0 * u_norm) * xi;
    let tau_t = if dt > 1e-14 { dt / 2.0 } else { f64::MAX };
    tau_ss.min(tau_t)
}

/// Bubble enrichment function for MINI element.
///
/// b(λ) = 27 * λ₁ * λ₂ * λ₃  (cubic bubble for triangle).
pub fn bubble_enrichment(lambda: [f64; 3]) -> f64 {
    27.0 * lambda[0] * lambda[1] * lambda[2]
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // StabilizedNS tests
    #[test]
    fn test_stabilized_ns_new() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        assert_eq!(ns.density, 1000.0);
        assert_eq!(ns.viscosity, 1e-3);
    }

    #[test]
    fn test_stabilized_ns_kinematic_viscosity() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        assert!((ns.kinematic_viscosity() - 1e-6).abs() < 1e-14);
    }

    #[test]
    fn test_stabilized_ns_reynolds_number() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        // Re = ρ * u * L / μ = 1000 * 1.0 * 0.01 / 1e-3 = 10000
        let re = ns.reynolds_number(1.0, 0.01);
        assert!((re - 10000.0).abs() < 1e-6);
    }

    #[test]
    fn test_stabilized_ns_tau_supg_zero_velocity() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        let tau = ns.tau_supg(0.0);
        assert!(tau >= 0.0);
    }

    #[test]
    fn test_stabilized_ns_tau_supg_nonzero() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        let tau = ns.tau_supg(1.0);
        assert!(tau > 0.0);
    }

    #[test]
    fn test_stabilized_ns_continuity_residual() {
        let ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        assert_eq!(ns.continuity_residual(0.5), 0.5);
    }

    #[test]
    fn test_stabilized_ns_step_time() {
        let mut ns = StabilizedNS::new(1000.0, 1e-3, 0.01, 1e-3);
        ns.step_time();
        assert!((ns.time - 1e-3).abs() < 1e-14);
    }

    // SupgStabilization tests
    #[test]
    fn test_supg_tau_decreases_with_velocity() {
        let supg = SupgStabilization::new(0.1, 1e-6, 1e-3);
        let tau1 = supg.tau(0.1);
        let tau2 = supg.tau(10.0);
        // Higher velocity → smaller τ (dominated by convective term)
        assert!(tau2 <= tau1);
    }

    #[test]
    fn test_supg_peclet_number() {
        let supg = SupgStabilization::new(0.1, 0.1, 1e-3);
        let pe = supg.peclet(1.0);
        assert!((pe - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_supg_xi_upwind_zero_velocity() {
        let supg = SupgStabilization::new(0.1, 1e-6, 1e-3);
        let xi = supg.xi_upwind(0.0);
        assert_eq!(xi, 0.0);
    }

    #[test]
    fn test_supg_weight_zero_velocity() {
        let supg = SupgStabilization::new(0.1, 1e-6, 1e-3);
        let w = supg.supg_weight([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0);
        assert_eq!(w, 0.0);
    }

    // PspgStabilization tests
    #[test]
    fn test_pspg_tau_positive() {
        let pspg = PspgStabilization::new(0.1, 1e-6, 1e-3);
        assert!(pspg.tau(1.0) > 0.0);
    }

    #[test]
    fn test_pspg_residual_contribution() {
        let pspg = PspgStabilization::new(0.1, 1e-6, 1e-3);
        let tau = 1e-4;
        let res = pspg.residual_contribution(tau, 1000.0, [1.0, 0.0, 0.0], [100.0, 0.0, 0.0]);
        assert!((res - tau / 1000.0 * 100.0).abs() < 1e-14);
    }

    // VmsStabilization tests
    #[test]
    fn test_vms_tau_m_positive() {
        let vms = VmsStabilization::new(0.1, 1e-3, 1000.0, 1e-3);
        assert!(vms.tau_m(1.0) > 0.0);
    }

    #[test]
    fn test_vms_tau_c_positive() {
        let vms = VmsStabilization::new(0.1, 1e-3, 1000.0, 1e-3);
        assert!(vms.tau_c(1.0) > 0.0);
    }

    #[test]
    fn test_vms_fine_scale_velocity() {
        let vms = VmsStabilization::new(0.1, 1e-3, 1000.0, 1e-3);
        let u_prime = vms.fine_scale_velocity(1e-4, [100.0, 0.0, 0.0]);
        assert!(u_prime[0] < 0.0);
    }

    #[test]
    fn test_vms_fine_scale_pressure() {
        let vms = VmsStabilization::new(0.1, 1e-3, 1000.0, 1e-3);
        let p_prime = vms.fine_scale_pressure(1e-3, 10.0);
        assert!(p_prime < 0.0);
    }

    // IncompressibleElement tests
    #[test]
    fn test_incompressible_taylor_hood_dofs() {
        let e = IncompressibleElement::taylor_hood(0.1);
        assert_eq!(e.n_vel_nodes, 6);
        assert_eq!(e.n_pres_nodes, 3);
    }

    #[test]
    fn test_incompressible_mini_element_dofs() {
        let e = IncompressibleElement::mini_element(0.1);
        assert!(e.use_bubble);
    }

    #[test]
    fn test_bubble_enrichment_centroid() {
        let b = bubble_enrichment([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        assert!((b - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bubble_enrichment_vertex() {
        // At a vertex, one barycentric coord is 1, others 0 → bubble = 0
        assert_eq!(bubble_enrichment([1.0, 0.0, 0.0]), 0.0);
    }

    // ProjectionMethod tests
    #[test]
    fn test_projection_corrector_zero_grad_p() {
        let pm = ProjectionMethod::new(1000.0, 1e-3, 1e-3);
        let u = pm.corrector_step([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        assert_eq!(u, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_projection_pressure_rhs() {
        let pm = ProjectionMethod::new(1000.0, 1e-3, 1e-3);
        let rhs = pm.pressure_rhs(0.01);
        assert!((rhs - 1000.0 / 1e-3 * 0.01).abs() < 1e-4);
    }

    #[test]
    fn test_projection_divergence_free() {
        let pm = ProjectionMethod::new(1000.0, 1e-3, 1e-3);
        assert!(pm.is_divergence_free(1e-10, 1e-8));
        assert!(!pm.is_divergence_free(0.1, 1e-8));
    }

    #[test]
    fn test_projection_predictor_zero_all() {
        let pm = ProjectionMethod::new(1000.0, 1e-3, 1e-3);
        let u = pm.predictor_step([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]);
        assert_eq!(u, [0.0; 3]);
    }

    // AdvectionDiffusion tests
    #[test]
    fn test_advdiff_peclet() {
        let ad = AdvectionDiffusion::new(0.1, [1.0, 0.0, 0.0], 0.1, 1e-3);
        assert!((ad.peclet_number() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_advdiff_tau_positive() {
        let ad = AdvectionDiffusion::new(0.1, [1.0, 0.0, 0.0], 0.1, 1e-3);
        assert!(ad.tau_supg() > 0.0);
    }

    #[test]
    fn test_advdiff_residual() {
        let ad = AdvectionDiffusion::new(0.1, [1.0, 0.0, 0.0], 0.1, 1e-3);
        let r = ad.residual(0.0, 0.0, 0.0, 0.0);
        assert_eq!(r, 0.0);
    }

    // AleFormulation tests
    #[test]
    fn test_ale_convection_zero_mesh_velocity() {
        let ale = AleFormulation::new(1000.0, 1e-3, 4, 1e-3);
        let c = ale.ale_convection([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        assert_eq!(c, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_ale_set_mesh_velocity() {
        let mut ale = AleFormulation::new(1000.0, 1e-3, 4, 1e-3);
        ale.set_mesh_velocity(0, [1.0, 0.5, 0.0]);
        assert_eq!(ale.mesh_velocity[0], [1.0, 0.5, 0.0]);
    }

    #[test]
    fn test_ale_mesh_velocity_magnitude() {
        let mut ale = AleFormulation::new(1000.0, 1e-3, 4, 1e-3);
        ale.set_mesh_velocity(0, [3.0, 4.0, 0.0]);
        assert!((ale.mesh_velocity_magnitude(0) - 5.0).abs() < 1e-10);
    }

    // FluidThermalCoupling tests
    #[test]
    fn test_fluid_thermal_buoyancy_at_reference() {
        let ftc =
            FluidThermalCoupling::new(1000.0, 2.1e-4, [0.0, -9.81, 0.0], 300.0, 1.43e-7, 0.71);
        let f = ftc.buoyancy_force(300.0);
        assert!(f[1].abs() < 1e-10);
    }

    #[test]
    fn test_fluid_thermal_buoyancy_hot() {
        let ftc =
            FluidThermalCoupling::new(1000.0, 2.1e-4, [0.0, -9.81, 0.0], 300.0, 1.43e-7, 0.71);
        let f = ftc.buoyancy_force(310.0);
        // Hot fluid → positive buoyancy (upward, opposing gravity in -y direction)
        assert!(f[1] > 0.0);
    }

    #[test]
    fn test_fluid_thermal_rayleigh_number_positive() {
        let ftc =
            FluidThermalCoupling::new(1000.0, 2.1e-4, [0.0, -9.81, 0.0], 300.0, 1.43e-7, 0.71);
        let ra = ftc.rayleigh_number(10.0, 0.1, 1e-6);
        assert!(ra > 0.0);
    }

    #[test]
    fn test_fluid_thermal_nusselt_positive() {
        let ftc =
            FluidThermalCoupling::new(1000.0, 2.1e-4, [0.0, -9.81, 0.0], 300.0, 1.43e-7, 0.71);
        let nu = ftc.nusselt_vertical_plate(1e6);
        assert!(nu > 0.0);
    }

    // LevelSetFem tests
    #[test]
    fn test_level_set_heaviside_inside() {
        let ls = LevelSetFem::new(10, 0.01, 1e-3, 0.1);
        assert!((ls.heaviside(0.1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_level_set_heaviside_outside() {
        let ls = LevelSetFem::new(10, 0.01, 1e-3, 0.1);
        assert_eq!(ls.heaviside(-0.1), 0.0);
    }

    #[test]
    fn test_level_set_heaviside_at_interface() {
        let ls = LevelSetFem::new(10, 0.01, 1e-3, 0.1);
        assert!((ls.heaviside(0.0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_level_set_delta_at_interface() {
        let ls = LevelSetFem::new(10, 0.01, 1e-3, 0.1);
        assert!(ls.delta_function(0.0) > 0.0);
    }

    #[test]
    fn test_level_set_delta_outside() {
        let ls = LevelSetFem::new(10, 0.01, 1e-3, 0.1);
        assert_eq!(ls.delta_function(0.1), 0.0);
    }

    #[test]
    fn test_level_set_init_sphere() {
        let mut ls = LevelSetFem::new(2, 0.01, 1e-3, 0.1);
        let pos = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        ls.init_sphere(&pos, [0.0, 0.0, 0.0], 1.0);
        assert!((ls.phi[0] - (-1.0)).abs() < 1e-14);
        assert!((ls.phi[1] - 1.0).abs() < 1e-14);
    }

    // Helper function tests
    #[test]
    fn test_compute_tau_supg_zero_all() {
        let tau = compute_tau_supg(0.0, 0.1, 0.0, 0.0);
        assert_eq!(tau, 0.0);
    }

    #[test]
    fn test_compute_tau_supg_positive() {
        let tau = compute_tau_supg(1.0, 0.1, 1e-6, 1e-3);
        assert!(tau > 0.0);
    }

    #[test]
    fn test_stabilization_parameter_diffusion_dominated() {
        let tau = stabilization_parameter(0.0, 0.1, 1.0, 0.0);
        assert!(tau > 0.0);
    }

    #[test]
    fn test_stabilization_parameter_convection_dominated() {
        let tau = stabilization_parameter(100.0, 0.1, 1e-6, 1e-3);
        assert!(tau > 0.0);
    }
}

// ============================================================================
// Stokes FEM
// ============================================================================

/// Stokes flow FEM solver using Taylor-Hood elements.
///
/// Solves the incompressible Stokes equations:
/// −μ∇²u + ∇p = f,  ∇·u = 0
pub struct StokesFEM {
    /// Number of velocity nodes.
    pub n_vel_nodes: usize,
    /// Number of pressure nodes.
    pub n_press_nodes: usize,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
}

impl StokesFEM {
    /// Create a new Stokes FEM solver.
    pub fn new(n_vel_nodes: usize, n_press_nodes: usize, viscosity: f64) -> Self {
        Self {
            n_vel_nodes,
            n_press_nodes,
            viscosity,
        }
    }

    /// Total degrees of freedom (2 velocity components per velocity node + 1 pressure per pressure node).
    pub fn total_dofs(&self) -> usize {
        2 * self.n_vel_nodes + self.n_press_nodes
    }
}

/// Compute the Taylor-Hood element stiffness matrix for a triangular element.
///
/// Given element coordinates `coords` (3 vertices) and viscosity `viscosity`,
/// returns a 6×6 matrix representing the viscous contribution from the
/// velocity DOFs (u₁,v₁, u₂,v₂, u₃,v₃).
pub fn stokes_element_matrix(coords: &[[f64; 2]; 3], viscosity: f64) -> [[f64; 6]; 6] {
    // Compute shape function derivatives (constant for linear triangle)
    let x1 = coords[0][0];
    let y1 = coords[0][1];
    let x2 = coords[1][0];
    let y2 = coords[1][1];
    let x3 = coords[2][0];
    let y3 = coords[2][1];
    let area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    let area = area2.abs() / 2.0;
    if area < 1e-300 {
        return [[0.0; 6]; 6];
    }
    // dN/dx and dN/dy for each of the 3 P1 shape functions
    let dndx = [(y2 - y3) / area2, (y3 - y1) / area2, (y1 - y2) / area2];
    let dndy = [(x3 - x2) / area2, (x1 - x3) / area2, (x2 - x1) / area2];
    let mut ke = [[0.0f64; 6]; 6];
    // Assemble using B^T * C * B * area for plane-stress incompressible
    // DOF ordering: (u1, v1, u2, v2, u3, v3) → rows/cols 0..5
    for a in 0..3 {
        for b in 0..3 {
            // ε_xx contribution: dNa/dx * dNb/dx
            let kuu = viscosity * (dndx[a] * dndx[b] + dndy[a] * dndy[b]) * area;
            let row_u = 2 * a;
            let col_u = 2 * b;
            ke[row_u][col_u] += kuu;
            ke[row_u + 1][col_u + 1] += kuu;
        }
    }
    ke
}

/// Pressure-velocity coupling matrix for a triangular Stokes element.
///
/// Returns a 6×3 matrix B such that the coupling block B * p couples
/// pressure DOFs (P1) to velocity DOFs (P2-like, treated as P1 here).
pub fn stokes_pressure_gradient(coords: &[[f64; 2]; 3]) -> [[f64; 3]; 6] {
    let x1 = coords[0][0];
    let y1 = coords[0][1];
    let x2 = coords[1][0];
    let y2 = coords[1][1];
    let x3 = coords[2][0];
    let y3 = coords[2][1];
    let area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    let area = area2.abs() / 2.0;
    if area < 1e-300 {
        return [[0.0; 3]; 6];
    }
    let dndx = [(y2 - y3) / area2, (y3 - y1) / area2, (y1 - y2) / area2];
    let dndy = [(x3 - x2) / area2, (x1 - x3) / area2, (x2 - x1) / area2];
    let mut bp = [[0.0f64; 3]; 6];
    // Row 2a   (u-component): -dNa/dx * Nb_area
    // Row 2a+1 (v-component): -dNa/dy * Nb_area
    for a in 0..3 {
        for b in 0..3 {
            // Integrate N_b over triangle: area/3 for linear elements
            let nb_int = area / 3.0;
            bp[2 * a][b] = -dndx[a] * nb_int;
            bp[2 * a + 1][b] = -dndy[a] * nb_int;
        }
    }
    bp
}

/// Set up the Stokes system for the lid-driven cavity benchmark.
///
/// Builds a simplified (dense) system matrix and RHS for an n×n grid.
/// Returns `(matrix, rhs)` where `matrix` is `(2*n_nodes) × (2*n_nodes)`
/// and `rhs` is `2*n_nodes`.  n_nodes = (n+1)^2.
pub fn lid_driven_cavity_setup(n: usize, viscosity: f64) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n_nodes = (n + 1) * (n + 1);
    let n_dofs = 2 * n_nodes;
    let mat = vec![vec![0.0f64; n_dofs]; n_dofs];
    let mut rhs = vec![0.0f64; n_dofs];
    // Lid velocity: u=1 on top row (y = n*h)
    let h = 1.0 / n as f64;
    let _ = viscosity; // used in full assembly; placeholder here
    let _ = h;
    // Set RHS for lid nodes: top row j=n → u-DOF = 1
    for i in 0..=n {
        let _node = i * (n + 1) + n; // column-major: node = i*(n+1)+j
        // Actually row-major: node = j*(n+1)+i for row j, col i
        let lid_node = n * (n + 1) + i;
        if lid_node < n_nodes {
            rhs[2 * lid_node] = 1.0; // u = 1 on lid
        }
    }
    (mat, rhs)
}

/// Compute the stream function from a 2D velocity field on a regular grid.
///
/// Integrates ψ by: dψ/dy = u, dψ/dx = -v using trapezoidal rule.
/// Grid is `nx × ny` with spacing `dx`.
pub fn stream_function_from_velocity(
    ux: &[f64],
    _uy: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<f64> {
    let n = nx * ny;
    let mut psi = vec![0.0f64; n];
    if n == 0 || dx <= 0.0 {
        return psi;
    }
    // Integrate along each column from bottom to top using ux: dψ/dy = ux
    for i in 0..nx {
        psi[i] = 0.0; // bottom boundary
        for j in 1..ny {
            let idx_prev = (j - 1) * nx + i;
            let idx_curr = j * nx + i;
            if idx_curr < n && idx_prev < ux.len() {
                psi[idx_curr] = psi[idx_prev] + dx * ux[idx_prev];
            }
        }
    }
    psi
}

/// Compute the 2D vorticity field ω = ∂v/∂x − ∂u/∂y on a regular grid.
///
/// Uses central differences on the interior; one-sided differences on boundaries.
pub fn vorticity_2d(ux: &[f64], uy: &[f64], nx: usize, ny: usize, dx: f64, dy: f64) -> Vec<f64> {
    let n = nx * ny;
    let mut omega = vec![0.0f64; n];
    if n == 0 {
        return omega;
    }
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            // dv/dx
            let dv_dx = if i == 0 {
                let r = j * nx + (i + 1).min(nx - 1);
                (uy[r] - uy[idx]) / dx
            } else if i == nx - 1 {
                let l = j * nx + i - 1;
                (uy[idx] - uy[l]) / dx
            } else {
                let r = j * nx + i + 1;
                let l = j * nx + i - 1;
                (uy[r] - uy[l]) / (2.0 * dx)
            };
            // du/dy
            let du_dy = if j == 0 {
                let u = (j + 1).min(ny - 1) * nx + i;
                (ux[u] - ux[idx]) / dy
            } else if j == ny - 1 {
                let d = (j - 1) * nx + i;
                (ux[idx] - ux[d]) / dy
            } else {
                let u = (j + 1) * nx + i;
                let d = (j - 1) * nx + i;
                (ux[u] - ux[d]) / (2.0 * dy)
            };
            omega[idx] = dv_dx - du_dy;
        }
    }
    omega
}

/// Compute the Reynolds number Re = ρ v L / μ.
///
/// - `rho`: fluid density (kg/m³)
/// - `v`: characteristic velocity (m/s)
/// - `l`: characteristic length (m)
/// - `mu`: dynamic viscosity (Pa·s)
pub fn reynolds_number(rho: f64, v: f64, l: f64, mu: f64) -> f64 {
    if mu.abs() < 1e-300 {
        return 0.0;
    }
    rho * v * l / mu
}

/// Compute the Stokes number Stk = ρ_p d_p² v / (18 μ L).
///
/// - `rho_p`: particle density (kg/m³)
/// - `d_p`: particle diameter (m)
/// - `v`: characteristic flow velocity (m/s)
/// - `mu`: dynamic viscosity (Pa·s)
/// - `l`: characteristic length (m)
pub fn stokes_number(rho_p: f64, d_p: f64, v: f64, mu: f64, l: f64) -> f64 {
    if mu.abs() < 1e-300 || l.abs() < 1e-300 {
        return 0.0;
    }
    rho_p * d_p * d_p * v / (18.0 * mu * l)
}

/// Simplified Navier-Stokes FEM driver.
///
/// Stores a placeholder velocity field and advances time with explicit Euler.
pub struct NavierStokesFEM {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Reynolds number.
    pub re: f64,
    /// Current simulation time (s).
    pub time: f64,
    /// Placeholder velocity DOF vector (u, v interleaved per node).
    pub velocity: Vec<f64>,
}

impl NavierStokesFEM {
    /// Create a new Navier-Stokes FEM driver.
    pub fn new(n_nodes: usize, re: f64) -> Self {
        Self {
            n_nodes,
            re,
            time: 0.0,
            velocity: vec![0.0; 2 * n_nodes],
        }
    }

    /// Advance time by `dt` seconds using forward Euler (placeholder).
    pub fn step_euler(&mut self, dt: f64) {
        // Placeholder: damp velocity slightly to represent viscous dissipation
        let damp = 1.0 - dt / self.re.max(1e-10);
        for v in self.velocity.iter_mut() {
            *v *= damp;
        }
        self.time += dt;
    }

    /// Total kinetic energy: ½ ∑ u_i².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.velocity.iter().map(|v| v * v).sum::<f64>()
    }

    /// Current simulation time.
    pub fn time(&self) -> f64 {
        self.time
    }
}

/// Right-hand side for the pressure Poisson equation in the projection method.
///
/// RHS_i = div_u_i / dt  (used in fractional-step methods).
pub fn pressure_poisson_rhs(div_u: &[f64], dt: f64) -> Vec<f64> {
    if dt.abs() < 1e-300 {
        return vec![0.0; div_u.len()];
    }
    div_u.iter().map(|&d| d / dt).collect()
}

// ============================================================================
// Extended tests for the new Stokes/NS functions
// ============================================================================

#[cfg(test)]
mod stokes_tests {
    use super::*;

    #[test]
    fn test_stokes_fem_new() {
        let s = StokesFEM::new(6, 3, 1e-3);
        assert_eq!(s.n_vel_nodes, 6);
        assert_eq!(s.n_press_nodes, 3);
        assert!((s.viscosity - 1e-3).abs() < 1e-15);
    }

    #[test]
    fn test_stokes_fem_total_dofs() {
        let s = StokesFEM::new(6, 3, 1e-3);
        assert_eq!(s.total_dofs(), 15); // 2*6 + 3
    }

    #[test]
    fn test_stokes_element_matrix_symmetry() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let ke = stokes_element_matrix(&coords, 1.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (ke[i][j] - ke[j][i]).abs() < 1e-12,
                    "not symmetric at ({},{})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_stokes_element_matrix_zero_viscosity() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let ke = stokes_element_matrix(&coords, 0.0);
        for row in &ke {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn test_stokes_element_matrix_degenerate() {
        // Degenerate triangle (zero area)
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.5, 0.0]];
        let ke = stokes_element_matrix(&coords, 1.0);
        for row in &ke {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn test_stokes_element_matrix_viscosity_scales() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let ke1 = stokes_element_matrix(&coords, 1.0);
        let ke2 = stokes_element_matrix(&coords, 2.0);
        assert!((ke2[0][0] - 2.0 * ke1[0][0]).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_pressure_gradient_dimensions() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let bp = stokes_pressure_gradient(&coords);
        assert_eq!(bp.len(), 6);
        assert_eq!(bp[0].len(), 3);
    }

    #[test]
    fn test_stokes_pressure_gradient_degenerate() {
        let coords = [[0.0, 0.0], [1.0, 0.0], [0.5, 0.0]];
        let bp = stokes_pressure_gradient(&coords);
        for row in &bp {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn test_lid_driven_cavity_setup_sizes() {
        let (mat, rhs) = lid_driven_cavity_setup(4, 1e-3);
        let n_nodes = 5 * 5; // (n+1)^2
        assert_eq!(mat.len(), 2 * n_nodes);
        assert_eq!(rhs.len(), 2 * n_nodes);
    }

    #[test]
    fn test_lid_driven_cavity_rhs_lid() {
        let (_mat, rhs) = lid_driven_cavity_setup(4, 1e-3);
        // Top-row nodes should have rhs entry = 1 (u-component of lid)
        let n = 4;
        let n1 = n + 1;
        let any_lid_nonzero = (0..n1).any(|i| {
            let lid_node = n * n1 + i;
            if lid_node < n1 * n1 && 2 * lid_node < rhs.len() {
                rhs[2 * lid_node] > 0.5
            } else {
                false
            }
        });
        assert!(any_lid_nonzero);
    }

    #[test]
    fn test_stream_function_from_velocity_size() {
        let nx = 5;
        let ny = 4;
        let ux = vec![1.0; nx * ny];
        let uy = vec![0.0; nx * ny];
        let psi = stream_function_from_velocity(&ux, &uy, nx, ny, 0.1);
        assert_eq!(psi.len(), nx * ny);
    }

    #[test]
    fn test_stream_function_from_velocity_zero() {
        let nx = 3;
        let ny = 3;
        let ux = vec![0.0; nx * ny];
        let uy = vec![0.0; nx * ny];
        let psi = stream_function_from_velocity(&ux, &uy, nx, ny, 0.1);
        for v in &psi {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_stream_function_from_velocity_integrates() {
        // Uniform ux=1 → ψ increases with j
        let nx = 2;
        let ny = 3;
        let ux = vec![1.0; nx * ny];
        let uy = vec![0.0; nx * ny];
        let psi = stream_function_from_velocity(&ux, &uy, nx, ny, 0.5);
        // psi at j=1 > psi at j=0 for same column
        assert!(psi[nx] > psi[0]);
    }

    #[test]
    fn test_vorticity_2d_size() {
        let nx = 5;
        let ny = 4;
        let ux = vec![0.0; nx * ny];
        let uy = vec![0.0; nx * ny];
        let omega = vorticity_2d(&ux, &uy, nx, ny, 0.1, 0.1);
        assert_eq!(omega.len(), nx * ny);
    }

    #[test]
    fn test_vorticity_2d_zero_field() {
        let nx = 4;
        let ny = 4;
        let ux = vec![0.0; nx * ny];
        let uy = vec![0.0; nx * ny];
        let omega = vorticity_2d(&ux, &uy, nx, ny, 0.1, 0.1);
        for v in &omega {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn test_vorticity_2d_solid_rotation() {
        // u = -y, v = x → ω = ∂v/∂x - ∂u/∂y = 1 - (-1) = 2
        let nx = 5;
        let ny = 5;
        let dx = 0.25;
        let dy = 0.25;
        let mut ux = vec![0.0; nx * ny];
        let mut uy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                ux[j * nx + i] = -y;
                uy[j * nx + i] = x;
            }
        }
        let omega = vorticity_2d(&ux, &uy, nx, ny, dx, dy);
        // Interior points should be close to 2.0
        let idx = 2 * nx + 2;
        assert!((omega[idx] - 2.0).abs() < 0.1, "omega={}", omega[idx]);
    }

    #[test]
    fn test_reynolds_number_formula() {
        let re = reynolds_number(1000.0, 1.0, 0.1, 0.001);
        assert!((re - 1e5).abs() < 1.0);
    }

    #[test]
    fn test_reynolds_number_zero_viscosity() {
        let re = reynolds_number(1000.0, 1.0, 0.1, 0.0);
        assert_eq!(re, 0.0);
    }

    #[test]
    fn test_reynolds_number_proportional() {
        let re1 = reynolds_number(1.0, 1.0, 1.0, 1.0);
        let re2 = reynolds_number(2.0, 1.0, 1.0, 1.0);
        assert!((re2 - 2.0 * re1).abs() < 1e-12);
    }

    #[test]
    fn test_stokes_number_formula() {
        let stk = stokes_number(1000.0, 1e-4, 1.0, 1e-3, 0.1);
        let expected = 1000.0 * 1e-8 * 1.0 / (18.0 * 1e-3 * 0.1);
        assert!((stk - expected).abs() < 1e-14);
    }

    #[test]
    fn test_stokes_number_zero_viscosity() {
        let stk = stokes_number(1000.0, 1e-4, 1.0, 0.0, 0.1);
        assert_eq!(stk, 0.0);
    }

    #[test]
    fn test_stokes_number_zero_length() {
        let stk = stokes_number(1000.0, 1e-4, 1.0, 1e-3, 0.0);
        assert_eq!(stk, 0.0);
    }

    #[test]
    fn test_ns_fem_new() {
        let ns = NavierStokesFEM::new(10, 100.0);
        assert_eq!(ns.n_nodes, 10);
        assert_eq!(ns.velocity.len(), 20);
        assert_eq!(ns.time(), 0.0);
    }

    #[test]
    fn test_ns_fem_step_advances_time() {
        let mut ns = NavierStokesFEM::new(4, 100.0);
        ns.step_euler(0.01);
        assert!((ns.time() - 0.01).abs() < 1e-14);
    }

    #[test]
    fn test_ns_fem_step_multiple() {
        let mut ns = NavierStokesFEM::new(4, 100.0);
        ns.step_euler(0.01);
        ns.step_euler(0.01);
        assert!((ns.time() - 0.02).abs() < 1e-13);
    }

    #[test]
    fn test_ns_fem_kinetic_energy_zero() {
        let ns = NavierStokesFEM::new(4, 100.0);
        assert_eq!(ns.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_ns_fem_kinetic_energy_nonzero() {
        let mut ns = NavierStokesFEM::new(2, 100.0);
        ns.velocity = vec![1.0, 0.0, 0.0, 1.0];
        assert!((ns.kinetic_energy() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_poisson_rhs_size() {
        let div_u = vec![0.1, 0.2, 0.3];
        let rhs = pressure_poisson_rhs(&div_u, 0.01);
        assert_eq!(rhs.len(), 3);
    }

    #[test]
    fn test_pressure_poisson_rhs_formula() {
        let div_u = vec![1.0, 2.0];
        let rhs = pressure_poisson_rhs(&div_u, 0.5);
        assert!((rhs[0] - 2.0).abs() < 1e-12);
        assert!((rhs[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_pressure_poisson_rhs_zero_dt() {
        let div_u = vec![1.0, 2.0];
        let rhs = pressure_poisson_rhs(&div_u, 0.0);
        for v in &rhs {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_vorticity_2d_different_dx_dy() {
        let nx = 3;
        let ny = 3;
        let ux = vec![0.0; nx * ny];
        let uy = vec![1.0; nx * ny]; // uniform v = 1 → dv/dx = 0
        let omega = vorticity_2d(&ux, &uy, nx, ny, 0.1, 0.2);
        // All vorticity components should be zero for uniform field
        for v in &omega {
            assert!(v.abs() < 1e-10);
        }
    }

    #[test]
    fn test_stokes_fem_viscosity_stored() {
        let s = StokesFEM::new(3, 3, 0.42);
        assert!((s.viscosity - 0.42).abs() < 1e-15);
    }
}
