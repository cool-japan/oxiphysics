//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{CS2, CX, CY, NQ, W};

/// Hele-Shaw flow parameters in a thin-gap device.
#[derive(Debug, Clone)]
pub struct HeleShawParams {
    /// Gap height (m).
    pub gap_height: f64,
    /// Lateral width (m).
    pub width: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
}
impl HeleShawParams {
    /// Compute Hele-Shaw permeability K = h² / 12.
    pub fn permeability(&self) -> f64 {
        self.gap_height * self.gap_height / 12.0
    }
    /// Darcy velocity u = −K/μ · ∇P.
    pub fn darcy_velocity(&self, pressure_gradient: f64) -> f64 {
        -self.permeability() / self.viscosity * pressure_gradient
    }
    /// Flow resistance per unit length R = 12 μ / h².
    pub fn flow_resistance(&self) -> f64 {
        12.0 * self.viscosity / (self.gap_height * self.gap_height)
    }
}
/// Thrombus growth model coupling platelet activation and coagulation cascade.
///
/// Thrombus volume fraction increases with local fibrin deposition and adhered platelets.
pub struct ThrombosisModel {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Thrombus volume fraction per cell ∈ \[0, 1\].
    pub thrombus_fraction: Vec<f64>,
    /// Per-cell coagulation state.
    pub coag: Vec<CoagulationState>,
    /// Platelet collection.
    pub platelets: Vec<Platelet>,
    /// Platelet activation rate constant (1/(Pa·s)).
    pub k_act: f64,
    /// Shear activation threshold (Pa).
    pub tau_thresh: f64,
    /// Thrombus growth rate per fibrin (1/nM·s).
    pub growth_rate: f64,
}
impl ThrombosisModel {
    /// Construct a new thrombosis model on a `nx` × `ny` grid with `n_platelets` platelets.
    pub fn new(nx: usize, ny: usize, n_platelets: usize, k_act: f64, tau_thresh: f64) -> Self {
        let n = nx * ny;
        let platelets: Vec<Platelet> = (0..n_platelets)
            .map(|i| {
                let x = (i % nx) as f64 + 0.5;
                let y = (i / nx % ny) as f64 + 0.5;
                Platelet::new([x, y])
            })
            .collect();
        Self {
            nx,
            ny,
            thrombus_fraction: vec![0.0; n],
            coag: vec![CoagulationState::physiological(); n],
            platelets,
            k_act,
            tau_thresh,
            growth_rate: 1e-5,
        }
    }
    /// Advance thrombosis model by one step with given local shear stress field.
    pub fn step(&mut self, shear_stress: &[f64], dt: f64) {
        let nx = self.nx;
        for p in self.platelets.iter_mut() {
            let ix = clamp(p.pos[0], 0.0, (nx - 1) as f64) as usize;
            let iy = clamp(p.pos[1], 0.0, (self.ny - 1) as f64) as usize;
            let idx = ix + iy * nx;
            let tau = if idx < shear_stress.len() {
                shear_stress[idx]
            } else {
                0.0
            };
            p.update_activation(tau, dt, self.k_act, self.tau_thresh);
        }
        let adhered: Vec<usize> = self
            .platelets
            .iter()
            .filter(|p| p.state == PlateletState::Adhered)
            .map(|p| {
                let ix = clamp(p.pos[0], 0.0, (nx - 1) as f64) as usize;
                let iy = clamp(p.pos[1], 0.0, (self.ny - 1) as f64) as usize;
                ix + iy * nx
            })
            .collect();
        for &idx in &adhered {
            let act = 1.0;
            self.coag[idx].step(act, dt);
            let fibrin = self.coag[idx].fibrin;
            self.thrombus_fraction[idx] = clamp(
                self.thrombus_fraction[idx] + self.growth_rate * fibrin * dt,
                0.0,
                1.0,
            );
        }
    }
    /// Total thrombus volume (sum of fractions × cell volume).
    pub fn total_thrombus_volume(&self, cell_volume: f64) -> f64 {
        self.thrombus_fraction.iter().sum::<f64>() * cell_volume
    }
    /// Number of adhered platelets.
    pub fn adhered_count(&self) -> usize {
        self.platelets
            .iter()
            .filter(|p| p.state == PlateletState::Adhered)
            .count()
    }
}
/// Parameters governing red blood cell membrane mechanics.
#[derive(Debug, Clone)]
pub struct RbcMembraneParams {
    /// Shear modulus of the RBC membrane (N/m).
    pub shear_modulus: f64,
    /// Bending modulus (N·m).
    pub bending_modulus: f64,
    /// Area-incompressibility modulus (N/m).
    pub area_modulus: f64,
    /// Rest length of each membrane spring (m).
    pub rest_length: f64,
    /// Immersed boundary coupling coefficient.
    pub ib_stiffness: f64,
}
/// 2-D LBM microfluidics simulation supporting Hele-Shaw, Dean flow, and droplet dynamics.
pub struct MicrofluidicsLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Carrier fluid distribution functions \[nx * ny * NQ\].
    pub f: Vec<f64>,
    /// Droplet phase distribution functions \[nx * ny * NQ\].
    pub g: Vec<f64>,
    /// Density of carrier phase.
    pub rho_c: Vec<f64>,
    /// Density of droplet phase.
    pub rho_d: Vec<f64>,
    /// x-velocity.
    pub ux: Vec<f64>,
    /// y-velocity.
    pub uy: Vec<f64>,
    /// Relaxation time for carrier phase.
    pub tau_c: f64,
    /// Relaxation time for droplet phase.
    pub tau_d: f64,
    /// Surface tension parameter (Shan-Chen coupling).
    pub g_sc: f64,
    /// Collection of droplets (Lagrangian tracking).
    pub droplets: Vec<Droplet>,
    /// Hele-Shaw parameters.
    pub hele_shaw: HeleShawParams,
    /// Dean flow parameters.
    pub dean_flow: DeanFlowParams,
    /// Solid mask.
    pub is_wall: Vec<bool>,
}
impl MicrofluidicsLbm {
    /// Create a new microfluidics LBM simulation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        tau_c: f64,
        tau_d: f64,
        g_sc: f64,
        hele_shaw: HeleShawParams,
        dean_flow: DeanFlowParams,
    ) -> Self {
        let n = nx * ny;
        let is_wall: Vec<bool> = (0..n)
            .map(|idx| {
                let iy = idx / nx;
                iy == 0 || iy == ny - 1
            })
            .collect();
        let init_f: Vec<f64> = (0..n * NQ).map(|k| W[k % NQ]).collect();
        let init_g: Vec<f64> = (0..n * NQ).map(|k| W[k % NQ] * 0.1).collect();
        Self {
            nx,
            ny,
            f: init_f,
            g: init_g,
            rho_c: vec![1.0; n],
            rho_d: vec![0.1; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            tau_c,
            tau_d,
            g_sc,
            droplets: Vec::new(),
            hele_shaw,
            dean_flow,
            is_wall,
        }
    }
    /// Add a droplet to the simulation.
    pub fn add_droplet(&mut self, droplet: Droplet) {
        self.droplets.push(droplet);
    }
    /// Compute Hele-Shaw effective body force for current pressure gradient.
    pub fn hele_shaw_force(&self, dp_dx: f64) -> f64 {
        self.hele_shaw.darcy_velocity(dp_dx)
    }
    /// Perform one LBM step (BGK, Shan-Chen pseudo-potential for multiphase).
    pub fn step(&mut self, pressure_gradient: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        for idx in 0..n {
            let mut rc = 0.0f64;
            let mut rd = 0.0f64;
            let mut ux = 0.0f64;
            let mut uy = 0.0f64;
            for q in 0..NQ {
                rc += self.f[idx * NQ + q];
                rd += self.g[idx * NQ + q];
                ux += (self.f[idx * NQ + q] + self.g[idx * NQ + q]) * CX[q];
                uy += (self.f[idx * NQ + q] + self.g[idx * NQ + q]) * CY[q];
            }
            self.rho_c[idx] = rc;
            self.rho_d[idx] = rd;
            let rho_tot = rc + rd;
            self.ux[idx] = if rho_tot > 1e-12 { ux / rho_tot } else { 0.0 };
            self.uy[idx] = if rho_tot > 1e-12 { uy / rho_tot } else { 0.0 };
        }
        let omega_c = 1.0 / self.tau_c;
        let omega_d = 1.0 / self.tau_d;
        let g_sc = self.g_sc;
        let mut f_new = self.f.clone();
        let mut g_new = self.g.clone();
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = ix + iy * nx;
                if self.is_wall[idx] {
                    let opp = [0usize, 3, 4, 1, 2, 7, 8, 5, 6];
                    for q in 0..NQ {
                        f_new[idx * NQ + opp[q]] = self.f[idx * NQ + q];
                        g_new[idx * NQ + opp[q]] = self.g[idx * NQ + q];
                    }
                    continue;
                }
                let psi_c = self.rho_c[idx].clamp(1e-6, 1.0).ln();
                let psi_d = self.rho_d[idx].clamp(1e-6, 1.0).ln();
                let mut fx_sc_c = 0.0f64;
                let mut fy_sc_c = 0.0f64;
                let mut fx_sc_d = 0.0f64;
                let mut fy_sc_d = 0.0f64;
                for q in 1..NQ {
                    let nx2 = ((ix as isize + CX[q] as isize).rem_euclid(nx as isize)) as usize;
                    let ny2 = ((iy as isize + CY[q] as isize).rem_euclid(ny as isize)) as usize;
                    let idx2 = nx2 + ny2 * nx;
                    let psi_c2 = self.rho_c[idx2].clamp(1e-6, 1.0).ln();
                    let psi_d2 = self.rho_d[idx2].clamp(1e-6, 1.0).ln();
                    fx_sc_c += W[q] * CX[q] * psi_c2;
                    fy_sc_c += W[q] * CY[q] * psi_c2;
                    fx_sc_d += W[q] * CX[q] * psi_d2;
                    fy_sc_d += W[q] * CY[q] * psi_d2;
                    let _ = psi_c2;
                    let _ = psi_d2;
                }
                let fx_c = -g_sc * psi_c * fx_sc_d + pressure_gradient;
                let fy_c = -g_sc * psi_c * fy_sc_d;
                let fx_d = -g_sc * psi_d * fx_sc_c;
                let fy_d = -g_sc * psi_d * fy_sc_c;
                let ux_eq_c = self.ux[idx] + self.tau_c * fx_c / self.rho_c[idx].max(1e-12);
                let uy_eq_c = self.uy[idx] + self.tau_c * fy_c / self.rho_c[idx].max(1e-12);
                let ux_eq_d = self.ux[idx] + self.tau_d * fx_d / self.rho_d[idx].max(1e-12);
                let uy_eq_d = self.uy[idx] + self.tau_d * fy_d / self.rho_d[idx].max(1e-12);
                let feq_c = feq(self.rho_c[idx], ux_eq_c, uy_eq_c);
                let feq_d = feq(self.rho_d[idx], ux_eq_d, uy_eq_d);
                for q in 0..NQ {
                    f_new[idx * NQ + q] =
                        self.f[idx * NQ + q] * (1.0 - omega_c) + omega_c * feq_c[q];
                    g_new[idx * NQ + q] =
                        self.g[idx * NQ + q] * (1.0 - omega_d) + omega_d * feq_d[q];
                }
            }
        }
        self.f = f_new;
        self.g = g_new;
        for d in self.droplets.iter_mut() {
            let ix = clamp(d.pos[0], 0.0, (nx - 1) as f64) as usize;
            let iy = clamp(d.pos[1], 0.0, (ny - 1) as f64) as usize;
            let idx = ix + iy * nx;
            d.vel[0] = self.ux[idx];
            d.vel[1] = self.uy[idx];
            d.advance(1.0);
        }
    }
    /// Run `n_steps` steps.
    pub fn run(&mut self, n_steps: usize, pressure_gradient: f64) {
        for _ in 0..n_steps {
            self.step(pressure_gradient);
        }
    }
    /// Number of active droplets.
    pub fn droplet_count(&self) -> usize {
        self.droplets.len()
    }
}
/// Parameters for the Carreau-Yasuda non-Newtonian viscosity model.
///
/// The apparent viscosity is:
/// η = η_∞ + (η_0 − η_∞) \[1 + (λ γ̇)^a\]^((n-1)/a)
#[derive(Debug, Clone)]
pub struct CarreauYasudaParams {
    /// Zero-shear-rate viscosity (Pa·s).
    pub eta_0: f64,
    /// Infinite-shear-rate viscosity (Pa·s).
    pub eta_inf: f64,
    /// Relaxation time (s).
    pub lambda: f64,
    /// Power-law index (dimensionless).
    pub n: f64,
    /// Yasuda exponent (dimensionless); a = 2 → Carreau model.
    pub a: f64,
}
impl CarreauYasudaParams {
    /// Blood rheology parameters (Cho & Kensey 1991).
    pub fn blood() -> Self {
        Self {
            eta_0: 0.056,
            eta_inf: 0.00345,
            lambda: 3.313,
            n: 0.3568,
            a: 2.0,
        }
    }
    /// Compute apparent dynamic viscosity for a given shear rate γ̇ (s⁻¹).
    pub fn viscosity(&self, gamma_dot: f64) -> f64 {
        let gamma_dot = gamma_dot.abs();
        let lambda_gamma = self.lambda * gamma_dot;
        let base = 1.0 + lambda_gamma.powf(self.a);
        self.eta_inf + (self.eta_0 - self.eta_inf) * base.powf((self.n - 1.0) / self.a)
    }
}
/// A single platelet in the thrombosis model.
#[derive(Debug, Clone)]
pub struct Platelet {
    /// Position (x, y) in lattice units.
    pub pos: [f64; 2],
    /// Velocity (x, y) in lattice units.
    pub vel: [f64; 2],
    /// Activation level ∈ \[0, 1\].
    pub activation: f64,
    /// Current state.
    pub state: PlateletState,
    /// Accumulated shear exposure (Pa·s).
    pub shear_history: f64,
}
impl Platelet {
    /// Create a new resting platelet at position `pos`.
    pub fn new(pos: [f64; 2]) -> Self {
        Self {
            pos,
            vel: [0.0; 2],
            activation: 0.0,
            state: PlateletState::Resting,
            shear_history: 0.0,
        }
    }
    /// Update platelet state based on local shear stress and time step.
    ///
    /// Uses the Soares model: dA/dt = k_act · τ · (1 − A).
    pub fn update_activation(&mut self, shear_stress: f64, dt: f64, k_act: f64, tau_thresh: f64) {
        if self.state == PlateletState::Adhered {
            return;
        }
        self.shear_history += shear_stress * dt;
        let rate = k_act * (shear_stress - tau_thresh).max(0.0) * (1.0 - self.activation);
        self.activation = clamp(self.activation + rate * dt, 0.0, 1.0);
        self.state = if self.activation > 0.9 {
            PlateletState::Adhered
        } else if self.activation > 0.1 {
            PlateletState::Activated
        } else {
            PlateletState::Resting
        };
    }
}
/// Immersed-boundary RBC model embedded in a 2-D LBM fluid.
///
/// The RBC membrane is discretised as a ring of Lagrangian marker points
/// connected by springs (shear elasticity) plus bending potentials.
/// The immersed boundary method couples Lagrangian forces to the Eulerian grid.
pub struct RedBloodCellLbm {
    /// Number of membrane marker points.
    pub n_markers: usize,
    /// Marker positions (x, y) in lattice units.
    pub marker_pos: Vec<[f64; 2]>,
    /// Marker velocities (x, y) in lattice units.
    pub marker_vel: Vec<[f64; 2]>,
    /// Membrane parameters.
    pub params: RbcMembraneParams,
    /// Eulerian force density x-component \[nx * ny\].
    pub force_x: Vec<f64>,
    /// Eulerian force density y-component \[nx * ny\].
    pub force_y: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Resting area (lattice units²).
    pub rest_area: f64,
}
impl RedBloodCellLbm {
    /// Create a circular RBC membrane with `n_markers` points centred at (`cx`, `cy`),
    /// radius `r`, embedded in a grid of size `nx` × `ny`.
    pub fn new_circle(
        n_markers: usize,
        cx: f64,
        cy: f64,
        r: f64,
        nx: usize,
        ny: usize,
        params: RbcMembraneParams,
    ) -> Self {
        let marker_pos: Vec<[f64; 2]> = (0..n_markers)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n_markers as f64;
                [cx + r * theta.cos(), cy + r * theta.sin()]
            })
            .collect();
        let rest_area = PI * r * r;
        Self {
            n_markers,
            marker_pos,
            marker_vel: vec![[0.0; 2]; n_markers],
            params,
            force_x: vec![0.0; nx * ny],
            force_y: vec![0.0; nx * ny],
            nx,
            ny,
            rest_area,
        }
    }
    /// Compute shear spring force on marker `i` from its two neighbours.
    pub fn shear_force(&self, i: usize) -> [f64; 2] {
        let nm = self.n_markers;
        let prev = (i + nm - 1) % nm;
        let next = (i + 1) % nm;
        let pi = self.marker_pos[i];
        let pp = self.marker_pos[prev];
        let pn = self.marker_pos[next];
        let dp = sub2(pi, pp);
        let dn = sub2(pi, pn);
        let lp = len2(dp);
        let ln = len2(dn);
        let kappa = self.params.shear_modulus;
        let l0 = self.params.rest_length;
        let fp = if lp > 1e-12 {
            scale2(dp, -kappa * (lp - l0) / lp)
        } else {
            [0.0; 2]
        };
        let fn_ = if ln > 1e-12 {
            scale2(dn, -kappa * (ln - l0) / ln)
        } else {
            [0.0; 2]
        };
        add2(fp, fn_)
    }
    /// Compute bending force on marker `i` (discrete curvature penalty).
    pub fn bending_force(&self, i: usize) -> [f64; 2] {
        let nm = self.n_markers;
        let prev = (i + nm - 1) % nm;
        let next = (i + 1) % nm;
        let pi = self.marker_pos[i];
        let pp = self.marker_pos[prev];
        let pn = self.marker_pos[next];
        let laplacian = [pp[0] - 2.0 * pi[0] + pn[0], pp[1] - 2.0 * pi[1] + pn[1]];
        scale2(laplacian, self.params.bending_modulus)
    }
    /// Compute current enclosed area using the shoelace formula.
    pub fn current_area(&self) -> f64 {
        let nm = self.n_markers;
        let mut area = 0.0f64;
        for i in 0..nm {
            let j = (i + 1) % nm;
            area += self.marker_pos[i][0] * self.marker_pos[j][1];
            area -= self.marker_pos[j][0] * self.marker_pos[i][1];
        }
        area.abs() * 0.5
    }
    /// Compute area-conservation force on marker `i`.
    pub fn area_force(&self, i: usize) -> [f64; 2] {
        let area = self.current_area();
        let nm = self.n_markers;
        let prev = (i + nm - 1) % nm;
        let next = (i + 1) % nm;
        let chord = sub2(self.marker_pos[next], self.marker_pos[prev]);
        let normal = [chord[1], -chord[0]];
        let delta_a = area - self.rest_area;
        scale2(normal, -self.params.area_modulus * delta_a / (nm as f64))
    }
    /// Spread Lagrangian forces to Eulerian grid using delta function (IB method).
    pub fn spread_forces(&mut self, fluid_ux: &[f64], fluid_uy: &[f64]) {
        for v in self.force_x.iter_mut() {
            *v = 0.0;
        }
        for v in self.force_y.iter_mut() {
            *v = 0.0;
        }
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..self.n_markers {
            let f_shear = self.shear_force(i);
            let f_bend = self.bending_force(i);
            let f_area = self.area_force(i);
            let fx_lag = f_shear[0] + f_bend[0] + f_area[0];
            let fy_lag = f_shear[1] + f_bend[1] + f_area[1];
            let xm = self.marker_pos[i][0];
            let ym = self.marker_pos[i][1];
            let ix0 = xm.floor() as isize;
            let iy0 = ym.floor() as isize;
            let frac_x = xm - ix0 as f64;
            let frac_y = ym - iy0 as f64;
            let mut u_interp = [0.0f64; 2];
            for (djx, djy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                let ix = (ix0 + djx).rem_euclid(nx as isize) as usize;
                let iy = (iy0 + djy).rem_euclid(ny as isize) as usize;
                let w = (if djx == 0 { 1.0 - frac_x } else { frac_x })
                    * (if djy == 0 { 1.0 - frac_y } else { frac_y });
                let idx = ix + iy * nx;
                u_interp[0] += w * fluid_ux[idx];
                u_interp[1] += w * fluid_uy[idx];
            }
            let ib_fx = self.params.ib_stiffness * (u_interp[0] - self.marker_vel[i][0]) + fx_lag;
            let ib_fy = self.params.ib_stiffness * (u_interp[1] - self.marker_vel[i][1]) + fy_lag;
            for (djx, djy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                let ix = (ix0 + djx).rem_euclid(nx as isize) as usize;
                let iy = (iy0 + djy).rem_euclid(ny as isize) as usize;
                let w = (if djx == 0 { 1.0 - frac_x } else { frac_x })
                    * (if djy == 0 { 1.0 - frac_y } else { frac_y });
                let idx = ix + iy * nx;
                self.force_x[idx] += w * ib_fx;
                self.force_y[idx] += w * ib_fy;
            }
            self.marker_vel[i] = u_interp;
        }
    }
    /// Advance marker positions by `dt` using current marker velocities.
    pub fn advance_markers(&mut self, dt: f64) {
        for i in 0..self.n_markers {
            self.marker_pos[i][0] += self.marker_vel[i][0] * dt;
            self.marker_pos[i][1] += self.marker_vel[i][1] * dt;
        }
    }
    /// Compute membrane perimeter.
    pub fn perimeter(&self) -> f64 {
        let nm = self.n_markers;
        let mut perim = 0.0f64;
        for i in 0..nm {
            let j = (i + 1) % nm;
            perim += len2(sub2(self.marker_pos[j], self.marker_pos[i]));
        }
        perim
    }
    /// Compute circularity index: 4π A / P².
    pub fn circularity(&self) -> f64 {
        let a = self.current_area();
        let p = self.perimeter();
        if p < 1e-12 {
            return 0.0;
        }
        4.0 * PI * a / (p * p)
    }
}
/// Activation state of a platelet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlateletState {
    /// Resting (unactivated) platelet.
    Resting,
    /// Partially activated.
    Activated,
    /// Fully activated and adhered.
    Adhered,
}
/// Three-element Windkessel model for arterial pressure-flow coupling.
///
/// The 3-element Windkessel consists of:
/// - Characteristic impedance R₁ (proximal resistance, Pa·s/m³)
/// - Compliance C (arterial compliance, m³/Pa)
/// - Peripheral resistance R₂ (distal resistance, Pa·s/m³)
///
/// Governing ODE: C dP/dt = Q − (P − P_out) / R₂ − R₁ C dQ/dt
pub struct ArteryModel {
    /// Proximal resistance (Pa·s/m³).
    pub r1: f64,
    /// Arterial compliance (m³/Pa).
    pub c: f64,
    /// Peripheral resistance (Pa·s/m³).
    pub r2: f64,
    /// Venous (outflow) pressure (Pa).
    pub p_out: f64,
    /// Current arterial pressure (Pa).
    pub pressure: f64,
    /// Wall compliance (m²/Pa) – radial distension.
    pub wall_compliance: f64,
    /// Pulse wave velocity (m/s).
    pub pwv: f64,
    /// Previous flow rate (m³/s) for derivative approximation.
    pub q_prev: f64,
}
impl ArteryModel {
    /// Construct a new Windkessel artery model.
    ///
    /// `pwv` is the pulse wave velocity: c = √(Eh / (2ρR)) (Moens-Korteweg).
    pub fn new(r1: f64, c: f64, r2: f64, p_out: f64, pwv: f64, wall_compliance: f64) -> Self {
        Self {
            r1,
            c,
            r2,
            p_out,
            pressure: p_out + 10000.0,
            wall_compliance,
            pwv,
            q_prev: 0.0,
        }
    }
    /// Advance the Windkessel ODE by one time step `dt` given inflow `q` (m³/s).
    ///
    /// Uses implicit Euler for the RC term to ensure unconditional stability:
    /// P^(n+1) = (P^n + dt/C * (Q - P_out/R2 - R1 dQ/dt)) / (1 + dt/(C R2))
    pub fn step(&mut self, q: f64, dt: f64) {
        if self.c <= 0.0 {
            return;
        }
        let dq_dt = (q - self.q_prev) / dt.max(1e-30);
        let rhs = q / self.c
            - (self.pressure - self.p_out) / (self.c * self.r2.max(1e-30))
            - self.r1 / self.c * dq_dt;
        let denom = 1.0 + dt / (self.c * self.r2.max(1e-30));
        let dp_dt = rhs / denom;
        self.pressure += dp_dt * dt;
        self.q_prev = q;
    }
    /// Compute radial displacement of arterial wall due to pressure (m).
    ///
    /// δr = C_w · (P − P_out) / (2π r₀ L) where C_w = wall_compliance.
    pub fn wall_displacement(&self, r0: f64, length: f64) -> f64 {
        let delta_p = self.pressure - self.p_out;
        if r0 < 1e-12 || length < 1e-12 {
            return 0.0;
        }
        self.wall_compliance * delta_p / (2.0 * PI * r0 * length)
    }
    /// Moens-Korteweg pulse wave velocity: c = √(E h / (2 ρ R)).
    pub fn moens_korteweg(
        youngs_modulus: f64,
        wall_thickness: f64,
        density: f64,
        radius: f64,
    ) -> f64 {
        if density <= 0.0 || radius <= 0.0 {
            return 0.0;
        }
        (youngs_modulus * wall_thickness / (2.0 * density * radius)).sqrt()
    }
    /// Current pressure (Pa).
    pub fn current_pressure(&self) -> f64 {
        self.pressure
    }
}
/// Single cell in the biofiltration LBM simulation.
#[derive(Debug, Clone)]
pub struct BiofiltrationCell {
    /// Biofilm volume fraction ∈ \[0, 1\].
    pub biofilm_fraction: f64,
    /// Nutrient concentration (g/m³).
    pub nutrient: f64,
    /// Porosity ε = 1 − biofilm_fraction.
    pub porosity: f64,
    /// Local permeability (m²) via Kozeny-Carman.
    pub permeability: f64,
}
impl BiofiltrationCell {
    /// Create a new empty cell with full porosity.
    pub fn new(nutrient: f64) -> Self {
        Self {
            biofilm_fraction: 0.0,
            nutrient,
            porosity: 1.0,
            permeability: 1.0,
        }
    }
    /// Update porosity and permeability from current biofilm fraction.
    pub fn update_porous_properties(&mut self, grain_diameter: f64) {
        self.porosity = clamp(1.0 - self.biofilm_fraction, 0.05, 1.0);
        let eps = self.porosity;
        let one_m_eps = 1.0 - eps;
        self.permeability = if one_m_eps < 1e-12 {
            0.0
        } else {
            eps * eps * eps * grain_diameter * grain_diameter / (180.0 * one_m_eps * one_m_eps)
        };
    }
}
/// Coagulation cascade state (simplified 4-species model).
#[derive(Debug, Clone)]
pub struct CoagulationState {
    /// Prothrombin concentration (nM).
    pub prothrombin: f64,
    /// Thrombin concentration (nM).
    pub thrombin: f64,
    /// Fibrinogen concentration (nM).
    pub fibrinogen: f64,
    /// Fibrin concentration (nM).
    pub fibrin: f64,
}
impl CoagulationState {
    /// Physiological initial concentrations.
    pub fn physiological() -> Self {
        Self {
            prothrombin: 1400.0,
            thrombin: 0.0,
            fibrinogen: 8800.0,
            fibrin: 0.0,
        }
    }
    /// Advance coagulation cascade by `dt` seconds with platelet activation `act` ∈ \[0,1\].
    ///
    /// Simplified kinetics: thrombin generation from prothrombin driven by activation.
    pub fn step(&mut self, act: f64, dt: f64) {
        let k_pt = 0.5 * act;
        let k_fi = 0.1;
        let k_deg = 0.01;
        let d_thrombin = k_pt * self.prothrombin * dt - k_deg * self.thrombin * dt;
        let d_prothrombin = -k_pt * self.prothrombin * dt;
        let d_fibrin = k_fi * self.thrombin * self.fibrinogen * dt / 1000.0;
        let d_fibrinogen = -d_fibrin;
        self.thrombin = (self.thrombin + d_thrombin).max(0.0);
        self.prothrombin = (self.prothrombin + d_prothrombin).max(0.0);
        self.fibrin = (self.fibrin + d_fibrin).max(0.0);
        self.fibrinogen = (self.fibrinogen + d_fibrinogen).max(0.0);
    }
}
/// 2-D Lattice Boltzmann blood-flow simulator with Carreau-Yasuda viscosity and
/// pulsatile (Womersley) forcing.
///
/// # Overview
/// - D2Q9 lattice
/// - BGK collision with locally varying relaxation time τ(x,y) computed from
///   the Carreau-Yasuda model and the local shear rate.
/// - External body force via Guo's scheme.
/// - Pulsatile pressure gradient: G(t) = G_mean + G_amp · sin(2π f_heart t).
pub struct BloodFlowLbm {
    /// Grid width (lattice units).
    pub nx: usize,
    /// Grid height (lattice units).
    pub ny: usize,
    /// Distribution functions \[nx * ny * NQ\].
    pub f: Vec<f64>,
    /// Density field \[nx * ny\].
    pub rho: Vec<f64>,
    /// x-velocity field \[nx * ny\].
    pub ux: Vec<f64>,
    /// y-velocity field \[nx * ny\].
    pub uy: Vec<f64>,
    /// Carreau-Yasuda parameters.
    pub cy_params: CarreauYasudaParams,
    /// Mean pressure gradient (lattice units).
    pub g_mean: f64,
    /// Amplitude of pulsatile pressure gradient (lattice units).
    pub g_amp: f64,
    /// Heart-rate frequency (lattice units, cycles per timestep).
    pub f_heart: f64,
    /// Current time step.
    pub time_step: usize,
    /// Womersley number (dimensionless).
    pub womersley: f64,
    /// Solid mask: true = wall cell.
    pub is_wall: Vec<bool>,
}
impl BloodFlowLbm {
    /// Construct a new blood-flow LBM simulation.
    ///
    /// `womersley` = R √(ω ρ / μ_0) where R is vessel radius, ω = 2π f_heart,
    /// ρ is density and μ_0 is zero-shear viscosity.
    pub fn new(
        nx: usize,
        ny: usize,
        cy_params: CarreauYasudaParams,
        g_mean: f64,
        g_amp: f64,
        f_heart: f64,
        womersley: f64,
    ) -> Self {
        let n = nx * ny;
        let rho0 = 1.0f64;
        let f: Vec<f64> = (0..n * NQ)
            .map(|k| {
                let q = k % NQ;
                W[q] * rho0
            })
            .collect();
        let is_wall: Vec<bool> = (0..n)
            .map(|idx| {
                let iy = idx / nx;
                iy == 0 || iy == ny - 1
            })
            .collect();
        Self {
            nx,
            ny,
            f,
            rho: vec![rho0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            cy_params,
            g_mean,
            g_amp,
            f_heart,
            time_step: 0,
            womersley,
            is_wall,
        }
    }
    /// Compute instantaneous pressure gradient at current time.
    pub fn pressure_gradient(&self) -> f64 {
        let t = self.time_step as f64;
        self.g_mean + self.g_amp * (2.0 * PI * self.f_heart * t).sin()
    }
    /// Compute local shear rate magnitude at cell (ix, iy) from velocity field.
    pub fn shear_rate(&self, ix: usize, iy: usize) -> f64 {
        let nx = self.nx;
        let xp = ((ix + 1) % nx) + iy * nx;
        let xm = ix.wrapping_sub(1).min(nx - 1) + iy * nx;
        let yp = ix + ((iy + 1).min(self.ny - 1)) * nx;
        let ym = ix + iy.saturating_sub(1) * nx;
        let dux_dx = (self.ux[xp] - self.ux[xm]) * 0.5;
        let duy_dy = (self.uy[yp] - self.uy[ym]) * 0.5;
        let dux_dy = (self.ux[yp] - self.ux[ym]) * 0.5;
        let duy_dx = (self.uy[xp] - self.uy[xm]) * 0.5;
        let s_xx = dux_dx;
        let s_yy = duy_dy;
        let s_xy = 0.5 * (dux_dy + duy_dx);
        (2.0 * (s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy)).sqrt()
    }
    /// Compute local relaxation time τ from local shear rate.
    pub fn local_tau(&self, gamma_dot: f64) -> f64 {
        let nu = self.cy_params.viscosity(gamma_dot);
        3.0 * nu + 0.5
    }
    /// Perform one LBM time step (collision + streaming + macroscopic update).
    pub fn step(&mut self) {
        let g = self.pressure_gradient();
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut f_new = self.f.clone();
        for idx in 0..n {
            if self.is_wall[idx] {
                continue;
            }
            let ix = idx % nx;
            let iy = idx / nx;
            let gamma = self.shear_rate(ix, iy);
            let tau = self.local_tau(gamma);
            let omega = 1.0 / tau;
            let rho = self.rho[idx];
            let ux = self.ux[idx];
            let uy = self.uy[idx];
            let f_eq = feq(rho, ux, uy);
            for q in 0..NQ {
                let f_force = W[q]
                    * rho
                    * ((CX[q] - ux) / CS2 + CX[q] * (CX[q] * ux + CY[q] * uy) / (CS2 * CS2))
                    * g;
                f_new[idx * NQ + q] = self.f[idx * NQ + q] * (1.0 - omega)
                    + omega * f_eq[q]
                    + (1.0 - 0.5 * omega) * f_force;
            }
        }
        let mut f_stream = vec![0.0f64; n * NQ];
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = ix + iy * nx;
                if self.is_wall[idx] {
                    for q in 0..NQ {
                        let opp = [0, 3, 4, 1, 2, 7, 8, 5, 6][q];
                        f_stream[idx * NQ + opp] = f_new[idx * NQ + q];
                    }
                    continue;
                }
                for q in 0..NQ {
                    let ex = CX[q] as isize;
                    let ey = CY[q] as isize;
                    let nx_dest = (ix as isize + ex).rem_euclid(nx as isize) as usize;
                    let ny_dest = (iy as isize + ey).rem_euclid(ny as isize) as usize;
                    let dest = nx_dest + ny_dest * nx;
                    f_stream[dest * NQ + q] = f_new[idx * NQ + q];
                }
            }
        }
        self.f = f_stream;
        for idx in 0..n {
            if self.is_wall[idx] {
                continue;
            }
            let mut rho = 0.0f64;
            let mut ux = 0.0f64;
            let mut uy = 0.0f64;
            for q in 0..NQ {
                let fi = self.f[idx * NQ + q];
                rho += fi;
                ux += fi * CX[q];
                uy += fi * CY[q];
            }
            self.rho[idx] = rho;
            self.ux[idx] = ux / rho;
            self.uy[idx] = uy / rho;
        }
        self.time_step += 1;
    }
    /// Run `n_steps` time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }
    /// Compute mean x-velocity (flow rate proxy) over interior cells.
    pub fn mean_velocity(&self) -> f64 {
        let count: usize = self.is_wall.iter().filter(|&&w| !w).count();
        if count == 0 {
            return 0.0;
        }
        let sum: f64 = self
            .ux
            .iter()
            .zip(self.is_wall.iter())
            .filter(|&(_, w)| !w)
            .map(|(u, _)| u)
            .sum();
        sum / count as f64
    }
    /// Womersley number stored at construction time.
    pub fn womersley_number(&self) -> f64 {
        self.womersley
    }
}
/// Dean flow parameters for flow in a curved channel.
#[derive(Debug, Clone)]
pub struct DeanFlowParams {
    /// Channel hydraulic radius (m).
    pub radius: f64,
    /// Radius of curvature of the channel centreline (m).
    pub curvature_radius: f64,
    /// Mean flow velocity (m/s).
    pub mean_velocity: f64,
    /// Kinematic viscosity (m²/s).
    pub kinematic_viscosity: f64,
}
impl DeanFlowParams {
    /// Dean number: De = Re √(a / R_c).
    pub fn dean_number(&self) -> f64 {
        let re = self.mean_velocity * self.radius / self.kinematic_viscosity;
        re * (self.radius / self.curvature_radius).sqrt()
    }
    /// Secondary flow velocity scale u_sec ≈ 1.8 De^(1.5) ν / a (approximate).
    pub fn secondary_velocity(&self) -> f64 {
        let de = self.dean_number();
        if de < 1e-12 {
            return 0.0;
        }
        1.8 * de.powf(1.5) * self.kinematic_viscosity / self.radius
    }
    /// Pressure drop enhancement factor due to curvature (empirical).
    pub fn pressure_drop_factor(&self) -> f64 {
        let de = self.dean_number();
        if de > 100.0 {
            0.316 * de.powf(-0.25) * (1.0 + 0.033 * (de.log10()).powf(4.0))
        } else {
            1.0 + 0.033 * (de.log10().max(0.001)).powf(4.0)
        }
    }
}
/// Droplet in a microfluidic channel.
#[derive(Debug, Clone)]
pub struct Droplet {
    /// Position (x, y) (lattice units).
    pub pos: [f64; 2],
    /// Velocity (x, y) (lattice units/step).
    pub vel: [f64; 2],
    /// Droplet radius (lattice units).
    pub radius: f64,
    /// Internal viscosity ratio λ = μ_d / μ_c.
    pub viscosity_ratio: f64,
    /// Surface tension coefficient (lattice units).
    pub surface_tension: f64,
    /// Capillary number Ca = μ U / σ.
    pub capillary_number: f64,
}
impl Droplet {
    /// Create a new droplet.
    pub fn new(pos: [f64; 2], radius: f64, viscosity_ratio: f64, surface_tension: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 2],
            radius,
            viscosity_ratio,
            surface_tension,
            capillary_number: 0.0,
        }
    }
    /// Update capillary number given carrier velocity `u` and viscosity `mu`.
    pub fn update_capillary_number(&mut self, u: f64, mu: f64) {
        self.capillary_number = if self.surface_tension > 1e-30 {
            mu * u / self.surface_tension
        } else {
            0.0
        };
    }
    /// Critical capillary number above which droplet breaks (Taylor 1934, λ = 1).
    pub fn critical_capillary_number(&self) -> f64 {
        let lam = self.viscosity_ratio;
        0.5 * (1.0 + lam) / (1.0 + 1.5 * lam) * (19.0 * lam + 16.0) / (16.0 * lam + 16.0)
    }
    /// Check whether droplet will break under current Ca.
    pub fn will_break(&self) -> bool {
        self.capillary_number > self.critical_capillary_number()
    }
    /// Advance droplet position.
    pub fn advance(&mut self, dt: f64) {
        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
    }
}
/// 2-D Lattice Boltzmann biofiltration simulator.
///
/// Couples:
/// - LBM fluid flow through porous media (Brinkman equation via body force).
/// - Reaction-diffusion nutrient transport.
/// - Monod-kinetics biofilm growth.
pub struct BiofiltrationLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions for fluid \[nx * ny * NQ\].
    pub f: Vec<f64>,
    /// Nutrient distribution functions \[nx * ny * NQ\].
    pub g: Vec<f64>,
    /// Fluid density.
    pub rho: Vec<f64>,
    /// x-velocity.
    pub ux: Vec<f64>,
    /// y-velocity.
    pub uy: Vec<f64>,
    /// Nutrient concentration field.
    pub nutrient: Vec<f64>,
    /// Per-cell biofiltration state.
    pub cells: Vec<BiofiltrationCell>,
    /// Fluid relaxation time τ.
    pub tau_fluid: f64,
    /// Nutrient diffusion relaxation time.
    pub tau_nutrient: f64,
    /// Grain diameter for Kozeny-Carman (m).
    pub grain_diameter: f64,
    /// Monod half-saturation constant (g/m³).
    pub ks: f64,
    /// Maximum growth rate (1/s).
    pub mu_max: f64,
    /// Biofilm growth rate coefficient (per step).
    pub growth_coeff: f64,
    /// Inlet nutrient concentration (g/m³).
    pub c_inlet: f64,
}
impl BiofiltrationLbm {
    /// Construct a new biofiltration simulation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        tau_fluid: f64,
        tau_nutrient: f64,
        grain_diameter: f64,
        ks: f64,
        mu_max: f64,
        c_inlet: f64,
    ) -> Self {
        let n = nx * ny;
        let init_f: Vec<f64> = (0..n * NQ).map(|k| W[k % NQ]).collect();
        let init_g: Vec<f64> = (0..n * NQ).map(|k| W[k % NQ] * c_inlet).collect();
        let cells: Vec<BiofiltrationCell> = (0..n)
            .map(|idx| {
                let ix = idx % nx;
                let nutrient_init = if ix == 0 { c_inlet } else { 0.0 };
                BiofiltrationCell::new(nutrient_init)
            })
            .collect();
        Self {
            nx,
            ny,
            f: init_f,
            g: init_g,
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            nutrient: vec![0.0; n],
            cells,
            tau_fluid,
            tau_nutrient,
            grain_diameter,
            ks,
            mu_max,
            growth_coeff: 1e-5,
            c_inlet,
        }
    }
    /// Monod growth rate at concentration `c`.
    pub fn monod_rate(&self, c: f64) -> f64 {
        if c <= 0.0 || self.ks <= 0.0 {
            return 0.0;
        }
        self.mu_max * c / (self.ks + c)
    }
    /// Perform one simulation step.
    pub fn step(&mut self, pressure_gradient: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let omega_f = 1.0 / self.tau_fluid;
        let omega_n = 1.0 / self.tau_nutrient;
        for idx in 0..n {
            let mut rho = 0.0f64;
            let mut ux = 0.0f64;
            let mut uy = 0.0f64;
            let mut nutrient = 0.0f64;
            for q in 0..NQ {
                let fi = self.f[idx * NQ + q];
                let gi = self.g[idx * NQ + q];
                rho += fi;
                ux += fi * CX[q];
                uy += fi * CY[q];
                nutrient += gi;
            }
            self.rho[idx] = rho;
            self.ux[idx] = if rho > 1e-12 { ux / rho } else { 0.0 };
            self.uy[idx] = if rho > 1e-12 { uy / rho } else { 0.0 };
            self.nutrient[idx] = nutrient;
        }
        let mut f_new = vec![0.0f64; n * NQ];
        let mut g_new = vec![0.0f64; n * NQ];
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = ix + iy * nx;
                let eps = self.cells[idx].porosity;
                let perm = self.cells[idx].permeability;
                let nu = CS2 * (self.tau_fluid - 0.5);
                let brinkman_x = if perm > 1e-30 {
                    -nu / perm * self.ux[idx] * self.rho[idx] * eps
                } else {
                    0.0
                };
                let f_total_x = pressure_gradient + brinkman_x;
                let feq_f = feq(self.rho[idx], self.ux[idx], self.uy[idx]);
                let nutrient_c = self.nutrient[idx];
                let feq_n = feq(1.0, self.ux[idx], self.uy[idx]);
                let c = nutrient_c.max(0.0);
                let biofilm = self.cells[idx].biofilm_fraction;
                let consumption = self.monod_rate(c) * biofilm;
                for q in 0..NQ {
                    let f_force = W[q]
                        * self.rho[idx]
                        * ((CX[q] - self.ux[idx]) / CS2
                            + CX[q] * (CX[q] * self.ux[idx] + CY[q] * self.uy[idx]) / (CS2 * CS2))
                        * f_total_x;
                    let fi_coll = self.f[idx * NQ + q] * (1.0 - omega_f)
                        + omega_f * feq_f[q]
                        + (1.0 - 0.5 * omega_f) * f_force;
                    let gi_eq = feq_n[q] * c;
                    let gi_coll = self.g[idx * NQ + q] * (1.0 - omega_n) + omega_n * gi_eq
                        - W[q] * consumption;
                    let ex = CX[q] as isize;
                    let ey = CY[q] as isize;
                    let iy_wall_top = iy == ny - 1 && ey > 0;
                    let iy_wall_bot = iy == 0 && ey < 0;
                    if iy_wall_top || iy_wall_bot {
                        let opp = [0, 3, 4, 1, 2, 7, 8, 5, 6][q];
                        f_new[idx * NQ + opp] += fi_coll;
                        g_new[idx * NQ + opp] += gi_coll;
                    } else {
                        let nx2 = ((ix as isize + ex).rem_euclid(nx as isize)) as usize;
                        let ny2 = ((iy as isize + ey).rem_euclid(ny as isize)) as usize;
                        let dest = nx2 + ny2 * nx;
                        f_new[dest * NQ + q] += fi_coll;
                        g_new[dest * NQ + q] += gi_coll;
                    }
                }
            }
        }
        for iy in 0..ny {
            let idx = iy * nx;
            for q in 0..NQ {
                g_new[idx * NQ + q] = W[q] * self.c_inlet;
            }
        }
        self.f = f_new;
        self.g = g_new;
        for idx in 0..n {
            let c = self.nutrient[idx].max(0.0);
            let growth = self.growth_coeff * self.monod_rate(c);
            self.cells[idx].biofilm_fraction =
                clamp(self.cells[idx].biofilm_fraction + growth, 0.0, 0.95);
            self.cells[idx].update_porous_properties(self.grain_diameter);
        }
    }
    /// Run `n_steps` time steps.
    pub fn run(&mut self, n_steps: usize, pressure_gradient: f64) {
        for _ in 0..n_steps {
            self.step(pressure_gradient);
        }
    }
    /// Compute mean nutrient concentration across the grid.
    pub fn mean_nutrient(&self) -> f64 {
        if self.nutrient.is_empty() {
            return 0.0;
        }
        self.nutrient.iter().sum::<f64>() / self.nutrient.len() as f64
    }
    /// Compute total biofilm volume fraction (sum over all cells).
    pub fn total_biofilm(&self) -> f64 {
        self.cells.iter().map(|c| c.biofilm_fraction).sum()
    }
    /// Compute pressure drop across the filter (Darcy law estimate).
    pub fn pressure_drop(&self, length: f64, flow_velocity: f64) -> f64 {
        let mean_perm: f64 =
            self.cells.iter().map(|c| c.permeability).sum::<f64>() / self.cells.len() as f64;
        if mean_perm < 1e-30 {
            return f64::INFINITY;
        }
        let nu = CS2 * (self.tau_fluid - 0.5);
        nu * flow_velocity * length / mean_perm
    }
}
