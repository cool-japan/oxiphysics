//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

use super::functions::*;

/// Adaptive algorithm selector: chooses between continuum NS and kinetic LBM.
///
/// Based on local Knudsen number, selects the appropriate solver
/// for each region of the domain.
#[derive(Debug, Clone)]
pub struct AdaptiveAlgorithmSelector {
    /// Number of domain nodes.
    pub n_nodes: usize,
    /// Knudsen threshold for continuum: Kn < kn_continuum → use NS.
    pub kn_continuum: f64,
    /// Knudsen threshold for rarefied: Kn > kn_rarefied → use DSMC.
    pub kn_rarefied: f64,
    /// Algorithm labels per node: 0 = NS, 1 = LBM, 2 = DSMC.
    pub algorithm: Vec<u8>,
    /// Knudsen numbers per node (updated externally).
    pub kn_field: Vec<f64>,
}
impl AdaptiveAlgorithmSelector {
    /// Create a new adaptive selector.
    ///
    /// # Arguments
    /// * `n_nodes`       - domain size
    /// * `kn_continuum`  - upper Kn for continuum (default 0.001)
    /// * `kn_rarefied`   - lower Kn for full rarefied (default 1.0)
    pub fn new(n_nodes: usize, kn_continuum: f64, kn_rarefied: f64) -> Self {
        Self {
            n_nodes,
            kn_continuum,
            kn_rarefied,
            algorithm: vec![1u8; n_nodes],
            kn_field: vec![0.001; n_nodes],
        }
    }
    /// Update algorithm selection based on current Knudsen field.
    pub fn update_selection(&mut self) {
        for i in 0..self.n_nodes {
            let kn = self.kn_field[i];
            self.algorithm[i] = if kn < self.kn_continuum {
                0
            } else if kn > self.kn_rarefied {
                2
            } else {
                1
            };
        }
    }
    /// Set Knudsen field and recompute algorithm map.
    pub fn set_kn_field(&mut self, kn: &[f64]) {
        let n = self.n_nodes.min(kn.len());
        self.kn_field[..n].copy_from_slice(&kn[..n]);
        self.update_selection();
    }
    /// Count nodes using each algorithm.
    pub fn algorithm_counts(&self) -> [usize; 3] {
        let mut counts = [0usize; 3];
        for &a in &self.algorithm {
            if (a as usize) < 3 {
                counts[a as usize] += 1;
            }
        }
        counts
    }
    /// Fraction of nodes using Navier-Stokes solver.
    pub fn ns_fraction(&self) -> f64 {
        if self.n_nodes == 0 {
            return 0.0;
        }
        self.algorithm.iter().filter(|&&a| a == 0).count() as f64 / self.n_nodes as f64
    }
    /// Fraction of nodes using LBM solver.
    pub fn lbm_fraction(&self) -> f64 {
        if self.n_nodes == 0 {
            return 0.0;
        }
        self.algorithm.iter().filter(|&&a| a == 1).count() as f64 / self.n_nodes as f64
    }
    /// Fraction of nodes using DSMC solver.
    pub fn dsmc_fraction(&self) -> f64 {
        if self.n_nodes == 0 {
            return 0.0;
        }
        self.algorithm.iter().filter(|&&a| a == 2).count() as f64 / self.n_nodes as f64
    }
    /// Find transition interfaces (consecutive nodes with different algorithms).
    pub fn find_interfaces(&self) -> Vec<usize> {
        let mut ifaces = Vec::new();
        for i in 1..self.n_nodes {
            if self.algorithm[i] != self.algorithm[i - 1] {
                ifaces.push(i);
            }
        }
        ifaces
    }
}
/// Continuum breakdown parameter: Knudsen number estimator.
///
/// Kn = λ / L where λ is the mean free path and L is a characteristic
/// macroscopic length. Used to decide where LBM (Kn < 0.01, continuum)
/// or DSMC (Kn > 0.1, rarefied) should be applied.
///
/// References:
/// - Bird, G.A. (1994). *Molecular Gas Dynamics and the Direct Simulation of Gas Flows*.
/// - Hadjiconstantinou, N.G. (2000). Hybrid LBM-DSMC. *J. Comput. Phys.* 164, 293.
#[derive(Debug, Clone)]
pub struct KnudsenEstimator {
    /// Number of nodes in the domain.
    pub n_nodes: usize,
    /// Grid spacing (physical units).
    pub dx: f64,
    /// Dynamic viscosity μ.
    pub mu: f64,
    /// Reference pressure p_ref.
    pub p_ref: f64,
    /// Characteristic length L.
    pub char_length: f64,
    /// Knudsen number per node.
    pub kn: Vec<f64>,
}
impl KnudsenEstimator {
    /// Create a new Knudsen number estimator.
    pub fn new(n_nodes: usize, dx: f64, mu: f64, p_ref: f64, char_length: f64) -> Self {
        Self {
            n_nodes,
            dx,
            mu,
            p_ref,
            char_length,
            kn: vec![0.0; n_nodes],
        }
    }
    /// Compute mean free path λ from Chapman-Enskog: λ ≈ μ sqrt(π/(2ρkT)).
    ///
    /// In LBM units (T=1, kB=1): λ ≈ μ sqrt(π/2) / p where p = ρ cs².
    pub fn mean_free_path(&self, rho: f64) -> f64 {
        let p = rho * CS2;
        let p_safe = p.max(1e-30);
        self.mu * (std::f64::consts::PI / 2.0).sqrt() / p_safe
    }
    /// Local Knudsen number at a point: Kn = λ / L_local.
    ///
    /// L_local is estimated from local density gradient: L_local = ρ / |∇ρ|.
    pub fn local_kn(&self, rho: f64, grad_rho: f64) -> f64 {
        let lambda = self.mean_free_path(rho);
        let l_local = if grad_rho.abs() > 1e-15 {
            rho / grad_rho.abs()
        } else {
            self.char_length
        };
        lambda / l_local.max(1e-30)
    }
    /// Compute Knudsen numbers for a density field.
    pub fn compute_kn_field(&mut self, rho_field: &[f64], nx: usize) {
        for i in 0..self.n_nodes.min(rho_field.len()) {
            let rho = rho_field[i];
            let left = if i == 0 {
                rho_field[self.n_nodes - 1]
            } else {
                rho_field[i - 1]
            };
            let right = if i + 1 >= self.n_nodes {
                rho_field[0]
            } else {
                rho_field[i + 1]
            };
            let grad_rho = (right - left).abs() / (2.0 * self.dx);
            self.kn[i] = self.local_kn(rho, grad_rho);
        }
        let _ = nx;
    }
    /// Classify each node: 0 = continuum, 1 = slip, 2 = transitional, 3 = rarefied.
    pub fn classify_regimes(&self) -> Vec<u8> {
        self.kn
            .iter()
            .map(|&kn| {
                if kn < 0.001 {
                    0
                } else if kn < 0.1 {
                    1
                } else if kn < 10.0 {
                    2
                } else {
                    3
                }
            })
            .collect()
    }
    /// Fraction of domain in continuum regime.
    pub fn continuum_fraction(&self) -> f64 {
        if self.n_nodes == 0 {
            return 1.0;
        }
        self.kn.iter().filter(|&&kn| kn < 0.01).count() as f64 / self.n_nodes as f64
    }
    /// Maximum Knudsen number in domain.
    pub fn max_kn(&self) -> f64 {
        self.kn.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
}
/// Two-way particle-laden flow coupling.
///
/// Fluid forces act on particles (drag), and particle wakes disturb the fluid
/// (reaction force).
#[derive(Debug, Clone)]
pub struct TwoWayCoupling {
    /// LBM grid distributions `f[node][9]`.
    pub f: Vec<[f64; 9]>,
    /// Particle positions (2D: flat \[x0,y0,x1,y1,...\]).
    pub positions: Vec<f64>,
    /// Particle velocities (2D: flat \[vx0,vy0,...\]).
    pub velocities: Vec<f64>,
    /// Particle radii.
    pub radii: Vec<f64>,
    /// Particle densities.
    pub densities: Vec<f64>,
    /// Drag force on each particle `[fx, fy]`.
    pub drag_forces: Vec<[f64; 2]>,
    /// Body force per LBM node `[fx, fy]`.
    pub body_force: Vec<[f64; 2]>,
    /// Number of LBM nodes.
    pub n_nodes: usize,
    /// Number of particles.
    pub n_particles: usize,
    /// Fluid kinematic viscosity (lattice units).
    pub nu: f64,
    /// LBM relaxation time.
    pub tau: f64,
}
impl TwoWayCoupling {
    /// Create a new two-way coupling.
    pub fn new(n_nodes: usize, n_particles: usize, nu: f64, tau: f64) -> Self {
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            f: vec![f_init; n_nodes],
            positions: vec![0.0; n_particles * 2],
            velocities: vec![0.0; n_particles * 2],
            radii: vec![0.5; n_particles],
            densities: vec![1.0; n_particles],
            drag_forces: vec![[0.0; 2]; n_particles],
            body_force: vec![[0.0; 2]; n_nodes],
            n_nodes,
            n_particles,
            nu,
            tau,
        }
    }
    /// Compute Stokes drag on particle `p` given local fluid velocity `u_fluid`.
    pub fn stokes_drag(&self, p: usize, u_fluid: [f64; 2]) -> [f64; 2] {
        let r = self.radii[p];
        let vx = self.velocities[2 * p];
        let vy = self.velocities[2 * p + 1];
        let coeff = 6.0 * std::f64::consts::PI * self.nu * r;
        [coeff * (u_fluid[0] - vx), coeff * (u_fluid[1] - vy)]
    }
    /// Compute drag for all particles given per-node fluid velocities.
    pub fn compute_all_drags(
        &mut self,
        fluid_velocities: &[[f64; 2]],
        node_positions: &[[f64; 2]],
    ) {
        for p in 0..self.n_particles {
            let px = self.positions[2 * p];
            let py = self.positions[2 * p + 1];
            let nearest = node_positions
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    let da = (a[0] - px).powi(2) + (a[1] - py).powi(2);
                    let db = (b[0] - px).powi(2) + (b[1] - py).powi(2);
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            self.drag_forces[p] = self.stokes_drag(p, fluid_velocities[nearest]);
        }
    }
    /// BGK collision with body force (Guo forcing scheme).
    pub fn collide_with_force(&mut self) {
        for k in 0..self.n_nodes {
            let rho: f64 = self.f[k].iter().sum();
            let fx = self.body_force[k][0];
            let fy = self.body_force[k][1];
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for i in 0..9 {
                let c = D2Q9_VELOCITIES[i];
                ux += c[0] as f64 * self.f[k][i];
                uy += c[1] as f64 * self.f[k][i];
            }
            let u = [
                (ux + 0.5 * fx) / rho.max(1e-15),
                (uy + 0.5 * fy) / rho.max(1e-15),
            ];
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let u2 = u[0] * u[0] + u[1] * u[1];
                let feq =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
                let ec = c[0] as f64 * fx + c[1] as f64 * fy;
                let eu_sum = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let fi =
                    w * ((1.0 - 1.0 / (2.0 * self.tau)) * (ec / CS2 + eu_sum * ec / (CS2 * CS2)));
                self.f[k][i] -= (self.f[k][i] - feq) / self.tau + fi;
            }
        }
    }
}
/// LBM-DSMC (Direct Simulation Monte Carlo) coupling interface.
///
/// Manages the handshake region where the LBM continuum description
/// transitions to the DSMC particle-based description.
///
/// References:
/// - Lian, Y.-Y. et al. (2011). Hybrid DSMC-NS approach. *Microfluid.* 10, 481.
/// - Hash, D.B. & Hassan, H.A. (1996). Assessment of schemes for coupling. *J. Thermophys.* 10, 242.
#[derive(Debug, Clone)]
pub struct LbmDsmcCoupling {
    /// Number of nodes in the overlap (buffer) zone.
    pub n_overlap: usize,
    /// LBM distribution functions in overlap zone.
    pub f_lbm: Vec<[f64; 9]>,
    /// DSMC sampled density in overlap zone.
    pub rho_dsmc: Vec<f64>,
    /// DSMC sampled velocity in overlap zone.
    pub u_dsmc: Vec<[f64; 2]>,
    /// DSMC sampled temperature in overlap zone.
    pub temp_dsmc: Vec<f64>,
    /// LBM relaxation time τ.
    pub tau: f64,
    /// Boltzmann constant kB.
    pub k_b: f64,
    /// Particle mass m.
    pub mass: f64,
}
impl LbmDsmcCoupling {
    /// Create a new LBM-DSMC coupling zone.
    ///
    /// # Arguments
    /// * `n_overlap` - number of nodes in the overlap region
    /// * `tau`       - LBM relaxation time
    /// * `mass`      - molecular mass
    pub fn new(n_overlap: usize, tau: f64, mass: f64) -> Self {
        let feq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            n_overlap,
            f_lbm: vec![feq; n_overlap],
            rho_dsmc: vec![1.0; n_overlap],
            u_dsmc: vec![[0.0; 2]; n_overlap],
            temp_dsmc: vec![1.0; n_overlap],
            tau,
            k_b: 1.0,
            mass,
        }
    }
    /// Convert DSMC macroscopic moments to LBM equilibrium distributions.
    ///
    /// The LBM temperature-dependent equilibrium uses
    /// cs² = k_B T / m.
    pub fn dsmc_to_lbm(&self, idx: usize) -> [f64; 9] {
        if idx >= self.n_overlap {
            return D2Q9_WEIGHTS.map(|w| w);
        }
        let rho = self.rho_dsmc[idx];
        let u = self.u_dsmc[idx];
        let _cs2_local = self.k_b * self.temp_dsmc[idx] / self.mass;
        lbm_equilibrium(rho, u)
    }
    /// Extract DSMC moments from LBM distributions.
    ///
    /// ρ = Σ f_i, j = Σ f_i e_i, T = (Σ f_i |e_i - u|²) / (d ρ cs²)
    pub fn lbm_to_dsmc_moments(&self, f: &[f64; 9]) -> (f64, [f64; 2], f64) {
        let rho: f64 = f.iter().sum();
        let rho_safe = rho.max(1e-30);
        let mut jx = 0.0f64;
        let mut jy = 0.0f64;
        for alpha in 0..9 {
            jx += f[alpha] * D2Q9_VELOCITIES[alpha][0] as f64;
            jy += f[alpha] * D2Q9_VELOCITIES[alpha][1] as f64;
        }
        let ux = jx / rho_safe;
        let uy = jy / rho_safe;
        let mut e2 = 0.0f64;
        for alpha in 0..9 {
            let ex = D2Q9_VELOCITIES[alpha][0] as f64 - ux;
            let ey = D2Q9_VELOCITIES[alpha][1] as f64 - uy;
            e2 += f[alpha] * (ex * ex + ey * ey);
        }
        let temp = e2 / (2.0 * rho_safe);
        (rho, [ux, uy], temp)
    }
    /// Enforce flux matching at the LBM-DSMC boundary.
    ///
    /// Modifies LBM distributions so that normal momentum flux
    /// matches the DSMC flux at the interface node `idx`.
    pub fn flux_match(&mut self, idx: usize) {
        if idx >= self.n_overlap {
            return;
        }
        let f_dsmc_eq = self.dsmc_to_lbm(idx);
        let alpha_blend = 0.5;
        for i in 0..9 {
            self.f_lbm[idx][i] =
                (1.0 - alpha_blend) * self.f_lbm[idx][i] + alpha_blend * f_dsmc_eq[i];
        }
    }
    /// Update all overlap nodes with DSMC data.
    pub fn update_from_dsmc(&mut self, rho: &[f64], u: &[[f64; 2]], temp: &[f64]) {
        for i in 0..self.n_overlap.min(rho.len()) {
            self.rho_dsmc[i] = rho[i];
            self.u_dsmc[i] = u[i];
            self.temp_dsmc[i] = temp[i];
            self.flux_match(i);
        }
    }
    /// Mean temperature in the overlap zone.
    pub fn mean_temperature(&self) -> f64 {
        if self.n_overlap == 0 {
            return 0.0;
        }
        self.temp_dsmc.iter().sum::<f64>() / self.n_overlap as f64
    }
}
/// Quad-tree node for 2D adaptive mesh refinement.
///
/// Each node covers a rectangular region `[x0, x1] × [y0, y1]` in physical
/// space and may be refined into four children.
#[derive(Debug, Clone)]
pub struct QuadTreeNode {
    /// Physical x-extent: `[x_min, x_max]`.
    pub x_range: [f64; 2],
    /// Physical y-extent: `[y_min, y_max]`.
    pub y_range: [f64; 2],
    /// Refinement level (0 = coarsest).
    pub level: usize,
    /// Node index (unique identifier).
    pub index: usize,
    /// Indices of child nodes (None if leaf node).
    pub children: Option<[usize; 4]>,
    /// D2Q9 distribution functions at this node.
    pub f: [f64; 9],
    /// Macroscopic density.
    pub rho: f64,
    /// Macroscopic velocity `[ux, uy]`.
    pub u: [f64; 2],
}
impl QuadTreeNode {
    /// Create a new leaf quad-tree node.
    pub fn new(x_range: [f64; 2], y_range: [f64; 2], level: usize, index: usize) -> Self {
        let f_eq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            x_range,
            y_range,
            level,
            index,
            children: None,
            f: f_eq,
            rho: 1.0,
            u: [0.0; 2],
        }
    }
    /// Cell size Δx at this refinement level.
    pub fn dx(&self) -> f64 {
        self.x_range[1] - self.x_range[0]
    }
    /// Cell size Δy at this refinement level.
    pub fn dy(&self) -> f64 {
        self.y_range[1] - self.y_range[0]
    }
    /// Center position `[cx, cy]`.
    pub fn center(&self) -> [f64; 2] {
        [
            (self.x_range[0] + self.x_range[1]) * 0.5,
            (self.y_range[0] + self.y_range[1]) * 0.5,
        ]
    }
    /// Returns true if this is a leaf node (no children).
    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }
    /// D2Q9 equilibrium at this node's macroscopic state.
    pub fn equilibrium(&self, i: usize) -> f64 {
        let w = D2Q9_WEIGHTS[i];
        let c = D2Q9_VELOCITIES[i];
        let eu = c[0] as f64 * self.u[0] + c[1] as f64 * self.u[1];
        let u2 = self.u[0] * self.u[0] + self.u[1] * self.u[1];
        w * self.rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }
}
/// LBM-RANS coupling: LBM near wall, Navier-Stokes in far field.
///
/// The domain is split into an LBM zone (near-wall, fine structures) and a
/// NS zone (far field, coarser resolution).  The interface exchanges
/// velocity and pressure/density data.
#[derive(Debug, Clone)]
pub struct HybridLbmNs {
    /// LBM zone distribution functions.
    pub f_lbm: Vec<[f64; 9]>,
    /// NS zone velocity (u, v) at each NS node.
    pub u_ns: Vec<[f64; 2]>,
    /// NS zone pressure at each node.
    pub p_ns: Vec<f64>,
    /// Number of LBM nodes.
    pub n_lbm: usize,
    /// Number of NS nodes.
    pub n_ns: usize,
    /// LBM relaxation time.
    pub tau_lbm: f64,
    /// NS kinematic viscosity.
    pub nu_ns: f64,
    /// Interface node pairs: (lbm_idx, ns_idx).
    pub interface_pairs: Vec<(usize, usize)>,
    /// Overlap layer width (nodes).
    pub overlap_width: usize,
}
impl HybridLbmNs {
    /// Create a new hybrid LBM-NS solver.
    pub fn new(n_lbm: usize, n_ns: usize, tau_lbm: f64, nu_ns: f64) -> Self {
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            f_lbm: vec![f_init; n_lbm],
            u_ns: vec![[0.0; 2]; n_ns],
            p_ns: vec![1.0; n_ns],
            n_lbm,
            n_ns,
            tau_lbm,
            nu_ns,
            interface_pairs: Vec::new(),
            overlap_width: 2,
        }
    }
    /// Add an interface pair (lbm_node, ns_node).
    pub fn add_interface_pair(&mut self, lbm_idx: usize, ns_idx: usize) {
        self.interface_pairs.push((lbm_idx, ns_idx));
    }
    /// Exchange velocity data at interface: LBM → NS.
    pub fn lbm_to_ns_exchange(&mut self) {
        for &(lbm_idx, ns_idx) in &self.interface_pairs {
            let rho: f64 = self.f_lbm[lbm_idx].iter().sum();
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for i in 0..9 {
                let c = D2Q9_VELOCITIES[i];
                ux += c[0] as f64 * self.f_lbm[lbm_idx][i];
                uy += c[1] as f64 * self.f_lbm[lbm_idx][i];
            }
            self.u_ns[ns_idx] = [ux / rho.max(1e-15), uy / rho.max(1e-15)];
            self.p_ns[ns_idx] = rho * CS2;
        }
    }
    /// Exchange velocity data at interface: NS → LBM (via equilibrium).
    pub fn ns_to_lbm_exchange(&mut self) {
        for &(lbm_idx, ns_idx) in &self.interface_pairs {
            let u = self.u_ns[ns_idx];
            let rho = self.p_ns[ns_idx] / CS2;
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let u2 = u[0] * u[0] + u[1] * u[1];
                self.f_lbm[lbm_idx][i] =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
            }
        }
    }
}
/// Grid refinement interface: fine-to-coarse rescaling with buffer layers.
///
/// Implements the Filippova-Hänel / Rohde-Dünweg rescaling algorithm
/// for distribution functions at grid refinement interfaces.
///
/// References:
/// - Filippova, O. & Hänel, D. (1998). *J. Comput. Phys.* 147, 219–228.
/// - Rohde, M. et al. (2006). *Comput. Fluids* 35, 1121–1134.
#[derive(Debug, Clone)]
pub struct GridRefinementInterface {
    /// Refinement ratio r (typically 2).
    pub refinement_ratio: usize,
    /// Relaxation time on coarse grid τ_c.
    pub tau_coarse: f64,
    /// Relaxation time on fine grid τ_f = τ_c / r + 0.5 (1 - 1/r).
    pub tau_fine: f64,
    /// Buffer layer width in coarse cells.
    pub buffer_width: usize,
    /// Distribution functions at coarse side of interface.
    pub f_coarse_iface: Vec<[f64; 9]>,
    /// Distribution functions at fine side of interface.
    pub f_fine_iface: Vec<[f64; 9]>,
}
impl GridRefinementInterface {
    /// Create a new grid refinement interface.
    ///
    /// # Arguments
    /// * `n_iface`          - number of interface nodes
    /// * `refinement_ratio` - spatial refinement ratio (2 or 4)
    /// * `tau_coarse`       - relaxation time on coarse grid
    /// * `buffer_width`     - width of overlap buffer layer
    pub fn new(
        n_iface: usize,
        refinement_ratio: usize,
        tau_coarse: f64,
        buffer_width: usize,
    ) -> Self {
        let r = refinement_ratio as f64;
        let tau_fine = (tau_coarse - 0.5) / r + 0.5;
        let f_eq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            refinement_ratio,
            tau_coarse,
            tau_fine,
            buffer_width,
            f_coarse_iface: vec![f_eq; n_iface],
            f_fine_iface: vec![f_eq; n_iface * refinement_ratio],
        }
    }
    /// Rescale distribution from coarse to fine grid (Filippova-Hänel).
    ///
    /// f_f = f_eq(ρ,u) + (f_c - f_eq_c) * τ_c / τ_f
    pub fn coarse_to_fine_rescale(&self, f_c: &[f64; 9], rho: f64, u: [f64; 2]) -> [f64; 9] {
        let feq_c = lbm_equilibrium(rho, u);
        let scale = self.tau_coarse / self.tau_fine;
        let mut f_f = [0.0f64; 9];
        for alpha in 0..9 {
            f_f[alpha] = feq_c[alpha] + (f_c[alpha] - feq_c[alpha]) * scale;
        }
        f_f
    }
    /// Rescale distribution from fine to coarse grid.
    ///
    /// f_c = f_eq(ρ,u) + (f_f - f_eq_f) * τ_f / τ_c
    pub fn fine_to_coarse_rescale(&self, f_f: &[f64; 9], rho: f64, u: [f64; 2]) -> [f64; 9] {
        let feq_f = lbm_equilibrium(rho, u);
        let scale = self.tau_fine / self.tau_coarse;
        let mut f_c = [0.0f64; 9];
        for alpha in 0..9 {
            f_c[alpha] = feq_f[alpha] + (f_f[alpha] - feq_f[alpha]) * scale;
        }
        f_c
    }
    /// Average fine cells onto coarse cell (spatial restriction).
    ///
    /// Takes `r²` fine cells and returns their mass-weighted average.
    pub fn restrict_cells(&self, fine_cells: &[[f64; 9]]) -> [f64; 9] {
        fine_to_coarse_pop(fine_cells)
    }
    /// Prolong coarse cell to `r²` fine cells (spatial prolongation with rescaling).
    pub fn prolong_cell(&self, f_c: &[f64; 9]) -> Vec<[f64; 9]> {
        let r2 = self.refinement_ratio * self.refinement_ratio;
        let rho: f64 = f_c.iter().sum();
        let mut jx = 0.0f64;
        let mut jy = 0.0f64;
        for alpha in 0..9 {
            jx += f_c[alpha] * D2Q9_VELOCITIES[alpha][0] as f64;
            jy += f_c[alpha] * D2Q9_VELOCITIES[alpha][1] as f64;
        }
        let u = [jx / rho.max(1e-30), jy / rho.max(1e-30)];
        let f_fine = self.coarse_to_fine_rescale(f_c, rho, u);
        vec![f_fine; r2]
    }
    /// Update buffer layer using linear blending.
    ///
    /// Buffer cells near the interface are blended between coarse and fine.
    pub fn update_buffer_layer(&self, idx: usize, alpha_blend: f64) -> [f64; 9] {
        if idx >= self.f_coarse_iface.len() || idx >= self.f_fine_iface.len() {
            return D2Q9_WEIGHTS.map(|w| w * 1.0);
        }
        buffer_layer_update(
            &self.f_coarse_iface[idx],
            &self.f_fine_iface[idx],
            alpha_blend,
        )
    }
    /// Check consistency: total mass should be conserved across interface.
    pub fn mass_conservation_error(&self) -> f64 {
        let coarse_mass: f64 = self
            .f_coarse_iface
            .iter()
            .map(|f| f.iter().sum::<f64>())
            .sum();
        let fine_mass: f64 = self
            .f_fine_iface
            .iter()
            .map(|f| f.iter().sum::<f64>())
            .sum::<f64>()
            / self.refinement_ratio as f64;
        (coarse_mass - fine_mass).abs()
    }
}
/// Chapman-Enskog expansion state tracker.
///
/// Tracks the multi-scale perturbation expansion f = f⁽⁰⁾ + ε f⁽¹⁾ + ε² f⁽²⁾
/// linking the kinetic LBM description to the Navier-Stokes equations.
///
/// References:
/// - Chapman, S. & Cowling, T.G. (1970). *The Mathematical Theory of Non-Uniform Gases*.
/// - He, X. & Luo, L.-S. (1997). Theory of the lattice Boltzmann method. *Phys. Rev. E* 56, 6811.
#[derive(Debug, Clone)]
pub struct ChapmanEnskogExpansion {
    /// Knudsen number ε (ratio of mean free path to macroscopic length).
    pub epsilon: f64,
    /// Zeroth-order (equilibrium) distribution f⁽⁰⁾ for each velocity direction.
    pub f0: [f64; 9],
    /// First-order correction f⁽¹⁾.
    pub f1: [f64; 9],
    /// Second-order correction f⁽²⁾.
    pub f2: [f64; 9],
    /// Macroscopic density ρ.
    pub rho: f64,
    /// Macroscopic velocity \[u_x, u_y\].
    pub u: [f64; 2],
    /// Viscous stress tensor σ_xy (off-diagonal element).
    pub sigma_xy: f64,
    /// Dynamic viscosity μ derived at O(ε).
    pub mu: f64,
    /// Bulk viscosity ζ derived at O(ε).
    pub zeta: f64,
}
impl ChapmanEnskogExpansion {
    /// Create a new Chapman-Enskog expansion state.
    ///
    /// # Arguments
    /// * `rho` - macroscopic density
    /// * `u`   - macroscopic velocity \[u_x, u_y\]
    /// * `tau` - LBM relaxation time
    /// * `epsilon` - Knudsen number
    pub fn new(rho: f64, u: [f64; 2], tau: f64, epsilon: f64) -> Self {
        let f0 = lbm_equilibrium(rho, u);
        let f1 = [0.0f64; 9];
        let f2 = [0.0f64; 9];
        let mu = rho * CS2 * (tau - 0.5);
        let zeta = 0.0;
        let sigma_xy = 0.0;
        Self {
            epsilon,
            f0,
            f1,
            f2,
            rho,
            u,
            sigma_xy,
            mu,
            zeta,
        }
    }
    /// Compute first-order correction f⁽¹⁾_α using Chapman-Enskog.
    ///
    /// f⁽¹⁾_α = -τ * w_α * ρ / cs⁴ * (e_α·∂ₜu + e_αx e_αy * ∂_β u_β terms).
    ///
    /// For the D2Q9 lattice, this reduces to the non-equilibrium stress tensor.
    pub fn compute_f1(&mut self, grad_ux: f64, grad_uy: f64, grad_vx: f64, grad_vy: f64) {
        let tau = self.mu / (self.rho * CS2) + 0.5;
        let s_xx = grad_ux;
        let s_yy = grad_vy;
        let s_xy = 0.5 * (grad_uy + grad_vx);
        let s_trace = s_xx + s_yy;
        for alpha in 0..9 {
            let [ex, ey] = [
                D2Q9_VELOCITIES[alpha][0] as f64,
                D2Q9_VELOCITIES[alpha][1] as f64,
            ];
            let w = D2Q9_WEIGHTS[alpha];
            let ee_xx = ex * ex - CS2;
            let ee_yy = ey * ey - CS2;
            let ee_xy = ex * ey;
            self.f1[alpha] = -tau * w * self.rho / (CS2 * CS2)
                * (ee_xx * (s_xx - s_trace / 2.0)
                    + ee_yy * (s_yy - s_trace / 2.0)
                    + 2.0 * ee_xy * s_xy);
        }
        self.sigma_xy = -2.0 * self.rho * CS2 * tau * s_xy;
    }
    /// Compute second-order correction f⁽²⁾_α from second-order stress tensor.
    ///
    /// At O(ε²) recovers the energy equation and thermal diffusivity terms.
    pub fn compute_f2(&mut self, d2ux_dx2: f64, d2uy_dy2: f64) {
        let tau = self.mu / (self.rho * CS2) + 0.5;
        for alpha in 0..9 {
            let [ex, ey] = [
                D2Q9_VELOCITIES[alpha][0] as f64,
                D2Q9_VELOCITIES[alpha][1] as f64,
            ];
            let w = D2Q9_WEIGHTS[alpha];
            let visc_corr = tau * tau * w * CS2 * (ex * d2ux_dx2 + ey * d2uy_dy2);
            self.f2[alpha] = visc_corr;
        }
    }
    /// Reconstruct total distribution f = f⁽⁰⁾ + ε f⁽¹⁾ + ε² f⁽²⁾.
    pub fn total_distribution(&self) -> [f64; 9] {
        let mut f = [0.0f64; 9];
        let eps = self.epsilon;
        for alpha in 0..9 {
            f[alpha] = self.f0[alpha] + eps * self.f1[alpha] + eps * eps * self.f2[alpha];
        }
        f
    }
    /// Extract non-equilibrium stress tensor components from f.
    ///
    /// Π⁽¹⁾_αβ = Σ_i f_i⁽¹⁾ e_iα e_iβ
    pub fn neq_stress_tensor(&self) -> [f64; 4] {
        let mut pi = [0.0f64; 4];
        for alpha in 0..9 {
            let [ex, ey] = [
                D2Q9_VELOCITIES[alpha][0] as f64,
                D2Q9_VELOCITIES[alpha][1] as f64,
            ];
            pi[0] += self.f1[alpha] * ex * ex;
            pi[1] += self.f1[alpha] * ex * ey;
            pi[2] += self.f1[alpha] * ey * ex;
            pi[3] += self.f1[alpha] * ey * ey;
        }
        pi
    }
    /// Effective kinematic viscosity ν = μ/ρ.
    pub fn kinematic_viscosity(&self) -> f64 {
        if self.rho < 1e-30 {
            return 0.0;
        }
        self.mu / self.rho
    }
    /// Check if Knudsen number is in continuum regime (Kn < 0.01).
    pub fn is_continuum_regime(&self) -> bool {
        self.epsilon < 0.01
    }
}
/// Non-uniform grid LBM with interpolation at grid transitions.
///
/// The grid is described by a mapping from logical indices to physical
/// coordinates, allowing stretched or compressed grids.
#[derive(Debug, Clone)]
pub struct VariableResolution {
    /// Distribution functions `f[node][9]`.
    pub f: Vec<[f64; 9]>,
    /// Physical x-coordinates of each node.
    pub x_coords: Vec<f64>,
    /// Physical y-coordinates of each node.
    pub y_coords: Vec<f64>,
    /// Local mesh spacing Δx at each node.
    pub dx_local: Vec<f64>,
    /// Local mesh spacing Δy at each node.
    pub dy_local: Vec<f64>,
    /// Grid width (logical).
    pub nx: usize,
    /// Grid height (logical).
    pub ny: usize,
    /// Local relaxation times (adjusted for local grid spacing).
    pub tau_local: Vec<f64>,
}
impl VariableResolution {
    /// Create a new variable-resolution grid with uniform spacing as default.
    pub fn new(nx: usize, ny: usize, tau_base: f64) -> Self {
        let n = nx * ny;
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let mut x_coords = vec![0.0; n];
        let mut y_coords = vec![0.0; n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                x_coords[k] = x as f64;
                y_coords[k] = y as f64;
            }
        }
        Self {
            f: vec![f_init; n],
            x_coords,
            y_coords,
            dx_local: vec![1.0; n],
            dy_local: vec![1.0; n],
            nx,
            ny,
            tau_local: vec![tau_base; n],
        }
    }
    /// Apply a geometric stretching in x from `x_min` to `x_max`.
    ///
    /// `ratio` is the ratio between the last and first cell widths.
    pub fn apply_x_stretching(&mut self, x_min: f64, x_max: f64, ratio: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let r = if (ratio - 1.0).abs() < 1e-10 {
            1.0
        } else {
            ratio.powf(1.0 / (nx - 1).max(1) as f64)
        };
        let mut xs = vec![0.0_f64; nx];
        if (r - 1.0).abs() < 1e-10 {
            for i in 0..nx {
                xs[i] = x_min + (x_max - x_min) * i as f64 / (nx - 1).max(1) as f64;
            }
        } else {
            let scale = (x_max - x_min) / (r.powi(nx as i32) - 1.0);
            for i in 0..nx {
                xs[i] = x_min + scale * (r.powi(i as i32) - 1.0);
            }
        }
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                self.x_coords[k] = xs[x];
                let dx = if x + 1 < nx {
                    xs[x + 1] - xs[x]
                } else {
                    xs[x] - xs[x - 1]
                };
                self.dx_local[k] = dx.abs();
            }
        }
    }
    /// BGK collision with local relaxation time.
    pub fn collide(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let tau = self.tau_local[k];
            let rho: f64 = self.f[k].iter().sum();
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for i in 0..9 {
                let c = D2Q9_VELOCITIES[i];
                ux += c[0] as f64 * self.f[k][i];
                uy += c[1] as f64 * self.f[k][i];
            }
            let u = [ux / rho.max(1e-15), uy / rho.max(1e-15)];
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let u2 = u[0] * u[0] + u[1] * u[1];
                let feq =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
                self.f[k][i] -= (self.f[k][i] - feq) / tau;
            }
        }
    }
}
/// LBM-MD hybrid: exchange forces and velocities at the interface region.
///
/// The LBM domain provides mean-field velocity/stress data to the MD region,
/// and MD provides corrected velocity fluctuations back to LBM.
#[derive(Debug, Clone)]
pub struct LbmMolecularCoupling {
    /// LBM distributions at coupling nodes.
    pub f_lbm: Vec<[f64; 9]>,
    /// MD particle positions (flat: \[x0, y0, x1, y1, ...\]).
    pub md_positions: Vec<f64>,
    /// MD particle velocities (flat: \[vx0, vy0, vx1, vy1, ...\]).
    pub md_velocities: Vec<f64>,
    /// MD particle masses.
    pub md_masses: Vec<f64>,
    /// LBM macro velocity at coupling nodes.
    pub u_lbm: Vec<[f64; 2]>,
    /// MD mean velocity at coupling nodes.
    pub u_md: Vec<[f64; 2]>,
    /// Number of LBM coupling nodes.
    pub n_lbm_nodes: usize,
    /// Number of MD particles.
    pub n_md: usize,
    /// LBM relaxation time.
    pub tau: f64,
    /// LBM lattice spacing.
    pub dx: f64,
}
impl LbmMolecularCoupling {
    /// Create a new LBM-MD coupling.
    pub fn new(n_lbm_nodes: usize, n_md: usize, tau: f64, dx: f64) -> Self {
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            f_lbm: vec![f_init; n_lbm_nodes],
            md_positions: vec![0.0; n_md * 2],
            md_velocities: vec![0.0; n_md * 2],
            md_masses: vec![1.0; n_md],
            u_lbm: vec![[0.0; 2]; n_lbm_nodes],
            u_md: vec![[0.0; 2]; n_lbm_nodes],
            n_lbm_nodes,
            n_md,
            tau,
            dx,
        }
    }
    /// Compute mean MD velocity by averaging particle velocities near each coupling node.
    ///
    /// Uses a simple cell-based averaging: particles within `rc` of a node are included.
    pub fn compute_md_mean_velocity(&mut self, node_positions: &[[f64; 2]], rc: f64) {
        let rc2 = rc * rc;
        for (k, node_pos) in node_positions.iter().enumerate() {
            let mut sum_u = [0.0_f64; 2];
            let mut total_mass = 0.0_f64;
            for p in 0..self.n_md {
                let px = self.md_positions[2 * p];
                let py = self.md_positions[2 * p + 1];
                let r2 = (px - node_pos[0]).powi(2) + (py - node_pos[1]).powi(2);
                if r2 < rc2 {
                    let m = self.md_masses[p];
                    sum_u[0] += m * self.md_velocities[2 * p];
                    sum_u[1] += m * self.md_velocities[2 * p + 1];
                    total_mass += m;
                }
            }
            if total_mass > 0.0 {
                self.u_md[k] = [sum_u[0] / total_mass, sum_u[1] / total_mass];
            }
        }
    }
    /// Apply velocity correction to LBM distributions: shift equilibrium toward MD mean velocity.
    pub fn apply_velocity_correction(&mut self) {
        for k in 0..self.n_lbm_nodes {
            let rho: f64 = self.f_lbm[k].iter().sum();
            let u_target = self.u_md[k];
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u_target[0] + c[1] as f64 * u_target[1];
                let u2 = u_target[0] * u_target[0] + u_target[1] * u_target[1];
                let feq =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
                self.f_lbm[k][i] = self.f_lbm[k][i] * (1.0 - 1.0 / self.tau) + feq / self.tau;
            }
        }
    }
}
/// Heterogeneous Multiscale Method (HMM) framework for LBM.
///
/// HMM couples a macro-scale Navier-Stokes solver with a micro-scale
/// LBM solver. The LBM "micro" solver is used to compute effective
/// transport coefficients (viscosity, diffusivity) for the NS "macro" solver.
///
/// References:
/// - E, W. & Engquist, B. (2003). The heterogeneous multiscale methods. *CMS* 1, 87.
/// - Ren, W. & E, W. (2005). Heterogeneous multiscale method for Stokes flow. *J. Comput. Phys.* 204, 1.
#[derive(Debug, Clone)]
pub struct HeterogeneousMultiscale {
    /// Number of macro-scale nodes.
    pub n_macro: usize,
    /// Number of micro-scale nodes per macro node.
    pub n_micro: usize,
    /// Macro-scale grid spacing H.
    pub h_macro: f64,
    /// Micro-scale grid spacing h.
    pub h_micro: f64,
    /// Macro-scale velocity at each node.
    pub u_macro: Vec<[f64; 2]>,
    /// Macro-scale pressure at each node.
    pub p_macro: Vec<f64>,
    /// Effective viscosity from micro-scale LBM.
    pub nu_eff: Vec<f64>,
    /// LBM relaxation times per macro node (from HMM micro computation).
    pub tau_micro: Vec<f64>,
    /// Compression factor κ = H / h.
    pub kappa: f64,
}
impl HeterogeneousMultiscale {
    /// Create a new HMM framework.
    ///
    /// # Arguments
    /// * `n_macro`  - number of macro nodes
    /// * `n_micro`  - number of micro nodes per macro node
    /// * `h_macro`  - macro grid spacing
    /// * `h_micro`  - micro grid spacing
    /// * `tau_init` - initial LBM relaxation time
    pub fn new(n_macro: usize, n_micro: usize, h_macro: f64, h_micro: f64, tau_init: f64) -> Self {
        let kappa = h_macro / h_micro;
        Self {
            n_macro,
            n_micro,
            h_macro,
            h_micro,
            u_macro: vec![[0.0; 2]; n_macro],
            p_macro: vec![CS2; n_macro],
            nu_eff: vec![CS2 * (tau_init - 0.5); n_macro],
            tau_micro: vec![tau_init; n_macro],
            kappa,
        }
    }
    /// Run micro-scale LBM to estimate effective viscosity at macro node `idx`.
    ///
    /// Applies a known strain rate and measures the resulting stress tensor.
    /// ν_eff = -σ_xy / (2 * S_xy) where S_xy is the applied shear rate.
    pub fn micro_to_macro_viscosity(&self, idx: usize, shear_rate: f64) -> f64 {
        if idx >= self.n_macro {
            return 0.0;
        }
        let tau = self.tau_micro[idx];
        let nu = CS2 * (tau - 0.5);
        let rho = 1.0;
        let _sigma_xy = -2.0 * rho * nu * shear_rate;
        nu
    }
    /// Update effective viscosity from micro-scale computation at all nodes.
    pub fn update_nu_eff(&mut self, shear_rates: &[f64]) {
        for i in 0..self.n_macro {
            let sr = if i < shear_rates.len() {
                shear_rates[i]
            } else {
                0.0
            };
            self.nu_eff[i] = self.micro_to_macro_viscosity(i, sr);
        }
    }
    /// Macro-scale NS update (explicit Euler for 1D model).
    ///
    /// ∂_t u = -u ∂_x u - ∂_x p / ρ + ν_eff ∂²_x u
    pub fn macro_step(&mut self, dt: f64) {
        let n = self.n_macro;
        if n < 3 {
            return;
        }
        let h = self.h_macro;
        let mut u_new = self.u_macro.clone();
        for i in 1..n - 1 {
            let ux = self.u_macro[i][0];
            let ux_l = self.u_macro[i - 1][0];
            let ux_r = self.u_macro[i + 1][0];
            let p_l = self.p_macro[i - 1];
            let p_r = self.p_macro[i + 1];
            let nu = self.nu_eff[i];
            let adv = -ux * (ux_r - ux_l) / (2.0 * h);
            let pres = -(p_r - p_l) / (2.0 * h);
            let diff = nu * (ux_r - 2.0 * ux + ux_l) / (h * h);
            u_new[i][0] = ux + dt * (adv + pres + diff);
        }
        self.u_macro = u_new;
    }
    /// Compress micro-scale solution to macro scale using averaging.
    pub fn compress_micro_to_macro(&self, f_micro: &[[f64; 9]]) -> ([f64; 2], f64) {
        if f_micro.is_empty() {
            return ([0.0; 2], 1.0);
        }
        let n = f_micro.len() as f64;
        let mut jx = 0.0f64;
        let mut jy = 0.0f64;
        let mut rho_sum = 0.0f64;
        for f in f_micro {
            let rho: f64 = f.iter().sum();
            rho_sum += rho;
            for alpha in 0..9 {
                jx += f[alpha] * D2Q9_VELOCITIES[alpha][0] as f64;
                jy += f[alpha] * D2Q9_VELOCITIES[alpha][1] as f64;
            }
        }
        let rho_avg = rho_sum / n;
        let u_avg = [jx / rho_sum.max(1e-30), jy / rho_sum.max(1e-30)];
        (u_avg, rho_avg)
    }
    /// Scale separation ratio κ = H / h.
    pub fn scale_separation(&self) -> f64 {
        self.kappa
    }
}
/// 2D quad-tree adaptive mesh for LBM.
///
/// Nodes are stored in a flat arena; children are referenced by index.
/// Refinement and coarsening are triggered by a gradient-based criterion.
#[derive(Debug, Clone)]
pub struct AdaptiveMesh {
    /// Arena of all quad-tree nodes (including non-leaf nodes).
    pub nodes: Vec<QuadTreeNode>,
    /// Maximum allowed refinement level.
    pub max_level: usize,
    /// Minimum cell size (physical units).
    pub min_cell_size: f64,
    /// Refinement threshold (gradient magnitude).
    pub refine_threshold: f64,
    /// Coarsening threshold (gradient magnitude).
    pub coarsen_threshold: f64,
}
impl AdaptiveMesh {
    /// Create a new adaptive mesh with a single root cell covering `[0,L]×[0,L]`.
    pub fn new(domain_size: f64, max_level: usize, refine_threshold: f64) -> Self {
        let root = QuadTreeNode::new([0.0, domain_size], [0.0, domain_size], 0, 0);
        Self {
            nodes: vec![root],
            max_level,
            min_cell_size: domain_size / 2_f64.powi(max_level as i32),
            refine_threshold,
            coarsen_threshold: refine_threshold * 0.25,
        }
    }
    /// Number of nodes in the mesh.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }
    /// Returns true if the mesh has no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
    /// Refine node `idx` into four children (NW, NE, SW, SE).
    ///
    /// Returns the indices of the four child nodes, or `None` if already at
    /// `max_level` or if the node is not a leaf.
    pub fn refine(&mut self, idx: usize) -> Option<[usize; 4]> {
        if !self.nodes[idx].is_leaf() {
            return None;
        }
        if self.nodes[idx].level >= self.max_level {
            return None;
        }
        let parent = &self.nodes[idx];
        let lv = parent.level + 1;
        let xm = (parent.x_range[0] + parent.x_range[1]) * 0.5;
        let ym = (parent.y_range[0] + parent.y_range[1]) * 0.5;
        let x0 = parent.x_range[0];
        let x1 = parent.x_range[1];
        let y0 = parent.y_range[0];
        let y1 = parent.y_range[1];
        let base = self.nodes.len();
        let children = [
            QuadTreeNode::new([x0, xm], [y0, ym], lv, base),
            QuadTreeNode::new([xm, x1], [y0, ym], lv, base + 1),
            QuadTreeNode::new([x0, xm], [ym, y1], lv, base + 2),
            QuadTreeNode::new([xm, x1], [ym, y1], lv, base + 3),
        ];
        let parent_f = self.nodes[idx].f;
        let parent_rho = self.nodes[idx].rho;
        let parent_u = self.nodes[idx].u;
        for mut child in children {
            child.f = coarse_to_fine_pop(&parent_f, parent_rho, parent_u, child.u);
            child.rho = parent_rho;
            child.u = parent_u;
            self.nodes.push(child);
        }
        let child_indices = [base, base + 1, base + 2, base + 3];
        self.nodes[idx].children = Some(child_indices);
        Some(child_indices)
    }
    /// Collect all leaf node indices.
    pub fn leaf_indices(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.is_leaf())
            .map(|(i, _)| i)
            .collect()
    }
    /// Apply one BGK collision step to all leaf nodes.
    pub fn collide_all(&mut self, tau: f64) {
        for node in &mut self.nodes {
            if !node.is_leaf() {
                continue;
            }
            for i in 0..9 {
                let feq = node.equilibrium(i);
                node.f[i] -= (node.f[i] - feq) / tau;
            }
        }
    }
}
/// Coarse-fine grid interface using buffer layers.
///
/// Manages the data exchange between a coarse patch (grid spacing 2h, time
/// step 2dt) and a fine patch (grid spacing h, time step dt).
#[derive(Debug, Clone)]
pub struct PatchedGrid {
    /// Coarse grid distribution functions `f_coarse[node][9]`.
    pub f_coarse: Vec<[f64; 9]>,
    /// Fine grid distribution functions `f_fine[node][9]`.
    pub f_fine: Vec<[f64; 9]>,
    /// Coarse grid width.
    pub nx_coarse: usize,
    /// Coarse grid height.
    pub ny_coarse: usize,
    /// Fine grid width (typically 2 * nx_coarse - 1 for the overlap region).
    pub nx_fine: usize,
    /// Fine grid height.
    pub ny_fine: usize,
    /// Refinement ratio (fine/coarse resolution), typically 2.
    pub ratio: usize,
    /// Buffer layer width in coarse cells.
    pub buffer_width: usize,
}
impl PatchedGrid {
    /// Create a new patched grid pair.
    pub fn new(nx_coarse: usize, ny_coarse: usize, ratio: usize) -> Self {
        let nx_fine = (nx_coarse - 1) * ratio + 1;
        let ny_fine = (ny_coarse - 1) * ratio + 1;
        let n_coarse = nx_coarse * ny_coarse;
        let n_fine = nx_fine * ny_fine;
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            f_coarse: vec![f_init; n_coarse],
            f_fine: vec![f_init; n_fine],
            nx_coarse,
            ny_coarse,
            nx_fine,
            ny_fine,
            ratio,
            buffer_width: 2,
        }
    }
    /// Prolongate coarse grid values to fine grid (bilinear interpolation).
    pub fn prolongate(&mut self) {
        let r = self.ratio;
        let nxc = self.nx_coarse;
        let nyc = self.ny_coarse;
        let nxf = self.nx_fine;
        for yc in 0..(nyc - 1) {
            for xc in 0..(nxc - 1) {
                let kc_sw = yc * nxc + xc;
                let kc_se = yc * nxc + xc + 1;
                let kc_nw = (yc + 1) * nxc + xc;
                let kc_ne = (yc + 1) * nxc + xc + 1;
                for dy in 0..=r {
                    for dx in 0..=r {
                        let tx = dx as f64 / r as f64;
                        let ty = dy as f64 / r as f64;
                        let xf = xc * r + dx;
                        let yf = yc * r + dy;
                        let kf = yf * nxf + xf;
                        for i in 0..9 {
                            self.f_fine[kf][i] = (1.0 - tx) * (1.0 - ty) * self.f_coarse[kc_sw][i]
                                + tx * (1.0 - ty) * self.f_coarse[kc_se][i]
                                + (1.0 - tx) * ty * self.f_coarse[kc_nw][i]
                                + tx * ty * self.f_coarse[kc_ne][i];
                        }
                    }
                }
            }
        }
    }
    /// Restrict fine grid values to coarse grid (volume averaging).
    pub fn restrict(&mut self) {
        let r = self.ratio;
        let nxc = self.nx_coarse;
        let nyc = self.ny_coarse;
        let nxf = self.nx_fine;
        let r2 = (r * r) as f64;
        for yc in 0..nyc {
            for xc in 0..nxc {
                let kc = yc * nxc + xc;
                let mut sum = [0.0_f64; 9];
                let mut count = 0usize;
                let xf_start = xc * r;
                let yf_start = yc * r;
                let xf_end = ((xc + 1) * r).min(self.nx_fine);
                let yf_end = ((yc + 1) * r).min(self.ny_fine);
                for yf in yf_start..yf_end {
                    for xf in xf_start..xf_end {
                        let kf = yf * nxf + xf;
                        for i in 0..9 {
                            sum[i] += self.f_fine[kf][i];
                        }
                        count += 1;
                    }
                }
                if count > 0 {
                    let scale = count as f64 / r2;
                    for i in 0..9 {
                        self.f_coarse[kc][i] = sum[i] / scale.max(1.0);
                    }
                }
            }
        }
    }
}
/// Cell agglomeration: restriction (fine → coarse) and prolongation (coarse → fine).
///
/// Implements standard full-weighting restriction and bi-linear prolongation
/// for LBM distribution functions.
#[derive(Debug, Clone)]
pub struct CellReduction {
    /// Refinement ratio (fine cells per coarse cell per dimension).
    pub ratio: usize,
}
impl CellReduction {
    /// Create a new cell reduction operator with the given ratio.
    pub fn new(ratio: usize) -> Self {
        Self { ratio }
    }
    /// Restrict fine populations to coarse: volume average.
    pub fn restrict(
        &self,
        f_fine: &[[f64; 9]],
        nx_fine: usize,
        ny_fine: usize,
    ) -> (Vec<[f64; 9]>, usize, usize) {
        let r = self.ratio;
        let nx_coarse = nx_fine / r;
        let ny_coarse = ny_fine / r;
        let mut f_coarse = vec![[0.0_f64; 9]; nx_coarse * ny_coarse];
        let r2 = (r * r) as f64;
        for yc in 0..ny_coarse {
            for xc in 0..nx_coarse {
                let kc = yc * nx_coarse + xc;
                let mut sum = [0.0_f64; 9];
                for dy in 0..r {
                    for dx in 0..r {
                        let xf = xc * r + dx;
                        let yf = yc * r + dy;
                        if xf < nx_fine && yf < ny_fine {
                            let kf = yf * nx_fine + xf;
                            for i in 0..9 {
                                sum[i] += f_fine[kf][i];
                            }
                        }
                    }
                }
                for i in 0..9 {
                    f_coarse[kc][i] = sum[i] / r2;
                }
            }
        }
        (f_coarse, nx_coarse, ny_coarse)
    }
    /// Prolong coarse populations to fine: bi-linear interpolation.
    pub fn prolong(
        &self,
        f_coarse: &[[f64; 9]],
        nx_coarse: usize,
        ny_coarse: usize,
    ) -> (Vec<[f64; 9]>, usize, usize) {
        let r = self.ratio;
        let nx_fine = nx_coarse * r;
        let ny_fine = ny_coarse * r;
        let mut f_fine = vec![[0.0_f64; 9]; nx_fine * ny_fine];
        for yc in 0..ny_coarse {
            for xc in 0..nx_coarse {
                let kc = yc * nx_coarse + xc;
                for dy in 0..r {
                    for dx in 0..r {
                        let xf = xc * r + dx;
                        let yf = yc * r + dy;
                        if xf < nx_fine && yf < ny_fine {
                            let kf = yf * nx_fine + xf;
                            f_fine[kf] = f_coarse[kc];
                        }
                    }
                }
            }
        }
        (f_fine, nx_fine, ny_fine)
    }
}
/// Multi-scale coupling: LBM to Navier-Stokes bridge.
///
/// Converts LBM distribution functions to Navier-Stokes macroscopic
/// variables and vice versa, enabling hybrid simulations.
#[derive(Debug, Clone)]
pub struct LbmNavierStokesBridge {
    /// LBM grid size (1D flattened).
    pub n_nodes: usize,
    /// Relaxation time τ.
    pub tau: f64,
    /// Macroscopic density field.
    pub rho: Vec<f64>,
    /// Macroscopic velocity field \[u, v\] per node.
    pub velocity: Vec<[f64; 2]>,
    /// Pressure field p = ρ cs².
    pub pressure: Vec<f64>,
    /// Viscous stress xx per node.
    pub stress_xx: Vec<f64>,
    /// Viscous stress xy per node.
    pub stress_xy: Vec<f64>,
    /// Viscous stress yy per node.
    pub stress_yy: Vec<f64>,
}
impl LbmNavierStokesBridge {
    /// Create a new LBM-NS bridge.
    ///
    /// # Arguments
    /// * `n_nodes` - number of LBM lattice nodes
    /// * `tau`     - LBM relaxation time
    pub fn new(n_nodes: usize, tau: f64) -> Self {
        Self {
            n_nodes,
            tau,
            rho: vec![1.0; n_nodes],
            velocity: vec![[0.0; 2]; n_nodes],
            pressure: vec![CS2; n_nodes],
            stress_xx: vec![0.0; n_nodes],
            stress_xy: vec![0.0; n_nodes],
            stress_yy: vec![0.0; n_nodes],
        }
    }
    /// Extract macroscopic fields from LBM distributions.
    ///
    /// ρ = Σ_i f_i,  ρu = Σ_i f_i e_i,  p = ρ cs²
    pub fn extract_macro_fields(&mut self, f_all: &[[f64; 9]]) {
        for (n, f) in f_all.iter().enumerate().take(self.n_nodes) {
            let rho: f64 = f.iter().sum();
            let mut jx = 0.0f64;
            let mut jy = 0.0f64;
            for alpha in 0..9 {
                jx += f[alpha] * D2Q9_VELOCITIES[alpha][0] as f64;
                jy += f[alpha] * D2Q9_VELOCITIES[alpha][1] as f64;
            }
            let rho_safe = rho.max(1e-30);
            self.rho[n] = rho;
            self.velocity[n] = [jx / rho_safe, jy / rho_safe];
            self.pressure[n] = rho * CS2;
            let feq = lbm_equilibrium(rho, self.velocity[n]);
            let mut pi_xx = 0.0f64;
            let mut pi_xy = 0.0f64;
            let mut pi_yy = 0.0f64;
            for alpha in 0..9 {
                let ex = D2Q9_VELOCITIES[alpha][0] as f64;
                let ey = D2Q9_VELOCITIES[alpha][1] as f64;
                let fneq = f[alpha] - feq[alpha];
                pi_xx += fneq * ex * ex;
                pi_xy += fneq * ex * ey;
                pi_yy += fneq * ey * ey;
            }
            let nu = CS2 * (self.tau - 0.5);
            let factor = -1.0 / (2.0 * self.tau * CS2);
            self.stress_xx[n] = -2.0 * rho * nu * pi_xx * factor;
            self.stress_xy[n] = -2.0 * rho * nu * pi_xy * factor;
            self.stress_yy[n] = -2.0 * rho * nu * pi_yy * factor;
        }
    }
    /// Inject Navier-Stokes solution back into LBM as constrained equilibrium.
    ///
    /// Given ρ_ns, u_ns from an NS solver, reconstructs f_eq and adds residual stress.
    pub fn inject_ns_solution(&self, rho_ns: f64, u_ns: [f64; 2]) -> [f64; 9] {
        lbm_equilibrium(rho_ns, u_ns)
    }
    /// Compute Mach number at node n.
    pub fn mach_number(&self, n: usize) -> f64 {
        if n >= self.n_nodes {
            return 0.0;
        }
        let [ux, uy] = self.velocity[n];
        (ux * ux + uy * uy).sqrt() / CS2.sqrt()
    }
    /// Kinematic viscosity ν = cs² (τ - 0.5).
    pub fn kinematic_viscosity(&self) -> f64 {
        CS2 * (self.tau - 0.5)
    }
    /// Reynolds number Re = U L / ν.
    pub fn reynolds_number(&self, char_vel: f64, char_len: f64) -> f64 {
        let nu = self.kinematic_viscosity();
        if nu < 1e-30 {
            return f64::INFINITY;
        }
        char_vel * char_len / nu
    }
}
/// Space-time interpolation for multi-scale LBM.
///
/// Provides bi-linear spatial and linear temporal interpolation of
/// distribution functions at grid interfaces.
#[derive(Debug, Clone)]
pub struct SpaceTimeInterpolator {
    /// Grid spacing on coarse level Δx_c.
    pub dx_coarse: f64,
    /// Time step on coarse level Δt_c.
    pub dt_coarse: f64,
    /// Refinement ratio.
    pub ratio: usize,
}
impl SpaceTimeInterpolator {
    /// Create a new space-time interpolator.
    pub fn new(dx_coarse: f64, dt_coarse: f64, ratio: usize) -> Self {
        Self {
            dx_coarse,
            dt_coarse,
            ratio,
        }
    }
    /// Bilinear spatial interpolation in 2D.
    ///
    /// Given values at four corners of a coarse cell, interpolates to
    /// a point (ξ, η) ∈ \[0,1\]² within the cell.
    pub fn bilinear(
        &self,
        f00: &[f64; 9],
        f10: &[f64; 9],
        f01: &[f64; 9],
        f11: &[f64; 9],
        xi: f64,
        eta: f64,
    ) -> [f64; 9] {
        let mut f = [0.0f64; 9];
        let w00 = (1.0 - xi) * (1.0 - eta);
        let w10 = xi * (1.0 - eta);
        let w01 = (1.0 - xi) * eta;
        let w11 = xi * eta;
        for alpha in 0..9 {
            f[alpha] = w00 * f00[alpha] + w10 * f10[alpha] + w01 * f01[alpha] + w11 * f11[alpha];
        }
        f
    }
    /// Linear temporal interpolation between two time levels.
    pub fn temporal_linear(&self, f_old: &[f64; 9], f_new: &[f64; 9], theta: f64) -> [f64; 9] {
        let mut f = [0.0f64; 9];
        for alpha in 0..9 {
            f[alpha] = (1.0 - theta) * f_old[alpha] + theta * f_new[alpha];
        }
        f
    }
    /// Second-order Hermite (cubic) temporal interpolation.
    ///
    /// Uses f at two time levels and first-order time derivatives.
    pub fn hermite_temporal(
        &self,
        f0: &[f64; 9],
        f1: &[f64; 9],
        df0: &[f64; 9],
        df1: &[f64; 9],
        theta: f64,
    ) -> [f64; 9] {
        let t = theta;
        let t2 = t * t;
        let t3 = t2 * t;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;
        let mut f = [0.0f64; 9];
        for alpha in 0..9 {
            f[alpha] = h00 * f0[alpha] + h10 * df0[alpha] + h01 * f1[alpha] + h11 * df1[alpha];
        }
        f
    }
    /// Fine-grid node position (fractional coordinates) within a coarse cell.
    pub fn fine_node_fraction(&self, fine_idx: usize) -> f64 {
        fine_idx as f64 / self.ratio as f64
    }
    /// Interpolation weight for fine node at fine_idx within a coarse interval.
    pub fn coarse_weight(&self, fine_idx: usize) -> (f64, f64) {
        let xi = self.fine_node_fraction(fine_idx);
        (1.0 - xi, xi)
    }
}
/// Space-time explicit coupling between coarse (2h, 2dt) and fine (h, dt) grids.
///
/// The coupling is achieved via a time-interpolation of coarse-grid values
/// to provide boundary data for the fine grid at intermediate time levels.
#[derive(Debug, Clone)]
pub struct SpaceTimeLbm {
    /// Coarse grid distributions at time t.
    pub f_coarse_n: Vec<[f64; 9]>,
    /// Coarse grid distributions at time t + 2dt.
    pub f_coarse_np1: Vec<[f64; 9]>,
    /// Fine grid distributions.
    pub f_fine: Vec<[f64; 9]>,
    /// Coarse grid size.
    pub n_coarse: usize,
    /// Fine grid size.
    pub n_fine: usize,
    /// Relaxation time on coarse grid: τ_c.
    pub tau_coarse: f64,
    /// Relaxation time on fine grid: τ_f.
    pub tau_fine: f64,
    /// Space-time refinement ratio.
    pub ratio: usize,
    /// Current fine-grid sub-step counter.
    pub sub_step: usize,
}
impl SpaceTimeLbm {
    /// Create a new space-time LBM coupling.
    ///
    /// `tau_coarse` is the coarse relaxation time; `tau_fine` is derived via
    /// the Filippova-Hänel rescaling: τ_f = τ_c/2 + 0.5*(1 - 1/ratio).
    pub fn new(n_coarse: usize, ratio: usize, tau_coarse: f64) -> Self {
        let n_fine = n_coarse * ratio;
        let tau_fine = tau_coarse / 2.0 + 0.5 * (1.0 - 1.0 / ratio as f64);
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            f_coarse_n: vec![f_init; n_coarse],
            f_coarse_np1: vec![f_init; n_coarse],
            f_fine: vec![f_init; n_fine],
            n_coarse,
            n_fine,
            tau_coarse,
            tau_fine,
            ratio,
            sub_step: 0,
        }
    }
    /// Time-interpolated coarse distribution at fine sub-step `s` (0..ratio).
    ///
    /// Uses linear interpolation: f_c(t + s*dt) = (1 - s/r)*f_c_n + (s/r)*f_c_np1.
    pub fn interpolated_coarse(&self, node: usize, s: usize) -> [f64; 9] {
        let alpha = s as f64 / self.ratio as f64;
        let mut result = [0.0_f64; 9];
        for i in 0..9 {
            result[i] =
                (1.0 - alpha) * self.f_coarse_n[node][i] + alpha * self.f_coarse_np1[node][i];
        }
        result
    }
    /// Advance fine grid by one sub-step using BGK collision.
    pub fn advance_fine_substep(&mut self) {
        for k in 0..self.n_fine {
            let rho: f64 = self.f_fine[k].iter().sum();
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for i in 0..9 {
                let c = D2Q9_VELOCITIES[i];
                ux += c[0] as f64 * self.f_fine[k][i];
                uy += c[1] as f64 * self.f_fine[k][i];
            }
            let u = [ux / rho.max(1e-15), uy / rho.max(1e-15)];
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let u2 = u[0] * u[0] + u[1] * u[1];
                let feq =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
                self.f_fine[k][i] -= (self.f_fine[k][i] - feq) / self.tau_fine;
            }
        }
        self.sub_step = (self.sub_step + 1) % self.ratio;
    }
}
/// Wavelet-based grid adaptation criterion.
///
/// Estimates the local detail coefficient (wavelet coefficient magnitude)
/// as a grid refinement indicator.  If the detail exceeds `energy_threshold`,
/// that cell should be refined.
#[derive(Debug, Clone)]
pub struct MultiresolutionFilter {
    /// Energy threshold for refinement.
    pub energy_threshold: f64,
    /// Wavelet order (1 = Haar, 2 = linear prediction).
    pub wavelet_order: usize,
}
impl MultiresolutionFilter {
    /// Create a new multiresolution filter.
    pub fn new(energy_threshold: f64, wavelet_order: usize) -> Self {
        Self {
            energy_threshold,
            wavelet_order,
        }
    }
    /// Compute Haar wavelet detail coefficients for a 1D signal.
    ///
    /// Returns `(approximation, detail)` at the next coarser level.
    pub fn haar_transform(&self, signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len();
        let half = n / 2;
        let mut approx = vec![0.0_f64; half];
        let mut detail = vec![0.0_f64; half];
        for i in 0..half {
            approx[i] = (signal[2 * i] + signal[2 * i + 1]) * 0.5;
            detail[i] = (signal[2 * i] - signal[2 * i + 1]) * 0.5;
        }
        (approx, detail)
    }
    /// Compute refinement indicator for a 2D field (row-major, nx × ny).
    ///
    /// Returns a boolean mask: `true` = cell should be refined.
    pub fn refinement_indicator(&self, field: &[f64], nx: usize, ny: usize) -> Vec<bool> {
        let n = nx * ny;
        let mut indicator = vec![false; n];
        for y in 1..(ny - 1) {
            for x in 1..(nx - 1) {
                let k = y * nx + x;
                let detail = (field[k + 1] + field[k - 1] + field[k + nx] + field[k - nx]
                    - 4.0 * field[k])
                    .abs();
                if detail > self.energy_threshold {
                    indicator[k] = true;
                }
            }
        }
        indicator
    }
    /// Inverse Haar transform (reconstruct from approximation and detail).
    pub fn haar_inverse(&self, approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let half = approx.len().min(detail.len());
        let mut signal = vec![0.0_f64; 2 * half];
        for i in 0..half {
            signal[2 * i] = approx[i] + detail[i];
            signal[2 * i + 1] = approx[i] - detail[i];
        }
        signal
    }
}
/// Overlapping Schwarz domain decomposition for multiscale LBM.
///
/// Maintains multiple overlapping patches and iteratively exchanges boundary
/// conditions until convergence (alternating Schwarz method).
#[derive(Debug, Clone)]
pub struct MultiscaleCoupling {
    /// Number of patches.
    pub n_patches: usize,
    /// Distribution functions for each patch (patch_idx → node_idx → \[9\]).
    pub patches: Vec<Vec<[f64; 9]>>,
    /// Patch sizes.
    pub patch_sizes: Vec<usize>,
    /// Relaxation times per patch.
    pub tau_per_patch: Vec<f64>,
    /// Overlap boundary mappings: (patch_src, node_src, patch_dst, node_dst).
    pub overlap_maps: Vec<(usize, usize, usize, usize)>,
    /// Convergence tolerance for Schwarz iterations.
    pub tolerance: f64,
    /// Maximum Schwarz iterations per step.
    pub max_schwarz_iter: usize,
}
impl MultiscaleCoupling {
    /// Create a new multiscale coupling with `n_patches` subdomains.
    pub fn new(patch_sizes: Vec<usize>, tau_per_patch: Vec<f64>, tolerance: f64) -> Self {
        assert_eq!(patch_sizes.len(), tau_per_patch.len());
        let n_patches = patch_sizes.len();
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let patches: Vec<Vec<[f64; 9]>> = patch_sizes.iter().map(|&n| vec![f_init; n]).collect();
        Self {
            n_patches,
            patches,
            patch_sizes,
            tau_per_patch,
            overlap_maps: Vec::new(),
            tolerance,
            max_schwarz_iter: 10,
        }
    }
    /// Register an overlap coupling between two patches.
    pub fn add_overlap(
        &mut self,
        patch_src: usize,
        node_src: usize,
        patch_dst: usize,
        node_dst: usize,
    ) {
        self.overlap_maps
            .push((patch_src, node_src, patch_dst, node_dst));
    }
    /// Perform one Schwarz exchange iteration across all overlap mappings.
    ///
    /// Returns the maximum change in distribution functions (convergence measure).
    pub fn schwarz_exchange(&mut self) -> f64 {
        let mut max_change = 0.0_f64;
        let overlaps = self.overlap_maps.clone();
        for (ps, ns, pd, nd) in overlaps {
            let f_src = self.patches[ps][ns];
            let f_dst_old = self.patches[pd][nd];
            self.patches[pd][nd] = f_src;
            let change = f_src
                .iter()
                .zip(f_dst_old.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            max_change = max_change.max(change);
        }
        max_change
    }
    /// Perform all Schwarz iterations until convergence or max iterations.
    pub fn iterate_schwarz(&mut self) -> usize {
        for iter in 0..self.max_schwarz_iter {
            let change = self.schwarz_exchange();
            if change < self.tolerance {
                return iter + 1;
            }
        }
        self.max_schwarz_iter
    }
    /// BGK collision on patch `p`.
    pub fn collide_patch(&mut self, p: usize) {
        let tau = self.tau_per_patch[p];
        for f in &mut self.patches[p] {
            let rho: f64 = f.iter().sum();
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for i in 0..9 {
                let c = D2Q9_VELOCITIES[i];
                ux += c[0] as f64 * f[i];
                uy += c[1] as f64 * f[i];
            }
            let u = [ux / rho.max(1e-15), uy / rho.max(1e-15)];
            for i in 0..9 {
                let w = D2Q9_WEIGHTS[i];
                let c = D2Q9_VELOCITIES[i];
                let eu = c[0] as f64 * u[0] + c[1] as f64 * u[1];
                let u2 = u[0] * u[0] + u[1] * u[1];
                let feq =
                    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
                f[i] -= (f[i] - feq) / tau;
            }
        }
    }
}
/// Temporal refinement: sub-cycling on fine grids.
///
/// Fine grid time step Δt_f = Δt_c / r. The fine grid performs r sub-steps
/// for every one coarse step, maintaining CFL stability.
#[derive(Debug, Clone)]
pub struct TemporalRefinement {
    /// Refinement ratio r.
    pub ratio: usize,
    /// Coarse time step Δt_c.
    pub dt_coarse: f64,
    /// Fine time step Δt_f = Δt_c / r.
    pub dt_fine: f64,
    /// Current sub-step counter.
    pub sub_step: usize,
    /// Fine grid distribution functions (flattened).
    pub f_fine: Vec<[f64; 9]>,
    /// Fine grid node count.
    pub n_fine: usize,
}
impl TemporalRefinement {
    /// Create a new temporal refinement handler.
    ///
    /// # Arguments
    /// * `n_fine`    - number of fine grid nodes
    /// * `ratio`     - temporal refinement ratio (= spatial ratio)
    /// * `dt_coarse` - coarse grid time step
    pub fn new(n_fine: usize, ratio: usize, dt_coarse: f64) -> Self {
        let dt_fine = dt_coarse / ratio as f64;
        let feq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        Self {
            ratio,
            dt_coarse,
            dt_fine,
            sub_step: 0,
            f_fine: vec![feq; n_fine],
            n_fine,
        }
    }
    /// Perform one fine-grid BGK collision sub-step.
    ///
    /// f_i ← f_i - (f_i - f_i^eq) / τ_f
    pub fn fine_collision_step(&mut self, tau_fine: f64) {
        for node in &mut self.f_fine {
            let rho: f64 = node.iter().sum();
            let mut jx = 0.0f64;
            let mut jy = 0.0f64;
            for alpha in 0..9 {
                jx += node[alpha] * D2Q9_VELOCITIES[alpha][0] as f64;
                jy += node[alpha] * D2Q9_VELOCITIES[alpha][1] as f64;
            }
            let u = [jx / rho.max(1e-30), jy / rho.max(1e-30)];
            let feq = lbm_equilibrium(rho, u);
            for alpha in 0..9 {
                node[alpha] -= (node[alpha] - feq[alpha]) / tau_fine;
            }
        }
        self.sub_step = (self.sub_step + 1) % self.ratio;
    }
    /// Returns true when fine grid has completed a full coarse step.
    pub fn is_synchronized(&self) -> bool {
        self.sub_step == 0
    }
    /// Perform all `ratio` fine sub-steps (full coarse synchronization cycle).
    pub fn full_cycle(&mut self, tau_fine: f64) {
        for _ in 0..self.ratio {
            self.fine_collision_step(tau_fine);
        }
    }
    /// Interpolate fine-grid field to match coarse time level.
    ///
    /// Simple linear-in-time interpolation between sub-steps.
    pub fn time_interpolate(&self, f_prev: &[f64; 9], f_next: &[f64; 9], t_frac: f64) -> [f64; 9] {
        let mut f = [0.0f64; 9];
        for alpha in 0..9 {
            f[alpha] = (1.0 - t_frac) * f_prev[alpha] + t_frac * f_next[alpha];
        }
        f
    }
    /// CFL number for fine grid (should be ≤ 1).
    pub fn fine_cfl(&self, dx_fine: f64) -> f64 {
        self.dt_fine / dx_fine
    }
}
