//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::boundary::{Boundary, apply_boundaries_2d};
use crate::collision::{bgk_collide_2d, bgk_collide_3d};
use crate::grid::{LbmGrid2D, LbmGrid3D, equilibrium_2d};
use crate::lattice::{CS2, LatticeType};
use crate::streaming::{stream_2d, stream_3d};
use crate::turbulence::smagorinsky_omega;

/// Result of a single ensemble member run.
#[derive(Debug, Clone)]
pub struct EnsembleMemberResult {
    /// Index of this ensemble member (0-based).
    pub member_idx: usize,
    /// Final step count.
    pub steps: usize,
    /// Whether the run converged.
    pub converged: bool,
    /// Final mean velocity.
    pub mean_velocity: f64,
    /// Final mean density.
    pub mean_density: f64,
    /// Final Mach number.
    pub mach: f64,
}
/// Guo body-force scheme parameters.
#[derive(Debug, Clone)]
pub struct GuoForceParams {
    /// Body force vector `[fx, fy, fz]`.
    pub force: [f64; 3],
    /// Relaxation frequency omega (used to compute the correction factor).
    pub omega: f64,
}
impl GuoForceParams {
    /// Create new Guo force parameters.
    pub fn new(force: [f64; 3], omega: f64) -> Self {
        Self { force, omega }
    }
    /// Correction prefactor `(1 - omega/2)`.
    pub fn prefactor(&self) -> f64 {
        1.0 - 0.5 * self.omega
    }
    /// Compute the force term for direction `i` of D2Q9 at macroscopic velocity `(ux, uy)`.
    ///
    /// `Fi = wi * (1 - omega/2) * [ (ei/cs2 + (ei·u)ei/cs4) · F ]`
    pub fn d2q9_force_term(&self, i: usize, ux: f64, uy: f64) -> f64 {
        use crate::lattice::CS2;
        let weights = [
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
        let ex = [0i32, 1, 0, -1, 0, 1, -1, -1, 1];
        let ey = [0i32, 0, 1, 0, -1, 1, 1, -1, -1];
        let cx = ex[i] as f64;
        let cy = ey[i] as f64;
        let eu = cx * ux + cy * uy;
        let pref = self.prefactor();
        let fx = self.force[0];
        let fy = self.force[1];
        weights[i]
            * pref
            * ((cx / CS2 + eu * cx / (CS2 * CS2)) * fx + (cy / CS2 + eu * cy / (CS2 * CS2)) * fy)
    }
}
/// Snapshot of a simulation state for restart purposes.
#[derive(Debug, Clone)]
pub struct SimulationSnapshot {
    /// Step count at snapshot time.
    pub step_count: usize,
    /// Density field.
    pub density: Vec<f64>,
    /// x-velocity field.
    pub ux: Vec<f64>,
    /// y-velocity field.
    pub uy: Vec<f64>,
}
impl SimulationSnapshot {
    /// Create a snapshot from the current simulation state.
    pub fn from_simulation(sim: &LbmSimulation) -> Self {
        let vel = sim.get_velocity_field();
        Self {
            step_count: sim.step_count,
            density: sim.get_density_field(),
            ux: vel.iter().map(|v| v[0]).collect(),
            uy: vel.iter().map(|v| v[1]).collect(),
        }
    }
    /// Number of cells in the snapshot.
    pub fn n_cells(&self) -> usize {
        self.density.len()
    }
}
impl SimulationSnapshot {
    /// Compute the relative change in ux between two snapshots.
    pub fn velocity_change(&self, other: &Self) -> f64 {
        assert_eq!(self.ux.len(), other.ux.len(), "Snapshot size mismatch");
        let mut diff_sq = 0.0_f64;
        let mut ref_sq = 0.0_f64;
        for (a, b) in self.ux.iter().zip(other.ux.iter()) {
            diff_sq += (a - b).powi(2);
            ref_sq += b * b;
        }
        if ref_sq < 1e-30 {
            return diff_sq.sqrt();
        }
        (diff_sq / ref_sq).sqrt()
    }
    /// Check whether the density field is physically valid (all > 0).
    pub fn is_density_valid(&self) -> bool {
        self.density.iter().all(|&r| r > 0.0)
    }
    /// Check whether the velocity field is bounded (|u| < threshold).
    pub fn is_velocity_bounded(&self, threshold: f64) -> bool {
        self.ux
            .iter()
            .zip(self.uy.iter())
            .all(|(ux, uy)| (ux * ux + uy * uy).sqrt() < threshold)
    }
    /// Mean density of the snapshot.
    pub fn mean_density(&self) -> f64 {
        if self.density.is_empty() {
            return 0.0;
        }
        self.density.iter().sum::<f64>() / self.density.len() as f64
    }
}
/// MRT relaxation parameter set for D3Q19 lattice.
#[derive(Debug, Clone)]
pub struct MrtParamsD3Q19 {
    /// Relaxation rates for all 19 moments.
    pub s: [f64; 19],
}
impl MrtParamsD3Q19 {
    /// Construct with uniform viscous rate omega and bulk rate s_bulk.
    pub fn uniform(omega: f64, s_bulk: f64) -> Self {
        let mut s = [s_bulk; 19];
        s[0] = 0.0;
        s[3] = 0.0;
        s[5] = 0.0;
        s[7] = 0.0;
        s[9] = omega;
        s[11] = omega;
        s[13] = omega;
        s[14] = omega;
        s[15] = omega;
        Self { s }
    }
}
/// Fully-featured LBM simulation runner.
///
/// Supports D2Q9, D3Q19, and D3Q27 lattices with optional body forcing
/// and Smagorinsky turbulence.  Use [`LbmSimulation::new`] with an
/// [`LbmConfig`] to create an instance.
///
/// # Example
/// ```no_run
/// use oxiphysics_lbm::simulation::{LbmConfig, LbmSimulation};
///
/// let mut cfg = LbmConfig::d2q9(20, 10, 1.0 / 6.0);
/// cfg.body_force = Some([1e-5, 0.0, 0.0]);
/// let mut sim = LbmSimulation::new(cfg);
/// sim.step_n(100);
/// let vel = sim.get_velocity_field();
/// assert_eq!(vel.len(), 20 * 10);
/// ```
#[derive(Debug, Clone)]
pub struct LbmSimulation {
    /// Simulation configuration.
    pub config: LbmConfig,
    /// Internal grid storage.
    pub(super) grid: SimGrid,
    /// 2D boundary conditions (used only for 2D grids).
    pub(super) boundaries_2d: Vec<Boundary>,
    /// Current step counter.
    pub step_count: usize,
}
impl LbmSimulation {
    /// Create a new simulation from a configuration.
    pub fn new(config: LbmConfig) -> Self {
        let grid = match config.lattice_type {
            LatticeType::D2Q9 => {
                SimGrid::Grid2D(LbmGrid2D::new(config.nx, config.ny, LatticeType::D2Q9))
            }
            LatticeType::D3Q19 => SimGrid::Grid3D(LbmGrid3D::new(
                config.nx,
                config.ny,
                config.nz,
                LatticeType::D3Q19,
            )),
            LatticeType::D3Q27 => SimGrid::Grid3D(LbmGrid3D::new(
                config.nx,
                config.ny,
                config.nz,
                LatticeType::D3Q27,
            )),
        };
        Self {
            config,
            grid,
            boundaries_2d: Vec::new(),
            step_count: 0,
        }
    }
    /// Set 2D boundary conditions (only applicable to D2Q9 simulations).
    pub fn set_boundaries(&mut self, boundaries: Vec<Boundary>) {
        self.boundaries_2d = boundaries;
    }
    /// Perform a single complete LBM step.
    ///
    /// The step sequence is:
    /// 1. BGK collision (with optional Smagorinsky turbulence model for 2D)
    /// 2. Optional body-force application (simple velocity shift for 2D)
    /// 3. Streaming
    /// 4. Boundary conditions
    /// 5. Macroscopic recomputation
    pub fn step(&mut self, _dt: f64) {
        let omega = self.config.omega();
        let force = self.config.body_force;
        let smag_cs = self.config.smagorinsky_cs;
        match &mut self.grid {
            SimGrid::Grid2D(g) => {
                if let Some(cs) = smag_cs {
                    let nx = g.nx;
                    let ny = g.ny;
                    let q = g.lattice.q();
                    g.compute_macroscopic();
                    for y in 0..ny {
                        for x in 0..nx {
                            let local_omega = smagorinsky_omega(g, cs, omega, x, y);
                            let k = g.idx(x, y);
                            let rho_k = g.rho[k];
                            let ux_k = g.ux[k];
                            let uy_k = g.uy[k];
                            for i in 0..q {
                                let w = g.lattice.weight(i);
                                let c = g.lattice.velocity_2d(i);
                                let feq =
                                    equilibrium_2d(w, rho_k, ux_k, uy_k, c[0] as f64, c[1] as f64);
                                g.f[i][k] -= local_omega * (g.f[i][k] - feq);
                            }
                        }
                    }
                } else {
                    bgk_collide_2d(g, omega);
                }
                if let Some(f) = force {
                    let nx = g.nx;
                    let ny = g.ny;
                    let q = g.lattice.q();
                    g.compute_macroscopic();
                    for y in 0..ny {
                        for x in 0..nx {
                            let k = g.idx(x, y);
                            let rho = g.rho[k];
                            g.ux[k] += f[0] * 0.5 / rho;
                            g.uy[k] += f[1] * 0.5 / rho;
                            let ux_k = g.ux[k];
                            let uy_k = g.uy[k];
                            for i in 0..q {
                                let w = g.lattice.weight(i);
                                let c = g.lattice.velocity_2d(i);
                                let cx = c[0] as f64;
                                let cy = c[1] as f64;
                                let tau = 1.0 / omega;
                                let eu = cx * ux_k + cy * uy_k;
                                let fi_term = w
                                    * (1.0 - 0.5 * omega)
                                    * ((cx / CS2 + eu * cx / (CS2 * CS2)) * f[0]
                                        + (cy / CS2 + eu * cy / (CS2 * CS2)) * f[1]);
                                let _ = tau;
                                g.f[i][k] += fi_term;
                            }
                        }
                    }
                }
                stream_2d(g);
                let boundaries = self.boundaries_2d.clone();
                apply_boundaries_2d(g, &boundaries);
                g.compute_macroscopic();
            }
            SimGrid::Grid3D(g) => {
                bgk_collide_3d(g, omega);
                stream_3d(g);
                g.compute_macroscopic();
            }
        }
        self.step_count += 1;
    }
    /// Perform `n` complete LBM steps.
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step(1.0);
        }
    }
    /// Return the velocity field as a `Vec<[f64; 3]>` (uz = 0 for 2D).
    ///
    /// The flat index matches `z * ny * nx + y * nx + x`.
    pub fn get_velocity_field(&self) -> Vec<[f64; 3]> {
        match &self.grid {
            SimGrid::Grid2D(g) => {
                let n = g.nx * g.ny;
                (0..n).map(|k| [g.ux[k], g.uy[k], 0.0]).collect()
            }
            SimGrid::Grid3D(g) => {
                let n = g.nx * g.ny * g.nz;
                (0..n).map(|k| [g.ux[k], g.uy[k], g.uz[k]]).collect()
            }
        }
    }
    /// Return the pressure field `p = cs² * rho` for every cell.
    pub fn get_pressure_field(&self) -> Vec<f64> {
        match &self.grid {
            SimGrid::Grid2D(g) => g.rho.iter().map(|&r| r * CS2).collect(),
            SimGrid::Grid3D(g) => g.rho.iter().map(|&r| r * CS2).collect(),
        }
    }
    /// Return the density field for every cell.
    pub fn get_density_field(&self) -> Vec<f64> {
        match &self.grid {
            SimGrid::Grid2D(g) => g.rho.clone(),
            SimGrid::Grid3D(g) => g.rho.clone(),
        }
    }
    /// Compute the relative L2 error of the x-velocity profile against the
    /// analytical Poiseuille parabola.
    ///
    /// Only meaningful for 2D channel flow with walls at y = 0 and y = ny-1.
    /// Requires that the flow has been driven long enough to reach steady state.
    ///
    /// Returns the relative L2 error `||u_num - u_ana||₂ / ||u_ana||₂`.
    pub fn poiseuille_flow_error(&self) -> f64 {
        match &self.grid {
            SimGrid::Grid2D(g) => {
                let nx = g.nx;
                let ny = g.ny;
                let h = (ny - 2) as f64;
                let mid_x = nx / 2;
                let mut u_max = 0.0_f64;
                for y in 1..(ny - 1) {
                    let k = g.idx(mid_x, y);
                    u_max = u_max.max(g.ux[k]);
                }
                if u_max < 1e-15 {
                    return f64::INFINITY;
                }
                let nu = self.config.viscosity;
                let dp_dx = u_max * 8.0 * nu / (h * h);
                let mut sq_num = 0.0_f64;
                let mut sq_ana = 0.0_f64;
                for y in 1..(ny - 1) {
                    let yf = y as f64;
                    let u_ana = dp_dx / (2.0 * nu) * (yf - 0.5) * ((ny as f64 - 1.5) - yf);
                    for x in 0..nx {
                        let k = g.idx(x, y);
                        let diff = g.ux[k] - u_ana;
                        sq_num += diff * diff;
                        sq_ana += u_ana * u_ana;
                    }
                }
                if sq_ana < 1e-30 {
                    return f64::INFINITY;
                }
                (sq_num / sq_ana).sqrt()
            }
            SimGrid::Grid3D(_) => f64::NAN,
        }
    }
    /// Borrow the inner 2D grid (panics for 3D simulations).
    pub fn grid_2d(&self) -> &LbmGrid2D {
        match &self.grid {
            SimGrid::Grid2D(g) => g,
            SimGrid::Grid3D(_) => panic!("grid_2d() called on a 3D simulation"),
        }
    }
    /// Borrow the inner 3D grid (panics for 2D simulations).
    pub fn grid_3d(&self) -> &LbmGrid3D {
        match &self.grid {
            SimGrid::Grid2D(_) => panic!("grid_3d() called on a 2D simulation"),
            SimGrid::Grid3D(g) => g,
        }
    }
}
/// Statistics collected during a simulation run.
#[derive(Debug, Clone)]
pub struct SimulationStatistics {
    /// Mean velocity magnitude over time.
    pub mean_velocity_history: Vec<f64>,
    /// Mean density over time.
    pub mean_density_history: Vec<f64>,
    /// Maximum velocity magnitude over time.
    pub max_velocity_history: Vec<f64>,
    /// Mach number history.
    pub mach_history: Vec<f64>,
}
impl SimulationStatistics {
    /// Create new empty statistics.
    pub fn new() -> Self {
        Self {
            mean_velocity_history: Vec::new(),
            mean_density_history: Vec::new(),
            max_velocity_history: Vec::new(),
            mach_history: Vec::new(),
        }
    }
    /// Record statistics from the current simulation state.
    pub fn record(&mut self, sim: &LbmSimulation) {
        let vel = sim.get_velocity_field();
        let rho = sim.get_density_field();
        let mut sum_vel = 0.0;
        let mut max_vel = 0.0_f64;
        for v in &vel {
            let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            sum_vel += mag;
            max_vel = max_vel.max(mag);
        }
        let n = vel.len() as f64;
        self.mean_velocity_history.push(sum_vel / n);
        self.max_velocity_history.push(max_vel);
        let mean_rho: f64 = rho.iter().sum::<f64>() / rho.len() as f64;
        self.mean_density_history.push(mean_rho);
        let cs = (1.0 / 3.0_f64).sqrt();
        self.mach_history.push(max_vel / cs);
    }
    /// Return the latest mean velocity.
    pub fn latest_mean_velocity(&self) -> f64 {
        self.mean_velocity_history.last().copied().unwrap_or(0.0)
    }
    /// Return the latest Mach number.
    pub fn latest_mach(&self) -> f64 {
        self.mach_history.last().copied().unwrap_or(0.0)
    }
    /// Check if the simulation is stable (Ma < 0.3).
    pub fn is_stable(&self) -> bool {
        self.latest_mach() < 0.3
    }
}
impl SimulationStatistics {
    /// Return the minimum Mach number recorded.
    pub fn min_mach(&self) -> f64 {
        self.mach_history.iter().cloned().fold(f64::MAX, f64::min)
    }
    /// Return the maximum Mach number recorded.
    pub fn max_mach(&self) -> f64 {
        self.mach_history.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Return the mean of a slice.
    fn mean_of(v: &[f64]) -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.iter().sum::<f64>() / v.len() as f64
    }
    /// Mean velocity across the entire history.
    pub fn mean_velocity_overall(&self) -> f64 {
        Self::mean_of(&self.mean_velocity_history)
    }
    /// Mean density across the entire history.
    pub fn mean_density_overall(&self) -> f64 {
        Self::mean_of(&self.mean_density_history)
    }
    /// Number of statistics records collected.
    pub fn n_records(&self) -> usize {
        self.mean_velocity_history.len()
    }
}
/// In-memory checkpoint store for simulation state.
pub struct CheckpointStore {
    /// List of saved snapshots.
    pub snapshots: Vec<SimulationSnapshot>,
    /// Maximum number of checkpoints to retain.
    pub capacity: usize,
}
impl CheckpointStore {
    /// Create a new checkpoint store with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            snapshots: Vec::new(),
            capacity,
        }
    }
    /// Save a checkpoint of the current simulation state.
    pub fn save(&mut self, sim: &LbmSimulation) {
        if self.snapshots.len() >= self.capacity {
            self.snapshots.remove(0);
        }
        self.snapshots
            .push(SimulationSnapshot::from_simulation(sim));
    }
    /// Return the most recent checkpoint, if any.
    pub fn latest(&self) -> Option<&SimulationSnapshot> {
        self.snapshots.last()
    }
    /// Number of stored checkpoints.
    pub fn len(&self) -> usize {
        self.snapshots.len()
    }
    /// Whether the store is empty.
    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }
    /// Clear all stored checkpoints.
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }
}
/// A 2D LBM simulation combining grid, boundaries, and solver parameters.
#[derive(Debug, Clone)]
pub struct LbmSimulation2D {
    /// The computational grid.
    pub grid: LbmGrid2D,
    /// Boundary conditions.
    pub boundaries: Vec<Boundary>,
    /// Relaxation frequency omega = 1/tau.
    pub omega: f64,
    /// Current time step counter.
    pub step_count: usize,
}
impl LbmSimulation2D {
    /// Create a new simulation with the given grid size and relaxation frequency.
    pub fn new(nx: usize, ny: usize, omega: f64) -> Self {
        Self {
            grid: LbmGrid2D::new(nx, ny, LatticeType::D2Q9),
            boundaries: Vec::new(),
            omega,
            step_count: 0,
        }
    }
    /// Create a simulation from an existing grid.
    pub fn from_grid(grid: LbmGrid2D, omega: f64) -> Self {
        Self {
            grid,
            boundaries: Vec::new(),
            omega,
            step_count: 0,
        }
    }
    /// Set boundary conditions.
    pub fn set_boundaries(&mut self, boundaries: Vec<Boundary>) {
        self.boundaries = boundaries;
    }
    /// Perform a single time step: collide → stream → apply boundaries.
    pub fn step(&mut self) {
        bgk_collide_2d(&mut self.grid, self.omega);
        stream_2d(&mut self.grid);
        apply_boundaries_2d(&mut self.grid, &self.boundaries);
        self.grid.compute_macroscopic();
        self.step_count += 1;
    }
    /// Run the simulation for `n_steps` time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }
}
/// MRT relaxation parameter set for D2Q9 lattice.
///
/// The 9 relaxation rates correspond to the 9 moment equations:
/// rho, energy, energy-square, jx, qx, jy, qy, pxx, pxy.
#[derive(Debug, Clone)]
pub struct MrtParamsD2Q9 {
    /// Relaxation rates s\[0..9\]. Conserved modes have s=0.
    pub s: [f64; 9],
}
impl MrtParamsD2Q9 {
    /// Standard MRT rates for D2Q9 (Lallemand & Luo 2000).
    ///
    /// Conserved moments (density, momentum) have s=0.
    /// Viscous moments use the BGK-equivalent rate omega = 1/tau.
    pub fn standard(omega: f64) -> Self {
        Self {
            s: [0.0, 1.63, 1.14, 0.0, 1.92, 0.0, 1.92, omega, omega],
        }
    }
    /// Maximally stable MRT: all non-conserved rates = 1.8.
    pub fn max_stable() -> Self {
        Self {
            s: [0.0, 1.8, 1.8, 0.0, 1.8, 0.0, 1.8, 1.8, 1.8],
        }
    }
    /// Return the effective viscosity rate (s\[7\] which controls shear stress).
    pub fn shear_omega(&self) -> f64 {
        self.s[7]
    }
}
/// Ensemble runner that executes multiple independent simulations and
/// averages the results.
///
/// Each member shares the same `LbmConfig` but may have a different initial
/// perturbation applied via `initial_ux_perturbation`.
pub struct EnsembleRunner {
    /// Shared simulation configuration.
    pub config: LbmConfig,
    /// Number of ensemble members.
    pub n_members: usize,
    /// Steps per member run.
    pub steps_per_member: usize,
    /// Small initial perturbation amplitude applied to ux (lattice units).
    pub perturbation_amplitude: f64,
    /// Collected results from all members.
    pub results: Vec<EnsembleMemberResult>,
}
impl EnsembleRunner {
    /// Create a new ensemble runner.
    pub fn new(config: LbmConfig, n_members: usize, steps_per_member: usize) -> Self {
        Self {
            config,
            n_members,
            steps_per_member,
            perturbation_amplitude: 1e-6,
            results: Vec::new(),
        }
    }
    /// Run all ensemble members and store individual results.
    ///
    /// Each member is initialised fresh from `config`.  A deterministic
    /// perturbation proportional to the member index is added to ux so the
    /// members are not identical.
    pub fn run_all(&mut self) {
        self.results.clear();
        for idx in 0..self.n_members {
            let result = self.run_member(idx);
            self.results.push(result);
        }
    }
    /// Run a single ensemble member and return its result.
    fn run_member(&self, member_idx: usize) -> EnsembleMemberResult {
        let mut sim = LbmSimulation::new(self.config.clone());
        let perturb = self.perturbation_amplitude * (member_idx as f64 + 1.0);
        if let SimGrid::Grid2D(ref mut g) = sim.grid {
            for k in 0..g.nx * g.ny {
                g.ux[k] += perturb * ((k % 7) as f64 - 3.0);
            }
        }
        let mut stats = SimulationStatistics::new();
        let n_cells = sim.get_velocity_field().len();
        let mut monitor = ConvergenceMonitor::new(n_cells, 1e-6, 100);
        let mut converged = false;
        for step in 0..self.steps_per_member {
            sim.step(1.0);
            if step % 100 == 99 {
                let vel = sim.get_velocity_field();
                let ux: Vec<f64> = vel.iter().map(|v| v[0]).collect();
                let uy: Vec<f64> = vel.iter().map(|v| v[1]).collect();
                converged = monitor.check(&ux, &uy);
                if converged {
                    break;
                }
            }
        }
        stats.record(&sim);
        EnsembleMemberResult {
            member_idx,
            steps: sim.step_count,
            converged,
            mean_velocity: stats.latest_mean_velocity(),
            mean_density: stats.mean_density_overall(),
            mach: stats.latest_mach(),
        }
    }
    /// Compute the ensemble-averaged mean velocity across all members.
    pub fn ensemble_mean_velocity(&self) -> f64 {
        if self.results.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.results.iter().map(|r| r.mean_velocity).sum();
        sum / self.results.len() as f64
    }
    /// Compute the ensemble-averaged Mach number.
    pub fn ensemble_mean_mach(&self) -> f64 {
        if self.results.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.results.iter().map(|r| r.mach).sum();
        sum / self.results.len() as f64
    }
    /// Count how many members converged.
    pub fn n_converged(&self) -> usize {
        self.results.iter().filter(|r| r.converged).count()
    }
    /// Compute variance of mean velocity across members.
    pub fn velocity_variance(&self) -> f64 {
        if self.results.len() < 2 {
            return 0.0;
        }
        let mean = self.ensemble_mean_velocity();
        let var: f64 = self
            .results
            .iter()
            .map(|r| (r.mean_velocity - mean).powi(2))
            .sum::<f64>()
            / (self.results.len() - 1) as f64;
        var
    }
}
/// Adaptive step controller for LBM simulations.
///
/// Adjusts the effective relaxation parameter based on stability criteria.
pub struct StepController {
    /// Minimum allowed viscosity.
    pub nu_min: f64,
    /// Maximum allowed viscosity.
    pub nu_max: f64,
    /// Current effective viscosity.
    pub nu_current: f64,
    /// Safety factor for stability (0 < factor <= 1).
    pub safety_factor: f64,
}
impl StepController {
    /// Create a new step controller.
    pub fn new(nu_min: f64, nu_max: f64) -> Self {
        Self {
            nu_min,
            nu_max,
            nu_current: nu_min,
            safety_factor: 0.9,
        }
    }
    /// Adjust viscosity based on maximum velocity in the domain.
    ///
    /// Ensures Ma = u_max / cs < 0.3 for stability.
    pub fn adjust_viscosity(&mut self, u_max: f64) {
        let cs = (1.0 / 3.0_f64).sqrt();
        let ma = u_max / cs;
        if ma > 0.3 {
            self.nu_current = (self.nu_current * 1.1).min(self.nu_max);
        } else if ma < 0.1 {
            self.nu_current = (self.nu_current * 0.95).max(self.nu_min);
        }
    }
    /// Return the current relaxation frequency.
    pub fn omega(&self) -> f64 {
        let tau = self.nu_current / CS2 + 0.5;
        1.0 / tau
    }
}
/// Monitors convergence of an LBM simulation by tracking velocity residuals.
pub struct ConvergenceMonitor {
    /// Previous velocity field snapshot (flat).
    pub(super) prev_ux: Vec<f64>,
    /// Previous velocity field snapshot (flat).
    pub(super) prev_uy: Vec<f64>,
    /// History of L2 residuals.
    pub residual_history: Vec<f64>,
    /// Convergence threshold.
    pub tolerance: f64,
    /// Check interval (every N steps).
    pub check_interval: usize,
}
impl ConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(n_cells: usize, tolerance: f64, check_interval: usize) -> Self {
        Self {
            prev_ux: vec![0.0; n_cells],
            prev_uy: vec![0.0; n_cells],
            residual_history: Vec::new(),
            tolerance,
            check_interval,
        }
    }
    /// Update the monitor with current velocities. Returns true if converged.
    pub fn check(&mut self, ux: &[f64], uy: &[f64]) -> bool {
        let n = ux.len();
        let mut diff_sq = 0.0;
        let mut ref_sq = 0.0;
        for i in 0..n {
            let dx = ux[i] - self.prev_ux[i];
            let dy = uy[i] - self.prev_uy[i];
            diff_sq += dx * dx + dy * dy;
            ref_sq += ux[i] * ux[i] + uy[i] * uy[i];
        }
        let residual = if ref_sq > 1e-30 {
            (diff_sq / ref_sq).sqrt()
        } else {
            diff_sq.sqrt()
        };
        self.residual_history.push(residual);
        self.prev_ux.copy_from_slice(ux);
        self.prev_uy.copy_from_slice(uy);
        residual < self.tolerance
    }
    /// Return the latest residual.
    pub fn latest_residual(&self) -> f64 {
        self.residual_history.last().copied().unwrap_or(f64::MAX)
    }
    /// Check if residuals are decreasing (last 3 values).
    pub fn is_decreasing(&self) -> bool {
        let n = self.residual_history.len();
        if n < 3 {
            return true;
        }
        self.residual_history[n - 1] < self.residual_history[n - 2]
            && self.residual_history[n - 2] < self.residual_history[n - 3]
    }
}
impl ConvergenceMonitor {
    /// Return the minimum residual ever recorded.
    pub fn min_residual(&self) -> f64 {
        self.residual_history
            .iter()
            .cloned()
            .fold(f64::MAX, f64::min)
    }
    /// Return the maximum residual ever recorded.
    pub fn max_residual(&self) -> f64 {
        self.residual_history
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }
    /// Number of residuals recorded so far.
    pub fn n_records(&self) -> usize {
        self.residual_history.len()
    }
    /// Check if the residual trend is stagnating (last 5 changes < 1e-4 relative).
    pub fn is_stagnating(&self) -> bool {
        let n = self.residual_history.len();
        if n < 6 {
            return false;
        }
        let last = self.residual_history[n - 1];
        let older = self.residual_history[n - 6];
        if older.abs() < 1e-30 {
            return true;
        }
        ((last - older) / older).abs() < 1e-4
    }
}
/// Configuration for a general LBM simulation.
///
/// Specifies grid dimensions, fluid properties, and optional physics modules.
#[derive(Debug, Clone)]
pub struct LbmConfig {
    /// Grid size in x-direction.
    pub nx: usize,
    /// Grid size in y-direction.
    pub ny: usize,
    /// Grid size in z-direction (set to 1 for 2D simulations).
    pub nz: usize,
    /// Kinematic viscosity in lattice units.
    pub viscosity: f64,
    /// Maximum inlet/reference velocity in lattice units.
    pub max_velocity: f64,
    /// Lattice type: D2Q9, D3Q19, or D3Q27.
    pub lattice_type: LatticeType,
    /// Optional body-force vector `[Fx, Fy, Fz]` applied every step.
    pub body_force: Option<[f64; 3]>,
    /// Optional Smagorinsky constant for turbulence modelling.
    /// When `None`, plain BGK collision is used.
    pub smagorinsky_cs: Option<f64>,
}
impl LbmConfig {
    /// Create a D2Q9 configuration with sensible defaults.
    pub fn d2q9(nx: usize, ny: usize, viscosity: f64) -> Self {
        Self {
            nx,
            ny,
            nz: 1,
            viscosity,
            max_velocity: 0.1,
            lattice_type: LatticeType::D2Q9,
            body_force: None,
            smagorinsky_cs: None,
        }
    }
    /// Create a D3Q19 configuration.
    pub fn d3q19(nx: usize, ny: usize, nz: usize, viscosity: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            viscosity,
            max_velocity: 0.1,
            lattice_type: LatticeType::D3Q19,
            body_force: None,
            smagorinsky_cs: None,
        }
    }
    /// Create a D3Q27 configuration.
    pub fn d3q27(nx: usize, ny: usize, nz: usize, viscosity: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            viscosity,
            max_velocity: 0.1,
            lattice_type: LatticeType::D3Q27,
            body_force: None,
            smagorinsky_cs: None,
        }
    }
    /// Compute the BGK relaxation frequency from viscosity.
    ///
    /// `ω = 1 / τ`,  `τ = ν / cs² + 0.5 = 3ν + 0.5`
    pub fn omega(&self) -> f64 {
        let tau = self.viscosity / CS2 + 0.5;
        1.0 / tau
    }
}
impl LbmConfig {
    /// Kinematic viscosity from omega: ν = cs²(1/ω - 0.5)
    pub fn viscosity_from_omega(omega: f64) -> f64 {
        CS2 * (1.0 / omega - 0.5)
    }
    /// Reynolds number: Re = U * L / ν
    pub fn reynolds_number(&self, u_ref: f64, l_ref: f64) -> f64 {
        u_ref * l_ref / self.viscosity.max(1e-30)
    }
    /// Mach number at the reference velocity.
    pub fn mach_number(&self) -> f64 {
        let cs = (1.0_f64 / 3.0).sqrt();
        self.max_velocity / cs
    }
    /// Check if the Mach number is within the incompressible limit (Ma < 0.3).
    pub fn is_incompressible(&self) -> bool {
        self.mach_number() < 0.3
    }
}
/// Adaptive dt controller that adjusts the physical time-step based on the
/// current CFL (Courant–Friedrichs–Lewy) condition and a target Mach number.
///
/// In LBM the lattice spacing `Δx` and lattice time-step `Δt_lat` are both
/// unity by convention.  The *physical* time-step is related to the physical
/// velocity by:
///
/// `Δt_phys = Ma_target * cs_lat / u_ref_phys * Δt_lat`
///
/// This controller recomputes the physical dt every `adjust_interval` steps
/// so that the maximum Mach number in the simulation stays close to
/// `ma_target`.
#[derive(Debug, Clone)]
pub struct AdaptiveDtController {
    /// Target Mach number (dimensionless).
    pub ma_target: f64,
    /// Minimum allowed Mach number (clamp from below).
    pub ma_min: f64,
    /// Maximum allowed Mach number (clamp from above, for stability).
    pub ma_max: f64,
    /// How many simulation steps between dt adjustments.
    pub adjust_interval: usize,
    /// Safety factor (0 < safety <= 1) to stay below stability limit.
    pub safety: f64,
    /// Current physical time-step (arbitrary units, starts at 1).
    pub dt_current: f64,
    /// History of dt values after each adjustment.
    pub dt_history: Vec<f64>,
}
impl AdaptiveDtController {
    /// Create a new adaptive dt controller with sensible defaults.
    ///
    /// # Arguments
    /// * `ma_target`       - target Ma, e.g. 0.1 for incompressible LBM
    /// * `adjust_interval` - steps between adjustments
    pub fn new(ma_target: f64, adjust_interval: usize) -> Self {
        Self {
            ma_target,
            ma_min: 0.01,
            ma_max: 0.3,
            adjust_interval,
            safety: 0.95,
            dt_current: 1.0,
            dt_history: Vec::new(),
        }
    }
    /// Update dt given the current maximum velocity in the domain.
    ///
    /// Returns the new dt.
    pub fn update(&mut self, u_max: f64) -> f64 {
        let cs = (1.0_f64 / 3.0).sqrt();
        if u_max < 1e-30 {
            self.dt_history.push(self.dt_current);
            return self.dt_current;
        }
        let ma_current = u_max / cs;
        let scale = (self.ma_target / ma_current) * self.safety;
        let dt_new = (self.dt_current * scale)
            .max(self.dt_current * (self.ma_min / ma_current.max(1e-30)))
            .min(self.dt_current * (self.ma_max / ma_current.max(1e-30)));
        self.dt_current = dt_new;
        self.dt_history.push(dt_new);
        dt_new
    }
    /// Return the latest recorded dt.
    pub fn latest_dt(&self) -> f64 {
        self.dt_history.last().copied().unwrap_or(self.dt_current)
    }
    /// Number of adjustments made so far.
    pub fn n_adjustments(&self) -> usize {
        self.dt_history.len()
    }
    /// Check whether the current Ma is within acceptable limits.
    pub fn is_stable(&self, u_max: f64) -> bool {
        let cs = (1.0_f64 / 3.0).sqrt();
        let ma = u_max / cs;
        ma < self.ma_max
    }
}
/// Full-featured simulation loop with checkpoint saving, convergence
/// monitoring, adaptive dt, and user-supplied output callbacks.
pub struct FullStepLoop {
    /// Maximum steps to run.
    pub max_steps: usize,
    /// Convergence tolerance for residual.
    pub tolerance: f64,
    /// Check convergence every `check_interval` steps.
    pub check_interval: usize,
    /// Save checkpoint every `checkpoint_interval` steps (0 = disabled).
    pub checkpoint_interval: usize,
    /// Output callback every `output_interval` steps (0 = disabled).
    pub output_interval: usize,
    /// Internal checkpoint store.
    pub checkpoint_store: CheckpointStore,
    /// Adaptive dt controller (optional).
    pub adaptive_dt: Option<AdaptiveDtController>,
    /// Residuals recorded during the run.
    pub residuals: Vec<f64>,
    /// Statistics collected during the run.
    pub stats: SimulationStatistics,
}
impl FullStepLoop {
    /// Create a new `FullStepLoop` with the given settings.
    pub fn new(
        max_steps: usize,
        tolerance: f64,
        check_interval: usize,
        checkpoint_interval: usize,
    ) -> Self {
        Self {
            max_steps,
            tolerance,
            check_interval,
            checkpoint_interval,
            output_interval: 0,
            checkpoint_store: CheckpointStore::new(10),
            adaptive_dt: None,
            residuals: Vec::new(),
            stats: SimulationStatistics::new(),
        }
    }
    /// Enable adaptive time-step control with the given target Mach number.
    pub fn enable_adaptive_dt(&mut self, ma_target: f64) {
        self.adaptive_dt = Some(AdaptiveDtController::new(ma_target, self.check_interval));
    }
    /// Run the simulation, returning `(steps_taken, converged)`.
    pub fn run(&mut self, sim: &mut LbmSimulation) -> (usize, bool) {
        let n_cells = sim.get_velocity_field().len();
        let mut monitor = ConvergenceMonitor::new(n_cells, self.tolerance, self.check_interval);
        let mut converged = false;
        for step in 0..self.max_steps {
            let dt = if let Some(ref adc) = self.adaptive_dt {
                adc.dt_current
            } else {
                1.0
            };
            sim.step(dt);
            if let Some(ref mut adc) = self.adaptive_dt
                && step % adc.adjust_interval == 0
            {
                let vel = sim.get_velocity_field();
                let u_max = vel
                    .iter()
                    .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
                    .fold(0.0_f64, f64::max);
                adc.update(u_max);
            }
            if step % self.check_interval == 0 {
                let vel = sim.get_velocity_field();
                let ux: Vec<f64> = vel.iter().map(|v| v[0]).collect();
                let uy: Vec<f64> = vel.iter().map(|v| v[1]).collect();
                converged = monitor.check(&ux, &uy);
                self.residuals.push(monitor.latest_residual());
                if step % (5 * self.check_interval) == 0 {
                    self.stats.record(sim);
                }
                if converged {
                    if self.checkpoint_interval > 0 {
                        self.checkpoint_store.save(sim);
                    }
                    return (step + 1, true);
                }
            }
            if self.checkpoint_interval > 0 && step % self.checkpoint_interval == 0 && step > 0 {
                self.checkpoint_store.save(sim);
            }
        }
        (self.max_steps, converged)
    }
    /// Return the last residual from the run.
    pub fn last_residual(&self) -> f64 {
        self.residuals.last().copied().unwrap_or(f64::MAX)
    }
    /// Whether the residuals are monotonically decreasing over the last window.
    pub fn is_converging(&self) -> bool {
        let n = self.residuals.len();
        if n < 3 {
            return true;
        }
        self.residuals[n - 1] < self.residuals[n - 3]
    }
}
/// High-level run loop that drives a simulation until convergence or max steps.
pub struct RunLoop {
    /// Maximum number of steps to run.
    pub max_steps: usize,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Interval between convergence checks.
    pub check_interval: usize,
    /// Statistics collected during the run.
    pub stats: SimulationStatistics,
    /// Convergence residuals recorded during the run.
    pub residuals: Vec<f64>,
}
impl RunLoop {
    /// Create a new run loop.
    pub fn new(max_steps: usize, tolerance: f64, check_interval: usize) -> Self {
        Self {
            max_steps,
            tolerance,
            check_interval,
            stats: SimulationStatistics::new(),
            residuals: Vec::new(),
        }
    }
    /// Run until convergence or max steps.
    ///
    /// Returns the number of steps actually taken and whether the simulation
    /// converged.
    pub fn run(&mut self, sim: &mut LbmSimulation) -> (usize, bool) {
        let n_cells = sim.get_velocity_field().len();
        let mut monitor = ConvergenceMonitor::new(n_cells, self.tolerance, self.check_interval);
        let mut converged = false;
        for step in 0..self.max_steps {
            sim.step(1.0);
            if step % self.check_interval == 0 {
                let vel = sim.get_velocity_field();
                let ux: Vec<f64> = vel.iter().map(|v| v[0]).collect();
                let uy: Vec<f64> = vel.iter().map(|v| v[1]).collect();
                converged = monitor.check(&ux, &uy);
                self.residuals.push(monitor.latest_residual());
                if step % (10 * self.check_interval) == 0 {
                    self.stats.record(sim);
                }
                if converged {
                    return (step + 1, true);
                }
            }
        }
        (self.max_steps, converged)
    }
    /// Return the last recorded residual.
    pub fn last_residual(&self) -> f64 {
        self.residuals.last().copied().unwrap_or(f64::MAX)
    }
}
/// Post-processing utilities for LBM velocity and pressure fields.
pub struct FlowFieldAnalysis;
impl FlowFieldAnalysis {
    /// Compute the L2 norm of the velocity field.
    pub fn velocity_l2_norm(vel: &[[f64; 3]]) -> f64 {
        let sum: f64 = vel
            .iter()
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum();
        (sum / vel.len().max(1) as f64).sqrt()
    }
    /// Compute the maximum velocity magnitude.
    pub fn max_velocity(vel: &[[f64; 3]]) -> f64 {
        vel.iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }
    /// Compute the mean velocity magnitude.
    pub fn mean_velocity(vel: &[[f64; 3]]) -> f64 {
        if vel.is_empty() {
            return 0.0;
        }
        let sum: f64 = vel
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .sum();
        sum / vel.len() as f64
    }
    /// Compute the relative L2 error between two velocity fields.
    pub fn relative_l2_error(vel_num: &[[f64; 3]], vel_ref: &[[f64; 3]]) -> f64 {
        assert_eq!(vel_num.len(), vel_ref.len(), "Field size mismatch");
        let mut err_sq = 0.0_f64;
        let mut ref_sq = 0.0_f64;
        for (vn, vr) in vel_num.iter().zip(vel_ref.iter()) {
            let dx = vn[0] - vr[0];
            let dy = vn[1] - vr[1];
            let dz = vn[2] - vr[2];
            err_sq += dx * dx + dy * dy + dz * dz;
            ref_sq += vr[0] * vr[0] + vr[1] * vr[1] + vr[2] * vr[2];
        }
        if ref_sq < 1e-30 {
            return err_sq.sqrt();
        }
        (err_sq / ref_sq).sqrt()
    }
    /// Compute the vorticity (z-component) at each interior cell of a 2D grid
    /// using central differences.
    ///
    /// `omega_z = duy/dx - dux/dy`
    ///
    /// Returns a flat vector of length `nx * ny`; boundary cells are set to 0.
    pub fn vorticity_z(vel: &[[f64; 3]], nx: usize, ny: usize) -> Vec<f64> {
        let mut vort = vec![0.0_f64; nx * ny];
        for y in 1..(ny - 1) {
            for x in 1..(nx - 1) {
                let k = y * nx + x;
                let ke = y * nx + (x + 1);
                let kw = y * nx + (x - 1);
                let kn = (y + 1) * nx + x;
                let ks = (y - 1) * nx + x;
                let duy_dx = (vel[ke][1] - vel[kw][1]) / 2.0;
                let dux_dy = (vel[kn][0] - vel[ks][0]) / 2.0;
                vort[k] = duy_dx - dux_dy;
            }
        }
        vort
    }
    /// Compute the divergence at each interior cell of a 2D grid.
    ///
    /// `div = dux/dx + duy/dy`
    pub fn divergence_2d(vel: &[[f64; 3]], nx: usize, ny: usize) -> Vec<f64> {
        let mut div = vec![0.0_f64; nx * ny];
        for y in 1..(ny - 1) {
            for x in 1..(nx - 1) {
                let k = y * nx + x;
                let dux_dx = (vel[y * nx + x + 1][0] - vel[y * nx + x - 1][0]) / 2.0;
                let duy_dy = (vel[(y + 1) * nx + x][1] - vel[(y - 1) * nx + x][1]) / 2.0;
                div[k] = dux_dx + duy_dy;
            }
        }
        div
    }
    /// Compute the kinetic energy per cell: `KE = 0.5 * rho * u^2`.
    pub fn kinetic_energy_field(vel: &[[f64; 3]], rho: &[f64]) -> Vec<f64> {
        vel.iter()
            .zip(rho.iter())
            .map(|(v, &r)| 0.5 * r * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .collect()
    }
    /// Compute the total kinetic energy (sum over all cells).
    pub fn total_kinetic_energy(vel: &[[f64; 3]], rho: &[f64]) -> f64 {
        Self::kinetic_energy_field(vel, rho).iter().sum()
    }
    /// Compute the mean pressure.
    pub fn mean_pressure(pressure: &[f64]) -> f64 {
        if pressure.is_empty() {
            return 0.0;
        }
        pressure.iter().sum::<f64>() / pressure.len() as f64
    }
    /// Compute the pressure variance.
    pub fn pressure_variance(pressure: &[f64]) -> f64 {
        if pressure.len() < 2 {
            return 0.0;
        }
        let mean = Self::mean_pressure(pressure);
        pressure.iter().map(|&p| (p - mean).powi(2)).sum::<f64>() / (pressure.len() - 1) as f64
    }
}
/// Enum holding either a 2D or a 3D computational grid.
#[derive(Debug, Clone)]
pub(super) enum SimGrid {
    Grid2D(LbmGrid2D),
    Grid3D(LbmGrid3D),
}
