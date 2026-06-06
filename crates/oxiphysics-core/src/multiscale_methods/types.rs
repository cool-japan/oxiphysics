//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Slow-fast system splitter for stiff ODEs.
///
/// Models `dy/dt = f_slow(y, z) + f_fast(y, z)` where `z = z(y)` at slow
/// equilibrium.  Implements the heterogeneous multiscale method for ODEs
/// (HMM-ODE): run the fast system to equilibrium, then advance the slow
/// variables.
pub struct SlowFastSplitter {
    /// Slow time step.
    pub dt_slow: f64,
    /// Fast time step.
    pub dt_fast: f64,
    /// Number of fast steps per slow step.
    pub fast_steps: usize,
}
impl SlowFastSplitter {
    /// Create a slow-fast splitter.
    pub fn new(dt_slow: f64, dt_fast: f64, fast_steps: usize) -> Self {
        Self {
            dt_slow,
            dt_fast,
            fast_steps,
        }
    }
    /// Advance one slow step using averaged fast dynamics.
    ///
    /// * `slow`      – Slow variable state.
    /// * `fast`      – Fast variable state (input/output).
    /// * `slow_rhs`  – Returns `f_slow(slow, fast)`.
    /// * `fast_rhs`  – Returns `f_fast(slow, fast)`.
    ///
    /// Returns updated slow state.
    pub fn step<FS, FF>(
        &self,
        slow: &[f64],
        fast: &mut [f64],
        slow_rhs: FS,
        fast_rhs: FF,
    ) -> Vec<f64>
    where
        FS: Fn(&[f64], &[f64]) -> Vec<f64>,
        FF: Fn(&[f64], &[f64]) -> Vec<f64>,
    {
        let mut f_avg = vec![0.0f64; slow.len()];
        for step_i in 0..self.fast_steps {
            let ff = fast_rhs(slow, fast);
            for (fi, ffi) in fast.iter_mut().zip(ff.iter()) {
                *fi += self.dt_fast * ffi;
            }
            if step_i >= self.fast_steps / 2 {
                let fs = slow_rhs(slow, fast);
                for (a, si) in f_avg.iter_mut().zip(fs.iter()) {
                    *a += si;
                }
            }
        }
        let n_avg = (self.fast_steps - self.fast_steps / 2).max(1) as f64;
        slow.iter()
            .zip(f_avg.iter())
            .map(|(s, fa)| s + self.dt_slow * fa / n_avg)
            .collect()
    }
}
/// Coarse Projective Integration with RK4 extrapolation.
///
/// Instead of simple Euler extrapolation, this uses a Runge-Kutta 4
/// step on the coarse time derivative estimated from the micro solver.
pub struct CpiRk4 {
    /// Micro time step.
    pub dt_micro: f64,
    /// Number of micro steps for estimation.
    pub estim_steps: usize,
    /// Coarse (projective) time step.
    pub dt_coarse: f64,
}
impl CpiRk4 {
    /// Create a CPI-RK4 integrator.
    pub fn new(dt_micro: f64, estim_steps: usize, dt_coarse: f64) -> Self {
        Self {
            dt_micro,
            estim_steps,
            dt_coarse,
        }
    }
    /// Estimate the coarse time derivative by running the micro solver.
    fn estimate_derivative<M, R>(
        &self,
        coarse: &[f64],
        lift: &impl Fn(&[f64]) -> Vec<f64>,
        micro_step: &M,
        restrict: &R,
    ) -> Vec<f64>
    where
        M: Fn(&[f64]) -> Vec<f64>,
        R: Fn(&[f64]) -> Vec<f64>,
    {
        let mut fine = lift(coarse);
        let c0 = restrict(&fine);
        for _ in 0..self.estim_steps {
            fine = micro_step(&fine);
        }
        let c1 = restrict(&fine);
        let dt = self.dt_micro * self.estim_steps as f64;
        c0.iter()
            .zip(c1.iter())
            .map(|(a, b)| (b - a) / dt)
            .collect()
    }
    /// Take one coarse RK4 step.
    ///
    /// * `coarse`    – Current coarse state.
    /// * `lift`      – Lift operator.
    /// * `micro_step` – Fine-scale stepper.
    /// * `restrict`  – Restrict operator.
    ///
    /// Returns updated coarse state.
    pub fn step<M, R>(
        &self,
        coarse: &[f64],
        lift: &impl Fn(&[f64]) -> Vec<f64>,
        micro_step: &M,
        restrict: &R,
    ) -> Vec<f64>
    where
        M: Fn(&[f64]) -> Vec<f64>,
        R: Fn(&[f64]) -> Vec<f64>,
    {
        let h = self.dt_coarse;
        let k1 = self.estimate_derivative(coarse, lift, micro_step, restrict);
        let s2: Vec<f64> = coarse
            .iter()
            .zip(k1.iter())
            .map(|(c, k)| c + 0.5 * h * k)
            .collect();
        let k2 = self.estimate_derivative(&s2, lift, micro_step, restrict);
        let s3: Vec<f64> = coarse
            .iter()
            .zip(k2.iter())
            .map(|(c, k)| c + 0.5 * h * k)
            .collect();
        let k3 = self.estimate_derivative(&s3, lift, micro_step, restrict);
        let s4: Vec<f64> = coarse
            .iter()
            .zip(k3.iter())
            .map(|(c, k)| c + h * k)
            .collect();
        let k4 = self.estimate_derivative(&s4, lift, micro_step, restrict);
        coarse
            .iter()
            .enumerate()
            .map(|(i, &c)| c + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            .collect()
    }
}
/// Result of a homogenization cell-problem solve.
#[derive(Debug, Clone)]
pub struct HomogenizationResult {
    /// Effective isotropic property (scalar for simplicity).
    pub effective_property: f64,
    /// Voigt upper bound.
    pub voigt_bound: f64,
    /// Reuss lower bound.
    pub reuss_bound: f64,
    /// Hashin-Shtrikman upper bound (two-phase).
    pub hs_upper: f64,
    /// Hashin-Shtrikman lower bound (two-phase).
    pub hs_lower: f64,
    /// Volume fraction of phase 1.
    pub volume_fraction: f64,
}
/// A coarse grid for sequential multiscale upscaling.
///
/// Each coarse cell aggregates a block of fine-grid cells.
pub struct CoarseGrid {
    /// Number of coarse cells.
    pub n_coarse: usize,
    /// Upscaling ratio (number of fine cells per coarse cell).
    pub ratio: usize,
    /// Fine-grid cell values.
    pub fine_values: Vec<f64>,
    /// Coarse-grid cell values (computed by upscaling).
    pub coarse_values: Vec<f64>,
}
impl CoarseGrid {
    /// Create a coarse grid and immediately compute coarse values by arithmetic upscaling.
    pub fn from_fine(fine_values: Vec<f64>, ratio: usize) -> Self {
        let n_fine = fine_values.len();
        let n_coarse = n_fine.div_ceil(ratio);
        let mut coarse_values = vec![0.0f64; n_coarse];
        for (ci, chunk) in fine_values.chunks(ratio).enumerate() {
            coarse_values[ci] = chunk.iter().sum::<f64>() / chunk.len() as f64;
        }
        Self {
            n_coarse,
            ratio,
            fine_values,
            coarse_values,
        }
    }
    /// Downscale coarse values back to the fine grid using piecewise-constant interpolation.
    pub fn downscale(&self) -> Vec<f64> {
        let mut fine = Vec::with_capacity(self.fine_values.len());
        for (ci, &cv) in self.coarse_values.iter().enumerate() {
            let start = ci * self.ratio;
            let end = (start + self.ratio).min(self.fine_values.len());
            for _ in start..end {
                fine.push(cv);
            }
        }
        fine
    }
    /// Downscale with linear interpolation between coarse-cell centres.
    pub fn downscale_linear(&self) -> Vec<f64> {
        let n_fine = self.fine_values.len();
        let mut fine = vec![0.0f64; n_fine];
        for (fi, fv) in fine.iter_mut().enumerate() {
            let ci_f = fi as f64 / self.ratio as f64;
            let ci0 = (ci_f as usize).min(self.n_coarse - 1);
            let ci1 = (ci0 + 1).min(self.n_coarse - 1);
            let alpha = ci_f - ci0 as f64;
            *fv = (1.0 - alpha) * self.coarse_values[ci0] + alpha * self.coarse_values[ci1];
        }
        fine
    }
}
/// A rectangular subregion used in concurrent multiscale domain decomposition.
#[derive(Debug, Clone)]
pub struct SubDomain {
    /// Domain index.
    pub id: usize,
    /// Lower-left corner.
    pub min: [f64; 3],
    /// Upper-right corner.
    pub max: [f64; 3],
    /// True if this domain uses the fine (micro) solver.
    pub is_fine: bool,
}
impl SubDomain {
    /// Create a new subdomain.
    pub fn new(id: usize, min: [f64; 3], max: [f64; 3], is_fine: bool) -> Self {
        Self {
            id,
            min,
            max,
            is_fine,
        }
    }
    /// Volume of the subdomain.
    pub fn volume(&self) -> f64 {
        let d = sub3(self.max, self.min);
        d[0] * d[1] * d[2]
    }
    /// Centre of the subdomain.
    pub fn centre(&self) -> [f64; 3] {
        scale3(add3(self.min, self.max), 0.5)
    }
    /// Return `true` if the point `p` lies inside this subdomain.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
    /// Return `true` if this subdomain overlaps `other`.
    pub fn overlaps(&self, other: &SubDomain) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}
/// Resolution level for adaptive multiscale switching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolutionLevel {
    /// Coarse (macro-scale) solver.
    Coarse,
    /// Intermediate (meso-scale) solver.
    Meso,
    /// Fine (micro-scale) solver.
    Fine,
}
/// Criterion for triggering a resolution change.
#[derive(Debug, Clone)]
pub struct AdaptiveCriterion {
    /// Error threshold above which we switch to a finer level.
    pub refine_threshold: f64,
    /// Error threshold below which we switch to a coarser level.
    pub coarsen_threshold: f64,
    /// Maximum gradient magnitude before forced refinement.
    pub gradient_threshold: f64,
}
impl AdaptiveCriterion {
    /// Create a new adaptive criterion with given thresholds.
    pub fn new(refine_threshold: f64, coarsen_threshold: f64, gradient_threshold: f64) -> Self {
        Self {
            refine_threshold,
            coarsen_threshold,
            gradient_threshold,
        }
    }
    /// Determine the appropriate resolution for the given local error and gradient.
    pub fn decide(&self, error: f64, gradient: f64, current: ResolutionLevel) -> ResolutionLevel {
        if error > self.refine_threshold || gradient > self.gradient_threshold {
            match current {
                ResolutionLevel::Coarse => ResolutionLevel::Meso,
                ResolutionLevel::Meso => ResolutionLevel::Fine,
                ResolutionLevel::Fine => ResolutionLevel::Fine,
            }
        } else if error < self.coarsen_threshold && gradient < self.coarsen_threshold {
            match current {
                ResolutionLevel::Fine => ResolutionLevel::Meso,
                ResolutionLevel::Meso => ResolutionLevel::Coarse,
                ResolutionLevel::Coarse => ResolutionLevel::Coarse,
            }
        } else {
            current
        }
    }
}
/// Adaptive multiscale mesh node.
#[derive(Debug, Clone)]
pub struct AdaptiveNode {
    /// Position of this node.
    pub position: [f64; 3],
    /// Current field value at this node.
    pub value: f64,
    /// Current resolution level.
    pub level: ResolutionLevel,
    /// Local error estimate.
    pub error: f64,
}
impl AdaptiveNode {
    /// Create a new adaptive node.
    pub fn new(position: [f64; 3], value: f64) -> Self {
        Self {
            position,
            value,
            level: ResolutionLevel::Coarse,
            error: 0.0,
        }
    }
    /// Update the error estimate (simple: |Δvalue| / max(|value|, 1e-10)).
    pub fn update_error(&mut self, prev_value: f64) {
        self.error = (self.value - prev_value).abs() / self.value.abs().max(1e-10);
    }
}
/// Configuration for the Heterogeneous Multiscale Method.
#[derive(Debug, Clone)]
pub struct HmmConfig {
    /// Number of macro quadrature points.
    pub n_macro_points: usize,
    /// Number of micro time steps per macro step.
    pub micro_steps: usize,
    /// Micro time step size.
    pub dt_micro: f64,
    /// Macro time step size.
    pub dt_macro: f64,
    /// Size of the micro domain (unit-cell side length).
    pub micro_domain_size: f64,
}
/// Discretised periodic unit-cell used in homogenization theory.
///
/// The cell is divided into `nx × ny × nz` voxels.  Each voxel carries a
/// local material tensor (isotropic here, stored as a scalar conductivity /
/// Young's modulus).
pub struct UnitCell {
    /// Number of voxels along x.
    pub nx: usize,
    /// Number of voxels along y.
    pub ny: usize,
    /// Number of voxels along z.
    pub nz: usize,
    /// Physical side length of one voxel.
    pub h: f64,
    /// Local isotropic property (e.g., conductivity) at each voxel.
    pub property: Vec<f64>,
}
impl UnitCell {
    /// Create a uniform unit cell where every voxel has `value`.
    pub fn uniform(nx: usize, ny: usize, nz: usize, h: f64, value: f64) -> Self {
        let n = nx * ny * nz;
        Self {
            nx,
            ny,
            nz,
            h,
            property: vec![value; n],
        }
    }
    /// Create a two-phase checkerboard unit cell.
    ///
    /// Even-sum indices get `v0`, odd-sum indices get `v1`.
    pub fn checkerboard(nx: usize, ny: usize, nz: usize, h: f64, v0: f64, v1: f64) -> Self {
        let mut prop = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    prop.push(if (i + j + k) % 2 == 0 { v0 } else { v1 });
                }
            }
        }
        Self {
            nx,
            ny,
            nz,
            h,
            property: prop,
        }
    }
    /// Return the linear voxel index.
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }
    /// Arithmetic mean of the local property (Voigt upper bound for
    /// conductivity homogenization).
    pub fn voigt_average(&self) -> f64 {
        let sum: f64 = self.property.iter().sum();
        sum / (self.property.len() as f64)
    }
    /// Harmonic mean of the local property (Reuss lower bound).
    pub fn reuss_average(&self) -> f64 {
        let sum: f64 = self.property.iter().map(|&p| 1.0 / p.max(1e-300)).sum();
        (self.property.len() as f64) / sum
    }
    /// Geometric mean of the local property.
    pub fn geometric_average(&self) -> f64 {
        let log_sum: f64 = self.property.iter().map(|&p| p.max(1e-300).ln()).sum();
        (log_sum / self.property.len() as f64).exp()
    }
}
/// Multi-level error estimator that tracks errors across scale levels.
pub struct MultilevelErrorEstimator {
    /// Error tolerance per level.
    pub tolerances: Vec<f64>,
    /// Recorded errors per level.
    pub errors: Vec<f64>,
}
impl MultilevelErrorEstimator {
    /// Create a new multilevel estimator with uniform tolerance `tol` for `n_levels` levels.
    pub fn new(n_levels: usize, tol: f64) -> Self {
        Self {
            tolerances: vec![tol; n_levels],
            errors: vec![0.0; n_levels],
        }
    }
    /// Record the error at a given level.
    pub fn record(&mut self, level: usize, error: f64) {
        if level < self.errors.len() {
            self.errors[level] = error;
        }
    }
    /// Return `true` if all errors are within tolerance.
    pub fn all_converged(&self) -> bool {
        self.errors
            .iter()
            .zip(self.tolerances.iter())
            .all(|(e, t)| e <= t)
    }
    /// Total error as the Euclidean norm of per-level errors.
    pub fn total_error(&self) -> f64 {
        self.errors.iter().map(|e| e * e).sum::<f64>().sqrt()
    }
}
/// Equation-free coarse-grained projective integration driver.
///
/// This implements the basic equation-free scheme:
/// 1. **Lift**: map coarse state to fine-scale initial condition.
/// 2. **Run**: advance the fine-scale model for a healing + estimation window.
/// 3. **Restrict**: coarsen the fine-scale state.
/// 4. **Project**: extrapolate coarse state over a large time step.
pub struct EquationFreeIntegrator {
    /// Time step for the fine-scale (microscale) solver.
    pub dt_micro: f64,
    /// Number of micro steps for the healing phase.
    pub heal_steps: usize,
    /// Number of micro steps used to estimate the coarse time derivative.
    pub estim_steps: usize,
    /// Projective time step (coarse step).
    pub dt_project: f64,
}
impl EquationFreeIntegrator {
    /// Create a new equation-free integrator.
    pub fn new(dt_micro: f64, heal_steps: usize, estim_steps: usize, dt_project: f64) -> Self {
        Self {
            dt_micro,
            heal_steps,
            estim_steps,
            dt_project,
        }
    }
    /// Perform one projective integration step.
    ///
    /// * `coarse_state` – Current coarse-grained state vector.
    /// * `lift`         – Closure that converts coarse state to fine state.
    /// * `micro_step`   – Closure that advances the fine state by `dt_micro`.
    /// * `restrict`     – Closure that converts fine state to coarse state.
    ///
    /// Returns the new coarse state after the projective step.
    pub fn step<F, M, R>(
        &self,
        coarse_state: &[f64],
        lift: F,
        micro_step: M,
        restrict: R,
    ) -> Vec<f64>
    where
        F: Fn(&[f64]) -> Vec<f64>,
        M: Fn(&[f64]) -> Vec<f64>,
        R: Fn(&[f64]) -> Vec<f64>,
    {
        let mut fine = lift(coarse_state);
        for _ in 0..self.heal_steps {
            fine = micro_step(&fine);
        }
        let c0 = restrict(&fine);
        for _ in 0..self.estim_steps {
            fine = micro_step(&fine);
        }
        let c1 = restrict(&fine);
        let dt_estim = self.dt_micro * self.estim_steps as f64;
        let dcoarse: Vec<f64> = c0
            .iter()
            .zip(c1.iter())
            .map(|(a, b)| (b - a) / dt_estim)
            .collect();
        c1.iter()
            .zip(dcoarse.iter())
            .map(|(c, dc)| c + dc * self.dt_project)
            .collect()
    }
}
/// A node on the nudged elastic band (NEB) for minimum energy path finding.
#[derive(Debug, Clone)]
pub struct NebImage {
    /// Configuration coordinates.
    pub coords: Vec<f64>,
    /// Energy at this image.
    pub energy: f64,
    /// Force on this image.
    pub force: Vec<f64>,
}
impl NebImage {
    /// Create a new NEB image.
    pub fn new(coords: Vec<f64>) -> Self {
        let n = coords.len();
        Self {
            coords,
            energy: 0.0,
            force: vec![0.0; n],
        }
    }
}
/// Domain decomposition for concurrent multiscale simulation.
pub struct ConcurrentDomainDecomposition {
    /// List of all subdomains.
    pub domains: Vec<SubDomain>,
}
impl ConcurrentDomainDecomposition {
    /// Create an empty domain decomposition.
    pub fn new() -> Self {
        Self {
            domains: Vec::new(),
        }
    }
    /// Add a subdomain.
    pub fn add_domain(&mut self, domain: SubDomain) {
        self.domains.push(domain);
    }
    /// Find which subdomain contains point `p`.  Returns `None` if no match.
    pub fn find_domain(&self, p: [f64; 3]) -> Option<usize> {
        self.domains.iter().position(|d| d.contains(p))
    }
    /// Return the overlap (handshake) pairs – subdomains that share an interface.
    pub fn overlap_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..self.domains.len() {
            for j in (i + 1)..self.domains.len() {
                if self.domains[i].overlaps(&self.domains[j]) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
    /// Interpolate a field value at position `p` from the coarse domain,
    /// using inverse-distance weighting from subdomain centres.
    pub fn interpolate_field(&self, p: [f64; 3], field_values: &[f64]) -> f64 {
        let mut weight_sum = 0.0f64;
        let mut value_sum = 0.0f64;
        for (i, dom) in self.domains.iter().enumerate() {
            if i >= field_values.len() {
                break;
            }
            let c = dom.centre();
            let d = norm3(sub3(p, c)).max(1e-12);
            let w = 1.0 / d;
            weight_sum += w;
            value_sum += w * field_values[i];
        }
        if weight_sum < 1e-300 {
            0.0
        } else {
            value_sum / weight_sum
        }
    }
}
/// A macro-scale quadrature point carrying effective data computed by
/// the micro solver.
#[derive(Debug, Clone)]
pub struct MacroPoint {
    /// Position in macro domain.
    pub position: f64,
    /// Effective flux at this point.
    pub effective_flux: f64,
    /// Effective stiffness tensor (isotropic 1D: scalar).
    pub effective_stiffness: f64,
}
/// Atom in the atomistic region.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Current position.
    pub pos: [f64; 3],
    /// Current velocity.
    pub vel: [f64; 3],
    /// Mass.
    pub mass: f64,
}
impl Atom {
    /// Create a new atom.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64) -> Self {
        Self { pos, vel, mass }
    }
    /// Kinetic energy of this atom.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
}
/// Phase-field model parameters for Allen-Cahn microstructure evolution.
#[derive(Debug, Clone)]
pub struct PhaseFieldParams {
    /// Interfacial energy parameter ε².
    pub epsilon_sq: f64,
    /// Relaxation coefficient M.
    pub mobility: f64,
    /// Double-well height W.
    pub well_height: f64,
}
impl PhaseFieldParams {
    /// Create a new set of phase-field parameters.
    pub fn new(epsilon_sq: f64, mobility: f64, well_height: f64) -> Self {
        Self {
            epsilon_sq,
            mobility,
            well_height,
        }
    }
    /// Interfacial width `η = ε √(8/W)`.
    pub fn interface_width(&self) -> f64 {
        self.epsilon_sq.sqrt() * (8.0 / self.well_height.max(1e-300)).sqrt()
    }
    /// Interface velocity scale `V₀ = M W`.
    pub fn velocity_scale(&self) -> f64 {
        self.mobility * self.well_height
    }
}
/// A coarse-grained bead representing a group of atoms.
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Position.
    pub pos: [f64; 3],
    /// Velocity.
    pub vel: [f64; 3],
    /// Force accumulator.
    pub force: [f64; 3],
    /// Effective mass.
    pub mass: f64,
}
impl CgBead {
    /// Create a new coarse-grained bead.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            mass,
        }
    }
    /// Reset the force accumulator to zero.
    pub fn reset_force(&mut self) {
        self.force = [0.0; 3];
    }
    /// Velocity-Verlet position half-step.
    pub fn vv_pos_step(&mut self, dt: f64) {
        self.pos = add3(self.pos, scale3(self.vel, dt));
    }
    /// Velocity-Verlet velocity update.
    pub fn vv_vel_update(&mut self, force_new: [f64; 3], dt: f64) {
        let acc_old = scale3(self.force, 1.0 / self.mass.max(1e-300));
        let acc_new = scale3(force_new, 1.0 / self.mass.max(1e-300));
        self.vel = add3(self.vel, scale3(add3(acc_old, acc_new), 0.5 * dt));
        self.force = force_new;
    }
}
/// Look-up table for pre-computed micro simulation data.
///
/// In information-passing multiscale methods the micro simulator is run
/// off-line for a set of macro inputs, and results are tabulated.
pub struct MicroDataTable {
    /// Macro inputs (e.g., strain values).
    pub inputs: Vec<f64>,
    /// Corresponding micro outputs (e.g., effective stress).
    pub outputs: Vec<f64>,
}
impl MicroDataTable {
    /// Create a new table from paired input/output vectors.
    pub fn new(inputs: Vec<f64>, outputs: Vec<f64>) -> Self {
        Self { inputs, outputs }
    }
    /// Query the table at macro input `x` using linear interpolation.
    pub fn query(&self, x: f64) -> f64 {
        let n = self.inputs.len().min(self.outputs.len());
        if n == 0 {
            return 0.0;
        }
        if n == 1 {
            return self.outputs[0];
        }
        if x <= self.inputs[0] {
            return self.outputs[0];
        }
        if x >= self.inputs[n - 1] {
            return self.outputs[n - 1];
        }
        let idx = self.inputs.partition_point(|&v| v <= x);
        let i0 = idx.saturating_sub(1);
        let i1 = i0 + 1;
        let alpha = (x - self.inputs[i0]) / (self.inputs[i1] - self.inputs[i0]).max(1e-300);
        (1.0 - alpha) * self.outputs[i0] + alpha * self.outputs[i1]
    }
    /// Build a table by sampling a micro model closure at `n` points in `[a, b]`.
    pub fn build<F: Fn(f64) -> f64>(a: f64, b: f64, n: usize, micro_fn: F) -> Self {
        let inputs: Vec<f64> = (0..n)
            .map(|i| a + (b - a) * i as f64 / (n - 1).max(1) as f64)
            .collect();
        let outputs: Vec<f64> = inputs.iter().map(|&x| micro_fn(x)).collect();
        Self::new(inputs, outputs)
    }
}
/// Richardson extrapolation error estimator.
///
/// Given solutions on a coarse grid (`h`) and a fine grid (`h/2`),
/// estimates the leading-order error and extrapolated solution.
#[derive(Debug, Clone)]
pub struct RichardsonEstimator {
    /// Order of the numerical method (e.g., 2 for second-order).
    pub order: f64,
}
impl RichardsonEstimator {
    /// Create a Richardson estimator for a method of given `order`.
    pub fn new(order: f64) -> Self {
        Self { order }
    }
    /// Compute Richardson-extrapolated value from coarse (`u_c`) and fine (`u_f`) solutions.
    ///
    /// `u_exact ≈ u_f + (u_f - u_c) / (2^order - 1)`
    pub fn extrapolate(&self, u_coarse: f64, u_fine: f64) -> f64 {
        let r = 2.0_f64.powf(self.order);
        u_fine + (u_fine - u_coarse) / (r - 1.0)
    }
    /// Estimate the error in the fine solution.
    pub fn error_estimate(&self, u_coarse: f64, u_fine: f64) -> f64 {
        let r = 2.0_f64.powf(self.order);
        (u_fine - u_coarse).abs() / (r - 1.0)
    }
    /// Estimate the effective order of convergence from three grid levels.
    ///
    /// Returns `p ≈ log2((u_c - u_m) / (u_m - u_f))`.
    pub fn effective_order(&self, u_coarse: f64, u_meso: f64, u_fine: f64) -> f64 {
        let num = (u_coarse - u_meso).abs();
        let den = (u_meso - u_fine).abs();
        if den < 1e-300 {
            self.order
        } else {
            (num / den).log2()
        }
    }
}
