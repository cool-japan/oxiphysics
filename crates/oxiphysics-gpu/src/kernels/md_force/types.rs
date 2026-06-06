//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{compute_lj_force, erfc_approx};

/// Mock GPU Ewald real-space short-range kernel.
///
/// Computes the real-space contribution to Ewald sum for all pairs within
/// `params.r_cutoff`. In a real GPU implementation this would be a WGSL
/// compute shader; here we provide a CPU reference implementation.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n, periodic wrapped)
///   - `inputs[1]`: charges `[q1, q2, ...]` (n values)
///   - `inputs[2]`: `[alpha, r_cutoff, box_length]` (3 values)
///
/// **Outputs:**
///   - `outputs[0]`: real-space forces `[fx, fy, fz, ...]` (3n)
///   - `outputs[1]`: `[real_space_energy]`
pub struct EwaldRealSpaceKernel;
/// Mock GPU virial stress tensor accumulation kernel.
///
/// Computes the pairwise virial stress tensor W_αβ = -Σ r_α F_β for all
/// LJ pairs within the cutoff. Returns all 6 independent tensor components.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: `[epsilon, sigma, cutoff]`
///
/// **Outputs:**
///   - `outputs[0]`: `[Wxx, Wyy, Wzz, Wxy, Wxz, Wyz]` (6 values)
pub struct VirialStressTensorKernel;
/// Kernel that computes Lennard-Jones pair forces and potential energy.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: `[epsilon, sigma, cutoff]` (3 values)
///
/// **Outputs:**
///   - `outputs[0]`: forces `[fx, fy, fz, ...]` (3n)
///   - `outputs[1]`: `[total_potential_energy]` (1 value)
pub struct LennardJonesKernel;
/// A pairwise force kernel parameterised by a Lennard-Jones potential, a
/// cutoff distance, and an optional energy shift.
pub struct PairForceKernel {
    /// Lennard-Jones parameters.
    pub lj: LjPotential,
    /// Cutoff distance (pairs beyond this are skipped).
    pub cutoff: f64,
    /// If `true`, shift the potential so V(cutoff) = 0.
    pub shift: bool,
}
impl PairForceKernel {
    /// Construct a `PairForceKernel`.
    pub fn new(lj: LjPotential, cutoff: f64, shift: bool) -> Self {
        Self { lj, cutoff, shift }
    }
    /// Evaluate the shifted potential correction at the cutoff.
    fn shift_energy(&self) -> f64 {
        if self.shift {
            let (e_cut, _) = compute_lj_force(self.cutoff, &self.lj);
            e_cut
        } else {
            0.0
        }
    }
    /// Compute force and energy for a pair at distance `r`.
    pub fn evaluate(&self, r: f64) -> (f64, f64) {
        if r >= self.cutoff || r < 1e-30 {
            return (0.0, 0.0);
        }
        let (e, f) = compute_lj_force(r, &self.lj);
        let e_shifted = if self.shift {
            e - self.shift_energy()
        } else {
            e
        };
        (e_shifted, f)
    }
}
/// A harmonic angle between three atoms i-j-k.
///
/// `V(θ) = 0.5 · k_theta · (θ - θ0)²`
#[derive(Debug, Clone, Copy)]
pub struct HarmonicAngle {
    /// Force constant (energy/rad²).
    pub k_theta: f64,
    /// Equilibrium angle in radians.
    pub theta0: f64,
    /// Atom index i (first end).
    pub atom_i: usize,
    /// Atom index j (central apex).
    pub atom_j: usize,
    /// Atom index k (second end).
    pub atom_k: usize,
}
impl HarmonicAngle {
    /// Create a new harmonic angle.
    pub fn new(atom_i: usize, atom_j: usize, atom_k: usize, k_theta: f64, theta0: f64) -> Self {
        Self {
            k_theta,
            theta0,
            atom_i,
            atom_j,
            atom_k,
        }
    }
}
/// Kernel that computes Coulomb pair forces and potential energy.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: charges `[q1, q2, ...]` (n values)
///   - `inputs[2]`: `[k_e, cutoff]` (2 values)
///
/// **Outputs:**
///   - `outputs[0]`: forces `[fx, fy, fz, ...]` (3n)
///   - `outputs[1]`: `[total_potential_energy]` (1 value)
pub struct CoulombKernel;
/// A harmonic bond between two atoms.
///
/// `V(r) = 0.5 · k · (r - r0)²`
/// `F = -k · (r - r0) · r̂`
#[derive(Debug, Clone, Copy)]
pub struct HarmonicBond {
    /// Spring constant (energy/length²).
    pub k: f64,
    /// Equilibrium bond length.
    pub r0: f64,
    /// Index of atom i.
    pub atom_i: usize,
    /// Index of atom j.
    pub atom_j: usize,
}
impl HarmonicBond {
    /// Create a new harmonic bond.
    pub fn new(atom_i: usize, atom_j: usize, k: f64, r0: f64) -> Self {
        Self {
            k,
            r0,
            atom_i,
            atom_j,
        }
    }
}
/// Mock GPU pair energy accumulation kernel.
///
/// Given a pre-built neighbor list (as pairs `[i, j, ...]`), computes
/// LJ pair energies and accumulates them per-particle.
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x, y, z, ...]` (3n)
///   - `inputs[1]`: pair list `[i, j, i, j, ...]` (2 * num_pairs)
///   - `inputs[2]`: `[epsilon, sigma, cutoff]`
///
/// **Outputs:**
///   - `outputs[0]`: per-particle energy (n values)
///   - `outputs[1]`: `[total_energy]`
pub struct PairEnergyAccumulateKernel;
/// Coulomb potential parameters.
///
/// V(r) = k_e · q_i · q_j / r
/// F(r) = k_e · q_i · q_j / r^2  (magnitude, positive = repulsive for same-sign charges)
#[derive(Debug, Clone, Copy)]
pub struct CoulombPotential {
    /// Coulomb constant (e.g. 332.06 for kcal/mol·Å·e² in real units).
    pub k_e: f64,
}
impl CoulombPotential {
    /// Construct a new `CoulombPotential`.
    pub fn new(k_e: f64) -> Self {
        Self { k_e }
    }
    /// Compute Coulomb energy and force magnitude for charges `qi`, `qj` at distance `r`.
    ///
    /// Returns `(energy, force_magnitude)`.
    /// Force is positive for repulsion (same-sign charges).
    pub fn compute(&self, qi: f64, qj: f64, r: f64) -> (f64, f64) {
        if r < 1e-30 {
            return (f64::INFINITY, f64::INFINITY);
        }
        let energy = self.k_e * qi * qj / r;
        let force_mag = self.k_e * qi * qj / (r * r);
        (energy, force_mag)
    }
}
/// Lennard-Jones potential parameters.
///
/// The LJ potential is: V(r) = 4·ε·\[(σ/r)^12 − (σ/r)^6\]
///
/// The equilibrium (minimum) distance is r_min = 2^(1/6)·σ.
#[derive(Debug, Clone, Copy)]
pub struct LjPotential {
    /// Well depth (energy units).
    pub epsilon: f64,
    /// Finite distance at which V(r) = 0 (length units).
    pub sigma: f64,
}
impl LjPotential {
    /// Construct a new `LjPotential`.
    pub fn new(epsilon: f64, sigma: f64) -> Self {
        Self { epsilon, sigma }
    }
    /// Distance at which the potential is minimum (force = 0).
    ///
    /// r_min = 2^(1/6) · σ
    pub fn r_min(&self) -> f64 {
        2.0_f64.powf(1.0 / 6.0) * self.sigma
    }
    /// Potential energy at the minimum: V(r_min) = -ε.
    pub fn well_depth(&self) -> f64 {
        -self.epsilon
    }
}
/// Cutoff scheme for pairwise interactions.
#[derive(Debug, Clone, Copy)]
pub enum CutoffScheme {
    /// Hard cutoff: force and energy are zero beyond cutoff.
    Hard {
        /// Cutoff distance.
        cutoff: f64,
    },
    /// Shifted potential: V(r) - V(rc).
    Shifted {
        /// Cutoff distance.
        cutoff: f64,
    },
    /// Switched: smooth force switch between `r_switch` and `r_cutoff`.
    Switched {
        /// Inner switch radius.
        r_switch: f64,
        /// Outer cutoff radius.
        r_cutoff: f64,
    },
}
impl CutoffScheme {
    /// Returns the cutoff distance for this scheme.
    pub fn cutoff_distance(&self) -> f64 {
        match *self {
            CutoffScheme::Hard { cutoff } => cutoff,
            CutoffScheme::Shifted { cutoff } => cutoff,
            CutoffScheme::Switched { r_cutoff, .. } => r_cutoff,
        }
    }
    /// Compute the switching function value at distance `r`.
    /// Returns 1.0 inside, 0.0 outside, smooth transition for Switched.
    pub fn switch_value(&self, r: f64) -> f64 {
        match *self {
            CutoffScheme::Hard { cutoff } => {
                if r < cutoff {
                    1.0
                } else {
                    0.0
                }
            }
            CutoffScheme::Shifted { cutoff } => {
                if r < cutoff {
                    1.0
                } else {
                    0.0
                }
            }
            CutoffScheme::Switched { r_switch, r_cutoff } => {
                if r <= r_switch {
                    1.0
                } else if r >= r_cutoff {
                    0.0
                } else {
                    let t = (r - r_switch) / (r_cutoff - r_switch);
                    1.0 - t * t * (3.0 - 2.0 * t)
                }
            }
        }
    }
}
/// A force buffer that accumulates forces for each particle.
pub struct ForceBuffer {
    /// Forces per particle: \[fx, fy, fz\] for each.
    pub forces: Vec<[f64; 3]>,
    /// Potential energy contributions per particle.
    pub energies: Vec<f64>,
    /// Virial tensor diagonal components per particle (for pressure computation).
    pub virial: Vec<[f64; 3]>,
}
impl ForceBuffer {
    /// Create a new force buffer for `n` particles, initialized to zero.
    pub fn new(n: usize) -> Self {
        Self {
            forces: vec![[0.0; 3]; n],
            energies: vec![0.0; n],
            virial: vec![[0.0; 3]; n],
        }
    }
    /// Zero out all forces, energies, and virials.
    pub fn clear(&mut self) {
        for f in &mut self.forces {
            *f = [0.0; 3];
        }
        for e in &mut self.energies {
            *e = 0.0;
        }
        for v in &mut self.virial {
            *v = [0.0; 3];
        }
    }
    /// Add a pairwise force contribution between particles i and j.
    ///
    /// `f_ij` is the force on particle i due to particle j (Newton III: j gets -f_ij).
    /// `e_ij` is the pairwise potential energy (split equally).
    /// `dx` is the displacement vector from j to i.
    pub fn add_pair(&mut self, i: usize, j: usize, f_ij: [f64; 3], e_ij: f64, dx: [f64; 3]) {
        for k in 0..3 {
            self.forces[i][k] += f_ij[k];
            self.forces[j][k] -= f_ij[k];
            let w = f_ij[k] * dx[k];
            self.virial[i][k] += w * 0.5;
            self.virial[j][k] += w * 0.5;
        }
        self.energies[i] += e_ij * 0.5;
        self.energies[j] += e_ij * 0.5;
    }
    /// Total potential energy.
    pub fn total_energy(&self) -> f64 {
        self.energies.iter().sum()
    }
    /// Total force (should be zero for an isolated system by Newton III).
    pub fn total_force(&self) -> [f64; 3] {
        let mut total = [0.0; 3];
        for f in &self.forces {
            total[0] += f[0];
            total[1] += f[1];
            total[2] += f[2];
        }
        total
    }
    /// Total virial (trace of virial tensor).
    pub fn total_virial(&self) -> f64 {
        self.virial.iter().map(|v| v[0] + v[1] + v[2]).sum()
    }
    /// Reduce forces across multiple buffers (for parallel accumulation).
    pub fn reduce_from(&mut self, others: &[ForceBuffer]) {
        for other in others {
            for i in 0..self.forces.len().min(other.forces.len()) {
                for k in 0..3 {
                    self.forces[i][k] += other.forces[i][k];
                    self.virial[i][k] += other.virial[i][k];
                }
                self.energies[i] += other.energies[i];
            }
        }
    }
}
/// Mock GPU PPPM charge-assignment kernel (NGP: nearest grid point).
///
/// Maps particle charges onto a 3-D mesh using the nearest grid point rule.
/// In production this would be a GPU shader; here we provide a CPU mock.
pub struct PppmChargeAssignKernel;
/// The full virial stress tensor (6 independent components for symmetric tensor).
///
/// Components are stored as `[xx, yy, zz, xy, xz, yz]`.
#[derive(Debug, Clone, Copy)]
pub struct VirialTensor {
    /// Components `[Wxx, Wyy, Wzz, Wxy, Wxz, Wyz]`.
    pub components: [f64; 6],
}
impl VirialTensor {
    /// Zero tensor.
    pub fn zero() -> Self {
        VirialTensor {
            components: [0.0; 6],
        }
    }
    /// Trace (Wxx + Wyy + Wzz).
    pub fn trace(&self) -> f64 {
        self.components[0] + self.components[1] + self.components[2]
    }
    /// Pressure contribution: -trace / (3 * volume).
    pub fn pressure_contribution(&self, volume: f64) -> f64 {
        if volume < 1e-30 {
            return 0.0;
        }
        -self.trace() / (3.0 * volume)
    }
    /// Add another virial tensor.
    pub fn add(&self, other: &VirialTensor) -> VirialTensor {
        let mut c = self.components;
        for (ck, &ok) in c.iter_mut().zip(other.components.iter()) {
            *ck += ok;
        }
        VirialTensor { components: c }
    }
}
/// GPU-style temperature scaling kernel (CPU mock).
///
/// **Inputs:**
///   - `inputs[0]`: velocities `[vx,vy,vz, ...]` (3n)
///   - `inputs[1]`: masses (n values)
///   - `inputs[2]`: `[t_target, k_boltzmann]`
///
/// **Outputs:**
///   - `outputs[0]`: scaled velocities `[vx,vy,vz, ...]` (3n)
///   - `outputs[1]`: `[t_before, t_after]`
pub struct TemperatureScaleKernel;
/// Simple cell-list based neighbor list for MD simulations.
pub struct NeighborList {
    /// For each particle, a list of neighbor particle indices.
    pub neighbors: Vec<Vec<usize>>,
    /// Cutoff distance used to build the list.
    pub cutoff: f64,
    /// Skin distance (list is valid while particles move less than skin/2).
    pub skin: f64,
}
impl NeighborList {
    /// Build a neighbor list from positions using a brute-force O(n²) scan.
    pub fn build_brute_force(positions: &[[f64; 3]], cutoff: f64, skin: f64) -> Self {
        let n = positions.len();
        let r_list = cutoff + skin;
        let r_list2 = r_list * r_list;
        let mut neighbors = vec![Vec::new(); n];
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < r_list2 {
                    neighbors[i].push(j);
                    neighbors[j].push(i);
                }
            }
        }
        Self {
            neighbors,
            cutoff,
            skin,
        }
    }
    /// Number of particles.
    pub fn num_particles(&self) -> usize {
        self.neighbors.len()
    }
    /// Total number of neighbor pairs (each counted once).
    pub fn num_pairs(&self) -> usize {
        let total: usize = self.neighbors.iter().map(|n| n.len()).sum();
        total / 2
    }
    /// Check if the list needs rebuilding given maximum particle displacement.
    pub fn needs_rebuild(&self, max_displacement: f64) -> bool {
        max_displacement > self.skin * 0.5
    }
}
/// GPU-style bond force kernel (CPU mock).
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x0,y0,z0, x1,y1,z1, ...]` (3n)
///   - `inputs[1]`: bond table `[i0, j0, k0, r0_0, i1, j1, k1, r0_1, ...]`
///     (4 values per bond: atom_i, atom_j, k, r0)
///
/// **Outputs:**
///   - `outputs[0]`: forces `[fx,fy,fz, ...]` (3n)
///   - `outputs[1]`: per-bond energies (one per bond)
pub struct BondForceKernel;
/// GPU neighbor list update status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NlistUpdateStatus {
    /// List is still valid; no rebuild required.
    Valid,
    /// List was rebuilt.
    Rebuilt,
}
/// Mock GPU neighbor list update kernel.
///
/// Checks if maximum displacement exceeds `skin/2`; if so, rebuilds the list.
///
/// **Inputs:**
///   - `inputs[0]`: current positions (3n)
///   - `inputs[1]`: reference positions at last rebuild (3n)
///   - `inputs[2]`: `[cutoff, skin]`
///
/// **Outputs:**
///   - `outputs[0]`: neighbor list as pairs `[i, j, i, j, ...]` (2*num_pairs)
///   - `outputs[1]`: `[status, num_pairs]` (2 values; status: 0=Valid, 1=Rebuilt)
pub struct NlistUpdateKernel;
/// GPU-style angle force kernel (CPU mock).
///
/// **Inputs:**
///   - `inputs[0]`: positions `[x0,y0,z0, ...]` (3n)
///   - `inputs[1]`: angle table `[i0,j0,k0,k_theta0,theta0_0, i1,j1,k1,k_theta1,theta0_1, ...]`
///     (5 values per angle)
///
/// **Outputs:**
///   - `outputs[0]`: forces `[fx,fy,fz, ...]` (3n)
///   - `outputs[1]`: per-angle energies
pub struct AngleForceKernel;
/// Parameters for Ewald summation of long-range electrostatics.
///
/// Splits Coulomb interaction into short-range (real space) and long-range
/// (reciprocal space) parts via a Gaussian splitting parameter α.
#[derive(Debug, Clone, Copy)]
pub struct EwaldParams {
    /// Splitting parameter α (1/length units).
    pub alpha: f64,
    /// Real-space cutoff.
    pub r_cutoff: f64,
    /// K-space cutoff squared (|k|² ≤ k_cutoff_sq).
    pub k_cutoff_sq: f64,
    /// Box length (cubic box assumed).
    pub box_length: f64,
}
impl EwaldParams {
    /// Construct `EwaldParams`.
    pub fn new(alpha: f64, r_cutoff: f64, k_cutoff_sq: f64, box_length: f64) -> Self {
        Self {
            alpha,
            r_cutoff,
            k_cutoff_sq,
            box_length,
        }
    }
    /// Estimate the real-space accuracy parameter erfc(α·r_cutoff).
    pub fn real_space_accuracy(&self) -> f64 {
        let x = self.alpha * self.r_cutoff;
        erfc_approx(x)
    }
}
/// Particle-mesh grid configuration for PPPM.
#[derive(Debug, Clone, Copy)]
pub struct PppmGrid {
    /// Number of grid points in each dimension.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Number of grid points in z.
    pub nz: usize,
    /// Box length.
    pub box_length: f64,
    /// Charge assignment order (1=NGP, 2=CIC, 3=TSC).
    pub order: usize,
}
impl PppmGrid {
    /// Create a new `PppmGrid`.
    pub fn new(nx: usize, ny: usize, nz: usize, box_length: f64, order: usize) -> Self {
        Self {
            nx,
            ny,
            nz,
            box_length,
            order,
        }
    }
    /// Grid spacing in x.
    pub fn dx(&self) -> f64 {
        self.box_length / self.nx as f64
    }
    /// Grid spacing in y.
    pub fn dy(&self) -> f64 {
        self.box_length / self.ny as f64
    }
    /// Grid spacing in z.
    pub fn dz(&self) -> f64 {
        self.box_length / self.nz as f64
    }
    /// Total number of grid points.
    pub fn total_points(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}
