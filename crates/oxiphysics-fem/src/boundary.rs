#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Boundary condition types and application methods.
//!
//! Supports Dirichlet (prescribed displacement), Neumann (applied force),
//! Robin (spring/elastic), periodic, symmetry, time-dependent, mixed, and
//! contact boundary conditions for finite element systems.
//! Dense-matrix variants are provided for small self-contained systems;
//! the sparse (`CsrMatrix`) variants are used in the main pipeline.

use crate::sparse::CsrMatrix;

// ── DirichletBc ───────────────────────────────────────────────────────────────

/// A Dirichlet (essential) boundary condition.
///
/// Prescribes the displacement at a specific degree of freedom.
#[derive(Debug, Clone)]
pub struct DirichletBc {
    /// Node index.
    pub node: usize,
    /// Degree of freedom: 0 = x, 1 = y, 2 = z.
    pub dof: usize,
    /// Prescribed displacement value.
    pub value: f64,
}

impl DirichletBc {
    /// Create a new Dirichlet boundary condition.
    pub fn new(node: usize, dof: usize, value: f64) -> Self {
        Self { node, dof, value }
    }

    /// Return the global DOF index for this boundary condition.
    pub fn global_dof(&self) -> usize {
        self.node * 3 + self.dof
    }
}

// ── NeumannBc ─────────────────────────────────────────────────────────────────

/// A Neumann (natural) boundary condition.
///
/// Applies a force vector at a specific node.
#[derive(Debug, Clone)]
pub struct NeumannBc {
    /// Node index.
    pub node: usize,
    /// Force vector components \[fx, fy, fz\].
    pub force: [f64; 3],
}

impl NeumannBc {
    /// Create a new Neumann boundary condition.
    pub fn new(node: usize, force: [f64; 3]) -> Self {
        Self { node, force }
    }
}

// ── DirichletBC (scalar DOF form) ─────────────────────────────────────────────

/// Dirichlet BC targeting an arbitrary scalar DOF index (not restricted to ×3 layout).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DirichletBC {
    /// Global node id.
    pub node_id: usize,
    /// DOF index within the node.
    pub dof: usize,
    /// Prescribed value.
    pub value: f64,
}

#[allow(dead_code)]
impl DirichletBC {
    /// Create a new scalar-DOF Dirichlet BC.
    pub fn new(node_id: usize, dof: usize, value: f64) -> Self {
        Self {
            node_id,
            dof,
            value,
        }
    }
}

// ── NeumannBC (scalar DOF form) ───────────────────────────────────────────────

/// Neumann BC applying a scalar force to a single DOF.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NeumannBC {
    /// Global node id.
    pub node_id: usize,
    /// DOF index within the node.
    pub dof: usize,
    /// Applied force.
    pub force: f64,
}

#[allow(dead_code)]
impl NeumannBC {
    /// Create a new scalar-DOF Neumann BC.
    pub fn new(node_id: usize, dof: usize, force: f64) -> Self {
        Self {
            node_id,
            dof,
            force,
        }
    }
}

// ── RobinBC ───────────────────────────────────────────────────────────────────

/// Robin (spring/elastic) boundary condition.
///
/// Adds a spring stiffness `stiffness` to the diagonal and a corresponding
/// force contribution `stiffness * value` to the RHS.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RobinBC {
    /// Global node id.
    pub node_id: usize,
    /// DOF index within the node.
    pub dof: usize,
    /// Spring stiffness.
    pub stiffness: f64,
    /// Prescribed reference value.
    pub value: f64,
}

#[allow(dead_code)]
impl RobinBC {
    /// Create a new Robin BC.
    pub fn new(node_id: usize, dof: usize, stiffness: f64, value: f64) -> Self {
        Self {
            node_id,
            dof,
            stiffness,
            value,
        }
    }

    /// Global DOF index.
    pub fn global_dof(&self, dofs_per_node: usize) -> usize {
        self.node_id * dofs_per_node + self.dof
    }
}

// ── RobinBCImproved ───────────────────────────────────────────────────────────

/// Improved Robin BC with separate convective and radiative components.
///
/// Models the general Robin condition:
///   α * u + β * du/dn = g
///
/// Applied to the system as:
///   K_diag\[i\] += α,  f\[i\] += g
///
/// The improved variant also supports a damping (velocity-dependent) term.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RobinBCImproved {
    /// Global node id.
    pub node_id: usize,
    /// DOF index within the node.
    pub dof: usize,
    /// Coefficient α (stiffness-like, multiplies u).
    pub alpha: f64,
    /// Coefficient β (flux-like, multiplies du/dn).
    pub beta: f64,
    /// RHS value g.
    pub g: f64,
    /// Optional damping coefficient (multiplies velocity).
    pub damping: f64,
}

#[allow(dead_code)]
impl RobinBCImproved {
    /// Create a new improved Robin BC.
    pub fn new(node_id: usize, dof: usize, alpha: f64, beta: f64, g: f64) -> Self {
        Self {
            node_id,
            dof,
            alpha,
            beta,
            g,
            damping: 0.0,
        }
    }

    /// Create with damping.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.damping = damping;
        self
    }

    /// Global DOF index.
    pub fn global_dof(&self, dofs_per_node: usize) -> usize {
        self.node_id * dofs_per_node + self.dof
    }

    /// Apply this improved Robin BC to a dense stiffness matrix and force vector.
    pub fn apply_dense(&self, k: &mut [f64], f: &mut [f64], n_dofs: usize, dofs_per_node: usize) {
        let i = self.global_dof(dofs_per_node);
        assert!(i < n_dofs, "Robin DOF {i} out of range {n_dofs}");
        k[i * n_dofs + i] += self.alpha;
        f[i] += self.g;
    }

    /// Apply damping contribution to a damping matrix diagonal.
    pub fn apply_damping(&self, c_diag: &mut [f64], dofs_per_node: usize) {
        let i = self.global_dof(dofs_per_node);
        if i < c_diag.len() {
            c_diag[i] += self.damping;
        }
    }
}

// ── PeriodicBC ────────────────────────────────────────────────────────────────

/// Periodic boundary condition linking two DOFs.
///
/// Enforces u\[master\] == u\[slave\] (optionally with an offset).
/// Applied via the penalty method or elimination.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PeriodicBC {
    /// Master node id.
    pub master_node: usize,
    /// Slave node id.
    pub slave_node: usize,
    /// DOF index within the node.
    pub dof: usize,
    /// Optional constant offset: u_slave = u_master + offset.
    pub offset: f64,
}

#[allow(dead_code)]
impl PeriodicBC {
    /// Create a new periodic BC with zero offset.
    pub fn new(master_node: usize, slave_node: usize, dof: usize) -> Self {
        Self {
            master_node,
            slave_node,
            dof,
            offset: 0.0,
        }
    }

    /// Create a periodic BC with a specified offset.
    pub fn with_offset(master_node: usize, slave_node: usize, dof: usize, offset: f64) -> Self {
        Self {
            master_node,
            slave_node,
            dof,
            offset,
        }
    }

    /// Master global DOF (assumes 3 DOF per node).
    pub fn master_global_dof(&self) -> usize {
        self.master_node * 3 + self.dof
    }

    /// Slave global DOF (assumes 3 DOF per node).
    pub fn slave_global_dof(&self) -> usize {
        self.slave_node * 3 + self.dof
    }
}

/// Apply periodic BCs to a dense stiffness matrix and force vector using penalty method.
///
/// Enforces u_slave = u_master + offset by adding penalty terms:
///   K\[s,s\] += penalty,  K\[m,m\] += penalty,  K\[s,m\] -= penalty,  K\[m,s\] -= penalty
///   f\[s\] += penalty * offset,  f\[m\] -= penalty * offset
#[allow(dead_code)]
pub fn apply_periodic_dense_penalty(
    k: &mut [f64],
    f: &mut [f64],
    n_dofs: usize,
    bcs: &[PeriodicBC],
    penalty: f64,
) {
    for bc in bcs {
        let m = bc.master_global_dof();
        let s = bc.slave_global_dof();
        assert!(m < n_dofs && s < n_dofs);

        k[s * n_dofs + s] += penalty;
        k[m * n_dofs + m] += penalty;
        k[s * n_dofs + m] -= penalty;
        k[m * n_dofs + s] -= penalty;
        f[s] += penalty * bc.offset;
        f[m] -= penalty * bc.offset;
    }
}

/// Apply periodic BCs by elimination (tie slave DOF to master DOF).
///
/// Row and column of slave DOF are modified so that the slave equation
/// becomes: u_slave - u_master = offset.
#[allow(dead_code)]
pub fn apply_periodic_dense_elimination(
    k: &mut [f64],
    f: &mut [f64],
    n_dofs: usize,
    bcs: &[PeriodicBC],
) {
    for bc in bcs {
        let m = bc.master_global_dof();
        let s = bc.slave_global_dof();
        assert!(m < n_dofs && s < n_dofs);

        // Transfer slave contributions to master
        for j in 0..n_dofs {
            if j != s {
                k[m * n_dofs + j] += k[s * n_dofs + j];
                k[j * n_dofs + m] += k[j * n_dofs + s];
            }
        }
        k[m * n_dofs + m] += k[s * n_dofs + s];
        f[m] += f[s];

        // Zero slave row and column, set constraint equation
        for j in 0..n_dofs {
            k[s * n_dofs + j] = 0.0;
            k[j * n_dofs + s] = 0.0;
        }
        k[s * n_dofs + s] = 1.0;
        k[s * n_dofs + m] = -1.0;
        f[s] = bc.offset;
    }
}

// ── SymmetryBC ────────────────────────────────────────────────────────────────

/// Symmetry boundary condition.
///
/// Enforces zero displacement normal to the symmetry plane.
/// The plane is defined by specifying which DOF is constrained to zero.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SymmetryBC {
    /// Node indices on the symmetry plane.
    pub nodes: Vec<usize>,
    /// DOF to constrain (0=x, 1=y, 2=z — normal to the symmetry plane).
    pub constrained_dof: usize,
}

#[allow(dead_code)]
impl SymmetryBC {
    /// Create a new symmetry BC.
    pub fn new(nodes: Vec<usize>, constrained_dof: usize) -> Self {
        assert!(constrained_dof < 3, "DOF must be 0, 1, or 2");
        Self {
            nodes,
            constrained_dof,
        }
    }

    /// Convert to a list of Dirichlet BCs (zero displacement on the constrained DOF).
    pub fn to_dirichlet_bcs(&self) -> Vec<DirichletBC> {
        self.nodes
            .iter()
            .map(|&node| DirichletBC::new(node, self.constrained_dof, 0.0))
            .collect()
    }

    /// Apply symmetry BC directly to a dense system.
    pub fn apply_dense(&self, k: &mut [f64], f: &mut [f64], n_dofs: usize) {
        let bcs = self.to_dirichlet_bcs();
        apply_dirichlet_dense(k, f, n_dofs, &bcs);
    }
}

// ── TimeDependentBC ───────────────────────────────────────────────────────────

/// A time-dependent Dirichlet boundary condition.
///
/// The prescribed value is a function of time: u(t) = amplitude * f(t)
/// where f(t) is one of several built-in time functions.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum TimeFunction {
    /// Constant value.
    Constant(f64),
    /// Linear ramp: f(t) = t / ramp_time, clamped to \[0, 1\].
    Ramp {
        /// Time over which the ramp reaches its full value.
        ramp_time: f64,
    },
    /// Sinusoidal: f(t) = sin(2π * frequency * t).
    Sinusoidal {
        /// Oscillation frequency in Hz.
        frequency: f64,
    },
    /// Step function: f(t) = 0 for t < step_time, 1 for t >= step_time.
    Step {
        /// Time at which the step occurs.
        step_time: f64,
    },
    /// Triangular wave with given period.
    Triangular {
        /// Period of the triangular wave.
        period: f64,
    },
    /// Exponential decay: f(t) = exp(-decay_rate * t).
    ExponentialDecay {
        /// Rate of exponential decay.
        decay_rate: f64,
    },
}

#[allow(dead_code)]
impl TimeFunction {
    /// Evaluate the time function at time `t`.
    pub fn evaluate(&self, t: f64) -> f64 {
        match self {
            TimeFunction::Constant(c) => *c,
            TimeFunction::Ramp { ramp_time } => {
                if *ramp_time <= 0.0 {
                    1.0
                } else {
                    (t / ramp_time).clamp(0.0, 1.0)
                }
            }
            TimeFunction::Sinusoidal { frequency } => {
                (2.0 * std::f64::consts::PI * frequency * t).sin()
            }
            TimeFunction::Step { step_time } => {
                if t >= *step_time {
                    1.0
                } else {
                    0.0
                }
            }
            TimeFunction::Triangular { period } => {
                if *period <= 0.0 {
                    return 0.0;
                }
                let phase = (t / period) % 1.0;
                if phase < 0.5 {
                    4.0 * phase - 1.0
                } else {
                    3.0 - 4.0 * phase
                }
            }
            TimeFunction::ExponentialDecay { decay_rate } => (-decay_rate * t).exp(),
        }
    }
}

/// Time-dependent Dirichlet BC.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TimeDependentBC {
    /// Node id.
    pub node_id: usize,
    /// DOF index.
    pub dof: usize,
    /// Amplitude multiplier.
    pub amplitude: f64,
    /// Time function.
    pub time_function: TimeFunction,
}

#[allow(dead_code)]
impl TimeDependentBC {
    /// Create a new time-dependent BC.
    pub fn new(node_id: usize, dof: usize, amplitude: f64, time_function: TimeFunction) -> Self {
        Self {
            node_id,
            dof,
            amplitude,
            time_function,
        }
    }

    /// Evaluate the prescribed value at time `t`.
    pub fn value_at(&self, t: f64) -> f64 {
        self.amplitude * self.time_function.evaluate(t)
    }

    /// Convert to a standard DirichletBC at a given time.
    pub fn to_dirichlet_bc(&self, t: f64) -> DirichletBC {
        DirichletBC::new(self.node_id, self.dof, self.value_at(t))
    }
}

/// Apply time-dependent Dirichlet BCs at a given time instant.
#[allow(dead_code)]
pub fn apply_time_dependent_dense(
    k: &mut [f64],
    f: &mut [f64],
    n_dofs: usize,
    bcs: &[TimeDependentBC],
    t: f64,
) {
    let dirichlet_bcs: Vec<DirichletBC> = bcs.iter().map(|bc| bc.to_dirichlet_bc(t)).collect();
    apply_dirichlet_dense(k, f, n_dofs, &dirichlet_bcs);
}

/// Apply time-dependent Dirichlet BCs to a sparse system.
#[allow(dead_code)]
pub fn apply_time_dependent_sparse(
    k: &mut CsrMatrix,
    f: &mut [f64],
    bcs: &[TimeDependentBC],
    t: f64,
) {
    let dirichlet_bcs: Vec<DirichletBc> = bcs
        .iter()
        .map(|bc| DirichletBc::new(bc.node_id, bc.dof, bc.value_at(t)))
        .collect();
    apply_dirichlet(k, f, &dirichlet_bcs);
}

// ── ContactBC ─────────────────────────────────────────────────────────────────

/// Penalty-based contact boundary condition between pairs of nodes.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactBC {
    /// Pairs of (node_a, node_b) in contact.
    pub node_pairs: Vec<(usize, usize)>,
    /// Initial gap distance (negative = penetration).
    pub gap: f64,
    /// Penalty stiffness.
    pub penalty: f64,
}

#[allow(dead_code)]
impl ContactBC {
    /// Create a new contact BC.
    pub fn new(node_pairs: Vec<(usize, usize)>, gap: f64, penalty: f64) -> Self {
        Self {
            node_pairs,
            gap,
            penalty,
        }
    }
}

// ── MixedBC ──────────────────────────────────────────────────────────────────

/// Mixed boundary condition combining Dirichlet on some DOFs and Neumann on others.
///
/// For example, a roller support constrains displacement in one direction
/// while allowing free movement (zero force) in others.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MixedBC {
    /// Node id.
    pub node_id: usize,
    /// Dirichlet DOFs and their prescribed values.
    pub dirichlet_dofs: Vec<(usize, f64)>,
    /// Neumann DOFs and their prescribed forces.
    pub neumann_dofs: Vec<(usize, f64)>,
}

#[allow(dead_code)]
impl MixedBC {
    /// Create a new mixed BC.
    pub fn new(
        node_id: usize,
        dirichlet_dofs: Vec<(usize, f64)>,
        neumann_dofs: Vec<(usize, f64)>,
    ) -> Self {
        Self {
            node_id,
            dirichlet_dofs,
            neumann_dofs,
        }
    }

    /// Create a roller support (fix one DOF, free others).
    pub fn roller(node_id: usize, fixed_dof: usize, fixed_value: f64) -> Self {
        let dirichlet_dofs = vec![(fixed_dof, fixed_value)];
        Self {
            node_id,
            dirichlet_dofs,
            neumann_dofs: Vec::new(),
        }
    }

    /// Create a pin support (fix two DOFs, free the third).
    pub fn pin_2d(node_id: usize, free_dof: usize) -> Self {
        let mut dirichlet_dofs = Vec::new();
        for d in 0..3 {
            if d != free_dof {
                dirichlet_dofs.push((d, 0.0));
            }
        }
        Self {
            node_id,
            dirichlet_dofs,
            neumann_dofs: Vec::new(),
        }
    }

    /// Convert to separate Dirichlet and Neumann BC lists.
    pub fn to_separate_bcs(&self) -> (Vec<DirichletBC>, Vec<NeumannBC>) {
        let dirichlet: Vec<DirichletBC> = self
            .dirichlet_dofs
            .iter()
            .map(|&(dof, val)| DirichletBC::new(self.node_id, dof, val))
            .collect();
        let neumann: Vec<NeumannBC> = self
            .neumann_dofs
            .iter()
            .map(|&(dof, force)| NeumannBC::new(self.node_id, dof, force))
            .collect();
        (dirichlet, neumann)
    }

    /// Apply this mixed BC to a dense system.
    pub fn apply_dense(&self, k: &mut [f64], f: &mut [f64], n_dofs: usize) {
        let (dirichlet, neumann) = self.to_separate_bcs();
        apply_neumann_scalar(f, n_dofs, &neumann);
        apply_dirichlet_dense(k, f, n_dofs, &dirichlet);
    }
}

// ── Dense matrix BC helpers ───────────────────────────────────────────────────

/// Apply Dirichlet BCs to a **dense** `n_dofs × n_dofs` row-major stiffness
/// matrix `K` (length `n_dofs * n_dofs`) and force vector `f` using the
/// penalty method.
///
/// For each prescribed DOF `i` with value `u_i`:
/// - `K[i*n + i]` is increased by `penalty`
/// - `f[i]` is increased by `penalty * u_i`
///
/// A large `penalty` (e.g. 1e30 × max diagonal) enforces the BC while keeping
/// the matrix non-singular.
#[allow(dead_code)]
pub fn apply_dirichlet_dense_penalty(
    k: &mut Vec<f64>,
    f: &mut Vec<f64>,
    n_dofs: usize,
    bcs: &[DirichletBC],
    penalty: f64,
) {
    for bc in bcs {
        let i = bc.node_id * 3 + bc.dof; // assumes 3 DOF/node; generalize if needed
        assert!(i < n_dofs, "BC DOF {i} out of range {n_dofs}");
        k[i * n_dofs + i] += penalty;
        f[i] += penalty * bc.value;
    }
}

/// Apply Dirichlet BCs to a dense matrix using the **elimination** (row/column
/// zeroing) method.
///
/// For each prescribed DOF `i` with value `u_i`:
/// - `f[j] -= K[j,i] * u_i` for all `j != i`
/// - Row `i` is zeroed and `K[i,i]` set to 1
/// - Column `i` is zeroed
/// - `f[i] = u_i`
#[allow(dead_code)]
pub fn apply_dirichlet_dense(k: &mut [f64], f: &mut [f64], n_dofs: usize, bcs: &[DirichletBC]) {
    for bc in bcs {
        let dof = bc.node_id * 3 + bc.dof;
        assert!(dof < n_dofs);

        // Adjust RHS
        if bc.value.abs() > 0.0 {
            for j in 0..n_dofs {
                if j != dof {
                    f[j] -= k[j * n_dofs + dof] * bc.value;
                }
            }
        }

        // Zero row and column, set diagonal to 1
        for j in 0..n_dofs {
            k[dof * n_dofs + j] = 0.0;
            k[j * n_dofs + dof] = 0.0;
        }
        k[dof * n_dofs + dof] = 1.0;
        f[dof] = bc.value;
    }
}

/// Apply Neumann BCs to a dense force vector.
#[allow(dead_code)]
pub fn apply_neumann_scalar(f: &mut [f64], n_dofs: usize, bcs: &[NeumannBC]) {
    for bc in bcs {
        let i = bc.node_id * 3 + bc.dof;
        assert!(i < n_dofs, "Neumann DOF {i} out of range");
        f[i] += bc.force;
    }
}

/// Apply Robin BCs to a dense stiffness matrix and force vector.
#[allow(dead_code)]
pub fn apply_robin(k: &mut [f64], f: &mut [f64], n_dofs: usize, bcs: &[RobinBC]) {
    for bc in bcs {
        let i = bc.global_dof(3);
        assert!(i < n_dofs, "Robin DOF {i} out of range");
        k[i * n_dofs + i] += bc.stiffness;
        f[i] += bc.stiffness * bc.value;
    }
}

/// Apply improved Robin BCs to a dense stiffness matrix and force vector.
#[allow(dead_code)]
pub fn apply_robin_improved(
    k: &mut [f64],
    f: &mut [f64],
    n_dofs: usize,
    bcs: &[RobinBCImproved],
    dofs_per_node: usize,
) {
    for bc in bcs {
        bc.apply_dense(k, f, n_dofs, dofs_per_node);
    }
}

// ── Sparse BC application ─────────────────────────────────────────────────────

/// Apply Dirichlet BCs to a sparse CSR matrix using the penalty method.
///
/// For each prescribed DOF `i` with value `u_i`:
/// - The diagonal entry `K[i,i]` is increased by `penalty`
/// - `f[i]` is increased by `penalty * u_i`
#[allow(dead_code)]
pub fn apply_dirichlet_sparse_penalty(
    k: &mut CsrMatrix,
    f: &mut [f64],
    bcs: &[DirichletBc],
    penalty: f64,
) {
    for bc in bcs {
        let dof = bc.global_dof();
        assert!(dof < k.nrows, "BC DOF {dof} out of range");

        // Find diagonal entry and add penalty
        let start = k.row_ptr[dof];
        let end = k.row_ptr[dof + 1];
        for idx in start..end {
            if k.col_indices[idx] == dof {
                k.values[idx] += penalty;
                break;
            }
        }
        f[dof] += penalty * bc.value;
    }
}

/// Apply Robin BCs to a sparse CSR matrix and force vector.
#[allow(dead_code)]
pub fn apply_robin_sparse(k: &mut CsrMatrix, f: &mut [f64], bcs: &[RobinBC], dofs_per_node: usize) {
    for bc in bcs {
        let i = bc.global_dof(dofs_per_node);
        assert!(i < k.nrows, "Robin DOF {i} out of range");

        // Find diagonal entry and add stiffness
        let start = k.row_ptr[i];
        let end = k.row_ptr[i + 1];
        for idx in start..end {
            if k.col_indices[idx] == i {
                k.values[idx] += bc.stiffness;
                break;
            }
        }
        f[i] += bc.stiffness * bc.value;
    }
}

/// Apply Neumann BCs to a sparse system force vector (scalar DOF form).
#[allow(dead_code)]
pub fn apply_neumann_sparse(f: &mut [f64], bcs: &[NeumannBC]) {
    for bc in bcs {
        let i = bc.node_id * 3 + bc.dof;
        if i < f.len() {
            f[i] += bc.force;
        }
    }
}

// ── BoundarySet ───────────────────────────────────────────────────────────────

/// Collection of all boundary conditions for a FEM system.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct BoundarySet {
    dirichlet: Vec<DirichletBC>,
    neumann: Vec<NeumannBC>,
    robin: Vec<RobinBC>,
    periodic: Vec<PeriodicBC>,
    symmetry: Vec<SymmetryBC>,
    time_dependent: Vec<TimeDependentBC>,
    mixed: Vec<MixedBC>,
}

#[allow(dead_code)]
impl BoundarySet {
    /// Create an empty boundary set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a Dirichlet BC.
    pub fn add_dirichlet(&mut self, bc: DirichletBC) {
        self.dirichlet.push(bc);
    }

    /// Add a Neumann BC.
    pub fn add_neumann(&mut self, bc: NeumannBC) {
        self.neumann.push(bc);
    }

    /// Add a Robin BC.
    pub fn add_robin(&mut self, bc: RobinBC) {
        self.robin.push(bc);
    }

    /// Add a periodic BC.
    pub fn add_periodic(&mut self, bc: PeriodicBC) {
        self.periodic.push(bc);
    }

    /// Add a symmetry BC.
    pub fn add_symmetry(&mut self, bc: SymmetryBC) {
        self.symmetry.push(bc);
    }

    /// Add a time-dependent BC.
    pub fn add_time_dependent(&mut self, bc: TimeDependentBC) {
        self.time_dependent.push(bc);
    }

    /// Add a mixed BC.
    pub fn add_mixed(&mut self, bc: MixedBC) {
        self.mixed.push(bc);
    }

    /// Total number of boundary conditions.
    pub fn count(&self) -> usize {
        self.dirichlet.len()
            + self.neumann.len()
            + self.robin.len()
            + self.periodic.len()
            + self.symmetry.iter().map(|s| s.nodes.len()).sum::<usize>()
            + self.time_dependent.len()
            + self.mixed.len()
    }

    /// Apply all BCs to the dense stiffness matrix and force vector.
    ///
    /// Order: Neumann first (adds forces), then Robin (adds stiffness and
    /// forces), then Dirichlet (elimination, overrides other contributions).
    pub fn apply_all(&self, k: &mut Vec<f64>, f: &mut Vec<f64>, n_dofs: usize, penalty: f64) {
        // Mixed BCs (Neumann part)
        for bc in &self.mixed {
            let (_, neumann) = bc.to_separate_bcs();
            apply_neumann_scalar(f, n_dofs, &neumann);
        }

        apply_neumann_scalar(f, n_dofs, &self.neumann);
        apply_robin(k, f, n_dofs, &self.robin);

        // Periodic BCs
        if !self.periodic.is_empty() {
            apply_periodic_dense_penalty(k, f, n_dofs, &self.periodic, penalty);
        }

        // Symmetry BCs
        for sym in &self.symmetry {
            let sym_bcs = sym.to_dirichlet_bcs();
            apply_dirichlet_dense_penalty(k, f, n_dofs, &sym_bcs, penalty);
        }

        // Mixed BCs (Dirichlet part)
        for bc in &self.mixed {
            let (dirichlet, _) = bc.to_separate_bcs();
            apply_dirichlet_dense_penalty(k, f, n_dofs, &dirichlet, penalty);
        }

        apply_dirichlet_dense_penalty(k, f, n_dofs, &self.dirichlet, penalty);
    }

    /// Apply all BCs including time-dependent ones at a given time.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_all_at_time(
        &self,
        k: &mut Vec<f64>,
        f: &mut Vec<f64>,
        n_dofs: usize,
        penalty: f64,
        t: f64,
    ) {
        self.apply_all(k, f, n_dofs, penalty);

        // Time-dependent BCs override static Dirichlet
        if !self.time_dependent.is_empty() {
            let td_bcs: Vec<DirichletBC> = self
                .time_dependent
                .iter()
                .map(|bc| bc.to_dirichlet_bc(t))
                .collect();
            apply_dirichlet_dense_penalty(k, f, n_dofs, &td_bcs, penalty);
        }
    }
}

// ── Sparse CsrMatrix variants (original API) ──────────────────────────────────

/// Apply Dirichlet boundary conditions by row/column zeroing (CSR).
///
/// For each prescribed DOF `i` with value `u_i`:
/// - Off-diagonal entries in row and column `i` are zeroed
/// - `K[i][i]` is set to 1.0
/// - `f[i]` is set to `u_i`
/// - The load vector is adjusted: `f[j] -= K[j][i] * u_i` for all other DOFs `j`
///
/// This approach preserves symmetry and avoids ill-conditioning.
pub fn apply_dirichlet(k: &mut CsrMatrix, f: &mut [f64], bcs: &[DirichletBc]) {
    for bc in bcs {
        let dof = bc.global_dof();
        assert!(dof < k.nrows, "BC DOF {dof} out of range");

        if bc.value.abs() > 0.0 {
            for (row, f_row) in f.iter_mut().enumerate().take(k.nrows) {
                if row == dof {
                    continue;
                }
                let rs = k.row_ptr[row];
                let re = k.row_ptr[row + 1];
                for idx in rs..re {
                    if k.col_indices[idx] == dof {
                        *f_row -= k.values[idx] * bc.value;
                    }
                }
            }
        }

        let start = k.row_ptr[dof];
        let end = k.row_ptr[dof + 1];
        for idx in start..end {
            if k.col_indices[idx] == dof {
                k.values[idx] = 1.0;
            } else {
                k.values[idx] = 0.0;
            }
        }

        for row in 0..k.nrows {
            if row == dof {
                continue;
            }
            let rs = k.row_ptr[row];
            let re = k.row_ptr[row + 1];
            for idx in rs..re {
                if k.col_indices[idx] == dof {
                    k.values[idx] = 0.0;
                }
            }
        }

        f[dof] = bc.value;
    }
}

/// Apply Neumann boundary conditions by adding forces to the load vector (CSR).
pub fn apply_neumann(f: &mut [f64], bcs: &[NeumannBc]) {
    for bc in bcs {
        let base = bc.node * 3;
        f[base] += bc.force[0];
        f[base + 1] += bc.force[1];
        f[base + 2] += bc.force[2];
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── original sparse tests ─────────────────────────────────────────────

    #[test]
    fn test_dirichlet_bc_zeros_row_col() {
        let triplets: Vec<(usize, usize, f64)> = vec![
            (0, 0, 10.0),
            (0, 1, 1.0),
            (0, 2, 2.0),
            (1, 0, 1.0),
            (1, 1, 20.0),
            (1, 2, 3.0),
            (2, 0, 2.0),
            (2, 1, 3.0),
            (2, 2, 30.0),
        ];
        let mut k = CsrMatrix::from_triplets(3, 3, &triplets);
        let mut f = vec![100.0, 200.0, 300.0];

        let bcs = vec![DirichletBc {
            node: 0,
            dof: 1,
            value: 0.0,
        }];
        apply_dirichlet(&mut k, &mut f, &bcs);

        assert_eq!(k.get(1, 0), 0.0);
        assert!((k.get(1, 1) - 1.0).abs() < 1e-12);
        assert_eq!(k.get(1, 2), 0.0);
        assert_eq!(k.get(0, 1), 0.0);
        assert_eq!(k.get(2, 1), 0.0);
    }

    #[test]
    fn test_dirichlet_bc_applied() {
        let triplets: Vec<(usize, usize, f64)> =
            vec![(0, 0, 8.0), (0, 1, 3.0), (1, 0, 3.0), (1, 1, 5.0)];
        let mut k = CsrMatrix::from_triplets(2, 2, &triplets);
        let mut f = vec![16.0, 10.0];

        let bcs = vec![DirichletBc {
            node: 0,
            dof: 0,
            value: 0.0,
        }];
        apply_dirichlet(&mut k, &mut f, &bcs);

        assert!((k.get(0, 0) - 1.0).abs() < 1e-12, "diagonal should be 1");
        assert!(
            k.get(0, 1).abs() < 1e-12,
            "off-diagonal in row 0 should be 0"
        );
        assert!(
            k.get(1, 0).abs() < 1e-12,
            "off-diagonal in col 0 should be 0"
        );
        assert!(
            (f[0] - 0.0).abs() < 1e-12,
            "f[0] should equal prescribed value"
        );
    }

    #[test]
    fn test_neumann_bc() {
        let mut f = vec![0.0; 9]; // 3 nodes
        let bcs = vec![NeumannBc::new(2, [0.0, -1000.0, 0.0])];
        apply_neumann(&mut f, &bcs);
        assert_eq!(f[6], 0.0);
        assert_eq!(f[7], -1000.0);
        assert_eq!(f[8], 0.0);
    }

    // ── dense Dirichlet tests ─────────────────────────────────────────────

    #[test]
    fn dense_dirichlet_elimination_zeros_row_col() {
        let mut k = vec![4.0, 2.0, 2.0, 3.0];
        let mut f = vec![8.0, 6.0];
        let bcs = vec![DirichletBC::new(0, 0, 0.0)];
        apply_dirichlet_dense(&mut k, &mut f, 2, &bcs);
        assert!((k[0] - 1.0).abs() < 1e-12);
        assert!(k[1].abs() < 1e-12);
        assert!(k[2].abs() < 1e-12);
        assert!((f[0]).abs() < 1e-12);
    }

    #[test]
    fn dense_dirichlet_nonzero_value_adjusts_rhs() {
        let mut k = vec![4.0, 2.0, 2.0, 3.0];
        let mut f = vec![0.0, 0.0];
        let bcs = vec![DirichletBC::new(0, 0, 2.0)];
        apply_dirichlet_dense(&mut k, &mut f, 2, &bcs);
        assert!((f[1] - (-4.0)).abs() < 1e-12, "f[1]={}", f[1]);
        assert!((f[0] - 2.0).abs() < 1e-12, "f[0]={}", f[0]);
    }

    // ── Neumann scalar tests ──────────────────────────────────────────────

    #[test]
    fn neumann_scalar_adds_force() {
        let mut f = vec![0.0; 6]; // 2 nodes, 3 DOF each
        let bcs = vec![NeumannBC::new(1, 1, -500.0)]; // node 1, y-DOF
        apply_neumann_scalar(&mut f, 6, &bcs);
        assert_eq!(f[4], -500.0); // node 1, dof 1 → index 4
    }

    // ── Robin BC tests ────────────────────────────────────────────────────

    #[test]
    fn robin_bc_adds_stiffness_to_diagonal() {
        let mut k = vec![1.0, 0.0, 0.0, 1.0]; // 2×2 identity
        let mut f = vec![0.0, 0.0];
        let bcs = vec![RobinBC::new(0, 0, 100.0, 0.5)];
        apply_robin(&mut k, &mut f, 2, &bcs);
        assert!((k[0] - 101.0).abs() < 1e-12);
        assert!((f[0] - 50.0).abs() < 1e-12);
    }

    // ── BoundarySet combined apply test ───────────────────────────────────

    #[test]
    fn boundary_set_apply_all_combined() {
        let mut set = BoundarySet::new();
        set.add_neumann(NeumannBC::new(0, 2, 10.0)); // node 0, z-DOF += 10
        set.add_robin(RobinBC::new(0, 0, 50.0, 1.0)); // K[0,0] += 50, f[0] += 50
        set.add_dirichlet(DirichletBC::new(0, 1, 0.0)); // fix node 0, y-DOF

        let n = 6; // 2 nodes
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 1.0;
        }
        let mut f = vec![0.0; n];

        set.apply_all(&mut k, &mut f, n, 1e12);

        assert!((f[2] - 10.0).abs() < 1e-9, "f[2]={}", f[2]);
        assert!((k[0] - 51.0).abs() < 1e-9, "K[0,0]={}", k[0]);
        assert!(k[n + 1] > 1e11, "K[1,1]={}", k[n + 1]);
    }

    // ── Robin BC Improved tests ───────────────────────────────────────────

    #[test]
    fn robin_improved_applies_alpha_and_g() {
        let bc = RobinBCImproved::new(0, 0, 200.0, 1.0, 50.0);
        let n = 3;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 1.0;
        }
        let mut f = vec![0.0; n];
        bc.apply_dense(&mut k, &mut f, n, 3);
        assert!((k[0] - 201.0).abs() < 1e-12, "K[0,0]={}", k[0]);
        assert!((f[0] - 50.0).abs() < 1e-12, "f[0]={}", f[0]);
    }

    #[test]
    fn robin_improved_damping() {
        let bc = RobinBCImproved::new(0, 1, 100.0, 1.0, 0.0).with_damping(5.0);
        let mut c_diag = vec![0.0; 6];
        bc.apply_damping(&mut c_diag, 3);
        assert!((c_diag[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn apply_robin_improved_batch() {
        let bcs = vec![
            RobinBCImproved::new(0, 0, 10.0, 1.0, 5.0),
            RobinBCImproved::new(0, 1, 20.0, 1.0, 3.0),
        ];
        let n = 6;
        let mut k = vec![0.0; n * n];
        let mut f = vec![0.0; n];
        apply_robin_improved(&mut k, &mut f, n, &bcs, 3);
        assert!((k[0] - 10.0).abs() < 1e-12);
        assert!((k[n + 1] - 20.0).abs() < 1e-12);
        assert!((f[0] - 5.0).abs() < 1e-12);
        assert!((f[1] - 3.0).abs() < 1e-12);
    }

    // ── Periodic BC tests ────────────────────────────────────────────────

    #[test]
    fn periodic_bc_penalty_ties_dofs() {
        // 2 nodes (6 DOFs), tie node 0 x-DOF to node 1 x-DOF
        let n = 6;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 10.0;
        }
        let mut f = vec![0.0; n];
        let bcs = vec![PeriodicBC::new(0, 1, 0)]; // tie DOF 0 to DOF 3
        let penalty = 1e6;
        apply_periodic_dense_penalty(&mut k, &mut f, n, &bcs, penalty);

        // K[0,0] and K[3,3] should increase by penalty
        assert!((k[0] - (10.0 + penalty)).abs() < 1e-6);
        assert!((k[3 * n + 3] - (10.0 + penalty)).abs() < 1e-6);
        // K[0,3] and K[3,0] should decrease by penalty
        assert!((k[3] - (-penalty)).abs() < 1e-6);
        assert!((k[3 * n] - (-penalty)).abs() < 1e-6);
    }

    #[test]
    fn periodic_bc_with_offset() {
        let n = 6;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 10.0;
        }
        let mut f = vec![0.0; n];
        let bcs = vec![PeriodicBC::with_offset(0, 1, 0, 2.0)];
        let penalty = 1e6;
        apply_periodic_dense_penalty(&mut k, &mut f, n, &bcs, penalty);

        // f[slave] += penalty * offset = 1e6 * 2.0
        assert!((f[3] - penalty * 2.0).abs() < 1e-6);
        // f[master] -= penalty * offset
        assert!((f[0] - (-penalty * 2.0)).abs() < 1e-6);
    }

    #[test]
    fn periodic_bc_elimination() {
        let n = 6;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 10.0;
        }
        let mut f = vec![1.0; n];
        let bcs = vec![PeriodicBC::new(0, 1, 0)];
        apply_periodic_dense_elimination(&mut k, &mut f, n, &bcs);

        // Slave row: k[s,s] = 1, k[s,m] = -1, rest zero
        assert!((k[3 * n + 3] - 1.0).abs() < 1e-12);
        assert!((k[3 * n] - (-1.0)).abs() < 1e-12);
        assert!((f[3]).abs() < 1e-12); // offset = 0
    }

    // ── Symmetry BC tests ────────────────────────────────────────────────

    #[test]
    fn symmetry_bc_fixes_normal_dof() {
        let sym = SymmetryBC::new(vec![0, 1], 1); // fix y-DOF for nodes 0 and 1
        let bcs = sym.to_dirichlet_bcs();
        assert_eq!(bcs.len(), 2);
        assert_eq!(bcs[0].node_id, 0);
        assert_eq!(bcs[0].dof, 1);
        assert_eq!(bcs[0].value, 0.0);
        assert_eq!(bcs[1].node_id, 1);
    }

    #[test]
    fn symmetry_bc_apply_dense() {
        let sym = SymmetryBC::new(vec![0], 0); // fix x-DOF for node 0
        let n = 3;
        let mut k = vec![5.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 5.0];
        let mut f = vec![10.0, 20.0, 30.0];
        sym.apply_dense(&mut k, &mut f, n);
        // DOF 0 should be eliminated
        assert!((k[0] - 1.0).abs() < 1e-12);
        assert!(k[1].abs() < 1e-12);
        assert!(k[2].abs() < 1e-12);
        assert!((f[0]).abs() < 1e-12);
    }

    // ── Time-dependent BC tests ──────────────────────────────────────────

    #[test]
    fn time_function_constant() {
        let tf = TimeFunction::Constant(3.125);
        assert!((tf.evaluate(0.0) - 3.125).abs() < 1e-12);
        assert!((tf.evaluate(100.0) - 3.125).abs() < 1e-12);
    }

    #[test]
    fn time_function_ramp() {
        let tf = TimeFunction::Ramp { ramp_time: 2.0 };
        assert!((tf.evaluate(0.0)).abs() < 1e-12);
        assert!((tf.evaluate(1.0) - 0.5).abs() < 1e-12);
        assert!((tf.evaluate(2.0) - 1.0).abs() < 1e-12);
        assert!((tf.evaluate(5.0) - 1.0).abs() < 1e-12); // clamped
    }

    #[test]
    fn time_function_sinusoidal() {
        let tf = TimeFunction::Sinusoidal { frequency: 1.0 };
        assert!(tf.evaluate(0.0).abs() < 1e-12);
        assert!((tf.evaluate(0.25) - 1.0).abs() < 1e-12); // sin(π/2) = 1
    }

    #[test]
    fn time_function_step() {
        let tf = TimeFunction::Step { step_time: 1.0 };
        assert!((tf.evaluate(0.5)).abs() < 1e-12);
        assert!((tf.evaluate(1.0) - 1.0).abs() < 1e-12);
        assert!((tf.evaluate(2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn time_function_triangular() {
        let tf = TimeFunction::Triangular { period: 1.0 };
        // At t=0: phase=0, f = 4*0 - 1 = -1
        assert!((tf.evaluate(0.0) - (-1.0)).abs() < 1e-12);
        // At t=0.25: phase=0.25, f = 4*0.25 - 1 = 0
        assert!(tf.evaluate(0.25).abs() < 1e-12);
        // At t=0.5: phase=0.5, f = 3 - 4*0.5 = 1
        assert!((tf.evaluate(0.5) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn time_function_exponential_decay() {
        let tf = TimeFunction::ExponentialDecay { decay_rate: 1.0 };
        assert!((tf.evaluate(0.0) - 1.0).abs() < 1e-12);
        assert!((tf.evaluate(1.0) - (-1.0_f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn time_dependent_bc_value_at() {
        let bc = TimeDependentBC::new(0, 0, 5.0, TimeFunction::Ramp { ramp_time: 1.0 });
        assert!((bc.value_at(0.0)).abs() < 1e-12);
        assert!((bc.value_at(0.5) - 2.5).abs() < 1e-12);
        assert!((bc.value_at(1.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn time_dependent_bc_to_dirichlet() {
        let bc = TimeDependentBC::new(2, 1, 10.0, TimeFunction::Constant(1.0));
        let dbc = bc.to_dirichlet_bc(0.0);
        assert_eq!(dbc.node_id, 2);
        assert_eq!(dbc.dof, 1);
        assert!((dbc.value - 10.0).abs() < 1e-12);
    }

    #[test]
    fn apply_time_dependent_dense_test() {
        let bcs = vec![TimeDependentBC::new(
            0,
            0,
            1.0,
            TimeFunction::Step { step_time: 0.5 },
        )];
        let n = 3;
        let mut k = vec![5.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 5.0];
        let mut f = vec![10.0, 20.0, 30.0];

        // Before step time: value = 0
        apply_time_dependent_dense(&mut k, &mut f, n, &bcs, 0.0);
        assert!((f[0]).abs() < 1e-12);

        // Reset and test after step time
        let mut k2 = vec![5.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 5.0];
        let mut f2 = vec![10.0, 20.0, 30.0];
        apply_time_dependent_dense(&mut k2, &mut f2, n, &bcs, 1.0);
        assert!((f2[0] - 1.0).abs() < 1e-12);
    }

    // ── Mixed BC tests ──────────────────────────────────────────────────

    #[test]
    fn mixed_bc_roller() {
        let bc = MixedBC::roller(0, 1, 0.0); // fix y, free x and z
        let (dirichlet, neumann) = bc.to_separate_bcs();
        assert_eq!(dirichlet.len(), 1);
        assert_eq!(dirichlet[0].dof, 1);
        assert!(neumann.is_empty());
    }

    #[test]
    fn mixed_bc_pin_2d() {
        let bc = MixedBC::pin_2d(0, 2); // fix x and y, free z
        let (dirichlet, _) = bc.to_separate_bcs();
        assert_eq!(dirichlet.len(), 2);
        let dofs: Vec<usize> = dirichlet.iter().map(|d| d.dof).collect();
        assert!(dofs.contains(&0));
        assert!(dofs.contains(&1));
        assert!(!dofs.contains(&2));
    }

    #[test]
    fn mixed_bc_apply_dense_test() {
        let bc = MixedBC::new(
            0,
            vec![(0, 0.0)],    // fix x-DOF
            vec![(1, -100.0)], // apply force on y-DOF
        );
        let n = 3;
        let mut k = vec![5.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 5.0];
        let mut f = vec![0.0; n];
        bc.apply_dense(&mut k, &mut f, n);

        // y-DOF force applied
        assert!((f[1] - (-100.0)).abs() < 1e-12);
        // x-DOF eliminated
        assert!((k[0] - 1.0).abs() < 1e-12);
        assert!((f[0]).abs() < 1e-12);
    }

    // ── Sparse BC application tests ──────────────────────────────────────

    #[test]
    fn dirichlet_sparse_penalty_adds_to_diagonal() {
        let triplets: Vec<(usize, usize, f64)> =
            vec![(0, 0, 10.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 20.0)];
        let mut k = CsrMatrix::from_triplets(2, 2, &triplets);
        let mut f = vec![0.0; 2];
        let bcs = vec![DirichletBc::new(0, 0, 3.0)];
        let penalty = 1e10;
        apply_dirichlet_sparse_penalty(&mut k, &mut f, &bcs, penalty);

        assert!((k.get(0, 0) - (10.0 + penalty)).abs() < 1.0);
        assert!((f[0] - penalty * 3.0).abs() < 1.0);
        // Off-diagonals unchanged
        assert!((k.get(0, 1) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn robin_sparse_test() {
        let triplets: Vec<(usize, usize, f64)> = vec![(0, 0, 10.0), (1, 1, 20.0), (2, 2, 30.0)];
        let mut k = CsrMatrix::from_triplets(3, 3, &triplets);
        let mut f = vec![0.0; 3];
        let bcs = vec![RobinBC::new(0, 0, 50.0, 2.0)];
        apply_robin_sparse(&mut k, &mut f, &bcs, 3);
        assert!((k.get(0, 0) - 60.0).abs() < 1e-12);
        assert!((f[0] - 100.0).abs() < 1e-12);
    }

    #[test]
    fn neumann_sparse_adds_force() {
        let mut f = vec![0.0; 6];
        let bcs = vec![NeumannBC::new(1, 2, 42.0)];
        apply_neumann_sparse(&mut f, &bcs);
        assert!((f[5] - 42.0).abs() < 1e-12);
    }

    // ── BoundarySet count and time-dependent tests ──────────────────────

    #[test]
    fn boundary_set_count() {
        let mut set = BoundarySet::new();
        assert_eq!(set.count(), 0);
        set.add_dirichlet(DirichletBC::new(0, 0, 0.0));
        set.add_neumann(NeumannBC::new(0, 1, 1.0));
        set.add_robin(RobinBC::new(0, 2, 10.0, 0.0));
        set.add_periodic(PeriodicBC::new(0, 1, 0));
        set.add_symmetry(SymmetryBC::new(vec![0, 1, 2], 1));
        set.add_time_dependent(TimeDependentBC::new(0, 0, 1.0, TimeFunction::Constant(1.0)));
        set.add_mixed(MixedBC::roller(0, 0, 0.0));
        assert_eq!(set.count(), 9); // 1+1+1+1+3+1+1 (symmetry has 3 nodes)
    }

    #[test]
    fn boundary_set_apply_all_at_time() {
        let mut set = BoundarySet::new();
        set.add_time_dependent(TimeDependentBC::new(
            0,
            0,
            5.0,
            TimeFunction::Ramp { ramp_time: 1.0 },
        ));

        let n = 3;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 1.0;
        }
        let mut f = vec![0.0; n];

        set.apply_all_at_time(&mut k, &mut f, n, 1e12, 0.5);
        // At t=0.5, value = 5.0 * 0.5 = 2.5
        // Penalty method: f[0] += penalty * 2.5
        assert!(
            f[0] > 1e11,
            "f[0] should be large due to penalty, got {}",
            f[0]
        );
    }

    #[test]
    fn boundary_set_with_mixed_and_periodic() {
        let mut set = BoundarySet::new();
        set.add_mixed(MixedBC::new(0, vec![(0, 0.0)], vec![(1, 50.0)]));
        set.add_periodic(PeriodicBC::new(0, 1, 2)); // tie z-DOF

        let n = 6;
        let mut k = vec![0.0; n * n];
        for i in 0..n {
            k[i * n + i] = 1.0;
        }
        let mut f = vec![0.0; n];
        let penalty = 1e6;

        set.apply_all(&mut k, &mut f, n, penalty);

        // Neumann part of mixed: f[1] += 50
        assert!((f[1] - 50.0).abs() < 1e-6, "f[1]={}", f[1]);
    }
}
