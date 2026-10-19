// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM (Finite Element Method) Assembly and Solver.
//!
//! [`PyFemAssembly`] analyses a 3-D pin-jointed truss: nodes with coordinates,
//! two-node bars of a given axial stiffness `EA`, supports (fixed DOFs) and
//! nodal loads. The bar stiffness and axial force come from the engine's
//! [`TrussElement3D`](oxiphysics::fem::truss_frame::TrussElement3D); the
//! assembly scatters the `6 x 6` element matrices into the global system and
//! solves the free DOFs with a dense direct solver.

use oxiphysics::fem::truss_frame::TrussElement3D;
use pyo3::prelude::*;

// ===========================================================================
// FEM (Finite Element Method) Assembly and Solver
// ===========================================================================

/// Largest number of nodes of an assembly.
const MAX_NODES: usize = 1 << 20;

/// Largest number of free DOFs the dense direct solver accepts.
///
/// The reduced stiffness matrix is stored densely (`8 * n^2` bytes) and
/// factorised in `O(n^3)`; the limit keeps a mistyped node count from aborting
/// the interpreter. Bigger models belong to the sparse `PyFemSolver`.
const MAX_DENSE_DOFS: usize = 3000;

/// A DOF whose stiffness is below this fraction of the largest diagonal
/// stiffness counts as unstiffened.
const ZERO_STIFFNESS_FRACTION: f64 = 1.0e-12;

/// Relative tolerance of the declared bar length against the node distance.
const LENGTH_TOLERANCE: f64 = 1.0e-6;

/// Shortest bar the assembly accepts (m).
const MIN_BAR_LENGTH: f64 = 1.0e-12;

/// Reasons a [`PyFemAssembly`] rejects its input or cannot solve.
///
/// An out-of-range node or DOF surfaces in Python as an `IndexError`, every
/// other error as a `ValueError`.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FemError {
    /// A node index is not smaller than the number of nodes.
    #[error("node index {node} is out of range for an assembly with {node_count} nodes")]
    NodeOutOfRange {
        /// The offending node index.
        node: usize,
        /// Number of nodes of the assembly.
        node_count: usize,
    },
    /// A DOF index is not smaller than the number of DOFs.
    #[error("DOF index {dof} is out of range for an assembly with {n_dofs} DOFs")]
    DofOutOfRange {
        /// The offending DOF index.
        dof: usize,
        /// Number of DOFs of the assembly.
        n_dofs: usize,
    },
    /// More nodes were requested than an assembly can hold.
    #[error("{requested} nodes exceed the supported maximum of {limit}")]
    TooManyNodes {
        /// Requested node count.
        requested: usize,
        /// Largest supported node count.
        limit: usize,
    },
    /// A number that must be finite (and positive where stated) is not.
    #[error("{what} must be {requirement}, got {value}")]
    InvalidValue {
        /// What the value describes.
        what: &'static str,
        /// What the value has to satisfy.
        requirement: &'static str,
        /// The offending value.
        value: f64,
    },
    /// Both ends of a bar are the same node.
    #[error("bar {element} joins node {node} to itself")]
    BarToSameNode {
        /// Index of the bar (the index the bar would get).
        element: usize,
        /// The node used twice.
        node: usize,
    },
    /// A bar ends at a node whose coordinates were never set.
    #[error("node {node} has no coordinates (set them with set_node or from_coords)")]
    MissingNodeCoordinates {
        /// The node without coordinates.
        node: usize,
    },
    /// The two end nodes of a bar coincide, so it has no direction.
    #[error("bar {element} has zero length: its nodes {node_a} and {node_b} coincide")]
    DegenerateBar {
        /// Index of the bar.
        element: usize,
        /// First end node.
        node_a: usize,
        /// Second end node.
        node_b: usize,
    },
    /// The length given for a bar differs from the distance of its nodes.
    #[error("bar {element} declares the length {declared} m but its nodes are {geometric} m apart")]
    InconsistentLength {
        /// Index of the bar.
        element: usize,
        /// The length given when the bar was added.
        declared: f64,
        /// The distance of the end nodes.
        geometric: f64,
    },
    /// A loaded DOF is not stiffened by any bar.
    #[error("DOF {dof} carries a load of {load} N but no bar has stiffness along it")]
    LoadOnUnstiffenedDof {
        /// The loaded DOF.
        dof: usize,
        /// The load on it (N).
        load: f64,
    },
    /// The stiffness matrix of the free DOFs is singular.
    #[error(
        "the structure is unstable: its stiffness matrix is singular at DOF {dof} \
         (missing supports or a mechanism)"
    )]
    UnstableStructure {
        /// A DOF of the singular mode.
        dof: usize,
    },
    /// Too many free DOFs for the dense solver.
    #[error("{free_dofs} free DOFs exceed the dense solver limit of {limit}")]
    TooLargeForDenseSolver {
        /// Number of free DOFs.
        free_dofs: usize,
        /// Largest number the dense solver accepts.
        limit: usize,
    },
}

impl From<FemError> for PyErr {
    fn from(err: FemError) -> PyErr {
        let message = err.to_string();
        match err {
            FemError::NodeOutOfRange { .. } | FemError::DofOutOfRange { .. } => {
                pyo3::exceptions::PyIndexError::new_err(message)
            }
            _ => pyo3::exceptions::PyValueError::new_err(message),
        }
    }
}

/// A single 2-node bar (truss) element for FEM assembly.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct FemBarElement {
    /// Global indices of the two end nodes.
    #[pyo3(get, set)]
    pub nodes: [usize; 2],
    /// Young's modulus times cross-sectional area (EA).
    #[pyo3(get, set)]
    pub ea: f64,
    /// Length given when the bar was added, if any (m).
    ///
    /// The length of a bar is the distance of its end nodes; a declared length
    /// is only checked against it when the assembly is built.
    #[pyo3(get, set)]
    pub length: Option<f64>,
}

/// A bar of the assembled model: its DOFs and the engine's element.
#[derive(Debug, Clone)]
struct AssembledBar {
    /// Global DOFs `[u_x, u_y, u_z]` of the first node, then of the second.
    dofs: [usize; 6],
    /// Global `6 x 6` element stiffness matrix.
    stiffness: [[f64; 6]; 6],
    /// The engine element (axial force evaluation).
    element: TrussElement3D,
}

/// FEM assembly for a 3-D truss structure.
///
/// Each of the `n_nodes` nodes has three DOFs `u_x, u_y, u_z` (DOF
/// `3 * node + axis`). The workflow is: give the nodes coordinates
/// ([`from_coords`](Self::from_coords) or [`set_node`](Self::set_node)), add
/// bars between them, fix the supported DOFs, apply nodal loads and
/// [`solve`](Self::solve). Bars may point in any direction: the element
/// stiffness is `EA/L * [[c cᵀ, -c cᵀ], [-c cᵀ, c cᵀ]]` for the unit direction
/// `c` of the bar, so loads and stiffness act in all three directions.
///
/// The loads are never altered by `assemble` or `solve`; Dirichlet conditions
/// are imposed by eliminating the fixed DOFs.
#[pyclass(from_py_object)]
#[derive(Debug, Clone)]
pub struct PyFemAssembly {
    /// Number of degrees of freedom (nodes × 3 for 3-D).
    #[pyo3(get)]
    pub n_dofs: usize,
    /// Node coordinates; `None` until set.
    nodes: Vec<Option<[f64; 3]>>,
    /// Bar elements.
    elements: Vec<FemBarElement>,
    /// External force vector (the loads as the user applied them).
    forces: Vec<f64>,
    /// Fixed DOFs (Dirichlet BCs: displacement = 0).
    fixed: Vec<bool>,
    /// Validated bars, built by `assemble`.
    assembled_bars: Vec<AssembledBar>,
    /// Displacement solution vector.
    displacements: Vec<f64>,
    /// Whether `assembled_bars` matches the nodes and elements.
    assembled: bool,
    /// Whether `displacements` is the solution of the current model.
    solved: bool,
}

#[pymethods]
impl PyFemAssembly {
    /// Create a new FEM assembly with `n_nodes` 3-D nodes whose coordinates are
    /// still to be set.
    #[new]
    pub fn new(n_nodes: usize) -> Result<Self, FemError> {
        if n_nodes > MAX_NODES {
            return Err(FemError::TooManyNodes {
                requested: n_nodes,
                limit: MAX_NODES,
            });
        }
        let n_dofs = n_nodes * 3;
        Ok(Self {
            n_dofs,
            nodes: vec![None; n_nodes],
            elements: Vec::new(),
            forces: vec![0.0; n_dofs],
            fixed: vec![false; n_dofs],
            assembled_bars: Vec::new(),
            displacements: vec![0.0; n_dofs],
            assembled: false,
            solved: false,
        })
    }

    /// Create an assembly with one node per coordinate triple `[x, y, z]`.
    #[staticmethod]
    pub fn from_coords(coords: Vec<[f64; 3]>) -> Result<Self, FemError> {
        let mut assembly = Self::new(coords.len())?;
        for (node, position) in coords.into_iter().enumerate() {
            assembly.set_node(node, position)?;
        }
        Ok(assembly)
    }

    /// Set the coordinates `[x, y, z]` of `node`.
    pub fn set_node(&mut self, node: usize, coords: [f64; 3]) -> Result<(), FemError> {
        if let Some(&bad) = coords.iter().find(|c| !c.is_finite()) {
            return Err(FemError::InvalidValue {
                what: "node coordinate",
                requirement: "finite",
                value: bad,
            });
        }
        let node_count = self.nodes.len();
        let slot = self
            .nodes
            .get_mut(node)
            .ok_or(FemError::NodeOutOfRange { node, node_count })?;
        *slot = Some(coords);
        self.invalidate_model();
        Ok(())
    }

    /// Coordinates of `node`, or `None` if the node does not exist or has no
    /// coordinates yet.
    pub fn node(&self, node: usize) -> Option<[f64; 3]> {
        self.nodes.get(node).copied().flatten()
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of bar elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }

    /// Add a bar element between two nodes with the axial stiffness `ea`
    /// (Young's modulus times cross-section) and return its index.
    ///
    /// The bar spans the distance of its nodes. `length` may state that
    /// distance as well; it is only checked against the node coordinates when
    /// the assembly is built.
    #[pyo3(signature = (node_a, node_b, ea, length=None))]
    pub fn add_bar_element(
        &mut self,
        node_a: usize,
        node_b: usize,
        ea: f64,
        length: Option<f64>,
    ) -> Result<usize, FemError> {
        let node_count = self.nodes.len();
        for node in [node_a, node_b] {
            if node >= node_count {
                return Err(FemError::NodeOutOfRange { node, node_count });
            }
        }
        let element = self.elements.len();
        if node_a == node_b {
            return Err(FemError::BarToSameNode {
                element,
                node: node_a,
            });
        }
        if !(ea.is_finite() && ea > 0.0) {
            return Err(FemError::InvalidValue {
                what: "axial stiffness EA",
                requirement: "finite and positive",
                value: ea,
            });
        }
        if let Some(declared) = length
            && !(declared.is_finite() && declared > 0.0)
        {
            return Err(FemError::InvalidValue {
                what: "bar length",
                requirement: "finite and positive",
                value: declared,
            });
        }
        self.elements.push(FemBarElement {
            nodes: [node_a, node_b],
            ea,
            length,
        });
        self.invalidate_model();
        Ok(element)
    }

    /// Apply an external force to a DOF (node_index * 3 + direction); forces
    /// applied to the same DOF add up.
    pub fn apply_force(&mut self, dof: usize, force: f64) -> Result<(), FemError> {
        if !force.is_finite() {
            return Err(FemError::InvalidValue {
                what: "force",
                requirement: "finite",
                value: force,
            });
        }
        let n_dofs = self.n_dofs;
        let slot = self
            .forces
            .get_mut(dof)
            .ok_or(FemError::DofOutOfRange { dof, n_dofs })?;
        *slot += force;
        self.invalidate_results();
        Ok(())
    }

    /// The external force currently applied to `dof` (N), or `0.0` if the DOF
    /// does not exist.
    pub fn force(&self, dof: usize) -> f64 {
        self.forces.get(dof).copied().unwrap_or(0.0)
    }

    /// Remove all applied forces.
    pub fn clear_forces(&mut self) {
        self.forces.fill(0.0);
        self.invalidate_results();
    }

    /// Fix a DOF (zero displacement Dirichlet BC).
    pub fn fix_dof(&mut self, dof: usize) -> Result<(), FemError> {
        let n_dofs = self.n_dofs;
        let slot = self
            .fixed
            .get_mut(dof)
            .ok_or(FemError::DofOutOfRange { dof, n_dofs })?;
        *slot = true;
        self.invalidate_results();
        Ok(())
    }

    /// Build the bar elements of the model.
    ///
    /// Checks that every bar has two distinct end nodes with coordinates, that
    /// it has a direction (a non-zero length) and that a declared length agrees
    /// with the node distance, and computes the `6 x 6` global stiffness matrix
    /// of every bar (`EA/L * [[c cᵀ, -c cᵀ], [-c cᵀ, c cᵀ]]` for the unit
    /// direction `c`). `solve` calls this when the model has changed.
    pub fn assemble(&mut self) -> Result<(), FemError> {
        let mut bars = Vec::with_capacity(self.elements.len());
        for (index, element) in self.elements.iter().enumerate() {
            let [node_a, node_b] = element.nodes;
            let position = |node: usize| {
                self.nodes
                    .get(node)
                    .copied()
                    .flatten()
                    .ok_or(FemError::MissingNodeCoordinates { node })
            };
            let (start, end) = (position(node_a)?, position(node_b)?);
            // Only the product EA matters: it is carried by the modulus.
            let truss = TrussElement3D {
                elastic_modulus: element.ea,
                area: 1.0,
                node_i: start,
                node_j: end,
            };
            let length = truss.length();
            if length <= MIN_BAR_LENGTH {
                return Err(FemError::DegenerateBar {
                    element: index,
                    node_a,
                    node_b,
                });
            }
            if let Some(declared) = element.length
                && (declared - length).abs() > LENGTH_TOLERANCE * declared.max(length)
            {
                return Err(FemError::InconsistentLength {
                    element: index,
                    declared,
                    geometric: length,
                });
            }
            bars.push(AssembledBar {
                dofs: [
                    3 * node_a,
                    3 * node_a + 1,
                    3 * node_a + 2,
                    3 * node_b,
                    3 * node_b + 1,
                    3 * node_b + 2,
                ],
                stiffness: truss.global_stiffness(),
                element: truss,
            });
        }
        self.assembled_bars = bars;
        self.assembled = true;
        self.invalidate_results();
        Ok(())
    }

    /// Solve the linear system K·u = f for the free DOFs.
    ///
    /// Fixed DOFs have zero displacement; their support reactions are not part
    /// of the system, so loads applied to them are carried by the support. The
    /// system is solved by Gaussian elimination with partial pivoting; the
    /// applied loads are left untouched, so the model can be modified and
    /// solved again. Returns `True` on success, `False` for an assembly
    /// without nodes.
    ///
    /// # Errors
    ///
    /// `ValueError` if the model is invalid (see [`assemble`](Self::assemble)),
    /// if a loaded free DOF has no stiffness (for example a `y` load on a node
    /// whose bars all run along `x`), if the structure is unstable (singular
    /// stiffness matrix: missing supports or a mechanism) or if it has more
    /// than 3000 free DOFs.
    pub fn solve(&mut self) -> Result<bool, FemError> {
        if self.n_dofs == 0 {
            return Ok(false);
        }
        if !self.assembled {
            self.assemble()?;
        }
        self.solved = false;

        // Classify the DOFs: fixed, unstiffened (no bar acts along them) or free.
        let diagonal = self.diagonal_stiffness();
        let largest = diagonal.iter().copied().fold(0.0_f64, f64::max);
        let threshold = ZERO_STIFFNESS_FRACTION * largest;
        let mut free_index = vec![usize::MAX; self.n_dofs];
        let mut free_dofs: Vec<usize> = Vec::new();
        for (dof, &stiffness) in diagonal.iter().enumerate() {
            if self.fixed[dof] {
                continue;
            }
            if stiffness <= threshold {
                if self.forces[dof] != 0.0 {
                    return Err(FemError::LoadOnUnstiffenedDof {
                        dof,
                        load: self.forces[dof],
                    });
                }
                continue;
            }
            free_index[dof] = free_dofs.len();
            free_dofs.push(dof);
        }
        let n_free = free_dofs.len();
        if n_free > MAX_DENSE_DOFS {
            return Err(FemError::TooLargeForDenseSolver {
                free_dofs: n_free,
                limit: MAX_DENSE_DOFS,
            });
        }

        // Reduced system K_ff u_f = f_f, built into local storage.
        let mut matrix = vec![0.0; n_free * n_free];
        let mut rhs: Vec<f64> = free_dofs.iter().map(|&dof| self.forces[dof]).collect();
        for bar in &self.assembled_bars {
            for (row_local, &row_dof) in bar.dofs.iter().enumerate() {
                let row = free_index[row_dof];
                if row == usize::MAX {
                    continue;
                }
                for (col_local, &col_dof) in bar.dofs.iter().enumerate() {
                    let col = free_index[col_dof];
                    if col != usize::MAX {
                        matrix[row * n_free + col] += bar.stiffness[row_local][col_local];
                    }
                }
            }
        }
        if let Err(pivot) = solve_dense(&mut matrix, &mut rhs, n_free) {
            let dof = free_dofs.get(pivot).copied().unwrap_or(0);
            return Err(FemError::UnstableStructure { dof });
        }

        let mut displacements = vec![0.0; self.n_dofs];
        for (&dof, &value) in free_dofs.iter().zip(&rhs) {
            displacements[dof] = value;
        }
        self.displacements = displacements;
        self.solved = true;
        Ok(true)
    }

    /// Get the displacement at DOF `dof`, or `0.0` if out of bounds (or if the
    /// model has not been solved since it last changed).
    pub fn displacement(&self, dof: usize) -> f64 {
        self.displacements.get(dof).copied().unwrap_or(0.0)
    }

    /// Get all displacements as an owned `Vec<f64>`.
    pub fn all_displacements(&self) -> Vec<f64> {
        self.displacements.clone()
    }

    /// Compute the axial force in element `elem_idx`, positive in tension.
    ///
    /// It is `EA/L * c · (u_b - u_a)` for the unit direction `c` of the bar and
    /// the displacements of its end nodes. Returns `None` if the index is out
    /// of bounds or the model has not been solved since it last changed.
    pub fn element_force(&self, elem_idx: usize) -> Option<f64> {
        if !self.solved {
            return None;
        }
        let bar = self.assembled_bars.get(elem_idx)?;
        let mut end_displacements = [0.0; 6];
        for (slot, &dof) in end_displacements.iter_mut().zip(&bar.dofs) {
            *slot = self.displacements.get(dof).copied().unwrap_or(0.0);
        }
        Some(bar.element.axial_force(end_displacements))
    }
}

impl PyFemAssembly {
    /// The geometry or the bars changed: the assembled bars and the results are stale.
    fn invalidate_model(&mut self) {
        self.assembled = false;
        self.assembled_bars.clear();
        self.invalidate_results();
    }

    /// The loads, supports or model changed: the displacements are stale.
    fn invalidate_results(&mut self) {
        self.solved = false;
        self.displacements.fill(0.0);
    }

    /// Stiffness `K_ii` of every DOF (diagonal of the global matrix).
    fn diagonal_stiffness(&self) -> Vec<f64> {
        let mut diagonal = vec![0.0; self.n_dofs];
        for bar in &self.assembled_bars {
            for (local, &dof) in bar.dofs.iter().enumerate() {
                diagonal[dof] += bar.stiffness[local][local];
            }
        }
        diagonal
    }
}

/// Solve the dense system `a x = b` (`a` row-major `n x n`) in place by
/// Gaussian elimination with partial pivoting; the solution replaces `b`.
///
/// Returns the column at which the matrix turned out to be singular, i.e. its
/// pivot was not larger than `n * eps * max|a_ij|`.
fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> Result<(), usize> {
    let scale = a.iter().fold(0.0_f64, |largest, v| largest.max(v.abs()));
    let tolerance = n as f64 * f64::EPSILON * scale;
    for col in 0..n {
        let mut pivot_row = col;
        let mut best = a[col * n + col].abs();
        for row in (col + 1)..n {
            let candidate = a[row * n + col].abs();
            if candidate > best {
                best = candidate;
                pivot_row = row;
            }
        }
        if best.is_nan() || best <= tolerance {
            return Err(col);
        }
        if pivot_row != col {
            for j in 0..n {
                a.swap(col * n + j, pivot_row * n + j);
            }
            b.swap(col, pivot_row);
        }
        let pivot = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            if factor == 0.0 {
                continue;
            }
            for j in col..n {
                let delta = factor * a[col * n + j];
                a[row * n + j] -= delta;
            }
            let delta = factor * b[col];
            b[row] -= delta;
        }
    }
    for i in (0..n).rev() {
        let row = &a[i * n..(i + 1) * n];
        let tail: f64 = row[i + 1..]
            .iter()
            .zip(&b[i + 1..])
            .map(|(c, x)| c * x)
            .sum();
        b[i] = (b[i] - tail) / row[i];
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Axial stiffness of every bar of the tests (N).
    const EA: f64 = 1.0e6;

    fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{what}: {actual} vs {expected}"
        );
    }

    /// DOF of `node` along `axis` (0 = x, 1 = y, 2 = z).
    fn dof(node: usize, axis: usize) -> usize {
        3 * node + axis
    }

    /// Two bars at 45 degrees from the supports (-1, 0, 0) and (1, 0, 0) up to the
    /// apex (0, 1, 0); both supports are fixed in all directions.
    fn two_bar_truss() -> PyFemAssembly {
        let mut assembly =
            PyFemAssembly::from_coords(vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
                .expect("three nodes");
        assembly.add_bar_element(0, 2, EA, None).expect("left bar");
        assembly.add_bar_element(1, 2, EA, None).expect("right bar");
        for support_dof in 0..6 {
            assembly.fix_dof(support_dof).expect("support DOF");
        }
        assembly
    }

    #[test]
    fn test_symmetric_two_bar_truss_carries_an_apex_load_in_equal_compression() {
        let load = 1000.0;
        let mut assembly = two_bar_truss();
        assembly.apply_force(dof(2, 1), -load).expect("apex y load");
        assert!(assembly.solve().expect("a stable truss"));
        let member_force = -load / 2.0_f64.sqrt();
        for bar in 0..2 {
            let force = assembly.element_force(bar).expect("bar exists");
            assert_close(force, member_force, 1e-9, "member force");
        }
        // The apex sinks by P L / (EA sin^2 45 deg * 2) = P sqrt(2) / EA and does not sway.
        assert_close(
            assembly.displacement(dof(2, 1)),
            -load * 2.0_f64.sqrt() / EA,
            1e-12,
            "apex deflection",
        );
        assert_close(assembly.displacement(dof(2, 0)), 0.0, 1e-15, "apex sway");
        assert_close(assembly.displacement(dof(2, 2)), 0.0, 1e-15, "apex z");
        // The supports do not move.
        for support_dof in 0..6 {
            assert_eq!(assembly.displacement(support_dof), 0.0);
        }
    }

    #[test]
    fn test_a_sideways_load_goes_into_one_tension_and_one_compression_bar() {
        let mut assembly = two_bar_truss();
        assembly.apply_force(dof(2, 0), 500.0).expect("apex x load");
        assert!(assembly.solve().expect("a stable truss"));
        // Equilibrium of the apex: (N1 - N2) / sqrt(2) = P in x, (N1 + N2) / sqrt(2) = 0 in y.
        let (left, right) = (
            assembly.element_force(0).expect("left bar"),
            assembly.element_force(1).expect("right bar"),
        );
        assert_close(left, 500.0 / 2.0_f64.sqrt(), 1e-9, "left bar");
        assert_close(right, -500.0 / 2.0_f64.sqrt(), 1e-9, "right bar");
    }

    #[test]
    fn test_loads_survive_assembly_and_solving() {
        let mut assembly = two_bar_truss();
        assembly.apply_force(dof(2, 1), -1000.0).expect("y load");
        assembly.apply_force(dof(2, 0), 250.0).expect("x load");
        assembly.assemble().expect("a valid model");
        let first = {
            assembly.solve().expect("a stable truss");
            assembly.all_displacements()
        };
        // Neither assembling nor solving touches the applied loads ...
        assert_eq!(assembly.force(dof(2, 1)), -1000.0);
        assert_eq!(assembly.force(dof(2, 0)), 250.0);
        // ... so solving again, also after re-assembling, gives the same answer.
        assembly.solve().expect("a stable truss");
        assert_eq!(assembly.all_displacements(), first);
        assembly.assemble().expect("a valid model");
        assembly.solve().expect("a stable truss");
        assert_eq!(assembly.all_displacements(), first);
        assert_eq!(assembly.force(dof(2, 1)), -1000.0);
        // The model is linear: doubling the loads doubles the displacements.
        assembly.apply_force(dof(2, 1), -1000.0).expect("y load");
        assembly.apply_force(dof(2, 0), 250.0).expect("x load");
        assembly.solve().expect("a stable truss");
        for (doubled, single) in assembly.all_displacements().iter().zip(&first) {
            assert_close(*doubled, 2.0 * single, 1e-15, "superposition");
        }
        // Clearing the loads leaves the structure unloaded.
        assembly.clear_forces();
        assert_eq!(assembly.force(dof(2, 1)), 0.0);
        assembly.solve().expect("a stable truss");
        assert!(assembly.all_displacements().iter().all(|u| *u == 0.0));
    }

    #[test]
    fn test_loads_on_supported_dofs_are_carried_by_the_supports() {
        let mut assembly = two_bar_truss();
        assembly.apply_force(0, 12345.0).expect("load on a support");
        assembly.apply_force(dof(2, 1), -1000.0).expect("y load");
        assembly.solve().expect("a stable truss");
        assert_eq!(assembly.force(0), 12345.0);
        assert_close(
            assembly.element_force(0).expect("left bar"),
            -1000.0 / 2.0_f64.sqrt(),
            1e-9,
            "member force",
        );
    }

    #[test]
    fn test_tripod_under_a_general_load_satisfies_statics() {
        let apex = [0.3, 0.2, 1.5];
        let feet = [[1.0, 0.0, 0.0], [-0.6, 0.9, 0.0], [-0.4, -0.8, 0.1]];
        let stiffness = [1.0e6, 2.0e6, 1.5e6];
        let load = [200.0, -350.0, -900.0];
        let mut coords = vec![apex];
        coords.extend(feet);
        let mut assembly = PyFemAssembly::from_coords(coords).expect("four nodes");
        for (foot, ea) in stiffness.iter().enumerate() {
            assembly
                .add_bar_element(0, foot + 1, *ea, None)
                .expect("leg");
            for axis in 0..3 {
                assembly.fix_dof(dof(foot + 1, axis)).expect("foot DOF");
            }
        }
        for (axis, component) in load.iter().enumerate() {
            assembly
                .apply_force(dof(0, axis), *component)
                .expect("apex load");
        }
        assert!(assembly.solve().expect("a stable tripod"));
        // Unit vectors from the apex to the feet: the legs pull the apex along them.
        let directions: Vec<nalgebra::Vector3<f64>> = feet
            .iter()
            .map(|foot| {
                nalgebra::Vector3::new(foot[0] - apex[0], foot[1] - apex[1], foot[2] - apex[2])
                    .normalize()
            })
            .collect();
        let forces: Vec<f64> = (0..3)
            .map(|leg| assembly.element_force(leg).expect("leg force"))
            .collect();
        let resultant = directions
            .iter()
            .zip(&forces)
            .fold(nalgebra::Vector3::from(load), |sum, (direction, force)| {
                sum + direction * *force
            });
        assert!(
            resultant.norm() < 1e-8,
            "apex out of equilibrium: {resultant:?}"
        );
        // Statically determinate: the forces solve U N = -F for U = [u1 u2 u3].
        let directions_matrix =
            nalgebra::Matrix3::from_columns(&[directions[0], directions[1], directions[2]]);
        let expected = directions_matrix.try_inverse().expect("independent legs")
            * -nalgebra::Vector3::from(load);
        for (leg, force) in forces.iter().enumerate() {
            assert_close(*force, expected[leg], 1e-8, "leg force");
        }
        // The displacement of the apex follows from the leg forces: u . c_i = N_i L_i / EA_i.
        let displacement = nalgebra::Vector3::new(
            assembly.displacement(dof(0, 0)),
            assembly.displacement(dof(0, 1)),
            assembly.displacement(dof(0, 2)),
        );
        for (leg, foot) in feet.iter().enumerate() {
            let length = (nalgebra::Vector3::from(*foot) - nalgebra::Vector3::from(apex)).norm();
            // The leg runs from the foot (fixed) to the apex: its direction is -u_i.
            let extension = -directions[leg].dot(&displacement);
            assert_close(
                extension,
                forces[leg] * length / stiffness[leg],
                1e-12,
                "leg extension",
            );
        }
    }

    #[test]
    fn test_x_aligned_bar_matches_the_analytic_elongation() {
        let mut assembly =
            PyFemAssembly::from_coords(vec![[0.0; 3], [2.0, 0.0, 0.0]]).expect("two nodes");
        assembly.add_bar_element(0, 1, EA, Some(2.0)).expect("bar");
        assembly.fix_dof(0).expect("DOF 0");
        assembly.apply_force(dof(1, 0), 800.0).expect("axial load");
        assert!(assembly.solve().expect("a stable bar"));
        assert_close(
            assembly.displacement(dof(1, 0)),
            800.0 * 2.0 / EA,
            1e-15,
            "elongation",
        );
        assert_close(
            assembly.element_force(0).expect("bar"),
            800.0,
            1e-9,
            "axial force",
        );
        // The directions without stiffness stay at rest.
        assert_eq!(assembly.displacement(dof(1, 1)), 0.0);
        assert_eq!(assembly.displacement(dof(1, 2)), 0.0);
    }

    #[test]
    fn test_a_load_without_stiffness_is_rejected_instead_of_erased() {
        let mut assembly =
            PyFemAssembly::from_coords(vec![[0.0; 3], [1.0, 0.0, 0.0]]).expect("two nodes");
        assembly.add_bar_element(0, 1, EA, None).expect("bar");
        for axis in 0..3 {
            assembly.fix_dof(dof(0, axis)).expect("support DOF");
        }
        // A y load on a node whose only bar runs along x has nothing to resist it ...
        assembly.apply_force(dof(1, 1), 10.0).expect("y load");
        assert_eq!(
            assembly.solve(),
            Err(FemError::LoadOnUnstiffenedDof {
                dof: dof(1, 1),
                load: 10.0
            })
        );
        // ... and the load is still there, the error does not swallow it.
        assert_eq!(assembly.force(dof(1, 1)), 10.0);
        // Without the offending load the same structure solves.
        assembly.clear_forces();
        assembly.apply_force(dof(1, 0), 10.0).expect("x load");
        assert!(assembly.solve().expect("a stable bar"));
    }

    #[test]
    fn test_unstable_structures_are_reported() {
        // A single inclined bar: the free end can swing sideways.
        let mut assembly =
            PyFemAssembly::from_coords(vec![[0.0; 3], [1.0, 2.0, 2.0]]).expect("two nodes");
        assembly.add_bar_element(0, 1, EA, Some(3.0)).expect("bar");
        for axis in 0..3 {
            assembly.fix_dof(dof(0, axis)).expect("support DOF");
        }
        assembly.apply_force(dof(1, 0), 100.0).expect("load");
        assert!(matches!(
            assembly.solve(),
            Err(FemError::UnstableStructure { .. })
        ));
        // Without any support the structure is free to move as a whole.
        let mut free = PyFemAssembly::from_coords(vec![[0.0; 3], [1.0, 1.0, 1.0]]).expect("nodes");
        free.add_bar_element(0, 1, EA, None).expect("bar");
        assert!(matches!(
            free.solve(),
            Err(FemError::UnstableStructure { .. })
        ));
        // A failed solve leaves no results behind.
        assert!(assembly.element_force(0).is_none());
    }

    #[test]
    fn test_invalid_models_are_reported_when_assembling() {
        let mut assembly = PyFemAssembly::new(2).expect("two nodes");
        assembly.add_bar_element(0, 1, EA, None).expect("bar");
        assembly.fix_dof(0).expect("DOF 0");
        // No coordinates at all.
        assert_eq!(
            assembly.assemble(),
            Err(FemError::MissingNodeCoordinates { node: 0 })
        );
        assert_eq!(
            assembly.solve(),
            Err(FemError::MissingNodeCoordinates { node: 0 })
        );
        assembly.set_node(0, [1.0, 1.0, 1.0]).expect("node 0");
        assert_eq!(
            assembly.assemble(),
            Err(FemError::MissingNodeCoordinates { node: 1 })
        );
        // Coincident end nodes.
        assembly.set_node(1, [1.0, 1.0, 1.0]).expect("node 1");
        assert_eq!(
            assembly.assemble(),
            Err(FemError::DegenerateBar {
                element: 0,
                node_a: 0,
                node_b: 1
            })
        );
        // A declared length that disagrees with the geometry.
        assembly.set_node(1, [2.0, 1.0, 1.0]).expect("node 1");
        assembly
            .add_bar_element(0, 1, EA, Some(2.0))
            .expect("second bar");
        assert_eq!(
            assembly.assemble(),
            Err(FemError::InconsistentLength {
                element: 1,
                declared: 2.0,
                geometric: 1.0
            })
        );
    }

    #[test]
    fn test_input_validation() {
        let mut assembly = PyFemAssembly::new(2).expect("two nodes");
        assert_eq!(
            assembly.set_node(2, [0.0; 3]),
            Err(FemError::NodeOutOfRange {
                node: 2,
                node_count: 2
            })
        );
        assert!(matches!(
            assembly.set_node(0, [0.0, f64::NAN, 0.0]),
            Err(FemError::InvalidValue { .. })
        ));
        assert_eq!(
            assembly.add_bar_element(0, 5, EA, None),
            Err(FemError::NodeOutOfRange {
                node: 5,
                node_count: 2
            })
        );
        assert_eq!(
            assembly.add_bar_element(1, 1, EA, None),
            Err(FemError::BarToSameNode {
                element: 0,
                node: 1
            })
        );
        for ea in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                assembly.add_bar_element(0, 1, ea, None),
                Err(FemError::InvalidValue { .. })
            ));
        }
        for length in [0.0, -2.0, f64::NAN] {
            assert!(matches!(
                assembly.add_bar_element(0, 1, EA, Some(length)),
                Err(FemError::InvalidValue { .. })
            ));
        }
        assert_eq!(assembly.element_count(), 0);
        assert_eq!(
            assembly.apply_force(6, 1.0),
            Err(FemError::DofOutOfRange { dof: 6, n_dofs: 6 })
        );
        assert!(matches!(
            assembly.apply_force(0, f64::INFINITY),
            Err(FemError::InvalidValue { .. })
        ));
        assert_eq!(
            assembly.fix_dof(7),
            Err(FemError::DofOutOfRange { dof: 7, n_dofs: 6 })
        );
        assert_eq!(assembly.force(99), 0.0);
        assert_eq!(assembly.displacement(99), 0.0);
        assert!(assembly.element_force(0).is_none());
        assert_eq!(assembly.node(5), None);
        assert_eq!(assembly.node(0), None);
        // An assembly without nodes has nothing to solve.
        let mut empty = PyFemAssembly::new(0).expect("no nodes");
        assert_eq!(empty.solve(), Ok(false));
    }

    #[test]
    fn test_unreasonable_sizes_are_rejected_without_allocating() {
        for nodes in [MAX_NODES + 1, usize::MAX] {
            assert_eq!(
                PyFemAssembly::new(nodes).err(),
                Some(FemError::TooManyNodes {
                    requested: nodes,
                    limit: MAX_NODES
                })
            );
        }
        // More free DOFs than the dense solver takes: a spiral of bars along which
        // every node is stiffened in all three directions.
        let coords: Vec<[f64; 3]> = (0..1100)
            .map(|i| {
                let t = i as f64;
                [t, (0.7 * t).sin(), (1.3 * t).cos()]
            })
            .collect();
        let mut assembly = PyFemAssembly::from_coords(coords).expect("1100 nodes");
        for node in 0..1098 {
            assembly
                .add_bar_element(node, node + 1, EA, None)
                .expect("bar");
            assembly
                .add_bar_element(node, node + 2, EA, None)
                .expect("skip bar");
        }
        assembly
            .add_bar_element(1098, 1099, EA, None)
            .expect("last bar");
        for axis in 0..3 {
            assembly.fix_dof(dof(0, axis)).expect("support DOF");
        }
        assert!(matches!(
            assembly.solve(),
            Err(FemError::TooLargeForDenseSolver {
                limit: MAX_DENSE_DOFS,
                ..
            })
        ));
    }

    #[test]
    fn test_results_are_dropped_when_the_model_changes() {
        let mut assembly = two_bar_truss();
        assembly.apply_force(dof(2, 1), -1000.0).expect("load");
        assert!(assembly.element_force(0).is_none(), "nothing solved yet");
        assembly.solve().expect("a stable truss");
        assert!(assembly.element_force(0).is_some());
        assert!(assembly.displacement(dof(2, 1)) < 0.0);
        // Any change of loads, supports, geometry or bars makes the results stale.
        assembly.apply_force(dof(2, 1), -1.0).expect("load");
        assert!(assembly.element_force(0).is_none());
        assert_eq!(assembly.displacement(dof(2, 1)), 0.0);
        assembly.solve().expect("a stable truss");
        assembly.set_node(2, [0.0, 2.0, 0.0]).expect("node");
        assert!(assembly.element_force(0).is_none());
        assembly.solve().expect("a stable truss");
        assert!(assembly.element_force(0).is_some());
        assembly.add_bar_element(0, 1, EA, None).expect("tie");
        assert!(assembly.element_force(0).is_none());
        assembly.solve().expect("a stable truss");
        // The new tie between the supports carries nothing; its index is 2.
        assert_close(
            assembly.element_force(2).expect("tie"),
            0.0,
            1e-9,
            "tie force",
        );
        assert!(assembly.element_force(3).is_none());
    }

    #[test]
    fn test_taller_apex_carries_more_axial_force_per_load() {
        // The member force of the symmetric truss is P / (2 sin(angle)).
        for height in [0.5, 1.0, 2.0, 5.0] {
            let mut assembly = two_bar_truss();
            assembly.set_node(2, [0.0, height, 0.0]).expect("apex");
            assembly.apply_force(dof(2, 1), -1000.0).expect("load");
            assembly.solve().expect("a stable truss");
            let sine = height / (1.0 + height * height).sqrt();
            assert_close(
                assembly.element_force(0).expect("bar"),
                -1000.0 / (2.0 * sine),
                1e-8,
                "member force",
            );
        }
    }

    #[test]
    fn test_solve_dense_handles_row_swaps_and_reports_singular_matrices() {
        // [[0, 2, 1], [1, 1, 1], [2, 1, 3]] x = (7, 6, 13) has the solution (1, 2, 3).
        let mut a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0];
        let mut b = vec![7.0, 6.0, 13.0];
        assert_eq!(solve_dense(&mut a, &mut b, 3), Ok(()));
        for (value, expected) in b.iter().zip([1.0, 2.0, 3.0]) {
            assert_close(*value, expected, 1e-12, "solution");
        }
        // Proportional rows make the matrix singular.
        let mut singular = vec![1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 0.0, 1.0];
        let mut rhs = vec![1.0, 2.0, 3.0];
        assert_eq!(solve_dense(&mut singular, &mut rhs, 3), Err(2));
        // The empty system is solved trivially; a NaN entry is reported.
        assert_eq!(solve_dense(&mut [], &mut [], 0), Ok(()));
        let mut nan = vec![f64::NAN];
        let mut rhs = vec![1.0];
        assert_eq!(solve_dense(&mut nan, &mut rhs, 1), Err(0));
    }
}
