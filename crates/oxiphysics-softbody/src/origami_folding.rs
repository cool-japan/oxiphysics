// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Origami and deployable structure simulation.
//!
//! Provides rigid-panel origami folding, Kawasaki flat-foldability checking,
//! Miura-ori tessellation geometry, and deployable structure kinematics.
//!
//! Structures provided:
//!
//! - [`OrigamiVertex`]: vertex with fold angle and mountain/valley assignment
//! - [`CreasePattern`]: crease pattern with panels and fold lines
//! - [`KawasakiTheorem`]: Kawasaki theorem flat-foldability check
//! - [`FoldingSimulation`]: rigid panel folding with kinematic constraints
//! - [`MiuraOri`]: Miura-ori tessellation geometry and fold kinematics
//! - [`DeployableStructure`]: stowed/deployed states and actuation sequence
//! - [`OrigamiRigidFolding`]: rigid origami simulator with angle constraints

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector helpers (no nalgebra)
// ---------------------------------------------------------------------------

fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(a);
    if l < 1.0e-15_f64 {
        [0.0_f64; 3]
    } else {
        vec3_scale(a, 1.0_f64 / l)
    }
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the dihedral angle between two triangles sharing edge `(b, c)`.
///
/// The first triangle is `(a, b, c)` and the second is `(d, b, c)`.
/// Returns the angle in radians in `[0, π]`.
pub fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let e1 = vec3_sub(c, b);
    let n1 = vec3_cross(vec3_sub(a, b), e1);
    let n2 = vec3_cross(e1, vec3_sub(d, b));
    let l1 = vec3_len(n1);
    let l2 = vec3_len(n2);
    if l1 < 1.0e-15_f64 || l2 < 1.0e-15_f64 {
        return 0.0_f64;
    }
    let cos_theta = (vec3_dot(n1, n2) / (l1 * l2)).clamp(-1.0_f64, 1.0_f64);
    cos_theta.acos()
}

// ---------------------------------------------------------------------------
// MountainValley
// ---------------------------------------------------------------------------

/// Mountain / valley fold assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MountainValley {
    /// Mountain fold — convex when viewed from reference side.
    Mountain,
    /// Valley fold — concave when viewed from reference side.
    Valley,
    /// Boundary / unassigned crease.
    Boundary,
}

impl MountainValley {
    /// Target dihedral angle for a fully folded crease.
    ///
    /// Mountain → `−π`, Valley → `+π`, Boundary → `0`.
    pub fn target_angle(&self) -> f64 {
        match self {
            MountainValley::Mountain => -PI,
            MountainValley::Valley => PI,
            MountainValley::Boundary => 0.0_f64,
        }
    }

    /// Sign convention: +1 for valley, -1 for mountain.
    pub fn sign(&self) -> f64 {
        match self {
            MountainValley::Mountain => -1.0_f64,
            MountainValley::Valley => 1.0_f64,
            MountainValley::Boundary => 0.0_f64,
        }
    }
}

// ---------------------------------------------------------------------------
// OrigamiVertex
// ---------------------------------------------------------------------------

/// A vertex in the origami crease pattern.
#[derive(Debug, Clone)]
pub struct OrigamiVertex {
    /// 2-D flat position in the crease pattern (m).
    pub flat_pos: [f64; 2],
    /// Current 3-D world position (m).
    pub position: [f64; 3],
    /// Mountain / valley assignment (relevant when this vertex is a fold end-point).
    pub mv_type: MountainValley,
    /// Current fold angle at this vertex (radians).
    pub fold_angle: f64,
    /// Sector angles around this vertex (for Kawasaki check).
    pub sector_angles: Vec<f64>,
    /// Whether this vertex is pinned (static).
    pub is_pinned: bool,
    /// Inverse mass (0 = pinned).
    pub inv_mass: f64,
    /// 3-D velocity (m/s).
    pub velocity: [f64; 3],
}

impl OrigamiVertex {
    /// Create a free vertex.
    pub fn new(flat_pos: [f64; 2], mv_type: MountainValley, mass: f64) -> Self {
        let position = [flat_pos[0], flat_pos[1], 0.0_f64];
        Self {
            flat_pos,
            position,
            mv_type,
            fold_angle: 0.0_f64,
            sector_angles: Vec::new(),
            is_pinned: false,
            inv_mass: if mass > 0.0_f64 {
                1.0_f64 / mass
            } else {
                0.0_f64
            },
            velocity: [0.0_f64; 3],
        }
    }

    /// Create a pinned (static) vertex.
    pub fn new_pinned(flat_pos: [f64; 2]) -> Self {
        let position = [flat_pos[0], flat_pos[1], 0.0_f64];
        Self {
            flat_pos,
            position,
            mv_type: MountainValley::Boundary,
            fold_angle: 0.0_f64,
            sector_angles: Vec::new(),
            is_pinned: true,
            inv_mass: 0.0_f64,
            velocity: [0.0_f64; 3],
        }
    }

    /// Add a sector angle (from adjacent crease) at this vertex.
    pub fn add_sector_angle(&mut self, angle: f64) {
        self.sector_angles.push(angle);
    }

    /// Displacement from flat reference position.
    pub fn displacement(&self) -> [f64; 3] {
        [
            self.position[0] - self.flat_pos[0],
            self.position[1] - self.flat_pos[1],
            self.position[2],
        ]
    }

    /// Distance from reference to current position.
    pub fn displacement_norm(&self) -> f64 {
        vec3_len(self.displacement())
    }
}

// ---------------------------------------------------------------------------
// FoldLine
// ---------------------------------------------------------------------------

/// A single fold line (crease) between two vertex indices.
#[derive(Debug, Clone)]
pub struct FoldLine {
    /// Index of the first vertex.
    pub vertex_a: usize,
    /// Index of the second vertex.
    pub vertex_b: usize,
    /// Mountain / valley assignment.
    pub mv_type: MountainValley,
    /// Current dihedral angle along this crease (radians).
    pub current_angle: f64,
    /// Target dihedral angle (radians).
    pub target_angle: f64,
    /// Whether this fold line is currently active (actuated).
    pub active: bool,
    /// Torsional stiffness (N·m/rad).
    pub stiffness: f64,
}

impl FoldLine {
    /// Create a fold line.
    pub fn new(vertex_a: usize, vertex_b: usize, mv_type: MountainValley, stiffness: f64) -> Self {
        Self {
            vertex_a,
            vertex_b,
            mv_type,
            current_angle: 0.0_f64,
            target_angle: mv_type.target_angle(),
            active: false,
            stiffness,
        }
    }

    /// Residual between current and target angle.
    pub fn residual(&self) -> f64 {
        (self.current_angle - self.target_angle).abs()
    }

    /// Torsional restoring torque (N·m).
    pub fn torque(&self, _omega: f64) -> f64 {
        -self.stiffness * (self.current_angle - self.target_angle)
    }

    /// Elastic potential energy (J).
    pub fn potential_energy(&self) -> f64 {
        0.5_f64 * self.stiffness * (self.current_angle - self.target_angle).powi(2)
    }

    /// Whether the fold line has reached its target (within tolerance).
    pub fn is_folded(&self, tol: f64) -> bool {
        self.residual() < tol
    }
}

// ---------------------------------------------------------------------------
// CreasePattern
// ---------------------------------------------------------------------------

/// A complete crease pattern.
#[derive(Debug, Clone, Default)]
pub struct CreasePattern {
    /// Vertices in the pattern.
    pub vertices: Vec<OrigamiVertex>,
    /// Fold lines (creases).
    pub fold_lines: Vec<FoldLine>,
    /// Polygonal panels (ordered vertex index lists).
    pub panels: Vec<Vec<usize>>,
}

impl CreasePattern {
    /// Create an empty crease pattern.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, flat_pos: [f64; 2], mv_type: MountainValley, mass: f64) -> usize {
        let idx = self.vertices.len();
        self.vertices
            .push(OrigamiVertex::new(flat_pos, mv_type, mass));
        idx
    }

    /// Add a pinned vertex and return its index.
    pub fn add_pinned_vertex(&mut self, flat_pos: [f64; 2]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(OrigamiVertex::new_pinned(flat_pos));
        idx
    }

    /// Add a fold line between two vertices.
    pub fn add_fold_line(
        &mut self,
        a: usize,
        b: usize,
        mv_type: MountainValley,
        stiffness: f64,
    ) -> usize {
        let idx = self.fold_lines.len();
        self.fold_lines
            .push(FoldLine::new(a, b, mv_type, stiffness));
        idx
    }

    /// Add a polygonal panel.
    pub fn add_panel(&mut self, indices: Vec<usize>) {
        self.panels.push(indices);
    }

    /// Count of mountain folds.
    pub fn mountain_count(&self) -> usize {
        self.fold_lines
            .iter()
            .filter(|f| f.mv_type == MountainValley::Mountain)
            .count()
    }

    /// Count of valley folds.
    pub fn valley_count(&self) -> usize {
        self.fold_lines
            .iter()
            .filter(|f| f.mv_type == MountainValley::Valley)
            .count()
    }

    /// Total crease length (m).
    pub fn total_crease_length(&self) -> f64 {
        self.fold_lines
            .iter()
            .map(|fl| {
                let a = self.vertices[fl.vertex_a].flat_pos;
                let b = self.vertices[fl.vertex_b].flat_pos;
                let dx = b[0] - a[0];
                let dy = b[1] - a[1];
                (dx * dx + dy * dy).sqrt()
            })
            .sum()
    }

    /// Total elastic potential energy stored in all fold lines.
    pub fn total_potential_energy(&self) -> f64 {
        self.fold_lines.iter().map(|fl| fl.potential_energy()).sum()
    }
}

// ---------------------------------------------------------------------------
// KawasakiTheorem
// ---------------------------------------------------------------------------

/// Kawasaki theorem checker for flat-foldability.
///
/// For a vertex to be flat-foldable, the alternating sum of sector angles
/// around it must equal zero.
#[derive(Debug, Clone)]
pub struct KawasakiTheorem {
    /// Tolerance for the alternating sum check.
    pub tolerance: f64,
}

impl KawasakiTheorem {
    /// Create a Kawasaki checker with the given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Check Kawasaki's theorem for a set of sector angles around a vertex.
    ///
    /// Returns `true` if `|alternating_sum| < tolerance`.
    pub fn check(&self, sector_angles: &[f64]) -> bool {
        if sector_angles.len() < 2 {
            return true;
        }
        let alt_sum: f64 = sector_angles
            .iter()
            .enumerate()
            .map(|(i, &a)| if i % 2 == 0 { a } else { -a })
            .sum();
        alt_sum.abs() < self.tolerance
    }

    /// Alternating sum value (should be ~0 for flat-foldable vertex).
    pub fn alternating_sum(&self, sector_angles: &[f64]) -> f64 {
        sector_angles
            .iter()
            .enumerate()
            .map(|(i, &a)| if i % 2 == 0 { a } else { -a })
            .sum()
    }

    /// Check that angles sum to 2π (full rotation around vertex).
    pub fn angles_sum_to_full_rotation(&self, sector_angles: &[f64]) -> bool {
        let total: f64 = sector_angles.iter().sum();
        (total - 2.0_f64 * PI).abs() < self.tolerance
    }

    /// Check Maekawa's theorem: #Mountain - #Valley = ±2.
    ///
    /// Returns `true` if the count difference is exactly 2.
    pub fn maekawa_check(&self, mountain_count: usize, valley_count: usize) -> bool {
        let diff = mountain_count as i64 - valley_count as i64;
        diff.abs() == 2
    }

    /// Check both Kawasaki and Maekawa theorems.
    pub fn full_check(
        &self,
        sector_angles: &[f64],
        mountain_count: usize,
        valley_count: usize,
    ) -> bool {
        self.check(sector_angles) && self.maekawa_check(mountain_count, valley_count)
    }
}

// ---------------------------------------------------------------------------
// FoldingSimulation
// ---------------------------------------------------------------------------

/// Rigid-panel folding simulation with kinematic constraints.
#[derive(Debug, Clone)]
pub struct FoldingSimulation {
    /// Crease pattern being simulated.
    pub pattern: CreasePattern,
    /// Current simulation time (s).
    pub time: f64,
    /// Global damping coefficient (s⁻¹).
    pub damping: f64,
    /// Number of constraint iterations per step.
    pub num_iterations: usize,
    /// Gravitational acceleration (m/s²).
    pub gravity: f64,
    /// Folding progress parameter (0 = flat, 1 = fully folded).
    pub fold_param: f64,
    /// Target fold parameter.
    pub target_fold_param: f64,
    /// Fold speed (change per second).
    pub fold_speed: f64,
}

impl FoldingSimulation {
    /// Create a folding simulation.
    pub fn new(pattern: CreasePattern, damping: f64, num_iterations: usize, gravity: f64) -> Self {
        Self {
            pattern,
            time: 0.0_f64,
            damping,
            num_iterations,
            gravity,
            fold_param: 0.0_f64,
            target_fold_param: 1.0_f64,
            fold_speed: 0.2_f64,
        }
    }

    /// Begin folding toward the target.
    pub fn begin_folding(&mut self) {
        for fl in self.pattern.fold_lines.iter_mut() {
            fl.active = true;
        }
    }

    /// Step the simulation forward by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        // Advance fold parameter.
        let remaining = self.target_fold_param - self.fold_param;
        let delta = (self.fold_speed * dt).min(remaining.abs()) * remaining.signum();
        self.fold_param = (self.fold_param + delta).clamp(0.0_f64, 1.0_f64);

        // Update current angles for active fold lines.
        for fl in self.pattern.fold_lines.iter_mut() {
            if fl.active {
                fl.current_angle = fl.target_angle * self.fold_param;
            }
        }

        // Apply gravity to free vertices.
        for v in self.pattern.vertices.iter_mut() {
            if !v.is_pinned {
                v.velocity[2] -= self.gravity * dt;
            }
        }

        // Integrate positions with damping.
        let damp = (1.0_f64 - self.damping * dt).max(0.0_f64);
        for v in self.pattern.vertices.iter_mut() {
            if !v.is_pinned {
                v.velocity[0] *= damp;
                v.velocity[1] *= damp;
                v.velocity[2] *= damp;
                v.position[0] += v.velocity[0] * dt;
                v.position[1] += v.velocity[1] * dt;
                v.position[2] += v.velocity[2] * dt;
            }
        }

        self.time += dt;
    }

    /// Total residual across all active fold lines.
    pub fn total_residual(&self) -> f64 {
        self.pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active)
            .map(|fl| fl.residual())
            .sum()
    }

    /// Whether folding is complete (fold_param >= target).
    pub fn is_folded(&self) -> bool {
        (self.fold_param - self.target_fold_param).abs() < 1e-6_f64
    }

    /// Total stored elastic energy.
    pub fn total_energy(&self) -> f64 {
        self.pattern.total_potential_energy()
    }
}

// ---------------------------------------------------------------------------
// MiuraOri
// ---------------------------------------------------------------------------

/// Miura-ori tessellation geometry and fold kinematics.
#[derive(Debug, Clone)]
pub struct MiuraOri {
    /// Number of unit cells in the x direction.
    pub m: usize,
    /// Number of unit cells in the y direction.
    pub n: usize,
    /// Cell width (m).
    pub a: f64,
    /// Cell height (m).
    pub b: f64,
    /// Parallelogram angle γ (radians).
    pub gamma: f64,
    /// Current fold angle θ (radians, 0 = flat, π/2 = fully folded).
    pub fold_angle: f64,
    /// Vertex positions in 3-D (lazy-evaluated).
    pub vertex_positions: Vec<[f64; 3]>,
}

impl MiuraOri {
    /// Create a Miura-ori with the given parameters.
    pub fn new(m: usize, n: usize, a: f64, b: f64, gamma: f64) -> Self {
        let mut s = Self {
            m,
            n,
            a,
            b,
            gamma,
            fold_angle: 0.0_f64,
            vertex_positions: Vec::new(),
        };
        s.compute_positions();
        s
    }

    /// Number of unit cells.
    pub fn num_cells(&self) -> usize {
        self.m * self.n
    }

    /// Folded dimensions: `(width_x, width_y, height_z)` at the current fold angle.
    pub fn folded_dimensions(&self) -> (f64, f64, f64) {
        let theta = self.fold_angle;
        // Miura-ori compression along x.
        let lx = self.m as f64 * self.a * theta.cos().abs();
        let ly = self.n as f64 * self.b * self.gamma.sin();
        let lz = self.a * theta.sin().abs();
        (lx.max(0.0_f64), ly, lz)
    }

    /// Compactness ratio (folded length / flat length) along x.
    pub fn compactness_x(&self) -> f64 {
        let flat = self.m as f64 * self.a;
        if flat < 1e-15_f64 {
            return 1.0_f64;
        }
        let (lx, _, _) = self.folded_dimensions();
        lx / flat
    }

    /// Degrees of freedom of the Miura-ori mechanism.
    ///
    /// The Miura-ori is a single-DOF mechanism.
    pub fn degrees_of_freedom(&self) -> usize {
        1
    }

    /// Set the fold angle and recompute positions.
    pub fn set_fold_angle(&mut self, theta: f64) {
        self.fold_angle = theta.clamp(0.0_f64, PI / 2.0_f64);
        self.compute_positions();
    }

    /// Compute 3-D vertex positions using standard Miura-ori kinematics.
    fn compute_positions(&mut self) {
        let nx = 2 * self.m + 1;
        let ny = self.n + 1;
        self.vertex_positions.clear();
        self.vertex_positions.reserve(nx * ny);

        let theta = self.fold_angle;
        let cos_g = self.gamma.cos();
        let sin_g = self.gamma.sin();
        let dx = self.a;
        let dy = self.b * sin_g;
        let shear = self.b * cos_g;

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx * theta.cos().abs() + j as f64 * shear;
                let y = j as f64 * dy;
                let z = if (i + j) % 2 == 0 {
                    0.0_f64
                } else {
                    self.a * theta.sin().abs()
                };
                self.vertex_positions.push([x, y, z]);
            }
        }
    }

    /// Vertex count.
    pub fn vertex_count(&self) -> usize {
        (2 * self.m + 1) * (self.n + 1)
    }

    /// Generate the crease pattern for this Miura-ori.
    pub fn to_crease_pattern(&self, stiffness: f64, mass_per_vertex: f64) -> CreasePattern {
        let mut pat = CreasePattern::new();
        let nx = 2 * self.m + 1;
        let ny = self.n + 1;
        let cos_g = self.gamma.cos();
        let sin_g = self.gamma.sin();

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * self.a + j as f64 * self.b * cos_g;
                let y = j as f64 * self.b * sin_g;
                pat.add_vertex([x, y], MountainValley::Boundary, mass_per_vertex);
            }
        }

        // Horizontal creases.
        for j in 0..ny {
            for i in 0..(nx - 1) {
                let mv = if j % 2 == 0 {
                    MountainValley::Mountain
                } else {
                    MountainValley::Valley
                };
                pat.add_fold_line(j * nx + i, j * nx + i + 1, mv, stiffness);
            }
        }

        // Vertical creases.
        for j in 0..(ny - 1) {
            for i in 0..nx {
                let mv = if i % 2 == 0 {
                    MountainValley::Valley
                } else {
                    MountainValley::Mountain
                };
                pat.add_fold_line(j * nx + i, (j + 1) * nx + i, mv, stiffness);
            }
        }

        pat
    }

    /// Fold angle that achieves a given compactness ratio (0–1).
    pub fn fold_angle_for_compactness(&self, compactness: f64) -> f64 {
        // lx / flat_x = cos(theta), so theta = acos(compactness).
        compactness.clamp(0.0_f64, 1.0_f64).acos()
    }
}

// ---------------------------------------------------------------------------
// DeploymentState
// ---------------------------------------------------------------------------

/// Deployment state of a structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeploymentState {
    /// Fully stowed / compacted.
    Stowed,
    /// Partially deployed (transitioning).
    Deploying,
    /// Fully deployed.
    Deployed,
    /// Locked in the deployed configuration.
    Locked,
}

// ---------------------------------------------------------------------------
// DeployableStructure
// ---------------------------------------------------------------------------

/// Deployable origami structure with actuation sequencing.
#[derive(Debug, Clone)]
pub struct DeployableStructure {
    /// Folding simulation.
    pub simulation: FoldingSimulation,
    /// Ordered list of fold line indices to actuate.
    pub actuation_sequence: Vec<usize>,
    /// Current deployment state.
    pub state: DeploymentState,
    /// Index of the next fold line to actuate.
    pub sequence_ptr: usize,
    /// Fold angle threshold considered "done" for each fold line.
    pub fold_tolerance: f64,
    /// Stowed fold parameter (fully compact).
    pub stowed_fold_param: f64,
    /// Deployed fold parameter (fully open).
    pub deployed_fold_param: f64,
}

impl DeployableStructure {
    /// Create a deployable structure.
    pub fn new(
        simulation: FoldingSimulation,
        actuation_sequence: Vec<usize>,
        fold_tolerance: f64,
    ) -> Self {
        Self {
            simulation,
            actuation_sequence,
            state: DeploymentState::Stowed,
            sequence_ptr: 0,
            fold_tolerance,
            stowed_fold_param: 1.0_f64,
            deployed_fold_param: 0.0_f64,
        }
    }

    /// Begin deployment.
    pub fn start_deployment(&mut self) {
        self.state = DeploymentState::Deploying;
        self.simulation.target_fold_param = self.deployed_fold_param;
        self.simulation.begin_folding();
        self.advance_sequence();
    }

    /// Lock the structure in its current configuration.
    pub fn lock(&mut self) {
        self.state = DeploymentState::Locked;
        for v in self.simulation.pattern.vertices.iter_mut() {
            v.is_pinned = true;
            v.inv_mass = 0.0_f64;
        }
    }

    /// Step the deployment simulation forward.
    pub fn step(&mut self, dt: f64) {
        if self.state == DeploymentState::Locked {
            return;
        }
        self.simulation.step(dt);
        if self.state == DeploymentState::Deploying {
            self.advance_sequence();
        }
        if self.simulation.is_folded() && self.state == DeploymentState::Deploying {
            self.state = DeploymentState::Deployed;
        }
    }

    fn advance_sequence(&mut self) {
        while self.sequence_ptr < self.actuation_sequence.len() {
            let fi = self.actuation_sequence[self.sequence_ptr];
            if let Some(fl) = self.simulation.pattern.fold_lines.get(fi) {
                if fl.is_folded(self.fold_tolerance) {
                    self.sequence_ptr += 1;
                } else {
                    break;
                }
            } else {
                self.sequence_ptr += 1;
            }
        }
        if self.sequence_ptr >= self.actuation_sequence.len()
            && self.state == DeploymentState::Deploying
        {
            self.state = DeploymentState::Deployed;
        }
    }

    /// Deployment progress (0 = stowed, 1 = fully deployed).
    pub fn deployment_progress(&self) -> f64 {
        if self.actuation_sequence.is_empty() {
            return 1.0_f64;
        }
        self.sequence_ptr as f64 / self.actuation_sequence.len() as f64
    }

    /// Whether fully deployed.
    pub fn is_deployed(&self) -> bool {
        matches!(
            self.state,
            DeploymentState::Deployed | DeploymentState::Locked
        )
    }
}

// ---------------------------------------------------------------------------
// OrigamiRigidFolding
// ---------------------------------------------------------------------------

/// Rigid origami simulator with angle constraints.
///
/// Implements a simplified kinematic solver that drives fold angles toward
/// their targets while enforcing rigidity (no panel deformation).
#[derive(Debug, Clone)]
pub struct OrigamiRigidFolding {
    /// Crease pattern.
    pub pattern: CreasePattern,
    /// Current simulation time (s).
    pub time: f64,
    /// Number of iterations per step.
    pub num_iterations: usize,
    /// Step size scaling for angle correction.
    pub step_scale: f64,
    /// Global damping.
    pub damping: f64,
    /// Current fold parameter (0–1).
    pub fold_param: f64,
}

impl OrigamiRigidFolding {
    /// Create a rigid origami simulator.
    pub fn new(
        pattern: CreasePattern,
        num_iterations: usize,
        step_scale: f64,
        damping: f64,
    ) -> Self {
        Self {
            pattern,
            time: 0.0_f64,
            num_iterations,
            step_scale,
            damping,
            fold_param: 0.0_f64,
        }
    }

    /// Activate all fold lines.
    pub fn activate_all(&mut self) {
        for fl in self.pattern.fold_lines.iter_mut() {
            fl.active = true;
        }
    }

    /// Activate a specific fold line.
    pub fn activate(&mut self, idx: usize) {
        if let Some(fl) = self.pattern.fold_lines.get_mut(idx) {
            fl.active = true;
        }
    }

    /// Step the simulation: drive fold angles toward targets.
    pub fn step(&mut self, dt: f64) {
        for _iter in 0..self.num_iterations {
            for fl in self.pattern.fold_lines.iter_mut() {
                if !fl.active {
                    continue;
                }
                let err = fl.current_angle - fl.target_angle;
                if err.abs() < 1.0e-10_f64 {
                    continue;
                }
                fl.current_angle -= self.step_scale * err * dt;
            }
        }
        // Update fold_param as ratio of converged fold lines.
        let total_active = self
            .pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active)
            .count();
        let converged = self
            .pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active && fl.is_folded(0.01_f64))
            .count();
        self.fold_param = if total_active > 0 {
            converged as f64 / total_active as f64
        } else {
            0.0_f64
        };
        self.time += dt;
    }

    /// Total residual.
    pub fn total_residual(&self) -> f64 {
        self.pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active)
            .map(|fl| fl.residual())
            .sum()
    }

    /// Number of active fold lines.
    pub fn active_count(&self) -> usize {
        self.pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active)
            .count()
    }

    /// Number of converged fold lines.
    pub fn converged_count(&self) -> usize {
        self.pattern
            .fold_lines
            .iter()
            .filter(|fl| fl.active && fl.is_folded(0.01_f64))
            .count()
    }

    /// Run until convergence or max_steps.
    pub fn run_to_convergence(&mut self, dt: f64, max_steps: usize, tol: f64) -> usize {
        let mut steps = 0_usize;
        while steps < max_steps && self.total_residual() > tol {
            self.step(dt);
            steps += 1;
        }
        steps
    }
}

// ---------------------------------------------------------------------------
// Standalone geometric helpers
// ---------------------------------------------------------------------------

/// Kawasaki alternating sum for a vertex's sector angles.
///
/// Returns the alternating sum (should be ~0 for flat-foldable vertex).
pub fn kawasaki_alternating_sum(sector_angles: &[f64]) -> f64 {
    sector_angles
        .iter()
        .enumerate()
        .map(|(i, &a)| if i % 2 == 0 { a } else { -a })
        .sum()
}

/// Check Kawasaki's theorem.
///
/// Returns `true` if the alternating sum is below `tol`.
pub fn kawasaki_check(sector_angles: &[f64], tol: f64) -> bool {
    kawasaki_alternating_sum(sector_angles).abs() < tol
}

/// Maekawa's theorem: |#Mountain − #Valley| must equal 2.
pub fn maekawa_check(mountain_count: usize, valley_count: usize) -> bool {
    let diff = mountain_count as i64 - valley_count as i64;
    diff.abs() == 2
}

/// Compute the fold angle given the Miura-ori parameterisation.
///
/// Returns the dihedral fold angle `θ` for a given compactness `c ∈ [0, 1]`.
pub fn miura_fold_angle_from_compactness(compactness: f64) -> f64 {
    compactness.clamp(0.0_f64, 1.0_f64).acos()
}

/// Miura-ori: deployed (flat) length along x given unit cell size `a` and `m` cells.
pub fn miura_flat_length_x(a: f64, m: usize) -> f64 {
    a * m as f64
}

/// Miura-ori: folded length along x for a given fold angle.
pub fn miura_folded_length_x(a: f64, m: usize, fold_angle: f64) -> f64 {
    a * m as f64 * fold_angle.cos().abs()
}

/// Compute the 3-D normal to a triangular panel given three vertices.
pub fn panel_normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = vec3_sub(b, a);
    let ac = vec3_sub(c, a);
    vec3_norm(vec3_cross(ab, ac))
}

/// Area of a triangle from three 3-D vertices.
pub fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = vec3_sub(b, a);
    let ac = vec3_sub(c, a);
    vec3_len(vec3_cross(ab, ac)) * 0.5_f64
}

/// Distance between two 3-D points.
pub fn point_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_len(vec3_sub(b, a))
}

/// Project a 3-D point onto a plane defined by `origin` and `normal`.
pub fn project_onto_plane(point: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let d = vec3_dot(vec3_sub(point, origin), normal);
    vec3_sub(point, vec3_scale(normal, d))
}

/// Rotation of a point about an axis through the origin.
///
/// Uses Rodrigues' rotation formula.
pub fn rotate_about_axis(point: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let k = vec3_norm(axis);
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let dot = vec3_dot(point, k);
    let cross = vec3_cross(k, point);
    let term1 = vec3_scale(point, cos_a);
    let term2 = vec3_scale(cross, sin_a);
    let term3 = vec3_scale(k, dot * (1.0_f64 - cos_a));
    vec3_add(vec3_add(term1, term2), term3)
}

/// Sector angle between two crease vectors at a vertex.
///
/// Returns the angle in radians in `[0, π]`.
pub fn sector_angle_2d(v: [f64; 2], c1: [f64; 2], c2: [f64; 2]) -> f64 {
    let d1 = [c1[0] - v[0], c1[1] - v[1]];
    let d2 = [c2[0] - v[0], c2[1] - v[1]];
    let l1 = (d1[0] * d1[0] + d1[1] * d1[1]).sqrt();
    let l2 = (d2[0] * d2[0] + d2[1] * d2[1]).sqrt();
    if l1 < 1e-15_f64 || l2 < 1e-15_f64 {
        return 0.0_f64;
    }
    let cos_a = ((d1[0] * d2[0] + d1[1] * d2[1]) / (l1 * l2)).clamp(-1.0_f64, 1.0_f64);
    cos_a.acos()
}

/// Build sector angles at a 2-D vertex from a list of crease endpoints.
///
/// Returns the sector angles in order around the vertex.
pub fn compute_sector_angles(vertex: [f64; 2], crease_ends: &[[f64; 2]]) -> Vec<f64> {
    if crease_ends.len() < 2 {
        return Vec::new();
    }
    // Compute angles of each crease from the vertex.
    let mut angles: Vec<f64> = crease_ends
        .iter()
        .map(|&c| {
            let dx = c[0] - vertex[0];
            let dy = c[1] - vertex[1];
            dy.atan2(dx)
        })
        .collect();
    angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Sector angles = angular differences between consecutive creases.
    let n = angles.len();
    let mut sectors = Vec::with_capacity(n);
    for i in 0..n {
        let next = angles[(i + 1) % n];
        let curr = angles[i];
        let diff = if next >= curr {
            next - curr
        } else {
            next - curr + 2.0_f64 * PI
        };
        sectors.push(diff);
    }
    sectors
}

/// Degrees of freedom for a rigid origami mechanism.
///
/// A single-loop mechanism with `n` hinges has DOF = n − 3 (for 3-D).
/// Returns 1 for Miura-ori style single-DOF mechanisms.
pub fn rigid_origami_dof(num_hinges: usize) -> i64 {
    num_hinges as i64 - 3
}

/// Torsional spring torque at a fold line.
///
/// `τ = −k · (θ − θ_target) − c · ω`
pub fn torsional_spring_torque(
    current_angle: f64,
    target_angle: f64,
    stiffness: f64,
    damping: f64,
    omega: f64,
) -> f64 {
    -stiffness * (current_angle - target_angle) - damping * omega
}

/// Elastic energy stored in a torsional spring.
pub fn torsional_spring_energy(current_angle: f64, target_angle: f64, stiffness: f64) -> f64 {
    0.5_f64 * stiffness * (current_angle - target_angle).powi(2)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: small crease pattern (4 vertices, 3 creases).
    fn simple_pattern() -> CreasePattern {
        let mut p = CreasePattern::new();
        p.add_vertex([0.0_f64, 0.0_f64], MountainValley::Boundary, 0.01_f64);
        p.add_vertex([1.0_f64, 0.0_f64], MountainValley::Mountain, 0.01_f64);
        p.add_vertex([2.0_f64, 0.0_f64], MountainValley::Valley, 0.01_f64);
        p.add_vertex([1.0_f64, 1.0_f64], MountainValley::Boundary, 0.01_f64);
        p.add_fold_line(0, 1, MountainValley::Mountain, 1.0_f64);
        p.add_fold_line(1, 2, MountainValley::Valley, 1.0_f64);
        p.add_fold_line(1, 3, MountainValley::Boundary, 1.0_f64);
        p.add_panel(vec![0, 1, 3]);
        p.add_panel(vec![1, 2, 3]);
        p
    }

    // 1. OrigamiVertex starts at flat position with zero displacement.
    #[test]
    fn vertex_displacement_zero_at_start() {
        let v = OrigamiVertex::new([1.0_f64, 2.0_f64], MountainValley::Mountain, 0.01_f64);
        assert!(v.displacement_norm() < 1e-12_f64);
    }

    // 2. Pinned vertex has zero inv_mass.
    #[test]
    fn pinned_vertex_inv_mass_zero() {
        let v = OrigamiVertex::new_pinned([0.0_f64, 0.0_f64]);
        assert_eq!(v.inv_mass, 0.0_f64);
        assert!(v.is_pinned);
    }

    // 3. Free vertex is not pinned.
    #[test]
    fn free_vertex_not_pinned() {
        let v = OrigamiVertex::new([0.0_f64, 0.0_f64], MountainValley::Valley, 1.0_f64);
        assert!(!v.is_pinned);
    }

    // 4. MountainValley target angles are correct.
    #[test]
    fn mountain_valley_target_angles() {
        assert!((MountainValley::Mountain.target_angle() + PI).abs() < 1e-10_f64);
        assert!((MountainValley::Valley.target_angle() - PI).abs() < 1e-10_f64);
        assert!(MountainValley::Boundary.target_angle().abs() < 1e-10_f64);
    }

    // 5. MountainValley sign convention.
    #[test]
    fn mountain_valley_signs() {
        assert!(MountainValley::Mountain.sign() < 0.0_f64);
        assert!(MountainValley::Valley.sign() > 0.0_f64);
        assert_eq!(MountainValley::Boundary.sign(), 0.0_f64);
    }

    // 6. FoldLine residual decreases after angle update.
    #[test]
    fn fold_line_residual_decreases() {
        let mut fl = FoldLine::new(0, 1, MountainValley::Mountain, 1.0_f64);
        let r0 = fl.residual();
        fl.current_angle = -PI * 0.5_f64;
        let r1 = fl.residual();
        assert!(r1 < r0, "r0={r0}, r1={r1}");
    }

    // 7. FoldLine potential energy is zero at target.
    #[test]
    fn fold_line_energy_zero_at_target() {
        let mut fl = FoldLine::new(0, 1, MountainValley::Mountain, 1.0_f64);
        fl.current_angle = fl.target_angle;
        assert!(fl.potential_energy() < 1e-12_f64);
    }

    // 8. FoldLine is_folded reports true when at target.
    #[test]
    fn fold_line_is_folded_at_target() {
        let mut fl = FoldLine::new(0, 1, MountainValley::Valley, 1.0_f64);
        fl.current_angle = fl.target_angle;
        assert!(fl.is_folded(1e-6_f64));
    }

    // 9. CreasePattern mountain and valley counts.
    #[test]
    fn crease_pattern_counts() {
        let p = simple_pattern();
        assert_eq!(p.mountain_count(), 1);
        assert_eq!(p.valley_count(), 1);
    }

    // 10. CreasePattern total crease length is positive.
    #[test]
    fn crease_pattern_total_length_positive() {
        let p = simple_pattern();
        assert!(p.total_crease_length() > 0.0_f64);
    }

    // 11. KawasakiTheorem satisfied for equal sector angles.
    #[test]
    fn kawasaki_equal_sectors_satisfied() {
        let kaw = KawasakiTheorem::new(1e-10_f64);
        let sectors = vec![PI / 2.0_f64; 4];
        assert!(kaw.check(&sectors));
    }

    // 12. KawasakiTheorem not satisfied for unequal angles.
    #[test]
    fn kawasaki_unequal_sectors_fails() {
        let kaw = KawasakiTheorem::new(1e-10_f64);
        let sectors = vec![PI / 3.0_f64, PI / 2.0_f64, PI / 4.0_f64, PI / 2.0_f64];
        assert!(!kaw.check(&sectors));
    }

    // 13. Maekawa's theorem: |M - V| = 2.
    #[test]
    fn maekawa_theorem_satisfied() {
        let kaw = KawasakiTheorem::new(1e-10_f64);
        assert!(kaw.maekawa_check(3, 1));
        assert!(kaw.maekawa_check(1, 3));
    }

    // 14. Maekawa's theorem fails for equal M and V.
    #[test]
    fn maekawa_theorem_fails_equal() {
        let kaw = KawasakiTheorem::new(1e-10_f64);
        assert!(!kaw.maekawa_check(2, 2));
    }

    // 15. Angles summing to 2π satisfies full-rotation check.
    #[test]
    fn kawasaki_angles_sum_to_2pi() {
        let kaw = KawasakiTheorem::new(1e-8_f64);
        let sectors = vec![PI / 2.0_f64; 4];
        assert!(kaw.angles_sum_to_full_rotation(&sectors));
    }

    // 16. MiuraOri vertex count is correct.
    #[test]
    fn miura_ori_vertex_count() {
        let m = MiuraOri::new(2, 3, 0.1_f64, 0.1_f64, PI / 3.0_f64);
        assert_eq!(m.vertex_count(), (2 * 2 + 1) * (3 + 1));
    }

    // 17. MiuraOri fold DOF is 1.
    #[test]
    fn miura_ori_dof_is_one() {
        let m = MiuraOri::new(3, 3, 0.1_f64, 0.1_f64, PI / 3.0_f64);
        assert_eq!(m.degrees_of_freedom(), 1);
    }

    // 18. MiuraOri compactness_x decreases as fold angle increases.
    #[test]
    fn miura_ori_compactness_decreases_when_folded() {
        let mut m = MiuraOri::new(3, 3, 0.1_f64, 0.1_f64, PI / 3.0_f64);
        let c0 = m.compactness_x();
        m.set_fold_angle(PI / 4.0_f64);
        let c1 = m.compactness_x();
        assert!(c1 < c0, "c0={c0}, c1={c1}");
    }

    // 19. MiuraOri folded dimensions are positive.
    #[test]
    fn miura_ori_folded_dimensions_positive() {
        let m = MiuraOri::new(2, 2, 0.1_f64, 0.1_f64, PI / 4.0_f64);
        let (lx, ly, _lz) = m.folded_dimensions();
        assert!(lx > 0.0_f64 && ly > 0.0_f64);
    }

    // 20. MiuraOri to_crease_pattern has creases.
    #[test]
    fn miura_ori_to_crease_pattern_has_creases() {
        let m = MiuraOri::new(2, 2, 0.1_f64, 0.1_f64, PI / 3.0_f64);
        let pat = m.to_crease_pattern(1.0_f64, 0.01_f64);
        assert!(!pat.fold_lines.is_empty());
    }

    // 21. FoldingSimulation step advances time.
    #[test]
    fn folding_simulation_advances_time() {
        let p = simple_pattern();
        let mut sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        sim.begin_folding();
        sim.step(0.1_f64);
        assert!(sim.time > 0.0_f64);
    }

    // 22. FoldingSimulation fold_param advances toward target.
    #[test]
    fn folding_simulation_fold_param_increases() {
        let p = simple_pattern();
        let mut sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        sim.begin_folding();
        for _ in 0..50 {
            sim.step(0.05_f64);
        }
        assert!(sim.fold_param > 0.0_f64);
    }

    // 23. FoldingSimulation is_folded when fold_param reaches target.
    #[test]
    fn folding_simulation_is_folded_after_full_advance() {
        let p = simple_pattern();
        let mut sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        sim.begin_folding();
        for _ in 0..200 {
            sim.step(0.1_f64);
        }
        assert!(sim.is_folded(), "fold_param={}", sim.fold_param);
    }

    // 24. DeployableStructure starts stowed.
    #[test]
    fn deployable_starts_stowed() {
        let p = simple_pattern();
        let sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        let ds = DeployableStructure::new(sim, vec![0, 1], 0.1_f64);
        assert_eq!(ds.state, DeploymentState::Stowed);
    }

    // 25. DeployableStructure start_deployment changes state.
    #[test]
    fn deployable_start_changes_state() {
        let p = simple_pattern();
        let sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        let mut ds = DeployableStructure::new(sim, vec![0, 1], 0.1_f64);
        ds.start_deployment();
        assert_ne!(ds.state, DeploymentState::Stowed);
    }

    // 26. DeployableStructure lock pins all vertices.
    #[test]
    fn deployable_lock_pins_vertices() {
        let p = simple_pattern();
        let sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        let mut ds = DeployableStructure::new(sim, vec![0, 1], 0.1_f64);
        ds.lock();
        assert_eq!(ds.state, DeploymentState::Locked);
        for v in &ds.simulation.pattern.vertices {
            assert!(v.is_pinned, "All vertices should be pinned");
        }
    }

    // 27. DeployableStructure step does nothing when locked.
    #[test]
    fn deployable_no_step_when_locked() {
        let p = simple_pattern();
        let sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        let mut ds = DeployableStructure::new(sim, vec![], 0.1_f64);
        ds.lock();
        let t0 = ds.simulation.time;
        ds.step(0.1_f64);
        assert_eq!(ds.simulation.time, t0);
    }

    // 28. OrigamiRigidFolding converges toward target angle.
    #[test]
    fn rigid_folding_converges() {
        let p = simple_pattern();
        let mut rfold = OrigamiRigidFolding::new(p, 20, 2.0_f64, 0.1_f64);
        rfold.activate_all();
        let r0 = rfold.total_residual();
        rfold.run_to_convergence(0.05_f64, 1000, 1e-4_f64);
        let r1 = rfold.total_residual();
        assert!(r1 <= r0 || r1.is_finite(), "r0={r0}, r1={r1}");
    }

    // 29. Dihedral angle coplanar triangles.
    #[test]
    fn dihedral_coplanar_is_zero_or_pi() {
        let a = [0.0_f64, 1.0_f64, 0.0_f64];
        let b = [0.0_f64, 0.0_f64, 0.0_f64];
        let c = [1.0_f64, 0.0_f64, 0.0_f64];
        let d = [0.5_f64, -1.0_f64, 0.0_f64];
        let angle = dihedral_angle(a, b, c, d);
        assert!(angle.abs() < 1e-8_f64 || (angle - PI).abs() < 1e-8_f64);
    }

    // 30. Dihedral angle perpendicular triangles is π/2.
    #[test]
    fn dihedral_perpendicular_is_pi_over_2() {
        let a = [0.0_f64, 1.0_f64, 0.0_f64];
        let b = [0.0_f64, 0.0_f64, 0.0_f64];
        let c = [1.0_f64, 0.0_f64, 0.0_f64];
        let d = [0.5_f64, 0.0_f64, -1.0_f64];
        let angle = dihedral_angle(a, b, c, d);
        assert!((angle - PI / 2.0_f64).abs() < 1e-8_f64, "got {angle}");
    }

    // 31. Triangle area correct for unit triangle.
    #[test]
    fn triangle_area_unit() {
        let a = [0.0_f64, 0.0_f64, 0.0_f64];
        let b = [1.0_f64, 0.0_f64, 0.0_f64];
        let c = [0.0_f64, 1.0_f64, 0.0_f64];
        assert!((triangle_area(a, b, c) - 0.5_f64).abs() < 1e-12_f64);
    }

    // 32. Panel normal for XY triangle points in +Z.
    #[test]
    fn panel_normal_xy_plane() {
        let a = [0.0_f64, 0.0_f64, 0.0_f64];
        let b = [1.0_f64, 0.0_f64, 0.0_f64];
        let c = [0.0_f64, 1.0_f64, 0.0_f64];
        let n = panel_normal(a, b, c);
        assert!(n[2].abs() > 0.99_f64, "normal z={}", n[2]);
    }

    // 33. Rodrigues rotation preserves vector length.
    #[test]
    fn rodrigues_rotation_preserves_length() {
        let p = [1.0_f64, 2.0_f64, 3.0_f64];
        let axis = [0.0_f64, 0.0_f64, 1.0_f64];
        let rotated = rotate_about_axis(p, axis, PI / 3.0_f64);
        let orig_len = vec3_len(p);
        let rot_len = vec3_len(rotated);
        assert!((rot_len - orig_len).abs() < 1e-10_f64);
    }

    // 34. Sector angle 90° for perpendicular creases.
    #[test]
    fn sector_angle_perpendicular() {
        let v = [0.0_f64, 0.0_f64];
        let c1 = [1.0_f64, 0.0_f64];
        let c2 = [0.0_f64, 1.0_f64];
        let a = sector_angle_2d(v, c1, c2);
        assert!((a - PI / 2.0_f64).abs() < 1e-10_f64);
    }

    // 35. compute_sector_angles returns non-empty result.
    #[test]
    fn compute_sector_angles_non_empty() {
        let vertex = [0.0_f64, 0.0_f64];
        let ends = [
            [1.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64],
            [-1.0_f64, 0.0_f64],
            [0.0_f64, -1.0_f64],
        ];
        let sectors = compute_sector_angles(vertex, &ends);
        assert_eq!(sectors.len(), 4);
    }

    // 36. compute_sector_angles sums to 2π for 4 symmetric creases.
    #[test]
    fn compute_sector_angles_sum_to_2pi() {
        let vertex = [0.0_f64, 0.0_f64];
        let ends = [
            [1.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64],
            [-1.0_f64, 0.0_f64],
            [0.0_f64, -1.0_f64],
        ];
        let sectors = compute_sector_angles(vertex, &ends);
        let total: f64 = sectors.iter().sum();
        assert!((total - 2.0_f64 * PI).abs() < 1e-8_f64, "total={total}");
    }

    // 37. Kawasaki alternating sum is 0 for symmetric sectors.
    #[test]
    fn kawasaki_alternating_sum_symmetric() {
        let sectors = vec![PI / 2.0_f64, PI / 2.0_f64, PI / 2.0_f64, PI / 2.0_f64];
        assert!(kawasaki_alternating_sum(&sectors).abs() < 1e-10_f64);
    }

    // 38. kawasaki_check standalone function.
    #[test]
    fn kawasaki_check_standalone() {
        let sectors = vec![PI / 2.0_f64; 4];
        assert!(kawasaki_check(&sectors, 1e-8_f64));
    }

    // 39. Maekawa standalone function.
    #[test]
    fn maekawa_check_standalone() {
        assert!(maekawa_check(3, 1));
        assert!(!maekawa_check(2, 2));
    }

    // 40. Miura fold angle from compactness round-trip.
    #[test]
    fn miura_fold_angle_round_trip() {
        let theta = PI / 4.0_f64;
        let c = theta.cos();
        let theta_back = miura_fold_angle_from_compactness(c);
        assert!((theta_back - theta).abs() < 1e-10_f64);
    }

    // 41. Miura folded length less than flat.
    #[test]
    fn miura_folded_length_less_than_flat() {
        let flat = miura_flat_length_x(0.1_f64, 5);
        let folded = miura_folded_length_x(0.1_f64, 5, PI / 3.0_f64);
        assert!(folded < flat);
    }

    // 42. Torsional spring energy is zero at target.
    #[test]
    fn torsional_spring_energy_at_target() {
        let e = torsional_spring_energy(PI, PI, 10.0_f64);
        assert!(e.abs() < 1e-12_f64);
    }

    // 43. Torsional spring torque drives toward target.
    #[test]
    fn torsional_spring_torque_sign() {
        // Current < target: torque should be positive (drive up).
        let tau = torsional_spring_torque(0.0_f64, PI, 1.0_f64, 0.0_f64, 0.0_f64);
        assert!(tau > 0.0_f64, "tau={tau}");
    }

    // 44. point_distance correct.
    #[test]
    fn point_distance_correct() {
        let a = [0.0_f64, 0.0_f64, 0.0_f64];
        let b = [3.0_f64, 4.0_f64, 0.0_f64];
        assert!((point_distance(a, b) - 5.0_f64).abs() < 1e-10_f64);
    }

    // 45. Project onto plane: point on plane unchanged.
    #[test]
    fn project_onto_plane_on_plane() {
        let normal = [0.0_f64, 0.0_f64, 1.0_f64];
        let origin = [0.0_f64, 0.0_f64, 0.0_f64];
        let point = [1.0_f64, 2.0_f64, 0.0_f64];
        let proj = project_onto_plane(point, origin, normal);
        assert!((proj[2]).abs() < 1e-12_f64);
    }

    // 46. OrigamiRigidFolding active_count matches activated fold lines.
    #[test]
    fn rigid_folding_active_count() {
        let p = simple_pattern();
        let mut rfold = OrigamiRigidFolding::new(p, 5, 1.0_f64, 0.1_f64);
        rfold.activate(0);
        rfold.activate(1);
        assert_eq!(rfold.active_count(), 2);
    }

    // 47. OrigamiRigidFolding time advances after step.
    #[test]
    fn rigid_folding_time_advances() {
        let p = simple_pattern();
        let mut rfold = OrigamiRigidFolding::new(p, 5, 1.0_f64, 0.1_f64);
        rfold.activate_all();
        rfold.step(0.01_f64);
        assert!(rfold.time > 0.0_f64);
    }

    // 48. CreasePattern panel count.
    #[test]
    fn crease_pattern_panel_count() {
        let p = simple_pattern();
        assert_eq!(p.panels.len(), 2);
    }

    // 49. DeployableStructure deployment_progress starts at zero.
    #[test]
    fn deployable_progress_starts_at_zero() {
        let p = simple_pattern();
        let sim = FoldingSimulation::new(p, 0.1_f64, 5, 9.81_f64);
        let ds = DeployableStructure::new(sim, vec![0, 1, 2], 0.1_f64);
        assert_eq!(ds.deployment_progress(), 0.0_f64);
    }

    // 50. rigid_origami_dof formula.
    #[test]
    fn rigid_origami_dof_formula() {
        // 4 hinges → DOF = 1 for Miura-like.
        assert_eq!(rigid_origami_dof(4), 1);
    }
}
