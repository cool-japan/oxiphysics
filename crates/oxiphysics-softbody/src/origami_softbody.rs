// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Origami and foldable structure simulation.
//!
//! Provides:
//!
//! - **Crease pattern representation** with vertices, edges and faces
//! - **Mountain / valley folds** classification and constraint enforcement
//! - **Rigid origami kinematics** using fold angle interpolation
//! - **Fold angle interpolation** (linear, smooth ease-in-out)
//! - **Miura-ori pattern** generation
//! - **Yoshizawa-Randlett** diagram primitives
//! - **Waterbomb** base pattern generation
//! - **Flat-foldability** conditions (Kawasaki's theorem, Maekawa's theorem)
//! - **Origami stiffness** (crease stiffness model)
//! - **Self-folding actuation** (thermal / stimulus-driven fold angle evolution)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector helpers ([f64; 3] — no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-D vectors.
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract.
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale.
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

/// Normalise; returns zero vector if degenerate.
fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = v3_norm(a);
    if l < 1.0e-14 {
        return [0.0, 0.0, 0.0];
    }
    v3_scale(a, 1.0 / l)
}

/// Zero vector.
fn v3_zero() -> [f64; 3] {
    [0.0, 0.0, 0.0]
}

/// Distance between two points.
fn v3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    v3_norm(v3_sub(a, b))
}

// ---------------------------------------------------------------------------
// Crease / fold types
// ---------------------------------------------------------------------------

/// Type of a crease line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CreaseType {
    /// Mountain fold (folds away from viewer).
    Mountain,
    /// Valley fold (folds toward viewer).
    Valley,
    /// Flat / boundary edge (no fold).
    Flat,
    /// Boundary edge of the paper.
    Boundary,
}

impl CreaseType {
    /// Sign convention: mountain = +1, valley = -1, flat/boundary = 0.
    pub fn sign(&self) -> i32 {
        match self {
            CreaseType::Mountain => 1,
            CreaseType::Valley => -1,
            CreaseType::Flat | CreaseType::Boundary => 0,
        }
    }

    /// Whether this crease is foldable (mountain or valley).
    pub fn is_foldable(&self) -> bool {
        matches!(self, CreaseType::Mountain | CreaseType::Valley)
    }
}

/// A crease edge in a crease pattern.
#[derive(Clone, Debug)]
pub struct CreaseEdge {
    /// Index of start vertex.
    pub v0: usize,
    /// Index of end vertex.
    pub v1: usize,
    /// Type of crease.
    pub crease_type: CreaseType,
    /// Target fold angle \[radians\]. Positive = mountain, negative = valley.
    pub target_angle: f64,
    /// Current fold angle \[radians\].
    pub current_angle: f64,
    /// Stiffness of this crease \[N*m/rad\].
    pub stiffness: f64,
    /// Left face index (if any).
    pub face_left: Option<usize>,
    /// Right face index (if any).
    pub face_right: Option<usize>,
}

impl CreaseEdge {
    /// Create a new crease edge.
    pub fn new(v0: usize, v1: usize, crease_type: CreaseType, target_angle: f64) -> Self {
        Self {
            v0,
            v1,
            crease_type,
            target_angle,
            current_angle: 0.0,
            stiffness: 1.0,
            face_left: None,
            face_right: None,
        }
    }

    /// Angle error (difference between current and target).
    pub fn angle_error(&self) -> f64 {
        self.target_angle - self.current_angle
    }

    /// Torque exerted by the crease stiffness.
    pub fn torque(&self) -> f64 {
        self.stiffness * self.angle_error()
    }
}

/// A face (panel) in the crease pattern.
#[derive(Clone, Debug)]
pub struct CreaseFace {
    /// Vertex indices forming this face (polygon, CCW order).
    pub vertices: Vec<usize>,
    /// Edge indices bordering this face.
    pub edges: Vec<usize>,
}

impl CreaseFace {
    /// Create a new face.
    pub fn new(vertices: Vec<usize>, edges: Vec<usize>) -> Self {
        Self { vertices, edges }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Compute area given vertex positions.
    pub fn area(&self, positions: &[[f64; 3]]) -> f64 {
        if self.vertices.len() < 3 {
            return 0.0;
        }
        let mut total = 0.0;
        let a = positions[self.vertices[0]];
        for i in 1..self.vertices.len() - 1 {
            let b = positions[self.vertices[i]];
            let c = positions[self.vertices[i + 1]];
            let ab = v3_sub(b, a);
            let ac = v3_sub(c, a);
            total += v3_norm(v3_cross(ab, ac)) * 0.5;
        }
        total
    }

    /// Compute face normal (average).
    pub fn normal(&self, positions: &[[f64; 3]]) -> [f64; 3] {
        if self.vertices.len() < 3 {
            return [0.0, 0.0, 1.0];
        }
        let a = positions[self.vertices[0]];
        let b = positions[self.vertices[1]];
        let c = positions[self.vertices[2]];
        let ab = v3_sub(b, a);
        let ac = v3_sub(c, a);
        v3_normalize(v3_cross(ab, ac))
    }
}

// ---------------------------------------------------------------------------
// Crease pattern
// ---------------------------------------------------------------------------

/// A complete crease pattern.
#[derive(Clone, Debug)]
pub struct CreasePattern {
    /// Vertex positions (initially flat, in 3-D).
    pub vertices: Vec<[f64; 3]>,
    /// Crease edges.
    pub edges: Vec<CreaseEdge>,
    /// Faces (panels).
    pub faces: Vec<CreaseFace>,
}

impl CreasePattern {
    /// Create an empty crease pattern.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
        }
    }

    /// Add a vertex, return its index.
    pub fn add_vertex(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(pos);
        idx
    }

    /// Add a crease edge, return its index.
    pub fn add_edge(
        &mut self,
        v0: usize,
        v1: usize,
        crease_type: CreaseType,
        target_angle: f64,
    ) -> usize {
        let idx = self.edges.len();
        self.edges
            .push(CreaseEdge::new(v0, v1, crease_type, target_angle));
        idx
    }

    /// Add a face, return its index.
    pub fn add_face(&mut self, vertices: Vec<usize>, edges: Vec<usize>) -> usize {
        let idx = self.faces.len();
        self.faces.push(CreaseFace::new(vertices, edges));
        idx
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    /// Count mountain creases.
    pub fn count_mountains(&self) -> usize {
        self.edges
            .iter()
            .filter(|e| e.crease_type == CreaseType::Mountain)
            .count()
    }

    /// Count valley creases.
    pub fn count_valleys(&self) -> usize {
        self.edges
            .iter()
            .filter(|e| e.crease_type == CreaseType::Valley)
            .count()
    }

    /// Total crease length.
    pub fn total_crease_length(&self) -> f64 {
        self.edges
            .iter()
            .filter(|e| e.crease_type.is_foldable())
            .map(|e| v3_dist(self.vertices[e.v0], self.vertices[e.v1]))
            .sum()
    }

    /// Set all crease stiffnesses.
    pub fn set_all_stiffness(&mut self, stiffness: f64) {
        for e in &mut self.edges {
            e.stiffness = stiffness;
        }
    }

    /// Get edges incident to a vertex.
    pub fn vertex_edges(&self, v: usize) -> Vec<usize> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.v0 == v || e.v1 == v)
            .map(|(i, _)| i)
            .collect()
    }

    /// Get crease types around a vertex.
    pub fn vertex_crease_types(&self, v: usize) -> Vec<CreaseType> {
        self.vertex_edges(v)
            .iter()
            .map(|&ei| self.edges[ei].crease_type)
            .collect()
    }
}

impl Default for CreasePattern {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Flat-foldability conditions
// ---------------------------------------------------------------------------

/// Kawasaki's theorem: at an interior vertex, the alternating sum of sector
/// angles must equal zero (equivalently, opposite angles sum to pi).
///
/// Returns the residual (should be near zero for flat-foldable).
pub fn kawasaki_condition(sector_angles: &[f64]) -> f64 {
    let n = sector_angles.len();
    if n < 4 || !n.is_multiple_of(2) {
        return f64::INFINITY;
    }
    let mut sum_even = 0.0;
    let mut sum_odd = 0.0;
    for (i, &a) in sector_angles.iter().enumerate() {
        if i % 2 == 0 {
            sum_even += a;
        } else {
            sum_odd += a;
        }
    }
    (sum_even - sum_odd).abs()
}

/// Check Kawasaki's theorem: alternating angle sums equal pi.
pub fn is_kawasaki_satisfied(sector_angles: &[f64], tolerance: f64) -> bool {
    kawasaki_condition(sector_angles) < tolerance
}

/// Maekawa's theorem: at any interior vertex, M - V = +/- 2.
///
/// `m` = number of mountain creases, `v` = number of valley creases.
/// Returns true if Maekawa's condition holds.
pub fn maekawa_condition(mountains: usize, valleys: usize) -> bool {
    let diff = (mountains as i64 - valleys as i64).unsigned_abs();
    diff == 2
}

/// Check if a vertex satisfies Maekawa's theorem given crease types.
pub fn check_maekawa(crease_types: &[CreaseType]) -> bool {
    let m = crease_types
        .iter()
        .filter(|&&t| t == CreaseType::Mountain)
        .count();
    let v = crease_types
        .iter()
        .filter(|&&t| t == CreaseType::Valley)
        .count();
    maekawa_condition(m, v)
}

/// Total number of creases at a flat-foldable vertex must be even.
pub fn is_even_degree(crease_types: &[CreaseType]) -> bool {
    let foldable: Vec<_> = crease_types.iter().filter(|t| t.is_foldable()).collect();
    foldable.len() % 2 == 0
}

// ---------------------------------------------------------------------------
// Fold angle interpolation
// ---------------------------------------------------------------------------

/// Linear fold angle interpolation.
pub fn fold_angle_linear(start: f64, end: f64, t: f64) -> f64 {
    start + (end - start) * t.clamp(0.0, 1.0)
}

/// Smooth (ease-in-out) fold angle interpolation using smoothstep.
pub fn fold_angle_smooth(start: f64, end: f64, t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    let s = t * t * (3.0 - 2.0 * t); // smoothstep
    start + (end - start) * s
}

/// Quintic smooth interpolation (smoother ease-in-out).
pub fn fold_angle_quintic(start: f64, end: f64, t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    let s = t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
    start + (end - start) * s
}

/// Sinusoidal fold angle interpolation.
pub fn fold_angle_sinusoidal(start: f64, end: f64, t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    let s = 0.5 * (1.0 - (PI * t).cos());
    start + (end - start) * s
}

// ---------------------------------------------------------------------------
// Rigid origami kinematics
// ---------------------------------------------------------------------------

/// Compute the dihedral angle between two face normals sharing an edge.
pub fn dihedral_angle(n1: [f64; 3], n2: [f64; 3], edge_dir: [f64; 3]) -> f64 {
    let cos_a = v3_dot(n1, n2).clamp(-1.0, 1.0);
    let cross = v3_cross(n1, n2);
    let sin_a = v3_dot(cross, v3_normalize(edge_dir));
    sin_a.atan2(cos_a)
}

/// Rotate a point around an axis (Rodrigues' rotation formula).
pub fn rotate_around_axis(
    point: [f64; 3],
    axis_origin: [f64; 3],
    axis_dir: [f64; 3],
    angle: f64,
) -> [f64; 3] {
    let k = v3_normalize(axis_dir);
    let p = v3_sub(point, axis_origin);
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    // p' = p*cos(a) + (k x p)*sin(a) + k*(k.p)*(1-cos(a))
    let kxp = v3_cross(k, p);
    let kdp = v3_dot(k, p);
    let rotated = v3_add(
        v3_add(v3_scale(p, cos_a), v3_scale(kxp, sin_a)),
        v3_scale(k, kdp * (1.0 - cos_a)),
    );
    v3_add(rotated, axis_origin)
}

/// Fold a set of vertices around a crease edge by a given angle.
///
/// `vertices_to_fold` are indices that should be rotated.
pub fn fold_vertices(
    positions: &mut [[f64; 3]],
    edge_v0: usize,
    edge_v1: usize,
    angle: f64,
    vertices_to_fold: &[usize],
) {
    let origin = positions[edge_v0];
    let axis = v3_sub(positions[edge_v1], positions[edge_v0]);
    for &vi in vertices_to_fold {
        positions[vi] = rotate_around_axis(positions[vi], origin, axis, angle);
    }
}

// ---------------------------------------------------------------------------
// Pattern generators
// ---------------------------------------------------------------------------

/// Generate a Miura-ori crease pattern.
///
/// `nx` x `ny` grid of parallelogram cells, with `angle` controlling the
/// zigzag angle. Returns positions in the flat (unfolded) state.
pub fn miura_ori_pattern(
    nx: usize,
    ny: usize,
    cell_width: f64,
    cell_height: f64,
    zigzag_angle: f64,
) -> CreasePattern {
    let mut cp = CreasePattern::new();
    let cols = nx + 1;
    let rows = ny + 1;

    // Generate vertices
    for j in 0..rows {
        for i in 0..cols {
            let x_offset = if j % 2 == 1 {
                cell_width * zigzag_angle.tan() * 0.5
            } else {
                0.0
            };
            let x = i as f64 * cell_width + x_offset;
            let y = j as f64 * cell_height;
            cp.add_vertex([x, y, 0.0]);
        }
    }

    // Horizontal edges (alternating mountain/valley)
    for j in 0..rows {
        for i in 0..nx {
            let v0 = j * cols + i;
            let v1 = j * cols + i + 1;
            let ct = if j % 2 == 0 {
                CreaseType::Mountain
            } else {
                CreaseType::Valley
            };
            let target = if ct == CreaseType::Mountain {
                PI / 3.0
            } else {
                -PI / 3.0
            };
            cp.add_edge(v0, v1, ct, target);
        }
    }

    // Vertical edges (alternating)
    for j in 0..ny {
        for i in 0..cols {
            let v0 = j * cols + i;
            let v1 = (j + 1) * cols + i;
            let ct = if i % 2 == 0 {
                CreaseType::Mountain
            } else {
                CreaseType::Valley
            };
            let target = if ct == CreaseType::Mountain {
                PI / 3.0
            } else {
                -PI / 3.0
            };
            cp.add_edge(v0, v1, ct, target);
        }
    }

    // Faces
    for j in 0..ny {
        for i in 0..nx {
            let v0 = j * cols + i;
            let v1 = j * cols + i + 1;
            let v2 = (j + 1) * cols + i + 1;
            let v3 = (j + 1) * cols + i;
            cp.add_face(vec![v0, v1, v2, v3], vec![]);
        }
    }

    cp
}

/// Generate a waterbomb base crease pattern.
///
/// A classic square with diagonal and horizontal/vertical creases.
pub fn waterbomb_pattern(size: f64) -> CreasePattern {
    let mut cp = CreasePattern::new();
    let s = size / 2.0;

    // 5 vertices: 4 corners + center
    let v0 = cp.add_vertex([-s, -s, 0.0]); // bottom-left
    let v1 = cp.add_vertex([s, -s, 0.0]); // bottom-right
    let v2 = cp.add_vertex([s, s, 0.0]); // top-right
    let v3 = cp.add_vertex([-s, s, 0.0]); // top-left
    let vc = cp.add_vertex([0.0, 0.0, 0.0]); // center

    // Boundary edges
    cp.add_edge(v0, v1, CreaseType::Boundary, 0.0);
    cp.add_edge(v1, v2, CreaseType::Boundary, 0.0);
    cp.add_edge(v2, v3, CreaseType::Boundary, 0.0);
    cp.add_edge(v3, v0, CreaseType::Boundary, 0.0);

    // Diagonal creases (valleys)
    cp.add_edge(v0, vc, CreaseType::Valley, -PI / 2.0);
    cp.add_edge(v1, vc, CreaseType::Valley, -PI / 2.0);
    cp.add_edge(v2, vc, CreaseType::Valley, -PI / 2.0);
    cp.add_edge(v3, vc, CreaseType::Valley, -PI / 2.0);

    // Mid-edge creases (mountains) - connect midpoints through center
    // For simplicity, add midpoints
    let m01 = cp.add_vertex([0.0, -s, 0.0]);
    let m12 = cp.add_vertex([s, 0.0, 0.0]);
    let m23 = cp.add_vertex([0.0, s, 0.0]);
    let m30 = cp.add_vertex([-s, 0.0, 0.0]);

    cp.add_edge(m01, vc, CreaseType::Mountain, PI / 2.0);
    cp.add_edge(m12, vc, CreaseType::Mountain, PI / 2.0);
    cp.add_edge(m23, vc, CreaseType::Mountain, PI / 2.0);
    cp.add_edge(m30, vc, CreaseType::Mountain, PI / 2.0);

    // Faces (8 triangles)
    cp.add_face(vec![v0, m01, vc], vec![]);
    cp.add_face(vec![m01, v1, vc], vec![]);
    cp.add_face(vec![v1, m12, vc], vec![]);
    cp.add_face(vec![m12, v2, vc], vec![]);
    cp.add_face(vec![v2, m23, vc], vec![]);
    cp.add_face(vec![m23, v3, vc], vec![]);
    cp.add_face(vec![v3, m30, vc], vec![]);
    cp.add_face(vec![m30, v0, vc], vec![]);

    cp
}

/// Generate a Yoshizawa-Randlett diagram basic fold sequence.
///
/// Returns a simple square paper with a single valley fold along the middle.
pub fn yoshizawa_valley_fold(size: f64) -> CreasePattern {
    let mut cp = CreasePattern::new();
    let s = size / 2.0;

    let v0 = cp.add_vertex([-s, -s, 0.0]);
    let v1 = cp.add_vertex([s, -s, 0.0]);
    let v2 = cp.add_vertex([s, s, 0.0]);
    let v3 = cp.add_vertex([-s, s, 0.0]);
    let m_left = cp.add_vertex([-s, 0.0, 0.0]);
    let m_right = cp.add_vertex([s, 0.0, 0.0]);

    // Boundary
    cp.add_edge(v0, v1, CreaseType::Boundary, 0.0);
    cp.add_edge(v1, m_right, CreaseType::Boundary, 0.0);
    cp.add_edge(m_right, v2, CreaseType::Boundary, 0.0);
    cp.add_edge(v2, v3, CreaseType::Boundary, 0.0);
    cp.add_edge(v3, m_left, CreaseType::Boundary, 0.0);
    cp.add_edge(m_left, v0, CreaseType::Boundary, 0.0);

    // Valley fold
    cp.add_edge(m_left, m_right, CreaseType::Valley, -PI);

    // Two faces
    cp.add_face(vec![v0, v1, m_right, m_left], vec![]);
    cp.add_face(vec![m_left, m_right, v2, v3], vec![]);

    cp
}

/// Generate a Yoshizawa-Randlett mountain fold sequence.
pub fn yoshizawa_mountain_fold(size: f64) -> CreasePattern {
    let mut cp = yoshizawa_valley_fold(size);
    // Convert the valley fold to mountain
    for edge in &mut cp.edges {
        if edge.crease_type == CreaseType::Valley {
            edge.crease_type = CreaseType::Mountain;
            edge.target_angle = PI;
        }
    }
    cp
}

// ---------------------------------------------------------------------------
// Origami stiffness model
// ---------------------------------------------------------------------------

/// Crease stiffness model: torque as a function of fold angle.
#[derive(Clone, Debug)]
pub struct CreaseStiffnessModel {
    /// Linear stiffness coefficient \[N*m/rad\].
    pub linear_stiffness: f64,
    /// Nonlinear (cubic) stiffness coefficient.
    pub cubic_stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
    /// Rest angle.
    pub rest_angle: f64,
}

impl CreaseStiffnessModel {
    /// Create a linear stiffness model.
    pub fn linear(stiffness: f64) -> Self {
        Self {
            linear_stiffness: stiffness,
            cubic_stiffness: 0.0,
            damping: 0.0,
            rest_angle: 0.0,
        }
    }

    /// Create a nonlinear stiffness model.
    pub fn nonlinear(linear: f64, cubic: f64, damping: f64) -> Self {
        Self {
            linear_stiffness: linear,
            cubic_stiffness: cubic,
            damping,
            rest_angle: 0.0,
        }
    }

    /// Compute torque for a given angle and angular velocity.
    pub fn torque(&self, angle: f64, angular_velocity: f64) -> f64 {
        let da = angle - self.rest_angle;
        let linear = self.linear_stiffness * da;
        let cubic = self.cubic_stiffness * da * da * da;
        let damp = self.damping * angular_velocity;
        -(linear + cubic + damp)
    }

    /// Potential energy stored in the crease.
    pub fn potential_energy(&self, angle: f64) -> f64 {
        let da = angle - self.rest_angle;
        0.5 * self.linear_stiffness * da * da + 0.25 * self.cubic_stiffness * da.powi(4)
    }
}

// ---------------------------------------------------------------------------
// Self-folding actuation
// ---------------------------------------------------------------------------

/// Self-folding actuator model.
///
/// Drives fold angle toward a target via a stimulus (e.g. temperature).
#[derive(Clone, Debug)]
pub struct SelfFoldingActuator {
    /// Target fold angle.
    pub target_angle: f64,
    /// Actuation rate \[rad/s per unit stimulus\].
    pub rate: f64,
    /// Current stimulus level (0 = off, 1 = fully on).
    pub stimulus: f64,
    /// Maximum angular velocity \[rad/s\].
    pub max_velocity: f64,
}

impl SelfFoldingActuator {
    /// Create a new actuator.
    pub fn new(target_angle: f64, rate: f64) -> Self {
        Self {
            target_angle,
            rate,
            stimulus: 0.0,
            max_velocity: 10.0,
        }
    }

    /// Set stimulus level (clamped to \[0, 1\]).
    pub fn set_stimulus(&mut self, level: f64) {
        self.stimulus = level.clamp(0.0, 1.0);
    }

    /// Compute the desired angular velocity given current angle.
    pub fn desired_velocity(&self, current_angle: f64) -> f64 {
        let error = self.target_angle - current_angle;
        let vel = self.stimulus * self.rate * error;
        vel.clamp(-self.max_velocity, self.max_velocity)
    }

    /// Update the fold angle by one time step.
    pub fn update(&self, current_angle: f64, dt: f64) -> f64 {
        let vel = self.desired_velocity(current_angle);
        current_angle + vel * dt
    }
}

// ---------------------------------------------------------------------------
// Origami simulator
// ---------------------------------------------------------------------------

/// Simple origami folding simulator.
#[derive(Clone, Debug)]
pub struct OrigamiSimulator {
    /// The crease pattern.
    pub pattern: CreasePattern,
    /// Fold angle velocities.
    pub angular_velocities: Vec<f64>,
    /// Stiffness model for each crease.
    pub stiffness_models: Vec<CreaseStiffnessModel>,
    /// Optional self-folding actuators.
    pub actuators: Vec<Option<SelfFoldingActuator>>,
    /// Current simulation time.
    pub time: f64,
    /// Global damping.
    pub damping: f64,
}

impl OrigamiSimulator {
    /// Create a new simulator from a crease pattern.
    pub fn new(pattern: CreasePattern, stiffness: f64) -> Self {
        let n = pattern.edges.len();
        let stiffness_models: Vec<_> = (0..n)
            .map(|_| CreaseStiffnessModel::linear(stiffness))
            .collect();
        Self {
            angular_velocities: vec![0.0; n],
            actuators: vec![None; n],
            stiffness_models,
            pattern,
            time: 0.0,
            damping: 0.5,
        }
    }

    /// Set an actuator on a specific edge.
    pub fn set_actuator(&mut self, edge_idx: usize, actuator: SelfFoldingActuator) {
        if edge_idx < self.actuators.len() {
            self.actuators[edge_idx] = Some(actuator);
        }
    }

    /// Step the simulation forward by `dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.pattern.edges.len();
        for i in 0..n {
            let edge = &self.pattern.edges[i];
            if !edge.crease_type.is_foldable() {
                continue;
            }

            // Stiffness torque toward target
            let torque =
                self.stiffness_models[i].torque(edge.current_angle, self.angular_velocities[i]);

            // Actuator contribution
            let actuator_vel = if let Some(ref act) = self.actuators[i] {
                act.desired_velocity(edge.current_angle)
            } else {
                0.0
            };

            // Simple integration
            let total_torque = torque - self.damping * self.angular_velocities[i];
            self.angular_velocities[i] += total_torque * dt + actuator_vel * dt;
            self.pattern.edges[i].current_angle += self.angular_velocities[i] * dt;
        }
        self.time += dt;
    }

    /// Run simulation for `n_steps`.
    pub fn run(&mut self, n_steps: usize, dt: f64) {
        for _ in 0..n_steps {
            self.step(dt);
        }
    }

    /// Compute total potential energy in all creases.
    pub fn total_energy(&self) -> f64 {
        let mut total = 0.0;
        for (i, edge) in self.pattern.edges.iter().enumerate() {
            if edge.crease_type.is_foldable() {
                total += self.stiffness_models[i].potential_energy(edge.current_angle);
            }
        }
        total
    }

    /// Get maximum angle error across all foldable creases.
    pub fn max_angle_error(&self) -> f64 {
        self.pattern
            .edges
            .iter()
            .filter(|e| e.crease_type.is_foldable())
            .map(|e| e.angle_error().abs())
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Fold sequence
// ---------------------------------------------------------------------------

/// A single fold step in a fold sequence.
#[derive(Clone, Debug)]
pub struct FoldStep {
    /// Edge index to fold.
    pub edge_index: usize,
    /// Target angle for this step.
    pub target_angle: f64,
    /// Vertices to rotate.
    pub affected_vertices: Vec<usize>,
}

/// Execute a fold sequence on a crease pattern.
pub fn execute_fold_sequence(pattern: &mut CreasePattern, steps: &[FoldStep]) {
    for step in steps {
        let edge = &pattern.edges[step.edge_index];
        let angle_delta = step.target_angle - edge.current_angle;
        fold_vertices(
            &mut pattern.vertices,
            edge.v0,
            edge.v1,
            angle_delta,
            &step.affected_vertices,
        );
        pattern.edges[step.edge_index].current_angle = step.target_angle;
    }
}

// ---------------------------------------------------------------------------
// Utility: measure folded state
// ---------------------------------------------------------------------------

/// Bounding box of the origami in its current configuration.
pub fn bounding_box(positions: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    if positions.is_empty() {
        return (v3_zero(), v3_zero());
    }
    let mut min_pt = positions[0];
    let mut max_pt = positions[0];
    for p in positions {
        for k in 0..3 {
            if p[k] < min_pt[k] {
                min_pt[k] = p[k];
            }
            if p[k] > max_pt[k] {
                max_pt[k] = p[k];
            }
        }
    }
    (min_pt, max_pt)
}

/// Compute the compaction ratio (folded bbox volume / flat bbox volume).
pub fn compaction_ratio(flat_positions: &[[f64; 3]], folded_positions: &[[f64; 3]]) -> f64 {
    let (fmin, fmax) = bounding_box(flat_positions);
    let (cmin, cmax) = bounding_box(folded_positions);

    let flat_vol = (fmax[0] - fmin[0]).max(1e-12)
        * (fmax[1] - fmin[1]).max(1e-12)
        * (fmax[2] - fmin[2]).max(1e-12);
    let fold_vol = (cmax[0] - cmin[0]).max(1e-12)
        * (cmax[1] - cmin[1]).max(1e-12)
        * (cmax[2] - cmin[2]).max(1e-12);

    fold_vol / flat_vol
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // 1. CreaseType signs
    #[test]
    fn test_crease_type_signs() {
        assert_eq!(CreaseType::Mountain.sign(), 1);
        assert_eq!(CreaseType::Valley.sign(), -1);
        assert_eq!(CreaseType::Flat.sign(), 0);
        assert_eq!(CreaseType::Boundary.sign(), 0);
    }

    // 2. CreaseType foldable
    #[test]
    fn test_crease_foldable() {
        assert!(CreaseType::Mountain.is_foldable());
        assert!(CreaseType::Valley.is_foldable());
        assert!(!CreaseType::Flat.is_foldable());
        assert!(!CreaseType::Boundary.is_foldable());
    }

    // 3. Kawasaki condition: four 90-degree sectors
    #[test]
    fn test_kawasaki_90_deg() {
        let angles = vec![PI / 2.0, PI / 2.0, PI / 2.0, PI / 2.0];
        let residual = kawasaki_condition(&angles);
        assert!(
            residual < 1.0e-10,
            "Four equal sectors should satisfy Kawasaki, got {:.6}",
            residual
        );
    }

    // 4. Kawasaki condition: non-flat-foldable
    #[test]
    fn test_kawasaki_fails() {
        let angles = vec![PI / 3.0, PI / 2.0, PI / 3.0, PI / 2.0 + PI / 6.0];
        let residual = kawasaki_condition(&angles);
        // sum_even = pi/3 + pi/3 = 2pi/3, sum_odd = pi/2 + pi/2+pi/6 = pi + pi/6
        // |2pi/3 - 7pi/6| != 0
        assert!(residual > 0.1, "Should fail Kawasaki");
    }

    // 5. Maekawa condition: 3M + 1V
    #[test]
    fn test_maekawa_3m_1v() {
        assert!(maekawa_condition(3, 1));
    }

    // 6. Maekawa condition: 2M + 2V fails
    #[test]
    fn test_maekawa_2m_2v_fails() {
        assert!(!maekawa_condition(2, 2));
    }

    // 7. Even degree
    #[test]
    fn test_even_degree() {
        let types = vec![
            CreaseType::Mountain,
            CreaseType::Valley,
            CreaseType::Mountain,
            CreaseType::Valley,
        ];
        assert!(is_even_degree(&types));
    }

    // 8. Linear fold angle interpolation
    #[test]
    fn test_fold_linear() {
        assert!((fold_angle_linear(0.0, PI, 0.0)).abs() < EPS);
        assert!((fold_angle_linear(0.0, PI, 1.0) - PI).abs() < EPS);
        assert!((fold_angle_linear(0.0, PI, 0.5) - PI / 2.0).abs() < EPS);
    }

    // 9. Smooth interpolation endpoints
    #[test]
    fn test_fold_smooth_endpoints() {
        assert!((fold_angle_smooth(0.0, 1.0, 0.0)).abs() < EPS);
        assert!((fold_angle_smooth(0.0, 1.0, 1.0) - 1.0).abs() < EPS);
    }

    // 10. Smooth interpolation is monotonic
    #[test]
    fn test_fold_smooth_monotonic() {
        let mut prev = fold_angle_smooth(0.0, 1.0, 0.0);
        for i in 1..=20 {
            let t = i as f64 / 20.0;
            let val = fold_angle_smooth(0.0, 1.0, t);
            assert!(val >= prev - EPS, "Smooth should be monotonic");
            prev = val;
        }
    }

    // 11. Quintic endpoints
    #[test]
    fn test_quintic_endpoints() {
        assert!((fold_angle_quintic(0.0, 2.0, 0.0)).abs() < EPS);
        assert!((fold_angle_quintic(0.0, 2.0, 1.0) - 2.0).abs() < EPS);
    }

    // 12. Sinusoidal interpolation midpoint
    #[test]
    fn test_sinusoidal_midpoint() {
        let mid = fold_angle_sinusoidal(0.0, 2.0, 0.5);
        assert!(
            (mid - 1.0).abs() < 1.0e-6,
            "Midpoint should be 1.0, got {:.6}",
            mid
        );
    }

    // 13. Rotation: identity (angle=0)
    #[test]
    fn test_rotate_identity() {
        let p = [1.0, 2.0, 3.0];
        let r = rotate_around_axis(p, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0);
        for k in 0..3 {
            assert!((r[k] - p[k]).abs() < 1.0e-10, "Should be identity");
        }
    }

    // 14. Rotation: 90 degrees around Z
    #[test]
    fn test_rotate_90_z() {
        let p = [1.0, 0.0, 0.0];
        let r = rotate_around_axis(p, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], PI / 2.0);
        assert!((r[0]).abs() < 1.0e-10, "x should be ~0, got {:.6}", r[0]);
        assert!(
            (r[1] - 1.0).abs() < 1.0e-10,
            "y should be ~1, got {:.6}",
            r[1]
        );
    }

    // 15. Miura-ori pattern topology
    #[test]
    fn test_miura_ori_topology() {
        let cp = miura_ori_pattern(3, 2, 1.0, 1.0, 0.3);
        assert_eq!(cp.num_vertices(), 4 * 3); // (3+1)*(2+1) = 12
        assert_eq!(cp.num_faces(), 6); // 3*2 = 6
        assert!(cp.num_edges() > 0);
    }

    // 16. Miura-ori has mountains and valleys
    #[test]
    fn test_miura_ori_crease_types() {
        let cp = miura_ori_pattern(3, 2, 1.0, 1.0, 0.3);
        assert!(cp.count_mountains() > 0);
        assert!(cp.count_valleys() > 0);
    }

    // 17. Waterbomb pattern topology
    #[test]
    fn test_waterbomb_topology() {
        let cp = waterbomb_pattern(2.0);
        assert_eq!(
            cp.num_faces(),
            8,
            "Waterbomb should have 8 triangular faces"
        );
        assert!(cp.count_mountains() > 0);
        assert!(cp.count_valleys() > 0);
    }

    // 18. Yoshizawa valley fold
    #[test]
    fn test_yoshizawa_valley() {
        let cp = yoshizawa_valley_fold(2.0);
        assert_eq!(cp.num_faces(), 2);
        assert_eq!(cp.count_valleys(), 1);
    }

    // 19. Yoshizawa mountain fold
    #[test]
    fn test_yoshizawa_mountain() {
        let cp = yoshizawa_mountain_fold(2.0);
        assert_eq!(cp.count_mountains(), 1);
        assert_eq!(cp.count_valleys(), 0);
    }

    // 20. Crease stiffness: zero torque at rest
    #[test]
    fn test_stiffness_zero_at_rest() {
        let model = CreaseStiffnessModel::linear(10.0);
        let t = model.torque(0.0, 0.0);
        assert!(t.abs() < EPS, "Torque at rest should be zero");
    }

    // 21. Crease stiffness: nonzero torque when displaced
    #[test]
    fn test_stiffness_nonzero_displaced() {
        let model = CreaseStiffnessModel::linear(10.0);
        let t = model.torque(0.5, 0.0);
        assert!(t.abs() > EPS, "Torque should be nonzero when displaced");
    }

    // 22. Potential energy at rest is zero
    #[test]
    fn test_potential_energy_rest() {
        let model = CreaseStiffnessModel::linear(10.0);
        assert!(model.potential_energy(0.0) < EPS);
    }

    // 23. Self-folding actuator: zero stimulus => zero velocity
    #[test]
    fn test_actuator_zero_stimulus() {
        let act = SelfFoldingActuator::new(PI, 1.0);
        let vel = act.desired_velocity(0.0);
        assert!(vel.abs() < EPS, "Zero stimulus should give zero velocity");
    }

    // 24. Self-folding actuator: full stimulus drives toward target
    #[test]
    fn test_actuator_full_stimulus() {
        let mut act = SelfFoldingActuator::new(PI, 1.0);
        act.set_stimulus(1.0);
        let vel = act.desired_velocity(0.0);
        assert!(vel > 0.0, "Should drive toward positive target");
    }

    // 25. Simulator doesn't crash
    #[test]
    fn test_simulator_basic() {
        let cp = yoshizawa_valley_fold(1.0);
        let mut sim = OrigamiSimulator::new(cp, 5.0);
        sim.run(100, 0.001);
        assert!(sim.time > 0.0);
    }

    // 26. Bounding box of unit square
    #[test]
    fn test_bounding_box() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let (mn, mx) = bounding_box(&pts);
        assert!((mn[0]).abs() < EPS && (mn[1]).abs() < EPS);
        assert!((mx[0] - 1.0).abs() < EPS && (mx[1] - 1.0).abs() < EPS);
    }

    // 27. Compaction ratio of flat = 1
    #[test]
    fn test_compaction_flat() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let ratio = compaction_ratio(&pts, &pts);
        assert!(
            (ratio - 1.0).abs() < 1.0e-6,
            "Flat compaction should be 1.0, got {:.6}",
            ratio
        );
    }

    // 28. Dihedral angle between parallel normals is 0
    #[test]
    fn test_dihedral_parallel() {
        let n1 = [0.0, 0.0, 1.0];
        let n2 = [0.0, 0.0, 1.0];
        let edge = [1.0, 0.0, 0.0];
        let a = dihedral_angle(n1, n2, edge);
        assert!(
            a.abs() < 1.0e-10,
            "Parallel normals => dihedral 0, got {:.6}",
            a
        );
    }

    // 29. CreaseEdge torque at target is zero
    #[test]
    fn test_crease_edge_torque_at_target() {
        let mut e = CreaseEdge::new(0, 1, CreaseType::Mountain, PI / 2.0);
        e.current_angle = PI / 2.0;
        assert!(e.torque().abs() < EPS);
    }

    // 30. Face area computation
    #[test]
    fn test_face_area() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let face = CreaseFace::new(vec![0, 1, 2], vec![]);
        let area = face.area(&positions);
        assert!(
            (area - 0.5).abs() < 1.0e-10,
            "Triangle area should be 0.5, got {:.6}",
            area
        );
    }

    // 31. Total crease length
    #[test]
    fn test_total_crease_length() {
        let mut cp = CreasePattern::new();
        cp.add_vertex([0.0, 0.0, 0.0]);
        cp.add_vertex([1.0, 0.0, 0.0]);
        cp.add_vertex([0.0, 1.0, 0.0]);
        cp.add_edge(0, 1, CreaseType::Mountain, PI / 2.0);
        cp.add_edge(1, 2, CreaseType::Valley, -PI / 2.0);
        cp.add_edge(0, 2, CreaseType::Boundary, 0.0);
        let length = cp.total_crease_length();
        // Mountain edge: 1.0, Valley edge: sqrt(2) ~ 1.414
        let expected = 1.0 + 2.0_f64.sqrt();
        assert!(
            (length - expected).abs() < 1.0e-6,
            "Total crease length should be {:.6}, got {:.6}",
            expected,
            length
        );
    }
}
