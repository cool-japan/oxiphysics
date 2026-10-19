// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Origami mechanics simulation module.
//!
//! Implements crease pattern graphs, Miura-ori fold simulation, rigid origami
//! kinematics, compliant fold hinge mechanics, waterbomb base, Yoshimura
//! patterns, deployable structures, bistable origami (snap-through),
//! folded-sheet stiffness, and SMA-driven origami actuators.
//!
//! All vectors are plain `[f64; 3]` — no nalgebra dependency.

use std::f64::consts::PI;

use crate::{Error, Result};

// ── vector helpers ────────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-vector. Returns zero vector if norm is too small.
#[inline]
pub fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

/// Compute the dihedral angle between two triangles sharing an edge.
///
/// The edge is defined by vertices `e0`..`e1`. `p0` is a vertex on the first
/// triangle (not on the edge) and `p1` is a vertex on the second triangle.
/// Returns the angle in radians in \[0, π\].
pub fn dihedral_angle(e0: [f64; 3], e1: [f64; 3], p0: [f64; 3], p1: [f64; 3]) -> f64 {
    let edge = vec3_sub(e1, e0);
    let n0 = vec3_cross(vec3_sub(p0, e0), edge);
    let n1 = vec3_cross(edge, vec3_sub(p1, e0));
    let cos_a = vec3_dot(n0, n1) / (vec3_norm(n0) * vec3_norm(n1) + 1e-15);
    // Convention: flat (coplanar) panels give dihedral angle = π; fully folded = 0.
    PI - cos_a.clamp(-1.0, 1.0).acos()
}

// ── crease pattern ────────────────────────────────────────────────────────────

/// Assignment of a crease in an origami pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CreaseAssignment {
    /// Mountain fold (convex when viewed from above).
    Mountain,
    /// Valley fold (concave when viewed from above).
    Valley,
    /// Flat crease or boundary.
    Flat,
    /// Border/boundary edge.
    Border,
}

/// A single crease (fold line) in a crease pattern.
#[derive(Debug, Clone)]
pub struct Crease {
    /// Index of the start vertex.
    pub v0: usize,
    /// Index of the end vertex.
    pub v1: usize,
    /// Fold assignment.
    pub assignment: CreaseAssignment,
    /// Rest (target) fold angle in radians.
    pub rest_angle: f64,
}

impl Crease {
    /// Create a new crease.
    pub fn new(v0: usize, v1: usize, assignment: CreaseAssignment, rest_angle: f64) -> Self {
        Self {
            v0,
            v1,
            assignment,
            rest_angle,
        }
    }
}

/// A crease pattern graph: vertices in 2-D (x, y) and a list of creases.
#[derive(Debug, Clone)]
pub struct CreasePattern {
    /// Flat (x, y) positions of pattern vertices.
    pub vertices: Vec<[f64; 2]>,
    /// Crease edges.
    pub creases: Vec<Crease>,
}

impl CreasePattern {
    /// Create an empty crease pattern.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            creases: Vec::new(),
        }
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, x: f64, y: f64) -> usize {
        let idx = self.vertices.len();
        self.vertices.push([x, y]);
        idx
    }

    /// Add a crease between two vertex indices.
    pub fn add_crease(
        &mut self,
        v0: usize,
        v1: usize,
        assignment: CreaseAssignment,
        rest_angle: f64,
    ) {
        self.creases
            .push(Crease::new(v0, v1, assignment, rest_angle));
    }

    /// Count mountain folds.
    pub fn count_mountain(&self) -> usize {
        self.creases
            .iter()
            .filter(|c| c.assignment == CreaseAssignment::Mountain)
            .count()
    }

    /// Count valley folds.
    pub fn count_valley(&self) -> usize {
        self.creases
            .iter()
            .filter(|c| c.assignment == CreaseAssignment::Valley)
            .count()
    }

    /// Verify Maekawa's theorem at a given interior vertex: |M - V| == 2.
    ///
    /// Returns `true` if the theorem holds, `false` otherwise.
    pub fn check_maekawa(&self, vertex: usize) -> bool {
        let mut m = 0i32;
        let mut v = 0i32;
        for c in &self.creases {
            if c.v0 == vertex || c.v1 == vertex {
                match c.assignment {
                    CreaseAssignment::Mountain => m += 1,
                    CreaseAssignment::Valley => v += 1,
                    _ => {}
                }
            }
        }
        (m - v).abs() == 2
    }
}

impl Default for CreasePattern {
    fn default() -> Self {
        Self::new()
    }
}

// ── Miura-ori ─────────────────────────────────────────────────────────────────

/// Dimensions of the folded Miura-ori cell (Schenk & Guest, PNAS 2013).
///
/// Around one vertex the four creases of a rigidly folded Miura-ori are
///
/// * `c1 = ( L, 0, H)` and `c3 = (-L, 0, H)`: the two halves of the straight
///   crease line (edge length `a`), which rise out of the reference plane,
/// * `c2 = ( V,  S, 0)` and `c4 = ( V, -S, 0)`: the two halves of the zigzag
///   crease line (edge length `b`), which stay in the reference plane,
///
/// so `L`, `S` and `H` are the extents of the cell along x, y and z and `V`
/// is the zigzag amplitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MiuraCellDimensions {
    /// `S`: y-extent (advance along the zigzag) of one `b`-edge.
    pub s: f64,
    /// `L`: x-extent of one `a`-edge.
    pub l: f64,
    /// `H`: height of the folded cell, the z-extent of one `a`-edge.
    pub h: f64,
    /// `V`: x-extent (zigzag amplitude) of one `b`-edge.
    pub v: f64,
}

/// Miura-ori unit cell descriptor.
///
/// A Miura-ori sheet tiles the plane with parallelogram facets of side
/// lengths `a` (along the straight creases, the x-direction) and `b` (along
/// the zigzag creases, the y-direction) and acute sector angle `gamma`.
/// One unit cell is one period of the folded pattern: `2 × 2` facets.
///
/// `fold_angle` is the angle `θ` between the facets and the reference plane
/// (`θ = 0` flat, `θ = π/2` flat-folded).  With `φ` the elevation of the
/// straight creases (`sin φ = sin θ sin γ`) the cell dimensions are
///
/// ```text
/// H = a sin θ sin γ
/// S = b cos θ tan γ / sqrt(1 + cos²θ tan²γ)
/// L = a sqrt(1 - sin²θ sin²γ)
/// V = b / sqrt(1 + cos²θ tan²γ)
/// ```
///
/// (Schenk & Guest, "Geometry of Miura-folded metamaterials", PNAS 2013),
/// and the projected cell is `2L × 2S`.
#[derive(Debug, Clone)]
pub struct MiuraOriUnit {
    /// Side length `a` of a facet along the x-direction (straight creases).
    pub a: f64,
    /// Side length `b` of a facet along the zigzag (y-direction) creases.
    pub b: f64,
    /// Acute angle of the parallelogram facets (radians).
    pub gamma: f64,
    /// Current fold angle (θ), radians, in \[0, π/2\].  0 = flat, π/2 = flat-folded.
    pub fold_angle: f64,
}

impl MiuraOriUnit {
    /// Create a new Miura-ori unit cell.
    pub fn new(a: f64, b: f64, gamma: f64) -> Self {
        Self {
            a,
            b,
            gamma,
            fold_angle: 0.0,
        }
    }

    /// Dimensions `S`, `L`, `H`, `V` of the folded cell at the current fold angle.
    pub fn dimensions(&self) -> MiuraCellDimensions {
        let (sin_t, cos_t) = self.fold_angle.sin_cos();
        let (sin_g, cos_g) = self.gamma.sin_cos();
        // Elevation of the straight creases: sin φ = sin θ sin γ.
        let l = self.a * (1.0 - sin_t * sin_t * sin_g * sin_g).max(0.0).sqrt();
        let h = self.a * sin_t * sin_g;
        // 1 + cos²θ tan²γ = (cos²γ + cos²θ sin²γ) / cos²γ, which also stays
        // well defined for γ = π/2 (where tan γ is infinite).
        let planar = cos_g * cos_g + cos_t * cos_t * sin_g * sin_g;
        let norm = planar.sqrt();
        // Only |cos γ| enters: a facet angle γ and its supplement describe the
        // same parallelogram.  For rectangular facets (γ = π/2) folded flat
        // (θ = π/2) the expression is 0/0; its limit is the unfolded zigzag
        // with S = b and V = 0.
        let (s, v) = if norm > 0.0 {
            (self.b * cos_t * sin_g / norm, self.b * cos_g.abs() / norm)
        } else {
            (self.b, 0.0)
        };
        MiuraCellDimensions { s, l, h, v }
    }

    /// Projected width of the unit cell (x-extent, one period `2L`).
    ///
    /// `W = 2 a sqrt(1 - sin²θ sin²γ)`; it equals `2a` for the flat sheet.
    pub fn projected_width(&self) -> f64 {
        2.0 * self.dimensions().l
    }

    /// Projected height of the unit cell (y-extent, one zigzag period `2S`).
    ///
    /// `2S = 2 b cos θ tan γ / sqrt(1 + cos²θ tan²γ)`; it equals `2 b sin γ`
    /// for the flat sheet and vanishes for the flat-folded state `θ = π/2`.
    pub fn projected_height(&self) -> f64 {
        2.0 * self.dimensions().s
    }

    /// Out-of-plane height `H = a sin θ sin γ` of the folded cell.
    pub fn height(&self) -> f64 {
        self.dimensions().h
    }

    /// In-plane Poisson ratio `ν = -(dS/S) / (dL/L)` of the Miura-ori sheet
    /// (strain of the y-extent over strain of the x-extent).
    ///
    /// From `d ln S / dθ = -tan θ / (1 + cos²θ tan²γ)` and
    /// `d ln L / dθ = -sin θ cos θ sin²γ / (1 - sin²θ sin²γ)` the θ-derivatives
    /// combine, using `1 - sin²θ sin²γ = cos²γ + cos²θ sin²γ`, to
    ///
    /// ```text
    /// ν = -cos²γ / (cos²θ sin²γ) = -1 / (cos²θ tan²γ)
    /// ```
    ///
    /// which is negative for every fold angle (the sheet is auxetic) and
    /// independent of the edge lengths `a` and `b`.  The value is `NaN` where
    /// the strains are undefined: for the flat-folded state `θ = π/2` and for
    /// a degenerate facet angle (`sin γ = 0`).
    pub fn poisson_ratio(&self) -> f64 {
        let (sin_g, cos_g) = self.gamma.sin_cos();
        let cos_t = self.fold_angle.cos();
        let denominator = cos_t * cos_t * sin_g * sin_g;
        if denominator < 1e-24 {
            return f64::NAN;
        }
        -(cos_g * cos_g) / denominator
    }

    /// Advance fold by `delta_theta` radians, clamping to \[0, π/2\].
    pub fn fold_step(&mut self, delta_theta: f64) {
        self.fold_angle = (self.fold_angle + delta_theta).clamp(0.0, PI / 2.0);
    }
}

/// A tiled Miura-ori sheet of `nx × ny` unit cells.
#[derive(Debug, Clone)]
pub struct MiuraOriSheet {
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Unit cell geometry.
    pub unit: MiuraOriUnit,
}

impl MiuraOriSheet {
    /// Create a new Miura-ori sheet.
    pub fn new(nx: usize, ny: usize, a: f64, b: f64, gamma: f64) -> Self {
        Self {
            nx,
            ny,
            unit: MiuraOriUnit::new(a, b, gamma),
        }
    }

    /// Total projected width of the sheet: `nx` unit cells of width `2L`.
    ///
    /// The sideways overhang `V` of the zigzag boundary is not included.
    pub fn total_width(&self) -> f64 {
        self.nx as f64 * self.unit.projected_width()
    }

    /// Total projected height of the sheet: `ny` unit cells of height `2S`.
    pub fn total_height(&self) -> f64 {
        self.ny as f64 * self.unit.projected_height()
    }

    /// Apply a uniform fold angle (clamped to \[0, π/2\]) to the entire sheet.
    pub fn set_fold_angle(&mut self, theta: f64) {
        self.unit.fold_angle = theta.clamp(0.0, PI / 2.0);
    }

    /// Fold uniformly by `delta_theta` radians.
    pub fn fold_step(&mut self, delta_theta: f64) {
        self.unit.fold_step(delta_theta);
    }
}

// ── rigid origami panel-hinge system ─────────────────────────────────────────

/// A rigid panel in a rigid-origami kinematic model.
#[derive(Debug, Clone)]
pub struct RigidPanel {
    /// Panel index.
    pub id: usize,
    /// Corner vertices in 3-D world coordinates.
    pub vertices: Vec<[f64; 3]>,
    /// Current orientation as a 3×3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Translation of the panel's reference point.
    pub translation: [f64; 3],
}

impl RigidPanel {
    /// Create a flat panel from a list of rest vertices (coplanar).
    pub fn new(id: usize, vertices: Vec<[f64; 3]>) -> Self {
        let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        Self {
            id,
            vertices,
            rotation: rot,
            translation: [0.0; 3],
        }
    }

    /// Compute the panel's centroid.
    pub fn centroid(&self) -> [f64; 3] {
        let n = self.vertices.len() as f64;
        let mut c = [0.0f64; 3];
        for v in &self.vertices {
            c[0] += v[0];
            c[1] += v[1];
            c[2] += v[2];
        }
        [c[0] / n, c[1] / n, c[2] / n]
    }

    /// Apply a rotation matrix `r` to all vertices around the centroid.
    pub fn apply_rotation(&mut self, r: [[f64; 3]; 3]) {
        let cen = self.centroid();
        for v in &mut self.vertices {
            let local = vec3_sub(*v, cen);
            let rotated = mat3_mul_vec(r, local);
            *v = vec3_add(rotated, cen);
        }
    }
}

/// Apply a 3×3 matrix (row-major) to a 3-vector.
fn mat3_mul_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Build a rotation matrix around an arbitrary axis `axis` by angle `theta`.
///
/// Uses Rodrigues' rotation formula.
pub fn rotation_matrix_axis_angle(axis: [f64; 3], theta: f64) -> [[f64; 3]; 3] {
    let [kx, ky, kz] = vec3_normalize(axis);
    let c = theta.cos();
    let s = theta.sin();
    let t = 1.0 - c;
    [
        [t * kx * kx + c, t * kx * ky - s * kz, t * kx * kz + s * ky],
        [t * kx * ky + s * kz, t * ky * ky + c, t * ky * kz - s * kx],
        [t * kx * kz - s * ky, t * ky * kz + s * kx, t * kz * kz + c],
    ]
}

/// A hinge connecting two panels along a shared fold edge.
#[derive(Debug, Clone)]
pub struct FoldHinge {
    /// Index of the first panel.
    pub panel_a: usize,
    /// Index of the second panel.
    pub panel_b: usize,
    /// Fold axis in world coordinates (unit vector along the shared edge).
    pub axis: [f64; 3],
    /// Current fold angle (dihedral angle between the two panels), radians.
    pub angle: f64,
    /// Rest (target) fold angle, radians.
    pub rest_angle: f64,
    /// Hinge stiffness (N·m/rad).
    pub stiffness: f64,
}

impl FoldHinge {
    /// Create a new fold hinge.
    pub fn new(
        panel_a: usize,
        panel_b: usize,
        axis: [f64; 3],
        rest_angle: f64,
        stiffness: f64,
    ) -> Self {
        Self {
            panel_a,
            panel_b,
            axis: vec3_normalize(axis),
            angle: 0.0,
            rest_angle,
            stiffness,
        }
    }

    /// Elastic restoring torque: τ = k * (rest - angle).
    pub fn restoring_torque(&self) -> f64 {
        self.stiffness * (self.rest_angle - self.angle)
    }

    /// Advance the hinge angle by `delta` radians.
    pub fn rotate(&mut self, delta: f64) {
        self.angle += delta;
    }
}

// ── compliant fold hinge ──────────────────────────────────────────────────────

/// Compliant fold hinge model based on a linear torsional spring.
///
/// Models the fold crease as a torsional spring connecting two rigid panels.
#[derive(Debug, Clone)]
pub struct CompliantHinge {
    /// Torsional stiffness (N·m/rad).
    pub k_torsion: f64,
    /// Damping coefficient (N·m·s/rad).
    pub c_damp: f64,
    /// Current angle (rad).
    pub angle: f64,
    /// Angular velocity (rad/s).
    pub angular_velocity: f64,
    /// Rest angle (rad).
    pub rest_angle: f64,
}

impl CompliantHinge {
    /// Create a new compliant hinge.
    pub fn new(k_torsion: f64, c_damp: f64, rest_angle: f64) -> Self {
        Self {
            k_torsion,
            c_damp,
            angle: 0.0,
            angular_velocity: 0.0,
            rest_angle,
        }
    }

    /// Applied torque from an external source (e.g., SMA actuator).
    pub fn torque_elastic(&self) -> f64 {
        -self.k_torsion * (self.angle - self.rest_angle)
    }

    /// Damping torque.
    pub fn torque_damping(&self) -> f64 {
        -self.c_damp * self.angular_velocity
    }

    /// Net torque on the hinge.
    pub fn net_torque(&self) -> f64 {
        self.torque_elastic() + self.torque_damping()
    }

    /// Integrate the hinge by one step `dt` with inertia `I` (kg·m²).
    pub fn step(&mut self, dt: f64, inertia: f64, external_torque: f64) {
        let total_torque = self.net_torque() + external_torque;
        let alpha = total_torque / (inertia + 1e-15);
        self.angular_velocity += alpha * dt;
        self.angle += self.angular_velocity * dt;
    }
}

// ── waterbomb base ────────────────────────────────────────────────────────────

/// Waterbomb base origami vertex.
///
/// The waterbomb base consists of 6 creases meeting at a central vertex:
/// alternating mountain and valley folds arranged radially.
#[derive(Debug, Clone)]
pub struct WaterbombBase {
    /// Half-size of the square sheet.
    pub half_size: f64,
    /// Fold angle θ (0 = flat, π/2 = fully folded).
    pub fold_angle: f64,
}

impl WaterbombBase {
    /// Create a waterbomb base from a square sheet of given `side`.
    pub fn new(side: f64) -> Self {
        Self {
            half_size: side * 0.5,
            fold_angle: 0.0,
        }
    }

    /// Height of the central apex above the base plane.
    pub fn apex_height(&self) -> f64 {
        self.half_size * self.fold_angle.sin()
    }

    /// Projected footprint radius at current fold angle.
    pub fn footprint_radius(&self) -> f64 {
        self.half_size * self.fold_angle.cos()
    }

    /// Set fold to specific angle in \[0, π/2\].
    pub fn set_fold_angle(&mut self, theta: f64) {
        self.fold_angle = theta.clamp(0.0, PI / 2.0);
    }

    /// Fold fraction in \[0, 1\].
    pub fn fold_fraction(&self) -> f64 {
        self.fold_angle / (PI / 2.0)
    }
}

// ── Yoshimura pattern ─────────────────────────────────────────────────────────

/// Yoshimura (diamond) buckling pattern parameters for a cylindrical shell.
#[derive(Debug, Clone)]
pub struct YoshimuraPattern {
    /// Radius of the cylinder.
    pub radius: f64,
    /// Total height of the cylinder.
    pub height: f64,
    /// Number of circumferential diamond columns.
    pub n_cols: usize,
    /// Number of axial diamond rows.
    pub n_rows: usize,
    /// Current axial compression ratio (0 = uncompressed, 1 = fully folded).
    pub compression: f64,
}

impl YoshimuraPattern {
    /// Create a Yoshimura pattern.
    pub fn new(radius: f64, height: f64, n_cols: usize, n_rows: usize) -> Self {
        Self {
            radius,
            height,
            n_cols,
            n_rows,
            compression: 0.0,
        }
    }

    /// Effective compressed height.
    pub fn compressed_height(&self) -> f64 {
        self.height * (1.0 - self.compression)
    }

    /// Diamond half-angle in the circumferential direction.
    pub fn diamond_half_angle(&self) -> f64 {
        PI / self.n_cols as f64
    }

    /// Compressive strain.
    pub fn compressive_strain(&self) -> f64 {
        self.compression
    }

    /// Apply compression increment, clamp to \[0, 1).
    pub fn compress(&mut self, delta: f64) {
        self.compression = (self.compression + delta).clamp(0.0, 0.999);
    }
}

// ── deployable structures ─────────────────────────────────────────────────────

/// State of a deployable structure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DeployState {
    /// Fully stowed/packed.
    Stowed,
    /// Partially deployed.
    Deploying,
    /// Fully deployed.
    Deployed,
}

/// An origami-inspired deployable panel structure.
///
/// Stores a set of panels and fold hinges. Deployment is driven by advancing
/// all hinge angles toward their rest angles.
#[derive(Debug, Clone)]
pub struct DeployableStructure {
    /// Rigid panels.
    pub panels: Vec<RigidPanel>,
    /// Fold hinges.
    pub hinges: Vec<FoldHinge>,
    /// Current deployment state.
    pub state: DeployState,
    /// Deployment progress in \[0, 1\].
    pub progress: f64,
}

impl DeployableStructure {
    /// Create a new deployable structure.
    pub fn new(panels: Vec<RigidPanel>, hinges: Vec<FoldHinge>) -> Self {
        Self {
            panels,
            hinges,
            state: DeployState::Stowed,
            progress: 0.0,
        }
    }

    /// Advance deployment by `step` (fraction of full deployment).
    pub fn deploy_step(&mut self, step: f64) {
        self.progress = (self.progress + step).clamp(0.0, 1.0);
        for hinge in &mut self.hinges {
            hinge.angle = hinge.rest_angle * self.progress;
        }
        self.state = if self.progress >= 1.0 - 1e-6 {
            DeployState::Deployed
        } else if self.progress <= 1e-6 {
            DeployState::Stowed
        } else {
            DeployState::Deploying
        };
    }

    /// Reset deployment to stowed state.
    pub fn reset(&mut self) {
        self.progress = 0.0;
        for hinge in &mut self.hinges {
            hinge.angle = 0.0;
        }
        self.state = DeployState::Stowed;
    }

    /// Check if all hinges have reached their rest angles (within tolerance).
    pub fn is_fully_deployed(&self) -> bool {
        self.hinges
            .iter()
            .all(|h| (h.angle - h.rest_angle).abs() < 1e-3)
    }
}

// ── bistable origami ──────────────────────────────────────────────────────────

/// Bistable snap-through origami model.
///
/// Models a bistable crease with a double-well potential:
/// V(θ) = k_snap * (θ² - θ_c²)² where θ_c is the critical angle.
#[derive(Debug, Clone)]
pub struct BistableOrigami {
    /// Double-well snap stiffness (N·m/rad³).
    pub k_snap: f64,
    /// Critical angle θ_c (rad) — half the barrier position.
    pub theta_c: f64,
    /// Current fold angle (rad).
    pub angle: f64,
    /// Angular velocity (rad/s).
    pub angular_velocity: f64,
}

impl BistableOrigami {
    /// Create a bistable origami element.
    pub fn new(k_snap: f64, theta_c: f64) -> Self {
        Self {
            k_snap,
            theta_c,
            angle: -theta_c * 0.99,
            angular_velocity: 0.0,
        }
    }

    /// Potential energy V(θ).
    pub fn potential_energy(&self) -> f64 {
        let x = self.angle * self.angle - self.theta_c * self.theta_c;
        self.k_snap * x * x
    }

    /// Restoring torque: dV/dθ with negation (force = -dV/dθ).
    pub fn restoring_torque(&self) -> f64 {
        -4.0 * self.k_snap * self.angle * (self.angle * self.angle - self.theta_c * self.theta_c)
    }

    /// Check which stable state the crease is in.
    ///
    /// Returns `1` for the positive well, `-1` for the negative well, `0` near the barrier.
    pub fn stable_state(&self) -> i32 {
        if self.angle > self.theta_c * 0.1 {
            1
        } else if self.angle < -self.theta_c * 0.1 {
            -1
        } else {
            0
        }
    }

    /// Integrate by one time step `dt` with inertia `I` and external torque.
    pub fn step(&mut self, dt: f64, inertia: f64, external_torque: f64) {
        let torque = self.restoring_torque() + external_torque;
        let alpha = torque / (inertia + 1e-15);
        self.angular_velocity += alpha * dt;
        self.angle += self.angular_velocity * dt;
    }
}

// ── folded-sheet stiffness ────────────────────────────────────────────────────

/// Compute the effective longitudinal stiffness of a Miura-ori sheet.
///
/// Based on the homogenized stiffness model from Schenk & Guest (2013).
///
/// # Arguments
/// * `t` – sheet thickness (m)
/// * `e_mat` – material Young's modulus (Pa)
/// * `a` – unit cell dimension along x (m)
/// * `b` – unit cell dimension along y (m)
/// * `gamma` – acute parallelogram angle (rad)
/// * `theta` – fold angle (rad)
///
/// Returns the effective stiffness per unit width (N/m).
pub fn miura_longitudinal_stiffness(
    t: f64,
    e_mat: f64,
    a: f64,
    b: f64,
    gamma: f64,
    theta: f64,
) -> f64 {
    // Simplified analytic approximation.
    // E_eff ∝ E·t³ / (a·b) * f(γ, θ)
    let f_geom = gamma.sin().powi(2) * theta.cos().powi(2)
        / (a * b
            * (1.0 + theta.sin().powi(2) * gamma.cos().powi(2) / gamma.sin().powi(2)).max(1e-10));
    e_mat * t.powi(3) * f_geom
}

/// Compute the effective transverse stiffness of a Miura-ori sheet (N/m).
pub fn miura_transverse_stiffness(
    t: f64,
    e_mat: f64,
    a: f64,
    b: f64,
    gamma: f64,
    theta: f64,
) -> f64 {
    let f_geom = gamma.cos().powi(2) * theta.sin().powi(2) / (a * b + 1e-15);
    e_mat * t.powi(3) * f_geom
}

// ── SMA origami actuator ──────────────────────────────────────────────────────

/// Phase of a shape memory alloy wire.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SmaPhase {
    /// Martensitic (low-temperature) phase — extended/soft.
    Martensite,
    /// Austenitic (high-temperature) phase — contracted/stiff.
    Austenite,
}

/// Shape Memory Alloy (SMA) actuator driving an origami fold hinge.
///
/// Uses a simplified two-phase model: the SMA wire contracts when heated
/// above `T_austenite` and extends when cooled below `T_martensite`.
#[derive(Debug, Clone)]
pub struct SmaActuator {
    /// Wire cross-sectional area (m²).
    pub area: f64,
    /// Wire length at martensite (rest) state (m).
    pub l_martensite: f64,
    /// Wire length at austenite (contracted) state (m).
    pub l_austenite: f64,
    /// Austenite finish temperature (K).
    pub t_af: f64,
    /// Martensite finish temperature (K).
    pub t_mf: f64,
    /// Current wire temperature (K).
    pub temperature: f64,
    /// Current phase.
    pub phase: SmaPhase,
    /// Recovery stress (Pa) when contracting.
    pub recovery_stress: f64,
}

impl SmaActuator {
    /// Create a new SMA actuator.
    pub fn new(
        area: f64,
        l_martensite: f64,
        l_austenite: f64,
        t_af: f64,
        t_mf: f64,
        recovery_stress: f64,
    ) -> Self {
        Self {
            area,
            l_martensite,
            l_austenite,
            t_af,
            t_mf,
            temperature: t_mf - 10.0,
            phase: SmaPhase::Martensite,
            recovery_stress,
        }
    }

    /// Update phase based on current temperature.
    pub fn update_phase(&mut self) {
        if self.phase == SmaPhase::Martensite && self.temperature >= self.t_af {
            self.phase = SmaPhase::Austenite;
        } else if self.phase == SmaPhase::Austenite && self.temperature <= self.t_mf {
            self.phase = SmaPhase::Martensite;
        }
    }

    /// Current wire length (interpolated by temperature).
    pub fn current_length(&self) -> f64 {
        match self.phase {
            SmaPhase::Martensite => self.l_martensite,
            SmaPhase::Austenite => {
                let frac = ((self.temperature - self.t_mf) / (self.t_af - self.t_mf + 1e-15))
                    .clamp(0.0, 1.0);
                self.l_martensite + frac * (self.l_austenite - self.l_martensite)
            }
        }
    }

    /// Actuator force (N) = stress × area (only during contraction).
    pub fn actuation_force(&self) -> f64 {
        if self.phase == SmaPhase::Austenite {
            self.recovery_stress * self.area
        } else {
            0.0
        }
    }

    /// Heat the actuator by `delta_t` Kelvin and update phase.
    pub fn heat(&mut self, delta_t: f64) {
        self.temperature += delta_t;
        self.update_phase();
    }

    /// Cool the actuator by `delta_t` Kelvin and update phase.
    pub fn cool(&mut self, delta_t: f64) {
        self.temperature -= delta_t;
        self.update_phase();
    }

    /// Contraction ratio (0 = no contraction, 1 = full contraction).
    pub fn contraction_ratio(&self) -> f64 {
        let contraction = self.l_martensite - self.current_length();
        (contraction / (self.l_martensite - self.l_austenite + 1e-15)).clamp(0.0, 1.0)
    }
}

// ── SMA-driven fold system ────────────────────────────────────────────────────

/// An origami fold hinge driven by an SMA wire actuator.
#[derive(Debug, Clone)]
pub struct SmaFoldSystem {
    /// Underlying compliant fold hinge.
    pub hinge: CompliantHinge,
    /// SMA wire actuator.
    pub actuator: SmaActuator,
    /// Mechanical advantage (moment arm) between wire contraction and fold angle (m/rad).
    pub mechanical_advantage: f64,
}

impl SmaFoldSystem {
    /// Create a new SMA-driven fold system.
    pub fn new(hinge: CompliantHinge, actuator: SmaActuator, mechanical_advantage: f64) -> Self {
        Self {
            hinge,
            actuator,
            mechanical_advantage,
        }
    }

    /// Compute the torque generated by the SMA wire on the hinge.
    pub fn sma_torque(&self) -> f64 {
        self.actuator.actuation_force() * self.mechanical_advantage
    }

    /// Integrate by one time step `dt` with panel inertia `I` (kg·m²).
    pub fn step(&mut self, dt: f64, inertia: f64) {
        let external_torque = self.sma_torque();
        self.hinge.step(dt, inertia, external_torque);
    }

    /// Activate actuator by heating `delta_t` Kelvin.
    pub fn actuate(&mut self, delta_t: f64) {
        self.actuator.heat(delta_t);
    }

    /// Deactivate actuator by cooling `delta_t` Kelvin.
    pub fn deactivate(&mut self, delta_t: f64) {
        self.actuator.cool(delta_t);
    }

    /// Current fold angle.
    pub fn fold_angle(&self) -> f64 {
        self.hinge.angle
    }
}

// ── snap-through force model ──────────────────────────────────────────────────

/// Snap-through force curve for a bistable arch or crease.
///
/// Computes the characteristic non-monotonic force-displacement relationship
/// used to characterise snap-through behaviour.
///
/// F(x) = k * x * (x² - x_c²) · (-1)  — cubic restoring force.
pub fn snap_through_force(k: f64, x: f64, x_c: f64) -> f64 {
    -k * x * (x * x - x_c * x_c)
}

/// Find the snap-through critical force magnitude (peak of the S-curve).
///
/// F_max = k * (2/3√3) * x_c³
pub fn snap_through_critical_force(k: f64, x_c: f64) -> f64 {
    k * (2.0 / (3.0 * 3.0_f64.sqrt())) * x_c.powi(3)
}

// ── vertex fold angle solver ──────────────────────────────────────────────────

/// Folding branch of a flat-foldable degree-4 vertex.
///
/// A degree-4 vertex whose sector angles satisfy `α0 + α2 = α1 + α3 = π` can
/// be folded rigidly along two different one-parameter motions that both
/// start from the flat state.  The branch selects the relation between the
/// fold angles of the opposite creases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Degree4Branch {
    /// `ρ2 = ρ0` and `ρ3 = -ρ1`, with
    /// `tan(ρ1/2) = -cos((α0+α1)/2) / cos((α0-α1)/2) · tan(ρ0/2)`.
    A,
    /// `ρ2 = -ρ0` and `ρ3 = ρ1`, with
    /// `tan(ρ1/2) = sin((α0+α1)/2) / sin((α0-α1)/2) · tan(ρ0/2)`.
    B,
}

/// Product of two 3×3 matrices (row-major).
fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for (i, out_row) in out.iter_mut().enumerate() {
        for (j, entry) in out_row.iter_mut().enumerate() {
            *entry = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

/// Closure matrix of a degree-4 origami vertex.
///
/// Walking once around the vertex, each crease `i` contributes a rotation
/// about the crease direction by its fold angle `ρ_i`, and each sector `i`
/// (between crease `i` and crease `i + 1`) a rotation about the sheet normal
/// by its sector angle `α_i`:
///
/// ```text
/// M = R_x(ρ0) R_z(α0) R_x(ρ1) R_z(α1) R_x(ρ2) R_z(α2) R_x(ρ3) R_z(α3)
/// ```
///
/// The fold angles describe a rigid folding of the vertex exactly when `M` is
/// the identity (Belcastro & Hull, 2002).
pub fn degree4_vertex_closure(alpha: [f64; 4], rho: [f64; 4]) -> [[f64; 3]; 3] {
    let mut m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for (&sector, &fold) in alpha.iter().zip(rho.iter()) {
        m = mat3_mul(m, rotation_matrix_axis_angle([1.0, 0.0, 0.0], fold));
        m = mat3_mul(m, rotation_matrix_axis_angle([0.0, 0.0, 1.0], sector));
    }
    m
}

/// Solve for the fold angles of a flat-foldable degree-4 origami vertex.
///
/// `alpha = [α0, α1, α2, α3]` are the sector angles in order around the
/// vertex (`α_i` lies between crease `i` and crease `i + 1`), `rho0` is the
/// driving fold angle of crease 0, and the returned `[ρ0, ρ1, ρ2, ρ3]` are the
/// fold angles of the creases on the selected rigid-folding `branch`
/// (see [`Degree4Branch`]); they satisfy the closure condition
/// `degree4_vertex_closure(alpha, ρ) = I`.
///
/// Flat-foldability requires every sector angle to lie in `(0, π)`, the
/// sector angles to sum to `2π` (developable vertex) and Kawasaki's condition
/// `α0 + α2 = π`.  For such a vertex the closure condition reduces to the
/// half-angle relations of Huffman (1976) documented on [`Degree4Branch`].
///
/// # Errors
///
/// Returns [`Error::InvalidParameter`] if `rho0` is not finite or outside
/// `[-π, π]`, a sector angle is outside `(0, π)`, the sector angles do not sum
/// to `2π`, or Kawasaki's condition is violated (tolerance `1e-9`).
pub fn degree4_vertex_fold_angles(
    rho0: f64,
    alpha: [f64; 4],
    branch: Degree4Branch,
) -> Result<[f64; 4]> {
    const TOLERANCE: f64 = 1e-9;
    if !rho0.is_finite() || rho0.abs() > PI {
        return Err(Error::InvalidParameter(format!(
            "degree-4 vertex: driving fold angle must be finite and within [-π, π], got {rho0}"
        )));
    }
    if let Some(bad) = alpha
        .iter()
        .find(|&&sector| !(sector.is_finite() && sector > 0.0 && sector < PI))
    {
        return Err(Error::InvalidParameter(format!(
            "degree-4 vertex: sector angles must lie in (0, π), got {bad}"
        )));
    }
    let total: f64 = alpha.iter().sum();
    if (total - 2.0 * PI).abs() > TOLERANCE {
        return Err(Error::InvalidParameter(format!(
            "degree-4 vertex: sector angles must sum to 2π (developable vertex), got {total}"
        )));
    }
    let kawasaki = alpha[0] + alpha[2];
    if (kawasaki - PI).abs() > TOLERANCE {
        return Err(Error::InvalidParameter(format!(
            "degree-4 vertex is not flat-foldable: Kawasaki's condition α0 + α2 = π is violated (α0 + α2 = {kawasaki})"
        )));
    }
    let half_tan = (0.5 * rho0).tan();
    let half_sum = 0.5 * (alpha[0] + alpha[1]);
    let half_diff = 0.5 * (alpha[0] - alpha[1]);
    match branch {
        Degree4Branch::A => {
            // |α0 - α1| < π, so cos((α0 - α1)/2) > 0.
            let rho1 = 2.0 * (-half_sum.cos() / half_diff.cos() * half_tan).atan();
            Ok([rho0, rho1, rho0, -rho1])
        }
        Degree4Branch::B => {
            // sin((α0 - α1)/2) vanishes for α0 = α1; atan2 with a non-negative
            // denominator evaluates the quotient without dividing by zero.
            let numerator = half_sum.sin() * half_tan;
            let denominator = half_diff.sin();
            let rho1 = if denominator < 0.0 {
                2.0 * (-numerator).atan2(-denominator)
            } else {
                2.0 * numerator.atan2(denominator)
            };
            Ok([rho0, rho1, -rho0, rho1])
        }
    }
}

// ── panel normal ──────────────────────────────────────────────────────────────

/// Compute the unit normal of a planar panel from its first three vertices.
pub fn panel_normal(v: &[[f64; 3]]) -> [f64; 3] {
    if v.len() < 3 {
        return [0.0, 0.0, 1.0];
    }
    let e1 = vec3_sub(v[1], v[0]);
    let e2 = vec3_sub(v[2], v[0]);
    vec3_normalize(vec3_cross(e1, e2))
}

// ── energy model for fold pattern ────────────────────────────────────────────

/// Elastic energy stored in a set of fold hinges.
pub fn total_fold_energy(hinges: &[FoldHinge]) -> f64 {
    hinges
        .iter()
        .map(|h| 0.5 * h.stiffness * (h.angle - h.rest_angle).powi(2))
        .sum()
}

/// Total elastic energy of compliant hinges.
pub fn total_compliant_hinge_energy(hinges: &[CompliantHinge]) -> f64 {
    hinges
        .iter()
        .map(|h| 0.5 * h.k_torsion * (h.angle - h.rest_angle).powi(2))
        .sum()
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// Crease vectors `[c1, c2, c3, c4]` of one vertex of the folded Miura-ori
    /// in counter-clockwise order of the flat sheet.
    fn miura_vertex_creases(unit: &MiuraOriUnit) -> [[f64; 3]; 4] {
        let d = unit.dimensions();
        [
            [d.l, 0.0, d.h],
            [d.v, d.s, 0.0],
            [-d.l, 0.0, d.h],
            [d.v, -d.s, 0.0],
        ]
    }

    // 1. vec3_add basic
    #[test]
    fn test_vec3_add() {
        let r = vec3_add([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((r[0] - 5.0).abs() < EPS);
        assert!((r[1] - 7.0).abs() < EPS);
        assert!((r[2] - 9.0).abs() < EPS);
    }

    // 2. vec3_cross orthogonality
    #[test]
    fn test_vec3_cross_orthogonal() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[0]).abs() < EPS);
        assert!((z[1]).abs() < EPS);
        assert!((z[2] - 1.0).abs() < EPS);
    }

    // 3. vec3_normalize unit length
    #[test]
    fn test_vec3_normalize_unit() {
        let v = vec3_normalize([3.0, 4.0, 0.0]);
        assert!((vec3_norm(v) - 1.0).abs() < EPS);
    }

    // 4. vec3_normalize zero vector
    #[test]
    fn test_vec3_normalize_zero() {
        let v = vec3_normalize([0.0, 0.0, 0.0]);
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }

    // 5. dihedral_angle flat panels
    #[test]
    fn test_dihedral_angle_flat() {
        // Two triangles in the XY plane sharing edge (0,0,0)-(1,0,0)
        let angle = dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        );
        // Should be π (flat)
        assert!(
            (angle - PI).abs() < 1e-6,
            "flat dihedral should be π, got {angle}"
        );
    }

    // 6. dihedral_angle perpendicular panels
    #[test]
    fn test_dihedral_angle_perpendicular() {
        let angle = dihedral_angle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.0, 1.0],
        );
        assert!(
            (angle - PI / 2.0).abs() < 1e-6,
            "perpendicular angle should be π/2, got {angle}"
        );
    }

    // 7. CreasePattern Maekawa theorem
    #[test]
    fn test_maekawa_theorem() {
        // Maekawa's theorem: |M - V| = 2 for any flat-foldable interior vertex.
        // Use 3 Mountain + 1 Valley: |3 - 1| = 2.
        let mut cp = CreasePattern::new();
        let center = cp.add_vertex(0.0, 0.0);
        let _v1 = cp.add_vertex(1.0, 0.0);
        let _v2 = cp.add_vertex(0.0, 1.0);
        let _v3 = cp.add_vertex(-1.0, 0.0);
        let _v4 = cp.add_vertex(0.0, -1.0);
        cp.add_crease(center, 1, CreaseAssignment::Mountain, PI);
        cp.add_crease(center, 2, CreaseAssignment::Mountain, PI);
        cp.add_crease(center, 3, CreaseAssignment::Mountain, PI);
        cp.add_crease(center, 4, CreaseAssignment::Valley, 0.0);
        assert!(
            cp.check_maekawa(center),
            "Maekawa theorem should hold for 3M+1V"
        );
    }

    // 8. CreasePattern count mountain/valley
    #[test]
    fn test_crease_counts() {
        let mut cp = CreasePattern::new();
        let a = cp.add_vertex(0.0, 0.0);
        let b = cp.add_vertex(1.0, 0.0);
        let c = cp.add_vertex(0.0, 1.0);
        cp.add_crease(a, b, CreaseAssignment::Mountain, PI);
        cp.add_crease(a, c, CreaseAssignment::Valley, 0.0);
        cp.add_crease(b, c, CreaseAssignment::Mountain, PI);
        assert_eq!(cp.count_mountain(), 2);
        assert_eq!(cp.count_valley(), 1);
    }

    // 9. MiuraOriUnit projected dimensions: the flat sheet has its maximal extent
    #[test]
    fn test_miura_projected_dimensions() {
        let unit = MiuraOriUnit::new(1.0, 1.0, PI / 4.0);
        // Flat (θ = 0): width = 2a, height = 2 b sin γ, no out-of-plane height.
        assert!((unit.projected_width() - 2.0).abs() < EPS);
        assert!((unit.projected_height() - 2.0 * (PI / 4.0).sin()).abs() < EPS);
        assert!(unit.height().abs() < EPS);
    }

    // 10. MiuraOriUnit fold changes dimensions: both in-plane extents shrink
    #[test]
    fn test_miura_fold_step() {
        let mut unit = MiuraOriUnit::new(1.0, 1.0, PI / 4.0);
        let (w0, h0) = (unit.projected_width(), unit.projected_height());
        unit.fold_step(PI / 4.0);
        assert!(unit.projected_width() > 0.0);
        assert!(unit.projected_width() < w0);
        assert!(unit.projected_height() > 0.0);
        assert!(unit.projected_height() < h0);
        assert!(unit.height() > 0.0);
    }

    // 11. MiuraOriSheet total dimensions scale with cell count
    #[test]
    fn test_miura_sheet_dimensions() {
        let mut sheet = MiuraOriSheet::new(3, 4, 1.0, 1.0, PI / 4.0);
        sheet.set_fold_angle(PI / 4.0);
        let w1 = sheet.total_width();
        let h1 = sheet.total_height();
        // Double the cells
        let mut sheet2 = MiuraOriSheet::new(6, 4, 1.0, 1.0, PI / 4.0);
        sheet2.set_fold_angle(PI / 4.0);
        assert!((sheet2.total_width() - 2.0 * w1).abs() < 1e-10);
        assert!((sheet2.total_height() - h1).abs() < 1e-10);
    }

    // 12. rotation_matrix_axis_angle: rotation of z around z by 90 degrees
    #[test]
    fn test_rotation_matrix_z90() {
        let r = rotation_matrix_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let x_rot = mat3_mul_vec(r, [1.0, 0.0, 0.0]);
        // Should become [0, 1, 0]
        assert!((x_rot[0]).abs() < 1e-10);
        assert!((x_rot[1] - 1.0).abs() < 1e-10);
        assert!((x_rot[2]).abs() < 1e-10);
    }

    // 13. FoldHinge restoring torque
    #[test]
    fn test_fold_hinge_restoring_torque() {
        let hinge = FoldHinge::new(0, 1, [0.0, 0.0, 1.0], PI / 2.0, 10.0);
        // angle = 0, rest = π/2 → torque = 10 * π/2
        let expected = 10.0 * PI / 2.0;
        assert!((hinge.restoring_torque() - expected).abs() < 1e-10);
    }

    // 14. CompliantHinge step changes angle
    #[test]
    fn test_compliant_hinge_step() {
        let mut hinge = CompliantHinge::new(1.0, 0.0, PI / 4.0);
        // Apply no external torque; elastic force drives it toward rest
        for _ in 0..1000 {
            hinge.step(0.01, 0.001, 0.0);
        }
        // angle should have moved toward rest_angle
        assert!(
            hinge.angle > 0.0,
            "hinge should have moved toward positive rest angle"
        );
    }

    // 15. WaterbombBase apex height
    #[test]
    fn test_waterbomb_apex_height() {
        let mut wb = WaterbombBase::new(2.0);
        wb.set_fold_angle(PI / 2.0);
        // half_size = 1.0, sin(π/2) = 1 → height = 1.0
        assert!((wb.apex_height() - 1.0).abs() < EPS);
    }

    // 16. WaterbombBase footprint radius
    #[test]
    fn test_waterbomb_footprint() {
        let mut wb = WaterbombBase::new(2.0);
        wb.set_fold_angle(PI / 4.0);
        let expected = 1.0 * (PI / 4.0).cos();
        assert!((wb.footprint_radius() - expected).abs() < EPS);
    }

    // 17. YoshimuraPattern compressed height
    #[test]
    fn test_yoshimura_compression() {
        let mut yp = YoshimuraPattern::new(0.1, 1.0, 8, 4);
        yp.compress(0.3);
        assert!((yp.compressed_height() - 0.7).abs() < EPS);
    }

    // 18. DeployableStructure deploy_step
    #[test]
    fn test_deployable_deploy_step() {
        let panels = vec![RigidPanel::new(0, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])];
        let hinges = vec![FoldHinge::new(0, 0, [0.0, 0.0, 1.0], PI / 2.0, 1.0)];
        let mut ds = DeployableStructure::new(panels, hinges);
        ds.deploy_step(0.5);
        assert_eq!(ds.state, DeployState::Deploying);
        assert!((ds.progress - 0.5).abs() < EPS);
        ds.deploy_step(0.5);
        assert_eq!(ds.state, DeployState::Deployed);
    }

    // 19. DeployableStructure reset
    #[test]
    fn test_deployable_reset() {
        let panels = vec![];
        let hinges = vec![FoldHinge::new(0, 1, [1.0, 0.0, 0.0], PI, 1.0)];
        let mut ds = DeployableStructure::new(panels, hinges);
        ds.deploy_step(1.0);
        ds.reset();
        assert_eq!(ds.state, DeployState::Stowed);
        assert!(ds.progress.abs() < EPS);
    }

    // 20. BistableOrigami potential energy wells
    #[test]
    fn test_bistable_potential_wells() {
        let bs = BistableOrigami::new(1.0, 1.0);
        // At θ = ±θ_c: V = 0
        let mut bs_pos = bs.clone();
        bs_pos.angle = 1.0;
        assert!(bs_pos.potential_energy().abs() < 1e-10);
        let mut bs_neg = bs.clone();
        bs_neg.angle = -1.0;
        assert!(bs_neg.potential_energy().abs() < 1e-10);
    }

    // 21. BistableOrigami snap-through stable states
    #[test]
    fn test_bistable_stable_states() {
        let mut bs = BistableOrigami::new(1.0, 1.0);
        bs.angle = -0.99;
        assert_eq!(bs.stable_state(), -1);
        bs.angle = 0.99;
        assert_eq!(bs.stable_state(), 1);
    }

    // 22. BistableOrigami step integrates correctly
    #[test]
    fn test_bistable_step_integration() {
        let mut bs = BistableOrigami::new(1.0, 1.0);
        bs.angle = 1.5; // past the positive well
        let initial_v = bs.angular_velocity;
        bs.step(0.001, 0.01, 0.0);
        // Force should drive it back toward θ_c = 1.0
        assert!(bs.angular_velocity != initial_v);
    }

    // 23. SmaActuator phase transition
    #[test]
    fn test_sma_phase_transition() {
        let mut sma = SmaActuator::new(1e-6, 0.1, 0.09, 350.0, 300.0, 200e6);
        assert_eq!(sma.phase, SmaPhase::Martensite);
        sma.heat(60.0); // now at 350 K = t_af
        assert_eq!(sma.phase, SmaPhase::Austenite);
    }

    // 24. SmaActuator contraction ratio
    #[test]
    fn test_sma_contraction_ratio() {
        let mut sma = SmaActuator::new(1e-6, 0.1, 0.09, 350.0, 300.0, 200e6);
        sma.heat(60.0); // fully austenitic
        let cr = sma.contraction_ratio();
        assert!(cr > 0.0 && cr <= 1.0);
    }

    // 25. SmaActuator actuation force
    #[test]
    fn test_sma_actuation_force() {
        let mut sma = SmaActuator::new(1e-6, 0.1, 0.09, 350.0, 300.0, 200e6);
        assert!(sma.actuation_force().abs() < EPS, "no force in martensite");
        sma.heat(60.0);
        assert!(sma.actuation_force() > 0.0, "force in austenite");
    }

    // 26. SmaFoldSystem torque drives fold angle
    #[test]
    fn test_sma_fold_system_activation() {
        let hinge = CompliantHinge::new(0.01, 0.001, PI / 2.0);
        let sma = SmaActuator::new(1e-5, 0.1, 0.08, 350.0, 300.0, 200e6);
        let mut sys = SmaFoldSystem::new(hinge, sma, 0.01);
        sys.actuate(60.0); // trigger austenite
        let initial_angle = sys.fold_angle();
        for _ in 0..500 {
            sys.step(0.001, 0.001);
        }
        assert!(
            sys.fold_angle() > initial_angle,
            "fold should increase under SMA actuation"
        );
    }

    // 27. miura_longitudinal_stiffness positive
    #[test]
    fn test_miura_stiffness_positive() {
        let s = miura_longitudinal_stiffness(0.001, 70e9, 0.05, 0.05, PI / 4.0, PI / 4.0);
        assert!(s > 0.0, "stiffness should be positive, got {s}");
    }

    // 28. snap_through_force zero at stable points
    #[test]
    fn test_snap_through_zero_at_equilibria() {
        let k = 1.0;
        let x_c = 1.0;
        // Force zero at x=0, x=±x_c
        assert!(snap_through_force(k, 0.0, x_c).abs() < EPS);
        assert!(snap_through_force(k, x_c, x_c).abs() < EPS);
        assert!(snap_through_force(k, -x_c, x_c).abs() < EPS);
    }

    // 29. total_fold_energy sum
    #[test]
    fn test_total_fold_energy() {
        let hinges = vec![
            FoldHinge {
                panel_a: 0,
                panel_b: 1,
                axis: [1.0, 0.0, 0.0],
                angle: 1.0,
                rest_angle: 0.0,
                stiffness: 2.0,
            },
            FoldHinge {
                panel_a: 1,
                panel_b: 2,
                axis: [0.0, 1.0, 0.0],
                angle: 0.5,
                rest_angle: 0.0,
                stiffness: 4.0,
            },
        ];
        // 0.5 * 2 * 1² + 0.5 * 4 * 0.25 = 1.0 + 0.5 = 1.5
        assert!((total_fold_energy(&hinges) - 1.5).abs() < EPS);
    }

    // 30. panel_normal correct
    #[test]
    fn test_panel_normal_xy_plane() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let n = panel_normal(&verts);
        // Normal should be (0, 0, 1) or (0, 0, -1)
        assert!(n[0].abs() < EPS);
        assert!(n[1].abs() < EPS);
        assert!((n[2].abs() - 1.0).abs() < EPS);
    }

    // 31. RigidPanel centroid
    #[test]
    fn test_rigid_panel_centroid() {
        let panel = RigidPanel::new(
            0,
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
            ],
        );
        let c = panel.centroid();
        assert!((c[0] - 1.0).abs() < EPS);
        assert!((c[1] - 1.0).abs() < EPS);
    }

    // 32. degree4_vertex_fold_angles: both branches close the vertex
    #[test]
    fn test_degree4_vertex_closure_both_branches() {
        let cases = [(1.1, 0.7), (0.7, 1.1), (1.0, 1.0), (2.2, 0.4), (0.3, 2.9)];
        for &(a0, a1) in &cases {
            let alpha = [a0, a1, PI - a0, PI - a1];
            for &rho0 in &[0.3, 1.2, 2.0, -0.8, 3.0, -PI] {
                for branch in [Degree4Branch::A, Degree4Branch::B] {
                    let rho = degree4_vertex_fold_angles(rho0, alpha, branch)
                        .expect("flat-foldable vertex");
                    assert_eq!(rho[0], rho0);
                    let m = degree4_vertex_closure(alpha, rho);
                    for (i, row) in m.iter().enumerate() {
                        for (j, &value) in row.iter().enumerate() {
                            let expected = if i == j { 1.0 } else { 0.0 };
                            assert!(
                                (value - expected).abs() < 1e-10,
                                "{branch:?} a=({a0},{a1}) rho0={rho0}: M[{i}][{j}] = {value}"
                            );
                        }
                    }
                }
            }
        }
    }

    // 32a. closure holds on a deterministic random sweep of sector angles and
    // driving fold angles, on both branches
    #[test]
    fn test_degree4_vertex_closure_random_sweep() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next_unit = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..2000 {
            let a0 = 0.05 + (PI - 0.1) * next_unit();
            let a1 = 0.05 + (PI - 0.1) * next_unit();
            let alpha = [a0, a1, PI - a0, PI - a1];
            let rho0 = (2.0 * next_unit() - 1.0) * (PI - 0.01);
            for branch in [Degree4Branch::A, Degree4Branch::B] {
                let rho = degree4_vertex_fold_angles(rho0, alpha, branch).expect("flat-foldable");
                let m = degree4_vertex_closure(alpha, rho);
                for (i, row) in m.iter().enumerate() {
                    for (j, &value) in row.iter().enumerate() {
                        let expected = if i == j { 1.0 } else { 0.0 };
                        assert!(
                            (value - expected).abs() < 1e-9,
                            "{branch:?} a=({a0},{a1}) rho0={rho0}: M[{i}][{j}] = {value}"
                        );
                    }
                }
            }
        }
    }

    // 32b. fold angles agree with an independent numerical solve of the closure
    #[test]
    fn test_degree4_vertex_matches_numerical_closure_solution() {
        let alpha = [1.1, 0.7, PI - 1.1, PI - 0.7];
        // (rho0, branch A rho1, branch B rho1) from a Newton solve of
        // R_x(ρ0) R_z(α0) ... R_x(ρ3) R_z(α3) = I.
        let reference = [
            (0.3, -0.191132, 1.074808),
            (1.2, -0.818796, 2.431568),
            (2.0, -1.558512, 2.818729),
        ];
        for &(rho0, rho1_a, rho1_b) in &reference {
            let a = degree4_vertex_fold_angles(rho0, alpha, Degree4Branch::A).unwrap();
            assert!((a[1] - rho1_a).abs() < 2e-6, "A: {} vs {rho1_a}", a[1]);
            assert!((a[2] - rho0).abs() < EPS && (a[3] + a[1]).abs() < EPS);
            let b = degree4_vertex_fold_angles(rho0, alpha, Degree4Branch::B).unwrap();
            assert!((b[1] - rho1_b).abs() < 2e-6, "B: {} vs {rho1_b}", b[1]);
            assert!((b[2] + rho0).abs() < EPS && (b[3] - b[1]).abs() < EPS);
        }
    }

    // 32c. the flat state is the common start of both branches
    #[test]
    fn test_degree4_vertex_flat_state() {
        let alpha = [1.1, 0.7, PI - 1.1, PI - 0.7];
        for branch in [Degree4Branch::A, Degree4Branch::B] {
            let rho = degree4_vertex_fold_angles(0.0, alpha, branch).unwrap();
            for r in rho {
                assert!(
                    r.abs() < EPS,
                    "{branch:?}: flat vertex must have zero folds"
                );
            }
        }
    }

    // 32d. Maekawa: three creases fold one way, one the other way
    #[test]
    fn test_degree4_vertex_maekawa_pattern() {
        // α0 + α1 < π: branch A gives signs (+, -, +, +) for ρ0 > 0.
        let alpha = [1.1, 0.7, PI - 1.1, PI - 0.7];
        let rho = degree4_vertex_fold_angles(0.9, alpha, Degree4Branch::A).unwrap();
        let positive = rho.iter().filter(|&&r| r > 0.0).count();
        let negative = rho.iter().filter(|&&r| r < 0.0).count();
        assert_eq!((positive, negative), (3, 1));
        // Branch B (α0 + α1 < π, α0 > α1): (+, +, -, +).
        let rho = degree4_vertex_fold_angles(0.9, alpha, Degree4Branch::B).unwrap();
        let positive = rho.iter().filter(|&&r| r > 0.0).count();
        let negative = rho.iter().filter(|&&r| r < 0.0).count();
        assert_eq!((positive, negative), (3, 1));
    }

    // 32e. non-flat-foldable or malformed input is rejected
    #[test]
    fn test_degree4_vertex_rejects_invalid_input() {
        // Sectors sum to π, not 2π (the former "symmetric" test vertex).
        let alpha = [PI / 4.0; 4];
        assert!(degree4_vertex_fold_angles(PI / 6.0, alpha, Degree4Branch::A).is_err());
        // Sum is 2π but α0 + α2 ≠ π: developable, not flat-foldable.
        let alpha = [1.0, 2.0, 1.5, 2.0 * PI - 4.5];
        assert!(degree4_vertex_fold_angles(0.5, alpha, Degree4Branch::B).is_err());
        // Sector angle of 0 or π (sum 2π and Kawasaki's condition still hold).
        let alpha = [0.0, 1.0, PI, PI - 1.0];
        assert!(degree4_vertex_fold_angles(0.5, alpha, Degree4Branch::A).is_err());
        let alpha = [PI, 1.0, 0.0, PI - 1.0];
        assert!(degree4_vertex_fold_angles(0.5, alpha, Degree4Branch::A).is_err());
        // Driving angle not finite or outside [-π, π].
        let alpha = [1.1, 0.7, PI - 1.1, PI - 0.7];
        for bad in [f64::NAN, f64::INFINITY, 3.5, -4.0] {
            assert!(degree4_vertex_fold_angles(bad, alpha, Degree4Branch::A).is_err());
        }
        match degree4_vertex_fold_angles(0.5, [PI / 4.0; 4], Degree4Branch::A) {
            Err(Error::InvalidParameter(message)) => {
                assert!(message.contains("2π"), "message: {message}")
            }
            other => panic!("expected InvalidParameter, got {other:?}"),
        }
    }

    // 32f. the old formula broke at ρ0 = π/2 (tan diverges); the new one does not
    #[test]
    fn test_degree4_vertex_is_regular_at_right_angle_fold() {
        let alpha = [1.1, 0.7, PI - 1.1, PI - 0.7];
        for branch in [Degree4Branch::A, Degree4Branch::B] {
            let rho = degree4_vertex_fold_angles(PI / 2.0, alpha, branch).unwrap();
            assert!(rho.iter().all(|r| r.is_finite()));
            let m = degree4_vertex_closure(alpha, rho);
            for (i, row) in m.iter().enumerate() {
                for (j, &value) in row.iter().enumerate() {
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!((value - expected).abs() < 1e-10);
                }
            }
        }
    }

    // 33. MiuraOriUnit Poisson ratio is negative (auxetic) for every fold angle
    #[test]
    fn test_miura_poisson_ratio() {
        let mut unit = MiuraOriUnit::new(1.0, 1.0, PI / 4.0);
        for k in 1..=29 {
            unit.fold_angle = 0.05 + (1.5 - 0.05) * k as f64 / 29.0;
            let nu = unit.poisson_ratio();
            assert!(nu < 0.0, "Miura-ori must be auxetic, got ν = {nu}");
        }
        // ν = -1 / (cos²θ tan²γ): at γ = π/4, θ = π/4 this is -2.
        unit.fold_angle = PI / 4.0;
        assert!((unit.poisson_ratio() + 2.0).abs() < 1e-12);
    }

    // 33b. ν = -(dS/S)/(dL/L) against central differences of the cell dimensions
    #[test]
    fn test_miura_poisson_ratio_matches_finite_differences() {
        let h = 1e-6;
        for &(a, b, gamma) in &[(1.0, 1.0, PI / 4.0), (2.0, 0.7, 1.1), (0.3, 1.9, 0.5)] {
            for &theta in &[0.1, 0.4, 0.9, 1.3] {
                let mut unit = MiuraOriUnit::new(a, b, gamma);
                unit.fold_angle = theta;
                let mut plus = unit.clone();
                plus.fold_angle = theta + h;
                let mut minus = unit.clone();
                minus.fold_angle = theta - h;
                let d_ln_s =
                    (plus.projected_height().ln() - minus.projected_height().ln()) / (2.0 * h);
                let d_ln_l =
                    (plus.projected_width().ln() - minus.projected_width().ln()) / (2.0 * h);
                let numeric = -d_ln_s / d_ln_l;
                let analytic = unit.poisson_ratio();
                assert!(
                    (numeric - analytic).abs() < 1e-5 * (1.0 + analytic.abs()),
                    "a={a} b={b} γ={gamma} θ={theta}: numeric {numeric} vs analytic {analytic}"
                );
                // ν does not depend on the edge lengths.
                let reference = MiuraOriUnit {
                    a: 1.0,
                    b: 1.0,
                    gamma,
                    fold_angle: theta,
                };
                assert!((reference.poisson_ratio() - analytic).abs() < 1e-12);
            }
        }
    }

    // 33c. degenerate states give NaN instead of a bogus finite value
    #[test]
    fn test_miura_poisson_ratio_degenerate_states() {
        let mut unit = MiuraOriUnit::new(1.0, 1.0, PI / 4.0);
        unit.fold_angle = PI / 2.0;
        assert!(unit.poisson_ratio().is_nan(), "flat-folded state");
        let degenerate = MiuraOriUnit::new(1.0, 1.0, 0.0);
        assert!(degenerate.poisson_ratio().is_nan(), "sin γ = 0");
        // Rectangular facets (γ = π/2) fold like a plain accordion: ν = 0.
        let mut rect = MiuraOriUnit::new(1.0, 1.0, PI / 2.0);
        rect.fold_angle = 0.8;
        assert!(rect.poisson_ratio().abs() < 1e-12);
    }

    // 33d. dimensions at the flat and flat-folded limits
    #[test]
    fn test_miura_dimension_limits() {
        let (a, b, gamma) = (1.3, 0.8, 1.0);
        let mut unit = MiuraOriUnit::new(a, b, gamma);
        let flat = unit.dimensions();
        assert!((flat.l - a).abs() < EPS);
        assert!((flat.s - b * gamma.sin()).abs() < EPS);
        assert!((flat.v - b * gamma.cos()).abs() < EPS);
        assert!(flat.h.abs() < EPS);
        unit.fold_angle = PI / 2.0;
        let folded = unit.dimensions();
        assert!((folded.l - a * gamma.cos()).abs() < EPS);
        assert!(folded.s.abs() < EPS);
        assert!((folded.h - a * gamma.sin()).abs() < EPS);
        assert!((folded.v - b).abs() < EPS);
        // fold_step saturates at the flat-folded state.
        unit.fold_angle = 0.0;
        unit.fold_step(10.0);
        assert!((unit.fold_angle - PI / 2.0).abs() < EPS);
        unit.fold_step(-10.0);
        assert!(unit.fold_angle.abs() < EPS);
    }

    // 33e. the folded cell is a rigid folding of the flat vertex with facet tilt θ
    #[test]
    fn test_miura_vertex_is_rigid_and_facets_tilt_by_theta() {
        let (a, b, gamma) = (1.3, 0.8, 1.0);
        for &theta in &[0.0, 0.2, 0.7, 1.2, 1.5] {
            let mut unit = MiuraOriUnit::new(a, b, gamma);
            unit.fold_angle = theta;
            let creases = miura_vertex_creases(&unit);
            let lengths = [a, b, a, b];
            for (crease, &length) in creases.iter().zip(lengths.iter()) {
                assert!((vec3_norm(*crease) - length).abs() < 1e-12);
            }
            // Sector angles of the flat sheet are preserved: γ, π-γ, π-γ, γ.
            let sectors = [gamma, PI - gamma, PI - gamma, gamma];
            for i in 0..4 {
                let u = vec3_normalize(creases[i]);
                let v = vec3_normalize(creases[(i + 1) % 4]);
                let angle = vec3_dot(u, v).clamp(-1.0, 1.0).acos();
                assert!((angle - sectors[i]).abs() < 1e-10, "θ={theta} sector {i}");
                // Every facet is inclined by θ against the reference plane.
                let normal = vec3_normalize(vec3_cross(creases[i], creases[(i + 1) % 4]));
                assert!(
                    (normal[2] - theta.cos()).abs() < 1e-10,
                    "θ={theta} facet {i}"
                );
            }
        }
    }

    // 33f. closed-form branch B reproduces the dihedral angles of the folded Miura vertex
    #[test]
    fn test_degree4_branch_b_matches_miura_geometry() {
        let (a, b, gamma) = (1.3, 0.8, 1.0);
        let alpha = [gamma, PI - gamma, PI - gamma, gamma];
        for &theta in &[0.1, 0.5, 0.9, 1.3] {
            let mut unit = MiuraOriUnit::new(a, b, gamma);
            unit.fold_angle = theta;
            let creases = miura_vertex_creases(&unit);
            // Oriented facet normals: facet i lies between crease i and crease i+1.
            let normals: Vec<[f64; 3]> = (0..4)
                .map(|i| vec3_normalize(vec3_cross(creases[i], creases[(i + 1) % 4])))
                .collect();
            // |fold angle| of crease i = angle between the facets on either side.
            let geometric: Vec<f64> = (0..4)
                .map(|i| {
                    let before = normals[(i + 3) % 4];
                    vec3_dot(before, normals[i]).clamp(-1.0, 1.0).acos()
                })
                .collect();
            // The facets on either side of a zigzag crease are tilted by +θ and
            // -θ against the reference plane, so the zigzag creases fold by 2θ.
            assert!((geometric[1] - 2.0 * theta).abs() < 1e-9);
            assert!((geometric[3] - 2.0 * theta).abs() < 1e-9);
            let solved = degree4_vertex_fold_angles(geometric[0], alpha, Degree4Branch::B).unwrap();
            for i in 0..4 {
                assert!(
                    (solved[i].abs() - geometric[i]).abs() < 1e-9,
                    "θ={theta} crease {i}: solved {} vs geometric {}",
                    solved[i].abs(),
                    geometric[i]
                );
            }
        }
    }

    // 34. total_compliant_hinge_energy
    #[test]
    fn test_total_compliant_hinge_energy() {
        let hinges = vec![
            CompliantHinge {
                k_torsion: 10.0,
                c_damp: 0.1,
                angle: 1.0,
                angular_velocity: 0.0,
                rest_angle: 0.0,
            },
            CompliantHinge {
                k_torsion: 5.0,
                c_damp: 0.1,
                angle: 2.0,
                angular_velocity: 0.0,
                rest_angle: 0.0,
            },
        ];
        // 0.5*10*1 + 0.5*5*4 = 5 + 10 = 15
        assert!((total_compliant_hinge_energy(&hinges) - 15.0).abs() < EPS);
    }

    // 35. snap_through_critical_force positive
    #[test]
    fn test_snap_through_critical_force() {
        let f = snap_through_critical_force(1.0, 1.0);
        assert!(f > 0.0, "critical force should be positive, got {f}");
    }
}
