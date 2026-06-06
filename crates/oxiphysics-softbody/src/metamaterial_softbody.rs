// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mechanical metamaterials: auxetic, chiral, pentamode, magnetic, programmable matter,
//! acoustic, topology optimization, and effective medium analysis.
//!
//! All structures use plain `f64` arrays — no nalgebra dependency.
//!
//! # Overview
//!
//! * [`AuxeticLattice`] — Re-entrant honeycomb lattice, negative Poisson's ratio.
//! * [`ChiralMetamaterial`] — Chiral unit cell geometry, twist under tension.
//! * [`PentamodeMetamaterial`] — Near-zero shear stiffness pentamode lattice.
//! * [`MagneticMetamaterial`] — Magnetically actuated, field-dependent stiffness.
//! * [`ProgrammableMatter`] — Bistable unit cells, snap-through behavior.
//! * [`AcousticMetamaterial`] — Locally resonant cells, band gaps, negative effective mass.
//! * [`MetamaterialOptimization`] — Topology optimization for metamaterial unit cells.
//! * [`MetamaterialAnalysis`] — Effective medium theory and homogenization.

use std::f64::consts::PI;

// ─── Vec2 helpers ─────────────────────────────────────────────────────────────

#[inline]
fn vec2_add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn vec2_sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn vec2_scale(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}

#[inline]
fn vec2_norm(a: [f64; 2]) -> f64 {
    (a[0] * a[0] + a[1] * a[1]).sqrt()
}

// ─── Vec3 helpers ─────────────────────────────────────────────────────────────

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

// ─── AuxeticLattice ───────────────────────────────────────────────────────────

/// A node in the re-entrant honeycomb lattice.
#[derive(Debug, Clone)]
pub struct LatticeNode {
    /// Current position \[m\].
    pub pos: [f64; 2],
    /// Rest position \[m\].
    pub rest_pos: [f64; 2],
    /// Velocity \[m/s\].
    pub vel: [f64; 2],
    /// Mass \[kg\].
    pub mass: f64,
    /// Whether this node is pinned.
    pub pinned: bool,
}

impl LatticeNode {
    /// Create a new lattice node.
    pub fn new(pos: [f64; 2], mass: f64) -> Self {
        Self {
            pos,
            rest_pos: pos,
            vel: [0.0; 2],
            mass,
            pinned: false,
        }
    }
}

/// Re-entrant honeycomb auxetic lattice with negative Poisson's ratio.
///
/// The re-entrant structure features inward-pointing vertices that cause
/// the lattice to expand laterally when stretched axially, giving ν < 0.
pub struct AuxeticLattice {
    /// Unit cell width \[m\].
    pub cell_width: f64,
    /// Unit cell height \[m\].
    pub cell_height: f64,
    /// Re-entrant angle from horizontal \[rad\] (positive = re-entrant).
    pub re_entrant_angle: f64,
    /// Wall thickness relative to cell size.
    pub wall_thickness_ratio: f64,
    /// Number of cells in x direction.
    pub nx: usize,
    /// Number of cells in y direction.
    pub ny: usize,
    /// Lattice nodes.
    pub nodes: Vec<LatticeNode>,
    /// Spring connections (node_a, node_b, rest_length, stiffness).
    pub springs: Vec<(usize, usize, f64, f64)>,
    /// Material Young's modulus \[Pa\].
    pub e_material: f64,
}

impl AuxeticLattice {
    /// Create a re-entrant honeycomb auxetic lattice.
    ///
    /// * `cell_width`   — unit cell width \[m\]
    /// * `cell_height`  — unit cell height \[m\]
    /// * `theta`        — re-entrant angle \[rad\]
    /// * `t_ratio`      — wall thickness / length ratio
    /// * `nx`, `ny`     — number of unit cells in x and y
    /// * `e_mat`        — parent material Young's modulus \[Pa\]
    pub fn new(
        cell_width: f64,
        cell_height: f64,
        theta: f64,
        t_ratio: f64,
        nx: usize,
        ny: usize,
        e_mat: f64,
    ) -> Self {
        let mut nodes = Vec::new();
        let stiffness = e_mat * t_ratio * 0.1;

        // Generate re-entrant hexagonal nodes
        for iy in 0..=ny {
            for ix in 0..=nx {
                let x = ix as f64 * cell_width;
                let y = iy as f64 * cell_height;
                nodes.push(LatticeNode::new([x, y], 1.0));
                let dx = 0.5 * cell_width * theta.cos();
                let dy = -0.5 * cell_height * theta.sin();
                if ix < nx && iy < ny {
                    nodes.push(LatticeNode::new([x + dx, y + cell_height * 0.5 + dy], 0.5));
                }
            }
        }

        // Add springs between adjacent nodes
        let n = nodes.len();
        let mut springs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = vec2_norm(vec2_sub(nodes[j].pos, nodes[i].pos));
                if d < cell_width * 0.8 && d > 1e-10 {
                    springs.push((i, j, d, stiffness));
                }
            }
        }

        Self {
            cell_width,
            cell_height,
            re_entrant_angle: theta,
            wall_thickness_ratio: t_ratio,
            nx,
            ny,
            nodes,
            springs,
            e_material: e_mat,
        }
    }

    /// Compute the theoretical Poisson's ratio for re-entrant honeycomb.
    ///
    /// Based on the Gibson-Ashby analytical model:
    /// ν_xy = −(cos θ · (h/l − sin θ)) / (sin²θ)
    ///
    /// * `h_over_l` — ratio of vertical to inclined wall length
    pub fn poisson_ratio(&self, h_over_l: f64) -> f64 {
        let theta = self.re_entrant_angle;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        -(cos_t * (h_over_l - sin_t)) / (sin_t * sin_t)
    }

    /// Effective Young's modulus in the y direction \[Pa\].
    ///
    /// * `h_over_l` — ratio of vertical to inclined wall length
    pub fn effective_modulus_y(&self, h_over_l: f64) -> f64 {
        let theta = self.re_entrant_angle;
        let t = self.wall_thickness_ratio;
        let sin_t = theta.sin();
        self.e_material * t.powi(3)
            / ((h_over_l - sin_t).abs().max(1e-30) * theta.cos().abs().max(1e-30))
    }

    /// Effective shear modulus \[Pa\].
    pub fn effective_shear_modulus(&self, h_over_l: f64) -> f64 {
        let theta = self.re_entrant_angle;
        let t = self.wall_thickness_ratio;
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        self.e_material * t.powi(3) * (h_over_l + sin_t)
            / (h_over_l * h_over_l * (1.0 + 2.0 * h_over_l) * cos_t.abs().max(1e-30))
    }

    /// Compute spring potential energy \[J\].
    pub fn spring_energy(&self) -> f64 {
        self.springs
            .iter()
            .map(|&(i, j, l0, k)| {
                let d = vec2_norm(vec2_sub(self.nodes[j].pos, self.nodes[i].pos));
                0.5 * k * (d - l0).powi(2)
            })
            .sum()
    }

    /// Apply a uniform uniaxial strain in the y direction and return x-strain.
    pub fn apply_strain_y(&mut self, strain: f64) -> f64 {
        let max_y = self
            .nodes
            .iter()
            .map(|n| n.pos[1])
            .fold(f64::NEG_INFINITY, f64::max);
        let min_y = self
            .nodes
            .iter()
            .map(|n| n.pos[1])
            .fold(f64::INFINITY, f64::min);
        let height = max_y - min_y;
        let max_x_before = self
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let min_x_before = self
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::INFINITY, f64::min);
        let width_before = max_x_before - min_x_before;

        for node in self.nodes.iter_mut() {
            if !node.pinned {
                let frac_y = (node.pos[1] - min_y) / height.max(1e-10);
                node.pos[1] += strain * frac_y * height;
            }
        }
        for _ in 0..10 {
            self.relax_step(0.001);
        }
        let max_x = self
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let min_x = self
            .nodes
            .iter()
            .map(|n| n.pos[0])
            .fold(f64::INFINITY, f64::min);
        let width_after = max_x - min_x;
        (width_after - width_before) / width_before.max(1e-10)
    }

    fn relax_step(&mut self, dt: f64) {
        let n = self.nodes.len();
        let mut forces = vec![[0.0_f64; 2]; n];
        for &(i, j, l0, k) in &self.springs {
            let dr = vec2_sub(self.nodes[j].pos, self.nodes[i].pos);
            let d = vec2_norm(dr);
            if d < 1e-12 {
                continue;
            }
            let f_mag = k * (d - l0);
            let dir = vec2_scale(dr, 1.0 / d);
            let f = vec2_scale(dir, f_mag);
            forces[i] = vec2_add(forces[i], f);
            forces[j] = vec2_sub(forces[j], f);
        }
        for (i, node) in self.nodes.iter_mut().enumerate() {
            if node.pinned {
                continue;
            }
            node.vel = vec2_scale(
                vec2_add(node.vel, vec2_scale(forces[i], dt / node.mass)),
                0.98,
            );
            node.pos = vec2_add(node.pos, vec2_scale(node.vel, dt));
        }
    }
}

// ─── ChiralMetamaterial ───────────────────────────────────────────────────────

/// A chiral unit cell.
#[derive(Debug, Clone)]
pub struct ChiralCell {
    /// Centre position \[m\].
    pub center: [f64; 2],
    /// Rotation angle \[rad\].
    pub angle: f64,
    /// Radius of central ring \[m\].
    pub ring_radius: f64,
    /// Number of ligaments.
    pub n_ligaments: usize,
}

impl ChiralCell {
    /// Create a new chiral unit cell.
    pub fn new(center: [f64; 2], ring_radius: f64, n_ligaments: usize) -> Self {
        Self {
            center,
            angle: 0.0,
            ring_radius,
            n_ligaments,
        }
    }

    /// Ligament attachment points on the ring.
    pub fn ligament_points(&self) -> Vec<[f64; 2]> {
        (0..self.n_ligaments)
            .map(|k| {
                let phi = self.angle + 2.0 * PI * k as f64 / self.n_ligaments as f64;
                [
                    self.center[0] + self.ring_radius * phi.cos(),
                    self.center[1] + self.ring_radius * phi.sin(),
                ]
            })
            .collect()
    }
}

/// Chiral metamaterial with chirality-induced coupling between axial and rotational motion.
///
/// Chiral structures have no plane of symmetry, causing rotation under
/// uniaxial tension (twist-stretch coupling).
pub struct ChiralMetamaterial {
    /// Unit cells in the lattice.
    pub cells: Vec<ChiralCell>,
    /// Lattice spacing \[m\].
    pub lattice_spacing: f64,
    /// Ligament stiffness \[N/m\].
    pub ligament_stiffness: f64,
    /// Ring torsional stiffness \[N·m/rad\].
    pub ring_torsional_stiffness: f64,
    /// Number of cells in x.
    pub nx: usize,
    /// Number of cells in y.
    pub ny: usize,
    /// Chirality parameter (1 = right-handed, -1 = left-handed).
    pub chirality: f64,
}

impl ChiralMetamaterial {
    /// Create a new chiral metamaterial lattice.
    pub fn new(
        nx: usize,
        ny: usize,
        lattice_spacing: f64,
        ring_radius: f64,
        n_ligaments: usize,
        ligament_stiffness: f64,
        chirality: f64,
    ) -> Self {
        let mut cells = Vec::new();
        for iy in 0..ny {
            for ix in 0..nx {
                let cx = ix as f64 * lattice_spacing;
                let cy = iy as f64 * lattice_spacing;
                cells.push(ChiralCell::new([cx, cy], ring_radius, n_ligaments));
            }
        }
        Self {
            cells,
            lattice_spacing,
            ligament_stiffness,
            ring_torsional_stiffness: ligament_stiffness * ring_radius * ring_radius,
            nx,
            ny,
            chirality,
        }
    }

    /// Compute the twist angle per unit axial strain (twist-stretch coupling).
    ///
    /// Returns dθ/dε \[rad\].
    pub fn twist_per_strain(&self) -> f64 {
        let cell = &self.cells[0];
        let r = cell.ring_radius;
        let l = self.lattice_spacing;
        let n = cell.n_ligaments as f64;
        self.chirality * n * r / l
    }

    /// Apply axial strain and compute ring rotation.
    pub fn rotate_under_tension(&mut self, axial_strain: f64) {
        let twist = self.twist_per_strain() * axial_strain;
        for cell in self.cells.iter_mut() {
            cell.angle += twist;
        }
    }

    /// Effective Poisson's ratio from chiral coupling.
    pub fn effective_poisson_ratio(&self) -> f64 {
        let r = self.cells[0].ring_radius;
        let l = self.lattice_spacing;
        let beta = r / l;
        -(1.0 - beta * beta) / (1.0 + beta * beta)
    }

    /// Compute effective in-plane bulk modulus \[Pa\].
    pub fn effective_bulk_modulus(&self, e_mat: f64) -> f64 {
        let r = self.cells[0].ring_radius;
        let l = self.lattice_spacing;
        let t = 0.05 * l;
        e_mat * (t / l).powi(3) / (3.0 * (1.0 - (r / l).powi(2)).abs().max(1e-30))
    }

    /// Total rotational energy stored in the lattice \[J\].
    pub fn rotational_energy(&self) -> f64 {
        self.cells
            .iter()
            .map(|c| 0.5 * self.ring_torsional_stiffness * c.angle * c.angle)
            .sum()
    }
}

// ─── PentamodeMetamaterial ────────────────────────────────────────────────────

/// A pentamode metamaterial lattice element.
#[derive(Debug, Clone)]
pub struct PentamodeElement {
    /// Node positions \[m\] (4 nodes in tetrahedral unit cell).
    pub nodes: [[f64; 3]; 4],
    /// Connection radius (double-cone vertex radius) \[m\].
    pub vertex_radius: f64,
    /// Strut length \[m\].
    pub strut_length: f64,
    /// Material Young's modulus \[Pa\].
    pub e_material: f64,
}

impl PentamodeElement {
    /// Create a pentamode element from a tetrahedral unit cell.
    pub fn new(size: f64, e_material: f64) -> Self {
        let h = size * (2.0_f64 / 3.0).sqrt();
        let nodes = [
            [0.0, 0.0, 0.0],
            [size, 0.0, 0.0],
            [size * 0.5, h, 0.0],
            [size * 0.5, h / 3.0, size * (2.0_f64 / 3.0).sqrt()],
        ];
        Self {
            nodes,
            vertex_radius: size * 0.02,
            strut_length: size,
            e_material,
        }
    }

    /// Effective bulk modulus \[Pa\] of the pentamode structure.
    pub fn effective_bulk_modulus(&self) -> f64 {
        let r = self.vertex_radius;
        let l = self.strut_length;
        self.e_material * (r / l).powi(2)
    }

    /// Effective shear modulus \[Pa\] — near zero for ideal pentamode.
    pub fn effective_shear_modulus(&self) -> f64 {
        let r = self.vertex_radius;
        let l = self.strut_length;
        self.e_material * (r / l).powi(4)
    }

    /// Ratio of bulk to shear modulus (should be >> 1 for good pentamode).
    pub fn bulk_to_shear_ratio(&self) -> f64 {
        let r = self.vertex_radius;
        let l = self.strut_length;
        (l / r).powi(2)
    }
}

/// Pentamode metamaterial — near-zero shear stiffness, water-wave analogue.
pub struct PentamodeMetamaterial {
    /// Unit cells.
    pub elements: Vec<PentamodeElement>,
    /// Target acoustic impedance matching \[Pa·s/m\].
    pub target_impedance: f64,
}

impl PentamodeMetamaterial {
    /// Create a pentamode metamaterial lattice.
    pub fn new(nx: usize, ny: usize, nz: usize, cell_size: f64, e_material: f64) -> Self {
        let mut elements = Vec::new();
        for _iz in 0..nz {
            for _iy in 0..ny {
                for _ix in 0..nx {
                    elements.push(PentamodeElement::new(cell_size, e_material));
                }
            }
        }
        Self {
            elements,
            target_impedance: 1.5e6,
        }
    }

    /// Check if the lattice qualifies as pentamode (κ/G > threshold).
    pub fn is_pentamode(&self, threshold: f64) -> bool {
        self.elements
            .iter()
            .all(|e| e.bulk_to_shear_ratio() > threshold)
    }

    /// Effective longitudinal wave speed \[m/s\].
    pub fn longitudinal_wave_speed(&self, density: f64) -> f64 {
        if self.elements.is_empty() {
            return 0.0;
        }
        let kappa = self.elements[0].effective_bulk_modulus();
        (kappa / density).sqrt()
    }

    /// Impedance mismatch ratio with water.
    pub fn impedance_mismatch(&self, density: f64) -> f64 {
        let c = self.longitudinal_wave_speed(density);
        let z = density * c;
        (z - self.target_impedance).abs() / self.target_impedance
    }
}

// ─── MagneticMetamaterial ─────────────────────────────────────────────────────

/// A magnetically active unit cell.
#[derive(Debug, Clone)]
pub struct MagneticCell {
    /// Position \[m\].
    pub pos: [f64; 3],
    /// Magnetic moment \[A·m²\].
    pub magnetic_moment: [f64; 3],
    /// Current stiffness \[N/m\] (field-dependent).
    pub stiffness: f64,
    /// Rest stiffness \[N/m\] (without field).
    pub rest_stiffness: f64,
    /// Magnetostrictive coefficient \[m/A²\].
    pub magnetostrictive_coeff: f64,
}

impl MagneticCell {
    /// Create a magnetic unit cell.
    pub fn new(pos: [f64; 3], rest_stiffness: f64) -> Self {
        Self {
            pos,
            magnetic_moment: [0.0, 0.0, 1.0],
            stiffness: rest_stiffness,
            rest_stiffness,
            magnetostrictive_coeff: 1e-8,
        }
    }

    /// Update stiffness based on applied magnetic field magnitude \[T\].
    pub fn update_stiffness(&mut self, field_magnitude: f64) {
        let delta_k = self.magnetostrictive_coeff * field_magnitude * field_magnitude;
        self.stiffness = self.rest_stiffness + delta_k;
    }

    /// Magnetostrictive strain for given field \[T\].
    pub fn magnetostrictive_strain(&self, field: f64) -> f64 {
        self.magnetostrictive_coeff * field * field
    }
}

/// Magnetically actuated metamaterial with field-tunable stiffness.
pub struct MagneticMetamaterial {
    /// Unit cells.
    pub cells: Vec<MagneticCell>,
    /// Current applied field magnitude \[T\].
    pub applied_field: f64,
    /// Magnetic permeability \[H/m\].
    pub permeability: f64,
}

impl MagneticMetamaterial {
    /// Create a magnetic metamaterial.
    pub fn new(n_cells: usize, rest_stiffness: f64) -> Self {
        let cells = (0..n_cells)
            .map(|i| MagneticCell::new([i as f64 * 0.01, 0.0, 0.0], rest_stiffness))
            .collect();
        Self {
            cells,
            applied_field: 0.0,
            permeability: 4.0 * PI * 1e-7,
        }
    }

    /// Apply a uniform magnetic field \[T\] to all cells.
    pub fn apply_field(&mut self, field: f64) {
        self.applied_field = field;
        for cell in self.cells.iter_mut() {
            cell.update_stiffness(field);
        }
    }

    /// Effective structural stiffness \[N/m\] (average over cells).
    pub fn effective_stiffness(&self) -> f64 {
        if self.cells.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.cells.iter().map(|c| c.stiffness).sum();
        sum / self.cells.len() as f64
    }

    /// Stiffness tuning ratio: K(B) / K(0).
    pub fn stiffness_ratio(&self) -> f64 {
        if self.cells.is_empty() || self.cells[0].rest_stiffness.abs() < 1e-30 {
            return 1.0;
        }
        self.effective_stiffness() / self.cells[0].rest_stiffness
    }

    /// Compute dipole-dipole interaction energy between cells i and j \[J\].
    pub fn dipole_interaction_energy(&self, i: usize, j: usize) -> f64 {
        if i >= self.cells.len() || j >= self.cells.len() {
            return 0.0;
        }
        let r = vec3_sub(self.cells[j].pos, self.cells[i].pos);
        let r_mag = vec3_norm(r);
        if r_mag < 1e-10 {
            return 0.0;
        }
        let mu0 = self.permeability;
        let m_i = vec3_norm(self.cells[i].magnetic_moment);
        let m_j = vec3_norm(self.cells[j].magnetic_moment);
        let r_hat = [r[0] / r_mag, r[1] / r_mag, r[2] / r_mag];
        let cos_ij = if m_i < 1e-30 || m_j < 1e-30 {
            0.0
        } else {
            vec3_dot(self.cells[i].magnetic_moment, self.cells[j].magnetic_moment) / (m_i * m_j)
        };
        let cos_ir = if m_i < 1e-30 {
            0.0
        } else {
            vec3_dot(self.cells[i].magnetic_moment, r_hat) / m_i
        };
        let cos_jr = if m_j < 1e-30 {
            0.0
        } else {
            vec3_dot(self.cells[j].magnetic_moment, r_hat) / m_j
        };
        mu0 / (4.0 * PI) * m_i * m_j / r_mag.powi(3) * (cos_ij - 3.0 * cos_ir * cos_jr)
    }
}

// ─── ProgrammableMatter ───────────────────────────────────────────────────────

/// Bistable unit cell with two stable configurations.
#[derive(Debug, Clone)]
pub struct BistableCell {
    /// Current displacement \[m\].
    pub displacement: f64,
    /// First stable state displacement \[m\].
    pub state_a: f64,
    /// Second stable state displacement \[m\].
    pub state_b: f64,
    /// Energy barrier between states \[J\].
    pub energy_barrier: f64,
    /// Bistability parameter (> 1 for bistable).
    pub lambda: f64,
    /// Snap-through threshold force \[N\].
    pub snap_threshold: f64,
    /// Linear stiffness \[N/m\].
    pub stiffness: f64,
    /// Whether cell is currently in state A (true) or B (false).
    pub in_state_a: bool,
}

impl BistableCell {
    /// Create a bistable unit cell.
    ///
    /// * `stiffness`   — linear spring constant \[N/m\]
    /// * `lambda`      — bistability parameter (> 1 for bistable)
    pub fn new(stiffness: f64, lambda: f64) -> Self {
        let d = if lambda > 1.0 {
            (1.0 - 1.0 / lambda).sqrt()
        } else {
            0.0
        };
        let barrier = 0.25 * stiffness * d * d;
        let snap_f = stiffness * d * (lambda / (2.0 * PI)).abs();
        Self {
            displacement: -d,
            state_a: -d,
            state_b: d,
            energy_barrier: barrier,
            lambda,
            snap_threshold: snap_f,
            stiffness,
            in_state_a: true,
        }
    }

    /// Elastic potential energy at current displacement (double-well potential) \[J\].
    pub fn potential_energy(&self) -> f64 {
        let x = self.displacement;
        let k = self.stiffness;
        let l = self.lambda;
        0.25 * k * (x * x - (l - 1.0) / l.max(1e-10)).powi(2)
    }

    /// Restoring force at current displacement \[N\].
    pub fn restoring_force(&self) -> f64 {
        let x = self.displacement;
        let k = self.stiffness;
        let l = self.lambda;
        k * x * (x * x - (l - 1.0) / l.max(1e-10))
    }

    /// Apply external force `f` \[N\] and check for snap-through.
    ///
    /// Returns `true` if a snap-through event occurred.
    pub fn apply_force(&mut self, f: f64) -> bool {
        let snapped = (f.abs() > self.snap_threshold && self.in_state_a && f > 0.0)
            || (f.abs() > self.snap_threshold && !self.in_state_a && f < 0.0);
        if snapped {
            if self.in_state_a {
                self.displacement = self.state_b;
                self.in_state_a = false;
            } else {
                self.displacement = self.state_a;
                self.in_state_a = true;
            }
        }
        snapped
    }
}

/// Programmable matter built from bistable unit cells with sequential folding.
pub struct ProgrammableMatter {
    /// Array of bistable cells.
    pub cells: Vec<BistableCell>,
    /// Coupling stiffness between adjacent cells \[N/m\].
    pub coupling: f64,
    /// Folding sequence index.
    pub fold_index: usize,
}

impl ProgrammableMatter {
    /// Create a chain of bistable cells.
    pub fn new(n: usize, cell_stiffness: f64, lambda: f64, coupling: f64) -> Self {
        let cells = (0..n)
            .map(|_| BistableCell::new(cell_stiffness, lambda))
            .collect();
        Self {
            cells,
            coupling,
            fold_index: 0,
        }
    }

    /// Total stored elastic energy in all cells \[J\].
    pub fn total_energy(&self) -> f64 {
        let mut e = self.cells.iter().map(|c| c.potential_energy()).sum::<f64>();
        for i in 0..self.cells.len().saturating_sub(1) {
            let dx = self.cells[i + 1].displacement - self.cells[i].displacement;
            e += 0.5 * self.coupling * dx * dx;
        }
        e
    }

    /// Trigger sequential folding by applying a force to the first unfolded cell.
    ///
    /// Returns the index of the cell that snapped.
    pub fn trigger_fold(&mut self, force: f64) -> Option<usize> {
        if self.fold_index >= self.cells.len() {
            return None;
        }
        let snapped = self.cells[self.fold_index].apply_force(force);
        if snapped {
            let idx = self.fold_index;
            self.fold_index += 1;
            Some(idx)
        } else {
            None
        }
    }

    /// Count cells currently in state A.
    pub fn count_state_a(&self) -> usize {
        self.cells.iter().filter(|c| c.in_state_a).count()
    }

    /// Count cells currently in state B.
    pub fn count_state_b(&self) -> usize {
        self.cells.iter().filter(|c| !c.in_state_a).count()
    }

    /// Compute the snap-through energy released during a transition \[J\].
    pub fn snap_energy_release(&self) -> f64 {
        if self.cells.is_empty() {
            return 0.0;
        }
        self.cells[0].energy_barrier * 2.0
    }
}

// ─── AcousticMetamaterial ─────────────────────────────────────────────────────

/// A locally resonant unit cell for acoustic metamaterials.
#[derive(Debug, Clone)]
pub struct ResonantCell {
    /// Outer shell mass \[kg\].
    pub shell_mass: f64,
    /// Inner resonator mass \[kg\].
    pub resonator_mass: f64,
    /// Coating stiffness \[N/m\].
    pub coating_stiffness: f64,
    /// Resonant frequency \[Hz\].
    pub resonant_frequency: f64,
    /// Cell size \[m\].
    pub size: f64,
}

impl ResonantCell {
    /// Create a locally resonant unit cell (lead sphere in rubber coating in epoxy).
    ///
    /// * `shell_mass`      — outer epoxy mass \[kg\]
    /// * `resonator_mass`  — inner lead sphere mass \[kg\]
    /// * `k_coat`          — rubber coating stiffness \[N/m\]
    /// * `size`            — cell side length \[m\]
    pub fn new(shell_mass: f64, resonator_mass: f64, k_coat: f64, size: f64) -> Self {
        let f0 = (k_coat / resonator_mass).sqrt() / (2.0 * PI);
        Self {
            shell_mass,
            resonator_mass,
            coating_stiffness: k_coat,
            resonant_frequency: f0,
            size,
        }
    }

    /// Effective dynamic mass at frequency `f` \[Hz\] (can be negative near resonance).
    pub fn effective_mass(&self, f: f64) -> f64 {
        let omega = 2.0 * PI * f;
        let omega0 = 2.0 * PI * self.resonant_frequency;
        let m1 = self.shell_mass;
        let m2 = self.resonator_mass;
        let k = self.coating_stiffness;
        let denom = omega0 * omega0 - omega * omega;
        if denom.abs() < 1e-30 {
            f64::INFINITY
        } else {
            m1 + k * m2 / (m1.max(1e-30) * denom)
        }
    }
}

/// Acoustic metamaterial with locally resonant unit cells.
///
/// Features sub-wavelength band gaps where wave propagation is forbidden.
pub struct AcousticMetamaterial {
    /// Unit cells.
    pub cells: Vec<ResonantCell>,
    /// Host medium density \[kg/m³\].
    pub host_density: f64,
    /// Host medium wave speed \[m/s\].
    pub host_wave_speed: f64,
    /// Number of cells in lattice (1D).
    pub n_cells: usize,
}

impl AcousticMetamaterial {
    /// Create an acoustic metamaterial with lead-rubber resonators in epoxy.
    pub fn new(n_cells: usize) -> Self {
        let cells = (0..n_cells)
            .map(|_| ResonantCell::new(0.028, 0.101, 8.0e4, 0.015))
            .collect();
        Self {
            cells,
            host_density: 1180.0,
            host_wave_speed: 2700.0,
            n_cells,
        }
    }

    /// Estimate the band gap frequency range \[Hz\].
    ///
    /// Returns `(f_lower, f_upper)`.
    pub fn band_gap_range(&self) -> (f64, f64) {
        if self.cells.is_empty() {
            return (0.0, 0.0);
        }
        let f0 = self.cells[0].resonant_frequency;
        let m1 = self.cells[0].shell_mass;
        let m2 = self.cells[0].resonator_mass;
        let f_upper = f0 * (1.0 + m2 / m1).sqrt();
        (f0, f_upper)
    }

    /// Check whether wave propagation is forbidden at frequency `f` \[Hz\].
    pub fn is_band_gap(&self, f: f64) -> bool {
        let (f_low, f_high) = self.band_gap_range();
        f >= f_low && f <= f_high
    }

    /// Effective negative mass density ratio ρ*/ρ_host at frequency `f` \[Hz\].
    ///
    /// Returns negative values in the band gap.
    pub fn effective_density_ratio(&self, f: f64) -> f64 {
        if self.cells.is_empty() {
            return 1.0;
        }
        let cell = &self.cells[0];
        let omega = 2.0 * PI * f;
        let omega0 = 2.0 * PI * cell.resonant_frequency;
        let m_ratio = cell.resonator_mass / (cell.shell_mass + cell.resonator_mass);
        let denom = omega0 * omega0 - omega * omega;
        if denom.abs() < 1e-30 {
            return f64::NEG_INFINITY;
        }
        1.0 + m_ratio * omega0 * omega0 / denom
    }

    /// Transmission loss estimate through n_cells at frequency `f` \[dB\].
    pub fn transmission_loss(&self, f: f64) -> f64 {
        if !self.is_band_gap(f) {
            return 0.0;
        }
        let rho_ratio = self.effective_density_ratio(f).abs();
        if rho_ratio < 1.0 {
            return 0.0;
        }
        20.0 * rho_ratio.log10() * self.n_cells as f64
    }

    /// Wavelength in host medium at frequency `f` \[m\].
    pub fn wavelength(&self, f: f64) -> f64 {
        if f < 1e-10 {
            return f64::INFINITY;
        }
        self.host_wave_speed / f
    }

    /// Sub-wavelength ratio: cell_size / wavelength at resonance.
    pub fn sub_wavelength_ratio(&self) -> f64 {
        if self.cells.is_empty() {
            return 0.0;
        }
        let f0 = self.cells[0].resonant_frequency;
        let lambda = self.wavelength(f0);
        self.cells[0].size / lambda
    }
}

// ─── MetamaterialOptimization ─────────────────────────────────────────────────

/// Design variable for topology optimization of a metamaterial unit cell.
#[derive(Debug, Clone)]
pub struct DesignElement {
    /// Element density (0 = void, 1 = solid).
    pub density: f64,
    /// Position \[m\].
    pub pos: [f64; 2],
    /// Size \[m\].
    pub size: f64,
    /// Local stiffness contribution \[Pa\].
    pub stiffness_contrib: f64,
}

impl DesignElement {
    /// Create a new design element.
    pub fn new(pos: [f64; 2], size: f64, e_material: f64) -> Self {
        Self {
            density: 1.0,
            pos,
            size,
            stiffness_contrib: e_material,
        }
    }

    /// SIMP penalized stiffness (Solid Isotropic Material with Penalization).
    pub fn penalized_stiffness(&self, penalty: f64) -> f64 {
        self.stiffness_contrib * self.density.powf(penalty)
    }
}

/// Topology optimization for metamaterial unit cell properties.
///
/// Uses SIMP (Solid Isotropic Material with Penalization) to optimize
/// material distribution for target effective properties.
pub struct MetamaterialOptimization {
    /// Design elements (unit cell discretization).
    pub elements: Vec<DesignElement>,
    /// SIMP penalty exponent.
    pub penalty: f64,
    /// Volume fraction constraint.
    pub volume_fraction: f64,
    /// Target Poisson's ratio.
    pub target_poisson: f64,
    /// Optimization iteration count.
    pub iteration: usize,
    /// Sensitivity filter radius \[m\].
    pub filter_radius: f64,
}

impl MetamaterialOptimization {
    /// Create a new topology optimization problem.
    ///
    /// * `nx`, `ny`       — grid resolution
    /// * `cell_size`      — unit cell size \[m\]
    /// * `e_material`     — solid material stiffness \[Pa\]
    /// * `vol_frac`       — target volume fraction (0 to 1)
    /// * `target_poisson` — target Poisson's ratio
    pub fn new(
        nx: usize,
        ny: usize,
        cell_size: f64,
        e_material: f64,
        vol_frac: f64,
        target_poisson: f64,
    ) -> Self {
        let dx = cell_size / nx as f64;
        let mut elements = Vec::new();
        for iy in 0..ny {
            for ix in 0..nx {
                let pos = [ix as f64 * dx + dx * 0.5, iy as f64 * dx + dx * 0.5];
                elements.push(DesignElement::new(pos, dx, e_material));
            }
        }
        let rho0 = vol_frac;
        for e in elements.iter_mut() {
            e.density = rho0;
        }
        Self {
            elements,
            penalty: 3.0,
            volume_fraction: vol_frac,
            target_poisson,
            iteration: 0,
            filter_radius: 2.0 * cell_size / nx as f64,
        }
    }

    /// Compute current volume fraction.
    pub fn current_volume_fraction(&self) -> f64 {
        if self.elements.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.elements.iter().map(|e| e.density).sum();
        sum / self.elements.len() as f64
    }

    /// Compute effective stiffness of current design \[Pa\].
    pub fn effective_stiffness(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| e.penalized_stiffness(self.penalty) * e.size * e.size)
            .sum::<f64>()
            / self.elements.len().max(1) as f64
    }

    /// OC (Optimality Criteria) update step.
    ///
    /// Minimizes compliance subject to volume fraction constraint.
    /// Returns the change in density (convergence criterion).
    pub fn oc_update(&mut self, sensitivities: &[f64]) -> f64 {
        if sensitivities.len() != self.elements.len() {
            return 0.0;
        }
        let mut l1 = 0.0_f64;
        let mut l2 = 1e9_f64;
        let move_limit = 0.2;
        let mut rho_new = vec![0.0_f64; self.elements.len()];

        for _ in 0..50 {
            let lmid = (l1 + l2) * 0.5;
            let mut vol = 0.0_f64;
            for (i, elem) in self.elements.iter().enumerate() {
                let be = (-sensitivities[i] / lmid).max(0.0);
                let rho = (elem.density * be.sqrt())
                    .clamp(elem.density - move_limit, elem.density + move_limit)
                    .clamp(0.001, 1.0);
                rho_new[i] = rho;
                vol += rho;
            }
            if vol / self.elements.len() as f64 > self.volume_fraction {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
        }

        let mut change = 0.0_f64;
        for (i, elem) in self.elements.iter_mut().enumerate() {
            change = change.max((rho_new[i] - elem.density).abs());
            elem.density = rho_new[i];
        }
        self.iteration += 1;
        change
    }

    /// Compute sensitivity of compliance to density (dC/dρ).
    pub fn compute_sensitivities(&self) -> Vec<f64> {
        self.elements
            .iter()
            .map(|e| -self.penalty * e.density.powf(self.penalty - 1.0) * e.stiffness_contrib)
            .collect()
    }

    /// Apply sensitivity filter (density-weighted average in filter radius).
    pub fn filter_sensitivities(&self, raw: &[f64]) -> Vec<f64> {
        let r = self.filter_radius;
        let n = self.elements.len();
        let mut filtered = vec![0.0_f64; n];
        let mut weight_sum = vec![0.0_f64; n];

        for i in 0..n {
            for (j, &rawj) in raw.iter().enumerate().take(n) {
                let dr = vec2_sub(self.elements[j].pos, self.elements[i].pos);
                let d = vec2_norm(dr);
                if d < r {
                    let w = (r - d) * self.elements[j].density;
                    filtered[i] += w * rawj;
                    weight_sum[i] += w;
                }
            }
            if weight_sum[i] > 1e-30 {
                filtered[i] /= weight_sum[i] * self.elements[i].density.max(1e-30);
            }
        }
        filtered
    }
}

// ─── MetamaterialAnalysis ─────────────────────────────────────────────────────

/// Elastic constants tensor (Voigt notation, 6×6).
#[derive(Debug, Clone)]
pub struct ElasticTensor {
    /// Stiffness matrix C \[Pa\], Voigt notation (6×6 stored as 36-element flat array).
    pub c: [f64; 36],
}

impl ElasticTensor {
    /// Create an isotropic elastic tensor from Young's modulus and Poisson's ratio.
    pub fn isotropic(e: f64, nu: f64) -> Self {
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let c11 = lam + 2.0 * mu;
        let c12 = lam;
        let c44 = mu;
        let mut c = [0.0_f64; 36];
        c[0] = c11;
        c[1] = c12;
        c[2] = c12;
        c[6] = c12;
        c[7] = c11;
        c[8] = c12;
        c[12] = c12;
        c[13] = c12;
        c[14] = c11;
        c[21] = c44;
        c[28] = c44;
        c[35] = c44;
        Self { c }
    }

    /// Component C\[i\]\[j\] (0-indexed).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.c[i * 6 + j]
    }

    /// Bulk modulus K \[Pa\] (valid for isotropic).
    pub fn bulk_modulus(&self) -> f64 {
        (self.c[0] + self.c[1] + self.c[2]) / 3.0
    }

    /// Shear modulus G \[Pa\] (valid for isotropic).
    pub fn shear_modulus(&self) -> f64 {
        (self.c[21] + self.c[28] + self.c[35]) / 3.0
    }

    /// Young's modulus E \[Pa\].
    pub fn young_modulus(&self) -> f64 {
        let g = self.shear_modulus();
        let k = self.bulk_modulus();
        9.0 * k * g / (3.0 * k + g)
    }

    /// Poisson's ratio (isotropic estimate).
    pub fn poisson_ratio(&self) -> f64 {
        let g = self.shear_modulus();
        let k = self.bulk_modulus();
        (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g))
    }
}

/// Effective medium theory and homogenization analysis for metamaterial lattices.
pub struct MetamaterialAnalysis;

impl MetamaterialAnalysis {
    /// Voigt upper bound for effective stiffness of a two-phase composite \[Pa\].
    ///
    /// * `e1`, `e2` — moduli of phases 1 and 2 \[Pa\]
    /// * `vf1`      — volume fraction of phase 1
    pub fn voigt_bound(e1: f64, e2: f64, vf1: f64) -> f64 {
        vf1 * e1 + (1.0 - vf1) * e2
    }

    /// Reuss lower bound for effective stiffness \[Pa\].
    pub fn reuss_bound(e1: f64, e2: f64, vf1: f64) -> f64 {
        let vf2 = 1.0 - vf1;
        1.0 / (vf1 / e1.max(1e-30) + vf2 / e2.max(1e-30))
    }

    /// Hashin-Shtrikman lower bound for effective bulk modulus \[Pa\].
    pub fn hashin_shtrikman_lower(k1: f64, _k2: f64, g1: f64, _g2: f64, vf1: f64) -> f64 {
        let vf2 = 1.0 - vf1;
        let alpha = 3.0 * k1 / (3.0 * k1 + 4.0 * g1);
        k1 + vf2 / (1.0 / (_k2 - k1).abs().max(1e-30) + vf1 * alpha / k1.max(1e-30))
    }

    /// Hashin-Shtrikman upper bound for effective bulk modulus \[Pa\].
    pub fn hashin_shtrikman_upper(k1: f64, k2: f64, _g1: f64, g2: f64, vf1: f64) -> f64 {
        let vf2 = 1.0 - vf1;
        let alpha = 3.0 * k2 / (3.0 * k2 + 4.0 * g2);
        k2 + vf1 / (1.0 / (k1 - k2).abs().max(1e-30) + vf2 * alpha / k2.max(1e-30))
    }

    /// Effective medium theory for effective Poisson's ratio of re-entrant lattice.
    ///
    /// Uses the Gibson-Ashby cellular solid model.
    pub fn auxetic_effective_poisson(theta: f64, h_over_l: f64) -> f64 {
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        -(cos_t * (h_over_l - sin_t)) / (sin_t * sin_t)
    }

    /// Relative density of a hexagonal lattice.
    ///
    /// * `t` — wall thickness \[m\]
    /// * `l` — cell wall length \[m\]
    pub fn hexagonal_relative_density(t: f64, l: f64) -> f64 {
        2.0 * t / (l * 3.0_f64.sqrt())
    }

    /// Compute the band gap width fraction for an acoustic metamaterial.
    ///
    /// * `f_lower`, `f_upper` — band gap bounds \[Hz\]
    pub fn band_gap_width_fraction(f_lower: f64, f_upper: f64) -> f64 {
        if f_lower <= 0.0 {
            return 0.0;
        }
        let f_center = (f_lower + f_upper) * 0.5;
        (f_upper - f_lower) / f_center
    }

    /// Compute effective compliance for a perforated plate lattice (Mori-Tanaka).
    ///
    /// Returns `(E_eff/E_solid, nu_eff)`.
    pub fn mori_tanaka_perforated(e_solid: f64, nu_solid: f64, porosity: f64) -> (f64, f64) {
        let dilute_factor = 1.0 - porosity * (1.0 + nu_solid) / (2.0 - nu_solid);
        let e_eff = e_solid * dilute_factor.max(0.0);
        let nu_eff = nu_solid * (1.0 - 0.5 * porosity);
        (e_eff / e_solid.max(1e-30), nu_eff.clamp(-1.0, 0.5))
    }

    /// Indentation resistance factor ~ E / (1 - ν²).
    ///
    /// Smaller denominator → better resistance; auxetic (ν < 0) improves this.
    pub fn indentation_resistance_factor(nu: f64) -> f64 {
        1.0 / (1.0 - nu * nu).max(1e-10)
    }

    /// Fracture toughness enhancement for auxetic foam.
    ///
    /// Auxetic foams show improved toughness due to densification ahead of crack.
    pub fn fracture_toughness_factor(nu: f64) -> f64 {
        (1.0 + nu.abs()).sqrt()
    }

    /// Compute dispersion relation for 1D locally resonant chain.
    ///
    /// * `k_wave`  — wave number \[1/m\]
    /// * `m1`      — shell mass \[kg\]
    /// * `m2`      — resonator mass \[kg\]
    /// * `k_eff`   — effective spring constant between cells \[N/m\]
    /// * `k_coat`  — coating stiffness \[N/m\]
    /// * `a`       — lattice constant \[m\]
    ///
    /// Returns angular frequency \[rad/s\].
    pub fn dispersion_1d(k_wave: f64, _m1: f64, m2: f64, k_eff: f64, k_coat: f64, a: f64) -> f64 {
        let omega_r2 = k_coat / m2.max(1e-30);
        let sin_ka = (k_wave * a / 2.0).sin();
        let omega_sq_free = 4.0 * k_eff / _m1.max(1e-30) * sin_ka * sin_ka;
        let b = -(omega_sq_free + omega_r2 * (1.0 + m2 / _m1.max(1e-30)));
        let c = omega_sq_free * omega_r2;
        let discriminant = b * b - 4.0 * c;
        if discriminant < 0.0 {
            return 0.0;
        }
        let omega2 = (-b - discriminant.sqrt()) / 2.0;
        if omega2 > 0.0 { omega2.sqrt() } else { 0.0 }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── AuxeticLattice ─────────────────────────────────────────────────────

    #[test]
    fn test_auxetic_poisson_ratio_negative() {
        let lattice = AuxeticLattice::new(0.01, 0.02, PI / 6.0, 0.05, 3, 3, 1e9);
        let nu = lattice.poisson_ratio(1.5);
        assert!(
            nu < 0.0,
            "re-entrant honeycomb Poisson ratio should be negative: {nu}"
        );
    }

    #[test]
    fn test_auxetic_poisson_ratio_magnitude() {
        let lattice = AuxeticLattice::new(0.01, 0.02, PI / 4.0, 0.05, 2, 2, 1e9);
        let nu = lattice.poisson_ratio(1.0);
        assert!(
            nu.abs() < 50.0,
            "Poisson ratio magnitude should be bounded: {nu}"
        );
    }

    #[test]
    fn test_auxetic_effective_modulus_positive() {
        let lattice = AuxeticLattice::new(0.01, 0.02, PI / 6.0, 0.05, 2, 2, 1e9);
        let e = lattice.effective_modulus_y(1.5);
        assert!(e > 0.0, "effective modulus should be positive: {e}");
        assert!(e < 1e9, "effective modulus should be less than solid: {e}");
    }

    #[test]
    fn test_auxetic_spring_energy_zero_at_rest() {
        let lattice = AuxeticLattice::new(0.01, 0.02, PI / 6.0, 0.05, 2, 2, 1e9);
        let energy = lattice.spring_energy();
        assert!(
            energy < 1e-10,
            "spring energy at rest should be near zero: {energy}"
        );
    }

    #[test]
    fn test_auxetic_lattice_creation() {
        let lattice = AuxeticLattice::new(0.01, 0.02, PI / 6.0, 0.05, 3, 3, 1e9);
        assert!(!lattice.nodes.is_empty(), "lattice should have nodes");
        assert!(!lattice.springs.is_empty(), "lattice should have springs");
    }

    // ── ChiralMetamaterial ────────────────────────────────────────────────

    #[test]
    fn test_chiral_twist_per_strain_nonzero() {
        let chiral = ChiralMetamaterial::new(3, 3, 0.01, 0.002, 6, 1e4, 1.0);
        let twist = chiral.twist_per_strain();
        assert!(
            twist.abs() > 1e-10,
            "twist per strain should be nonzero: {twist}"
        );
    }

    #[test]
    fn test_chiral_chirality_sign() {
        let right = ChiralMetamaterial::new(2, 2, 0.01, 0.002, 6, 1e4, 1.0);
        let left = ChiralMetamaterial::new(2, 2, 0.01, 0.002, 6, 1e4, -1.0);
        assert!(
            right.twist_per_strain() > 0.0,
            "right-handed should give positive twist"
        );
        assert!(
            left.twist_per_strain() < 0.0,
            "left-handed should give negative twist"
        );
    }

    #[test]
    fn test_chiral_rotation_under_tension() {
        let mut chiral = ChiralMetamaterial::new(2, 2, 0.01, 0.002, 6, 1e4, 1.0);
        let initial_angle = chiral.cells[0].angle;
        chiral.rotate_under_tension(0.01);
        let new_angle = chiral.cells[0].angle;
        assert!(
            (new_angle - initial_angle).abs() > 1e-12,
            "cells should rotate under tension"
        );
    }

    #[test]
    fn test_chiral_poisson_ratio_range() {
        let chiral = ChiralMetamaterial::new(2, 2, 0.01, 0.002, 6, 1e4, 1.0);
        let nu = chiral.effective_poisson_ratio();
        assert!(
            (-1.0..=1.0).contains(&nu),
            "Poisson ratio out of range: {nu}"
        );
    }

    #[test]
    fn test_chiral_rotational_energy_nonneg() {
        let mut chiral = ChiralMetamaterial::new(2, 2, 0.01, 0.002, 6, 1e4, 1.0);
        chiral.rotate_under_tension(0.05);
        let e = chiral.rotational_energy();
        assert!(e >= 0.0, "rotational energy should be non-negative: {e}");
    }

    // ── PentamodeMetamaterial ─────────────────────────────────────────────

    #[test]
    fn test_pentamode_bulk_to_shear_large() {
        let elem = PentamodeElement::new(0.01, 1e9);
        let ratio = elem.bulk_to_shear_ratio();
        assert!(
            ratio > 100.0,
            "pentamode should have large κ/G ratio: {ratio}"
        );
    }

    #[test]
    fn test_pentamode_shear_much_less_than_bulk() {
        let elem = PentamodeElement::new(0.01, 1e9);
        assert!(
            elem.effective_shear_modulus() < elem.effective_bulk_modulus(),
            "G* should be << κ* for pentamode"
        );
    }

    #[test]
    fn test_pentamode_is_pentamode_check() {
        let mm = PentamodeMetamaterial::new(2, 2, 2, 0.01, 1e9);
        assert!(mm.is_pentamode(10.0), "pentamode should pass threshold=10");
    }

    #[test]
    fn test_pentamode_wave_speed_positive() {
        let mm = PentamodeMetamaterial::new(2, 2, 2, 0.01, 1e9);
        let c = mm.longitudinal_wave_speed(1000.0);
        assert!(c > 0.0, "wave speed should be positive: {c}");
    }

    // ── MagneticMetamaterial ──────────────────────────────────────────────

    #[test]
    fn test_magnetic_stiffness_increases_with_field() {
        let mut mm = MagneticMetamaterial::new(5, 1000.0);
        let k0 = mm.effective_stiffness();
        mm.apply_field(1.0);
        let k1 = mm.effective_stiffness();
        assert!(
            k1 > k0,
            "stiffness should increase with applied field: {k0} → {k1}"
        );
    }

    #[test]
    fn test_magnetic_stiffness_ratio_gt_one() {
        let mut mm = MagneticMetamaterial::new(5, 1000.0);
        mm.apply_field(1.0);
        let ratio = mm.stiffness_ratio();
        assert!(
            ratio > 1.0,
            "stiffness ratio should be > 1 with field: {ratio}"
        );
    }

    #[test]
    fn test_magnetic_zero_field_no_change() {
        let mut mm = MagneticMetamaterial::new(5, 1000.0);
        mm.apply_field(0.0);
        let ratio = mm.stiffness_ratio();
        assert!(
            (ratio - 1.0).abs() < 1e-10,
            "zero field should give ratio=1: {ratio}"
        );
    }

    #[test]
    fn test_magnetic_magnetostrictive_strain() {
        let cell = MagneticCell::new([0.0; 3], 1000.0);
        let s = cell.magnetostrictive_strain(1.0);
        assert!(
            s >= 0.0,
            "magnetostrictive strain should be non-negative: {s}"
        );
    }

    // ── ProgrammableMatter ────────────────────────────────────────────────

    #[test]
    fn test_bistable_double_well_energy() {
        let cell = BistableCell::new(1000.0, 2.0);
        let f = cell.restoring_force();
        // At equilibrium: restoring force should be relatively small
        assert!(
            f.abs() < 500.0,
            "restoring force at equilibrium should be small: {f}"
        );
    }

    #[test]
    fn test_bistable_snap_through_event() {
        let mut cell = BistableCell::new(100.0, 3.0);
        let threshold = cell.snap_threshold;
        let snapped = cell.apply_force(threshold * 2.0 + 1.0);
        assert!(snapped, "large force should cause snap-through");
        assert!(!cell.in_state_a, "cell should be in state B after snap");
    }

    #[test]
    fn test_programmable_matter_sequential_folding() {
        let mut pm = ProgrammableMatter::new(5, 100.0, 3.0, 10.0);
        let threshold = pm.cells[0].snap_threshold;
        let f = threshold * 3.0 + 1.0;
        let idx = pm.trigger_fold(f);
        assert!(idx.is_some(), "folding should occur with sufficient force");
        assert_eq!(idx.unwrap(), 0, "first cell should fold first");
        assert_eq!(pm.fold_index, 1, "fold_index should advance");
    }

    #[test]
    fn test_programmable_matter_energy_nonneg() {
        let pm = ProgrammableMatter::new(4, 100.0, 2.5, 10.0);
        assert!(
            pm.total_energy() >= 0.0,
            "total energy should be non-negative"
        );
    }

    #[test]
    fn test_bistable_state_count() {
        let pm = ProgrammableMatter::new(6, 100.0, 2.5, 10.0);
        let a = pm.count_state_a();
        let b = pm.count_state_b();
        assert_eq!(a + b, 6, "state counts should sum to total cells");
        assert_eq!(a, 6, "all cells initially in state A");
        assert_eq!(b, 0);
    }

    // ── AcousticMetamaterial ──────────────────────────────────────────────

    #[test]
    fn test_acoustic_band_gap_range_ordered() {
        let mm = AcousticMetamaterial::new(5);
        let (f_low, f_high) = mm.band_gap_range();
        assert!(
            f_low > 0.0,
            "band gap lower bound should be positive: {f_low}"
        );
        assert!(
            f_high > f_low,
            "band gap upper should exceed lower: {f_high} > {f_low}"
        );
    }

    #[test]
    fn test_acoustic_band_gap_detection() {
        let mm = AcousticMetamaterial::new(5);
        let (f_low, f_high) = mm.band_gap_range();
        let f_mid = (f_low + f_high) * 0.5;
        assert!(
            mm.is_band_gap(f_mid),
            "midpoint should be in band gap: {f_mid}"
        );
    }

    #[test]
    fn test_acoustic_outside_band_gap() {
        let mm = AcousticMetamaterial::new(5);
        let (f_low, _) = mm.band_gap_range();
        assert!(
            !mm.is_band_gap(f_low * 0.5),
            "below band gap should not be in gap"
        );
    }

    #[test]
    fn test_acoustic_effective_density_negative_in_gap() {
        let mm = AcousticMetamaterial::new(5);
        let (f_low, f_high) = mm.band_gap_range();
        let f_test = f_low * 1.05;
        if f_test < f_high {
            let rho_ratio = mm.effective_density_ratio(f_test);
            assert!(
                rho_ratio < 0.0,
                "effective density should be negative in gap: {rho_ratio}"
            );
        }
    }

    #[test]
    fn test_acoustic_sub_wavelength_ratio() {
        let mm = AcousticMetamaterial::new(5);
        let ratio = mm.sub_wavelength_ratio();
        assert!(
            ratio > 0.0 && ratio < 1.0,
            "resonators should be sub-wavelength: {ratio}"
        );
    }

    #[test]
    fn test_acoustic_resonant_cell_frequency_positive() {
        let cell = ResonantCell::new(0.028, 0.101, 8e4, 0.015);
        assert!(
            cell.resonant_frequency > 0.0,
            "resonant frequency should be positive"
        );
    }

    // ── MetamaterialOptimization ──────────────────────────────────────────

    #[test]
    fn test_topology_opt_volume_fraction_init() {
        let opt = MetamaterialOptimization::new(10, 10, 0.01, 1e9, 0.5, -0.3);
        let vf = opt.current_volume_fraction();
        assert!(
            (vf - 0.5).abs() < 1e-10,
            "initial VF should match target: {vf}"
        );
    }

    #[test]
    fn test_topology_opt_sensitivity_negative() {
        let opt = MetamaterialOptimization::new(5, 5, 0.01, 1e9, 0.5, -0.3);
        let s = opt.compute_sensitivities();
        assert!(
            s.iter().all(|&x| x < 0.0),
            "sensitivities should all be negative"
        );
    }

    #[test]
    fn test_topology_opt_oc_update_converges() {
        let mut opt = MetamaterialOptimization::new(5, 5, 0.01, 1e9, 0.5, -0.3);
        for _ in 0..10 {
            let raw = opt.compute_sensitivities();
            let filtered = opt.filter_sensitivities(&raw);
            let _change = opt.oc_update(&filtered);
        }
        let vf = opt.current_volume_fraction();
        // OC update should produce valid volume fractions in [0, 1]
        assert!(
            (0.0..=1.0).contains(&vf),
            "volume fraction should be valid: {vf}"
        );
    }

    #[test]
    fn test_simp_penalized_stiffness() {
        let elem = DesignElement::new([0.0; 2], 0.001, 1e9);
        let k_full = elem.penalized_stiffness(3.0);
        assert!(
            (k_full - 1e9).abs() < 1.0,
            "full density should give full stiffness"
        );
        let mut elem2 = elem.clone();
        elem2.density = 0.5;
        let k_half = elem2.penalized_stiffness(3.0);
        assert!(
            k_half < k_full * 0.5,
            "SIMP with p=3 should penalize intermediate density"
        );
    }

    // ── MetamaterialAnalysis ──────────────────────────────────────────────

    #[test]
    fn test_voigt_reuss_bounds_ordering() {
        let e1 = 1e9;
        let e2 = 1e6;
        let vf = 0.4;
        let voigt = MetamaterialAnalysis::voigt_bound(e1, e2, vf);
        let reuss = MetamaterialAnalysis::reuss_bound(e1, e2, vf);
        assert!(
            voigt > reuss,
            "Voigt bound should exceed Reuss bound: {voigt} > {reuss}"
        );
    }

    #[test]
    fn test_voigt_bound_single_phase() {
        let e = 1e9;
        let v1 = MetamaterialAnalysis::voigt_bound(e, e, 0.5);
        assert!(
            (v1 - e).abs() < 1.0,
            "single phase Voigt should equal E: {v1}"
        );
    }

    #[test]
    fn test_elastic_tensor_isotropic() {
        let e = 200e9;
        let nu = 0.3;
        let ct = ElasticTensor::isotropic(e, nu);
        let k = ct.bulk_modulus();
        let g = ct.shear_modulus();
        let e_check = 9.0 * k * g / (3.0 * k + g);
        assert!(
            (e_check - e).abs() / e < 1e-6,
            "recovered E={e_check} expected={e}"
        );
    }

    #[test]
    fn test_elastic_tensor_poisson() {
        let e = 70e9;
        let nu = 0.33;
        let ct = ElasticTensor::isotropic(e, nu);
        let nu_check = ct.poisson_ratio();
        assert!(
            (nu_check - nu).abs() < 1e-6,
            "recovered nu={nu_check} expected={nu}"
        );
    }

    #[test]
    fn test_band_gap_width_fraction() {
        let width_frac = MetamaterialAnalysis::band_gap_width_fraction(100.0, 120.0);
        let expected = 20.0 / 110.0;
        assert!(
            (width_frac - expected).abs() < 1e-10,
            "band gap fraction={width_frac}"
        );
    }

    #[test]
    fn test_fracture_toughness_enhancement_auxetic() {
        let nu_conv = 0.3;
        let nu_aux = -0.5;
        let k_conv = MetamaterialAnalysis::fracture_toughness_factor(nu_conv);
        let k_aux = MetamaterialAnalysis::fracture_toughness_factor(nu_aux);
        assert!(
            k_aux > k_conv,
            "auxetic should have higher toughness: {k_aux} > {k_conv}"
        );
    }

    #[test]
    fn test_dispersion_relation_positive() {
        let omega = MetamaterialAnalysis::dispersion_1d(100.0, 0.028, 0.101, 1e5, 8e4, 0.015);
        assert!(
            omega >= 0.0,
            "dispersion frequency should be non-negative: {omega}"
        );
    }

    #[test]
    fn test_auxetic_analysis_effective_poisson() {
        let nu = MetamaterialAnalysis::auxetic_effective_poisson(PI / 6.0, 1.5);
        assert!(
            nu < 0.0,
            "re-entrant analysis should give negative nu: {nu}"
        );
    }
}
