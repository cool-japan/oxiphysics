// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Membrane mechanics for soft body simulation.
//!
//! This module provides thin shell/membrane elements, tension-only behaviour,
//! wrinkling detection (tension field theory), pneumatic pressure loading,
//! membrane inflation, form-finding (force density method), cable-membrane
//! interaction, and pre-stress models.
//!
//! All quantities are in SI units (Pa, m, N, ...) unless otherwise stated.
//! Uses `[f64; 3]` arrays for vectors (no nalgebra dependency).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector helpers (plain [f64; 3], no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-vectors.
#[inline]
fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract b from a.
#[inline]
fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a vector.
#[inline]
fn vec_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
fn vec_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length of a 3-vector.
#[inline]
fn vec_len(a: [f64; 3]) -> f64 {
    vec_dot(a, a).sqrt()
}

/// Normalize a 3-vector, returning zero if degenerate.
#[inline]
fn vec_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec_len(a);
    if len < 1e-15 {
        [0.0; 3]
    } else {
        vec_scale(a, 1.0 / len)
    }
}

/// Negate a 3-vector.
#[inline]
fn vec_neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

/// Add scaled vector: a + s * b.
#[inline]
fn vec_add_scaled(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
}

/// Squared length of a 3-vector.
#[inline]
fn vec_len_sq(a: [f64; 3]) -> f64 {
    vec_dot(a, a)
}

/// Zero 3-vector.
#[inline]
fn vec_zero() -> [f64; 3] {
    [0.0; 3]
}

// ---------------------------------------------------------------------------
// Membrane node
// ---------------------------------------------------------------------------

/// A node (vertex) in the membrane mesh.
#[derive(Debug, Clone)]
pub struct MembraneNode {
    /// Current position (m).
    pub position: [f64; 3],
    /// Previous position (for Verlet integration).
    pub prev_position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// External force accumulator (N).
    pub force: [f64; 3],
    /// Inverse mass (1/kg). Zero for fixed nodes.
    pub inv_mass: f64,
}

impl MembraneNode {
    /// Create a new membrane node with given position and mass.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            prev_position: position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass: if mass > 0.0 { 1.0 / mass } else { 0.0 },
        }
    }

    /// Create a fixed (static) membrane node.
    pub fn new_fixed(position: [f64; 3]) -> Self {
        Self {
            position,
            prev_position: position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass: 0.0,
        }
    }

    /// Returns true if the node is fixed (infinite mass).
    pub fn is_fixed(&self) -> bool {
        self.inv_mass == 0.0
    }
}

// ---------------------------------------------------------------------------
// Wrinkling state (tension field theory)
// ---------------------------------------------------------------------------

/// Wrinkling state of a membrane element based on tension field theory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WrinklingState {
    /// Both principal stresses are tensile -- fully taut.
    Taut,
    /// One principal stress is compressive -- wrinkled.
    Wrinkled,
    /// Both principal stresses are zero or compressive -- slack.
    Slack,
}

/// Determine wrinkling state from the two in-plane principal stresses.
///
/// # Arguments
/// * `sigma_1` -- first principal stress (Pa).
/// * `sigma_2` -- second principal stress (Pa).
///
/// Returns the [`WrinklingState`].
pub fn wrinkling_state(sigma_1: f64, sigma_2: f64) -> WrinklingState {
    if sigma_1 > 0.0 && sigma_2 > 0.0 {
        WrinklingState::Taut
    } else if sigma_1 > 0.0 || sigma_2 > 0.0 {
        WrinklingState::Wrinkled
    } else {
        WrinklingState::Slack
    }
}

/// Compute wrinkle wavelength for a thin membrane under uniaxial tension.
///
/// `lambda = 2 * pi * (D / N)^0.25`
///
/// where D is bending stiffness and N is membrane tension per unit width.
///
/// # Arguments
/// * `bending_stiffness` -- D = E*t^3 / (12*(1-nu^2)) (N m).
/// * `tension_per_width` -- N (N/m).
///
/// Returns wrinkle wavelength (m).
pub fn wrinkle_wavelength(bending_stiffness: f64, tension_per_width: f64) -> f64 {
    if tension_per_width <= 0.0 {
        return 0.0;
    }
    2.0 * PI * (bending_stiffness / tension_per_width).sqrt().sqrt()
}

/// Compute wrinkle amplitude from membrane strain and wavelength.
///
/// `A = lambda / pi * sqrt(epsilon_comp)`
///
/// where epsilon_comp is the compressive strain in the wrinkle direction.
///
/// # Arguments
/// * `wavelength` -- wrinkle wavelength (m).
/// * `compressive_strain` -- magnitude of compressive strain (positive value).
///
/// Returns wrinkle amplitude (m).
pub fn wrinkle_amplitude(wavelength: f64, compressive_strain: f64) -> f64 {
    if compressive_strain <= 0.0 || wavelength <= 0.0 {
        return 0.0;
    }
    wavelength / PI * compressive_strain.sqrt()
}

// ---------------------------------------------------------------------------
// Membrane triangle element
// ---------------------------------------------------------------------------

/// A constant-strain triangle (CST) membrane element.
///
/// Supports tension-only behaviour: compressive stresses are zeroed out
/// to enforce the membrane (no bending stiffness) assumption.
#[derive(Debug, Clone)]
pub struct MembraneTriangle {
    /// Node indices (into a global node array).
    pub nodes: [usize; 3],
    /// Thickness (m).
    pub thickness: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poissons_ratio: f64,
    /// Pre-stress in local x direction (Pa).
    pub prestress_x: f64,
    /// Pre-stress in local y direction (Pa).
    pub prestress_y: f64,
    /// Reference (undeformed) area (m^2).
    pub ref_area: f64,
    /// Whether to enforce tension-only behaviour.
    pub tension_only: bool,
}

impl MembraneTriangle {
    /// Create a new membrane triangle element.
    pub fn new(
        nodes: [usize; 3],
        thickness: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
        positions: &[[f64; 3]],
    ) -> Self {
        let area = triangle_area(
            positions[nodes[0]],
            positions[nodes[1]],
            positions[nodes[2]],
        );
        Self {
            nodes,
            thickness,
            youngs_modulus,
            poissons_ratio,
            prestress_x: 0.0,
            prestress_y: 0.0,
            ref_area: area,
            tension_only: true,
        }
    }

    /// Set pre-stress values for the element.
    pub fn set_prestress(&mut self, px: f64, py: f64) {
        self.prestress_x = px;
        self.prestress_y = py;
    }

    /// Compute the current area of the triangle.
    pub fn current_area(&self, positions: &[[f64; 3]]) -> f64 {
        triangle_area(
            positions[self.nodes[0]],
            positions[self.nodes[1]],
            positions[self.nodes[2]],
        )
    }

    /// Compute the area ratio (current / reference).
    pub fn area_ratio(&self, positions: &[[f64; 3]]) -> f64 {
        if self.ref_area < 1e-30 {
            return 1.0;
        }
        self.current_area(positions) / self.ref_area
    }

    /// Compute unit normal of the triangle in current configuration.
    pub fn normal(&self, positions: &[[f64; 3]]) -> [f64; 3] {
        let p0 = positions[self.nodes[0]];
        let p1 = positions[self.nodes[1]];
        let p2 = positions[self.nodes[2]];
        let e1 = vec_sub(p1, p0);
        let e2 = vec_sub(p2, p0);
        vec_normalize(vec_cross(e1, e2))
    }

    /// Compute membrane tension force per unit length (N/m).
    ///
    /// Uses a simplified isotropic plane-stress model with area strain.
    ///
    /// Returns the scalar membrane tension (N/m) = stress * thickness.
    pub fn membrane_tension(&self, positions: &[[f64; 3]]) -> f64 {
        let ratio = self.area_ratio(positions);
        let strain = ratio - 1.0;
        let factor = self.youngs_modulus / (1.0 - self.poissons_ratio * self.poissons_ratio);
        let stress = factor * strain + self.prestress_x;
        if self.tension_only && stress < 0.0 {
            0.0
        } else {
            stress * self.thickness
        }
    }

    /// Compute the wrinkling state of this element.
    pub fn element_wrinkling_state(&self, positions: &[[f64; 3]]) -> WrinklingState {
        let ratio = self.area_ratio(positions);
        let strain = ratio - 1.0;
        let factor = self.youngs_modulus / (1.0 - self.poissons_ratio * self.poissons_ratio);
        // Simplified: assume equal biaxial strain for the isotropic approximation.
        let sigma_1 = factor * strain + self.prestress_x;
        let sigma_2 = factor * strain + self.prestress_y;
        wrinkling_state(sigma_1, sigma_2)
    }

    /// Apply internal membrane forces to nodes (distribute equally to 3 vertices).
    ///
    /// Uses a simplified approach: forces are applied along the element normal
    /// proportional to the area change (like a pressure).
    pub fn apply_forces(&self, nodes: &mut [MembraneNode]) {
        let positions: Vec<[f64; 3]> = nodes.iter().map(|n| n.position).collect();
        let tension = self.membrane_tension(&positions);
        if tension.abs() < 1e-30 {
            return;
        }
        let n = self.normal(&positions);
        let area = self.current_area(&positions);
        // Force magnitude distributed to 3 nodes
        let f_mag = tension * area / 3.0;
        for &idx in &self.nodes {
            if nodes[idx].inv_mass > 0.0 {
                let f = vec_scale(n, -f_mag * (self.area_ratio(&positions) - 1.0).signum());
                nodes[idx].force = vec_add(nodes[idx].force, f);
            }
        }
    }
}

/// Compute the area of a triangle defined by three points.
pub fn triangle_area(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let e1 = vec_sub(p1, p0);
    let e2 = vec_sub(p2, p0);
    0.5 * vec_len(vec_cross(e1, e2))
}

/// Compute the normal of a triangle defined by three points.
pub fn triangle_normal(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> [f64; 3] {
    let e1 = vec_sub(p1, p0);
    let e2 = vec_sub(p2, p0);
    vec_normalize(vec_cross(e1, e2))
}

// ---------------------------------------------------------------------------
// Pneumatic pressure loading
// ---------------------------------------------------------------------------

/// Apply uniform pneumatic pressure to a triangulated membrane surface.
///
/// Distributes pressure force (F = p * A * n) equally to the three vertices
/// of each triangle.
///
/// # Arguments
/// * `nodes` -- mutable slice of membrane nodes.
/// * `triangles` -- slice of triangle index triples.
/// * `pressure` -- internal gauge pressure (Pa).
pub fn apply_pneumatic_pressure(
    nodes: &mut [MembraneNode],
    triangles: &[[usize; 3]],
    pressure: f64,
) {
    for tri in triangles {
        let p0 = nodes[tri[0]].position;
        let p1 = nodes[tri[1]].position;
        let p2 = nodes[tri[2]].position;
        let n = triangle_normal(p0, p1, p2);
        let area = triangle_area(p0, p1, p2);
        let f_per_node = vec_scale(n, pressure * area / 3.0);
        for &idx in tri {
            if nodes[idx].inv_mass > 0.0 {
                nodes[idx].force = vec_add(nodes[idx].force, f_per_node);
            }
        }
    }
}

/// Compute the enclosed volume of a closed triangulated membrane.
///
/// Uses the divergence theorem: V = (1/6) * sum(n_i dot p0_i * area_i * 2).
///
/// # Arguments
/// * `positions` -- node positions.
/// * `triangles` -- triangle index triples.
///
/// Returns the enclosed volume (m^3), sign depends on winding order.
pub fn enclosed_volume(positions: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    let mut vol = 0.0;
    for tri in triangles {
        let p0 = positions[tri[0]];
        let p1 = positions[tri[1]];
        let p2 = positions[tri[2]];
        // Signed volume of tetrahedron with origin: (1/6) * p0 . (p1 x p2)
        vol += vec_dot(p0, vec_cross(p1, p2));
    }
    vol / 6.0
}

/// Compute the total surface area of a triangulated membrane.
///
/// # Arguments
/// * `positions` -- node positions.
/// * `triangles` -- triangle index triples.
///
/// Returns the total surface area (m^2).
pub fn total_surface_area(positions: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    let mut area = 0.0;
    for tri in triangles {
        area += triangle_area(positions[tri[0]], positions[tri[1]], positions[tri[2]]);
    }
    area
}

// ---------------------------------------------------------------------------
// Membrane inflation model
// ---------------------------------------------------------------------------

/// Spherical membrane inflation model (Mooney-Rivlin rubber balloon).
///
/// For a thin spherical shell of initial radius R0, thickness t0, under
/// internal pressure p, the circumferential stretch ratio lambda satisfies:
///
/// `p = (2 * t0 / R0) * (C1 + C2 / lambda^2) * (lambda - 1/lambda^5)`
///
/// (simplified from the full Mooney-Rivlin model for biaxial extension).
///
/// # Arguments
/// * `r0` -- initial (undeformed) radius (m).
/// * `t0` -- initial wall thickness (m).
/// * `c1` -- Mooney-Rivlin constant C1 (Pa).
/// * `c2` -- Mooney-Rivlin constant C2 (Pa).
/// * `stretch` -- circumferential stretch ratio lambda (dimensionless, >= 1).
///
/// Returns the internal pressure (Pa).
pub fn spherical_inflation_pressure(r0: f64, t0: f64, c1: f64, c2: f64, stretch: f64) -> f64 {
    let lam = stretch;
    let lam2 = lam * lam;
    let lam5 = lam2 * lam2 * lam;
    (2.0 * t0 / r0) * (c1 + c2 / lam2) * (lam - 1.0 / lam5)
}

/// Compute the current radius of a spherical balloon at given stretch.
///
/// `R = R0 * lambda`
///
/// # Arguments
/// * `r0` -- initial radius (m).
/// * `stretch` -- circumferential stretch ratio.
///
/// Returns the current radius (m).
pub fn balloon_current_radius(r0: f64, stretch: f64) -> f64 {
    r0 * stretch
}

/// Compute the current thickness of a spherical balloon at given stretch.
///
/// `t = t0 / lambda^2` (incompressibility assumption).
///
/// # Arguments
/// * `t0` -- initial thickness (m).
/// * `stretch` -- circumferential stretch ratio.
///
/// Returns the current thickness (m).
pub fn balloon_current_thickness(t0: f64, stretch: f64) -> f64 {
    t0 / (stretch * stretch)
}

/// Compute the volume of a spherical balloon at given stretch.
///
/// `V = (4/3) * pi * (R0 * lambda)^3`
///
/// # Arguments
/// * `r0` -- initial radius (m).
/// * `stretch` -- circumferential stretch ratio.
///
/// Returns volume (m^3).
pub fn balloon_volume(r0: f64, stretch: f64) -> f64 {
    let r = r0 * stretch;
    4.0 / 3.0 * PI * r * r * r
}

// ---------------------------------------------------------------------------
// Form-finding: force density method
// ---------------------------------------------------------------------------

/// Force density method for membrane form-finding.
///
/// Given a cable/membrane network with `n` free nodes, `n_fixed` fixed nodes,
/// and connectivity, finds the equilibrium shape under prescribed force densities.
///
/// The force density q_e = force_per_length / length_e for each element.
///
/// The equilibrium equations are: D^T Q D x = p_ext - D^T Q D_f x_f
/// where D is the connectivity matrix and Q = diag(q_e).
///
/// This simplified solver handles a single-axis (e.g., z-coordinate) for
/// a planar boundary problem.
#[derive(Debug, Clone)]
pub struct ForceDensityMethod {
    /// Number of free nodes.
    pub n_free: usize,
    /// Number of fixed (boundary) nodes.
    pub n_fixed: usize,
    /// Element connectivity: (free_or_fixed_node_i, free_or_fixed_node_j).
    pub elements: Vec<[usize; 2]>,
    /// Force density for each element (N/m^2).
    pub force_densities: Vec<f64>,
    /// Fixed node coordinates (one component, e.g. z).
    pub fixed_coords: Vec<f64>,
    /// External load on each free node (one component).
    pub external_loads: Vec<f64>,
}

impl ForceDensityMethod {
    /// Create a new force density method problem.
    pub fn new(
        n_free: usize,
        n_fixed: usize,
        elements: Vec<[usize; 2]>,
        force_densities: Vec<f64>,
        fixed_coords: Vec<f64>,
        external_loads: Vec<f64>,
    ) -> Self {
        Self {
            n_free,
            n_fixed,
            elements,
            force_densities,
            fixed_coords,
            external_loads,
        }
    }

    /// Solve for the equilibrium coordinates of free nodes using Gauss-Seidel iteration.
    ///
    /// # Arguments
    /// * `max_iterations` -- maximum number of Gauss-Seidel iterations.
    /// * `tolerance` -- convergence tolerance on the coordinate change.
    ///
    /// Returns the equilibrium coordinates of the `n_free` nodes.
    pub fn solve(&self, max_iterations: usize, tolerance: f64) -> Vec<f64> {
        let n = self.n_free;
        let total = n + self.n_fixed;
        let mut coords = vec![0.0; total];
        // Set fixed node coordinates
        coords[n..n + self.n_fixed].copy_from_slice(&self.fixed_coords[..self.n_fixed]);

        // Build stiffness contributions per node
        // K_ii = sum of q_e for elements incident on node i
        // K_ij = -q_e for element connecting i and j
        let mut diag = vec![0.0; n];
        for (e, elem) in self.elements.iter().enumerate() {
            let q = self.force_densities[e];
            let i = elem[0];
            let j = elem[1];
            if i < n {
                diag[i] += q;
            }
            if j < n {
                diag[j] += q;
            }
        }

        // Gauss-Seidel iteration
        for _iter in 0..max_iterations {
            let mut max_change = 0.0_f64;
            for i in 0..n {
                if diag[i].abs() < 1e-30 {
                    continue;
                }
                let mut rhs = self.external_loads[i];
                for (e, elem) in self.elements.iter().enumerate() {
                    let q = self.force_densities[e];
                    if elem[0] == i {
                        rhs += q * coords[elem[1]];
                    } else if elem[1] == i {
                        rhs += q * coords[elem[0]];
                    }
                }
                let new_val = rhs / diag[i];
                max_change = max_change.max((new_val - coords[i]).abs());
                coords[i] = new_val;
            }
            if max_change < tolerance {
                break;
            }
        }

        coords[..n].to_vec()
    }
}

// ---------------------------------------------------------------------------
// Cable element for cable-membrane interaction
// ---------------------------------------------------------------------------

/// A cable element for cable-membrane interaction.
///
/// Cables carry tension only (no compression, no bending).
#[derive(Debug, Clone)]
pub struct CableElement {
    /// Node indices (two endpoints).
    pub nodes: [usize; 2],
    /// Rest length (m).
    pub rest_length: f64,
    /// Axial stiffness EA (N).
    pub axial_stiffness: f64,
    /// Pre-tension force (N).
    pub pretension: f64,
}

impl CableElement {
    /// Create a new cable element.
    pub fn new(nodes: [usize; 2], rest_length: f64, axial_stiffness: f64, pretension: f64) -> Self {
        Self {
            nodes,
            rest_length,
            axial_stiffness,
            pretension,
        }
    }

    /// Create a cable element from node positions with zero pretension.
    pub fn from_positions(nodes: [usize; 2], positions: &[[f64; 3]], axial_stiffness: f64) -> Self {
        let rest = vec_len(vec_sub(positions[nodes[1]], positions[nodes[0]]));
        Self::new(nodes, rest, axial_stiffness, 0.0)
    }

    /// Compute the current length of the cable.
    pub fn current_length(&self, positions: &[[f64; 3]]) -> f64 {
        vec_len(vec_sub(positions[self.nodes[1]], positions[self.nodes[0]]))
    }

    /// Compute the cable tension force (N). Returns 0 if in compression.
    pub fn tension(&self, positions: &[[f64; 3]]) -> f64 {
        let l = self.current_length(positions);
        if l < 1e-15 {
            return self.pretension.max(0.0);
        }
        let strain = (l - self.rest_length) / self.rest_length;
        let force = self.axial_stiffness * strain + self.pretension;
        force.max(0.0) // tension only
    }

    /// Apply cable forces to the nodes.
    pub fn apply_forces(&self, nodes: &mut [MembraneNode]) {
        let p0 = nodes[self.nodes[0]].position;
        let p1 = nodes[self.nodes[1]].position;
        let diff = vec_sub(p1, p0);
        let l = vec_len(diff);
        if l < 1e-15 {
            return;
        }
        let positions: Vec<[f64; 3]> = nodes.iter().map(|n| n.position).collect();
        let t = self.tension(&positions);
        if t < 1e-30 {
            return;
        }
        let dir = vec_scale(diff, 1.0 / l);
        let f = vec_scale(dir, t);
        if nodes[self.nodes[0]].inv_mass > 0.0 {
            nodes[self.nodes[0]].force = vec_add(nodes[self.nodes[0]].force, f);
        }
        if nodes[self.nodes[1]].inv_mass > 0.0 {
            nodes[self.nodes[1]].force = vec_add(nodes[self.nodes[1]].force, vec_neg(f));
        }
    }
}

// ---------------------------------------------------------------------------
// Pre-stress model
// ---------------------------------------------------------------------------

/// Isotropic pre-stress state for a membrane.
#[derive(Debug, Clone)]
pub struct MembranePreStress {
    /// Pre-stress in the first principal direction (Pa).
    pub sigma_1: f64,
    /// Pre-stress in the second principal direction (Pa).
    pub sigma_2: f64,
    /// Orientation angle of the first principal direction (rad).
    pub theta: f64,
}

impl MembranePreStress {
    /// Create a new pre-stress state.
    pub fn new(sigma_1: f64, sigma_2: f64, theta: f64) -> Self {
        Self {
            sigma_1,
            sigma_2,
            theta,
        }
    }

    /// Create an isotropic (equal biaxial) pre-stress.
    pub fn isotropic(sigma: f64) -> Self {
        Self::new(sigma, sigma, 0.0)
    }

    /// Compute the stress components in global x-y coordinates.
    ///
    /// Returns (sigma_xx, sigma_yy, sigma_xy).
    pub fn global_components(&self) -> (f64, f64, f64) {
        let c = self.theta.cos();
        let s = self.theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;
        let sxx = self.sigma_1 * c2 + self.sigma_2 * s2;
        let syy = self.sigma_1 * s2 + self.sigma_2 * c2;
        let sxy = (self.sigma_1 - self.sigma_2) * cs;
        (sxx, syy, sxy)
    }

    /// Compute the mean pre-stress (Pa).
    pub fn mean_stress(&self) -> f64 {
        0.5 * (self.sigma_1 + self.sigma_2)
    }

    /// Returns true if the pre-stress is purely tensile.
    pub fn is_tensile(&self) -> bool {
        self.sigma_1 >= 0.0 && self.sigma_2 >= 0.0
    }
}

// ---------------------------------------------------------------------------
// Membrane mesh
// ---------------------------------------------------------------------------

/// A membrane mesh consisting of nodes, triangles, and optional cables.
#[derive(Debug, Clone)]
pub struct MembraneMesh {
    /// Nodes in the mesh.
    pub nodes: Vec<MembraneNode>,
    /// Triangle connectivity.
    pub triangles: Vec<MembraneTriangle>,
    /// Cable elements.
    pub cables: Vec<CableElement>,
    /// Internal gauge pressure (Pa).
    pub pressure: f64,
    /// Gravity acceleration vector (m/s^2).
    pub gravity: [f64; 3],
    /// Damping coefficient for dynamic relaxation.
    pub damping: f64,
}

impl Default for MembraneMesh {
    fn default() -> Self {
        Self::new()
    }
}

impl MembraneMesh {
    /// Create a new empty membrane mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            triangles: Vec::new(),
            cables: Vec::new(),
            pressure: 0.0,
            gravity: [0.0, -9.81, 0.0],
            damping: 0.01,
        }
    }

    /// Create a flat rectangular membrane mesh in the XZ plane.
    ///
    /// # Arguments
    /// * `nx` -- number of nodes in x direction.
    /// * `nz` -- number of nodes in z direction.
    /// * `lx` -- total length in x (m).
    /// * `lz` -- total length in z (m).
    /// * `mass_per_node` -- mass assigned to each node (kg).
    /// * `thickness` -- membrane thickness (m).
    /// * `youngs_modulus` -- E (Pa).
    /// * `poissons_ratio` -- nu.
    pub fn create_rectangular(
        nx: usize,
        nz: usize,
        lx: f64,
        lz: f64,
        mass_per_node: f64,
        thickness: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
    ) -> Self {
        let mut mesh = Self::new();
        let dx = if nx > 1 { lx / (nx - 1) as f64 } else { 0.0 };
        let dz = if nz > 1 { lz / (nz - 1) as f64 } else { 0.0 };

        // Create nodes
        for iz in 0..nz {
            for ix in 0..nx {
                let pos = [ix as f64 * dx, 0.0, iz as f64 * dz];
                mesh.nodes.push(MembraneNode::new(pos, mass_per_node));
            }
        }

        // Create triangles (two per quad)
        let positions: Vec<[f64; 3]> = mesh.nodes.iter().map(|n| n.position).collect();
        for iz in 0..(nz - 1) {
            for ix in 0..(nx - 1) {
                let i00 = iz * nx + ix;
                let i10 = iz * nx + ix + 1;
                let i01 = (iz + 1) * nx + ix;
                let i11 = (iz + 1) * nx + ix + 1;
                // Lower-left triangle
                mesh.triangles.push(MembraneTriangle::new(
                    [i00, i10, i01],
                    thickness,
                    youngs_modulus,
                    poissons_ratio,
                    &positions,
                ));
                // Upper-right triangle
                mesh.triangles.push(MembraneTriangle::new(
                    [i10, i11, i01],
                    thickness,
                    youngs_modulus,
                    poissons_ratio,
                    &positions,
                ));
            }
        }

        mesh
    }

    /// Fix the boundary nodes of a rectangular mesh.
    pub fn fix_boundary(&mut self, nx: usize, nz: usize) {
        for iz in 0..nz {
            for ix in 0..nx {
                if ix == 0 || ix == nx - 1 || iz == 0 || iz == nz - 1 {
                    let idx = iz * nx + ix;
                    self.nodes[idx].inv_mass = 0.0;
                }
            }
        }
    }

    /// Get positions as a flat vector of \[f64; 3\].
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.nodes.iter().map(|n| n.position).collect()
    }

    /// Number of nodes in the mesh.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of triangles in the mesh.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Total surface area of the membrane.
    pub fn surface_area(&self) -> f64 {
        let pos = self.positions();
        let tris: Vec<[usize; 3]> = self.triangles.iter().map(|t| t.nodes).collect();
        total_surface_area(&pos, &tris)
    }

    /// Clear all force accumulators.
    pub fn clear_forces(&mut self) {
        for node in &mut self.nodes {
            node.force = vec_zero();
        }
    }

    /// Apply gravity forces to all nodes.
    pub fn apply_gravity(&mut self) {
        let g = self.gravity;
        for node in &mut self.nodes {
            if node.inv_mass > 0.0 {
                let mass = 1.0 / node.inv_mass;
                let f = vec_scale(g, mass);
                node.force = vec_add(node.force, f);
            }
        }
    }

    /// Apply pneumatic pressure to the membrane surface.
    pub fn apply_pressure(&mut self) {
        let p = self.pressure;
        if p.abs() < 1e-30 {
            return;
        }
        let tris: Vec<[usize; 3]> = self.triangles.iter().map(|t| t.nodes).collect();
        apply_pneumatic_pressure(&mut self.nodes, &tris, p);
    }

    /// Apply cable forces.
    pub fn apply_cable_forces(&mut self) {
        let cables = self.cables.clone();
        for cable in &cables {
            cable.apply_forces(&mut self.nodes);
        }
    }

    /// Perform one time step using Verlet integration with damping.
    ///
    /// # Arguments
    /// * `dt` -- time step (s).
    pub fn step(&mut self, dt: f64) {
        self.clear_forces();
        self.apply_gravity();
        self.apply_pressure();
        self.apply_cable_forces();

        let damping = self.damping;
        let dt2 = dt * dt;

        for node in &mut self.nodes {
            if node.inv_mass == 0.0 {
                continue;
            }
            let acc = vec_scale(node.force, node.inv_mass);
            // Verlet integration with damping
            let new_pos = [
                node.position[0] * (2.0 - damping) - node.prev_position[0] * (1.0 - damping)
                    + acc[0] * dt2,
                node.position[1] * (2.0 - damping) - node.prev_position[1] * (1.0 - damping)
                    + acc[1] * dt2,
                node.position[2] * (2.0 - damping) - node.prev_position[2] * (1.0 - damping)
                    + acc[2] * dt2,
            ];
            node.velocity = vec_scale(vec_sub(new_pos, node.position), 1.0 / dt);
            node.prev_position = node.position;
            node.position = new_pos;
        }
    }

    /// Run multiple simulation steps.
    pub fn simulate(&mut self, dt: f64, steps: usize) {
        for _ in 0..steps {
            self.step(dt);
        }
    }
}

// ---------------------------------------------------------------------------
// Membrane bending stiffness
// ---------------------------------------------------------------------------

/// Compute membrane/shell bending stiffness D.
///
/// `D = E * t^3 / (12 * (1 - nu^2))`
///
/// # Arguments
/// * `youngs_modulus` -- E (Pa).
/// * `thickness` -- t (m).
/// * `poissons_ratio` -- nu.
///
/// Returns bending stiffness (N m).
pub fn bending_stiffness(youngs_modulus: f64, thickness: f64, poissons_ratio: f64) -> f64 {
    youngs_modulus * thickness.powi(3) / (12.0 * (1.0 - poissons_ratio * poissons_ratio))
}

/// Compute critical buckling stress for a circular membrane under compression.
///
/// `sigma_cr = k * D / (R^2 * t)`
///
/// where k is a buckling coefficient (k ~ 3.6 for clamped edges).
///
/// # Arguments
/// * `stiffness_d` -- bending stiffness D (N m).
/// * `radius` -- membrane radius R (m).
/// * `thickness` -- t (m).
/// * `buckling_coeff` -- k (dimensionless, ~3.6 for clamped, ~1.22 for simply-supported).
///
/// Returns critical buckling stress (Pa).
pub fn critical_buckling_stress(
    stiffness_d: f64,
    radius: f64,
    thickness: f64,
    buckling_coeff: f64,
) -> f64 {
    buckling_coeff * stiffness_d / (radius * radius * thickness)
}

// ---------------------------------------------------------------------------
// Membrane stress from pressure (thin-walled vessels)
// ---------------------------------------------------------------------------

/// Hoop stress in a thin-walled cylindrical pressure vessel.
///
/// `sigma_hoop = p * R / t`
///
/// # Arguments
/// * `pressure` -- internal gauge pressure (Pa).
/// * `radius` -- cylinder radius (m).
/// * `thickness` -- wall thickness (m).
///
/// Returns hoop stress (Pa).
pub fn cylindrical_hoop_stress(pressure: f64, radius: f64, thickness: f64) -> f64 {
    pressure * radius / thickness
}

/// Hoop stress in a thin-walled spherical pressure vessel.
///
/// `sigma = p * R / (2 * t)`
///
/// # Arguments
/// * `pressure` -- internal gauge pressure (Pa).
/// * `radius` -- sphere radius (m).
/// * `thickness` -- wall thickness (m).
///
/// Returns hoop stress (Pa).
pub fn spherical_hoop_stress(pressure: f64, radius: f64, thickness: f64) -> f64 {
    pressure * radius / (2.0 * thickness)
}

/// Axial stress in a thin-walled cylindrical pressure vessel.
///
/// `sigma_axial = p * R / (2 * t)`
///
/// # Arguments
/// * `pressure` -- internal gauge pressure (Pa).
/// * `radius` -- cylinder radius (m).
/// * `thickness` -- wall thickness (m).
///
/// Returns axial stress (Pa).
pub fn cylindrical_axial_stress(pressure: f64, radius: f64, thickness: f64) -> f64 {
    pressure * radius / (2.0 * thickness)
}

// ---------------------------------------------------------------------------
// Membrane strain energy
// ---------------------------------------------------------------------------

/// Compute membrane strain energy density for isotropic plane stress.
///
/// `U = (E * t / (2 * (1 - nu^2))) * (eps_x^2 + eps_y^2 + 2*nu*eps_x*eps_y + ((1-nu)/2)*gamma_xy^2)`
///
/// # Arguments
/// * `youngs_modulus` -- E (Pa).
/// * `thickness` -- t (m).
/// * `poissons_ratio` -- nu.
/// * `eps_x` -- normal strain in x.
/// * `eps_y` -- normal strain in y.
/// * `gamma_xy` -- shear strain.
///
/// Returns strain energy per unit area (J/m^2).
pub fn membrane_strain_energy_density(
    youngs_modulus: f64,
    thickness: f64,
    poissons_ratio: f64,
    eps_x: f64,
    eps_y: f64,
    gamma_xy: f64,
) -> f64 {
    let nu = poissons_ratio;
    let factor = youngs_modulus * thickness / (2.0 * (1.0 - nu * nu));
    factor
        * (eps_x * eps_x
            + eps_y * eps_y
            + 2.0 * nu * eps_x * eps_y
            + (1.0 - nu) / 2.0 * gamma_xy * gamma_xy)
}

// ---------------------------------------------------------------------------
// Airbag deployment model
// ---------------------------------------------------------------------------

/// Simple airbag deployment model: pressure-volume relationship.
///
/// Given an ideal gas inflator:
/// `p * V = n * R_gas * T`
///
/// As volume increases from inflation, pressure drops. This function
/// computes the pressure given inflated volume and gas parameters.
///
/// # Arguments
/// * `moles` -- amount of gas n (mol).
/// * `temperature` -- gas temperature T (K).
/// * `volume` -- current volume V (m^3).
///
/// Returns internal gauge pressure (Pa, above atmospheric).
pub fn airbag_pressure(moles: f64, temperature: f64, volume: f64) -> f64 {
    const R_GAS_CONST: f64 = 8.314;
    const P_ATM: f64 = 101_325.0;
    let p_abs = moles * R_GAS_CONST * temperature / volume;
    (p_abs - P_ATM).max(0.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- wrinkling state ---

    #[test]
    fn test_wrinkling_taut() {
        assert_eq!(wrinkling_state(100.0, 50.0), WrinklingState::Taut);
    }

    #[test]
    fn test_wrinkling_wrinkled() {
        assert_eq!(wrinkling_state(100.0, -10.0), WrinklingState::Wrinkled);
    }

    #[test]
    fn test_wrinkling_slack() {
        assert_eq!(wrinkling_state(-5.0, -10.0), WrinklingState::Slack);
    }

    #[test]
    fn test_wrinkling_zero_stresses_slack() {
        assert_eq!(wrinkling_state(0.0, 0.0), WrinklingState::Slack);
    }

    // --- wrinkle wavelength and amplitude ---

    #[test]
    fn test_wrinkle_wavelength_positive() {
        let lam = wrinkle_wavelength(1e-6, 100.0);
        assert!(
            lam > 0.0,
            "wavelength must be positive for positive tension"
        );
    }

    #[test]
    fn test_wrinkle_wavelength_zero_tension() {
        let lam = wrinkle_wavelength(1e-6, 0.0);
        assert!((lam).abs() < EPS, "wavelength is zero for zero tension");
    }

    #[test]
    fn test_wrinkle_amplitude_positive() {
        let amp = wrinkle_amplitude(0.01, 0.001);
        assert!(amp > 0.0, "amplitude must be positive for positive strain");
    }

    #[test]
    fn test_wrinkle_amplitude_zero_strain() {
        let amp = wrinkle_amplitude(0.01, 0.0);
        assert!(amp.abs() < EPS, "amplitude is zero for zero strain");
    }

    // --- triangle geometry ---

    #[test]
    fn test_triangle_area_unit() {
        let area = triangle_area([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((area - 0.5).abs() < EPS, "unit right triangle area = 0.5");
    }

    #[test]
    fn test_triangle_normal_unit() {
        let n = triangle_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((n[2] - 1.0).abs() < EPS, "normal should be +Z");
    }

    // --- pneumatic pressure ---

    #[test]
    fn test_pneumatic_pressure_applies_force() {
        let mut nodes = vec![
            MembraneNode::new([0.0, 0.0, 0.0], 1.0),
            MembraneNode::new([1.0, 0.0, 0.0], 1.0),
            MembraneNode::new([0.0, 0.0, 1.0], 1.0),
        ];
        let tris = vec![[0, 1, 2]];
        apply_pneumatic_pressure(&mut nodes, &tris, 1000.0);
        // Triangle [0,0,0]->[1,0,0]->[0,0,1] has normal in -Y direction
        // so positive pressure pushes in -Y. Check that force magnitude is nonzero.
        let total_fy: f64 = nodes.iter().map(|n| n.force[1]).sum();
        assert!(
            total_fy.abs() > 0.0,
            "pressure should produce force, got {:.6}",
            total_fy
        );
    }

    #[test]
    fn test_pneumatic_pressure_zero_pressure() {
        let mut nodes = vec![
            MembraneNode::new([0.0, 0.0, 0.0], 1.0),
            MembraneNode::new([1.0, 0.0, 0.0], 1.0),
            MembraneNode::new([0.0, 0.0, 1.0], 1.0),
        ];
        let tris = vec![[0, 1, 2]];
        apply_pneumatic_pressure(&mut nodes, &tris, 0.0);
        for n in &nodes {
            assert!(vec_len(n.force) < EPS, "zero pressure -> zero force");
        }
    }

    // --- enclosed volume ---

    #[test]
    fn test_enclosed_volume_tetrahedron() {
        // Unit tetrahedron with one vertex at origin
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // Closed tetrahedron has 4 triangles (outward-facing)
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let vol = enclosed_volume(&positions, &tris).abs();
        assert!(
            (vol - 1.0 / 6.0).abs() < 1e-10,
            "unit tet volume = 1/6, got {:.6}",
            vol
        );
    }

    // --- surface area ---

    #[test]
    fn test_total_surface_area() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let area = total_surface_area(&positions, &tris);
        assert!((area - 0.5).abs() < EPS, "single triangle area");
    }

    // --- spherical inflation ---

    #[test]
    fn test_spherical_inflation_zero_stretch_zero_pressure() {
        let p = spherical_inflation_pressure(0.1, 0.001, 1e5, 0.0, 1.0);
        assert!(
            p.abs() < 1e-6,
            "at lambda=1, pressure should be ~0, got {:.6}",
            p
        );
    }

    #[test]
    fn test_spherical_inflation_pressure_positive_for_stretch() {
        let p = spherical_inflation_pressure(0.1, 0.001, 1e5, 5e4, 1.5);
        assert!(p > 0.0, "stretched balloon -> positive pressure");
    }

    #[test]
    fn test_balloon_volume_increases_with_stretch() {
        let v1 = balloon_volume(0.1, 1.0);
        let v2 = balloon_volume(0.1, 2.0);
        assert!(v2 > v1, "larger stretch -> larger volume");
    }

    #[test]
    fn test_balloon_thickness_decreases_with_stretch() {
        let t1 = balloon_current_thickness(0.001, 1.0);
        let t2 = balloon_current_thickness(0.001, 2.0);
        assert!(t2 < t1, "stretch thins the wall");
    }

    // --- force density method ---

    #[test]
    fn test_force_density_single_free_node() {
        // One free node connected to two fixed nodes at z=0 and z=1
        // With equal force densities and no load, solution should be z=0.5
        let fdm = ForceDensityMethod::new(
            1,
            2,
            vec![[0, 1], [0, 2]], // free node 0 connected to fixed nodes 1 and 2
            vec![1.0, 1.0],       // equal force densities
            vec![0.0, 1.0],       // fixed at z=0 and z=1
            vec![0.0],            // no external load
        );
        let result = fdm.solve(100, 1e-12);
        assert!(
            (result[0] - 0.5).abs() < 1e-10,
            "midpoint should be at 0.5, got {:.6}",
            result[0]
        );
    }

    #[test]
    fn test_force_density_with_load() {
        // Single free node with downward load, should shift below midpoint
        let fdm = ForceDensityMethod::new(
            1,
            2,
            vec![[0, 1], [0, 2]],
            vec![1.0, 1.0],
            vec![0.0, 1.0],
            vec![-0.5], // downward load
        );
        let result = fdm.solve(100, 1e-12);
        assert!(
            result[0] < 0.5,
            "loaded node should be below midpoint, got {:.6}",
            result[0]
        );
    }

    // --- cable element ---

    #[test]
    fn test_cable_tension_stretched() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let cable = CableElement::new([0, 1], 1.0, 1000.0, 0.0);
        let t = cable.tension(&positions);
        assert!(t > 0.0, "stretched cable has positive tension");
    }

    #[test]
    fn test_cable_tension_compressed_is_zero() {
        let positions = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let cable = CableElement::new([0, 1], 1.0, 1000.0, 0.0);
        let t = cable.tension(&positions);
        assert!(t.abs() < EPS, "compressed cable has zero tension");
    }

    #[test]
    fn test_cable_pretension() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let cable = CableElement::new([0, 1], 1.0, 1000.0, 50.0);
        let t = cable.tension(&positions);
        assert!(
            (t - 50.0).abs() < EPS,
            "at rest length, tension = pretension, got {:.6}",
            t
        );
    }

    // --- pre-stress ---

    #[test]
    fn test_prestress_isotropic_components() {
        let ps = MembranePreStress::isotropic(100e3);
        let (sxx, syy, sxy) = ps.global_components();
        assert!((sxx - 100e3).abs() < EPS, "isotropic sxx");
        assert!((syy - 100e3).abs() < EPS, "isotropic syy");
        assert!(sxy.abs() < EPS, "isotropic sxy = 0");
    }

    #[test]
    fn test_prestress_mean() {
        let ps = MembranePreStress::new(100e3, 200e3, 0.0);
        assert!((ps.mean_stress() - 150e3).abs() < EPS, "mean stress");
    }

    #[test]
    fn test_prestress_is_tensile() {
        let ps = MembranePreStress::new(100e3, 200e3, 0.0);
        assert!(ps.is_tensile(), "both positive -> tensile");
        let ps2 = MembranePreStress::new(100e3, -50e3, 0.0);
        assert!(!ps2.is_tensile(), "one negative -> not tensile");
    }

    // --- bending stiffness and buckling ---

    #[test]
    fn test_bending_stiffness_formula() {
        let d = bending_stiffness(200e9, 0.001, 0.3);
        let expected = 200e9 * 0.001_f64.powi(3) / (12.0 * (1.0 - 0.09));
        assert!((d - expected).abs() < 1e-6, "D formula");
    }

    #[test]
    fn test_critical_buckling_stress_positive() {
        let d = bending_stiffness(200e9, 0.001, 0.3);
        let sigma_cr = critical_buckling_stress(d, 0.1, 0.001, 3.6);
        assert!(sigma_cr > 0.0, "critical stress must be positive");
    }

    // --- thin-wall stresses ---

    #[test]
    fn test_cylindrical_hoop_stress_formula() {
        let sigma = cylindrical_hoop_stress(1e6, 0.5, 0.01);
        assert!((sigma - 50e6).abs() < 1.0, "hoop = pR/t");
    }

    #[test]
    fn test_spherical_hoop_half_cylindrical() {
        let p = 1e6;
        let r = 0.5;
        let t = 0.01;
        let s_cyl = cylindrical_hoop_stress(p, r, t);
        let s_sph = spherical_hoop_stress(p, r, t);
        assert!(
            (s_sph - s_cyl / 2.0).abs() < 1.0,
            "spherical hoop = cylindrical hoop / 2"
        );
    }

    // --- strain energy density ---

    #[test]
    fn test_strain_energy_zero_strain() {
        let u = membrane_strain_energy_density(200e9, 0.001, 0.3, 0.0, 0.0, 0.0);
        assert!(u.abs() < EPS, "zero strain -> zero energy");
    }

    #[test]
    fn test_strain_energy_positive_for_nonzero_strain() {
        let u = membrane_strain_energy_density(200e9, 0.001, 0.3, 0.01, 0.0, 0.0);
        assert!(u > 0.0, "nonzero strain -> positive energy");
    }

    // --- membrane mesh ---

    #[test]
    fn test_rectangular_mesh_topology() {
        let mesh = MembraneMesh::create_rectangular(4, 3, 3.0, 2.0, 1.0, 0.001, 1e6, 0.3);
        assert_eq!(mesh.num_nodes(), 12, "4*3 = 12 nodes");
        assert_eq!(mesh.num_triangles(), 12, "3*2*2 = 12 triangles");
    }

    #[test]
    fn test_rectangular_mesh_surface_area() {
        let mesh = MembraneMesh::create_rectangular(3, 3, 2.0, 2.0, 1.0, 0.001, 1e6, 0.3);
        let area = mesh.surface_area();
        assert!(
            (area - 4.0).abs() < 1e-10,
            "2x2 flat membrane area = 4 m^2, got {:.6}",
            area
        );
    }

    #[test]
    fn test_mesh_step_moves_nodes_under_gravity() {
        let mut mesh = MembraneMesh::create_rectangular(3, 3, 1.0, 1.0, 1.0, 0.001, 1e6, 0.3);
        // Fix corners only
        mesh.nodes[0].inv_mass = 0.0;
        mesh.nodes[2].inv_mass = 0.0;
        mesh.nodes[6].inv_mass = 0.0;
        mesh.nodes[8].inv_mass = 0.0;
        let y_before = mesh.nodes[4].position[1]; // center node
        mesh.simulate(1.0 / 60.0, 10);
        let y_after = mesh.nodes[4].position[1];
        assert!(
            y_after < y_before,
            "center node should sag under gravity, before={:.6} after={:.6}",
            y_before,
            y_after
        );
    }

    // --- airbag model ---

    #[test]
    fn test_airbag_pressure_positive_for_small_volume() {
        let p = airbag_pressure(1.0, 300.0, 0.001);
        assert!(
            p > 0.0,
            "high moles in small volume -> positive gauge pressure"
        );
    }

    #[test]
    fn test_airbag_pressure_decreases_with_volume() {
        // Use enough moles that both are above atmospheric
        let p1 = airbag_pressure(1.0, 300.0, 0.01);
        let p2 = airbag_pressure(1.0, 300.0, 0.1);
        assert!(p1 > p2, "larger volume -> lower pressure");
    }

    // --- membrane element ---

    #[test]
    fn test_membrane_triangle_ref_area() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tri = MembraneTriangle::new([0, 1, 2], 0.001, 1e6, 0.3, &positions);
        assert!(
            (tri.ref_area - 0.5).abs() < EPS,
            "ref area of unit right triangle"
        );
    }

    #[test]
    fn test_membrane_triangle_tension_only() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0], // compressed (original was at 1.0)
            [0.0, 0.5, 0.0], // compressed
        ];
        let ref_positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut tri = MembraneTriangle::new([0, 1, 2], 0.001, 1e6, 0.3, &ref_positions);
        tri.tension_only = true;
        let tension = tri.membrane_tension(&positions);
        assert!(
            tension.abs() < EPS,
            "compressed membrane has zero tension (tension-only)"
        );
    }
}
