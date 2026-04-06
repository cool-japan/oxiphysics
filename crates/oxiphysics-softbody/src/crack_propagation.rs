// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crack propagation and remeshing for soft body simulation.
//!
//! Implements fracture mechanics concepts including crack tip tracking,
//! stress intensity factor computation, and mesh splitting.

#![allow(dead_code)]

/// A crack tip in the simulation mesh.
#[derive(Debug, Clone)]
pub struct CrackTip {
    /// World-space position of the crack tip.
    pub position: [f64; 3],
    /// Unit direction of crack propagation.
    pub direction: [f64; 3],
    /// Mode-I stress intensity factor (Pa·√m).
    pub stress_intensity: f64,
    /// Current propagation velocity (m/s).
    pub propagation_velocity: f64,
}

impl CrackTip {
    /// Create a new crack tip.
    pub fn new(position: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            position,
            direction,
            stress_intensity: 0.0,
            propagation_velocity: 0.0,
        }
    }
}

/// Material fracture criteria.
#[derive(Debug, Clone)]
pub struct CrackCriteria {
    /// Mode-I fracture toughness (MPa·√m).
    pub k_ic: f64,
    /// Maximum crack propagation speed (m/s).
    pub max_propagation_speed: f64,
    /// Minimum element size for remeshing (m).
    pub min_element_size: f64,
}

impl CrackCriteria {
    /// Fracture criteria for glass (K_Ic ≈ 0.7 MPa·√m).
    pub fn glass() -> Self {
        Self {
            k_ic: 0.7,
            max_propagation_speed: 1500.0,
            min_element_size: 1e-4,
        }
    }

    /// Fracture criteria for steel (K_Ic ≈ 50 MPa·√m).
    pub fn steel() -> Self {
        Self {
            k_ic: 50.0,
            max_propagation_speed: 500.0,
            min_element_size: 1e-3,
        }
    }

    /// Returns true if the applied stress intensity exceeds fracture toughness.
    pub fn should_propagate(&self, k_i: f64) -> bool {
        k_i > self.k_ic
    }

    /// Compute crack propagation speed: v = v_max * (1 - K_Ic/K_I), clamped to \[0, v_max\].
    pub fn propagation_speed(&self, k_i: f64) -> f64 {
        if k_i <= self.k_ic {
            return 0.0;
        }
        let v = self.max_propagation_speed * (1.0 - self.k_ic / k_i);
        v.clamp(0.0, self.max_propagation_speed)
    }
}

/// A node (particle) in the soft body mesh.
#[derive(Debug, Clone)]
pub struct ParticleNode {
    /// World-space position.
    pub position: [f64; 3],
    /// Velocity.
    pub velocity: [f64; 3],
    /// Node mass (kg).
    pub mass: f64,
    /// Unique node identifier.
    pub id: usize,
}

impl ParticleNode {
    /// Create a new particle node.
    pub fn new(position: [f64; 3], mass: f64, id: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            id,
        }
    }
}

/// A triangulated soft body mesh.
#[derive(Debug, Clone)]
pub struct SoftBodyMesh {
    /// All particle nodes.
    pub nodes: Vec<ParticleNode>,
    /// Edges as pairs of node indices.
    pub edges: Vec<[usize; 2]>,
    /// Triangles as triples of node indices.
    pub triangles: Vec<[usize; 3]>,
}

impl SoftBodyMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Add a node and return its index.
    pub fn add_node(&mut self, pos: [f64; 3], mass: f64) -> usize {
        let id = self.nodes.len();
        self.nodes.push(ParticleNode::new(pos, mass, id));
        id
    }

    /// Add an edge between nodes i and j.
    pub fn add_edge(&mut self, i: usize, j: usize) {
        self.edges.push([i, j]);
    }

    /// Add a triangle with corners i, j, k.
    pub fn add_triangle(&mut self, i: usize, j: usize, k: usize) {
        self.triangles.push([i, j, k]);
    }

    /// Compute the length of edge at `edge_idx`.
    pub fn edge_length(&self, edge_idx: usize) -> f64 {
        let [a, b] = self.edges[edge_idx];
        let pa = &self.nodes[a].position;
        let pb = &self.nodes[b].position;
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Number of nodes in the mesh.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
}

impl Default for SoftBodyMesh {
    fn default() -> Self {
        Self::new()
    }
}

/// A propagating crack with history tracking.
#[derive(Debug, Clone)]
pub struct CrackPath {
    /// Active crack tips.
    pub tips: Vec<CrackTip>,
    /// History of crack tip positions.
    pub history: Vec<[f64; 3]>,
}

impl CrackPath {
    /// Create a new crack path starting at `start` propagating in `direction`.
    pub fn new(start: [f64; 3], direction: [f64; 3]) -> Self {
        let tip = CrackTip::new(start, direction);
        Self {
            tips: vec![tip],
            history: vec![start],
        }
    }

    /// Advance the primary crack tip by `speed * dt` along its direction.
    pub fn advance_tip(&mut self, dt: f64, speed: f64) {
        if let Some(tip) = self.tips.first_mut() {
            tip.propagation_velocity = speed;
            let step = speed * dt;
            tip.position[0] += tip.direction[0] * step;
            tip.position[1] += tip.direction[1] * step;
            tip.position[2] += tip.direction[2] * step;
            self.history.push(tip.position);
        }
    }

    /// Return a reference to the primary (first) crack tip, if any.
    pub fn current_tip(&self) -> Option<&CrackTip> {
        self.tips.first()
    }

    /// Total arc length of the crack path history.
    pub fn path_length(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        self.history
            .windows(2)
            .map(|w| {
                let a = w[0];
                let b = w[1];
                let dx = b[0] - a[0];
                let dy = b[1] - a[1];
                let dz = b[2] - a[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .sum()
    }
}

/// Split edge `edge_idx` at its midpoint, inserting a new node.
///
/// The original edge is replaced by two new edges sharing the midpoint.
/// Returns the index of the new midpoint node.
pub fn split_edge(mesh: &mut SoftBodyMesh, edge_idx: usize) -> usize {
    let [a, b] = mesh.edges[edge_idx];
    let pa = mesh.nodes[a].position;
    let pb = mesh.nodes[b].position;
    let mid = [
        (pa[0] + pb[0]) * 0.5,
        (pa[1] + pb[1]) * 0.5,
        (pa[2] + pb[2]) * 0.5,
    ];
    let avg_mass = (mesh.nodes[a].mass + mesh.nodes[b].mass) * 0.5;
    let mid_idx = mesh.add_node(mid, avg_mass);

    // Replace original edge with [a, mid].
    mesh.edges[edge_idx] = [a, mid_idx];
    // Add the second half.
    mesh.edges.push([mid_idx, b]);

    mid_idx
}

/// Propagate the crack one step if `K_I > K_Ic`.
///
/// Advances the crack tip, finds and splits the nearest mesh edge, and
/// returns `true` if propagation occurred.
pub fn propagate_crack(
    mesh: &mut SoftBodyMesh,
    crack: &mut CrackPath,
    criteria: &CrackCriteria,
    k_i: f64,
    dt: f64,
) -> bool {
    if !criteria.should_propagate(k_i) {
        return false;
    }
    let speed = criteria.propagation_speed(k_i);
    crack.advance_tip(dt, speed);

    if let Some(tip_pos) = crack.current_tip().map(|t| t.position)
        && let Some(edge_idx) = nearest_edge_to_point(mesh, tip_pos)
    {
        split_edge(mesh, edge_idx);
    }
    true
}

/// Compute a simplified mode-I stress intensity factor near the crack tip.
///
/// K_I ≈ E·√(π·a) · mean(|u| / √r) for nodes near the tip,
/// where `a` is the current crack half-length and `r` is the node–tip distance.
/// E (Young's modulus) is taken as 1.0 (normalised) in this simplified form.
pub fn compute_stress_intensity(
    crack: &CrackPath,
    mesh: &SoftBodyMesh,
    displacements: &[[f64; 3]],
) -> f64 {
    let tip = match crack.current_tip() {
        Some(t) => t,
        None => return 0.0,
    };

    let a = crack.path_length().max(1e-12);
    let pi_a = std::f64::consts::PI * a;

    let mut sum = 0.0_f64;
    let mut count = 0usize;

    let search_radius = a * 10.0_f64.max(1e-3);

    for (i, node) in mesh.nodes.iter().enumerate() {
        if i >= displacements.len() {
            break;
        }
        let dx = node.position[0] - tip.position[0];
        let dy = node.position[1] - tip.position[1];
        let dz = node.position[2] - tip.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-12 || r > search_radius {
            continue;
        }
        let u = &displacements[i];
        let u_mag = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        sum += u_mag / r.sqrt();
        count += 1;
    }

    if count == 0 {
        return 0.0;
    }
    let avg = sum / count as f64;
    pi_a.sqrt() * avg
}

/// Find the index of the mesh edge whose midpoint is closest to `point`.
pub fn nearest_edge_to_point(mesh: &SoftBodyMesh, point: [f64; 3]) -> Option<usize> {
    if mesh.edges.is_empty() {
        return None;
    }
    let mut best_idx = 0;
    let mut best_dist2 = f64::MAX;

    for (i, &[a, b]) in mesh.edges.iter().enumerate() {
        let pa = &mesh.nodes[a].position;
        let pb = &mesh.nodes[b].position;
        let mx = (pa[0] + pb[0]) * 0.5 - point[0];
        let my = (pa[1] + pb[1]) * 0.5 - point[1];
        let mz = (pa[2] + pb[2]) * 0.5 - point[2];
        let d2 = mx * mx + my * my + mz * mz;
        if d2 < best_dist2 {
            best_dist2 = d2;
            best_idx = i;
        }
    }
    Some(best_idx)
}

// ---------------------------------------------------------------------------
// Crack tip tracking
// ---------------------------------------------------------------------------

/// Tracker that monitors multiple crack tips and their velocities.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrackTracker {
    /// All active crack paths.
    pub cracks: Vec<CrackPath>,
    /// Total energy dissipated by crack propagation (J).
    pub total_energy_dissipated: f64,
    /// Energy release rate G (J/m^2) for tracking.
    pub energy_release_rate: f64,
}

#[allow(dead_code)]
impl CrackTracker {
    /// Create a new tracker with no cracks.
    pub fn new() -> Self {
        Self {
            cracks: Vec::new(),
            total_energy_dissipated: 0.0,
            energy_release_rate: 0.0,
        }
    }

    /// Add a new crack starting at `position` propagating in `direction`.
    pub fn add_crack(&mut self, position: [f64; 3], direction: [f64; 3]) {
        self.cracks.push(CrackPath::new(position, direction));
    }

    /// Number of active cracks.
    pub fn crack_count(&self) -> usize {
        self.cracks.len()
    }

    /// Total path length across all cracks.
    pub fn total_crack_length(&self) -> f64 {
        self.cracks.iter().map(|c| c.path_length()).sum()
    }

    /// Number of active tips across all cracks.
    pub fn total_tip_count(&self) -> usize {
        self.cracks.iter().map(|c| c.tips.len()).sum()
    }
}

impl Default for CrackTracker {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Crack branching
// ---------------------------------------------------------------------------

impl CrackPath {
    /// Branch the primary crack tip into two tips at an angle `branch_angle`
    /// (radians) from the current propagation direction.
    ///
    /// After branching, `tips` will contain 2 new tips replacing the primary.
    /// The new tips diverge symmetrically in the XY plane.
    #[allow(dead_code)]
    pub fn branch(&mut self, branch_angle: f64) {
        if let Some(tip) = self.tips.first().cloned() {
            // Compute two new directions by rotating the propagation direction
            let cos_a = branch_angle.cos();
            let sin_a = branch_angle.sin();
            let d = tip.direction;

            // Rotate in XY plane (simplified 2D branching)
            let dir1 = [
                d[0] * cos_a - d[1] * sin_a,
                d[0] * sin_a + d[1] * cos_a,
                d[2],
            ];
            let dir2 = [
                d[0] * cos_a + d[1] * sin_a,
                -d[0] * sin_a + d[1] * cos_a,
                d[2],
            ];

            // Normalize
            let len1 = (dir1[0] * dir1[0] + dir1[1] * dir1[1] + dir1[2] * dir1[2]).sqrt();
            let len2 = (dir2[0] * dir2[0] + dir2[1] * dir2[1] + dir2[2] * dir2[2]).sqrt();

            let norm_dir1 = if len1 > 1e-15 {
                [dir1[0] / len1, dir1[1] / len1, dir1[2] / len1]
            } else {
                d
            };
            let norm_dir2 = if len2 > 1e-15 {
                [dir2[0] / len2, dir2[1] / len2, dir2[2] / len2]
            } else {
                d
            };

            let mut tip1 = CrackTip::new(tip.position, norm_dir1);
            tip1.stress_intensity = tip.stress_intensity;
            let mut tip2 = CrackTip::new(tip.position, norm_dir2);
            tip2.stress_intensity = tip.stress_intensity;

            self.tips.clear();
            self.tips.push(tip1);
            self.tips.push(tip2);
        }
    }

    /// Number of tips (1 for unbranched, 2+ after branching).
    #[allow(dead_code)]
    pub fn tip_count(&self) -> usize {
        self.tips.len()
    }
}

// ---------------------------------------------------------------------------
// Crack arrest
// ---------------------------------------------------------------------------

/// Criteria for crack arrest (stopping).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrackArrestCriteria {
    /// Minimum stress intensity below which the crack arrests (MPa*sqrt(m)).
    pub k_arrest: f64,
    /// Minimum propagation speed below which the crack arrests (m/s).
    pub min_speed: f64,
    /// Maximum crack length before arrest (m).
    pub max_length: f64,
}

#[allow(dead_code)]
impl CrackArrestCriteria {
    /// Create arrest criteria.
    pub fn new(k_arrest: f64, min_speed: f64, max_length: f64) -> Self {
        Self {
            k_arrest,
            min_speed,
            max_length,
        }
    }

    /// Check if a crack should be arrested.
    pub fn should_arrest(&self, k_i: f64, speed: f64, length: f64) -> bool {
        k_i < self.k_arrest || speed < self.min_speed || length > self.max_length
    }
}

/// Check arrest and deactivate crack tips if criteria are met.
#[allow(dead_code)]
pub fn check_arrest(crack: &mut CrackPath, criteria: &CrackArrestCriteria, k_i: f64) -> bool {
    let speed = crack
        .current_tip()
        .map(|t| t.propagation_velocity)
        .unwrap_or(0.0);
    let length = crack.path_length();

    if criteria.should_arrest(k_i, speed, length) {
        // Set all tip velocities to zero
        for tip in &mut crack.tips {
            tip.propagation_velocity = 0.0;
        }
        true
    } else {
        false
    }
}

// ---------------------------------------------------------------------------
// Dynamic fracture with energy tracking
// ---------------------------------------------------------------------------

/// Dynamic fracture state tracking.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DynamicFracture {
    /// Fracture criteria.
    pub criteria: CrackCriteria,
    /// Arrest criteria.
    pub arrest: CrackArrestCriteria,
    /// Critical energy release rate G_c (J/m^2).
    pub g_c: f64,
    /// Accumulated energy dissipated by fracture (J).
    pub energy_dissipated: f64,
    /// History of energy dissipation per step.
    pub energy_history: Vec<f64>,
}

#[allow(dead_code)]
impl DynamicFracture {
    /// Create a new dynamic fracture model.
    pub fn new(criteria: CrackCriteria, g_c: f64) -> Self {
        Self {
            arrest: CrackArrestCriteria::new(criteria.k_ic * 0.5, 0.1, 10.0),
            criteria,
            g_c,
            energy_dissipated: 0.0,
            energy_history: Vec::new(),
        }
    }

    /// Compute energy release rate from stress intensity: G = K_I^2 / E.
    ///
    /// `youngs_modulus` is in Pa.
    pub fn energy_release_rate(&self, k_i: f64, youngs_modulus: f64) -> f64 {
        if youngs_modulus.abs() < 1e-30 {
            return 0.0;
        }
        k_i * k_i / youngs_modulus
    }

    /// Step the dynamic fracture: propagate if G > G_c, track energy.
    ///
    /// Returns the crack advance distance in this step (0 if no propagation).
    pub fn step(
        &mut self,
        crack: &mut CrackPath,
        mesh: &mut SoftBodyMesh,
        k_i: f64,
        youngs_modulus: f64,
        dt: f64,
    ) -> f64 {
        let g = self.energy_release_rate(k_i, youngs_modulus);

        if g < self.g_c {
            self.energy_history.push(0.0);
            return 0.0;
        }

        // Propagate
        let speed = self.criteria.propagation_speed(k_i);
        let advance = speed * dt;
        crack.advance_tip(dt, speed);

        // Split nearest edge
        if let Some(tip_pos) = crack.current_tip().map(|t| t.position)
            && let Some(edge_idx) = nearest_edge_to_point(mesh, tip_pos)
        {
            split_edge(mesh, edge_idx);
        }

        // Energy dissipated = G_c * new_area (simplified as G_c * advance * 1.0)
        let energy = self.g_c * advance;
        self.energy_dissipated += energy;
        self.energy_history.push(energy);

        advance
    }

    /// Total energy dissipated so far.
    pub fn total_energy_dissipated(&self) -> f64 {
        self.energy_dissipated
    }
}

/// Compute per-node displacement magnitudes.
#[allow(dead_code)]
pub fn displacement_magnitudes(displacements: &[[f64; 3]]) -> Vec<f64> {
    displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .collect()
}

/// Compute the maximum displacement magnitude.
#[allow(dead_code)]
pub fn max_displacement(displacements: &[[f64; 3]]) -> f64 {
    displacement_magnitudes(displacements)
        .into_iter()
        .fold(0.0_f64, f64::max)
}

// ---------------------------------------------------------------------------
// Stress intensity criterion (K-field evaluation)
// ---------------------------------------------------------------------------

/// Mode-I, II, III stress intensity factors for linear elastic fracture.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct StressIntensityFactors {
    /// Mode-I (opening) stress intensity factor (Pa·√m).
    pub k1: f64,
    /// Mode-II (sliding) stress intensity factor (Pa·√m).
    pub k2: f64,
    /// Mode-III (tearing) stress intensity factor (Pa·√m).
    pub k3: f64,
}

#[allow(dead_code)]
impl StressIntensityFactors {
    /// Create with only Mode-I loading.
    pub fn mode1(k1: f64) -> Self {
        Self {
            k1,
            k2: 0.0,
            k3: 0.0,
        }
    }

    /// Effective mixed-mode stress intensity (Irwin criterion).
    pub fn effective(&self) -> f64 {
        (self.k1 * self.k1 + self.k2 * self.k2 + self.k3 * self.k3).sqrt()
    }

    /// Maximum circumferential stress criterion: propagation angle (radians)
    /// from the current crack direction for mixed-mode loading.
    ///
    /// Returns the kink angle θ using the max-σ_θθ criterion:
    ///   θ = 2 * atan((K1 - sqrt(K1^2 + 8*K2^2)) / (4*K2))
    pub fn kink_angle(&self) -> f64 {
        if self.k2.abs() < 1e-30 {
            return 0.0;
        }
        let discriminant = self.k1 * self.k1 + 8.0 * self.k2 * self.k2;
        let num = self.k1 - discriminant.sqrt();
        let denom = 4.0 * self.k2;
        2.0 * (num / denom).atan()
    }

    /// Strain energy release rate G = (K1^2 + K2^2) / E + K3^2 / (2*mu)
    /// for plane-strain conditions.
    ///
    /// `youngs_modulus` E in Pa, `poisson_ratio` ν.
    pub fn energy_release_rate(&self, youngs_modulus: f64, poisson_ratio: f64) -> f64 {
        if youngs_modulus.abs() < 1e-30 {
            return 0.0;
        }
        let e_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
        let mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
        (self.k1 * self.k1 + self.k2 * self.k2) / e_prime
            + (if mu.abs() > 1e-30 {
                self.k3 * self.k3 / (2.0 * mu)
            } else {
                0.0
            })
    }
}

// ---------------------------------------------------------------------------
// Cohesive zone element
// ---------------------------------------------------------------------------

/// A cohesive zone element connecting two node pairs.
///
/// Models progressive damage across the interface using a bilinear traction-separation law.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CohesiveElement {
    /// Indices of the top-face node pair `[a, b]`.
    pub top_nodes: [usize; 2],
    /// Indices of the bottom-face node pair `[c, d]`.
    pub bottom_nodes: [usize; 2],
    /// Peak traction (Pa).
    pub t_max: f64,
    /// Critical separation (m) at which damage = 1.
    pub delta_c: f64,
    /// Current damage variable D ∈ \[0, 1\].
    pub damage: f64,
    /// Whether the element has fully fractured (D = 1).
    pub is_fractured: bool,
}

#[allow(dead_code)]
impl CohesiveElement {
    /// Create a new cohesive element.
    pub fn new(top_nodes: [usize; 2], bottom_nodes: [usize; 2], t_max: f64, delta_c: f64) -> Self {
        Self {
            top_nodes,
            bottom_nodes,
            t_max,
            delta_c,
            damage: 0.0,
            is_fractured: false,
        }
    }

    /// Compute opening displacement between mid-points of top and bottom node pairs.
    pub fn opening_displacement(&self, mesh: &SoftBodyMesh) -> f64 {
        let n_len = mesh.nodes.len();
        let [a, b] = self.top_nodes;
        let [c, d] = self.bottom_nodes;
        if a >= n_len || b >= n_len || c >= n_len || d >= n_len {
            return 0.0;
        }
        let top_mid = [
            (mesh.nodes[a].position[0] + mesh.nodes[b].position[0]) * 0.5,
            (mesh.nodes[a].position[1] + mesh.nodes[b].position[1]) * 0.5,
            (mesh.nodes[a].position[2] + mesh.nodes[b].position[2]) * 0.5,
        ];
        let bot_mid = [
            (mesh.nodes[c].position[0] + mesh.nodes[d].position[0]) * 0.5,
            (mesh.nodes[c].position[1] + mesh.nodes[d].position[1]) * 0.5,
            (mesh.nodes[c].position[2] + mesh.nodes[d].position[2]) * 0.5,
        ];
        let diff = [
            top_mid[0] - bot_mid[0],
            top_mid[1] - bot_mid[1],
            top_mid[2] - bot_mid[2],
        ];
        (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt()
    }

    /// Update damage using a bilinear traction-separation law.
    ///
    /// Returns the current traction magnitude (Pa).
    pub fn update_damage(&mut self, delta: f64) -> f64 {
        if self.is_fractured {
            return 0.0;
        }
        let delta_0 = self.delta_c * 0.1; // onset of damage
        if delta <= delta_0 {
            // Linear elastic regime
            return self.t_max * delta / delta_0;
        }
        // Damage growth
        let d_new = (delta - delta_0) / (self.delta_c - delta_0 + 1e-30);
        self.damage = d_new.clamp(0.0, 1.0);
        if self.damage >= 1.0 {
            self.is_fractured = true;
            return 0.0;
        }
        self.t_max * (1.0 - self.damage) * (self.delta_c - delta) / (self.delta_c - delta_0 + 1e-30)
    }
}

/// Insert cohesive zone elements along the edges on both sides of a crack plane.
///
/// For each edge in `edge_indices`, inserts a cohesive element whose top/bottom
/// faces are the edge endpoints. Returns the newly created cohesive elements.
#[allow(dead_code)]
pub fn insert_cohesive_elements(
    mesh: &SoftBodyMesh,
    edge_indices: &[usize],
    t_max: f64,
    delta_c: f64,
) -> Vec<CohesiveElement> {
    edge_indices
        .iter()
        .filter_map(|&ei| {
            if ei >= mesh.edges.len() {
                return None;
            }
            let [a, b] = mesh.edges[ei];
            // Create a degenerate cohesive element where top/bottom share the same two nodes
            // (placeholder until full mesh duplication is implemented)
            Some(CohesiveElement::new([a, b], [a, b], t_max, delta_c))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Crack front tracking
// ---------------------------------------------------------------------------

/// Snap-shot of a crack front at a given simulation time.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrackFrontSnapshot {
    /// Simulation time when this snapshot was taken.
    pub time: f64,
    /// Positions of all crack tips at this time.
    pub tip_positions: Vec<[f64; 3]>,
    /// Corresponding propagation velocities.
    pub tip_velocities: Vec<f64>,
}

/// Tracks the evolution of a crack front over time.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrackFrontTracker {
    /// Time-ordered list of crack front snapshots.
    pub snapshots: Vec<CrackFrontSnapshot>,
    /// Maximum number of snapshots to retain.
    pub max_snapshots: usize,
}

#[allow(dead_code)]
impl CrackFrontTracker {
    /// Create a new tracker with the given snapshot capacity.
    pub fn new(max_snapshots: usize) -> Self {
        Self {
            snapshots: Vec::new(),
            max_snapshots: max_snapshots.max(1),
        }
    }

    /// Record a snapshot from the given crack path.
    pub fn record(&mut self, crack: &CrackPath, time: f64) {
        let tip_positions: Vec<[f64; 3]> = crack.tips.iter().map(|t| t.position).collect();
        let tip_velocities: Vec<f64> = crack.tips.iter().map(|t| t.propagation_velocity).collect();
        self.snapshots.push(CrackFrontSnapshot {
            time,
            tip_positions,
            tip_velocities,
        });
        // Trim oldest if over capacity
        if self.snapshots.len() > self.max_snapshots {
            self.snapshots.remove(0);
        }
    }

    /// Number of recorded snapshots.
    pub fn snapshot_count(&self) -> usize {
        self.snapshots.len()
    }

    /// Maximum crack propagation velocity ever recorded.
    pub fn peak_velocity(&self) -> f64 {
        self.snapshots
            .iter()
            .flat_map(|s| s.tip_velocities.iter())
            .cloned()
            .fold(0.0_f64, f64::max)
    }

    /// Total crack area grown, estimated as `sum of step distances × 1 m` (unit thickness).
    pub fn total_area_grown(&self) -> f64 {
        if self.snapshots.len() < 2 {
            return 0.0;
        }
        let mut area = 0.0_f64;
        for w in self.snapshots.windows(2) {
            let a = &w[0];
            let b = &w[1];
            let n = a.tip_positions.len().min(b.tip_positions.len());
            for i in 0..n {
                let pa = a.tip_positions[i];
                let pb = b.tip_positions[i];
                let dx = pb[0] - pa[0];
                let dy = pb[1] - pa[1];
                let dz = pb[2] - pa[2];
                area += (dx * dx + dy * dy + dz * dz).sqrt();
            }
        }
        area
    }
}

// ---------------------------------------------------------------------------
// Fragment identification
// ---------------------------------------------------------------------------

/// A fragment is a set of node indices that form a connected component
/// after crack propagation has severed some edges.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Node indices belonging to this fragment.
    pub node_indices: Vec<usize>,
    /// Fragment index (assigned after identification).
    pub id: usize,
}

/// Identify connected fragments in a mesh by building a node-adjacency graph
/// from the remaining (unbroken) edges.
///
/// Returns a list of fragments. Each fragment contains the node indices of
/// one connected component.
#[allow(dead_code)]
pub fn identify_fragments(mesh: &SoftBodyMesh) -> Vec<Fragment> {
    let n = mesh.node_count();
    if n == 0 {
        return Vec::new();
    }
    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &[a, b] in &mesh.edges {
        if a < n && b < n {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    // BFS/DFS to find connected components
    let mut visited = vec![false; n];
    let mut fragments = Vec::new();
    let mut frag_id = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        while let Some(node) = stack.pop() {
            if visited[node] {
                continue;
            }
            visited[node] = true;
            component.push(node);
            for &nb in &adj[node] {
                if !visited[nb] {
                    stack.push(nb);
                }
            }
        }
        fragments.push(Fragment {
            node_indices: component,
            id: frag_id,
        });
        frag_id += 1;
    }
    fragments
}

/// Compute the centre of mass of a fragment.
#[allow(dead_code)]
pub fn fragment_centre_of_mass(fragment: &Fragment, mesh: &SoftBodyMesh) -> [f64; 3] {
    if fragment.node_indices.is_empty() {
        return [0.0; 3];
    }
    let mut sum = [0.0_f64; 3];
    let mut total_mass = 0.0_f64;
    for &idx in &fragment.node_indices {
        if idx >= mesh.nodes.len() {
            continue;
        }
        let m = mesh.nodes[idx].mass;
        sum[0] += mesh.nodes[idx].position[0] * m;
        sum[1] += mesh.nodes[idx].position[1] * m;
        sum[2] += mesh.nodes[idx].position[2] * m;
        total_mass += m;
    }
    if total_mass < 1e-30 {
        return [0.0; 3];
    }
    [
        sum[0] / total_mass,
        sum[1] / total_mass,
        sum[2] / total_mass,
    ]
}

// ---------------------------------------------------------------------------
// Branching heuristic
// ---------------------------------------------------------------------------

/// Evaluate whether a crack should branch based on energy and velocity criteria.
///
/// Returns `true` when branching is recommended.
#[allow(dead_code)]
pub fn should_branch(
    k_i: f64,
    k_ic: f64,
    velocity: f64,
    terminal_velocity: f64,
    branch_ratio: f64,
) -> bool {
    // Branch when K exceeds `branch_ratio × K_Ic` AND velocity > 40% of terminal
    k_i > branch_ratio * k_ic && velocity > 0.4 * terminal_velocity
}

/// Compute the branching angle (radians) from the Yoffe model:
///   branch_angle ≈ arccos(1 - 0.14*(v/c)^2) where c is Rayleigh wave speed.
///
/// Clamps to \[0, π/3\].
#[allow(dead_code)]
pub fn branching_angle(velocity: f64, rayleigh_speed: f64) -> f64 {
    if rayleigh_speed < 1e-12 {
        return 0.0;
    }
    let ratio = (velocity / rayleigh_speed).min(1.0);

    (1.0 - 0.14 * ratio * ratio)
        .acos()
        .min(std::f64::consts::FRAC_PI_3)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. CrackCriteria::glass has correct K_Ic.
    #[test]
    fn test_crack_criteria_glass_k_ic() {
        let g = CrackCriteria::glass();
        assert!((g.k_ic - 0.7).abs() < 1e-10);
    }

    // 2. CrackCriteria::steel has correct K_Ic.
    #[test]
    fn test_crack_criteria_steel_k_ic() {
        let s = CrackCriteria::steel();
        assert!((s.k_ic - 50.0).abs() < 1e-10);
    }

    // 3. should_propagate returns false when K_I <= K_Ic.
    #[test]
    fn test_should_not_propagate_below_threshold() {
        let g = CrackCriteria::glass();
        assert!(!g.should_propagate(0.5));
        assert!(!g.should_propagate(0.7));
    }

    // 4. should_propagate returns true when K_I > K_Ic.
    #[test]
    fn test_should_propagate_above_threshold() {
        let g = CrackCriteria::glass();
        assert!(g.should_propagate(1.0));
    }

    // 5. propagation_speed returns 0 at or below K_Ic.
    #[test]
    fn test_propagation_speed_zero_at_threshold() {
        let g = CrackCriteria::glass();
        assert_eq!(g.propagation_speed(0.5), 0.0);
        assert_eq!(g.propagation_speed(0.7), 0.0);
    }

    // 6. propagation_speed is positive above K_Ic.
    #[test]
    fn test_propagation_speed_positive_above_threshold() {
        let g = CrackCriteria::glass();
        let v = g.propagation_speed(1.4); // 2 * K_Ic → v = v_max * 0.5
        assert!((v - g.max_propagation_speed * 0.5).abs() < 1e-6);
    }

    // 7. propagation_speed is clamped to v_max.
    #[test]
    fn test_propagation_speed_clamped() {
        let g = CrackCriteria::glass();
        let v = g.propagation_speed(1e12);
        assert!(v <= g.max_propagation_speed);
    }

    // 8. SoftBodyMesh::new starts empty.
    #[test]
    fn test_soft_body_mesh_new_empty() {
        let m = SoftBodyMesh::new();
        assert_eq!(m.node_count(), 0);
        assert_eq!(m.triangle_count(), 0);
        assert!(m.edges.is_empty());
    }

    // 9. add_node increments node count.
    #[test]
    fn test_add_node_increments_count() {
        let mut m = SoftBodyMesh::new();
        let i0 = m.add_node([0.0, 0.0, 0.0], 1.0);
        let i1 = m.add_node([1.0, 0.0, 0.0], 1.0);
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(m.node_count(), 2);
    }

    // 10. edge_length computes correct distance.
    #[test]
    fn test_edge_length_correct() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([3.0, 4.0, 0.0], 1.0);
        m.add_edge(0, 1);
        assert!((m.edge_length(0) - 5.0).abs() < 1e-10);
    }

    // 11. add_triangle increments triangle count.
    #[test]
    fn test_add_triangle_increments_count() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_node([0.0, 1.0, 0.0], 1.0);
        m.add_triangle(0, 1, 2);
        assert_eq!(m.triangle_count(), 1);
    }

    // 12. CrackPath::new initialises with one tip and one history point.
    #[test]
    fn test_crack_path_new() {
        let start = [1.0, 2.0, 3.0];
        let dir = [1.0, 0.0, 0.0];
        let cp = CrackPath::new(start, dir);
        assert_eq!(cp.tips.len(), 1);
        assert_eq!(cp.history.len(), 1);
    }

    // 13. advance_tip moves the crack tip in the correct direction.
    #[test]
    fn test_advance_tip_moves_position() {
        let mut cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cp.advance_tip(1.0, 2.0);
        let pos = cp.current_tip().unwrap().position;
        assert!((pos[0] - 2.0).abs() < 1e-10);
        assert!(pos[1].abs() < 1e-10);
    }

    // 14. path_length accumulates correctly.
    #[test]
    fn test_path_length_accumulates() {
        let mut cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cp.advance_tip(1.0, 3.0); // +3 in x
        cp.advance_tip(1.0, 4.0); // +4 in x
        let len = cp.path_length();
        assert!((len - 7.0).abs() < 1e-8);
    }

    // 15. path_length is zero for a single-point history.
    #[test]
    fn test_path_length_zero_single_point() {
        let cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(cp.path_length(), 0.0);
    }

    // 16. split_edge increases node count by 1 and edge count by 1.
    #[test]
    fn test_split_edge_increases_counts() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([2.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let before_nodes = m.node_count();
        let before_edges = m.edges.len();
        split_edge(&mut m, 0);
        assert_eq!(m.node_count(), before_nodes + 1);
        assert_eq!(m.edges.len(), before_edges + 1);
    }

    // 17. split_edge places midpoint correctly.
    #[test]
    fn test_split_edge_midpoint_position() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([4.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mid_idx = split_edge(&mut m, 0);
        let mid_pos = m.nodes[mid_idx].position;
        assert!((mid_pos[0] - 2.0).abs() < 1e-10);
        assert!(mid_pos[1].abs() < 1e-10);
    }

    // 18. nearest_edge_to_point returns None for empty mesh.
    #[test]
    fn test_nearest_edge_empty_mesh() {
        let m = SoftBodyMesh::new();
        assert!(nearest_edge_to_point(&m, [0.0, 0.0, 0.0]).is_none());
    }

    // 19. nearest_edge_to_point finds the closest edge.
    #[test]
    fn test_nearest_edge_finds_closest() {
        let mut m = SoftBodyMesh::new();
        // Edge 0: midpoint at (0.5, 0, 0)
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        // Edge 1: midpoint at (10.5, 0, 0)
        m.add_node([10.0, 0.0, 0.0], 1.0);
        m.add_node([11.0, 0.0, 0.0], 1.0);
        m.add_edge(2, 3);
        // Query close to edge 0's midpoint.
        let result = nearest_edge_to_point(&m, [0.6, 0.0, 0.0]);
        assert_eq!(result, Some(0));
    }

    // 20. propagate_crack returns false when K_I < K_Ic.
    #[test]
    fn test_propagate_crack_no_propagation_below_threshold() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mut crack = CrackPath::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let criteria = CrackCriteria::glass();
        let propagated = propagate_crack(&mut m, &mut crack, &criteria, 0.3, 0.01);
        assert!(!propagated);
    }

    // 21. propagate_crack returns true and modifies mesh when K_I > K_Ic.
    #[test]
    fn test_propagate_crack_propagates_above_threshold() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mut crack = CrackPath::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let criteria = CrackCriteria::glass();
        let before_nodes = m.node_count();
        let propagated = propagate_crack(&mut m, &mut crack, &criteria, 2.0, 0.01);
        assert!(propagated);
        assert!(m.node_count() > before_nodes);
    }

    // 22. compute_stress_intensity returns 0 with no nodes near tip.
    #[test]
    fn test_compute_stress_intensity_empty() {
        let m = SoftBodyMesh::new();
        let crack = CrackPath::new([100.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let displacements: Vec<[f64; 3]> = vec![];
        let k = compute_stress_intensity(&crack, &m, &displacements);
        assert_eq!(k, 0.0);
    }

    // 23. compute_stress_intensity increases with larger displacements.
    #[test]
    fn test_compute_stress_intensity_increases_with_displacement() {
        let mut m = SoftBodyMesh::new();
        m.add_node([0.1, 0.0, 0.0], 1.0);
        m.add_node([-0.1, 0.0, 0.0], 1.0);

        let mut crack = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        // Give the crack a nonzero path length.
        crack.advance_tip(1.0, 0.5);

        let small_disp: Vec<[f64; 3]> = vec![[0.001, 0.0, 0.0], [0.001, 0.0, 0.0]];
        let large_disp: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

        let k_small = compute_stress_intensity(&crack, &m, &small_disp);
        let k_large = compute_stress_intensity(&crack, &m, &large_disp);
        assert!(k_large > k_small, "k_large={k_large}, k_small={k_small}");
    }

    // 24. CrackTip stores correct initial values.
    #[test]
    fn test_crack_tip_initial_values() {
        let ct = CrackTip::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0]);
        assert_eq!(ct.stress_intensity, 0.0);
        assert_eq!(ct.propagation_velocity, 0.0);
        assert_eq!(ct.position, [1.0, 2.0, 3.0]);
    }

    // 25. Default trait for SoftBodyMesh works identically to new().
    #[test]
    fn test_soft_body_mesh_default() {
        let m: SoftBodyMesh = Default::default();
        assert_eq!(m.node_count(), 0);
    }

    // 26. CrackTracker starts empty.
    #[test]
    fn test_crack_tracker_new() {
        let t = CrackTracker::new();
        assert_eq!(t.crack_count(), 0);
        assert_eq!(t.total_tip_count(), 0);
        assert!(t.total_crack_length() < 1e-15);
    }

    // 27. CrackTracker add_crack increments count.
    #[test]
    fn test_crack_tracker_add_crack() {
        let mut t = CrackTracker::new();
        t.add_crack([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        t.add_crack([5.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(t.crack_count(), 2);
        assert_eq!(t.total_tip_count(), 2);
    }

    // 28. CrackTracker total_crack_length after advancing.
    #[test]
    fn test_crack_tracker_total_length() {
        let mut t = CrackTracker::new();
        t.add_crack([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        t.cracks[0].advance_tip(1.0, 5.0); // 5m advance
        assert!((t.total_crack_length() - 5.0).abs() < 1e-8);
    }

    // 29. Crack branching creates two tips.
    #[test]
    fn test_crack_branching() {
        let mut cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(cp.tip_count(), 1);
        cp.branch(0.3);
        assert_eq!(cp.tip_count(), 2);
    }

    // 30. Branched tips have different directions.
    #[test]
    fn test_branched_tips_diverge() {
        let mut cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cp.branch(0.5);
        let d1 = cp.tips[0].direction;
        let d2 = cp.tips[1].direction;
        // Directions should differ
        let diff = (d1[0] - d2[0]).abs() + (d1[1] - d2[1]).abs();
        assert!(diff > 0.01, "Branched directions should differ");
    }

    // 31. CrackArrestCriteria should_arrest below k_arrest.
    #[test]
    fn test_arrest_below_k() {
        let ac = CrackArrestCriteria::new(1.0, 0.1, 100.0);
        assert!(ac.should_arrest(0.5, 10.0, 1.0));
    }

    // 32. CrackArrestCriteria should not arrest above all thresholds.
    #[test]
    fn test_no_arrest_above_thresholds() {
        let ac = CrackArrestCriteria::new(1.0, 0.1, 100.0);
        assert!(!ac.should_arrest(2.0, 10.0, 5.0));
    }

    // 33. CrackArrestCriteria arrest on max_length.
    #[test]
    fn test_arrest_max_length() {
        let ac = CrackArrestCriteria::new(0.1, 0.01, 10.0);
        assert!(ac.should_arrest(5.0, 100.0, 15.0));
    }

    // 34. check_arrest sets velocity to zero.
    #[test]
    fn test_check_arrest_zeros_velocity() {
        let mut cp = CrackPath::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        cp.advance_tip(1.0, 100.0);
        let ac = CrackArrestCriteria::new(10.0, 0.1, 100.0);
        let arrested = check_arrest(&mut cp, &ac, 0.5); // k_i < k_arrest
        assert!(arrested);
        assert_eq!(cp.tips[0].propagation_velocity, 0.0);
    }

    // 35. DynamicFracture energy_release_rate formula.
    #[test]
    fn test_energy_release_rate() {
        let criteria = CrackCriteria::glass();
        let df = DynamicFracture::new(criteria, 10.0);
        let g = df.energy_release_rate(1.0, 70e9); // K_I = 1 MPa*sqrt(m), E = 70 GPa
        let expected = 1.0 / 70e9;
        assert!((g - expected).abs() < 1e-20);
    }

    // 36. DynamicFracture step tracks energy.
    #[test]
    fn test_dynamic_fracture_step() {
        let criteria = CrackCriteria::glass();
        let mut df = DynamicFracture::new(criteria, 1.0); // G_c = 1.0
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mut crack = CrackPath::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
        // Need G = K_I^2 / E > G_c = 1.0, so K_I > sqrt(E * G_c)
        // For E = 70e9: K_I > sqrt(70e9) ≈ 264575, use 300000
        let advance = df.step(&mut crack, &mut m, 3e5, 70e9, 0.001);
        assert!(advance > 0.0, "Should propagate");
        assert!(df.total_energy_dissipated() > 0.0, "Should track energy");
    }

    // 37. DynamicFracture no propagation below G_c.
    #[test]
    fn test_dynamic_fracture_no_propagation() {
        let criteria = CrackCriteria::steel();
        let mut df = DynamicFracture::new(criteria, 1e10); // Very high G_c
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mut crack = CrackPath::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let advance = df.step(&mut crack, &mut m, 1.0, 200e9, 0.001);
        assert_eq!(advance, 0.0);
    }

    // 38. displacement_magnitudes.
    #[test]
    fn test_displacement_magnitudes() {
        let disps = vec![[3.0, 4.0, 0.0], [0.0, 0.0, 1.0]];
        let mags = displacement_magnitudes(&disps);
        assert!((mags[0] - 5.0).abs() < 1e-10);
        assert!((mags[1] - 1.0).abs() < 1e-10);
    }

    // 39. max_displacement.
    #[test]
    fn test_max_displacement() {
        let disps = vec![[3.0, 4.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((max_displacement(&disps) - 5.0).abs() < 1e-10);
    }

    // 40. CrackTracker default.
    #[test]
    fn test_crack_tracker_default() {
        let t: CrackTracker = Default::default();
        assert_eq!(t.crack_count(), 0);
    }

    // 41. ParticleNode velocity starts at zero.
    #[test]
    fn test_particle_node_velocity_zero() {
        let n = ParticleNode::new([1.0, 2.0, 3.0], 1.0, 0);
        assert_eq!(n.velocity, [0.0, 0.0, 0.0]);
    }

    // 42. DynamicFracture energy history tracks steps.
    #[test]
    fn test_dynamic_fracture_energy_history() {
        let criteria = CrackCriteria::glass();
        let mut df = DynamicFracture::new(criteria, 1.0);
        let mut m = SoftBodyMesh::new();
        m.add_node([0.0, 0.0, 0.0], 1.0);
        m.add_node([1.0, 0.0, 0.0], 1.0);
        m.add_edge(0, 1);
        let mut crack = CrackPath::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]);
        df.step(&mut crack, &mut m, 2.0, 70e9, 0.001);
        df.step(&mut crack, &mut m, 0.1, 70e9, 0.001); // Below threshold
        assert_eq!(df.energy_history.len(), 2);
    }

    // --- StressIntensityFactors ---

    // 43. Mode-I only effective equals K1.
    #[test]
    fn test_sif_mode1_effective() {
        let sif = StressIntensityFactors::mode1(2.5);
        assert!((sif.effective() - 2.5).abs() < 1e-12);
    }

    // 44. Mixed mode effective is RMS.
    #[test]
    fn test_sif_mixed_mode_effective() {
        let sif = StressIntensityFactors {
            k1: 3.0,
            k2: 4.0,
            k3: 0.0,
        };
        assert!((sif.effective() - 5.0).abs() < 1e-10);
    }

    // 45. Kink angle is zero for pure Mode-I.
    #[test]
    fn test_sif_kink_angle_mode1() {
        let sif = StressIntensityFactors::mode1(1.0);
        assert!(sif.kink_angle().abs() < 1e-12);
    }

    // 46. Kink angle is nonzero for mixed mode.
    #[test]
    fn test_sif_kink_angle_mixed() {
        let sif = StressIntensityFactors {
            k1: 1.0,
            k2: 0.5,
            k3: 0.0,
        };
        let angle = sif.kink_angle();
        assert!(
            angle.abs() > 0.01,
            "mixed-mode kink angle should be nonzero"
        );
    }

    // 47. Energy release rate is zero for zero K.
    #[test]
    fn test_sif_energy_release_rate_zero() {
        let sif = StressIntensityFactors {
            k1: 0.0,
            k2: 0.0,
            k3: 0.0,
        };
        assert!((sif.energy_release_rate(70e9, 0.3)).abs() < 1e-30);
    }

    // 48. Energy release rate is positive for nonzero K.
    #[test]
    fn test_sif_energy_release_rate_positive() {
        let sif = StressIntensityFactors::mode1(1e6);
        let g = sif.energy_release_rate(200e9, 0.3);
        assert!(g > 0.0, "G should be positive, got {g}");
    }

    // --- CohesiveElement ---

    // 49. CohesiveElement starts with zero damage.
    #[test]
    fn test_cohesive_element_initial_damage() {
        let ce = CohesiveElement::new([0, 1], [2, 3], 1e6, 1e-4);
        assert!((ce.damage).abs() < 1e-12);
        assert!(!ce.is_fractured);
    }

    // 50. CohesiveElement update_damage in elastic regime.
    #[test]
    fn test_cohesive_element_elastic_regime() {
        let mut ce = CohesiveElement::new([0, 1], [2, 3], 1e6, 1e-4);
        let t = ce.update_damage(0.0); // zero opening
        assert!(t.abs() < 1e-12, "zero displacement → zero traction");
        assert!(!ce.is_fractured);
    }

    // 51. CohesiveElement fully fractures at delta_c.
    #[test]
    fn test_cohesive_element_fracture() {
        let mut ce = CohesiveElement::new([0, 1], [2, 3], 1e6, 1e-4);
        ce.update_damage(2e-4); // beyond delta_c
        assert!(ce.is_fractured, "element should fracture beyond delta_c");
    }

    // 52. insert_cohesive_elements returns correct count.
    #[test]
    fn test_insert_cohesive_elements_count() {
        let mut mesh = SoftBodyMesh::new();
        mesh.add_node([0.0; 3], 1.0);
        mesh.add_node([1.0, 0.0, 0.0], 1.0);
        mesh.add_node([2.0, 0.0, 0.0], 1.0);
        mesh.add_edge(0, 1);
        mesh.add_edge(1, 2);
        let cohesive = insert_cohesive_elements(&mesh, &[0, 1], 1e6, 1e-4);
        assert_eq!(cohesive.len(), 2);
    }

    // --- CrackFrontTracker ---

    // 53. CrackFrontTracker starts empty.
    #[test]
    fn test_crack_front_tracker_empty() {
        let t = CrackFrontTracker::new(10);
        assert_eq!(t.snapshot_count(), 0);
        assert!((t.peak_velocity()).abs() < 1e-12);
    }

    // 54. Recording a snapshot increments count.
    #[test]
    fn test_crack_front_tracker_record() {
        let mut t = CrackFrontTracker::new(10);
        let crack = CrackPath::new([0.0; 3], [1.0, 0.0, 0.0]);
        t.record(&crack, 0.1);
        assert_eq!(t.snapshot_count(), 1);
    }

    // 55. Tracker respects max_snapshots capacity.
    #[test]
    fn test_crack_front_tracker_cap() {
        let mut t = CrackFrontTracker::new(3);
        let crack = CrackPath::new([0.0; 3], [1.0, 0.0, 0.0]);
        for i in 0..10 {
            t.record(&crack, i as f64 * 0.1);
        }
        assert_eq!(t.snapshot_count(), 3, "should cap at max_snapshots");
    }

    // 56. Peak velocity is tracked.
    #[test]
    fn test_crack_front_tracker_peak_velocity() {
        let mut t = CrackFrontTracker::new(10);
        let mut crack = CrackPath::new([0.0; 3], [1.0, 0.0, 0.0]);
        crack.advance_tip(0.1, 100.0);
        t.record(&crack, 0.1);
        assert!((t.peak_velocity() - 100.0).abs() < 1e-6);
    }

    // 57. total_area_grown with two snapshots.
    #[test]
    fn test_crack_front_tracker_area_grown() {
        let mut t = CrackFrontTracker::new(10);
        let mut crack = CrackPath::new([0.0; 3], [1.0, 0.0, 0.0]);
        t.record(&crack, 0.0);
        crack.advance_tip(1.0, 5.0);
        t.record(&crack, 1.0);
        let area = t.total_area_grown();
        assert!(area > 0.0, "area grown should be positive, got {area}");
    }

    // --- Fragment identification ---

    // 58. Single-node mesh gives one fragment.
    #[test]
    fn test_fragment_single_node() {
        let mut mesh = SoftBodyMesh::new();
        mesh.add_node([0.0; 3], 1.0);
        let frags = identify_fragments(&mesh);
        assert_eq!(frags.len(), 1);
        assert_eq!(frags[0].node_indices.len(), 1);
    }

    // 59. Connected mesh gives one fragment.
    #[test]
    fn test_fragment_connected_mesh() {
        let mut mesh = SoftBodyMesh::new();
        mesh.add_node([0.0; 3], 1.0);
        mesh.add_node([1.0, 0.0, 0.0], 1.0);
        mesh.add_node([2.0, 0.0, 0.0], 1.0);
        mesh.add_edge(0, 1);
        mesh.add_edge(1, 2);
        let frags = identify_fragments(&mesh);
        assert_eq!(frags.len(), 1);
    }

    // 60. Disconnected mesh gives multiple fragments.
    #[test]
    fn test_fragment_disconnected_mesh() {
        let mut mesh = SoftBodyMesh::new();
        mesh.add_node([0.0; 3], 1.0); // 0
        mesh.add_node([1.0, 0.0, 0.0], 1.0); // 1
        mesh.add_node([10.0, 0.0, 0.0], 1.0); // 2 (disconnected)
        mesh.add_edge(0, 1);
        // No edge connecting 2 to rest
        let frags = identify_fragments(&mesh);
        assert_eq!(frags.len(), 2, "should have 2 disconnected fragments");
    }

    // 61. Empty mesh gives no fragments.
    #[test]
    fn test_fragment_empty_mesh() {
        let mesh = SoftBodyMesh::new();
        let frags = identify_fragments(&mesh);
        assert!(frags.is_empty());
    }

    // 62. Fragment centre of mass.
    #[test]
    fn test_fragment_centre_of_mass() {
        let mut mesh = SoftBodyMesh::new();
        mesh.add_node([0.0; 3], 1.0);
        mesh.add_node([2.0, 0.0, 0.0], 1.0);
        mesh.add_edge(0, 1);
        let frags = identify_fragments(&mesh);
        assert_eq!(frags.len(), 1);
        let com = fragment_centre_of_mass(&frags[0], &mesh);
        assert!((com[0] - 1.0).abs() < 1e-10, "CoM should be at midpoint");
    }

    // --- Branching heuristic ---

    // 63. No branching below thresholds.
    #[test]
    fn test_should_branch_below_threshold() {
        assert!(!should_branch(1.0, 1.0, 100.0, 2000.0, 2.0));
    }

    // 64. Branching above both thresholds.
    #[test]
    fn test_should_branch_above_threshold() {
        assert!(should_branch(3.0, 1.0, 1000.0, 2000.0, 2.0));
    }

    // 65. Branching angle is zero at zero velocity.
    #[test]
    fn test_branching_angle_zero_velocity() {
        let angle = branching_angle(0.0, 3000.0);
        // At v=0, cos(angle) = 1 - 0.14*0 = 1.0 → acos(1.0) = 0.0
        assert!(
            angle.abs() < 1e-10,
            "angle at zero velocity should be 0, got {angle}"
        );
    }

    // 66. Branching angle increases with velocity.
    #[test]
    fn test_branching_angle_increases_with_velocity() {
        let a1 = branching_angle(500.0, 3000.0);
        let a2 = branching_angle(1000.0, 3000.0);
        assert!(a2 >= a1, "angle should increase with velocity");
    }

    // 67. Branching angle is capped at PI/3.
    #[test]
    fn test_branching_angle_capped() {
        let angle = branching_angle(1e10, 3000.0);
        assert!(angle <= std::f64::consts::FRAC_PI_3 + 1e-10);
    }
}
