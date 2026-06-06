//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{hoop_stress_angle, make_perp_fracture, make_perpendicular, norm3};

/// Result of splitting a node at a crack path.
pub struct SplitResult {
    /// Index of the original node.
    pub original_node: usize,
    /// Index of the newly created duplicate node.
    pub new_node: usize,
    /// Edges that were reassigned to use the new node.
    pub reassigned_edges: Vec<usize>,
}
/// A single vertex on the crack front, plus the local normal to the crack plane.
#[derive(Debug, Clone, Copy)]
pub struct CrackFrontVertex {
    /// Position of this front vertex in world space.
    pub position: [f64; 3],
    /// Unit normal to the crack plane at this vertex.
    pub normal: [f64; 3],
    /// Current Mode-I stress intensity factor at this vertex.
    pub k1: f64,
    /// Current Mode-II stress intensity factor at this vertex.
    pub k2: f64,
    /// Whether this vertex has been arrested.
    pub arrested: bool,
}
impl CrackFrontVertex {
    /// Create a new front vertex.
    pub fn new(position: [f64; 3], normal: [f64; 3], k1: f64, k2: f64) -> Self {
        let n = norm3(normal);
        let normal = if n > 1e-12 {
            [normal[0] / n, normal[1] / n, normal[2] / n]
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            position,
            normal,
            k1,
            k2,
            arrested: false,
        }
    }
    /// Maximum hoop stress criterion: propagation direction angle θ_c.
    ///
    /// From Williams (1957): `θ_c = 2 * atan(K_I/(4*K_II) - sign(K_II)*sqrt((K_I/(4*K_II))^2 + 0.5))`
    ///
    /// Returns 0 when K_II = 0 (pure Mode-I).
    pub fn max_hoop_stress_angle(&self) -> f64 {
        if self.k2.abs() < 1e-30 {
            return 0.0;
        }
        let ratio = self.k1 / (4.0 * self.k2);
        let sign = if self.k2 >= 0.0 { 1.0_f64 } else { -1.0_f64 };
        let inner = ratio * ratio + 0.5;
        2.0 * (ratio - sign * inner.sqrt()).atan()
    }
    /// Equivalent (effective) stress intensity factor:
    /// `K_eq = cos(θ/2) * (K_I * cos²(θ/2) - 1.5 * K_II * sin(θ))`
    pub fn equivalent_sif(&self) -> f64 {
        let theta = self.max_hoop_stress_angle();
        let ch = (theta / 2.0).cos();
        let s = theta.sin();
        ch * (self.k1 * ch * ch - 1.5 * self.k2 * s)
    }
}
/// Tracks the full crack front: a closed or open loop of [`CrackFrontVertex`] entries.
#[derive(Debug, Clone)]
pub struct CrackFront {
    /// Ordered list of vertices forming the crack front.
    pub vertices: Vec<CrackFrontVertex>,
    /// Fracture toughness K_Ic (Pa√m).
    pub k1c: f64,
    /// Branching threshold: when K_eq > k_branch the front may branch.
    pub k_branch: f64,
    /// Minimum propagation increment (m).
    pub da_min: f64,
    /// Maximum propagation increment (m).
    pub da_max: f64,
    /// Indices of tetrahedral elements that have been deleted due to fracture.
    pub deleted_elements: Vec<usize>,
}
impl CrackFront {
    /// Create a new crack front.
    pub fn new(k1c: f64, k_branch: f64, da_min: f64, da_max: f64) -> Self {
        Self {
            vertices: Vec::new(),
            k1c,
            k_branch,
            da_min,
            da_max,
            deleted_elements: Vec::new(),
        }
    }
    /// Add a vertex to the crack front.
    pub fn add_vertex(&mut self, v: CrackFrontVertex) {
        self.vertices.push(v);
    }
    /// Estimate stress intensity factors at each front vertex from a nodal
    /// stress field.
    ///
    /// `stress_voigt[i]` is the stress tensor at node `i` in Voigt notation
    /// `[s11, s22, s33, s12, s23, s13]`.
    /// `node_positions[i]` is the 3-D position of node `i`.
    ///
    /// For each front vertex the method finds the nearest stress sample and
    /// projects the stress onto the local crack-front coordinate frame to
    /// obtain K_I and K_II via the near-tip asymptotic formula
    /// `K = σ * sqrt(2π r)` where `r` is the distance to the nearest sample.
    pub fn estimate_sif_from_stress_field(
        &mut self,
        stress_voigt: &[[f64; 6]],
        node_positions: &[[f64; 3]],
    ) {
        for fv in &mut self.vertices {
            if fv.arrested {
                continue;
            }
            let mut best_dist = f64::MAX;
            let mut best_idx = 0;
            for (i, &pos) in node_positions.iter().enumerate() {
                let dx = pos[0] - fv.position[0];
                let dy = pos[1] - fv.position[1];
                let dz = pos[2] - fv.position[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < best_dist {
                    best_dist = d;
                    best_idx = i;
                }
            }
            if best_dist < 1e-30 || best_idx >= stress_voigt.len() {
                continue;
            }
            let s = stress_voigt[best_idx];
            let n = fv.normal;
            let s_nn = s[0] * n[0] * n[0]
                + s[1] * n[1] * n[1]
                + s[2] * n[2] * n[2]
                + 2.0 * s[3] * n[0] * n[1]
                + 2.0 * s[4] * n[1] * n[2]
                + 2.0 * s[5] * n[0] * n[2];
            let s_shear = (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]).sqrt();
            let r = best_dist.max(1e-10);
            let sq = (2.0 * std::f64::consts::PI * r).sqrt();
            fv.k1 = s_nn * sq;
            fv.k2 = s_shear * sq;
        }
    }
    /// Propagate each active front vertex in the direction given by the
    /// maximum hoop stress criterion.
    ///
    /// The propagation increment is clamped to `[da_min, da_max]` and scales
    /// as `da = da_max * (K_eq / K_Ic - 1)` when `K_eq > K_Ic`.
    ///
    /// Returns the number of vertices that propagated.
    pub fn propagate(&mut self) -> usize {
        let mut count = 0;
        for fv in &mut self.vertices {
            if fv.arrested {
                continue;
            }
            let k_eq = fv.equivalent_sif();
            if k_eq <= self.k1c {
                continue;
            }
            let da = (self.da_max * (k_eq / self.k1c - 1.0)).clamp(self.da_min, self.da_max);
            let theta = fv.max_hoop_stress_angle();
            let n = fv.normal;
            let t = make_perpendicular(n);
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            let dir = [
                cos_t * n[0] + sin_t * t[0],
                cos_t * n[1] + sin_t * t[1],
                cos_t * n[2] + sin_t * t[2],
            ];
            fv.position[0] += da * dir[0];
            fv.position[1] += da * dir[1];
            fv.position[2] += da * dir[2];
            let dl = norm3(dir);
            if dl > 1e-12 {
                fv.normal = [dir[0] / dl, dir[1] / dl, dir[2] / dl];
            }
            count += 1;
        }
        count
    }
    /// Check for branching: vertices with `K_eq > k_branch` spawn a new
    /// branched vertex at ±45° of the current propagation direction.
    ///
    /// Returns the newly created branch vertices (not yet inserted).
    pub fn check_branching(&self) -> Vec<CrackFrontVertex> {
        let mut branches = Vec::new();
        for fv in &self.vertices {
            if fv.arrested {
                continue;
            }
            let k_eq = fv.equivalent_sif();
            if k_eq <= self.k_branch {
                continue;
            }
            let n = fv.normal;
            let t = make_perpendicular(n);
            for sign in [1.0_f64, -1.0] {
                let angle: f64 = sign * std::f64::consts::FRAC_PI_4;
                let c = angle.cos();
                let s = angle.sin();
                let dir = [
                    c * n[0] + s * t[0],
                    c * n[1] + s * t[1],
                    c * n[2] + s * t[2],
                ];
                let dl = norm3(dir);
                let branch_normal = if dl > 1e-12 {
                    [dir[0] / dl, dir[1] / dl, dir[2] / dl]
                } else {
                    n
                };
                branches.push(CrackFrontVertex {
                    position: fv.position,
                    normal: branch_normal,
                    k1: fv.k1 * 0.5,
                    k2: fv.k2 * 0.5,
                    arrested: false,
                });
            }
        }
        branches
    }
    /// Insert branch vertices into the front (call after `check_branching`).
    pub fn apply_branches(&mut self, branches: Vec<CrackFrontVertex>) {
        self.vertices.extend(branches);
    }
    /// Element deletion: mark tetrahedral elements as deleted when any of
    /// their four nodes falls within `deletion_radius` of any front vertex.
    ///
    /// `tet_indices[e]` is the `[i0, i1, i2, i3]` node-index tuple for element `e`.
    /// `positions[i]` is the world position of node `i`.
    ///
    /// Returns the number of newly deleted elements.
    pub fn delete_fractured_elements(
        &mut self,
        tet_indices: &[[usize; 4]],
        positions: &[[f64; 3]],
        deletion_radius: f64,
    ) -> usize {
        let r2 = deletion_radius * deletion_radius;
        let mut count = 0;
        for (e, tet) in tet_indices.iter().enumerate() {
            if self.deleted_elements.contains(&e) {
                continue;
            }
            let mut should_delete = false;
            'outer: for &node_idx in tet.iter() {
                if node_idx >= positions.len() {
                    continue;
                }
                let p = positions[node_idx];
                for fv in &self.vertices {
                    let dx = p[0] - fv.position[0];
                    let dy = p[1] - fv.position[1];
                    let dz = p[2] - fv.position[2];
                    if dx * dx + dy * dy + dz * dz <= r2 {
                        should_delete = true;
                        break 'outer;
                    }
                }
            }
            if should_delete {
                self.deleted_elements.push(e);
                count += 1;
            }
        }
        count
    }
    /// Returns `true` if element `e` has been deleted.
    pub fn is_element_deleted(&self, e: usize) -> bool {
        self.deleted_elements.contains(&e)
    }
    /// Number of active (non-arrested) front vertices.
    pub fn active_vertex_count(&self) -> usize {
        self.vertices.iter().filter(|v| !v.arrested).count()
    }
    /// Arrest all vertices whose K_eq is below K_Ic.
    pub fn arrest_sub_critical(&mut self) {
        for fv in &mut self.vertices {
            if fv.equivalent_sif() < self.k1c {
                fv.arrested = true;
            }
        }
    }
}
/// Configuration for crack-branching detection.
///
/// A crack may bifurcate into two branches when the mode-mixity angle or
/// the crack speed exceeds a threshold.
#[derive(Debug, Clone, Copy)]
pub struct BranchingCriterion {
    /// Angle (rad) above which branching occurs.
    ///
    /// Typically `π/8` (22.5°).
    pub angle_threshold: f64,
    /// Minimum K_eq / K_Ic ratio at which branching is possible.
    pub k_ratio_threshold: f64,
    /// Half-angle of the two branches relative to the main crack direction (rad).
    pub branch_half_angle: f64,
}
impl BranchingCriterion {
    /// Construct with default parameters (angle = π/8, k_ratio = 2.0, branch = π/8).
    pub fn default_params() -> Self {
        Self {
            angle_threshold: std::f64::consts::PI / 8.0,
            k_ratio_threshold: 2.0,
            branch_half_angle: std::f64::consts::PI / 8.0,
        }
    }
    /// Determine whether the crack at `tip` should branch given `k1`, `k2`, `k1c`.
    ///
    /// Returns `None` if branching does not occur.
    /// Returns `Some((branch_a, branch_b))` with the two new crack tips otherwise.
    pub fn crack_branching(
        &self,
        tip: &CrackTip3D,
        k1: f64,
        k2: f64,
        k1c: f64,
    ) -> Option<(CrackTip3D, CrackTip3D)> {
        let k_eq = (k1 * k1 + k2 * k2).sqrt();
        if k_eq < self.k_ratio_threshold * k1c {
            return None;
        }
        let mix_angle = hoop_stress_angle(k1, k2).abs();
        if mix_angle < self.angle_threshold {
            return None;
        }
        let d = tip.direction;
        let t = make_perp_fracture(d);
        let alpha = self.branch_half_angle;
        let make_branch_dir = |sign: f64| -> [f64; 3] {
            let c = alpha.cos();
            let s = alpha.sin() * sign;
            let raw = [
                c * d[0] + s * t[0],
                c * d[1] + s * t[1],
                c * d[2] + s * t[2],
            ];
            let n = norm3(raw);
            if n > 1e-12 {
                [raw[0] / n, raw[1] / n, raw[2] / n]
            } else {
                d
            }
        };
        let dir_a = make_branch_dir(1.0);
        let dir_b = make_branch_dir(-1.0);
        let branch_a = CrackTip3D::new(tip.position, dir_a);
        let branch_b = CrackTip3D::new(tip.position, dir_b);
        Some((branch_a, branch_b))
    }
}
/// Fractureable mesh -- collection of nodes plus edges that can break.
pub struct FractureMesh {
    /// Node positions (world space).
    pub positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Node masses.
    pub masses: Vec<f64>,
    /// All edges (may become broken over time).
    pub edges: Vec<FractureEdge>,
    /// Fracture model applied to all edges.
    pub fracture_model: FractureModel,
    /// Which nodes are fixed (pinned boundary conditions).
    pub fixed_nodes: Vec<bool>,
}
impl FractureMesh {
    /// Create an empty mesh with the given fracture model.
    pub fn new(fracture_model: FractureModel) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            masses: Vec::new(),
            edges: Vec::new(),
            fracture_model,
            fixed_nodes: Vec::new(),
        }
    }
    /// Add a node and return its index.
    pub fn add_node(&mut self, pos: [f64; 3], mass: f64) -> usize {
        let idx = self.positions.len();
        self.positions.push(pos);
        self.velocities.push([0.0; 3]);
        self.masses.push(mass);
        self.fixed_nodes.push(false);
        idx
    }
    /// Add an edge between nodes `a` and `b`.
    ///
    /// The rest length is set to the current Euclidean distance between them.
    /// Returns the edge index.
    pub fn add_edge(&mut self, a: usize, b: usize, stiffness: f64) -> usize {
        let pa = self.positions[a];
        let pb = self.positions[b];
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        let rest_length = (dx * dx + dy * dy + dz * dz).sqrt();
        let idx = self.edges.len();
        self.edges
            .push(FractureEdge::new(a, b, rest_length, stiffness));
        idx
    }
    /// Fix (pin) a node so it is unaffected by forces.
    pub fn pin_node(&mut self, idx: usize) {
        self.fixed_nodes[idx] = true;
    }
    /// Compute the net spring force on each node from all non-broken edges.
    pub fn compute_forces(&self) -> Vec<[f64; 3]> {
        let n = self.positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        for edge in &self.edges {
            let (fa, fb) = edge.spring_force(&self.positions);
            let a = edge.node_a;
            let b = edge.node_b;
            forces[a][0] += fa[0];
            forces[a][1] += fa[1];
            forces[a][2] += fa[2];
            forces[b][0] += fb[0];
            forces[b][1] += fb[1];
            forces[b][2] += fb[2];
        }
        forces
    }
    /// Check every edge against the fracture model.
    ///
    /// Returns the indices of edges that were newly broken in this call.
    pub fn check_fractures(&mut self) -> Vec<usize> {
        let mut newly_broken = Vec::new();
        let positions = &self.positions.clone();
        for (i, edge) in self.edges.iter_mut().enumerate() {
            if edge.check_fracture(positions, &self.fracture_model) {
                newly_broken.push(i);
            }
        }
        newly_broken
    }
    /// Total elastic energy stored in all intact edges.
    pub fn total_elastic_energy(&self) -> f64 {
        self.edges
            .iter()
            .map(|e| e.elastic_energy(&self.positions))
            .sum()
    }
    /// Total kinetic energy of all non-fixed nodes.
    pub fn total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.positions.len() {
            if self.fixed_nodes[i] {
                continue;
            }
            let v = &self.velocities[i];
            let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            ke += 0.5 * self.masses[i] * v_sq;
        }
        ke
    }
    /// Maximum strain across all intact edges.
    pub fn max_strain(&self) -> f64 {
        self.edges
            .iter()
            .filter(|e| !e.broken)
            .map(|e| e.current_strain(&self.positions))
            .fold(0.0_f64, f64::max)
    }
    /// Semi-implicit Euler integration step.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let forces = self.compute_forces();
        let n = self.positions.len();
        for (i, (vel, (force, (mass, fixed)))) in self
            .velocities
            .iter_mut()
            .zip(
                forces
                    .iter()
                    .zip(self.masses.iter().zip(self.fixed_nodes.iter())),
            )
            .enumerate()
        {
            let _ = i;
            if *fixed {
                continue;
            }
            let inv_m = 1.0 / mass;
            for (v, (f, g)) in vel.iter_mut().zip(force.iter().zip(gravity.iter())) {
                *v += dt * (f * inv_m + g);
            }
        }
        for i in 0..n {
            if self.fixed_nodes[i] {
                continue;
            }
            self.positions[i][0] += dt * self.velocities[i][0];
            self.positions[i][1] += dt * self.velocities[i][1];
            self.positions[i][2] += dt * self.velocities[i][2];
        }
        self.check_fractures();
    }
    /// Step with velocity damping.
    pub fn step_damped(&mut self, dt: f64, gravity: [f64; 3], damping: f64) {
        let forces = self.compute_forces();
        let n = self.positions.len();
        for (i, (vel, (force, (mass, fixed)))) in self
            .velocities
            .iter_mut()
            .zip(
                forces
                    .iter()
                    .zip(self.masses.iter().zip(self.fixed_nodes.iter())),
            )
            .enumerate()
        {
            let _ = i;
            if *fixed {
                continue;
            }
            let inv_m = 1.0 / mass;
            for (v, (f, g)) in vel.iter_mut().zip(force.iter().zip(gravity.iter())) {
                *v += dt * (f * inv_m + g);
                *v *= 1.0 - damping;
            }
        }
        for i in 0..n {
            if self.fixed_nodes[i] {
                continue;
            }
            self.positions[i][0] += dt * self.velocities[i][0];
            self.positions[i][1] += dt * self.velocities[i][1];
            self.positions[i][2] += dt * self.velocities[i][2];
        }
        self.check_fractures();
    }
    /// Number of broken edges.
    pub fn n_broken(&self) -> usize {
        self.edges.iter().filter(|e| e.broken).count()
    }
    /// Number of intact (non-broken) edges.
    pub fn active_edge_count(&self) -> usize {
        self.edges.iter().filter(|e| !e.broken).count()
    }
    /// Find connected components of nodes using BFS over non-broken edges.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let n = self.positions.len();
        let mut visited = vec![false; n];
        let mut components: Vec<Vec<usize>> = Vec::new();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for edge in &self.edges {
            if !edge.broken {
                adj[edge.node_a].push(edge.node_b);
                adj[edge.node_b].push(edge.node_a);
            }
        }
        for start in 0..n {
            if visited[start] {
                continue;
            }
            let mut component = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            while let Some(node) = queue.pop_front() {
                component.push(node);
                for &neighbor in &adj[node] {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        queue.push_back(neighbor);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        components
    }
}
/// A 3-D crack tip state used during propagation.
#[derive(Debug, Clone)]
pub struct CrackTip3D {
    /// Current position of the crack tip in 3-D space.
    pub position: [f64; 3],
    /// Current propagation direction (unit vector).
    pub direction: [f64; 3],
    /// Total crack length accumulated so far (m).
    pub length: f64,
    /// Whether the crack has arrested (stopped propagating).
    pub arrested: bool,
}
impl CrackTip3D {
    /// Create a new crack tip at `position` propagating in `direction`.
    pub fn new(position: [f64; 3], direction: [f64; 3]) -> Self {
        let n = norm3(direction);
        let dir = if n > f64::EPSILON {
            [direction[0] / n, direction[1] / n, direction[2] / n]
        } else {
            [1.0, 0.0, 0.0]
        };
        Self {
            position,
            direction: dir,
            length: 0.0,
            arrested: false,
        }
    }
    /// Advance the crack tip by `da` in the current direction.
    pub fn advance(&mut self, da: f64) {
        if self.arrested {
            return;
        }
        self.position[0] += self.direction[0] * da;
        self.position[1] += self.direction[1] * da;
        self.position[2] += self.direction[2] * da;
        self.length += da;
    }
    /// Steer the crack in a new direction (re-normalised).
    pub fn steer(&mut self, new_direction: [f64; 3]) {
        let n = norm3(new_direction);
        if n > f64::EPSILON {
            self.direction = [
                new_direction[0] / n,
                new_direction[1] / n,
                new_direction[2] / n,
            ];
        }
    }
}
/// Which criterion to use when deciding how a crack propagates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropagationCriterion {
    /// Maximum hoop stress criterion (Williams 1957).
    ///
    /// The crack extends in the direction that maximises the circumferential
    /// stress.  Direction angle: `θ_c` from `CrackFrontVertex::max_hoop_stress_angle`.
    MaxHoopStress,
    /// Minimum strain energy density criterion (Sih 1974).
    ///
    /// The crack propagates toward the minimum of the strain energy density S.
    MinStrainEnergy,
    /// Cohesive zone model.
    ///
    /// Propagation is governed by the cohesive traction-separation law;
    /// the crack advances when the crack-opening displacement exceeds the
    /// critical value `delta_c`.
    CohesiveZone,
}
/// A crack represented as a sequence of points.
pub struct CrackPath {
    /// Points along the crack path.
    pub points: Vec<[f64; 3]>,
    /// Current crack tip direction.
    pub direction: [f64; 3],
    /// Current crack half-length.
    pub half_length: f64,
}
impl CrackPath {
    /// Create a new crack at the given origin with initial direction.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        let dir = if len > 1e-12 {
            [direction[0] / len, direction[1] / len, direction[2] / len]
        } else {
            [1.0, 0.0, 0.0]
        };
        Self {
            points: vec![origin],
            direction: dir,
            half_length: 0.0,
        }
    }
    /// Propagate the crack by a given increment in the current direction.
    pub fn propagate(&mut self, increment: f64) {
        let tip = self.tip();
        let new_tip = [
            tip[0] + self.direction[0] * increment,
            tip[1] + self.direction[1] * increment,
            tip[2] + self.direction[2] * increment,
        ];
        self.points.push(new_tip);
        self.half_length += increment;
    }
    /// Propagate the crack in a new direction.
    pub fn propagate_with_direction(&mut self, increment: f64, new_direction: [f64; 3]) {
        let len = (new_direction[0] * new_direction[0]
            + new_direction[1] * new_direction[1]
            + new_direction[2] * new_direction[2])
            .sqrt();
        if len > 1e-12 {
            self.direction = [
                new_direction[0] / len,
                new_direction[1] / len,
                new_direction[2] / len,
            ];
        }
        self.propagate(increment);
    }
    /// Get the current crack tip position.
    pub fn tip(&self) -> [f64; 3] {
        *self.points.last().unwrap_or(&[0.0; 3])
    }
    /// Total crack path length.
    pub fn total_length(&self) -> f64 {
        let mut length = 0.0;
        for i in 1..self.points.len() {
            let dx = self.points[i][0] - self.points[i - 1][0];
            let dy = self.points[i][1] - self.points[i - 1][1];
            let dz = self.points[i][2] - self.points[i - 1][2];
            length += (dx * dx + dy * dy + dz * dz).sqrt();
        }
        length
    }
    /// Number of segments in the crack path.
    pub fn segment_count(&self) -> usize {
        if self.points.len() > 1 {
            self.points.len() - 1
        } else {
            0
        }
    }
}
/// Bilinear traction-separation law for a cohesive zone.
///
/// The model tracks the opening displacement δ of an interface:
/// - For δ ≤ δ₀:  T = (T_c / δ₀) * δ  (linear loading)
/// - For δ₀ < δ < δ_c:  T = T_c * (δ_c - δ) / (δ_c - δ₀)  (softening)
/// - For δ ≥ δ_c:  T = 0  (fully separated)
///
/// Energy of fracture: G_c = 0.5 * T_c * δ_c
#[derive(Debug, Clone)]
pub struct CohesiveZone {
    /// Peak traction T_c (Pa).
    pub peak_traction: f64,
    /// Critical opening at peak traction δ₀ (m).
    pub delta_0: f64,
    /// Critical opening displacement δ_c (m, fully separated).
    pub delta_c: f64,
    /// Maximum historic opening displacement (irreversibility).
    pub delta_max: f64,
    /// Whether this interface has fully separated.
    pub separated: bool,
}
impl CohesiveZone {
    /// Create a new cohesive zone element.
    pub fn new(peak_traction: f64, delta_0: f64, delta_c: f64) -> Self {
        assert!(delta_0 > 0.0 && delta_c > delta_0, "require 0 < δ₀ < δ_c");
        Self {
            peak_traction,
            delta_0,
            delta_c,
            delta_max: 0.0,
            separated: false,
        }
    }
    /// Evaluate the traction for a given opening displacement `delta`.
    ///
    /// Updates the irreversible maximum opening history.
    pub fn traction(&mut self, delta: f64) -> f64 {
        if self.separated {
            return 0.0;
        }
        if delta < 0.0 {
            return 0.0;
        }
        if delta > self.delta_max {
            self.delta_max = delta;
        }
        if self.delta_max >= self.delta_c {
            self.separated = true;
            return 0.0;
        }
        if self.delta_max <= self.delta_0 {
            self.peak_traction * delta / self.delta_0
        } else {
            let d_eff = self.delta_max.min(self.delta_c);
            let t_eff = self.peak_traction * (self.delta_c - d_eff) / (self.delta_c - self.delta_0);
            if self.delta_max > 1e-30 {
                t_eff * delta / d_eff
            } else {
                0.0
            }
        }
    }
    /// Fracture energy per unit area: G_c = 0.5 * T_c * δ_c.
    pub fn fracture_energy(&self) -> f64 {
        0.5 * self.peak_traction * self.delta_c
    }
    /// Current secant stiffness T(δ) / δ.
    pub fn secant_stiffness(&self, delta: f64) -> f64 {
        if delta < 1e-30 || self.separated {
            return 0.0;
        }
        let t = if delta <= self.delta_0 {
            self.peak_traction * delta / self.delta_0
        } else if delta < self.delta_c {
            self.peak_traction * (self.delta_c - delta) / (self.delta_c - self.delta_0)
        } else {
            0.0
        };
        t / delta
    }
    /// Whether the interface has fully separated (δ ≥ δ_c).
    pub fn is_separated(&self) -> bool {
        self.separated
    }
    /// Reset the cohesive zone (for reversible test scenarios).
    pub fn reset(&mut self) {
        self.delta_max = 0.0;
        self.separated = false;
    }
}
/// Which fracture criterion to use.
pub enum FractureMode {
    /// Break when (length - rest_length) / rest_length > strain_threshold.
    StrainBased,
    /// Break when constraint force magnitude > stress_threshold.
    StressBased,
    /// Break when either criterion is met.
    Combined,
}
/// A single fragment generated by fracture.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Node indices belonging to this fragment.
    pub node_indices: Vec<usize>,
    /// Centre of mass position (computed from node positions).
    pub center_of_mass: [f64; 3],
    /// Linear velocity of the fragment's centre of mass.
    pub velocity: [f64; 3],
    /// Mass of the fragment.
    pub mass: f64,
}
impl Fragment {
    /// Build a fragment from a list of node indices.
    ///
    /// The centre of mass, velocity, and mass are computed from the mesh data.
    pub fn from_nodes(
        node_indices: Vec<usize>,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
    ) -> Self {
        let total_mass: f64 = node_indices.iter().map(|&i| masses[i]).sum();
        let mut com = [0.0; 3];
        let mut vel = [0.0; 3];
        if total_mass > 1e-30 {
            for &i in &node_indices {
                let m = masses[i];
                com[0] += m * positions[i][0];
                com[1] += m * positions[i][1];
                com[2] += m * positions[i][2];
                vel[0] += m * velocities[i][0];
                vel[1] += m * velocities[i][1];
                vel[2] += m * velocities[i][2];
            }
            let inv_m = 1.0 / total_mass;
            com = [com[0] * inv_m, com[1] * inv_m, com[2] * inv_m];
            vel = [vel[0] * inv_m, vel[1] * inv_m, vel[2] * inv_m];
        }
        Self {
            node_indices,
            center_of_mass: com,
            velocity: vel,
            mass: total_mass,
        }
    }
    /// Kinetic energy of the fragment.
    pub fn kinetic_energy(&self) -> f64 {
        let v = self.velocity;
        0.5 * self.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    }
}
/// A progressive damage model where damage accumulates over load cycles.
pub struct DamageModel {
    /// Current damage level in \[0, 1\]. 0 = intact, 1 = fully broken.
    pub damage: f64,
    /// Damage increment per unit of excess strain beyond threshold.
    pub damage_rate: f64,
    /// Strain threshold below which no damage accumulates.
    pub threshold_strain: f64,
    /// Stiffness degradation factor: effective_stiffness = (1 - damage) * original.
    pub original_stiffness: f64,
}
impl DamageModel {
    /// Create a new damage model.
    pub fn new(damage_rate: f64, threshold_strain: f64, original_stiffness: f64) -> Self {
        Self {
            damage: 0.0,
            damage_rate,
            threshold_strain,
            original_stiffness,
        }
    }
    /// Accumulate damage based on current strain.
    ///
    /// Returns true if the element is fully broken (damage >= 1).
    pub fn accumulate(&mut self, strain: f64) -> bool {
        if strain > self.threshold_strain {
            let excess = strain - self.threshold_strain;
            self.damage = (self.damage + self.damage_rate * excess).min(1.0);
        }
        self.damage >= 1.0
    }
    /// Current effective stiffness accounting for damage.
    pub fn effective_stiffness(&self) -> f64 {
        (1.0 - self.damage) * self.original_stiffness
    }
    /// Whether the element is fully broken.
    pub fn is_broken(&self) -> bool {
        self.damage >= 1.0
    }
    /// Reset damage to zero.
    pub fn reset(&mut self) {
        self.damage = 0.0;
    }
}
/// Fracture threshold -- when constraint/edge strain exceeds this, it breaks.
pub struct FractureModel {
    /// Relative elongation at which edge breaks (e.g. 0.3 = 30 %).
    pub strain_threshold: f64,
    /// Force threshold for breaking.
    pub stress_threshold: f64,
    /// Which criterion to apply.
    pub mode: FractureMode,
}
/// A single edge/constraint that can fracture.
pub struct FractureEdge {
    /// Index of the first node.
    pub node_a: usize,
    /// Index of the second node.
    pub node_b: usize,
    /// Natural (rest) length of this edge.
    pub rest_length: f64,
    /// Spring stiffness coefficient.
    pub stiffness: f64,
    /// Whether this edge has been permanently broken.
    pub broken: bool,
    /// Accumulated damage in \[0, 1\].
    pub damage: f64,
}
impl FractureEdge {
    /// Create a new intact edge.
    pub fn new(node_a: usize, node_b: usize, rest_length: f64, stiffness: f64) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            stiffness,
            broken: false,
            damage: 0.0,
        }
    }
    /// Compute (|r_b - r_a| - rest_length) / rest_length.
    pub fn current_strain(&self, positions: &[[f64; 3]]) -> f64 {
        let pa = positions[self.node_a];
        let pb = positions[self.node_b];
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if self.rest_length.abs() < f64::EPSILON {
            return 0.0;
        }
        (len - self.rest_length) / self.rest_length
    }
    /// Effective stiffness accounting for damage.
    pub fn effective_stiffness(&self) -> f64 {
        (1.0 - self.damage) * self.stiffness
    }
    /// Hooke's law spring force along the edge.
    ///
    /// Returns `(force_on_a, force_on_b)`. Returns `([0,0,0], [0,0,0])` if
    /// the edge is broken.
    pub fn spring_force(&self, positions: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
        if self.broken {
            return ([0.0; 3], [0.0; 3]);
        }
        let pa = positions[self.node_a];
        let pb = positions[self.node_b];
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < f64::EPSILON {
            return ([0.0; 3], [0.0; 3]);
        }
        let extension = len - self.rest_length;
        let mag = self.effective_stiffness() * extension;
        let nx = dx / len;
        let ny = dy / len;
        let nz = dz / len;
        let fa = [mag * nx, mag * ny, mag * nz];
        let fb = [-mag * nx, -mag * ny, -mag * nz];
        (fa, fb)
    }
    /// Compute the elastic energy stored in this edge.
    ///
    /// `E = 0.5 * k * (L - L0)^2`
    pub fn elastic_energy(&self, positions: &[[f64; 3]]) -> f64 {
        if self.broken {
            return 0.0;
        }
        let pa = positions[self.node_a];
        let pb = positions[self.node_b];
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        let extension = len - self.rest_length;
        0.5 * self.effective_stiffness() * extension * extension
    }
    /// Check fracture criteria and mark edge broken if exceeded.
    ///
    /// Returns `true` on the first fracture event (transition intact -> broken).
    pub fn check_fracture(&mut self, positions: &[[f64; 3]], model: &FractureModel) -> bool {
        if self.broken {
            return false;
        }
        let strain = self.current_strain(positions);
        let (fa, _) = self.spring_force(positions);
        let force_mag = (fa[0] * fa[0] + fa[1] * fa[1] + fa[2] * fa[2]).sqrt();
        let strain_exceeded = strain > model.strain_threshold;
        let stress_exceeded = force_mag > model.stress_threshold;
        let should_break = match model.mode {
            FractureMode::StrainBased => strain_exceeded,
            FractureMode::StressBased => stress_exceeded,
            FractureMode::Combined => strain_exceeded || stress_exceeded,
        };
        if should_break {
            self.broken = true;
            true
        } else {
            false
        }
    }
    /// Accumulate damage on this edge.
    ///
    /// Returns true if damage reaches 1.0 and the edge breaks.
    pub fn accumulate_damage(
        &mut self,
        positions: &[[f64; 3]],
        damage_rate: f64,
        threshold: f64,
    ) -> bool {
        if self.broken {
            return false;
        }
        let strain = self.current_strain(positions);
        if strain > threshold {
            let excess = strain - threshold;
            self.damage = (self.damage + damage_rate * excess).min(1.0);
            if self.damage >= 1.0 {
                self.broken = true;
                return true;
            }
        }
        false
    }
}
