//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// Dynamic fracture toughness model (Freund's model).
///
/// The dynamic fracture toughness decreases with crack velocity:
/// `K_Id = K_Ic * (1 - v/cR)^m`
///
/// where `v` is crack velocity, `cR` is Rayleigh wave speed, and `m` is a
/// material constant (typically 0.5 for glass).
#[derive(Debug, Clone)]
pub struct DynamicFractureToughness {
    /// Static fracture toughness K_Ic (Pa·√m).
    pub k_ic: f64,
    /// Rayleigh wave speed cR (m/s).
    pub rayleigh_speed: f64,
    /// Velocity exponent m (material constant).
    pub exponent: f64,
}
impl DynamicFractureToughness {
    /// Create a new dynamic fracture toughness model.
    pub fn new(k_ic: f64, rayleigh_speed: f64, exponent: f64) -> Self {
        Self {
            k_ic,
            rayleigh_speed,
            exponent,
        }
    }
    /// Compute K_Id at the given crack velocity.
    ///
    /// `K_Id = K_Ic * (1 - v/cR)^m`  for v < cR
    pub fn k_id(&self, crack_velocity: f64) -> f64 {
        if self.rayleigh_speed < 1e-15 {
            return self.k_ic;
        }
        let v_norm = (crack_velocity / self.rayleigh_speed).clamp(0.0, 1.0 - 1e-9);
        self.k_ic * (1.0 - v_norm).powf(self.exponent)
    }
    /// Maximum crack velocity below the Rayleigh limit (Freund's terminal velocity).
    ///
    /// In practice, crack velocity asymptotically approaches cR but never reaches it.
    pub fn terminal_velocity_fraction(&self) -> f64 {
        0.95
    }
    /// Check if the given stress intensity K_I will drive crack growth.
    pub fn will_propagate(&self, k_i: f64, crack_velocity: f64) -> bool {
        k_i >= self.k_id(crack_velocity)
    }
}
/// A bond (spring) connecting two nodes in the fracture mesh.
pub struct FractureBond {
    /// Index of the first endpoint node.
    pub node_a: usize,
    /// Index of the second endpoint node.
    pub node_b: usize,
    /// Rest length (m).
    pub rest_length: f64,
    /// Spring stiffness (N/m).
    pub stiffness: f64,
    /// Maximum force magnitude before the bond breaks (N).
    pub strength: f64,
    /// Energy released per unit area when the bond fractures (J/m²).
    pub fracture_energy: f64,
    /// Whether this bond has been broken.
    pub broken: bool,
    /// Normalised stretch: |current_length − rest_length| / rest_length.
    pub current_stretch: f64,
}
impl FractureBond {
    /// Create a new intact bond with the given parameters.
    pub fn new(
        node_a: usize,
        node_b: usize,
        rest_length: f64,
        stiffness: f64,
        strength: f64,
        fracture_energy: f64,
    ) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            stiffness,
            strength,
            fracture_energy,
            broken: false,
            current_stretch: 0.0,
        }
    }
}
/// A bond-based soft-body mesh that supports dynamic fracture propagation.
pub struct FractureMesh {
    /// All nodes in the mesh.
    pub nodes: Vec<FractureNode>,
    /// All bonds in the mesh.
    pub bonds: Vec<FractureBond>,
    /// Material properties.
    pub material: FractureMaterial,
    /// Active crack tips.
    pub crack_tips: Vec<CrackTip>,
    /// Gravitational acceleration vector (m/s²).
    pub gravity: [f64; 3],
}
impl FractureMesh {
    /// Create a 2-D (XY) grid of nodes with `nx × ny` vertices separated by
    /// spacing `dx`.  Bonds are added for horizontal, vertical, and diagonal
    /// neighbours.  Bond stiffness is derived from the material's Young's
    /// modulus.
    pub fn new_grid(nx: usize, ny: usize, dx: f64, mat: FractureMaterial) -> Self {
        assert!(nx >= 2 && ny >= 2, "Grid must be at least 2×2");
        let node_mass = mat.density * dx * dx;
        let stiffness = mat.youngs_modulus * dx;
        let strength = mat.tensile_strength * dx * dx;
        let g_c = mat.critical_energy_release_rate();
        let mut nodes = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let pos = [i as f64 * dx, j as f64 * dx, 0.0];
                nodes.push(FractureNode::new(pos, node_mass));
            }
        }
        let idx = |i: usize, j: usize| j * nx + i;
        let mut bonds = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                let a = idx(i, j);
                if i + 1 < nx {
                    let b = idx(i + 1, j);
                    bonds.push(FractureBond::new(a, b, dx, stiffness, strength, g_c));
                }
                if j + 1 < ny {
                    let b = idx(i, j + 1);
                    bonds.push(FractureBond::new(a, b, dx, stiffness, strength, g_c));
                }
                if i + 1 < nx && j + 1 < ny {
                    let b = idx(i + 1, j + 1);
                    let diag = dx * std::f64::consts::SQRT_2;
                    bonds.push(FractureBond::new(a, b, diag, stiffness, strength, g_c));
                }
                if i > 0 && j + 1 < ny {
                    let b = idx(i - 1, j + 1);
                    let diag = dx * std::f64::consts::SQRT_2;
                    bonds.push(FractureBond::new(a, b, diag, stiffness, strength, g_c));
                }
            }
        }
        Self {
            nodes,
            bonds,
            material: mat,
            crack_tips: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
        }
    }
    /// Compute the elastic bond forces on every node.  Returns a `Vec` with
    /// one `[f64; 3]` per node containing the accumulated force (N).
    pub fn compute_bond_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0f64; 3]; self.nodes.len()];
        for bond in &self.bonds {
            if bond.broken {
                continue;
            }
            let pa = self.nodes[bond.node_a].pos;
            let pb = self.nodes[bond.node_b].pos;
            let delta = vec_sub(pb, pa);
            let current_len = vec_length(delta);
            if current_len < 1e-15 {
                continue;
            }
            let extension = current_len - bond.rest_length;
            let force_mag = bond.stiffness * extension;
            let dir = vec_normalize(delta);
            let f = vec_scale(dir, force_mag);
            forces[bond.node_a] = vec_add(forces[bond.node_a], f);
            forces[bond.node_b] = vec_sub(forces[bond.node_b], f);
        }
        forces
    }
    /// Check every bond for fracture.  A bond is broken if its normalised
    /// stretch exceeds the critical stretch *or* the force magnitude exceeds
    /// the bond strength.  Broken bonds increment the damage of their
    /// endpoint nodes and may spawn or update crack tips.
    pub fn check_fracture(&mut self) {
        let eps_c = self.material.critical_stretch();
        let positions: Vec<[f64; 3]> = self.nodes.iter().map(|n| n.pos).collect();
        let mut newly_broken: Vec<(usize, usize, [f64; 3])> = Vec::new();
        for bond in &mut self.bonds {
            if bond.broken {
                continue;
            }
            let pa = positions[bond.node_a];
            let pb = positions[bond.node_b];
            let delta = vec_sub(pb, pa);
            let current_len = vec_length(delta);
            let stretch = if bond.rest_length > 1e-15 {
                (current_len - bond.rest_length).abs() / bond.rest_length
            } else {
                0.0
            };
            bond.current_stretch = stretch;
            let extension = current_len - bond.rest_length;
            let force_mag = (bond.stiffness * extension).abs();
            if stretch > eps_c || force_mag > bond.strength {
                bond.broken = true;
                let mid = vec_scale(vec_add(pa, pb), 0.5);
                let bond_dir = vec_normalize(delta);
                let crack_dir = [-bond_dir[1], bond_dir[0], 0.0];
                newly_broken.push((bond.node_a, bond.node_b, crack_dir));
                let _ = mid;
            }
        }
        for (na, nb, crack_dir) in newly_broken {
            let total_bonds_a = self
                .bonds
                .iter()
                .filter(|b| b.node_a == na || b.node_b == na)
                .count()
                .max(1);
            let total_bonds_b = self
                .bonds
                .iter()
                .filter(|b| b.node_a == nb || b.node_b == nb)
                .count()
                .max(1);
            self.nodes[na].damage = (self.nodes[na].damage + 1.0 / total_bonds_a as f64).min(1.0);
            self.nodes[nb].damage = (self.nodes[nb].damage + 1.0 / total_bonds_b as f64).min(1.0);
            let tip_node = if self.nodes[na].damage >= self.nodes[nb].damage {
                na
            } else {
                nb
            };
            let already_tracked = self.crack_tips.iter().any(|ct| ct.node_idx == tip_node);
            if !already_tracked {
                self.crack_tips.push(CrackTip::new(tip_node, crack_dir));
            }
        }
    }
    /// Propagate each crack tip by breaking nearby bonds whose stress
    /// intensity exceeds K_Ic.
    ///
    /// The stress intensity factor K_I is estimated from the local elastic
    /// energy density of bonds connected to the crack-tip node.
    pub fn propagate_crack(&mut self, dt: f64) {
        let k_ic = self.material.fracture_toughness;
        let forces = self.compute_bond_forces();
        let positions: Vec<[f64; 3]> = self.nodes.iter().map(|n| n.pos).collect();
        let n_tips = self.crack_tips.len();
        for t in 0..n_tips {
            let tip_node = self.crack_tips[t].node_idx;
            let f_tip = forces[tip_node];
            let f_mag = vec_length(f_tip);
            let connected_rest_lengths: Vec<f64> = self
                .bonds
                .iter()
                .filter(|b| !b.broken && (b.node_a == tip_node || b.node_b == tip_node))
                .map(|b| b.rest_length)
                .collect();
            if connected_rest_lengths.is_empty() {
                self.crack_tips[t].stress_intensity = 0.0;
                continue;
            }
            let r_avg: f64 =
                connected_rest_lengths.iter().sum::<f64>() / connected_rest_lengths.len() as f64;
            let k_i = f_mag / (std::f64::consts::PI * r_avg).sqrt();
            self.crack_tips[t].stress_intensity = k_i;
            if k_i <= k_ic {
                continue;
            }
            let crack_dir = self.crack_tips[t].direction;
            let mut best_bond: Option<usize> = None;
            let mut best_cos: f64 = -1.0;
            for (bi, bond) in self.bonds.iter().enumerate() {
                if bond.broken {
                    continue;
                }
                let (other, _self_node) = if bond.node_a == tip_node {
                    (bond.node_b, bond.node_a)
                } else if bond.node_b == tip_node {
                    (bond.node_a, bond.node_b)
                } else {
                    continue;
                };
                let d = vec_normalize(vec_sub(positions[other], positions[tip_node]));
                let cos_angle = vec_dot(d, crack_dir);
                if cos_angle > best_cos {
                    best_cos = cos_angle;
                    best_bond = Some(bi);
                }
            }
            if let Some(bi) = best_bond {
                let (na, nb) = (self.bonds[bi].node_a, self.bonds[bi].node_b);
                self.bonds[bi].broken = true;
                let tb_a = self
                    .bonds
                    .iter()
                    .filter(|b| b.node_a == na || b.node_b == na)
                    .count()
                    .max(1);
                let tb_b = self
                    .bonds
                    .iter()
                    .filter(|b| b.node_a == nb || b.node_b == nb)
                    .count()
                    .max(1);
                self.nodes[na].damage = (self.nodes[na].damage + 1.0 / tb_a as f64).min(1.0);
                self.nodes[nb].damage = (self.nodes[nb].damage + 1.0 / tb_b as f64).min(1.0);
                let next_node = if self.bonds[bi].node_a == tip_node {
                    nb
                } else {
                    na
                };
                let advance_dist = vec_length(vec_sub(positions[next_node], positions[tip_node]));
                self.crack_tips[t].total_length += advance_dist;
                let speed = advance_dist / dt.max(1e-15);
                self.crack_tips[t].velocity = vec_scale(crack_dir, speed);
                self.crack_tips[t].node_idx = next_node;
            }
        }
    }
    /// Semi-implicit Euler integration step.
    ///
    /// 1. Accumulate bond forces and gravity.
    /// 2. Update velocities.
    /// 3. Update positions.
    /// 4. Run fracture check.
    pub fn integrate(&mut self, dt: f64) {
        let bond_forces = self.compute_bond_forces();
        let g = self.gravity;
        for (i, node) in self.nodes.iter_mut().enumerate() {
            if node.fixed {
                continue;
            }
            let fg = vec_scale(g, node.mass);
            let f_total = vec_add(fg, bond_forces[i]);
            let accel = vec_scale(f_total, 1.0 / node.mass);
            node.vel = vec_add(node.vel, vec_scale(accel, dt));
            node.pos = vec_add(node.pos, vec_scale(node.vel, dt));
        }
        self.check_fracture();
    }
    /// Advance the simulation by one timestep `dt` (integration + crack propagation).
    pub fn step(&mut self, dt: f64) {
        self.integrate(dt);
        self.propagate_crack(dt);
    }
    /// Number of broken bonds.
    pub fn broken_bond_count(&self) -> usize {
        self.bonds.iter().filter(|b| b.broken).count()
    }
    /// Total fracture energy released (sum of `fracture_energy` for broken bonds).
    pub fn total_fracture_energy(&self) -> f64 {
        self.bonds
            .iter()
            .filter(|b| b.broken)
            .map(|b| b.fracture_energy)
            .sum()
    }
    /// Damage value for every node.
    pub fn damage_field(&self) -> Vec<f64> {
        self.nodes.iter().map(|n| n.damage).collect()
    }
    /// Sum of `total_length` across all crack tips.
    pub fn crack_length(&self) -> f64 {
        self.crack_tips.iter().map(|ct| ct.total_length).sum()
    }
}
impl FractureMesh {
    /// Apply a velocity damping factor to all non-fixed nodes.
    ///
    /// `damping` should be in `[0, 1]`: 1 = no damping, 0 = freeze all motion.
    pub fn apply_damping(&mut self, damping: f64) {
        for node in &mut self.nodes {
            if !node.fixed {
                node.vel = vec_scale(node.vel, damping);
            }
        }
    }
    /// Perform multiple integration sub-steps of size `dt / sub_steps`.
    pub fn step_substeps(&mut self, dt: f64, sub_steps: usize) {
        let sub_dt = dt / sub_steps.max(1) as f64;
        for _ in 0..sub_steps {
            self.step(sub_dt);
        }
    }
    /// Apply a concentrated external force to node `idx`.
    pub fn apply_nodal_force(&mut self, idx: usize, force: [f64; 3], dt: f64) {
        if idx < self.nodes.len() && !self.nodes[idx].fixed {
            let accel = vec_scale(force, 1.0 / self.nodes[idx].mass);
            self.nodes[idx].vel = vec_add(self.nodes[idx].vel, vec_scale(accel, dt));
        }
    }
    /// Attempt to branch every crack tip whose speed exceeds the Yoffe threshold.
    ///
    /// New crack tips are spawned at `±branch_angle_deg` from the current
    /// direction.  The original tip is retained.
    pub fn try_branch_cracks(&mut self, branch_angle_deg: f64, threshold_fraction: f64) {
        let c_r = rayleigh_wave_speed(&self.material);
        let mut new_tips: Vec<CrackTip> = Vec::new();
        for tip in &self.crack_tips {
            let speed = vec_length(tip.velocity);
            if should_branch(speed, c_r, threshold_fraction) {
                let (d_pos, d_neg) = branch_directions(tip.direction, branch_angle_deg);
                new_tips.push(CrackTip::new(tip.node_idx, d_pos));
                new_tips.push(CrackTip::new(tip.node_idx, d_neg));
            }
        }
        for nt in new_tips {
            let already = self.crack_tips.iter().any(|ct| {
                ct.node_idx == nt.node_idx
                    && (vec_dot(ct.direction, nt.direction) - 1.0).abs() < 1e-6
            });
            if !already {
                self.crack_tips.push(nt);
            }
        }
    }
    /// Update dynamic stress intensity factors for all crack tips using the
    /// Freund velocity correction.
    pub fn update_dynamic_sif(&mut self) {
        let c_r = rayleigh_wave_speed(&self.material);
        for tip in &mut self.crack_tips {
            let speed = vec_length(tip.velocity);
            let k_dyn = dynamic_stress_intensity(tip.stress_intensity, speed, c_r);
            tip.stress_intensity = k_dyn;
        }
    }
    /// Return the Rayleigh wave speed for this mesh's material.
    pub fn rayleigh_speed(&self) -> f64 {
        rayleigh_wave_speed(&self.material)
    }
    /// Return the longitudinal wave speed for this mesh's material.
    pub fn longitudinal_speed(&self) -> f64 {
        longitudinal_wave_speed(&self.material)
    }
    /// Compute total elastic strain energy: Σ ½ k (extension)².
    pub fn total_strain_energy(&self) -> f64 {
        self.bonds
            .iter()
            .filter(|b| !b.broken)
            .map(|b| {
                let pa = self.nodes[b.node_a].pos;
                let pb = self.nodes[b.node_b].pos;
                let ext = vec_length(vec_sub(pb, pa)) - b.rest_length;
                0.5 * b.stiffness * ext * ext
            })
            .sum()
    }
    /// Count intact bonds.
    pub fn intact_bond_count(&self) -> usize {
        self.bonds.iter().filter(|b| !b.broken).count()
    }
    /// Maximum stretch across all intact bonds.
    pub fn max_stretch(&self) -> f64 {
        self.bonds
            .iter()
            .filter(|b| !b.broken)
            .map(|b| b.current_stretch)
            .fold(0.0_f64, f64::max)
    }
    /// Average damage across all nodes.
    pub fn average_damage(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        self.nodes.iter().map(|n| n.damage).sum::<f64>() / self.nodes.len() as f64
    }
    /// Maximum speed of any non-fixed node (m/s).
    pub fn max_node_speed(&self) -> f64 {
        self.nodes
            .iter()
            .filter(|n| !n.fixed)
            .map(|n| vec_length(n.vel))
            .fold(0.0_f64, f64::max)
    }
    /// Kinetic energy of the mesh: Σ ½ m |v|².
    pub fn kinetic_energy(&self) -> f64 {
        self.nodes
            .iter()
            .filter(|n| !n.fixed)
            .map(|n| 0.5 * n.mass * vec_dot(n.vel, n.vel))
            .sum()
    }
    /// Centre of mass of the mesh.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.nodes.iter().map(|n| n.mass).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for n in &self.nodes {
            com[0] += n.mass * n.pos[0];
            com[1] += n.mass * n.pos[1];
            com[2] += n.mass * n.pos[2];
        }
        vec_scale(com, 1.0 / total_mass)
    }
    /// Apply a uniform tensile load in the Y direction: top-row nodes are
    /// displaced by `displacement` and bottom-row nodes are fixed.
    ///
    /// The mesh must have been built with `new_grid`.
    pub fn apply_tensile_load_y(&mut self, displacement: f64) {
        let y_min = self
            .nodes
            .iter()
            .map(|n| n.pos[1])
            .fold(f64::INFINITY, f64::min);
        let y_max = self
            .nodes
            .iter()
            .map(|n| n.pos[1])
            .fold(f64::NEG_INFINITY, f64::max);
        let tol = 1e-10;
        for n in &mut self.nodes {
            if (n.pos[1] - y_min).abs() < tol {
                n.fixed = true;
            } else if (n.pos[1] - y_max).abs() < tol {
                n.fixed = true;
                n.pos[1] += displacement;
            }
        }
    }
}
/// Crack branching model based on the maximum tangential stress criterion
/// (Erdogan & Sih 1963).
///
/// For a mixed-mode crack (mode I + mode II), the branching angle θ satisfies:
/// ```text
/// K_I sin(θ) + K_II (3 cos(θ) - 1) = 0
/// ```
/// Solved analytically: `tan(θ/2) = K_I / (4*K_II) - sign(K_II) * sqrt((K_I/(4*K_II))² + 0.5)`
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct CrackBranching {
    /// Mode-I stress intensity factor (Pa·√m).
    pub k_i: f64,
    /// Mode-II stress intensity factor (Pa·√m).
    pub k_ii: f64,
}
#[allow(dead_code)]
impl CrackBranching {
    /// Create a new branching model.
    pub fn new(k_i: f64, k_ii: f64) -> Self {
        Self { k_i, k_ii }
    }
    /// Compute the branching angle θ (radians) from the current crack plane.
    ///
    /// For pure mode-I (K_II = 0) → θ = 0 (straight ahead).
    /// Returns angle in \[-π/2, π/2\].
    pub fn compute_branching_angle(&self) -> f64 {
        if self.k_ii.abs() < 1e-30 {
            return 0.0;
        }
        let ki = self.k_i;
        let kii = self.k_ii;
        let disc = (ki * ki + 8.0 * kii * kii).sqrt();
        let theta = 2.0 * ((ki - disc) / (4.0 * kii)).atan();
        theta.clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2)
    }
    /// Compute the branching angle in degrees.
    pub fn compute_branching_angle_deg(&self) -> f64 {
        self.compute_branching_angle().to_degrees()
    }
    /// Check whether mixed-mode conditions favour branching.
    ///
    /// Branching is favoured when |K_II / K_I| > threshold (typically ~ 0.2).
    pub fn is_branching_favoured(&self, threshold_ratio: f64) -> bool {
        if self.k_i.abs() < 1e-30 {
            return false;
        }
        (self.k_ii / self.k_i).abs() >= threshold_ratio
    }
}
/// Represents a propagating crack tip in the mesh.
pub struct CrackTip {
    /// Index of the node at the current crack-tip location.
    pub node_idx: usize,
    /// Unit vector pointing in the direction of crack propagation.
    pub direction: [f64; 3],
    /// Crack propagation velocity vector (m/s).
    pub velocity: [f64; 3],
    /// Mode-I stress intensity factor K_I (Pa·√m).
    pub stress_intensity: f64,
    /// Total crack length accumulated so far (m).
    pub total_length: f64,
}
impl CrackTip {
    /// Create a new crack tip.
    pub fn new(node_idx: usize, direction: [f64; 3]) -> Self {
        Self {
            node_idx,
            direction,
            velocity: [0.0; 3],
            stress_intensity: 0.0,
            total_length: 0.0,
        }
    }
}
/// A single node in the fracture mesh.
pub struct FractureNode {
    /// Current position (m).
    pub pos: [f64; 3],
    /// Current velocity (m/s).
    pub vel: [f64; 3],
    /// Nodal mass (kg).
    pub mass: f64,
    /// If `true` the node is kinematically fixed (Dirichlet BC).
    pub fixed: bool,
    /// Local damage indicator in \[0, 1\].  0 = intact, 1 = fully damaged.
    pub damage: f64,
}
impl FractureNode {
    /// Create a new node at `pos` with the given mass.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            fixed: false,
            damage: 0.0,
        }
    }
}
/// Linear-elastic material properties relevant to brittle fracture.
pub struct FractureMaterial {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν (dimensionless).
    pub poisson_ratio: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
    /// Tensile strength σ_t (Pa).
    pub tensile_strength: f64,
    /// Mode-I fracture toughness K_Ic (Pa·√m).
    pub fracture_toughness: f64,
}
impl FractureMaterial {
    /// Create a new material.
    pub fn new(
        youngs_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        tensile_strength: f64,
        fracture_toughness: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
            density,
            tensile_strength,
            fracture_toughness,
        }
    }
    /// Critical stretch at which a bond breaks: ε_c = σ_t / E.
    pub fn critical_stretch(&self) -> f64 {
        self.tensile_strength / self.youngs_modulus
    }
    /// Critical energy release rate (plane stress): G_c = K_Ic² / E.
    pub fn critical_energy_release_rate(&self) -> f64 {
        self.fracture_toughness * self.fracture_toughness / self.youngs_modulus
    }
}
/// Dynamic fracture energy balance tracker.
///
/// Tracks kinetic energy, strain energy, surface energy created by fracture,
/// and external work to verify the energy balance during dynamic fracture.
#[derive(Debug, Clone)]
pub struct DynamicEnergyBalance {
    /// Kinetic energy at current time.
    pub kinetic_energy: f64,
    /// Elastic strain energy at current time.
    pub strain_energy: f64,
    /// Total fracture surface energy (sum of Gc * area of all broken surfaces).
    pub fracture_surface_energy: f64,
    /// Cumulative external work done on the system.
    pub external_work: f64,
    /// History of total energy (kinetic + strain + fracture) per time step.
    pub energy_history: Vec<f64>,
}
impl DynamicEnergyBalance {
    /// Create a new energy balance tracker.
    pub fn new() -> Self {
        Self {
            kinetic_energy: 0.0,
            strain_energy: 0.0,
            fracture_surface_energy: 0.0,
            external_work: 0.0,
            energy_history: Vec::new(),
        }
    }
    /// Update kinetic energy from nodal velocities and masses.
    pub fn update_kinetic(&mut self, nodes: &[FractureNode]) {
        self.kinetic_energy = nodes
            .iter()
            .map(|n| {
                let v2 = n.vel[0] * n.vel[0] + n.vel[1] * n.vel[1] + n.vel[2] * n.vel[2];
                0.5 * n.mass * v2
            })
            .sum();
    }
    /// Update strain energy from bonds.
    ///
    /// Strain energy per bond: U = 0.5 * k * (stretch * rest_length)^2.
    pub fn update_strain(&mut self, bonds: &[FractureBond]) {
        self.strain_energy = bonds
            .iter()
            .filter(|b| !b.broken)
            .map(|b| {
                let delta = b.current_stretch * b.rest_length;
                0.5 * b.stiffness * delta * delta
            })
            .sum();
    }
    /// Add fracture surface energy when a bond breaks.
    ///
    /// `bond_area` is the representative area associated with the bond (m^2).
    pub fn add_fracture_event(&mut self, gc: f64, bond_area: f64) {
        self.fracture_surface_energy += gc * bond_area;
    }
    /// Record current total energy to history.
    pub fn record(&mut self) {
        let total = self.kinetic_energy + self.strain_energy + self.fracture_surface_energy;
        self.energy_history.push(total);
    }
    /// Return the total mechanical + fracture energy.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy + self.strain_energy + self.fracture_surface_energy
    }
    /// Compute relative energy drift from initial value.
    pub fn relative_drift(&self) -> f64 {
        if self.energy_history.len() < 2 {
            return 0.0;
        }
        let e0 = self.energy_history[0];
        if e0.abs() < 1e-60 {
            return 0.0;
        }
        let e_last = *self
            .energy_history
            .last()
            .expect("collection should not be empty");
        (e_last - e0).abs() / e0.abs()
    }
}
/// Cohesive zone traction-separation law (linear softening).
///
/// Models the process zone ahead of a crack tip as a cohesive interface
/// that transmits traction up to a maximum value, then softens linearly
/// to zero at the critical opening displacement.
#[derive(Debug, Clone)]
pub struct CohesiveZone {
    /// Maximum traction (cohesive strength) T_max (Pa).
    pub t_max: f64,
    /// Critical separation δ_c (m) at which traction = 0.
    pub delta_c: f64,
    /// Current opening displacement δ (m).
    pub delta: f64,
    /// Maximum opening displacement reached (for irreversibility).
    pub delta_max: f64,
    /// Whether the interface has fully failed (δ ≥ δ_c).
    pub failed: bool,
}
impl CohesiveZone {
    /// Create a new cohesive zone with given strength and critical separation.
    pub fn new(t_max: f64, delta_c: f64) -> Self {
        Self {
            t_max,
            delta_c,
            delta: 0.0,
            delta_max: 0.0,
            failed: false,
        }
    }
    /// Compute the traction for the current opening displacement.
    ///
    /// Linear softening:
    /// - `T = T_max * δ / δ_c`  if loading from virgin state
    /// - Linear softening after peak: `T = T_max * (1 - δ/δ_c)`
    /// - `T = 0` if failed
    pub fn traction(&self, delta: f64) -> f64 {
        if self.failed || delta <= 0.0 {
            return 0.0;
        }
        if delta >= self.delta_c {
            return 0.0;
        }
        self.t_max * (1.0 - delta / self.delta_c)
    }
    /// Update the state given a new opening displacement.
    ///
    /// Returns the traction at the new opening.
    pub fn update(&mut self, new_delta: f64) -> f64 {
        if self.failed {
            return 0.0;
        }
        if new_delta >= self.delta_c {
            self.failed = true;
            self.delta = self.delta_c;
            return 0.0;
        }
        if new_delta > self.delta_max {
            self.delta_max = new_delta;
        }
        self.delta = new_delta.max(0.0);
        self.traction(self.delta)
    }
    /// Fracture energy: area under the traction-separation curve.
    ///
    /// For linear softening: `G_c = 0.5 * T_max * δ_c`.
    pub fn fracture_energy(&self) -> f64 {
        0.5 * self.t_max * self.delta_c
    }
    /// Damage variable D in \[0, 1\] (0 = intact, 1 = fully failed).
    pub fn damage(&self) -> f64 {
        if self.failed {
            1.0
        } else {
            (self.delta_max / self.delta_c).clamp(0.0, 1.0)
        }
    }
}
/// Tracks crack opening displacement along a crack path.
///
/// Stores (position along crack, opening displacement) pairs.
#[derive(Debug, Clone)]
pub struct CodProfile {
    /// Normalized positions along the crack (0 = tip, 1 = mouth).
    pub s: Vec<f64>,
    /// Opening displacement at each position (m).
    pub opening: Vec<f64>,
}
impl CodProfile {
    /// Create an empty COD profile.
    pub fn new() -> Self {
        Self {
            s: Vec::new(),
            opening: Vec::new(),
        }
    }
    /// Add a measurement point.
    pub fn add_point(&mut self, s: f64, opening: f64) {
        self.s.push(s);
        self.opening.push(opening);
    }
    /// Maximum opening displacement (CTOD at crack mouth).
    pub fn max_opening(&self) -> f64 {
        self.opening.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Crack tip opening displacement (minimum opening at tip end, s≈0).
    pub fn ctod(&self) -> f64 {
        if self.s.is_empty() {
            return 0.0;
        }
        self.s
            .iter()
            .zip(self.opening.iter())
            .min_by(|(s1, _), (s2, _)| s1.partial_cmp(s2).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(_, &o)| o)
            .unwrap_or(0.0)
    }
    /// Interpolate opening at position `s` using linear interpolation.
    pub fn interpolate(&self, s: f64) -> f64 {
        if self.s.len() < 2 {
            return self.opening.first().copied().unwrap_or(0.0);
        }
        for i in 0..self.s.len() - 1 {
            let s0 = self.s[i];
            let s1 = self.s[i + 1];
            if s >= s0 && s <= s1 {
                let t = (s - s0) / (s1 - s0).max(1e-15);
                return self.opening[i] + t * (self.opening[i + 1] - self.opening[i]);
            }
        }
        *self.opening.last().expect("collection should not be empty")
    }
}
/// Crack propagation model using the dynamic fracture toughness K_ID.
///
/// Provides the KID-based crack-tip velocity formula from Freund (1990):
/// ```text
/// v_tip = c_R * ( 1 - K_Id(v) / K_I )   when K_I > K_Ic, else 0
/// ```
/// Solved iteratively to self-consistency.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CrackPropagation {
    /// Dynamic fracture toughness model.
    pub toughness: DynamicFractureToughness,
    /// Applied (static) stress intensity factor K_I (Pa·√m).
    pub k_i_applied: f64,
}
#[allow(dead_code)]
impl CrackPropagation {
    /// Create a new propagation model.
    pub fn new(k_ic: f64, rayleigh_speed: f64, exponent: f64, k_i_applied: f64) -> Self {
        Self {
            toughness: DynamicFractureToughness::new(k_ic, rayleigh_speed, exponent),
            k_i_applied,
        }
    }
    /// Compute the self-consistent crack-tip velocity (m/s).
    ///
    /// Iterates the fixed-point equation:
    /// `v = c_R * max(0, 1 - K_Id(v) / K_I)`
    /// until convergence or `max_iter` iterations.
    pub fn compute_crack_tip_velocity(&self, max_iter: usize) -> f64 {
        let c_r = self.toughness.rayleigh_speed;
        let k_i = self.k_i_applied;
        let k_ic = self.toughness.k_ic;
        if k_i <= k_ic || c_r < 1e-15 {
            return 0.0;
        }
        let mut v = c_r * 0.1;
        for _ in 0..max_iter {
            let k_id = self.toughness.k_id(v);
            let v_new = c_r * (1.0 - k_id / k_i).max(0.0);
            if (v_new - v).abs() < 1e-6 * c_r {
                return v_new;
            }
            v = v_new;
        }
        v
    }
    /// Compute the crack-tip velocity directly without iteration (single-step estimate).
    ///
    /// Uses the static K_Id at v = 0 to get an upper-bound estimate.
    pub fn crack_tip_velocity_static_estimate(&self) -> f64 {
        let k_id_static = self.toughness.k_id(0.0);
        let k_i = self.k_i_applied;
        let c_r = self.toughness.rayleigh_speed;
        c_r * (1.0 - k_id_static / k_i).max(0.0)
    }
}
/// Fragment size distribution model based on the Mott-Grady theory for
/// dynamic fragmentation of a ring or rod under rapid expansion.
///
/// The characteristic fragment length (mean fragment size) is:
/// ```text
/// L_f = (24 * K_Ic² / (ρ * ε̇² * E))^(1/3)
/// ```
/// where:
/// * K_Ic – fracture toughness (Pa·√m)
/// * ρ    – material density (kg/m³)
/// * ε̇   – strain rate (1/s)
/// * E    – Young's modulus (Pa)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FragmentDistribution {
    /// Fracture toughness K_Ic (Pa·√m).
    pub fracture_toughness: f64,
    /// Material density ρ (kg/m³).
    pub density: f64,
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
}
#[allow(dead_code)]
impl FragmentDistribution {
    /// Create a new fragment distribution model.
    pub fn new(fracture_toughness: f64, density: f64, youngs_modulus: f64) -> Self {
        Self {
            fracture_toughness,
            density,
            youngs_modulus,
        }
    }
    /// Compute the Mott-Grady mean fragment size (m) at the given strain rate (1/s).
    ///
    /// Returns 0 if the strain rate is zero or negative.
    pub fn compute_mott_grady(&self, strain_rate: f64) -> f64 {
        if strain_rate < 1e-30 || self.density < 1e-30 || self.youngs_modulus < 1e-30 {
            return 0.0;
        }
        let kic = self.fracture_toughness;
        let rho = self.density;
        let edot = strain_rate;
        let e = self.youngs_modulus;
        let numerator = 24.0 * kic * kic;
        let denominator = rho * edot * edot * e;
        (numerator / denominator).cbrt()
    }
    /// Estimate the number of fragments per unit volume.
    ///
    /// N_f ≈ 1 / L_f³
    pub fn fragment_number_density(&self, strain_rate: f64) -> f64 {
        let lf = self.compute_mott_grady(strain_rate);
        if lf < 1e-30 {
            return 0.0;
        }
        1.0 / (lf * lf * lf)
    }
    /// Compute the fragment kinetic energy per unit volume assuming fragments
    /// fly apart at their characteristic velocity v = ε̇ * L_f / 2.
    pub fn fragment_kinetic_energy_density(&self, strain_rate: f64) -> f64 {
        let lf = self.compute_mott_grady(strain_rate);
        let v_frag = strain_rate * lf / 2.0;
        0.5 * self.density * v_frag * v_frag
    }
}
