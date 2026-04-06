//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// A single tissue layer with distinct mechanical properties.
#[derive(Debug, Clone)]
pub struct TissueLayer {
    /// Layer name (e.g. "epidermis", "dermis").
    pub name: String,
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Layer thickness \[m\].
    pub thickness: f64,
    /// Density \[kg/m³\].
    pub density: f64,
    /// Fracture toughness K_Ic \[Pa·√m\].
    pub fracture_toughness: f64,
}
impl TissueLayer {
    /// Create a new tissue layer.
    pub fn new(
        name: impl Into<String>,
        young_modulus: f64,
        poisson_ratio: f64,
        thickness: f64,
        density: f64,
        fracture_toughness: f64,
    ) -> Self {
        Self {
            name: name.into(),
            young_modulus,
            poisson_ratio,
            thickness,
            density,
            fracture_toughness,
        }
    }
    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Bulk modulus K = E / (3(1−2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
    /// First Lamé parameter λ = ν E / ((1+ν)(1−2ν)).
    pub fn lame_lambda(&self) -> f64 {
        let nu = self.poisson_ratio;
        self.young_modulus * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }
}
/// Geometry type of a surgical tool tip.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToolGeometry {
    /// Spherical tip (radius stored in `tip_radius`).
    Sphere,
    /// Cylindrical needle (radius stored in `tip_radius`).
    Needle,
    /// Flat blade (width stored in `tip_radius`).
    Blade,
}
/// Haptic feedback state for a surgical tool: force, torque, energy.
#[derive(Debug, Clone)]
pub struct HapticFeedback {
    /// Resultant force at tool \[N\].
    pub force: [f64; 3],
    /// Resultant torque at tool \[N·m\].
    pub torque: [f64; 3],
    /// Contact deformation energy \[J\].
    pub contact_energy: f64,
    /// Cutting energy expended this step \[J\].
    pub cutting_energy: f64,
    /// Average deformation gradient determinant (volume change) of contacted elements.
    pub mean_det_f: f64,
}
impl HapticFeedback {
    /// Create a zeroed haptic feedback state.
    pub fn zero() -> Self {
        Self {
            force: [0.0; 3],
            torque: [0.0; 3],
            contact_energy: 0.0,
            cutting_energy: 0.0,
            mean_det_f: 1.0,
        }
    }
    /// Magnitude of feedback force \[N\].
    pub fn force_magnitude(&self) -> f64 {
        norm3(self.force)
    }
    /// Magnitude of feedback torque \[N·m\].
    pub fn torque_magnitude(&self) -> f64 {
        norm3(self.torque)
    }
    /// Add a contact contribution from a tissue node with penetration depth.
    ///
    /// F += k_contact · depth · n_contact
    pub fn add_contact(
        &mut self,
        normal: [f64; 3],
        depth: f64,
        stiffness: f64,
        moment_arm: [f64; 3],
    ) {
        let f = scale3(normal, stiffness * depth);
        self.force = add3(self.force, f);
        self.torque = add3(self.torque, cross3(moment_arm, f));
        self.contact_energy += 0.5 * stiffness * depth * depth;
    }
    /// Reset feedback to zero (call at the start of each step).
    pub fn reset(&mut self) {
        *self = Self::zero();
    }
}
/// A multi-layer soft tissue model with heterogeneous elastic properties and
/// pre-stress support.
#[derive(Debug, Clone)]
pub struct TissueModel {
    /// Tissue nodes (particle positions).
    pub nodes: Vec<TissueNode>,
    /// Tetrahedral element connectivity (node indices).
    pub elements: Vec<[usize; 4]>,
    /// Layer definitions.
    pub layers: Vec<TissueLayer>,
    /// Global pre-stress tensor component σ₀ \[Pa\] (isotropic, for simplicity).
    pub pre_stress: f64,
    /// Rayleigh damping coefficient α (mass-proportional).
    pub damping_alpha: f64,
    /// Rayleigh damping coefficient β (stiffness-proportional).
    pub damping_beta: f64,
}
impl TissueModel {
    /// Create a new tissue model.
    pub fn new(layers: Vec<TissueLayer>, pre_stress: f64) -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
            layers,
            pre_stress,
            damping_alpha: 0.02,
            damping_beta: 0.002,
        }
    }
    /// Add a node to the tissue.
    pub fn add_node(&mut self, node: TissueNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }
    /// Add a tetrahedral element.
    pub fn add_element(&mut self, indices: [usize; 4]) {
        self.elements.push(indices);
    }
    /// Compute the deformation gradient F for element `elem_idx` given the
    /// rest (reference) positions.
    ///
    /// F = (x_deformed − x0) · (X_rest − X0)⁻¹
    ///
    /// This returns F as a flat row-major 3×3 array.
    pub fn deformation_gradient(&self, elem_idx: usize, rest_positions: &[[f64; 3]]) -> [f64; 9] {
        let idx = self.elements[elem_idx];
        let x = [
            self.nodes[idx[0]].position,
            self.nodes[idx[1]].position,
            self.nodes[idx[2]].position,
            self.nodes[idx[3]].position,
        ];
        let x_ref = [
            rest_positions[idx[0]],
            rest_positions[idx[1]],
            rest_positions[idx[2]],
            rest_positions[idx[3]],
        ];
        let d0 = sub3(x_ref[1], x_ref[0]);
        let d1 = sub3(x_ref[2], x_ref[0]);
        let d2 = sub3(x_ref[3], x_ref[0]);
        let e0 = sub3(x[1], x[0]);
        let e1 = sub3(x[2], x[0]);
        let e2 = sub3(x[3], x[0]);
        let ds_inv = mat3_inv([
            d0[0], d1[0], d2[0], d0[1], d1[1], d2[1], d0[2], d1[2], d2[2],
        ]);
        let dm = [
            e0[0], e1[0], e2[0], e0[1], e1[1], e2[1], e0[2], e1[2], e2[2],
        ];
        mat3_mul(dm, ds_inv)
    }
    /// Linear elastic energy density W = λ/2 tr(ε)² + μ tr(ε²)
    /// where ε = (F + F^T)/2 − I (small strain approximation).
    pub fn elastic_energy_density(&self, elem_idx: usize, rest_positions: &[[f64; 3]]) -> f64 {
        let layer_idx = self.nodes[self.elements[elem_idx][0]].layer_index;
        let layer = &self.layers[layer_idx.min(self.layers.len() - 1)];
        let lam = layer.lame_lambda();
        let mu = layer.shear_modulus();
        let f = self.deformation_gradient(elem_idx, rest_positions);
        let eps = [
            f[0] - 1.0,
            0.5 * (f[1] + f[3]),
            0.5 * (f[2] + f[6]),
            0.5 * (f[3] + f[1]),
            f[4] - 1.0,
            0.5 * (f[5] + f[7]),
            0.5 * (f[6] + f[2]),
            0.5 * (f[7] + f[5]),
            f[8] - 1.0,
        ];
        let trace_eps = eps[0] + eps[4] + eps[8];
        let trace_eps2 = mat3_frobenius_sq(eps);
        0.5 * lam * trace_eps * trace_eps + mu * trace_eps2
    }
    /// Integrate one explicit time step for all free nodes.
    ///
    /// Uses forward Euler with Rayleigh damping.
    pub fn step(&mut self, dt: f64, external_force: [f64; 3]) {
        for node in &mut self.nodes {
            if node.is_static {
                continue;
            }
            let damp = scale3(node.velocity, -self.damping_alpha);
            let accel = scale3(
                add3(add3(node.force, damp), external_force),
                node.inv_mass(),
            );
            node.velocity = add3(node.velocity, scale3(accel, dt));
            node.position = add3(node.position, scale3(node.velocity, dt));
            node.force = [0.0; 3];
        }
    }
}
/// A rigid surgical tool with a pose, geometry, and contact detection.
#[derive(Debug, Clone)]
pub struct SurgicalTool {
    /// Tool identifier.
    pub name: String,
    /// Tool geometry type.
    pub geometry: ToolGeometry,
    /// Tip / blade radius or half-width \[m\].
    pub tip_radius: f64,
    /// Current tool-tip position \[m\].
    pub position: [f64; 3],
    /// Tool orientation as unit axis \[−\].
    pub orientation: [f64; 3],
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
}
impl SurgicalTool {
    /// Create a new surgical tool.
    pub fn new(
        name: impl Into<String>,
        geometry: ToolGeometry,
        tip_radius: f64,
        position: [f64; 3],
    ) -> Self {
        Self {
            name: name.into(),
            geometry,
            tip_radius,
            position,
            orientation: [0.0, 1.0, 0.0],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
        }
    }
    /// Move the tool tip by a translation vector.
    pub fn translate(&mut self, delta: [f64; 3]) {
        self.position = add3(self.position, delta);
        self.velocity = delta;
    }
    /// Detect contact with a tissue node.  Returns the penetration depth > 0
    /// if the node is inside the tool, or ≤ 0 if not.
    pub fn contact_depth(&self, node_pos: [f64; 3], node_radius: f64) -> f64 {
        let dist = norm3(sub3(node_pos, self.position));
        let threshold = self.tip_radius + node_radius;
        threshold - dist
    }
    /// Check if a tissue node is in contact with the tool.
    pub fn is_in_contact(&self, node_pos: [f64; 3], node_radius: f64) -> bool {
        self.contact_depth(node_pos, node_radius) > 0.0
    }
    /// Contact normal (from tool to node).
    pub fn contact_normal(&self, node_pos: [f64; 3]) -> [f64; 3] {
        normalize3(sub3(node_pos, self.position))
    }
    /// Advance the tool pose with a time step dt \[s\].
    pub fn integrate(&mut self, dt: f64) {
        self.position = add3(self.position, scale3(self.velocity, dt));
    }
}
/// Cutting model with crack-path tracking and energy criterion.
#[derive(Debug, Clone)]
pub struct CuttingModel {
    /// Tissue edges that can be cut.
    pub edges: Vec<TissueEdge>,
    /// Critical strain energy release rate G_c \[J/m²\].
    pub gc: f64,
    /// Blade sharpness factor s ∈ (0, 1]; 1 = perfectly sharp.
    pub sharpness: f64,
}
impl CuttingModel {
    /// Create a new cutting model.
    ///
    /// # Parameters
    /// - `gc`        : Critical energy release rate \[J/m²\].
    /// - `sharpness` : Blade sharpness ∈ (0, 1].
    pub fn new(gc: f64, sharpness: f64) -> Self {
        Self {
            edges: Vec::new(),
            gc,
            sharpness: sharpness.clamp(1e-6, 1.0),
        }
    }
    /// Register an edge as cuttable.
    pub fn add_edge(&mut self, edge: TissueEdge) {
        self.edges.push(edge);
    }
    /// Apply cutting force at a contact point, damaging nearby edges.
    ///
    /// Damage increment is proportional to applied force scaled by sharpness
    /// and divided by the critical energy release rate.
    ///
    /// # Parameters
    /// - `contact_pos`     : Tool contact position \[m\].
    /// - `cutting_force`   : Magnitude of cutting force \[N\].
    /// - `influence_radius`: Radius within which edges receive damage \[m\].
    /// - `nodes`           : Tissue node positions.
    pub fn apply_cut(
        &mut self,
        contact_pos: [f64; 3],
        cutting_force: f64,
        influence_radius: f64,
        nodes: &[[f64; 3]],
    ) {
        for edge in &mut self.edges {
            if edge.is_severed() {
                continue;
            }
            let mid = scale3(add3(nodes[edge.node_a], nodes[edge.node_b]), 0.5);
            let dist = norm3(sub3(mid, contact_pos));
            if dist < influence_radius {
                let attenuation = 1.0 - dist / influence_radius;
                let delta_d = self.sharpness * cutting_force * attenuation / self.gc;
                edge.apply_damage(delta_d);
            }
        }
    }
    /// Count severed edges.
    pub fn num_severed(&self) -> usize {
        self.edges.iter().filter(|e| e.is_severed()).count()
    }
    /// Crack path: ordered list of severed edge midpoints.
    ///
    /// Returns positions of midpoints of all severed edges given node positions.
    pub fn crack_path(&self, nodes: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.edges
            .iter()
            .filter(|e| e.is_severed())
            .map(|e| scale3(add3(nodes[e.node_a], nodes[e.node_b]), 0.5))
            .collect()
    }
}
/// A particle / node in the tissue mesh.
#[derive(Debug, Clone)]
pub struct TissueNode {
    /// World-space position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Accumulated force \[N\].
    pub force: [f64; 3],
    /// Node mass \[kg\].
    pub mass: f64,
    /// True if this node is pinned (zero inverse mass).
    pub is_static: bool,
    /// Layer index this node belongs to.
    pub layer_index: usize,
}
impl TissueNode {
    /// Create a dynamic tissue node.
    pub fn new(position: [f64; 3], mass: f64, layer_index: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            mass,
            is_static: false,
            layer_index,
        }
    }
    /// Create a static (pinned) tissue node.
    pub fn new_static(position: [f64; 3], layer_index: usize) -> Self {
        Self {
            is_static: true,
            mass: 0.0,
            ..Self::new(position, 0.0, layer_index)
        }
    }
    /// Inverse mass (0 for static nodes).
    pub fn inv_mass(&self) -> f64 {
        if self.is_static || self.mass <= 0.0 {
            0.0
        } else {
            1.0 / self.mass
        }
    }
}
/// Orchestrates one simulation step: tissue integration, tool contact, cutting.
#[derive(Debug, Clone)]
pub struct SimulationStep {
    /// Tissue model.
    pub tissue: TissueModel,
    /// Surgical tool.
    pub tool: SurgicalTool,
    /// Cutting model.
    pub cutting: CuttingModel,
    /// List of sutures.
    pub sutures: Vec<SuturModel>,
    /// Haptic feedback for the current step.
    pub haptic: HapticFeedback,
    /// Contact stiffness k_c \[N/m\].
    pub contact_stiffness: f64,
    /// Tissue node interaction radius \[m\].
    pub node_radius: f64,
    /// Rest node positions (used for deformation gradient computation).
    pub rest_positions: Vec<[f64; 3]>,
    /// Simulation time \[s\].
    pub time: f64,
}
impl SimulationStep {
    /// Create a new simulation step orchestrator.
    ///
    /// # Parameters
    /// - `tissue`           : Tissue model.
    /// - `tool`             : Surgical tool.
    /// - `cutting`          : Cutting model.
    /// - `contact_stiffness`: Contact stiffness \[N/m\].
    /// - `node_radius`      : Radius of each tissue node \[m\].
    pub fn new(
        tissue: TissueModel,
        tool: SurgicalTool,
        cutting: CuttingModel,
        contact_stiffness: f64,
        node_radius: f64,
    ) -> Self {
        let rest_positions = tissue.nodes.iter().map(|n| n.position).collect();
        Self {
            tissue,
            tool,
            cutting,
            sutures: Vec::new(),
            haptic: HapticFeedback::zero(),
            contact_stiffness,
            node_radius,
            rest_positions,
            time: 0.0,
        }
    }
    /// Add a suture to the simulation.
    pub fn add_suture(&mut self, suture: SuturModel) {
        self.sutures.push(suture);
    }
    /// Advance the simulation by time step dt \[s\].
    ///
    /// Steps performed:
    /// 1. Reset haptic feedback.
    /// 2. Detect and respond to tool-tissue contact.
    /// 3. Apply cutting forces.
    /// 4. Apply suture constraints.
    /// 5. Integrate tissue dynamics.
    /// 6. Advance healing.
    /// 7. Integrate tool pose.
    pub fn advance(&mut self, dt: f64, gravity: [f64; 3]) {
        self.haptic.reset();
        let tool_pos = self.tool.position;
        let tool_r = self.tool.tip_radius;
        let node_r = self.node_radius;
        let k_c = self.contact_stiffness;
        for node in &mut self.tissue.nodes {
            if node.is_static {
                continue;
            }
            let dist = norm3(sub3(node.position, tool_pos));
            let threshold = tool_r + node_r;
            if dist < threshold {
                let depth = threshold - dist;
                let normal = normalize3(sub3(node.position, tool_pos));
                let contact_force = scale3(normal, k_c * depth);
                node.force = add3(node.force, contact_force);
                let moment_arm = sub3(node.position, tool_pos);
                self.haptic
                    .add_contact(scale3(normal, -1.0), depth, k_c, moment_arm);
            }
        }
        let cutting_force = self.haptic.force_magnitude();
        let node_positions: Vec<[f64; 3]> = self.tissue.nodes.iter().map(|n| n.position).collect();
        self.cutting
            .apply_cut(tool_pos, cutting_force, tool_r * 2.0, &node_positions);
        self.haptic.cutting_energy += cutting_force * norm3(self.tool.velocity) * dt;
        for suture in &mut self.sutures {
            if !suture.active {
                continue;
            }
            let (force_a, _tension) = suture.tension_force(&self.tissue.nodes);
            let force_b = scale3(force_a, -1.0);
            let ia = suture.node_a;
            let ib = suture.node_b;
            self.tissue.nodes[ia].force = add3(self.tissue.nodes[ia].force, force_a);
            self.tissue.nodes[ib].force = add3(self.tissue.nodes[ib].force, force_b);
        }
        self.tissue.step(dt, gravity);
        for suture in &mut self.sutures {
            suture.advance_healing(dt);
        }
        self.tool.integrate(dt);
        self.time += dt;
    }
    /// Number of fully severed edges in the cutting model.
    pub fn num_cuts(&self) -> usize {
        self.cutting.num_severed()
    }
}
/// Neo-Hookean material constants for a soft tissue element.
#[derive(Debug, Clone)]
pub struct NeoHookeanParams {
    /// Shear modulus μ \[Pa\].
    pub mu: f64,
    /// Bulk modulus κ \[Pa\].
    pub kappa: f64,
    /// Viscous damping coefficient η \[Pa·s\].
    pub eta: f64,
}
impl NeoHookeanParams {
    /// Create neo-Hookean parameters from Young's modulus E and Poisson ratio ν,
    /// plus viscous damping coefficient η.
    pub fn from_young(young_modulus: f64, poisson_ratio: f64, eta: f64) -> Self {
        let nu = poisson_ratio;
        let mu = young_modulus / (2.0 * (1.0 + nu));
        let kappa = young_modulus / (3.0 * (1.0 - 2.0 * nu));
        Self { mu, kappa, eta }
    }
    /// Neo-Hookean strain energy density: W = μ/2 (I₁ - 3) - μ ln J + κ/2 (ln J)²
    ///
    /// where I₁ = tr(C) = tr(F^T F), J = det(F).
    pub fn strain_energy_density(&self, f: &[f64; 9]) -> f64 {
        let j = mat3_det(*f);
        if j <= 0.0 {
            return 0.0;
        }
        let ft_f = mat3_sym_product_fft(*f);
        let i1 = ft_f[0] + ft_f[4] + ft_f[8];
        let ln_j = j.ln();
        self.mu / 2.0 * (i1 - 3.0) - self.mu * ln_j + self.kappa / 2.0 * ln_j * ln_j
    }
    /// Viscous stress contribution: σ_visc = η * (F + F^T) / 2 (simplified).
    pub fn viscous_stress(&self, f: &[f64; 9]) -> f64 {
        let sym_trace = (f[0] + f[4] + f[8]) / 3.0;
        self.eta * sym_trace
    }
}
/// A neo-Hookean viscoelastic soft tissue simulation model.
///
/// Each tetrahedral element has its own neo-Hookean material properties and
/// a viscous damping term that resists rapid deformation.
#[derive(Debug, Clone)]
pub struct SoftTissueModel {
    /// Node positions (world space) \[m\].
    pub nodes: Vec<[f64; 3]>,
    /// Node velocities \[m/s\].
    pub velocities: Vec<[f64; 3]>,
    /// Node masses \[kg\].
    pub masses: Vec<f64>,
    /// Boolean static mask (true = fixed).
    pub fixed: Vec<bool>,
    /// Tetrahedral element connectivity.
    pub elements: Vec<[usize; 4]>,
    /// Neo-Hookean material per element.
    pub materials: Vec<NeoHookeanParams>,
    /// Rest (reference) positions \[m\].
    pub rest_positions: Vec<[f64; 3]>,
    /// Rayleigh mass damping α.
    pub alpha: f64,
}
impl SoftTissueModel {
    /// Create a new soft tissue model.
    pub fn new(alpha: f64) -> Self {
        Self {
            nodes: Vec::new(),
            velocities: Vec::new(),
            masses: Vec::new(),
            fixed: Vec::new(),
            elements: Vec::new(),
            materials: Vec::new(),
            rest_positions: Vec::new(),
            alpha,
        }
    }
    /// Add a node to the model.
    ///
    /// Returns the node index.
    pub fn add_node(&mut self, position: [f64; 3], mass: f64, is_fixed: bool) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(position);
        self.velocities.push([0.0; 3]);
        self.masses.push(mass);
        self.fixed.push(is_fixed);
        self.rest_positions.push(position);
        idx
    }
    /// Add a tetrahedral element with the given material.
    pub fn add_element(&mut self, indices: [usize; 4], material: NeoHookeanParams) {
        self.elements.push(indices);
        self.materials.push(material);
    }
    /// Compute the deformation gradient F for element `e`.
    pub fn deformation_gradient(&self, e: usize) -> [f64; 9] {
        let idx = self.elements[e];
        let x = [
            self.nodes[idx[0]],
            self.nodes[idx[1]],
            self.nodes[idx[2]],
            self.nodes[idx[3]],
        ];
        let x_ref = [
            self.rest_positions[idx[0]],
            self.rest_positions[idx[1]],
            self.rest_positions[idx[2]],
            self.rest_positions[idx[3]],
        ];
        let d0 = sub3(x_ref[1], x_ref[0]);
        let d1 = sub3(x_ref[2], x_ref[0]);
        let d2 = sub3(x_ref[3], x_ref[0]);
        let e0 = sub3(x[1], x[0]);
        let e1 = sub3(x[2], x[0]);
        let e2 = sub3(x[3], x[0]);
        let ds = [
            d0[0], d1[0], d2[0], d0[1], d1[1], d2[1], d0[2], d1[2], d2[2],
        ];
        let dm = [
            e0[0], e1[0], e2[0], e0[1], e1[1], e2[1], e0[2], e1[2], e2[2],
        ];
        mat3_mul(dm, mat3_inv(ds))
    }
    /// Strain energy density for element `e`.
    pub fn element_energy(&self, e: usize) -> f64 {
        let f = self.deformation_gradient(e);
        self.materials[e].strain_energy_density(&f)
    }
    /// Total strain energy across all elements.
    pub fn total_energy(&self) -> f64 {
        (0..self.elements.len())
            .map(|e| self.element_energy(e))
            .sum()
    }
    /// Volume of tetrahedral element `e` in deformed configuration.
    pub fn element_volume(&self, e: usize) -> f64 {
        let idx = self.elements[e];
        let p0 = self.nodes[idx[0]];
        let p1 = self.nodes[idx[1]];
        let p2 = self.nodes[idx[2]];
        let p3 = self.nodes[idx[3]];
        let a = sub3(p1, p0);
        let b = sub3(p2, p0);
        let c = sub3(p3, p0);
        let cr = cross3(b, c);
        (dot3(a, cr) / 6.0).abs()
    }
    /// Explicit Euler step integrating tissue dynamics.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        for i in 0..self.nodes.len() {
            if self.fixed[i] || self.masses[i] <= 0.0 {
                continue;
            }
            let damp = scale3(self.velocities[i], -self.alpha);
            let accel = scale3(add3(gravity, damp), 1.0);
            self.velocities[i] = add3(self.velocities[i], scale3(accel, dt));
            self.nodes[i] = add3(self.nodes[i], scale3(self.velocities[i], dt));
        }
    }
    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }
    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
}
/// Simulation of topological cuts using a damage-based approach.
///
/// Tracks a set of edges that can be progressively damaged until severed,
/// at which point they are removed from the active mesh connectivity.
#[derive(Debug, Clone)]
pub struct CuttingSimulation {
    /// Edges that participate in cutting.
    pub edges: Vec<MeshEdge>,
    /// Critical energy release rate \[J/m²\].
    pub gc: f64,
    /// Tool sharpness ∈ (0, 1].
    pub sharpness: f64,
    /// Severed edge indices (for remeshing).
    pub severed_indices: Vec<usize>,
}
impl CuttingSimulation {
    /// Create a new cutting simulation.
    pub fn new(gc: f64, sharpness: f64) -> Self {
        Self {
            edges: Vec::new(),
            gc,
            sharpness: sharpness.clamp(1e-6, 1.0),
            severed_indices: Vec::new(),
        }
    }
    /// Add a cuttable edge.
    pub fn add_edge(&mut self, edge: MeshEdge) {
        self.edges.push(edge);
    }
    /// Apply a cutting stroke at `contact_pos` with magnitude `force` \[N\].
    ///
    /// Edges whose midpoints lie within `radius` \[m\] receive damage.
    pub fn apply_cut(
        &mut self,
        contact_pos: [f64; 3],
        force: f64,
        radius: f64,
        nodes: &[[f64; 3]],
    ) {
        for (i, edge) in self.edges.iter_mut().enumerate() {
            if edge.severed {
                continue;
            }
            let mid = edge.midpoint(nodes);
            let dist = norm3(sub3(mid, contact_pos));
            if dist < radius {
                let atten = 1.0 - dist / radius;
                let delta = self.sharpness * force * atten / self.gc;
                edge.apply_damage(delta);
                if edge.severed
                    && !CuttingSimulation::already_recorded(&Self::severed_snapshot(&[i]), &[i])
                {
                }
            }
        }
        for (i, edge) in self.edges.iter().enumerate() {
            if edge.severed && !self.severed_indices.contains(&i) {
                self.severed_indices.push(i);
            }
        }
    }
    fn already_recorded(indices: &[usize], new: &[usize]) -> bool {
        new.iter().all(|n| indices.contains(n))
    }
    fn severed_snapshot(indices: &[usize]) -> Vec<usize> {
        indices.to_vec()
    }
    /// Number of severed edges.
    pub fn n_severed(&self) -> usize {
        self.severed_indices.len()
    }
    /// Crack path: midpoints of all severed edges.
    pub fn crack_path(&self, nodes: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.severed_indices
            .iter()
            .map(|&i| self.edges[i].midpoint(nodes))
            .collect()
    }
    /// Remaining intact edges count.
    pub fn n_intact(&self) -> usize {
        self.edges.iter().filter(|e| !e.severed).count()
    }
    /// Remove all severed edges and return them (simulating remeshing step).
    pub fn remesh(&mut self) -> Vec<MeshEdge> {
        let severed: Vec<MeshEdge> = self
            .severed_indices
            .iter()
            .map(|&i| self.edges[i].clone())
            .collect();
        self.edges.retain(|e| !e.severed);
        self.severed_indices.clear();
        severed
    }
}
/// Parameters for surgical tool-tissue contact with Coulomb friction.
#[derive(Debug, Clone)]
pub struct ContactParams {
    /// Normal contact stiffness \[N/m\].
    pub stiffness: f64,
    /// Normal contact damping \[N·s/m\].
    pub damping: f64,
    /// Coefficient of friction μ.
    pub friction: f64,
    /// Maximum static friction multiplier (>= 1).
    pub static_friction_multiplier: f64,
}
impl ContactParams {
    /// Create contact parameters.
    pub fn new(stiffness: f64, damping: f64, friction: f64) -> Self {
        Self {
            stiffness,
            damping,
            friction,
            static_friction_multiplier: 1.2_f64,
        }
    }
}
/// Surgical tool-tissue collision response with Coulomb friction.
#[derive(Debug, Clone)]
pub struct CollisionResponseSurgical {
    /// Contact parameters.
    pub params: ContactParams,
    /// Last computed contact normal force \[N\].
    pub normal_force: f64,
    /// Last computed friction force vector \[N\].
    pub friction_force: [f64; 3],
}
impl CollisionResponseSurgical {
    /// Create a new collision response handler.
    pub fn new(params: ContactParams) -> Self {
        Self {
            params,
            normal_force: 0.0,
            friction_force: [0.0; 3],
        }
    }
    /// Compute the contact response for a node penetrating the tool.
    ///
    /// Returns `(force_on_node, force_on_tool)` pair.
    pub fn compute_contact(
        &mut self,
        penetration_depth: f64,
        contact_normal: [f64; 3],
        relative_velocity: [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        if penetration_depth <= 0.0 {
            self.normal_force = 0.0;
            self.friction_force = [0.0; 3];
            return ([0.0; 3], [0.0; 3]);
        }
        let vn = dot3(relative_velocity, contact_normal);
        let fn_val =
            (self.params.stiffness * penetration_depth - self.params.damping * vn).max(0.0);
        self.normal_force = fn_val;
        let vt = sub3(relative_velocity, scale3(contact_normal, vn));
        let vt_norm = norm3(vt);
        let ft = if vt_norm > 1e-12 {
            let ft_max = self.params.friction * fn_val;
            let dir = normalize3(vt);
            scale3(dir, -ft_max.min(vt_norm * self.params.stiffness))
        } else {
            [0.0; 3]
        };
        self.friction_force = ft;
        let force_on_node = add3(scale3(contact_normal, fn_val), ft);
        let force_on_tool = scale3(force_on_node, -1.0);
        (force_on_node, force_on_tool)
    }
    /// Total contact force magnitude on the node.
    pub fn total_force_magnitude(&self) -> f64 {
        let total = add3(
            scale3([0.0, 0.0, 1.0], self.normal_force),
            self.friction_force,
        );
        norm3(total)
    }
}
/// A suture thread instance connecting two tissue nodes.
#[derive(Debug, Clone)]
pub struct SutureThread {
    /// Index of first tissue node.
    pub node_a: usize,
    /// Index of second tissue node.
    pub node_b: usize,
    /// Rest length of the thread \[m\] (target gap at closure).
    pub rest_length: f64,
    /// Material properties.
    pub material: SutureMaterial,
    /// Whether the thread is intact.
    pub intact: bool,
    /// Accumulated healing fraction ∈ \[0, 1\].
    pub healing: f64,
}
impl SutureThread {
    /// Create a new suture thread.
    pub fn new(node_a: usize, node_b: usize, rest_length: f64, material: SutureMaterial) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            material,
            intact: true,
            healing: 0.0,
        }
    }
    /// Compute the force vector applied to node_a (towards node_b).
    ///
    /// Returns (\[force on A\], tension magnitude).
    pub fn force(&self, nodes: &[[f64; 3]]) -> ([f64; 3], f64) {
        if !self.intact {
            return ([0.0; 3], 0.0);
        }
        let pa = nodes[self.node_a];
        let pb = nodes[self.node_b];
        let delta = sub3(pb, pa);
        let dist = norm3(delta);
        if dist < 1e-12 {
            return ([0.0; 3], 0.0);
        }
        let elongation = (dist - self.rest_length).max(0.0);
        let tension = self.material.elastic_force(elongation);
        let dir = normalize3(delta);
        (scale3(dir, tension), tension)
    }
    /// Advance healing by dt \[s\] with time constant tau \[s\].
    pub fn heal(&mut self, dt: f64, tau: f64) {
        if self.intact {
            self.healing = (self.healing + dt / tau).min(1.0);
        }
    }
    /// Cut the thread (set intact = false).
    pub fn cut(&mut self) {
        self.intact = false;
    }
    /// Check if the thread has fully healed.
    pub fn is_healed(&self) -> bool {
        self.healing >= 1.0
    }
}
/// A blood vessel that can rupture and bleed.
#[derive(Debug, Clone)]
pub struct BloodVessel {
    /// Vessel position (centre) in tissue \[m\].
    pub position: [f64; 3],
    /// Vessel radius \[m\].
    pub radius: f64,
    /// Wall tensile strength \[Pa\].
    pub wall_strength: f64,
    /// Mean arterial pressure \[Pa\].
    pub pressure: f64,
    /// Whether the vessel is ruptured.
    pub ruptured: bool,
}
impl BloodVessel {
    /// Create a new intact blood vessel.
    pub fn new(position: [f64; 3], radius: f64, wall_strength: f64, pressure: f64) -> Self {
        Self {
            position,
            radius,
            wall_strength,
            pressure,
            ruptured: false,
        }
    }
    /// Check if the vessel ruptures under applied stress \[Pa\].
    ///
    /// Uses a simple von-Mises criterion: stress > wall_strength → rupture.
    pub fn check_rupture(&mut self, applied_stress: f64) -> bool {
        if applied_stress > self.wall_strength {
            self.ruptured = true;
        }
        self.ruptured
    }
    /// Bleeding rate \[m³/s\] based on Poiseuille flow from the ruptured end.
    ///
    /// Returns 0.0 if not ruptured.
    pub fn bleeding_rate(&self, blood_viscosity: f64) -> f64 {
        if !self.ruptured || blood_viscosity < 1e-15 {
            return 0.0;
        }
        let l = self.radius.max(1e-6);
        std::f64::consts::PI * self.radius.powi(4) * self.pressure / (8.0 * blood_viscosity * l)
    }
    /// Volume of blood lost after `dt` \[s\].
    pub fn blood_loss(&self, dt: f64, blood_viscosity: f64) -> f64 {
        self.bleeding_rate(blood_viscosity) * dt
    }
}
/// A topological edge between two tissue nodes.
#[derive(Debug, Clone)]
pub struct TissueEdge {
    /// First node index.
    pub node_a: usize,
    /// Second node index.
    pub node_b: usize,
    /// Rest length \[m\].
    pub rest_length: f64,
    /// Cut status.
    pub status: EdgeStatus,
    /// Accumulated damage D ∈ \[0, 1\].
    pub damage: f64,
}
impl TissueEdge {
    /// Create a new intact edge.
    pub fn new(node_a: usize, node_b: usize, rest_length: f64) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            status: EdgeStatus::Intact,
            damage: 0.0,
        }
    }
    /// Apply damage and update edge status.
    ///
    /// `delta_damage` is the increment of damage applied in one step.
    pub fn apply_damage(&mut self, delta_damage: f64) {
        self.damage = (self.damage + delta_damage).min(1.0);
        self.status = if self.damage >= 1.0 {
            EdgeStatus::Severed
        } else if self.damage > 0.0 {
            EdgeStatus::Partial(self.damage)
        } else {
            EdgeStatus::Intact
        };
    }
    /// True if this edge has been fully severed.
    pub fn is_severed(&self) -> bool {
        matches!(self.status, EdgeStatus::Severed)
    }
}
/// A half-edge for topological cutting.
#[derive(Debug, Clone)]
pub struct MeshEdge {
    /// Node index A.
    pub node_a: usize,
    /// Node index B.
    pub node_b: usize,
    /// Whether this edge has been severed.
    pub severed: bool,
    /// Accumulated damage ∈ \[0, 1\].
    pub damage: f64,
    /// Critical damage threshold for severance.
    pub damage_threshold: f64,
}
impl MeshEdge {
    /// Create a new intact edge.
    pub fn new(node_a: usize, node_b: usize, damage_threshold: f64) -> Self {
        Self {
            node_a,
            node_b,
            severed: false,
            damage: 0.0,
            damage_threshold,
        }
    }
    /// Apply damage increment and check for severance.
    pub fn apply_damage(&mut self, delta: f64) {
        self.damage = (self.damage + delta).min(1.0);
        if self.damage >= self.damage_threshold {
            self.severed = true;
        }
    }
    /// Midpoint of the edge given node positions.
    pub fn midpoint(&self, nodes: &[[f64; 3]]) -> [f64; 3] {
        scale3(add3(nodes[self.node_a], nodes[self.node_b]), 0.5)
    }
}
/// Model tracking all vessels in a tissue region.
#[derive(Debug, Clone)]
pub struct BleedingModel {
    /// Blood vessels in the tissue.
    pub vessels: Vec<BloodVessel>,
    /// Blood dynamic viscosity \[Pa·s\] (water ≈ 1e-3, blood ≈ 3e-3).
    pub blood_viscosity: f64,
    /// Total blood volume lost \[m³\].
    pub total_blood_loss: f64,
}
impl BleedingModel {
    /// Create a new bleeding model.
    pub fn new(vessels: Vec<BloodVessel>, blood_viscosity: f64) -> Self {
        Self {
            vessels,
            blood_viscosity,
            total_blood_loss: 0.0,
        }
    }
    /// Apply a cutting force at `contact_pos` and check vessel ruptures.
    ///
    /// Vessels within `radius` \[m\] experience stress proportional to the force.
    pub fn apply_trauma(&mut self, contact_pos: [f64; 3], force: f64, radius: f64) -> usize {
        let mut new_ruptures = 0;
        for vessel in &mut self.vessels {
            if vessel.ruptured {
                continue;
            }
            let dist = norm3(sub3(vessel.position, contact_pos));
            if dist < radius {
                let stress = force / (std::f64::consts::PI * (radius * radius).max(1e-20))
                    * (1.0 - dist / radius);
                if vessel.check_rupture(stress) {
                    new_ruptures += 1;
                }
            }
        }
        new_ruptures
    }
    /// Advance the bleeding model by `dt` \[s\].
    pub fn step(&mut self, dt: f64) {
        for vessel in &self.vessels {
            self.total_blood_loss += vessel.blood_loss(dt, self.blood_viscosity);
        }
    }
    /// Number of ruptured vessels.
    pub fn n_ruptured(&self) -> usize {
        self.vessels.iter().filter(|v| v.ruptured).count()
    }
    /// Total bleeding rate \[m³/s\] across all ruptured vessels.
    pub fn total_bleeding_rate(&self) -> f64 {
        self.vessels
            .iter()
            .map(|v| v.bleeding_rate(self.blood_viscosity))
            .sum()
    }
}
/// Healing state of a suture.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HealingState {
    /// Suture is active, no healing.
    Active,
    /// Tissue is partially healed; healing ∈ (0, 1).
    Healing(f64),
    /// Tissue fully healed; suture can be removed.
    Healed,
}
/// Mechanical properties of a suture thread.
#[derive(Debug, Clone)]
pub struct SutureMaterial {
    /// Thread Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Thread cross-sectional area \[m²\].
    pub cross_section: f64,
    /// Maximum tensile stress before failure \[Pa\].
    pub ultimate_stress: f64,
    /// Knot efficiency factor κ ∈ (0, 1]; 1 = no knot weakening.
    pub knot_efficiency: f64,
    /// Pretension force applied during tying \[N\].
    pub pretension: f64,
}
impl SutureMaterial {
    /// Create a new suture material.
    pub fn new(
        young_modulus: f64,
        cross_section: f64,
        ultimate_stress: f64,
        knot_efficiency: f64,
        pretension: f64,
    ) -> Self {
        Self {
            young_modulus,
            cross_section,
            ultimate_stress,
            knot_efficiency: knot_efficiency.clamp(0.0, 1.0),
            pretension,
        }
    }
    /// Axial stiffness EA \[N\] of the thread.
    pub fn axial_stiffness(&self) -> f64 {
        self.young_modulus * self.cross_section
    }
    /// Maximum allowable force considering knot efficiency.
    pub fn max_force(&self) -> f64 {
        self.ultimate_stress * self.cross_section * self.knot_efficiency
    }
    /// Elastic force for a given elongation \[m\] above rest length.
    pub fn elastic_force(&self, elongation: f64) -> f64 {
        (self.pretension + self.axial_stiffness() * elongation.max(0.0)).min(self.max_force())
    }
    /// Check if the suture would fail under the given applied force.
    pub fn will_fail(&self, applied_force: f64) -> bool {
        applied_force > self.max_force()
    }
}
/// Status of a tissue edge with respect to cutting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeStatus {
    /// Edge is intact.
    Intact,
    /// Edge is partially cut (crack is progressing).
    Partial(f64),
    /// Edge is fully severed.
    Severed,
}
/// A suture constraint linking two tissue nodes.
#[derive(Debug, Clone)]
pub struct SuturModel {
    /// First tissue node index.
    pub node_a: usize,
    /// Second tissue node index.
    pub node_b: usize,
    /// Target separation distance \[m\] (usually zero for wound closure).
    pub target_gap: f64,
    /// Suture stiffness \[N/m\].
    pub stiffness: f64,
    /// Maximum tension before suture fails \[N\].
    pub max_tension: f64,
    /// Healing state.
    pub healing: HealingState,
    /// Accumulated healing time \[s\].
    pub heal_time: f64,
    /// Time constant for healing \[s\].
    pub heal_tau: f64,
    /// Whether the suture is still active.
    pub active: bool,
}
impl SuturModel {
    /// Create a new suture between two nodes.
    ///
    /// # Parameters
    /// - `node_a`, `node_b` : Tissue node indices being sutured.
    /// - `target_gap`       : Desired gap between nodes after closure \[m\].
    /// - `stiffness`        : Spring stiffness \[N/m\].
    /// - `max_tension`      : Failure tension \[N\].
    /// - `heal_tau`         : Healing time constant \[s\].
    pub fn new(
        node_a: usize,
        node_b: usize,
        target_gap: f64,
        stiffness: f64,
        max_tension: f64,
        heal_tau: f64,
    ) -> Self {
        Self {
            node_a,
            node_b,
            target_gap,
            stiffness,
            max_tension,
            healing: HealingState::Active,
            heal_time: 0.0,
            heal_tau,
            active: true,
        }
    }
    /// Compute the suture force vector applied to node_a.
    ///
    /// F = k (d − target_gap) · n̂
    ///
    /// where d is the current distance and n̂ is the unit vector from a to b.
    pub fn tension_force(&self, nodes: &[TissueNode]) -> ([f64; 3], f64) {
        let pa = nodes[self.node_a].position;
        let pb = nodes[self.node_b].position;
        let delta = sub3(pb, pa);
        let dist = norm3(delta);
        if dist < 1e-12 {
            return ([0.0; 3], 0.0);
        }
        let n = normalize3(delta);
        let tension = self.stiffness * (dist - self.target_gap);
        let force = scale3(n, tension);
        (force, tension.abs())
    }
    /// Advance healing by time dt \[s\].
    ///
    /// Healing fraction increases exponentially: h(t) = 1 − e^(−t/τ).
    pub fn advance_healing(&mut self, dt: f64) {
        if !self.active {
            return;
        }
        self.heal_time += dt;
        let fraction = 1.0 - (-self.heal_time / self.heal_tau).exp();
        self.healing = if fraction >= 0.99 {
            HealingState::Healed
        } else {
            HealingState::Healing(fraction)
        };
    }
    /// Check suture tension against failure criterion.
    ///
    /// Returns `true` if the suture has failed (broken).
    pub fn check_failure(&mut self, nodes: &[TissueNode]) -> bool {
        let (_force, tension) = self.tension_force(nodes);
        if tension > self.max_tension {
            self.active = false;
            true
        } else {
            false
        }
    }
    /// Healing fraction ∈ \[0, 1\].
    pub fn healing_fraction(&self) -> f64 {
        match self.healing {
            HealingState::Active => 0.0,
            HealingState::Healing(f) => f,
            HealingState::Healed => 1.0,
        }
    }
}
