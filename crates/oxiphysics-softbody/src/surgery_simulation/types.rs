//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::elasticity::{
    TetRest, explicit_stability_limit, linear_elastic_energy_density, linear_elastic_stress,
    mat3_adjugate, mat3_symmetric_part, mat3_transpose, stability_substeps,
};
use super::functions::*;
use std::collections::HashSet;
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
///
/// The tetrahedral elements are linear elastic (small strain). The Lamé
/// parameters of an element come from the layer of its first node, the
/// deformation is measured against [`TissueModel::rest_positions`], and the
/// isotropic [`TissueModel::pre_stress`] is added to the element stress.
/// Energy is dissipated by Rayleigh damping `C = α M + β K` with the
/// mass-proportional coefficient [`TissueModel::damping_alpha`] and the
/// stiffness-proportional coefficient [`TissueModel::damping_beta`].
///
/// Small-strain elasticity is geometrically linear: infinitesimal rigid
/// rotations of an element are force-free, and a large rotation is not resisted
/// as it should be. The model is therefore meant for tissue that is anchored
/// well enough (enough pinned nodes) for the deformation to stay small; use
/// [`SoftTissueModel`], whose neo-Hookean energy is rotation invariant, for
/// large deformations and rotations.
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
    /// Rayleigh damping coefficient α \[1/s\] (mass-proportional: force `-α m v`).
    pub damping_alpha: f64,
    /// Rayleigh damping coefficient β \[s\] (stiffness-proportional: force `-β K v`).
    pub damping_beta: f64,
    /// Rest (reference) node positions \[m\], recorded by [`TissueModel::add_node`].
    ///
    /// The element forces measure the deformation against these positions.
    pub rest_positions: Vec<[f64; 3]>,
    /// Severed flag of each element (parallel to `elements`).
    ///
    /// A severed element no longer carries load, see
    /// [`TissueModel::sever_element`].
    pub severed_elements: Vec<bool>,
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
            rest_positions: Vec::new(),
            severed_elements: Vec::new(),
        }
    }
    /// Add a node to the tissue.
    ///
    /// The position the node has when it is added becomes its rest position.
    pub fn add_node(&mut self, node: TissueNode) -> usize {
        let idx = self.nodes.len();
        self.rest_positions.push(node.position);
        self.nodes.push(node);
        idx
    }
    /// Add a tetrahedral element.
    pub fn add_element(&mut self, indices: [usize; 4]) {
        self.elements.push(indices);
        self.severed_elements.push(false);
    }
    /// `true` if the element has been severed and therefore carries no load.
    pub fn is_element_severed(&self, element: usize) -> bool {
        self.severed_elements.get(element).copied().unwrap_or(false)
    }
    /// Sever an element: it stops exerting elastic forces.
    ///
    /// Returns `false` (and changes nothing) if `element` is out of range.
    pub fn sever_element(&mut self, element: usize) -> bool {
        if element >= self.elements.len() {
            return false;
        }
        if self.severed_elements.len() < self.elements.len() {
            self.severed_elements.resize(self.elements.len(), false);
        }
        self.severed_elements[element] = true;
        true
    }
    /// Sever every element that contains both end nodes of one of the given
    /// (cut) edges.
    ///
    /// Returns the number of elements that were newly severed.
    pub fn sever_elements_by_edges(
        &mut self,
        edges: impl IntoIterator<Item = (usize, usize)>,
    ) -> usize {
        let cut: HashSet<(usize, usize)> = edges
            .into_iter()
            .map(|(a, b)| (a.min(b), a.max(b)))
            .collect();
        if cut.is_empty() {
            return 0;
        }
        let mut newly_severed = 0;
        for element in 0..self.elements.len() {
            if self.is_element_severed(element) {
                continue;
            }
            let nodes = self.elements[element];
            let touches_cut_edge = (0..4).any(|a| {
                ((a + 1)..4).any(|b| {
                    let (lo, hi) = (nodes[a].min(nodes[b]), nodes[a].max(nodes[b]));
                    cut.contains(&(lo, hi))
                })
            });
            if touches_cut_edge && self.sever_element(element) {
                newly_severed += 1;
            }
        }
        newly_severed
    }
    /// Layer that governs element `element`: the layer of its first node
    /// (clamped to the last layer).  `None` if the model has no layers or the
    /// element refers to a missing node.
    fn layer_of_element(&self, element: usize) -> Option<&TissueLayer> {
        let first_node = self.elements.get(element)?[0];
        let layer_index = self.nodes.get(first_node)?.layer_index;
        self.layers
            .get(layer_index.min(self.layers.len().saturating_sub(1)))
    }
    /// Elements that currently carry load, with their precomputed rest data.
    ///
    /// Severed elements, elements that refer to missing nodes, flat rest
    /// tetrahedra and elements whose layer has non-finite Lamé parameters
    /// (for example Poisson's ratio 0.5) are left out.
    fn load_bearing_elements(&self) -> Vec<TissueElement> {
        let known_nodes = self.nodes.len().min(self.rest_positions.len());
        let mut elements = Vec::with_capacity(self.elements.len());
        for (element, &nodes) in self.elements.iter().enumerate() {
            if self.is_element_severed(element) || nodes.iter().any(|&i| i >= known_nodes) {
                continue;
            }
            let Some(layer) = self.layer_of_element(element) else {
                continue;
            };
            let lame_lambda = layer.lame_lambda();
            let shear_modulus = layer.shear_modulus();
            if !lame_lambda.is_finite() || !shear_modulus.is_finite() {
                continue;
            }
            let Some(rest) = TetRest::from_positions(nodes.map(|i| self.rest_positions[i])) else {
                continue;
            };
            elements.push(TissueElement {
                nodes,
                rest,
                lame_lambda,
                shear_modulus,
            });
        }
        elements
    }
    /// Add the internal forces of `elements` to `forces` (indexed by node).
    ///
    /// The element force is `-V0 σ Dm⁻ᵀ` for the stress
    /// `σ = C : sym(∇(u + β v)) + σ₀ I`, where `u` is the displacement from the
    /// rest position and `v` the velocity: the `u` part is the elastic force
    /// `-K u`, the `β v` part the stiffness-proportional Rayleigh damping
    /// `-β K v`, and `σ₀ I` the isotropic pre-stress.
    fn accumulate_internal_forces(&self, elements: &[TissueElement], forces: &mut [[f64; 3]]) {
        for element in elements {
            let field = element.nodes.map(|i| {
                let node = &self.nodes[i];
                add3(
                    sub3(node.position, self.rest_positions[i]),
                    scale3(node.velocity, self.damping_beta),
                )
            });
            let gradient = element.rest.gradient(&field);
            let stress = linear_elastic_stress(
                &gradient,
                element.lame_lambda,
                element.shear_modulus,
                self.pre_stress,
            );
            let nodal = element.rest.nodal_forces(&stress);
            for (&node, force) in element.nodes.iter().zip(nodal.iter()) {
                forces[node] = add3(forces[node], *force);
            }
        }
    }
    /// Internal force on every node \[N\] for the current positions and
    /// velocities: linear-elastic element forces, pre-stress and
    /// stiffness-proportional damping.
    ///
    /// Severed elements and elements without a usable rest shape or material
    /// do not contribute. Applied forces ([`TissueNode::force`]) and gravity
    /// are not included.
    pub fn internal_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0; 3]; self.nodes.len()];
        let elements = self.load_bearing_elements();
        self.accumulate_internal_forces(&elements, &mut forces);
        forces
    }
    /// Strain energy stored in the load-bearing elements \[J\]:
    /// `Σ V0 (λ/2 tr(ε)² + μ tr(ε²) + σ₀ tr(ε))` with the small strain
    /// `ε = sym(∇u)` of the displacement from the rest positions.
    ///
    /// Without damping the internal forces are exactly the negative gradient of
    /// this energy with respect to the node positions.
    pub fn elastic_energy(&self) -> f64 {
        self.load_bearing_elements()
            .iter()
            .map(|element| {
                let displacement = element
                    .nodes
                    .map(|i| sub3(self.nodes[i].position, self.rest_positions[i]));
                let gradient = element.rest.gradient(&displacement);
                let strain = mat3_symmetric_part(gradient);
                let density = linear_elastic_energy_density(
                    &strain,
                    element.lame_lambda,
                    element.shear_modulus,
                );
                let pre_stress_work = self.pre_stress * (strain[0] + strain[4] + strain[8]);
                element.rest.volume * (density + pre_stress_work)
            })
            .sum()
    }
    /// Largest time step \[s\] for which explicit integration of this tissue
    /// is guaranteed to be stable.
    ///
    /// It follows from a conservative bound on the highest natural frequency of
    /// the linear-elastic model and the Rayleigh damping coefficients, see
    /// [`TissueModel::step`]. Infinite if nothing limits the step (no elastic
    /// element and no damping).
    pub fn critical_time_step(&self) -> f64 {
        self.stable_time_step(&self.load_bearing_elements())
    }
    fn stable_time_step(&self, elements: &[TissueElement]) -> f64 {
        let mut stiffness = vec![0.0; self.nodes.len()];
        for element in elements {
            let sums = element
                .rest
                .stiffness_row_sums(element.lame_lambda, element.shear_modulus);
            for (&node, sum) in element.nodes.iter().zip(sums) {
                stiffness[node] += sum;
            }
        }
        let omega_sq = self
            .nodes
            .iter()
            .zip(stiffness.iter())
            .filter(|(node, _)| !node.is_static && node.inv_mass() > 0.0)
            .map(|(node, &k)| k * node.inv_mass())
            .fold(0.0_f64, f64::max);
        let damping_rate = self.damping_alpha.max(0.0) + self.damping_beta.max(0.0) * omega_sq;
        explicit_stability_limit(omega_sq.sqrt(), damping_rate)
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
        let Some(layer) = self.layer_of_element(elem_idx) else {
            return 0.0;
        };
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
    /// Integrate the tissue forward by `dt` \[s\] with semi-implicit
    /// (symplectic) Euler and Rayleigh damping.
    ///
    /// `gravity` is an *acceleration* \[m/s²\] such as `[0.0, -9.81, 0.0]`; it
    /// acts on every free node independent of its mass. Forces accumulated in
    /// [`TissueNode::force`] \[N\] are applied for the whole step and cleared
    /// afterwards. Static nodes and nodes without positive mass do not move.
    ///
    /// The acceleration of a free node of mass `m` is
    ///
    /// `a = (f_applied + f_internal) / m + g - α v`,
    ///
    /// where `f_internal` collects the linear-elastic element forces of the
    /// load-bearing elements (see [`TissueModel::internal_forces`]) and
    /// `-α v` is the mass-proportional Rayleigh damping.
    ///
    /// Explicit integration is only conditionally stable. The step is
    /// therefore split into the smallest number of equal sub-steps that stays
    /// below half of [`TissueModel::critical_time_step`], capped at 4096
    /// sub-steps; beyond the cap the caller has to reduce `dt`, the stiffness
    /// or the stiffness-proportional damping `damping_beta`. For fine meshes
    /// the explicitly integrated damping term `-β K v` limits the step to about
    /// `2 / (β ω²)`, which is usually the limit that hits the cap first.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let elements = self.load_bearing_elements();
        let substeps = stability_substeps(dt, self.stable_time_step(&elements));
        let h = dt / substeps as f64;
        let alpha = self.damping_alpha;
        let mut internal = vec![[0.0; 3]; self.nodes.len()];
        for _ in 0..substeps {
            internal.fill([0.0; 3]);
            self.accumulate_internal_forces(&elements, &mut internal);
            for (node, internal_force) in self.nodes.iter_mut().zip(internal.iter()) {
                let inv_mass = node.inv_mass();
                if node.is_static || inv_mass == 0.0 {
                    continue;
                }
                let force = add3(node.force, *internal_force);
                let accel = sub3(
                    add3(scale3(force, inv_mass), gravity),
                    scale3(node.velocity, alpha),
                );
                node.velocity = add3(node.velocity, scale3(accel, h));
                node.position = add3(node.position, scale3(node.velocity, h));
            }
        }
        for node in &mut self.nodes {
            node.force = [0.0; 3];
        }
    }
}
/// Precomputed data of a load-bearing element of a [`TissueModel`].
struct TissueElement {
    nodes: [usize; 4],
    rest: TetRest,
    lame_lambda: f64,
    shear_modulus: f64,
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
    /// 3. Apply cutting forces and sever the tissue elements that contain a
    ///    fully cut edge (see [`TissueModel::sever_elements_by_edges`]).
    /// 4. Apply suture constraints.
    /// 5. Integrate tissue dynamics; `gravity` is an acceleration \[m/s²\].
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
        let cut_edges: Vec<(usize, usize)> = self
            .cutting
            .edges
            .iter()
            .filter(|edge| edge.is_severed())
            .map(|edge| (edge.node_a, edge.node_b))
            .collect();
        if !cut_edges.is_empty() {
            self.tissue.sever_elements_by_edges(cut_edges);
        }
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
    /// Mean normal component `η tr(f) / 3` of the Newtonian viscous stress
    /// `η sym(f)` for the velocity gradient `f` \[1/s\].
    ///
    /// This is a scalar summary only; [`SoftTissueModel::internal_forces`] uses
    /// the full viscous stress tensor of every element.
    pub fn viscous_stress(&self, f: &[f64; 9]) -> f64 {
        let sym_trace = (f[0] + f[4] + f[8]) / 3.0;
        self.eta * sym_trace
    }
}
/// A neo-Hookean viscoelastic soft tissue simulation model.
///
/// Each tetrahedral element has its own [`NeoHookeanParams`]. The elastic
/// nodal forces are the negative gradient of the neo-Hookean strain energy
/// `V0 W(F)` (see [`NeoHookeanParams::strain_energy_density`]) of the element
/// with rest volume `V0`; elements that are inverted (`det F <= 0`) exert no
/// force. Dissipation comes from two terms: a Newtonian viscous stress
/// `σ_v = η sym(L)` driven by the velocity gradient `L = Ḟ F⁻¹` of the
/// element (the viscosity `η` of its [`NeoHookeanParams`], rotation invariant
/// and never energy-producing) and the mass-proportional Rayleigh damping `α`.
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
    /// Strain energy density \[Pa\] of element `e` (energy per unit rest
    /// volume).
    pub fn element_energy(&self, e: usize) -> f64 {
        let f = self.deformation_gradient(e);
        self.materials[e].strain_energy_density(&f)
    }
    /// Sum of the strain energy *densities* \[Pa\] of all elements, see
    /// [`SoftTissueModel::element_energy`].
    ///
    /// This is not an energy in joules; weigh each element with its rest
    /// volume, [`SoftTissueModel::strain_energy`], for that.
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
    /// Number of nodes that have every per-node array (position, velocity,
    /// mass, fixed flag and rest position); elements may only refer to those.
    fn known_node_count(&self) -> usize {
        self.nodes
            .len()
            .min(self.velocities.len())
            .min(self.masses.len())
            .min(self.fixed.len())
            .min(self.rest_positions.len())
    }
    /// Elements that can exert forces, with their precomputed rest data.
    ///
    /// Elements without a material, with missing nodes, a flat rest shape or
    /// non-finite material constants are left out.
    fn load_bearing_elements(&self) -> Vec<SoftTissueElement> {
        let known_nodes = self.known_node_count();
        let mut elements = Vec::with_capacity(self.elements.len());
        for (&nodes, material) in self.elements.iter().zip(self.materials.iter()) {
            if nodes.iter().any(|&i| i >= known_nodes) {
                continue;
            }
            if !material.mu.is_finite() || !material.kappa.is_finite() || !material.eta.is_finite()
            {
                continue;
            }
            let Some(rest) = TetRest::from_positions(nodes.map(|i| self.rest_positions[i])) else {
                continue;
            };
            elements.push(SoftTissueElement {
                nodes,
                rest,
                mu: material.mu,
                kappa: material.kappa,
                eta: material.eta.max(0.0),
            });
        }
        elements
    }
    /// Nodal forces of one element for the current positions and velocities,
    /// or `None` when the element is inverted (`det F <= 0`).
    ///
    /// The first Piola-Kirchhoff stress is
    /// `P = μ (F - F⁻ᵀ) + κ ln J F⁻ᵀ + η J sym(Ḟ F⁻¹) F⁻ᵀ`: the first two terms
    /// derive from [`NeoHookeanParams::strain_energy_density`], the last is the
    /// Newtonian viscous stress `η sym(L)` pulled back to the reference
    /// configuration. The nodal forces are `-V0 P Dm⁻ᵀ`.
    fn element_forces(&self, element: &SoftTissueElement) -> Option<[[f64; 3]; 4]> {
        let positions = element.nodes.map(|i| self.nodes[i]);
        let f = element.rest.gradient(&positions);
        let j = mat3_det(f);
        if !j.is_finite() || j <= 0.0 {
            return None;
        }
        let inv_j = 1.0 / j;
        let mut f_inv = mat3_adjugate(f);
        for entry in &mut f_inv {
            *entry *= inv_j;
        }
        let f_inv_t = mat3_transpose(f_inv);
        let ln_j = j.ln();
        let mut stress = [0.0; 9];
        for k in 0..9 {
            stress[k] = element.mu * (f[k] - f_inv_t[k]) + element.kappa * ln_j * f_inv_t[k];
        }
        if element.eta > 0.0 {
            let velocities = element.nodes.map(|i| self.velocities[i]);
            let f_dot = element.rest.gradient(&velocities);
            let rate = mat3_symmetric_part(mat3_mul(f_dot, f_inv));
            let viscous = mat3_mul(rate, f_inv_t);
            for k in 0..9 {
                stress[k] += element.eta * j * viscous[k];
            }
        }
        Some(element.rest.nodal_forces(&stress))
    }
    /// Internal force on every node \[N\] for the current positions and
    /// velocities: the neo-Hookean elastic forces plus the viscous forces of
    /// the materials' viscosity `η`.
    ///
    /// Inverted elements and elements without a usable rest shape or material
    /// do not contribute. Gravity and the Rayleigh damping `-α m v` are not
    /// included.
    pub fn internal_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0; 3]; self.nodes.len()];
        for element in &self.load_bearing_elements() {
            if let Some(nodal) = self.element_forces(element) {
                for (&node, force) in element.nodes.iter().zip(nodal.iter()) {
                    forces[node] = add3(forces[node], *force);
                }
            }
        }
        forces
    }
    /// Rest volume \[m³\] of element `e`; `0` if the element does not exist,
    /// refers to a node without rest position or has a degenerate rest shape.
    pub fn rest_volume(&self, e: usize) -> f64 {
        let Some(&nodes) = self.elements.get(e) else {
            return 0.0;
        };
        let mut rest = [[0.0; 3]; 4];
        for (slot, &node) in rest.iter_mut().zip(nodes.iter()) {
            match self.rest_positions.get(node) {
                Some(position) => *slot = *position,
                None => return 0.0,
            }
        }
        TetRest::from_positions(rest).map_or(0.0, |tet| tet.volume)
    }
    /// Total strain energy \[J\]: the energy density of every element that can
    /// exert forces (see [`SoftTissueModel::internal_forces`]) weighted with its
    /// rest volume, `Σ V0 W(F)`.
    ///
    /// Without viscosity the internal forces are the negative gradient of this
    /// energy with respect to the node positions. Inverted elements store no
    /// energy, just as they exert no force; elements that cannot exert forces
    /// (missing nodes, flat rest shape, non-finite material) are left out.
    pub fn strain_energy(&self) -> f64 {
        self.load_bearing_elements()
            .iter()
            .map(|element| {
                let material = NeoHookeanParams {
                    mu: element.mu,
                    kappa: element.kappa,
                    eta: element.eta,
                };
                let positions = element.nodes.map(|i| self.nodes[i]);
                let f = element.rest.gradient(&positions);
                element.rest.volume * material.strain_energy_density(&f)
            })
            .sum()
    }
    /// Estimate of the largest time step \[s\] for which explicit integration
    /// of this tissue is stable for small strains.
    ///
    /// It combines a rigorous bound on the highest natural frequency of the
    /// linearised neo-Hookean model with a bound on the damping rate of the
    /// viscosity `η` and of `α`, see [`SoftTissueModel::step`]. The viscous
    /// damping matrix does not commute with the stiffness matrix, so the
    /// combined limit is a sound estimate rather than a proof;
    /// [`SoftTissueModel::step`] keeps half of it as a safety margin. Infinite
    /// if nothing limits the step.
    pub fn critical_time_step(&self) -> f64 {
        self.stable_time_step(&self.load_bearing_elements())
    }
    fn stable_time_step(&self, elements: &[SoftTissueElement]) -> f64 {
        let known_nodes = self.known_node_count();
        let mut stiffness = vec![0.0; known_nodes];
        let mut viscosity = vec![0.0; known_nodes];
        for element in elements {
            // Small-strain limit of the neo-Hookean energy: Lamé parameters
            // (κ, μ); the Newtonian viscosity η acts like λ = 0, μ = η / 2.
            let elastic = element.rest.stiffness_row_sums(element.kappa, element.mu);
            let viscous = element.rest.stiffness_row_sums(0.0, 0.5 * element.eta);
            for ((&node, elastic_sum), viscous_sum) in
                element.nodes.iter().zip(elastic).zip(viscous)
            {
                stiffness[node] += elastic_sum;
                viscosity[node] += viscous_sum;
            }
        }
        let mut omega_sq = 0.0_f64;
        let mut viscous_rate = 0.0_f64;
        for node in 0..known_nodes {
            if self.fixed[node] || self.masses[node] <= 0.0 {
                continue;
            }
            let inv_mass = 1.0 / self.masses[node];
            omega_sq = omega_sq.max(stiffness[node] * inv_mass);
            viscous_rate = viscous_rate.max(viscosity[node] * inv_mass);
        }
        explicit_stability_limit(omega_sq.sqrt(), self.alpha.max(0.0) + viscous_rate)
    }
    /// Integrate the tissue forward by `dt` \[s\] with semi-implicit
    /// (symplectic) Euler.
    ///
    /// `gravity` is an acceleration \[m/s²\]. Fixed nodes and nodes without
    /// positive mass do not move; the acceleration of a free node of mass `m`
    /// is
    ///
    /// `a = f_internal / m + g - α v`
    ///
    /// with the neo-Hookean and viscous element forces of
    /// [`SoftTissueModel::internal_forces`].
    ///
    /// Explicit integration is only conditionally stable. The step is
    /// therefore split into the smallest number of equal sub-steps that stays
    /// below half of [`SoftTissueModel::critical_time_step`], capped at 4096
    /// sub-steps; beyond the cap the caller has to reduce `dt`, the stiffness
    /// or the viscosity `η` (for fine meshes the explicitly integrated viscous
    /// force limits the step to about `2 / γ`, where `γ` is the damping rate
    /// of `α` and `η`). Strongly compressed elements are stiffer than the
    /// bound assumes, which the safety margin only partly covers.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let known_nodes = self.known_node_count();
        let elements = self.load_bearing_elements();
        let substeps = stability_substeps(dt, self.stable_time_step(&elements));
        let h = dt / substeps as f64;
        let mut forces = vec![[0.0; 3]; known_nodes];
        for _ in 0..substeps {
            forces.fill([0.0; 3]);
            for element in &elements {
                if let Some(nodal) = self.element_forces(element) {
                    for (&node, force) in element.nodes.iter().zip(nodal.iter()) {
                        forces[node] = add3(forces[node], *force);
                    }
                }
            }
            for (i, force) in forces.iter().enumerate() {
                if self.fixed[i] || self.masses[i] <= 0.0 {
                    continue;
                }
                let accel = sub3(
                    add3(scale3(*force, 1.0 / self.masses[i]), gravity),
                    scale3(self.velocities[i], self.alpha),
                );
                self.velocities[i] = add3(self.velocities[i], scale3(accel, h));
                self.nodes[i] = add3(self.nodes[i], scale3(self.velocities[i], h));
            }
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
/// Precomputed data of an element of a [`SoftTissueModel`] that can exert forces.
struct SoftTissueElement {
    nodes: [usize; 4],
    rest: TetRest,
    mu: f64,
    kappa: f64,
    eta: f64,
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
