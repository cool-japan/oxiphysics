//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
/// One-way coupling: drives deformable body kinematics from a rigid body.
///
/// The rigid body's position and velocity are prescribed; the deformable mesh
/// nodes are constrained to follow via a damped spring.
#[derive(Debug, Clone)]
pub struct KinematicDriver {
    /// Rigid body world position (m).
    pub rigid_position: [f64; 3],
    /// Rigid body world velocity (m/s).
    pub rigid_velocity: [f64; 3],
    /// Spring stiffness for the kinematic constraint (N/m).
    pub drive_stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub drive_damping: f64,
    /// Attachment offset from rigid COM for each driven node, in local frame (m).
    pub attachment_offsets: Vec<[f64; 3]>,
}
impl KinematicDriver {
    /// Construct a kinematic driver for `n` nodes at the rigid COM (zero offsets).
    pub fn new_at_com(n: usize, drive_stiffness: f64, drive_damping: f64) -> Self {
        Self {
            rigid_position: [0.0; 3],
            rigid_velocity: [0.0; 3],
            drive_stiffness,
            drive_damping,
            attachment_offsets: vec![[0.0; 3]; n],
        }
    }
    /// Target world position for node `i` (m).
    pub fn target_position(&self, i: usize) -> [f64; 3] {
        vec3_add(self.rigid_position, self.attachment_offsets[i])
    }
    /// Driving force on node `i` given its current position and velocity.
    pub fn drive_force(&self, i: usize, node_pos: [f64; 3], node_vel: [f64; 3]) -> [f64; 3] {
        let target = self.target_position(i);
        let disp = vec3_sub(target, node_pos);
        let vel_err = vec3_sub(self.rigid_velocity, node_vel);
        let spring = vec3_scale(disp, self.drive_stiffness);
        let damp = vec3_scale(vel_err, self.drive_damping);
        vec3_add(spring, damp)
    }
    /// Update rigid body state.
    pub fn set_rigid_state(&mut self, position: [f64; 3], velocity: [f64; 3]) {
        self.rigid_position = position;
        self.rigid_velocity = velocity;
    }
}
/// Attaches a specific deformable mesh node to a rigid body with a spring constraint.
#[derive(Debug, Clone)]
pub struct CouplingJoint {
    /// Index of the attached mesh node.
    pub node_index: usize,
    /// Attachment point on the rigid body (world frame, m).
    pub anchor_world: [f64; 3],
    /// Spring stiffness (N/m).
    pub spring_stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Maximum force before the joint breaks (N); `f64::INFINITY` = unbreakable.
    pub max_force: f64,
    /// Whether the joint is intact.
    pub intact: bool,
}
impl CouplingJoint {
    /// Construct a new joint.
    pub fn new(
        node_index: usize,
        anchor_world: [f64; 3],
        spring_stiffness: f64,
        damping: f64,
    ) -> Self {
        Self {
            node_index,
            anchor_world,
            spring_stiffness,
            damping,
            max_force: f64::INFINITY,
            intact: true,
        }
    }
    /// Spring force on the node (N).
    pub fn spring_force(&self, node_pos: [f64; 3], node_vel: [f64; 3]) -> [f64; 3] {
        if !self.intact {
            return [0.0; 3];
        }
        let disp = vec3_sub(self.anchor_world, node_pos);
        let spring = vec3_scale(disp, self.spring_stiffness);
        let damp = vec3_scale(node_vel, -self.damping);
        vec3_add(spring, damp)
    }
    /// Check if the spring force exceeds `max_force` and break the joint if so.
    pub fn check_breakage(&mut self, node_pos: [f64; 3], node_vel: [f64; 3]) {
        let f = self.spring_force(node_pos, node_vel);
        if vec3_norm(f) > self.max_force {
            self.intact = false;
        }
    }
}
/// Transfers an impulse between a rigid body and a set of deformable mesh nodes.
///
/// The impulse is distributed to the `k` nearest nodes weighted by inverse
/// distance, and Newton's third law ensures the rigid body receives the
/// equal-and-opposite impulse.
#[derive(Debug, Clone)]
pub struct ImpulseTransfer {
    /// Rigid body inverse mass (1/kg); zero for static.
    pub rigid_inv_mass: f64,
    /// Per-node inverse masses (1/kg).
    pub node_inv_masses: Vec<f64>,
    /// Positions of mesh nodes (m), used for distance weighting.
    pub node_positions: Vec<[f64; 3]>,
}
impl ImpulseTransfer {
    /// Construct with a uniform node mass.
    pub fn new_uniform(node_positions: Vec<[f64; 3]>, node_mass: f64, rigid_mass: f64) -> Self {
        let n = node_positions.len();
        Self {
            rigid_inv_mass: if rigid_mass > 1e-12 {
                1.0 / rigid_mass
            } else {
                0.0
            },
            node_inv_masses: vec![1.0 / node_mass.max(1e-12); n],
            node_positions,
        }
    }
    /// Distribute `impulse` (N·s) from a contact point `point` to up to `k` nearest nodes.
    ///
    /// Returns per-node velocity deltas (m/s).
    pub fn distribute_to_nodes(
        &self,
        impulse: [f64; 3],
        point: [f64; 3],
        k: usize,
    ) -> Vec<[f64; 3]> {
        let n = self.node_positions.len();
        if n == 0 {
            return Vec::new();
        }
        let mut indexed: Vec<(usize, f64)> = self
            .node_positions
            .iter()
            .enumerate()
            .map(|(i, p)| (i, vec3_norm(vec3_sub(*p, point))))
            .collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let k_actual = k.min(n);
        let nearest = &indexed[..k_actual];
        let eps = 1e-6;
        let weights: Vec<f64> = nearest.iter().map(|(_, d)| 1.0 / (d + eps)).collect();
        let weight_sum: f64 = weights.iter().sum();
        let mut result = vec![[0.0_f64; 3]; n];
        for ((idx, _), w) in nearest.iter().zip(weights.iter()) {
            let frac = w / weight_sum.max(1e-12);
            let node_impulse = vec3_scale(impulse, frac);
            result[*idx] = vec3_scale(node_impulse, self.node_inv_masses[*idx]);
        }
        result
    }
    /// Velocity delta applied to the rigid body (equal-and-opposite to total impulse).
    pub fn rigid_velocity_delta(&self, impulse: [f64; 3]) -> [f64; 3] {
        vec3_scale(vec3_scale(impulse, -1.0), self.rigid_inv_mass)
    }
}
/// Method used to enforce the rigid-deformable coupling constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplingMethod {
    /// Penalty spring coupling — simple and fast, allows small penetration.
    PenaltySpring,
    /// Lagrange multiplier — exact, increases system DOF.
    LagrangeMultiplier,
    /// Mortar method — projection-based, best for sliding interfaces.
    Mortar,
}
/// Coupled rigid-deformable simulation step.
#[derive(Debug, Clone)]
pub struct CoupledSimulation {
    /// Rigid body position (m).
    pub rigid_position: [f64; 3],
    /// Rigid body velocity (m/s).
    pub rigid_velocity: [f64; 3],
    /// Rigid body mass (kg).
    pub mass: f64,
    /// Coupling bond.
    pub bond: RigidDeformableBond,
    /// Deformation state.
    pub deform: DeformationState,
    /// Thermal coupling.
    pub thermal: ThermalExpansionCoupling,
    /// Fluid-structure interface.
    pub fsi: FluidStructureInterface,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}
impl CoupledSimulation {
    /// Construct a default coupled simulation (in air — zero fluid pressure).
    pub fn new_default() -> Self {
        let mut fsi = FluidStructureInterface::submerged_plate();
        fsi.fluid_pressure = 0.0;
        Self {
            rigid_position: [0.0; 3],
            rigid_velocity: [0.0; 3],
            mass: 100.0,
            bond: RigidDeformableBond::default_rubber_seal(),
            deform: DeformationState::new_undeformed(16),
            thermal: ThermalExpansionCoupling::aluminium_with_rubber_seal(),
            fsi,
            gravity: [0.0, -9.81, 0.0],
        }
    }
    /// Advance the simulation by one time step `dt` (s).
    ///
    /// Steps:
    /// 1. Compute coupling forces (penalty spring + FSI + thermal).
    /// 2. Integrate rigid body (semi-implicit Euler).
    /// 3. Update deformable mesh nodal forces.
    pub fn step(&mut self, dt: f64, external_force: [f64; 3]) {
        if !self.bond.active {
            return;
        }
        let thermal_f = self.thermal.thermal_force();
        let fsi_f = self.fsi.pressure_force();
        self.fsi.update_rigid_velocity(self.rigid_velocity);
        let max_disp = self.deform.max_displacement();
        let coupling_f = self.bond.penalty_stiffness * max_disp;
        let gravity_f = vec3_scale(self.gravity, self.mass);
        let total_f = vec3_add(
            vec3_add(external_force, gravity_f),
            vec3_add(fsi_f, [thermal_f, 0.0, coupling_f]),
        );
        let accel = vec3_scale(total_f, 1.0 / self.mass.max(1e-6));
        self.rigid_velocity = vec3_add(self.rigid_velocity, vec3_scale(accel, dt));
        self.rigid_position = vec3_add(self.rigid_position, vec3_scale(self.rigid_velocity, dt));
        let n = self.deform.nodal_forces.len();
        let nodal_f_val = coupling_f / n.max(1) as f64;
        for nf in &mut self.deform.nodal_forces {
            *nf = [nodal_f_val, 0.0, 0.0];
        }
        let force_mag = vec3_norm(total_f);
        if force_mag > self.bond.max_force {
            self.bond.debond();
        }
    }
    /// Kinetic energy of the rigid body (J).
    pub fn rigid_kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.rigid_velocity, self.rigid_velocity)
    }
}
/// Couples a rigid body to a softbody cluster via the shape-matching algorithm.
///
/// The optimal rigid transform (translation + rotation) that maps the rest
/// positions to the current positions is computed.  The rigid body is then
/// driven towards that extracted transform.
#[derive(Debug, Clone)]
pub struct ShapeMatchingCoupling {
    /// Rest-pose positions of the cluster nodes (m).
    pub rest_positions: Vec<[f64; 3]>,
    /// Per-node mass weights (kg).
    pub node_masses: Vec<f64>,
    /// Blending factor α ∈ \[0, 1\]: 0 = fully elastic, 1 = fully rigid.
    pub stiffness_alpha: f64,
    /// Cluster goal positions updated each step.
    pub goal_positions: Vec<[f64; 3]>,
}
impl ShapeMatchingCoupling {
    /// Construct with uniform node masses.
    pub fn new_uniform(
        rest_positions: Vec<[f64; 3]>,
        node_mass: f64,
        stiffness_alpha: f64,
    ) -> Self {
        let n = rest_positions.len();
        let goals = rest_positions.clone();
        Self {
            rest_positions,
            node_masses: vec![node_mass; n],
            stiffness_alpha,
            goal_positions: goals,
        }
    }
    /// Total cluster mass (kg).
    pub fn total_mass(&self) -> f64 {
        self.node_masses.iter().sum()
    }
    /// Centre of mass of the rest configuration (m).
    pub fn rest_com(&self) -> [f64; 3] {
        let total = self.total_mass().max(1e-12);
        let mut com = [0.0; 3];
        for (p, m) in self.rest_positions.iter().zip(self.node_masses.iter()) {
            com = vec3_add(com, vec3_scale(*p, *m));
        }
        vec3_scale(com, 1.0 / total)
    }
    /// Update goal positions by blending current toward rigid-driven targets.
    ///
    /// `current` — current node positions, `targets` — rigid-driven targets.
    pub fn update_goals(&mut self, current: &[[f64; 3]], targets: &[[f64; 3]]) {
        let a = self.stiffness_alpha;
        self.goal_positions = current
            .iter()
            .zip(targets.iter())
            .map(|(c, t)| vec3_add(vec3_scale(*c, 1.0 - a), vec3_scale(*t, a)))
            .collect();
    }
    /// Restoring force on each node towards its goal (N/m per unit stiffness).
    pub fn restoring_forces(&self, current: &[[f64; 3]], stiffness: f64) -> Vec<[f64; 3]> {
        current
            .iter()
            .zip(self.goal_positions.iter())
            .map(|(c, g)| vec3_scale(vec3_sub(*g, *c), stiffness))
            .collect()
    }
}
/// Hertz contact + plastic deformation: crater depth model.
#[derive(Debug, Clone)]
pub struct ImpactDeformation {
    /// Effective elastic modulus E* (Pa).
    pub effective_modulus: f64,
    /// Yield strength of the deformable body (Pa).
    pub yield_strength: f64,
    /// Impactor mass (kg).
    pub impactor_mass: f64,
    /// Impactor radius (m).
    pub impactor_radius: f64,
    /// Current plastic indentation depth (m).
    pub plastic_depth: f64,
}
impl ImpactDeformation {
    /// Default steel-on-steel impact.
    pub fn steel_on_steel() -> Self {
        Self {
            effective_modulus: 115e9,
            yield_strength: 250e6,
            impactor_mass: 0.5,
            impactor_radius: 0.01,
            plastic_depth: 0.0,
        }
    }
    /// Default rubber ball on aluminium plate.
    pub fn rubber_on_aluminium() -> Self {
        Self {
            effective_modulus: 3e9,
            yield_strength: 100e6,
            impactor_mass: 0.05,
            impactor_radius: 0.025,
            plastic_depth: 0.0,
        }
    }
    /// Peak Hertz contact force (N) for a given approach velocity (m/s).
    pub fn peak_hertz_force(&self, velocity: f64) -> f64 {
        let e = self.effective_modulus;
        let r = self.impactor_radius;
        let m = self.impactor_mass;
        (2.0 * m * e * r * velocity * velocity).sqrt()
    }
    /// Elastic contact radius (m) at peak force.
    pub fn hertz_contact_radius(&self, velocity: f64) -> f64 {
        let force = self.peak_hertz_force(velocity);
        let numer = 3.0 * force * self.impactor_radius;
        let denom = 4.0 * self.effective_modulus;
        (numer / denom).powf(1.0 / 3.0)
    }
    /// Plastic crater depth (m) for a given impact velocity (m/s).
    ///
    /// Uses simplified strain-energy balance.
    pub fn plastic_crater_depth(&mut self, velocity: f64) -> f64 {
        let fy = self.yield_strength;
        let ke = 0.5 * self.impactor_mass * velocity * velocity;
        let area = std::f64::consts::PI * self.impactor_radius * self.impactor_radius;
        let depth = ke / (fy * area).max(1e-12);
        self.plastic_depth = depth;
        depth
    }
    /// Returns `true` if the impact has caused yielding.
    pub fn caused_yielding(&self, velocity: f64) -> bool {
        let force = self.peak_hertz_force(velocity);
        let area = std::f64::consts::PI * self.impactor_radius * self.impactor_radius;
        force / area.max(1e-12) > self.yield_strength
    }
}
/// State of the deformable mesh at a coupling interface.
#[derive(Debug, Clone)]
pub struct DeformationState {
    /// Nodal displacements at the interface nodes (m).
    pub displacements: Vec<[f64; 3]>,
    /// Nodal velocities (m/s).
    pub velocities: Vec<[f64; 3]>,
    /// Equivalent nodal forces from coupling (N).
    pub nodal_forces: Vec<[f64; 3]>,
    /// Von Mises strain at each node (dimensionless).
    pub strains: Vec<f64>,
}
impl DeformationState {
    /// Create an undeformed state with `n` nodes.
    pub fn new_undeformed(n: usize) -> Self {
        Self {
            displacements: vec![[0.0; 3]; n],
            velocities: vec![[0.0; 3]; n],
            nodal_forces: vec![[0.0; 3]; n],
            strains: vec![0.0; n],
        }
    }
    /// Apply a uniform displacement to all nodes.
    pub fn apply_uniform_displacement(&mut self, disp: [f64; 3]) {
        for d in &mut self.displacements {
            *d = vec3_add(*d, disp);
        }
    }
    /// Maximum nodal displacement (m).
    pub fn max_displacement(&self) -> f64 {
        self.displacements
            .iter()
            .map(|d| vec3_norm(*d))
            .fold(0.0_f64, f64::max)
    }
    /// Update nodal forces from an external force vector.
    pub fn set_nodal_forces(&mut self, forces: Vec<[f64; 3]>) {
        self.nodal_forces = forces;
    }
    /// RMS nodal strain.
    pub fn rms_strain(&self) -> f64 {
        if self.strains.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = self.strains.iter().map(|s| s * s).sum();
        (sum_sq / self.strains.len() as f64).sqrt()
    }
}
/// Linear blend skinning (LBS) weight matrix for a deformable proxy.
///
/// Each row corresponds to a mesh node; each column to a rigid bone/joint.
/// `weights[node][bone]` ∈ \[0, 1\] and each row sums to 1.
#[derive(Debug, Clone)]
pub struct SkinningMatrix {
    /// Weight matrix: `weights[node_index][bone_index]`.
    pub weights: Vec<Vec<f64>>,
    /// Number of mesh nodes.
    pub num_nodes: usize,
    /// Number of bones (rigid bodies in a chain).
    pub num_bones: usize,
}
impl SkinningMatrix {
    /// Construct with all weight on the first bone (rigid attachment).
    pub fn new_rigid(num_nodes: usize, num_bones: usize) -> Self {
        let mut weights = vec![vec![0.0; num_bones]; num_nodes];
        for row in &mut weights {
            if num_bones > 0 {
                row[0] = 1.0;
            }
        }
        Self {
            weights,
            num_nodes,
            num_bones,
        }
    }
    /// Construct with uniform weights across all bones.
    pub fn new_uniform(num_nodes: usize, num_bones: usize) -> Self {
        let w = if num_bones > 0 {
            1.0 / num_bones as f64
        } else {
            0.0
        };
        Self {
            weights: vec![vec![w; num_bones]; num_nodes],
            num_nodes,
            num_bones,
        }
    }
    /// Normalise each row so weights sum to 1.
    pub fn normalise_rows(&mut self) {
        for row in &mut self.weights {
            let s: f64 = row.iter().sum();
            if s > 1e-12 {
                for w in row.iter_mut() {
                    *w /= s;
                }
            }
        }
    }
    /// Compute skinned position of node `i` given `bone_positions` (m).
    pub fn skinned_position(&self, node_index: usize, bone_positions: &[[f64; 3]]) -> [f64; 3] {
        let mut pos = [0.0; 3];
        for (b, &w) in self.weights[node_index].iter().enumerate() {
            if b < bone_positions.len() {
                pos = vec3_add(pos, vec3_scale(bone_positions[b], w));
            }
        }
        pos
    }
    /// Row sum for node `i` (should be ≈1 after normalisation).
    pub fn row_sum(&self, node_index: usize) -> f64 {
        self.weights[node_index].iter().sum()
    }
}
/// Bond parameters connecting a rigid body to a deformable mesh at a surface.
#[derive(Debug, Clone)]
pub struct RigidDeformableBond {
    /// Rigid body ID (external reference).
    pub rigid_body_id: usize,
    /// Coupling method.
    pub method: CouplingMethod,
    /// Penalty spring stiffness (N/m).
    pub penalty_stiffness: f64,
    /// Penalty damping (N·s/m).
    pub penalty_damping: f64,
    /// Maximum allowable coupling force before debond (N).
    pub max_force: f64,
    /// Whether the bond is active.
    pub active: bool,
    /// Number of contact points in this bond.
    pub num_contact_points: usize,
}
impl RigidDeformableBond {
    /// Default penalty coupling for a rubber seal on a steel plate.
    pub fn default_rubber_seal() -> Self {
        Self {
            rigid_body_id: 0,
            method: CouplingMethod::PenaltySpring,
            penalty_stiffness: 1e6,
            penalty_damping: 1e4,
            max_force: 50_000.0,
            active: true,
            num_contact_points: 16,
        }
    }
    /// Mortar coupling for structural steel.
    pub fn structural_mortar(rigid_body_id: usize) -> Self {
        Self {
            rigid_body_id,
            method: CouplingMethod::Mortar,
            penalty_stiffness: 5e8,
            penalty_damping: 5e5,
            max_force: 500_000.0,
            active: true,
            num_contact_points: 64,
        }
    }
    /// Deactivate (debond) the coupling.
    pub fn debond(&mut self) {
        self.active = false;
    }
}
/// Coupling force computed at a contact point.
#[derive(Debug, Clone)]
pub struct CouplingForce {
    /// Contact point world position (m).
    pub point: [f64; 3],
    /// Force vector applied on the rigid body (N).
    pub rigid_force: [f64; 3],
    /// Equal and opposite force applied to the deformable mesh node (N).
    pub deformable_force: [f64; 3],
    /// Coupling method used.
    pub method: CouplingMethod,
}
impl CouplingForce {
    /// Compute a penalty spring coupling force.
    ///
    /// `projection` — surface projection result.
    /// `penalty_k` — penalty stiffness (N/m).
    /// `damping_c` — penalty damping (N·s/m).
    /// `relative_velocity` — velocity of rigid point relative to mesh (m/s).
    pub fn penalty(
        projection: &SurfaceProjection,
        penalty_k: f64,
        damping_c: f64,
        relative_velocity: f64,
    ) -> Self {
        let penetration = (-projection.gap).max(0.0);
        let spring_force = penalty_k * penetration;
        let damping_force = damping_c * relative_velocity.abs();
        let magnitude = spring_force + damping_force;
        let f = vec3_scale(projection.normal, magnitude);
        CouplingForce {
            point: projection.rigid_point,
            rigid_force: f,
            deformable_force: vec3_scale(f, -1.0),
            method: CouplingMethod::PenaltySpring,
        }
    }
    /// Magnitude of the coupling force (N).
    pub fn magnitude(&self) -> f64 {
        vec3_norm(self.rigid_force)
    }
}
/// Projection of a rigid-body surface point onto a deformable mesh.
#[derive(Debug, Clone)]
pub struct SurfaceProjection {
    /// Rigid-body surface point in world coordinates (m).
    pub rigid_point: [f64; 3],
    /// Projected point on the deformable mesh surface (m).
    pub projected_point: [f64; 3],
    /// Index of the nearest mesh triangle.
    pub triangle_idx: usize,
    /// Barycentric coordinates within that triangle.
    pub barycentric: [f64; 3],
    /// Surface normal at the projected point (unit vector).
    pub normal: [f64; 3],
    /// Gap (signed distance, positive = not penetrating) (m).
    pub gap: f64,
}
impl SurfaceProjection {
    /// Compute a projection from a rigid point to a triangle defined by
    /// three vertices `v0`, `v1`, `v2`.
    pub fn project_to_triangle(
        rigid_point: [f64; 3],
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        triangle_idx: usize,
    ) -> Self {
        let e1 = vec3_sub(v1, v0);
        let e2 = vec3_sub(v2, v0);
        let normal_raw = [
            e2[1] * e1[2] - e2[2] * e1[1],
            e2[2] * e1[0] - e2[0] * e1[2],
            e2[0] * e1[1] - e2[1] * e1[0],
        ];
        let normal = vec3_normalise(normal_raw);
        let delta = vec3_sub(rigid_point, v0);
        let gap = vec3_dot(delta, normal);
        let projected = vec3_sub(rigid_point, vec3_scale(normal, gap));
        let area2 = vec3_norm(normal_raw);
        let bary = if area2 > 1e-12 {
            let c0 = [
                (v1[1] - v2[1]) * (projected[0] - v2[0]) + (v2[0] - v1[0]) * (projected[1] - v2[1]),
                (v2[1] - v0[1]) * (projected[0] - v0[0]) + (v0[0] - v2[0]) * (projected[1] - v0[1]),
                0.0,
            ];
            let w0 = c0[0] / area2;
            let w1 = c0[1] / area2;
            let w2 = 1.0 - w0 - w1;
            [w0.max(0.0), w1.max(0.0), w2.max(0.0)]
        } else {
            [1.0 / 3.0; 3]
        };
        SurfaceProjection {
            rigid_point,
            projected_point: projected,
            triangle_idx,
            barycentric: bary,
            normal,
            gap,
        }
    }
    /// Returns `true` if the rigid point is penetrating the mesh surface.
    pub fn is_penetrating(&self) -> bool {
        self.gap < 0.0
    }
}
/// Rigid body acting as a proxy for a deformable surface.
///
/// Stores the reference (rest) shape and a per-node deformation map so that
/// a rigid-body solver can drive (or be driven by) a soft-body simulation.
#[derive(Debug, Clone)]
pub struct DeformableProxy {
    /// Rest-pose node positions (m), one entry per surface node.
    pub reference_shape: Vec<[f64; 3]>,
    /// Current deformed positions of the surface nodes (m).
    pub deformed_shape: Vec<[f64; 3]>,
    /// Rigid-body centre of mass (m).
    pub rigid_com: [f64; 3],
    /// Rigid-body orientation as a 3×3 rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Proxy stiffness: restoring spring from deformed to rigid-driven position (N/m).
    pub proxy_stiffness: f64,
}
impl DeformableProxy {
    /// Construct a proxy from a reference shape with identity rotation.
    pub fn new(reference_shape: Vec<[f64; 3]>, proxy_stiffness: f64) -> Self {
        let deformed = reference_shape.clone();
        Self {
            reference_shape,
            deformed_shape: deformed,
            rigid_com: [0.0; 3],
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            proxy_stiffness,
        }
    }
    /// Apply the rigid-body transform (rotation R + translation t) to the
    /// reference shape to obtain the rigid-driven target positions.
    pub fn rigid_driven_positions(&self) -> Vec<[f64; 3]> {
        let r = &self.rotation;
        self.reference_shape
            .iter()
            .map(|p| {
                let rp = [
                    r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
                    r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
                    r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
                ];
                vec3_add(rp, self.rigid_com)
            })
            .collect()
    }
    /// Compute per-node spring forces pulling deformed nodes towards rigid targets (N).
    pub fn spring_forces(&self) -> Vec<[f64; 3]> {
        let targets = self.rigid_driven_positions();
        self.deformed_shape
            .iter()
            .zip(targets.iter())
            .map(|(d, t)| vec3_scale(vec3_sub(*t, *d), self.proxy_stiffness))
            .collect()
    }
    /// Number of surface nodes.
    pub fn num_nodes(&self) -> usize {
        self.reference_shape.len()
    }
    /// Maximum deviation between deformed and rest shape (m).
    pub fn max_deformation(&self) -> f64 {
        self.reference_shape
            .iter()
            .zip(self.deformed_shape.iter())
            .map(|(r, d)| vec3_norm(vec3_sub(*d, *r)))
            .fold(0.0_f64, f64::max)
    }
}
/// Mutual force exchange between a rigid body and a deformable mesh.
///
/// The deformable mesh exerts a reaction force on the rigid body equal to
/// the negative of the sum of all nodal coupling forces.
#[derive(Debug, Clone)]
pub struct TwoWayCoupling {
    /// Rigid body mass (kg).
    pub rigid_mass: f64,
    /// Rigid body velocity (m/s).
    pub rigid_velocity: [f64; 3],
    /// Spring stiffness of the coupling (N/m).
    pub coupling_stiffness: f64,
    /// Per-node coupling spring rest lengths (m).
    pub rest_lengths: Vec<f64>,
    /// Per-node attachment points on the rigid body in world frame (m).
    pub attachment_points: Vec<[f64; 3]>,
}
impl TwoWayCoupling {
    /// Construct for `n` nodes, all attached at the rigid COM with zero rest length.
    pub fn new(n: usize, rigid_mass: f64, coupling_stiffness: f64) -> Self {
        Self {
            rigid_mass,
            rigid_velocity: [0.0; 3],
            coupling_stiffness,
            rest_lengths: vec![0.0; n],
            attachment_points: vec![[0.0; 3]; n],
        }
    }
    /// Compute coupling force on a node at `node_pos` connected to attachment point `i` (N).
    pub fn node_force(&self, i: usize, node_pos: [f64; 3]) -> [f64; 3] {
        let attach = self.attachment_points[i];
        let delta = vec3_sub(attach, node_pos);
        let dist = vec3_norm(delta);
        let rest = self.rest_lengths[i];
        let extension = dist - rest;
        if dist < 1e-12 || extension <= 0.0 {
            return [0.0; 3];
        }
        vec3_scale(vec3_normalise(delta), extension * self.coupling_stiffness)
    }
    /// Reaction force on the rigid body: negative sum of all nodal forces (N).
    pub fn rigid_reaction_force(&self, node_positions: &[[f64; 3]]) -> [f64; 3] {
        let total_node_force: [f64; 3] = node_positions
            .iter()
            .enumerate()
            .fold([0.0; 3], |acc, (i, p)| {
                vec3_add(acc, self.node_force(i, *p))
            });
        vec3_scale(total_node_force, -1.0)
    }
    /// Integrate rigid body velocity under reaction force for time step `dt` (s).
    pub fn integrate_rigid(&mut self, node_positions: &[[f64; 3]], dt: f64) {
        let f = self.rigid_reaction_force(node_positions);
        let accel = vec3_scale(f, 1.0 / self.rigid_mass.max(1e-12));
        self.rigid_velocity = vec3_add(self.rigid_velocity, vec3_scale(accel, dt));
    }
}
/// Enforces velocity compatibility at the rigid-deformable interface.
///
/// Projects the deformable node velocity onto the rigid body velocity at the
/// attachment point (Baumgarte stabilisation).
#[derive(Debug, Clone)]
pub struct VelocityConstraint {
    /// Node index in the deformable mesh.
    pub node_index: usize,
    /// Rigid body velocity (m/s).
    pub rigid_velocity: [f64; 3],
    /// Constraint normal (unit vector, direction of constraint).
    pub normal: [f64; 3],
    /// Baumgarte stabilisation coefficient β ∈ \[0, 1\].
    pub baumgarte_beta: f64,
    /// Position error (m) — drift from the constraint manifold.
    pub position_error: f64,
}
impl VelocityConstraint {
    /// Construct a new velocity constraint.
    pub fn new(node_index: usize, normal: [f64; 3], baumgarte_beta: f64) -> Self {
        Self {
            node_index,
            rigid_velocity: [0.0; 3],
            normal: vec3_normalise(normal),
            baumgarte_beta,
            position_error: 0.0,
        }
    }
    /// Relative velocity (m/s) along the constraint normal.
    pub fn relative_velocity(&self, node_vel: [f64; 3]) -> f64 {
        vec3_dot(vec3_sub(node_vel, self.rigid_velocity), self.normal)
    }
    /// Constraint violation (scalar): velocity mismatch + Baumgarte correction.
    pub fn violation(&self, node_vel: [f64; 3], dt: f64) -> f64 {
        self.relative_velocity(node_vel) + self.baumgarte_beta * self.position_error / dt.max(1e-12)
    }
    /// Lagrange multiplier λ for the constraint force (N·s/m²).
    ///
    /// `node_inv_mass` — inverse mass of the node (1/kg).
    pub fn lagrange_multiplier(&self, node_vel: [f64; 3], node_inv_mass: f64, dt: f64) -> f64 {
        let c_dot = self.violation(node_vel, dt);
        -c_dot / node_inv_mass.max(1e-12)
    }
}
/// Thermal expansion coupling: temperature-dependent rigid body size change.
#[derive(Debug, Clone)]
pub struct ThermalExpansionCoupling {
    /// Coefficient of thermal expansion (1/K).
    pub cte: f64,
    /// Reference temperature (K).
    pub reference_temp: f64,
    /// Current temperature (K).
    pub current_temp: f64,
    /// Nominal length (m) of the rigid body in the constrained direction.
    pub nominal_length: f64,
    /// Stiffness of the deformable mesh at the interface (N/m).
    pub interface_stiffness: f64,
}
impl ThermalExpansionCoupling {
    /// Default aluminium alloy with rubber seal.
    pub fn aluminium_with_rubber_seal() -> Self {
        Self {
            cte: 23.1e-6,
            reference_temp: 293.15,
            current_temp: 293.15,
            nominal_length: 0.5,
            interface_stiffness: 1e7,
        }
    }
    /// Thermal strain (dimensionless).
    pub fn thermal_strain(&self) -> f64 {
        self.cte * (self.current_temp - self.reference_temp)
    }
    /// Thermal expansion displacement (m).
    pub fn expansion(&self) -> f64 {
        self.thermal_strain() * self.nominal_length
    }
    /// Interface force due to thermal expansion (N).
    pub fn thermal_force(&self) -> f64 {
        self.expansion() * self.interface_stiffness
    }
    /// Update temperature.
    pub fn set_temperature(&mut self, temp_k: f64) {
        self.current_temp = temp_k;
    }
}
/// Tracks energy transfer between the rigid body and the deformable subsystem.
#[derive(Debug, Clone)]
pub struct EnergyExchange {
    /// Cumulative energy injected into the deformable system from the rigid body (J).
    pub energy_to_deformable: f64,
    /// Cumulative energy returned from the deformable system to the rigid body (J).
    pub energy_to_rigid: f64,
    /// Current kinetic energy of the rigid body (J).
    pub rigid_ke: f64,
    /// Current elastic strain energy of the deformable body (J).
    pub deformable_strain_energy: f64,
}
impl EnergyExchange {
    /// Construct with zero energy state.
    pub fn new() -> Self {
        Self {
            energy_to_deformable: 0.0,
            energy_to_rigid: 0.0,
            rigid_ke: 0.0,
            deformable_strain_energy: 0.0,
        }
    }
    /// Record a coupling power exchange over time step `dt` (s).
    ///
    /// Positive `power_to_deformable` means energy flows rigid → deformable.
    pub fn record_exchange(&mut self, power_to_deformable: f64, dt: f64) {
        if power_to_deformable >= 0.0 {
            self.energy_to_deformable += power_to_deformable * dt;
        } else {
            self.energy_to_rigid += (-power_to_deformable) * dt;
        }
    }
    /// Update the instantaneous KE and strain energy values.
    pub fn update(&mut self, rigid_ke: f64, strain_energy: f64) {
        self.rigid_ke = rigid_ke;
        self.deformable_strain_energy = strain_energy;
    }
    /// Total energy in the coupled system (J).
    pub fn total_energy(&self) -> f64 {
        self.rigid_ke + self.deformable_strain_energy
    }
    /// Net energy transferred (positive = net flow to deformable) (J).
    pub fn net_transfer(&self) -> f64 {
        self.energy_to_deformable - self.energy_to_rigid
    }
}
/// Detects and resolves contacts between a rigid body surface and a deformable mesh.
#[derive(Debug, Clone)]
pub struct ContactBridge {
    /// Contact detection tolerance (m): nodes within this distance are in contact.
    pub contact_threshold: f64,
    /// Normal restitution coefficient e ∈ \[0, 1\].
    pub restitution: f64,
    /// Friction coefficient μ.
    pub friction: f64,
    /// Rigid body position (m).
    pub rigid_position: [f64; 3],
    /// Rigid body radius (simplified sphere shape) (m).
    pub rigid_radius: f64,
}
impl ContactBridge {
    /// Construct a default contact bridge.
    pub fn new(contact_threshold: f64, restitution: f64, friction: f64) -> Self {
        Self {
            contact_threshold,
            restitution,
            friction,
            rigid_position: [0.0; 3],
            rigid_radius: 0.5,
        }
    }
    /// Detect which mesh nodes are in contact with the rigid sphere.
    ///
    /// Returns indices and penetration depths (m, positive = penetrating).
    pub fn detect_contacts(&self, node_positions: &[[f64; 3]]) -> Vec<(usize, f64)> {
        node_positions
            .iter()
            .enumerate()
            .filter_map(|(i, p)| {
                let dist = vec3_norm(vec3_sub(*p, self.rigid_position));
                let penetration = self.rigid_radius + self.contact_threshold - dist;
                if penetration > 0.0 {
                    Some((i, penetration))
                } else {
                    None
                }
            })
            .collect()
    }
    /// Compute contact normal impulse for a contacting node.
    ///
    /// `v_rel` — relative normal velocity (m/s, positive = approaching).
    /// `node_inv_mass` — node inverse mass (1/kg).
    /// Returns impulse magnitude (N·s).
    pub fn contact_impulse(&self, v_rel: f64, node_inv_mass: f64) -> f64 {
        if v_rel <= 0.0 {
            return 0.0;
        }
        let e = self.restitution;
        let inv_m_rigid = 1.0_f64;
        let denom = inv_m_rigid + node_inv_mass;
        if denom < 1e-12 {
            return 0.0;
        }
        -(1.0 + e) * v_rel / denom
    }
    /// Number of contacts detected.
    pub fn contact_count(&self, node_positions: &[[f64; 3]]) -> usize {
        self.detect_contacts(node_positions).len()
    }
}
/// Transmission of rigid-body vibration to a structural frequency response.
#[derive(Debug, Clone)]
pub struct VibrationTransmission {
    /// Natural frequency of the structural mode (Hz).
    pub natural_freq: f64,
    /// Damping ratio of the structural mode.
    pub damping_ratio: f64,
    /// Coupling mass (modal mass of the interface) (kg).
    pub modal_mass: f64,
    /// Coupling stiffness (N/m).
    pub coupling_stiffness: f64,
}
impl VibrationTransmission {
    /// Default for a steel panel attached to a rigid machine frame.
    pub fn steel_panel() -> Self {
        Self {
            natural_freq: 150.0,
            damping_ratio: 0.02,
            modal_mass: 5.0,
            coupling_stiffness: 1e7,
        }
    }
    /// Frequency response function (FRF) amplitude at excitation frequency `f` (Hz).
    ///
    /// Returns the structural displacement per unit force amplitude (m/N).
    pub fn frf_amplitude(&self, f: f64) -> f64 {
        let fn_ = self.natural_freq.max(1e-6);
        let zeta = self.damping_ratio;
        let k = self.coupling_stiffness.max(1.0);
        let r = f / fn_;
        let denom = ((1.0 - r * r).powi(2) + (2.0 * zeta * r).powi(2)).sqrt();
        1.0 / (k * denom.max(1e-12))
    }
    /// Peak FRF amplitude at resonance (m/N).
    pub fn peak_frf(&self) -> f64 {
        let k = self.coupling_stiffness.max(1.0);
        let zeta = self.damping_ratio.max(1e-6);
        1.0 / (2.0 * k * zeta)
    }
    /// RMS response acceleration (m/s²) given input RMS force (N) at frequency `f`.
    pub fn rms_acceleration(&self, f: f64, rms_force: f64) -> f64 {
        let disp_frf = self.frf_amplitude(f);
        let omega = 2.0 * std::f64::consts::PI * f;
        disp_frf * rms_force * omega * omega
    }
}
/// Fluid-structure interface: pressure loads on rigid surface.
#[derive(Debug, Clone)]
pub struct FluidStructureInterface {
    /// Surface area of the interface (m²).
    pub surface_area: f64,
    /// Current fluid pressure at the interface (Pa).
    pub fluid_pressure: f64,
    /// Surface normal (outward from fluid) (unit vector).
    pub surface_normal: [f64; 3],
    /// Rigid body position (m) — updated from simulation.
    pub rigid_position: [f64; 3],
    /// Rigid body velocity (m/s).
    pub rigid_velocity: [f64; 3],
    /// Added-mass coefficient (fluid-induced inertia multiplier).
    pub added_mass_coeff: f64,
    /// Structural mass (kg).
    pub structural_mass: f64,
}
impl FluidStructureInterface {
    /// Default for a submerged rigid plate in water.
    pub fn submerged_plate() -> Self {
        Self {
            surface_area: 0.25,
            fluid_pressure: 101_325.0,
            surface_normal: [0.0, 1.0, 0.0],
            rigid_position: [0.0; 3],
            rigid_velocity: [0.0; 3],
            added_mass_coeff: 1.0,
            structural_mass: 20.0,
        }
    }
    /// Pressure force on the rigid surface (N) along the surface normal.
    pub fn pressure_force(&self) -> [f64; 3] {
        let magnitude = self.fluid_pressure * self.surface_area;
        vec3_scale(self.surface_normal, magnitude)
    }
    /// Added mass (kg) contributing to rigid-body inertia in the normal direction.
    pub fn added_mass(&self) -> f64 {
        let rho_water = 1000.0_f64;
        let volume_approx = self.surface_area * 0.1;
        rho_water * volume_approx * self.added_mass_coeff
    }
    /// Update rigid body velocity (feedback from solver).
    pub fn update_rigid_velocity(&mut self, velocity: [f64; 3]) {
        self.rigid_velocity = velocity;
    }
    /// Compute the kinematic feedback pressure correction to the fluid (Pa).
    ///
    /// When the rigid surface moves into the fluid, the pressure rises locally.
    pub fn kinematic_feedback_pressure(&self) -> f64 {
        let rho_water = 1000.0_f64;
        let vn = vec3_dot(self.rigid_velocity, self.surface_normal);
        0.5 * rho_water * vn * vn
    }
}
/// Adhesive layer model for soft attachment between rigid and deformable bodies.
#[derive(Debug, Clone)]
pub struct SoftAttachment {
    /// Normal (mode-I) stiffness per unit area (N/m³).
    pub normal_stiffness: f64,
    /// Shear (mode-II) stiffness per unit area (N/m³).
    pub shear_stiffness: f64,
    /// Normal debonding strength (Pa).
    pub normal_strength: f64,
    /// Shear debonding strength (Pa).
    pub shear_strength: f64,
    /// Contact area (m²).
    pub contact_area: f64,
    /// Whether the attachment has debonded.
    pub debonded: bool,
}
impl SoftAttachment {
    /// Default epoxy adhesive bond.
    pub fn epoxy_bond(contact_area: f64) -> Self {
        Self {
            normal_stiffness: 1e9,
            shear_stiffness: 5e8,
            normal_strength: 20e6,
            shear_strength: 15e6,
            contact_area,
            debonded: false,
        }
    }
    /// Default rubber gasket.
    pub fn rubber_gasket(contact_area: f64) -> Self {
        Self {
            normal_stiffness: 1e7,
            shear_stiffness: 5e6,
            normal_strength: 2e6,
            shear_strength: 1.5e6,
            contact_area,
            debonded: false,
        }
    }
    /// Compute normal and shear coupling forces given displacements (m).
    ///
    /// Returns `(normal_force, shear_force)` in N and checks for debonding.
    pub fn coupling_force(&mut self, normal_disp: f64, shear_disp: f64) -> (f64, f64) {
        if self.debonded {
            return (0.0, 0.0);
        }
        let fn_ = self.normal_stiffness * normal_disp * self.contact_area;
        let fs = self.shear_stiffness * shear_disp * self.contact_area;
        let sigma_n = fn_ / self.contact_area.max(1e-12);
        let tau = fs.abs() / self.contact_area.max(1e-12);
        if sigma_n > self.normal_strength || tau > self.shear_strength {
            self.debonded = true;
            return (0.0, 0.0);
        }
        (fn_, fs)
    }
}
