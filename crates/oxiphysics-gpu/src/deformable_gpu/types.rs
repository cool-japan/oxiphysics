//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    add3, cross3, dot3, len2_3, length3, normalize3, quat_rotate, scale3, sub3,
};

/// An XPBD distance constraint between two particles.
#[derive(Debug, Clone)]
pub struct XpbdDistConstraint {
    /// Index of particle A.
    pub i: usize,
    /// Index of particle B.
    pub j: usize,
    /// Rest length.
    pub rest_len: f64,
    /// Compliance (inverse stiffness) in m/N.
    pub compliance: f64,
}
impl XpbdDistConstraint {
    /// Create a new distance constraint.
    pub fn new(i: usize, j: usize, rest_len: f64, compliance: f64) -> Self {
        Self {
            i,
            j,
            rest_len,
            compliance,
        }
    }
}
/// Skinning kernel specification: choose between linear blend skinning (LBS)
/// and dual-quaternion skinning (DQS).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkinningMode {
    /// Linear blend skinning (fast, may have volume loss artefacts).
    LinearBlend,
    /// Dual-quaternion skinning (better arc preservation).
    DualQuaternion,
}
/// Volume preservation constraint for a tetrahedron defined by four particle indices.
#[derive(Debug, Clone)]
pub struct VolumeConstraint {
    /// Particle indices `[a, b, c, d]`.
    pub indices: [usize; 4],
    /// Rest volume of the tetrahedron.
    pub rest_volume: f64,
    /// Compliance (inverse stiffness) for volume preservation.
    pub compliance: f64,
}
impl VolumeConstraint {
    /// Create a new volume preservation constraint.
    pub fn new(indices: [usize; 4], rest_volume: f64, compliance: f64) -> Self {
        Self {
            indices,
            rest_volume,
            compliance,
        }
    }
}
/// Configuration for the shared GPU physics pipeline.
#[derive(Debug, Clone)]
pub struct PhysicsPipelineConfig {
    /// Gravity vector `[gx, gy, gz]` in m/s².
    pub gravity: [f64; 3],
    /// Time step size in seconds.
    pub dt: f64,
    /// Number of substeps per call to `step()`.
    pub substeps: u32,
    /// Enable collision response.
    pub collision_enabled: bool,
    /// Enable FEM forces.
    pub fem_enabled: bool,
    /// Enable XPBD constraints.
    pub xpbd_enabled: bool,
}
/// FEM GPU solver using corotational elements.
///
/// Assembles nodal forces from all elements and applies them to a
/// `DeformableGpuBuffer`.
#[derive(Debug)]
pub struct FEMGpuSolver {
    /// Corotational tetrahedral elements.
    pub elements: Vec<CorotationalElement>,
    /// Total number of nodes.
    pub num_nodes: usize,
}
impl FEMGpuSolver {
    /// Create a new solver for a mesh with `num_nodes` nodes.
    pub fn new(num_nodes: usize) -> Self {
        Self {
            elements: Vec::new(),
            num_nodes,
        }
    }
    /// Add a corotational element.
    pub fn add_element(&mut self, elem: CorotationalElement) {
        self.num_nodes = self
            .num_nodes
            .max(elem.nodes.iter().copied().max().unwrap_or(0) + 1);
        self.elements.push(elem);
    }
    /// Assemble and apply nodal forces to `buf.forces`.
    pub fn assemble_forces(&self, buf: &mut DeformableGpuBuffer) {
        for elem in &self.elements {
            let node_forces = elem.compute_nodal_forces(&buf.positions);
            for (k, &ni) in elem.nodes.iter().enumerate() {
                if ni < buf.len() {
                    buf.forces[ni] = add3(buf.forces[ni], node_forces[k]);
                }
            }
        }
    }
    /// Compute total strain energy (sum over elements).
    pub fn strain_energy(&self, buf: &DeformableGpuBuffer) -> f64 {
        self.elements
            .iter()
            .map(|elem| {
                let [ia, ib, ic, id] = elem.nodes;
                if ia >= buf.len() || ib >= buf.len() || ic >= buf.len() || id >= buf.len() {
                    return 0.0;
                }
                let (_lambda, mu) = elem.lame();
                let pairs = [
                    (buf.positions[ib], buf.positions[ia]),
                    (buf.positions[ic], buf.positions[ia]),
                    (buf.positions[id], buf.positions[ia]),
                ];
                pairs
                    .iter()
                    .zip(elem.dm_inv.iter())
                    .map(|((pb, pa), dm_row)| {
                        let ds = sub3(*pb, *pa);
                        let rest_dir_len =
                            (dm_row[0] * dm_row[0] + dm_row[1] * dm_row[1] + dm_row[2] * dm_row[2])
                                .sqrt();
                        let stretch = length3(ds) * rest_dir_len;
                        0.5 * mu * (stretch - 1.0) * (stretch - 1.0) * elem.rest_volume
                    })
                    .sum::<f64>()
            })
            .sum()
    }
}
/// A detected contact between a particle and a surface.
#[derive(Debug, Clone)]
pub struct ContactInfo {
    /// Particle index.
    pub particle_idx: usize,
    /// Contact normal (pointing away from the surface).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
}
/// GPU collision detection and response for deformable bodies.
///
/// Detects particle-plane, particle-sphere, and particle-face contacts
/// and resolves them with position projection or impulse.
#[derive(Debug)]
pub struct GPUCollisionResponse {
    /// Penalty stiffness for force-based resolution.
    pub penalty_stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
    /// Contact threshold distance.
    pub contact_threshold: f64,
    /// Coefficient of restitution `[0, 1]`.
    pub restitution: f64,
}
impl GPUCollisionResponse {
    /// Create a new collision response handler.
    pub fn new(
        penalty_stiffness: f64,
        damping: f64,
        contact_threshold: f64,
        restitution: f64,
    ) -> Self {
        Self {
            penalty_stiffness,
            damping,
            contact_threshold,
            restitution,
        }
    }
    /// Default collision parameters.
    pub fn default_params() -> Self {
        Self::new(1e5, 100.0, 1e-3, 0.3)
    }
    /// Detect contacts between particles and an infinite half-space plane.
    ///
    /// The half-space is defined by `normal` (pointing outward) and `offset`
    /// (distance from origin along normal).
    pub fn detect_plane_contacts(
        &self,
        positions: &[[f64; 3]],
        normal: [f64; 3],
        offset: f64,
    ) -> Vec<ContactInfo> {
        let n = normalize3(normal);
        positions
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| {
                let d = dot3(p, n) - offset;
                if d < self.contact_threshold {
                    let depth = self.contact_threshold - d;
                    Some(ContactInfo {
                        particle_idx: i,
                        normal: n,
                        depth,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
    /// Detect contacts between particles and a sphere.
    pub fn detect_sphere_contacts(
        &self,
        positions: &[[f64; 3]],
        centre: [f64; 3],
        radius: f64,
    ) -> Vec<ContactInfo> {
        positions
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| {
                let d = sub3(p, centre);
                let dist = length3(d);
                if dist < radius + self.contact_threshold && dist > 1e-15 {
                    let depth = radius - dist + self.contact_threshold;
                    let normal = scale3(d, 1.0 / dist);
                    Some(ContactInfo {
                        particle_idx: i,
                        normal,
                        depth,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
    /// Resolve contacts by projecting positions out of penetration.
    pub fn resolve_positions(&self, buf: &mut DeformableGpuBuffer, contacts: &[ContactInfo]) {
        for c in contacts {
            let i = c.particle_idx;
            if i >= buf.len() || buf.inv_masses[i] < 1e-15 {
                continue;
            }
            buf.positions[i] = add3(buf.positions[i], scale3(c.normal, c.depth));
        }
    }
    /// Resolve contacts by applying impulse-based velocity correction.
    pub fn resolve_velocities(&self, buf: &mut DeformableGpuBuffer, contacts: &[ContactInfo]) {
        for c in contacts {
            let i = c.particle_idx;
            if i >= buf.len() || buf.inv_masses[i] < 1e-15 {
                continue;
            }
            let vn = dot3(buf.velocities[i], c.normal);
            if vn < 0.0 {
                let impulse = -(1.0 + self.restitution) * vn;
                buf.velocities[i] = add3(buf.velocities[i], scale3(c.normal, impulse));
            }
        }
    }
    /// Apply penalty forces for all plane contacts.
    pub fn apply_plane_penalty_forces(
        &self,
        buf: &mut DeformableGpuBuffer,
        normal: [f64; 3],
        offset: f64,
    ) {
        let contacts = self.detect_plane_contacts(&buf.positions, normal, offset);
        for c in &contacts {
            let i = c.particle_idx;
            if i >= buf.len() || buf.inv_masses[i] < 1e-15 {
                continue;
            }
            let vn = dot3(buf.velocities[i], c.normal);
            let f_pen = scale3(c.normal, self.penalty_stiffness * c.depth);
            let f_damp = scale3(c.normal, -self.damping * vn.min(0.0));
            buf.forces[i] = add3(buf.forces[i], add3(f_pen, f_damp));
        }
    }
}
/// Full deformable GPU simulation pipeline: predict → solve → update.
///
/// Manages a `DeformableGpuBuffer` and drives the XPBD solver, FEM forces,
/// and collision response through each time step.
#[derive(Debug)]
pub struct DeformableGpuPipeline {
    /// Particle buffer (SoA layout).
    pub buf: DeformableGpuBuffer,
    /// XPBD constraint solver.
    pub xpbd: XPBDGpuSolver,
    /// FEM elastic force solver.
    pub fem: FEMGpuSolver,
    /// Collision response handler.
    pub collision: GPUCollisionResponse,
    /// Gravity vector `[gx, gy, gz]` m/s².
    pub gravity: [f64; 3],
    /// Time step in seconds.
    pub dt: f64,
    /// Number of XPBD sub-iterations per step.
    pub substeps: u32,
    /// Total elapsed simulation time (s).
    pub time: f64,
    /// Enable FEM forces.
    pub fem_enabled: bool,
    /// Enable XPBD constraints.
    pub xpbd_enabled: bool,
    /// Enable collision response.
    pub collision_enabled: bool,
}
impl DeformableGpuPipeline {
    /// Create a new pipeline with default settings.
    pub fn new(gravity: [f64; 3], dt: f64, substeps: u32) -> Self {
        let num_nodes = 0;
        Self {
            buf: DeformableGpuBuffer::new(),
            xpbd: XPBDGpuSolver::new(10),
            fem: FEMGpuSolver::new(num_nodes),
            collision: GPUCollisionResponse::default_params(),
            gravity,
            dt,
            substeps,
            time: 0.0,
            fem_enabled: true,
            xpbd_enabled: true,
            collision_enabled: true,
        }
    }
    /// Add a particle to the simulation and return its index.
    pub fn add_particle(&mut self, pos: [f64; 3], mass: f64) -> usize {
        self.buf.push(pos, [0.0; 3], mass)
    }
    /// Advance by one macro step (subdivided into `substeps` sub-steps).
    pub fn step(&mut self) {
        let sub_dt = self.dt / self.substeps as f64;
        for _ in 0..self.substeps {
            self.sub_step(sub_dt);
        }
        self.time += self.dt;
    }
    /// Execute one sub-step.
    fn sub_step(&mut self, dt: f64) {
        self.buf.reset_forces();
        self.buf.apply_gravity(self.gravity);
        if self.fem_enabled {
            self.fem.assemble_forces(&mut self.buf);
        }
        self.buf.integrate(dt);
        if self.xpbd_enabled {
            self.xpbd.solve(&mut self.buf, dt);
        }
        if self.collision_enabled {
            let contacts =
                self.collision
                    .detect_plane_contacts(&self.buf.positions, [0.0, 1.0, 0.0], 0.0);
            self.collision.resolve_positions(&mut self.buf, &contacts);
            self.collision.resolve_velocities(&mut self.buf, &contacts);
        }
    }
    /// Return total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.buf.total_kinetic_energy()
    }
    /// Return the centre of mass of all particles.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let n = self.buf.len();
        if n == 0 {
            return [0.0; 3];
        }
        let (total_m, total_pos) = (0..n).fold((0.0f64, [0.0f64; 3]), |(m, acc), i| {
            (
                m + self.buf.masses[i],
                add3(acc, scale3(self.buf.positions[i], self.buf.masses[i])),
            )
        });
        if total_m < 1e-15 {
            [0.0; 3]
        } else {
            scale3(total_pos, 1.0 / total_m)
        }
    }
}
/// GPU deformable mesh: vertex buffer, face buffer, skinning weights, and
/// blend shapes.
#[derive(Debug)]
pub struct DeformableGpuMesh {
    /// Vertex buffer.
    pub vertices: Vec<DeformableVertex>,
    /// Face (triangle) index buffer.
    pub faces: Vec<DeformableFace>,
    /// Blend shape (morph target) list.
    pub blend_shapes: Vec<BlendShape>,
}
impl DeformableGpuMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            faces: Vec::new(),
            blend_shapes: Vec::new(),
        }
    }
    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, v: DeformableVertex) -> u32 {
        let idx = self.vertices.len() as u32;
        self.vertices.push(v);
        idx
    }
    /// Add a triangular face.
    pub fn add_face(&mut self, a: u32, b: u32, c: u32) {
        self.faces.push(DeformableFace::new(a, b, c));
    }
    /// Add a blend shape and return its index.
    pub fn add_blend_shape(&mut self, bs: BlendShape) -> usize {
        let idx = self.blend_shapes.len();
        self.blend_shapes.push(bs);
        idx
    }
    /// Apply all active blend shapes to vertex positions (additive).
    pub fn apply_blend_shapes(&mut self) {
        for v in &mut self.vertices {
            v.position = v.rest_position;
        }
        for bs in &self.blend_shapes {
            if bs.weight.abs() < 1e-15 {
                continue;
            }
            for (i, v) in self.vertices.iter_mut().enumerate() {
                if i < bs.offsets.len() {
                    let off = scale3(bs.offsets[i], bs.weight);
                    v.position = add3(v.position, off);
                }
            }
        }
    }
    /// Compute axis-aligned bounding box of current vertex positions.
    ///
    /// Returns `([min_x, min_y, min_z], [max_x, max_y, max_z])`.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for k in 0..3 {
                lo[k] = lo[k].min(v.position[k]);
                hi[k] = hi[k].max(v.position[k]);
            }
        }
        (lo, hi)
    }
    /// Compute face normals for all triangles.
    pub fn compute_face_normals(&self) -> Vec<[f64; 3]> {
        self.faces
            .iter()
            .map(|f| {
                let a = self.vertices[f.a as usize].position;
                let b = self.vertices[f.b as usize].position;
                let c = self.vertices[f.c as usize].position;
                normalize3(cross3(sub3(b, a), sub3(c, a)))
            })
            .collect()
    }
}
/// A linear tetrahedral FEM element.
#[derive(Debug, Clone)]
pub struct TetElement {
    /// Vertex indices `[a, b, c, d]`.
    pub nodes: [usize; 4],
    /// Rest-pose volume (m³).
    pub rest_volume: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl TetElement {
    /// Create a new tetrahedral element.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nodes: [usize; 4],
        rest_volume: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            nodes,
            rest_volume,
            youngs_modulus,
            poisson_ratio,
        }
    }
    /// Compute Lamé parameters from Young's modulus and Poisson's ratio.
    pub fn lame_params(&self) -> (f64, f64) {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        (lambda, mu)
    }
}
/// A triangular face in the deformable mesh.
#[derive(Debug, Clone, Copy)]
pub struct DeformableFace {
    /// Vertex index A.
    pub a: u32,
    /// Vertex index B.
    pub b: u32,
    /// Vertex index C.
    pub c: u32,
}
impl DeformableFace {
    /// Construct a new face from three vertex indices.
    pub fn new(a: u32, b: u32, c: u32) -> Self {
        Self { a, b, c }
    }
}
/// A joint transform for skinning, represented as a 4×3 matrix (row-major).
///
/// The matrix encodes the joint-to-world affine transform: `M * v`.
#[derive(Debug, Clone, Copy)]
pub struct JointTransform {
    /// Row 0 of the 3×3 rotation part, plus translation x in column 3.
    pub row0: [f64; 4],
    /// Row 1 of the 3×3 rotation part, plus translation y in column 3.
    pub row1: [f64; 4],
    /// Row 2 of the 3×3 rotation part, plus translation z in column 3.
    pub row2: [f64; 4],
}
impl JointTransform {
    /// Identity joint transform.
    pub fn identity() -> Self {
        Self {
            row0: [1.0, 0.0, 0.0, 0.0],
            row1: [0.0, 1.0, 0.0, 0.0],
            row2: [0.0, 0.0, 1.0, 0.0],
        }
    }
    /// Apply this transform to a point (homogeneous w=1).
    pub fn transform_point(&self, p: [f64; 3]) -> [f64; 3] {
        [
            self.row0[0] * p[0] + self.row0[1] * p[1] + self.row0[2] * p[2] + self.row0[3],
            self.row1[0] * p[0] + self.row1[1] * p[1] + self.row1[2] * p[2] + self.row1[3],
            self.row2[0] * p[0] + self.row2[1] * p[1] + self.row2[2] * p[2] + self.row2[3],
        ]
    }
    /// Blend (add) another transform scaled by `w` into this one.
    #[allow(dead_code)]
    fn add_scaled(&self, other: &JointTransform, w: f64) -> JointTransform {
        let blend_row = |r: [f64; 4], o: [f64; 4]| {
            [
                r[0] + w * o[0],
                r[1] + w * o[1],
                r[2] + w * o[2],
                r[3] + w * o[3],
            ]
        };
        JointTransform {
            row0: blend_row(self.row0, other.row0),
            row1: blend_row(self.row1, other.row1),
            row2: blend_row(self.row2, other.row2),
        }
    }
}
/// FEM GPU kernel (CPU mock): element stiffness assembly and force computation.
#[derive(Debug)]
pub struct FEMGpuKernel {
    /// Tetrahedral elements.
    pub elements: Vec<TetElement>,
    /// Number of nodes.
    pub num_nodes: usize,
}
impl FEMGpuKernel {
    /// Create a new FEM kernel for a mesh with `num_nodes` nodes.
    pub fn new(num_nodes: usize) -> Self {
        Self {
            elements: Vec::new(),
            num_nodes,
        }
    }
    /// Add a tetrahedral element.
    pub fn add_element(&mut self, elem: TetElement) {
        self.num_nodes = self
            .num_nodes
            .max(elem.nodes.iter().copied().max().unwrap_or(0) + 1);
        self.elements.push(elem);
    }
    /// Compute internal elastic forces for all elements.
    ///
    /// `positions` – current node positions `[x, y, z]`.
    /// `rest_positions` – rest-pose node positions.
    ///
    /// Returns a per-node force vector.
    pub fn compute_forces(
        &self,
        positions: &[[f64; 3]],
        rest_positions: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0f64; 3]; self.num_nodes];
        for elem in &self.elements {
            let [ia, ib, ic, id] = elem.nodes;
            if ia >= positions.len()
                || ib >= positions.len()
                || ic >= positions.len()
                || id >= positions.len()
            {
                continue;
            }
            let pa = positions[ia];
            let pb = positions[ib];
            let pc = positions[ic];
            let pd = positions[id];
            let ra = rest_positions[ia];
            let rb = rest_positions[ib];
            let rc = rest_positions[ic];
            let rd = rest_positions[id];
            let e1 = sub3(pb, pa);
            let e2 = sub3(pc, pa);
            let e3 = sub3(pd, pa);
            let r1 = sub3(rb, ra);
            let r2 = sub3(rc, ra);
            let r3 = sub3(rd, ra);
            let (_lambda, mu) = elem.lame_params();
            let stretch_force = |e_cur: [f64; 3], e_rest: [f64; 3]| -> [f64; 3] {
                let len_cur = length3(e_cur);
                let len_rest = length3(e_rest);
                if len_rest < 1e-15 {
                    return [0.0; 3];
                }
                let stretch = (len_cur - len_rest) / len_rest;
                let dir = normalize3(e_cur);
                scale3(dir, -mu * stretch * elem.rest_volume)
            };
            let f1 = stretch_force(e1, r1);
            let f2 = stretch_force(e2, r2);
            let f3 = stretch_force(e3, r3);
            let fa = [
                -(f1[0] + f2[0] + f3[0]),
                -(f1[1] + f2[1] + f3[1]),
                -(f1[2] + f2[2] + f3[2]),
            ];
            forces[ia] = add3(forces[ia], fa);
            forces[ib] = add3(forces[ib], f1);
            forces[ic] = add3(forces[ic], f2);
            forces[id] = add3(forces[id], f3);
        }
        forces
    }
    /// Compute the total elastic potential energy of the mesh.
    pub fn compute_strain_energy(
        &self,
        positions: &[[f64; 3]],
        rest_positions: &[[f64; 3]],
    ) -> f64 {
        let mut energy = 0.0;
        for elem in &self.elements {
            let [ia, ib, ic, id] = elem.nodes;
            if ia >= positions.len()
                || ib >= positions.len()
                || ic >= positions.len()
                || id >= positions.len()
            {
                continue;
            }
            let (_lambda, mu) = elem.lame_params();
            for (cur, rest) in [
                (positions[ib], rest_positions[ib]),
                (positions[ic], rest_positions[ic]),
                (positions[id], rest_positions[id]),
            ] {
                let e_cur = sub3(cur, positions[ia]);
                let e_rest = sub3(rest, rest_positions[ia]);
                let len_cur = length3(e_cur);
                let len_rest = length3(e_rest);
                if len_rest > 1e-15 {
                    let stretch = (len_cur - len_rest) / len_rest;
                    energy += 0.5 * mu * stretch * stretch * elem.rest_volume;
                }
            }
        }
        energy
    }
}
/// Parameters for GPU collision response.
#[derive(Debug, Clone)]
pub struct CollisionResponseParams {
    /// Penalty stiffness coefficient (N/m).
    pub penalty_stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Contact threshold distance (m).
    pub contact_threshold: f64,
}
impl CollisionResponseParams {
    /// Create default collision parameters.
    pub fn new(penalty_stiffness: f64, damping: f64, contact_threshold: f64) -> Self {
        Self {
            penalty_stiffness,
            damping,
            contact_threshold,
        }
    }
}
/// A blend shape (morph target) applied to the mesh.
#[derive(Debug, Clone)]
pub struct BlendShape {
    /// Human-readable name of this blend shape (e.g. `"smile"`).
    pub name: String,
    /// Per-vertex displacement offsets `[dx, dy, dz]`.
    pub offsets: Vec<[f64; 3]>,
    /// Current blend weight in `[0, 1]`.
    pub weight: f64,
}
impl BlendShape {
    /// Create a new blend shape with given offsets and weight `0`.
    pub fn new(name: impl Into<String>, offsets: Vec<[f64; 3]>) -> Self {
        Self {
            name: name.into(),
            offsets,
            weight: 0.0,
        }
    }
}
/// Unified physics pipeline for softbody, rigid, and fluid simulations on
/// shared GPU resources (CPU mock).
///
/// Manages a deformable mesh alongside XPBD and FEM kernels, and runs the
/// full simulation loop per time step.
#[derive(Debug)]
pub struct PhysicsGpuPipeline {
    /// Pipeline configuration.
    pub config: PhysicsPipelineConfig,
    /// The deformable mesh being simulated.
    pub mesh: DeformableGpuMesh,
    /// XPBD constraint solver.
    pub xpbd: XPBDGpuKernel,
    /// FEM elastic force kernel.
    pub fem: FEMGpuKernel,
    /// Collision response kernel.
    pub collision: CollisionResponseGpu,
    /// Per-vertex inverse masses.
    pub inv_masses: Vec<f64>,
    /// Per-vertex velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Rest-pose positions for FEM.
    pub rest_positions: Vec<[f64; 3]>,
    /// Total simulation time elapsed (s).
    pub time: f64,
}
impl PhysicsGpuPipeline {
    /// Create a new pipeline with the given configuration and mesh.
    pub fn new(config: PhysicsPipelineConfig, mesh: DeformableGpuMesh) -> Self {
        let n = mesh.vertices.len();
        let inv_masses: Vec<f64> = mesh
            .vertices
            .iter()
            .map(|v| if v.mass < 1e-15 { 0.0 } else { 1.0 / v.mass })
            .collect();
        let velocities: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| v.velocity).collect();
        let rest_positions: Vec<[f64; 3]> = mesh.vertices.iter().map(|v| v.rest_position).collect();
        Self {
            xpbd: XPBDGpuKernel::new(8),
            fem: FEMGpuKernel::new(n),
            collision: CollisionResponseGpu::new(CollisionResponseParams::default()),
            config,
            mesh,
            inv_masses,
            velocities,
            rest_positions,
            time: 0.0,
        }
    }
    /// Advance the simulation by one macro time step (subdivided into
    /// `config.substeps` sub-steps).
    pub fn step(&mut self) {
        let sub_dt = self.config.dt / self.config.substeps as f64;
        for _ in 0..self.config.substeps {
            self.substep(sub_dt);
        }
        self.time += self.config.dt;
    }
    /// Execute a single sub-step.
    fn substep(&mut self, dt: f64) {
        let n = self.mesh.vertices.len();
        let mut positions: Vec<[f64; 3]> = self.mesh.vertices.iter().map(|v| v.position).collect();
        if self.config.xpbd_enabled {
            self.xpbd.integrate(
                &mut positions,
                &mut self.velocities,
                &self.inv_masses,
                self.config.gravity,
                dt,
            );
        }
        if self.config.fem_enabled && !self.rest_positions.is_empty() {
            let fem_forces = self.fem.compute_forces(&positions, &self.rest_positions);
            for i in 0..n.min(fem_forces.len()) {
                if self.inv_masses[i] < 1e-15 {
                    continue;
                }
                let a = scale3(fem_forces[i], self.inv_masses[i]);
                self.velocities[i] = add3(self.velocities[i], scale3(a, dt));
                positions[i] = add3(positions[i], scale3(self.velocities[i], dt));
            }
        }
        if self.config.xpbd_enabled {
            self.xpbd.solve_step(&mut positions, &self.inv_masses, dt);
        }
        for (i, v) in self.mesh.vertices.iter_mut().enumerate() {
            if i < positions.len() {
                if self.inv_masses[i] > 1e-15 {
                    let dp = sub3(positions[i], v.position);
                    v.velocity = scale3(dp, 1.0 / dt);
                }
                v.position = positions[i];
            }
        }
        for (i, v) in self.mesh.vertices.iter().enumerate() {
            if i < self.velocities.len() {
                self.velocities[i] = v.velocity;
            }
        }
    }
    /// Return the current kinetic energy of the mesh.
    pub fn kinetic_energy(&self) -> f64 {
        self.mesh
            .vertices
            .iter()
            .map(|v| {
                let v2 = len2_3(v.velocity);
                0.5 * v.mass * v2
            })
            .sum()
    }
    /// Return the current centre of mass of the mesh.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let (total_m, total_pos) = self
            .mesh
            .vertices
            .iter()
            .fold((0.0f64, [0.0f64; 3]), |(m, acc), v| {
                (m + v.mass, add3(acc, scale3(v.position, v.mass)))
            });
        if total_m < 1e-15 {
            [0.0; 3]
        } else {
            scale3(total_pos, 1.0 / total_m)
        }
    }
}
/// A dual quaternion for dual-quaternion skinning.
///
/// Stored as `(q_real, q_dual)` where both are `[w, x, y, z]` quaternions.
#[derive(Debug, Clone, Copy)]
pub struct DualQuat {
    /// Real part `[w, x, y, z]`.
    pub real: [f64; 4],
    /// Dual part `[w, x, y, z]`.
    pub dual: [f64; 4],
}
impl DualQuat {
    /// Construct identity dual quaternion (no rotation, no translation).
    pub fn identity() -> Self {
        Self {
            real: [1.0, 0.0, 0.0, 0.0],
            dual: [0.0, 0.0, 0.0, 0.0],
        }
    }
    /// Construct a dual quaternion from a pure translation `[tx, ty, tz]`.
    pub fn from_translation(t: [f64; 3]) -> Self {
        Self {
            real: [1.0, 0.0, 0.0, 0.0],
            dual: [0.0, 0.5 * t[0], 0.5 * t[1], 0.5 * t[2]],
        }
    }
    /// Add another dual quaternion scaled by `w` (for DQS blending).
    pub fn add_scaled(&self, other: &DualQuat, w: f64) -> DualQuat {
        DualQuat {
            real: [
                self.real[0] + w * other.real[0],
                self.real[1] + w * other.real[1],
                self.real[2] + w * other.real[2],
                self.real[3] + w * other.real[3],
            ],
            dual: [
                self.dual[0] + w * other.dual[0],
                self.dual[1] + w * other.dual[1],
                self.dual[2] + w * other.dual[2],
                self.dual[3] + w * other.dual[3],
            ],
        }
    }
    /// Normalise so that `||real|| = 1`.
    pub fn normalize(&self) -> DualQuat {
        let len = (self.real[0] * self.real[0]
            + self.real[1] * self.real[1]
            + self.real[2] * self.real[2]
            + self.real[3] * self.real[3])
            .sqrt();
        if len < 1e-15 {
            return *self;
        }
        DualQuat {
            real: [
                self.real[0] / len,
                self.real[1] / len,
                self.real[2] / len,
                self.real[3] / len,
            ],
            dual: [
                self.dual[0] / len,
                self.dual[1] / len,
                self.dual[2] / len,
                self.dual[3] / len,
            ],
        }
    }
    /// Transform a point using this (normalised) dual quaternion.
    pub fn transform_point(&self, p: [f64; 3]) -> [f64; 3] {
        let dq = self.normalize();
        let qr = dq.real;
        let qd = dq.dual;
        let tx = 2.0 * (qd[0] * (-qr[1]) + qd[1] * qr[0] + qd[2] * (-qr[3]) - qd[3] * (-qr[2]));
        let ty = 2.0 * (qd[0] * (-qr[2]) - qd[1] * (-qr[3]) + qd[2] * qr[0] + qd[3] * (-qr[1]));
        let tz = 2.0 * (qd[0] * (-qr[3]) + qd[1] * (-qr[2]) - qd[2] * (-qr[1]) + qd[3] * qr[0]);
        let rot_p = quat_rotate(qr, p);
        add3(rot_p, [tx, ty, tz])
    }
}
/// A single particle node in the deformable GPU simulation (SoA layout entry).
///
/// Stores position, velocity, mass, and accumulated force for one particle.
#[derive(Debug, Clone)]
pub struct DeformableGpuNode {
    /// Current position `[x, y, z]`.
    pub pos: [f64; 3],
    /// Current velocity `[vx, vy, vz]`.
    pub vel: [f64; 3],
    /// Particle mass in kg.
    pub mass: f64,
    /// Accumulated external force `[fx, fy, fz]` (reset each step).
    pub force: [f64; 3],
}
impl DeformableGpuNode {
    /// Create a new particle at `pos` with the given mass and zero velocity/force.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            force: [0.0; 3],
        }
    }
    /// Add an external force to this node.
    pub fn apply_force(&mut self, f: [f64; 3]) {
        self.force = add3(self.force, f);
    }
    /// Reset the accumulated force to zero.
    pub fn reset_force(&mut self) {
        self.force = [0.0; 3];
    }
    /// Compute the kinetic energy of this node (`0.5 m v²`).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
}
/// Corotational FEM element: stores rest-pose shape matrix inverse.
#[derive(Debug, Clone)]
pub struct CorotationalElement {
    /// Node indices `[a, b, c, d]`.
    pub nodes: [usize; 4],
    /// Rest-pose shape matrix inverse (3×3, row-major).
    pub dm_inv: [[f64; 3]; 3],
    /// Rest volume (m³).
    pub rest_volume: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl CorotationalElement {
    /// Construct a corotational element from four node positions.
    ///
    /// `rest` – rest-pose positions for nodes `[a, b, c, d]`.
    pub fn from_rest(nodes: [usize; 4], rest: &[[f64; 3]; 4], e: f64, nu: f64) -> Self {
        let ra = rest[0];
        let rb = rest[1];
        let rc = rest[2];
        let rd = rest[3];
        let c0 = sub3(rb, ra);
        let c1 = sub3(rc, ra);
        let c2 = sub3(rd, ra);
        let det = c0[0] * (c1[1] * c2[2] - c1[2] * c2[1]) - c0[1] * (c1[0] * c2[2] - c1[2] * c2[0])
            + c0[2] * (c1[0] * c2[1] - c1[1] * c2[0]);
        let rest_volume = det.abs() / 6.0;
        let dm_inv = if det.abs() > 1e-15 {
            let inv_det = 1.0 / det;
            [
                [
                    (c1[1] * c2[2] - c1[2] * c2[1]) * inv_det,
                    (c0[2] * c2[1] - c0[1] * c2[2]) * inv_det,
                    (c0[1] * c1[2] - c0[2] * c1[1]) * inv_det,
                ],
                [
                    (c1[2] * c2[0] - c1[0] * c2[2]) * inv_det,
                    (c0[0] * c2[2] - c0[2] * c2[0]) * inv_det,
                    (c0[2] * c1[0] - c0[0] * c1[2]) * inv_det,
                ],
                [
                    (c1[0] * c2[1] - c1[1] * c2[0]) * inv_det,
                    (c0[1] * c2[0] - c0[0] * c2[1]) * inv_det,
                    (c0[0] * c1[1] - c0[1] * c1[0]) * inv_det,
                ],
            ]
        } else {
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        };
        Self {
            nodes,
            dm_inv,
            rest_volume,
            youngs_modulus: e,
            poisson_ratio: nu,
        }
    }
    /// Compute Lamé parameters.
    pub fn lame(&self) -> (f64, f64) {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        (lambda, mu)
    }
    /// Compute nodal forces for this element given current positions.
    ///
    /// Returns forces for nodes `[a, b, c, d]` in order.
    pub fn compute_nodal_forces(&self, positions: &[[f64; 3]]) -> [[f64; 3]; 4] {
        let [ia, ib, ic, id] = self.nodes;
        if ia >= positions.len()
            || ib >= positions.len()
            || ic >= positions.len()
            || id >= positions.len()
        {
            return [[0.0; 3]; 4];
        }
        let pa = positions[ia];
        let pb = positions[ib];
        let pc = positions[ic];
        let pd = positions[id];
        let ds0 = sub3(pb, pa);
        let ds1 = sub3(pc, pa);
        let ds2 = sub3(pd, pa);
        let (_lambda, mu) = self.lame();
        let mut fa = [0.0f64; 3];
        let mut fb = [0.0f64; 3];
        let mut fc = [0.0f64; 3];
        let mut fd = [0.0f64; 3];
        let f_col = |ds: [f64; 3], dm_row: [f64; 3]| -> [f64; 3] {
            [
                ds[0] * dm_row[0] + ds[1] * dm_row[1] + ds[2] * dm_row[2],
                ds0[0] * dm_row[0] + ds0[1] * dm_row[1] + ds0[2] * dm_row[2],
                ds1[0] * dm_row[0] + ds1[1] * dm_row[1] + ds1[2] * dm_row[2],
            ]
        };
        let _ = f_col;
        for (ds, dm_row, fp) in [
            (ds0, self.dm_inv[0], &mut fb),
            (ds1, self.dm_inv[1], &mut fc),
            (ds2, self.dm_inv[2], &mut fd),
        ] {
            let rest_dir_len =
                (dm_row[0] * dm_row[0] + dm_row[1] * dm_row[1] + dm_row[2] * dm_row[2]).sqrt();
            let stretch = length3(ds) * rest_dir_len;
            let dir = normalize3(ds);
            let f_elastic = scale3(dir, -mu * (stretch - 1.0) * self.rest_volume);
            *fp = add3(*fp, f_elastic);
            fa = add3(fa, scale3(f_elastic, -1.0));
        }
        [fa, fb, fc, fd]
    }
}
/// Shape-matching constraint: pulls a group of particles toward a reference shape.
#[derive(Debug, Clone)]
pub struct ShapeMatchConstraint {
    /// Particle indices belonging to this shape.
    pub indices: Vec<usize>,
    /// Rest-pose positions for each particle in `indices`.
    pub rest_positions: Vec<[f64; 3]>,
    /// Stiffness factor `[0, 1]`.
    pub stiffness: f64,
}
impl ShapeMatchConstraint {
    /// Create a new shape-matching constraint.
    pub fn new(indices: Vec<usize>, rest_positions: Vec<[f64; 3]>, stiffness: f64) -> Self {
        Self {
            indices,
            rest_positions,
            stiffness,
        }
    }
}
/// XPBD GPU solver with distance, volume, and shape-matching constraints.
///
/// This is a higher-level wrapper around `XPBDGpuKernel` that also handles
/// volume preservation and shape matching.
#[derive(Debug)]
pub struct XPBDGpuSolver {
    /// Distance constraints.
    pub distance_constraints: Vec<XpbdDistConstraint>,
    /// Volume preservation constraints.
    pub volume_constraints: Vec<VolumeConstraint>,
    /// Shape matching constraints.
    pub shape_constraints: Vec<ShapeMatchConstraint>,
    /// Number of solver iterations per sub-step.
    pub iterations: u32,
}
impl XPBDGpuSolver {
    /// Create a new solver with the given number of iterations.
    pub fn new(iterations: u32) -> Self {
        Self {
            distance_constraints: Vec::new(),
            volume_constraints: Vec::new(),
            shape_constraints: Vec::new(),
            iterations,
        }
    }
    /// Add a distance constraint between particles `i` and `j`.
    pub fn add_distance(&mut self, i: usize, j: usize, rest_len: f64, compliance: f64) {
        self.distance_constraints
            .push(XpbdDistConstraint::new(i, j, rest_len, compliance));
    }
    /// Add a volume preservation constraint.
    pub fn add_volume(&mut self, vc: VolumeConstraint) {
        self.volume_constraints.push(vc);
    }
    /// Add a shape-matching constraint.
    pub fn add_shape(&mut self, sc: ShapeMatchConstraint) {
        self.shape_constraints.push(sc);
    }
    /// Solve one XPBD step on a `DeformableGpuBuffer`.
    pub fn solve(&self, buf: &mut DeformableGpuBuffer, dt: f64) {
        let alpha = 1.0 / (dt * dt);
        for _ in 0..self.iterations {
            for c in &self.distance_constraints {
                if c.i >= buf.len() || c.j >= buf.len() {
                    continue;
                }
                let pi = buf.positions[c.i];
                let pj = buf.positions[c.j];
                let d = sub3(pi, pj);
                let dist = length3(d);
                if dist < 1e-15 {
                    continue;
                }
                let constraint_val = dist - c.rest_len;
                let wi = buf.inv_masses[c.i];
                let wj = buf.inv_masses[c.j];
                let denom = wi + wj + c.compliance * alpha;
                if denom < 1e-15 {
                    continue;
                }
                let lambda = -constraint_val / denom;
                let grad = normalize3(d);
                buf.positions[c.i] = add3(pi, scale3(grad, lambda * wi));
                buf.positions[c.j] = add3(pj, scale3(grad, -lambda * wj));
            }
            for vc in &self.volume_constraints {
                let [ia, ib, ic, id] = vc.indices;
                if ia >= buf.len() || ib >= buf.len() || ic >= buf.len() || id >= buf.len() {
                    continue;
                }
                let pa = buf.positions[ia];
                let pb = buf.positions[ib];
                let pc = buf.positions[ic];
                let pd = buf.positions[id];
                let ab = sub3(pb, pa);
                let ac = sub3(pc, pa);
                let ad = sub3(pd, pa);
                let cur_vol = dot3(ab, cross3(ac, ad)).abs() / 6.0;
                let vol_err = cur_vol - vc.rest_volume;
                if vol_err.abs() < 1e-15 {
                    continue;
                }
                let centroid = scale3(add3(add3(pa, pb), add3(pc, pd)), 0.25);
                let scale_factor = if cur_vol > 1e-15 {
                    (vc.rest_volume / cur_vol).cbrt()
                } else {
                    1.0_f64
                };
                let correction_weight = (1.0 - scale_factor) * vc.compliance.min(1.0);
                for &idx in &[ia, ib, ic, id] {
                    if buf.inv_masses[idx] < 1e-15 {
                        continue;
                    }
                    let to_centroid = sub3(centroid, buf.positions[idx]);
                    buf.positions[idx] =
                        add3(buf.positions[idx], scale3(to_centroid, correction_weight));
                }
            }
            for sc in &self.shape_constraints {
                if sc.indices.is_empty() {
                    continue;
                }
                let n = sc.indices.len() as f64;
                let mut cur_centroid = [0.0f64; 3];
                for &idx in &sc.indices {
                    if idx < buf.len() {
                        cur_centroid = add3(cur_centroid, buf.positions[idx]);
                    }
                }
                cur_centroid = scale3(cur_centroid, 1.0 / n);
                let mut rest_centroid = [0.0f64; 3];
                for &rp in &sc.rest_positions {
                    rest_centroid = add3(rest_centroid, rp);
                }
                rest_centroid = scale3(rest_centroid, 1.0 / n);
                for (k, &idx) in sc.indices.iter().enumerate() {
                    if idx >= buf.len() || buf.inv_masses[idx] < 1e-15 {
                        continue;
                    }
                    let rest_offset = if k < sc.rest_positions.len() {
                        sub3(sc.rest_positions[k], rest_centroid)
                    } else {
                        [0.0; 3]
                    };
                    let target = add3(cur_centroid, rest_offset);
                    let delta = sub3(target, buf.positions[idx]);
                    buf.positions[idx] = add3(buf.positions[idx], scale3(delta, sc.stiffness));
                }
            }
        }
    }
}
/// A single vertex in the deformable mesh, including skinning metadata.
#[derive(Debug, Clone)]
pub struct DeformableVertex {
    /// Rest-pose position `[x, y, z]`.
    pub rest_position: [f64; 3],
    /// Current world-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// Current velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Vertex mass in kg.
    pub mass: f64,
    /// Skinning joint indices (up to 4 per vertex).
    pub joint_indices: [u32; 4],
    /// Skinning weights corresponding to `joint_indices` (must sum ≤ 1).
    pub joint_weights: [f64; 4],
}
impl DeformableVertex {
    /// Construct a new vertex at the given rest position with unit mass.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            rest_position: position,
            position,
            velocity: [0.0; 3],
            mass: 1.0,
            joint_indices: [0; 4],
            joint_weights: [1.0, 0.0, 0.0, 0.0],
        }
    }
}
/// GPU-side SoA (Structure-of-Arrays) buffer for deformable body particles.
///
/// Stores positions, velocities, masses, and forces as separate arrays to
/// enable cache-friendly access patterns in kernel dispatches.
#[derive(Debug, Default)]
pub struct DeformableGpuBuffer {
    /// Particle positions `[[x, y, z\], ...]`.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities `[[vx, vy, vz\], ...]`.
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses `[m0, m1, ...]`.
    pub masses: Vec<f64>,
    /// Accumulated forces `[[fx, fy, fz\], ...]`.
    pub forces: Vec<[f64; 3]>,
    /// Inverse masses `[1/m0, 1/m1, ...]` (0 for fixed particles).
    pub inv_masses: Vec<f64>,
}
impl DeformableGpuBuffer {
    /// Create an empty buffer.
    pub fn new() -> Self {
        Self::default()
    }
    /// Build a buffer from a list of `DeformableGpuNode` particles.
    pub fn from_nodes(nodes: &[DeformableGpuNode]) -> Self {
        let positions = nodes.iter().map(|n| n.pos).collect();
        let velocities = nodes.iter().map(|n| n.vel).collect();
        let masses: Vec<f64> = nodes.iter().map(|n| n.mass).collect();
        let forces = nodes.iter().map(|n| n.force).collect();
        let inv_masses = masses
            .iter()
            .map(|&m| if m > 1e-15 { 1.0 / m } else { 0.0 })
            .collect();
        Self {
            positions,
            velocities,
            masses,
            forces,
            inv_masses,
        }
    }
    /// Number of particles in the buffer.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Returns `true` if the buffer holds no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Add a particle and return its index.
    pub fn push(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64) -> usize {
        let idx = self.positions.len();
        self.positions.push(pos);
        self.velocities.push(vel);
        self.masses.push(mass);
        self.forces.push([0.0; 3]);
        self.inv_masses
            .push(if mass > 1e-15 { 1.0 / mass } else { 0.0 });
        idx
    }
    /// Reset all accumulated forces to zero.
    pub fn reset_forces(&mut self) {
        for f in &mut self.forces {
            *f = [0.0; 3];
        }
    }
    /// Apply gravity to all free particles (adds `m * g` to force).
    pub fn apply_gravity(&mut self, gravity: [f64; 3]) {
        for i in 0..self.len() {
            if self.inv_masses[i] > 1e-15 {
                let fg = scale3(gravity, self.masses[i]);
                self.forces[i] = add3(self.forces[i], fg);
            }
        }
    }
    /// Semi-implicit Euler integration over all free particles.
    pub fn integrate(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.inv_masses[i] < 1e-15 {
                continue;
            }
            let a = scale3(self.forces[i], self.inv_masses[i]);
            self.velocities[i] = add3(self.velocities[i], scale3(a, dt));
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }
    /// Compute total kinetic energy of the buffer.
    pub fn total_kinetic_energy(&self) -> f64 {
        (0..self.len())
            .map(|i| 0.5 * self.masses[i] * dot3(self.velocities[i], self.velocities[i]))
            .sum()
    }
    /// Copy positions into a flat `Vec`f64` (x0,y0,z0, x1,y1,z1, ...).
    pub fn flatten_positions(&self) -> Vec<f64> {
        self.positions.iter().flat_map(|&p| p).collect()
    }
}
/// GPU collision response kernel (CPU mock).
///
/// Detects vertex-face penetration and applies penalty forces.
#[derive(Debug)]
pub struct CollisionResponseGpu {
    /// Collision parameters.
    pub params: CollisionResponseParams,
}
impl CollisionResponseGpu {
    /// Create a new collision response kernel with the given parameters.
    pub fn new(params: CollisionResponseParams) -> Self {
        Self { params }
    }
    /// Compute the signed distance from a point `p` to the plane of triangle
    /// `(v0, v1, v2)`.
    ///
    /// Negative values mean `p` is behind the triangle.
    pub fn signed_distance_to_triangle(
        &self,
        p: [f64; 3],
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
    ) -> f64 {
        let n = normalize3(cross3(sub3(v2, v0), sub3(v1, v0)));
        dot3(sub3(p, v0), n)
    }
    /// Apply penalty forces for vertex-face collisions.
    ///
    /// For each vertex in `query_positions`, checks against faces defined by
    /// `face_vertices`. Adds penalty forces to `forces`.
    pub fn apply_penalty_forces(
        &self,
        query_positions: &[[f64; 3]],
        query_velocities: &[[f64; 3]],
        face_vertices: &[([f64; 3], [f64; 3], [f64; 3])],
        forces: &mut [[f64; 3]],
    ) {
        for (vi, &p) in query_positions.iter().enumerate() {
            if vi >= forces.len() {
                break;
            }
            let vel = if vi < query_velocities.len() {
                query_velocities[vi]
            } else {
                [0.0; 3]
            };
            for &(v0, v1, v2) in face_vertices {
                let sd = self.signed_distance_to_triangle(p, v0, v1, v2);
                if sd > 0.0 || sd < -self.params.contact_threshold {
                    continue;
                }
                let n = normalize3(cross3(sub3(v2, v0), sub3(v1, v0)));
                let pen = -sd;
                let f_pen = scale3(n, self.params.penalty_stiffness * pen);
                let vel_n = dot3(vel, n);
                let f_damp = scale3(n, -self.params.damping * vel_n.min(0.0));
                forces[vi] = add3(forces[vi], add3(f_pen, f_damp));
            }
        }
    }
    /// Apply sphere-vertex collision penalty forces.
    ///
    /// Pushes vertices outside a sphere of radius `r` centred at `centre`.
    pub fn apply_sphere_penalty(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        centre: [f64; 3],
        radius: f64,
        forces: &mut [[f64; 3]],
    ) {
        for (i, &p) in positions.iter().enumerate() {
            if i >= forces.len() {
                break;
            }
            let d = sub3(p, centre);
            let dist = length3(d);
            if dist >= radius || dist < 1e-15 {
                continue;
            }
            let pen = radius - dist;
            let n = scale3(d, 1.0 / dist);
            let f_pen = scale3(n, self.params.penalty_stiffness * pen);
            let vel = if i < velocities.len() {
                velocities[i]
            } else {
                [0.0; 3]
            };
            let vel_n = dot3(vel, n);
            let f_damp = scale3(n, -self.params.damping * vel_n.min(0.0));
            forces[i] = add3(forces[i], add3(f_pen, f_damp));
        }
    }
}
/// XPBD GPU kernel (CPU mock): iterative position-based constraint solve.
#[derive(Debug)]
pub struct XPBDGpuKernel {
    /// List of distance constraints.
    pub constraints: Vec<XpbdDistConstraint>,
    /// Number of solver iterations per time step.
    pub iterations: u32,
}
impl XPBDGpuKernel {
    /// Create a new XPBD kernel with the given number of iterations.
    pub fn new(iterations: u32) -> Self {
        Self {
            constraints: Vec::new(),
            iterations,
        }
    }
    /// Add a distance constraint.
    pub fn add_constraint(&mut self, c: XpbdDistConstraint) {
        self.constraints.push(c);
    }
    /// Solve one XPBD step.
    ///
    /// * `positions` – mutable positions `\[x, y, z\]` for each particle.
    /// * `inv_masses` – inverse mass for each particle (0 = fixed).
    /// * `dt` – time step in seconds.
    pub fn solve_step(&self, positions: &mut [[f64; 3]], inv_masses: &[f64], dt: f64) {
        let alpha = 1.0 / (dt * dt);
        for _ in 0..self.iterations {
            for c in &self.constraints {
                if c.i >= positions.len() || c.j >= positions.len() {
                    continue;
                }
                let pi = positions[c.i];
                let pj = positions[c.j];
                let d = sub3(pi, pj);
                let dist = length3(d);
                if dist < 1e-15 {
                    continue;
                }
                let constraint_val = dist - c.rest_len;
                let wi = inv_masses[c.i];
                let wj = inv_masses[c.j];
                let denom = wi + wj + c.compliance * alpha;
                if denom < 1e-15 {
                    continue;
                }
                let lambda = -constraint_val / denom;
                let grad = normalize3(d);
                positions[c.i] = add3(pi, scale3(grad, lambda * wi));
                positions[c.j] = add3(pj, scale3(grad, -lambda * wj));
            }
        }
    }
    /// Integrate particle velocities and positions under gravity.
    ///
    /// Classic semi-implicit Euler: apply gravity, update position.
    pub fn integrate(
        &self,
        positions: &mut [[f64; 3]],
        velocities: &mut [[f64; 3]],
        inv_masses: &[f64],
        gravity: [f64; 3],
        dt: f64,
    ) {
        for i in 0..positions.len() {
            if inv_masses[i] < 1e-15 {
                continue;
            }
            velocities[i] = add3(velocities[i], scale3(gravity, dt));
            positions[i] = add3(positions[i], scale3(velocities[i], dt));
        }
    }
}
/// Skinning kernel: applies joint transforms to mesh vertices.
#[derive(Debug)]
pub struct SkinningKernel {
    /// Skinning mode.
    pub mode: SkinningMode,
    /// Joint transforms (world-space).
    pub joints: Vec<JointTransform>,
    /// Dual quaternion joints for DQS mode.
    pub dq_joints: Vec<DualQuat>,
}
impl SkinningKernel {
    /// Create a skinning kernel with the given mode and number of joints.
    pub fn new(mode: SkinningMode, num_joints: usize) -> Self {
        Self {
            mode,
            joints: vec![JointTransform::identity(); num_joints],
            dq_joints: vec![DualQuat::identity(); num_joints],
        }
    }
    /// Apply linear blend skinning to a single vertex.
    pub fn apply_lbs(&self, v: &DeformableVertex) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        let mut total_w = 0.0f64;
        for k in 0..4 {
            let w = v.joint_weights[k];
            if w < 1e-15 {
                continue;
            }
            let ji = v.joint_indices[k] as usize;
            if ji >= self.joints.len() {
                continue;
            }
            let tp = self.joints[ji].transform_point(v.rest_position);
            out = add3(out, scale3(tp, w));
            total_w += w;
        }
        if total_w > 1e-15 {
            scale3(out, 1.0 / total_w)
        } else {
            v.rest_position
        }
    }
    /// Apply dual-quaternion skinning to a single vertex.
    pub fn apply_dqs(&self, v: &DeformableVertex) -> [f64; 3] {
        let mut blended = DualQuat {
            real: [0.0; 4],
            dual: [0.0; 4],
        };
        let mut total_w = 0.0f64;
        for k in 0..4 {
            let w = v.joint_weights[k];
            if w < 1e-15 {
                continue;
            }
            let ji = v.joint_indices[k] as usize;
            if ji >= self.dq_joints.len() {
                continue;
            }
            blended = blended.add_scaled(&self.dq_joints[ji], w);
            total_w += w;
        }
        if total_w < 1e-15 {
            return v.rest_position;
        }
        blended.transform_point(v.rest_position)
    }
    /// Dispatch skinning over the full mesh, updating `vertex.position`.
    ///
    /// This CPU-mock version iterates over all vertices sequentially.
    pub fn dispatch(&self, mesh: &mut DeformableGpuMesh) {
        for v in &mut mesh.vertices {
            v.position = match self.mode {
                SkinningMode::LinearBlend => self.apply_lbs(v),
                SkinningMode::DualQuaternion => self.apply_dqs(v),
            };
        }
    }
}
