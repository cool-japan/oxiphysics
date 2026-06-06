//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Kernel that integrates positions given velocities.
pub struct IntegratePositionKernel;
/// Projected Gauss-Seidel (PGS) constraint solver kernel.
///
/// Iteratively solves positional constraints using position-based dynamics
/// (PBD) style updates.
pub struct ConstraintSolverKernel;
impl ConstraintSolverKernel {
    /// Solve distance constraints for one iteration.
    ///
    /// Updates positions in-place to satisfy the constraints.
    /// `inv_masses[i]` is the inverse mass of body i.
    pub fn solve_distance_constraints(
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        constraints: &[DistanceConstraint],
        dt: f64,
    ) {
        let alpha = 1.0 / (dt * dt);
        for c in constraints {
            let a = c.body_a;
            let b = c.body_b;
            let w_a = inv_masses[a];
            let w_b = inv_masses[b];
            let w_total = w_a + w_b;
            if w_total < 1e-30 {
                continue;
            }
            let dx = [
                positions[b][0] - positions[a][0],
                positions[b][1] - positions[a][1],
                positions[b][2] - positions[a][2],
            ];
            let dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            if dist < 1e-30 {
                continue;
            }
            let err = dist - c.rest_length;
            let tilde_compliance = c.compliance * alpha;
            let lambda = -err / (w_total + tilde_compliance);
            let correction = lambda / dist;
            for k in 0..3 {
                positions[a][k] -= w_a * correction * dx[k];
                positions[b][k] += w_b * correction * dx[k];
            }
        }
    }
}
/// Sleeping parameters used by [`SleepTest`].
#[derive(Debug, Clone, Copy)]
pub struct SleepParams {
    /// Linear velocity threshold below which the body is considered dormant.
    pub linear_threshold: f64,
    /// Angular velocity threshold.
    pub angular_threshold: f64,
    /// Number of consecutive frames below threshold before the body sleeps.
    pub sleep_frames: u32,
}
/// Full rigid body state for one body.
///
/// Orientation is stored as a unit quaternion `[qx, qy, qz, qw]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidBodyState {
    /// Centre-of-mass position.
    pub position: [f64; 3],
    /// Linear velocity.
    pub velocity: [f64; 3],
    /// Orientation quaternion `[qx, qy, qz, qw]`.
    pub orientation: [f64; 4],
    /// Angular velocity in world frame.
    pub angular_velocity: [f64; 3],
    /// Inverse mass (0 = infinite mass / static body).
    pub inverse_mass: f64,
}
impl RigidBodyState {
    /// Construct a `RigidBodyState` at rest.
    pub fn at_rest(position: [f64; 3], inverse_mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            inverse_mass,
        }
    }
}
/// Sleeping state for a rigid body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SleepState {
    /// Body is fully awake and integrated every step.
    Awake,
    /// Body is sleeping and skipped during integration.
    Sleeping,
}
/// Processes a batch of contacts and accumulates impulses.
///
/// This is the CPU-side mock of the GPU contact-batch dispatch kernel.
/// In a real GPU implementation each contact would be handled by one thread;
/// atomic adds would accumulate the impulses into per-body registers.
pub struct ContactBatchProcessor {
    /// Coefficient of restitution for all contacts.
    pub restitution: f64,
    /// Positional correction fraction (Baumgarte stabilisation).
    pub baumgarte: f64,
}
impl ContactBatchProcessor {
    /// Create a new processor.
    pub fn new(restitution: f64, baumgarte: f64) -> Self {
        Self {
            restitution,
            baumgarte,
        }
    }
    /// Process all contacts and return accumulated impulses for each body.
    ///
    /// `n_bodies` is the total body count; the returned `Vec` has one entry
    /// per body, initialised to zero and populated by this function.
    pub fn process_contacts(
        &self,
        contacts: &[ContactPoint],
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        inv_masses: &[f64],
        n_bodies: usize,
        dt: f64,
    ) -> Vec<AccumulatedImpulse> {
        let mut impulses = vec![AccumulatedImpulse::default(); n_bodies];
        for c in contacts {
            let a = c.body_a;
            let b = c.body_b;
            let w_a = inv_masses[a];
            let w_b = inv_masses[b];
            let w_total = w_a + w_b;
            if w_total < 1e-30 {
                continue;
            }
            let corr = self.baumgarte * c.depth / (w_total * dt);
            let rel_vel: f64 = (0..3)
                .map(|k| (velocities[b][k] - velocities[a][k]) * c.normal[k])
                .sum();
            let j = if rel_vel < 0.0 {
                -(1.0 + self.restitution) * rel_vel / w_total
            } else {
                0.0
            };
            let impulse_a = [
                -w_a * (j + corr) * c.normal[0],
                -w_a * (j + corr) * c.normal[1],
                -w_a * (j + corr) * c.normal[2],
            ];
            let impulse_b = [
                w_b * (j + corr) * c.normal[0],
                w_b * (j + corr) * c.normal[1],
                w_b * (j + corr) * c.normal[2],
            ];
            for k in 0..3 {
                impulses[a].linear[k] += impulse_a[k];
                impulses[b].linear[k] += impulse_b[k];
            }
            let ra = [
                c.position[0] - positions[a][0],
                c.position[1] - positions[a][1],
                c.position[2] - positions[a][2],
            ];
            let rb = [
                c.position[0] - positions[b][0],
                c.position[1] - positions[b][1],
                c.position[2] - positions[b][2],
            ];
            let ang_a = [
                ra[1] * impulse_a[2] - ra[2] * impulse_a[1],
                ra[2] * impulse_a[0] - ra[0] * impulse_a[2],
                ra[0] * impulse_a[1] - ra[1] * impulse_a[0],
            ];
            let ang_b = [
                rb[1] * impulse_b[2] - rb[2] * impulse_b[1],
                rb[2] * impulse_b[0] - rb[0] * impulse_b[2],
                rb[0] * impulse_b[1] - rb[1] * impulse_b[0],
            ];
            impulses[a].add_angular(ang_a);
            impulses[b].add_angular(ang_b);
        }
        impulses
    }
}
/// Island (connected component) solver for the constraint graph.
///
/// Bodies connected by constraints form "islands" that can be solved
/// independently, enabling parallel execution.
pub struct IslandSolver {
    /// Island ID for each body. Bodies with the same ID are in the same island.
    pub island_ids: Vec<usize>,
    /// Number of islands.
    pub num_islands: usize,
}
impl IslandSolver {
    /// Build islands from constraints using union-find.
    pub fn build(num_bodies: usize, constraints: &[(usize, usize)]) -> Self {
        let mut parent: Vec<usize> = (0..num_bodies).collect();
        let mut rank = vec![0usize; num_bodies];
        let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        };
        for &(a, b) in constraints {
            let ra = find(&mut parent, a);
            let rb = find(&mut parent, b);
            if ra != rb {
                if rank[ra] < rank[rb] {
                    parent[ra] = rb;
                } else if rank[ra] > rank[rb] {
                    parent[rb] = ra;
                } else {
                    parent[rb] = ra;
                    rank[ra] += 1;
                }
            }
        }
        let mut island_map = std::collections::HashMap::new();
        let mut island_ids = vec![0usize; num_bodies];
        let mut next_id = 0usize;
        for (i, island_id) in island_ids.iter_mut().enumerate() {
            let root = find(&mut parent, i);
            let id = *island_map.entry(root).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            });
            *island_id = id;
        }
        Self {
            island_ids,
            num_islands: next_id,
        }
    }
    /// Get all body indices belonging to a specific island.
    pub fn bodies_in_island(&self, island_id: usize) -> Vec<usize> {
        self.island_ids
            .iter()
            .enumerate()
            .filter(|&(_, id)| *id == island_id)
            .map(|(i, _)| i)
            .collect()
    }
}
/// Tests whether a rigid body should transition to or from sleeping.
pub struct SleepTest {
    /// Per-body consecutive dormant frame counter.
    pub dormant_frames: Vec<u32>,
    /// Per-body current sleep state.
    pub sleep_states: Vec<SleepState>,
}
impl SleepTest {
    /// Create a new sleep test for `n` bodies (all awake initially).
    pub fn new(n: usize) -> Self {
        Self {
            dormant_frames: vec![0; n],
            sleep_states: vec![SleepState::Awake; n],
        }
    }
    /// Update sleep states for all bodies given the current states.
    ///
    /// Returns the number of bodies that transitioned to sleep this step.
    pub fn update(&mut self, states: &[RigidBodyState], params: &SleepParams) -> usize {
        let mut newly_sleeping = 0;
        for (i, s) in states.iter().enumerate() {
            if s.inverse_mass < 1e-30 {
                self.sleep_states[i] = SleepState::Sleeping;
                continue;
            }
            let lin_sq: f64 = s.velocity.iter().map(|v| v * v).sum();
            let ang_sq: f64 = s.angular_velocity.iter().map(|w| w * w).sum();
            let dormant = lin_sq < params.linear_threshold * params.linear_threshold
                && ang_sq < params.angular_threshold * params.angular_threshold;
            if dormant {
                self.dormant_frames[i] += 1;
                if self.dormant_frames[i] >= params.sleep_frames
                    && self.sleep_states[i] == SleepState::Awake
                {
                    self.sleep_states[i] = SleepState::Sleeping;
                    newly_sleeping += 1;
                }
            } else {
                self.dormant_frames[i] = 0;
                self.sleep_states[i] = SleepState::Awake;
            }
        }
        newly_sleeping
    }
    /// Wake all bodies.
    pub fn wake_all(&mut self) {
        for s in &mut self.sleep_states {
            *s = SleepState::Awake;
        }
        for f in &mut self.dormant_frames {
            *f = 0;
        }
    }
    /// Return the number of sleeping bodies.
    pub fn sleeping_count(&self) -> usize {
        self.sleep_states
            .iter()
            .filter(|&&s| s == SleepState::Sleeping)
            .count()
    }
}
/// Accumulated impulse for one body during one solve step.
#[derive(Debug, Clone, Copy, Default)]
pub struct AccumulatedImpulse {
    /// Linear impulse (Δ momentum).
    pub linear: [f64; 3],
    /// Angular impulse (Δ angular momentum).
    pub angular: [f64; 3],
}
impl AccumulatedImpulse {
    /// Add a linear impulse contribution.
    pub fn add_linear(&mut self, impulse: [f64; 3]) {
        for (l, &imp) in self.linear.iter_mut().zip(impulse.iter()) {
            *l += imp;
        }
    }
    /// Add an angular impulse contribution.
    pub fn add_angular(&mut self, impulse: [f64; 3]) {
        for (a, &imp) in self.angular.iter_mut().zip(impulse.iter()) {
            *a += imp;
        }
    }
    /// Apply the accumulated impulse to a body state.
    pub fn apply(&self, state: &mut RigidBodyState) {
        let im = state.inverse_mass;
        for k in 0..3 {
            state.velocity[k] += self.linear[k] * im;
            state.angular_velocity[k] += self.angular[k] * im;
        }
    }
    /// L2 norm of the linear impulse.
    pub fn linear_magnitude(&self) -> f64 {
        self.linear.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
    /// L2 norm of the angular impulse.
    pub fn angular_magnitude(&self) -> f64 {
        self.angular.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}
/// Kernel that integrates velocities given forces and masses.
pub struct IntegrateVelocityKernel;
/// A distance constraint between two bodies.
#[derive(Debug, Clone, Copy)]
pub struct DistanceConstraint {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Rest length of the constraint.
    pub rest_length: f64,
    /// Compliance (inverse stiffness). 0 = rigid.
    pub compliance: f64,
}
/// Rigid body data in Structure-of-Arrays (SoA) layout.
///
/// This layout maps directly to GPU buffer uploads where each quantity
/// (position, velocity, quaternion, etc.) occupies a contiguous memory region,
/// enabling vectorized (SIMD-width) processing in CUDA/WGSL kernels.
pub struct SoaRigidBody {
    /// Number of bodies.
    pub count: usize,
    pub pos_x: Vec<f64>,
    pub pos_y: Vec<f64>,
    pub pos_z: Vec<f64>,
    pub vel_x: Vec<f64>,
    pub vel_y: Vec<f64>,
    pub vel_z: Vec<f64>,
    pub quat_x: Vec<f64>,
    pub quat_y: Vec<f64>,
    pub quat_z: Vec<f64>,
    pub quat_w: Vec<f64>,
    pub omega_x: Vec<f64>,
    pub omega_y: Vec<f64>,
    pub omega_z: Vec<f64>,
    pub inv_mass: Vec<f64>,
}
impl SoaRigidBody {
    /// Build an SoA buffer from a slice of `RigidBodyState` values.
    pub fn from_slice(states: &[RigidBodyState]) -> Self {
        let n = states.len();
        let mut s = Self {
            count: n,
            pos_x: Vec::with_capacity(n),
            pos_y: Vec::with_capacity(n),
            pos_z: Vec::with_capacity(n),
            vel_x: Vec::with_capacity(n),
            vel_y: Vec::with_capacity(n),
            vel_z: Vec::with_capacity(n),
            quat_x: Vec::with_capacity(n),
            quat_y: Vec::with_capacity(n),
            quat_z: Vec::with_capacity(n),
            quat_w: Vec::with_capacity(n),
            omega_x: Vec::with_capacity(n),
            omega_y: Vec::with_capacity(n),
            omega_z: Vec::with_capacity(n),
            inv_mass: Vec::with_capacity(n),
        };
        for rb in states {
            s.pos_x.push(rb.position[0]);
            s.pos_y.push(rb.position[1]);
            s.pos_z.push(rb.position[2]);
            s.vel_x.push(rb.velocity[0]);
            s.vel_y.push(rb.velocity[1]);
            s.vel_z.push(rb.velocity[2]);
            s.quat_x.push(rb.orientation[0]);
            s.quat_y.push(rb.orientation[1]);
            s.quat_z.push(rb.orientation[2]);
            s.quat_w.push(rb.orientation[3]);
            s.omega_x.push(rb.angular_velocity[0]);
            s.omega_y.push(rb.angular_velocity[1]);
            s.omega_z.push(rb.angular_velocity[2]);
            s.inv_mass.push(rb.inverse_mass);
        }
        s
    }
    /// Convert back to a `Vec<RigidBodyState>`.
    pub fn to_vec(&self) -> Vec<RigidBodyState> {
        (0..self.count)
            .map(|i| RigidBodyState {
                position: [self.pos_x[i], self.pos_y[i], self.pos_z[i]],
                velocity: [self.vel_x[i], self.vel_y[i], self.vel_z[i]],
                orientation: [
                    self.quat_x[i],
                    self.quat_y[i],
                    self.quat_z[i],
                    self.quat_w[i],
                ],
                angular_velocity: [self.omega_x[i], self.omega_y[i], self.omega_z[i]],
                inverse_mass: self.inv_mass[i],
            })
            .collect()
    }
    /// Integrate all bodies using Euler stepping (SIMD-friendly loop).
    ///
    /// Processes all bodies in a single loop over the SoA arrays, mimicking
    /// how a GPU kernel would launch one thread per body.
    pub fn integrate_euler(&mut self, forces: &[[f64; 3]], torques: &[[f64; 3]], dt: f64) {
        for i in 0..self.count {
            let im = self.inv_mass[i];
            let ax = forces[i][0] * im;
            let ay = forces[i][1] * im;
            let az = forces[i][2] * im;
            self.vel_x[i] += ax * dt;
            self.vel_y[i] += ay * dt;
            self.vel_z[i] += az * dt;
            self.pos_x[i] += self.vel_x[i] * dt;
            self.pos_y[i] += self.vel_y[i] * dt;
            self.pos_z[i] += self.vel_z[i] * dt;
            let alpha_x = torques[i][0] * im;
            let alpha_y = torques[i][1] * im;
            let alpha_z = torques[i][2] * im;
            self.omega_x[i] += alpha_x * dt;
            self.omega_y[i] += alpha_y * dt;
            self.omega_z[i] += alpha_z * dt;
            let q = [
                self.quat_x[i],
                self.quat_y[i],
                self.quat_z[i],
                self.quat_w[i],
            ];
            let omega = [self.omega_x[i], self.omega_y[i], self.omega_z[i]];
            let dq = quat_derivative(q, omega);
            self.quat_x[i] += dq[0] * dt;
            self.quat_y[i] += dq[1] * dt;
            self.quat_z[i] += dq[2] * dt;
            self.quat_w[i] += dq[3] * dt;
            let qn = quat_normalise([
                self.quat_x[i],
                self.quat_y[i],
                self.quat_z[i],
                self.quat_w[i],
            ]);
            self.quat_x[i] = qn[0];
            self.quat_y[i] = qn[1];
            self.quat_z[i] = qn[2];
            self.quat_w[i] = qn[3];
        }
    }
}
/// Axis-aligned bounding box for broadphase.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: [f64; 3],
    pub max: [f64; 3],
}
impl Aabb {
    /// Create an AABB centred at `position` with `half_extents`.
    pub fn from_center(position: [f64; 3], half_extents: [f64; 3]) -> Self {
        Self {
            min: [
                position[0] - half_extents[0],
                position[1] - half_extents[1],
                position[2] - half_extents[2],
            ],
            max: [
                position[0] + half_extents[0],
                position[1] + half_extents[1],
                position[2] + half_extents[2],
            ],
        }
    }
    /// Check if two AABBs overlap.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Expand the AABB by a margin in all directions.
    pub fn expand(&mut self, margin: f64) {
        for k in 0..3 {
            self.min[k] -= margin;
            self.max[k] += margin;
        }
    }
    /// Compute the volume of the AABB.
    pub fn volume(&self) -> f64 {
        (self.max[0] - self.min[0]) * (self.max[1] - self.min[1]) * (self.max[2] - self.min[2])
    }
}
/// Broadphase update kernel: recomputes AABBs from body positions and
/// finds overlapping pairs using a brute-force sweep.
pub struct BroadphaseUpdateKernel;
impl BroadphaseUpdateKernel {
    /// Update AABBs from body states and find overlapping pairs.
    ///
    /// `half_extents[i]` is the half-extent of body i's AABB.
    /// Returns a list of overlapping body pairs `(i, j)` with `i < j`.
    pub fn find_overlapping_pairs(
        positions: &[[f64; 3]],
        half_extents: &[[f64; 3]],
        margin: f64,
    ) -> Vec<(usize, usize)> {
        let n = positions.len();
        let mut aabbs: Vec<Aabb> = Vec::with_capacity(n);
        for i in 0..n {
            let mut aabb = Aabb::from_center(positions[i], half_extents[i]);
            aabb.expand(margin);
            aabbs.push(aabb);
        }
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if aabbs[i].overlaps(&aabbs[j]) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
}
/// A contact point between two bodies.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// Contact position in world space.
    pub position: [f64; 3],
    /// Contact normal (pointing from body A to body B).
    pub normal: [f64; 3],
    /// Penetration depth (positive means overlap).
    pub depth: f64,
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
}
/// Semi-implicit Euler integration kernel.
///
/// Updates velocity first, then uses the new velocity to update position.
/// This is symplectic and better preserves energy than standard Euler.
pub struct SemiImplicitEulerKernel;
/// A kernel that normalizes a batch of quaternions to unit length.
///
/// In a real CUDA kernel each thread handles one quaternion; this mock
/// processes them in a sequential loop with the same memory access pattern.
pub struct QuaternionNormKernel;
impl QuaternionNormKernel {
    /// Normalize a batch of quaternions to unit length.
    pub fn normalize_batch(quats: &[[f64; 4]]) -> Vec<[f64; 4]> {
        quats.iter().map(|&q| quat_normalise(q)).collect()
    }
    /// In-place normalize a batch of quaternions.
    pub fn normalize_in_place(quats: &mut [[f64; 4]]) {
        for q in quats.iter_mut() {
            *q = quat_normalise(*q);
        }
    }
    /// Count how many quaternions deviate from unit length by more than `tol`.
    pub fn count_denormalized(quats: &[[f64; 4]], tol: f64) -> usize {
        quats
            .iter()
            .filter(|&&q| {
                let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
                (len - 1.0).abs() > tol
            })
            .count()
    }
}
/// Contact generation kernel: generates contact points between
/// sphere-shaped bodies.
pub struct ContactGenerationKernel;
impl ContactGenerationKernel {
    /// Generate contacts between sphere bodies.
    ///
    /// `positions[i]` is the position of body i.
    /// `radii[i]` is the radius of body i.
    /// `pairs` is the list of potentially overlapping pairs from broadphase.
    pub fn generate_sphere_contacts(
        positions: &[[f64; 3]],
        radii: &[f64],
        pairs: &[(usize, usize)],
    ) -> Vec<ContactPoint> {
        let mut contacts = Vec::new();
        for &(a, b) in pairs {
            let dx = [
                positions[b][0] - positions[a][0],
                positions[b][1] - positions[a][1],
                positions[b][2] - positions[a][2],
            ];
            let dist_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            let sum_r = radii[a] + radii[b];
            if dist_sq < sum_r * sum_r {
                let dist = dist_sq.sqrt();
                let depth = sum_r - dist;
                let normal = if dist > 1e-14 {
                    [dx[0] / dist, dx[1] / dist, dx[2] / dist]
                } else {
                    [0.0, 1.0, 0.0]
                };
                let contact_pos = [
                    positions[a][0] + normal[0] * radii[a],
                    positions[a][1] + normal[1] * radii[a],
                    positions[a][2] + normal[2] * radii[a],
                ];
                contacts.push(ContactPoint {
                    position: contact_pos,
                    normal,
                    depth,
                    body_a: a,
                    body_b: b,
                });
            }
        }
        contacts
    }
    /// Apply contact impulses to resolve penetration.
    ///
    /// Uses a simple position-correction approach.
    pub fn resolve_contacts(
        positions: &mut [[f64; 3]],
        velocities: &mut [[f64; 3]],
        inv_masses: &[f64],
        contacts: &[ContactPoint],
        restitution: f64,
    ) {
        for c in contacts {
            let a = c.body_a;
            let b = c.body_b;
            let w_a = inv_masses[a];
            let w_b = inv_masses[b];
            let w_total = w_a + w_b;
            if w_total < 1e-30 {
                continue;
            }
            let correction = c.depth / w_total;
            for (k, &n_k) in c.normal.iter().enumerate() {
                positions[a][k] -= w_a * correction * n_k;
                positions[b][k] += w_b * correction * n_k;
            }
            let rel_vel = [
                velocities[b][0] - velocities[a][0],
                velocities[b][1] - velocities[a][1],
                velocities[b][2] - velocities[a][2],
            ];
            let vel_along_normal =
                rel_vel[0] * c.normal[0] + rel_vel[1] * c.normal[1] + rel_vel[2] * c.normal[2];
            if vel_along_normal < 0.0 {
                let j = -(1.0 + restitution) * vel_along_normal / w_total;
                for (k, &n_k) in c.normal.iter().enumerate() {
                    velocities[a][k] -= w_a * j * n_k;
                    velocities[b][k] += w_b * j * n_k;
                }
            }
        }
    }
}
