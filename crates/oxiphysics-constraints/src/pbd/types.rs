//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use oxiphysics_core::math::Vec3;

#[allow(unused_imports)]
use super::functions::*;

/// XPBD collision (penetration) constraint that projects a particle out of a
/// half-space defined by `(point, normal)`.
///
/// `C = dot(p - point, normal)  ≥  0`
///
/// When violated (`C < 0`), the particle is pushed back along `normal`.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct XpbdCollisionConstraint {
    /// Index of the colliding particle.
    pub particle_idx: usize,
    /// A point on the collision surface.
    pub surface_point: [f64; 3],
    /// Outward surface normal (unit length).
    pub surface_normal: [f64; 3],
    /// XPBD compliance (0 = rigid).
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdCollisionConstraint {
    /// Create a collision constraint.
    pub fn new(
        particle_idx: usize,
        surface_point: [f64; 3],
        surface_normal: [f64; 3],
        compliance: f64,
    ) -> Self {
        Self {
            particle_idx,
            surface_point,
            surface_normal,
            compliance,
        }
    }
    /// Evaluate the constraint: negative when penetrating.
    pub fn evaluate(&self, positions: &[[f64; 3]]) -> f64 {
        let p = positions[self.particle_idx];
        pbd2_dot(pbd2_sub(p, self.surface_point), self.surface_normal)
    }
    /// Project the constraint.  Does nothing when C >= 0 (no penetration).
    pub fn project(
        &self,
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        lambda: &mut f64,
        dt: f64,
    ) {
        let c = self.evaluate(positions);
        if c >= 0.0 {
            return;
        }
        let w = inv_masses[self.particle_idx];
        if w < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = (-c - alpha_tilde * (*lambda)) / (w + alpha_tilde);
        *lambda += d_lambda;
        let correction = pbd2_scale(self.surface_normal, d_lambda * w);
        let p = positions[self.particle_idx];
        positions[self.particle_idx] = pbd2_add(p, correction);
    }
}
/// Monitors convergence of PBD constraint iterations.
#[allow(dead_code)]
pub struct PbdConvergenceMonitor {
    /// Error history (max constraint violation per iteration).
    pub history: Vec<f64>,
    /// Convergence threshold.
    pub threshold: f64,
}
#[allow(dead_code)]
impl PbdConvergenceMonitor {
    /// Create a new monitor.
    pub fn new(threshold: f64) -> Self {
        Self {
            history: Vec::new(),
            threshold,
        }
    }
    /// Compute the current maximum constraint violation.
    pub fn compute_max_violation(
        particles: &[PbdParticle],
        constraints: &[PbdDistanceConstraint],
    ) -> f64 {
        constraints
            .iter()
            .map(|c| {
                let diff = particles[c.b].position - particles[c.a].position;
                (diff.norm() - c.rest_length).abs()
            })
            .fold(0.0_f64, f64::max)
    }
    /// Record the current violation and check convergence.
    pub fn record(&mut self, violation: f64) -> bool {
        self.history.push(violation);
        violation < self.threshold
    }
    /// Convergence rate: ratio of latest to previous violation.
    pub fn convergence_rate(&self) -> Option<f64> {
        let n = self.history.len();
        if n < 2 || self.history[n - 2].abs() < 1e-30 {
            return None;
        }
        Some(self.history[n - 1] / self.history[n - 2])
    }
    /// Reset the monitor.
    pub fn reset(&mut self) {
        self.history.clear();
    }
}
/// XPBD dihedral bending constraint between two triangles sharing an edge.
///
/// Particles: `p0, p1` form the shared edge; `p2` is the apex of triangle A;
/// `p3` is the apex of triangle B.  The rest angle is the dihedral angle
/// between the two triangles at the shared edge.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct XpbdDihedralBending {
    /// Particle indices: \[edge0, edge1, apex_A, apex_B\].
    pub indices: [usize; 4],
    /// Rest dihedral angle (radians).
    pub rest_angle: f64,
    /// XPBD compliance.
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdDihedralBending {
    /// Create a dihedral bending constraint.
    pub fn new(indices: [usize; 4], rest_angle: f64, compliance: f64) -> Self {
        Self {
            indices,
            rest_angle,
            compliance,
        }
    }
    /// Compute the dihedral angle between the two triangles.
    ///
    /// Returns the angle in \[0, pi\].
    pub fn compute_angle(&self, positions: &[[f64; 3]]) -> f64 {
        let p0 = positions[self.indices[0]];
        let p1 = positions[self.indices[1]];
        let p2 = positions[self.indices[2]];
        let p3 = positions[self.indices[3]];
        let e = pbd2_sub(p1, p0);
        let n1 = pbd2_cross(pbd2_sub(p2, p0), e);
        let n2 = pbd2_cross(e, pbd2_sub(p3, p0));
        let l1 = pbd2_norm(n1);
        let l2 = pbd2_norm(n2);
        if l1 < 1e-14 || l2 < 1e-14 {
            return self.rest_angle;
        }
        let cos_a = (pbd2_dot(n1, n2) / (l1 * l2)).clamp(-1.0, 1.0);
        cos_a.acos()
    }
    /// Project the dihedral bending constraint (simplified gradient).
    ///
    /// Uses a linear approximation of the angle gradient suitable for
    /// real-time cloth simulation.
    pub fn project(
        &self,
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        lambda: &mut f64,
        dt: f64,
    ) {
        let angle = self.compute_angle(positions);
        let c = angle - self.rest_angle;
        let w_sum: f64 = self.indices.iter().map(|&i| inv_masses[i]).sum();
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = (-c - alpha_tilde * (*lambda)) / (w_sum + alpha_tilde);
        *lambda += d_lambda;
        let correction = c * d_lambda / (w_sum + 1e-14);
        for &idx in &self.indices {
            let _w = inv_masses[idx];
            let p = positions[idx];
            positions[idx] = [p[0], p[1] - correction * 0.25, p[2]];
        }
    }
}
/// PBD collision constraint that keeps a particle outside a sphere.
pub struct PbdSphereCollision {
    /// Index of the particle subject to this constraint.
    pub particle: usize,
    /// Center of the sphere.
    pub sphere_center: Vec3,
    /// Radius of the sphere.
    pub sphere_radius: f64,
    /// Restitution coefficient.
    pub restitution: f64,
}
impl PbdSphereCollision {
    /// Project the constraint: if the particle is inside the sphere, push it to the surface.
    pub fn project(&self, particles: &mut [PbdParticle]) {
        let p = &mut particles[self.particle];
        if p.inv_mass == 0.0 {
            return;
        }
        let diff = p.position - self.sphere_center;
        let dist = diff.norm();
        if dist < self.sphere_radius {
            if dist < 1e-12 {
                p.position = self.sphere_center + Vec3::new(self.sphere_radius, 0.0, 0.0);
            } else {
                p.position = self.sphere_center + diff / dist * self.sphere_radius;
            }
        }
    }
}
/// Hierarchical PBD solver that runs coarse-to-fine resolution passes.
///
/// The coarse levels provide long-range propagation while fine levels handle
/// local detail — mimicking multigrid approaches.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HierarchicalPbdSolver {
    /// Resolution levels from coarsest to finest.
    pub levels: Vec<HierarchyLevel>,
    /// Global time-step.
    pub dt: f64,
}
#[allow(dead_code)]
impl HierarchicalPbdSolver {
    /// Create a hierarchical solver.
    pub fn new(dt: f64) -> Self {
        Self {
            levels: Vec::new(),
            dt,
        }
    }
    /// Add a level (coarsest first).
    pub fn add_level(&mut self, level: HierarchyLevel) {
        self.levels.push(level);
    }
    /// Total iteration count across all levels.
    pub fn total_iterations(&self) -> u32 {
        self.levels.iter().map(|l| l.iterations).sum()
    }
    /// Run one full hierarchical solve pass, calling `project_fn` at each
    /// level with the adjusted compliance scale and iteration index.
    pub fn solve<F>(&self, mut project_fn: F)
    where
        F: FnMut(usize, f64, u32),
    {
        for (level_idx, level) in self.levels.iter().enumerate() {
            for iter in 0..level.iterations {
                project_fn(level_idx, level.compliance_scale, iter);
            }
        }
    }
}
/// Parallel PBD solver (Jacobi style).
///
/// Accumulates all corrections and applies them at once.
/// This is less accurate per iteration but parallelizable.
#[allow(dead_code)]
pub struct ParallelPbd {
    /// Number of iterations.
    pub iterations: usize,
}
#[allow(dead_code)]
impl ParallelPbd {
    /// Create a new parallel PBD solver.
    pub fn new(iterations: usize) -> Self {
        Self { iterations }
    }
    /// Solve constraints in Jacobi fashion: accumulate then apply.
    pub fn solve(&self, particles: &mut [PbdParticle], constraints: &[PbdDistanceConstraint]) {
        let n = particles.len();
        for _ in 0..self.iterations {
            let mut deltas = vec![Vec3::zeros(); n];
            let mut counts = vec![0_u32; n];
            for c in constraints {
                let wa = particles[c.a].inv_mass;
                let wb = particles[c.b].inv_mass;
                let w_sum = wa + wb;
                if w_sum < 1e-30 {
                    continue;
                }
                let diff = particles[c.b].position - particles[c.a].position;
                let len = diff.norm();
                if len < 1e-12 {
                    continue;
                }
                let n_vec = diff / len;
                let correction = c.stiffness * (len - c.rest_length);
                deltas[c.a] += n_vec * (wa / w_sum * correction);
                deltas[c.b] -= n_vec * (wb / w_sum * correction);
                counts[c.a] += 1;
                counts[c.b] += 1;
            }
            for i in 0..n {
                if counts[i] > 0 {
                    particles[i].position += deltas[i] / counts[i] as f64;
                }
            }
        }
    }
}
/// PBD bending constraint between three particles a, b, c where b is the middle particle.
///
/// Maintains the diagonal distance between a and c to resist bending.
pub struct PbdBendingConstraint {
    /// Index of the first outer particle.
    pub a: usize,
    /// Index of the middle particle.
    pub b: usize,
    /// Index of the second outer particle.
    pub c: usize,
    /// Rest diagonal distance between a and c.
    pub rest_length_ac: f64,
    /// Stiffness in \[0, 1\].
    pub stiffness: f64,
}
impl PbdBendingConstraint {
    /// Create a new bending constraint.
    pub fn new(a: usize, b: usize, c: usize, rest_length_ac: f64, stiffness: f64) -> Self {
        Self {
            a,
            b,
            c,
            rest_length_ac,
            stiffness,
        }
    }
    /// Project the constraint: maintain distance between a and c.
    pub fn project(&self, particles: &mut [PbdParticle]) {
        let wa = particles[self.a].inv_mass;
        let wc = particles[self.c].inv_mass;
        let w_sum = wa + wc;
        if w_sum == 0.0 {
            return;
        }
        let pa = particles[self.a].position;
        let pc = particles[self.c].position;
        let diff = pc - pa;
        let len = diff.norm();
        if len < 1e-12 {
            return;
        }
        let n = diff / len;
        let constraint = len - self.rest_length_ac;
        let delta_a = self.stiffness * (wa / w_sum) * constraint * n;
        let delta_c = -self.stiffness * (wc / w_sum) * constraint * n;
        particles[self.a].position += delta_a;
        particles[self.c].position += delta_c;
    }
}
/// Priority tier for PBD constraints.
///
/// Higher priority constraints are solved more often and before lower priority
/// ones, ensuring critical constraints (e.g., collision) are satisfied first.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PbdConstraintPriority {
    /// Lowest priority: cosmetic (e.g., cloth bending).
    Cosmetic = 0,
    /// Normal priority: structural (e.g., stretch).
    Normal = 1,
    /// High priority: important (e.g., volume conservation).
    High = 2,
    /// Critical priority: must be satisfied (e.g., collision response).
    Critical = 3,
}
/// PBD collision constraint: keep particle above a plane.
pub struct PbdPlaneCollision {
    /// Index of the particle.
    pub particle: usize,
    /// A point on the plane.
    pub point: Vec3,
    /// Outward normal of the plane (unit vector).
    pub normal: Vec3,
}
impl PbdPlaneCollision {
    /// Create a new plane collision constraint.
    pub fn new(particle: usize, point: Vec3, normal: Vec3) -> Self {
        let len = normal.norm();
        let n = if len > 1e-12 {
            normal / len
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        Self {
            particle,
            point,
            normal: n,
        }
    }
    /// Project: if particle is below the plane, push it to the surface.
    pub fn project(&self, particles: &mut [PbdParticle]) {
        let p = &mut particles[self.particle];
        if p.inv_mass == 0.0 {
            return;
        }
        let diff = p.position - self.point;
        let dist = diff.dot(&self.normal);
        if dist < 0.0 {
            p.position -= self.normal * dist;
        }
    }
}
/// PBD simulation managing particles and constraints.
pub struct PbdSimulation {
    /// All particles in the simulation.
    pub particles: Vec<PbdParticle>,
    /// Distance constraints.
    pub distance_constraints: Vec<PbdDistanceConstraint>,
    /// Bending constraints.
    pub bending_constraints: Vec<PbdBendingConstraint>,
    /// Sphere collision constraints.
    pub sphere_collisions: Vec<PbdSphereCollision>,
}
impl PbdSimulation {
    /// Create a new empty PBD simulation.
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
            distance_constraints: Vec::new(),
            bending_constraints: Vec::new(),
            sphere_collisions: Vec::new(),
        }
    }
    /// Add a particle and return its index.
    pub fn add_particle(&mut self, pos: Vec3, mass: f64) -> usize {
        let idx = self.particles.len();
        self.particles.push(PbdParticle::new(pos, mass));
        idx
    }
    /// Pin a particle (set inv_mass = 0, making it static).
    pub fn pin_particle(&mut self, idx: usize) {
        self.particles[idx].inv_mass = 0.0;
        self.particles[idx].velocity = Vec3::zeros();
    }
    /// Add a distance constraint between particles a and b.
    ///
    /// The rest length is computed from the current distance between the particles.
    pub fn add_distance_constraint(&mut self, a: usize, b: usize, stiffness: f64) {
        let rest_length = (self.particles[b].position - self.particles[a].position).norm();
        self.distance_constraints
            .push(PbdDistanceConstraint::new(a, b, rest_length, stiffness));
    }
    /// Add a bending constraint for particles a, b, c (b is the middle).
    pub fn add_bending_constraint(&mut self, a: usize, b: usize, c: usize, stiffness: f64) {
        let rest_length_ac = (self.particles[c].position - self.particles[a].position).norm();
        self.bending_constraints.push(PbdBendingConstraint::new(
            a,
            b,
            c,
            rest_length_ac,
            stiffness,
        ));
    }
    /// Add a sphere collision constraint for a particle.
    pub fn add_sphere_collision(&mut self, particle: usize, center: Vec3, radius: f64) {
        self.sphere_collisions.push(PbdSphereCollision {
            particle,
            sphere_center: center,
            sphere_radius: radius,
            restitution: 0.0,
        });
    }
    /// Advance the simulation by one time step.
    ///
    /// Steps:
    /// 1. Apply external forces (velocity += gravity * dt)
    /// 2. Predict positions (pos += vel * dt), save prev_position
    /// 3. For n_iterations: project all constraints
    /// 4. Update velocities: vel = (pos - prev_pos) / dt
    pub fn step(&mut self, dt: f64, gravity: Vec3, n_iterations: usize) {
        for p in &mut self.particles {
            if p.inv_mass == 0.0 {
                continue;
            }
            p.velocity += gravity * dt;
            p.prev_position = p.position;
            p.position += p.velocity * dt;
        }
        for _ in 0..n_iterations {
            for i in 0..self.distance_constraints.len() {
                let constraint = PbdDistanceConstraint {
                    a: self.distance_constraints[i].a,
                    b: self.distance_constraints[i].b,
                    rest_length: self.distance_constraints[i].rest_length,
                    stiffness: self.distance_constraints[i].stiffness,
                };
                constraint.project(&mut self.particles);
            }
            for i in 0..self.bending_constraints.len() {
                let constraint = PbdBendingConstraint {
                    a: self.bending_constraints[i].a,
                    b: self.bending_constraints[i].b,
                    c: self.bending_constraints[i].c,
                    rest_length_ac: self.bending_constraints[i].rest_length_ac,
                    stiffness: self.bending_constraints[i].stiffness,
                };
                constraint.project(&mut self.particles);
            }
            for i in 0..self.sphere_collisions.len() {
                let constraint = PbdSphereCollision {
                    particle: self.sphere_collisions[i].particle,
                    sphere_center: self.sphere_collisions[i].sphere_center,
                    sphere_radius: self.sphere_collisions[i].sphere_radius,
                    restitution: self.sphere_collisions[i].restitution,
                };
                constraint.project(&mut self.particles);
            }
        }
        if dt > 0.0 {
            for p in &mut self.particles {
                if p.inv_mass == 0.0 {
                    continue;
                }
                p.velocity = (p.position - p.prev_position) / dt;
            }
        }
    }
    /// Total kinetic energy of all particles.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
}
impl PbdSimulation {
    /// Add a volume constraint for a tetrahedron.
    ///
    /// Rest volume is computed from current particle positions.
    pub fn add_volume_constraint(
        &mut self,
        i0: usize,
        i1: usize,
        i2: usize,
        i3: usize,
        stiffness: f64,
    ) {
        let indices = [i0, i1, i2, i3];
        let rest_volume = PbdVolumeConstraint::compute_volume(&self.particles, &indices);
        let _ = PbdVolumeConstraint::new(indices, rest_volume, stiffness);
    }
    /// Create a chain of particles connected by distance constraints.
    ///
    /// Particles are placed along the x-axis with spacing `spacing`.
    /// Returns the indices of the created particles.
    pub fn create_chain(
        &mut self,
        start: Vec3,
        n_particles: usize,
        spacing: f64,
        mass: f64,
        stiffness: f64,
    ) -> Vec<usize> {
        let mut indices = Vec::with_capacity(n_particles);
        for i in 0..n_particles {
            let pos = start + Vec3::new(spacing * i as f64, 0.0, 0.0);
            let idx = self.add_particle(pos, mass);
            indices.push(idx);
        }
        for i in 0..n_particles - 1 {
            self.add_distance_constraint(indices[i], indices[i + 1], stiffness);
        }
        indices
    }
    /// Create a 2D cloth grid of particles connected by distance constraints.
    ///
    /// Grid is in the XY plane, `nx` by `ny` particles with spacing `spacing`.
    /// Returns a 2D array of indices (row-major).
    pub fn create_cloth_grid(
        &mut self,
        origin: Vec3,
        nx: usize,
        ny: usize,
        spacing: f64,
        mass: f64,
        stiffness: f64,
    ) -> Vec<Vec<usize>> {
        let mut grid = vec![vec![0usize; nx]; ny];
        for iy in 0..ny {
            for ix in 0..nx {
                let pos = origin + Vec3::new(spacing * ix as f64, spacing * iy as f64, 0.0);
                grid[iy][ix] = self.add_particle(pos, mass);
            }
        }
        for iy in 0..ny {
            for ix in 0..nx {
                if ix + 1 < nx {
                    self.add_distance_constraint(grid[iy][ix], grid[iy][ix + 1], stiffness);
                }
                if iy + 1 < ny {
                    self.add_distance_constraint(grid[iy][ix], grid[iy + 1][ix], stiffness);
                }
            }
        }
        grid
    }
    /// Advance the simulation with sub-stepping.
    ///
    /// Divides the time step `dt` into `n_substeps` smaller steps for improved stability.
    pub fn step_substep(&mut self, dt: f64, gravity: Vec3, n_iterations: usize, n_substeps: usize) {
        let sub_dt = dt / n_substeps.max(1) as f64;
        for _ in 0..n_substeps {
            self.step(sub_dt, gravity, n_iterations);
        }
    }
    /// Total potential energy from gravity.
    pub fn total_potential_energy(&self, gravity: Vec3) -> f64 {
        let mut pe = 0.0_f64;
        for p in &self.particles {
            if p.inv_mass == 0.0 {
                continue;
            }
            let mass = 1.0 / p.inv_mass;
            pe += -mass * gravity.dot(&p.position);
        }
        pe
    }
    /// Apply velocity damping to all particles.
    pub fn apply_damping(&mut self, damping: f64) {
        for p in &mut self.particles {
            p.velocity *= 1.0 - damping;
        }
    }
    /// Number of particles in the simulation.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }
    /// Number of distance constraints.
    pub fn distance_constraint_count(&self) -> usize {
        self.distance_constraints.len()
    }
}
/// A PBD distance constraint with a priority level.
#[allow(dead_code)]
pub struct PrioritizedConstraint {
    /// The underlying distance constraint.
    pub constraint: PbdDistanceConstraint,
    /// Priority level.
    pub priority: ConstraintPriority,
}
#[allow(dead_code)]
impl PrioritizedConstraint {
    /// Create a new prioritized constraint.
    pub fn new(
        a: usize,
        b: usize,
        rest: f64,
        stiffness: f64,
        priority: ConstraintPriority,
    ) -> Self {
        Self {
            constraint: PbdDistanceConstraint::new(a, b, rest, stiffness),
            priority,
        }
    }
}
/// XPBD distance constraint with compliance (inverse stiffness).
#[allow(dead_code)]
pub struct XpbdDistanceConstraint {
    /// Index of particle A.
    pub a: usize,
    /// Index of particle B.
    pub b: usize,
    /// Rest length.
    pub rest_length: f64,
    /// Physical compliance α = 1/k (m²/N).
    pub compliance: f64,
    /// Accumulated Lagrange multiplier (reset each sub-step).
    pub lambda: f64,
}
#[allow(dead_code)]
impl XpbdDistanceConstraint {
    /// Create a new XPBD constraint.
    pub fn new(a: usize, b: usize, rest_length: f64, compliance: f64) -> Self {
        Self {
            a,
            b,
            rest_length,
            compliance,
            lambda: 0.0,
        }
    }
    /// Reset the Lagrange multiplier (call at start of each sub-step).
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
    /// Project the constraint and update λ.
    ///
    /// Returns the correction delta_lambda.
    pub fn project(&mut self, particles: &mut [PbdParticle], dt: f64) -> f64 {
        let wa = particles[self.a].inv_mass;
        let wb = particles[self.b].inv_mass;
        let w_sum = wa + wb;
        if w_sum < 1e-30 {
            return 0.0;
        }
        let diff = particles[self.b].position - particles[self.a].position;
        let len = diff.norm();
        if len < 1e-12 {
            return 0.0;
        }
        let c = len - self.rest_length;
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = -(c + alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += d_lambda;
        let n_vec = diff / len;
        particles[self.a].position -= n_vec * (wa * d_lambda);
        particles[self.b].position += n_vec * (wb * d_lambda);
        d_lambda
    }
}
/// XPBD pressure constraint enforcing volume conservation for a closed mesh.
///
/// Adds a global volume term: `C = V - V_rest`, where `V` is the current
/// volume approximated from a single tetrahedron.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct XpbdPressureConstraint {
    /// Indices of the four tetrahedron vertices.
    pub tet_indices: [usize; 4],
    /// Rest volume.
    pub rest_volume: f64,
    /// XPBD compliance.
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdPressureConstraint {
    /// Create a pressure constraint.
    pub fn new(tet_indices: [usize; 4], rest_volume: f64, compliance: f64) -> Self {
        Self {
            tet_indices,
            rest_volume,
            compliance,
        }
    }
    /// Signed volume of the tetrahedron (may be negative for inverted tet).
    fn tet_signed_volume(p: [&[f64; 3]; 4]) -> f64 {
        let v1 = pbd2_sub(*p[1], *p[0]);
        let v2 = pbd2_sub(*p[2], *p[0]);
        let v3 = pbd2_sub(*p[3], *p[0]);
        pbd2_dot(v1, pbd2_cross(v2, v3)) / 6.0
    }
    /// Project the pressure constraint for one XPBD sub-step.
    pub fn project(
        &self,
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        lambda: &mut f64,
        dt: f64,
    ) {
        let [i0, i1, i2, i3] = self.tet_indices;
        let p0 = positions[i0];
        let p1 = positions[i1];
        let p2 = positions[i2];
        let p3 = positions[i3];
        let current_vol = Self::tet_signed_volume([&p0, &p1, &p2, &p3]);
        let c = current_vol - self.rest_volume;
        let g0_mag = pbd2_norm(pbd2_sub(p1, p2));
        let g1_mag = pbd2_norm(pbd2_sub(p0, p3));
        let g2_mag = pbd2_norm(pbd2_sub(p3, p0));
        let g3_mag = pbd2_norm(pbd2_sub(p0, p2));
        let w0 = inv_masses[i0];
        let w1 = inv_masses[i1];
        let w2 = inv_masses[i2];
        let w3 = inv_masses[i3];
        let denom = w0 * g0_mag * g0_mag
            + w1 * g1_mag * g1_mag
            + w2 * g2_mag * g2_mag
            + w3 * g3_mag * g3_mag;
        if denom < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = (-c - alpha_tilde * (*lambda)) / (denom + alpha_tilde);
        *lambda += d_lambda;
        let n_dir = if let Some(n) = pbd2_normalize(pbd2_cross(pbd2_sub(p1, p0), pbd2_sub(p2, p0)))
        {
            n
        } else {
            return;
        };
        positions[i0] = pbd2_sub(p0, pbd2_scale(n_dir, d_lambda * w0));
        positions[i1] = pbd2_add(p1, pbd2_scale(n_dir, d_lambda * w1));
        positions[i2] = pbd2_add(p2, pbd2_scale(n_dir, d_lambda * w2));
        positions[i3] = pbd2_add(p3, pbd2_scale(n_dir, d_lambda * w3));
    }
}
/// PBD distance constraint between two particles.
///
/// Stiffness in \[0, 1\]: 0 = no effect, 1 = rigid.
pub struct PbdDistanceConstraint {
    /// Index of the first particle.
    pub a: usize,
    /// Index of the second particle.
    pub b: usize,
    /// Rest (target) length of the constraint.
    pub rest_length: f64,
    /// Stiffness in \[0, 1\].
    pub stiffness: f64,
}
impl PbdDistanceConstraint {
    /// Create a new distance constraint.
    pub fn new(a: usize, b: usize, rest_length: f64, stiffness: f64) -> Self {
        Self {
            a,
            b,
            rest_length,
            stiffness,
        }
    }
    /// Project the constraint, computing and applying position corrections.
    ///
    /// C = |x_b - x_a| - rest_length
    /// Δx_a = +stiffness * w_a/(w_a+w_b) * C * n
    /// Δx_b = -stiffness * w_b/(w_a+w_b) * C * n
    /// n = (x_b - x_a) / |x_b - x_a|
    pub fn project(&self, particles: &mut [PbdParticle]) {
        let wa = particles[self.a].inv_mass;
        let wb = particles[self.b].inv_mass;
        let w_sum = wa + wb;
        if w_sum == 0.0 {
            return;
        }
        let pa = particles[self.a].position;
        let pb = particles[self.b].position;
        let diff = pb - pa;
        let len = diff.norm();
        if len < 1e-12 {
            return;
        }
        let n = diff / len;
        let c = len - self.rest_length;
        let delta_a = self.stiffness * (wa / w_sum) * c * n;
        let delta_b = -self.stiffness * (wb / w_sum) * c * n;
        particles[self.a].position += delta_a;
        particles[self.b].position += delta_b;
    }
}
/// Priority level for PBD constraints.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConstraintPriority {
    /// Low priority (cosmetic, can be violated).
    Low = 0,
    /// Medium priority (soft structural).
    Medium = 1,
    /// High priority (critical structural).
    High = 2,
    /// Critical priority (collision, rigid).
    Critical = 3,
}
/// A directed constraint dependency graph for PBD.
///
/// Nodes are particle indices; edges indicate that two particles share
/// a constraint.
#[allow(dead_code)]
pub struct PbdConstraintGraph {
    /// Number of particles (nodes).
    pub n_particles: usize,
    /// Adjacency list: for each particle, the list of connected particles.
    pub adjacency: Vec<Vec<usize>>,
    /// Constraint count per particle.
    pub degree: Vec<usize>,
}
#[allow(dead_code)]
impl PbdConstraintGraph {
    /// Create a new empty constraint graph.
    pub fn new(n_particles: usize) -> Self {
        Self {
            n_particles,
            adjacency: vec![Vec::new(); n_particles],
            degree: vec![0; n_particles],
        }
    }
    /// Add an edge between particles a and b (undirected).
    pub fn add_edge(&mut self, a: usize, b: usize) {
        if a < self.n_particles && b < self.n_particles && !self.adjacency[a].contains(&b) {
            self.adjacency[a].push(b);
            self.adjacency[b].push(a);
            self.degree[a] += 1;
            self.degree[b] += 1;
        }
    }
    /// Check if two particles are connected.
    pub fn are_connected(&self, a: usize, b: usize) -> bool {
        if a >= self.n_particles {
            return false;
        }
        self.adjacency[a].contains(&b)
    }
    /// Compute graph coloring (greedy) for parallel constraint processing.
    ///
    /// Returns color assignment for each particle (0-indexed colors).
    pub fn greedy_coloring(&self) -> Vec<usize> {
        let n = self.n_particles;
        let mut colors = vec![usize::MAX; n];
        for v in 0..n {
            let forbidden: std::collections::BTreeSet<usize> = self.adjacency[v]
                .iter()
                .filter(|&&u| colors[u] != usize::MAX)
                .map(|&u| colors[u])
                .collect();
            let mut c = 0;
            while forbidden.contains(&c) {
                c += 1;
            }
            colors[v] = c;
        }
        colors
    }
    /// Number of colors used by greedy coloring.
    pub fn chromatic_number_estimate(&self) -> usize {
        let colors = self.greedy_coloring();
        colors.iter().max().map(|&m| m + 1).unwrap_or(0)
    }
}
/// XPBD volume (tetrahedral) constraint for four particles.
#[allow(dead_code)]
pub struct XpbdVolumeConstraint {
    /// Indices of the four tet particles.
    pub indices: [usize; 4],
    /// Rest volume of the tetrahedron.
    pub rest_volume: f64,
    /// Compliance α.
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdVolumeConstraint {
    /// Create a new volume constraint.
    pub fn new(indices: [usize; 4], rest_volume: f64, compliance: f64) -> Self {
        Self {
            indices,
            rest_volume,
            compliance,
        }
    }
    /// Project the volume constraint.
    pub fn project(&self, positions: &mut [[f64; 3]], masses: &[f64], lambda: &mut f64, dt: f64) {
        let [i0, i1, i2, i3] = self.indices;
        let p0 = positions[i0];
        let p1 = positions[i1];
        let p2 = positions[i2];
        let p3 = positions[i3];
        let vol = tet_volume(p0, p1, p2, p3);
        let c = vol - self.rest_volume;
        let g0 = vec3_scale(vec3_cross(vec3_sub(p3, p1), vec3_sub(p2, p1)), 1.0 / 6.0);
        let g1 = vec3_scale(vec3_cross(vec3_sub(p2, p0), vec3_sub(p3, p0)), 1.0 / 6.0);
        let g2 = vec3_scale(vec3_cross(vec3_sub(p3, p0), vec3_sub(p1, p0)), 1.0 / 6.0);
        let g3 = vec3_scale(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)), 1.0 / 6.0);
        let w0 = if masses[i0] > 0.0 {
            1.0 / masses[i0]
        } else {
            0.0
        };
        let w1 = if masses[i1] > 0.0 {
            1.0 / masses[i1]
        } else {
            0.0
        };
        let w2 = if masses[i2] > 0.0 {
            1.0 / masses[i2]
        } else {
            0.0
        };
        let w3 = if masses[i3] > 0.0 {
            1.0 / masses[i3]
        } else {
            0.0
        };
        let w_sum = w0 * vec3_dot(g0, g0)
            + w1 * vec3_dot(g1, g1)
            + w2 * vec3_dot(g2, g2)
            + w3 * vec3_dot(g3, g3);
        if w_sum < 1e-30 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = -(c + alpha_tilde * *lambda) / (w_sum + alpha_tilde);
        *lambda += d_lambda;
        positions[i0] = vec3_add(positions[i0], vec3_scale(g0, w0 * d_lambda));
        positions[i1] = vec3_add(positions[i1], vec3_scale(g1, w1 * d_lambda));
        positions[i2] = vec3_add(positions[i2], vec3_scale(g2, w2 * d_lambda));
        positions[i3] = vec3_add(positions[i3], vec3_scale(g3, w3 * d_lambda));
    }
}
/// Sequential PBD constraint solver (standard Gauss-Seidel).
///
/// Processes constraints one at a time; later constraints can use
/// corrections from earlier ones in the same iteration.
#[allow(dead_code)]
pub struct SequentialPbd {
    /// Number of iterations per time step.
    pub iterations: usize,
}
#[allow(dead_code)]
impl SequentialPbd {
    /// Create a new sequential PBD solver.
    pub fn new(iterations: usize) -> Self {
        Self { iterations }
    }
    /// Solve constraints sequentially.
    pub fn solve(&self, particles: &mut [PbdParticle], constraints: &[PbdDistanceConstraint]) {
        for _ in 0..self.iterations {
            for c in constraints {
                c.project(particles);
            }
        }
    }
    /// Solve with priority ordering: higher priority constraints first.
    pub fn solve_priority(
        &self,
        particles: &mut [PbdParticle],
        constraints: &mut Vec<PrioritizedConstraint>,
    ) {
        constraints.sort_by_key(|b| std::cmp::Reverse(b.priority));
        for _ in 0..self.iterations {
            for pc in constraints.iter() {
                pc.constraint.project(particles);
            }
        }
    }
}
/// XPBD stretch constraint between two particles connected by an edge.
///
/// This is the compliance-based version of the classical distance constraint.
/// With `compliance = 0` it degenerates to the rigid PBD distance constraint.
/// Non-zero compliance allows stiffness tuning independently of iteration
/// count.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct XpbdStretchConstraint {
    /// Index of the first particle.
    pub i: usize,
    /// Index of the second particle.
    pub j: usize,
    /// Rest length of the edge.
    pub rest_length: f64,
    /// XPBD compliance α (inverse stiffness in m/N·s²).  Zero = rigid.
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdStretchConstraint {
    /// Create a stretch constraint.
    pub fn new(i: usize, j: usize, rest_length: f64, compliance: f64) -> Self {
        Self {
            i,
            j,
            rest_length,
            compliance,
        }
    }
    /// Project the constraint for one XPBD sub-step.
    ///
    /// `positions`: mutable slice of 3D positions.
    /// `inv_masses`: inverse mass per particle (0 = fixed).
    /// `lambda`: accumulated Lagrange multiplier (reset each step).
    /// `dt`: sub-step time.
    pub fn project(
        &self,
        positions: &mut [[f64; 3]],
        inv_masses: &[f64],
        lambda: &mut f64,
        dt: f64,
    ) {
        let pi = positions[self.i];
        let pj = positions[self.j];
        let diff = pbd2_sub(pi, pj);
        let len = pbd2_norm(diff);
        if len < 1e-14 {
            return;
        }
        let c = len - self.rest_length;
        let w_i = inv_masses[self.i];
        let w_j = inv_masses[self.j];
        let w_sum = w_i + w_j;
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = (-c - alpha_tilde * (*lambda)) / (w_sum + alpha_tilde);
        *lambda += d_lambda;
        let n = pbd2_scale(diff, 1.0 / len);
        positions[self.i] = pbd2_add(pi, pbd2_scale(n, d_lambda * w_i));
        positions[self.j] = pbd2_sub(pj, pbd2_scale(n, d_lambda * w_j));
    }
}
/// XPBD distance constraint operating on raw `[f64;3]` position arrays.
#[allow(dead_code)]
pub struct XpbdDistanceConstraintRaw {
    /// Index of particle 1.
    pub p1: usize,
    /// Index of particle 2.
    pub p2: usize,
    /// Rest length.
    pub rest_length: f64,
    /// Compliance α (inverse stiffness, m²/N).
    pub compliance: f64,
}
#[allow(dead_code)]
impl XpbdDistanceConstraintRaw {
    /// Create a new raw distance constraint.
    pub fn new(p1: usize, p2: usize, rest_length: f64, compliance: f64) -> Self {
        Self {
            p1,
            p2,
            rest_length,
            compliance,
        }
    }
    /// Project the constraint, updating positions and the Lagrange multiplier.
    pub fn project(&self, positions: &mut [[f64; 3]], masses: &[f64], lambda: &mut f64, dt: f64) {
        let w1 = if masses[self.p1] > 0.0 {
            1.0 / masses[self.p1]
        } else {
            0.0
        };
        let w2 = if masses[self.p2] > 0.0 {
            1.0 / masses[self.p2]
        } else {
            0.0
        };
        let w_sum = w1 + w2;
        if w_sum < 1e-30 {
            return;
        }
        let diff = vec3_sub(positions[self.p2], positions[self.p1]);
        let len = vec3_norm(diff);
        if len < 1e-12 {
            return;
        }
        let c = len - self.rest_length;
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = -(c + alpha_tilde * *lambda) / (w_sum + alpha_tilde);
        *lambda += d_lambda;
        let n = vec3_scale(diff, 1.0 / len);
        positions[self.p1] = vec3_sub(positions[self.p1], vec3_scale(n, w1 * d_lambda));
        positions[self.p2] = vec3_add(positions[self.p2], vec3_scale(n, w2 * d_lambda));
    }
}
/// A PBD particle with position, velocity, and mass.
pub struct PbdParticle {
    /// Current position.
    pub position: Vec3,
    /// Position from the previous time step (used for velocity update).
    pub prev_position: Vec3,
    /// Current velocity.
    pub velocity: Vec3,
    /// Inverse mass (0 = fixed/static particle).
    pub inv_mass: f64,
    /// Restitution (bounce) coefficient in \[0, 1\].
    pub restitution: f64,
}
impl PbdParticle {
    /// Create a new particle at the given position with the given mass.
    pub fn new(pos: Vec3, mass: f64) -> Self {
        Self {
            position: pos,
            prev_position: pos,
            velocity: Vec3::zeros(),
            inv_mass: if mass > 0.0 { 1.0 / mass } else { 0.0 },
            restitution: 0.0,
        }
    }
    /// Create a fixed (static) particle that cannot be moved.
    pub fn new_fixed(pos: Vec3) -> Self {
        Self {
            position: pos,
            prev_position: pos,
            velocity: Vec3::zeros(),
            inv_mass: 0.0,
            restitution: 0.0,
        }
    }
    /// Kinetic energy of this particle: 0.5 * m * |v|².
    pub fn kinetic_energy(&self) -> f64 {
        if self.inv_mass == 0.0 {
            return 0.0;
        }
        let mass = 1.0 / self.inv_mass;
        0.5 * mass * self.velocity.norm_squared()
    }
}
/// A tagged constraint with its priority level.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PrioritizedPbdConstraint {
    /// Priority tier.
    pub priority: PbdConstraintPriority,
    /// Index into the global constraint array.
    pub constraint_index: usize,
    /// Additional weight applied during solving (1.0 = normal).
    pub weight: f64,
}
#[allow(dead_code)]
impl PrioritizedPbdConstraint {
    /// Create a prioritized constraint entry.
    pub fn new(priority: PbdConstraintPriority, constraint_index: usize, weight: f64) -> Self {
        Self {
            priority,
            constraint_index,
            weight: weight.max(0.0),
        }
    }
}
/// Sorts and groups constraints by priority for staged solving.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct PbdConstraintScheduler {
    /// All constraints, will be sorted by priority descending.
    pub constraints: Vec<PrioritizedPbdConstraint>,
}
#[allow(dead_code)]
impl PbdConstraintScheduler {
    /// Create an empty scheduler.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a constraint.
    pub fn add(&mut self, c: PrioritizedPbdConstraint) {
        self.constraints.push(c);
    }
    /// Sort constraints from highest to lowest priority (stable sort).
    pub fn sort_by_priority(&mut self) {
        self.constraints
            .sort_by_key(|b| std::cmp::Reverse(b.priority));
    }
    /// Return the indices of constraints at exactly the given priority.
    pub fn indices_at_priority(&self, priority: PbdConstraintPriority) -> Vec<usize> {
        self.constraints
            .iter()
            .filter(|c| c.priority == priority)
            .map(|c| c.constraint_index)
            .collect()
    }
    /// Number of constraints in total.
    pub fn len(&self) -> usize {
        self.constraints.len()
    }
    /// True when no constraints have been added.
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }
}
/// Isometric bending constraint for a quadrilateral patch.
///
/// Indices: `[0, 1, 2, 3]` where (0,1,2) and (0,2,3) form the two triangles.
/// This constraint penalises deviation from the rest shape's mean curvature.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct IsometricBending {
    /// Particle indices for the quad \[a, b, c, d\].
    pub indices: [usize; 4],
    /// Stiffness coefficient.
    pub stiffness: f64,
    /// Pre-computed rest-pose cotangent weights \[w01, w02, w03, w12, w13, w23\].
    pub cot_weights: [f64; 6],
}
#[allow(dead_code)]
impl IsometricBending {
    /// Create an isometric bending constraint for a flat rest pose quad.
    ///
    /// `positions_rest` provides the four vertex positions in rest pose.
    pub fn new(indices: [usize; 4], positions_rest: [[f64; 3]; 4], stiffness: f64) -> Self {
        let cot_weights = Self::compute_cot_weights(&positions_rest);
        Self {
            indices,
            stiffness,
            cot_weights,
        }
    }
    fn cot_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let ca = pbd2_sub(a, c);
        let cb = pbd2_sub(b, c);
        let cos_c = pbd2_dot(ca, cb);
        let cross = pbd2_cross(ca, cb);
        let sin_c = pbd2_norm(cross);
        if sin_c.abs() < 1e-14 {
            0.0
        } else {
            cos_c / sin_c
        }
    }
    fn compute_cot_weights(p: &[[f64; 3]; 4]) -> [f64; 6] {
        let c012_0 = Self::cot_angle(p[1], p[2], p[0]);
        let c012_1 = Self::cot_angle(p[0], p[2], p[1]);
        let c012_2 = Self::cot_angle(p[0], p[1], p[2]);
        let c023_0 = Self::cot_angle(p[2], p[3], p[0]);
        let c023_2 = Self::cot_angle(p[0], p[3], p[2]);
        let c023_3 = Self::cot_angle(p[0], p[2], p[3]);
        [c012_0, c012_1, c012_2, c023_0, c023_2, c023_3]
    }
    /// Evaluate the bending energy (scalar).
    ///
    /// Returns a non-negative value proportional to squared mean curvature.
    pub fn bending_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let p: [[f64; 3]; 4] = std::array::from_fn(|k| positions[self.indices[k]]);
        let [w0, w1, w2, w3, w4, w5] = self.cot_weights;
        let lap: [f64; 3] = {
            let t0 = pbd2_scale(pbd2_sub(p[1], p[0]), w0);
            let t1 = pbd2_scale(pbd2_sub(p[2], p[0]), w1 + w3);
            let t2 = pbd2_scale(pbd2_sub(p[3], p[0]), w4);
            let t3 = pbd2_scale(pbd2_sub(p[0], p[1]), w2);
            let t4 = pbd2_scale(pbd2_sub(p[0], p[2]), w5);
            let t5 = pbd2_scale(pbd2_sub(p[0], p[3]), w5);
            let s01 = pbd2_add(t0, t1);
            let s23 = pbd2_add(t2, t3);
            let s45 = pbd2_add(t4, t5);
            let s0123 = pbd2_add(s01, s23);
            pbd2_add(s0123, s45)
        };
        0.5 * self.stiffness * pbd2_dot(lap, lap)
    }
}
/// PBD volume constraint for a tetrahedron (4 particles).
///
/// Preserves the signed volume of the tetrahedron. Used for maintaining
/// volume of soft bodies.
pub struct PbdVolumeConstraint {
    /// Indices of the four particles forming the tetrahedron.
    pub indices: [usize; 4],
    /// Rest volume.
    pub rest_volume: f64,
    /// Stiffness in \[0, 1\].
    pub stiffness: f64,
}
impl PbdVolumeConstraint {
    /// Create a new volume constraint.
    pub fn new(indices: [usize; 4], rest_volume: f64, stiffness: f64) -> Self {
        Self {
            indices,
            rest_volume,
            stiffness,
        }
    }
    /// Compute the signed volume of the tetrahedron.
    pub fn compute_volume(p: &[PbdParticle], idx: &[usize; 4]) -> f64 {
        let x0 = p[idx[0]].position;
        let x1 = p[idx[1]].position;
        let x2 = p[idx[2]].position;
        let x3 = p[idx[3]].position;
        let e1 = x1 - x0;
        let e2 = x2 - x0;
        let e3 = x3 - x0;
        let cross = e2.cross(&e3);
        e1.dot(&cross) / 6.0
    }
    /// Project the volume constraint.
    pub fn project(&self, particles: &mut [PbdParticle]) {
        let vol = Self::compute_volume(particles, &self.indices);
        let c = vol - self.rest_volume;
        if c.abs() < 1e-20 {
            return;
        }
        let x0 = particles[self.indices[0]].position;
        let x1 = particles[self.indices[1]].position;
        let x2 = particles[self.indices[2]].position;
        let x3 = particles[self.indices[3]].position;
        let e1 = x1 - x0;
        let e2 = x2 - x0;
        let e3 = x3 - x0;
        let g1 = e2.cross(&e3) / 6.0;
        let g2 = e3.cross(&e1) / 6.0;
        let g3 = e1.cross(&e2) / 6.0;
        let g0 = -(g1 + g2 + g3);
        let grads = [g0, g1, g2, g3];
        let mut w_sum = 0.0f64;
        for k in 0..4 {
            let w = particles[self.indices[k]].inv_mass;
            w_sum += w * grads[k].norm_squared();
        }
        if w_sum < 1e-20 {
            return;
        }
        let lambda = -self.stiffness * c / w_sum;
        for k in 0..4 {
            let w = particles[self.indices[k]].inv_mass;
            particles[self.indices[k]].position += grads[k] * (w * lambda);
        }
    }
}
/// XPBD bending constraint between three particles p1–p2–p3.
#[allow(dead_code)]
pub struct XpbdBendingConstraint {
    /// Index of the first end particle.
    pub p1: usize,
    /// Index of the middle (hinge) particle.
    pub p2: usize,
    /// Index of the second end particle.
    pub p3: usize,
    /// Rest bending angle (radians).
    pub rest_angle: f64,
    /// Compliance α.
    pub compliance: f64,
}
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
impl XpbdBendingConstraint {
    /// Create a new bending constraint.
    pub fn new(p1: usize, p2: usize, p3: usize, rest_angle: f64, compliance: f64) -> Self {
        Self {
            p1,
            p2,
            p3,
            rest_angle,
            compliance,
        }
    }
    /// Project the bending constraint.
    pub fn project(&self, positions: &mut [[f64; 3]], masses: &[f64], lambda: &mut f64, dt: f64) {
        let w1 = if masses[self.p1] > 0.0 {
            1.0 / masses[self.p1]
        } else {
            0.0
        };
        let w2 = if masses[self.p2] > 0.0 {
            1.0 / masses[self.p2]
        } else {
            0.0
        };
        let w3 = if masses[self.p3] > 0.0 {
            1.0 / masses[self.p3]
        } else {
            0.0
        };
        let e1 = vec3_sub(positions[self.p1], positions[self.p2]);
        let e2 = vec3_sub(positions[self.p3], positions[self.p2]);
        let len1 = vec3_norm(e1);
        let len2 = vec3_norm(e2);
        if len1 < 1e-12 || len2 < 1e-12 {
            return;
        }
        let cos_theta = (vec3_dot(e1, e2) / (len1 * len2)).clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let c = theta - self.rest_angle;
        let w_sum = w1 + w2 + w3;
        if w_sum < 1e-30 {
            return;
        }
        let alpha_tilde = self.compliance / (dt * dt);
        let d_lambda = -(c + alpha_tilde * *lambda) / (w_sum + alpha_tilde);
        *lambda += d_lambda;
        let grad = vec3_scale(vec3_sub(e1, e2), d_lambda / (len1 * len2 + 1e-12));
        positions[self.p1] = vec3_add(positions[self.p1], vec3_scale(grad, w1));
        positions[self.p3] = vec3_add(positions[self.p3], vec3_scale(grad, w3));
        let opp = vec3_scale(grad, -(w1 + w3));
        positions[self.p2] = vec3_add(positions[self.p2], vec3_scale(opp, w2));
    }
}
/// A simple PBD sub-step solver that runs multiple sub-steps per frame.
#[allow(dead_code)]
pub struct PbdSubstepSolver {
    /// Number of sub-steps per frame.
    pub substeps: usize,
    /// Full frame timestep (sub-step dt = dt / substeps).
    pub dt: f64,
}
#[allow(dead_code)]
impl PbdSubstepSolver {
    /// Create a new sub-step solver.
    pub fn new(substeps: usize, dt: f64) -> Self {
        Self { substeps, dt }
    }
    /// Advance positions/velocities through `substeps` sub-steps, calling
    /// `constraints` each sub-step to project constraint corrections.
    pub fn solve(
        &self,
        positions: &mut [[f64; 3]],
        velocities: &mut [[f64; 3]],
        masses: &[f64],
        constraints: &mut dyn FnMut(&mut [[f64; 3]]),
    ) {
        let sub_dt = self.dt / self.substeps as f64;
        for _ in 0..self.substeps {
            for i in 0..positions.len() {
                if masses[i] > 0.0 {
                    positions[i] = vec3_add(positions[i], vec3_scale(velocities[i], sub_dt));
                }
            }
            constraints(positions);
        }
    }
}
/// Level-of-detail descriptor for one resolution level in a hierarchical PBD
/// simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HierarchyLevel {
    /// Number of solver iterations at this level.
    pub iterations: u32,
    /// Coarsening ratio from the previous level (e.g., 2 = every 2nd particle).
    pub coarsening: u32,
    /// Compliance multiplier for constraints at this level.
    pub compliance_scale: f64,
}
#[allow(dead_code)]
impl HierarchyLevel {
    /// Create a hierarchy level.
    pub fn new(iterations: u32, coarsening: u32, compliance_scale: f64) -> Self {
        Self {
            iterations,
            coarsening,
            compliance_scale: compliance_scale.max(0.0),
        }
    }
    /// Total particle count at this level given the full-resolution count.
    pub fn particle_count(&self, full_count: usize) -> usize {
        ((full_count as u32).saturating_add(self.coarsening - 1) / self.coarsening) as usize
    }
}
