//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

/// SPH particle storage in Structure-of-Arrays (SoA) layout.
///
/// Keeps each field in its own contiguous Vec so that operations that touch
/// only one field (e.g. density summation) have maximal cache locality.
#[derive(Debug, Clone)]
pub struct ParticleSetSoA {
    /// Particle positions \[\[x, y, z\\]; N].
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities \[\[vx, vy, vz\\]; N].
    pub velocities: Vec<[f64; 3]>,
    /// Net force on each particle \[\[fx, fy, fz\\]; N].
    pub forces: Vec<[f64; 3]>,
    /// Particle masses \[kg\].
    pub masses: Vec<f64>,
    /// Particle densities \[kg/m³\].
    pub densities: Vec<f64>,
    /// Particle pressures \[Pa\].
    pub pressures: Vec<f64>,
}
impl ParticleSetSoA {
    /// Create an empty `ParticleSetSoA` with pre-allocated capacity for `n` particles.
    pub fn new(n: usize) -> Self {
        Self {
            positions: Vec::with_capacity(n),
            velocities: Vec::with_capacity(n),
            forces: Vec::with_capacity(n),
            masses: Vec::with_capacity(n),
            densities: Vec::with_capacity(n),
            pressures: Vec::with_capacity(n),
        }
    }
    /// Append a particle with given position, velocity, and mass.
    ///
    /// Forces, density, and pressure are initialised to zero.
    pub fn add_particle(&mut self, position: [f64; 3], velocity: [f64; 3], mass: f64) {
        self.positions.push(position);
        self.velocities.push(velocity);
        self.forces.push([0.0; 3]);
        self.masses.push(mass);
        self.densities.push(0.0);
        self.pressures.push(0.0);
    }
    /// Remove the particle at index `i` using swap-remove (O(1), does not preserve order).
    ///
    /// # Panics
    /// Panics if `i >= self.len()`.
    pub fn remove_particle(&mut self, i: usize) {
        let last = self.len() - 1;
        self.positions.swap(i, last);
        self.velocities.swap(i, last);
        self.forces.swap(i, last);
        self.masses.swap(i, last);
        self.densities.swap(i, last);
        self.pressures.swap(i, last);
        self.positions.pop();
        self.velocities.pop();
        self.forces.pop();
        self.masses.pop();
        self.densities.pop();
        self.pressures.pop();
    }
    /// Number of particles currently stored.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Returns `true` if there are no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Apply periodic boundary conditions.
    ///
    /// For each particle, wraps coordinates into `[box_min, box_max)` along
    /// each axis using the modulo operation.
    pub fn apply_pbc(&mut self, box_min: [f64; 3], box_max: [f64; 3]) {
        let size = [
            box_max[0] - box_min[0],
            box_max[1] - box_min[1],
            box_max[2] - box_min[2],
        ];
        for pos in &mut self.positions {
            for k in 0..3 {
                if size[k] < 1e-30 {
                    continue;
                }
                let shifted = pos[k] - box_min[k];
                let wrapped = shifted - (shifted / size[k]).floor() * size[k];
                pos[k] = box_min[k] + wrapped;
            }
        }
    }
    /// Forward-Euler integration: `v += (F/m)*dt`, `x += v*dt`.
    pub fn integrate_euler(&mut self, dt: f64) {
        for i in 0..self.len() {
            let inv_m = 1.0 / self.masses[i];
            for k in 0..3 {
                self.velocities[i][k] += self.forces[i][k] * inv_m * dt;
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }
    }
    /// Leapfrog velocity kick: `v += (F/m) * dt`.
    ///
    /// Call this *after* computing forces and *before* `integrate_leapfrog_drift`.
    pub fn integrate_leapfrog_kick(&mut self, dt: f64) {
        for i in 0..self.len() {
            let inv_m = 1.0 / self.masses[i];
            for k in 0..3 {
                self.velocities[i][k] += self.forces[i][k] * inv_m * dt;
            }
        }
    }
    /// Leapfrog position drift: `x += v * dt`.
    ///
    /// Call this *after* `integrate_leapfrog_kick`.
    pub fn integrate_leapfrog_drift(&mut self, dt: f64) {
        for i in 0..self.len() {
            for k in 0..3 {
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }
    }
    /// Compute the mass-weighted centre of mass position.
    ///
    /// Returns `[0, 0, 0]` for an empty set.
    pub fn center_of_mass(&self) -> [f64; 3] {
        if self.is_empty() {
            return [0.0; 3];
        }
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for i in 0..self.len() {
            for k in 0..3 {
                com[k] += self.masses[i] * self.positions[i][k];
            }
        }
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
        com
    }
    /// Total translational kinetic energy: `Σ 0.5 * m_i * |v_i|²`.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| {
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * m * v2
            })
            .sum()
    }
    /// Zero all forces in preparation for force accumulation.
    pub fn clear_forces(&mut self) {
        for f in &mut self.forces {
            *f = [0.0; 3];
        }
    }
    /// Apply a uniform gravitational acceleration `g` to all particles.
    pub fn apply_gravity_soa(&mut self, g: [f64; 3]) {
        for i in 0..self.len() {
            for k in 0..3 {
                self.forces[i][k] += self.masses[i] * g[k];
            }
        }
    }
    /// Total linear momentum: `Σ m_i * v_i`.
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut p = [0.0_f64; 3];
        for i in 0..self.len() {
            for k in 0..3 {
                p[k] += self.masses[i] * self.velocities[i][k];
            }
        }
        p
    }
}
/// A tracker that records the history of a single particle's position over time.
#[allow(dead_code)]
pub struct ParticleTracker {
    /// Index of the tracked particle in the particle set.
    pub particle_idx: usize,
    /// Recorded positions (one per step).
    pub trajectory: Vec<[f64; 3]>,
    /// Recorded simulation times.
    pub times: Vec<f64>,
}
impl ParticleTracker {
    /// Create a new tracker for the particle at `particle_idx`.
    pub fn new(particle_idx: usize) -> Self {
        Self {
            particle_idx,
            trajectory: Vec::new(),
            times: Vec::new(),
        }
    }
    /// Record the current position of the tracked particle.
    pub fn record(&mut self, particles: &SphParticleSet, time: f64) {
        if self.particle_idx < particles.len() {
            self.trajectory.push(particles.positions[self.particle_idx]);
            self.times.push(time);
        }
    }
    /// Total path length of the recorded trajectory.
    pub fn path_length(&self) -> f64 {
        if self.trajectory.len() < 2 {
            return 0.0;
        }
        self.trajectory
            .windows(2)
            .map(|w| {
                let dx = w[1][0] - w[0][0];
                let dy = w[1][1] - w[0][1];
                let dz = w[1][2] - w[0][2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .sum()
    }
    /// Mean speed over the recorded trajectory.
    pub fn mean_speed(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let total_time = self.times[self.times.len() - 1] - self.times[0];
        if total_time < 1e-30 {
            return 0.0;
        }
        self.path_length() / total_time
    }
    /// Clear all recorded data.
    pub fn reset(&mut self) {
        self.trajectory.clear();
        self.times.clear();
    }
}
/// A lazy view into an `SphParticleSet` selecting only particles that pass a
/// user-supplied predicate.
///
/// The filter stores indices into the owning set; iterating gives those
/// indices.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParticleFilter {
    /// Indices of the selected particles.
    pub indices: Vec<usize>,
}
#[allow(dead_code)]
impl ParticleFilter {
    /// Select particles from `ps` where `predicate(i) == true`.
    pub fn new<F>(ps: &SphParticleSet, predicate: F) -> Self
    where
        F: Fn(usize) -> bool,
    {
        let indices = (0..ps.len()).filter(|&i| predicate(i)).collect();
        Self { indices }
    }
    /// Number of selected particles.
    pub fn count(&self) -> usize {
        self.indices.len()
    }
    /// Whether the filter selects no particles.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
    /// Positions of the selected particles.
    pub fn positions<'a>(&self, ps: &'a SphParticleSet) -> Vec<&'a [f64; 3]> {
        self.indices.iter().map(|&i| &ps.positions[i]).collect()
    }
    /// Velocities of the selected particles.
    pub fn velocities<'a>(&self, ps: &'a SphParticleSet) -> Vec<&'a [f64; 3]> {
        self.indices.iter().map(|&i| &ps.velocities[i]).collect()
    }
    /// Total kinetic energy of the selected particles.
    pub fn kinetic_energy(&self, ps: &SphParticleSet) -> f64 {
        self.indices
            .iter()
            .map(|&i| {
                let v = &ps.velocities[i];
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * ps.masses[i] * v2
            })
            .sum()
    }
    /// Centre of mass of the selected particles.
    pub fn center_of_mass(&self, ps: &SphParticleSet) -> [f64; 3] {
        let total_m: f64 = self.indices.iter().map(|&i| ps.masses[i]).sum();
        if total_m == 0.0 {
            return [0.0; 3];
        }
        let mut com = [0.0_f64; 3];
        for &i in &self.indices {
            for k in 0..3 {
                com[k] += ps.masses[i] * ps.positions[i][k];
            }
        }
        for k in 0..3 {
            com[k] /= total_m;
        }
        com
    }
    /// Average density of the selected particles.
    pub fn average_density(&self, ps: &SphParticleSet) -> f64 {
        if self.indices.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.indices.iter().map(|&i| ps.densities[i]).sum();
        sum / self.indices.len() as f64
    }
}
/// Summary statistics computed across an ensemble of `SphParticleSet`
/// snapshots (e.g. multiple simulation replicas or consecutive time steps).
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct EnsembleStats {
    /// Number of snapshots accumulated.
    pub count: usize,
    /// Mean kinetic energy over snapshots.
    pub mean_ke: f64,
    /// Mean total momentum magnitude over snapshots.
    pub mean_mom_mag: f64,
    /// Maximum density error relative to `rho0` over all snapshots.
    pub max_density_error: f64,
    /// Rest density reference value.
    pub rho0: f64,
}
#[allow(dead_code)]
impl EnsembleStats {
    /// Create empty stats with given rest density reference.
    pub fn new(rho0: f64) -> Self {
        Self {
            rho0,
            ..Default::default()
        }
    }
    /// Accumulate statistics from a single snapshot.
    pub fn accumulate(&mut self, ps: &SphParticleSet) {
        self.count += 1;
        let ke = ps.kinetic_energy();
        let n = self.count as f64;
        self.mean_ke = self.mean_ke * (n - 1.0) / n + ke / n;
        let mom = ps.total_momentum();
        let mom_mag = (mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]).sqrt();
        self.mean_mom_mag = self.mean_mom_mag * (n - 1.0) / n + mom_mag / n;
        let rho0 = self.rho0;
        let local_max_err = ps
            .densities
            .iter()
            .map(|&d| (d - rho0).abs() / rho0)
            .fold(0.0_f64, f64::max);
        if local_max_err > self.max_density_error {
            self.max_density_error = local_max_err;
        }
    }
    /// Whether any snapshots have been accumulated.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}
/// A single SPH particle's full data (for initialization convenience).
#[derive(Debug, Clone)]
pub struct SphParticle {
    /// Position of the particle \[m\].
    pub position: Vec3,
    /// Velocity of the particle \[m/s\].
    pub velocity: Vec3,
    /// Acceleration accumulated this step \[m/s²\].
    pub acceleration: Vec3,
    /// Current density \[kg/m³\].
    pub density: f64,
    /// Current pressure \[Pa\].
    pub pressure: f64,
    /// Mass of the particle \[kg\].
    pub mass: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Whether this is a fixed boundary (ghost) particle.
    pub is_boundary: bool,
}
impl SphParticle {
    /// Create a new particle with given position, velocity, and mass.
    ///
    /// All other fields are initialised to zero / false.
    pub fn new(position: Vec3, velocity: Vec3, mass: f64) -> Self {
        Self {
            position,
            velocity,
            acceleration: Vec3::zeros(),
            density: 0.0,
            pressure: 0.0,
            mass,
            viscosity: 0.0,
            is_boundary: false,
        }
    }
    /// Create a boundary particle (fixed, no dynamic update).
    pub fn new_boundary(position: Vec3, mass: f64) -> Self {
        Self {
            position,
            velocity: Vec3::zeros(),
            acceleration: Vec3::zeros(),
            density: 0.0,
            pressure: 0.0,
            mass,
            viscosity: 0.0,
            is_boundary: true,
        }
    }
    /// Kinetic energy of this single particle: ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.velocity.x * self.velocity.x
            + self.velocity.y * self.velocity.y
            + self.velocity.z * self.velocity.z;
        0.5 * self.mass * v2
    }
}
/// Structure-of-arrays particle storage for cache-friendly SPH computations.
///
/// Each field is a `Vec`_`` aligned by particle index.  All vecs always have
/// the same length (invariant maintained by [`add`](Self::add) and
/// [`remove`](Self::remove)).
#[derive(Debug, Clone)]
pub struct SphParticleSet {
    /// Particle positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities \[m/s\].
    pub velocities: Vec<[f64; 3]>,
    /// Accumulated accelerations \[m/s²\].
    pub accelerations: Vec<[f64; 3]>,
    /// Current densities \[kg/m³\].
    pub densities: Vec<f64>,
    /// Current pressures \[Pa\].
    pub pressures: Vec<f64>,
    /// Particle masses \[kg\].
    pub masses: Vec<f64>,
    /// Dynamic viscosities \[Pa·s\].
    pub viscosities: Vec<f64>,
    /// Boundary flags (true = fixed ghost particle).
    pub is_boundary: Vec<bool>,
}
impl SphParticleSet {
    /// Create an empty particle set.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            accelerations: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            masses: Vec::new(),
            viscosities: Vec::new(),
            is_boundary: Vec::new(),
        }
    }
    /// Pre-allocate storage for `n` particles.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            positions: Vec::with_capacity(n),
            velocities: Vec::with_capacity(n),
            accelerations: Vec::with_capacity(n),
            densities: Vec::with_capacity(n),
            pressures: Vec::with_capacity(n),
            masses: Vec::with_capacity(n),
            viscosities: Vec::with_capacity(n),
            is_boundary: Vec::with_capacity(n),
        }
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether the set contains no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Add a single particle from an [`SphParticle`] descriptor.
    pub fn add(&mut self, particle: &SphParticle) {
        self.positions.push([
            particle.position.x,
            particle.position.y,
            particle.position.z,
        ]);
        self.velocities.push([
            particle.velocity.x,
            particle.velocity.y,
            particle.velocity.z,
        ]);
        self.accelerations.push([0.0, 0.0, 0.0]);
        self.densities.push(particle.density);
        self.pressures.push(particle.pressure);
        self.masses.push(particle.mass);
        self.viscosities.push(particle.viscosity);
        self.is_boundary.push(particle.is_boundary);
    }
    /// Convenience: add a particle by raw arrays (no [`SphParticle`] needed).
    pub fn add_raw(
        &mut self,
        pos: [f64; 3],
        vel: [f64; 3],
        mass: f64,
        viscosity: f64,
        boundary: bool,
    ) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.accelerations.push([0.0, 0.0, 0.0]);
        self.densities.push(0.0);
        self.pressures.push(0.0);
        self.masses.push(mass);
        self.viscosities.push(viscosity);
        self.is_boundary.push(boundary);
    }
    /// Remove a particle by index using swap-remove (O(1)).
    pub fn remove(&mut self, index: usize) {
        self.positions.swap_remove(index);
        self.velocities.swap_remove(index);
        self.accelerations.swap_remove(index);
        self.densities.swap_remove(index);
        self.pressures.swap_remove(index);
        self.masses.swap_remove(index);
        self.viscosities.swap_remove(index);
        self.is_boundary.swap_remove(index);
    }
    /// Reset all accelerations to zero.
    pub fn clear_accelerations(&mut self) {
        for a in &mut self.accelerations {
            *a = [0.0, 0.0, 0.0];
        }
    }
    /// Total kinetic energy: Σ ½ mᵢ |vᵢ|².
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| {
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * m * v2
            })
            .sum()
    }
    /// Centre of mass: Σ mᵢ rᵢ / Σ mᵢ.
    ///
    /// Returns `[0,0,0]` for an empty set.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass == 0.0 {
            return [0.0, 0.0, 0.0];
        }
        let mut com = [0.0_f64; 3];
        for (pos, &m) in self.positions.iter().zip(self.masses.iter()) {
            com[0] += m * pos[0];
            com[1] += m * pos[1];
            com[2] += m * pos[2];
        }
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
        com
    }
    /// Apply a body acceleration `g` (gravity) to every non-boundary particle.
    ///
    /// `aᵢ += g`  for fluid particles only.
    pub fn apply_gravity(&mut self, g: [f64; 3]) {
        for (a, &boundary) in self.accelerations.iter_mut().zip(self.is_boundary.iter()) {
            if !boundary {
                a[0] += g[0];
                a[1] += g[1];
                a[2] += g[2];
            }
        }
    }
    /// Euler integration: `v += a·dt`, `r += v·dt` for non-boundary particles.
    pub fn integrate(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            self.velocities[i][0] += self.accelerations[i][0] * dt;
            self.velocities[i][1] += self.accelerations[i][1] * dt;
            self.velocities[i][2] += self.accelerations[i][2] * dt;
            self.positions[i][0] += self.velocities[i][0] * dt;
            self.positions[i][1] += self.velocities[i][1] * dt;
            self.positions[i][2] += self.velocities[i][2] * dt;
        }
    }
    /// Apply periodic boundary conditions.
    ///
    /// Wraps each coordinate into `[box_lo, box_hi)`.
    pub fn apply_pbc(&mut self, box_lo: [f64; 3], box_hi: [f64; 3]) {
        for pos in &mut self.positions {
            for k in 0..3 {
                let len = box_hi[k] - box_lo[k];
                if len <= 0.0 {
                    continue;
                }
                while pos[k] < box_lo[k] {
                    pos[k] += len;
                }
                while pos[k] >= box_hi[k] {
                    pos[k] -= len;
                }
            }
        }
    }
    /// Total momentum: Σ mᵢ vᵢ.
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut p = [0.0_f64; 3];
        for (v, &m) in self.velocities.iter().zip(self.masses.iter()) {
            p[0] += m * v[0];
            p[1] += m * v[1];
            p[2] += m * v[2];
        }
        p
    }
    /// Maximum speed |vᵢ| over all particles.
    pub fn max_speed(&self) -> f64 {
        self.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }
    /// Symplectic Euler integration (kick-drift): `v += a·dt`, `r += v·dt`.
    ///
    /// Alias for [`integrate`](Self::integrate) — the name clarifies the
    /// integration scheme for documentation purposes.
    pub fn integrate_euler(&mut self, dt: f64) {
        self.integrate(dt);
    }
    /// Velocity Verlet integration (full step).
    ///
    /// Performs the position update `r += v·dt + 0.5·a·dt²` and the
    /// half-kick `v += 0.5·a·dt`.  The caller is responsible for
    /// computing forces at the new positions and calling
    /// [`finish_verlet`](Self::finish_verlet) to complete the velocity step.
    pub fn integrate_verlet_half(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                self.positions[i][k] +=
                    self.velocities[i][k] * dt + 0.5 * self.accelerations[i][k] * dt * dt;
            }
            for k in 0..3 {
                self.velocities[i][k] += 0.5 * self.accelerations[i][k] * dt;
            }
        }
    }
    /// Complete the Velocity Verlet step: `v += 0.5·a_new·dt`.
    ///
    /// Call this after re-computing accelerations at the new positions.
    pub fn finish_verlet(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                self.velocities[i][k] += 0.5 * self.accelerations[i][k] * dt;
            }
        }
    }
    /// Leapfrog (kick-drift-kick) integration.
    ///
    /// This is equivalent to Velocity Verlet but expressed in kick-drift form:
    /// 1. half-kick: `v += 0.5·a·dt`
    /// 2. drift:     `r += v·dt`
    /// 3. (caller recomputes accelerations)
    /// 4. half-kick: `v += 0.5·a·dt`
    ///
    /// This method performs steps 1–2.  Call [`finish_verlet`](Self::finish_verlet)
    /// after recomputing accelerations for step 4.
    pub fn leapfrog_kick_drift(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                self.velocities[i][k] += 0.5 * self.accelerations[i][k] * dt;
            }
            for k in 0..3 {
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }
    }
    /// Total potential energy under uniform gravity: Σ mᵢ g · rᵢ.
    ///
    /// `g` is the gravity vector (e.g. `[0, -9.81, 0]`).
    pub fn potential_energy(&self, g: [f64; 3]) -> f64 {
        self.positions
            .iter()
            .zip(self.masses.iter())
            .map(|(r, &m)| -m * (g[0] * r[0] + g[1] * r[1] + g[2] * r[2]))
            .sum()
    }
    /// Total energy (kinetic + gravitational potential).
    pub fn total_energy(&self, g: [f64; 3]) -> f64 {
        self.kinetic_energy() + self.potential_energy(g)
    }
    /// Angular momentum about the origin: Σ mᵢ (rᵢ × vᵢ).
    pub fn angular_momentum(&self) -> [f64; 3] {
        let mut l = [0.0_f64; 3];
        for i in 0..self.len() {
            let r = &self.positions[i];
            let v = &self.velocities[i];
            let m = self.masses[i];
            l[0] += m * (r[1] * v[2] - r[2] * v[1]);
            l[1] += m * (r[2] * v[0] - r[0] * v[2]);
            l[2] += m * (r[0] * v[1] - r[1] * v[0]);
        }
        l
    }
    /// Bounding box of all particle positions: `(min_corner, max_corner)`.
    ///
    /// Returns `([0,0,0], [0,0,0])` for an empty set.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = self.positions[0];
        let mut hi = self.positions[0];
        for pos in &self.positions[1..] {
            for k in 0..3 {
                if pos[k] < lo[k] {
                    lo[k] = pos[k];
                }
                if pos[k] > hi[k] {
                    hi[k] = pos[k];
                }
            }
        }
        (lo, hi)
    }
    /// Average density over all particles.
    pub fn average_density(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.densities.iter().sum();
        sum / self.len() as f64
    }
    /// Maximum density over all particles.
    pub fn max_density(&self) -> f64 {
        self.densities.iter().copied().fold(0.0_f64, f64::max)
    }
    /// Scale all velocities by a factor (useful for damping / rescaling).
    pub fn scale_velocities(&mut self, factor: f64) {
        for v in &mut self.velocities {
            v[0] *= factor;
            v[1] *= factor;
            v[2] *= factor;
        }
    }
    /// Number of fluid (non-boundary) particles.
    pub fn fluid_count(&self) -> usize {
        self.is_boundary.iter().filter(|&&b| !b).count()
    }
    /// Number of boundary particles.
    pub fn boundary_count(&self) -> usize {
        self.is_boundary.iter().filter(|&&b| b).count()
    }
    /// Total mass of all particles.
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }
    /// Clamp all velocities to a maximum magnitude.
    pub fn clamp_velocities(&mut self, max_speed: f64) {
        let max2 = max_speed * max_speed;
        for v in &mut self.velocities {
            let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if v2 > max2 {
                let scale = max_speed / v2.sqrt();
                v[0] *= scale;
                v[1] *= scale;
                v[2] *= scale;
            }
        }
    }
    /// Reflect particles off a wall at `wall_pos` along axis `axis` (0=x, 1=y, 2=z).
    ///
    /// Particles that have crossed below `wall_pos` on the given axis are
    /// reflected back and their velocity component is reversed with the given
    /// coefficient of restitution.
    pub fn reflect_wall(&mut self, axis: usize, wall_pos: f64, restitution: f64) {
        debug_assert!(axis < 3);
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            if self.positions[i][axis] < wall_pos {
                self.positions[i][axis] = 2.0 * wall_pos - self.positions[i][axis];
                self.velocities[i][axis] *= -restitution;
            }
        }
    }
}
impl SphParticleSet {
    /// Compute the mean inter-particle distance to the `k` nearest neighbours
    /// for each fluid particle.
    ///
    /// Returns a vector of `n` values (one per particle).  Boundary particles
    /// return `0.0`.  This is an O(n²) brute-force implementation suitable for
    /// testing and small systems.
    pub fn mean_neighbor_distance(&self, k: usize) -> Vec<f64> {
        let n = self.len();
        let k = k.min(n.saturating_sub(1));
        let mut result = vec![0.0_f64; n];
        if k == 0 {
            return result;
        }
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let pi = &self.positions[i];
            let mut dists: Vec<f64> = (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    let pj = &self.positions[j];
                    let dx = pi[0] - pj[0];
                    let dy = pi[1] - pj[1];
                    let dz = pi[2] - pj[2];
                    (dx * dx + dy * dy + dz * dz).sqrt()
                })
                .collect();
            dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let sum: f64 = dists.iter().take(k).sum();
            result[i] = sum / k as f64;
        }
        result
    }
    /// For each particle, count the number of particles within distance
    /// `radius` (excluding itself).
    pub fn neighbor_counts(&self, radius: f64) -> Vec<usize> {
        let n = self.len();
        let mut counts = vec![0usize; n];
        let r2 = radius * radius;
        for i in 0..n {
            let pi = &self.positions[i];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let pj = &self.positions[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    counts[i] += 1;
                }
            }
        }
        counts
    }
    /// Variance of densities across all particles.
    pub fn density_variance(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let mean = self.average_density();
        let var: f64 = self
            .densities
            .iter()
            .map(|&d| (d - mean).powi(2))
            .sum::<f64>()
            / self.len() as f64;
        var
    }
    /// Standard deviation of densities.
    pub fn density_std(&self) -> f64 {
        self.density_variance().sqrt()
    }
    /// Variance of speeds (|v_i|) across all fluid particles.
    pub fn speed_variance(&self) -> f64 {
        let fluid: Vec<f64> = self
            .velocities
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(v, _)| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .collect();
        if fluid.is_empty() {
            return 0.0;
        }
        let mean = fluid.iter().sum::<f64>() / fluid.len() as f64;
        fluid.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / fluid.len() as f64
    }
}
impl SphParticleSet {
    /// Compute the divergence of the velocity field at each particle using the
    /// SPH discrete sum.
    ///
    /// `h` is the SPH smoothing length.  Uses the cubic spline kernel gradient.
    ///
    /// Returns a vector of `n` divergence estimates.
    pub fn velocity_divergence(&self, h: f64) -> Vec<f64> {
        let n = self.len();
        let mut divv = vec![0.0_f64; n];
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let pi = &self.positions[i];
            let vi = &self.velocities[i];
            let mut d = 0.0_f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let pj = &self.positions[j];
                let vj = &self.velocities[j];
                let rx = pi[0] - pj[0];
                let ry = pi[1] - pj[1];
                let rz = pi[2] - pj[2];
                let r2 = rx * rx + ry * ry + rz * rz;
                let r = r2.sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let dw = cubic_spline_gradient(r, h);
                let dv_x = vi[0] - vj[0];
                let dv_y = vi[1] - vj[1];
                let dv_z = vi[2] - vj[2];
                let m_over_rho = self.masses[j] / self.densities[j].max(1e-20);
                d += m_over_rho * (dv_x * rx / r + dv_y * ry / r + dv_z * rz / r) * dw;
            }
            divv[i] = -d;
        }
        divv
    }
    /// Compute the vorticity magnitude at each particle.
    pub fn vorticity_magnitude(&self, h: f64) -> Vec<f64> {
        let n = self.len();
        let mut omega = vec![0.0_f64; n];
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let pi = &self.positions[i];
            let vi = &self.velocities[i];
            let mut curl = [0.0_f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let pj = &self.positions[j];
                let vj = &self.velocities[j];
                let rx = pi[0] - pj[0];
                let ry = pi[1] - pj[1];
                let rz = pi[2] - pj[2];
                let r = (rx * rx + ry * ry + rz * rz).sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let dw = cubic_spline_gradient(r, h);
                let m_over_rho = self.masses[j] / self.densities[j].max(1e-20);
                let dv = [vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]];
                let r_hat = [rx / r, ry / r, rz / r];
                let gw = [r_hat[0] * dw, r_hat[1] * dw, r_hat[2] * dw];
                curl[0] += m_over_rho * (dv[1] * gw[2] - dv[2] * gw[1]);
                curl[1] += m_over_rho * (dv[2] * gw[0] - dv[0] * gw[2]);
                curl[2] += m_over_rho * (dv[0] * gw[1] - dv[1] * gw[0]);
            }
            omega[i] = (curl[0] * curl[0] + curl[1] * curl[1] + curl[2] * curl[2]).sqrt();
        }
        omega
    }
}
impl SphParticleSet {
    /// Export positions to a flat CSV string: `x,y,z\n` per particle.
    pub fn export_positions_csv(&self) -> String {
        let mut s = String::from("x,y,z\n");
        for p in &self.positions {
            s.push_str(&format!("{},{},{}\n", p[0], p[1], p[2]));
        }
        s
    }
    /// Export full state to a flat CSV string.
    pub fn export_full_csv(&self) -> String {
        let mut s = String::from("x,y,z,vx,vy,vz,density,pressure,mass,is_boundary\n");
        for i in 0..self.len() {
            let p = &self.positions[i];
            let v = &self.velocities[i];
            s.push_str(&format!(
                "{},{},{},{},{},{},{},{},{},{}\n",
                p[0],
                p[1],
                p[2],
                v[0],
                v[1],
                v[2],
                self.densities[i],
                self.pressures[i],
                self.masses[i],
                self.is_boundary[i] as u8,
            ));
        }
        s
    }
    /// Export positions to a flat `Vec`f64` in XYZ order.
    pub fn export_positions_flat(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect()
    }
    /// Import positions from a flat `Vec`f64` in XYZ order.
    pub fn import_positions_flat(&mut self, flat: &[f64]) {
        assert!(
            flat.len().is_multiple_of(3),
            "flat positions must be divisible by 3"
        );
        let n = flat.len() / 3;
        assert_eq!(n, self.len(), "flat positions length mismatch");
        for i in 0..n {
            self.positions[i] = [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]];
        }
    }
    /// Compute per-particle smoothed density from the SPH kernel.
    ///
    /// Uses the cubic spline kernel W(r, h).
    pub fn compute_sph_densities(&mut self, h: f64) {
        let n = self.len();
        for i in 0..n {
            let pi = self.positions[i];
            let mut density = 0.0;
            for j in 0..n {
                let pj = self.positions[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                density += self.masses[j] * cubic_spline_kernel(r, h);
            }
            self.densities[i] = density;
        }
    }
}
impl SphParticleSet {
    /// Merge two particles `i` and `j` (by index) into a single particle at
    /// their combined centre of mass.  The higher-index particle is removed.
    ///
    /// Returns the index of the surviving (merged) particle.
    pub fn merge(&mut self, i: usize, j: usize) -> usize {
        assert!(i < self.len() && j < self.len() && i != j);
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        let m_lo = self.masses[lo];
        let m_hi = self.masses[hi];
        let total = m_lo + m_hi;
        let inv = if total > 1e-30 { 1.0 / total } else { 0.0 };
        for k in 0..3 {
            self.positions[lo][k] =
                (self.positions[lo][k] * m_lo + self.positions[hi][k] * m_hi) * inv;
            self.velocities[lo][k] =
                (self.velocities[lo][k] * m_lo + self.velocities[hi][k] * m_hi) * inv;
        }
        self.masses[lo] = total;
        self.densities[lo] = (self.densities[lo] + self.densities[hi]) * 0.5;
        self.pressures[lo] = (self.pressures[lo] + self.pressures[hi]) * 0.5;
        self.remove(hi);
        lo
    }
    /// Split particle `idx` into two child particles.
    ///
    /// Each child has half the mass of the parent and is displaced by ±`offset`
    /// along `split_axis` (0=x, 1=y, 2=z).
    pub fn split(&mut self, idx: usize, offset: f64, split_axis: usize) {
        assert!(idx < self.len());
        assert!(split_axis < 3);
        let half_mass = self.masses[idx] * 0.5;
        let mut pos_a = self.positions[idx];
        let mut pos_b = self.positions[idx];
        pos_a[split_axis] -= offset;
        pos_b[split_axis] += offset;
        let vel = self.velocities[idx];
        let dens = self.densities[idx];
        let pres = self.pressures[idx];
        let visc = self.viscosities[idx];
        let bound = self.is_boundary[idx];
        self.positions[idx] = pos_a;
        self.masses[idx] = half_mass;
        self.densities[idx] = dens;
        self.pressures[idx] = pres;
        self.positions.push(pos_b);
        self.velocities.push(vel);
        self.accelerations.push([0.0; 3]);
        self.densities.push(dens);
        self.pressures.push(pres);
        self.masses.push(half_mass);
        self.viscosities.push(visc);
        self.is_boundary.push(bound);
    }
}
impl SphParticleSet {
    /// Compute the Morton code for particle `i` given the particle set's
    /// bounding box.  Returns 0 for an out-of-range particle.
    #[allow(dead_code)]
    pub fn morton_code(&self, i: usize, box_lo: [f64; 3], box_hi: [f64; 3]) -> u64 {
        let p = &self.positions[i];
        let bits = (1u32 << 21) - 1;
        let mut coords = [0u32; 3];
        for k in 0..3 {
            let span = box_hi[k] - box_lo[k];
            if span <= 0.0 {
                continue;
            }
            let t = ((p[k] - box_lo[k]) / span).clamp(0.0, 1.0);
            coords[k] = (t * bits as f64) as u32;
        }
        morton_encode(coords[0], coords[1], coords[2])
    }
    /// Sort all particles by their Morton code (Z-order curve).
    ///
    /// This reorders the SoA arrays so that spatially close particles are also
    /// close in memory, improving cache performance for neighbour loops.
    #[allow(dead_code)]
    pub fn sort_by_morton(&mut self, box_lo: [f64; 3], box_hi: [f64; 3]) {
        let n = self.len();
        if n == 0 {
            return;
        }
        let mut codes: Vec<(u64, usize)> = (0..n)
            .map(|i| (self.morton_code(i, box_lo, box_hi), i))
            .collect();
        codes.sort_unstable_by_key(|&(c, _)| c);
        let perm: Vec<usize> = codes.iter().map(|&(_, i)| i).collect();
        fn permute<T: Clone>(v: &[T], perm: &[usize]) -> Vec<T> {
            perm.iter().map(|&i| v[i].clone()).collect()
        }
        self.positions = permute(&self.positions, &perm);
        self.velocities = permute(&self.velocities, &perm);
        self.accelerations = permute(&self.accelerations, &perm);
        self.densities = permute(&self.densities, &perm);
        self.pressures = permute(&self.pressures, &perm);
        self.masses = permute(&self.masses, &perm);
        self.viscosities = permute(&self.viscosities, &perm);
        self.is_boundary = permute(&self.is_boundary, &perm);
    }
    /// Iterate over (index, position) pairs for all fluid (non-boundary)
    /// particles as an immutable slice reference.
    ///
    /// Returns `(i, pos)` where `pos = &self.positions[i]`.
    #[allow(dead_code)]
    pub fn iter_fluid_positions(&self) -> impl Iterator<Item = (usize, &[f64; 3])> {
        self.positions
            .iter()
            .enumerate()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|((i, p), _)| (i, p))
    }
    /// Compute per-particle adaptive smoothing lengths using the standard
    /// SPH criterion `h_i = η * (m_i / ρ_i)^(1/3)`.
    ///
    /// `eta` is typically 1.2–1.5.  Returns a vector of smoothing lengths
    /// (one per particle); boundary particles receive the fallback value `h0`.
    #[allow(dead_code)]
    pub fn adaptive_smoothing_lengths(&self, eta: f64, h0: f64) -> Vec<f64> {
        self.masses
            .iter()
            .zip(self.densities.iter())
            .zip(self.is_boundary.iter())
            .map(|((&m, &rho), &boundary)| {
                if boundary || rho < 1e-20 {
                    h0
                } else {
                    eta * (m / rho).cbrt()
                }
            })
            .collect()
    }
    /// Leapfrog full-step integration (kick-drift-kick) for all fluid
    /// particles with a *single* consistent acceleration field.
    ///
    /// This is the recommended integrator for energy conservation in SPH:
    /// 1. half-kick: `v += 0.5 * a_old * dt`
    /// 2. drift:     `r += v * dt`
    /// 3. (caller recomputes `a` at new positions)
    /// 4. half-kick: `v += 0.5 * a_new * dt`  ← call `finish_verlet`
    ///
    /// This convenience wrapper performs steps 1–2 *and* the final half-kick
    /// using the *current* acceleration, making it suitable for first-order
    /// tests where `a` is held constant.
    #[allow(dead_code)]
    pub fn leapfrog_full(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                self.velocities[i][k] += 0.5 * self.accelerations[i][k] * dt;
            }
            for k in 0..3 {
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
            for k in 0..3 {
                self.velocities[i][k] += 0.5 * self.accelerations[i][k] * dt;
            }
        }
    }
    /// Velocity Verlet full step (constant acceleration variant).
    ///
    /// Performs `r += v*dt + 0.5*a*dt²` then `v += a*dt`.
    /// For variable acceleration problems use `integrate_verlet_half` +
    /// recompute + `finish_verlet` instead.
    #[allow(dead_code)]
    pub fn verlet_full(&mut self, dt: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                self.positions[i][k] +=
                    self.velocities[i][k] * dt + 0.5 * self.accelerations[i][k] * dt * dt;
            }
            for k in 0..3 {
                self.velocities[i][k] += self.accelerations[i][k] * dt;
            }
        }
    }
    /// Radius of gyration of the particle set about its centre of mass.
    ///
    /// `R_g = sqrt( Σ m_i |r_i - r_com|² / Σ m_i )`
    #[allow(dead_code)]
    pub fn radius_of_gyration(&self) -> f64 {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-30 {
            return 0.0;
        }
        let com = self.center_of_mass();
        let weighted_sq: f64 = self
            .positions
            .iter()
            .zip(self.masses.iter())
            .map(|(p, &m)| {
                let dx = p[0] - com[0];
                let dy = p[1] - com[1];
                let dz = p[2] - com[2];
                m * (dx * dx + dy * dy + dz * dz)
            })
            .sum();
        (weighted_sq / total_mass).sqrt()
    }
    /// Compute the temperature proxy: `T_i ∝ |v_i - v_cm|²`.
    ///
    /// Subtracts the centre-of-mass velocity and returns the mass-weighted
    /// mean squared peculiar speed (units: m²/s²).  Useful as a rough proxy
    /// for granular temperature in granular SPH.
    #[allow(dead_code)]
    pub fn granular_temperature(&self) -> f64 {
        let n = self.len();
        if n == 0 {
            return 0.0;
        }
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-30 {
            return 0.0;
        }
        let mut vcm = [0.0_f64; 3];
        for (v, &m) in self.velocities.iter().zip(self.masses.iter()) {
            vcm[0] += m * v[0];
            vcm[1] += m * v[1];
            vcm[2] += m * v[2];
        }
        vcm[0] /= total_mass;
        vcm[1] /= total_mass;
        vcm[2] /= total_mass;
        let mut sum = 0.0_f64;
        for (v, &m) in self.velocities.iter().zip(self.masses.iter()) {
            let dv0 = v[0] - vcm[0];
            let dv1 = v[1] - vcm[1];
            let dv2 = v[2] - vcm[2];
            sum += m * (dv0 * dv0 + dv1 * dv1 + dv2 * dv2);
        }
        sum / total_mass
    }
    /// Add a particle with position, velocity, mass, and smoothing length `h`.
    ///
    /// Convenience wrapper that delegates to `add_raw`, recording `h` via a
    /// returned index and setting viscosity to 0.
    #[allow(dead_code)]
    pub fn add_particle_h(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, _h: f64) -> usize {
        let idx = self.len();
        self.add_raw(pos, vel, mass, 0.0, false);
        idx
    }
    /// Remove a particle by index (explicit swap-remove alias for clarity).
    #[allow(dead_code)]
    pub fn remove_particle(&mut self, index: usize) {
        self.remove(index);
    }
    /// Iterator over all particle indices (fluids and boundaries).
    #[allow(dead_code)]
    pub fn iter_indices(&self) -> std::ops::Range<usize> {
        0..self.len()
    }
    /// Compute pairwise pressure force gradient using the symmetric SPH
    /// formulation:
    ///   `a_press_i += -Σ_j m_j (P_i/ρ_i² + P_j/ρ_j²) ∇W_ij`
    ///
    /// Accumulates into `self.accelerations` (does *not* clear first).
    #[allow(dead_code)]
    pub fn accumulate_pressure_acceleration(&mut self, h: f64) {
        let n = self.len();
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let pi = self.positions[i];
            let rho_i = self.densities[i].max(1e-20);
            let p_i = self.pressures[i];
            let m_i = self.masses[i];
            for j in 0..n {
                if j == i {
                    continue;
                }
                let pj = self.positions[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let rho_j = self.densities[j].max(1e-20);
                let p_j = self.pressures[j];
                let m_j = self.masses[j];
                let dw = cubic_spline_gradient(r, h);
                let factor = -m_j * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw;
                let inv_r = 1.0 / r;
                self.accelerations[i][0] += factor * dx * inv_r;
                self.accelerations[i][1] += factor * dy * inv_r;
                self.accelerations[i][2] += factor * dz * inv_r;
                if !self.is_boundary[j] {
                    let factor_j = m_i * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw;
                    self.accelerations[j][0] += factor_j * dx * inv_r;
                    self.accelerations[j][1] += factor_j * dy * inv_r;
                    self.accelerations[j][2] += factor_j * dz * inv_r;
                }
            }
        }
    }
    /// Morris et al. (1997) viscosity acceleration.
    ///
    /// `a_visc_i += Σ_j m_j * 4 * μ / ((ρ_i + ρ_j) * r²) * (v_ij · r_ij) * ∇W_ij`
    ///
    /// `mu` is the dynamic viscosity \[Pa·s\].  Accumulates into accelerations.
    #[allow(dead_code)]
    pub fn accumulate_morris_viscosity(&mut self, h: f64, mu: f64) {
        if mu.abs() < 1e-30 {
            return;
        }
        let n = self.len();
        for i in 0..n {
            if self.is_boundary[i] {
                continue;
            }
            let pi = self.positions[i];
            let vi = self.velocities[i];
            let rho_i = self.densities[i].max(1e-20);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let pj = self.positions[j];
                let vj = self.velocities[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let rho_j = self.densities[j].max(1e-20);
                let m_j = self.masses[j];
                let dw = cubic_spline_gradient(r, h);
                let dvx = vi[0] - vj[0];
                let dvy = vi[1] - vj[1];
                let dvz = vi[2] - vj[2];
                let coeff = 4.0 * mu / (rho_i + rho_j) * dw / r;
                let factor = m_j * coeff;
                self.accelerations[i][0] += factor * dvx;
                self.accelerations[i][1] += factor * dvy;
                self.accelerations[i][2] += factor * dvz;
            }
        }
    }
    /// WCSPH pressure update: `P_i = B * ((ρ_i / ρ₀)^γ - 1)`.
    ///
    /// `rho0` is the rest density, `B = rho0 * c0² / gamma` is the stiffness,
    /// `c0` is the numerical speed of sound.
    #[allow(dead_code)]
    pub fn update_wcsph_pressure(&mut self, rho0: f64, c0: f64, gamma: f64) {
        let b = rho0 * c0 * c0 / gamma;
        for i in 0..self.len() {
            let ratio = self.densities[i] / rho0;
            self.pressures[i] = b * (ratio.powf(gamma) - 1.0);
        }
    }
    /// Compute SPH density summation with per-particle smoothing length `h_vec`.
    ///
    /// Each particle `i` uses `h_vec[i]` as its kernel radius.
    #[allow(dead_code)]
    pub fn compute_density_variable_h(&mut self, h_vec: &[f64]) {
        let n = self.len();
        assert_eq!(h_vec.len(), n, "h_vec length must match particle count");
        for i in 0..n {
            let hi = h_vec[i];
            let pi = self.positions[i];
            let mut rho = self.masses[i] * cubic_spline_kernel(0.0, hi);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let pj = self.positions[j];
                let dx = pi[0] - pj[0];
                let dy = pi[1] - pj[1];
                let dz = pi[2] - pj[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                rho += self.masses[j] * cubic_spline_kernel(r, hi);
            }
            self.densities[i] = rho;
        }
    }
    /// Compute the CFL-limited time step for this particle set.
    ///
    /// `dt_cfl = cfl * h / (c0 + v_max)` where `v_max` is the maximum speed.
    #[allow(dead_code)]
    pub fn cfl_timestep(&self, h: f64, c0: f64, cfl: f64) -> f64 {
        let v_max = self.max_speed();
        let denom = c0 + v_max;
        if denom < 1e-30 {
            return f64::MAX;
        }
        cfl * h / denom
    }
    /// Viscous time step limit: `dt_visc = 0.125 * h² / nu`.
    #[allow(dead_code)]
    pub fn viscous_timestep(&self, h: f64, nu: f64) -> f64 {
        if nu < 1e-30 {
            return f64::MAX;
        }
        0.125 * h * h / nu
    }
    /// Combined adaptive time step: `min(dt_cfl, dt_visc)` clamped to `[dt_min, dt_max]`.
    #[allow(dead_code)]
    pub fn adaptive_timestep(
        &self,
        h: f64,
        c0: f64,
        nu: f64,
        cfl: f64,
        dt_min: f64,
        dt_max: f64,
    ) -> f64 {
        let dt_c = self.cfl_timestep(h, c0, cfl);
        let dt_v = self.viscous_timestep(h, nu);
        dt_c.min(dt_v).clamp(dt_min, dt_max)
    }
}
impl SphParticleSet {
    /// Compute the root-mean-square speed of all fluid particles.
    #[allow(dead_code)]
    pub fn rms_speed(&self) -> f64 {
        let fluid_count = self.fluid_count();
        if fluid_count == 0 {
            return 0.0;
        }
        let sum_v2: f64 = self
            .velocities
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(v, _)| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum();
        (sum_v2 / fluid_count as f64).sqrt()
    }
    /// Mean velocity vector of all fluid particles.
    #[allow(dead_code)]
    pub fn mean_velocity(&self) -> [f64; 3] {
        let fluid_count = self.fluid_count();
        if fluid_count == 0 {
            return [0.0; 3];
        }
        let mut sum = [0.0_f64; 3];
        for (v, &b) in self.velocities.iter().zip(self.is_boundary.iter()) {
            if !b {
                sum[0] += v[0];
                sum[1] += v[1];
                sum[2] += v[2];
            }
        }
        let n = fluid_count as f64;
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
    /// Mean pressure of all fluid particles.
    #[allow(dead_code)]
    pub fn mean_pressure(&self) -> f64 {
        let fluid_count = self.fluid_count();
        if fluid_count == 0 {
            return 0.0;
        }
        let sum: f64 = self
            .pressures
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(p, _)| p)
            .sum();
        sum / fluid_count as f64
    }
    /// Maximum pressure among all fluid particles.
    #[allow(dead_code)]
    pub fn max_pressure(&self) -> f64 {
        self.pressures
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(p, _)| *p)
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum pressure among all fluid particles.
    #[allow(dead_code)]
    pub fn min_pressure(&self) -> f64 {
        self.pressures
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(p, _)| *p)
            .fold(f64::INFINITY, f64::min)
    }
    /// Pressure standard deviation among all fluid particles.
    #[allow(dead_code)]
    pub fn pressure_std(&self) -> f64 {
        let fluid_count = self.fluid_count();
        if fluid_count == 0 {
            return 0.0;
        }
        let mean = self.mean_pressure();
        let variance: f64 = self
            .pressures
            .iter()
            .zip(self.is_boundary.iter())
            .filter(|&(_, b)| !b)
            .map(|(p, _)| (p - mean).powi(2))
            .sum::<f64>()
            / fluid_count as f64;
        variance.sqrt()
    }
    /// Compute the total (vector) angular momentum about an arbitrary pivot
    /// point `pivot`: Σ mᵢ (rᵢ - pivot) × vᵢ.
    #[allow(dead_code)]
    pub fn angular_momentum_about(&self, pivot: [f64; 3]) -> [f64; 3] {
        let mut l = [0.0_f64; 3];
        for i in 0..self.len() {
            let m = self.masses[i];
            let r = [
                self.positions[i][0] - pivot[0],
                self.positions[i][1] - pivot[1],
                self.positions[i][2] - pivot[2],
            ];
            let v = &self.velocities[i];
            l[0] += m * (r[1] * v[2] - r[2] * v[1]);
            l[1] += m * (r[2] * v[0] - r[0] * v[2]);
            l[2] += m * (r[0] * v[1] - r[1] * v[0]);
        }
        l
    }
    /// Translate all particles by a displacement vector `delta`.
    #[allow(dead_code)]
    pub fn translate(&mut self, delta: [f64; 3]) {
        for pos in &mut self.positions {
            pos[0] += delta[0];
            pos[1] += delta[1];
            pos[2] += delta[2];
        }
    }
    /// Add a uniform velocity offset to all fluid particles.
    #[allow(dead_code)]
    pub fn add_velocity_offset(&mut self, dv: [f64; 3]) {
        for (v, &b) in self.velocities.iter_mut().zip(self.is_boundary.iter()) {
            if !b {
                v[0] += dv[0];
                v[1] += dv[1];
                v[2] += dv[2];
            }
        }
    }
    /// Set all densities to a uniform value `rho0`.
    #[allow(dead_code)]
    pub fn set_uniform_density(&mut self, rho0: f64) {
        for d in &mut self.densities {
            *d = rho0;
        }
    }
    /// Set all pressures to zero.
    #[allow(dead_code)]
    pub fn zero_pressures(&mut self) {
        for p in &mut self.pressures {
            *p = 0.0;
        }
    }
    /// Apply a density-weighted Tait EOS pressure update: P = B*((ρ/ρ₀)^γ - 1).
    ///
    /// `b = ρ₀ * c₀² / γ`.
    #[allow(dead_code)]
    pub fn update_pressures_tait(&mut self, rho0: f64, c0: f64, gamma: f64) {
        let big_b = rho0 * c0 * c0 / gamma;
        for i in 0..self.len() {
            self.pressures[i] = big_b * ((self.densities[i] / rho0).powf(gamma) - 1.0);
        }
    }
    /// SPH density summation using a cubic-spline kernel.
    ///
    /// Updates `self.densities` in-place using the standard summation
    /// `ρᵢ = Σⱼ mⱼ W(|rᵢ - rⱼ|, h)` with cubic spline kernel.
    #[allow(dead_code)]
    pub fn compute_density_sph(&mut self, h: f64) {
        let n = self.len();
        for i in 0..n {
            let mut rho = 0.0_f64;
            for j in 0..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r <= 2.0 * h {
                    rho += self.masses[j] * cubic_spline_w(r, h);
                }
            }
            self.densities[i] = rho;
        }
    }
    /// Count particles in the positive-x half-space (`x > 0`).
    #[allow(dead_code)]
    pub fn count_positive_x_half(&self) -> usize {
        self.positions.iter().filter(|p| p[0] > 0.0).count()
    }
    /// Count particles in a rectangular sub-box `[lo, hi)`.
    #[allow(dead_code)]
    pub fn count_in_box(&self, lo: [f64; 3], hi: [f64; 3]) -> usize {
        self.positions
            .iter()
            .filter(|p| {
                p[0] >= lo[0]
                    && p[0] < hi[0]
                    && p[1] >= lo[1]
                    && p[1] < hi[1]
                    && p[2] >= lo[2]
                    && p[2] < hi[2]
            })
            .count()
    }
    /// Compute the pressure gradient force contribution from the Tait EOS
    /// at particle `i` summed over all neighbours within `2h`.
    ///
    /// Returns the force vector `[fx, fy, fz]`.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn pressure_gradient_force(&self, i: usize, h: f64) -> [f64; 3] {
        let n = self.len();
        let rho_i = self.densities[i].max(1e-20);
        let p_i = self.pressures[i];
        let mut f = [0.0_f64; 3];
        for j in 0..n {
            if j == i {
                continue;
            }
            let dx = self.positions[i][0] - self.positions[j][0];
            let dy = self.positions[i][1] - self.positions[j][1];
            let dz = self.positions[i][2] - self.positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            if r < 1e-14 || r > 2.0 * h {
                continue;
            }
            let rho_j = self.densities[j].max(1e-20);
            let p_j = self.pressures[j];
            let dw = cubic_spline_gradient(r, h);
            let factor = -self.masses[j] * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) * dw / r;
            f[0] += factor * dx;
            f[1] += factor * dy;
            f[2] += factor * dz;
        }
        [
            f[0] * self.masses[i],
            f[1] * self.masses[i],
            f[2] * self.masses[i],
        ]
    }
}
impl SphParticleSet {
    /// Position Verlet: x(t+dt) = 2*x(t) - x(t-dt) + a(t)*dt².
    ///
    /// Requires the previous position array `prev_positions` (same length).
    /// Updates `self.positions` in place; previous positions are also updated
    /// to `self.positions` (the current step becomes the previous for next time).
    #[allow(dead_code)]
    pub fn integrate_verlet_position(&mut self, prev_positions: &mut Vec<[f64; 3]>, dt: f64) {
        assert_eq!(
            prev_positions.len(),
            self.len(),
            "prev_positions length mismatch"
        );
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            let dt2 = dt * dt;
            let new_pos = [
                2.0 * self.positions[i][0] - prev_positions[i][0] + self.accelerations[i][0] * dt2,
                2.0 * self.positions[i][1] - prev_positions[i][1] + self.accelerations[i][1] * dt2,
                2.0 * self.positions[i][2] - prev_positions[i][2] + self.accelerations[i][2] * dt2,
            ];
            self.velocities[i] = [
                (new_pos[0] - prev_positions[i][0]) / (2.0 * dt),
                (new_pos[1] - prev_positions[i][1]) / (2.0 * dt),
                (new_pos[2] - prev_positions[i][2]) / (2.0 * dt),
            ];
            prev_positions[i] = self.positions[i];
            self.positions[i] = new_pos;
        }
    }
    /// Apply reflective (bounce) boundary conditions on all six faces of a box.
    ///
    /// Particles that cross a wall are reflected and their normal velocity
    /// component is negated with a coefficient of restitution `e ∈ [0, 1]`.
    #[allow(dead_code)]
    pub fn apply_box_reflection(&mut self, box_lo: [f64; 3], box_hi: [f64; 3], restitution: f64) {
        for i in 0..self.len() {
            if self.is_boundary[i] {
                continue;
            }
            for k in 0..3 {
                if self.positions[i][k] < box_lo[k] {
                    self.positions[i][k] = 2.0 * box_lo[k] - self.positions[i][k];
                    self.velocities[i][k] *= -restitution;
                } else if self.positions[i][k] > box_hi[k] {
                    self.positions[i][k] = 2.0 * box_hi[k] - self.positions[i][k];
                    self.velocities[i][k] *= -restitution;
                }
            }
        }
    }
    /// Enumerate particles sorted by distance from a reference point `ref_pos`.
    ///
    /// Returns a vector of `(index, distance)` pairs sorted ascending by
    /// distance.
    #[allow(dead_code)]
    pub fn sorted_by_distance_from(&self, ref_pos: [f64; 3]) -> Vec<(usize, f64)> {
        let mut pairs: Vec<(usize, f64)> = self
            .positions
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let dx = p[0] - ref_pos[0];
                let dy = p[1] - ref_pos[1];
                let dz = p[2] - ref_pos[2];
                (i, (dx * dx + dy * dy + dz * dz).sqrt())
            })
            .collect();
        pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        pairs
    }
    /// Number of particles closer than `r` to a given point `center`.
    #[allow(dead_code)]
    pub fn count_within_radius(&self, center: [f64; 3], r: f64) -> usize {
        let r2 = r * r;
        self.positions
            .iter()
            .filter(|p| {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                dx * dx + dy * dy + dz * dz < r2
            })
            .count()
    }
    /// Compute the density-weighted velocity field divergence at particle `i`
    /// using an approximate SPH discretization with cubic-spline kernel.
    #[allow(dead_code)]
    pub fn divergence_at(&self, i: usize, h: f64) -> f64 {
        let n = self.len();
        let pi = &self.positions[i];
        let vi = &self.velocities[i];
        let mut div = 0.0_f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = pi[0] - self.positions[j][0];
            let dy = pi[1] - self.positions[j][1];
            let dz = pi[2] - self.positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            if r < 1e-14 || r > 2.0 * h {
                continue;
            }
            let dw = cubic_spline_gradient(r, h);
            let dv = [
                self.velocities[j][0] - vi[0],
                self.velocities[j][1] - vi[1],
                self.velocities[j][2] - vi[2],
            ];
            let m_rho = self.masses[j] / self.densities[j].max(1e-20);
            let dot = dv[0] * dx / r + dv[1] * dy / r + dv[2] * dz / r;
            div += m_rho * dot * dw;
        }
        -div
    }
    /// Return a `Vec`usize` of particle indices sorted by the given field:
    /// `"density"`, `"pressure"`, `"speed"`, or `"mass"`.
    ///
    /// Returns an empty vec on unknown field names.
    #[allow(dead_code)]
    pub fn argsort_by_field(&self, field: &str) -> Vec<usize> {
        let n = self.len();
        let mut idx: Vec<usize> = (0..n).collect();
        match field {
            "density" => idx.sort_unstable_by(|&a, &b| {
                self.densities[a]
                    .partial_cmp(&self.densities[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
            "pressure" => idx.sort_unstable_by(|&a, &b| {
                self.pressures[a]
                    .partial_cmp(&self.pressures[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
            "speed" => idx.sort_unstable_by(|&a, &b| {
                let sa = (self.velocities[a][0].powi(2)
                    + self.velocities[a][1].powi(2)
                    + self.velocities[a][2].powi(2))
                .sqrt();
                let sb = (self.velocities[b][0].powi(2)
                    + self.velocities[b][1].powi(2)
                    + self.velocities[b][2].powi(2))
                .sqrt();
                sa.partial_cmp(&sb).unwrap_or(std::cmp::Ordering::Equal)
            }),
            "mass" => idx.sort_unstable_by(|&a, &b| {
                self.masses[a]
                    .partial_cmp(&self.masses[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
            _ => {}
        }
        idx
    }
    /// Verify SoA invariant: all arrays have the same length.
    ///
    /// Returns `Ok(())` if consistent, `Err(msg)` otherwise.
    #[allow(dead_code)]
    pub fn check_invariant(&self) -> Result<(), String> {
        let n = self.positions.len();
        macro_rules! check {
            ($field:ident) => {
                if self.$field.len() != n {
                    return Err(format!(
                        "SoA invariant violated: positions.len()={} but {}.len()={}",
                        n,
                        stringify!($field),
                        self.$field.len()
                    ));
                }
            };
        }
        check!(velocities);
        check!(accelerations);
        check!(densities);
        check!(pressures);
        check!(masses);
        check!(viscosities);
        check!(is_boundary);
        Ok(())
    }
}
/// Structure-of-arrays particle storage for cache-friendly SPH computations.
///
/// This legacy type uses [`Vec3`] for positions/velocities/forces and is
/// retained for compatibility with code that already uses it.
#[derive(Debug, Clone)]
pub struct ParticleSet {
    /// Particle positions.
    pub positions: Vec<Vec3>,
    /// Particle velocities.
    pub velocities: Vec<Vec3>,
    /// Particle densities (computed each step).
    pub densities: Vec<f64>,
    /// Particle pressures (computed from densities).
    pub pressures: Vec<f64>,
    /// Accumulated forces on each particle.
    pub forces: Vec<Vec3>,
    /// Particle masses.
    pub masses: Vec<f64>,
}
impl ParticleSet {
    /// Create an empty particle set.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            pressures: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
        }
    }
    /// Create a particle set with pre-allocated capacity.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            positions: Vec::with_capacity(n),
            velocities: Vec::with_capacity(n),
            densities: Vec::with_capacity(n),
            pressures: Vec::with_capacity(n),
            forces: Vec::with_capacity(n),
            masses: Vec::with_capacity(n),
        }
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether the particle set is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Add a single particle.
    pub fn add_particle(&mut self, particle: &SphParticle) {
        self.positions.push(particle.position);
        self.velocities.push(particle.velocity);
        self.densities.push(0.0);
        self.pressures.push(0.0);
        self.forces.push(Vec3::zeros());
        self.masses.push(particle.mass);
    }
    /// Remove a particle by index (swap-remove for O(1)).
    pub fn remove(&mut self, index: usize) {
        self.positions.swap_remove(index);
        self.velocities.swap_remove(index);
        self.densities.swap_remove(index);
        self.pressures.swap_remove(index);
        self.forces.swap_remove(index);
        self.masses.swap_remove(index);
    }
    /// Reset all forces to zero.
    pub fn clear_forces(&mut self) {
        for f in &mut self.forces {
            *f = Vec3::zeros();
        }
    }
}
