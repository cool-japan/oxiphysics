//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// A single particle used in volumetric PBD simulation.
#[derive(Debug, Clone)]
pub struct VolumeParticle {
    /// Current world-space position.
    pub position: [f64; 3],
    /// Position at the previous time step (used for Verlet integration).
    pub prev_position: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// When `true` the particle is pinned and never moved by the solver.
    pub fixed: bool,
}
/// PBD tetrahedral soft body.
///
/// Simulation loop:
/// 1. `step` integrates positions with Verlet and gravity.
/// 2. `solve_volume_constraints` enforces incompressibility.
/// 3. `solve_edge_constraints` enforces rest edge lengths.
#[derive(Debug, Clone)]
pub struct TetrahedralSoftBody {
    /// All particles in the mesh.
    pub particles: Vec<VolumeParticle>,
    /// Each entry contains four particle indices forming a tetrahedron.
    pub tetrahedra: Vec<[usize; 4]>,
    /// Rest volume for each tetrahedron (same order as `tetrahedra`).
    pub rest_volumes: Vec<f64>,
    /// Stiffness coefficient for the volume constraint (0 = fully compliant, 1 = rigid).
    pub volume_stiffness: f64,
    /// Stiffness coefficient for edge distance constraints.
    pub edge_stiffness: f64,
    /// `(i, j, rest_len)` — one entry per unique edge in the mesh.
    pub rest_lengths: Vec<(usize, usize, f64)>,
}
impl TetrahedralSoftBody {
    /// Build a new soft body.
    ///
    /// Rest volumes and edge rest lengths are computed automatically from the
    /// initial particle positions.
    pub fn new(particles: Vec<VolumeParticle>, tetrahedra: Vec<[usize; 4]>) -> Self {
        let rest_volumes: Vec<f64> = tetrahedra
            .iter()
            .map(|tet| {
                let v = tet_volume_signed(
                    particles[tet[0]].position,
                    particles[tet[1]].position,
                    particles[tet[2]].position,
                    particles[tet[3]].position,
                );
                v.abs()
            })
            .collect();
        use std::collections::HashSet;
        let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
        for tet in &tetrahedra {
            let edges = [
                (tet[0], tet[1]),
                (tet[0], tet[2]),
                (tet[0], tet[3]),
                (tet[1], tet[2]),
                (tet[1], tet[3]),
                (tet[2], tet[3]),
            ];
            for (a, b) in edges {
                let key = if a < b { (a, b) } else { (b, a) };
                edge_set.insert(key);
            }
        }
        let rest_lengths: Vec<(usize, usize, f64)> = edge_set
            .into_iter()
            .map(|(a, b)| {
                let d = len3(sub3(particles[a].position, particles[b].position));
                (a, b, d)
            })
            .collect();
        Self {
            particles,
            tetrahedra,
            rest_volumes,
            volume_stiffness: 0.5,
            edge_stiffness: 0.5,
            rest_lengths,
        }
    }
    /// Signed volume of a single tetrahedron given a particle slice.
    pub fn tet_volume(p: &[VolumeParticle], tet: &[usize; 4]) -> f64 {
        tet_volume_signed(
            p[tet[0]].position,
            p[tet[1]].position,
            p[tet[2]].position,
            p[tet[3]].position,
        )
        .abs()
    }
    /// Sum of all tetrahedral volumes in the current configuration.
    pub fn total_volume(&self) -> f64 {
        self.tetrahedra
            .iter()
            .map(|tet| Self::tet_volume(&self.particles, tet))
            .sum()
    }
    /// Project volume constraints for all tetrahedra.
    ///
    /// Constraint: `C = V_current - V_rest = 0`.
    ///
    /// Gradient (Green's theorem): each vertex gradient is `±(1/6) * (e_j × e_k)`.
    pub fn solve_volume_constraints(&mut self) {
        for (idx, &tet) in self.tetrahedra.iter().enumerate() {
            let v_rest = self.rest_volumes[idx];
            let p0 = self.particles[tet[0]].position;
            let p1 = self.particles[tet[1]].position;
            let p2 = self.particles[tet[2]].position;
            let p3 = self.particles[tet[3]].position;
            let v_cur = tet_volume_signed(p0, p1, p2, p3).abs();
            let c = v_cur - v_rest;
            if c.abs() < 1e-15 {
                continue;
            }
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);
            let e3 = sub3(p3, p0);
            let g0 = cross3(sub3(p2, p1), sub3(p3, p1));
            let g1 = cross3(e2, e3);
            let g2 = cross3(e3, e1);
            let g3 = cross3(e1, e2);
            let inv6 = 1.0 / 6.0;
            let g0 = scale3(g0, inv6);
            let g1 = scale3(g1, inv6);
            let g2 = scale3(g2, inv6);
            let g3 = scale3(g3, inv6);
            let w = [
                if self.particles[tet[0]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[0]].mass
                },
                if self.particles[tet[1]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[1]].mass
                },
                if self.particles[tet[2]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[2]].mass
                },
                if self.particles[tet[3]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[3]].mass
                },
            ];
            let denom = w[0] * dot3(g0, g0)
                + w[1] * dot3(g1, g1)
                + w[2] * dot3(g2, g2)
                + w[3] * dot3(g3, g3);
            if denom < 1e-15 {
                continue;
            }
            let lambda = -self.volume_stiffness * c / denom;
            let gs = [g0, g1, g2, g3];
            for (k, &ti) in tet.iter().enumerate() {
                if !self.particles[ti].fixed {
                    let dp = scale3(gs[k], lambda * w[k]);
                    self.particles[ti].position = add3(self.particles[ti].position, dp);
                }
            }
        }
    }
    /// Project PBD distance constraints for every edge in the mesh.
    pub fn solve_edge_constraints(&mut self) {
        for &(i, j, rest_len) in &self.rest_lengths {
            let pi = self.particles[i].position;
            let pj = self.particles[j].position;
            let diff = sub3(pi, pj);
            let cur_len = len3(diff);
            if cur_len < 1e-15 {
                continue;
            }
            let c = cur_len - rest_len;
            if c.abs() < 1e-15 {
                continue;
            }
            let wi = if self.particles[i].fixed {
                0.0
            } else {
                1.0 / self.particles[i].mass
            };
            let wj = if self.particles[j].fixed {
                0.0
            } else {
                1.0 / self.particles[j].mass
            };
            let denom = wi + wj;
            if denom < 1e-15 {
                continue;
            }
            let dir = scale3(diff, 1.0 / cur_len);
            let lambda = -self.edge_stiffness * c / denom;
            if !self.particles[i].fixed {
                let dp = scale3(dir, lambda * wi);
                self.particles[i].position = add3(self.particles[i].position, dp);
            }
            if !self.particles[j].fixed {
                let dp = scale3(dir, -lambda * wj);
                self.particles[j].position = add3(self.particles[j].position, dp);
            }
        }
    }
    /// Advance the simulation by one time step.
    ///
    /// Applies Verlet integration with external gravity, then projects both
    /// volume and edge constraints.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let dt2 = dt * dt;
        for p in &mut self.particles {
            if p.fixed {
                continue;
            }
            let vel = sub3(p.position, p.prev_position);
            let accel = gravity;
            let new_pos = add3(add3(p.position, vel), scale3(accel, dt2));
            p.prev_position = p.position;
            p.position = new_pos;
        }
        self.solve_volume_constraints();
        self.solve_edge_constraints();
    }
    /// Apply an external impulse `force * dt` to a single particle.
    ///
    /// Modifying `prev_position` injects velocity into the Verlet integrator.
    pub fn apply_external_force(&mut self, particle_idx: usize, force: [f64; 3], dt: f64) {
        let p = &mut self.particles[particle_idx];
        if p.fixed {
            return;
        }
        let dv = scale3(force, dt / p.mass);
        p.prev_position = sub3(p.prev_position, dv);
    }
}
/// Solve volume constraints with internal pressure feedback.
///
/// Like `solve_volume_constraints`, but additionally applies a pressure
/// force proportional to the volume deficit, pushing particles outward
/// when compressed.
impl TetrahedralSoftBody {
    /// Project volume constraints with an additional internal pressure correction.
    #[allow(dead_code)]
    pub fn solve_volume_with_pressure(&mut self, pressure: f64) {
        for (idx, &tet) in self.tetrahedra.iter().enumerate() {
            let v_rest = self.rest_volumes[idx];
            let p0 = self.particles[tet[0]].position;
            let p1 = self.particles[tet[1]].position;
            let p2 = self.particles[tet[2]].position;
            let p3 = self.particles[tet[3]].position;
            let v_cur = tet_volume_signed(p0, p1, p2, p3).abs();
            let c = v_cur - v_rest;
            let grads = tet_volume_gradient(p0, p1, p2, p3);
            let w = [
                if self.particles[tet[0]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[0]].mass
                },
                if self.particles[tet[1]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[1]].mass
                },
                if self.particles[tet[2]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[2]].mass
                },
                if self.particles[tet[3]].fixed {
                    0.0
                } else {
                    1.0 / self.particles[tet[3]].mass
                },
            ];
            let lambda = volume_constraint_lambda(&grads, &w, c, self.volume_stiffness);
            let pressure_correction = if v_cur < v_rest {
                pressure * (v_rest - v_cur) / v_rest
            } else {
                0.0
            };
            for (k, &ti) in tet.iter().enumerate() {
                if !self.particles[ti].fixed {
                    let dp = scale3(grads[k], (lambda + pressure_correction) * w[k]);
                    self.particles[ti].position = add3(self.particles[ti].position, dp);
                }
            }
        }
    }
    /// Compute the total volume using the centroid decomposition method.
    ///
    /// For each tet, computes the signed volume of 4 sub-tets formed with
    /// the centroid, then verifies consistency.
    #[allow(dead_code)]
    pub fn total_volume_centroid(&self) -> f64 {
        let mut total = 0.0_f64;
        for tet in &self.tetrahedra {
            let positions: [[f64; 3]; 4] = [
                self.particles[tet[0]].position,
                self.particles[tet[1]].position,
                self.particles[tet[2]].position,
                self.particles[tet[3]].position,
            ];
            let c = tet_centroid(&positions);
            let v0 = tet_volume_signed(c, positions[0], positions[1], positions[2]).abs();
            let v1 = tet_volume_signed(c, positions[0], positions[1], positions[3]).abs();
            let v2 = tet_volume_signed(c, positions[0], positions[2], positions[3]).abs();
            let v3 = tet_volume_signed(c, positions[1], positions[2], positions[3]).abs();
            total += v0 + v1 + v2 + v3;
        }
        total
    }
    /// Compute volume change ratio: V_current / V_rest for each tetrahedron.
    #[allow(dead_code)]
    pub fn volume_ratios(&self) -> Vec<f64> {
        self.tetrahedra
            .iter()
            .enumerate()
            .map(|(idx, tet)| {
                let v = Self::tet_volume(&self.particles, tet);
                let v_rest = self.rest_volumes[idx];
                if v_rest.abs() < 1e-15 {
                    1.0
                } else {
                    v / v_rest
                }
            })
            .collect()
    }
    /// Step with pressure-volume work.
    ///
    /// Uses the ideal gas law to compute internal pressure based on
    /// volume change, then applies it during constraint projection.
    #[allow(dead_code)]
    pub fn step_with_gas_law(&mut self, dt: f64, gravity: [f64; 3], reference_pressure: f64) {
        let dt2 = dt * dt;
        for p in &mut self.particles {
            if p.fixed {
                continue;
            }
            let vel = sub3(p.position, p.prev_position);
            let new_pos = add3(add3(p.position, vel), scale3(gravity, dt2));
            p.prev_position = p.position;
            p.position = new_pos;
        }
        let total_rest: f64 = self.rest_volumes.iter().sum();
        let total_cur = self.total_volume();
        let pressure = ideal_gas_pressure(reference_pressure, total_rest, total_cur);
        self.solve_volume_with_pressure(pressure);
        self.solve_edge_constraints();
    }
}
impl TetrahedralSoftBody {
    /// Apply volume preservation penalty forces via position correction.
    /// Moves particles to reduce volume error.
    #[allow(dead_code)]
    pub fn apply_volume_preservation(&mut self, rest_volumes: &[f64], stiffness: f64) {
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.position).collect();
        let grad =
            volume_preservation_gradient(&self.tetrahedra, rest_volumes, &positions, stiffness);
        for (i, p) in self.particles.iter_mut().enumerate() {
            if !p.fixed {
                let inv_m = 1.0 / p.mass;
                let scale = inv_m * 0.0001;
                p.position[0] -= grad[i][0] * scale;
                p.position[1] -= grad[i][1] * scale;
                p.position[2] -= grad[i][2] * scale;
            }
        }
    }
    /// Return the current volumes for each tetrahedron.
    #[allow(dead_code)]
    pub fn compute_current_volumes(&self) -> Vec<f64> {
        self.tetrahedra
            .iter()
            .map(|tet| {
                tet_volume_signed(
                    self.particles[tet[0]].position,
                    self.particles[tet[1]].position,
                    self.particles[tet[2]].position,
                    self.particles[tet[3]].position,
                )
                .abs()
            })
            .collect()
    }
    /// Count the number of free (non-fixed) particles.
    #[allow(dead_code)]
    pub fn free_particle_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.fixed).count()
    }
    /// Estimate kinetic energy using Verlet velocity approximation: v ≈ (pos - prev_pos) / dt.
    /// Uses dt = 1.0 as a unit time.
    #[allow(dead_code)]
    pub fn approx_kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                if p.fixed {
                    0.0
                } else {
                    let vel = sub3(p.position, p.prev_position);
                    let v2 = dot3(vel, vel);
                    0.5 * p.mass * v2
                }
            })
            .sum()
    }
    /// Apply plane collision to all particles (positional correction only).
    #[allow(dead_code)]
    pub fn apply_plane_collision(&mut self, plane: &CollisionPlane) {
        let n = plane.normal;
        let p0 = plane.point;
        for p in &mut self.particles {
            if p.fixed {
                continue;
            }
            let diff = sub3(p.position, p0);
            let dist = dot3(diff, n);
            if dist < 0.0 {
                let depth = -dist;
                p.position[0] += n[0] * depth;
                p.position[1] += n[1] * depth;
                p.position[2] += n[2] * depth;
                let vel = sub3(p.position, p.prev_position);
                let vn = dot3(vel, n);
                if vn < 0.0 {
                    p.prev_position[0] -= n[0] * vn * 1.5;
                    p.prev_position[1] -= n[1] * vn * 1.5;
                    p.prev_position[2] -= n[2] * vn * 1.5;
                }
            }
        }
    }
    /// Simulate one step with Boyle's law gas pressure applied to a surface mesh.
    /// The gas pressure modifies prev_position to inject velocity into Verlet.
    #[allow(dead_code)]
    pub fn step_with_pressure(
        &mut self,
        triangles: &[[usize; 3]],
        rest_pressure: f64,
        rest_volume: f64,
        gravity: [f64; 3],
        dt: f64,
    ) {
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.position).collect();
        let current_volume = signed_surface_volume(triangles, &positions).abs();
        let pressure = boyle_pressure(rest_pressure, rest_volume, current_volume);
        for tri in triangles {
            let v0 = positions[tri[0]];
            let v1 = positions[tri[1]];
            let v2 = positions[tri[2]];
            let e01 = sub3(v1, v0);
            let e02 = sub3(v2, v0);
            let area_normal = cross3(e01, e02);
            let area = len3(area_normal) * 0.5;
            let normal = if len3(area_normal) > 1e-30 {
                normalize3(area_normal)
            } else {
                [0.0, 0.0, 0.0]
            };
            let force_per_particle = pressure * area / 3.0;
            for &vi in tri {
                let p = &mut self.particles[vi];
                if !p.fixed {
                    let accel = scale3(normal, force_per_particle / p.mass);
                    p.prev_position[0] -= accel[0] * dt * dt;
                    p.prev_position[1] -= accel[1] * dt * dt;
                    p.prev_position[2] -= accel[2] * dt * dt;
                }
            }
        }
        self.step(dt, gravity);
    }
}
/// Sphere-particle soft collision: push particles out of a sphere.
#[allow(dead_code)]
pub struct CollisionSphere {
    /// Centre of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}
/// Plane collision represented by a point on the plane and an outward normal.
#[allow(dead_code)]
pub struct CollisionPlane {
    /// A point on the plane.
    pub point: [f64; 3],
    /// Unit outward normal.
    pub normal: [f64; 3],
}
