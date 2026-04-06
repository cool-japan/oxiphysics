//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    add, apply_wrinkling, compute_dihedral_angle, compute_lift_force, cross, dot, length,
    normalize, resolve_cloth_floor_collision, resolve_cloth_sphere_collision,
    resolve_self_collision, scale, sub,
};
use crate::constraint::{BendingConstraint, DistanceConstraint};
use crate::particle::{SoftBody, SoftParticle};
use oxiphysics_core::math::{Real, Vec3};

/// A collection of seam constraints acting on a flat position array.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct SeamSystem {
    /// All seam constraints.
    pub seams: Vec<SeamConstraint>,
}
#[allow(dead_code)]
impl SeamSystem {
    /// Create an empty seam system.
    pub fn new() -> Self {
        Self { seams: Vec::new() }
    }
    /// Add a seam.
    pub fn add_seam(&mut self, seam: SeamConstraint) {
        self.seams.push(seam);
    }
    /// Apply one PBD iteration of all seam constraints.
    ///
    /// `positions` and `inv_masses` are modified in place.
    pub fn solve(&mut self, positions: &mut Vec<[f64; 3]>, inv_masses: &[f64]) {
        let n = self.seams.len();
        for idx in 0..n {
            let (da, db) = {
                let seam = &mut self.seams[idx];
                seam.apply(positions, inv_masses)
            };
            let a = self.seams[idx].vertex_a;
            let b = self.seams[idx].vertex_b;
            if inv_masses[a] > 0.0 {
                positions[a] = add(positions[a], da);
            }
            if inv_masses[b] > 0.0 {
                positions[b] = add(positions[b], db);
            }
        }
    }
    /// Number of active (non-torn) seams.
    pub fn active_count(&self) -> usize {
        self.seams.iter().filter(|s| !s.torn).count()
    }
    /// Number of torn seams.
    pub fn torn_count(&self) -> usize {
        self.seams.iter().filter(|s| s.torn).count()
    }
}
/// A cloth mesh created from a regular grid, integrated with the XPBD solver.
#[derive(Debug, Clone)]
pub struct XpbdClothMesh {
    /// Soft body containing the cloth particles.
    pub body: SoftBody,
    /// Distance constraints (structural + shear).
    pub distance_constraints: Vec<DistanceConstraint>,
    /// Bending constraints (across two triangles sharing an edge).
    pub bending_constraints: Vec<BendingConstraint>,
    /// Number of vertices along the X axis.
    pub nx: usize,
    /// Number of vertices along the Y axis.
    pub ny: usize,
}
impl XpbdClothMesh {
    /// Create a new cloth mesh lying in the XZ plane.
    ///
    /// * `nx`, `ny` - number of vertices in each dimension (must be >= 2).
    /// * `width`, `height` - physical dimensions of the cloth.
    /// * `mass_per_particle` - mass assigned to each particle.
    /// * `compliance` - XPBD compliance for the generated constraints.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        width: Real,
        height: Real,
        mass_per_particle: Real,
        compliance: Real,
    ) -> Self {
        assert!(nx >= 2 && ny >= 2, "Cloth grid must be at least 2x2");
        let dx = width / (nx - 1) as Real;
        let dy = height / (ny - 1) as Real;
        let mut particles = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let pos = Vec3::new(i as Real * dx, 0.0, j as Real * dy);
                particles.push(SoftParticle::new(pos, mass_per_particle));
            }
        }
        let idx = |ix: usize, iy: usize| -> usize { iy * nx + ix };
        let mut triangles: Vec<[usize; 3]> = Vec::new();
        for j in 0..(ny - 1) {
            for i in 0..(nx - 1) {
                let a = idx(i, j);
                let b = idx(i + 1, j);
                let c = idx(i, j + 1);
                let d = idx(i + 1, j + 1);
                triangles.push([a, b, c]);
                triangles.push([b, d, c]);
            }
        }
        let mut body = SoftBody::from_particles(particles);
        body.triangles = triangles;
        let mut distance_constraints = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                if i + 1 < nx {
                    distance_constraints.push(DistanceConstraint::from_particles(
                        idx(i, j),
                        idx(i + 1, j),
                        &body.particles,
                        compliance,
                    ));
                }
                if j + 1 < ny {
                    distance_constraints.push(DistanceConstraint::from_particles(
                        idx(i, j),
                        idx(i, j + 1),
                        &body.particles,
                        compliance,
                    ));
                }
                if i + 1 < nx && j + 1 < ny {
                    distance_constraints.push(DistanceConstraint::from_particles(
                        idx(i, j),
                        idx(i + 1, j + 1),
                        &body.particles,
                        compliance,
                    ));
                    distance_constraints.push(DistanceConstraint::from_particles(
                        idx(i + 1, j),
                        idx(i, j + 1),
                        &body.particles,
                        compliance,
                    ));
                }
            }
        }
        let mut bending_constraints = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                if i + 2 < nx && j + 1 < ny {
                    bending_constraints.push(BendingConstraint::from_particles(
                        [idx(i, j), idx(i + 2, j), idx(i + 1, j), idx(i + 1, j + 1)],
                        &body.particles,
                        compliance * 10.0,
                    ));
                }
                if j + 2 < ny && i + 1 < nx {
                    bending_constraints.push(BendingConstraint::from_particles(
                        [idx(i, j), idx(i, j + 2), idx(i, j + 1), idx(i + 1, j + 1)],
                        &body.particles,
                        compliance * 10.0,
                    ));
                }
            }
        }
        Self {
            body,
            distance_constraints,
            bending_constraints,
            nx,
            ny,
        }
    }
    /// Total number of particles.
    pub fn num_particles(&self) -> usize {
        self.nx * self.ny
    }
    /// Total number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.body.triangles.len()
    }
    /// Pin (fix) the particle at grid position (`ix`, `iy`).
    pub fn pin(&mut self, ix: usize, iy: usize) {
        let idx = iy * self.nx + ix;
        self.body.particles[idx].inverse_mass = 0.0;
    }
}
/// A spring edge connecting two cloth vertices.
#[derive(Debug, Clone)]
pub struct ClothEdge {
    /// First vertex index.
    pub a: usize,
    /// Second vertex index.
    pub b: usize,
    /// Rest (natural) length of the spring.
    pub rest_length: f64,
    /// Spring stiffness.
    pub stiffness: f64,
    /// Whether this edge has been torn.
    pub torn: bool,
    /// Stretch ratio (|d - rest| / rest) before tearing.
    pub tear_threshold: f64,
}

impl ClothEdge {
    /// Create a new intact edge with default tear threshold.
    pub fn new(a: usize, b: usize, rest_length: f64, stiffness: f64) -> Self {
        Self {
            a,
            b,
            rest_length,
            stiffness,
            torn: false,
            tear_threshold: 1.0,
        }
    }

    /// Compute the strain of this edge: (current_length - rest_length) / rest_length.
    pub fn strain(&self, vertices: &[ClothVertex]) -> f64 {
        let pa = vertices[self.a].pos;
        let pb = vertices[self.b].pos;
        let dx = pb[0] - pa[0];
        let dy = pb[1] - pa[1];
        let dz = pb[2] - pa[2];
        let current_len = (dx * dx + dy * dy + dz * dz).sqrt();
        if self.rest_length.abs() < f64::EPSILON {
            return 0.0;
        }
        (current_len - self.rest_length) / self.rest_length
    }
}

/// Standalone cloth simulation mesh (XPBD integration + tearing + wind).
///
/// Uses raw `[f64; 3]` arrays — no external math crate required.
pub struct ClothMesh {
    /// All vertices of the cloth.
    pub vertices: Vec<ClothVertex>,
    /// All spring edges.
    pub edges: Vec<ClothEdge>,
    /// All triangular faces.
    pub triangles: Vec<ClothTriangle>,
    /// Gravity acceleration vector.
    pub gravity: [f64; 3],
    /// Air resistance damping coefficient.
    pub air_resistance: f64,
    /// Wind velocity vector.
    pub wind: [f64; 3],
}
impl ClothMesh {
    /// Create a `nx × ny` grid cloth lying in the XZ plane.
    ///
    /// Generates structural, shear, and bending edges.
    pub fn new_grid(nx: usize, ny: usize, width: f64, height: f64, mass: f64) -> Self {
        assert!(nx >= 2 && ny >= 2, "Grid must be at least 2x2");
        let dx = width / (nx - 1) as f64;
        let dz = height / (ny - 1) as f64;
        let mut vertices = Vec::with_capacity(nx * ny);
        for r in 0..ny {
            for c in 0..nx {
                let pos = [c as f64 * dx, 0.0, r as f64 * dz];
                vertices.push(ClothVertex {
                    pos,
                    prev_pos: pos,
                    vel: [0.0, 0.0, 0.0],
                    mass,
                    inv_mass: if mass > 0.0 { 1.0 / mass } else { 0.0 },
                    fixed: false,
                    normal: [0.0, 1.0, 0.0],
                });
            }
        }
        let vert_idx = |r: usize, c: usize| r * nx + c;
        let tri_area = |p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]| -> f64 {
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            length(cross(e1, e2)) * 0.5
        };
        let mut edges = Vec::new();
        let mut triangles = Vec::new();
        for r in 0..ny {
            for c in 0..nx {
                if c + 1 < nx {
                    let a = vert_idx(r, c);
                    let b = vert_idx(r, c + 1);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 1.0,
                        torn: false,
                        tear_threshold: 0.5,
                    });
                }
                if r + 1 < ny {
                    let a = vert_idx(r, c);
                    let b = vert_idx(r + 1, c);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 1.0,
                        torn: false,
                        tear_threshold: 0.5,
                    });
                }
                if c + 1 < nx && r + 1 < ny {
                    let a = vert_idx(r, c);
                    let b = vert_idx(r + 1, c + 1);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 0.8,
                        torn: false,
                        tear_threshold: 0.6,
                    });
                    let a = vert_idx(r, c + 1);
                    let b = vert_idx(r + 1, c);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 0.8,
                        torn: false,
                        tear_threshold: 0.6,
                    });
                }
                if c + 2 < nx {
                    let a = vert_idx(r, c);
                    let b = vert_idx(r, c + 2);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 0.3,
                        torn: false,
                        tear_threshold: 1.0,
                    });
                }
                if r + 2 < ny {
                    let a = vert_idx(r, c);
                    let b = vert_idx(r + 2, c);
                    let rl = length(sub(vertices[b].pos, vertices[a].pos));
                    edges.push(ClothEdge {
                        a,
                        b,
                        rest_length: rl,
                        stiffness: 0.3,
                        torn: false,
                        tear_threshold: 1.0,
                    });
                }
            }
        }
        for r in 0..(ny - 1) {
            for c in 0..(nx - 1) {
                let a = vert_idx(r, c);
                let b = vert_idx(r, c + 1);
                let cc = vert_idx(r + 1, c);
                let d = vert_idx(r + 1, c + 1);
                let ra = tri_area(vertices[a].pos, vertices[b].pos, vertices[cc].pos);
                triangles.push(ClothTriangle {
                    indices: [a, b, cc],
                    rest_area: ra,
                });
                let ra2 = tri_area(vertices[b].pos, vertices[d].pos, vertices[cc].pos);
                triangles.push(ClothTriangle {
                    indices: [b, d, cc],
                    rest_area: ra2,
                });
            }
        }
        Self {
            vertices,
            edges,
            triangles,
            gravity: [0.0, -9.81, 0.0],
            air_resistance: 0.1,
            wind: [0.0, 0.0, 0.0],
        }
    }
    /// Pin one of the 4 corners (0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right).
    pub fn pin_corner(&mut self, corner: usize) {
        let n = self.vertices.len();
        let mut nx = n;
        for i in 1..n {
            if self.vertices[i].pos[2] > 1e-12 {
                nx = i;
                break;
            }
        }
        let ny = n / nx;
        let idx = match corner {
            0 => 0,
            1 => nx - 1,
            2 => (ny - 1) * nx,
            3 => ny * nx - 1,
            _ => panic!("corner must be 0..=3"),
        };
        self.vertices[idx].fixed = true;
        self.vertices[idx].inv_mass = 0.0;
    }
    /// Fix all vertices on the top row (row 0).
    pub fn pin_top_edge(&mut self) {
        let n = self.vertices.len();
        let mut nx = n;
        for i in 1..n {
            if self.vertices[i].pos[2] > 1e-12 {
                nx = i;
                break;
            }
        }
        for c in 0..nx {
            self.vertices[c].fixed = true;
            self.vertices[c].inv_mass = 0.0;
        }
    }
    /// Apply gravity to all free vertices: vel += gravity * dt.
    pub fn apply_gravity(&mut self, dt: f64) {
        for v in self.vertices.iter_mut() {
            if !v.fixed {
                v.vel[0] += self.gravity[0] * dt;
                v.vel[1] += self.gravity[1] * dt;
                v.vel[2] += self.gravity[2] * dt;
            }
        }
    }
    /// Apply wind force to each triangle, distributing to 3 vertices.
    pub fn apply_wind(&mut self, dt: f64) {
        let air_density = 1.225_f64;
        let tris: Vec<[usize; 3]> = self.triangles.iter().map(|t| t.indices).collect();
        for tri in &tris {
            let p0 = self.vertices[tri[0]].pos;
            let p1 = self.vertices[tri[1]].pos;
            let p2 = self.vertices[tri[2]].pos;
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            let area_normal = cross(e1, e2);
            let area2 = length(area_normal);
            if area2 < 1e-12 {
                continue;
            }
            let n_hat = scale(area_normal, 1.0 / area2);
            let area = area2 * 0.5;
            let v0 = self.vertices[tri[0]].vel;
            let v1 = self.vertices[tri[1]].vel;
            let v2 = self.vertices[tri[2]].vel;
            let vel_tri = scale(add(add(v0, v1), v2), 1.0 / 3.0);
            let w = self.wind;
            let rel = sub(w, vel_tri);
            let dot_val = dot(n_hat, rel);
            let force_mag = air_density * dot_val * area / 3.0;
            let force = scale(n_hat, force_mag);
            for &vi in tri.iter() {
                if !self.vertices[vi].fixed {
                    let f = scale(force, dt / self.vertices[vi].mass.max(1e-12));
                    self.vertices[vi].vel = add(self.vertices[vi].vel, f);
                }
            }
        }
    }
    /// Apply air resistance damping: vel *= exp(-air_resistance * dt).
    pub fn apply_air_resistance(&mut self, dt: f64) {
        let factor = (-self.air_resistance * dt).exp();
        for v in self.vertices.iter_mut() {
            if !v.fixed {
                v.vel = scale(v.vel, factor);
            }
        }
    }
    /// Predict positions: pos_new = pos + vel * dt.
    pub fn predict_positions(&mut self, dt: f64) {
        for v in self.vertices.iter_mut() {
            if !v.fixed {
                v.prev_pos = v.pos;
                v.pos = add(v.pos, scale(v.vel, dt));
            }
        }
    }
    /// Update velocities from position change: vel = (pos - prev_pos) / dt.
    pub fn update_velocities(&mut self, dt: f64) {
        for v in self.vertices.iter_mut() {
            if !v.fixed {
                v.vel = scale(sub(v.pos, v.prev_pos), 1.0 / dt);
            }
        }
    }
    /// XPBD distance constraint solving for each non-torn edge.
    pub fn solve_stretch_constraints(&mut self, dt: f64) {
        let compliance = 1e-4_f64;
        let alpha = compliance / (dt * dt);
        let edge_data: Vec<(usize, usize, f64, bool)> = self
            .edges
            .iter()
            .map(|e| (e.a, e.b, e.rest_length, e.torn))
            .collect();
        for (a, b, rest, torn) in edge_data {
            if torn {
                continue;
            }
            let pa = self.vertices[a].pos;
            let pb = self.vertices[b].pos;
            let wa = self.vertices[a].inv_mass;
            let wb = self.vertices[b].inv_mass;
            let d_vec = sub(pb, pa);
            let d = length(d_vec);
            if d < 1e-12 {
                continue;
            }
            let dir = scale(d_vec, 1.0 / d);
            let c = d - rest;
            let denom = wa + wb + alpha;
            if denom.abs() < 1e-30 {
                continue;
            }
            let lambda = -c / denom;
            if !self.vertices[a].fixed {
                self.vertices[a].pos = add(pa, scale(dir, -lambda * wa));
            }
            if !self.vertices[b].fixed {
                self.vertices[b].pos = add(pb, scale(dir, lambda * wb));
            }
        }
    }
    /// Simple bending constraint: for bending edges (stiffness < 0.5),
    /// softly maintain rest length.
    pub fn solve_bending_constraints(&mut self) {
        let compliance = 1e-2_f64;
        let edge_data: Vec<(usize, usize, f64, f64, bool)> = self
            .edges
            .iter()
            .map(|e| (e.a, e.b, e.rest_length, e.stiffness, e.torn))
            .collect();
        for (a, b, rest, stiffness, torn) in edge_data {
            if stiffness >= 0.5 || torn {
                continue;
            }
            let pa = self.vertices[a].pos;
            let pb = self.vertices[b].pos;
            let wa = self.vertices[a].inv_mass;
            let wb = self.vertices[b].inv_mass;
            let d_vec = sub(pb, pa);
            let d = length(d_vec);
            if d < 1e-12 {
                continue;
            }
            let dir = scale(d_vec, 1.0 / d);
            let c = d - rest;
            let denom = wa + wb + compliance;
            if denom.abs() < 1e-30 {
                continue;
            }
            let lambda = -c / denom;
            if !self.vertices[a].fixed {
                self.vertices[a].pos = add(pa, scale(dir, -lambda * wa));
            }
            if !self.vertices[b].fixed {
                self.vertices[b].pos = add(pb, scale(dir, lambda * wb));
            }
        }
    }
    /// Mark edges as torn if stretch ratio exceeds tear_threshold.
    pub fn check_tearing(&mut self) {
        for edge in self.edges.iter_mut() {
            if edge.torn {
                continue;
            }
            let pa = self.vertices[edge.a].pos;
            let pb = self.vertices[edge.b].pos;
            let d = length(sub(pb, pa));
            let stretch = (d - edge.rest_length).abs() / edge.rest_length.max(1e-12);
            if stretch > edge.tear_threshold {
                edge.torn = true;
            }
        }
    }
    /// Run one simulation step: predict → solve stretch (×iterations) →
    /// bending → check tearing → update velocities.
    pub fn step(&mut self, dt: f64, iterations: usize) {
        self.apply_gravity(dt);
        self.apply_wind(dt);
        self.apply_air_resistance(dt);
        self.predict_positions(dt);
        for _ in 0..iterations {
            self.solve_stretch_constraints(dt);
        }
        self.solve_bending_constraints();
        self.check_tearing();
        self.update_velocities(dt);
    }
    /// Compute area-weighted per-vertex normals from triangle face normals.
    pub fn compute_vertex_normals(&mut self) {
        for v in self.vertices.iter_mut() {
            v.normal = [0.0, 0.0, 0.0];
        }
        let tris: Vec<[usize; 3]> = self.triangles.iter().map(|t| t.indices).collect();
        for tri in &tris {
            let p0 = self.vertices[tri[0]].pos;
            let p1 = self.vertices[tri[1]].pos;
            let p2 = self.vertices[tri[2]].pos;
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            let area_normal = cross(e1, e2);
            for &vi in tri.iter() {
                self.vertices[vi].normal = add(self.vertices[vi].normal, area_normal);
            }
        }
        for v in self.vertices.iter_mut() {
            v.normal = normalize(v.normal);
        }
    }
    /// Return the number of torn edges.
    pub fn total_torn_edges(&self) -> usize {
        self.edges.iter().filter(|e| e.torn).count()
    }
    /// Return the total kinetic energy: sum(0.5 * m * |v|^2).
    pub fn kinetic_energy(&self) -> f64 {
        self.vertices
            .iter()
            .map(|v| 0.5 * v.mass * dot(v.vel, v.vel))
            .sum()
    }
    /// Return the total elastic (spring) potential energy: sum(0.5 * k * (d - rest)^2).
    pub fn potential_energy(&self) -> f64 {
        self.edges
            .iter()
            .filter(|e| !e.torn)
            .map(|e| {
                let pa = self.vertices[e.a].pos;
                let pb = self.vertices[e.b].pos;
                let d = length(sub(pb, pa));
                let stretch = d - e.rest_length;
                0.5 * e.stiffness * stretch * stretch
            })
            .sum()
    }
}
impl ClothMesh {
    /// Resolve pairwise self-collision between all free vertices.
    ///
    /// Pairs closer than `min_dist` are pushed apart using XPBD-style
    /// positional correction.
    pub fn resolve_self_collision(&mut self, min_dist: f64) {
        resolve_self_collision(&mut self.vertices, min_dist);
    }
    /// Resolve collision of the cloth against a sphere obstacle.
    pub fn collide_with_sphere(&mut self, sphere: &RigidSphere) {
        resolve_cloth_sphere_collision(&mut self.vertices, sphere);
    }
    /// Resolve collision of the cloth against an infinite floor at `floor_y`.
    pub fn collide_with_floor(&mut self, floor_y: f64, restitution: f64) {
        resolve_cloth_floor_collision(&mut self.vertices, floor_y, restitution);
    }
    /// Apply aerodynamic lift in addition to drag.
    ///
    /// `lift_coeff` is the flat-plate lift coefficient.
    pub fn apply_lift(&mut self, dt: f64, lift_coeff: f64) {
        let air_density = 1.225_f64;
        let tris: Vec<[usize; 3]> = self.triangles.iter().map(|t| t.indices).collect();
        let wind = self.wind;
        for tri in &tris {
            let p0 = self.vertices[tri[0]].pos;
            let p1 = self.vertices[tri[1]].pos;
            let p2 = self.vertices[tri[2]].pos;
            let v0 = self.vertices[tri[0]].vel;
            let v1 = self.vertices[tri[1]].vel;
            let v2 = self.vertices[tri[2]].vel;
            let v_tri = scale(add(add(v0, v1), v2), 1.0 / 3.0);
            let lift = compute_lift_force(p0, p1, p2, v_tri, wind, air_density, lift_coeff);
            let f_per_vert = scale(lift, 1.0 / 3.0);
            for &vi in tri.iter() {
                if !self.vertices[vi].fixed {
                    let m = self.vertices[vi].mass.max(1e-12);
                    self.vertices[vi].vel = add(self.vertices[vi].vel, scale(f_per_vert, dt / m));
                }
            }
        }
    }
    /// Displace each vertex along its normal to simulate cloth wrinkling.
    ///
    /// `wrinkle_scale` controls the amplitude of the wrinkle displacement.
    pub fn apply_wrinkling(&mut self, wrinkle_scale: f64) {
        self.compute_vertex_normals();
        apply_wrinkling(&mut self.vertices, &self.edges, wrinkle_scale);
    }
    /// Return indices of all torn edges.
    pub fn torn_edge_indices(&self) -> Vec<usize> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.torn)
            .map(|(i, _)| i)
            .collect()
    }
    /// Repair (un-tear) all edges whose stretch is back below the repair ratio.
    ///
    /// An edge is repaired if `|d - rest| / rest < repair_ratio`.
    pub fn repair_tears(&mut self, repair_ratio: f64) {
        for edge in &mut self.edges {
            if !edge.torn {
                continue;
            }
            let pa = self.vertices[edge.a].pos;
            let pb = self.vertices[edge.b].pos;
            let d = length(sub(pb, pa));
            let stretch = (d - edge.rest_length).abs() / edge.rest_length.max(1e-12);
            if stretch < repair_ratio {
                edge.torn = false;
            }
        }
    }
    /// Force-tear all edges in the given index list.
    pub fn force_tear_edges(&mut self, indices: &[usize]) {
        for &i in indices {
            if i < self.edges.len() {
                self.edges[i].torn = true;
            }
        }
    }
    /// Full simulation step with optional self-collision and sphere obstacle.
    ///
    /// `sphere` may be `None` if no obstacle is present.
    pub fn step_with_collision(
        &mut self,
        dt: f64,
        iterations: usize,
        self_collision_dist: Option<f64>,
        sphere: Option<&RigidSphere>,
        floor_y: Option<f64>,
    ) {
        self.step(dt, iterations);
        if let Some(d) = self_collision_dist {
            self.resolve_self_collision(d);
        }
        if let Some(s) = sphere {
            self.collide_with_sphere(s);
        }
        if let Some(fy) = floor_y {
            self.collide_with_floor(fy, 0.3);
        }
    }
    /// Return the average triangle area of the cloth.
    pub fn average_triangle_area(&self) -> f64 {
        if self.triangles.is_empty() {
            return 0.0;
        }
        let total: f64 = self
            .triangles
            .iter()
            .map(|t| {
                let p0 = self.vertices[t.indices[0]].pos;
                let p1 = self.vertices[t.indices[1]].pos;
                let p2 = self.vertices[t.indices[2]].pos;
                length(cross(sub(p1, p0), sub(p2, p0))) * 0.5
            })
            .sum();
        total / self.triangles.len() as f64
    }
    /// Return the total area of the cloth (sum of all triangle areas).
    pub fn total_area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let p0 = self.vertices[t.indices[0]].pos;
                let p1 = self.vertices[t.indices[1]].pos;
                let p2 = self.vertices[t.indices[2]].pos;
                length(cross(sub(p1, p0), sub(p2, p0))) * 0.5
            })
            .sum()
    }
    /// Return the maximum velocity magnitude across all vertices.
    pub fn max_vertex_speed(&self) -> f64 {
        self.vertices
            .iter()
            .map(|v| length(v.vel))
            .fold(0.0_f64, f64::max)
    }
    /// Return the axis-aligned bounding box `([xmin,ymin,zmin], [xmax,ymax,zmax])`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        if self.vertices.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = self.vertices[0].pos;
        let mut hi = self.vertices[0].pos;
        for v in &self.vertices[1..] {
            for k in 0..3 {
                if v.pos[k] < lo[k] {
                    lo[k] = v.pos[k];
                }
                if v.pos[k] > hi[k] {
                    hi[k] = v.pos[k];
                }
            }
        }
        (lo, hi)
    }
}
/// A constraint on a PBD cloth simulation.
#[derive(Debug, Clone)]
pub enum ClothConstraint {
    /// Stretch (distance) constraint between two particles.
    Stretch {
        /// Index of the first particle.
        i: usize,
        /// Index of the second particle.
        j: usize,
        /// Rest (natural) length.
        rest: f64,
        /// Stiffness in \[0, 1\].
        k: f64,
    },
    /// Dihedral bending constraint across four particles forming two triangles.
    Bend {
        /// First particle (shared edge endpoint 0).
        i: usize,
        /// Second particle (shared edge endpoint 1).
        j: usize,
        /// Third particle (wing of triangle A).
        k: usize,
        /// Fourth particle (wing of triangle B).
        l: usize,
        /// Rest dihedral angle (radians).
        rest_angle: f64,
        /// Bending stiffness in \[0, 1\].
        stiffness: f64,
    },
}
/// A single cloth vertex.
#[derive(Debug, Clone)]
pub struct ClothVertex {
    /// Current world-space position.
    pub pos: [f64; 3],
    /// Position at the previous time-step.
    pub prev_pos: [f64; 3],
    /// Current velocity.
    pub vel: [f64; 3],
    /// Vertex mass (kg).
    pub mass: f64,
    /// Inverse mass (0 for fixed/pinned vertices).
    pub inv_mass: f64,
    /// If `true` the vertex is pinned and never moves.
    pub fixed: bool,
    /// Per-vertex normal for wind interaction.
    pub normal: [f64; 3],
}
/// A rigid sphere obstacle.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RigidSphere {
    /// Centre of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere (m).
    pub radius: f64,
    /// Restitution coefficient in \[0, 1\].
    pub restitution: f64,
}
impl RigidSphere {
    /// Create a new rigid sphere.
    pub fn new(center: [f64; 3], radius: f64, restitution: f64) -> Self {
        Self {
            center,
            radius,
            restitution,
        }
    }
}
/// PBD cloth mesh using flat position/velocity/mass arrays and a constraint list.
///
/// Triangles are stored as row-major pairs from a grid:
/// for a grid of `nx × ny`, triangles are \[(r,c),(r,c+1),(r+1,c)\] and \[(r,c+1),(r+1,c+1),(r+1,c)\].
#[derive(Debug, Clone)]
pub struct PbdClothMesh {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Stretch and bend constraints.
    pub constraints: Vec<ClothConstraint>,
    /// Indices of pinned (kinematically fixed) particles.
    pub pinned: Vec<usize>,
    /// Triangle connectivity (used for wind/area).
    pub triangles: Vec<[usize; 3]>,
    /// Grid width (nx).
    pub nx: usize,
    /// Grid height (ny).
    pub ny: usize,
}
impl PbdClothMesh {
    /// Create a flat `nx × ny` grid lying in the XZ plane.
    ///
    /// Generates stretch constraints (horizontal, vertical, shear) and bend
    /// constraints (dihedral across shared edges).  No particles are pinned.
    pub fn grid(nx: usize, ny: usize, spacing: f64, mass: f64) -> Self {
        assert!(nx >= 2 && ny >= 2, "Grid must be at least 2×2");
        let n = nx * ny;
        let mut positions = Vec::with_capacity(n);
        for r in 0..ny {
            for c in 0..nx {
                positions.push([c as f64 * spacing, 0.0, r as f64 * spacing]);
            }
        }
        let velocities = vec![[0.0_f64; 3]; n];
        let masses = vec![mass; n];
        let vidx = |r: usize, c: usize| r * nx + c;
        let mut constraints = Vec::new();
        for r in 0..ny {
            for c in 0..nx - 1 {
                let a = vidx(r, c);
                let b = vidx(r, c + 1);
                let rest = length(sub(positions[b], positions[a]));
                constraints.push(ClothConstraint::Stretch {
                    i: a,
                    j: b,
                    rest,
                    k: 1.0,
                });
            }
        }
        for r in 0..ny - 1 {
            for c in 0..nx {
                let a = vidx(r, c);
                let b = vidx(r + 1, c);
                let rest = length(sub(positions[b], positions[a]));
                constraints.push(ClothConstraint::Stretch {
                    i: a,
                    j: b,
                    rest,
                    k: 1.0,
                });
            }
        }
        for r in 0..ny - 1 {
            for c in 0..nx - 1 {
                let a = vidx(r, c);
                let b = vidx(r + 1, c + 1);
                let rest_ab = length(sub(positions[b], positions[a]));
                constraints.push(ClothConstraint::Stretch {
                    i: a,
                    j: b,
                    rest: rest_ab,
                    k: 0.8,
                });
                let c2 = vidx(r, c + 1);
                let d = vidx(r + 1, c);
                let rest_cd = length(sub(positions[d], positions[c2]));
                constraints.push(ClothConstraint::Stretch {
                    i: c2,
                    j: d,
                    rest: rest_cd,
                    k: 0.8,
                });
            }
        }
        for r in 1..ny - 1 {
            for c in 0..nx - 1 {
                let i = vidx(r, c);
                let j = vidx(r, c + 1);
                let k_idx = vidx(r - 1, c);
                let l_idx = vidx(r + 1, c);
                let rest_angle = compute_dihedral_angle(
                    positions[i],
                    positions[j],
                    positions[k_idx],
                    positions[l_idx],
                );
                constraints.push(ClothConstraint::Bend {
                    i,
                    j,
                    k: k_idx,
                    l: l_idx,
                    rest_angle,
                    stiffness: 0.3,
                });
            }
        }
        for r in 0..ny - 1 {
            for c in 1..nx - 1 {
                let i = vidx(r, c);
                let j = vidx(r + 1, c);
                let k_idx = vidx(r, c - 1);
                let l_idx = vidx(r, c + 1);
                let rest_angle = compute_dihedral_angle(
                    positions[i],
                    positions[j],
                    positions[k_idx],
                    positions[l_idx],
                );
                constraints.push(ClothConstraint::Bend {
                    i,
                    j,
                    k: k_idx,
                    l: l_idx,
                    rest_angle,
                    stiffness: 0.3,
                });
            }
        }
        let mut triangles = Vec::with_capacity(2 * (nx - 1) * (ny - 1));
        for r in 0..ny - 1 {
            for c in 0..nx - 1 {
                let a = vidx(r, c);
                let b = vidx(r, c + 1);
                let cc = vidx(r + 1, c);
                let d = vidx(r + 1, c + 1);
                triangles.push([a, b, cc]);
                triangles.push([b, d, cc]);
            }
        }
        Self {
            positions,
            velocities,
            masses,
            constraints,
            pinned: Vec::new(),
            triangles,
            nx,
            ny,
        }
    }
    /// Pin a particle by index (it will not move).
    pub fn pin(&mut self, idx: usize) {
        if !self.pinned.contains(&idx) {
            self.pinned.push(idx);
        }
    }
    #[inline]
    fn is_pinned(&self, idx: usize) -> bool {
        self.pinned.contains(&idx)
    }
    /// One XPBD timestep under gravity with `iterations` constraint passes.
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, dt: f64, gravity: [f64; 3], iterations: usize) {
        let n = self.positions.len();
        let mut pred: Vec<[f64; 3]> = Vec::with_capacity(n);
        for i in 0..n {
            if self.is_pinned(i) {
                pred.push(self.positions[i]);
            } else {
                let v = add(self.velocities[i], scale(gravity, dt));
                pred.push(add(self.positions[i], scale(v, dt)));
            }
        }
        let compliance = 1e-4_f64;
        let alpha = compliance / (dt * dt);
        let inv_masses: Vec<f64> = self
            .masses
            .iter()
            .map(|&m| if m > 1e-30 { 1.0 / m } else { 0.0 })
            .collect();
        let constraints_snap = self.constraints.clone();
        for _ in 0..iterations {
            for c in &constraints_snap {
                match c {
                    ClothConstraint::Stretch { i, j, rest, k } => {
                        let pi = pred[*i];
                        let pj = pred[*j];
                        let d_vec = sub(pj, pi);
                        let d = length(d_vec);
                        if d < 1e-12 {
                            continue;
                        }
                        let dir = scale(d_vec, 1.0 / d);
                        let constraint = d - rest;
                        let wi = if self.is_pinned(*i) {
                            0.0
                        } else {
                            inv_masses[*i]
                        };
                        let wj = if self.is_pinned(*j) {
                            0.0
                        } else {
                            inv_masses[*j]
                        };
                        let denom = wi + wj + alpha / k;
                        if denom < 1e-30 {
                            continue;
                        }
                        let lam = -constraint / denom;
                        if !self.is_pinned(*i) {
                            pred[*i] = add(pi, scale(dir, -lam * wi));
                        }
                        if !self.is_pinned(*j) {
                            pred[*j] = add(pj, scale(dir, lam * wj));
                        }
                    }
                    ClothConstraint::Bend {
                        i,
                        j,
                        k: ki,
                        l,
                        rest_angle,
                        stiffness: stiff,
                    } => {
                        let pi = pred[*i];
                        let pj = pred[*j];
                        let pk = pred[*ki];
                        let pl = pred[*l];
                        let cur = compute_dihedral_angle(pi, pj, pk, pl);
                        let err = cur - rest_angle;
                        if err.abs() < 1e-10 {
                            continue;
                        }
                        let wi = if self.is_pinned(*ki) {
                            0.0
                        } else {
                            inv_masses[*ki]
                        };
                        let wl = if self.is_pinned(*l) {
                            0.0
                        } else {
                            inv_masses[*l]
                        };
                        let total_w = wi + wl;
                        if total_w < 1e-30 {
                            continue;
                        }
                        let mid = scale(add(pi, pj), 0.5);
                        let dk = normalize(sub(pk, mid));
                        let dl = normalize(sub(pl, mid));
                        let corr_mag = err * stiff * 0.01;
                        if !self.is_pinned(*ki) {
                            pred[*ki] = sub(pk, scale(dk, corr_mag * wi / total_w));
                        }
                        if !self.is_pinned(*l) {
                            pred[*l] = sub(pl, scale(dl, corr_mag * wl / total_w));
                        }
                    }
                }
            }
        }
        for i in 0..n {
            if self.is_pinned(i) {
                self.velocities[i] = [0.0; 3];
                continue;
            }
            self.velocities[i] = scale(sub(pred[i], self.positions[i]), 1.0 / dt);
            self.positions[i] = pred[i];
        }
    }
    /// Apply aerodynamic wind/drag to each particle using triangle normals.
    pub fn apply_wind(&mut self, wind: [f64; 3]) {
        let air_density = 1.225_f64;
        let dt = 1.0 / 60.0;
        let tris = self.triangles.clone();
        for tri in &tris {
            let p0 = self.positions[tri[0]];
            let p1 = self.positions[tri[1]];
            let p2 = self.positions[tri[2]];
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            let area_normal = cross(e1, e2);
            let area2 = length(area_normal);
            if area2 < 1e-12 {
                continue;
            }
            let n_hat = scale(area_normal, 1.0 / area2);
            let area = area2 * 0.5;
            let v0 = self.velocities[tri[0]];
            let v1 = self.velocities[tri[1]];
            let v2 = self.velocities[tri[2]];
            let vel_tri = scale(add(add(v0, v1), v2), 1.0 / 3.0);
            let rel = sub(wind, vel_tri);
            let dot_val = dot(n_hat, rel);
            let force_mag = air_density * dot_val * area / 3.0;
            let force = scale(n_hat, force_mag);
            for &vi in tri.iter() {
                if !self.is_pinned(vi) {
                    let m = self.masses[vi].max(1e-12);
                    self.velocities[vi] = add(self.velocities[vi], scale(force, dt / m));
                }
            }
        }
    }
    /// Push all particles inside sphere of `center`/`r` to its surface.
    pub fn apply_sphere_collision(&mut self, center: [f64; 3], r: f64) {
        for i in 0..self.positions.len() {
            if self.is_pinned(i) {
                continue;
            }
            let d = sub(self.positions[i], center);
            let dist = length(d);
            if dist < r && dist > 1e-15 {
                let n = normalize(d);
                self.positions[i] = add(center, scale(n, r));
                let vn = dot(self.velocities[i], n);
                if vn < 0.0 {
                    self.velocities[i] = sub(self.velocities[i], scale(n, vn));
                }
            }
        }
    }
    /// Sum of triangle areas (current geometry).
    pub fn total_area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|tri| {
                let p0 = self.positions[tri[0]];
                let p1 = self.positions[tri[1]];
                let p2 = self.positions[tri[2]];
                length(cross(sub(p1, p0), sub(p2, p0))) * 0.5
            })
            .sum()
    }
}
/// A seam constraint that connects vertex `a` on one cloth piece to vertex `b`
/// on another (or on the same mesh) at a user-specified rest distance.
///
/// When `rest_distance` is 0 the constraint acts as a hard weld; non-zero
/// values model an elastic seam.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SeamConstraint {
    /// Index of the first vertex (in the mesh vertex array).
    pub vertex_a: usize,
    /// Index of the second vertex.
    pub vertex_b: usize,
    /// Rest distance (m).  0 = weld.
    pub rest_distance: f64,
    /// Stiffness in \[0, 1\].
    pub stiffness: f64,
    /// Whether this seam has been broken (torn open).
    pub torn: bool,
    /// Strain threshold above which the seam tears.
    pub tear_threshold: f64,
}
#[allow(dead_code)]
impl SeamConstraint {
    /// Create a seam welding `a` to `b` at their current distance.
    pub fn new_weld(a: usize, b: usize, stiffness: f64) -> Self {
        Self {
            vertex_a: a,
            vertex_b: b,
            rest_distance: 0.0,
            stiffness,
            torn: false,
            tear_threshold: f64::MAX,
        }
    }
    /// Create a seam with an explicit rest distance.
    pub fn new_elastic(
        a: usize,
        b: usize,
        rest_distance: f64,
        stiffness: f64,
        tear_threshold: f64,
    ) -> Self {
        Self {
            vertex_a: a,
            vertex_b: b,
            rest_distance,
            stiffness,
            torn: false,
            tear_threshold,
        }
    }
    /// Current stretch ratio: `(d - rest) / rest` (or `d / ε` when rest ≈ 0).
    pub fn stretch_ratio(&self, positions: &[[f64; 3]]) -> f64 {
        let pa = positions[self.vertex_a];
        let pb = positions[self.vertex_b];
        let d = length(sub(pb, pa));
        if self.rest_distance < 1e-12 {
            d
        } else {
            (d - self.rest_distance) / self.rest_distance
        }
    }
    /// Apply the seam correction to predicted positions (PBD style).
    ///
    /// `inv_masses[i]` = 1/m_i (0 for pinned particles).
    /// Returns `(delta_a, delta_b)` position corrections.
    pub fn apply(&mut self, positions: &[[f64; 3]], inv_masses: &[f64]) -> ([f64; 3], [f64; 3]) {
        if self.torn {
            return ([0.0; 3], [0.0; 3]);
        }
        let pa = positions[self.vertex_a];
        let pb = positions[self.vertex_b];
        let diff = sub(pb, pa);
        let d = length(diff);
        let strain = if self.rest_distance < 1e-12 {
            d
        } else {
            (d - self.rest_distance) / self.rest_distance
        };
        if strain > self.tear_threshold {
            self.torn = true;
            return ([0.0; 3], [0.0; 3]);
        }
        if d < 1e-14 {
            return ([0.0; 3], [0.0; 3]);
        }
        let dir = scale(diff, 1.0 / d);
        let c = d - self.rest_distance;
        let wa = inv_masses[self.vertex_a];
        let wb = inv_masses[self.vertex_b];
        let denom = wa + wb;
        if denom < 1e-30 {
            return ([0.0; 3], [0.0; 3]);
        }
        let lambda = -c * self.stiffness / denom;
        let da = scale(dir, -lambda * wa);
        let db = scale(dir, lambda * wb);
        (da, db)
    }
}
/// A triangular face of the cloth mesh.
#[derive(Debug, Clone)]
pub struct ClothTriangle {
    /// Vertex indices of the triangle.
    pub indices: [usize; 3],
    /// Rest area of the triangle.
    pub rest_area: f64,
}
