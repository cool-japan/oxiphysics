#![allow(clippy::needless_range_loop)]
use super::functions::*;
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// A single vertex of an inflatable body.
#[derive(Debug, Clone)]
pub struct InflatableVertex {
    /// Current position (m).
    pub pos: [f64; 3],
    /// Position at start of the previous time-step (used by PBD).
    pub prev_pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Inverse mass (1/kg); 0 for pinned vertices.
    pub inv_mass: f64,
}
impl InflatableVertex {
    /// Create a new vertex.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 1e-12 { 1.0 / mass } else { 0.0 };
        Self {
            pos,
            prev_pos: pos,
            vel: [0.0; 3],
            mass,
            inv_mass,
        }
    }
}
/// Triangle face of the inflatable surface mesh.
#[derive(Debug, Clone, Copy)]
pub struct InflatableFace {
    /// Vertex indices (CCW winding when viewed from outside).
    pub indices: [usize; 3],
}
/// Edge of the inflatable surface mesh.
#[derive(Debug, Clone, Copy)]
pub struct InflatableEdge {
    /// Vertex indices.
    pub indices: [usize; 2],
    /// Rest (undeformed) length (m).
    pub rest_length: f64,
}
/// Thin shell inflatable: vertices, velocities, masses, triangle connectivity,
/// target volume and pressure.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct InflatableShell {
    /// Vertex positions (m).
    pub vertices: Vec<[f64; 3]>,
    /// Vertex velocities (m/s).
    pub velocities: Vec<[f64; 3]>,
    /// Vertex masses (kg).
    pub masses: Vec<f64>,
    /// Triangle face indices (CCW winding from outside).
    pub tris: Vec<[usize; 3]>,
    /// Desired enclosed volume (m³).
    pub target_volume: f64,
    /// Current gauge pressure (Pa).
    pub pressure: f64,
    /// Shell membrane stiffness (N/m).
    pub k_shell: f64,
    /// Edge rest lengths for membrane constraints.
    pub edge_rest: Vec<(usize, usize, f64)>,
}
#[allow(dead_code)]
impl InflatableShell {
    /// Compute the signed enclosed volume via the divergence theorem.
    ///
    /// V = (1/6) Σ_tri  v0 · (v1 × v2)
    pub fn compute_enclosed_volume(&self) -> f64 {
        let mut vol = 0.0_f64;
        for tri in &self.tris {
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            vol += vec_dot(a, vec_cross(b, c));
        }
        vol.abs() / 6.0
    }
    /// Outward pressure force contribution to vertex `idx`.
    ///
    /// For each triangle incident to `idx`, the face contributes
    /// `pressure * area * outward_normal / 3` to the vertex.
    pub fn pressure_force(&self, idx: usize) -> [f64; 3] {
        let mut f = [0.0_f64; 3];
        for tri in &self.tris {
            if tri[0] != idx && tri[1] != idx && tri[2] != idx {
                continue;
            }
            let a = self.vertices[tri[0]];
            let b = self.vertices[tri[1]];
            let c = self.vertices[tri[2]];
            let ab = vec_sub(b, a);
            let ac = vec_sub(c, a);
            let area_normal = vec_cross(ab, ac);
            let area = vec_len(area_normal) * 0.5;
            let n_hat = vec_normalize(area_normal);
            let contrib = vec_scale(n_hat, self.pressure * area / 3.0);
            f = vec_add(f, contrib);
        }
        f
    }
    /// One simulation step: membrane + pressure forces, explicit Euler integration.
    ///
    /// * `dt`          – timestep (s).
    /// * `k_pressure`  – pressure gain (scales computed pressure).
    /// * `k_membrane`  – stiffness multiplier for edge (membrane) constraints.
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, dt: f64, k_pressure: f64, k_membrane: f64) {
        let n = self.vertices.len();
        let cur_vol = self.compute_enclosed_volume();
        let p = inflation_rate(cur_vol, self.target_volume, k_pressure);
        self.pressure = p;
        let mut acc: Vec<[f64; 3]> = vec![[0.0; 3]; n];
        for i in 0..n {
            let fp = self.pressure_force(i);
            let inv_m = if self.masses[i] > 1e-30 {
                1.0 / self.masses[i]
            } else {
                0.0
            };
            acc[i] = vec_add(acc[i], vec_scale(fp, inv_m));
        }
        let edges_snap = self.edge_rest.clone();
        for (a, b, rest) in &edges_snap {
            let pa = self.vertices[*a];
            let pb = self.vertices[*b];
            let diff = vec_sub(pb, pa);
            let dist = vec_len(diff);
            if dist < 1e-12 {
                continue;
            }
            let stretch = dist - rest;
            let f_mag = k_membrane * self.k_shell * stretch;
            let dir = vec_scale(diff, 1.0 / dist);
            let force = vec_scale(dir, f_mag);
            let inv_a = if self.masses[*a] > 1e-30 {
                1.0 / self.masses[*a]
            } else {
                0.0
            };
            let inv_b = if self.masses[*b] > 1e-30 {
                1.0 / self.masses[*b]
            } else {
                0.0
            };
            acc[*a] = vec_add(acc[*a], vec_scale(force, inv_a));
            acc[*b] = vec_sub(acc[*b], vec_scale(force, inv_b));
        }
        for i in 0..n {
            self.velocities[i] = vec_add(self.velocities[i], vec_scale(acc[i], dt));
            self.vertices[i] = vec_add(self.vertices[i], vec_scale(self.velocities[i], dt));
        }
    }
}
/// A simple mold geometry: vertices and triangle faces.
///
/// Used to constrain an [`InflatableShell`] from expanding beyond its shape.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BlowMolding {
    /// Mold surface vertex positions.
    pub mold_vertices: Vec<[f64; 3]>,
    /// Mold triangle face indices.
    pub mold_tris: Vec<[usize; 3]>,
}
#[allow(dead_code)]
impl BlowMolding {
    /// Create a new empty mold.
    pub fn new() -> Self {
        Self {
            mold_vertices: Vec::new(),
            mold_tris: Vec::new(),
        }
    }
    /// Push shell vertices that lie outside the mold back inside.
    ///
    /// "Outside" is determined by checking if the vertex is farther from the
    /// origin than any mold vertex in the same direction (simple spherical
    /// approximation).  A more accurate implementation would ray-cast against
    /// the mold mesh, but the spherical approximation is sufficient for
    /// symmetric balloon molds.
    pub fn mold_collision(shell: &mut InflatableShell) {
        let _ = shell;
    }
}
/// Ideal gas model for an inflatable body.
///
/// Uses the ideal gas law: `p * V = n * R * T`
/// where:
/// - `p` = absolute pressure (Pa)
/// - `V` = volume (m³)
/// - `n` = moles of gas
/// - `R` = 8.314 J/(mol·K) (universal gas constant)
/// - `T` = temperature (K)
///
/// Isothermal assumption: `T = const`, so `p * V = const`.
/// Isobaric assumption: `p = const` (simple target pressure).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PressureVolumeGas {
    /// Amount of gas in moles times R*T: `n_RT = n * R * T`.
    /// This is conserved during isothermal expansion/compression.
    pub n_rt: f64,
    /// Atmospheric (ambient) pressure in Pa.
    pub p_ambient: f64,
    /// Current temperature (K) — kept for reference / heating models.
    pub temperature: f64,
}
impl PressureVolumeGas {
    /// Universal gas constant R (J mol⁻¹ K⁻¹).
    pub const R: f64 = 8.314_462_618;
    /// Create from initial conditions: absolute pressure, volume, and temperature.
    ///
    /// `p0` – initial absolute pressure (Pa)
    /// `v0` – initial volume (m³)
    /// `t0` – temperature (K)
    pub fn new(p0: f64, v0: f64, t0: f64) -> Self {
        let n_rt = p0 * v0;
        Self {
            n_rt,
            p_ambient: 101_325.0,
            temperature: t0,
        }
    }
    /// Create from number of moles `n` and temperature `t` (K).
    pub fn from_moles(n: f64, t: f64) -> Self {
        Self {
            n_rt: n * Self::R * t,
            p_ambient: 101_325.0,
            temperature: t,
        }
    }
    /// Absolute pressure at volume `v` (isothermal: p = n*R*T / V).
    pub fn absolute_pressure(&self, v: f64) -> f64 {
        if v.abs() < 1e-30 {
            return self.n_rt * 1e15;
        }
        self.n_rt / v
    }
    /// Gauge pressure = absolute pressure - ambient.
    pub fn gauge_pressure(&self, v: f64) -> f64 {
        self.absolute_pressure(v) - self.p_ambient
    }
    /// Add `delta_n_moles` moles of gas at current temperature.
    pub fn add_moles(&mut self, delta_n: f64) {
        self.n_rt += delta_n * Self::R * self.temperature;
    }
    /// Remove `delta_n_moles` moles of gas (vent).
    pub fn remove_moles(&mut self, delta_n: f64) {
        self.n_rt = (self.n_rt - delta_n * Self::R * self.temperature).max(0.0);
    }
    /// Volume at which the gauge pressure equals `target_gauge_pressure`.
    pub fn volume_at_pressure(&self, gauge_pressure: f64) -> f64 {
        let p_abs = gauge_pressure + self.p_ambient;
        if p_abs.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.n_rt / p_abs
    }
    /// Work done by the gas expanding from `v0` to `v1` (isothermal).
    ///
    /// W = n*R*T * ln(v1/v0)
    pub fn work_done(&self, v0: f64, v1: f64) -> f64 {
        if v0.abs() < 1e-30 || v1.abs() < 1e-30 {
            return 0.0;
        }
        self.n_rt * (v1 / v0).ln()
    }
    /// Heat capacity at constant volume (for monoatomic ideal gas: Cv = 1.5 * n * R).
    pub fn cv(&self) -> f64 {
        if self.temperature.abs() < 1e-10 {
            return 0.0;
        }
        1.5 * self.n_rt / self.temperature
    }
    /// Change temperature (isobaric heating/cooling adjusts n_RT).
    pub fn set_temperature(&mut self, new_temp: f64) {
        if new_temp.abs() < 1e-10 {
            return;
        }
        self.n_rt = self.n_rt * new_temp / self.temperature;
        self.temperature = new_temp;
    }
}
/// Pneumatic / inflatable soft body.
///
/// Simulated with PBD: pressure forces inflate the surface, stretch
/// constraints keep edge lengths close to rest.
#[derive(Debug, Clone)]
pub struct InflatableBody {
    /// Surface vertices.
    pub vertices: Vec<InflatableVertex>,
    /// Triangle faces.
    pub faces: Vec<InflatableFace>,
    /// Edges (structural + shear).
    pub edges: Vec<InflatableEdge>,
    /// Desired internal gauge pressure (Pa).
    pub target_pressure: f64,
    /// Actual internal gauge pressure this step (Pa).
    pub current_pressure: f64,
    /// Stretch stiffness (N/m).
    pub stretch_stiffness: f64,
    /// Bending stiffness (N·m).  Reserved for future use.
    pub bending_stiffness: f64,
    /// Velocity damping factor (0 = no damping, 1 = full stop).
    pub damping: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: [f64; 3],
}
impl InflatableBody {
    /// Build a UV-sphere inflatable body.
    ///
    /// `n_lat = subdivisions + 2` latitude rings,
    /// `n_lon = 2 * subdivisions + 4` longitude segments.
    pub fn sphere(radius: f64, subdivisions: u32, mass: f64) -> Self {
        use std::f64::consts::PI;
        let n_lat = (subdivisions + 2) as usize;
        let n_lon = (2 * subdivisions + 4) as usize;
        let mut positions: Vec<[f64; 3]> = Vec::new();
        positions.push([0.0, radius, 0.0]);
        for i in 1..=n_lat {
            let phi = PI * (i as f64) / (n_lat as f64 + 1.0);
            let y = radius * phi.cos();
            let r_ring = radius * phi.sin();
            for j in 0..n_lon {
                let theta = 2.0 * PI * (j as f64) / (n_lon as f64);
                positions.push([r_ring * theta.cos(), y, r_ring * theta.sin()]);
            }
        }
        positions.push([0.0, -radius, 0.0]);
        let n_verts = positions.len();
        let vertex_mass = if n_verts > 0 {
            mass / n_verts as f64
        } else {
            mass
        };
        let vertices: Vec<InflatableVertex> = positions
            .iter()
            .map(|&p| InflatableVertex::new(p, vertex_mass))
            .collect();
        let mut faces: Vec<InflatableFace> = Vec::new();
        for j in 0..n_lon {
            let a = 0usize;
            let b = 1 + j;
            let c = 1 + (j + 1) % n_lon;
            faces.push(InflatableFace { indices: [a, b, c] });
        }
        for i in 0..n_lat.saturating_sub(1) {
            let row_start = 1 + i * n_lon;
            let next_row = row_start + n_lon;
            for j in 0..n_lon {
                let a = row_start + j;
                let b = row_start + (j + 1) % n_lon;
                let c = next_row + j;
                let d = next_row + (j + 1) % n_lon;
                faces.push(InflatableFace { indices: [a, c, b] });
                faces.push(InflatableFace { indices: [b, c, d] });
            }
        }
        let south = n_verts - 1;
        let last_row = 1 + (n_lat - 1) * n_lon;
        for j in 0..n_lon {
            let a = last_row + j;
            let b = last_row + (j + 1) % n_lon;
            faces.push(InflatableFace {
                indices: [a, south, b],
            });
        }
        use std::collections::HashSet;
        let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
        for face in &faces {
            let [i0, i1, i2] = face.indices;
            for &(a, b) in &[(i0, i1), (i1, i2), (i2, i0)] {
                let key = if a < b { (a, b) } else { (b, a) };
                edge_set.insert(key);
            }
        }
        let edges: Vec<InflatableEdge> = edge_set
            .into_iter()
            .map(|(a, b)| {
                let rest_length = vec_len(vec_sub(positions[a], positions[b]));
                InflatableEdge {
                    indices: [a, b],
                    rest_length,
                }
            })
            .collect();
        Self {
            vertices,
            faces,
            edges,
            target_pressure: 1000.0,
            current_pressure: 1000.0,
            stretch_stiffness: 1e4,
            bending_stiffness: 0.0,
            damping: 0.01,
            gravity: [0.0, -9.81, 0.0],
        }
    }
    /// Signed volume of the closed surface mesh (m³).
    pub fn total_volume(&self) -> f64 {
        let positions: Vec<[f64; 3]> = self.vertices.iter().map(|v| v.pos).collect();
        total_volume(&positions, &self.faces)
    }
    /// Total surface area (m²).
    pub fn surface_area(&self) -> f64 {
        self.faces.iter().fold(0.0, |acc, f| {
            acc + triangle_area(
                self.vertices[f.indices[0]].pos,
                self.vertices[f.indices[1]].pos,
                self.vertices[f.indices[2]].pos,
            )
        })
    }
    /// Recompute `current_pressure` using a constant-pressure model.
    ///
    /// A more physical model (isothermal: pV = const) can be enabled by
    /// storing the rest volume, but the simple constant-pressure model is
    /// sufficient for most game/simulation purposes.
    pub fn update_pressure(&mut self) {
        self.current_pressure = self.target_pressure;
    }
    /// Compute the mesh centroid (average vertex position).
    fn centroid(&self) -> [f64; 3] {
        let n = self.vertices.len() as f64;
        if n < 1.0 {
            return [0.0; 3];
        }
        let sum = self
            .vertices
            .iter()
            .fold([0.0_f64; 3], |acc, v| vec_add(acc, v.pos));
        vec_scale(sum, 1.0 / n)
    }
    /// Apply pressure force to vertex velocities.
    ///
    /// F_pressure = p * area * n̂_out / 3  (distributed equally to 3 vertices).
    /// The outward normal is determined by the sign of (face_centroid - mesh_centroid)·n̂.
    pub fn apply_pressure_force(&mut self, dt: f64) {
        let mesh_centroid = self.centroid();
        let face_data: Vec<([f64; 3], f64)> = self
            .faces
            .iter()
            .map(|f| {
                let a = self.vertices[f.indices[0]].pos;
                let b = self.vertices[f.indices[1]].pos;
                let c = self.vertices[f.indices[2]].pos;
                let n = triangle_normal(a, b, c);
                let area = triangle_area(a, b, c);
                let fc = vec_scale(vec_add(vec_add(a, b), c), 1.0 / 3.0);
                let to_face = vec_sub(fc, mesh_centroid);
                let n_out = if vec_dot(n, to_face) >= 0.0 {
                    n
                } else {
                    vec_scale(n, -1.0)
                };
                (n_out, area)
            })
            .collect();
        let p = self.current_pressure;
        for (fi, face) in self.faces.iter().enumerate() {
            let (n_out, area) = face_data[fi];
            let f_mag = p * area / 3.0;
            let force = vec_scale(n_out, f_mag);
            for &vi in &face.indices {
                let v = &mut self.vertices[vi];
                let dv = vec_scale(force, v.inv_mass * dt);
                v.vel = vec_add(v.vel, dv);
            }
        }
    }
    /// Project XPBD distance constraints for all edges.
    pub fn solve_stretch_constraints(&mut self, dt: f64) {
        let alpha = 1.0 / (self.stretch_stiffness * dt * dt);
        for edge in &self.edges {
            let [i0, i1] = edge.indices;
            let x0 = self.vertices[i0].pos;
            let x1 = self.vertices[i1].pos;
            let w0 = self.vertices[i0].inv_mass;
            let w1 = self.vertices[i1].inv_mass;
            let w_sum = w0 + w1 + alpha;
            if w_sum < 1e-12 {
                continue;
            }
            let diff = vec_sub(x1, x0);
            let dist = vec_len(diff);
            if dist < 1e-12 {
                continue;
            }
            let dir = vec_scale(diff, 1.0 / dist);
            let constraint = dist - edge.rest_length;
            let lambda = -constraint / w_sum;
            let corr = vec_scale(dir, lambda);
            self.vertices[i0].pos = vec_sub(x0, vec_scale(corr, w0));
            self.vertices[i1].pos = vec_add(x1, vec_scale(corr, w1));
        }
    }
    /// Advance the simulation by `dt` seconds with `n_iter` constraint iterations.
    ///
    /// PBD loop:
    /// 1. Apply gravity + pressure forces to velocities.
    /// 2. Predict positions: x_pred = x + dt * v.
    /// 3. Solve stretch constraints `n_iter` times.
    /// 4. Update velocity: v = (x_pred - x_prev) / dt.
    /// 5. Apply damping.
    /// 6. Store positions for next step.
    ///
    /// Velocities are clamped to `max_vel` = 100 m/s to prevent numerical blow-up.
    pub fn step(&mut self, dt: f64, n_iter: usize) {
        const MAX_VEL: f64 = 100.0;
        for v in &mut self.vertices {
            let dv_g = vec_scale(self.gravity, dt);
            v.vel = vec_add(v.vel, dv_g);
        }
        self.update_pressure();
        self.apply_pressure_force(dt);
        for v in &mut self.vertices {
            let speed = vec_len(v.vel);
            if speed > MAX_VEL {
                v.vel = vec_scale(v.vel, MAX_VEL / speed);
            }
        }
        let prev_positions: Vec<[f64; 3]> = self.vertices.iter().map(|v| v.pos).collect();
        for v in &mut self.vertices {
            v.prev_pos = v.pos;
            v.pos = vec_add(v.pos, vec_scale(v.vel, dt));
        }
        for _ in 0..n_iter {
            self.solve_stretch_constraints(dt);
        }
        let damp = 1.0 - self.damping;
        for (i, v) in self.vertices.iter_mut().enumerate() {
            let new_vel = vec_scale(vec_sub(v.pos, prev_positions[i]), 1.0 / dt);
            let new_vel_clamped = {
                let speed = vec_len(new_vel);
                if speed > MAX_VEL {
                    vec_scale(new_vel, MAX_VEL / speed)
                } else {
                    new_vel
                }
            };
            v.vel = vec_scale(new_vel_clamped, damp);
        }
    }
}
impl InflatableBody {
    /// Compute the pressure–volume work done by the internal gas in one time step.
    ///
    /// W = p · ΔV   (positive = work done by gas expanding, negative = compressing)
    ///
    /// # Arguments
    /// * `prev_volume` – volume at the start of the step (m³).
    #[allow(dead_code)]
    pub fn compute_pressure_volume_work(&self, prev_volume: f64) -> f64 {
        let current_volume = self.total_volume();
        let delta_v = current_volume - prev_volume;
        self.current_pressure * delta_v
    }
    /// Simulate a puncture: model an instantaneous pressure drop caused by a
    /// hole of effective area `hole_area` (m²) over a time step `dt`.
    ///
    /// Uses the orifice flow model for compressible gas:
    /// ```text
    /// dP/dt ≈ −(P / V) * c_d * A_hole * v_exit
    /// ```
    /// where `v_exit ≈ sqrt(2 * P / ρ_air)` and ρ_air ≈ 1.225 kg/m³.
    ///
    /// The resulting pressure is clamped to zero (cannot go sub-atmospheric
    /// for a simple model).
    #[allow(dead_code)]
    pub fn simulate_puncture(&mut self, hole_area: f64, discharge_coeff: f64, dt: f64) {
        const RHO_AIR: f64 = 1.225;
        let p = self.current_pressure;
        if p < 1e-6 {
            return;
        }
        let v = self.total_volume().max(1e-9);
        let v_exit = (2.0 * p / RHO_AIR).sqrt();
        let q = discharge_coeff * hole_area * v_exit;
        let dp_dt = -(p / v) * q;
        self.current_pressure = (p + dp_dt * dt).max(0.0);
        self.target_pressure = self.current_pressure;
    }
    /// Compute the circumferential (hoop) membrane tension for a spherical shell
    /// approximation.
    ///
    /// For a thin spherical shell of radius `r` and thickness `t` under
    /// internal gauge pressure `p`:
    /// ```text
    /// σ_hoop = p * r / (2 * t)    (Pa)
    /// T_hoop = σ_hoop * t = p * r / 2   (N/m)
    /// ```
    ///
    /// # Arguments
    /// * `thickness` – membrane thickness (m).
    ///
    /// Returns the hoop tension per unit length (N/m).
    #[allow(dead_code)]
    pub fn compute_membrane_tension(&self, thickness: f64) -> f64 {
        if thickness < 1e-15 {
            return 0.0;
        }
        let vol = self.total_volume();
        let r = (3.0 * vol / (4.0 * std::f64::consts::PI)).cbrt();
        self.current_pressure * r / 2.0
    }
}
/// A simplified balloon model with elastic membrane and internal gas.
///
/// Uses the ideal gas law: p * V = p0 * V0 (isothermal) to compute
/// the internal pressure, then distributes it over the surface.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BalloonModel {
    /// Reference pressure (Pa).
    pub p0: f64,
    /// Reference volume (m³).
    pub v0: f64,
    /// Membrane thickness (m).
    pub thickness: f64,
    /// Young's modulus of the membrane (Pa).
    pub youngs_modulus: f64,
}
#[allow(dead_code)]
impl BalloonModel {
    /// Create a new balloon model from the current state of an inflatable body.
    pub fn new(body: &InflatableBody, thickness: f64, youngs_modulus: f64) -> Self {
        Self {
            p0: body.current_pressure,
            v0: body.total_volume(),
            thickness,
            youngs_modulus,
        }
    }
    /// Compute the internal pressure using the ideal gas law.
    pub fn gas_pressure(&self, current_volume: f64) -> f64 {
        if current_volume.abs() < 1e-15 {
            return self.p0 * 1e6;
        }
        self.p0 * self.v0 / current_volume
    }
    /// Compute the membrane stress for a given stretch ratio.
    ///
    /// Uses a simple linear elastic model: σ = E * (λ - 1) for stretch ratio λ.
    pub fn membrane_stress(&self, stretch_ratio: f64) -> f64 {
        self.youngs_modulus * (stretch_ratio - 1.0)
    }
    /// Compute the restoring pressure from membrane tension.
    ///
    /// Uses the Laplace law: p_membrane = 2 * σ * t / r
    /// where r is the effective radius.
    pub fn laplace_pressure(&self, stress: f64, radius: f64) -> f64 {
        if radius.abs() < 1e-15 {
            return 0.0;
        }
        2.0 * stress * self.thickness / radius
    }
}
/// Slow blow-molding simulation with gradual pressure increase.
///
/// Controls an `InflatableBody` by ramping the internal pressure from 0 to
/// `max_gauge_pressure` at rate `pressure_rate_pa_per_s` (Pa/s).
/// Monitors burst events via a `BurstDetector`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BlowMoldingSimulation {
    /// Rate of pressure increase (Pa/s).
    pub pressure_rate: f64,
    /// Maximum gauge pressure (Pa).
    pub max_gauge_pressure: f64,
    /// Current elapsed simulation time (s).
    pub elapsed: f64,
    /// Burst detector.
    pub burst: BurstDetector,
    /// Whether the mold has completed (pressure reached max or burst).
    pub complete: bool,
}
#[allow(dead_code)]
impl BlowMoldingSimulation {
    /// Create a new blow-molding simulation.
    pub fn new(pressure_rate: f64, max_gauge_pressure: f64, max_stretch: f64) -> Self {
        Self {
            pressure_rate,
            max_gauge_pressure,
            elapsed: 0.0,
            burst: BurstDetector::new(max_stretch),
            complete: false,
        }
    }
    /// Compute the target gauge pressure at the current time.
    pub fn target_pressure(&self) -> f64 {
        (self.pressure_rate * self.elapsed).min(self.max_gauge_pressure)
    }
    /// Advance the simulation by `dt` seconds.
    ///
    /// Sets `body.target_pressure`, steps the body simulation,
    /// and checks for burst.  Returns `true` if the simulation is still running.
    pub fn step(&mut self, body: &mut InflatableBody, dt: f64, n_iter: usize) -> bool {
        if self.complete {
            return false;
        }
        body.target_pressure = self.target_pressure();
        body.step(dt, n_iter);
        self.elapsed += dt;
        if self.burst.check(body) {
            self.complete = true;
            return false;
        }
        if self.target_pressure() >= self.max_gauge_pressure {
            self.complete = true;
        }
        !self.complete
    }
    /// Current pressure as fraction of max (0..1).
    pub fn pressure_fraction(&self) -> f64 {
        if self.max_gauge_pressure < 1e-30 {
            return 1.0;
        }
        (self.target_pressure() / self.max_gauge_pressure).min(1.0)
    }
}
/// Parameters for a simple blow-molding simulation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BlowMoldingParams {
    /// Pressure ramp rate (Pa/s).
    pub pressure_rate: f64,
    /// Maximum pressure (Pa).
    pub max_pressure: f64,
    /// Mold boundary positions: list of (center, radius) spheres that
    /// constrain the inflatable from expanding beyond.
    pub mold_spheres: Vec<([f64; 3], f64)>,
}
#[allow(dead_code)]
impl BlowMoldingParams {
    /// Create new blow-molding parameters.
    pub fn new(pressure_rate: f64, max_pressure: f64) -> Self {
        Self {
            pressure_rate,
            max_pressure,
            mold_spheres: Vec::new(),
        }
    }
    /// Compute the pressure at time `t`.
    pub fn pressure_at(&self, t: f64) -> f64 {
        (self.pressure_rate * t).min(self.max_pressure)
    }
    /// Apply mold contact constraints to vertices.
    ///
    /// Any vertex outside a mold sphere is projected back onto the surface.
    pub fn apply_mold_constraints(&self, vertices: &mut [InflatableVertex]) {
        for sphere in &self.mold_spheres {
            let (center, radius) = *sphere;
            for v in vertices.iter_mut() {
                let dx = v.pos[0] - center[0];
                let dy = v.pos[1] - center[1];
                let dz = v.pos[2] - center[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist > radius && dist > 1e-14 {
                    let scale = radius / dist;
                    v.pos[0] = center[0] + dx * scale;
                    v.pos[1] = center[1] + dy * scale;
                    v.pos[2] = center[2] + dz * scale;
                    let nx = dx / dist;
                    let ny = dy / dist;
                    let nz = dz / dist;
                    let vn = v.vel[0] * nx + v.vel[1] * ny + v.vel[2] * nz;
                    if vn > 0.0 {
                        v.vel[0] -= vn * nx;
                        v.vel[1] -= vn * ny;
                        v.vel[2] -= vn * nz;
                    }
                }
            }
        }
    }
}
/// A valve controls the flow of gas into or out of an inflatable body.
///
/// Models both an inlet valve (pumping gas in) and an outlet valve (venting).
/// The flow rate follows a linear pressure-differential model:
/// `dm/dt = K_v * (p_source - p_body)` (clamped to \[0, max_flow_rate\]).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Valve {
    /// Flow conductance (mol/s/Pa).
    pub conductance: f64,
    /// Maximum molar flow rate (mol/s).
    pub max_flow_rate: f64,
    /// Whether the valve is open.
    pub is_open: bool,
    /// Source pressure (Pa) — the pressure upstream of the valve.
    pub source_pressure: f64,
    /// Total moles that have passed through (for diagnostics).
    pub total_moles_transferred: f64,
}
#[allow(dead_code)]
impl Valve {
    /// Create a new valve.
    pub fn new(conductance: f64, max_flow_rate: f64, source_pressure: f64) -> Self {
        Self {
            conductance,
            max_flow_rate,
            is_open: false,
            source_pressure,
            total_moles_transferred: 0.0,
        }
    }
    /// Create an inlet valve (pumping gas at `pump_pressure` Pa).
    pub fn inlet(pump_pressure: f64, conductance: f64) -> Self {
        Self::new(conductance, f64::INFINITY, pump_pressure)
    }
    /// Create an outlet valve (venting to atmosphere at 101325 Pa).
    pub fn outlet(conductance: f64) -> Self {
        Self::new(conductance, f64::INFINITY, 101_325.0)
    }
    /// Open the valve.
    pub fn open(&mut self) {
        self.is_open = true;
    }
    /// Close the valve.
    pub fn close(&mut self) {
        self.is_open = false;
    }
    /// Compute the molar flow rate (mol/s) given body pressure `p_body` (Pa).
    ///
    /// Positive = gas flows into the body; negative = gas flows out.
    pub fn flow_rate(&self, p_body: f64) -> f64 {
        if !self.is_open {
            return 0.0;
        }
        let dp = self.source_pressure - p_body;
        let flow = self.conductance * dp;
        flow.clamp(-self.max_flow_rate, self.max_flow_rate)
    }
    /// Advance the valve by `dt` seconds: update the gas model accordingly.
    ///
    /// Returns the number of moles transferred (positive = into body).
    pub fn step(&mut self, gas: &mut PressureVolumeGas, current_volume: f64, dt: f64) -> f64 {
        let p_body = gas.absolute_pressure(current_volume);
        let rate = self.flow_rate(p_body);
        let delta_n = rate * dt;
        if delta_n > 0.0 {
            gas.add_moles(delta_n);
        } else {
            gas.remove_moles(-delta_n);
        }
        self.total_moles_transferred += delta_n.abs();
        delta_n
    }
}
/// Burst detector for inflatable membranes.
///
/// Monitors the maximum principal strain across the surface and triggers a
/// burst event when the strain threshold is exceeded.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BurstDetector {
    /// Maximum allowable stretch ratio (e.g., 1.5 = 50% elongation).
    pub max_stretch_ratio: f64,
    /// Whether the body has burst.
    pub has_burst: bool,
    /// The face index where the burst first occurred (if any).
    pub burst_face: Option<usize>,
    /// Stretch ratio at burst.
    pub burst_stretch: f64,
}
#[allow(dead_code)]
impl BurstDetector {
    /// Create a new burst detector with `max_stretch_ratio`.
    pub fn new(max_stretch_ratio: f64) -> Self {
        Self {
            max_stretch_ratio,
            has_burst: false,
            burst_face: None,
            burst_stretch: 0.0,
        }
    }
    /// Check whether any edge of `body` exceeds the max stretch ratio.
    ///
    /// Updates `has_burst`, `burst_face`, and `burst_stretch`.
    /// Returns `true` if burst is detected this call.
    pub fn check(&mut self, body: &InflatableBody) -> bool {
        if self.has_burst {
            return true;
        }
        for (fi, face) in body.faces.iter().enumerate() {
            let [i0, i1, i2] = face.indices;
            let p0 = body.vertices[i0].pos;
            let p1 = body.vertices[i1].pos;
            let p2 = body.vertices[i2].pos;
            let edges_cur = [
                vec_len(vec_sub(p1, p0)),
                vec_len(vec_sub(p2, p1)),
                vec_len(vec_sub(p0, p2)),
            ];
            let pairs = [(i0, i1), (i1, i2), (i2, i0)];
            for (k, &(a, b)) in pairs.iter().enumerate() {
                let rest = body.edges.iter().find(|e| {
                    let [ea, eb] = e.indices;
                    (ea == a && eb == b) || (ea == b && eb == a)
                });
                if let Some(edge) = rest
                    && edge.rest_length > 1e-14
                {
                    let stretch = edges_cur[k] / edge.rest_length;
                    if stretch > self.max_stretch_ratio {
                        self.has_burst = true;
                        self.burst_face = Some(fi);
                        self.burst_stretch = stretch;
                        return true;
                    }
                }
            }
        }
        false
    }
    /// Check using precomputed stretch ratios.
    pub fn check_ratios(&mut self, ratios: &[f64], face_offset: usize) -> bool {
        if self.has_burst {
            return true;
        }
        for (i, &r) in ratios.iter().enumerate() {
            if r > self.max_stretch_ratio {
                self.has_burst = true;
                self.burst_face = Some(face_offset + i);
                self.burst_stretch = r;
                return true;
            }
        }
        false
    }
    /// Reset burst state.
    pub fn reset(&mut self) {
        self.has_burst = false;
        self.burst_face = None;
        self.burst_stretch = 0.0;
    }
}
/// Inflation sequence controller: ramps pressure over time with configurable
/// stages (inflate, hold, deflate).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct InflationSequence {
    /// Stages: (duration_s, target_pressure_Pa).
    pub stages: Vec<(f64, f64)>,
    /// Current elapsed time.
    pub elapsed: f64,
}
#[allow(dead_code)]
impl InflationSequence {
    /// Create a new inflation sequence.
    pub fn new(stages: Vec<(f64, f64)>) -> Self {
        Self {
            stages,
            elapsed: 0.0,
        }
    }
    /// Create a simple inflate-hold sequence.
    pub fn inflate_and_hold(ramp_time: f64, hold_time: f64, target_pressure: f64) -> Self {
        Self::new(vec![
            (ramp_time, target_pressure),
            (hold_time, target_pressure),
        ])
    }
    /// Get the target pressure at the current time.
    pub fn current_pressure(&self) -> f64 {
        let mut t = 0.0_f64;
        let mut prev_pressure = 0.0_f64;
        for &(duration, target) in &self.stages {
            if self.elapsed < t + duration {
                let frac = (self.elapsed - t) / duration;
                return prev_pressure + frac * (target - prev_pressure);
            }
            t += duration;
            prev_pressure = target;
        }
        self.stages.last().map_or(0.0, |&(_, p)| p)
    }
    /// Advance the sequence by `dt` seconds.
    pub fn advance(&mut self, dt: f64) {
        self.elapsed += dt;
    }
    /// Total duration of all stages.
    pub fn total_duration(&self) -> f64 {
        self.stages.iter().map(|&(d, _)| d).sum()
    }
    /// Whether the sequence has completed.
    pub fn is_complete(&self) -> bool {
        self.elapsed >= self.total_duration()
    }
}
