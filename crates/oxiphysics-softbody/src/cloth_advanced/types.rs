//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Hashin-type failure criterion for woven composites applied to fabric.
#[derive(Debug, Clone)]
pub struct FabricFailureCriterion {
    /// Tensile warp strength (Pa).
    pub sigma_warp_t: f64,
    /// Compressive warp strength (Pa).
    pub sigma_warp_c: f64,
    /// Tensile weft strength (Pa).
    pub sigma_weft_t: f64,
    /// Compressive weft strength (Pa).
    pub sigma_weft_c: f64,
    /// Shear strength (Pa).
    pub tau_s: f64,
}
impl FabricFailureCriterion {
    /// Create for cotton plain weave.
    pub fn cotton_plain() -> Self {
        Self {
            sigma_warp_t: 400e6,
            sigma_warp_c: 200e6,
            sigma_weft_t: 350e6,
            sigma_weft_c: 180e6,
            tau_s: 80e6,
        }
    }
    /// Hashin warp-tension failure index (< 1 = safe).
    pub fn warp_tension_index(&self, sigma_1: f64, tau_12: f64) -> f64 {
        (sigma_1 / self.sigma_warp_t).powi(2) + (tau_12 / self.tau_s).powi(2)
    }
    /// Hashin weft-tension failure index.
    pub fn weft_tension_index(&self, sigma_2: f64, tau_12: f64) -> f64 {
        (sigma_2 / self.sigma_weft_t).powi(2) + (tau_12 / self.tau_s).powi(2)
    }
    /// Maximum of all failure indices (> 1 = failure).
    pub fn max_failure_index(&self, sigma_1: f64, sigma_2: f64, tau_12: f64) -> f64 {
        let f1 = if sigma_1 >= 0.0 {
            self.warp_tension_index(sigma_1, tau_12)
        } else {
            (sigma_1 / self.sigma_warp_c).powi(2)
        };
        let f2 = if sigma_2 >= 0.0 {
            self.weft_tension_index(sigma_2, tau_12)
        } else {
            (sigma_2 / self.sigma_weft_c).powi(2)
        };
        let fs = (tau_12 / self.tau_s).powi(2);
        f1.max(f2).max(fs)
    }
    /// Check if failed.
    pub fn is_failed(&self, sigma_1: f64, sigma_2: f64, tau_12: f64) -> bool {
        self.max_failure_index(sigma_1, sigma_2, tau_12) >= 1.0
    }
}
/// Sewn cloth panels with seam constraints.
#[derive(Debug, Clone)]
pub struct ClothSewing {
    /// Particles in all panels (flattened).
    pub particles: Vec<ClothParticle>,
    /// Panel boundary: `panel_starts[i]` is the start index of panel `i`.
    pub panel_starts: Vec<usize>,
    /// Seam constraints.
    pub seams: Vec<SeamLine>,
}
impl ClothSewing {
    /// Create a new sewing system with `num_panels` of `particles_per_panel` particles.
    pub fn new(num_panels: usize, particles_per_panel: usize, mass: f64) -> Self {
        let total = num_panels * particles_per_panel;
        let particles = (0..total)
            .map(|i| ClothParticle::new([i as f64 * 0.01, 0.0, 0.0], mass))
            .collect();
        let panel_starts = (0..num_panels).map(|i| i * particles_per_panel).collect();
        Self {
            particles,
            panel_starts,
            seams: Vec::new(),
        }
    }
    /// Add a seam between panels.
    pub fn add_seam(&mut self, seam: SeamLine) {
        self.seams.push(seam);
    }
    /// Apply seam closure constraints (pulls seam particles together).
    pub fn apply_seam_constraints(&mut self, dt: f64) {
        for seam in &self.seams {
            let np = seam.num_pairs();
            for i in 0..np {
                let ia = seam.panel_a_indices[i];
                let ib = seam.panel_b_indices[i];
                if ia >= self.particles.len() || ib >= self.particles.len() {
                    continue;
                }
                let pa = self.particles[ia].position;
                let pb = self.particles[ib].position;
                let diff = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                let dist = norm3(diff);
                if dist < 1e-15 {
                    continue;
                }
                let dir = normalize3(diff);
                let correction = seam.stiffness * dist * dt;
                if !self.particles[ia].is_pinned() {
                    for (v, d) in self.particles[ia].velocity.iter_mut().zip(dir.iter()) {
                        *v += correction * d;
                    }
                }
                if !self.particles[ib].is_pinned() {
                    for (v, d) in self.particles[ib].velocity.iter_mut().zip(dir.iter()) {
                        *v -= correction * d;
                    }
                }
            }
        }
    }
}
/// Broad-phase cloth-cloth collision using sorted AABB sweep.
#[derive(Debug, Clone)]
pub struct ClothCollisionBroadphase {
    /// AABB leaves for each triangle.
    pub leaves: Vec<ClothBvhLeaf>,
    /// Collision pairs found.
    pub pairs: Vec<(usize, usize)>,
}
impl ClothCollisionBroadphase {
    /// Create broad-phase from triangle list.
    pub fn new() -> Self {
        Self {
            leaves: Vec::new(),
            pairs: Vec::new(),
        }
    }
    /// Add a triangle to the broad phase.
    pub fn add_triangle(&mut self, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], margin: f64) {
        let idx = self.leaves.len();
        self.leaves
            .push(ClothBvhLeaf::from_triangle(idx, p0, p1, p2, margin));
    }
    /// Find all overlapping AABB pairs (brute force for small sets).
    pub fn find_pairs(&mut self) {
        self.pairs.clear();
        let n = self.leaves.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if self.leaves[i].overlaps(&self.leaves[j]) {
                    self.pairs.push((i, j));
                }
            }
        }
    }
    /// Number of candidate collision pairs.
    pub fn pair_count(&self) -> usize {
        self.pairs.len()
    }
}
/// Yarn cross-section geometry.
#[derive(Debug, Clone)]
pub struct YarnSection {
    /// Radius of circular cross-section (m).
    pub radius: f64,
    /// Linear density (tex = g/km).
    pub linear_density: f64,
    /// Young's modulus of yarn (Pa).
    pub elastic_modulus: f64,
    /// Bending stiffness EI (N·m²).
    pub bending_stiffness: f64,
    /// Torsional stiffness GJ (N·m²).
    pub torsional_stiffness: f64,
}
impl YarnSection {
    /// Create from radius and elastic modulus (assumes circular cross-section).
    pub fn circular(radius: f64, e: f64, g: f64) -> Self {
        let i = std::f64::consts::PI * radius.powi(4) / 4.0;
        let j = 2.0 * i;
        Self {
            radius,
            linear_density: 1.0,
            elastic_modulus: e,
            bending_stiffness: e * i,
            torsional_stiffness: g * j,
        }
    }
    /// Axial stiffness EA.
    pub fn axial_stiffness(&self) -> f64 {
        let area = std::f64::consts::PI * self.radius * self.radius;
        self.elastic_modulus * area
    }
}
/// A 2D garment panel (pattern piece) before construction.
#[derive(Debug, Clone)]
pub struct GarmentPanel {
    /// Panel name.
    pub name: String,
    /// 2D vertices of the pattern outline (m).
    pub vertices: Vec<[f64; 2]>,
    /// Grainline direction (2D unit vector, typically along warp).
    pub grainline: [f64; 2],
    /// Seam allowance width (m).
    pub seam_allowance: f64,
}
impl GarmentPanel {
    /// Create a rectangular panel.
    pub fn rectangle(name: &str, width: f64, height: f64) -> Self {
        Self {
            name: name.to_string(),
            vertices: vec![[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]],
            grainline: [0.0, 1.0],
            seam_allowance: 0.015,
        }
    }
    /// Panel area (shoelace formula).
    pub fn area(&self) -> f64 {
        let n = self.vertices.len();
        if n < 3 {
            return 0.0;
        }
        let mut sum = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            let [xi, yi] = self.vertices[i];
            let [xj, yj] = self.vertices[j];
            sum += xi * yj - xj * yi;
        }
        0.5 * sum.abs()
    }
    /// Perimeter (sum of edge lengths).
    pub fn perimeter(&self) -> f64 {
        let n = self.vertices.len();
        let mut total = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            let [xi, yi] = self.vertices[i];
            let [xj, yj] = self.vertices[j];
            total += ((xj - xi).powi(2) + (yj - yi).powi(2)).sqrt();
        }
        total
    }
    /// Centroid of the panel.
    pub fn centroid(&self) -> [f64; 2] {
        if self.vertices.is_empty() {
            return [0.0, 0.0];
        }
        let n = self.vertices.len() as f64;
        let cx = self.vertices.iter().map(|v| v[0]).sum::<f64>() / n;
        let cy = self.vertices.iter().map(|v| v[1]).sum::<f64>() / n;
        [cx, cy]
    }
    /// Sewing length excluding seam allowances (usable seam).
    pub fn net_perimeter(&self) -> f64 {
        self.perimeter() - self.seam_allowance * self.vertices.len() as f64
    }
}
/// Thermal properties of fabric for heat transfer simulation.
#[derive(Debug, Clone)]
pub struct FabricThermal {
    /// Thermal conductivity (W/m·K).
    pub conductivity: f64,
    /// Specific heat capacity (J/kg·K).
    pub specific_heat: f64,
    /// Emissivity (0–1).
    pub emissivity: f64,
    /// Moisture vapor permeability (kg/m²·s·Pa).
    pub vapor_permeability: f64,
    /// Thickness (m).
    pub thickness: f64,
    /// Areal density (kg/m²).
    pub areal_density: f64,
}
impl FabricThermal {
    /// Cotton thermal properties.
    pub fn cotton() -> Self {
        Self {
            conductivity: 0.071,
            specific_heat: 1340.0,
            emissivity: 0.95,
            vapor_permeability: 1e-11,
            thickness: 0.3e-3,
            areal_density: 0.12,
        }
    }
    /// Polyester thermal properties.
    pub fn polyester() -> Self {
        Self {
            conductivity: 0.045,
            specific_heat: 1380.0,
            emissivity: 0.90,
            vapor_permeability: 5e-13,
            thickness: 0.25e-3,
            areal_density: 0.10,
        }
    }
    /// Thermal resistance (CLO unit = 0.155 m²K/W).
    pub fn thermal_resistance(&self) -> f64 {
        self.thickness / self.conductivity
    }
    /// CLO value (clothing insulation unit).
    pub fn clo_value(&self) -> f64 {
        self.thermal_resistance() / 0.155
    }
    /// Transient temperature distribution (1D slab, semi-infinite).
    ///
    /// Temperature at surface after heating to `t_surface` for time `t` (s).
    /// `x` — depth from surface (m).
    pub fn temperature_at_depth(&self, t_surface: f64, t_initial: f64, x: f64, t: f64) -> f64 {
        let density = self.areal_density / self.thickness;
        let alpha = self.conductivity / (density * self.specific_heat);
        if alpha <= 0.0 || t <= 0.0 {
            return t_initial;
        }
        let xi = x / (2.0 * (alpha * t).sqrt());
        let erfc_val = erfc_fabric(xi);
        t_initial + (t_surface - t_initial) * erfc_val
    }
}
/// Anisotropic elastic fabric model fitted from biaxial extension test data.
#[derive(Debug, Clone)]
pub struct ElasticFabric {
    /// Fitted stiffness in x direction (warp).
    pub kx: f64,
    /// Fitted stiffness in y direction (weft).
    pub ky: f64,
    /// Poisson ratio xy.
    pub nu_xy: f64,
    /// Shear stiffness.
    pub ks: f64,
    /// Rest area.
    pub rest_area: f64,
}
impl ElasticFabric {
    /// Create an elastic fabric model.
    pub fn new(kx: f64, ky: f64, nu_xy: f64, ks: f64, rest_area: f64) -> Self {
        Self {
            kx,
            ky,
            nu_xy,
            ks,
            rest_area,
        }
    }
    /// Compute the 2×2 stiffness tensor `C` for in-plane deformation.
    ///
    /// Returns `[[Cxx, Cxy\], [Cxy, Cyy]]`.
    pub fn stiffness_tensor(&self) -> [[f64; 2]; 2] {
        let denom = 1.0 - self.nu_xy * self.nu_xy * self.ky / self.kx;
        if denom.abs() < 1e-15 {
            return [[self.kx, 0.0], [0.0, self.ky]];
        }
        [
            [self.kx / denom, self.nu_xy * self.ky / denom],
            [self.nu_xy * self.ky / denom, self.ky / denom],
        ]
    }
    /// In-plane force given strain `[ex, ey]`.
    ///
    /// Returns `[Nx, Ny]` (N/m).
    pub fn in_plane_force(&self, strain: [f64; 2]) -> [f64; 2] {
        let c = self.stiffness_tensor();
        [
            c[0][0] * strain[0] + c[0][1] * strain[1],
            c[1][0] * strain[0] + c[1][1] * strain[1],
        ]
    }
    /// Compute the energy density for a given strain state.
    pub fn strain_energy_density(&self, strain: [f64; 2]) -> f64 {
        let f = self.in_plane_force(strain);
        0.5 * (f[0] * strain[0] + f[1] * strain[1])
    }
}
/// Woven fabric mechanical model.
#[derive(Debug, Clone)]
pub struct WovenFabric {
    /// Warp yarn section.
    pub warp: YarnSection,
    /// Weft yarn section.
    pub weft: YarnSection,
    /// Weave pattern.
    pub pattern: WeavePattern,
    /// Warp thread density (threads/m).
    pub warp_density: f64,
    /// Weft thread density (threads/m).
    pub weft_density: f64,
    /// Friction coefficient between yarns.
    pub yarn_friction: f64,
}
impl WovenFabric {
    /// Create a basic plain-weave cotton fabric.
    pub fn cotton_plain() -> Self {
        let yarn = YarnSection::circular(0.15e-3, 8e9, 3e9);
        Self {
            warp: yarn.clone(),
            weft: yarn,
            pattern: WeavePattern::Plain,
            warp_density: 3000.0,
            weft_density: 2800.0,
            yarn_friction: 0.3,
        }
    }
    /// Fabric areal density (kg/m²).
    pub fn areal_density(&self) -> f64 {
        let warp_mass = self.warp_density * self.warp.linear_density * 1e-6;
        let weft_mass = self.weft_density * self.weft.linear_density * 1e-6;
        warp_mass + weft_mass
    }
    /// Warp-direction tensile stiffness per unit width (N/m).
    pub fn warp_stiffness(&self) -> f64 {
        let ea = self.warp.axial_stiffness();
        let crimp = self.pattern.crimp_factor();
        ea * self.warp_density * (1.0 - crimp)
    }
    /// Weft-direction tensile stiffness per unit width.
    pub fn weft_stiffness(&self) -> f64 {
        let ea = self.weft.axial_stiffness();
        let crimp = self.pattern.crimp_factor();
        ea * self.weft_density * (1.0 - crimp)
    }
    /// Shear stiffness (Kawabata model approximation).
    pub fn shear_stiffness(&self) -> f64 {
        let n_contact =
            self.warp_density * self.weft_density * 4.0 * self.warp.radius * self.weft.radius;
        let pitch = 1.0 / self.warp_density.max(1.0);
        self.yarn_friction * n_contact * pitch
    }
    /// Thickness of fabric (approximate: 2 yarn diameters).
    pub fn thickness(&self) -> f64 {
        2.0 * (self.warp.radius + self.weft.radius)
    }
}
/// Water absorption simulation for cloth.
#[derive(Debug, Clone)]
pub struct ClothWetting {
    /// Water content per particle (kg/m²).
    pub water_content: Vec<f64>,
    /// Dry surface mass density (kg/m²).
    pub dry_density: f64,
    /// Maximum water capacity (kg/m²).
    pub max_capacity: f64,
    /// Absorption rate (m²/s).
    pub absorption_rate: f64,
    /// Evaporation rate (kg/m²/s).
    pub evaporation_rate: f64,
}
impl ClothWetting {
    /// Create a new wetting model for `n` particles.
    pub fn new(n: usize, dry_density: f64, max_capacity: f64) -> Self {
        Self {
            water_content: vec![0.0; n],
            dry_density,
            max_capacity,
            absorption_rate: 1e-3,
            evaporation_rate: 1e-4,
        }
    }
    /// Apply wetting from a water source to particle `i`.
    pub fn wet_particle(&mut self, i: usize, amount: f64) {
        if i < self.water_content.len() {
            self.water_content[i] = (self.water_content[i] + amount).min(self.max_capacity);
        }
    }
    /// Effective mass density at particle `i` (dry + water).
    pub fn effective_density(&self, i: usize) -> f64 {
        let w = if i < self.water_content.len() {
            self.water_content[i]
        } else {
            0.0
        };
        self.dry_density + w
    }
    /// Update evaporation for one time step.
    pub fn update_evaporation(&mut self, dt: f64) {
        for w in &mut self.water_content {
            *w = (*w - self.evaporation_rate * dt).max(0.0);
        }
    }
    /// Total water mass in the cloth.
    pub fn total_water_mass(&self) -> f64 {
        self.water_content.iter().sum()
    }
}
/// A cloth simulation particle.
#[derive(Debug, Clone)]
pub struct ClothParticle {
    /// Position in world space.
    pub position: [f64; 3],
    /// Velocity.
    pub velocity: [f64; 3],
    /// Accumulated external force.
    pub force: [f64; 3],
    /// Inverse mass (0 = static/pinned).
    pub inv_mass: f64,
}
impl ClothParticle {
    /// Create a new cloth particle.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 1e-15 { 1.0 / mass } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass,
        }
    }
    /// Create a pinned (static) particle.
    pub fn pinned(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass: 0.0,
        }
    }
    /// Returns `true` if this particle is pinned (static).
    pub fn is_pinned(&self) -> bool {
        self.inv_mass < 1e-15
    }
}
/// Cloth render mesh created by subdividing the simulation mesh.
#[derive(Debug, Clone)]
pub struct ClothRenderMesh {
    /// Render vertices.
    pub vertices: Vec<RenderVertex>,
    /// Triangle indices.
    pub indices: Vec<[u32; 3]>,
}
impl ClothRenderMesh {
    /// Build a render mesh from cloth particles and triangle connectivity.
    ///
    /// Computes smooth normals by averaging face normals at each vertex.
    pub fn from_particles(particles: &[ClothParticle], triangles: &[[usize; 3]]) -> Self {
        let n = particles.len();
        let mut normals = vec![[0.0f32; 3]; n];
        let mut norm_count = vec![0u32; n];
        for tri in triangles {
            let [i0, i1, i2] = *tri;
            if i0 >= n || i1 >= n || i2 >= n {
                continue;
            }
            let p0 = particles[i0].position;
            let p1 = particles[i1].position;
            let p2 = particles[i2].position;
            let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            let fn_ = normalize3(cross3(e1, e2));
            for &vi in &[i0, i1, i2] {
                normals[vi][0] += fn_[0] as f32;
                normals[vi][1] += fn_[1] as f32;
                normals[vi][2] += fn_[2] as f32;
                norm_count[vi] += 1;
            }
        }
        let vertices: Vec<RenderVertex> = particles
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let cnt = norm_count[i].max(1) as f32;
                let n_raw = [
                    normals[i][0] / cnt,
                    normals[i][1] / cnt,
                    normals[i][2] / cnt,
                ];
                let n_len =
                    (n_raw[0] * n_raw[0] + n_raw[1] * n_raw[1] + n_raw[2] * n_raw[2]).sqrt();
                let nm = if n_len > 1e-7 {
                    [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len]
                } else {
                    [0.0, 1.0, 0.0]
                };
                RenderVertex {
                    position: [
                        p.position[0] as f32,
                        p.position[1] as f32,
                        p.position[2] as f32,
                    ],
                    normal: nm,
                    uv: [i as f32 / n as f32, 0.0],
                }
            })
            .collect();
        let indices: Vec<[u32; 3]> = triangles
            .iter()
            .filter_map(|&[i0, i1, i2]| {
                if i0 < n && i1 < n && i2 < n {
                    Some([i0 as u32, i1 as u32, i2 as u32])
                } else {
                    None
                }
            })
            .collect();
        Self { vertices, indices }
    }
    /// Number of triangles.
    pub fn num_triangles(&self) -> usize {
        self.indices.len()
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
}
/// Edge element for tearing simulation.
#[derive(Debug, Clone)]
pub struct ClothEdge {
    /// Index of particle A.
    pub a: usize,
    /// Index of particle B.
    pub b: usize,
    /// Rest length.
    pub rest_length: f64,
    /// Maximum strain before rupture.
    pub rupture_strain: f64,
    /// Whether this edge has torn.
    pub torn: bool,
}
impl ClothEdge {
    /// Create a new cloth edge.
    pub fn new(a: usize, b: usize, rest_length: f64, rupture_strain: f64) -> Self {
        Self {
            a,
            b,
            rest_length,
            rupture_strain,
            torn: false,
        }
    }
    /// Compute current strain given particle positions.
    pub fn strain(&self, particles: &[ClothParticle]) -> f64 {
        if self.a >= particles.len() || self.b >= particles.len() {
            return 0.0;
        }
        let pa = particles[self.a].position;
        let pb = particles[self.b].position;
        let dist = norm3([pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]]);
        (dist - self.rest_length) / self.rest_length
    }
    /// Check and mark the edge as torn if strain exceeds rupture threshold.
    pub fn check_rupture(&mut self, particles: &[ClothParticle]) -> bool {
        if !self.torn && self.strain(particles) > self.rupture_strain {
            self.torn = true;
        }
        self.torn
    }
}
/// Mode of cloth interaction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InteractionMode {
    /// Pinch-grabbing a cloth vertex.
    Grab,
    /// Pulling a cloth vertex toward a target.
    Pull,
    /// Cutting along a line.
    Cut,
}
/// Bounding volume hierarchy leaf node for cloth collision.
#[derive(Debug, Clone)]
pub struct ClothBvhLeaf {
    /// Triangle index.
    pub tri_idx: usize,
    /// Axis-aligned bounding box: \[min_x, min_y, min_z, max_x, max_y, max_z\].
    pub aabb: [f64; 6],
}
impl ClothBvhLeaf {
    /// Create from triangle vertex positions.
    pub fn from_triangle(
        idx: usize,
        p0: [f64; 3],
        p1: [f64; 3],
        p2: [f64; 3],
        margin: f64,
    ) -> Self {
        let min_x = p0[0].min(p1[0]).min(p2[0]) - margin;
        let min_y = p0[1].min(p1[1]).min(p2[1]) - margin;
        let min_z = p0[2].min(p1[2]).min(p2[2]) - margin;
        let max_x = p0[0].max(p1[0]).max(p2[0]) + margin;
        let max_y = p0[1].max(p1[1]).max(p2[1]) + margin;
        let max_z = p0[2].max(p1[2]).max(p2[2]) + margin;
        Self {
            tri_idx: idx,
            aabb: [min_x, min_y, min_z, max_x, max_y, max_z],
        }
    }
    /// Check AABB overlap with another leaf.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.aabb[0] <= other.aabb[3]
            && self.aabb[3] >= other.aabb[0]
            && self.aabb[1] <= other.aabb[4]
            && self.aabb[4] >= other.aabb[1]
            && self.aabb[2] <= other.aabb[5]
            && self.aabb[5] >= other.aabb[2]
    }
}
/// Knitted fabric loop model.
#[derive(Debug, Clone)]
pub struct KnittedFabric {
    /// Yarn section.
    pub yarn: YarnSection,
    /// Loop length (m) — controls stitch size.
    pub loop_length: f64,
    /// Courses per unit length (1/m).
    pub courses: f64,
    /// Wales per unit width (1/m).
    pub wales: f64,
    /// Knit structure: true = jersey, false = rib.
    pub is_jersey: bool,
}
impl KnittedFabric {
    /// Create a typical jersey cotton knit.
    pub fn jersey_cotton() -> Self {
        let yarn = YarnSection::circular(0.1e-3, 5e9, 2e9);
        Self {
            yarn,
            loop_length: 4e-3,
            courses: 150.0,
            wales: 120.0,
            is_jersey: true,
        }
    }
    /// Stitch density (loops/m²).
    pub fn stitch_density(&self) -> f64 {
        self.courses * self.wales
    }
    /// Areal density (kg/m²).
    pub fn areal_density(&self) -> f64 {
        self.stitch_density() * self.loop_length * self.yarn.linear_density * 1e-6
    }
    /// Extensibility: normalized stretch at 50% of break force.
    pub fn extensibility(&self) -> f64 {
        if self.is_jersey { 1.5 } else { 0.8 }
    }
    /// Loop geometry: loop height from Peirce model.
    pub fn loop_height(&self) -> f64 {
        self.loop_length / std::f64::consts::PI
    }
    /// Bending rigidity of fabric (N·m²/m width).
    pub fn bending_rigidity(&self) -> f64 {
        self.yarn.bending_stiffness * self.wales
    }
}
/// Woven fabric anisotropy model (plain weave approximation).
#[derive(Debug, Clone)]
pub struct FabricWeave {
    /// Warp direction (unit vector).
    pub warp: [f64; 3],
    /// Weft direction (unit vector).
    pub weft: [f64; 3],
    /// Stiffness along warp direction (N/m).
    pub warp_stiffness: f64,
    /// Stiffness along weft direction (N/m).
    pub weft_stiffness: f64,
    /// Shear stiffness (N/m).
    pub shear_stiffness: f64,
}
impl FabricWeave {
    /// Create a new fabric weave model.
    pub fn new(
        warp: [f64; 3],
        weft: [f64; 3],
        warp_stiffness: f64,
        weft_stiffness: f64,
        shear_stiffness: f64,
    ) -> Self {
        Self {
            warp: normalize3(warp),
            weft: normalize3(weft),
            warp_stiffness,
            weft_stiffness,
            shear_stiffness,
        }
    }
    /// Effective stiffness in direction `d` (unit vector) using off-axis formula.
    ///
    /// Uses a simple trigonometric blending of warp and weft stiffnesses.
    pub fn stiffness_in_direction(&self, d: [f64; 3]) -> f64 {
        let dn = normalize3(d);
        let cos_warp = dot3(dn, self.warp).abs();
        let cos_weft = dot3(dn, self.weft).abs();
        self.warp_stiffness * cos_warp * cos_warp
            + self.weft_stiffness * cos_weft * cos_weft
            + self.shear_stiffness * (1.0 - cos_warp * cos_warp - cos_weft * cos_weft).abs()
    }
    /// Compute the force on a spring edge given its direction and stretch.
    pub fn edge_force(&self, edge_dir: [f64; 3], stretch: f64) -> f64 {
        let k = self.stiffness_in_direction(edge_dir);
        k * stretch
    }
}
/// Multi-layer cloth simulation with inter-layer friction.
#[derive(Debug, Clone)]
pub struct LayeredCloth {
    /// Particles for each layer.
    pub layers: Vec<Vec<ClothParticle>>,
    /// Inter-layer friction coefficient.
    pub inter_layer_friction: f64,
    /// Air pocket thickness between layers (m).
    pub air_pocket_thickness: f64,
    /// Interlocking spring stiffness.
    pub interlocking_stiffness: f64,
}
impl LayeredCloth {
    /// Create a new layered cloth with `num_layers` of `num_particles` each.
    pub fn new(
        num_layers: usize,
        num_particles: usize,
        mass_per_particle: f64,
        inter_layer_friction: f64,
        air_pocket_thickness: f64,
    ) -> Self {
        let layers = (0..num_layers)
            .map(|l| {
                (0..num_particles)
                    .map(|_| {
                        ClothParticle::new(
                            [0.0, l as f64 * air_pocket_thickness, 0.0],
                            mass_per_particle,
                        )
                    })
                    .collect()
            })
            .collect();
        Self {
            layers,
            inter_layer_friction,
            air_pocket_thickness,
            interlocking_stiffness: 100.0,
        }
    }
    /// Number of layers.
    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }
    /// Apply inter-layer friction forces (simplified, between adjacent layers).
    pub fn apply_inter_layer_friction(&mut self, dt: f64) {
        let nl = self.layers.len();
        if nl < 2 {
            return;
        }
        for layer_idx in 0..(nl - 1) {
            let n = self.layers[layer_idx]
                .len()
                .min(self.layers[layer_idx + 1].len());
            for i in 0..n {
                let v_a = self.layers[layer_idx][i].velocity;
                let v_b = self.layers[layer_idx + 1][i].velocity;
                let dv = [v_b[0] - v_a[0], v_b[1] - v_a[1], v_b[2] - v_a[2]];
                let friction_force = [
                    self.inter_layer_friction * dv[0],
                    self.inter_layer_friction * dv[1],
                    self.inter_layer_friction * dv[2],
                ];
                if !self.layers[layer_idx][i].is_pinned() {
                    let inv_m = self.layers[layer_idx][i].inv_mass;
                    self.layers[layer_idx][i].velocity[0] += friction_force[0] * inv_m * dt;
                    self.layers[layer_idx][i].velocity[1] += friction_force[1] * inv_m * dt;
                    self.layers[layer_idx][i].velocity[2] += friction_force[2] * inv_m * dt;
                }
                if !self.layers[layer_idx + 1][i].is_pinned() {
                    let inv_m = self.layers[layer_idx + 1][i].inv_mass;
                    self.layers[layer_idx + 1][i].velocity[0] -= friction_force[0] * inv_m * dt;
                    self.layers[layer_idx + 1][i].velocity[1] -= friction_force[1] * inv_m * dt;
                    self.layers[layer_idx + 1][i].velocity[2] -= friction_force[2] * inv_m * dt;
                }
            }
        }
    }
    /// Total kinetic energy across all layers.
    pub fn kinetic_energy(&self) -> f64 {
        self.layers
            .iter()
            .flat_map(|l| l.iter())
            .map(|p| {
                let m = if p.inv_mass > 1e-15 {
                    1.0 / p.inv_mass
                } else {
                    0.0
                };
                0.5 * m * dot3(p.velocity, p.velocity)
            })
            .sum()
    }
}
/// Drape simulation parameters for a fabric over an obstacle.
#[derive(Debug, Clone)]
pub struct DrapeSimParams {
    /// Number of drape petals (symmetry order).
    pub petal_count: usize,
    /// Fabric radius (m).
    pub fabric_radius: f64,
    /// Support pedestal radius (m).
    pub pedestal_radius: f64,
    /// Pedestal height (m).
    pub pedestal_height: f64,
    /// Fabric linear density (kg/m²).
    pub areal_density: f64,
    /// Bending rigidity (N·m).
    pub bending_rigidity: f64,
    /// Gravity (m/s²).
    pub gravity: f64,
}
impl DrapeSimParams {
    /// Create with typical muslin parameters.
    pub fn muslin() -> Self {
        Self {
            petal_count: 5,
            fabric_radius: 0.3,
            pedestal_radius: 0.1,
            pedestal_height: 0.1,
            areal_density: 0.12,
            bending_rigidity: 5e-6,
            gravity: 9.81,
        }
    }
    /// Drape coefficient (dimensionless) — fraction of area projected below pedestal.
    ///
    /// High drape coefficient → fabric hangs down; low → stiff.
    pub fn drape_coefficient(&self) -> f64 {
        let gravity_force = self.areal_density * self.gravity;
        let stiffness_char =
            self.bending_rigidity / (gravity_force * self.pedestal_radius * self.pedestal_radius);
        (1.0 - stiffness_char.min(1.0)).max(0.0)
    }
    /// Number of wrinkles/folds (approximated as petal_count for symmetric drape).
    pub fn fold_count(&self) -> usize {
        self.petal_count
    }
    /// Petal amplitude (height of draping fold tip).
    pub fn petal_amplitude(&self) -> f64 {
        let overhang = self.fabric_radius - self.pedestal_radius;
        overhang
            * (1.0
                - self.bending_rigidity
                    / (self.areal_density * self.gravity * overhang * overhang + 1e-30))
                .max(0.0)
    }
}
/// A seam line constraint connecting two panels.
#[derive(Debug, Clone)]
pub struct SeamLine {
    /// Indices into panel A particles.
    pub panel_a_indices: Vec<usize>,
    /// Indices into panel B particles.
    pub panel_b_indices: Vec<usize>,
    /// Seam stiffness.
    pub stiffness: f64,
    /// Current seam gap width (m).
    pub gap: f64,
}
impl SeamLine {
    /// Create a new seam line.
    pub fn new(panel_a_indices: Vec<usize>, panel_b_indices: Vec<usize>, stiffness: f64) -> Self {
        Self {
            panel_a_indices,
            panel_b_indices,
            stiffness,
            gap: 0.0,
        }
    }
    /// Number of seam pairs.
    pub fn num_pairs(&self) -> usize {
        self.panel_a_indices.len().min(self.panel_b_indices.len())
    }
}
/// A garment composed of multiple sewing panels with avatar collision.
#[derive(Debug, Clone)]
pub struct ClothGarment {
    /// Sewing system.
    pub sewing: ClothSewing,
    /// Avatar collision spheres: (centre, radius).
    pub avatar_spheres: Vec<([f64; 3], f64)>,
    /// Garment name.
    pub name: String,
}
impl ClothGarment {
    /// Create a new garment.
    pub fn new(
        name: impl Into<String>,
        num_panels: usize,
        particles_per_panel: usize,
        mass: f64,
    ) -> Self {
        Self {
            sewing: ClothSewing::new(num_panels, particles_per_panel, mass),
            avatar_spheres: Vec::new(),
            name: name.into(),
        }
    }
    /// Add a collision sphere to the avatar.
    pub fn add_avatar_sphere(&mut self, centre: [f64; 3], radius: f64) {
        self.avatar_spheres.push((centre, radius));
    }
    /// Resolve collisions with avatar spheres.
    pub fn resolve_avatar_collisions(&mut self) {
        for p in &mut self.sewing.particles {
            if p.is_pinned() {
                continue;
            }
            for &(centre, radius) in &self.avatar_spheres {
                let diff = [
                    p.position[0] - centre[0],
                    p.position[1] - centre[1],
                    p.position[2] - centre[2],
                ];
                let dist = norm3(diff);
                if dist < radius && dist > 1e-15 {
                    let n = normalize3(diff);
                    let corr = radius - dist;
                    p.position[0] += n[0] * corr;
                    p.position[1] += n[1] * corr;
                    p.position[2] += n[2] * corr;
                }
            }
        }
    }
    /// Total number of particles.
    pub fn num_particles(&self) -> usize {
        self.sewing.particles.len()
    }
}
/// Bernoulli-based wind pressure on cloth.
#[derive(Debug, Clone)]
pub struct WindDraping {
    /// Wind velocity vector (m/s).
    pub wind_velocity: [f64; 3],
    /// Air density (kg/m³).
    pub rho: f64,
    /// Drag coefficient.
    pub cd: f64,
}
impl WindDraping {
    /// Create a new wind model.
    pub fn new(wind_velocity: [f64; 3], rho: f64, cd: f64) -> Self {
        Self {
            wind_velocity,
            rho,
            cd,
        }
    }
    /// Compute wind force on a triangle defined by three particle positions.
    ///
    /// Returns force contributions to each vertex.
    pub fn triangle_force(&self, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> [[f64; 3]; 3] {
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let n_raw = cross3(e1, e2);
        let area = norm3(n_raw) * 0.5;
        if area < 1e-15 {
            return [[0.0; 3]; 3];
        }
        let n = normalize3(n_raw);
        let v_dot_n = dot3(self.wind_velocity, n);
        let dynamic_pressure = 0.5 * self.rho * v_dot_n * v_dot_n;
        let force_mag = self.cd * dynamic_pressure * area;
        let sign = if v_dot_n < 0.0 { -1.0 } else { 1.0 };
        let f = [
            n[0] * force_mag * sign / 3.0,
            n[1] * force_mag * sign / 3.0,
            n[2] * force_mag * sign / 3.0,
        ];
        [f, f, f]
    }
    /// Apply wind forces to all particles given triangle connectivity.
    pub fn apply(&self, particles: &mut [ClothParticle], triangles: &[[usize; 3]]) {
        for tri in triangles {
            let [i0, i1, i2] = *tri;
            if i0 >= particles.len() || i1 >= particles.len() || i2 >= particles.len() {
                continue;
            }
            let p0 = particles[i0].position;
            let p1 = particles[i1].position;
            let p2 = particles[i2].position;
            let forces = self.triangle_force(p0, p1, p2);
            for (pi, fi) in [(i0, forces[0]), (i1, forces[1]), (i2, forces[2])] {
                if !particles[pi].is_pinned() {
                    particles[pi].force[0] += fi[0];
                    particles[pi].force[1] += fi[1];
                    particles[pi].force[2] += fi[2];
                }
            }
        }
    }
}
/// Simulation mesh vertex for rendering.
#[derive(Debug, Clone)]
pub struct RenderVertex {
    /// World-space position.
    pub position: [f32; 3],
    /// Smooth normal.
    pub normal: [f32; 3],
    /// UV texture coordinates.
    pub uv: [f32; 2],
}
/// Woven fabric weave pattern type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WeavePattern {
    /// Plain weave: alternating over/under.
    Plain,
    /// Twill weave: diagonal lines.
    Twill(u32, u32),
    /// Satin weave: long floats.
    Satin(u32),
    /// Leno weave: twisted warp yarns.
    Leno,
}
impl WeavePattern {
    /// Float length (maximum consecutive overs) for weave type.
    pub fn float_length(&self) -> u32 {
        match self {
            WeavePattern::Plain => 1,
            WeavePattern::Twill(w, _) => *w,
            WeavePattern::Satin(n) => *n - 1,
            WeavePattern::Leno => 1,
        }
    }
    /// Weave crimp factor (≈ (float_length + 1) / (2 * repeat)).
    pub fn crimp_factor(&self) -> f64 {
        match self {
            WeavePattern::Plain => 0.5,
            WeavePattern::Twill(w, f) => *w as f64 / (*w + *f) as f64,
            WeavePattern::Satin(n) => (*n as f64 - 1.0) / (*n as f64),
            WeavePattern::Leno => 0.4,
        }
    }
}
/// Garment fitting quality metrics.
#[derive(Debug, Clone, Default)]
pub struct FitMetrics {
    /// Ease (slack) in circumference (m) — positive = comfortable.
    pub ease: f64,
    /// Average strain in fabric at wearing configuration.
    pub avg_strain: f64,
    /// Maximum strain.
    pub max_strain: f64,
    /// Regions with contact pressure > threshold.
    pub pressure_areas: usize,
    /// Overall fit score \[0..1\]: 1 = perfect.
    pub fit_score: f64,
}
impl FitMetrics {
    /// Compute fit score from ease and strain values.
    pub fn compute_score(&mut self) {
        let ease_score = (self.ease / 0.05).clamp(0.0, 1.0);
        let strain_score = 1.0 - (self.avg_strain / 0.3).min(1.0);
        self.fit_score = 0.6 * ease_score + 0.4 * strain_score;
    }
}
/// Cloth interaction controller for user manipulation.
#[derive(Debug, Clone)]
pub struct ClothInteraction {
    /// Currently selected particle index.
    pub selected_particle: Option<usize>,
    /// Target position for grab/pull operations.
    pub target_position: [f64; 3],
    /// Current interaction mode.
    pub mode: InteractionMode,
    /// Pull/grab stiffness.
    pub stiffness: f64,
    /// Particles that have been cut (removed from simulation).
    pub cut_particles: Vec<usize>,
}
impl ClothInteraction {
    /// Create a new cloth interaction controller.
    pub fn new(stiffness: f64) -> Self {
        Self {
            selected_particle: None,
            target_position: [0.0; 3],
            mode: InteractionMode::Grab,
            stiffness,
            cut_particles: Vec::new(),
        }
    }
    /// Select a particle for interaction.
    pub fn select(&mut self, particle_idx: usize, mode: InteractionMode) {
        self.selected_particle = Some(particle_idx);
        self.mode = mode;
    }
    /// Deselect the current particle.
    pub fn deselect(&mut self) {
        self.selected_particle = None;
    }
    /// Apply the interaction force to the selected particle.
    pub fn apply(&self, particles: &mut [ClothParticle], dt: f64) {
        let Some(idx) = self.selected_particle else {
            return;
        };
        if idx >= particles.len() {
            return;
        }
        if particles[idx].is_pinned() {
            return;
        }
        match self.mode {
            InteractionMode::Grab | InteractionMode::Pull => {
                let p = particles[idx].position;
                let diff = [
                    self.target_position[0] - p[0],
                    self.target_position[1] - p[1],
                    self.target_position[2] - p[2],
                ];
                let inv_m = particles[idx].inv_mass;
                for (v, d) in particles[idx].velocity.iter_mut().zip(diff.iter()) {
                    *v += self.stiffness * d * inv_m * dt;
                }
            }
            InteractionMode::Cut => {
                particles[idx].inv_mass = 0.0;
            }
        }
    }
    /// Perform a closest-particle query.
    pub fn find_closest_particle(
        &self,
        particles: &[ClothParticle],
        query: [f64; 3],
    ) -> Option<usize> {
        particles
            .iter()
            .enumerate()
            .filter(|(_, p)| !p.is_pinned())
            .min_by(|(_, pa), (_, pb)| {
                let da = norm3([
                    pa.position[0] - query[0],
                    pa.position[1] - query[1],
                    pa.position[2] - query[2],
                ]);
                let db = norm3([
                    pb.position[0] - query[0],
                    pb.position[1] - query[1],
                    pb.position[2] - query[2],
                ]);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
    }
}
/// Cloth tearing simulation.
#[derive(Debug, Clone)]
pub struct ClothTearing {
    /// Cloth particles.
    pub particles: Vec<ClothParticle>,
    /// Cloth edges with rupture criteria.
    pub edges: Vec<ClothEdge>,
}
impl ClothTearing {
    /// Create a new cloth tearing system.
    pub fn new(particles: Vec<ClothParticle>, edges: Vec<ClothEdge>) -> Self {
        Self { particles, edges }
    }
    /// Update rupture state for all edges and return count of newly torn edges.
    pub fn check_all_ruptures(&mut self) -> usize {
        let particles = self.particles.clone();
        let mut count = 0;
        for e in self.edges.iter_mut() {
            if !e.torn && e.check_rupture(&particles) {
                count += 1;
            }
        }
        count
    }
    /// Number of intact edges.
    pub fn intact_edge_count(&self) -> usize {
        self.edges.iter().filter(|e| !e.torn).count()
    }
    /// Number of torn edges.
    pub fn torn_edge_count(&self) -> usize {
        self.edges.iter().filter(|e| e.torn).count()
    }
}
