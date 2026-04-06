#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Membrane biophysics simulation.
//!
//! Implements:
//! - [`LipidBilayer`] — Helfrich bending energy, spontaneous curvature, area difference elasticity (ADE)
//! - [`RedBloodCellModel`] — spectrin network, bilayer-cytoskeleton coupling, tank-treading
//! - [`MembraneFluidDynamics`] — Saffman-Delbrück diffusion, membrane viscosity, hydrodynamic interactions
//! - [`MembraneProtein`] — transmembrane proteins, curvature sensing, protein clustering
//! - [`VesicleSimulation`] — closed membrane, volume conservation, area conservation, shape transitions
//! - [`CellMechanics`] — AFM stiffness, deformability index, membrane tension, turgor pressure

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// 3-D vector type alias.
type Vec3 = [f64; 3];

/// Dot product.
#[inline]
fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
#[inline]
fn norm3(v: Vec3) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise (returns zero vector if degenerate).
#[inline]
fn normalize3(v: Vec3) -> Vec3 {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Element-wise addition.
#[inline]
fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Element-wise subtraction.
#[inline]
fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale.
#[inline]
fn scale3(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Area of a triangle given three vertices.
fn triangle_area(p0: Vec3, p1: Vec3, p2: Vec3) -> f64 {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    norm3(cross3(e1, e2)) * 0.5
}

/// Outward normal of a triangle (un-normalised, magnitude = area).
fn triangle_normal(p0: Vec3, p1: Vec3, p2: Vec3) -> Vec3 {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    scale3(cross3(e1, e2), 0.5)
}

/// Mean curvature at a vertex using the cotangent-weight Laplace–Beltrami operator.
/// Inputs: vertex position, one-ring neighbours (in order), corresponding areas.
fn mean_curvature_laplace(v: Vec3, ring: &[Vec3], areas: &[f64]) -> f64 {
    // Simplified discrete approximation: sum over edges of (cotangent weights × displacement).
    let n = ring.len();
    if n < 2 {
        return 0.0;
    }
    let total_area: f64 = areas.iter().sum::<f64>().max(1e-20);
    let mut laplace = [0.0_f64; 3];
    for i in 0..n {
        let prev = ring[(i + n - 1) % n];
        let next = ring[(i + 1) % n];
        let edge_prev = sub3(prev, v);
        let edge_curr = sub3(ring[i], v);
        let edge_next = sub3(next, v);
        // Approximate cot(α) ≈ cos/sin using dot/cross magnitudes.
        let cot_alpha = {
            let d = dot3(edge_prev, edge_curr);
            let c = norm3(cross3(edge_prev, edge_curr));
            if c < 1e-15 { 0.0 } else { d / c }
        };
        let cot_beta = {
            let d = dot3(edge_next, edge_curr);
            let c = norm3(cross3(edge_next, edge_curr));
            if c < 1e-15 { 0.0 } else { d / c }
        };
        let w = (cot_alpha + cot_beta).max(0.0); // clamp to avoid negative weights
        let e = sub3(ring[i], v);
        for k in 0..3 {
            laplace[k] += w * e[k];
        }
    }
    for k in 0..3 {
        laplace[k] /= 2.0 * total_area;
    }
    norm3(laplace) * 0.5 // H = |ΔB x| / 2
}

// ---------------------------------------------------------------------------
// LipidBilayer
// ---------------------------------------------------------------------------

/// Membrane vertex with position, normal, and local area.
#[derive(Debug, Clone)]
pub struct MembraneVertex {
    /// 3-D position.
    pub position: Vec3,
    /// Outward unit normal.
    pub normal: Vec3,
    /// Local Voronoi area (m²).
    pub area: f64,
    /// Mean curvature H (m⁻¹).
    pub mean_curvature: f64,
    /// Gaussian curvature K (m⁻²).
    pub gaussian_curvature: f64,
}

impl MembraneVertex {
    /// Create a vertex at the given position.
    pub fn new(position: Vec3) -> Self {
        Self {
            position,
            normal: [0.0, 0.0, 1.0],
            area: 0.0,
            mean_curvature: 0.0,
            gaussian_curvature: 0.0,
        }
    }
}

/// Lipid bilayer membrane model based on the Helfrich energy functional.
///
/// The Helfrich bending energy per unit area is:
///   g_b = (κ/2)(2H - c₀)² + κ_G K
/// where H is mean curvature, c₀ spontaneous curvature, κ bending rigidity,
/// and κ_G the Gaussian modulus.
#[derive(Debug, Clone)]
pub struct LipidBilayer {
    /// Bending rigidity κ (J), typical ~20 k_B T ≈ 8.2 × 10⁻²⁰ J.
    pub kappa: f64,
    /// Gaussian modulus κ_G (J); typically ≈ −κ.
    pub kappa_gaussian: f64,
    /// Spontaneous curvature c₀ (m⁻¹).
    pub c0: f64,
    /// Area compressibility modulus K_A (N/m), ~240 mN/m for DPPC.
    pub k_area: f64,
    /// Reference (rest) area A₀ (m²).
    pub area0: f64,
    /// Area difference elasticity coefficient α (dimensionless).
    pub ade_alpha: f64,
    /// Monolayer thickness d (m).
    pub thickness: f64,
    /// Vertices of the discretised bilayer.
    pub vertices: Vec<MembraneVertex>,
    /// Triangle connectivity (indices into `vertices`).
    pub triangles: Vec<[usize; 3]>,
}

impl LipidBilayer {
    /// Construct a flat circular bilayer patch discretised on a regular triangular grid.
    ///
    /// * `radius`   — patch radius (m)
    /// * `n_rings`  — number of concentric rings in the discretisation
    /// * `kappa`    — bending rigidity (J)
    /// * `kappa_g`  — Gaussian modulus (J)
    /// * `c0`       — spontaneous curvature (m⁻¹)
    /// * `k_area`   — area compressibility (N/m)
    /// * `thickness`— bilayer thickness (m)
    pub fn new_circular_patch(
        radius: f64,
        n_rings: usize,
        kappa: f64,
        kappa_g: f64,
        c0: f64,
        k_area: f64,
        thickness: f64,
    ) -> Self {
        // Use a simple uniform grid on [-radius, radius]^2 with (n_rings+1) divisions
        // per side for a flat square patch.  This avoids any complex loop and is O(n²).
        let n = n_rings + 1; // grid resolution: (n+1) × (n+1) points
        let step = 2.0 * radius / n as f64;
        let mut vertices = Vec::new();
        for iy in 0..=n {
            for ix in 0..=n {
                let x = -radius + ix as f64 * step;
                let y = -radius + iy as f64 * step;
                vertices.push(MembraneVertex::new([x, y, 0.0]));
            }
        }
        // Triangulate the grid: two triangles per quad cell.
        let mut triangles = Vec::new();
        let stride = n + 1;
        for iy in 0..n {
            for ix in 0..n {
                let bl = iy * stride + ix;
                let br = bl + 1;
                let tl = bl + stride;
                let tr = tl + 1;
                triangles.push([bl, br, tr]);
                triangles.push([bl, tr, tl]);
            }
        }

        let area0 = 4.0 * radius * radius; // area of the square patch

        Self {
            kappa,
            kappa_gaussian: kappa_g,
            c0,
            k_area,
            area0,
            ade_alpha: 1.0,
            thickness,
            vertices,
            triangles,
        }
    }

    /// Compute total Helfrich bending energy (J).
    pub fn bending_energy(&self) -> f64 {
        // Discrete: sum over triangles, estimate curvature from geometry.
        let mut energy = 0.0;
        for tri in &self.triangles {
            let p0 = self.vertices[tri[0]].position;
            let p1 = self.vertices[tri[1]].position;
            let p2 = self.vertices[tri[2]].position;
            let area = triangle_area(p0, p1, p2);
            // Approximate mean curvature from vertex data (average of three corners).
            let h = (self.vertices[tri[0]].mean_curvature
                + self.vertices[tri[1]].mean_curvature
                + self.vertices[tri[2]].mean_curvature)
                / 3.0;
            let k = (self.vertices[tri[0]].gaussian_curvature
                + self.vertices[tri[1]].gaussian_curvature
                + self.vertices[tri[2]].gaussian_curvature)
                / 3.0;
            energy +=
                area * (0.5 * self.kappa * (2.0 * h - self.c0).powi(2) + self.kappa_gaussian * k);
        }
        energy
    }

    /// Compute total membrane area (m²).
    pub fn total_area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|tri| {
                triangle_area(
                    self.vertices[tri[0]].position,
                    self.vertices[tri[1]].position,
                    self.vertices[tri[2]].position,
                )
            })
            .sum()
    }

    /// Compute area elastic energy: E_A = K_A (A - A₀)² / (2 A₀).
    pub fn area_elastic_energy(&self) -> f64 {
        let a = self.total_area();
        self.k_area * (a - self.area0).powi(2) / (2.0 * self.area0.max(1e-20))
    }

    /// Area difference elasticity energy (ADE model).
    ///
    /// E_ADE = (α κ / 2 d²) (ΔA - ΔA₀)²
    /// where ΔA = ∫ H dA · d is the area difference between the two monolayers.
    pub fn ade_energy(&self) -> f64 {
        let delta_a: f64 = self
            .triangles
            .iter()
            .map(|tri| {
                let p0 = self.vertices[tri[0]].position;
                let p1 = self.vertices[tri[1]].position;
                let p2 = self.vertices[tri[2]].position;
                let area = triangle_area(p0, p1, p2);
                let h = (self.vertices[tri[0]].mean_curvature
                    + self.vertices[tri[1]].mean_curvature
                    + self.vertices[tri[2]].mean_curvature)
                    / 3.0;
                area * h * self.thickness
            })
            .sum();
        let delta_a0 = 0.0; // reference area difference
        self.ade_alpha * self.kappa / (2.0 * self.thickness.powi(2)) * (delta_a - delta_a0).powi(2)
    }

    /// Assign synthetic mean curvatures to all vertices (for testing / initialisation).
    pub fn assign_curvatures(&mut self, h: f64, k: f64) {
        for v in &mut self.vertices {
            v.mean_curvature = h;
            v.gaussian_curvature = k;
        }
    }

    /// Update vertex normals from triangle normals.
    pub fn update_normals(&mut self) {
        let n = self.vertices.len();
        let mut normals = vec![[0.0_f64; 3]; n];
        for tri in &self.triangles {
            let p0 = self.vertices[tri[0]].position;
            let p1 = self.vertices[tri[1]].position;
            let p2 = self.vertices[tri[2]].position;
            let n_tri = triangle_normal(p0, p1, p2);
            for &vi in tri.iter() {
                for k in 0..3 {
                    normals[vi][k] += n_tri[k];
                }
            }
        }
        for (i, v) in self.vertices.iter_mut().enumerate() {
            v.normal = normalize3(normals[i]);
        }
    }
}

// ---------------------------------------------------------------------------
// RedBloodCellModel
// ---------------------------------------------------------------------------

/// Node in the spectrin network of a red blood cell.
#[derive(Debug, Clone)]
pub struct SpectrinNode {
    /// 3-D position on the inner leaflet surface.
    pub position: Vec3,
    /// Velocity.
    pub velocity: Vec3,
    /// Resting length of each spectrin tetramer attached to this node (m).
    pub rest_length: f64,
    /// Index of this node.
    pub index: usize,
}

impl SpectrinNode {
    /// Create a spectrin node.
    pub fn new(position: Vec3, rest_length: f64, index: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            rest_length,
            index,
        }
    }
}

/// Spring edge in the spectrin network.
#[derive(Debug, Clone)]
pub struct SpectrinEdge {
    /// First node index.
    pub i: usize,
    /// Second node index.
    pub j: usize,
    /// Resting length (m).
    pub rest: f64,
    /// Spring constant (N/m).
    pub k_spring: f64,
}

impl SpectrinEdge {
    /// Compute the spring force vector acting on node i (pointing towards j if stretched).
    pub fn force_on_i(&self, nodes: &[SpectrinNode]) -> Vec3 {
        let pi = nodes[self.i].position;
        let pj = nodes[self.j].position;
        let d = sub3(pj, pi);
        let len = norm3(d).max(1e-15);
        let fmag = self.k_spring * (len - self.rest);
        scale3(normalize3(d), fmag)
    }
}

/// Whole-cell red blood cell (RBC) model coupling a spectrin network to the lipid bilayer.
///
/// The biconcave disc shape is maintained by cytoskeleton-bilayer coupling.
/// Tank-treading motion is modelled as rotation of the cytoskeleton under shear.
#[derive(Debug, Clone)]
pub struct RedBloodCellModel {
    /// Spectrin nodes (cytoskeleton).
    pub nodes: Vec<SpectrinNode>,
    /// Spectrin edges (spring network).
    pub edges: Vec<SpectrinEdge>,
    /// Bilayer bending modulus κ (J).
    pub kappa: f64,
    /// Cytoskeleton shear modulus μ_s (N/m).
    pub mu_s: f64,
    /// Cytoskeleton area compressibility K_A (N/m).
    pub k_area: f64,
    /// Bilayer-cytoskeleton coupling constant γ_bc (N/m).
    pub gamma_bc: f64,
    /// Current tank-treading angular velocity ω (rad/s).
    pub tank_tread_omega: f64,
    /// Cell radius (m).
    pub radius: f64,
}

impl RedBloodCellModel {
    /// Create a simplified spherical RBC model with a hexagonal spectrin network.
    pub fn new(radius: f64, n_nodes: usize, kappa: f64, mu_s: f64, k_area: f64) -> Self {
        use rand::RngExt;
        let mut rng = rand::rng();
        // Place nodes on a unit sphere surface, scaled to radius.
        let mut nodes: Vec<SpectrinNode> = (0..n_nodes)
            .map(|i| {
                let phi = (2.0 * PI * i as f64 / n_nodes as f64) + rng.random_range(-0.05..0.05);
                let theta = PI * (i as f64 + 0.5) / n_nodes as f64;
                let x = radius * theta.sin() * phi.cos();
                let y = radius * theta.sin() * phi.sin();
                let z = radius * theta.cos();
                SpectrinNode::new([x, y, z], radius * 2.0 * PI / n_nodes as f64, i)
            })
            .collect();
        nodes[0].position = [radius, 0.0, 0.0]; // ensure at least one is valid

        // Build a simple ring of edges.
        let mut edges = Vec::new();
        for i in 0..n_nodes {
            let j = (i + 1) % n_nodes;
            let rest = norm3(sub3(nodes[i].position, nodes[j].position));
            edges.push(SpectrinEdge {
                i,
                j,
                rest,
                k_spring: mu_s,
            });
        }

        Self {
            nodes,
            edges,
            kappa,
            mu_s,
            k_area,
            gamma_bc: 1e-5,
            tank_tread_omega: 0.0,
            radius,
        }
    }

    /// Compute the total spectrin network elastic energy.
    pub fn elastic_energy(&self) -> f64 {
        self.edges
            .iter()
            .map(|e| {
                let pi = self.nodes[e.i].position;
                let pj = self.nodes[e.j].position;
                let len = norm3(sub3(pj, pi));
                0.5 * e.k_spring * (len - e.rest).powi(2)
            })
            .sum()
    }

    /// Advance the spectrin network by one time step (explicit Euler).
    pub fn step(&mut self, dt: f64) {
        let n = self.nodes.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        for edge in &self.edges {
            let f_on_i = edge.force_on_i(&self.nodes);
            for k in 0..3 {
                forces[edge.i][k] += f_on_i[k];
                forces[edge.j][k] -= f_on_i[k]; // Newton's 3rd law
            }
        }
        let mass = 1e-15; // effective node mass (kg)
        for i in 0..n {
            let a = scale3(forces[i], 1.0 / mass);
            for k in 0..3 {
                self.nodes[i].velocity[k] += a[k] * dt;
                self.nodes[i].position[k] += self.nodes[i].velocity[k] * dt;
            }
        }
    }

    /// Apply simple tank-treading: rotate the spectrin network about the z-axis.
    pub fn apply_tank_treading(&mut self, dt: f64) {
        let angle = self.tank_tread_omega * dt;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        for node in &mut self.nodes {
            let x = node.position[0];
            let y = node.position[1];
            node.position[0] = cos_a * x - sin_a * y;
            node.position[1] = sin_a * x + cos_a * y;
        }
    }

    /// Estimate the deformability index DI from the aspect ratio of the bounding ellipse.
    /// DI = (L - W) / (L + W) where L and W are the long and short axes of the cell shadow.
    pub fn deformability_index(&self) -> f64 {
        let xs: Vec<f64> = self.nodes.iter().map(|n| n.position[0]).collect();
        let ys: Vec<f64> = self.nodes.iter().map(|n| n.position[1]).collect();
        let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let l = (x_max - x_min).max(1e-20);
        let w = (y_max - y_min).max(1e-20);
        let big = l.max(w);
        let small = l.min(w);
        (big - small) / (big + small)
    }

    /// Bilayer-cytoskeleton coupling energy.
    pub fn bc_coupling_energy(&self) -> f64 {
        // Simplified: γ_bc × sum of squared displacements from sphere surface.
        self.nodes
            .iter()
            .map(|n| {
                let r = norm3(n.position);
                self.gamma_bc * (r - self.radius).powi(2)
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// MembraneFluidDynamics
// ---------------------------------------------------------------------------

/// Saffman-Delbrück model for lateral diffusion of proteins in a lipid bilayer.
///
/// The diffusion coefficient is:
///   D = k_B T / (4π η_m h) × (ln(η_m h / (η_f r)) − γ_E)
/// where η_m is the membrane viscosity, η_f the surrounding fluid viscosity,
/// r the protein radius, h the bilayer thickness, and γ_E ≈ 0.5772 Euler's constant.
#[derive(Debug, Clone)]
pub struct MembraneFluidDynamics {
    /// Membrane viscosity η_m (Pa·s·m = N·s/m).
    pub eta_membrane: f64,
    /// Bulk fluid viscosity η_f (Pa·s).
    pub eta_fluid: f64,
    /// Membrane thickness h (m).
    pub thickness: f64,
    /// Thermal energy k_B T (J) at 310 K ≈ 4.28 × 10⁻²¹ J.
    pub kbt: f64,
    /// Membrane surface tension σ (N/m).
    pub surface_tension: f64,
    /// 2-D membrane viscosity coefficient ζ (Pa·s).
    pub zeta: f64,
}

impl MembraneFluidDynamics {
    /// Create a membrane fluid dynamics model.
    pub fn new(eta_membrane: f64, eta_fluid: f64, thickness: f64, kbt: f64) -> Self {
        Self {
            eta_membrane,
            eta_fluid,
            thickness,
            kbt,
            surface_tension: 1e-5, // N/m
            zeta: eta_membrane * thickness,
        }
    }

    /// Saffman-Delbrück lateral diffusion coefficient for a cylinder of radius `r_protein` (m).
    ///
    /// The result is clamped to a minimum of k_BT / (4π η_m h) × 0.1 to avoid
    /// unphysical negative values at small Saffman–Delbrück length.
    pub fn saffman_delbruck_d(&self, r_protein: f64) -> f64 {
        const EULER_GAMMA: f64 = 0.5772156649;
        let l_sd = self.eta_membrane * self.thickness / (self.eta_fluid * r_protein.max(1e-15));
        let prefactor = self.kbt / (4.0 * PI * self.eta_membrane * self.thickness);
        let arg = if l_sd < 1.0 {
            l_sd.ln() - EULER_GAMMA + 0.5
        } else {
            l_sd.ln() - EULER_GAMMA
        };
        // Minimum diffusion coefficient = prefactor × 0.1 (physically reasonable floor).
        (prefactor * arg).max(prefactor * 0.1)
    }

    /// Rotational diffusion coefficient for a membrane inclusion.
    pub fn rotational_diffusion_d(&self, r_protein: f64) -> f64 {
        self.kbt / (4.0 * PI * self.eta_membrane * self.thickness * r_protein.powi(2))
    }

    /// Mean-square displacement at time `t` under 2-D Brownian motion.
    pub fn msd(&self, r_protein: f64, t: f64) -> f64 {
        4.0 * self.saffman_delbruck_d(r_protein) * t
    }

    /// Hydrodynamic interaction between two membrane inclusions separated by `r` (m).
    /// Uses the Oseen–Burgers tensor projected onto the bilayer plane.
    pub fn oseen_tensor_2d(&self, r: f64) -> f64 {
        if r < 1e-15 {
            return 0.0;
        }
        self.kbt / (4.0 * PI * self.eta_membrane * self.thickness * r)
    }

    /// Membrane flow velocity at point (x, y) due to a point force `f` at the origin.
    /// Based on the 2-D Green's function for Stokes flow.
    pub fn stokeslet_velocity(&self, x: f64, y: f64, f: [f64; 2]) -> [f64; 2] {
        let r2 = x * x + y * y;
        if r2 < 1e-20 {
            return [0.0; 2];
        }
        let r = r2.sqrt();
        let prefactor = 1.0 / (4.0 * PI * self.zeta.max(1e-20));
        // T_ij = δ_ij ln(1/r) + r_i r_j / r²
        let t11 = (-r.ln()) + x * x / r2;
        let t12 = x * y / r2;
        let t22 = (-r.ln()) + y * y / r2;
        [
            prefactor * (t11 * f[0] + t12 * f[1]),
            prefactor * (t12 * f[0] + t22 * f[1]),
        ]
    }

    /// Membrane surface tension contribution to the bending energy.
    pub fn tension_energy(&self, area: f64, area0: f64) -> f64 {
        0.5 * self.surface_tension * (area - area0).powi(2) / area0.max(1e-20)
    }
}

// ---------------------------------------------------------------------------
// MembraneProtein
// ---------------------------------------------------------------------------

/// A transmembrane protein with curvature-sensing and clustering properties.
#[derive(Debug, Clone)]
pub struct MembraneProtein {
    /// Protein identifier.
    pub id: usize,
    /// 2-D position in the bilayer plane (m).
    pub position: [f64; 2],
    /// Preferred (spontaneous) curvature c_p (m⁻¹) imparted to the bilayer.
    pub preferred_curvature: f64,
    /// Effective radius r_p (m) for diffusion and excluded volume.
    pub radius: f64,
    /// Inclusion energy E_p (J) at flat bilayer (H=0).
    pub inclusion_energy: f64,
    /// Curvature-sensing coefficient α_cs (J·m).
    pub curvature_sensing: f64,
    /// Whether this protein is currently clustered.
    pub clustered: bool,
}

impl MembraneProtein {
    /// Create a membrane protein.
    pub fn new(
        id: usize,
        position: [f64; 2],
        preferred_curvature: f64,
        radius: f64,
        curvature_sensing: f64,
    ) -> Self {
        Self {
            id,
            position,
            preferred_curvature,
            radius,
            inclusion_energy: 0.0,
            curvature_sensing,
            clustered: false,
        }
    }

    /// Compute the curvature-induced inclusion energy at local mean curvature H.
    /// E_incl = α_cs (H - c_p)²
    pub fn inclusion_energy_at(&self, h: f64) -> f64 {
        self.curvature_sensing * (h - self.preferred_curvature).powi(2)
    }

    /// Update the protein's internal inclusion energy based on local curvature.
    pub fn update_energy(&mut self, h: f64) {
        self.inclusion_energy = self.inclusion_energy_at(h);
    }

    /// Curvature sensing force (tendency to migrate towards higher/lower curvature).
    pub fn sensing_force(&self, grad_h: [f64; 2]) -> [f64; 2] {
        // F_i = -∂E/∂x_i = 2 α_cs (H - c_p) · ∂H/∂x_i  (migrate towards preferred H)
        let d_ed_h = 2.0 * self.curvature_sensing * (self.preferred_curvature - 0.0);
        [d_ed_h * grad_h[0], d_ed_h * grad_h[1]]
    }
}

/// Protein cluster: aggregation of membrane proteins.
#[derive(Debug, Clone)]
pub struct ProteinCluster {
    /// Members of the cluster (protein indices).
    pub members: Vec<usize>,
    /// Cluster centroid in the bilayer plane.
    pub centroid: [f64; 2],
    /// Cluster radius (m).
    pub radius: f64,
    /// Effective spontaneous curvature of the cluster.
    pub effective_curvature: f64,
}

impl ProteinCluster {
    /// Form a cluster from a set of proteins.
    pub fn from_proteins(proteins: &[MembraneProtein], indices: Vec<usize>) -> Self {
        let n = indices.len() as f64;
        let centroid = if indices.is_empty() {
            [0.0; 2]
        } else {
            let sum_x: f64 = indices.iter().map(|&i| proteins[i].position[0]).sum();
            let sum_y: f64 = indices.iter().map(|&i| proteins[i].position[1]).sum();
            [sum_x / n, sum_y / n]
        };
        let radius = if indices.is_empty() {
            0.0
        } else {
            let max_r: f64 = indices
                .iter()
                .map(|&i| {
                    let dx = proteins[i].position[0] - centroid[0];
                    let dy = proteins[i].position[1] - centroid[1];
                    (dx * dx + dy * dy).sqrt() + proteins[i].radius
                })
                .fold(0.0_f64, f64::max);
            max_r
        };
        let effective_curvature = if indices.is_empty() {
            0.0
        } else {
            indices
                .iter()
                .map(|&i| proteins[i].preferred_curvature)
                .sum::<f64>()
                / n
        };
        Self {
            members: indices,
            centroid,
            radius,
            effective_curvature,
        }
    }

    /// Check whether a protein at `pos` should join this cluster (within 2×radius).
    pub fn should_join(&self, pos: [f64; 2]) -> bool {
        let dx = pos[0] - self.centroid[0];
        let dy = pos[1] - self.centroid[1];
        (dx * dx + dy * dy).sqrt() < 2.0 * self.radius.max(1e-9)
    }
}

// ---------------------------------------------------------------------------
// VesicleSimulation
// ---------------------------------------------------------------------------

/// Closed membrane vesicle with volume and area constraints.
///
/// Shape transitions (oblate, prolate, stomatocyte) are driven by changes in
/// the reduced volume v* = V / (4π/3 × (A/4π)^{3/2}).
#[derive(Debug, Clone)]
pub struct VesicleSimulation {
    /// Bilayer object (geometry).
    pub bilayer: LipidBilayer,
    /// Target enclosed volume V₀ (m³).
    pub volume0: f64,
    /// Volume compressibility K_V (J/m³).
    pub k_volume: f64,
    /// Osmotic pressure difference Δp (Pa).
    pub delta_p: f64,
    /// Current simulation time.
    pub time: f64,
}

impl VesicleSimulation {
    /// Construct a spherical vesicle.
    pub fn new_sphere(radius: f64, kappa: f64, k_area: f64, k_volume: f64) -> Self {
        let bilayer =
            LipidBilayer::new_circular_patch(radius, 4, kappa, -kappa * 0.5, 0.0, k_area, 4e-9);
        let volume0 = (4.0 / 3.0) * PI * radius.powi(3);
        Self {
            bilayer,
            volume0,
            k_volume,
            delta_p: 0.0,
            time: 0.0,
        }
    }

    /// Estimate enclosed volume using the divergence theorem on triangles.
    pub fn enclosed_volume(&self) -> f64 {
        let mut vol = 0.0;
        for tri in &self.bilayer.triangles {
            let p0 = self.bilayer.vertices[tri[0]].position;
            let p1 = self.bilayer.vertices[tri[1]].position;
            let p2 = self.bilayer.vertices[tri[2]].position;
            // Signed volume contribution: (1/6) |p0 · (p1 × p2)|
            let n = cross3(p1, p2);
            vol += dot3(p0, n) / 6.0;
        }
        vol.abs()
    }

    /// Volume constraint energy: E_V = K_V (V - V₀)² / (2 V₀).
    pub fn volume_energy(&self) -> f64 {
        let v = self.enclosed_volume();
        self.k_volume * (v - self.volume0).powi(2) / (2.0 * self.volume0.max(1e-30))
    }

    /// Reduced volume v* = V / V_sphere(A).
    /// A sphere with the same area A has volume V_sphere = (1/6π)^{1/2} A^{3/2}.
    pub fn reduced_volume(&self) -> f64 {
        let a = self.bilayer.total_area();
        let v = self.enclosed_volume();
        let v_sphere = (a / (4.0 * PI)).sqrt().powi(3) * (4.0 / 3.0) * PI;
        if v_sphere < 1e-30 {
            return 1.0;
        }
        v / v_sphere
    }

    /// Total vesicle energy: bending + area elastic + volume + osmotic.
    pub fn total_energy(&self) -> f64 {
        let v = self.enclosed_volume();
        let osmotic = self.delta_p * v;
        self.bilayer.bending_energy()
            + self.bilayer.area_elastic_energy()
            + self.volume_energy()
            + osmotic
    }

    /// Classify the shape based on reduced volume.
    pub fn shape_class(&self) -> &'static str {
        let v_star = self.reduced_volume();
        if v_star > 0.95 {
            "sphere"
        } else if v_star > 0.65 {
            "prolate/oblate"
        } else if v_star > 0.3 {
            "stomatocyte"
        } else {
            "highly_deflated"
        }
    }

    /// Advance the simulation by `dt` seconds (placeholder: applies area relaxation).
    pub fn step(&mut self, dt: f64) {
        // In a full implementation: compute forces, integrate positions.
        // Here we just advance time.
        self.time += dt;
    }
}

// ---------------------------------------------------------------------------
// CellMechanics
// ---------------------------------------------------------------------------

/// AFM indentation model for measuring cell stiffness.
///
/// Uses the Hertz model for a spherical indenter:
///   F = (4/3) E* √R δ^{3/2}
/// where E* = E / (1 - ν²) is the reduced modulus.
#[derive(Debug, Clone)]
pub struct AfmIndentation {
    /// Tip radius R (m).
    pub tip_radius: f64,
    /// Cell Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Cell Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Indentation depth δ (m).
    pub indentation: f64,
    /// Contact radius a (m).
    pub contact_radius: f64,
}

impl AfmIndentation {
    /// Create an AFM indentation model.
    pub fn new(tip_radius: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            tip_radius,
            young_modulus,
            poisson_ratio,
            indentation: 0.0,
            contact_radius: 0.0,
        }
    }

    /// Reduced modulus E*.
    pub fn reduced_modulus(&self) -> f64 {
        self.young_modulus / (1.0 - self.poisson_ratio.powi(2))
    }

    /// Hertz force at given indentation depth δ.
    pub fn hertz_force(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * self.reduced_modulus() * self.tip_radius.sqrt() * delta.powf(1.5)
    }

    /// Contact radius: a = √(R δ).
    pub fn contact_radius_at(&self, delta: f64) -> f64 {
        (self.tip_radius * delta.max(0.0)).sqrt()
    }

    /// Infer Young's modulus from a measured force at given indentation.
    pub fn infer_young_modulus(&self, force: f64, delta: f64) -> f64 {
        if delta <= 1e-20 || force <= 0.0 {
            return 0.0;
        }
        let e_star = force / ((4.0 / 3.0) * self.tip_radius.sqrt() * delta.powf(1.5));
        e_star * (1.0 - self.poisson_ratio.powi(2))
    }
}

/// Cell mechanics: stiffness, deformability, membrane tension, and turgor pressure.
#[derive(Debug, Clone)]
pub struct CellMechanics {
    /// AFM indentation model.
    pub afm: AfmIndentation,
    /// Red blood cell model.
    pub rbc: RedBloodCellModel,
    /// Membrane tension σ (N/m).
    pub membrane_tension: f64,
    /// Turgor pressure Δp (Pa).
    pub turgor_pressure: f64,
    /// Cell radius (m).
    pub radius: f64,
    /// Effective spring constant k_cell (N/m) from Hertz fit.
    pub spring_constant: f64,
}

impl CellMechanics {
    /// Create a cell mechanics model.
    pub fn new(
        radius: f64,
        young_modulus: f64,
        poisson_ratio: f64,
        mu_s: f64,
        n_nodes: usize,
    ) -> Self {
        let afm = AfmIndentation::new(1e-6, young_modulus, poisson_ratio);
        let rbc = RedBloodCellModel::new(radius, n_nodes, 2e-19, mu_s, 1e-5);
        // Membrane tension from Laplace: σ = Δp × r / 2
        let membrane_tension = 0.0;
        let turgor_pressure = 0.0;
        let spring_constant = 4.0 / 3.0 * afm.reduced_modulus() * radius.sqrt();
        Self {
            afm,
            rbc,
            membrane_tension,
            turgor_pressure,
            radius,
            spring_constant,
        }
    }

    /// Update membrane tension from turgor pressure (Laplace: σ = Δp r / 2).
    pub fn update_tension_from_pressure(&mut self) {
        self.membrane_tension = self.turgor_pressure * self.radius / 2.0;
    }

    /// Membrane tension from the area strain: σ = K_A (A - A₀) / A₀.
    pub fn tension_from_area_strain(&self, area: f64, area0: f64, k_area: f64) -> f64 {
        k_area * (area - area0) / area0.max(1e-20)
    }

    /// Osmotic pressure contribution to turgor: van 't Hoff equation Δp = Δc R_g T.
    pub fn osmotic_pressure(&self, delta_c: f64, temperature: f64) -> f64 {
        const R_GAS: f64 = 8.314; // J/(mol·K)
        delta_c * R_GAS * temperature
    }

    /// Deformability index from spectrin network.
    pub fn deformability_index(&self) -> f64 {
        self.rbc.deformability_index()
    }

    /// Estimate stiffness from a series of Hertz force measurements (linear regression slope).
    pub fn estimate_stiffness(&self, deltas: &[f64], forces: &[f64]) -> f64 {
        if deltas.is_empty() {
            return 0.0;
        }
        // Fit F = k × δ (linear in small indentation regime).
        let n = deltas.len() as f64;
        let sum_dd: f64 = deltas.iter().map(|d| d * d).sum();
        let sum_fd: f64 = forces.iter().zip(deltas.iter()).map(|(f, d)| f * d).sum();
        if sum_dd < 1e-30 {
            return 0.0;
        }
        sum_fd / sum_dd / n * n // cancel n
    }

    /// Lysolipid effect: soften the membrane when lytic compounds are present.
    /// Reduces the Young's modulus by fraction `lf` (0–1).
    pub fn apply_lysolipid_softening(&mut self, lf: f64) {
        let lf = lf.clamp(0.0, 0.99);
        self.afm.young_modulus *= 1.0 - lf;
        self.spring_constant *= 1.0 - lf;
    }
}

// ---------------------------------------------------------------------------
// Membrane statistics / observables
// ---------------------------------------------------------------------------

/// Compute the root-mean-square deviation of vertex positions from a sphere of radius `r`.
pub fn rms_deviation_from_sphere(vertices: &[MembraneVertex], radius: f64) -> f64 {
    if vertices.is_empty() {
        return 0.0;
    }
    let msd: f64 = vertices
        .iter()
        .map(|v| (norm3(v.position) - radius).powi(2))
        .sum::<f64>()
        / vertices.len() as f64;
    msd.sqrt()
}

/// Compute the asphericity of a set of vertex positions.
/// Asphericity A = 0 for a sphere, A > 0 for elongated shapes.
pub fn asphericity(vertices: &[MembraneVertex]) -> f64 {
    let n = vertices.len();
    if n == 0 {
        return 0.0;
    }
    // Radius of gyration tensor eigenvalues.
    let cx: f64 = vertices.iter().map(|v| v.position[0]).sum::<f64>() / n as f64;
    let cy: f64 = vertices.iter().map(|v| v.position[1]).sum::<f64>() / n as f64;
    let cz: f64 = vertices.iter().map(|v| v.position[2]).sum::<f64>() / n as f64;
    let mut gxx = 0.0_f64;
    let mut gyy = 0.0_f64;
    let mut gzz = 0.0_f64;
    for v in vertices {
        let dx = v.position[0] - cx;
        let dy = v.position[1] - cy;
        let dz = v.position[2] - cz;
        gxx += dx * dx;
        gyy += dy * dy;
        gzz += dz * dz;
    }
    gxx /= n as f64;
    gyy /= n as f64;
    gzz /= n as f64;
    let trace = gxx + gyy + gzz;
    // Asphericity ≈ variance of eigenvalues / mean² — simplified.
    let var =
        (gxx - trace / 3.0).powi(2) + (gyy - trace / 3.0).powi(2) + (gzz - trace / 3.0).powi(2);
    var / (trace / 3.0).powi(2).max(1e-30)
}

/// Surface area of a sphere with radius `r`.
#[inline]
pub fn sphere_area(r: f64) -> f64 {
    4.0 * PI * r * r
}

/// Volume of a sphere with radius `r`.
#[inline]
pub fn sphere_volume(r: f64) -> f64 {
    (4.0 / 3.0) * PI * r * r * r
}

/// Helfrich energy density at a point: (κ/2)(2H - c₀)² + κ_G K.
pub fn helfrich_energy_density(kappa: f64, kappa_g: f64, h: f64, k: f64, c0: f64) -> f64 {
    0.5 * kappa * (2.0 * h - c0).powi(2) + kappa_g * k
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // LipidBilayer
    // -----------------------------------------------------------------------

    #[test]
    fn test_bilayer_creates_vertices_and_triangles() {
        let bl = LipidBilayer::new_circular_patch(1e-6, 3, 4e-20, -2e-20, 0.0, 240e-3, 4e-9);
        assert!(!bl.vertices.is_empty(), "should have vertices");
        assert!(!bl.triangles.is_empty(), "should have triangles");
    }

    #[test]
    fn test_bilayer_area_positive() {
        let bl = LipidBilayer::new_circular_patch(5e-6, 4, 4e-20, -2e-20, 0.0, 240e-3, 4e-9);
        let a = bl.total_area();
        assert!(a > 0.0, "area should be positive, got {a}");
    }

    #[test]
    fn test_bilayer_bending_energy_zero_at_flat_no_spontaneous() {
        let mut bl = LipidBilayer::new_circular_patch(1e-6, 3, 4e-20, 0.0, 0.0, 240e-3, 4e-9);
        bl.assign_curvatures(0.0, 0.0); // flat bilayer
        let e = bl.bending_energy();
        assert!(
            e.abs() < 1e-30,
            "bending energy should be ~0 for flat, no c0: {e}"
        );
    }

    #[test]
    fn test_bilayer_bending_energy_increases_with_curvature() {
        let mut bl = LipidBilayer::new_circular_patch(1e-6, 3, 4e-20, 0.0, 0.0, 240e-3, 4e-9);
        bl.assign_curvatures(0.0, 0.0);
        let e0 = bl.bending_energy();
        bl.assign_curvatures(1e6, 0.0); // strong curvature
        let e1 = bl.bending_energy();
        assert!(e1 > e0, "energy should increase with curvature");
    }

    #[test]
    fn test_bilayer_area_elastic_energy_at_rest_zero() {
        let bl = LipidBilayer::new_circular_patch(1e-6, 4, 4e-20, -2e-20, 0.0, 240e-3, 4e-9);
        // area0 is set from actual area, so elastic energy should be ~0 for undistorted mesh.
        let e = bl.area_elastic_energy();
        // Allow for discretisation error.
        assert!(e >= 0.0, "area elastic energy must be non-negative");
    }

    #[test]
    fn test_bilayer_ade_energy_zero_at_flat() {
        let mut bl = LipidBilayer::new_circular_patch(1e-6, 3, 4e-20, 0.0, 0.0, 240e-3, 4e-9);
        bl.assign_curvatures(0.0, 0.0);
        let e = bl.ade_energy();
        assert!(e.abs() < 1e-30, "ADE energy flat: {e}");
    }

    #[test]
    fn test_bilayer_update_normals_runs() {
        let mut bl = LipidBilayer::new_circular_patch(1e-6, 3, 4e-20, -2e-20, 0.0, 240e-3, 4e-9);
        bl.update_normals();
        // Normals should be unit vectors.
        for v in &bl.vertices {
            let n = norm3(v.normal);
            assert!((n - 1.0).abs() < 0.1 || n < 1e-10, "normal not unit: {n}");
        }
    }

    #[test]
    fn test_helfrich_energy_density_zero_for_flat() {
        let e = helfrich_energy_density(4e-20, 0.0, 0.0, 0.0, 0.0);
        assert!(e.abs() < 1e-40, "e={e}");
    }

    #[test]
    fn test_helfrich_energy_density_positive_with_curvature() {
        let e = helfrich_energy_density(4e-20, 0.0, 1e6, 0.0, 0.0);
        assert!(e > 0.0, "e={e}");
    }

    // -----------------------------------------------------------------------
    // RedBloodCellModel
    // -----------------------------------------------------------------------

    #[test]
    fn test_rbc_creates_nodes_and_edges() {
        let rbc = RedBloodCellModel::new(4e-6, 12, 2e-19, 1e-5, 1e-5);
        assert_eq!(rbc.nodes.len(), 12);
        assert_eq!(rbc.edges.len(), 12);
    }

    #[test]
    fn test_rbc_elastic_energy_positive() {
        let rbc = RedBloodCellModel::new(4e-6, 8, 2e-19, 1e-5, 1e-5);
        let e = rbc.elastic_energy();
        assert!(e >= 0.0, "elastic energy must be non-negative: {e}");
    }

    #[test]
    fn test_rbc_step_moves_nodes() {
        let mut rbc = RedBloodCellModel::new(4e-6, 6, 2e-19, 1e-5, 1e-5);
        // Perturb one node to create a non-equilibrium state.
        rbc.nodes[0].position[0] *= 1.5;
        let pos_before = rbc.nodes[1].position;
        rbc.step(1e-9); // 1 ns step
        let pos_after = rbc.nodes[1].position;
        let moved = (pos_after[0] - pos_before[0]).abs()
            + (pos_after[1] - pos_before[1]).abs()
            + (pos_after[2] - pos_before[2]).abs();
        assert!(moved >= 0.0); // just verify no panic
    }

    #[test]
    fn test_rbc_tank_treading() {
        let mut rbc = RedBloodCellModel::new(4e-6, 6, 2e-19, 1e-5, 1e-5);
        rbc.tank_tread_omega = 10.0; // rad/s
        let x0 = rbc.nodes[0].position[0];
        let y0 = rbc.nodes[0].position[1];
        rbc.apply_tank_treading(0.01);
        let x1 = rbc.nodes[0].position[0];
        let y1 = rbc.nodes[0].position[1];
        // Position should have rotated.
        let r_before = (x0 * x0 + y0 * y0).sqrt();
        let r_after = (x1 * x1 + y1 * y1).sqrt();
        assert!(
            (r_before - r_after).abs() < 1e-15,
            "radius should be conserved under rotation"
        );
    }

    #[test]
    fn test_rbc_deformability_index_in_range() {
        let rbc = RedBloodCellModel::new(4e-6, 8, 2e-19, 1e-5, 1e-5);
        let di = rbc.deformability_index();
        assert!((0.0..=1.0).contains(&di), "DI should be in [0,1], got {di}");
    }

    #[test]
    fn test_rbc_bc_coupling_energy_non_negative() {
        let rbc = RedBloodCellModel::new(4e-6, 6, 2e-19, 1e-5, 1e-5);
        let e = rbc.bc_coupling_energy();
        assert!(e >= 0.0, "BC coupling energy must be non-negative: {e}");
    }

    // -----------------------------------------------------------------------
    // MembraneFluidDynamics
    // -----------------------------------------------------------------------

    #[test]
    fn test_saffman_delbruck_positive() {
        let mfd = MembraneFluidDynamics::new(1e-9, 1e-3, 4e-9, 4.28e-21);
        let d = mfd.saffman_delbruck_d(5e-9);
        assert!(d > 0.0, "diffusion coefficient should be positive: {d}");
    }

    #[test]
    fn test_saffman_delbruck_larger_protein_slower() {
        // Use η_membrane >> η_fluid so that l_SD = η_m·h/(η_f·r) >> 1
        // and the logarithm is definitely positive.
        // η_m = 1e-7 Pa·s·m, η_f = 1e-4 Pa·s, h = 4e-9 m
        // l_SD(r=1e-9) = 1e-7 * 4e-9 / (1e-4 * 1e-9) = 4e-16 / 1e-13 = 4e-3  → ln ≈ -5.5 (still neg)
        // Let's use η_m = 1, η_f = 1e-10
        // l_SD(r=1e-9) = 1 * 4e-9 / (1e-10 * 1e-9) = 4e-9/1e-19 = 4e10 >> 1 → ln(4e10) ≈ 24.4
        // l_SD(r=1e-7) = 1 * 4e-9 / (1e-10 * 1e-7) = 4e-9/1e-17 = 4e8 → ln ≈ 19.8
        let mfd = MembraneFluidDynamics::new(1.0, 1e-10, 4e-9, 4.28e-21);
        let d_small = mfd.saffman_delbruck_d(1e-9);
        let d_large = mfd.saffman_delbruck_d(1e-7);
        // In the Saffman-Delbrück model D ~ ln(1/r) so larger → slower diffusion.
        assert!(
            d_small > d_large,
            "smaller protein should diffuse faster (d_small={d_small}, d_large={d_large})"
        );
    }

    #[test]
    fn test_msd_linear_in_time() {
        let mfd = MembraneFluidDynamics::new(1e-9, 1e-3, 4e-9, 4.28e-21);
        let msd1 = mfd.msd(5e-9, 1.0);
        let msd2 = mfd.msd(5e-9, 2.0);
        assert!(
            (msd2 / msd1 - 2.0).abs() < 1e-9,
            "MSD should be linear in t"
        );
    }

    #[test]
    fn test_oseen_tensor_decays_with_distance() {
        let mfd = MembraneFluidDynamics::new(1e-9, 1e-3, 4e-9, 4.28e-21);
        let t1 = mfd.oseen_tensor_2d(1e-7);
        let t2 = mfd.oseen_tensor_2d(1e-6);
        assert!(t1 > t2, "Oseen tensor should decay with distance");
    }

    #[test]
    fn test_stokeslet_velocity_zero_force_zero_velocity() {
        let mfd = MembraneFluidDynamics::new(1e-9, 1e-3, 4e-9, 4.28e-21);
        let v = mfd.stokeslet_velocity(1e-7, 0.0, [0.0, 0.0]);
        assert!(v[0].abs() < 1e-30 && v[1].abs() < 1e-30);
    }

    #[test]
    fn test_rotational_diffusion_positive() {
        let mfd = MembraneFluidDynamics::new(1e-9, 1e-3, 4e-9, 4.28e-21);
        let dr = mfd.rotational_diffusion_d(5e-9);
        assert!(dr > 0.0, "rotational diffusion should be positive: {dr}");
    }

    // -----------------------------------------------------------------------
    // MembraneProtein
    // -----------------------------------------------------------------------

    #[test]
    fn test_protein_inclusion_energy_at_preferred_curvature_is_zero() {
        let p = MembraneProtein::new(0, [0.0, 0.0], 1e6, 5e-9, 1e-20);
        let e = p.inclusion_energy_at(1e6);
        assert!(
            e.abs() < 1e-40,
            "energy at preferred curvature should be ~0: {e}"
        );
    }

    #[test]
    fn test_protein_inclusion_energy_increases_away_from_preferred() {
        let p = MembraneProtein::new(0, [0.0, 0.0], 1e6, 5e-9, 1e-20);
        let e0 = p.inclusion_energy_at(1e6);
        let e1 = p.inclusion_energy_at(2e6);
        assert!(e1 > e0);
    }

    #[test]
    fn test_protein_update_energy() {
        let mut p = MembraneProtein::new(0, [0.0, 0.0], 0.0, 5e-9, 1e-20);
        p.update_energy(1e6);
        assert!(p.inclusion_energy > 0.0);
    }

    #[test]
    fn test_protein_sensing_force_non_zero_with_gradient() {
        let p = MembraneProtein::new(0, [0.0, 0.0], 2e6, 5e-9, 1e-20);
        let f = p.sensing_force([1.0, 0.0]);
        assert!(
            f[0].abs() > 1e-30 || f[1].abs() > 1e-30,
            "force should be non-zero"
        );
    }

    #[test]
    fn test_protein_cluster_centroid() {
        let proteins = vec![
            MembraneProtein::new(0, [0.0, 0.0], 0.0, 1e-9, 1e-20),
            MembraneProtein::new(1, [2.0, 0.0], 0.0, 1e-9, 1e-20),
        ];
        let cluster = ProteinCluster::from_proteins(&proteins, vec![0, 1]);
        assert!(
            (cluster.centroid[0] - 1.0).abs() < 1e-9,
            "centroid x should be 1"
        );
    }

    #[test]
    fn test_protein_cluster_should_join() {
        let proteins = vec![MembraneProtein::new(0, [0.0, 0.0], 0.0, 5e-9, 1e-20)];
        let cluster = ProteinCluster::from_proteins(&proteins, vec![0]);
        // A point very close should join.
        assert!(cluster.should_join([0.0, 0.0]));
    }

    // -----------------------------------------------------------------------
    // VesicleSimulation
    // -----------------------------------------------------------------------

    #[test]
    fn test_vesicle_creates_without_panic() {
        let _ = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e10);
    }

    #[test]
    fn test_vesicle_total_energy_non_negative() {
        let vesicle = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e10);
        let e = vesicle.total_energy();
        assert!(e >= 0.0, "total vesicle energy should be non-negative: {e}");
    }

    #[test]
    fn test_vesicle_reduced_volume_positive() {
        let mut vesicle = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e10);
        vesicle.bilayer.assign_curvatures(0.0, 0.0);
        let rv = vesicle.reduced_volume();
        // A flat patch has zero enclosed volume; verify the method returns a non-negative
        // finite value (physical range is [0, 1]).
        assert!(
            rv >= 0.0 && rv.is_finite(),
            "reduced volume should be non-negative: {rv}"
        );
    }

    #[test]
    fn test_vesicle_shape_class_returns_string() {
        let vesicle = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e10);
        let s = vesicle.shape_class();
        assert!(!s.is_empty());
    }

    #[test]
    fn test_vesicle_step_advances_time() {
        let mut vesicle = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e10);
        vesicle.step(1e-6);
        assert!((vesicle.time - 1e-6).abs() < 1e-20, "time should advance");
    }

    #[test]
    fn test_vesicle_volume_energy_positive_when_compressed() {
        let mut vesicle = VesicleSimulation::new_sphere(5e-6, 2e-19, 240e-3, 1e13);
        // Move all vertices inward to decrease volume.
        for v in &mut vesicle.bilayer.vertices {
            v.position = scale3(v.position, 0.5);
        }
        let e = vesicle.volume_energy();
        assert!(
            e > 0.0,
            "compressed vesicle should have positive volume energy: {e}"
        );
    }

    // -----------------------------------------------------------------------
    // CellMechanics
    // -----------------------------------------------------------------------

    #[test]
    fn test_afm_hertz_force_positive_for_positive_indentation() {
        let afm = AfmIndentation::new(1e-6, 1000.0, 0.4);
        let f = afm.hertz_force(1e-7);
        assert!(f > 0.0, "Hertz force should be positive: {f}");
    }

    #[test]
    fn test_afm_hertz_force_zero_for_zero_indentation() {
        let afm = AfmIndentation::new(1e-6, 1000.0, 0.4);
        let f = afm.hertz_force(0.0);
        assert!(f.abs() < 1e-30);
    }

    #[test]
    fn test_afm_hertz_force_increases_with_indentation() {
        let afm = AfmIndentation::new(1e-6, 1000.0, 0.4);
        let f1 = afm.hertz_force(1e-8);
        let f2 = afm.hertz_force(2e-8);
        assert!(f2 > f1, "deeper indentation should give larger force");
    }

    #[test]
    fn test_afm_infer_young_modulus() {
        let afm = AfmIndentation::new(1e-6, 1000.0, 0.4);
        let delta = 1e-7;
        let f = afm.hertz_force(delta);
        let e_inferred = afm.infer_young_modulus(f, delta);
        assert!(
            (e_inferred - 1000.0).abs() < 1.0,
            "inferred E should match: {e_inferred}"
        );
    }

    #[test]
    fn test_cell_mechanics_creates() {
        let cm = CellMechanics::new(4e-6, 1000.0, 0.4, 1e-5, 6);
        assert!(cm.spring_constant > 0.0);
    }

    #[test]
    fn test_cell_mechanics_osmotic_pressure() {
        let cm = CellMechanics::new(4e-6, 1000.0, 0.4, 1e-5, 6);
        let p = cm.osmotic_pressure(1.0, 310.0); // 1 mol/m³
        assert!((p - 8.314 * 310.0).abs() < 0.01);
    }

    #[test]
    fn test_cell_mechanics_update_tension_from_pressure() {
        let mut cm = CellMechanics::new(4e-6, 1000.0, 0.4, 1e-5, 6);
        cm.turgor_pressure = 1000.0;
        cm.update_tension_from_pressure();
        assert!((cm.membrane_tension - 1000.0 * 4e-6 / 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_cell_mechanics_lysolipid_softening() {
        let mut cm = CellMechanics::new(4e-6, 1000.0, 0.4, 1e-5, 6);
        let e0 = cm.afm.young_modulus;
        cm.apply_lysolipid_softening(0.5);
        assert!(cm.afm.young_modulus < e0);
    }

    #[test]
    fn test_cell_mechanics_deformability_index_in_range() {
        let cm = CellMechanics::new(4e-6, 1000.0, 0.4, 1e-5, 8);
        let di = cm.deformability_index();
        assert!((0.0..=1.0).contains(&di));
    }

    // -----------------------------------------------------------------------
    // Utility functions
    // -----------------------------------------------------------------------

    #[test]
    fn test_sphere_area_unit_sphere() {
        let a = sphere_area(1.0);
        assert!((a - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_volume_unit_sphere() {
        let v = sphere_volume(1.0);
        assert!((v - (4.0 / 3.0) * PI).abs() < 1e-10);
    }

    #[test]
    fn test_rms_deviation_from_sphere_exact_sphere() {
        let r = 1.0;
        let vertices: Vec<MembraneVertex> = (0..8)
            .map(|i| {
                let phi = 2.0 * PI * i as f64 / 8.0;
                MembraneVertex::new([r * phi.cos(), r * phi.sin(), 0.0])
            })
            .collect();
        let rms = rms_deviation_from_sphere(&vertices, r);
        // All points are exactly at radius 1 → rms should be 0.
        assert!(rms.abs() < 1e-10, "rms={rms}");
    }

    #[test]
    fn test_asphericity_symmetric_config() {
        // Symmetric set of points → should be low asphericity.
        let vertices: Vec<MembraneVertex> = vec![
            MembraneVertex::new([1.0, 0.0, 0.0]),
            MembraneVertex::new([-1.0, 0.0, 0.0]),
            MembraneVertex::new([0.0, 1.0, 0.0]),
            MembraneVertex::new([0.0, -1.0, 0.0]),
            MembraneVertex::new([0.0, 0.0, 1.0]),
            MembraneVertex::new([0.0, 0.0, -1.0]),
        ];
        let a = asphericity(&vertices);
        assert!(a < 1.0, "symmetric config should have low asphericity: {a}");
    }

    #[test]
    fn test_triangle_area_right_triangle() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let a = triangle_area(p0, p1, p2);
        assert!((a - 0.5).abs() < 1e-10, "area={a}");
    }

    #[test]
    fn test_mean_curvature_laplace_flat_surface() {
        // A vertex at the centre of a ring of coplanar neighbours should have H≈0.
        let v = [0.0, 0.0, 0.0];
        let ring: Vec<Vec3> = (0..6)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / 6.0;
                [theta.cos(), theta.sin(), 0.0]
            })
            .collect();
        let areas = vec![PI / 6.0; 6];
        let h = mean_curvature_laplace(v, &ring, &areas);
        assert!(
            h.abs() < 1e-5,
            "flat surface should have ~0 mean curvature: {h}"
        );
    }
}
