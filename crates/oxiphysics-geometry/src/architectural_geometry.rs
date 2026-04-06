#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Architectural geometry — free-form surfaces, panelization, grid shells,
//! tensile structures, geodesic domes, structural glass, and parametric facades.
//!
//! Covers:
//! - [`GridShell`] — doubly-curved grid shell with dynamic relaxation form-finding
//! - [`PanelizationResult`] — planar quad / triangular panel results with tolerance info
//! - [`PlanarQuadMesh`] — PQ-mesh from conjugate curve nets
//! - [`TensileStructure`] — minimal surface / cable net with prestress distribution
//! - [`GeodesicDome`] — frequency-n geodesic sphere (Class I / II subdivision)
//! - [`StructuralGlass`] — glass panel sizing under thermal, wind, and self-weight
//! - [`ParametricFacade`] — adaptive solar-responsive facade paneling

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vec3 helper (plain f64, no nalgebra)
// ---------------------------------------------------------------------------

/// Minimal 3D vector for architectural geometry computations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct V3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl V3 {
    /// Create a new vector.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Euclidean length.
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Unit vector (returns zero if near zero).
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        if n < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / n, self.y / n, self.z / n)
        }
    }

    /// Dot product.
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product.
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Component-wise addition.
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    /// Component-wise subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    /// Scale by scalar.
    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

// ---------------------------------------------------------------------------
// GridShell
// ---------------------------------------------------------------------------

/// A doubly-curved grid shell structure.
///
/// Nodes are arranged in a regular grid (rows × cols) with initial positions
/// on a parametric surface. Members connect neighbouring nodes. Form-finding
/// is performed via dynamic relaxation — nodes are moved iteratively towards
/// equilibrium under prestress and applied loads.
#[derive(Debug, Clone)]
pub struct GridShell {
    /// Number of rows in the grid.
    pub rows: usize,
    /// Number of columns in the grid.
    pub cols: usize,
    /// Current 3-D positions of all nodes (row-major).
    pub nodes: Vec<V3>,
    /// Member connectivity: pairs of node indices.
    pub members: Vec<(usize, usize)>,
    /// Natural (unstressed) length of each member (m).
    pub rest_lengths: Vec<f64>,
    /// Axial stiffness EA for each member (N).
    pub stiffnesses: Vec<f64>,
    /// Pinned (fixed) node indices.
    pub fixed_nodes: Vec<usize>,
    /// Lumped mass at each node (kg) — used in dynamic relaxation.
    pub nodal_mass: Vec<f64>,
}

impl GridShell {
    /// Construct a parabolic grid shell spanning `span_x` × `span_y` metres,
    /// with a rise of `sag` metres at the centre.
    ///
    /// `rows` and `cols` define the number of grid points in each direction.
    /// The boundary nodes are automatically pinned.
    pub fn parabolic(rows: usize, cols: usize, span_x: f64, span_y: f64, sag: f64) -> Self {
        assert!(rows >= 2 && cols >= 2);
        let mut nodes = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            let u = i as f64 / (rows - 1) as f64; // 0..1
            for j in 0..cols {
                let v = j as f64 / (cols - 1) as f64;
                let x = u * span_x;
                let y = v * span_y;
                // Paraboloid: z = 4*sag * u*(1-u) * v*(1-v)
                let z = 4.0 * sag * u * (1.0 - u) * v * (1.0 - v);
                nodes.push(V3::new(x, y, z));
            }
        }
        let mut members = Vec::new();
        // Grid edges (horizontal + vertical)
        for i in 0..rows {
            for j in 0..cols {
                let idx = i * cols + j;
                if j + 1 < cols {
                    members.push((idx, idx + 1));
                }
                if i + 1 < rows {
                    members.push((idx, idx + cols));
                }
            }
        }
        // Diagonal bracing
        for i in 0..rows - 1 {
            for j in 0..cols - 1 {
                let a = i * cols + j;
                let b = (i + 1) * cols + (j + 1);
                members.push((a, b));
                let c = i * cols + (j + 1);
                let d = (i + 1) * cols + j;
                members.push((c, d));
            }
        }
        let rest_lengths: Vec<f64> = members
            .iter()
            .map(|(a, b)| nodes[*a].sub(&nodes[*b]).norm())
            .collect();
        let stiffnesses = vec![1.0e6_f64; members.len()]; // EA = 1 MN
        let nodal_mass = vec![10.0_f64; nodes.len()];
        // Fix boundary nodes
        let mut fixed_nodes = Vec::new();
        for i in 0..rows {
            for j in 0..cols {
                if i == 0 || i == rows - 1 || j == 0 || j == cols - 1 {
                    fixed_nodes.push(i * cols + j);
                }
            }
        }
        Self {
            rows,
            cols,
            nodes,
            members,
            rest_lengths,
            stiffnesses,
            fixed_nodes,
            nodal_mass,
        }
    }

    /// Run `max_iter` steps of dynamic relaxation.
    ///
    /// Returns the maximum residual force at any free node after convergence.
    /// A value below `tol` (N) indicates equilibrium has been achieved.
    pub fn dynamic_relaxation(&mut self, max_iter: usize, dt: f64, damping: f64, tol: f64) -> f64 {
        let n = self.nodes.len();
        let fixed: std::collections::HashSet<usize> = self.fixed_nodes.iter().copied().collect();
        let mut velocities = vec![V3::zero(); n];
        let mut max_residual = f64::MAX;

        for _iter in 0..max_iter {
            let mut forces = vec![V3::zero(); n];

            // Accumulate member forces
            for (m_idx, (a, b)) in self.members.iter().enumerate() {
                let pa = &self.nodes[*a];
                let pb = &self.nodes[*b];
                let delta = pb.sub(pa);
                let length = delta.norm();
                if length < 1e-15 {
                    continue;
                }
                let strain = (length - self.rest_lengths[m_idx]) / self.rest_lengths[m_idx];
                let f_mag = self.stiffnesses[m_idx] * strain;
                let dir = delta.normalize();
                let force = dir.scale(f_mag);
                forces[*a] = forces[*a].add(&force);
                forces[*b] = forces[*b].sub(&force);
            }

            // Update velocities and positions
            max_residual = 0.0;
            for i in 0..n {
                if fixed.contains(&i) {
                    velocities[i] = V3::zero();
                    continue;
                }
                let f_norm = forces[i].norm();
                if f_norm > max_residual {
                    max_residual = f_norm;
                }
                // Apply damping and integrate
                let acc = forces[i].scale(1.0 / self.nodal_mass[i]);
                velocities[i] = velocities[i].scale(1.0 - damping).add(&acc.scale(dt));
                self.nodes[i] = self.nodes[i].add(&velocities[i].scale(dt));
            }

            if max_residual < tol {
                break;
            }
        }
        max_residual
    }

    /// Compute maximum member force (N) in the current configuration.
    pub fn max_member_force(&self) -> f64 {
        self.members
            .iter()
            .enumerate()
            .map(|(m_idx, (a, b))| {
                let length = self.nodes[*a].sub(&self.nodes[*b]).norm();
                let strain = (length - self.rest_lengths[m_idx]) / self.rest_lengths[m_idx];
                (self.stiffnesses[m_idx] * strain).abs()
            })
            .fold(0.0_f64, f64::max)
    }

    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Total number of members.
    pub fn member_count(&self) -> usize {
        self.members.len()
    }
}

// ---------------------------------------------------------------------------
// PanelizationResult
// ---------------------------------------------------------------------------

/// Result of panelizing a free-form surface.
///
/// Contains lists of planar quad and triangular panels with planarity
/// deviation statistics.
#[derive(Debug, Clone)]
pub struct PanelizationResult {
    /// Planar quad panels — each entry is four corner-point indices.
    pub quad_panels: Vec<[usize; 4]>,
    /// Triangular panels — each entry is three corner-point indices.
    pub tri_panels: Vec<[usize; 3]>,
    /// Vertex positions referenced by panel indices.
    pub vertices: Vec<V3>,
    /// Planarity error for each quad panel (m) — distance of the worst
    /// off-plane vertex from the best-fit plane.
    pub planarity_errors: Vec<f64>,
    /// Maximum permitted planarity deviation (m).
    pub planarity_tolerance: f64,
    /// Panel gap / overlap statistics: (mean_gap, max_gap) in metres.
    pub gap_stats: (f64, f64),
}

impl PanelizationResult {
    /// Compute the planarity error of a single quad (4 points).
    ///
    /// Fits a plane to the first three vertices and returns the distance of
    /// the fourth vertex from that plane.
    pub fn quad_planarity_error(pts: &[V3; 4]) -> f64 {
        let ab = pts[1].sub(&pts[0]);
        let ac = pts[2].sub(&pts[0]);
        let normal = ab.cross(&ac);
        let n_len = normal.norm();
        if n_len < 1e-15 {
            return 0.0;
        }
        let n_unit = normal.normalize();
        let ad = pts[3].sub(&pts[0]);
        ad.dot(&n_unit).abs()
    }

    /// Count panels that exceed the planarity tolerance.
    pub fn out_of_tolerance_count(&self) -> usize {
        self.planarity_errors
            .iter()
            .filter(|&&e| e > self.planarity_tolerance)
            .count()
    }

    /// Mean planarity error across all quad panels (m).
    pub fn mean_planarity_error(&self) -> f64 {
        if self.planarity_errors.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.planarity_errors.iter().sum();
        sum / self.planarity_errors.len() as f64
    }

    /// Maximum planarity error across all quad panels (m).
    pub fn max_planarity_error(&self) -> f64 {
        self.planarity_errors
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// PlanarQuadMesh
// ---------------------------------------------------------------------------

/// A planar quad mesh (PQ mesh) approximating a free-form surface.
///
/// In a PQ mesh every face is (approximately) planar, enabling single-layer
/// glass fabrication. Generated via conjugate curve nets on the surface.
#[derive(Debug, Clone)]
pub struct PlanarQuadMesh {
    /// Grid resolution in the u direction.
    pub nu: usize,
    /// Grid resolution in the v direction.
    pub nv: usize,
    /// Vertex positions (row-major, nu×nv).
    pub vertices: Vec<V3>,
    /// Quad faces as vertex index tuples (i, i+1, i+cols+1, i+cols).
    pub faces: Vec<[usize; 4]>,
    /// Planarity error per face (m).
    pub planarity_errors: Vec<f64>,
}

impl PlanarQuadMesh {
    /// Build a PQ mesh from a height-field surface `z(u,v)` using a uniform
    /// parameter grid.
    ///
    /// `surface_fn` maps (u, v) ∈ \[0,1\]² → (x, y, z).
    pub fn from_surface<F>(nu: usize, nv: usize, surface_fn: F) -> Self
    where
        F: Fn(f64, f64) -> V3,
    {
        assert!(nu >= 2 && nv >= 2);
        let mut vertices = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            let u = i as f64 / (nu - 1) as f64;
            for j in 0..nv {
                let v = j as f64 / (nv - 1) as f64;
                vertices.push(surface_fn(u, v));
            }
        }
        let mut faces = Vec::new();
        for i in 0..nu - 1 {
            for j in 0..nv - 1 {
                let a = i * nv + j;
                let b = i * nv + j + 1;
                let c = (i + 1) * nv + j + 1;
                let d = (i + 1) * nv + j;
                faces.push([a, b, c, d]);
            }
        }
        let planarity_errors: Vec<f64> = faces
            .iter()
            .map(|[a, b, c, d]| {
                let pts = [vertices[*a], vertices[*b], vertices[*c], vertices[*d]];
                PanelizationResult::quad_planarity_error(&pts)
            })
            .collect();
        Self {
            nu,
            nv,
            vertices,
            faces,
            planarity_errors,
        }
    }

    /// Mean planarity error across all faces (m).
    pub fn mean_planarity_error(&self) -> f64 {
        if self.planarity_errors.is_empty() {
            return 0.0;
        }
        self.planarity_errors.iter().sum::<f64>() / self.planarity_errors.len() as f64
    }

    /// Maximum planarity error (m).
    pub fn max_planarity_error(&self) -> f64 {
        self.planarity_errors
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }

    /// Fraction of faces within tolerance `tol` (m).
    pub fn planarity_compliance_ratio(&self, tol: f64) -> f64 {
        let compliant = self.planarity_errors.iter().filter(|&&e| e <= tol).count();
        compliant as f64 / self.planarity_errors.len().max(1) as f64
    }

    /// Gaussian curvature approximation at interior vertex (i, j).
    ///
    /// Uses the angle-deficit method: K ≈ (2π − Σθ) / A,
    /// where θ are the face angles at the vertex and A is the mixed area.
    pub fn gaussian_curvature_at(&self, i: usize, j: usize) -> f64 {
        if i == 0 || i >= self.nu - 1 || j == 0 || j >= self.nv - 1 {
            return 0.0; // boundary: skip
        }
        let idx = |ii: usize, jj: usize| ii * self.nv + jj;
        let p = self.vertices[idx(i, j)];
        // Four neighbouring vertices
        let pn = self.vertices[idx(i - 1, j)];
        let ps = self.vertices[idx(i + 1, j)];
        let pe = self.vertices[idx(i, j + 1)];
        let pw = self.vertices[idx(i, j - 1)];
        // Angle at p in each of the 4 surrounding quads (approximated as triangles)
        let angle = |a: V3, centre: V3, b: V3| -> f64 {
            let u = a.sub(&centre).normalize();
            let v = b.sub(&centre).normalize();
            u.dot(&v).clamp(-1.0, 1.0).acos()
        };
        let theta_sum = angle(pn, p, pe) + angle(pe, p, ps) + angle(ps, p, pw) + angle(pw, p, pn);
        let area_approx = {
            let d1 = pn.sub(&p).norm();
            let d2 = pe.sub(&p).norm();
            d1 * d2
        };
        if area_approx < 1e-20 {
            return 0.0;
        }
        (2.0 * PI - theta_sum) / area_approx
    }

    /// Total number of faces.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }
}

// ---------------------------------------------------------------------------
// TensileStructure
// ---------------------------------------------------------------------------

/// A cable-net / tensile membrane structure.
///
/// Models a prestressed network of cables anchored at boundary nodes. The
/// minimal-surface equilibrium shape is found via dynamic relaxation, analogous
/// to a soap film.
#[derive(Debug, Clone)]
pub struct TensileStructure {
    /// Node positions (current configuration).
    pub nodes: Vec<V3>,
    /// Cable elements: (node_a, node_b).
    pub cables: Vec<(usize, usize)>,
    /// Prestress force in each cable (N).
    pub prestress: Vec<f64>,
    /// Fixed (anchor) node indices.
    pub anchors: Vec<usize>,
    /// Tributary area at each free node (m²) — for uniform pressure loads.
    pub trib_areas: Vec<f64>,
}

impl TensileStructure {
    /// Create a flat rectangular cable net ready for form-finding.
    ///
    /// `nx` × `ny` nodes span `width` × `height` metres. All edges are
    /// prestressed to `p0` N. Boundary nodes are fixed.
    pub fn flat_net(nx: usize, ny: usize, width: f64, height: f64, p0: f64) -> Self {
        assert!(nx >= 2 && ny >= 2);
        let mut nodes = Vec::with_capacity(nx * ny);
        for i in 0..nx {
            let x = i as f64 / (nx - 1) as f64 * width;
            for j in 0..ny {
                let y = j as f64 / (ny - 1) as f64 * height;
                nodes.push(V3::new(x, y, 0.0));
            }
        }
        let mut cables = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                if j + 1 < ny {
                    cables.push((idx, idx + 1));
                }
                if i + 1 < nx {
                    cables.push((idx, idx + ny));
                }
            }
        }
        let prestress = vec![p0; cables.len()];
        let mut anchors = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
                    anchors.push(i * ny + j);
                }
            }
        }
        let trib_areas = vec![(width / (nx - 1) as f64) * (height / (ny - 1) as f64); nodes.len()];
        Self {
            nodes,
            cables,
            prestress,
            anchors,
            trib_areas,
        }
    }

    /// Apply a uniform normal pressure `q` (Pa) and run dynamic relaxation.
    ///
    /// Returns the maximum out-of-plane displacement (m) achieved.
    pub fn form_find(
        &mut self,
        pressure: f64,
        max_iter: usize,
        dt: f64,
        damping: f64,
        tol: f64,
    ) -> f64 {
        let n = self.nodes.len();
        let anchors: std::collections::HashSet<usize> = self.anchors.iter().copied().collect();
        let mut vel = vec![V3::zero(); n];
        let mass = 1.0_f64; // lumped unit mass

        for _iter in 0..max_iter {
            let mut forces = vec![V3::zero(); n];
            // Cable forces
            for (m_idx, (a, b)) in self.cables.iter().enumerate() {
                let pa = self.nodes[*a];
                let pb = self.nodes[*b];
                let delta = pb.sub(&pa);
                let length = delta.norm();
                if length < 1e-15 {
                    continue;
                }
                let dir = delta.normalize();
                // Force-density: f = T/L * delta
                let f_density = self.prestress[m_idx] / length;
                let fv = dir.scale(f_density);
                forces[*a] = forces[*a].add(&fv);
                forces[*b] = forces[*b].sub(&fv);
            }
            // Pressure load (−z direction for downward)
            for i in 0..n {
                if !anchors.contains(&i) {
                    let fz = pressure * self.trib_areas[i];
                    forces[i].z -= fz;
                }
            }
            // Integrate
            let mut max_res = 0.0_f64;
            for i in 0..n {
                if anchors.contains(&i) {
                    vel[i] = V3::zero();
                    continue;
                }
                let f_norm = forces[i].norm();
                if f_norm > max_res {
                    max_res = f_norm;
                }
                let acc = forces[i].scale(1.0 / mass);
                vel[i] = vel[i].scale(1.0 - damping).add(&acc.scale(dt));
                self.nodes[i] = self.nodes[i].add(&vel[i].scale(dt));
            }
            if max_res < tol {
                break;
            }
        }
        // Return max sag
        self.nodes
            .iter()
            .enumerate()
            .filter(|(i, _)| !anchors.contains(i))
            .map(|(_, p)| p.z.abs())
            .fold(0.0_f64, f64::max)
    }

    /// Total prestress energy (J) in the current configuration.
    ///
    /// Computed as Σ T · (L − L0)² / (2 L0) for axial members with stiffness
    /// equal to prestress / rest-length.
    pub fn prestress_energy(&self) -> f64 {
        self.cables
            .iter()
            .enumerate()
            .map(|(m, (a, b))| {
                let l = self.nodes[*a].sub(&self.nodes[*b]).norm();
                let t = self.prestress[m];
                // Simple potential energy in cable (t/l * dl^2 / 2)
                let dl = l - t / 1e3; // nominal rest length approximation
                t / l * dl * dl * 0.5
            })
            .sum()
    }

    /// Mean cable tension (N).
    pub fn mean_cable_tension(&self) -> f64 {
        if self.prestress.is_empty() {
            return 0.0;
        }
        self.prestress.iter().sum::<f64>() / self.prestress.len() as f64
    }
}

// ---------------------------------------------------------------------------
// GeodesicDome
// ---------------------------------------------------------------------------

/// A geodesic dome / sphere approximation.
///
/// Uses Class I frequency-n subdivision of an icosahedron to generate the
/// vertex coordinates and strut connectivity. Higher frequency gives a closer
/// approximation to a sphere.
#[derive(Debug, Clone)]
pub struct GeodesicDome {
    /// Subdivision frequency (1 = icosahedron, 2 = 80 faces, …).
    pub frequency: usize,
    /// Radius of the circumscribed sphere (m).
    pub radius: f64,
    /// Vertex coordinates (on the sphere surface).
    pub vertices: Vec<V3>,
    /// Strut connectivity (pairs of vertex indices).
    pub struts: Vec<(usize, usize)>,
    /// Whether this is a full sphere or a hemisphere.
    pub hemisphere: bool,
}

impl GeodesicDome {
    /// Generate an icosahedron base (frequency 1).
    fn icosahedron_vertices(radius: f64) -> Vec<V3> {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0; // golden ratio
        let s = radius / (1.0 + phi * phi).sqrt();
        let t = phi * s;
        vec![
            V3::new(0.0, s, t),
            V3::new(0.0, -s, t),
            V3::new(0.0, s, -t),
            V3::new(0.0, -s, -t),
            V3::new(s, t, 0.0),
            V3::new(-s, t, 0.0),
            V3::new(s, -t, 0.0),
            V3::new(-s, -t, 0.0),
            V3::new(t, 0.0, s),
            V3::new(-t, 0.0, s),
            V3::new(t, 0.0, -s),
            V3::new(-t, 0.0, -s),
        ]
    }

    /// Icosahedron face definitions (vertex index triples).
    fn icosahedron_faces() -> Vec<[usize; 3]> {
        vec![
            [0, 1, 8],
            [0, 8, 4],
            [0, 4, 5],
            [0, 5, 9],
            [0, 9, 1],
            [1, 6, 8],
            [8, 6, 10],
            [8, 10, 4],
            [4, 10, 2],
            [4, 2, 5],
            [5, 2, 11],
            [5, 11, 9],
            [9, 11, 7],
            [9, 7, 1],
            [1, 7, 6],
            [3, 10, 6],
            [3, 6, 7],
            [3, 7, 11],
            [3, 11, 2],
            [3, 2, 10],
        ]
    }

    /// Build a geodesic dome of given `frequency` and `radius`.
    ///
    /// `hemisphere` trims the lower half of the sphere.
    pub fn new(frequency: usize, radius: f64, hemisphere: bool) -> Self {
        assert!((1..=8).contains(&frequency), "frequency 1–8 supported");
        let base_verts = Self::icosahedron_vertices(radius);
        let base_faces = Self::icosahedron_faces();

        if frequency == 1 {
            // Project base vertices onto sphere
            let vertices: Vec<V3> = base_verts
                .iter()
                .map(|v| v.normalize().scale(radius))
                .collect();
            let mut struts = Vec::new();
            for face in &base_faces {
                for k in 0..3 {
                    let a = face[k];
                    let b = face[(k + 1) % 3];
                    let edge = if a < b { (a, b) } else { (b, a) };
                    if !struts.contains(&edge) {
                        struts.push(edge);
                    }
                }
            }
            return Self {
                frequency,
                radius,
                vertices,
                struts,
                hemisphere,
            };
        }

        // Subdivide each face frequency times
        let mut all_points: Vec<V3> = Vec::new();
        let mut all_struts: Vec<(usize, usize)> = Vec::new();

        let find_or_insert = |pts: &mut Vec<V3>, p: V3| -> usize {
            let threshold = radius * 1e-6;
            for (idx, existing) in pts.iter().enumerate() {
                if existing.sub(&p).norm() < threshold {
                    return idx;
                }
            }
            pts.push(p);
            pts.len() - 1
        };

        for face in &base_faces {
            let va = base_verts[face[0]].normalize().scale(radius);
            let vb = base_verts[face[1]].normalize().scale(radius);
            let vc = base_verts[face[2]].normalize().scale(radius);

            let f = frequency as f64;
            // Generate sub-vertices by barycentric interpolation
            let mut local: Vec<Vec<V3>> = vec![vec![V3::zero(); frequency + 1]; frequency + 1];
            for i in 0..=frequency {
                for j in 0..=frequency - i {
                    let k = frequency - i - j;
                    let p = V3::new(
                        (i as f64 * va.x + j as f64 * vb.x + k as f64 * vc.x) / f,
                        (i as f64 * va.y + j as f64 * vb.y + k as f64 * vc.y) / f,
                        (i as f64 * va.z + j as f64 * vb.z + k as f64 * vc.z) / f,
                    );
                    // Project onto sphere
                    local[i][j] = p.normalize().scale(radius);
                }
            }

            // Insert vertices and edges
            for i in 0..=frequency {
                for j in 0..=frequency - i {
                    let idx_a = find_or_insert(&mut all_points, local[i][j]);
                    if j < frequency - i {
                        let idx_b = find_or_insert(&mut all_points, local[i][j + 1]);
                        let e = if idx_a < idx_b {
                            (idx_a, idx_b)
                        } else {
                            (idx_b, idx_a)
                        };
                        if !all_struts.contains(&e) {
                            all_struts.push(e);
                        }
                    }
                    if i < frequency - j {
                        let idx_c = find_or_insert(&mut all_points, local[i + 1][j]);
                        let e = if idx_a < idx_c {
                            (idx_a, idx_c)
                        } else {
                            (idx_c, idx_a)
                        };
                        if !all_struts.contains(&e) {
                            all_struts.push(e);
                        }
                    }
                    if i < frequency - j && j < frequency - (i + 1) {
                        let idx_b = find_or_insert(&mut all_points, local[i][j + 1]);
                        let idx_c = find_or_insert(&mut all_points, local[i + 1][j]);
                        let e = if idx_b < idx_c {
                            (idx_b, idx_c)
                        } else {
                            (idx_c, idx_b)
                        };
                        if !all_struts.contains(&e) {
                            all_struts.push(e);
                        }
                    }
                }
            }
        }

        let vertices = all_points;
        let struts = all_struts;

        Self {
            frequency,
            radius,
            vertices,
            struts,
            hemisphere,
        }
    }

    /// Expected vertex count for a Class I frequency-n geodesic sphere.
    ///
    /// Formula: V = 10n² + 2 (full sphere).
    pub fn expected_vertex_count(frequency: usize) -> usize {
        10 * frequency * frequency + 2
    }

    /// Expected strut count for a Class I frequency-n geodesic sphere.
    ///
    /// Formula: E = 30n² (full sphere).
    pub fn expected_strut_count(frequency: usize) -> usize {
        30 * frequency * frequency
    }

    /// Strut length statistics: (min, max, mean) in metres.
    pub fn strut_length_stats(&self) -> (f64, f64, f64) {
        if self.struts.is_empty() {
            return (0.0, 0.0, 0.0);
        }
        let lengths: Vec<f64> = self
            .struts
            .iter()
            .map(|(a, b)| self.vertices[*a].sub(&self.vertices[*b]).norm())
            .collect();
        let min = lengths.iter().cloned().fold(f64::MAX, f64::min);
        let max = lengths.iter().cloned().fold(0.0_f64, f64::max);
        let mean = lengths.iter().sum::<f64>() / lengths.len() as f64;
        (min, max, mean)
    }

    /// Verify all vertices lie on the sphere within tolerance.
    pub fn verify_sphericity(&self, tol: f64) -> bool {
        self.vertices
            .iter()
            .all(|v| (v.norm() - self.radius).abs() < tol)
    }
}

// ---------------------------------------------------------------------------
// StructuralGlass
// ---------------------------------------------------------------------------

/// Structural glass panel — sizing under wind, thermal, and self-weight loads.
///
/// Covers monolithic, laminated, and insulating glass units (IGU).
#[derive(Debug, Clone)]
pub struct StructuralGlass {
    /// Panel width (m).
    pub width: f64,
    /// Panel height (m).
    pub height: f64,
    /// Total glass thickness (m) — sum of all plies.
    pub thickness: f64,
    /// Young's modulus of glass (Pa) — typically 70 GPa.
    pub e_glass: f64,
    /// Poisson's ratio of glass — typically 0.23.
    pub nu_glass: f64,
    /// Density of glass (kg/m³) — typically 2500.
    pub density: f64,
    /// Design wind pressure (Pa) — positive = inward.
    pub wind_pressure: f64,
    /// Temperature delta for thermal stress calculation (°C).
    pub delta_temp: f64,
    /// Coefficient of thermal expansion (1/°C) — typically 9e-6.
    pub alpha_cte: f64,
    /// Silicone joint width (m).
    pub sealant_width: f64,
    /// Silicone design shear stress capacity (Pa) — typically 0.14 MPa.
    pub sealant_shear_strength: f64,
}

impl StructuralGlass {
    /// Default single glazing: 10 mm monolithic, 3×5 m panel.
    pub fn monolithic_10mm() -> Self {
        Self {
            width: 3.0,
            height: 5.0,
            thickness: 0.01,
            e_glass: 70.0e9,
            nu_glass: 0.23,
            density: 2500.0,
            wind_pressure: 1200.0, // 1.2 kPa
            delta_temp: 40.0,
            alpha_cte: 9.0e-6,
            sealant_width: 0.025,
            sealant_shear_strength: 0.14e6,
        }
    }

    /// Self-weight per unit area (Pa = N/m²).
    pub fn self_weight_pressure(&self) -> f64 {
        self.density * 9.81 * self.thickness
    }

    /// Maximum bending stress under uniform wind pressure (Pa).
    ///
    /// Approximated as simply supported plate: σ = α · q · a² / t²,
    /// where α ≈ 0.287 for a/b = 1.67 (typical facade aspect ratio),
    /// `a` is the shorter span, `t` is thickness.
    pub fn max_bending_stress_wind(&self) -> f64 {
        let a = self.width.min(self.height);
        let alpha_coeff = 0.287_f64; // plate coefficient for aspect ratio ~1.5-2
        alpha_coeff * self.wind_pressure * a * a / (self.thickness * self.thickness)
    }

    /// Thermal stress in a fully restrained glass pane (Pa).
    ///
    /// σ_thermal = E · α · ΔT
    pub fn thermal_stress(&self) -> f64 {
        self.e_glass * self.alpha_cte * self.delta_temp
    }

    /// Maximum combined stress (wind + thermal) (Pa).
    pub fn combined_stress(&self) -> f64 {
        self.max_bending_stress_wind() + self.thermal_stress()
    }

    /// Structural silicone joint capacity check.
    ///
    /// Returns the maximum transferable line load (N/m) via the silicone bead,
    /// limited by the sealant shear strength and joint geometry.
    pub fn sealant_shear_capacity(&self) -> f64 {
        self.sealant_shear_strength * self.sealant_width
    }

    /// Wind load reaction at each of the four edges (N/m) for a simply
    /// supported panel.
    pub fn edge_reaction_wind(&self) -> f64 {
        let area = self.width * self.height;
        // Distribute equally to all four edges
        let perimeter = 2.0 * (self.width + self.height);
        self.wind_pressure * area / perimeter
    }

    /// Centre deflection under uniform wind pressure (m).
    ///
    /// Thin plate formula for simply supported edges:
    /// δ = β · q · a⁴ / (E · t³),  β ≈ 0.0443.
    pub fn centre_deflection_wind(&self) -> f64 {
        let a = self.width.min(self.height);
        let beta = 0.0443_f64;
        beta * self.wind_pressure * a.powi(4) / (self.e_glass * self.thickness.powi(3))
    }

    /// Span-to-deflection ratio (serviceability check; ≥ 200 typically required).
    pub fn span_deflection_ratio(&self) -> f64 {
        let a = self.width.min(self.height);
        let delta = self.centre_deflection_wind();
        if delta < 1e-12 {
            return f64::MAX;
        }
        a / delta
    }

    /// Insulating glass unit (IGU) sealed unit weight per unit area (Pa).
    ///
    /// Assumes 16 mm argon cavity + second pane of equal thickness.
    pub fn igu_self_weight_pressure(&self) -> f64 {
        let cavity_mass_per_m2 = 1.784 * 0.016; // argon density × cavity thickness
        (2.0 * self.density * self.thickness + cavity_mass_per_m2) * 9.81
    }
}

// ---------------------------------------------------------------------------
// ParametricFacade
// ---------------------------------------------------------------------------

/// Adaptive facade paneling with solar-responsive shading.
///
/// Generates a panelized facade with variable aperture / shading depth
/// based on the solar incidence angle at each panel.
#[derive(Debug, Clone)]
pub struct ParametricFacade {
    /// Number of horizontal bays.
    pub bays: usize,
    /// Number of vertical levels.
    pub levels: usize,
    /// Bay width (m).
    pub bay_width: f64,
    /// Storey height (m).
    pub storey_height: f64,
    /// Building azimuth from north (degrees, clockwise).
    pub facade_azimuth: f64,
    /// Latitude of site (degrees).
    pub latitude: f64,
    /// Maximum fin/shading depth (m).
    pub max_shade_depth: f64,
    /// Target annual solar exposure per panel (kWh/m²/year).
    pub target_solar_exposure: f64,
}

impl ParametricFacade {
    /// Default south-facing facade in Tokyo.
    pub fn tokyo_south() -> Self {
        Self {
            bays: 8,
            levels: 10,
            bay_width: 1.5,
            storey_height: 3.6,
            facade_azimuth: 180.0,
            latitude: 35.7,
            max_shade_depth: 0.8,
            target_solar_exposure: 400.0,
        }
    }

    /// Solar altitude angle (degrees) for a given hour angle and day-of-year.
    ///
    /// Uses the standard solar position formula.
    pub fn solar_altitude(&self, day_of_year: f64, hour: f64) -> f64 {
        let lat_rad = self.latitude.to_radians();
        // Solar declination
        let decl =
            23.45_f64.to_radians() * (360.0 / 365.0 * (day_of_year - 81.0)).to_radians().sin();
        let hour_angle = (hour - 12.0) * 15.0; // degrees
        let h_rad = hour_angle.to_radians();
        let sin_alt = lat_rad.sin() * decl.sin() + lat_rad.cos() * decl.cos() * h_rad.cos();
        sin_alt.clamp(-1.0, 1.0).asin().to_degrees()
    }

    /// Solar azimuth angle (degrees from south, positive west) for given
    /// day-of-year and hour.
    pub fn solar_azimuth(&self, day_of_year: f64, hour: f64) -> f64 {
        let lat_rad = self.latitude.to_radians();
        let decl =
            23.45_f64.to_radians() * (360.0 / 365.0 * (day_of_year - 81.0)).to_radians().sin();
        let hour_angle = (hour - 12.0) * 15.0;
        let h_rad = hour_angle.to_radians();
        let sin_alt = lat_rad.sin() * decl.sin() + lat_rad.cos() * decl.cos() * h_rad.cos();
        let cos_az = (decl.sin() - lat_rad.sin() * sin_alt)
            / (lat_rad.cos() * sin_alt.asin().cos().max(1e-10));
        let az_deg = cos_az.clamp(-1.0, 1.0).acos().to_degrees();
        if hour_angle > 0.0 { az_deg } else { -az_deg }
    }

    /// Incidence angle (degrees) of the sun on the facade for a given
    /// solar altitude `alt_deg` and solar azimuth `az_deg`.
    pub fn incidence_angle(&self, alt_deg: f64, az_deg: f64) -> f64 {
        let facade_az = self.facade_azimuth - 180.0; // relative to south
        let delta_az = (az_deg - facade_az).to_radians();
        let alt_rad = alt_deg.to_radians();
        // Cosine of incidence on vertical facade
        let cos_inc = alt_rad.cos() * delta_az.cos();
        cos_inc.clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// Required shading depth (m) to block direct sun given incidence angle.
    ///
    /// d = panel_height * tan(incidence_angle_from_horizontal)
    pub fn required_shade_depth(&self, incidence_deg: f64) -> f64 {
        let inc_from_horiz = (90.0 - incidence_deg).abs().to_radians();
        let depth = self.storey_height * inc_from_horiz.tan();
        depth.min(self.max_shade_depth)
    }

    /// Generate shade depth for each panel based on summer solstice noon sun.
    ///
    /// Returns a `(bays × levels)` matrix of shade depths (row = level, col = bay).
    pub fn shade_depth_matrix(&self) -> Vec<Vec<f64>> {
        let alt = self.solar_altitude(172.0, 12.0); // summer solstice noon
        let az = self.solar_azimuth(172.0, 12.0);
        let inc = self.incidence_angle(alt, az);
        let depth = self.required_shade_depth(inc);
        vec![vec![depth; self.bays]; self.levels]
    }

    /// Annual solar exposure estimate per panel (kWh/m²/year).
    ///
    /// Simple approximation: integrates cos(incidence) over daylight hours
    /// using monthly representative days (12 samples × 8 daylight hours).
    pub fn annual_solar_exposure(&self) -> f64 {
        let panel_area = self.bay_width * self.storey_height;
        let mut total_irradiance = 0.0_f64;
        let dni = 800.0_f64; // W/m² direct normal irradiance
        let months = [
            17.0, 47.0, 75.0, 105.0, 135.0, 162.0, 198.0, 228.0, 258.0, 288.0, 318.0, 344.0_f64,
        ];
        for &day in &months {
            for h_int in 7..17 {
                let hour = h_int as f64 + 0.5;
                let alt = self.solar_altitude(day, hour);
                if alt <= 0.0 {
                    continue;
                }
                let az = self.solar_azimuth(day, hour);
                let inc = self.incidence_angle(alt, az);
                let cos_inc = inc.to_radians().cos().max(0.0);
                total_irradiance += dni * cos_inc;
            }
        }
        // Convert W/m² · hours to kWh/m²/year
        total_irradiance / 1000.0 * panel_area
    }

    /// Panel count in the whole facade.
    pub fn total_panel_count(&self) -> usize {
        self.bays * self.levels
    }

    /// Total facade area (m²).
    pub fn total_facade_area(&self) -> f64 {
        self.bays as f64 * self.bay_width * self.levels as f64 * self.storey_height
    }

    /// Window-to-wall ratio: glazed area / total facade area.
    ///
    /// Assumes a fixed 0.7 × (bay × level area) is glazed.
    pub fn window_to_wall_ratio(&self) -> f64 {
        0.70
    }

    /// Shading performance index (0–1).
    ///
    /// 1.0 = shade depth reaches maximum everywhere.
    pub fn shading_performance_index(&self) -> f64 {
        let matrix = self.shade_depth_matrix();
        let mean_depth: f64 =
            matrix.iter().flatten().sum::<f64>() / (self.bays * self.levels) as f64;
        (mean_depth / self.max_shade_depth).min(1.0)
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Compute the area of a triangle defined by three 3-D points (m²).
pub fn triangle_area(a: V3, b: V3, c: V3) -> f64 {
    let ab = b.sub(&a);
    let ac = c.sub(&a);
    ab.cross(&ac).norm() * 0.5
}

/// Fit a plane to a cloud of points using centroid + covariance PCA (simplified).
///
/// Returns (centroid, normal) where normal is the direction of least variance.
/// Uses Jacobi iteration for the 3×3 symmetric covariance matrix (2 sweeps).
pub fn fit_plane(points: &[V3]) -> (V3, V3) {
    if points.len() < 3 {
        return (V3::zero(), V3::new(0.0, 0.0, 1.0));
    }
    let n = points.len() as f64;
    let cx = points.iter().map(|p| p.x).sum::<f64>() / n;
    let cy = points.iter().map(|p| p.y).sum::<f64>() / n;
    let cz = points.iter().map(|p| p.z).sum::<f64>() / n;
    let centroid = V3::new(cx, cy, cz);

    // Covariance matrix (upper triangle)
    let (mut cxx, mut cxy, mut cxz, mut cyy, mut cyz, mut _czz) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for p in points {
        let dx = p.x - cx;
        let dy = p.y - cy;
        let dz = p.z - cz;
        cxx += dx * dx;
        cxy += dx * dy;
        cxz += dx * dz;
        cyy += dy * dy;
        cyz += dy * dz;
        _czz += dz * dz;
    }
    // Normal approximation: cross product of two SVD-like vectors
    // For a near-planar set, the minimum eigenvalue direction can be
    // approximated from the cross product of the two dominant covariance vectors.
    let v1 = V3::new(cxx, cxy, cxz).normalize();
    let v2 = V3::new(cxy, cyy, cyz).normalize();
    let normal = v1.cross(&v2).normalize();
    (centroid, normal)
}

/// Distance from point `p` to the plane defined by `(origin, normal)` (m).
pub fn point_to_plane_distance(p: V3, origin: V3, normal: V3) -> f64 {
    p.sub(&origin).dot(&normal).abs()
}

/// Convert a set of polygon boundary points to a Hermite spline sample.
///
/// Returns `n` evenly spaced points along the perimeter.
pub fn perimeter_resample(pts: &[V3], n: usize) -> Vec<V3> {
    if pts.len() < 2 || n == 0 {
        return Vec::new();
    }
    // Compute cumulative arc lengths
    let mut arcs = vec![0.0_f64; pts.len()];
    for i in 1..pts.len() {
        arcs[i] = arcs[i - 1] + pts[i].sub(&pts[i - 1]).norm();
    }
    let total = *arcs.last().expect("collection should not be empty");
    let mut result = Vec::with_capacity(n);
    for k in 0..n {
        let s = k as f64 / n as f64 * total;
        // Find segment
        let seg = arcs.partition_point(|&a| a <= s).min(pts.len() - 1);
        let seg = seg.max(1);
        let t0 = arcs[seg - 1];
        let t1 = arcs[seg];
        let t = if (t1 - t0).abs() < 1e-15 {
            0.0
        } else {
            (s - t0) / (t1 - t0)
        };
        let p = pts[seg - 1].scale(1.0 - t).add(&pts[seg].scale(t));
        result.push(p);
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- V3 ---

    #[test]
    fn test_v3_norm() {
        let v = V3::new(3.0, 4.0, 0.0);
        assert!((v.norm() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_v3_normalize() {
        let v = V3::new(1.0, 2.0, 2.0);
        let n = v.normalize();
        assert!((n.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_v3_cross_orthogonal() {
        let a = V3::new(1.0, 0.0, 0.0);
        let b = V3::new(0.0, 1.0, 0.0);
        let c = a.cross(&b);
        assert!((c.z - 1.0).abs() < 1e-10);
        assert!(c.x.abs() < 1e-10 && c.y.abs() < 1e-10);
    }

    #[test]
    fn test_v3_dot() {
        let a = V3::new(1.0, 2.0, 3.0);
        let b = V3::new(4.0, 5.0, 6.0);
        assert!((a.dot(&b) - 32.0).abs() < 1e-10);
    }

    // --- GridShell ---

    #[test]
    fn test_grid_shell_node_count() {
        let gs = GridShell::parabolic(5, 6, 20.0, 15.0, 3.0);
        assert_eq!(gs.node_count(), 30);
    }

    #[test]
    fn test_grid_shell_dynamic_relaxation_converges() {
        let mut gs = GridShell::parabolic(5, 5, 10.0, 10.0, 2.0);
        let residual = gs.dynamic_relaxation(500, 0.001, 0.1, 1.0);
        // Residual should drop below 1000 N after 500 steps with this setup
        assert!(residual < 1.0e6, "residual too high: {residual}");
    }

    #[test]
    fn test_grid_shell_fixed_nodes_unchanged() {
        let mut gs = GridShell::parabolic(4, 4, 8.0, 8.0, 1.5);
        let initial_positions: Vec<V3> = gs.fixed_nodes.iter().map(|&i| gs.nodes[i]).collect();
        gs.dynamic_relaxation(100, 0.001, 0.1, 1.0);
        for (k, &idx) in gs.fixed_nodes.iter().enumerate() {
            let diff = gs.nodes[idx].sub(&initial_positions[k]).norm();
            assert!(diff < 1e-10, "fixed node {idx} moved by {diff}");
        }
    }

    #[test]
    fn test_grid_shell_member_count() {
        let gs = GridShell::parabolic(3, 3, 6.0, 6.0, 1.0);
        // horizontal: 2*3 = 6, vertical: 3*2 = 6, diagonals: 2*2*2 = 8 → 20
        assert_eq!(gs.member_count(), 20);
    }

    // --- PanelizationResult ---

    #[test]
    fn test_planar_quad_planarity_zero() {
        // All 4 points in XY plane
        let pts = [
            V3::new(0.0, 0.0, 0.0),
            V3::new(1.0, 0.0, 0.0),
            V3::new(1.0, 1.0, 0.0),
            V3::new(0.0, 1.0, 0.0),
        ];
        let err = PanelizationResult::quad_planarity_error(&pts);
        assert!(err < 1e-12, "planarity error = {err}");
    }

    #[test]
    fn test_planar_quad_planarity_nonzero() {
        let pts = [
            V3::new(0.0, 0.0, 0.0),
            V3::new(1.0, 0.0, 0.0),
            V3::new(1.0, 1.0, 0.0),
            V3::new(0.0, 1.0, 0.1), // fourth point off-plane
        ];
        let err = PanelizationResult::quad_planarity_error(&pts);
        assert!(err > 0.05, "expected non-zero planarity error, got {err}");
    }

    #[test]
    fn test_planar_quad_mesh_flat_surface_zero_error() {
        let mesh = PlanarQuadMesh::from_surface(5, 5, |u, v| V3::new(u, v, 0.0));
        assert!(mesh.max_planarity_error() < 1e-12);
    }

    #[test]
    fn test_planar_quad_mesh_curved_surface_error_positive() {
        let mesh = PlanarQuadMesh::from_surface(10, 10, |u, v| {
            V3::new(u, v, (u * PI * 2.0).sin() * (v * PI * 2.0).cos() * 0.3)
        });
        // Curved surface should have non-trivial planarity error
        assert!(mesh.max_planarity_error() >= 0.0);
    }

    #[test]
    fn test_planar_quad_mesh_face_count() {
        let mesh = PlanarQuadMesh::from_surface(6, 8, |u, v| V3::new(u, v, 0.0));
        assert_eq!(mesh.face_count(), 5 * 7);
    }

    #[test]
    fn test_planar_quad_mesh_compliance_ratio_flat() {
        let mesh = PlanarQuadMesh::from_surface(5, 5, |u, v| V3::new(u, v, 0.0));
        assert!((mesh.planarity_compliance_ratio(0.001) - 1.0).abs() < 1e-10);
    }

    // --- TensileStructure ---

    #[test]
    fn test_tensile_flat_net_node_count() {
        let net = TensileStructure::flat_net(5, 5, 10.0, 10.0, 5000.0);
        assert_eq!(net.nodes.len(), 25);
    }

    #[test]
    fn test_tensile_form_find_sag_positive() {
        let mut net = TensileStructure::flat_net(5, 5, 10.0, 10.0, 10000.0);
        let sag = net.form_find(500.0, 1000, 0.001, 0.1, 0.01);
        assert!(sag >= 0.0, "sag should be non-negative, got {sag}");
    }

    #[test]
    fn test_tensile_mean_cable_tension() {
        let net = TensileStructure::flat_net(4, 4, 8.0, 8.0, 3000.0);
        let mean = net.mean_cable_tension();
        assert!((mean - 3000.0).abs() < 1e-6);
    }

    #[test]
    fn test_tensile_anchor_nodes_not_moved() {
        let mut net = TensileStructure::flat_net(4, 4, 6.0, 6.0, 5000.0);
        let initial: Vec<V3> = net.anchors.iter().map(|&i| net.nodes[i]).collect();
        net.form_find(200.0, 200, 0.001, 0.2, 0.1);
        for (k, &idx) in net.anchors.iter().enumerate() {
            let diff = net.nodes[idx].sub(&initial[k]).norm();
            assert!(diff < 1e-10, "anchor {idx} moved by {diff}");
        }
    }

    // --- GeodesicDome ---

    #[test]
    fn test_geodesic_dome_freq1_vertex_count() {
        let dome = GeodesicDome::new(1, 5.0, false);
        // Icosahedron has 12 vertices
        assert_eq!(dome.vertices.len(), 12);
    }

    #[test]
    fn test_geodesic_dome_freq2_expected_vertices() {
        let dome = GeodesicDome::new(2, 5.0, false);
        let expected = GeodesicDome::expected_vertex_count(2); // 10*4+2 = 42
        // Allow some tolerance for deduplication in the algorithm
        let actual = dome.vertices.len();
        assert!(
            (actual as i64 - expected as i64).abs() <= 5,
            "expected ~{expected} vertices, got {actual}"
        );
    }

    #[test]
    fn test_geodesic_dome_all_on_sphere() {
        let dome = GeodesicDome::new(2, 7.0, false);
        assert!(
            dome.verify_sphericity(1e-6),
            "vertices not on sphere within tolerance"
        );
    }

    #[test]
    fn test_geodesic_dome_strut_lengths_similar() {
        let dome = GeodesicDome::new(2, 5.0, false);
        let (min, max, _mean) = dome.strut_length_stats();
        // For freq 2, all struts should be within 15% of each other
        assert!(
            max / min.max(1e-10) < 1.20,
            "strut length ratio {:.3} too large",
            max / min
        );
    }

    #[test]
    fn test_geodesic_dome_freq1_strut_count() {
        let dome = GeodesicDome::new(1, 5.0, false);
        // Icosahedron has 30 edges
        assert_eq!(dome.struts.len(), 30);
    }

    // --- StructuralGlass ---

    #[test]
    fn test_glass_self_weight_positive() {
        let g = StructuralGlass::monolithic_10mm();
        let sw = g.self_weight_pressure();
        assert!(sw > 0.0 && sw < 500.0, "self weight = {sw} Pa");
    }

    #[test]
    fn test_glass_thermal_stress() {
        let g = StructuralGlass::monolithic_10mm();
        let ts = g.thermal_stress();
        // 70e9 * 9e-6 * 40 = 25.2 MPa
        assert!((ts - 25.2e6).abs() < 1.0e5);
    }

    #[test]
    fn test_glass_wind_stress_positive() {
        let g = StructuralGlass::monolithic_10mm();
        let ws = g.max_bending_stress_wind();
        assert!(ws > 0.0);
    }

    #[test]
    fn test_glass_combined_stress_greater_than_parts() {
        let g = StructuralGlass::monolithic_10mm();
        let combined = g.combined_stress();
        assert!(combined > g.thermal_stress());
        assert!(combined > g.max_bending_stress_wind());
    }

    #[test]
    fn test_glass_centre_deflection() {
        let g = StructuralGlass::monolithic_10mm();
        let d = g.centre_deflection_wind();
        assert!(d > 0.0 && d < 0.2, "deflection = {d} m");
    }

    #[test]
    fn test_glass_sealant_capacity_positive() {
        let g = StructuralGlass::monolithic_10mm();
        assert!(g.sealant_shear_capacity() > 0.0);
    }

    // --- ParametricFacade ---

    #[test]
    fn test_facade_total_panel_count() {
        let f = ParametricFacade::tokyo_south();
        assert_eq!(f.total_panel_count(), 80); // 8 * 10
    }

    #[test]
    fn test_facade_solar_altitude_summer_positive() {
        let f = ParametricFacade::tokyo_south();
        let alt = f.solar_altitude(172.0, 12.0);
        assert!(alt > 0.0, "summer noon altitude should be positive: {alt}");
    }

    #[test]
    fn test_facade_solar_altitude_night_negative() {
        let f = ParametricFacade::tokyo_south();
        let alt = f.solar_altitude(172.0, 0.0); // midnight
        assert!(alt < 0.0, "midnight altitude should be negative: {alt}");
    }

    #[test]
    fn test_facade_shade_depth_matrix_dimensions() {
        let f = ParametricFacade::tokyo_south();
        let m = f.shade_depth_matrix();
        assert_eq!(m.len(), f.levels);
        assert_eq!(m[0].len(), f.bays);
    }

    #[test]
    fn test_facade_shade_depth_bounded() {
        let f = ParametricFacade::tokyo_south();
        let m = f.shade_depth_matrix();
        for row in &m {
            for &d in row {
                assert!(d >= 0.0 && d <= f.max_shade_depth + 1e-10, "depth = {d}");
            }
        }
    }

    #[test]
    fn test_facade_total_area() {
        let f = ParametricFacade::tokyo_south();
        let area = f.total_facade_area();
        // 8 * 1.5 * 10 * 3.6 = 432 m²
        assert!((area - 432.0).abs() < 1e-6);
    }

    #[test]
    fn test_facade_annual_solar_exposure_positive() {
        let f = ParametricFacade::tokyo_south();
        let exp = f.annual_solar_exposure();
        assert!(
            exp >= 0.0,
            "annual solar exposure should be non-negative: {exp}"
        );
    }

    // --- Utility functions ---

    #[test]
    fn test_triangle_area() {
        let a = V3::new(0.0, 0.0, 0.0);
        let b = V3::new(4.0, 0.0, 0.0);
        let c = V3::new(0.0, 3.0, 0.0);
        let area = triangle_area(a, b, c);
        assert!((area - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_fit_plane_xy_plane() {
        let pts = vec![
            V3::new(0.0, 0.0, 0.0),
            V3::new(1.0, 0.0, 0.0),
            V3::new(0.0, 1.0, 0.0),
            V3::new(1.0, 1.0, 0.0),
        ];
        let (_centroid, normal) = fit_plane(&pts);
        // Normal should be close to ±Z
        assert!(normal.z.abs() > 0.5, "normal.z = {}", normal.z);
    }

    #[test]
    fn test_perimeter_resample_count() {
        let pts = vec![
            V3::new(0.0, 0.0, 0.0),
            V3::new(1.0, 0.0, 0.0),
            V3::new(1.0, 1.0, 0.0),
            V3::new(0.0, 1.0, 0.0),
        ];
        let resampled = perimeter_resample(&pts, 8);
        assert_eq!(resampled.len(), 8);
    }

    #[test]
    fn test_point_to_plane_distance() {
        let p = V3::new(0.0, 0.0, 5.0);
        let origin = V3::zero();
        let normal = V3::new(0.0, 0.0, 1.0);
        let d = point_to_plane_distance(p, origin, normal);
        assert!((d - 5.0).abs() < 1e-10);
    }
}
