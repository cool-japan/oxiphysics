// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cosserat rod / beam soft-body simulation.
//!
//! Provides a geometrically exact Cosserat rod model including:
//!
//! - Cosserat rod theory: centroid curve + material frame (director triad)
//! - Darboux vector, intrinsic curvature and twist
//! - Kirchhoff rod kinematics (inextensible/unshearable approximation)
//! - Simo–Reissner rod (extensible, shearable)
//! - Discrete Cosserat rod: chain of rigid frames connected by elastic joints
//! - Strain measures: stretching, shearing, bending, twisting
//! - Elastic energy (Kirchhoff model: bending + twist stiffness)
//! - Geometric exactness: finite rotations via rotation matrices
//! - Buckling load estimation (Euler buckling extended to rods)
//! - Tendon-driven rod actuation
//! - Active Cosserat rod (muscle-like distributed moment)
//! - Cosserat plate: two-dimensional generalisation
//! - Wrinkle simulation via wrinkling criterion

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers (no nalgebra – plain [f64; 3] and [f64; 9])
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the cross product a × b.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector.
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Normalise a 3-vector; returns the zero vector if nearly zero.
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}

/// Add two 3-vectors.
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract b from a.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by scalar s.
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// 3×3 identity matrix stored column-major as \[f64; 9\].
fn mat3_identity() -> [f64; 9] {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}

/// Matrix-vector product R·v for R stored row-major as \[f64; 9\].
fn mat3_mul_vec(r: [f64; 9], v: [f64; 3]) -> [f64; 3] {
    [
        r[0] * v[0] + r[1] * v[1] + r[2] * v[2],
        r[3] * v[0] + r[4] * v[1] + r[5] * v[2],
        r[6] * v[0] + r[7] * v[1] + r[8] * v[2],
    ]
}

/// Matrix-matrix product A·B (both row-major \[f64; 9\]).
fn mat3_mul(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}

/// Transpose of a 3×3 matrix.
fn mat3_transpose(a: [f64; 9]) -> [f64; 9] {
    [a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8]]
}

/// Rodrigues rotation matrix: rotate by axis `axis` (unit vector) by angle `theta` \[rad\].
///
/// R = I + sin θ · \[n\]× + (1 − cos θ) · \[n\]×²
pub fn rodrigues(axis: [f64; 3], theta: f64) -> [f64; 9] {
    let n = normalize3(axis);
    let (nx, ny, nz) = (n[0], n[1], n[2]);
    let (s, c) = (theta.sin(), theta.cos());
    let t = 1.0 - c;
    [
        c + nx * nx * t,
        nx * ny * t - nz * s,
        nx * nz * t + ny * s,
        ny * nx * t + nz * s,
        c + ny * ny * t,
        ny * nz * t - nx * s,
        nz * nx * t - ny * s,
        nz * ny * t + nx * s,
        c + nz * nz * t,
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// DirectorTriad: material frame (d1, d2, d3)
// ─────────────────────────────────────────────────────────────────────────────

/// Material frame (director triad) of a Cosserat rod cross-section.
///
/// * `d1`, `d2` – Cross-section principal axes (orthogonal unit vectors).
/// * `d3` – Tangent vector (unit vector along rod centreline).
#[derive(Debug, Clone, Copy)]
pub struct DirectorTriad {
    /// First cross-section principal axis d₁.
    pub d1: [f64; 3],
    /// Second cross-section principal axis d₂.
    pub d2: [f64; 3],
    /// Tangent director d₃ (along centreline).
    pub d3: [f64; 3],
}

impl DirectorTriad {
    /// Create a new [`DirectorTriad`] from three orthogonal unit vectors.
    pub fn new(d1: [f64; 3], d2: [f64; 3], d3: [f64; 3]) -> Self {
        Self { d1, d2, d3 }
    }

    /// Canonical frame aligned with the z-axis (d3 = \[0,0,1\]).
    pub fn canonical() -> Self {
        Self {
            d1: [1.0, 0.0, 0.0],
            d2: [0.0, 1.0, 0.0],
            d3: [0.0, 0.0, 1.0],
        }
    }

    /// Convert to rotation matrix R whose columns are d1, d2, d3.
    pub fn to_rotation_matrix(&self) -> [f64; 9] {
        [
            self.d1[0], self.d2[0], self.d3[0], self.d1[1], self.d2[1], self.d3[1], self.d1[2],
            self.d2[2], self.d3[2],
        ]
    }

    /// Apply a Rodrigues rotation `R` to all three directors.
    pub fn rotate(&self, r: [f64; 9]) -> Self {
        Self {
            d1: mat3_mul_vec(r, self.d1),
            d2: mat3_mul_vec(r, self.d2),
            d3: mat3_mul_vec(r, self.d3),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Darboux vector and strain measures
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Darboux vector κ = (κ₁, κ₂, τ) from adjacent director triads.
///
/// The Darboux vector encodes the curvatures κ₁, κ₂ and twist τ of the rod.
/// Given two adjacent frames at arc-length positions s and s+Δs:
///
/// ω ≈ (R(s)ᵀ · (R(s+Δs) − R(s))) / Δs  mapped to its axial vector.
///
/// Returns `[kappa1, kappa2, tau]` in the material frame.
///
/// # Arguments
/// * `frame_a` – Director triad at s.
/// * `frame_b` – Director triad at s+Δs.
/// * `ds` – Arc-length increment Δs \[m\].
pub fn darboux_vector(frame_a: &DirectorTriad, frame_b: &DirectorTriad, ds: f64) -> [f64; 3] {
    // Relative rotation dR = Ra^T * Rb.
    let ra = frame_a.to_rotation_matrix();
    let rb = frame_b.to_rotation_matrix();
    let ra_t = mat3_transpose(ra);
    let dr = mat3_mul(ra_t, rb);
    // Extract axial vector of skew-symmetric part of dR: (dR - I).
    // ω × = (dR - I) / ds  → axial: [dr[7]-dr[5], dr[2]-dr[6], dr[3]-dr[1]] / (2*ds)
    let wx = (dr[7] - dr[5]) / (2.0 * ds);
    let wy = (dr[2] - dr[6]) / (2.0 * ds);
    let wz = (dr[3] - dr[1]) / (2.0 * ds);
    [wx, wy, wz]
}

/// Stretch-shear strain vector γ = Rᵀ · ∂r/∂s − e₃.
///
/// For an inextensible, unshearable rod (Kirchhoff), γ = 0.
/// For a general Cosserat rod, γ ≠ 0 when there is stretching (γ₃ ≠ 0)
/// or shearing (γ₁, γ₂ ≠ 0).
///
/// # Arguments
/// * `frame` – Material frame at this point.
/// * `dr_ds` – Derivative of centroid position w.r.t. arc-length ∂r/∂s.
pub fn stretch_shear_strain(frame: &DirectorTriad, dr_ds: [f64; 3]) -> [f64; 3] {
    let ra = frame.to_rotation_matrix();
    let ra_t = mat3_transpose(ra);
    let gamma = mat3_mul_vec(ra_t, dr_ds);
    // Subtract reference e3 = [0,0,1].
    [gamma[0], gamma[1], gamma[2] - 1.0]
}

// ─────────────────────────────────────────────────────────────────────────────
// CosseratRodProperties
// ─────────────────────────────────────────────────────────────────────────────

/// Material and geometric properties of a Cosserat rod.
#[derive(Debug, Clone, Copy)]
pub struct CosseratRodProperties {
    /// Bending stiffness EI₁ about d₁ \[N·m²\].
    pub ei1: f64,
    /// Bending stiffness EI₂ about d₂ \[N·m²\].
    pub ei2: f64,
    /// Torsional stiffness GJ \[N·m²\].
    pub gj: f64,
    /// Axial stiffness EA \[N\].
    pub ea: f64,
    /// Shear stiffness GA₁ \[N\].
    pub ga1: f64,
    /// Shear stiffness GA₂ \[N\].
    pub ga2: f64,
    /// Linear mass density \[kg/m\].
    pub rho_a: f64,
}

impl CosseratRodProperties {
    /// Create rod properties for a circular cross-section rod.
    ///
    /// # Arguments
    /// * `radius` – Cross-section radius \[m\].
    /// * `young_modulus` – Young's modulus E \[Pa\].
    /// * `shear_modulus` – Shear modulus G \[Pa\].
    /// * `density` – Material density \[kg/m³\].
    pub fn circular(radius: f64, young_modulus: f64, shear_modulus: f64, density: f64) -> Self {
        let r2 = radius * radius;
        let area = PI * r2;
        let i_bending = PI * r2 * r2 / 4.0;
        let j_polar = PI * r2 * r2 / 2.0;
        Self {
            ei1: young_modulus * i_bending,
            ei2: young_modulus * i_bending,
            gj: shear_modulus * j_polar,
            ea: young_modulus * area,
            ga1: shear_modulus * area,
            ga2: shear_modulus * area,
            rho_a: density * area,
        }
    }

    /// Create rod properties for a rectangular cross-section rod.
    ///
    /// # Arguments
    /// * `width` – Cross-section width b \[m\].
    /// * `height` – Cross-section height h \[m\].
    /// * `young_modulus` – Young's modulus E \[Pa\].
    /// * `shear_modulus` – Shear modulus G \[Pa\].
    /// * `density` – Material density \[kg/m³\].
    pub fn rectangular(
        width: f64,
        height: f64,
        young_modulus: f64,
        shear_modulus: f64,
        density: f64,
    ) -> Self {
        let area = width * height;
        let i1 = width * height * height * height / 12.0;
        let i2 = height * width * width * width / 12.0;
        // Approximate torsion constant for rectangle (Saint-Venant).
        let (a, b) = if width >= height {
            (width, height)
        } else {
            (height, width)
        };
        let j = a
            * b
            * b
            * b
            * (1.0 / 3.0 - 0.21 * b / a * (1.0 - b * b * b * b / (12.0 * a * a * a * a)));
        Self {
            ei1: young_modulus * i1,
            ei2: young_modulus * i2,
            gj: shear_modulus * j,
            ea: young_modulus * area,
            ga1: shear_modulus * area,
            ga2: shear_modulus * area,
            rho_a: density * area,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Elastic energy: Kirchhoff rod model
// ─────────────────────────────────────────────────────────────────────────────

/// Kirchhoff rod elastic energy density per unit arc-length.
///
/// E_elastic = ½·(EI₁·(κ₁ − κ₁₀)² + EI₂·(κ₂ − κ₂₀)² + GJ·(τ − τ₀)²)
///
/// where (κ₁₀, κ₂₀, τ₀) are the intrinsic (rest) curvatures/twist.
///
/// # Arguments
/// * `props` – Rod material properties.
/// * `kappa` – Current Darboux vector \[κ₁, κ₂, τ\].
/// * `kappa0` – Intrinsic Darboux vector (rest curvature/twist).
pub fn kirchhoff_energy_density(
    props: &CosseratRodProperties,
    kappa: [f64; 3],
    kappa0: [f64; 3],
) -> f64 {
    let dk1 = kappa[0] - kappa0[0];
    let dk2 = kappa[1] - kappa0[1];
    let dtau = kappa[2] - kappa0[2];
    0.5 * (props.ei1 * dk1 * dk1 + props.ei2 * dk2 * dk2 + props.gj * dtau * dtau)
}

/// Simo–Reissner rod elastic energy density (extensible + shearable).
///
/// E = ½·(EA·γ₃² + GA₁·γ₁² + GA₂·γ₂²) + ½·(EI₁·κ₁² + EI₂·κ₂² + GJ·τ²)
///
/// # Arguments
/// * `props` – Rod material properties.
/// * `gamma` – Stretch-shear strain \[γ₁, γ₂, γ₃\].
/// * `kappa` – Bending-twist strain \[κ₁, κ₂, τ\].
pub fn simo_reissner_energy_density(
    props: &CosseratRodProperties,
    gamma: [f64; 3],
    kappa: [f64; 3],
) -> f64 {
    let e_shear_stretch = 0.5
        * (props.ga1 * gamma[0] * gamma[0]
            + props.ga2 * gamma[1] * gamma[1]
            + props.ea * gamma[2] * gamma[2]);
    let e_bend_twist = 0.5
        * (props.ei1 * kappa[0] * kappa[0]
            + props.ei2 * kappa[1] * kappa[1]
            + props.gj * kappa[2] * kappa[2]);
    e_shear_stretch + e_bend_twist
}

/// Internal moment vector M from Darboux vector and stiffness.
///
/// M = (EI₁·κ₁, EI₂·κ₂, GJ·τ)
///
/// # Arguments
/// * `props` – Rod material properties.
/// * `kappa` – Darboux vector \[κ₁, κ₂, τ\].
/// * `kappa0` – Intrinsic Darboux vector.
pub fn internal_moment(
    props: &CosseratRodProperties,
    kappa: [f64; 3],
    kappa0: [f64; 3],
) -> [f64; 3] {
    [
        props.ei1 * (kappa[0] - kappa0[0]),
        props.ei2 * (kappa[1] - kappa0[1]),
        props.gj * (kappa[2] - kappa0[2]),
    ]
}

/// Internal force vector N from stretch-shear strain.
///
/// N = (GA₁·γ₁, GA₂·γ₂, EA·γ₃)
///
/// # Arguments
/// * `props` – Rod material properties.
/// * `gamma` – Stretch-shear strain.
pub fn internal_force(props: &CosseratRodProperties, gamma: [f64; 3]) -> [f64; 3] {
    [
        props.ga1 * gamma[0],
        props.ga2 * gamma[1],
        props.ea * gamma[2],
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// DiscreteCosseratNode: one segment of the discrete Cosserat rod
// ─────────────────────────────────────────────────────────────────────────────

/// A single node in the discrete Cosserat rod model.
///
/// Each node stores the centroid position and material frame of one segment,
/// along with dynamic state (linear and angular velocities) for simulation.
#[derive(Debug, Clone)]
pub struct CosseratNode {
    /// Centroid position \[m\].
    pub position: [f64; 3],
    /// Material frame (director triad).
    pub frame: DirectorTriad,
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\] in global frame.
    pub omega: [f64; 3],
    /// Is this node pinned (static)?
    pub is_pinned: bool,
    /// Accumulated external force \[N\].
    pub ext_force: [f64; 3],
    /// Accumulated external moment \[N·m\].
    pub ext_moment: [f64; 3],
}

impl CosseratNode {
    /// Create a new free [`CosseratNode`] at position with canonical frame.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            frame: DirectorTriad::canonical(),
            velocity: [0.0; 3],
            omega: [0.0; 3],
            is_pinned: false,
            ext_force: [0.0; 3],
            ext_moment: [0.0; 3],
        }
    }

    /// Create a pinned (static) [`CosseratNode`].
    pub fn new_pinned(position: [f64; 3]) -> Self {
        let mut node = Self::new(position);
        node.is_pinned = true;
        node
    }

    /// Apply an external force to this node.
    pub fn apply_force(&mut self, force: [f64; 3]) {
        self.ext_force = add3(self.ext_force, force);
    }

    /// Apply an external moment to this node.
    pub fn apply_moment(&mut self, moment: [f64; 3]) {
        self.ext_moment = add3(self.ext_moment, moment);
    }

    /// Clear accumulated forces and moments.
    pub fn clear_loads(&mut self) {
        self.ext_force = [0.0; 3];
        self.ext_moment = [0.0; 3];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DiscreteCosseratRod: chain of rigid frames
// ─────────────────────────────────────────────────────────────────────────────

/// Discrete Cosserat rod: a chain of rigid frames connected by elastic joints.
///
/// Each segment is defined by two adjacent nodes. The elastic restoring moments
/// and forces follow the Kirchhoff (or Simo–Reissner) constitutive law.
#[derive(Debug, Clone)]
pub struct DiscreteCosseratRod {
    /// Nodes of the discrete rod (N+1 nodes for N segments).
    pub nodes: Vec<CosseratNode>,
    /// Intrinsic Darboux vector at each segment (rest curvature/twist).
    pub intrinsic_kappa: Vec<[f64; 3]>,
    /// Rest segment lengths \[m\].
    pub rest_lengths: Vec<f64>,
    /// Uniform material properties.
    pub props: CosseratRodProperties,
}

impl DiscreteCosseratRod {
    /// Create a straight discrete Cosserat rod along the z-axis.
    ///
    /// # Arguments
    /// * `num_segments` – Number of segments N (creates N+1 nodes).
    /// * `total_length` – Total rest length L \[m\].
    /// * `props` – Rod material properties.
    pub fn straight(num_segments: usize, total_length: f64, props: CosseratRodProperties) -> Self {
        let seg_len = total_length / num_segments as f64;
        let mut nodes = Vec::with_capacity(num_segments + 1);
        for i in 0..=num_segments {
            let z = i as f64 * seg_len;
            nodes.push(CosseratNode::new([0.0, 0.0, z]));
        }
        // Pin the first node.
        nodes[0].is_pinned = true;
        let intrinsic_kappa = vec![[0.0; 3]; num_segments];
        let rest_lengths = vec![seg_len; num_segments];
        Self {
            nodes,
            intrinsic_kappa,
            rest_lengths,
            props,
        }
    }

    /// Number of segments.
    pub fn num_segments(&self) -> usize {
        self.nodes.len().saturating_sub(1)
    }

    /// Compute Darboux vector for segment `i`.
    pub fn darboux(&self, i: usize) -> [f64; 3] {
        let ds = self.rest_lengths[i];
        darboux_vector(&self.nodes[i].frame, &self.nodes[i + 1].frame, ds)
    }

    /// Compute elastic energy of the full rod.
    pub fn elastic_energy(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.num_segments() {
            let kappa = self.darboux(i);
            let kappa0 = self.intrinsic_kappa[i];
            let ds = self.rest_lengths[i];
            energy += kirchhoff_energy_density(&self.props, kappa, kappa0) * ds;
        }
        energy
    }

    /// Apply gravity to all free nodes.
    ///
    /// # Arguments
    /// * `g` – Gravitational acceleration vector \[m/s²\].
    pub fn apply_gravity(&mut self, g: [f64; 3]) {
        let n_seg = self.num_segments() as f64;
        let total_len: f64 = self.rest_lengths.iter().sum();
        let seg_mass = self.props.rho_a * total_len / n_seg;
        for node in &mut self.nodes {
            if !node.is_pinned {
                let gravity_force = scale3(g, seg_mass);
                node.apply_force(gravity_force);
            }
        }
    }

    /// Semi-implicit Euler integration step.
    ///
    /// Computes elastic restoring forces/moments from curvature mismatch
    /// and integrates the equations of motion.
    ///
    /// # Arguments
    /// * `dt` – Time step \[s\].
    /// * `damping` – Velocity damping coefficient (0=none, 1=critical).
    pub fn step(&mut self, dt: f64, damping: f64) {
        let n_seg = self.num_segments();

        // Compute elastic restoring moments at each joint.
        let mut elastic_moments = vec![[0.0_f64; 3]; self.nodes.len()];
        let mut elastic_forces = vec![[0.0_f64; 3]; self.nodes.len()];

        for i in 0..n_seg {
            let kappa = self.darboux(i);
            let kappa0 = self.intrinsic_kappa[i];
            let ds = self.rest_lengths[i];
            let moment_local = internal_moment(&self.props, kappa, kappa0);
            // Rotate to global frame using node i's frame.
            let r_i = self.nodes[i].frame.to_rotation_matrix();
            let moment_global = mat3_mul_vec(r_i, moment_local);
            // Distribute moment to adjacent nodes.
            let half_m = scale3(moment_global, 0.5);
            elastic_moments[i] = add3(elastic_moments[i], half_m);
            elastic_moments[i + 1] = add3(elastic_moments[i + 1], half_m);

            // Compute centroid stretch (restoring force for Simo–Reissner).
            let pos_i = self.nodes[i].position;
            let pos_j = self.nodes[i + 1].position;
            let dr = sub3(pos_j, pos_i);
            let dr_norm = norm3(dr);
            if dr_norm > 1e-30 {
                let stretch = (dr_norm - ds) / ds;
                let f_mag = self.props.ea * stretch;
                let f_dir = scale3(normalize3(dr), f_mag);
                elastic_forces[i] = add3(elastic_forces[i], f_dir);
                elastic_forces[i + 1] = sub3(elastic_forces[i + 1], f_dir);
            }
        }

        // Integrate velocities and positions.
        let inv_mass =
            1.0 / (self.props.rho_a * self.rest_lengths.iter().sum::<f64>() / n_seg as f64);

        for i in 0..self.nodes.len() {
            if self.nodes[i].is_pinned {
                continue;
            }
            let total_force = add3(add3(self.nodes[i].ext_force, elastic_forces[i]), [0.0; 3]);
            let accel = scale3(total_force, inv_mass);
            // Velocity update with damping.
            let v_new = scale3(
                add3(self.nodes[i].velocity, scale3(accel, dt)),
                1.0 - damping * dt,
            );
            self.nodes[i].velocity = v_new;
            self.nodes[i].position = add3(self.nodes[i].position, scale3(v_new, dt));

            // Angular velocity and frame update.
            let total_moment = add3(self.nodes[i].ext_moment, elastic_moments[i]);
            let alpha = scale3(total_moment, inv_mass);
            let omega_new = scale3(
                add3(self.nodes[i].omega, scale3(alpha, dt)),
                1.0 - damping * dt,
            );
            self.nodes[i].omega = omega_new;
            let omega_mag = norm3(omega_new);
            if omega_mag > 1e-30 {
                let rot = rodrigues(omega_new, omega_mag * dt);
                self.nodes[i].frame = self.nodes[i].frame.rotate(rot);
            }

            self.nodes[i].clear_loads();
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Buckling of rods
// ─────────────────────────────────────────────────────────────────────────────

/// Euler critical buckling load for a Cosserat rod.
///
/// P_cr = π² · EI / (K·L)²
///
/// # Arguments
/// * `ei` – Bending stiffness EI \[N·m²\].
/// * `length` – Rod length L \[m\].
/// * `effective_length_factor` – Effective length factor K
///   (1.0 = pinned-pinned, 0.5 = fixed-fixed, 0.7 = fixed-pinned, 2.0 = fixed-free).
pub fn euler_buckling_load(ei: f64, length: f64, effective_length_factor: f64) -> f64 {
    let le = effective_length_factor * length;
    PI * PI * ei / (le * le)
}

/// Slenderness ratio λ = L_eff / r_gyration.
///
/// A rod buckles elastically (Euler) for λ > critical slenderness.
///
/// # Arguments
/// * `effective_length` – Effective buckling length K·L \[m\].
/// * `radius_of_gyration` – r = sqrt(I/A) \[m\].
pub fn slenderness_ratio(effective_length: f64, radius_of_gyration: f64) -> f64 {
    effective_length / radius_of_gyration
}

/// Rank of a rod cross-section: radius of gyration r = sqrt(I/A).
///
/// # Arguments
/// * `second_moment_area` – I \[m⁴\].
/// * `area` – A \[m²\].
pub fn radius_of_gyration(second_moment_area: f64, area: f64) -> f64 {
    (second_moment_area / area).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tendon-driven rod actuation
// ─────────────────────────────────────────────────────────────────────────────

/// Tendon-driven rod: compute bending moment from tendon tension and eccentricity.
///
/// A tendon running parallel to the rod axis at eccentricity `e` from the neutral axis
/// creates a bending moment M = T · e when tensioned.
///
/// # Arguments
/// * `tension` – Tendon tension T \[N\].
/// * `eccentricity` – Perpendicular distance e from neutral axis \[m\].
pub fn tendon_bending_moment(tension: f64, eccentricity: f64) -> f64 {
    tension * eccentricity
}

/// Effective curvature change κ_act induced by tendon actuation.
///
/// κ_act = M_act / EI = (T · e) / EI
///
/// # Arguments
/// * `tension` – Tendon tension T \[N\].
/// * `eccentricity` – Eccentricity e \[m\].
/// * `ei` – Bending stiffness EI \[N·m²\].
pub fn tendon_induced_curvature(tension: f64, eccentricity: f64, ei: f64) -> f64 {
    tension * eccentricity / ei
}

/// TendonActuatedRod: a rod with a set of tendons for position control.
#[derive(Debug, Clone)]
pub struct TendonActuatedRod {
    /// Underlying discrete Cosserat rod.
    pub rod: DiscreteCosseratRod,
    /// Tendon eccentricities (d1 direction, d2 direction) for each tendon.
    pub tendon_eccentricities: Vec<[f64; 2]>,
    /// Current tendon tensions \[N\].
    pub tendon_tensions: Vec<f64>,
}

impl TendonActuatedRod {
    /// Create a new [`TendonActuatedRod`] from a discrete rod and tendon layout.
    pub fn new(rod: DiscreteCosseratRod, tendon_eccentricities: Vec<[f64; 2]>) -> Self {
        let n = tendon_eccentricities.len();
        Self {
            rod,
            tendon_eccentricities,
            tendon_tensions: vec![0.0; n],
        }
    }

    /// Set tendon tension for tendon `idx`.
    pub fn set_tension(&mut self, idx: usize, tension: f64) {
        if idx < self.tendon_tensions.len() {
            self.tendon_tensions[idx] = tension;
        }
    }

    /// Compute the net actuation bending moment per unit length at each node.
    ///
    /// Returns `[m_d1, m_d2]` moments (about d1 and d2 axes respectively).
    pub fn actuation_moment(&self) -> [f64; 2] {
        let mut m1 = 0.0;
        let mut m2 = 0.0;
        for (i, ecc) in self.tendon_eccentricities.iter().enumerate() {
            let t = self.tendon_tensions[i];
            m1 += t * ecc[1]; // tension × e2 → bending about d1
            m2 += t * ecc[0]; // tension × e1 → bending about d2
        }
        [m1, m2]
    }

    /// Apply tendon actuation moments to all rod nodes and step.
    ///
    /// # Arguments
    /// * `dt` – Time step \[s\].
    /// * `damping` – Damping coefficient.
    pub fn step(&mut self, dt: f64, damping: f64) {
        let act_m = self.actuation_moment();
        let n_nodes = self.rod.nodes.len();
        for i in 0..n_nodes {
            if !self.rod.nodes[i].is_pinned {
                // Project actuation moments into global frame via current node frame.
                let r = self.rod.nodes[i].frame.to_rotation_matrix();
                // m_d1 acts about d1, m_d2 about d2.
                let m_local = [act_m[0], act_m[1], 0.0];
                let m_global = mat3_mul_vec(r, m_local);
                self.rod.nodes[i].apply_moment(m_global);
            }
        }
        self.rod.step(dt, damping);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Active Cosserat rod (muscle-like distributed moment)
// ─────────────────────────────────────────────────────────────────────────────

/// Active Cosserat rod: a rod with intrinsic curvature that can be varied
/// to simulate muscle-like actuation.
///
/// By changing the intrinsic (rest) curvature κ₀(s, t) over time, the rod
/// bends like a biological muscle or soft actuator.
#[derive(Debug, Clone)]
pub struct ActiveCosseratRod {
    /// Underlying discrete Cosserat rod.
    pub rod: DiscreteCosseratRod,
    /// Target intrinsic curvature profiles κ₀_target for each segment.
    pub target_kappa: Vec<[f64; 3]>,
    /// Activation time constant τ_act \[s\].
    pub tau_activation: f64,
}

impl ActiveCosseratRod {
    /// Create a new [`ActiveCosseratRod`].
    pub fn new(rod: DiscreteCosseratRod, tau_activation: f64) -> Self {
        let n_seg = rod.num_segments();
        let target_kappa = rod.intrinsic_kappa.clone();
        Self {
            rod,
            target_kappa: target_kappa.into_iter().take(n_seg).collect(),
            tau_activation,
        }
    }

    /// Set the target curvature for segment `idx`.
    pub fn set_target_curvature(&mut self, idx: usize, kappa_target: [f64; 3]) {
        if idx < self.target_kappa.len() {
            self.target_kappa[idx] = kappa_target;
        }
    }

    /// Update intrinsic curvature towards target with first-order dynamics.
    ///
    /// κ₀(t + Δt) = κ₀(t) + (κ_target − κ₀(t)) / τ_act · Δt
    pub fn update_activation(&mut self, dt: f64) {
        for i in 0..self.rod.intrinsic_kappa.len() {
            let kappa0 = self.rod.intrinsic_kappa[i];
            let target = self.target_kappa[i];
            let dkappa = scale3(sub3(target, kappa0), dt / self.tau_activation);
            self.rod.intrinsic_kappa[i] = add3(kappa0, dkappa);
        }
    }

    /// Step: update activation and integrate rod dynamics.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3], damping: f64) {
        self.update_activation(dt);
        self.rod.apply_gravity(gravity);
        self.rod.step(dt, damping);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cosserat plate theory
// ─────────────────────────────────────────────────────────────────────────────

/// Cosserat plate: a two-dimensional Cosserat body with in-plane and bending DOF.
///
/// Stores an M×N grid of position + director nodes.
#[derive(Debug, Clone)]
pub struct CosseratPlate {
    /// Grid dimensions (rows × cols).
    pub rows: usize,
    /// Grid columns.
    pub cols: usize,
    /// Node positions \[m\].
    pub positions: Vec<[f64; 3]>,
    /// Node director triads.
    pub frames: Vec<DirectorTriad>,
    /// Rest spacing in row direction \[m\].
    pub dx: f64,
    /// Rest spacing in column direction \[m\].
    pub dy: f64,
    /// Bending stiffness D = E·h³ / (12·(1−ν²)) \[N·m\].
    pub bending_stiffness: f64,
    /// Membrane stiffness K = E·h / (1−ν²) \[N/m\].
    pub membrane_stiffness: f64,
}

impl CosseratPlate {
    /// Create a flat Cosserat plate in the XY plane.
    ///
    /// # Arguments
    /// * `rows` – Number of nodes in the y direction.
    /// * `cols` – Number of nodes in the x direction.
    /// * `lx` – Total plate length in x \[m\].
    /// * `ly` – Total plate length in y \[m\].
    /// * `bending_stiffness` – D \[N·m\].
    /// * `membrane_stiffness` – K \[N/m\].
    pub fn flat(
        rows: usize,
        cols: usize,
        lx: f64,
        ly: f64,
        bending_stiffness: f64,
        membrane_stiffness: f64,
    ) -> Self {
        let dx = lx / (cols.saturating_sub(1).max(1)) as f64;
        let dy = ly / (rows.saturating_sub(1).max(1)) as f64;
        let mut positions = Vec::with_capacity(rows * cols);
        let frames = vec![DirectorTriad::canonical(); rows * cols];
        for j in 0..rows {
            for i in 0..cols {
                positions.push([i as f64 * dx, j as f64 * dy, 0.0]);
            }
        }
        Self {
            rows,
            cols,
            positions,
            frames,
            dx,
            dy,
            bending_stiffness,
            membrane_stiffness,
        }
    }

    /// Index into the flat node array.
    pub fn idx(&self, row: usize, col: usize) -> usize {
        row * self.cols + col
    }

    /// Estimate bending energy from out-of-plane curvature (finite difference).
    pub fn bending_energy(&self) -> f64 {
        let mut energy = 0.0;
        for j in 1..self.rows.saturating_sub(1) {
            for i in 1..self.cols.saturating_sub(1) {
                let p = self.positions[self.idx(j, i)];
                let px = self.positions[self.idx(j, i + 1)];
                let pm = self.positions[self.idx(j, i - 1)];
                let py = self.positions[self.idx(j + 1, i)];
                let pn = self.positions[self.idx(j - 1, i)];
                // Laplacian ∇²w ≈ (w_px + w_mx + w_py + w_ny - 4w_c) / h²
                let kappa_xx = (px[2] - 2.0 * p[2] + pm[2]) / (self.dx * self.dx);
                let kappa_yy = (py[2] - 2.0 * p[2] + pn[2]) / (self.dy * self.dy);
                let kappa_sum = kappa_xx + kappa_yy;
                energy += 0.5 * self.bending_stiffness * kappa_sum * kappa_sum * self.dx * self.dy;
            }
        }
        energy
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Wrinkle simulation
// ─────────────────────────────────────────────────────────────────────────────

/// Wrinkling criterion: Stein-Hedgepeth tension field theory.
///
/// A membrane wrinkles when the minor principal stress σ₂ < 0.
/// Wrinkling state is classified by the principal stresses.
///
/// Returns `WrinkleState`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WrinkleState {
    /// Both principal stresses are non-negative: taut membrane.
    Taut,
    /// Minor principal stress < 0: wrinkled membrane.
    Wrinkled,
    /// Both principal stresses are negative: slack membrane.
    Slack,
}

/// Classify wrinkling state from principal stresses.
///
/// # Arguments
/// * `sigma1` – Major principal stress \[Pa\].
/// * `sigma2` – Minor principal stress \[Pa\].
pub fn wrinkling_state(sigma1: f64, sigma2: f64) -> WrinkleState {
    if sigma1 < 0.0 && sigma2 < 0.0 {
        WrinkleState::Slack
    } else if sigma2 < 0.0 {
        WrinkleState::Wrinkled
    } else {
        WrinkleState::Taut
    }
}

/// Estimate wrinkle wavelength λ_w for a thin elastic membrane.
///
/// λ_w = π · (D / N_w)^(1/4)
///
/// where D = E·h³/(12(1−ν²)) is the bending stiffness and
/// N_w is the compressive membrane force per unit width.
///
/// # Arguments
/// * `bending_stiffness_d` – Bending stiffness D \[N·m\].
/// * `compressive_force_per_width` – |N_w| \[N/m\].
pub fn wrinkle_wavelength(bending_stiffness_d: f64, compressive_force_per_width: f64) -> f64 {
    if compressive_force_per_width < 1e-30 {
        f64::INFINITY
    } else {
        PI * (bending_stiffness_d / compressive_force_per_width).powf(0.25)
    }
}

/// Estimate wrinkle amplitude from stretch ratio and plate geometry.
///
/// A_w ≈ sqrt(ε · L²  / (π²·t)) where ε is excess arc-length strain.
///
/// # Arguments
/// * `excess_strain` – Excess compressive strain ε (>0).
/// * `wavelength` – Wrinkle wavelength λ_w \[m\].
pub fn wrinkle_amplitude(excess_strain: f64, wavelength: f64) -> f64 {
    if excess_strain <= 0.0 {
        0.0
    } else {
        wavelength * (excess_strain / PI).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Intrinsic curvature helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Constant curvature arc: position along arc parameterised by arc-length s.
///
/// For a planar arc in the XY plane with curvature κ:
/// r(s) = (sin(κs)/κ, (1-cos(κs))/κ, 0)
///
/// # Arguments
/// * `kappa` – Curvature κ \[rad/m\] (0 = straight).
/// * `s` – Arc-length parameter \[m\].
pub fn constant_curvature_position(kappa: f64, s: f64) -> [f64; 3] {
    if kappa.abs() < 1e-30 {
        [s, 0.0, 0.0]
    } else {
        [
            (kappa * s).sin() / kappa,
            (1.0 - (kappa * s).cos()) / kappa,
            0.0,
        ]
    }
}

/// Helix centreline: position along a helix with curvature κ and torsion τ.
///
/// r(s) = (r·cos(s/ρ), r·sin(s/ρ), p·s/(2π·ρ))
/// where ρ = 1/κ is the radius of curvature.
///
/// # Arguments
/// * `kappa` – Curvature κ \[rad/m\].
/// * `tau_helix` – Torsion τ \[rad/m\].
/// * `s` – Arc-length \[m\].
pub fn helix_position(kappa: f64, tau_helix: f64, s: f64) -> [f64; 3] {
    if kappa.abs() < 1e-30 {
        return [0.0, 0.0, s];
    }
    let _rho = 1.0 / kappa;
    // For a helix: r = kappa/(kappa²+tau²) (curvature sets radius)
    let denom = kappa * kappa + tau_helix * tau_helix;
    let r_helix = kappa / denom;
    let pitch_over_2pi = tau_helix / denom;
    let theta = s * (kappa * kappa + tau_helix * tau_helix).sqrt();
    [
        r_helix * theta.cos(),
        r_helix * theta.sin(),
        pitch_over_2pi * theta,
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Geometric exactness helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Rotation vector θ (axis-angle) from rotation matrix R.
///
/// φ = acos((tr(R) − 1)/2), axis = (1/(2 sin φ))·(R − Rᵀ) axial vector.
/// Returns \[axis_x·φ, axis_y·φ, axis_z·φ\].
pub fn rotation_vector_from_matrix(r: [f64; 9]) -> [f64; 3] {
    let trace = r[0] + r[4] + r[8];
    let cos_phi = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
    let phi = cos_phi.acos();
    if phi.abs() < 1e-12 {
        return [0.0; 3];
    }
    let s = 1.0 / (2.0 * phi.sin());
    let wx = (r[7] - r[5]) * s;
    let wy = (r[2] - r[6]) * s;
    let wz = (r[3] - r[1]) * s;
    [wx * phi, wy * phi, wz * phi]
}

/// Compose two rotations represented as rotation vectors (axis-angle).
///
/// Uses Rodrigues formula: first apply θ₁, then θ₂.
/// Returns the combined rotation vector.
pub fn compose_rotation_vectors(theta1: [f64; 3], theta2: [f64; 3]) -> [f64; 3] {
    let phi1 = norm3(theta1);
    let phi2 = norm3(theta2);
    let r1 = if phi1 < 1e-30 {
        mat3_identity()
    } else {
        rodrigues(theta1, phi1)
    };
    let r2 = if phi2 < 1e-30 {
        mat3_identity()
    } else {
        rodrigues(theta2, phi2)
    };
    let r_combined = mat3_mul(r2, r1);
    rotation_vector_from_matrix(r_combined)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. DirectorTriad canonical frame has orthonormal directors.
    #[test]
    fn test_canonical_frame_orthonormal() {
        let frame = DirectorTriad::canonical();
        assert!((dot3(frame.d1, frame.d2)).abs() < EPS, "d1 ⊥ d2");
        assert!((dot3(frame.d2, frame.d3)).abs() < EPS, "d2 ⊥ d3");
        assert!((dot3(frame.d1, frame.d3)).abs() < EPS, "d1 ⊥ d3");
        assert!((norm3(frame.d1) - 1.0).abs() < EPS, "d1 unit");
        assert!((norm3(frame.d3) - 1.0).abs() < EPS, "d3 unit");
    }

    // 2. Rodrigues rotation with angle 0 gives identity.
    #[test]
    fn test_rodrigues_zero_angle() {
        let r = rodrigues([0.0, 0.0, 1.0], 0.0);
        let v = [1.0, 2.0, 3.0];
        let rv = mat3_mul_vec(r, v);
        for i in 0..3 {
            assert!(
                (rv[i] - v[i]).abs() < EPS,
                "Rodrigues(θ=0) should give identity: {rv:?} vs {v:?}"
            );
        }
    }

    // 3. Rodrigues rotation 90° about z rotates x to y.
    #[test]
    fn test_rodrigues_90_about_z() {
        let r = rodrigues([0.0, 0.0, 1.0], PI / 2.0);
        let v = [1.0, 0.0, 0.0];
        let rv = mat3_mul_vec(r, v);
        assert!(
            (rv[0]).abs() < 1e-10,
            "x→y: x component should be 0, got {}",
            rv[0]
        );
        assert!(
            (rv[1] - 1.0).abs() < 1e-10,
            "x→y: y component should be 1, got {}",
            rv[1]
        );
        assert!(
            (rv[2]).abs() < 1e-10,
            "x→y: z component should be 0, got {}",
            rv[2]
        );
    }

    // 4. Darboux vector is zero for adjacent identical frames.
    #[test]
    fn test_darboux_zero_for_same_frames() {
        let frame = DirectorTriad::canonical();
        let kappa = darboux_vector(&frame, &frame, 0.1);
        for i in 0..3 {
            assert!(
                kappa[i].abs() < EPS,
                "Darboux vector should be zero for identical frames, got {kappa:?}"
            );
        }
    }

    // 5. Stretch-shear strain is zero for canonical frame with dr/ds = e3.
    #[test]
    fn test_stretch_shear_zero_for_straight() {
        let frame = DirectorTriad::canonical();
        let dr_ds = [0.0, 0.0, 1.0]; // d3 direction
        let gamma = stretch_shear_strain(&frame, dr_ds);
        for i in 0..3 {
            assert!(
                gamma[i].abs() < EPS,
                "No strain for undeformed straight rod: gamma={gamma:?}"
            );
        }
    }

    // 6. Kirchhoff energy is zero at rest (kappa = kappa0).
    #[test]
    fn test_kirchhoff_energy_zero_at_rest() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let kappa0 = [0.0, 0.0, 0.0];
        let energy = kirchhoff_energy_density(&props, kappa0, kappa0);
        assert!(
            energy.abs() < EPS,
            "Kirchhoff energy at rest should be zero, got {energy}"
        );
    }

    // 7. Kirchhoff energy increases with curvature departure from rest.
    #[test]
    fn test_kirchhoff_energy_increases_with_deformation() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let kappa0 = [0.0; 3];
        let kappa1 = [1.0, 0.0, 0.0];
        let kappa2 = [2.0, 0.0, 0.0];
        let e1 = kirchhoff_energy_density(&props, kappa1, kappa0);
        let e2 = kirchhoff_energy_density(&props, kappa2, kappa0);
        assert!(
            e2 > e1,
            "Higher curvature → higher elastic energy: {e1} vs {e2}"
        );
    }

    // 8. Internal moment is zero at rest.
    #[test]
    fn test_internal_moment_zero_at_rest() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let m = internal_moment(&props, [0.0; 3], [0.0; 3]);
        for i in 0..3 {
            assert!(m[i].abs() < EPS, "Moment should be zero at rest, got {m:?}");
        }
    }

    // 9. Internal force has correct sign for tension.
    #[test]
    fn test_internal_force_sign() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let gamma_tension = [0.0, 0.0, 0.1]; // positive stretch → tension
        let f = internal_force(&props, gamma_tension);
        assert!(
            f[2] > 0.0,
            "Positive stretch → positive axial force, got {}",
            f[2]
        );
    }

    // 10. Discrete rod straight creation has correct number of nodes.
    #[test]
    fn test_discrete_rod_creation() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let rod = DiscreteCosseratRod::straight(10, 1.0, props);
        assert_eq!(rod.nodes.len(), 11, "10 segments → 11 nodes");
        assert_eq!(rod.num_segments(), 10);
    }

    // 11. Discrete rod first node is pinned.
    #[test]
    fn test_discrete_rod_first_pinned() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let rod = DiscreteCosseratRod::straight(5, 0.5, props);
        assert!(rod.nodes[0].is_pinned, "First node should be pinned");
    }

    // 12. Euler buckling load decreases with length.
    #[test]
    fn test_buckling_load_decreases_with_length() {
        let ei = 1000.0;
        let p1 = euler_buckling_load(ei, 1.0, 1.0);
        let p2 = euler_buckling_load(ei, 2.0, 1.0);
        assert!(p2 < p1, "Longer rod buckles at lower load: {p1} vs {p2}");
    }

    // 13. Euler buckling load for fixed-free (K=2) is quarter of pinned-pinned (K=1).
    #[test]
    fn test_buckling_load_boundary_conditions() {
        let ei = 1000.0;
        let l = 1.0;
        let p_pinned = euler_buckling_load(ei, l, 1.0);
        let p_cantilever = euler_buckling_load(ei, l, 2.0);
        assert!(
            (p_pinned - 4.0 * p_cantilever).abs() < 1e-6,
            "Cantilever buckling = pinned/4: {p_pinned} vs {p_cantilever}"
        );
    }

    // 14. Tendon bending moment scales with eccentricity.
    #[test]
    fn test_tendon_moment_scales_with_eccentricity() {
        let m1 = tendon_bending_moment(100.0, 0.01);
        let m2 = tendon_bending_moment(100.0, 0.02);
        assert!(
            (m2 - 2.0 * m1).abs() < EPS,
            "Tendon moment scales with eccentricity"
        );
    }

    // 15. Tendon induced curvature is positive for positive tension.
    #[test]
    fn test_tendon_curvature_positive() {
        let kappa = tendon_induced_curvature(100.0, 0.01, 1000.0);
        assert!(
            kappa > 0.0,
            "Positive tension → positive curvature, got {kappa}"
        );
    }

    // 16. Wrinkling state: taut for both positive stresses.
    #[test]
    fn test_wrinkling_taut() {
        assert_eq!(wrinkling_state(100.0, 50.0), WrinkleState::Taut);
    }

    // 17. Wrinkling state: wrinkled for negative minor stress.
    #[test]
    fn test_wrinkling_wrinkled() {
        assert_eq!(wrinkling_state(100.0, -10.0), WrinkleState::Wrinkled);
    }

    // 18. Wrinkling state: slack for both negative stresses.
    #[test]
    fn test_wrinkling_slack() {
        assert_eq!(wrinkling_state(-10.0, -50.0), WrinkleState::Slack);
    }

    // 19. Wrinkle wavelength increases with bending stiffness.
    #[test]
    fn test_wrinkle_wavelength_increases_with_stiffness() {
        let lam1 = wrinkle_wavelength(1e-4, 1.0);
        let lam2 = wrinkle_wavelength(1e-3, 1.0);
        assert!(
            lam2 > lam1,
            "Higher D → longer wrinkle wavelength: {lam1} vs {lam2}"
        );
    }

    // 20. Wrinkle amplitude is zero for non-positive excess strain.
    #[test]
    fn test_wrinkle_amplitude_zero_no_strain() {
        let a = wrinkle_amplitude(0.0, 0.01);
        assert!(
            a.abs() < EPS,
            "Zero excess strain → zero wrinkle amplitude, got {a}"
        );
    }

    // 21. Constant curvature: at s=0 gives origin.
    #[test]
    fn test_constant_curvature_origin() {
        let pos = constant_curvature_position(1.0, 0.0);
        for i in 0..3 {
            assert!(
                pos[i].abs() < EPS,
                "Arc at s=0 should be at origin, got {pos:?}"
            );
        }
    }

    // 22. Constant curvature: zero curvature gives straight line.
    #[test]
    fn test_constant_curvature_straight() {
        let pos = constant_curvature_position(0.0, 2.5);
        assert!(
            (pos[0] - 2.5).abs() < EPS,
            "Zero curvature → straight line at x=s: {pos:?}"
        );
    }

    // 23. Cosserat plate flat creation has correct node count.
    #[test]
    fn test_cosserat_plate_node_count() {
        let plate = CosseratPlate::flat(4, 5, 1.0, 0.8, 1e-3, 1e5);
        assert_eq!(plate.positions.len(), 20, "4×5 plate has 20 nodes");
        assert_eq!(plate.rows, 4);
        assert_eq!(plate.cols, 5);
    }

    // 24. Plate bending energy is zero for flat plate.
    #[test]
    fn test_plate_bending_energy_flat() {
        let plate = CosseratPlate::flat(4, 4, 1.0, 1.0, 1e-3, 1e5);
        let e = plate.bending_energy();
        assert!(
            e.abs() < EPS,
            "Flat plate bending energy should be zero, got {e}"
        );
    }

    // 25. Active rod activation converges toward target curvature.
    #[test]
    fn test_active_rod_activation_converges() {
        let props = CosseratRodProperties::circular(0.005, 100e9, 40e9, 1000.0);
        let rod = DiscreteCosseratRod::straight(5, 0.5, props);
        let mut active = ActiveCosseratRod::new(rod, 0.1);
        let target = [1.0, 0.0, 0.0];
        active.set_target_curvature(2, target);
        for _ in 0..100 {
            active.update_activation(0.01);
        }
        let kappa = active.rod.intrinsic_kappa[2];
        assert!(
            (kappa[0] - 1.0).abs() < 0.01,
            "Active curvature should converge to target: {kappa:?}"
        );
    }

    // 26. Rotation vector from identity matrix is zero vector.
    #[test]
    fn test_rotation_vector_from_identity() {
        let id = mat3_identity();
        let rv = rotation_vector_from_matrix(id);
        for i in 0..3 {
            assert!(
                rv[i].abs() < EPS,
                "Identity rotation vector should be zero, got {rv:?}"
            );
        }
    }

    // 27. Compose rotations: 90°+90° about z = 180° about z.
    #[test]
    fn test_compose_rotations_90_90() {
        let theta90 = [0.0, 0.0, PI / 2.0];
        let composed = compose_rotation_vectors(theta90, theta90);
        // Should be ≈ [0, 0, π].
        assert!(
            (composed[2].abs() - PI).abs() < 1e-10,
            "90°+90° about z should give 180°, got {composed:?}"
        );
    }

    // 28. CosseratRodProperties circular: EI1 = EI2 (isotropic circular cross-section).
    #[test]
    fn test_circular_rod_isotropic() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        assert!(
            (props.ei1 - props.ei2).abs() < EPS,
            "Circular rod should have EI1=EI2: {} vs {}",
            props.ei1,
            props.ei2
        );
    }

    // 29. Discrete rod elastic energy is zero when straight (no deformation).
    #[test]
    fn test_elastic_energy_zero_straight() {
        let props = CosseratRodProperties::circular(0.01, 200e9, 80e9, 7800.0);
        let rod = DiscreteCosseratRod::straight(5, 1.0, props);
        let energy = rod.elastic_energy();
        assert!(
            energy.abs() < EPS,
            "Straight rod at rest should have zero elastic energy, got {energy}"
        );
    }

    // 30. DiscreteCosseratRod step: free end node moves under gravity.
    #[test]
    fn test_rod_step_gravity() {
        let props = CosseratRodProperties::circular(0.005, 1e9, 0.4e9, 1000.0);
        let mut rod = DiscreteCosseratRod::straight(4, 0.4, props);
        let initial_pos = rod.nodes[4].position;
        let g = [0.0, -9.81, 0.0];
        let dt = 1e-3;
        for _ in 0..50 {
            rod.apply_gravity(g);
            rod.step(dt, 0.01);
        }
        let final_y = rod.nodes[4].position[1];
        assert!(
            final_y < initial_pos[1],
            "Free end should fall under gravity: initial_y={}, final_y={final_y}",
            initial_pos[1]
        );
    }
}
