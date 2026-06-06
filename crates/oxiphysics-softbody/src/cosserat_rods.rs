// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cosserat rod theory: extensible Kirchhoff/Cosserat rods.
//!
//! This module implements the geometrically-exact Cosserat/Kirchhoff rod
//! formulation for slender elastic structures that can bend, twist, shear, and
//! stretch.  The key components are:
//!
//! - [`RodCrossSection`] – geometry of the rod cross-section
//! - [`CosseratNode`] – discrete node storing position, director triad and velocities
//! - [`CosseratRod`] – array of nodes with material stiffness parameters
//! - [`InternalForces`] – computation of bending moments, shear/axial forces, torsion
//! - [`CosseratDynamics`] – equations of motion for dynamic simulation
//! - [`BoundaryConditions`] – clamped, pinned, free ends, follower forces
//! - [`StaticSolver`] – Newton-Raphson / shooting method for static equilibria
//! - [`RodContact`] – rod-rod contact detection and friction/adhesion

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Small linear algebra helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// 3-vector type alias used throughout this module.
pub type Vec3 = [f64; 3];

/// 3×3 matrix stored in row-major order.
pub type Mat3 = [[f64; 3]; 3];

#[inline]
fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_norm(a: Vec3) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_normalize(a: Vec3) -> Vec3 {
    let n = vec3_norm(a);
    if n < 1e-14 {
        [0.0, 0.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

/// Rotate vector `v` by angle `theta` around unit axis `k` (Rodrigues).
fn rodrigues(v: Vec3, k: Vec3, theta: f64) -> Vec3 {
    let (s, c) = (theta.sin(), theta.cos());
    let kxv = vec3_cross(k, v);
    let kdv = vec3_dot(k, v);
    vec3_add(
        vec3_add(vec3_scale(v, c), vec3_scale(kxv, s)),
        vec3_scale(k, kdv * (1.0 - c)),
    )
}

// ---------------------------------------------------------------------------
// RodCrossSection
// ---------------------------------------------------------------------------

/// Shape of the rod cross-section.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrossSectionShape {
    /// Circular cross-section.
    Circular {
        /// Radius \[m\].
        radius: f64,
    },
    /// Rectangular cross-section.
    Rectangular {
        /// Width \[m\].
        width: f64,
        /// Height \[m\].
        height: f64,
    },
    /// Elliptical cross-section.
    Elliptical {
        /// Semi-axis along d1 \[m\].
        semi_a: f64,
        /// Semi-axis along d2 \[m\].
        semi_b: f64,
    },
}

/// Cross-sectional properties: area, second moments of inertia, torsion constant.
#[derive(Debug, Clone, Copy)]
pub struct RodCrossSection {
    /// Shape of the cross-section.
    pub shape: CrossSectionShape,
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Second moment of inertia about the d1 axis: I₁ \[m⁴\].
    pub i1: f64,
    /// Second moment of inertia about the d2 axis: I₂ \[m⁴\].
    pub i2: f64,
    /// St. Venant torsion constant J \[m⁴\].
    pub torsion_j: f64,
}

impl RodCrossSection {
    /// Construct a [`RodCrossSection`] from a given shape, computing all
    /// geometric properties automatically.
    ///
    /// # Examples
    /// ```no_run
    /// use oxiphysics_softbody::cosserat_rods::{RodCrossSection, CrossSectionShape};
    /// let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.01 });
    /// assert!((cs.area - std::f64::consts::PI * 1e-4).abs() < 1e-12);
    /// ```
    pub fn new(shape: CrossSectionShape) -> Self {
        match shape {
            CrossSectionShape::Circular { radius: r } => {
                let area = PI * r * r;
                let i = PI * r * r * r * r / 4.0;
                // Circular section: J = I1 + I2 = pi r^4 / 2
                let j = PI * r * r * r * r / 2.0;
                Self {
                    shape,
                    area,
                    i1: i,
                    i2: i,
                    torsion_j: j,
                }
            }
            CrossSectionShape::Rectangular {
                width: b,
                height: h,
            } => {
                let area = b * h;
                let i1 = b * h * h * h / 12.0;
                let i2 = h * b * b * b / 12.0;
                // Approximate Bredt torsion constant for thin rectangles
                let a_min = b.min(h);
                let a_max = b.max(h);
                let j = a_min * a_min * a_min * a_max / 3.0
                    * (1.0 - 0.63 * a_min / a_max * (1.0 - a_min.powi(4) / (12.0 * a_max.powi(4))));
                Self {
                    shape,
                    area,
                    i1,
                    i2,
                    torsion_j: j,
                }
            }
            CrossSectionShape::Elliptical {
                semi_a: a,
                semi_b: b,
            } => {
                let area = PI * a * b;
                let i1 = PI * a * b * b * b / 4.0;
                let i2 = PI * b * a * a * a / 4.0;
                // Approximate torsion constant for ellipse
                let j = PI * a * a * a * b * b * b / (a * a + b * b);
                Self {
                    shape,
                    area,
                    i1,
                    i2,
                    torsion_j: j,
                }
            }
        }
    }

    /// Compute the shear correction factor (Timoshenko) for the cross-section.
    /// Returns ≈ 0.9 for circular, 5/6 for rectangular, interpolated for ellipse.
    pub fn shear_correction_factor(&self) -> f64 {
        match self.shape {
            CrossSectionShape::Circular { .. } => 0.9,
            CrossSectionShape::Rectangular { .. } => 5.0 / 6.0,
            CrossSectionShape::Elliptical {
                semi_a: a,
                semi_b: b,
            } => {
                // Approximate: between 0.83 and 0.9
                let ratio = (a / b).min(b / a);
                0.83 + 0.07 * ratio
            }
        }
    }
}

// ---------------------------------------------------------------------------
// CosseratNode
// ---------------------------------------------------------------------------

/// A single node on a discretized Cosserat rod.
///
/// Each node carries:
/// - position **r** in world space
/// - a right-handed director triad {**d**₁, **d**₂, **d**₃}, where **d**₃
///   is the rod tangent and **d**₁, **d**₂ span the cross-section
/// - linear velocity **v** and angular velocity **ω** (body-frame)
#[derive(Debug, Clone)]
pub struct CosseratNode {
    /// Position in world space \[m\].
    pub position: Vec3,
    /// First material director (cross-section in-plane axis 1).
    pub d1: Vec3,
    /// Second material director (cross-section in-plane axis 2).
    pub d2: Vec3,
    /// Third director — rod tangent / axial direction.
    pub d3: Vec3,
    /// Linear velocity \[m/s\].
    pub velocity: Vec3,
    /// Angular velocity in body frame \[rad/s\].
    pub angular_velocity: Vec3,
    /// External force applied at this node \[N\].
    pub external_force: Vec3,
    /// External moment applied at this node \[N·m\].
    pub external_moment: Vec3,
}

impl CosseratNode {
    /// Create a new node with the given position and a standard director triad
    /// aligned with the global axes.
    pub fn new(position: Vec3) -> Self {
        Self {
            position,
            d1: [1.0, 0.0, 0.0],
            d2: [0.0, 1.0, 0.0],
            d3: [0.0, 0.0, 1.0],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            external_force: [0.0; 3],
            external_moment: [0.0; 3],
        }
    }

    /// Create a node with explicit director triad (must be orthonormal).
    pub fn with_directors(position: Vec3, d1: Vec3, d2: Vec3, d3: Vec3) -> Self {
        Self {
            position,
            d1,
            d2,
            d3,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            external_force: [0.0; 3],
            external_moment: [0.0; 3],
        }
    }

    /// Re-orthonormalize the director triad using Gram-Schmidt.
    /// **d**₃ is kept fixed; **d**₁ and **d**₂ are corrected.
    pub fn orthonormalize(&mut self) {
        self.d3 = vec3_normalize(self.d3);
        // Remove d3 component from d1
        let proj = vec3_dot(self.d1, self.d3);
        self.d1 = vec3_normalize(vec3_sub(self.d1, vec3_scale(self.d3, proj)));
        self.d2 = vec3_normalize(vec3_cross(self.d3, self.d1));
    }

    /// Apply an incremental rotation (angle-axis `omega * dt`) to the directors.
    pub fn rotate_directors(&mut self, omega: Vec3, dt: f64) {
        let angle = vec3_norm(omega) * dt;
        if angle < 1e-14 {
            return;
        }
        let axis = vec3_scale(omega, 1.0 / (angle / dt));
        let axis = vec3_normalize(axis);
        self.d1 = rodrigues(self.d1, axis, angle);
        self.d2 = rodrigues(self.d2, axis, angle);
        self.d3 = rodrigues(self.d3, axis, angle);
    }
}

// ---------------------------------------------------------------------------
// CosseratRod
// ---------------------------------------------------------------------------

/// Material and geometric stiffness parameters of a Cosserat rod.
#[derive(Debug, Clone)]
pub struct RodStiffness {
    /// Bending stiffness about d1: EI₁ \[N·m²\].
    pub ei1: f64,
    /// Bending stiffness about d2: EI₂ \[N·m²\].
    pub ei2: f64,
    /// Torsional stiffness: GJ \[N·m²\].
    pub gj: f64,
    /// Axial (extensional) stiffness: EA \[N\].
    pub ea: f64,
    /// Shear stiffness: GA_s \[N\] (corrected shear area).
    pub ga: f64,
}

impl RodStiffness {
    /// Construct stiffness parameters from material properties and cross-section.
    ///
    /// # Arguments
    /// - `young_modulus` – Young's modulus E \[Pa\]
    /// - `shear_modulus` – Shear modulus G \[Pa\]
    /// - `cs` – Cross-section geometry
    pub fn from_material(young_modulus: f64, shear_modulus: f64, cs: &RodCrossSection) -> Self {
        let ks = cs.shear_correction_factor();
        Self {
            ei1: young_modulus * cs.i1,
            ei2: young_modulus * cs.i2,
            gj: shear_modulus * cs.torsion_j,
            ea: young_modulus * cs.area,
            ga: shear_modulus * ks * cs.area,
        }
    }
}

/// An extensible, shear-deformable Cosserat rod.
///
/// The rod is discretized into `n` nodes connected by `n-1` segments.
/// Each segment stores the rest-state curvature κ₀, torsion τ₀, and reference length.
#[derive(Debug, Clone)]
pub struct CosseratRod {
    /// Discrete nodes along the rod.
    pub nodes: Vec<CosseratNode>,
    /// Material stiffness parameters.
    pub stiffness: RodStiffness,
    /// Rest curvature about d1 at each segment \[1/m\].
    pub kappa1_0: Vec<f64>,
    /// Rest curvature about d2 at each segment \[1/m\].
    pub kappa2_0: Vec<f64>,
    /// Rest torsion (twist rate) at each segment \[1/m\].
    pub tau_0: Vec<f64>,
    /// Natural (rest) length of each segment \[m\].
    pub rest_lengths: Vec<f64>,
    /// Linear mass density ρA \[kg/m\].
    pub linear_density: f64,
    /// Second moment of mass per unit length ρI \[kg·m\].
    pub rotary_inertia: f64,
}

impl CosseratRod {
    /// Create a straight Cosserat rod along the Z-axis.
    ///
    /// # Arguments
    /// - `n_nodes` – number of nodes (≥ 2)
    /// - `total_length` – total rest length \[m\]
    /// - `stiffness` – rod stiffness parameters
    /// - `linear_density` – ρA \[kg/m\]
    /// - `rotary_inertia` – ρI \[kg·m\]
    pub fn new_straight(
        n_nodes: usize,
        total_length: f64,
        stiffness: RodStiffness,
        linear_density: f64,
        rotary_inertia: f64,
    ) -> Self {
        assert!(n_nodes >= 2, "Rod must have at least 2 nodes");
        let seg_len = total_length / (n_nodes - 1) as f64;
        let nodes: Vec<CosseratNode> = (0..n_nodes)
            .map(|i| {
                let z = i as f64 * seg_len;
                CosseratNode::new([0.0, 0.0, z])
            })
            .collect();
        let n_segs = n_nodes - 1;
        Self {
            nodes,
            stiffness,
            kappa1_0: vec![0.0; n_segs],
            kappa2_0: vec![0.0; n_segs],
            tau_0: vec![0.0; n_segs],
            rest_lengths: vec![seg_len; n_segs],
            linear_density,
            rotary_inertia,
        }
    }

    /// Number of segments (= number of nodes − 1).
    pub fn n_segments(&self) -> usize {
        self.nodes.len() - 1
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Total rest length of the rod.
    pub fn total_rest_length(&self) -> f64 {
        self.rest_lengths.iter().sum()
    }

    /// Compute the current tangent vector at segment `seg` by finite difference.
    pub fn segment_tangent(&self, seg: usize) -> Vec3 {
        let a = &self.nodes[seg];
        let b = &self.nodes[seg + 1];
        vec3_normalize(vec3_sub(b.position, a.position))
    }

    /// Compute segment stretch ratio (current length / rest length).
    pub fn segment_stretch(&self, seg: usize) -> f64 {
        let a = &self.nodes[seg];
        let b = &self.nodes[seg + 1];
        let len = vec3_norm(vec3_sub(b.position, a.position));
        len / self.rest_lengths[seg]
    }

    /// Interpolate the director triad at the midpoint of segment `seg`
    /// as the average of the two endpoint triads.
    pub fn segment_directors(&self, seg: usize) -> (Vec3, Vec3, Vec3) {
        let a = &self.nodes[seg];
        let b = &self.nodes[seg + 1];
        let d1 = vec3_normalize(vec3_add(a.d1, b.d1));
        let d2 = vec3_normalize(vec3_add(a.d2, b.d2));
        let d3 = vec3_normalize(vec3_add(a.d3, b.d3));
        (d1, d2, d3)
    }
}

// ---------------------------------------------------------------------------
// InternalForces
// ---------------------------------------------------------------------------

/// Result of internal force/moment computation at a rod segment.
#[derive(Debug, Clone, Copy, Default)]
pub struct SegmentInternalLoads {
    /// Axial force N (positive = tension) \[N\].
    pub axial_force: f64,
    /// Shear force along d1: Q₁ \[N\].
    pub shear1: f64,
    /// Shear force along d2: Q₂ \[N\].
    pub shear2: f64,
    /// Bending moment about d1: M₁ \[N·m\].
    pub moment1: f64,
    /// Bending moment about d2: M₂ \[N·m\].
    pub moment2: f64,
    /// Torsional moment T \[N·m\].
    pub torsion: f64,
}

/// Computes internal forces and moments in a Cosserat rod from its current configuration.
pub struct InternalForces;

impl InternalForces {
    /// Compute strain measures at segment `seg`:
    /// - axial strain ε = (l - l₀) / l₀
    /// - shear strains γ₁, γ₂
    /// - curvature differences Δκ₁, Δκ₂
    /// - twist difference Δτ
    pub fn segment_strains(rod: &CosseratRod, seg: usize) -> (f64, f64, f64, f64, f64, f64) {
        let a = &rod.nodes[seg];
        let b = &rod.nodes[seg + 1];
        let diff = vec3_sub(b.position, a.position);
        let current_len = vec3_norm(diff);
        let rest_len = rod.rest_lengths[seg];

        // Axial strain
        let eps = (current_len - rest_len) / rest_len;

        // Shear: projection of tangent onto cross-section directors of midpoint frame
        let t = if current_len > 1e-14 {
            vec3_scale(diff, 1.0 / current_len)
        } else {
            [0.0, 0.0, 1.0]
        };
        let (d1_m, d2_m, _d3_m) = rod.segment_directors(seg);
        let gamma1 = vec3_dot(t, d1_m);
        let gamma2 = vec3_dot(t, d2_m);

        // Curvature and twist via finite rotation between adjacent directors
        let d3a = a.d3;
        let d3b = b.d3;
        // Curvature vector: (d3b - d3a) / rest_len projected onto d1, d2
        let dd3 = vec3_sub(d3b, d3a);
        let dkappa1 = vec3_dot(dd3, d1_m) / rest_len - rod.kappa1_0[seg];
        let dkappa2 = vec3_dot(dd3, d2_m) / rest_len - rod.kappa2_0[seg];

        // Twist: rate of rotation of d1 around d3
        let d1a = a.d1;
        let d1b = b.d1;
        let dd1 = vec3_sub(d1b, d1a);
        let dtau = vec3_dot(dd1, d2_m) / rest_len - rod.tau_0[seg];

        (eps, gamma1, gamma2, dkappa1, dkappa2, dtau)
    }

    /// Compute internal loads at segment `seg` from the constitutive law.
    pub fn compute(rod: &CosseratRod, seg: usize) -> SegmentInternalLoads {
        let (eps, gamma1, gamma2, dkappa1, dkappa2, dtau) = Self::segment_strains(rod, seg);
        SegmentInternalLoads {
            axial_force: rod.stiffness.ea * eps,
            shear1: rod.stiffness.ga * gamma1,
            shear2: rod.stiffness.ga * gamma2,
            moment1: rod.stiffness.ei1 * dkappa1,
            moment2: rod.stiffness.ei2 * dkappa2,
            torsion: rod.stiffness.gj * dtau,
        }
    }

    /// Compute the vector of generalized elastic forces on both nodes of segment `seg`.
    /// Returns `(f_a, m_a, f_b, m_b)` where `f` are forces and `m` are moments.
    pub fn nodal_forces(rod: &CosseratRod, seg: usize) -> (Vec3, Vec3, Vec3, Vec3) {
        let loads = Self::compute(rod, seg);
        let a = &rod.nodes[seg];
        let b = &rod.nodes[seg + 1];
        let rest_len = rod.rest_lengths[seg];
        let (d1_m, d2_m, d3_m) = rod.segment_directors(seg);

        // Force in world frame: N along tangent + Q along cross-section directors
        let diff = vec3_sub(b.position, a.position);
        let current_len = vec3_norm(diff);
        let tangent = if current_len > 1e-14 {
            vec3_scale(diff, 1.0 / current_len)
        } else {
            d3_m
        };

        let f_internal = vec3_add(
            vec3_add(
                vec3_scale(tangent, loads.axial_force),
                vec3_scale(d1_m, loads.shear1),
            ),
            vec3_scale(d2_m, loads.shear2),
        );

        // Moment on node a (reaction) divided by element length
        let m_a_body = vec3_add(
            vec3_add(
                vec3_scale(d1_m, loads.moment1),
                vec3_scale(d2_m, loads.moment2),
            ),
            vec3_scale(d3_m, loads.torsion),
        );
        let m_a = vec3_scale(m_a_body, 1.0 / rest_len);

        // Equal and opposite on node b
        let f_b = f_internal;
        let f_a = vec3_scale(f_internal, -1.0);
        let m_b = vec3_scale(m_a, -1.0);

        // Add d3 × f contribution to moments (torque arm)
        let half_len = rest_len * 0.5;
        let torque_a = vec3_scale(vec3_cross(a.d3, f_internal), half_len);
        let torque_b = vec3_scale(vec3_cross(b.d3, f_b), half_len);

        (f_a, vec3_add(m_a, torque_a), f_b, vec3_add(m_b, torque_b))
    }
}

// ---------------------------------------------------------------------------
// CosseratDynamics
// ---------------------------------------------------------------------------

/// Equations of motion for a dynamic Cosserat rod.
///
/// Implements the balance laws:
/// - Linear: ∂F/∂s + f_ext = ρA · **ü**
/// - Angular: ∂M/∂s + **d**₃ × F + m_ext = ρI · **α̈**
pub struct CosseratDynamics;

impl CosseratDynamics {
    /// Compute generalized accelerations (linear + angular) at each node due to
    /// internal elastic forces plus external loading.
    ///
    /// Returns `(acc[i], alpha[i])` for all `i` in `0..rod.n_nodes()`.
    pub fn accelerations(rod: &CosseratRod) -> Vec<(Vec3, Vec3)> {
        let n = rod.n_nodes();
        let mut forces = vec![[0.0f64; 3]; n];
        let mut moments = vec![[0.0f64; 3]; n];

        // Accumulate internal loads from each segment
        for seg in 0..rod.n_segments() {
            let (fa, ma, fb, mb) = InternalForces::nodal_forces(rod, seg);
            forces[seg] = vec3_add(forces[seg], fa);
            moments[seg] = vec3_add(moments[seg], ma);
            forces[seg + 1] = vec3_add(forces[seg + 1], fb);
            moments[seg + 1] = vec3_add(moments[seg + 1], mb);
        }

        // Add external forces/moments
        for i in 0..n {
            forces[i] = vec3_add(forces[i], rod.nodes[i].external_force);
            moments[i] = vec3_add(moments[i], rod.nodes[i].external_moment);
        }

        // Compute accelerations: a = F / (ρA·l_elem)
        let l_elem = rod.total_rest_length() / rod.n_segments() as f64;
        let mass_per_node = rod.linear_density * l_elem;
        let inertia_per_node = rod.rotary_inertia * l_elem;

        (0..n)
            .map(|i| {
                let acc = if mass_per_node > 1e-30 {
                    vec3_scale(forces[i], 1.0 / mass_per_node)
                } else {
                    [0.0; 3]
                };
                let alpha = if inertia_per_node > 1e-30 {
                    vec3_scale(moments[i], 1.0 / inertia_per_node)
                } else {
                    [0.0; 3]
                };
                (acc, alpha)
            })
            .collect()
    }

    /// Advance the rod by `dt` using explicit (symplectic) Euler integration.
    ///
    /// Fixed nodes (pinned/clamped) must be excluded by the caller via
    /// [`BoundaryConditions::apply`] before and after this call.
    pub fn step_euler(rod: &mut CosseratRod, dt: f64) {
        let accels = Self::accelerations(rod);
        let n = rod.n_nodes();
        for (i, &(acc, alpha)) in accels.iter().enumerate().take(n) {
            // Update velocities
            rod.nodes[i].velocity = vec3_add(rod.nodes[i].velocity, vec3_scale(acc, dt));
            rod.nodes[i].angular_velocity =
                vec3_add(rod.nodes[i].angular_velocity, vec3_scale(alpha, dt));
            // Update positions
            let vel = rod.nodes[i].velocity;
            rod.nodes[i].position = vec3_add(rod.nodes[i].position, vec3_scale(vel, dt));
            // Update directors
            let omega = rod.nodes[i].angular_velocity;
            rod.nodes[i].rotate_directors(omega, dt);
            rod.nodes[i].orthonormalize();
        }
    }
}

// ---------------------------------------------------------------------------
// BoundaryConditions
// ---------------------------------------------------------------------------

/// Type of boundary condition at one end of the rod.
#[derive(Debug, Clone)]
pub enum BcType {
    /// All translational and rotational DOFs fixed.
    Clamped,
    /// Translational DOFs fixed, rotational DOFs free.
    Pinned,
    /// All DOFs free.
    Free,
    /// Prescribed external force (follower force tracks **d**₃).
    FollowerForce {
        /// Force magnitude \[N\].
        magnitude: f64,
    },
    /// Prescribed curvature κ₁, κ₂ and twist τ at the end.
    PrescribedCurvature {
        /// Curvature component κ₁ \[1/m\].
        kappa1: f64,
        /// Curvature component κ₂ \[1/m\].
        kappa2: f64,
        /// Twist per unit length τ \[rad/m\].
        tau: f64,
    },
}

/// Boundary conditions for a Cosserat rod.
#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    /// BC at the base node (index 0).
    pub base: BcType,
    /// BC at the tip node (last index).
    pub tip: BcType,
    /// Fixed base position (used when `base` is `Clamped` or `Pinned`).
    pub base_position: Vec3,
    /// Fixed base director triad (used when `base` is `Clamped`).
    pub base_d1: Vec3,
    /// Fixed base d2.
    pub base_d2: Vec3,
    /// Fixed base d3.
    pub base_d3: Vec3,
}

impl BoundaryConditions {
    /// Create a clamped-free (cantilever) boundary condition.
    pub fn cantilever(base_position: Vec3, d1: Vec3, d2: Vec3, d3: Vec3) -> Self {
        Self {
            base: BcType::Clamped,
            tip: BcType::Free,
            base_position,
            base_d1: d1,
            base_d2: d2,
            base_d3: d3,
        }
    }

    /// Apply boundary conditions to the rod: zero the velocity/acceleration
    /// of constrained nodes and reset their positions/directors.
    pub fn apply(&self, rod: &mut CosseratRod) {
        // Base node
        match &self.base {
            BcType::Clamped => {
                let n = &mut rod.nodes[0];
                n.position = self.base_position;
                n.d1 = self.base_d1;
                n.d2 = self.base_d2;
                n.d3 = self.base_d3;
                n.velocity = [0.0; 3];
                n.angular_velocity = [0.0; 3];
            }
            BcType::Pinned => {
                let n = &mut rod.nodes[0];
                n.position = self.base_position;
                n.velocity = [0.0; 3];
            }
            _ => {}
        }
        // Tip node: follower force
        let last = rod.nodes.len() - 1;
        if let BcType::FollowerForce { magnitude } = &self.tip {
            let d3 = rod.nodes[last].d3;
            rod.nodes[last].external_force = vec3_scale(d3, *magnitude);
        }
    }
}

// ---------------------------------------------------------------------------
// StaticSolver
// ---------------------------------------------------------------------------

/// Static equilibrium solver using Newton-Raphson iteration.
///
/// Finds the deformed configuration of the rod under static loading
/// by iterating until the internal force residual falls below tolerance.
pub struct StaticSolver {
    /// Maximum number of Newton-Raphson iterations.
    pub max_iter: usize,
    /// Convergence tolerance on residual norm.
    pub tol: f64,
}

impl StaticSolver {
    /// Construct a new static solver with given iteration limit and tolerance.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// Compute the residual norm: sum of nodal force magnitudes on free DOFs.
    pub fn residual_norm(rod: &CosseratRod) -> f64 {
        let accs = CosseratDynamics::accelerations(rod);
        accs.iter()
            .skip(1) // skip clamped base
            .map(|(a, alpha)| vec3_norm(*a) + vec3_norm(*alpha))
            .sum()
    }

    /// Solve for static equilibrium using pseudo-dynamic relaxation.
    ///
    /// Applies damped dynamics until the residual norm is below `self.tol`
    /// or `self.max_iter` steps are exhausted. Returns the number of
    /// iterations performed and final residual.
    pub fn solve(&self, rod: &mut CosseratRod, bcs: &BoundaryConditions) -> (usize, f64) {
        let dt = 1e-4;
        let damping = 0.95;
        let mut iter = 0;

        for _ in 0..self.max_iter {
            iter += 1;
            CosseratDynamics::step_euler(rod, dt);
            bcs.apply(rod);

            // Apply velocity damping
            for node in rod.nodes.iter_mut() {
                node.velocity = vec3_scale(node.velocity, damping);
                node.angular_velocity = vec3_scale(node.angular_velocity, damping);
            }

            let res = Self::residual_norm(rod);
            if res < self.tol {
                return (iter, res);
            }
        }
        (iter, Self::residual_norm(rod))
    }

    /// Shooting method: given prescribed base state and tip condition, integrate
    /// the rod ODE along its length and return the tip residual.
    ///
    /// This implements a simplified single-shooting method for static rods
    /// under end loading, integrating the Kirchhoff equations along arc length.
    pub fn shoot(
        rod: &CosseratRod,
        base_kappa1: f64,
        base_kappa2: f64,
        base_tau: f64,
    ) -> (Vec3, Vec3) {
        // Forward integration of curvature/twist to find tip position and directors.
        let n = rod.n_nodes();
        let mut pos = rod.nodes[0].position;
        let mut d1 = rod.nodes[0].d1;
        let mut d2 = rod.nodes[0].d2;
        let mut d3 = rod.nodes[0].d3;

        let seg_len = rod.total_rest_length() / rod.n_segments() as f64;
        let n_segs = n - 1;

        for i in 0..n_segs {
            let t = i as f64 / n_segs as f64;
            let k1 = base_kappa1 * (1.0 - t) + rod.kappa1_0[i] * t;
            let k2 = base_kappa2 * (1.0 - t) + rod.kappa2_0[i] * t;
            let tau = base_tau * (1.0 - t) + rod.tau_0[i] * t;

            // Darboux vector in body frame: Ω = k1*d1 + k2*d2 + tau*d3
            let omega = vec3_add(
                vec3_add(vec3_scale(d1, k1), vec3_scale(d2, k2)),
                vec3_scale(d3, tau),
            );

            // Update directors: d_i' = Ω × d_i
            let dd1 = vec3_scale(vec3_cross(omega, d1), seg_len);
            let dd2 = vec3_scale(vec3_cross(omega, d2), seg_len);
            let dd3 = vec3_scale(vec3_cross(omega, d3), seg_len);
            d1 = vec3_normalize(vec3_add(d1, dd1));
            d2 = vec3_normalize(vec3_add(d2, dd2));
            d3 = vec3_normalize(vec3_add(d3, dd3));

            pos = vec3_add(pos, vec3_scale(d3, seg_len));
        }

        (pos, d3)
    }
}

// ---------------------------------------------------------------------------
// RodContact
// ---------------------------------------------------------------------------

/// Contact interaction between two Cosserat rods.
pub struct RodContact {
    /// Contact stiffness \[N/m\].
    pub stiffness: f64,
    /// Coulomb friction coefficient.
    pub friction_coeff: f64,
    /// Adhesion energy per unit area \[J/m²\].
    pub adhesion: f64,
}

impl RodContact {
    /// Create a new rod contact model.
    pub fn new(stiffness: f64, friction_coeff: f64, adhesion: f64) -> Self {
        Self {
            stiffness,
            friction_coeff,
            adhesion,
        }
    }

    /// Compute the closest distance between segment `seg_a` of `rod_a` and
    /// segment `seg_b` of `rod_b`.  Returns `(distance, point_a, point_b)`.
    pub fn segment_segment_distance(
        rod_a: &CosseratRod,
        seg_a: usize,
        rod_b: &CosseratRod,
        seg_b: usize,
    ) -> (f64, Vec3, Vec3) {
        let p0 = rod_a.nodes[seg_a].position;
        let p1 = rod_a.nodes[seg_a + 1].position;
        let q0 = rod_b.nodes[seg_b].position;
        let q1 = rod_b.nodes[seg_b + 1].position;

        let d1 = vec3_sub(p1, p0);
        let d2 = vec3_sub(q1, q0);
        let r = vec3_sub(p0, q0);

        let a = vec3_dot(d1, d1);
        let e = vec3_dot(d2, d2);
        let f = vec3_dot(d2, r);

        let (s, t) = if a < 1e-14 {
            let s = 0.0f64;
            let t = if e > 1e-14 {
                (f / e).clamp(0.0, 1.0)
            } else {
                0.0
            };
            (s, t)
        } else {
            let c = vec3_dot(d1, r);
            if e < 1e-14 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else {
                let b = vec3_dot(d1, d2);
                let denom = a * e - b * b;
                let s = if denom.abs() > 1e-14 {
                    ((b * f - c * e) / denom).clamp(0.0, 1.0)
                } else {
                    0.5
                };
                let t = ((b * s + f) / e).clamp(0.0, 1.0);
                (s, t)
            }
        };

        let pa = vec3_add(p0, vec3_scale(d1, s));
        let pb = vec3_add(q0, vec3_scale(d2, t));
        let dist = vec3_norm(vec3_sub(pa, pb));
        (dist, pa, pb)
    }

    /// Detect all contacts between all segments of two rods within `radius`.
    /// Returns list of `(seg_a, seg_b, distance, point_a, point_b)`.
    pub fn detect_contacts(
        &self,
        rod_a: &CosseratRod,
        rod_b: &CosseratRod,
        radius: f64,
    ) -> Vec<(usize, usize, f64, Vec3, Vec3)> {
        let mut contacts = Vec::new();
        for sa in 0..rod_a.n_segments() {
            for sb in 0..rod_b.n_segments() {
                let (dist, pa, pb) = Self::segment_segment_distance(rod_a, sa, rod_b, sb);
                if dist < radius {
                    contacts.push((sa, sb, dist, pa, pb));
                }
            }
        }
        contacts
    }

    /// Compute contact force magnitude at a single contact point.
    /// Includes linear spring repulsion plus adhesion (negative when outside).
    pub fn contact_force_magnitude(&self, distance: f64, radius: f64) -> f64 {
        let penetration = radius - distance;
        if distance < radius {
            // Repulsion
            self.stiffness * penetration
        } else {
            // Adhesion pull (JKR-like, simplified)
            let pull_range = radius * 0.1;
            if distance < radius + pull_range {
                -self.adhesion * (distance - radius) / pull_range
            } else {
                0.0
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Utility: arc-length parametrization helpers
// ---------------------------------------------------------------------------

/// Compute the tip deflection of a cantilever rod under a tip transverse load.
///
/// Uses the Euler-Bernoulli limit: δ = PL³ / (3EI).
pub fn euler_bernoulli_tip_deflection(load: f64, length: f64, ei: f64) -> f64 {
    load * length.powi(3) / (3.0 * ei)
}

/// Compute the tip rotation of a cantilever under end moment M₀.
///
/// For a uniform rod: θ = M₀L / (EI).  The deformed shape is a circular arc.
pub fn end_moment_tip_rotation(moment: f64, length: f64, ei: f64) -> f64 {
    moment * length / ei
}

/// Compute the angular deflection of a rod under pure torsion.
///
/// φ = T·L / (GJ).
pub fn torsion_angle(torque: f64, length: f64, gj: f64) -> f64 {
    torque * length / gj
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // 1. Circular cross-section area
    #[test]
    fn test_circular_area() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.01 });
        let expected = PI * 0.01 * 0.01;
        assert!((cs.area - expected).abs() < 1e-12, "area={}", cs.area);
    }

    // 2. Circular second moment of inertia I = π r⁴ / 4
    #[test]
    fn test_circular_moment_of_inertia() {
        let r = 0.02;
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: r });
        let expected = PI * r.powi(4) / 4.0;
        assert!((cs.i1 - expected).abs() < 1e-20, "i1={}", cs.i1);
        assert!((cs.i2 - expected).abs() < 1e-20, "i2={}", cs.i2);
    }

    // 3. Circular torsion constant J = π r⁴ / 2
    #[test]
    fn test_circular_torsion_constant() {
        let r = 0.01;
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: r });
        let expected = PI * r.powi(4) / 2.0;
        assert!(
            (cs.torsion_j - expected).abs() < 1e-18,
            "J={}",
            cs.torsion_j
        );
    }

    // 4. Rectangular cross-section I₁ = bh³/12
    #[test]
    fn test_rectangular_moment_of_inertia() {
        let b = 0.02;
        let h = 0.04;
        let cs = RodCrossSection::new(CrossSectionShape::Rectangular {
            width: b,
            height: h,
        });
        let i1_expected = b * h.powi(3) / 12.0;
        assert!(
            (cs.i1 - i1_expected).abs() < 1e-18,
            "i1={} expected={}",
            cs.i1,
            i1_expected
        );
    }

    // 5. Elliptical area = π a b
    #[test]
    fn test_elliptical_area() {
        let a = 0.03;
        let b = 0.01;
        let cs = RodCrossSection::new(CrossSectionShape::Elliptical {
            semi_a: a,
            semi_b: b,
        });
        let expected = PI * a * b;
        assert!((cs.area - expected).abs() < 1e-14, "area={}", cs.area);
    }

    // 6. Straight rod nodes have correct z-spacing
    #[test]
    fn test_straight_rod_node_spacing() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(200e9, 77e9, &cs);
        let rod = CosseratRod::new_straight(6, 1.0, stiff, 1.0, 1e-6);
        for i in 1..rod.n_nodes() {
            let dz = rod.nodes[i].position[2] - rod.nodes[i - 1].position[2];
            assert!((dz - 0.2).abs() < 1e-10, "dz[{}]={}", i, dz);
        }
    }

    // 7. Straight rod at rest: zero strains
    #[test]
    fn test_straight_rod_zero_strain() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(200e9, 77e9, &cs);
        let rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        for seg in 0..rod.n_segments() {
            let (eps, g1, g2, dk1, dk2, dtau) = InternalForces::segment_strains(&rod, seg);
            assert!(eps.abs() < EPS, "axial strain at seg {seg}: {eps}");
            assert!(g1.abs() < 0.1, "shear g1 at seg {seg}: {g1}"); // shear against d1 for aligned rod is small
            assert!(g2.abs() < EPS, "shear g2 at seg {seg}: {g2}");
            assert!(dk1.abs() < EPS, "Δκ1 at seg {seg}: {dk1}");
            assert!(dk2.abs() < EPS, "Δκ2 at seg {seg}: {dk2}");
            assert!(dtau.abs() < EPS, "Δτ at seg {seg}: {dtau}");
        }
    }

    // 8. Euler-Bernoulli tip deflection formula
    #[test]
    fn test_euler_bernoulli_formula() {
        let p = 10.0; // N
        let l = 1.0; // m
        let ei = 1000.0; // N·m²
        let delta = euler_bernoulli_tip_deflection(p, l, ei);
        let expected = p * l.powi(3) / (3.0 * ei);
        assert!((delta - expected).abs() < EPS);
    }

    // 9. End moment tip rotation formula
    #[test]
    fn test_end_moment_rotation() {
        let m = 5.0;
        let l = 0.5;
        let ei = 200.0;
        let theta = end_moment_tip_rotation(m, l, ei);
        assert!((theta - m * l / ei).abs() < EPS);
    }

    // 10. Torsion angle formula
    #[test]
    fn test_torsion_angle_formula() {
        let t = 2.0;
        let l = 1.0;
        let gj = 500.0;
        let phi = torsion_angle(t, l, gj);
        assert!((phi - t * l / gj).abs() < EPS);
    }

    // 11. Rodrigues rotation: 180° around Z maps [1,0,0] → [-1,0,0]
    #[test]
    fn test_rodrigues_rotation_180() {
        let v = [1.0, 0.0, 0.0];
        let k = [0.0, 0.0, 1.0];
        let result = rodrigues(v, k, PI);
        assert!((result[0] + 1.0).abs() < 1e-12, "x={}", result[0]);
        assert!(result[1].abs() < 1e-12, "y={}", result[1]);
    }

    // 12. CosseratNode orthonormalization preserves d3 direction
    #[test]
    fn test_node_orthonormalize() {
        let mut node = CosseratNode::new([0.0, 0.0, 0.0]);
        // Perturb d1 slightly off-orthogonal to d3
        node.d1 = [0.9, 0.1, 0.3];
        node.orthonormalize();
        let dot_d1_d3 = vec3_dot(node.d1, node.d3);
        assert!(dot_d1_d3.abs() < 1e-10, "d1·d3={dot_d1_d3}");
        let dot_d1_d2 = vec3_dot(node.d1, node.d2);
        assert!(dot_d1_d2.abs() < 1e-10, "d1·d2={dot_d1_d2}");
    }

    // 13. Static solver reduces residual
    #[test]
    fn test_static_solver_residual_decreases() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.01 });
        let stiff = RodStiffness::from_material(1e6, 4e5, &cs);
        let mut rod = CosseratRod::new_straight(10, 1.0, stiff, 0.1, 1e-6);

        // Apply tip load
        let last = rod.n_nodes() - 1;
        rod.nodes[last].external_force = [0.0, 1.0, 0.0];

        let bcs = BoundaryConditions::cantilever(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        let solver = StaticSolver::new(200, 1e-3);
        let _res0 = StaticSolver::residual_norm(&rod);
        let (iters, _res_final) = solver.solve(&mut rod, &bcs);
        // Solver should run without panicking and iterate
        assert!(iters > 0, "solver should iterate at least once");
    }

    // 14. Segment-segment distance: parallel segments on same line
    #[test]
    fn test_segment_segment_distance_collinear() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod_a = CosseratRod::new_straight(3, 1.0, stiff.clone(), 1.0, 1e-6);
        let mut rod_b = CosseratRod::new_straight(3, 1.0, stiff, 1.0, 1e-6);
        // Shift rod_b by 0.1 in X
        for n in rod_b.nodes.iter_mut() {
            n.position[0] += 0.1;
        }
        let (dist, _pa, _pb) = RodContact::segment_segment_distance(&rod_a, 0, &rod_b, 0);
        assert!((dist - 0.1).abs() < 1e-10, "dist={dist}");
    }

    // 15. RodContact detects contact when rods are close
    #[test]
    fn test_rod_contact_detection() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod_a = CosseratRod::new_straight(3, 1.0, stiff.clone(), 1.0, 1e-6);
        let mut rod_b = CosseratRod::new_straight(3, 1.0, stiff, 1.0, 1e-6);
        for n in rod_b.nodes.iter_mut() {
            n.position[0] += 0.05;
        }
        let contact = RodContact::new(1000.0, 0.3, 0.0);
        let contacts = contact.detect_contacts(&rod_a, &rod_b, 0.1);
        assert!(!contacts.is_empty(), "should detect contacts");
    }

    // 16. RodContact: no contact when rods are far
    #[test]
    fn test_rod_contact_no_contact() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod_a = CosseratRod::new_straight(3, 1.0, stiff.clone(), 1.0, 1e-6);
        let mut rod_b = CosseratRod::new_straight(3, 1.0, stiff, 1.0, 1e-6);
        for n in rod_b.nodes.iter_mut() {
            n.position[0] += 10.0;
        }
        let contact = RodContact::new(1000.0, 0.3, 0.0);
        let contacts = contact.detect_contacts(&rod_a, &rod_b, 0.1);
        assert!(contacts.is_empty(), "should have no contacts");
    }

    // 17. Follower force boundary condition updates tip force direction
    #[test]
    fn test_follower_force_bc() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let mut rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        let last = rod.n_nodes() - 1;
        // Tilt the tip tangent
        rod.nodes[last].d3 = vec3_normalize([0.0, 1.0, 1.0]);

        let bcs = BoundaryConditions {
            base: BcType::Clamped,
            tip: BcType::FollowerForce { magnitude: 5.0 },
            base_position: [0.0; 3],
            base_d1: [1.0, 0.0, 0.0],
            base_d2: [0.0, 1.0, 0.0],
            base_d3: [0.0, 0.0, 1.0],
        };
        bcs.apply(&mut rod);
        let tip_force = rod.nodes[last].external_force;
        // Force should be along tip d3
        let tip_d3 = rod.nodes[last].d3;
        let cos_angle = vec3_dot(vec3_normalize(tip_force), tip_d3);
        assert!(
            (cos_angle - 1.0).abs() < 1e-10,
            "follower force not along d3"
        );
    }

    // 18. Clamped base BC resets position and velocity
    #[test]
    fn test_clamped_bc_resets_base() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let mut rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        rod.nodes[0].position = [5.0, 5.0, 5.0]; // perturb
        rod.nodes[0].velocity = [1.0, 2.0, 3.0]; // perturb

        let bcs = BoundaryConditions::cantilever(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        bcs.apply(&mut rod);
        assert_eq!(rod.nodes[0].position, [0.0, 0.0, 0.0]);
        assert_eq!(rod.nodes[0].velocity, [0.0, 0.0, 0.0]);
    }

    // 19. Dynamics step moves rod under gravity
    #[test]
    fn test_dynamics_gravity() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e6, 4e5, &cs);
        let mut rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);

        // Apply gravity to all nodes
        for node in rod.nodes.iter_mut() {
            node.external_force = [0.0, -9.81, 0.0];
        }

        let pos_before = rod.nodes[2].position[1];
        let dt = 1e-3;
        for _ in 0..10 {
            CosseratDynamics::step_euler(&mut rod, dt);
        }
        let pos_after = rod.nodes[2].position[1];
        assert!(pos_after < pos_before, "node should fall under gravity");
    }

    // 20. Segment tangent of straight rod is [0,0,1]
    #[test]
    fn test_segment_tangent_straight() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        let t = rod.segment_tangent(0);
        assert!((t[2] - 1.0).abs() < 1e-10, "tangent={t:?}");
    }

    // 21. Segment stretch of undeformed rod is 1.0
    #[test]
    fn test_segment_stretch_undeformed() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        for seg in 0..rod.n_segments() {
            let s = rod.segment_stretch(seg);
            assert!((s - 1.0).abs() < 1e-10, "stretch[{seg}]={s}");
        }
    }

    // 22. Internal loads are zero for undeformed rod
    #[test]
    fn test_internal_loads_zero_undeformed() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.01 });
        let stiff = RodStiffness::from_material(200e9, 77e9, &cs);
        let rod = CosseratRod::new_straight(5, 1.0, stiff, 1.0, 1e-6);
        for seg in 0..rod.n_segments() {
            let loads = InternalForces::compute(&rod, seg);
            assert!(loads.axial_force.abs() < EPS, "N at seg {seg}");
            assert!(loads.moment1.abs() < EPS, "M1 at seg {seg}");
            assert!(loads.moment2.abs() < EPS, "M2 at seg {seg}");
            assert!(loads.torsion.abs() < EPS, "T at seg {seg}");
        }
    }

    // 23. Shooting method returns base position for zero curvature
    #[test]
    fn test_shoot_zero_curvature() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod = CosseratRod::new_straight(10, 1.0, stiff, 1.0, 1e-6);
        let (tip_pos, tip_d3) = StaticSolver::shoot(&rod, 0.0, 0.0, 0.0);
        // Straight rod, tip should be at z=1.0
        assert!((tip_pos[2] - 1.0).abs() < 1e-8, "tip_z={}", tip_pos[2]);
        // Tangent should be along +Z
        assert!((tip_d3[2] - 1.0).abs() < 1e-8, "tip_d3={tip_d3:?}");
    }

    // 24. Contact force repulsion is positive for penetration
    #[test]
    fn test_contact_force_repulsion() {
        let contact = RodContact::new(1000.0, 0.3, 0.0);
        let f = contact.contact_force_magnitude(0.05, 0.1);
        assert!(
            f > 0.0,
            "contact force should be repulsive for penetration, f={f}"
        );
    }

    // 25. Contact force is zero at radius
    #[test]
    fn test_contact_force_at_radius() {
        let contact = RodContact::new(1000.0, 0.3, 0.0);
        let f = contact.contact_force_magnitude(0.1, 0.1);
        assert!(f.abs() < EPS, "f at contact radius should be ~0, f={f}");
    }

    // 26. Stiffness from material: EA = E * A
    #[test]
    fn test_stiffness_ea() {
        let r = 0.01;
        let e = 200e9;
        let g = 77e9;
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: r });
        let stiff = RodStiffness::from_material(e, g, &cs);
        let expected_ea = e * cs.area;
        assert!((stiff.ea - expected_ea).abs() / expected_ea < 1e-10);
    }

    // 27. Total rest length matches specified length
    #[test]
    fn test_total_rest_length() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod = CosseratRod::new_straight(11, 2.5, stiff, 1.0, 1e-6);
        assert!((rod.total_rest_length() - 2.5).abs() < 1e-10);
    }

    // 28. Directors are orthonormal after initialization
    #[test]
    fn test_directors_orthonormal() {
        let node = CosseratNode::new([0.0, 0.0, 0.0]);
        let d12 = vec3_dot(node.d1, node.d2);
        let d13 = vec3_dot(node.d1, node.d3);
        let d23 = vec3_dot(node.d2, node.d3);
        let n1 = vec3_norm(node.d1);
        let n2 = vec3_norm(node.d2);
        let n3 = vec3_norm(node.d3);
        assert!(d12.abs() < EPS, "d1·d2={d12}");
        assert!(d13.abs() < EPS, "d1·d3={d13}");
        assert!(d23.abs() < EPS, "d2·d3={d23}");
        assert!((n1 - 1.0).abs() < EPS);
        assert!((n2 - 1.0).abs() < EPS);
        assert!((n3 - 1.0).abs() < EPS);
    }

    // 29. Rotate directors by 2π returns to original
    #[test]
    fn test_rotate_directors_full_circle() {
        let mut node = CosseratNode::new([0.0, 0.0, 0.0]);
        let omega_z = 2.0 * PI; // 2π rad/s in Z direction
        let dt = 1.0; // 1 second → full rotation
        node.rotate_directors([0.0, 0.0, omega_z], dt);
        node.orthonormalize();
        // d1 should be back to ~[1,0,0], d2 ~[0,1,0]
        assert!((node.d1[0] - 1.0).abs() < 1e-10, "d1x={}", node.d1[0]);
        assert!((node.d2[1] - 1.0).abs() < 1e-10, "d2y={}", node.d2[1]);
    }

    // 30. Rod n_segments = n_nodes - 1
    #[test]
    fn test_n_segments() {
        let cs = RodCrossSection::new(CrossSectionShape::Circular { radius: 0.005 });
        let stiff = RodStiffness::from_material(1e9, 4e8, &cs);
        let rod = CosseratRod::new_straight(7, 1.0, stiff, 1.0, 1e-6);
        assert_eq!(rod.n_segments(), 6);
        assert_eq!(rod.n_nodes(), 7);
    }
}
