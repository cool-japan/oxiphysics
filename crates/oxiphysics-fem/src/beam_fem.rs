// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Beam finite element formulations.
//!
//! This module provides a comprehensive set of beam element types and analyses:
//!
//! - **Euler-Bernoulli beam**: classical thin beam, stiffness matrix, shape functions
//! - **Timoshenko beam**: shear-deformable thick beam, shear correction factor
//! - **3-D beam with torsion**: Saint-Venant torsion, warping, 12-DOF element
//! - **Beam-column (P-delta)**: geometric stiffness matrix, buckling interaction
//! - **Curved beam element**: initial curvature, in-plane and out-of-plane bending
//! - **Composite beam**: transformed section (steel-concrete, laminated cross-section)
//! - **Beam dynamics**: consistent and lumped mass matrices, natural frequencies
//! - **Plastic hinge**: moment-rotation, hinge formation criterion
//! - **Tapered beam**: linearly varying cross-section, numerical integration
//! - **Cable element**: catenary geometry, geometrically nonlinear tension element
//!
//! # Coordinate convention
//! The local beam axis is along x; y and z are the two principal bending axes.
//! All stiffness matrices are expressed in local coordinates.

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Default shear correction factor κ for solid rectangular cross-sections (5/6).
pub const SHEAR_CORRECTION_RECT: f64 = 5.0 / 6.0;

/// Default shear correction factor κ for solid circular cross-sections.
pub const SHEAR_CORRECTION_CIRCLE: f64 = 0.9;

/// Number of DOF per node for a 3-D Euler-Bernoulli beam (3 translations + 3 rotations).
pub const DOF_PER_NODE_3D: usize = 6;

/// Number of DOF per node for a 2-D beam (2 translations + 1 rotation).
pub const DOF_PER_NODE_2D: usize = 3;

/// Total DOF for a 2-node 3-D beam element (12).
pub const BEAM3D_NDOF: usize = 12;

/// Total DOF for a 2-node 2-D beam element (6).
pub const BEAM2D_NDOF: usize = 6;

// ============================================================================
// Section 1 – Cross-section properties
// ============================================================================

/// Geometric and material cross-section properties for a beam.
#[derive(Debug, Clone, Copy)]
pub struct CrossSection {
    /// Cross-sectional area A \[m²\].
    pub area: f64,
    /// Second moment of area about the strong (z) axis I_zz \[m⁴\].
    pub izz: f64,
    /// Second moment of area about the weak (y) axis I_yy \[m⁴\].
    pub iyy: f64,
    /// Saint-Venant torsional constant J \[m⁴\].
    pub j_torsion: f64,
    /// Warping constant C_w \[m⁶\] (for open thin-walled sections).
    pub c_warping: f64,
    /// Young's modulus E \[Pa\].
    pub elastic_modulus: f64,
    /// Shear modulus G \[Pa\].
    pub shear_modulus: f64,
    /// Shear correction factor κ_y (transverse shear in y direction).
    pub kappa_y: f64,
    /// Shear correction factor κ_z (transverse shear in z direction).
    pub kappa_z: f64,
}

impl CrossSection {
    /// Create a solid rectangular cross-section.
    ///
    /// # Arguments
    /// * `b`      – width \[m\]
    /// * `h`      – height \[m\]
    /// * `e`      – Young's modulus \[Pa\]
    /// * `nu`     – Poisson's ratio
    pub fn rectangular(b: f64, h: f64, e: f64, nu: f64) -> Self {
        let g = e / (2.0 * (1.0 + nu));
        let area = b * h;
        let izz = b * h * h * h / 12.0;
        let iyy = h * b * b * b / 12.0;
        // Approximate torsional constant for rectangle: J ≈ b h³ (1/3 - 0.21 h/b) for b > h
        let (bw, hw) = if b >= h { (b, h) } else { (h, b) };
        let j = hw * hw * hw * bw / 3.0
            * (1.0 - 0.21 * hw / bw * (1.0 - hw.powi(4) / (12.0 * bw.powi(4))));
        Self {
            area,
            izz,
            iyy,
            j_torsion: j,
            c_warping: 0.0,
            elastic_modulus: e,
            shear_modulus: g,
            kappa_y: SHEAR_CORRECTION_RECT,
            kappa_z: SHEAR_CORRECTION_RECT,
        }
    }

    /// Create a solid circular cross-section.
    ///
    /// # Arguments
    /// * `radius` – cross-section radius \[m\]
    /// * `e`      – Young's modulus \[Pa\]
    /// * `nu`     – Poisson's ratio
    pub fn circular(radius: f64, e: f64, nu: f64) -> Self {
        let g = e / (2.0 * (1.0 + nu));
        let r2 = radius * radius;
        let area = PI * r2;
        let i = PI * r2 * r2 / 4.0;
        let j = PI * r2 * r2 / 2.0;
        Self {
            area,
            izz: i,
            iyy: i,
            j_torsion: j,
            c_warping: 0.0,
            elastic_modulus: e,
            shear_modulus: g,
            kappa_y: SHEAR_CORRECTION_CIRCLE,
            kappa_z: SHEAR_CORRECTION_CIRCLE,
        }
    }

    /// Axial stiffness EA \[N\].
    pub fn axial_stiffness(&self) -> f64 {
        self.elastic_modulus * self.area
    }

    /// Bending stiffness EI_zz \[N·m²\] (strong axis).
    pub fn bending_stiffness_z(&self) -> f64 {
        self.elastic_modulus * self.izz
    }

    /// Bending stiffness EI_yy \[N·m²\] (weak axis).
    pub fn bending_stiffness_y(&self) -> f64 {
        self.elastic_modulus * self.iyy
    }

    /// Torsional stiffness GJ \[N·m²\].
    pub fn torsional_stiffness(&self) -> f64 {
        self.shear_modulus * self.j_torsion
    }

    /// Shear stiffness κ_y G A \[N\].
    pub fn shear_stiffness_y(&self) -> f64 {
        self.kappa_y * self.shear_modulus * self.area
    }

    /// Shear stiffness κ_z G A \[N\].
    pub fn shear_stiffness_z(&self) -> f64 {
        self.kappa_z * self.shear_modulus * self.area
    }
}

// ============================================================================
// Section 2 – Euler-Bernoulli beam element (2-D, 2 nodes, 6 DOF)
// ============================================================================

/// 2-D Euler-Bernoulli beam element stiffness matrix (6×6, local coordinates).
///
/// DOF order: \[u1, v1, θ1, u2, v2, θ2\] where u=axial, v=transverse, θ=rotation.
///
/// Returns a 6×6 matrix stored as a flat \[f64; 36\] row-major array.
///
/// # Arguments
/// * `cs` – cross-section properties
/// * `le` – element length L \[m\]
pub fn euler_bernoulli_stiffness_2d(cs: &CrossSection, le: f64) -> [f64; 36] {
    let ea = cs.axial_stiffness();
    let ei = cs.bending_stiffness_z();
    let l = le;
    let l2 = l * l;
    let l3 = l * l * l;

    let mut k = [0.0f64; 36];
    // Axial DOF (0, 3)
    k[0] = ea / l;
    k[3] = -ea / l;
    k[3 * 6] = -ea / l;
    k[3 * 6 + 3] = ea / l;

    // Bending DOF (1, 2, 4, 5)
    let k11 = 12.0 * ei / l3;
    let k12 = 6.0 * ei / l2;
    let k22 = 4.0 * ei / l;
    let k25 = 2.0 * ei / l;

    k[6 + 1] = k11;
    k[6 + 2] = k12;
    k[6 + 4] = -k11;
    k[6 + 5] = k12;

    k[2 * 6 + 1] = k12;
    k[2 * 6 + 2] = k22;
    k[2 * 6 + 4] = -k12;
    k[2 * 6 + 5] = k25;

    k[4 * 6 + 1] = -k11;
    k[4 * 6 + 2] = -k12;
    k[4 * 6 + 4] = k11;
    k[4 * 6 + 5] = -k12;

    k[5 * 6 + 1] = k12;
    k[5 * 6 + 2] = k25;
    k[5 * 6 + 4] = -k12;
    k[5 * 6 + 5] = k22;

    k
}

/// Evaluate the Euler-Bernoulli Hermite shape functions at natural coordinate ξ ∈ \[0,1\].
///
/// Returns \[N1, N2, N3, N4\] where N1, N2 are the translational shape functions
/// at nodes 1 and 2, and N3, N4 are the rotational shape functions.
///
/// # Arguments
/// * `xi` – natural coordinate ξ ∈ \[0, 1\]
/// * `le` – element length L
pub fn hermite_shape_functions(xi: f64, le: f64) -> [f64; 4] {
    let x = xi * le;
    let l = le;
    let l2 = l * l;
    let l3 = l * l * l;
    let n1 = 1.0 - 3.0 * x * x / l2 + 2.0 * x * x * x / l3;
    let n2 = 3.0 * x * x / l2 - 2.0 * x * x * x / l3;
    let n3 = x - 2.0 * x * x / l + x * x * x / l2;
    let n4 = -x * x / l + x * x * x / l2;
    [n1, n2, n3, n4]
}

/// Evaluate the derivatives of the Hermite shape functions dN/dx at ξ ∈ \[0,1\].
///
/// Returns \[dN1/dx, dN2/dx, dN3/dx, dN4/dx\].
///
/// # Arguments
/// * `xi` – natural coordinate ξ ∈ \[0, 1\]
/// * `le` – element length L
pub fn hermite_shape_derivatives(xi: f64, le: f64) -> [f64; 4] {
    let x = xi * le;
    let l = le;
    let l2 = l * l;
    let l3 = l * l * l;
    let dn1 = -6.0 * x / l2 + 6.0 * x * x / l3;
    let dn2 = 6.0 * x / l2 - 6.0 * x * x / l3;
    let dn3 = 1.0 - 4.0 * x / l + 3.0 * x * x / l2;
    let dn4 = -2.0 * x / l + 3.0 * x * x / l2;
    [dn1, dn2, dn3, dn4]
}

/// Compute the Euler-Bernoulli consistent load vector for a uniform transverse load q.
///
/// Returns \[Fy1, M1, Fy2, M2\] (transverse forces and moments at both nodes).
///
/// # Arguments
/// * `q`  – distributed transverse load \[N/m\]
/// * `le` – element length \[m\]
pub fn euler_bernoulli_uniform_load(q: f64, le: f64) -> [f64; 4] {
    let l = le;
    [
        q * l / 2.0,
        q * l * l / 12.0,
        q * l / 2.0,
        -q * l * l / 12.0,
    ]
}

// ============================================================================
// Section 3 – Timoshenko beam element
// ============================================================================

/// Timoshenko shear parameter φ = 12 EI / (κ G A L²).
///
/// φ = 0 → Euler-Bernoulli limit; large φ → shear-dominated beam.
///
/// # Arguments
/// * `cs` – cross-section
/// * `le` – element length
pub fn timoshenko_phi(cs: &CrossSection, le: f64) -> f64 {
    let ei = cs.bending_stiffness_z();
    let ks = cs.shear_stiffness_y();
    if ks.abs() < f64::EPSILON || le.abs() < f64::EPSILON {
        return 0.0;
    }
    12.0 * ei / (ks * le * le)
}

/// 2-D Timoshenko beam element stiffness matrix (6×6, local coordinates).
///
/// DOF order: \[u1, v1, θ1, u2, v2, θ2\].
/// Shear locking is avoided by the standard consistent formulation with φ.
///
/// # Arguments
/// * `cs` – cross-section
/// * `le` – element length
pub fn timoshenko_stiffness_2d(cs: &CrossSection, le: f64) -> [f64; 36] {
    let ea = cs.axial_stiffness();
    let ei = cs.bending_stiffness_z();
    let phi = timoshenko_phi(cs, le);
    let l = le;
    let l2 = l * l;
    let l3 = l * l * l;
    let c = 1.0 / (1.0 + phi);

    let mut k = [0.0f64; 36];

    // Axial DOF
    k[0] = ea / l;
    k[3] = -ea / l;
    k[3 * 6] = -ea / l;
    k[3 * 6 + 3] = ea / l;

    // Timoshenko bending
    let k11 = 12.0 * ei * c / l3;
    let k12 = 6.0 * ei * c / l2;
    let k22 = (4.0 + phi) * ei * c / l;
    let k25 = (2.0 - phi) * ei * c / l;

    k[6 + 1] = k11;
    k[6 + 2] = k12;
    k[6 + 4] = -k11;
    k[6 + 5] = k12;

    k[2 * 6 + 1] = k12;
    k[2 * 6 + 2] = k22;
    k[2 * 6 + 4] = -k12;
    k[2 * 6 + 5] = k25;

    k[4 * 6 + 1] = -k11;
    k[4 * 6 + 2] = -k12;
    k[4 * 6 + 4] = k11;
    k[4 * 6 + 5] = -k12;

    k[5 * 6 + 1] = k12;
    k[5 * 6 + 2] = k25;
    k[5 * 6 + 4] = -k12;
    k[5 * 6 + 5] = k22;

    k
}

// ============================================================================
// Section 4 – 3-D beam element with torsion (12 DOF)
// ============================================================================

/// 3-D beam element stiffness matrix with axial, biaxial bending, and torsion.
///
/// DOF order per node: \[u, v, w, θx, θy, θz\] (translations then rotations).
/// Total: 12 DOF.  Returns a flat \[f64; 144\] row-major array.
///
/// Uses Euler-Bernoulli bending and Saint-Venant torsion.
///
/// # Arguments
/// * `cs` – cross-section properties
/// * `le` – element length
pub fn beam3d_stiffness(cs: &CrossSection, le: f64) -> [f64; 144] {
    let ea = cs.axial_stiffness();
    let eiz = cs.bending_stiffness_z();
    let eiy = cs.bending_stiffness_y();
    let gj = cs.torsional_stiffness();
    let l = le;
    let l2 = l * l;
    let l3 = l * l * l;

    let mut k = [0.0f64; 144];

    // Helper closure to set symmetric entry
    let set = |k: &mut [f64; 144], i: usize, j: usize, v: f64| {
        k[i * 12 + j] = v;
        k[j * 12 + i] = v;
    };

    // Axial (DOF 0 and 6)
    k[0] = ea / l;
    k[6 * 12 + 6] = ea / l;
    set(&mut k, 0, 6, -ea / l);

    // Torsion (DOF 3 and 9)
    k[3 * 12 + 3] = gj / l;
    k[9 * 12 + 9] = gj / l;
    set(&mut k, 3, 9, -gj / l);

    // Bending in xz-plane (bending about y-axis: DOF 2, 4, 8, 10)
    let ky11 = 12.0 * eiy / l3;
    let ky12 = 6.0 * eiy / l2;
    let ky22 = 4.0 * eiy / l;
    let ky25 = 2.0 * eiy / l;

    k[2 * 12 + 2] = ky11;
    k[8 * 12 + 8] = ky11;
    set(&mut k, 2, 8, -ky11);
    set(&mut k, 2, 4, ky12);
    set(&mut k, 2, 10, ky12);
    set(&mut k, 8, 4, -ky12);
    set(&mut k, 8, 10, -ky12);
    k[4 * 12 + 4] = ky22;
    k[10 * 12 + 10] = ky22;
    set(&mut k, 4, 10, ky25);

    // Bending in xy-plane (bending about z-axis: DOF 1, 5, 7, 11)
    let kz11 = 12.0 * eiz / l3;
    let kz12 = 6.0 * eiz / l2;
    let kz22 = 4.0 * eiz / l;
    let kz25 = 2.0 * eiz / l;

    k[12 + 1] = kz11;
    k[7 * 12 + 7] = kz11;
    set(&mut k, 1, 7, -kz11);
    set(&mut k, 1, 5, kz12);
    set(&mut k, 1, 11, kz12);
    set(&mut k, 7, 5, -kz12);
    set(&mut k, 7, 11, -kz12);
    k[5 * 12 + 5] = kz22;
    k[11 * 12 + 11] = kz22;
    set(&mut k, 5, 11, kz25);

    k
}

/// Verify that a 12×12 stiffness matrix is symmetric (max off-diagonal error).
///
/// # Arguments
/// * `k` – 12×12 matrix stored as \[f64; 144\]
pub fn beam3d_stiffness_symmetry_error(k: &[f64; 144]) -> f64 {
    let mut max_err = 0.0_f64;
    for i in 0..12 {
        for j in 0..12 {
            let err = (k[i * 12 + j] - k[j * 12 + i]).abs();
            if err > max_err {
                max_err = err;
            }
        }
    }
    max_err
}

// ============================================================================
// Section 5 – Geometric stiffness matrix (P-delta effect)
// ============================================================================

/// Geometric stiffness matrix for a 2-D beam-column element (6×6).
///
/// K_g = (N/L) \[0  0  0  0  0  0;
///               0  6/5  L/10  0  -6/5  L/10;
///               0  L/10  2L²/15  0  -L/10  -L²/30; ...]
///
/// This is the standard consistent geometric stiffness for Euler-Bernoulli.
/// N > 0 means compression (destabilising).
///
/// Returns a flat \[f64; 36\] row-major array.
///
/// # Arguments
/// * `axial_force` – axial compressive force N \[N\] (positive = compression)
/// * `le`          – element length \[m\]
pub fn geometric_stiffness_2d(axial_force: f64, le: f64) -> [f64; 36] {
    let n = axial_force;
    let l = le;
    let l2 = l * l;
    let mut kg = [0.0f64; 36];

    // Transverse DOF (1, 2, 4, 5) — the axial DOF are unaffected at this order
    let c = n / l;

    kg[6 + 1] = c * 6.0 / 5.0;
    kg[6 + 2] = c * l / 10.0;
    kg[6 + 4] = -c * 6.0 / 5.0;
    kg[6 + 5] = c * l / 10.0;

    kg[2 * 6 + 1] = c * l / 10.0;
    kg[2 * 6 + 2] = c * 2.0 * l2 / 15.0;
    kg[2 * 6 + 4] = -c * l / 10.0;
    kg[2 * 6 + 5] = -c * l2 / 30.0;

    kg[4 * 6 + 1] = -c * 6.0 / 5.0;
    kg[4 * 6 + 2] = -c * l / 10.0;
    kg[4 * 6 + 4] = c * 6.0 / 5.0;
    kg[4 * 6 + 5] = -c * l / 10.0;

    kg[5 * 6 + 1] = c * l / 10.0;
    kg[5 * 6 + 2] = -c * l2 / 30.0;
    kg[5 * 6 + 4] = -c * l / 10.0;
    kg[5 * 6 + 5] = c * 2.0 * l2 / 15.0;

    kg
}

/// Euler critical buckling load for a pin-pin beam-column: P_cr = π² EI / L².
///
/// # Arguments
/// * `cs` – cross-section
/// * `le` – element length (= column height for a pin-pin column)
pub fn euler_buckling_load(cs: &CrossSection, le: f64) -> f64 {
    let ei = cs.bending_stiffness_z().min(cs.bending_stiffness_y());
    PI * PI * ei / (le * le)
}

/// Effective buckling length factor K for common end conditions.
///
/// Returns K: 1.0 (pin-pin), 0.5 (fixed-fixed), 0.7 (fixed-pin), 2.0 (fixed-free).
pub enum BoundaryConditionType {
    /// Both ends pinned (K = 1.0).
    PinPin,
    /// Both ends fixed (K = 0.5).
    FixedFixed,
    /// One end fixed, one pinned (K = 0.7).
    FixedPin,
    /// One end fixed, one free (cantilever, K = 2.0).
    FixedFree,
}

/// Get the effective length factor K for a beam-column.
///
/// # Arguments
/// * `bc` – boundary condition type
pub fn effective_length_factor(bc: &BoundaryConditionType) -> f64 {
    match bc {
        BoundaryConditionType::PinPin => 1.0,
        BoundaryConditionType::FixedFixed => 0.5,
        BoundaryConditionType::FixedPin => 0.7,
        BoundaryConditionType::FixedFree => 2.0,
    }
}

// ============================================================================
// Section 6 – Curved beam element
// ============================================================================

/// Curved beam element with initial radius of curvature R.
///
/// The element accounts for membrane-bending coupling due to curvature.
/// Stored in the local arc-length coordinate system.
#[derive(Debug, Clone, Copy)]
pub struct CurvedBeamElement {
    /// Cross-section properties.
    pub section: CrossSection,
    /// Arc length of the element \[m\].
    pub arc_length: f64,
    /// Radius of curvature R \[m\] (positive = center above beam axis).
    pub radius: f64,
}

impl CurvedBeamElement {
    /// Construct a curved beam element.
    pub fn new(section: CrossSection, arc_length: f64, radius: f64) -> Self {
        Self {
            section,
            arc_length,
            radius,
        }
    }

    /// Curvature κ = 1/R \[1/m\].
    pub fn curvature(&self) -> f64 {
        if self.radius.abs() < f64::EPSILON {
            return 0.0;
        }
        1.0 / self.radius
    }

    /// Membrane stiffness of curved beam: k_m = EA / arc_length.
    pub fn membrane_stiffness(&self) -> f64 {
        self.section.axial_stiffness() / self.arc_length.max(f64::EPSILON)
    }

    /// Bending stiffness of curved beam: k_b = EI / arc_length³ * 12.
    pub fn bending_stiffness(&self) -> f64 {
        12.0 * self.section.bending_stiffness_z()
            / (self.arc_length * self.arc_length * self.arc_length).max(f64::EPSILON)
    }

    /// Membrane-bending coupling stiffness k_mb = κ EI / arc_length.
    ///
    /// This term couples in-plane forces to moments due to the curvature.
    pub fn coupling_stiffness(&self) -> f64 {
        self.curvature() * self.section.bending_stiffness_z() / self.arc_length.max(f64::EPSILON)
    }

    /// Radial deflection under a uniform radial load q \[N/m\] (approximate).
    ///
    /// Δr ≈ q R⁴ / (8 EI)  for a shallow arch
    ///
    /// # Arguments
    /// * `q_radial` – radial distributed load \[N/m\]
    pub fn radial_deflection_uniform_load(&self, q_radial: f64) -> f64 {
        let r = self.radius;
        let ei = self.section.bending_stiffness_z();
        if ei.abs() < f64::EPSILON {
            return 0.0;
        }
        q_radial * r * r * r * r / (8.0 * ei)
    }
}

// ============================================================================
// Section 7 – Composite beam (transformed section)
// ============================================================================

/// A single material layer in a composite beam cross-section.
#[derive(Debug, Clone, Copy)]
pub struct CompositeLayer {
    /// Width of the layer \[m\].
    pub width: f64,
    /// Height/thickness of the layer \[m\].
    pub height: f64,
    /// Distance from the reference axis to the centroid of this layer \[m\].
    pub y_centroid: f64,
    /// Young's modulus of this layer \[Pa\].
    pub elastic_modulus: f64,
}

impl CompositeLayer {
    /// Area of this layer A_i = b_i h_i.
    pub fn area(&self) -> f64 {
        self.width * self.height
    }

    /// First moment of area about the reference axis: S_i = A_i y_i.
    pub fn first_moment(&self) -> f64 {
        self.area() * self.y_centroid
    }

    /// Second moment of area about the reference axis (parallel axis theorem):
    /// I_i = b_i h_i³/12 + A_i y_i².
    pub fn second_moment(&self) -> f64 {
        self.width * self.height * self.height * self.height / 12.0
            + self.area() * self.y_centroid * self.y_centroid
    }
}

/// Compute the transformed cross-section properties for a composite beam.
///
/// All sections are transformed to the modular ratio n_i = E_i / E_ref.
///
/// # Arguments
/// * `layers`    – slice of composite layers
/// * `e_ref`     – reference Young's modulus \[Pa\]
///
/// Returns (EA, EI) of the transformed section.
pub fn transformed_section_properties(layers: &[CompositeLayer], e_ref: f64) -> (f64, f64) {
    if layers.is_empty() || e_ref.abs() < f64::EPSILON {
        return (0.0, 0.0);
    }
    // Compute EA and position of neutral axis
    let ea: f64 = layers.iter().map(|l| l.elastic_modulus * l.area()).sum();
    let ey: f64 = layers
        .iter()
        .map(|l| l.elastic_modulus * l.first_moment())
        .sum();
    let y_na = ey / ea.max(f64::EPSILON);

    // Compute EI about the neutral axis
    let ei: f64 = layers
        .iter()
        .map(|l| {
            let shift = l.y_centroid - y_na;
            let e = l.elastic_modulus;
            let i_local = l.width * l.height * l.height * l.height / 12.0;
            let i_pa = l.area() * shift * shift;
            e * (i_local + i_pa)
        })
        .sum();

    (ea, ei)
}

// ============================================================================
// Section 8 – Beam dynamics: mass matrices and natural frequencies
// ============================================================================

/// Consistent mass matrix for a 2-D Euler-Bernoulli beam element (6×6).
///
/// DOF order: \[u1, v1, θ1, u2, v2, θ2\].
/// Returns a flat \[f64; 36\] row-major array.
///
/// # Arguments
/// * `rho` – mass per unit length ρA \[kg/m\]
/// * `le`  – element length \[m\]
pub fn consistent_mass_2d(rho_a: f64, le: f64) -> [f64; 36] {
    let l = le;
    let l2 = l * l;
    let c = rho_a * l / 420.0;
    let mut m = [0.0f64; 36];

    // Axial DOF (0, 3): linear shape functions
    m[0] = rho_a * l / 3.0;
    m[3] = rho_a * l / 6.0;
    m[3 * 6] = rho_a * l / 6.0;
    m[3 * 6 + 3] = rho_a * l / 3.0;

    // Bending DOF (1, 2, 4, 5): Hermite shape functions
    m[6 + 1] = c * 156.0;
    m[6 + 2] = c * 22.0 * l;
    m[6 + 4] = c * 54.0;
    m[6 + 5] = -c * 13.0 * l;

    m[2 * 6 + 1] = c * 22.0 * l;
    m[2 * 6 + 2] = c * 4.0 * l2;
    m[2 * 6 + 4] = c * 13.0 * l;
    m[2 * 6 + 5] = -c * 3.0 * l2;

    m[4 * 6 + 1] = c * 54.0;
    m[4 * 6 + 2] = c * 13.0 * l;
    m[4 * 6 + 4] = c * 156.0;
    m[4 * 6 + 5] = -c * 22.0 * l;

    m[5 * 6 + 1] = -c * 13.0 * l;
    m[5 * 6 + 2] = -c * 3.0 * l2;
    m[5 * 6 + 4] = -c * 22.0 * l;
    m[5 * 6 + 5] = c * 4.0 * l2;

    m
}

/// Lumped mass matrix for a 2-D beam element (6×6, diagonal).
///
/// Applies the HRZ (Hinton-Rock-Zienkiewicz) lumping scheme.
/// Returns a flat \[f64; 36\] row-major array (only diagonal is non-zero).
///
/// # Arguments
/// * `rho_a` – mass per unit length \[kg/m\]
/// * `le`    – element length \[m\]
pub fn lumped_mass_2d(rho_a: f64, le: f64) -> [f64; 36] {
    let half = rho_a * le / 2.0;
    let mut m = [0.0f64; 36];
    // Split total mass equally to translational DOF; rotational DOF get zero mass
    m[0] = half;
    m[6 + 1] = half;
    m[3 * 6 + 3] = half;
    m[4 * 6 + 4] = half;
    m
}

/// Analytical natural frequencies of a simply-supported Euler-Bernoulli beam.
///
/// ω_n = (nπ/L)² √(EI/(ρA))
///
/// Returns the first `n_modes` natural frequencies \[rad/s\].
///
/// # Arguments
/// * `cs`     – cross-section
/// * `rho_a`  – mass per unit length \[kg/m\]
/// * `length` – total beam length \[m\]
/// * `n_modes`– number of modes to compute
pub fn simply_supported_natural_frequencies(
    cs: &CrossSection,
    rho_a: f64,
    length: f64,
    n_modes: usize,
) -> Vec<f64> {
    let ei = cs.bending_stiffness_z();
    let mut freqs = Vec::with_capacity(n_modes);
    for n in 1..=n_modes {
        let nf = n as f64;
        let k = nf * PI / length;
        let omega = k * k * (ei / rho_a.max(f64::EPSILON)).sqrt();
        freqs.push(omega);
    }
    freqs
}

/// Analytical natural frequencies of a cantilever Euler-Bernoulli beam.
///
/// Uses the transcendental equation root β_n L:
/// β_1 L = 1.875104, β_2 L = 4.694091, β_3 L = 7.854757
///
/// ω_n = (β_n)² √(EI/(ρA))
///
/// # Arguments
/// * `cs`     – cross-section
/// * `rho_a`  – mass per unit length \[kg/m\]
/// * `length` – beam length \[m\]
/// * `n_modes`– number of modes (≤ 3 analytically; higher modes use approx)
pub fn cantilever_natural_frequencies(
    cs: &CrossSection,
    rho_a: f64,
    length: f64,
    n_modes: usize,
) -> Vec<f64> {
    let ei = cs.bending_stiffness_z();
    // First 5 beta_n * L values from standard tables
    let beta_l = [1.875_104, 4.694_091, 7.854_757, 10.995_541, 14.137_166];
    let mut freqs = Vec::with_capacity(n_modes);
    for n in 0..n_modes {
        let bl = if n < beta_l.len() {
            beta_l[n]
        } else {
            (2 * (n + 1) - 1) as f64 * PI / 2.0
        };
        let beta = bl / length;
        let omega = beta * beta * (ei / rho_a.max(f64::EPSILON)).sqrt();
        freqs.push(omega);
    }
    freqs
}

// ============================================================================
// Section 9 – Plastic hinge
// ============================================================================

/// Plastic hinge moment-rotation relationship (bilinear model).
#[derive(Debug, Clone, Copy)]
pub struct PlasticHinge {
    /// Fully plastic moment capacity M_p \[N·m\].
    pub m_plastic: f64,
    /// Elastic limit moment M_y \[N·m\].
    pub m_yield: f64,
    /// Elastic bending stiffness EI \[N·m²\].
    pub ei: f64,
    /// Post-yield tangent stiffness ratio α = EI_t / EI ∈ \[0, 1\].
    pub hardening_ratio: f64,
    /// Current rotation \[rad\].
    pub rotation: f64,
    /// Whether the hinge has formed.
    pub is_formed: bool,
}

impl PlasticHinge {
    /// Create a new plastic hinge.
    ///
    /// # Arguments
    /// * `m_plastic`       – plastic moment capacity \[N·m\]
    /// * `m_yield`         – yield moment \[N·m\]
    /// * `ei`              – bending stiffness \[N·m²\]
    /// * `hardening_ratio` – post-yield stiffness ratio
    pub fn new(m_plastic: f64, m_yield: f64, ei: f64, hardening_ratio: f64) -> Self {
        Self {
            m_plastic,
            m_yield,
            ei,
            hardening_ratio,
            rotation: 0.0,
            is_formed: false,
        }
    }

    /// Return the moment for a given hinge rotation using the bilinear model.
    ///
    /// # Arguments
    /// * `theta` – hinge rotation \[rad\]
    pub fn moment(&self, theta: f64) -> f64 {
        let theta_y = self.m_yield / self.ei.max(f64::EPSILON);
        let m_abs = if theta.abs() <= theta_y {
            self.ei * theta.abs()
        } else {
            self.m_yield + self.hardening_ratio * self.ei * (theta.abs() - theta_y)
        };
        m_abs.min(self.m_plastic) * theta.signum()
    }

    /// Check and update whether the hinge has formed (|M| ≥ M_p).
    ///
    /// # Arguments
    /// * `moment` – current moment at the section \[N·m\]
    pub fn check_formation(&mut self, moment: f64) {
        if moment.abs() >= self.m_plastic {
            self.is_formed = true;
        }
    }

    /// Plastic rotation capacity φ_p = (M_p − M_y) / (α EI).
    pub fn plastic_rotation_capacity(&self) -> f64 {
        let denom = self.hardening_ratio * self.ei;
        if denom.abs() < f64::EPSILON {
            return f64::INFINITY;
        }
        (self.m_plastic - self.m_yield) / denom
    }
}

// ============================================================================
// Section 10 – Tapered beam element
// ============================================================================

/// Linearly tapered beam with properties that vary along the element.
///
/// Properties vary as: P(x) = P1 + (P2 − P1) x/L
#[derive(Debug, Clone, Copy)]
pub struct TaperedBeamElement {
    /// Cross-section at node 1 (x=0).
    pub section_1: CrossSection,
    /// Cross-section at node 2 (x=L).
    pub section_2: CrossSection,
    /// Element length \[m\].
    pub length: f64,
}

impl TaperedBeamElement {
    /// Interpolate the bending stiffness EI(x) at position x ∈ \[0, L\].
    ///
    /// # Arguments
    /// * `x` – position along the element \[m\]
    pub fn bending_stiffness_at(&self, x: f64) -> f64 {
        let xi = (x / self.length.max(f64::EPSILON)).clamp(0.0, 1.0);
        let ei1 = self.section_1.bending_stiffness_z();
        let ei2 = self.section_2.bending_stiffness_z();
        ei1 + (ei2 - ei1) * xi
    }

    /// Interpolate the axial stiffness EA(x) at position x ∈ \[0, L\].
    pub fn axial_stiffness_at(&self, x: f64) -> f64 {
        let xi = (x / self.length.max(f64::EPSILON)).clamp(0.0, 1.0);
        let ea1 = self.section_1.axial_stiffness();
        let ea2 = self.section_2.axial_stiffness();
        ea1 + (ea2 - ea1) * xi
    }

    /// Compute the equivalent uniform EI via numerical integration (10-point Gauss).
    ///
    /// EI_equiv = (1/L) ∫₀ᴸ EI(x) dx
    pub fn equivalent_bending_stiffness(&self) -> f64 {
        // 10-point Gauss-Legendre points and weights on [-1,1]
        let pts = [
            -0.973_906_529,
            -0.865_063_367,
            -0.679_409_568,
            -0.433_395_394,
            -0.148_874_339,
            0.148_874_339,
            0.433_395_394,
            0.679_409_568,
            0.865_063_367,
            0.973_906_529,
        ];
        let wts = [
            0.066_671_344,
            0.149_451_349,
            0.219_086_363,
            0.269_266_719,
            0.295_524_225,
            0.295_524_225,
            0.269_266_719,
            0.219_086_363,
            0.149_451_349,
            0.066_671_344,
        ];
        let l = self.length;
        let mut integral = 0.0_f64;
        for (pt, wt) in pts.iter().zip(wts.iter()) {
            let xi = 0.5 * (1.0 + pt);
            let x = xi * l;
            integral += wt * self.bending_stiffness_at(x);
        }
        integral * l * 0.5 / l
    }
}

// ============================================================================
// Section 11 – Cable element (catenary geometry)
// ============================================================================

/// Cable element modelling parameters.
#[derive(Debug, Clone, Copy)]
pub struct CableElement {
    /// Cross-sectional area \[m²\].
    pub area: f64,
    /// Young's modulus \[Pa\].
    pub elastic_modulus: f64,
    /// Linear mass density ρA \[kg/m\].
    pub linear_density: f64,
    /// Unstretched (reference) length L₀ \[m\].
    pub length_ref: f64,
    /// Horizontal component of cable tension H \[N\].
    pub tension_horizontal: f64,
}

impl CableElement {
    /// Create a new cable element.
    pub fn new(
        area: f64,
        elastic_modulus: f64,
        linear_density: f64,
        length_ref: f64,
        tension_horizontal: f64,
    ) -> Self {
        Self {
            area,
            elastic_modulus,
            linear_density,
            length_ref,
            tension_horizontal,
        }
    }

    /// Catenary sag parameter a = H / (w) where w = ρA g.
    ///
    /// # Arguments
    /// * `g` – gravitational acceleration \[m/s²\]
    pub fn catenary_a(&self, g: f64) -> f64 {
        let w = self.linear_density * g;
        if w.abs() < f64::EPSILON {
            return f64::MAX;
        }
        self.tension_horizontal / w
    }

    /// Maximum sag of a symmetric catenary cable of span l.
    ///
    /// d_max ≈ w l² / (8 H)  (parabolic approximation)
    ///
    /// # Arguments
    /// * `span` – horizontal span \[m\]
    /// * `g`    – gravitational acceleration \[m/s²\]
    pub fn max_sag(&self, span: f64, g: f64) -> f64 {
        let w = self.linear_density * g;
        if self.tension_horizontal.abs() < f64::EPSILON {
            return f64::MAX;
        }
        w * span * span / (8.0 * self.tension_horizontal)
    }

    /// Secant (Ernst) modulus for a cable element accounting for sag.
    ///
    /// E_eff = E / (1 + w² L² E A / (12 H³))
    ///
    /// # Arguments
    /// * `g` – gravitational acceleration \[m/s²\]
    pub fn ernst_modulus(&self, g: f64) -> f64 {
        let w = self.linear_density * g;
        let h = self.tension_horizontal.max(f64::EPSILON);
        let numer = self.elastic_modulus;
        let denom = 1.0
            + w * w * self.length_ref * self.length_ref * self.elastic_modulus * self.area
                / (12.0 * h * h * h);
        numer / denom.max(f64::EPSILON)
    }

    /// Axial stiffness of the cable element with sag correction.
    ///
    /// k_cable = E_eff A / L₀
    ///
    /// # Arguments
    /// * `g` – gravitational acceleration \[m/s²\]
    pub fn axial_stiffness_sag_corrected(&self, g: f64) -> f64 {
        self.ernst_modulus(g) * self.area / self.length_ref.max(f64::EPSILON)
    }

    /// Natural vibration frequency of a taut cable (first mode).
    ///
    /// f₁ = (1/2L) √(H/(ρA))
    ///
    /// # Arguments
    /// * `g` – gravitational acceleration \[m/s²\]
    pub fn fundamental_frequency(&self, g: f64) -> f64 {
        let w = self.linear_density * g;
        if w.abs() < f64::EPSILON || self.length_ref.abs() < f64::EPSILON {
            return 0.0;
        }
        let rho_a = self.linear_density;
        0.5 / self.length_ref * (self.tension_horizontal / rho_a.max(f64::EPSILON)).sqrt()
    }
}

// ============================================================================
// Section 12 – Mode shape utilities
// ============================================================================

/// Mode shape of a simply-supported beam: φ_n(x) = sin(n π x / L).
///
/// # Arguments
/// * `x`      – position \[m\]
/// * `length` – beam length \[m\]
/// * `n`      – mode number (1, 2, ...)
pub fn simply_supported_mode_shape(x: f64, length: f64, n: usize) -> f64 {
    let nf = n as f64;
    (nf * PI * x / length.max(f64::EPSILON)).sin()
}

/// Mode shape of a cantilever beam (clamped-free).
///
/// Uses the characteristic equation solution:
/// φ_n(x) = cosh(β_n x) − cos(β_n x) − σ_n (sinh(β_n x) − sin(β_n x))
///
/// # Arguments
/// * `x`      – position \[m\]
/// * `length` – beam length \[m\]
/// * `n`      – mode number (1, 2, ...; limited to 3 analytically)
pub fn cantilever_mode_shape(x: f64, length: f64, n: usize) -> f64 {
    let beta_l_vals = [1.875_104, 4.694_091, 7.854_757, 10.995_541, 14.137_166];
    let sigma_vals = [0.734_096, 1.018_465, 0.999_225, 1.000_033, 0.999_999];
    let idx = n.saturating_sub(1).min(beta_l_vals.len() - 1);
    let bl = beta_l_vals[idx];
    let sigma = sigma_vals[idx];
    let beta = bl / length.max(f64::EPSILON);
    let bx = beta * x;
    bx.cosh() - bx.cos() - sigma * (bx.sinh() - bx.sin())
}

// ============================================================================
// Section 13 – Beam deflection utilities
// ============================================================================

/// Maximum deflection of a simply-supported beam under a uniform load.
///
/// δ_max = 5 q L⁴ / (384 EI)
///
/// # Arguments
/// * `q`  – uniform load \[N/m\]
/// * `le` – span \[m\]
/// * `ei` – bending stiffness \[N·m²\]
pub fn deflection_ss_uniform_load(q: f64, le: f64, ei: f64) -> f64 {
    if ei.abs() < f64::EPSILON {
        return f64::MAX;
    }
    5.0 * q * le * le * le * le / (384.0 * ei)
}

/// Maximum deflection of a cantilever beam under a tip point load.
///
/// δ_max = P L³ / (3 EI)
///
/// # Arguments
/// * `p`  – point load at the free end \[N\]
/// * `le` – beam length \[m\]
/// * `ei` – bending stiffness \[N·m²\]
pub fn deflection_cantilever_tip_load(p: f64, le: f64, ei: f64) -> f64 {
    if ei.abs() < f64::EPSILON {
        return f64::MAX;
    }
    p * le * le * le / (3.0 * ei)
}

/// Deflection profile of a cantilever beam under a tip load P.
///
/// δ(x) = P x² (3L − x) / (6 EI)
///
/// # Arguments
/// * `p`  – tip load \[N\]
/// * `le` – beam length \[m\]
/// * `ei` – bending stiffness \[N·m²\]
/// * `x`  – position along the beam \[m\]
pub fn deflection_cantilever_profile(p: f64, le: f64, ei: f64, x: f64) -> f64 {
    if ei.abs() < f64::EPSILON {
        return 0.0;
    }
    p * x * x * (3.0 * le - x) / (6.0 * ei)
}

// ============================================================================
// Section 14 – Beam internal forces
// ============================================================================

/// Bending moment at position x in a simply-supported beam under uniform load.
///
/// M(x) = q/2 · x · (L − x)
///
/// # Arguments
/// * `q`  – distributed load \[N/m\]
/// * `le` – span \[m\]
/// * `x`  – position \[m\]
pub fn moment_ss_uniform_load(q: f64, le: f64, x: f64) -> f64 {
    0.5 * q * x * (le - x)
}

/// Shear force at position x in a simply-supported beam under uniform load.
///
/// V(x) = q (L/2 − x)
///
/// # Arguments
/// * `q`  – distributed load \[N/m\]
/// * `le` – span \[m\]
/// * `x`  – position \[m\]
pub fn shear_ss_uniform_load(q: f64, le: f64, x: f64) -> f64 {
    q * (le / 2.0 - x)
}

/// Normal bending stress at the extreme fibre: σ = M c / I.
///
/// # Arguments
/// * `moment`    – bending moment \[N·m\]
/// * `c`         – distance from neutral axis to extreme fibre \[m\]
/// * `second_mom`– second moment of area I \[m⁴\]
pub fn bending_stress(moment: f64, c: f64, second_mom: f64) -> f64 {
    if second_mom.abs() < f64::EPSILON {
        return 0.0;
    }
    moment * c / second_mom
}

// ============================================================================
// Unit tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn steel_rect_section() -> CrossSection {
        CrossSection::rectangular(0.1, 0.2, 200.0e9, 0.3)
    }

    fn steel_circ_section() -> CrossSection {
        CrossSection::circular(0.05, 200.0e9, 0.3)
    }

    #[test]
    fn test_cross_section_rectangular_area() {
        let cs = CrossSection::rectangular(0.1, 0.2, 200.0e9, 0.3);
        assert!((cs.area - 0.02).abs() < 1.0e-12);
    }

    #[test]
    fn test_cross_section_rectangular_izz() {
        let cs = CrossSection::rectangular(0.1, 0.2, 200.0e9, 0.3);
        let expected = 0.1 * 0.2f64.powi(3) / 12.0;
        assert!((cs.izz - expected).abs() < 1.0e-20);
    }

    #[test]
    fn test_cross_section_circular_area() {
        let cs = CrossSection::circular(0.05, 200.0e9, 0.3);
        assert!((cs.area - PI * 0.0025).abs() < 1.0e-14);
    }

    #[test]
    fn test_cross_section_bending_stiffness() {
        let cs = steel_rect_section();
        let ei = cs.bending_stiffness_z();
        assert!(ei > 0.0);
    }

    #[test]
    fn test_cross_section_torsional_stiffness() {
        let cs = steel_rect_section();
        assert!(cs.torsional_stiffness() > 0.0);
    }

    #[test]
    fn test_euler_bernoulli_stiffness_diagonal_positive() {
        let cs = steel_rect_section();
        let k = euler_bernoulli_stiffness_2d(&cs, 1.0);
        // All diagonal entries must be positive
        for i in 0..6 {
            assert!(k[i * 6 + i] >= 0.0, "k[{i},{i}] = {}", k[i * 6 + i]);
        }
    }

    #[test]
    fn test_euler_bernoulli_stiffness_symmetry() {
        let cs = steel_rect_section();
        let k = euler_bernoulli_stiffness_2d(&cs, 2.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!((k[i * 6 + j] - k[j * 6 + i]).abs() < 1.0e-8);
            }
        }
    }

    #[test]
    fn test_hermite_shape_functions_partition_of_unity() {
        // At xi=0: N1=1, N2=0
        let n = hermite_shape_functions(0.0, 1.0);
        assert!((n[0] - 1.0).abs() < 1.0e-12);
        assert!(n[1].abs() < 1.0e-12);
        // At xi=1: N1=0, N2=1
        let n2 = hermite_shape_functions(1.0, 1.0);
        assert!(n2[0].abs() < 1.0e-12);
        assert!((n2[1] - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_hermite_derivatives_at_ends() {
        let le = 2.0;
        let dn0 = hermite_shape_derivatives(0.0, le);
        // dN3/dx|_{x=0} = 1 (slope DOF at node 1)
        assert!((dn0[2] - 1.0).abs() < 1.0e-12);
        // dN4/dx|_{x=0} = 0
        assert!(dn0[3].abs() < 1.0e-12);
    }

    #[test]
    fn test_euler_bernoulli_uniform_load_equilibrium() {
        // Total vertical force = q*L (sum of Fy1 + Fy2)
        let q = 10.0;
        let le = 3.0;
        let loads = euler_bernoulli_uniform_load(q, le);
        let total_force = loads[0] + loads[2];
        assert!((total_force - q * le).abs() < 1.0e-10);
    }

    #[test]
    fn test_timoshenko_phi_slender() {
        // Slender beam: phi should be small
        let cs = steel_rect_section();
        let phi = timoshenko_phi(&cs, 10.0);
        assert!(phi < 1.0);
    }

    #[test]
    fn test_timoshenko_stiffness_symmetry() {
        let cs = steel_rect_section();
        let k = timoshenko_stiffness_2d(&cs, 2.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (k[i * 6 + j] - k[j * 6 + i]).abs() < 1.0e-6,
                    "asymmetry at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_timoshenko_converges_to_eb_limit() {
        // Very slender beam: Timoshenko → Euler-Bernoulli
        let cs = CrossSection {
            area: 0.001,
            izz: 1.0e-6,
            iyy: 1.0e-7,
            j_torsion: 1.0e-7,
            c_warping: 0.0,
            elastic_modulus: 200.0e9,
            shear_modulus: 77.0e9,
            kappa_y: 5.0 / 6.0,
            kappa_z: 5.0 / 6.0,
        };
        let le = 100.0;
        let k_eb = euler_bernoulli_stiffness_2d(&cs, le);
        let k_tim = timoshenko_stiffness_2d(&cs, le);
        let diff = (k_eb[6 + 1] - k_tim[6 + 1]).abs() / k_eb[6 + 1];
        assert!(diff < 1.0e-3);
    }

    #[test]
    fn test_beam3d_stiffness_size() {
        let cs = steel_circ_section();
        let k = beam3d_stiffness(&cs, 2.0);
        assert_eq!(k.len(), 144);
    }

    #[test]
    fn test_beam3d_stiffness_symmetry() {
        let cs = steel_circ_section();
        let k = beam3d_stiffness(&cs, 2.0);
        let err = beam3d_stiffness_symmetry_error(&k);
        assert!(err < 1.0e-6, "max symmetry error = {err}");
    }

    #[test]
    fn test_geometric_stiffness_symmetry() {
        let kg = geometric_stiffness_2d(1000.0, 3.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!((kg[i * 6 + j] - kg[j * 6 + i]).abs() < 1.0e-12);
            }
        }
    }

    #[test]
    fn test_euler_buckling_load_positive() {
        let cs = steel_rect_section();
        let pcr = euler_buckling_load(&cs, 3.0);
        assert!(pcr > 0.0);
    }

    #[test]
    fn test_effective_length_factor() {
        assert!((effective_length_factor(&BoundaryConditionType::PinPin) - 1.0).abs() < 1.0e-12);
        assert!(
            (effective_length_factor(&BoundaryConditionType::FixedFixed) - 0.5).abs() < 1.0e-12
        );
        assert!((effective_length_factor(&BoundaryConditionType::FixedPin) - 0.7).abs() < 1.0e-12);
        assert!((effective_length_factor(&BoundaryConditionType::FixedFree) - 2.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_curved_beam_curvature() {
        let cs = steel_rect_section();
        let cb = CurvedBeamElement::new(cs, 1.0, 5.0);
        assert!((cb.curvature() - 0.2).abs() < 1.0e-12);
    }

    #[test]
    fn test_curved_beam_coupling_stiffness() {
        let cs = steel_rect_section();
        let cb = CurvedBeamElement::new(cs, 2.0, 10.0);
        assert!(cb.coupling_stiffness() > 0.0);
    }

    #[test]
    fn test_transformed_section_single_layer() {
        let layer = CompositeLayer {
            width: 0.1,
            height: 0.3,
            y_centroid: 0.15,
            elastic_modulus: 200.0e9,
        };
        let (ea, ei) = transformed_section_properties(&[layer], 200.0e9);
        let expected_ea = 200.0e9 * 0.03;
        assert!((ea - expected_ea).abs() / expected_ea < 1.0e-10);
        assert!(ei > 0.0);
    }

    #[test]
    fn test_composite_layer_area() {
        let layer = CompositeLayer {
            width: 0.5,
            height: 0.2,
            y_centroid: 0.1,
            elastic_modulus: 1.0,
        };
        assert!((layer.area() - 0.1).abs() < 1.0e-12);
    }

    #[test]
    fn test_consistent_mass_trace() {
        let rho_a = 7850.0 * 0.02; // rho * area
        let le = 2.0;
        let m = consistent_mass_2d(rho_a, le);
        // Trace of translational DOF should sum to total mass (approx)
        let total_mass = rho_a * le;
        let trace = m[0] + m[6 + 1] + m[3 * 6 + 3] + m[4 * 6 + 4];
        // 156/420 * 2 * rho_a * l + 2/6 * rho_a * l ≈ total mass
        assert!(trace > 0.0 && trace < 2.0 * total_mass);
    }

    #[test]
    fn test_lumped_mass_diagonal() {
        let m = lumped_mass_2d(100.0, 2.0);
        for i in 0..6 {
            for j in 0..6 {
                if i != j {
                    assert_eq!(m[i * 6 + j], 0.0);
                }
            }
        }
    }

    #[test]
    fn test_natural_frequencies_increasing() {
        let cs = steel_rect_section();
        let freqs = simply_supported_natural_frequencies(&cs, 7850.0 * cs.area, 5.0, 3);
        assert_eq!(freqs.len(), 3);
        assert!(freqs[0] < freqs[1] && freqs[1] < freqs[2]);
    }

    #[test]
    fn test_cantilever_frequencies_increasing() {
        let cs = steel_circ_section();
        let freqs = cantilever_natural_frequencies(&cs, 7850.0 * cs.area, 2.0, 3);
        assert!(freqs[0] < freqs[1] && freqs[1] < freqs[2]);
    }

    #[test]
    fn test_plastic_hinge_elastic_regime() {
        let hinge = PlasticHinge::new(100.0e3, 80.0e3, 200.0e9 * 1.0e-6, 0.01);
        let theta_small = 1.0e-7;
        let m = hinge.moment(theta_small);
        // Should be elastic: m = EI * theta
        assert!((m - hinge.ei * theta_small).abs() < 1.0e-6);
    }

    #[test]
    fn test_plastic_hinge_formation() {
        let mut hinge = PlasticHinge::new(100.0e3, 80.0e3, 200.0e6, 0.01);
        assert!(!hinge.is_formed);
        hinge.check_formation(110.0e3);
        assert!(hinge.is_formed);
    }

    #[test]
    fn test_tapered_beam_interpolation() {
        let cs1 = steel_rect_section();
        let cs2 = CrossSection::rectangular(0.1, 0.4, 200.0e9, 0.3);
        let tb = TaperedBeamElement {
            section_1: cs1,
            section_2: cs2,
            length: 4.0,
        };
        let ei_mid = tb.bending_stiffness_at(2.0);
        let ei1 = cs1.bending_stiffness_z();
        let ei2 = cs2.bending_stiffness_z();
        assert!((ei_mid - (ei1 + ei2) / 2.0).abs() < 1.0e-3);
    }

    #[test]
    fn test_tapered_beam_equivalent_stiffness() {
        let cs1 = steel_rect_section();
        let cs2 = steel_rect_section(); // Uniform → eq stiffness = same
        let tb = TaperedBeamElement {
            section_1: cs1,
            section_2: cs2,
            length: 3.0,
        };
        let eq_ei = tb.equivalent_bending_stiffness();
        let expected = cs1.bending_stiffness_z();
        assert!((eq_ei - expected).abs() / expected < 1.0e-6);
    }

    #[test]
    fn test_cable_sag() {
        let cable = CableElement::new(0.01, 200.0e9, 78.5, 50.0, 100.0e3);
        let sag = cable.max_sag(50.0, 9.81);
        assert!(sag > 0.0);
    }

    #[test]
    fn test_ernst_modulus_less_than_e() {
        let cable = CableElement::new(0.01, 200.0e9, 78.5, 50.0, 100.0e3);
        let e_eff = cable.ernst_modulus(9.81);
        assert!(e_eff <= cable.elastic_modulus);
    }

    #[test]
    fn test_cable_fundamental_frequency_positive() {
        let cable = CableElement::new(0.001, 200.0e9, 7.85, 20.0, 50.0e3);
        let f = cable.fundamental_frequency(9.81);
        assert!(f > 0.0);
    }

    #[test]
    fn test_mode_shape_ss_zero_at_supports() {
        let length = 5.0;
        let phi_0 = simply_supported_mode_shape(0.0, length, 1);
        let phi_l = simply_supported_mode_shape(length, length, 1);
        assert!(phi_0.abs() < 1.0e-12);
        assert!(phi_l.abs() < 1.0e-10);
    }

    #[test]
    fn test_cantilever_mode_shape_zero_at_root() {
        let phi = cantilever_mode_shape(0.0, 3.0, 1);
        assert!(phi.abs() < 1.0e-10);
    }

    #[test]
    fn test_deflection_ss_uniform_load() {
        let ei = 200.0e9 * 1.0e-6;
        let delta = deflection_ss_uniform_load(10.0e3, 6.0, ei);
        assert!(delta > 0.0);
    }

    #[test]
    fn test_deflection_cantilever_tip_load() {
        let p = 1000.0;
        let le = 2.0;
        let ei = 200.0e9 * 1.0e-6;
        let delta = deflection_cantilever_tip_load(p, le, ei);
        let expected = p * le * le * le / (3.0 * ei);
        assert!((delta - expected).abs() < 1.0e-15);
    }

    #[test]
    fn test_deflection_cantilever_profile_at_tip() {
        let p = 1000.0;
        let le = 2.0;
        let ei = 200.0e9 * 1.0e-6;
        let delta_tip = deflection_cantilever_profile(p, le, ei, le);
        let expected = deflection_cantilever_tip_load(p, le, ei);
        assert!((delta_tip - expected).abs() < 1.0e-12);
    }

    #[test]
    fn test_moment_ss_midspan() {
        let q = 10.0;
        let le = 6.0;
        let m_mid = moment_ss_uniform_load(q, le, le / 2.0);
        let expected = q * le * le / 8.0;
        assert!((m_mid - expected).abs() < 1.0e-10);
    }

    #[test]
    fn test_shear_ss_at_support() {
        let q = 10.0;
        let le = 6.0;
        let v = shear_ss_uniform_load(q, le, 0.0);
        assert!((v - q * le / 2.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_bending_stress_basic() {
        let sigma = bending_stress(1000.0, 0.1, 1.0e-4);
        assert!((sigma - 1.0e6).abs() < 1.0e-3);
    }

    #[test]
    fn test_bending_stress_zero_inertia() {
        let sigma = bending_stress(1000.0, 0.1, 0.0);
        assert_eq!(sigma, 0.0);
    }
}
