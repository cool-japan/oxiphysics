// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shell finite element formulations.
//!
//! This module provides:
//!
//! - **Kirchhoff shell element**: thin shell, curvature-based bending
//! - **Mindlin-Reissner shell**: thick shell with transverse shear DOF
//! - **Membrane-bending coupling**: eccentricity-based coupling matrix
//! - **Drilling degrees of freedom**: in-plane rotation stabilisation
//! - **Shear locking remedies**: reduced integration, MITC formulation
//! - **Shell normals and thickness interpolation**
//! - **Composite shell laminate (CLT)**: \[A\], \[B\], \[D\] matrices
//! - **Shell buckling analysis**: geometric stiffness, critical load
//! - **Thermal shell analysis**: thermal load vector
//! - **Shell contact**: gap function, penalty stiffness

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Constants
// ============================================================================

/// Default shear correction factor κ for Mindlin-Reissner plates (5/6).
pub const SHEAR_CORRECTION_FACTOR: f64 = 5.0 / 6.0;

/// Number of DOF per node for a Mindlin-Reissner shell (6: u,v,w,θx,θy,θz).
pub const MINDLIN_DOF_PER_NODE: usize = 6;

/// Number of DOF per node for a Kirchhoff thin shell (5: u,v,w,θx,θy).
pub const KIRCHHOFF_DOF_PER_NODE: usize = 5;

// ============================================================================
// 1. Shell Element Geometry
// ============================================================================

/// 3-D shell node: position, normal vector, and thickness.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ShellNode {
    /// Nodal coordinates (x, y, z) \[m\].
    pub coords: [f64; 3],
    /// Unit normal to the shell mid-surface at this node.
    pub normal: [f64; 3],
    /// Shell thickness at this node \[m\].
    pub thickness: f64,
}

impl ShellNode {
    /// Create a new shell node.
    pub fn new(coords: [f64; 3], normal: [f64; 3], thickness: f64) -> Self {
        let norm = vec3_norm(normal);
        let n = if norm > 1.0e-14 {
            vec3_scale(normal, 1.0 / norm)
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            coords,
            normal: n,
            thickness,
        }
    }

    /// Top surface position: mid + (t/2) * n.
    pub fn top_surface(&self) -> [f64; 3] {
        vec3_add(self.coords, vec3_scale(self.normal, 0.5 * self.thickness))
    }

    /// Bottom surface position: mid - (t/2) * n.
    pub fn bottom_surface(&self) -> [f64; 3] {
        vec3_sub(self.coords, vec3_scale(self.normal, 0.5 * self.thickness))
    }
}

/// Compute the unit normal to a flat triangle defined by three nodes.
///
/// Uses the cross product of two edge vectors.  Returns the unit normal.
#[allow(dead_code)]
pub fn triangle_normal(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> [f64; 3] {
    let e1 = vec3_sub(p2, p1);
    let e2 = vec3_sub(p3, p1);
    let n = vec3_cross(e1, e2);
    let norm = vec3_norm(n);
    if norm < 1.0e-14 {
        [0.0, 0.0, 1.0]
    } else {
        vec3_scale(n, 1.0 / norm)
    }
}

/// Compute the area of a triangle from three 3-D points.
#[allow(dead_code)]
pub fn triangle_area(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let e1 = vec3_sub(p2, p1);
    let e2 = vec3_sub(p3, p1);
    let cross = vec3_cross(e1, e2);
    0.5 * vec3_norm(cross)
}

// ============================================================================
// 2. Kirchhoff Thin Shell Element
// ============================================================================

/// Kirchhoff thin shell element properties.
///
/// Based on classical thin-shell theory (no transverse shear deformation).
/// Suitable for shells with t/L < 1/20.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct KirchhoffShell {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Shell thickness \[m\].
    pub thickness: f64,
    /// Mass density \[kg/m³\].
    pub density: f64,
}

impl KirchhoffShell {
    /// Create a new Kirchhoff shell element.
    pub fn new(young_modulus: f64, poisson_ratio: f64, thickness: f64, density: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            thickness,
            density,
        }
    }

    /// Flexural rigidity D = E t³ / (12 (1 - ν²)).
    pub fn flexural_rigidity(&self) -> f64 {
        let t = self.thickness;
        self.young_modulus * t * t * t / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }

    /// Membrane stiffness A = E t / (1 - ν²).
    pub fn membrane_stiffness(&self) -> f64 {
        self.young_modulus * self.thickness / (1.0 - self.poisson_ratio * self.poisson_ratio)
    }

    /// Bending constitutive matrix \[D_b\] (3×3).
    ///
    /// Relates bending moments {Mx, My, Mxy} to curvatures {κx, κy, κxy}.
    pub fn bending_matrix(&self) -> [[f64; 3]; 3] {
        let d = self.flexural_rigidity();
        let nu = self.poisson_ratio;
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }

    /// Membrane constitutive matrix \[D_m\] (3×3).
    ///
    /// Relates membrane forces {Nx, Ny, Nxy} to membrane strains {εx, εy, γxy}.
    pub fn membrane_matrix(&self) -> [[f64; 3]; 3] {
        let a = self.membrane_stiffness();
        let nu = self.poisson_ratio;
        [
            [a, a * nu, 0.0],
            [a * nu, a, 0.0],
            [0.0, 0.0, a * (1.0 - nu) / 2.0],
        ]
    }

    /// Surface mass per unit area: ρ t.
    pub fn surface_mass(&self) -> f64 {
        self.density * self.thickness
    }

    /// Critical buckling pressure for a simply-supported rectangular plate a×b.
    ///
    /// Approximation using the lowest buckling mode (m=n=1):
    /// p_cr = D * π² / a² * (1 + (a/b)²)².
    pub fn critical_buckling_pressure(&self, a: f64, b: f64) -> f64 {
        let d = self.flexural_rigidity();
        let ratio = a / b;
        d * PI * PI / (a * a) * (1.0 + ratio * ratio).powi(2)
    }
}

// ============================================================================
// 3. Mindlin-Reissner Shell Element
// ============================================================================

/// Mindlin-Reissner shell element (thick shell with transverse shear).
///
/// Includes 5 DOF per node: (u, v, w, θx, θy).
/// Transverse shear deformation is included using a shear correction factor.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MindlinShell {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Shell thickness \[m\].
    pub thickness: f64,
    /// Mass density \[kg/m³\].
    pub density: f64,
    /// Shear correction factor κ (default 5/6).
    pub kappa: f64,
}

impl MindlinShell {
    /// Create a new Mindlin-Reissner shell element.
    pub fn new(young_modulus: f64, poisson_ratio: f64, thickness: f64, density: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            thickness,
            density,
            kappa: SHEAR_CORRECTION_FACTOR,
        }
    }

    /// Create with custom shear correction factor.
    pub fn with_kappa(mut self, kappa: f64) -> Self {
        self.kappa = kappa;
        self
    }

    /// Shear modulus: G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Transverse shear stiffness: Gs = κ G t.
    pub fn transverse_shear_stiffness(&self) -> f64 {
        self.kappa * self.shear_modulus() * self.thickness
    }

    /// Flexural rigidity D = E t³ / (12 (1 - ν²)).
    pub fn flexural_rigidity(&self) -> f64 {
        let t = self.thickness;
        self.young_modulus * t * t * t / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }

    /// Shear constitutive matrix \[D_s\] (2×2).
    ///
    /// Relates transverse shear forces {Qx, Qy} to shear strains {γxz, γyz}.
    pub fn shear_matrix(&self) -> [[f64; 2]; 2] {
        let gs = self.transverse_shear_stiffness();
        [[gs, 0.0], [0.0, gs]]
    }

    /// Bending constitutive matrix \[D_b\] (3×3).
    pub fn bending_matrix(&self) -> [[f64; 3]; 3] {
        let d = self.flexural_rigidity();
        let nu = self.poisson_ratio;
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }

    /// Rotary inertia per unit area: I_rot = ρ t³ / 12.
    pub fn rotary_inertia(&self) -> f64 {
        self.density * self.thickness.powi(3) / 12.0
    }

    /// Check shear locking risk: returns true if t/L > 1/10 (thick regime).
    pub fn is_thick(&self, element_size: f64) -> bool {
        self.thickness / element_size > 0.1
    }
}

// ============================================================================
// 4. Membrane-Bending Coupling
// ============================================================================

/// Coupling matrix \[B\] relating membrane forces and bending moments to
/// the mid-surface strains via eccentricity `e` (offset from reference surface).
///
/// For a shell with eccentricity `e`:
/// {N} = \[A\] {ε0} + \[B\] {κ}
/// {M} = \[B\] {ε0} + \[D\] {κ}
///
/// Returns the coupling matrix B (3×3) = e * \[A_m\].
#[allow(dead_code)]
pub fn membrane_bending_coupling(shell: &KirchhoffShell, eccentricity: f64) -> [[f64; 3]; 3] {
    let dm = shell.membrane_matrix();
    let e = eccentricity;
    let mut b = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] = e * dm[i][j];
        }
    }
    b
}

/// Compute the membrane strain from in-plane displacements (u, v) and
/// their gradients (∂u/∂x, ∂u/∂y, ∂v/∂x, ∂v/∂y).
///
/// Returns {εx, εy, γxy}.
#[allow(dead_code)]
pub fn membrane_strain(du_dx: f64, dv_dy: f64, du_dy: f64, dv_dx: f64) -> [f64; 3] {
    [du_dx, dv_dy, du_dy + dv_dx]
}

/// Compute bending curvatures from rotation gradients (Mindlin formulation).
///
/// κx = -∂θx/∂x, κy = -∂θy/∂y, κxy = -(∂θx/∂y + ∂θy/∂x).
#[allow(dead_code)]
pub fn bending_curvature(dthx_dx: f64, dthy_dy: f64, dthx_dy: f64, dthy_dx: f64) -> [f64; 3] {
    [-dthx_dx, -dthy_dy, -(dthx_dy + dthy_dx)]
}

// ============================================================================
// 5. Drilling Degrees of Freedom
// ============================================================================

/// Drilling rotation stabilisation stiffness for a flat shell element.
///
/// The drilling DOF θz is stabilised with a penalty term:
/// k_drill = α * G * t * A  where α is a small parameter (≈ 0.001).
///
/// `area` is the element area, `alpha` is the penalty parameter.
#[allow(dead_code)]
pub fn drilling_dof_stiffness(shell: &MindlinShell, area: f64, alpha: f64) -> f64 {
    alpha * shell.shear_modulus() * shell.thickness * area
}

/// Compute the in-plane rotation (drilling) strain from displacement gradients.
///
/// ω_z = (∂v/∂x - ∂u/∂y) / 2  (antisymmetric part of in-plane strain).
#[allow(dead_code)]
pub fn drilling_strain(du_dy: f64, dv_dx: f64) -> f64 {
    (dv_dx - du_dy) / 2.0
}

// ============================================================================
// 6. Shear Locking Remedies
// ============================================================================

/// Compute the assumed transverse shear strain at a 4-node MITC4 element.
///
/// MITC (Mixed Interpolation of Tensorial Components) evaluates shear strains
/// at tying points on the element edges.
///
/// `gamma_ab` = shear strain at the 4 tying points (A, B, C, D) on edges.
/// Returns the interpolated shear strain at parametric point (xi, eta).
#[allow(dead_code)]
pub fn mitc4_shear_strain(
    gamma_a: f64,
    gamma_b: f64,
    gamma_c: f64,
    gamma_d: f64,
    xi: f64,
    eta: f64,
) -> [f64; 2] {
    // MITC4 interpolation (simplified):
    // γxz = (1-eta)/2 * γA + (1+eta)/2 * γB
    // γyz = (1-xi)/2  * γC + (1+xi)/2  * γD
    let gamma_xz = 0.5 * (1.0 - eta) * gamma_a + 0.5 * (1.0 + eta) * gamma_b;
    let gamma_yz = 0.5 * (1.0 - xi) * gamma_c + 0.5 * (1.0 + xi) * gamma_d;
    [gamma_xz, gamma_yz]
}

/// Reduced integration weight for 1-point quadrature on a 2-D element.
///
/// Returns (xi, eta, weight) = (0, 0, 4) for the single Gauss point.
#[allow(dead_code)]
pub fn reduced_integration_point() -> (f64, f64, f64) {
    (0.0, 0.0, 4.0)
}

/// Full 2×2 Gauss quadrature points and weights on \[-1,1\]².
///
/// Returns array of (xi, eta, weight) tuples.
#[allow(dead_code)]
pub fn gauss_2x2() -> [(f64, f64, f64); 4] {
    let g = 1.0 / 3.0_f64.sqrt();
    [(-g, -g, 1.0), (g, -g, 1.0), (g, g, 1.0), (-g, g, 1.0)]
}

/// Compute the shear locking ratio for a Mindlin plate element.
///
/// The shear locking parameter β = (G * k * t² ) / (D * L²).
/// Values β >> 1 indicate severe shear locking risk.
#[allow(dead_code)]
pub fn shear_locking_parameter(shell: &MindlinShell, element_size: f64) -> f64 {
    let d = shell.flexural_rigidity();
    let gs = shell.transverse_shear_stiffness();
    gs * element_size * element_size / d
}

// ============================================================================
// 7. Composite Shell Laminate (Classical Lamination Theory — CLT)
// ============================================================================

/// A single laminate ply with material properties and orientation.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct Ply {
    /// Ply thickness \[m\].
    pub thickness: f64,
    /// Fibre orientation angle θ (radians from x-axis).
    pub angle_rad: f64,
    /// Longitudinal modulus E1 \[Pa\].
    pub e1: f64,
    /// Transverse modulus E2 \[Pa\].
    pub e2: f64,
    /// Shear modulus G12 \[Pa\].
    pub g12: f64,
    /// Major Poisson's ratio ν12.
    pub nu12: f64,
}

impl Ply {
    /// Create a new laminate ply.
    pub fn new(thickness: f64, angle_deg: f64, e1: f64, e2: f64, g12: f64, nu12: f64) -> Self {
        Self {
            thickness,
            angle_rad: angle_deg * PI / 180.0,
            e1,
            e2,
            g12,
            nu12,
        }
    }

    /// Minor Poisson's ratio: ν21 = ν12 * E2/E1.
    pub fn nu21(&self) -> f64 {
        self.nu12 * self.e2 / self.e1
    }

    /// On-axis reduced stiffness matrix Q (3×3) in material coordinates.
    pub fn q_matrix(&self) -> [[f64; 3]; 3] {
        let denom = 1.0 - self.nu12 * self.nu21();
        let q11 = self.e1 / denom;
        let q22 = self.e2 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q66 = self.g12;
        [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
    }

    /// Transformed reduced stiffness matrix Q̄ (3×3) at angle θ.
    pub fn q_bar_matrix(&self) -> [[f64; 3]; 3] {
        let q = self.q_matrix();
        let c = self.angle_rad.cos();
        let s = self.angle_rad.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let cs = c * s;
        let c2s2 = c2 * s2;
        [
            [
                q[0][0] * c4 + 2.0 * (q[0][1] + 2.0 * q[2][2]) * c2s2 + q[1][1] * s4,
                (q[0][0] + q[1][1] - 4.0 * q[2][2]) * c2s2 + q[0][1] * (c4 + s4),
                (q[0][0] - q[0][1] - 2.0 * q[2][2]) * cs * c2
                    - (q[1][1] - q[0][1] - 2.0 * q[2][2]) * cs * s2,
            ],
            [
                (q[0][0] + q[1][1] - 4.0 * q[2][2]) * c2s2 + q[0][1] * (c4 + s4),
                q[0][0] * s4 + 2.0 * (q[0][1] + 2.0 * q[2][2]) * c2s2 + q[1][1] * c4,
                (q[0][0] - q[0][1] - 2.0 * q[2][2]) * cs * s2
                    - (q[1][1] - q[0][1] - 2.0 * q[2][2]) * cs * c2,
            ],
            [
                (q[0][0] - q[0][1] - 2.0 * q[2][2]) * cs * c2
                    - (q[1][1] - q[0][1] - 2.0 * q[2][2]) * cs * s2,
                (q[0][0] - q[0][1] - 2.0 * q[2][2]) * cs * s2
                    - (q[1][1] - q[0][1] - 2.0 * q[2][2]) * cs * c2,
                (q[0][0] + q[1][1] - 2.0 * q[0][1] - 2.0 * q[2][2]) * c2s2 + q[2][2] * (c4 + s4),
            ],
        ]
    }
}

/// Laminate stacking: multiple plies forming a composite shell.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Laminate {
    /// Ordered list of plies from bottom (z_0) to top (z_n).
    pub plies: Vec<Ply>,
}

impl Laminate {
    /// Create a new laminate.
    pub fn new(plies: Vec<Ply>) -> Self {
        Self { plies }
    }

    /// Total thickness of the laminate.
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }

    /// Z-coordinates of ply interfaces (bottom = -t/2, top = +t/2).
    ///
    /// Returns n+1 values for n plies.
    pub fn z_coords(&self) -> Vec<f64> {
        let h = self.total_thickness();
        let mut z = Vec::with_capacity(self.plies.len() + 1);
        let mut zi = -0.5 * h;
        z.push(zi);
        for ply in &self.plies {
            zi += ply.thickness;
            z.push(zi);
        }
        z
    }

    /// CLT extensional stiffness matrix \[A\] (3×3).
    ///
    /// A_ij = Σ_k Q̄_ij^(k) * (z_{k+1} - z_k).
    pub fn a_matrix(&self) -> [[f64; 3]; 3] {
        let z = self.z_coords();
        let mut a = [[0.0f64; 3]; 3];
        for (k, ply) in self.plies.iter().enumerate() {
            let qb = ply.q_bar_matrix();
            let dz = z[k + 1] - z[k];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] += qb[i][j] * dz;
                }
            }
        }
        a
    }

    /// CLT coupling stiffness matrix \[B\] (3×3).
    ///
    /// B_ij = (1/2) Σ_k Q̄_ij^(k) * (z_{k+1}² - z_k²).
    pub fn b_matrix(&self) -> [[f64; 3]; 3] {
        let z = self.z_coords();
        let mut b = [[0.0f64; 3]; 3];
        for (k, ply) in self.plies.iter().enumerate() {
            let qb = ply.q_bar_matrix();
            let dz2 = z[k + 1] * z[k + 1] - z[k] * z[k];
            for i in 0..3 {
                for j in 0..3 {
                    b[i][j] += 0.5 * qb[i][j] * dz2;
                }
            }
        }
        b
    }

    /// CLT bending stiffness matrix \[D\] (3×3).
    ///
    /// D_ij = (1/3) Σ_k Q̄_ij^(k) * (z_{k+1}³ - z_k³).
    pub fn d_matrix(&self) -> [[f64; 3]; 3] {
        let z = self.z_coords();
        let mut d = [[0.0f64; 3]; 3];
        for (k, ply) in self.plies.iter().enumerate() {
            let qb = ply.q_bar_matrix();
            let dz3 = z[k + 1].powi(3) - z[k].powi(3);
            for i in 0..3 {
                for j in 0..3 {
                    d[i][j] += (1.0 / 3.0) * qb[i][j] * dz3;
                }
            }
        }
        d
    }

    /// Check if the laminate is symmetric (B ≈ 0).
    pub fn is_symmetric(&self) -> bool {
        let b = self.b_matrix();
        b.iter()
            .flat_map(|row| row.iter())
            .all(|&v| v.abs() < 1.0e-10)
    }

    /// Effective longitudinal modulus of the laminate: Ex = (A11 - A12²/A22) / h.
    pub fn effective_ex(&self) -> f64 {
        let a = self.a_matrix();
        let h = self.total_thickness();
        if h < 1.0e-15 || a[1][1].abs() < 1.0e-30 {
            return 0.0;
        }
        (a[0][0] - a[0][1] * a[0][1] / a[1][1]) / h
    }
}

// ============================================================================
// 8. Shell Buckling Analysis
// ============================================================================

/// Geometric stiffness matrix entry for membrane pre-stress.
///
/// For a shell element under in-plane stress resultants (Nx, Ny, Nxy),
/// the geometric stiffness integrand at a Gauss point is:
/// k_geo = Nx * (∂w/∂x)² + 2 Ny * (∂w/∂y)² + 2 Nxy * (∂w/∂x)(∂w/∂y)
///
/// Returns the scalar geometric stiffness contribution.
#[allow(dead_code)]
pub fn geometric_stiffness_integrand(nx: f64, ny: f64, nxy: f64, dw_dx: f64, dw_dy: f64) -> f64 {
    nx * dw_dx * dw_dx + ny * dw_dy * dw_dy + 2.0 * nxy * dw_dx * dw_dy
}

/// Estimate the critical buckling load factor for a simply-supported
/// flat rectangular plate under uniform in-plane compression Nx.
///
/// Uses the Euler buckling formula:
/// λ_cr = π² D / (Nx * a²) * (m + n² * a²/b²)²
/// for mode (m, n) = (1, 1).
#[allow(dead_code)]
pub fn plate_buckling_factor(shell: &KirchhoffShell, nx_applied: f64, a: f64, b: f64) -> f64 {
    let d = shell.flexural_rigidity();
    let k_cr = (1.0 + (a / b).powi(2)).powi(2);
    PI * PI * d * k_cr / (nx_applied * a * a)
}

/// Compute the buckling pressure for a cylindrical shell (Donnell theory, simplified).
///
/// p_cr = 2E / (sqrt(3(1-ν²))) * (t/R)² for long cylinders.
///
/// `radius` is the shell radius, `t` the thickness.
#[allow(dead_code)]
pub fn cylinder_buckling_pressure(e: f64, nu: f64, t: f64, radius: f64) -> f64 {
    2.0 * e / (3.0_f64 * (1.0 - nu * nu)).sqrt() * (t / radius).powi(2)
}

/// Compute eigenvalue (critical load multiplier) by inverse iteration.
///
/// Given stiffness `k` and geometric stiffness `kg` (both scalars for 1-D),
/// returns λ_cr = k / kg.
#[allow(dead_code)]
pub fn critical_load_multiplier(k: f64, kg: f64) -> f64 {
    if kg.abs() < 1.0e-30 {
        f64::INFINITY
    } else {
        k / kg
    }
}

// ============================================================================
// 9. Thermal Shell Analysis
// ============================================================================

/// Thermal load vector for a Kirchhoff shell element.
///
/// The thermal membrane resultant:
/// NT = A_11 * α * ΔT * t  (simplified isotropic).
///
/// The thermal bending resultant due to temperature gradient ΔT through thickness:
/// MT = D_11 * α * ΔT_z  where ΔT_z = ΔT / t.
///
/// Returns (NT, MT) pair (membrane and bending thermal loads).
#[allow(dead_code)]
pub fn thermal_load_resultants(
    shell: &KirchhoffShell,
    delta_t_mid: f64,
    delta_t_gradient: f64,
    alpha_thermal: f64,
) -> (f64, f64) {
    let a = shell.membrane_stiffness();
    let d = shell.flexural_rigidity();
    let nu = shell.poisson_ratio;
    let nt = a * (1.0 + nu) * alpha_thermal * delta_t_mid;
    let mt = d * (1.0 + nu) * alpha_thermal * delta_t_gradient / shell.thickness;
    (nt, mt)
}

/// Compute the thermal strain in a shell ply at temperature change ΔT.
///
/// Returns {αΔT, αΔT, 0} for an isotropic material.
#[allow(dead_code)]
pub fn thermal_strain_isotropic(alpha: f64, delta_t: f64) -> [f64; 3] {
    [alpha * delta_t, alpha * delta_t, 0.0]
}

/// Compute the thermal stress resultant vector {Nth_x, Nth_y, 0} for a
/// laminate with uniform temperature change ΔT.
///
/// N_th = \[A\] {ε_th}.
#[allow(dead_code)]
pub fn laminate_thermal_force(laminate: &Laminate, alpha: f64, delta_t: f64) -> [f64; 3] {
    let a = laminate.a_matrix();
    let eps_th = thermal_strain_isotropic(alpha, delta_t);
    let mut n = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            n[i] += a[i][j] * eps_th[j];
        }
    }
    n
}

// ============================================================================
// 10. Shell Contact
// ============================================================================

/// Shell contact gap function between two shell mid-surfaces.
///
/// Returns the signed gap distance g = (x2 - x1) · n_contact.
/// Positive gap means no contact; negative or zero means penetration.
#[allow(dead_code)]
pub fn shell_contact_gap(x1: [f64; 3], x2: [f64; 3], n_contact: [f64; 3]) -> f64 {
    let diff = vec3_sub(x2, x1);
    vec3_dot(diff, n_contact)
}

/// Penalty contact force for shell contact.
///
/// f_contact = -k_penalty * min(g, 0) * n_contact.
/// If gap g > 0, no contact force; if g < 0, apply penalty.
#[allow(dead_code)]
pub fn shell_contact_force(gap: f64, k_penalty: f64, n_contact: [f64; 3]) -> [f64; 3] {
    if gap >= 0.0 {
        [0.0; 3]
    } else {
        vec3_scale(n_contact, -k_penalty * gap)
    }
}

/// Hertz contact stiffness estimate for two curved shells.
///
/// k_Hertz = (4/3) * E* * sqrt(R_eff)
/// where E* = E / (2(1-ν²)) is the combined modulus and
/// R_eff = R1*R2/(R1+R2) is the effective radius.
#[allow(dead_code)]
pub fn hertz_contact_stiffness(e: f64, nu: f64, r1: f64, r2: f64) -> f64 {
    let e_star = e / (2.0 * (1.0 - nu * nu));
    let r_eff = if (r1 + r2).abs() < 1.0e-30 {
        return 0.0;
    } else {
        r1 * r2 / (r1 + r2)
    };
    (4.0 / 3.0) * e_star * r_eff.sqrt()
}

/// Augmented Lagrangian contact update for shell contact.
///
/// Updates the Lagrange multiplier: λ_new = max(0, λ + k_penalty * gap).
#[allow(dead_code)]
pub fn augmented_lagrangian_update(lambda: f64, k_penalty: f64, gap: f64) -> f64 {
    (lambda + k_penalty * gap).max(0.0)
}

// ============================================================================
// Helper vector utilities (no nalgebra)
// ============================================================================

/// Add two 3-D vectors.
#[allow(dead_code)]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors: a - b.
#[allow(dead_code)]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector by a scalar.
#[allow(dead_code)]
pub fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-D vectors.
#[allow(dead_code)]
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-D vectors.
#[allow(dead_code)]
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-D vector.
#[allow(dead_code)]
pub fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Matrix-vector product: (3×3) matrix times (3) vector.
#[allow(dead_code)]
pub fn mat3x3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    let mut result = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i] += m[i][j] * v[j];
        }
    }
    result
}

/// Matrix product: (3×3) × (3×3).
#[allow(dead_code)]
pub fn mat3x3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Trace (sum of diagonal) of a 3×3 matrix.
#[allow(dead_code)]
pub fn mat3x3_trace(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}

/// Compute the Frobenius norm of a 3×3 matrix.
#[allow(dead_code)]
pub fn mat3x3_frobenius(m: [[f64; 3]; 3]) -> f64 {
    m.iter()
        .flat_map(|row| row.iter())
        .map(|&v| v * v)
        .sum::<f64>()
        .sqrt()
}

// ============================================================================
// 11. Additional Utilities
// ============================================================================

/// Interpolate shell thickness between two nodes using linear shape functions.
///
/// N1 = (1 - xi) / 2,  N2 = (1 + xi) / 2.
#[allow(dead_code)]
pub fn interpolate_thickness(t1: f64, t2: f64, xi: f64) -> f64 {
    0.5 * (1.0 - xi) * t1 + 0.5 * (1.0 + xi) * t2
}

/// Compute the shear strain energy in a Mindlin shell element.
///
/// U_s = (1/2) * κ * G * t * (γxz² + γyz²) * area.
#[allow(dead_code)]
pub fn shear_strain_energy(shell: &MindlinShell, gamma_xz: f64, gamma_yz: f64, area: f64) -> f64 {
    let gs = shell.transverse_shear_stiffness();
    0.5 * gs * (gamma_xz * gamma_xz + gamma_yz * gamma_yz) * area
}

/// Compute the bending strain energy for a Kirchhoff shell.
///
/// U_b = (1/2) {κ}^T \[D_b\] {κ} * area.
#[allow(dead_code)]
pub fn bending_strain_energy(shell: &KirchhoffShell, curvature: [f64; 3], area: f64) -> f64 {
    let db = shell.bending_matrix();
    let db_kappa = mat3x3_vec(db, curvature);
    let dot: f64 = curvature
        .iter()
        .zip(db_kappa.iter())
        .map(|(&k, &dk)| k * dk)
        .sum();
    0.5 * dot * area
}

/// Compute the membrane strain energy.
///
/// U_m = (1/2) {ε}^T \[D_m\] {ε} * area.
#[allow(dead_code)]
pub fn membrane_strain_energy(shell: &KirchhoffShell, strain: [f64; 3], area: f64) -> f64 {
    let dm = shell.membrane_matrix();
    let dm_eps = mat3x3_vec(dm, strain);
    let dot: f64 = strain
        .iter()
        .zip(dm_eps.iter())
        .map(|(&e, &de)| e * de)
        .sum();
    0.5 * dot * area
}

/// Compute the natural frequency of a simply-supported rectangular Kirchhoff plate.
///
/// ω_mn = π² * sqrt(D / (ρ * t)) * (m²/a² + n²/b²).
#[allow(dead_code)]
pub fn plate_natural_frequency(shell: &KirchhoffShell, a: f64, b: f64, m: u32, n: u32) -> f64 {
    let d = shell.flexural_rigidity();
    let rho_t = shell.density * shell.thickness;
    PI * PI * (d / rho_t).sqrt() * (m as f64 * m as f64 / (a * a) + n as f64 * n as f64 / (b * b))
}

/// Compute the out-of-plane deflection at the center of a simply-supported
/// uniformly-loaded plate (Navier series, single term approximation).
///
/// w_max ≈ 16 q0 a⁴ / (π⁶ D * (1 + (a/b)²)²) for square plate.
#[allow(dead_code)]
pub fn plate_center_deflection(shell: &KirchhoffShell, q0: f64, a: f64, b: f64) -> f64 {
    let d = shell.flexural_rigidity();
    let ratio = a / b;
    16.0 * q0 * a.powi(4) / (PI.powi(6) * d * (1.0 + ratio * ratio).powi(2))
}

/// Compute the effective in-plane stiffness of a symmetric laminate \[A11_eff\].
///
/// Returns the effective A11/h for use in membrane analysis.
#[allow(dead_code)]
pub fn effective_a11(laminate: &Laminate) -> f64 {
    let h = laminate.total_thickness();
    if h < 1.0e-14 {
        return 0.0;
    }
    laminate.a_matrix()[0][0] / h
}

/// Compute the thermal expansion coefficient of a symmetric laminate in x-direction.
///
/// α_eff_x ≈ Σ_k Q̄11^k * α_k * t_k / (A11 * h).
#[allow(dead_code)]
pub fn laminate_effective_alpha(laminate: &Laminate, alpha_plies: &[f64]) -> f64 {
    if laminate.plies.len() != alpha_plies.len() {
        return 0.0;
    }
    let a = laminate.a_matrix();
    let h = laminate.total_thickness();
    if a[0][0].abs() < 1.0e-30 || h < 1.0e-14 {
        return 0.0;
    }
    let numerator: f64 = laminate
        .plies
        .iter()
        .zip(alpha_plies.iter())
        .map(|(ply, &alpha)| ply.q_bar_matrix()[0][0] * alpha * ply.thickness)
        .sum();
    numerator / (a[0][0])
}

/// Compute the force resultant norm |{N}| = sqrt(Nx² + Ny² + Nxy²).
#[allow(dead_code)]
pub fn force_resultant_norm(n: [f64; 3]) -> f64 {
    (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Transform a stress resultant vector from global to ply local coordinates.
///
/// Uses rotation matrix T for angle θ.
#[allow(dead_code)]
pub fn transform_stress_to_ply(n_global: [f64; 3], angle_rad: f64) -> [f64; 3] {
    let c = angle_rad.cos();
    let s = angle_rad.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;
    let t = [[c2, s2, 2.0 * cs], [s2, c2, -2.0 * cs], [-cs, cs, c2 - s2]];
    mat3x3_vec(t, n_global)
}

/// Compute the failure index for a ply using the Tsai-Wu criterion (2-D simplified).
///
/// FI = (σ1/F1t)² + (σ2/F2t)² - (σ1*σ2)/(F1t*F2t) + (τ12/F6)²
/// Values FI ≥ 1 indicate failure.
#[allow(dead_code)]
pub fn tsai_wu_failure_index(
    sigma1: f64,
    sigma2: f64,
    tau12: f64,
    f1t: f64,
    f2t: f64,
    f6: f64,
) -> f64 {
    (sigma1 / f1t).powi(2) + (sigma2 / f2t).powi(2) - sigma1 * sigma2 / (f1t * f2t)
        + (tau12 / f6).powi(2)
}

// ============================================================================
// Unit Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-10;

    fn default_kirchhoff() -> KirchhoffShell {
        KirchhoffShell::new(200.0e9, 0.3, 0.01, 7850.0)
    }

    fn default_mindlin() -> MindlinShell {
        MindlinShell::new(200.0e9, 0.3, 0.01, 7850.0)
    }

    fn isotropic_ply(angle_deg: f64) -> Ply {
        Ply::new(0.001, angle_deg, 200.0e9, 200.0e9, 76.9e9, 0.3)
    }

    #[test]
    fn test_kirchhoff_flexural_rigidity() {
        let s = default_kirchhoff();
        let d_expected = 200.0e9 * 0.01_f64.powi(3) / (12.0 * (1.0 - 0.09));
        assert!((s.flexural_rigidity() - d_expected).abs() / d_expected < 1.0e-12);
    }

    #[test]
    fn test_kirchhoff_bending_matrix_symmetry() {
        let s = default_kirchhoff();
        let db = s.bending_matrix();
        assert!((db[0][1] - db[1][0]).abs() < EPS);
    }

    #[test]
    fn test_kirchhoff_membrane_matrix_symmetry() {
        let s = default_kirchhoff();
        let dm = s.membrane_matrix();
        assert!((dm[0][1] - dm[1][0]).abs() < EPS);
    }

    #[test]
    fn test_kirchhoff_critical_buckling_positive() {
        let s = default_kirchhoff();
        let p = s.critical_buckling_pressure(1.0, 1.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_mindlin_shear_modulus() {
        let s = default_mindlin();
        let g_expected = 200.0e9 / (2.0 * 1.3);
        assert!((s.shear_modulus() - g_expected).abs() / g_expected < 1.0e-12);
    }

    #[test]
    fn test_mindlin_transverse_shear_stiffness_positive() {
        let s = default_mindlin();
        assert!(s.transverse_shear_stiffness() > 0.0);
    }

    #[test]
    fn test_mindlin_is_thick() {
        let s = default_mindlin();
        // t=0.01, L=0.05: ratio = 0.2 > 0.1 => thick
        assert!(s.is_thick(0.05));
        // t=0.01, L=0.5: ratio = 0.02 < 0.1 => thin
        assert!(!s.is_thick(0.5));
    }

    #[test]
    fn test_mindlin_shear_locking_parameter() {
        let s = default_mindlin();
        let beta = shear_locking_parameter(&s, 0.1);
        assert!(beta > 0.0);
    }

    #[test]
    fn test_shell_node_normal_normalized() {
        let n = ShellNode::new([0.0, 0.0, 0.0], [3.0, 0.0, 4.0], 0.01);
        let norm = vec3_norm(n.normal);
        assert!((norm - 1.0).abs() < EPS);
    }

    #[test]
    fn test_shell_node_top_bottom_symmetry() {
        let node = ShellNode::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.02);
        let top = node.top_surface();
        let bot = node.bottom_surface();
        assert!((top[2] - 0.01).abs() < EPS);
        assert!((bot[2] + 0.01).abs() < EPS);
    }

    #[test]
    fn test_triangle_normal_flat_xy() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let n = triangle_normal(p1, p2, p3);
        // Normal should be [0, 0, 1]
        assert!((n[2] - 1.0).abs() < EPS);
        assert!(n[0].abs() < EPS);
        assert!(n[1].abs() < EPS);
    }

    #[test]
    fn test_triangle_area() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let area = triangle_area(p1, p2, p3);
        assert!((area - 0.5).abs() < EPS);
    }

    #[test]
    fn test_membrane_strain_pure_extension() {
        let eps = membrane_strain(0.001, 0.0, 0.0, 0.0);
        assert!((eps[0] - 0.001).abs() < EPS);
        assert!(eps[1].abs() < EPS);
        assert!(eps[2].abs() < EPS);
    }

    #[test]
    fn test_bending_curvature_zero() {
        let k = bending_curvature(0.0, 0.0, 0.0, 0.0);
        assert_eq!(k, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_drilling_dof_stiffness_positive() {
        let s = default_mindlin();
        let k_d = drilling_dof_stiffness(&s, 0.01, 0.001);
        assert!(k_d > 0.0);
    }

    #[test]
    fn test_drilling_strain_antisymmetric() {
        // du/dy = dv/dx => ω_z = 0
        let omega = drilling_strain(0.5, 0.5);
        assert_eq!(omega, 0.0);
    }

    #[test]
    fn test_mitc4_shear_strain_center() {
        // At xi=0, eta=0: average of the 4 tying point values
        let gs = mitc4_shear_strain(0.1, 0.2, 0.3, 0.4, 0.0, 0.0);
        assert!((gs[0] - 0.15).abs() < EPS);
        assert!((gs[1] - 0.35).abs() < EPS);
    }

    #[test]
    fn test_gauss_2x2_weight_sum() {
        let pts = gauss_2x2();
        let w_sum: f64 = pts.iter().map(|&(_, _, w)| w).sum();
        assert!((w_sum - 4.0).abs() < EPS);
    }

    #[test]
    fn test_ply_q_matrix_symmetry() {
        let ply = isotropic_ply(0.0);
        let q = ply.q_matrix();
        assert!((q[0][1] - q[1][0]).abs() < EPS);
    }

    #[test]
    fn test_laminate_total_thickness() {
        let plies = vec![isotropic_ply(0.0), isotropic_ply(90.0)];
        let lam = Laminate::new(plies);
        assert!((lam.total_thickness() - 0.002).abs() < EPS);
    }

    #[test]
    fn test_laminate_a_matrix_symmetry() {
        let plies = vec![isotropic_ply(0.0), isotropic_ply(45.0), isotropic_ply(90.0)];
        let lam = Laminate::new(plies);
        let a = lam.a_matrix();
        assert!((a[0][1] - a[1][0]).abs() < 1.0e-6);
    }

    #[test]
    fn test_laminate_b_matrix_symmetric_layup() {
        // Symmetric [0/90]s laminate should have B ≈ 0
        let plies = vec![
            isotropic_ply(0.0),
            isotropic_ply(90.0),
            isotropic_ply(90.0),
            isotropic_ply(0.0),
        ];
        let lam = Laminate::new(plies);
        assert!(lam.is_symmetric(), "symmetric laminate should have B ≈ 0");
    }

    #[test]
    fn test_laminate_d_matrix_positive_definite_diagonal() {
        let plies = vec![isotropic_ply(0.0), isotropic_ply(90.0)];
        let lam = Laminate::new(plies);
        let d = lam.d_matrix();
        assert!(d[0][0] > 0.0);
        assert!(d[1][1] > 0.0);
        assert!(d[2][2] > 0.0);
    }

    #[test]
    fn test_geometric_stiffness_integrand_uniaxial() {
        // Pure Nx compression, no Ny or Nxy
        let ks = geometric_stiffness_integrand(1.0, 0.0, 0.0, 1.0, 0.0);
        assert!((ks - 1.0).abs() < EPS);
    }

    #[test]
    fn test_plate_buckling_factor_positive() {
        let s = default_kirchhoff();
        let lambda = plate_buckling_factor(&s, 1.0e6, 1.0, 1.0);
        assert!(lambda > 0.0);
    }

    #[test]
    fn test_cylinder_buckling_pressure_positive() {
        let p = cylinder_buckling_pressure(200.0e9, 0.3, 0.01, 1.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_critical_load_multiplier_zero_geometric_stiffness() {
        let lambda = critical_load_multiplier(1000.0, 0.0);
        assert!(lambda.is_infinite());
    }

    #[test]
    fn test_thermal_load_resultants_nonzero() {
        let s = default_kirchhoff();
        let (nt, mt) = thermal_load_resultants(&s, 100.0, 0.0, 12.0e-6);
        assert!(nt.abs() > 0.0);
        assert!(mt.abs() < EPS); // zero gradient => zero bending moment
    }

    #[test]
    fn test_laminate_thermal_force_isotropic() {
        let plies = vec![isotropic_ply(0.0)];
        let lam = Laminate::new(plies);
        let n = laminate_thermal_force(&lam, 12.0e-6, 100.0);
        // For isotropic single ply, Nx = Ny, Nxy ≈ 0
        assert!((n[0] - n[1]).abs() / n[0].abs() < 1.0e-10);
    }

    #[test]
    fn test_shell_contact_gap_no_penetration() {
        let x1 = [0.0, 0.0, 0.0];
        let x2 = [0.0, 0.0, 1.0];
        let n = [0.0, 0.0, 1.0];
        let g = shell_contact_gap(x1, x2, n);
        assert!((g - 1.0).abs() < EPS);
    }

    #[test]
    fn test_shell_contact_force_no_contact() {
        let f = shell_contact_force(0.5, 1.0e6, [0.0, 0.0, 1.0]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_shell_contact_force_penetration() {
        let f = shell_contact_force(-0.001, 1.0e6, [0.0, 0.0, 1.0]);
        // Force should be positive (pushing apart)
        assert!(f[2] > 0.0);
    }

    #[test]
    fn test_hertz_contact_stiffness_positive() {
        let k = hertz_contact_stiffness(200.0e9, 0.3, 0.01, 0.01);
        assert!(k > 0.0);
    }

    #[test]
    fn test_plate_natural_frequency_positive() {
        let s = default_kirchhoff();
        let w = plate_natural_frequency(&s, 1.0, 1.0, 1, 1);
        assert!(w > 0.0);
    }

    #[test]
    fn test_tsai_wu_failure_below_one() {
        // Low stresses should give FI < 1
        let fi = tsai_wu_failure_index(100.0e6, 50.0e6, 20.0e6, 1500.0e6, 600.0e6, 80.0e6);
        assert!(fi < 1.0);
    }

    #[test]
    fn test_tsai_wu_failure_above_one() {
        // Stresses at strength limit: σ1 = F1t => FI ≈ 1
        let fi = tsai_wu_failure_index(1500.0e6, 0.0, 0.0, 1500.0e6, 600.0e6, 80.0e6);
        assert!((fi - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_bending_strain_energy_zero_curvature() {
        let s = default_kirchhoff();
        let u = bending_strain_energy(&s, [0.0, 0.0, 0.0], 1.0);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_membrane_strain_energy_positive() {
        let s = default_kirchhoff();
        let u = membrane_strain_energy(&s, [0.001, 0.0, 0.0], 0.01);
        assert!(u > 0.0);
    }

    #[test]
    fn test_vec3_cross_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[2] - 1.0).abs() < EPS);
        assert!(c[0].abs() < EPS);
        assert!(c[1].abs() < EPS);
    }

    #[test]
    fn test_interpolate_thickness_midpoint() {
        let t = interpolate_thickness(0.01, 0.02, 0.0);
        assert!((t - 0.015).abs() < EPS);
    }

    #[test]
    fn test_shear_strain_energy_zero() {
        let s = default_mindlin();
        let u = shear_strain_energy(&s, 0.0, 0.0, 0.01);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_augmented_lagrangian_no_penetration() {
        // gap > 0 (no penetration): lambda stays >= 0
        let lam = augmented_lagrangian_update(0.0, 1.0e6, 0.001);
        assert!(lam >= 0.0);
    }
}
