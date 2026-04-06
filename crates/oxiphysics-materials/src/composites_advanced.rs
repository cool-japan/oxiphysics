// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced composite material models: laminate theory, fiber-matrix, damage.
//!
//! Implements Classical Laminate Theory (CLT), Halpin-Tsai micromechanics,
//! Mori-Tanaka homogenization, failure criteria (Tsai-Wu, Puck, Max Stress),
//! progressive damage, and delamination models.

#![warn(missing_docs)]

/// Fiber material properties.
#[derive(Debug, Clone)]
pub struct FiberProperties {
    /// Axial (fiber direction) Young's modulus \[Pa\].
    pub e1: f64,
    /// Transverse Young's modulus \[Pa\].
    pub e2: f64,
    /// In-plane shear modulus \[Pa\].
    pub g12: f64,
    /// Major Poisson's ratio.
    pub nu12: f64,
    /// Fiber volume fraction (0-1).
    pub volume_fraction: f64,
    /// Fiber aspect ratio (length/diameter).
    pub aspect_ratio: f64,
}

impl FiberProperties {
    /// Create new fiber properties.
    pub fn new(
        e1: f64,
        e2: f64,
        g12: f64,
        nu12: f64,
        volume_fraction: f64,
        aspect_ratio: f64,
    ) -> Self {
        Self {
            e1,
            e2,
            g12,
            nu12,
            volume_fraction,
            aspect_ratio,
        }
    }

    /// Create carbon fiber (T300) properties.
    pub fn carbon_t300() -> Self {
        Self {
            e1: 230e9,
            e2: 15e9,
            g12: 15e9,
            nu12: 0.27,
            volume_fraction: 0.6,
            aspect_ratio: 1000.0,
        }
    }

    /// Create glass fiber (E-glass) properties.
    pub fn glass_e() -> Self {
        Self {
            e1: 72e9,
            e2: 72e9,
            g12: 30e9,
            nu12: 0.22,
            volume_fraction: 0.5,
            aspect_ratio: 500.0,
        }
    }
}

/// Matrix material properties.
#[derive(Debug, Clone)]
pub struct MatrixProperties {
    /// Young's modulus \[Pa\].
    pub modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Tensile yield stress \[Pa\].
    pub yield_stress: f64,
    /// Fracture toughness K_Ic \[Pa·√m\].
    pub fracture_toughness: f64,
}

impl MatrixProperties {
    /// Create new matrix properties.
    pub fn new(
        modulus: f64,
        poisson_ratio: f64,
        yield_stress: f64,
        fracture_toughness: f64,
    ) -> Self {
        Self {
            modulus,
            poisson_ratio,
            yield_stress,
            fracture_toughness,
        }
    }

    /// Create epoxy resin matrix properties.
    pub fn epoxy() -> Self {
        Self {
            modulus: 3.5e9,
            poisson_ratio: 0.35,
            yield_stress: 80e6,
            fracture_toughness: 0.5e6,
        }
    }

    /// Shear modulus G = E / (2*(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
}

/// Halpin-Tsai micromechanics model.
///
/// Predicts elastic constants of a composite from fiber + matrix properties.
#[derive(Debug, Clone)]
pub struct HalpinTsai {
    /// Fiber properties.
    pub fiber: FiberProperties,
    /// Matrix properties.
    pub matrix: MatrixProperties,
}

impl HalpinTsai {
    /// Create a new Halpin-Tsai model.
    pub fn new(fiber: FiberProperties, matrix: MatrixProperties) -> Self {
        Self { fiber, matrix }
    }

    fn vf(&self) -> f64 {
        self.fiber.volume_fraction
    }
    fn vm(&self) -> f64 {
        1.0 - self.fiber.volume_fraction
    }

    /// Longitudinal modulus E11 (rule of mixtures).
    pub fn e11(&self) -> f64 {
        self.fiber.e1 * self.vf() + self.matrix.modulus * self.vm()
    }

    /// Transverse modulus E22 using Halpin-Tsai.
    pub fn e22(&self) -> f64 {
        let xi = 2.0; // geometry parameter for circular fibers
        let ef = self.fiber.e2;
        let em = self.matrix.modulus;
        let eta = (ef / em - 1.0) / (ef / em + xi);
        em * (1.0 + xi * eta * self.vf()) / (1.0 - eta * self.vf())
    }

    /// In-plane shear modulus G12 using Halpin-Tsai.
    pub fn g12(&self) -> f64 {
        let xi = 1.0;
        let gf = self.fiber.g12;
        let gm = self.matrix.shear_modulus();
        let eta = (gf / gm - 1.0) / (gf / gm + xi);
        gm * (1.0 + xi * eta * self.vf()) / (1.0 - eta * self.vf())
    }

    /// Major Poisson's ratio ν12 (rule of mixtures).
    pub fn nu12(&self) -> f64 {
        self.fiber.nu12 * self.vf() + self.matrix.poisson_ratio * self.vm()
    }

    /// Minor Poisson's ratio ν21 = ν12 * E22 / E11.
    pub fn nu21(&self) -> f64 {
        self.nu12() * self.e22() / self.e11()
    }

    /// Voigt upper bound for E11.
    pub fn voigt_e11(&self) -> f64 {
        self.fiber.e1 * self.vf() + self.matrix.modulus * self.vm()
    }

    /// Reuss lower bound for E22.
    pub fn reuss_e22(&self) -> f64 {
        1.0 / (self.vf() / self.fiber.e2 + self.vm() / self.matrix.modulus)
    }
}

/// Mori-Tanaka homogenization for ellipsoidal inclusions.
#[derive(Debug, Clone)]
pub struct MoriTanakaComposite {
    /// Fiber (inclusion) properties.
    pub fiber: FiberProperties,
    /// Matrix properties.
    pub matrix: MatrixProperties,
}

impl MoriTanakaComposite {
    /// Create a new Mori-Tanaka model.
    pub fn new(fiber: FiberProperties, matrix: MatrixProperties) -> Self {
        Self { fiber, matrix }
    }

    /// Effective bulk modulus (simplified spherical inclusion).
    pub fn effective_bulk_modulus(&self) -> f64 {
        let km = self.matrix.modulus / (3.0 * (1.0 - 2.0 * self.matrix.poisson_ratio));
        let gm = self.matrix.shear_modulus();
        // Eshelby factor for spheres
        let kf = self.fiber.e1 / (3.0 * (1.0 - 2.0 * self.fiber.nu12));
        let vf = self.fiber.volume_fraction;
        let vm = 1.0 - vf;
        let alpha = 3.0 * km / (3.0 * km + 4.0 * gm);
        km + vf * (kf - km) / (1.0 + vm * alpha * (kf - km) / km)
    }

    /// Effective shear modulus (simplified).
    pub fn effective_shear_modulus(&self) -> f64 {
        let gm = self.matrix.shear_modulus();
        let gf = self.fiber.g12;
        let vf = self.fiber.volume_fraction;
        let vm = 1.0 - vf;
        let km = self.matrix.modulus / (3.0 * (1.0 - 2.0 * self.matrix.poisson_ratio));
        let beta = 6.0 * (km + 2.0 * gm) / (5.0 * (3.0 * km + 4.0 * gm));
        gm + vf * (gf - gm) / (1.0 + vm * beta * (gf - gm) / gm)
    }
}

/// Single lamina (ply) in a laminate stack.
#[derive(Debug, Clone)]
pub struct Lamina {
    /// Longitudinal modulus E1 \[Pa\].
    pub e1: f64,
    /// Transverse modulus E2 \[Pa\].
    pub e2: f64,
    /// In-plane shear modulus G12 \[Pa\].
    pub g12: f64,
    /// Major Poisson's ratio.
    pub nu12: f64,
    /// Ply thickness \[m\].
    pub thickness: f64,
    /// Fiber angle from laminate x-axis \[radians\].
    pub angle: f64,
}

impl Lamina {
    /// Create a new lamina.
    pub fn new(e1: f64, e2: f64, g12: f64, nu12: f64, thickness: f64, angle_deg: f64) -> Self {
        Self {
            e1,
            e2,
            g12,
            nu12,
            thickness,
            angle: angle_deg.to_radians(),
        }
    }

    /// Minor Poisson's ratio.
    pub fn nu21(&self) -> f64 {
        self.nu12 * self.e2 / self.e1
    }

    /// Reduced stiffness matrix Q in principal material coordinates \[Q11, Q12, Q22, Q66\].
    pub fn q_matrix(&self) -> [f64; 4] {
        let denom = 1.0 - self.nu12 * self.nu21();
        [
            self.e1 / denom,             // Q11
            self.nu12 * self.e2 / denom, // Q12
            self.e2 / denom,             // Q22
            self.g12,                    // Q66
        ]
    }

    /// Transformed reduced stiffness matrix Q_bar (in laminate coords).
    /// Returns \[Q11, Q12, Q16, Q22, Q26, Q66\].
    pub fn q_bar(&self) -> [f64; 6] {
        let q = self.q_matrix();
        let c = self.angle.cos();
        let s = self.angle.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let c2s2 = c2 * s2;
        [
            q[0] * c4 + 2.0 * (q[1] + 2.0 * q[3]) * c2s2 + q[2] * s4,
            (q[0] + q[2] - 4.0 * q[3]) * c2s2 + q[1] * (c4 + s4),
            (q[0] - q[1] - 2.0 * q[3]) * c * s * c2 - (q[2] - q[1] - 2.0 * q[3]) * c * s * s2,
            q[0] * s4 + 2.0 * (q[1] + 2.0 * q[3]) * c2s2 + q[2] * c4,
            (q[0] - q[1] - 2.0 * q[3]) * c * s * s2 - (q[2] - q[1] - 2.0 * q[3]) * c * s * c2,
            (q[0] + q[2] - 2.0 * q[1] - 2.0 * q[3]) * c2s2 + q[3] * (c4 + s4),
        ]
    }
}

/// Classical Laminate Theory stiffness matrices.
///
/// Computes \[A|B|D\] matrices for a laminate stack.
#[derive(Debug, Clone)]
pub struct ClassicalLaminateTheory {
    /// Laminae in the stack (bottom to top).
    pub laminae: Vec<Lamina>,
}

impl ClassicalLaminateTheory {
    /// Create CLT from a list of laminae.
    pub fn new(laminae: Vec<Lamina>) -> Self {
        Self { laminae }
    }

    /// Total thickness \[m\].
    pub fn total_thickness(&self) -> f64 {
        self.laminae.iter().map(|l| l.thickness).sum()
    }

    /// Z-coordinates of ply interfaces (bottom to top).
    pub fn z_coords(&self) -> Vec<f64> {
        let h = self.total_thickness();
        let mut z = vec![-h / 2.0];
        for l in &self.laminae {
            let top = *z.last().expect("collection should not be empty") + l.thickness;
            z.push(top);
        }
        z
    }

    /// A matrix (extensional stiffness), 3x3 as \[A11, A12, A16, A22, A26, A66\].
    pub fn a_matrix(&self) -> [f64; 6] {
        let z = self.z_coords();
        let mut a = [0.0f64; 6];
        for (i, l) in self.laminae.iter().enumerate() {
            let qb = l.q_bar();
            let dz = z[i + 1] - z[i];
            for j in 0..6 {
                a[j] += qb[j] * dz;
            }
        }
        a
    }

    /// B matrix (coupling stiffness) \[B11, B12, B16, B22, B26, B66\].
    pub fn b_matrix(&self) -> [f64; 6] {
        let z = self.z_coords();
        let mut b = [0.0f64; 6];
        for (i, l) in self.laminae.iter().enumerate() {
            let qb = l.q_bar();
            let dz2 = (z[i + 1].powi(2) - z[i].powi(2)) / 2.0;
            for j in 0..6 {
                b[j] += qb[j] * dz2;
            }
        }
        b
    }

    /// D matrix (bending stiffness) \[D11, D12, D16, D22, D26, D66\].
    pub fn d_matrix(&self) -> [f64; 6] {
        let z = self.z_coords();
        let mut d = [0.0f64; 6];
        for (i, l) in self.laminae.iter().enumerate() {
            let qb = l.q_bar();
            let dz3 = (z[i + 1].powi(3) - z[i].powi(3)) / 3.0;
            for j in 0..6 {
                d[j] += qb[j] * dz3;
            }
        }
        d
    }

    /// Check if laminate is symmetric (B ≈ 0).
    pub fn is_symmetric(&self) -> bool {
        let b = self.b_matrix();
        b.iter().all(|&v| v.abs() < 1e-6)
    }

    /// Check if laminate is balanced (A16 = A26 ≈ 0).
    pub fn is_balanced(&self) -> bool {
        let a = self.a_matrix();
        a[2].abs() < 1e-6 && a[4].abs() < 1e-6
    }
}

/// Laminate response under in-plane loads and moments.
#[derive(Debug, Clone)]
pub struct LaminateResponse {
    /// Applied in-plane force resultants \[N/m\]: \[Nx, Ny, Nxy\].
    pub n: [f64; 3],
    /// Applied moment resultants \[N\]: \[Mx, My, Mxy\].
    pub m: [f64; 3],
    /// Midplane strains \[ε_x0, ε_y0, γ_xy0\].
    pub midplane_strains: [f64; 3],
    /// Midplane curvatures \[κ_x, κ_y, κ_xy\].
    pub curvatures: [f64; 3],
}

impl LaminateResponse {
    /// Compute midplane strains for a symmetric laminate (B=0) under pure N.
    /// Uses a_inv (3x3 inverse A in Voigt notation order \[11,12,22,16,26,66\]).
    pub fn from_n_symmetric(n: [f64; 3], a: [f64; 6]) -> Self {
        // Simplified: compute midplane strains from [A]{ε} = {N}
        // For unidirectional symmetric: approximate 3x3 solve
        let a11 = a[0];
        let a12 = a[1];
        let a22 = a[3];
        let det = a11 * a22 - a12 * a12;
        let eps_x = (a22 * n[0] - a12 * n[1]) / det;
        let eps_y = (-a12 * n[0] + a11 * n[1]) / det;
        let gamma_xy = if a[5] > 0.0 { n[2] / a[5] } else { 0.0 };
        Self {
            n,
            m: [0.0; 3],
            midplane_strains: [eps_x, eps_y, gamma_xy],
            curvatures: [0.0; 3],
        }
    }
}

/// Tsai-Wu quadratic failure criterion.
#[derive(Debug, Clone)]
pub struct TsaiWuCriterion {
    /// Longitudinal tensile strength Xt \[Pa\].
    pub xt: f64,
    /// Longitudinal compressive strength Xc \[Pa\].
    pub xc: f64,
    /// Transverse tensile strength Yt \[Pa\].
    pub yt: f64,
    /// Transverse compressive strength Yc \[Pa\].
    pub yc: f64,
    /// In-plane shear strength S12 \[Pa\].
    pub s12: f64,
    /// Interaction term F12* (typically -0.5 ≤ F12* ≤ 0).
    pub f12_star: f64,
}

impl TsaiWuCriterion {
    /// Create Tsai-Wu criterion.
    pub fn new(xt: f64, xc: f64, yt: f64, yc: f64, s12: f64, f12_star: f64) -> Self {
        Self {
            xt,
            xc,
            yt,
            yc,
            s12,
            f12_star,
        }
    }

    /// F1 coefficient.
    pub fn f1(&self) -> f64 {
        1.0 / self.xt - 1.0 / self.xc
    }

    /// F2 coefficient.
    pub fn f2(&self) -> f64 {
        1.0 / self.yt - 1.0 / self.yc
    }

    /// F11 coefficient.
    pub fn f11(&self) -> f64 {
        1.0 / (self.xt * self.xc)
    }

    /// F22 coefficient.
    pub fn f22(&self) -> f64 {
        1.0 / (self.yt * self.yc)
    }

    /// F66 coefficient.
    pub fn f66(&self) -> f64 {
        1.0 / (self.s12 * self.s12)
    }

    /// F12 interaction coefficient.
    pub fn f12(&self) -> f64 {
        self.f12_star / (self.f11() * self.f22()).sqrt()
    }

    /// Tsai-Wu failure index (FI < 1 = safe).
    pub fn failure_index(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        self.f1() * sigma1
            + self.f2() * sigma2
            + self.f11() * sigma1 * sigma1
            + self.f22() * sigma2 * sigma2
            + self.f66() * tau12 * tau12
            + 2.0 * self.f12() * sigma1 * sigma2
    }

    /// Check if stress state fails (FI ≥ 1).
    pub fn fails(&self, sigma1: f64, sigma2: f64, tau12: f64) -> bool {
        self.failure_index(sigma1, sigma2, tau12) >= 1.0
    }

    /// Check convexity condition: |F12*| < 1 (equivalently F11*F22 > F12²).
    ///
    /// Using the normalized form: convexity holds when `f12_star² < 1`.
    pub fn is_convex(&self) -> bool {
        self.f12_star * self.f12_star < 1.0
    }

    /// Strength ratio R: scale factor to failure. Solves quadratic a*R²+b*R-1=0.
    pub fn strength_ratio(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let a = self.f11() * sigma1 * sigma1
            + self.f22() * sigma2 * sigma2
            + self.f66() * tau12 * tau12
            + 2.0 * self.f12() * sigma1 * sigma2;
        let b = self.f1() * sigma1 + self.f2() * sigma2;
        if a < 1e-30 {
            if b.abs() < 1e-30 {
                return f64::INFINITY;
            }
            return 1.0 / b;
        }
        let disc = b * b + 4.0 * a;
        if disc < 0.0 {
            return f64::INFINITY;
        }
        (-b + disc.sqrt()) / (2.0 * a)
    }
}

/// Maximum stress failure criterion.
#[derive(Debug, Clone)]
pub struct MaxStressCriterion {
    /// Longitudinal tensile strength Xt \[Pa\].
    pub xt: f64,
    /// Longitudinal compressive strength Xc \[Pa\].
    pub xc: f64,
    /// Transverse tensile strength Yt \[Pa\].
    pub yt: f64,
    /// Transverse compressive strength Yc \[Pa\].
    pub yc: f64,
    /// In-plane shear strength S12 \[Pa\].
    pub s12: f64,
}

impl MaxStressCriterion {
    /// Create max-stress criterion.
    pub fn new(xt: f64, xc: f64, yt: f64, yc: f64, s12: f64) -> Self {
        Self {
            xt,
            xc,
            yt,
            yc,
            s12,
        }
    }

    /// Check failure.
    pub fn fails(&self, sigma1: f64, sigma2: f64, tau12: f64) -> bool {
        if sigma1 > 0.0 && sigma1 >= self.xt {
            return true;
        }
        if sigma1 < 0.0 && sigma1.abs() >= self.xc {
            return true;
        }
        if sigma2 > 0.0 && sigma2 >= self.yt {
            return true;
        }
        if sigma2 < 0.0 && sigma2.abs() >= self.yc {
            return true;
        }
        if tau12.abs() >= self.s12 {
            return true;
        }
        false
    }

    /// Maximum stress failure index (max ratio to strength).
    pub fn failure_index(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let f1 = if sigma1 > 0.0 {
            sigma1 / self.xt
        } else {
            sigma1.abs() / self.xc
        };
        let f2 = if sigma2 > 0.0 {
            sigma2 / self.yt
        } else {
            sigma2.abs() / self.yc
        };
        let f6 = tau12.abs() / self.s12;
        f1.max(f2).max(f6)
    }
}

/// Puck failure criterion (fiber fracture and inter-fiber fracture).
#[derive(Debug, Clone)]
pub struct PuckCriterion {
    /// Longitudinal tensile strength Xt \[Pa\].
    pub xt: f64,
    /// Longitudinal compressive strength Xc \[Pa\].
    pub xc: f64,
    /// Transverse tensile strength Yt \[Pa\].
    pub yt: f64,
    /// Transverse compressive strength Yc \[Pa\].
    pub yc: f64,
    /// Shear strength S12 \[Pa\].
    pub s12: f64,
    /// Master fracture plane inclination (typically 53°).
    pub theta_fp_deg: f64,
}

impl PuckCriterion {
    /// Create Puck criterion.
    pub fn new(xt: f64, xc: f64, yt: f64, yc: f64, s12: f64) -> Self {
        Self {
            xt,
            xc,
            yt,
            yc,
            s12,
            theta_fp_deg: 53.0,
        }
    }

    /// Fiber fracture failure index.
    pub fn fiber_fracture_index(&self, sigma1: f64) -> f64 {
        if sigma1 >= 0.0 {
            sigma1 / self.xt
        } else {
            sigma1.abs() / self.xc
        }
    }

    /// Inter-fiber fracture failure index (mode A: transverse tension).
    pub fn iff_mode_a(&self, sigma2: f64, tau12: f64) -> f64 {
        if sigma2 < 0.0 {
            return 0.0;
        }
        let p12 = 0.27; // friction coefficient (typical)
        ((tau12 / self.s12).powi(2)
            + (1.0 - p12 * self.yt / self.s12).powi(2) * (sigma2 / self.yt).powi(2))
        .sqrt()
            + p12 * sigma2 / self.s12
    }

    /// Inter-fiber fracture failure index (mode B: transverse compression).
    pub fn iff_mode_b(&self, sigma2: f64, tau12: f64) -> f64 {
        if sigma2 >= 0.0 {
            return 0.0;
        }
        let p12c = 0.32;
        let ra = self.s12 / (2.0 * p12c);
        let term1 = (tau12 / (2.0 * (1.0 + p12c) * self.s12 * ra)).powi(2);
        let term2 = (sigma2 / self.yc).powi(2);
        (term1 + term2).sqrt() + p12c * sigma2 / self.yc
    }
}

/// Progressive damage model (ply-by-ply degradation).
#[derive(Debug, Clone)]
pub struct ProgressiveDamage {
    /// Current stiffness degradation factor per ply (1 = intact).
    pub degradation: Vec<f64>,
    /// Failure criterion (Tsai-Wu).
    pub criterion: TsaiWuCriterion,
    /// Stiffness reduction factor upon failure.
    pub reduction_factor: f64,
}

impl ProgressiveDamage {
    /// Create progressive damage model for n plies.
    pub fn new(n_plies: usize, criterion: TsaiWuCriterion, reduction_factor: f64) -> Self {
        Self {
            degradation: vec![1.0; n_plies],
            criterion,
            reduction_factor,
        }
    }

    /// Apply load step and degrade failed plies.
    /// Returns (new_degradation_flags, ply failure count).
    pub fn apply_load_step(&mut self, stresses: &[[f64; 3]]) -> usize {
        let mut count = 0;
        for (i, s) in stresses.iter().enumerate() {
            if self.degradation[i] > 1e-6 {
                let fi = self.criterion.failure_index(s[0], s[1], s[2]);
                if fi >= 1.0 {
                    self.degradation[i] *= self.reduction_factor;
                    count += 1;
                }
            }
        }
        count
    }

    /// Check if all plies have failed (last-ply-failure).
    pub fn all_plies_failed(&self) -> bool {
        self.degradation.iter().all(|&d| d < 0.01)
    }

    /// Number of failed plies.
    pub fn failed_count(&self) -> usize {
        self.degradation.iter().filter(|&&d| d < 0.99).count()
    }
}

/// Delamination model (fracture mechanics approach).
#[derive(Debug, Clone)]
pub struct DelaminationModel {
    /// Mode-I interlaminar fracture toughness G_Ic \[J/m²\].
    pub g_ic: f64,
    /// Mode-II interlaminar fracture toughness G_IIc \[J/m²\].
    pub g_iic: f64,
    /// BK criterion exponent.
    pub bk_exponent: f64,
}

impl DelaminationModel {
    /// Create delamination model.
    pub fn new(g_ic: f64, g_iic: f64, bk_exponent: f64) -> Self {
        Self {
            g_ic,
            g_iic,
            bk_exponent,
        }
    }

    /// Benzeggagh-Kenane (BK) mixed-mode toughness.
    pub fn bk_toughness(&self, g_i: f64, g_ii: f64) -> f64 {
        let g_total = g_i + g_ii;
        if g_total < 1e-30 {
            return self.g_ic;
        }
        let mode_ratio = g_ii / g_total;
        self.g_ic + (self.g_iic - self.g_ic) * mode_ratio.powf(self.bk_exponent)
    }

    /// Check if delamination occurs under (G_I, G_II) loading.
    pub fn delamination_occurs(&self, g_i: f64, g_ii: f64) -> bool {
        let g_total = g_i + g_ii;
        let g_c = self.bk_toughness(g_i, g_ii);
        g_total >= g_c
    }

    /// Failure index for mixed-mode delamination (quadratic criterion).
    pub fn failure_index_quadratic(&self, g_i: f64, g_ii: f64) -> f64 {
        (g_i / self.g_ic).powi(2) + (g_ii / self.g_iic).powi(2)
    }
}

/// Woven composite geometry parameters.
#[derive(Debug, Clone)]
pub struct WovenComposite {
    /// Weave style (Plain, Twill, Satin).
    pub weave_style: WeaveStyle,
    /// Crimp angle \[degrees\].
    pub crimp_angle_deg: f64,
    /// Fiber volume fraction.
    pub vf: f64,
    /// Warp/weft ratio.
    pub warp_weft_ratio: f64,
}

/// Weave style for woven composites.
#[derive(Debug, Clone, PartialEq)]
pub enum WeaveStyle {
    /// Plain weave (1/1 interlacing).
    Plain,
    /// Twill weave (2/2 interlacing).
    Twill,
    /// Satin weave (4H, 5H, 8H).
    Satin {
        /// Number of harnesses (typically 4, 5, or 8)
        harness: u32,
    },
}

impl WovenComposite {
    /// Create a plain weave composite.
    pub fn plain(crimp_angle_deg: f64, vf: f64) -> Self {
        Self {
            weave_style: WeaveStyle::Plain,
            crimp_angle_deg,
            vf,
            warp_weft_ratio: 1.0,
        }
    }

    /// Effective modulus reduction due to crimp (cos⁴θ approximation).
    pub fn crimp_modulus_factor(&self) -> f64 {
        let theta = self.crimp_angle_deg.to_radians();
        theta.cos().powi(4)
    }

    /// Approximate in-plane modulus (Rule-of-mixtures with crimp).
    pub fn approximate_modulus(&self, fiber_e: f64, matrix_e: f64) -> f64 {
        let rom = fiber_e * self.vf + matrix_e * (1.0 - self.vf);
        rom * self.crimp_modulus_factor()
    }

    /// Interlacement density (number of crossings per unit cell).
    pub fn interlacement_density(&self) -> f64 {
        match &self.weave_style {
            WeaveStyle::Plain => 1.0,
            WeaveStyle::Twill => 0.5,
            WeaveStyle::Satin { harness } => 1.0 / (*harness as f64),
        }
    }
}

/// Homogenization scheme for composite effective properties.
#[derive(Debug, Clone)]
pub struct HomogenizationScheme;

impl HomogenizationScheme {
    /// Voigt (iso-strain) upper bound for modulus.
    pub fn voigt(e_fiber: f64, e_matrix: f64, vf: f64) -> f64 {
        e_fiber * vf + e_matrix * (1.0 - vf)
    }

    /// Reuss (iso-stress) lower bound for modulus.
    pub fn reuss(e_fiber: f64, e_matrix: f64, vf: f64) -> f64 {
        1.0 / (vf / e_fiber + (1.0 - vf) / e_matrix)
    }

    /// Hill average = (Voigt + Reuss) / 2.
    pub fn hill(e_fiber: f64, e_matrix: f64, vf: f64) -> f64 {
        (Self::voigt(e_fiber, e_matrix, vf) + Self::reuss(e_fiber, e_matrix, vf)) / 2.0
    }

    /// Hashin-Shtrikman upper bound (spherical inclusions).
    pub fn hashin_shtrikman_upper(k1: f64, g1: f64, k2: f64, vf: f64) -> f64 {
        let vm = 1.0 - vf;
        k1 + vm / (1.0 / (k2 - k1) + 3.0 * vf / (3.0 * k1 + 4.0 * g1))
    }

    /// Verify Voigt ≥ Reuss (fundamental bound).
    pub fn verify_bounds(e_fiber: f64, e_matrix: f64, vf: f64) -> bool {
        Self::voigt(e_fiber, e_matrix, vf) >= Self::reuss(e_fiber, e_matrix, vf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon_epoxy() -> (FiberProperties, MatrixProperties) {
        (FiberProperties::carbon_t300(), MatrixProperties::epoxy())
    }

    #[test]
    fn test_fiber_properties_carbon() {
        let f = FiberProperties::carbon_t300();
        assert!(f.e1 > 200e9);
        assert!(f.volume_fraction > 0.0 && f.volume_fraction < 1.0);
    }

    #[test]
    fn test_matrix_shear_modulus() {
        let m = MatrixProperties::epoxy();
        let g = m.shear_modulus();
        assert!((g - m.modulus / (2.0 * (1.0 + m.poisson_ratio))).abs() < 1.0);
    }

    #[test]
    fn test_halpin_tsai_e11_rom() {
        let (f, m) = carbon_epoxy();
        let ht = HalpinTsai::new(f.clone(), m.clone());
        // E11 is rule of mixtures
        let rom = f.e1 * f.volume_fraction + m.modulus * (1.0 - f.volume_fraction);
        assert!((ht.e11() - rom).abs() < 1.0);
    }

    #[test]
    fn test_halpin_tsai_e22_greater_reuss() {
        let (f, m) = carbon_epoxy();
        let ht = HalpinTsai::new(f, m);
        // H-T gives better (higher) E22 than Reuss lower bound
        assert!(ht.e22() >= ht.reuss_e22() - 1e6);
    }

    #[test]
    fn test_halpin_tsai_nu12_positive() {
        let (f, m) = carbon_epoxy();
        let ht = HalpinTsai::new(f, m);
        assert!(ht.nu12() > 0.0 && ht.nu12() < 0.5);
    }

    #[test]
    fn test_halpin_tsai_nu21_reciprocal() {
        let (f, m) = carbon_epoxy();
        let ht = HalpinTsai::new(f, m);
        // nu21/E22 = nu12/E11 (reciprocal relation)
        let lhs = ht.nu21() / ht.e22();
        let rhs = ht.nu12() / ht.e11();
        assert!((lhs - rhs).abs() < 1e-15);
    }

    #[test]
    fn test_voigt_reuss_bounds_fiber() {
        let (f, m) = carbon_epoxy();
        let ht = HalpinTsai::new(f, m);
        assert!(ht.voigt_e11() >= ht.reuss_e22());
    }

    #[test]
    fn test_mori_tanaka_bulk_positive() {
        let (f, m) = carbon_epoxy();
        let mt = MoriTanakaComposite::new(f, m);
        assert!(mt.effective_bulk_modulus() > 0.0);
    }

    #[test]
    fn test_mori_tanaka_shear_positive() {
        let (f, m) = carbon_epoxy();
        let mt = MoriTanakaComposite::new(f, m);
        assert!(mt.effective_shear_modulus() > 0.0);
    }

    #[test]
    fn test_lamina_q_matrix_positive_definite() {
        let l = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, 0.0);
        let q = l.q_matrix();
        assert!(q[0] > 0.0 && q[2] > 0.0 && q[3] > 0.0);
    }

    #[test]
    fn test_lamina_nu21_reciprocal() {
        let l = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, 0.0);
        let nu21 = l.nu21();
        // nu21 = nu12 * E2/E1
        let expected = 0.28 * 8e9 / 130e9;
        assert!((nu21 - expected).abs() < 1e-15);
    }

    #[test]
    fn test_q_bar_zero_angle_equals_q() {
        let l = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, 0.0);
        let q = l.q_matrix();
        let qb = l.q_bar();
        assert!((qb[0] - q[0]).abs() < 1e3); // Q11_bar ≈ Q11
        assert!((qb[3] - q[2]).abs() < 1e3); // Q22_bar ≈ Q22
    }

    #[test]
    fn test_clt_total_thickness() {
        let l1 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 0.0);
        let l2 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 90.0);
        let clt = ClassicalLaminateTheory::new(vec![l1, l2]);
        assert!((clt.total_thickness() - 0.5e-3).abs() < 1e-10);
    }

    #[test]
    fn test_clt_z_coords_symmetric() {
        let l1 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 0.0);
        let l2 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 90.0);
        let clt = ClassicalLaminateTheory::new(vec![l1, l2]);
        let z = clt.z_coords();
        assert_eq!(z.len(), 3);
        assert!((z[0] + z[2]).abs() < 1e-15); // symmetric about 0
    }

    #[test]
    fn test_clt_a_matrix_nonzero() {
        let l1 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 0.0);
        let l2 = Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.25e-3, 90.0);
        let clt = ClassicalLaminateTheory::new(vec![l1, l2]);
        let a = clt.a_matrix();
        assert!(a[0] > 0.0 && a[3] > 0.0);
    }

    #[test]
    fn test_clt_symmetric_laminate_b_zero() {
        // [0/90/90/0] is symmetric -> B=0
        let mk_l = |angle: f64| Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, angle);
        let clt = ClassicalLaminateTheory::new(vec![mk_l(0.0), mk_l(90.0), mk_l(90.0), mk_l(0.0)]);
        assert!(clt.is_symmetric());
    }

    #[test]
    fn test_clt_balanced_laminate() {
        // [0/+45/-45/90] balanced -> A16=A26≈0
        let mk_l = |angle: f64| Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, angle);
        let clt = ClassicalLaminateTheory::new(vec![
            mk_l(0.0),
            mk_l(45.0),
            mk_l(-45.0),
            mk_l(90.0),
            mk_l(90.0),
            mk_l(-45.0),
            mk_l(45.0),
            mk_l(0.0),
        ]);
        // Symmetric so B=0; check A16/A26 small
        let a = clt.a_matrix();
        assert!(a[2].abs() < a[0] * 1e-10);
    }

    #[test]
    fn test_tsai_wu_convexity() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        assert!(tw.is_convex());
    }

    #[test]
    fn test_tsai_wu_origin_safe() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        assert!(!tw.fails(0.0, 0.0, 0.0));
        assert!(tw.failure_index(0.0, 0.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_tsai_wu_longitudinal_tensile_failure() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        // Just above Xt should fail
        assert!(tw.fails(1600e6, 0.0, 0.0));
    }

    #[test]
    fn test_tsai_wu_strength_ratio_positive() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        let r = tw.strength_ratio(500e6, 10e6, 5e6);
        assert!(r > 0.0);
    }

    #[test]
    fn test_max_stress_safe_origin() {
        let ms = MaxStressCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        assert!(!ms.fails(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_max_stress_fiber_failure() {
        let ms = MaxStressCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        assert!(ms.fails(1600e6, 0.0, 0.0));
    }

    #[test]
    fn test_max_stress_shear_failure() {
        let ms = MaxStressCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        assert!(ms.fails(0.0, 0.0, 80e6));
    }

    #[test]
    fn test_max_stress_failure_index_range() {
        let ms = MaxStressCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        let fi = ms.failure_index(500e6, 20e6, 30e6);
        assert!(fi > 0.0);
    }

    #[test]
    fn test_puck_fiber_fracture_tension() {
        let p = PuckCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        let fi = p.fiber_fracture_index(1600e6);
        assert!(fi > 1.0);
    }

    #[test]
    fn test_puck_iff_mode_a_transverse_tension() {
        let p = PuckCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6);
        // Large transverse tension should give fi > 1
        let fi = p.iff_mode_a(60e6, 0.0);
        assert!(fi > 1.0);
    }

    #[test]
    fn test_progressive_damage_count() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        let mut pd = ProgressiveDamage::new(3, tw, 0.01);
        let stresses = [
            [1600e6_f64, 0.0, 0.0],
            [100e6, 10e6, 5e6],
            [1700e6, 0.0, 0.0],
        ];
        let failed = pd.apply_load_step(&stresses);
        assert_eq!(failed, 2);
    }

    #[test]
    fn test_progressive_damage_not_all_failed() {
        let tw = TsaiWuCriterion::new(1500e6, 1200e6, 50e6, 200e6, 70e6, -0.5);
        let pd = ProgressiveDamage::new(4, tw, 0.5);
        assert!(!pd.all_plies_failed());
    }

    #[test]
    fn test_delamination_bk_toughness_mode1() {
        let dm = DelaminationModel::new(200.0, 800.0, 2.284);
        let gc = dm.bk_toughness(200.0, 0.0);
        assert!((gc - 200.0).abs() < 1e-6);
    }

    #[test]
    fn test_delamination_bk_toughness_mode2() {
        let dm = DelaminationModel::new(200.0, 800.0, 2.284);
        let gc = dm.bk_toughness(0.0, 800.0);
        assert!((gc - 800.0).abs() < 1.0);
    }

    #[test]
    fn test_delamination_failure_detection() {
        let dm = DelaminationModel::new(200.0, 800.0, 2.284);
        assert!(dm.delamination_occurs(300.0, 0.0)); // G_I > G_Ic
        assert!(!dm.delamination_occurs(100.0, 0.0)); // G_I < G_Ic
    }

    #[test]
    fn test_woven_composite_crimp_factor_leq1() {
        let wc = WovenComposite::plain(5.0, 0.5);
        let factor = wc.crimp_modulus_factor();
        assert!(factor > 0.0 && factor <= 1.0);
    }

    #[test]
    fn test_woven_composite_modulus_less_than_rom() {
        let wc = WovenComposite::plain(5.0, 0.5);
        let rom = 200e9 * 0.5 + 3e9 * 0.5;
        let woven_mod = wc.approximate_modulus(200e9, 3e9);
        assert!(woven_mod < rom);
    }

    #[test]
    fn test_woven_satin_lower_interlacement() {
        let plain = WovenComposite::plain(5.0, 0.5);
        let satin = WovenComposite {
            weave_style: WeaveStyle::Satin { harness: 8 },
            crimp_angle_deg: 5.0,
            vf: 0.5,
            warp_weft_ratio: 1.0,
        };
        assert!(plain.interlacement_density() > satin.interlacement_density());
    }

    #[test]
    fn test_homogenization_voigt_reuss_bounds() {
        let vf = 0.6;
        let ef = 230e9;
        let em = 3.5e9;
        assert!(HomogenizationScheme::verify_bounds(ef, em, vf));
    }

    #[test]
    fn test_homogenization_voigt_gt_hill_gt_reuss() {
        let vf = 0.6;
        let ef = 230e9;
        let em = 3.5e9;
        let v = HomogenizationScheme::voigt(ef, em, vf);
        let h = HomogenizationScheme::hill(ef, em, vf);
        let r = HomogenizationScheme::reuss(ef, em, vf);
        assert!(v >= h);
        assert!(h >= r);
    }

    #[test]
    fn test_homogenization_vf0_gives_matrix() {
        let em = 3.5e9;
        let ef = 230e9;
        let v = HomogenizationScheme::voigt(ef, em, 0.0);
        let r = HomogenizationScheme::reuss(ef, em, 0.0);
        assert!((v - em).abs() < 1.0);
        assert!((r - em).abs() < 1.0);
    }

    #[test]
    fn test_homogenization_vf1_gives_fiber() {
        let em = 3.5e9;
        let ef = 230e9;
        let v = HomogenizationScheme::voigt(ef, em, 1.0);
        let r = HomogenizationScheme::reuss(ef, em, 1.0);
        assert!((v - ef).abs() < 1.0);
        assert!((r - ef).abs() < 1.0);
    }

    #[test]
    fn test_laminate_response_symmetric() {
        let mk_l = |angle: f64| Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, angle);
        let clt = ClassicalLaminateTheory::new(vec![mk_l(0.0), mk_l(90.0), mk_l(90.0), mk_l(0.0)]);
        let a = clt.a_matrix();
        let resp = LaminateResponse::from_n_symmetric([1000.0, 0.0, 0.0], a);
        assert!(resp.midplane_strains[0].abs() > 0.0);
    }

    #[test]
    fn test_clt_d_matrix_positive_diagonal() {
        let mk_l = |angle: f64| Lamina::new(130e9, 8e9, 4.5e9, 0.28, 0.125e-3, angle);
        let clt = ClassicalLaminateTheory::new(vec![mk_l(0.0), mk_l(90.0)]);
        let d = clt.d_matrix();
        assert!(d[0] > 0.0 && d[3] > 0.0 && d[5] > 0.0);
    }

    #[test]
    fn test_mori_tanaka_bounds() {
        let (f, m) = carbon_epoxy();
        let em = m.modulus;
        let ef = f.e1;
        let mt = MoriTanakaComposite::new(f, m);
        let k_eff = mt.effective_bulk_modulus();
        let km = em / (3.0 * (1.0 - 2.0 * 0.35));
        let kf = ef / (3.0 * (1.0 - 2.0 * 0.27));
        // Effective should be between matrix and fiber bulk moduli
        assert!(k_eff >= km.min(kf) - 1e9);
        assert!(k_eff <= km.max(kf) + 1e9);
    }
}
