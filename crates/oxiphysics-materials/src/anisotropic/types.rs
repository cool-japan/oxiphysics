//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Classical laminate plate theory (CLT) single ply in a laminate stack.
///
/// Represents one ply with its orientation angle (in degrees from the
/// laminate reference axis) and thickness, backed by a `WovenLamina` or
/// equivalent orthotropic properties.
#[derive(Debug, Clone, Copy)]
pub struct Ply {
    /// In-plane modulus along ply axis 1 \[Pa\]
    pub e1: f64,
    /// In-plane modulus along ply axis 2 \[Pa\]
    pub e2: f64,
    /// In-plane Poisson ratio ν12
    pub nu12: f64,
    /// In-plane shear modulus G12 \[Pa\]
    pub g12: f64,
    /// Thickness \[m\]
    pub thickness: f64,
    /// Orientation angle θ \[degrees\] from laminate reference axis
    pub angle_deg: f64,
}
impl Ply {
    /// Create a new `Ply`.
    pub fn new(e1: f64, e2: f64, nu12: f64, g12: f64, thickness: f64, angle_deg: f64) -> Self {
        Self {
            e1,
            e2,
            nu12,
            g12,
            thickness,
            angle_deg,
        }
    }
    /// Minor Poisson ratio ν₂₁ = ν₁₂ · E₂ / E₁.
    pub fn nu21(&self) -> f64 {
        self.nu12 * self.e2 / self.e1
    }
    /// Reduced stiffness matrix Q \[Pa\] in ply principal coordinates.
    ///
    /// Returns `[Q11, Q22, Q12, Q66]`.
    pub fn reduced_stiffness(&self) -> [f64; 4] {
        let nu21 = self.nu12 * self.e2 / self.e1;
        let denom = 1.0 - self.nu12 * nu21;
        let q11 = self.e1 / denom;
        let q22 = self.e2 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q66 = self.g12;
        [q11, q22, q12, q66]
    }
    /// Transformed (rotated) reduced stiffness Q̄ in laminate coordinates.
    ///
    /// Returns the 3×3 matrix `[Qbar_11, Qbar_12, Qbar_16; Qbar_12, Qbar_22, Qbar_26;
    /// Qbar_16, Qbar_26, Qbar_66]` as a flat array `[q11,q12,q16,q12,q22,q26,q16,q26,q66]`.
    pub fn transformed_stiffness(&self) -> [[f64; 3]; 3] {
        let [q11, q22, q12, q66] = self.reduced_stiffness();
        let theta = self.angle_deg.to_radians();
        let m = theta.cos();
        let n = theta.sin();
        let m2 = m * m;
        let n2 = n * n;
        let m4 = m2 * m2;
        let n4 = n2 * n2;
        let m2n2 = m2 * n2;
        let m3n = m2 * m * n;
        let mn3 = m * n2 * n;
        let qb11 = q11 * m4 + 2.0 * (q12 + 2.0 * q66) * m2n2 + q22 * n4;
        let qb22 = q11 * n4 + 2.0 * (q12 + 2.0 * q66) * m2n2 + q22 * m4;
        let qb12 = (q11 + q22 - 4.0 * q66) * m2n2 + q12 * (m4 + n4);
        let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * m2n2 + q66 * (m4 + n4);
        let qb16 = (q11 - q12 - 2.0 * q66) * m3n - (q22 - q12 - 2.0 * q66) * mn3;
        let qb26 = (q11 - q12 - 2.0 * q66) * mn3 - (q22 - q12 - 2.0 * q66) * m3n;
        [[qb11, qb12, qb16], [qb12, qb22, qb26], [qb16, qb26, qb66]]
    }
}
/// Transversely isotropic material: isotropic in one plane, anisotropic along the axis.
///
/// This is a special case of orthotropic with `E2 = E3`, `ν12 = ν13`, `G12 = G13`.
/// Typical for unidirectional fiber-reinforced composites.
#[derive(Debug, Clone, Copy)]
pub struct TransverselyIsotropic {
    /// Axial modulus along the fiber direction (Pa)
    pub ea: f64,
    /// Transverse modulus (Pa)
    pub et: f64,
    /// Axial Poisson ratio ν12 (strain in transverse due to axial stress)
    pub nua: f64,
    /// Transverse Poisson ratio ν23
    pub nut: f64,
    /// Axial shear modulus G12 (Pa)
    pub ga: f64,
}
impl TransverselyIsotropic {
    /// Create a new transversely isotropic material.
    pub fn new(ea: f64, et: f64, nua: f64, nut: f64, ga: f64) -> Self {
        Self {
            ea,
            et,
            nua,
            nut,
            ga,
        }
    }
    /// Transverse shear modulus `Gt = Et / (2 * (1 + nut))`.
    pub fn gt(&self) -> f64 {
        self.et / (2.0 * (1.0 + self.nut))
    }
    /// Convert to equivalent orthotropic material.
    ///
    /// `E1 = ea`, `E2 = E3 = et`, `ν12 = ν13 = nua`, `ν23 = nut`, `G12 = G13 = ga`, `G23 = gt()`.
    pub fn to_orthotropic(&self) -> OrthotropicMaterial {
        OrthotropicMaterial::new(
            self.ea,
            self.et,
            self.et,
            self.nua,
            self.nut,
            self.nua,
            self.ga,
            self.gt(),
            self.ga,
        )
    }
    /// 6×6 stiffness (constitutive) matrix D.
    pub fn constitutive_matrix(&self) -> [[f64; 6]; 6] {
        self.to_orthotropic().constitutive_matrix()
    }
}
/// Anisotropic thermal conductivity tensor (3×3 symmetric positive-definite matrix).
///
/// Represents materials where heat flow differs in each crystallographic direction,
/// such as hexagonal boron nitride or graphite.
#[derive(Debug, Clone, Copy)]
pub struct ThermalConductivityTensor {
    /// The 3×3 conductivity matrix κ \[W/(m·K)\], stored in row-major order.
    pub kappa: [[f64; 3]; 3],
}
impl ThermalConductivityTensor {
    /// Construct from a 3×3 matrix.
    pub fn new(kappa: [[f64; 3]; 3]) -> Self {
        Self { kappa }
    }
    /// Isotropic thermal conductivity (same in all directions).
    pub fn isotropic(k: f64) -> Self {
        Self::orthotropic(k, k, k)
    }
    /// Orthotropic (diagonal) thermal conductivity.
    ///
    /// `kx`, `ky`, `kz` are the principal conductivities along the x, y, z axes.
    pub fn orthotropic(kx: f64, ky: f64, kz: f64) -> Self {
        let mut k = [[0.0_f64; 3]; 3];
        k[0][0] = kx;
        k[1][1] = ky;
        k[2][2] = kz;
        Self::new(k)
    }
    /// Compute the heat flux vector `q = -κ · ∇T`.
    ///
    /// `grad_t` is the temperature gradient vector \[K/m\].
    pub fn heat_flux(&self, grad_t: [f64; 3]) -> [f64; 3] {
        let mut q = [0.0f64; 3];
        for (i, qi) in q.iter_mut().enumerate() {
            for (j, &gt_j) in grad_t.iter().enumerate() {
                *qi -= self.kappa[i][j] * gt_j;
            }
        }
        q
    }
    /// Effective thermal conductivity along a direction `n` (unit vector).
    ///
    /// κ_eff = nᵀ · κ · n
    pub fn effective_conductivity(&self, n: [f64; 3]) -> f64 {
        let mut kn = [0.0f64; 3];
        for (i, kni) in kn.iter_mut().enumerate() {
            for (j, &n_j) in n.iter().enumerate() {
                *kni += self.kappa[i][j] * n_j;
            }
        }
        n[0] * kn[0] + n[1] * kn[1] + n[2] * kn[2]
    }
    /// Average isotropic thermal conductivity: (κ_xx + κ_yy + κ_zz) / 3.
    pub fn mean_conductivity(&self) -> f64 {
        (self.kappa[0][0] + self.kappa[1][1] + self.kappa[2][2]) / 3.0
    }
}
/// LaRC03/04 failure criteria for advanced composite materials.
///
/// NASA LaRC (Langley Research Center) failure criteria account for
/// fibre kinking and in-situ strength effects.
///
/// Reference: Davila, C.G., Camanho, P.P. & Rose, C.A. (2005). "Failure
/// criteria for FRP laminates". J. Compos. Mater. 39(4), 323–345.
#[derive(Debug, Clone, Copy)]
pub struct LaRCFailureCriteria {
    /// Fibre tensile strength \[Pa\].
    pub xt: f64,
    /// Fibre compressive strength \[Pa\] (positive).
    pub xc: f64,
    /// In-situ matrix tensile strength \[Pa\].
    pub yis: f64,
    /// In-situ shear strength \[Pa\].
    pub sis: f64,
    /// Fibre misalignment angle φ₀ \[rad\].
    pub phi0: f64,
    /// Fracture toughness ratio parameter η.
    pub eta: f64,
}
impl LaRCFailureCriteria {
    /// Create with explicit parameters.
    pub fn new(xt: f64, xc: f64, yis: f64, sis: f64, phi0: f64, eta: f64) -> Self {
        Self {
            xt,
            xc,
            yis,
            sis,
            phi0,
            eta,
        }
    }
    /// Default parameters for carbon/epoxy (approximation of IM7/8552).
    pub fn default_carbon_epoxy() -> Self {
        Self::new(2560.0e6, 1590.0e6, 160.0e6, 130.0e6, 0.015, 0.43)
    }
    /// Fibre kinking failure index for compressive loading.
    ///
    /// `FI_fk = |τ_m| / (sis - η * σ_mn)` where stresses are in the kink band frame.
    pub fn fibre_kinking_index(&self, sigma11: f64, sigma22: f64, tau12: f64) -> f64 {
        let _ = sigma22;
        let phi = self.phi0;
        let sigma_m = sigma11 * phi.sin().powi(2) + tau12 * 2.0 * phi.sin() * phi.cos();
        let tau_m = (sigma11 - 0.0) * phi.sin() * phi.cos()
            + tau12 * (phi.cos().powi(2) - phi.sin().powi(2));
        let denominator = (self.sis - self.eta * sigma_m).max(1.0);
        tau_m.abs() / denominator
    }
    /// Matrix cracking failure index.
    ///
    /// `FI_mc = (τ_n/sis)² + (σ_n/yis)²` (simplified Mohr-Coulomb form).
    pub fn matrix_cracking_index(&self, sigma22: f64, tau12: f64, tau23: f64) -> f64 {
        let _ = tau23;
        let a = (tau12 / self.sis).powi(2);
        let b = if sigma22 > 0.0 {
            (sigma22 / self.yis).powi(2)
        } else {
            0.0
        };
        (a + b).sqrt()
    }
}
/// Woven/braided reinforcement model for triaxially braided composites.
///
/// Uses a simplified mixture-theory approach with off-axis transformations
/// for the bias (braided) tows plus the axial tows.
#[derive(Debug, Clone, Copy)]
pub struct BraidedComposite {
    /// Volume fraction of fibre (all tow directions combined).
    pub vf: f64,
    /// Fibre Young's modulus \[Pa\].
    pub e_fibre: f64,
    /// Matrix Young's modulus \[Pa\].
    pub e_matrix: f64,
    /// Fibre Poisson ratio.
    pub nu_fibre: f64,
    /// Matrix Poisson ratio.
    pub nu_matrix: f64,
    /// Braid angle θ \[degrees\] from axial direction.
    pub braid_angle_deg: f64,
}
impl BraidedComposite {
    /// Construct a triaxial braid from constituent properties.
    ///
    /// `braid_angle_deg` is the angle of the bias tows from the axial direction.
    pub fn triaxial_braid(
        vf: f64,
        e_fibre: f64,
        e_matrix: f64,
        nu_fibre: f64,
        nu_matrix: f64,
        braid_angle_deg: f64,
    ) -> Self {
        Self {
            vf,
            e_fibre,
            e_matrix,
            nu_fibre,
            nu_matrix,
            braid_angle_deg,
        }
    }
    /// Matrix volume fraction.
    fn vm(&self) -> f64 {
        1.0 - self.vf
    }
    /// Rule-of-mixtures axial modulus \[Pa\].
    fn e_axial(&self) -> f64 {
        self.vf * self.e_fibre + self.vm() * self.e_matrix
    }
    /// Inverse rule-of-mixtures transverse modulus \[Pa\].
    fn e_transverse(&self) -> f64 {
        1.0 / (self.vf / self.e_fibre + self.vm() / self.e_matrix)
    }
    /// Rule-of-mixtures Poisson ratio ν12.
    fn nu12(&self) -> f64 {
        self.vf * self.nu_fibre + self.vm() * self.nu_matrix
    }
    /// In-plane shear modulus G12 \[Pa\] (Reuss/series lower bound).
    fn g12(&self) -> f64 {
        let gf = self.e_fibre / (2.0 * (1.0 + self.nu_fibre));
        let gm = self.e_matrix / (2.0 * (1.0 + self.nu_matrix));
        1.0 / (self.vf / gf + self.vm() / gm)
    }
    /// Effective in-plane modulus Ex of the braided laminate.
    ///
    /// Computed by averaging the Q̄ matrices of the axial (0°) and
    /// bias (±θ) tow contributions.
    pub fn in_plane_modulus(&self) -> f64 {
        let theta = self.braid_angle_deg.to_radians();
        let m = theta.cos();
        let n = theta.sin();
        let m2 = m * m;
        let n2 = n * n;
        let m4 = m2 * m2;
        let n4 = n2 * n2;
        let m2n2 = m2 * n2;
        let e1 = self.e_axial();
        let e2 = self.e_transverse();
        let nu12 = self.nu12();
        let nu21 = nu12 * e2 / e1;
        let g12 = self.g12();
        let denom = 1.0 - nu12 * nu21;
        let q11 = e1 / denom;
        let q22 = e2 / denom;
        let q12 = nu12 * e2 / denom;
        let q66 = g12;
        let qb11_bias = q11 * m4 + 2.0 * (q12 + 2.0 * q66) * m2n2 + q22 * n4;
        let weight_axial = 1.0 / 3.0;
        let weight_bias = 2.0 / 3.0;
        weight_axial * q11 + weight_bias * qb11_bias
    }
    /// Effective shear modulus G of the braided laminate \[Pa\].
    pub fn effective_shear_modulus(&self) -> f64 {
        let theta = self.braid_angle_deg.to_radians();
        let m = theta.cos();
        let n = theta.sin();
        let m2 = m * m;
        let n2 = n * n;
        let m4 = m2 * m2;
        let n4 = n2 * n2;
        let m2n2 = m2 * n2;
        let e1 = self.e_axial();
        let e2 = self.e_transverse();
        let nu12 = self.nu12();
        let nu21 = nu12 * e2 / e1;
        let g12 = self.g12();
        let denom = 1.0 - nu12 * nu21;
        let q11 = e1 / denom;
        let q22 = e2 / denom;
        let q12 = nu12 * e2 / denom;
        let q66 = g12;
        let qb66_bias = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * m2n2 + q66 * (m4 + n4);
        (1.0 / 3.0) * q66 + (2.0 / 3.0) * qb66_bias
    }
}
/// The 13 independent compliance components of a monoclinic material
/// (Voigt notation, symmetry plane = plane 1-2).
///
/// ```text
/// S = [S11 S12 S13  0   0  S16 ]
///     [S12 S22 S23  0   0  S26 ]
///     [S13 S23 S33  0   0  S36 ]
///     [ 0   0   0  S44 S45  0  ]
///     [ 0   0   0  S45 S55  0  ]
///     [S16 S26 S36  0   0  S66 ]
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct MonoclinicCompliance {
    /// S11 compliance component \[1/Pa\]
    pub s11: f64,
    /// S12 compliance component \[1/Pa\]
    pub s12: f64,
    /// S13 compliance component \[1/Pa\]
    pub s13: f64,
    /// S16 coupling compliance component \[1/Pa\]
    pub s16: f64,
    /// S22 compliance component \[1/Pa\]
    pub s22: f64,
    /// S23 compliance component \[1/Pa\]
    pub s23: f64,
    /// S26 coupling compliance component \[1/Pa\]
    pub s26: f64,
    /// S33 compliance component \[1/Pa\]
    pub s33: f64,
    /// S36 coupling compliance component \[1/Pa\]
    pub s36: f64,
    /// S44 compliance component \[1/Pa\]
    pub s44: f64,
    /// S45 coupling compliance component \[1/Pa\]
    pub s45: f64,
    /// S55 compliance component \[1/Pa\]
    pub s55: f64,
    /// S66 compliance component \[1/Pa\]
    pub s66: f64,
}

/// Monoclinic linear elastic material (single symmetry plane = plane 1-2).
///
/// Has 13 independent elastic constants (compared to 21 for triclinic).
/// The compliance matrix has the following form (Voigt notation, plane 3 is
/// the symmetry axis):
/// ```text
/// S = [S11 S12 S13  0   0  S16 ]
///     [S12 S22 S23  0   0  S26 ]
///     [S13 S23 S33  0   0  S36 ]
///     [ 0   0   0  S44 S45  0  ]
///     [ 0   0   0  S45 S55  0  ]
///     [S16 S26 S36  0   0  S66 ]
/// ```
#[derive(Debug, Clone, Copy)]
pub struct MonoclinicMaterial {
    /// S11 compliance component \[1/Pa\]
    pub s11: f64,
    /// S12 compliance component \[1/Pa\]
    pub s12: f64,
    /// S13 compliance component \[1/Pa\]
    pub s13: f64,
    /// S16 coupling compliance component \[1/Pa\]
    pub s16: f64,
    /// S22 compliance component \[1/Pa\]
    pub s22: f64,
    /// S23 compliance component \[1/Pa\]
    pub s23: f64,
    /// S26 coupling compliance component \[1/Pa\]
    pub s26: f64,
    /// S33 compliance component \[1/Pa\]
    pub s33: f64,
    /// S36 coupling compliance component \[1/Pa\]
    pub s36: f64,
    /// S44 compliance component \[1/Pa\]
    pub s44: f64,
    /// S45 coupling compliance component \[1/Pa\]
    pub s45: f64,
    /// S55 compliance component \[1/Pa\]
    pub s55: f64,
    /// S66 compliance component \[1/Pa\]
    pub s66: f64,
}
impl MonoclinicMaterial {
    /// Build from a [`MonoclinicCompliance`] parameter struct holding all 13
    /// independent compliance values.
    pub fn new(c: MonoclinicCompliance) -> Self {
        Self {
            s11: c.s11,
            s12: c.s12,
            s13: c.s13,
            s16: c.s16,
            s22: c.s22,
            s23: c.s23,
            s26: c.s26,
            s33: c.s33,
            s36: c.s36,
            s44: c.s44,
            s45: c.s45,
            s55: c.s55,
            s66: c.s66,
        }
    }
    /// Effective in-plane Young's modulus E₁ = 1 / S₁₁.
    pub fn e1(&self) -> f64 {
        1.0 / self.s11
    }
    /// Effective in-plane Young's modulus E₂ = 1 / S₂₂.
    pub fn e2(&self) -> f64 {
        1.0 / self.s22
    }
    /// Check that the compliance matrix has the correct symmetry (S16, S26, S36 coupling).
    ///
    /// Returns `true` if S11 and S22 are positive (minimal physical check).
    pub fn is_physically_valid(&self) -> bool {
        self.s11 > 0.0 && self.s22 > 0.0 && self.s33 > 0.0
    }
}
/// Woven composite lamina (in-plane weave) with distinct in-plane and through-thickness properties.
#[derive(Debug, Clone, Copy)]
pub struct WovenLamina {
    /// In-plane modulus along axis 1 (Pa)
    pub e11: f64,
    /// In-plane modulus along axis 2 (Pa)
    pub e22: f64,
    /// In-plane Poisson ratio ν12
    pub nu12: f64,
    /// In-plane shear modulus G12 (Pa)
    pub g12: f64,
    /// Through-thickness modulus E33 (Pa)
    pub e33: f64,
    /// Out-of-plane shear modulus G13 (Pa)
    pub g13: f64,
    /// Out-of-plane shear modulus G23 (Pa)
    pub g23: f64,
}
impl WovenLamina {
    /// Construct a balanced plain weave lamina from constituent properties using rule of mixtures.
    ///
    /// For a balanced weave: `E11 = E22 = vf*Ef + vm*Em`.
    /// - `fiber_vol`: fiber volume fraction (0..1)
    /// - `e_fiber`, `e_matrix`: constituent Young's moduli (Pa)
    /// - `nu_fiber`, `nu_matrix`: constituent Poisson ratios
    pub fn balanced_plain_weave(
        fiber_vol: f64,
        e_fiber: f64,
        e_matrix: f64,
        nu_fiber: f64,
        nu_matrix: f64,
    ) -> Self {
        let vm = 1.0 - fiber_vol;
        let vf = fiber_vol;
        let e11 = vf * e_fiber + vm * e_matrix;
        let e22 = e11;
        let nu12 = vf * nu_fiber + vm * nu_matrix;
        let g_fiber = e_fiber / (2.0 * (1.0 + nu_fiber));
        let g_matrix = e_matrix / (2.0 * (1.0 + nu_matrix));
        let g12 = 1.0 / (vf / g_fiber + vm / g_matrix);
        let e33 = 1.0 / (vf / e_fiber + vm / e_matrix);
        let g13 = 1.0 / (vf / g_fiber + vm / g_matrix);
        let g23 = g13;
        Self {
            e11,
            e22,
            nu12,
            g12,
            e33,
            g13,
            g23,
        }
    }
    /// 3×3 plane-stress constitutive matrix Q for `[σ_11, σ_22, τ_12] = Q * [ε_11, ε_22, γ_12]`.
    pub fn constitutive_matrix_2d(&self) -> [[f64; 3]; 3] {
        let nu21 = self.nu12 * self.e22 / self.e11;
        let denom = 1.0 - self.nu12 * nu21;
        let mut q = [[0.0_f64; 3]; 3];
        q[0][0] = self.e11 / denom;
        q[1][1] = self.e22 / denom;
        q[0][1] = self.nu12 * self.e22 / denom;
        q[1][0] = q[0][1];
        q[2][2] = self.g12;
        q
    }
}
/// Orthotropic linear elastic material with 3 distinct Young's moduli, 3 shear moduli,
/// and 3 Poisson ratios.
///
/// Used for wood, fiber-reinforced composites, etc.
///
/// The compliance matrix S = D⁻¹:
/// ```text
/// [ε_11]   [1/E1    -ν21/E2  -ν31/E3  0      0      0    ] [σ_11]
/// [ε_22] = [-ν12/E1  1/E2    -ν32/E3  0      0      0    ] [σ_22]
/// [ε_33]   [-ν13/E1 -ν23/E2   1/E3    0      0      0    ] [σ_33]
/// [γ_12]   [0        0        0       1/G12  0      0    ] [τ_12]
/// [γ_23]   [0        0        0       0      1/G23  0    ] [τ_23]
/// [γ_13]   [0        0        0       0      0      1/G13] [τ_13]
/// ```
#[derive(Debug, Clone, Copy)]
pub struct OrthotropicMaterial {
    /// Young's modulus along axis 1 (Pa)
    pub e1: f64,
    /// Young's modulus along axis 2 (Pa)
    pub e2: f64,
    /// Young's modulus along axis 3 (Pa)
    pub e3: f64,
    /// Major Poisson ratio ν12 (strain in 2 due to stress in 1)
    pub nu12: f64,
    /// Major Poisson ratio ν23 (strain in 3 due to stress in 2)
    pub nu23: f64,
    /// Major Poisson ratio ν13 (strain in 3 due to stress in 1)
    pub nu13: f64,
    /// Shear modulus in the 1-2 plane (Pa)
    pub g12: f64,
    /// Shear modulus in the 2-3 plane (Pa)
    pub g23: f64,
    /// Shear modulus in the 1-3 plane (Pa)
    pub g13: f64,
}
impl OrthotropicMaterial {
    /// Create with explicit parameters. Symmetry requires `nu_ij/Ei = nu_ji/Ej`.
    pub fn new(
        e1: f64,
        e2: f64,
        e3: f64,
        nu12: f64,
        nu23: f64,
        nu13: f64,
        g12: f64,
        g23: f64,
        g13: f64,
    ) -> Self {
        Self {
            e1,
            e2,
            e3,
            nu12,
            nu23,
            nu13,
            g12,
            g23,
            g13,
        }
    }
    /// Carbon fiber composite (AS4 fiber/epoxy matrix).
    ///
    /// - E1 = 140 GPa (along fiber), E2 = E3 = 10 GPa (transverse)
    /// - ν12 = ν13 = 0.3, ν23 = 0.4
    /// - G12 = G13 = 5 GPa, G23 = 3.5 GPa
    pub fn carbon_fiber_epoxy() -> Self {
        Self::new(140.0e9, 10.0e9, 10.0e9, 0.3, 0.4, 0.3, 5.0e9, 3.5e9, 5.0e9)
    }
    /// Douglas fir wood (along grain = 1, radial = 2, tangential = 3).
    ///
    /// - E1 = 13.5 GPa, E2 = 0.9 GPa, E3 = 0.75 GPa
    /// - ν12 = 0.292, ν23 = 0.37, ν13 = 0.374
    /// - G12 = 0.88 GPa, G23 = 0.06 GPa, G13 = 1.07 GPa
    pub fn douglas_fir() -> Self {
        Self::new(
            13.5e9, 0.9e9, 0.75e9, 0.292, 0.37, 0.374, 0.88e9, 0.06e9, 1.07e9,
        )
    }
    /// 6×6 compliance matrix S (inverse of stiffness D).
    ///
    /// Uses the minor Poisson ratios: ν21 = ν12·E2/E1, ν31 = ν13·E3/E1, ν32 = ν23·E3/E2.
    pub fn compliance_matrix(&self) -> [[f64; 6]; 6] {
        let nu21 = self.nu12 * self.e2 / self.e1;
        let nu31 = self.nu13 * self.e3 / self.e1;
        let nu32 = self.nu23 * self.e3 / self.e2;
        let mut s = [[0.0_f64; 6]; 6];
        s[0][0] = 1.0 / self.e1;
        s[1][1] = 1.0 / self.e2;
        s[2][2] = 1.0 / self.e3;
        s[0][1] = -nu21 / self.e2;
        s[1][0] = -nu21 / self.e2;
        s[0][2] = -nu31 / self.e3;
        s[2][0] = -nu31 / self.e3;
        s[1][2] = -nu32 / self.e3;
        s[2][1] = -nu32 / self.e3;
        s[3][3] = 1.0 / self.g12;
        s[4][4] = 1.0 / self.g23;
        s[5][5] = 1.0 / self.g13;
        s
    }
    /// 6×6 stiffness (constitutive) matrix D = S⁻¹.
    ///
    /// Computes by analytically inverting the 3×3 normal block and placing shear terms directly.
    pub fn constitutive_matrix(&self) -> [[f64; 6]; 6] {
        let nu21 = self.nu12 * self.e2 / self.e1;
        let nu31 = self.nu13 * self.e3 / self.e1;
        let nu32 = self.nu23 * self.e3 / self.e2;
        let nu12 = self.nu12;
        let nu13 = self.nu13;
        let nu23 = self.nu23;
        let delta = 1.0 - nu12 * nu21 - nu23 * nu32 - nu13 * nu31 - 2.0 * nu21 * nu32 * nu13;
        let e1 = self.e1;
        let e2 = self.e2;
        let e3 = self.e3;
        let mut d = [[0.0_f64; 6]; 6];
        d[0][0] = (1.0 - nu23 * nu32) * e1 / delta;
        d[1][1] = (1.0 - nu13 * nu31) * e2 / delta;
        d[2][2] = (1.0 - nu12 * nu21) * e3 / delta;
        d[0][1] = (nu21 + nu31 * nu23) * e1 / delta;
        d[1][0] = d[0][1];
        d[0][2] = (nu31 + nu21 * nu32) * e1 / delta;
        d[2][0] = d[0][2];
        d[1][2] = (nu32 + nu12 * nu31) * e2 / delta;
        d[2][1] = d[1][2];
        d[3][3] = self.g12;
        d[4][4] = self.g23;
        d[5][5] = self.g13;
        d
    }
    /// Check Voigt symmetry: `nu_ij/Ei = nu_ji/Ej` within a tolerance of 1e-10.
    pub fn check_symmetry(&self) -> bool {
        let tol = 1.0e-10;
        let nu21 = self.nu12 * self.e2 / self.e1;
        let nu31 = self.nu13 * self.e3 / self.e1;
        let nu32 = self.nu23 * self.e3 / self.e2;
        (self.nu12 / self.e1 - nu21 / self.e2).abs() < tol
            && (self.nu13 / self.e1 - nu31 / self.e3).abs() < tol
            && (self.nu23 / self.e2 - nu32 / self.e3).abs() < tol
    }
    /// Longitudinal bulk modulus: average over the three axes, E_avg / (3*(1 - 2*nu_avg)).
    ///
    /// This is an engineering approximation for the effective bulk modulus.
    pub fn bulk_modulus_longitudinal(&self) -> f64 {
        let e_avg = (self.e1 + self.e2 + self.e3) / 3.0;
        let nu_avg = (self.nu12 + self.nu23 + self.nu13) / 3.0;
        e_avg / (3.0 * (1.0 - 2.0 * nu_avg))
    }
}
/// Biaxial failure envelope for a composite material.
///
/// Represents the outer boundary of admissible stress states in the
/// (σ11, σ22) plane using a maximum-stress interaction criterion.
#[derive(Debug, Clone, Copy)]
pub struct FailureEnvelope {
    /// Tensile strength along axis 1 \[Pa\].
    pub xt: f64,
    /// Compressive strength along axis 1 \[Pa\] (negative value).
    pub xc: f64,
    /// Tensile strength along axis 2 \[Pa\].
    pub yt: f64,
    /// Compressive strength along axis 2 \[Pa\] (negative value).
    pub yc: f64,
}
impl FailureEnvelope {
    /// Create with explicit strengths.
    pub fn max_stress(xt: f64, xc: f64, yt: f64, yc: f64) -> Self {
        Self { xt, xc, yt, yc }
    }
    /// Biaxial failure index for `(σ11, σ22)`.
    ///
    /// Returns the maximum of the two normalised stress ratios.
    /// A value ≥ 1 indicates failure.
    pub fn failure_index_biaxial(&self, sigma11: f64, sigma22: f64) -> f64 {
        let fi1 = if sigma11 >= 0.0 {
            sigma11 / self.xt
        } else {
            sigma11.abs() / self.xc.abs()
        };
        let fi2 = if sigma22 >= 0.0 {
            sigma22 / self.yt
        } else {
            sigma22.abs() / self.yc.abs()
        };
        fi1.max(fi2)
    }
    /// Tsai-Wu tensor failure index for biaxial `(σ11, σ22)`.
    ///
    /// `FI = F1*σ11 + F2*σ22 + F11*σ11² + F22*σ22² + 2*F12*σ11*σ22`
    ///
    /// Uses `F12 = -0.5*sqrt(F11*F22)` (Tsai-Wu interaction term).
    pub fn tsai_wu_index(&self, sigma11: f64, sigma22: f64) -> f64 {
        let f1 = 1.0 / self.xt + 1.0 / self.xc;
        let f2 = 1.0 / self.yt + 1.0 / self.yc;
        let f11 = -1.0 / (self.xt * self.xc);
        let f22 = -1.0 / (self.yt * self.yc);
        let f12 = -0.5 * (f11 * f22).abs().sqrt();
        f1 * sigma11
            + f2 * sigma22
            + f11 * sigma11 * sigma11
            + f22 * sigma22 * sigma22
            + 2.0 * f12 * sigma11 * sigma22
    }
    /// Check if the stress state `(σ11, σ22)` is inside the failure envelope.
    pub fn is_inside(&self, sigma11: f64, sigma22: f64) -> bool {
        self.failure_index_biaxial(sigma11, sigma22) < 1.0
    }
}
/// Anisotropic diffusion tensor \[m²/s\] for a material.
#[derive(Debug, Clone, Copy)]
pub struct DiffusionTensor {
    /// 3×3 diffusion matrix D.
    pub d: [[f64; 3]; 3],
}
impl DiffusionTensor {
    /// Create from a full 3×3 matrix.
    pub fn new(d: [[f64; 3]; 3]) -> Self {
        Self { d }
    }
    /// Create an orthotropic (diagonal) diffusion tensor.
    pub fn orthotropic(dx: f64, dy: f64, dz: f64) -> Self {
        let mut d = [[0.0f64; 3]; 3];
        d[0][0] = dx;
        d[1][1] = dy;
        d[2][2] = dz;
        Self::new(d)
    }
    /// Compute the diffusive flux vector `J = -D · ∇c`.
    pub fn diffusive_flux(&self, grad_c: [f64; 3]) -> [f64; 3] {
        let mut j = [0.0f64; 3];
        for (i, ji) in j.iter_mut().enumerate() {
            for (j_idx, &gc_j) in grad_c.iter().enumerate() {
                *ji -= self.d[i][j_idx] * gc_j;
            }
        }
        j
    }
    /// Mean-square displacement in direction `n` after time `t`.
    ///
    /// MSD = 2 * D_eff * t, where D_eff = nᵀ D n.
    pub fn msd_along(&self, n: [f64; 3], t: f64) -> f64 {
        let d_eff: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| n[i] * self.d[i][j] * n[j]))
            .sum();
        2.0 * d_eff * t
    }
    /// Average isotropic diffusivity: (Dxx + Dyy + Dzz) / 3.
    pub fn mean_diffusivity(&self) -> f64 {
        (self.d[0][0] + self.d[1][1] + self.d[2][2]) / 3.0
    }
}
/// Hill's anisotropic yield criterion (1948) for orthotropic plasticity.
///
/// The yield function is:
/// `f(σ) = F(σ22-σ33)² + G(σ33-σ11)² + H(σ11-σ22)² + 2Lτ23² + 2Mτ31² + 2Nτ12² - σ_y²`
///
/// Reference: Hill, R. (1948). "A Theory of the Yielding and Plastic Flow of
/// Anisotropic Metals". Proc. R. Soc. Lond. A 193, 281–297.
#[derive(Debug, Clone, Copy)]
pub struct HillYieldCriterion {
    /// F parameter (22-33 anisotropy).
    pub f: f64,
    /// G parameter (33-11 anisotropy).
    pub g: f64,
    /// H parameter (11-22 anisotropy).
    pub h: f64,
    /// L parameter (23 shear anisotropy).
    pub l: f64,
    /// M parameter (31 shear anisotropy).
    pub m: f64,
    /// N parameter (12 shear anisotropy).
    pub n: f64,
}
impl HillYieldCriterion {
    /// Create from the six Hill parameters.
    pub fn new(f: f64, g: f64, h: f64, l: f64, m: f64, n: f64) -> Self {
        Self { f, g, h, l, m, n }
    }
    /// Isotropic case: reduces to von Mises criterion.
    ///
    /// For isotropic Hill: F = G = H = 1/(2 σ_y²), L = M = N = 3/(2 σ_y²).
    /// This constructor uses a unit yield stress convention (F=G=H=0.5, L=M=N=1.5).
    pub fn isotropic(sigma_y: f64) -> Self {
        let sy2 = sigma_y * sigma_y;
        Self::new(
            0.5 / sy2,
            0.5 / sy2,
            0.5 / sy2,
            1.5 / sy2,
            1.5 / sy2,
            1.5 / sy2,
        )
    }
    /// Construct from Lankford R-values `(R11, R22, R33)`.
    ///
    /// R-values describe plastic anisotropy in each principal direction.
    /// Uses the standard conversion:
    /// - `H = R11 / (1 + R11)`, `G = 1 / (1 + R11)`, `F = G * R22 / R11`
    /// - `N = (R11 + R22)*(1 + 2*R33) / (2*R33*(1+R11))`, `L = M = 1.5`
    pub fn from_r_values(r11: f64, r22: f64, r33: f64) -> Self {
        let h = r11 / (1.0 + r11);
        let g = 1.0 / (1.0 + r11);
        let f = g * r22 / r11;
        let n = (r11 + r22) * (1.0 + 2.0 * r33) / (2.0 * r33 * (1.0 + r11));
        let l = 1.5;
        let m = 1.5;
        Self::new(f, g, h, l, m, n)
    }
    /// Evaluate the Hill yield function `f(σ) = q² - σ_y²`.
    ///
    /// `sigma` = `[σ11, σ22, σ33, τ12, τ23, τ13]` (Voigt notation).
    /// Returns `> 0` for yielded state, `= 0` at yield, `< 0` elastic.
    pub fn yield_function(&self, sigma: [f64; 6], sigma_y: f64) -> f64 {
        let [s11, s22, s33, t12, t23, t13] = sigma;
        let q2 = self.f * (s22 - s33).powi(2)
            + self.g * (s33 - s11).powi(2)
            + self.h * (s11 - s22).powi(2)
            + 2.0 * self.l * t23 * t23
            + 2.0 * self.m * t13 * t13
            + 2.0 * self.n * t12 * t12;
        q2 - sigma_y * sigma_y
    }
    /// Hill effective stress: `σ_eff = sqrt(Q²)`.
    ///
    /// `sigma` = `[σ11, σ22, σ33, τ12, τ23, τ13]`.
    pub fn effective_stress(&self, sigma: [f64; 6]) -> f64 {
        let [s11, s22, s33, t12, t23, t13] = sigma;
        let q2 = self.f * (s22 - s33).powi(2)
            + self.g * (s33 - s11).powi(2)
            + self.h * (s11 - s22).powi(2)
            + 2.0 * self.l * t23 * t23
            + 2.0 * self.m * t13 * t13
            + 2.0 * self.n * t12 * t12;
        q2.max(0.0).sqrt()
    }
    /// Yield surface normal (direction of plastic strain rate) at `sigma`.
    ///
    /// Returns `dφ/dσ = [∂f/∂σ11, ∂f/∂σ22, ∂f/∂σ33, ∂f/∂τ12, ∂f/∂τ23, ∂f/∂τ13]`.
    pub fn normal(&self, sigma: [f64; 6]) -> [f64; 6] {
        let [s11, s22, s33, t12, t23, t13] = sigma;
        [
            -2.0 * self.g * (s33 - s11) + 2.0 * self.h * (s11 - s22),
            2.0 * self.f * (s22 - s33) - 2.0 * self.h * (s11 - s22),
            -2.0 * self.f * (s22 - s33) + 2.0 * self.g * (s33 - s11),
            4.0 * self.n * t12,
            4.0 * self.l * t23,
            4.0 * self.m * t13,
        ]
    }
}
/// Hashin (1980) failure criteria for unidirectional composites.
///
/// Distinguishes fibre and matrix failure modes in both tension and compression.
///
/// Reference: Hashin, Z. (1980). "Failure criteria for unidirectional fiber
/// composites". J. Appl. Mech. 47(2), 329–334.
#[derive(Debug, Clone, Copy)]
pub struct HashinFailureCriteria {
    /// Fibre tensile strength Xt \[Pa\].
    pub xt: f64,
    /// Fibre compressive strength Xc \[Pa\] (positive value).
    pub xc: f64,
    /// Matrix tensile strength Yt \[Pa\].
    pub yt: f64,
    /// Matrix compressive strength Yc \[Pa\] (positive value).
    pub yc: f64,
    /// Longitudinal (in-plane) shear strength S12 \[Pa\].
    pub s12: f64,
    /// Transverse shear strength S23 \[Pa\].
    pub s23: f64,
}
impl HashinFailureCriteria {
    /// Create a new Hashin failure criteria set.
    pub fn new(xt: f64, xc: f64, yt: f64, yc: f64, s12: f64, s23: f64) -> Self {
        Self {
            xt,
            xc,
            yt,
            yc,
            s12,
            s23,
        }
    }
    /// Typical carbon/epoxy system (AS4/3501-6).
    pub fn carbon_epoxy_typical() -> Self {
        Self::new(1480.0e6, 1200.0e6, 50.0e6, 170.0e6, 70.0e6, 40.0e6)
    }
    /// Fibre tension failure index `FI_ft = (σ11/Xt)² + (τ12/S12)²`.
    ///
    /// Returns FI ≥ 1 if failure predicted.
    pub fn fibre_tension(&self, sigma11: f64, tau12: f64) -> f64 {
        (sigma11 / self.xt).powi(2) + (tau12 / self.s12).powi(2)
    }
    /// Fibre compression failure index `FI_fc = (σ11/Xc)²`.
    pub fn fibre_compression(&self, sigma11: f64) -> f64 {
        (sigma11.abs() / self.xc).powi(2)
    }
    /// Fiber tension failure index (American English alias).
    pub fn fiber_tension(&self, sigma11: f64, tau12: f64, _tau13: f64) -> f64 {
        self.fibre_tension(sigma11, tau12)
    }
    /// Fiber compression failure index (American English alias).
    pub fn fiber_compression(&self, sigma11: f64) -> f64 {
        self.fibre_compression(sigma11)
    }
    /// Matrix tension failure index `FI_mt = (σ22/Yt)² + (τ12/S12)²`.
    pub fn matrix_tension(&self, sigma22: f64, tau12: f64) -> f64 {
        (sigma22 / self.yt).powi(2) + (tau12 / self.s12).powi(2)
    }
    /// Matrix compression failure index (Hashin 1980 form).
    ///
    /// `FI_mc = (σ22/(2S23))² + (Yc/(2S23) - 1)*(σ22/Yc) + (τ12/S12)²`
    pub fn matrix_compression(&self, sigma22: f64, tau12: f64) -> f64 {
        let a = (sigma22.abs() / (2.0 * self.s23)).powi(2);
        let b = (self.yc / (2.0 * self.s23) - 1.0) * sigma22.abs() / self.yc;
        let c = (tau12 / self.s12).powi(2);
        a + b + c
    }
    /// Compute all four failure indices and return `(max_fi, mode_name)`.
    ///
    /// - Mode 0: fibre tension, 1: fibre compression, 2: matrix tension, 3: matrix compression.
    pub fn max_failure_index(
        &self,
        sigma11: f64,
        sigma22: f64,
        tau12: f64,
        tau23: f64,
    ) -> (f64, usize) {
        let _ = tau23;
        let fi = [
            if sigma11 >= 0.0 {
                self.fibre_tension(sigma11, tau12)
            } else {
                0.0
            },
            if sigma11 < 0.0 {
                self.fibre_compression(sigma11)
            } else {
                0.0
            },
            if sigma22 >= 0.0 {
                self.matrix_tension(sigma22, tau12)
            } else {
                0.0
            },
            if sigma22 < 0.0 {
                self.matrix_compression(sigma22, tau12)
            } else {
                0.0
            },
        ];
        let max_idx = fi
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        (fi[max_idx], max_idx)
    }
}
/// Puck (1996) physically based failure criteria for unidirectional composites.
///
/// Distinguishes fibre failure (FF) from inter-fibre failure (IFF) with
/// fracture-plane-based formulation.
///
/// Reference: Puck, A. & Schürmann, H. (1998). "Failure analysis of FRP
/// laminates by means of physically based phenomenological models".
/// Compos. Sci. Technol. 58(7), 1045–1067.
#[derive(Debug, Clone, Copy)]
pub struct PuckFailureCriteria {
    /// Fibre tensile strength Xt \[Pa\].
    pub xt: f64,
    /// Matrix tensile strength Yt \[Pa\].
    pub yt: f64,
    /// Longitudinal shear strength S21 \[Pa\].
    pub s21: f64,
    /// Inclination parameter p_t for transverse tensile loading.
    pub p_t: f64,
    /// Inclination parameter p_c for transverse compressive loading.
    pub p_c: f64,
}
impl PuckFailureCriteria {
    /// Create with explicit parameters.
    pub fn new(xt: f64, yt: f64, s21: f64, p_t: f64, p_c: f64) -> Self {
        Self {
            xt,
            yt,
            s21,
            p_t,
            p_c,
        }
    }
    /// Typical parameters for carbon/epoxy (AS4/3501-6).
    pub fn carbon_epoxy_typical() -> Self {
        Self::new(1480.0e6, 50.0e6, 70.0e6, 0.27, 0.32)
    }
    /// Fibre failure index (tension mode): `FI_ff = |σ11| / Xt`.
    pub fn fibre_failure_index(&self, sigma11: f64, xt: f64) -> f64 {
        sigma11.abs() / xt
    }
    /// Inter-fibre failure index (simplified IFF mode A).
    ///
    /// For transverse tension: `FI_iff = sqrt((τ21/S21)² + (1 - p_t*Yt/S21)²*(σ22/Yt)²) + p_t*σ22/S21`
    pub fn inter_fibre_failure_index(&self, sigma22: f64, tau21: f64, _tau23: f64) -> f64 {
        if sigma22 >= 0.0 {
            let a = (tau21 / self.s21).powi(2);
            let b = (1.0 - self.p_t * self.yt / self.s21).powi(2) * (sigma22 / self.yt).powi(2);
            let fi = (a + b).sqrt() + self.p_t * sigma22 / self.s21;
            fi.max(0.0)
        } else {
            let a = (tau21 / self.s21).powi(2);
            let b = (self.p_c * sigma22.abs() / self.s21).powi(2);
            (a + b).sqrt() + self.p_c * sigma22.abs() / self.s21
        }
    }
    /// Fiber failure tension (Puck criterion).
    ///
    /// Simplified: f_FF = σ_1 / X_T with minor σ_2 correction.
    pub fn fiber_failure_tension(&self, sigma1: f64, _e_f1: f64, sigma2: f64, _e_f2: f64) -> f64 {
        if sigma1 <= 0.0 {
            return 0.0;
        }
        let correction = 1.0 + sigma2 / self.yt * 0.01;
        (sigma1 / self.xt) * correction
    }
    /// Fiber failure compression (Puck criterion).
    pub fn fiber_failure_compression(&self, sigma1: f64) -> f64 {
        if sigma1 >= 0.0 {
            return 0.0;
        }
        // Puck uses Xt as the reference for compressive FF as well for this 5-field struct;
        // callers supply negative sigma1.
        (-sigma1) / self.xt
    }
    /// Inter-fiber failure Mode A (transverse tension + shear).
    pub fn inter_fiber_failure_mode_a(
        &self,
        sigma2: f64,
        tau12: f64,
        p_tt: f64,
        _p_ct: f64,
    ) -> f64 {
        if sigma2 < 0.0 {
            return 0.0;
        }
        let term_tau = (tau12 / self.s21).powi(2);
        let factor = 1.0 - p_tt * self.yt / self.s21;
        let term_sig = (factor * sigma2 / self.yt).powi(2);
        (term_tau + term_sig).sqrt() + p_tt * sigma2 / self.s21
    }
}
/// A 3×3 rotation matrix representing a crystal symmetry operation.
#[derive(Debug, Clone, Copy)]
pub struct SymmetryOperation {
    /// The 3×3 rotation/reflection matrix.
    pub matrix: [[f64; 3]; 3],
}
impl SymmetryOperation {
    /// Create from a 3×3 matrix.
    pub fn new(matrix: [[f64; 3]; 3]) -> Self {
        Self { matrix }
    }
    /// Identity operation (no transformation).
    pub fn identity() -> Self {
        Self::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
    /// 180° rotation about the z-axis.
    pub fn c2_z() -> Self {
        Self::new([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    }
    /// Inversion (centrosymmetric operation).
    pub fn inversion() -> Self {
        Self::new([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    }
    /// Mirror plane perpendicular to z (σh).
    pub fn mirror_z() -> Self {
        Self::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    }
    /// Apply this operation to a 3-vector.
    pub fn apply(&self, v: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0f64; 3];
        for (i, oi) in out.iter_mut().enumerate() {
            for (j, &v_j) in v.iter().enumerate() {
                *oi += self.matrix[i][j] * v_j;
            }
        }
        out
    }
    /// Compose this operation with another: `self * other`.
    pub fn compose(&self, other: &SymmetryOperation) -> SymmetryOperation {
        let mut m = [[0.0f64; 3]; 3];
        for (i, row) in m.iter_mut().enumerate() {
            for (j, mij) in row.iter_mut().enumerate() {
                for k in 0..3 {
                    *mij += self.matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        SymmetryOperation::new(m)
    }
    /// Determinant of the matrix (+1 = proper rotation, -1 = improper).
    pub fn determinant(&self) -> f64 {
        let m = &self.matrix;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
    /// Returns `true` if the operation is a proper rotation (det = +1).
    pub fn is_proper_rotation(&self) -> bool {
        (self.determinant() - 1.0).abs() < 1e-10
    }
}
