//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

use std::f64::consts::PI;

/// A water molecule with positions for O, H1, H2 (and optionally M for 4-site).
///
/// All positions are in Angstroms.
#[derive(Debug, Clone)]
pub struct WaterMolecule {
    /// Position of the oxygen atom \[x, y, z\] in Å.
    pub oxygen: [f64; 3],
    /// Position of the first hydrogen atom \[x, y, z\] in Å.
    pub hydrogen1: [f64; 3],
    /// Position of the second hydrogen atom \[x, y, z\] in Å.
    pub hydrogen2: [f64; 3],
    /// Position of the virtual M-site (TIP4P). None for 3-site models.
    pub m_site: Option<[f64; 3]>,
}
impl WaterMolecule {
    /// O-H bond length: 0.9572 Å (experimental geometry).
    pub fn bond_length() -> f64 {
        0.9572
    }
    /// H-O-H bond angle: 104.52° (experimental geometry).
    pub fn bond_angle_deg() -> f64 {
        104.52
    }
    /// Create a 3-site water molecule with the oxygen at the given position,
    /// using the default TIP3P geometry.
    pub fn new(oxygen: [f64; 3]) -> Self {
        Self::with_geometry(oxygen, &WaterGeometry::tip3p())
    }
    /// Create a water molecule using explicit geometry parameters.
    pub fn with_geometry(oxygen: [f64; 3], geom: &WaterGeometry) -> Self {
        let r = geom.r_oh;
        let half_angle_rad = geom.half_angle_rad();
        let h1 = [
            oxygen[0] + r * half_angle_rad.sin(),
            oxygen[1] + r * half_angle_rad.cos(),
            oxygen[2],
        ];
        let h2 = [
            oxygen[0] - r * half_angle_rad.sin(),
            oxygen[1] + r * half_angle_rad.cos(),
            oxygen[2],
        ];
        let m_site = if geom.has_m_site() {
            Some([oxygen[0], oxygen[1] + geom.r_om, oxygen[2]])
        } else {
            None
        };
        Self {
            oxygen,
            hydrogen1: h1,
            hydrogen2: h2,
            m_site,
        }
    }
    /// Create a TIP4P water molecule with the M-site placed automatically.
    pub fn tip4p(oxygen: [f64; 3]) -> Self {
        Self::with_geometry(oxygen, &WaterGeometry::tip4p())
    }
    /// Create an SPC water molecule.
    pub fn spc(oxygen: [f64; 3]) -> Self {
        Self::with_geometry(oxygen, &WaterGeometry::spc())
    }
    /// Create an SPC/E water molecule.
    pub fn spce(oxygen: [f64; 3]) -> Self {
        Self::with_geometry(oxygen, &WaterGeometry::spce())
    }
    /// TIP3P parameters for this molecule (convenience wrapper).
    pub fn tip3p() -> WaterParams {
        WaterParams::tip3p()
    }
    /// Center of mass of the water molecule.
    pub fn center_of_mass(&self, params: &WaterParams) -> [f64; 3] {
        let total = params.mass_o + 2.0 * params.mass_h;
        [
            (params.mass_o * self.oxygen[0]
                + params.mass_h * self.hydrogen1[0]
                + params.mass_h * self.hydrogen2[0])
                / total,
            (params.mass_o * self.oxygen[1]
                + params.mass_h * self.hydrogen1[1]
                + params.mass_h * self.hydrogen2[1])
                / total,
            (params.mass_o * self.oxygen[2]
                + params.mass_h * self.hydrogen1[2]
                + params.mass_h * self.hydrogen2[2])
                / total,
        ]
    }
    /// Dipole moment vector of the water molecule.
    ///
    /// For 3-site models: p = q_H*r_H1 + q_H*r_H2 + q_O*r_O
    /// For 4-site models: p = q_H*r_H1 + q_H*r_H2 + q_M*r_M  (q_O = 0)
    ///
    /// Units: e*Å (multiply by 4.803 to convert to Debye).
    pub fn dipole_moment(&self, params: &WaterParams) -> [f64; 3] {
        if let Some(m) = self.m_site {
            [
                params.q_h * self.hydrogen1[0] + params.q_h * self.hydrogen2[0] + params.q_m * m[0],
                params.q_h * self.hydrogen1[1] + params.q_h * self.hydrogen2[1] + params.q_m * m[1],
                params.q_h * self.hydrogen1[2] + params.q_h * self.hydrogen2[2] + params.q_m * m[2],
            ]
        } else {
            [
                params.q_h * self.hydrogen1[0]
                    + params.q_h * self.hydrogen2[0]
                    + params.q_o * self.oxygen[0],
                params.q_h * self.hydrogen1[1]
                    + params.q_h * self.hydrogen2[1]
                    + params.q_o * self.oxygen[1],
                params.q_h * self.hydrogen1[2]
                    + params.q_h * self.hydrogen2[2]
                    + params.q_o * self.oxygen[2],
            ]
        }
    }
    /// Dipole moment magnitude in Debye.
    pub fn dipole_magnitude_debye(&self, params: &WaterParams) -> f64 {
        let p = self.dipole_moment(params);
        let mag_ea = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        mag_ea * 4.803
    }
    /// Update the M-site position based on current O, H1, H2 positions.
    /// The M-site lies along the bisector of H1-O-H2 at distance r_OM from O.
    pub fn update_m_site(&mut self, geom: &WaterGeometry) {
        if !geom.has_m_site() {
            self.m_site = None;
            return;
        }
        let mut d1 = [0.0; 3];
        let mut d2 = [0.0; 3];
        for i in 0..3 {
            d1[i] = self.hydrogen1[i] - self.oxygen[i];
            d2[i] = self.hydrogen2[i] - self.oxygen[i];
        }
        let n1 = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        let n2 = (d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]).sqrt();
        if n1 < 1e-12 || n2 < 1e-12 {
            return;
        }
        let mut bisector = [0.0; 3];
        for i in 0..3 {
            bisector[i] = d1[i] / n1 + d2[i] / n2;
        }
        let nb =
            (bisector[0] * bisector[0] + bisector[1] * bisector[1] + bisector[2] * bisector[2])
                .sqrt();
        if nb < 1e-12 {
            return;
        }
        let mut m = [0.0; 3];
        for i in 0..3 {
            m[i] = self.oxygen[i] + geom.r_om * bisector[i] / nb;
        }
        self.m_site = Some(m);
    }
    /// Verify that bond lengths match constraints to within tolerance.
    pub fn check_constraints(&self, geom: &WaterGeometry, tol: f64) -> bool {
        let r1 = dist(&self.oxygen, &self.hydrogen1);
        let r2 = dist(&self.oxygen, &self.hydrogen2);
        let r_hh = dist(&self.hydrogen1, &self.hydrogen2);
        (r1 - geom.r_oh).abs() < tol
            && (r2 - geom.r_oh).abs() < tol
            && (r_hh - geom.r_hh()).abs() < tol
    }
}
/// Geometric parameters for a water model.
#[derive(Debug, Clone)]
pub struct WaterGeometry {
    /// O-H bond length (Å).
    pub r_oh: f64,
    /// H-O-H bond angle (degrees).
    pub angle_hoh_deg: f64,
    /// Distance from O to the virtual M-site along bisector (Å). 0 for 3-site models.
    pub r_om: f64,
}
impl WaterGeometry {
    /// TIP3P geometry: r_OH = 0.9572 Å, angle = 104.52°, no M-site.
    pub fn tip3p() -> Self {
        Self {
            r_oh: 0.9572,
            angle_hoh_deg: 104.52,
            r_om: 0.0,
        }
    }
    /// SPC geometry: r_OH = 1.0 Å, angle = 109.47° (tetrahedral), no M-site.
    pub fn spc() -> Self {
        Self {
            r_oh: 1.0,
            angle_hoh_deg: 109.47,
            r_om: 0.0,
        }
    }
    /// SPC/E geometry (same as SPC): r_OH = 1.0 Å, angle = 109.47°, no M-site.
    pub fn spce() -> Self {
        Self::spc()
    }
    /// TIP4P geometry: r_OH = 0.9572 Å, angle = 104.52°, r_OM = 0.15 Å.
    pub fn tip4p() -> Self {
        Self {
            r_oh: 0.9572,
            angle_hoh_deg: 104.52,
            r_om: 0.15,
        }
    }
    /// H-H distance derived from r_OH and the bond angle.
    pub fn r_hh(&self) -> f64 {
        let half_angle = self.angle_hoh_deg / 2.0 * PI / 180.0;
        2.0 * self.r_oh * half_angle.sin()
    }
    /// Half bond angle in radians.
    pub fn half_angle_rad(&self) -> f64 {
        self.angle_hoh_deg / 2.0 * PI / 180.0
    }
    /// Whether this geometry has a virtual M-site (4-site model).
    pub fn has_m_site(&self) -> bool {
        self.r_om > 0.0
    }
}
/// Ice Ih (hexagonal ice) unit cell parameters (Angstroms).
///
/// - `a` = 4.513 Å, `c` = 7.352 Å (space group P6₃/mmc).
/// - 4 molecules per unit cell.
#[allow(dead_code)]
pub struct IceIhCell {
    /// Lattice parameter a (Å).
    pub a: f64,
    /// Lattice parameter c (Å).
    pub c: f64,
    /// Number of molecules per unit cell (4 for ice Ih).
    pub n_per_cell: usize,
}
impl IceIhCell {
    /// Default ice Ih parameters at 0 °C.
    pub fn new() -> Self {
        Self {
            a: 4.513,
            c: 7.352,
            n_per_cell: 4,
        }
    }
    /// Unit cell volume (Å³).
    pub fn cell_volume(&self) -> f64 {
        (3.0_f64.sqrt() / 2.0) * self.a * self.a * self.c
    }
    /// Mass density of ice Ih in g/cm³.
    ///
    /// Uses M_water = 18.015 g/mol and Avogadro's number N_A = 6.022e23 mol⁻¹.
    pub fn density_gcc(&self) -> f64 {
        let m_water_g = 18.015 / 6.022e23_f64;
        let n_mol = self.n_per_cell as f64;
        let vol_cm3 = self.cell_volume() * 1e-24_f64;
        n_mol * m_water_g / vol_cm3
    }
    /// Number density (molecules per Å³).
    pub fn number_density(&self) -> f64 {
        self.n_per_cell as f64 / self.cell_volume()
    }
}
/// Force field parameters for a water model.
#[derive(Debug, Clone)]
pub struct WaterParams {
    /// Model name (e.g. "TIP3P", "SPC", "SPC/E", "TIP4P").
    pub name: String,
    /// Model type enum.
    pub model_type: WaterModelType,
    /// Partial charge on oxygen (elementary charge units, e).
    pub q_o: f64,
    /// Partial charge on each hydrogen (elementary charge units, e).
    pub q_h: f64,
    /// Lennard-Jones σ for oxygen-oxygen interactions (Å).
    pub sigma_o: f64,
    /// Lennard-Jones ε for oxygen-oxygen interactions (kJ/mol).
    pub epsilon_o: f64,
    /// Mass of oxygen (atomic mass units, u).
    pub mass_o: f64,
    /// Mass of hydrogen (atomic mass units, u).
    pub mass_h: f64,
    /// Geometry for this model.
    pub geometry: WaterGeometry,
    /// Charge on M-site (for 4-site models). 0 for 3-site models.
    pub q_m: f64,
}
impl WaterParams {
    /// TIP3P water model parameters.
    ///
    /// Charges: q_O = -0.834 e, q_H = +0.417 e.
    /// LJ: sigma = 3.1507 Å, epsilon = 0.636 kJ/mol.
    pub fn tip3p() -> Self {
        Self {
            name: "TIP3P".to_string(),
            model_type: WaterModelType::Tip3p,
            q_o: -0.834,
            q_h: 0.417,
            sigma_o: 3.1507,
            epsilon_o: 0.636,
            mass_o: 15.999,
            mass_h: 1.008,
            geometry: WaterGeometry::tip3p(),
            q_m: 0.0,
        }
    }
    /// SPC (Simple Point Charge) water model parameters.
    ///
    /// Charges: q_O = -0.82 e, q_H = +0.41 e.
    /// LJ: sigma = 3.166 Å, epsilon = 0.6502 kJ/mol.
    pub fn spc() -> Self {
        Self {
            name: "SPC".to_string(),
            model_type: WaterModelType::Spc,
            q_o: -0.82,
            q_h: 0.41,
            sigma_o: 3.166,
            epsilon_o: 0.6502,
            mass_o: 15.999,
            mass_h: 1.008,
            geometry: WaterGeometry::spc(),
            q_m: 0.0,
        }
    }
    /// SPC/E (Extended Simple Point Charge) water model parameters.
    ///
    /// Charges: q_O = -0.8476 e, q_H = +0.4238 e.
    /// LJ: sigma = 3.166 Å, epsilon = 0.6502 kJ/mol.
    pub fn spce() -> Self {
        Self {
            name: "SPC/E".to_string(),
            model_type: WaterModelType::Spce,
            q_o: -0.8476,
            q_h: 0.4238,
            sigma_o: 3.166,
            epsilon_o: 0.6502,
            mass_o: 15.999,
            mass_h: 1.008,
            geometry: WaterGeometry::spce(),
            q_m: 0.0,
        }
    }
    /// TIP4P water model parameters.
    ///
    /// In TIP4P the negative charge sits on a virtual M-site along the
    /// H-O-H bisector rather than on the oxygen.
    ///
    /// Charges: q_O = 0, q_H = +0.52 e, q_M = -1.04 e.
    /// LJ (on O): sigma = 3.1540 Å, epsilon = 0.6480 kJ/mol.
    pub fn tip4p() -> Self {
        Self {
            name: "TIP4P".to_string(),
            model_type: WaterModelType::Tip4p,
            q_o: 0.0,
            q_h: 0.52,
            sigma_o: 3.1540,
            epsilon_o: 0.6480,
            mass_o: 15.999,
            mass_h: 1.008,
            geometry: WaterGeometry::tip4p(),
            q_m: -1.04,
        }
    }
    /// Total molecular charge (should be zero for a neutral water model).
    pub fn total_charge(&self) -> f64 {
        self.q_o + 2.0 * self.q_h + self.q_m
    }
    /// Constraint bond lengths for this model: (r_OH, r_HH).
    pub fn constraint_distances(&self) -> (f64, f64) {
        (self.geometry.r_oh, self.geometry.r_hh())
    }
    /// LJ well depth in kJ/mol.
    pub fn lj_well_depth(&self) -> f64 {
        self.epsilon_o
    }
    /// LJ sigma in Å.
    pub fn lj_sigma(&self) -> f64 {
        self.sigma_o
    }
    /// Get parameters by model type.
    pub fn from_model_type(model_type: WaterModelType) -> Self {
        match model_type {
            WaterModelType::Tip3p => Self::tip3p(),
            WaterModelType::Spc => Self::spc(),
            WaterModelType::Spce => Self::spce(),
            WaterModelType::Tip4p => Self::tip4p(),
        }
    }
}
/// A TIP5P water molecule with O, H1, H2, LP1, LP2 positions (Å).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Tip5pMolecule {
    /// Oxygen position.
    pub oxygen: [f64; 3],
    /// First hydrogen.
    pub hydrogen1: [f64; 3],
    /// Second hydrogen.
    pub hydrogen2: [f64; 3],
    /// First lone-pair site.
    pub lp1: [f64; 3],
    /// Second lone-pair site.
    pub lp2: [f64; 3],
}
impl Tip5pMolecule {
    /// Build a TIP5P molecule with O at `oxygen` using default geometry.
    pub fn new(oxygen: [f64; 3]) -> Self {
        let p = Tip5pParams::new();
        let hoh_half = p.angle_hoh_deg / 2.0 * PI / 180.0;
        let lp_half = p.angle_lp_olp_deg / 2.0 * PI / 180.0;
        let h1 = [
            oxygen[0] + p.r_oh * hoh_half.sin(),
            oxygen[1] + p.r_oh * hoh_half.cos(),
            oxygen[2],
        ];
        let h2 = [
            oxygen[0] - p.r_oh * hoh_half.sin(),
            oxygen[1] + p.r_oh * hoh_half.cos(),
            oxygen[2],
        ];
        let lp1 = [
            oxygen[0] + p.r_olp * lp_half.sin(),
            oxygen[1] - p.r_olp * lp_half.cos(),
            oxygen[2],
        ];
        let lp2 = [
            oxygen[0] - p.r_olp * lp_half.sin(),
            oxygen[1] - p.r_olp * lp_half.cos(),
            oxygen[2],
        ];
        Self {
            oxygen,
            hydrogen1: h1,
            hydrogen2: h2,
            lp1,
            lp2,
        }
    }
    /// Check that all O-H bonds are at the target length within tolerance.
    pub fn check_bond_lengths(&self, tol: f64) -> bool {
        let p = Tip5pParams::new();
        let r1 = dist(&self.oxygen, &self.hydrogen1);
        let r2 = dist(&self.oxygen, &self.hydrogen2);
        (r1 - p.r_oh).abs() < tol && (r2 - p.r_oh).abs() < tol
    }
    /// O-LP1 distance.
    pub fn r_olp1(&self) -> f64 {
        dist(&self.oxygen, &self.lp1)
    }
}
/// Enum of supported water model types.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WaterModelType {
    /// TIP3P — Jorgensen et al. (1983)
    Tip3p,
    /// SPC — Berendsen et al. (1981)
    Spc,
    /// SPC/E — Berendsen et al. (1987)
    Spce,
    /// TIP4P — Jorgensen et al. (1983), 4-site model
    Tip4p,
}
/// Flexible SPC (f-SPC) water model.
///
/// Adds harmonic bond stretching and angle bending to the rigid SPC model.
///
/// Harmonic bond potential:  V_bond = (K_r/2) * (r - r_0)²
/// Harmonic angle potential: V_angle = (K_theta/2) * (theta - theta_0)²
///
/// Reference: Toukan & Rahman, Phys. Rev. B 31, 2643 (1985).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlexibleSpcParams {
    /// Equilibrium O-H bond length r_0 (Å).
    pub r_oh: f64,
    /// Equilibrium H-O-H angle theta_0 (degrees).
    pub angle_hoh_deg: f64,
    /// Bond force constant K_r (kJ/mol/Å²).
    pub k_bond: f64,
    /// Angle force constant K_theta (kJ/mol/rad²).
    pub k_angle: f64,
    /// Electrostatic charges (same as SPC).
    pub q_o: f64,
    /// Hydrogen partial charge (e).
    pub q_h: f64,
    /// LJ parameters.
    pub sigma_o: f64,
    /// LJ well depth for oxygen (kJ/mol).
    pub epsilon_o: f64,
}
impl FlexibleSpcParams {
    /// Default f-SPC parameters.
    pub fn new() -> Self {
        Self {
            r_oh: 1.0,
            angle_hoh_deg: 109.47,
            k_bond: 4637.0,
            k_angle: 383.0,
            q_o: -0.82,
            q_h: 0.41,
            sigma_o: 3.166,
            epsilon_o: 0.6502,
        }
    }
    /// Bond stretching energy for a single O-H bond.
    pub fn bond_energy(&self, r_oh: f64) -> f64 {
        let dr = r_oh - self.r_oh;
        0.5 * self.k_bond * dr * dr
    }
    /// Angle bending energy.
    pub fn angle_energy(&self, angle_rad: f64) -> f64 {
        let theta0 = self.angle_hoh_deg * PI / 180.0;
        let dtheta = angle_rad - theta0;
        0.5 * self.k_angle * dtheta * dtheta
    }
    /// Total intramolecular potential energy of a water molecule.
    pub fn intramolecular_energy(&self, mol: &WaterMolecule) -> f64 {
        let r1 = dist(&mol.oxygen, &mol.hydrogen1);
        let r2 = dist(&mol.oxygen, &mol.hydrogen2);
        let d1 = [
            mol.hydrogen1[0] - mol.oxygen[0],
            mol.hydrogen1[1] - mol.oxygen[1],
            mol.hydrogen1[2] - mol.oxygen[2],
        ];
        let d2 = [
            mol.hydrogen2[0] - mol.oxygen[0],
            mol.hydrogen2[1] - mol.oxygen[1],
            mol.hydrogen2[2] - mol.oxygen[2],
        ];
        let n1 = (d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]).sqrt();
        let n2 = (d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]).sqrt();
        let cos_theta = if n1 < 1e-12 || n2 < 1e-12 {
            -1.0
        } else {
            (d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]) / (n1 * n2)
        };
        let theta = cos_theta.clamp(-1.0, 1.0).acos();
        self.bond_energy(r1) + self.bond_energy(r2) + self.angle_energy(theta)
    }
}
/// SPC / SPC-E water-model helper: harmonic H-O-H angle energy.
///
/// Even though SPC treats water as rigid, one can evaluate the energy penalty
/// when the geometry deviates from the ideal tetrahedral angle.  This is
/// useful in flexible SPC variants or for testing.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Spc {
    /// Harmonic force constant for the H-O-H angle (kJ mol⁻¹ rad⁻²).
    pub k_theta: f64,
    /// Equilibrium H-O-H angle (radians).
    pub theta_0: f64,
}
impl Spc {
    /// SPC parameters: θ₀ = 109.47°, k_θ = 383.0 kJ mol⁻¹ rad⁻².
    #[allow(dead_code)]
    #[allow(clippy::self_named_constructors)]
    pub fn spc() -> Self {
        Self {
            k_theta: 383.0,
            theta_0: 109.47_f64.to_radians(),
        }
    }
    /// SPC/E parameters — same geometry as SPC.
    #[allow(dead_code)]
    pub fn spce() -> Self {
        Self::spc()
    }
    /// Harmonic angle energy: E = ½ k_θ (θ − θ₀)².
    ///
    /// # Arguments
    /// * `r_o`  — oxygen position (Å)
    /// * `r_h1` — hydrogen-1 position (Å)
    /// * `r_h2` — hydrogen-2 position (Å)
    ///
    /// Returns energy in kJ mol⁻¹.
    #[allow(dead_code)]
    pub fn compute_angle_energy(&self, r_o: [f64; 3], r_h1: [f64; 3], r_h2: [f64; 3]) -> f64 {
        let v1 = [r_h1[0] - r_o[0], r_h1[1] - r_o[1], r_h1[2] - r_o[2]];
        let v2 = [r_h2[0] - r_o[0], r_h2[1] - r_o[1], r_h2[2] - r_o[2]];
        let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if len1 < 1e-15 || len2 < 1e-15 {
            return 0.0;
        }
        let cos_theta = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (len1 * len2);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let d = theta - self.theta_0;
        0.5 * self.k_theta * d * d
    }
    /// Compute the current H-O-H angle (radians) for the given geometry.
    #[allow(dead_code)]
    pub fn compute_angle(r_o: [f64; 3], r_h1: [f64; 3], r_h2: [f64; 3]) -> f64 {
        let v1 = [r_h1[0] - r_o[0], r_h1[1] - r_o[1], r_h1[2] - r_o[2]];
        let v2 = [r_h2[0] - r_o[0], r_h2[1] - r_o[1], r_h2[2] - r_o[2]];
        let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if len1 < 1e-15 || len2 < 1e-15 {
            return 0.0;
        }
        let cos_theta = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (len1 * len2);
        cos_theta.clamp(-1.0, 1.0).acos()
    }
}
/// TIP4P virtual-site (M-site) parameters and force redistribution.
///
/// In TIP4P the negative charge sits on a virtual M-site located along the
/// H-O-H bisector at distance `r_om` from oxygen.  During MD the force on the
/// M-site must be redistributed to the three real atoms (O, H1, H2) before
/// integrating equations of motion.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Tip4p {
    /// O-M distance along bisector (Å).  Default 0.15 Å.
    pub r_om: f64,
    /// Weighting factor γ = r_OM / r_OL where r_OL is the projection onto each O-H bond.
    pub gamma: f64,
}
impl Tip4p {
    /// Default TIP4P parameters (r_OM = 0.15 Å, HOH = 104.52°).
    #[allow(dead_code)]
    pub fn default_params() -> Self {
        let r_om = 0.15_f64;
        let r_oh = 0.9572_f64;
        let half_angle = 52.26_f64.to_radians();
        let gamma = r_om / (r_oh * half_angle.cos());
        Self { r_om, gamma }
    }
    /// Redistribute the force on the virtual M-site to the real atoms.
    ///
    /// The M-site is located at: **r_M = r_O + γ*(r_H1 + r_H2 - 2*r_O)**
    ///
    /// By the chain rule the force contributions to O, H1, H2 are:
    /// - `f_O += (1 - 2γ) * f_M`
    /// - `f_H1 += γ * f_M`
    /// - `f_H2 += γ * f_M`
    ///
    /// Returns `(delta_f_o, delta_f_h1, delta_f_h2)` — increments to add.
    #[allow(dead_code)]
    pub fn compute_virtual_site_force_redistribution(
        &self,
        f_m: [f64; 3],
    ) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let g = self.gamma;
        let one_minus_2g = 1.0 - 2.0 * g;
        let df_o = [
            one_minus_2g * f_m[0],
            one_minus_2g * f_m[1],
            one_minus_2g * f_m[2],
        ];
        let df_h1 = [g * f_m[0], g * f_m[1], g * f_m[2]];
        let df_h2 = [g * f_m[0], g * f_m[1], g * f_m[2]];
        (df_o, df_h1, df_h2)
    }
    /// Compute the position of the M-site given the three real-atom positions.
    ///
    /// **r_M = r_O + γ*(r_H1 + r_H2 - 2*r_O)**
    #[allow(dead_code)]
    pub fn m_site_position(&self, r_o: [f64; 3], r_h1: [f64; 3], r_h2: [f64; 3]) -> [f64; 3] {
        let g = self.gamma;
        [
            r_o[0] + g * (r_h1[0] + r_h2[0] - 2.0 * r_o[0]),
            r_o[1] + g * (r_h1[1] + r_h2[1] - 2.0 * r_o[1]),
            r_o[2] + g * (r_h1[2] + r_h2[2] - 2.0 * r_o[2]),
        ]
    }
    /// Verify that force redistribution conserves total force.
    ///
    /// Sum of redistributed forces must equal the original M-site force.
    #[allow(dead_code)]
    pub fn force_conservation_check(&self, f_m: [f64; 3]) -> f64 {
        let (df_o, df_h1, df_h2) = self.compute_virtual_site_force_redistribution(f_m);
        let mut max_err = 0.0_f64;
        for k in 0..3 {
            let err = (df_o[k] + df_h1[k] + df_h2[k] - f_m[k]).abs();
            if err > max_err {
                max_err = err;
            }
        }
        max_err
    }
}
/// A cluster of water molecules with hydrogen-bond network utilities.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WaterCluster {
    /// The water molecules in the cluster.
    pub molecules: Vec<WaterMolecule>,
    /// Hydrogen-bond distance cutoff O···O (Å).
    pub r_oo_cut: f64,
    /// Hydrogen-bond angle cutoff ∠D–H···A (degrees); bonds with angle < this are excluded.
    pub angle_cut_deg: f64,
}
impl WaterCluster {
    /// Create a new cluster with standard geometric thresholds.
    ///
    /// Defaults: r_OO < 3.5 Å, ∠OHO > 150°.
    #[allow(dead_code)]
    pub fn new(molecules: Vec<WaterMolecule>) -> Self {
        Self {
            molecules,
            r_oo_cut: 3.5,
            angle_cut_deg: 150.0,
        }
    }
    /// Create a cluster with custom thresholds.
    #[allow(dead_code)]
    pub fn with_thresholds(
        molecules: Vec<WaterMolecule>,
        r_oo_cut: f64,
        angle_cut_deg: f64,
    ) -> Self {
        Self {
            molecules,
            r_oo_cut,
            angle_cut_deg,
        }
    }
    /// Count the number of hydrogen bonds in the cluster using a full geometric criterion.
    ///
    /// A hydrogen bond D–H···A exists if:
    /// 1. r(O_D … O_A) < `r_oo_cut`
    /// 2. ∠(O_D – H – O_A) > `angle_cut_deg`
    ///
    /// Each qualifying D–H···A triple counts as one hydrogen bond.
    /// Returns the total H-bond count (integer cast to f64).
    #[allow(dead_code)]
    pub fn compute_hydrogen_bond_count(&self) -> f64 {
        let mols = &self.molecules;
        let n = mols.len();
        let angle_min_cos = self.angle_cut_deg.to_radians().cos();
        let mut count = 0usize;
        for i in 0..n {
            for h_pos in [mols[i].hydrogen1, mols[i].hydrogen2] {
                let o_d = mols[i].oxygen;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let o_a = mols[j].oxygen;
                    let r_oo = dist3(o_d, o_a);
                    if r_oo >= self.r_oo_cut {
                        continue;
                    }
                    let v_hd = [o_d[0] - h_pos[0], o_d[1] - h_pos[1], o_d[2] - h_pos[2]];
                    let v_ha = [o_a[0] - h_pos[0], o_a[1] - h_pos[1], o_a[2] - h_pos[2]];
                    let len_hd = (v_hd[0] * v_hd[0] + v_hd[1] * v_hd[1] + v_hd[2] * v_hd[2]).sqrt();
                    let len_ha = (v_ha[0] * v_ha[0] + v_ha[1] * v_ha[1] + v_ha[2] * v_ha[2]).sqrt();
                    if len_hd < 1e-15 || len_ha < 1e-15 {
                        continue;
                    }
                    let cos_ang = (v_hd[0] * v_ha[0] + v_hd[1] * v_ha[1] + v_hd[2] * v_ha[2])
                        / (len_hd * len_ha);
                    if cos_ang <= angle_min_cos {
                        count += 1;
                    }
                }
            }
        }
        count as f64
    }
    /// Average number of H-bonds per water molecule.
    #[allow(dead_code)]
    pub fn average_hbonds_per_molecule(&self) -> f64 {
        let n = self.molecules.len();
        if n == 0 {
            return 0.0;
        }
        self.compute_hydrogen_bond_count() / n as f64
    }
    /// Returns the adjacency list of H-bond pairs (donor molecule index, acceptor index).
    #[allow(dead_code)]
    pub fn hbond_network(&self) -> Vec<(usize, usize)> {
        let mols = &self.molecules;
        let n = mols.len();
        let angle_min_cos = self.angle_cut_deg.to_radians().cos();
        let mut pairs = Vec::new();
        for i in 0..n {
            for h_pos in [mols[i].hydrogen1, mols[i].hydrogen2] {
                let o_d = mols[i].oxygen;
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let o_a = mols[j].oxygen;
                    let r_oo = dist3(o_d, o_a);
                    if r_oo >= self.r_oo_cut {
                        continue;
                    }
                    let v_hd = [o_d[0] - h_pos[0], o_d[1] - h_pos[1], o_d[2] - h_pos[2]];
                    let v_ha = [o_a[0] - h_pos[0], o_a[1] - h_pos[1], o_a[2] - h_pos[2]];
                    let len_hd = (v_hd[0] * v_hd[0] + v_hd[1] * v_hd[1] + v_hd[2] * v_hd[2]).sqrt();
                    let len_ha = (v_ha[0] * v_ha[0] + v_ha[1] * v_ha[1] + v_ha[2] * v_ha[2]).sqrt();
                    if len_hd < 1e-15 || len_ha < 1e-15 {
                        continue;
                    }
                    let cos_ang = (v_hd[0] * v_ha[0] + v_hd[1] * v_ha[1] + v_hd[2] * v_ha[2])
                        / (len_hd * len_ha);
                    if cos_ang <= angle_min_cos {
                        pairs.push((i, j));
                    }
                }
            }
        }
        pairs
    }
}
/// Summary of computed properties for a water model.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WaterModelSummary {
    /// Model name.
    pub name: String,
    /// Dipole moment magnitude (Debye).
    pub dipole_debye: f64,
    /// H-H distance (Å).
    pub r_hh: f64,
    /// LJ σ (Å).
    pub sigma: f64,
    /// LJ ε (kJ/mol).
    pub epsilon: f64,
    /// Whether the model has a virtual site.
    pub has_virtual_site: bool,
}
/// TIP5P water model parameters.
///
/// TIP5P has 5 interaction sites: O, H1, H2, and two lone-pair (LP) sites.
/// Charges: q_O = 0, q_H = +0.241 e, q_LP = −0.241 e each.
/// LJ on oxygen: σ = 3.1200 Å, ε = 0.6694 kJ/mol.
///
/// Reference: Mahoney & Jorgensen, J. Chem. Phys. 112, 8910 (2000).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Tip5pParams {
    /// Partial charge on each hydrogen (e).
    pub q_h: f64,
    /// Partial charge on each lone-pair site (e).
    pub q_lp: f64,
    /// LJ σ for O-O (Å).
    pub sigma_o: f64,
    /// LJ ε for O-O (kJ/mol).
    pub epsilon_o: f64,
    /// O-H bond length (Å).
    pub r_oh: f64,
    /// H-O-H angle (degrees).
    pub angle_hoh_deg: f64,
    /// O-LP distance (Å).
    pub r_olp: f64,
    /// LP-O-LP angle (degrees).
    pub angle_lp_olp_deg: f64,
}
impl Tip5pParams {
    /// TIP5P parameters as published.
    pub fn new() -> Self {
        Self {
            q_h: 0.241,
            q_lp: -0.241,
            sigma_o: 3.120,
            epsilon_o: 0.6694,
            r_oh: 0.9572,
            angle_hoh_deg: 104.52,
            r_olp: 0.70,
            angle_lp_olp_deg: 109.47,
        }
    }
    /// Total charge per molecule (should be 0).
    pub fn total_charge(&self) -> f64 {
        2.0 * self.q_h + 2.0 * self.q_lp
    }
}
