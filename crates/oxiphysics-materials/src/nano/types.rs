//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// Models the matrix-fiber interfacial transition zone (ITZ).
///
/// Accounts for debonding criterion and van der Waals adhesion forces.
#[derive(Debug, Clone)]
pub struct InterfacialZone {
    /// Interface fracture energy Gc (J/m²).
    pub fracture_energy: f64,
    /// Interface strength (Pa).
    pub strength: f64,
    /// van der Waals adhesion energy (J/m²).
    pub vdw_energy: f64,
    /// Interface thickness (m).
    pub thickness: f64,
    /// Debonding flag.
    pub debonded: bool,
}
impl InterfacialZone {
    /// Creates a new interfacial zone model.
    pub fn new(fracture_energy: f64, strength: f64, vdw_energy: f64, thickness: f64) -> Self {
        Self {
            fracture_energy,
            strength,
            vdw_energy,
            thickness,
            debonded: false,
        }
    }
    /// Check debonding: debond if applied traction exceeds strength.
    pub fn check_debond(&mut self, traction: f64) {
        if traction >= self.strength {
            self.debonded = true;
        }
    }
    /// van der Waals force per unit area (Lennard-Jones approximation).
    ///
    /// F/A = A_H / (6 * pi * z^3) where A_H is Hamaker constant.
    pub fn vdw_force_per_area(&self, separation_nm: f64) -> f64 {
        let a_h = 1e-19;
        let z = separation_nm * 1e-9;
        a_h / (6.0 * std::f64::consts::PI * z.powi(3))
    }
    /// Cohesive zone traction for bilinear law.
    pub fn cohesive_traction(&self, separation: f64, delta_c: f64) -> f64 {
        let delta_ratio = (separation / delta_c).clamp(0.0, 1.0);
        self.strength * (1.0 - delta_ratio)
    }
}
/// Classical Laminate Theory: computes the ABD stiffness matrix.
///
/// The ABD matrix relates mid-plane strains/curvatures to in-plane forces/moments.
#[derive(Debug, Clone)]
pub struct CompositeLaminate {
    /// Stack of plies (bottom to top).
    pub plies: Vec<Ply>,
}
impl CompositeLaminate {
    /// Creates a new composite laminate from a list of plies.
    pub fn new(plies: Vec<Ply>) -> Self {
        Self { plies }
    }
    /// Total laminate thickness.
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }
    /// Computes the A11 membrane stiffness component.
    ///
    /// A11 = Σ Qbar11_k * (z_k - z_{k-1}).
    pub fn a11(&self) -> f64 {
        self.plies
            .iter()
            .map(|p| p.transformed_q11() * p.thickness)
            .sum()
    }
    /// Computes the D11 bending stiffness component (simplified: assumes symmetric laminate).
    ///
    /// D11 = (1/12) * Σ Qbar11_k * t_k^3.
    pub fn d11(&self) -> f64 {
        self.plies
            .iter()
            .map(|p| p.transformed_q11() * p.thickness.powi(3) / 12.0)
            .sum()
    }
    /// Checks Tsai-Wu failure criterion for the first ply.
    ///
    /// Returns the safety factor (>1 = safe).
    pub fn tsai_wu_safety_factor(&self, n11: f64, xt: f64, xc: f64) -> f64 {
        let f1 = 1.0 / xt - 1.0 / xc;
        let f11 = 1.0 / (xt * xc);
        let a_coef = f11 * n11 * n11;
        let b_coef = f1 * n11;
        if a_coef + b_coef <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / (a_coef + b_coef)
    }
}
/// Virial stress tensor computation bridging molecular to continuum mechanics.
pub struct VirialStress {
    /// Volume of the representative volume element (m³).
    pub volume_m3: f64,
    /// Temperature (K) for kinetic contribution.
    pub temperature: f64,
    /// Number of atoms in the RVE.
    pub n_atoms: usize,
    /// Atom mass (kg).
    pub atom_mass: f64,
}
impl VirialStress {
    /// Silicon RVE for nanoscale stress calculation.
    pub fn silicon_rve(side_nm: f64) -> Self {
        let side_m = side_nm * 1.0e-9;
        VirialStress {
            volume_m3: side_m.powi(3),
            temperature: 300.0,
            n_atoms: (side_nm / 0.543 * 8.0).round() as usize,
            atom_mass: 28.0 * 1.660e-27,
        }
    }
    /// Virial kinetic stress (Pa) = N*m*kT / V per component.
    pub fn kinetic_stress_pa(&self) -> f64 {
        const KB: f64 = 1.380649e-23;
        self.n_atoms as f64 * self.atom_mass * KB * self.temperature / self.volume_m3
    }
    /// Virial potential stress (Pa) from pair force-distance list.
    /// `pairs`: list of (force_x: f64, force_y: f64, force_z: f64, dx: f64, dy: f64, dz: f64).
    pub fn potential_stress_xx(&self, pairs: &[(f64, f64, f64, f64, f64, f64)]) -> f64 {
        let sum: f64 = pairs
            .iter()
            .map(|(fx, _fy, _fz, rx, _ry, _rz)| fx * rx)
            .sum();
        -sum / self.volume_m3
    }
    /// Cauchy-Born rule: mapping atomistic strain to continuum deformation gradient.
    /// Given deformation gradient F\[3x3\] as flat array \[F11,F12,F13,F21,...\],
    /// returns the Green-Lagrange strain E = 0.5*(F^T*F - I) as 6-component Voigt.
    #[allow(clippy::too_many_arguments)]
    pub fn cauchy_born_strain(f: &[f64; 9]) -> [f64; 6] {
        let ftf11 = f[0] * f[0] + f[3] * f[3] + f[6] * f[6];
        let ftf22 = f[1] * f[1] + f[4] * f[4] + f[7] * f[7];
        let ftf33 = f[2] * f[2] + f[5] * f[5] + f[8] * f[8];
        let ftf12 = f[0] * f[1] + f[3] * f[4] + f[6] * f[7];
        let ftf13 = f[0] * f[2] + f[3] * f[5] + f[6] * f[8];
        let ftf23 = f[1] * f[2] + f[4] * f[5] + f[7] * f[8];
        [
            0.5 * (ftf11 - 1.0),
            0.5 * (ftf22 - 1.0),
            0.5 * (ftf33 - 1.0),
            ftf12,
            ftf13,
            ftf23,
        ]
    }
    /// Atomistic-to-continuum stress via Hardy formulation (simplified).
    pub fn hardy_stress_xx(&self, velocities: &[[f64; 3]], forces: &[[f64; 3]]) -> f64 {
        const KB: f64 = 1.380649e-23;
        let kinetic: f64 = velocities
            .iter()
            .map(|v| self.atom_mass * v[0] * v[0])
            .sum();
        let potential: f64 = forces
            .iter()
            .enumerate()
            .map(|(i, f_arr)| {
                if i < velocities.len() {
                    f_arr[0] * velocities[i][0] * 1.0e-9
                } else {
                    0.0
                }
            })
            .sum();
        let _ = KB;
        -(kinetic + potential) / self.volume_m3
    }
}
/// Quantum dot confinement energy using particle-in-a-box model.
pub struct QuantumDotConfinement {
    /// Dot radius (nm).
    pub radius_nm: f64,
    /// Bulk bandgap Eg (eV).
    pub bulk_bandgap_ev: f64,
    /// Effective electron mass ratio m*_e/m_0.
    pub eff_mass_electron: f64,
    /// Effective hole mass ratio m*_h/m_0.
    pub eff_mass_hole: f64,
    /// Dielectric constant ε_r.
    pub dielectric_constant: f64,
    /// Material designation.
    pub designation: String,
}
impl QuantumDotConfinement {
    /// Boltzmann constant (eV/K).
    #[allow(dead_code)]
    const KB_EV: f64 = 8.617e-5;
    /// Planck constant (eV·s).
    #[allow(dead_code)]
    const H_EV_S: f64 = 4.136e-15;
    /// Free electron mass (kg).
    const M0: f64 = 9.109e-31;
    /// CdSe quantum dot.
    pub fn cdse(radius_nm: f64) -> Self {
        QuantumDotConfinement {
            radius_nm,
            bulk_bandgap_ev: 1.74,
            eff_mass_electron: 0.13,
            eff_mass_hole: 0.45,
            dielectric_constant: 9.4,
            designation: "CdSe".to_string(),
        }
    }
    /// InP quantum dot.
    pub fn inp(radius_nm: f64) -> Self {
        QuantumDotConfinement {
            radius_nm,
            bulk_bandgap_ev: 1.35,
            eff_mass_electron: 0.077,
            eff_mass_hole: 0.64,
            dielectric_constant: 12.4,
            designation: "InP".to_string(),
        }
    }
    /// Confinement energy (eV) from effective mass approximation.
    /// ΔE = ℏ²π²/(2 * μ * R²) where 1/μ = 1/m*_e + 1/m*_h.
    pub fn confinement_energy_ev(&self) -> f64 {
        const HBAR_J_S: f64 = 1.0546e-34;
        const EV_TO_J: f64 = 1.602e-19;
        let r = self.radius_nm * 1.0e-9;
        let mu = 1.0 / (1.0 / self.eff_mass_electron + 1.0 / self.eff_mass_hole);
        let mu_kg = mu * Self::M0;
        let e_j = std::f64::consts::PI.powi(2) * HBAR_J_S.powi(2) / (2.0 * mu_kg * r.powi(2));
        e_j / EV_TO_J
    }
    /// Coulomb interaction correction (eV) — first-order electrostatic.
    pub fn coulomb_correction_ev(&self) -> f64 {
        const E_SQ_OVER_4PI_EPS0: f64 = 14.4;
        let r_angstrom = self.radius_nm * 10.0;
        -1.786 * E_SQ_OVER_4PI_EPS0 / (self.dielectric_constant * r_angstrom)
    }
    /// Size-dependent bandgap Eg(R) = Eg_bulk + ΔE_confinement + ΔE_coulomb.
    pub fn effective_bandgap_ev(&self) -> f64 {
        self.bulk_bandgap_ev + self.confinement_energy_ev() + self.coulomb_correction_ev()
    }
    /// Emission wavelength (nm) corresponding to the bandgap.
    pub fn emission_wavelength_nm(&self) -> f64 {
        let eg = self.effective_bandgap_ev().max(0.01);
        1239.8 / eg
    }
    /// Bohr radius aB (nm) = ε_r * a0 / μ.
    pub fn bohr_radius_nm(&self) -> f64 {
        let a0 = 0.0529;
        let mu = 1.0 / (1.0 / self.eff_mass_electron + 1.0 / self.eff_mass_hole);
        self.dielectric_constant * a0 / mu
    }
}
/// Effective stiffness of a matrix with dispersed nanoparticles.
///
/// Uses Mori-Tanaka homogenization for spherical inclusions.
#[derive(Debug, Clone)]
pub struct NanoparticleDispersion {
    /// Particle volume fraction (0 to 1).
    pub volume_fraction: f64,
    /// Particle Young's modulus (Pa).
    pub ep: f64,
    /// Matrix Young's modulus (Pa).
    pub em: f64,
    /// Particle Poisson's ratio.
    pub nu_p: f64,
    /// Matrix Poisson's ratio.
    pub nu_m: f64,
}
impl NanoparticleDispersion {
    /// Creates a new nanoparticle dispersion.
    pub fn new(volume_fraction: f64, ep: f64, em: f64, nu_p: f64, nu_m: f64) -> Self {
        Self {
            volume_fraction: volume_fraction.clamp(0.0, 1.0),
            ep,
            em,
            nu_p,
            nu_m,
        }
    }
    /// Bulk modulus K from E and ν.
    fn bulk_modulus(e: f64, nu: f64) -> f64 {
        e / (3.0 * (1.0 - 2.0 * nu))
    }
    /// Shear modulus G from E and ν.
    fn shear_modulus(e: f64, nu: f64) -> f64 {
        e / (2.0 * (1.0 + nu))
    }
    /// Effective bulk modulus from Mori-Tanaka.
    pub fn effective_bulk_modulus(&self) -> f64 {
        let km = Self::bulk_modulus(self.em, self.nu_m);
        let kp = Self::bulk_modulus(self.ep, self.nu_p);
        let gm = Self::shear_modulus(self.em, self.nu_m);
        let f = self.volume_fraction;
        let alpha = 3.0 * km / (3.0 * km + 4.0 * gm);
        let dk = kp - km;
        km + f * dk / (1.0 + (1.0 - f) * dk / (km + 4.0 * gm / 3.0) * alpha)
    }
    /// Effective Young's modulus from Mori-Tanaka bulk and shear moduli.
    pub fn effective_youngs_modulus(&self) -> f64 {
        let km = Self::bulk_modulus(self.em, self.nu_m);
        let kp = Self::bulk_modulus(self.ep, self.nu_p);
        let gm = Self::shear_modulus(self.em, self.nu_m);
        let gp = Self::shear_modulus(self.ep, self.nu_p);
        let f = self.volume_fraction;
        let alpha_k = 3.0 * km / (3.0 * km + 4.0 * gm);
        let beta_g = 6.0 * (km + 2.0 * gm) / (5.0 * (3.0 * km + 4.0 * gm));
        let k_eff = km + f * (kp - km) / (1.0 + (1.0 - f) * (kp - km) / km * alpha_k);
        let g_eff = gm + f * (gp - gm) / (1.0 + (1.0 - f) * (gp - gm) / gm * beta_g);
        9.0 * k_eff * g_eff / (3.0 * k_eff + g_eff)
    }
}
/// Graphene in-plane elastic stiffness tensor and derived properties.
pub struct GrapheneElasticStiffness {
    /// In-plane stiffness C11 (N/m, 2D modulus × thickness).
    pub c11: f64,
    /// In-plane stiffness C12 (N/m).
    pub c12: f64,
    /// In-plane stiffness C44 (N/m).
    pub c44: f64,
    /// Graphene layer thickness (nm, typically 0.335 nm for bulk conversion).
    pub layer_thickness_nm: f64,
}
impl GrapheneElasticStiffness {
    /// Perfect monolayer graphene elastic constants (from DFT/experiments).
    pub fn monolayer() -> Self {
        GrapheneElasticStiffness {
            c11: 352.0,
            c12: 60.0,
            c44: 146.0,
            layer_thickness_nm: 0.335,
        }
    }
    /// Graphene oxide (reduced GO) with ~15% defect density.
    pub fn reduced_graphene_oxide() -> Self {
        GrapheneElasticStiffness {
            c11: 250.0,
            c12: 45.0,
            c44: 100.0,
            layer_thickness_nm: 0.335,
        }
    }
    /// In-plane Young's modulus E2D (N/m) = C11 - C12²/C11.
    pub fn youngs_modulus_2d(&self) -> f64 {
        self.c11 - self.c12 * self.c12 / self.c11
    }
    /// In-plane Poisson's ratio ν = C12 / C11.
    pub fn poisson_ratio(&self) -> f64 {
        self.c12 / self.c11
    }
    /// 3D Young's modulus (GPa) using layer thickness for bulk conversion.
    pub fn youngs_modulus_3d_gpa(&self) -> f64 {
        let t = self.layer_thickness_nm * 1.0e-9;
        (self.youngs_modulus_2d() / t) / 1.0e9
    }
    /// Bending stiffness κ (eV) ≈ 1.4 eV for monolayer graphene.
    pub fn bending_stiffness_ev(&self) -> f64 {
        let t = self.layer_thickness_nm * 1.0e-9;
        let nu = self.poisson_ratio();
        let kappa_nm = self.c11 * 1e9 * t.powi(3) / (12.0 * (1.0 - nu * nu));
        kappa_nm * 6.241e18
    }
}
/// Carbon nanotube (CNT) mechanical and thermal properties.
///
/// Models single-walled CNTs with armchair, zigzag, or chiral geometry.
#[derive(Debug, Clone)]
pub struct NanotubeProperties {
    /// Chirality type.
    pub chirality: CntChirality,
    /// Chiral indices (n, m).
    pub n: u32,
    /// Chiral index m.
    pub m: u32,
    /// Tube diameter (nm).
    pub diameter_nm: f64,
    /// Young's modulus (TPa, typically ~1 TPa for SWCNTs).
    pub youngs_modulus_tpa: f64,
    /// Tensile strength (GPa).
    pub tensile_strength_gpa: f64,
    /// Thermal conductivity (W/m/K).
    pub thermal_conductivity: f64,
}
impl NanotubeProperties {
    /// Creates an armchair (n,n) CNT.
    pub fn armchair(n: u32) -> Self {
        let a_cc = 0.142;
        let diameter = a_cc * (3.0_f64.sqrt()) * n as f64 / std::f64::consts::PI;
        Self {
            chirality: CntChirality::Armchair,
            n,
            m: n,
            diameter_nm: diameter,
            youngs_modulus_tpa: 1.0,
            tensile_strength_gpa: 130.0,
            thermal_conductivity: 3500.0,
        }
    }
    /// Creates a zigzag (n, 0) CNT.
    pub fn zigzag(n: u32) -> Self {
        let a_cc = 0.142;
        let diameter = a_cc * (3.0_f64.sqrt()) * n as f64 / std::f64::consts::PI;
        Self {
            chirality: CntChirality::Zigzag,
            n,
            m: 0,
            diameter_nm: diameter,
            youngs_modulus_tpa: 0.97,
            tensile_strength_gpa: 120.0,
            thermal_conductivity: 2000.0,
        }
    }
    /// Creates a chiral (n, m) CNT.
    pub fn chiral(n: u32, m: u32) -> Self {
        let a_cc = 0.142;
        let diameter = a_cc * ((n * n + n * m + m * m) as f64).sqrt() / std::f64::consts::PI;
        Self {
            chirality: CntChirality::Chiral,
            n,
            m,
            diameter_nm: diameter,
            youngs_modulus_tpa: 1.0,
            tensile_strength_gpa: 100.0,
            thermal_conductivity: 1500.0,
        }
    }
    /// Axial stiffness: EA = E * A, where A = π * d * t (t = 0.34 nm wall thickness).
    pub fn axial_stiffness_n(&self) -> f64 {
        let wall_thickness_nm = 0.34;
        let area_nm2 = std::f64::consts::PI * self.diameter_nm * wall_thickness_nm;
        let area_m2 = area_nm2 * 1e-18;
        let e_pa = self.youngs_modulus_tpa * 1e12;
        e_pa * area_m2
    }
    /// Returns true if the CNT is metallic (armchair or n-m divisible by 3).
    pub fn is_metallic(&self) -> bool {
        match self.chirality {
            CntChirality::Armchair => true,
            _ => (self.n as i32 - self.m as i32).abs() % 3 == 0,
        }
    }
}
/// Type of carbon nanotube chirality.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CntChirality {
    /// Armchair (n, n) — metallic.
    Armchair,
    /// Zigzag (n, 0) — can be semiconducting or metallic.
    Zigzag,
    /// Chiral (n, m) — general case.
    Chiral,
}
/// Thermoelectric material: Seebeck, Peltier, and ZT.
#[derive(Debug, Clone)]
pub struct ThermoelectricMaterial {
    /// Seebeck coefficient S (V/K).
    pub seebeck: f64,
    /// Electrical conductivity σ (S/m).
    pub electrical_conductivity: f64,
    /// Thermal conductivity κ (W/m/K).
    pub thermal_conductivity: f64,
    /// Operating temperature T (K).
    pub temperature: f64,
}
impl ThermoelectricMaterial {
    /// Creates a Bi2Te3-like thermoelectric (typical room-temperature TE).
    pub fn bismuth_telluride() -> Self {
        Self {
            seebeck: 200e-6,
            electrical_conductivity: 1e5,
            thermal_conductivity: 1.5,
            temperature: 300.0,
        }
    }
    /// Power factor PF = S² * σ (W/m/K²).
    pub fn power_factor(&self) -> f64 {
        self.seebeck * self.seebeck * self.electrical_conductivity
    }
    /// Figure of merit ZT = S² * σ * T / κ (dimensionless).
    pub fn zt(&self) -> f64 {
        self.power_factor() * self.temperature / self.thermal_conductivity
    }
    /// Peltier coefficient Π = S * T (V).
    pub fn peltier_coefficient(&self) -> f64 {
        self.seebeck * self.temperature
    }
    /// Carnot efficiency at temperature difference ΔT.
    pub fn device_efficiency(&self, delta_t: f64) -> f64 {
        let zt = self.zt();
        let t_h = self.temperature + delta_t;
        let carnot = delta_t / t_h;
        let sqrt_factor = (1.0 + zt).sqrt();
        let t_ratio = self.temperature / t_h;
        carnot * (sqrt_factor - 1.0) / (sqrt_factor + t_ratio)
    }
}
/// Dislocation mechanics: Peierls stress, Orowan strengthening, dislocation density.
pub struct DislocationMechanics {
    /// Shear modulus G (GPa).
    pub shear_modulus: f64,
    /// Burgers vector magnitude b (nm).
    pub burgers_nm: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
    /// Lattice parameter a (nm).
    pub lattice_parameter_nm: f64,
    /// Misfit parameter δ for precipitate.
    pub misfit: f64,
}
impl DislocationMechanics {
    /// FCC copper dislocation mechanics.
    pub fn fcc_copper() -> Self {
        DislocationMechanics {
            shear_modulus: 42.0,
            burgers_nm: 0.256,
            nu: 0.34,
            lattice_parameter_nm: 0.362,
            misfit: 0.005,
        }
    }
    /// BCC iron dislocation mechanics.
    pub fn bcc_iron() -> Self {
        DislocationMechanics {
            shear_modulus: 82.0,
            burgers_nm: 0.248,
            nu: 0.29,
            lattice_parameter_nm: 0.286,
            misfit: 0.01,
        }
    }
    /// Peierls-Nabarro stress τ_P (GPa).
    /// τ_P = 2G/(1-ν) * exp(-2π*w/b) where w = a/(2*(1-ν)).
    pub fn peierls_nabarro_stress(&self) -> f64 {
        let a = self.lattice_parameter_nm;
        let b = self.burgers_nm;
        let w = a / (2.0 * (1.0 - self.nu));

        2.0 * self.shear_modulus / (1.0 - self.nu) * (-2.0 * std::f64::consts::PI * w / b).exp()
    }
    /// Orowan strengthening stress Δτ (MPa) for obstacles of spacing λ.
    /// Δτ = M * G * b / λ where M = Taylor factor.
    pub fn orowan_stress_mpa(&self, obstacle_spacing_nm: f64) -> f64 {
        const M: f64 = 3.06;
        M * self.shear_modulus * self.burgers_nm / obstacle_spacing_nm * 1000.0
    }
    /// Taylor hardening: Δσ = M * α * G * b * sqrt(ρ).
    pub fn taylor_hardening_mpa(&self, dislocation_density_m2: f64) -> f64 {
        const M: f64 = 3.06;
        const ALPHA: f64 = 0.3;
        M * ALPHA
            * self.shear_modulus
            * self.burgers_nm
            * 1.0e-9
            * dislocation_density_m2.sqrt()
            * 1.0e9
    }
    /// Dislocation line energy per unit length (J/m).
    pub fn line_energy(&self) -> f64 {
        let b = self.burgers_nm * 1.0e-9;
        let g = self.shear_modulus * 1.0e9;
        g * b.powi(2) / 2.0
    }
    /// Critical resolved shear stress (Schmid's law).
    pub fn schmid_crss(&self, sigma_applied: f64, phi: f64, lambda: f64) -> f64 {
        sigma_applied * phi.cos() * lambda.cos()
    }
}
/// Graphene thermal conductivity model: size, defects, substrate coupling.
pub struct NanoGrapheneThermal {
    /// Ideal bulk in-plane conductivity (W/m·K), ≈5000 for suspended.
    pub k_bulk: f64,
    /// Substrate coupling parameter Γ (GW/m³·K) per unit area.
    pub substrate_coupling: f64,
    /// Defect scattering mean free path (nm).
    pub defect_mfp_nm: f64,
    /// Grain size in polycrystalline graphene (nm).
    pub grain_size_nm: f64,
}
impl NanoGrapheneThermal {
    /// Suspended graphene thermal model.
    pub fn suspended() -> Self {
        NanoGrapheneThermal {
            k_bulk: 5000.0,
            substrate_coupling: 0.0,
            defect_mfp_nm: 1000.0,
            grain_size_nm: 1000.0,
        }
    }
    /// Graphene on SiO₂ substrate.
    pub fn on_sio2() -> Self {
        NanoGrapheneThermal {
            k_bulk: 5000.0,
            substrate_coupling: 50.0,
            defect_mfp_nm: 500.0,
            grain_size_nm: 500.0,
        }
    }
    /// Effective phonon mean free path (nm) with defect and grain boundary scattering.
    pub fn effective_mfp_nm(&self, intrinsic_mfp: f64) -> f64 {
        1.0 / (1.0 / intrinsic_mfp + 1.0 / self.defect_mfp_nm + 1.0 / self.grain_size_nm)
    }
    /// Effective in-plane thermal conductivity (W/m·K).
    pub fn effective_conductivity(&self, intrinsic_mfp: f64) -> f64 {
        let mfp_eff = self.effective_mfp_nm(intrinsic_mfp);
        self.k_bulk * mfp_eff / intrinsic_mfp
    }
    /// Out-of-plane thermal resistance due to substrate (K·m²/W) for film thickness hf (nm).
    pub fn substrate_thermal_resistance(&self, hf_nm: f64) -> f64 {
        if self.substrate_coupling < 1.0e-15 {
            return f64::INFINITY;
        }
        let hf = hf_nm * 1.0e-9;
        1.0 / (self.substrate_coupling * 1.0e9 * hf)
    }
}
/// Molecular crystal elastic properties from bond stiffness.
#[derive(Debug, Clone)]
pub struct MolecularCrystal {
    /// Unit cell lattice parameters \[a, b, c\] in Angstroms.
    pub lattice: [f64; 3],
    /// Bond stiffness k (N/m).
    pub bond_stiffness: f64,
    /// Number of bonds per unit cell.
    pub bonds_per_cell: u32,
    /// Unit cell volume (m³).
    pub cell_volume_m3: f64,
}
impl MolecularCrystal {
    /// Creates a new molecular crystal model.
    pub fn new(lattice: [f64; 3], bond_stiffness: f64, bonds_per_cell: u32) -> Self {
        let vol = lattice[0] * lattice[1] * lattice[2] * 1e-30;
        Self {
            lattice,
            bond_stiffness,
            bonds_per_cell,
            cell_volume_m3: vol,
        }
    }
    /// Estimated Young's modulus from bond stiffness: E ≈ k * n / V^(1/3).
    pub fn youngs_modulus(&self) -> f64 {
        let vol_13 = self.cell_volume_m3.cbrt();
        self.bond_stiffness * self.bonds_per_cell as f64 / vol_13
    }
    /// Compressibility β = V / (dV/dp) ≈ 1 / (n * k / V).
    pub fn compressibility(&self) -> f64 {
        let e = self.youngs_modulus();
        if e < 1e-15 {
            return f64::INFINITY;
        }
        1.0 / e
    }
}
/// Size-dependent material strengthening and scaling models.
///
/// Covers Hall-Petch grain size strengthening and Weibull size scaling.
#[derive(Debug, Clone)]
pub struct SizeEffect {
    /// Hall-Petch coefficient k_y (Pa·m^0.5).
    pub hall_petch_k: f64,
    /// Friction stress σ0 (Pa).
    pub friction_stress: f64,
    /// Weibull modulus m (dimensionless).
    pub weibull_m: f64,
    /// Reference volume V0 (m³).
    pub reference_volume: f64,
    /// Reference strength σ0_w at V0 (Pa).
    pub reference_strength: f64,
}
impl SizeEffect {
    /// Creates a new size effect model.
    pub fn new(
        hall_petch_k: f64,
        friction_stress: f64,
        weibull_m: f64,
        reference_volume: f64,
        reference_strength: f64,
    ) -> Self {
        Self {
            hall_petch_k,
            friction_stress,
            weibull_m,
            reference_volume,
            reference_strength,
        }
    }
    /// Hall-Petch yield strength: σ_y = σ0 + k_y / sqrt(d).
    pub fn hall_petch_strength(&self, grain_diameter_m: f64) -> f64 {
        self.friction_stress + self.hall_petch_k / grain_diameter_m.sqrt()
    }
    /// Weibull size-scaled strength: σ_v = σ0 * (V0/V)^(1/m).
    pub fn weibull_strength(&self, volume_m3: f64) -> f64 {
        self.reference_strength * (self.reference_volume / volume_m3).powf(1.0 / self.weibull_m)
    }
    /// Weibull survival probability P_s = exp(-(V/V0) * (σ/σ0)^m).
    pub fn weibull_survival_probability(&self, stress: f64, volume: f64) -> f64 {
        let ratio = stress / self.reference_strength;
        (-(volume / self.reference_volume) * ratio.powf(self.weibull_m)).exp()
    }
}
/// Freely jointed chain (FJC) model for polymer elasticity.
///
/// Provides the Langevin-based force-extension relation and Flory-Huggins
/// mixing free energy.
#[derive(Debug, Clone)]
pub struct PolymerChain {
    /// Number of statistical segments N.
    pub n_segments: u32,
    /// Kuhn segment length b (m).
    pub segment_length: f64,
    /// Persistence length l_p (m).
    pub persistence_length: f64,
    /// Temperature (K).
    pub temperature: f64,
}
impl PolymerChain {
    /// Boltzmann constant (J/K).
    const KB: f64 = 1.380_649e-23;
    /// Creates a new polymer chain model.
    pub fn new(
        n_segments: u32,
        segment_length: f64,
        persistence_length: f64,
        temperature: f64,
    ) -> Self {
        Self {
            n_segments,
            segment_length,
            persistence_length,
            temperature,
        }
    }
    /// Contour length L0 = N * b.
    pub fn contour_length(&self) -> f64 {
        self.n_segments as f64 * self.segment_length
    }
    /// Root-mean-square end-to-end distance R_rms = b * sqrt(N).
    pub fn rms_end_to_end(&self) -> f64 {
        self.segment_length * (self.n_segments as f64).sqrt()
    }
    /// Langevin entropic force-extension relation F(r) ≈ (kT/b) * L^{-1}(r/L0).
    ///
    /// Uses Pade approximation for inverse Langevin: L^{-1}(x) ≈ x(3-x^2)/(1-x^2).
    pub fn entropic_force(&self, extension: f64) -> f64 {
        let l0 = self.contour_length();
        let x = (extension / l0).clamp(-0.999, 0.999);
        let inv_lang = x * (3.0 - x * x) / (1.0 - x * x);
        Self::KB * self.temperature / self.segment_length * inv_lang
    }
    /// Flory-Huggins free energy of mixing per lattice site.
    ///
    /// ΔF_mix = kT * \[φ/N * ln(φ) + (1-φ) * ln(1-φ) + χ * φ * (1-φ)\].
    pub fn flory_huggins_free_energy(&self, phi: f64, chi: f64) -> f64 {
        let phi = phi.clamp(1e-10, 1.0 - 1e-10);
        let n = self.n_segments as f64;
        Self::KB
            * self.temperature
            * (phi / n * phi.ln() + (1.0 - phi) * (1.0 - phi).ln() + chi * phi * (1.0 - phi))
    }
}
/// Fiber-matrix composite properties via micromechanics.
///
/// Implements the rule of mixtures and Halpin-Tsai equations.
#[derive(Debug, Clone)]
pub struct FiberMatrix {
    /// Fiber volume fraction Vf (0 to 1).
    pub vf: f64,
    /// Fiber longitudinal modulus Ef (Pa).
    pub ef: f64,
    /// Matrix modulus Em (Pa).
    pub em: f64,
    /// Fiber Poisson's ratio νf.
    pub nu_f: f64,
    /// Matrix Poisson's ratio νm.
    pub nu_m: f64,
}
impl FiberMatrix {
    /// Creates a new fiber-matrix system.
    pub fn new(vf: f64, ef: f64, em: f64, nu_f: f64, nu_m: f64) -> Self {
        let vf = vf.clamp(0.0, 1.0);
        Self {
            vf,
            ef,
            em,
            nu_f,
            nu_m,
        }
    }
    /// Matrix volume fraction Vm = 1 - Vf.
    pub fn vm(&self) -> f64 {
        1.0 - self.vf
    }
    /// Longitudinal modulus E1 (rule of mixtures).
    pub fn e1(&self) -> f64 {
        self.ef * self.vf + self.em * self.vm()
    }
    /// Transverse modulus E2 (inverse rule of mixtures).
    pub fn e2_rom(&self) -> f64 {
        self.ef * self.em / (self.ef * self.vm() + self.em * self.vf)
    }
    /// Transverse modulus E2 from Halpin-Tsai equations.
    ///
    /// Uses ξ = 2 for circular fibers.
    pub fn e2_halpin_tsai(&self) -> f64 {
        let xi = 2.0;
        let eta = (self.ef / self.em - 1.0) / (self.ef / self.em + xi);
        self.em * (1.0 + xi * eta * self.vf) / (1.0 - eta * self.vf)
    }
    /// Longitudinal Poisson's ratio ν12 (rule of mixtures).
    pub fn nu12(&self) -> f64 {
        self.nu_f * self.vf + self.nu_m * self.vm()
    }
    /// Shear modulus G12 (Halpin-Tsai).
    pub fn g12_halpin_tsai(&self) -> f64 {
        let gf = self.ef / (2.0 * (1.0 + self.nu_f));
        let gm = self.em / (2.0 * (1.0 + self.nu_m));
        let xi = 1.0;
        let eta = (gf / gm - 1.0) / (gf / gm + xi);
        gm * (1.0 + xi * eta * self.vf) / (1.0 - eta * self.vf)
    }
}
/// Carbon nanotube buckling analysis under axial and bending loads.
pub struct CntBuckling {
    /// CNT diameter d (nm).
    pub diameter_nm: f64,
    /// Wall thickness t (nm), typically 0.34 nm.
    pub wall_thickness_nm: f64,
    /// Young's modulus E (TPa).
    pub youngs_modulus_tpa: f64,
    /// CNT length L (nm).
    pub length_nm: f64,
}
impl CntBuckling {
    /// Create a standard SWCNT buckling model.
    pub fn swcnt(diameter_nm: f64, length_nm: f64) -> Self {
        CntBuckling {
            diameter_nm,
            wall_thickness_nm: 0.34,
            youngs_modulus_tpa: 1.0,
            length_nm,
        }
    }
    /// Moment of inertia I (nm⁴) for thin-walled cylinder.
    pub fn moment_of_inertia(&self) -> f64 {
        let r = self.diameter_nm / 2.0;
        let t = self.wall_thickness_nm;
        std::f64::consts::PI * r.powi(3) * t
    }
    /// Cross-sectional area A (nm²).
    pub fn cross_section_area(&self) -> f64 {
        std::f64::consts::PI * self.diameter_nm * self.wall_thickness_nm
    }
    /// Euler critical buckling load P_cr (nN) for pin-pin boundary conditions.
    pub fn euler_buckling_load(&self) -> f64 {
        let e_gpa = self.youngs_modulus_tpa * 1000.0;
        let i_nm4 = self.moment_of_inertia();
        let l = self.length_nm;
        std::f64::consts::PI.powi(2) * e_gpa * i_nm4 / (l * l)
    }
    /// Shell buckling critical stress (MPa) for thin cylindrical shell.
    /// σ_cr = E*t / (r * sqrt(3*(1-ν²)))
    pub fn shell_buckling_stress(&self) -> f64 {
        let r = self.diameter_nm / 2.0;
        let t = self.wall_thickness_nm;
        let e_gpa = self.youngs_modulus_tpa * 1000.0;
        let nu: f64 = 0.19;
        e_gpa * t / (r * (3.0 * (1.0 - nu * nu)).sqrt())
    }
    /// Slenderness ratio L/r for the CNT.
    pub fn slenderness_ratio(&self) -> f64 {
        let r_gyration = (self.moment_of_inertia() / self.cross_section_area()).sqrt();
        self.length_nm / r_gyration
    }
    /// Critical strain ε_cr for axial buckling.
    pub fn critical_strain(&self) -> f64 {
        let r = self.diameter_nm / 2.0;
        let t = self.wall_thickness_nm;
        let c = (3.0 * (1.0 - 0.19_f64.powi(2))).sqrt();
        t / (r * c)
    }
}
/// Biological material model: collagen fibers, bone composite, tissue.
#[derive(Debug, Clone)]
pub struct BioMaterial {
    /// Collagen fiber volume fraction (0–1).
    pub collagen_fraction: f64,
    /// Mineral (hydroxyapatite) volume fraction (0–1).
    pub mineral_fraction: f64,
    /// Collagen fiber modulus (Pa).
    pub collagen_modulus: f64,
    /// Mineral modulus (Pa).
    pub mineral_modulus: f64,
    /// Tissue Mooney-Rivlin parameter C1 (Pa).
    pub c1: f64,
    /// Tissue Mooney-Rivlin parameter C2 (Pa).
    pub c2: f64,
    /// Remodeling rate constant (s^-1).
    pub remodeling_rate: f64,
}
impl BioMaterial {
    /// Creates a cortical bone material model.
    pub fn cortical_bone() -> Self {
        Self {
            collagen_fraction: 0.35,
            mineral_fraction: 0.45,
            collagen_modulus: 1.5e9,
            mineral_modulus: 114e9,
            c1: 1e5,
            c2: 1e4,
            remodeling_rate: 1e-8,
        }
    }
    /// Creates a soft tissue (cartilage-like) model.
    pub fn soft_tissue() -> Self {
        Self {
            collagen_fraction: 0.15,
            mineral_fraction: 0.0,
            collagen_modulus: 0.5e9,
            mineral_modulus: 0.0,
            c1: 1e4,
            c2: 5e3,
            remodeling_rate: 1e-7,
        }
    }
    /// Effective modulus via Voigt rule of mixtures.
    pub fn effective_modulus(&self) -> f64 {
        self.collagen_fraction * self.collagen_modulus
            + self.mineral_fraction * self.mineral_modulus
    }
    /// Mooney-Rivlin strain energy W = C1*(I1-3) + C2*(I2-3).
    ///
    /// I1 = λ1² + λ2² + λ3² (first invariant).
    pub fn mooney_rivlin_energy(&self, i1: f64, i2: f64) -> f64 {
        self.c1 * (i1 - 3.0) + self.c2 * (i2 - 3.0)
    }
    /// Remodeling stimulus: if strain energy W exceeds threshold, increase rate.
    pub fn remodeling_stimulus(&self, strain_energy: f64, threshold: f64) -> f64 {
        if strain_energy > threshold {
            self.remodeling_rate * (strain_energy - threshold)
        } else {
            0.0
        }
    }
}
/// Surface-to-volume ratio effects: Hall-Petch breakdown and surface stress.
pub struct SurfaceToVolumeEffects {
    /// Hall-Petch slope k_HP (MPa·μm^0.5).
    pub k_hp: f64,
    /// Friction stress σ_0 (MPa).
    pub sigma_0: f64,
    /// Hall-Petch breakdown grain size d_0 (nm).
    pub d_breakdown_nm: f64,
    /// Surface energy γ (J/m²).
    pub surface_energy: f64,
    /// Surface stress τ_s (N/m).
    pub surface_stress: f64,
    /// Atom diameter d_atom (nm).
    pub d_atom_nm: f64,
}
impl SurfaceToVolumeEffects {
    /// Nanocrystalline copper model.
    pub fn nc_copper() -> Self {
        SurfaceToVolumeEffects {
            k_hp: 145.0,
            sigma_0: 25.0,
            d_breakdown_nm: 15.0,
            surface_energy: 1.7,
            surface_stress: 1.5,
            d_atom_nm: 0.256,
        }
    }
    /// Nanocrystalline iron model.
    pub fn nc_iron() -> Self {
        SurfaceToVolumeEffects {
            k_hp: 600.0,
            sigma_0: 100.0,
            d_breakdown_nm: 20.0,
            surface_energy: 2.4,
            surface_stress: 2.0,
            d_atom_nm: 0.248,
        }
    }
    /// Hall-Petch yield strength (MPa). Below d_breakdown, returns softening estimate.
    pub fn yield_strength_mpa(&self, grain_size_nm: f64) -> f64 {
        if grain_size_nm >= self.d_breakdown_nm {
            self.sigma_0 + self.k_hp / grain_size_nm.sqrt()
        } else {
            let sigma_peak = self.sigma_0 + self.k_hp / self.d_breakdown_nm.sqrt();
            sigma_peak * (grain_size_nm / self.d_breakdown_nm)
        }
    }
    /// Surface-to-volume ratio for a sphere of diameter d (nm): S/V = 6/d.
    pub fn surface_to_volume_ratio(&self, diameter_nm: f64) -> f64 {
        6.0 / diameter_nm
    }
    /// Fraction of surface atoms in a spherical nanoparticle.
    pub fn surface_atom_fraction(&self, diameter_nm: f64) -> f64 {
        let n_surface_fraction = self.d_atom_nm / diameter_nm * 6.0;
        n_surface_fraction.min(1.0)
    }
    /// Capillary pressure (MPa) inside a spherical nanoparticle.
    pub fn capillary_pressure_mpa(&self, radius_nm: f64) -> f64 {
        2.0 * self.surface_stress / (radius_nm * 1.0e-9) / 1.0e6
    }
    /// Melting point depression ΔTm (K) via Gibbs-Thomson: ΔTm = 4*γ*Tm/(ΔHf*ρ*d).
    pub fn melting_point_depression(&self, t_m_bulk: f64, h_f: f64, rho: f64, d_nm: f64) -> f64 {
        4.0 * self.surface_energy * t_m_bulk / (h_f * rho * d_nm * 1.0e-9)
    }
}
/// Grain boundary mechanics: Read-Shockley energy, misorientation-dependent strength.
pub struct GrainBoundaryMechanics {
    /// Surface energy scale factor E_0 (J/m²).
    pub e0: f64,
    /// Critical misorientation angle θ_m (degrees) — typically 15°.
    pub theta_m: f64,
    /// Grain boundary diffusivity pre-exponential D_0gb (m²/s).
    pub d0_gb: f64,
    /// GB diffusion activation energy Q_gb (J/mol).
    pub q_gb: f64,
    /// Grain boundary width δ (nm).
    pub delta_nm: f64,
}
impl GrainBoundaryMechanics {
    /// Aluminum grain boundary model.
    pub fn aluminum_gb() -> Self {
        GrainBoundaryMechanics {
            e0: 0.32,
            theta_m: 15.0,
            d0_gb: 2.0e-7,
            q_gb: 84_000.0,
            delta_nm: 0.5,
        }
    }
    /// Nickel grain boundary model.
    pub fn nickel_gb() -> Self {
        GrainBoundaryMechanics {
            e0: 0.69,
            theta_m: 15.0,
            d0_gb: 1.7e-9,
            q_gb: 115_000.0,
            delta_nm: 0.5,
        }
    }
    /// Read-Shockley grain boundary energy (J/m²).
    /// E_gb = E_0 * θ/θ_m * (1 - ln(θ/θ_m)) for θ < θ_m.
    pub fn read_shockley_energy(&self, misorientation_deg: f64) -> f64 {
        let theta_m_rad = self.theta_m.to_radians();
        let theta_rad = misorientation_deg.to_radians().min(theta_m_rad);
        let ratio = theta_rad / theta_m_rad;
        if ratio < 1.0e-10 {
            return 0.0;
        }
        self.e0 * ratio * (1.0 - ratio.ln())
    }
    /// High-angle GB energy (≈ constant for θ > θ_m) (J/m²).
    pub fn high_angle_energy(&self) -> f64 {
        self.e0
    }
    /// Grain boundary diffusivity at temperature T (K).
    pub fn gb_diffusivity(&self, t_kelvin: f64) -> f64 {
        const R: f64 = 8.314;
        self.d0_gb * (-self.q_gb / (R * t_kelvin)).exp()
    }
    /// Effective GB diffusion flux (m²/s·m) = D_gb * δ.
    pub fn gb_diffusion_flux(&self, t_kelvin: f64) -> f64 {
        self.gb_diffusivity(t_kelvin) * self.delta_nm * 1.0e-9
    }
    /// Coble creep rate ε̇ (s⁻¹) for grain size d (m).
    pub fn coble_creep_rate(&self, sigma: f64, d_m: f64, t_kelvin: f64, omega: f64) -> f64 {
        const R: f64 = 8.314;
        let d_gb = self.gb_diffusivity(t_kelvin);
        let delta = self.delta_nm * 1.0e-9;
        148.0 * sigma * d_gb * delta * omega / (R * t_kelvin * d_m.powi(3))
    }
}
/// Thin film mechanics: biaxial stress, Stoney equation, delamination.
pub struct ThinFilmMechanics {
    /// Film Young's modulus Ef (GPa).
    pub ef: f64,
    /// Film Poisson's ratio νf.
    pub nu_f: f64,
    /// Film thickness hf (nm).
    pub hf_nm: f64,
    /// Film CTE αf (1/K).
    pub cte_film: f64,
    /// Substrate Young's modulus Es (GPa).
    pub es: f64,
    /// Substrate Poisson's ratio νs.
    pub nu_s: f64,
    /// Substrate thickness hs (mm).
    pub hs_mm: f64,
    /// Substrate CTE αs (1/K).
    pub cte_substrate: f64,
    /// Film-substrate interface fracture toughness Gc (J/m²).
    pub gc: f64,
}
impl ThinFilmMechanics {
    /// TiN film on silicon substrate.
    pub fn tin_on_silicon(hf_nm: f64) -> Self {
        ThinFilmMechanics {
            ef: 450.0,
            nu_f: 0.25,
            hf_nm,
            cte_film: 9.4e-6,
            es: 130.0,
            nu_s: 0.28,
            hs_mm: 0.725,
            cte_substrate: 2.6e-6,
            gc: 5.0,
        }
    }
    /// Cu film on SiO₂/Si for interconnect.
    pub fn cu_on_sio2(hf_nm: f64) -> Self {
        ThinFilmMechanics {
            ef: 130.0,
            nu_f: 0.34,
            hf_nm,
            cte_film: 16.5e-6,
            es: 73.0,
            nu_s: 0.17,
            hs_mm: 0.725,
            cte_substrate: 0.55e-6,
            gc: 2.0,
        }
    }
    /// Biaxial modulus M_f = Ef / (1 - νf).
    pub fn biaxial_modulus(&self) -> f64 {
        self.ef / (1.0 - self.nu_f)
    }
    /// Biaxial thermal mismatch stress σ_f (MPa) for temperature change ΔT.
    pub fn thermal_mismatch_stress_mpa(&self, delta_t: f64) -> f64 {
        let delta_cte = self.cte_film - self.cte_substrate;
        -self.biaxial_modulus() * delta_cte * delta_t * 1000.0
    }
    /// Stoney equation: substrate curvature κ (1/m) from film stress.
    pub fn stoney_curvature(&self, sigma_f_mpa: f64) -> f64 {
        let hf = self.hf_nm * 1.0e-9;
        let hs = self.hs_mm * 1.0e-3;
        let m_s = self.es / (1.0 - self.nu_s);
        let sigma_f_pa = sigma_f_mpa * 1.0e6;
        6.0 * sigma_f_pa * hf / (m_s * 1.0e9 * hs.powi(2))
    }
    /// Wafer bow (μm) from Stoney curvature for wafer radius r_wafer (mm).
    pub fn wafer_bow_um(&self, sigma_f_mpa: f64, r_wafer_mm: f64) -> f64 {
        let kappa = self.stoney_curvature(sigma_f_mpa);
        let r = r_wafer_mm * 1.0e-3;
        kappa * r.powi(2) / 2.0 * 1.0e6
    }
    /// Energy release rate G for channel cracking (Beuth solution).
    /// G ≈ Z * σ² * hf / Ef where Z ≈ 2 for channel crack.
    pub fn channel_crack_erg(&self, sigma_mpa: f64) -> f64 {
        const Z: f64 = 2.0;
        let sigma_pa = sigma_mpa * 1.0e6;
        let hf = self.hf_nm * 1.0e-9;
        let ef_pa = self.ef * 1.0e9;
        Z * sigma_pa.powi(2) * hf / ef_pa
    }
    /// Delamination driving force G_del (J/m²) for thin film buckling.
    pub fn delamination_energy(&self, sigma_mpa: f64) -> f64 {
        let sigma_pa = sigma_mpa * 1.0e6;
        let hf = self.hf_nm * 1.0e-9;
        let ef_pa = self.ef * 1.0e9;
        (1.0 - self.nu_f) * sigma_pa.powi(2) * hf / ef_pa
    }
}
/// Nanoindentation analysis using Oliver-Pharr method with pile-up correction.
pub struct NanoindentationOliverPharr {
    /// Indenter tip radius R (nm) for Berkovich tip.
    pub tip_radius_nm: f64,
    /// Indenter elastic modulus E_i (GPa), diamond = 1141 GPa.
    pub e_indenter: f64,
    /// Indenter Poisson's ratio ν_i.
    pub nu_indenter: f64,
    /// Berkovich geometric constant C = 24.5.
    pub c_geom: f64,
}
impl NanoindentationOliverPharr {
    /// Standard Berkovich diamond indenter.
    pub fn berkovich_diamond() -> Self {
        NanoindentationOliverPharr {
            tip_radius_nm: 50.0,
            e_indenter: 1141.0,
            nu_indenter: 0.07,
            c_geom: 24.5,
        }
    }
    /// Cube-corner indenter.
    pub fn cube_corner() -> Self {
        NanoindentationOliverPharr {
            tip_radius_nm: 40.0,
            e_indenter: 1141.0,
            nu_indenter: 0.07,
            c_geom: 2.598,
        }
    }
    /// Reduced modulus E_r from measured stiffness S and contact area A.
    /// E_r = S * sqrt(π) / (2 * sqrt(A)).
    pub fn reduced_modulus(&self, stiffness_n_per_nm: f64, contact_area_nm2: f64) -> f64 {
        stiffness_n_per_nm * std::f64::consts::PI.sqrt() / (2.0 * contact_area_nm2.sqrt())
    }
    /// Sample Young's modulus from reduced modulus.
    /// 1/E_r = (1-ν²)/E + (1-ν_i²)/E_i
    pub fn sample_youngs_modulus(&self, e_r_gpa: f64, nu_sample: f64) -> f64 {
        let inv_er = 1.0 / e_r_gpa;
        let inv_ei = (1.0 - self.nu_indenter.powi(2)) / self.e_indenter;
        (1.0 - nu_sample.powi(2)) / (inv_er - inv_ei)
    }
    /// Projected contact area from displacement (Oliver-Pharr).
    /// A = C * hc² where hc = hmax - 0.75 * Pmax / S.
    pub fn contact_area_nm2(&self, h_max_nm: f64, p_max_n: f64, stiffness_n_per_nm: f64) -> f64 {
        let h_c = h_max_nm - 0.75 * p_max_n / stiffness_n_per_nm;
        self.c_geom * h_c.max(0.0).powi(2)
    }
    /// Hardness H (GPa) = Pmax / A.
    pub fn hardness_gpa(&self, p_max_n: f64, contact_area_nm2: f64) -> f64 {
        if contact_area_nm2 < 1.0e-30 {
            return 0.0;
        }
        p_max_n / (contact_area_nm2 * 1.0e-18) / 1.0e9
    }
    /// Pile-up correction factor f_pu using elastic-perfectly-plastic assumption.
    /// f_pu = 1 + (E_r / H * 0.1) (simplified).
    pub fn pile_up_correction(&self, e_r_gpa: f64, h_gpa: f64) -> f64 {
        if h_gpa < 1.0e-10 {
            return 1.0;
        }
        1.0 + 0.1 * e_r_gpa / h_gpa
    }
}
/// Graphene sheet properties (monolayer or bilayer).
#[derive(Debug, Clone)]
pub struct GrapheneSheet {
    /// Number of layers (1 = monolayer, 2 = bilayer, etc.).
    pub n_layers: u32,
    /// In-plane Young's modulus (TPa).
    pub youngs_modulus_tpa: f64,
    /// In-plane tensile strength (GPa).
    pub tensile_strength_gpa: f64,
    /// Thermal conductivity (W/m/K).
    pub thermal_conductivity: f64,
    /// Vacancy defect fraction (0 = perfect, 1 = fully defected).
    pub vacancy_fraction: f64,
}
impl GrapheneSheet {
    /// Creates a perfect monolayer graphene sheet.
    pub fn monolayer() -> Self {
        Self {
            n_layers: 1,
            youngs_modulus_tpa: 1.0,
            tensile_strength_gpa: 130.0,
            thermal_conductivity: 5000.0,
            vacancy_fraction: 0.0,
        }
    }
    /// Creates a bilayer graphene sheet (AB stacking).
    pub fn bilayer() -> Self {
        Self {
            n_layers: 2,
            youngs_modulus_tpa: 0.98,
            tensile_strength_gpa: 120.0,
            thermal_conductivity: 4000.0,
            vacancy_fraction: 0.0,
        }
    }
    /// Effective Young's modulus accounting for vacancy defects.
    ///
    /// Uses a simple rule: E_eff = E0 * (1 - alpha * vacancy_fraction)^2.
    pub fn effective_youngs_modulus(&self) -> f64 {
        let alpha = 2.5;
        let frac = self.vacancy_fraction.clamp(0.0, 1.0);
        self.youngs_modulus_tpa * (1.0 - alpha * frac).max(0.0).powi(2)
    }
    /// Fracture strain reduction from vacancies (empirical model).
    pub fn fracture_strain(&self) -> f64 {
        let base_strain = 0.25;
        base_strain * (1.0 - self.vacancy_fraction.clamp(0.0, 1.0)).powi(2)
    }
    /// Thermal conductivity (simple layer scaling).
    pub fn effective_thermal_conductivity(&self) -> f64 {
        self.thermal_conductivity / (1.0 + 0.1 * (self.n_layers as f64 - 1.0))
    }
}
/// A single ply in a composite laminate.
#[derive(Debug, Clone)]
pub struct Ply {
    /// Ply thickness (m).
    pub thickness: f64,
    /// Fiber orientation angle (radians).
    pub angle: f64,
    /// Longitudinal Young's modulus E1 (Pa).
    pub e1: f64,
    /// Transverse Young's modulus E2 (Pa).
    pub e2: f64,
    /// In-plane shear modulus G12 (Pa).
    pub g12: f64,
    /// Major Poisson's ratio ν12.
    pub nu12: f64,
}
impl Ply {
    /// Creates a new ply.
    pub fn new(thickness: f64, angle_deg: f64, e1: f64, e2: f64, g12: f64, nu12: f64) -> Self {
        Self {
            thickness,
            angle: angle_deg.to_radians(),
            e1,
            e2,
            g12,
            nu12,
        }
    }
    /// Minor Poisson's ratio ν21 = ν12 * E2 / E1.
    pub fn nu21(&self) -> f64 {
        self.nu12 * self.e2 / self.e1
    }
    /// Reduced stiffness matrix Q11, Q12, Q22, Q66 in principal axes.
    pub fn reduced_stiffness(&self) -> [f64; 4] {
        let denom = 1.0 - self.nu12 * self.nu21();
        let q11 = self.e1 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q22 = self.e2 / denom;
        let q66 = self.g12;
        [q11, q12, q22, q66]
    }
    /// Transformed stiffness Qbar11 in the laminate frame.
    pub fn transformed_q11(&self) -> f64 {
        let [q11, q12, q22, q66] = self.reduced_stiffness();
        let c = self.angle.cos();
        let s = self.angle.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        q11 * c4 + 2.0 * (q12 + 2.0 * q66) * s2 * c2 + q22 * s4
    }
}
/// Nanoscale phonon heat transport: Fourier breakdown, Boltzmann transport.
pub struct PhononTransport {
    /// Phonon mean free path Λ (nm) at 300 K.
    pub mean_free_path_nm: f64,
    /// Debye temperature θ_D (K).
    pub debye_temp: f64,
    /// Room-temperature thermal conductivity k_bulk (W/m·K).
    pub k_bulk: f64,
    /// Phonon group velocity v_g (m/s).
    pub v_group: f64,
    /// Volumetric heat capacity ρ*Cv (J/m³·K).
    pub rho_cv: f64,
}
impl PhononTransport {
    /// Silicon phonon transport properties.
    pub fn silicon() -> Self {
        PhononTransport {
            mean_free_path_nm: 300.0,
            debye_temp: 645.0,
            k_bulk: 148.0,
            v_group: 6400.0,
            rho_cv: 1.63e6,
        }
    }
    /// Germanium phonon transport properties.
    pub fn germanium() -> Self {
        PhononTransport {
            mean_free_path_nm: 200.0,
            debye_temp: 360.0,
            k_bulk: 60.0,
            v_group: 3900.0,
            rho_cv: 1.66e6,
        }
    }
    /// Knudsen number Kn = Λ / L for characteristic length L (nm).
    pub fn knudsen_number(&self, length_nm: f64) -> f64 {
        self.mean_free_path_nm / length_nm
    }
    /// Effective thermal conductivity accounting for ballistic-diffusive crossover.
    /// k_eff = k_bulk / (1 + Kn) (Matthiessen's rule approximation).
    pub fn effective_conductivity(&self, length_nm: f64) -> f64 {
        let kn = self.knudsen_number(length_nm);
        self.k_bulk / (1.0 + kn)
    }
    /// Interface thermal resistance (Kapitza resistance) R_K (m²·K/W).
    /// Using diffuse mismatch model approximation.
    pub fn kapitza_resistance_dmm(&self, k2: f64) -> f64 {
        4.0 * (1.0 / self.k_bulk + 1.0 / k2) / self.v_group
    }
    /// Ballistic phonon heat flux (W/m²) in quasi-ballistic regime.
    /// q_bal = 0.5 * ρ*Cv * v_g * ΔT.
    pub fn ballistic_heat_flux(&self, delta_t: f64) -> f64 {
        0.5 * self.rho_cv * self.v_group * delta_t
    }
    /// Temperature-dependent conductivity k(T) ∝ 1/T above Debye temperature.
    pub fn conductivity_at_temp(&self, t_k: f64) -> f64 {
        if t_k < self.debye_temp {
            self.k_bulk * (self.debye_temp / t_k).powf(0.5)
        } else {
            self.k_bulk * self.debye_temp / t_k
        }
    }
}
