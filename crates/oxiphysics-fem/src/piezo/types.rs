//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::functions::{GAUSS2_ABSCISSAE, GAUSS2_WEIGHTS};

/// Parameters for a 1D piezoelectric transducer resonance model.
pub struct TransducerResonance {
    /// Elastic stiffness c_33 (Pa).
    pub c33: f64,
    /// Piezoelectric coefficient e_33 (C/m^2).
    pub e33: f64,
    /// Dielectric permittivity epsilon_33 (F/m).
    pub eps33: f64,
    /// Mass density (kg/m^3).
    pub density: f64,
    /// Thickness (m).
    pub thickness: f64,
}
impl TransducerResonance {
    /// Create a new transducer resonance model.
    pub fn new(c33: f64, e33: f64, eps33: f64, density: f64, thickness: f64) -> Self {
        Self {
            c33,
            e33,
            eps33,
            density,
            thickness,
        }
    }
    /// Compute the uncoupled thickness resonance frequency (Hz).
    ///
    /// f_r = (1 / 2t) * sqrt(c_33 / rho)
    pub fn resonance_frequency_uncoupled(&self) -> f64 {
        if self.density <= 0.0 || self.c33 <= 0.0 {
            return 0.0;
        }
        let v_sound = (self.c33 / self.density).sqrt();
        v_sound / (2.0 * self.thickness)
    }
    /// Compute the electromechanically stiffened resonance frequency (Hz).
    ///
    /// c_33_D = c_33_E + e_33^2 / epsilon_33
    pub fn resonance_frequency_stiffened(&self) -> f64 {
        if self.density <= 0.0 || self.c33 <= 0.0 {
            return 0.0;
        }
        let c33_stiffened = if self.eps33.abs() > 1e-60 {
            self.c33 + self.e33 * self.e33 / self.eps33
        } else {
            self.c33
        };
        let v_sound = (c33_stiffened / self.density).sqrt();
        v_sound / (2.0 * self.thickness)
    }
    /// Compute the thickness-mode electromechanical coupling factor k_t.
    ///
    /// k_t^2 = e_33^2 / (c_33 * epsilon_33)
    pub fn kt_coupling_factor(&self) -> f64 {
        if self.c33.abs() < 1e-60 || self.eps33.abs() < 1e-60 {
            return 0.0;
        }
        let kt2 = self.e33 * self.e33 / (self.c33 * self.eps33);
        kt2.sqrt()
    }
    /// Compute the anti-resonance frequency (Hz).
    ///
    /// Using the approximate formula: f_a ≈ f_r * sqrt(1 + k_t^2 * (pi^2 / 8))
    pub fn anti_resonance_frequency(&self) -> f64 {
        let fr = self.resonance_frequency_uncoupled();
        let kt = self.kt_coupling_factor();
        fr * (1.0 + kt * kt * std::f64::consts::PI * std::f64::consts::PI / 8.0).sqrt()
    }
    /// Compute the motional impedance magnitude at frequency `omega` (rad/s).
    ///
    /// Simplified 1D model: |Z_m| = omega * m / sqrt((1 - omega^2/omega_r^2)^2 + eta^2)
    /// where eta is mechanical loss factor.
    pub fn motional_impedance(&self, omega: f64, loss_factor: f64) -> f64 {
        let fr = self.resonance_frequency_uncoupled();
        let omega_r = 2.0 * std::f64::consts::PI * fr;
        let area = 1e-4;
        let mass = self.density * self.thickness * area;
        let denominator_sq = (1.0 - (omega / omega_r).powi(2)).powi(2) + loss_factor * loss_factor;
        if denominator_sq < 1e-60 {
            return f64::INFINITY;
        }
        omega * mass / denominator_sq.sqrt()
    }
}
/// Result of an electromechanical FRF computation.
pub struct EmFrfResult {
    /// Excitation frequency (rad/s).
    pub omega: f64,
    /// Mechanical response magnitude.
    pub mechanical_magnitude: f64,
    /// Electrical response (voltage or charge) magnitude.
    pub electrical_magnitude: f64,
}
/// Single finite element with coupled piezoelectric physics.
pub struct PiezoElement {
    /// Piezoelectric material of this element.
    pub material: PiezoMaterial,
    /// Element volume (m^3).
    pub volume: f64,
}
impl PiezoElement {
    /// Creates a new `PiezoElement`.
    pub fn new(material: PiezoMaterial, volume: f64) -> Self {
        PiezoElement { material, volume }
    }
    /// Returns the electrical coupling force: **F** = **e** * **epsilon** (C/m^2).
    pub fn coupling_force(&self, strain: &[f64; 6]) -> [f64; 3] {
        e_times_strain(&self.material.e_matrix, strain)
    }
    /// Returns the mechanically induced strain for an applied voltage.
    pub fn induced_strain(&self, voltage: f64, thickness: f64) -> [f64; 6] {
        const D33: f64 = 400e-12;
        let eps33 = D33 * voltage / thickness;
        [0.0, 0.0, eps33, 0.0, 0.0, 0.0]
    }
    /// Compute the element coupling matrix K_ue.
    pub fn element_coupling_matrix(&self) -> [[f64; 3]; 6] {
        self.material.coupling_matrix(self.volume)
    }
    /// Compute the element dielectric stiffness K_ee.
    pub fn element_dielectric_stiffness(&self) -> [[f64; 3]; 3] {
        self.material.dielectric_stiffness(self.volume)
    }
    /// Compute the electromechanical coupling coefficient k_33.
    ///
    /// The coupling coefficient k_33 quantifies the conversion efficiency
    /// between mechanical and electrical energy along the poling direction (3).
    ///
    /// k_33^2 = e_33^2 / (c_33 * epsilon_33)
    ///
    /// where:
    /// - e_33 is the piezoelectric stress coefficient (C/m^2)
    /// - c_33 is the elastic stiffness at constant electric field (Pa)
    /// - epsilon_33 is the permittivity at constant strain (F/m)
    pub fn compute_electromechanical_coupling(&self) -> f64 {
        let e33 = self.material.e_matrix[2][2];
        let c33 = self.material.c_matrix[2][2];
        let eps33 = self.material.epsilon_matrix[2][2];
        if c33 * eps33 < 1e-60 {
            return 0.0;
        }
        let k33_sq = (e33 * e33) / (c33 * eps33);
        k33_sq.sqrt()
    }
    /// Compute the resonance frequency of a piezoelectric element.
    ///
    /// For a thin piezoelectric disk/bar vibrating in the thickness mode,
    /// the fundamental resonance frequency is:
    ///
    /// f_r = (1 / 2L) * sqrt(c_33^D / rho)
    ///
    /// where:
    /// - L is the element length (m) along the poling direction
    /// - c_33^D is the stiffness at constant electric displacement
    ///   c_33^D = c_33^E + e_33^2 / epsilon_33
    /// - rho is the mass density (kg/m^3)
    ///
    /// Returns resonance frequency in Hz.
    pub fn compute_resonance_frequency(&self, length: f64, density: f64) -> f64 {
        if length < 1e-20 || density < 1e-20 {
            return 0.0;
        }
        let c33_e = self.material.c_matrix[2][2];
        let e33 = self.material.e_matrix[2][2];
        let eps33 = self.material.epsilon_matrix[2][2];
        let c33_d = if eps33 > 1e-30 {
            c33_e + e33 * e33 / eps33
        } else {
            c33_e
        };
        (1.0 / (2.0 * length)) * (c33_d / density).sqrt()
    }
    /// Compute the effective piezoelectric charge coefficient d_33.
    ///
    /// The d_33 coefficient relates strain to applied electric field:
    /// epsilon_3 = d_33 * E_3
    ///
    /// For a piezoelectric element, d_33 is derived from:
    /// d_33 = e_33 * s_33   (where s_33 is the elastic compliance)
    ///
    /// Alternatively, from the constitutive relations:
    /// d_33 = e_33 / c_33^E
    ///
    /// Returns d_33 in C/N (= m/V).
    pub fn compute_effective_piezo_charge(&self) -> f64 {
        let e33 = self.material.e_matrix[2][2];
        let c33 = self.material.c_matrix[2][2];
        if c33.abs() < 1e-60 {
            return 0.0;
        }
        e33 / c33
    }
    /// Compute the blocked force (force generated at zero strain).
    ///
    /// F_blocked = e_33 * V * area / thickness
    ///
    /// where V is applied voltage, area is electrode area, thickness is element length.
    pub fn compute_blocked_force(&self, voltage: f64, area: f64, thickness: f64) -> f64 {
        let e33 = self.material.e_matrix[2][2];
        if thickness.abs() < 1e-20 {
            return 0.0;
        }
        e33 * voltage * area / thickness
    }
    /// Compute the free displacement (displacement at zero force).
    ///
    /// delta_free = d_33 * V
    ///
    /// where d_33 = e_33 / c_33^E and V is applied voltage.
    pub fn compute_free_displacement(&self, voltage: f64) -> f64 {
        let d33 = self.compute_effective_piezo_charge();
        d33 * voltage
    }
}
/// Piezoelectric energy harvester parameters and efficiency calculations.
pub struct PiezoHarvester {
    /// Electromechanical coupling coefficient k^2.
    pub k_sq: f64,
    /// Mechanical quality factor Q_m.
    pub q_m: f64,
    /// Resonance frequency (rad/s).
    pub omega_r: f64,
    /// Optimum electrical load resistance (Ohm).
    pub r_opt: f64,
}
impl PiezoHarvester {
    /// Create a new piezo energy harvester.
    pub fn new(k_sq: f64, q_m: f64, omega_r: f64, r_opt: f64) -> Self {
        Self {
            k_sq,
            q_m,
            omega_r,
            r_opt,
        }
    }
    /// Maximum power conversion efficiency at resonance.
    ///
    /// eta_max = k^2 * Q_m / (1 + k^2 * Q_m)  (approximate)
    pub fn max_efficiency(&self) -> f64 {
        let x = self.k_sq * self.q_m;
        if (1.0 + x).abs() < 1e-30 {
            return 0.0;
        }
        x / (1.0 + x)
    }
    /// Power output at resonance for base acceleration amplitude a_0 and seismic mass m.
    ///
    /// P_max = m^2 * a_0^2 * Q_m / (2 * omega_r)
    pub fn max_power_at_resonance(&self, mass: f64, accel_amplitude: f64) -> f64 {
        if self.omega_r.abs() < 1e-30 {
            return 0.0;
        }
        let f_sq = mass * mass * accel_amplitude * accel_amplitude;
        f_sq * self.q_m / (2.0 * self.omega_r)
    }
    /// Figure of merit: FOM = k^2 * Q_m.
    pub fn figure_of_merit(&self) -> f64 {
        self.k_sq * self.q_m
    }
    /// Bandwidth (half-power, Hz): BW = f_r / Q_m.
    pub fn bandwidth_hz(&self) -> f64 {
        let f_r = self.omega_r / (2.0 * std::f64::consts::PI);
        if self.q_m.abs() < 1e-30 {
            return f64::INFINITY;
        }
        f_r / self.q_m
    }
}
/// Piezoelectric actuator model.
///
/// Computes the force and strain produced by an applied voltage.
pub struct PiezoActuator {
    /// Actuator material.
    pub material: PiezoMaterial,
    /// Actuator area (m^2).
    pub area: f64,
    /// Actuator thickness (m).
    pub thickness: f64,
}
impl PiezoActuator {
    /// Create a new piezoelectric actuator.
    pub fn new(material: PiezoMaterial, area: f64, thickness: f64) -> Self {
        Self {
            material,
            area,
            thickness,
        }
    }
    /// Compute the free strain produced by an applied voltage.
    ///
    /// Uses d-constants (d = e * C^{-1}): epsilon_free = d * E_field
    /// Simplified: epsilon_33 = d_33 * V / t
    pub fn free_strain(&self, voltage: f64) -> [f64; 6] {
        const D33: f64 = 400e-12;
        const D31: f64 = -175e-12;
        let e_field = voltage / self.thickness;
        [D31 * e_field, D31 * e_field, D33 * e_field, 0.0, 0.0, 0.0]
    }
    /// Compute the blocking force in the poling direction.
    ///
    /// F_block = d_33 * V / t * c_33 * A (clamped actuator)
    pub fn blocking_force(&self, voltage: f64) -> f64 {
        let strain = self.free_strain(voltage);
        let c33 = self.material.c_matrix[2][2];
        strain[2] * c33 * self.area
    }
    /// Compute the electrical energy stored in the actuator.
    ///
    /// U_e = 0.5 * epsilon_s * (V/t)^2 * Volume
    pub fn stored_energy(&self, voltage: f64) -> f64 {
        let eps33 = self.material.epsilon_matrix[2][2];
        let e_field = voltage / self.thickness;
        0.5 * eps33 * e_field * e_field * self.area * self.thickness
    }
}
/// 8-node hexahedral piezoelectric finite element.
///
/// Stores nodal coordinates, material, and assembles element stiffness matrices
/// using 2×2×2 Gauss quadrature.
pub struct PiezoHex8Element {
    /// Nodal coordinates `[8][3]` in physical space (m).
    pub coords: [[f64; 3]; 8],
    /// Piezoelectric material.
    pub material: PiezoMaterial,
}
impl PiezoHex8Element {
    /// Create a new 8-node hex piezoelectric element.
    pub fn new(coords: [[f64; 3]; 8], material: PiezoMaterial) -> Self {
        Self { coords, material }
    }
    /// Assemble the mechanical stiffness matrix K_uu (24×24).
    ///
    /// Uses 2×2×2 Gauss quadrature.
    /// Returns a flat row-major 24×24 matrix.
    pub fn stiffness_matrix(&self) -> Vec<f64> {
        let mut k = vec![0.0_f64; 24 * 24];
        for &xi in &GAUSS2_ABSCISSAE {
            for &eta in &GAUSS2_ABSCISSAE {
                for &zeta in &GAUSS2_ABSCISSAE {
                    let wt = GAUSS2_WEIGHTS[0] * GAUSS2_WEIGHTS[0] * GAUSS2_WEIGHTS[0];
                    let j = hex8_jacobian(&self.coords, xi, eta, zeta);
                    let det_j = det3(&j);
                    if det_j.abs() < 1e-30 {
                        continue;
                    }
                    let j_inv = inv3(&j);
                    let dn_nat = hex8_shape_derivatives(xi, eta, zeta);
                    let dndx = hex8_dndx(&j_inv, &dn_nat);
                    let b = hex8_b_matrix(&dndx);
                    let factor = wt * det_j;
                    for i in 0..24 {
                        for l in 0..24 {
                            let mut val = 0.0;
                            for r in 0..6 {
                                for s in 0..6 {
                                    val += b[r][i] * self.material.c_matrix[r][s] * b[s][l];
                                }
                            }
                            k[i * 24 + l] += factor * val;
                        }
                    }
                }
            }
        }
        k
    }
    /// Assemble the piezoelectric coupling matrix K_ue (24×3).
    ///
    /// Relates mechanical DOFs (24) to 3 electric DOFs (Ex, Ey, Ez).
    /// Returns flat row-major 24×3 matrix.
    pub fn coupling_matrix_full(&self) -> Vec<f64> {
        let mut k_ue = vec![0.0_f64; 24 * 3];
        for &xi in &GAUSS2_ABSCISSAE {
            for &eta in &GAUSS2_ABSCISSAE {
                for &zeta in &GAUSS2_ABSCISSAE {
                    let wt = GAUSS2_WEIGHTS[0].powi(3);
                    let j = hex8_jacobian(&self.coords, xi, eta, zeta);
                    let det_j = det3(&j);
                    if det_j.abs() < 1e-30 {
                        continue;
                    }
                    let j_inv = inv3(&j);
                    let dn_nat = hex8_shape_derivatives(xi, eta, zeta);
                    let dndx = hex8_dndx(&j_inv, &dn_nat);
                    let b = hex8_b_matrix(&dndx);
                    let factor = wt * det_j;
                    for i in 0..24 {
                        for q in 0..3 {
                            let mut val = 0.0;
                            for (r, b_row) in b.iter().enumerate() {
                                val += b_row[i] * self.material.e_matrix[r][q];
                            }
                            k_ue[i * 3 + q] += factor * val;
                        }
                    }
                }
            }
        }
        k_ue
    }
    /// Compute the element volume by integrating det(J) over the element.
    pub fn volume(&self) -> f64 {
        let mut vol = 0.0;
        for &xi in &GAUSS2_ABSCISSAE {
            for &eta in &GAUSS2_ABSCISSAE {
                for &zeta in &GAUSS2_ABSCISSAE {
                    let wt = GAUSS2_WEIGHTS[0].powi(3);
                    let j = hex8_jacobian(&self.coords, xi, eta, zeta);
                    vol += wt * det3(&j);
                }
            }
        }
        vol
    }
}
/// Piezoelectric material in d-form (strain-charge form).
///
/// The d-form constitutive equations are:
/// - `[epsilon] = [s_E][sigma] + [d]^T [E]` (mechanical strain)
/// - `[D] = [d][sigma] + [epsilon_T][E]`     (electric displacement)
///
/// This is complementary to the e-form (`e = d * c_E`).
pub struct PiezoMaterialDForm {
    /// Elastic compliance matrix `[s_E]` \[6×6\] (Pa^{-1}).
    pub s_matrix: [[f64; 6]; 6],
    /// Piezoelectric strain coefficient matrix `[d]` \[3×6\] (m/V or C/N).
    pub d_matrix: [[f64; 6]; 3],
    /// Dielectric permittivity at constant stress `[epsilon_T]` \[3×3\] (F/m).
    pub epsilon_t: [[f64; 3]; 3],
    /// Human-readable name.
    pub name: String,
}
impl PiezoMaterialDForm {
    /// Approximate PZT-5A in d-form.
    pub fn pzt5a_d_form() -> Self {
        let s11 = 16.4e-12_f64;
        let s12 = -5.74e-12_f64;
        let s13 = -7.22e-12_f64;
        let s33 = 18.8e-12_f64;
        let s44 = 47.5e-12_f64;
        let s66 = 2.0 * (s11 - s12);
        #[rustfmt::skip]
        let s_matrix = [
            [s11, s12, s13, 0.0, 0.0, 0.0],
            [s12, s11, s13, 0.0, 0.0, 0.0],
            [s13, s13, s33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, s44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, s44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, s66],
        ];
        let d31 = -171e-12_f64;
        let d33 = 374e-12_f64;
        let d15 = 584e-12_f64;
        #[rustfmt::skip]
        let d_matrix: [[f64; 6]; 3] = [
            [0.0, 0.0, 0.0, d15, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, d15, 0.0],
            [d31, d31, d33, 0.0, 0.0, 0.0],
        ];
        let eps0 = 8.854187817e-12_f64;
        let eps11 = 1700.0 * eps0;
        let eps33 = 1470.0 * eps0;
        let epsilon_t = [[eps11, 0.0, 0.0], [0.0, eps11, 0.0], [0.0, 0.0, eps33]];
        PiezoMaterialDForm {
            s_matrix,
            d_matrix,
            epsilon_t,
            name: "PZT-5A (d-form)".to_string(),
        }
    }
    /// Compute strain from stress and electric field: `epsilon = s_E sigma + d^T E`.
    pub fn strain(&self, stress: &[f64; 6], e_field: &[f64; 3]) -> [f64; 6] {
        let s_sigma = mat6_vec6_mul(&self.s_matrix, stress);
        let mut dt_e = [0.0_f64; 6];
        for (j, dte_j) in dt_e.iter_mut().enumerate() {
            for (k, &ek) in e_field.iter().enumerate() {
                *dte_j += self.d_matrix[k][j] * ek;
            }
        }
        let mut eps = [0.0_f64; 6];
        for (eps_i, (&ss_i, &dte_i)) in eps.iter_mut().zip(s_sigma.iter().zip(dt_e.iter())) {
            *eps_i = ss_i + dte_i;
        }
        eps
    }
    /// Compute electric displacement: `D = d sigma + epsilon_T E`.
    pub fn electric_displacement_d(&self, stress: &[f64; 6], e_field: &[f64; 3]) -> [f64; 3] {
        let mut d = [0.0_f64; 3];
        for (i, di) in d.iter_mut().enumerate() {
            for (j, &sj) in stress.iter().enumerate() {
                *di += self.d_matrix[i][j] * sj;
            }
            for (j, &ej) in e_field.iter().enumerate() {
                *di += self.epsilon_t[i][j] * ej;
            }
        }
        d
    }
    /// Compute the piezoelectric coupling coefficient k_33 (thickness coupling).
    ///
    /// k_33 = d_33 / sqrt(s33_E * eps33_T)
    pub fn k33_coupling_factor(&self) -> f64 {
        let d33 = self.d_matrix[2][2];
        let s33 = self.s_matrix[2][2];
        let eps33 = self.epsilon_t[2][2];
        let denom = (s33 * eps33).sqrt();
        if denom.abs() < 1e-60 {
            return 0.0;
        }
        d33 / denom
    }
}
/// PVDF (polyvinylidene fluoride) piezoelectric film material constants.
///
/// PVDF is a semi-crystalline polymer with piezoelectric properties.
/// Much more flexible than PZT, suitable for wearable / large-area sensors.
pub struct PvdfFilm {
    /// d31 piezoelectric coefficient (m/V or C/N).
    pub d31: f64,
    /// d33 piezoelectric coefficient (m/V or C/N).
    pub d33: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Relative dielectric permittivity.
    pub eps_r: f64,
    /// Film thickness (m).
    pub thickness: f64,
    /// Film area (m^2).
    pub area: f64,
}
impl PvdfFilm {
    /// Typical PVDF film parameters.
    pub fn typical() -> Self {
        Self {
            d31: 23e-12,
            d33: -33e-12,
            youngs_modulus: 2.0e9,
            eps_r: 12.0,
            thickness: 28e-6,
            area: 1e-4,
        }
    }
    /// Vacuum permittivity (F/m).
    const EPS0: f64 = 8.854187817e-12;
    /// Absolute dielectric permittivity eps = eps_r * eps_0 (F/m).
    pub fn permittivity(&self) -> f64 {
        self.eps_r * Self::EPS0
    }
    /// Capacitance C = eps_r * eps_0 * A / t (F).
    pub fn capacitance(&self) -> f64 {
        self.permittivity() * self.area / self.thickness
    }
    /// Open-circuit voltage for applied stress sigma (Pa) — 1-D d31 mode.
    ///
    /// V = d31 * sigma * t / (eps_r * eps_0)
    pub fn voltage_from_stress_d31(&self, sigma: f64) -> f64 {
        let eps = self.permittivity();
        if eps.abs() < 1e-60 {
            return 0.0;
        }
        self.d31 * sigma * self.thickness / eps
    }
    /// Open-circuit voltage for applied stress sigma (Pa) — d33 mode (thickness).
    ///
    /// V = d33 * sigma * t / (eps_r * eps_0)
    pub fn voltage_from_stress_d33(&self, sigma: f64) -> f64 {
        let eps = self.permittivity();
        if eps.abs() < 1e-60 {
            return 0.0;
        }
        self.d33 * sigma * self.thickness / eps
    }
    /// Charge generated by a force F applied over the film area (d31 mode).
    ///
    /// Q = d31 * F   (C)
    pub fn charge_from_force_d31(&self, force: f64) -> f64 {
        let sigma = force / self.area;
        self.d31 * sigma * self.area
    }
    /// Actuation strain for applied voltage V (d31 mode).
    ///
    /// epsilon = d31 * E_field = d31 * V / t
    pub fn actuation_strain_d31(&self, voltage: f64) -> f64 {
        self.d31 * voltage / self.thickness
    }
    /// Actuation strain for applied voltage V (d33 mode).
    pub fn actuation_strain_d33(&self, voltage: f64) -> f64 {
        self.d33 * voltage / self.thickness
    }
    /// Electromechanical coupling factor k31.
    ///
    /// k31^2 = d31^2 * E / eps_T
    pub fn k31_coupling_factor(&self) -> f64 {
        let eps = self.permittivity();
        if eps.abs() < 1e-60 || self.youngs_modulus.abs() < 1e-60 {
            return 0.0;
        }
        let k31_sq = self.d31 * self.d31 * self.youngs_modulus / eps;
        k31_sq.sqrt()
    }
    /// Bending stiffness of the PVDF film per unit width: EI = E * t^3 / 12.
    pub fn bending_stiffness_per_unit_width(&self) -> f64 {
        self.youngs_modulus * self.thickness.powi(3) / 12.0
    }
}
/// Piezoelectric sensor model.
///
/// Computes the output voltage from a sensor patch given the strain field.
pub struct PiezoSensor {
    /// Sensor material.
    pub material: PiezoMaterial,
    /// Sensor area (m^2).
    pub area: f64,
    /// Sensor thickness (m).
    pub thickness: f64,
}
impl PiezoSensor {
    /// Create a new piezoelectric sensor.
    pub fn new(material: PiezoMaterial, area: f64, thickness: f64) -> Self {
        Self {
            material,
            area,
            thickness,
        }
    }
    /// Compute the open-circuit voltage generated by a given strain.
    ///
    /// V = -e_33 * epsilon_33 * thickness / epsilon_s_33
    /// (simplified for uniaxial case along poling direction)
    pub fn open_circuit_voltage(&self, strain: &[f64; 6]) -> f64 {
        let _e33 = self.material.e_matrix[2][2];
        let eps33 = self.material.epsilon_matrix[2][2];
        if eps33.abs() < 1e-60 {
            return 0.0;
        }
        let d3: f64 = (0..6)
            .map(|j| self.material.e_matrix[j][2] * strain[j])
            .sum();
        -d3 * self.thickness / eps33
    }
    /// Compute the charge generated by a given strain.
    ///
    /// Q = D_3 * Area = (e * epsilon)_3 * Area
    pub fn charge_output(&self, strain: &[f64; 6]) -> f64 {
        let d = e_times_strain(&self.material.e_matrix, strain);
        d[2] * self.area
    }
}
/// Piezoelectric constitutive material.
pub struct PiezoMaterial {
    /// Elastic stiffness matrix \[6x6\] in Voigt notation (Pa).
    pub c_matrix: [[f64; 6]; 6],
    /// Piezoelectric coupling matrix stored as `[6][3]` (C/m^2).
    pub e_matrix: [[f64; 3]; 6],
    /// Dielectric permittivity matrix \[3x3\] (F/m).
    pub epsilon_matrix: [[f64; 3]; 3],
    /// Human-readable material name.
    pub name: String,
}
impl PiezoMaterial {
    /// Returns approximate PZT-5A piezoelectric material constants.
    pub fn pzt5a() -> Self {
        let c11 = 121.0e9_f64;
        let c12 = 75.4e9_f64;
        let c13 = 75.2e9_f64;
        let c33 = 111.0e9_f64;
        let c44 = 21.1e9_f64;
        let c66 = (c11 - c12) / 2.0;
        #[rustfmt::skip]
        let c_matrix: [[f64; 6]; 6] = [
            [c11, c12, c13, 0.0, 0.0, 0.0],
            [c12, c11, c13, 0.0, 0.0, 0.0],
            [c13, c13, c33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, c44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, c44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, c66],
        ];
        let e31 = -5.4_f64;
        let e33 = 15.8_f64;
        let e15 = 12.3_f64;
        #[rustfmt::skip]
        let e_matrix: [[f64; 3]; 6] = [
            [0.0, 0.0, e31],
            [0.0, 0.0, e31],
            [0.0, 0.0, e33],
            [e15, 0.0, 0.0],
            [0.0, e15, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let eps0 = 8.854187817e-12_f64;
        let eps11 = 1700.0 * eps0;
        let eps33 = 1470.0 * eps0;
        #[rustfmt::skip]
        let epsilon_matrix: [[f64; 3]; 3] = [
            [eps11, 0.0, 0.0],
            [0.0, eps11, 0.0],
            [0.0, 0.0, eps33],
        ];
        PiezoMaterial {
            c_matrix,
            e_matrix,
            epsilon_matrix,
            name: "PZT-5A".to_string(),
        }
    }
    /// Returns approximate PZT-5H piezoelectric material constants.
    pub fn pzt5h() -> Self {
        let c11 = 126.0e9_f64;
        let c12 = 79.5e9_f64;
        let c13 = 84.1e9_f64;
        let c33 = 117.0e9_f64;
        let c44 = 23.0e9_f64;
        let c66 = (c11 - c12) / 2.0;
        #[rustfmt::skip]
        let c_matrix: [[f64; 6]; 6] = [
            [c11, c12, c13, 0.0, 0.0, 0.0],
            [c12, c11, c13, 0.0, 0.0, 0.0],
            [c13, c13, c33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, c44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, c44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, c66],
        ];
        let e31 = -6.5_f64;
        let e33 = 23.3_f64;
        let e15 = 17.0_f64;
        #[rustfmt::skip]
        let e_matrix: [[f64; 3]; 6] = [
            [0.0, 0.0, e31],
            [0.0, 0.0, e31],
            [0.0, 0.0, e33],
            [e15, 0.0, 0.0],
            [0.0, e15, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let eps0 = 8.854187817e-12_f64;
        let eps11 = 3130.0 * eps0;
        let eps33 = 3400.0 * eps0;
        #[rustfmt::skip]
        let epsilon_matrix: [[f64; 3]; 3] = [
            [eps11, 0.0, 0.0],
            [0.0, eps11, 0.0],
            [0.0, 0.0, eps33],
        ];
        PiezoMaterial {
            c_matrix,
            e_matrix,
            epsilon_matrix,
            name: "PZT-5H".to_string(),
        }
    }
    /// Computes mechanical stress `sigma = c * epsilon - e^T * E`.
    pub fn mechanical_stress(&self, strain: &[f64; 6], e_field: &[f64; 3]) -> [f64; 6] {
        let c_eps = mat6_vec6_mul(&self.c_matrix, strain);
        let mut et_e = [0.0_f64; 6];
        for (j, ete_j) in et_e.iter_mut().enumerate() {
            for (i, &ei) in e_field.iter().enumerate() {
                *ete_j += self.e_matrix[j][i] * ei;
            }
        }
        let mut sigma = [0.0_f64; 6];
        for (sig_k, (&ceps_k, &ete_k)) in sigma.iter_mut().zip(c_eps.iter().zip(et_e.iter())) {
            *sig_k = ceps_k - ete_k;
        }
        sigma
    }
    /// Computes electric displacement `D = e * epsilon + epsilon_s * E`.
    pub fn electric_displacement(&self, strain: &[f64; 6], e_field: &[f64; 3]) -> [f64; 3] {
        let e_eps = e_times_strain(&self.e_matrix, strain);
        let eps_e = mat3_vec3_mul(&self.epsilon_matrix, e_field);
        let mut d = [0.0_f64; 3];
        for (di, (&ee_i, &eps_i)) in d.iter_mut().zip(e_eps.iter().zip(eps_e.iter())) {
            *di = ee_i + eps_i;
        }
        d
    }
    /// Compute the piezoelectric coupling stiffness matrix K_ue (mechanical-electrical coupling).
    ///
    /// For a single element with uniform strain, the coupling matrix relates
    /// mechanical DOFs to electrical DOFs: K_ue = integral B^T e^T dV.
    ///
    /// Simplified: returns e^T * volume as a 6x3 matrix.
    pub fn coupling_matrix(&self, volume: f64) -> [[f64; 3]; 6] {
        let mut k_ue = [[0.0; 3]; 6];
        for (j, kue_j) in k_ue.iter_mut().enumerate() {
            for (i, kue_ji) in kue_j.iter_mut().enumerate() {
                *kue_ji = self.e_matrix[j][i] * volume;
            }
        }
        k_ue
    }
    /// Compute the dielectric stiffness matrix K_ee (electrical DOFs only).
    ///
    /// K_ee = integral N_phi^T epsilon_s N_phi dV.
    ///
    /// Simplified: returns epsilon_s * volume as a 3x3 matrix.
    pub fn dielectric_stiffness(&self, volume: f64) -> [[f64; 3]; 3] {
        let mut k_ee = [[0.0; 3]; 3];
        for (i, kee_i) in k_ee.iter_mut().enumerate() {
            for (j, kee_ij) in kee_i.iter_mut().enumerate() {
                *kee_ij = self.epsilon_matrix[i][j] * volume;
            }
        }
        k_ee
    }
}
/// Coupled mechanical-electrical system with voltage DOFs.
///
/// Stores the assembled system:
/// ```text
/// [K_uu   K_ue ] [u  ]   [F_mech ]
/// [K_ue^T K_ee ] [phi] = [Q_elec ]
/// ```
pub struct CoupledPiezoSystem {
    /// Number of mechanical DOFs.
    pub n_mech_dofs: usize,
    /// Number of electrical (voltage) DOFs.
    pub n_elec_dofs: usize,
    /// Mechanical stiffness K_uu (n_mech x n_mech), stored row-major.
    pub k_uu: Vec<f64>,
    /// Coupling matrix K_ue (n_mech x n_elec), stored row-major.
    pub k_ue: Vec<f64>,
    /// Dielectric stiffness K_ee (n_elec x n_elec), stored row-major.
    pub k_ee: Vec<f64>,
}
impl CoupledPiezoSystem {
    /// Create a new coupled system with zero matrices.
    pub fn new(n_mech: usize, n_elec: usize) -> Self {
        Self {
            n_mech_dofs: n_mech,
            n_elec_dofs: n_elec,
            k_uu: vec![0.0; n_mech * n_mech],
            k_ue: vec![0.0; n_mech * n_elec],
            k_ee: vec![0.0; n_elec * n_elec],
        }
    }
    /// Add a mechanical stiffness entry.
    pub fn add_k_uu(&mut self, i: usize, j: usize, val: f64) {
        self.k_uu[i * self.n_mech_dofs + j] += val;
    }
    /// Add a coupling entry.
    pub fn add_k_ue(&mut self, i: usize, j: usize, val: f64) {
        self.k_ue[i * self.n_elec_dofs + j] += val;
    }
    /// Add a dielectric stiffness entry.
    pub fn add_k_ee(&mut self, i: usize, j: usize, val: f64) {
        self.k_ee[i * self.n_elec_dofs + j] += val;
    }
    /// Compute the total size of the coupled system.
    pub fn total_dofs(&self) -> usize {
        self.n_mech_dofs + self.n_elec_dofs
    }
    /// Compute the mechanical response for prescribed voltages (static condensation).
    ///
    /// u = K_uu^{-1} (F - K_ue * phi)
    /// Simplified: for a single DOF each, returns u = (F - K_ue * phi) / K_uu
    pub fn static_condensation_1dof(&self, f_mech: f64, phi: f64) -> f64 {
        if self.n_mech_dofs != 1 || self.n_elec_dofs != 1 {
            return 0.0;
        }
        let k_uu = self.k_uu[0];
        let k_ue = self.k_ue[0];
        if k_uu.abs() < 1e-60 {
            return 0.0;
        }
        (f_mech - k_ue * phi) / k_uu
    }
}
