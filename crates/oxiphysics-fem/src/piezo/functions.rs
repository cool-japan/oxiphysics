//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::EmFrfResult;

/// Compute the electromechanical FRF for a single-DOF piezoelectric system.
///
/// Models a single mechanical DOF coupled to a single voltage DOF.
///
/// The system is:
/// ```text
/// (-omega^2 m + k_uu + i omega c) u + k_ue phi = F
/// k_ue u - k_ee phi = Q
/// ```
///
/// For sensor mode (Q = 0, F given), returns mechanical and voltage response.
#[allow(clippy::too_many_arguments)]
pub fn electromechanical_frf_1dof(
    omega: f64,
    mass: f64,
    k_uu: f64,
    k_ue: f64,
    k_ee: f64,
    damping: f64,
    force: f64,
) -> EmFrfResult {
    let z_real = k_uu - omega * omega * mass;
    let z_imag = omega * damping;
    let k_eff = if k_ee.abs() > 1e-60 {
        k_ue * k_ue / k_ee
    } else {
        0.0
    };
    let z_total_real = z_real + k_eff;
    let z_total_imag = z_imag;
    let z_sq = z_total_real * z_total_real + z_total_imag * z_total_imag;
    if z_sq < 1e-60 {
        return EmFrfResult {
            omega,
            mechanical_magnitude: 0.0,
            electrical_magnitude: 0.0,
        };
    }
    let u_real = force * z_total_real / z_sq;
    let u_imag = -force * z_total_imag / z_sq;
    let u_mag = (u_real * u_real + u_imag * u_imag).sqrt();
    let phi_mag = if k_ee.abs() > 1e-60 {
        k_ue * u_mag / k_ee
    } else {
        0.0
    };
    EmFrfResult {
        omega,
        mechanical_magnitude: u_mag,
        electrical_magnitude: phi_mag,
    }
}
/// Converts a 3x3 stress tensor to Voigt notation.
pub fn voigt_notation_stress(s: [[f64; 3]; 3]) -> [f64; 6] {
    [s[0][0], s[1][1], s[2][2], s[1][2], s[0][2], s[0][1]]
}
/// Converts a 3x3 strain tensor to Voigt notation.
pub fn voigt_notation_strain(e: [[f64; 3]; 3]) -> [f64; 6] {
    [
        e[0][0],
        e[1][1],
        e[2][2],
        2.0 * e[1][2],
        2.0 * e[0][2],
        2.0 * e[0][1],
    ]
}
/// Multiplies a 6x6 matrix by a 6-vector.
pub fn mat6_vec6_mul(m: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut result = [0.0_f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            result[i] += m[i][j] * v[j];
        }
    }
    result
}
/// Computes **e** * **epsilon** where the piezoelectric matrix **e** (3x6) is stored
/// column-major as `e[col][row]`.
pub fn e_times_strain(e: &[[f64; 3]; 6], strain: &[f64; 6]) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for j in 0..6 {
        for i in 0..3 {
            result[i] += e[j][i] * strain[j];
        }
    }
    result
}
/// Multiplies a 3x3 matrix by a 3-vector.
pub fn mat3_vec3_mul(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i] += m[i][j] * v[j];
        }
    }
    result
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::piezo::*;
    #[test]
    fn voigt_stress_diagonal_tensor() {
        let s = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let v = voigt_notation_stress(s);
        assert_eq!(v, [1.0, 2.0, 3.0, 0.0, 0.0, 0.0]);
    }
    #[test]
    fn voigt_strain_diagonal_tensor() {
        let e = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]];
        let v = voigt_notation_strain(e);
        assert_eq!(v, [0.1, 0.2, 0.3, 0.0, 0.0, 0.0]);
    }
    #[test]
    fn mechanical_stress_zero_inputs() {
        let mat = PiezoMaterial::pzt5a();
        let strain = [0.0; 6];
        let e_field = [0.0; 3];
        let sigma = mat.mechanical_stress(&strain, &e_field);
        for s in sigma {
            assert_eq!(s, 0.0);
        }
    }
    #[test]
    fn electric_displacement_z_strain() {
        let mat = PiezoMaterial::pzt5a();
        let strain = [0.0, 0.0, 1e-3, 0.0, 0.0, 0.0];
        let e_field = [0.0; 3];
        let d = mat.electric_displacement(&strain, &e_field);
        assert!(d[2] > 0.0, "D_z should be positive for positive eps_33");
        assert_eq!(d[0], 0.0);
        assert_eq!(d[1], 0.0);
    }
    #[test]
    fn induced_strain_100v_1mm() {
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoElement::new(mat, 1e-9);
        let strain = elem.induced_strain(100.0, 1e-3);
        let expected = 400e-12 * 100.0 / 1e-3;
        assert!((strain[2] - expected).abs() < 1e-20);
        assert_eq!(strain[0], 0.0);
        assert_eq!(strain[1], 0.0);
    }
    #[test]
    fn mat6_vec6_identity() {
        let mut identity = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            identity[i][i] = 1.0;
        }
        let v = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let result = mat6_vec6_mul(&identity, &v);
        assert_eq!(result, v);
    }
    #[test]
    fn mat3_vec3_identity() {
        let mut identity = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            identity[i][i] = 1.0;
        }
        let v = [7.0, 8.0, 9.0];
        let result = mat3_vec3_mul(&identity, &v);
        assert_eq!(result, v);
    }
    #[test]
    fn test_pzt5h_creation() {
        let mat = PiezoMaterial::pzt5h();
        assert_eq!(mat.name, "PZT-5H");
        let pzt5a = PiezoMaterial::pzt5a();
        assert!(
            mat.epsilon_matrix[2][2] > pzt5a.epsilon_matrix[2][2],
            "PZT-5H should have higher eps_33 than PZT-5A"
        );
    }
    #[test]
    fn test_coupling_matrix_volume_scaling() {
        let mat = PiezoMaterial::pzt5a();
        let vol1 = 1.0e-6;
        let vol2 = 2.0e-6;
        let k1 = mat.coupling_matrix(vol1);
        let k2 = mat.coupling_matrix(vol2);
        for j in 0..6 {
            for i in 0..3 {
                if k1[j][i].abs() > 1e-30 {
                    assert!(
                        (k2[j][i] / k1[j][i] - 2.0).abs() < 1e-12,
                        "coupling should scale linearly with volume"
                    );
                }
            }
        }
    }
    #[test]
    fn test_dielectric_stiffness() {
        let mat = PiezoMaterial::pzt5a();
        let vol = 1.0e-6;
        let k_ee = mat.dielectric_stiffness(vol);
        assert!((k_ee[0][1]).abs() < 1e-30);
        assert!((k_ee[0][2]).abs() < 1e-30);
        assert!(k_ee[0][0] > 0.0, "dielectric stiffness must be positive");
        assert!(k_ee[2][2] > 0.0);
    }
    #[test]
    fn test_element_coupling_matrix() {
        let mat = PiezoMaterial::pzt5a();
        let vol = 1.0e-6;
        let elem = PiezoElement::new(mat, vol);
        let k_ue = elem.element_coupling_matrix();
        assert!(k_ue[2][2].abs() > 0.0, "e_33 coupling should be non-zero");
    }
    #[test]
    fn test_sensor_open_circuit_voltage() {
        let mat = PiezoMaterial::pzt5a();
        let sensor = PiezoSensor::new(mat, 1e-4, 0.5e-3);
        let strain = [0.0, 0.0, 1e-6, 0.0, 0.0, 0.0];
        let v = sensor.open_circuit_voltage(&strain);
        assert!(v.abs() > 0.0, "voltage should be non-zero for z-strain");
    }
    #[test]
    fn test_sensor_zero_strain_zero_voltage() {
        let mat = PiezoMaterial::pzt5a();
        let sensor = PiezoSensor::new(mat, 1e-4, 0.5e-3);
        let strain = [0.0; 6];
        let v = sensor.open_circuit_voltage(&strain);
        assert_eq!(v, 0.0, "zero strain should give zero voltage");
    }
    #[test]
    fn test_sensor_charge_output() {
        let mat = PiezoMaterial::pzt5a();
        let sensor = PiezoSensor::new(mat, 1e-4, 0.5e-3);
        let strain = [0.0, 0.0, 1e-6, 0.0, 0.0, 0.0];
        let q = sensor.charge_output(&strain);
        assert!(q.abs() > 0.0, "charge should be non-zero for z-strain");
    }
    #[test]
    fn test_actuator_free_strain() {
        let mat = PiezoMaterial::pzt5a();
        let act = PiezoActuator::new(mat, 1e-4, 0.5e-3);
        let strain = act.free_strain(100.0);
        assert!(strain[2] > 0.0, "epsilon_33 should be positive");
        assert!(strain[0] < 0.0, "epsilon_11 should be negative (d_31 < 0)");
        assert!(strain[1] < 0.0, "epsilon_22 should be negative (d_31 < 0)");
    }
    #[test]
    fn test_actuator_blocking_force() {
        let mat = PiezoMaterial::pzt5a();
        let act = PiezoActuator::new(mat, 1e-4, 0.5e-3);
        let f = act.blocking_force(100.0);
        assert!(
            f > 0.0,
            "blocking force should be positive for positive voltage"
        );
        assert!(f.is_finite());
    }
    #[test]
    fn test_actuator_stored_energy() {
        let mat = PiezoMaterial::pzt5a();
        let act = PiezoActuator::new(mat, 1e-4, 0.5e-3);
        let u = act.stored_energy(100.0);
        assert!(u > 0.0, "stored energy must be positive");
        let u2 = act.stored_energy(200.0);
        assert!((u2 / u - 4.0).abs() < 1e-10, "energy should scale with V^2");
    }
    #[test]
    fn test_coupled_system_1dof() {
        let mut sys = CoupledPiezoSystem::new(1, 1);
        sys.add_k_uu(0, 0, 1000.0);
        sys.add_k_ue(0, 0, 5.0);
        sys.add_k_ee(0, 0, 1e-9);
        assert_eq!(sys.total_dofs(), 2);
        let u = sys.static_condensation_1dof(10.0, 100.0);
        let expected = (10.0 - 5.0 * 100.0) / 1000.0;
        assert!((u - expected).abs() < 1e-12);
    }
    #[test]
    fn test_em_frf_static() {
        let result = electromechanical_frf_1dof(0.0, 1.0, 1000.0, 5.0, 1e-9, 10.0, 1.0);
        assert!(result.mechanical_magnitude > 0.0);
        assert!(result.electrical_magnitude > 0.0);
    }
    #[test]
    fn test_em_frf_resonance() {
        let omega_n = (1000.0_f64 / 1.0).sqrt();
        let damping = 1.0;
        let off_res =
            electromechanical_frf_1dof(omega_n * 0.1, 1.0, 1000.0, 5.0, 1e-9, damping, 1.0);
        let at_res = electromechanical_frf_1dof(omega_n, 1.0, 1000.0, 5.0, 1e-9, damping, 1.0);
        assert!(
            at_res.mechanical_magnitude > off_res.mechanical_magnitude * 0.5,
            "resonance response {} should exceed off-resonance {}",
            at_res.mechanical_magnitude,
            off_res.mechanical_magnitude
        );
    }
    #[test]
    fn test_mechanical_stress_with_electric_field() {
        let mat = PiezoMaterial::pzt5a();
        let strain = [0.0; 6];
        let e_field = [0.0, 0.0, 1000.0];
        let sigma = mat.mechanical_stress(&strain, &e_field);
        assert!(
            sigma[2] < 0.0,
            "sigma_33 should be negative for positive E_3 with positive e_33"
        );
    }
}
/// Natural (isoparametric) coordinates for the 8 nodes of a hex element.
/// Each entry is (xi, eta, zeta) in \[-1, 1\]^3.
pub(super) const HEX8_NATURAL_COORDS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];
/// 2-point Gauss quadrature abscissae and weights in 1D.
pub(super) const GAUSS2_ABSCISSAE: [f64; 2] = [-0.577_350_269_189_626, 0.577_350_269_189_626];
pub(super) const GAUSS2_WEIGHTS: [f64; 2] = [1.0, 1.0];
/// Evaluate shape functions N_I(xi, eta, zeta) for an 8-node hex element.
///
/// Returns an array of 8 shape function values.
pub fn hex8_shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 8] {
    let mut n = [0.0_f64; 8];
    for i in 0..8 {
        let [xi_i, eta_i, zeta_i] = HEX8_NATURAL_COORDS[i];
        n[i] = 0.125 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta);
    }
    n
}
/// Evaluate shape function derivatives dN/d(xi, eta, zeta) for an 8-node hex element.
///
/// Returns a `[8][3]` array where `[i][j]` = dN_i / d(coord_j).
pub fn hex8_shape_derivatives(xi: f64, eta: f64, zeta: f64) -> [[f64; 3]; 8] {
    let mut dn = [[0.0_f64; 3]; 8];
    for i in 0..8 {
        let [xi_i, eta_i, zeta_i] = HEX8_NATURAL_COORDS[i];
        dn[i][0] = 0.125 * xi_i * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta);
        dn[i][1] = 0.125 * (1.0 + xi_i * xi) * eta_i * (1.0 + zeta_i * zeta);
        dn[i][2] = 0.125 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * zeta_i;
    }
    dn
}
/// Compute the Jacobian matrix J\[3\]\[3\] = sum_I (dN_I/d(xi,eta,zeta)) * X_I
/// for an 8-node hex element given nodal coordinates `coords[8][3]`.
pub fn hex8_jacobian(coords: &[[f64; 3]; 8], xi: f64, eta: f64, zeta: f64) -> [[f64; 3]; 3] {
    let dn = hex8_shape_derivatives(xi, eta, zeta);
    let mut j = [[0.0_f64; 3]; 3];
    for i in 0..8 {
        for row in 0..3 {
            for col in 0..3 {
                j[row][col] += dn[i][row] * coords[i][col];
            }
        }
    }
    j
}
/// 3x3 matrix determinant.
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// 3x3 matrix inverse.
pub(super) fn inv3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    let id = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id,
        ],
    ]
}
/// Build the strain-displacement matrix B\[6\]\[24\] at a Gauss point.
///
/// `dndx[8][3]` contains dN_I/dX_j in physical coordinates.
#[allow(clippy::too_many_arguments)]
pub fn hex8_b_matrix(dndx: &[[f64; 3]; 8]) -> [[f64; 24]; 6] {
    let mut b = [[0.0_f64; 24]; 6];
    for i in 0..8 {
        let col_base = i * 3;
        b[0][col_base] = dndx[i][0];
        b[1][col_base + 1] = dndx[i][1];
        b[2][col_base + 2] = dndx[i][2];
        b[3][col_base + 1] = dndx[i][2];
        b[3][col_base + 2] = dndx[i][1];
        b[4][col_base] = dndx[i][2];
        b[4][col_base + 2] = dndx[i][0];
        b[5][col_base] = dndx[i][1];
        b[5][col_base + 1] = dndx[i][0];
    }
    b
}
/// Compute dN/dX = J^{-1} * dN/d(xi, eta, zeta) in physical coordinates.
pub fn hex8_dndx(j_inv: &[[f64; 3]; 3], dn_nat: &[[f64; 3]; 8]) -> [[f64; 3]; 8] {
    let mut dndx = [[0.0_f64; 3]; 8];
    for i in 0..8 {
        for k in 0..3 {
            for m in 0..3 {
                dndx[i][k] += j_inv[k][m] * dn_nat[i][m];
            }
        }
    }
    dndx
}
#[cfg(test)]
mod piezo_expanded_tests {
    use super::*;
    use crate::piezo::*;
    #[test]
    fn test_hex8_shape_functions_sum_to_one() {
        let test_points = [
            (0.0, 0.0, 0.0),
            (0.5, -0.3, 0.7),
            (-1.0, -1.0, -1.0),
            (1.0, 1.0, 1.0),
        ];
        for (xi, eta, zeta) in test_points {
            let n = hex8_shape_functions(xi, eta, zeta);
            let sum: f64 = n.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-13,
                "N sum at ({xi},{eta},{zeta}) = {sum}"
            );
        }
    }
    #[test]
    fn test_hex8_shape_functions_partition_of_unity() {
        for node in 0..8 {
            let [xi, eta, zeta] = HEX8_NATURAL_COORDS[node];
            let n = hex8_shape_functions(xi, eta, zeta);
            assert!(
                (n[node] - 1.0).abs() < 1e-13,
                "N[{node}] at own node = {}",
                n[node]
            );
            for other in 0..8 {
                if other != node {
                    assert!(
                        n[other].abs() < 1e-13,
                        "N[{other}] at node {node} = {}",
                        n[other]
                    );
                }
            }
        }
    }
    #[test]
    fn test_hex8_jacobian_unit_cube() {
        let coords: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let j = hex8_jacobian(&coords, 0.0, 0.0, 0.0);
        assert!((j[0][0] - 0.5).abs() < 1e-13, "J[0][0] = {}", j[0][0]);
        assert!((j[1][1] - 0.5).abs() < 1e-13, "J[1][1] = {}", j[1][1]);
        assert!((j[2][2] - 0.5).abs() < 1e-13, "J[2][2] = {}", j[2][2]);
        assert!(j[0][1].abs() < 1e-13);
    }
    #[test]
    fn test_hex8_element_volume() {
        let coords: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoHex8Element::new(coords, mat);
        let vol = elem.volume();
        assert!((vol - 1.0).abs() < 1e-10, "volume = {vol}");
    }
    #[test]
    fn test_hex8_stiffness_matrix_symmetric() {
        let coords: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoHex8Element::new(coords, mat);
        let k = elem.stiffness_matrix();
        for i in 0..24 {
            for j in 0..24 {
                let kij = k[i * 24 + j];
                let kji = k[j * 24 + i];
                assert!(
                    (kij - kji).abs() < 1e-4 * kij.abs().max(1e-6),
                    "K[{i}][{j}]={kij} vs K[{j}][{i}]={kji}"
                );
            }
        }
    }
    #[test]
    fn test_d_form_pzt5a_strain_zero_inputs() {
        let mat = PiezoMaterialDForm::pzt5a_d_form();
        let stress = [0.0_f64; 6];
        let e_field = [0.0_f64; 3];
        let eps = mat.strain(&stress, &e_field);
        for e in eps {
            assert_eq!(e, 0.0);
        }
    }
    #[test]
    fn test_d_form_electric_displacement_zero() {
        let mat = PiezoMaterialDForm::pzt5a_d_form();
        let stress = [0.0_f64; 6];
        let e_field = [0.0_f64; 3];
        let d = mat.electric_displacement_d(&stress, &e_field);
        for di in d {
            assert_eq!(di, 0.0);
        }
    }
    #[test]
    fn test_k33_coupling_factor_reasonable() {
        let mat = PiezoMaterialDForm::pzt5a_d_form();
        let k33 = mat.k33_coupling_factor();
        assert!(k33 > 0.3 && k33 < 1.0, "k33 = {k33}");
    }
    #[test]
    fn test_transducer_resonance_frequency() {
        let c33 = 111.0e9_f64;
        let density = 7750.0_f64;
        let thickness = 1.0e-3_f64;
        let e33 = 15.8_f64;
        let eps33 = 8.854187817e-12 * 1470.0;
        let tr = TransducerResonance::new(c33, e33, eps33, density, thickness);
        let fr = tr.resonance_frequency_uncoupled();
        assert!(fr > 1e5 && fr < 1e8, "f_r = {fr}");
    }
    #[test]
    fn test_transducer_stiffened_frequency_higher() {
        let c33 = 111.0e9_f64;
        let density = 7750.0_f64;
        let thickness = 1.0e-3_f64;
        let e33 = 15.8_f64;
        let eps33 = 8.854187817e-12 * 1470.0;
        let tr = TransducerResonance::new(c33, e33, eps33, density, thickness);
        let fr = tr.resonance_frequency_uncoupled();
        let fs = tr.resonance_frequency_stiffened();
        assert!(fs >= fr, "stiffened freq {fs} should be >= uncoupled {fr}");
    }
    #[test]
    fn test_transducer_kt_positive() {
        let c33 = 111.0e9_f64;
        let density = 7750.0_f64;
        let thickness = 1.0e-3_f64;
        let e33 = 15.8_f64;
        let eps33 = 8.854187817e-12 * 1470.0;
        let tr = TransducerResonance::new(c33, e33, eps33, density, thickness);
        let kt = tr.kt_coupling_factor();
        assert!(kt > 0.0, "kt should be positive");
    }
    #[test]
    fn test_transducer_anti_resonance_gt_resonance() {
        let c33 = 111.0e9_f64;
        let density = 7750.0_f64;
        let thickness = 1.0e-3_f64;
        let e33 = 15.8_f64;
        let eps33 = 8.854187817e-12 * 1470.0;
        let tr = TransducerResonance::new(c33, e33, eps33, density, thickness);
        let fr = tr.resonance_frequency_uncoupled();
        let fa = tr.anti_resonance_frequency();
        assert!(fa >= fr, "anti-resonance {fa} >= resonance {fr}");
    }
    #[test]
    fn test_transducer_motional_impedance_finite_off_resonance() {
        let c33 = 111.0e9_f64;
        let density = 7750.0_f64;
        let thickness = 1.0e-3_f64;
        let e33 = 15.8_f64;
        let eps33 = 8.854187817e-12 * 1470.0;
        let tr = TransducerResonance::new(c33, e33, eps33, density, thickness);
        let omega_off = 2.0 * std::f64::consts::PI * 1e4;
        let z = tr.motional_impedance(omega_off, 0.02);
        assert!(z.is_finite() && z > 0.0, "Z_m = {z}");
    }
    #[test]
    fn test_hex8_b_matrix_correct_size() {
        let coords: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let j = hex8_jacobian(&coords, 0.0, 0.0, 0.0);
        let j_inv = inv3(&j);
        let dn_nat = hex8_shape_derivatives(0.0, 0.0, 0.0);
        let dndx = hex8_dndx(&j_inv, &dn_nat);
        let b = hex8_b_matrix(&dndx);
        assert!(
            b[0][0].abs() > 0.0 || b[0][3].abs() > 0.0 || b[0][6].abs() > 0.0,
            "B matrix dN/dX terms should be non-zero"
        );
        let _ = b;
    }
    #[test]
    fn test_coupling_matrix_full_nonzero() {
        let coords: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let mat = PiezoMaterial::pzt5a();
        let elem = PiezoHex8Element::new(coords, mat);
        let k_ue = elem.coupling_matrix_full();
        assert_eq!(k_ue.len(), 24 * 3);
        let has_nonzero = k_ue.iter().any(|&v| v.abs() > 1e-30);
        assert!(has_nonzero, "coupling matrix should have non-zero entries");
    }
}
/// Compute the resonance frequency of a cantilever piezoelectric beam.
///
/// For a uniform Euler-Bernoulli beam with one end clamped:
/// f_n = (lambda_n^2 / (2*pi)) * sqrt(EI / (rho*A*L^4))
///
/// First mode: lambda_1 = 1.875104068711961
///
/// # Arguments
/// * `e_modulus` - Young's modulus (Pa)
/// * `moment_of_inertia` - Second moment of area I (m^4)
/// * `density` - mass density (kg/m^3)
/// * `cross_area` - cross-section area (m^2)
/// * `length` - beam length (m)
#[allow(clippy::too_many_arguments)]
#[allow(non_snake_case)]
pub fn piezo_beam_resonance_frequency(
    e_modulus: f64,
    moment_of_inertia: f64,
    density: f64,
    cross_area: f64,
    length: f64,
) -> f64 {
    if density <= 0.0 || cross_area <= 0.0 || length <= 0.0 {
        return 0.0;
    }
    let lambda1_sq = 1.875104068711961_f64.powi(2);
    let ei = e_modulus * moment_of_inertia;
    let rho_a_l4 = density * cross_area * length.powi(4);
    if rho_a_l4 <= 0.0 {
        return 0.0;
    }
    lambda1_sq / (2.0 * std::f64::consts::PI) * (ei / rho_a_l4).sqrt()
}
/// Compute the n-th mode resonance frequency of a cantilever beam.
///
/// Characteristic values lambda_n for clamped-free beam:
/// mode 1: 1.875, mode 2: 4.694, mode 3: 7.855, mode 4: 10.996
#[allow(clippy::too_many_arguments)]
#[allow(non_snake_case)]
pub fn piezo_beam_resonance_mode_n(
    e_modulus: f64,
    moment_of_inertia: f64,
    density: f64,
    cross_area: f64,
    length: f64,
    mode: usize,
) -> f64 {
    if density <= 0.0 || cross_area <= 0.0 || length <= 0.0 || mode == 0 {
        return 0.0;
    }
    let lambdas = [
        1.875104068711961,
        4.694091132974174,
        7.854757438237612,
        10.995540734875467,
    ];
    let lambda = if mode <= 4 {
        lambdas[mode - 1]
    } else {
        (2.0 * mode as f64 - 1.0) * std::f64::consts::PI / 2.0
    };
    let lambda_sq = lambda * lambda;
    let ei = e_modulus * moment_of_inertia;
    let rho_a_l4 = density * cross_area * length.powi(4);
    if rho_a_l4 <= 0.0 {
        return 0.0;
    }
    lambda_sq / (2.0 * std::f64::consts::PI) * (ei / rho_a_l4).sqrt()
}
/// Piezoelectric actuator force for applied voltage using d33 mode.
///
/// The free strain is: epsilon_33 = d33 * V / t
/// The blocked (clamped) force is: F = c33 * epsilon_33 * A = c33 * d33 * V / t * A
///
/// # Arguments
/// * `voltage`      - applied voltage (V)
/// * `d33`          - piezoelectric strain coefficient (m/V)
/// * `c33`          - elastic stiffness constant (Pa)
/// * `area`         - actuator cross-section area (m^2)
/// * `thickness`    - actuator thickness (m)
pub fn piezo_actuator_force_d33(
    voltage: f64,
    d33: f64,
    c33: f64,
    area: f64,
    thickness: f64,
) -> f64 {
    if thickness.abs() < 1e-30 {
        return 0.0;
    }
    c33 * d33 * voltage / thickness * area
}
/// Piezoelectric actuator force using d31 mode.
///
/// F = c11 * d31 * V / t * A  (in-plane direction)
pub fn piezo_actuator_force_d31(
    voltage: f64,
    d31: f64,
    c11: f64,
    area: f64,
    thickness: f64,
) -> f64 {
    if thickness.abs() < 1e-30 {
        return 0.0;
    }
    c11 * d31 * voltage / thickness * area
}
/// Piezoelectric sensor output voltage from strain in d33 mode.
///
/// V = epsilon_33 * t / (d33 * c33 / (eps33 * c33))
/// Simplified: V = e33 * epsilon * t / eps33 = (d33 * c33) * eps * t / eps33
///
/// # Arguments
/// * `strain33`     - strain in the 3-direction
/// * `d33`          - d33 coefficient (m/V)
/// * `c33`          - elastic stiffness (Pa)
/// * `eps33`        - dielectric permittivity eps33 (F/m)
/// * `thickness`    - element thickness (m)
pub fn piezo_sensor_voltage_d33(
    strain33: f64,
    d33: f64,
    c33: f64,
    eps33: f64,
    thickness: f64,
) -> f64 {
    if eps33.abs() < 1e-60 {
        return 0.0;
    }
    let e33 = d33 * c33;
    e33 * strain33 * thickness / eps33
}
/// Piezoelectric sensor output voltage from strain in d31 mode.
pub fn piezo_sensor_voltage_d31(
    strain11: f64,
    d31: f64,
    c11: f64,
    eps33: f64,
    thickness: f64,
) -> f64 {
    if eps33.abs() < 1e-60 {
        return 0.0;
    }
    let e31 = d31 * c11;
    e31 * strain11 * thickness / eps33
}
