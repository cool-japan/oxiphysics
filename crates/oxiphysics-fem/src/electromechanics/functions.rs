//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
/// Compute the electric field **E = −∇φ** inside a linear tetrahedral element.
///
/// `phi` – the four nodal potential values `[φ1, φ2, φ3, φ4]`.
/// `nodes` – the four node coordinates `[[x,y,z\]; 4]`.
/// `_volume` – element volume (not needed for constant-field tet, kept for API
///             consistency).
pub fn compute_electric_field_from_potential(
    phi: &[f64],
    nodes: &[[f64; 3]; 4],
    _volume: f64,
) -> [f64; 3] {
    let [p0, p1, p2, p3] = nodes;
    let j = [
        [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
        [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
        [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
    ];
    let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
        - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
        + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
    let inv_det = 1.0 / det;
    let ji = [
        [
            (j[1][1] * j[2][2] - j[1][2] * j[2][1]) * inv_det,
            (j[0][2] * j[2][1] - j[0][1] * j[2][2]) * inv_det,
            (j[0][1] * j[1][2] - j[0][2] * j[1][1]) * inv_det,
        ],
        [
            (j[1][2] * j[2][0] - j[1][0] * j[2][2]) * inv_det,
            (j[0][0] * j[2][2] - j[0][2] * j[2][0]) * inv_det,
            (j[0][2] * j[1][0] - j[0][0] * j[1][2]) * inv_det,
        ],
        [
            (j[1][0] * j[2][1] - j[1][1] * j[2][0]) * inv_det,
            (j[0][1] * j[2][0] - j[0][0] * j[2][1]) * inv_det,
            (j[0][0] * j[1][1] - j[0][1] * j[1][0]) * inv_det,
        ],
    ];
    let grad_ref: [[f64; 3]; 4] = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut grad_phi = [0.0_f64; 3];
    for i in 0..4 {
        let mut dn = [0.0_f64; 3];
        for row in 0..3 {
            let mut s = 0.0;
            for col in 0..3 {
                s += ji[col][row] * grad_ref[i][col];
            }
            dn[row] = s;
        }
        for k in 0..3 {
            grad_phi[k] += dn[k] * phi[i];
        }
    }
    [-grad_phi[0], -grad_phi[1], -grad_phi[2]]
}
#[cfg(test)]
mod tests {

    use crate::electromechanics::*;
    /// The unit tetrahedron with vertices (0,0,0),(1,0,0),(0,1,0),(0,0,1) has
    /// volume 1/6.
    #[test]
    fn test_unit_tet_volume() {
        let nodes: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let vol = PiezoFemAssembler::compute_element_volume(&nodes);
        assert!((vol - 1.0 / 6.0).abs() < 1e-12, "volume = {vol}");
    }
    /// `blocked_force` must be positive when voltage > 0 and e33 > 0.
    #[test]
    fn test_blocked_force_positive() {
        let mat = PiezoMaterial::isotropic_with_piezo(67e9, 0.31, 1200.0, 15.1, -5.2, 12.3);
        let actuator = PiezoActuator {
            length: 0.01,
            area: 1e-4,
            material: mat,
            voltage: 100.0,
        };
        let f = actuator.blocked_force();
        assert!(f > 0.0, "blocked_force = {f}");
    }
    /// The Gauss–Seidel solver must converge in fewer than `max_iter`
    /// iterations on a small 10×10 grid with boundary voltages set.
    #[test]
    fn test_electrostatic_solver_converges() {
        let nx = 10;
        let ny = 10;
        let mut solver = ElectrostaticSolver::new(nx, ny, 1e-3);
        for ix in 0..nx {
            solver.phi[(ny - 1) * nx + ix] = 1.0;
        }
        let max_iter = 50_000;
        let iters = solver.solve_laplace_gauss_seidel(max_iter, 1e-8);
        assert!(
            iters < max_iter,
            "did not converge in {max_iter} iterations (ran {iters})"
        );
    }
    /// `isotropic_with_piezo` should produce a symmetric elastic stiffness.
    #[test]
    fn test_elastic_stiffness_symmetry() {
        let mat = PiezoMaterial::isotropic_with_piezo(200e9, 0.3, 1.0, 0.0, 0.0, 0.0);
        let c = mat.elastic_stiffness;
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c[i][j] - c[j][i]).abs() < 1e-6,
                    "C[{i}][{j}]={} != C[{j}][{i}]={}",
                    c[i][j],
                    c[j][i]
                );
            }
        }
    }
    /// `compute_b_matrix` must return a 6×12 matrix (6 strain components,
    /// 12 mechanical DOFs).
    #[test]
    fn test_b_matrix_dimensions() {
        let nodes: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let b = PiezoFemAssembler::compute_b_matrix(&nodes);
        assert_eq!(b.len(), 6);
        for row in &b {
            assert_eq!(row.len(), 12);
        }
    }
}
/// Compute the capacitance matrix between electrode groups.
///
/// The global stiffness matrix K of the Laplace problem gives capacitances
/// via C = K (boundary DOF sub-matrix).  This function extracts the
/// sub-matrix for the given electrode node sets.
///
/// # Arguments
/// * `k_global` – global n×n electrostatic stiffness (row-major flat Vec)
/// * `n` – number of DOFs
/// * `electrode_a` – indices of electrode A nodes
/// * `electrode_b` – indices of electrode B nodes
///
/// Returns the 2×2 capacitance matrix \[\[C_aa, C_ab\\], \[C_ba, C_bb\]].
pub fn capacitance_matrix(
    k_global: &[Vec<f64>],
    electrode_a: &[usize],
    electrode_b: &[usize],
) -> [[f64; 2]; 2] {
    let sum_kab = |rows: &[usize], cols: &[usize]| -> f64 {
        let mut s = 0.0;
        for &r in rows {
            for &c in cols {
                if r < k_global.len() && c < k_global[r].len() {
                    s += k_global[r][c];
                }
            }
        }
        s
    };
    let kaa = sum_kab(electrode_a, electrode_a);
    let kab = sum_kab(electrode_a, electrode_b);
    let kba = sum_kab(electrode_b, electrode_a);
    let kbb = sum_kab(electrode_b, electrode_b);
    [[kaa, kab], [kba, kbb]]
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::electromechanics::*;
    #[test]
    fn test_electrostatic_element_stiffness_positive_definite() {
        let elem = ElectrostaticElement {
            coords: [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            permittivity: 8.854e-12,
        };
        let k = elem.element_stiffness();
        for i in 0..4 {
            assert!(k[i][i] > 0.0, "k[{i}][{i}] = {}", k[i][i]);
        }
        for i in 0..4 {
            for j in 0..4 {
                assert!((k[i][j] - k[j][i]).abs() < 1e-25, "asymmetry at ({i},{j})");
            }
        }
    }
    #[test]
    fn test_electrostatic_element_electric_field() {
        let elem = ElectrostaticElement {
            coords: [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            permittivity: 8.854e-12,
        };
        let v = 100.0;
        let nodal_phi = [0.0, 0.0, 0.0, v];
        let e = elem.electric_field(&nodal_phi);
        assert!(e[2].abs() > 0.0, "Ez should be non-zero");
        assert!((e[2] - (-v)).abs() < 1e-8, "Ez = {}", e[2]);
    }
    #[test]
    fn test_capacitance_matrix_symmetry() {
        let k: Vec<Vec<f64>> = vec![vec![2.0, -1.0], vec![-1.0, 2.0]];
        let c = capacitance_matrix(&k, &[0], &[1]);
        assert!((c[0][1] - c[1][0]).abs() < 1e-12, "C should be symmetric");
    }
    #[test]
    fn test_dielectric_breakdown_below_threshold() {
        let bd = DielectricBreakdown::new(3e7, 3.5);
        assert!(!bd.is_broken_down(1e7), "should not break down at 10 MV/m");
        assert!(bd.is_broken_down(4e7), "should break down at 40 MV/m");
    }
    #[test]
    fn test_dielectric_breakdown_safety_margin() {
        let bd = DielectricBreakdown::new(3e7, 3.5);
        let margin = bd.safety_margin(1e7);
        assert!((margin - 3.0).abs() < 1e-10, "safety margin = {margin}");
    }
    #[test]
    fn test_dielectric_breakdown_energy_density() {
        let bd = DielectricBreakdown::new(3e7, 3.5);
        let u = bd.energy_density_at_breakdown();
        assert!(u > 0.0 && u.is_finite(), "energy density = {u}");
    }
    #[test]
    fn test_mems_electrostatic_force_increases_with_voltage() {
        let mems = MemsActuator::new(1e-6, 2e-6, 100.0, 1.0);
        let f1 = mems.electrostatic_force(10.0, 2e-6);
        let f2 = mems.electrostatic_force(20.0, 2e-6);
        assert!(f2 > f1, "force should increase with voltage");
    }
    #[test]
    fn test_mems_pull_in_voltage_positive() {
        let mems = MemsActuator::new(1e-6, 2e-6, 100.0, 1.0);
        let v_pi = mems.pull_in_voltage();
        assert!(v_pi > 0.0 && v_pi.is_finite(), "V_pi = {v_pi}");
    }
    #[test]
    fn test_mems_equilibrium_below_pull_in() {
        let mems = MemsActuator::new(1e-6, 2e-6, 100.0, 1.0);
        let v_pi = mems.pull_in_voltage();
        let gap = mems.equilibrium_gap(v_pi * 0.5);
        assert!(gap.is_some(), "should have equilibrium below pull-in");
        let g = gap.unwrap();
        assert!(g > 0.0 && g < mems.rest_gap, "gap = {g}");
    }
    #[test]
    fn test_mems_capacitance() {
        let mems = MemsActuator::new(1e-6, 2e-6, 100.0, 1.0);
        let c1 = mems.capacitance(2e-6);
        let c2 = mems.capacitance(1e-6);
        assert!(c2 > c1, "smaller gap should give larger capacitance");
    }
    #[test]
    fn test_piezo_coupling_d33_positive() {
        let pc = PiezoCoupling::new(67e9, 0.31, 15.1, -5.2, 12.3);
        let d33 = pc.d33();
        assert!(d33 > 0.0, "d33 should be positive: {d33}");
    }
    #[test]
    fn test_piezo_coupling_induced_strain() {
        let pc = PiezoCoupling::new(67e9, 0.31, 15.1, -5.2, 12.3);
        let eps = pc.induced_strain_33(1e6);
        assert!(eps > 0.0, "induced strain should be positive: {eps}");
    }
    #[test]
    fn test_piezo_coupling_k33_range() {
        pub(super) const EPS0: f64 = 8.854_187_817e-12;
        let pc = PiezoCoupling::new(67e9, 0.31, 15.1, -5.2, 12.3);
        let eps33 = 1200.0 * EPS0;
        let k33 = pc.k33(eps33);
        assert!(k33 > 0.0, "k33 should be positive: {k33}");
    }
    #[test]
    fn test_electric_field_from_potential_uniform_field() {
        let nodes: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.001],
        ];
        let vol = PiezoFemAssembler::compute_element_volume(&nodes);
        let phi = [0.0, 0.0, 0.0, -1.0];
        let e = compute_electric_field_from_potential(&phi, &nodes, vol);
        assert!(e[2].abs() > 0.0, "Ez should be nonzero");
    }
    #[test]
    fn test_electrostatic_solver_potential_gradient() {
        let nx = 5;
        let ny = 5;
        let mut solver = ElectrostaticSolver::new(nx, ny, 1e-3);
        for ix in 0..nx {
            solver.phi[(ny - 1) * nx + ix] = 1.0;
        }
        solver.solve_laplace_gauss_seidel(10000, 1e-6);
        for iy in 1..ny - 1 {
            for ix in 1..nx - 1 {
                let phi = solver.phi[iy * nx + ix];
                assert!(phi > -1e-6 && phi < 1.0 + 1e-6, "phi({ix},{iy}) = {phi}");
            }
        }
    }
}
/// Compute the Maxwell stress tensor for a dielectric with permittivity `eps`.
///
/// T_ij = ε (E_i E_j - δ_ij |E|²/2)
#[allow(dead_code)]
pub fn maxwell_stress_tensor(e_field: &[f64; 3], eps: f64) -> [[f64; 3]; 3] {
    let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            t[i][j] = eps * (e_field[i] * e_field[j] - 0.5 * delta * e_sq);
        }
    }
    t
}
/// Compute the Maxwell stress tensor trace.
///
/// tr(T_Maxwell) = ε (|E|² - 3/2 |E|²) = -ε |E|²/2
#[allow(dead_code)]
pub fn maxwell_stress_trace(e_field: &[f64; 3], eps: f64) -> f64 {
    let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
    -0.5 * eps * e_sq
}
/// Electrostatic energy density: U = ε |E|² / 2.
#[allow(dead_code)]
pub fn electrostatic_energy_density(e_field: &[f64; 3], eps: f64) -> f64 {
    let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
    0.5 * eps * e_sq
}
/// Electric force on a dielectric volume from Maxwell stress tensor.
///
/// For a uniform field, the net force from T_Maxwell on a surface element n dA.
/// F_i = Σ_j T_ij n_j * area
#[allow(dead_code)]
pub fn maxwell_force(e_field: &[f64; 3], normal: &[f64; 3], area: f64, eps: f64) -> [f64; 3] {
    let t = maxwell_stress_tensor(e_field, eps);
    let mut f = [0.0_f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            f[i] += t[i][j] * normal[j];
        }
        f[i] *= area;
    }
    f
}
#[cfg(test)]
mod tests_electromechanics_extended {
    use super::*;
    use crate::electromechanics::*;
    fn unit_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    /// d-form strain should be zero under zero stress and field.
    #[test]
    fn test_piezo_d_zero_input() {
        let pd = PiezoD::from_pzt(67e9, 0.31, 1200.0, 400e-12, -185e-12, 580e-12);
        let s = pd.strain(&[0.0; 6], &[0.0; 3]);
        for (i, &si) in s.iter().enumerate() {
            assert!(si.abs() < 1e-30, "strain[{i}] = {si} should be 0");
        }
    }
    /// d-form compliance matrix should be symmetric.
    #[test]
    fn test_piezo_d_compliance_symmetry() {
        let pd = PiezoD::from_pzt(67e9, 0.31, 1200.0, 400e-12, -185e-12, 580e-12);
        for i in 0..6 {
            for j in 0..6 {
                let diff = (pd.compliance[i][j] - pd.compliance[j][i]).abs();
                assert!(diff < 1e-30, "compliance[{i}][{j}] ≠ compliance[{j}][{i}]");
            }
        }
    }
    /// Electric displacement under zero inputs is zero.
    #[test]
    fn test_piezo_d_zero_displacement() {
        let pd = PiezoD::from_pzt(67e9, 0.31, 1200.0, 400e-12, -185e-12, 580e-12);
        let d = pd.electric_displacement(&[0.0; 6], &[0.0; 3]);
        for (i, &di) in d.iter().enumerate() {
            assert!(di.abs() < 1e-30, "D[{i}] = {di} should be 0");
        }
    }
    /// Electric displacement is non-zero under applied stress.
    #[test]
    fn test_piezo_d_displacement_nonzero() {
        let pd = PiezoD::from_pzt(67e9, 0.31, 1200.0, 400e-12, -185e-12, 580e-12);
        let stress = [1e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let d = pd.electric_displacement(&stress, &[0.0; 3]);
        let d_mag: f64 = d.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(
            d_mag > 0.0,
            "electric displacement should be nonzero under stress: {d_mag}"
        );
    }
    /// Electrostrictive strain is zero for zero field.
    #[test]
    fn test_electrostriction_zero_field() {
        let m = ElectrostrictiveMaterial::new(1e-18, 1e-19);
        let eps = m.strain_tensor(&[0.0; 3]);
        for i in 0..3 {
            for j in 0..3 {
                assert!(eps[i][j].abs() < 1e-30);
            }
        }
    }
    /// Electrostrictive strain tensor is symmetric.
    #[test]
    fn test_electrostriction_symmetry() {
        let m = ElectrostrictiveMaterial::new(1e-18, 1e-19);
        let e = [1e6, 0.5e6, 0.0];
        let eps = m.strain_tensor(&e);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (eps[i][j] - eps[j][i]).abs() < 1e-30,
                    "eps[{i}][{j}] ≠ eps[{j}][{i}]"
                );
            }
        }
    }
    /// Electrostrictive energy density is positive for nonzero field.
    #[test]
    fn test_electrostriction_energy_positive() {
        let m = ElectrostrictiveMaterial::new(1e-18, 1e-19);
        let e = [1e6, 0.0, 0.0];
        let u = m.energy_density(&e);
        assert!(u > 0.0, "energy density = {u}");
    }
    /// Maxwell stress is symmetric.
    #[test]
    fn test_dea_maxwell_stress_symmetry() {
        let dea = DielectricElastomer::new(1e6, 5e6, 4.0, 100e-6);
        let e = [1e6, 0.5e6, 0.0];
        let t = dea.maxwell_stress_tensor(&e);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (t[i][j] - t[j][i]).abs() < 1e-20,
                    "Maxwell stress not symmetric at ({i},{j})"
                );
            }
        }
    }
    /// Maxwell pressure increases with voltage.
    #[test]
    fn test_dea_maxwell_pressure_increases_with_voltage() {
        let dea = DielectricElastomer::new(1e6, 5e6, 4.0, 100e-6);
        let p1 = dea.maxwell_pressure(100.0);
        let p2 = dea.maxwell_pressure(200.0);
        assert!(p2 > p1, "Maxwell pressure should increase with voltage");
    }
    /// Actuation strain is negative (thinning).
    #[test]
    fn test_dea_actuation_strain_negative() {
        let dea = DielectricElastomer::new(1e6, 5e6, 4.0, 100e-6);
        let eps = dea.actuation_strain(1000.0);
        assert!(
            eps < 0.0,
            "actuation strain should be negative (thinning): {eps}"
        );
    }
    /// Electrostatic energy density is positive.
    #[test]
    fn test_dea_electrostatic_energy_positive() {
        let dea = DielectricElastomer::new(1e6, 5e6, 4.0, 100e-6);
        let u = dea.electrostatic_energy_density(500.0);
        assert!(u > 0.0, "electrostatic energy density = {u}");
    }
    /// Maxwell stress tensor is symmetric.
    #[test]
    fn test_maxwell_stress_tensor_symmetry() {
        let e = [1.0e5, 2.0e5, 0.5e5];
        let t = maxwell_stress_tensor(&e, 8.854e-12);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (t[i][j] - t[j][i]).abs() < 1e-30,
                    "not symmetric at ({i},{j})"
                );
            }
        }
    }
    /// Maxwell stress trace is negative for non-zero field.
    #[test]
    fn test_maxwell_stress_trace_negative() {
        let e = [0.0, 0.0, 1e5];
        let tr = maxwell_stress_trace(&e, 8.854e-12);
        assert!(tr < 0.0, "Maxwell trace should be negative: {tr}");
    }
    /// Electrostatic energy density is positive.
    #[test]
    fn test_electrostatic_energy_density_positive() {
        let e = [1e3, 2e3, 3e3];
        let u = electrostatic_energy_density(&e, 8.854e-12);
        assert!(u > 0.0, "energy density = {u}");
    }
    /// Maxwell force is non-zero for non-zero field and normal.
    #[test]
    fn test_maxwell_force_nonzero() {
        let e = [0.0, 0.0, 1e5];
        let normal = [0.0, 0.0, 1.0];
        let f = maxwell_force(&e, &normal, 1e-4, 8.854e-12);
        let f_mag: f64 = f.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(f_mag > 0.0, "force magnitude = {f_mag}");
    }
    /// Global assembler produces symmetric matrix.
    #[test]
    fn test_global_assembler_symmetry() {
        let mut local = PiezoFemAssembler::new();
        local.nodes = vec![
            ElectromechanicalNode::new(0.0, 0.0, 0.0),
            ElectromechanicalNode::new(1.0, 0.0, 0.0),
            ElectromechanicalNode::new(0.0, 1.0, 0.0),
            ElectromechanicalNode::new(0.0, 0.0, 1.0),
        ];
        let mat = PiezoMaterial::isotropic_with_piezo(67e9, 0.31, 1200.0, 15.1, -5.2, 12.3);
        local.materials = vec![mat];
        local.elements = vec![ElectromechanicalElement {
            nodes: [0, 1, 2, 3],
            material_id: 0,
        }];
        let mut global = GlobalPiezoAssembler::new(local);
        global.assemble();
        let n = global.n_dof;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (global.k_global[i][j] - global.k_global[j][i]).abs() < 1e-8,
                    "K_global[{i}][{j}] ≠ K_global[{j}][{i}]"
                );
            }
        }
    }
    /// Extracted K_uu sub-block has correct size.
    #[test]
    fn test_global_assembler_kuu_size() {
        let mut local = PiezoFemAssembler::new();
        local.nodes = vec![
            ElectromechanicalNode::new(0.0, 0.0, 0.0),
            ElectromechanicalNode::new(1.0, 0.0, 0.0),
            ElectromechanicalNode::new(0.0, 1.0, 0.0),
            ElectromechanicalNode::new(0.0, 0.0, 1.0),
        ];
        let mat = PiezoMaterial::isotropic_with_piezo(67e9, 0.31, 1200.0, 0.0, 0.0, 0.0);
        local.materials = vec![mat];
        local.elements = vec![ElectromechanicalElement {
            nodes: [0, 1, 2, 3],
            material_id: 0,
        }];
        let mut global = GlobalPiezoAssembler::new(local);
        global.assemble();
        let kuu = global.extract_k_uu();
        assert_eq!(kuu.len(), 12);
        assert_eq!(kuu[0].len(), 12);
    }
    /// Extracted K_pp sub-block has correct size.
    #[test]
    fn test_global_assembler_kpp_size() {
        let mut local = PiezoFemAssembler::new();
        local.nodes = vec![
            ElectromechanicalNode::new(0.0, 0.0, 0.0),
            ElectromechanicalNode::new(1.0, 0.0, 0.0),
            ElectromechanicalNode::new(0.0, 1.0, 0.0),
            ElectromechanicalNode::new(0.0, 0.0, 1.0),
        ];
        let mat = PiezoMaterial::isotropic_with_piezo(67e9, 0.31, 1200.0, 0.0, 0.0, 0.0);
        local.materials = vec![mat];
        local.elements = vec![ElectromechanicalElement {
            nodes: [0, 1, 2, 3],
            material_id: 0,
        }];
        let mut global = GlobalPiezoAssembler::new(local);
        global.assemble();
        let kpp = global.extract_k_pp();
        assert_eq!(kpp.len(), 4);
        assert_eq!(kpp[0].len(), 4);
    }
    /// Sensor short-circuit charge is positive for positive stress.
    #[test]
    fn test_sensor_short_circuit_charge() {
        let sensor = PiezoSensor::new(67e9, 0.31, 15.1, -5.2, 12.3, 1e-4, 1e-3, 1200.0, 1e6);
        let q = sensor.short_circuit_charge(1e6);
        assert!(q > 0.0, "short-circuit charge = {q}");
    }
    /// Sensor open-circuit voltage is positive for positive stress.
    #[test]
    fn test_sensor_open_circuit_voltage() {
        let sensor = PiezoSensor::new(67e9, 0.31, 15.1, -5.2, 12.3, 1e-4, 1e-3, 1200.0, 1e6);
        let v = sensor.open_circuit_voltage(1e6);
        assert!(v > 0.0, "open-circuit voltage = {v}");
    }
    /// Sensor capacitance is positive.
    #[test]
    fn test_sensor_capacitance() {
        let sensor = PiezoSensor::new(67e9, 0.31, 15.1, -5.2, 12.3, 1e-4, 1e-3, 1200.0, 1e6);
        let c = sensor.capacitance();
        assert!(c > 0.0 && c.is_finite(), "capacitance = {c}");
    }
    /// Power output increases with stress amplitude.
    #[test]
    fn test_sensor_power_increases_with_stress() {
        let sensor = PiezoSensor::new(67e9, 0.31, 15.1, -5.2, 12.3, 1e-4, 1e-3, 1200.0, 1e3);
        let p1 = sensor.power_output(1e6, 100.0);
        let p2 = sensor.power_output(2e6, 100.0);
        assert!(
            p2 > p1,
            "power output should increase with stress: p1={p1}, p2={p2}"
        );
    }
    /// Open-circuit voltage is proportional to stress and length.
    #[test]
    fn test_piezo_coupling_open_circuit_voltage_proportional() {
        let pc = PiezoCoupling::new(67e9, 0.31, 15.1, -5.2, 12.3);
        pub(super) const EPS0: f64 = 8.854_187_817e-12;
        let eps33 = 1200.0 * EPS0;
        let v1 = pc.open_circuit_voltage(1e6, 0.01, eps33);
        let v2 = pc.open_circuit_voltage(2e6, 0.01, eps33);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-10,
            "Voc should double when stress doubles: v1={v1}, v2={v2}"
        );
    }
    /// k33 is in physically reasonable range for PZT.
    #[test]
    fn test_piezo_coupling_k33_reasonable() {
        pub(super) const EPS0: f64 = 8.854_187_817e-12;
        let pc = PiezoCoupling::new(67e9, 0.31, 15.1, -5.2, 12.3);
        let eps33 = 1200.0 * EPS0;
        let k = pc.k33(eps33);
        assert!(k > 0.0 && k < 1.0, "k33 = {k}, expected (0,1)");
    }
    /// Electric field magnitude function.
    #[test]
    fn test_electric_field_magnitude_zero() {
        let solver = ElectrostaticSolver::new(4, 4, 1e-3);
        let mag = solver.electric_field_magnitude();
        for &m in &mag {
            assert!(m.abs() < 1e-30, "mag = {m}");
        }
    }
    /// Electric field solver converges with Dirichlet boundary.
    #[test]
    fn test_electrostatic_solver_dirichlet() {
        let mut solver = ElectrostaticSolver::new(5, 5, 1e-3);
        for ix in 0..5 {
            solver.phi[4 * 5 + ix] = 1.0;
        }
        let iters = solver.solve_laplace_gauss_seidel(50_000, 1e-6);
        assert!(iters <= 50_000, "should converge");
        let center = solver.phi[2 * 5 + 2];
        assert!(center > 0.1 && center < 0.9, "center phi = {center}");
    }
    /// Default assembler has empty collections.
    #[test]
    fn test_piezo_assembler_empty() {
        let asm = PiezoFemAssembler::new();
        assert!(asm.nodes.is_empty());
        assert!(asm.elements.is_empty());
        assert!(asm.materials.is_empty());
    }
    /// Element volume of a scaled tet.
    #[test]
    fn test_element_volume_scaled() {
        let nodes: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        let vol = PiezoFemAssembler::compute_element_volume(&nodes);
        assert!((vol - 8.0 / 6.0).abs() < 1e-10, "vol = {vol}");
    }
    /// Breakdown at exactly the threshold.
    #[test]
    fn test_dielectric_breakdown_at_threshold() {
        let bd = DielectricBreakdown::new(3e7, 3.5);
        assert!(bd.is_broken_down(3e7), "should break at exactly threshold");
    }
    /// Safety margin at zero field is infinite.
    #[test]
    fn test_dielectric_breakdown_zero_field() {
        let bd = DielectricBreakdown::new(3e7, 3.5);
        let m = bd.safety_margin(0.0);
        assert!(m.is_infinite(), "safety margin at zero field = {m}");
    }
    /// MEMS capacitance at rest gap equals nominal value.
    #[test]
    fn test_mems_capacitance_at_rest() {
        pub(super) const EPS0: f64 = 8.854_187_817e-12;
        let area = 1e-6;
        let gap = 2e-6;
        let eps_r = 1.0;
        let mems = MemsActuator::new(area, gap, 100.0, eps_r);
        let c = mems.capacitance(gap);
        let expected = eps_r * EPS0 * area / gap;
        assert!((c - expected).abs() < 1e-30, "C = {c}, expected {expected}");
    }
    /// Pull-in voltage for vacuum gap.
    #[test]
    fn test_mems_pull_in_analytical() {
        pub(super) const EPS0: f64 = 8.854_187_817e-12;
        let k = 100.0;
        let g0 = 2e-6;
        let area = 1e-6;
        let eps_r = 1.0;
        let mems = MemsActuator::new(area, g0, k, eps_r);
        let v_pi = mems.pull_in_voltage();
        let expected = ((8.0 * k * g0.powi(3)) / (27.0 * eps_r * EPS0 * area)).sqrt();
        assert!(
            (v_pi - expected).abs() / expected < 1e-6,
            "V_pi = {v_pi}, expected {expected}"
        );
    }
    /// B-matrix rows sum consistent with constant-field assumption.
    #[test]
    fn test_b_matrix_unit_tet() {
        let nodes = unit_tet();
        let b = PiezoFemAssembler::compute_b_matrix(&nodes);
        let sum_row0: f64 = b[0].iter().enumerate().step_by(3).map(|(_, &v)| v).sum();
        assert!(sum_row0.abs() < 1e-10, "B matrix row 0 sum = {sum_row0}");
    }
}
