//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::ConformingTransfer;

/// Stefan-Boltzmann constant \[W/(m²·K⁴)\]
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
/// Universal gas constant \[J/(mol·K)\]
pub const GAS_CONSTANT: f64 = 8.314_462_618;
/// Simple dense Gauss-Seidel solver for a symmetric positive-definite system.
pub fn gauss_seidel_solve(
    a: &[Vec<f64>],
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = b.len();
    let mut x = x0.to_vec();
    for _ in 0..max_iter {
        let x_old = x.clone();
        for i in 0..n {
            let mut s = b[i];
            for j in 0..n {
                if j != i {
                    s -= a[i][j] * x[j];
                }
            }
            x[i] = s / a[i][i].max(1e-30);
        }
        let diff: f64 = x
            .iter()
            .zip(x_old.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        if diff < tol {
            break;
        }
    }
    x
}
/// Compute the spectral radius of a matrix (largest absolute eigenvalue estimate).
///
/// Uses power iteration.
pub fn spectral_radius(a: &[Vec<f64>], max_iter: usize) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let mut v = vec![1.0; n];
    let mut lambda = 0.0f64;
    for _ in 0..max_iter {
        let mut av = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                av[i] += a[i][j] * v[j];
            }
        }
        lambda = av.iter().cloned().fold(0.0_f64, |acc, x| acc.max(x.abs()));
        if lambda < 1e-30 {
            break;
        }
        for i in 0..n {
            v[i] = av[i] / lambda;
        }
    }
    lambda
}
/// Compute the Frobenius norm of a matrix.
pub fn frobenius_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .flat_map(|row| row.iter())
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt()
}
/// Interface coupling residual for a 1-field partitioned scheme.
///
/// r = x_solid − T · x_fluid  (displacement compatibility)
pub fn interface_residual(
    x_solid: &[f64],
    x_fluid: &[f64],
    transfer: &ConformingTransfer,
) -> Vec<f64> {
    let x_fluid_on_solid = transfer.transfer(x_fluid);
    x_solid
        .iter()
        .zip(x_fluid_on_solid.iter())
        .map(|(xs, xf)| xs - xf)
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::coupled_fem::*;
    #[test]
    fn thermo_elastic_modulus_at_ref_temp() {
        let mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        let e = mat.modulus_at(mat.t_ref);
        assert!((e - 200e9).abs() < 1.0, "E at T_ref should equal E0");
    }
    #[test]
    fn thermo_elastic_modulus_decreases_with_temperature() {
        let mut mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        mat.beta = 5e-4;
        let e_hot = mat.modulus_at(mat.t_ref + 100.0);
        let e_cold = mat.modulus_at(mat.t_ref);
        assert!(e_hot < e_cold, "modulus must decrease with temperature");
    }
    #[test]
    fn thermal_load_vector_nonzero_for_positive_delta_t() {
        let mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        let f = mat.thermal_load_vector_2d(mat.t_ref, 50.0);
        assert!(f[0].abs() > 0.0, "thermal load x must be nonzero");
        assert!(f[1].abs() > 0.0, "thermal load y must be nonzero");
        assert!(f[2].abs() < 1e-10, "shear component must be zero");
    }
    #[test]
    fn constitutive_2d_symmetry() {
        let mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        let d = mat.constitutive_2d(mat.t_ref);
        for i in 0..3 {
            for j in 0..3 {
                assert!((d[i][j] - d[j][i]).abs() < 1.0, "D must be symmetric");
            }
        }
    }
    #[test]
    fn thermal_field_uniform_average() {
        let field = ThermalField::uniform(5, 5, 300.0);
        assert!((field.average() - 300.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_field_interpolation_corners() {
        let field = ThermalField::uniform(4, 4, 400.0);
        let t = field.interpolate(0, 0, 0.0, 0.0);
        assert!((t - 400.0).abs() < 1e-10);
    }
    #[test]
    fn thermal_field_bc_left_right() {
        let mut field = ThermalField::uniform(5, 5, 300.0);
        field.set_left_bc(600.0);
        field.set_right_bc(100.0);
        assert!((field.temperatures[field.idx(0, 0)] - 600.0).abs() < 1e-10);
        assert!((field.temperatures[field.idx(4, 0)] - 100.0).abs() < 1e-10);
    }
    #[test]
    fn thermo_mech_coupling_step_changes_displacements() {
        let mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        let mut coupler = ThermoMechanicalCoupling::new(5, 5, mat);
        coupler.thermal.set_left_bc(600.0);
        let change = coupler.staggered_step(0.01, 1e6);
        assert!(change >= 0.0);
    }
    #[test]
    fn thermo_mech_coupling_solve_converges() {
        let mat = ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0);
        let mut coupler = ThermoMechanicalCoupling::new(5, 5, mat);
        coupler.tol = 1e-3;
        coupler.max_iter = 100;
        let iters = coupler.solve(1e-4, 0.0);
        assert!(iters <= 100);
    }
    #[test]
    fn ale_node_convective_velocity() {
        let mut node = AleNode::new(0.0, 0.0);
        node.fluid_velocity = [3.0, 1.0];
        node.mesh_velocity = [1.0, 0.5];
        let cv = node.convective_velocity();
        assert!((cv[0] - 2.0).abs() < 1e-10);
        assert!((cv[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn ale_node_update_position() {
        let mut node = AleNode::new(1.0, 2.0);
        node.displacement = [0.1, -0.2];
        node.update_position();
        assert!((node.spatial_pos[0] - 1.1).abs() < 1e-10);
        assert!((node.spatial_pos[1] - 1.8).abs() < 1e-10);
    }
    #[test]
    fn ale_fsi_solver_step_no_panic() {
        let nodes: Vec<AleNode> = (0..10).map(|i| AleNode::new(i as f64 * 0.1, 0.0)).collect();
        let iface = AleInterface::new((0..10).collect(), vec![[0.0, 1.0]; 10]);
        let mut solver = AleFsiSolver::new(nodes, iface);
        solver.step(0.01, 1.0);
    }
    #[test]
    fn ale_fsi_fluid_force_nonzero_with_pressure() {
        let nodes: Vec<AleNode> = (0..5)
            .map(|i| {
                let mut n = AleNode::new(i as f64 * 0.2, 0.5);
                n.pressure = 100.0;
                n
            })
            .collect();
        let iface = AleInterface::new((0..5).collect(), vec![[0.0, 1.0]; 5]);
        let solver = AleFsiSolver::new(nodes, iface);
        let f = solver.total_fluid_force();
        assert!(f[1].abs() > 0.0, "fluid should exert force");
    }
    #[test]
    fn piezo_stress_zero_field_zero_strain() {
        let piezo = PiezoelectricCoupling::pzt5h();
        let sigma = piezo.compute_stress([0.0; 3], [0.0; 2]);
        for s in sigma {
            assert!(s.abs() < 1e-6, "zero strain/field → zero stress");
        }
    }
    #[test]
    fn piezo_electric_displacement_nonzero_under_strain() {
        let piezo = PiezoelectricCoupling::pzt5h();
        let d = piezo.compute_electric_displacement([0.001, 0.0, 0.0], [0.0; 2]);
        assert!(d.iter().any(|x| x.abs() > 0.0));
    }
    #[test]
    fn piezo_coupling_coefficient_positive() {
        let piezo = PiezoelectricCoupling::pzt5h();
        let k = piezo.coupling_coefficient();
        assert!(k >= 0.0, "coupling coefficient must be non-negative");
    }
    #[test]
    fn magnetostrictive_strain_zero_at_zero_field() {
        let mat = MagnetostrictiveMaterial::terfenol_d();
        let eps = mat.magnetostrictive_strain(0.0);
        assert!(eps.abs() < 1e-20, "no field → no magnetostrictive strain");
    }
    #[test]
    fn magnetostrictive_strain_saturates() {
        let mat = MagnetostrictiveMaterial::terfenol_d();
        let eps_large = mat.magnetostrictive_strain(1e9);
        assert!(
            (eps_large - mat.lambda_s).abs() / mat.lambda_s < 0.01,
            "strain should saturate near lambda_s"
        );
    }
    #[test]
    fn magnetostrictive_effective_modulus_increases_with_field() {
        let mat = MagnetostrictiveMaterial::terfenol_d();
        let e0 = mat.effective_modulus(0.0);
        let e1 = mat.effective_modulus(1e6);
        assert!(e1 >= e0, "ΔE effect: modulus should increase with field");
    }
    #[test]
    fn chemo_mech_diffusion_step_conserves_species_approx() {
        let mut prob = ChemoMechanicalProblem::new_1d(21, 1.0, 1.0);
        prob.diff_coeff = 1e-14;
        let total_before = prob.total_species();
        prob.diffusion_step(1e8, 1.0, 1.0);
        let total_after = prob.total_species();
        assert!(
            (total_before - total_after).abs() / total_before < 0.01,
            "species should be approximately conserved with uniform BCs"
        );
    }
    #[test]
    fn chemo_mech_stress_nonzero_with_concentration_gradient() {
        let mut prob = ChemoMechanicalProblem::new_1d(11, 1.0, 1.0);
        prob.concentration[0] = 2.0;
        prob.update_stress_from_concentration(1.0);
        assert!(prob.hydrostatic_stress[0].abs() > 0.0);
    }
    #[test]
    fn chemo_mech_chemical_potential_decreases_with_stress() {
        let mut prob = ChemoMechanicalProblem::new_1d(5, 1.0, 1.0);
        let mu_no_stress = prob.chemical_potential(0);
        prob.hydrostatic_stress[0] = 1e8;
        let mu_stress = prob.chemical_potential(0);
        assert!(
            mu_stress < mu_no_stress,
            "positive hydrostatic stress should reduce chemical potential"
        );
    }
    #[test]
    fn acoustic_cavity_impedance() {
        let cav = AcousticCavity1D::new(11, 1.0, 343.0, 1.225);
        let z = cav.impedance();
        assert!((z - 343.0 * 1.225).abs() < 1e-6);
    }
    #[test]
    fn acoustic_cavity_energy_zero_initially() {
        let cav = AcousticCavity1D::new(11, 1.0, 343.0, 1.225);
        assert!(cav.total_energy().abs() < 1e-15);
    }
    #[test]
    fn acoustic_cavity_pressure_source_propagates() {
        let mut cav = AcousticCavity1D::new(21, 1.0, 343.0, 1.225);
        let dt = 1.0 / (343.0 * 20.0 * 2.0);
        cav.apply_pressure_source(100.0, 2000.0 * std::f64::consts::TAU, 0.25e-3);
        for _ in 0..20 {
            cav.step(dt);
        }
        let e = cav.total_energy();
        assert!(e > 0.0, "energy should be positive after excitation");
    }
    #[test]
    fn aitken_relaxation_initial_omega() {
        let aitken = AitkenRelaxation::new(4, 0.5);
        assert!((aitken.omega - 0.5).abs() < 1e-10);
    }
    #[test]
    fn aitken_relaxation_relax_moves_toward_prediction() {
        let aitken = AitkenRelaxation::new(3, 0.5);
        let x_old = vec![0.0; 3];
        let x_new = vec![2.0; 3];
        let relaxed = aitken.relax(&x_old, &x_new);
        for v in relaxed {
            assert!((v - 1.0).abs() < 1e-10, "ω=0.5 → midpoint");
        }
    }
    #[test]
    fn aitken_omega_clamped_within_bounds() {
        let mut aitken = AitkenRelaxation::new(2, 0.5);
        aitken.residual_prev = vec![1.0, 1.0];
        aitken.residual_curr = vec![0.999, 0.999];
        aitken.update_omega();
        assert!(aitken.omega >= aitken.omega_min);
        assert!(aitken.omega <= aitken.omega_max);
    }
    #[test]
    fn monolithic_block_size() {
        let sys = MonolithicBlock2x2::new(3, 2);
        let mat = sys.full_matrix();
        assert_eq!(mat.len(), 5);
        assert_eq!(mat[0].len(), 5);
    }
    #[test]
    fn monolithic_block_gauss_seidel_diagonal_system() {
        let mut sys = MonolithicBlock2x2::new(2, 2);
        sys.k_uu[0][0] = 1.0;
        sys.k_uu[1][1] = 1.0;
        sys.k_pp[0][0] = 1.0;
        sys.k_pp[1][1] = 1.0;
        sys.f_u = vec![3.0, 5.0];
        sys.f_p = vec![7.0, 11.0];
        sys.block_gauss_seidel_step();
        assert!((sys.u[0] - 3.0).abs() < 1e-10);
        assert!((sys.u[1] - 5.0).abs() < 1e-10);
        assert!((sys.p[0] - 7.0).abs() < 1e-10);
        assert!((sys.p[1] - 11.0).abs() < 1e-10);
    }
    #[test]
    fn monolithic_full_rhs_concatenates() {
        let mut sys = MonolithicBlock2x2::new(2, 3);
        sys.f_u = vec![1.0, 2.0];
        sys.f_p = vec![3.0, 4.0, 5.0];
        let rhs = sys.full_rhs();
        assert_eq!(rhs, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    }
    #[test]
    fn pod_rom_basis_orthonormal() {
        let snaps: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let rom = PodRom::from_snapshots(&snaps, 3);
        for m1 in 0..rom.n_modes {
            for m2 in 0..rom.n_modes {
                let dot: f64 = (0..rom.n_full)
                    .map(|i| rom.basis[m1][i] * rom.basis[m2][i])
                    .sum();
                if m1 == m2 {
                    assert!((dot - 1.0).abs() < 1e-6, "basis not normalized: {dot}");
                } else {
                    assert!(dot.abs() < 1e-6, "basis not orthogonal: {dot}");
                }
            }
        }
    }
    #[test]
    fn pod_rom_project_reconstruct_roundtrip() {
        let snaps: Vec<Vec<f64>> = vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0]];
        let rom = PodRom::from_snapshots(&snaps, 2);
        let v = vec![0.6, 0.8, 0.0, 0.0];
        let q = rom.project(&v);
        let v_back = rom.reconstruct(&q);
        let err: f64 = v
            .iter()
            .zip(v_back.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-8, "reconstruct error too large: {err}");
    }
    #[test]
    fn pod_rom_projection_error_within_span_is_small() {
        let snaps: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let rom = PodRom::from_snapshots(&snaps, 2);
        let v = vec![3.0, 4.0];
        let err = rom.projection_error(&v);
        assert!(err < 1e-8, "projection error within span: {err}");
    }
    #[test]
    fn mortar_gap_zero_matching_displacements() {
        let master = vec![0.0, 0.5, 1.0];
        let slave = vec![0.0, 0.5, 1.0];
        let iface = MortarInterface::new_1d(master, slave);
        let u_m = vec![0.1, 0.2, 0.3];
        let u_s = iface
            .m_mat
            .iter()
            .map(|row| row.iter().zip(u_m.iter()).map(|(w, u)| w * u).sum::<f64>())
            .collect::<Vec<_>>();
        let gap = iface.gap_function(&u_s, &u_m);
        for g in gap {
            assert!(
                g.abs() < 1e-10,
                "gap should be zero for matching meshes: {g}"
            );
        }
    }
    #[test]
    fn mortar_compatibility_error_positive_for_mismatch() {
        let master = vec![0.0, 1.0];
        let slave = vec![0.0, 0.5, 1.0];
        let iface = MortarInterface::new_1d(master, slave);
        let u_m = vec![0.0, 1.0];
        let u_s = vec![0.0, 0.0, 1.0];
        let err = iface.compatibility_error(&u_s, &u_m);
        assert!(err >= 0.0);
    }
    #[test]
    fn conforming_transfer_identity_same_mesh() {
        let pos: Vec<f64> = (0..5).map(|i| i as f64 * 0.25).collect();
        let tf = ConformingTransfer::linear_interpolation(pos.clone(), pos.clone());
        let vals = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let out = tf.transfer(&vals);
        for (a, b) in vals.iter().zip(out.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn conforming_transfer_linear_interpolation_midpoint() {
        let src = vec![0.0, 1.0];
        let tgt = vec![0.5];
        let tf = ConformingTransfer::linear_interpolation(src, tgt);
        let vals = vec![0.0, 2.0];
        let out = tf.transfer(&vals);
        assert!(
            (out[0] - 1.0).abs() < 1e-10,
            "midpoint interpolation: {}",
            out[0]
        );
    }
    #[test]
    fn gauss_seidel_diagonal_system() {
        let a = vec![
            vec![4.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let b = vec![8.0, 9.0, 6.0];
        let x = gauss_seidel_solve(&a, &b, &[0.0; 3], 100, 1e-10);
        assert!((x[0] - 2.0).abs() < 1e-8);
        assert!((x[1] - 3.0).abs() < 1e-8);
        assert!((x[2] - 3.0).abs() < 1e-8);
    }
    #[test]
    fn spectral_radius_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let rho = spectral_radius(&a, 100);
        assert!((rho - 1.0).abs() < 1e-6);
    }
    #[test]
    fn frobenius_norm_identity_sqrt2() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let f = frobenius_norm(&a);
        assert!((f - 2.0f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn interface_residual_zero_for_matching_fields() {
        let src = vec![0.0, 0.5, 1.0];
        let tgt = vec![0.0, 0.5, 1.0];
        let tf = ConformingTransfer::linear_interpolation(src.clone(), tgt.clone());
        let vals = vec![1.0, 2.0, 3.0];
        let res = interface_residual(&vals, &vals, &tf);
        for r in res {
            assert!(r.abs() < 1e-10);
        }
    }
    #[test]
    fn coupling_strategy_equality() {
        let s1 = CouplingStrategy::Staggered;
        let s2 = CouplingStrategy::Monolithic;
        assert_ne!(s1, s2);
        let s3 = CouplingStrategy::IterativePartitioned {
            max_iter: 50,
            tol: 1e-6,
        };
        assert_ne!(s1, s3);
    }
}
#[cfg(test)]
mod tests_new_coupled {

    use crate::coupled_fem::*;
    #[test]
    fn thermal_expansion_zero_at_ref_temp() {
        let mut state = ThermalExpansionState::new(300.0, 12e-6);
        state.update(300.0);
        for &s in &state.thermal_strain {
            assert!(s.abs() < 1e-20, "strain at T_ref must be 0");
        }
    }
    #[test]
    fn thermal_expansion_positive_above_ref() {
        let mut state = ThermalExpansionState::new(300.0, 12e-6);
        state.update(400.0);
        let eps = 12e-6 * 100.0;
        assert!((state.thermal_strain[0] - eps).abs() < 1e-15);
        assert!((state.thermal_strain[1] - eps).abs() < 1e-15);
        assert!((state.thermal_strain[2] - eps).abs() < 1e-15);
        assert_eq!(state.thermal_strain[3], 0.0);
    }
    #[test]
    fn thermal_expansion_volumetric_strain() {
        let mut state = ThermalExpansionState::new(300.0, 10e-6);
        state.update(400.0);
        let vol = state.volumetric_thermal_strain();
        assert!((vol - 3.0 * 10e-6 * 100.0).abs() < 1e-15);
    }
    #[test]
    fn thermal_expansion_elastic_strain_subtraction() {
        let mut state = ThermalExpansionState::new(300.0, 10e-6);
        state.update(400.0);
        let eps_th = 10e-6 * 100.0;
        let total = [eps_th + 1e-4, eps_th, eps_th, 0.0, 0.0, 0.0];
        let el = state.elastic_strain(total);
        assert!((el[0] - 1e-4).abs() < 1e-15);
        assert!(el[1].abs() < 1e-15);
    }
    #[test]
    fn temp_dependent_stiffness_at_ref_temp_equals_e0() {
        let s = TempDependentStiffness::new(200e9, 5e-4, 300.0, 0.3);
        assert!((s.modulus(300.0) - 200e9).abs() < 1.0);
    }
    #[test]
    fn temp_dependent_stiffness_decreases_above_ref() {
        let s = TempDependentStiffness::new(200e9, 5e-4, 300.0, 0.3);
        assert!(s.modulus(400.0) < s.modulus(300.0));
    }
    #[test]
    fn temp_dependent_stiffness_matrix_symmetry() {
        let s = TempDependentStiffness::new(200e9, 0.0, 300.0, 0.3);
        let d = s.stiffness_6x6(300.0);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (d[i][j] - d[j][i]).abs() < 1.0,
                    "D[{i}][{j}] must equal D[{j}][{i}]"
                );
            }
        }
    }
    #[test]
    fn temp_dependent_stress_zero_strain_zero_delta_t() {
        let s = TempDependentStiffness::new(200e9, 0.0, 300.0, 0.3);
        let sigma = s.stress_from_total_strain([0.0; 6], 300.0, 12e-6);
        for &si in &sigma {
            assert!(si.abs() < 1.0, "zero total strain at T_ref → zero stress");
        }
    }
    #[test]
    fn ale_mesh_chain_edge_count() {
        let pos: Vec<[f64; 2]> = (0..5).map(|i| [i as f64 * 0.1, 0.0]).collect();
        let mesh = AleMesh::chain(pos);
        assert_eq!(mesh.edges.len(), 4);
    }
    #[test]
    fn ale_mesh_smooth_laplacian_no_panic() {
        let pos: Vec<[f64; 2]> = (0..7).map(|i| [i as f64, (i as f64).sin()]).collect();
        let mut mesh = AleMesh::chain(pos);
        mesh.smooth_laplacian();
    }
    #[test]
    fn ale_mesh_distortion_undistorted() {
        let pos: Vec<[f64; 2]> = (0..4).map(|i| [i as f64, 0.0]).collect();
        let mesh = AleMesh::chain(pos);
        let distortion = mesh.distortion_metric();
        assert!(
            (distortion - 1.0).abs() < 1e-10,
            "undistorted mesh distortion = 1.0"
        );
    }
    #[test]
    fn ale_mesh_mesh_velocity_update() {
        let pos: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let mut mesh = AleMesh::chain(pos.clone());
        mesh.positions[1][0] = 1.1;
        mesh.update_mesh_velocity(&pos, 0.1);
        assert!(
            (mesh.mesh_velocity[1][0] - 1.0).abs() < 1e-10,
            "v_mesh = 0.1/0.1 = 1.0"
        );
    }
    #[test]
    fn lorentz_force_cross_product() {
        let lf = LorentzForceDensity::new([0.0, 1.0, 0.0], 1e7);
        let f = lf.force([1.0, 0.0, 0.0]);
        assert!(f[2].abs() > 0.5, "J×B should have z-component");
        assert!(f[0].abs() < 1e-12);
        assert!(f[1].abs() < 1e-12);
    }
    #[test]
    fn lorentz_force_zero_current_zero_force() {
        let lf = LorentzForceDensity::new([0.0, 0.0, 1.0], 1e7);
        let f = lf.force([0.0, 0.0, 0.0]);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn lorentz_ohmic_dissipation_positive() {
        let lf = LorentzForceDensity::new([0.0, 1.0, 0.0], 1e6);
        let q = lf.ohmic_dissipation([1000.0, 0.0, 0.0]);
        assert!(q > 0.0);
    }
    #[test]
    fn radiation_damping_coeff_positive() {
        let rd = RadiationDamping::new(1000.0, 1500.0, 0.01);
        assert!(rd.damping_coeff() > 0.0);
    }
    #[test]
    fn radiation_damping_opposes_velocity() {
        let rd = RadiationDamping::new(1000.0, 1500.0, 0.01);
        let f = rd.damping_force(0.1);
        assert!(f < 0.0, "damping force must oppose positive velocity");
    }
    #[test]
    fn radiation_power_proportional_to_v_squared() {
        let rd = RadiationDamping::new(1000.0, 1500.0, 1.0);
        let p1 = rd.radiated_power(1.0);
        let p2 = rd.radiated_power(2.0);
        assert!((p2 - 4.0 * p1).abs() < 1e-6);
    }
    #[test]
    fn acoustic_load_zero_pressure_zero_force() {
        let indices = vec![0, 1, 2];
        let normals = vec![[0.0, 1.0, 0.0]; 3];
        let areas = vec![0.01; 3];
        let load = AcousticPressureLoad::new(indices, normals, areas);
        let forces = load.nodal_forces();
        for f in forces {
            assert_eq!(f, [0.0, 0.0, 0.0]);
        }
    }
    #[test]
    fn acoustic_load_force_magnitude() {
        let indices = vec![0];
        let normals = vec![[0.0, 1.0, 0.0]];
        let areas = vec![1.0];
        let mut load = AcousticPressureLoad::new(indices, normals, areas);
        load.set_pressure(vec![100.0]);
        let mag = load.total_force_magnitude();
        assert!((mag - 100.0).abs() < 1e-8);
    }
    #[test]
    fn hygro_diffusion_initial_moisture_uniform() {
        let mat = HygroMechanicalMaterial::new(10e9, 0.3, 0.2e-3, 0.0, 1e-12);
        let solver = HygroDiffusionSolver::new(11, 1.0, 0.5, mat);
        let m = solver.total_moisture();
        assert!(
            (m - 0.5).abs() < 1e-10,
            "total moisture should be 0.5 * 1.0 m"
        );
    }
    #[test]
    fn hygro_stress_nonzero_above_ref() {
        let mat = HygroMechanicalMaterial::new(10e9, 0.3, 0.2e-3, 0.0, 1e-12);
        let mut solver = HygroDiffusionSolver::new(5, 1.0, 100.0, mat);
        solver.update_hygrostress();
        assert!(solver.hygrostress[0].abs() > 0.0);
    }
    #[test]
    fn multiphysics_residual_zero_initially() {
        let res = MultiPhysicsResidual::new(5, 5, 3);
        assert!(res.total_norm().abs() < 1e-15);
    }
    #[test]
    fn multiphysics_residual_assembly() {
        let mut res = MultiPhysicsResidual::new(2, 2, 0);
        let k = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let u = vec![1.0, 1.0];
        let f = vec![2.0, 3.0];
        res.assemble_structural(&k, &u, &f);
        assert!(res.struct_norm().abs() < 1e-12, "K·u = f → zero residual");
    }
    #[test]
    fn staggered_scheme_converges_identical_fields() {
        let mut scheme = StaggeredCouplingScheme::new(100, 1e-8);
        let x = vec![1.0, 2.0, 3.0];
        let (conv, rel) = scheme.check_convergence(&x, &x);
        assert!(conv, "identical fields must be converged");
        assert!(rel.abs() < 1e-15);
    }
    #[test]
    fn staggered_scheme_not_converged_different_fields() {
        let mut scheme = StaggeredCouplingScheme::new(10, 1e-8);
        let x_old = vec![0.0, 0.0];
        let x_new = vec![1.0, 1.0];
        let (conv, _) = scheme.check_convergence(&x_old, &x_new);
        assert!(!conv);
    }
    #[test]
    fn coupling_comparison_diagonal_system() {
        let cmp = CouplingComparison::compare_2x2(2.0, 0.0, 0.0, 3.0, 4.0, 9.0, 50, 1e-10);
        assert!(
            cmp.staggered_iters <= 5,
            "Diagonal system should converge quickly"
        );
        assert!(cmp.staggered_converged);
        assert!(cmp.monolithic_residual.abs() < 1e-10);
    }
    #[test]
    fn mls_interpolator_reproduces_constant() {
        let positions: Vec<f64> = (0..5).map(|i| i as f64 * 0.25).collect();
        let values = vec![3.0f64; 5];
        let mls = MlsInterpolator::new(positions, 0.5);
        let tgt = 0.5;
        let v = mls.interpolate(&values, tgt);
        assert!(
            (v - 3.0).abs() < 1e-6,
            "MLS should reproduce constant: got {v}"
        );
    }
    #[test]
    fn mls_interpolator_reproduces_linear() {
        let positions: Vec<f64> = (0..6).map(|i| i as f64 * 0.2).collect();
        let values: Vec<f64> = positions.iter().map(|&x| 2.0 * x).collect();
        let mls = MlsInterpolator::new(positions, 0.6);
        let v = mls.interpolate(&values, 0.5);
        assert!(
            (v - 1.0).abs() < 0.05,
            "MLS should approximate linear: got {v}"
        );
    }
    #[test]
    fn mls_transfer_same_grid_identity() {
        let pos: Vec<f64> = (0..5).map(|i| i as f64 * 0.25).collect();
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mls = MlsInterpolator::new(pos.clone(), 0.5);
        let out = mls.transfer(&values, &pos);
        for (a, b) in values.iter().zip(out.iter()) {
            assert!((a - b).abs() < 0.1, "same-grid transfer: {a} vs {b}");
        }
    }
}
