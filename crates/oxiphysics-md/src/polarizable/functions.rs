//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{InducedDipole, InducibleDipoleForce};

#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn norm2(a: [f64; 3]) -> f64 {
    dot(a, a)
}
#[inline]
pub(super) fn norm(a: [f64; 3]) -> f64 {
    norm2(a).sqrt()
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale(s: f64, a: [f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Electric field at position `r` due to a point dipole `dipole`.
///
/// E(r) = (1/r^3) * (3*(p·r̂)*r̂ - p)
pub fn dipole_field_at(dipole: &InducedDipole, r: [f64; 3]) -> [f64; 3] {
    InducibleDipoleForce::dipole_electric_field(dipole.dipole, sub(r, dipole.position))
}
/// Polarization energy: U_pol = -0.5 * sum_i (p_i · p_i / alpha_i).
///
/// For linear response (p = alpha * E_ind), this equals -0.5 * sum alpha_i * |E_ind_i|^2.
pub fn polarization_energy(dipoles: &[[f64; 3]], polarizabilities: &[f64]) -> f64 {
    let mut u = 0.0f64;
    for (i, (&p, &alpha)) in dipoles.iter().zip(polarizabilities.iter()).enumerate() {
        let _ = i;
        if alpha > 1e-30 {
            u -= 0.5 * dot(p, p) / alpha;
        }
    }
    u
}
/// Self-consistent field (SCF) loop to solve for induced dipoles.
///
/// Iterates until the RMS dipole change is less than `tol` or `max_iter`
/// steps have been taken.  Dipole-dipole interactions between sites are
/// included; Coulomb fields from fixed charges use a simple 1/r^2 model.
///
/// Returns the final induced dipole moments for each site.
#[allow(clippy::too_many_arguments)]
pub fn self_consistent_dipoles(
    positions: &[[f64; 3]],
    polarizabilities: &[f64],
    charges: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut dipoles = vec![[0.0f64; 3]; n];
    for _iter in 0..max_iter {
        let old_dipoles = dipoles.clone();
        let mut new_dipoles = vec![[0.0f64; 3]; n];
        for i in 0..n {
            let mut e_field = [0.0f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r_vec = sub(positions[i], positions[j]);
                let r2 = norm2(r_vec);
                if r2 < 1e-24 {
                    continue;
                }
                let r3 = r2 * r2.sqrt();
                let q = charges[j];
                for d in 0..3 {
                    e_field[d] += q * r_vec[d] / r3;
                }
                let e_dip = InducibleDipoleForce::dipole_electric_field(old_dipoles[j], r_vec);
                e_field = add(e_field, e_dip);
            }
            new_dipoles[i] = scale(polarizabilities[i], e_field);
        }
        let mut rms = 0.0f64;
        for i in 0..n {
            let d = sub(new_dipoles[i], old_dipoles[i]);
            rms += norm2(d);
        }
        rms = (rms / (3 * n) as f64).sqrt();
        dipoles = new_dipoles;
        if rms < tol {
            break;
        }
    }
    dipoles
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_spring_energy_zero_at_rest() {
        let mut model = DrudeModel::new(500.0, 2.6);
        model.add_atom([0.0, 0.0, 0.0], 0.0, 12.0, 1.5);
        model.add_atom([3.0, 0.0, 0.0], -0.5, 16.0, 0.8);
        let e = model.spring_energy();
        assert!(
            e.abs() < 1e-15,
            "spring energy at rest should be 0, got {e}"
        );
    }
    #[test]
    fn test_drude_charge_from_alpha() {
        let k = 500.0_f64;
        let alpha = 1.5_f64;
        let q = DrudeModel::drude_charge_from_alpha(alpha, k);
        let expected = (k * alpha).sqrt();
        assert!((q - expected).abs() < 1e-12, "q = {q}, expected {expected}");
        assert!(((q * q) - k * alpha).abs() < 1e-10);
    }
    #[test]
    fn test_thole_screening_at_zero() {
        let s = TholeDamping::thole_screening(0.0, 2.6);
        assert!(s.abs() < 1e-15, "s(0) should be 0, got {s}");
    }
    #[test]
    fn test_thole_screening_large_u() {
        let s = TholeDamping::thole_screening(100.0, 2.6);
        assert!((s - 1.0).abs() < 1e-10, "s(large) should be ~1, got {s}");
    }
    #[test]
    fn test_point_dipole_single_site_converges_one_iter() {
        let mut model = PointDipoleModel::new();
        model.add_site([0.0, 0.0, 0.0], 0.0, 1.0);
        let e_ext = vec![[1.0, 0.0, 0.0]];
        let iters = model.solve_self_consistent(&e_ext, 100, 1e-12);
        assert!(
            iters <= 2,
            "single site should converge in at most 2 iterations, got {iters}"
        );
        let mu = model.dipoles[0];
        assert!(
            (mu[0] - 1.0).abs() < 1e-12,
            "mu_x should be alpha * E_x = 1.0, got {}",
            mu[0]
        );
    }
    #[test]
    fn test_polarization_energy_negative_aligned() {
        let mut model = PointDipoleModel::new();
        model.add_site([0.0, 0.0, 0.0], 0.0, 2.0);
        let e_ext = vec![[1.0, 0.0, 0.0]];
        model.solve_self_consistent(&e_ext, 100, 1e-12);
        let u = model.polarization_energy(&e_ext);
        assert!(u < 0.0, "polarization energy should be negative, got {u}");
    }
    #[test]
    fn test_dipole_electric_field_falloff() {
        let mu = [0.0, 0.0, 1.0];
        let r1 = [0.0, 0.0, 2.0_f64];
        let r2 = [0.0, 0.0, 4.0_f64];
        let e1 = InducibleDipoleForce::dipole_electric_field(mu, r1);
        let e2 = InducibleDipoleForce::dipole_electric_field(mu, r2);
        let mag1 = norm(e1);
        let mag2 = norm(e2);
        let ratio = mag2 / mag1;
        assert!(
            (ratio - 0.125).abs() < 1e-10,
            "field ratio should be 0.125 (1/r^3 falloff), got {ratio}"
        );
    }
    #[test]
    fn test_drude_displace_and_spring_energy() {
        let mut model = DrudeModel::new(100.0, 2.6);
        model.add_atom([0.0, 0.0, 0.0], 0.0, 12.0, 1.0);
        model.displace_drude(0, [0.1, 0.0, 0.0]);
        let e = model.spring_energy();
        assert!(
            (e - 0.5).abs() < 1e-12,
            "spring energy after displacement: {e}"
        );
    }
    #[test]
    fn test_drude_compute_dipoles() {
        let mut model = DrudeModel::new(100.0, 2.6);
        model.add_atom([0.0, 0.0, 0.0], 0.0, 12.0, 1.0);
        model.displace_drude(0, [0.5, 0.0, 0.0]);
        model.compute_dipoles();
        let mu = model.atoms[0].dipole;
        assert!((mu[0] - 5.0).abs() < 1e-10, "mu_x: {}", mu[0]);
    }
    #[test]
    fn test_drude_total_dipole() {
        let mut model = DrudeModel::new(100.0, 2.6);
        model.add_atom([0.0, 0.0, 0.0], 0.0, 12.0, 1.0);
        model.add_atom([3.0, 0.0, 0.0], 0.0, 12.0, 1.0);
        model.displace_drude(0, [0.1, 0.0, 0.0]);
        model.displace_drude(1, [-0.1, 0.0, 0.0]);
        model.compute_dipoles();
        let total = model.total_dipole();
        assert!(total[0].abs() < 1e-10, "total dipole x should be ~0");
    }
    #[test]
    fn test_drude_minimize() {
        let mut model = DrudeModel::new(100.0, 2.6);
        model.add_atom([0.0, 0.0, 0.0], 0.0, 12.0, 1.0);
        let e_ext = [[1.0, 0.0, 0.0]];
        let iters = model.minimize_drude(&e_ext, 1000, 0.001, 1e-6);
        assert!(iters < 1000, "should converge, took {iters}");
        let d = sub(
            model.atoms[0].drude.as_ref().unwrap().pos,
            model.atoms[0].pos,
        );
        assert!((d[0] - 0.1).abs() < 0.01, "displacement x: {}", d[0]);
    }
    #[test]
    fn test_shell_model_spring_energy_zero() {
        let mut model = ShellModel::new();
        model.add_atom([0.0, 0.0, 0.0], 1.0, -1.0, 16.0, 100.0);
        let e = model.spring_energy();
        assert!(e.abs() < 1e-15, "spring energy at rest: {e}");
    }
    #[test]
    fn test_shell_model_polarizability() {
        let mut model = ShellModel::new();
        model.add_atom([0.0, 0.0, 0.0], 1.0, -2.0, 16.0, 100.0);
        assert!(
            (model.atoms[0].polarizability - 0.04).abs() < 1e-12,
            "polarizability: {}",
            model.atoms[0].polarizability
        );
    }
    #[test]
    fn test_shell_model_induced_dipole() {
        let mut model = ShellModel::new();
        model.add_atom([0.0, 0.0, 0.0], 1.0, -2.0, 16.0, 100.0);
        model.atoms[0].shell_pos = [0.1, 0.0, 0.0];
        let mu = model.induced_dipole(0);
        assert!((mu[0] - (-0.2)).abs() < 1e-12, "mu_x: {}", mu[0]);
    }
    #[test]
    fn test_shell_model_relax() {
        let mut model = ShellModel::new();
        model.add_atom([0.0, 0.0, 0.0], 1.0, -2.0, 16.0, 100.0);
        let e_ext = [[10.0, 0.0, 0.0]];
        let iters = model.relax_shells(&e_ext, 1000, 0.001, 1e-6);
        assert!(iters < 1000, "should converge");
        let d = model.displacement(0);
        assert!((d[0] - (-0.2)).abs() < 0.01, "displacement: {}", d[0]);
    }
    #[test]
    fn test_cos_spring_energy_zero() {
        let mut cos = ChargeOnSpring::new(1.0);
        cos.add_site([0.0, 0.0, 0.0], 1.0, 100.0);
        assert!(cos.spring_energy().abs() < 1e-15);
    }
    #[test]
    fn test_cos_effective_polarizability() {
        let mut cos = ChargeOnSpring::new(1.0);
        cos.add_site([0.0, 0.0, 0.0], 5.0, 50.0);
        let alpha = cos.effective_polarizability(0);
        assert!((alpha - 0.5).abs() < 1e-12, "alpha: {alpha}");
    }
    #[test]
    fn test_cos_update_positions() {
        let mut cos = ChargeOnSpring::new(1.0);
        cos.add_site([0.0, 0.0, 0.0], 2.0, 100.0);
        let e_ext = [[1.0, 0.0, 0.0]];
        cos.update_positions(&e_ext);
        let d = sub(cos.spring_positions[0], cos.positions[0]);
        assert!((d[0] - 0.02).abs() < 1e-12, "d_x: {}", d[0]);
    }
    #[test]
    fn test_sor_single_site() {
        let mut model = PointDipoleModel::new();
        model.add_site([0.0, 0.0, 0.0], 0.0, 1.0);
        let e_ext = vec![[1.0, 0.0, 0.0]];
        let iters = model.solve_sor(&e_ext, 100, 1e-12, 1.0);
        assert!(
            iters <= 2,
            "SOR single site should converge quickly, got {iters}"
        );
        assert!(
            (model.dipoles[0][0] - 1.0).abs() < 1e-12,
            "mu_x: {}",
            model.dipoles[0][0]
        );
    }
    #[test]
    fn test_total_dipole_model() {
        let mut model = PointDipoleModel::new();
        model.add_site([0.0, 0.0, 0.0], 0.0, 1.0);
        model.add_site([5.0, 0.0, 0.0], 0.0, 1.0);
        let e_ext = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        model.solve_self_consistent(&e_ext, 100, 1e-12);
        let total = model.total_dipole();
        assert!(total[0] > 0.0, "total dipole should be positive");
    }
    #[test]
    fn test_dipole_dipole_force_symmetry() {
        let mu_i = [1.0, 0.0, 0.0];
        let mu_j = [0.0, 1.0, 0.0];
        let r = [2.0, 0.0, 0.0];
        let f_ij = InducibleDipoleForce::dipole_dipole_force(mu_i, mu_j, r);
        let f_ji = InducibleDipoleForce::dipole_dipole_force(mu_j, mu_i, scale(-1.0, r));
        for k in 0..3 {
            assert!(
                (f_ij[k] + f_ji[k]).abs() < 1e-10,
                "Newton's 3rd law violated: f_ij[{k}]={}, f_ji[{k}]={}",
                f_ij[k],
                f_ji[k]
            );
        }
    }
    #[test]
    fn test_thole_linear_at_zero() {
        let s = TholeDamping::thole_screening_linear(0.0, 2.6);
        assert!(s.abs() < 1e-15, "linear screening at 0: {s}");
    }
    #[test]
    fn test_thole_linear_large_u() {
        let s = TholeDamping::thole_screening_linear(100.0, 2.6);
        assert!((s - 1.0).abs() < 1e-10, "linear screening at large u: {s}");
    }
    #[test]
    fn test_charge_dipole_interaction_sign() {
        let q = 1.0;
        let mu = [1.0, 0.0, 0.0];
        let r = [2.0, 0.0, 0.0];
        let e = InducibleDipoleForce::charge_dipole_interaction(q, mu, r);
        assert!(
            e < 0.0,
            "charge-dipole should be negative when aligned: {e}"
        );
    }
    #[test]
    fn test_drucker_spring_zero_at_core() {
        let model = DruckerModel::new(1.5, 500.0);
        let r_core = [1.0, 2.0, 3.0];
        let r_drude = [1.0, 2.0, 3.0];
        let f = model.drude_spring_force(r_core, r_drude);
        for d in 0..3 {
            assert!(
                f[d].abs() < 1e-15,
                "spring force should be zero at core, f[{d}]={}",
                f[d]
            );
        }
    }
    #[test]
    fn test_drucker_dipole_linear_with_field() {
        let alpha = 2.0;
        let model = DruckerModel::new(alpha, 500.0);
        let e1 = [1.0, 0.0, 0.0];
        let e2 = [3.0, 0.0, 0.0];
        let p1 = model.dipole_moment(e1);
        let p2 = model.dipole_moment(e2);
        assert!(
            (p2[0] - 3.0 * p1[0]).abs() < 1e-12,
            "dipole should scale linearly with field"
        );
        assert!((p1[0] - alpha * e1[0]).abs() < 1e-12, "dipole = alpha * E");
    }
    #[test]
    fn test_scf_dipoles_converges() {
        let positions = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let alphas = vec![1.0, 1.0];
        let charges = vec![1.0, -1.0];
        let dipoles = self_consistent_dipoles(&positions, &alphas, &charges, 200, 1e-10);
        assert_eq!(dipoles.len(), 2, "should return one dipole per site");
        for dp in &dipoles {
            for d in 0..3 {
                assert!(dp[d].is_finite(), "dipole component should be finite");
            }
        }
    }
    #[test]
    fn test_polarization_energy_negative() {
        let dipoles = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let alphas = vec![2.0, 2.0];
        let u = polarization_energy(&dipoles, &alphas);
        assert!(
            u <= 0.0,
            "polarization energy should be non-positive, got {u}"
        );
    }
}
#[cfg(test)]
mod tests_polarizable_ext {
    use super::super::types::*;

    #[test]
    fn test_qeq_solver_builds_system() {
        let atoms = vec![QeqAtom::new(4.533, 6.919), QeqAtom::new(5.341, 10.126)];
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let solver = QeqSolver::new();
        let (a, b) = solver.build_system(&atoms, &positions, 0.0);
        assert_eq!(a.len(), 3, "matrix should be 3×3 for 2 atoms");
        assert!((b[2] - 0.0).abs() < 1e-15, "last rhs entry = Q_total = 0");
    }
    #[test]
    fn test_qeq_solve_charge_neutral() {
        let atoms = vec![QeqAtom::new(4.533, 6.919), QeqAtom::new(5.341, 10.126)];
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let solver = QeqSolver::new();
        let charges = solver
            .solve(&atoms, &positions, 0.0)
            .expect("solve should succeed");
        let sum: f64 = charges.iter().sum();
        assert!(sum.abs() < 1e-10, "total charge should be 0, got {sum}");
    }
    #[test]
    fn test_fluctuating_charge_kinetic_energy() {
        let n = 3;
        let mut fq = FluctuatingCharge::new(n, 1.0, vec![4.5, 5.3, 6.1], vec![7.0, 10.0, 8.5]);
        fq.charge_velocities = vec![1.0, 2.0, 3.0];
        let ke = fq.charge_kinetic_energy();
        assert!(
            (ke - 7.0).abs() < 1e-12,
            "kinetic energy should be 7.0, got {ke}"
        );
    }
    #[test]
    fn test_thole_damping_ext_exponential_zero() {
        let td = TholeDampingExt::new(TholeModel::Exponential, 1.3);
        let f = td.damping_fn(0.0);
        assert!(f.abs() < 1e-15, "damping at u=0 should be 0, got {f}");
    }
    #[test]
    fn test_thole_damping_ext_exponential_large_u() {
        let td = TholeDampingExt::new(TholeModel::Exponential, 1.3);
        let f = td.damping_fn(10.0);
        assert!(
            (f - 1.0).abs() < 1e-6,
            "damping at large u should approach 1"
        );
    }
    #[test]
    fn test_charmm_drude_params_cg331() {
        let p = CharmmDrudeParams::cg331();
        assert_eq!(p.atom_type, "CG331");
        assert!(p.polarizability > 0.0);
        assert!(p.k_drude > 0.0);
    }
    #[test]
    fn test_charmm_drude_param_table_lookup() {
        let table = CharmmDrudeParamTable::default_charmm();
        let p = table
            .get("OG311")
            .expect("OG311 should be in default table");
        assert_eq!(p.atom_type, "OG311");
    }
    #[test]
    fn test_swm4ndp_derived_polarizability() {
        let params = Swm4NdpParams::new();
        let alpha = params.derived_polarizability();
        assert!(alpha > 0.0, "derived polarizability should be positive");
    }
    #[test]
    fn test_amoeba_induction_energy_negative_for_aligned() {
        let _n = 2;
        let alpha = vec![1.0, 1.0];
        let mut pol = AmoebaPolarization::new(alpha, 0.39);
        pol.mu_ind = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let e_perm = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let u = pol.induction_energy(&e_perm);
        assert!(
            u < 0.0,
            "induction energy should be negative for aligned dipoles/fields"
        );
    }
    #[test]
    fn test_drude_ext_spring_forces_zero_at_coincidence() {
        let core_pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let drude_pos = core_pos.clone();
        let state = DrudeExtLagrangian::new(
            core_pos,
            drude_pos,
            vec![12.0, 16.0],
            vec![0.4, 0.4],
            1000.0,
            1.0,
            5.0,
        );
        let fs = state.spring_forces();
        for f in &fs {
            for d in 0..3 {
                assert!(
                    f[d].abs() < 1e-15,
                    "spring force should be zero at coincidence"
                );
            }
        }
    }
    #[test]
    fn test_drude_ext_spring_forces_nonzero_displaced() {
        let core_pos = vec![[0.0, 0.0, 0.0]];
        let drude_pos = vec![[0.1, 0.0, 0.0]];
        let state =
            DrudeExtLagrangian::new(core_pos, drude_pos, vec![12.0], vec![0.4], 1000.0, 1.0, 5.0);
        let fs = state.spring_forces();
        assert!(
            (fs[0][0] - (-100.0)).abs() < 1e-10,
            "spring force x should be -100"
        );
    }
    #[test]
    fn test_drude_self_energy_zero_displacement() {
        let e = Drude::compute_self_energy(0.0, 0.0, 0.0, 1000.0);
        assert!(
            (e - 0.0).abs() < 1e-15,
            "self energy should be 0 at coincidence, got {e}"
        );
    }
    #[test]
    fn test_drude_self_energy_positive() {
        let e = Drude::compute_self_energy(0.1, 0.0, 0.0, 1000.0);
        assert!(
            e > 0.0,
            "self energy must be positive for displaced Drude, got {e}"
        );
        let expected = 0.5 * 1000.0 * 0.01;
        assert!(
            (e - expected).abs() < 1e-10,
            "self energy mismatch: {e} vs {expected}"
        );
    }
    #[test]
    fn test_drude_self_energy_3d() {
        let e = Drude::compute_self_energy(0.1, 0.2, 0.3, 500.0);
        let expected = 0.5 * 500.0 * 0.14;
        assert!(
            (e - expected).abs() < 1e-10,
            "3D self energy mismatch: {e} vs {expected}"
        );
    }
    #[test]
    fn test_induced_dipole_model_energy_single_site() {
        let model = InducedDipoleModel::new(vec![2.0], vec![[1.0, 0.0, 0.0]]);
        let e = model.compute_polarization_energy();
        assert!(
            (e - (-0.25)).abs() < 1e-12,
            "polarization energy mismatch: {e}"
        );
    }
    #[test]
    fn test_induced_dipole_model_energy_zero_dipole() {
        let model = InducedDipoleModel::new(vec![1.5], vec![[0.0, 0.0, 0.0]]);
        let e = model.compute_polarization_energy();
        assert!(
            (e - 0.0).abs() < 1e-15,
            "zero dipole energy should be 0, got {e}"
        );
    }
    #[test]
    fn test_induced_dipole_model_energy_two_sites() {
        let model = InducedDipoleModel::new(vec![1.0, 2.0], vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let e = model.compute_polarization_energy();
        assert!((e - (-0.75)).abs() < 1e-12, "two-site energy mismatch: {e}");
    }
    #[test]
    fn test_thole_smeared_charge_large_r() {
        let q = 1.0;
        let alpha_i = 1.0;
        let alpha_j = 1.0;
        let a = 2.0;
        let r = 100.0;
        let q_smeared = TholeDamping::compute_smeared_charge(q, r, alpha_i, alpha_j, a);
        assert!(
            (q_smeared - q).abs() < 1e-6,
            "smeared charge should equal bare charge at large r: {q_smeared}"
        );
    }
    #[test]
    fn test_thole_smeared_charge_zero_r() {
        let q_smeared = TholeDamping::compute_smeared_charge(1.0, 0.0, 1.0, 1.0, 2.0);
        assert!(
            q_smeared.abs() < 1e-10,
            "smeared charge at r=0 should be 0, got {q_smeared}"
        );
    }
    #[test]
    fn test_thole_smeared_charge_sign_preserved() {
        let q_pos = TholeDamping::compute_smeared_charge(1.0, 3.0, 1.0, 1.0, 2.0);
        let q_neg = TholeDamping::compute_smeared_charge(-1.0, 3.0, 1.0, 1.0, 2.0);
        assert!(
            q_pos > 0.0,
            "positive charge should remain positive: {q_pos}"
        );
        assert!(
            q_neg < 0.0,
            "negative charge should remain negative: {q_neg}"
        );
        assert!(
            (q_pos + q_neg).abs() < 1e-14,
            "equal magnitude, opposite signs should cancel"
        );
    }
}
