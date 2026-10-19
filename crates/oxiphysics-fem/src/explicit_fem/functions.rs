//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Helper: linearly interpolate a (time, value) table.
pub(super) fn interpolate_table(table: &[(f64, f64)], t: f64) -> f64 {
    if table.is_empty() {
        return 0.0;
    }
    if t <= table[0].0 {
        return table[0].1;
    }
    if t >= table[table.len() - 1].0 {
        return table[table.len() - 1].1;
    }
    for i in 0..table.len() - 1 {
        if (table[i].0..=table[i + 1].0).contains(&t) {
            let dt = table[i + 1].0 - table[i].0;
            if dt < 1e-30 {
                return table[i].1;
            }
            let frac = (t - table[i].0) / dt;
            return table[i].1 + frac * (table[i + 1].1 - table[i].1);
        }
    }
    0.0
}
/// Compute the consistent mass matrix for a 4-node linear tetrahedron (12×12).
///
/// `M_e = rho * V / 20 * block_matrix`
pub fn tetrahedral_consistent_mass(nodes: &[[f64; 3]; 4], density: f64) -> [[f64; 12]; 12] {
    let v01 = [
        nodes[1][0] - nodes[0][0],
        nodes[1][1] - nodes[0][1],
        nodes[1][2] - nodes[0][2],
    ];
    let v02 = [
        nodes[2][0] - nodes[0][0],
        nodes[2][1] - nodes[0][1],
        nodes[2][2] - nodes[0][2],
    ];
    let v03 = [
        nodes[3][0] - nodes[0][0],
        nodes[3][1] - nodes[0][1],
        nodes[3][2] - nodes[0][2],
    ];
    let cross = [
        v01[1] * v02[2] - v01[2] * v02[1],
        v01[2] * v02[0] - v01[0] * v02[2],
        v01[0] * v02[1] - v01[1] * v02[0],
    ];
    let vol = (cross[0] * v03[0] + cross[1] * v03[1] + cross[2] * v03[2]).abs() / 6.0;
    let coeff = density * vol / 20.0;
    let mut m = [[0.0f64; 12]; 12];
    for i in 0..4 {
        for j in 0..4 {
            let val = if i == j { 2.0 } else { 1.0 } * coeff;
            for d in 0..3 {
                m[i * 3 + d][j * 3 + d] = val;
            }
        }
    }
    m
}
/// Row-sum lumped mass vector for a 4-node linear tetrahedron.
pub fn tetrahedral_lumped_mass_row_sum(nodes: &[[f64; 3]; 4], density: f64) -> [f64; 12] {
    let m = tetrahedral_consistent_mass(nodes, density);
    let mut lumped = [0.0f64; 12];
    for (i, lumped_i) in lumped.iter_mut().enumerate() {
        for &mij in &m[i] {
            *lumped_i += mij;
        }
    }
    lumped
}
/// HRZ lumped mass for a 4-node linear tetrahedron.
///
/// Scales each diagonal entry so that the sum of diagonal entries equals the
/// total element mass (sum of all entries of the consistent mass matrix / 3,
/// since there are 3 uncoupled direction blocks).
pub fn tetrahedral_lumped_mass_hrz(nodes: &[[f64; 3]; 4], density: f64) -> [f64; 12] {
    let m = tetrahedral_consistent_mass(nodes, density);
    let total_mass: f64 = m.iter().map(|row| row.iter().sum::<f64>()).sum::<f64>() / 3.0;
    let diag_sum: f64 = (0..12).map(|i| m[i][i]).sum();
    let scale = if diag_sum > 1e-30 {
        3.0 * total_mass / diag_sum
    } else {
        1.0
    };
    let mut lumped = [0.0f64; 12];
    for i in 0..12 {
        lumped[i] = m[i][i] * scale;
    }
    lumped
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::explicit_fem::*;
    #[test]
    fn test_explicit_solver_creation() {
        let config = ExplicitFemConfig::default();
        let solver = ExplicitFemSolver::new(12, config);
        assert_eq!(solver.state.n_dof(), 12);
        assert_eq!(solver.mass.len(), 12);
    }
    #[test]
    fn test_critical_time_step() {
        let config = ExplicitFemConfig {
            cfl_factor: 0.9,
            wave_speed: 5000.0,
            min_element_length: 0.01,
            ..Default::default()
        };
        let solver = ExplicitFemSolver::new(6, config);
        let dt = solver.critical_time_step();
        assert!((dt - 0.9 * 0.01 / 5000.0).abs() < 1e-15);
    }
    #[test]
    fn test_explicit_step_free_vibration() {
        let m: f64 = 1.0;
        let k = 100.0;
        let x0 = 0.01;
        let config = ExplicitFemConfig {
            cfl_factor: 0.5,
            wave_speed: (k / m).sqrt(),
            min_element_length: 0.01,
            check_energy: false,
            ..Default::default()
        };
        let mut solver = ExplicitFemSolver::new(1, config);
        solver.mass[0] = m;
        solver.state.displacement[0] = x0;
        let dt = solver.critical_time_step();
        solver.step(dt, |u| vec![k * u[0]]);
        assert!(solver.state.displacement[0] != x0 || solver.state.velocity[0].abs() > 0.0);
    }
    #[test]
    fn test_fix_dofs() {
        let config = ExplicitFemConfig::default();
        let mut solver = ExplicitFemSolver::new(6, config);
        solver.fix_dofs(&[0, 1, 2]);
        assert!(solver.fixed_dof[0]);
        assert!(solver.fixed_dof[1]);
        assert!(!solver.fixed_dof[3]);
    }
    #[test]
    fn test_kinetic_energy_zero_initial() {
        let config = ExplicitFemConfig::default();
        let solver = ExplicitFemSolver::new(6, config);
        assert_eq!(solver.kinetic_energy(), 0.0);
    }
    #[test]
    fn test_kinetic_energy_nonzero() {
        let config = ExplicitFemConfig::default();
        let mut solver = ExplicitFemSolver::new(3, config);
        solver.mass = vec![2.0, 2.0, 2.0];
        solver.state.velocity = vec![1.0, 0.0, 0.0];
        assert!((solver.kinetic_energy() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_assemble_lumped_mass() {
        let config = ExplicitFemConfig::default();
        let mut solver = ExplicitFemSolver::new(4, config);
        let elem = vec![(vec![0, 1, 2, 3], 4.0)];
        solver.assemble_lumped_mass(&elem);
        for i in 0..4 {
            assert!((solver.mass[i] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_estimate_stable_dt() {
        let config = ExplicitFemConfig::default();
        let solver = ExplicitFemSolver::new(3, config);
        let dt = solver.estimate_stable_dt_from_stiffness(1e6, 0.1);
        assert!(dt > 0.0);
    }
    #[test]
    fn test_set_initial_velocity() {
        let config = ExplicitFemConfig::default();
        let mut solver = ExplicitFemSolver::new(6, config);
        solver.set_initial_velocity(3, 5.0);
        assert_eq!(solver.state.velocity[3], 5.0);
    }
    #[test]
    fn test_set_external_force() {
        let config = ExplicitFemConfig::default();
        let mut solver = ExplicitFemSolver::new(6, config);
        solver.set_external_force(2, 100.0);
        assert_eq!(solver.external_force[2], 100.0);
    }
    #[test]
    fn test_run_until() {
        let config = ExplicitFemConfig {
            check_energy: false,
            ..Default::default()
        };
        let mut solver = ExplicitFemSolver::new(3, config);
        solver.mass = vec![1.0, 1.0, 1.0];
        solver.run_until(1e-6, |_u| vec![0.0, 0.0, 0.0]);
        assert!(solver.state.time >= 1e-6 - 1e-15);
    }
    #[test]
    fn test_hourglass_control_creation() {
        let hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 10);
        assert_eq!(hg.n_elements, 10);
        assert_eq!(hg.hourglass_energy.len(), 10);
    }
    #[test]
    fn test_hourglass_mode_amplitudes_zero_displacement() {
        let hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 1);
        let disp = [[0.0f64; 3]; 8];
        let amps = hg.compute_mode_amplitudes(&disp);
        for mode in &amps {
            for &v in mode {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn test_hourglass_mode_amplitudes_nonzero() {
        let hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 1);
        let mut disp = [[0.0f64; 3]; 8];
        for (n, &g) in HourglassControl::HOURGLASS_GAMMA[0].iter().enumerate() {
            disp[n][0] = g * 0.01;
        }
        let amps = hg.compute_mode_amplitudes(&disp);
        assert!(amps[0][0].abs() > 0.0);
    }
    #[test]
    fn test_stiffness_forces_zero_amplitude() {
        let hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 1);
        let q = [[0.0f64; 3]; 4];
        let forces = hg.stiffness_forces(&q, 200e9, 1e-6);
        for f in &forces {
            for &v in f {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn test_viscous_forces_zero() {
        let hg = HourglassControl::new(HourglassType::Viscous, 0.05, 0.1, 1);
        let q_dot = [[0.0f64; 3]; 4];
        let forces = hg.viscous_forces(&q_dot, 7800.0, 5000.0, 1e-6);
        for f in &forces {
            for &v in f {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn test_combined_forces() {
        let hg = HourglassControl::new(HourglassType::Combined, 0.05, 0.1, 1);
        let q = [[0.001f64; 3]; 4];
        let q_dot = [[0.01f64; 3]; 4];
        let forces = hg.combined_forces(&q, &q_dot, 200e9, 7800.0, 5000.0, 1e-6);
        let sum: f64 = forces.iter().flat_map(|f| f.iter()).map(|v| v.abs()).sum();
        assert!(sum > 0.0);
    }
    #[test]
    fn test_total_hourglass_energy() {
        let mut hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 3);
        hg.hourglass_energy = vec![1.0, 2.0, 3.0];
        assert!((hg.total_hourglass_energy() - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_hourglass_acceptable() {
        let mut hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 2);
        hg.hourglass_energy = vec![0.5, 0.5];
        assert!(hg.is_hourglass_acceptable(100.0, 0.05));
        assert!(!hg.is_hourglass_acceptable(100.0, 0.005));
    }
    #[test]
    fn test_update_hourglass_energy() {
        let mut hg = HourglassControl::new(HourglassType::Stiffness, 0.05, 0.1, 1);
        let q = [[1.0f64; 3]; 4];
        hg.update_hourglass_energy(0, &q, 1000.0);
        assert!((hg.hourglass_energy[0] - 6000.0).abs() < 1e-9);
    }
    #[test]
    fn test_contact_detection_creation() {
        let cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        assert_eq!(cd.pairs.len(), 0);
    }
    #[test]
    fn test_node_to_segment_gap_2d_no_contact() {
        let slave = [0.0, 2.0];
        let ma = [-1.0, 0.0];
        let mb = [1.0, 0.0];
        let (gap, t, _normal) = ContactDetectionFem::node_to_segment_gap_2d(slave, ma, mb);
        assert!(gap > 0.0, "gap should be positive (no contact): {}", gap);
        assert!((t - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_node_to_segment_gap_2d_contact() {
        let slave = [0.0, -0.01];
        let ma = [-1.0, 0.0];
        let mb = [1.0, 0.0];
        let (gap, _t, _normal) = ContactDetectionFem::node_to_segment_gap_2d(slave, ma, mb);
        assert!(gap < 0.0, "should be penetration: {}", gap);
    }
    #[test]
    fn test_penalty_force_no_contact() {
        let cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        assert_eq!(cd.penalty_force(0.01), 0.0);
    }
    #[test]
    fn test_penalty_force_contact() {
        let cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        let f = cd.penalty_force(-0.001);
        assert!((f - 1000.0).abs() < 1e-9);
    }
    #[test]
    fn test_augmented_lagrangian_force() {
        let cd = ContactDetectionFem::new(ContactEnforcement::AugmentedLagrangian, 1e4, 0.3);
        let f = cd.augmented_lagrangian_force(-0.001, 0.0);
        assert!(f > 0.0);
    }
    #[test]
    fn test_detect_no_pairs_when_far() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        cd.search_radius = 1e-4;
        let positions = vec![0.0, 10.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0];
        cd.detect_node_to_segment(&positions, &[[1, 2]], &[[0.0, 1.0, 0.0]], &[0])
            .unwrap();
        assert_eq!(cd.pairs.len(), 0);
    }
    #[test]
    fn test_contact_pair_creation() {
        let pair = ContactPair::new(5, [3, 4]);
        assert_eq!(pair.slave_node, 5);
        assert_eq!(pair.master_segment, [3, 4]);
        assert_eq!(pair.gap, 0.0);
    }
    #[test]
    fn test_assemble_contact_forces_empty() {
        let cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        let mut force = vec![0.0f64; 9];
        cd.assemble_contact_forces(&mut force);
        assert!(force.iter().all(|&f| f == 0.0));
    }
    #[test]
    fn test_update_lagrange_multipliers() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Lagrange, 1e6, 0.3);
        let mut pair = ContactPair::new(0, [1, 2]);
        pair.gap = -0.001;
        pair.lambda = 0.0;
        cd.pairs.push(pair);
        cd.update_lagrange_multipliers();
        assert!(cd.pairs[0].lambda >= 0.0);
    }
    /// Slave node `0` and master segment `1-2` along x (outward normal +z).
    fn wall_positions(slave: [f64; 3]) -> Vec<f64> {
        vec![slave[0], slave[1], slave[2], -1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    }
    #[test]
    fn test_node_to_segment_gap_3d_is_signed() {
        let a = [-1.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let n = [0.0, 0.0, 1.0];
        // 1e-3 below the segment (inside the master body): penetration.
        let (gap, t, normal) =
            ContactDetectionFem::node_to_segment_gap_3d([0.0, 0.0, -1e-3], a, b, n).unwrap();
        assert!((gap + 1e-3).abs() < 1e-15, "gap = {gap}");
        assert!((t - 0.5).abs() < 1e-15);
        assert_eq!(normal, [0.0, 0.0, 1.0]);
        // 2e-3 above: separation.
        let (gap, t, _) =
            ContactDetectionFem::node_to_segment_gap_3d([0.5, 0.0, 2e-3], a, b, n).unwrap();
        assert!((gap - 2e-3).abs() < 1e-15, "gap = {gap}");
        assert!((t - 0.75).abs() < 1e-15);
        // The normal is normalised, and its sign decides the side.
        let (gap, _, normal) =
            ContactDetectionFem::node_to_segment_gap_3d([0.0, 0.0, 2e-3], a, b, [0.0, 0.0, -4.0])
                .unwrap();
        assert!((gap + 2e-3).abs() < 1e-15, "gap = {gap}");
        assert_eq!(normal, [0.0, 0.0, -1.0]);
    }
    #[test]
    fn test_node_to_segment_gap_3d_rejects_degenerate_normal() {
        let a = [-1.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        for bad in [
            [0.0, 0.0, 0.0],
            [f64::NAN, 0.0, 1.0],
            [f64::INFINITY, 0.0, 0.0],
        ] {
            assert!(matches!(
                ContactDetectionFem::node_to_segment_gap_3d([0.0, 0.0, 0.0], a, b, bad),
                Err(crate::Error::InvalidInput(_))
            ));
        }
    }
    #[test]
    fn test_penetrating_slave_gets_penalty_force_and_zero_net_force() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        cd.search_radius = 1e-2;
        let positions = wall_positions([0.25, 0.0, -1e-3]);
        cd.detect_node_to_segment(&positions, &[[1, 2]], &[[0.0, 0.0, 1.0]], &[0])
            .unwrap();
        assert_eq!(cd.pairs.len(), 1);
        let pair = &cd.pairs[0];
        assert!((pair.gap + 1e-3).abs() < 1e-15, "gap = {}", pair.gap);
        assert!((pair.t - 0.625).abs() < 1e-15, "t = {}", pair.t);
        assert_eq!(pair.normal, [0.0, 0.0, 1.0]);
        let mut force = vec![0.0_f64; 9];
        cd.assemble_contact_forces(&mut force);
        // k_pen * 1e-3 = 1000 N pushes the slave out along +z ...
        assert!((force[2] - 1000.0).abs() < 1e-9, "slave fz = {}", force[2]);
        // ... and the reaction is split (1 - t) : t over the master nodes.
        assert!(
            (force[5] + 375.0).abs() < 1e-9,
            "master 1 fz = {}",
            force[5]
        );
        assert!(
            (force[8] + 625.0).abs() < 1e-9,
            "master 2 fz = {}",
            force[8]
        );
        // Force and moment (about y) balance.
        let sum: f64 = (0..3).map(|i| force[i * 3 + 2]).sum();
        assert!(sum.abs() < 1e-9, "net force {sum}");
        let xs = [0.25, -1.0, 1.0];
        let moment: f64 = (0..3).map(|i| xs[i] * force[i * 3 + 2]).sum();
        assert!(moment.abs() < 1e-9, "net moment {moment}");
        // No spurious lateral components.
        for i in 0..3 {
            assert_eq!(force[i * 3], 0.0);
            assert_eq!(force[i * 3 + 1], 0.0);
        }
    }
    #[test]
    fn test_separated_slave_has_positive_gap_and_no_force() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        cd.search_radius = 1e-2;
        let positions = wall_positions([0.0, 0.0, 5e-3]);
        cd.detect_node_to_segment(&positions, &[[1, 2]], &[[0.0, 0.0, 1.0]], &[0])
            .unwrap();
        assert_eq!(cd.pairs.len(), 1);
        assert!((cd.pairs[0].gap - 5e-3).abs() < 1e-15);
        let mut force = vec![0.0_f64; 9];
        cd.assemble_contact_forces(&mut force);
        assert!(force.iter().all(|&f| f == 0.0));
    }
    #[test]
    fn test_detection_ignores_nodes_beyond_the_segment_ends_and_deep_tunnelling() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        cd.search_radius = 1.0;
        let normals = [[0.0, 0.0, 1.0]];
        // Projects beyond the segment end (t = 1.25) although within the radius.
        let beyond = wall_positions([1.5, 0.0, -1e-4]);
        cd.detect_node_to_segment(&beyond, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert!(cd.pairs.is_empty(), "node beyond the end was paired");
        // Behind the segment but farther than the capture radius.
        cd.search_radius = 1e-2;
        let deep = wall_positions([0.0, 0.0, -0.5]);
        cd.detect_node_to_segment(&deep, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert!(cd.pairs.is_empty(), "tunnelled node was paired");
        // A node that is an end node of the segment never contacts it.
        let own = wall_positions([0.0, 0.0, 0.0]);
        cd.detect_node_to_segment(&own, &[[0, 2]], &normals, &[0])
            .unwrap();
        assert!(cd.pairs.is_empty());
    }
    #[test]
    fn test_detect_node_to_segment_validates_input() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Penalty, 1e6, 0.3);
        let positions = wall_positions([0.0, 0.0, -1e-4]);
        let mut seeded = ContactPair::new(0, [1, 2]);
        seeded.lambda = 3.0;
        cd.pairs.push(seeded);
        let unit_z = [[0.0, 0.0, 1.0]];
        assert!(matches!(
            cd.detect_node_to_segment(&positions, &[[1, 2]], &[], &[0]),
            Err(crate::Error::InvalidInput(_))
        ));
        assert!(matches!(
            cd.detect_node_to_segment(&positions, &[[1, 2]], &[[0.0, 0.0, 0.0]], &[0]),
            Err(crate::Error::InvalidInput(_))
        ));
        assert!(matches!(
            cd.detect_node_to_segment(&positions, &[[1, 7]], &unit_z, &[0]),
            Err(crate::Error::InvalidInput(_))
        ));
        assert!(matches!(
            cd.detect_node_to_segment(&positions, &[[1, 2]], &unit_z, &[9]),
            Err(crate::Error::InvalidInput(_))
        ));
        // A rejected request leaves the existing pairs alone.
        assert_eq!(cd.pairs.len(), 1);
        assert_eq!(cd.pairs[0].lambda, 3.0);
    }
    #[test]
    fn test_detection_carries_multipliers_over() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::AugmentedLagrangian, 1e6, 0.3);
        cd.search_radius = 1e-2;
        let normals = [[0.0, 0.0, 1.0]];
        let positions = wall_positions([0.0, 0.0, -1e-3]);
        cd.detect_node_to_segment(&positions, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert_eq!(cd.pairs[0].lambda, 0.0);
        cd.pairs[0].lambda = 7.0;
        // Same pair detected again: the multiplier survives.
        cd.detect_node_to_segment(&positions, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert_eq!(cd.pairs[0].lambda, 7.0);
        // A different segment starts from zero.
        cd.detect_node_to_segment(
            &positions,
            &[[1, 2], [2, 1]],
            &[normals[0], [0.0, 0.0, -1.0]],
            &[0],
        )
        .unwrap();
        let by_segment = |s: [usize; 2]| cd.pairs.iter().find(|p| p.master_segment == s);
        assert_eq!(by_segment([1, 2]).map(|p| p.lambda), Some(7.0));
        assert_eq!(by_segment([2, 1]).map(|p| p.lambda), Some(0.0));
        // The pair leaves the candidate set: its multiplier is dropped.
        let far = wall_positions([0.0, 0.0, 1.0]);
        cd.detect_node_to_segment(&far, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert!(cd.pairs.is_empty());
        cd.detect_node_to_segment(&positions, &[[1, 2]], &normals, &[0])
            .unwrap();
        assert_eq!(cd.pairs[0].lambda, 0.0);
    }
    #[test]
    fn test_augmented_lagrangian_force_is_compressive_positive() {
        let cd = ContactDetectionFem::new(ContactEnforcement::AugmentedLagrangian, 1e4, 0.3);
        // eps = 1e4, gap = -1e-3: the penalty term alone gives 10.
        assert!((cd.augmented_lagrangian_force(-1e-3, 0.0) - 10.0).abs() < 1e-12);
        // A pushing multiplier adds to it.
        assert!((cd.augmented_lagrangian_force(-1e-3, 5.0) - 15.0).abs() < 1e-12);
        // Separation relaxes the multiplier; the force never becomes tensile.
        assert!((cd.augmented_lagrangian_force(1e-3, 20.0) - 10.0).abs() < 1e-12);
        assert_eq!(cd.augmented_lagrangian_force(1e-3, 5.0), 0.0);
        assert_eq!(cd.augmented_lagrangian_force(1.0, 0.0), 0.0);
    }
    #[test]
    fn test_multiplier_update_grows_under_penetration_and_relaxes_when_open() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Lagrange, 1e6, 0.3);
        cd.augmentation = 1e4;
        cd.pairs.push(ContactPair::new(0, [1, 2]));
        cd.pairs[0].gap = -1e-3;
        cd.update_lagrange_multipliers();
        assert!((cd.pairs[0].lambda - 10.0).abs() < 1e-12);
        cd.update_lagrange_multipliers();
        assert!((cd.pairs[0].lambda - 20.0).abs() < 1e-12);
        // The pair opens: lambda decays and never goes negative.
        cd.pairs[0].gap = 1e-3;
        cd.update_lagrange_multipliers();
        assert!((cd.pairs[0].lambda - 10.0).abs() < 1e-12);
        cd.update_lagrange_multipliers();
        assert_eq!(cd.pairs[0].lambda, 0.0);
        cd.update_lagrange_multipliers();
        assert_eq!(cd.pairs[0].lambda, 0.0);
    }
    #[test]
    fn test_lagrange_enforcement_force_is_the_multiplier() {
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Lagrange, 1e6, 0.3);
        let mut pair = ContactPair::new(0, [1, 2]);
        pair.normal = [0.0, 0.0, 1.0];
        pair.gap = -1e-3;
        pair.lambda = 7.0;
        cd.pairs.push(pair);
        let mut force = vec![0.0_f64; 9];
        cd.assemble_contact_forces(&mut force);
        assert!((force[2] - 7.0).abs() < 1e-12, "slave fz = {}", force[2]);
        assert!((force[5] + 3.5).abs() < 1e-12 && (force[8] + 3.5).abs() < 1e-12);
        // Zero multiplier: no force even if the pair is penetrating.
        cd.pairs[0].lambda = 0.0;
        let mut force = vec![0.0_f64; 9];
        cd.assemble_contact_forces(&mut force);
        assert!(force.iter().all(|&f| f == 0.0));
    }
    /// One mass on a spring (stiffness `k`) pressed by `load` into a rigid wall
    /// at zero clearance.  The contact is solved with the augmented Lagrangian
    /// (Uzawa) iteration implemented by `ContactDetectionFem`.
    #[test]
    fn test_augmented_lagrangian_drives_gap_to_zero_on_rigid_wall() {
        let (k, load) = (1.0e3, 50.0);
        let mut cd = ContactDetectionFem::new(ContactEnforcement::AugmentedLagrangian, 1e6, 0.0);
        cd.augmentation = 4.0 * k;
        let mut pair = ContactPair::new(0, [1, 2]);
        pair.normal = [0.0, 0.0, 1.0];
        cd.pairs.push(pair);
        let eps = cd.augmentation;
        let mut gaps = Vec::new();
        for _ in 0..=20 {
            // Equilibrium load = k u + max(0, lambda + eps u) for the current lambda.
            let u = ((load - cd.pairs[0].lambda) / (k + eps)).max(0.0);
            cd.pairs[0].gap = -u;
            let contact = cd.augmented_lagrangian_force(cd.pairs[0].gap, cd.pairs[0].lambda);
            assert!(
                (k * u + contact - load).abs() < 1e-9 * load,
                "equilibrium violated: {}",
                k * u + contact - load
            );
            gaps.push(cd.pairs[0].gap);
            cd.update_lagrange_multipliers();
        }
        let initial = gaps[0].abs();
        assert!(initial > 0.0);
        // The error contracts by k / (k + eps) = 0.2 per update (checked while
        // the gap is still far above the rounding level of `load - lambda`).
        for w in gaps.windows(2).take(8) {
            assert!(
                (w[1] / w[0] - 0.2).abs() < 1e-8,
                "contraction {}",
                w[1] / w[0]
            );
        }
        assert!(gaps[20].abs() < 1e-8 * initial, "gap {}", gaps[20]);
        assert!(
            (cd.pairs[0].lambda - load).abs() < 1e-8 * load,
            "lambda = {}",
            cd.pairs[0].lambda
        );
    }
    /// Same problem with pure Lagrange enforcement: the contact force is the
    /// multiplier itself, which the Uzawa update drives to the applied load.
    #[test]
    fn test_lagrange_enforcement_converges_to_applied_load() {
        let (k, load) = (1.0e3, 50.0);
        let mut cd = ContactDetectionFem::new(ContactEnforcement::Lagrange, 1e6, 0.0);
        cd.augmentation = 0.5 * k;
        let mut pair = ContactPair::new(0, [1, 2]);
        pair.normal = [0.0, 0.0, 1.0];
        cd.pairs.push(pair);
        for _ in 0..40 {
            // Mass equilibrium: k u + lambda = load, penetration u = -gap.
            let u = (load - cd.pairs[0].lambda) / k;
            cd.pairs[0].gap = -u;
            cd.update_lagrange_multipliers();
        }
        assert!(
            (cd.pairs[0].lambda - load).abs() < 1e-9 * load,
            "lambda = {}",
            cd.pairs[0].lambda
        );
        assert!(cd.pairs[0].gap.abs() < 1e-12);
        let mut force = vec![0.0_f64; 9];
        cd.assemble_contact_forces(&mut force);
        assert!(
            (force[2] - load).abs() < 1e-8 * load,
            "slave fz = {}",
            force[2]
        );
    }
    #[test]
    fn test_erosion_creation() {
        let ea = ErosionAlgorithm::new(FailureCriterion::PlasticStrain(1.0), 10, true);
        assert_eq!(ea.states.len(), 10);
        assert_eq!(ea.n_eroded, 0);
    }
    #[test]
    fn test_should_erode_plastic_strain() {
        let state = ElementErosionState::new(1.0);
        let mut s = state;
        s.plastic_strain = 1.5;
        assert!(s.should_erode(&FailureCriterion::PlasticStrain(1.0)));
        s.plastic_strain = 0.5;
        assert!(!s.should_erode(&FailureCriterion::PlasticStrain(1.0)));
    }
    #[test]
    fn test_should_erode_damage() {
        let mut state = ElementErosionState::new(1.0);
        state.damage = 0.95;
        assert!(state.should_erode(&FailureCriterion::Damage(0.9)));
    }
    #[test]
    fn test_should_erode_principal_stress() {
        let mut state = ElementErosionState::new(1.0);
        state.max_principal_stress = 600e6;
        assert!(state.should_erode(&FailureCriterion::PrincipalStress(500e6)));
    }
    #[test]
    fn test_should_erode_volume_fraction() {
        let mut state = ElementErosionState::new(1.0);
        state.volume_ratio = 0.05;
        assert!(state.should_erode(&FailureCriterion::VolumeFraction(0.1)));
    }
    #[test]
    fn test_update_element_erodes() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::PlasticStrain(0.5), 5, true);
        ea.set_element_mass(2, 10.0);
        let eroded = ea.update_element(2, 0.0, 0.8, 0.0, 1.0, 1);
        assert!(eroded);
        assert!(ea.states[2].eroded);
        assert_eq!(ea.n_eroded, 1);
        assert!((ea.total_eroded_mass - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_update_element_no_erosion() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::PlasticStrain(1.0), 5, true);
        let eroded = ea.update_element(0, 0.0, 0.5, 0.0, 1.0, 1);
        assert!(!eroded);
    }
    #[test]
    fn test_eroded_elements_list() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::Damage(0.9), 4, true);
        ea.update_element(1, 0.95, 0.0, 0.0, 1.0, 1);
        ea.update_element(3, 0.95, 0.0, 0.0, 1.0, 2);
        let eroded = ea.eroded_elements();
        assert_eq!(eroded, vec![1, 3]);
    }
    #[test]
    fn test_active_elements_list() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::Damage(0.9), 4, true);
        ea.update_element(0, 0.95, 0.0, 0.0, 1.0, 1);
        let active = ea.active_elements();
        assert_eq!(active, vec![1, 2, 3]);
    }
    #[test]
    fn test_erosion_fraction() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::PlasticStrain(0.5), 4, true);
        ea.update_element(0, 0.0, 1.0, 0.0, 1.0, 1);
        ea.update_element(2, 0.0, 1.0, 0.0, 1.0, 2);
        assert!((ea.erosion_fraction() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_redistribute_mass() {
        let mut ea = ErosionAlgorithm::new(FailureCriterion::PlasticStrain(0.5), 2, true);
        ea.set_element_mass(0, 4.0);
        ea.update_element(0, 0.0, 1.0, 0.0, 1.0, 1);
        let mut node_masses = vec![0.0f64; 4];
        ea.redistribute_mass_to_nodes(0, &mut node_masses, &[0, 1, 2, 3]);
        let total: f64 = node_masses.iter().sum();
        assert!((total - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_p_wave_speed() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let cp = wp.p_wave_speed();
        assert!(cp > 0.0);
        assert!(cp > 5000.0 && cp < 7000.0);
    }
    #[test]
    fn test_s_wave_speed() {
        let wp = WavePropagationFem::new(WaveType::Transverse, 200e9, 7800.0, 0.3, 0.01);
        let cs = wp.s_wave_speed();
        assert!(cs > 3000.0 && cs < 4000.0);
    }
    #[test]
    fn test_bar_wave_speed() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let cb = wp.bar_wave_speed();
        assert!(cb > 4000.0 && cb < 6000.0);
    }
    #[test]
    fn test_rayleigh_wave_speed() {
        let wp = WavePropagationFem::new(WaveType::Rayleigh, 200e9, 7800.0, 0.3, 0.01);
        let cr = wp.rayleigh_wave_speed();
        let cs = wp.s_wave_speed();
        assert!(cr < cs);
        assert!(cr > 0.0);
    }
    #[test]
    fn test_dispersion_curve_longitudinal() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let freqs: Vec<f64> = (1..=5).map(|i| i as f64 * 1000.0).collect();
        let result = wp.dispersion_curve(&freqs);
        assert_eq!(result.frequencies.len(), 5);
        assert!(result.phase_velocities.iter().all(|&c| c > 0.0));
    }
    #[test]
    fn test_dispersion_curve_flexural() {
        let wp = WavePropagationFem::new(WaveType::Flexural, 200e9, 7800.0, 0.3, 0.01);
        let freqs = vec![100.0, 1000.0, 10000.0];
        let result = wp.dispersion_curve(&freqs);
        for (cp, cg) in result
            .phase_velocities
            .iter()
            .zip(result.group_velocities.iter())
        {
            assert!((cg / cp - 2.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_reflection_transmission_same_material() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let c = wp.bar_wave_speed();
        let (r, t) = wp.reflection_transmission(7800.0, c);
        assert!(r.abs() < 1e-12, "R should be 0 for same material: {}", r);
        assert!(
            (t - 1.0).abs() < 1e-12,
            "T should be 1 for same material: {}",
            t
        );
    }
    #[test]
    fn test_arrival_time() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let t = wp.arrival_time(5.064);
        assert!(t > 0.0);
    }
    #[test]
    fn test_fem_dispersion() {
        let ratio = WavePropagationFem::fem_dispersion(0.001, 5000.0, 100.0);
        assert!(ratio > 0.0);
        assert!((ratio - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_hopkinson_pulse() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        assert_eq!(wp.hopkinson_pulse(0.0, 1e6, 1e-3), 1e6);
        assert_eq!(wp.hopkinson_pulse(2e-3, 1e6, 1e-3), 0.0);
    }
    #[test]
    fn test_blast_friedlander_at_zero() {
        let params = BlastParameters {
            peak_pressure: 1e6,
            duration: 5e-3,
            decay_coeff: 1.0,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::Friedlander, params);
        let p = blast.pressure_at(0.0);
        assert!((p - 1e6).abs() < 1e-6);
    }
    #[test]
    fn test_blast_friedlander_zero_after_duration() {
        let params = BlastParameters {
            peak_pressure: 1e6,
            duration: 5e-3,
            decay_coeff: 1.0,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::Friedlander, params);
        assert_eq!(blast.pressure_at(10e-3), 0.0);
    }
    #[test]
    fn test_blast_triangular() {
        let params = BlastParameters {
            peak_pressure: 2e6,
            duration: 10e-3,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::Triangular, params);
        assert!((blast.pressure_at(0.0) - 2e6).abs() < 1e-6);
        assert!((blast.pressure_at(5e-3) - 1e6).abs() < 1e-6);
        assert_eq!(blast.pressure_at(15e-3), 0.0);
    }
    #[test]
    fn test_blast_conwep_positive_pressure() {
        let params = BlastParameters {
            standoff: 2.0,
            charge_mass: 1.0,
            peak_pressure: 1e5,
            duration: 5e-3,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::ConWep, params);
        // The shock reaches 2 m from 1 kg of TNT after about 1.7 ms: nothing
        // earlier, a decaying pressure afterwards.
        assert_eq!(blast.pressure_at(1e-3), 0.0);
        let wave = BlastWave::surface_burst(2.0, 1.0).unwrap();
        let p = blast.pressure_at(wave.arrival_time + 0.1 * wave.positive_duration);
        assert!(p > 0.0 && p < wave.peak_overpressure, "p = {p}");
    }
    #[test]
    fn test_blast_conwep_follows_kingery_bulmash_pulse() {
        // 1 t of TNT at 50 m: UN IATG 01.80 Table 5 gives 0.43 bar and 82.4 ms.
        let params = BlastParameters {
            standoff: 50.0,
            charge_mass: 1000.0,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::ConWep, params);
        let wave = BlastWave::surface_burst(50.0, 1000.0).unwrap();
        assert!((wave.peak_overpressure - 0.43e5).abs() < 0.015 * 0.43e5);
        assert!((wave.arrival_time - 82.4e-3).abs() < 0.015 * 82.4e-3);
        assert_eq!(blast.pressure_at(0.5 * wave.arrival_time), 0.0);
        let at_front = blast.pressure_at(wave.arrival_time);
        assert!((at_front - wave.peak_overpressure).abs() < 1e-9 * wave.peak_overpressure);
        let mid = wave.arrival_time + 0.5 * wave.positive_duration;
        let p_mid = blast.pressure_at(mid);
        assert!(
            p_mid > 0.0 && p_mid < 0.5 * wave.peak_overpressure,
            "p_mid = {p_mid}"
        );
        let after = wave.arrival_time + 1.01 * wave.positive_duration;
        assert_eq!(blast.pressure_at(after), 0.0);
    }
    #[test]
    fn test_blast_conwep_ignores_the_friedlander_inputs() {
        let base = BlastParameters {
            standoff: 8.0,
            charge_mass: 25.0,
            ..Default::default()
        };
        let other = BlastParameters {
            peak_pressure: 7.0e9,
            duration: 1.0,
            decay_coeff: 42.0,
            ..base.clone()
        };
        let a = BlastLoadingFem::new(BlastModel::ConWep, base);
        let b = BlastLoadingFem::new(BlastModel::ConWep, other);
        let wave = BlastWave::surface_burst(8.0, 25.0).unwrap();
        for k in 0..=12 {
            let t = wave.arrival_time + wave.positive_duration * f64::from(k) / 10.0;
            assert_eq!(a.pressure_at(t), b.pressure_at(t));
        }
    }
    #[test]
    fn test_blast_conwep_without_charge_applies_no_load() {
        for mass in [0.0, -5.0, f64::NAN, f64::INFINITY] {
            let params = BlastParameters {
                standoff: 3.0,
                charge_mass: mass,
                ..Default::default()
            };
            let blast = BlastLoadingFem::new(BlastModel::ConWep, params);
            for t in [0.0, 1e-4, 1e-3, 1e-2, 1.0] {
                assert_eq!(blast.pressure_at(t), 0.0, "mass {mass}, t = {t}");
            }
        }
    }
    #[test]
    fn test_blast_conwep_peak_decays_to_zero_with_standoff() {
        // The previous fit tended to 0.84 MPa (about 8 atm) at infinite standoff.
        let mut previous = f64::INFINITY;
        for standoff in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1.0e3, 1.0e5] {
            let wave = BlastWave::surface_burst(standoff, 1.0).unwrap();
            let params = BlastParameters {
                standoff,
                charge_mass: 1.0,
                ..Default::default()
            };
            let blast = BlastLoadingFem::new(BlastModel::ConWep, params);
            let peak = blast.pressure_at(wave.arrival_time);
            assert!(peak < previous, "R = {standoff}: {peak} >= {previous}");
            previous = peak;
        }
        // 40 m from 1 kg of TNT: below 5 kPa, and essentially nothing at 100 km.
        let params = BlastParameters {
            standoff: 40.0,
            charge_mass: 1.0,
            ..Default::default()
        };
        let near = BlastLoadingFem::new(BlastModel::ConWep, params);
        let wave = BlastWave::surface_burst(40.0, 1.0).unwrap();
        assert!(near.pressure_at(wave.arrival_time) < 5.0e3);
        assert!(previous < 10.0, "peak at 100 km: {previous} Pa");
    }
    #[test]
    fn test_blast_conwep_nodal_forces_follow_the_pulse() {
        let params = BlastParameters {
            standoff: 6.0,
            charge_mass: 50.0,
            ..Default::default()
        };
        let mut blast = BlastLoadingFem::new(BlastModel::ConWep, params);
        blast.surface_nodes = vec![0];
        blast.node_normals = vec![[0.0, 0.0, -1.0]];
        blast.tributary_areas = vec![0.25];
        let wave = BlastWave::surface_burst(6.0, 50.0).unwrap();
        let mut early = vec![0.0f64; 3];
        blast.assemble_blast_forces(0.5 * wave.arrival_time, &mut early);
        assert!(early.iter().all(|&f| f == 0.0));
        let mut front = vec![0.0f64; 3];
        blast.assemble_blast_forces(wave.arrival_time, &mut front);
        let expected = -wave.peak_overpressure * 0.25;
        assert!(
            (front[2] - expected).abs() < 1e-9 * expected.abs(),
            "fz = {}",
            front[2]
        );
        assert_eq!((front[0], front[1]), (0.0, 0.0));
    }
    #[test]
    fn test_blast_user_defined_interpolation() {
        let table = vec![(0.0, 0.0), (1.0, 100.0), (2.0, 0.0)];
        let params = BlastParameters {
            table,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::UserDefined, params);
        assert!((blast.pressure_at(0.5) - 50.0).abs() < 1e-10);
        assert!((blast.pressure_at(1.5) - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_assemble_blast_forces() {
        let params = BlastParameters {
            peak_pressure: 1e6,
            duration: 1.0,
            ..Default::default()
        };
        let mut blast = BlastLoadingFem::new(BlastModel::Triangular, params);
        blast.surface_nodes = vec![0];
        blast.node_normals = vec![[0.0, 1.0, 0.0]];
        blast.tributary_areas = vec![1.0];
        let mut force = vec![0.0f64; 3];
        blast.assemble_blast_forces(0.0, &mut force);
        assert!((force[1] - 1e6).abs() < 1.0);
    }
    #[test]
    fn test_compute_impulse() {
        let params = BlastParameters::default();
        let mut blast = BlastLoadingFem::new(BlastModel::Triangular, params);
        blast.pressure_history = vec![(0.0, 100.0), (1.0, 50.0), (2.0, 0.0)];
        let impulse = blast.compute_impulse();
        assert!((impulse - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_dynamic_load_factor() {
        let params = BlastParameters::default();
        let mut blast = BlastLoadingFem::new(BlastModel::Triangular, params);
        blast.response_history = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.01)];
        let dlf = blast.dynamic_load_factor(0.01);
        assert!((dlf - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_reflected_pressure() {
        let params = BlastParameters {
            peak_pressure: 1e5,
            duration: 1.0,
            ..Default::default()
        };
        let blast = BlastLoadingFem::new(BlastModel::Triangular, params);
        let pr = blast.reflected_pressure(0.0, 101325.0);
        assert!(pr > blast.pressure_at(0.0));
    }
    #[test]
    fn test_tetrahedral_consistent_mass_positive() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let m = tetrahedral_consistent_mass(&nodes, 1000.0);
        for (i, row) in m.iter().enumerate() {
            assert!(row[i] > 0.0);
        }
    }
    #[test]
    fn test_tetrahedral_consistent_mass_symmetric() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let m = tetrahedral_consistent_mass(&nodes, 7800.0);
        for (i, row) in m.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - m[j][i]).abs() < 1e-15);
            }
        }
    }
    #[test]
    fn test_lumped_mass_row_sum_positive() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let lumped = tetrahedral_lumped_mass_row_sum(&nodes, 1000.0);
        for &v in &lumped {
            assert!(v > 0.0);
        }
    }
    #[test]
    fn test_lumped_mass_hrz_preserves_total_mass() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let density = 1000.0;
        let hrz = tetrahedral_lumped_mass_hrz(&nodes, density);
        let m = tetrahedral_consistent_mass(&nodes, density);
        let total_consistent: f64 = m.iter().map(|row| row.iter().sum::<f64>()).sum::<f64>() / 3.0;
        let total_hrz: f64 = hrz.iter().sum::<f64>() / 3.0;
        assert!(
            (total_hrz - total_consistent).abs() / total_consistent < 1e-8,
            "HRZ mass={}, consistent mass={}",
            total_hrz,
            total_consistent
        );
        for &v in &hrz {
            assert!(v > 0.0);
        }
    }
    #[test]
    fn test_lumped_mass_total_conserved() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        let density = 7800.0;
        let m = tetrahedral_consistent_mass(&nodes, density);
        let total_consistent: f64 = m.iter().map(|row| row.iter().sum::<f64>()).sum::<f64>() / 3.0;
        let lumped = tetrahedral_lumped_mass_row_sum(&nodes, density);
        let total_lumped: f64 = lumped.iter().sum::<f64>() / 3.0;
        assert!((total_consistent - total_lumped).abs() / total_consistent < 1e-10);
    }
    #[test]
    fn test_interpolate_table_boundaries() {
        let table = vec![(0.0, 10.0), (1.0, 20.0)];
        assert_eq!(interpolate_table(&table, -1.0), 10.0);
        assert_eq!(interpolate_table(&table, 2.0), 20.0);
        assert!((interpolate_table(&table, 0.5) - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_table_empty() {
        let table: Vec<(f64, f64)> = Vec::new();
        assert_eq!(interpolate_table(&table, 0.5), 0.0);
    }
    #[test]
    fn test_wave_reflection_transmission_conservation() {
        let wp = WavePropagationFem::new(WaveType::Longitudinal, 200e9, 7800.0, 0.3, 0.01);
        let z1 = 7800.0 * wp.bar_wave_speed();
        let density2 = 2700.0;
        let c2 = (70e9 / 2700.0_f64).sqrt();
        let z2 = density2 * c2;
        let (r, t) = wp.reflection_transmission(density2, c2);
        let power_check = r * r + (z1 / z2) * t * t;
        assert!((power_check - 1.0).abs() < 1e-10);
    }
}
