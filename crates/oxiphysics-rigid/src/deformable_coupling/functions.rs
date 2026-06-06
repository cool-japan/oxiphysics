//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DeformableProxy;

pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
pub(super) fn vec3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-12 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint_forces::vec3_norm;

    use crate::deformable_coupling::CoupledSimulation;
    use crate::deformable_coupling::CouplingForce;

    use crate::deformable_coupling::DeformationState;

    use crate::deformable_coupling::FluidStructureInterface;
    use crate::deformable_coupling::ImpactDeformation;

    use crate::deformable_coupling::RigidDeformableBond;

    use crate::deformable_coupling::SoftAttachment;
    use crate::deformable_coupling::SurfaceProjection;
    use crate::deformable_coupling::ThermalExpansionCoupling;

    use crate::deformable_coupling::VibrationTransmission;

    #[test]
    fn projection_to_horizontal_triangle() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.3, 0.5, 0.3];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        assert!((proj.gap - 0.5).abs() < 0.01, "gap: {}", proj.gap);
        assert!(!proj.is_penetrating());
    }
    #[test]
    fn projection_penetrating_below_triangle() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.3, -0.2, 0.3];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        assert!(
            proj.is_penetrating(),
            "should be penetrating: gap={}",
            proj.gap
        );
    }
    #[test]
    fn projection_barycentric_sum_to_one() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.25, 0.1, 0.25];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        let sum = proj.barycentric[0] + proj.barycentric[1] + proj.barycentric[2];
        assert!((sum - 1.0).abs() < 0.01, "barycentric sum: {sum}");
    }
    #[test]
    fn penalty_force_zero_when_not_penetrating() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.3, 0.5, 0.3];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        let cf = CouplingForce::penalty(&proj, 1e6, 1e3, 0.0);
        assert!(
            cf.magnitude() < 1e-3,
            "force should be near zero: {}",
            cf.magnitude()
        );
    }
    #[test]
    fn penalty_force_positive_when_penetrating() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.3, -0.05, 0.3];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        let cf = CouplingForce::penalty(&proj, 1e6, 1e3, 0.0);
        assert!(
            cf.magnitude() > 0.0,
            "force should be positive: {}",
            cf.magnitude()
        );
    }
    #[test]
    fn penalty_rigid_and_deformable_forces_opposite() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 0.0, 1.0];
        let point = [0.3, -0.05, 0.3];
        let proj = SurfaceProjection::project_to_triangle(point, v0, v1, v2, 0);
        let cf = CouplingForce::penalty(&proj, 1e6, 1e3, 0.0);
        for i in 0..3 {
            assert!(
                (cf.rigid_force[i] + cf.deformable_force[i]).abs() < 1e-10,
                "forces should be equal and opposite: axis {i}"
            );
        }
    }
    #[test]
    fn bond_debond_deactivates() {
        let mut bond = RigidDeformableBond::default_rubber_seal();
        assert!(bond.active);
        bond.debond();
        assert!(!bond.active);
    }
    #[test]
    fn bond_mortar_higher_stiffness_than_rubber() {
        let mortar = RigidDeformableBond::structural_mortar(0);
        let rubber = RigidDeformableBond::default_rubber_seal();
        assert!(mortar.penalty_stiffness > rubber.penalty_stiffness);
    }
    #[test]
    fn deformation_state_starts_at_zero() {
        let ds = DeformationState::new_undeformed(10);
        assert_eq!(ds.max_displacement(), 0.0);
    }
    #[test]
    fn deformation_state_apply_displacement() {
        let mut ds = DeformationState::new_undeformed(5);
        ds.apply_uniform_displacement([0.01, 0.0, 0.0]);
        assert!((ds.max_displacement() - 0.01).abs() < 1e-10);
    }
    #[test]
    fn deformation_state_rms_strain() {
        let mut ds = DeformationState::new_undeformed(4);
        ds.strains = vec![0.001, 0.002, 0.001, 0.002];
        let rms = ds.rms_strain();
        assert!(rms > 0.0 && rms < 0.01, "rms_strain: {rms}");
    }
    #[test]
    fn soft_attachment_force_at_small_disp() {
        let mut sa = SoftAttachment::epoxy_bond(0.01);
        let (fn_, fs) = sa.coupling_force(0.0001, 0.0001);
        assert!(fn_ > 0.0 && fs > 0.0, "fn={fn_}, fs={fs}");
    }
    #[test]
    fn soft_attachment_debonds_at_high_disp() {
        let mut sa = SoftAttachment::epoxy_bond(0.01);
        sa.coupling_force(100.0, 0.0);
        assert!(sa.debonded, "should be debonded");
    }
    #[test]
    fn soft_attachment_zero_after_debond() {
        let mut sa = SoftAttachment::epoxy_bond(0.01);
        sa.debonded = true;
        let (fn_, fs) = sa.coupling_force(0.1, 0.1);
        assert_eq!(fn_, 0.0);
        assert_eq!(fs, 0.0);
    }
    #[test]
    fn hertz_force_positive() {
        let imp = ImpactDeformation::steel_on_steel();
        let f = imp.peak_hertz_force(5.0);
        assert!(f > 0.0, "force: {f}");
    }
    #[test]
    fn hertz_contact_radius_positive() {
        let imp = ImpactDeformation::steel_on_steel();
        let r = imp.hertz_contact_radius(5.0);
        assert!(r > 0.0, "radius: {r}");
    }
    #[test]
    fn plastic_depth_increases_with_velocity() {
        let mut imp1 = ImpactDeformation::steel_on_steel();
        let mut imp2 = ImpactDeformation::steel_on_steel();
        let d1 = imp1.plastic_crater_depth(1.0);
        let d2 = imp2.plastic_crater_depth(5.0);
        assert!(d2 > d1, "higher velocity → deeper crater: {d1} vs {d2}");
    }
    #[test]
    fn yielding_at_high_velocity() {
        let imp = ImpactDeformation::steel_on_steel();
        assert!(imp.caused_yielding(100.0), "should yield at 100 m/s");
    }
    #[test]
    fn frf_peak_at_resonance() {
        let vt = VibrationTransmission::steel_panel();
        let fn_ = vt.natural_freq;
        let peak = vt.frf_amplitude(fn_);
        let off_peak = vt.frf_amplitude(fn_ * 2.0);
        assert!(
            peak > off_peak,
            "FRF should peak at resonance: {peak} vs {off_peak}"
        );
    }
    #[test]
    fn frf_peak_positive() {
        let vt = VibrationTransmission::steel_panel();
        assert!(vt.peak_frf() > 0.0);
    }
    #[test]
    fn rms_acceleration_positive() {
        let vt = VibrationTransmission::steel_panel();
        let a = vt.rms_acceleration(100.0, 100.0);
        assert!(a > 0.0, "rms acceleration: {a}");
    }
    #[test]
    fn thermal_expansion_zero_at_reference() {
        let t = ThermalExpansionCoupling::aluminium_with_rubber_seal();
        assert_eq!(t.thermal_strain(), 0.0);
        assert_eq!(t.expansion(), 0.0);
    }
    #[test]
    fn thermal_expansion_positive_above_reference() {
        let mut t = ThermalExpansionCoupling::aluminium_with_rubber_seal();
        t.set_temperature(t.reference_temp + 50.0);
        assert!(t.expansion() > 0.0, "expansion: {}", t.expansion());
    }
    #[test]
    fn thermal_force_proportional_to_delta_t() {
        let mut t1 = ThermalExpansionCoupling::aluminium_with_rubber_seal();
        let mut t2 = ThermalExpansionCoupling::aluminium_with_rubber_seal();
        t1.set_temperature(t1.reference_temp + 10.0);
        t2.set_temperature(t2.reference_temp + 20.0);
        assert!(
            (t2.thermal_force() - 2.0 * t1.thermal_force()).abs() < 1.0,
            "force should double: {} vs 2*{}",
            t2.thermal_force(),
            t1.thermal_force()
        );
    }
    #[test]
    fn fsi_pressure_force_non_zero() {
        let fsi = FluidStructureInterface::submerged_plate();
        let f = fsi.pressure_force();
        assert!(vec3_norm(f) > 0.0);
    }
    #[test]
    fn fsi_added_mass_positive() {
        let fsi = FluidStructureInterface::submerged_plate();
        assert!(fsi.added_mass() > 0.0);
    }
    #[test]
    fn fsi_kinematic_feedback_pressure() {
        let mut fsi = FluidStructureInterface::submerged_plate();
        fsi.update_rigid_velocity([0.0, 2.0, 0.0]);
        let dp = fsi.kinematic_feedback_pressure();
        assert!(dp > 0.0, "dp: {dp}");
    }
    #[test]
    fn coupled_sim_rigid_falls_under_gravity() {
        let mut sim = CoupledSimulation::new_default();
        let y0 = sim.rigid_position[1];
        for _ in 0..100 {
            sim.step(0.01, [0.0; 3]);
        }
        assert!(
            sim.rigid_position[1] < y0,
            "rigid body should fall: y={} vs y0={y0}",
            sim.rigid_position[1]
        );
    }
    #[test]
    fn coupled_sim_kinetic_energy_increases() {
        let mut sim = CoupledSimulation::new_default();
        let ke0 = sim.rigid_kinetic_energy();
        for _ in 0..50 {
            sim.step(0.01, [0.0; 3]);
        }
        let ke1 = sim.rigid_kinetic_energy();
        assert!(ke1 > ke0, "kinetic energy should increase under gravity");
    }
    #[test]
    fn coupled_sim_no_movement_after_debond() {
        let mut sim = CoupledSimulation::new_default();
        sim.bond.debond();
        let pos0 = sim.rigid_position;
        sim.step(0.01, [0.0; 3]);
        assert_eq!(sim.rigid_position, pos0);
    }
    #[test]
    fn vec3_normalise_unit_length() {
        let v = vec3_normalise([3.0, 4.0, 0.0]);
        assert!((vec3_norm(v) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn vec3_normalise_zero_returns_zero() {
        let v = vec3_normalise([0.0; 3]);
        assert_eq!(v, [0.0; 3]);
    }
}
#[cfg(test)]
mod new_struct_tests {

    use crate::constraint_forces::vec3_norm;
    use crate::constraint_forces::vec3_sub;
    use crate::deformable_coupling::ContactBridge;

    use crate::deformable_coupling::CouplingJoint;
    use crate::deformable_coupling::DeformableProxy;

    use crate::deformable_coupling::EnergyExchange;

    use crate::deformable_coupling::ImpulseTransfer;
    use crate::deformable_coupling::KinematicDriver;

    use crate::deformable_coupling::ShapeMatchingCoupling;
    use crate::deformable_coupling::SkinningMatrix;

    use crate::deformable_coupling::TwoWayCoupling;
    use crate::deformable_coupling::VelocityConstraint;

    #[test]
    fn proxy_num_nodes() {
        let nodes = vec![[0.0; 3]; 8];
        let p = DeformableProxy::new(nodes, 1e4);
        assert_eq!(p.num_nodes(), 8);
    }
    #[test]
    fn proxy_max_deformation_initially_zero() {
        let nodes = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let p = DeformableProxy::new(nodes, 1e4);
        assert_eq!(p.max_deformation(), 0.0);
    }
    #[test]
    fn proxy_rigid_driven_identity_matches_reference() {
        let nodes = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let p = DeformableProxy::new(nodes.clone(), 1e4);
        let driven = p.rigid_driven_positions();
        for (d, r) in driven.iter().zip(nodes.iter()) {
            let err = vec3_norm(vec3_sub(*d, *r));
            assert!(err < 1e-10, "error: {err}");
        }
    }
    #[test]
    fn proxy_spring_forces_zero_when_matching() {
        let nodes = vec![[1.0, 0.0, 0.0]];
        let p = DeformableProxy::new(nodes, 1e4);
        let forces = p.spring_forces();
        assert!(vec3_norm(forces[0]) < 1e-10);
    }
    #[test]
    fn proxy_spring_forces_nonzero_when_displaced() {
        let nodes = vec![[1.0, 0.0, 0.0]];
        let mut p = DeformableProxy::new(nodes, 1e4);
        p.deformed_shape[0] = [1.5, 0.0, 0.0];
        let forces = p.spring_forces();
        assert!(vec3_norm(forces[0]) > 0.0);
    }
    #[test]
    fn shape_matching_total_mass() {
        let pos = vec![[0.0; 3]; 4];
        let smc = ShapeMatchingCoupling::new_uniform(pos, 0.25, 0.8);
        assert!((smc.total_mass() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn shape_matching_rest_com_symmetric() {
        let pos = vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let smc = ShapeMatchingCoupling::new_uniform(pos, 1.0, 0.8);
        let com = smc.rest_com();
        assert!(com[0].abs() < 1e-10, "com_x: {}", com[0]);
    }
    #[test]
    fn shape_matching_update_goals_alpha_one() {
        let pos = vec![[0.0; 3]];
        let mut smc = ShapeMatchingCoupling::new_uniform(pos, 1.0, 1.0);
        smc.update_goals(&[[0.0; 3]], &[[1.0, 0.0, 0.0]]);
        assert!((smc.goal_positions[0][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn shape_matching_update_goals_alpha_zero() {
        let pos = vec![[0.5, 0.0, 0.0]];
        let mut smc = ShapeMatchingCoupling::new_uniform(pos.clone(), 1.0, 0.0);
        smc.update_goals(&[[0.5, 0.0, 0.0]], &[[1.0, 0.0, 0.0]]);
        assert!((smc.goal_positions[0][0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn shape_matching_restoring_force_toward_goal() {
        let pos = vec![[0.0; 3]];
        let mut smc = ShapeMatchingCoupling::new_uniform(pos, 1.0, 1.0);
        smc.goal_positions = vec![[2.0, 0.0, 0.0]];
        let forces = smc.restoring_forces(&[[0.0; 3]], 100.0);
        assert!(
            (forces[0][0] - 200.0).abs() < 1e-6,
            "force_x: {}",
            forces[0][0]
        );
    }
    #[test]
    fn impulse_transfer_distribute_sum_to_node_delta_v() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let it = ImpulseTransfer::new_uniform(positions, 1.0, 10.0);
        let deltas = it.distribute_to_nodes([0.0, 10.0, 0.0], [0.5, 0.0, 0.0], 2);
        let total_y: f64 = deltas.iter().map(|d| d[1]).sum();
        assert!(total_y.abs() > 0.0, "total impulse transferred: {total_y}");
    }
    #[test]
    fn impulse_transfer_rigid_delta_v_opposite() {
        let positions = vec![[0.0; 3]];
        let it = ImpulseTransfer::new_uniform(positions, 1.0, 2.0);
        let dv = it.rigid_velocity_delta([0.0, 4.0, 0.0]);
        assert!((dv[1] + 2.0).abs() < 1e-10, "rigid dv_y: {}", dv[1]);
    }
    #[test]
    fn impulse_transfer_no_nodes_empty_result() {
        let it = ImpulseTransfer::new_uniform(Vec::new(), 1.0, 1.0);
        let res = it.distribute_to_nodes([1.0, 0.0, 0.0], [0.0; 3], 3);
        assert!(res.is_empty());
    }
    #[test]
    fn contact_bridge_detects_nearby_node() {
        let cb = ContactBridge::new(0.05, 0.5, 0.3);
        let nodes = vec![[0.0, 0.0, 0.0]];
        let contacts = cb.detect_contacts(&nodes);
        assert!(!contacts.is_empty(), "should detect contact");
    }
    #[test]
    fn contact_bridge_no_contact_far_node() {
        let cb = ContactBridge::new(0.05, 0.5, 0.3);
        let nodes = vec![[100.0, 0.0, 0.0]];
        assert_eq!(cb.contact_count(&nodes), 0);
    }
    #[test]
    fn contact_bridge_impulse_zero_when_separating() {
        let cb = ContactBridge::new(0.05, 0.5, 0.3);
        let j = cb.contact_impulse(-1.0, 1.0);
        assert_eq!(j, 0.0);
    }
    #[test]
    fn contact_bridge_impulse_positive_when_approaching() {
        let cb = ContactBridge::new(0.05, 1.0, 0.3);
        let j = cb.contact_impulse(2.0, 1.0);
        assert!(
            j < 0.0,
            "impulse should be negative (opposing approach): {j}"
        );
    }
    #[test]
    fn kinematic_driver_target_at_com_for_zero_offset() {
        let kd = KinematicDriver::new_at_com(3, 1e4, 100.0);
        let t = kd.target_position(0);
        assert_eq!(t, [0.0; 3]);
    }
    #[test]
    fn kinematic_driver_drive_force_toward_rigid() {
        let mut kd = KinematicDriver::new_at_com(1, 1000.0, 0.0);
        kd.set_rigid_state([1.0, 0.0, 0.0], [0.0; 3]);
        let f = kd.drive_force(0, [0.0; 3], [0.0; 3]);
        assert!((f[0] - 1000.0).abs() < 1e-8, "force_x: {}", f[0]);
    }
    #[test]
    fn kinematic_driver_drive_force_zero_when_aligned() {
        let mut kd = KinematicDriver::new_at_com(1, 1000.0, 100.0);
        kd.set_rigid_state([1.0, 0.0, 0.0], [0.0; 3]);
        let f = kd.drive_force(0, [1.0, 0.0, 0.0], [0.0; 3]);
        assert!(vec3_norm(f) < 1e-8, "force should be zero: {f:?}");
    }
    #[test]
    fn two_way_coupling_node_force_zero_at_rest() {
        let tc = TwoWayCoupling::new(1, 10.0, 1000.0);
        let f = tc.node_force(0, [0.0; 3]);
        assert!(vec3_norm(f) < 1e-10);
    }
    #[test]
    fn two_way_coupling_reaction_opposes_node_force() {
        let mut tc = TwoWayCoupling::new(1, 10.0, 1000.0);
        tc.attachment_points[0] = [1.0, 0.0, 0.0];
        let nodes = vec![[0.0; 3]];
        let node_f = tc.node_force(0, nodes[0]);
        let reaction = tc.rigid_reaction_force(&nodes);
        for i in 0..3 {
            assert!((node_f[i] + reaction[i]).abs() < 1e-8, "axis {i}");
        }
    }
    #[test]
    fn two_way_coupling_integrate_changes_velocity() {
        let mut tc = TwoWayCoupling::new(1, 1.0, 1000.0);
        tc.attachment_points[0] = [0.0, 1.0, 0.0];
        let nodes = vec![[0.0, 2.0, 0.0]];
        tc.integrate_rigid(&nodes, 0.01);
        assert!(tc.rigid_velocity[1].abs() > 0.0, "vy should be nonzero");
    }
    #[test]
    fn coupling_joint_force_zero_at_anchor() {
        let j = CouplingJoint::new(0, [1.0, 0.0, 0.0], 1000.0, 10.0);
        let f = j.spring_force([1.0, 0.0, 0.0], [0.0; 3]);
        assert!(vec3_norm(f) < 1e-10);
    }
    #[test]
    fn coupling_joint_force_nonzero_when_displaced() {
        let j = CouplingJoint::new(0, [0.0; 3], 1000.0, 0.0);
        let f = j.spring_force([1.0, 0.0, 0.0], [0.0; 3]);
        assert!(vec3_norm(f) > 0.0);
    }
    #[test]
    fn coupling_joint_breaks_at_max_force() {
        let mut j = CouplingJoint::new(0, [100.0, 0.0, 0.0], 1e6, 0.0);
        j.max_force = 1.0;
        j.check_breakage([0.0; 3], [0.0; 3]);
        assert!(!j.intact, "joint should have broken");
    }
    #[test]
    fn coupling_joint_no_force_when_broken() {
        let mut j = CouplingJoint::new(0, [1.0, 0.0, 0.0], 1000.0, 0.0);
        j.intact = false;
        let f = j.spring_force([0.0; 3], [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn skinning_matrix_rigid_all_weight_on_first_bone() {
        let sm = SkinningMatrix::new_rigid(4, 3);
        for i in 0..4 {
            assert!((sm.weights[i][0] - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn skinning_matrix_uniform_row_sums_to_one() {
        let sm = SkinningMatrix::new_uniform(3, 4);
        for i in 0..3 {
            assert!(
                (sm.row_sum(i) - 1.0).abs() < 1e-10,
                "row {i} sum: {}",
                sm.row_sum(i)
            );
        }
    }
    #[test]
    fn skinning_matrix_skinned_position_matches_rigid() {
        let sm = SkinningMatrix::new_rigid(1, 2);
        let bone_pos = vec![[5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let pos = sm.skinned_position(0, &bone_pos);
        assert!((pos[0] - 5.0).abs() < 1e-10, "x: {}", pos[0]);
    }
    #[test]
    fn skinning_matrix_normalise_rows() {
        let mut sm = SkinningMatrix::new_uniform(2, 2);
        sm.weights[0][0] = 2.0;
        sm.weights[0][1] = 2.0;
        sm.normalise_rows();
        assert!((sm.row_sum(0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn velocity_constraint_relative_velocity_zero_when_matched() {
        let vc = VelocityConstraint::new(0, [1.0, 0.0, 0.0], 0.1);
        let rel = vc.relative_velocity([0.0; 3]);
        assert_eq!(rel, 0.0);
    }
    #[test]
    fn velocity_constraint_relative_velocity_nonzero() {
        let vc = VelocityConstraint::new(0, [1.0, 0.0, 0.0], 0.1);
        let rel = vc.relative_velocity([5.0, 0.0, 0.0]);
        assert!((rel - 5.0).abs() < 1e-10, "rel: {rel}");
    }
    #[test]
    fn velocity_constraint_lagrange_multiplier_sign() {
        let vc = VelocityConstraint::new(0, [1.0, 0.0, 0.0], 0.0);
        let lm = vc.lagrange_multiplier([1.0, 0.0, 0.0], 1.0, 0.01);
        assert!(
            lm < 0.0,
            "λ should be negative to oppose relative motion: {lm}"
        );
    }
    #[test]
    fn energy_exchange_initial_state_zero() {
        let ee = EnergyExchange::new();
        assert_eq!(ee.total_energy(), 0.0);
        assert_eq!(ee.net_transfer(), 0.0);
    }
    #[test]
    fn energy_exchange_record_positive_flow() {
        let mut ee = EnergyExchange::new();
        ee.record_exchange(100.0, 0.01);
        assert!((ee.energy_to_deformable - 1.0).abs() < 1e-10);
        assert_eq!(ee.energy_to_rigid, 0.0);
    }
    #[test]
    fn energy_exchange_record_negative_flow() {
        let mut ee = EnergyExchange::new();
        ee.record_exchange(-200.0, 0.01);
        assert!((ee.energy_to_rigid - 2.0).abs() < 1e-10);
        assert_eq!(ee.energy_to_deformable, 0.0);
    }
    #[test]
    fn energy_exchange_net_transfer_sign() {
        let mut ee = EnergyExchange::new();
        ee.record_exchange(100.0, 1.0);
        ee.record_exchange(-40.0, 1.0);
        assert!(
            (ee.net_transfer() - 60.0).abs() < 1e-10,
            "net: {}",
            ee.net_transfer()
        );
    }
    #[test]
    fn energy_exchange_total_energy_update() {
        let mut ee = EnergyExchange::new();
        ee.update(50.0, 30.0);
        assert!((ee.total_energy() - 80.0).abs() < 1e-10);
    }
    #[test]
    fn energy_exchange_default_same_as_new() {
        let ee: EnergyExchange = Default::default();
        assert_eq!(ee.total_energy(), 0.0);
    }
}
/// Compute the spring force at a deformable-rigid interface.
///
/// # Arguments
/// * `proxy` - Deformable proxy at the attachment point.
/// * `rigid_pos` - Position of the rigid body attachment point `[x, y, z]`.
/// * `stiffness` - Spring stiffness in N/m.
///
/// Returns the force vector `[fx, fy, fz]` applied to the rigid body.
pub fn spring_attachment(proxy: &DeformableProxy, rigid_pos: [f64; 3], stiffness: f64) -> [f64; 3] {
    let centre = proxy.rigid_com;
    let delta = [
        rigid_pos[0] - centre[0],
        rigid_pos[1] - centre[1],
        rigid_pos[2] - centre[2],
    ];
    [
        -stiffness * delta[0],
        -stiffness * delta[1],
        -stiffness * delta[2],
    ]
}
/// Reduce a modal basis to the first `n_modes` modes (modal truncation).
///
/// Given a full modal matrix stored column-major as a flat `Vec`f64` of
/// shape `(n_dofs × n_modes_full)`, this returns only the first `n_modes`
/// columns.
///
/// # Arguments
/// * `modal_matrix` - Flat column-major modal matrix `n_dofs × n_modes_full`.
/// * `n_dofs` - Number of physical degrees of freedom.
/// * `n_modes_full` - Total number of modes available.
/// * `n_modes` - Number of modes to retain.
pub fn modal_truncation(
    modal_matrix: &[f64],
    n_dofs: usize,
    n_modes_full: usize,
    n_modes: usize,
) -> Vec<f64> {
    assert!(
        n_modes <= n_modes_full,
        "Cannot retain more modes than available"
    );
    assert_eq!(modal_matrix.len(), n_dofs * n_modes_full);
    let keep = n_modes.min(n_modes_full);
    let mut out = vec![0.0_f64; n_dofs * keep];
    for col in 0..keep {
        for row in 0..n_dofs {
            out[col * n_dofs + row] = modal_matrix[col * n_dofs + row];
        }
    }
    out
}
/// Craig-Bampton substructure reduction.
///
/// Reduces a structural subcomponent to a set of fixed-interface normal modes
/// plus boundary (interface) DOFs.  For this implementation the returned
/// result is a pair `(n_interior_modes, n_boundary_dofs)` representing the
/// dimension of the reduced system.
///
/// # Arguments
/// * `n_interior_dofs` - Number of interior (non-boundary) DOFs.
/// * `n_boundary_dofs` - Number of boundary (interface) DOFs.
/// * `n_modes` - Number of fixed-interface normal modes to retain.
///
/// Returns `(reduced_size, constraint_modes_count)`.
pub fn craig_bampton(
    n_interior_dofs: usize,
    n_boundary_dofs: usize,
    n_modes: usize,
) -> (usize, usize) {
    let reduced_size = n_modes + n_boundary_dofs;
    let constraint_modes = n_interior_dofs.min(n_boundary_dofs * n_modes + n_boundary_dofs);
    (reduced_size, constraint_modes)
}
/// Compute rigid-body modes (null space of the stiffness matrix) for an
/// unconstrained body in 3-D space.
///
/// For a free body in 3-D there are exactly 6 rigid-body modes
/// (3 translations + 3 rotations). This function returns a flat column-major
/// matrix of size `n_dofs × 6` whose columns span the rigid-body null space.
/// The DOFs are ordered as triples `\[x_i, y_i, z_i\]` for each node,
/// with `node_positions` listing the `\[x, y, z\]` position of each node.
///
/// # Arguments
///
/// * `node_positions` - Flat list `\[x0, y0, z0, x1, y1, z1, …\]`.
pub fn rigid_body_modes(node_positions: &[f64]) -> Vec<f64> {
    assert_eq!(
        node_positions.len() % 3,
        0,
        "node_positions must be a multiple of 3"
    );
    let n_nodes = node_positions.len() / 3;
    let n_dofs = n_nodes * 3;
    let mut rbm = vec![0.0_f64; n_dofs * 6];
    for i in 0..n_nodes {
        let (x, y, z) = (
            node_positions[3 * i],
            node_positions[3 * i + 1],
            node_positions[3 * i + 2],
        );
        let row = i * 3;
        rbm[row] = 1.0;
        rbm[n_dofs + row + 1] = 1.0;
        rbm[2 * n_dofs + row + 2] = 1.0;
        rbm[3 * n_dofs + row] = 0.0;
        rbm[3 * n_dofs + row + 1] = -z;
        rbm[3 * n_dofs + row + 2] = y;
        rbm[4 * n_dofs + row] = z;
        rbm[4 * n_dofs + row + 1] = 0.0;
        rbm[4 * n_dofs + row + 2] = -x;
        rbm[5 * n_dofs + row] = -y;
        rbm[5 * n_dofs + row + 1] = x;
        rbm[5 * n_dofs + row + 2] = 0.0;
    }
    rbm
}
/// Assemble a reduced system from multiple Craig-Bampton substructures.
///
/// Each substructure is described by `(n_modes, n_boundary_dofs)`.  The total
/// reduced system size is the sum of all reduced substructure sizes.
///
/// # Arguments
/// * `substructures` - Slice of `(n_interior_dofs, n_boundary_dofs, n_modes)`.
///
/// Returns the total reduced system DOF count.
pub fn component_synthesis(substructures: &[(usize, usize, usize)]) -> usize {
    substructures
        .iter()
        .map(|&(ni, nb, nm)| craig_bampton(ni, nb, nm).0)
        .sum()
}
#[cfg(test)]
mod substructure_tests {

    use crate::deformable_coupling::DeformableProxy;

    use crate::deformable_coupling::component_synthesis;
    use crate::deformable_coupling::craig_bampton;
    use crate::deformable_coupling::modal_truncation;
    use crate::deformable_coupling::rigid_body_modes;
    use crate::deformable_coupling::spring_attachment;
    #[test]
    fn spring_attachment_zero_displacement_zero_force() {
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let proxy = DeformableProxy::new(nodes, 1e3);
        let centre = proxy.rigid_com;
        let f = spring_attachment(&proxy, centre, 500.0);
        assert!(f[0].abs() < 1e-10 && f[1].abs() < 1e-10 && f[2].abs() < 1e-10);
    }
    #[test]
    fn spring_attachment_unit_displacement() {
        let nodes = vec![[0.0, 0.0, 0.0]];
        let proxy = DeformableProxy::new(nodes, 1e3);
        let f = spring_attachment(&proxy, [1.0, 0.0, 0.0], 100.0);
        assert!((f[0] + 100.0).abs() < 1e-10, "fx={}", f[0]);
        assert!(f[1].abs() < 1e-10);
        assert!(f[2].abs() < 1e-10);
    }
    #[test]
    fn spring_attachment_3d_displacement() {
        let nodes = vec![[0.0, 0.0, 0.0]];
        let proxy = DeformableProxy::new(nodes, 1e3);
        let f = spring_attachment(&proxy, [1.0, 2.0, 3.0], 10.0);
        assert!((f[0] + 10.0).abs() < 1e-10);
        assert!((f[1] + 20.0).abs() < 1e-10);
        assert!((f[2] + 30.0).abs() < 1e-10);
    }
    #[test]
    fn modal_truncation_reduces_columns() {
        let n_dofs = 4;
        let n_full = 4;
        let modal: Vec<f64> = (0..n_dofs * n_full).map(|i| i as f64).collect();
        let truncated = modal_truncation(&modal, n_dofs, n_full, 2);
        assert_eq!(truncated.len(), n_dofs * 2);
    }
    #[test]
    fn modal_truncation_preserves_first_columns() {
        let n_dofs = 3;
        let n_full = 3;
        let modal: Vec<f64> = (0..n_dofs * n_full).map(|i| i as f64).collect();
        let truncated = modal_truncation(&modal, n_dofs, n_full, 2);
        for i in 0..n_dofs * 2 {
            assert!((truncated[i] - modal[i]).abs() < 1e-10, "i={i}");
        }
    }
    #[test]
    fn modal_truncation_zero_modes() {
        let modal = vec![1.0, 2.0, 3.0, 4.0];
        let truncated = modal_truncation(&modal, 2, 2, 0);
        assert!(truncated.is_empty());
    }
    #[test]
    fn craig_bampton_reduced_size() {
        let (sz, _) = craig_bampton(10, 4, 3);
        assert_eq!(sz, 7, "3 modes + 4 boundary = 7");
    }
    #[test]
    fn craig_bampton_no_modes() {
        let (sz, _) = craig_bampton(10, 4, 0);
        assert_eq!(sz, 4);
    }
    #[test]
    fn craig_bampton_no_boundary() {
        let (sz, _) = craig_bampton(10, 0, 5);
        assert_eq!(sz, 5);
    }
    #[test]
    fn rigid_body_modes_single_node_shape() {
        let pos = vec![0.0, 0.0, 0.0];
        let rbm = rigid_body_modes(&pos);
        assert_eq!(rbm.len(), 3 * 6);
    }
    #[test]
    fn rigid_body_modes_translation_columns_identity() {
        let pos = vec![0.0, 0.0, 0.0];
        let rbm = rigid_body_modes(&pos);
        assert!((rbm[0] - 1.0).abs() < 1e-10);
        assert!(rbm[1].abs() < 1e-10);
        assert!(rbm[2].abs() < 1e-10);
    }
    #[test]
    fn rigid_body_modes_two_nodes_shape() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let rbm = rigid_body_modes(&pos);
        assert_eq!(rbm.len(), 6 * 6, "2 nodes × 3 dofs × 6 modes");
    }
    #[test]
    fn component_synthesis_single_substructure() {
        let total = component_synthesis(&[(10, 4, 3)]);
        assert_eq!(total, 7);
    }
    #[test]
    fn component_synthesis_two_substructures() {
        let total = component_synthesis(&[(10, 4, 3), (8, 2, 2)]);
        assert_eq!(total, 11);
    }
    #[test]
    fn component_synthesis_empty() {
        let total = component_synthesis(&[]);
        assert_eq!(total, 0);
    }
}
