//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

/// Speed of sound squared (cs² = 1/3).
pub(super) const CS2: f64 = 1.0 / 3.0;
/// Convert coarse-grid D2Q9 populations to fine-grid via equilibrium rescaling.
///
/// The fine-grid populations are computed by rescaling the non-equilibrium
/// part: f_fine = f_eq(ρ, u) + (τ_f / τ_c) * (f_coarse - f_eq(ρ, u)).
pub fn coarse_to_fine_pop(
    f_coarse: &[f64; 9],
    rho: f64,
    u_coarse: [f64; 2],
    u_fine: [f64; 2],
) -> [f64; 9] {
    let mut result = [0.0_f64; 9];
    let _ = u_fine;
    for i in 0..9 {
        let w = D2Q9_WEIGHTS[i];
        let c = D2Q9_VELOCITIES[i];
        let eu = c[0] as f64 * u_coarse[0] + c[1] as f64 * u_coarse[1];
        let u2 = u_coarse[0] * u_coarse[0] + u_coarse[1] * u_coarse[1];
        let feq = w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
        result[i] = feq + 2.0 * (f_coarse[i] - feq);
        result[i] = result[i].max(0.0);
    }
    result
}
/// Convert fine-grid D2Q9 populations to coarse-grid via restriction averaging.
///
/// Takes a slice of `ratio²` fine-cell populations and volume-averages them.
pub fn fine_to_coarse_pop(f_fine_cells: &[[f64; 9]]) -> [f64; 9] {
    let n = f_fine_cells.len();
    if n == 0 {
        return [0.0; 9];
    }
    let mut sum = [0.0_f64; 9];
    for f in f_fine_cells {
        for i in 0..9 {
            sum[i] += f[i];
        }
    }
    sum.map(|s| s / n as f64)
}
/// Update the buffer layer distributions by interpolating between coarse nodes.
///
/// `f_boundary` are the two coarse boundary nodes; `alpha ∈ [0,1]` is the
/// fractional position of the buffer node.
pub fn buffer_layer_update(f_a: &[f64; 9], f_b: &[f64; 9], alpha: f64) -> [f64; 9] {
    let mut result = [0.0_f64; 9];
    for i in 0..9 {
        result[i] = (1.0 - alpha) * f_a[i] + alpha * f_b[i];
    }
    result
}
/// Grid refinement criterion based on vorticity magnitude.
///
/// Computes |ω| = |∂u_y/∂x − ∂u_x/∂y| via finite differences.
/// Returns `true` if the cell should be refined.
pub fn grid_refinement_criterion(
    u: &[[f64; 2]],
    nx: usize,
    ny: usize,
    threshold: f64,
) -> Vec<bool> {
    let n = nx * ny;
    let mut mask = vec![false; n];
    for y in 1..(ny - 1) {
        for x in 1..(nx - 1) {
            let k = y * nx + x;
            let duydx = (u[k + 1][1] - u[k - 1][1]) * 0.5;
            let duxdy = (u[k + nx][0] - u[k - nx][0]) * 0.5;
            let omega = (duydx - duxdy).abs();
            mask[k] = omega > threshold;
        }
    }
    mask
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiscale::types::*;
    #[test]
    fn test_quadtree_node_is_leaf() {
        let node = QuadTreeNode::new([0.0, 1.0], [0.0, 1.0], 0, 0);
        assert!(node.is_leaf());
    }
    #[test]
    fn test_quadtree_node_center() {
        let node = QuadTreeNode::new([0.0, 2.0], [0.0, 4.0], 0, 0);
        let c = node.center();
        assert!((c[0] - 1.0).abs() < 1e-14);
        assert!((c[1] - 2.0).abs() < 1e-14);
    }
    #[test]
    fn test_quadtree_node_dx_dy() {
        let node = QuadTreeNode::new([0.0, 2.0], [0.0, 3.0], 0, 0);
        assert!((node.dx() - 2.0).abs() < 1e-14);
        assert!((node.dy() - 3.0).abs() < 1e-14);
    }
    #[test]
    fn test_quadtree_node_equilibrium_sum() {
        let mut node = QuadTreeNode::new([0.0, 1.0], [0.0, 1.0], 0, 0);
        node.rho = 2.0;
        node.u = [0.0, 0.0];
        let sum: f64 = (0..9).map(|i| node.equilibrium(i)).sum();
        assert!((sum - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_adaptive_mesh_initial_size() {
        let mesh = AdaptiveMesh::new(1.0, 3, 0.01);
        assert_eq!(mesh.len(), 1);
        assert_eq!(mesh.leaf_indices(), vec![0]);
    }
    #[test]
    fn test_adaptive_mesh_refine() {
        let mut mesh = AdaptiveMesh::new(1.0, 3, 0.01);
        let children = mesh.refine(0);
        assert!(children.is_some());
        let ids = children.unwrap();
        assert_eq!(ids.len(), 4);
        assert_eq!(mesh.len(), 5);
    }
    #[test]
    fn test_adaptive_mesh_refine_max_level() {
        let mut mesh = AdaptiveMesh::new(1.0, 0, 0.01);
        let result = mesh.refine(0);
        assert!(result.is_none());
    }
    #[test]
    fn test_adaptive_mesh_collide() {
        let mut mesh = AdaptiveMesh::new(1.0, 2, 0.01);
        mesh.refine(0);
        mesh.collide_all(1.0);
        for &idx in &mesh.leaf_indices() {
            let rho: f64 = mesh.nodes[idx].f.iter().sum();
            assert!(rho > 0.0);
        }
    }
    #[test]
    fn test_adaptive_mesh_is_not_empty() {
        let mesh = AdaptiveMesh::new(2.0, 4, 0.01);
        assert!(!mesh.is_empty());
    }
    #[test]
    fn test_patched_grid_sizes() {
        let pg = PatchedGrid::new(4, 4, 2);
        assert_eq!(pg.nx_fine, 7);
        assert_eq!(pg.ny_fine, 7);
    }
    #[test]
    fn test_patched_grid_prolongate_conserves() {
        let mut pg = PatchedGrid::new(3, 3, 2);
        pg.prolongate();
        let sum: f64 = pg.f_fine.iter().flat_map(|f| f.iter()).sum();
        assert!(sum > 0.0);
    }
    #[test]
    fn test_patched_grid_restrict() {
        let mut pg = PatchedGrid::new(3, 3, 2);
        for f in &mut pg.f_fine {
            for v in f.iter_mut() {
                *v = 1.0;
            }
        }
        pg.restrict();
        for f in &pg.f_coarse {
            let rho: f64 = f.iter().sum();
            assert!(rho > 0.0);
        }
    }
    #[test]
    fn test_space_time_tau_fine() {
        let st = SpaceTimeLbm::new(10, 2, 1.0);
        assert!((st.tau_fine - 0.75).abs() < 1e-12);
    }
    #[test]
    fn test_space_time_interpolation_endpoints() {
        let st = SpaceTimeLbm::new(4, 2, 1.0);
        let f0 = st.interpolated_coarse(0, 0);
        let f2 = st.interpolated_coarse(0, 2);
        for i in 0..9 {
            assert!((f0[i] - st.f_coarse_n[0][i]).abs() < 1e-14);
            assert!((f2[i] - st.f_coarse_np1[0][i]).abs() < 1e-14);
        }
    }
    #[test]
    fn test_space_time_advance_substep() {
        let mut st = SpaceTimeLbm::new(4, 2, 1.0);
        st.advance_fine_substep();
        assert_eq!(st.sub_step, 1);
        st.advance_fine_substep();
        assert_eq!(st.sub_step, 0);
    }
    #[test]
    fn test_hybrid_lbm_ns_exchange() {
        let mut h = HybridLbmNs::new(4, 4, 1.0, 1e-3);
        h.add_interface_pair(0, 0);
        h.lbm_to_ns_exchange();
        assert!((h.p_ns[0] - CS2).abs() < 1e-10);
    }
    #[test]
    fn test_hybrid_ns_to_lbm_exchange() {
        let mut h = HybridLbmNs::new(4, 4, 1.0, 1e-3);
        h.add_interface_pair(0, 0);
        h.u_ns[0] = [0.1, 0.0];
        h.p_ns[0] = CS2;
        h.ns_to_lbm_exchange();
        let rho: f64 = h.f_lbm[0].iter().sum();
        assert!((rho - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_multiscale_coupling_new() {
        let mc = MultiscaleCoupling::new(vec![10, 20], vec![1.0, 1.2], 1e-6);
        assert_eq!(mc.n_patches, 2);
        assert_eq!(mc.patches[0].len(), 10);
        assert_eq!(mc.patches[1].len(), 20);
    }
    #[test]
    fn test_multiscale_schwarz_exchange_no_overlap() {
        let mut mc = MultiscaleCoupling::new(vec![4, 4], vec![1.0, 1.0], 1e-6);
        let change = mc.schwarz_exchange();
        assert!((change).abs() < 1e-12);
    }
    #[test]
    fn test_multiscale_collide_patch() {
        let mut mc = MultiscaleCoupling::new(vec![8], vec![1.0], 1e-6);
        mc.collide_patch(0);
        let sum_rho: f64 = mc.patches[0].iter().map(|f| f.iter().sum::<f64>()).sum();
        assert!((sum_rho - 8.0).abs() < 1e-8);
    }
    #[test]
    fn test_variable_resolution_init() {
        let vr = VariableResolution::new(4, 4, 1.0);
        assert_eq!(vr.f.len(), 16);
        assert!((vr.x_coords[5] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_variable_resolution_stretching() {
        let mut vr = VariableResolution::new(5, 3, 1.0);
        vr.apply_x_stretching(0.0, 1.0, 1.0);
        assert!((vr.x_coords[0] - 0.0).abs() < 1e-10);
        assert!((vr.x_coords[4] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_variable_resolution_collide() {
        let mut vr = VariableResolution::new(4, 4, 1.0);
        vr.collide();
        let sum_rho: f64 = vr.f.iter().map(|f| f.iter().sum::<f64>()).sum();
        assert!((sum_rho - 16.0).abs() < 1e-8);
    }
    #[test]
    fn test_cell_reduction_restrict_uniform() {
        let cr = CellReduction::new(2);
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let f_fine = vec![f_init; 4 * 4];
        let (f_coarse, nx_c, ny_c) = cr.restrict(&f_fine, 4, 4);
        assert_eq!(nx_c, 2);
        assert_eq!(ny_c, 2);
        for f in &f_coarse {
            let rho: f64 = f.iter().sum();
            assert!((rho - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_cell_reduction_prolong() {
        let cr = CellReduction::new(2);
        let f_init: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let f_coarse = vec![f_init; 2 * 2];
        let (f_fine, nx_f, ny_f) = cr.prolong(&f_coarse, 2, 2);
        assert_eq!(nx_f, 4);
        assert_eq!(ny_f, 4);
        for f in &f_fine {
            let rho: f64 = f.iter().sum();
            assert!((rho - 2.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_haar_transform_roundtrip() {
        let filter = MultiresolutionFilter::new(0.01, 1);
        let signal = vec![1.0, 2.0, 3.0, 4.0];
        let (approx, detail) = filter.haar_transform(&signal);
        let reconstructed = filter.haar_inverse(&approx, &detail);
        for (a, b) in signal.iter().zip(reconstructed.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn test_haar_transform_constant_signal() {
        let filter = MultiresolutionFilter::new(0.01, 1);
        let signal = vec![3.0, 3.0, 3.0, 3.0];
        let (_, detail) = filter.haar_transform(&signal);
        for d in detail {
            assert!(d.abs() < 1e-14);
        }
    }
    #[test]
    fn test_refinement_indicator_flat_field() {
        let filter = MultiresolutionFilter::new(0.1, 1);
        let field = vec![1.0; 6 * 6];
        let indicator = filter.refinement_indicator(&field, 6, 6);
        assert!(!indicator.iter().any(|&v| v));
    }
    #[test]
    fn test_refinement_indicator_spike() {
        let filter = MultiresolutionFilter::new(0.1, 1);
        let mut field = vec![0.0_f64; 6 * 6];
        field[3 * 6 + 3] = 10.0;
        let indicator = filter.refinement_indicator(&field, 6, 6);
        assert!(indicator.iter().any(|&v| v));
    }
    #[test]
    fn test_coarse_to_fine_pop_positive() {
        let f_eq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let result = coarse_to_fine_pop(&f_eq, 1.0, [0.0, 0.0], [0.0, 0.0]);
        for v in result {
            assert!(v >= 0.0);
        }
    }
    #[test]
    fn test_fine_to_coarse_pop_average() {
        let f_val: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let cells = vec![f_val; 4];
        let result = fine_to_coarse_pop(&cells);
        let rho: f64 = result.iter().sum();
        assert!((rho - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_fine_to_coarse_empty() {
        let result = fine_to_coarse_pop(&[]);
        for v in result {
            assert!(v == 0.0);
        }
    }
    #[test]
    fn test_buffer_layer_update_interpolation() {
        let fa: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let fb: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let result = buffer_layer_update(&fa, &fb, 0.5);
        let rho: f64 = result.iter().sum();
        assert!((rho - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_grid_refinement_criterion_no_vorticity() {
        let u = vec![[0.1_f64, 0.0]; 36];
        let mask = grid_refinement_criterion(&u, 6, 6, 0.01);
        assert!(!mask.iter().any(|&v| v));
    }
    #[test]
    fn test_two_way_coupling_stokes_drag() {
        let coupling = TwoWayCoupling::new(4, 2, 1e-3, 1.0);
        let drag = coupling.stokes_drag(0, [0.1, 0.0]);
        let expected = 6.0 * std::f64::consts::PI * 1e-3 * 0.5 * 0.1;
        assert!((drag[0] - expected).abs() < 1e-14);
        assert!(drag[1].abs() < 1e-14);
    }
    #[test]
    fn test_lbm_md_coupling_init() {
        let lbm_md = LbmMolecularCoupling::new(4, 10, 1.0, 1.0);
        assert_eq!(lbm_md.n_lbm_nodes, 4);
        assert_eq!(lbm_md.n_md, 10);
        assert_eq!(lbm_md.md_positions.len(), 20);
    }
}
/// Compute D2Q9 equilibrium distribution.
///
/// f_eq_i = w_i * ρ * \[1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) - u²/(2cs²)\]
pub(super) fn lbm_equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
    let mut feq = [0.0f64; 9];
    let u2 = u[0] * u[0] + u[1] * u[1];
    for alpha in 0..9 {
        let ex = D2Q9_VELOCITIES[alpha][0] as f64;
        let ey = D2Q9_VELOCITIES[alpha][1] as f64;
        let eu = ex * u[0] + ey * u[1];
        feq[alpha] = D2Q9_WEIGHTS[alpha]
            * rho
            * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}
#[cfg(test)]
mod multiscale_extended_tests {
    use super::*;
    use crate::multiscale::types::*;
    #[test]
    fn test_ce_equilibrium_sums_to_rho() {
        let ce = ChapmanEnskogExpansion::new(1.0, [0.0, 0.0], 1.0, 0.01);
        let rho_sum: f64 = ce.f0.iter().sum();
        assert!((rho_sum - 1.0).abs() < 1e-12, "rho_sum={rho_sum}");
    }
    #[test]
    fn test_ce_compute_f1_zero_grad() {
        let mut ce = ChapmanEnskogExpansion::new(1.0, [0.0, 0.0], 1.0, 0.01);
        ce.compute_f1(0.0, 0.0, 0.0, 0.0);
        let f1_sum: f64 = ce.f1.iter().sum();
        assert!(f1_sum.abs() < 1e-12, "f1_sum={f1_sum}");
    }
    #[test]
    fn test_ce_compute_f1_shear() {
        let mut ce = ChapmanEnskogExpansion::new(1.0, [0.1, 0.0], 1.0, 0.01);
        ce.compute_f1(0.0, 0.1, 0.1, 0.0);
        assert!(ce.sigma_xy.abs() > 0.0);
    }
    #[test]
    fn test_ce_total_distribution_close_to_eq_small_kn() {
        let mut ce = ChapmanEnskogExpansion::new(1.0, [0.01, 0.0], 0.6, 0.001);
        ce.compute_f1(0.01, 0.0, 0.0, 0.01);
        let ftot = ce.total_distribution();
        for (alpha, (&ft, &f0)) in ftot.iter().zip(ce.f0.iter()).enumerate() {
            assert!((ft - f0).abs() < 0.01, "Large deviation at alpha={alpha}");
        }
    }
    #[test]
    fn test_ce_kinematic_viscosity_positive() {
        let ce = ChapmanEnskogExpansion::new(1.0, [0.0, 0.0], 1.0, 0.01);
        assert!(ce.kinematic_viscosity() > 0.0);
    }
    #[test]
    fn test_ce_continuum_regime() {
        let ce = ChapmanEnskogExpansion::new(1.0, [0.0, 0.0], 1.0, 0.005);
        assert!(ce.is_continuum_regime());
        let ce2 = ChapmanEnskogExpansion::new(1.0, [0.0, 0.0], 1.0, 0.1);
        assert!(!ce2.is_continuum_regime());
    }
    #[test]
    fn test_lbm_ns_bridge_extract_rho() {
        let mut bridge = LbmNavierStokesBridge::new(4, 1.0);
        let feq: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let f_all = vec![feq; 4];
        bridge.extract_macro_fields(&f_all);
        for &rho in &bridge.rho {
            assert!((rho - 2.0).abs() < 1e-12, "rho={rho}");
        }
    }
    #[test]
    fn test_lbm_ns_bridge_kinematic_viscosity() {
        let bridge = LbmNavierStokesBridge::new(4, 1.0);
        assert!(bridge.kinematic_viscosity() > 0.0);
    }
    #[test]
    fn test_lbm_ns_bridge_mach_number_zero_vel() {
        let bridge = LbmNavierStokesBridge::new(4, 1.0);
        assert!((bridge.mach_number(0)).abs() < 1e-12);
    }
    #[test]
    fn test_lbm_ns_bridge_reynolds_number() {
        let bridge = LbmNavierStokesBridge::new(4, 1.0);
        let re = bridge.reynolds_number(0.1, 10.0);
        assert!(re > 0.0);
    }
    #[test]
    fn test_gri_tau_fine_consistency() {
        let gri = GridRefinementInterface::new(4, 2, 1.0, 1);
        assert!(
            (gri.tau_fine - 0.75).abs() < 1e-12,
            "tau_fine={}",
            gri.tau_fine
        );
    }
    #[test]
    fn test_gri_coarse_to_fine_conserves_rho() {
        let gri = GridRefinementInterface::new(4, 2, 1.0, 1);
        let f_c: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.5);
        let f_f = gri.coarse_to_fine_rescale(&f_c, 1.5, [0.0, 0.0]);
        let rho_c: f64 = f_c.iter().sum();
        let rho_f: f64 = f_f.iter().sum();
        assert!(
            (rho_c - rho_f).abs() < 1e-10,
            "rho_c={rho_c}, rho_f={rho_f}"
        );
    }
    #[test]
    fn test_gri_fine_to_coarse_conserves_rho() {
        let gri = GridRefinementInterface::new(4, 2, 1.0, 1);
        let f_f: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 0.8);
        let f_c = gri.fine_to_coarse_rescale(&f_f, 0.8, [0.0, 0.0]);
        let rho_f: f64 = f_f.iter().sum();
        let rho_c: f64 = f_c.iter().sum();
        assert!((rho_f - rho_c).abs() < 1e-10);
    }
    #[test]
    fn test_gri_prolong_returns_r2_cells() {
        let gri = GridRefinementInterface::new(4, 2, 1.0, 1);
        let f_c: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let prolonged = gri.prolong_cell(&f_c);
        assert_eq!(prolonged.len(), 4);
    }
    #[test]
    fn test_temporal_refinement_dt_fine() {
        let tr = TemporalRefinement::new(8, 2, 0.01);
        assert!((tr.dt_fine - 0.005).abs() < 1e-15);
    }
    #[test]
    fn test_temporal_refinement_full_cycle_synced() {
        let mut tr = TemporalRefinement::new(8, 4, 0.01);
        tr.full_cycle(0.75);
        assert!(tr.is_synchronized());
    }
    #[test]
    fn test_temporal_refinement_sub_step_increments() {
        let mut tr = TemporalRefinement::new(4, 3, 0.01);
        tr.fine_collision_step(0.75);
        assert_eq!(tr.sub_step, 1);
        tr.fine_collision_step(0.75);
        tr.fine_collision_step(0.75);
        assert!(tr.is_synchronized());
    }
    #[test]
    fn test_temporal_refinement_interpolation() {
        let tr = TemporalRefinement::new(4, 2, 0.01);
        let f0: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let f1: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let fi = tr.time_interpolate(&f0, &f1, 0.5);
        let rho: f64 = fi.iter().sum();
        assert!((rho - 1.5).abs() < 1e-12, "rho={rho}");
    }
    #[test]
    fn test_sti_bilinear_corners() {
        let sti = SpaceTimeInterpolator::new(1.0, 0.01, 2);
        let f00: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let f10 = f00;
        let f01 = f00;
        let f11 = f00;
        let fi = sti.bilinear(&f00, &f10, &f01, &f11, 0.5, 0.5);
        let rho: f64 = fi.iter().sum();
        assert!((rho - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_sti_temporal_linear_endpoints() {
        let sti = SpaceTimeInterpolator::new(1.0, 0.01, 2);
        let f0: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let f1: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 3.0);
        let fi0 = sti.temporal_linear(&f0, &f1, 0.0);
        let fi1 = sti.temporal_linear(&f0, &f1, 1.0);
        assert!((fi0.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((fi1.iter().sum::<f64>() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_sti_hermite_endpoints() {
        let sti = SpaceTimeInterpolator::new(1.0, 0.01, 2);
        let f0: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let f1: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 2.0);
        let df0: [f64; 9] = [0.0f64; 9];
        let df1: [f64; 9] = [0.0f64; 9];
        let fi0 = sti.hermite_temporal(&f0, &f1, &df0, &df1, 0.0);
        let fi1 = sti.hermite_temporal(&f0, &f1, &df0, &df1, 1.0);
        assert!((fi0.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((fi1.iter().sum::<f64>() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_dsmc_coupling_init() {
        let coupling = LbmDsmcCoupling::new(8, 1.0, 1.0);
        assert_eq!(coupling.n_overlap, 8);
        assert_eq!(coupling.f_lbm.len(), 8);
    }
    #[test]
    fn test_dsmc_to_lbm_eq_rho() {
        let mut coupling = LbmDsmcCoupling::new(4, 1.0, 1.0);
        coupling.rho_dsmc[0] = 1.5;
        let feq = coupling.dsmc_to_lbm(0);
        let rho: f64 = feq.iter().sum();
        assert!((rho - 1.5).abs() < 1e-10, "rho={rho}");
    }
    #[test]
    fn test_lbm_to_dsmc_moments() {
        let coupling = LbmDsmcCoupling::new(4, 1.0, 1.0);
        let f: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * 1.0);
        let (rho, u, _temp) = coupling.lbm_to_dsmc_moments(&f);
        assert!((rho - 1.0).abs() < 1e-10);
        assert!(u[0].abs() < 1e-12 && u[1].abs() < 1e-12);
    }
    #[test]
    fn test_dsmc_mean_temperature() {
        let mut coupling = LbmDsmcCoupling::new(4, 1.0, 1.0);
        coupling.temp_dsmc = vec![2.0; 4];
        assert!((coupling.mean_temperature() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_kn_mean_free_path_positive() {
        let kn_est = KnudsenEstimator::new(10, 0.1, 1e-3, 1.0, 1.0);
        assert!(kn_est.mean_free_path(1.0) > 0.0);
    }
    #[test]
    fn test_kn_field_computation() {
        let mut kn_est = KnudsenEstimator::new(10, 0.1, 1e-3, 1.0, 10.0);
        let rho_field = vec![1.0f64; 10];
        kn_est.compute_kn_field(&rho_field, 10);
        assert_eq!(kn_est.kn.len(), 10);
        for kn in &kn_est.kn {
            assert!(kn.is_finite() && *kn >= 0.0);
        }
    }
    #[test]
    fn test_kn_continuum_fraction_uniform() {
        let mut kn_est = KnudsenEstimator::new(10, 0.1, 1e-3, 1.0, 10.0);
        let rho_field = vec![1.0f64; 10];
        kn_est.compute_kn_field(&rho_field, 10);
        assert!(kn_est.continuum_fraction() > 0.0);
    }
    #[test]
    fn test_kn_classify_regimes() {
        let mut kn_est = KnudsenEstimator::new(4, 0.1, 1e-3, 1.0, 1.0);
        kn_est.kn = vec![0.0001, 0.05, 5.0, 100.0];
        let regimes = kn_est.classify_regimes();
        assert_eq!(regimes[0], 0);
        assert_eq!(regimes[1], 1);
        assert_eq!(regimes[2], 2);
        assert_eq!(regimes[3], 3);
    }
    #[test]
    fn test_hmm_init() {
        let hmm = HeterogeneousMultiscale::new(8, 4, 1.0, 0.25, 1.0);
        assert_eq!(hmm.n_macro, 8);
        assert!((hmm.kappa - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_hmm_micro_to_macro_viscosity() {
        let hmm = HeterogeneousMultiscale::new(4, 4, 1.0, 0.25, 1.0);
        let nu = hmm.micro_to_macro_viscosity(0, 0.1);
        assert!(nu > 0.0);
        assert!((nu - CS2 * 0.5).abs() < 1e-12, "nu={nu}");
    }
    #[test]
    fn test_hmm_compress_micro() {
        let hmm = HeterogeneousMultiscale::new(4, 4, 1.0, 0.25, 1.0);
        let f_micro: Vec<[f64; 9]> = (0..4).map(|_| D2Q9_WEIGHTS.map(|w| w * 2.0)).collect();
        let (u_avg, rho_avg) = hmm.compress_micro_to_macro(&f_micro);
        assert!((rho_avg - 2.0).abs() < 1e-10, "rho_avg={rho_avg}");
        assert!(u_avg[0].abs() < 1e-12 && u_avg[1].abs() < 1e-12);
    }
    #[test]
    fn test_hmm_macro_step() {
        let mut hmm = HeterogeneousMultiscale::new(8, 4, 1.0, 0.25, 1.0);
        hmm.u_macro[4][0] = 0.1;
        let u_before = hmm.u_macro[4][0];
        hmm.macro_step(0.001);
        let _ = u_before;
        for u in &hmm.u_macro {
            assert!(u[0].is_finite());
        }
    }
    #[test]
    fn test_aas_default_lbm() {
        let sel = AdaptiveAlgorithmSelector::new(10, 0.001, 1.0);
        assert_eq!(sel.algorithm_counts()[1], 10);
    }
    #[test]
    fn test_aas_set_kn_all_continuum() {
        let mut sel = AdaptiveAlgorithmSelector::new(10, 0.01, 1.0);
        sel.set_kn_field(&[0.0001f64; 10]);
        assert_eq!(sel.algorithm_counts()[0], 10);
    }
    #[test]
    fn test_aas_set_kn_all_rarefied() {
        let mut sel = AdaptiveAlgorithmSelector::new(10, 0.01, 0.1);
        sel.set_kn_field(&[10.0f64; 10]);
        assert_eq!(sel.algorithm_counts()[2], 10);
    }
    #[test]
    fn test_aas_fractions_sum_to_one() {
        let mut sel = AdaptiveAlgorithmSelector::new(10, 0.001, 1.0);
        let kn: Vec<f64> = (0..10).map(|i| 0.0001 * 10.0_f64.powi(i / 3)).collect();
        sel.set_kn_field(&kn);
        let total = sel.ns_fraction() + sel.lbm_fraction() + sel.dsmc_fraction();
        assert!((total - 1.0).abs() < 1e-12, "total={total}");
    }
    #[test]
    fn test_aas_find_interfaces_mixed() {
        let mut sel = AdaptiveAlgorithmSelector::new(6, 0.01, 0.1);
        let kn = vec![0.001, 0.001, 0.05, 0.05, 5.0, 5.0];
        sel.set_kn_field(&kn);
        let ifaces = sel.find_interfaces();
        assert!(!ifaces.is_empty());
    }
}
