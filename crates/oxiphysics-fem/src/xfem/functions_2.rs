//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Reinitialise a level set field to a signed distance function by a
/// single-pass correction.
///
/// This is a simplified iterative regularisation: each nodal value is updated
/// to `±min(|phi_i|, min_neighbor + h)` where `h` is the mesh spacing.
/// For accurate results, a full fast-marching method should be used.
///
/// # Arguments
/// * `phi` – mutable reference to the nodal level set values
/// * `connectivity` – slice of edge pairs (node_a, node_b)
/// * `h` – characteristic element size
pub fn reinitialise_level_set(phi: &mut [f64], connectivity: &[(usize, usize)], h: f64) {
    let n = phi.len();
    let mut new_phi = phi.to_vec();
    for &(a, b) in connectivity {
        if a >= n || b >= n {
            continue;
        }
        let da = phi[a];
        let db = phi[b];
        let sign_a = if da >= 0.0 { 1.0 } else { -1.0 };
        let sign_b = if db >= 0.0 { 1.0 } else { -1.0 };
        let proposed_a = sign_a * (db.abs() + h).min(da.abs());
        let proposed_b = sign_b * (da.abs() + h).min(db.abs());
        if proposed_a.abs() < new_phi[a].abs() {
            new_phi[a] = proposed_a;
        }
        if proposed_b.abs() < new_phi[b].abs() {
            new_phi[b] = proposed_b;
        }
    }
    phi.copy_from_slice(&new_phi);
}
#[cfg(test)]
mod tests_xfem_expanded {
    use super::*;
    use crate::xfem::*;
    #[test]
    fn test_sub_triangle_area() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let tri = SubTriangle::new(verts, 1);
        assert!((tri.area() - 0.5).abs() < 1e-12, "area = {}", tri.area());
    }
    #[test]
    fn test_sub_triangle_centroid() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let tri = SubTriangle::new(verts, 1);
        let c = tri.centroid();
        assert!((c[0] - 1.0 / 3.0).abs() < 1e-12, "cx = {}", c[0]);
        assert!((c[1] - 1.0 / 3.0).abs() < 1e-12, "cy = {}", c[1]);
    }
    #[test]
    fn test_sub_triangle_gauss_points_weight_sum() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let tri = SubTriangle::new(verts, 1);
        let gp = tri.gauss_points();
        let total_w: f64 = gp.iter().map(|&(_, _, w)| w).sum();
        assert!((total_w - 0.5).abs() < 1e-12, "weight sum = {total_w}");
    }
    #[test]
    fn test_sub_triangulate_no_cut() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let phi = [1.0, 1.0, 1.0];
        let tris = sub_triangulate_triangle(&verts, &phi);
        assert_eq!(tris.len(), 1, "no cut: should return 1 sub-triangle");
        assert_eq!(tris[0].side, 1);
    }
    #[test]
    fn test_sub_triangulate_cut_area_conservation() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let phi = [-1.0, 1.0, 0.0];
        let tris = sub_triangulate_triangle(&verts, &phi);
        let original_area = {
            let a = [verts[1][0] - verts[0][0], verts[1][1] - verts[0][1]];
            let b = [verts[2][0] - verts[0][0], verts[2][1] - verts[0][1]];
            0.5 * (a[0] * b[1] - a[1] * b[0]).abs()
        };
        let sub_area: f64 = tris.iter().map(|t| t.area()).sum();
        assert!(
            (sub_area - original_area).abs() < 1e-10,
            "area conservation: {sub_area} vs {original_area}"
        );
    }
    #[test]
    fn test_total_xfem_dofs() {
        let n = total_xfem_dofs(10, 3, 4, 2);
        assert_eq!(n, 66);
    }
    #[test]
    fn test_build_enriched_dof_map_none_for_unenriched() {
        let map = build_enriched_dof_map(5, 3, &[1, 2], &[]);
        assert!(map[0].is_none(), "node 0 not enriched");
        assert!(map[3].is_none(), "node 3 not enriched");
    }
    #[test]
    fn test_build_enriched_dof_map_heaviside_offset() {
        let map = build_enriched_dof_map(4, 3, &[0], &[]);
        assert_eq!(
            map[0],
            Some(12),
            "Heaviside offset for node 0 = {}",
            map[0].unwrap_or(0)
        );
    }
    #[test]
    fn test_locate_crack_tip_basic() {
        let crack = Crack::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let nodes = vec![[0.5, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let ls = LevelSet::from_crack(&nodes, &crack);
        let result = locate_crack_tip_from_level_set(&nodes, &ls);
        assert!(result.is_some());
        let (idx, _) = result.unwrap();
        assert!(idx < 2, "tip should be near crack tip, got node {idx}");
    }
    #[test]
    fn test_cmod_basic() {
        let u_upper = [0.0, 1e-4, 0.0];
        let u_lower = [0.0, -1e-4, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let cmod = crack_mouth_opening_displacement(u_upper, u_lower, normal);
        assert!((cmod - 2e-4).abs() < 1e-15, "CMOD = {cmod}");
    }
    #[test]
    fn test_ctod_irwin_positive() {
        let ctod = ctod_irwin(30e6, 200e9, 0.3, 1e-3);
        assert!(ctod > 0.0 && ctod.is_finite(), "CTOD = {ctod}");
    }
    #[test]
    fn test_ctod_irwin_scales_with_k1() {
        let c1 = ctod_irwin(30e6, 200e9, 0.3, 1e-3);
        let c2 = ctod_irwin(60e6, 200e9, 0.3, 1e-3);
        assert!((c2 / c1 - 2.0).abs() < 1e-10, "CTOD should scale with K_I");
    }
    #[test]
    fn test_interaction_integral_trapezoidal_constant() {
        let n = 100;
        let vals = vec![1.0; n];
        let r = 1.0;
        let i_val = interaction_integral_trapezoidal(&vals, r);
        assert!(
            (i_val - 2.0 * std::f64::consts::PI).abs() < 0.01,
            "integral = {i_val}"
        );
    }
    #[test]
    fn test_reinitialise_level_set_no_change_far_from_interface() {
        let mut phi = vec![10.0, 9.0, 8.0];
        let connectivity = vec![(0usize, 1usize), (1, 2)];
        let phi_before = phi.clone();
        reinitialise_level_set(&mut phi, &connectivity, 1.0);
        for (i, (before, after)) in phi_before.iter().zip(phi.iter()).enumerate() {
            assert!(
                before.signum() == after.signum(),
                "sign flipped at node {i}"
            );
        }
    }
    #[test]
    fn test_sub_triangle_gauss_points_inside() {
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let tri = SubTriangle::new(verts, 1);
        for (xi, eta, _w) in &tri.gauss_points() {
            assert!(*xi >= 0.0 - 1e-10, "xi = {xi}");
            assert!(*eta >= 0.0 - 1e-10, "eta = {eta}");
            assert!(*xi + *eta <= 1.0 + 1e-10, "xi+eta = {}", xi + eta);
        }
    }
}
