//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::FoamDict;

/// Estimate mesh Reynolds number: `Re = U * L / nu`.
pub fn mesh_reynolds_number(u_ref: f64, length_scale: f64, nu: f64) -> f64 {
    if nu < 1e-30 {
        return f64::INFINITY;
    }
    u_ref * length_scale / nu
}
/// Courant-Friedrichs-Lewy (CFL) number for an explicit time scheme.
///
/// `CFL = U * dt / dx`
pub fn cfl_number(u: f64, dt: f64, dx: f64) -> f64 {
    if dx < 1e-30 {
        return f64::INFINITY;
    }
    u * dt / dx
}
/// Compute the maximum allowable time step for explicit time integration
/// given a CFL limit.
///
/// `dt_max = CFL_limit * dx / U`
pub fn max_dt_from_cfl(cfl_limit: f64, dx: f64, u: f64) -> f64 {
    if u < 1e-30 {
        return f64::INFINITY;
    }
    cfl_limit * dx / u
}
/// Merge two `FoamDict` objects, with entries from `other` overwriting same-key
/// entries in `base`.
pub fn merge_foam_dicts(base: &FoamDict, other: &FoamDict) -> FoamDict {
    let mut result = base.clone();
    for (key, val) in &other.entries {
        result.entries.retain(|(k, _)| k != key);
        result.entries.push((key.clone(), val.clone()));
    }
    result
}
#[cfg(test)]
mod tests_openfoam_mesh {
    use crate::openfoam::*;
    #[test]
    fn test_face_area_unit_square() {
        let verts = vec![
            [0.0, 0.0, 0.0_f64],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let area = face_area(&verts);
        assert!((area - 1.0).abs() < 1e-10, "area={}", area);
    }
    #[test]
    fn test_face_area_triangle() {
        let verts = vec![[0.0, 0.0, 0.0_f64], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let area = face_area(&verts);
        assert!((area - 6.0).abs() < 1e-10, "area={}", area);
    }
    #[test]
    fn test_face_area_degenerate() {
        let verts = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0]];
        assert!((face_area(&verts)).abs() < 1e-10);
    }
    #[test]
    fn test_face_centroid_quad() {
        let verts = vec![
            [0.0, 0.0, 0.0_f64],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
        ];
        let c = face_centroid(&verts);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
        assert!((c[2]).abs() < 1e-10);
    }
    #[test]
    fn test_face_normal_xy_plane() {
        let verts = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let n = face_normal(&verts);
        assert!((n[2].abs() - 1.0).abs() < 1e-10, "n={:?}", n);
    }
    #[test]
    fn test_non_orthogonality_zero_for_aligned() {
        let face_norm = [1.0, 0.0, 0.0];
        let owner_centre = [0.0, 0.0, 0.0];
        let face_centre = [1.0, 0.0, 0.0];
        let noo = non_orthogonality(face_norm, owner_centre, face_centre);
        assert!(noo.abs() < 1e-8, "noo={}", noo);
    }
    #[test]
    fn test_non_orthogonality_ninety_degrees() {
        let face_norm = [0.0, 1.0, 0.0];
        let owner_centre = [0.0, 0.0, 0.0];
        let face_centre = [1.0, 0.0, 0.0];
        let noo = non_orthogonality(face_norm, owner_centre, face_centre);
        assert!((noo - 90.0).abs() < 1e-6, "noo={}", noo);
    }
    #[test]
    fn test_average_cell_volume_unit_cube() {
        let vol = average_cell_volume(1.0, 1.0, 1.0, 1, 1, 1);
        assert!((vol - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_average_cell_volume_subdivided() {
        let vol = average_cell_volume(2.0, 2.0, 2.0, 2, 2, 2);
        assert!((vol - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quad_aspect_ratio_square() {
        let ar = quad_face_aspect_ratio(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!((ar - 1.0).abs() < 1e-10, "ar={}", ar);
    }
    #[test]
    fn test_quad_aspect_ratio_rectangle() {
        let ar = quad_face_aspect_ratio(
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!((ar - 4.0).abs() < 1e-10, "ar={}", ar);
    }
    #[test]
    fn test_k_epsilon_initial_values() {
        let (k, eps) = k_epsilon_initial(1.0, 0.05, 0.1);
        assert!((k - 0.00375).abs() < 1e-8, "k={}", k);
        assert!(eps > 0.0, "epsilon must be positive");
    }
    #[test]
    fn test_k_omega_initial_values() {
        let (k, omega) = k_omega_initial(1.0, 0.05, 0.1);
        assert!(k > 0.0);
        assert!(omega > 0.0);
    }
    #[test]
    fn test_turbulent_viscosity_ke() {
        let nu_t = turbulent_viscosity_ke(0.00375, 0.1);
        let expected = 0.09 * 0.00375_f64.powi(2) / 0.1;
        assert!((nu_t - expected).abs() < 1e-12, "nu_t={}", nu_t);
    }
    #[test]
    fn test_turbulent_viscosity_ko() {
        let nu_t = turbulent_viscosity_ko(1.0, 2.0);
        assert!((nu_t - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_find_cell_containing_point_center() {
        let cell = find_cell_containing_point([1.5, 1.5, 1.5], 2.0, 2.0, 2.0, 2, 2, 2);
        assert!(cell.is_some(), "should find a cell");
    }
    #[test]
    fn test_find_cell_containing_point_outside() {
        let cell = find_cell_containing_point([3.0, 0.0, 0.0], 2.0, 2.0, 2.0, 2, 2, 2);
        assert!(cell.is_none());
    }
    #[test]
    fn test_probe_scalar_basic() {
        let vals = vec![10.0, 20.0, 30.0];
        assert_eq!(probe_scalar(&vals, 1), Some(20.0));
        assert_eq!(probe_scalar(&vals, 5), None);
    }
    #[test]
    fn test_probe_vector_basic() {
        let vals: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let v = probe_vector(&vals, 0).unwrap();
        assert!((v[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cfl_number_basic() {
        let cfl = cfl_number(1.0, 0.01, 0.1);
        assert!((cfl - 0.1).abs() < 1e-10, "cfl={}", cfl);
    }
    #[test]
    fn test_max_dt_from_cfl() {
        let dt = max_dt_from_cfl(0.5, 0.1, 1.0);
        assert!((dt - 0.05).abs() < 1e-10, "dt={}", dt);
    }
    #[test]
    fn test_mesh_reynolds_number() {
        let re = mesh_reynolds_number(1.0, 1.0, 1e-6);
        assert!((re - 1e6).abs() < 1.0, "re={}", re);
    }
    #[test]
    fn test_merge_foam_dicts_basic() {
        let mut base = FoamDict::new();
        base.insert("a", FoamValue::Scalar(1.0));
        base.insert("b", FoamValue::Scalar(2.0));
        let mut other = FoamDict::new();
        other.insert("b", FoamValue::Scalar(99.0));
        other.insert("c", FoamValue::Scalar(3.0));
        let merged = merge_foam_dicts(&base, &other);
        assert_eq!(merged.get_scalar("a"), Some(1.0));
        assert_eq!(merged.get_scalar("b"), Some(99.0));
        assert_eq!(merged.get_scalar("c"), Some(3.0));
    }
    #[test]
    fn test_merge_foam_dicts_empty_other() {
        let mut base = FoamDict::new();
        base.insert("x", FoamValue::Scalar(5.0));
        let other = FoamDict::new();
        let merged = merge_foam_dicts(&base, &other);
        assert_eq!(merged.get_scalar("x"), Some(5.0));
    }
    #[test]
    fn test_hex_cell_volume_unit_cube() {
        let corners: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let vol = hex_cell_volume(&corners);
        assert!((vol - 1.0).abs() < 1e-8, "vol={}", vol);
    }
    #[test]
    fn test_cell_centroid_single_face() {
        let points: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let faces: Vec<Vec<usize>> = vec![vec![0, 1, 2, 3]];
        let c = cell_centroid_from_faces(&[0], &faces, &points);
        assert!((c[0] - 0.5).abs() < 1e-10);
        assert!((c[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_green_gauss_gradient_constant_field() {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let n_cells = mesh.n_cells;
        let values = vec![1.0_f64; n_cells];
        let cell_centres: Vec<[f64; 3]> = (0..n_cells)
            .map(|i| {
                let ix = i % 2;
                let iy = (i / 2) % 2;
                let iz = i / 4;
                [
                    (ix as f64 + 0.5) * 0.5,
                    (iy as f64 + 0.5) * 0.5,
                    (iz as f64 + 0.5) * 0.5,
                ]
            })
            .collect();
        let grad = green_gauss_gradient(
            &values,
            &cell_centres,
            &mesh.faces,
            &mesh.points,
            &mesh.owner,
            &mesh.neighbour,
            n_cells,
        );
        assert_eq!(grad.len(), n_cells);
        for g in &grad {
            assert!(
                g[0].is_finite() && g[1].is_finite() && g[2].is_finite(),
                "gradient contains non-finite value: {:?}",
                g
            );
        }
    }
}
