//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_extended {
    use crate::csg::CsgSmoothUnion;
    use crate::csg::CsgTree;
    use crate::csg::MarchingCell;
    use crate::csg::SdfBox;
    use crate::csg::SdfCappedCylinder;
    use crate::csg::SdfCapsule;
    use crate::csg::SdfCone;
    use crate::csg::SdfRoundedBox;
    use crate::csg::SdfSphere;
    use crate::csg::*;
    #[test]
    fn capped_cylinder_inside_is_negative() {
        let c = SdfCappedCylinder::new([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(c.sdf([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn capped_cylinder_outside_is_positive() {
        let c = SdfCappedCylinder::new([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(c.sdf([3.0, 0.0, 0.0]) > 0.0);
        assert!(c.sdf([0.0, 5.0, 0.0]) > 0.0);
    }
    #[test]
    fn capped_cylinder_surface_near_zero() {
        let c = SdfCappedCylinder::new([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(c.sdf([1.0, 0.0, 0.0]).abs() < 1e-10);
        assert!(c.sdf([0.0, 2.0, 0.0]).abs() < 1e-10);
    }
    #[test]
    fn capsule_inside_is_negative() {
        let c = SdfCapsule::new([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(c.sdf([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn capsule_on_surface_near_zero() {
        let c = SdfCapsule::new([0.0, 0.0, 0.0], 1.0, 1.0);
        assert!(c.sdf([1.0, 0.0, 0.0]).abs() < 1e-10);
    }
    #[test]
    fn capsule_outside_is_positive() {
        let c = SdfCapsule::new([0.0, 0.0, 0.0], 0.5, 1.0);
        assert!(c.sdf([5.0, 0.0, 0.0]) > 0.0);
    }
    #[test]
    fn rounded_box_inside_is_negative() {
        let b = SdfRoundedBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.1);
        assert!(b.sdf([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn rounded_box_outside_is_positive() {
        let b = SdfRoundedBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.1);
        assert!(b.sdf([5.0, 0.0, 0.0]) > 0.0);
    }
    #[test]
    fn rounded_box_smaller_than_sharp_box() {
        let sharp = SdfBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let rounded = SdfRoundedBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.2);
        let p = [0.0, 0.0, 0.0];
        assert!(
            rounded.sdf(p) < sharp.sdf(p) + 1e-9,
            "rounded box sdf ({}) should be <= sharp box sdf ({}) at center",
            rounded.sdf(p),
            sharp.sdf(p)
        );
    }
    fn sphere_bounds_fn(s: &dyn ImplicitSurface) -> ([f64; 3], [f64; 3]) {
        let _ = s;
        ([-2.0; 3], [2.0; 3])
    }
    #[test]
    fn csg_tree_aabb_union_expands() {
        let a = CsgTree::leaf(SdfSphere::new([-3.0, 0.0, 0.0], 1.0));
        let b = CsgTree::leaf(SdfSphere::new([3.0, 0.0, 0.0], 1.0));
        let u = CsgTree::union(a, b);
        let (mn, mx) = csg_tree_aabb(&u, &sphere_bounds_fn);
        for k in 0..3 {
            assert!((mn[k] + 2.0).abs() < 1e-9, "mn[{k}]={}", mn[k]);
            assert!((mx[k] - 2.0).abs() < 1e-9, "mx[{k}]={}", mx[k]);
        }
    }
    #[test]
    fn csg_tree_aabb_leaf_returns_bounds() {
        let tree = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let (mn, mx) = csg_tree_aabb(&tree, &sphere_bounds_fn);
        assert_eq!(mn, [-2.0; 3]);
        assert_eq!(mx, [2.0; 3]);
    }
    #[test]
    fn csg_tree_aabb_difference_uses_left() {
        let a = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 3.0));
        let b = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.5));
        let d = CsgTree::difference(a, b);
        let (mn, mx) = csg_tree_aabb(&d, &sphere_bounds_fn);
        assert_eq!(mn, [-2.0; 3]);
        assert_eq!(mx, [2.0; 3]);
    }
    #[test]
    fn marching_cell_has_surface_mixed_signs() {
        let mut vals = [-1.0f64; 8];
        vals[7] = 1.0;
        let cell = MarchingCell {
            min: [0.0; 3],
            step: 1.0,
            corner_values: vals,
        };
        assert!(cell.has_surface());
    }
    #[test]
    fn marching_cell_no_surface_all_negative() {
        let cell = MarchingCell {
            min: [0.0; 3],
            step: 1.0,
            corner_values: [-1.0; 8],
        };
        assert!(!cell.has_surface());
    }
    #[test]
    fn marching_cell_no_surface_all_positive() {
        let cell = MarchingCell {
            min: [0.0; 3],
            step: 1.0,
            corner_values: [1.0; 8],
        };
        assert!(!cell.has_surface());
    }
    #[test]
    fn marching_cell_extract_vertices_nonempty_on_crossing() {
        let mut vals = [-1.0f64; 8];
        vals[0] = 1.0;
        let cell = MarchingCell {
            min: [0.0; 3],
            step: 1.0,
            corner_values: vals,
        };
        let verts = cell.extract_vertices();
        assert!(!verts.is_empty(), "should have crossing vertices");
    }
    #[test]
    fn marching_cubes_sample_sphere_cells_nonempty() {
        let s = SdfSphere::new([0.0, 0.0, 0.0], 1.0);
        let cells = marching_cubes_sample(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], 8);
        assert!(!cells.is_empty(), "should find surface cells for sphere");
    }
    #[test]
    fn marching_cubes_vertices_nonempty_for_sphere() {
        let s = SdfSphere::new([0.0, 0.0, 0.0], 1.0);
        let verts = marching_cubes_vertices(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], 8);
        assert!(!verts.is_empty(), "should find vertices on sphere surface");
    }
    #[test]
    fn marching_cubes_all_inside_gives_no_surface() {
        let s = SdfSphere::new([0.0, 0.0, 0.0], 100.0);
        let cells = marching_cubes_sample(&s, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], 4);
        assert!(
            cells.is_empty(),
            "all-inside grid should have no surface cells"
        );
    }
    #[test]
    fn cone_inside_is_negative() {
        let c = SdfCone::new([0.0, 0.0, 0.0], std::f64::consts::PI / 4.0, 2.0);
        let d = c.sdf([0.0, 1.0, 0.0]);
        assert!(
            d < 0.0,
            "point on axis inside cone should be negative, got {d}"
        );
    }
    #[test]
    fn cone_above_base_is_positive() {
        let c = SdfCone::new([0.0, 0.0, 0.0], std::f64::consts::PI / 4.0, 1.0);
        let d = c.sdf([0.0, 5.0, 0.0]);
        assert!(d > 0.0, "point above cone base should be positive, got {d}");
    }
    #[test]
    fn csg_union_sphere_and_cylinder() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let c = CsgTree::leaf(SdfCappedCylinder::new([0.0, 0.0, 0.0], 0.5, 3.0));
        let u = CsgTree::union(s, c);
        assert!(
            csg_sdf(&u, [0.0, 2.5, 0.0]) < 0.0,
            "inside cylinder should be inside union"
        );
        assert!(
            csg_sdf(&u, [5.0, 0.0, 0.0]) > 0.0,
            "outside both should be outside union"
        );
    }
    #[test]
    fn csg_intersection_sphere_and_box() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 2.0));
        let b = CsgTree::leaf(SdfBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]));
        let i = CsgTree::intersection(s, b);
        assert!(
            csg_sdf(&i, [0.0, 0.0, 0.0]) < 0.0,
            "center is in intersection"
        );
        assert!(
            csg_sdf(&i, [1.5, 0.0, 0.0]) > 0.0,
            "outside box is outside intersection"
        );
    }
    #[test]
    fn smooth_union_vs_sharp_at_midpoint() {
        let a = Box::new(SdfSphere::new([-1.0, 0.0, 0.0], 1.0));
        let b = Box::new(SdfSphere::new([1.0, 0.0, 0.0], 1.0));
        let su = CsgSmoothUnion::new(a, b, 0.5);
        let d = su.sdf([0.0, 0.0, 0.0]);
        assert!(
            d <= 0.5,
            "smooth union at midpoint should be blended: d={d}"
        );
    }
    #[test]
    fn csg_surface_area_sphere_positive() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let area = csg_surface_area(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], 10);
        assert!(area > 0.0, "sphere surface area should be > 0, got {area}");
    }
    #[test]
    fn csg_surface_area_increases_with_resolution() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let a_low = csg_surface_area(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], 8);
        let a_high = csg_surface_area(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5], 16);
        assert!(
            a_low > 0.0 && a_high > 0.0,
            "areas should be positive: low={a_low} high={a_high}"
        );
    }
    #[test]
    fn csg_node_surface_area_sphere_approx() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let area = csg_node_surface_area(&s, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5]);
        assert!(
            area > 5.0 && area < 30.0,
            "sphere SA estimate should be in [5,30], got {area}"
        );
    }
    #[test]
    fn csg_offset_sdf_outward_expands_shape() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let p = [1.3, 0.0, 0.0];
        let d_original = csg_sdf(&s, p);
        let d_offset = csg_offset_sdf(&s, p, 0.5);
        assert!(d_original > 0.0, "original should be outside: {d_original}");
        assert!(
            d_offset < 0.0,
            "offset=0.5 outward should include point: {d_offset}"
        );
    }
    #[test]
    fn csg_offset_sdf_inward_shrinks_shape() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let p = [0.0, 0.0, 0.0];
        let d_offset = csg_offset_sdf(&s, p, -0.5);
        assert!(
            d_offset < 0.0,
            "inward offset should keep origin inside: {d_offset}"
        );
    }
    #[test]
    fn csg_offset_contains_outward() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let p = [1.1, 0.0, 0.0];
        assert!(
            csg_offset_contains(&s, p, 0.2),
            "point at distance 0.1 should be inside offset 0.2 surface"
        );
    }
    #[test]
    fn csg_simplify_leaf_always_true() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        assert!(
            csg_tree_simplify_check(&s, [0.0, 0.0, 0.0]),
            "leaf node is always simplified"
        );
    }
    #[test]
    fn csg_simplify_union_inside_one_branch() {
        let a = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 2.0));
        let b = CsgTree::leaf(SdfSphere::new([10.0, 0.0, 0.0], 1.0));
        let u = CsgTree::union(a, b);
        assert!(
            csg_tree_simplify_check(&u, [0.0, 0.0, 0.0]),
            "union: probe inside left branch → simplified"
        );
    }
    #[test]
    fn csg_offset_gradient_is_unit_length() {
        let s = CsgTree::leaf(SdfSphere::new([0.0, 0.0, 0.0], 1.0));
        let p = [1.0, 0.0, 0.0];
        let g = csg_offset_gradient(&s, p);
        let len = (g[0].powi(2) + g[1].powi(2) + g[2].powi(2)).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-3,
            "offset gradient should be unit length, got {len}"
        );
    }
}
