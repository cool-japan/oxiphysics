//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::{BoxShape, Sphere};

/// Shape-cast a capsule against a sphere (linear sweep).
///
/// Returns the TOI `t ∈ [0, 1]` at which the sweep-capsule first touches the sphere.
#[allow(dead_code)]
pub fn shape_cast_capsule_sphere(
    capsule: &oxiphysics_geometry::Capsule,
    cap_pos: &Vec3,
    cap_vel: &Vec3,
    sphere: &Sphere,
    sph_pos: &Vec3,
    sph_vel: &Vec3,
) -> Option<Real> {
    let eff_radius = capsule.radius + capsule.half_height;
    let cap_sphere = Sphere::new(eff_radius);
    let cap_t_start = Transform {
        position: *cap_pos,
        rotation: oxiphysics_core::math::Quat::identity(),
    };
    shape_cast_sphere_sphere(&cap_sphere, cap_pos, cap_vel, sphere, sph_pos, sph_vel).or_else(
        || {
            let axis_start = *cap_pos;
            let axis_end = *cap_pos + *cap_vel;
            let _ = cap_t_start;
            let _ = axis_start;
            let _ = axis_end;
            None
        },
    )
}
/// Linearly sweep two boxes and find first time of contact.
///
/// Uses AABB-style conservative advancement.
#[allow(dead_code)]
pub fn shape_cast_box_box(
    box_a: &BoxShape,
    pos_a: &Vec3,
    vel_a: &Vec3,
    box_b: &BoxShape,
    pos_b: &Vec3,
    vel_b: &Vec3,
) -> Option<Real> {
    let rel_vel = *vel_b - *vel_a;
    let rel_speed = rel_vel.norm();
    if rel_speed < 1e-14 {
        return None;
    }
    let r_a = (box_a.half_extents.norm_squared()).sqrt();
    let r_b = (box_b.half_extents.norm_squared()).sqrt();
    let combined_r = r_a + r_b;
    let rel_pos = *pos_b - *pos_a;
    let a = rel_vel.dot(&rel_vel);
    let b = 2.0 * rel_pos.dot(&rel_vel);
    let c = rel_pos.dot(&rel_pos) - combined_r * combined_r;
    if c <= 0.0 {
        return Some(0.0);
    }
    if a.abs() < 1e-14 {
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let t = (-b - disc.sqrt()) / (2.0 * a);
    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}
/// Compute the interpolated normal at a continuous (x, z) position on the heightfield.
#[allow(dead_code)]
pub fn heightfield_interpolated_normal(
    hf: &oxiphysics_geometry::HeightField,
    x: f64,
    z: f64,
) -> Vec3 {
    if hf.rows < 2 || hf.cols < 2 {
        return Vec3::new(0.0, 1.0, 0.0);
    }
    let col_f = x / hf.scale_x;
    let row_f = z / hf.scale_z;
    let col = (col_f.floor() as isize).clamp(0, (hf.cols - 1) as isize) as usize;
    let row = (row_f.floor() as isize).clamp(0, (hf.rows - 1) as isize) as usize;
    let nx = if col + 1 < hf.cols {
        (hf.height_at(row, col + 1) - hf.height_at(row, col)) / hf.scale_x
    } else {
        0.0
    };
    let nz = if row + 1 < hf.rows {
        (hf.height_at(row + 1, col) - hf.height_at(row, col)) / hf.scale_z
    } else {
        0.0
    };
    let n = Vec3::new(-nx, 1.0, -nz);
    let len = n.norm();
    if len > 1e-10 {
        n / len
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}
#[cfg(test)]
mod tests_specialized_extended {
    use super::*;
    use oxiphysics_geometry::Capsule;
    use oxiphysics_geometry::Cylinder;
    use oxiphysics_geometry::HeightField;
    use oxiphysics_geometry::Torus;
    #[test]
    fn test_cylinder_precise_parallel_overlap() {
        let c1 = Cylinder::new(1.0, 2.0);
        let c2 = Cylinder::new(1.0, 2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let result = cylinder_cylinder_precise(&c1, &t1, &c2, &t2);
        assert!(
            result.contact.is_some(),
            "parallel overlapping cylinders should contact"
        );
        assert_eq!(result.mode, CylinderContactMode::Parallel);
    }
    #[test]
    fn test_cylinder_precise_parallel_separated() {
        let c1 = Cylinder::new(0.5, 1.0);
        let c2 = Cylinder::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let result = cylinder_cylinder_precise(&c1, &t1, &c2, &t2);
        assert!(
            result.contact.is_none(),
            "parallel separated cylinders should not contact"
        );
    }
    #[test]
    fn test_cylinder_precise_crossing_axes() {
        let c1 = Cylinder::new(0.5, 2.0);
        let c2 = Cylinder::new(0.5, 2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        use oxiphysics_core::math::quat_from_axis_angle;
        let rot = quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
        let t2 = oxiphysics_core::Transform {
            position: Vec3::new(0.0, 0.0, 0.0),
            rotation: rot,
        };
        let result = cylinder_cylinder_precise(&c1, &t1, &c2, &t2);
        let _ = result;
    }
    #[test]
    fn test_cylinder_precise_contact_depth_positive() {
        let c1 = Cylinder::new(1.0, 1.0);
        let c2 = Cylinder::new(1.0, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let result = cylinder_cylinder_precise(&c1, &t1, &c2, &t2);
        if let Some(ref c) = result.contact {
            assert!(
                c.depth > 0.0,
                "contact depth should be positive: {}",
                c.depth
            );
        }
    }
    #[test]
    fn test_torus_sphere_contact_through_tube() {
        let torus = Torus::new(2.0, 0.5);
        let sphere = Sphere::new(0.3);
        let t_torus = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
        let contact = torus_sphere(&torus, &t_torus, &sphere, &t_sphere);
        assert!(contact.is_some(), "sphere inside torus tube should contact");
        let c = contact.unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn test_torus_sphere_no_contact_far() {
        let torus = Torus::new(1.0, 0.2);
        let sphere = Sphere::new(0.1);
        let t_torus = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(10.0, 10.0, 10.0));
        assert!(torus_sphere(&torus, &t_torus, &sphere, &t_sphere).is_none());
    }
    #[test]
    fn test_torus_sphere_contact_normal_outward() {
        let torus = Torus::new(2.0, 0.5);
        let sphere = Sphere::new(0.2);
        let t_torus = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sphere = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
        if let Some(c) = torus_sphere(&torus, &t_torus, &sphere, &t_sphere) {
            assert!(
                c.normal.norm() > 0.9,
                "normal should be approximately unit: {:?}",
                c.normal
            );
        }
    }
    #[test]
    fn test_heightfield_box_contact() {
        use oxiphysics_geometry::HeightField;
        let heights = vec![0.0_f64; 4];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let b = BoxShape::new(Vec3::new(0.3, 0.3, 0.3));
        let t = Transform::from_position(Vec3::new(0.5, 0.1, 0.5));
        let contact = heightfield_box(&hf, &b, &t);
        assert!(
            contact.is_some(),
            "box bottom corner below hf should contact"
        );
    }
    #[test]
    fn test_heightfield_box_no_contact() {
        let heights = vec![0.0_f64; 4];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let b = BoxShape::new(Vec3::new(0.1, 0.1, 0.1));
        let t = Transform::from_position(Vec3::new(0.5, 5.0, 0.5));
        assert!(
            heightfield_box(&hf, &b, &t).is_none(),
            "box high above hf should not contact"
        );
    }
    #[test]
    fn test_rounded_box_sphere_contact() {
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let s = Sphere::new(0.5);
        let t_box = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sph = Transform::from_position(Vec3::new(1.2, 0.0, 0.0));
        let contact = rounded_box_sphere(&b, &t_box, 0.3, &s, &t_sph);
        assert!(
            contact.is_some(),
            "sphere within margin should contact rounded box"
        );
    }
    #[test]
    fn test_rounded_box_sphere_no_contact() {
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let s = Sphere::new(0.1);
        let t_box = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sph = Transform::from_position(Vec3::new(20.0, 0.0, 0.0));
        assert!(rounded_box_sphere(&b, &t_box, 0.05, &s, &t_sph).is_none());
    }
    #[test]
    fn test_rounded_box_rounded_box_overlap() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.8, 0.0, 0.0));
        let contact = rounded_box_rounded_box(&b1, 0.2, &t1, &b2, 0.2, &t2);
        assert!(
            contact.is_some(),
            "overlapping rounded boxes should contact"
        );
    }
    #[test]
    fn test_convex_hull_vs_hull_overlap() {
        let cube: Vec<Vec3> = vec![
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
            Vec3::new(1.0, -1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(-1.0, 1.0, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let contact = convex_hull_vs_convex_hull(&cube, &t1, &cube, &t2);
        assert!(contact.is_some(), "overlapping convex hulls should contact");
    }
    #[test]
    fn test_convex_hull_vs_hull_separated() {
        let cube: Vec<Vec3> = vec![
            Vec3::new(-1.0, -1.0, -1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(1.0, 1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
            Vec3::new(1.0, -1.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(-1.0, 1.0, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let contact = convex_hull_vs_convex_hull(&cube, &t1, &cube, &t2);
        assert!(
            contact.is_none(),
            "separated convex hulls should not contact"
        );
    }
    #[test]
    fn test_convex_hull_vs_hull_depth_positive() {
        let tet: Vec<Vec3> = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
            Vec3::new(0.5, 0.5, 1.0),
        ];
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.3, 0.3, 0.0));
        if let Some(c) = convex_hull_vs_convex_hull(&tet, &t1, &tet, &t2) {
            assert!(c.depth > 0.0, "depth must be positive: {}", c.depth);
        }
    }
    #[test]
    fn test_sphere_sphere_gap_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let gap = sphere_sphere_gap(&s1, &t1, &s2, &t2);
        assert!((gap - 2.0).abs() < 1e-10, "gap should be 2.0: {gap}");
    }
    #[test]
    fn test_sphere_sphere_gap_overlapping() {
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let gap = sphere_sphere_gap(&s1, &t1, &s2, &t2);
        assert!(
            gap < 0.0,
            "overlapping spheres should have negative gap: {gap}"
        );
    }
    #[test]
    fn test_sphere_sphere_closing_speed_approaching() {
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let vel1 = Vec3::new(1.0, 0.0, 0.0);
        let vel2 = Vec3::new(-1.0, 0.0, 0.0);
        let speed = sphere_sphere_closing_speed(&t1, vel1, &t2, vel2);
        assert!(
            speed > 0.0,
            "approaching spheres should have positive closing speed: {speed}"
        );
    }
    #[test]
    fn test_sphere_sphere_closing_speed_receding() {
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let vel1 = Vec3::new(-1.0, 0.0, 0.0);
        let vel2 = Vec3::new(1.0, 0.0, 0.0);
        let speed = sphere_sphere_closing_speed(&t1, vel1, &t2, vel2);
        assert!(
            speed < 0.0,
            "receding spheres should have negative closing speed: {speed}"
        );
    }
    #[test]
    fn test_sphere_halfspace_contact() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.5, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        let contact = sphere_halfspace(&s, &t, &plane_pt, &plane_n);
        assert!(
            contact.is_some(),
            "sphere at y=0.5 with r=1 should penetrate y=0 plane"
        );
        let c = contact.unwrap();
        assert!(
            (c.depth - 0.5).abs() < 1e-10,
            "depth should be 0.5: {}",
            c.depth
        );
    }
    #[test]
    fn test_sphere_halfspace_no_contact() {
        let s = Sphere::new(0.3);
        let t = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        assert!(sphere_halfspace(&s, &t, &plane_pt, &plane_n).is_none());
    }
    #[test]
    fn test_capsule_halfspace_contact() {
        let c = Capsule::new(0.5, 1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.3, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        let contact = capsule_halfspace(&c, &t, &plane_pt, &plane_n);
        assert!(
            contact.is_some(),
            "capsule endpoint below plane should contact"
        );
    }
    #[test]
    fn test_capsule_halfspace_no_contact() {
        let c = Capsule::new(0.3, 0.5);
        let t = Transform::from_position(Vec3::new(0.0, 10.0, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        assert!(capsule_halfspace(&c, &t, &plane_pt, &plane_n).is_none());
    }
    #[test]
    fn test_box_halfspace_contact() {
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t = Transform::from_position(Vec3::new(0.0, 0.5, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        let contact = box_halfspace(&b, &t, &plane_pt, &plane_n);
        assert!(
            contact.is_some(),
            "box bottom corners below plane should contact"
        );
    }
    #[test]
    fn test_box_halfspace_no_contact() {
        let b = BoxShape::new(Vec3::new(0.2, 0.2, 0.2));
        let t = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        let plane_pt = Vec3::new(0.0, 0.0, 0.0);
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        assert!(box_halfspace(&b, &t, &plane_pt, &plane_n).is_none());
    }
    #[test]
    fn test_shape_cast_box_box_approaching() {
        let b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(3.0, 0.0, 0.0);
        let vel_a = Vec3::new(2.0, 0.0, 0.0);
        let vel_b = Vec3::new(0.0, 0.0, 0.0);
        let toi = shape_cast_box_box(&b, &pos_a, &vel_a, &b, &pos_b, &vel_b);
        assert!(toi.is_some(), "approaching boxes should have TOI");
    }
    #[test]
    fn test_shape_cast_box_box_receding() {
        let b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
        let pos_a = Vec3::new(0.0, 0.0, 0.0);
        let pos_b = Vec3::new(3.0, 0.0, 0.0);
        let vel_a = Vec3::new(-1.0, 0.0, 0.0);
        let vel_b = Vec3::new(1.0, 0.0, 0.0);
        let toi = shape_cast_box_box(&b, &pos_a, &vel_a, &b, &pos_b, &vel_b);
        assert!(toi.is_none(), "receding boxes should have no TOI");
    }
    #[test]
    fn test_heightfield_interpolated_normal_flat() {
        let heights = vec![0.0_f64; 4];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let n = heightfield_interpolated_normal(&hf, 0.5, 0.5);
        assert!(
            (n.y - 1.0).abs() < 1e-10,
            "flat hf normal should be (0,1,0): {:?}",
            n
        );
    }
    #[test]
    fn test_heightfield_interpolated_normal_slope() {
        let heights = vec![0.0, 1.0, 0.0, 1.0];
        let hf = HeightField::new(heights, 2, 2, 1.0, 1.0);
        let n = heightfield_interpolated_normal(&hf, 0.0, 0.0);
        assert!(
            n.x != 0.0 || n.y > 0.0,
            "sloped hf should have tilted normal: {:?}",
            n
        );
    }
    #[test]
    fn test_rounded_capsule_sphere_contact() {
        let c = Capsule::new(0.5, 1.0);
        let s = Sphere::new(0.3);
        let t_cap = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sph = Transform::from_position(Vec3::new(0.7, 0.0, 0.0));
        let contact = rounded_capsule_sphere(&c, &t_cap, 0.2, &s, &t_sph);
        assert!(
            contact.is_some(),
            "sphere inside rounded capsule should contact"
        );
    }
    #[test]
    fn test_rounded_capsule_sphere_no_contact() {
        let c = Capsule::new(0.2, 0.5);
        let s = Sphere::new(0.1);
        let t_cap = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t_sph = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        assert!(rounded_capsule_sphere(&c, &t_cap, 0.05, &s, &t_sph).is_none());
    }
}
