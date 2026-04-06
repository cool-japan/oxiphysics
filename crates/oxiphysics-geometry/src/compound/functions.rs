//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::ChildShapeKind;

/// Compute the 3×3 inertia tensor of a sphere about its own center.
#[allow(dead_code)]
pub fn sphere_inertia(mass: f64, r: f64) -> [[f64; 3]; 3] {
    let i = 2.0 / 5.0 * mass * r * r;
    [[i, 0.0, 0.0], [0.0, i, 0.0], [0.0, 0.0, i]]
}
/// Compute the 3×3 inertia tensor of an axis-aligned box about its own center.
///
/// `hx`, `hy`, `hz` are the half-extents.
#[allow(dead_code)]
pub fn box_inertia(mass: f64, hx: f64, hy: f64, hz: f64) -> [[f64; 3]; 3] {
    let i_xx = mass / 3.0 * (hy * hy + hz * hz);
    let i_yy = mass / 3.0 * (hx * hx + hz * hz);
    let i_zz = mass / 3.0 * (hx * hx + hy * hy);
    [[i_xx, 0.0, 0.0], [0.0, i_yy, 0.0], [0.0, 0.0, i_zz]]
}
/// Ray-sphere intersection. Returns (toi, normal).
#[allow(dead_code)]
pub(super) fn ray_sphere(
    origin: [f64; 3],
    dir: [f64; 3],
    radius: f64,
    max_toi: f64,
) -> Option<(f64, [f64; 3])> {
    let a = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
    let b = 2.0 * (origin[0] * dir[0] + origin[1] * dir[1] + origin[2] * dir[2]);
    let c = origin[0] * origin[0] + origin[1] * origin[1] + origin[2] * origin[2] - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    let t = if t1 >= 0.0 { t1 } else { t2 };
    if t < 0.0 || t > max_toi {
        return None;
    }
    let p = [
        origin[0] + dir[0] * t,
        origin[1] + dir[1] * t,
        origin[2] + dir[2] * t,
    ];
    let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    let n = if len > 1e-12 {
        [p[0] / len, p[1] / len, p[2] / len]
    } else {
        [0.0, 1.0, 0.0]
    };
    Some((t, n))
}
/// Ray-box intersection (axis-aligned, centered at origin). Returns (toi, normal).
#[allow(dead_code)]
pub(super) fn ray_box(
    origin: [f64; 3],
    dir: [f64; 3],
    half_extents: [f64; 3],
    max_toi: f64,
) -> Option<(f64, [f64; 3])> {
    let mut tmin = f64::NEG_INFINITY;
    let mut tmax = f64::INFINITY;
    let mut normal = [0.0; 3];
    for i in 0..3 {
        if dir[i].abs() < 1e-12 {
            if origin[i] < -half_extents[i] || origin[i] > half_extents[i] {
                return None;
            }
        } else {
            let t1 = (-half_extents[i] - origin[i]) / dir[i];
            let t2 = (half_extents[i] - origin[i]) / dir[i];
            let (ta, tb, sign) = if t1 < t2 {
                (t1, t2, -1.0)
            } else {
                (t2, t1, 1.0)
            };
            if ta > tmin {
                tmin = ta;
                normal = [0.0; 3];
                normal[i] = sign;
            }
            tmax = tmax.min(tb);
            if tmin > tmax {
                return None;
            }
        }
    }
    if tmin < 0.0 {
        tmin = 0.0;
    }
    if tmin > max_toi {
        return None;
    }
    Some((tmin, normal))
}
/// Ray-capsule intersection (Y-axis aligned, centered at origin).
#[allow(dead_code)]
pub(super) fn ray_capsule(
    origin: [f64; 3],
    dir: [f64; 3],
    radius: f64,
    half_height: f64,
    max_toi: f64,
) -> Option<(f64, [f64; 3])> {
    let mut best: Option<(f64, [f64; 3])> = None;
    let top_o = [origin[0], origin[1] - half_height, origin[2]];
    if let Some((t, n)) = ray_sphere(top_o, dir, radius, max_toi) {
        let hit_y = origin[1] + dir[1] * t;
        if hit_y >= half_height && best.as_ref().is_none_or(|(bt, _)| t < *bt) {
            best = Some((t, n));
        }
    }
    let bot_o = [origin[0], origin[1] + half_height, origin[2]];
    if let Some((t, n)) = ray_sphere(bot_o, dir, radius, max_toi) {
        let hit_y = origin[1] + dir[1] * t;
        if hit_y <= -half_height && best.as_ref().is_none_or(|(bt, _)| t < *bt) {
            best = Some((t, n));
        }
    }
    let a = dir[0] * dir[0] + dir[2] * dir[2];
    let b = 2.0 * (origin[0] * dir[0] + origin[2] * dir[2]);
    let c = origin[0] * origin[0] + origin[2] * origin[2] - radius * radius;
    let disc = b * b - 4.0 * a * c;
    if a > 1e-12 && disc >= 0.0 {
        let sqrt_disc = disc.sqrt();
        for &t in &[(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)] {
            if t >= 0.0 && t <= max_toi {
                let hit_y = origin[1] + dir[1] * t;
                if hit_y.abs() <= half_height {
                    let hx = origin[0] + dir[0] * t;
                    let hz = origin[2] + dir[2] * t;
                    let len = (hx * hx + hz * hz).sqrt();
                    let n = if len > 1e-12 {
                        [hx / len, 0.0, hz / len]
                    } else {
                        [1.0, 0.0, 0.0]
                    };
                    if best.as_ref().is_none_or(|(bt, _)| t < *bt) {
                        best = Some((t, n));
                    }
                }
            }
        }
    }
    best
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Compound;
    use crate::Shape;
    use crate::box_shape::BoxShape;
    use crate::compound::ChildShapeKind;

    use crate::compound::CompoundShape;

    use crate::sphere::Sphere;
    use oxiphysics_core::Real;
    use oxiphysics_core::Transform;
    use oxiphysics_core::math::Vec3;
    use std::f64::consts::PI;
    use std::sync::Arc;
    #[test]
    fn test_compound_volume() {
        let s1: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let s2: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let compound = Compound::new(vec![
            (Transform::default(), s1.clone()),
            (
                Transform::from_position(Vec3::new(5.0, 0.0, 0.0)),
                s2.clone(),
            ),
        ]);
        assert!((compound.volume() - 2.0 * s1.volume()).abs() < 1e-10);
    }
    /// Compound of two unit spheres → volume = 2 * (4/3)π
    #[test]
    fn test_compound_two_spheres_volume() {
        let s1: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let s2: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let compound = Compound::new(vec![
            (Transform::default(), s1),
            (Transform::from_position(Vec3::new(4.0, 0.0, 0.0)), s2),
        ]);
        let expected = 2.0 * (4.0 / 3.0) * PI;
        assert!(
            (compound.volume() - expected).abs() < 1e-6,
            "volume={} expected={}",
            compound.volume(),
            expected
        );
    }
    /// Two equal spheres on x-axis at ±d → I_yy = 2*(I_sphere + m_child*d²)
    #[test]
    fn test_compound_inertia_parallel_axis() {
        let r = 1.0_f64;
        let d = 3.0_f64;
        let total_mass = 10.0_f64;
        let s1: Arc<dyn Shape> = Arc::new(Sphere::new(r));
        let s2: Arc<dyn Shape> = Arc::new(Sphere::new(r));
        let compound = Compound::new(vec![
            (Transform::from_position(Vec3::new(d, 0.0, 0.0)), s1),
            (Transform::from_position(Vec3::new(-d, 0.0, 0.0)), s2),
        ]);
        let m_child = total_mass / 2.0;
        let i_sphere = 2.0 / 5.0 * m_child * r * r;
        let expected_iyy = 2.0 * (i_sphere + m_child * d * d);
        let inertia = compound.inertia_tensor(total_mass);
        let iyy = inertia[(1, 1)];
        assert!(
            (iyy - expected_iyy).abs() < 1e-2,
            "I_yy={} expected={}",
            iyy,
            expected_iyy
        );
    }
    /// Compound with box at origin and sphere at (5,0,0); ray from (10,0,0) toward -x → hits sphere first (TOI≈4)
    #[test]
    fn test_compound_raycast_hits_child() {
        let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
        let sphere: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let compound = Compound::new(vec![
            (Transform::default(), box_shape),
            (Transform::from_position(Vec3::new(5.0, 0.0, 0.0)), sphere),
        ]);
        let origin = Vec3::new(10.0, 0.0, 0.0);
        let direction = Vec3::new(-1.0, 0.0, 0.0);
        let hit = compound.ray_cast(&origin, &direction, 100.0);
        assert!(hit.is_some(), "ray should hit the compound shape");
        let hit = hit.unwrap();
        assert!((hit.toi - 4.0).abs() < 1e-2, "toi={} expected≈4.0", hit.toi);
    }
    /// Compound of N shapes → total mass = sum of child masses
    #[test]
    fn test_compound_mass_properties_additive() {
        let density = 500.0_f64;
        let shapes: Vec<Arc<dyn Shape>> = vec![
            Arc::new(Sphere::new(1.0)),
            Arc::new(Sphere::new(0.5)),
            Arc::new(BoxShape::new(Vec3::new(1.0, 1.0, 1.0))),
        ];
        let offsets = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(-3.0, 0.0, 0.0),
        ];
        let expected_total_mass: Real = shapes.iter().map(|s| density * s.volume()).sum();
        let children: Vec<(Transform, Arc<dyn Shape>)> = offsets
            .iter()
            .zip(shapes.iter())
            .map(|(pos, s)| (Transform::from_position(*pos), s.clone()))
            .collect();
        let compound = Compound::new(children);
        let props = compound.mass_properties(density);
        assert!(
            (props.mass - expected_total_mass).abs() < 1e-6,
            "compound mass={} expected={}",
            props.mass,
            expected_total_mass
        );
    }
    #[test]
    fn test_compound_shape_single_sphere_volume() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 2.0);
        let expected = (4.0 / 3.0) * PI * 8.0;
        assert!(
            (cs.total_volume() - expected).abs() < 1e-6,
            "vol={} expected={}",
            cs.total_volume(),
            expected
        );
    }
    #[test]
    fn test_compound_shape_box_volume() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        assert!(
            (cs.total_volume() - 48.0).abs() < 1e-10,
            "vol={}",
            cs.total_volume()
        );
    }
    #[test]
    fn test_compound_shape_aabb_single_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([1.0, 2.0, 3.0], 0.5);
        let (min, max) = cs.aabb();
        assert!((min[0] - 0.5).abs() < 1e-10);
        assert!((min[1] - 1.5).abs() < 1e-10);
        assert!((min[2] - 2.5).abs() < 1e-10);
        assert!((max[0] - 1.5).abs() < 1e-10);
        assert!((max[1] - 2.5).abs() < 1e-10);
        assert!((max[2] - 3.5).abs() < 1e-10);
    }
    #[test]
    fn test_compound_shape_aabb_multiple() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-5.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 1.0);
        let (min, max) = cs.aabb();
        assert!((min[0] - (-6.0)).abs() < 1e-10);
        assert!((max[0] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_compound_shape_com_symmetric() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-3.0, 0.0, 0.0], 1.0);
        cs.add_sphere([3.0, 0.0, 0.0], 1.0);
        let com = cs.center_of_mass();
        assert!((com[0]).abs() < 1e-10, "com_x={}", com[0]);
        assert!((com[1]).abs() < 1e-10, "com_y={}", com[1]);
        assert!((com[2]).abs() < 1e-10, "com_z={}", com[2]);
    }
    #[test]
    fn test_compound_shape_com_weighted() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let com = cs.center_of_mass();
        assert!(
            (com[0] - 5.0).abs() < 1e-10,
            "com_x={} expected 5.0",
            com[0]
        );
    }
    #[test]
    fn test_compound_shape_ray_through_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let origin = [-5.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let hit = cs.ray_cast(origin, dir, 100.0);
        assert!(hit.is_some(), "should hit sphere");
        let (toi, normal, idx) = hit.unwrap();
        assert_eq!(idx, 0);
        assert!((toi - 4.0).abs() < 1e-6, "toi={} expected 4.0", toi);
        assert!((normal[0] - (-1.0)).abs() < 1e-6, "normal_x={}", normal[0]);
    }
    #[test]
    fn test_compound_shape_ray_misses() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let origin = [-5.0, 0.0, 0.0];
        let dir = [-1.0, 0.0, 0.0];
        let hit = cs.ray_cast(origin, dir, 100.0);
        assert!(hit.is_none(), "should not hit sphere from behind");
    }
    #[test]
    fn test_compound_shape_ray_closest_child() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([5.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let origin = [0.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let hit = cs.ray_cast(origin, dir, 100.0);
        assert!(hit.is_some());
        let (_, _, idx) = hit.unwrap();
        assert_eq!(idx, 0, "should hit closer sphere first");
    }
    #[test]
    fn test_compound_shape_contains_point_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(cs.contains_point([0.0, 0.0, 0.0]));
        assert!(cs.contains_point([1.0, 0.0, 0.0]));
        assert!(!cs.contains_point([3.0, 0.0, 0.0]));
    }
    #[test]
    fn test_compound_shape_contains_point_box() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(cs.contains_point([0.5, 0.5, 0.5]));
        assert!(!cs.contains_point([1.5, 0.0, 0.0]));
    }
    #[test]
    fn test_compound_shape_child_count() {
        let mut cs = CompoundShape::new();
        assert_eq!(cs.child_count(), 0);
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_box([1.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        cs.add_capsule([2.0, 0.0, 0.0], 0.5, 1.0);
        assert_eq!(cs.child_count(), 3);
    }
    #[test]
    fn test_compound_shape_capsule_volume() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 2.0);
        let sphere_vol = (4.0 / 3.0) * PI;
        let cyl_vol = PI * 4.0;
        let expected = sphere_vol + cyl_vol;
        assert!(
            (cs.total_volume() - expected).abs() < 1e-6,
            "vol={} expected={}",
            cs.total_volume(),
            expected
        );
    }
    #[test]
    fn test_compound_shape_ray_cast_box() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let origin = [-5.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let hit = cs.ray_cast(origin, dir, 100.0);
        assert!(hit.is_some());
        let (toi, normal, _) = hit.unwrap();
        assert!((toi - 4.0).abs() < 1e-6, "toi={} expected 4.0", toi);
        assert!((normal[0] - (-1.0)).abs() < 1e-6);
    }
    #[test]
    fn test_compound_shape_inertia_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let density = 1.0;
        let inertia = cs.inertia_tensor(density);
        let vol = cs.total_volume();
        let mass = density * vol;
        let expected_i = 2.0 / 5.0 * mass;
        assert!(
            (inertia[0][0] - expected_i).abs() < 1e-6,
            "I_xx={}",
            inertia[0][0]
        );
        assert!(
            (inertia[1][1] - expected_i).abs() < 1e-6,
            "I_yy={}",
            inertia[1][1]
        );
        assert!(
            (inertia[2][2] - expected_i).abs() < 1e-6,
            "I_zz={}",
            inertia[2][2]
        );
    }
    #[test]
    fn test_compound_shape_inertia_off_diagonal_zero_for_symmetric() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let inertia = cs.inertia_tensor(1.0);
        assert!(inertia[0][1].abs() < 1e-10);
        assert!(inertia[0][2].abs() < 1e-10);
        assert!(inertia[1][2].abs() < 1e-10);
    }
    #[test]
    fn test_bounding_sphere_single_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 2.0);
        let (center, radius) = cs.bounding_sphere();
        assert!((center[0]).abs() < 1e-10);
        assert!(radius >= 2.0);
    }
    #[test]
    fn test_bounding_sphere_two_spheres() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-3.0, 0.0, 0.0], 1.0);
        cs.add_sphere([3.0, 0.0, 0.0], 1.0);
        let (_, radius) = cs.bounding_sphere();
        assert!(radius >= 4.0, "radius={}", radius);
    }
    #[test]
    fn test_merge_with_adds_children() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 1.0);
        let mut cs2 = CompoundShape::new();
        cs2.add_box([2.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let merged = cs1.merge_with(&cs2);
        assert_eq!(merged.child_count(), 2);
    }
    #[test]
    fn test_scale_doubles_radius() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([1.0, 0.0, 0.0], 1.0);
        cs.scale(2.0);
        match cs.children[0].shape_kind {
            ChildShapeKind::Sphere { radius } => assert!((radius - 2.0).abs() < 1e-10),
            _ => panic!("expected sphere"),
        }
        assert!((cs.children[0].center[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_translate_shifts_centers() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.translate([1.0, 2.0, 3.0]);
        assert!((cs.children[0].center[0] - 1.0).abs() < 1e-10);
        assert!((cs.children[0].center[1] - 2.0).abs() < 1e-10);
        assert!((cs.children[0].center[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_overlaps_sphere_hit() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(cs.overlaps_sphere([0.5, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_overlaps_sphere_miss() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 0.5);
        assert!(!cs.overlaps_sphere([5.0, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_ray_cast_all_returns_both() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([2.0, 0.0, 0.0], 0.5);
        cs.add_sphere([5.0, 0.0, 0.0], 0.5);
        let hits = cs.ray_cast_all([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert_eq!(hits.len(), 2, "should hit both spheres");
        assert!(hits[0].0 < hits[1].0, "first hit should be closer");
    }
    #[test]
    fn test_merged_aabb_contains_children() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-5.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 1.0);
        let (mn, mx) = cs.merged_aabb();
        assert!(mn[0] <= -6.0 + 1e-9, "min_x={}", mn[0]);
        assert!(mx[0] >= 6.0 - 1e-9, "max_x={}", mx[0]);
    }
    #[test]
    fn test_compound_aabb_struct_all_aabbs() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_box([3.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let all: Vec<_> = cs
            .children
            .iter()
            .map(CompoundShape::child_aabb_public)
            .collect();
        assert_eq!(all.len(), 2);
        assert!((all[0].0[0] - (-1.0)).abs() < 1e-10);
        assert!((all[1].0[0] - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_raycast_returns_t_and_index() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let result = cs.raycast([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some());
        let (t, idx) = result.unwrap();
        assert_eq!(idx, 0);
        assert!((t - 4.0).abs() < 1e-6, "t={}", t);
    }
    #[test]
    fn test_volume_sum() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 2.0);
        let v1 = (4.0 / 3.0) * PI;
        let v2 = (4.0 / 3.0) * PI * 8.0;
        assert!(
            (cs.volume() - (v1 + v2)).abs() < 1e-6,
            "vol={}",
            cs.volume()
        );
    }
    #[test]
    fn test_center_of_mass_weighted_average() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let com = cs.center_of_mass_weighted(&[1.0, 1.0]);
        assert!((com[0] - 5.0).abs() < 1e-10, "com_x={}", com[0]);
        let com2 = cs.center_of_mass_weighted(&[1.0, 3.0]);
        assert!((com2[0] - 7.5).abs() < 1e-10, "com2_x={}", com2[0]);
    }
    #[test]
    fn test_inertia_tensor_from_masses_symmetric() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([4.0, 0.0, 0.0], 1.0);
        let masses = [1.0, 1.0];
        let i = cs.inertia_tensor_from_masses(&masses);
        assert!(
            (i[0][1] - i[1][0]).abs() < 1e-10,
            "not symmetric [0][1] vs [1][0]"
        );
        assert!(
            (i[0][2] - i[2][0]).abs() < 1e-10,
            "not symmetric [0][2] vs [2][0]"
        );
        assert!(
            (i[1][2] - i[2][1]).abs() < 1e-10,
            "not symmetric [1][2] vs [2][1]"
        );
    }
    #[test]
    fn test_closest_point_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let (cp, idx) = cs.closest_point([3.0, 0.0, 0.0]);
        assert_eq!(idx, 0);
        assert!((cp[0] - 1.0).abs() < 1e-9, "cp_x={}", cp[0]);
        assert!(cp[1].abs() < 1e-9);
        assert!(cp[2].abs() < 1e-9);
    }
    #[test]
    fn test_closest_point_selects_nearest_child() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-10.0, 0.0, 0.0], 1.0);
        cs.add_sphere([2.0, 0.0, 0.0], 1.0);
        let (_cp, idx) = cs.closest_point([4.0, 0.0, 0.0]);
        assert_eq!(idx, 1, "should select nearer child");
    }
    #[test]
    fn test_sphere_inertia_helper() {
        let i = sphere_inertia(5.0, 2.0);
        for k in 0..3 {
            assert!((i[k][k] - 8.0).abs() < 1e-10, "I[{k}][{k}]={}", i[k][k]);
        }
        assert!(i[0][1].abs() < 1e-15);
    }
    #[test]
    fn test_box_inertia_helper() {
        let i = box_inertia(3.0, 1.0, 2.0, 3.0);
        assert!((i[0][0] - 13.0).abs() < 1e-10, "I_xx={}", i[0][0]);
        assert!(i[0][1].abs() < 1e-15);
    }
}
/// Check if a point in local space is inside a `ChildShapeKind`.
pub(super) fn child_kind_contains(kind: &ChildShapeKind, p: [f64; 3]) -> bool {
    match kind {
        ChildShapeKind::Sphere { radius } => {
            p[0] * p[0] + p[1] * p[1] + p[2] * p[2] <= radius * radius
        }
        ChildShapeKind::Box { half_extents } => {
            p[0].abs() <= half_extents[0]
                && p[1].abs() <= half_extents[1]
                && p[2].abs() <= half_extents[2]
        }
        ChildShapeKind::Capsule {
            radius,
            half_height,
        } => {
            let clamped_y = p[1].clamp(-half_height, *half_height);
            let ry = p[1] - clamped_y;
            p[0] * p[0] + ry * ry + p[2] * p[2] <= radius * radius
        }
    }
}
/// Ray cast against a `ChildShapeKind` at the origin.
pub(super) fn ray_cast_kind(
    kind: &ChildShapeKind,
    origin: [f64; 3],
    dir: [f64; 3],
    max_toi: f64,
) -> Option<(f64, [f64; 3])> {
    match kind {
        ChildShapeKind::Sphere { radius } => ray_sphere(origin, dir, *radius, max_toi),
        ChildShapeKind::Box { half_extents } => ray_box(origin, dir, *half_extents, max_toi),
        ChildShapeKind::Capsule {
            radius,
            half_height,
        } => ray_capsule(origin, dir, *radius, *half_height, max_toi),
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::compound::ChildShapeKind;

    use crate::compound::CompoundShapeEx;
    use crate::compound::LocalTransform;
    use crate::compound::child_kind_contains;

    use std::f64::consts::PI;

    #[test]
    fn test_local_to_world_identity() {
        let t = LocalTransform::identity();
        let p = [1.0, 2.0, 3.0];
        let w = t.local_to_world(p);
        assert!((w[0] - 1.0).abs() < 1e-12);
        assert!((w[1] - 2.0).abs() < 1e-12);
        assert!((w[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_world_to_local_identity() {
        let t = LocalTransform::identity();
        let p = [5.0, -3.0, 7.0];
        let l = t.world_to_local(p);
        assert!((l[0] - 5.0).abs() < 1e-12);
        assert!((l[1] - (-3.0)).abs() < 1e-12);
        assert!((l[2] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_local_to_world_translation() {
        let t = LocalTransform::from_translation([10.0, 20.0, 30.0]);
        let p = [1.0, 0.0, 0.0];
        let w = t.local_to_world(p);
        assert!((w[0] - 11.0).abs() < 1e-12);
        assert!((w[1] - 20.0).abs() < 1e-12);
        assert!((w[2] - 30.0).abs() < 1e-12);
    }
    #[test]
    fn test_world_to_local_translation_roundtrip() {
        let t = LocalTransform::from_translation([3.0, -1.0, 5.0]);
        let world_p = [7.0, 4.0, 8.0];
        let local_p = t.world_to_local(world_p);
        let back = t.local_to_world(local_p);
        for i in 0..3 {
            assert!(
                (back[i] - world_p[i]).abs() < 1e-10,
                "axis {i}: {} != {}",
                back[i],
                world_p[i]
            );
        }
    }
    #[test]
    fn test_local_to_world_rotation_90_deg_y() {
        let t = LocalTransform {
            translation: [0.0; 3],
            rot: [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        };
        let p = [1.0, 0.0, 0.0];
        let w = t.local_to_world(p);
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "rotation should preserve length");
    }
    #[test]
    fn test_compound_ex_aabb_single_sphere_identity() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::identity(), 2.0);
        let (mn, mx) = cs.aabb();
        assert!((mn[0] - (-2.0)).abs() < 1e-10, "min_x={}", mn[0]);
        assert!((mx[0] - 2.0).abs() < 1e-10, "max_x={}", mx[0]);
    }
    #[test]
    fn test_compound_ex_aabb_translated_sphere() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::from_translation([5.0, 0.0, 0.0]), 1.0);
        let (mn, mx) = cs.aabb();
        assert!((mn[0] - 4.0).abs() < 1e-10, "min_x={}", mn[0]);
        assert!((mx[0] - 6.0).abs() < 1e-10, "max_x={}", mx[0]);
    }
    #[test]
    fn test_compound_ex_contains_point_sphere() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::from_translation([3.0, 0.0, 0.0]), 1.5);
        assert!(
            cs.contains_point([3.0, 0.0, 0.0]),
            "center should be inside"
        );
        assert!(!cs.contains_point([0.0, 0.0, 0.0]), "origin is outside");
    }
    #[test]
    fn test_compound_ex_contains_point_box() {
        let mut cs = CompoundShapeEx::new();
        cs.add_box(LocalTransform::identity(), [1.0, 2.0, 3.0]);
        assert!(cs.contains_point([0.5, 1.0, 2.0]));
        assert!(!cs.contains_point([1.5, 0.0, 0.0]));
    }
    #[test]
    fn test_compound_ex_ray_cast_sphere() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::from_translation([5.0, 0.0, 0.0]), 1.0);
        let hit = cs.ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_some(), "should hit sphere");
        let (t, _) = hit.unwrap();
        assert!((t - 4.0).abs() < 1e-6, "t={t} expected 4.0");
    }
    #[test]
    fn test_compound_ex_volume_sphere() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::identity(), 1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((cs.volume() - expected).abs() < 1e-9, "vol={}", cs.volume());
    }
    #[test]
    fn test_compound_ex_inertia_tensor_sphere_at_origin() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::identity(), 1.0);
        let density = 1.0;
        let i = cs.inertia_tensor(density);
        let vol = (4.0 / 3.0) * PI;
        let mass = density * vol;
        let expected = 2.0 / 5.0 * mass;
        assert!((i[0][0] - expected).abs() < 1e-9, "I_xx={}", i[0][0]);
        assert!((i[1][1] - expected).abs() < 1e-9, "I_yy={}", i[1][1]);
        assert!((i[2][2] - expected).abs() < 1e-9, "I_zz={}", i[2][2]);
    }
    #[test]
    fn test_compound_ex_inertia_tensor_parallel_axis() {
        let d = 3.0_f64;
        let density = 1.0;
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::from_translation([d, 0.0, 0.0]), 1.0);
        let i = cs.inertia_tensor(density);
        let vol = (4.0 / 3.0) * PI;
        let mass = density * vol;
        let expected_iyy = 2.0 / 5.0 * mass + mass * d * d;
        assert!(
            (i[1][1] - expected_iyy).abs() < 1e-6,
            "I_yy={} expected={}",
            i[1][1],
            expected_iyy
        );
    }
    #[test]
    fn test_child_kind_contains_sphere() {
        let k = ChildShapeKind::Sphere { radius: 2.0 };
        assert!(child_kind_contains(&k, [0.0, 0.0, 0.0]));
        assert!(child_kind_contains(&k, [1.9, 0.0, 0.0]));
        assert!(!child_kind_contains(&k, [2.1, 0.0, 0.0]));
    }
    #[test]
    fn test_child_kind_contains_box() {
        let k = ChildShapeKind::Box {
            half_extents: [1.0, 2.0, 3.0],
        };
        assert!(child_kind_contains(&k, [0.5, 1.5, 2.5]));
        assert!(!child_kind_contains(&k, [1.5, 0.0, 0.0]));
    }
    #[test]
    fn test_child_kind_contains_capsule() {
        let k = ChildShapeKind::Capsule {
            radius: 1.0,
            half_height: 2.0,
        };
        assert!(child_kind_contains(&k, [0.5, 1.0, 0.5]));
        assert!(!child_kind_contains(&k, [2.0, 0.0, 0.0]));
    }
}
#[cfg(test)]
mod tests_extended2 {

    use crate::compound::ChildShapeKind;
    use crate::compound::CompoundChild;
    use crate::compound::CompoundShape;

    use std::f64::consts::PI;

    #[test]
    fn test_remove_child_decreases_count() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 1.0);
        assert_eq!(cs.child_count(), 2);
        cs.remove_child(0);
        assert_eq!(cs.child_count(), 1);
    }
    #[test]
    fn test_swap_remove_child() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 2.0);
        cs.add_sphere([10.0, 0.0, 0.0], 3.0);
        cs.swap_remove_child(0);
        assert_eq!(cs.child_count(), 2);
    }
    #[test]
    fn test_replace_with_sphere_changes_kind() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        cs.replace_with_sphere(0, [1.0, 2.0, 3.0], 0.5);
        match cs.children[0].shape_kind {
            ChildShapeKind::Sphere { radius } => {
                assert!((radius - 0.5).abs() < 1e-10, "radius should be 0.5");
            }
            _ => panic!("expected Sphere after replace"),
        }
        assert!((cs.children[0].center[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_replace_with_box_changes_kind() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.replace_with_box(0, [2.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        match cs.children[0].shape_kind {
            ChildShapeKind::Box { half_extents } => {
                assert!((half_extents[0] - 0.5).abs() < 1e-10);
            }
            _ => panic!("expected Box after replace"),
        }
    }
    #[test]
    fn test_is_empty_and_clear() {
        let mut cs = CompoundShape::new();
        assert!(cs.is_empty(), "newly created should be empty");
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(!cs.is_empty());
        cs.clear();
        assert!(cs.is_empty(), "should be empty after clear");
    }
    #[test]
    fn test_closest_point_with_dist2_sphere() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let (cp, d2, idx) = cs.closest_point_with_dist2([3.0, 0.0, 0.0]);
        assert_eq!(idx, 0);
        assert!((cp[0] - 1.0).abs() < 1e-9);
        assert!((d2 - 4.0).abs() < 1e-9, "dist2 should be 4, got {d2}");
    }
    #[test]
    fn test_broad_phase_pairs_overlap() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 2.0);
        let mut cs2 = CompoundShape::new();
        cs2.add_sphere([1.0, 0.0, 0.0], 2.0);
        let pairs = cs1.broad_phase_pairs(&cs2);
        assert!(
            !pairs.is_empty(),
            "overlapping spheres should give broad-phase pair"
        );
    }
    #[test]
    fn test_broad_phase_pairs_no_overlap() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 0.5);
        let mut cs2 = CompoundShape::new();
        cs2.add_sphere([100.0, 0.0, 0.0], 0.5);
        let pairs = cs1.broad_phase_pairs(&cs2);
        assert!(
            pairs.is_empty(),
            "distant spheres should not produce broad-phase pairs"
        );
    }
    #[test]
    fn test_overlaps_compound_hit() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 1.5);
        let mut cs2 = CompoundShape::new();
        cs2.add_sphere([1.0, 0.0, 0.0], 1.5);
        assert!(cs1.overlaps_compound(&cs2));
    }
    #[test]
    fn test_overlaps_compound_miss() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 0.5);
        let mut cs2 = CompoundShape::new();
        cs2.add_sphere([50.0, 0.0, 0.0], 0.5);
        assert!(!cs1.overlaps_compound(&cs2));
    }
    #[test]
    fn test_centroid_with_densities_equal() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let c = cs.centroid_with_densities(&[1.0, 1.0]);
        assert!((c[0] - 5.0).abs() < 1e-9, "centroid_x={}", c[0]);
    }
    #[test]
    fn test_centroid_with_densities_unequal() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let c = cs.centroid_with_densities(&[1.0, 9.0]);
        assert!(
            c[0] > 5.0,
            "centroid should be above 5 with heavier right child"
        );
    }
    #[test]
    fn test_penetration_depth_sphere_penetrates() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 2.0);
        let result = cs.penetration_depth_sphere([0.5, 0.0, 0.0], 2.0);
        assert!(result.is_some(), "should detect penetration");
        let (depth, idx) = result.unwrap();
        assert_eq!(idx, 0);
        assert!(
            depth < 0.0,
            "depth should be negative for penetration, got {depth}"
        );
    }
    #[test]
    fn test_penetration_depth_sphere_no_penetration() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let result = cs.penetration_depth_sphere([10.0, 0.0, 0.0], 0.5);
        assert!(result.is_none(), "no penetration for distant sphere");
    }
    #[test]
    fn test_child_masses_sum() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 2.0);
        let density = 3.0;
        let masses = cs.child_masses(density);
        let expected: f64 = cs.total_mass(density);
        let actual: f64 = masses.iter().sum();
        assert!((actual - expected).abs() < 1e-9, "mass sum mismatch");
    }
    #[test]
    fn test_total_mass() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let density = 2.0;
        let expected = density * (4.0 / 3.0) * PI;
        let m = cs.total_mass(density);
        assert!(
            (m - expected).abs() < 1e-9,
            "total_mass={m} expected={expected}"
        );
    }
    #[test]
    fn test_child_aabb_public() {
        let child = CompoundChild {
            center: [1.0, 2.0, 3.0],
            shape_kind: ChildShapeKind::Sphere { radius: 1.0 },
        };
        let (mn, mx) = CompoundShape::child_aabb_public(&child);
        assert!((mn[0] - 0.0).abs() < 1e-10);
        assert!((mx[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_expanded_aabb_increases_size() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let (mn0, mx0) = cs.aabb();
        let (mn1, mx1) = cs.expanded_aabb(0.5);
        for i in 0..3 {
            assert!(mn1[i] < mn0[i], "expanded min should be smaller");
            assert!(mx1[i] > mx0[i], "expanded max should be larger");
        }
    }
    #[test]
    fn test_sphere_overlaps_aabb_inside() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(cs.sphere_overlaps_aabb([0.0, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_sphere_overlaps_aabb_outside() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(!cs.sphere_overlaps_aabb([100.0, 0.0, 0.0], 0.5));
    }
}
#[cfg(test)]
mod tests_extra {

    use crate::compound::ChildShapeKind;

    use crate::compound::CompoundShape;
    use crate::compound::CompoundShapeEx;
    use crate::compound::LocalTransform;

    use std::f64::consts::PI;

    #[test]
    fn test_aabb_union_three_spheres() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([-5.0, 0.0, 0.0], 1.0);
        cs.add_sphere([0.0, 5.0, 0.0], 1.0);
        cs.add_sphere([0.0, 0.0, 5.0], 1.0);
        let (mn, mx) = cs.aabb();
        assert!(mn[0] <= -6.0 + 1e-9);
        assert!(mx[1] >= 6.0 - 1e-9);
        assert!(mx[2] >= 6.0 - 1e-9);
    }
    #[test]
    fn test_aabb_union_empty_compound() {
        let cs = CompoundShape::new();
        let (mn, mx) = cs.aabb();
        for i in 0..3 {
            assert!((mn[i]).abs() < 1e-12);
            assert!((mx[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_aabb_tight_for_single_box() {
        let mut cs = CompoundShape::new();
        cs.add_box([3.0, 1.0, -2.0], [2.0, 1.0, 0.5]);
        let (mn, mx) = cs.aabb();
        assert!((mn[0] - 1.0).abs() < 1e-9, "min_x={}", mn[0]);
        assert!((mx[0] - 5.0).abs() < 1e-9, "max_x={}", mx[0]);
        assert!((mn[1] - 0.0).abs() < 1e-9, "min_y={}", mn[1]);
        assert!((mx[1] - 2.0).abs() < 1e-9, "max_y={}", mx[1]);
    }
    #[test]
    fn test_aabb_capsule() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 3.0);
        let (mn, mx) = cs.aabb();
        assert!((mn[0] - (-1.0)).abs() < 1e-9);
        assert!(
            (mx[1] - 4.0).abs() < 1e-9,
            "max_y for capsule with hh=3, r=1 should be 4, got {}",
            mx[1]
        );
    }
    #[test]
    fn test_raycast_selects_first_among_three_children() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([2.0, 0.0, 0.0], 0.5);
        cs.add_sphere([5.0, 0.0, 0.0], 0.5);
        cs.add_sphere([8.0, 0.0, 0.0], 0.5);
        let result = cs.raycast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_some());
        let (_, idx) = result.unwrap();
        assert_eq!(idx, 0, "should hit first (closest) sphere, got idx={idx}");
    }
    #[test]
    fn test_raycast_all_sorts_by_toi() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([3.0, 0.0, 0.0], 0.5);
        cs.add_sphere([7.0, 0.0, 0.0], 0.5);
        let hits = cs.ray_cast_all([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert_eq!(hits.len(), 2, "should hit both spheres");
        assert!(hits[0].0 < hits[1].0, "hits should be sorted by toi");
    }
    #[test]
    fn test_ray_cast_none_when_all_behind_origin() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let result = cs.ray_cast([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(result.is_none(), "should not hit sphere behind ray origin");
    }
    #[test]
    fn test_ray_cast_box_all_six_faces() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let h = cs.ray_cast([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 100.0);
        assert!(h.is_some(), "+X face miss");
        let h = cs.ray_cast([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0);
        assert!(h.is_some(), "-X face miss");
        let h = cs.ray_cast([0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 100.0);
        assert!(h.is_some(), "+Y face miss");
        let h = cs.ray_cast([0.0, 0.0, -5.0], [0.0, 0.0, 1.0], 100.0);
        assert!(h.is_some(), "-Z face miss");
    }
    #[test]
    fn test_total_volume_scales_with_sphere_radius() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 1.0);
        let mut cs2 = CompoundShape::new();
        cs2.add_sphere([0.0, 0.0, 0.0], 2.0);
        let ratio = cs2.total_volume() / cs1.total_volume();
        assert!(
            (ratio - 8.0).abs() < 1e-6,
            "volume ratio should be 8, got {ratio}"
        );
    }
    #[test]
    fn test_total_volume_box_scales_with_half_extents() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let expected = 8.0 * 2.0 * 3.0 * 4.0;
        assert!(
            (cs.total_volume() - expected).abs() < 1e-9,
            "box volume={}, expected={expected}",
            cs.total_volume()
        );
    }
    #[test]
    fn test_total_mass_proportional_to_density() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let m1 = cs.total_mass(1.0);
        let m2 = cs.total_mass(3.0);
        assert!(
            (m2 / m1 - 3.0).abs() < 1e-9,
            "mass should scale with density"
        );
    }
    #[test]
    fn test_child_masses_length_matches_children() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_box([5.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        cs.add_capsule([10.0, 0.0, 0.0], 0.5, 1.0);
        let masses = cs.child_masses(2.0);
        assert_eq!(masses.len(), 3);
        for &m in &masses {
            assert!(m > 0.0, "all child masses should be positive");
        }
    }
    #[test]
    fn test_centroid_with_densities_single_child() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([3.0, 4.0, 5.0], 1.0);
        let c = cs.centroid_with_densities(&[1.0]);
        assert!((c[0] - 3.0).abs() < 1e-9, "cx={}", c[0]);
        assert!((c[1] - 4.0).abs() < 1e-9, "cy={}", c[1]);
        assert!((c[2] - 5.0).abs() < 1e-9, "cz={}", c[2]);
    }
    #[test]
    fn test_centroid_uses_default_density_for_missing_entries() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let c = cs.centroid_with_densities(&[1.0]);
        assert!((c[0] - 5.0).abs() < 1e-9, "cx={}", c[0]);
    }
    #[test]
    fn test_inertia_tensor_single_sphere_at_origin_all_diagonal() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        let i = cs.inertia_tensor(1.0);
        assert!(i[0][1].abs() < 1e-10);
        assert!(i[0][2].abs() < 1e-10);
        assert!(i[1][2].abs() < 1e-10);
    }
    #[test]
    fn test_inertia_tensor_from_masses_single_box() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let i = cs.inertia_tensor_from_masses(&[10.0]);
        let expected_xx = 10.0 / 3.0 * (4.0 + 9.0);
        assert!(
            (i[0][0] - expected_xx).abs() < 1e-6,
            "I_xx={}, expected={expected_xx}",
            i[0][0]
        );
    }
    #[test]
    fn test_inertia_tensor_two_equal_spheres_parallel_axis() {
        let d = 2.0_f64;
        let mut cs = CompoundShape::new();
        cs.add_sphere([-d, 0.0, 0.0], 1.0);
        cs.add_sphere([d, 0.0, 0.0], 1.0);
        let density = 1.0;
        let i = cs.inertia_tensor(density);
        let vol_sphere = (4.0 / 3.0) * PI;
        let mass_each = density * vol_sphere;
        let i_sphere_y = 2.0 / 5.0 * mass_each;
        let expected_iyy = 2.0 * (i_sphere_y + mass_each * d * d);
        assert!(
            (i[1][1] - expected_iyy).abs() < 1e-6,
            "I_yy={}, expected={expected_iyy}",
            i[1][1]
        );
    }
    #[test]
    fn test_contains_point_capsule_inside_cylinder() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(
            cs.contains_point([0.5, 1.0, 0.0]),
            "should be inside capsule cylinder"
        );
    }
    #[test]
    fn test_contains_point_capsule_inside_hemisphere() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(
            cs.contains_point([0.0, 2.5, 0.0]),
            "should be inside top hemisphere"
        );
    }
    #[test]
    fn test_contains_point_capsule_outside() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(
            !cs.contains_point([0.0, 5.0, 0.0]),
            "should be outside capsule"
        );
    }
    #[test]
    fn test_contains_point_multiple_shapes_uses_union() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        assert!(
            cs.contains_point([0.5, 0.0, 0.0]),
            "should be inside first sphere"
        );
        assert!(
            cs.contains_point([10.5, 0.0, 0.0]),
            "should be inside second sphere"
        );
        assert!(
            !cs.contains_point([5.0, 0.0, 0.0]),
            "should be outside both spheres"
        );
    }
    #[test]
    fn test_closest_point_on_box_surface() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let (cp, _idx) = cs.closest_point([5.0, 0.0, 0.0]);
        assert!(
            (cp[0] - 1.0).abs() < 1e-9,
            "closest x should be at surface 1.0, got {}",
            cp[0]
        );
    }
    #[test]
    fn test_closest_point_on_capsule_cylinder_part() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([0.0, 0.0, 0.0], 1.0, 2.0);
        let (cp, _idx) = cs.closest_point([5.0, 0.0, 0.0]);
        assert!(
            (cp[0] - 1.0).abs() < 1e-9,
            "closest x on capsule should be 1.0, got {}",
            cp[0]
        );
    }
    #[test]
    fn test_merge_preserves_all_shapes() {
        let mut cs1 = CompoundShape::new();
        cs1.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs1.add_box([2.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        let mut cs2 = CompoundShape::new();
        cs2.add_capsule([4.0, 0.0, 0.0], 0.5, 1.0);
        let merged = cs1.merge_with(&cs2);
        assert_eq!(merged.child_count(), 3);
        let vol = merged.total_volume();
        assert!(vol > cs1.total_volume(), "merged volume should be larger");
    }
    #[test]
    fn test_clear_then_add_works() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([5.0, 0.0, 0.0], 2.0);
        cs.clear();
        assert!(cs.is_empty());
        cs.add_sphere([1.0, 0.0, 0.0], 0.5);
        assert_eq!(cs.child_count(), 1);
    }
    #[test]
    fn test_scale_capsule() {
        let mut cs = CompoundShape::new();
        cs.add_capsule([1.0, 0.0, 0.0], 1.0, 2.0);
        cs.scale(3.0);
        match cs.children[0].shape_kind {
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => {
                assert!((radius - 3.0).abs() < 1e-10);
                assert!((half_height - 6.0).abs() < 1e-10);
            }
            _ => panic!("expected Capsule"),
        }
        assert!((cs.children[0].center[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_translate_multiple_children() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([1.0, 2.0, 3.0], 1.0);
        cs.add_box([4.0, 5.0, 6.0], [1.0, 1.0, 1.0]);
        cs.translate([10.0, 0.0, -5.0]);
        assert!((cs.children[0].center[0] - 11.0).abs() < 1e-10);
        assert!((cs.children[0].center[2] - (-2.0)).abs() < 1e-10);
        assert!((cs.children[1].center[0] - 14.0).abs() < 1e-10);
    }
    #[test]
    fn test_compound_ex_volume_two_shapes() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::identity(), 1.0);
        cs.add_box(
            LocalTransform::from_translation([5.0, 0.0, 0.0]),
            [1.0, 1.0, 1.0],
        );
        let v_sphere = (4.0 / 3.0) * PI;
        let v_box = 8.0;
        assert!(
            (cs.volume() - (v_sphere + v_box)).abs() < 1e-9,
            "total volume mismatch"
        );
    }
    #[test]
    fn test_compound_ex_ray_cast_misses() {
        let mut cs = CompoundShapeEx::new();
        cs.add_sphere(LocalTransform::from_translation([0.0, 10.0, 0.0]), 1.0);
        let hit = cs.ray_cast([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_none(), "ray along X should not hit sphere at y=10");
    }
    #[test]
    fn test_compound_ex_contains_point_capsule() {
        let mut cs = CompoundShapeEx::new();
        cs.add_capsule(LocalTransform::from_translation([0.0, 5.0, 0.0]), 1.0, 3.0);
        assert!(
            cs.contains_point([0.0, 5.0, 0.0]),
            "center of capsule should be inside"
        );
        assert!(
            !cs.contains_point([0.0, 0.0, 0.0]),
            "origin should be outside"
        );
    }
    #[test]
    fn test_local_transform_direction_roundtrip() {
        let t = LocalTransform {
            translation: [3.0, -1.0, 2.0],
            rot: [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        };
        let v = [1.0, 0.0, 0.0];
        let world = t.local_to_world_dir(v);
        let back = t.world_to_local_dir(world);
        for i in 0..3 {
            assert!(
                (back[i] - v[i]).abs() < 1e-10,
                "dir roundtrip failed at index {i}"
            );
        }
    }
    #[test]
    fn test_overlaps_sphere_at_boundary() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(
            cs.overlaps_sphere([1.5, 0.0, 0.0], 1.0),
            "overlapping spheres should be detected"
        );
    }
    #[test]
    fn test_sphere_overlaps_aabb_boundary() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(cs.sphere_overlaps_aabb([2.5, 0.0, 0.0], 1.5));
    }
    #[test]
    fn test_penetration_depth_box_penetrates() {
        let mut cs = CompoundShape::new();
        cs.add_box([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let result = cs.penetration_depth_sphere([1.0, 0.0, 0.0], 2.0);
        assert!(result.is_some(), "sphere inside box should penetrate");
        let (depth, _idx) = result.unwrap();
        assert!(depth < 0.0, "penetration depth should be negative");
    }
    #[test]
    fn test_closest_point_with_dist2_multiple_children() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_sphere([10.0, 0.0, 0.0], 1.0);
        let (_cp, _d2, idx) = cs.closest_point_with_dist2([3.0, 0.0, 0.0]);
        assert_eq!(idx, 0, "child 0 is closer");
    }
    #[test]
    fn test_inertia_tensor_is_symmetric_three_mixed_shapes() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([1.0, 0.0, 0.0], 1.0);
        cs.add_box([-1.0, 2.0, 0.0], [0.5, 1.0, 0.5]);
        cs.add_capsule([0.0, -3.0, 1.0], 0.8, 1.5);
        let i = cs.inertia_tensor(2.0);
        assert!((i[0][1] - i[1][0]).abs() < 1e-9, "[0][1] vs [1][0]");
        assert!((i[0][2] - i[2][0]).abs() < 1e-9, "[0][2] vs [2][0]");
        assert!((i[1][2] - i[2][1]).abs() < 1e-9, "[1][2] vs [2][1]");
    }
    #[test]
    fn test_inertia_tensor_diagonal_positive() {
        let mut cs = CompoundShape::new();
        cs.add_sphere([0.0, 0.0, 0.0], 1.0);
        cs.add_box([5.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let i = cs.inertia_tensor(1.0);
        assert!(i[0][0] > 0.0, "I_xx should be positive");
        assert!(i[1][1] > 0.0, "I_yy should be positive");
        assert!(i[2][2] > 0.0, "I_zz should be positive");
    }
}
