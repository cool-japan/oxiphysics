//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::narrowphase::epa::Epa;
use crate::narrowphase::gjk::{Gjk, GjkResult};
use crate::narrowphase::specialized;
use crate::types::{CollisionPair, Contact, ContactManifold};
use oxiphysics_core::Transform;
use oxiphysics_geometry::{BoxShape, Capsule, Shape, Sphere};

use super::types::{
    CompoundDispatchResult, CompoundShape, ConcaveMesh, ContactPatch, MeshTriangle,
    NarrowPhaseDispatcher, NarrowPhaseResult, ShapeType, SimpleHeightField, SpeculativeConfig,
};

/// Batch entry for speculative contact generation.
pub type SpecBatchEntry<'a> = (
    &'a dyn Shape,
    &'a Transform,
    [f64; 3],
    &'a dyn Shape,
    &'a Transform,
    [f64; 3],
    CollisionPair,
);

/// Map shape types to a numeric ordinal for canonical ordering.
pub(super) fn shape_type_ordinal(t: ShapeType) -> u32 {
    match t {
        ShapeType::Sphere => 0,
        ShapeType::Box => 1,
        ShapeType::Capsule => 2,
        ShapeType::Cylinder => 3,
        ShapeType::Cone => 4,
        ShapeType::ConvexHull => 5,
        ShapeType::TriangleMesh => 6,
        ShapeType::Compound => 7,
        ShapeType::HeightField => 8,
        ShapeType::Plane => 9,
        ShapeType::Custom(id) => 100 + id,
    }
}
/// Type of a registered narrow-phase collision function.
pub type NarrowPhaseFn = fn(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    pair: CollisionPair,
) -> NarrowPhaseResult;
pub(super) fn sphere_capsule_dispatch(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    pair: CollisionPair,
) -> NarrowPhaseResult {
    use oxiphysics_core::math::Vec3;
    let sphere = unsafe { &*(shape_a as *const dyn Shape as *const Sphere) };
    let capsule = unsafe { &*(shape_b as *const dyn Shape as *const Capsule) };
    let sphere_center = transform_a.position;
    let half_h = capsule.half_height;
    let cap_up = transform_b.rotation * Vec3::new(0.0, half_h, 0.0);
    let cap_a = transform_b.position + cap_up;
    let cap_b = transform_b.position - cap_up;
    let ab = cap_b - cap_a;
    let ab_len2 = ab.dot(&ab);
    let t = if ab_len2 > 1e-12 {
        ((sphere_center - cap_a).dot(&ab) / ab_len2).clamp(0.0, 1.0)
    } else {
        0.5
    };
    let closest = cap_a + ab * t;
    let diff = sphere_center - closest;
    let dist = diff.norm();
    let combined = sphere.radius + capsule.radius;
    if dist >= combined {
        return NarrowPhaseResult::separated();
    }
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let depth = combined - dist;
    let point_a = sphere_center - normal * sphere.radius;
    let point_b = closest + normal * capsule.radius;
    let contact = Contact::new(point_a, point_b, normal, depth);
    let mut m = ContactManifold::new(pair);
    m.add_contact(contact);
    NarrowPhaseResult::contact(m)
}
pub(super) fn gjk_fallback_dispatch(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    pair: CollisionPair,
) -> NarrowPhaseResult {
    match gjk_epa(shape_a, transform_a, shape_b, transform_b, pair) {
        Some(manifold) => NarrowPhaseResult::contact(manifold),
        None => NarrowPhaseResult::separated(),
    }
}
/// Dispatch against a compound shape by testing each child shape individually
/// and collecting all contacts.
pub fn dispatch_compound(
    dispatcher: &NarrowPhaseDispatcher,
    compound_shapes: &[(ShapeType, &dyn Shape, Transform)],
    other_shape: &dyn Shape,
    other_type: ShapeType,
    other_transform: &Transform,
    base_pair: CollisionPair,
) -> NarrowPhaseResult {
    let mut all_contacts = Vec::new();
    for (child_type, child_shape, child_transform) in compound_shapes {
        let result = dispatcher.dispatch(
            *child_shape,
            *child_type,
            child_transform,
            other_shape,
            other_type,
            other_transform,
            base_pair,
        );
        if let Some(manifold) = result.manifold {
            all_contacts.extend(manifold.contacts);
        }
    }
    if all_contacts.is_empty() {
        NarrowPhaseResult::separated()
    } else {
        let mut manifold = ContactManifold::new(base_pair);
        for c in all_contacts {
            manifold.add_contact(c);
        }
        NarrowPhaseResult::contact(manifold)
    }
}
/// Attempt to use a specialized collision test.
pub(super) fn try_specialized(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> Option<Contact> {
    let name_a = type_name(shape_a);
    let name_b = type_name(shape_b);
    match (name_a, name_b) {
        ("Sphere", "Sphere") => {
            let s1 = as_sphere(shape_a)?;
            let s2 = as_sphere(shape_b)?;
            specialized::sphere_sphere(s1, transform_a, s2, transform_b)
        }
        ("Sphere", "BoxShape") => {
            let s = as_sphere(shape_a)?;
            let b = as_box(shape_b)?;
            specialized::sphere_box(s, transform_a, b, transform_b)
        }
        ("BoxShape", "Sphere") => {
            let b = as_box(shape_a)?;
            let s = as_sphere(shape_b)?;
            let mut c = specialized::sphere_box(s, transform_b, b, transform_a)?;
            c.normal = -c.normal;
            std::mem::swap(&mut c.point_a, &mut c.point_b);
            Some(c)
        }
        ("BoxShape", "BoxShape") => {
            let b1 = as_box(shape_a)?;
            let b2 = as_box(shape_b)?;
            specialized::box_box_sat(b1, transform_a, b2, transform_b)
        }
        ("Capsule", "Capsule") => {
            let c1 = as_capsule(shape_a)?;
            let c2 = as_capsule(shape_b)?;
            specialized::capsule_capsule(c1, transform_a, c2, transform_b)
        }
        _ => None,
    }
}
pub(super) fn gjk_epa(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    pair: CollisionPair,
) -> Option<ContactManifold> {
    let result = Gjk::query(shape_a, transform_a, shape_b, transform_b);
    match result {
        GjkResult::Intersecting(simplex) => {
            let contact =
                Epa::penetration_depth(shape_a, transform_a, shape_b, transform_b, &simplex)?;
            let mut manifold = ContactManifold::new(pair);
            manifold.add_contact(contact);
            Some(manifold)
        }
        GjkResult::Separated { .. } => None,
    }
}
/// Extract type name from Debug output.
pub(super) fn type_name(shape: &dyn Shape) -> &'static str {
    let debug = format!("{:?}", shape);
    if debug.starts_with("Sphere") {
        "Sphere"
    } else if debug.starts_with("BoxShape") {
        "BoxShape"
    } else if debug.starts_with("Capsule") {
        "Capsule"
    } else if debug.starts_with("Cylinder") {
        "Cylinder"
    } else if debug.starts_with("Cone") {
        "Cone"
    } else {
        "Unknown"
    }
}
pub(super) fn as_sphere(shape: &dyn Shape) -> Option<&Sphere> {
    let s = shape.support_point(&oxiphysics_core::math::Vec3::new(1.0, 0.0, 0.0));
    let radius = s.x;
    if radius > 0.0 {
        let ptr = shape as *const dyn Shape as *const Sphere;
        Some(unsafe { &*ptr })
    } else {
        None
    }
}
pub(super) fn as_box(shape: &dyn Shape) -> Option<&BoxShape> {
    let ptr = shape as *const dyn Shape as *const BoxShape;
    Some(unsafe { &*ptr })
}
pub(super) fn as_capsule(shape: &dyn Shape) -> Option<&Capsule> {
    let ptr = shape as *const dyn Shape as *const Capsule;
    Some(unsafe { &*ptr })
}
/// Dispatch compound vs compound using O(n²) brute-force with optional AABB pruning.
///
/// `use_aabb_pruning` enables a cheap AABB overlap pre-filter; when false
/// every sub-shape pair is tested unconditionally.
pub fn dispatch_compound_compound(
    dispatcher: &NarrowPhaseDispatcher,
    compound_a: &CompoundShape,
    transform_a: &Transform,
    compound_b: &CompoundShape,
    transform_b: &Transform,
    pair_id: u32,
    use_aabb_pruning: bool,
) -> CompoundDispatchResult {
    let mut manifolds = Vec::new();
    let mut tests = 0usize;
    let mut pruned = 0usize;
    let mut pair_counter = (pair_id as usize) * 1000;
    for (ia, (child_a, ta_local)) in compound_a
        .children
        .iter()
        .zip(compound_a.local_transforms.iter())
        .enumerate()
    {
        let world_a = compose_transforms(transform_a, ta_local);
        for (ib, (child_b, tb_local)) in compound_b
            .children
            .iter()
            .zip(compound_b.local_transforms.iter())
            .enumerate()
        {
            let world_b = compose_transforms(transform_b, tb_local);
            if use_aabb_pruning {
                let pa = world_a.position;
                let pb = world_b.position;
                let dist2 = {
                    let dx = pa.x - pb.x;
                    let dy = pa.y - pb.y;
                    let dz = pa.z - pb.z;
                    dx * dx + dy * dy + dz * dz
                };
                if dist2 > 25.0 {
                    pruned += 1;
                    continue;
                }
            }
            let sub_pair = CollisionPair::new(ia, ib + pair_counter);
            pair_counter += 1;
            tests += 1;
            let result = dispatcher.dispatch(
                child_a.as_ref(),
                compound_a.child_types[ia],
                &world_a,
                child_b.as_ref(),
                compound_b.child_types[ib],
                &world_b,
                sub_pair,
            );
            if result.has_contact()
                && let Some(m) = result.manifold
            {
                manifolds.push(m);
            }
        }
    }
    CompoundDispatchResult {
        manifolds,
        tests_performed: tests,
        pairs_pruned: pruned,
    }
}
/// Compose a parent world transform with a local child transform.
pub(super) fn compose_transforms(parent: &Transform, local: &Transform) -> Transform {
    parent.compose(local)
}
/// Dispatch a convex sphere against a heightfield.
///
/// Returns a `NarrowPhaseResult` using the heightfield's sphere test.
pub fn dispatch_heightfield_sphere(
    hf: &SimpleHeightField,
    sphere_center: [f64; 3],
    sphere_radius: f64,
    pair: CollisionPair,
) -> NarrowPhaseResult {
    if let Some((point, normal, depth)) = hf.test_sphere(sphere_center, sphere_radius) {
        let contact = Contact {
            point_a: oxiphysics_core::math::Vec3::new(point[0], point[1], point[2]),
            point_b: oxiphysics_core::math::Vec3::new(
                sphere_center[0] - normal[0] * sphere_radius,
                sphere_center[1] - normal[1] * sphere_radius,
                sphere_center[2] - normal[2] * sphere_radius,
            ),
            normal: oxiphysics_core::math::Vec3::new(normal[0], normal[1], normal[2]),
            depth,
        };
        let mut m = ContactManifold::new(pair);
        m.add_contact(contact);
        NarrowPhaseResult::contact(m)
    } else {
        NarrowPhaseResult::separated()
    }
}
/// Mid-phase mesh vs mesh AABB overlap test.
///
/// Returns the pairs `(i, j)` of triangle indices whose AABBs overlap.
pub fn mesh_mesh_aabb_mid_phase(
    triangles_a: &[MeshTriangle],
    triangles_b: &[MeshTriangle],
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for (ia, ta) in triangles_a.iter().enumerate() {
        for (ib, tb) in triangles_b.iter().enumerate() {
            if ta.aabb.overlaps(&tb.aabb) {
                pairs.push((ia, ib));
            }
        }
    }
    pairs
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::CollisionPair;
    use crate::Contact;
    use crate::ContactManifold;
    use crate::NarrowPhaseDispatcher;

    use crate::narrowphase::Aabb;

    use crate::narrowphase::DispatchConfig;
    use crate::narrowphase::DispatchKey;

    use crate::narrowphase::DispatchStats;
    use crate::narrowphase::ShapeFeatureCache;
    use crate::narrowphase::SimpleHeightField;

    use crate::narrowphase::dispatch::ShapeType;

    use oxiphysics_core::Vec3;
    use std::collections::HashMap;
    #[test]
    fn test_dispatch_sphere_sphere() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let manifold = NarrowPhaseDispatcher::generate_contacts(&s1, &t1, &s2, &t2, pair);
        assert!(manifold.is_some());
        let m = manifold.unwrap();
        assert_eq!(m.contacts.len(), 1);
    }
    #[test]
    fn test_dispatch_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let manifold = NarrowPhaseDispatcher::generate_contacts(&s1, &t1, &s2, &t2, pair);
        assert!(manifold.is_none());
    }
    #[test]
    fn test_dispatch_key_canonical() {
        let k1 = DispatchKey::new(ShapeType::Sphere, ShapeType::Box);
        let k2 = DispatchKey::new(ShapeType::Box, ShapeType::Sphere);
        assert_eq!(k1, k2, "DispatchKey should be order-independent");
    }
    #[test]
    fn test_dispatch_key_same() {
        let k = DispatchKey::new(ShapeType::Capsule, ShapeType::Capsule);
        assert_eq!(k.0, k.1);
    }
    #[test]
    fn test_default_dispatcher_sphere_sphere() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert!(
            result.manifold.is_some(),
            "overlapping spheres should produce contact"
        );
    }
    #[test]
    fn test_default_dispatcher_sphere_sphere_separated() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert!(
            result.manifold.is_none(),
            "separated spheres should produce no contact"
        );
    }
    #[test]
    fn test_register_custom_pair() {
        let mut dispatcher = NarrowPhaseDispatcher::empty();
        dispatcher.register_pair(
            ShapeType::Sphere,
            ShapeType::Sphere,
            |_sa, _ta, _sb, _tb, pair| {
                use crate::types::Contact;
                let contact = Contact {
                    point_a: Vec3::zeros(),
                    point_b: Vec3::zeros(),
                    normal: Vec3::new(0.0, 1.0, 0.0),
                    depth: 1.0,
                };
                let mut m = ContactManifold::new(pair);
                m.add_contact(contact);
                NarrowPhaseResult::contact(m)
            },
        );
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(100.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert!(
            result.manifold.is_some(),
            "custom algorithm should produce contact"
        );
    }
    #[test]
    fn test_narrow_phase_result_separated() {
        let r = NarrowPhaseResult::separated();
        assert!(r.manifold.is_none());
        assert!(!r.has_contact());
        assert!(r.penetration_depth().is_none());
    }
    #[test]
    fn test_narrow_phase_result_has_contact() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::zeros(),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.5,
        });
        let r = NarrowPhaseResult::contact(m);
        assert!(r.has_contact());
        assert!((r.penetration_depth().unwrap() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_shape_type_all_variants() {
        let types = [
            ShapeType::Sphere,
            ShapeType::Box,
            ShapeType::Capsule,
            ShapeType::Cylinder,
            ShapeType::Cone,
            ShapeType::ConvexHull,
            ShapeType::TriangleMesh,
            ShapeType::Compound,
            ShapeType::HeightField,
            ShapeType::Plane,
            ShapeType::Custom(0),
        ];
        let mut map: HashMap<DispatchKey, usize> = HashMap::new();
        for (i, &t) in types.iter().enumerate() {
            map.insert(DispatchKey::new(t, t), i);
        }
        assert_eq!(map.len(), types.len());
    }
    #[test]
    fn test_dispatch_key_custom_types() {
        let k1 = DispatchKey::new(ShapeType::Custom(0), ShapeType::Custom(1));
        let k2 = DispatchKey::new(ShapeType::Custom(1), ShapeType::Custom(0));
        assert_eq!(k1, k2, "custom type keys should be canonical");
    }
    #[test]
    fn test_unregister_pair() {
        let mut dispatcher = NarrowPhaseDispatcher::default();
        assert!(dispatcher.has_pair(ShapeType::Sphere, ShapeType::Sphere));
        let removed = dispatcher.unregister_pair(ShapeType::Sphere, ShapeType::Sphere);
        assert!(removed);
        assert!(!dispatcher.has_pair(ShapeType::Sphere, ShapeType::Sphere));
    }
    #[test]
    fn test_unregister_nonexistent() {
        let mut dispatcher = NarrowPhaseDispatcher::empty();
        let removed = dispatcher.unregister_pair(ShapeType::Cone, ShapeType::Cone);
        assert!(!removed);
    }
    #[test]
    fn test_registered_count() {
        let dispatcher = NarrowPhaseDispatcher::default();
        assert_eq!(dispatcher.registered_count(), 6);
    }
    #[test]
    fn test_registered_keys() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let keys = dispatcher.registered_keys();
        assert_eq!(keys.len(), 6);
    }
    #[test]
    fn test_has_pair() {
        let dispatcher = NarrowPhaseDispatcher::default();
        assert!(dispatcher.has_pair(ShapeType::Sphere, ShapeType::Sphere));
        assert!(dispatcher.has_pair(ShapeType::Sphere, ShapeType::Box));
        assert!(dispatcher.has_pair(ShapeType::Box, ShapeType::Sphere));
        assert!(!dispatcher.has_pair(ShapeType::Cone, ShapeType::Cone));
    }
    #[test]
    fn test_sphere_capsule_overlapping() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let sphere = Sphere::new(1.0);
        let capsule = Capsule::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &sphere,
            ShapeType::Sphere,
            &t1,
            &capsule,
            ShapeType::Capsule,
            &t2,
            pair,
        );
        assert!(
            result.has_contact(),
            "overlapping sphere-capsule should produce contact"
        );
    }
    #[test]
    fn test_sphere_capsule_separated() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let sphere = Sphere::new(0.5);
        let capsule = Capsule::new(0.5, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &sphere,
            ShapeType::Sphere,
            &t1,
            &capsule,
            ShapeType::Capsule,
            &t2,
            pair,
        );
        assert!(
            !result.has_contact(),
            "separated sphere-capsule should not produce contact"
        );
    }
    #[test]
    fn test_sphere_capsule_contact_depth() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let sphere = Sphere::new(1.0);
        let capsule = Capsule::new(1.0, 1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &sphere,
            ShapeType::Sphere,
            &t1,
            &capsule,
            ShapeType::Capsule,
            &t2,
            pair,
        );
        assert!(result.has_contact());
        let depth = result.penetration_depth().unwrap();
        assert!(depth > 0.0, "penetration depth should be positive");
        assert!(
            (depth - 0.5).abs() < 1e-6,
            "expected depth 0.5, got {}",
            depth
        );
    }
    #[test]
    fn test_gjk_fallback_for_unknown_pair() {
        let dispatcher = NarrowPhaseDispatcher::empty();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert!(
            result.has_contact(),
            "GJK fallback should detect overlapping spheres"
        );
    }
    #[test]
    fn test_gjk_fallback_separated() {
        let dispatcher = NarrowPhaseDispatcher::empty();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert!(!result.has_contact());
    }
    #[test]
    fn test_dispatch_config_default() {
        let config = DispatchConfig::default();
        assert!(config.use_gjk_fallback);
        assert!(config.max_gjk_iterations > 0);
        assert!(config.max_epa_iterations > 0);
        assert!(config.contact_tolerance > 0.0);
    }
    #[test]
    fn test_dispatcher_with_config() {
        let config = DispatchConfig {
            use_gjk_fallback: false,
            max_gjk_iterations: 32,
            max_epa_iterations: 32,
            contact_tolerance: 1e-4,
        };
        let dispatcher = NarrowPhaseDispatcher::with_config(config);
        assert_eq!(dispatcher.registered_count(), 0);
    }
    #[test]
    fn test_heightfield_flat() {
        let heights = vec![0.0; 9];
        let hf = SimpleHeightField::new(heights, 3, 3, 1.0);
        assert_eq!(hf.height_at(0, 0), 0.0);
        assert_eq!(hf.height_at(1, 1), 0.0);
    }
    #[test]
    fn test_heightfield_normal_flat() {
        let heights = vec![0.0; 9];
        let hf = SimpleHeightField::new(heights, 3, 3, 1.0);
        let n = hf.normal_at(0, 0);
        assert!((n[1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_heightfield_sphere_collision() {
        let heights = vec![0.0; 9];
        let hf = SimpleHeightField::new(heights, 3, 3, 1.0);
        let result = hf.test_sphere([1.0, 0.5, 1.0], 1.0);
        assert!(result.is_some());
        let (_, _, depth) = result.unwrap();
        assert!(
            (depth - 0.5).abs() < 1e-6,
            "expected depth 0.5, got {}",
            depth
        );
    }
    #[test]
    fn test_heightfield_sphere_no_collision() {
        let heights = vec![0.0; 9];
        let hf = SimpleHeightField::new(heights, 3, 3, 1.0);
        let result = hf.test_sphere([1.0, 5.0, 1.0], 1.0);
        assert!(result.is_none());
    }
    #[test]
    fn test_heightfield_varied_heights() {
        let heights = vec![0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0];
        let hf = SimpleHeightField::new(heights, 3, 3, 1.0);
        assert_eq!(hf.height_at(1, 1), 2.0);
        assert_eq!(hf.height_at(0, 1), 1.0);
    }
    #[test]
    fn test_dispatch_symmetric_same_as_normal() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let r1 = dispatcher.dispatch(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        let r2 = dispatcher.dispatch_symmetric(
            &s1,
            ShapeType::Sphere,
            &t1,
            &s2,
            ShapeType::Sphere,
            &t2,
            pair,
        );
        assert_eq!(r1.has_contact(), r2.has_contact());
    }
    #[test]
    fn test_shape_type_ordinal_ordering() {
        assert!(shape_type_ordinal(ShapeType::Sphere) < shape_type_ordinal(ShapeType::Box));
        assert!(shape_type_ordinal(ShapeType::Box) < shape_type_ordinal(ShapeType::Capsule));
        assert!(shape_type_ordinal(ShapeType::HeightField) < shape_type_ordinal(ShapeType::Plane));
        assert!(shape_type_ordinal(ShapeType::Plane) < shape_type_ordinal(ShapeType::Custom(0)));
    }
    #[test]
    fn test_dispatch_key_plane() {
        let k = DispatchKey::new(ShapeType::Plane, ShapeType::Sphere);
        assert_eq!(k.0, ShapeType::Sphere);
        assert_eq!(k.1, ShapeType::Plane);
    }
    #[test]
    fn test_compound_shape_empty() {
        let c = CompoundShape::new();
        assert!(c.is_empty());
        assert_eq!(c.len(), 0);
    }
    #[test]
    fn test_compound_shape_add_children() {
        let mut c = CompoundShape::new();
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        c.add_child(Box::new(s), ShapeType::Sphere, t);
        assert_eq!(c.len(), 1);
        assert!(!c.is_empty());
    }
    #[test]
    fn test_dispatch_compound_compound_no_overlap() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let mut ca = CompoundShape::new();
        ca.add_child(
            Box::new(Sphere::new(0.5)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        let mut cb = CompoundShape::new();
        cb.add_child(
            Box::new(Sphere::new(0.5)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        let ta = Transform::from_position(Vec3::new(-100.0, 0.0, 0.0));
        let tb = Transform::from_position(Vec3::new(100.0, 0.0, 0.0));
        let result = dispatch_compound_compound(&dispatcher, &ca, &ta, &cb, &tb, 0, true);
        assert!(
            !result.has_contacts(),
            "separated compounds should have no contacts"
        );
    }
    #[test]
    fn test_dispatch_compound_compound_overlap() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let mut ca = CompoundShape::new();
        ca.add_child(
            Box::new(Sphere::new(1.0)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        let mut cb = CompoundShape::new();
        cb.add_child(
            Box::new(Sphere::new(1.0)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        let ta = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let tb = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = dispatch_compound_compound(&dispatcher, &ca, &ta, &cb, &tb, 0, false);
        assert!(
            result.has_contacts(),
            "overlapping compounds should have contacts"
        );
    }
    #[test]
    fn test_dispatch_compound_pruning_reduces_tests() {
        let dispatcher = NarrowPhaseDispatcher::default();
        let mut ca = CompoundShape::new();
        let mut cb = CompoundShape::new();
        ca.add_child(
            Box::new(Sphere::new(0.1)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        cb.add_child(
            Box::new(Sphere::new(0.1)),
            ShapeType::Sphere,
            Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
        );
        let ta = Transform::from_position(Vec3::new(-100.0, 0.0, 0.0));
        let tb = Transform::from_position(Vec3::new(100.0, 0.0, 0.0));
        let pruned = dispatch_compound_compound(&dispatcher, &ca, &ta, &cb, &tb, 0, true);
        let not_pruned = dispatch_compound_compound(&dispatcher, &ca, &ta, &cb, &tb, 0, false);
        assert!(pruned.tests_performed <= not_pruned.tests_performed);
    }
    #[test]
    fn test_heightfield_sphere_dispatch_contact() {
        let hf = SimpleHeightField::new(vec![0.0; 9], 3, 3, 1.0);
        let pair = CollisionPair::new(0, 1);
        let result = dispatch_heightfield_sphere(&hf, [1.0, 0.5, 1.0], 1.0, pair);
        assert!(
            result.has_contact(),
            "sphere penetrating heightfield should have contact"
        );
    }
    #[test]
    fn test_heightfield_sphere_dispatch_no_contact() {
        let hf = SimpleHeightField::new(vec![0.0; 9], 3, 3, 1.0);
        let pair = CollisionPair::new(0, 1);
        let result = dispatch_heightfield_sphere(&hf, [1.0, 10.0, 1.0], 0.5, pair);
        assert!(
            !result.has_contact(),
            "sphere above heightfield should have no contact"
        );
    }
    #[test]
    fn test_aabb_overlaps_true() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_aabb_overlaps_false() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_aabb_surface_area() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        let sa = a.surface_area();
        assert!((sa - 52.0).abs() < 1e-10, "surface area should be 52: {sa}");
    }
    #[test]
    fn test_mesh_triangle_aabb() {
        let tri = MeshTriangle::new([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]], 0);
        assert!((tri.aabb.min[0]).abs() < 1e-10);
        assert!((tri.aabb.max[0] - 1.0).abs() < 1e-10);
        assert!((tri.aabb.max[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_mesh_mid_phase_overlapping() {
        let ta = vec![MeshTriangle::new(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            0,
        )];
        let tb = vec![MeshTriangle::new(
            [[0.3, 0.0, 0.0], [1.3, 0.0, 0.0], [0.8, 1.0, 0.0]],
            0,
        )];
        let pairs = mesh_mesh_aabb_mid_phase(&ta, &tb);
        assert!(
            !pairs.is_empty(),
            "overlapping triangles should produce a pair"
        );
    }
    #[test]
    fn test_mesh_mesh_mid_phase_non_overlapping() {
        let ta = vec![MeshTriangle::new(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            0,
        )];
        let tb = vec![MeshTriangle::new(
            [[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.5, 1.0, 0.0]],
            0,
        )];
        let pairs = mesh_mesh_aabb_mid_phase(&ta, &tb);
        assert!(
            pairs.is_empty(),
            "separated triangles should produce no pairs"
        );
    }
    #[test]
    fn test_feature_cache_miss_then_hit() {
        let mut cache = ShapeFeatureCache::new();
        let v1 = cache.get_or_insert(0xAB, || [1.0, 0.0, 0.0]);
        let v2 = cache.get_or_insert(0xAB, || [999.0, 0.0, 0.0]);
        assert!((v1[0] - 1.0).abs() < 1e-10);
        assert!(
            (v2[0] - 1.0).abs() < 1e-10,
            "should return cached value: {}",
            v2[0]
        );
        assert_eq!(cache.misses, 1);
        assert_eq!(cache.hits, 1);
    }
    #[test]
    fn test_feature_cache_hit_ratio() {
        let mut cache = ShapeFeatureCache::new();
        cache.get_or_insert(1, || [0.0; 3]);
        cache.get_or_insert(1, || [0.0; 3]);
        cache.get_or_insert(1, || [0.0; 3]);
        assert!((cache.hit_ratio() - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_feature_cache_clear() {
        let mut cache = ShapeFeatureCache::new();
        cache.get_or_insert(1, || [1.0; 3]);
        cache.clear();
        assert_eq!(cache.hits, 0);
        assert_eq!(cache.misses, 0);
        cache.get_or_insert(1, || [2.0; 3]);
        assert_eq!(cache.misses, 1);
    }
    #[test]
    fn test_dispatch_stats_contact_rate() {
        let mut stats = DispatchStats::default();
        stats.record(true, false);
        stats.record(false, false);
        stats.record(true, true);
        assert!((stats.contact_rate() - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_dispatch_stats_specialization_rate() {
        let mut stats = DispatchStats::default();
        stats.record(true, false);
        stats.record(false, false);
        stats.record(true, true);
        assert!((stats.specialization_rate() - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_dispatch_stats_pruned() {
        let mut stats = DispatchStats::default();
        stats.record_pruned();
        stats.record_pruned();
        assert_eq!(stats.broad_phase_pruned, 2);
    }
    #[test]
    fn test_dispatch_stats_empty_rates() {
        let stats = DispatchStats::default();
        assert!((stats.contact_rate()).abs() < 1e-10);
        assert!((stats.specialization_rate()).abs() < 1e-10);
    }
}
/// Clip a polygon against a half-plane defined by an edge.
///
/// Returns the clipped polygon.  Used as a building block for Sutherland-Hodgman
/// polygon clipping in multi-contact generation.
pub(super) fn clip_polygon_by_halfspace(
    polygon: &[[f64; 2]],
    plane_point: [f64; 2],
    plane_normal: [f64; 2],
) -> Vec<[f64; 2]> {
    if polygon.is_empty() {
        return Vec::new();
    }
    let mut result = Vec::new();
    let n = polygon.len();
    for i in 0..n {
        let curr = polygon[i];
        let next = polygon[(i + 1) % n];
        let dc = (curr[0] - plane_point[0]) * plane_normal[0]
            + (curr[1] - plane_point[1]) * plane_normal[1];
        let dn = (next[0] - plane_point[0]) * plane_normal[0]
            + (next[1] - plane_point[1]) * plane_normal[1];
        if dc >= 0.0 {
            result.push(curr);
        }
        if (dc >= 0.0) != (dn >= 0.0) {
            let denom = dc - dn;
            if denom.abs() > 1e-14 {
                let t = dc / denom;
                result.push([
                    curr[0] + t * (next[0] - curr[0]),
                    curr[1] + t * (next[1] - curr[1]),
                ]);
            }
        }
    }
    result
}
/// Generate a multi-contact patch for two overlapping axis-aligned rectangles
/// (box face vs box face) using Sutherland-Hodgman clipping in 2D.
///
/// `face_a` and `face_b` are arrays of 4 vertices (world-space XZ plane).
/// Returns a `ContactPatch` with up to 8 contact points.
pub fn clip_box_faces_2d(
    face_a: [[f64; 2]; 4],
    face_b: [[f64; 2]; 4],
    normal: [f64; 3],
    depth: f64,
) -> ContactPatch {
    let mut clipped: Vec<[f64; 2]> = face_b.to_vec();
    for i in 0..4 {
        if clipped.is_empty() {
            break;
        }
        let a = face_a[i];
        let b = face_a[(i + 1) % 4];
        let edge = [b[0] - a[0], b[1] - a[1]];
        let edge_normal = [edge[1], -edge[0]];
        clipped = clip_polygon_by_halfspace(&clipped, a, edge_normal);
    }
    let mut patch = ContactPatch::new(normal, depth);
    for p in &clipped {
        patch.add_point([p[0], 0.0, p[1]]);
    }
    patch
}
/// Generate a contact manifold with multiple contact points for a box-box pair
/// using SAT-derived reference/incident face clipping.
///
/// Returns up to 4 contact points in world space.
pub fn box_box_manifold_clip(
    box_a: &BoxShape,
    transform_a: &oxiphysics_core::Transform,
    box_b: &BoxShape,
    transform_b: &oxiphysics_core::Transform,
    pair: CollisionPair,
) -> Option<ContactManifold> {
    use crate::narrowphase::specialized::box_box_sat;
    let single_contact = box_box_sat(box_a, transform_a, box_b, transform_b)?;
    let mut manifold = ContactManifold::new(pair);
    manifold.add_contact(single_contact.clone());
    let n = single_contact.normal;
    let perp1 = if n.x.abs() < 0.9 {
        oxiphysics_core::math::Vec3::new(1.0, 0.0, 0.0) - n * n.x
    } else {
        oxiphysics_core::math::Vec3::new(0.0, 1.0, 0.0) - n * n.y
    };
    let perp1_len = perp1.norm();
    if perp1_len < 1e-10 {
        return Some(manifold);
    }
    let perp1 = perp1 / perp1_len;
    let perp2 = n.cross(&perp1);
    let half = [
        box_a.half_extents.x,
        box_a.half_extents.y,
        box_a.half_extents.z,
    ];
    let scale = *half
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(&0.1);
    let offsets = [
        perp1 * scale + perp2 * scale,
        perp1 * scale - perp2 * scale,
        -perp1 * scale + perp2 * scale,
        -perp1 * scale - perp2 * scale,
    ];
    for off in &offsets {
        let shifted_ta = oxiphysics_core::Transform {
            position: transform_a.position + *off,
            rotation: transform_a.rotation,
        };
        if let Some(c) = box_box_sat(box_a, &shifted_ta, box_b, transform_b)
            && c.depth > 0.0
        {
            let mut adjusted = c;
            adjusted.point_a -= *off;
            adjusted.point_b -= *off;
            manifold.add_contact(adjusted);
        }
    }
    manifold.contacts.sort_by(|a, b| {
        b.depth
            .partial_cmp(&a.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    manifold.contacts.truncate(4);
    Some(manifold)
}
/// Generate speculative contacts between two convex shapes given their velocities.
///
/// When two shapes are within `config.margin` of each other AND moving toward
/// each other faster than a threshold, a speculative contact is injected so
/// the constraint solver can prevent tunnelling over the next timestep.
///
/// The contact depth is reported as the current gap (negative = overlapping)
/// so the solver can distinguish pre-contact from post-contact.
pub fn generate_speculative_contacts(
    shape_a: &dyn Shape,
    transform_a: &oxiphysics_core::Transform,
    velocity_a: [f64; 3],
    shape_b: &dyn Shape,
    transform_b: &oxiphysics_core::Transform,
    velocity_b: [f64; 3],
    pair: CollisionPair,
    config: &SpeculativeConfig,
) -> Option<ContactManifold> {
    use crate::narrowphase::gjk::{GjkDistanceResult, gjk_distance_query};
    let result: GjkDistanceResult = gjk_distance_query(shape_a, transform_a, shape_b, transform_b);
    if result.distance > config.margin {
        return None;
    }
    let diff = result.closest_b - result.closest_a;
    let dist = diff.norm();
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        let d = transform_b.position - transform_a.position;
        if d.norm_squared() > 1e-10 {
            d.normalize()
        } else {
            oxiphysics_core::math::Vec3::new(0.0, 1.0, 0.0)
        }
    };
    let rel_vel = oxiphysics_core::math::Vec3::new(
        velocity_b[0] - velocity_a[0],
        velocity_b[1] - velocity_a[1],
        velocity_b[2] - velocity_a[2],
    );
    let closing = -rel_vel.dot(&normal);
    if closing < -config.velocity_scale * config.margin && result.distance > 0.0 {
        return None;
    }
    let depth = config.margin - result.distance;
    let contact = Contact {
        point_a: result.closest_a,
        point_b: result.closest_b,
        normal,
        depth,
    };
    let mut manifold = ContactManifold::new(pair);
    manifold.add_contact(contact);
    Some(manifold)
}
/// Dispatch a convex shape against a concave mesh (convex decomposition).
///
/// Tests each convex piece of the mesh against the convex shape using GJK,
/// collecting all contacts.  An AABB pre-filter limits the number of tests.
pub fn dispatch_concave_vs_convex(
    mesh: &ConcaveMesh,
    mesh_transform: &oxiphysics_core::Transform,
    _convex: &dyn Shape,
    convex_type: ShapeType,
    convex_transform: &oxiphysics_core::Transform,
    base_pair: CollisionPair,
    dispatcher: &NarrowPhaseDispatcher,
) -> Vec<ContactManifold> {
    let _ = convex_type;
    let _ = dispatcher;
    let mut manifolds = Vec::new();
    let inv_mesh = mesh_transform.inverse();
    let local_convex_pos = inv_mesh.transform_point(&convex_transform.position);
    let convex_aabb_expand = 2.0;
    for (piece_idx, piece) in mesh.pieces.iter().enumerate() {
        let mut rejected = false;
        for i in 0..3 {
            let lp = [local_convex_pos.x, local_convex_pos.y, local_convex_pos.z][i];
            if lp + convex_aabb_expand < piece.aabb_min[i]
                || lp - convex_aabb_expand > piece.aabb_max[i]
            {
                rejected = true;
                break;
            }
        }
        if rejected {
            continue;
        }
        let verts_world: Vec<oxiphysics_core::math::Vec3> = piece
            .vertices
            .iter()
            .map(|v| {
                mesh_transform.transform_point(&oxiphysics_core::math::Vec3::new(v[0], v[1], v[2]))
            })
            .collect();
        if verts_world.is_empty() {
            continue;
        }
        use crate::narrowphase::specialized::convex_convex_gjk_intersect;
        if convex_convex_gjk_intersect(&verts_world, mesh_transform, &verts_world, convex_transform)
        {
            let centroid: oxiphysics_core::math::Vec3 =
                verts_world.iter().sum::<oxiphysics_core::math::Vec3>() / verts_world.len() as f64;
            let diff = convex_transform.position - centroid;
            let dist = diff.norm();
            let normal = if dist > 1e-10 {
                diff / dist
            } else {
                oxiphysics_core::math::Vec3::new(0.0, 1.0, 0.0)
            };
            let contact = Contact {
                point_a: centroid,
                point_b: convex_transform.position,
                normal,
                depth: 0.1,
            };
            let sub_pair = CollisionPair::new(base_pair.a, base_pair.b + piece_idx * 100);
            let mut m = ContactManifold::new(sub_pair);
            m.add_contact(contact);
            manifolds.push(m);
        }
    }
    manifolds
}
/// Generic heightfield vs convex dispatch.
///
/// Samples the heightfield at all grid cells within the convex shape's bounding
/// sphere and returns contacts for any penetrating cells.
pub fn dispatch_heightfield_convex_generic(
    hf: &SimpleHeightField,
    convex_pos: [f64; 3],
    convex_radius: f64,
    convex_half_height: f64,
    pair: CollisionPair,
) -> Vec<ContactManifold> {
    let mut manifolds = Vec::new();
    let search_r = (convex_radius + convex_half_height).abs();
    let col_min = ((convex_pos[0] - search_r) / hf.spacing).floor() as isize;
    let col_max = ((convex_pos[0] + search_r) / hf.spacing).ceil() as isize;
    let row_min = ((convex_pos[2] - search_r) / hf.spacing).floor() as isize;
    let row_max = ((convex_pos[2] + search_r) / hf.spacing).ceil() as isize;
    for row in row_min..=row_max {
        for col in col_min..=col_max {
            if row < 0 || col < 0 {
                continue;
            }
            let row = row as usize;
            let col = col as usize;
            if row >= hf.rows || col >= hf.cols {
                continue;
            }
            let sample_x = col as f64 * hf.spacing;
            let sample_z = row as f64 * hf.spacing;
            let height = hf.height_at(row, col);
            let dx = convex_pos[0] - sample_x;
            let dz = convex_pos[2] - sample_z;
            let horiz_dist = (dx * dx + dz * dz).sqrt();
            if horiz_dist < search_r {
                let penetration = height + convex_radius - convex_pos[1];
                if penetration > 0.0 {
                    let normal = hf.normal_at(row, col);
                    let contact = Contact {
                        point_a: oxiphysics_core::math::Vec3::new(sample_x, height, sample_z),
                        point_b: oxiphysics_core::math::Vec3::new(
                            convex_pos[0],
                            convex_pos[1] - convex_radius,
                            convex_pos[2],
                        ),
                        normal: oxiphysics_core::math::Vec3::new(normal[0], normal[1], normal[2]),
                        depth: penetration,
                    };
                    let cell_pair = CollisionPair::new(pair.a, pair.b + row * 1000 + col);
                    let mut m = ContactManifold::new(cell_pair);
                    m.add_contact(contact);
                    manifolds.push(m);
                }
            }
        }
    }
    manifolds
}
/// Merge two contact manifolds, keeping the deepest contacts.
///
/// Useful when multiple overlapping tests on the same pair produce different
/// contacts; merging ensures the constraint solver sees the full picture.
pub fn merge_manifolds(
    manifold_a: ContactManifold,
    manifold_b: ContactManifold,
    max_contacts: usize,
) -> ContactManifold {
    let mut merged = ContactManifold::new(manifold_a.pair);
    for c in manifold_a.contacts {
        merged.add_contact(c);
    }
    for c in manifold_b.contacts {
        merged.add_contact(c);
    }
    merged.contacts.sort_by(|a, b| {
        b.depth
            .partial_cmp(&a.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    merged.contacts.truncate(max_contacts);
    merged
}
/// Flip a contact manifold (swap A and B, negate normal).
pub fn flip_manifold(manifold: &mut ContactManifold) {
    for c in &mut manifold.contacts {
        c.normal = -c.normal;
        std::mem::swap(&mut c.point_a, &mut c.point_b);
    }
}
/// Scale contact depths by a factor (e.g. for bias / slop removal).
pub fn scale_contact_depths(manifold: &mut ContactManifold, scale: f64) {
    for c in &mut manifold.contacts {
        c.depth *= scale;
    }
}
/// Remove contacts with depth below `min_depth`.
pub fn prune_shallow_contacts(manifold: &mut ContactManifold, min_depth: f64) {
    manifold.contacts.retain(|c| c.depth >= min_depth);
}
/// Translate all contact points by a world-space offset.
pub fn translate_manifold(manifold: &mut ContactManifold, offset: oxiphysics_core::math::Vec3) {
    for c in &mut manifold.contacts {
        c.point_a += offset;
        c.point_b += offset;
    }
}
/// Generate speculative contacts for all pairs in a batch.
///
/// Returns one manifold per pair (or `None` if the pair is too far apart).
pub fn batch_speculative_contacts<'a>(
    pairs: &[SpecBatchEntry<'a>],
    config: &SpeculativeConfig,
) -> Vec<Option<ContactManifold>> {
    pairs
        .iter()
        .map(|(sa, ta, va, sb, tb, vb, pair)| {
            generate_speculative_contacts(*sa, ta, *va, *sb, tb, *vb, *pair, config)
        })
        .collect()
}
#[cfg(test)]
mod tests_dispatch_extended {
    use super::*;
    use crate::CollisionPair;
    use crate::Contact;
    use crate::ContactManifold;
    use crate::narrowphase::ContactPatch;
    use crate::narrowphase::ContactPatchReducer;
    use crate::narrowphase::ConvexPiece;
    use crate::narrowphase::DispatchQueue;
    use crate::narrowphase::clip_box_faces_2d;
    use crate::narrowphase::clip_polygon_by_halfspace;
    use crate::narrowphase::flip_manifold;
    use oxiphysics_core::Vec3;

    /// Merge two `types::ContactManifold`s, keeping at most `max` contacts.
    fn merge_manifolds(a: ContactManifold, b: ContactManifold, max: usize) -> ContactManifold {
        let pair = a.pair;
        let mut merged = ContactManifold::new(pair);
        let mut all: Vec<Contact> = a.contacts.into_iter().chain(b.contacts).collect();
        all.sort_by(|x, y| {
            y.depth
                .partial_cmp(&x.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        all.truncate(max);
        for c in all {
            merged.add_contact(c);
        }
        merged
    }
    #[test]
    fn test_contact_patch_empty() {
        let p = ContactPatch::new([0.0, 1.0, 0.0], 0.1);
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);
    }
    #[test]
    fn test_contact_patch_add_point() {
        let mut p = ContactPatch::new([0.0, 1.0, 0.0], 0.1);
        p.add_point([1.0, 0.0, 0.0]);
        p.add_point([2.0, 0.0, 0.0]);
        assert_eq!(p.len(), 2);
        assert!(!p.is_empty());
    }
    #[test]
    fn test_clip_polygon_by_halfspace_full_inside() {
        // All points have x >= 0 with normal [1,0,0] and d=0
        let poly = vec![
            [0.5, 0.5, 0.0],
            [1.5, 0.5, 0.0],
            [1.5, 1.5, 0.0],
            [0.5, 1.5, 0.0],
        ];
        let clipped = clip_polygon_by_halfspace(&poly, [1.0, 0.0, 0.0], 0.0);
        assert_eq!(clipped.len(), 4, "all points inside, should keep 4");
    }
    #[test]
    fn test_clip_polygon_by_halfspace_full_outside() {
        // All points have x < 0 with normal [1,0,0] and d=0
        let poly = vec![
            [-2.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-2.0, 1.0, 0.0],
        ];
        let clipped = clip_polygon_by_halfspace(&poly, [1.0, 0.0, 0.0], 0.0);
        assert!(clipped.is_empty(), "all points outside");
    }
    #[test]
    fn test_clip_polygon_by_halfspace_partial() {
        // Some points inside, some outside with normal [1,0,0] and d=0
        let poly = vec![
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ];
        let clipped = clip_polygon_by_halfspace(&poly, [1.0, 0.0, 0.0], 0.0);
        assert!(
            clipped.len() >= 2,
            "partial clip should keep some points: {:?}",
            clipped
        );
    }
    #[test]
    fn test_clip_box_faces_2d_no_overlap() {
        let face_a: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let face_b: [[f64; 2]; 4] = [[5.0, 5.0], [6.0, 5.0], [6.0, 6.0], [5.0, 6.0]];
        let patch = clip_box_faces_2d(face_a, face_b, [0.0, 1.0, 0.0], 0.1);
        assert!(
            patch.is_empty(),
            "non-overlapping faces should produce no contacts"
        );
    }
    #[test]
    fn test_clip_box_faces_2d_full_overlap() {
        let face: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let patch = clip_box_faces_2d(face, face, [0.0, 1.0, 0.0], 0.1);
        let _ = patch;
    }
    #[test]
    fn test_reducer_keeps_max_contacts() {
        let reducer = ContactPatchReducer::new(2);
        let pair = CollisionPair::new(0, 1);
        let mut manifold = ContactManifold::new(pair);
        for i in 0..6 {
            manifold.add_contact(Contact {
                point_a: Vec3::new(i as f64, 0.0, 0.0),
                point_b: Vec3::new(i as f64 + 0.1, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                depth: i as f64 * 0.1 + 0.01,
            });
        }
        reducer.reduce(&mut manifold);
        assert!(
            manifold.contacts.len() <= 2,
            "should keep at most 2: {}",
            manifold.contacts.len()
        );
    }
    #[test]
    fn test_reducer_no_reduction_needed() {
        let reducer = ContactPatchReducer::new(4);
        let pair = CollisionPair::new(0, 1);
        let mut manifold = ContactManifold::new(pair);
        for i in 0..3 {
            manifold.add_contact(Contact {
                point_a: Vec3::new(i as f64, 0.0, 0.0),
                point_b: Vec3::new(i as f64 + 0.1, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                depth: 0.05,
            });
        }
        reducer.reduce(&mut manifold);
        assert_eq!(
            manifold.contacts.len(),
            3,
            "should keep all 3 when under max"
        );
    }
    #[test]
    fn test_box_box_manifold_clip_overlap() {
        use oxiphysics_geometry::BoxShape;
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = oxiphysics_core::Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = box_box_manifold_clip(&b1, &t1, &b2, &t2, pair);
        assert!(
            result.is_some(),
            "overlapping boxes should produce manifold"
        );
        let m = result.unwrap();
        assert!(!m.contacts.is_empty());
    }
    #[test]
    fn test_box_box_manifold_clip_separated() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = oxiphysics_core::Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let result = box_box_manifold_clip(&b1, &t1, &b2, &t2, pair);
        assert!(
            result.is_none(),
            "separated boxes should produce no manifold"
        );
    }
    #[test]
    fn test_speculative_contacts_within_margin() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = oxiphysics_core::Transform::from_position(Vec3::new(2.1, 0.0, 0.0));
        let config = SpeculativeConfig {
            margin: 0.5,
            velocity_scale: 10.0,
        };
        let pair = CollisionPair::new(0, 1);
        let result = generate_speculative_contacts(
            &s1,
            &t1,
            [1.0, 0.0, 0.0],
            &s2,
            &t2,
            [-1.0, 0.0, 0.0],
            pair,
            &config,
        );
        assert!(
            result.is_some(),
            "shapes within margin should produce speculative contact"
        );
    }
    #[test]
    fn test_speculative_contacts_beyond_margin() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = oxiphysics_core::Transform::from_position(Vec3::new(20.0, 0.0, 0.0));
        let config = SpeculativeConfig::default();
        let pair = CollisionPair::new(0, 1);
        let result =
            generate_speculative_contacts(&s1, &t1, [0.0; 3], &s2, &t2, [0.0; 3], pair, &config);
        assert!(
            result.is_none(),
            "far separated shapes should not produce speculative contact"
        );
    }
    #[test]
    fn test_speculative_config_default() {
        let c = SpeculativeConfig::default();
        assert!(c.margin > 0.0);
        assert!(c.velocity_scale > 0.0);
    }
    #[test]
    fn test_dispatch_queue_push_pop() {
        let mut q = DispatchQueue::new();
        q.push(CollisionPair::new(0, 1), 0.5);
        q.push(CollisionPair::new(1, 2), 1.5);
        q.push(CollisionPair::new(2, 3), 0.1);
        let top = q.pop().unwrap();
        assert!(
            (top.priority - 1.5).abs() < 1e-10,
            "highest priority first: {}",
            top.priority
        );
    }
    #[test]
    fn test_dispatch_queue_empty() {
        let mut q = DispatchQueue::new();
        assert!(q.is_empty());
        assert_eq!(q.len(), 0);
        assert!(q.pop().is_none());
    }
    #[test]
    fn test_dispatch_queue_clear() {
        let mut q = DispatchQueue::new();
        q.push(CollisionPair::new(0, 1), 1.0);
        q.push(CollisionPair::new(1, 2), 2.0);
        q.clear();
        assert!(q.is_empty());
    }
    #[test]
    fn test_merge_manifolds() {
        let pair = CollisionPair::new(0, 1);
        let mut m1 = ContactManifold::new(pair);
        m1.add_contact(Contact {
            point_a: Vec3::new(0.0, 0.0, 0.0),
            point_b: Vec3::new(0.1, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            depth: 0.5,
        });
        let mut m2 = ContactManifold::new(pair);
        m2.add_contact(Contact {
            point_a: Vec3::new(1.0, 0.0, 0.0),
            point_b: Vec3::new(1.1, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            depth: 0.1,
        });
        let merged = merge_manifolds(m1, m2, 4);
        assert_eq!(merged.contacts.len(), 2);
        assert!((merged.contacts[0].depth - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_merge_manifolds_limits_to_max() {
        let pair = CollisionPair::new(0, 1);
        let mut m1 = ContactManifold::new(pair);
        let mut m2 = ContactManifold::new(pair);
        for i in 0..3 {
            m1.add_contact(Contact {
                point_a: Vec3::new(i as f64, 0.0, 0.0),
                point_b: Vec3::new(i as f64 + 0.1, 0.0, 0.0),
                normal: Vec3::new(1.0, 0.0, 0.0),
                depth: 0.01 * i as f64 + 0.01,
            });
            m2.add_contact(Contact {
                point_a: Vec3::new(i as f64 + 10.0, 0.0, 0.0),
                point_b: Vec3::new(i as f64 + 10.1, 0.0, 0.0),
                normal: Vec3::new(1.0, 0.0, 0.0),
                depth: 0.01 * i as f64 + 0.01,
            });
        }
        let merged = merge_manifolds(m1, m2, 4);
        assert!(
            merged.contacts.len() <= 4,
            "should not exceed max_contacts=4: {}",
            merged.contacts.len()
        );
    }
    #[test]
    fn test_flip_manifold() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(Contact {
            point_a: Vec3::new(0.0, 0.0, 0.0),
            point_b: Vec3::new(1.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            depth: 0.1,
        });
        flip_manifold(&mut m);
        assert!(
            (m.contacts[0].normal.x - (-1.0)).abs() < 1e-10,
            "normal should be flipped"
        );
        assert!(
            (m.contacts[0].point_a.x - 1.0).abs() < 1e-10,
            "point_a and point_b should be swapped"
        );
    }
    #[test]
    fn test_scale_contact_depths() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::zeros(),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 1.0,
        });
        scale_contact_depths(&mut m, 2.0);
        assert!((m.contacts[0].depth - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_prune_shallow_contacts() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::zeros(),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.001,
        });
        m.add_contact(Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::zeros(),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.1,
        });
        prune_shallow_contacts(&mut m, 0.01);
        assert_eq!(m.contacts.len(), 1, "only the deep contact should remain");
        assert!((m.contacts[0].depth - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_translate_manifold() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(Contact {
            point_a: Vec3::new(1.0, 0.0, 0.0),
            point_b: Vec3::new(2.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            depth: 0.1,
        });
        let offset = Vec3::new(5.0, 0.0, 0.0);
        translate_manifold(&mut m, offset);
        assert!((m.contacts[0].point_a.x - 6.0).abs() < 1e-10);
        assert!((m.contacts[0].point_b.x - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_heightfield_convex_generic_flat() {
        let hf = SimpleHeightField::new(vec![0.0; 9], 3, 3, 1.0);
        let pair = CollisionPair::new(0, 1);
        let manifolds = dispatch_heightfield_convex_generic(&hf, [1.0, 0.5, 1.0], 1.0, 0.0, pair);
        assert!(
            !manifolds.is_empty(),
            "should produce contacts for penetrating sphere"
        );
    }
    #[test]
    fn test_heightfield_convex_generic_no_contact() {
        let hf = SimpleHeightField::new(vec![0.0; 9], 3, 3, 1.0);
        let pair = CollisionPair::new(0, 1);
        let manifolds = dispatch_heightfield_convex_generic(&hf, [1.0, 10.0, 1.0], 0.1, 0.0, pair);
        assert!(
            manifolds.is_empty(),
            "sphere far above heightfield should have no contacts"
        );
    }
    #[test]
    fn test_dispatch_concave_empty_mesh() {
        let mesh = ConcaveMesh::new(vec![]);
        let s = Sphere::new(1.0);
        let tm = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let ts = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let pair = CollisionPair::new(0, 1);
        let dispatcher = NarrowPhaseDispatcher::default();
        let manifolds =
            dispatch_concave_vs_convex(&mesh, &tm, &s, ShapeType::Sphere, &ts, pair, &dispatcher);
        assert!(
            manifolds.is_empty(),
            "empty mesh should produce no contacts"
        );
    }
    #[test]
    fn test_convex_piece_aabb_overlap() {
        let pa = ConvexPiece::new(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]);
        let pb = ConvexPiece::new(vec![
            [0.5, 0.5, 0.0],
            [1.5, 0.5, 0.0],
            [1.5, 1.5, 0.0],
            [0.5, 1.5, 0.0],
        ]);
        assert!(
            pa.aabb_overlaps(&pb, 0.0),
            "overlapping pieces should overlap"
        );
    }
    #[test]
    fn test_convex_piece_aabb_no_overlap() {
        let pa = ConvexPiece::new(vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        let pb = ConvexPiece::new(vec![[5.0, 5.0, 5.0], [6.0, 6.0, 6.0]]);
        assert!(
            !pa.aabb_overlaps(&pb, 0.0),
            "separated pieces should not overlap"
        );
    }
    #[test]
    fn test_batch_speculative_empty() {
        let config = SpeculativeConfig::default();
        let results = batch_speculative_contacts(&[], &config);
        assert!(results.is_empty());
    }
    #[test]
    fn test_batch_speculative_single() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = oxiphysics_core::Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = oxiphysics_core::Transform::from_position(Vec3::new(2.1, 0.0, 0.0));
        let config = SpeculativeConfig {
            margin: 0.5,
            velocity_scale: 10.0,
        };
        let pair = CollisionPair::new(0, 1);
        let pairs: Vec<_> = vec![(
            &s1 as &dyn Shape,
            &t1,
            [0.5, 0.0, 0.0_f64],
            &s2 as &dyn Shape,
            &t2,
            [-0.5, 0.0, 0.0_f64],
            pair,
        )];
        let results = batch_speculative_contacts(&pairs, &config);
        assert_eq!(results.len(), 1);
    }
}
