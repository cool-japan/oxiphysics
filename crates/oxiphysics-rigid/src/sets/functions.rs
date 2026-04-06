//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::RigidBodySet;

/// Alias for `RigidBodySet` used in contexts that don't need the "Rigid" prefix.
pub type BodySet = RigidBodySet;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Aabb3;
    use crate::BodyInterpolator;
    use crate::BodyPairCache;
    use crate::BodyState;
    use crate::BvhTree;
    use crate::Collider;
    use crate::ColliderSet;
    use crate::ConstraintGraph;
    use crate::RigidBody;
    use crate::SpatialHash;
    use oxiphysics_core::BodyHandle;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_geometry::{Shape, Sphere};
    use std::sync::Arc;
    #[test]
    fn test_body_set_insert_get() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        assert_eq!(set.len(), 1);
        assert!(set.get(h).is_some());
    }
    #[test]
    fn test_body_set_remove() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.remove(h);
        assert_eq!(set.len(), 0);
        assert!(set.get(h).is_none());
    }
    #[test]
    fn test_body_set_generation() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.remove(h1);
        let h2 = set.insert(RigidBody::new(2.0));
        assert!(set.get(h1).is_none());
        assert!(set.get(h2).is_some());
        assert_eq!(h2.index, h1.index);
        assert_ne!(h2.generation, h1.generation);
    }
    #[test]
    fn test_collider_set() {
        let mut set = ColliderSet::new();
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let h = set.insert(Collider::new(shape));
        assert_eq!(set.len(), 1);
        assert!(set.get(h).is_some());
    }
    #[test]
    fn test_body_set_iter() {
        let mut set = RigidBodySet::new();
        let _h1 = set.insert(RigidBody::new(1.0));
        let _h2 = set.insert(RigidBody::new(2.0));
        let bodies: Vec<_> = set.iter().collect();
        assert_eq!(bodies.len(), 2);
    }
    #[test]
    fn test_body_set_iter_mut() {
        let mut set = RigidBodySet::new();
        let _h1 = set.insert(RigidBody::new(1.0));
        for (_, body) in set.iter_mut() {
            body.velocity = Vec3::new(1.0, 0.0, 0.0);
        }
        let body = set.iter().next().unwrap().1;
        assert!((body.velocity.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_attach_collider_updates_body() {
        let mut body_set = RigidBodySet::new();
        let body_handle = body_set.insert(RigidBody::new(1.0));
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        let collider = Collider::new(shape).with_density(1000.0);
        body_set.attach_collider_with_auto_inertia(body_handle, &collider);
        let body = body_set.get(body_handle).unwrap();
        let expected_mass = collider.mass_properties().mass;
        assert!(
            (body.mass - expected_mass).abs() < 1e-6,
            "body mass={} expected={}",
            body.mass,
            expected_mass
        );
        assert!(body.mass > 100.0);
    }
    #[test]
    fn apply_gravity_increases_force_accumulator() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(2.0));
        set.apply_gravity([0.0, -9.81, 0.0], 0.01);
        let fy = set.get(h).unwrap().force_accumulator.y;
        assert!((fy - (-19.62)).abs() < 1e-6, "fy={fy}");
    }
    #[test]
    fn apply_gravity_skips_static_bodies() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new_static());
        set.apply_gravity([0.0, -9.81, 0.0], 0.01);
        let fy = set.get(h).unwrap().force_accumulator.y;
        assert!(
            fy.abs() < 1e-12,
            "static body should not accumulate gravity force"
        );
    }
    #[test]
    fn step_velocities_integrates_force() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        {
            let b = set.get_mut(h).unwrap();
            b.force_accumulator = Vec3::new(10.0, 0.0, 0.0);
        }
        set.step_velocities(0.1);
        let vx = set.get(h).unwrap().velocity.x;
        assert!((vx - 1.0).abs() < 1e-10, "vx={vx}");
    }
    #[test]
    fn clear_forces_zeroes_accumulators() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        {
            let b = set.get_mut(h).unwrap();
            b.force_accumulator = Vec3::new(5.0, 5.0, 5.0);
        }
        set.clear_forces();
        let f = set.get(h).unwrap().force_accumulator;
        assert!(f.norm() < 1e-12);
    }
    #[test]
    fn active_sleeping_iterators_partition() {
        use crate::body::BodyState;
        let mut set = RigidBodySet::new();
        let _h1 = set.insert(RigidBody::new(1.0));
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h2).unwrap().state = BodyState::Sleeping;
        let active: Vec<_> = set.active_bodies().collect();
        let sleeping: Vec<_> = set.sleeping_bodies().collect();
        assert_eq!(active.len(), 1);
        assert_eq!(sleeping.len(), 1);
    }
    #[test]
    fn constraint_graph_single_island() {
        let mut g = ConstraintGraph::new();
        let h0 = BodyHandle::new(0, 0);
        let h1 = BodyHandle::new(1, 0);
        let h2 = BodyHandle::new(2, 0);
        g.add_constraint(h0, h1);
        g.add_constraint(h1, h2);
        let islands = g.find_islands(&[h0, h1, h2]);
        assert_eq!(islands.len(), 1);
        assert_eq!(islands[0].len(), 3);
    }
    #[test]
    fn constraint_graph_two_islands() {
        let mut g = ConstraintGraph::new();
        let h0 = BodyHandle::new(0, 0);
        let h1 = BodyHandle::new(1, 0);
        let h2 = BodyHandle::new(2, 0);
        let h3 = BodyHandle::new(3, 0);
        g.add_constraint(h0, h1);
        g.add_constraint(h2, h3);
        let islands = g.find_islands(&[h0, h1, h2, h3]);
        assert_eq!(islands.len(), 2);
    }
    #[test]
    fn constraint_graph_singletons_without_constraints() {
        let g = ConstraintGraph::new();
        let handles: Vec<BodyHandle> = (0..4).map(|i| BodyHandle::new(i, 0)).collect();
        let islands = g.find_islands(&handles);
        assert_eq!(
            islands.len(),
            4,
            "each unconstrained body is its own island"
        );
    }
    #[test]
    fn test_dynamic_static_kinematic_iterators() {
        use crate::body::BodyType;
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new_static());
        {
            let mut b = RigidBody::new(1.0);
            b.body_type = BodyType::Kinematic;
            set.insert(b);
        }
        assert_eq!(set.dynamic_bodies().count(), 1);
        assert_eq!(set.static_bodies().count(), 1);
        assert_eq!(set.kinematic_bodies().count(), 1);
    }
    #[test]
    fn test_capacity_and_free_count() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        let _h2 = set.insert(RigidBody::new(2.0));
        set.remove(h1);
        assert_eq!(set.capacity(), 2);
        assert_eq!(set.free_count(), 1);
    }
    #[test]
    fn test_bodies_by_mass_range() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(5.0));
        set.insert(RigidBody::new(10.0));
        let found = set.bodies_by_mass_range(4.0, 6.0);
        assert_eq!(found.len(), 1);
        assert!((found[0].1.mass - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_heaviest_body() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(100.0));
        set.insert(RigidBody::new(5.0));
        let (_, body) = set.heaviest_body().unwrap();
        assert!((body.mass - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_lightest_dynamic_body() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(10.0));
        set.insert(RigidBody::new(2.0));
        set.insert(RigidBody::new(50.0));
        let (_, body) = set.lightest_dynamic_body().unwrap();
        assert!((body.mass - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_fastest_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(10.0, 0.0, 0.0);
        set.insert(RigidBody::new(1.0));
        let (fh, _) = set.fastest_body().unwrap();
        assert_eq!(fh.index, h.index);
    }
    #[test]
    fn test_bodies_within_radius() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h2).unwrap().transform.position = Vec3::new(100.0, 0.0, 0.0);
        let found = set.bodies_within_radius(Vec3::zeros(), 5.0);
        assert_eq!(found.len(), 1);
    }
    #[test]
    fn test_bodies_in_aabb() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(0.5, 0.5, 0.5);
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h2).unwrap().transform.position = Vec3::new(10.0, 10.0, 10.0);
        let found = set.bodies_in_aabb(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(found.len(), 1);
    }
    #[test]
    fn test_nearest_body() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(1.0, 0.0, 0.0);
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h2).unwrap().transform.position = Vec3::new(10.0, 0.0, 0.0);
        let (nh, _, dist) = set.nearest_body(Vec3::zeros()).unwrap();
        assert_eq!(nh.index, h1.index);
        assert!((dist - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_freeze_all() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(5.0, 5.0, 5.0);
        set.freeze_all();
        assert!(set.get(h).unwrap().velocity.norm() < 1e-12);
    }
    #[test]
    fn test_remove_where() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(100.0));
        set.insert(RigidBody::new(2.0));
        let removed = set.remove_where(|b| b.mass > 50.0);
        assert_eq!(removed.len(), 1);
        assert_eq!(set.len(), 2);
    }
    #[test]
    fn test_for_each_mut() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(2.0));
        set.for_each_mut(|_| true, |b| b.velocity = Vec3::new(1.0, 0.0, 0.0));
        for (_, b) in set.iter() {
            assert!((b.velocity.x - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_statistics() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(3.0));
        set.insert(RigidBody::new_static());
        let h = set.insert(RigidBody::new(2.0));
        set.get_mut(h).unwrap().state = BodyState::Sleeping;
        let stats = set.statistics();
        assert_eq!(stats.total_count, 4);
        assert_eq!(stats.dynamic_count, 3);
        assert_eq!(stats.static_count, 1);
        assert_eq!(stats.active_count, 3);
        assert_eq!(stats.sleeping_count, 1);
    }
    #[test]
    fn test_export_import_positions() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(1.0, 2.0, 3.0);
        set.get_mut(h2).unwrap().transform.position = Vec3::new(4.0, 5.0, 6.0);
        let exported = set.export_positions();
        assert_eq!(exported.len(), 6);
        assert!((exported[0] - 1.0).abs() < 1e-12);
        assert!((exported[3] - 4.0).abs() < 1e-12);
        let new_pos = vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
        set.import_positions(&new_pos);
        assert!((set.get(h1).unwrap().transform.position.x - 10.0).abs() < 1e-12);
        assert!((set.get(h2).unwrap().transform.position.x - 40.0).abs() < 1e-12);
    }
    #[test]
    fn test_export_velocities() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(1.0, 2.0, 3.0);
        let vel = set.export_velocities();
        assert_eq!(vel.len(), 3);
        assert!((vel[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_export_masses() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(3.0));
        set.insert(RigidBody::new(7.0));
        let masses = set.export_masses();
        assert_eq!(masses.len(), 2);
    }
    #[test]
    fn test_all_handles() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        set.insert(RigidBody::new(2.0));
        set.insert(RigidBody::new(3.0));
        assert_eq!(set.all_handles().len(), 3);
    }
    #[test]
    fn test_active_handles() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new(1.0));
        let h2 = set.insert(RigidBody::new(2.0));
        set.get_mut(h2).unwrap().state = BodyState::Sleeping;
        assert_eq!(set.active_handles().len(), 1);
    }
    #[test]
    fn test_constraint_graph_degree() {
        let mut g = ConstraintGraph::new();
        let h0 = BodyHandle::new(0, 0);
        let h1 = BodyHandle::new(1, 0);
        let h2 = BodyHandle::new(2, 0);
        g.add_constraint(h0, h1);
        g.add_constraint(h0, h2);
        assert_eq!(g.degree(h0), 2);
        assert_eq!(g.degree(h1), 1);
    }
    #[test]
    fn test_constraint_graph_are_connected() {
        let mut g = ConstraintGraph::new();
        let h0 = BodyHandle::new(0, 0);
        let h1 = BodyHandle::new(1, 0);
        let h2 = BodyHandle::new(2, 0);
        g.add_constraint(h0, h1);
        assert!(g.are_connected(h0, h1));
        assert!(g.are_connected(h1, h0));
        assert!(!g.are_connected(h0, h2));
    }
    #[test]
    fn test_constraint_graph_node_count() {
        let mut g = ConstraintGraph::new();
        let h0 = BodyHandle::new(0, 0);
        let h1 = BodyHandle::new(1, 0);
        g.add_constraint(h0, h1);
        assert_eq!(g.node_count(), 2);
    }
    #[test]
    fn test_import_velocities() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.import_velocities(&[7.0, 8.0, 9.0]);
        let v = set.get(h).unwrap().velocity;
        assert!((v.x - 7.0).abs() < 1e-12);
        assert!((v.y - 8.0).abs() < 1e-12);
        assert!((v.z - 9.0).abs() < 1e-12);
    }
    #[test]
    fn spatial_hash_insert_and_query() {
        let mut hash: SpatialHash = SpatialHash::new(2.0);
        hash.insert(0, [0.5, 0.5, 0.5]);
        hash.insert(1, [1.5, 0.5, 0.5]);
        hash.insert(2, [100.0, 100.0, 100.0]);
        let near = hash.query_cell([0.5, 0.5, 0.5]);
        assert!(near.contains(&0));
        assert!(near.contains(&1));
        assert!(!near.contains(&2));
    }
    #[test]
    fn spatial_hash_clear() {
        let mut hash: SpatialHash = SpatialHash::new(1.0);
        hash.insert(0, [0.0, 0.0, 0.0]);
        hash.clear();
        let near = hash.query_cell([0.0, 0.0, 0.0]);
        assert!(near.is_empty());
    }
    #[test]
    fn spatial_hash_radius_query() {
        let mut hash: SpatialHash = SpatialHash::new(1.0);
        hash.insert(0, [0.0, 0.0, 0.0]);
        hash.insert(1, [0.9, 0.0, 0.0]);
        hash.insert(2, [5.0, 0.0, 0.0]);
        let near = hash.query_radius([0.0, 0.0, 0.0], 1.5);
        assert!(near.contains(&0));
        assert!(near.contains(&1));
        assert!(!near.contains(&2));
    }
    #[test]
    fn spatial_hash_rebuild_from_body_set() {
        let mut set = RigidBodySet::new();
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let h2 = set.insert(RigidBody::new(1.0));
        set.get_mut(h2).unwrap().transform.position = Vec3::new(50.0, 0.0, 0.0);
        let mut hash = SpatialHash::new(5.0);
        for (hndl, body) in set.iter() {
            let p = body.transform.position;
            hash.insert(hndl.index as usize, [p.x, p.y, p.z]);
        }
        let near = hash.query_radius([0.0, 0.0, 0.0], 4.0);
        assert_eq!(near.len(), 1);
        assert_eq!(near[0], h1.index as usize);
    }
    #[test]
    fn aabb3_contains_point() {
        let aabb = Aabb3::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(aabb.contains([0.0, 0.0, 0.0]));
        assert!(!aabb.contains([2.0, 0.0, 0.0]));
    }
    #[test]
    fn aabb3_overlap() {
        let a = Aabb3::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Aabb3::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.overlaps(&b));
        let c = Aabb3::new([3.0, 3.0, 3.0], [5.0, 5.0, 5.0]);
        assert!(!a.overlaps(&c));
    }
    #[test]
    fn aabb3_surface_area() {
        let a = Aabb3::new([0.0; 3], [1.0, 2.0, 3.0]);
        assert!((a.surface_area() - 22.0).abs() < 1e-10);
    }
    #[test]
    fn aabb3_merge() {
        let a = Aabb3::new([0.0; 3], [1.0; 3]);
        let b = Aabb3::new([-1.0; 3], [2.0; 3]);
        let merged = a.merge(&b);
        assert_eq!(merged.min, [-1.0, -1.0, -1.0]);
        assert_eq!(merged.max, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn bvh_build_and_query() {
        let aabbs: Vec<(usize, Aabb3)> = vec![
            (0, Aabb3::new([0.0; 3], [1.0; 3])),
            (1, Aabb3::new([2.0; 3], [3.0; 3])),
            (2, Aabb3::new([10.0; 3], [11.0; 3])),
        ];
        let bvh = BvhTree::build(&aabbs);
        let query = Aabb3::new([-0.5; 3], [1.5; 3]);
        let hits = bvh.query(&query);
        assert!(hits.contains(&0), "leaf 0 should be hit");
        assert!(!hits.contains(&2), "leaf 2 should not be hit");
    }
    #[test]
    fn bvh_empty() {
        let bvh = BvhTree::build(&[]);
        let query = Aabb3::new([-1.0; 3], [1.0; 3]);
        assert!(bvh.query(&query).is_empty());
    }
    #[test]
    fn interpolator_stores_and_retrieves_snapshot() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(2.0, 0.0, 0.0);
        set.get_mut(h).unwrap().transform.position = Vec3::new(1.0, 0.0, 0.0);
        let mut interp = BodyInterpolator::default();
        interp.snapshot_previous(&set);
        set.get_mut(h).unwrap().transform.position = Vec3::new(3.0, 0.0, 0.0);
        let pos = interp.interpolate_position(h, set.get(h).unwrap(), 0.5);
        assert!((pos.x - 2.0).abs() < 1e-10, "pos.x={}", pos.x);
    }
    #[test]
    fn interpolator_alpha_zero_returns_previous() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let mut interp = BodyInterpolator::default();
        interp.snapshot_previous(&set);
        set.get_mut(h).unwrap().transform.position = Vec3::new(10.0, 0.0, 0.0);
        let pos = interp.interpolate_position(h, set.get(h).unwrap(), 0.0);
        assert!((pos.x - 0.0).abs() < 1e-10, "pos.x={}", pos.x);
    }
    #[test]
    fn interpolator_alpha_one_returns_current() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let mut interp = BodyInterpolator::default();
        interp.snapshot_previous(&set);
        set.get_mut(h).unwrap().transform.position = Vec3::new(10.0, 0.0, 0.0);
        let pos = interp.interpolate_position(h, set.get(h).unwrap(), 1.0);
        assert!((pos.x - 10.0).abs() < 1e-10, "pos.x={}", pos.x);
    }
    #[test]
    fn interpolator_velocity_blend() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(0.0, 0.0, 0.0);
        let mut interp = BodyInterpolator::default();
        interp.snapshot_previous(&set);
        set.get_mut(h).unwrap().velocity = Vec3::new(10.0, 0.0, 0.0);
        let v = interp.interpolate_velocity(h, set.get(h).unwrap(), 0.3);
        assert!((v.x - 3.0).abs() < 1e-10, "v.x={}", v.x);
    }
    #[test]
    fn interpolator_clear_and_fallback() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(5.0, 0.0, 0.0);
        let mut interp = BodyInterpolator::default();
        interp.snapshot_previous(&set);
        interp.clear();
        assert!(interp.is_empty());
        let pos = interp.interpolate_position(h, set.get(h).unwrap(), 0.5);
        assert!((pos.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn body_pair_cache_set_get() {
        let mut cache = BodyPairCache::new();
        let h0 = oxiphysics_core::BodyHandle::new(0, 0);
        let h1 = oxiphysics_core::BodyHandle::new(1, 0);
        cache.set(h0, h1, 42.0);
        assert_eq!(cache.get(h0, h1), Some(42.0));
        assert_eq!(cache.get(h1, h0), Some(42.0));
    }
    #[test]
    fn body_pair_cache_remove() {
        let mut cache = BodyPairCache::new();
        let h0 = oxiphysics_core::BodyHandle::new(0, 0);
        let h1 = oxiphysics_core::BodyHandle::new(1, 0);
        cache.set(h0, h1, 1.0);
        cache.remove(h0, h1);
        assert_eq!(cache.get(h0, h1), None);
    }
    #[test]
    fn body_pair_cache_evict_stale() {
        let mut set = RigidBodySet::new();
        let h0 = set.insert(RigidBody::new(1.0));
        let h1 = set.insert(RigidBody::new(1.0));
        let h2_phantom = oxiphysics_core::BodyHandle::new(99, 0);
        let mut cache = BodyPairCache::new();
        cache.set(h0, h1, 1.0);
        cache.set(h0, h2_phantom, 2.0);
        cache.evict_stale(&set);
        assert!(cache.get(h0, h1).is_some());
        assert_eq!(cache.get(h0, h2_phantom), None);
    }
    #[test]
    fn body_pair_cache_len() {
        let mut cache = BodyPairCache::new();
        assert!(cache.is_empty());
        let h0 = oxiphysics_core::BodyHandle::new(0, 0);
        let h1 = oxiphysics_core::BodyHandle::new(1, 0);
        cache.set(h0, h1, 0.0);
        assert_eq!(cache.len(), 1);
        cache.clear();
        assert!(cache.is_empty());
    }
    #[test]
    fn aabb3_from_points() {
        let pts = vec![[1.0, 2.0, 3.0], [-1.0, 0.0, 4.0], [2.0, -1.0, 1.0]];
        let aabb = Aabb3::from_points(&pts).unwrap();
        assert_eq!(aabb.min, [-1.0, -1.0, 1.0]);
        assert_eq!(aabb.max, [2.0, 2.0, 4.0]);
    }
    #[test]
    fn aabb3_expand() {
        let a = Aabb3::new([0.0; 3], [1.0; 3]);
        let b = a.expand(0.5);
        assert_eq!(b.min, [-0.5; 3]);
        assert_eq!(b.max, [1.5; 3]);
    }
    #[test]
    fn aabb3_sq_dist_inside() {
        let a = Aabb3::new([-1.0; 3], [1.0; 3]);
        assert_eq!(a.sq_dist_to_point([0.0, 0.0, 0.0]), 0.0);
    }
    #[test]
    fn aabb3_sq_dist_outside() {
        let a = Aabb3::new([0.0; 3], [1.0; 3]);
        let d2 = a.sq_dist_to_point([2.0, 0.0, 0.0]);
        assert!((d2 - 1.0).abs() < 1e-10, "d2={d2}");
    }
    #[test]
    fn aabb3_volume() {
        let a = Aabb3::new([0.0; 3], [2.0, 3.0, 4.0]);
        assert!((a.volume() - 24.0).abs() < 1e-10);
    }
    #[test]
    fn bvh_rebuild_from_body_set() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let bvh = BvhTree::rebuild_from_body_set(&set, 0.5);
        assert_eq!(bvh.node_count(), 1);
        let q = Aabb3::new([-0.1; 3], [0.1; 3]);
        let hits = bvh.query(&q);
        assert_eq!(hits.len(), 1);
    }
    #[test]
    fn bvh_multiple_primitives_all_hit_by_large_query() {
        let aabbs: Vec<(usize, Aabb3)> = (0..5)
            .map(|i| {
                let x = i as f64 * 2.0;
                (i, Aabb3::new([x, 0.0, 0.0], [x + 1.0, 1.0, 1.0]))
            })
            .collect();
        let bvh = BvhTree::build(&aabbs);
        let big_query = Aabb3::new([-1.0; 3], [20.0; 3]);
        let hits = bvh.query(&big_query);
        assert_eq!(hits.len(), 5, "all 5 primitives should be hit");
    }
    #[test]
    fn bvh_no_hits_outside_query() {
        let aabbs = vec![(0usize, Aabb3::new([0.0; 3], [1.0; 3]))];
        let bvh = BvhTree::build(&aabbs);
        let q = Aabb3::new([5.0; 3], [6.0; 3]);
        assert!(bvh.query(&q).is_empty());
    }
    #[test]
    fn test_total_momentum_single_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(2.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(3.0, 0.0, 0.0);
        let (p, _) = set.compute_total_momentum();
        assert!((p[0] - 6.0).abs() < 1e-10, "px={}", p[0]);
        assert!(p[1].abs() < 1e-10);
        assert!(p[2].abs() < 1e-10);
    }
    #[test]
    fn test_total_momentum_static_body_excluded() {
        let mut set = RigidBodySet::new();
        set.insert(RigidBody::new_static());
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().velocity = Vec3::new(5.0, 0.0, 0.0);
        let (p, _) = set.compute_total_momentum();
        assert!((p[0] - 5.0).abs() < 1e-10, "px={}", p[0]);
    }
    #[test]
    fn test_total_momentum_two_equal_bodies_opposite_velocity() {
        let mut set = RigidBodySet::new();
        let h0 = set.insert(RigidBody::new(1.0));
        set.get_mut(h0).unwrap().velocity = Vec3::new(1.0, 0.0, 0.0);
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().velocity = Vec3::new(-1.0, 0.0, 0.0);
        let (p, _) = set.compute_total_momentum();
        assert!(p[0].abs() < 1e-10, "px should cancel: {}", p[0]);
    }
    #[test]
    fn test_com_single_body() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(4.0, 5.0, 6.0);
        let com = set.compute_center_of_mass().unwrap();
        assert!((com[0] - 4.0).abs() < 1e-10);
        assert!((com[1] - 5.0).abs() < 1e-10);
        assert!((com[2] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_com_two_equal_mass_bodies() {
        let mut set = RigidBodySet::new();
        let h0 = set.insert(RigidBody::new(1.0));
        set.get_mut(h0).unwrap().transform.position = Vec3::new(0.0, 0.0, 0.0);
        let h1 = set.insert(RigidBody::new(1.0));
        set.get_mut(h1).unwrap().transform.position = Vec3::new(2.0, 0.0, 0.0);
        let com = set.compute_center_of_mass().unwrap();
        assert!((com[0] - 1.0).abs() < 1e-10, "com_x={}", com[0]);
    }
    #[test]
    fn test_com_empty_set_returns_none() {
        let set = RigidBodySet::new();
        assert!(set.compute_center_of_mass().is_none());
    }
    #[test]
    fn test_explosion_impulse_body_in_range_gains_velocity() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(1.0, 0.0, 0.0);
        set.apply_explosion_impulse([0.0; 3], 5.0, 100.0);
        let v = set.get(h).unwrap().velocity;
        assert!(
            v.x > 0.0,
            "body should gain velocity from explosion: vx={}",
            v.x
        );
    }
    #[test]
    fn test_explosion_impulse_body_outside_range_unaffected() {
        let mut set = RigidBodySet::new();
        let h = set.insert(RigidBody::new(1.0));
        set.get_mut(h).unwrap().transform.position = Vec3::new(100.0, 0.0, 0.0);
        set.apply_explosion_impulse([0.0; 3], 5.0, 100.0);
        let v = set.get(h).unwrap().velocity;
        assert!(
            v.x.abs() < 1e-12,
            "body outside blast radius should be unaffected: vx={}",
            v.x
        );
    }
    #[test]
    fn test_batch_update_transforms_updates_position_and_velocity() {
        let mut set = RigidBodySet::new();
        let h0 = set.insert(RigidBody::new(1.0));
        let h1 = set.insert(RigidBody::new(1.0));
        let transforms = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        set.batch_update_transforms(&[h0, h1], &transforms);
        let b0 = set.get(h0).unwrap();
        assert!((b0.transform.position.x - 1.0).abs() < 1e-10);
        assert!((b0.velocity.x - 4.0).abs() < 1e-10);
        let b1 = set.get(h1).unwrap();
        assert!((b1.transform.position.x - 7.0).abs() < 1e-10);
        assert!((b1.velocity.y - 11.0).abs() < 1e-10);
    }
}
