#![allow(clippy::ptr_arg)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::math::Real;

use super::types::{BodyArena, BodyMode, Contact};

#[cfg(test)]
mod tests {

    use crate::Body;
    use crate::BodyArena;
    use crate::BodyHandle;

    use crate::BodyVelocity;
    use crate::PhysicsEvent;
    use crate::PhysicsWorld;
    use crate::TimeStep;
    #[test]
    fn test_body_velocity_linear_speed_sq() {
        let v = BodyVelocity::new([3.0, 4.0, 0.0], [0.0; 3]);
        assert!((v.linear_speed_sq() - 25.0).abs() < 1e-12);
    }
    #[test]
    fn test_body_velocity_damp() {
        let mut v = BodyVelocity::new([10.0, 0.0, 0.0], [0.0; 3]);
        v.damp(1.0, 1.0, 0.5);
        assert!((v.linear[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_body_velocity_apply_linear_impulse() {
        let mut v = BodyVelocity::default();
        v.apply_linear_impulse([10.0, 0.0, 0.0], 0.5);
        assert!((v.linear[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_body_apply_force_static_noop() {
        let mut body = Body::new_static([0.0; 3]);
        body.apply_force([100.0, 0.0, 0.0]);
        assert_eq!(
            body.force[0], 0.0,
            "static bodies should not accumulate force"
        );
    }
    #[test]
    fn test_body_integrate_velocities() {
        let mut body = Body::new_dynamic([0.0; 3], 2.0);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        body.integrate_velocities([0.0, -10.0, 0.0], 1.0);
        assert!((body.velocity.linear[1] - (-10.0)).abs() < 1e-10);
    }
    #[test]
    fn test_body_integrate_positions() {
        let mut body = Body::new_dynamic([0.0; 3], 1.0);
        body.velocity.linear = [1.0, 0.0, 0.0];
        body.integrate_positions(2.0);
        assert!((body.transform.position.x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_body_clear_forces() {
        let mut body = Body::new_dynamic([0.0; 3], 1.0);
        body.apply_force([5.0, 5.0, 5.0]);
        body.clear_forces();
        assert_eq!(body.force, [0.0; 3]);
    }
    #[test]
    fn test_body_arena_insert_remove() {
        let mut arena = BodyArena::new();
        let handle = arena.insert(Body::new_dynamic([0.0; 3], 1.0));
        assert_eq!(arena.len(), 1);
        let removed = arena.remove(handle);
        assert!(removed.is_some());
        assert_eq!(arena.len(), 0);
    }
    #[test]
    fn test_body_arena_stale_handle() {
        let mut arena = BodyArena::new();
        let h = arena.insert(Body::new_dynamic([0.0; 3], 1.0));
        arena.remove(h);
        assert!(arena.get(h).is_none());
    }
    #[test]
    fn test_body_arena_generation_reuse() {
        let mut arena = BodyArena::new();
        let h1 = arena.insert(Body::new_dynamic([0.0; 3], 1.0));
        arena.remove(h1);
        let h2 = arena.insert(Body::new_dynamic([1.0, 0.0, 0.0], 2.0));
        assert_eq!(h1.index, h2.index);
        assert_ne!(h1.generation, h2.generation);
        assert!(arena.get(h2).is_some());
        assert!(arena.get(h1).is_none());
    }
    #[test]
    fn test_world_add_remove_body() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        assert_eq!(world.body_count(), 1);
        world.remove_body(h);
        assert_eq!(world.body_count(), 0);
    }
    #[test]
    fn test_world_step_accumulates_time() {
        let mut world = PhysicsWorld::new();
        let ts = TimeStep::new(1.0 / 60.0);
        world.step(&ts);
        assert!((world.time - 1.0 / 60.0).abs() < 1e-12);
    }
    #[test]
    fn test_world_step_with_accumulator_substeps() {
        let mut world = PhysicsWorld::new();
        world.fixed_dt = 1.0 / 60.0;
        let steps = world.step_with_accumulator(1.0 / 30.0);
        assert_eq!(steps, 2);
    }
    #[test]
    fn test_world_gravity_applies() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.body_mut(h).unwrap().linear_damping = 0.0;
        world.step(&TimeStep::new(1.0));
        let vy = world.body(h).unwrap().velocity.linear[1];
        assert!((vy - (-9.81)).abs() < 1e-6, "vy={}", vy);
    }
    #[test]
    fn test_world_bodies_in_radius() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.add_body(Body::new_dynamic([100.0, 0.0, 0.0], 1.0));
        let nearby = world.bodies_in_radius([0.0; 3], 5.0);
        assert_eq!(nearby.len(), 1);
    }
    #[test]
    fn test_world_kinetic_energy_at_rest() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_static([0.0; 3]));
        assert_eq!(world.total_kinetic_energy(), 0.0);
    }
    #[test]
    fn test_world_compute_islands_single_body() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.compute_islands();
        assert_eq!(world.islands.len(), 1);
        assert_eq!(world.islands[0].size(), 1);
    }
    #[test]
    fn test_world_event_queue() {
        let mut world = PhysicsWorld::new();
        let h = BodyHandle::new(0, 0);
        world.push_event(PhysicsEvent::BodySlept { body: h });
        let events = world.drain_events();
        assert_eq!(events.len(), 1);
        assert!(world.events.is_empty());
    }
    #[test]
    fn test_world_total_momentum_at_rest() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_static([0.0; 3]));
        let p = world.total_momentum();
        assert_eq!(p, [0.0; 3]);
    }
}
/// Resolve a single contact impulse between body A and optionally body B.
///
/// Modifies the velocity fields of the bodies in the arena in-place.
/// This is a single-iteration impulse solver (Dantzig / sequential impulse).
pub fn resolve_contact(bodies: &mut BodyArena, contact: &Contact) {
    let pa = {
        if let Some(b) = bodies.get(contact.body_a) {
            [
                b.transform.position.x,
                b.transform.position.y,
                b.transform.position.z,
            ]
        } else {
            return;
        }
    };
    let (va, inv_ma, inv_ia) = {
        if let Some(b) = bodies.get(contact.body_a) {
            (
                b.velocity.linear,
                b.inv_mass,
                if b.inertia > 0.0 {
                    1.0 / b.inertia
                } else {
                    0.0
                },
            )
        } else {
            return;
        }
    };
    let (vb, inv_mb, inv_ib, pb) = if let Some(bh) = contact.body_b {
        if let Some(b) = bodies.get(bh) {
            let pb = [
                b.transform.position.x,
                b.transform.position.y,
                b.transform.position.z,
            ];
            (
                b.velocity.linear,
                b.inv_mass,
                if b.inertia > 0.0 {
                    1.0 / b.inertia
                } else {
                    0.0
                },
                pb,
            )
        } else {
            return;
        }
    } else {
        ([0.0; 3], 0.0, 0.0, contact.point)
    };
    let vrel_n = {
        let dvx = va[0] - vb[0];
        let dvy = va[1] - vb[1];
        let dvz = va[2] - vb[2];
        dvx * contact.normal[0] + dvy * contact.normal[1] + dvz * contact.normal[2]
    };
    if vrel_n >= 0.0 {
        return;
    }
    let ra = [
        contact.point[0] - pa[0],
        contact.point[1] - pa[1],
        contact.point[2] - pa[2],
    ];
    let rb = [
        contact.point[0] - pb[0],
        contact.point[1] - pb[1],
        contact.point[2] - pb[2],
    ];
    fn cross_n_sq(r: [Real; 3], n: [Real; 3], inv_i: Real) -> Real {
        let cx = r[1] * n[2] - r[2] * n[1];
        let cy = r[2] * n[0] - r[0] * n[2];
        let cz = r[0] * n[1] - r[1] * n[0];
        inv_i * (cx * cx + cy * cy + cz * cz)
    }
    let k = inv_ma
        + inv_mb
        + cross_n_sq(ra, contact.normal, inv_ia)
        + cross_n_sq(rb, contact.normal, inv_ib);
    if k < 1e-12 {
        return;
    }
    let e = contact.restitution;
    let j = -(1.0 + e) * vrel_n / k;
    if let Some(b) = bodies.get_mut(contact.body_a)
        && b.mode == BodyMode::Dynamic
    {
        b.velocity.linear[0] += j * inv_ma * contact.normal[0];
        b.velocity.linear[1] += j * inv_ma * contact.normal[1];
        b.velocity.linear[2] += j * inv_ma * contact.normal[2];
        let ax = ra[1] * contact.normal[2] - ra[2] * contact.normal[1];
        let ay = ra[2] * contact.normal[0] - ra[0] * contact.normal[2];
        let az = ra[0] * contact.normal[1] - ra[1] * contact.normal[0];
        b.velocity.angular[0] += inv_ia * ax * j;
        b.velocity.angular[1] += inv_ia * ay * j;
        b.velocity.angular[2] += inv_ia * az * j;
    }
    if let Some(bh) = contact.body_b
        && let Some(b) = bodies.get_mut(bh)
        && b.mode == BodyMode::Dynamic
    {
        b.velocity.linear[0] -= j * inv_mb * contact.normal[0];
        b.velocity.linear[1] -= j * inv_mb * contact.normal[1];
        b.velocity.linear[2] -= j * inv_mb * contact.normal[2];
        let bx = rb[1] * contact.normal[2] - rb[2] * contact.normal[1];
        let by = rb[2] * contact.normal[0] - rb[0] * contact.normal[2];
        let bz = rb[0] * contact.normal[1] - rb[1] * contact.normal[0];
        b.velocity.angular[0] -= inv_ib * bx * j;
        b.velocity.angular[1] -= inv_ib * by * j;
        b.velocity.angular[2] -= inv_ib * bz * j;
    }
}
/// Continuous collision detection between two moving spheres.
///
/// Returns the time of impact `t ∈ [0, 1]` (fraction of the step), or
/// `None` if no collision occurs during the step.
pub fn ccd_sphere_sphere(
    pos_a: [Real; 3],
    vel_a: [Real; 3],
    radius_a: Real,
    pos_b: [Real; 3],
    vel_b: [Real; 3],
    radius_b: Real,
) -> Option<Real> {
    let dv = [
        vel_b[0] - vel_a[0],
        vel_b[1] - vel_a[1],
        vel_b[2] - vel_a[2],
    ];
    let dp = [
        pos_b[0] - pos_a[0],
        pos_b[1] - pos_a[1],
        pos_b[2] - pos_a[2],
    ];
    let r = radius_a + radius_b;
    let a = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    let b = 2.0 * (dp[0] * dv[0] + dp[1] * dv[1] + dp[2] * dv[2]);
    let c = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - r * r;
    if a < 1e-20 {
        if c <= 0.0 {
            return Some(0.0);
        }
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t = (-b - sqrt_disc) / (2.0 * a);
    if (0.0..=1.0).contains(&t) {
        return Some(t);
    }
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    if (0.0..=1.0).contains(&t2) {
        return Some(t2);
    }
    None
}
/// CCD between a moving sphere and a static plane.
///
/// Returns the time of impact `t ∈ [0, 1]` or `None`.
pub fn ccd_sphere_plane(
    pos: [Real; 3],
    vel: [Real; 3],
    radius: Real,
    plane_normal: [Real; 3],
    plane_offset: Real,
) -> Option<Real> {
    let d0 = pos[0] * plane_normal[0] + pos[1] * plane_normal[1] + pos[2] * plane_normal[2]
        - plane_offset;
    let dv = vel[0] * plane_normal[0] + vel[1] * plane_normal[1] + vel[2] * plane_normal[2];
    if dv.abs() < 1e-12 {
        return if d0 <= radius { Some(0.0) } else { None };
    }
    let t = (radius - d0) / dv;
    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}
/// Analyse the constraint graph and return body groups (islands) sorted by size.
///
/// This is a standalone utility that works on a slice of constraints without
/// requiring access to the full `PhysicsWorld`.
pub fn constraint_islands(n_bodies: usize, constraints: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n_bodies).collect();
    fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for &(a, b) in constraints {
        if a < n_bodies && b < n_bodies {
            let ra = find(&mut parent, a);
            let rb = find(&mut parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
    }
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for i in 0..n_bodies {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }
    let mut result: Vec<Vec<usize>> = groups.into_values().collect();
    result.sort_by_key(|v: &Vec<usize>| std::cmp::Reverse(v.len()));
    result
}
#[cfg(test)]
mod world_expanded_tests {

    use crate::Body;
    use crate::BodyArena;
    use crate::BodyHandle;
    use crate::BodyMode;

    use crate::PhysicsWorld;
    use crate::Real;
    use crate::TimeStep;
    use crate::world::BodySnapshot;
    use crate::world::BvhNode;
    use crate::world::Contact;
    use crate::world::ForceField;
    use crate::world::SpatialHashGrid;
    use crate::world::ccd_sphere_plane;
    use crate::world::ccd_sphere_sphere;
    use crate::world::constraint_islands;
    use crate::world::resolve_contact;
    use crate::world::types::Frustum;
    #[test]
    fn test_spatial_hash_grid_candidate_pairs() {
        let mut grid = SpatialHashGrid::new(1.0);
        grid.insert(1, [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        grid.insert(2, [0.2, 0.2, 0.2], [0.7, 0.7, 0.7]);
        let pairs = grid.candidate_pairs();
        assert!(!pairs.is_empty(), "should find at least one candidate pair");
    }
    #[test]
    fn test_spatial_hash_grid_no_pairs() {
        let mut grid = SpatialHashGrid::new(1.0);
        grid.insert(1, [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]);
        grid.insert(2, [5.0, 5.0, 5.0], [5.5, 5.5, 5.5]);
        let pairs = grid.candidate_pairs();
        assert!(pairs.is_empty(), "should find no pairs for distant bodies");
    }
    #[test]
    fn test_spatial_hash_grid_clear() {
        let mut grid = SpatialHashGrid::new(1.0);
        grid.insert(1, [0.0; 3], [1.0; 3]);
        grid.clear();
        assert!(grid.candidate_pairs().is_empty());
    }
    #[test]
    fn test_bvh_build_and_query() {
        use crate::math::Vec3;
        use crate::types::Aabb;
        let aabbs: Vec<(Aabb, usize)> = (0..8)
            .map(|i| {
                let x = i as f64;
                (
                    Aabb::new(Vec3::new(x, 0.0, 0.0), Vec3::new(x + 0.5, 0.5, 0.5)),
                    i,
                )
            })
            .collect();
        let mut items: Vec<(Aabb, usize)> = aabbs;
        let bvh = BvhNode::build(&mut items).unwrap();
        let query = crate::types::Aabb::new(Vec3::new(2.0, 0.0, 0.0), Vec3::new(4.0, 1.0, 1.0));
        let mut results = Vec::new();
        bvh.query_aabb(&query, &mut results);
        assert!(!results.is_empty(), "BVH query should return results");
    }
    #[test]
    fn test_bvh_empty_returns_none() {
        let items: &mut [(crate::types::Aabb, usize)] = &mut [];
        assert!(BvhNode::build(items).is_none());
    }
    #[test]
    fn test_constant_force_field() {
        let ff = ForceField::constant([0.0, -9.81, 0.0]);
        let f = ff.force_on([0.0, 5.0, 0.0], [0.0; 3], 1.0);
        assert!((f[1] + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_drag_force_field_opposes_velocity() {
        let ff = ForceField::drag(0.5);
        let f = ff.force_on([0.0; 3], [2.0, 0.0, 0.0], 1.0);
        assert!(f[0] < 0.0, "drag should oppose motion");
    }
    #[test]
    fn test_point_gravity_attracts() {
        let ff = ForceField::point_gravity([0.0, 0.0, 0.0], 1.0, 0.01);
        let f = ff.force_on([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(f[0] < 0.0, "should attract toward origin: fx={}", f[0]);
    }
    #[test]
    fn test_inactive_force_field_returns_zero() {
        let mut ff = ForceField::constant([0.0, -9.81, 0.0]);
        ff.active = false;
        let f = ff.force_on([0.0; 3], [0.0; 3], 1.0);
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn test_resolve_contact_separates_bodies() {
        let mut arena = BodyArena::new();
        let mut b1 = Body::new_dynamic([0.0; 3], 1.0);
        let mut b2 = Body::new_dynamic([0.5, 0.0, 0.0], 1.0);
        b1.velocity.linear = [1.0, 0.0, 0.0];
        b2.velocity.linear = [-1.0, 0.0, 0.0];
        b1.linear_damping = 0.0;
        b2.linear_damping = 0.0;
        let h1 = arena.insert(b1);
        let h2 = arena.insert(b2);
        let contact = Contact::new(
            h1,
            Some(h2),
            [0.25, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.1,
            0.5,
            0.3,
        );
        resolve_contact(&mut arena, &contact);
        let va = arena.get(h1).unwrap().velocity.linear[0];
        assert!(
            va < 1.0,
            "body A velocity should decrease after contact: va={va}"
        );
    }
    #[test]
    fn test_ccd_sphere_sphere_approaching() {
        let t = ccd_sphere_sphere(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [3.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
        );
        assert!(t.is_some(), "approaching spheres should collide");
        let toi = t.unwrap();
        assert!((0.0..=1.0).contains(&toi), "toi out of [0,1]: {toi}");
    }
    #[test]
    fn test_ccd_sphere_sphere_not_colliding() {
        let t = ccd_sphere_sphere(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(
            t.is_none(),
            "parallel moving spheres far apart should not collide"
        );
    }
    #[test]
    fn test_ccd_sphere_plane() {
        let t = ccd_sphere_plane(
            [0.0, 1.4, 0.0],
            [0.0, -10.0, 0.0],
            0.5,
            [0.0, 1.0, 0.0],
            0.0,
        );
        assert!(
            t.is_some(),
            "sphere moving fast toward plane should collide within step"
        );
        let toi = t.unwrap();
        assert!((toi - 0.09).abs() < 1e-10, "toi={toi}");
    }
    #[test]
    fn test_ccd_sphere_plane_close() {
        let t = ccd_sphere_plane([0.0, 1.4, 0.0], [0.0, -2.0, 0.0], 0.5, [0.0, 1.0, 0.0], 0.0);
        assert!(t.is_some(), "sphere should hit plane within step");
        let toi = t.unwrap();
        assert!((toi - 0.45).abs() < 1e-10, "toi={toi}");
    }
    #[test]
    fn test_constraint_islands_two_groups() {
        let islands = constraint_islands(6, &[(0, 1), (1, 2), (3, 4)]);
        assert_eq!(
            islands.len(),
            3,
            "should have 3 islands (2 connected + 1 isolated)"
        );
    }
    #[test]
    fn test_constraint_islands_all_connected() {
        let islands = constraint_islands(4, &[(0, 1), (1, 2), (2, 3)]);
        assert_eq!(islands.len(), 1);
        assert_eq!(islands[0].len(), 4);
    }
    #[test]
    fn test_world_snapshot_restore() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        let snap = world.snapshot();
        world.step(&TimeStep::new(1.0 / 60.0));
        world.restore_snapshot(&snap);
        assert!((world.time - snap.time).abs() < 1e-12);
        let pos = world.body(h).unwrap().transform.position;
        assert!(
            pos.norm() < 1e-10,
            "position should be restored to zero: {:?}",
            pos
        );
    }
    #[test]
    fn test_world_apply_force_fields() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.body_mut(h).unwrap().linear_damping = 0.0;
        let fields = vec![ForceField::constant([10.0, 0.0, 0.0])];
        world.apply_force_fields(&fields);
        let f = world.body(h).unwrap().force[0];
        assert!((f - 10.0).abs() < 1e-10, "force should be 10.0: {f}");
    }
    #[test]
    fn test_world_damp_region() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.body_mut(h).unwrap().velocity.linear = [10.0, 0.0, 0.0];
        world.damp_region([0.0; 3], 5.0, 0.5);
        let v = world.body(h).unwrap().velocity.linear[0];
        assert!((v - 5.0).abs() < 1e-10, "velocity should be halved: {v}");
    }
    #[test]
    fn test_world_most_energetic_body() {
        let mut world = PhysicsWorld::new();
        let _h1 = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        let h2 = world.add_body(Body::new_dynamic([1.0, 0.0, 0.0], 1.0));
        world.body_mut(h2).unwrap().velocity.linear = [5.0, 0.0, 0.0];
        let best = world.most_energetic_body();
        assert_eq!(best, Some(h2), "h2 should be most energetic");
    }
    #[test]
    fn test_world_wake_all() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.body_mut(h).unwrap().mode = BodyMode::Sleeping;
        world.wake_all();
        assert_eq!(world.body(h).unwrap().mode, BodyMode::Dynamic);
    }
    #[test]
    fn test_world_diagnostics() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.add_body(Body::new_static([5.0, 0.0, 0.0]));
        let diag = world.diagnostics();
        assert_eq!(diag.total_bodies, 2);
        assert_eq!(diag.n_dynamic, 1);
        assert_eq!(diag.n_static, 1);
    }
    #[test]
    fn test_world_set_velocity() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.set_velocity(h, [3.0, 4.0, 0.0], [0.0; 3]);
        let v = &world.body(h).unwrap().velocity;
        assert!((v.linear[0] - 3.0).abs() < 1e-10);
        assert!((v.linear[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_world_set_position() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.set_position(h, [7.0, 8.0, 9.0]);
        let p = world.body(h).unwrap().transform.position;
        assert!((p.x - 7.0).abs() < 1e-10);
        assert!((p.y - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_world_angular_momentum() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        let l = world.total_angular_momentum();
        assert_eq!(l, [0.0; 3]);
    }
    #[test]
    fn test_world_apply_position_correction() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.apply_position_correction(h, [1.0, 2.0, 3.0]);
        let p = world.body(h).unwrap().transform.position;
        assert!((p.x - 1.0).abs() < 1e-10);
        assert!((p.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_world_body_snapshot() {
        let body = Body::new_dynamic([1.0, 2.0, 3.0], 5.0);
        let h = BodyHandle::new(0, 0);
        let snap = BodySnapshot::from_body(h, &body);
        assert!((snap.transform.position.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_query_sphere_finds_nearby() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([1.0, 0.0, 0.0], 1.0));
        let result = world.query_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(
            result.contains(&h),
            "body at distance 1 should be inside radius 2"
        );
    }
    #[test]
    fn test_query_sphere_excludes_distant() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([10.0, 0.0, 0.0], 1.0));
        let result = world.query_sphere([0.0, 0.0, 0.0], 2.0);
        assert!(
            !result.contains(&h),
            "body at distance 10 should be outside radius 2"
        );
    }
    #[test]
    fn test_query_sphere_empty_world() {
        let world = PhysicsWorld::new();
        let result = world.query_sphere([0.0, 0.0, 0.0], 100.0);
        assert!(result.is_empty());
    }
    #[test]
    fn test_query_frustum_includes_origin_body() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([0.0, 0.0, 0.0], 1.0));
        let big = 10.0f64;
        let planes: [[Real; 4]; 6] = [
            [1.0, 0.0, 0.0, big],
            [-1.0, 0.0, 0.0, big],
            [0.0, 1.0, 0.0, big],
            [0.0, -1.0, 0.0, big],
            [0.0, 0.0, 1.0, big],
            [0.0, 0.0, -1.0, big],
        ];
        let frustum = Frustum::new(planes);
        let result = world.query_frustum(&frustum);
        assert!(
            result.contains(&h),
            "origin body should be inside large frustum"
        );
    }
    #[test]
    fn test_query_frustum_excludes_far_body() {
        let mut world = PhysicsWorld::new();
        let h = world.add_body(Body::new_dynamic([100.0, 100.0, 100.0], 1.0));
        let s = 1.0f64;
        let planes: [[Real; 4]; 6] = [
            [1.0, 0.0, 0.0, s],
            [-1.0, 0.0, 0.0, s],
            [0.0, 1.0, 0.0, s],
            [0.0, -1.0, 0.0, s],
            [0.0, 0.0, 1.0, s],
            [0.0, 0.0, -1.0, s],
        ];
        let frustum = Frustum::new(planes);
        let result = world.query_frustum(&frustum);
        assert!(
            !result.contains(&h),
            "far body should be outside tight frustum"
        );
    }
    #[test]
    fn test_compute_stats_body_count() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([0.0; 3], 1.0));
        world.add_body(Body::new_static([1.0, 0.0, 0.0]));
        let stats = world.compute_stats();
        assert_eq!(stats.body_count, 2);
    }
    #[test]
    fn test_compute_stats_zero_time_on_new_world() {
        let world = PhysicsWorld::new();
        let stats = world.compute_stats();
        assert_eq!(stats.simulation_time, 0.0);
    }
    #[test]
    fn test_serialize_snapshot_contains_time() {
        let world = PhysicsWorld::new();
        let snap = world.serialize_snapshot();
        assert!(
            snap.lines[0].starts_with("time="),
            "first line should start with 'time='"
        );
    }
    #[test]
    fn test_serialize_snapshot_body_entry() {
        let mut world = PhysicsWorld::new();
        world.add_body(Body::new_dynamic([3.0, 0.0, 0.0], 2.0));
        let snap = world.serialize_snapshot();
        let has_body_line = snap.lines.iter().any(|l| l.starts_with("body "));
        assert!(has_body_line, "snapshot should contain a body line");
    }
}
