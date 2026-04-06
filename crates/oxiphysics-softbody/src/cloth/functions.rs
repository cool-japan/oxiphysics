//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use crate::particle::SoftParticle;
use oxiphysics_core::math::Vec3;

use super::types::{ClothEdge, ClothVertex, RigidSphere};

/// Apply aerodynamic wind force to particles based on triangle normals and
/// areas.
///
/// For each triangle the wind force is proportional to the projected area
/// facing the wind direction, distributed equally among the three vertices.
pub fn apply_wind(particles: &mut [SoftParticle], wind: &Vec3, triangles: &[[usize; 3]]) {
    for tri in triangles {
        let p0 = particles[tri[0]].position;
        let p1 = particles[tri[1]].position;
        let p2 = particles[tri[2]].position;
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let normal = e1.cross(&e2);
        let area_x2 = normal.norm();
        if area_x2 < 1e-12 {
            continue;
        }
        let n_hat = normal / area_x2;
        let area = area_x2 * 0.5;
        let cos_angle = n_hat.dot(wind);
        let force = n_hat * (cos_angle * area / 3.0);
        for &vi in tri {
            if !particles[vi].is_static() {
                particles[vi].external_force += force;
            }
        }
    }
}
#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = length(a);
    if l < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / l)
    }
}
/// Resolve self-collision for cloth vertices using a simple repulsion sphere.
///
/// For every pair of vertices closer than `min_dist`, applies a positional
/// correction that pushes them apart to `min_dist`.  Only non-fixed vertices
/// are corrected.
#[allow(dead_code)]
pub fn resolve_self_collision(vertices: &mut Vec<ClothVertex>, min_dist: f64) {
    let n = vertices.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let pi = vertices[i].pos;
            let pj = vertices[j].pos;
            let d = sub(pj, pi);
            let dist = length(d);
            if dist < min_dist && dist > 1e-14 {
                let overlap = min_dist - dist;
                let dir = scale(d, 1.0 / dist);
                let wi = vertices[i].inv_mass;
                let wj = vertices[j].inv_mass;
                let w_sum = wi + wj;
                if w_sum < 1e-30 {
                    continue;
                }
                let corr = scale(dir, overlap / w_sum);
                if !vertices[i].fixed {
                    vertices[i].pos = sub(vertices[i].pos, scale(corr, wi));
                }
                if !vertices[j].fixed {
                    vertices[j].pos = add(vertices[j].pos, scale(corr, wj));
                }
            }
        }
    }
}
/// Resolve collision between all cloth vertices and a rigid sphere.
///
/// Vertices that penetrate the sphere are pushed to its surface and their
/// inward velocity components are removed (or reflected with restitution).
#[allow(dead_code)]
pub fn resolve_cloth_sphere_collision(vertices: &mut Vec<ClothVertex>, sphere: &RigidSphere) {
    for v in vertices.iter_mut() {
        if v.fixed {
            continue;
        }
        let d = sub(v.pos, sphere.center);
        let dist = length(d);
        if dist < sphere.radius && dist > 1e-14 {
            let n_hat = scale(d, 1.0 / dist);
            v.pos = add(sphere.center, scale(n_hat, sphere.radius + 1e-6));
            let vn = dot(v.vel, n_hat);
            if vn < 0.0 {
                v.vel = sub(v.vel, scale(n_hat, (1.0 + sphere.restitution) * vn));
            }
        }
    }
}
/// Resolve collision between all cloth vertices and a rigid axis-aligned box.
///
/// Vertices below `floor_y` are pushed back up with restitution applied to the
/// y-velocity.
#[allow(dead_code)]
pub fn resolve_cloth_floor_collision(
    vertices: &mut Vec<ClothVertex>,
    floor_y: f64,
    restitution: f64,
) {
    for v in vertices.iter_mut() {
        if v.fixed {
            continue;
        }
        if v.pos[1] < floor_y {
            v.pos[1] = floor_y;
            if v.vel[1] < 0.0 {
                v.vel[1] *= -restitution;
            }
        }
    }
}
/// Lift force model for cloth triangles.
///
/// Computes an aerodynamic lift force using a flat-plate approximation:
/// `F_lift = 0.5 * ρ_air * C_L * A * |v_rel|² * n̂`
/// where `v_rel` is the relative wind velocity and `C_L` depends on the
/// angle of attack.
#[allow(dead_code)]
pub fn compute_lift_force(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    v_tri: [f64; 3],
    wind: [f64; 3],
    air_density: f64,
    lift_coeff: f64,
) -> [f64; 3] {
    let e1 = sub(p1, p0);
    let e2 = sub(p2, p0);
    let area_normal = cross(e1, e2);
    let area2 = length(area_normal);
    if area2 < 1e-14 {
        return [0.0; 3];
    }
    let n_hat = scale(area_normal, 1.0 / area2);
    let area = area2 * 0.5;
    let v_rel = sub(wind, v_tri);
    let v_rel_sq = dot(v_rel, v_rel);
    let sin_alpha = dot(normalize(v_rel), n_hat).abs().min(1.0);
    let f_mag = 0.5 * air_density * lift_coeff * area * v_rel_sq * sin_alpha;
    scale(n_hat, f_mag)
}
/// Wrinkle amplitude at a cloth vertex.
///
/// Estimates the wrinkling intensity from the compressive strain at the
/// vertex, expressed as the ratio of compressive edge stretch to a reference
/// rest length.
#[allow(dead_code)]
pub fn vertex_wrinkle_amplitude(
    vertex_idx: usize,
    edges: &[ClothEdge],
    vertices: &[ClothVertex],
    wrinkle_scale: f64,
) -> f64 {
    let mut compressive_strain = 0.0_f64;
    let mut count = 0usize;
    for e in edges {
        if e.torn {
            continue;
        }
        if e.a == vertex_idx || e.b == vertex_idx {
            let pa = vertices[e.a].pos;
            let pb = vertices[e.b].pos;
            let cur = length(sub(pb, pa));
            let strain = (cur - e.rest_length) / e.rest_length.max(1e-12);
            if strain < 0.0 {
                compressive_strain += (-strain).abs();
                count += 1;
            }
        }
    }
    if count == 0 {
        return 0.0;
    }
    wrinkle_scale * compressive_strain / count as f64
}
/// Displace cloth vertices along their normals by a wrinkle amplitude.
///
/// This is a post-processing wrinkling enhancement: each vertex is displaced
/// by `amplitude(v) * normal(v)`.
#[allow(dead_code)]
pub fn apply_wrinkling(vertices: &mut Vec<ClothVertex>, edges: &[ClothEdge], wrinkle_scale: f64) {
    let n = vertices.len();
    let mut amplitudes = vec![0.0_f64; n];
    for i in 0..n {
        amplitudes[i] = vertex_wrinkle_amplitude(i, edges, vertices, wrinkle_scale);
    }
    for (i, v) in vertices.iter_mut().enumerate() {
        if v.fixed {
            continue;
        }
        let amp = amplitudes[i];
        v.pos = add(v.pos, scale(v.normal, amp));
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ClothMesh;

    use crate::XpbdClothMesh;

    use crate::cloth::length;
    use crate::cloth::sub;

    #[test]
    fn test_new_grid_vertex_count() {
        let nx = 5;
        let ny = 4;
        let cloth = ClothMesh::new_grid(nx, ny, 3.0, 2.0, 1.0);
        assert_eq!(cloth.vertices.len(), nx * ny);
    }
    #[test]
    fn test_pin_corner_fixes_vertex() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.pin_corner(0);
        assert!(cloth.vertices[0].fixed);
        assert_eq!(cloth.vertices[0].inv_mass, 0.0);
    }
    #[test]
    fn test_step_no_panic() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        for _ in 0..10 {
            cloth.step(1.0 / 60.0, 5);
        }
    }
    #[test]
    fn test_gravity_free_vertex_falls() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            cloth.step(dt, 5);
        }
        let y = cloth.vertices[3].pos[1];
        assert!(
            y < -0.01,
            "Free vertex should fall under gravity, got y={y}"
        );
    }
    #[test]
    fn test_tearing_overstretched_edge() {
        let mut cloth = ClothMesh::new_grid(2, 2, 1.0, 1.0, 1.0);
        cloth.edges[0].rest_length = 0.01;
        cloth.edges[0].tear_threshold = 0.1;
        let b = cloth.edges[0].b;
        cloth.vertices[b].pos = [100.0, 0.0, 0.0];
        cloth.check_tearing();
        assert!(cloth.edges[0].torn, "Overstretched edge should be torn");
    }
    #[test]
    fn test_compute_vertex_normals_no_panic() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.compute_vertex_normals();
        for v in &cloth.vertices {
            let l = length(v.normal);
            assert!(l < 1.0 + 1e-6, "Normal length should be <= 1, got {l}");
        }
    }
    #[test]
    fn test_kinetic_energy_non_negative() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        for _ in 0..10 {
            cloth.step(1.0 / 60.0, 3);
        }
        let ke = cloth.kinetic_energy();
        assert!(ke >= 0.0, "Kinetic energy must be non-negative, got {ke}");
    }
    #[test]
    fn test_xpbd_cloth_mesh_topology() {
        let cloth = XpbdClothMesh::new(4, 3, 3.0, 2.0, 1.0, 0.001);
        assert_eq!(cloth.num_particles(), 12);
        assert_eq!(cloth.num_triangles(), 12);
        assert!(!cloth.distance_constraints.is_empty());
        assert!(!cloth.bending_constraints.is_empty());
    }
    #[test]
    fn test_self_collision_pushes_apart() {
        let mut cloth = ClothMesh::new_grid(2, 2, 1.0, 1.0, 1.0);
        cloth.vertices[0].pos = [0.0, 0.0, 0.0];
        cloth.vertices[1].pos = [0.001, 0.0, 0.0];
        let min_dist = 0.1;
        cloth.resolve_self_collision(min_dist);
        let d = length(sub(cloth.vertices[1].pos, cloth.vertices[0].pos));
        assert!(
            d >= min_dist - 1e-10,
            "After self-collision, distance should be >= min_dist: {d}"
        );
    }
    #[test]
    fn test_sphere_collision_pushes_out() {
        let mut cloth = ClothMesh::new_grid(2, 2, 1.0, 1.0, 1.0);
        let sphere = RigidSphere::new([0.0, 0.0, 0.0], 1.0, 0.5);
        cloth.vertices[0].pos = [0.1, 0.0, 0.0];
        cloth.collide_with_sphere(&sphere);
        let d = length(sub(cloth.vertices[0].pos, sphere.center));
        assert!(
            d >= sphere.radius - 1e-10,
            "Vertex should be pushed outside sphere, d={d}"
        );
    }
    #[test]
    fn test_floor_collision_pushes_up() {
        let mut cloth = ClothMesh::new_grid(2, 2, 1.0, 1.0, 1.0);
        cloth.vertices[0].pos[1] = -0.5;
        cloth.vertices[0].vel[1] = -1.0;
        cloth.collide_with_floor(0.0, 0.5);
        assert!(
            cloth.vertices[0].pos[1] >= 0.0 - 1e-14,
            "Floor collision should move vertex above floor"
        );
        assert!(
            cloth.vertices[0].vel[1] >= 0.0,
            "Floor collision should reverse downward velocity"
        );
    }
    #[test]
    fn test_wrinkling_no_panic() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        for _ in 0..10 {
            cloth.step(1.0 / 60.0, 3);
        }
        cloth.apply_wrinkling(0.01);
    }
    #[test]
    fn test_repair_tears() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.edges[0].torn = true;
        cloth.repair_tears(100.0);
        assert!(
            !cloth.edges[0].torn,
            "Edge should be repaired at high repair_ratio"
        );
    }
    #[test]
    fn test_force_tear_edges() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        assert!(!cloth.edges[0].torn);
        cloth.force_tear_edges(&[0, 1]);
        assert!(cloth.edges[0].torn);
        assert!(cloth.edges[1].torn);
    }
    #[test]
    fn test_torn_edge_indices() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.edges[2].torn = true;
        cloth.edges[4].torn = true;
        let torn = cloth.torn_edge_indices();
        assert_eq!(torn, vec![2, 4]);
    }
    #[test]
    fn test_total_area_positive() {
        let cloth = ClothMesh::new_grid(4, 4, 2.0, 2.0, 1.0);
        let area = cloth.total_area();
        assert!(area > 0.0, "Total area should be positive, got {area}");
    }
    #[test]
    fn test_bounding_box_bounds_vertices() {
        let cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        let (lo, hi) = cloth.bounding_box();
        for v in &cloth.vertices {
            for k in 0..3 {
                assert!(
                    v.pos[k] >= lo[k] - 1e-12 && v.pos[k] <= hi[k] + 1e-12,
                    "Vertex pos[{k}]={} outside bounding box [{}, {}]",
                    v.pos[k],
                    lo[k],
                    hi[k]
                );
            }
        }
    }
    #[test]
    fn test_step_with_collision_no_panic() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        let sphere = RigidSphere::new([0.5, -0.5, 0.5], 0.3, 0.5);
        for _ in 0..10 {
            cloth.step_with_collision(1.0 / 60.0, 3, Some(0.05), Some(&sphere), Some(-2.0));
        }
    }
    #[test]
    fn test_lift_no_panic() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.wind = [0.0, 0.0, 5.0];
        cloth.apply_lift(1.0 / 60.0, 1.2);
    }
    #[test]
    fn test_max_vertex_speed_non_negative() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        for _ in 0..5 {
            cloth.step(1.0 / 60.0, 3);
        }
        assert!(cloth.max_vertex_speed() >= 0.0);
    }
}
/// Compute the dihedral angle between two triangles sharing edge (p0, p1).
///
/// `pk` and `pl` are the "wing" vertices of each triangle.  Returns the angle
/// in radians in \[0, π\].
#[inline]
pub(super) fn compute_dihedral_angle(
    p0: [f64; 3],
    p1: [f64; 3],
    pk: [f64; 3],
    pl: [f64; 3],
) -> f64 {
    let e = sub(p1, p0);
    let e_len = length(e);
    if e_len < 1e-14 {
        return 0.0;
    }
    let e_hat = scale(e, 1.0 / e_len);
    let ek = sub(pk, p0);
    let el = sub(pl, p0);
    let nk = cross(e_hat, normalize(ek));
    let nl = cross(e_hat, normalize(el));
    let cos_theta = dot(nk, nl).clamp(-1.0, 1.0);
    cos_theta.acos()
}
#[cfg(test)]
mod pbd_cloth_tests {

    use crate::PbdClothMesh;

    use crate::cloth::length;
    use crate::cloth::sub;

    #[test]
    fn test_pbd_grid_vertex_count() {
        let nx = 5;
        let ny = 4;
        let cloth = PbdClothMesh::grid(nx, ny, 0.1, 1.0);
        assert_eq!(cloth.positions.len(), nx * ny);
        assert_eq!(cloth.velocities.len(), nx * ny);
        assert_eq!(cloth.masses.len(), nx * ny);
    }
    #[test]
    fn test_pbd_pinned_particle_does_not_move() {
        let mut cloth = PbdClothMesh::grid(4, 4, 0.1, 1.0);
        cloth.pin(0);
        let pos_before = cloth.positions[0];
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..60 {
            cloth.step(1.0 / 60.0, gravity, 3);
        }
        let dx = length(sub(cloth.positions[0], pos_before));
        assert!(dx < 1e-12, "Pinned particle should not move, moved by {dx}");
    }
    #[test]
    fn test_pbd_sphere_collision_pushes_out() {
        let mut cloth = PbdClothMesh::grid(3, 3, 1.0, 1.0);
        cloth.positions[4] = [0.1, 0.0, 0.0];
        let center = [0.0_f64; 3];
        let r = 0.5;
        cloth.apply_sphere_collision(center, r);
        let dist = length(sub(cloth.positions[4], center));
        assert!(
            dist >= r - 1e-9,
            "Particle should be outside sphere after collision, dist={dist}"
        );
    }
    #[test]
    fn test_pbd_step_moves_free_particle_under_gravity() {
        let mut cloth = PbdClothMesh::grid(3, 3, 1.0, 1.0);
        for c in 0..3 {
            cloth.pin(c);
        }
        let init_y = cloth.positions[3][1];
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..30 {
            cloth.step(1.0 / 60.0, gravity, 3);
        }
        let final_y = cloth.positions[3][1];
        assert!(final_y < init_y, "Free particle should fall under gravity");
    }
    #[test]
    fn test_pbd_total_area_positive() {
        let cloth = PbdClothMesh::grid(4, 4, 0.5, 1.0);
        let area = cloth.total_area();
        assert!(area > 0.0, "Total area should be positive, got {area}");
    }
}
/// Green–Lagrange strain tensor components from triangle deformation.
///
/// Given a triangle in 3-D defined by vertices `p0`, `p1`, `p2` and its
/// reference (rest) shape `q0`, `q1`, `q2` lying in the XY plane, returns
/// the Green–Lagrange strain `E = 0.5 * (F^T F - I)` in Voigt notation
/// `[E11, E22, E12]`.
///
/// The deformation gradient `F` maps the reference triangle to the current
/// one.
#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn green_lagrange_strain(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    q0: [f64; 2],
    q1: [f64; 2],
    q2: [f64; 2],
) -> [f64; 3] {
    let dq1 = [q1[0] - q0[0], q1[1] - q0[1]];
    let dq2 = [q2[0] - q0[0], q2[1] - q0[1]];
    let det = dq1[0] * dq2[1] - dq1[1] * dq2[0];
    if det.abs() < 1e-20 {
        return [0.0; 3];
    }
    let inv_det = 1.0 / det;
    let inv_dq = [
        [dq2[1] * inv_det, -dq2[0] * inv_det],
        [-dq1[1] * inv_det, dq1[0] * inv_det],
    ];
    let dp1 = sub(p1, p0);
    let dp2 = sub(p2, p0);
    let f0 = [
        dp1[0] * inv_dq[0][0] + dp2[0] * inv_dq[1][0],
        dp1[0] * inv_dq[0][1] + dp2[0] * inv_dq[1][1],
    ];
    let f1 = [
        dp1[1] * inv_dq[0][0] + dp2[1] * inv_dq[1][0],
        dp1[1] * inv_dq[0][1] + dp2[1] * inv_dq[1][1],
    ];
    let f2 = [
        dp1[2] * inv_dq[0][0] + dp2[2] * inv_dq[1][0],
        dp1[2] * inv_dq[0][1] + dp2[2] * inv_dq[1][1],
    ];
    let ftf = [
        f0[0] * f0[0] + f1[0] * f1[0] + f2[0] * f2[0],
        f0[1] * f0[1] + f1[1] * f1[1] + f2[1] * f2[1],
        f0[0] * f0[1] + f1[0] * f1[1] + f2[0] * f2[1],
    ];
    [0.5 * (ftf[0] - 1.0), 0.5 * (ftf[1] - 1.0), 0.5 * ftf[2]]
}
/// St. Venant–Kirchhoff strain energy for a triangle.
///
/// `W = λ/2 * tr(E)^2 + μ * tr(E^2)`
///
/// where λ and μ are the Lamé parameters and E is the Green–Lagrange strain.
///
/// `strain_voigt` = `[E11, E22, E12]`.
#[allow(dead_code)]
pub fn svk_energy(strain_voigt: [f64; 3], lambda: f64, mu: f64) -> f64 {
    let e11 = strain_voigt[0];
    let e22 = strain_voigt[1];
    let e12 = strain_voigt[2];
    let tr_e = e11 + e22;
    let tr_e2 = e11 * e11 + e22 * e22 + 2.0 * e12 * e12;
    0.5 * lambda * tr_e * tr_e + mu * tr_e2
}
/// St. Venant–Kirchhoff second Piola–Kirchhoff stress (Voigt notation).
///
/// `S = λ * tr(E) * I + 2μ * E`
///
/// Returns `[S11, S22, S12]`.
#[allow(dead_code)]
pub fn svk_stress(strain_voigt: [f64; 3], lambda: f64, mu: f64) -> [f64; 3] {
    let e11 = strain_voigt[0];
    let e22 = strain_voigt[1];
    let e12 = strain_voigt[2];
    let tr_e = e11 + e22;
    [
        lambda * tr_e + 2.0 * mu * e11,
        lambda * tr_e + 2.0 * mu * e22,
        2.0 * mu * e12,
    ]
}
/// Convert Young's modulus `E` and Poisson's ratio `nu` to 2-D plane-stress
/// Lamé parameters.
#[allow(dead_code)]
pub fn lame_from_young_poisson(young: f64, nu: f64) -> (f64, f64) {
    let mu = young / (2.0 * (1.0 + nu));
    let lambda = young * nu / ((1.0 + nu) * (1.0 - nu));
    (lambda, mu)
}
/// Detect compressive regions in a cloth mesh.
///
/// Returns a list of `(vertex_index, compression_measure)` pairs for all
/// vertices where the average edge strain is negative (compressive).
///
/// The compression measure is the magnitude of the average compressive strain.
#[allow(dead_code)]
pub fn detect_wrinkle_regions(
    vertices: &[ClothVertex],
    edges: &[ClothEdge],
    compression_threshold: f64,
) -> Vec<(usize, f64)> {
    let n = vertices.len();
    let mut sums = vec![0.0_f64; n];
    let mut counts = vec![0usize; n];
    for e in edges {
        if e.torn {
            continue;
        }
        let pa = vertices[e.a].pos;
        let pb = vertices[e.b].pos;
        let d = length(sub(pb, pa));
        let strain = (d - e.rest_length) / e.rest_length.max(1e-12);
        if strain < 0.0 {
            sums[e.a] += strain.abs();
            sums[e.b] += strain.abs();
            counts[e.a] += 1;
            counts[e.b] += 1;
        }
    }
    let mut result = Vec::new();
    for i in 0..n {
        if counts[i] > 0 {
            let avg = sums[i] / counts[i] as f64;
            if avg >= compression_threshold {
                result.push((i, avg));
            }
        }
    }
    result
}
/// Propagate a tear through a cloth mesh starting from a given torn edge.
///
/// Starting from torn edges, the algorithm looks at adjacent (sharing a vertex)
/// edges and tears them if their current strain exceeds `propagation_threshold`.
///
/// Returns the number of newly torn edges.
#[allow(dead_code)]
pub fn propagate_tears(
    edges: &mut Vec<ClothEdge>,
    vertices: &[ClothVertex],
    propagation_threshold: f64,
    max_propagation_steps: usize,
) -> usize {
    let mut newly_torn = 0;
    for _ in 0..max_propagation_steps {
        let mut torn_this_step = 0;
        let n_edges = edges.len();
        let mut torn_vertices: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for e in edges.iter() {
            if e.torn {
                torn_vertices.insert(e.a);
                torn_vertices.insert(e.b);
            }
        }
        for i in 0..n_edges {
            if edges[i].torn {
                continue;
            }
            let a = edges[i].a;
            let b = edges[i].b;
            if !torn_vertices.contains(&a) && !torn_vertices.contains(&b) {
                continue;
            }
            let pa = vertices[a].pos;
            let pb = vertices[b].pos;
            let d = length(sub(pb, pa));
            let strain = (d - edges[i].rest_length) / edges[i].rest_length.max(1e-12);
            if strain > propagation_threshold {
                edges[i].torn = true;
                torn_this_step += 1;
            }
        }
        newly_torn += torn_this_step;
        if torn_this_step == 0 {
            break;
        }
    }
    newly_torn
}
#[cfg(test)]
mod cloth_extended_tests {
    use super::*;
    use crate::ClothMesh;

    use crate::SeamConstraint;
    use crate::SeamSystem;

    use crate::cloth::add;
    use crate::cloth::length;
    use crate::cloth::sub;
    use crate::detect_wrinkle_regions;
    use crate::lame_from_young_poisson;
    use crate::propagate_tears;
    use crate::svk_energy;
    use crate::svk_stress;
    #[test]
    fn test_seam_weld_pulls_vertices_together() {
        let mut seam = SeamConstraint::new_weld(0, 1, 1.0);
        let mut positions = vec![[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        seam.apply(&positions, &inv_masses);
        let (da, db) = seam.apply(&positions, &inv_masses);
        assert!(
            da[0] > 0.0 || da[0].abs() < 1e-10,
            "vertex 0 correction should be toward vertex 1"
        );
        let _ = (da, db);
        positions[0] = add(positions[0], da);
        positions[1] = add(positions[1], db);
        let d = length(sub(positions[1], positions[0]));
        assert!(
            d < 1.0 - 1e-10,
            "vertices should be closer after seam constraint, d={d}"
        );
    }
    #[test]
    fn test_seam_elastic_stretch_ratio() {
        let seam = SeamConstraint::new_elastic(0, 1, 1.0, 1.0, 0.5);
        let positions = [[0.0_f64; 3], [1.5, 0.0, 0.0]];
        let ratio = seam.stretch_ratio(&positions);
        assert!(
            (ratio - 0.5).abs() < 1e-10,
            "stretch ratio should be 0.5, got {ratio}"
        );
    }
    #[test]
    fn test_seam_tears_when_threshold_exceeded() {
        let mut seam = SeamConstraint::new_elastic(0, 1, 1.0, 1.0, 0.1);
        let positions = [[0.0_f64; 3], [2.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let _ = seam.apply(&positions, &inv_masses);
        assert!(seam.torn, "seam should tear when strain > threshold");
    }
    #[test]
    fn test_seam_system_active_count() {
        let mut sys = SeamSystem::new();
        sys.add_seam(SeamConstraint::new_weld(0, 1, 1.0));
        sys.add_seam(SeamConstraint::new_weld(1, 2, 1.0));
        assert_eq!(sys.active_count(), 2);
        sys.seams[0].torn = true;
        assert_eq!(sys.active_count(), 1);
        assert_eq!(sys.torn_count(), 1);
    }
    #[test]
    fn test_seam_system_solve_moves_vertices() {
        let mut sys = SeamSystem::new();
        sys.add_seam(SeamConstraint::new_weld(0, 1, 1.0));
        let mut positions = vec![[0.0_f64; 3], [2.0, 0.0, 0.0]];
        let inv_masses = vec![1.0, 1.0];
        let d_before = length(sub(positions[1], positions[0]));
        sys.solve(&mut positions, &inv_masses);
        let d_after = length(sub(positions[1], positions[0]));
        assert!(
            d_after < d_before,
            "seam system should pull vertices together"
        );
    }
    #[test]
    fn test_seam_weld_both_pinned_no_correction() {
        let mut seam = SeamConstraint::new_weld(0, 1, 1.0);
        let positions = [[0.0_f64; 3], [1.0, 0.0, 0.0]];
        let inv_masses = [0.0, 0.0];
        let (da, db) = seam.apply(&positions, &inv_masses);
        assert_eq!(da, [0.0; 3], "no correction when both pinned");
        assert_eq!(db, [0.0; 3], "no correction when both pinned");
    }
    #[test]
    fn test_svk_energy_zero_at_rest() {
        let strain = [0.0, 0.0, 0.0];
        let (lambda, mu) = lame_from_young_poisson(1e6, 0.3);
        let e = svk_energy(strain, lambda, mu);
        assert!(e.abs() < 1e-20, "zero strain → zero SVK energy, got {e}");
    }
    #[test]
    fn test_svk_energy_positive_under_stretch() {
        let (lambda, mu) = lame_from_young_poisson(1e6, 0.3);
        let strain = [0.1, 0.0, 0.0];
        let e = svk_energy(strain, lambda, mu);
        assert!(
            e > 0.0,
            "SVK energy should be positive under stretch, got {e}"
        );
    }
    #[test]
    fn test_svk_stress_zero_at_zero_strain() {
        let (lambda, mu) = lame_from_young_poisson(1e6, 0.3);
        let stress = svk_stress([0.0; 3], lambda, mu);
        assert!(stress[0].abs() < 1e-20);
        assert!(stress[1].abs() < 1e-20);
        assert!(stress[2].abs() < 1e-20);
    }
    #[test]
    fn test_svk_stress_positive_under_tension() {
        let (lambda, mu) = lame_from_young_poisson(1e6, 0.3);
        let strain = [0.05, 0.0, 0.0];
        let stress = svk_stress(strain, lambda, mu);
        assert!(
            stress[0] > 0.0,
            "S11 should be positive under tension, got {}",
            stress[0]
        );
    }
    #[test]
    fn test_green_lagrange_strain_identity_deformation() {
        let q0 = [0.0_f64, 0.0];
        let q1 = [1.0_f64, 0.0];
        let q2 = [0.0_f64, 1.0];
        let p0 = [0.0, 0.0, 0.0_f64];
        let p1 = [1.0, 0.0, 0.0_f64];
        let p2 = [0.0, 1.0, 0.0_f64];
        let strain = green_lagrange_strain(p0, p1, p2, q0, q1, q2);
        for (k, &s) in strain.iter().enumerate() {
            assert!(
                s.abs() < 1e-10,
                "strain[{k}] should be 0 for identity, got {s}"
            );
        }
    }
    #[test]
    fn test_green_lagrange_strain_uniaxial_stretch() {
        let q0 = [0.0_f64, 0.0];
        let q1 = [1.0_f64, 0.0];
        let q2 = [0.0_f64, 1.0];
        let p0 = [0.0, 0.0, 0.0_f64];
        let p1 = [1.5, 0.0, 0.0_f64];
        let p2 = [0.0, 1.0, 0.0_f64];
        let strain = green_lagrange_strain(p0, p1, p2, q0, q1, q2);
        assert!(
            (strain[0] - 0.625).abs() < 1e-10,
            "E11 should be 0.625 for 1.5× stretch, got {}",
            strain[0]
        );
    }
    #[test]
    fn test_lame_from_young_poisson_steel() {
        let (lambda, mu) = lame_from_young_poisson(200e9, 0.3);
        assert!(
            (mu - 76.923e9).abs() / 76.923e9 < 0.01,
            "mu for steel should be ~76.9 GPa, got {mu}"
        );
        assert!(lambda > 0.0, "lambda should be positive");
    }
    #[test]
    fn test_detect_wrinkle_regions_no_wrinkles_when_stretched() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        let wrinkles = detect_wrinkle_regions(&cloth.vertices, &cloth.edges, 0.001);
        assert!(
            wrinkles.is_empty(),
            "no wrinkles when cloth is at rest or stretched"
        );
        cloth.pin_top_edge();
        for _ in 0..5 {
            cloth.step(1.0 / 60.0, 3);
        }
        let _wrinkles2 = detect_wrinkle_regions(&cloth.vertices, &cloth.edges, 0.001);
    }
    #[test]
    fn test_detect_wrinkle_regions_finds_compressed_edges() {
        let mut cloth = ClothMesh::new_grid(2, 2, 1.0, 1.0, 1.0);
        let b_idx = cloth.edges[0].b;
        cloth.vertices[b_idx].pos = cloth.vertices[cloth.edges[0].a].pos;
        cloth.vertices[b_idx].pos[0] += 0.01;
        let wrinkles = detect_wrinkle_regions(&cloth.vertices, &cloth.edges, 0.0);
        assert!(
            !wrinkles.is_empty(),
            "compressed edge should produce wrinkle region"
        );
    }
    #[test]
    fn test_propagate_tears_from_initial_tear() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.edges[0].torn = true;
        let b1 = cloth.edges[1].b;
        cloth.vertices[b1].pos[0] += 5.0;
        let newly_torn = propagate_tears(&mut cloth.edges, &cloth.vertices, 0.01, 3);
        assert!(
            newly_torn > 0 || cloth.edges[1].torn,
            "tear should propagate to adjacent overstretched edges"
        );
    }
    #[test]
    fn test_propagate_tears_no_propagation_when_unstretched() {
        let mut cloth = ClothMesh::new_grid(4, 4, 1.0, 1.0, 1.0);
        cloth.edges[0].torn = true;
        let count = propagate_tears(&mut cloth.edges, &cloth.vertices, 100.0, 3);
        assert_eq!(count, 0, "no propagation for high threshold, got {count}");
    }
    #[test]
    fn test_propagate_tears_stops_after_max_steps() {
        let mut cloth = ClothMesh::new_grid(5, 5, 1.0, 1.0, 1.0);
        cloth.edges[0].torn = true;
        for i in 1..cloth.edges.len() {
            let b_idx = cloth.edges[i].b;
            cloth.vertices[b_idx].pos[0] += 100.0;
        }
        let _count = propagate_tears(&mut cloth.edges, &cloth.vertices, 0.001, 2);
    }
    #[test]
    fn test_cloth_seam_and_step_combined() {
        let mut cloth = ClothMesh::new_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.pin_top_edge();
        for _ in 0..20 {
            cloth.step(1.0 / 60.0, 3);
        }
        let bot_y = cloth.vertices[cloth.vertices.len() - 1].pos[1];
        assert!(
            bot_y < -0.01,
            "cloth bottom should have fallen, got {bot_y}"
        );
    }
    #[test]
    fn test_svk_energy_symmetry() {
        let (lambda, mu) = lame_from_young_poisson(1e6, 0.3);
        let e1 = svk_energy([0.1, 0.05, 0.0], lambda, mu);
        let e2 = svk_energy([0.05, 0.1, 0.0], lambda, mu);
        assert!(
            (e1 - e2).abs() < 1e-10,
            "SVK energy should be symmetric under axis swap, {e1} vs {e2}"
        );
    }
}
