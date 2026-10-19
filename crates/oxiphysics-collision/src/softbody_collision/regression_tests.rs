// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the box query, the contact response and the exact
//! edge-edge CCD of [`super`].

use super::*;

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn vec_close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    len3(sub3(a, b)) <= tol
}

// ---------------------------------------------------------------------------
// PrimitiveVsDeformable::query_box
// ---------------------------------------------------------------------------

const BOX_MIN: [f64; 3] = [-1.0, -1.0, -1.0];
const BOX_MAX: [f64; 3] = [1.0, 1.0, 1.0];

#[test]
fn query_box_inside_vertex_near_positive_x_face() {
    let contacts = PrimitiveVsDeformable::query_box(&[[0.9, 0.0, 0.0]], BOX_MIN, BOX_MAX);
    assert_eq!(contacts.len(), 1);
    let c = &contacts[0];
    assert_eq!(c.vertex_index, 0);
    assert!(approx(c.depth, 0.1, 1e-12), "depth={}", c.depth);
    assert_eq!(c.normal, [1.0, 0.0, 0.0]);
    assert!(vec_close(c.point, [1.0, 0.0, 0.0], 1e-12), "{:?}", c.point);
}

#[test]
fn query_box_inside_vertex_near_negative_y_face() {
    let contacts = PrimitiveVsDeformable::query_box(&[[0.0, -0.95, 0.0]], BOX_MIN, BOX_MAX);
    assert_eq!(contacts.len(), 1);
    let c = &contacts[0];
    assert!(approx(c.depth, 0.05, 1e-12), "depth={}", c.depth);
    assert_eq!(c.normal, [0.0, -1.0, 0.0]);
    assert!(vec_close(c.point, [0.0, -1.0, 0.0], 1e-12), "{:?}", c.point);
}

#[test]
fn query_box_picks_nearest_face_on_every_axis() {
    let min = [-2.0, -3.0, -4.0];
    let max = [2.0, 3.0, 4.0];
    let cases: [([f64; 3], [f64; 3], f64); 6] = [
        ([-1.9, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.1),
        ([1.7, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3),
        ([0.0, -2.5, 0.0], [0.0, -1.0, 0.0], 0.5),
        ([0.0, 2.9, 0.0], [0.0, 1.0, 0.0], 0.1),
        ([0.0, 0.0, -3.2], [0.0, 0.0, -1.0], 0.8),
        ([0.0, 0.0, 3.6], [0.0, 0.0, 1.0], 0.4),
    ];
    for (vertex, normal, depth) in cases {
        let contacts = PrimitiveVsDeformable::query_box(&[vertex], min, max);
        assert_eq!(contacts.len(), 1, "vertex {vertex:?}");
        let c = &contacts[0];
        assert_eq!(c.normal, normal, "vertex {vertex:?}");
        assert!(approx(c.depth, depth, 1e-12), "vertex {vertex:?}");
        // The contact point lies on the face and the vertex is `depth` behind it.
        let moved = add3(c.point, scale3(c.normal, -c.depth));
        assert!(vec_close(moved, vertex, 1e-12), "vertex {vertex:?}");
    }
}

#[test]
fn query_box_depth_resolves_the_contact_exactly() {
    // Pushing the vertex along the normal by the depth puts it on the surface.
    let vertex = [0.3, -0.2, 0.85];
    let contacts = PrimitiveVsDeformable::query_box(&[vertex], BOX_MIN, BOX_MAX);
    assert_eq!(contacts.len(), 1);
    let c = &contacts[0];
    assert_eq!(c.normal, [0.0, 0.0, 1.0]);
    let resolved = add3(vertex, scale3(c.normal, c.depth));
    assert!(vec_close(resolved, [0.3, -0.2, 1.0], 1e-12), "{resolved:?}");
}

#[test]
fn query_box_reports_vertex_indices_and_skips_outside_and_surface_vertices() {
    let vertices = [
        [5.0, 0.0, 0.0],   // outside
        [0.0, 0.0, 0.5],   // inside
        [1.0, 0.0, 0.0],   // exactly on the surface
        [-0.5, 0.25, 0.0], // inside
        [0.0, f64::NAN, 0.0],
    ];
    let contacts = PrimitiveVsDeformable::query_box(&vertices, BOX_MIN, BOX_MAX);
    let indices: Vec<usize> = contacts.iter().map(|c| c.vertex_index).collect();
    assert_eq!(indices, vec![1, 3]);
    assert!(
        contacts
            .iter()
            .all(|c| c.depth > 0.0 && c.depth.is_finite())
    );
}

#[test]
fn query_box_tie_between_faces_is_resolved_deterministically() {
    // The centre of a cube is equidistant from all six faces: -X wins.
    let contacts = PrimitiveVsDeformable::query_box(&[[0.0, 0.0, 0.0]], BOX_MIN, BOX_MAX);
    assert_eq!(contacts.len(), 1);
    assert_eq!(contacts[0].normal, [-1.0, 0.0, 0.0]);
    assert!(approx(contacts[0].depth, 1.0, 1e-12));
}

#[test]
fn query_box_inverted_or_nan_box_is_empty() {
    let vertices = [[0.0, 0.0, 0.0]];
    let inverted = PrimitiveVsDeformable::query_box(&vertices, [1.0, -1.0, -1.0], [-1.0, 1.0, 1.0]);
    assert!(inverted.is_empty());
    let nan_box =
        PrimitiveVsDeformable::query_box(&vertices, [f64::NAN, -1.0, -1.0], [1.0, 1.0, 1.0]);
    assert!(nan_box.is_empty());
}

#[test]
fn query_box_with_thickness_reports_the_shell_outside_the_box() {
    // 0.03 outside the +X face with a 0.05 shell.
    let contacts = PrimitiveVsDeformable::query_box_with_thickness(
        &[[1.03, 0.0, 0.0], [1.06, 0.0, 0.0]],
        BOX_MIN,
        BOX_MAX,
        0.05,
    );
    assert_eq!(contacts.len(), 1);
    let c = &contacts[0];
    assert_eq!(c.vertex_index, 0);
    assert!(approx(c.depth, 0.02, 1e-12), "depth={}", c.depth);
    assert!(vec_close(c.normal, [1.0, 0.0, 0.0], 1e-12));
    assert!(vec_close(c.point, [1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn query_box_with_thickness_uses_the_corner_direction_outside_a_corner() {
    let contacts = PrimitiveVsDeformable::query_box_with_thickness(
        &[[1.03, 1.03, 0.0]],
        BOX_MIN,
        BOX_MAX,
        0.05,
    );
    assert_eq!(contacts.len(), 1);
    let c = &contacts[0];
    let dist = (2.0_f64 * 0.03 * 0.03).sqrt();
    assert!(approx(c.depth, 0.05 - dist, 1e-12), "depth={}", c.depth);
    let s = std::f64::consts::FRAC_1_SQRT_2;
    assert!(vec_close(c.normal, [s, s, 0.0], 1e-12), "{:?}", c.normal);
    assert!(vec_close(c.point, [1.0, 1.0, 0.0], 1e-12), "{:?}", c.point);
}

#[test]
fn query_box_with_thickness_adds_the_shell_to_the_inside_depth() {
    let contacts =
        PrimitiveVsDeformable::query_box_with_thickness(&[[0.9, 0.0, 0.0]], BOX_MIN, BOX_MAX, 0.05);
    assert_eq!(contacts.len(), 1);
    assert!(approx(contacts[0].depth, 0.15, 1e-12));
    assert_eq!(contacts[0].normal, [1.0, 0.0, 0.0]);
    assert!(vec_close(contacts[0].point, [1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn query_box_with_zero_thickness_matches_query_box() {
    let vertices = [[0.2, 0.1, -0.4], [3.0, 0.0, 0.0], [0.0, 0.0, 0.99]];
    let a = PrimitiveVsDeformable::query_box(&vertices, BOX_MIN, BOX_MAX);
    let b = PrimitiveVsDeformable::query_box_with_thickness(&vertices, BOX_MIN, BOX_MAX, -1.0);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.vertex_index, y.vertex_index);
        assert_eq!(x.depth, y.depth);
        assert_eq!(x.normal, y.normal);
        assert_eq!(x.point, y.point);
    }
}

// ---------------------------------------------------------------------------
// ContactResponse
// ---------------------------------------------------------------------------

fn deform_contact(vertex_index: usize, depth: f64) -> DeformContact {
    DeformContact {
        vertex_index,
        tri_index: 0,
        depth,
        normal: [0.0, 1.0, 0.0],
        point_on_rigid: [0.0; 3],
    }
}

fn edge_result(distance: f64) -> EdgeEdgeResult {
    EdgeEdgeResult {
        edge_a: 0,
        edge_b: 0,
        distance,
        s: 0.5,
        t: 0.5,
        normal: [0.0, 1.0, 0.0],
    }
}

#[test]
fn contact_response_new_has_zero_thickness_and_builder_sets_it() {
    let cr = ContactResponse::new(0.3, 0.6);
    assert_eq!(cr.thickness, 0.0);
    let cr = cr.with_thickness(0.05);
    assert_eq!(cr.thickness, 0.05);
    assert_eq!(cr.restitution, 0.3);
    assert_eq!(cr.friction, 0.6);
}

#[test]
fn deform_rigid_separating_vertex_gets_no_impulse() {
    let cr = ContactResponse::new(0.5, 0.5);
    let velocities = [[0.0, 1.0, 0.0]];
    let masses = [1.0];
    // Not penetrating and separating: nothing to do at all.
    let none = cr.solve_deform_rigid(&[deform_contact(0, 0.0)], &velocities, &masses);
    assert!(none.is_empty());
    // Penetrating but separating: position correction only, never a pull-back impulse.
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.02)], &velocities, &masses);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].impulse, [0.0; 3]);
    assert!(vec_close(out[0].correction, [0.0, 0.02, 0.0], 1e-15));
}

#[test]
fn deform_rigid_resting_vertex_gets_no_normal_impulse() {
    let cr = ContactResponse::new(1.0, 0.5);
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.01)], &[[2.0, 0.0, 0.0]], &[1.0]);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].impulse, [0.0; 3]);
}

#[test]
fn deform_rigid_approaching_vertex_gets_the_restitution_impulse() {
    let cr = ContactResponse::new(0.5, 0.0);
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.02)], &[[0.0, -1.0, 0.0]], &[2.0]);
    assert_eq!(out.len(), 1);
    // j_n = (1 + e) * |v_n| * m = 1.5 * 1 * 2.
    assert!(vec_close(out[0].impulse, [0.0, 3.0, 0.0], 1e-12));
    assert!(vec_close(out[0].correction, [0.0, 0.02, 0.0], 1e-15));
}

#[test]
fn deform_rigid_friction_never_reverses_the_tangential_velocity() {
    // Sliding at 0.1 m/s while hitting the surface at 2 m/s: Coulomb friction
    // mu * j_n = 2 would far exceed the momentum that merely stops the sliding.
    let velocity = [0.1, -2.0, 0.0];
    let mass = 1.0;
    let cr = ContactResponse::new(0.0, 1.0);
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.0)], &[velocity], &[mass]);
    assert_eq!(out.len(), 1);
    let new_vt = velocity[0] + out[0].impulse[0] / mass;
    assert!(
        new_vt.abs() < 1e-12,
        "tangential velocity after friction: {new_vt}"
    );
    assert!(approx(out[0].impulse[1], 2.0, 1e-12));

    // Weak friction only slows the sliding down: mu * j_n = 0.02 < m |v_t|.
    let cr = ContactResponse::new(0.0, 0.01);
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.0)], &[velocity], &[mass]);
    let new_vt = velocity[0] + out[0].impulse[0] / mass;
    assert!(approx(new_vt, 0.08, 1e-12), "tangential velocity: {new_vt}");
}

#[test]
fn deform_rigid_friction_clamp_scales_with_mass() {
    let velocity = [0.5, -4.0, 0.0];
    let mass = 3.0;
    let cr = ContactResponse::new(0.0, 5.0);
    let out = cr.solve_deform_rigid(&[deform_contact(0, 0.0)], &[velocity], &[mass]);
    assert_eq!(out.len(), 1);
    // The tangential impulse is exactly the momentum that stops the sliding.
    assert!(approx(out[0].impulse[0], -mass * 0.5, 1e-12));
}

#[test]
fn deform_rigid_skips_invalid_indices_and_masses() {
    let cr = ContactResponse::new(0.0, 0.0);
    let velocities = [[0.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, -1.0, 0.0]];
    let masses = [1.0, 0.0, f64::NAN];
    let contacts = [
        deform_contact(0, 0.0),
        deform_contact(1, 0.0), // zero mass
        deform_contact(2, 0.0), // NaN mass
        deform_contact(9, 0.0), // out of range
    ];
    let out = cr.solve_deform_rigid(&contacts, &velocities, &masses);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].vertex_index, 0);
    // A mass slice shorter than the velocity slice must not panic either.
    let short = cr.solve_deform_rigid(&[deform_contact(2, 0.0)], &velocities, &masses[..2]);
    assert!(short.is_empty());
}

#[test]
fn edge_edge_depth_comes_from_the_thickness_not_the_friction() {
    let velocities = [[0.0, -1.0, 0.0], [0.0, -1.0, 0.0]];
    let masses = [1.0, 1.0];
    let edges = [(0, 1)];
    let results = [edge_result(0.02)];
    for friction in [0.0, 0.2, 0.9] {
        let cr = ContactResponse::new(0.0, friction).with_thickness(0.05);
        let out = cr.solve_edge_edge(&results, &velocities, &edges, &masses);
        assert_eq!(out.len(), 2, "friction {friction}");
        for response in &out {
            assert!(
                vec_close(response.correction, [0.0, 0.03, 0.0], 1e-12),
                "friction {friction}: {:?}",
                response.correction
            );
            assert!(approx(response.impulse[1], 1.0, 1e-12));
        }
    }
}

#[test]
fn edge_edge_default_thickness_produces_no_position_correction() {
    let cr = ContactResponse::new(0.0, 0.9);
    let out = cr.solve_edge_edge(
        &[edge_result(0.02)],
        &[[0.0, -1.0, 0.0], [0.0, -1.0, 0.0]],
        &[(0, 1)],
        &[1.0, 1.0],
    );
    assert_eq!(out.len(), 2);
    for response in &out {
        assert_eq!(response.correction, [0.0; 3]);
    }
}

#[test]
fn edge_edge_separating_vertices_get_no_impulse() {
    let cr = ContactResponse::new(0.5, 0.5).with_thickness(0.05);
    let velocities = [[0.0, 1.0, 0.0], [0.0, 3.0, 0.0]];
    let out = cr.solve_edge_edge(&[edge_result(0.02)], &velocities, &[(0, 1)], &[1.0, 1.0]);
    assert_eq!(out.len(), 2);
    for response in &out {
        assert_eq!(response.impulse, [0.0; 3]);
        assert!(vec_close(response.correction, [0.0, 0.03, 0.0], 1e-12));
    }
    // Outside the shell and separating: nothing to report.
    let far = cr.solve_edge_edge(&[edge_result(0.2)], &velocities, &[(0, 1)], &[1.0, 1.0]);
    assert!(far.is_empty());
}

#[test]
fn edge_edge_friction_is_clamped_like_the_vertex_contact() {
    let cr = ContactResponse::new(0.0, 1.0).with_thickness(0.05);
    let velocities = [[0.1, -2.0, 0.0], [0.1, -2.0, 0.0]];
    let out = cr.solve_edge_edge(&[edge_result(0.0)], &velocities, &[(0, 1)], &[1.0, 1.0]);
    assert_eq!(out.len(), 2);
    for response in &out {
        let new_vt = 0.1 + response.impulse[0];
        assert!(new_vt.abs() < 1e-12, "tangential velocity: {new_vt}");
    }
}

#[test]
fn edge_edge_skips_missing_edges_vertices_and_masses() {
    let cr = ContactResponse::new(0.0, 0.0).with_thickness(0.05);
    let velocities = [[0.0, -1.0, 0.0], [0.0, -1.0, 0.0]];
    let masses = [1.0, 1.0];
    let mut missing_edge = edge_result(0.0);
    missing_edge.edge_a = 4;
    assert!(
        cr.solve_edge_edge(&[missing_edge], &velocities, &[(0, 1)], &masses)
            .is_empty()
    );
    // Vertex 7 does not exist: only vertex 0 responds.
    let out = cr.solve_edge_edge(&[edge_result(0.0)], &velocities, &[(0, 7)], &masses);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].vertex_index, 0);
    // Zero mass is skipped.
    let out = cr.solve_edge_edge(&[edge_result(0.0)], &velocities, &[(0, 1)], &[1.0, 0.0]);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].vertex_index, 0);
}

// ---------------------------------------------------------------------------
// CcdDeformable::edge_edge_ccd
// ---------------------------------------------------------------------------

/// Static edge A along +X from the origin, and an edge B along Y at `x` that
/// sweeps down in Z from `z0` to `z1` during the step.
struct SweepScene {
    edges_a: Vec<(usize, usize)>,
    verts_a: Vec<[f64; 3]>,
    edges_b: Vec<(usize, usize)>,
    verts0_b: Vec<[f64; 3]>,
    verts1_b: Vec<[f64; 3]>,
}

fn sweep_scene(x: f64, z0: f64, z1: f64) -> SweepScene {
    SweepScene {
        edges_a: vec![(0, 1)],
        verts_a: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        edges_b: vec![(0, 1)],
        verts0_b: vec![[x, -0.5, z0], [x, 0.5, z0]],
        verts1_b: vec![[x, -0.5, z1], [x, 0.5, z1]],
    }
}

fn run_ccd(ccd: &CcdDeformable, scene: &SweepScene) -> Vec<EdgeCcdResult> {
    ccd.edge_edge_ccd(
        &scene.edges_a,
        &scene.verts_a,
        &scene.verts_a,
        &scene.edges_b,
        &scene.verts0_b,
        &scene.verts1_b,
    )
}

#[test]
fn edge_edge_ccd_detects_x_crossing_near_the_edge_ends() {
    // B crosses A at x = 0.98, 2 cm from A's end point. The midpoints of the
    // two edges stay at least 0.48 apart, so the old midpoint approximation
    // never reported this crossing.
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(0.98, 1.0, -1.0);
    let results = run_ccd(&ccd, &scene);
    assert_eq!(results.len(), 1, "{results:?}");
    let r = &results[0];
    assert_eq!((r.edge_a, r.edge_b), (0, 0));
    assert!(approx(r.toi, 0.5, 1e-6), "toi={}", r.toi);
    // B starts above A, so the normal from B to A points down.
    assert!(
        vec_close(r.normal, [0.0, 0.0, -1.0], 1e-9),
        "{:?}",
        r.normal
    );
}

#[test]
fn edge_edge_ccd_normal_is_oriented_by_the_side_edge_a_started_on() {
    // Same crossing but B starts below A and rises: the normal from B to A flips.
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(0.98, -1.0, 1.0);
    let results = run_ccd(&ccd, &scene);
    assert_eq!(results.len(), 1);
    assert!(approx(results[0].toi, 0.5, 1e-6));
    assert!(vec_close(results[0].normal, [0.0, 0.0, 1.0], 1e-9));
}

#[test]
fn edge_edge_ccd_misses_a_crossing_beyond_the_end_of_the_edge() {
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(1.5, 1.0, -1.0);
    assert!(run_ccd(&ccd, &scene).is_empty());
}

#[test]
fn edge_edge_ccd_reports_the_first_shell_entry_of_a_graze() {
    // B never reaches A (it stops 0.03 above), but it enters the 0.05 shell at
    // z(t) = 1 - 0.97 t = 0.05, i.e. t = 0.95 / 0.97.
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(0.5, 1.0, 0.03);
    let results = run_ccd(&ccd, &scene);
    assert_eq!(results.len(), 1);
    assert!(
        approx(results[0].toi, 0.95 / 0.97, 1e-6),
        "toi={}",
        results[0].toi
    );
    assert!(vec_close(results[0].normal, [0.0, 0.0, -1.0], 1e-9));
}

#[test]
fn edge_edge_ccd_reports_initial_contact_at_time_zero() {
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(0.5, 0.03, 0.03);
    let results = run_ccd(&ccd, &scene);
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].toi, 0.0);
    assert!(vec_close(results[0].normal, [0.0, 0.0, -1.0], 1e-9));
}

#[test]
fn edge_edge_ccd_reports_indices_of_the_colliding_pair_only() {
    // Polyline A along X with two edges, two B edges: only A1 x B0 collide.
    let ccd = CcdDeformable::new(0.05, 0.05);
    let verts_a = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let edges_a = vec![(0, 1), (1, 2)];
    let verts0_b = vec![
        [1.98, -0.5, 1.0],
        [1.98, 0.5, 1.0],
        [50.0, -0.5, 1.0],
        [50.0, 0.5, 1.0],
    ];
    let verts1_b = vec![
        [1.98, -0.5, -1.0],
        [1.98, 0.5, -1.0],
        [50.0, -0.5, -1.0],
        [50.0, 0.5, -1.0],
    ];
    let edges_b = vec![(0, 1), (2, 3)];
    let results = ccd.edge_edge_ccd(&edges_a, &verts_a, &verts_a, &edges_b, &verts0_b, &verts1_b);
    assert_eq!(results.len(), 1, "{results:?}");
    assert_eq!((results[0].edge_a, results[0].edge_b), (1, 0));
    assert!(approx(results[0].toi, 0.5, 1e-6));
}

#[test]
fn edge_edge_ccd_skips_edges_with_out_of_range_vertices() {
    let ccd = CcdDeformable::new(0.05, 0.05);
    let mut scene = sweep_scene(0.98, 1.0, -1.0);
    // A has an invalid edge first; B has one invalid edge too.
    scene.edges_a.insert(0, (0, 9));
    scene.edges_b.push((1, 4));
    let results = run_ccd(&ccd, &scene);
    assert_eq!(results.len(), 1);
    // The valid A edge is now at index 1.
    assert_eq!((results[0].edge_a, results[0].edge_b), (1, 0));

    // End positions missing for a vertex: that edge is skipped, not a panic.
    let short = ccd.edge_edge_ccd(
        &scene.edges_a,
        &scene.verts_a,
        &scene.verts_a,
        &scene.edges_b,
        &scene.verts0_b,
        &scene.verts1_b[..1],
    );
    assert!(short.is_empty());
}

#[test]
fn edge_edge_ccd_handles_parallel_edges_approaching_each_other() {
    // Two parallel edges along X; B descends onto A and ends 0.02 above it.
    let ccd = CcdDeformable::new(0.05, 0.05);
    let verts_a = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let verts0_b = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]];
    let verts1_b = vec![[0.0, 0.0, 0.02], [1.0, 0.0, 0.02]];
    let results = ccd.edge_edge_ccd(
        &[(0, 1)],
        &verts_a,
        &verts_a,
        &[(0, 1)],
        &verts0_b,
        &verts1_b,
    );
    assert_eq!(results.len(), 1);
    // z(t) = 1 - 0.98 t reaches the 0.05 shell at t = 0.95 / 0.98.
    assert!(
        approx(results[0].toi, 0.95 / 0.98, 1e-6),
        "{}",
        results[0].toi
    );
    assert!(vec_close(results[0].normal, [0.0, 0.0, -1.0], 1e-9));
}

#[test]
fn edge_edge_ccd_returns_nothing_for_far_apart_edges() {
    let ccd = CcdDeformable::new(0.05, 0.05);
    let scene = sweep_scene(0.5, 10.0, 8.0);
    assert!(run_ccd(&ccd, &scene).is_empty());
}

// ---------------------------------------------------------------------------
// CcdDeformable::edge_edge_ccd with a contact radius of zero
// ---------------------------------------------------------------------------

/// End points of edge A at the crossing time of the zero-radius scenes.
const CROSSING_A0: [f64; 3] = [0.1, -0.2, 0.3];
const CROSSING_A1: [f64; 3] = [1.7, 0.9, -0.4];
/// Direction of edge B at the crossing time of the zero-radius scenes.
const CROSSING_DIR_B: [f64; 3] = [0.3, -0.5, 0.9];
/// Velocity of edge B in the translating scenes: mostly along the common
/// normal of the two edges.
const TRANSLATION: [f64; 3] = [0.84, -1.95, -0.68];

/// Four points `[a0, a1, b0, b1]` that move linearly over one step: edge A joins
/// the first two, edge B the last two.
struct CrossingPair {
    start: [[f64; 3]; 4],
    end: [[f64; 3]; 4],
}

/// Build four points with the given constant velocities such that, at the
/// normalised time `t_cross`, edge B crosses edge A at parameter `s_cross` of A
/// (and at parameter `0.6 / 1.4` of B), all in generic, non-round positions.
fn crossing_pair(s_cross: f64, t_cross: f64, velocities: [[f64; 3]; 4]) -> CrossingPair {
    let crossing_point = add3(CROSSING_A0, scale3(sub3(CROSSING_A1, CROSSING_A0), s_cross));
    let at_cross = [
        CROSSING_A0,
        CROSSING_A1,
        sub3(crossing_point, scale3(CROSSING_DIR_B, 0.6)),
        add3(crossing_point, scale3(CROSSING_DIR_B, 0.8)),
    ];
    let mut start = at_cross;
    let mut end = at_cross;
    for (k, velocity) in velocities.iter().enumerate() {
        start[k] = sub3(at_cross[k], scale3(*velocity, t_cross));
        end[k] = add3(at_cross[k], scale3(*velocity, 1.0 - t_cross));
    }
    CrossingPair { start, end }
}

/// Velocities of the four points when edge A rests and edge B translates.
fn translating_b(velocity: [f64; 3]) -> [[f64; 3]; 4] {
    [[0.0; 3], [0.0; 3], velocity, velocity]
}

/// The pair with every position multiplied by `factor` (the crossing time and
/// the crossing geometry are unchanged).
fn scaled(pair: &CrossingPair, factor: f64) -> CrossingPair {
    CrossingPair {
        start: pair.start.map(|p| scale3(p, factor)),
        end: pair.end.map(|p| scale3(p, factor)),
    }
}

fn run_crossing_pair(pair: &CrossingPair, radius: f64) -> Vec<EdgeCcdResult> {
    CcdDeformable::new(radius, radius).edge_edge_ccd(
        &[(0, 1)],
        &pair.start[..2],
        &pair.end[..2],
        &[(0, 1)],
        &pair.start[2..],
        &pair.end[2..],
    )
}

/// Unit common normal of the two edges at the crossing.
fn crossing_normal() -> [f64; 3] {
    normalize3(cross3(sub3(CROSSING_A1, CROSSING_A0), CROSSING_DIR_B))
}

#[test]
fn edge_edge_ccd_with_zero_radius_detects_a_generic_exact_crossing() {
    // B passes through A at t = 0.4. The closest-point distance of the edges at
    // the rounded root time is about 1e-16, not 0, so a contact radius of
    // exactly zero used to reject the root and report nothing.
    let pair = crossing_pair(0.37, 0.4, translating_b(TRANSLATION));
    let results = run_crossing_pair(&pair, 0.0);
    assert_eq!(results.len(), 1, "{results:?}");
    let r = &results[0];
    assert_eq!((r.edge_a, r.edge_b), (0, 0));
    assert!(approx(r.toi, 0.4, 1e-9), "toi={}", r.toi);
    // B moves along +n through A, so it starts on the -n side of A and the
    // normal from B towards A is +n.
    assert!(approx(len3(r.normal), 1.0, 1e-12), "{:?}", r.normal);
    assert!(
        dot3(r.normal, crossing_normal()) > 1.0 - 1e-9,
        "{:?} vs {:?}",
        r.normal,
        crossing_normal()
    );
}

#[test]
fn edge_edge_ccd_with_zero_radius_orients_the_normal_by_the_starting_side() {
    // The reversed motion crosses at the same time but B now starts on the +n
    // side of A, so the normal flips.
    let reversed = scale3(TRANSLATION, -1.0);
    let pair = crossing_pair(0.37, 0.4, translating_b(reversed));
    let results = run_crossing_pair(&pair, 0.0);
    assert_eq!(results.len(), 1, "{results:?}");
    assert!(approx(results[0].toi, 0.4, 1e-9), "toi={}", results[0].toi);
    assert!(
        dot3(results[0].normal, crossing_normal()) < -(1.0 - 1e-9),
        "{:?}",
        results[0].normal
    );
}

#[test]
fn edge_edge_ccd_with_zero_radius_solves_a_true_cubic() {
    // All four end points move with different velocities, so the coplanarity
    // condition is a genuine cubic in time rather than a linear equation.
    let velocities = [
        [0.30, -0.20, 0.50],
        [-0.40, 0.60, 0.10],
        [0.70, -1.10, -0.90],
        [-0.50, 0.90, 1.40],
    ];
    let pair = crossing_pair(0.37, 0.55, velocities);
    let results = run_crossing_pair(&pair, 0.0);
    assert_eq!(results.len(), 1, "{results:?}");
    assert!(approx(results[0].toi, 0.55, 1e-9), "toi={}", results[0].toi);
    assert!(
        dot3(results[0].normal, crossing_normal()) > 1.0 - 1e-9,
        "{:?}",
        results[0].normal
    );
}

#[test]
fn edge_edge_ccd_with_zero_radius_misses_a_crossing_beyond_the_end_of_the_edge() {
    // The line of B crosses the line of A 0.79 beyond A's end point, so the two
    // segments never touch (they stay 0.78 apart at their closest).
    let pair = crossing_pair(1.4, 0.4, translating_b(TRANSLATION));
    assert!(run_crossing_pair(&pair, 0.0).is_empty());
}

#[test]
fn edge_edge_ccd_with_zero_radius_does_not_depend_on_the_coordinate_scale() {
    // The rounding floor of the radius grows with the magnitude of the
    // coordinates, so the same crossing is found at every scale.
    let base = crossing_pair(0.37, 0.4, translating_b(TRANSLATION));
    for factor in [0.1, 1.0, 1e3] {
        let results = run_crossing_pair(&scaled(&base, factor), 0.0);
        assert_eq!(results.len(), 1, "factor {factor}: {results:?}");
        assert!(
            approx(results[0].toi, 0.4, 1e-9),
            "factor {factor}: toi={}",
            results[0].toi
        );
        assert!(
            dot3(results[0].normal, crossing_normal()) > 1.0 - 1e-9,
            "factor {factor}: {:?}",
            results[0].normal
        );
    }
}

#[test]
fn edge_edge_ccd_reports_nothing_for_a_negative_or_nan_radius() {
    let pair = crossing_pair(0.37, 0.4, translating_b(TRANSLATION));
    assert_eq!(run_crossing_pair(&pair, 0.0).len(), 1);
    assert!(run_crossing_pair(&pair, -0.01).is_empty());
    assert!(run_crossing_pair(&pair, f64::NAN).is_empty());
}

// ---------------------------------------------------------------------------
// CcdDeformable::edge_edge_ccd does not depend on the units of the mesh
// ---------------------------------------------------------------------------

/// The scene with every position multiplied by `factor`.
fn scaled_scene(scene: &SweepScene, factor: f64) -> SweepScene {
    let scale_all = |points: &[[f64; 3]]| -> Vec<[f64; 3]> {
        points.iter().map(|&p| scale3(p, factor)).collect()
    };
    SweepScene {
        edges_a: scene.edges_a.clone(),
        verts_a: scale_all(&scene.verts_a),
        edges_b: scene.edges_b.clone(),
        verts0_b: scale_all(&scene.verts0_b),
        verts1_b: scale_all(&scene.verts1_b),
    }
}

#[test]
fn edge_edge_ccd_is_invariant_under_uniform_scaling_of_mesh_and_radius() {
    // A crossing near the end of A, a graze of the radius shell and an initial
    // contact, each with the expected time of impact.
    let scenes = [
        (sweep_scene(0.98, 1.0, -1.0), 0.5),
        (sweep_scene(0.5, 1.0, 0.03), 0.95 / 0.97),
        (sweep_scene(0.5, 0.03, 0.03), 0.0),
    ];
    for (scene, expected_toi) in &scenes {
        for factor in [1e-4, 1e-2, 1.0, 1e2] {
            let ccd = CcdDeformable::new(0.05 * factor, 0.05 * factor);
            let results = run_ccd(&ccd, &scaled_scene(scene, factor));
            assert_eq!(results.len(), 1, "factor {factor}: {results:?}");
            assert!(
                approx(results[0].toi, *expected_toi, 1e-6),
                "factor {factor}: toi={} expected {expected_toi}",
                results[0].toi
            );
            assert!(
                vec_close(results[0].normal, [0.0, 0.0, -1.0], 1e-9),
                "factor {factor}: {:?}",
                results[0].normal
            );
        }
    }
}

#[test]
fn edge_edge_ccd_detects_a_shallow_crossing_of_short_edges() {
    // Two edges of about 1 cm cross at an angle of one degree: B slides down
    // through A in the plane z = 0.005 - 0.01 t and meets it at t = 0.5. The
    // squared lengths of the edges are so small that their parallelism, judged
    // with an absolute tolerance, would make the closest points meaningless.
    let angle = 1.0_f64.to_radians();
    let half = scale3([angle.cos(), angle.sin(), 0.0], 0.004);
    let verts_a = [[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]];
    let centre0 = [0.005, 0.0, 0.005];
    let centre1 = [0.005, 0.0, -0.005];
    let verts0_b = [sub3(centre0, half), add3(centre0, half)];
    let verts1_b = [sub3(centre1, half), add3(centre1, half)];
    // On the raw coordinates the shared routine would place the closest points
    // about 7e-5 apart at the crossing, i.e. outside both of these radii.
    for radius in [0.0, 1e-5] {
        let results = CcdDeformable::new(radius, radius).edge_edge_ccd(
            &[(0, 1)],
            &verts_a,
            &verts_a,
            &[(0, 1)],
            &verts0_b,
            &verts1_b,
        );
        assert_eq!(results.len(), 1, "radius {radius}: {results:?}");
        assert!(
            approx(results[0].toi, 0.5, 1e-9),
            "radius {radius}: toi={}",
            results[0].toi
        );
        assert!(
            vec_close(results[0].normal, [0.0, 0.0, -1.0], 1e-9),
            "radius {radius}: {:?}",
            results[0].normal
        );
    }
}
