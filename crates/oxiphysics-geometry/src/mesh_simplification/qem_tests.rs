//! Regression tests for QEM simplification.
//!
//! The earlier simplifier kept stale entries in its heap, so a position and cost
//! computed for one edge were applied to another; scanned every triangle on every
//! collapse; accepted collapses that pinch or fold the surface; and read only
//! three of the six options of [`QemConfig`] (`check_normal_deviation` always
//! said yes).

use super::*;
use crate::test_meshes::{flat_grid, icosphere};
use std::collections::HashMap;

type Mesh = (Vec<[f64; 3]>, Vec<[usize; 3]>);

fn config(target: usize) -> QemConfig {
    QemConfig {
        target_triangles: target,
        max_error: 0.0,
        preserve_boundary: true,
        max_normal_deviation: 0.0,
        preserve_texture_seams: false,
        boundary_weight: 100.0,
    }
}

/// Number of triangles that use each undirected edge.
fn edge_usage(triangles: &[[usize; 3]]) -> HashMap<(usize, usize), usize> {
    let mut usage = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            *usage.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    usage
}

/// Asserts that the mesh is a closed, valid 2-manifold of genus 0 without
/// unreferenced vertices.
fn assert_closed_sphere_topology(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) {
    for tri in triangles {
        assert!(
            tri.iter().all(|&v| v < vertices.len()),
            "index out of range"
        );
        assert!(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
    }
    let usage = edge_usage(triangles);
    for (edge, count) in &usage {
        assert_eq!(*count, 2, "edge {edge:?} is used by {count} triangles");
    }
    let mut used = vec![false; vertices.len()];
    for tri in triangles {
        for &v in tri {
            used[v] = true;
        }
    }
    assert!(used.iter().all(|&u| u), "unreferenced vertex");
    assert_eq!(
        euler_characteristic(vertices.len(), triangles),
        2,
        "not a topological sphere"
    );
    let validation = validate_topology(vertices, triangles);
    assert!(validation.is_valid, "{validation:?}");
    assert_eq!(validation.isolated_vertices, 0);
}

/// Points on the surface: a barycentric lattice of `n + 1` points per triangle edge.
fn surface_samples(vertices: &[[f64; 3]], triangles: &[[usize; 3]], n: usize) -> Vec<[f64; 3]> {
    let mut samples = Vec::new();
    for tri in triangles {
        let (a, b, c) = (vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
        for i in 0..=n {
            for j in 0..=(n - i) {
                let (u, v) = (i as f64 / n as f64, j as f64 / n as f64);
                let w = 1.0 - u - v;
                samples.push([
                    w * a[0] + u * b[0] + v * c[0],
                    w * a[1] + u * b[1] + v * c[1],
                    w * a[2] + u * b[2] + v * c[2],
                ]);
            }
        }
    }
    samples
}

fn one_sided_hausdorff(from: &Mesh, to: &Mesh) -> f64 {
    surface_samples(&from.0, &from.1, 4)
        .iter()
        .map(|&p| {
            to.1.iter()
                .map(|tri| point_to_triangle_dist(p, to.0[tri[0]], to.0[tri[1]], to.0[tri[2]]))
                .fold(f64::MAX, f64::min)
        })
        .fold(0.0, f64::max)
}

/// Symmetric Hausdorff distance between the surfaces, estimated on dense samples.
fn hausdorff(a: &Mesh, b: &Mesh) -> f64 {
    one_sided_hausdorff(a, b).max(one_sided_hausdorff(b, a))
}

fn radius(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

fn sorted_positions(points: &[[f64; 3]]) -> Vec<[u64; 3]> {
    let mut keys: Vec<[u64; 3]> = points
        .iter()
        .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
        .collect();
    keys.sort_unstable();
    keys
}

fn cube() -> Mesh {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let triangles = vec![
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [3, 6, 2],
        [3, 7, 6],
        [0, 7, 3],
        [0, 4, 7],
        [1, 2, 6],
        [1, 6, 5],
    ];
    (vertices, triangles)
}

// ---------------------------------------------------------------------------
// Closed surfaces
// ---------------------------------------------------------------------------

#[test]
fn icosphere_at_half_is_a_valid_manifold_close_to_the_sphere() {
    let original = icosphere(2);
    assert_eq!(original.1.len(), 320);
    let simplified = qem_simplify(&original.0, &original.1, &config(160));
    // Every collapse of a closed surface removes two triangles.
    assert_eq!(simplified.1.len(), 160);
    assert_closed_sphere_topology(&simplified.0, &simplified.1);
    // The surface stays close to the sphere: no spikes, no collapse of the shape.
    // (The minimiser of the quadric of chord planes can lie a little outside the
    // sphere, hence the upper bound above one.)
    for &p in &simplified.0 {
        let r = radius(p);
        assert!((0.9..=1.05).contains(&r), "vertex at radius {r}");
    }
    let distance = hausdorff(&original, &simplified);
    assert!(distance < 0.05, "Hausdorff distance {distance}");
    let volume = mesh_signed_volume(&simplified.0, &simplified.1);
    let reference = mesh_signed_volume(&original.0, &original.1);
    assert!(
        volume > 0.85 * reference && volume < 1.01 * reference,
        "volume {volume} of {reference}"
    );
}

#[test]
fn simplified_faces_still_point_outwards() {
    let original = icosphere(3);
    let (vertices, triangles) = qem_simplify(&original.0, &original.1, &config(320));
    assert_eq!(triangles.len(), 320);
    for tri in &triangles {
        let (a, b, c) = (vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
        let normal = cross3(sub3(b, a), sub3(c, a));
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        let outward = dot3(normal, centroid);
        assert!(outward > 0.0, "folded triangle {tri:?}");
    }
    assert_closed_sphere_topology(&vertices, &triangles);
}

#[test]
fn the_target_is_reached_down_to_a_coarse_mesh() {
    let original = icosphere(2);
    for target in [300, 250, 200, 120, 80, 60, 40] {
        let (vertices, triangles) = qem_simplify(&original.0, &original.1, &config(target));
        assert_eq!(triangles.len(), target, "target {target}");
        assert_closed_sphere_topology(&vertices, &triangles);
    }
    // An odd target is met from above in steps of two triangles.
    let (_, triangles) = qem_simplify(&original.0, &original.1, &config(101));
    assert_eq!(triangles.len(), 100);
}

#[test]
fn a_closed_surface_is_never_reduced_below_a_tetrahedron() {
    let original = icosphere(1);
    let (vertices, triangles) = qem_simplify(&original.0, &original.1, &config(1));
    assert!(triangles.len() >= 4, "{} triangles", triangles.len());
    assert_closed_sphere_topology(&vertices, &triangles);
    let (cube_vertices, cube_triangles) = cube();
    let mut free = config(1);
    free.preserve_boundary = false;
    let (vertices, triangles) = qem_simplify(&cube_vertices, &cube_triangles, &free);
    assert!(triangles.len() >= 4);
    assert_closed_sphere_topology(&vertices, &triangles);
}

#[test]
fn the_result_does_not_depend_on_anything_but_the_input() {
    let original = icosphere(2);
    let first = qem_simplify(&original.0, &original.1, &config(90));
    let second = qem_simplify(&original.0, &original.1, &config(90));
    assert_eq!(first, second);
}

#[test]
fn an_ellipsoid_is_simplified_along_its_curvature() {
    let (vertices, triangles) = icosphere(3);
    let stretched: Vec<[f64; 3]> = vertices
        .iter()
        .map(|p| [2.0 * p[0], p[1], 0.5 * p[2]])
        .collect();
    let simplified = qem_simplify(&stretched, &triangles, &config(320));
    assert_eq!(simplified.1.len(), 320);
    assert_closed_sphere_topology(&simplified.0, &simplified.1);
    let distance = hausdorff_symmetric(
        &stretched,
        &triangles,
        &simplified.0,
        &simplified.1,
        usize::MAX,
    );
    assert!(distance < 0.08, "Hausdorff distance {distance}");
}

#[test]
fn a_large_mesh_is_simplified_to_a_tenth() {
    let original = icosphere(4);
    assert_eq!(original.1.len(), 5120);
    let (vertices, triangles) = qem_simplify(&original.0, &original.1, &config(512));
    assert_eq!(triangles.len(), 512);
    assert_closed_sphere_topology(&vertices, &triangles);
    for &p in &vertices {
        let r = radius(p);
        assert!((0.95..=1.05).contains(&r), "vertex at radius {r}");
    }
    // Vertex-to-surface distances in both directions (the dense estimate of the
    // other tests is too slow for this size).
    let distance = hausdorff_symmetric(&original.0, &original.1, &vertices, &triangles, usize::MAX);
    assert!(distance < 0.05, "Hausdorff distance {distance}");
}

#[test]
fn the_default_configuration_halves_a_sphere() {
    let original = icosphere(2);
    let config = QemConfig::default_config(original.1.len());
    assert_eq!(config.target_triangles, 160);
    let simplified = qem_simplify(&original.0, &original.1, &config);
    assert_eq!(simplified.1.len(), 160);
    assert_closed_sphere_topology(&simplified.0, &simplified.1);
    let distance = hausdorff(&original, &simplified);
    assert!(distance < 0.05, "Hausdorff distance {distance}");
}

#[test]
fn simplify_to_ratio_halves_a_sphere() {
    let original = icosphere(2);
    let (vertices, triangles) = simplify_to_ratio(&original.0, &original.1, 0.5);
    assert_eq!(triangles.len(), 160);
    assert_closed_sphere_topology(&vertices, &triangles);
}

// ---------------------------------------------------------------------------
// max_error and max_normal_deviation
// ---------------------------------------------------------------------------

#[test]
fn the_error_threshold_limits_the_simplification() {
    let original = icosphere(2);
    let with_error = |max_error: f64| QemConfig {
        max_error,
        ..config(0)
    };
    let untouched = qem_simplify(&original.0, &original.1, &with_error(1.0e-12));
    assert_eq!(untouched.1.len(), 320);

    let mut previous = untouched.1.len();
    for max_error in [1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 10.0] {
        let (vertices, triangles) = qem_simplify(&original.0, &original.1, &with_error(max_error));
        assert!(
            triangles.len() <= previous,
            "max_error {max_error}: {} triangles after {previous}",
            triangles.len()
        );
        assert_closed_sphere_topology(&vertices, &triangles);
        previous = triangles.len();
    }
    assert!(previous < 320, "a large threshold must allow collapses");
}

#[test]
fn collapsing_everything_allowed_leaves_a_valid_small_mesh() {
    let original = icosphere(2);
    let (vertices, triangles) = qem_simplify(&original.0, &original.1, &config(0));
    assert!(
        triangles.len() >= 4 && triangles.len() < 40,
        "{}",
        triangles.len()
    );
    assert_closed_sphere_topology(&vertices, &triangles);
}

#[test]
fn the_normal_deviation_limit_blocks_collapses_that_tilt_faces() {
    let original = icosphere(2);
    let with_limit = |max_normal_deviation: f64| QemConfig {
        max_normal_deviation,
        ..config(160)
    };
    // Neighbouring faces of the sphere differ by far more than 0.01 rad, so each
    // collapse tilts some face by more than that.
    let blocked = qem_simplify(&original.0, &original.1, &with_limit(0.01));
    assert_eq!(blocked.1.len(), 320);
    assert_eq!(blocked.0.len(), original.0.len());
    // A generous limit allows the whole reduction.
    let free = qem_simplify(&original.0, &original.1, &with_limit(1.2));
    assert_eq!(free.1.len(), 160);
    assert_closed_sphere_topology(&free.0, &free.1);
    // A moderate limit lies in between.
    let moderate = qem_simplify(&original.0, &original.1, &with_limit(0.3));
    assert!(moderate.1.len() >= 160 && moderate.1.len() <= 320);
    assert_closed_sphere_topology(&moderate.0, &moderate.1);
}

#[test]
fn the_cube_is_only_simplified_when_tilting_faces_is_allowed() {
    let (vertices, triangles) = cube();
    let loose = QemConfig {
        preserve_boundary: false,
        boundary_weight: 0.0,
        ..config(6)
    };
    let simplified = qem_simplify(&vertices, &triangles, &loose);
    assert!(simplified.1.len() < 12, "{} triangles", simplified.1.len());
    assert_closed_sphere_topology(&simplified.0, &simplified.1);

    // Every collapse of the cube moves a corner by at least a third of an edge,
    // which tilts a neighbouring face far beyond 0.1 rad.
    let tight = QemConfig {
        max_normal_deviation: 0.1,
        ..loose
    };
    let kept = qem_simplify(&vertices, &triangles, &tight);
    assert_eq!(kept.1.len(), 12);
    assert_eq!(kept.0.len(), 8);
}

// ---------------------------------------------------------------------------
// Boundaries
// ---------------------------------------------------------------------------

/// Boundary edges of a mesh as pairs of positions.
fn boundary_segments(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<([u64; 3], [u64; 3])> {
    let key = |p: [f64; 3]| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
    let mut segments: Vec<([u64; 3], [u64; 3])> = find_boundary_edges(triangles)
        .into_iter()
        .map(|(a, b)| {
            let (ka, kb) = (key(vertices[a]), key(vertices[b]));
            (ka.min(kb), ka.max(kb))
        })
        .collect();
    segments.sort_unstable();
    segments
}

fn signed_area_z(vertices: &[[f64; 3]], tri: &[usize; 3]) -> f64 {
    let (a, b, c) = (vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

#[test]
fn preserving_the_boundary_keeps_every_boundary_vertex_and_edge() {
    let (vertices, triangles) = flat_grid(8, 8, 1.0);
    let simplified = qem_simplify(&vertices, &triangles, &config(40));
    assert!(simplified.1.len() <= 41, "{} triangles", simplified.1.len());
    assert!(simplified.1.len() < triangles.len());

    // The boundary polyline is exactly the original one.
    assert_eq!(
        boundary_segments(&simplified.0, &simplified.1),
        boundary_segments(&vertices, &triangles)
    );
    let boundary_positions = |mesh: &Mesh| -> Vec<[f64; 3]> {
        find_boundary_vertices(&mesh.1)
            .into_iter()
            .map(|v| mesh.0[v])
            .collect()
    };
    assert_eq!(
        sorted_positions(&boundary_positions(&simplified)),
        sorted_positions(&boundary_positions(&(vertices.clone(), triangles.clone())))
    );
    // Nothing folded over and the surface stays flat.
    for tri in &simplified.1 {
        assert!(
            signed_area_z(&simplified.0, tri) > 0.0,
            "folded triangle {tri:?}"
        );
    }
    assert!(simplified.0.iter().all(|p| p[2] == 0.0));
    let area = mesh_surface_area(&simplified.0, &simplified.1);
    assert!((area - 64.0).abs() < 1e-9, "area {area}");
    assert!(validate_topology(&simplified.0, &simplified.1).is_valid);
}

#[test]
fn interior_vertices_are_absorbed_by_pinned_boundary_vertices() {
    let (vertices, triangles) = flat_grid(6, 6, 1.0);
    let simplified = qem_simplify(&vertices, &triangles, &config(1));
    // 24 boundary vertices stay; only interior vertices can disappear.
    let boundary = find_boundary_vertices(&triangles).len();
    assert_eq!(boundary, 24);
    assert!(simplified.0.len() >= boundary);
    assert_eq!(
        boundary_segments(&simplified.0, &simplified.1),
        boundary_segments(&vertices, &triangles)
    );
    assert!(simplified.1.len() < triangles.len() / 2);
    assert_eq!(
        euler_characteristic(simplified.0.len(), &simplified.1),
        1,
        "a disc"
    );
}

#[test]
fn the_boundary_penalty_keeps_corners_and_the_perimeter() {
    let (vertices, triangles) = flat_grid(6, 6, 1.0);
    let corners = [
        [0.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [0.0, 6.0, 0.0],
        [6.0, 6.0, 0.0],
    ];
    let penalised = QemConfig {
        preserve_boundary: false,
        boundary_weight: 1000.0,
        ..config(20)
    };
    let (kept_vertices, kept_triangles) = qem_simplify(&vertices, &triangles, &penalised);
    assert!(kept_triangles.len() <= 21);
    for corner in corners {
        assert!(
            kept_vertices.contains(&corner),
            "corner {corner:?} was lost"
        );
    }
    for v in find_boundary_vertices(&kept_triangles) {
        let p = kept_vertices[v];
        let on_perimeter = [p[0], p[1]]
            .iter()
            .any(|&x| x.abs() < 1e-9 || (x - 6.0).abs() < 1e-9);
        assert!(on_perimeter, "boundary vertex {p:?} left the perimeter");
    }
    assert!(validate_topology(&kept_vertices, &kept_triangles).is_valid);
    for tri in &kept_triangles {
        assert!(signed_area_z(&kept_vertices, tri) > 0.0);
    }

    // Without the penalty a corner is merged to the middle of an edge.
    let free = QemConfig {
        boundary_weight: 0.0,
        ..penalised
    };
    let (free_vertices, _) = qem_simplify(&vertices, &triangles, &free);
    let corners_kept = corners
        .iter()
        .filter(|corner| free_vertices.contains(corner))
        .count();
    assert!(corners_kept < 4, "all corners survive without a penalty");
}

// ---------------------------------------------------------------------------
// Texture seams
// ---------------------------------------------------------------------------

/// Seam along the grid line `u = 3` of a 6 x 6 grid.
fn seam_line() -> TextureSeamFlags {
    let mut seams = TextureSeamFlags::new();
    for j in 0..6 {
        seams.mark_seam(j * 7 + 3, (j + 1) * 7 + 3);
    }
    seams
}

fn seam_config() -> QemConfig {
    QemConfig {
        preserve_boundary: false,
        preserve_texture_seams: true,
        boundary_weight: 0.0,
        ..config(18)
    }
}

#[test]
fn texture_seams_survive_simplification() {
    let (vertices, triangles) = flat_grid(6, 6, 1.0);
    let seams = seam_line();
    let (sv, st) = qem_simplify_with_seams(&vertices, &triangles, &seam_config(), Some(&seams));
    assert!(st.len() <= 19, "{} triangles", st.len());
    let at = |x: f64, y: f64| sv.iter().position(|p| *p == [x, y, 0.0]);
    let edges = edge_usage(&st);
    for j in 0..6 {
        let a = at(3.0, j as f64);
        let b = at(3.0, (j + 1) as f64);
        assert!(
            a.is_some() && b.is_some(),
            "seam vertex of row {j} was lost"
        );
        if let (Some(a), Some(b)) = (a, b) {
            assert!(
                edges.contains_key(&(a.min(b), a.max(b))),
                "seam edge {j} was lost"
            );
        }
    }
    assert!(validate_topology(&sv, &st).is_valid);
    for tri in &st {
        assert!(signed_area_z(&sv, tri) > 0.0);
    }
}

#[test]
fn seams_matter_only_when_they_are_preserved() {
    let (vertices, triangles) = flat_grid(6, 6, 1.0);
    let seams = seam_line();
    let preserved = qem_simplify_with_seams(&vertices, &triangles, &seam_config(), Some(&seams));
    let ignored = qem_simplify(&vertices, &triangles, &seam_config());
    assert_ne!(preserved, ignored, "the seam had no effect");
    // Seams given but not requested are ignored.
    let not_requested = QemConfig {
        preserve_texture_seams: false,
        ..seam_config()
    };
    assert_eq!(
        qem_simplify_with_seams(&vertices, &triangles, &not_requested, Some(&seams)),
        qem_simplify(&vertices, &triangles, &not_requested)
    );
    // Without the seam some seam vertex disappears.
    let lost = (0..=6).any(|j| !ignored.0.contains(&[3.0, j as f64, 0.0]));
    assert!(lost, "the unconstrained run keeps the whole seam by chance");
}

#[test]
fn seam_edges_that_are_not_in_the_mesh_are_ignored() {
    let (vertices, triangles) = flat_grid(4, 4, 1.0);
    let mut seams = TextureSeamFlags::new();
    seams.mark_seam(0, 24); // opposite corners: not an edge
    seams.mark_seam(2, 1000); // vertex out of range
    let config = QemConfig {
        preserve_boundary: false,
        preserve_texture_seams: true,
        boundary_weight: 0.0,
        ..config(10)
    };
    let with = qem_simplify_with_seams(&vertices, &triangles, &config, Some(&seams));
    let without = qem_simplify_with_seams(&vertices, &triangles, &config, None);
    assert_eq!(with, without);
}

// ---------------------------------------------------------------------------
// Robustness
// ---------------------------------------------------------------------------

#[test]
fn invalid_triangles_are_ignored() {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [f64::NAN, 0.0, 0.0],
        [5.0, 5.0, 5.0],
    ];
    let triangles = vec![
        [0, 1, 2],
        [1, 3, 2],
        [0, 0, 1],  // repeated vertex
        [0, 1, 99], // missing vertex
        [0, 1, 4],  // not finite
    ];
    let (sv, st) = qem_simplify(&vertices, &triangles, &config(0));
    assert_eq!(st.len(), 2);
    // Vertices 4 and 5 are not used by a valid triangle.
    assert_eq!(sv.len(), 4);
    assert!(qem_simplify(&[], &[], &config(0)).0.is_empty());
    assert!(qem_simplify(&vertices, &[], &config(0)).1.is_empty());
}

#[test]
fn a_target_above_the_triangle_count_changes_nothing() {
    let (vertices, triangles) = flat_grid(3, 3, 1.0);
    let (sv, st) = qem_simplify(&vertices, &triangles, &config(1000));
    assert_eq!(sv, vertices);
    assert_eq!(st, triangles);
}

#[test]
fn non_manifold_edges_are_left_alone() {
    // Three triangles share the edge (0, 1).
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, -1.0, 0.0],
        [0.5, 0.0, 1.0],
    ];
    let triangles = vec![[0, 1, 2], [0, 1, 3], [0, 1, 4]];
    let free = QemConfig {
        preserve_boundary: false,
        boundary_weight: 0.0,
        ..config(1)
    };
    let (sv, st) = qem_simplify(&vertices, &triangles, &free);
    assert_eq!(st.len(), 3);
    assert_eq!(sv.len(), 5);
}

#[test]
fn a_single_triangle_is_never_collapsed_away() {
    let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let triangles = vec![[0, 1, 2]];
    let free = QemConfig {
        preserve_boundary: false,
        boundary_weight: 0.0,
        ..config(0)
    };
    let (sv, st) = qem_simplify(&vertices, &triangles, &free);
    assert_eq!(st, triangles);
    assert_eq!(sv, vertices);
}

// ---------------------------------------------------------------------------
// Quadrics
// ---------------------------------------------------------------------------

#[test]
fn the_optimal_point_does_not_depend_on_the_scale_of_the_quadric() {
    // Three orthogonal planes meet at (1, 2, 3); scaling the quadric must not turn
    // a regular system into a "singular" one.
    for scale in [1.0, 1.0e-4, 1.0e-8, 1.0e6] {
        let q = SymMat4::from_plane(1.0, 0.0, 0.0, -1.0)
            .add(&SymMat4::from_plane(0.0, 1.0, 0.0, -2.0))
            .add(&SymMat4::from_plane(0.0, 0.0, 1.0, -3.0))
            .scale(scale);
        let p = q
            .optimal_point()
            .unwrap_or_else(|| panic!("scale {scale}: no optimum"));
        assert!((p[0] - 1.0).abs() < 1e-9, "{p:?}");
        assert!((p[1] - 2.0).abs() < 1e-9, "{p:?}");
        assert!((p[2] - 3.0).abs() < 1e-9, "{p:?}");
    }
}

#[test]
fn oblique_degenerate_quadrics_have_no_optimal_point() {
    // One plane (rank 1) and two planes (rank 2) with oblique normals are singular
    // although rounding errors make their determinant slightly non-zero.
    let n = [1.0 / 3.0_f64.sqrt(); 3];
    let plane = SymMat4::from_plane(n[0], n[1], n[2], -0.5);
    let mut one = SymMat4::zero();
    for _ in 0..6 {
        one = one.add(&plane);
    }
    assert!(one.optimal_point().is_none());

    let m = [
        2.0 / 21.0_f64.sqrt(),
        1.0 / 21.0_f64.sqrt(),
        4.0 / 21.0_f64.sqrt(),
    ];
    let two = one.add(&SymMat4::from_plane(m[0], m[1], m[2], 0.25));
    assert!(two.optimal_point().is_none());
}

// ---------------------------------------------------------------------------
// check_normal_deviation and the half-edge mesh
// ---------------------------------------------------------------------------

fn grid_half_edge_mesh() -> (HalfEdgeMesh, Vec<[usize; 3]>) {
    let (vertices, triangles) = flat_grid(4, 4, 1.0);
    (
        HalfEdgeMesh::from_triangles(vertices, &triangles),
        triangles,
    )
}

#[test]
fn vertex_faces_returns_every_face_of_every_vertex() {
    for (vertices, triangles) in [
        flat_grid(4, 4, 1.0),
        flat_grid(1, 1, 1.0),
        cube(),
        icosphere(1),
    ] {
        let mesh = HalfEdgeMesh::from_triangles(vertices.clone(), &triangles);
        for v in 0..vertices.len() {
            let mut found = mesh.vertex_faces(v);
            found.sort_unstable();
            let expected: Vec<usize> = triangles
                .iter()
                .enumerate()
                .filter(|(_, tri)| tri.contains(&v))
                .map(|(f, _)| f)
                .collect();
            assert_eq!(found, expected, "faces of vertex {v}");

            let mut neighbours = mesh.vertex_neighbors(v);
            neighbours.sort_unstable();
            let mut expected_neighbours: Vec<usize> = expected
                .iter()
                .flat_map(|&f| triangles[f])
                .filter(|&u| u != v)
                .collect();
            expected_neighbours.sort_unstable();
            expected_neighbours.dedup();
            assert_eq!(neighbours, expected_neighbours, "neighbours of vertex {v}");
        }
    }
}

#[test]
fn deviation_check_accepts_a_collapse_that_stays_in_the_plane() {
    let (mesh, _) = grid_half_edge_mesh();
    // Vertices 12 = (2, 2) and 13 = (3, 2); the middle of the edge.
    assert!(check_normal_deviation(&mesh, 12, 13, [2.5, 2.0, 0.0], 0.1));
    assert!(check_normal_deviation(&mesh, 12, 13, [2.5, 2.0, 0.0], 1.5));
}

#[test]
fn deviation_check_rejects_a_collapse_that_tilts_faces() {
    let (mesh, _) = grid_half_edge_mesh();
    // Lifting the merged vertex by 0.2 over faces of size one tilts them by more
    // than 0.1 rad but by less than one radian.
    let lifted = [2.5, 2.0, 0.2];
    assert!(!check_normal_deviation(&mesh, 12, 13, lifted, 0.1));
    assert!(check_normal_deviation(&mesh, 12, 13, lifted, 1.0));
}

#[test]
fn deviation_check_rejects_folded_and_degenerate_faces() {
    let (mesh, _) = grid_half_edge_mesh();
    // Moving the merged vertex far to the left folds the face (1, 2), (2, 2),
    // (2, 3) over; the angular limit of 3 rad does not hide a fold.
    assert!(!check_normal_deviation(
        &mesh,
        12,
        13,
        [-3.0, 2.0, 0.0],
        3.0
    ));
    // On the line through (2, 1) and (3, 1) the face (2, 1), (3, 1), (3, 2)
    // degenerates to a segment.
    assert!(!check_normal_deviation(&mesh, 12, 13, [2.7, 1.0, 0.0], 3.0));
}

#[test]
fn deviation_check_without_a_limit_accepts_everything() {
    let (mesh, _) = grid_half_edge_mesh();
    for limit in [0.0, -1.0, f64::NAN] {
        assert!(check_normal_deviation(
            &mesh,
            12,
            13,
            [-30.0, 2.0, 7.0],
            limit
        ));
    }
}

#[test]
fn deviation_check_rejects_unusable_vertices() {
    let (mut mesh, _) = grid_half_edge_mesh();
    assert!(!check_normal_deviation(
        &mesh,
        12,
        999,
        [2.5, 2.0, 0.0],
        0.5
    ));
    mesh.vertex_deleted[13] = true;
    assert!(!check_normal_deviation(&mesh, 12, 13, [2.5, 2.0, 0.0], 0.5));
}

#[test]
fn deviation_check_sees_the_faces_on_both_sides_of_a_boundary_vertex() {
    let (vertices, triangles) = flat_grid(4, 4, 1.0);
    let mesh = HalfEdgeMesh::from_triangles(vertices, &triangles);
    // Vertex 2 = (2, 0) is on the boundary and has faces on both sides. Moving it
    // upwards tilts all of them; the check must find every face, whichever
    // outgoing half-edge is stored for the vertex.
    for v in [2, 10, 5, 14, 0, 24] {
        let faces = mesh.vertex_faces(v).len();
        let expected = triangles.iter().filter(|tri| tri.contains(&v)).count();
        assert_eq!(faces, expected, "vertex {v}");
    }
    let neighbour = 3; // (3, 0)
    assert!(!check_normal_deviation(
        &mesh,
        2,
        neighbour,
        [2.5, 0.0, 0.5],
        0.1
    ));
    // A move to the left folds only the faces on the left of vertex 2.
    assert!(!check_normal_deviation(
        &mesh,
        2,
        neighbour,
        [0.2, 0.0, 0.0],
        3.0
    ));
}
