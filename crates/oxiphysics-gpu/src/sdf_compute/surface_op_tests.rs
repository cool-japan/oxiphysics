// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Behavioural tests of the SDF surface operators: marching cubes, morphological
//! opening/closing (with distance re-initialisation) and Laplacian smoothing /
//! sharpening.  The table-level tests live in `mc_tables.rs` and the Eikonal
//! solver tests in `eikonal.rs`.

use std::collections::HashMap;
use std::f64::consts::PI;

use super::*;

// ── helpers ──────────────────────────────────────────────────────────────────

fn sphere_grid(n: usize, dx: f64, center: [f64; 3], radius: f64) -> SdfGrid {
    let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    grid.compute_sphere_sdf(center, radius);
    grid
}

fn box_grid(n: usize, dx: f64, center: [f64; 3], half: [f64; 3]) -> SdfGrid {
    let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    grid.compute_box_sdf(center, half);
    grid
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Vertex key after quantising to `1e-9`, so bit-level noise cannot split a vertex.
fn quantise(p: [f64; 3]) -> [i64; 3] {
    let q = |x: f64| (x / 1e-9).round() as i64;
    [q(p[0]), q(p[1]), q(p[2])]
}

/// Topology and geometry of a triangle soup.
struct MeshStats {
    vertices: usize,
    edges: usize,
    faces: usize,
    /// Smallest / largest number of triangles sharing an undirected edge.
    min_edge_use: usize,
    max_edge_use: usize,
    /// Largest number of times one *directed* edge occurs (1 for an oriented manifold).
    max_directed_use: usize,
    degenerate: usize,
    components: usize,
    area: f64,
    /// Signed volume by the divergence theorem (positive for outward-facing triangles).
    volume: f64,
}

fn mesh_stats(triangles: &[Triangle]) -> MeshStats {
    let mut ids: HashMap<[i64; 3], usize> = HashMap::new();
    let mut parent: Vec<usize> = Vec::new();
    let mut undirected: HashMap<(usize, usize), usize> = HashMap::new();
    let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
    let (mut area, mut volume, mut degenerate) = (0.0, 0.0, 0);

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for triangle in triangles {
        let mut corner = [0_usize; 3];
        for (slot, &p) in corner.iter_mut().zip(&triangle.v) {
            let next = ids.len();
            let id = *ids.entry(quantise(p)).or_insert(next);
            if id == parent.len() {
                parent.push(id);
            }
            *slot = id;
        }
        if corner[0] == corner[1] || corner[1] == corner[2] || corner[0] == corner[2] {
            degenerate += 1;
            continue;
        }
        let [a, b, c] = triangle.v;
        area += 0.5 * norm(cross(sub(b, a), sub(c, a)));
        volume += dot(a, cross(b, c)) / 6.0;
        for s in 0..3 {
            let (u, v) = (corner[s], corner[(s + 1) % 3]);
            *directed.entry((u, v)).or_insert(0) += 1;
            *undirected.entry((u.min(v), u.max(v))).or_insert(0) += 1;
            let (ru, rv) = (find(&mut parent, u), find(&mut parent, v));
            if ru != rv {
                parent[ru] = rv;
            }
        }
    }
    let mut roots = std::collections::HashSet::new();
    for id in 0..parent.len() {
        roots.insert(find(&mut parent, id));
    }
    MeshStats {
        vertices: ids.len(),
        edges: undirected.len(),
        faces: triangles.len() - degenerate,
        min_edge_use: undirected.values().copied().min().unwrap_or(0),
        max_edge_use: undirected.values().copied().max().unwrap_or(0),
        max_directed_use: directed.values().copied().max().unwrap_or(0),
        degenerate,
        components: roots.len(),
        area,
        volume,
    }
}

fn euler_characteristic(stats: &MeshStats) -> i64 {
    stats.vertices as i64 - stats.edges as i64 + stats.faces as i64
}

fn assert_closed_oriented_manifold(stats: &MeshStats) {
    assert_eq!(stats.degenerate, 0, "degenerate triangles");
    assert_eq!(
        stats.min_edge_use, 2,
        "an edge is used by fewer than two triangles"
    );
    assert_eq!(
        stats.max_edge_use, 2,
        "an edge is used by more than two triangles"
    );
    assert_eq!(
        stats.max_directed_use, 1,
        "neighbouring triangles disagree on orientation"
    );
}

// ── marching cubes ───────────────────────────────────────────────────────────

#[test]
fn marching_cubes_sphere_is_a_watertight_oriented_surface_with_the_right_size() {
    let radius = 1.2;
    let grid = sphere_grid(34, 0.1, [1.713, 1.679, 1.707], radius);
    let triangles = marching_cubes(&grid, 0.0);
    let stats = mesh_stats(&triangles);
    assert_closed_oriented_manifold(&stats);
    assert_eq!(
        euler_characteristic(&stats),
        2,
        "a sphere has Euler characteristic 2"
    );
    assert_eq!(stats.components, 1);
    let area = 4.0 * PI * radius * radius;
    let volume = 4.0 / 3.0 * PI * radius.powi(3);
    assert!(
        (stats.area / area - 1.0).abs() < 0.02,
        "area {} vs {area}",
        stats.area
    );
    assert!(
        stats.volume > 0.0,
        "triangles must face outwards, volume {}",
        stats.volume
    );
    assert!(
        (stats.volume / volume - 1.0).abs() < 0.03,
        "volume {} vs {volume}",
        stats.volume
    );
    assert_eq!(mesh_triangle_count(&grid, 0.0), triangles.len());
}

#[test]
fn marching_cubes_vertices_lie_on_the_level_set_and_normals_follow_the_gradient() {
    let center = [1.713, 1.679, 1.707];
    let radius = 1.2;
    let dx = 0.1;
    let triangles = marching_cubes(&sphere_grid(34, dx, center, radius), 0.0);
    assert!(!triangles.is_empty());
    for triangle in &triangles {
        for &p in &triangle.v {
            let d = norm(sub(p, center)) - radius;
            assert!(
                d.abs() < 0.03 * dx,
                "vertex {p:?} is {d} away from the sphere"
            );
        }
        let n = cross(
            sub(triangle.v[1], triangle.v[0]),
            sub(triangle.v[2], triangle.v[0]),
        );
        let centroid = [
            (triangle.v[0][0] + triangle.v[1][0] + triangle.v[2][0]) / 3.0,
            (triangle.v[0][1] + triangle.v[1][1] + triangle.v[2][1]) / 3.0,
            (triangle.v[0][2] + triangle.v[1][2] + triangle.v[2][2]) / 3.0,
        ];
        let outward = sub(centroid, center);
        let cos = dot(n, outward) / (norm(n) * norm(outward));
        assert!(cos > 0.0, "normal points inwards (cos {cos})");
    }
}

#[test]
fn marching_cubes_torus_is_a_closed_surface_of_genus_one() {
    let (major, minor) = (0.8, 0.3);
    let center = [1.213, 1.187, 1.203];
    let mut grid = SdfGrid::new(40, 40, 40, 0.06, [0.0; 3]);
    grid.compute_torus_sdf(center, major, minor);
    let stats = mesh_stats(&marching_cubes(&grid, 0.0));
    assert_closed_oriented_manifold(&stats);
    assert_eq!(
        euler_characteristic(&stats),
        0,
        "a torus has Euler characteristic 0"
    );
    assert_eq!(stats.components, 1);
    let volume = 2.0 * PI * PI * major * minor * minor;
    assert!(
        (stats.volume / volume - 1.0).abs() < 0.05,
        "volume {} vs {volume}",
        stats.volume
    );
    let area = 4.0 * PI * PI * major * minor;
    assert!(
        (stats.area / area - 1.0).abs() < 0.04,
        "area {} vs {area}",
        stats.area
    );
}

#[test]
fn marching_cubes_two_separate_blobs_form_two_closed_components() {
    let dx = 0.1;
    let a = sphere_grid(40, dx, [1.013, 1.987, 1.969], 0.7);
    let b = sphere_grid(40, dx, [2.987, 2.013, 2.011], 0.9);
    let stats = mesh_stats(&marching_cubes(&union_sdf(&a, &b), 0.0));
    assert_closed_oriented_manifold(&stats);
    assert_eq!(stats.components, 2);
    assert_eq!(euler_characteristic(&stats), 4, "two spheres: 2 + 2");
    let volume = 4.0 / 3.0 * PI * (0.7_f64.powi(3) + 0.9_f64.powi(3));
    assert!((stats.volume / volume - 1.0).abs() < 0.05);
}

#[test]
fn marching_cubes_merged_blobs_with_a_neck_form_one_component() {
    let dx = 0.1;
    let a = sphere_grid(40, dx, [1.513, 2.011, 1.989], 0.9);
    let b = sphere_grid(40, dx, [2.487, 1.989, 2.011], 0.9);
    let stats = mesh_stats(&marching_cubes(&union_sdf(&a, &b), 0.0));
    assert_closed_oriented_manifold(&stats);
    assert_eq!(stats.components, 1);
    assert_eq!(euler_characteristic(&stats), 2);
}

#[test]
fn marching_cubes_planar_field_is_a_flat_sheet_with_exact_area() {
    let (n, dx, x0) = (9, 0.25, 0.9);
    let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                grid.set(i, j, k, grid.world_pos(i, j, k)[0] - x0);
            }
        }
    }
    let triangles = marching_cubes(&grid, 0.0);
    let cells_per_side = n - 1;
    assert_eq!(triangles.len(), 2 * cells_per_side * cells_per_side);
    let side = cells_per_side as f64 * dx;
    let stats = mesh_stats(&triangles);
    assert!(
        (stats.area - side * side).abs() < 1e-12,
        "area {}",
        stats.area
    );
    for triangle in &triangles {
        for p in triangle.v {
            assert!((p[0] - x0).abs() < 1e-12, "vertex off the plane: {p:?}");
        }
        let n = cross(
            sub(triangle.v[1], triangle.v[0]),
            sub(triangle.v[2], triangle.v[0]),
        );
        let unit = [n[0] / norm(n), n[1] / norm(n), n[2] / norm(n)];
        assert!(
            (unit[0] - 1.0).abs() < 1e-12 && unit[1].abs() < 1e-12 && unit[2].abs() < 1e-12,
            "normal must point to increasing x, got {unit:?}"
        );
    }
}

#[test]
fn marching_cubes_single_inside_corner_cuts_exactly_one_triangle() {
    let mut grid = SdfGrid::new(2, 2, 2, 0.5, [0.0; 3]);
    grid.values.iter_mut().for_each(|v| *v = 0.5);
    grid.set(0, 0, 0, -0.1);
    let triangles = marching_cubes(&grid, 0.0);
    assert_eq!(triangles.len(), 1);
    // The crossings sit at t = 0.1 / 0.6 = 1/6 of the edge from the inside corner.
    let corner = grid.world_pos(0, 0, 0);
    let step = 0.5 / 6.0;
    let expected = [
        [corner[0] + step, corner[1], corner[2]],
        [corner[0], corner[1] + step, corner[2]],
        [corner[0], corner[1], corner[2] + step],
    ];
    for want in expected {
        assert!(
            triangles[0].v.iter().any(|&p| norm(sub(p, want)) < 1e-12),
            "missing vertex {want:?} in {:?}",
            triangles[0].v
        );
    }
    let normal = cross(
        sub(triangles[0].v[1], triangles[0].v[0]),
        sub(triangles[0].v[2], triangles[0].v[0]),
    );
    assert!(
        dot(normal, [1.0, 1.0, 1.0]) > 0.0,
        "the triangle must face away from the inside corner, normal {normal:?}"
    );
}

#[test]
fn marching_cubes_extracts_the_requested_isovalue() {
    let dx = 0.1;
    let (center, radius, iso) = ([1.713, 1.679, 1.707], 0.9, 0.3);
    let stats = mesh_stats(&marching_cubes(&sphere_grid(34, dx, center, radius), iso));
    assert_closed_oriented_manifold(&stats);
    let area = 4.0 * PI * (radius + iso).powi(2);
    assert!(
        (stats.area / area - 1.0).abs() < 0.02,
        "area {} vs {area}",
        stats.area
    );
}

#[test]
fn marching_cubes_of_uniform_or_too_small_grids_is_empty() {
    let mut grid = SdfGrid::new(4, 4, 4, 0.1, [0.0; 3]);
    grid.values.iter_mut().for_each(|v| *v = -1.0);
    assert!(marching_cubes(&grid, 0.0).is_empty());
    grid.values.iter_mut().for_each(|v| *v = 1.0);
    assert!(marching_cubes(&grid, 0.0).is_empty());
    let flat = SdfGrid::new(5, 5, 1, 0.1, [0.0; 3]);
    assert!(marching_cubes(&flat, 0.0).is_empty());
}

#[test]
fn marching_cubes_output_is_deterministic() {
    let grid = sphere_grid(24, 0.1, [1.213, 1.187, 1.203], 0.8);
    let first = marching_cubes(&grid, 0.0);
    let second = marching_cubes(&grid, 0.0);
    assert_eq!(first.len(), second.len());
    for (a, b) in first.iter().zip(&second) {
        assert_eq!(a.v, b.v);
    }
}

// ── redistancing ─────────────────────────────────────────────────────────────

#[test]
fn redistance_recovers_a_unit_gradient_sphere_from_a_distorted_field() {
    let (n, dx) = (48, 0.05);
    let center = [1.213, 1.187, 1.2];
    let exact = sphere_grid(n, dx, center, 0.7);
    // Same zero level set, wrong metric: |grad| = 3 near the surface.
    let mut distorted = sphere_grid(n, dx, center, 0.7);
    distorted
        .values
        .iter_mut()
        .for_each(|v| *v = 3.0 * *v + 2.0 * v.powi(3));
    let fixed = sdf_redistance(&distorted);
    let mut worst_near = 0.0_f64;
    for (got, want) in fixed.values.iter().zip(&exact.values) {
        assert_eq!(*got < 0.0, *want < 0.0, "sign changed");
        if want.abs() < 4.0 * dx {
            worst_near = worst_near.max((got - want).abs());
        }
    }
    assert!(
        worst_near < 0.2 * dx,
        "error {worst_near} within four cells of the surface"
    );
    let mut gradient_range = (f64::MAX, 0.0_f64);
    for i in 1..n - 1 {
        for j in 1..n - 1 {
            for k in 1..n - 1 {
                let depth = exact.get(i, j, k).abs();
                if depth > 2.0 * dx && depth < 8.0 * dx {
                    let g = norm(fixed.gradient_at(i, j, k));
                    gradient_range = (gradient_range.0.min(g), gradient_range.1.max(g));
                }
            }
        }
    }
    assert!(
        gradient_range.0 > 0.85 && gradient_range.1 < 1.15,
        "gradient magnitude left [0.85, 1.15]: {gradient_range:?}"
    );
}

#[test]
fn redistance_is_idempotent_on_its_own_output() {
    let (n, dx) = (40, 0.05);
    let fixed = sdf_redistance(&sphere_grid(n, dx, [1.013, 0.987, 1.01], 0.6));
    let again = sdf_redistance(&fixed);
    for (a, b) in again.values.iter().zip(&fixed.values) {
        assert!(
            (a - b).abs() < 0.05 * dx,
            "redistancing moved a value by {}",
            (a - b).abs()
        );
    }
}

#[test]
fn redistance_makes_both_sides_of_a_two_blob_field_consistent() {
    let dx = 0.05;
    let a = sphere_grid(48, dx, [0.613, 1.187, 1.2], 0.35);
    let b = sphere_grid(48, dx, [1.787, 1.213, 1.2], 0.35);
    let union = union_sdf(&a, &b);
    let mut scrambled = sdf_redistance(&union);
    // Destroy every magnitude but keep the sign pattern.
    scrambled
        .values
        .iter_mut()
        .for_each(|v| *v = v.signum() * 5.0);
    let fixed = sdf_redistance(&scrambled);
    // The mid-plane between the blobs is as far from either as the exact union says.
    let (mi, mj, mk) = (24, 24, 24);
    assert!(
        (fixed.get(mi, mj, mk) - union.get(mi, mj, mk)).abs() < 0.5 * dx,
        "{} vs {}",
        fixed.get(mi, mj, mk),
        union.get(mi, mj, mk)
    );
}

// ── morphology ───────────────────────────────────────────────────────────────

const N_MORPH: usize = 64;
const DX_MORPH: f64 = 0.025;

/// Radius of the test bodies: not a multiple of the cell size, so that no grid
/// point lies exactly on a surface (which would make marching cubes emit
/// zero-area triangles).
const BODY_RADIUS: f64 = 0.3011;

/// Sphere of radius [`BODY_RADIUS`] with a thin spike (3 cells wide) along +x and
/// a thick protrusion (10 cells wide) along -x.
fn spiked_sphere() -> (SdfGrid, [f64; 3]) {
    let probe = SdfGrid::new(N_MORPH, N_MORPH, N_MORPH, DX_MORPH, [0.0; 3]);
    let center = probe.world_pos(32, 32, 32);
    let body = sphere_grid(N_MORPH, DX_MORPH, center, BODY_RADIUS);
    let thin = box_grid(
        N_MORPH,
        DX_MORPH,
        [center[0] + 0.275, center[1], center[2]],
        [0.275, 1.5 * DX_MORPH, 1.5 * DX_MORPH],
    );
    let thick = box_grid(
        N_MORPH,
        DX_MORPH,
        [center[0] - 0.3, center[1], center[2]],
        [0.3, 5.0 * DX_MORPH, 5.0 * DX_MORPH],
    );
    (union_sdf(&union_sdf(&body, &thin), &thick), center)
}

#[test]
fn open_removes_a_thin_spike_but_keeps_the_body_and_a_thick_protrusion() {
    let (grid, center) = spiked_sphere();
    let offset = 3.0 * DX_MORPH;
    let opened = sdf_open(&grid, offset);
    let at = |g: &SdfGrid, di: i32| {
        let i = usize::try_from(32 + di).expect("index stays inside the grid");
        g.get(i, 32, 32)
    };

    // Spike cells further than the body radius plus a cell are inside before, outside after.
    let mut probed = 0;
    for di in 0..22_usize {
        let i = 32 + di;
        let distance_from_centre = grid.world_pos(i, 32, 32)[0] - center[0];
        if distance_from_centre > BODY_RADIUS + DX_MORPH {
            assert!(grid.get(i, 32, 32) < 0.0, "cell {i} is part of the spike");
            assert!(
                opened.get(i, 32, 32) > 0.0,
                "cell {i} of the thin spike survived the opening: {}",
                opened.get(i, 32, 32)
            );
            probed += 1;
        }
    }
    assert!(
        probed >= 8,
        "the test must probe the spike, got {probed} cells"
    );
    // Body centre stays inside; the output genuinely differs from the input.
    assert!(at(&opened, 0) < 0.0);
    let changed = grid
        .values
        .iter()
        .zip(&opened.values)
        .filter(|(a, b)| (**a - **b).abs() > DX_MORPH)
        .count();
    assert!(
        changed > 20,
        "opening left the field essentially unchanged ({changed} cells)"
    );
    // The thick protrusion (half-width 5 cells > offset) keeps its axis inside.
    for di in [-6, -10, -14] {
        assert!(
            at(&opened, di) < 0.0,
            "thick protrusion lost at offset {di}"
        );
    }
    // Anti-extensive: nothing outside the shape gets inside.  Outside the union the
    // min() of the component distances is the exact distance, so the opened
    // distance may not be smaller; inside it is not (the union field is too
    // shallow where the nearest boundary of a component is no boundary of the
    // union), and only the sign is compared.
    for (before, after) in grid.values.iter().zip(&opened.values) {
        if *before > 0.0 {
            assert!(
                *after >= *before - 0.75 * DX_MORPH,
                "opening grew the shape: {before} -> {after}"
            );
        }
        if *before > 0.2 * DX_MORPH {
            assert!(
                *after > 0.0,
                "an outside cell ({before}) became inside after opening"
            );
        }
    }
}

#[test]
fn open_reproduces_a_shape_that_the_structuring_element_fits_into() {
    let center = [0.8125; 3];
    let grid = sphere_grid(N_MORPH, DX_MORPH, center, 0.4);
    let opened = sdf_open(&grid, 3.0 * DX_MORPH);
    let closed = sdf_close(&grid, 3.0 * DX_MORPH);
    for result in [&opened, &closed] {
        for (before, after) in grid.values.iter().zip(&result.values) {
            if before.abs() < 4.0 * DX_MORPH {
                assert!(
                    (after - before).abs() < 0.25 * DX_MORPH,
                    "a convex ball must survive the operation, {before} -> {after}"
                );
            }
        }
        let a = count_interior_cells(&grid) as f64;
        let b = count_interior_cells(result) as f64;
        assert!(
            (a - b).abs() / a < 0.01,
            "volume changed from {a} to {b} cells"
        );
    }
}

#[test]
fn open_is_idempotent() {
    let (grid, _) = spiked_sphere();
    let offset = 3.0 * DX_MORPH;
    let once = sdf_open(&grid, offset);
    let twice = sdf_open(&once, offset);
    for (a, b) in once.values.iter().zip(&twice.values) {
        // First-order Eikonal accuracy: a few tenths of a cell near ridges.
        if a.abs() < 3.0 * DX_MORPH {
            assert!(
                (a - b).abs() < 0.5 * DX_MORPH,
                "second opening moved the surface: {a} -> {b}"
            );
        }
        if a.abs() > 0.5 * DX_MORPH {
            assert_eq!(
                *a < 0.0,
                *b < 0.0,
                "second opening flipped a cell ({a} -> {b})"
            );
        }
    }
}

/// Two 12-cell-thick slabs facing each other across a gap of `gap_cells` cells.
fn two_slabs(gap_cells: usize) -> SdfGrid {
    let first_end = 20;
    let second_start = first_end + gap_cells;
    let low = 8.0 * DX_MORPH;
    let node_x = |index: usize| (index as f64 + 0.5) * DX_MORPH;
    // Slab faces sit halfway between cell centres so no sample lies on a face.
    let a_hi = 0.5 * (node_x(first_end - 1) + node_x(first_end));
    let b_lo = 0.5 * (node_x(second_start - 1) + node_x(second_start));
    let b_hi = b_lo + 12.0 * DX_MORPH;
    let half_yz = 0.5 * (N_MORPH as f64 * DX_MORPH) - 4.0 * DX_MORPH;
    let mid_yz = 0.5 * N_MORPH as f64 * DX_MORPH;
    let a = box_grid(
        N_MORPH,
        DX_MORPH,
        [0.5 * (low + a_hi), mid_yz, mid_yz],
        [0.5 * (a_hi - low), half_yz, half_yz],
    );
    let b = box_grid(
        N_MORPH,
        DX_MORPH,
        [0.5 * (b_lo + b_hi), mid_yz, mid_yz],
        [0.5 * (b_hi - b_lo), half_yz, half_yz],
    );
    union_sdf(&a, &b)
}

#[test]
fn close_fills_a_gap_narrower_than_twice_the_offset() {
    let offset = 3.0 * DX_MORPH;
    let grid = two_slabs(3);
    let closed = sdf_close(&grid, offset);
    // Gap cells are x-indices 20, 21, 22.
    for i in 20..23 {
        assert!(grid.get(i, 32, 32) > 0.0, "gap cell {i} must start outside");
        assert!(
            closed.get(i, 32, 32) < 0.0,
            "gap cell {i} must be filled by closing, got {}",
            closed.get(i, 32, 32)
        );
    }
    // Extensive: nothing that was inside is lost, up to resolution.
    for (before, after) in grid.values.iter().zip(&closed.values) {
        assert!(
            *after <= *before + 0.75 * DX_MORPH,
            "closing shrank the shape: {before} -> {after}"
        );
        // The Eikonal solve is first order: the eight sharp slab corners, which sit on a
        // three-way ridge of the dilated shape, come out about half a cell shallower.
        if *before < -DX_MORPH {
            assert!(
                *after < 0.0,
                "an inside cell ({before}) became outside after closing"
            );
        }
    }
    // Far-away empty space stays empty and slab interiors stay solid.
    assert!(closed.get(40, 32, 32) > 0.0);
    assert!(closed.get(32, 2, 2) > 0.0);
    assert!(closed.get(14, 32, 32) < 0.0 && closed.get(28, 32, 32) < 0.0);
}

#[test]
fn close_keeps_a_gap_wider_than_twice_the_offset() {
    let offset = 3.0 * DX_MORPH;
    let grid = two_slabs(10);
    let closed = sdf_close(&grid, offset);
    // Gap cells are x-indices 20..30; the middle ones are further than the offset from both slabs.
    for i in 23..27 {
        assert!(
            closed.get(i, 32, 32) > 0.0,
            "gap cell {i} must stay open, got {}",
            closed.get(i, 32, 32)
        );
    }
}

#[test]
fn open_removes_a_bridge_narrower_than_twice_the_offset() {
    // Two blobs joined by a bridge: open() with a large enough offset separates them.
    let center_a = [0.5125, 0.8125, 0.8125];
    let center_b = [1.1125, 0.8125, 0.8125];
    let a = sphere_grid(N_MORPH, DX_MORPH, center_a, 0.2511);
    let b = sphere_grid(N_MORPH, DX_MORPH, center_b, 0.2511);
    let bridge = box_grid(
        N_MORPH,
        DX_MORPH,
        [0.8125, 0.8125, 0.8125],
        [0.3, 1.5 * DX_MORPH, 1.5 * DX_MORPH],
    );
    let joined = union_sdf(&union_sdf(&a, &b), &bridge);
    let before = mesh_stats(&marching_cubes(&joined, 0.0));
    assert_eq!(before.components, 1);
    let opened = sdf_open(&joined, 3.0 * DX_MORPH);
    let after = mesh_stats(&marching_cubes(&opened, 0.0));
    assert_eq!(after.components, 2, "the bridge must be cut");
    assert_closed_oriented_manifold(&after);
    assert!(
        opened.get(32, 32, 32) > 0.0,
        "bridge midpoint must be outside after opening"
    );
}

#[test]
fn open_of_a_shape_smaller_than_the_structuring_element_is_empty() {
    let grid = sphere_grid(32, 0.05, [0.8, 0.8, 0.8], 0.1);
    assert!(count_interior_cells(&grid) > 0);
    let opened = sdf_open(&grid, 0.15);
    assert_eq!(count_interior_cells(&opened), 0);
    assert!(opened.values.iter().all(|&v| v == f64::MAX));
}

#[test]
fn close_of_a_grid_filled_by_the_dilation_stays_filled() {
    let grid = sphere_grid(16, 0.1, [0.8, 0.8, 0.8], 0.5);
    let closed = sdf_close(&grid, 3.0);
    assert!(closed.values.iter().all(|&v| v == -f64::MAX));
}

#[test]
fn a_negative_offset_exchanges_opening_and_closing() {
    let (grid, _) = spiked_sphere();
    let offset = 2.0 * DX_MORPH;
    let open_negative = sdf_open(&grid, -offset);
    let close_positive = sdf_close(&grid, offset);
    assert_eq!(open_negative.values, close_positive.values);
    let close_negative = sdf_close(&grid, -offset);
    let open_positive = sdf_open(&grid, offset);
    assert_eq!(close_negative.values, open_positive.values);
}

// ── Laplacian smoothing / sharpening ─────────────────────────────────────────

/// Constant-zero grid with one negative spike in the interior.
fn spike_grid(n: usize, dx: f64) -> (SdfGrid, (usize, usize, usize)) {
    let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    grid.values.iter_mut().for_each(|v| *v = 0.0);
    let c = n / 2;
    grid.set(c, c, c, -1.0);
    (grid, (c, c, c))
}

fn max_abs(grid: &SdfGrid) -> f64 {
    grid.values.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

#[test]
fn smoothing_shrinks_a_spike_and_sharpening_grows_it() {
    let (grid, (i, j, k)) = spike_grid(9, 0.1);
    let smoothed = sdf_laplacian_smooth(&grid, 1, 0.001);
    assert!(smoothed.get(i, j, k) > -1.0 && smoothed.get(i, j, k) < 0.0);
    assert!(
        smoothed.get(i + 1, j, k) < 0.0,
        "the spike must diffuse into its neighbours"
    );
    let sharpened = sdf_laplacian_sharpen(&grid, 0.001);
    assert!(
        sharpened.get(i, j, k) < -1.0,
        "sharpening must deepen the spike"
    );
    assert!(
        sharpened.get(i + 1, j, k) > 0.0,
        "sharpening must overshoot the neighbours"
    );
    // Exact one-step values: the Laplacian at the spike is 6 / dx^2 = 600, and
    // -1 / dx^2 = -100 at a neighbour; smoothing adds dt * laplacian, sharpening
    // subtracts amount * laplacian.
    assert!((smoothed.get(i, j, k) - (-1.0 + 0.001 * 600.0)).abs() < 1e-12);
    assert!((smoothed.get(i + 1, j, k) - (0.001 * -100.0)).abs() < 1e-12);
    assert!((sharpened.get(i, j, k) - (-1.0 - 0.001 * 600.0)).abs() < 1e-12);
    assert!((sharpened.get(i + 1, j, k) - (-0.001 * -100.0)).abs() < 1e-12);
}

#[test]
fn smoothing_never_increases_the_maximum_even_for_a_large_time_step() {
    // dt / dx^2 of 1.0, 0.7 and 3.3: up to twenty times the explicit stability
    // limit of one step.  An unstable scheme multiplies the checkerboard mode by
    // 1 - 12 dt / dx^2 per step.
    for (n, dx, dt) in [(11, 0.1, 0.01), (11, 0.1, 0.007), (11, 0.1, 0.033)] {
        let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
        let mut rng = 0x2545_F491_u32;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    rng = rng.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    let noise = f64::from(rng >> 8) / f64::from(1_u32 << 24) * 2.0 - 1.0;
                    let checker = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                    grid.set(i, j, k, 0.5 * checker + 0.5 * noise);
                }
            }
        }
        let mut previous = max_abs(&grid);
        let mut current = sdf_laplacian_smooth(&grid, 1, dt);
        for _ in 0..4 {
            let now = max_abs(&current);
            assert!(
                now <= previous + 1e-12,
                "dt = {dt}: smoothing amplified {previous} -> {now}"
            );
            previous = now;
            current = sdf_laplacian_smooth(&current, 1, dt);
        }
    }
}

#[test]
fn smoothing_damps_the_checkerboard_mode_at_the_exact_substep_rate() {
    // dt / dx^2 = 0.7 -> m = ceil(4.2) = 5 sub-steps of coefficient 0.14, each
    // multiplying an interior checkerboard by 1 - 12 * 0.14 = -0.68.  Three
    // iterations are 15 sub-steps; the centre of a 41^3 grid is 20 cells from the
    // boundary, which the information cannot reach in 15 sub-steps.
    let (n, dx, dt) = (41, 0.1, 0.007);
    let mut grid = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                grid.set(i, j, k, if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 });
            }
        }
    }
    let smoothed = sdf_laplacian_smooth(&grid, 3, dt);
    let centre = n / 2;
    let start = grid.get(centre, centre, centre);
    let expected = start * (-0.68_f64).powi(15);
    assert!(
        (smoothed.get(centre, centre, centre) - expected).abs() < 1e-9,
        "centre {} vs {expected}",
        smoothed.get(centre, centre, centre)
    );
    assert!(expected.abs() < 5e-3);
}

#[test]
fn smoothing_sub_steps_match_manual_small_steps() {
    let (spike, _) = spike_grid(9, 0.1);
    // dt / dx^2 = 0.75 -> m = ceil(4.5) = 5 sub-steps of dt / 5 (ratio 0.15).
    let big = sdf_laplacian_smooth(&spike, 1, 0.0075);
    // Five iterations of dt / 5 are stable on their own (6 * 0.15 < 1) and equal.
    let manual = sdf_laplacian_smooth(&spike, 5, 0.0015);
    for (a, b) in big.values.iter().zip(&manual.values) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }
    // The spike has diffused into its surroundings: strictly between untouched and gone.
    assert!(
        big.get(4, 4, 4) > -0.1 && big.get(4, 4, 4) < -0.01,
        "{}",
        big.get(4, 4, 4)
    );
}

#[test]
fn smoothing_preserves_linear_fields_and_adds_the_laplacian_of_a_quadratic() {
    let (n, dx) = (9, 0.1);
    let mut linear = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    let mut quadratic = SdfGrid::new(n, n, n, dx, [0.0; 3]);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let [x, y, z] = linear.world_pos(i, j, k);
                linear.set(i, j, k, x + 2.0 * y - 3.0 * z);
                quadratic.set(i, j, k, x * x + y * y + z * z);
            }
        }
    }
    let smoothed = sdf_laplacian_smooth(&linear, 5, 0.004);
    for (a, b) in smoothed.values.iter().zip(&linear.values) {
        assert!(
            (a - b).abs() < 1e-12,
            "a linear field has zero Laplacian: {a} vs {b}"
        );
    }
    // One stable sub-step (dt / dx^2 = 0.1 < 1/6) adds dt * 6 in the interior and
    // dt * 4 on a face (one axis skipped).
    let dt = 0.001;
    let once = sdf_laplacian_smooth(&quadratic, 1, dt);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let axes_with_two_neighbours =
                    [i, j, k].iter().filter(|&&c| c > 0 && c < n - 1).count();
                let want = quadratic.get(i, j, k) + dt * 2.0 * axes_with_two_neighbours as f64;
                assert!(
                    (once.get(i, j, k) - want).abs() < 1e-12,
                    "cell ({i},{j},{k}): {} vs {want}",
                    once.get(i, j, k)
                );
            }
        }
    }
}

#[test]
fn mean_curvature_smoothing_shrinks_a_sphere_at_the_analytic_rate() {
    // For phi = r - R the flow phi_t = laplacian(phi) = 2 / r gives phi(r, t) = r - R + 2t / r,
    // so the zero level set moves to r = (R + sqrt(R^2 - 8t)) / 2.
    let (n, dx, radius, t) = (48, 0.05, 0.5, 0.00625);
    let center = SdfGrid::new(n, n, n, dx, [0.0; 3]).world_pos(24, 24, 24);
    let grid = sphere_grid(n, dx, center, radius);
    let smoothed = sdf_mean_curvature_smooth(&grid, t);
    // Find the zero crossing along the +x axis through the centre.
    let c = 24;
    let mut crossing = None;
    for i in c..n - 1 {
        let (a, b) = (smoothed.get(i, c, c), smoothed.get(i + 1, c, c));
        if a < 0.0 && b >= 0.0 {
            let xa = smoothed.world_pos(i, c, c)[0] - center[0];
            crossing = Some(xa + dx * a / (a - b));
            break;
        }
    }
    let crossing = crossing.expect("the shrunken sphere must still cross zero");
    let expected = 0.5 * (radius + (radius * radius - 8.0 * t).sqrt());
    assert!(
        (crossing - expected).abs() < 0.1 * dx,
        "zero crossing at {crossing}, analytic {expected}, started at {radius}"
    );
    assert!(crossing < radius - 0.3 * dx, "the sphere must shrink");
}

#[test]
fn smoothing_and_sharpening_ignore_invalid_parameters() {
    let (spike, _) = spike_grid(5, 0.1);
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            sdf_laplacian_smooth(&spike, 2, dt).values,
            spike.values,
            "dt = {dt}"
        );
    }
    assert_eq!(sdf_laplacian_smooth(&spike, 0, 0.01).values, spike.values);
    let mut broken = SdfGrid::new(3, 3, 3, 0.1, [0.0; 3]);
    broken.values.truncate(5);
    assert_eq!(sdf_laplacian_smooth(&broken, 1, 0.001).values.len(), 5);
    assert_eq!(sdf_laplacian_sharpen(&broken, 0.001).values.len(), 5);
    let zero_spacing = SdfGrid::new(3, 3, 3, 0.0, [0.0; 3]);
    assert_eq!(
        sdf_laplacian_smooth(&zero_spacing, 1, 0.001).values,
        zero_spacing.values
    );
}

#[test]
fn smoothing_a_sphere_keeps_it_a_sphere() {
    let grid = sphere_grid(24, 0.1, [1.213, 1.187, 1.2], 0.8);
    let smoothed = sdf_laplacian_smooth(&grid, 2, 0.002);
    assert!(max_abs(&smoothed) <= max_abs(&grid) + 1e-12);
    for (a, b) in smoothed.values.iter().zip(&grid.values) {
        // The surface moves by 2 t / R = 0.01 at most; only cells that close can flip.
        if b.abs() > 0.03 {
            assert_eq!(
                *a < 0.0,
                *b < 0.0,
                "light smoothing must not flip a cell ({b} -> {a})"
            );
        }
    }
}

// ── module layout ────────────────────────────────────────────────────────────

#[test]
fn fast_sweeping_update_is_reachable_through_every_path() {
    let seeded = || {
        let mut grid = SdfGrid::new(7, 7, 7, 0.5, [0.0; 3]);
        grid.set(3, 3, 3, 0.0);
        grid
    };
    let mut via_functions = seeded();
    let mut via_eikonal = seeded();
    let mut via_reexport = seeded();
    functions::fast_sweeping_update(&mut via_functions);
    eikonal::fast_sweeping_update(&mut via_eikonal);
    fast_sweeping_update(&mut via_reexport);
    assert_eq!(via_functions.values, via_eikonal.values);
    assert_eq!(via_functions.values, via_reexport.values);
    assert!(
        (via_functions.get(3, 3, 6) - 1.5).abs() < 1e-12,
        "three cells along an axis are 1.5 away"
    );
}
