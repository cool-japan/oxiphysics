// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the grid-accelerated vertex-face self-collision
//! resolution.
//!
//! The exhaustive scan that the grid replaced cost
//! `iterations * vertices * triangles` narrow-phase tests (and cloned the
//! triangle list every iteration). The grid must find exactly the same contacts
//! and apply exactly the same corrections, so the tests compare against a
//! verbatim copy of the exhaustive algorithm bit for bit, and check the work
//! done through [`SelfCollisionStats`] instead of wall-clock time.

use super::triangle_grid::{REBUILD_BUDGET_BASE, REBUILD_BUDGET_PARTICLES, TriangleGrid};
use super::{FabricCollision, FabricParticle, SelfCollisionStats, create_fabric_mesh};
use oxiphysics_core::math::{Real, Vec3};

/// Particle spacing of the generated sheets \[m\].
const SPACING: Real = 0.05;

/// Deterministic SplitMix64 generator: reproducible scenarios without
/// depending on an external random number crate.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> Real {
        (self.next_u64() >> 11) as Real / (1_u64 << 53) as Real
    }

    /// Uniform in `[lo, hi)`.
    fn range(&mut self, lo: Real, hi: Real) -> Real {
        lo + (hi - lo) * self.unit()
    }
}

/// `layers` parallel `n x n` sheets stacked `layer_gap` apart, each slid
/// sideways by a fraction of a cell and jittered, with random masses and about
/// one pinned particle in ten. A gap below the thickness makes the sheets
/// penetrate each other all over.
fn stacked_sheets(
    n: usize,
    layers: usize,
    layer_gap: Real,
    thickness: Real,
    seed: u64,
) -> FabricCollision {
    let mut rng = SplitMix64(seed);
    let extent = SPACING * (n - 1) as Real;
    let mut particles: Vec<FabricParticle> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for layer in 0..layers {
        let (sheet, sheet_triangles) = create_fabric_mesh(n, n, extent, extent, 1.0);
        let base = particles.len();
        let slide = 0.37 * SPACING * layer as Real;
        for p in &sheet {
            let position = Vec3::new(
                p.position.x + slide + rng.range(-0.1, 0.1) * SPACING,
                layer as Real * layer_gap + rng.range(-0.25, 0.25) * layer_gap,
                p.position.z + 0.5 * slide + rng.range(-0.1, 0.1) * SPACING,
            );
            particles.push(if rng.unit() < 0.1 {
                FabricParticle::pinned(position)
            } else {
                FabricParticle::new(position, rng.range(0.01, 0.05))
            });
        }
        triangles.extend(
            sheet_triangles
                .iter()
                .map(|t| [t[0] + base, t[1] + base, t[2] + base]),
        );
    }
    FabricCollision::new(particles, triangles, thickness, 0.3)
}

/// A single flat `n x n` sheet in the XZ plane: nothing touches anything.
fn flat_sheet(n: usize, thickness: Real) -> FabricCollision {
    let extent = SPACING * (n - 1) as Real;
    let (particles, triangles) = create_fabric_mesh(n, n, extent, extent, 1.0);
    FabricCollision::new(particles, triangles, thickness, 0.3)
}

/// Verbatim copy of the exhaustive algorithm that the grid replaced: every
/// vertex against every triangle, `collision_iterations` times.
fn resolve_exhaustively(fc: &mut FabricCollision) {
    for _ in 0..fc.collision_iterations {
        let n_verts = fc.particles.len();
        let tris = fc.triangles.clone();
        for vi in 0..n_verts {
            for tri in &tris {
                if let Some((depth, normal)) = fc.vertex_face_test(vi, *tri) {
                    let w_v = fc.particles[vi].inv_mass;
                    let w_tri: Real = tri.iter().map(|&i| fc.particles[i].inv_mass).sum();
                    let total_w = w_v + w_tri / 3.0;
                    if total_w < 1e-12 {
                        continue;
                    }
                    let correction = depth / total_w;
                    if !fc.particles[vi].is_pinned() {
                        fc.particles[vi].position += normal * (correction * w_v);
                    }
                    for &ti in tri.iter() {
                        if !fc.particles[ti].is_pinned() {
                            let w_i = fc.particles[ti].inv_mass;
                            fc.particles[ti].position -= normal * (correction * w_i / 3.0);
                        }
                    }
                }
            }
        }
    }
}

fn bits(v: &Vec3) -> [u64; 3] {
    [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
}

/// Assert that both systems hold bit-identical particle positions.
fn assert_same_positions(label: &str, got: &FabricCollision, expected: &FabricCollision) {
    assert_eq!(got.particles.len(), expected.particles.len(), "{label}");
    for (i, (g, e)) in got
        .particles
        .iter()
        .zip(expected.particles.iter())
        .enumerate()
    {
        assert_eq!(
            bits(&g.position),
            bits(&e.position),
            "{label}: particle {i} differs from the exhaustive scan: {:?} vs {:?}",
            g.position,
            e.position
        );
    }
}

/// Number of (vertex, triangle) pairs that currently pass the narrow phase.
fn contact_pairs(fc: &FabricCollision) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for vertex in 0..fc.particles.len() {
        for (t, tri) in fc.triangles.iter().enumerate() {
            if fc.vertex_face_test(vertex, *tri).is_some() {
                pairs.push((vertex, t));
            }
        }
    }
    pairs
}

// ---------------------------------------------------------------------------
// Equivalence with the exhaustive scan
// ---------------------------------------------------------------------------

#[test]
fn grid_candidates_are_ascending_and_cover_every_contact() {
    for seed in 1..=4 {
        let fc = stacked_sheets(10, 3, 0.008, 0.02, seed);
        let grid = TriangleGrid::build(&fc.particles, &fc.triangles, fc.collision_thickness);
        assert!(grid.is_indexed(), "seed {seed}: the sheets are finite");
        for vertex in 0..fc.particles.len() {
            let list = grid.candidates(vertex);
            assert!(
                list.windows(2).all(|w| w[0] < w[1]),
                "seed {seed}: candidates of vertex {vertex} must be strictly ascending"
            );
        }
        let pairs = contact_pairs(&fc);
        assert!(
            pairs.len() > 50,
            "seed {seed}: the scenario must contain many contacts, got {}",
            pairs.len()
        );
        for (vertex, triangle) in pairs {
            assert!(
                grid.candidates(vertex).contains(&triangle),
                "seed {seed}: contact of vertex {vertex} with triangle {triangle} is not a candidate"
            );
        }
    }
}

#[test]
fn grid_resolution_is_bit_identical_to_the_exhaustive_scan() {
    for (n, layers, gap, thickness, seed) in [
        (8, 2, 0.008, 0.02, 11),
        (10, 3, 0.008, 0.02, 12),
        (12, 3, 0.004, 0.012, 13),
        (9, 4, 0.006, 0.025, 14),
        (14, 2, 0.01, 0.03, 15),
    ] {
        let mut expected = stacked_sheets(n, layers, gap, thickness, seed);
        let mut got = expected.clone();
        resolve_exhaustively(&mut expected);
        let stats = got.resolve_self_collisions_with_stats();
        assert!(
            stats.contacts_resolved > 20,
            "n={n} layers={layers} seed={seed}: scenario must exercise the corrections, got {stats:?}"
        );
        assert_same_positions(
            &format!("n={n} layers={layers} seed={seed}"),
            &got,
            &expected,
        );
    }
}

#[test]
fn resolve_self_collisions_matches_the_stats_variant() {
    let mut plain = stacked_sheets(10, 3, 0.008, 0.02, 21);
    let mut with_stats = plain.clone();
    plain.resolve_self_collisions();
    let stats = with_stats.resolve_self_collisions_with_stats();
    assert!(stats.contacts_resolved > 0);
    assert_same_positions("plain vs stats", &plain, &with_stats);
}

#[test]
fn contacts_are_resolved_for_every_iteration_count() {
    for iterations in [1, 2, 3, 5, 8] {
        let mut expected = stacked_sheets(9, 3, 0.008, 0.02, 31);
        expected.collision_iterations = iterations;
        let mut got = expected.clone();
        resolve_exhaustively(&mut expected);
        got.resolve_self_collisions();
        assert_same_positions(&format!("{iterations} iterations"), &got, &expected);
    }
}

#[test]
fn dense_stack_rebuilds_the_grid_and_stays_exact() {
    // Seven sheets packed 0.1 * thickness apart: every particle is pushed
    // many times, so some move farther than the grid's safety margin and the
    // grid has to be rebuilt while a vertex is being processed.
    let thickness = 0.02;
    let mut expected = stacked_sheets(9, 7, 0.1 * thickness, thickness, 41);
    expected.collision_iterations = 8;
    let mut got = expected.clone();
    resolve_exhaustively(&mut expected);
    let stats = got.resolve_self_collisions_with_stats();
    assert!(
        stats.grid_rebuilds > 0,
        "the dense stack must outgrow the safety margin, got {stats:?}"
    );
    assert_same_positions("dense stack", &got, &expected);
}

#[test]
fn exhausting_the_rebuild_budget_keeps_the_result_exact() {
    // Twelve tiny sheets packed a tenth of the thickness apart: so many
    // corrections that the grid is rebuilt until its budget is used up and the
    // exhaustive scan finishes the call.
    let thickness = 0.02;
    let mut expected = stacked_sheets(3, 12, 0.1 * thickness, thickness, 121);
    expected.collision_iterations = 8;
    let mut got = expected.clone();
    let budget = REBUILD_BUDGET_BASE + got.particles.len() / REBUILD_BUDGET_PARTICLES;
    resolve_exhaustively(&mut expected);
    let stats = got.resolve_self_collisions_with_stats();
    assert_eq!(
        stats.grid_rebuilds, budget,
        "the stack must use up the rebuild budget, got {stats:?}"
    );
    assert_same_positions("rebuild budget exhausted", &got, &expected);
}

#[test]
fn contact_across_a_cell_boundary_is_found() {
    // A vertical triangle just below a cell boundary and a vertex just above
    // it: the contact crosses the boundary, so the vertex looks into a cell
    // other than the one the triangle's own box starts in.
    let thickness = 0.01;
    let corner = |x: Real, y: Real, z: Real| FabricParticle::new(Vec3::new(x, y, z), 1.0);
    let build = |x_plane: Real, x_vertex: Real| {
        FabricCollision::new(
            vec![
                corner(x_plane, 0.0, 0.0),
                corner(x_plane, 1.0, 0.0),
                corner(x_plane, 0.0, 1.0),
                corner(x_vertex, 0.25, 0.25),
            ],
            vec![[0, 1, 2]],
            thickness,
            0.3,
        )
    };
    let probe = build(0.0, 0.5);
    let cell_size = TriangleGrid::build(&probe.particles, &probe.triangles, thickness)
        .cell_size()
        .expect("finite sheet");
    for boundary in [1.0, 7.0, -3.0] {
        let x_boundary = boundary * cell_size;
        let mut fc = build(x_boundary - 0.2 * thickness, x_boundary + 0.3 * thickness);
        let mut expected = fc.clone();
        resolve_exhaustively(&mut expected);
        let stats = fc.resolve_self_collisions_with_stats();
        assert!(
            stats.contacts_resolved >= 1,
            "boundary {boundary}: the contact must be found, got {stats:?}"
        );
        assert!(
            fc.particles[3].position.x > x_boundary + 0.3 * thickness,
            "boundary {boundary}: the vertex must be pushed away from the plane"
        );
        assert_same_positions(&format!("boundary {boundary}"), &fc, &expected);
    }
}

#[test]
fn resolution_far_from_the_origin_is_exact() {
    let mut expected = stacked_sheets(10, 3, 0.008, 0.02, 51);
    let offset = Vec3::new(1.0e6, -2.0e6, 3.0e6);
    for p in expected.particles.iter_mut() {
        p.position += offset;
        p.prev_position = p.position;
    }
    let mut got = expected.clone();
    resolve_exhaustively(&mut expected);
    let stats = got.resolve_self_collisions_with_stats();
    assert!(stats.contacts_resolved > 0, "{stats:?}");
    assert_same_positions("far from the origin", &got, &expected);
}

#[test]
fn resolution_is_deterministic() {
    let mut a = stacked_sheets(10, 3, 0.008, 0.02, 61);
    let mut b = a.clone();
    let stats_a = a.resolve_self_collisions_with_stats();
    let stats_b = b.resolve_self_collisions_with_stats();
    assert_eq!(stats_a, stats_b);
    assert_same_positions("repeat run", &a, &b);
}

// ---------------------------------------------------------------------------
// Cost
// ---------------------------------------------------------------------------

#[test]
fn untouched_sheet_costs_a_fraction_of_the_exhaustive_scan() {
    let n = 40;
    let mut fc = flat_sheet(n, 0.005);
    let vertices = fc.particles.len();
    let triangles = fc.triangles.len();
    let stats = fc.resolve_self_collisions_with_stats();
    assert_eq!(stats.contacts_resolved, 0, "a flat sheet has no contacts");
    assert_eq!(stats.grid_rebuilds, 0);
    let exhaustive_pass = vertices * triangles;
    assert!(
        stats.pairs_tested * 25 < exhaustive_pass,
        "{} narrow-phase tests for {vertices} vertices and {triangles} triangles \
         (exhaustive pass: {exhaustive_pass})",
        stats.pairs_tested
    );
}

#[test]
fn work_grows_linearly_with_the_mesh_size() {
    let work = |n: usize| {
        let mut fc = flat_sheet(n, 0.005);
        fc.resolve_self_collisions_with_stats().pairs_tested as Real
    };
    let small = work(20);
    let large = work(40);
    // Four times the vertices: about four times the work, where the exhaustive
    // scan would take sixteen times as long.
    let ratio = large / small;
    assert!(
        ratio > 2.5 && ratio < 6.0,
        "doubling the resolution in both directions should about quadruple the work, got {ratio} \
         ({small} -> {large})"
    );
}

#[test]
fn stacked_sheets_cost_far_less_than_the_exhaustive_scan() {
    let mut fc = stacked_sheets(30, 3, 0.004, 0.012, 71);
    let exhaustive_pass = fc.particles.len() * fc.triangles.len();
    let iterations = fc.collision_iterations;
    let stats = fc.resolve_self_collisions_with_stats();
    assert!(stats.contacts_resolved > 100, "{stats:?}");
    assert!(
        stats.pairs_tested * 8 < exhaustive_pass * iterations,
        "{} narrow-phase tests, the exhaustive scan would run {} ({stats:?})",
        stats.pairs_tested,
        exhaustive_pass * iterations
    );
}

#[test]
fn iterations_stop_once_nothing_is_left_to_resolve() {
    // Without contacts the second to fifth iteration would repeat the first
    // one verbatim; the work equals one pass over the candidate lists.
    let mut fc = flat_sheet(16, 0.005);
    assert_eq!(fc.collision_iterations, 5);
    let grid = TriangleGrid::build(&fc.particles, &fc.triangles, fc.collision_thickness);
    let one_pass: usize = (0..fc.particles.len())
        .map(|vertex| grid.candidates(vertex).len())
        .sum();
    let stats = fc.resolve_self_collisions_with_stats();
    assert_eq!(stats.pairs_tested, one_pass);
    assert_eq!(stats.contacts_resolved, 0);
}

// ---------------------------------------------------------------------------
// Grid maintenance
// ---------------------------------------------------------------------------

#[test]
fn grid_reports_drift_beyond_the_safety_margin_and_rebuilds() {
    let thickness = 0.02;
    let mut fc = flat_sheet(10, thickness);
    // A free vertex (it belongs to no triangle) hovering over the sheet.
    let mover = fc.particles.len();
    fc.particles
        .push(FabricParticle::new(Vec3::new(0.22, 0.01, 0.22), 1.0));
    let mut grid = TriangleGrid::build(&fc.particles, &fc.triangles, thickness);
    assert!(grid.is_indexed());
    let hovering_candidates = grid.candidates(mover).to_vec();
    assert!(!hovering_candidates.is_empty());
    let corner_candidates = grid.candidates(0).to_vec();
    assert!(!grid.has_drifted(&fc.particles, &[mover, 0]));

    // The triangles are small compared with the thickness, so the safety margin
    // is the thickness: a move inside it is tolerated.
    fc.particles[mover].position.y += 0.9 * thickness;
    assert!(!grid.has_drifted(&fc.particles, &[mover]));
    // Beyond the margin the lists may be incomplete.
    fc.particles[mover].position.y += 0.2 * thickness;
    assert!(grid.has_drifted(&fc.particles, &[mover]));
    assert!(
        !grid.has_drifted(&fc.particles, &[0, 1, 2]),
        "only the moved particle drifted"
    );

    // Teleport the particle far away: after the rebuild nothing is near it,
    // it is no longer drifted, and everybody else keeps their candidates.
    fc.particles[mover].position = Vec3::new(50.0, 50.0, 50.0);
    grid.rebuild(&fc.particles, &fc.triangles, thickness);
    assert!(grid.is_indexed());
    assert!(!grid.has_drifted(&fc.particles, &[mover]));
    assert!(grid.candidates(mover).is_empty());
    assert_eq!(grid.candidates(0), corner_candidates.as_slice());
}

#[test]
fn grid_margin_grows_with_the_triangle_size_for_thin_collision_layers() {
    // Spacing 0.05 and thickness 0.005: the margin is a quarter of the triangle
    // size (0.0125), not the thickness, so contact-rich meshes rarely rebuild.
    let thickness = 0.005;
    let mut fc = flat_sheet(10, thickness);
    let mover = fc.particles.len();
    fc.particles
        .push(FabricParticle::new(Vec3::new(0.22, 0.003, 0.22), 1.0));
    let grid = TriangleGrid::build(&fc.particles, &fc.triangles, thickness);
    assert!(grid.is_indexed());

    fc.particles[mover].position.y += 0.011;
    assert!(
        !grid.has_drifted(&fc.particles, &[mover]),
        "a move below a quarter of the triangle size must be tolerated"
    );
    fc.particles[mover].position.y += 0.003;
    assert!(
        grid.has_drifted(&fc.particles, &[mover]),
        "a move beyond a quarter of the triangle size must be reported"
    );
}

#[test]
fn grid_gives_up_after_its_rebuild_budget_and_scans_everything() {
    let thickness = 0.01;
    let fc = flat_sheet(4, thickness);
    let mut grid = TriangleGrid::build(&fc.particles, &fc.triangles, thickness);
    let budget = REBUILD_BUDGET_BASE + fc.particles.len() / REBUILD_BUDGET_PARTICLES;
    assert!(budget >= 1);
    for rebuild in 1..=budget {
        assert!(
            grid.rebuild(&fc.particles, &fc.triangles, thickness),
            "rebuild {rebuild} of {budget} is within the budget"
        );
        assert!(grid.is_indexed());
    }
    assert!(
        !grid.rebuild(&fc.particles, &fc.triangles, thickness),
        "the budget is used up"
    );
    assert!(!grid.is_indexed());
    let everything: Vec<usize> = (0..fc.triangles.len()).collect();
    for vertex in 0..fc.particles.len() {
        assert_eq!(grid.candidates(vertex), everything.as_slice());
    }
    assert!(
        !grid.has_drifted(&fc.particles, &[0]),
        "nothing can drift without a grid"
    );
    // Once given up the grid stays given up.
    assert!(!grid.rebuild(&fc.particles, &fc.triangles, thickness));
    assert!(!grid.is_indexed());
}

// ---------------------------------------------------------------------------
// Degenerate input
// ---------------------------------------------------------------------------

#[test]
fn degenerate_parameters_leave_the_particles_untouched() {
    let reference = stacked_sheets(8, 3, 0.008, 0.02, 91);
    let nothing = SelfCollisionStats::default();

    for thickness in [0.0, -0.5, Real::NAN] {
        let mut fc = reference.clone();
        fc.collision_thickness = thickness;
        assert_eq!(
            fc.resolve_self_collisions_with_stats(),
            nothing,
            "thickness {thickness}"
        );
        assert_same_positions(&format!("thickness {thickness}"), &fc, &reference);
    }

    let mut no_iterations = reference.clone();
    no_iterations.collision_iterations = 0;
    assert_eq!(no_iterations.resolve_self_collisions_with_stats(), nothing);
    assert_same_positions("zero iterations", &no_iterations, &reference);

    let mut no_triangles = reference.clone();
    no_triangles.triangles.clear();
    assert_eq!(no_triangles.resolve_self_collisions_with_stats(), nothing);
    assert_same_positions("no triangles", &no_triangles, &reference);

    let mut empty = FabricCollision::new(Vec::new(), Vec::new(), 0.01, 0.3);
    assert_eq!(empty.resolve_self_collisions_with_stats(), nothing);
}

#[test]
fn infinite_thickness_falls_back_to_the_exhaustive_scan() {
    let mut expected = stacked_sheets(6, 2, 0.008, Real::INFINITY, 101);
    expected.collision_iterations = 1;
    let mut got = expected.clone();
    let grid = TriangleGrid::build(&got.particles, &got.triangles, got.collision_thickness);
    assert!(!grid.is_indexed());
    resolve_exhaustively(&mut expected);
    let stats = got.resolve_self_collisions_with_stats();
    assert_eq!(
        stats.pairs_tested,
        got.particles.len() * got.triangles.len()
    );
    assert_same_positions("infinite thickness", &got, &expected);
}

#[test]
fn non_finite_particles_fall_back_to_the_exhaustive_scan() {
    for poison in [Real::NAN, Real::INFINITY, Real::NEG_INFINITY] {
        let mut expected = stacked_sheets(8, 3, 0.008, 0.02, 111);
        // A stray particle that belongs to no triangle.
        expected
            .particles
            .push(FabricParticle::new(Vec3::new(poison, 0.3, 0.1), 1.0));
        let mut got = expected.clone();
        let grid = TriangleGrid::build(&got.particles, &got.triangles, got.collision_thickness);
        assert!(!grid.is_indexed(), "poison {poison}");

        resolve_exhaustively(&mut expected);
        let stats = got.resolve_self_collisions_with_stats();
        assert!(stats.contacts_resolved > 0, "poison {poison}: {stats:?}");
        assert_same_positions(&format!("poison {poison}"), &got, &expected);
    }
}

#[test]
fn triangle_with_a_missing_particle_is_not_indexed() {
    let particles = vec![
        FabricParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
        FabricParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
    ];
    let triangles = vec![[0, 1, 7]];
    let grid = TriangleGrid::build(&particles, &triangles, 0.01);
    assert!(!grid.is_indexed());
    assert_eq!(grid.candidates(0), &[0]);
}

#[test]
fn pinned_geometry_is_never_moved() {
    let thickness = 0.01;
    let mut fc = FabricCollision::new(
        vec![
            FabricParticle::pinned(Vec3::new(0.0, 0.0, 0.0)),
            FabricParticle::pinned(Vec3::new(1.0, 0.0, 0.0)),
            FabricParticle::pinned(Vec3::new(0.0, 0.0, 1.0)),
            FabricParticle::pinned(Vec3::new(0.25, 0.002, 0.25)),
        ],
        vec![[0, 1, 2]],
        thickness,
        0.3,
    );
    let reference = fc.clone();
    let stats = fc.resolve_self_collisions_with_stats();
    assert_eq!(stats.contacts_resolved, 0, "nothing can move: {stats:?}");
    assert!(stats.pairs_tested >= 1);
    assert_same_positions("all pinned", &fc, &reference);
}

#[test]
fn single_contact_pushes_vertex_and_triangle_apart_by_inverse_mass() {
    // Unit masses: the vertex moves by c, each corner by c / 3 the other way,
    // with c = depth / (1 + 1) for depth = thickness - distance.
    let thickness = 0.01;
    let mut fc = FabricCollision::new(
        vec![
            FabricParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            FabricParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
            FabricParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            FabricParticle::new(Vec3::new(0.25, 0.002, 0.25), 1.0),
        ],
        vec![[0, 1, 2]],
        thickness,
        0.3,
    );
    fc.collision_iterations = 1;
    let stats = fc.resolve_self_collisions_with_stats();
    assert_eq!(stats.contacts_resolved, 1);
    let depth = thickness - 0.002;
    let correction = depth / 2.0;
    let vertex_y = fc.particles[3].position.y;
    assert!(
        (vertex_y - (0.002 + correction)).abs() < 1.0e-15,
        "{vertex_y}"
    );
    for corner in 0..3 {
        let y = fc.particles[corner].position.y;
        assert!(
            (y + correction / 3.0).abs() < 1.0e-15,
            "corner {corner}: {y}"
        );
    }
}
