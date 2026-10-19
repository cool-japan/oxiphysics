// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bounded Voronoi cells.
//!
//! The Voronoi cell of site `s_i` is the intersection of the half-planes
//! `{x : |x - s_i| <= |x - s_j|}` over all other sites `s_j`.  The cells of
//! the sites on the convex hull of the input are unbounded, so a polygon made
//! of circumcentres alone cannot describe them.  The cells built here start
//! from the clipping rectangle and are cut by those bisector half-planes, so
//! they are always closed polygons and tile the rectangle exactly.

use super::functions::{bowyer_watson, voronoi_neighbors};

/// Relative tolerance (in units of the clipping rectangle's larger side) used
/// to decide whether a vertex lies on the wrong side of a bisector and to
/// merge coincident polygon vertices.
const RELATIVE_TOLERANCE: f64 = 1e-12;

/// Clip a convex polygon to the half-plane `{x : (x - p) . n <= 0}`
/// (Sutherland-Hodgman).
///
/// The vertices of `poly` keep their orientation.  Vertices on the boundary
/// of the half-plane count as inside.  The result is empty when the polygon
/// lies completely outside.
pub(super) fn clip_polygon_halfplane(poly: &[[f64; 2]], p: [f64; 2], n: [f64; 2]) -> Vec<[f64; 2]> {
    let count = poly.len();
    let mut output: Vec<[f64; 2]> = Vec::with_capacity(count + 1);
    for k in 0..count {
        let cur = poly[k];
        let nxt = poly[(k + 1) % count];
        let d_cur = (cur[0] - p[0]) * n[0] + (cur[1] - p[1]) * n[1];
        let d_nxt = (nxt[0] - p[0]) * n[0] + (nxt[1] - p[1]) * n[1];
        let cur_inside = d_cur <= 0.0;
        let nxt_inside = d_nxt <= 0.0;
        if cur_inside {
            output.push(cur);
        }
        if cur_inside != nxt_inside {
            // The classes differ, so d_cur != d_nxt and the division is safe.
            let t = d_cur / (d_cur - d_nxt);
            output.push([
                cur[0] + t * (nxt[0] - cur[0]),
                cur[1] + t * (nxt[1] - cur[1]),
            ]);
        }
    }
    output
}

/// Remove consecutive (cyclically) coincident vertices.
fn remove_duplicate_vertices(poly: &mut Vec<[f64; 2]>, tolerance: f64) {
    let near = |a: [f64; 2], b: [f64; 2]| {
        (a[0] - b[0]).abs() <= tolerance && (a[1] - b[1]).abs() <= tolerance
    };
    poly.dedup_by(|later, earlier| near(*earlier, *later));
    while poly.len() > 1 && near(poly[0], poly[poly.len() - 1]) {
        poly.pop();
    }
}

/// Whether any vertex of `poly` lies strictly on the far side of the plane
/// `(x - p) . n = 0`, by more than `tolerance` (in units of `n`'s length).
fn has_vertex_beyond(poly: &[[f64; 2]], p: [f64; 2], n: [f64; 2], tolerance: f64) -> bool {
    let n_len = (n[0] * n[0] + n[1] * n[1]).sqrt();
    let limit = tolerance * n_len;
    poly.iter()
        .any(|v| (v[0] - p[0]) * n[0] + (v[1] - p[1]) * n[1] > limit)
}

/// Voronoi cell of `sites[i]` inside `rect`, a counter-clockwise rectangle.
///
/// `hints` are sites expected to be Voronoi neighbours of `sites[i]`; clipping
/// by them first shrinks the cell quickly.  The result does not depend on the
/// hints: afterwards every other site whose bisector can still reach the cell
/// is tested, so a wrong or missing hint never changes the cell.
fn bounded_cell(
    sites: &[[f64; 2]],
    i: usize,
    rect: &[[f64; 2]],
    hints: &[usize],
    tolerance: f64,
) -> Vec<[f64; 2]> {
    let si = sites[i];
    let mut poly = rect.to_vec();

    // Clip by the bisector of (s_i, s_j).  `None` means s_j owns the whole of
    // s_i's region because both sites coincide and j has the lower index.
    let clip = |poly: &mut Vec<[f64; 2]>, j: usize| -> Option<()> {
        let sj = sites[j];
        let n = [sj[0] - si[0], sj[1] - si[1]];
        if n[0] == 0.0 && n[1] == 0.0 {
            // Coincident sites share one region: the lower index keeps it.
            return if j < i { None } else { Some(()) };
        }
        let mid = [(si[0] + sj[0]) * 0.5, (si[1] + sj[1]) * 0.5];
        if has_vertex_beyond(poly, mid, n, tolerance) {
            let mut clipped = clip_polygon_halfplane(poly, mid, n);
            remove_duplicate_vertices(&mut clipped, tolerance);
            *poly = clipped;
        }
        Some(())
    };

    for &j in hints {
        if j == i || j >= sites.len() {
            continue;
        }
        if clip(&mut poly, j).is_none() || poly.is_empty() {
            return Vec::new();
        }
    }

    if poly.is_empty() {
        return poly;
    }

    // Exact verification: a site at distance d from s_i has its bisector at
    // distance d / 2, so it can only cut a cell contained in the disc of
    // radius `radius` around s_i when d < 2 * radius.  Clipping only shrinks
    // the cell, which makes the skip test valid for the rest of the pass.
    let mut radius_sq = max_distance_sq(&poly, si);
    for (j, &sj) in sites.iter().enumerate() {
        if j == i {
            continue;
        }
        let dx = sj[0] - si[0];
        let dy = sj[1] - si[1];
        let d_sq = dx * dx + dy * dy;
        if d_sq > 4.0 * radius_sq * (1.0 + 1e-9) {
            continue;
        }
        if clip(&mut poly, j).is_none() || poly.is_empty() {
            return Vec::new();
        }
        radius_sq = max_distance_sq(&poly, si);
    }
    poly
}

/// Largest squared distance from `centre` to a vertex of `poly`.
fn max_distance_sq(poly: &[[f64; 2]], centre: [f64; 2]) -> f64 {
    poly.iter()
        .map(|v| (v[0] - centre[0]).powi(2) + (v[1] - centre[1]).powi(2))
        .fold(0.0, f64::max)
}

/// Voronoi cells of `sites`, clipped to the rectangle
/// `bounds = [xmin, xmax, ymin, ymax]`.
///
/// Entry `i` is the polygon (counter-clockwise vertices) of the part of the
/// rectangle that is closer to `sites[i]` than to any other site.  Unlike the
/// circumcentre polygons of [`delaunay_to_voronoi`](super::delaunay_to_voronoi),
/// the cells of sites on the convex hull of the input are closed by the
/// rectangle, so the cells of all sites tile the rectangle: their areas add up
/// to the area of the rectangle.
///
/// Every cell is the rectangle cut by the bisector half-plane of each other
/// site.  The Delaunay neighbours are used as a first guess to cut quickly,
/// and a distance test over all sites then guarantees the exact result even
/// for collinear sites, sites without Delaunay triangles, and triangles lost
/// to the finite super-triangle of the triangulation.
///
/// Sites that coincide share one region, which goes to the lowest index; the
/// others get an empty cell.  A site with a non-finite coordinate has an empty
/// cell and does not constrain the others.  When the rectangle is empty or
/// non-finite, all cells are empty.
pub fn bounded_voronoi_cells(sites: &[[f64; 2]], bounds: [f64; 4]) -> Vec<Vec<[f64; 2]>> {
    let count = sites.len();
    let [xmin, xmax, ymin, ymax] = bounds;
    let valid_bounds = bounds.iter().all(|b| b.is_finite()) && xmax > xmin && ymax > ymin;
    if count == 0 || !valid_bounds {
        return vec![Vec::new(); count];
    }
    let finite: Vec<usize> = (0..count)
        .filter(|&i| sites[i][0].is_finite() && sites[i][1].is_finite())
        .collect();
    let finite_sites: Vec<[f64; 2]> = finite.iter().map(|&i| sites[i]).collect();
    let triangles = bowyer_watson(&finite_sites);
    let adjacency = voronoi_neighbors(finite_sites.len(), &triangles);

    let rect = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]];
    let tolerance = RELATIVE_TOLERANCE * (xmax - xmin).max(ymax - ymin);
    let mut cells: Vec<Vec<[f64; 2]>> = vec![Vec::new(); count];
    for (k, &original) in finite.iter().enumerate() {
        cells[original] = bounded_cell(&finite_sites, k, &rect, &adjacency[k], tolerance);
    }
    cells
}

#[cfg(test)]
mod tests {
    use super::super::{VoronoiDiagram, lloyd_relaxation};
    use super::*;

    /// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
    fn lcg(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn polygon_area(poly: &[[f64; 2]]) -> f64 {
        let n = poly.len();
        if n < 3 {
            return 0.0;
        }
        let mut twice = 0.0;
        for k in 0..n {
            let a = poly[k];
            let b = poly[(k + 1) % n];
            twice += a[0] * b[1] - b[0] * a[1];
        }
        twice * 0.5
    }

    fn contains(poly: &[[f64; 2]], p: [f64; 2], tolerance: f64) -> bool {
        let n = poly.len();
        if n < 3 {
            return false;
        }
        (0..n).all(|k| {
            let a = poly[k];
            let b = poly[(k + 1) % n];
            let cross = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
            cross >= -tolerance
        })
    }

    fn nearest(sites: &[[f64; 2]], p: [f64; 2]) -> usize {
        let mut best = 0;
        let mut best_d = f64::INFINITY;
        for (i, s) in sites.iter().enumerate() {
            let d = (s[0] - p[0]).powi(2) + (s[1] - p[1]).powi(2);
            if d < best_d {
                best_d = d;
                best = i;
            }
        }
        best
    }

    fn total_area(cells: &[Vec<[f64; 2]>]) -> f64 {
        cells.iter().map(|c| polygon_area(c)).sum()
    }

    #[test]
    fn clip_halfplane_cuts_unit_square_in_half() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        // Keep x <= 0.5.
        let clipped = clip_polygon_halfplane(&square, [0.5, 0.0], [1.0, 0.0]);
        assert_eq!(clipped.len(), 4);
        assert!((polygon_area(&clipped) - 0.5).abs() < 1e-15);
        assert!(clipped.iter().all(|v| v[0] <= 0.5 + 1e-15));
    }

    #[test]
    fn clip_halfplane_keeps_or_discards_whole_polygon() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let kept = clip_polygon_halfplane(&square, [2.0, 0.0], [1.0, 0.0]);
        assert_eq!(kept.len(), 4);
        let dropped = clip_polygon_halfplane(&square, [-1.0, 0.0], [1.0, 0.0]);
        assert!(dropped.is_empty());
    }

    #[test]
    fn clip_halfplane_diagonal_cut_makes_triangle() {
        let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        // Keep x + y <= 1: the lower-left triangle of area 1/2.
        let clipped = clip_polygon_halfplane(&square, [0.5, 0.5], [1.0, 1.0]);
        assert!((polygon_area(&clipped) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn twenty_random_sites_tile_the_bounding_box() {
        let mut next = lcg(0x5eed_0001);
        let bounds = [0.0, 3.0, -1.0, 1.0];
        let sites: Vec<[f64; 2]> = (0..20)
            .map(|_| [3.0 * next(), -1.0 + 2.0 * next()])
            .collect();
        let cells = bounded_voronoi_cells(&sites, bounds);
        assert_eq!(cells.len(), 20);
        let bbox_area = 3.0 * 2.0;
        let sum = total_area(&cells);
        assert!(
            ((sum - bbox_area) / bbox_area).abs() < 1e-9,
            "cell areas {sum} must add up to the bounding box area {bbox_area}"
        );
        for (i, cell) in cells.iter().enumerate() {
            assert!(cell.len() >= 3, "cell {i} must be a polygon");
            assert!(
                polygon_area(cell) > 0.0,
                "cell {i} must be counter-clockwise"
            );
            assert!(
                contains(cell, sites[i], 1e-12),
                "cell {i} must contain its site"
            );
        }
    }

    #[test]
    fn every_sample_point_lies_in_the_cell_of_its_nearest_site() {
        let mut next = lcg(0xfeed_beef);
        let bounds = [-2.0, 2.0, -2.0, 2.0];
        let sites: Vec<[f64; 2]> = (0..25)
            .map(|_| [-2.0 + 4.0 * next(), -2.0 + 4.0 * next()])
            .collect();
        let cells = bounded_voronoi_cells(&sites, bounds);
        let mut checked = 0;
        for gx in 0..40 {
            for gy in 0..40 {
                let p = [
                    -2.0 + (gx as f64 + 0.5) * 0.1,
                    -2.0 + (gy as f64 + 0.5) * 0.1,
                ];
                let owner = nearest(&sites, p);
                assert!(
                    contains(&cells[owner], p, 1e-9),
                    "point {p:?} must lie in the cell of its nearest site {owner}"
                );
                for (j, cell) in cells.iter().enumerate() {
                    if j != owner {
                        assert!(
                            !contains(cell, p, -1e-9),
                            "point {p:?} must lie strictly outside the cell of site {j}"
                        );
                    }
                }
                checked += 1;
            }
        }
        assert_eq!(checked, 1600);
    }

    #[test]
    fn hull_sites_get_closed_quadrant_cells() {
        // Four sites in convex position: every cell is unbounded, and the
        // circumcentre polygons of the Delaunay triangulation degenerate to a
        // single point.  The bounded cells are the four quadrants.
        let sites = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]];
        let cells = bounded_voronoi_cells(&sites, [0.0, 1.0, 0.0, 1.0]);
        for (i, cell) in cells.iter().enumerate() {
            assert_eq!(cell.len(), 4, "quadrant cell {i} has four corners");
            assert!((polygon_area(cell) - 0.25).abs() < 1e-12);
        }
        assert!((total_area(&cells) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn collinear_sites_get_parallel_strips() {
        // No Delaunay triangle exists; the verification pass must still find
        // the exact cells.
        let sites = [[0.1, 0.5], [0.4, 0.5], [0.9, 0.5]];
        let cells = bounded_voronoi_cells(&sites, [0.0, 1.0, 0.0, 1.0]);
        let expected = [0.25, 0.4, 0.35];
        for (cell, &area) in cells.iter().zip(expected.iter()) {
            assert!((polygon_area(cell) - area).abs() < 1e-12);
        }
    }

    #[test]
    fn single_site_owns_the_whole_rectangle() {
        let cells = bounded_voronoi_cells(&[[0.3, 0.3]], [0.0, 2.0, 0.0, 1.0]);
        assert_eq!(cells.len(), 1);
        assert!((polygon_area(&cells[0]) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn coincident_sites_share_one_region() {
        let sites = [[0.25, 0.5], [0.75, 0.5], [0.25, 0.5]];
        let cells = bounded_voronoi_cells(&sites, [0.0, 1.0, 0.0, 1.0]);
        assert!((polygon_area(&cells[0]) - 0.5).abs() < 1e-12);
        assert!((polygon_area(&cells[1]) - 0.5).abs() < 1e-12);
        assert!(
            cells[2].is_empty(),
            "the higher duplicate gets an empty cell"
        );
        assert!((total_area(&cells) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sites_outside_the_rectangle_still_tile_it() {
        let sites = [[-1.0, 0.5], [0.5, 0.5], [3.0, 0.5], [0.5, 4.0]];
        let cells = bounded_voronoi_cells(&sites, [0.0, 1.0, 0.0, 1.0]);
        assert!((total_area(&cells) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn non_finite_site_gets_an_empty_cell_and_constrains_nobody() {
        let sites = [[0.25, 0.5], [f64::NAN, 0.5], [0.75, 0.5]];
        let cells = bounded_voronoi_cells(&sites, [0.0, 1.0, 0.0, 1.0]);
        assert!(cells[1].is_empty());
        assert!((polygon_area(&cells[0]) - 0.5).abs() < 1e-12);
        assert!((polygon_area(&cells[2]) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn empty_input_and_empty_rectangle() {
        assert!(bounded_voronoi_cells(&[], [0.0, 1.0, 0.0, 1.0]).is_empty());
        let cells = bounded_voronoi_cells(&[[0.5, 0.5]], [1.0, 1.0, 0.0, 1.0]);
        assert_eq!(cells.len(), 1);
        assert!(cells[0].is_empty());
        let cells = bounded_voronoi_cells(&[[0.5, 0.5]], [0.0, f64::INFINITY, 0.0, 1.0]);
        assert!(cells[0].is_empty());
    }

    #[test]
    fn regular_lattice_of_sites_gets_unit_square_cells() {
        // On a lattice every four neighbouring sites are cocircular, so the
        // Delaunay triangulation is not unique; the cells must still be the
        // unit squares around the sites.
        let sites: Vec<[f64; 2]> = (0..4)
            .flat_map(|i| (0..4).map(move |j| [i as f64 + 0.5, j as f64 + 0.5]))
            .collect();
        let cells = bounded_voronoi_cells(&sites, [0.0, 4.0, 0.0, 4.0]);
        for (i, cell) in cells.iter().enumerate() {
            assert_eq!(cell.len(), 4, "cell {i} is a square");
            assert!((polygon_area(cell) - 1.0).abs() < 1e-12, "cell {i}");
            assert!(contains(cell, sites[i], 1e-12));
        }
    }

    #[test]
    fn voronoi_diagram_cell_areas_add_up_to_the_bounding_box() {
        let mut next = lcg(0x0d1a_9a77);
        let bbox = [-1.0, 3.0, 0.5, 2.5];
        let pts: Vec<[f64; 2]> = (0..20)
            .map(|_| [-1.0 + 4.0 * next(), 0.5 + 2.0 * next()])
            .collect();
        let bbox_area = 4.0 * 2.0;
        let diagram = VoronoiDiagram::from_sites(&pts, bbox);
        let total: f64 = (0..pts.len()).map(|i| diagram.cell_area(i)).sum();
        assert!(
            ((total - bbox_area) / bbox_area).abs() < 1e-9,
            "cell areas {total} must add up to {bbox_area}"
        );
        assert!((0..pts.len()).all(|i| diagram.cell_area(i) > 0.0));
        // The diagram of the relaxed sites tiles the box as well.
        let relaxed = diagram.lloyd_relaxation(3);
        let total: f64 = (0..pts.len()).map(|i| relaxed.cell_area(i)).sum();
        assert!(((total - bbox_area) / bbox_area).abs() < 1e-9);
    }

    #[test]
    fn lloyd_step_moves_hull_sites_to_their_quadrant_centroids() {
        // All four sites are on the convex hull of the input, so the polygon of
        // circumcentres of such a site is not its cell and a relaxation built
        // on it left them in place.  Their bounded cells are the four quadrants.
        let bounds = [0.0, 1.0, 0.0, 1.0];
        let sites = [[0.1, 0.1], [0.9, 0.1], [0.9, 0.9], [0.1, 0.9]];
        let relaxed = lloyd_relaxation(&sites, bounds, 1);
        let expected = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]];
        for (site, want) in relaxed.iter().zip(expected) {
            assert!(
                (site[0] - want[0]).abs() < 1e-12 && (site[1] - want[1]).abs() < 1e-12,
                "got {relaxed:?}"
            );
        }
        // The quadrant centroids are a fixed point of the relaxation.
        let again = lloyd_relaxation(&relaxed, bounds, 1);
        for (site, want) in again.iter().zip(expected) {
            assert!((site[0] - want[0]).abs() < 1e-12 && (site[1] - want[1]).abs() < 1e-12);
        }
    }

    /// Integral of `|x - site|^2` over a convex polygon: the exact triangle
    /// formula `A / 12 (sum |v|^2 + |sum v|^2)` on a fan, with the vertices taken
    /// relative to `site`.
    fn quantisation_energy(poly: &[[f64; 2]], site: [f64; 2]) -> f64 {
        let n = poly.len();
        if n < 3 {
            return 0.0;
        }
        let rel = |k: usize| [poly[k][0] - site[0], poly[k][1] - site[1]];
        let first = rel(0);
        let mut energy = 0.0;
        for k in 1..n - 1 {
            let (b, c) = (rel(k), rel(k + 1));
            let area = 0.5
                * ((b[0] - first[0]) * (c[1] - first[1]) - (b[1] - first[1]) * (c[0] - first[0]));
            let sum = [first[0] + b[0] + c[0], first[1] + b[1] + c[1]];
            let squares = first[0] * first[0]
                + first[1] * first[1]
                + b[0] * b[0]
                + b[1] * b[1]
                + c[0] * c[0]
                + c[1] * c[1];
            energy += area.abs() / 12.0 * (squares + sum[0] * sum[0] + sum[1] * sum[1]);
        }
        energy
    }

    #[test]
    fn lloyd_relaxation_never_increases_the_quantisation_energy() {
        // Moving every site to the centroid of its cell minimises the energy
        // for that partition, and re-partitioning lowers it further.  The
        // property needs the true bounded cells of the hull sites as well.
        let mut next = lcg(0x11c0_ffee);
        let bounds = [0.0, 2.0, 0.0, 1.0];
        let mut sites: Vec<[f64; 2]> = (0..30).map(|_| [2.0 * next(), next()]).collect();
        let energy = |sites: &[[f64; 2]]| -> f64 {
            bounded_voronoi_cells(sites, bounds)
                .iter()
                .zip(sites)
                .map(|(cell, &site)| quantisation_energy(cell, site))
                .sum()
        };
        let mut previous = energy(&sites);
        for step in 0..8 {
            sites = lloyd_relaxation(&sites, bounds, 1);
            let current = energy(&sites);
            assert!(
                current <= previous + 1e-12,
                "step {step}: energy rose from {previous} to {current}"
            );
            previous = current;
        }
    }

    #[test]
    fn quantisation_energy_of_a_square_matches_the_closed_form() {
        // Integral of x^2 + y^2 over the unit square centred on the site: 1/6.
        let square = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]];
        assert!((quantisation_energy(&square, [0.0, 0.0]) - 1.0 / 6.0).abs() < 1e-15);
    }
}
