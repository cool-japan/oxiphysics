//! Voronoi-pole approximation of the medial axis of a surface sample.
//!
//! The Voronoi vertices of a dense sample of a surface are the circumcentres of
//! its Delaunay tetrahedra. For every sample the *pole* is the farthest Voronoi
//! vertex of its cell (Amenta and Bern, "Surface reconstruction by Voronoi
//! filtering", 1999): it lies near the medial axis, and its distance to the
//! sample approximates the radius of the maximal ball touching the surface at
//! that sample.

use super::functions::{add3, cross3, dot3, len3, scale3, sub3};
use super::types::MedialAxisPoint;
use crate::computational_geometry::delaunay_3d::{Tetrahedralization, delaunay_3d};
use crate::error::{Error, Result};

/// Circumsphere `(centre, radius)` of the tetrahedron `(a, b, c, d)`, or `None`
/// when it is not usable as a medial ball.
///
/// The centre solves `2 (p - a) . x = |p - a|^2` for `p` in `b, c, d`, which
/// Cramer's rule writes as `x = (h1 r2 x r3 + h2 r3 x r1 + h3 r1 x r2) / det`
/// with `r_i = p_i - a`, `h_i = |r_i|^2 / 2` and `det = r1 . (r2 x r3)`. A
/// nearly flat tetrahedron (`|det|` below `1e-10 |r1| |r2| |r3|`) has a
/// numerically meaningless circumsphere and is rejected, as is any non-finite
/// result.
fn circumsphere(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<([f64; 3], f64)> {
    let r1 = sub3(b, a);
    let r2 = sub3(c, a);
    let r3 = sub3(d, a);
    let r2_cross_r3 = cross3(r2, r3);
    let det = dot3(r1, r2_cross_r3);
    let scale = len3(r1) * len3(r2) * len3(r3);
    if !det.is_finite() || !scale.is_finite() || scale <= 0.0 || det.abs() <= 1e-10 * scale {
        return None;
    }
    let numerator = add3(
        add3(
            scale3(r2_cross_r3, 0.5 * dot3(r1, r1)),
            scale3(cross3(r3, r1), 0.5 * dot3(r2, r2)),
        ),
        scale3(cross3(r1, r2), 0.5 * dot3(r3, r3)),
    );
    let offset = scale3(numerator, 1.0 / det);
    let radius = len3(offset);
    if !radius.is_finite() {
        return None;
    }
    Some((add3(a, offset), radius))
}

/// Area-weighted outward normal of every vertex that lies on the boundary of
/// the tetrahedralization (`None` for interior samples).
///
/// A boundary face is a tetrahedron face without a neighbour. Its outward sense
/// is fixed against the tetrahedron's opposite vertex instead of trusting the
/// face winding.
fn boundary_vertex_normals(tetra: &Tetrahedralization) -> Vec<Option<[f64; 3]>> {
    let points = &tetra.points;
    let mut normals: Vec<Option<[f64; 3]>> = vec![None; points.len()];
    for tet in &tetra.tets {
        for (opposite, &neighbour) in tet.neigh.iter().enumerate() {
            if neighbour >= 0 {
                continue;
            }
            let face = tet.face(opposite);
            let origin = points[face[0]];
            let mut normal = cross3(sub3(points[face[1]], origin), sub3(points[face[2]], origin));
            if dot3(normal, sub3(points[tet.v[opposite]], origin)) > 0.0 {
                normal = scale3(normal, -1.0);
            }
            for &vertex in &face {
                normals[vertex] = Some(add3(normals[vertex].unwrap_or([0.0; 3]), normal));
            }
        }
    }
    normals
}

/// Diagonal length of the axis-aligned bounding box of `points`.
fn bounding_box_diagonal(points: &[[f64; 3]]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    len3(sub3(hi, lo))
}

/// Compute an approximate medial axis of a surface sample from its Voronoi
/// poles.
///
/// The surface samples are tetrahedralized with the exact 3-D Delaunay
/// triangulation; the circumcentre of every tetrahedron is a Voronoi vertex.
/// For each sample the pole is the farthest Voronoi vertex among its incident
/// tetrahedra, and the returned [`MedialAxisPoint`] carries
///
/// * `position`: the pole (a medial-axis candidate),
/// * `radius`: the distance from the pole to the sample, i.e. the radius of the
///   empty ball centred at the pole that touches the sample,
/// * `closest_verts`: the sample itself and the other vertex of the pole's
///   tetrahedron nearest to the pole (all four vertices are equidistant up to
///   rounding).
///
/// Samples on the boundary of the tetrahedralization (the convex hull) have an
/// unbounded Voronoi cell whose farthest finite vertex can lie outside the
/// shape, so their pole is restricted to Voronoi vertices on the inner side of
/// the sample's outward hull normal. Circumspheres of nearly flat tetrahedra
/// and spheres larger than the bounding-box diagonal cannot be medial balls
/// and are ignored; a sample without any usable incident tetrahedron (or an
/// exact duplicate of an earlier sample) yields no point.
///
/// At most `num_candidates` points are returned; when more samples have a pole
/// they are thinned to an evenly spaced subset (in sample order) so the cap
/// does not bias the result towards the start of the input. The output is
/// deterministic.
///
/// # Errors
///
/// * [`Error::TooFewPoints`] for fewer than four vertices,
/// * [`Error::DegenerateGeometry`] for non-finite coordinates or when all
///   vertices are coplanar (no 3-D tetrahedralization, hence no 3-D medial
///   axis, exists).
pub fn compute_medial_axis_approx(
    surface_vertices: &[[f64; 3]],
    num_candidates: usize,
) -> Result<Vec<MedialAxisPoint>> {
    let n = surface_vertices.len();
    if n < 4 {
        return Err(Error::too_few_points(4, n));
    }
    if surface_vertices
        .iter()
        .any(|p| p.iter().any(|x| !x.is_finite()))
    {
        return Err(Error::degenerate_geometry(
            "medial axis: surface vertex has a non-finite coordinate",
        ));
    }
    let tetra = delaunay_3d(surface_vertices).ok_or_else(|| {
        Error::degenerate_geometry(
            "medial axis: surface vertices are coplanar, so no 3-D Delaunay tetrahedralization exists",
        )
    })?;
    let diagonal = bounding_box_diagonal(surface_vertices);
    let spheres: Vec<Option<([f64; 3], f64)>> = tetra
        .tets
        .iter()
        .map(|t| {
            circumsphere(
                surface_vertices[t.v[0]],
                surface_vertices[t.v[1]],
                surface_vertices[t.v[2]],
                surface_vertices[t.v[3]],
            )
            .filter(|&(_, radius)| radius <= diagonal)
        })
        .collect();
    let hull_normals = boundary_vertex_normals(&tetra);

    // Pole of every sample: the incident tet with the largest usable circumsphere.
    let mut pole_tet: Vec<Option<(usize, f64)>> = vec![None; n];
    for (ti, tet) in tetra.tets.iter().enumerate() {
        let Some((centre, radius)) = spheres[ti] else {
            continue;
        };
        for &v in &tet.v {
            if let Some(normal) = hull_normals[v]
                && dot3(sub3(centre, surface_vertices[v]), normal) >= 0.0
            {
                continue;
            }
            if pole_tet[v].is_none_or(|(_, best)| radius > best) {
                pole_tet[v] = Some((ti, radius));
            }
        }
    }

    let mut points = Vec::new();
    for (sample, pole) in pole_tet.iter().enumerate() {
        let Some(&(ti, _)) = pole.as_ref() else {
            continue;
        };
        let Some((centre, _)) = spheres[ti] else {
            continue;
        };
        let distance_to = |vertex: usize| len3(sub3(surface_vertices[vertex], centre));
        let mut partner: Option<(usize, f64)> = None;
        for &vertex in &tetra.tets[ti].v {
            if vertex == sample {
                continue;
            }
            let d = distance_to(vertex);
            if partner.is_none_or(|(best_vertex, best_d)| {
                d < best_d || (d == best_d && vertex < best_vertex)
            }) {
                partner = Some((vertex, d));
            }
        }
        let Some((partner_vertex, _)) = partner else {
            continue;
        };
        points.push(MedialAxisPoint {
            position: centre,
            radius: distance_to(sample),
            closest_verts: [sample, partner_vertex],
        });
    }

    if points.len() > num_candidates {
        let total = points.len();
        points = (0..num_candidates)
            .map(|k| points[k * total / num_candidates].clone())
            .collect();
    }
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `n` points spread over a sphere (Fibonacci lattice).
    fn fibonacci_sphere(n: usize, radius: f64, centre: [f64; 3]) -> Vec<[f64; 3]> {
        let golden_angle = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());
        (0..n)
            .map(|i| {
                let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
                let ring = (1.0 - z * z).sqrt();
                let theta = golden_angle * i as f64;
                [
                    centre[0] + radius * ring * theta.cos(),
                    centre[1] + radius * ring * theta.sin(),
                    centre[2] + radius * z,
                ]
            })
            .collect()
    }

    fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
        len3(sub3(a, b))
    }

    #[test]
    fn sphere_poles_sit_at_the_centre_with_the_sphere_radius() {
        let samples = fibonacci_sphere(200, 2.0, [0.0; 3]);
        let medial = compute_medial_axis_approx(&samples, 1000).expect("sphere samples");
        assert_eq!(medial.len(), 200, "one pole per sample");
        for (i, m) in medial.iter().enumerate() {
            assert!(
                len3(m.position) < 0.1,
                "pole {i} at {:?} is not at the centre",
                m.position
            );
            assert!(
                (m.radius - 2.0).abs() < 0.1,
                "pole {i} has radius {}",
                m.radius
            );
        }
    }

    #[test]
    fn translated_sphere_poles_follow_the_centre() {
        let centre = [1.0, -2.0, 3.0];
        let samples = fibonacci_sphere(150, 1.5, centre);
        let medial = compute_medial_axis_approx(&samples, 1000).expect("sphere samples");
        assert_eq!(medial.len(), 150);
        for m in &medial {
            assert!(distance(m.position, centre) < 0.15, "pole {:?}", m.position);
            assert!((m.radius - 1.5).abs() < 0.15, "radius {}", m.radius);
        }
    }

    #[test]
    fn prolate_ellipsoid_poles_lie_on_the_focal_segment() {
        // Semi-axes (3, 1, 1): the medial axis is the segment on the x axis with
        // |x| <= (a^2 - b^2) / a = 8/3, and balls shrink from radius 1 to 1/3.
        let samples: Vec<[f64; 3]> = fibonacci_sphere(600, 1.0, [0.0; 3])
            .into_iter()
            .map(|p| [3.0 * p[0], p[1], p[2]])
            .collect();
        let medial = compute_medial_axis_approx(&samples, 10_000).expect("ellipsoid samples");
        assert!(medial.len() >= 590, "got {} poles", medial.len());
        for m in &medial {
            let off_axis = m.position[1].hypot(m.position[2]);
            assert!(
                off_axis < 0.1,
                "pole {:?} is {off_axis} off the axis",
                m.position
            );
            assert!(
                m.position[0].abs() < 2.75,
                "pole {:?} beyond the segment",
                m.position
            );
            assert!(
                m.radius > 0.25 && m.radius < 1.15,
                "pole {:?} has radius {}",
                m.position,
                m.radius
            );
        }
    }

    #[test]
    fn cube_corner_poles_sit_at_the_cube_centre() {
        let mut corners = Vec::new();
        for x in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                for z in [0.0, 1.0] {
                    corners.push([x, y, z]);
                }
            }
        }
        let medial = compute_medial_axis_approx(&corners, 100).expect("cube corners");
        assert_eq!(medial.len(), 8);
        for m in &medial {
            assert!(
                distance(m.position, [0.5; 3]) < 1e-12,
                "pole {:?}",
                m.position
            );
            assert!((m.radius - 0.75f64.sqrt()).abs() < 1e-12);
        }
    }

    #[test]
    fn closest_vertices_are_equidistant_from_the_pole() {
        let samples = fibonacci_sphere(120, 2.0, [0.5, 0.0, -1.0]);
        let medial = compute_medial_axis_approx(&samples, 1000).expect("sphere samples");
        for m in &medial {
            let [first, second] = m.closest_verts;
            assert_ne!(first, second);
            let d0 = distance(m.position, samples[first]);
            let d1 = distance(m.position, samples[second]);
            assert!((d0 - m.radius).abs() < 1e-9 * m.radius.max(1.0));
            assert!(
                (d1 - m.radius).abs() < 1e-6 * m.radius.max(1.0),
                "d1 = {d1}"
            );
        }
    }

    #[test]
    fn coplanar_samples_are_rejected() {
        let circle: Vec<[f64; 3]> = (0..16)
            .map(|i| {
                let t = i as f64 / 16.0 * std::f64::consts::TAU;
                [t.cos(), t.sin(), 0.0]
            })
            .collect();
        let err = compute_medial_axis_approx(&circle, 4).expect_err("coplanar input");
        assert!(err.is_degenerate(), "unexpected error: {err}");
    }

    #[test]
    fn too_few_or_non_finite_samples_are_rejected() {
        let three = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(matches!(
            compute_medial_axis_approx(&three, 4),
            Err(Error::TooFewPoints {
                required: 4,
                actual: 3
            })
        ));
        let mut samples = fibonacci_sphere(20, 1.0, [0.0; 3]);
        samples[7][1] = f64::NAN;
        let err = compute_medial_axis_approx(&samples, 4).expect_err("NaN coordinate");
        assert!(err.is_degenerate(), "unexpected error: {err}");
    }

    #[test]
    fn candidate_cap_thins_the_output_evenly() {
        let samples = fibonacci_sphere(200, 2.0, [0.0; 3]);
        let capped = compute_medial_axis_approx(&samples, 10).expect("sphere samples");
        assert_eq!(capped.len(), 10);
        let picked: Vec<usize> = capped.iter().map(|m| m.closest_verts[0]).collect();
        assert_eq!(picked, (0..10).map(|k| k * 20).collect::<Vec<_>>());
        assert!(
            compute_medial_axis_approx(&samples, 0)
                .expect("zero cap")
                .is_empty()
        );
        assert_eq!(
            compute_medial_axis_approx(&samples, 5000)
                .expect("large cap")
                .len(),
            200
        );
    }

    #[test]
    fn medial_axis_output_is_deterministic() {
        let samples = fibonacci_sphere(100, 1.0, [0.0; 3]);
        let first = compute_medial_axis_approx(&samples, 100).expect("first run");
        let second = compute_medial_axis_approx(&samples, 100).expect("second run");
        assert_eq!(first.len(), second.len());
        for (a, b) in first.iter().zip(&second) {
            assert_eq!(a.position, b.position);
            assert_eq!(a.radius, b.radius);
            assert_eq!(a.closest_verts, b.closest_verts);
        }
    }
}
