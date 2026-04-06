//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
use super::super::gjk::{Simplex, SupportPoint};
#[allow(unused_imports)]
use super::functions_2::*;
use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

#[allow(unused_imports)]
use super::functions::*;

/// A face of the EPA polytope (raw array version).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EpaFaceRaw {
    /// Indices of the three vertices forming this face.
    pub vertices: [usize; 3],
    /// Outward-pointing normal (unit vector).
    pub normal: [f64; 3],
    /// Distance from the origin to this face along the normal.
    pub distance: f64,
}
/// Stateful EPA solver providing minimum-penetration-axis and polytope-refinement
/// operations as explicit, named entry points.
///
/// The underlying algorithm is the same expanding-polytope approach used by
/// `epa_penetration`, but wrapped in a struct to allow step-by-step refinement
/// and integration testing.
#[allow(dead_code)]
pub struct EpaSolver {
    /// Current polytope being refined.
    pub(super) polytope: Option<EpaPolytope>,
    /// Number of iterations performed so far.
    pub(super) iterations: usize,
    /// Last convergence tolerance used.
    pub(super) tolerance: f64,
}
#[allow(dead_code)]
impl EpaSolver {
    /// Create a new solver seeded from a GJK tetrahedron simplex.
    pub fn from_simplex(simplex: &[[f64; 3]; 4]) -> Self {
        Self {
            polytope: Some(EpaPolytope::from_gjk_simplex(simplex)),
            iterations: 0,
            tolerance: 1e-6,
        }
    }
    /// Override the convergence tolerance (default 1e-6).
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }
    /// Compute the *minimum penetration axis* – the contact normal and depth
    /// for the current polytope state without expanding further.
    ///
    /// Returns `(normal, depth)` where `normal` is a unit vector pointing from
    /// shape B to shape A and `depth` is the penetration depth.
    ///
    /// If the polytope is empty or unset, returns `([0,1,0], 0.0)`.
    pub fn compute_penetration_axis(&self) -> ([f64; 3], f64) {
        let poly = match &self.polytope {
            Some(p) => p,
            None => return ([0.0, 1.0, 0.0], 0.0),
        };
        if poly.faces.is_empty() {
            return ([0.0, 1.0, 0.0], 0.0);
        }
        let (_idx, dist, normal) = poly.closest_face();
        (normal, dist)
    }
    /// Perform a single polytope refinement step.
    ///
    /// `support_fn(dir)` should return the support point of the Minkowski
    /// difference in direction `dir`.
    ///
    /// Returns `true` if the polytope was successfully expanded (i.e. the
    /// algorithm has not yet converged), or `false` if it converged or failed.
    pub fn refine_polytope(&mut self, support_fn: &mut dyn FnMut([f64; 3]) -> [f64; 3]) -> bool {
        let poly = match self.polytope.as_mut() {
            Some(p) => p,
            None => return false,
        };
        if poly.faces.is_empty() {
            return false;
        }
        self.iterations += 1;
        let closest_idx = poly.find_closest_face();
        let closest = poly.faces[closest_idx].clone();
        let new_point = support_fn(closest.normal);
        let new_dist = epa_dot3(new_point, closest.normal);
        if (new_dist - closest.distance).abs() < self.tolerance {
            return false;
        }
        poly.expand(new_point)
    }
    /// Run the solver to completion, performing at most `max_iter` refinement
    /// steps and returning the final penetration result.
    ///
    /// Equivalent to calling `refine_polytope` in a loop and then
    /// `compute_penetration_axis`.
    pub fn solve(
        &mut self,
        support_fn: &mut dyn FnMut([f64; 3]) -> [f64; 3],
        max_iter: usize,
    ) -> Option<EpaPenetration> {
        for _ in 0..max_iter {
            if !self.refine_polytope(support_fn) {
                break;
            }
        }
        let (normal, depth) = self.compute_penetration_axis();
        if depth <= 0.0 && normal == [0.0, 1.0, 0.0] {
            return None;
        }
        let poly = self.polytope.as_ref()?;
        if poly.faces.is_empty() {
            return None;
        }
        let closest_idx = poly.find_closest_face();
        let closest = &poly.faces[closest_idx];
        let witness_a = poly.vertices[closest.vertices[0]];
        let witness_b = epa_sub3(witness_a, epa_scale3(normal, depth));
        Some(EpaPenetration {
            normal,
            depth,
            witness_a,
            witness_b,
        })
    }
    /// Number of refinement steps performed so far.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
    /// Current number of faces in the polytope.
    pub fn face_count(&self) -> usize {
        self.polytope.as_ref().map(|p| p.faces.len()).unwrap_or(0)
    }
    /// Current number of vertices in the polytope.
    pub fn vertex_count(&self) -> usize {
        self.polytope
            .as_ref()
            .map(|p| p.vertices.len())
            .unwrap_or(0)
    }
}
/// EPA polytope storing vertices and faces.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EpaPolytope {
    /// Vertices of the polytope (Minkowski-difference points).
    pub vertices: Vec<[f64; 3]>,
    /// Faces of the polytope.
    pub faces: Vec<EpaFaceRaw>,
}
#[allow(dead_code)]
impl EpaPolytope {
    /// Initialize the polytope from a GJK tetrahedron (4 vertices).
    pub fn from_gjk_simplex(simplex: &[[f64; 3]; 4]) -> Self {
        let vertices = simplex.to_vec();
        let face_indices: [[usize; 3]; 4] = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
        let mut faces = Vec::new();
        for idx in &face_indices {
            let a = vertices[idx[0]];
            let b = vertices[idx[1]];
            let c = vertices[idx[2]];
            let normal = epa_face_normal(a, b, c);
            let mut distance = epa_dot3(normal, a);
            let mut n = normal;
            if distance < 0.0 {
                n = epa_negate3(n);
                distance = -distance;
            }
            faces.push(EpaFaceRaw {
                vertices: *idx,
                normal: n,
                distance,
            });
        }
        Self { vertices, faces }
    }
    /// Find the index of the closest face to the origin.
    pub fn find_closest_face(&self) -> usize {
        self.faces
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.distance
                    .partial_cmp(&b.distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Return `(face_index, distance, normal)` for the face closest to the origin.
    pub fn closest_face(&self) -> (usize, f64, [f64; 3]) {
        let idx = self.find_closest_face();
        if self.faces.is_empty() {
            return (0, 0.0, [0.0, 1.0, 0.0]);
        }
        let f = &self.faces[idx];
        (idx, f.distance, f.normal)
    }
    /// Expand the polytope toward a new support point.
    ///
    /// Returns `true` if the polytope was successfully expanded, `false` if the
    /// support point did not make progress.
    pub fn expand(&mut self, support_point: [f64; 3]) -> bool {
        let new_idx = self.vertices.len();
        let mut edges: Vec<(usize, usize)> = Vec::new();
        let mut i = 0;
        let mut any_removed = false;
        while i < self.faces.len() {
            let face = &self.faces[i];
            let v = self.vertices[face.vertices[0]];
            let to_point = epa_sub3(support_point, v);
            if epa_dot3(face.normal, to_point) > 1e-10 {
                let [a, b, c] = face.vertices;
                epa_add_edge(&mut edges, a, b);
                epa_add_edge(&mut edges, b, c);
                epa_add_edge(&mut edges, c, a);
                self.faces.swap_remove(i);
                any_removed = true;
            } else {
                i += 1;
            }
        }
        if !any_removed || edges.is_empty() {
            return false;
        }
        self.vertices.push(support_point);
        for &(ea, eb) in &edges {
            let a = self.vertices[ea];
            let b = self.vertices[eb];
            let c = support_point;
            let normal = epa_face_normal(a, b, c);
            let mut distance = epa_dot3(normal, a);
            let mut n = normal;
            if distance < 0.0 {
                n = epa_negate3(n);
                distance = -distance;
            }
            self.faces.push(EpaFaceRaw {
                vertices: [ea, eb, new_idx],
                normal: n,
                distance,
            });
        }
        true
    }
}
/// Configuration for EPA termination.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EpaConfig {
    /// Maximum number of expansion iterations.
    pub max_iter: usize,
    /// Convergence tolerance (distance improvement threshold).
    pub tolerance: f64,
    /// Maximum number of faces allowed in the polytope.
    pub max_faces: usize,
}
#[allow(dead_code)]
impl EpaConfig {
    /// Default EPA configuration.
    pub fn default() -> Self {
        Self {
            max_iter: 64,
            tolerance: 1e-6,
            max_faces: 256,
        }
    }
    /// Create a high-precision configuration.
    pub fn high_precision() -> Self {
        Self {
            max_iter: 128,
            tolerance: 1e-9,
            max_faces: 512,
        }
    }
}
/// Summary statistics about an EPA polytope's geometry.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct EpaPolytopeStats {
    /// Number of faces.
    pub face_count: usize,
    /// Number of vertices.
    pub vertex_count: usize,
    /// Minimum face distance from origin.
    pub min_distance: f64,
    /// Maximum face distance from origin.
    pub max_distance: f64,
    /// Average face distance from origin.
    pub avg_distance: f64,
    /// Minimum face area.
    pub min_area: f64,
    /// Maximum face area.
    pub max_area: f64,
}
/// Witness point data: contact positions on each shape.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EpaWitness {
    /// Contact point on shape A.
    pub point_a: [f64; 3],
    /// Contact point on shape B.
    pub point_b: [f64; 3],
    /// Contact normal (unit vector, from B to A).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
}
/// EPA algorithm for computing penetration depth and contact normal.
pub struct Epa;
impl Epa {
    /// Compute penetration depth and contact information from a GJK simplex
    /// that has been determined to contain the origin.
    pub fn penetration_depth(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
        simplex: &Simplex,
    ) -> Option<Contact> {
        if simplex.len() < 4 {
            return Self::fallback_contact(simplex);
        }
        let mut vertices: Vec<SupportPoint> = simplex.points.clone();
        let mut faces = build_initial_faces(&vertices);
        for _ in 0..MAX_ITERATIONS {
            let closest_idx = match faces.iter().enumerate().min_by(|(_, a), (_, b)| {
                a.distance
                    .partial_cmp(&b.distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }) {
                Some((idx, _)) => idx,
                None => return None,
            };
            let closest = faces[closest_idx].clone();
            let search_dir = closest.normal;
            let new_point =
                support_minkowski(shape_a, transform_a, shape_b, transform_b, &search_dir);
            let new_dist = new_point.point.dot(&search_dir);
            if (new_dist - closest.distance).abs() < TOLERANCE {
                return Some(build_contact(&closest, &vertices));
            }
            let new_idx = vertices.len();
            vertices.push(new_point);
            let mut edges: Vec<(usize, usize)> = Vec::new();
            let mut i = 0;
            while i < faces.len() {
                let face = &faces[i];
                let v = vertices[face.indices[0]].point;
                if face.normal.dot(&(new_point.point - v)) > 0.0 {
                    let [a, b, c] = face.indices;
                    add_edge(&mut edges, a, b);
                    add_edge(&mut edges, b, c);
                    add_edge(&mut edges, c, a);
                    faces.swap_remove(i);
                } else {
                    i += 1;
                }
            }
            for &(ea, eb) in &edges {
                let normal = compute_face_normal(
                    &vertices[ea].point,
                    &vertices[eb].point,
                    &vertices[new_idx].point,
                );
                let distance = normal.dot(&vertices[ea].point);
                let (normal, distance) = if distance < 0.0 {
                    (-normal, -distance)
                } else {
                    (normal, distance)
                };
                faces.push(Face {
                    indices: [ea, eb, new_idx],
                    normal,
                    distance,
                });
            }
        }
        faces
            .iter()
            .min_by(|a, b| {
                a.distance
                    .partial_cmp(&b.distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|face| build_contact(face, &vertices))
    }
    fn fallback_contact(simplex: &Simplex) -> Option<Contact> {
        if simplex.is_empty() {
            return None;
        }
        let (closest, _) = simplex.closest_point_to_origin();
        let depth = closest.norm();
        let normal = if depth > 1e-10 {
            -closest.normalize()
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let point_a: Vec3 =
            simplex.points.iter().map(|p| p.support_a).sum::<Vec3>() / simplex.len() as Real;
        let point_b: Vec3 =
            simplex.points.iter().map(|p| p.support_b).sum::<Vec3>() / simplex.len() as Real;
        Some(Contact::new(point_a, point_b, normal, depth))
    }
}
/// Statistics collected during an EPA run.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct EpaStats {
    /// Number of expansion iterations performed.
    pub iterations: usize,
    /// Final number of faces in the polytope.
    pub face_count: usize,
    /// Final penetration depth.
    pub final_depth: f64,
    /// Whether the algorithm converged within the tolerance.
    pub converged: bool,
}
/// A triangle face of the EPA polytope.
#[derive(Debug, Clone)]
pub(super) struct Face {
    pub(super) indices: [usize; 3],
    pub(super) normal: Vec3,
    pub(super) distance: Real,
}
/// A priority-queue-like structure that tracks the EPA face with minimum distance.
///
/// In practice, EPA typically scans all faces each iteration (O(n)), but this
/// wrapper exposes a clean interface to find the best face.
#[allow(dead_code)]
pub struct EpaFaceQueue<'a> {
    pub(super) polytope: &'a EpaPolytope,
}
#[allow(dead_code)]
impl<'a> EpaFaceQueue<'a> {
    /// Create a new face queue backed by the given polytope.
    pub fn new(polytope: &'a EpaPolytope) -> Self {
        Self { polytope }
    }
    /// Return the index of the face with minimum distance from origin.
    pub fn pop_min(&self) -> Option<usize> {
        if self.polytope.faces.is_empty() {
            return None;
        }
        Some(self.polytope.find_closest_face())
    }
    /// Return `(index, distance, normal)` for the closest face.
    pub fn peek_min(&self) -> Option<(usize, f64, [f64; 3])> {
        let idx = self.pop_min()?;
        let f = &self.polytope.faces[idx];
        Some((idx, f.distance, f.normal))
    }
    /// Number of faces in the backing polytope.
    pub fn len(&self) -> usize {
        self.polytope.faces.len()
    }
    /// Whether the backing polytope has no faces.
    pub fn is_empty(&self) -> bool {
        self.polytope.faces.is_empty()
    }
}
/// Penetration result from the EPA algorithm.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EpaPenetration {
    /// Penetration normal (unit vector, points from B to A).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Witness point on shape A.
    pub witness_a: [f64; 3],
    /// Witness point on shape B.
    pub witness_b: [f64; 3],
}
