// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Discrete Helmholtz-Hodge decomposition of tangent vector fields on triangle meshes;
//! see [`HodgeDecomposition`] for the method.

use super::support::{Components, EdgeTable, ScaledSpd};
use super::{
    DiscreteMesh, FaceVectorField, VertexVectorField, add3, cross3, dot3, len3, scale3, sub3,
};
use crate::error::{Error, Result};

/// Result of the discrete Helmholtz-Hodge decomposition of a tangent vector field.
///
/// A tangent vector field `X`, constant on every triangle, is split into three
/// mutually `L2`-orthogonal parts (Polthier & Preuss, *Identifying vector field
/// singularities using a discrete Hodge decomposition*, 2003):
///
/// ```text
/// X = grad f  +  N x grad psi  +  h
/// ```
///
/// * the **exact** part `grad f` is the gradient of a continuous piecewise linear
///   function `f` (one value per vertex),
/// * the **co-exact** part `N x grad psi` is the gradient of a *non-conforming*
///   piecewise linear function `psi` rotated by 90 degrees about the face normal `N`;
///   `psi` is a Crouzeix-Raviart function, continuous at the midpoints of the edges
///   and described by one value per edge,
/// * the **harmonic** part `h` is what remains.  It is orthogonal to both other
///   parts, hence discretely divergence free and curl free.
///
/// `f` and `psi` minimise `sum_t A_t |grad f - X_t|^2` and `sum_t A_t |N x grad psi -
/// X_t|^2`, i.e. they solve the symmetric positive definite stiffness systems
/// `S f = b` and `C psi = c` whose right-hand sides are the integrated divergence of
/// `X` and of the rotated field `-N x X`.
///
/// The pairing "conforming `f`, Crouzeix-Raviart `psi`" is what makes the
/// decomposition **complete**: the face-constant tangent fields have dimension `2F`
/// and the dimensions of the three parts add up to exactly that, because by
/// Euler's formula `(V - 1) + (E - 1) + 2 genus = 2F` on a closed surface.  The
/// harmonic part is therefore the *topological* part only: it vanishes identically
/// on every sphere, disc and other surface of genus zero, for **every** input
/// field, and spans a space of dimension `2 genus` on a closed surface
/// (`2 genus + loops - 1` on a connected surface with that many boundary loops).
/// With two conforming spaces the dimensions would only add up to `2V - 2`, about
/// half of `2F`, and most of any field would end up in the "harmonic" remainder.
///
/// **Boundaries.**  The potential `f` takes the natural (Neumann) boundary
/// condition; the stream function takes `psi = 0` at the edges of the boundary, so
/// the co-exact part is tangent to the boundary.  This is the usual splitting of a
/// flow into a potential part and a divergence-free part that does not cross the
/// boundary.  The constant of a potential is free where no condition fixes it (`f`
/// always, `psi` on closed components); it is chosen so that the area mean is zero.
///
/// The triangles must form a manifold with consistently oriented triangles
/// (counter-clockwise seen from the side their normal points to).  Triangles of
/// zero area carry no field and are ignored.
///
/// `exact_faces[t] + co_exact_faces[t] + harmonic_faces[t]` is the tangential part
/// of the input field on triangle `t`.
///
/// Computed by [`hodge_decompose_faces`] and [`hodge_decompose`].
#[derive(Debug, Clone)]
pub struct HodgeDecomposition {
    /// Potential `f` of the exact part, one value per vertex.
    pub scalar_potential: Vec<f64>,
    /// Stream function `psi` of the co-exact part: the value of the
    /// Crouzeix-Raviart function at the midpoint of every edge of
    /// [`edges`](Self::edges).
    pub stream_function: Vec<f64>,
    /// The edges the stream function is defined on, as `(low, high)` vertex pairs in
    /// lexicographic order.  Edges of zero-area triangles are not included.
    pub edges: Vec<(usize, usize)>,
    /// Exact part `grad f` on every triangle.
    pub exact_faces: Vec<[f64; 3]>,
    /// Co-exact part `N x grad psi` on every triangle.
    pub co_exact_faces: Vec<[f64; 3]>,
    /// Harmonic part on every triangle.
    pub harmonic_faces: Vec<[f64; 3]>,
    /// Harmonic part per vertex: the area-weighted mean of
    /// [`harmonic_faces`](Self::harmonic_faces) over the triangles around the vertex.
    pub harmonic: Vec<[f64; 3]>,
}

/// Normal, area and hat-function gradients of a triangle that has an area.
struct HatBasis {
    /// Unit normal (right-hand rule on the vertex order).
    normal: [f64; 3],
    /// Area.
    area: f64,
    /// Gradients of the three linear hat functions of the corners.
    grad: [[f64; 3]; 3],
}

impl HatBasis {
    /// With `e_c` the counter-clockwise edge opposite corner `c`, the gradient of the
    /// hat function of corner `c` is `N x e_c / (2 A)`.  `None` for a (numerically)
    /// degenerate triangle.
    fn new(vertices: &[[f64; 3]], tri: &[usize; 3]) -> Option<Self> {
        let p = [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]];
        let opposite = [sub3(p[2], p[1]), sub3(p[0], p[2]), sub3(p[1], p[0])];
        let n = cross3(sub3(p[1], p[0]), sub3(p[2], p[0]));
        let twice_area = len3(n);
        let longest_sq = opposite
            .iter()
            .map(|e| dot3(*e, *e))
            .fold(0.0_f64, f64::max);
        let usable = twice_area.is_finite() && twice_area > 1e-14 * longest_sq;
        if !usable {
            return None;
        }
        let normal = scale3(n, 1.0 / twice_area);
        Some(Self {
            normal,
            area: 0.5 * twice_area,
            grad: std::array::from_fn(|c| scale3(cross3(normal, opposite[c]), 1.0 / twice_area)),
        })
    }
}

/// A triangle that has an area.
struct ActiveFace {
    /// Index of the triangle in the mesh.
    index: usize,
    /// The three vertices.
    corners: [usize; 3],
    /// The edge opposite each corner, as an index into the edge table.
    edges: [usize; 3],
    /// Geometry of the triangle.
    basis: HatBasis,
}

/// The usable triangles of a mesh with their edge table and components.
struct Surface<'a> {
    mesh: &'a DiscreteMesh,
    faces: Vec<ActiveFace>,
    table: EdgeTable,
    components: Components,
}

/// Reject non-finite vertices and triangles that reference missing vertices.
fn validate_mesh(mesh: &DiscreteMesh) -> Result<()> {
    let n = mesh.num_vertices();
    for (i, v) in mesh.vertices.iter().enumerate() {
        if !v.iter().all(|c| c.is_finite()) {
            return Err(Error::invalid_mesh(format!(
                "vertex {i} has a non-finite coordinate"
            )));
        }
    }
    for (f, tri) in mesh.triangles.iter().enumerate() {
        if let Some(&bad) = tri.iter().find(|&&v| v >= n) {
            return Err(Error::invalid_mesh(format!(
                "triangle {f} references vertex {bad} but the mesh has {n} vertices"
            )));
        }
    }
    Ok(())
}

impl<'a> Surface<'a> {
    /// # Errors
    /// [`Error::InvalidMesh`] for a malformed mesh, a non-manifold edge, or triangles
    /// that are not consistently oriented.
    fn new(mesh: &'a DiscreteMesh) -> Result<Self> {
        validate_mesh(mesh)?;
        let n = mesh.num_vertices();
        let mut indices: Vec<usize> = Vec::new();
        let mut bases: Vec<HatBasis> = Vec::new();
        for (index, tri) in mesh.triangles.iter().enumerate() {
            if let Some(basis) = HatBasis::new(&mesh.vertices, tri) {
                indices.push(index);
                bases.push(basis);
            }
        }
        let corner_lists: Vec<[usize; 3]> = indices.iter().map(|&i| mesh.triangles[i]).collect();
        let table = EdgeTable::build(n, &corner_lists)?;
        if !table.consistently_oriented {
            return Err(Error::invalid_mesh(
                "the triangles are not consistently oriented (two triangles traverse a shared \
                 edge in the same direction); the co-exact part needs a global orientation",
            ));
        }
        let components = Components::label(n, &corner_lists);
        let faces = indices
            .into_iter()
            .zip(bases)
            .zip(corner_lists.iter().zip(table.face_edges.iter()))
            .map(|((index, basis), (corners, edges))| ActiveFace {
                index,
                corners: *corners,
                edges: *edges,
                basis,
            })
            .collect();
        Ok(Self {
            mesh,
            faces,
            table,
            components,
        })
    }

    /// Vertices pinned to remove the constant of the potential: the first vertex of
    /// every connected component.
    fn pinned_vertices(&self) -> Vec<bool> {
        let mut seen = vec![false; self.components.count];
        let mut pinned = vec![false; self.mesh.num_vertices()];
        for (v, &c) in self.components.id.iter().enumerate() {
            if !seen[c] {
                seen[c] = true;
                pinned[v] = true;
            }
        }
        pinned
    }

    /// Whether a component has a boundary edge, i.e. an edge of a single triangle.
    fn component_has_boundary(&self) -> Vec<bool> {
        let mut has_boundary = vec![false; self.components.count];
        for (&(low, _), &count) in self.table.edges.iter().zip(self.table.face_count.iter()) {
            if count == 1 {
                has_boundary[self.components.id[low]] = true;
            }
        }
        has_boundary
    }

    /// Edges pinned to zero for the stream function: all boundary edges, and the
    /// first edge of every component without boundary (to remove the constant).
    fn pinned_edges(&self) -> Vec<bool> {
        let has_boundary = self.component_has_boundary();
        let mut closed_seen = vec![false; self.components.count];
        self.table
            .edges
            .iter()
            .zip(self.table.face_count.iter())
            .map(|(&(low, _), &count)| {
                let c = self.components.id[low];
                let first_of_closed = !has_boundary[c] && !closed_seen[c];
                if first_of_closed {
                    closed_seen[c] = true;
                }
                count == 1 || first_of_closed
            })
            .collect()
    }

    /// Per component: the integral over its triangles of a function whose mean over a
    /// triangle is `mean_over_face`, and the area of the component.
    fn integrals(&self, mean_over_face: impl Fn(&ActiveFace) -> f64) -> (Vec<f64>, Vec<f64>) {
        let mut integral = vec![0.0_f64; self.components.count];
        let mut area = vec![0.0_f64; self.components.count];
        for face in &self.faces {
            let c = self.components.id[face.corners[0]];
            integral[c] += face.basis.area * mean_over_face(face);
            area[c] += face.basis.area;
        }
        (integral, area)
    }

    /// Decompose the face-constant field `field` (one vector per triangle of the
    /// mesh; normal components are discarded).
    fn decompose(&self, field: &[[f64; 3]]) -> Result<HodgeDecomposition> {
        let n = self.mesh.num_vertices();
        let edge_count = self.table.edges.len();

        // Tangential parts and the two least-squares problems.
        let mut tangent = vec![[0.0_f64; 3]; self.mesh.triangles.len()];
        let mut potential_entries = Vec::with_capacity(9 * self.faces.len());
        let mut stream_entries = Vec::with_capacity(9 * self.faces.len());
        let mut potential_rhs = vec![0.0_f64; n];
        let mut stream_rhs = vec![0.0_f64; edge_count];
        for face in &self.faces {
            let basis = &face.basis;
            let x = field[face.index];
            let x = sub3(x, scale3(basis.normal, dot3(x, basis.normal)));
            tangent[face.index] = x;
            // The rotated field whose divergence drives the stream function: -N x X.
            let rotated = scale3(cross3(basis.normal, x), -1.0);
            for (a, grad_a) in basis.grad.iter().enumerate() {
                potential_rhs[face.corners[a]] += basis.area * dot3(x, *grad_a);
                // The Crouzeix-Raviart basis function of the edge opposite corner a is
                // 1 - 2 phi_a, whose gradient is -2 grad phi_a.
                stream_rhs[face.edges[a]] += -2.0 * basis.area * dot3(rotated, *grad_a);
                for (b, grad_b) in basis.grad.iter().enumerate() {
                    let stiffness = basis.area * dot3(*grad_a, *grad_b);
                    potential_entries.push((face.corners[a], face.corners[b], stiffness));
                    stream_entries.push((face.edges[a], face.edges[b], 4.0 * stiffness));
                }
            }
        }

        let mut potential = ScaledSpd::new(n, &potential_entries, &self.pinned_vertices())
            .solve(&potential_rhs, "Hodge decomposition: scalar potential")?;
        let mut stream = ScaledSpd::new(edge_count, &stream_entries, &self.pinned_edges())
            .solve(&stream_rhs, "Hodge decomposition: stream function")?;

        // Fix the free constants: zero area mean per component (for the stream function
        // only on components without boundary, where it is free).
        let (integral, area) =
            self.integrals(|face| face.corners.iter().map(|&v| potential[v]).sum::<f64>() / 3.0);
        for (v, value) in potential.iter_mut().enumerate() {
            let c = self.components.id[v];
            if area[c] > 0.0 {
                *value -= integral[c] / area[c];
            }
        }
        let has_boundary = self.component_has_boundary();
        let (integral, area) =
            self.integrals(|face| face.edges.iter().map(|&e| stream[e]).sum::<f64>() / 3.0);
        for (&(low, _), value) in self.table.edges.iter().zip(stream.iter_mut()) {
            let c = self.components.id[low];
            if !has_boundary[c] && area[c] > 0.0 {
                *value -= integral[c] / area[c];
            }
        }

        // The three parts on every triangle.
        let triangle_count = self.mesh.triangles.len();
        let mut exact_faces = vec![[0.0_f64; 3]; triangle_count];
        let mut co_exact_faces = vec![[0.0_f64; 3]; triangle_count];
        let mut harmonic_faces = vec![[0.0_f64; 3]; triangle_count];
        let mut harmonic = vec![[0.0_f64; 3]; n];
        let mut vertex_area = vec![0.0_f64; n];
        for face in &self.faces {
            let basis = &face.basis;
            let mut grad_f = [0.0_f64; 3];
            let mut grad_psi = [0.0_f64; 3];
            for (a, grad_a) in basis.grad.iter().enumerate() {
                grad_f = add3(grad_f, scale3(*grad_a, potential[face.corners[a]]));
                grad_psi = add3(grad_psi, scale3(*grad_a, -2.0 * stream[face.edges[a]]));
            }
            let co_exact = cross3(basis.normal, grad_psi);
            let rest = sub3(sub3(tangent[face.index], grad_f), co_exact);
            exact_faces[face.index] = grad_f;
            co_exact_faces[face.index] = co_exact;
            harmonic_faces[face.index] = rest;
            for &v in &face.corners {
                harmonic[v] = add3(harmonic[v], scale3(rest, basis.area));
                vertex_area[v] += basis.area;
            }
        }
        for (h, &a) in harmonic.iter_mut().zip(vertex_area.iter()) {
            if a > 0.0 {
                *h = scale3(*h, 1.0 / a);
            }
        }
        Ok(HodgeDecomposition {
            scalar_potential: potential,
            stream_function: stream,
            edges: self.table.edges.clone(),
            exact_faces,
            co_exact_faces,
            harmonic_faces,
            harmonic,
        })
    }
}

/// Hodge decomposition of a tangent vector field with one vector per triangle.
///
/// The normal component of the input vectors is discarded.  See
/// [`HodgeDecomposition`] for the discretisation, the boundary conditions and the
/// meaning of the three parts.
///
/// # Errors
/// * [`Error::LengthMismatch`] if `field` does not have one vector per triangle;
/// * [`Error::InvalidMesh`] for a non-finite vertex, a triangle that references a
///   missing vertex, an edge shared by more than two triangles, or triangles that
///   are not consistently oriented;
/// * [`Error::General`] for a non-finite field vector;
/// * [`Error::ConvergenceFailure`] if a linear solve does not converge.
pub fn hodge_decompose_faces(
    mesh: &DiscreteMesh,
    field: &FaceVectorField,
) -> Result<HodgeDecomposition> {
    if field.vectors.len() != mesh.num_faces() {
        return Err(Error::length_mismatch(
            mesh.num_faces(),
            field.vectors.len(),
        ));
    }
    if let Some(bad) = field
        .vectors
        .iter()
        .position(|v| !v.iter().all(|c| c.is_finite()))
    {
        return Err(Error::general(format!(
            "the vector of triangle {bad} is not finite"
        )));
    }
    Surface::new(mesh)?.decompose(&field.vectors)
}

/// Hodge decomposition of a tangent vector field given at the vertices.
///
/// The field of a triangle is the mean of the vectors at its three vertices (the
/// mean value of their linear interpolation); then [`hodge_decompose_faces`]
/// applies.  The vectors need not be exactly tangent: the normal component is
/// discarded.
///
/// # Errors
/// [`Error::LengthMismatch`] if `field` does not have one vector per vertex, and
/// the errors of [`hodge_decompose_faces`].
pub fn hodge_decompose(
    mesh: &DiscreteMesh,
    field: &VertexVectorField,
) -> Result<HodgeDecomposition> {
    let n = mesh.num_vertices();
    if field.vectors.len() != n {
        return Err(Error::length_mismatch(n, field.vectors.len()));
    }
    validate_mesh(mesh)?;
    let per_face: Vec<[f64; 3]> = mesh
        .triangles
        .iter()
        .map(|t| {
            let sum = add3(
                add3(field.vectors[t[0]], field.vectors[t[1]]),
                field.vectors[t[2]],
            );
            scale3(sum, 1.0 / 3.0)
        })
        .collect();
    hodge_decompose_faces(mesh, &FaceVectorField::from_vec(per_face))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discrete_geometry::test_support::{XorShift, discrete, open_cylinder, torus};
    use crate::test_meshes::{flat_grid, icosphere};

    fn random_field(mesh: &DiscreteMesh, seed: u64) -> FaceVectorField {
        let mut rng = XorShift::new(seed);
        FaceVectorField::from_vec(
            (0..mesh.num_faces())
                .map(|_| [rng.signed(), rng.signed(), rng.signed()])
                .collect(),
        )
    }

    /// The tangential part of a face field.
    fn tangent_of(mesh: &DiscreteMesh, field: &[[f64; 3]]) -> Vec<[f64; 3]> {
        field
            .iter()
            .enumerate()
            .map(|(f, x)| {
                let n = mesh.face_normal(f);
                sub3(*x, scale3(n, dot3(*x, n)))
            })
            .collect()
    }

    fn triangle_area(mesh: &DiscreteMesh, t: &[usize; 3]) -> f64 {
        0.5 * len3(cross3(
            sub3(mesh.vertices[t[1]], mesh.vertices[t[0]]),
            sub3(mesh.vertices[t[2]], mesh.vertices[t[0]]),
        ))
    }

    /// The L2 inner product of two face fields.
    fn inner(mesh: &DiscreteMesh, a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
        a.iter()
            .zip(b.iter())
            .zip(mesh.triangles.iter())
            .map(|((x, y), t)| triangle_area(mesh, t) * dot3(*x, *y))
            .sum()
    }

    fn norm(mesh: &DiscreteMesh, a: &[[f64; 3]]) -> f64 {
        inner(mesh, a, a).sqrt()
    }

    /// Gradient on every triangle of the linear interpolation of vertex values.
    fn gradient_of(mesh: &DiscreteMesh, values: &[f64]) -> Vec<[f64; 3]> {
        mesh.triangles
            .iter()
            .map(|tri| {
                let basis = HatBasis::new(&mesh.vertices, tri).expect("non-degenerate triangle");
                let mut g = [0.0; 3];
                for (a, grad) in basis.grad.iter().enumerate() {
                    g = add3(g, scale3(*grad, values[tri[a]]));
                }
                g
            })
            .collect()
    }

    fn rotated(mesh: &DiscreteMesh, field: &[[f64; 3]]) -> Vec<[f64; 3]> {
        field
            .iter()
            .enumerate()
            .map(|(f, x)| cross3(mesh.face_normal(f), *x))
            .collect()
    }

    fn random_values(n: usize, seed: u64) -> Vec<f64> {
        let mut rng = XorShift::new(seed);
        (0..n).map(|_| rng.signed()).collect()
    }

    /// Area-weighted mean of a vertex function.
    fn area_mean(mesh: &DiscreteMesh, values: &[f64]) -> f64 {
        let (mut integral, mut area) = (0.0, 0.0);
        for t in &mesh.triangles {
            let a = triangle_area(mesh, t);
            integral += a * (values[t[0]] + values[t[1]] + values[t[2]]) / 3.0;
            area += a;
        }
        integral / area
    }

    fn boundary_flags(mesh: &DiscreteMesh) -> Vec<bool> {
        EdgeTable::build(mesh.num_vertices(), &mesh.triangles)
            .expect("manifold mesh")
            .boundary_vertices(mesh.num_vertices())
    }

    fn test_surfaces() -> Vec<(&'static str, DiscreteMesh)> {
        vec![
            ("sphere", discrete(icosphere(2))),
            ("torus", discrete(torus(24, 12, 3.0, 1.0))),
            ("disc", discrete(flat_grid(8, 8, 1.0))),
            ("annulus", discrete(open_cylinder(16, 5, 1.0, 2.0))),
        ]
    }

    /// Number of linearly independent vectors among `vectors`, counting a direction
    /// only if the part of a vector that is orthogonal to the previous directions is
    /// longer than `threshold`.
    fn numerical_rank(vectors: &[Vec<f64>], threshold: f64) -> usize {
        let mut basis: Vec<Vec<f64>> = Vec::new();
        for v in vectors {
            let mut r = v.clone();
            for b in &basis {
                let c: f64 = r.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
                for (ri, bi) in r.iter_mut().zip(b.iter()) {
                    *ri -= c * bi;
                }
            }
            let len = r.iter().map(|x| x * x).sum::<f64>().sqrt();
            if len > threshold {
                basis.push(r.iter().map(|x| x / len).collect());
            }
        }
        basis.len()
    }

    fn flatten(field: &[[f64; 3]]) -> Vec<f64> {
        field.iter().flat_map(|v| v.iter().copied()).collect()
    }

    #[test]
    fn random_field_on_a_sphere_has_no_harmonic_part() {
        // On a genus-zero surface the exact and co-exact parts exhaust every face field.
        let m = discrete(icosphere(2));
        let field = random_field(&m, 1);
        let d = hodge_decompose_faces(&m, &field).expect("valid mesh");
        let x = tangent_of(&m, &field.vectors);
        let scale = norm(&m, &x);
        assert!(
            norm(&m, &d.harmonic_faces) < 1e-8 * scale,
            "harmonic part {} of {scale}",
            norm(&m, &d.harmonic_faces)
        );
        assert!(norm(&m, &d.exact_faces) > 0.2 * scale);
        assert!(norm(&m, &d.co_exact_faces) > 0.2 * scale);
        let worst = d.harmonic.iter().map(|h| len3(*h)).fold(0.0_f64, f64::max);
        assert!(worst < 1e-8, "vertex harmonic field {worst}");
    }

    #[test]
    fn parts_add_up_and_are_orthogonal_on_every_kind_of_surface() {
        for (name, m) in test_surfaces() {
            let field = random_field(&m, 5);
            let d = hodge_decompose_faces(&m, &field).expect("valid mesh");
            let x = tangent_of(&m, &field.vectors);
            for (f, x_f) in x.iter().enumerate() {
                let sum = add3(
                    add3(d.exact_faces[f], d.co_exact_faces[f]),
                    d.harmonic_faces[f],
                );
                assert!(len3(sub3(sum, *x_f)) < 1e-12, "{name}: face {f}");
            }
            let scale = inner(&m, &x, &x);
            let (e, c, h) = (&d.exact_faces, &d.co_exact_faces, &d.harmonic_faces);
            for (label, value) in [
                ("exact.co-exact", inner(&m, e, c)),
                ("exact.harmonic", inner(&m, e, h)),
                ("co-exact.harmonic", inner(&m, c, h)),
            ] {
                assert!(
                    value.abs() < 1e-8 * scale,
                    "{name}: {label} = {value} (scale {scale})"
                );
            }
            let parts = inner(&m, e, e) + inner(&m, c, c) + inner(&m, h, h);
            assert!((parts - scale).abs() < 1e-8 * scale, "{name}: Pythagoras");
            assert_eq!(d.scalar_potential.len(), m.num_vertices());
            assert_eq!(d.stream_function.len(), d.edges.len());
            assert_eq!(d.exact_faces.len(), m.num_faces());
            assert_eq!(d.harmonic.len(), m.num_vertices());
        }
    }

    #[test]
    fn gradient_of_a_vertex_function_is_recovered_with_its_potential() {
        for (name, m) in test_surfaces() {
            let values = random_values(m.num_vertices(), 9);
            let x = gradient_of(&m, &values);
            let d = hodge_decompose_faces(&m, &FaceVectorField::from_vec(x.clone()))
                .expect("valid mesh");
            let scale = norm(&m, &x);
            assert!(
                norm(&m, &d.co_exact_faces) < 1e-8 * scale,
                "{name}: co-exact"
            );
            assert!(
                norm(&m, &d.harmonic_faces) < 1e-8 * scale,
                "{name}: harmonic"
            );
            // f0 up to a constant: the potential has zero area mean.
            let mean = area_mean(&m, &values);
            for (v, (&f, &f0)) in d.scalar_potential.iter().zip(values.iter()).enumerate() {
                assert!(
                    (f - (f0 - mean)).abs() < 1e-8,
                    "{name}: vertex {v}: {f} vs {}",
                    f0 - mean
                );
            }
            assert!(area_mean(&m, &d.scalar_potential).abs() < 1e-12);
        }
    }

    #[test]
    fn rotated_gradient_of_a_vertex_function_is_purely_co_exact() {
        for (name, m) in test_surfaces() {
            // On surfaces with boundary the stream function vanishes on the boundary.
            let boundary = boundary_flags(&m);
            let values: Vec<f64> = random_values(m.num_vertices(), 21)
                .iter()
                .zip(boundary.iter())
                .map(|(&g, &b)| if b { 0.0 } else { g })
                .collect();
            let x = rotated(&m, &gradient_of(&m, &values));
            let d = hodge_decompose_faces(&m, &FaceVectorField::from_vec(x.clone()))
                .expect("valid mesh");
            let scale = norm(&m, &x);
            assert!(norm(&m, &d.exact_faces) < 1e-8 * scale, "{name}: exact");
            assert!(
                norm(&m, &d.harmonic_faces) < 1e-8 * scale,
                "{name}: harmonic"
            );
            assert!(
                d.scalar_potential.iter().all(|f| f.abs() < 1e-8),
                "{name}: the potential of a divergence-free field is constant"
            );
            let mut worst = 0.0_f64;
            for (f, x_f) in x.iter().enumerate() {
                worst = worst.max(len3(sub3(d.co_exact_faces[f], *x_f)));
                // The old code subtracted a vector along the vertex normal; the real
                // co-exact part is tangent to its face.
                assert!(dot3(d.co_exact_faces[f], m.face_normal(f)).abs() < 1e-10);
            }
            assert!(
                worst < 1e-8,
                "{name}: co-exact part differs from the input by {worst}"
            );
            // psi at an edge midpoint is the mean of the end values of a linear function.
            let offsets: Vec<f64> = d
                .edges
                .iter()
                .zip(d.stream_function.iter())
                .map(|(&(a, b), &psi)| psi - 0.5 * (values[a] + values[b]))
                .collect();
            let low = offsets.iter().copied().fold(f64::INFINITY, f64::min);
            let high = offsets.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(
                high - low < 1e-8,
                "{name}: stream function offset spread {}",
                high - low
            );
        }
    }

    #[test]
    fn harmonic_part_has_the_dimension_of_the_first_betti_number() {
        // Sphere: 0; torus: 2; disc: 0; annulus (two boundary loops): 1.
        let expected = [
            ("sphere", 0usize),
            ("torus", 2),
            ("disc", 0),
            ("annulus", 1),
        ];
        for ((name, m), (expected_name, rank)) in test_surfaces().into_iter().zip(expected) {
            assert_eq!(name, expected_name);
            let mut parts: Vec<Vec<f64>> = Vec::new();
            let mut field_size = 0.0_f64;
            for seed in 0..5 {
                let field = random_field(&m, 100 + seed);
                let d = hodge_decompose_faces(&m, &field).expect("valid mesh");
                parts.push(flatten(&d.harmonic_faces));
                let x = flatten(&tangent_of(&m, &field.vectors));
                field_size = field_size.max(x.iter().map(|c| c * c).sum::<f64>().sqrt());
            }
            // A direction counts if it is longer than 1e-7 of the size of the input.
            assert_eq!(numerical_rank(&parts, 1e-7 * field_size), rank, "{name}");
        }
    }

    #[test]
    fn harmonic_part_is_its_own_decomposition() {
        let m = discrete(torus(24, 12, 3.0, 1.0));
        let d = hodge_decompose_faces(&m, &random_field(&m, 3)).expect("valid mesh");
        let again = hodge_decompose_faces(&m, &FaceVectorField::from_vec(d.harmonic_faces.clone()))
            .expect("valid mesh");
        let scale = norm(&m, &d.harmonic_faces);
        assert!(scale > 1e-3);
        assert!(norm(&m, &again.exact_faces) < 1e-8 * scale);
        assert!(norm(&m, &again.co_exact_faces) < 1e-8 * scale);
        for f in 0..m.num_faces() {
            assert!(len3(sub3(again.harmonic_faces[f], d.harmonic_faces[f])) < 1e-8);
        }
    }

    /// Relative deviation of the harmonic part from the azimuthal field `d(phi)` of a
    /// torus of revolution with `nu x nv` vertices.
    fn azimuthal_deviation(nu: usize, nv: usize) -> f64 {
        let m = discrete(torus(nu, nv, 3.0, 1.0));
        let mut field = VertexVectorField::zeros(m.num_vertices());
        for (v, p) in m.vertices.iter().enumerate() {
            let rho_sq = p[0] * p[0] + p[1] * p[1];
            field.vectors[v] = [-p[1] / rho_sq, p[0] / rho_sq, 0.0];
        }
        let d = hodge_decompose(&m, &field).expect("valid mesh");
        let per_face: Vec<[f64; 3]> = m
            .triangles
            .iter()
            .map(|t| {
                scale3(
                    add3(
                        add3(field.vectors[t[0]], field.vectors[t[1]]),
                        field.vectors[t[2]],
                    ),
                    1.0 / 3.0,
                )
            })
            .collect();
        let x = tangent_of(&m, &per_face);
        let missing: Vec<[f64; 3]> = x
            .iter()
            .zip(d.harmonic_faces.iter())
            .map(|(a, b)| sub3(*a, *b))
            .collect();
        norm(&m, &missing) / norm(&m, &x)
    }

    #[test]
    fn azimuthal_field_of_a_torus_is_harmonic() {
        // The differential of the azimuth around the axis is closed and co-closed on
        // a torus of revolution: all of it is harmonic, up to the discretisation,
        // which halves with every refinement of the mesh.
        let coarse = azimuthal_deviation(24, 12);
        let fine = azimuthal_deviation(48, 24);
        assert!(fine < 0.05, "relative deviation {fine}");
        assert!(fine < 0.6 * coarse, "no convergence: {coarse} -> {fine}");
    }

    #[test]
    fn smooth_field_on_a_sphere_recovers_its_potentials() {
        let m = discrete(icosphere(4));
        // X = grad_S f0 + p x e_x with f0 = z + x y: the second term is the rotated
        // surface gradient of g0 = x on the unit sphere.
        let mut field = VertexVectorField::zeros(m.num_vertices());
        for (v, p) in m.vertices.iter().enumerate() {
            let g = [p[1], p[0], 1.0];
            let grad_f0 = sub3(g, scale3(*p, dot3(g, *p)));
            let rot_g0 = cross3(*p, [1.0, 0.0, 0.0]);
            field.vectors[v] = add3(grad_f0, rot_g0);
        }
        let d = hodge_decompose(&m, &field).expect("valid mesh");
        assert!(
            d.harmonic.iter().map(|h| len3(*h)).fold(0.0, f64::max) < 1e-8,
            "a sphere has no harmonic fields"
        );
        let f0: Vec<f64> = m.vertices.iter().map(|p| p[2] + p[0] * p[1]).collect();
        let mean = area_mean(&m, &f0);
        let (mut diff, mut reference) = (0.0, 0.0);
        for (&f, &exact) in d.scalar_potential.iter().zip(f0.iter()) {
            diff += (f - (exact - mean)).powi(2);
            reference += (exact - mean).powi(2);
        }
        let potential_error = (diff / reference).sqrt();
        assert!(potential_error < 0.005, "potential error {potential_error}");
        // psi at an edge midpoint approximates g0 = x there.
        let (mut diff, mut reference) = (0.0, 0.0);
        for (&(a, b), &psi) in d.edges.iter().zip(d.stream_function.iter()) {
            let g0 = 0.5 * (m.vertices[a][0] + m.vertices[b][0]);
            diff += (psi - g0).powi(2);
            reference += g0 * g0;
        }
        let stream_error = (diff / reference).sqrt();
        assert!(stream_error < 0.005, "stream function error {stream_error}");
    }

    #[test]
    fn constant_field_on_a_flat_patch_is_a_pure_gradient() {
        let m = discrete(flat_grid(8, 6, 1.0));
        let x = vec![[1.0, 2.0, 0.0]; m.num_faces()];
        let d = hodge_decompose_faces(&m, &FaceVectorField::from_vec(x)).expect("valid mesh");
        let linear: Vec<f64> = m.vertices.iter().map(|p| p[0] + 2.0 * p[1]).collect();
        let mean = area_mean(&m, &linear);
        for (v, &f) in d.scalar_potential.iter().enumerate() {
            assert!((f - (linear[v] - mean)).abs() < 1e-9, "vertex {v}");
        }
        assert!(norm(&m, &d.co_exact_faces) < 1e-9);
        assert!(norm(&m, &d.harmonic_faces) < 1e-9);
    }

    #[test]
    fn normal_component_of_the_input_is_discarded() {
        let m = discrete(icosphere(2));
        let field = random_field(&m, 8);
        let tilted = FaceVectorField::from_vec(
            field
                .vectors
                .iter()
                .enumerate()
                .map(|(f, x)| add3(*x, scale3(m.face_normal(f), 3.0 + f as f64)))
                .collect(),
        );
        let a = hodge_decompose_faces(&m, &field).expect("valid mesh");
        let b = hodge_decompose_faces(&m, &tilted).expect("valid mesh");
        for f in 0..m.num_faces() {
            assert!(len3(sub3(a.exact_faces[f], b.exact_faces[f])) < 1e-9);
            assert!(len3(sub3(a.co_exact_faces[f], b.co_exact_faces[f])) < 1e-9);
        }
    }

    #[test]
    fn connected_components_are_decomposed_independently() {
        let (mut vertices, mut faces) = icosphere(1);
        let offset = vertices.len();
        let first_faces = faces.len();
        let (more_vertices, more_faces) = icosphere(1);
        vertices.extend(more_vertices.iter().map(|p| [p[0] + 10.0, p[1], p[2]]));
        faces.extend(
            more_faces
                .iter()
                .map(|t| [t[0] + offset, t[1] + offset, t[2] + offset]),
        );
        let m = DiscreteMesh::new(vertices, faces);
        // A gradient field on the first sphere only.
        let mut values = random_values(m.num_vertices(), 4);
        for value in &mut values[offset..] {
            *value = 0.0;
        }
        let x = gradient_of(&m, &values);
        let d = hodge_decompose_faces(&m, &FaceVectorField::from_vec(x)).expect("valid mesh");
        assert!(d.scalar_potential[offset..].iter().all(|f| f.abs() < 1e-9));
        let first = DiscreteMesh::new(
            m.vertices[..offset].to_vec(),
            m.triangles[..first_faces].to_vec(),
        );
        let second = DiscreteMesh::new(
            m.vertices[offset..].to_vec(),
            m.triangles[first_faces..]
                .iter()
                .map(|t| [t[0] - offset, t[1] - offset, t[2] - offset])
                .collect(),
        );
        let mean = area_mean(&first, &values[..offset]);
        for (v, (&f, &value)) in d.scalar_potential[..offset]
            .iter()
            .zip(values.iter())
            .enumerate()
        {
            assert!((f - (value - mean)).abs() < 1e-8, "vertex {v}");
        }
        // Each component has its own constant: zero mean on both.
        assert!(area_mean(&first, &d.scalar_potential[..offset]).abs() < 1e-12);
        assert!(area_mean(&second, &d.scalar_potential[offset..]).abs() < 1e-12);
    }

    #[test]
    fn zero_field_gives_zero_parts() {
        let m = discrete(icosphere(1));
        let d =
            hodge_decompose(&m, &VertexVectorField::zeros(m.num_vertices())).expect("valid mesh");
        assert_eq!(d.scalar_potential.len(), m.num_vertices());
        assert_eq!(d.harmonic.len(), m.num_vertices());
        assert!(d.scalar_potential.iter().all(|f| *f == 0.0));
        assert!(d.stream_function.iter().all(|f| *f == 0.0));
        assert!(d.harmonic.iter().all(|h| *h == [0.0; 3]));
    }

    #[test]
    fn meshes_without_triangles_decompose_to_nothing() {
        let none = DiscreteMesh::new(Vec::new(), Vec::new());
        let d = hodge_decompose_faces(&none, &FaceVectorField::from_vec(Vec::new()))
            .expect("an empty mesh is valid");
        assert!(d.scalar_potential.is_empty() && d.edges.is_empty() && d.harmonic.is_empty());
        let points = DiscreteMesh::new(vec![[0.0; 3], [1.0, 0.0, 0.0]], Vec::new());
        let d = hodge_decompose(&points, &VertexVectorField::zeros(2)).expect("valid mesh");
        assert_eq!(d.scalar_potential, vec![0.0, 0.0]);
        assert!(d.exact_faces.is_empty() && d.stream_function.is_empty());
    }

    #[test]
    fn invalid_input_is_reported() {
        let m = discrete(icosphere(1));
        let short = FaceVectorField::from_vec(vec![[0.0; 3]; m.num_faces() - 1]);
        assert!(matches!(
            hodge_decompose_faces(&m, &short),
            Err(Error::LengthMismatch { .. })
        ));
        assert!(matches!(
            hodge_decompose(&m, &VertexVectorField::zeros(m.num_vertices() + 1)),
            Err(Error::LengthMismatch { .. })
        ));
        let mut nan_field = random_field(&m, 1);
        nan_field.vectors[4][1] = f64::NAN;
        assert!(matches!(
            hodge_decompose_faces(&m, &nan_field),
            Err(Error::General(_))
        ));
        let good = random_field(&m, 1);
        let mut bad = m.clone();
        bad.triangles[2][1] = m.num_vertices() + 3;
        assert!(matches!(
            hodge_decompose_faces(&bad, &good),
            Err(Error::InvalidMesh { .. })
        ));
        let mut bad = m.clone();
        bad.vertices[5][2] = f64::INFINITY;
        assert!(matches!(
            hodge_decompose_faces(&bad, &good),
            Err(Error::InvalidMesh { .. })
        ));
        // One flipped triangle: the orientation is no longer consistent.
        let mut bad = m.clone();
        bad.triangles[7].swap(0, 1);
        assert!(matches!(
            hodge_decompose_faces(&bad, &good),
            Err(Error::InvalidMesh { .. })
        ));
        // A third triangle on an existing edge.
        let mut bad = m.clone();
        let edge = (bad.triangles[0][0], bad.triangles[0][1]);
        bad.vertices.push([9.0, 9.0, 9.0]);
        let apex = bad.vertices.len() - 1;
        bad.triangles.push([edge.0, edge.1, apex]);
        let more = FaceVectorField::from_vec(vec![[0.0; 3]; bad.num_faces()]);
        assert!(matches!(
            hodge_decompose_faces(&bad, &more),
            Err(Error::InvalidMesh { .. })
        ));
    }

    #[test]
    fn degenerate_triangles_and_isolated_vertices_are_ignored() {
        let base = discrete(icosphere(1));
        let n = base.num_vertices();
        let mut vertices = base.vertices.clone();
        vertices.push([40.0, 0.0, 0.0]); // isolated
        vertices.extend([[20.0, 0.0, 0.0], [21.0, 0.0, 0.0], [22.0, 0.0, 0.0]]);
        let mut triangles = base.triangles.clone();
        triangles.push([n + 1, n + 2, n + 3]); // zero area
        let m = DiscreteMesh::new(vertices, triangles);
        let field = random_field(&base, 2);
        let mut padded = field.vectors.clone();
        padded.push([1.0, 1.0, 1.0]);
        let with_junk = hodge_decompose_faces(&m, &FaceVectorField::from_vec(padded))
            .expect("degenerate triangles are skipped");
        let plain = hodge_decompose_faces(&base, &field).expect("valid mesh");
        assert_eq!(with_junk.edges, plain.edges);
        for f in 0..base.num_faces() {
            assert!(len3(sub3(with_junk.exact_faces[f], plain.exact_faces[f])) < 1e-10);
        }
        assert_eq!(with_junk.exact_faces[base.num_faces()], [0.0; 3]);
        assert_eq!(with_junk.harmonic_faces[base.num_faces()], [0.0; 3]);
        assert_eq!(with_junk.scalar_potential[n], 0.0);
        assert_eq!(with_junk.harmonic[n], [0.0; 3]);
    }

    #[test]
    fn decomposition_is_deterministic() {
        let m = discrete(torus(16, 8, 3.0, 1.0));
        let field = random_field(&m, 6);
        let a = hodge_decompose_faces(&m, &field).expect("valid mesh");
        let b = hodge_decompose_faces(&m, &field).expect("valid mesh");
        let bits = |v: &[[f64; 3]]| v.iter().flatten().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.exact_faces), bits(&b.exact_faces));
        assert_eq!(bits(&a.co_exact_faces), bits(&b.co_exact_faces));
    }

    #[test]
    fn ten_thousand_vertices_decompose_quickly() {
        let m = discrete(icosphere(5));
        assert_eq!(m.num_vertices(), 10242);
        let field = random_field(&m, 12);
        let start = std::time::Instant::now();
        let d = hodge_decompose_faces(&m, &field).expect("valid mesh");
        let elapsed = start.elapsed().as_secs_f64();
        assert!(
            elapsed < 60.0,
            "10k-vertex decomposition took {elapsed:.2} s"
        );
        let x = tangent_of(&m, &field.vectors);
        assert!(norm(&m, &d.harmonic_faces) < 1e-7 * norm(&m, &x));
    }
}
