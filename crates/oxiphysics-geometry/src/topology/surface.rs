//! Validated index-form triangle surface used by the surgery operators.
//!
//! The half-edge structure stores every adjacency twice (`next`, `prev`,
//! `twin`, `face`, `vertex_he`), so editing it in place makes it easy to leave
//! a stale pointer behind. The surgery operators therefore edit the plain
//! `vertices` / `triangles` form of the surface, check that the edited surface
//! is still a valid oriented 2-manifold, and only then rebuild the half-edge
//! mesh from it, which leaves the mesh untouched when an operation is
//! rejected.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use super::functions::{len3, sub3};
use super::types::HalfEdgeMesh;
use crate::error::{Error, Result};

/// Marker for the "vertex at infinity" that closes every boundary loop in the
/// link-condition test (a cone over the boundary turns a surface with boundary
/// into a closed one).
pub(super) const INFINITE_VERTEX: usize = usize::MAX;

/// Map from a directed edge `(u, v)` to the `(triangle, corner)` whose corner
/// `corner` starts the edge `u -> v` (`triangles[triangle][corner] == u`).
pub(super) type DirectedEdges = HashMap<(usize, usize), (usize, usize)>;

/// Location of a half-edge: the face it belongs to and its corner in that face.
pub(super) type HalfEdgeLocator = Vec<Option<(usize, usize)>>;

/// Disjoint-set forest over vertex indices (iterative path compression).
pub(super) struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    pub(super) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    pub(super) fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    pub(super) fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
            self.parent[hi] = lo;
        }
    }
}

/// An oriented triangle surface in index form.
#[derive(Debug, Clone)]
pub(super) struct Surface {
    /// Vertex positions.
    pub(super) vertices: Vec<[f64; 3]>,
    /// Counter-clockwise vertex triples.
    pub(super) triangles: Vec<[usize; 3]>,
}

impl Surface {
    /// Extract and validate the surface stored in `mesh`.
    ///
    /// Besides the surface this returns, for every half-edge index, the face and
    /// corner it corresponds to (`None` for half-edges that no face loop visits).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidMesh`] when a face loop is not a triangle, a pointer is out
    /// of range, or the surface is not a valid oriented 2-manifold (see
    /// [`Surface::validate`]).
    pub(super) fn from_mesh(mesh: &HalfEdgeMesh) -> Result<(Self, HalfEdgeLocator)> {
        let vertex_count = mesh.vertices.len();
        let mut triangles = Vec::with_capacity(mesh.faces.len());
        let mut locator: HalfEdgeLocator = vec![None; mesh.half_edges.len()];
        for (fi, face) in mesh.faces.iter().enumerate() {
            let mut loop_he = [0usize; 3];
            let mut cur = face.start_he;
            for slot in &mut loop_he {
                let he = mesh.half_edges.get(cur).ok_or_else(|| {
                    Error::invalid_mesh(format!("face {fi}: half-edge index {cur} is out of range"))
                })?;
                *slot = cur;
                cur = he.next;
            }
            if cur != face.start_he {
                return Err(Error::invalid_mesh(format!(
                    "face {fi} is not a triangle (its half-edge loop does not close after 3 steps)"
                )));
            }
            let mut triangle = [0usize; 3];
            for (corner, &he_idx) in loop_he.iter().enumerate() {
                let origin = mesh.half_edges[he_idx].origin;
                if origin >= vertex_count {
                    return Err(Error::invalid_mesh(format!(
                        "half-edge {he_idx} starts at vertex {origin} but the mesh has {vertex_count} vertices"
                    )));
                }
                if locator[he_idx].is_some() {
                    return Err(Error::invalid_mesh(format!(
                        "half-edge {he_idx} belongs to more than one face"
                    )));
                }
                locator[he_idx] = Some((fi, corner));
                triangle[corner] = origin;
            }
            triangles.push(triangle);
        }
        let surface = Self {
            vertices: mesh.vertices.clone(),
            triangles,
        };
        surface.validate()?;
        Ok((surface, locator))
    }

    /// Rebuild a half-edge mesh from this surface.
    pub(super) fn into_mesh(self) -> HalfEdgeMesh {
        HalfEdgeMesh::from_triangles(self.vertices, &self.triangles)
    }

    /// Check that the surface is a valid oriented 2-manifold triangulation.
    ///
    /// The surface must have at least one face; every face must have three
    /// distinct, in-range vertices; no two faces may share the same vertex set;
    /// every directed edge may be used by at most one face (which also bounds
    /// every undirected edge to two faces with opposite directions, i.e. a
    /// consistent orientation); and the faces around each vertex must form a
    /// single fan (an open fan at a boundary vertex, a closed one otherwise).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidMesh`] describing the first violation found.
    pub(super) fn validate(&self) -> Result<()> {
        if self.triangles.is_empty() {
            return Err(Error::invalid_mesh("the mesh has no faces"));
        }
        let vertex_count = self.vertices.len();
        let mut seen_edges: HashSet<(usize, usize)> =
            HashSet::with_capacity(self.triangles.len() * 3);
        let mut seen_faces: HashSet<[usize; 3]> = HashSet::with_capacity(self.triangles.len());
        for (fi, tri) in self.triangles.iter().enumerate() {
            if let Some(&bad) = tri.iter().find(|&&v| v >= vertex_count) {
                return Err(Error::invalid_mesh(format!(
                    "face {fi} references vertex {bad} but the mesh has {vertex_count} vertices"
                )));
            }
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                return Err(Error::invalid_mesh(format!(
                    "face {fi} is degenerate: ({}, {}, {})",
                    tri[0], tri[1], tri[2]
                )));
            }
            let mut key = *tri;
            key.sort_unstable();
            if !seen_faces.insert(key) {
                return Err(Error::invalid_mesh(format!(
                    "face {fi} repeats the vertex set {key:?} of another face"
                )));
            }
            for k in 0..3 {
                let (u, v) = (tri[k], tri[(k + 1) % 3]);
                if !seen_edges.insert((u, v)) {
                    return Err(Error::invalid_mesh(format!(
                        "directed edge {u}->{v} is used by two faces (non-manifold edge or inconsistent orientation)"
                    )));
                }
            }
        }
        self.check_vertex_fans()
    }

    /// The faces around every vertex must form one fan.
    ///
    /// The link of a vertex `v` consists of the edge `a -> b` of every face
    /// `(v, a, b)`. Because directed edges are unique, every link vertex has at
    /// most one outgoing and one incoming link edge, so the link is a disjoint
    /// union of directed paths and cycles; `v` is a manifold vertex exactly when
    /// there is a single component.
    fn check_vertex_fans(&self) -> Result<()> {
        let mut links: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for tri in &self.triangles {
            for k in 0..3 {
                links
                    .entry(tri[k])
                    .or_default()
                    .push((tri[(k + 1) % 3], tri[(k + 2) % 3]));
            }
        }
        for (v, edges) in &links {
            let successor: HashMap<usize, usize> = edges.iter().copied().collect();
            let heads: HashSet<usize> = edges.iter().map(|&(_, head)| head).collect();
            let mut path_starts = edges
                .iter()
                .map(|&(tail, _)| tail)
                .filter(|tail| !heads.contains(tail));
            let first = path_starts.next();
            if path_starts.next().is_some() {
                return Err(Error::invalid_mesh(format!(
                    "vertex {v} is non-manifold: its faces form several separate fans"
                )));
            }
            let start = first.unwrap_or(edges[0].0);
            let mut walked = 0usize;
            let mut cur = start;
            while let Some(&next) = successor.get(&cur) {
                walked += 1;
                cur = next;
                if cur == start || walked > edges.len() {
                    break;
                }
            }
            if walked != edges.len() {
                return Err(Error::invalid_mesh(format!(
                    "vertex {v} is non-manifold: its faces form several separate fans"
                )));
            }
        }
        Ok(())
    }

    /// Directed-edge index of the surface.
    pub(super) fn directed_edges(&self) -> DirectedEdges {
        let mut map = HashMap::with_capacity(self.triangles.len() * 3);
        for (ti, tri) in self.triangles.iter().enumerate() {
            for k in 0..3 {
                map.insert((tri[k], tri[(k + 1) % 3]), (ti, k));
            }
        }
        map
    }

    /// Set of undirected edges as `(min, max)` pairs.
    pub(super) fn undirected_edges(&self) -> BTreeSet<(usize, usize)> {
        let mut edges = BTreeSet::new();
        for tri in &self.triangles {
            for k in 0..3 {
                let (u, v) = (tri[k], tri[(k + 1) % 3]);
                edges.insert((u.min(v), u.max(v)));
            }
        }
        edges
    }

    /// Number of vertices that belong to at least one face.
    pub(super) fn referenced_vertex_count(&self) -> usize {
        self.triangles
            .iter()
            .flatten()
            .copied()
            .collect::<HashSet<usize>>()
            .len()
    }

    /// Euler characteristic `V - E + F` over the referenced vertices.
    pub(super) fn euler_characteristic(&self) -> i64 {
        self.referenced_vertex_count() as i64 - self.undirected_edges().len() as i64
            + self.triangles.len() as i64
    }

    /// Connected-component label (smallest vertex index of the component) of
    /// every vertex; vertices without faces are their own components.
    pub(super) fn component_labels(&self) -> Vec<usize> {
        let mut sets = DisjointSet::new(self.vertices.len());
        for tri in &self.triangles {
            sets.union(tri[0], tri[1]);
            sets.union(tri[1], tri[2]);
        }
        (0..self.vertices.len()).map(|v| sets.find(v)).collect()
    }

    /// Number of connected components that contain at least one face.
    pub(super) fn component_count(&self) -> usize {
        let labels = self.component_labels();
        self.triangles
            .iter()
            .map(|tri| labels[tri[0]])
            .collect::<HashSet<usize>>()
            .len()
    }

    /// Boundary loops as vertex cycles, each in boundary direction (the surface
    /// lies to the left of every boundary edge `cycle[i] -> cycle[i + 1]`).
    ///
    /// The loops are ordered by the `(triangle, corner)` of their first boundary
    /// edge, and each loop starts at the tail of that edge, so the result is
    /// deterministic. The surface must be a valid manifold (one outgoing
    /// boundary edge per boundary vertex).
    pub(super) fn boundary_loops(&self) -> Vec<Vec<usize>> {
        let directed = self.directed_edges();
        let next_boundary: HashMap<usize, usize> = directed
            .keys()
            .filter(|&&(u, v)| !directed.contains_key(&(v, u)))
            .copied()
            .collect();
        let mut visited: HashSet<usize> = HashSet::new();
        let mut loops = Vec::new();
        for tri in &self.triangles {
            for k in 0..3 {
                let (u, v) = (tri[k], tri[(k + 1) % 3]);
                if next_boundary.get(&u) != Some(&v) || visited.contains(&u) {
                    continue;
                }
                let mut cycle = vec![u];
                visited.insert(u);
                let mut cur = v;
                while cur != u {
                    cycle.push(cur);
                    visited.insert(cur);
                    match next_boundary.get(&cur) {
                        Some(&next) => cur = next,
                        None => break,
                    }
                }
                loops.push(cycle);
            }
        }
        loops
    }

    /// Whether vertex `v` lies on a boundary edge.
    pub(super) fn is_boundary_vertex(&self, directed: &DirectedEdges, v: usize) -> bool {
        directed
            .keys()
            .any(|&(a, b)| (a == v || b == v) && !directed.contains_key(&(b, a)))
    }

    /// Neighbour vertices of `v` (adjacent through an edge).
    pub(super) fn neighbours(&self, v: usize) -> BTreeSet<usize> {
        let mut result = BTreeSet::new();
        for tri in &self.triangles {
            if tri.contains(&v) {
                result.extend(tri.iter().copied().filter(|&w| w != v));
            }
        }
        result
    }

    /// Check the link condition (Dey, Edelsbrunner, Guha and Nekhayev, "Topology
    /// preserving edge contraction") for collapsing the edge `a -- b`.
    ///
    /// Collapsing preserves the topology of a triangulated 2-manifold exactly when
    /// the vertices adjacent to both endpoints are the vertices opposite the edge
    /// in its one or two faces, where every boundary loop is closed by a vertex at
    /// infinity (so two boundary vertices may only be joined by a boundary edge).
    /// The test additionally rejects a collapse that would merge two faces onto
    /// the same vertex set (for example any edge of a tetrahedron, or an edge
    /// whose opposite vertex has valence three).
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] explaining which condition is violated.
    pub(super) fn check_collapse(&self, a: usize, b: usize) -> Result<()> {
        let directed = self.directed_edges();
        let mut link: BTreeSet<usize> = BTreeSet::new();
        let mut faces_on_edge = 0usize;
        for tri in &self.triangles {
            if tri.contains(&a) && tri.contains(&b) {
                faces_on_edge += 1;
                link.extend(tri.iter().copied().filter(|&w| w != a && w != b));
            }
        }
        if faces_on_edge == 0 {
            return Err(Error::invalid_mesh(format!(
                "no face contains both {a} and {b}"
            )));
        }
        if faces_on_edge == 1 {
            link.insert(INFINITE_VERTEX);
        }
        let mut common: BTreeSet<usize> = self
            .neighbours(a)
            .intersection(&self.neighbours(b))
            .copied()
            .collect();
        if self.is_boundary_vertex(&directed, a) && self.is_boundary_vertex(&directed, b) {
            common.insert(INFINITE_VERTEX);
        }
        if common != link {
            let extra: Vec<usize> = common.difference(&link).copied().collect();
            let detail = if extra == [INFINITE_VERTEX] {
                "the boundary (both endpoints are boundary vertices but the edge is interior)"
                    .to_string()
            } else {
                format!("the neighbours {extra:?}")
            };
            return Err(Error::unsupported(
                "edge_collapse",
                format!(
                    "link condition violated for edge {a}-{b}: the endpoints also share {detail}"
                ),
            ));
        }
        let faces: HashSet<[usize; 3]> = self
            .triangles
            .iter()
            .map(|tri| {
                let mut key = *tri;
                key.sort_unstable();
                key
            })
            .collect();
        for tri in &self.triangles {
            if tri.contains(&a) && !tri.contains(&b) {
                let mut merged = [0usize; 3];
                for (slot, &v) in merged.iter_mut().zip(tri.iter()) {
                    *slot = if v == a { b } else { v };
                }
                merged.sort_unstable();
                if faces.contains(&merged) {
                    return Err(Error::unsupported(
                        "edge_collapse",
                        format!(
                            "collapsing edge {a}-{b} would create two faces on the vertex set {merged:?}"
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Euclidean length of the edge `u -- v` (`f64::MAX` if it is not finite).
    pub(super) fn edge_length(&self, u: usize, v: usize) -> f64 {
        let length = len3(sub3(self.vertices[u], self.vertices[v]));
        if length.is_finite() { length } else { f64::MAX }
    }
}
