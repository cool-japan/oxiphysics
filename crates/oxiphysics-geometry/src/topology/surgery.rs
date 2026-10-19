//! Topological surgery on a half-edge mesh: edge split, edge collapse, handle
//! attachment and handle deletion.
//!
//! Every operation works on the plain `vertices` / `triangles` form of the
//! surface (the `surface` module): the half-edge mesh is read into that form,
//! edited, validated as an oriented 2-manifold, and only then rebuilt. An
//! operation that fails therefore leaves the mesh exactly as it was.

use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet, VecDeque};

use super::functions::{add3, scale3};
use super::surface::{DirectedEdges, HalfEdgeLocator, Surface};
use super::types::HalfEdgeMesh;
use crate::error::{Error, Result};

/// Topological surgery operations on a half-edge mesh.
///
/// The mesh must be a valid oriented triangle 2-manifold: every face a
/// triangle, every directed edge used once, and the faces around every vertex
/// forming a single fan. An operation reports [`Error::InvalidMesh`] otherwise,
/// and [`Error::Unsupported`] when the request is well-formed but the surgery
/// cannot be performed (for example a collapse that violates the link
/// condition). On any error the mesh is left untouched.
///
/// After a successful operation the half-edge arrays are rebuilt, so half-edge
/// and face indices obtained before the call are no longer meaningful; the
/// documentation of each operation states what happens to vertex indices.
pub struct TopologicalSurgery<'a> {
    /// Reference to the mesh being operated on.
    pub mesh: &'a mut HalfEdgeMesh,
}

/// Face, endpoints and opposite vertex of one half-edge.
struct HalfEdgeRef {
    face: usize,
    origin: usize,
    destination: usize,
    opposite: usize,
}

impl<'a> TopologicalSurgery<'a> {
    /// Create a surgery operator.
    pub fn new(mesh: &'a mut HalfEdgeMesh) -> Self {
        Self { mesh }
    }

    /// Boundary loops of the mesh as vertex cycles.
    ///
    /// Each loop lists its vertices in boundary direction (the surface lies to
    /// the left of every edge `loop[i] -> loop[i + 1]`). Loops are ordered by
    /// the first boundary edge that occurs when the faces are scanned in index
    /// order, and each loop starts at the origin of that edge; the loop indices
    /// accepted by [`TopologicalSurgery::attach_handle`] refer to this order.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidMesh`] if the mesh is not a valid oriented 2-manifold.
    pub fn boundary_loops(&self) -> Result<Vec<Vec<usize>>> {
        let (surface, _) = Surface::from_mesh(self.mesh)?;
        Ok(surface.boundary_loops())
    }

    /// Attach a handle: join two boundary loops of the same connected component
    /// with a tube of triangles, increasing the genus of that component by one.
    ///
    /// `boundary_a` and `boundary_b` index [`TopologicalSurgery::boundary_loops`].
    /// The tube pairs the vertices of the two loops in opposite cyclic order (so
    /// that the result stays orientable) and is triangulated as a strip of
    /// `len_a + len_b` triangles with no new vertices; among all cyclic
    /// alignments and two strip rules (greedy shortest rung, and balanced
    /// advance) the valid strip with the shortest total rung length is used.
    /// Loops may have different lengths. The Euler characteristic is unchanged
    /// and the number of boundary loops drops by two, so a sphere with two holes
    /// becomes a torus. Vertex indices do not change.
    ///
    /// # Errors
    ///
    /// * [`Error::IndexOutOfBounds`] if a loop index does not exist,
    /// * [`Error::Unsupported`] if both indices are equal, if the loops lie in
    ///   different connected components (joining them would merge components
    ///   instead of adding a handle), or if every possible strip would duplicate
    ///   an existing edge,
    /// * [`Error::InvalidMesh`] if the mesh is not a valid oriented 2-manifold.
    pub fn attach_handle(&mut self, boundary_a: usize, boundary_b: usize) -> Result<()> {
        let (mut surface, _) = Surface::from_mesh(self.mesh)?;
        let loops = surface.boundary_loops();
        for index in [boundary_a, boundary_b] {
            if index >= loops.len() {
                return Err(Error::index_out_of_bounds(index, loops.len()));
            }
        }
        if boundary_a == boundary_b {
            return Err(Error::unsupported(
                "attach_handle",
                "a handle joins two distinct boundary loops",
            ));
        }
        let loop_a = &loops[boundary_a];
        let loop_b = &loops[boundary_b];
        let labels = surface.component_labels();
        if labels[loop_a[0]] != labels[loop_b[0]] {
            return Err(Error::unsupported(
                "attach_handle",
                "the two boundary loops lie in different connected components",
            ));
        }
        let band = best_band(&surface, loop_a, loop_b)?;
        surface.triangles.extend(band);
        surface.validate()?;
        *self.mesh = surface.into_mesh();
        Ok(())
    }

    /// Delete a handle: cut a closed component of positive genus along a
    /// non-separating loop of edges, then cap the two resulting holes with
    /// fans, decreasing the genus of that component by one.
    ///
    /// The first closed component (by lowest vertex index) with Euler
    /// characteristic `<= 0` is used. The loop comes from a tree-cotree
    /// decomposition (Eppstein 2003): a shortest-path spanning tree of the
    /// vertex graph, a spanning tree of the dual graph among the remaining edges,
    /// and for every leftover edge the fundamental cycle it closes in the tree.
    /// Such cycles generate the first homology of the surface; the shortest one
    /// whose removal keeps the surface connected is cut. Every loop vertex is
    /// duplicated (the faces on one side of the loop move to the copy) and each of
    /// the two holes receives a cone vertex at the centroid of its boundary.
    /// Existing vertex indices are preserved; the loop-vertex copies and the two
    /// cone vertices are appended. The Euler characteristic rises by two.
    ///
    /// # Errors
    ///
    /// * [`Error::Unsupported`] if no closed component has positive genus or no
    ///   loop produces a valid result,
    /// * [`Error::InvalidMesh`] if the mesh is not a valid oriented 2-manifold.
    pub fn delete_handle(&mut self) -> Result<()> {
        let (surface, _) = Surface::from_mesh(self.mesh)?;
        let directed = surface.directed_edges();
        let labels = surface.component_labels();
        let component = first_closed_component_with_handle(&surface, &directed, &labels)
            .ok_or_else(|| {
                Error::unsupported(
                    "delete_handle",
                    "no closed connected component has a handle (genus > 0)",
                )
            })?;
        let expected_components = surface.component_count();
        let expected_euler = surface.euler_characteristic() + 2;
        for cycle in handle_loop_candidates(&surface, &labels, component) {
            let Ok(result) = cut_and_cap(&surface, &directed, &cycle) else {
                continue;
            };
            if result.component_count() == expected_components
                && result.euler_characteristic() == expected_euler
            {
                *self.mesh = result.into_mesh();
                return Ok(());
            }
        }
        Err(Error::unsupported(
            "delete_handle",
            "no non-separating loop of edges could be cut and capped",
        ))
    }

    /// Edge collapse: contract the edge of half-edge `he_idx` (from `v_start` to
    /// `v_end`) into one vertex, placed at the edge midpoint.
    ///
    /// The two faces incident to the edge (one for a boundary edge) disappear and
    /// the edge's three edges per face merge in pairs, so `V`, `E` and `F` drop by
    /// `1`, `3` and `2` (`1`, `2` and `1` on the boundary) and the Euler
    /// characteristic is preserved. The collapse is only performed when it keeps
    /// the surface a 2-manifold of the same topology: the *link condition* must
    /// hold (the endpoints' common neighbours are exactly the vertices opposite
    /// the edge, boundary loops being closed by a vertex at infinity) and no two
    /// faces may end up on the same vertex set (so no edge of a tetrahedron, and
    /// no edge next to a valence-3 vertex, can be collapsed).
    ///
    /// `v_start` is removed from the vertex array, so every vertex index above it
    /// shifts down by one. Returns the new index of the surviving vertex
    /// (`v_end`, or `v_end - 1` when `v_end > v_start`).
    ///
    /// # Errors
    ///
    /// * [`Error::IndexOutOfBounds`] if `he_idx` is not a half-edge of the mesh,
    /// * [`Error::Unsupported`] if the link condition is violated,
    /// * [`Error::InvalidMesh`] if the mesh is not a valid oriented 2-manifold.
    pub fn edge_collapse(&mut self, he_idx: usize) -> Result<usize> {
        let (mut surface, locator) = Surface::from_mesh(self.mesh)?;
        let edge = half_edge_ref(self.mesh, &locator, &surface, he_idx)?;
        let (removed, survivor) = (edge.origin, edge.destination);
        surface.check_collapse(removed, survivor)?;
        let midpoint = scale3(
            add3(surface.vertices[removed], surface.vertices[survivor]),
            0.5,
        );
        surface
            .triangles
            .retain(|tri| !(tri.contains(&removed) && tri.contains(&survivor)));
        for corner in surface.triangles.iter_mut().flatten() {
            if *corner == removed {
                *corner = survivor;
            }
            if *corner > removed {
                *corner -= 1;
            }
        }
        surface.vertices[survivor] = midpoint;
        surface.vertices.remove(removed);
        surface.validate()?;
        *self.mesh = surface.into_mesh();
        Ok(if survivor > removed {
            survivor - 1
        } else {
            survivor
        })
    }

    /// Edge split: insert a new vertex at the midpoint of the edge of half-edge
    /// `he_idx` and split each incident triangle in two, so that `V`, `E` and `F`
    /// grow by `1`, `3` and `2` for an interior edge (`1`, `2` and `1` on the
    /// boundary) and the Euler characteristic is preserved.
    ///
    /// Existing vertex and face indices keep their meaning (the split faces are
    /// rewritten in place and the new faces appended). Returns the index of the
    /// new vertex.
    ///
    /// # Errors
    ///
    /// * [`Error::IndexOutOfBounds`] if `he_idx` is not a half-edge of the mesh,
    /// * [`Error::InvalidMesh`] if the mesh is not a valid oriented 2-manifold.
    pub fn edge_split(&mut self, he_idx: usize) -> Result<usize> {
        let (mut surface, locator) = Surface::from_mesh(self.mesh)?;
        let edge = half_edge_ref(self.mesh, &locator, &surface, he_idx)?;
        let (a, b, c) = (edge.origin, edge.destination, edge.opposite);
        let directed = surface.directed_edges();
        let new_vertex = surface.vertices.len();
        surface
            .vertices
            .push(scale3(add3(surface.vertices[a], surface.vertices[b]), 0.5));
        surface.triangles[edge.face] = [a, new_vertex, c];
        surface.triangles.push([new_vertex, b, c]);
        if let Some(&(twin_face, corner)) = directed.get(&(b, a)) {
            let d = surface.triangles[twin_face][(corner + 2) % 3];
            surface.triangles[twin_face] = [b, new_vertex, d];
            surface.triangles.push([new_vertex, a, d]);
        }
        surface.validate()?;
        *self.mesh = surface.into_mesh();
        Ok(new_vertex)
    }
}

/// Resolve half-edge `he_idx` to its face and vertices.
fn half_edge_ref(
    mesh: &HalfEdgeMesh,
    locator: &HalfEdgeLocator,
    surface: &Surface,
    he_idx: usize,
) -> Result<HalfEdgeRef> {
    if he_idx >= mesh.half_edges.len() {
        return Err(Error::index_out_of_bounds(he_idx, mesh.half_edges.len()));
    }
    let Some((face, corner)) = locator[he_idx] else {
        return Err(Error::invalid_mesh(format!(
            "half-edge {he_idx} does not belong to any face"
        )));
    };
    let tri = surface.triangles[face];
    Ok(HalfEdgeRef {
        face,
        origin: tri[corner],
        destination: tri[(corner + 1) % 3],
        opposite: tri[(corner + 2) % 3],
    })
}

// ---------------------------------------------------------------------------
// Handle attachment: stitching a tube between two boundary loops
// ---------------------------------------------------------------------------

/// How the strip between the two loops advances at each step.
#[derive(Clone, Copy)]
enum StripRule {
    /// Advance along the loop that gives the shorter new rung.
    ShortestRung,
    /// Advance the loop that is proportionally behind, keeping both in step.
    Balanced,
}

/// The cheapest valid strip over all cyclic alignments and strip rules.
fn best_band(surface: &Surface, loop_a: &[usize], loop_b: &[usize]) -> Result<Vec<[usize; 3]>> {
    let existing = surface.undirected_edges();
    let mut best: Option<(f64, Vec<[usize; 3]>)> = None;
    for offset in 0..loop_b.len() {
        for rule in [StripRule::ShortestRung, StripRule::Balanced] {
            let Some((cost, band)) =
                strip_between(surface, loop_a, loop_b, offset, rule, &existing)
            else {
                continue;
            };
            if best.as_ref().is_none_or(|(best_cost, _)| cost < *best_cost) {
                best = Some((cost, band));
            }
        }
    }
    best.map(|(_, band)| band).ok_or_else(|| {
        Error::unsupported(
            "attach_handle",
            "every way of stitching the two loops would duplicate an existing edge",
        )
    })
}

/// Triangulate the tube between `loop_a` and `loop_b`.
///
/// `loop_b` is walked backwards starting at `loop_b[offset]` so that the tube
/// stays orientable. A strip step either advances along loop A with the
/// triangle `(a[i + 1], a[i], b[j])` (whose edge `a[i + 1] -> a[i]` is the twin
/// of boundary edge `a[i] -> a[i + 1]`) or along loop B with `(a[i], b[j],
/// b[j + 1])`; the lattice path from `(0, 0)` to `(len_a, len_b)` yields
/// `len_a + len_b` triangles. Returns the total rung length and the triangles,
/// or `None` if a rung would coincide with an existing edge or with another
/// rung of the strip.
fn strip_between(
    surface: &Surface,
    loop_a: &[usize],
    loop_b: &[usize],
    offset: usize,
    rule: StripRule,
    existing: &BTreeSet<(usize, usize)>,
) -> Option<(f64, Vec<[usize; 3]>)> {
    let (n, m) = (loop_a.len(), loop_b.len());
    let a_at = |i: usize| loop_a[i % n];
    let b_at = |j: usize| loop_b[(offset + m - j % m) % m];
    let rung_length = |i: usize, j: usize| surface.edge_length(a_at(i), b_at(j));
    let mut rungs: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut claim = |i: usize, j: usize| {
        let (u, v) = (a_at(i), b_at(j));
        let key = (u.min(v), u.max(v));
        !existing.contains(&key) && rungs.insert(key)
    };
    if !claim(0, 0) {
        return None;
    }
    let mut cost = rung_length(0, 0);
    let mut band = Vec::with_capacity(n + m);
    let (mut i, mut j) = (0usize, 0usize);
    while i < n || j < m {
        let advance_a = if i == n {
            false
        } else if j == m {
            true
        } else {
            match rule {
                StripRule::ShortestRung => rung_length(i + 1, j) <= rung_length(i, j + 1),
                StripRule::Balanced => i * m <= j * n,
            }
        };
        if advance_a {
            band.push([a_at(i + 1), a_at(i), b_at(j)]);
            i += 1;
        } else {
            band.push([a_at(i), b_at(j), b_at(j + 1)]);
            j += 1;
        }
        if i < n || j < m {
            if !claim(i, j) {
                return None;
            }
            cost += rung_length(i, j);
        }
    }
    Some((cost, band))
}

// ---------------------------------------------------------------------------
// Handle deletion: tree-cotree loop, cut, cap
// ---------------------------------------------------------------------------

/// Smallest vertex label of the first component that is closed and has Euler
/// characteristic `<= 0` (genus >= 1).
fn first_closed_component_with_handle(
    surface: &Surface,
    directed: &DirectedEdges,
    labels: &[usize],
) -> Option<usize> {
    let bordered: HashSet<usize> = directed
        .keys()
        .filter(|&&(u, v)| !directed.contains_key(&(v, u)))
        .map(|&(u, _)| labels[u])
        .collect();
    let mut vertices: BTreeMap<usize, HashSet<usize>> = BTreeMap::new();
    let mut edges: BTreeMap<usize, HashSet<(usize, usize)>> = BTreeMap::new();
    let mut faces: BTreeMap<usize, usize> = BTreeMap::new();
    for tri in &surface.triangles {
        let label = labels[tri[0]];
        vertices
            .entry(label)
            .or_default()
            .extend(tri.iter().copied());
        *faces.entry(label).or_default() += 1;
        for k in 0..3 {
            let (u, v) = (tri[k], tri[(k + 1) % 3]);
            edges.entry(label).or_default().insert((u.min(v), u.max(v)));
        }
    }
    faces.keys().copied().find(|label| {
        if bordered.contains(label) {
            return false;
        }
        let v = vertices.get(label).map_or(0, HashSet::len) as i64;
        let e = edges.get(label).map_or(0, HashSet::len) as i64;
        let f = faces.get(label).copied().unwrap_or(0) as i64;
        v - e + f <= 0
    })
}

/// Candidate handle loops of one closed component, shortest first.
///
/// Builds a shortest-path spanning tree `T` of the vertex graph rooted at the
/// lowest vertex, then a spanning tree `C` of the dual graph that crosses only
/// edges outside `T`. The `2g` edges in neither tree each close a simple cycle
/// with their tree path; these cycles are a basis of the first homology.
fn handle_loop_candidates(
    surface: &Surface,
    labels: &[usize],
    component: usize,
) -> Vec<Vec<usize>> {
    let faces: Vec<usize> = (0..surface.triangles.len())
        .filter(|&t| labels[surface.triangles[t][0]] == component)
        .collect();
    let mut adjacency: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    let mut edge_faces: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for &t in &faces {
        let tri = surface.triangles[t];
        for k in 0..3 {
            let (u, v) = (tri[k], tri[(k + 1) % 3]);
            adjacency.entry(u).or_default().insert(v);
            adjacency.entry(v).or_default().insert(u);
            edge_faces.entry((u.min(v), u.max(v))).or_default().push(t);
        }
    }
    let Some(&root) = adjacency.keys().next() else {
        return Vec::new();
    };
    let parent = shortest_path_tree(surface, &adjacency, root);
    let tree_edges: BTreeSet<(usize, usize)> = parent
        .iter()
        .map(|(&child, &up)| (child.min(up), child.max(up)))
        .collect();

    // Spanning tree of the dual graph that never crosses a primal tree edge.
    let mut cotree_edges: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut visited: HashSet<usize> = HashSet::new();
    let mut queue: VecDeque<usize> = VecDeque::new();
    if let Some(&first) = faces.first() {
        visited.insert(first);
        queue.push_back(first);
    }
    while let Some(t) = queue.pop_front() {
        let tri = surface.triangles[t];
        for k in 0..3 {
            let (u, v) = (tri[k], tri[(k + 1) % 3]);
            let key = (u.min(v), u.max(v));
            if tree_edges.contains(&key) {
                continue;
            }
            for &other in edge_faces.get(&key).into_iter().flatten() {
                if visited.insert(other) {
                    cotree_edges.insert(key);
                    queue.push_back(other);
                }
            }
        }
    }

    let mut candidates: Vec<(f64, (usize, usize), Vec<usize>)> = edge_faces
        .keys()
        .filter(|&&key| !tree_edges.contains(&key) && !cotree_edges.contains(&key))
        .filter_map(|&(u, v)| {
            let cycle = fundamental_cycle(&parent, u, v);
            if cycle.len() < 3 {
                return None;
            }
            let length: f64 = (0..cycle.len())
                .map(|i| surface.edge_length(cycle[i], cycle[(i + 1) % cycle.len()]))
                .sum();
            Some((length, (u, v), cycle))
        })
        .collect();
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    candidates.into_iter().map(|(_, _, cycle)| cycle).collect()
}

/// Dijkstra shortest-path tree over the vertex graph; returns the parent of
/// every non-root vertex. Ties are broken by vertex index, so the tree is
/// deterministic.
fn shortest_path_tree(
    surface: &Surface,
    adjacency: &BTreeMap<usize, BTreeSet<usize>>,
    root: usize,
) -> HashMap<usize, usize> {
    let mut distance: HashMap<usize, f64> = HashMap::new();
    let mut parent: HashMap<usize, usize> = HashMap::new();
    // Non-negative f64 bit patterns order like the numbers themselves.
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    distance.insert(root, 0.0);
    heap.push(Reverse((0.0f64.to_bits(), root)));
    while let Some(Reverse((bits, u))) = heap.pop() {
        let d = f64::from_bits(bits);
        if d > distance.get(&u).copied().unwrap_or(f64::INFINITY) {
            continue;
        }
        for &w in adjacency.get(&u).into_iter().flatten() {
            let candidate = d + surface.edge_length(u, w);
            if candidate < distance.get(&w).copied().unwrap_or(f64::INFINITY) {
                distance.insert(w, candidate);
                parent.insert(w, u);
                heap.push(Reverse((candidate.to_bits(), w)));
            }
        }
    }
    parent
}

/// The simple cycle `u -> .. -> lca -> .. -> v` formed by the tree path between
/// `u` and `v` (the closing edge `v -- u` is implicit).
fn fundamental_cycle(parent: &HashMap<usize, usize>, u: usize, v: usize) -> Vec<usize> {
    let climb = |start: usize| {
        let mut path = vec![start];
        let mut cur = start;
        while let Some(&up) = parent.get(&cur) {
            path.push(up);
            cur = up;
        }
        path
    };
    let from_u = climb(u);
    let from_v = climb(v);
    let position_in_v: HashMap<usize, usize> =
        from_v.iter().enumerate().map(|(i, &x)| (x, i)).collect();
    for (iu, x) in from_u.iter().enumerate() {
        if let Some(&iv) = position_in_v.get(x) {
            let mut cycle: Vec<usize> = from_u[..=iu].to_vec();
            cycle.extend(from_v[..iv].iter().rev());
            return cycle;
        }
    }
    Vec::new()
}

/// Cut the surface along the closed loop `cycle` and cap both new holes.
///
/// Rotating counter-clockwise around a loop vertex from its predecessor to its
/// successor sweeps the faces on the right of the loop; those faces are moved to
/// a copy of the vertex (a face may use the copy at one loop vertex and the
/// original at another, which is exactly what a non-separating cut needs). The
/// two boundary loops of the cut (the originals on the left and the copies on
/// the right) are each closed with a cone vertex at their centroid.
fn cut_and_cap(surface: &Surface, directed: &DirectedEdges, cycle: &[usize]) -> Result<Surface> {
    let k = cycle.len();
    let base = surface.vertices.len();
    let mut right_side: HashMap<(usize, usize), usize> = HashMap::new();
    for (i, &v) in cycle.iter().enumerate() {
        let previous = cycle[(i + k - 1) % k];
        let next = cycle[(i + 1) % k];
        let mut around = previous;
        let mut steps = 0usize;
        loop {
            let &(face, corner) = directed.get(&(v, around)).ok_or_else(|| {
                Error::unsupported("delete_handle", "the loop passes through a boundary vertex")
            })?;
            right_side.insert((face, v), base + i);
            let following = surface.triangles[face][(corner + 2) % 3];
            if following == next {
                break;
            }
            around = following;
            steps += 1;
            if steps > surface.triangles.len() {
                return Err(Error::invalid_mesh(format!(
                    "the faces around vertex {v} do not form a closed fan"
                )));
            }
        }
    }
    let mut triangles = surface.triangles.clone();
    for (face, tri) in triangles.iter_mut().enumerate() {
        for corner in tri.iter_mut() {
            if let Some(&copy) = right_side.get(&(face, *corner)) {
                *corner = copy;
            }
        }
    }
    let mut vertices = surface.vertices.clone();
    vertices.extend(cycle.iter().map(|&v| surface.vertices[v]));
    let cut = Surface {
        vertices,
        triangles,
    };
    cut.validate()?;

    let touched: HashSet<usize> = cycle.iter().copied().chain(base..base + k).collect();
    let holes: Vec<Vec<usize>> = cut
        .boundary_loops()
        .into_iter()
        .filter(|hole| hole.iter().any(|v| touched.contains(v)))
        .collect();
    if holes.len() != 2 {
        return Err(Error::unsupported(
            "delete_handle",
            format!(
                "cutting along the loop produced {} holes instead of 2",
                holes.len()
            ),
        ));
    }
    let mut capped = cut;
    for hole in &holes {
        let apex = capped.vertices.len();
        let sum = hole
            .iter()
            .fold([0.0f64; 3], |acc, &v| add3(acc, capped.vertices[v]));
        capped.vertices.push(scale3(sum, 1.0 / hole.len() as f64));
        for (idx, &u) in hole.iter().enumerate() {
            let v = hole[(idx + 1) % hole.len()];
            capped.triangles.push([v, u, apex]);
        }
    }
    capped.validate()?;
    Ok(capped)
}
