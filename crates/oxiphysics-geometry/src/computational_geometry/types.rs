//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{FaceId, HalfEdgeId, Point3, VertexId};
#[allow(unused_imports)]
use super::functions_2::*;

/// A Voronoi cell: the site point and its circumcenter-based Voronoi vertices.
#[derive(Debug, Clone)]
pub struct VoronoiCell2D {
    /// The generating site index.
    pub site: usize,
    /// Circumcenters of Delaunay triangles that form the Voronoi polygon vertices.
    pub circumcenters: Vec<Point2>,
}
/// A face in a line arrangement (a convex polygon cell bounded by arrangement edges).
#[derive(Debug, Clone)]
pub struct ArrangementFace {
    /// Indices of bounding arrangement vertices (in order).
    pub vertex_indices: Vec<usize>,
}
/// An arrangement of lines in the plane.
///
/// Stores the vertices (pairwise intersections) and adjacency information.
#[derive(Debug, Clone, Default)]
pub struct LineArrangement {
    /// The input lines.
    pub lines: Vec<Line2D>,
    /// All vertices of the arrangement (pairwise intersections).
    pub vertices: Vec<ArrangementVertex>,
}
impl LineArrangement {
    /// Construct the arrangement from a set of lines.
    ///
    /// Computes all O(n²) pairwise intersections, merging near-coincident ones.
    pub fn build(lines: &[Line2D]) -> Self {
        let n = lines.len();
        let mut vertices: Vec<ArrangementVertex> = Vec::new();
        let merge_eps = 1e-9;
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(pt) = line_intersect_2d(&lines[i], &lines[j]) {
                    let existing = vertices
                        .iter_mut()
                        .find(|v| v.point.dist_sq(pt) < merge_eps);
                    if let Some(v) = existing {
                        if !v.line_indices.contains(&i) {
                            v.line_indices.push(i);
                        }
                        if !v.line_indices.contains(&j) {
                            v.line_indices.push(j);
                        }
                    } else {
                        vertices.push(ArrangementVertex {
                            point: pt,
                            line_indices: vec![i, j],
                        });
                    }
                }
            }
        }
        Self {
            lines: lines.to_vec(),
            vertices,
        }
    }
    /// Number of vertices (intersection points).
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// By Euler's formula V - E + F = 2, compute the number of faces
    /// (including the unbounded outer face) from V and E.
    ///
    /// For n lines in general position: V = n(n-1)/2, E = n(n+1), F = n(n-1)/2 + n + 1.
    pub fn euler_face_count(&self) -> usize {
        let n = self.lines.len();
        n * (n - 1) / 2 + n + 1
    }
    /// Sort the vertices along a given line by arc-length parameter.
    pub fn vertices_on_line(&self, line_idx: usize) -> Vec<usize> {
        let mut idxs: Vec<usize> = self
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.line_indices.contains(&line_idx))
            .map(|(i, _)| i)
            .collect();
        let l = &self.lines[line_idx];
        let dir = Point2::new(-l.b, l.a);
        idxs.sort_by(|&a, &b| {
            let ta = dir.dot(self.vertices[a].point);
            let tb = dir.dot(self.vertices[b].point);
            ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
        });
        idxs
    }
}
/// A face in the half-edge mesh, storing one of its bounding half-edges.
#[derive(Debug, Clone)]
pub struct HEFace {
    /// One half-edge on the boundary of this face.
    pub half_edge: HalfEdgeId,
}
/// A single directed half-edge.
#[derive(Debug, Clone)]
pub struct HalfEdge {
    /// Origin vertex of this half-edge.
    pub origin: VertexId,
    /// The twin (opposite) half-edge.
    pub twin: HalfEdgeId,
    /// The next half-edge around the face (counter-clockwise).
    pub next: HalfEdgeId,
    /// The previous half-edge around the face.
    pub prev: HalfEdgeId,
    /// The face to the left of this half-edge.
    pub face: Option<FaceId>,
}
/// A visibility graph for path planning among 2D polygon obstacles.
///
/// Nodes are the vertices of all obstacles plus optional source/target points.
/// Edges connect pairs of nodes that can see each other (the open segment
/// between them does not cross any obstacle edge).
#[derive(Debug, Clone, Default)]
pub struct VisibilityGraph {
    /// All nodes (vertices + queries).
    pub nodes: Vec<Point2>,
    /// Adjacency list: `edges[i]` = list of node indices visible from node `i`.
    pub edges: Vec<Vec<usize>>,
}
impl VisibilityGraph {
    /// Build the visibility graph from a collection of obstacles.
    pub fn build(obstacles: &[Obstacle2D]) -> Self {
        let mut nodes: Vec<Point2> = Vec::new();
        for obs in obstacles {
            for &v in &obs.vertices {
                nodes.push(v);
            }
        }
        let all_edges: Vec<(Point2, Point2)> = obstacles.iter().flat_map(|o| o.edges()).collect();
        let n = nodes.len();
        let mut edges = vec![Vec::new(); n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                if segment_visible(nodes[i], nodes[j], &all_edges) {
                    edges[i].push(j);
                }
            }
        }
        Self { nodes, edges }
    }
    /// Add a query point (e.g., start or goal) and compute its visibility to existing nodes.
    pub fn add_query_point(&mut self, pt: Point2, obstacles: &[Obstacle2D]) {
        let all_edges: Vec<(Point2, Point2)> = obstacles.iter().flat_map(|o| o.edges()).collect();
        let n = self.nodes.len();
        let new_id = n;
        self.nodes.push(pt);
        self.edges.push(Vec::new());
        for i in 0..n {
            if segment_visible(pt, self.nodes[i], &all_edges) {
                self.edges[new_id].push(i);
                self.edges[i].push(new_id);
            }
        }
    }
    /// Number of nodes in the graph.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Shortest path (by Euclidean distance) from `src` to `dst` using Dijkstra.
    pub fn shortest_path(&self, src: usize, dst: usize) -> Option<Vec<usize>> {
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;
        let n = self.nodes.len();
        let mut dist = vec![f64::INFINITY; n];
        let mut prev = vec![usize::MAX; n];
        dist[src] = 0.0;
        let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
        heap.push(Reverse((0, src)));
        while let Some(Reverse((d_raw, u))) = heap.pop() {
            let d = d_raw as f64 / 1e9;
            if d > dist[u] + 1e-12 {
                continue;
            }
            if u == dst {
                break;
            }
            for &v in &self.edges[u] {
                let nd = dist[u] + self.nodes[u].dist(self.nodes[v]);
                if nd < dist[v] - 1e-12 {
                    dist[v] = nd;
                    prev[v] = u;
                    heap.push(Reverse(((nd * 1e9) as u64, v)));
                }
            }
        }
        if dist[dst].is_infinite() {
            return None;
        }
        let mut path = Vec::new();
        let mut cur = dst;
        while cur != usize::MAX {
            path.push(cur);
            cur = prev[cur];
        }
        path.reverse();
        Some(path)
    }
}
/// A vertex in the half-edge mesh, storing a representative outgoing half-edge.
#[derive(Debug, Clone)]
pub struct HEVertex {
    /// Position of this vertex.
    pub position: Point3,
    /// An outgoing half-edge from this vertex.
    pub half_edge: Option<HalfEdgeId>,
}
/// A line in the plane represented in the form ax + by = c.
#[derive(Debug, Clone, Copy)]
pub struct Line2D {
    /// Coefficient a.
    pub a: f64,
    /// Coefficient b.
    pub b: f64,
    /// Right-hand side c.
    pub c: f64,
}
impl Line2D {
    /// Construct a line through two points.
    pub fn through(p: Point2, q: Point2) -> Self {
        let a = q.y - p.y;
        let b = p.x - q.x;
        let c = a * p.x + b * p.y;
        Self { a, b, c }
    }
    /// Signed distance of point p from this line (positive on the a·x+b·y>c side).
    pub fn signed_dist(&self, p: Point2) -> f64 {
        let norm = (self.a * self.a + self.b * self.b).sqrt();
        if norm < 1e-15 {
            return 0.0;
        }
        (self.a * p.x + self.b * p.y - self.c) / norm
    }
    /// Evaluate which side of the line the point is on.
    /// Returns positive if on the left, negative if on the right, zero if on the line.
    pub fn side(&self, p: Point2) -> f64 {
        self.a * p.x + self.b * p.y - self.c
    }
}
/// A trapezoid in the trapezoidal decomposition of a planar subdivision.
///
/// Bounded by left and right endpoints and top/bottom line segments.
#[derive(Debug, Clone)]
pub struct Trapezoid {
    /// Left x-boundary.
    pub x_left: f64,
    /// Right x-boundary.
    pub x_right: f64,
    /// Bottom edge index (index into segment list, or `usize::MAX` for floor).
    pub bottom_seg: usize,
    /// Top edge index (index into segment list, or `usize::MAX` for ceiling).
    pub top_seg: usize,
    /// Face label (e.g., triangle index).
    pub face: Option<usize>,
}
/// A triangle in a 2D Delaunay triangulation (vertex indices into a point list).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DelaunayTri {
    /// First vertex index.
    pub a: usize,
    /// Second vertex index.
    pub b: usize,
    /// Third vertex index.
    pub c: usize,
}
impl DelaunayTri {
    /// Create a new Delaunay triangle from vertex indices.
    pub fn new(a: usize, b: usize, c: usize) -> Self {
        Self { a, b, c }
    }
}
/// A line segment in 2D for slab-based point location.
#[derive(Debug, Clone, Copy)]
pub struct Segment2D {
    /// Left endpoint.
    pub left: Point2,
    /// Right endpoint.
    pub right: Point2,
}
impl Segment2D {
    /// Create a new segment, normalising so `left.x <= right.x`.
    pub fn new(a: Point2, b: Point2) -> Self {
        if a.x <= b.x {
            Self { left: a, right: b }
        } else {
            Self { left: b, right: a }
        }
    }
    /// Y-value on the segment at x (by linear interpolation).
    pub fn y_at(&self, x: f64) -> Option<f64> {
        let dx = self.right.x - self.left.x;
        if dx.abs() < 1e-15 {
            return None;
        }
        let t = (x - self.left.x) / dx;
        if !(0.0..=1.0).contains(&t) {
            return None;
        }
        Some(self.left.y + t * (self.right.y - self.left.y))
    }
}
/// Slab-based point location structure.
///
/// Partitions the plane into vertical slabs at the x-coordinates of all
/// segment endpoints, then stores the ordered set of segments in each slab.
/// Query time: O(log n) slabs + O(log k) binary search within a slab.
#[derive(Debug, Clone, Default)]
pub struct SlabPointLocator {
    /// Sorted distinct x-coordinates of all endpoints.
    pub slab_xs: Vec<f64>,
    /// For each slab i (between slab_xs\[i\] and slab_xs\[i+1\]),
    /// the list of segment indices that cross that slab, sorted top-to-bottom
    /// at the slab midpoint.
    pub slab_segments: Vec<Vec<usize>>,
    /// The segments.
    pub segments: Vec<Segment2D>,
    /// Face labels per trapezoid cell: slab_faces\[slab\]\[k\] is the face below
    /// segments\[slab_segments\[slab\\]\[k\]].
    pub slab_faces: Vec<Vec<Option<usize>>>,
}
impl SlabPointLocator {
    /// Build the slab structure from a list of segments and optional face labels.
    ///
    /// `face_labels[i]` is the face above `segments[i]` (optional).
    pub fn build(segments: &[Segment2D], _face_labels: &[Option<usize>]) -> Self {
        let mut xs: Vec<f64> = Vec::new();
        for s in segments {
            xs.push(s.left.x);
            xs.push(s.right.x);
        }
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        let n_slabs = if xs.len() >= 2 { xs.len() - 1 } else { 0 };
        let mut slab_segments = vec![Vec::new(); n_slabs];
        let mut slab_faces = vec![Vec::new(); n_slabs];
        for slab_i in 0..n_slabs {
            let x_mid = (xs[slab_i] + xs[slab_i + 1]) / 2.0;
            let mut active: Vec<(f64, usize)> = segments
                .iter()
                .enumerate()
                .filter_map(|(idx, s)| s.y_at(x_mid).map(|y| (y, idx)))
                .collect();
            active.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            slab_segments[slab_i] = active.iter().map(|(_, idx)| *idx).collect();
            slab_faces[slab_i] = active.iter().map(|_| None).collect();
        }
        Self {
            slab_xs: xs,
            slab_segments,
            segments: segments.to_vec(),
            slab_faces,
        }
    }
    /// Locate a query point: returns the index of the first segment strictly
    /// above the point in the slab, or `None` if the point is above all segments.
    pub fn locate(&self, p: Point2) -> Option<usize> {
        let slab_i = self
            .slab_xs
            .partition_point(|&x| x <= p.x)
            .saturating_sub(1);
        if slab_i >= self.slab_segments.len() {
            return None;
        }
        for &seg_idx in &self.slab_segments[slab_i] {
            if let Some(y) = self.segments[seg_idx].y_at(p.x)
                && y >= p.y
            {
                return Some(seg_idx);
            }
        }
        None
    }
}
/// A face of the 3D convex hull: three vertex indices with outward normal.
#[derive(Debug, Clone)]
pub struct ConvexFace3D {
    /// Indices of the three vertices (counter-clockwise when viewed from outside).
    pub verts: [usize; 3],
    /// Outward-facing unit normal.
    pub normal: Point3,
}
/// A vertex in a line arrangement: the intersection of two or more lines.
#[derive(Debug, Clone)]
pub struct ArrangementVertex {
    /// The intersection point.
    pub point: Point2,
    /// Indices of lines passing through this vertex.
    pub line_indices: Vec<usize>,
}
/// A 2D polygon obstacle (simple, non-self-intersecting).
#[derive(Debug, Clone)]
pub struct Obstacle2D {
    /// Vertices of the obstacle polygon (in order).
    pub vertices: Vec<Point2>,
}
impl Obstacle2D {
    /// Construct an obstacle from a vertex list.
    pub fn new(vertices: Vec<Point2>) -> Self {
        Self { vertices }
    }
    /// Collect all edges as (Point2, Point2) pairs.
    pub fn edges(&self) -> Vec<(Point2, Point2)> {
        let n = self.vertices.len();
        (0..n)
            .map(|i| (self.vertices[i], self.vertices[(i + 1) % n]))
            .collect()
    }
}
/// A point in 2D space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}
impl Point2 {
    /// Construct a new 2D point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    /// Subtract another point, yielding a vector.
    pub fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y)
    }
    /// Add another point / vector.
    pub fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y)
    }
    /// Scale by a scalar.
    pub fn scale(self, t: f64) -> Self {
        Self::new(self.x * t, self.y * t)
    }
    /// 2D cross product (scalar z-component of the 3D cross product).
    pub fn cross(self, other: Self) -> f64 {
        self.x * other.y - self.y * other.x
    }
    /// Dot product.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }
    /// Euclidean distance squared to another point.
    pub fn dist_sq(self, other: Self) -> f64 {
        let d = self.sub(other);
        d.dot(d)
    }
    /// Euclidean distance to another point.
    pub fn dist(self, other: Self) -> f64 {
        self.dist_sq(other).sqrt()
    }
    /// 2D cross product of vectors (p1-p0) and (p2-p0).
    pub fn cross2(p0: Self, p1: Self, p2: Self) -> f64 {
        p1.sub(p0).cross(p2.sub(p0))
    }
}
/// A manifold polygon mesh represented with the half-edge data structure.
///
/// Supports O(1) adjacency queries: vertex–face, face–face, and edge–edge.
#[derive(Debug, Clone, Default)]
pub struct HalfEdgeMesh {
    /// All half-edges.
    pub half_edges: Vec<HalfEdge>,
    /// All vertices.
    pub vertices: Vec<HEVertex>,
    /// All faces.
    pub faces: Vec<HEFace>,
}
impl HalfEdgeMesh {
    /// Create an empty half-edge mesh.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a vertex at the given position and return its [`VertexId`].
    pub fn add_vertex(&mut self, pos: Point3) -> VertexId {
        let id = self.vertices.len();
        self.vertices.push(HEVertex {
            position: pos,
            half_edge: None,
        });
        id
    }
    /// Add a triangular face from three vertex indices (counter-clockwise winding).
    ///
    /// Allocates three half-edges and one face record.  Twin links are resolved
    /// lazily by `build_twin_links`.
    pub fn add_triangle(&mut self, v0: VertexId, v1: VertexId, v2: VertexId) -> FaceId {
        let face_id = self.faces.len();
        let he0 = self.half_edges.len();
        let he1 = he0 + 1;
        let he2 = he0 + 2;
        self.half_edges.push(HalfEdge {
            origin: v0,
            twin: usize::MAX,
            next: he1,
            prev: he2,
            face: Some(face_id),
        });
        self.half_edges.push(HalfEdge {
            origin: v1,
            twin: usize::MAX,
            next: he2,
            prev: he0,
            face: Some(face_id),
        });
        self.half_edges.push(HalfEdge {
            origin: v2,
            twin: usize::MAX,
            next: he0,
            prev: he1,
            face: Some(face_id),
        });
        self.faces.push(HEFace { half_edge: he0 });
        if self.vertices[v0].half_edge.is_none() {
            self.vertices[v0].half_edge = Some(he0);
        }
        if self.vertices[v1].half_edge.is_none() {
            self.vertices[v1].half_edge = Some(he1);
        }
        if self.vertices[v2].half_edge.is_none() {
            self.vertices[v2].half_edge = Some(he2);
        }
        face_id
    }
    /// Walk around a face and collect the vertex indices.
    pub fn face_vertices(&self, fid: FaceId) -> Vec<VertexId> {
        let start = self.faces[fid].half_edge;
        let mut verts = Vec::new();
        let mut cur = start;
        loop {
            verts.push(self.half_edges[cur].origin);
            cur = self.half_edges[cur].next;
            if cur == start {
                break;
            }
        }
        verts
    }
    /// Collect all faces that share vertex `vid`.
    pub fn vertex_faces(&self, vid: VertexId) -> Vec<FaceId> {
        let start_he = match self.vertices[vid].half_edge {
            Some(h) => h,
            None => return vec![],
        };
        let mut faces = Vec::new();
        let mut cur = start_he;
        loop {
            if let Some(f) = self.half_edges[cur].face {
                faces.push(f);
            }
            let twin = self.half_edges[cur].twin;
            if twin == usize::MAX {
                break;
            }
            cur = self.half_edges[twin].next;
            if cur == start_he {
                break;
            }
        }
        faces
    }
    /// Resolve twin half-edge links by matching directed edges (v_a, v_b) with (v_b, v_a).
    pub fn build_twin_links(&mut self) {
        use std::collections::HashMap;
        let n = self.half_edges.len();
        let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
        for i in 0..n {
            let origin = self.half_edges[i].origin;
            let dest = self.half_edges[self.half_edges[i].next].origin;
            edge_map.insert((origin, dest), i);
        }
        for i in 0..n {
            let origin = self.half_edges[i].origin;
            let dest = self.half_edges[self.half_edges[i].next].origin;
            if let Some(&twin) = edge_map.get(&(dest, origin)) {
                self.half_edges[i].twin = twin;
            }
        }
    }
    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }
    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }
    /// Compute the face normal (assumes planar convex polygon).
    pub fn face_normal(&self, fid: FaceId) -> Point3 {
        let verts = self.face_vertices(fid);
        if verts.len() < 3 {
            return [0.0; 3];
        }
        let p0 = self.vertices[verts[0]].position;
        let p1 = self.vertices[verts[1]].position;
        let p2 = self.vertices[verts[2]].position;
        let e1 = sub3(p1, p0);
        let e2 = sub3(p2, p0);
        normalize3(cross3(e1, e2))
    }
    /// Compute the centroid of a face.
    pub fn face_centroid(&self, fid: FaceId) -> Point3 {
        let verts = self.face_vertices(fid);
        let n = verts.len() as f64;
        let mut acc = [0.0f64; 3];
        for vid in &verts {
            let p = self.vertices[*vid].position;
            acc[0] += p[0];
            acc[1] += p[1];
            acc[2] += p[2];
        }
        [acc[0] / n, acc[1] / n, acc[2] / n]
    }
}
/// Result of an incremental 3D convex hull computation.
#[derive(Debug, Clone)]
pub struct ConvexHull3D {
    /// Vertices on the hull (subset of input points).
    pub vertices: Vec<Point3>,
    /// Triangular faces.
    pub faces: Vec<ConvexFace3D>,
}
impl ConvexHull3D {
    /// Compute the volume of the convex hull using the divergence theorem.
    pub fn volume(&self) -> f64 {
        let mut vol = 0.0f64;
        for face in &self.faces {
            let a = self.vertices[face.verts[0]];
            let b = self.vertices[face.verts[1]];
            let c = self.vertices[face.verts[2]];
            vol += dot3(a, cross3(b, c));
        }
        (vol / 6.0).abs()
    }
    /// Compute the surface area of the convex hull.
    pub fn surface_area(&self) -> f64 {
        let mut area = 0.0f64;
        for face in &self.faces {
            let a = self.vertices[face.verts[0]];
            let b = self.vertices[face.verts[1]];
            let c = self.vertices[face.verts[2]];
            let ab = sub3(b, a);
            let ac = sub3(c, a);
            area += 0.5 * mag3(cross3(ab, ac));
        }
        area
    }
}
/// Result of the art gallery approximation.
#[derive(Debug, Clone)]
pub struct ArtGalleryResult {
    /// Indices of the chosen guard vertices.
    pub guards: Vec<usize>,
    /// Coverage: for each vertex, `covered[i] = true` if vertex `i` is visible
    /// from at least one guard.
    pub covered: Vec<bool>,
}
