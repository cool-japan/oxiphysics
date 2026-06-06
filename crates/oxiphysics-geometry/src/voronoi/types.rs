//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{bowyer_watson, is_delaunay};
use super::functions::{
    circumcircle, circumsphere_3d, clip_polygon_to_bounds, delaunay_to_voronoi,
    flip_edge_if_needed, lloyd_relaxation, nearest_site, polygon_area_pt2d, polygon_centroid_pt2d,
};

/// A Delaunay triangulation in the new OO-style API.
#[derive(Debug, Clone)]
pub struct DelaunayTriangulation {
    /// Input points.
    pub points: Vec<Point2D>,
    /// Triangles as triples of indices into `points`.
    pub triangles: Vec<[usize; 3]>,
}
impl DelaunayTriangulation {
    /// Build a Delaunay triangulation from a set of 2-D points using
    /// the Bowyer-Watson incremental algorithm.
    pub fn from_points(pts: &[[f64; 2]]) -> Self {
        let points: Vec<Point2D> = pts.iter().map(|&p| Point2D::from_array(p)).collect();
        let raw = bowyer_watson(pts);
        let triangles: Vec<[usize; 3]> = raw.iter().map(|t| [t.a, t.b, t.c]).collect();
        Self { points, triangles }
    }
    /// Compute the circumscribed circle of a triangle given by three vertex indices.
    ///
    /// Returns `([cx, cy], radius)` or `([0,0], 0)` if degenerate.
    pub fn circumscribed_circle(&self, tri: [usize; 3]) -> ([f64; 2], f64) {
        let pa = self.points[tri[0]].to_array();
        let pb = self.points[tri[1]].to_array();
        let pc = self.points[tri[2]].to_array();
        match circumcircle(pa, pb, pc) {
            Some((center, r2)) => (center, r2.sqrt()),
            None => ([0.0, 0.0], 0.0),
        }
    }
    /// Check whether the triangulation satisfies the Delaunay condition.
    pub fn is_delaunay(&self) -> bool {
        let pts: Vec<[f64; 2]> = self.points.iter().map(|p| p.to_array()).collect();
        let tris: Vec<DelaunayTriangle> = self
            .triangles
            .iter()
            .map(|&t| DelaunayTriangle::new(t[0], t[1], t[2]))
            .collect();
        is_delaunay(&pts, &tris)
    }
    /// Attempt to flip edges that violate the Delaunay condition.
    ///
    /// Performs a single pass over all triangle pairs.
    pub fn flip_edge(&mut self) {
        let pts: Vec<[f64; 2]> = self.points.iter().map(|p| p.to_array()).collect();
        let mut tris: Vec<DelaunayTriangle> = self
            .triangles
            .iter()
            .map(|&t| DelaunayTriangle::new(t[0], t[1], t[2]))
            .collect();
        let n = tris.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if j < tris.len() {
                    flip_edge_if_needed(&mut tris, &pts, i, j);
                }
            }
        }
        self.triangles = tris.iter().map(|t| [t.a, t.b, t.c]).collect();
    }
    /// Compute the dual Voronoi diagram of this triangulation.
    ///
    /// Uses a tight bounding box around the point set with 20% margin.
    pub fn dual_voronoi(&self) -> VoronoiDiagram {
        let pts: Vec<[f64; 2]> = self.points.iter().map(|p| p.to_array()).collect();
        if pts.is_empty() {
            return VoronoiDiagram {
                sites: vec![],
                cells: vec![],
                bbox: [0.0, 1.0, 0.0, 1.0],
            };
        }
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for &p in &pts {
            if p[0] < xmin {
                xmin = p[0];
            }
            if p[0] > xmax {
                xmax = p[0];
            }
            if p[1] < ymin {
                ymin = p[1];
            }
            if p[1] > ymax {
                ymax = p[1];
            }
        }
        let margin_x = (xmax - xmin).max(1.0) * 0.2;
        let margin_y = (ymax - ymin).max(1.0) * 0.2;
        let bbox = [
            xmin - margin_x,
            xmax + margin_x,
            ymin - margin_y,
            ymax + margin_y,
        ];
        VoronoiDiagram::from_sites(&pts, bbox)
    }
}
/// A complete Voronoi diagram with sites, cells, and bounding box.
///
/// The bounding box is `[xmin, xmax, ymin, ymax]`.
#[derive(Debug, Clone)]
pub struct VoronoiDiagram {
    /// Generating sites.
    pub sites: Vec<VoronoiSite>,
    /// Voronoi cells (one per site).
    pub cells: Vec<VoronoiCell>,
    /// Bounding box `[xmin, xmax, ymin, ymax]`.
    pub bbox: [f64; 4],
}
impl VoronoiDiagram {
    /// Build a Voronoi diagram from a set of points and a bounding box.
    ///
    /// Uses the Bowyer-Watson algorithm + dual graph construction, with
    /// Sutherland-Hodgman clipping to the bounding box.
    ///
    /// # Parameters
    /// * `pts` – array of `[x, y]` positions.
    /// * `bbox` – `[xmin, xmax, ymin, ymax]`.
    pub fn from_sites(pts: &[[f64; 2]], bbox: [f64; 4]) -> Self {
        let sites: Vec<VoronoiSite> = pts
            .iter()
            .enumerate()
            .map(|(i, &p)| VoronoiSite::new(Point2D::from_array(p), i))
            .collect();
        if pts.is_empty() {
            return Self {
                sites,
                cells: vec![],
                bbox,
            };
        }
        let triangles = bowyer_watson(pts);
        let legacy_cells = delaunay_to_voronoi(pts, &triangles);
        let cells: Vec<VoronoiCell> = legacy_cells
            .iter()
            .enumerate()
            .map(|(i, lc)| {
                let clipped = clip_polygon_to_bounds(&lc.vertices, bbox);
                let verts: Vec<Point2D> = clipped.iter().map(|&v| Point2D::from_array(v)).collect();
                VoronoiCell::new(i, verts)
            })
            .collect();
        Self { sites, cells, bbox }
    }
    /// Find the index of the nearest site to point `p`.
    pub fn nearest_site(&self, p: [f64; 2]) -> usize {
        if self.sites.is_empty() {
            return 0;
        }
        let pts: Vec<[f64; 2]> = self.sites.iter().map(|s| s.pos.to_array()).collect();
        nearest_site(p, &pts)
    }
    /// Return the area of cell `i` (pre-computed from the cell polygon).
    pub fn cell_area(&self, i: usize) -> f64 {
        if i < self.cells.len() {
            self.cells[i].area
        } else {
            0.0
        }
    }
    /// Compute the centroid of cell `i`.
    pub fn centroid(&self, i: usize) -> [f64; 2] {
        if i < self.cells.len() {
            let c = self.cells[i].centroid();
            [c.x, c.y]
        } else {
            [0.0, 0.0]
        }
    }
    /// Perform Lloyd relaxation: move each site to the centroid of its cell,
    /// recompute the diagram, and repeat for `n_iter` iterations.
    pub fn lloyd_relaxation(&self, n_iter: usize) -> Self {
        let pts_arr: Vec<[f64; 2]> = self.sites.iter().map(|s| s.pos.to_array()).collect();
        let relaxed = lloyd_relaxation(&pts_arr, self.bbox, n_iter);
        Self::from_sites(&relaxed, self.bbox)
    }
}
/// A 2-D point with named x/y coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}
impl Point2D {
    /// Create a new 2-D point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    /// Squared distance to another point.
    #[inline]
    pub fn dist2(&self, other: &Point2D) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2)
    }
    /// Euclidean distance to another point.
    #[inline]
    pub fn dist(&self, other: &Point2D) -> f64 {
        self.dist2(other).sqrt()
    }
    /// Convert to array form.
    #[inline]
    pub fn to_array(&self) -> [f64; 2] {
        [self.x, self.y]
    }
    /// Create from array.
    #[inline]
    pub fn from_array(a: [f64; 2]) -> Self {
        Self { x: a[0], y: a[1] }
    }
}
/// A triangle in the Delaunay triangulation, storing indices into the point array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DelaunayTriangle {
    /// Index of first vertex.
    pub a: usize,
    /// Index of second vertex.
    pub b: usize,
    /// Index of third vertex.
    pub c: usize,
}
impl DelaunayTriangle {
    /// Create a new Delaunay triangle from three vertex indices.
    pub fn new(a: usize, b: usize, c: usize) -> Self {
        Self { a, b, c }
    }
    /// Returns all three vertex indices as an array.
    pub fn indices(&self) -> [usize; 3] {
        [self.a, self.b, self.c]
    }
    /// Check if this triangle contains any of the super-triangle vertices (indices >= offset).
    pub fn contains_super_vertex(&self, offset: usize) -> bool {
        self.a >= offset || self.b >= offset || self.c >= offset
    }
}
/// A Voronoi cell in the new API.
///
/// Contains the site index, polygon vertices (ordered), and pre-computed area.
#[derive(Debug, Clone)]
pub struct VoronoiCell {
    /// Index of the generating site in the parent diagram.
    pub site_index: usize,
    /// Ordered polygon vertices of the cell boundary.
    pub vertices: Vec<Point2D>,
    /// Pre-computed area of the cell.
    pub area: f64,
}
impl VoronoiCell {
    /// Create a new Voronoi cell and compute its area.
    pub fn new(site_index: usize, vertices: Vec<Point2D>) -> Self {
        let area = polygon_area_pt2d(&vertices);
        Self {
            site_index,
            vertices,
            area,
        }
    }
    /// Compute the centroid of this cell's polygon.
    pub fn centroid(&self) -> Point2D {
        polygon_centroid_pt2d(&self.vertices)
    }
    /// Perimeter of the cell polygon.
    pub fn perimeter(&self) -> f64 {
        let n = self.vertices.len();
        if n < 2 {
            return 0.0;
        }
        (0..n)
            .map(|i| self.vertices[i].dist(&self.vertices[(i + 1) % n]))
            .sum()
    }
}
/// A legacy Voronoi cell: the region closest to `site`.
///
/// `vertices` is the (possibly empty for unbounded) polygon of circumcenters.
#[derive(Debug, Clone)]
pub struct LegacyVoronoiCell {
    /// The generating site (input point).
    pub site: [f64; 2],
    /// Ordered polygon vertices (circumcenters of adjacent triangles).
    pub vertices: Vec<[f64; 2]>,
}
impl LegacyVoronoiCell {
    /// Create a new legacy Voronoi cell.
    pub fn new(site: [f64; 2], vertices: Vec<[f64; 2]>) -> Self {
        Self { site, vertices }
    }
}
/// A tetrahedron in a 3-D Delaunay tetrahedralisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DelaunayTetrahedron {
    /// Vertex indices (into the point array).
    pub v: [usize; 4],
}
impl DelaunayTetrahedron {
    /// Create from four vertex indices.
    pub fn new(a: usize, b: usize, c: usize, d: usize) -> Self {
        Self { v: [a, b, c, d] }
    }
    /// Return `true` if any vertex index is ≥ `offset` (super-tetrahedron detection).
    pub fn contains_super(&self, offset: usize) -> bool {
        self.v.iter().any(|&i| i >= offset)
    }
    /// Compute the circumsphere of this tetrahedron.
    pub fn circumsphere(&self, pts: &[[f64; 3]]) -> Option<([f64; 3], f64)> {
        circumsphere_3d(
            pts[self.v[0]],
            pts[self.v[1]],
            pts[self.v[2]],
            pts[self.v[3]],
        )
    }
}
/// A weighted site for the power diagram (Laguerre diagram).
///
/// The power distance from a query point `q` to a weighted site `(s, w)` is
/// `|q - s|² - w`.
#[derive(Debug, Clone, Copy)]
pub struct WeightedSite {
    /// 2-D position.
    pub position: [f64; 2],
    /// Weight (radius² for additive-weighted Voronoi / power diagram).
    pub weight: f64,
}
impl WeightedSite {
    /// Create a new weighted site.
    pub fn new(position: [f64; 2], weight: f64) -> Self {
        Self { position, weight }
    }
    /// Power distance from this site to query point `q`.
    pub fn power_distance(&self, q: [f64; 2]) -> f64 {
        let dx = q[0] - self.position[0];
        let dy = q[1] - self.position[1];
        dx * dx + dy * dy - self.weight
    }
}
/// A Voronoi generating site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoronoiSite {
    /// Position of the site.
    pub pos: Point2D,
    /// Index of this site in the diagram's site list.
    pub index: usize,
}
impl VoronoiSite {
    /// Create a new Voronoi site.
    pub fn new(pos: Point2D, index: usize) -> Self {
        Self { pos, index }
    }
}
