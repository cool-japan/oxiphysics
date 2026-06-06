//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{
    add, cross, dot, len, normalize, ray_aabb, ray_sphere, ray_triangle, reflect_ray, scale, sub,
};

/// The result of a successful ray-primitive intersection test.
#[derive(Debug, Clone, Copy)]
pub struct RayHit {
    /// Parameter `t` along the ray at the hit point.
    pub t: f64,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Surface normal at the hit point (unit length, outward).
    pub normal: [f64; 3],
    /// UV texture coordinates at the hit point `(u, v)`.
    pub uv: [f64; 2],
    /// Index of the hit primitive (face / mesh index).
    pub prim_id: usize,
}
impl RayHit {
    /// Construct a new ray hit record.
    pub fn new(t: f64, point: [f64; 3], normal: [f64; 3], uv: [f64; 2], prim_id: usize) -> Self {
        Self {
            t,
            point,
            normal,
            uv,
            prim_id,
        }
    }
    /// Return `true` if this hit is closer than `other`.
    pub fn closer_than(&self, other: &RayHit) -> bool {
        self.t < other.t
    }
}
/// Result of a batch ray-casting pass.
#[derive(Debug, Clone)]
pub struct BatchRayResult {
    /// Per-ray closest hit, or `None` if the ray missed.
    pub hits: Vec<Option<RayHit>>,
}
/// A node in a recursive ray tree.
///
/// Each ray may spawn up to two children: a reflected ray and/or a
/// refracted ray.
#[derive(Debug, Clone)]
pub struct RayTreeNode {
    /// The ray for this node.
    pub ray: Ray,
    /// Closest hit, if any.
    pub hit: Option<RayHit>,
    /// Attenuation factor for this ray's contribution.
    pub weight: f64,
    /// Reflected child (index into the tree node array).
    pub reflected_child: Option<usize>,
    /// Refracted child (index into the tree node array).
    pub refracted_child: Option<usize>,
}
impl RayTreeNode {
    /// Create a new ray tree node.
    pub fn new(ray: Ray, weight: f64) -> Self {
        Self {
            ray,
            hit: None,
            weight,
            reflected_child: None,
            refracted_child: None,
        }
    }
}
/// Capsule: a cylinder with hemispherical caps.
#[derive(Debug, Clone, Copy)]
pub struct Capsule {
    /// Centre of the first (bottom) hemisphere.
    pub a: [f64; 3],
    /// Centre of the second (top) hemisphere.
    pub b: [f64; 3],
    /// Capsule radius.
    pub radius: f64,
}
impl Capsule {
    /// Create a new capsule.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }
    /// Axis vector (unnormalised) from `a` to `b`.
    pub fn axis(&self) -> [f64; 3] {
        sub(self.b, self.a)
    }
    /// Length of the capsule (centre-to-centre).
    pub fn length(&self) -> f64 {
        len(self.axis())
    }
}
/// Oriented Bounding Box.
///
/// Stored as a centre, a rotation matrix (columns are OBB axes), and
/// half-extents along each local axis.
#[derive(Debug, Clone, Copy)]
pub struct Obb {
    /// Centre of the OBB.
    pub centre: [f64; 3],
    /// Rotation matrix – rows are the local X, Y, Z axes in world space.
    pub axes: [[f64; 3]; 3],
    /// Half-extents along each local axis.
    pub half_extents: [f64; 3],
}
impl Obb {
    /// Create an OBB from centre, axes and half-extents.
    pub fn new(centre: [f64; 3], axes: [[f64; 3]; 3], half_extents: [f64; 3]) -> Self {
        Self {
            centre,
            axes,
            half_extents,
        }
    }
    /// Create an axis-aligned OBB (convenience wrapper).
    pub fn from_aabb(aabb: &Aabb) -> Self {
        Self::new(
            aabb.centre(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            aabb.half_extents(),
        )
    }
}
/// Sphere primitive defined by a centre and radius.
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    /// Centre of the sphere.
    pub centre: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}
impl Sphere {
    /// Create a new sphere.
    pub fn new(centre: [f64; 3], radius: f64) -> Self {
        Self { centre, radius }
    }
    /// Compute the outward unit normal at a surface point.
    pub fn normal_at(&self, p: [f64; 3]) -> [f64; 3] {
        normalize(sub(p, self.centre))
    }
    /// Compute UV spherical coordinates `(u, v)` at a surface point.
    ///
    /// `u ∈ [0, 1]` maps azimuth, `v ∈ [0, 1]` maps elevation.
    pub fn uv_at(&self, p: [f64; 3]) -> [f64; 2] {
        let n = self.normal_at(p);
        let u = 1.0 - (n[0].atan2(n[2]) + PI) / (2.0 * PI);
        let v = (n[1].asin() + PI / 2.0) / PI;
        [u, v]
    }
}
/// A simple scene containing spheres and triangle meshes for ray casting.
pub struct Scene {
    /// Spheres in the scene.
    pub spheres: Vec<Sphere>,
    /// Triangle vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle index triples.
    pub triangles: Vec<(usize, usize, usize)>,
}
impl Scene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self {
            spheres: Vec::new(),
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }
    /// Add a sphere to the scene.
    pub fn add_sphere(&mut self, sphere: Sphere) {
        self.spheres.push(sphere);
    }
    /// Add a triangle by vertex positions.
    pub fn add_triangle(&mut self, v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) {
        let base = self.vertices.len();
        self.vertices.push(v0);
        self.vertices.push(v1);
        self.vertices.push(v2);
        self.triangles.push((base, base + 1, base + 2));
    }
    /// Cast a ray against all primitives in the scene and return the closest hit.
    pub fn cast(&self, ray: &Ray) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        for sphere in &self.spheres {
            if let Some(hit) = ray_sphere(ray, sphere) {
                best = Some(match best {
                    None => hit,
                    Some(prev) => {
                        if hit.t < prev.t {
                            hit
                        } else {
                            prev
                        }
                    }
                });
            }
        }
        for &(i, j, k) in &self.triangles {
            if let Some(hit) =
                ray_triangle(ray, self.vertices[i], self.vertices[j], self.vertices[k], i)
            {
                best = Some(match best {
                    None => hit,
                    Some(prev) => {
                        if hit.t < prev.t {
                            hit
                        } else {
                            prev
                        }
                    }
                });
            }
        }
        best
    }
}
/// A flat BVH built over a list of triangle primitives.
pub struct Bvh {
    /// Flat node array; node 0 is the root.
    pub nodes: Vec<BvhNode>,
    /// Primitives sorted according to the BVH build order.
    pub prim_indices: Vec<usize>,
    /// Vertex array shared with the scene.
    pub vertices: Vec<[f64; 3]>,
    /// Index triples (one per triangle).
    pub indices: Vec<(usize, usize, usize)>,
}
impl Bvh {
    /// Build a BVH over a triangle mesh using a simple median-split heuristic.
    pub fn build(vertices: Vec<[f64; 3]>, indices: Vec<(usize, usize, usize)>) -> Self {
        let n = indices.len();
        let mut prim_indices: Vec<usize> = (0..n).collect();
        let mut nodes = Vec::new();
        Self::build_recursive(&vertices, &indices, &mut prim_indices, 0, n, &mut nodes);
        Self {
            nodes,
            prim_indices,
            vertices,
            indices,
        }
    }
    fn triangle_centroid(vertices: &[[f64; 3]], tri: (usize, usize, usize)) -> [f64; 3] {
        let (i, j, k) = tri;
        let v0 = vertices[i];
        let v1 = vertices[j];
        let v2 = vertices[k];
        scale(add(add(v0, v1), v2), 1.0 / 3.0)
    }
    fn triangle_aabb(vertices: &[[f64; 3]], tri: (usize, usize, usize)) -> Aabb {
        let (i, j, k) = tri;
        let v0 = vertices[i];
        let v1 = vertices[j];
        let v2 = vertices[k];
        let min = [
            v0[0].min(v1[0]).min(v2[0]),
            v0[1].min(v1[1]).min(v2[1]),
            v0[2].min(v1[2]).min(v2[2]),
        ];
        let max = [
            v0[0].max(v1[0]).max(v2[0]),
            v0[1].max(v1[1]).max(v2[1]),
            v0[2].max(v1[2]).max(v2[2]),
        ];
        Aabb::new(min, max)
    }
    fn build_recursive(
        vertices: &[[f64; 3]],
        indices: &[(usize, usize, usize)],
        prim_indices: &mut Vec<usize>,
        start: usize,
        end: usize,
        nodes: &mut Vec<BvhNode>,
    ) -> usize {
        let node_idx = nodes.len();
        nodes.push(BvhNode::leaf(Aabb::new([0.0; 3], [0.0; 3]), 0, 0));
        let mut aabb = Aabb::new([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for &pi in &prim_indices[start..end] {
            let tri_bb = Self::triangle_aabb(vertices, indices[pi]);
            aabb = Aabb::merge(&aabb, &tri_bb);
        }
        let n_prims = end - start;
        if n_prims <= 2 {
            nodes[node_idx] = BvhNode::leaf(aabb, start, n_prims);
            return node_idx;
        }
        let d = sub(aabb.max, aabb.min);
        let axis = if d[0] >= d[1] && d[0] >= d[2] {
            0
        } else if d[1] >= d[2] {
            1
        } else {
            2
        };
        prim_indices[start..end].sort_by(|&a, &b| {
            let ca = Self::triangle_centroid(vertices, indices[a])[axis];
            let cb = Self::triangle_centroid(vertices, indices[b])[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = (start + end) / 2;
        let left = Self::build_recursive(vertices, indices, prim_indices, start, mid, nodes);
        let right = Self::build_recursive(vertices, indices, prim_indices, mid, end, nodes);
        nodes[node_idx] = BvhNode::internal(aabb, left, right);
        node_idx
    }
    /// Traverse the BVH and return the closest triangle hit.
    pub fn intersect(&self, ray: &Ray) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        let mut stack = vec![0usize];
        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            if ray_aabb(ray, &node.aabb).is_none() {
                continue;
            }
            if node.is_leaf() {
                for i in node.prim_start..(node.prim_start + node.prim_count) {
                    let pi = self.prim_indices[i];
                    let (vi, vj, vk) = self.indices[pi];
                    if vi >= self.vertices.len()
                        || vj >= self.vertices.len()
                        || vk >= self.vertices.len()
                    {
                        continue;
                    }
                    if let Some(hit) = ray_triangle(
                        ray,
                        self.vertices[vi],
                        self.vertices[vj],
                        self.vertices[vk],
                        pi,
                    ) {
                        best = Some(match best {
                            None => hit,
                            Some(prev) => {
                                if hit.t < prev.t {
                                    hit
                                } else {
                                    prev
                                }
                            }
                        });
                    }
                }
            } else {
                stack.push(node.left);
                stack.push(node.right);
            }
        }
        best
    }
}
/// A convex mesh represented as a list of vertices.
///
/// The support function finds the vertex that is most extreme along a
/// given direction.
#[derive(Debug, Clone)]
pub struct ConvexMesh {
    /// Vertex list in world space.
    pub vertices: Vec<[f64; 3]>,
}
impl ConvexMesh {
    /// Create a convex mesh from a vertex list.
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        Self { vertices }
    }
    /// Support function: returns the vertex most extreme in direction `dir`.
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        self.vertices
            .iter()
            .copied()
            .fold(self.vertices[0], |best, v| {
                if dot(v, dir) > dot(best, dir) {
                    v
                } else {
                    best
                }
            })
    }
    /// Compute the AABB of the mesh.
    pub fn aabb(&self) -> Aabb {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for &v in &self.vertices {
            for i in 0..3 {
                min[i] = min[i].min(v[i]);
                max[i] = max[i].max(v[i]);
            }
        }
        Aabb::new(min, max)
    }
}
/// Axis-Aligned Bounding Box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl Aabb {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Return the AABB centre.
    pub fn centre(&self) -> [f64; 3] {
        scale(add(self.min, self.max), 0.5)
    }
    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        scale(sub(self.max, self.min), 0.5)
    }
    /// Surface area of the AABB (used in SAH BVH cost).
    pub fn surface_area(&self) -> f64 {
        let d = sub(self.max, self.min);
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }
    /// Merge two AABBs into the smallest containing AABB.
    pub fn merge(a: &Aabb, b: &Aabb) -> Aabb {
        Aabb::new(
            [
                a.min[0].min(b.min[0]),
                a.min[1].min(b.min[1]),
                a.min[2].min(b.min[2]),
            ],
            [
                a.max[0].max(b.max[0]),
                a.max[1].max(b.max[1]),
                a.max[2].max(b.max[2]),
            ],
        )
    }
    /// Expand the AABB to include point `p`.
    pub fn expand(&self, p: [f64; 3]) -> Aabb {
        Aabb::new(
            [
                self.min[0].min(p[0]),
                self.min[1].min(p[1]),
                self.min[2].min(p[2]),
            ],
            [
                self.max[0].max(p[0]),
                self.max[1].max(p[1]),
                self.max[2].max(p[2]),
            ],
        )
    }
}
/// A flat recursive ray tree for a single primary ray.
///
/// Supports up to `max_depth` bounces with reflection and refraction.
pub struct RayTree {
    /// Flat array of ray tree nodes.
    pub nodes: Vec<RayTreeNode>,
    /// Maximum allowed bounce depth.
    pub max_depth: usize,
}
impl RayTree {
    /// Create an empty ray tree.
    pub fn new(max_depth: usize) -> Self {
        Self {
            nodes: Vec::new(),
            max_depth,
        }
    }
    /// Build the ray tree by tracing against a list of spheres.
    ///
    /// Recursively spawns reflection rays when hitting a sphere.
    pub fn trace(&mut self, primary: Ray, spheres: &[Sphere]) {
        self.nodes.clear();
        self.nodes.push(RayTreeNode::new(primary, 1.0));
        let mut i = 0;
        while i < self.nodes.len() {
            let depth = self.compute_depth(i);
            if depth >= self.max_depth {
                i += 1;
                continue;
            }
            let ray = self.nodes[i].ray;
            let weight = self.nodes[i].weight;
            let mut best: Option<RayHit> = None;
            for sphere in spheres {
                if let Some(hit) = ray_sphere(&ray, sphere) {
                    best = Some(match best {
                        None => hit,
                        Some(prev) => {
                            if hit.t < prev.t {
                                hit
                            } else {
                                prev
                            }
                        }
                    });
                }
            }
            self.nodes[i].hit = best;
            if let Some(hit) = best
                && weight > 0.01
            {
                let r_dir = reflect_ray(ray.direction, hit.normal);
                let r_ray = ray.spawn_offset(hit.point, r_dir);
                let child_idx = self.nodes.len();
                self.nodes.push(RayTreeNode::new(r_ray, weight * 0.8));
                self.nodes[i].reflected_child = Some(child_idx);
            }
            i += 1;
        }
    }
    /// Compute approximate depth of node `i` (linear search for parent).
    fn compute_depth(&self, target: usize) -> usize {
        let mut depth = 0usize;
        let mut current = target;
        'outer: loop {
            for (parent, node) in self.nodes.iter().enumerate() {
                if node.reflected_child == Some(current) || node.refracted_child == Some(current) {
                    depth += 1;
                    current = parent;
                    continue 'outer;
                }
            }
            break;
        }
        depth
    }
}
/// A ray with associated differentials, used for texture anti-aliasing
/// (ray differentials method, Igehy 1999).
#[derive(Debug, Clone, Copy)]
pub struct RayDifferential {
    /// Primary ray.
    pub ray: Ray,
    /// Offset ray slightly in the screen-X direction.
    pub rx: Ray,
    /// Offset ray slightly in the screen-Y direction.
    pub ry: Ray,
    /// Pixel footprint scale factor.
    pub pixel_spread: f64,
}
impl RayDifferential {
    /// Create a ray differential from a primary ray and two offset rays.
    pub fn new(ray: Ray, rx: Ray, ry: Ray) -> Self {
        Self {
            ray,
            rx,
            ry,
            pixel_spread: 1.0,
        }
    }
    /// Compute the area of the projected pixel footprint at parameter `t`.
    ///
    /// Returns an estimate of the screen-space pixel footprint on the surface.
    pub fn footprint_at(&self, t: f64) -> f64 {
        let p = self.ray.at(t);
        let px = self.rx.at(t);
        let py = self.ry.at(t);
        let dx = sub(px, p);
        let dy = sub(py, p);
        let c = cross(dx, dy);
        len(c)
    }
    /// Transfer differentials to a surface hit point.
    ///
    /// Given the surface normal at the hit, propagates the differentials
    /// using the plane transfer equations.
    pub fn transfer_to_surface(&self, t: f64, normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let n = normalize(normal);
        let p = self.ray.at(t);
        let px_world = self.rx.at(t);
        let py_world = self.ry.at(t);
        let transfer_point = |q: [f64; 3]| -> [f64; 3] {
            let offset = dot(sub(q, p), n);
            sub(q, scale(n, offset))
        };
        (transfer_point(px_world), transfer_point(py_world))
    }
}
/// A node in a flat BVH.
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// AABB of this node.
    pub aabb: Aabb,
    /// Index of the left child, or `usize::MAX` for a leaf.
    pub left: usize,
    /// Index of the right child, or `usize::MAX` for a leaf.
    pub right: usize,
    /// Leaf: starting index in the primitive array.
    pub prim_start: usize,
    /// Leaf: number of primitives.
    pub prim_count: usize,
}
impl BvhNode {
    /// Create an internal node.
    pub fn internal(aabb: Aabb, left: usize, right: usize) -> Self {
        Self {
            aabb,
            left,
            right,
            prim_start: 0,
            prim_count: 0,
        }
    }
    /// Create a leaf node.
    pub fn leaf(aabb: Aabb, prim_start: usize, prim_count: usize) -> Self {
        Self {
            aabb,
            left: usize::MAX,
            right: usize::MAX,
            prim_start,
            prim_count,
        }
    }
    /// Return `true` if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left == usize::MAX
    }
}
/// A height field stored as a grid of height values.
///
/// The grid covers `[0, width] × [0, height]` in the XZ plane.
pub struct Heightfield {
    /// Height values at grid cells (row-major: `heights[row * cols + col]`).
    pub heights: Vec<f64>,
    /// Number of columns (X direction).
    pub cols: usize,
    /// Number of rows (Z direction).
    pub rows: usize,
    /// Grid spacing (same in both X and Z).
    pub cell_size: f64,
}
impl Heightfield {
    /// Create a new heightfield.
    pub fn new(heights: Vec<f64>, cols: usize, rows: usize, cell_size: f64) -> Self {
        Self {
            heights,
            cols,
            rows,
            cell_size,
        }
    }
    /// Bilinear interpolation of height at `(x, z)`.
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let xi = (x / self.cell_size).max(0.0);
        let zi = (z / self.cell_size).max(0.0);
        let ix = (xi as usize).min(self.cols - 2);
        let iz = (zi as usize).min(self.rows - 2);
        let fx = xi - ix as f64;
        let fz = zi - iz as f64;
        let h00 = self.heights[iz * self.cols + ix];
        let h10 = self.heights[iz * self.cols + ix + 1];
        let h01 = self.heights[(iz + 1) * self.cols + ix];
        let h11 = self.heights[(iz + 1) * self.cols + ix + 1];
        h00 * (1.0 - fx) * (1.0 - fz)
            + h10 * fx * (1.0 - fz)
            + h01 * (1.0 - fx) * fz
            + h11 * fx * fz
    }
    /// Total width of the heightfield (X direction).
    pub fn width(&self) -> f64 {
        (self.cols - 1) as f64 * self.cell_size
    }
    /// Total depth of the heightfield (Z direction).
    pub fn depth(&self) -> f64 {
        (self.rows - 1) as f64 * self.cell_size
    }
}
/// An infinite ray parameterised as `origin + t * direction`.
///
/// The direction should be a unit vector for correct distance computations;
/// use [`Ray::new_normalised`] to construct one automatically.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Ray direction (unit length).
    pub direction: [f64; 3],
    /// Minimum valid parameter value (clipping near plane).
    pub t_min: f64,
    /// Maximum valid parameter value (clipping far plane).
    pub t_max: f64,
}
impl Ray {
    /// Create a ray with pre-normalised direction and default `[0, ∞)` range.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction,
            t_min: 0.0,
            t_max: f64::INFINITY,
        }
    }
    /// Create a ray and automatically normalise the direction vector.
    pub fn new_normalised(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self::new(origin, normalize(direction))
    }
    /// Evaluate the ray at parameter `t`: `P = origin + t * direction`.
    pub fn at(&self, t: f64) -> [f64; 3] {
        add(self.origin, scale(self.direction, t))
    }
    /// Return `true` if `t` lies within the valid range `[t_min, t_max]`.
    pub fn valid_t(&self, t: f64) -> bool {
        t >= self.t_min && t <= self.t_max
    }
    /// Spawn a secondary ray offset slightly to avoid self-intersection.
    pub fn spawn_offset(&self, origin: [f64; 3], direction: [f64; 3]) -> Ray {
        let eps = 1e-6;
        Ray::new(add(origin, scale(direction, eps)), normalize(direction))
    }
}
