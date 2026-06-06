//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// A material definition for ray tracing.
#[derive(Debug, Clone)]
pub struct Material {
    /// Material type.
    pub mat_type: MaterialType,
    /// Albedo / base color.
    pub albedo: [f64; 3],
    /// Roughness (0=smooth, 1=rough) for metal/PBR.
    pub roughness: f64,
    /// Metallic factor (0=dielectric, 1=metal) for PBR.
    pub metallic: f64,
    /// Index of refraction for dielectric.
    pub ior: f64,
    /// Emission color and strength.
    pub emission: [f64; 3],
    /// Specular exponent for Phong.
    pub shininess: f64,
    /// Ambient occlusion factor.
    pub ao: f64,
}
impl Material {
    /// Create a Lambertian diffuse material.
    pub fn diffuse(albedo: [f64; 3]) -> Self {
        Self {
            mat_type: MaterialType::Diffuse,
            albedo,
            roughness: 1.0,
            metallic: 0.0,
            ior: 1.0,
            emission: [0.0; 3],
            shininess: 0.0,
            ao: 1.0,
        }
    }
    /// Create a metal material.
    pub fn metal(albedo: [f64; 3], roughness: f64) -> Self {
        Self {
            mat_type: MaterialType::Metal,
            albedo,
            roughness: roughness.clamp(0.0, 1.0),
            metallic: 1.0,
            ior: 1.0,
            emission: [0.0; 3],
            shininess: 0.0,
            ao: 1.0,
        }
    }
    /// Create a glass (dielectric) material.
    pub fn glass(ior: f64) -> Self {
        Self {
            mat_type: MaterialType::Dielectric,
            albedo: [1.0; 3],
            roughness: 0.0,
            metallic: 0.0,
            ior,
            emission: [0.0; 3],
            shininess: 0.0,
            ao: 1.0,
        }
    }
    /// Create an emissive material.
    pub fn emissive(color: [f64; 3], strength: f64) -> Self {
        Self {
            mat_type: MaterialType::Emissive,
            albedo: color,
            roughness: 1.0,
            metallic: 0.0,
            ior: 1.0,
            emission: scale3(color, strength),
            shininess: 0.0,
            ao: 1.0,
        }
    }
    /// Create a PBR material.
    pub fn pbr(albedo: [f64; 3], metallic: f64, roughness: f64, ao: f64) -> Self {
        Self {
            mat_type: MaterialType::Pbr,
            albedo,
            roughness: roughness.clamp(0.0, 1.0),
            metallic: metallic.clamp(0.0, 1.0),
            ior: 1.5,
            emission: [0.0; 3],
            shininess: 0.0,
            ao,
        }
    }
}
/// Path tracer state for a single sample path.
#[derive(Debug, Clone)]
pub struct PathState {
    /// Current ray.
    pub ray: Ray,
    /// Accumulated throughput (product of BRDFs).
    pub throughput: [f64; 3],
    /// Accumulated radiance.
    pub radiance: [f64; 3],
    /// Current bounce depth.
    pub depth: u32,
    /// Maximum bounce depth.
    pub max_depth: u32,
}
impl PathState {
    /// Create a new path state.
    pub fn new(ray: Ray, max_depth: u32) -> Self {
        Self {
            ray,
            throughput: [1.0; 3],
            radiance: [0.0; 3],
            depth: 0,
            max_depth,
        }
    }
    /// Returns true if path tracing should continue.
    pub fn should_continue(&self) -> bool {
        self.depth < self.max_depth
            && (self.throughput[0] + self.throughput[1] + self.throughput[2]) > 1e-6
    }
    /// Apply Russian roulette termination. Returns false if path is terminated.
    pub fn russian_roulette(&mut self, survival_prob: f64) -> bool {
        if survival_prob >= 1.0 {
            return true;
        }
        let luminance =
            0.2126 * self.throughput[0] + 0.7152 * self.throughput[1] + 0.0722 * self.throughput[2];
        if luminance < survival_prob {
            return false;
        }
        self.throughput = scale3(self.throughput, 1.0 / survival_prob);
        true
    }
}
/// Result of a ray-scene intersection.
#[derive(Debug, Clone, Copy)]
pub struct HitRecord {
    /// Hit distance along the ray.
    pub t: f64,
    /// World-space hit position.
    pub position: [f64; 3],
    /// Outward-facing surface normal (normalized).
    pub normal: [f64; 3],
    /// UV texture coordinates.
    pub uv: [f64; 2],
    /// Index of the intersected triangle/primitive.
    pub prim_id: u32,
    /// Whether the ray hit the front face.
    pub front_face: bool,
    /// Index of the intersected material.
    pub material_id: u32,
}
impl HitRecord {
    /// Create a new hit record and set the face normal based on ray direction.
    pub fn new(
        t: f64,
        position: [f64; 3],
        outward_normal: [f64; 3],
        ray_dir: [f64; 3],
        uv: [f64; 2],
        prim_id: u32,
        material_id: u32,
    ) -> Self {
        let front_face = dot3(ray_dir, outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            scale3(outward_normal, -1.0)
        };
        Self {
            t,
            position,
            normal,
            uv,
            prim_id,
            front_face,
            material_id,
        }
    }
}
/// A node in a Bounding Volume Hierarchy.
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// Bounding box of this node.
    pub bounds: Aabb,
    /// For leaf: index of first primitive. For internal: left child index.
    pub left_or_first: u32,
    /// For leaf: primitive count (>0). For internal: 0.
    pub prim_count: u32,
}
impl BvhNode {
    /// Is this node a leaf?
    pub fn is_leaf(&self) -> bool {
        self.prim_count > 0
    }
}
/// A simple ray-traced scene.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    /// Scene triangles.
    pub triangles: Vec<Triangle>,
    /// Scene materials.
    pub materials: Vec<Material>,
    /// Point lights.
    pub lights: Vec<PointLight>,
    /// Area lights.
    pub area_lights: Vec<AreaLight>,
    /// Prebuilt BVH (None until built).
    pub bvh: Option<Bvh>,
}
impl Scene {
    /// Create a new empty scene.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a material and return its index.
    pub fn add_material(&mut self, mat: Material) -> u32 {
        let idx = self.materials.len() as u32;
        self.materials.push(mat);
        idx
    }
    /// Add a triangle.
    pub fn add_triangle(&mut self, tri: Triangle) {
        self.triangles.push(tri);
    }
    /// Add a point light.
    pub fn add_light(&mut self, light: PointLight) {
        self.lights.push(light);
    }
    /// Build the BVH over the current triangles.
    pub fn build_bvh(&mut self) {
        self.bvh = Some(Bvh::build(&self.triangles));
    }
    /// Add a unit box (12 triangles) centered at `center` with half-size `hs`.
    pub fn add_box(&mut self, center: [f64; 3], hs: [f64; 3], material_id: u32) {
        let [cx, cy, cz] = center;
        let [hx, hy, hz] = hs;
        let v = [
            [cx - hx, cy - hy, cz - hz],
            [cx + hx, cy - hy, cz - hz],
            [cx + hx, cy + hy, cz - hz],
            [cx - hx, cy + hy, cz - hz],
            [cx - hx, cy - hy, cz + hz],
            [cx + hx, cy - hy, cz + hz],
            [cx + hx, cy + hy, cz + hz],
            [cx - hx, cy + hy, cz + hz],
        ];
        let normals = [
            [0.0f64, 0.0, -1.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let faces: [[usize; 4]; 6] = [
            [0, 1, 2, 3],
            [5, 4, 7, 6],
            [4, 0, 3, 7],
            [1, 5, 6, 2],
            [4, 5, 1, 0],
            [3, 2, 6, 7],
        ];
        let uv_quad = [[0.0f64; 2], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        for (fi, face) in faces.iter().enumerate() {
            let n = normals[fi];
            let t0 = Triangle::new(
                v[face[0]],
                v[face[1]],
                v[face[2]],
                n,
                n,
                n,
                uv_quad[0],
                uv_quad[1],
                uv_quad[2],
                material_id,
            );
            let t1 = Triangle::new(
                v[face[0]],
                v[face[2]],
                v[face[3]],
                n,
                n,
                n,
                uv_quad[0],
                uv_quad[2],
                uv_quad[3],
                material_id,
            );
            self.add_triangle(t0);
            self.add_triangle(t1);
        }
    }
    /// Add a planar quad (2 triangles).
    pub fn add_quad(
        &mut self,
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        v3: [f64; 3],
        material_id: u32,
    ) {
        let n = normalize3(cross3(sub3(v1, v0), sub3(v2, v0)));
        let t0 = Triangle::new(
            v0,
            v1,
            v2,
            n,
            n,
            n,
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            material_id,
        );
        let t1 = Triangle::new(
            v0,
            v2,
            v3,
            n,
            n,
            n,
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            material_id,
        );
        self.add_triangle(t0);
        self.add_triangle(t1);
    }
    /// Intersect the scene with a ray using the prebuilt BVH.
    pub fn intersect<'a>(&'a self, ray: &Ray) -> Option<(HitRecord, &'a Triangle)> {
        match &self.bvh {
            Some(bvh) => bvh.intersect(ray, &self.triangles),
            None => {
                let mut best_t = ray.t_max;
                let mut best: Option<(HitRecord, usize)> = None;
                for (i, tri) in self.triangles.iter().enumerate() {
                    if let Some((t, u, v)) = tri.intersect_full(ray)
                        && t < best_t
                    {
                        best_t = t;
                        let pos = ray.at(t);
                        let norm = tri.interpolate_normal(u, v);
                        let uv = tri.interpolate_uv(u, v);
                        let hit = HitRecord::new(
                            t,
                            pos,
                            norm,
                            ray.direction,
                            uv,
                            i as u32,
                            tri.material_id,
                        );
                        best = Some((hit, i));
                    }
                }
                best.map(|(hit, i)| (hit, &self.triangles[i]))
            }
        }
    }
}
/// A triangle primitive for ray tracing.
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// Vertex positions.
    pub v: [[f64; 3]; 3],
    /// Per-vertex normals.
    pub n: [[f64; 3]; 3],
    /// Per-vertex UV coordinates.
    pub uv: [[f64; 2]; 3],
    /// Material index.
    pub material_id: u32,
}
impl Triangle {
    /// Create a new triangle.
    pub fn new(
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        n0: [f64; 3],
        n1: [f64; 3],
        n2: [f64; 3],
        uv0: [f64; 2],
        uv1: [f64; 2],
        uv2: [f64; 2],
        material_id: u32,
    ) -> Self {
        Self {
            v: [v0, v1, v2],
            n: [n0, n1, n2],
            uv: [uv0, uv1, uv2],
            material_id,
        }
    }
    /// Compute the geometric normal from vertex positions.
    pub fn geometric_normal(&self) -> [f64; 3] {
        let e1 = sub3(self.v[1], self.v[0]);
        let e2 = sub3(self.v[2], self.v[0]);
        normalize3(cross3(e1, e2))
    }
    /// Möller–Trumbore ray-triangle intersection.
    ///
    /// Returns `Some(t)` if the ray intersects the triangle, within `[t_min, t_max]`.
    pub fn intersect(&self, ray: &Ray) -> Option<f64> {
        let e1 = sub3(self.v[1], self.v[0]);
        let e2 = sub3(self.v[2], self.v[0]);
        let h = cross3(ray.direction, e2);
        let a = dot3(e1, h);
        if a.abs() < 1e-15 {
            return None;
        }
        let f = 1.0 / a;
        let s = sub3(ray.origin, self.v[0]);
        let u = f * dot3(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross3(s, e1);
        let v = f * dot3(ray.direction, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = f * dot3(e2, q);
        if t >= ray.t_min && t <= ray.t_max {
            Some(t)
        } else {
            None
        }
    }
    /// Full intersection returning barycentric coordinates.
    pub fn intersect_full(&self, ray: &Ray) -> Option<(f64, f64, f64)> {
        let e1 = sub3(self.v[1], self.v[0]);
        let e2 = sub3(self.v[2], self.v[0]);
        let h = cross3(ray.direction, e2);
        let a = dot3(e1, h);
        if a.abs() < 1e-15 {
            return None;
        }
        let f = 1.0 / a;
        let s = sub3(ray.origin, self.v[0]);
        let u = f * dot3(s, h);
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross3(s, e1);
        let v = f * dot3(ray.direction, q);
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let t = f * dot3(e2, q);
        if t >= ray.t_min && t <= ray.t_max {
            Some((t, u, v))
        } else {
            None
        }
    }
    /// Interpolate the shading normal at barycentric (u,v).
    pub fn interpolate_normal(&self, u: f64, v: f64) -> [f64; 3] {
        let w = 1.0 - u - v;
        let n = [
            w * self.n[0][0] + u * self.n[1][0] + v * self.n[2][0],
            w * self.n[0][1] + u * self.n[1][1] + v * self.n[2][1],
            w * self.n[0][2] + u * self.n[1][2] + v * self.n[2][2],
        ];
        normalize3(n)
    }
    /// Interpolate UV coordinates at barycentric (u,v).
    pub fn interpolate_uv(&self, u: f64, v: f64) -> [f64; 2] {
        let w = 1.0 - u - v;
        [
            w * self.uv[0][0] + u * self.uv[1][0] + v * self.uv[2][0],
            w * self.uv[0][1] + u * self.uv[1][1] + v * self.uv[2][1],
        ]
    }
    /// Axis-aligned bounding box of the triangle.
    pub fn aabb(&self) -> Aabb {
        let min_x = self.v[0][0].min(self.v[1][0]).min(self.v[2][0]);
        let min_y = self.v[0][1].min(self.v[1][1]).min(self.v[2][1]);
        let min_z = self.v[0][2].min(self.v[1][2]).min(self.v[2][2]);
        let max_x = self.v[0][0].max(self.v[1][0]).max(self.v[2][0]);
        let max_y = self.v[0][1].max(self.v[1][1]).max(self.v[2][1]);
        let max_z = self.v[0][2].max(self.v[1][2]).max(self.v[2][2]);
        Aabb {
            min: [min_x, min_y, min_z],
            max: [max_x, max_y, max_z],
        }
    }
    /// Centroid of the triangle.
    pub fn centroid(&self) -> [f64; 3] {
        [
            (self.v[0][0] + self.v[1][0] + self.v[2][0]) / 3.0,
            (self.v[0][1] + self.v[1][1] + self.v[2][1]) / 3.0,
            (self.v[0][2] + self.v[1][2] + self.v[2][2]) / 3.0,
        ]
    }
    /// Surface area of the triangle.
    pub fn area(&self) -> f64 {
        let e1 = sub3(self.v[1], self.v[0]);
        let e2 = sub3(self.v[2], self.v[0]);
        length3(cross3(e1, e2)) * 0.5
    }
}
/// Axis-aligned bounding box.
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
    /// Create an empty AABB (inverted extents).
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }
    /// Expand AABB to include a point.
    pub fn expand_point(&self, p: [f64; 3]) -> Self {
        Self {
            min: [
                self.min[0].min(p[0]),
                self.min[1].min(p[1]),
                self.min[2].min(p[2]),
            ],
            max: [
                self.max[0].max(p[0]),
                self.max[1].max(p[1]),
                self.max[2].max(p[2]),
            ],
        }
    }
    /// Merge two AABBs.
    pub fn merge(&self, other: &Self) -> Self {
        Self {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }
    /// Centroid of the AABB.
    pub fn centroid(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Surface area of the AABB.
    pub fn surface_area(&self) -> f64 {
        let d = [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ];
        2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
    }
    /// Longest axis (0=X, 1=Y, 2=Z).
    pub fn longest_axis(&self) -> usize {
        let d = [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ];
        if d[0] >= d[1] && d[0] >= d[2] {
            0
        } else if d[1] >= d[2] {
            1
        } else {
            2
        }
    }
    /// Ray-AABB intersection using slab method. Returns (t_near, t_far).
    pub fn intersect_ray(&self, ray: &Ray) -> Option<(f64, f64)> {
        let mut t_near = ray.t_min;
        let mut t_far = ray.t_max;
        for (&dir, (&min, (&max, &origin))) in ray
            .direction
            .iter()
            .zip(self.min.iter().zip(self.max.iter().zip(&ray.origin)))
        {
            let inv_d = if dir.abs() < 1e-15 {
                f64::INFINITY
            } else {
                1.0 / dir
            };
            let t0 = (min - origin) * inv_d;
            let t1 = (max - origin) * inv_d;
            let (t0, t1) = if inv_d < 0.0 { (t1, t0) } else { (t0, t1) };
            t_near = t_near.max(t0);
            t_far = t_far.min(t1);
            if t_far < t_near {
                return None;
            }
        }
        Some((t_near, t_far))
    }
}
/// Surface Area Heuristic (SAH) BVH over triangles.
#[derive(Debug, Clone)]
pub struct Bvh {
    /// All BVH nodes.
    pub nodes: Vec<BvhNode>,
    /// Primitive indices (reordered during build).
    pub prim_indices: Vec<u32>,
    /// Number of primitives.
    pub prim_count: usize,
}
impl Bvh {
    /// Build a BVH from a list of triangles using SAH.
    pub fn build(triangles: &[Triangle]) -> Self {
        let n = triangles.len();
        if n == 0 {
            return Self {
                nodes: Vec::new(),
                prim_indices: Vec::new(),
                prim_count: 0,
            };
        }
        let mut prim_indices: Vec<u32> = (0..n as u32).collect();
        let centroids: Vec<[f64; 3]> = triangles.iter().map(|t| t.centroid()).collect();
        let aabbs: Vec<Aabb> = triangles.iter().map(|t| t.aabb()).collect();
        let mut nodes = Vec::with_capacity(2 * n);
        let root_bounds = aabbs.iter().fold(Aabb::empty(), |acc, b| acc.merge(b));
        nodes.push(BvhNode {
            bounds: root_bounds,
            left_or_first: 0,
            prim_count: n as u32,
        });
        let mut stack = vec![0usize];
        while let Some(node_idx) = stack.pop() {
            let first = nodes[node_idx].left_or_first as usize;
            let count = nodes[node_idx].prim_count as usize;
            if count <= 4 {
                continue;
            }
            let parent_sa = nodes[node_idx].bounds.surface_area();
            let mut best_cost = f64::INFINITY;
            let mut best_axis = 0usize;
            let mut best_split = 0.0f64;
            for axis in [0usize, 1, 2] {
                let slice = &mut prim_indices[first..first + count];
                slice.sort_unstable_by(|&a, &b| {
                    centroids[a as usize][axis]
                        .partial_cmp(&centroids[b as usize][axis])
                        .expect("operation should succeed")
                });
                let mut left_bounds = Aabb::empty();
                let mut left_areas = Vec::with_capacity(count);
                for i in 0..count - 1 {
                    left_bounds = left_bounds.merge(&aabbs[slice[i] as usize]);
                    left_areas.push(left_bounds.surface_area());
                }
                let mut right_bounds = Aabb::empty();
                for i in (1..count).rev() {
                    right_bounds = right_bounds.merge(&aabbs[slice[i] as usize]);
                    let left_count = i;
                    let right_count = count - i;
                    let cost = (left_areas[i - 1] * left_count as f64
                        + right_bounds.surface_area() * right_count as f64)
                        / parent_sa;
                    if cost < best_cost {
                        best_cost = cost;
                        best_axis = axis;
                        best_split = centroids[slice[i] as usize][axis];
                    }
                }
            }
            let slice = &mut prim_indices[first..first + count];
            slice.sort_unstable_by(|&a, &b| {
                centroids[a as usize][best_axis]
                    .partial_cmp(&centroids[b as usize][best_axis])
                    .expect("operation should succeed")
            });
            let split_pos =
                slice.partition_point(|&idx| centroids[idx as usize][best_axis] < best_split);
            let split_pos = split_pos.clamp(1, count - 1);
            let left_count = split_pos;
            let right_count = count - split_pos;
            let left_bounds = prim_indices[first..first + left_count]
                .iter()
                .fold(Aabb::empty(), |acc, &i| acc.merge(&aabbs[i as usize]));
            let right_bounds = prim_indices[first + left_count..first + count]
                .iter()
                .fold(Aabb::empty(), |acc, &i| acc.merge(&aabbs[i as usize]));
            let left_child_idx = nodes.len();
            nodes.push(BvhNode {
                bounds: left_bounds,
                left_or_first: first as u32,
                prim_count: left_count as u32,
            });
            let right_child_idx = nodes.len();
            nodes.push(BvhNode {
                bounds: right_bounds,
                left_or_first: (first + left_count) as u32,
                prim_count: right_count as u32,
            });
            nodes[node_idx].left_or_first = left_child_idx as u32;
            nodes[node_idx].prim_count = 0;
            stack.push(left_child_idx);
            stack.push(right_child_idx);
        }
        Self {
            nodes,
            prim_indices,
            prim_count: n,
        }
    }
    /// Traverse the BVH and find the nearest intersection.
    pub fn intersect<'a>(
        &self,
        ray: &Ray,
        triangles: &'a [Triangle],
    ) -> Option<(HitRecord, &'a Triangle)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = Vec::with_capacity(64);
        stack.push(0usize);
        let mut best_t = ray.t_max;
        let mut best_hit: Option<(HitRecord, usize)> = None;
        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            let mut test_ray = *ray;
            test_ray.t_max = best_t;
            if node.bounds.intersect_ray(&test_ray).is_none() {
                continue;
            }
            if node.is_leaf() {
                let first = node.left_or_first as usize;
                let count = node.prim_count as usize;
                for i in first..first + count {
                    let tri_idx = self.prim_indices[i] as usize;
                    let tri = &triangles[tri_idx];
                    if let Some((t, u, v)) = tri.intersect_full(ray)
                        && t < best_t
                    {
                        best_t = t;
                        let pos = ray.at(t);
                        let norm = tri.interpolate_normal(u, v);
                        let uv = tri.interpolate_uv(u, v);
                        let hit = HitRecord::new(
                            t,
                            pos,
                            norm,
                            ray.direction,
                            uv,
                            tri_idx as u32,
                            tri.material_id,
                        );
                        best_hit = Some((hit, tri_idx));
                    }
                }
            } else {
                let left = node.left_or_first as usize;
                let right = left + 1;
                stack.push(left);
                stack.push(right);
            }
        }
        best_hit.map(|(hit, tri_idx)| (hit, &triangles[tri_idx]))
    }
    /// Test if any intersection exists (shadow ray query).
    pub fn intersect_any(&self, ray: &Ray, triangles: &[Triangle]) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        let mut stack = Vec::with_capacity(64);
        stack.push(0usize);
        while let Some(node_idx) = stack.pop() {
            let node = &self.nodes[node_idx];
            if node.bounds.intersect_ray(ray).is_none() {
                continue;
            }
            if node.is_leaf() {
                let first = node.left_or_first as usize;
                let count = node.prim_count as usize;
                for i in first..first + count {
                    let tri_idx = self.prim_indices[i] as usize;
                    if triangles[tri_idx].intersect(ray).is_some() {
                        return true;
                    }
                }
            } else {
                let left = node.left_or_first as usize;
                let right = left + 1;
                stack.push(left);
                stack.push(right);
            }
        }
        false
    }
}
/// Area light for soft shadow computation.
#[derive(Debug, Clone, Copy)]
pub struct AreaLight {
    /// Center position in world space.
    pub position: [f64; 3],
    /// Light color.
    pub color: [f64; 3],
    /// Intensity.
    pub intensity: f64,
    /// Light tangent (half-size along u axis).
    pub u_axis: [f64; 3],
    /// Light bitangent (half-size along v axis).
    pub v_axis: [f64; 3],
}
impl AreaLight {
    /// Create a new area light.
    pub fn new(
        position: [f64; 3],
        color: [f64; 3],
        intensity: f64,
        u_axis: [f64; 3],
        v_axis: [f64; 3],
    ) -> Self {
        Self {
            position,
            color,
            intensity,
            u_axis,
            v_axis,
        }
    }
    /// Sample a point on the light given stratified (su, sv) in \[-1,1\].
    pub fn sample_point(&self, su: f64, sv: f64) -> [f64; 3] {
        add3(
            add3(self.position, scale3(self.u_axis, su)),
            scale3(self.v_axis, sv),
        )
    }
}
/// A ray defined by an origin and normalized direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f64; 3],
    /// Normalized ray direction.
    pub direction: [f64; 3],
    /// Minimum valid hit distance.
    pub t_min: f64,
    /// Maximum valid hit distance.
    pub t_max: f64,
}
impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            origin,
            direction: normalize3(direction),
            t_min: 1e-4,
            t_max: f64::INFINITY,
        }
    }
    /// Evaluate the ray at parameter t: origin + t * direction.
    pub fn at(&self, t: f64) -> [f64; 3] {
        add3(self.origin, scale3(self.direction, t))
    }
}
/// Perspective camera model for primary ray generation.
#[derive(Debug, Clone)]
pub struct Camera {
    /// Camera position in world space.
    pub position: [f64; 3],
    /// Normalized forward direction.
    pub forward: [f64; 3],
    /// Normalized right direction.
    pub right: [f64; 3],
    /// Normalized up direction.
    pub up: [f64; 3],
    /// Vertical field of view in radians.
    pub fov_y: f64,
    /// Image aspect ratio (width/height).
    pub aspect: f64,
    /// Distance to the near plane.
    pub near: f64,
    /// Lens aperture radius (0 = pinhole).
    pub aperture: f64,
    /// Focus distance.
    pub focus_dist: f64,
}
impl Camera {
    /// Create a new perspective camera with look-at setup.
    pub fn look_at(
        eye: [f64; 3],
        target: [f64; 3],
        world_up: [f64; 3],
        fov_y_deg: f64,
        aspect: f64,
        aperture: f64,
        focus_dist: f64,
    ) -> Self {
        let forward = normalize3(sub3(target, eye));
        let right = normalize3(cross3(forward, world_up));
        let up = cross3(right, forward);
        Self {
            position: eye,
            forward,
            right,
            up,
            fov_y: fov_y_deg * PI / 180.0,
            aspect,
            near: 0.001,
            aperture,
            focus_dist,
        }
    }
    /// Generate a primary ray for pixel (px, py) in \[0,W) x \[0,H).
    pub fn generate_ray(&self, px: f64, py: f64, width: f64, height: f64) -> Ray {
        let half_h = (self.fov_y * 0.5).tan();
        let half_w = self.aspect * half_h;
        let ndc_x = (2.0 * (px + 0.5) / width - 1.0) * half_w;
        let ndc_y = (1.0 - 2.0 * (py + 0.5) / height) * half_h;
        let dir = normalize3(add3(
            add3(self.forward, scale3(self.right, ndc_x)),
            scale3(self.up, ndc_y),
        ));
        Ray::new(self.position, dir)
    }
    /// Generate a depth-of-field ray with lens sampling.
    pub fn generate_dof_ray(
        &self,
        px: f64,
        py: f64,
        width: f64,
        height: f64,
        lens_u: f64,
        lens_v: f64,
    ) -> Ray {
        let half_h = (self.fov_y * 0.5).tan();
        let half_w = self.aspect * half_h;
        let ndc_x = (2.0 * (px + 0.5) / width - 1.0) * half_w;
        let ndc_y = (1.0 - 2.0 * (py + 0.5) / height) * half_h;
        let focus_dir = normalize3(add3(
            add3(self.forward, scale3(self.right, ndc_x)),
            scale3(self.up, ndc_y),
        ));
        let focus_point = add3(self.position, scale3(focus_dir, self.focus_dist));
        let lens_offset = add3(
            scale3(self.right, lens_u * self.aperture),
            scale3(self.up, lens_v * self.aperture),
        );
        let origin = add3(self.position, lens_offset);
        let direction = normalize3(sub3(focus_point, origin));
        Ray::new(origin, direction)
    }
}
/// Material type enumeration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialType {
    /// Lambertian diffuse.
    Diffuse,
    /// Metal with roughness.
    Metal,
    /// Dielectric (glass/water).
    Dielectric,
    /// Emissive light source.
    Emissive,
    /// PBR material.
    Pbr,
}
/// A point light source.
#[derive(Debug, Clone, Copy)]
pub struct PointLight {
    /// Light position in world space.
    pub position: [f64; 3],
    /// Light color (linear RGB).
    pub color: [f64; 3],
    /// Light intensity.
    pub intensity: f64,
    /// Attenuation coefficients (constant, linear, quadratic).
    pub attenuation: [f64; 3],
}
impl PointLight {
    /// Create a new point light.
    pub fn new(position: [f64; 3], color: [f64; 3], intensity: f64) -> Self {
        Self {
            position,
            color,
            intensity,
            attenuation: [1.0, 0.0, 0.1],
        }
    }
    /// Compute attenuation at distance d.
    pub fn attenuate(&self, d: f64) -> f64 {
        let [c, l, q] = self.attenuation;
        1.0 / (c + l * d + q * d * d)
    }
}
/// Configuration for a full ray-trace render pass.
#[derive(Debug, Clone)]
pub struct RenderConfig {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Number of samples per pixel.
    pub spp: u32,
    /// Maximum path depth.
    pub max_depth: u32,
    /// Enable soft shadows.
    pub soft_shadows: bool,
    /// Enable ambient occlusion.
    pub ambient_occlusion: bool,
    /// Enable depth of field.
    pub depth_of_field: bool,
    /// Number of AO samples.
    pub ao_samples: u32,
    /// Number of shadow samples.
    pub shadow_samples: u32,
    /// Background color.
    pub background: [f64; 3],
    /// Ambient light color.
    pub ambient: [f64; 3],
    /// Tone mapping mode (0=none, 1=reinhard, 2=filmic, 3=ACES).
    pub tonemap: u32,
}
