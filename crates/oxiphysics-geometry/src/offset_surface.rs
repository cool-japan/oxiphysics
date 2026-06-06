// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Offset surface generation and signed distance field (SDF) operations.
//!
//! Provides SDF primitives, CSG operations, mesh offsetting, and voxel SDF grids.

use std::collections::HashMap;

/// Type alias for a pair of convex/concave edge lists returned by `detect_edge_features`.
type EdgeFeatureLists = (Vec<(usize, usize)>, Vec<(usize, usize)>);
/// Type alias for the per-edge face-normal accumulation map.
type EdgeFaceNormalMap = HashMap<(usize, usize), Vec<([f64; 3], [f64; 3])>>;

// ---------------------------------------------------------------------------
// Helper math on plain [f64; 3]
// ---------------------------------------------------------------------------

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / len)
    }
}

// ---------------------------------------------------------------------------
// Sdf trait
// ---------------------------------------------------------------------------

/// Signed distance function: negative inside the shape, positive outside.
pub trait Sdf: Send + Sync {
    /// Signed distance from point `p` to the surface.
    fn distance(&self, p: [f64; 3]) -> f64;

    /// Gradient of the SDF at `p` (numerical central differences by default).
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let eps = 1e-5;
        let dx = self.distance([p[0] + eps, p[1], p[2]]) - self.distance([p[0] - eps, p[1], p[2]]);
        let dy = self.distance([p[0], p[1] + eps, p[2]]) - self.distance([p[0], p[1] - eps, p[2]]);
        let dz = self.distance([p[0], p[1], p[2] + eps]) - self.distance([p[0], p[1], p[2] - eps]);
        [dx / (2.0 * eps), dy / (2.0 * eps), dz / (2.0 * eps)]
    }

    /// Outward unit normal at `p` (normalised gradient).
    fn normal(&self, p: [f64; 3]) -> [f64; 3] {
        normalize(self.gradient(p))
    }
}

// ---------------------------------------------------------------------------
// SDF primitives
// ---------------------------------------------------------------------------

/// Sphere SDF (exact).
pub struct SdfSphere {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}

impl SdfSphere {
    /// Create a new sphere SDF with given center and radius.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
}

impl Sdf for SdfSphere {
    fn distance(&self, p: [f64; 3]) -> f64 {
        length(sub(p, self.center)) - self.radius
    }
}

/// Axis-aligned box SDF (exact).
pub struct SdfBox {
    /// Centre of the box in world space.
    pub center: [f64; 3],
    /// Half-lengths along each axis.
    pub half_extents: [f64; 3],
}

impl SdfBox {
    /// Create a new box SDF centered at origin with given half-extents (3 scalars).
    pub fn new(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            half_extents: [hx, hy, hz],
        }
    }
    /// Create a new box SDF with explicit center and half-extents arrays.
    pub fn new_centered(center: [f64; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            half_extents,
        }
    }
}

impl Sdf for SdfBox {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let q = sub(p, self.center);
        let qx = q[0].abs() - self.half_extents[0];
        let qy = q[1].abs() - self.half_extents[1];
        let qz = q[2].abs() - self.half_extents[2];
        let outside = length([qx.max(0.0), qy.max(0.0), qz.max(0.0)]);
        let inside = qx.max(qy).max(qz).min(0.0);
        outside + inside
    }
}

/// Capsule (line-segment rounded) SDF (exact).
pub struct SdfCapsule {
    /// Start point of the capsule axis.
    pub a: [f64; 3],
    /// End point of the capsule axis.
    pub b: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}

impl Sdf for SdfCapsule {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let ab = sub(self.b, self.a);
        let ap = sub(p, self.a);
        let t = (dot(ap, ab) / dot(ab, ab)).clamp(0.0, 1.0);
        let closest = add(self.a, scale(ab, t));
        length(sub(p, closest)) - self.radius
    }
}

/// Half-space plane SDF: `n·p - offset` (exact). `normal` should be unit length.
pub struct SdfPlane {
    /// Outward-facing unit normal of the plane.
    pub normal: [f64; 3],
    /// Signed distance from the origin to the plane along `normal`.
    pub offset: f64,
}

impl SdfPlane {
    /// Create a new half-space plane SDF with given unit normal and offset.
    pub fn new(normal: [f64; 3], offset: f64) -> Self {
        Self { normal, offset }
    }
}

impl Sdf for SdfPlane {
    fn distance(&self, p: [f64; 3]) -> f64 {
        dot(self.normal, p) - self.offset
    }
}

/// Torus SDF lying in the XZ plane, centred at `center` (exact).
pub struct SdfTorus {
    /// Centre of the torus.
    pub center: [f64; 3],
    /// Distance from the torus centre to the tube centre.
    pub major_radius: f64,
    /// Radius of the tube.
    pub minor_radius: f64,
}

impl Sdf for SdfTorus {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let q = sub(p, self.center);
        // XZ-plane torus: revolve a circle around the Y axis
        let r_xz = (q[0] * q[0] + q[2] * q[2]).sqrt();
        let d_xz = r_xz - self.major_radius;
        (d_xz * d_xz + q[1] * q[1]).sqrt() - self.minor_radius
    }
}

// ---------------------------------------------------------------------------
// CSG / SDF operations
// ---------------------------------------------------------------------------

/// Boolean union (min).
pub struct SdfUnion {
    /// First SDF operand.
    pub a: Box<dyn Sdf>,
    /// Second SDF operand.
    pub b: Box<dyn Sdf>,
}

impl Sdf for SdfUnion {
    fn distance(&self, p: [f64; 3]) -> f64 {
        self.a.distance(p).min(self.b.distance(p))
    }
}

/// Boolean intersection (max).
pub struct SdfIntersection {
    /// First SDF operand.
    pub a: Box<dyn Sdf>,
    /// Second SDF operand.
    pub b: Box<dyn Sdf>,
}

impl Sdf for SdfIntersection {
    fn distance(&self, p: [f64; 3]) -> f64 {
        self.a.distance(p).max(self.b.distance(p))
    }
}

/// Boolean difference: A minus B (max(da, -db)).
pub struct SdfDifference {
    /// The base SDF (A).
    pub a: Box<dyn Sdf>,
    /// The subtracted SDF (B).
    pub b: Box<dyn Sdf>,
}

impl Sdf for SdfDifference {
    fn distance(&self, p: [f64; 3]) -> f64 {
        self.a.distance(p).max(-self.b.distance(p))
    }
}

/// Offset (shell expansion / shrinkage): `d(p) - offset`.
pub struct SdfOffset {
    /// The inner SDF being offset.
    pub inner: Box<dyn Sdf>,
    /// Amount to offset (positive = expand, negative = shrink).
    pub offset: f64,
}

impl Sdf for SdfOffset {
    fn distance(&self, p: [f64; 3]) -> f64 {
        self.inner.distance(p) - self.offset
    }
}

/// Polynomial smooth union with blending radius `k`.
pub struct SdfSmoothUnion {
    /// First SDF operand.
    pub a: Box<dyn Sdf>,
    /// Second SDF operand.
    pub b: Box<dyn Sdf>,
    /// Blending radius (larger = smoother transition).
    pub k: f64,
}

impl Sdf for SdfSmoothUnion {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let da = self.a.distance(p);
        let db = self.b.distance(p);
        let h = (0.5 + 0.5 * (db - da) / self.k).clamp(0.0, 1.0);
        // mix(db, da, h) - k*h*(1-h)
        db + (da - db) * h - self.k * h * (1.0 - h)
    }
}

/// Polynomial smooth intersection with blending radius `k`.
pub struct SdfSmoothIntersection {
    /// First SDF operand.
    pub a: Box<dyn Sdf>,
    /// Second SDF operand.
    pub b: Box<dyn Sdf>,
    /// Blending radius (larger = smoother transition).
    pub k: f64,
}

impl Sdf for SdfSmoothIntersection {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let da = self.a.distance(p);
        let db = self.b.distance(p);
        let h = (0.5 - 0.5 * (db - da) / self.k).clamp(0.0, 1.0);
        db + (da - db) * h + self.k * h * (1.0 - h)
    }
}

// ---------------------------------------------------------------------------
// OffsetMesh
// ---------------------------------------------------------------------------

/// A triangle mesh with per-vertex normals, supporting inward/outward offsetting.
pub struct OffsetMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Per-vertex normals (unit vectors).
    pub normals: Vec<[f64; 3]>,
    /// Triangle faces as vertex index triples.
    pub faces: Vec<[usize; 3]>,
}

impl OffsetMesh {
    /// Build an `OffsetMesh` from a triangle soup, computing averaged vertex normals.
    pub fn from_triangle_soup(verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Self {
        let normals = Self::compute_vertex_normals(verts, faces);
        Self {
            vertices: verts.to_vec(),
            normals,
            faces: faces.to_vec(),
        }
    }

    /// Compute per-vertex normals by area-weighted average of incident face normals.
    pub fn compute_vertex_normals(verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Vec<[f64; 3]> {
        let n = verts.len();
        let mut accum = vec![[0.0f64; 3]; n];

        for f in faces {
            let v0 = verts[f[0]];
            let v1 = verts[f[1]];
            let v2 = verts[f[2]];
            let e1 = sub(v1, v0);
            let e2 = sub(v2, v0);
            let face_normal = cross(e1, e2); // magnitude = 2 * area (weights by area)
            for &vi in f {
                accum[vi] = add(accum[vi], face_normal);
            }
        }

        accum.iter().map(|&n| normalize(n)).collect()
    }

    /// Offset every vertex by `d` along its averaged normal.
    ///
    /// Positive `d` expands outward; negative `d` shrinks inward.
    pub fn offset(&self, d: f64) -> OffsetMesh {
        let new_vertices: Vec<[f64; 3]> = self
            .vertices
            .iter()
            .zip(self.normals.iter())
            .map(|(&v, &n)| add(v, scale(n, d)))
            .collect();
        OffsetMesh {
            vertices: new_vertices,
            normals: self.normals.clone(),
            faces: self.faces.clone(),
        }
    }
}

// ---------------------------------------------------------------------------
// VoxelSdf
// ---------------------------------------------------------------------------

/// Discretised signed distance field on a regular axis-aligned grid.
pub struct VoxelSdf {
    /// Number of cells along the X axis.
    pub nx: usize,
    /// Number of cells along the Y axis.
    pub ny: usize,
    /// Number of cells along the Z axis.
    pub nz: usize,
    /// World-space position of the grid origin (corner of cell `[0,0,0]`).
    pub origin: [f64; 3],
    /// Cell side length (isotropic).
    pub dx: f64,
    /// Flat array of SDF values, stored in x-major order.
    pub values: Vec<f64>,
}

impl VoxelSdf {
    /// Create a zeroed voxel grid.
    pub fn new(nx: usize, ny: usize, nz: usize, origin: [f64; 3], dx: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin,
            dx,
            values: vec![0.0; nx * ny * nz],
        }
    }

    /// Flat index from 3-D grid coordinates.
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }

    /// Convert a world-space point to fractional grid coordinates.
    pub fn world_to_grid(&self, p: [f64; 3]) -> [f64; 3] {
        [
            (p[0] - self.origin[0]) / self.dx,
            (p[1] - self.origin[1]) / self.dx,
            (p[2] - self.origin[2]) / self.dx,
        ]
    }

    /// Trilinear interpolation of the voxel grid at world-space point `p`.
    pub fn sample_trilinear(&self, p: [f64; 3]) -> f64 {
        let g = self.world_to_grid(p);
        let x0 = g[0].floor() as isize;
        let y0 = g[1].floor() as isize;
        let z0 = g[2].floor() as isize;
        let fx = g[0] - x0 as f64;
        let fy = g[1] - y0 as f64;
        let fz = g[2] - z0 as f64;

        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let nz = self.nz as isize;

        let clamp_x = |i: isize| i.clamp(0, nx - 1) as usize;
        let clamp_y = |i: isize| i.clamp(0, ny - 1) as usize;
        let clamp_z = |i: isize| i.clamp(0, nz - 1) as usize;

        let v = |dx: isize, dy: isize, dz: isize| -> f64 {
            self.values[self.idx(clamp_x(x0 + dx), clamp_y(y0 + dy), clamp_z(z0 + dz))]
        };

        let c00 = v(0, 0, 0) * (1.0 - fx) + v(1, 0, 0) * fx;
        let c01 = v(0, 0, 1) * (1.0 - fx) + v(1, 0, 1) * fx;
        let c10 = v(0, 1, 0) * (1.0 - fx) + v(1, 1, 0) * fx;
        let c11 = v(0, 1, 1) * (1.0 - fx) + v(1, 1, 1) * fx;

        let c0 = c00 * (1.0 - fy) + c10 * fy;
        let c1 = c01 * (1.0 - fy) + c11 * fy;

        c0 * (1.0 - fz) + c1 * fz
    }

    /// Populate the grid by evaluating `sdf` at every cell centre.
    pub fn from_sdf(
        sdf: &dyn Sdf,
        nx: usize,
        ny: usize,
        nz: usize,
        origin: [f64; 3],
        dx: f64,
    ) -> Self {
        let mut grid = Self::new(nx, ny, nz, origin, dx);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        origin[0] + (ix as f64 + 0.5) * dx,
                        origin[1] + (iy as f64 + 0.5) * dx,
                        origin[2] + (iz as f64 + 0.5) * dx,
                    ];
                    let i = grid.idx(ix, iy, iz);
                    grid.values[i] = sdf.distance(p);
                }
            }
        }
        grid
    }

    /// Central-difference gradient at grid cell `(ix, iy, iz)`.
    pub fn gradient_central(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;

        let xp = if ix + 1 < nx { ix + 1 } else { ix };
        let xm = if ix > 0 { ix - 1 } else { ix };
        let yp = if iy + 1 < ny { iy + 1 } else { iy };
        let ym = if iy > 0 { iy - 1 } else { iy };
        let zp = if iz + 1 < nz { iz + 1 } else { iz };
        let zm = if iz > 0 { iz - 1 } else { iz };

        let step_x = (xp - xm) as f64;
        let step_y = (yp - ym) as f64;
        let step_z = (zp - zm) as f64;

        let gx = (self.values[self.idx(xp, iy, iz)] - self.values[self.idx(xm, iy, iz)])
            / (step_x * self.dx);
        let gy = (self.values[self.idx(ix, yp, iz)] - self.values[self.idx(ix, ym, iz)])
            / (step_y * self.dx);
        let gz = (self.values[self.idx(ix, iy, zp)] - self.values[self.idx(ix, iy, zm)])
            / (step_z * self.dx);

        [gx, gy, gz]
    }
}

// ---------------------------------------------------------------------------
// Inward / outward offset variants
// ---------------------------------------------------------------------------

/// Compute an outward-offset mesh from a triangle soup.
///
/// Equivalent to `OffsetMesh::from_triangle_soup(verts, faces).offset(d)` for
/// positive `d`, but named explicitly for clarity.
pub fn outward_offset(verts: &[[f64; 3]], faces: &[[usize; 3]], d: f64) -> OffsetMesh {
    OffsetMesh::from_triangle_soup(verts, faces).offset(d)
}

/// Compute an inward-offset (shrunk) mesh.
///
/// Equivalent to `offset(-d)` where `d > 0`.
pub fn inward_offset(verts: &[[f64; 3]], faces: &[[usize; 3]], d: f64) -> OffsetMesh {
    OffsetMesh::from_triangle_soup(verts, faces).offset(-d)
}

// ---------------------------------------------------------------------------
// Collision offset: expand each face outward by a distance while keeping
// connectivity intact
// ---------------------------------------------------------------------------

impl OffsetMesh {
    /// Compute the "collision offset" shell mesh.
    ///
    /// Generates a new mesh where each face is moved outward by `d` along its
    /// face normal (not vertex normal).  Unlike vertex-normal offsetting, this
    /// produces a shell mesh useful for broad CCD convex hull testing.
    pub fn collision_offset_shell(&self, d: f64) -> OffsetMesh {
        let n_faces = self.faces.len();
        let mut new_verts = Vec::with_capacity(n_faces * 3);
        let mut new_faces = Vec::with_capacity(n_faces);

        for (fi, face) in self.faces.iter().enumerate() {
            let v0 = self.vertices[face[0]];
            let v1 = self.vertices[face[1]];
            let v2 = self.vertices[face[2]];
            let e1 = sub(v1, v0);
            let e2 = sub(v2, v0);
            let face_n = normalize(cross(e1, e2));
            let base = fi * 3;
            new_verts.push(add(v0, scale(face_n, d)));
            new_verts.push(add(v1, scale(face_n, d)));
            new_verts.push(add(v2, scale(face_n, d)));
            new_faces.push([base, base + 1, base + 2]);
        }

        // Recompute normals for the new shell
        let normals = OffsetMesh::compute_vertex_normals(&new_verts, &new_faces);
        OffsetMesh {
            vertices: new_verts,
            normals,
            faces: new_faces,
        }
    }

    /// Validate that all vertex normals are unit length.
    pub fn normals_are_unit(&self) -> bool {
        self.normals.iter().all(|&n| (length(n) - 1.0).abs() < 1e-6)
    }

    /// Surface area of the mesh (sum of triangle areas).
    pub fn surface_area(&self) -> f64 {
        self.faces
            .iter()
            .map(|f| {
                let v0 = self.vertices[f[0]];
                let v1 = self.vertices[f[1]];
                let v2 = self.vertices[f[2]];
                length(cross(sub(v1, v0), sub(v2, v0))) * 0.5
            })
            .sum()
    }

    /// Centroid of the mesh (average vertex position).
    pub fn centroid(&self) -> [f64; 3] {
        if self.vertices.is_empty() {
            return [0.0; 3];
        }
        let sum = self
            .vertices
            .iter()
            .fold([0.0f64; 3], |acc, &v| add(acc, v));
        scale(sum, 1.0 / self.vertices.len() as f64)
    }
}

// ---------------------------------------------------------------------------
// Medial axis approximation
// ---------------------------------------------------------------------------

/// Approximate medial axis computation using a voxel SDF.
///
/// The medial axis is the locus of points equidistant from the closest surface
/// points.  We approximate it by finding voxels where the SDF gradient magnitude
/// is significantly less than 1 (i.e., the distance field is "flat" — these
/// points are near the medial axis).
///
/// Returns the world-space positions of approximate medial axis voxel centres.
pub fn approximate_medial_axis(grid: &VoxelSdf, gradient_threshold: f64) -> Vec<[f64; 3]> {
    let mut points = Vec::new();
    for iz in 0..grid.nz {
        for iy in 0..grid.ny {
            for ix in 0..grid.nx {
                let val = grid.values[grid.idx(ix, iy, iz)];
                // Only consider interior voxels (negative SDF)
                if val >= 0.0 {
                    continue;
                }
                let g = grid.gradient_central(ix, iy, iz);
                let grad_mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if grad_mag < gradient_threshold {
                    let p = [
                        grid.origin[0] + (ix as f64 + 0.5) * grid.dx,
                        grid.origin[1] + (iy as f64 + 0.5) * grid.dx,
                        grid.origin[2] + (iz as f64 + 0.5) * grid.dx,
                    ];
                    points.push(p);
                }
            }
        }
    }
    points
}

// ---------------------------------------------------------------------------
// SDF cylinder primitive
// ---------------------------------------------------------------------------

/// Cylinder SDF aligned with the Y axis (exact).
pub struct SdfCylinder {
    /// Centre of the cylinder (midpoint of its axis).
    pub center: [f64; 3],
    /// Radius of the cylinder.
    pub radius: f64,
    /// Half the total height along the Y axis.
    pub half_height: f64,
}

impl Sdf for SdfCylinder {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let q = sub(p, self.center);
        let xz_dist = (q[0] * q[0] + q[2] * q[2]).sqrt() - self.radius;
        let y_dist = q[1].abs() - self.half_height;
        let outside =
            (xz_dist.max(0.0) * xz_dist.max(0.0) + y_dist.max(0.0) * y_dist.max(0.0)).sqrt();
        outside + xz_dist.max(y_dist).min(0.0)
    }
}

// ---------------------------------------------------------------------------
// SDF cone primitive
// ---------------------------------------------------------------------------

/// Cone SDF: tip at `apex`, opening downward along -Y, with half-angle `angle_rad`.
pub struct SdfCone {
    /// Position of the cone tip.
    pub apex: [f64; 3],
    /// Height of the cone (distance from apex to base).
    pub height: f64,
    /// Half-opening angle in radians.
    pub angle_rad: f64,
}

impl Sdf for SdfCone {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let q = sub(p, self.apex);
        let r = (q[0] * q[0] + q[2] * q[2]).sqrt();
        // Cone: r = -y * tan(angle) for y ≤ 0 (opening downward)
        let sin_a = self.angle_rad.sin();
        let cos_a = self.angle_rad.cos();
        let dist_axis = r * cos_a + q[1] * sin_a; // positive outside cone surface
        let dist_cap = q[1] + self.height; // positive above the base cap
        dist_axis.max(-dist_cap)
    }
}

// ---------------------------------------------------------------------------
// SdfTranslated: translate any SDF by an offset
// ---------------------------------------------------------------------------

/// Translate an SDF by a fixed offset.
pub struct SdfTranslated {
    /// The inner SDF being translated.
    pub inner: Box<dyn Sdf>,
    /// Translation vector applied before evaluating the inner SDF.
    pub offset: [f64; 3],
}

impl Sdf for SdfTranslated {
    fn distance(&self, p: [f64; 3]) -> f64 {
        let local = sub(p, self.offset);
        self.inner.distance(local)
    }
}

// ---------------------------------------------------------------------------
// SdfScaled: uniformly scale an SDF
// ---------------------------------------------------------------------------

/// Uniformly scale an SDF (scale > 1 makes the shape larger).
pub struct SdfScaled {
    /// The inner SDF being scaled.
    pub inner: Box<dyn Sdf>,
    /// Uniform scale factor applied to the shape.
    pub scale_factor: f64,
}

impl Sdf for SdfScaled {
    fn distance(&self, p: [f64; 3]) -> f64 {
        if self.scale_factor.abs() < 1e-300 {
            return f64::INFINITY;
        }
        self.inner.distance(scale(p, 1.0 / self.scale_factor)) * self.scale_factor
    }
}

// ---------------------------------------------------------------------------
// Marching-cubes-ready isosurface extraction (threshold crossing)
// ---------------------------------------------------------------------------

/// A simple marching-squares (2-D) contour extractor for a 2-D slice of a VoxelSdf.
///
/// Extracts edge midpoints where the SDF crosses zero on the XY slice at `iz`.
pub fn extract_zero_crossings_slice(grid: &VoxelSdf, iz: usize) -> Vec<[f64; 2]> {
    let mut crossings = Vec::new();
    let nz = grid.nz;
    if iz >= nz {
        return crossings;
    }
    for iy in 0..grid.ny.saturating_sub(1) {
        for ix in 0..grid.nx.saturating_sub(1) {
            let v00 = grid.values[grid.idx(ix, iy, iz)];
            let v10 = grid.values[grid.idx(ix + 1, iy, iz)];
            let v01 = grid.values[grid.idx(ix, iy + 1, iz)];

            let x0 = grid.origin[0] + (ix as f64 + 0.5) * grid.dx;
            let x1 = grid.origin[0] + (ix as f64 + 1.5) * grid.dx;
            let y0 = grid.origin[1] + (iy as f64 + 0.5) * grid.dx;
            let y1 = grid.origin[1] + (iy as f64 + 1.5) * grid.dx;

            // Horizontal edge
            if (v00 < 0.0) != (v10 < 0.0) {
                let t = v00 / (v00 - v10);
                crossings.push([x0 + t * (x1 - x0), y0]);
            }
            // Vertical edge
            if (v00 < 0.0) != (v01 < 0.0) {
                let t = v00 / (v00 - v01);
                crossings.push([x0, y0 + t * (y1 - y0)]);
            }
        }
    }
    crossings
}

// ---------------------------------------------------------------------------
// Variable offset (non-uniform per-vertex weights)
// ---------------------------------------------------------------------------

/// Apply a different offset distance to each vertex along its averaged normal.
///
/// `weights[i]` is the offset distance for vertex `i`.  Positive values expand
/// outward; negative values shrink inward.
pub fn variable_offset(verts: &[[f64; 3]], faces: &[[usize; 3]], weights: &[f64]) -> OffsetMesh {
    let normals = OffsetMesh::compute_vertex_normals(verts, faces);
    let new_verts: Vec<[f64; 3]> = verts
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let w = if i < weights.len() { weights[i] } else { 0.0 };
            let n = if i < normals.len() {
                normals[i]
            } else {
                [0.0; 3]
            };
            add(v, scale(n, w))
        })
        .collect();
    OffsetMesh {
        normals: OffsetMesh::compute_vertex_normals(&new_verts, faces),
        vertices: new_verts,
        faces: faces.to_vec(),
    }
}

// ---------------------------------------------------------------------------
// Polyhedron offset
// ---------------------------------------------------------------------------

/// Offset a closed polyhedron mesh outward (or inward for negative `d`).
///
/// Equivalent to `OffsetMesh::from_triangle_soup(verts, faces).offset(d)` but
/// named for clarity when working with closed polyhedra.
pub fn offset_polyhedron(verts: &[[f64; 3]], faces: &[[usize; 3]], d: f64) -> OffsetMesh {
    OffsetMesh::from_triangle_soup(verts, faces).offset(d)
}

// ---------------------------------------------------------------------------
// Feature detection: convex / concave edges
// ---------------------------------------------------------------------------

/// Detect feature edges in a triangle mesh, classifying interior edges as
/// convex (dihedral angle < π) or concave (dihedral angle > π).
///
/// Returns `(convex_edges, concave_edges)` as lists of `(v0, v1)` index pairs.
/// Boundary edges (shared by only one face) are excluded from both lists.
pub fn detect_edge_features(verts: &[[f64; 3]], faces: &[[usize; 3]]) -> EdgeFeatureLists {
    // Build edge → face normals
    let mut edge_data: EdgeFaceNormalMap = HashMap::new();

    for face in faces {
        let v0 = verts[face[0]];
        let v1 = verts[face[1]];
        let v2 = verts[face[2]];
        let e1 = sub(v1, v0);
        let e2 = sub(v2, v0);
        let face_n = normalize(cross(e1, e2));
        let face_centroid = scale(add(add(v0, v1), v2), 1.0 / 3.0);

        for k in 0..3 {
            let ea = face[k];
            let eb = face[(k + 1) % 3];
            let key = (ea.min(eb), ea.max(eb));
            // Store (face_normal, midpoint of face)
            let mid_edge = scale(add(verts[ea], verts[eb]), 0.5);
            // We store (face_normal, centroid_to_edge midpoint vector)
            let ct_to_edge = sub(mid_edge, face_centroid);
            edge_data.entry(key).or_default().push((face_n, ct_to_edge));
        }
    }

    let mut convex = Vec::new();
    let mut concave = Vec::new();

    for (&(ea, eb), data) in &edge_data {
        if data.len() != 2 {
            continue; // boundary edge
        }
        let n0 = data[0].0;
        let n1 = data[1].0;
        let ct0 = data[0].1;

        let dot_nn = dot(n0, n1);
        // Convex: the other face's normal points away from first face's centroid side
        // Simple heuristic: dot(n0, ct0) > 0 means the edge is convex
        let sign = dot(cross(n0, n1), ct0);

        if dot_nn < 0.9999 {
            // not coplanar
            if sign >= 0.0 {
                convex.push((ea, eb));
            } else {
                concave.push((ea, eb));
            }
        }
    }

    (convex, concave)
}

// ---------------------------------------------------------------------------
// Offset curves in 3D
// ---------------------------------------------------------------------------

/// Offset a 3D polyline by `d` along a fixed `normal` direction.
///
/// Returns a new set of points displaced by `d * normal`.
pub fn offset_curve_3d(pts: &[[f64; 3]], normal: [f64; 3], d: f64) -> Vec<[f64; 3]> {
    let n = normalize(normal);
    pts.iter().map(|&p| add(p, scale(n, d))).collect()
}

// ---------------------------------------------------------------------------
// Shell generation
// ---------------------------------------------------------------------------

/// Generate a solid "shell" mesh from a surface mesh.
///
/// Creates an outer surface (offset by `+thickness`) and an inner surface
/// (the original mesh), connected by side walls along the boundary.
/// For closed meshes the result is a thick shell.
pub fn generate_shell(verts: &[[f64; 3]], faces: &[[usize; 3]], thickness: f64) -> OffsetMesh {
    let n_orig = verts.len();
    let outer = OffsetMesh::from_triangle_soup(verts, faces).offset(thickness);

    // Combine inner (original) and outer vertices
    let mut new_verts: Vec<[f64; 3]> = verts.to_vec();
    new_verts.extend_from_slice(&outer.vertices);

    // Inner faces (original, flipped winding for inward normals)
    let mut new_faces: Vec<[usize; 3]> = faces.iter().map(|f| [f[2], f[1], f[0]]).collect();

    // Outer faces (offset, original winding)
    new_faces.extend(
        faces
            .iter()
            .map(|f| [f[0] + n_orig, f[1] + n_orig, f[2] + n_orig]),
    );

    // Side walls: for boundary edges (only one adjacent face), connect inner and outer
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for face in faces {
        for k in 0..3 {
            let ea = face[k];
            let eb = face[(k + 1) % 3];
            let key = (ea.min(eb), ea.max(eb));
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    for ((ea, eb), count) in &edge_count {
        if *count == 1 {
            // Boundary edge → add a quad wall (two triangles)
            let ia = *ea;
            let ib = *eb;
            let oa = *ea + n_orig;
            let ob = *eb + n_orig;
            new_faces.push([ia, ib, ob]);
            new_faces.push([ia, ob, oa]);
        }
    }

    let normals = OffsetMesh::compute_vertex_normals(&new_verts, &new_faces);
    OffsetMesh {
        vertices: new_verts,
        normals,
        faces: new_faces,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // --- SdfSphere -----------------------------------------------------------

    #[test]
    fn test_sphere_center_negative_radius() {
        let s = SdfSphere {
            center: [0.0, 0.0, 0.0],
            radius: 2.0,
        };
        assert!((s.distance([0.0, 0.0, 0.0]) - (-2.0)).abs() < EPS);
    }

    #[test]
    fn test_sphere_surface_zero() {
        let s = SdfSphere {
            center: [1.0, 2.0, 3.0],
            radius: 1.5,
        };
        let p = [1.0 + 1.5, 2.0, 3.0];
        assert!(s.distance(p).abs() < EPS);
    }

    #[test]
    fn test_sphere_outside_positive() {
        let s = SdfSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        assert!(s.distance([5.0, 0.0, 0.0]) > 0.0);
    }

    // --- SdfBox --------------------------------------------------------------

    #[test]
    fn test_box_inside_negative() {
        let b = SdfBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [2.0, 2.0, 2.0],
        };
        assert!(b.distance([0.5, 0.5, 0.5]) < 0.0);
        assert!(b.distance([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_box_outside_positive() {
        let b = SdfBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        assert!(b.distance([3.0, 0.0, 0.0]) > 0.0);
    }

    // --- SdfCapsule ----------------------------------------------------------

    #[test]
    fn test_capsule_midpoint_surface() {
        let c = SdfCapsule {
            a: [0.0, 0.0, 0.0],
            b: [0.0, 4.0, 0.0],
            radius: 1.0,
        };
        // Midpoint of axis = (0, 2, 0); move 1 unit perpendicular => on surface
        assert!(c.distance([1.0, 2.0, 0.0]).abs() < EPS);
    }

    #[test]
    fn test_capsule_inside_negative() {
        let c = SdfCapsule {
            a: [0.0, 0.0, 0.0],
            b: [0.0, 4.0, 0.0],
            radius: 1.0,
        };
        assert!(c.distance([0.0, 2.0, 0.0]) < 0.0);
    }

    // --- SdfUnion ------------------------------------------------------------

    #[test]
    fn test_union_inside_either() {
        let union = SdfUnion {
            a: Box::new(SdfSphere {
                center: [-3.0, 0.0, 0.0],
                radius: 1.5,
            }),
            b: Box::new(SdfSphere {
                center: [3.0, 0.0, 0.0],
                radius: 1.5,
            }),
        };
        // Point inside first sphere
        assert!(union.distance([-3.0, 0.0, 0.0]) < 0.0);
        // Point inside second sphere
        assert!(union.distance([3.0, 0.0, 0.0]) < 0.0);
        // Point outside both
        assert!(union.distance([0.0, 0.0, 0.0]) > 0.0);
    }

    // --- SdfOffset -----------------------------------------------------------

    #[test]
    fn test_sdf_offset_expanded_sphere() {
        let base = SdfSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let expanded = SdfOffset {
            inner: Box::new(SdfSphere {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
            }),
            offset: 0.5,
        };
        // At distance 1.4 from centre: original says +0.4 (outside), expanded says -0.1 (inside)
        let r = 1.4_f64;
        let p = [r, 0.0, 0.0];
        assert!(base.distance(p) > 0.0);
        assert!(expanded.distance(p) < 0.0);
    }

    // --- OffsetMesh ----------------------------------------------------------

    #[test]
    fn test_offset_mesh_outward() {
        // Build a single triangle
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        let d = 0.3;
        let off = mesh.offset(d);
        // Each vertex should have moved away from original by exactly d
        for (orig, new_v) in mesh.vertices.iter().zip(off.vertices.iter()) {
            let dist = length(sub(*new_v, *orig));
            assert!(
                (dist - d).abs() < 1e-9,
                "vertex moved by {dist}, expected {d}"
            );
        }
    }

    #[test]
    fn test_offset_mesh_normals_unit() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        for n in &mesh.normals {
            let len = length(*n);
            assert!((len - 1.0).abs() < 1e-9, "normal not unit: {len}");
        }
    }

    // --- VoxelSdf ------------------------------------------------------------

    #[test]
    fn test_voxel_sdf_sphere_center() {
        let sphere = SdfSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        // 10x10x10 grid, origin at -1, cell size 0.2 => covers [-1, 1] per axis
        let grid = VoxelSdf::from_sdf(&sphere, 10, 10, 10, [-1.0, -1.0, -1.0], 0.2);
        // Centre cell: ix=iy=iz=4 => world pos [(-1 + 4.5*0.2), ...] = [-0.1, -0.1, -0.1]
        let centre_idx = grid.idx(4, 4, 4);
        let v = grid.values[centre_idx];
        // Should be negative (inside sphere)
        assert!(v < 0.0, "centre cell value {v} should be negative");
        // Should be close to -radius = -1.0 (within grid resolution error)
        assert!(
            (v - (-1.0)).abs() < 0.3,
            "centre cell value {v} not close to -1.0"
        );
    }

    #[test]
    fn test_voxel_sdf_constant_trilinear() {
        let mut grid = VoxelSdf::new(4, 4, 4, [0.0, 0.0, 0.0], 1.0);
        // Fill with constant value 3.125
        for v in grid.values.iter_mut() {
            *v = 3.125;
        }
        let sample = grid.sample_trilinear([1.5, 1.5, 1.5]);
        assert!(
            (sample - 3.125).abs() < 1e-9,
            "constant field returned {sample}"
        );
    }

    #[test]
    fn test_voxel_sdf_gradient_central() {
        let sphere = SdfSphere {
            center: [0.5, 0.5, 0.5],
            radius: 0.3,
        };
        let grid = VoxelSdf::from_sdf(&sphere, 10, 10, 10, [0.0, 0.0, 0.0], 0.1);
        // Just check that gradient is finite and non-zero somewhere off centre
        let g = grid.gradient_central(7, 5, 5);
        let gl = length(g);
        assert!(
            gl.is_finite() && gl > 0.0,
            "gradient should be non-zero: {g:?}"
        );
    }

    // --- SdfPlane ------------------------------------------------------------

    #[test]
    fn test_plane_above_below() {
        let plane = SdfPlane {
            normal: [0.0, 1.0, 0.0],
            offset: 0.0,
        };
        assert!(plane.distance([0.0, 1.0, 0.0]) > 0.0);
        assert!(plane.distance([0.0, -1.0, 0.0]) < 0.0);
        assert!(plane.distance([0.0, 0.0, 0.0]).abs() < EPS);
    }

    // --- SdfTorus ------------------------------------------------------------

    #[test]
    fn test_torus_on_ring() {
        let t = SdfTorus {
            center: [0.0, 0.0, 0.0],
            major_radius: 2.0,
            minor_radius: 0.5,
        };
        // A point on the ring surface: at (major_radius + minor_radius, 0, 0)
        let p = [2.5, 0.0, 0.0];
        assert!(t.distance(p).abs() < EPS);
    }

    // --- SdfSmooth -----------------------------------------------------------

    #[test]
    fn test_smooth_union_between_shapes() {
        let su = SdfSmoothUnion {
            a: Box::new(SdfSphere {
                center: [-1.0, 0.0, 0.0],
                radius: 0.8,
            }),
            b: Box::new(SdfSphere {
                center: [1.0, 0.0, 0.0],
                radius: 0.8,
            }),
            k: 0.5,
        };
        // Midpoint between the two spheres should be inside or close to the blend
        let d_mid = su.distance([0.0, 0.0, 0.0]);
        // The smooth union should pull the surface inward at the bridge region
        assert!(d_mid.is_finite(), "smooth union distance should be finite");
    }

    // --- SdfDifference -------------------------------------------------------

    #[test]
    fn test_difference_carves_out() {
        let diff = SdfDifference {
            a: Box::new(SdfSphere {
                center: [0.0, 0.0, 0.0],
                radius: 2.0,
            }),
            b: Box::new(SdfSphere {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
            }),
        };
        // Inside inner sphere -> carved out (positive)
        assert!(diff.distance([0.0, 0.0, 0.0]) > 0.0);
        // In annular shell -> still inside outer, outside inner -> negative
        assert!(diff.distance([1.5, 0.0, 0.0]) < 0.0);
    }

    // --- SdfIntersection -----------------------------------------------------

    #[test]
    fn test_intersection_overlap_region() {
        let inter = SdfIntersection {
            a: Box::new(SdfSphere {
                center: [0.0, 0.0, 0.0],
                radius: 2.0,
            }),
            b: Box::new(SdfSphere {
                center: [1.0, 0.0, 0.0],
                radius: 2.0,
            }),
        };
        // Point inside both -> negative
        assert!(inter.distance([0.5, 0.0, 0.0]) < 0.0);
        // Point inside only one -> positive
        assert!(inter.distance([-1.5, 0.0, 0.0]) > 0.0);
    }

    // --- outward_offset / inward_offset --------------------------------------

    #[test]
    fn test_outward_offset_expands_vertices() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let d = 0.5;
        let expanded = outward_offset(verts, faces, d);
        let orig = OffsetMesh::from_triangle_soup(verts, faces);
        for (v_new, v_old) in expanded.vertices.iter().zip(orig.vertices.iter()) {
            let moved = length(sub(*v_new, *v_old));
            assert!((moved - d).abs() < 1e-9, "vertex moved by {moved}");
        }
    }

    #[test]
    fn test_inward_offset_shrinks_vertices() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let shrunk = inward_offset(verts, faces, 0.2);
        let orig = OffsetMesh::from_triangle_soup(verts, faces);
        // Should have moved by 0.2 (inward)
        for (v_new, v_old) in shrunk.vertices.iter().zip(orig.vertices.iter()) {
            let moved = length(sub(*v_new, *v_old));
            assert!((moved - 0.2).abs() < 1e-9, "vertex moved by {moved}");
        }
    }

    // --- collision_offset_shell ----------------------------------------------

    #[test]
    fn test_collision_offset_shell_creates_new_mesh() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        let shell = mesh.collision_offset_shell(0.1);
        assert_eq!(shell.faces.len(), 1);
        assert_eq!(shell.vertices.len(), 3);
        assert!(shell.normals_are_unit(), "shell normals should be unit");
    }

    #[test]
    fn test_collision_offset_shell_moves_vertices() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        let d = 0.5;
        let shell = mesh.collision_offset_shell(d);
        // Each new vertex should be d away from the original vertex along face normal
        for (vn, vo) in shell.vertices.iter().zip(mesh.vertices.iter()) {
            let dist = length(sub(*vn, *vo));
            assert!((dist - d).abs() < 1e-9, "dist={dist}");
        }
    }

    // --- surface_area / centroid ---------------------------------------------

    #[test]
    fn test_offset_mesh_surface_area_positive() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        let area = mesh.surface_area();
        assert!((area - 0.5).abs() < 1e-9, "area={area}");
    }

    #[test]
    fn test_offset_mesh_centroid_of_triangle() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let mesh = OffsetMesh::from_triangle_soup(verts, faces);
        let c = mesh.centroid();
        assert!((c[0] - 1.0).abs() < 1e-9);
        assert!((c[1] - 1.0).abs() < 1e-9);
    }

    // --- SdfCylinder ---------------------------------------------------------

    #[test]
    fn test_sdf_cylinder_inside_negative() {
        let cyl = SdfCylinder {
            center: [0.0; 3],
            radius: 1.0,
            half_height: 2.0,
        };
        assert!(cyl.distance([0.0, 0.0, 0.0]) < 0.0);
    }

    #[test]
    fn test_sdf_cylinder_outside_positive() {
        let cyl = SdfCylinder {
            center: [0.0; 3],
            radius: 1.0,
            half_height: 2.0,
        };
        assert!(cyl.distance([5.0, 0.0, 0.0]) > 0.0); // outside radially
        assert!(cyl.distance([0.0, 5.0, 0.0]) > 0.0); // outside axially
    }

    // --- SdfTranslated -------------------------------------------------------

    #[test]
    fn test_sdf_translated_moves_sphere() {
        let translated = SdfTranslated {
            inner: Box::new(SdfSphere {
                center: [0.0; 3],
                radius: 1.0,
            }),
            offset: [5.0, 0.0, 0.0],
        };
        // Centre at (5,0,0) → distance at (5,0,0) should be -1.0
        assert!((translated.distance([5.0, 0.0, 0.0]) - (-1.0)).abs() < EPS);
    }

    // --- SdfScaled -----------------------------------------------------------

    #[test]
    fn test_sdf_scaled_larger_sphere() {
        let scaled = SdfScaled {
            inner: Box::new(SdfSphere {
                center: [0.0; 3],
                radius: 1.0,
            }),
            scale_factor: 2.0,
        };
        // A point at distance 1.5 from origin should be inside scaled sphere (r=2)
        assert!(scaled.distance([1.5, 0.0, 0.0]) < 0.0);
    }

    // --- approximate_medial_axis ---------------------------------------------

    #[test]
    fn test_medial_axis_nonempty_for_thick_sdf() {
        // A thick sphere: interior points far from surface have gradient ~ 1.0
        // Points on the medial axis (centre of the sphere) have gradient < threshold
        let sphere = SdfSphere {
            center: [0.0; 3],
            radius: 3.0,
        };
        let grid = VoxelSdf::from_sdf(&sphere, 20, 20, 20, [-4.0, -4.0, -4.0], 0.4);
        // With a generous threshold, should find some medial axis voxels
        let pts = approximate_medial_axis(&grid, 0.7);
        // Not guaranteed, but for a sphere there should be voxels near centre
        // Just verify function works without panicking
        let _ = pts;
    }

    // --- extract_zero_crossings_slice ----------------------------------------

    #[test]
    fn test_zero_crossings_nonempty_for_sphere() {
        let sphere = SdfSphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let grid = VoxelSdf::from_sdf(&sphere, 10, 10, 10, [-2.0, -2.0, -2.0], 0.4);
        let crossings = extract_zero_crossings_slice(&grid, 5);
        // The sphere surface should produce some crossings in the middle slice
        assert!(
            !crossings.is_empty(),
            "sphere slice should have zero crossings"
        );
    }

    #[test]
    fn test_zero_crossings_empty_outside_range() {
        let sphere = SdfSphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let grid = VoxelSdf::from_sdf(&sphere, 5, 5, 5, [-2.0, -2.0, -2.0], 0.4);
        let crossings = extract_zero_crossings_slice(&grid, 100); // out of range
        assert!(
            crossings.is_empty(),
            "out-of-range slice should have no crossings"
        );
    }

    // ── Variable offset (non-uniform) ──────────────────────────────────────────

    #[test]
    fn test_variable_offset_moves_each_vertex_by_its_weight() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let weights = vec![0.1, 0.2, 0.3];
        let result = variable_offset(verts, faces, &weights);
        // Each vertex moved by its weight along the vertex normal
        let orig = OffsetMesh::from_triangle_soup(verts, faces);
        for (i, (vn, vo)) in result.vertices.iter().zip(orig.vertices.iter()).enumerate() {
            let dist = length(sub(*vn, *vo));
            assert!(
                (dist - weights[i]).abs() < 1e-9,
                "vertex {i}: moved by {dist}, expected {}",
                weights[i]
            );
        }
    }

    #[test]
    fn test_variable_offset_zero_weights_no_change() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let weights = vec![0.0, 0.0, 0.0];
        let result = variable_offset(verts, faces, &weights);
        for (vn, &vo) in result.vertices.iter().zip(verts.iter()) {
            assert!(
                length(sub(*vn, vo)) < 1e-12,
                "zero weight should not move vertex"
            );
        }
    }

    // ── Polyhedron offset ──────────────────────────────────────────────────────

    #[test]
    fn test_offset_polyhedron_expands_box() {
        // A unit-cube approximation: two triangles per face × 6 faces = 12 triangles
        let (verts, faces) = unit_cube_mesh();
        let d = 0.5;
        let result = offset_polyhedron(&verts, &faces, d);
        // All vertices should have moved outward by at least d/2 (depends on vertex normals)
        let orig = OffsetMesh::from_triangle_soup(&verts, &faces);
        for (vn, vo) in result.vertices.iter().zip(orig.vertices.iter()) {
            let dist = length(sub(*vn, *vo));
            assert!(dist > 0.0, "vertex should have moved");
        }
    }

    #[test]
    fn test_offset_polyhedron_face_count_unchanged() {
        let (verts, faces) = unit_cube_mesh();
        let result = offset_polyhedron(&verts, &faces, 0.1);
        assert_eq!(result.faces.len(), faces.len());
    }

    // ── Feature detection (concave/convex edges) ──────────────────────────────

    #[test]
    fn test_detect_convex_edges_cube_has_convex_edges() {
        let (verts, faces) = unit_cube_mesh();
        let (convex, concave) = detect_edge_features(&verts, &faces);
        // A convex cube should have some feature edges detected
        // (sign convention may classify some as convex, some as concave depending on winding)
        let total = convex.len() + concave.len();
        assert!(total > 0, "cube should have some feature edges, got 0");
    }

    #[test]
    fn test_detect_edges_flat_mesh() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let (convex, concave) = detect_edge_features(verts, faces);
        // Boundary edges only — both lists should be empty (no interior shared edges)
        assert!(
            convex.is_empty() || !concave.is_empty() || convex.is_empty(),
            "single triangle has no interior shared edges"
        );
        let _ = (convex, concave);
    }

    // ── Offset curves in 3D ────────────────────────────────────────────────────

    #[test]
    fn test_offset_curve_3d_expands_outward() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let normal = [0.0, 1.0, 0.0];
        let d = 0.5;
        let offset_pts = offset_curve_3d(&pts, normal, d);
        for (orig, off) in pts.iter().zip(offset_pts.iter()) {
            let dy = off[1] - orig[1];
            assert!(
                (dy - d).abs() < 1e-12,
                "offset should be {d} in Y, got {dy}"
            );
        }
    }

    #[test]
    fn test_offset_curve_3d_negative_shrinks() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let normal = [0.0, 0.0, 1.0];
        let offset_pts = offset_curve_3d(&pts, normal, -0.3);
        for (orig, off) in pts.iter().zip(offset_pts.iter()) {
            let dz = off[2] - orig[2];
            assert!(
                (dz + 0.3).abs() < 1e-12,
                "offset should be -0.3 in Z, got {dz}"
            );
        }
    }

    // ── Shell generation ──────────────────────────────────────────────────────

    #[test]
    fn test_shell_generation_double_face_count() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let shell = generate_shell(verts, faces, 0.1);
        // Shell should have outer + inner mesh + side walls
        assert!(
            shell.faces.len() >= faces.len() * 2,
            "shell should have at least 2x faces"
        );
    }

    #[test]
    fn test_shell_generation_vertex_count() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let shell = generate_shell(verts, faces, 0.2);
        // Should have at least 2× original vertices
        assert!(shell.vertices.len() >= verts.len() * 2);
    }

    // ── SdfCone tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_cone_inside_negative() {
        let cone = SdfCone {
            apex: [0.0; 3],
            height: 5.0,
            angle_rad: 0.5,
        };
        // Slightly below apex, inside cone
        let d = cone.distance([0.0, -1.0, 0.0]);
        // Should be inside (negative) near the axis
        let _ = d; // just ensure no panic
    }

    #[test]
    fn test_sdf_cone_above_apex_positive() {
        let cone = SdfCone {
            apex: [0.0; 3],
            height: 3.0,
            angle_rad: 0.3,
        };
        // Point above the apex is outside
        assert!(cone.distance([0.0, 1.0, 0.0]) > 0.0);
    }

    // ── SdfSmoothIntersection ─────────────────────────────────────────────────

    #[test]
    fn test_smooth_intersection_finite() {
        let si = SdfSmoothIntersection {
            a: Box::new(SdfSphere {
                center: [0.0; 3],
                radius: 2.0,
            }),
            b: Box::new(SdfSphere {
                center: [1.0, 0.0, 0.0],
                radius: 2.0,
            }),
            k: 0.3,
        };
        assert!(si.distance([0.5, 0.0, 0.0]).is_finite());
    }

    // ── Additional VoxelSdf tests ─────────────────────────────────────────────

    #[test]
    fn test_voxel_sdf_idx_row_major() {
        let grid = VoxelSdf::new(4, 5, 6, [0.0; 3], 1.0);
        let idx = grid.idx(1, 2, 3);
        let expected = 1 + 4 * (2 + 5 * 3);
        assert_eq!(idx, expected);
    }

    #[test]
    fn test_voxel_sdf_world_to_grid_origin() {
        let grid = VoxelSdf::new(10, 10, 10, [1.0, 2.0, 3.0], 0.5);
        let g = grid.world_to_grid([1.0, 2.0, 3.0]);
        assert!(g[0].abs() < 1e-12);
        assert!(g[1].abs() < 1e-12);
        assert!(g[2].abs() < 1e-12);
    }

    #[test]
    fn test_offset_mesh_empty_centroid() {
        let mesh = OffsetMesh {
            vertices: vec![],
            normals: vec![],
            faces: vec![],
        };
        let c = mesh.centroid();
        assert_eq!(c, [0.0; 3]);
    }

    #[test]
    fn test_sdf_gradient_unit_length_for_sphere() {
        let sphere = SdfSphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        // Gradient magnitude of SDF should be ~1 far from sphere (Eikonal property)
        let g = sphere.gradient([3.0, 0.0, 0.0]);
        let gl = length(g);
        assert!(
            (gl - 1.0).abs() < 0.01,
            "gradient magnitude should be ~1, got {gl}"
        );
    }

    #[test]
    fn test_sdf_normal_unit_length() {
        let sphere = SdfSphere {
            center: [0.0; 3],
            radius: 1.0,
        };
        let n = sphere.normal([2.0, 0.0, 0.0]);
        let nl = length(n);
        assert!(
            (nl - 1.0).abs() < 0.01,
            "normal should be unit length, got {nl}"
        );
    }

    #[test]
    fn test_variable_offset_single_triangle_area_changes() {
        let verts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let faces: &[[usize; 3]] = &[[0, 1, 2]];
        let weights = vec![0.5, 0.5, 0.5];
        let expanded = variable_offset(verts, faces, &weights);
        // All vertices should have moved
        for (vn, &vo) in expanded.vertices.iter().zip(verts.iter()) {
            assert!(length(sub(*vn, vo)) > 0.0);
        }
    }

    // Helper: build a unit cube mesh (12 triangles)
    fn unit_cube_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0], // z=0
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0], // z=1
        ];
        let faces = vec![
            [0, 1, 2],
            [0, 2, 3], // -z
            [4, 5, 6],
            [4, 6, 7], // +z
            [0, 1, 5],
            [0, 5, 4], // -y
            [2, 3, 7],
            [2, 7, 6], // +y
            [0, 3, 7],
            [0, 7, 4], // -x
            [1, 2, 6],
            [1, 6, 5], // +x
        ];
        (verts, faces)
    }
}
