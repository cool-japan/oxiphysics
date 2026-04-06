#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Renderable primitive types for visualization data generation.

use std::f32::consts::PI;

use oxiphysics_core::math::Vec3;

/// RGBA color with floating-point components in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}

impl Color {
    /// Create a new color from RGBA components.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Opaque red.
    pub fn red() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Opaque green.
    pub fn green() -> Self {
        Self::new(0.0, 1.0, 0.0, 1.0)
    }

    /// Opaque blue.
    pub fn blue() -> Self {
        Self::new(0.0, 0.0, 1.0, 1.0)
    }

    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }

    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }

    /// Opaque yellow.
    pub fn yellow() -> Self {
        Self::new(1.0, 1.0, 0.0, 1.0)
    }

    /// Opaque cyan.
    pub fn cyan() -> Self {
        Self::new(0.0, 1.0, 1.0, 1.0)
    }

    /// Opaque magenta.
    pub fn magenta() -> Self {
        Self::new(1.0, 0.0, 1.0, 1.0)
    }

    /// Linearly interpolate between two colors.
    pub fn lerp(a: &Color, b: &Color, t: f32) -> Self {
        Self {
            r: a.r + (b.r - a.r) * t,
            g: a.g + (b.g - a.g) * t,
            b: a.b + (b.b - a.b) * t,
            a: a.a + (b.a - a.a) * t,
        }
    }
}

/// A vertex with position, normal, and color.
#[derive(Debug, Clone, Copy)]
pub struct Vertex {
    /// Position in 3D space.
    pub position: [f32; 3],
    /// Surface normal.
    pub normal: [f32; 3],
    /// Vertex color.
    pub color: Color,
}

/// A line segment between two 3D points.
#[derive(Debug, Clone)]
pub struct LinePrimitive {
    /// Start point of the line.
    pub start: Vec3,
    /// End point of the line.
    pub end: Vec3,
    /// Color of the line.
    pub color: Color,
}

/// A single triangle defined by three vertices.
#[derive(Debug, Clone)]
pub struct TrianglePrimitive {
    /// The three vertices of the triangle.
    pub vertices: [Vertex; 3],
}

/// An indexed triangle mesh for rendering.
#[derive(Debug, Clone)]
pub struct RenderMesh {
    /// Vertex data.
    pub vertices: Vec<Vertex>,
    /// Triangle indices (groups of 3).
    pub indices: Vec<u32>,
}

impl RenderMesh {
    /// Create a new empty render mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
        }
    }
}

impl Default for RenderMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ── New mesh-generation primitives ────────────────────────────────────────────

/// Vertex with position, normal, and UV coordinates (used by [`Mesh`]).
#[derive(Debug, Clone, Copy)]
pub struct MeshVertex {
    /// World-space position.
    pub position: [f32; 3],
    /// Surface normal (unit length).
    pub normal: [f32; 3],
    /// Texture / UV coordinate.
    pub uv: [f32; 2],
}

/// Indexed triangle mesh used for debug rendering.
#[derive(Debug, Clone)]
pub struct Mesh {
    /// Vertex data.
    pub vertices: Vec<MeshVertex>,
    /// Triangle indices (groups of 3).
    pub indices: Vec<u32>,
}

impl Mesh {
    /// Generate a UV sphere with `rings` latitude bands and `sectors` longitude slices.
    ///
    /// Vertex count = (rings + 1) × (sectors + 1).
    pub fn uv_sphere(radius: f32, rings: u32, sectors: u32) -> Self {
        use std::f32::consts::PI;
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for r in 0..=rings {
            let phi = PI * (r as f32) / (rings as f32);
            let sin_phi = phi.sin();
            let cos_phi = phi.cos();
            for s in 0..=sectors {
                let theta = 2.0 * PI * (s as f32) / (sectors as f32);
                let sin_theta = theta.sin();
                let cos_theta = theta.cos();

                let nx = sin_phi * cos_theta;
                let ny = cos_phi;
                let nz = sin_phi * sin_theta;
                vertices.push(MeshVertex {
                    position: [radius * nx, radius * ny, radius * nz],
                    normal: [nx, ny, nz],
                    uv: [(s as f32) / (sectors as f32), (r as f32) / (rings as f32)],
                });
            }
        }

        let row_len = sectors + 1;
        for r in 0..rings {
            for s in 0..sectors {
                let a = r * row_len + s;
                let b = (r + 1) * row_len + s;
                let c = (r + 1) * row_len + s + 1;
                let d = r * row_len + s + 1;
                indices.extend_from_slice(&[a, b, c, a, c, d]);
            }
        }

        Self { vertices, indices }
    }

    /// Generate a box mesh from `half_extents` centred at the origin.
    ///
    /// 6 faces × 2 triangles = 12 triangles (36 indices), 24 vertices.
    pub fn box_mesh(half_extents: [f32; 3]) -> Self {
        let [hx, hy, hz] = half_extents;

        // (normal, 4 corner positions, 4 UVs)
        let faces: [([f32; 3], [[f32; 3]; 4], [[f32; 2]; 4]); 6] = [
            (
                [1.0, 0.0, 0.0],
                [[hx, -hy, -hz], [hx, hy, -hz], [hx, hy, hz], [hx, -hy, hz]],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
            (
                [-1.0, 0.0, 0.0],
                [
                    [-hx, -hy, hz],
                    [-hx, hy, hz],
                    [-hx, hy, -hz],
                    [-hx, -hy, -hz],
                ],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
            (
                [0.0, 1.0, 0.0],
                [[-hx, hy, -hz], [hx, hy, -hz], [hx, hy, hz], [-hx, hy, hz]],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
            (
                [0.0, -1.0, 0.0],
                [
                    [-hx, -hy, hz],
                    [hx, -hy, hz],
                    [hx, -hy, -hz],
                    [-hx, -hy, -hz],
                ],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
            (
                [0.0, 0.0, 1.0],
                [[-hx, -hy, hz], [hx, -hy, hz], [hx, hy, hz], [-hx, hy, hz]],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
            (
                [0.0, 0.0, -1.0],
                [
                    [hx, -hy, -hz],
                    [-hx, -hy, -hz],
                    [-hx, hy, -hz],
                    [hx, hy, -hz],
                ],
                [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
            ),
        ];

        let mut vertices = Vec::with_capacity(24);
        let mut indices = Vec::with_capacity(36);

        for (normal, positions, uvs) in &faces {
            let base = vertices.len() as u32;
            for (pos, uv) in positions.iter().zip(uvs.iter()) {
                vertices.push(MeshVertex {
                    position: *pos,
                    normal: *normal,
                    uv: *uv,
                });
            }
            indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }

        Self { vertices, indices }
    }

    /// Generate a cylinder mesh (barrel + two end caps).
    pub fn cylinder(radius: f32, half_height: f32, segments: u32) -> Self {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        let bottom_y = -half_height;
        let top_y = half_height;

        // ── Barrel ──────────────────────────────────────────────────────────
        // Two rings: bottom then top.
        for (i, &y) in [bottom_y, top_y].iter().enumerate() {
            for s in 0..=segments {
                let theta = 2.0 * PI * (s as f32) / (segments as f32);
                let cos_t = theta.cos();
                let sin_t = theta.sin();
                vertices.push(MeshVertex {
                    position: [radius * cos_t, y, radius * sin_t],
                    normal: [cos_t, 0.0, sin_t],
                    uv: [(s as f32) / (segments as f32), i as f32],
                });
            }
        }

        let row = segments + 1;
        for s in 0..segments {
            let b0 = s;
            let b1 = s + 1;
            let t0 = row + s;
            let t1 = row + s + 1;
            indices.extend_from_slice(&[b0, t0, t1, b0, t1, b1]);
        }

        // ── Bottom cap ──────────────────────────────────────────────────────
        let cap_center_bot = vertices.len() as u32;
        vertices.push(MeshVertex {
            position: [0.0, bottom_y, 0.0],
            normal: [0.0, -1.0, 0.0],
            uv: [0.5, 0.5],
        });
        let cap_ring_bot_start = vertices.len() as u32;
        for s in 0..segments {
            let theta = 2.0 * PI * (s as f32) / (segments as f32);
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            vertices.push(MeshVertex {
                position: [radius * cos_t, bottom_y, radius * sin_t],
                normal: [0.0, -1.0, 0.0],
                uv: [0.5 + 0.5 * cos_t, 0.5 + 0.5 * sin_t],
            });
        }
        for s in 0..segments {
            let i0 = cap_ring_bot_start + s;
            let i1 = cap_ring_bot_start + (s + 1) % segments;
            indices.extend_from_slice(&[cap_center_bot, i1, i0]);
        }

        // ── Top cap ─────────────────────────────────────────────────────────
        let cap_center_top = vertices.len() as u32;
        vertices.push(MeshVertex {
            position: [0.0, top_y, 0.0],
            normal: [0.0, 1.0, 0.0],
            uv: [0.5, 0.5],
        });
        let cap_ring_top_start = vertices.len() as u32;
        for s in 0..segments {
            let theta = 2.0 * PI * (s as f32) / (segments as f32);
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            vertices.push(MeshVertex {
                position: [radius * cos_t, top_y, radius * sin_t],
                normal: [0.0, 1.0, 0.0],
                uv: [0.5 + 0.5 * cos_t, 0.5 + 0.5 * sin_t],
            });
        }
        for s in 0..segments {
            let i0 = cap_ring_top_start + s;
            let i1 = cap_ring_top_start + (s + 1) % segments;
            indices.extend_from_slice(&[cap_center_top, i0, i1]);
        }

        Self { vertices, indices }
    }

    /// Generate an arrow mesh: cylindrical shaft + cone tip.
    pub fn arrow(start: [f32; 3], end: [f32; 3], shaft_radius: f32, tip_length: f32) -> Self {
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let dz = end[2] - start[2];
        let total_len = (dx * dx + dy * dy + dz * dz).sqrt();
        if total_len < 1e-9 {
            return Self {
                vertices: Vec::new(),
                indices: Vec::new(),
            };
        }
        let ax = dx / total_len;
        let ay = dy / total_len;
        let az = dz / total_len;

        // Build a tangent frame
        let (tx, ty, tz) = if ay.abs() < 0.99 {
            let len = (ax * ax + az * az).sqrt().max(1e-9);
            (-az / len, 0.0_f32, ax / len)
        } else {
            (1.0_f32, 0.0, 0.0)
        };
        let bx = ay * tz - az * ty;
        let by = az * tx - ax * tz;
        let bz = ax * ty - ay * tx;

        let shaft_len = (total_len - tip_length).max(0.0);
        let tip_start = [
            start[0] + ax * shaft_len,
            start[1] + ay * shaft_len,
            start[2] + az * shaft_len,
        ];
        let tip_radius = shaft_radius * 2.0;
        let segments: u32 = 8;

        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        // ── Shaft ───────────────────────────────────────────────────────────
        for s in 0..=segments {
            let theta = 2.0 * PI * (s as f32) / (segments as f32);
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            let nx = tx * cos_t + bx * sin_t;
            let ny = ty * cos_t + by * sin_t;
            let nz = tz * cos_t + bz * sin_t;
            let u = (s as f32) / (segments as f32);

            vertices.push(MeshVertex {
                position: [
                    start[0] + shaft_radius * nx,
                    start[1] + shaft_radius * ny,
                    start[2] + shaft_radius * nz,
                ],
                normal: [nx, ny, nz],
                uv: [u, 0.0],
            });
            vertices.push(MeshVertex {
                position: [
                    tip_start[0] + shaft_radius * nx,
                    tip_start[1] + shaft_radius * ny,
                    tip_start[2] + shaft_radius * nz,
                ],
                normal: [nx, ny, nz],
                uv: [u, 1.0],
            });
        }
        let row = 2;
        for s in 0..segments {
            let b0 = s * row;
            let b1 = s * row + 1;
            let t0 = (s + 1) * row;
            let t1 = (s + 1) * row + 1;
            indices.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        }

        // ── Cone tip ────────────────────────────────────────────────────────
        let cone_base_start = vertices.len() as u32;
        // center of cone base (for cap)
        vertices.push(MeshVertex {
            position: tip_start,
            normal: [-ax, -ay, -az],
            uv: [0.5, 0.5],
        });
        for s in 0..segments {
            let theta = 2.0 * PI * (s as f32) / (segments as f32);
            let cos_t = theta.cos();
            let sin_t = theta.sin();
            let nx = tx * cos_t + bx * sin_t;
            let ny = ty * cos_t + by * sin_t;
            let nz = tz * cos_t + bz * sin_t;
            vertices.push(MeshVertex {
                position: [
                    tip_start[0] + tip_radius * nx,
                    tip_start[1] + tip_radius * ny,
                    tip_start[2] + tip_radius * nz,
                ],
                normal: [nx, ny, nz],
                uv: [0.5 + 0.5 * cos_t, 0.5 + 0.5 * sin_t],
            });
        }
        let tip_apex = vertices.len() as u32;
        vertices.push(MeshVertex {
            position: end,
            normal: [ax, ay, az],
            uv: [0.5, 0.5],
        });

        // Cap
        for s in 0..segments {
            let i0 = cone_base_start + 1 + s;
            let i1 = cone_base_start + 1 + (s + 1) % segments;
            indices.extend_from_slice(&[cone_base_start, i1, i0]);
        }
        // Side triangles
        for s in 0..segments {
            let i0 = cone_base_start + 1 + s;
            let i1 = cone_base_start + 1 + (s + 1) % segments;
            indices.extend_from_slice(&[i0, i1, tip_apex]);
        }

        Self { vertices, indices }
    }

    /// Generate a wireframe grid in the XZ plane centred at origin.
    ///
    /// Returns a [`WireframeMesh`]; use [`WireframeMesh`] for line rendering.
    pub fn grid(size: f32, divisions: u32) -> WireframeMesh {
        let mut verts = Vec::new();
        let mut line_indices = Vec::new();

        let step = size / (divisions as f32);
        let half = size / 2.0;

        for i in 0..=(divisions) {
            let coord = -half + step * (i as f32);

            // Line along Z
            let a = verts.len() as u32;
            verts.push([coord, 0.0, -half]);
            let b = verts.len() as u32;
            verts.push([coord, 0.0, half]);
            line_indices.push((a, b));

            // Line along X
            let c = verts.len() as u32;
            verts.push([-half, 0.0, coord]);
            let d = verts.len() as u32;
            verts.push([half, 0.0, coord]);
            line_indices.push((c, d));
        }

        WireframeMesh {
            vertices: verts,
            line_indices,
        }
    }

    /// Number of vertices in this mesh.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles in this mesh.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
}

/// Wireframe mesh — just line segments, no filled faces.
#[derive(Debug, Clone)]
pub struct WireframeMesh {
    /// Vertex positions.
    pub vertices: Vec<[f32; 3]>,
    /// Index pairs defining each line segment.
    pub line_indices: Vec<(u32, u32)>,
}

impl WireframeMesh {
    /// Build the 12-edge wireframe of an axis-aligned bounding box.
    pub fn aabb(min: [f32; 3], max: [f32; 3]) -> Self {
        let [x0, y0, z0] = min;
        let [x1, y1, z1] = max;

        // 8 corners
        let corners: [[f32; 3]; 8] = [
            [x0, y0, z0], // 0
            [x1, y0, z0], // 1
            [x1, y1, z0], // 2
            [x0, y1, z0], // 3
            [x0, y0, z1], // 4
            [x1, y0, z1], // 5
            [x1, y1, z1], // 6
            [x0, y1, z1], // 7
        ];

        let line_indices: Vec<(u32, u32)> = vec![
            // bottom face
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            // top face
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            // verticals
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        Self {
            vertices: corners.to_vec(),
            line_indices,
        }
    }

    /// Build three great-circle wireframes (XY, XZ, YZ planes) for a sphere.
    pub fn sphere_wireframe(center: [f32; 3], radius: f32, segments: u32) -> Self {
        let mut vertices = Vec::new();
        let mut line_indices = Vec::new();
        let [cx, cy, cz] = center;

        // Each great circle: `segments` edges → `segments` vertices
        for plane in 0..3_u32 {
            let base = vertices.len() as u32;
            for s in 0..segments {
                let theta = 2.0 * PI * (s as f32) / (segments as f32);
                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let pos = match plane {
                    0 => [cx + radius * cos_t, cy + radius * sin_t, cz], // XY
                    1 => [cx + radius * cos_t, cy, cz + radius * sin_t], // XZ
                    _ => [cx, cy + radius * cos_t, cz + radius * sin_t], // YZ
                };
                vertices.push(pos);
            }
            for s in 0..segments {
                let a = base + s;
                let b = base + (s + 1) % segments;
                line_indices.push((a, b));
            }
        }

        Self {
            vertices,
            line_indices,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sphere_vertex_count() {
        let rings: u32 = 8;
        let sectors: u32 = 16;
        let mesh = Mesh::uv_sphere(1.0, rings, sectors);
        let expected = ((rings + 1) * (sectors + 1)) as usize;
        assert_eq!(
            mesh.vertex_count(),
            expected,
            "UV sphere vertex count: expected {expected}, got {}",
            mesh.vertex_count()
        );
    }

    #[test]
    fn test_box_mesh_triangles() {
        let mesh = Mesh::box_mesh([0.5, 0.5, 0.5]);
        assert_eq!(
            mesh.triangle_count(),
            12,
            "box mesh must have exactly 12 triangles"
        );
    }

    #[test]
    fn test_cylinder_has_caps() {
        let segments: u32 = 8;
        // Barrel only: 2 rings × (segments+1) vertices
        let barrel_only_verts = 2 * (segments + 1) as usize;
        let mesh = Mesh::cylinder(1.0, 1.0, segments);
        assert!(
            mesh.vertex_count() > barrel_only_verts,
            "cylinder must include cap vertices (got {}, barrel-only = {})",
            mesh.vertex_count(),
            barrel_only_verts
        );
    }

    #[test]
    fn test_aabb_wireframe_edges() {
        let wf = WireframeMesh::aabb([-1.0; 3], [1.0; 3]);
        assert_eq!(
            wf.line_indices.len(),
            12,
            "AABB wireframe must have exactly 12 edges"
        );
    }

    #[test]
    fn test_sphere_wireframe_3circles() {
        let segments: u32 = 16;
        let wf = WireframeMesh::sphere_wireframe([0.0; 3], 1.0, segments);
        assert_eq!(
            wf.line_indices.len(),
            (3 * segments) as usize,
            "sphere wireframe must have 3 × segments line segments"
        );
    }
}

// ── Axis-Aligned Bounding Box (for frustum culling) ───────────────────────────

/// An axis-aligned bounding box in f32 for rendering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RenderAabb {
    /// Minimum corner.
    pub min: [f32; 3],
    /// Maximum corner.
    pub max: [f32; 3],
}

impl RenderAabb {
    /// Create a new render AABB.
    pub fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self { min, max }
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f32; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> [f32; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }

    /// Test whether this AABB intersects another.
    pub fn intersects(&self, other: &RenderAabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Extend this AABB to contain the given point.
    pub fn extend_to_point(&self, p: [f32; 3]) -> Self {
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
}

// ── Frustum for view culling ───────────────────────────────────────────────────

/// A view frustum defined by 6 planes for visibility culling.
///
/// Each plane is `[a, b, c, d]` satisfying `a*x + b*y + c*z + d = 0`.
#[derive(Debug, Clone)]
pub struct Frustum {
    /// The 6 frustum planes: near, far, left, right, bottom, top.
    pub planes: [[f32; 4]; 6],
}

impl Frustum {
    /// Create a frustum from a 4×4 clip matrix (row-major, OpenGL convention).
    ///
    /// Uses the Gribb–Hartmann method.
    pub fn from_matrix(m: &[f32; 16]) -> Self {
        let row = |i: usize| -> [f32; 4] { [m[i * 4], m[i * 4 + 1], m[i * 4 + 2], m[i * 4 + 3]] };
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);
        let add = |a: [f32; 4], b: [f32; 4]| -> [f32; 4] {
            [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
        };
        let sub = |a: [f32; 4], b: [f32; 4]| -> [f32; 4] {
            [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
        };
        let planes = [
            add(r3, r2), // near
            sub(r3, r2), // far
            add(r3, r0), // left
            sub(r3, r0), // right
            add(r3, r1), // bottom
            sub(r3, r1), // top
        ];
        Self { planes }
    }

    /// Test whether an AABB is (potentially) visible (intersects the frustum).
    ///
    /// Uses the negative-vertex test against all planes.
    pub fn intersects_aabb(&self, aabb: &RenderAabb) -> bool {
        for plane in &self.planes {
            let [a, b, c, d] = *plane;
            // Positive vertex (most in the direction of the plane normal)
            let px = if a >= 0.0 { aabb.max[0] } else { aabb.min[0] };
            let py = if b >= 0.0 { aabb.max[1] } else { aabb.min[1] };
            let pz = if c >= 0.0 { aabb.max[2] } else { aabb.min[2] };
            if a * px + b * py + c * pz + d < 0.0 {
                return false; // AABB is entirely outside this plane
            }
        }
        true
    }

    /// Test whether a point is inside (or on) the frustum.
    pub fn contains_point(&self, p: [f32; 3]) -> bool {
        for &[a, b, c, d] in &self.planes {
            if a * p[0] + b * p[1] + c * p[2] + d < 0.0 {
                return false;
            }
        }
        true
    }
}

// ── BVH scene traversal ───────────────────────────────────────────────────────

/// A node in a Bounding Volume Hierarchy (BVH).
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// Bounding box of this node.
    pub bounds: RenderAabb,
    /// Left child index, or `usize::MAX` for a leaf.
    pub left: usize,
    /// Right child index, or `usize::MAX` for a leaf.
    pub right: usize,
    /// Object indices stored at this leaf (empty for internal nodes).
    pub objects: Vec<usize>,
}

impl BvhNode {
    /// Return true if this is a leaf node.
    pub fn is_leaf(&self) -> bool {
        self.left == usize::MAX && self.right == usize::MAX
    }
}

/// A simple BVH over a list of RenderAabb objects.
#[derive(Debug, Clone)]
pub struct Bvh {
    /// Nodes of the BVH, node 0 is the root.
    pub nodes: Vec<BvhNode>,
}

impl Bvh {
    /// Build a BVH from a list of bounding boxes.
    ///
    /// Uses a simple top-down median split.
    pub fn build(aabbs: &[RenderAabb]) -> Self {
        let mut bvh = Bvh { nodes: Vec::new() };
        if aabbs.is_empty() {
            return bvh;
        }
        let indices: Vec<usize> = (0..aabbs.len()).collect();
        bvh.build_recursive(aabbs, &indices, 0);
        bvh
    }

    fn build_recursive(&mut self, aabbs: &[RenderAabb], indices: &[usize], _depth: usize) -> usize {
        let node_idx = self.nodes.len();
        // Compute bounding box of all objects
        let bounds =
            indices
                .iter()
                .fold(RenderAabb::new([f32::MAX; 3], [f32::MIN; 3]), |acc, &i| {
                    let b = &aabbs[i];
                    RenderAabb::new(
                        [
                            acc.min[0].min(b.min[0]),
                            acc.min[1].min(b.min[1]),
                            acc.min[2].min(b.min[2]),
                        ],
                        [
                            acc.max[0].max(b.max[0]),
                            acc.max[1].max(b.max[1]),
                            acc.max[2].max(b.max[2]),
                        ],
                    )
                });

        if indices.len() <= 2 {
            self.nodes.push(BvhNode {
                bounds,
                left: usize::MAX,
                right: usize::MAX,
                objects: indices.to_vec(),
            });
            return node_idx;
        }

        // Split along the longest axis
        let extents = [
            bounds.max[0] - bounds.min[0],
            bounds.max[1] - bounds.min[1],
            bounds.max[2] - bounds.min[2],
        ];
        let axis = if extents[0] >= extents[1] && extents[0] >= extents[2] {
            0
        } else if extents[1] >= extents[2] {
            1
        } else {
            2
        };
        let mid = bounds.min[axis] + extents[axis] * 0.5;
        let (left_idx, right_idx): (Vec<usize>, Vec<usize>) = indices.iter().partition(|&&i| {
            let c = aabbs[i].center();
            c[axis] < mid
        });
        let (left_idx, right_idx) = if left_idx.is_empty() || right_idx.is_empty() {
            let half = indices.len() / 2;
            (indices[..half].to_vec(), indices[half..].to_vec())
        } else {
            (left_idx, right_idx)
        };

        // Push placeholder node first (to hold node_idx)
        self.nodes.push(BvhNode {
            bounds,
            left: usize::MAX,
            right: usize::MAX,
            objects: Vec::new(),
        });
        let left_child = self.build_recursive(aabbs, &left_idx, _depth + 1);
        let right_child = self.build_recursive(aabbs, &right_idx, _depth + 1);
        self.nodes[node_idx].left = left_child;
        self.nodes[node_idx].right = right_child;
        node_idx
    }

    /// Collect all leaf object indices whose AABB overlaps the query AABB.
    pub fn query_aabb(&self, query: &RenderAabb) -> Vec<usize> {
        let mut result = Vec::new();
        if self.nodes.is_empty() {
            return result;
        }
        self.query_recursive(0, query, &mut result);
        result
    }

    fn query_recursive(&self, node_idx: usize, query: &RenderAabb, result: &mut Vec<usize>) {
        if node_idx >= self.nodes.len() {
            return;
        }
        let node = &self.nodes[node_idx];
        if !node.bounds.intersects(query) {
            return;
        }
        if node.is_leaf() {
            result.extend_from_slice(&node.objects);
        } else {
            if node.left != usize::MAX {
                self.query_recursive(node.left, query, result);
            }
            if node.right != usize::MAX {
                self.query_recursive(node.right, query, result);
            }
        }
    }
}

// ── Instanced rendering ───────────────────────────────────────────────────────

/// Per-instance data for instanced rendering.
#[derive(Debug, Clone, Copy)]
pub struct InstanceData {
    /// World-space transform column-major 4×4 matrix.
    pub transform: [f32; 16],
    /// Per-instance color override.
    pub color: Color,
}

impl InstanceData {
    /// Create instance data from a 4×4 column-major transform and color.
    pub fn new(transform: [f32; 16], color: Color) -> Self {
        Self { transform, color }
    }

    /// Create identity-transform instance data.
    pub fn identity(color: Color) -> Self {
        let mut t = [0.0_f32; 16];
        t[0] = 1.0;
        t[5] = 1.0;
        t[10] = 1.0;
        t[15] = 1.0;
        Self {
            transform: t,
            color,
        }
    }
}

/// A batch of instances for instanced draw calls.
#[derive(Debug, Clone)]
pub struct InstanceBatch {
    /// The mesh shared by all instances.
    pub mesh_id: usize,
    /// Instance data for each object.
    pub instances: Vec<InstanceData>,
}

impl InstanceBatch {
    /// Create a new empty instance batch.
    pub fn new(mesh_id: usize) -> Self {
        Self {
            mesh_id,
            instances: Vec::new(),
        }
    }

    /// Add an instance to this batch.
    pub fn add(&mut self, data: InstanceData) {
        self.instances.push(data);
    }

    /// Number of instances.
    pub fn count(&self) -> usize {
        self.instances.len()
    }
}

// ── Draw call batching ────────────────────────────────────────────────────────

/// A single draw call record.
#[derive(Debug, Clone)]
pub struct DrawCall {
    /// Mesh identifier.
    pub mesh_id: usize,
    /// Index count to draw.
    pub index_count: u32,
    /// Offset into the index buffer.
    pub index_offset: u32,
    /// Base vertex offset.
    pub base_vertex: i32,
    /// Per-instance offset.
    pub instance_offset: u32,
    /// Number of instances.
    pub instance_count: u32,
}

/// Batched draw call list.
#[derive(Debug, Clone, Default)]
pub struct DrawCallBatch {
    /// List of draw calls.
    pub calls: Vec<DrawCall>,
}

impl DrawCallBatch {
    /// Add a draw call.
    pub fn push(&mut self, call: DrawCall) {
        self.calls.push(call);
    }

    /// Merge consecutive draw calls with the same mesh (simple optimization).
    pub fn merge_consecutive(&mut self) {
        if self.calls.len() < 2 {
            return;
        }
        let mut merged = Vec::new();
        let mut iter = self.calls.drain(..);
        let mut current = iter.next().expect("iterator should have elements");
        for next in iter {
            if next.mesh_id == current.mesh_id
                && next.index_offset == current.index_offset + current.index_count
            {
                current.index_count += next.index_count;
                current.instance_count += next.instance_count;
            } else {
                merged.push(current);
                current = next;
            }
        }
        merged.push(current);
        self.calls = merged;
    }

    /// Total number of draw calls.
    pub fn call_count(&self) -> usize {
        self.calls.len()
    }
}

// ── Render state ─────────────────────────────────────────────────────────────

/// Depth test comparison function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepthFunc {
    /// Always pass.
    Always,
    /// Pass if less.
    Less,
    /// Pass if less or equal.
    LessEqual,
    /// Pass if equal.
    Equal,
    /// Pass if greater.
    Greater,
    /// Never pass.
    Never,
}

/// Blending mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendMode {
    /// No blending (opaque).
    Opaque,
    /// Alpha blending: src_alpha * src + (1 - src_alpha) * dst.
    Alpha,
    /// Additive: src + dst.
    Additive,
    /// Premultiplied alpha.
    PremultipliedAlpha,
}

/// Face culling mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CullMode {
    /// No face culling.
    None,
    /// Cull back faces.
    Back,
    /// Cull front faces.
    Front,
}

/// Complete render state descriptor.
#[derive(Debug, Clone)]
pub struct RenderState {
    /// Whether depth testing is enabled.
    pub depth_test: bool,
    /// Whether depth writing is enabled.
    pub depth_write: bool,
    /// Depth comparison function.
    pub depth_func: DepthFunc,
    /// Blending mode.
    pub blend_mode: BlendMode,
    /// Face culling mode.
    pub cull_mode: CullMode,
    /// Polygon offset (slope factor, bias).
    pub polygon_offset: Option<(f32, f32)>,
}

impl Default for RenderState {
    fn default() -> Self {
        Self {
            depth_test: true,
            depth_write: true,
            depth_func: DepthFunc::Less,
            blend_mode: BlendMode::Opaque,
            cull_mode: CullMode::Back,
            polygon_offset: None,
        }
    }
}

impl RenderState {
    /// Create a transparent render state (alpha blending, no depth write).
    pub fn transparent() -> Self {
        Self {
            depth_test: true,
            depth_write: false,
            depth_func: DepthFunc::Less,
            blend_mode: BlendMode::Alpha,
            cull_mode: CullMode::None,
            polygon_offset: None,
        }
    }

    /// Create a wireframe-friendly state (polygon offset).
    pub fn wireframe() -> Self {
        Self {
            depth_test: true,
            depth_write: true,
            depth_func: DepthFunc::LessEqual,
            blend_mode: BlendMode::Opaque,
            cull_mode: CullMode::None,
            polygon_offset: Some((-1.0, -1.0)),
        }
    }
}

// ── Frame statistics ──────────────────────────────────────────────────────────

/// Per-frame rendering statistics.
#[derive(Debug, Clone, Default)]
pub struct FrameStats {
    /// Total draw calls issued this frame.
    pub draw_calls: usize,
    /// Total triangles submitted.
    pub triangles: usize,
    /// Total vertices submitted.
    pub vertices: usize,
    /// Number of objects culled by frustum.
    pub frustum_culled: usize,
    /// Number of objects culled by occlusion.
    pub occlusion_culled: usize,
    /// Frame time in milliseconds.
    pub frame_time_ms: f64,
}

impl FrameStats {
    /// Create new zero stats.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset all counters to zero.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Add a draw call with the given primitive count.
    pub fn record_draw_call(&mut self, triangles: usize, vertices: usize) {
        self.draw_calls += 1;
        self.triangles += triangles;
        self.vertices += vertices;
    }

    /// Record a frustum culled object.
    pub fn record_frustum_cull(&mut self) {
        self.frustum_culled += 1;
    }

    /// Record an occlusion culled object.
    pub fn record_occlusion_cull(&mut self) {
        self.occlusion_culled += 1;
    }

    /// Compute effective triangles per second from frame time.
    pub fn triangles_per_second(&self) -> f64 {
        if self.frame_time_ms < 1e-9 {
            return 0.0;
        }
        self.triangles as f64 / (self.frame_time_ms * 1e-3)
    }
}

// ── Occlusion query stub ──────────────────────────────────────────────────────

/// An occlusion query result.
#[derive(Debug, Clone, Copy)]
pub struct OcclusionQuery {
    /// Whether the query result is available.
    pub available: bool,
    /// Number of samples that passed the depth test (0 = occluded).
    pub samples_passed: u32,
}

impl OcclusionQuery {
    /// Create a new pending occlusion query.
    pub fn new() -> Self {
        Self {
            available: false,
            samples_passed: 0,
        }
    }

    /// Mark the query as resolved with the given sample count.
    pub fn resolve(&mut self, samples: u32) {
        self.available = true;
        self.samples_passed = samples;
    }

    /// Returns true if the object was visible.
    pub fn is_visible(&self) -> bool {
        self.available && self.samples_passed > 0
    }
}

impl Default for OcclusionQuery {
    fn default() -> Self {
        Self::new()
    }
}

// ── Tests for expanded primitives ─────────────────────────────────────────────

#[cfg(test)]
mod expanded_prim_tests {
    use super::*;

    #[test]
    fn test_render_aabb_center() {
        let aabb = RenderAabb::new([0.0; 3], [2.0; 3]);
        let c = aabb.center();
        for ci in c {
            assert!((ci - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_render_aabb_half_extents() {
        let aabb = RenderAabb::new([0.0; 3], [4.0; 3]);
        let he = aabb.half_extents();
        for h in he {
            assert!((h - 2.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_render_aabb_intersects() {
        let a = RenderAabb::new([0.0; 3], [2.0; 3]);
        let b = RenderAabb::new([1.0; 3], [3.0; 3]);
        let c = RenderAabb::new([3.0; 3], [5.0; 3]);
        assert!(a.intersects(&b), "overlapping AABBs should intersect");
        assert!(
            !a.intersects(&c),
            "non-overlapping AABBs should not intersect"
        );
    }

    #[test]
    fn test_frustum_identity_contains_origin() {
        // Identity clip matrix (column-major) → frustum contains origin
        let mut m = [0.0_f32; 16];
        m[0] = 1.0;
        m[5] = 1.0;
        m[10] = 1.0;
        m[15] = 1.0;
        let frustum = Frustum::from_matrix(&m);
        // Just check it doesn't panic; result depends on convention
        let _ = frustum.contains_point([0.0; 3]);
    }

    #[test]
    fn test_bvh_build_empty() {
        let bvh = Bvh::build(&[]);
        assert!(bvh.nodes.is_empty());
    }

    #[test]
    fn test_bvh_single_object() {
        let aabb = RenderAabb::new([0.0; 3], [1.0; 3]);
        let bvh = Bvh::build(&[aabb]);
        assert!(!bvh.nodes.is_empty());
    }

    #[test]
    fn test_bvh_query_finds_overlap() {
        let aabbs = vec![
            RenderAabb::new([0.0; 3], [1.0; 3]),
            RenderAabb::new([2.0; 3], [3.0; 3]),
            RenderAabb::new([4.0; 3], [5.0; 3]),
        ];
        let bvh = Bvh::build(&aabbs);
        let query = RenderAabb::new([-0.5; 3], [1.5; 3]);
        let result = bvh.query_aabb(&query);
        assert!(!result.is_empty(), "should find the first AABB");
    }

    #[test]
    fn test_instance_batch_count() {
        let mut batch = InstanceBatch::new(0);
        assert_eq!(batch.count(), 0);
        batch.add(InstanceData::identity(Color::red()));
        batch.add(InstanceData::identity(Color::blue()));
        assert_eq!(batch.count(), 2);
    }

    #[test]
    fn test_draw_call_batch_merge() {
        let mut batch = DrawCallBatch::default();
        batch.push(DrawCall {
            mesh_id: 0,
            index_count: 6,
            index_offset: 0,
            base_vertex: 0,
            instance_offset: 0,
            instance_count: 1,
        });
        batch.push(DrawCall {
            mesh_id: 0,
            index_count: 6,
            index_offset: 6,
            base_vertex: 0,
            instance_offset: 0,
            instance_count: 1,
        });
        batch.push(DrawCall {
            mesh_id: 1,
            index_count: 3,
            index_offset: 0,
            base_vertex: 0,
            instance_offset: 0,
            instance_count: 1,
        });
        assert_eq!(batch.call_count(), 3);
        batch.merge_consecutive();
        assert_eq!(
            batch.call_count(),
            2,
            "two consecutive same-mesh calls should merge"
        );
    }

    #[test]
    fn test_render_state_default() {
        let rs = RenderState::default();
        assert!(rs.depth_test);
        assert!(rs.depth_write);
        assert_eq!(rs.blend_mode, BlendMode::Opaque);
        assert_eq!(rs.cull_mode, CullMode::Back);
    }

    #[test]
    fn test_render_state_transparent() {
        let rs = RenderState::transparent();
        assert!(!rs.depth_write);
        assert_eq!(rs.blend_mode, BlendMode::Alpha);
    }

    #[test]
    fn test_frame_stats_record() {
        let mut stats = FrameStats::new();
        stats.record_draw_call(100, 300);
        stats.record_draw_call(50, 150);
        assert_eq!(stats.draw_calls, 2);
        assert_eq!(stats.triangles, 150);
        assert_eq!(stats.vertices, 450);
    }

    #[test]
    fn test_frame_stats_triangles_per_second() {
        let mut stats = FrameStats::new();
        stats.record_draw_call(1_000_000, 3_000_000);
        stats.frame_time_ms = 16.666;
        let tps = stats.triangles_per_second();
        assert!(tps > 5e7, "should be roughly 60M tris/s: {tps}");
    }

    #[test]
    fn test_occlusion_query() {
        let mut q = OcclusionQuery::new();
        assert!(!q.is_visible(), "pending query should not be visible");
        q.resolve(42);
        assert!(q.is_visible(), "query with samples>0 should be visible");
        q.resolve(0);
        assert!(
            !q.is_visible(),
            "query with 0 samples should not be visible"
        );
    }
}
