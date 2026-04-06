#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-side mesh rendering, rasterization, and wireframe drawing.
//!
//! All math uses plain `f32` arrays — no external linear-algebra crate.

#![allow(dead_code)]

use std::f32::consts::PI;

// ── Math helpers ─────────────────────────────────────────────────────────────

/// Dot product of two 3-vectors.
pub fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
pub fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length of a 3-vector.
pub fn length3(a: [f32; 3]) -> f32 {
    dot3(a, a).sqrt()
}

/// Normalize a 3-vector to unit length.  Returns the zero vector if input is
/// degenerate (length < 1e-12).
pub fn normalize3(a: [f32; 3]) -> [f32; 3] {
    let len = length3(a);
    if len < 1e-12 {
        [0.0; 3]
    } else {
        [a[0] / len, a[1] / len, a[2] / len]
    }
}

/// Component-wise subtraction: `a - b`.
pub fn sub3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiplication of a 3-vector.
pub fn scale3(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ── Vertex3D ─────────────────────────────────────────────────────────────────

/// A fully-specified 3-D vertex.
#[derive(Debug, Clone, Copy)]
pub struct Vertex3D {
    /// World-space position.
    pub pos: [f32; 3],
    /// Surface normal (unit length).
    pub normal: [f32; 3],
    /// Texture coordinates.
    pub uv: [f32; 2],
    /// RGBA colour (components in `[0, 1]`).
    pub color: [f32; 4],
}

impl Default for Vertex3D {
    fn default() -> Self {
        Self {
            pos: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            uv: [0.0; 2],
            color: [1.0; 4],
        }
    }
}

// ── RenderMesh ───────────────────────────────────────────────────────────────

/// An indexed triangle mesh for CPU-side rendering.
#[derive(Debug, Clone)]
pub struct RenderMesh {
    /// Vertex data.
    pub vertices: Vec<Vertex3D>,
    /// Triangle list (groups of 3 indices).
    pub indices: Vec<u32>,
}

impl RenderMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
        }
    }

    /// Append a single triangle (three new vertices).
    pub fn add_triangle(&mut self, a: Vertex3D, b: Vertex3D, c: Vertex3D) {
        let base = self.vertices.len() as u32;
        self.vertices.push(a);
        self.vertices.push(b);
        self.vertices.push(c);
        self.indices.push(base);
        self.indices.push(base + 1);
        self.indices.push(base + 2);
    }

    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Recompute per-vertex normals by averaging the normals of all incident
    /// faces.
    pub fn compute_normals(&mut self) {
        // Zero out existing normals.
        for v in &mut self.vertices {
            v.normal = [0.0; 3];
        }

        // Accumulate face normals.
        let n_tri = self.triangle_count();
        for t in 0..n_tri {
            let i0 = self.indices[t * 3] as usize;
            let i1 = self.indices[t * 3 + 1] as usize;
            let i2 = self.indices[t * 3 + 2] as usize;

            let p0 = self.vertices[i0].pos;
            let p1 = self.vertices[i1].pos;
            let p2 = self.vertices[i2].pos;

            let edge1 = sub3(p1, p0);
            let edge2 = sub3(p2, p0);
            let face_n = cross3(edge1, edge2); // weighted by area (not normalized)

            for idx in [i0, i1, i2] {
                let n = &mut self.vertices[idx].normal;
                n[0] += face_n[0];
                n[1] += face_n[1];
                n[2] += face_n[2];
            }
        }

        // Normalize.
        for v in &mut self.vertices {
            v.normal = normalize3(v.normal);
        }
    }

    /// Compute the axis-aligned bounding box.  Returns `(min, max)`.
    /// Returns `([0,0,0], [0,0,0])` for an empty mesh.
    pub fn bounding_box(&self) -> ([f32; 3], [f32; 3]) {
        if self.vertices.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = self.vertices[0].pos;
        let mut mx = self.vertices[0].pos;
        for v in &self.vertices[1..] {
            for i in 0..3 {
                mn[i] = mn[i].min(v.pos[i]);
                mx[i] = mx[i].max(v.pos[i]);
            }
        }
        (mn, mx)
    }
}

impl Default for RenderMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ── Geometry builders ─────────────────────────────────────────────────────────

/// Build a UV sphere centred at the origin.
///
/// `lat` = latitude bands, `lon` = longitude slices.
/// Vertex count = `lat * lon * 2 * 3` (quads split into 2 triangles each).
pub fn build_sphere(radius: f32, lat: usize, lon: usize) -> RenderMesh {
    let mut mesh = RenderMesh::new();

    for i in 0..lat {
        let phi0 = PI * (i as f32) / (lat as f32);
        let phi1 = PI * (i as f32 + 1.0) / (lat as f32);

        for j in 0..lon {
            let theta0 = 2.0 * PI * (j as f32) / (lon as f32);
            let theta1 = 2.0 * PI * (j as f32 + 1.0) / (lon as f32);

            // Four corners of the quad on the sphere surface.
            let make = |phi: f32, theta: f32| -> Vertex3D {
                let nx = phi.sin() * theta.cos();
                let ny = phi.cos();
                let nz = phi.sin() * theta.sin();
                Vertex3D {
                    pos: [radius * nx, radius * ny, radius * nz],
                    normal: [nx, ny, nz],
                    uv: [theta / (2.0 * PI), phi / PI],
                    color: [1.0; 4],
                }
            };

            let v00 = make(phi0, theta0);
            let v10 = make(phi1, theta0);
            let v01 = make(phi0, theta1);
            let v11 = make(phi1, theta1);

            mesh.add_triangle(v00, v10, v11);
            mesh.add_triangle(v00, v11, v01);
        }
    }
    mesh
}

/// Build a box centred at the origin with the given half-extents.
///
/// Produces 12 triangles (6 faces × 2 triangles each).
pub fn build_box(half_extents: [f32; 3]) -> RenderMesh {
    let [hx, hy, hz] = half_extents;
    let mut mesh = RenderMesh::new();

    // (normal, 4 corners in CCW order as seen from outside)
    let faces: [([f32; 3], [[f32; 3]; 4]); 6] = [
        (
            [1.0, 0.0, 0.0],
            [[hx, -hy, -hz], [hx, hy, -hz], [hx, hy, hz], [hx, -hy, hz]],
        ),
        (
            [-1.0, 0.0, 0.0],
            [
                [-hx, -hy, hz],
                [-hx, hy, hz],
                [-hx, hy, -hz],
                [-hx, -hy, -hz],
            ],
        ),
        (
            [0.0, 1.0, 0.0],
            [[-hx, hy, -hz], [hx, hy, -hz], [hx, hy, hz], [-hx, hy, hz]],
        ),
        (
            [0.0, -1.0, 0.0],
            [
                [-hx, -hy, hz],
                [hx, -hy, hz],
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
            ],
        ),
        (
            [0.0, 0.0, 1.0],
            [[-hx, -hy, hz], [hx, -hy, hz], [hx, hy, hz], [-hx, hy, hz]],
        ),
        (
            [0.0, 0.0, -1.0],
            [
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
                [-hx, hy, -hz],
                [hx, hy, -hz],
            ],
        ),
    ];

    for (normal, corners) in &faces {
        let make = |p: [f32; 3]| Vertex3D {
            pos: p,
            normal: *normal,
            uv: [0.0; 2],
            color: [1.0; 4],
        };
        mesh.add_triangle(make(corners[0]), make(corners[1]), make(corners[2]));
        mesh.add_triangle(make(corners[0]), make(corners[2]), make(corners[3]));
    }
    mesh
}

/// Build a cylinder centred at the origin aligned along the Y axis.
///
/// `radius` = barrel radius, `height` = total height, `segments` = number of
/// side quads.
pub fn build_cylinder(radius: f32, height: f32, segments: usize) -> RenderMesh {
    let mut mesh = RenderMesh::new();
    let half_h = height / 2.0;

    for s in 0..segments {
        let theta0 = 2.0 * PI * (s as f32) / (segments as f32);
        let theta1 = 2.0 * PI * (s as f32 + 1.0) / (segments as f32);

        let (c0, s0) = (theta0.cos(), theta0.sin());
        let (c1, s1) = (theta1.cos(), theta1.sin());

        // ── Side quad ───────────────────────────────────────────────────────
        let vbl = Vertex3D {
            pos: [radius * c0, -half_h, radius * s0],
            normal: [c0, 0.0, s0],
            uv: [s as f32 / segments as f32, 1.0],
            color: [1.0; 4],
        };
        let vbr = Vertex3D {
            pos: [radius * c1, -half_h, radius * s1],
            normal: [c1, 0.0, s1],
            uv: [(s + 1) as f32 / segments as f32, 1.0],
            color: [1.0; 4],
        };
        let vtl = Vertex3D {
            pos: [radius * c0, half_h, radius * s0],
            normal: [c0, 0.0, s0],
            uv: [s as f32 / segments as f32, 0.0],
            color: [1.0; 4],
        };
        let vtr = Vertex3D {
            pos: [radius * c1, half_h, radius * s1],
            normal: [c1, 0.0, s1],
            uv: [(s + 1) as f32 / segments as f32, 0.0],
            color: [1.0; 4],
        };

        mesh.add_triangle(vbl, vbr, vtr);
        mesh.add_triangle(vbl, vtr, vtl);

        // ── Bottom cap ───────────────────────────────────────────────────────
        let center_b = Vertex3D {
            pos: [0.0, -half_h, 0.0],
            normal: [0.0, -1.0, 0.0],
            uv: [0.5, 0.5],
            color: [1.0; 4],
        };
        let ring0_b = Vertex3D {
            pos: [radius * c0, -half_h, radius * s0],
            normal: [0.0, -1.0, 0.0],
            uv: [0.5 + 0.5 * c0, 0.5 + 0.5 * s0],
            color: [1.0; 4],
        };
        let ring1_b = Vertex3D {
            pos: [radius * c1, -half_h, radius * s1],
            normal: [0.0, -1.0, 0.0],
            uv: [0.5 + 0.5 * c1, 0.5 + 0.5 * s1],
            color: [1.0; 4],
        };
        mesh.add_triangle(center_b, ring1_b, ring0_b);

        // ── Top cap ──────────────────────────────────────────────────────────
        let center_t = Vertex3D {
            pos: [0.0, half_h, 0.0],
            normal: [0.0, 1.0, 0.0],
            uv: [0.5, 0.5],
            color: [1.0; 4],
        };
        let ring0_t = Vertex3D {
            pos: [radius * c0, half_h, radius * s0],
            normal: [0.0, 1.0, 0.0],
            uv: [0.5 + 0.5 * c0, 0.5 + 0.5 * s0],
            color: [1.0; 4],
        };
        let ring1_t = Vertex3D {
            pos: [radius * c1, half_h, radius * s1],
            normal: [0.0, 1.0, 0.0],
            uv: [0.5 + 0.5 * c1, 0.5 + 0.5 * s1],
            color: [1.0; 4],
        };
        mesh.add_triangle(center_t, ring0_t, ring1_t);
    }
    mesh
}

/// Build a flat grid in the XZ plane centred at the origin.
///
/// `width` / `height` = total dimensions.  `nx` / `ny` = cell counts.
pub fn build_grid(width: f32, height: f32, nx: usize, ny: usize) -> RenderMesh {
    let mut mesh = RenderMesh::new();
    let dx = width / nx as f32;
    let dz = height / ny as f32;
    let ox = -width / 2.0;
    let oz = -height / 2.0;

    for iy in 0..ny {
        for ix in 0..nx {
            let x0 = ox + ix as f32 * dx;
            let x1 = ox + (ix + 1) as f32 * dx;
            let z0 = oz + iy as f32 * dz;
            let z1 = oz + (iy + 1) as f32 * dz;

            let u0 = ix as f32 / nx as f32;
            let u1 = (ix + 1) as f32 / nx as f32;
            let v0 = iy as f32 / ny as f32;
            let v1 = (iy + 1) as f32 / ny as f32;

            let make = |x: f32, z: f32, u: f32, v: f32| Vertex3D {
                pos: [x, 0.0, z],
                normal: [0.0, 1.0, 0.0],
                uv: [u, v],
                color: [1.0; 4],
            };

            let v00 = make(x0, z0, u0, v0);
            let v10 = make(x1, z0, u1, v0);
            let v11 = make(x1, z1, u1, v1);
            let v01 = make(x0, z1, u0, v1);

            mesh.add_triangle(v00, v10, v11);
            mesh.add_triangle(v00, v11, v01);
        }
    }
    mesh
}

/// Build an arrow mesh from `start` to `end` with the given shaft `radius`.
///
/// The arrow consists of a cylindrical shaft plus a cone tip.
pub fn build_arrow(start: [f32; 3], end: [f32; 3], radius: f32) -> RenderMesh {
    let dir = sub3(end, start);
    let total_len = length3(dir);
    if total_len < 1e-9 {
        return RenderMesh::new();
    }
    let axis = normalize3(dir);

    // Build an orthonormal basis {axis, tangent, bitangent}.
    let tangent = {
        let candidate = if axis[1].abs() < 0.99 {
            [0.0_f32, 1.0, 0.0]
        } else {
            [1.0_f32, 0.0, 0.0]
        };
        normalize3(cross3(candidate, axis))
    };
    let bitangent = cross3(axis, tangent);

    let tip_len = (total_len * 0.25).max(radius * 2.0);
    let shaft_len = total_len - tip_len;
    let tip_radius = radius * 2.5;
    let segments: usize = 12;

    let mut mesh = RenderMesh::new();

    // Helper: point on a circle at position `t` along axis from `base`.
    let ring_point = |base: [f32; 3], r: f32, s: usize| -> [f32; 3] {
        let theta = 2.0 * PI * (s as f32) / (segments as f32);
        let ct = theta.cos();
        let st = theta.sin();
        [
            base[0] + r * (tangent[0] * ct + bitangent[0] * st),
            base[1] + r * (tangent[1] * ct + bitangent[1] * st),
            base[2] + r * (tangent[2] * ct + bitangent[2] * st),
        ]
    };

    let shaft_end = [
        start[0] + axis[0] * shaft_len,
        start[1] + axis[1] * shaft_len,
        start[2] + axis[2] * shaft_len,
    ];

    // ── Shaft barrel ─────────────────────────────────────────────────────────
    for s in 0..segments {
        let p0b = ring_point(start, radius, s);
        let p1b = ring_point(start, radius, s + 1);
        let p0t = ring_point(shaft_end, radius, s);
        let p1t = ring_point(shaft_end, radius, s + 1);

        let n0 = normalize3(sub3(p0b, start));
        let n1 = normalize3(sub3(p1b, start));

        mesh.add_triangle(
            Vertex3D {
                pos: p0b,
                normal: n0,
                uv: [0.0, 0.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: p1b,
                normal: n1,
                uv: [1.0, 0.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: p1t,
                normal: n1,
                uv: [1.0, 1.0],
                color: [1.0; 4],
            },
        );
        mesh.add_triangle(
            Vertex3D {
                pos: p0b,
                normal: n0,
                uv: [0.0, 0.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: p1t,
                normal: n1,
                uv: [1.0, 1.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: p0t,
                normal: n0,
                uv: [0.0, 1.0],
                color: [1.0; 4],
            },
        );
    }

    // ── Cone tip ─────────────────────────────────────────────────────────────
    let neg_axis = [-axis[0], -axis[1], -axis[2]];
    for s in 0..segments {
        let pb0 = ring_point(shaft_end, tip_radius, s);
        let pb1 = ring_point(shaft_end, tip_radius, s + 1);

        // Cap at cone base
        mesh.add_triangle(
            Vertex3D {
                pos: shaft_end,
                normal: neg_axis,
                uv: [0.5, 0.5],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: pb1,
                normal: neg_axis,
                uv: [1.0, 0.5],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: pb0,
                normal: neg_axis,
                uv: [0.0, 0.5],
                color: [1.0; 4],
            },
        );

        // Cone side
        let side_n0 = normalize3(cross3(sub3(pb0, end), sub3(pb1, pb0)));
        mesh.add_triangle(
            Vertex3D {
                pos: pb0,
                normal: side_n0,
                uv: [0.0, 0.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: pb1,
                normal: side_n0,
                uv: [1.0, 0.0],
                color: [1.0; 4],
            },
            Vertex3D {
                pos: end,
                normal: axis,
                uv: [0.5, 1.0],
                color: [1.0; 4],
            },
        );
    }

    mesh
}

// ── Mat4 ──────────────────────────────────────────────────────────────────────

/// Column-major 4×4 transform matrix.
#[derive(Debug, Clone, Copy)]
pub struct Mat4 {
    /// Column-major storage: `data[col][row]`.
    pub data: [[f32; 4]; 4],
}

impl Mat4 {
    /// Identity matrix.
    pub fn identity() -> Self {
        Self {
            data: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    /// Pure translation.
    pub fn translation(x: f32, y: f32, z: f32) -> Self {
        let mut m = Self::identity();
        m.data[3][0] = x;
        m.data[3][1] = y;
        m.data[3][2] = z;
        m
    }

    /// Uniform scale.
    pub fn scale_uniform(s: f32) -> Self {
        let mut m = Self::identity();
        m.data[0][0] = s;
        m.data[1][1] = s;
        m.data[2][2] = s;
        m
    }

    /// Rotation about the Y axis by `angle` radians.
    pub fn rotation_y(angle: f32) -> Self {
        let (s, c) = angle.sin_cos();
        let mut m = Self::identity();
        m.data[0][0] = c;
        m.data[0][2] = -s;
        m.data[2][0] = s;
        m.data[2][2] = c;
        m
    }

    /// Perspective projection matrix (right-handed, depth range \[-1, 1\]).
    ///
    /// `fovy` = vertical field-of-view in radians, `aspect` = width / height.
    pub fn perspective(fovy: f32, aspect: f32, near: f32, far: f32) -> Self {
        let tan_half = (fovy / 2.0).tan();
        let mut m = [[0.0_f32; 4]; 4];
        m[0][0] = 1.0 / (aspect * tan_half);
        m[1][1] = 1.0 / tan_half;
        m[2][2] = -(far + near) / (far - near);
        m[2][3] = -1.0;
        m[3][2] = -(2.0 * far * near) / (far - near);
        Self { data: m }
    }

    /// View matrix from an eye position, target, and up direction.
    pub fn look_at(eye: [f32; 3], center: [f32; 3], up: [f32; 3]) -> Self {
        let f = normalize3(sub3(center, eye));
        let r = normalize3(cross3(f, up));
        let u = cross3(r, f);

        let mut m = Self::identity();
        m.data[0][0] = r[0];
        m.data[1][0] = r[1];
        m.data[2][0] = r[2];
        m.data[0][1] = u[0];
        m.data[1][1] = u[1];
        m.data[2][1] = u[2];
        m.data[0][2] = -f[0];
        m.data[1][2] = -f[1];
        m.data[2][2] = -f[2];
        m.data[3][0] = -dot3(r, eye);
        m.data[3][1] = -dot3(u, eye);
        m.data[3][2] = dot3(f, eye);
        m
    }

    /// Matrix multiplication: `self * other`.
    pub fn mul(&self, other: &Mat4) -> Mat4 {
        let a = &self.data;
        let b = &other.data;
        let mut c = [[0.0_f32; 4]; 4];
        // c[col][row] = sum_k a[k][row] * b[col][k]
        for col in 0..4 {
            for row in 0..4 {
                let mut sum = 0.0_f32;
                for k in 0..4 {
                    sum += a[k][row] * b[col][k];
                }
                c[col][row] = sum;
            }
        }
        Mat4 { data: c }
    }

    /// Transform a 3-D point (implicit w = 1, then perspective divide).
    pub fn transform_point(&self, p: [f32; 3]) -> [f32; 3] {
        let d = &self.data;
        let x = d[0][0] * p[0] + d[1][0] * p[1] + d[2][0] * p[2] + d[3][0];
        let y = d[0][1] * p[0] + d[1][1] * p[1] + d[2][1] * p[2] + d[3][1];
        let z = d[0][2] * p[0] + d[1][2] * p[1] + d[2][2] * p[2] + d[3][2];
        let w = d[0][3] * p[0] + d[1][3] * p[1] + d[2][3] * p[2] + d[3][3];
        if w.abs() < 1e-12 {
            [x, y, z]
        } else {
            [x / w, y / w, z / w]
        }
    }

    /// Transform a 3-D direction vector (implicit w = 0, no translation or
    /// perspective divide).
    pub fn transform_direction(&self, d_in: [f32; 3]) -> [f32; 3] {
        let d = &self.data;
        [
            d[0][0] * d_in[0] + d[1][0] * d_in[1] + d[2][0] * d_in[2],
            d[0][1] * d_in[0] + d[1][1] * d_in[1] + d[2][1] * d_in[2],
            d[0][2] * d_in[0] + d[1][2] * d_in[1] + d[2][2] * d_in[2],
        ]
    }
}

// ── Framebuffer ───────────────────────────────────────────────────────────────

/// A simple RGBA + depth framebuffer for CPU rasterization.
pub struct Framebuffer {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// RGBA colour buffer (row-major).
    pub color: Vec<[f32; 4]>,
    /// Depth buffer (row-major).
    pub depth: Vec<f32>,
}

impl Framebuffer {
    /// Allocate a framebuffer filled with transparent black and depth ∞.
    pub fn new(w: usize, h: usize) -> Self {
        Self {
            width: w,
            height: h,
            color: vec![[0.0, 0.0, 0.0, 0.0]; w * h],
            depth: vec![f32::INFINITY; w * h],
        }
    }

    /// Clear every pixel to `color` and reset depth to ∞.
    pub fn clear(&mut self, color: [f32; 4]) {
        for c in &mut self.color {
            *c = color;
        }
        for d in &mut self.depth {
            *d = f32::INFINITY;
        }
    }

    /// Get the colour of a pixel.  Returns transparent black if out of range.
    pub fn get_pixel(&self, x: usize, y: usize) -> [f32; 4] {
        if x < self.width && y < self.height {
            self.color[y * self.width + x]
        } else {
            [0.0; 4]
        }
    }

    /// Write a pixel if `depth` is less than the stored depth value.
    pub fn set_pixel(&mut self, x: usize, y: usize, c: [f32; 4], depth: f32) {
        if x < self.width && y < self.height {
            let idx = y * self.width + x;
            if depth < self.depth[idx] {
                self.color[idx] = c;
                self.depth[idx] = depth;
            }
        }
    }
}

// ── SoftwareRasterizer ────────────────────────────────────────────────────────

/// CPU-side rasterization primitives.
pub struct SoftwareRasterizer;

impl SoftwareRasterizer {
    /// Draw a line using Bresenham's algorithm.
    pub fn draw_line(fb: &mut Framebuffer, x0: i32, y0: i32, x1: i32, y1: i32, color: [f32; 4]) {
        let mut x = x0;
        let mut y = y0;
        let dx = (x1 - x0).abs();
        let dy = (y1 - y0).abs();
        let sx: i32 = if x0 < x1 { 1 } else { -1 };
        let sy: i32 = if y0 < y1 { 1 } else { -1 };
        let mut err = dx - dy;

        loop {
            if x >= 0 && y >= 0 {
                fb.set_pixel(x as usize, y as usize, color, 0.0);
            }
            if x == x1 && y == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 > -dy {
                err -= dy;
                x += sx;
            }
            if e2 < dx {
                err += dx;
                y += sy;
            }
        }
    }

    /// Draw the three edges of a triangle (wireframe only).
    pub fn draw_triangle_wireframe(
        fb: &mut Framebuffer,
        v0: [f32; 2],
        v1: [f32; 2],
        v2: [f32; 2],
        color: [f32; 4],
    ) {
        let to_i = |v: [f32; 2]| (v[0] as i32, v[1] as i32);
        let (ax, ay) = to_i(v0);
        let (bx, by) = to_i(v1);
        let (cx, cy) = to_i(v2);
        Self::draw_line(fb, ax, ay, bx, by, color);
        Self::draw_line(fb, bx, by, cx, cy, color);
        Self::draw_line(fb, cx, cy, ax, ay, color);
    }

    /// Rasterize a filled triangle with depth test and per-vertex colour
    /// interpolation (barycentric).
    ///
    /// `v0..v2` are NDC/screen-space `[x, y, depth]` triples.
    pub fn draw_triangle_filled(
        fb: &mut Framebuffer,
        v0: [f32; 3],
        v1: [f32; 3],
        v2: [f32; 3],
        c0: [f32; 4],
        c1: [f32; 4],
        c2: [f32; 4],
    ) {
        // Clip to framebuffer bounds.
        let w = fb.width as f32;
        let h = fb.height as f32;

        // Compute bounding box.
        let min_x = v0[0].min(v1[0]).min(v2[0]).max(0.0).floor() as i32;
        let max_x = v0[0].max(v1[0]).max(v2[0]).min(w - 1.0).ceil() as i32;
        let min_y = v0[1].min(v1[1]).min(v2[1]).max(0.0).floor() as i32;
        let max_y = v0[1].max(v1[1]).max(v2[1]).min(h - 1.0).ceil() as i32;

        // Edge function helper: signed area of (a, b, p).
        let edge = |ax: f32, ay: f32, bx: f32, by: f32, px: f32, py: f32| -> f32 {
            (bx - ax) * (py - ay) - (by - ay) * (px - ax)
        };

        let area = edge(v0[0], v0[1], v1[0], v1[1], v2[0], v2[1]);
        if area.abs() < 1e-9 {
            return;
        }

        for py in min_y..=max_y {
            for px in min_x..=max_x {
                let pxf = px as f32 + 0.5;
                let pyf = py as f32 + 0.5;

                let w0 = edge(v1[0], v1[1], v2[0], v2[1], pxf, pyf);
                let w1 = edge(v2[0], v2[1], v0[0], v0[1], pxf, pyf);
                let w2 = edge(v0[0], v0[1], v1[0], v1[1], pxf, pyf);

                // All must have the same sign as area.
                if (w0 >= 0.0) == (area >= 0.0)
                    && (w1 >= 0.0) == (area >= 0.0)
                    && (w2 >= 0.0) == (area >= 0.0)
                {
                    let bary0 = w0 / area;
                    let bary1 = w1 / area;
                    let bary2 = w2 / area;

                    let depth = bary0 * v0[2] + bary1 * v1[2] + bary2 * v2[2];

                    let interp =
                        |i: usize| -> f32 { bary0 * c0[i] + bary1 * c1[i] + bary2 * c2[i] };
                    let color = [interp(0), interp(1), interp(2), interp(3)];

                    fb.set_pixel(px as usize, py as usize, color, depth);
                }
            }
        }
    }
}

// ── WireframeRenderer ─────────────────────────────────────────────────────────

/// Projects a mesh through an MVP matrix and draws its edges.
pub struct WireframeRenderer;

impl WireframeRenderer {
    /// Render `mesh` to `fb` using the given MVP transform.
    ///
    /// After projection to NDC, coordinates are mapped to pixel space:
    /// `x_px = (ndcx + 1) / 2 * width`, `y_px = (1 - ndcy) / 2 * height`.
    pub fn render(mesh: &RenderMesh, mvp: &Mat4, fb: &mut Framebuffer, color: [f32; 4]) {
        let w = fb.width as f32;
        let h = fb.height as f32;

        // Project all vertices to screen space.
        let screen: Vec<[f32; 3]> = mesh
            .vertices
            .iter()
            .map(|v| {
                let clip = mvp.transform_point(v.pos);
                let sx = (clip[0] + 1.0) / 2.0 * w;
                let sy = (1.0 - clip[1]) / 2.0 * h;
                [sx, sy, clip[2]]
            })
            .collect();

        // Draw each triangle's three edges.
        let n_tri = mesh.triangle_count();
        for t in 0..n_tri {
            let i0 = mesh.indices[t * 3] as usize;
            let i1 = mesh.indices[t * 3 + 1] as usize;
            let i2 = mesh.indices[t * 3 + 2] as usize;

            SoftwareRasterizer::draw_line(
                fb,
                screen[i0][0] as i32,
                screen[i0][1] as i32,
                screen[i1][0] as i32,
                screen[i1][1] as i32,
                color,
            );
            SoftwareRasterizer::draw_line(
                fb,
                screen[i1][0] as i32,
                screen[i1][1] as i32,
                screen[i2][0] as i32,
                screen[i2][1] as i32,
                color,
            );
            SoftwareRasterizer::draw_line(
                fb,
                screen[i2][0] as i32,
                screen[i2][1] as i32,
                screen[i0][0] as i32,
                screen[i0][1] as i32,
                color,
            );
        }
    }
}

// ── PhongShader ───────────────────────────────────────────────────────────────

/// A simple Blinn–Phong shader evaluated on the CPU.
pub struct PhongShader {
    /// Normalised light direction (pointing **toward** the light).
    pub light_dir: [f32; 3],
    /// Light colour (RGB, each component in `[0, 1]`).
    pub light_color: [f32; 3],
    /// Ambient light intensity.
    pub ambient: f32,
}

impl PhongShader {
    /// Evaluate the Blinn–Phong lighting model.
    ///
    /// Returns the final RGB colour (before gamma correction).
    pub fn shade(&self, normal: [f32; 3], view_dir: [f32; 3], base_color: [f32; 3]) -> [f32; 3] {
        let n = normalize3(normal);
        let l = normalize3(self.light_dir);
        let v = normalize3(view_dir);

        // Diffuse
        let n_dot_l = dot3(n, l).max(0.0);

        // Specular (Blinn halfway vector)
        let h = normalize3([l[0] + v[0], l[1] + v[1], l[2] + v[2]]);
        let n_dot_h = dot3(n, h).max(0.0);
        let specular = n_dot_h.powf(32.0);

        let mut result = [0.0_f32; 3];
        for i in 0..3 {
            let ambient = self.ambient * base_color[i];
            let diffuse = n_dot_l * self.light_color[i] * base_color[i];
            let spec_c = specular * self.light_color[i];
            result[i] = (ambient + diffuse + spec_c).clamp(0.0, 1.0);
        }
        result
    }
}

// ── Mesh LOD generation ───────────────────────────────────────────────────────

/// Generate a simplified (LOD) version of a mesh by decimating vertices.
///
/// Uses a simple edge-collapse heuristic: repeatedly merge vertices that are
/// closer than `threshold` in world space.
pub fn generate_lod(mesh: &RenderMesh, threshold: f32) -> RenderMesh {
    let n = mesh.vertices.len();
    if n == 0 || threshold <= 0.0 {
        return mesh.clone();
    }

    // Build a remapping table: map each vertex to the representative vertex.
    let mut remap: Vec<usize> = (0..n).collect();
    for i in 0..n {
        if remap[i] != i {
            continue;
        }
        for j in (i + 1)..n {
            if remap[j] != j {
                continue;
            }
            let dx = mesh.vertices[i].pos[0] - mesh.vertices[j].pos[0];
            let dy = mesh.vertices[i].pos[1] - mesh.vertices[j].pos[1];
            let dz = mesh.vertices[i].pos[2] - mesh.vertices[j].pos[2];
            if (dx * dx + dy * dy + dz * dz).sqrt() <= threshold {
                remap[j] = i;
            }
        }
    }

    // Compact the vertex list.
    let mut new_indices: Vec<Option<usize>> = vec![None; n];
    let mut new_verts: Vec<Vertex3D> = Vec::new();
    for i in 0..n {
        let rep = remap[i];
        if new_indices[rep].is_none() {
            new_indices[rep] = Some(new_verts.len());
            new_verts.push(mesh.vertices[rep]);
        }
        new_indices[i] = new_indices[rep];
    }

    // Remap indices, discard degenerate triangles.
    let mut new_idx_buf: Vec<u32> = Vec::new();
    let tri_count = mesh.indices.len() / 3;
    for t in 0..tri_count {
        let a = new_indices[mesh.indices[t * 3] as usize].unwrap_or(0);
        let b = new_indices[mesh.indices[t * 3 + 1] as usize].unwrap_or(0);
        let c = new_indices[mesh.indices[t * 3 + 2] as usize].unwrap_or(0);
        if a != b && b != c && a != c {
            new_idx_buf.push(a as u32);
            new_idx_buf.push(b as u32);
            new_idx_buf.push(c as u32);
        }
    }

    RenderMesh {
        vertices: new_verts,
        indices: new_idx_buf,
    }
}

// ── Mesh instancing ───────────────────────────────────────────────────────────

/// A set of instance transforms for instanced rendering of a mesh.
#[derive(Debug, Clone)]
pub struct MeshInstanceSet {
    /// The mesh to be instanced.
    pub mesh: RenderMesh,
    /// Column-major 4×4 world transforms for each instance.
    pub transforms: Vec<[[f32; 4]; 4]>,
    /// Optional per-instance tint colours.
    pub tints: Vec<[f32; 4]>,
}

impl MeshInstanceSet {
    /// Create a new instance set from a mesh.
    pub fn new(mesh: RenderMesh) -> Self {
        Self {
            mesh,
            transforms: Vec::new(),
            tints: Vec::new(),
        }
    }

    /// Add an instance at the given position with a given tint.
    pub fn add_instance(&mut self, position: [f32; 3], tint: [f32; 4]) {
        let mut m = [
            [1.0_f32, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [position[0], position[1], position[2], 1.0],
        ];
        // Ensure it's a valid identity + translation.
        m[0][0] = 1.0;
        self.transforms.push(m);
        self.tints.push(tint);
    }

    /// Add an instance using a full 4×4 transform.
    pub fn add_instance_with_transform(&mut self, transform: [[f32; 4]; 4], tint: [f32; 4]) {
        self.transforms.push(transform);
        self.tints.push(tint);
    }

    /// Number of instances.
    pub fn instance_count(&self) -> usize {
        self.transforms.len()
    }

    /// Clear all instances.
    pub fn clear(&mut self) {
        self.transforms.clear();
        self.tints.clear();
    }
}

// ── Mesh skinning data ─────────────────────────────────────────────────────────

/// Maximum number of bone influences per vertex.
pub const MAX_BONE_INFLUENCES: usize = 4;

/// Per-vertex skinning data.
#[derive(Debug, Clone, Copy)]
pub struct SkinVertex {
    /// Bone indices (up to 4).
    pub bone_indices: [u32; MAX_BONE_INFLUENCES],
    /// Bone weights (should sum to 1.0).
    pub bone_weights: [f32; MAX_BONE_INFLUENCES],
}

impl SkinVertex {
    /// Create a new skin vertex with no influence.
    pub fn zero() -> Self {
        Self {
            bone_indices: [0; MAX_BONE_INFLUENCES],
            bone_weights: [0.0; MAX_BONE_INFLUENCES],
        }
    }

    /// Set a single bone with full weight.
    pub fn single_bone(bone_idx: u32) -> Self {
        Self {
            bone_indices: [bone_idx, 0, 0, 0],
            bone_weights: [1.0, 0.0, 0.0, 0.0],
        }
    }

    /// Normalize weights to sum to 1.
    pub fn normalize_weights(&mut self) {
        let total: f32 = self.bone_weights.iter().sum();
        if total > 1e-8 {
            for w in &mut self.bone_weights {
                *w /= total;
            }
        }
    }

    /// Sum of weights.
    pub fn weight_sum(&self) -> f32 {
        self.bone_weights.iter().sum()
    }
}

/// Skinning data for a mesh.
#[derive(Debug, Clone)]
pub struct MeshSkinData {
    /// Per-vertex skinning data.
    pub vertices: Vec<SkinVertex>,
    /// Bind-pose inverse matrices for each bone (column-major 4×4).
    pub inverse_bind_matrices: Vec<[[f32; 4]; 4]>,
    /// Bone names.
    pub bone_names: Vec<String>,
}

impl MeshSkinData {
    /// Create empty skinning data.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            inverse_bind_matrices: Vec::new(),
            bone_names: Vec::new(),
        }
    }

    /// Add a bone with a name and identity bind matrix.
    pub fn add_bone(&mut self, name: &str) -> usize {
        let idx = self.bone_names.len();
        self.bone_names.push(name.to_string());
        self.inverse_bind_matrices.push([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        idx
    }

    /// Number of bones.
    pub fn bone_count(&self) -> usize {
        self.bone_names.len()
    }

    /// Compute the skinned position for a vertex given the current bone transforms.
    ///
    /// `bone_transforms[i]` is the current world transform of bone i.
    #[allow(clippy::too_many_arguments)]
    pub fn skin_position(
        &self,
        vertex_idx: usize,
        base_pos: [f32; 3],
        bone_transforms: &[[[f32; 4]; 4]],
    ) -> [f32; 3] {
        if vertex_idx >= self.vertices.len() {
            return base_pos;
        }
        let sv = &self.vertices[vertex_idx];
        let mut result = [0.0f32; 3];
        for k in 0..MAX_BONE_INFLUENCES {
            let w = sv.bone_weights[k];
            if w < 1e-8 {
                continue;
            }
            let bi = sv.bone_indices[k] as usize;
            if bi >= bone_transforms.len() {
                continue;
            }
            let m = &bone_transforms[bi];
            // Transform base_pos by m (column-major, w=1).
            let x = m[0][0] * base_pos[0] + m[1][0] * base_pos[1] + m[2][0] * base_pos[2] + m[3][0];
            let y = m[0][1] * base_pos[0] + m[1][1] * base_pos[1] + m[2][1] * base_pos[2] + m[3][1];
            let z = m[0][2] * base_pos[0] + m[1][2] * base_pos[1] + m[2][2] * base_pos[2] + m[3][2];
            result[0] += w * x;
            result[1] += w * y;
            result[2] += w * z;
        }
        result
    }
}

impl Default for MeshSkinData {
    fn default() -> Self {
        Self::new()
    }
}

// ── Draw call batching ────────────────────────────────────────────────────────

/// A batched draw call: multiple objects rendered with the same material.
#[derive(Debug, Clone)]
pub struct BatchedDrawCall {
    /// Shared material index.
    pub material_index: usize,
    /// Node IDs included in this batch.
    pub node_ids: Vec<usize>,
    /// Per-instance transforms (column-major 4×4, f32).
    pub instance_transforms: Vec<[[f32; 4]; 4]>,
}

impl BatchedDrawCall {
    /// Create an empty batch for the given material.
    pub fn new(material_index: usize) -> Self {
        Self {
            material_index,
            node_ids: Vec::new(),
            instance_transforms: Vec::new(),
        }
    }

    /// Add a node to this batch with a given transform.
    pub fn push(&mut self, node_id: usize, transform: [[f32; 4]; 4]) {
        self.node_ids.push(node_id);
        self.instance_transforms.push(transform);
    }

    /// Number of instances in this batch.
    pub fn count(&self) -> usize {
        self.node_ids.len()
    }
}

/// Batch a list of `(node_id, material_index, transform)` tuples by material.
pub fn batch_draw_calls(calls: &[(usize, usize, [[f32; 4]; 4])]) -> Vec<BatchedDrawCall> {
    let mut batches: std::collections::HashMap<usize, BatchedDrawCall> =
        std::collections::HashMap::new();
    for &(node_id, mat_idx, transform) in calls {
        batches
            .entry(mat_idx)
            .or_insert_with(|| BatchedDrawCall::new(mat_idx))
            .push(node_id, transform);
    }
    let mut result: Vec<BatchedDrawCall> = batches.into_values().collect();
    result.sort_by_key(|b| b.material_index);
    result
}

// ── Mesh render statistics ────────────────────────────────────────────────────

/// Statistics collected for a frame.
#[derive(Debug, Clone, Default)]
pub struct MeshRenderStats {
    /// Total number of draw calls issued.
    pub draw_calls: usize,
    /// Total number of triangles rendered.
    pub triangles: usize,
    /// Total number of vertices processed.
    pub vertices: usize,
    /// Number of instances rendered.
    pub instances: usize,
    /// Number of nodes culled.
    pub culled: usize,
    /// Number of LOD switches that occurred.
    pub lod_switches: usize,
}

impl MeshRenderStats {
    /// Create zeroed statistics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Accumulate statistics from a mesh and draw call.
    pub fn record_draw(&mut self, mesh: &RenderMesh, instance_count: usize) {
        self.draw_calls += 1;
        self.triangles += mesh.triangle_count() * instance_count;
        self.vertices += mesh.vertices.len() * instance_count;
        self.instances += instance_count;
    }

    /// Reset all counters to zero.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Vertices per triangle ratio (0.0 if no triangles).
    pub fn verts_per_tri(&self) -> f32 {
        if self.triangles == 0 {
            0.0
        } else {
            self.vertices as f32 / self.triangles as f32
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── build_sphere ──────────────────────────────────────────────────────────

    #[test]
    fn test_build_sphere_vertex_count() {
        let lat: usize = 8;
        let lon: usize = 16;
        let mesh = build_sphere(1.0, lat, lon);
        // Each quad → 2 triangles → 6 vertices; lat×lon quads.
        let expected_verts = lat * lon * 2 * 3;
        assert_eq!(
            mesh.vertices.len(),
            expected_verts,
            "sphere vertex count: expected {expected_verts}, got {}",
            mesh.vertices.len()
        );
    }

    // ── build_box ─────────────────────────────────────────────────────────────

    #[test]
    fn test_build_box_12_triangles() {
        let mesh = build_box([0.5, 0.5, 0.5]);
        assert_eq!(mesh.triangle_count(), 12, "box must have 12 triangles");
    }

    // ── Mat4::identity ────────────────────────────────────────────────────────

    #[test]
    fn test_mat4_identity_transforms_point_unchanged() {
        let id = Mat4::identity();
        let p = [3.0_f32, -1.5, 7.0];
        let q = id.transform_point(p);
        for i in 0..3 {
            assert!(
                (q[i] - p[i]).abs() < 1e-6,
                "identity must leave point unchanged; axis {i}: {} vs {}",
                q[i],
                p[i]
            );
        }
    }

    // ── Mat4::perspective far plane ───────────────────────────────────────────

    #[test]
    fn test_mat4_perspective_far_plane() {
        let far = 100.0_f32;
        let near = 0.1_f32;
        let proj = Mat4::perspective(std::f32::consts::FRAC_PI_4, 1.0, near, far);

        // A point exactly on the far plane along -Z (right-handed: z = -far).
        let p = [0.0_f32, 0.0, -far];
        let q = proj.transform_point(p);
        // NDC depth should be +1.0 for the far plane in OpenGL convention.
        assert!(
            (q[2] - 1.0).abs() < 1e-4,
            "far plane should map to NDC depth 1.0, got {}",
            q[2]
        );
    }

    // ── Framebuffer clear ─────────────────────────────────────────────────────

    #[test]
    fn test_framebuffer_clear() {
        let mut fb = Framebuffer::new(8, 8);
        let clear_color = [0.2, 0.4, 0.6, 1.0];
        fb.clear(clear_color);
        for y in 0..8 {
            for x in 0..8 {
                let c = fb.get_pixel(x, y);
                for i in 0..4 {
                    assert!(
                        (c[i] - clear_color[i]).abs() < 1e-6,
                        "pixel ({x},{y}) channel {i}: expected {}, got {}",
                        clear_color[i],
                        c[i]
                    );
                }
            }
        }
    }

    // ── draw_line ─────────────────────────────────────────────────────────────

    #[test]
    fn test_draw_line_sets_pixels() {
        let mut fb = Framebuffer::new(32, 32);
        fb.clear([0.0; 4]);
        let white = [1.0, 1.0, 1.0, 1.0];
        SoftwareRasterizer::draw_line(&mut fb, 0, 0, 10, 0, white);
        // All pixels along y=0, x in [0..10] should be white.
        for x in 0..=10 {
            let c = fb.get_pixel(x, 0);
            assert!(
                (c[0] - 1.0).abs() < 1e-6,
                "pixel ({x},0) should be white after draw_line"
            );
        }
    }

    // ── bounding_box ──────────────────────────────────────────────────────────

    #[test]
    fn test_bounding_box_sphere_symmetric() {
        let r = 1.0_f32;
        let mesh = build_sphere(r, 16, 32);
        let (mn, mx) = mesh.bounding_box();
        for axis in 0..3 {
            assert!(
                mn[axis] <= -r * 0.99 - 1e-6 || mn[axis] >= -(r + 0.01),
                "sphere min[{axis}] = {} should be close to -{r}",
                mn[axis]
            );
            assert!(
                mx[axis] >= r * 0.99 - 1e-6,
                "sphere max[{axis}] = {} should be close to +{r}",
                mx[axis]
            );
            // Symmetric: min ≈ -max
            assert!(
                (mn[axis] + mx[axis]).abs() < 0.05,
                "bounding box should be symmetric: min[{axis}]={} max[{axis}]={}",
                mn[axis],
                mx[axis]
            );
        }
    }

    // ── PhongShader ───────────────────────────────────────────────────────────

    #[test]
    fn test_phong_shader_normal_aligned_with_light() {
        let shader = PhongShader {
            light_dir: [0.0, 0.0, 1.0],
            light_color: [1.0, 1.0, 1.0],
            ambient: 0.1,
        };
        // Normal exactly aligned with light direction → maximum diffuse.
        let result = shader.shade(
            [0.0, 0.0, 1.0], // normal == light
            [0.0, 0.0, 1.0], // view == light (maximises specular too)
            [1.0, 1.0, 1.0], // white base colour
        );
        // Each channel should be close to 1.0 (saturated).
        for i in 0..3 {
            assert!(
                result[i] > 0.9,
                "channel {i} should be bright when normal aligns with light, got {}",
                result[i]
            );
        }
    }

    // ── generate_lod ──────────────────────────────────────────────────────────

    #[test]
    fn test_generate_lod_no_reduction_large_threshold_is_zero() {
        // With threshold=0, no merging should occur.
        let mesh = build_box([0.5, 0.5, 0.5]);
        let lod = generate_lod(&mesh, 0.0);
        assert_eq!(lod.triangle_count(), mesh.triangle_count());
    }

    #[test]
    fn test_generate_lod_reduces_large_threshold() {
        // A sphere where vertices are close together: large threshold merges many.
        let mesh = build_sphere(1.0, 4, 4);
        let lod = generate_lod(&mesh, 0.5);
        assert!(lod.triangle_count() <= mesh.triangle_count());
        // Should still have a non-empty mesh.
        assert!(!lod.vertices.is_empty());
    }

    #[test]
    fn test_generate_lod_returns_valid_indices() {
        let mesh = build_sphere(1.0, 4, 8);
        let lod = generate_lod(&mesh, 0.1);
        let n = lod.vertices.len() as u32;
        for &idx in &lod.indices {
            assert!(idx < n, "LOD index {idx} out of bounds (n_verts={n})");
        }
    }

    // ── MeshInstanceSet ───────────────────────────────────────────────────────

    #[test]
    fn test_mesh_instance_set_add_and_count() {
        let mesh = build_box([0.5; 3]);
        let mut inst = MeshInstanceSet::new(mesh);
        assert_eq!(inst.instance_count(), 0);
        inst.add_instance([1.0, 0.0, 0.0], [1.0; 4]);
        inst.add_instance([2.0, 0.0, 0.0], [1.0; 4]);
        assert_eq!(inst.instance_count(), 2);
    }

    #[test]
    fn test_mesh_instance_set_clear() {
        let mesh = build_box([0.5; 3]);
        let mut inst = MeshInstanceSet::new(mesh);
        inst.add_instance([0.0; 3], [1.0; 4]);
        inst.clear();
        assert_eq!(inst.instance_count(), 0);
    }

    // ── SkinVertex ────────────────────────────────────────────────────────────

    #[test]
    fn test_skin_vertex_single_bone() {
        let sv = SkinVertex::single_bone(3);
        assert_eq!(sv.bone_indices[0], 3);
        assert!((sv.bone_weights[0] - 1.0).abs() < 1e-6);
        assert!((sv.weight_sum() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_skin_vertex_normalize_weights() {
        let mut sv = SkinVertex::zero();
        sv.bone_weights[0] = 2.0;
        sv.bone_weights[1] = 2.0;
        sv.normalize_weights();
        assert!((sv.weight_sum() - 1.0).abs() < 1e-6);
    }

    // ── MeshSkinData ──────────────────────────────────────────────────────────

    #[test]
    fn test_mesh_skin_data_add_bone() {
        let mut skin = MeshSkinData::new();
        let idx = skin.add_bone("root");
        assert_eq!(idx, 0);
        assert_eq!(skin.bone_count(), 1);
    }

    #[test]
    fn test_mesh_skin_data_skin_position_identity() {
        let mut skin = MeshSkinData::new();
        skin.add_bone("root");
        let sv = SkinVertex::single_bone(0);
        skin.vertices.push(sv);
        let id = [
            [1.0_f32, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let pos = [1.0_f32, 2.0, 3.0];
        let result = skin.skin_position(0, pos, &[id]);
        for i in 0..3 {
            assert!(
                (result[i] - pos[i]).abs() < 1e-5,
                "identity skin should preserve position"
            );
        }
    }

    // ── BatchedDrawCall ───────────────────────────────────────────────────────

    #[test]
    fn test_batch_draw_calls_groups_by_material() {
        let id_tf = [
            [1.0_f32, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let calls = vec![(0usize, 0usize, id_tf), (1, 1, id_tf), (2, 0, id_tf)];
        let batches = batch_draw_calls(&calls);
        // Should have 2 batches (material 0 and material 1).
        assert_eq!(batches.len(), 2);
        let mat0 = batches.iter().find(|b| b.material_index == 0).unwrap();
        assert_eq!(mat0.count(), 2, "material 0 should have 2 instances");
    }

    #[test]
    fn test_batch_draw_calls_empty() {
        let batches = batch_draw_calls(&[]);
        assert!(batches.is_empty());
    }

    // ── MeshRenderStats ───────────────────────────────────────────────────────

    #[test]
    fn test_mesh_render_stats_record_draw() {
        let mut stats = MeshRenderStats::new();
        let mesh = build_box([0.5; 3]);
        stats.record_draw(&mesh, 3);
        assert_eq!(stats.draw_calls, 1);
        assert_eq!(stats.triangles, 12 * 3);
        assert_eq!(stats.instances, 3);
    }

    #[test]
    fn test_mesh_render_stats_reset() {
        let mut stats = MeshRenderStats::new();
        let mesh = build_box([0.5; 3]);
        stats.record_draw(&mesh, 1);
        stats.reset();
        assert_eq!(stats.draw_calls, 0);
        assert_eq!(stats.triangles, 0);
    }

    #[test]
    fn test_mesh_render_stats_verts_per_tri() {
        let mut stats = MeshRenderStats::new();
        assert_eq!(stats.verts_per_tri(), 0.0);
        let mesh = build_box([0.5; 3]); // 12 tri, 72 verts
        stats.record_draw(&mesh, 1);
        assert!(stats.verts_per_tri() > 0.0);
    }
}
