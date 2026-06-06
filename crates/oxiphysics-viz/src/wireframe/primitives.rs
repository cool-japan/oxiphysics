// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Basic wireframe shape generators, WireframeMesh, DebugDraw, and contact normals.

use oxiphysics_core::Aabb;
use oxiphysics_core::math::Vec3;

use crate::primitives::{Color, LinePrimitive};

/// Generate 12 line edges for a wireframe box.
pub fn wireframe_box(center: Vec3, half_extents: Vec3, color: Color) -> Vec<LinePrimitive> {
    let hx = half_extents.x;
    let hy = half_extents.y;
    let hz = half_extents.z;
    let corners = [
        center + Vec3::new(-hx, -hy, -hz),
        center + Vec3::new(hx, -hy, -hz),
        center + Vec3::new(hx, hy, -hz),
        center + Vec3::new(-hx, hy, -hz),
        center + Vec3::new(-hx, -hy, hz),
        center + Vec3::new(hx, -hy, hz),
        center + Vec3::new(hx, hy, hz),
        center + Vec3::new(-hx, hy, hz),
    ];
    let edges: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    edges
        .iter()
        .map(|&(a, b)| LinePrimitive {
            start: corners[a],
            end: corners[b],
            color,
        })
        .collect()
}

/// Generate wireframe circles for a sphere (3 great circles: XY, XZ, YZ).
pub fn wireframe_sphere(
    center: Vec3,
    radius: f64,
    segments: usize,
    color: Color,
) -> Vec<LinePrimitive> {
    let mut lines = Vec::new();
    let seg = segments.max(3);
    for i in 0..seg {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((i + 1) as f64) / (seg as f64);
        lines.push(LinePrimitive {
            start: center + Vec3::new(radius * a0.cos(), radius * a0.sin(), 0.0),
            end: center + Vec3::new(radius * a1.cos(), radius * a1.sin(), 0.0),
            color,
        });
    }
    for i in 0..seg {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((i + 1) as f64) / (seg as f64);
        lines.push(LinePrimitive {
            start: center + Vec3::new(radius * a0.cos(), 0.0, radius * a0.sin()),
            end: center + Vec3::new(radius * a1.cos(), 0.0, radius * a1.sin()),
            color,
        });
    }
    for i in 0..seg {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((i + 1) as f64) / (seg as f64);
        lines.push(LinePrimitive {
            start: center + Vec3::new(0.0, radius * a0.cos(), radius * a0.sin()),
            end: center + Vec3::new(0.0, radius * a1.cos(), radius * a1.sin()),
            color,
        });
    }
    lines
}

/// Generate wireframe for a capsule.
pub fn wireframe_capsule(
    center: Vec3,
    half_height: f64,
    radius: f64,
    color: Color,
) -> Vec<LinePrimitive> {
    let mut lines = Vec::new();
    let segments = 16;
    let top_center = center + Vec3::new(0.0, half_height, 0.0);
    let bot_center = center - Vec3::new(0.0, half_height, 0.0);
    for i in 0..segments {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (segments as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((i + 1) as f64) / (segments as f64);
        let (dx0, dz0) = (radius * a0.cos(), radius * a0.sin());
        let (dx1, dz1) = (radius * a1.cos(), radius * a1.sin());
        lines.push(LinePrimitive {
            start: top_center + Vec3::new(dx0, 0.0, dz0),
            end: top_center + Vec3::new(dx1, 0.0, dz1),
            color,
        });
        lines.push(LinePrimitive {
            start: bot_center + Vec3::new(dx0, 0.0, dz0),
            end: bot_center + Vec3::new(dx1, 0.0, dz1),
            color,
        });
    }
    for k in 0..4 {
        let angle = std::f64::consts::FRAC_PI_2 * (k as f64);
        let dx = radius * angle.cos();
        let dz = radius * angle.sin();
        lines.push(LinePrimitive {
            start: top_center + Vec3::new(dx, 0.0, dz),
            end: bot_center + Vec3::new(dx, 0.0, dz),
            color,
        });
    }
    let half_seg = segments / 2;
    for i in 0..half_seg {
        let a0 = std::f64::consts::PI * (i as f64) / (half_seg as f64);
        let a1 = std::f64::consts::PI * ((i + 1) as f64) / (half_seg as f64);
        lines.push(LinePrimitive {
            start: top_center + Vec3::new(radius * a0.cos(), radius * a0.sin(), 0.0),
            end: top_center + Vec3::new(radius * a1.cos(), radius * a1.sin(), 0.0),
            color,
        });
        lines.push(LinePrimitive {
            start: top_center + Vec3::new(0.0, radius * a0.sin(), radius * a0.cos()),
            end: top_center + Vec3::new(0.0, radius * a1.sin(), radius * a1.cos()),
            color,
        });
    }
    for i in 0..half_seg {
        let a0 = std::f64::consts::PI + std::f64::consts::PI * (i as f64) / (half_seg as f64);
        let a1 = std::f64::consts::PI + std::f64::consts::PI * ((i + 1) as f64) / (half_seg as f64);
        lines.push(LinePrimitive {
            start: bot_center + Vec3::new(radius * a0.cos(), radius * a0.sin(), 0.0),
            end: bot_center + Vec3::new(radius * a1.cos(), radius * a1.sin(), 0.0),
            color,
        });
        lines.push(LinePrimitive {
            start: bot_center + Vec3::new(0.0, radius * a0.sin(), radius * a0.cos()),
            end: bot_center + Vec3::new(0.0, radius * a1.sin(), radius * a1.cos()),
            color,
        });
    }
    lines
}

/// Generate wireframe for an axis-aligned bounding box.
pub fn wireframe_aabb(aabb: &Aabb, color: Color) -> Vec<LinePrimitive> {
    let center = aabb.center();
    let half_extents = aabb.half_extents();
    wireframe_box(center, half_extents, color)
}

// ─── Wireframe mesh types ────────────────────────────────────────────────────

/// A vertex in a wireframe mesh.
#[derive(Debug, Clone, Copy)]
pub struct WireframeVertex {
    /// Position in world space.
    pub position: [f64; 3],
    /// RGB color.
    pub color: [f64; 3],
}

/// A line segment in a wireframe mesh, referencing vertex indices.
#[derive(Debug, Clone, Copy)]
pub struct WireframeLine {
    /// Start vertex index.
    pub start: usize,
    /// End vertex index.
    pub end: usize,
}

/// An indexed wireframe mesh.
#[derive(Debug, Clone)]
pub struct WireframeMesh {
    /// Vertices.
    pub vertices: Vec<WireframeVertex>,
    /// Line segments.
    pub lines: Vec<WireframeLine>,
}

impl WireframeMesh {
    /// Create a wireframe box from AABB min/max corners.
    pub fn from_aabb(min: [f64; 3], max: [f64; 3]) -> Self {
        let color = [1.0, 1.0, 1.0];
        let corners: Vec<WireframeVertex> = [
            [min[0], min[1], min[2]],
            [max[0], min[1], min[2]],
            [max[0], max[1], min[2]],
            [min[0], max[1], min[2]],
            [min[0], min[1], max[2]],
            [max[0], min[1], max[2]],
            [max[0], max[1], max[2]],
            [min[0], max[1], max[2]],
        ]
        .iter()
        .map(|&p| WireframeVertex { position: p, color })
        .collect();
        let edge_pairs: [(usize, usize); 12] = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let lines = edge_pairs
            .iter()
            .map(|&(a, b)| WireframeLine { start: a, end: b })
            .collect();
        Self {
            vertices: corners,
            lines,
        }
    }

    /// Create a wireframe sphere from great circles.
    pub fn from_sphere(center: [f64; 3], radius: f64, rings: usize, segments: usize) -> Self {
        let mut vertices = Vec::new();
        let mut lines = Vec::new();
        let color = [1.0, 1.0, 1.0];
        let seg = segments.max(3);
        for ring in 0..rings {
            let phi = std::f64::consts::PI * (ring as f64 + 1.0) / (rings as f64 + 1.0);
            let r = radius * phi.sin();
            let y = center[1] + radius * phi.cos();
            let base = vertices.len();
            for i in 0..seg {
                let theta = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
                vertices.push(WireframeVertex {
                    position: [center[0] + r * theta.cos(), y, center[2] + r * theta.sin()],
                    color,
                });
                lines.push(WireframeLine {
                    start: base + i,
                    end: base + (i + 1) % seg,
                });
            }
        }
        Self { vertices, lines }
    }

    /// Create a wireframe capsule between two points.
    pub fn from_capsule(p0: [f64; 3], p1: [f64; 3], radius: f64, segments: usize) -> Self {
        let seg = segments.max(4);
        let mut vertices = Vec::new();
        let mut lines = Vec::new();
        let color = [1.0, 1.0, 1.0];
        for center in [p0, p1] {
            let base = vertices.len();
            for i in 0..seg {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
                vertices.push(WireframeVertex {
                    position: [
                        center[0] + radius * angle.cos(),
                        center[1],
                        center[2] + radius * angle.sin(),
                    ],
                    color,
                });
                lines.push(WireframeLine {
                    start: base + i,
                    end: base + (i + 1) % seg,
                });
            }
        }
        for k in 0..4 {
            let idx = k * seg / 4;
            lines.push(WireframeLine {
                start: idx,
                end: seg + idx,
            });
        }
        Self { vertices, lines }
    }

    /// Create a wireframe ground grid.
    pub fn from_grid(nx: usize, nz: usize, spacing: f64) -> Self {
        let color = [0.5, 0.5, 0.5];
        let mut vertices = Vec::new();
        let mut lines = Vec::new();
        let half_x = (nx as f64 * spacing) / 2.0;
        let half_z = (nz as f64 * spacing) / 2.0;
        for iz in 0..=nz {
            let z = -half_z + iz as f64 * spacing;
            let a = vertices.len();
            vertices.push(WireframeVertex {
                position: [-half_x, 0.0, z],
                color,
            });
            vertices.push(WireframeVertex {
                position: [half_x, 0.0, z],
                color,
            });
            lines.push(WireframeLine {
                start: a,
                end: a + 1,
            });
        }
        for ix in 0..=nx {
            let x = -half_x + ix as f64 * spacing;
            let a = vertices.len();
            vertices.push(WireframeVertex {
                position: [x, 0.0, -half_z],
                color,
            });
            vertices.push(WireframeVertex {
                position: [x, 0.0, half_z],
                color,
            });
            lines.push(WireframeLine {
                start: a,
                end: a + 1,
            });
        }
        Self { vertices, lines }
    }

    /// Merge another wireframe mesh into this one, returning a new mesh.
    pub fn merge(&self, other: &Self) -> Self {
        let offset = self.vertices.len();
        let mut vertices = self.vertices.clone();
        vertices.extend_from_slice(&other.vertices);
        let mut lines = self.lines.clone();
        for l in &other.lines {
            lines.push(WireframeLine {
                start: l.start + offset,
                end: l.end + offset,
            });
        }
        Self { vertices, lines }
    }
}

// ─── DebugDraw ──────────────────────────────────────────────────────────────

/// Immediate-mode debug line drawing helper.
#[derive(Debug, Clone)]
pub struct DebugDraw {
    /// Collected line segments.
    pub lines: Vec<(WireframeVertex, WireframeVertex)>,
}

impl Default for DebugDraw {
    fn default() -> Self {
        Self::new()
    }
}

impl DebugDraw {
    /// Create a new empty debug drawer.
    pub fn new() -> Self {
        Self { lines: Vec::new() }
    }

    /// Draw a line from `a` to `b` with the given color.
    pub fn draw_line(&mut self, a: [f64; 3], b: [f64; 3], color: [f64; 3]) {
        self.lines.push((
            WireframeVertex { position: a, color },
            WireframeVertex { position: b, color },
        ));
    }

    /// Draw an arrow from `origin` in direction `dir` with the given `length`.
    pub fn draw_arrow(&mut self, origin: [f64; 3], dir: [f64; 3], length: f64, color: [f64; 3]) {
        let mag = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if mag < 1e-30 {
            return;
        }
        let n = [dir[0] / mag, dir[1] / mag, dir[2] / mag];
        let tip = [
            origin[0] + n[0] * length,
            origin[1] + n[1] * length,
            origin[2] + n[2] * length,
        ];
        self.draw_line(origin, tip, color);
        let barb_len = length * 0.2;
        let barb1 = [
            tip[0] - n[0] * barb_len + n[1] * barb_len,
            tip[1] - n[1] * barb_len - n[0] * barb_len,
            tip[2] - n[2] * barb_len,
        ];
        self.draw_line(tip, barb1, color);
    }

    /// Draw a cross (3-axis) at `center` with the given `size`.
    pub fn draw_cross(&mut self, center: [f64; 3], size: f64, color: [f64; 3]) {
        let h = size * 0.5;
        self.draw_line(
            [center[0] - h, center[1], center[2]],
            [center[0] + h, center[1], center[2]],
            color,
        );
        self.draw_line(
            [center[0], center[1] - h, center[2]],
            [center[0], center[1] + h, center[2]],
            color,
        );
        self.draw_line(
            [center[0], center[1], center[2] - h],
            [center[0], center[1], center[2] + h],
            color,
        );
    }

    /// Draw a circle in 3-D space.
    pub fn draw_circle(
        &mut self,
        center: [f64; 3],
        normal: [f64; 3],
        radius: f64,
        segments: usize,
        color: [f64; 3],
    ) {
        let seg = segments.max(3);
        let n_mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if n_mag < 1e-30 {
            return;
        }
        let n = [normal[0] / n_mag, normal[1] / n_mag, normal[2] / n_mag];
        let up = if n[1].abs() < 0.99 {
            [0.0, 1.0, 0.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let t = [
            n[1] * up[2] - n[2] * up[1],
            n[2] * up[0] - n[0] * up[2],
            n[0] * up[1] - n[1] * up[0],
        ];
        let t_mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        let t = [t[0] / t_mag, t[1] / t_mag, t[2] / t_mag];
        let b = [
            n[1] * t[2] - n[2] * t[1],
            n[2] * t[0] - n[0] * t[2],
            n[0] * t[1] - n[1] * t[0],
        ];
        let mut prev = [
            center[0] + radius * t[0],
            center[1] + radius * t[1],
            center[2] + radius * t[2],
        ];
        for i in 1..=seg {
            let angle = 2.0 * std::f64::consts::PI * (i as f64) / (seg as f64);
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            let curr = [
                center[0] + radius * (cos_a * t[0] + sin_a * b[0]),
                center[1] + radius * (cos_a * t[1] + sin_a * b[1]),
                center[2] + radius * (cos_a * t[2] + sin_a * b[2]),
            ];
            self.draw_line(prev, curr, color);
            prev = curr;
        }
    }

    /// Clear all lines.
    pub fn clear(&mut self) {
        self.lines.clear();
    }
    /// Return the number of accumulated lines.
    pub fn line_count(&self) -> usize {
        self.lines.len()
    }
}

/// Generate normal arrows at contact points.
pub fn contact_normal_lines(contacts: &[(Vec3, Vec3, f64)], color: Color) -> Vec<LinePrimitive> {
    contacts
        .iter()
        .map(|&(point, normal, depth)| LinePrimitive {
            start: point,
            end: point + normal * depth,
            color,
        })
        .collect()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wireframe_box_has_12_lines() {
        let lines = wireframe_box(Vec3::zeros(), Vec3::new(1.0, 1.0, 1.0), Color::white());
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn test_wireframe_sphere_line_count() {
        let segments = 8;
        let lines = wireframe_sphere(Vec3::zeros(), 1.0, segments, Color::white());
        assert_eq!(lines.len(), segments * 3);
    }

    #[test]
    fn test_wireframe_aabb() {
        let aabb = Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
        let lines = wireframe_aabb(&aabb, Color::yellow());
        assert_eq!(lines.len(), 12);
    }

    #[test]
    fn test_wireframe_mesh_from_aabb_12_edges() {
        let mesh = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.lines.len(), 12);
    }

    #[test]
    fn test_wireframe_mesh_from_sphere() {
        let mesh = WireframeMesh::from_sphere([0.0; 3], 1.0, 3, 8);
        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.lines.len(), 24);
    }

    #[test]
    fn test_wireframe_mesh_merge() {
        let a = WireframeMesh::from_aabb([0.0; 3], [1.0; 3]);
        let b = WireframeMesh::from_aabb([2.0; 3], [3.0; 3]);
        let merged = a.merge(&b);
        assert_eq!(merged.vertices.len(), 16);
        assert_eq!(merged.lines.len(), 24);
        for l in &merged.lines {
            assert!(l.start < merged.vertices.len());
            assert!(l.end < merged.vertices.len());
        }
    }

    #[test]
    fn test_debug_draw_line() {
        let mut dd = DebugDraw::new();
        dd.draw_line([0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(dd.line_count(), 1);
    }

    #[test]
    fn test_debug_draw_cross() {
        let mut dd = DebugDraw::new();
        dd.draw_cross([0.0; 3], 1.0, [1.0, 1.0, 1.0]);
        assert_eq!(dd.line_count(), 3);
    }

    #[test]
    fn test_debug_draw_circle() {
        let mut dd = DebugDraw::new();
        dd.draw_circle([0.0; 3], [0.0, 1.0, 0.0], 1.0, 16, [0.0, 1.0, 0.0]);
        assert_eq!(dd.line_count(), 16);
    }

    #[test]
    fn test_debug_draw_clear() {
        let mut dd = DebugDraw::new();
        dd.draw_line([0.0; 3], [1.0; 3], [1.0; 3]);
        dd.draw_line([0.0; 3], [2.0; 3], [1.0; 3]);
        dd.clear();
        assert_eq!(dd.line_count(), 0);
    }

    #[test]
    fn test_debug_draw_arrow() {
        let mut dd = DebugDraw::new();
        dd.draw_arrow([0.0; 3], [1.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);
        assert_eq!(dd.line_count(), 2);
    }

    #[test]
    fn test_debug_draw_arrow_zero_dir() {
        let mut dd = DebugDraw::new();
        dd.draw_arrow([0.0; 3], [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]);
        assert_eq!(dd.line_count(), 0);
    }

    #[test]
    fn test_wireframe_mesh_from_capsule() {
        let mesh = WireframeMesh::from_capsule([0.0; 3], [0.0, 2.0, 0.0], 0.5, 8);
        assert_eq!(mesh.vertices.len(), 16);
        assert_eq!(mesh.lines.len(), 20);
    }

    #[test]
    fn test_wireframe_mesh_from_grid() {
        let mesh = WireframeMesh::from_grid(4, 4, 1.0);
        assert_eq!(mesh.lines.len(), 10);
    }

    #[test]
    fn test_contact_normal_lines_count() {
        let contacts = vec![
            (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.1),
            (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 0.2),
        ];
        let lines = contact_normal_lines(&contacts, Color::new(1.0, 0.0, 0.0, 1.0));
        assert_eq!(lines.len(), 2);
    }

    #[test]
    fn test_contact_normal_lines_empty() {
        let lines = contact_normal_lines(&[], Color::new(1.0, 0.0, 0.0, 1.0));
        assert!(lines.is_empty());
    }
}
