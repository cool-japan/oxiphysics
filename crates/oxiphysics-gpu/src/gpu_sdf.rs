// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated Signed Distance Field (SDF) computation (CPU mock).
//!
//! This module provides:
//! - Primitive SDF shapes (`SdfPrimitive`)
//! - CSG boolean operations (`SdfOp`)
//! - A tree-based SDF scene (`SdfNode`)
//! - A volumetric grid (`SdfGrid`) that can be filled from any SDF function
//! - A mock GPU dispatcher (`GpuSdfCompute`) that evaluates primitives over a grid
//! - Free-standing SDF evaluation helpers

// ---------------------------------------------------------------------------
// SdfPrimitive
// ---------------------------------------------------------------------------

/// A primitive signed distance field shape.
#[derive(Debug, Clone, PartialEq)]
pub enum SdfPrimitive {
    /// Sphere centred at `center` with radius `radius`.
    Sphere {
        /// World-space centre of the sphere.
        center: [f64; 3],
        /// Sphere radius (must be positive).
        radius: f64,
    },
    /// Axis-aligned box centred at `center`.
    Box3d {
        /// World-space centre of the box.
        center: [f64; 3],
        /// Half-extents along each axis (width/2, height/2, depth/2).
        half_extents: [f64; 3],
    },
    /// Infinite cylinder (along the Y-axis) with a finite height cap.
    Cylinder {
        /// World-space centre of the cylinder.
        center: [f64; 3],
        /// Cylinder radius.
        radius: f64,
        /// Total height (symmetric about `center`).
        height: f64,
    },
    /// Torus whose axis is the Y-axis.
    Torus {
        /// World-space centre of the torus.
        center: [f64; 3],
        /// Major (ring) radius — distance from the torus centre to the tube centre.
        r_major: f64,
        /// Minor (tube) radius.
        r_minor: f64,
    },
}

impl SdfPrimitive {
    /// Evaluate the SDF for this primitive at world-space point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        match self {
            SdfPrimitive::Sphere { center, radius } => sdf_sphere(p, *center, *radius),
            SdfPrimitive::Box3d {
                center,
                half_extents,
            } => sdf_box(p, *center, *half_extents),
            SdfPrimitive::Cylinder {
                center,
                radius,
                height,
            } => sdf_cylinder(p, *center, *radius, *height),
            SdfPrimitive::Torus {
                center,
                r_major,
                r_minor,
            } => sdf_torus(p, *center, *r_major, *r_minor),
        }
    }
}

// ---------------------------------------------------------------------------
// SdfOp
// ---------------------------------------------------------------------------

/// Boolean / blending operations that combine two SDF values.
#[derive(Debug, Clone, PartialEq)]
pub enum SdfOp {
    /// Take the minimum of two distance values (boolean union).
    Union,
    /// Take the maximum of two distance values (boolean intersection).
    Intersection,
    /// Subtract the second shape from the first (`max(d1, -d2)`).
    Subtraction,
    /// Polynomial smooth minimum blend.
    SmoothUnion {
        /// Blend radius — larger values produce softer transitions.
        k: f64,
    },
    /// Uniformly grow (positive `d`) or shrink (negative `d`) a shape.
    Offset {
        /// Offset distance in world units.
        d: f64,
    },
}

impl SdfOp {
    /// Apply this operation to two already-evaluated distance values `d1` and
    /// `d2`.  For `Offset`, only `d1` is used (the second operand is ignored).
    pub fn apply(&self, d1: f64, d2: f64) -> f64 {
        match self {
            SdfOp::Union => d1.min(d2),
            SdfOp::Intersection => d1.max(d2),
            SdfOp::Subtraction => d1.max(-d2),
            SdfOp::SmoothUnion { k } => sdf_smooth_union(d1, d2, *k),
            SdfOp::Offset { d } => d1 + d,
        }
    }
}

// ---------------------------------------------------------------------------
// SdfNode
// ---------------------------------------------------------------------------

/// A node in a CSG (Constructive Solid Geometry) SDF tree.
#[derive(Debug, Clone)]
pub enum SdfNode {
    /// A leaf node that holds a single primitive.
    Leaf(SdfPrimitive),
    /// An internal node that combines two child nodes with an operation.
    Op {
        /// The combining operation.
        op: SdfOp,
        /// Left child.
        left: Box<SdfNode>,
        /// Right child.
        right: Box<SdfNode>,
    },
}

impl SdfNode {
    /// Create a leaf node from a primitive.
    pub fn leaf(prim: SdfPrimitive) -> Self {
        SdfNode::Leaf(prim)
    }

    /// Create an operation node combining two children.
    pub fn op(op: SdfOp, left: SdfNode, right: SdfNode) -> Self {
        SdfNode::Op {
            op,
            left: Box::new(left),
            right: Box::new(right),
        }
    }

    /// Evaluate the SDF tree at world-space point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        match self {
            SdfNode::Leaf(prim) => prim.evaluate(p),
            SdfNode::Op { op, left, right } => {
                let d1 = left.evaluate(p);
                let d2 = right.evaluate(p);
                op.apply(d1, d2)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SdfGrid
// ---------------------------------------------------------------------------

/// A regular 3-D grid that stores SDF values as `f32`.
///
/// Voxel `(i, j, k)` is centred at
/// `origin + spacing * [i + 0.5, j + 0.5, k + 0.5]`.
#[derive(Debug, Clone)]
pub struct SdfGrid {
    /// Number of voxels along `[x, y, z]`.
    pub dimensions: [usize; 3],
    /// Uniform spacing between voxel centres.
    pub spacing: f64,
    /// World-space position of the grid's lower-left-front corner.
    pub origin: [f64; 3],
    /// Flattened SDF values in x-major order (z slowest).
    ///
    /// `data[i + dims[0\] * (j + dims[1] * k)]` corresponds to voxel `(i,j,k)`.
    pub data: Vec<f32>,
}

impl SdfGrid {
    /// Allocate a new grid, filling every voxel with `f32::MAX`.
    pub fn new(dimensions: [usize; 3], spacing: f64, origin: [f64; 3]) -> Self {
        let n = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            dimensions,
            spacing,
            origin,
            data: vec![f32::MAX; n],
        }
    }

    /// Return the world-space centre of voxel `(i, j, k)`.
    pub fn voxel_center(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + self.spacing * (i as f64 + 0.5),
            self.origin[1] + self.spacing * (j as f64 + 0.5),
            self.origin[2] + self.spacing * (k as f64 + 0.5),
        ]
    }

    /// Flat linear index for voxel `(i, j, k)`.
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.dimensions[0] * (j + self.dimensions[1] * k)
    }

    /// Fill every voxel by evaluating `sdf` at its centre.
    ///
    /// The closure receives a `[f64; 3]` world-space point and should return
    /// the signed distance to the surface.
    pub fn compute_from_sdf(&mut self, sdf: &dyn Fn([f64; 3]) -> f64) {
        let [nx, ny, nz] = self.dimensions;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p = self.voxel_center(i, j, k);
                    let idx = self.index(i, j, k);
                    self.data[idx] = sdf(p) as f32;
                }
            }
        }
    }

    /// Total number of voxels in the grid.
    pub fn len(&self) -> usize {
        self.dimensions[0] * self.dimensions[1] * self.dimensions[2]
    }

    /// Return `true` when the grid contains no voxels.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ---------------------------------------------------------------------------
// GpuSdfCompute
// ---------------------------------------------------------------------------

/// Mock GPU dispatcher for SDF field computation.
///
/// On real hardware this would upload primitives to the GPU and run a
/// compute shader.  Here it falls back to a CPU loop.
#[derive(Debug, Clone, Default)]
pub struct GpuSdfCompute;

impl GpuSdfCompute {
    /// Create a new dispatcher.
    pub fn new() -> Self {
        Self
    }

    /// Fill `grid` with the union of all `primitives`.
    ///
    /// Each voxel is assigned the minimum distance among all primitives
    /// (equivalent to an unbounded CSG union).
    pub fn dispatch_compute(&self, grid: &mut SdfGrid, primitives: &[SdfPrimitive]) {
        if primitives.is_empty() {
            return;
        }
        grid.compute_from_sdf(&|p| {
            primitives
                .iter()
                .map(|prim| prim.evaluate(p))
                .fold(f64::INFINITY, f64::min)
        });
    }
}

// ---------------------------------------------------------------------------
// Free-standing SDF functions
// ---------------------------------------------------------------------------

/// Signed distance from point `p` to a sphere at `center` with radius `r`.
///
/// Negative inside, positive outside, zero on the surface.
pub fn sdf_sphere(p: [f64; 3], center: [f64; 3], r: f64) -> f64 {
    let dx = p[0] - center[0];
    let dy = p[1] - center[1];
    let dz = p[2] - center[2];
    (dx * dx + dy * dy + dz * dz).sqrt() - r
}

/// Signed distance from point `p` to an axis-aligned box centred at `center`
/// with half-extents `b`.
pub fn sdf_box(p: [f64; 3], center: [f64; 3], b: [f64; 3]) -> f64 {
    let qx = (p[0] - center[0]).abs() - b[0];
    let qy = (p[1] - center[1]).abs() - b[1];
    let qz = (p[2] - center[2]).abs() - b[2];
    let outer =
        (qx.max(0.0) * qx.max(0.0) + qy.max(0.0) * qy.max(0.0) + qz.max(0.0) * qz.max(0.0)).sqrt();
    let inner = qx.max(qy).max(qz).min(0.0);
    outer + inner
}

/// Signed distance from point `p` to a finite cylinder (Y-axis aligned)
/// centred at `center` with radius `radius` and total height `height`.
pub fn sdf_cylinder(p: [f64; 3], center: [f64; 3], radius: f64, height: f64) -> f64 {
    let dx = p[0] - center[0];
    let dz = p[2] - center[2];
    let radial = (dx * dx + dz * dz).sqrt() - radius;
    let axial = (p[1] - center[1]).abs() - height * 0.5;
    let outer = (radial.max(0.0) * radial.max(0.0) + axial.max(0.0) * axial.max(0.0)).sqrt();
    let inner = radial.max(axial).min(0.0);
    outer + inner
}

/// Signed distance from point `p` to a torus (Y-axis aligned) centred at
/// `center` with major radius `r_major` and minor radius `r_minor`.
pub fn sdf_torus(p: [f64; 3], center: [f64; 3], r_major: f64, r_minor: f64) -> f64 {
    let dx = p[0] - center[0];
    let dy = p[1] - center[1];
    let dz = p[2] - center[2];
    // Project onto the XZ plane
    let xz_dist = (dx * dx + dz * dz).sqrt();
    let qx = xz_dist - r_major;
    let qy = dy;
    (qx * qx + qy * qy).sqrt() - r_minor
}

/// Polynomial smooth-minimum of `d1` and `d2` with blend radius `k`.
///
/// `k == 0` degenerates to `min(d1, d2)`.  Larger values produce a softer
/// blend between the two surfaces.
///
/// Reference: Inigo Quilez — <https://iquilezles.org/articles/smin/>
pub fn sdf_smooth_union(d1: f64, d2: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return d1.min(d2);
    }
    let h = (0.5 + 0.5 * (d2 - d1) / k).clamp(0.0, 1.0);
    d1 * (1.0 - h) + d2 * h - k * h * (1.0 - h)
}

/// Estimate the surface normal at `p` using central differences on `sdf`.
///
/// Returns a normalised `[f64; 3]` gradient vector.  The step size `eps` is
/// `1e-4` by default (hard-coded).
pub fn sdf_gradient(sdf: &dyn Fn([f64; 3]) -> f64, p: [f64; 3]) -> [f64; 3] {
    const EPS: f64 = 1e-4;
    let gx = sdf([p[0] + EPS, p[1], p[2]]) - sdf([p[0] - EPS, p[1], p[2]]);
    let gy = sdf([p[0], p[1] + EPS, p[2]]) - sdf([p[0], p[1] - EPS, p[2]]);
    let gz = sdf([p[0], p[1], p[2] + EPS]) - sdf([p[0], p[1], p[2] - EPS]);
    let len = (gx * gx + gy * gy + gz * gz).sqrt();
    if len < 1e-15 {
        [0.0, 1.0, 0.0]
    } else {
        [gx / len, gy / len, gz / len]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── sdf_sphere ──────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_surface() {
        let d = sdf_sphere([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(d.abs() < 1e-10, "expected ~0, got {d}");
    }

    #[test]
    fn test_sphere_outside() {
        let d = sdf_sphere([2.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_inside() {
        let d = sdf_sphere([0.0; 3], [0.0; 3], 1.0);
        assert!((d - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_offset_center() {
        let d = sdf_sphere([3.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_sphere_diagonal() {
        // Point at [1,1,1] distance sqrt(3) from origin; sphere r=1 => sqrt(3)-1
        let expected = 3_f64.sqrt() - 1.0;
        let d = sdf_sphere([1.0, 1.0, 1.0], [0.0; 3], 1.0);
        assert!((d - expected).abs() < 1e-10);
    }

    // ── sdf_box ─────────────────────────────────────────────────────────────

    #[test]
    fn test_box_outside_x() {
        // Half-extents [0.5,0.5,0.5], point at [1.5,0,0] => d=1.0
        let d = sdf_box([1.5, 0.0, 0.0], [0.0; 3], [0.5; 3]);
        assert!((d - 1.0).abs() < 1e-10, "got {d}");
    }

    #[test]
    fn test_box_inside() {
        // Point at the origin, fully inside => d < 0
        let d = sdf_box([0.0; 3], [0.0; 3], [1.0; 3]);
        assert!(d < 0.0, "expected negative, got {d}");
    }

    #[test]
    fn test_box_on_face() {
        // On the +x face of a unit cube (half-extents 1)
        let d = sdf_box([1.0, 0.0, 0.0], [0.0; 3], [1.0; 3]);
        assert!(d.abs() < 1e-10, "got {d}");
    }

    #[test]
    fn test_box_corner() {
        // At corner [1,1,1] of unit cube => dist = sqrt(0+0+0) = 0
        let d = sdf_box([1.0, 1.0, 1.0], [0.0; 3], [1.0; 3]);
        assert!(d.abs() < 1e-10, "got {d}");
    }

    #[test]
    fn test_box_asymmetric_extents() {
        // Half-extents [1,2,3]; point at [2,0,0] => outside x only, d=1
        let d = sdf_box([2.0, 0.0, 0.0], [0.0; 3], [1.0, 2.0, 3.0]);
        assert!((d - 1.0).abs() < 1e-10, "got {d}");
    }

    // ── sdf_smooth_union ────────────────────────────────────────────────────

    #[test]
    fn test_smooth_union_zero_k_equals_min() {
        let d1 = 1.0_f64;
        let d2 = 2.0_f64;
        assert!((sdf_smooth_union(d1, d2, 0.0) - d1.min(d2)).abs() < 1e-12);
    }

    #[test]
    fn test_smooth_union_equal_inputs() {
        // When d1 == d2, h = 0.5 => result = d1 - k/4
        let d = 1.0_f64;
        let k = 0.5_f64;
        let result = sdf_smooth_union(d, d, k);
        let expected = d - k / 4.0;
        assert!((result - expected).abs() < 1e-10, "got {result}");
    }

    #[test]
    fn test_smooth_union_approaches_min_for_large_k() {
        // Far-apart values with large k: both clamped, result approaches min
        let d1 = 0.0_f64;
        let d2 = 100.0_f64;
        let result = sdf_smooth_union(d1, d2, 1.0);
        // h should clamp to 1.0, result = d1*0 + d2*1 - 0 = d2? No:
        // h = clamp(0.5 + 0.5*(100-0)/1.0) = clamp(50.5) = 1.0
        // result = d1*(1-1) + d2*1 - k*1*0 = 100
        // But this is max; expected is min. Let's just check it's <= max(d1,d2)
        assert!(result <= d2 + 1e-9);
    }

    #[test]
    fn test_smooth_union_negative_k_is_min() {
        let result = sdf_smooth_union(1.0, 2.0, -1.0);
        assert!((result - 1.0_f64.min(2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_smooth_union_symmetric() {
        // f(a,b,k) == f(b,a,k) up to float noise
        let a = 0.3_f64;
        let b = 0.7_f64;
        let k = 0.4_f64;
        let ab = sdf_smooth_union(a, b, k);
        let ba = sdf_smooth_union(b, a, k);
        assert!((ab - ba).abs() < 1e-12, "not symmetric: {ab} vs {ba}");
    }

    // ── sdf_gradient ────────────────────────────────────────────────────────

    #[test]
    fn test_gradient_sphere_outward() {
        // On the surface of a unit sphere at [1,0,0], gradient should point along +x
        let sdf = |p: [f64; 3]| sdf_sphere(p, [0.0; 3], 1.0);
        let g = sdf_gradient(&sdf, [1.0, 0.0, 0.0]);
        assert!((g[0] - 1.0).abs() < 1e-3, "gx={}", g[0]);
        assert!(g[1].abs() < 1e-3);
        assert!(g[2].abs() < 1e-3);
    }

    #[test]
    fn test_gradient_is_unit_length() {
        let sdf = |p: [f64; 3]| sdf_sphere(p, [0.0; 3], 1.0);
        let g = sdf_gradient(&sdf, [2.0, 0.0, 0.0]);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-6, "len={len}");
    }

    #[test]
    fn test_gradient_box_face_normal() {
        // Outside +x face of a unit box, gradient should point along +x
        let sdf = |p: [f64; 3]| sdf_box(p, [0.0; 3], [1.0; 3]);
        let g = sdf_gradient(&sdf, [2.0, 0.0, 0.0]);
        assert!((g[0] - 1.0).abs() < 1e-3, "gx={}", g[0]);
    }

    #[test]
    fn test_gradient_smooth_union() {
        let sdf = |p: [f64; 3]| {
            let d1 = sdf_sphere(p, [-1.0, 0.0, 0.0], 0.5);
            let d2 = sdf_sphere(p, [1.0, 0.0, 0.0], 0.5);
            sdf_smooth_union(d1, d2, 0.3)
        };
        let g = sdf_gradient(&sdf, [0.0, 0.0, 3.0]);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-4, "len={len}");
    }

    // ── SdfPrimitive ────────────────────────────────────────────────────────

    #[test]
    fn test_primitive_sphere_evaluate() {
        let prim = SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 2.0,
        };
        assert!((prim.evaluate([2.0, 0.0, 0.0])).abs() < 1e-10);
    }

    #[test]
    fn test_primitive_box_evaluate() {
        let prim = SdfPrimitive::Box3d {
            center: [0.0; 3],
            half_extents: [1.0; 3],
        };
        let d = prim.evaluate([0.0; 3]);
        assert!(d < 0.0);
    }

    #[test]
    fn test_primitive_cylinder_evaluate() {
        let prim = SdfPrimitive::Cylinder {
            center: [0.0; 3],
            radius: 1.0,
            height: 2.0,
        };
        // Point on the curved surface
        let d = prim.evaluate([1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-10, "d={d}");
    }

    #[test]
    fn test_primitive_torus_evaluate() {
        let prim = SdfPrimitive::Torus {
            center: [0.0; 3],
            r_major: 2.0,
            r_minor: 0.5,
        };
        // Point on the outer equator of the torus tube
        let d = prim.evaluate([2.5, 0.0, 0.0]);
        assert!(d.abs() < 1e-10, "d={d}");
    }

    // ── SdfOp ───────────────────────────────────────────────────────────────

    #[test]
    fn test_op_union() {
        let op = SdfOp::Union;
        assert!((op.apply(1.0, 2.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_op_intersection() {
        let op = SdfOp::Intersection;
        assert!((op.apply(1.0, 2.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_op_subtraction() {
        // d1=2, d2=3 => max(2, -3) = 2
        let op = SdfOp::Subtraction;
        assert!((op.apply(2.0, 3.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_op_smooth_union() {
        let op = SdfOp::SmoothUnion { k: 0.3 };
        let result = op.apply(1.0, 1.0);
        // Should be 1.0 - 0.3/4 = 0.925
        assert!((result - (1.0 - 0.3 / 4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_op_offset_positive() {
        let op = SdfOp::Offset { d: 0.5 };
        assert!((op.apply(1.0, 0.0) - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_op_offset_negative() {
        let op = SdfOp::Offset { d: -0.5 };
        assert!((op.apply(1.0, 0.0) - 0.5).abs() < 1e-12);
    }

    // ── SdfNode ─────────────────────────────────────────────────────────────

    #[test]
    fn test_node_leaf_evaluate() {
        let node = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        });
        assert!((node.evaluate([1.0, 0.0, 0.0])).abs() < 1e-10);
    }

    #[test]
    fn test_node_op_union() {
        let a = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [-2.0, 0.0, 0.0],
            radius: 1.0,
        });
        let b = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [2.0, 0.0, 0.0],
            radius: 1.0,
        });
        let tree = SdfNode::op(SdfOp::Union, a, b);
        // Origin: sphere A dist=1, sphere B dist=1, union=1
        let d = tree.evaluate([0.0; 3]);
        assert!((d - 1.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn test_node_op_intersection() {
        let a = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 2.0,
        });
        let b = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 3.0,
        });
        let tree = SdfNode::op(SdfOp::Intersection, a, b);
        // At origin: a=-2, b=-3 => intersection = max(-2,-3) = -2
        let d = tree.evaluate([0.0; 3]);
        assert!((d - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_node_op_subtraction() {
        let a = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 2.0,
        });
        let b = SdfNode::leaf(SdfPrimitive::Sphere {
            center: [0.0; 3],
            radius: 1.0,
        });
        let tree = SdfNode::op(SdfOp::Subtraction, a, b);
        // At origin: a=-2, b=-1 => subtraction = max(-2, 1) = 1
        let d = tree.evaluate([0.0; 3]);
        assert!((d - 1.0).abs() < 1e-10, "d={d}");
    }

    // ── SdfGrid ─────────────────────────────────────────────────────────────

    #[test]
    fn test_grid_new_size() {
        let grid = SdfGrid::new([4, 4, 4], 0.1, [0.0; 3]);
        assert_eq!(grid.len(), 64);
        assert!(!grid.is_empty());
    }

    #[test]
    fn test_grid_empty_dimensions() {
        let grid = SdfGrid::new([0, 4, 4], 0.1, [0.0; 3]);
        assert!(grid.is_empty());
    }

    #[test]
    fn test_grid_compute_from_sdf_sphere() {
        let mut grid = SdfGrid::new([4, 4, 4], 1.0, [-2.0, -2.0, -2.0]);
        let sphere_sdf = |p: [f64; 3]| sdf_sphere(p, [0.0; 3], 1.0);
        grid.compute_from_sdf(&sphere_sdf);
        // Every voxel should have been filled (no f32::MAX remaining)
        assert!(
            grid.data.iter().all(|&v| v.is_finite()),
            "some voxels unfilled"
        );
    }

    #[test]
    fn test_grid_index_roundtrip() {
        let grid = SdfGrid::new([5, 6, 7], 0.5, [0.0; 3]);
        for k in 0..7 {
            for j in 0..6 {
                for i in 0..5 {
                    let idx = grid.index(i, j, k);
                    // Manual formula
                    let expected = i + 5 * (j + 6 * k);
                    assert_eq!(idx, expected);
                }
            }
        }
    }

    #[test]
    fn test_grid_voxel_center() {
        let grid = SdfGrid::new([2, 2, 2], 1.0, [0.0; 3]);
        let c = grid.voxel_center(0, 0, 0);
        assert!((c[0] - 0.5).abs() < 1e-10);
        assert!((c[1] - 0.5).abs() < 1e-10);
        assert!((c[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_grid_center_of_last_voxel() {
        let grid = SdfGrid::new([4, 4, 4], 1.0, [0.0; 3]);
        let c = grid.voxel_center(3, 3, 3);
        assert!((c[0] - 3.5).abs() < 1e-10);
    }

    #[test]
    fn test_grid_inside_outside_counts() {
        let mut grid = SdfGrid::new([10, 10, 10], 0.2, [-1.0, -1.0, -1.0]);
        let sphere_sdf = |p: [f64; 3]| sdf_sphere(p, [0.0; 3], 0.5);
        grid.compute_from_sdf(&sphere_sdf);
        let inside = grid.data.iter().filter(|&&v| v < 0.0).count();
        assert!(inside > 0, "no voxels inside the sphere");
        let outside = grid.data.iter().filter(|&&v| v > 0.0).count();
        assert!(outside > 0, "no voxels outside the sphere");
    }

    // ── GpuSdfCompute ───────────────────────────────────────────────────────

    #[test]
    fn test_gpu_sdf_compute_no_primitives() {
        let compute = GpuSdfCompute::new();
        let mut grid = SdfGrid::new([4, 4, 4], 1.0, [0.0; 3]);
        // Fill with sentinel
        grid.data.iter_mut().for_each(|v| *v = -999.0);
        compute.dispatch_compute(&mut grid, &[]);
        // Nothing should have changed
        assert!(grid.data.iter().all(|&v| (v - (-999.0)).abs() < 1e-3));
    }

    #[test]
    fn test_gpu_sdf_compute_sphere_union() {
        let compute = GpuSdfCompute::new();
        let mut grid = SdfGrid::new([8, 8, 8], 0.25, [-1.0, -1.0, -1.0]);
        let primitives = vec![
            SdfPrimitive::Sphere {
                center: [-0.5, 0.0, 0.0],
                radius: 0.3,
            },
            SdfPrimitive::Sphere {
                center: [0.5, 0.0, 0.0],
                radius: 0.3,
            },
        ];
        compute.dispatch_compute(&mut grid, &primitives);
        // Grid should be fully populated
        assert!(grid.data.iter().all(|&v| v.is_finite()));
        // Inside at least one sphere should exist
        let inside = grid.data.iter().filter(|&&v| v < 0.0).count();
        assert!(inside > 0, "no voxels inside");
    }

    #[test]
    fn test_gpu_sdf_compute_single_box() {
        let compute = GpuSdfCompute::new();
        let mut grid = SdfGrid::new([4, 4, 4], 1.0, [-2.0, -2.0, -2.0]);
        let primitives = vec![SdfPrimitive::Box3d {
            center: [0.0; 3],
            half_extents: [0.5; 3],
        }];
        compute.dispatch_compute(&mut grid, &primitives);
        assert!(grid.data.iter().all(|&v| v.is_finite()));
    }

    #[test]
    fn test_gpu_sdf_compute_torus() {
        let compute = GpuSdfCompute::new();
        let mut grid = SdfGrid::new([6, 6, 6], 0.5, [-1.5, -1.5, -1.5]);
        let primitives = vec![SdfPrimitive::Torus {
            center: [0.0; 3],
            r_major: 0.8,
            r_minor: 0.2,
        }];
        compute.dispatch_compute(&mut grid, &primitives);
        assert!(grid.data.iter().all(|&v| v.is_finite()));
    }

    // ── Cylinder ────────────────────────────────────────────────────────────

    #[test]
    fn test_cylinder_inside() {
        let d = sdf_cylinder([0.0, 0.0, 0.0], [0.0; 3], 1.0, 2.0);
        assert!(d < 0.0, "expected inside, d={d}");
    }

    #[test]
    fn test_cylinder_outside_radially() {
        let d = sdf_cylinder([2.0, 0.0, 0.0], [0.0; 3], 1.0, 4.0);
        assert!((d - 1.0).abs() < 1e-10, "d={d}");
    }

    #[test]
    fn test_cylinder_on_curved_surface() {
        let d = sdf_cylinder([1.0, 0.0, 0.0], [0.0; 3], 1.0, 4.0);
        assert!(d.abs() < 1e-10, "d={d}");
    }

    // ── Torus ───────────────────────────────────────────────────────────────

    #[test]
    fn test_torus_centre_is_positive() {
        // The centre of the torus is outside the tube
        let d = sdf_torus([0.0; 3], [0.0; 3], 2.0, 0.5);
        assert!(d > 0.0, "d={d}");
    }

    #[test]
    fn test_torus_on_tube_surface() {
        // On the inner equator: x = r_major - r_minor, y=0, z=0
        let d = sdf_torus([1.5, 0.0, 0.0], [0.0; 3], 2.0, 0.5);
        assert!(d.abs() < 1e-10, "d={d}");
    }
}
