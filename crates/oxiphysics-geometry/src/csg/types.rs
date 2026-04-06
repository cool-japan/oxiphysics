//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ImplicitSurface;
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// An infinite plane: `dot(p, normal) - d = 0`.
pub struct SdfPlane {
    /// Unit outward normal.
    pub normal: [f64; 3],
    /// Plane offset: points with `dot(p, normal) > d` are outside.
    pub d: f64,
}
impl SdfPlane {
    /// Create a new plane.
    pub fn new(normal: [f64; 3], d: f64) -> Self {
        Self {
            normal: normalize(normal),
            d,
        }
    }
}
/// Smooth CSG difference.
///
/// SDF = smooth_max(sdf_a, -sdf_b, k).
pub struct CsgSmoothDifference {
    /// The base shape.
    pub a: Box<dyn ImplicitSurface>,
    /// The shape to subtract.
    pub b: Box<dyn ImplicitSurface>,
    /// Smoothing factor.
    pub k: f64,
}
impl CsgSmoothDifference {
    /// Create a new smooth difference.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// A cone (infinite half-cone) opening upward from the origin.
///
/// Defined by its half-angle in radians. SDF is exact on the surface.
pub struct SdfCone {
    /// Apex position.
    pub apex: [f64; 3],
    /// Half-angle in radians.
    pub half_angle: f64,
    /// Height of the finite cone.
    pub height: f64,
}
impl SdfCone {
    /// Create a new cone.
    pub fn new(apex: [f64; 3], half_angle: f64, height: f64) -> Self {
        Self {
            apex,
            half_angle,
            height,
        }
    }
}
/// A sphere defined by centre and radius.
pub struct SdfSphere {
    /// World-space centre.
    pub center: [f64; 3],
    /// Sphere radius (must be positive).
    pub radius: f64,
}
impl SdfSphere {
    /// Create a new sphere.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
}
/// CSG intersection: the region inside both `a` and `b`.
///
/// SDF = max(sdf_a, sdf_b).
pub struct CsgIntersection {
    /// First operand.
    pub a: Box<dyn ImplicitSurface>,
    /// Second operand.
    pub b: Box<dyn ImplicitSurface>,
}
impl CsgIntersection {
    /// Create a new intersection.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>) -> Self {
        Self { a, b }
    }
}
/// A node in a CSG tree.
///
/// Each node is either a leaf (wrapping an `ImplicitSurface`) or an interior
/// node combining two subtrees with a `CsgOp`.
pub enum CsgTree {
    /// A leaf primitive.
    Leaf(Box<dyn ImplicitSurface>),
    /// An interior combination node.
    Node {
        /// The boolean operation to apply.
        op: CsgOp,
        /// Left operand.
        left: Box<CsgTree>,
        /// Right operand.
        right: Box<CsgTree>,
    },
}
impl CsgTree {
    /// Construct a leaf node from any `ImplicitSurface`.
    pub fn leaf(s: impl ImplicitSurface + 'static) -> Self {
        CsgTree::Leaf(Box::new(s))
    }
    /// Construct a union node.
    pub fn union(left: CsgTree, right: CsgTree) -> Self {
        CsgTree::Node {
            op: CsgOp::Union,
            left: Box::new(left),
            right: Box::new(right),
        }
    }
    /// Construct an intersection node.
    pub fn intersection(left: CsgTree, right: CsgTree) -> Self {
        CsgTree::Node {
            op: CsgOp::Intersection,
            left: Box::new(left),
            right: Box::new(right),
        }
    }
    /// Construct a difference node (left minus right).
    pub fn difference(left: CsgTree, right: CsgTree) -> Self {
        CsgTree::Node {
            op: CsgOp::Difference,
            left: Box::new(left),
            right: Box::new(right),
        }
    }
}
/// CSG union: the region inside either `a` or `b`.
///
/// SDF = min(sdf_a, sdf_b).
pub struct CsgUnion {
    /// First operand.
    pub a: Box<dyn ImplicitSurface>,
    /// Second operand.
    pub b: Box<dyn ImplicitSurface>,
}
impl CsgUnion {
    /// Create a new union.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>) -> Self {
        Self { a, b }
    }
}
/// A torus centred at the origin lying in the XZ plane.
pub struct SdfTorus {
    /// Major radius (distance from centre to tube centre).
    pub major_radius: f64,
    /// Minor radius (tube radius).
    pub minor_radius: f64,
}
impl SdfTorus {
    /// Create a new torus.
    pub fn new(major_radius: f64, minor_radius: f64) -> Self {
        Self {
            major_radius,
            minor_radius,
        }
    }
}
/// CSG boolean operation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CsgOp {
    /// Union: region inside either child.
    Union,
    /// Intersection: region inside both children.
    Intersection,
    /// Difference: region inside the left child but outside the right child.
    Difference,
}
/// A capsule aligned with the Y axis.
pub struct SdfCapsule {
    /// Center of the capsule.
    pub center: [f64; 3],
    /// Capsule radius.
    pub radius: f64,
    /// Half-height of the cylindrical section.
    pub half_height: f64,
}
impl SdfCapsule {
    /// Create a new capsule.
    pub fn new(center: [f64; 3], radius: f64, half_height: f64) -> Self {
        Self {
            center,
            radius,
            half_height,
        }
    }
}
/// Offset a `CsgTree` surface inward (negative `offset`) or outward (positive).
///
/// Returns a wrapper `SDF f(p) = original_sdf(p) - offset` so that the
/// zero-level set is shifted by `offset` in the normal direction.
///
/// The returned object implements `ImplicitSurface`.
pub struct CsgOffsetSurface {
    /// The original SDF tree (as an evaluated function table over a grid is
    /// impractical, so we store the sampler function as a closure via callback).
    pub(super) inner_sdf: f64,
    /// Offset value (positive = expand outward, negative = shrink inward).
    pub offset: f64,
    /// Bounding box for the original shape (used for gradient estimation).
    pub bbox_min: [f64; 3],
    /// Maximum bounding extent for gradient computation.
    pub bbox_max: [f64; 3],
}
impl CsgOffsetSurface {
    /// Create an offset wrapper (placeholder; use `csg_offset_sdf` for actual evaluation).
    pub fn new(offset: f64, bbox_min: [f64; 3], bbox_max: [f64; 3]) -> Self {
        Self {
            inner_sdf: 0.0,
            offset,
            bbox_min,
            bbox_max,
        }
    }
}
/// A cell in a marching-cubes grid.
#[derive(Debug, Clone)]
pub struct MarchingCell {
    /// Min corner of the cell.
    pub min: [f64; 3],
    /// Cell side length.
    pub step: f64,
    /// SDF values at the 8 corners (index = ix*4 + iy*2 + iz).
    pub corner_values: [f64; 8],
}
impl MarchingCell {
    /// Return `true` if any corner is inside (negative SDF) and any is outside.
    pub fn has_surface(&self) -> bool {
        let has_inside = self.corner_values.iter().any(|&v| v < 0.0);
        let has_outside = self.corner_values.iter().any(|&v| v >= 0.0);
        has_inside && has_outside
    }
    /// Linear interpolation of position along edge between corner `i0` and `i1`.
    fn edge_vertex(&self, i0: usize, i1: usize) -> [f64; 3] {
        let offsets: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [self.step, 0.0, 0.0],
            [self.step, self.step, 0.0],
            [0.0, self.step, 0.0],
            [0.0, 0.0, self.step],
            [self.step, 0.0, self.step],
            [self.step, self.step, self.step],
            [0.0, self.step, self.step],
        ];
        let v0 = self.corner_values[i0];
        let v1 = self.corner_values[i1];
        let t = if (v1 - v0).abs() > 1e-14 {
            v0 / (v0 - v1)
        } else {
            0.5
        };
        let p0 = add(self.min, offsets[i0]);
        let p1 = add(self.min, offsets[i1]);
        lerp3(p0, p1, t)
    }
    /// Extract surface vertices from this cell (one per crossing edge).
    pub fn extract_vertices(&self) -> Vec<[f64; 3]> {
        if !self.has_surface() {
            return vec![];
        }
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
            .filter_map(|&(i0, i1)| {
                let v0 = self.corner_values[i0];
                let v1 = self.corner_values[i1];
                if v0.signum() != v1.signum() {
                    Some(self.edge_vertex(i0, i1))
                } else {
                    None
                }
            })
            .collect()
    }
}
/// An axis-aligned box defined by centre and half-extents.
pub struct SdfBox {
    /// World-space centre.
    pub center: [f64; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}
impl SdfBox {
    /// Create a new box.
    pub fn new(center: [f64; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            half_extents,
        }
    }
}
/// An infinite cylinder along the Y axis.
pub struct SdfCylinder {
    /// World-space centre (on the axis).
    pub center: [f64; 3],
    /// Cylinder radius.
    pub radius: f64,
}
impl SdfCylinder {
    /// Create a new cylinder.
    pub fn new(center: [f64; 3], radius: f64) -> Self {
        Self { center, radius }
    }
}
/// CSG difference: region inside `a` but outside `b`.
///
/// SDF = max(sdf_a, -sdf_b).
pub struct CsgDifference {
    /// The base shape.
    pub a: Box<dyn ImplicitSurface>,
    /// The shape to subtract.
    pub b: Box<dyn ImplicitSurface>,
}
impl CsgDifference {
    /// Create a new difference.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>) -> Self {
        Self { a, b }
    }
}
/// Smooth CSG union using the polynomial smooth-min kernel.
///
/// SDF = smooth_min(sdf_a, sdf_b, k).
pub struct CsgSmoothUnion {
    /// First operand.
    pub a: Box<dyn ImplicitSurface>,
    /// Second operand.
    pub b: Box<dyn ImplicitSurface>,
    /// Smoothing factor (larger -> smoother blend).
    pub k: f64,
}
impl CsgSmoothUnion {
    /// Create a new smooth union.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// Smooth CSG intersection.
///
/// SDF = smooth_max(sdf_a, sdf_b, k).
pub struct CsgSmoothIntersection {
    /// First operand.
    pub a: Box<dyn ImplicitSurface>,
    /// Second operand.
    pub b: Box<dyn ImplicitSurface>,
    /// Smoothing factor.
    pub k: f64,
}
impl CsgSmoothIntersection {
    /// Create a new smooth intersection.
    pub fn new(a: Box<dyn ImplicitSurface>, b: Box<dyn ImplicitSurface>, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// Classification of a point relative to a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlaneSide {
    /// Point is on the front (positive) side.
    Front,
    /// Point is on the back (negative) side.
    Back,
    /// Point is on the plane (within tolerance).
    OnPlane,
}
/// A finite capped cylinder aligned with the Y axis.
pub struct SdfCappedCylinder {
    /// World-space center.
    pub center: [f64; 3],
    /// Cylinder radius.
    pub radius: f64,
    /// Half-height (distance from centre to cap).
    pub half_height: f64,
}
impl SdfCappedCylinder {
    /// Create a new capped cylinder.
    pub fn new(center: [f64; 3], radius: f64, half_height: f64) -> Self {
        Self {
            center,
            radius,
            half_height,
        }
    }
}
/// Rounded box: axis-aligned box with rounded edges.
///
/// SDF = SdfBox.sdf(p) - rounding_radius
pub struct SdfRoundedBox {
    /// World-space centre.
    pub center: [f64; 3],
    /// Half-extents (before rounding).
    pub half_extents: [f64; 3],
    /// Corner rounding radius.
    pub radius: f64,
}
impl SdfRoundedBox {
    /// Create a new rounded box.
    pub fn new(center: [f64; 3], half_extents: [f64; 3], radius: f64) -> Self {
        Self {
            center,
            half_extents,
            radius,
        }
    }
}
