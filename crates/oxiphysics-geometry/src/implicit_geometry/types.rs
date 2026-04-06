//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

/// Smooth intersection using polynomial blend.
#[derive(Debug, Clone)]
pub struct SdfSmoothIntersection<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Blend radius.
    pub k: f64,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfSmoothIntersection<A, B> {
    /// Create a smooth intersection.
    pub fn new(a: A, b: B, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// Smooth union using the polynomial C¹ blend of Quilez.
///
/// `k` controls the blend radius; `k = 0` degenerates to hard union.
#[derive(Debug, Clone)]
pub struct SdfSmoothUnion<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Blend radius (≥ 0).
    pub k: f64,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfSmoothUnion<A, B> {
    /// Create a smooth union.
    pub fn new(a: A, b: B, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// SDF for a hexagonal prism aligned with the Y axis.
///
/// `radius` is the circumscribed circle radius; `half_height` is along Y.
#[derive(Debug, Clone, Copy)]
pub struct SdfHexagonalPrism {
    /// Circumscribed radius of the hexagonal cross-section \[m\].
    pub radius: f64,
    /// Half-height along Y \[m\].
    pub half_height: f64,
}
#[allow(dead_code)]
impl SdfHexagonalPrism {
    /// Create a hexagonal prism SDF.
    pub fn new(radius: f64, half_height: f64) -> Self {
        Self {
            radius,
            half_height,
        }
    }
}
/// Boolean operation type for the SDF expression tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp {
    /// Hard union.
    Union,
    /// Hard intersection.
    Intersection,
    /// Difference (left minus right).
    Difference,
    /// Smooth union (blend radius stored separately).
    SmoothUnion,
    /// Smooth intersection.
    SmoothIntersection,
    /// Smooth difference.
    SmoothDifference,
}
/// A single contact point between two bodies.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// Contact position in world space.
    pub position: [f64; 3],
    /// Contact normal (from B toward A).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
}
/// Translate an SDF by a constant offset.
#[derive(Debug, Clone)]
pub struct SdfTranslate<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Translation vector.
    pub offset: [f64; 3],
}
#[allow(dead_code)]
impl<S: Sdf> SdfTranslate<S> {
    /// Create a translated SDF.
    pub fn new(inner: S, offset: [f64; 3]) -> Self {
        Self { inner, offset }
    }
}
/// Infinite repetition of an SDF along all three axes.
#[derive(Debug, Clone)]
pub struct SdfRepeat<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Cell period in X.
    pub period_x: f64,
    /// Cell period in Y.
    pub period_y: f64,
    /// Cell period in Z.
    pub period_z: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfRepeat<S> {
    /// Create an infinitely repeated SDF.
    pub fn new(inner: S, period_x: f64, period_y: f64, period_z: f64) -> Self {
        Self {
            inner,
            period_x,
            period_y,
            period_z,
        }
    }
}
/// Scale an SDF uniformly (preserves exact signed distance).
#[derive(Debug, Clone)]
pub struct SdfScale<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Uniform scale factor (> 0).
    pub factor: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfScale<S> {
    /// Create a scaled SDF.
    pub fn new(inner: S, factor: f64) -> Self {
        Self { inner, factor }
    }
}
/// SDF for a torus in the XZ plane centred at the origin.
#[derive(Debug, Clone, Copy)]
pub struct SdfTorus {
    /// Major radius (centre to tube centre) \[m\].
    pub major: f64,
    /// Minor radius (tube radius) \[m\].
    pub minor: f64,
}
#[allow(dead_code)]
impl SdfTorus {
    /// Create a torus SDF.
    pub fn new(major: f64, minor: f64) -> Self {
        Self { major, minor }
    }
}
/// SDF for an infinite plane defined by `n·p = d` where `n` is the unit normal.
///
/// Points satisfying `n·p > d` are *outside* (positive side).
#[derive(Debug, Clone, Copy)]
pub struct SdfPlane {
    /// Unit outward normal.
    pub normal: [f64; 3],
    /// Signed distance of the plane from the origin along `normal`.
    pub offset: f64,
}
#[allow(dead_code)]
impl SdfPlane {
    /// Create a plane SDF from a (not necessarily unit) normal and offset.
    pub fn new(normal: [f64; 3], offset: f64) -> Self {
        Self {
            normal: norm(normal),
            offset,
        }
    }
}
/// Smooth difference (subtract `B` from `A`) with polynomial blend.
#[derive(Debug, Clone)]
pub struct SdfSmoothDifference<A, B> {
    /// Shape to subtract from.
    pub a: A,
    /// Shape to subtract.
    pub b: B,
    /// Blend radius.
    pub k: f64,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfSmoothDifference<A, B> {
    /// Create a smooth difference.
    pub fn new(a: A, b: B, k: f64) -> Self {
        Self { a, b, k }
    }
}
/// Apply torsion (twist around Y axis) to an SDF.
///
/// `strength` is twist angle per unit of Y \[radians/m\].
#[derive(Debug, Clone)]
pub struct SdfTwist<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Twist rate \[rad/m\].
    pub strength: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfTwist<S> {
    /// Create a twisted SDF.
    pub fn new(inner: S, strength: f64) -> Self {
        Self { inner, strength }
    }
}
/// A dense 3-D grid of SDF samples.
#[derive(Debug, Clone)]
pub struct SdfGrid {
    /// Number of grid cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z.
    pub nz: usize,
    /// World-space coordinates of cell (0,0,0).
    pub origin: [f64; 3],
    /// Uniform voxel spacing.
    pub spacing: f64,
    /// Flat SDF samples in (z, y, x) major order.
    pub data: Vec<f64>,
}
#[allow(dead_code)]
impl SdfGrid {
    /// Create an empty grid of given resolution and spacing.
    pub fn new(nx: usize, ny: usize, nz: usize, origin: [f64; 3], spacing: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin,
            spacing,
            data: vec![f64::INFINITY; nx * ny * nz],
        }
    }
    /// Sample an SDF and fill the grid.
    pub fn from_sdf<S: Sdf>(
        sdf: &S,
        nx: usize,
        ny: usize,
        nz: usize,
        origin: [f64; 3],
        spacing: f64,
    ) -> Self {
        let mut grid = Self::new(nx, ny, nz, origin, spacing);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = grid.cell_center(ix, iy, iz);
                    grid.data[iz * ny * nx + iy * nx + ix] = sdf.dist(p);
                }
            }
        }
        grid
    }
    /// World-space position of cell centre `(ix, iy, iz)`.
    pub fn cell_center(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.spacing,
            self.origin[1] + iy as f64 * self.spacing,
            self.origin[2] + iz as f64 * self.spacing,
        ]
    }
    /// Get the SDF value at integer index `(ix, iy, iz)`.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[iz * self.ny * self.nx + iy * self.nx + ix]
    }
    /// Set the SDF value at integer index `(ix, iy, iz)`.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, v: f64) {
        self.data[iz * self.ny * self.nx + iy * self.nx + ix] = v;
    }
    /// Trilinearly interpolate the SDF at world-space point `p`.
    pub fn interpolate(&self, p: [f64; 3]) -> f64 {
        let fx = (p[0] - self.origin[0]) / self.spacing;
        let fy = (p[1] - self.origin[1]) / self.spacing;
        let fz = (p[2] - self.origin[2]) / self.spacing;
        let ix = (fx as usize).min(self.nx.saturating_sub(2));
        let iy = (fy as usize).min(self.ny.saturating_sub(2));
        let iz = (fz as usize).min(self.nz.saturating_sub(2));
        let tx = (fx - ix as f64).clamp(0.0, 1.0);
        let ty = (fy - iy as f64).clamp(0.0, 1.0);
        let tz = (fz - iz as f64).clamp(0.0, 1.0);
        let v000 = self.get(ix, iy, iz);
        let v100 = self.get(ix + 1, iy, iz);
        let v010 = self.get(ix, iy + 1, iz);
        let v110 = self.get(ix + 1, iy + 1, iz);
        let v001 = self.get(ix, iy, iz + 1);
        let v101 = self.get(ix + 1, iy, iz + 1);
        let v011 = self.get(ix, iy + 1, iz + 1);
        let v111 = self.get(ix + 1, iy + 1, iz + 1);
        let c00 = v000 * (1.0 - tx) + v100 * tx;
        let c10 = v010 * (1.0 - tx) + v110 * tx;
        let c01 = v001 * (1.0 - tx) + v101 * tx;
        let c11 = v011 * (1.0 - tx) + v111 * tx;
        let c0 = c00 * (1.0 - ty) + c10 * ty;
        let c1 = c01 * (1.0 - ty) + c11 * ty;
        c0 * (1.0 - tz) + c1 * tz
    }
    /// Apply morphological dilation (grow surface by `delta`).
    pub fn dilate(&mut self, delta: f64) {
        for v in self.data.iter_mut() {
            *v -= delta;
        }
    }
    /// Apply morphological erosion (shrink surface by `delta`).
    pub fn erode(&mut self, delta: f64) {
        for v in self.data.iter_mut() {
            *v += delta;
        }
    }
}
/// SDF for a rounded cylinder aligned with the Y axis.
///
/// Cylinder of `radius` and `half_height`, with hemisphere caps of the same
/// radius (i.e. a capsule-like shape but with flat lateral surface).
#[derive(Debug, Clone, Copy)]
pub struct SdfRoundedCylinder {
    /// Cylinder radius \[m\].
    pub radius: f64,
    /// Half-height of the cylindrical body \[m\].
    pub half_height: f64,
    /// Corner rounding radius \[m\].
    pub rounding: f64,
}
#[allow(dead_code)]
impl SdfRoundedCylinder {
    /// Create a rounded cylinder SDF.
    pub fn new(radius: f64, half_height: f64, rounding: f64) -> Self {
        Self {
            radius,
            half_height,
            rounding,
        }
    }
}
/// SDF intersection of two shapes.
#[derive(Debug, Clone)]
pub struct SdfIntersection<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfIntersection<A, B> {
    /// Create an intersection.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}
/// Extrude a 2-D XZ-plane profile along the Y axis.
///
/// The wrapped function takes `[x, z]` and returns the 2-D signed distance.
pub struct SdfExtrude<F> {
    /// 2-D profile function: `f([x,z]) -> signed_distance`.
    pub profile: F,
    /// Half-extent along the Y axis.
    pub half_height: f64,
}
#[allow(dead_code)]
impl<F: Fn([f64; 2]) -> f64 + Send + Sync> SdfExtrude<F> {
    /// Create an extruded SDF.
    pub fn new(profile: F, half_height: f64) -> Self {
        Self {
            profile,
            half_height,
        }
    }
}
/// Revolve a 2-D profile curve (given as a function of r and y) around the Y axis.
///
/// The wrapped function takes `(r, y)` where `r = sqrt(x²+z²)` and returns
/// the signed distance in the 2-D profile plane.
pub struct SdfRevolution<F> {
    /// 2-D profile function: `f(r, y) -> signed_distance`.
    pub profile: F,
}
#[allow(dead_code)]
impl<F: Fn(f64, f64) -> f64 + Send + Sync> SdfRevolution<F> {
    /// Create a solid-of-revolution SDF.
    pub fn new(profile: F) -> Self {
        Self { profile }
    }
}
/// SDF union of two shapes.
#[derive(Debug, Clone)]
pub struct SdfUnion<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfUnion<A, B> {
    /// Create a union.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}
/// SDF for a finite cone with apex at origin pointing along +Y.
///
/// `half_angle` is the half-apex angle in radians; `height` is the cone height.
#[derive(Debug, Clone, Copy)]
pub struct SdfCone {
    /// Half-apex angle \[rad\].
    pub half_angle: f64,
    /// Cone height (distance from apex to base) \[m\].
    pub height: f64,
}
#[allow(dead_code)]
impl SdfCone {
    /// Create a cone SDF.
    pub fn new(half_angle: f64, height: f64) -> Self {
        Self { half_angle, height }
    }
}
/// Expand or contract a shape by a constant offset.
///
/// Positive `offset` grows the shape; negative shrinks it.
#[derive(Debug, Clone)]
pub struct SdfOffset<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Offset distance \[m\].
    pub offset: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfOffset<S> {
    /// Create an offset SDF.
    pub fn new(inner: S, offset: f64) -> Self {
        Self { inner, offset }
    }
}
/// Isosurface mesh extracted from a grid.
#[derive(Debug, Clone, Default)]
pub struct IsoMeshResult {
    /// Triangle list (each entry is one triangle).
    pub triangles: Vec<Triangle>,
    /// Total vertex count (3 × triangle count without sharing).
    pub vertex_count: usize,
}
/// SDF for a square-base pyramid with apex at `(0, height, 0)`.
///
/// Base is centred at origin in the XZ plane with half-side `half_base`.
#[derive(Debug, Clone, Copy)]
pub struct SdfPyramid {
    /// Half-side of the square base \[m\].
    pub half_base: f64,
    /// Height of the pyramid \[m\].
    pub height: f64,
}
#[allow(dead_code)]
impl SdfPyramid {
    /// Create a pyramid SDF.
    pub fn new(half_base: f64, height: f64) -> Self {
        Self { half_base, height }
    }
}
/// SDF difference: `A` minus `B`.
#[derive(Debug, Clone)]
pub struct SdfDifference<A, B> {
    /// Shape to subtract from.
    pub a: A,
    /// Shape to subtract.
    pub b: B,
}
#[allow(dead_code)]
impl<A: Sdf, B: Sdf> SdfDifference<A, B> {
    /// Create a difference.
    pub fn new(a: A, b: B) -> Self {
        Self { a, b }
    }
}
/// Displace an SDF surface with procedural value noise.
///
/// The displacement `amplitude` controls how far the surface is pushed;
/// `scale` controls the noise spatial frequency.
#[derive(Debug, Clone)]
pub struct SdfNoiseDisplace<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Noise spatial frequency scale.
    pub noise_scale: f64,
    /// Maximum displacement amplitude.
    pub amplitude: f64,
    /// Number of fBm octaves.
    pub octaves: u32,
}
#[allow(dead_code)]
impl<S: Sdf> SdfNoiseDisplace<S> {
    /// Create a noise-displaced SDF.
    pub fn new(inner: S, noise_scale: f64, amplitude: f64, octaves: u32) -> Self {
        Self {
            inner,
            noise_scale,
            amplitude,
            octaves,
        }
    }
}
/// SDF approximation for an axis-aligned ellipsoid with semi-axes `radii`.
///
/// Uses the Quilez approximation: exact on the boundary, approximate elsewhere.
#[derive(Debug, Clone, Copy)]
pub struct SdfEllipsoid {
    /// Semi-axis lengths (rx, ry, rz).
    pub radii: [f64; 3],
}
#[allow(dead_code)]
impl SdfEllipsoid {
    /// Create an ellipsoid SDF.
    pub fn new(rx: f64, ry: f64, rz: f64) -> Self {
        Self {
            radii: [rx, ry, rz],
        }
    }
}
/// A triangle in 3-D space (three vertices).
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// First vertex.
    pub v0: [f64; 3],
    /// Second vertex.
    pub v1: [f64; 3],
    /// Third vertex.
    pub v2: [f64; 3],
}
impl Triangle {
    /// Compute the triangle's outward normal (unnormalized).
    pub fn normal(&self) -> [f64; 3] {
        cross(sub(self.v1, self.v0), sub(self.v2, self.v0))
    }
}
/// Dynamic SDF expression tree node.
///
/// Allows building complex SDF scenes at runtime without generic type nesting.
pub enum SdfNode {
    /// A leaf: arbitrary boxed SDF.
    Leaf(Box<dyn Sdf>),
    /// A binary combination of two nodes.
    Binary {
        /// Binary operation.
        op: BoolOp,
        /// Left child.
        left: Box<SdfNode>,
        /// Right child.
        right: Box<SdfNode>,
        /// Blend radius (used for smooth ops).
        k: f64,
    },
}
#[allow(dead_code)]
impl SdfNode {
    /// Create a leaf node from any `Sdf` implementor.
    pub fn leaf<S: Sdf + 'static>(s: S) -> Self {
        Self::Leaf(Box::new(s))
    }
    /// Combine two nodes with a boolean operation.
    pub fn combine(op: BoolOp, left: Self, right: Self, k: f64) -> Self {
        Self::Binary {
            op,
            left: Box::new(left),
            right: Box::new(right),
            k,
        }
    }
    /// Evaluate the SDF tree at point `p`.
    pub fn eval(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Leaf(s) => s.dist(p),
            Self::Binary { op, left, right, k } => {
                let a = left.eval(p);
                let b = right.eval(p);
                match op {
                    BoolOp::Union => a.min(b),
                    BoolOp::Intersection => a.max(b),
                    BoolOp::Difference => a.max(-b),
                    BoolOp::SmoothUnion => sdf_smooth_union(a, b, *k),
                    BoolOp::SmoothIntersection => sdf_smooth_intersection(a, b, *k),
                    BoolOp::SmoothDifference => sdf_smooth_difference(a, b, *k),
                }
            }
        }
    }
}
/// SDF for a triangular prism (equilateral cross-section) aligned with Y axis.
///
/// `side` is the equilateral triangle side length; `half_height` is along Y.
#[derive(Debug, Clone, Copy)]
pub struct SdfTriangularPrism {
    /// Triangle side length \[m\].
    pub side: f64,
    /// Half-height along Y \[m\].
    pub half_height: f64,
}
#[allow(dead_code)]
impl SdfTriangularPrism {
    /// Create a triangular prism SDF.
    pub fn new(side: f64, half_height: f64) -> Self {
        Self { side, half_height }
    }
}
/// SDF approximation for the Schoen gyroid minimal surface.
///
/// The gyroid is an infinitely periodic surface: `sin x cos y + sin y cos z + sin z cos x = 0`.
/// This is an approximation; the exact gyroid has no closed-form SDF.
#[derive(Debug, Clone, Copy)]
pub struct SdfGyroid {
    /// Spatial frequency scale.
    pub scale: f64,
    /// Thickness of the gyroid sheet.
    pub thickness: f64,
}
#[allow(dead_code)]
impl SdfGyroid {
    /// Create a gyroid SDF.
    pub fn new(scale: f64, thickness: f64) -> Self {
        Self { scale, thickness }
    }
}
/// Hollow (shell) version of an SDF: the region between two offset surfaces.
///
/// Produces a thin shell of thickness `thickness` around the surface.
#[derive(Debug, Clone)]
pub struct SdfShell<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Shell thickness \[m\] (half on each side of the zero isosurface).
    pub thickness: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfShell<S> {
    /// Create a shell SDF of given `thickness`.
    pub fn new(inner: S, thickness: f64) -> Self {
        Self { inner, thickness }
    }
}
/// Apply a bend deformation around the Y axis to an SDF.
///
/// `strength` controls how much bending occurs per unit of X.
#[derive(Debug, Clone)]
pub struct SdfBend<S> {
    /// Wrapped SDF.
    pub inner: S,
    /// Bend rate \[rad/m\].
    pub strength: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfBend<S> {
    /// Create a bent SDF.
    pub fn new(inner: S, strength: f64) -> Self {
        Self { inner, strength }
    }
}
/// SDF for a line-segment tube: a capsule with flat ends (cylinder).
///
/// Computes distance to the axis segment `[a, b]` minus `radius`.
#[derive(Debug, Clone, Copy)]
pub struct SdfLineSegment {
    /// Start point of the segment.
    pub a: [f64; 3],
    /// End point of the segment.
    pub b: [f64; 3],
    /// Tube radius \[m\].
    pub radius: f64,
}
#[allow(dead_code)]
impl SdfLineSegment {
    /// Create a line-segment SDF.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }
}
/// Result of a ray-march query.
#[derive(Debug, Clone, Copy)]
pub struct RayMarchResult {
    /// Whether the ray hit the surface.
    pub hit: bool,
    /// Distance along the ray to the hit point (if `hit` is true).
    pub t: f64,
    /// Hit point in world space.
    pub point: [f64; 3],
    /// Number of steps taken.
    pub steps: u32,
}
/// Physics collision proxy that first checks a bounding sphere, then the detailed SDF.
///
/// Avoids expensive SDF evaluation for clearly non-colliding points.
#[derive(Debug, Clone)]
pub struct SdfBoundedProxy<S> {
    /// The detailed SDF.
    pub inner: S,
    /// Centre of the bounding sphere.
    pub bsphere_center: [f64; 3],
    /// Radius of the bounding sphere.
    pub bsphere_radius: f64,
}
#[allow(dead_code)]
impl<S: Sdf> SdfBoundedProxy<S> {
    /// Create a bounded proxy.
    pub fn new(inner: S, bsphere_center: [f64; 3], bsphere_radius: f64) -> Self {
        Self {
            inner,
            bsphere_center,
            bsphere_radius,
        }
    }
    /// Returns true if point `p` is definitely outside the bounding sphere.
    pub fn outside_bounds(&self, p: [f64; 3]) -> bool {
        len(sub(p, self.bsphere_center)) > self.bsphere_radius
    }
}
/// SDF for an axis-aligned box centred at the origin.
///
/// This is a standalone, self-contained version with `half_extents` per axis.
#[derive(Debug, Clone, Copy)]
pub struct SdfBox {
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}
#[allow(dead_code)]
impl SdfBox {
    /// Create a box SDF.
    pub fn new(hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            half_extents: [hx, hy, hz],
        }
    }
}
/// Result of an SDF collision query.
#[derive(Debug, Clone, Copy)]
pub struct SdfCollisionResult {
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Contact normal pointing from B toward A (unit vector).
    pub normal: [f64; 3],
    /// Contact point (on or near surface).
    pub contact_point: [f64; 3],
}
/// SDF for a sphere of given `radius` centred at the origin.
#[derive(Debug, Clone, Copy)]
pub struct SdfSphere {
    /// Sphere radius \[m\].
    pub radius: f64,
}
#[allow(dead_code)]
impl SdfSphere {
    /// Create a new sphere SDF with `radius`.
    pub fn new(radius: f64) -> Self {
        Self { radius }
    }
}
