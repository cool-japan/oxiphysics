//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::linear_algebra::{Mat4, Real, Vec3, Vec4};

/// A ray defined by an origin and direction.
#[derive(Debug, Clone)]
pub struct Ray {
    /// Origin of the ray.
    pub origin: Vec3,
    /// Direction of the ray (should be normalized).
    pub direction: Vec3,
}
impl Ray {
    /// Create a new ray.
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        Self { origin, direction }
    }
    /// Get the point along the ray at parameter `t`.
    pub fn point_at(&self, t: Real) -> Vec3 {
        self.origin + self.direction * t
    }
    /// Intersect this ray with a plane, returning the parameter `t` if they intersect.
    ///
    /// Returns `None` if the ray is parallel to (or lies in) the plane.
    pub fn intersect_plane(&self, plane: &Plane) -> Option<Real> {
        let denom = plane.normal.dot(&self.direction);
        if denom.abs() < 1e-10 {
            return None;
        }
        let t = (plane.distance - plane.normal.dot(&self.origin)) / denom;
        Some(t)
    }
}
/// A view frustum defined by 6 planes (near, far, left, right, top, bottom).
#[derive(Debug, Clone)]
pub struct Frustum {
    /// The six bounding planes of the frustum.
    pub planes: [Plane; 6],
}
impl Frustum {
    /// Extract frustum planes from a combined view-projection matrix.
    ///
    /// The resulting planes have inward-facing normals.
    pub fn from_view_projection(vp: &Mat4) -> Self {
        let row = |i: usize| -> Vec4 { Vec4::new(vp[(i, 0)], vp[(i, 1)], vp[(i, 2)], vp[(i, 3)]) };
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);
        let extract = |v: Vec4| -> Plane {
            let n = Vec3::new(v.x, v.y, v.z);
            let len = n.norm();
            if len < 1e-10 {
                Plane::new(Vec3::new(1.0, 0.0, 0.0), 0.0)
            } else {
                Plane::new(n / len, v.w / len)
            }
        };
        Frustum {
            planes: [
                extract(r3 + r2),
                extract(r3 - r2),
                extract(r3 + r0),
                extract(r3 - r0),
                extract(r3 + r1),
                extract(r3 - r1),
            ],
        }
    }
    /// Test whether a point lies inside (or on) the frustum.
    pub fn contains_point(&self, point: &Vec3) -> bool {
        for plane in &self.planes {
            if plane.signed_distance(point) < -1e-10 {
                return false;
            }
        }
        true
    }
    /// Test whether a sphere intersects the frustum.
    pub fn intersects_sphere(&self, center: &Vec3, radius: Real) -> bool {
        for plane in &self.planes {
            if plane.signed_distance(center) < -radius {
                return false;
            }
        }
        true
    }
}
impl Frustum {
    /// Test whether a sphere is fully contained inside the frustum.
    ///
    /// Returns `true` if every plane's signed distance to `center` is ≥ `radius`.
    #[allow(dead_code)]
    pub fn contains_sphere(&self, center: &Vec3, radius: Real) -> bool {
        self.planes
            .iter()
            .all(|p| p.signed_distance(center) >= radius)
    }
    /// Test whether an axis-aligned bounding box (AABB) intersects the frustum.
    ///
    /// Uses the separating-axis theorem: for each frustum plane, if all 8
    /// corners of the AABB are on the outside half-space, the AABB is
    /// rejected.
    #[allow(dead_code)]
    pub fn intersects_aabb(&self, aabb: &Aabb) -> bool {
        for plane in &self.planes {
            let p = Vec3::new(
                if plane.normal.x >= 0.0 {
                    aabb.max[0]
                } else {
                    aabb.min[0]
                },
                if plane.normal.y >= 0.0 {
                    aabb.max[1]
                } else {
                    aabb.min[1]
                },
                if plane.normal.z >= 0.0 {
                    aabb.max[2]
                } else {
                    aabb.min[2]
                },
            );
            if plane.signed_distance(&p) < 0.0 {
                return false;
            }
        }
        true
    }
    /// Build a frustum directly from a combined view-projection matrix.
    ///
    /// Alias for [`Frustum::from_view_projection`] with a more explicit name.
    #[allow(dead_code)]
    pub fn extract_from_view_proj(vp: &Mat4) -> Self {
        Self::from_view_projection(vp)
    }
}
/// A plane defined by a normal and signed distance from origin.
#[derive(Debug, Clone)]
pub struct Plane {
    /// Outward-facing normal (should be normalized).
    pub normal: Vec3,
    /// Signed distance from the origin along the normal.
    pub distance: Real,
}
impl Plane {
    /// Create a new plane.
    pub fn new(normal: Vec3, distance: Real) -> Self {
        Self { normal, distance }
    }
    /// Signed distance from a point to this plane.
    pub fn signed_distance(&self, point: &Vec3) -> Real {
        self.normal.dot(point) - self.distance
    }
    /// Create a plane from three non-collinear points.
    ///
    /// The normal is computed as `(b - a).cross(c - a)` normalized.
    /// Returns `None` if points are collinear.
    pub fn from_points(a: &Vec3, b: &Vec3, c: &Vec3) -> Option<Self> {
        let ab = b - a;
        let ac = c - a;
        let n = ab.cross(&ac);
        let len = n.norm();
        if len < 1e-10 {
            return None;
        }
        let normal = n / len;
        let distance = normal.dot(a);
        Some(Self { normal, distance })
    }
    /// Project a point onto this plane.
    pub fn project_point(&self, point: &Vec3) -> Vec3 {
        let d = self.signed_distance(point);
        point - self.normal * d
    }
    /// Reflect a point across this plane.
    pub fn reflect_point(&self, point: &Vec3) -> Vec3 {
        let d = self.signed_distance(point);
        point - self.normal * (2.0 * d)
    }
    /// Classify a point relative to this plane.
    ///
    /// Returns positive if the point is on the front side, negative on the
    /// back side, or zero if it lies on the plane (within tolerance).
    pub fn classify_point(&self, point: &Vec3, tolerance: Real) -> i32 {
        let d = self.signed_distance(point);
        if d > tolerance {
            1
        } else if d < -tolerance {
            -1
        } else {
            0
        }
    }
}
impl Plane {
    /// Ray–plane intersection returning the hit `t` parameter.
    ///
    /// Identical to [`Ray::intersect_plane`] but callable directly on a `Plane`.
    #[allow(dead_code)]
    pub fn intersect_ray(&self, ray: &Ray) -> Option<Real> {
        ray.intersect_plane(self)
    }
}
/// Axis-aligned bounding box in 3-D.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
#[allow(dead_code)]
impl Aabb {
    /// Create an AABB from its min and max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Create an AABB from a slice of points (returns `None` if empty).
    pub fn from_points(pts: &[[f64; 3]]) -> Option<Self> {
        if pts.is_empty() {
            return None;
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for &p in &pts[1..] {
            for i in 0..3 {
                if p[i] < mn[i] {
                    mn[i] = p[i];
                }
                if p[i] > mx[i] {
                    mx[i] = p[i];
                }
            }
        }
        Some(Self { min: mn, max: mx })
    }
    /// Test whether a point is inside (or on the boundary of) this AABB.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] <= self.max[i])
    }
    /// Smallest AABB that contains both `self` and `other`.
    pub fn union(&self, other: &Aabb) -> Aabb {
        let mut mn = self.min;
        let mut mx = self.max;
        for i in 0..3 {
            if other.min[i] < mn[i] {
                mn[i] = other.min[i];
            }
            if other.max[i] > mx[i] {
                mx[i] = other.max[i];
            }
        }
        Aabb { min: mn, max: mx }
    }
    /// Intersection of two AABBs.  Returns `None` if they do not overlap.
    pub fn intersect(&self, other: &Aabb) -> Option<Aabb> {
        let mut mn = [0.0_f64; 3];
        let mut mx = [0.0_f64; 3];
        for i in 0..3 {
            mn[i] = self.min[i].max(other.min[i]);
            mx[i] = self.max[i].min(other.max[i]);
            if mn[i] > mx[i] {
                return None;
            }
        }
        Some(Aabb { min: mn, max: mx })
    }
    /// Expand (or contract for negative `amount`) all sides by `amount`.
    pub fn expand(&self, amount: f64) -> Aabb {
        Aabb {
            min: [
                self.min[0] - amount,
                self.min[1] - amount,
                self.min[2] - amount,
            ],
            max: [
                self.max[0] + amount,
                self.max[1] + amount,
                self.max[2] + amount,
            ],
        }
    }
    /// Surface area of this AABB.
    pub fn surface_area(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }
    /// Volume of this AABB.
    pub fn volume(&self) -> f64 {
        (self.max[0] - self.min[0]) * (self.max[1] - self.min[1]) * (self.max[2] - self.min[2])
    }
    /// Centre of this AABB.
    pub fn centre(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Half-extents of this AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }
    /// Ray–AABB intersection test (slab method).
    ///
    /// `origin` and `dir` define the ray; `dir` need not be normalised.
    /// Returns `Some((t_enter, t_exit))` if the ray hits, otherwise `None`.
    /// Only positive-t intersections (in front of the origin) are reported.
    pub fn intersect_ray(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<(f64, f64)> {
        let mut t_min = f64::NEG_INFINITY;
        let mut t_max = f64::INFINITY;
        for i in 0..3 {
            if dir[i].abs() < 1e-300 {
                if origin[i] < self.min[i] || origin[i] > self.max[i] {
                    return None;
                }
            } else {
                let inv_d = 1.0 / dir[i];
                let t1 = (self.min[i] - origin[i]) * inv_d;
                let t2 = (self.max[i] - origin[i]) * inv_d;
                let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
                t_min = t_min.max(ta);
                t_max = t_max.min(tb);
                if t_min > t_max {
                    return None;
                }
            }
        }
        if t_max < 0.0 {
            return None;
        }
        Some((t_min, t_max))
    }
    /// Test whether another AABB overlaps this one.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.intersect(other).is_some()
    }
}
impl Aabb {
    /// Merge two AABBs into the smallest AABB containing both.
    ///
    /// Alias for [`Aabb::union`].
    #[allow(dead_code)]
    pub fn merge(&self, other: &Aabb) -> Aabb {
        self.union(other)
    }
    /// Expand this AABB by `amount` in every direction.
    ///
    /// Alias for [`Aabb::expand`] following the task spec naming.
    #[allow(dead_code)]
    pub fn expand_by(&self, amount: f64) -> Aabb {
        self.expand(amount)
    }
    /// Return the closest point on (or inside) this AABB to `p`.
    #[allow(dead_code)]
    pub fn closest_point(&self, p: [f64; 3]) -> [f64; 3] {
        [
            p[0].clamp(self.min[0], self.max[0]),
            p[1].clamp(self.min[1], self.max[1]),
            p[2].clamp(self.min[2], self.max[2]),
        ]
    }
    /// Test whether `p` is strictly inside (not on the boundary) of this AABB.
    #[allow(dead_code)]
    pub fn contains_point_strict(&self, p: [f64; 3]) -> bool {
        (0..3).all(|i| p[i] > self.min[i] && p[i] < self.max[i])
    }
}
/// A dual number `a + b * ε` where `ε² = 0`.
///
/// Dual numbers support forward-mode automatic differentiation: evaluating
/// `f(a + 1*ε)` gives `f(a) + f'(a)*ε`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// Real part.
    pub re: Real,
    /// Dual (infinitesimal) part.
    pub du: Real,
}
impl Dual {
    /// Construct a dual number.
    #[allow(dead_code)]
    pub fn new(re: Real, du: Real) -> Self {
        Self { re, du }
    }
    /// Constant (derivative = 0).
    #[allow(dead_code)]
    pub fn constant(re: Real) -> Self {
        Self { re, du: 0.0 }
    }
    /// Variable (derivative = 1).
    #[allow(dead_code)]
    pub fn variable(re: Real) -> Self {
        Self { re, du: 1.0 }
    }
    /// Square root: √(a + b*ε) = √a + b/(2√a) * ε.
    #[allow(dead_code)]
    pub fn sqrt(self) -> Self {
        let r = self.re.sqrt();
        Self {
            re: r,
            du: self.du / (2.0 * r),
        }
    }
    /// Sine: sin(a + b*ε) = sin(a) + b*cos(a)*ε.
    #[allow(dead_code)]
    pub fn sin(self) -> Self {
        Self {
            re: self.re.sin(),
            du: self.du * self.re.cos(),
        }
    }
    /// Cosine: cos(a + b*ε) = cos(a) - b*sin(a)*ε.
    #[allow(dead_code)]
    pub fn cos(self) -> Self {
        Self {
            re: self.re.cos(),
            du: -self.du * self.re.sin(),
        }
    }
    /// Exponential: exp(a + b*ε) = exp(a) + b*exp(a)*ε.
    #[allow(dead_code)]
    pub fn exp(self) -> Self {
        let e = self.re.exp();
        Self {
            re: e,
            du: self.du * e,
        }
    }
    /// Natural log: ln(a + b*ε) = ln(a) + b/a * ε.
    #[allow(dead_code)]
    pub fn ln(self) -> Self {
        Self {
            re: self.re.ln(),
            du: self.du / self.re,
        }
    }
    /// Raise to integer power: (a+b*ε)^n = a^n + n*b*a^{n-1}*ε.
    #[allow(dead_code)]
    pub fn powi(self, n: i32) -> Self {
        let re = self.re.powi(n);
        let du = n as Real * self.re.powi(n - 1) * self.du;
        Self { re, du }
    }
}
