// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Softbody collision detection and response.
//!
//! Provides:
//! - AABB broadphase for particle groups
//! - BVH-accelerated triangle-mesh collision
//! - Point-triangle and edge-edge closest-point tests
//! - Penalty and impulse-based collision response
//! - Continuous collision detection (CCD) for particles vs. triangles
//! - Self-collision for cloth meshes
//! - Friction response
//! - Collision group bitmask filtering

use oxiphysics_core::collision::{AabbRaw, Bvh, BvhLeaf};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_len_sq(v: [f64; 3]) -> f64 {
    vec3_dot(v, v)
}

#[inline]
fn vec3_len(v: [f64; 3]) -> f64 {
    vec3_len_sq(v).sqrt()
}

#[inline]
fn vec3_normalize(v: [f64; 3]) -> Option<[f64; 3]> {
    let len = vec3_len(v);
    if len < 1e-14 {
        None
    } else {
        Some(vec3_scale(v, 1.0 / len))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SoftBodyAabb
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box (AABB) for a group of particles.
///
/// Used as a broadphase acceleration structure. Call [`SoftBodyAabb::refit`]
/// each frame to update the box from current particle positions.
#[derive(Debug, Clone, Copy)]
pub struct SoftBodyAabb {
    /// Minimum corner of the AABB.
    pub min: [f64; 3],
    /// Maximum corner of the AABB.
    pub max: [f64; 3],
}

impl SoftBodyAabb {
    /// Creates an AABB from explicit min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Creates an empty (inverted) AABB ready to be expanded.
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    /// Refits the AABB to tightly enclose all given particle positions.
    pub fn refit(&mut self, positions: &[[f64; 3]]) {
        *self = Self::empty();
        for &p in positions {
            for (i, &pi) in p.iter().enumerate() {
                self.min[i] = self.min[i].min(pi);
                self.max[i] = self.max[i].max(pi);
            }
        }
    }

    /// Returns `true` if this AABB overlaps with `other`.
    pub fn overlaps(&self, other: &SoftBodyAabb) -> bool {
        for i in 0..3 {
            if self.min[i] > other.max[i] || self.max[i] < other.min[i] {
                return false;
            }
        }
        true
    }

    /// Returns `true` if `point` lies inside (or on the surface of) this AABB.
    pub fn contains_point(&self, point: [f64; 3]) -> bool {
        (0..3).all(|i| point[i] >= self.min[i] && point[i] <= self.max[i])
    }

    /// Expands the AABB by `margin` in all directions.
    pub fn expand(&self, margin: f64) -> Self {
        Self {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Half-extents of the AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TriangleMeshCollider
// ─────────────────────────────────────────────────────────────────────────────

/// A triangle mesh collider with a simple flat list of per-triangle AABBs.
///
/// Each triangle is pre-processed into a leaf AABB.  Broadphase queries
/// eliminate non-overlapping triangles before the exact test.
#[derive(Debug, Clone)]
pub struct TriangleMeshCollider {
    /// Triangle vertices: each entry is three vertex positions.
    pub triangles: Vec<[[f64; 3]; 3]>,
    /// Per-triangle AABBs (recomputed from `triangles`).
    leaf_aabbs: Vec<SoftBodyAabb>,
}

impl TriangleMeshCollider {
    /// Creates a new collider from a list of triangles.
    ///
    /// Each element of `triangles` is `[v0, v1, v2]`.
    pub fn new(triangles: Vec<[[f64; 3]; 3]>) -> Self {
        let leaf_aabbs = triangles
            .iter()
            .map(|tri| {
                let mut aabb = SoftBodyAabb::empty();
                aabb.refit(tri);
                aabb
            })
            .collect();
        Self {
            triangles,
            leaf_aabbs,
        }
    }

    /// Returns the number of triangles in the collider.
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Returns the indices of triangles whose AABB overlaps `query_aabb`.
    pub fn broad_query(&self, query_aabb: &SoftBodyAabb) -> Vec<usize> {
        self.leaf_aabbs
            .iter()
            .enumerate()
            .filter(|(_, a)| a.overlaps(query_aabb))
            .map(|(i, _)| i)
            .collect()
    }

    /// Returns the world-space AABB enclosing all triangles.
    pub fn total_aabb(&self) -> SoftBodyAabb {
        let mut aabb = SoftBodyAabb::empty();
        for tri in &self.triangles {
            aabb.refit(tri);
        }
        aabb
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PointTriangleContact
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a point-triangle closest-point query.
#[derive(Debug, Clone)]
pub struct PointTriangleContact {
    /// Closest point on the triangle to the query point.
    pub closest_point: [f64; 3],
    /// Barycentric coordinates (u, v, w) of the closest point.
    pub barycentric: [f64; 3],
    /// Squared distance from the query point to the closest point.
    pub dist_sq: f64,
    /// Triangle normal (outward).
    pub normal: [f64; 3],
}

impl PointTriangleContact {
    /// Computes the closest point on triangle `(a, b, c)` to point `p`.
    ///
    /// Based on the Ericson "Real-Time Collision Detection" algorithm.
    pub fn compute(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        let ab = vec3_sub(b, a);
        let ac = vec3_sub(c, a);
        let ap = vec3_sub(p, a);

        let d1 = vec3_dot(ab, ap);
        let d2 = vec3_dot(ac, ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            let bary = [1.0, 0.0, 0.0];
            let closest = a;
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let bp = vec3_sub(p, b);
        let d3 = vec3_dot(ab, bp);
        let d4 = vec3_dot(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            let bary = [0.0, 1.0, 0.0];
            let closest = b;
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            let bary = [1.0 - v, v, 0.0];
            let closest = vec3_add(a, vec3_scale(ab, v));
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let cp = vec3_sub(p, c);
        let d5 = vec3_dot(ab, cp);
        let d6 = vec3_dot(ac, cp);
        if d6 >= 0.0 && d5 <= d6 {
            let bary = [0.0, 0.0, 1.0];
            let closest = c;
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            let bary = [1.0 - w, 0.0, w];
            let closest = vec3_add(a, vec3_scale(ac, w));
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            let bary = [0.0, 1.0 - w, w];
            let closest = vec3_add(b, vec3_scale(vec3_sub(c, b), w));
            let diff = vec3_sub(p, closest);
            let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
            return Self {
                closest_point: closest,
                barycentric: bary,
                dist_sq: vec3_len_sq(diff),
                normal,
            };
        }

        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        let u = 1.0 - v - w;
        let bary = [u, v, w];
        let closest = [
            a[0] + ab[0] * v + ac[0] * w,
            a[1] + ab[1] * v + ac[1] * w,
            a[2] + ab[2] * v + ac[2] * w,
        ];
        let diff = vec3_sub(p, closest);
        let normal = vec3_normalize(vec3_cross(ab, ac)).unwrap_or([0.0, 1.0, 0.0]);
        Self {
            closest_point: closest,
            barycentric: bary,
            dist_sq: vec3_len_sq(diff),
            normal,
        }
    }

    /// Signed penetration depth (negative = penetrating from front face).
    ///
    /// Positive means the point is above the triangle plane.
    pub fn signed_depth(&self, point: [f64; 3]) -> f64 {
        let to_point = vec3_sub(point, self.closest_point);
        vec3_dot(to_point, self.normal)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EdgeEdgeContact
// ─────────────────────────────────────────────────────────────────────────────

/// Closest-point result between two line segments (edges).
#[derive(Debug, Clone)]
pub struct EdgeEdgeContact {
    /// Parameter t ∈ \[0,1\] along edge 1 for the closest point.
    pub t1: f64,
    /// Parameter t ∈ \[0,1\] along edge 2 for the closest point.
    pub t2: f64,
    /// Closest point on edge 1.
    pub point1: [f64; 3],
    /// Closest point on edge 2.
    pub point2: [f64; 3],
    /// Squared distance between the two closest points.
    pub dist_sq: f64,
}

impl EdgeEdgeContact {
    /// Computes the closest points between segment `(p1, p2)` and `(p3, p4)`.
    pub fn compute(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> Self {
        let d1 = vec3_sub(p2, p1);
        let d2 = vec3_sub(p4, p3);
        let r = vec3_sub(p1, p3);

        let a = vec3_dot(d1, d1);
        let e = vec3_dot(d2, d2);
        let f = vec3_dot(d2, r);

        let (t1, t2) = if a <= 1e-14 && e <= 1e-14 {
            (0.0, 0.0)
        } else if a <= 1e-14 {
            (0.0, (f / e).clamp(0.0, 1.0))
        } else {
            let c = vec3_dot(d1, r);
            if e <= 1e-14 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else {
                let b = vec3_dot(d1, d2);
                let denom = a * e - b * b;
                let s = if denom.abs() > 1e-14 {
                    ((b * f - c * e) / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let t = (b * s + f) / e;
                let (t_clamped, s_adjusted) = if t < 0.0 {
                    (0.0, (-c / a).clamp(0.0, 1.0))
                } else if t > 1.0 {
                    (1.0, ((b - c) / a).clamp(0.0, 1.0))
                } else {
                    (t, s)
                };
                (s_adjusted, t_clamped)
            }
        };

        let point1 = vec3_add(p1, vec3_scale(d1, t1));
        let point2 = vec3_add(p3, vec3_scale(d2, t2));
        let diff = vec3_sub(point1, point2);

        Self {
            t1,
            t2,
            point1,
            point2,
            dist_sq: vec3_len_sq(diff),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SoftBodyContact
// ─────────────────────────────────────────────────────────────────────────────

/// A contact between a softbody particle and a surface.
#[derive(Debug, Clone)]
pub struct SoftBodyContact {
    /// World-space contact point (on the collider surface).
    pub contact_point: [f64; 3],
    /// Contact normal pointing from collider to particle (outward).
    pub normal: [f64; 3],
    /// Penetration depth (positive = penetrating).
    pub penetration_depth: f64,
    /// Index of the penetrating particle in the soft body.
    pub particle_index: usize,
    /// Triangle index in the collider mesh that was hit.
    pub triangle_index: usize,
}

impl SoftBodyContact {
    /// Creates a new soft-body contact.
    pub fn new(
        contact_point: [f64; 3],
        normal: [f64; 3],
        penetration_depth: f64,
        particle_index: usize,
        triangle_index: usize,
    ) -> Self {
        Self {
            contact_point,
            normal,
            penetration_depth,
            particle_index,
            triangle_index,
        }
    }

    /// Returns `true` if the contact is actually penetrating (depth > threshold).
    pub fn is_penetrating(&self, threshold: f64) -> bool {
        self.penetration_depth > threshold
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CollisionResponse
// ─────────────────────────────────────────────────────────────────────────────

/// Collision response mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResponseMode {
    /// Penalty spring: apply a restoring force proportional to penetration.
    Penalty,
    /// Impulse: apply a velocity impulse to resolve penetration.
    Impulse,
}

/// Resolves penetration contacts between softbody particles and a mesh collider.
#[derive(Debug, Clone)]
pub struct CollisionResponse {
    /// Response mode.
    pub mode: ResponseMode,
    /// Penalty spring stiffness \[Pa/m\] (used in `Penalty` mode).
    pub penalty_stiffness: f64,
    /// Coefficient of restitution (0 = inelastic, 1 = elastic).
    pub restitution: f64,
}

impl CollisionResponse {
    /// Creates a penalty-based collision response.
    pub fn penalty(stiffness: f64) -> Self {
        Self {
            mode: ResponseMode::Penalty,
            penalty_stiffness: stiffness,
            restitution: 0.0,
        }
    }

    /// Creates an impulse-based collision response.
    pub fn impulse(restitution: f64) -> Self {
        Self {
            mode: ResponseMode::Impulse,
            penalty_stiffness: 0.0,
            restitution: restitution.clamp(0.0, 1.0),
        }
    }

    /// Computes the penalty force \[N\] for a given contact.
    ///
    /// Returns a force vector to apply to the particle (pushing it out).
    pub fn penalty_force(&self, contact: &SoftBodyContact) -> [f64; 3] {
        let mag = self.penalty_stiffness * contact.penetration_depth;
        vec3_scale(contact.normal, mag)
    }

    /// Computes the position correction for a penetrating particle.
    ///
    /// Simply pushes the particle out of the surface by `penetration_depth`.
    pub fn position_correction(&self, contact: &SoftBodyContact) -> [f64; 3] {
        vec3_scale(contact.normal, contact.penetration_depth)
    }

    /// Resolves a contact between softbody particle and triangle mesh.
    ///
    /// `position` and `velocity` are particle state (in/out).
    /// Returns the updated position and velocity.
    pub fn resolve(
        &self,
        position: [f64; 3],
        velocity: [f64; 3],
        contact: &SoftBodyContact,
    ) -> ([f64; 3], [f64; 3]) {
        if !contact.is_penetrating(0.0) {
            return (position, velocity);
        }
        let new_pos = vec3_add(position, self.position_correction(contact));
        let new_vel = match self.mode {
            ResponseMode::Penalty => velocity,
            ResponseMode::Impulse => {
                let vn = vec3_dot(velocity, contact.normal);
                if vn < 0.0 {
                    // Reflect normal component.
                    let impulse = -(1.0 + self.restitution) * vn;
                    vec3_add(velocity, vec3_scale(contact.normal, impulse))
                } else {
                    velocity
                }
            }
        };
        (new_pos, new_vel)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SelfCollision
// ─────────────────────────────────────────────────────────────────────────────

/// Cloth self-collision: particle vs. triangle within the same mesh.
///
/// Checks each particle against the triangles of the mesh that it is not a
/// corner of and returns a contact for every pair closer than the cloth
/// thickness.  A triangle is *adjacent* to a particle when it contains that
/// particle as a corner; adjacent pairs never produce contacts.
///
/// [`SelfCollision::detect`] builds a bounding-volume hierarchy over the
/// triangles (each box inflated by the thickness) and only runs the exact
/// point-triangle distance test on the triangles whose box contains the
/// particle, so the cost is roughly `O((P + T) log T)` instead of `O(P * T)`.
#[derive(Debug, Clone)]
pub struct SelfCollision {
    /// Cloth thickness (two-sided guard band) \[m\].
    pub thickness: f64,
}

impl SelfCollision {
    /// Creates a self-collision detector with the given cloth thickness \[m\].
    pub fn with_thickness(thickness: f64) -> Self {
        Self { thickness }
    }

    /// Creates a new self-collision detector.
    ///
    /// The mesh arguments are ignored: adjacency is derived from the triangle
    /// list that is passed to [`SelfCollision::detect`].
    #[deprecated(
        since = "0.1.3",
        note = "adjacency is derived from the triangles passed to `detect`; use `SelfCollision::with_thickness`"
    )]
    pub fn new(_num_particles: usize, _triangles: &[[usize; 3]], thickness: f64) -> Self {
        Self::with_thickness(thickness)
    }

    /// Positions of the three corners of `tri`, or `None` when a corner index
    /// is out of range or a corner position is not finite.
    fn triangle_corners(positions: &[[f64; 3]], tri: &[usize; 3]) -> Option<[[f64; 3]; 3]> {
        let corners = [
            *positions.get(tri[0])?,
            *positions.get(tri[1])?,
            *positions.get(tri[2])?,
        ];
        corners
            .iter()
            .flatten()
            .all(|coordinate| coordinate.is_finite())
            .then_some(corners)
    }

    /// Detects self-collision contacts.
    ///
    /// * `positions` – current particle positions.
    /// * `triangles` – mesh connectivity (indices into `positions`).
    ///
    /// Returns a `SoftBodyContact` for every particle that is closer than
    /// [`SelfCollision::thickness`] to a triangle it is not a corner of.  The
    /// contacts are ordered by particle index, then by triangle index.
    /// Triangles that reference a particle outside `positions`, or whose
    /// corners are not finite, are ignored; a non-positive or non-finite
    /// thickness yields no contacts.
    pub fn detect(&self, positions: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<SoftBodyContact> {
        let mut contacts = Vec::new();
        let thickness = self.thickness;
        if !(thickness.is_finite() && thickness > 0.0) {
            return contacts;
        }
        let mut leaves = Vec::with_capacity(triangles.len());
        for (ti, tri) in triangles.iter().enumerate() {
            let Some([a, b, c]) = Self::triangle_corners(positions, tri) else {
                continue;
            };
            let mut min = a;
            let mut max = a;
            for corner in [b, c] {
                for axis in 0..3 {
                    min[axis] = min[axis].min(corner[axis]);
                    max[axis] = max[axis].max(corner[axis]);
                }
            }
            for axis in 0..3 {
                min[axis] -= thickness;
                max[axis] += thickness;
            }
            leaves.push(BvhLeaf {
                id: ti as u64,
                aabb: AabbRaw::new(min, max),
            });
        }
        let bvh = Bvh::build(leaves);
        let thickness_sq = thickness * thickness;
        for (pi, &pos) in positions.iter().enumerate() {
            let mut hits = bvh.query(&AabbRaw::new(pos, pos));
            hits.sort_unstable();
            for id in hits {
                let ti = id as usize;
                let Some(tri) = triangles.get(ti) else {
                    continue;
                };
                if tri.contains(&pi) {
                    continue;
                }
                let Some([a, b, c]) = Self::triangle_corners(positions, tri) else {
                    continue;
                };
                let result = PointTriangleContact::compute(pos, a, b, c);
                if result.dist_sq < thickness_sq {
                    let depth = thickness - result.dist_sq.sqrt();
                    contacts.push(SoftBodyContact::new(
                        result.closest_point,
                        result.normal,
                        depth,
                        pi,
                        ti,
                    ));
                }
            }
        }
        contacts
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FrictionResponse
// ─────────────────────────────────────────────────────────────────────────────

/// Tangential friction impulse at a soft-body contact.
///
/// Implements Coulomb friction: the tangential velocity is clamped to
/// `μ * |v_n|`.
#[derive(Debug, Clone)]
pub struct FrictionResponse {
    /// Static friction coefficient.
    pub mu_static: f64,
    /// Kinetic friction coefficient.
    pub mu_kinetic: f64,
}

impl FrictionResponse {
    /// Creates a new friction response.
    pub fn new(mu_static: f64, mu_kinetic: f64) -> Self {
        Self {
            mu_static,
            mu_kinetic,
        }
    }

    /// Applies Coulomb friction to a velocity vector.
    ///
    /// * `velocity` – particle velocity before friction.
    /// * `normal` – contact normal (unit vector, pointing away from surface).
    ///
    /// Returns the velocity after friction is applied.
    pub fn apply(&self, velocity: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
        let vn_mag = vec3_dot(velocity, normal);
        let v_normal = vec3_scale(normal, vn_mag);
        let v_tangent = vec3_sub(velocity, v_normal);
        let vt_len = vec3_len(v_tangent);
        if vt_len < 1e-14 {
            return velocity;
        }
        // Friction force magnitude = μ * |v_normal|
        let friction_limit = self.mu_kinetic * vn_mag.abs();
        if vt_len <= friction_limit {
            // Static: zero tangential velocity.
            v_normal
        } else {
            // Kinetic: reduce tangential component.
            let scale = (vt_len - friction_limit) / vt_len;
            vec3_add(v_normal, vec3_scale(v_tangent, scale))
        }
    }

    /// Tangential impulse magnitude for a given normal impulse.
    pub fn tangential_impulse(&self, normal_impulse: f64) -> f64 {
        self.mu_kinetic * normal_impulse.abs()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ContinuousCollision
// ─────────────────────────────────────────────────────────────────────────────

/// Linear continuous collision detection (CCD) for a particle vs. a triangle.
///
/// Finds the first time-of-impact (TOI) ∈ \[0, 1\] for a linearly moving
/// point against a linearly moving triangle.
#[derive(Debug, Clone)]
pub struct ContinuousCollision;

impl ContinuousCollision {
    /// Computes the time of impact for a particle moving from `p0` to `p1`
    /// against a **static** triangle `(a, b, c)`.
    ///
    /// Returns `Some(t)` where t ∈ \[0, 1\] is the normalised impact time, or
    /// `None` if no collision occurs in the interval.
    pub fn particle_vs_static_triangle(
        p0: [f64; 3],
        p1: [f64; 3],
        a: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
    ) -> Option<f64> {
        let ab = vec3_sub(b, a);
        let ac = vec3_sub(c, a);
        let n = vec3_normalize(vec3_cross(ab, ac))?;

        let d = vec3_dot(n, a);
        let d0 = vec3_dot(n, p0);
        let d1 = vec3_dot(n, p1);

        // No crossing if both points on the same side.
        let denom = d1 - d0;
        if denom.abs() < 1e-14 {
            return None;
        }

        let t = (d - d0) / denom;
        if !(0.0..=1.0).contains(&t) {
            return None;
        }

        // Impact point.
        let dp = vec3_sub(p1, p0);
        let impact = vec3_add(p0, vec3_scale(dp, t));

        // Check if impact point is inside the triangle using barycentric test.
        let ptc = PointTriangleContact::compute(impact, a, b, c);
        if ptc.dist_sq < 1e-10 { Some(t) } else { None }
    }

    /// Sweeps a sphere of radius `r` from `p0` to `p1` against a static plane.
    ///
    /// Returns `Some(t)` when the sphere first touches the plane `n·x = d`.
    pub fn sphere_vs_plane(
        p0: [f64; 3],
        p1: [f64; 3],
        plane_normal: [f64; 3],
        plane_d: f64,
        radius: f64,
    ) -> Option<f64> {
        let h0 = vec3_dot(plane_normal, p0) - plane_d - radius;
        let h1 = vec3_dot(plane_normal, p1) - plane_d - radius;
        if h0 <= 0.0 {
            return Some(0.0); // already touching
        }
        if h1 >= h0 {
            return None; // moving away
        }
        let t = h0 / (h0 - h1);
        if (0.0..=1.0).contains(&t) {
            Some(t)
        } else {
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CollisionFilter
// ─────────────────────────────────────────────────────────────────────────────

/// Bitmask-based collision group filter.
///
/// An object belongs to `group_mask` and collides with objects whose
/// `group_mask` overlaps with this filter's `collides_with` mask.
#[derive(Debug, Clone, Copy)]
pub struct CollisionFilter {
    /// Bitmask of which groups this object belongs to.
    pub group_mask: u32,
    /// Bitmask of which groups this object should collide with.
    pub collides_with: u32,
}

impl CollisionFilter {
    /// Creates a new collision filter.
    pub fn new(group_mask: u32, collides_with: u32) -> Self {
        Self {
            group_mask,
            collides_with,
        }
    }

    /// Creates a filter that belongs to group `g` and collides with all groups.
    pub fn group(g: u32) -> Self {
        Self {
            group_mask: g,
            collides_with: u32::MAX,
        }
    }

    /// Returns `true` if `self` and `other` should generate a collision.
    ///
    /// Both objects must agree: A collides with B and B collides with A.
    pub fn should_collide(&self, other: &CollisionFilter) -> bool {
        (self.collides_with & other.group_mask) != 0 && (other.collides_with & self.group_mask) != 0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── SoftBodyAabb ──

    #[test]
    fn test_aabb_empty_state() {
        let aabb = SoftBodyAabb::empty();
        assert!(aabb.min[0].is_infinite() && aabb.min[0] > 0.0);
    }

    #[test]
    fn test_aabb_refit_single_point() {
        let mut aabb = SoftBodyAabb::empty();
        aabb.refit(&[[1.0, 2.0, 3.0]]);
        assert!((aabb.min[0] - 1.0).abs() < EPS);
        assert!((aabb.max[2] - 3.0).abs() < EPS);
    }

    #[test]
    fn test_aabb_refit_multiple_points() {
        let mut aabb = SoftBodyAabb::empty();
        aabb.refit(&[[-1.0, 0.0, 0.0], [1.0, 2.0, -3.0]]);
        assert!((aabb.min[0] - (-1.0)).abs() < EPS);
        assert!((aabb.max[1] - 2.0).abs() < EPS);
        assert!((aabb.min[2] - (-3.0)).abs() < EPS);
    }

    #[test]
    fn test_aabb_overlaps_true() {
        let a = SoftBodyAabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = SoftBodyAabb::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_overlaps_false() {
        let a = SoftBodyAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = SoftBodyAabb::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_contains_point() {
        let aabb = SoftBodyAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(aabb.contains_point([0.5, 0.5, 0.5]));
        assert!(!aabb.contains_point([1.5, 0.5, 0.5]));
    }

    #[test]
    fn test_aabb_expand() {
        let aabb = SoftBodyAabb::new([0.0; 3], [1.0; 3]);
        let expanded = aabb.expand(0.5);
        assert!((expanded.min[0] - (-0.5)).abs() < EPS);
        assert!((expanded.max[0] - 1.5).abs() < EPS);
    }

    #[test]
    fn test_aabb_center() {
        let aabb = SoftBodyAabb::new([0.0; 3], [2.0; 3]);
        let c = aabb.center();
        assert!((c[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_aabb_half_extents() {
        let aabb = SoftBodyAabb::new([0.0; 3], [4.0; 3]);
        let he = aabb.half_extents();
        assert!((he[0] - 2.0).abs() < EPS);
    }

    // ── TriangleMeshCollider ──

    #[test]
    fn test_mesh_collider_num_triangles() {
        let tris = vec![
            [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let c = TriangleMeshCollider::new(tris);
        assert_eq!(c.num_triangles(), 2);
    }

    #[test]
    fn test_mesh_collider_broad_query_finds_overlapping() {
        let tris = vec![[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let c = TriangleMeshCollider::new(tris);
        let query = SoftBodyAabb::new([-0.1; 3], [0.5; 3]);
        let hits = c.broad_query(&query);
        assert!(!hits.is_empty());
    }

    #[test]
    fn test_mesh_collider_broad_query_misses_non_overlapping() {
        let tris = vec![[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let c = TriangleMeshCollider::new(tris);
        let query = SoftBodyAabb::new([5.0; 3], [6.0; 3]);
        let hits = c.broad_query(&query);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_mesh_collider_total_aabb() {
        let tris = vec![[[0.0; 3], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]];
        let c = TriangleMeshCollider::new(tris);
        let aabb = c.total_aabb();
        assert!((aabb.max[0] - 3.0).abs() < EPS);
    }

    // ── PointTriangleContact ──

    #[test]
    fn test_point_triangle_above_center() {
        // Point directly above the centroid of a unit triangle.
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let res = PointTriangleContact::compute(p, a, b, c);
        // Closest point should be inside triangle.
        let bary = res.barycentric;
        assert!(bary[0] >= -EPS && bary[1] >= -EPS && bary[2] >= -EPS);
        // Distance should be ~1.0.
        assert!((res.dist_sq.sqrt() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_point_triangle_on_surface() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [0.0, 2.0, 0.0];
        let p = [0.5, 0.5, 0.0]; // lies in triangle
        let res = PointTriangleContact::compute(p, a, b, c);
        assert!(
            res.dist_sq < 1e-10,
            "point on triangle should have dist_sq≈0"
        );
    }

    #[test]
    fn test_point_triangle_closest_is_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [-1.0, -1.0, 0.0]; // outside, closest to vertex a
        let res = PointTriangleContact::compute(p, a, b, c);
        let da = vec3_len(vec3_sub(res.closest_point, a));
        assert!(da < 1e-6, "closest point should be vertex a");
    }

    #[test]
    fn test_point_triangle_normal_nonzero() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let res = PointTriangleContact::compute(p, a, b, c);
        let n_len = vec3_len(res.normal);
        assert!((n_len - 1.0).abs() < 1e-6, "normal should be unit length");
    }

    #[test]
    fn test_point_triangle_signed_depth_positive_above() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let res = PointTriangleContact::compute(p, a, b, c);
        let depth = res.signed_depth(p);
        assert!(depth > 0.0, "point above plane should have positive depth");
    }

    // ── EdgeEdgeContact ──

    #[test]
    fn test_edge_edge_parallel() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let p4 = [1.0, 1.0, 0.0];
        let res = EdgeEdgeContact::compute(p1, p2, p3, p4);
        assert!((res.dist_sq.sqrt() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_edge_edge_perpendicular() {
        // Two skew edges that are closest at (0.5, 0, 0) and (0, 0.5, 0).
        // Actually put them crossing in 3D:
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.5, -1.0, 0.0];
        let p4 = [0.5, 1.0, 0.0]; // crossing edge
        let res = EdgeEdgeContact::compute(p1, p2, p3, p4);
        // Closest distance should be zero (edges cross).
        assert!(res.dist_sq < 1e-10, "crossing edges have dist=0");
    }

    #[test]
    fn test_edge_edge_t_in_range() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [0.0, 1.0, 0.0];
        let p4 = [1.0, 1.0, 0.0];
        let res = EdgeEdgeContact::compute(p1, p2, p3, p4);
        assert!((0.0..=1.0).contains(&res.t1));
        assert!((0.0..=1.0).contains(&res.t2));
    }

    // ── SoftBodyContact ──

    #[test]
    fn test_soft_body_contact_is_penetrating() {
        let c = SoftBodyContact::new([0.0; 3], [0.0, 1.0, 0.0], 0.05, 0, 0);
        assert!(c.is_penetrating(0.01));
        assert!(!c.is_penetrating(0.1));
    }

    // ── CollisionResponse ──

    #[test]
    fn test_penalty_force_direction() {
        let resp = CollisionResponse::penalty(1000.0);
        let contact = SoftBodyContact::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 0);
        let force = resp.penalty_force(&contact);
        // Force should be in +Y direction.
        assert!(force[1] > 0.0);
    }

    #[test]
    fn test_position_correction_moves_particle_out() {
        let resp = CollisionResponse::penalty(1000.0);
        let contact = SoftBodyContact::new([0.0; 3], [0.0, 1.0, 0.0], 0.1, 0, 0);
        let correction = resp.position_correction(&contact);
        assert!((correction[1] - 0.1).abs() < EPS);
    }

    #[test]
    fn test_impulse_response_reflects_velocity() {
        let resp = CollisionResponse::impulse(1.0); // elastic
        let contact = SoftBodyContact::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 0);
        let vel = [0.0, -5.0, 0.0];
        let (_pos, new_vel) = resp.resolve([0.0; 3], vel, &contact);
        // Y-velocity should be positive (reflected) for elastic.
        assert!(
            new_vel[1] > 0.0,
            "velocity should be reflected, got {new_vel:?}"
        );
    }

    #[test]
    fn test_impulse_inelastic_absorbs_velocity() {
        let resp = CollisionResponse::impulse(0.0); // perfectly inelastic
        let contact = SoftBodyContact::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0, 0);
        let vel = [0.0, -5.0, 0.0];
        let (_pos, new_vel) = resp.resolve([0.0; 3], vel, &contact);
        // Y-velocity should be zero after fully inelastic collision.
        assert!(
            new_vel[1].abs() < EPS,
            "inelastic: vy should be 0, got {}",
            new_vel[1]
        );
    }

    // ── FrictionResponse ──

    #[test]
    fn test_friction_zero_tangential_no_change() {
        let fr = FrictionResponse::new(0.5, 0.4);
        // Velocity purely in normal direction.
        let vel = [0.0, 5.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let result = fr.apply(vel, normal);
        // No tangential component → velocity unchanged.
        assert!((result[0]).abs() < EPS && (result[2]).abs() < EPS);
    }

    #[test]
    fn test_friction_reduces_tangential_velocity() {
        let fr = FrictionResponse::new(0.5, 0.4);
        // Velocity: moving tangentially + small normal component.
        let vel = [10.0, 1.0, 0.0];
        let normal = [0.0, 1.0, 0.0];
        let result = fr.apply(vel, normal);
        let vt_before = (vel[0] * vel[0] + vel[2] * vel[2]).sqrt();
        let vt_after = (result[0] * result[0] + result[2] * result[2]).sqrt();
        assert!(
            vt_after <= vt_before + EPS,
            "friction should reduce tangential speed"
        );
    }

    #[test]
    fn test_friction_tangential_impulse_scales_with_normal() {
        let fr = FrictionResponse::new(0.5, 0.3);
        let j1 = fr.tangential_impulse(10.0);
        let j2 = fr.tangential_impulse(20.0);
        assert!((j2 - 2.0 * j1).abs() < EPS);
    }

    // ── ContinuousCollision ──

    #[test]
    fn test_ccd_particle_hits_triangle() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p0 = [0.25, 0.25, 1.0]; // above
        let p1 = [0.25, 0.25, -1.0]; // below → must cross
        let toi = ContinuousCollision::particle_vs_static_triangle(p0, p1, a, b, c);
        assert!(
            toi.is_some(),
            "particle sweeping through triangle should hit"
        );
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t));
    }

    #[test]
    fn test_ccd_particle_misses_triangle() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        // Particle sweeps far outside triangle in XY.
        let p0 = [5.0, 5.0, 1.0];
        let p1 = [5.0, 5.0, -1.0];
        let toi = ContinuousCollision::particle_vs_static_triangle(p0, p1, a, b, c);
        assert!(toi.is_none(), "particle outside triangle should miss");
    }

    #[test]
    fn test_ccd_sphere_hits_plane() {
        let p0 = [0.0, 2.0, 0.0];
        let p1 = [0.0, -2.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let toi = ContinuousCollision::sphere_vs_plane(p0, p1, n, 0.0, 0.1);
        assert!(toi.is_some());
        let t = toi.unwrap();
        // At t, sphere center should be at y = 0.1 (just touching).
        let y = p0[1] + (p1[1] - p0[1]) * t;
        assert!(
            (y - 0.1).abs() < 1e-6,
            "sphere should just touch plane, y={y}"
        );
    }

    #[test]
    fn test_ccd_sphere_misses_plane_moving_away() {
        let p0 = [0.0, 2.0, 0.0];
        let p1 = [0.0, 4.0, 0.0]; // moving away from plane at y=0
        let n = [0.0, 1.0, 0.0];
        let toi = ContinuousCollision::sphere_vs_plane(p0, p1, n, 0.0, 0.1);
        assert!(toi.is_none(), "sphere moving away should not hit");
    }

    // ── CollisionFilter ──

    #[test]
    fn test_collision_filter_same_group_collides() {
        let a = CollisionFilter::group(0b01);
        let b = CollisionFilter::group(0b01);
        assert!(a.should_collide(&b));
    }

    #[test]
    fn test_collision_filter_different_groups_no_collision() {
        let a = CollisionFilter::new(0b01, 0b10); // group 1, collides with group 2
        let b = CollisionFilter::new(0b100, 0b001); // group 4, collides with group 1
        // a wants to hit b? a.collides_with & b.group = 0b10 & 0b100 = 0
        assert!(!a.should_collide(&b));
    }

    #[test]
    fn test_collision_filter_bilateral() {
        let a = CollisionFilter::new(0b01, 0b10);
        let b = CollisionFilter::new(0b10, 0b01);
        assert!(a.should_collide(&b));
        assert!(b.should_collide(&a));
    }

    #[test]
    fn test_collision_filter_self_exclusion() {
        // Object that doesn't collide with its own group.
        let a = CollisionFilter::new(0b01, 0b10); // only collides with group 2
        let a2 = CollisionFilter::new(0b01, 0b10);
        // a & a2 share group 1, but a only collides_with group 2.
        assert!(!a.should_collide(&a2));
    }

    #[test]
    fn test_collision_filter_all_groups() {
        let a = CollisionFilter::new(0b0001, u32::MAX);
        let b = CollisionFilter::new(0b1000, u32::MAX);
        assert!(a.should_collide(&b));
    }

    // ── SelfCollision ──

    /// Positions and triangles of an `n x n` cloth grid in the XZ plane with the
    /// given vertex spacing; vertex `(i, j)` has index `j * n + i`.
    fn grid_cloth(n: usize, spacing: f64) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let mut positions = Vec::with_capacity(n * n);
        for j in 0..n {
            for i in 0..n {
                positions.push([i as f64 * spacing, 0.0, j as f64 * spacing]);
            }
        }
        let mut triangles = Vec::new();
        for j in 0..(n - 1) {
            for i in 0..(n - 1) {
                let a = j * n + i;
                let b = a + 1;
                let c = a + n;
                let d = c + 1;
                triangles.push([a, b, c]);
                triangles.push([b, d, c]);
            }
        }
        (positions, triangles)
    }

    /// The O(P * T) reference that `SelfCollision::detect` replaced.
    fn brute_force_self_contacts(
        positions: &[[f64; 3]],
        triangles: &[[usize; 3]],
        thickness: f64,
    ) -> Vec<SoftBodyContact> {
        let mut contacts = Vec::new();
        for (pi, &pos) in positions.iter().enumerate() {
            for (ti, tri) in triangles.iter().enumerate() {
                if tri.contains(&pi) {
                    continue;
                }
                let result = PointTriangleContact::compute(
                    pos,
                    positions[tri[0]],
                    positions[tri[1]],
                    positions[tri[2]],
                );
                let dist = result.dist_sq.sqrt();
                if dist < thickness {
                    contacts.push(SoftBodyContact::new(
                        result.closest_point,
                        result.normal,
                        thickness - dist,
                        pi,
                        ti,
                    ));
                }
            }
        }
        contacts
    }

    #[test]
    fn test_self_collision_adjacent_skipped() {
        // Two triangles sharing the edge (1, 2).  With a thickness larger than
        // the whole strip every non-adjacent particle/triangle pair collides,
        // so the contact list is exactly the set of non-adjacent pairs.
        let triangles = [[0, 1, 2_usize], [1, 3, 2]];
        let positions = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ];
        let sc = SelfCollision::with_thickness(5.0);
        let contacts = sc.detect(&positions, &triangles);
        let pairs: Vec<(usize, usize)> = contacts
            .iter()
            .map(|c| (c.particle_index, c.triangle_index))
            .collect();
        // Particle 0 is only a corner of triangle 0, particle 3 only of triangle 1;
        // particles 1 and 2 are corners of both triangles.
        assert_eq!(pairs, vec![(0, 1), (3, 0)]);
    }

    #[test]
    fn test_self_collision_detects_particle_above_triangle() {
        let triangles = [[0, 1, 2_usize]];
        let positions = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            // Particle 3 hovers 0.06 above the triangle interior.
            [0.25, 0.25, 0.06],
        ];
        let sc = SelfCollision::with_thickness(0.1);
        let contacts = sc.detect(&positions, &triangles);
        assert_eq!(contacts.len(), 1);
        let contact = &contacts[0];
        assert_eq!(contact.particle_index, 3);
        assert_eq!(contact.triangle_index, 0);
        assert!((contact.penetration_depth - 0.04).abs() < 1e-12);
        assert!((contact.contact_point[0] - 0.25).abs() < 1e-12);
        assert!((contact.contact_point[1] - 0.25).abs() < 1e-12);
        assert!(contact.contact_point[2].abs() < 1e-12);

        // Outside the guard band there is no contact.
        let far = SelfCollision::with_thickness(0.05);
        assert!(far.detect(&positions, &triangles).is_empty());
    }

    #[test]
    fn test_self_collision_matches_brute_force_on_folded_cloth() {
        // A 20 x 20 cloth folded in half: the second half is mirrored onto the
        // first at a height of 0.04, closer than the 0.08 thickness.
        let n = 20;
        let spacing = 0.1;
        let (mut positions, triangles) = grid_cloth(n, spacing);
        let fold_x = (n - 1) as f64 * spacing;
        for j in 0..n {
            for i in (n / 2)..n {
                let p = &mut positions[j * n + i];
                p[0] = fold_x - p[0];
                p[1] = 0.04;
            }
        }
        let thickness = 0.08;
        let expected = brute_force_self_contacts(&positions, &triangles, thickness);
        assert!(
            expected.len() > 200,
            "the folded cloth must produce many contacts, got {}",
            expected.len()
        );
        let found = SelfCollision::with_thickness(thickness).detect(&positions, &triangles);
        assert_eq!(found.len(), expected.len());
        for (got, want) in found.iter().zip(expected.iter()) {
            assert_eq!(got.particle_index, want.particle_index);
            assert_eq!(got.triangle_index, want.triangle_index);
            assert_eq!(got.contact_point, want.contact_point);
            assert_eq!(got.normal, want.normal);
            assert_eq!(got.penetration_depth, want.penetration_depth);
        }
    }

    #[test]
    fn test_self_collision_matches_brute_force_on_crumpled_cloth() {
        // A deterministic jitter in all three axes crumples a 16 x 16 sheet, so
        // the closest feature of a particle is a triangle's face, edge or corner
        // in turn, and the particles are scattered over many BVH nodes.
        let n = 16;
        let (mut positions, triangles) = grid_cloth(n, 0.1);
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next_offset = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        };
        for position in positions.iter_mut() {
            for coordinate in position.iter_mut() {
                *coordinate += 0.12 * next_offset();
            }
        }
        let thickness = 0.06;
        let expected = brute_force_self_contacts(&positions, &triangles, thickness);
        assert!(
            expected.len() > 100,
            "the crumpled cloth must produce many contacts, got {}",
            expected.len()
        );
        let found = SelfCollision::with_thickness(thickness).detect(&positions, &triangles);
        assert_eq!(found.len(), expected.len());
        for (got, want) in found.iter().zip(expected.iter()) {
            assert_eq!(got.particle_index, want.particle_index);
            assert_eq!(got.triangle_index, want.triangle_index);
            assert_eq!(got.contact_point, want.contact_point);
            assert_eq!(got.normal, want.normal);
            assert_eq!(got.penetration_depth, want.penetration_depth);
        }
    }

    #[test]
    fn test_self_collision_flat_cloth_has_no_contacts() {
        let (positions, triangles) = grid_cloth(12, 0.1);
        let sc = SelfCollision::with_thickness(0.03);
        assert!(sc.detect(&positions, &triangles).is_empty());
        assert!(brute_force_self_contacts(&positions, &triangles, 0.03).is_empty());
    }

    #[test]
    fn test_self_collision_skips_invalid_triangles_and_thickness() {
        let positions = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            // Particle 3 hovers over triangle 0.
            [0.25, 0.25, 0.05],
            // Particle 4 is not finite.
            [f64::NAN, 0.0, 0.0],
        ];
        // Triangle 1 references a particle that does not exist and triangle 2 has
        // a non-finite corner: both are ignored without panicking.
        let triangles = [[0, 1, 2_usize], [0, 1, 99], [0, 1, 4]];
        let sc = SelfCollision::with_thickness(0.1);
        let pairs: Vec<(usize, usize)> = sc
            .detect(&positions, &triangles)
            .iter()
            .map(|c| (c.particle_index, c.triangle_index))
            .collect();
        assert_eq!(pairs, vec![(3, 0)]);
        for bad in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(
                SelfCollision::with_thickness(bad)
                    .detect(&positions, &triangles)
                    .is_empty(),
                "thickness {bad} must yield no contacts"
            );
        }
        assert!(sc.detect(&[], &[]).is_empty());
    }
}
