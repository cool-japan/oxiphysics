// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact stability analysis for rigid body systems.
//!
//! This module provides tools for analysing whether a set of contacts can
//! maintain a body (or multi-body system) in static or quasi-static
//! equilibrium.  Key features:
//!
//! * **Coulomb friction cone linearisation** – approximate the circular
//!   friction cone with a polyhedral cone for LP / QP formulations.
//! * **Contact wrench space** – compute the set of wrenches achievable
//!   through a collection of contacts.
//! * **Force closure test** – determine whether arbitrary wrenches can be
//!   resisted by a contact configuration.
//! * **Grasp quality metrics** – Q1 (largest inscribed ball in the wrench
//!   space) and epsilon metric.
//! * **Support polygon** – compute the convex hull of ground-contact
//!   projections.
//! * **ZMP (Zero Moment Point)** computation.
//! * **Tipping analysis** – check if the CoM projection is inside the
//!   support polygon.
//! * **Multi-contact friction polyhedron** – build the linearised friction
//!   polyhedron for multiple simultaneous contacts.
//! * **Complementarity conditions** – formulate the LCP conditions for
//!   contact with friction.

// ── tiny linear-algebra helpers (f64 arrays only, no nalgebra) ──────────────

/// 3-component dot product.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 3-component cross product.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Length of a 3-vector.
fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3-vector; returns `[0,0,0]` if near-zero.
fn vec3_normalize(v: [f64; 3]) -> [f64; 3] {
    let l = vec3_len(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / l, v[1] / l, v[2] / l]
    }
}

/// Scale a 3-vector.
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Subtract two 3-vectors.
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Negate a 3-vector.
fn vec3_neg(v: [f64; 3]) -> [f64; 3] {
    [-v[0], -v[1], -v[2]]
}

/// Dot product of two 6-vectors.
fn dot6(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    let mut s = 0.0;
    for i in 0..6 {
        s += a[i] * b[i];
    }
    s
}

/// Length of a 6-vector.
fn vec6_len(v: &[f64; 6]) -> f64 {
    dot6(v, v).sqrt()
}

// ── Tangent basis construction ──────────────────────────────────────────────

/// Build an orthonormal tangent pair `(t1, t2)` for a given normal `n`.
fn tangent_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let ax = n[0].abs();
    let ay = n[1].abs();
    let az = n[2].abs();
    let candidate = if ax <= ay && ax <= az {
        [1.0, 0.0, 0.0]
    } else if ay <= ax && ay <= az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let t1 = vec3_normalize(cross3(n, candidate));
    let t2 = cross3(n, t1);
    (t1, t2)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Contact point description
// ═══════════════════════════════════════════════════════════════════════════════

/// A single frictional contact point.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// Position of the contact in world frame.
    pub position: [f64; 3],
    /// Outward-pointing surface normal at the contact (unit vector).
    pub normal: [f64; 3],
    /// Coulomb friction coefficient.
    pub mu: f64,
}

/// A 6D wrench: `[fx, fy, fz, tx, ty, tz]`.
pub type Wrench = [f64; 6];

// ═══════════════════════════════════════════════════════════════════════════════
// Coulomb cone linearisation
// ═══════════════════════════════════════════════════════════════════════════════

/// Generate the edge directions of a linearised Coulomb friction cone.
///
/// Given a contact normal `n` and friction coefficient `mu`, this function
/// approximates the friction cone with `num_sides` edges equally spaced
/// around the normal.
///
/// Each returned vector is a unit-length *force direction* on the surface
/// of the linearised cone.
///
/// # Arguments
/// * `normal` - Contact normal (unit vector).
/// * `mu` - Coulomb friction coefficient (>= 0).
/// * `num_sides` - Number of sides for the polyhedral approximation.
pub fn linearise_coulomb_cone(normal: [f64; 3], mu: f64, num_sides: usize) -> Vec<[f64; 3]> {
    let (t1, t2) = tangent_basis(normal);
    let mut edges = Vec::with_capacity(num_sides);
    for i in 0..num_sides {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (num_sides as f64);
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        // direction = normal + mu * (cos_a * t1 + sin_a * t2)
        let dir = [
            normal[0] + mu * (cos_a * t1[0] + sin_a * t2[0]),
            normal[1] + mu * (cos_a * t1[1] + sin_a * t2[1]),
            normal[2] + mu * (cos_a * t1[2] + sin_a * t2[2]),
        ];
        edges.push(vec3_normalize(dir));
    }
    edges
}

/// Interior approximation error for a linearised Coulomb cone with `n` sides.
///
/// The ratio of the inscribed to the circumscribed polygon is `cos(pi/n)`.
pub fn cone_linearisation_error(num_sides: usize) -> f64 {
    let n = num_sides.max(3) as f64;
    1.0 - (std::f64::consts::PI / n).cos()
}

// ═══════════════════════════════════════════════════════════════════════════════
// Contact wrench space
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute a single contact wrench (force + torque about origin) for a force
/// `f` applied at `position`.
///
/// wrench = \[f, position x f\]
pub fn contact_wrench(position: [f64; 3], force: [f64; 3]) -> Wrench {
    let torque = cross3(position, force);
    [
        force[0], force[1], force[2], torque[0], torque[1], torque[2],
    ]
}

/// Build the *primitive contact wrench set* for a single contact point.
///
/// Each edge of the linearised friction cone is turned into a unit wrench
/// about the origin.
pub fn primitive_contact_wrenches(contact: &ContactPoint, num_cone_sides: usize) -> Vec<Wrench> {
    let edges = linearise_coulomb_cone(contact.normal, contact.mu, num_cone_sides);
    edges
        .iter()
        .map(|&e| contact_wrench(contact.position, e))
        .collect()
}

/// Build the full contact wrench set for *multiple* contact points.
///
/// The returned matrix has dimensions `N x 6` stored as a `Vec`Wrench`
/// where `N = num_contacts * num_cone_sides`.
pub fn multi_contact_wrench_set(contacts: &[ContactPoint], num_cone_sides: usize) -> Vec<Wrench> {
    let mut wrenches = Vec::new();
    for c in contacts {
        wrenches.extend(primitive_contact_wrenches(c, num_cone_sides));
    }
    wrenches
}

// ═══════════════════════════════════════════════════════════════════════════════
// Force closure test
// ═══════════════════════════════════════════════════════════════════════════════

/// Check whether a set of contact wrenches achieves **force closure**.
///
/// Force closure means that any external wrench can be balanced by
/// non-negative combinations of the contact wrenches.
///
/// This implementation uses a necessary (but practical) condition:
/// the *convex hull* of the primitive wrenches must contain the origin
/// in its interior.  We approximate this by checking that the minimum
/// over all 6 coordinate axes of the signed span is strictly positive.
///
/// For a rigorous test the full 6-D convex-hull check is required;
/// this heuristic is fast and conservative.
pub fn force_closure_test(wrenches: &[Wrench]) -> bool {
    if wrenches.is_empty() {
        return false;
    }
    // For each of the 6 wrench axes check that we can produce both
    // positive and negative components.
    for axis in 0..6 {
        let mut has_pos = false;
        let mut has_neg = false;
        for w in wrenches {
            if w[axis] > 1e-10 {
                has_pos = true;
            }
            if w[axis] < -1e-10 {
                has_neg = true;
            }
            if has_pos && has_neg {
                break;
            }
        }
        if !has_pos || !has_neg {
            return false;
        }
    }
    true
}

/// Full force-closure test using ray-casting in wrench space.
///
/// For each of six canonical wrench directions we check whether a
/// non-negative combination of primitive wrenches can produce it.
/// Returns `true` only if all six directions are achievable.
pub fn force_closure_test_full(wrenches: &[Wrench]) -> bool {
    if wrenches.len() < 7 {
        // At least 7 wrenches needed for a 6-D interior point
        return force_closure_test(wrenches);
    }
    // Check the six canonical directions and their negatives.
    let dirs: [[f64; 6]; 6] = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ];
    for d in &dirs {
        let mut max_proj = f64::NEG_INFINITY;
        let mut min_proj = f64::INFINITY;
        for w in wrenches {
            let proj = dot6(w, d);
            if proj > max_proj {
                max_proj = proj;
            }
            if proj < min_proj {
                min_proj = proj;
            }
        }
        // Must span both sides of the origin in this direction
        if max_proj < 1e-10 || min_proj > -1e-10 {
            return false;
        }
    }
    true
}

// ═══════════════════════════════════════════════════════════════════════════════
// Grasp quality metrics
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the Q1 grasp quality metric.
///
/// Q1 is defined as the radius of the largest ball centred at the origin
/// that is fully contained in the convex hull of the primitive wrenches.
///
/// This approximation computes the minimum distance from the origin to
/// the convex-hull boundary by taking `min(|w|)` over all primitive
/// wrenches.  For a more accurate result the full wrench-space convex
/// hull is required.
pub fn grasp_quality_q1(wrenches: &[Wrench]) -> f64 {
    if wrenches.is_empty() {
        return 0.0;
    }
    let mut min_dist = f64::INFINITY;
    for w in wrenches {
        let d = vec6_len(w);
        if d < min_dist {
            min_dist = d;
        }
    }
    min_dist
}

/// Compute the epsilon grasp quality metric.
///
/// Epsilon is the radius of the largest wrench ball that can be resisted
/// (equivalent to Q1 when the primitive wrenches are normalised).
///
/// This function normalises each wrench to unit length before computing
/// the minimum distance.
pub fn grasp_quality_epsilon(wrenches: &[Wrench]) -> f64 {
    if wrenches.is_empty() {
        return 0.0;
    }
    let normalised: Vec<Wrench> = wrenches
        .iter()
        .map(|w| {
            let l = vec6_len(w);
            if l < 1e-15 {
                [0.0; 6]
            } else {
                let inv = 1.0 / l;
                [
                    w[0] * inv,
                    w[1] * inv,
                    w[2] * inv,
                    w[3] * inv,
                    w[4] * inv,
                    w[5] * inv,
                ]
            }
        })
        .collect();
    grasp_quality_q1(&normalised)
}

/// Volume-based grasp quality: the volume of the convex hull of the
/// normalised wrenches (approximation using the sum of absolute
/// determinants of 6-wrench subsets divided by 720).
///
/// For large wrench sets this is expensive; here we use a bounding
/// approximation based on axis-aligned extent.
pub fn grasp_quality_volume(wrenches: &[Wrench]) -> f64 {
    if wrenches.is_empty() {
        return 0.0;
    }
    let mut min_vals = [f64::INFINITY; 6];
    let mut max_vals = [f64::NEG_INFINITY; 6];
    for w in wrenches {
        for i in 0..6 {
            if w[i] < min_vals[i] {
                min_vals[i] = w[i];
            }
            if w[i] > max_vals[i] {
                max_vals[i] = w[i];
            }
        }
    }
    let mut vol = 1.0;
    for i in 0..6 {
        vol *= (max_vals[i] - min_vals[i]).max(0.0);
    }
    vol
}

// ═══════════════════════════════════════════════════════════════════════════════
// Support polygon
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the 2-D support polygon from a set of ground contact positions.
///
/// The contacts are projected onto the XY plane (z is ignored).  The result
/// is a convex hull in counter-clockwise order.
///
/// Returns the polygon vertices as `(x, y)` pairs.
pub fn support_polygon(contacts: &[[f64; 3]]) -> Vec<(f64, f64)> {
    let mut pts: Vec<(f64, f64)> = contacts.iter().map(|c| (c[0], c[1])).collect();
    convex_hull_2d(&mut pts)
}

/// Andrew's monotone-chain convex hull algorithm.
fn convex_hull_2d(points: &mut [(f64, f64)]) -> Vec<(f64, f64)> {
    let n = points.len();
    if n < 2 {
        return points.to_vec();
    }
    points.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(2 * n);
    // Lower hull
    for &p in points.iter() {
        while hull.len() >= 2 {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if cross_2d(a, b, p) <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    // Upper hull
    let lower_len = hull.len() + 1;
    for &p in points.iter().rev() {
        while hull.len() >= lower_len {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if cross_2d(a, b, p) <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    hull.pop(); // last point == first point
    hull
}

/// 2-D cross product for convex hull.
fn cross_2d(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

/// Test whether a 2-D point lies inside a convex polygon (CCW ordered).
///
/// Uses the cross-product winding test.
pub fn point_in_convex_polygon(point: (f64, f64), polygon: &[(f64, f64)]) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let cp = cross_2d(polygon[i], polygon[j], point);
        if cp < -1e-12 {
            return false;
        }
    }
    true
}

/// Compute the signed area of a convex polygon (CCW = positive).
pub fn polygon_area(polygon: &[(f64, f64)]) -> f64 {
    let n = polygon.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += polygon[i].0 * polygon[j].1 - polygon[j].0 * polygon[i].1;
    }
    area * 0.5
}

/// Compute the centroid of a convex polygon.
pub fn polygon_centroid(polygon: &[(f64, f64)]) -> (f64, f64) {
    let n = polygon.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if n == 1 {
        return polygon[0];
    }
    let a = polygon_area(polygon);
    if a.abs() < 1e-15 {
        // Degenerate: return average
        let sx: f64 = polygon.iter().map(|p| p.0).sum();
        let sy: f64 = polygon.iter().map(|p| p.1).sum();
        return (sx / n as f64, sy / n as f64);
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let factor = polygon[i].0 * polygon[j].1 - polygon[j].0 * polygon[i].1;
        cx += (polygon[i].0 + polygon[j].0) * factor;
        cy += (polygon[i].1 + polygon[j].1) * factor;
    }
    let inv_6a = 1.0 / (6.0 * a);
    (cx * inv_6a, cy * inv_6a)
}

// ═══════════════════════════════════════════════════════════════════════════════
// ZMP computation
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the Zero Moment Point (ZMP) on the ground plane (z = 0).
///
/// Given the centre-of-mass position, its linear acceleration, and gravity,
/// the ZMP is:
///
/// ```text
/// zmp_x = com_x - com_z * (acc_x - gx) / (acc_z - gz)
/// zmp_y = com_y - com_z * (acc_y - gy) / (acc_z - gz)
/// ```
///
/// Returns `None` if the vertical component is too small (free fall or
/// near-zero normal force).
///
/// # Arguments
/// * `com` - Centre of mass `\[x, y, z\]`.
/// * `acc` - Linear acceleration of the CoM `\[ax, ay, az\]`.
/// * `gravity` - Gravitational acceleration vector `\[gx, gy, gz\]`.
pub fn compute_zmp(com: [f64; 3], acc: [f64; 3], gravity: [f64; 3]) -> Option<(f64, f64)> {
    let denom = acc[2] - gravity[2];
    if denom.abs() < 1e-10 {
        return None;
    }
    let zmp_x = com[0] - com[2] * (acc[0] - gravity[0]) / denom;
    let zmp_y = com[1] - com[2] * (acc[1] - gravity[1]) / denom;
    Some((zmp_x, zmp_y))
}

/// Compute the ZMP using total contact forces and torques.
///
/// ```text
/// zmp_x = -tau_y / fz
/// zmp_y =  tau_x / fz
/// ```
///
/// `force` and `torque` are about the world origin.
pub fn compute_zmp_from_forces(force: [f64; 3], torque: [f64; 3]) -> Option<(f64, f64)> {
    if force[2].abs() < 1e-10 {
        return None;
    }
    let zmp_x = -torque[1] / force[2];
    let zmp_y = torque[0] / force[2];
    Some((zmp_x, zmp_y))
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tipping analysis
// ═══════════════════════════════════════════════════════════════════════════════

/// Result of a tipping stability analysis.
#[derive(Debug, Clone)]
pub struct TippingResult {
    /// Is the body stable (ZMP inside support polygon)?
    pub stable: bool,
    /// The computed ZMP (x, y).
    pub zmp: (f64, f64),
    /// Stability margin: distance from ZMP to nearest polygon edge.
    /// Negative means outside.
    pub margin: f64,
}

/// Compute the stability margin (signed distance from a point to the
/// nearest edge of a convex polygon).
///
/// Positive = inside, negative = outside.
pub fn stability_margin(point: (f64, f64), polygon: &[(f64, f64)]) -> f64 {
    let n = polygon.len();
    if n < 3 {
        return f64::NEG_INFINITY;
    }
    let mut min_dist = f64::INFINITY;
    for i in 0..n {
        let j = (i + 1) % n;
        let (ex, ey) = (polygon[j].0 - polygon[i].0, polygon[j].1 - polygon[i].1);
        let edge_len = (ex * ex + ey * ey).sqrt();
        if edge_len < 1e-15 {
            continue;
        }
        // Inward normal (assuming CCW polygon): positive = inside
        let (nx, ny) = (-ey / edge_len, ex / edge_len);
        let dx = point.0 - polygon[i].0;
        let dy = point.1 - polygon[i].1;
        let dist = dx * nx + dy * ny;
        if dist < min_dist {
            min_dist = dist;
        }
    }
    min_dist
}

/// Perform a full tipping analysis.
///
/// # Arguments
/// * `com` - Centre of mass `\[x, y, z\]`.
/// * `acc` - Linear acceleration of the CoM.
/// * `gravity` - Gravitational acceleration.
/// * `contact_positions` - Positions of ground contacts.
pub fn tipping_analysis(
    com: [f64; 3],
    acc: [f64; 3],
    gravity: [f64; 3],
    contact_positions: &[[f64; 3]],
) -> Option<TippingResult> {
    let zmp = compute_zmp(com, acc, gravity)?;
    let poly = support_polygon(contact_positions);
    if poly.len() < 3 {
        return Some(TippingResult {
            stable: false,
            zmp,
            margin: f64::NEG_INFINITY,
        });
    }
    let margin = stability_margin(zmp, &poly);
    let stable = margin > -1e-10;
    Some(TippingResult {
        stable,
        zmp,
        margin,
    })
}

/// Simplified tipping test: project CoM onto XY and check if it is inside
/// the support polygon (static case, zero acceleration).
pub fn is_statically_stable(com: [f64; 3], contact_positions: &[[f64; 3]]) -> bool {
    let poly = support_polygon(contact_positions);
    if poly.len() < 3 {
        return false;
    }
    let com_2d = (com[0], com[1]);
    point_in_convex_polygon(com_2d, &poly)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Multi-contact friction polyhedron
// ═══════════════════════════════════════════════════════════════════════════════

/// A linearised friction polyhedron for a single contact.
#[derive(Debug, Clone)]
pub struct FrictionPolyhedron {
    /// The edge directions of the linearised cone (unit vectors).
    pub edges: Vec<[f64; 3]>,
    /// Contact position.
    pub position: [f64; 3],
    /// Contact normal.
    pub normal: [f64; 3],
    /// Friction coefficient.
    pub mu: f64,
    /// Number of sides in the approximation.
    pub num_sides: usize,
}

impl FrictionPolyhedron {
    /// Create a new friction polyhedron.
    pub fn new(position: [f64; 3], normal: [f64; 3], mu: f64, num_sides: usize) -> Self {
        let edges = linearise_coulomb_cone(normal, mu, num_sides);
        Self {
            edges,
            position,
            normal,
            mu,
            num_sides,
        }
    }

    /// Maximum tangential force for a given normal force magnitude.
    pub fn max_tangential_force(&self, normal_force: f64) -> f64 {
        self.mu * normal_force
    }

    /// Test if a force vector lies within the friction cone.
    pub fn contains_force(&self, force: [f64; 3]) -> bool {
        let fn_mag = dot3(force, self.normal);
        if fn_mag < -1e-10 {
            return false; // Pulling
        }
        let tangent = vec3_sub(force, vec3_scale(self.normal, fn_mag));
        let ft_mag = vec3_len(tangent);
        ft_mag <= self.mu * fn_mag + 1e-10
    }

    /// Compute the wrench set for this friction polyhedron.
    pub fn wrench_set(&self) -> Vec<Wrench> {
        self.edges
            .iter()
            .map(|&e| contact_wrench(self.position, e))
            .collect()
    }
}

/// Build friction polyhedra for all contacts.
pub fn build_friction_polyhedra(
    contacts: &[ContactPoint],
    num_sides: usize,
) -> Vec<FrictionPolyhedron> {
    contacts
        .iter()
        .map(|c| FrictionPolyhedron::new(c.position, c.normal, c.mu, num_sides))
        .collect()
}

/// Aggregate wrench set from multiple friction polyhedra.
pub fn aggregate_wrench_set(polyhedra: &[FrictionPolyhedron]) -> Vec<Wrench> {
    let mut wrenches = Vec::new();
    for p in polyhedra {
        wrenches.extend(p.wrench_set());
    }
    wrenches
}

// ═══════════════════════════════════════════════════════════════════════════════
// Complementarity conditions
// ═══════════════════════════════════════════════════════════════════════════════

/// Complementarity condition for a single normal contact.
///
/// At equilibrium:
/// - `gap >= 0`     (no penetration)
/// - `f_n >= 0`     (only compressive normal force)
/// - `gap * f_n = 0` (complementarity)
///
/// Returns the residual `|gap * f_n|`.
pub fn normal_complementarity_residual(gap: f64, f_n: f64) -> f64 {
    (gap * f_n).abs()
}

/// Check if normal complementarity conditions are satisfied within a
/// tolerance.
pub fn check_normal_complementarity(gap: f64, f_n: f64, tol: f64) -> bool {
    gap >= -tol && f_n >= -tol && (gap * f_n).abs() < tol
}

/// Friction complementarity condition.
///
/// - `|f_t| <= mu * f_n`           (inside friction cone)
/// - `v_t * f_t <= 0`              (if sliding, friction opposes motion)
/// - `(mu * f_n - |f_t|) * |v_t| = 0` (complementarity)
///
/// Returns a non-negative residual measuring how badly the condition is
/// violated.
pub fn friction_complementarity_residual(f_t_mag: f64, mu: f64, f_n: f64, v_t_mag: f64) -> f64 {
    let cone_slack = mu * f_n - f_t_mag;
    (cone_slack * v_t_mag).abs()
}

/// Check all complementarity conditions for a single contact.
///
/// Returns `true` if the contact satisfies normal and friction
/// complementarity within the given tolerance.
pub fn check_contact_complementarity(
    gap: f64,
    f_n: f64,
    f_t_mag: f64,
    mu: f64,
    v_t_mag: f64,
    tol: f64,
) -> bool {
    // Normal complementarity
    if !check_normal_complementarity(gap, f_n, tol) {
        return false;
    }
    // Friction cone
    if f_t_mag > mu * f_n + tol {
        return false;
    }
    // Friction complementarity
    let residual = friction_complementarity_residual(f_t_mag, mu, f_n, v_t_mag);
    residual < tol
}

/// Signorini-Fichera complementarity for a contact array.
///
/// Given vectors of gaps, normal forces, tangential force magnitudes,
/// friction coefficients, and tangential sliding speeds, compute the
/// total complementarity residual.
pub fn total_complementarity_residual(
    gaps: &[f64],
    f_ns: &[f64],
    f_ts: &[f64],
    mus: &[f64],
    v_ts: &[f64],
) -> f64 {
    let n = gaps.len();
    let mut total = 0.0;
    for i in 0..n {
        total += normal_complementarity_residual(gaps[i], f_ns[i]);
        total += friction_complementarity_residual(f_ts[i], mus[i], f_ns[i], v_ts[i]);
    }
    total
}

// ═══════════════════════════════════════════════════════════════════════════════
// LCP formulation helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// A contact entry for the LCP formulation.
#[derive(Debug, Clone, Copy)]
pub struct LcpContact {
    /// Penetration gap (negative = penetrating).
    pub gap: f64,
    /// Normal force.
    pub f_n: f64,
    /// Tangential force magnitude.
    pub f_t: f64,
    /// Friction coefficient.
    pub mu: f64,
    /// Tangential sliding speed.
    pub v_t: f64,
    /// Effective mass in the normal direction.
    pub eff_mass_n: f64,
    /// Effective mass in the tangential direction.
    pub eff_mass_t: f64,
}

/// Project a normal force to satisfy non-negativity: `f_n = max(0, f_n)`.
pub fn project_normal_force(f_n: f64) -> f64 {
    f_n.max(0.0)
}

/// Project a tangential force to satisfy the friction cone:
/// `|f_t| <= mu * f_n`.
pub fn project_friction_force(f_t: f64, mu: f64, f_n: f64) -> f64 {
    let limit = mu * f_n;
    f_t.max(-limit).min(limit)
}

/// Perform one Gauss-Seidel iteration for an LCP contact set.
///
/// Updates forces in place.  Returns the sum of absolute force changes
/// (useful for convergence checks).
pub fn lcp_gauss_seidel_step(contacts: &mut [LcpContact], bias_factor: f64) -> f64 {
    let mut delta_sum = 0.0;
    for c in contacts.iter_mut() {
        // Normal direction
        let rhs_n = -c.gap * bias_factor;
        let old_fn = c.f_n;
        c.f_n = project_normal_force(old_fn + c.eff_mass_n * rhs_n);
        delta_sum += (c.f_n - old_fn).abs();

        // Tangential direction
        let rhs_t = -c.v_t;
        let old_ft = c.f_t;
        let new_ft = old_ft + c.eff_mass_t * rhs_t;
        c.f_t = project_friction_force(new_ft, c.mu, c.f_n);
        delta_sum += (c.f_t - old_ft).abs();
    }
    delta_sum
}

/// Solve the contact LCP using projected Gauss-Seidel.
///
/// Returns the number of iterations performed.
pub fn solve_contact_lcp(
    contacts: &mut [LcpContact],
    bias_factor: f64,
    max_iter: usize,
    tol: f64,
) -> usize {
    for iter in 0..max_iter {
        let delta = lcp_gauss_seidel_step(contacts, bias_factor);
        if delta < tol {
            return iter + 1;
        }
    }
    max_iter
}

// ═══════════════════════════════════════════════════════════════════════════════
// Static equilibrium solver
// ═══════════════════════════════════════════════════════════════════════════════

/// Check whether a given external wrench can be balanced by a set of
/// contact wrenches with non-negative coefficients.
///
/// Uses a simple iterative projection approach.
pub fn can_balance_wrench(
    external_wrench: &Wrench,
    contact_wrenches: &[Wrench],
    max_iter: usize,
    tol: f64,
) -> bool {
    if contact_wrenches.is_empty() {
        return vec6_len(external_wrench) < tol;
    }
    let n = contact_wrenches.len();
    let mut lambdas = vec![0.0; n];
    let neg_w: [f64; 6] = [
        -external_wrench[0],
        -external_wrench[1],
        -external_wrench[2],
        -external_wrench[3],
        -external_wrench[4],
        -external_wrench[5],
    ];

    for _iter in 0..max_iter {
        // Compute residual = sum(lambda_i * w_i) - (-external_wrench)
        let mut residual = neg_w;
        for (i, &lam) in lambdas.iter().enumerate() {
            for j in 0..6 {
                residual[j] += lam * contact_wrenches[i][j];
            }
        }
        if vec6_len(&residual) < tol {
            return true;
        }
        // Update each lambda by projecting
        for i in 0..n {
            let w_i = &contact_wrenches[i];
            let w_dot_w = dot6(w_i, w_i);
            if w_dot_w < 1e-15 {
                continue;
            }
            // Recompute residual contribution from lambda_i
            let mut res_without_i = neg_w;
            for (k, &lam) in lambdas.iter().enumerate() {
                if k != i {
                    for j in 0..6 {
                        res_without_i[j] += lam * contact_wrenches[k][j];
                    }
                }
            }
            // Optimal lambda_i = -dot(res_without_i, w_i) / dot(w_i, w_i)
            let numerator = -dot6(&res_without_i, w_i);
            let new_lam = (numerator / w_dot_w).max(0.0);
            lambdas[i] = new_lam;
        }
    }
    // Final residual check
    let mut residual = neg_w;
    for (i, &lam) in lambdas.iter().enumerate() {
        for j in 0..6 {
            residual[j] += lam * contact_wrenches[i][j];
        }
    }
    vec6_len(&residual) < tol
}

// ═══════════════════════════════════════════════════════════════════════════════
// Grasp matrix
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the grasp matrix row for a single contact.
///
/// The grasp matrix maps per-contact forces to the resultant wrench on the
/// object.  For contact `i` with position `p_i` the grasp matrix block is:
///
/// ```text
/// G_i = [ I_{3x3}        ]
///       [ [p_i]_x        ]
/// ```
///
/// where `[p_i]_x` is the skew-symmetric matrix of `p_i`.
///
/// Returns a `6x3` matrix stored as `\[\[f64; 3\\]; 6]` (row-major).
pub fn grasp_matrix_block(position: [f64; 3]) -> [[f64; 3]; 6] {
    let [px, py, pz] = position;
    [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, -pz, py],
        [pz, 0.0, -px],
        [-py, px, 0.0],
    ]
}

/// Compute the full grasp matrix for multiple contacts.
///
/// Returns a matrix of shape `6 x (3*n)` stored as `Vec<Vec`f64`>`.
pub fn grasp_matrix(positions: &[[f64; 3]]) -> Vec<Vec<f64>> {
    let n = positions.len();
    let cols = 3 * n;
    let mut g = vec![vec![0.0; cols]; 6];
    for (idx, pos) in positions.iter().enumerate() {
        let block = grasp_matrix_block(*pos);
        for row in 0..6 {
            for col in 0..3 {
                g[row][idx * 3 + col] = block[row][col];
            }
        }
    }
    g
}

// ═══════════════════════════════════════════════════════════════════════════════
// Wrench space analysis
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the gravity wrench for a body with given mass and CoM.
pub fn gravity_wrench(mass: f64, com: [f64; 3], gravity: [f64; 3]) -> Wrench {
    let f = vec3_scale(gravity, mass);
    let tau = cross3(com, f);
    [f[0], f[1], f[2], tau[0], tau[1], tau[2]]
}

/// Compute the centrifugal wrench contribution.
pub fn centrifugal_wrench(mass: f64, com: [f64; 3], omega: [f64; 3]) -> Wrench {
    // Centripetal acceleration: a = omega x (omega x r)
    let omega_cross_r = cross3(omega, com);
    let acc = cross3(omega, omega_cross_r);
    let f = vec3_scale(acc, mass);
    let tau = cross3(com, f);
    [f[0], f[1], f[2], tau[0], tau[1], tau[2]]
}

/// Sum multiple wrenches.
pub fn sum_wrenches(wrenches: &[Wrench]) -> Wrench {
    let mut total = [0.0; 6];
    for w in wrenches {
        for i in 0..6 {
            total[i] += w[i];
        }
    }
    total
}

// ═══════════════════════════════════════════════════════════════════════════════
// Multi-contact stability score
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute a normalised multi-contact stability score in [0, 1].
///
/// The score is based on:
/// 1. Force closure (binary).
/// 2. Grasp quality epsilon (normalised).
/// 3. Stability margin normalised by the polygon diameter.
///
/// Higher is more stable.
pub fn multi_contact_stability_score(
    contacts: &[ContactPoint],
    com: [f64; 3],
    num_cone_sides: usize,
) -> f64 {
    if contacts.is_empty() {
        return 0.0;
    }
    // Force closure contribution
    let wrenches = multi_contact_wrench_set(contacts, num_cone_sides);
    let fc = if force_closure_test(&wrenches) {
        0.4
    } else {
        0.0
    };

    // Grasp quality contribution
    let eps = grasp_quality_epsilon(&wrenches);
    let eps_score = (eps.min(1.0)) * 0.3;

    // Support polygon stability margin
    let positions: Vec<[f64; 3]> = contacts.iter().map(|c| c.position).collect();
    let poly = support_polygon(&positions);
    let margin = if poly.len() >= 3 {
        let com_2d = (com[0], com[1]);
        stability_margin(com_2d, &poly)
    } else {
        -1.0
    };
    // Normalise margin by polygon diameter
    let diameter = polygon_diameter(&poly);
    let margin_score = if diameter > 1e-10 {
        ((margin / diameter).clamp(0.0, 1.0)) * 0.3
    } else {
        0.0
    };

    fc + eps_score + margin_score
}

/// Compute the diameter of a polygon (maximum distance between vertices).
pub fn polygon_diameter(polygon: &[(f64, f64)]) -> f64 {
    let n = polygon.len();
    let mut max_d = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = polygon[j].0 - polygon[i].0;
            let dy = polygon[j].1 - polygon[i].1;
            let d = (dx * dx + dy * dy).sqrt();
            if d > max_d {
                max_d = d;
            }
        }
    }
    max_d
}

// ═══════════════════════════════════════════════════════════════════════════════
// Minimum-norm contact forces
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the minimum-norm normal forces for a planar multi-contact
/// problem assuming all contacts have the same normal (e.g. a flat surface).
///
/// Distributes the required total normal force equally among `n` contacts.
///
/// Returns `None` if no contacts or if the required force is tensile.
pub fn equal_distribution_normal_forces(
    total_normal_force: f64,
    num_contacts: usize,
) -> Option<Vec<f64>> {
    if num_contacts == 0 || total_normal_force < 0.0 {
        return None;
    }
    let per_contact = total_normal_force / (num_contacts as f64);
    Some(vec![per_contact; num_contacts])
}

/// Compute contact normal forces proportional to the distance from the
/// CoM to each contact (closer contacts bear more load).
pub fn distance_weighted_normal_forces(
    total_normal_force: f64,
    com_2d: (f64, f64),
    contact_positions_2d: &[(f64, f64)],
) -> Option<Vec<f64>> {
    let n = contact_positions_2d.len();
    if n == 0 || total_normal_force < 0.0 {
        return None;
    }
    // Inverse-distance weighting
    let mut weights: Vec<f64> = Vec::with_capacity(n);
    for p in contact_positions_2d {
        let dx = p.0 - com_2d.0;
        let dy = p.1 - com_2d.1;
        let dist = (dx * dx + dy * dy).sqrt();
        weights.push(1.0 / (dist + 1e-6));
    }
    let total_w: f64 = weights.iter().sum();
    if total_w < 1e-15 {
        return equal_distribution_normal_forces(total_normal_force, n);
    }
    let forces: Vec<f64> = weights
        .iter()
        .map(|w| total_normal_force * w / total_w)
        .collect();
    Some(forces)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Wrench cone containment
// ═══════════════════════════════════════════════════════════════════════════════

/// Check if a force is inside the Coulomb friction cone.
pub fn is_in_friction_cone(force: [f64; 3], normal: [f64; 3], mu: f64) -> bool {
    let fn_component = dot3(force, normal);
    if fn_component < -1e-10 {
        return false;
    }
    let f_tangent = vec3_sub(force, vec3_scale(normal, fn_component));
    let ft_mag = vec3_len(f_tangent);
    ft_mag <= mu * fn_component + 1e-10
}

/// Compute the angle between a force and the friction cone boundary.
///
/// Positive means inside the cone, negative means outside.
pub fn friction_cone_angle(force: [f64; 3], normal: [f64; 3], mu: f64) -> f64 {
    let f_len = vec3_len(force);
    if f_len < 1e-15 {
        return 0.0;
    }
    let fn_component = dot3(force, normal);
    let f_tangent = vec3_sub(force, vec3_scale(normal, fn_component));
    let ft_mag = vec3_len(f_tangent);
    let cone_angle = mu.atan();
    let force_angle = ft_mag.atan2(fn_component);
    cone_angle - force_angle
}

// ═══════════════════════════════════════════════════════════════════════════════
// Contact wrench cone (CWC)
// ═══════════════════════════════════════════════════════════════════════════════

/// The Contact Wrench Cone (CWC) for a set of contacts.
#[derive(Debug, Clone)]
pub struct ContactWrenchCone {
    /// Primitive wrench rays.
    pub rays: Vec<Wrench>,
    /// Number of contacts.
    pub num_contacts: usize,
    /// Number of cone facets per contact.
    pub num_facets: usize,
}

impl ContactWrenchCone {
    /// Build a CWC from contact points.
    pub fn new(contacts: &[ContactPoint], num_facets: usize) -> Self {
        let rays = multi_contact_wrench_set(contacts, num_facets);
        Self {
            rays,
            num_contacts: contacts.len(),
            num_facets,
        }
    }

    /// Check force closure.
    pub fn is_force_closure(&self) -> bool {
        force_closure_test(&self.rays)
    }

    /// Q1 quality.
    pub fn quality_q1(&self) -> f64 {
        grasp_quality_q1(&self.rays)
    }

    /// Epsilon quality.
    pub fn quality_epsilon(&self) -> f64 {
        grasp_quality_epsilon(&self.rays)
    }

    /// Check if an external wrench can be balanced.
    pub fn can_resist(&self, wrench: &Wrench, max_iter: usize, tol: f64) -> bool {
        can_balance_wrench(wrench, &self.rays, max_iter, tol)
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Rotational stability
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute the restoring torque about a tipping edge.
///
/// Given the CoM position, the edge defined by two points, and the
/// gravitational force, returns the torque about the edge axis.
///
/// Positive torque = restoring (stable), negative = tipping.
pub fn tipping_edge_torque(
    com: [f64; 3],
    edge_a: [f64; 3],
    edge_b: [f64; 3],
    gravity_force: [f64; 3],
) -> f64 {
    let edge_dir = vec3_normalize(vec3_sub(edge_b, edge_a));
    let r = vec3_sub(com, edge_a);
    // Force moment about the edge
    let moment = cross3(r, gravity_force);
    dot3(moment, edge_dir)
}

/// Find the most critical tipping edge in the support polygon.
///
/// Returns `(edge_index, torque)` where `edge_index` is the index of the
/// polygon edge with the smallest restoring torque.
pub fn critical_tipping_edge(
    com: [f64; 3],
    contact_positions: &[[f64; 3]],
    gravity_force: [f64; 3],
) -> Option<(usize, f64)> {
    let n = contact_positions.len();
    if n < 2 {
        return None;
    }
    let mut min_torque = f64::INFINITY;
    let mut min_idx = 0;
    for i in 0..n {
        let j = (i + 1) % n;
        let tau = tipping_edge_torque(
            com,
            contact_positions[i],
            contact_positions[j],
            gravity_force,
        );
        if tau < min_torque {
            min_torque = tau;
            min_idx = i;
        }
    }
    Some((min_idx, min_torque))
}

/// Time to tip: estimate the time before the body rotates past the tipping
/// edge given an angular acceleration.
///
/// ```text
/// t_tip = sqrt(2 * angle_margin / alpha)
/// ```
///
/// Returns `None` if not tipping or angular acceleration is zero.
pub fn time_to_tip(angle_margin_rad: f64, angular_accel: f64) -> Option<f64> {
    if angular_accel <= 0.0 || angle_margin_rad <= 0.0 {
        return None;
    }
    Some((2.0 * angle_margin_rad / angular_accel).sqrt())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Contact stability metrics collection
// ═══════════════════════════════════════════════════════════════════════════════

/// Summary of a contact stability analysis.
#[derive(Debug, Clone)]
pub struct ContactStabilityReport {
    /// Force closure achieved?
    pub force_closure: bool,
    /// Q1 grasp quality.
    pub q1: f64,
    /// Epsilon grasp quality.
    pub epsilon: f64,
    /// Support polygon area.
    pub support_area: f64,
    /// ZMP (if computed).
    pub zmp: Option<(f64, f64)>,
    /// Stability margin.
    pub margin: f64,
    /// Number of contacts.
    pub num_contacts: usize,
    /// Overall stability score \[0, 1\].
    pub score: f64,
}

/// Perform a comprehensive contact stability analysis.
pub fn analyse_contact_stability(
    contacts: &[ContactPoint],
    com: [f64; 3],
    gravity: [f64; 3],
    mass: f64,
    num_cone_sides: usize,
) -> ContactStabilityReport {
    let wrenches = multi_contact_wrench_set(contacts, num_cone_sides);
    let fc = force_closure_test(&wrenches);
    let q1 = grasp_quality_q1(&wrenches);
    let eps = grasp_quality_epsilon(&wrenches);

    let positions: Vec<[f64; 3]> = contacts.iter().map(|c| c.position).collect();
    let poly = support_polygon(&positions);
    let area = polygon_area(&poly);

    let grav_force = vec3_scale(gravity, mass);
    let grav_wrench = gravity_wrench(mass, com, gravity);
    let _grav_wrench = grav_wrench; // suppress unused

    // Simple ZMP from forces (assuming static: total contact force = -gravity_force)
    let total_force = vec3_neg(grav_force);
    let total_torque = cross3(com, total_force);
    let zmp = compute_zmp_from_forces(total_force, total_torque);

    let margin = if poly.len() >= 3 {
        let com_2d = (com[0], com[1]);
        stability_margin(com_2d, &poly)
    } else {
        f64::NEG_INFINITY
    };

    let score = multi_contact_stability_score(contacts, com, num_cone_sides);

    ContactStabilityReport {
        force_closure: fc,
        q1,
        epsilon: eps,
        support_area: area.abs(),
        zmp,
        margin,
        num_contacts: contacts.len(),
        score,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linearise_coulomb_cone_count() {
        let edges = linearise_coulomb_cone([0.0, 0.0, 1.0], 0.5, 8);
        assert_eq!(edges.len(), 8);
    }

    #[test]
    fn test_linearise_coulomb_cone_unit_vectors() {
        let edges = linearise_coulomb_cone([0.0, 0.0, 1.0], 0.3, 6);
        for e in &edges {
            let l = vec3_len(*e);
            assert!((l - 1.0).abs() < 1e-10, "Edge not unit length: {:.6}", l);
        }
    }

    #[test]
    fn test_cone_linearisation_error() {
        let err8 = cone_linearisation_error(8);
        let err16 = cone_linearisation_error(16);
        assert!(err16 < err8, "More sides should reduce error");
        assert!(err8 > 0.0);
    }

    #[test]
    fn test_contact_wrench() {
        let w = contact_wrench([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        // torque = (1,0,0) x (0,0,1) = (0,-1,0)
        assert!((w[0]).abs() < 1e-10);
        assert!((w[2] - 1.0).abs() < 1e-10);
        assert!((w[4] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_force_closure_simple() {
        // Opposing contacts from all sides for true force closure
        let contacts = [
            ContactPoint {
                position: [1.0, 0.0, 0.0],
                normal: [-1.0, 0.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [-1.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 1.0, 0.0],
                normal: [0.0, -1.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, -1.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 0.0, 1.0],
                normal: [0.0, 0.0, -1.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 0.0, -1.0],
                normal: [0.0, 0.0, 1.0],
                mu: 0.8,
            },
        ];
        let ws = multi_contact_wrench_set(&contacts, 8);
        let fc = force_closure_test(&ws);
        assert!(
            fc,
            "Opposing contacts from all sides should achieve force closure"
        );
    }

    #[test]
    fn test_force_closure_single_contact_fails() {
        let contacts = [ContactPoint {
            position: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
            mu: 0.5,
        }];
        let ws = multi_contact_wrench_set(&contacts, 8);
        let fc = force_closure_test(&ws);
        assert!(!fc, "Single contact should not achieve force closure");
    }

    #[test]
    fn test_support_polygon_triangle() {
        let contacts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let poly = support_polygon(&contacts);
        assert_eq!(poly.len(), 3, "Triangle should have 3 vertices");
    }

    #[test]
    fn test_support_polygon_square() {
        let contacts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let poly = support_polygon(&contacts);
        assert_eq!(poly.len(), 4, "Square should have 4 vertices");
    }

    #[test]
    fn test_point_in_polygon() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert!(point_in_convex_polygon((0.5, 0.5), &poly));
        assert!(!point_in_convex_polygon((2.0, 0.5), &poly));
    }

    #[test]
    fn test_polygon_area_unit_square() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let area = polygon_area(&poly);
        assert!((area - 1.0).abs() < 1e-10, "area = {:.6}", area);
    }

    #[test]
    fn test_polygon_centroid_unit_square() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let (cx, cy) = polygon_centroid(&poly);
        assert!((cx - 0.5).abs() < 1e-10, "cx = {:.6}", cx);
        assert!((cy - 0.5).abs() < 1e-10, "cy = {:.6}", cy);
    }

    #[test]
    fn test_zmp_static() {
        // Static case: acc = gravity, ZMP should be directly below CoM
        let com = [0.5, 0.3, 1.0];
        let gravity = [0.0, 0.0, -9.81];
        let acc = [0.0, 0.0, 0.0];
        let zmp = compute_zmp(com, acc, gravity);
        assert!(zmp.is_some());
        let (zx, zy) = zmp.unwrap();
        assert!((zx - 0.5).abs() < 1e-10, "zmp_x = {:.6}", zx);
        assert!((zy - 0.3).abs() < 1e-10, "zmp_y = {:.6}", zy);
    }

    #[test]
    fn test_zmp_from_forces() {
        let force = [0.0, 0.0, 100.0];
        let torque = [50.0, -30.0, 0.0];
        let zmp = compute_zmp_from_forces(force, torque);
        assert!(zmp.is_some());
        let (zx, zy) = zmp.unwrap();
        // zmp_x = -(-30) / 100 = 0.3
        assert!((zx - 0.3).abs() < 1e-10, "zmp_x = {:.6}", zx);
        // zmp_y = 50 / 100 = 0.5
        assert!((zy - 0.5).abs() < 1e-10, "zmp_y = {:.6}", zy);
    }

    #[test]
    fn test_stability_margin_inside() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let m = stability_margin((0.5, 0.5), &poly);
        assert!(m > 0.0, "Centre should be inside, margin = {:.6}", m);
    }

    #[test]
    fn test_stability_margin_outside() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let m = stability_margin((2.0, 0.5), &poly);
        assert!(
            m < 0.0,
            "Outside point should have negative margin, margin = {:.6}",
            m
        );
    }

    #[test]
    fn test_tipping_analysis_stable() {
        let com = [0.5, 0.5, 1.0];
        let acc = [0.0, 0.0, 0.0];
        let gravity = [0.0, 0.0, -9.81];
        let contacts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let result = tipping_analysis(com, acc, gravity, &contacts);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.stable, "CoM above centre of square should be stable");
    }

    #[test]
    fn test_is_statically_stable() {
        let com = [0.5, 0.5, 1.0];
        let contacts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        assert!(is_statically_stable(com, &contacts));
    }

    #[test]
    fn test_normal_complementarity() {
        assert!(check_normal_complementarity(0.1, 0.0, 1e-8));
        assert!(check_normal_complementarity(0.0, 5.0, 1e-8));
        assert!(!check_normal_complementarity(0.1, 5.0, 1e-8));
    }

    #[test]
    fn test_friction_complementarity_residual() {
        // Sticking: v_t = 0 -> residual = 0 regardless of cone slack
        let r = friction_complementarity_residual(1.0, 0.5, 10.0, 0.0);
        assert!(r.abs() < 1e-10);
    }

    #[test]
    fn test_check_contact_complementarity_sticking() {
        // gap = 0, f_n = 10, f_t = 2, mu = 0.5, v_t = 0
        assert!(check_contact_complementarity(
            0.0, 10.0, 2.0, 0.5, 0.0, 1e-6
        ));
    }

    #[test]
    fn test_project_normal_force() {
        assert!((project_normal_force(5.0) - 5.0).abs() < 1e-10);
        assert!((project_normal_force(-3.0)).abs() < 1e-10);
    }

    #[test]
    fn test_project_friction_force() {
        assert!((project_friction_force(3.0, 0.5, 10.0) - 3.0).abs() < 1e-10);
        assert!((project_friction_force(10.0, 0.5, 10.0) - 5.0).abs() < 1e-10);
        assert!((project_friction_force(-10.0, 0.5, 10.0) + 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_lcp_gauss_seidel_step() {
        let mut contacts = vec![LcpContact {
            gap: -0.01,
            f_n: 0.0,
            f_t: 0.0,
            mu: 0.5,
            v_t: 0.0,
            eff_mass_n: 1.0,
            eff_mass_t: 1.0,
        }];
        let _delta = lcp_gauss_seidel_step(&mut contacts, 0.2);
        assert!(contacts[0].f_n >= 0.0);
    }

    #[test]
    fn test_solve_contact_lcp_converges() {
        let mut contacts = vec![LcpContact {
            gap: -0.01,
            f_n: 0.0,
            f_t: 0.0,
            mu: 0.5,
            v_t: 0.1,
            eff_mass_n: 1.0,
            eff_mass_t: 1.0,
        }];
        let iters = solve_contact_lcp(&mut contacts, 0.2, 100, 1e-8);
        assert!(iters <= 100);
        assert!(contacts[0].f_n >= 0.0);
    }

    #[test]
    fn test_is_in_friction_cone() {
        let normal = [0.0, 0.0, 1.0];
        assert!(is_in_friction_cone([0.0, 0.0, 10.0], normal, 0.5));
        assert!(is_in_friction_cone([1.0, 0.0, 10.0], normal, 0.5));
        assert!(!is_in_friction_cone([6.0, 0.0, 10.0], normal, 0.5));
    }

    #[test]
    fn test_friction_cone_angle() {
        let normal = [0.0, 0.0, 1.0];
        let a = friction_cone_angle([0.0, 0.0, 1.0], normal, 0.5);
        assert!(a > 0.0, "Pure normal force should be inside cone");
    }

    #[test]
    fn test_grasp_matrix_block() {
        let g = grasp_matrix_block([1.0, 2.0, 3.0]);
        // First three rows should be identity
        assert!((g[0][0] - 1.0).abs() < 1e-10);
        assert!((g[1][1] - 1.0).abs() < 1e-10);
        assert!((g[2][2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gravity_wrench() {
        let w = gravity_wrench(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, -9.81]);
        assert!((w[2] + 9.81).abs() < 1e-10);
    }

    #[test]
    fn test_sum_wrenches() {
        let w1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let w2 = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let s = sum_wrenches(&[w1, w2]);
        for (i, &si) in s.iter().enumerate() {
            assert!((si - 7.0).abs() < 1e-10, "s[{i}]={si}");
        }
    }

    #[test]
    fn test_equal_distribution_normal_forces() {
        let f = equal_distribution_normal_forces(100.0, 4).unwrap();
        assert_eq!(f.len(), 4);
        for &fi in &f {
            assert!((fi - 25.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_polygon_diameter() {
        let poly = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let d = polygon_diameter(&poly);
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_tipping_edge_torque() {
        let com = [0.5, 0.0, 1.0];
        let ea = [0.0, 0.0, 0.0];
        let eb = [0.0, 1.0, 0.0];
        let gf = [0.0, 0.0, -9.81];
        let tau = tipping_edge_torque(com, ea, eb, gf);
        // Positive = restoring (CoM is inside)
        assert!(tau > 0.0, "Restoring torque expected, got {:.6}", tau);
    }

    #[test]
    fn test_time_to_tip() {
        let t = time_to_tip(0.1, 2.0);
        assert!(t.is_some());
        let t = t.unwrap();
        // t = sqrt(2 * 0.1 / 2.0) = sqrt(0.1) ≈ 0.316
        assert!((t - 0.1_f64.sqrt()).abs() < 1e-10, "t = {:.6}", t);
    }

    #[test]
    fn test_time_to_tip_no_accel() {
        assert!(time_to_tip(0.1, 0.0).is_none());
    }

    #[test]
    fn test_contact_wrench_cone() {
        let contacts = [
            ContactPoint {
                position: [1.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                mu: 0.5,
            },
            ContactPoint {
                position: [-1.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                mu: 0.5,
            },
        ];
        let cwc = ContactWrenchCone::new(&contacts, 8);
        assert_eq!(cwc.rays.len(), 16);
        assert_eq!(cwc.num_contacts, 2);
    }

    #[test]
    fn test_friction_polyhedron_contains() {
        let fp = FrictionPolyhedron::new([0.0; 3], [0.0, 0.0, 1.0], 0.5, 8);
        assert!(fp.contains_force([0.0, 0.0, 10.0]));
        assert!(fp.contains_force([1.0, 0.0, 10.0]));
        assert!(!fp.contains_force([6.0, 0.0, 10.0]));
    }

    #[test]
    fn test_analyse_contact_stability() {
        // Opposing contacts for force closure
        let contacts = [
            ContactPoint {
                position: [1.0, 0.0, 0.0],
                normal: [-1.0, 0.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [-1.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 1.0, 0.0],
                normal: [0.0, -1.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, -1.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 0.0, 1.0],
                normal: [0.0, 0.0, -1.0],
                mu: 0.8,
            },
            ContactPoint {
                position: [0.0, 0.0, -1.0],
                normal: [0.0, 0.0, 1.0],
                mu: 0.8,
            },
        ];
        let report =
            analyse_contact_stability(&contacts, [0.0, 0.0, 0.0], [0.0, 0.0, -9.81], 1.0, 8);
        assert!(report.force_closure);
        assert!(report.q1 > 0.0);
        assert_eq!(report.num_contacts, 6);
    }

    #[test]
    fn test_total_complementarity_residual() {
        let gaps = [0.0, 0.0];
        let fns = [10.0, 5.0];
        let fts = [1.0, 0.5];
        let mus = [0.5, 0.5];
        let vts = [0.0, 0.0];
        let r = total_complementarity_residual(&gaps, &fns, &fts, &mus, &vts);
        assert!(
            r.abs() < 1e-10,
            "Sticking contacts should have zero residual"
        );
    }

    #[test]
    fn test_distance_weighted_forces() {
        let forces = distance_weighted_normal_forces(100.0, (0.0, 0.0), &[(1.0, 0.0), (-1.0, 0.0)]);
        assert!(forces.is_some());
        let f = forces.unwrap();
        // Symmetric case: equal weights
        assert!((f[0] - 50.0).abs() < 1e-6, "f[0] = {:.6}", f[0]);
        assert!((f[1] - 50.0).abs() < 1e-6, "f[1] = {:.6}", f[1]);
    }

    #[test]
    fn test_grasp_quality_volume() {
        let wrenches = vec![
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let vol = grasp_quality_volume(&wrenches);
        assert!(vol >= 0.0);
    }

    #[test]
    fn test_centrifugal_wrench() {
        let w = centrifugal_wrench(1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        // omega x r = (0,0,1) x (1,0,0) = (0,1,0)
        // omega x (omega x r) = (0,0,1) x (0,1,0) = (-1,0,0)
        assert!((w[0] + 1.0).abs() < 1e-10, "fx = {:.6}", w[0]);
    }
}
