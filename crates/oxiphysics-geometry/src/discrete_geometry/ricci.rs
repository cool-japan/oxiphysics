// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Discrete conformal deformation of the metric of a triangle mesh (discrete Ricci
//! flow); the method is described in the documentation of the parent module.

use super::DiscreteMesh;
use super::support::{Components, EdgeTable, ScaledSpd};
use crate::error::{Error, Result};
use std::f64::consts::PI;

/// Largest number of step halvings in the line search of one Newton iteration.
const MAX_HALVINGS: usize = 40;
/// Armijo constant of the line search on the 2-norm of the curvature error.
const SUFFICIENT_DECREASE: f64 = 1e-4;
/// Relative accuracy to which the target curvature must satisfy Gauss-Bonnet.
const GAUSS_BONNET_TOLERANCE: f64 = 1e-9;
/// A triangle whose quadrupled area is below this fraction of its squared longest
/// side contributes nothing to the Jacobian (its angles are locally constant).
const DEGENERATE_AREA: f64 = 1e-12;

/// Outcome of a successful [`ricci_flow`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RicciFlowSummary {
    /// Number of Newton iterations performed (zero if the input already satisfied the
    /// tolerance).
    pub iterations: usize,
    /// Final curvature error `max_i |K_i - Kbar_i|`.
    pub residual: f64,
}

/// Angles of a triangle given by its three side lengths.
struct Corners {
    /// Interior angle at corner `k`, the corner opposite side `k`.
    angles: [f64; 3],
    /// `2 l_a l_b cos(angle k) = l_a^2 + l_b^2 - l_k^2` for the two sides `a`, `b`
    /// adjacent to corner `k`, in units of the squared longest side.
    dots: [f64; 3],
    /// Four times the area in units of the squared longest side; zero if the
    /// triangle inequality fails.
    four_area: f64,
}

/// Angles from the side lengths `l[k]` of the sides opposite the corners `k`; the
/// lengths must be finite and positive.
///
/// The sides are first divided by the longest one (the angles do not depend on the
/// scale), which keeps every product below in range however large or small the
/// lengths are.  The area comes from Kahan's numerically stable form of Heron's
/// formula and the angles from `atan2(4 A, l_a^2 + l_b^2 - l_k^2)`, which stays
/// accurate for tiny angles and returns the degenerate angles `0` and `pi` when the
/// triangle inequality fails.
fn triangle_corners(lengths: [f64; 3]) -> Corners {
    let longest = lengths[0].max(lengths[1]).max(lengths[2]);
    let l: [f64; 3] = std::array::from_fn(|k| lengths[k] / longest);
    let sq = [l[0] * l[0], l[1] * l[1], l[2] * l[2]];
    let mut sorted = l;
    sorted.sort_unstable_by(|x, y| y.total_cmp(x));
    let (a, b, c) = (sorted[0], sorted[1], sorted[2]);
    let product = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
    let four_area = if product > 0.0 { product.sqrt() } else { 0.0 };
    let dots: [f64; 3] = std::array::from_fn(|k| sq[(k + 1) % 3] + sq[(k + 2) % 3] - sq[k]);
    Corners {
        angles: std::array::from_fn(|k| four_area.atan2(dots[k])),
        dots,
        four_area,
    }
}

/// A mesh prepared for the conformal metric computations.
struct RicciSurface {
    triangles: Vec<[usize; 3]>,
    /// Side lengths of every triangle in the metric of the mesh, opposite corners 0, 1, 2.
    base: Vec<[f64; 3]>,
    /// The angle sum that a vertex has when it is flat: `2 pi` inside, `pi` on the
    /// boundary, zero for vertices that belong to no triangle.
    flat_angle: Vec<f64>,
    /// Vertices that belong to at least one triangle.
    referenced: Vec<bool>,
    components: Components,
    /// Total curvature every component must have (Gauss-Bonnet), indexed like
    /// `components`; zero for components without triangles.
    total_curvature: Vec<f64>,
    /// The first vertex of every component with triangles; pinned in the Newton system.
    pinned: Vec<bool>,
}

impl RicciSurface {
    /// # Errors
    /// * [`Error::InvalidMesh`] for a non-finite vertex, a triangle that references a
    ///   missing vertex, or an edge shared by more than two triangles;
    /// * [`Error::DegenerateGeometry`] for a triangle with a side of zero length.
    fn new(mesh: &DiscreteMesh) -> Result<Self> {
        let n = mesh.num_vertices();
        for (i, v) in mesh.vertices.iter().enumerate() {
            if !v.iter().all(|c| c.is_finite()) {
                return Err(Error::invalid_mesh(format!(
                    "vertex {i} has a non-finite coordinate"
                )));
            }
        }
        let table = EdgeTable::build(n, &mesh.triangles)?;
        let mut base = Vec::with_capacity(mesh.triangles.len());
        for (f, tri) in mesh.triangles.iter().enumerate() {
            let side = |k: usize| {
                let (a, b) = (
                    mesh.vertices[tri[(k + 1) % 3]],
                    mesh.vertices[tri[(k + 2) % 3]],
                );
                super::dist3(a, b)
            };
            let lengths = [side(0), side(1), side(2)];
            if !lengths.iter().all(|l| l.is_finite() && *l > 0.0) {
                return Err(Error::degenerate_geometry(format!(
                    "triangle {f} has a side of zero length"
                )));
            }
            base.push(lengths);
        }
        let boundary = table.boundary_vertices(n);
        let mut referenced = vec![false; n];
        for tri in &mesh.triangles {
            for &v in tri {
                referenced[v] = true;
            }
        }
        let flat_angle: Vec<f64> = referenced
            .iter()
            .zip(boundary.iter())
            .map(|(&used, &on_boundary)| match (used, on_boundary) {
                (false, _) => 0.0,
                (true, true) => PI,
                (true, false) => 2.0 * PI,
            })
            .collect();
        let components = Components::label(n, &mesh.triangles);
        let mut total_curvature = vec![0.0_f64; components.count];
        for (v, &angle) in flat_angle.iter().enumerate() {
            total_curvature[components.id[v]] += angle;
        }
        for tri in &mesh.triangles {
            total_curvature[components.id[tri[0]]] -= PI;
        }
        let mut pinned = vec![false; n];
        let mut seen = vec![false; components.count];
        for (v, &used) in referenced.iter().enumerate() {
            let c = components.id[v];
            if used && !seen[c] {
                seen[c] = true;
                pinned[v] = true;
            }
        }
        Ok(Self {
            triangles: mesh.triangles.clone(),
            base,
            flat_angle,
            referenced,
            components,
            total_curvature,
            pinned,
        })
    }

    fn vertex_count(&self) -> usize {
        self.flat_angle.len()
    }

    /// Check the conformal factors: right length, finite.
    fn check_factors(&self, u: &[f64]) -> Result<()> {
        if u.len() != self.vertex_count() {
            return Err(Error::length_mismatch(self.vertex_count(), u.len()));
        }
        if let Some(bad) = u.iter().position(|x| !x.is_finite()) {
            return Err(Error::general(format!(
                "conformal factor {bad} is not finite"
            )));
        }
        Ok(())
    }

    /// Check the target curvature: right length, finite, and satisfying the
    /// Gauss-Bonnet theorem on every connected component.
    fn check_target(&self, k_target: &[f64]) -> Result<()> {
        if k_target.len() != self.vertex_count() {
            return Err(Error::length_mismatch(self.vertex_count(), k_target.len()));
        }
        if let Some(bad) = k_target.iter().position(|x| !x.is_finite()) {
            return Err(Error::general(format!(
                "target curvature {bad} is not finite"
            )));
        }
        let mut sum = vec![0.0_f64; self.components.count];
        let mut magnitude = vec![0.0_f64; self.components.count];
        for (v, &used) in self.referenced.iter().enumerate() {
            if used {
                sum[self.components.id[v]] += k_target[v];
                magnitude[self.components.id[v]] += k_target[v].abs();
            }
        }
        let mut has_triangles = vec![false; self.components.count];
        for tri in &self.triangles {
            has_triangles[self.components.id[tri[0]]] = true;
        }
        for c in 0..self.components.count {
            if !has_triangles[c] {
                continue;
            }
            let required = self.total_curvature[c];
            let tolerance = GAUSS_BONNET_TOLERANCE * (1.0 + magnitude[c] + required.abs());
            if (sum[c] - required).abs() > tolerance {
                return Err(Error::out_of_range(
                    "total target curvature of a connected component (Gauss-Bonnet)",
                    sum[c],
                    required - tolerance,
                    required + tolerance,
                ));
            }
        }
        Ok(())
    }

    /// `exp(u_i / 2)` for every vertex; the edge `ij` is scaled by the product of the
    /// two end values.
    fn edge_scales(&self, u: &[f64]) -> Result<Vec<f64>> {
        let scales: Vec<f64> = u.iter().map(|x| (0.5 * x).exp()).collect();
        if scales.iter().any(|s| !s.is_finite() || *s == 0.0) {
            return Err(Error::general(
                "the conformal factors over- or underflow the edge lengths",
            ));
        }
        Ok(scales)
    }

    /// Side lengths of triangle `f` in the deformed metric.
    fn sides(&self, f: usize, scales: &[f64]) -> [f64; 3] {
        let tri = &self.triangles[f];
        std::array::from_fn(|k| {
            self.base[f][k] * scales[tri[(k + 1) % 3]] * scales[tri[(k + 2) % 3]]
        })
    }

    /// [`sides`](Self::sides), failing if a length over- or underflows.
    fn checked_sides(&self, f: usize, scales: &[f64]) -> Result<[f64; 3]> {
        let sides = self.sides(f, scales);
        if sides.iter().all(|l| l.is_finite() && *l > 0.0) {
            Ok(sides)
        } else {
            Err(Error::general(format!(
                "the conformal factors over- or underflow the edge lengths of triangle {f}"
            )))
        }
    }

    /// Integrated curvature of the deformed metric at every vertex.
    fn curvature(&self, u: &[f64]) -> Result<Vec<f64>> {
        let scales = self.edge_scales(u)?;
        let mut curvature = self.flat_angle.clone();
        for (f, tri) in self.triangles.iter().enumerate() {
            let corners = triangle_corners(self.checked_sides(f, &scales)?);
            for (&v, angle) in tri.iter().zip(corners.angles.iter()) {
                curvature[v] -= angle;
            }
        }
        Ok(curvature)
    }

    /// `Kbar - K` at the vertices that belong to a triangle, zero elsewhere.
    fn error(&self, k_target: &[f64], k: &[f64]) -> Vec<f64> {
        self.referenced
            .iter()
            .zip(k_target.iter().zip(k.iter()))
            .map(|(&used, (t, c))| if used { t - c } else { 0.0 })
            .collect()
    }

    /// The Jacobian `dK/du` of the curvature as `(row, column, value)` entries: the
    /// cotangent Laplacian of the deformed metric.  Triangles that violate the
    /// triangle inequality, or are nearly flat, do not contribute: their angles are
    /// (locally) constant.
    fn jacobian_entries(&self, u: &[f64]) -> Result<Vec<(usize, usize, f64)>> {
        let scales = self.edge_scales(u)?;
        let mut entries = Vec::with_capacity(12 * self.triangles.len());
        for (f, tri) in self.triangles.iter().enumerate() {
            let corners = triangle_corners(self.checked_sides(f, &scales)?);
            if corners.four_area <= DEGENERATE_AREA {
                continue;
            }
            for (k, dot) in corners.dots.iter().enumerate() {
                // The edge opposite corner k carries the weight cot(angle k) / 2.
                let weight = 0.5 * dot / corners.four_area;
                let (a, b) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
                entries.push((a, a, weight));
                entries.push((b, b, weight));
                entries.push((a, b, -weight));
                entries.push((b, a, -weight));
            }
        }
        Ok(entries)
    }

    /// The Newton step `delta` with `J delta = error` (the Jacobian is singular only
    /// by the constants of the components, which are pinned).
    fn newton_direction(&self, u: &[f64], error: &[f64]) -> Result<Vec<f64>> {
        let entries = self.jacobian_entries(u)?;
        // The error sums to zero on every component by Gauss-Bonnet; remove the
        // rounding residue so that the pinned system is consistent.
        let mut sum = vec![0.0_f64; self.components.count];
        let mut count = vec![0usize; self.components.count];
        for (v, &used) in self.referenced.iter().enumerate() {
            if used {
                sum[self.components.id[v]] += error[v];
                count[self.components.id[v]] += 1;
            }
        }
        let rhs: Vec<f64> = error
            .iter()
            .enumerate()
            .map(|(v, &e)| {
                let c = self.components.id[v];
                if self.referenced[v] {
                    e - sum[c] / count[c] as f64
                } else {
                    0.0
                }
            })
            .collect();
        ScaledSpd::new(self.vertex_count(), &entries, &self.pinned)
            .solve(&rhs, "discrete Ricci flow: Newton step")
    }
}

/// Largest absolute value and Euclidean norm of a vector.
fn norms(values: &[f64]) -> (f64, f64) {
    let largest = values.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let euclid = values.iter().map(|x| x * x).sum::<f64>().sqrt();
    (largest, euclid)
}

/// Integrated curvature of the metric `|p_i - p_j| exp((u_i + u_j) / 2)` at every
/// vertex: `2 pi` (`pi` on the boundary) minus the sum of the corner angles.
///
/// With `u = 0` this is the angle defect of the mesh itself.  Vertices that belong
/// to no triangle get zero.  The curvature sums to `2 pi chi` over a connected mesh
/// for every `u` (Gauss-Bonnet); a uniform shift of `u` rescales the metric and
/// leaves the curvature unchanged.
///
/// # Errors
/// * [`Error::LengthMismatch`] if `u` does not have one entry per vertex;
/// * [`Error::InvalidMesh`] / [`Error::DegenerateGeometry`] for a malformed mesh (see
///   [`ricci_flow`]);
/// * [`Error::General`] for non-finite factors or factors so large that the edge
///   lengths overflow.
pub fn ricci_curvature(mesh: &DiscreteMesh, u: &[f64]) -> Result<Vec<f64>> {
    let surface = RicciSurface::new(mesh)?;
    surface.check_factors(u)?;
    surface.curvature(u)
}

/// Lengths of the edges of the deformed metric `|p_i - p_j| exp((u_i + u_j) / 2)`,
/// as `(low, high, length)` for every edge `(low, high)` with `low < high`, sorted
/// lexicographically.
///
/// # Errors
/// As for [`ricci_curvature`].
pub fn conformal_edge_lengths(mesh: &DiscreteMesh, u: &[f64]) -> Result<Vec<(usize, usize, f64)>> {
    let surface = RicciSurface::new(mesh)?;
    surface.check_factors(u)?;
    let scales = surface.edge_scales(u)?;
    let table = EdgeTable::build(mesh.num_vertices(), &mesh.triangles)?;
    Ok(table
        .edges
        .iter()
        .map(|&(a, b)| {
            let length = super::dist3(mesh.vertices[a], mesh.vertices[b]) * scales[a] * scales[b];
            (a, b, length)
        })
        .collect())
}

/// One explicit step of the discrete Ricci flow `du_i/dt = Kbar_i - K_i`.
///
/// The curvature `K` of the metric `|p_i - p_j| exp((u_i + u_j) / 2)` is evaluated
/// at the incoming `u` (see [`ricci_curvature`]) and then
/// `u_i += step * (k_target_i - K_i)` is applied to every vertex that belongs to a
/// triangle.  A vertex whose curvature is above its target is shrunk, which
/// increases its corner angles and lowers its curvature.  Returns the largest
/// curvature error `max_i |K_i - k_target_i|` **before** the update.
///
/// The explicit flow is stable for steps below `2 / lambda_max`, where `lambda_max`
/// is the largest eigenvalue of the cotangent Laplacian of the current metric.  It
/// is at most twice the largest diagonal entry of that matrix, which is `2 * 6 *
/// cot(60 degrees) = 6.9` for a regular mesh of valence six, so steps below `0.29`
/// are stable there.  Use [`ricci_flow`] to solve for the target curvature
/// reliably.
///
/// `k_target` must satisfy the Gauss-Bonnet theorem on every connected component of
/// the mesh: its sum over the vertices of a component equals `2 pi chi` (the sum of
/// `2 pi` per interior and `pi` per boundary vertex, minus `pi` per triangle).
///
/// # Errors
/// * [`Error::LengthMismatch`] if `u` or `k_target` does not have one entry per
///   vertex;
/// * [`Error::OutOfRange`] if `step` is not positive and finite, or if `k_target`
///   violates Gauss-Bonnet on a component;
/// * [`Error::General`] for non-finite values;
/// * [`Error::InvalidMesh`] / [`Error::DegenerateGeometry`] for a non-finite vertex,
///   a triangle that references a missing vertex, an edge shared by more than two
///   triangles, or a triangle with a zero-length side.
pub fn ricci_flow_step(
    mesh: &DiscreteMesh,
    u: &mut [f64],
    k_target: &[f64],
    step: f64,
) -> Result<f64> {
    if !(step.is_finite() && step > 0.0) {
        return Err(Error::out_of_range(
            "step",
            step,
            f64::MIN_POSITIVE,
            f64::MAX,
        ));
    }
    let surface = RicciSurface::new(mesh)?;
    surface.check_factors(u)?;
    surface.check_target(k_target)?;
    let k = surface.curvature(u)?;
    let error = surface.error(k_target, &k);
    let (largest, _) = norms(&error);
    for (ui, e) in u.iter_mut().zip(error.iter()) {
        *ui += step * e;
    }
    Ok(largest)
}

/// Drive the conformal factors `u` of the metric `|p_i - p_j| exp((u_i + u_j) / 2)`
/// to the target curvature: iterate until `max_i |K_i - k_target_i| <= tol`.
///
/// This solves the steady state of the Ricci flow of [`ricci_flow_step`] with a
/// damped Newton iteration (Gu, Luo, Sun & Wu 2018).  The Jacobian of the curvature
/// is the cotangent Laplacian of the current metric; it is singular exactly by the
/// uniform rescaling of each connected component, which is fixed by keeping the
/// factor of one vertex per component.  The Newton step is halved until the
/// Euclidean norm of the curvature error decreases.  `u` holds the initial factors
/// and receives the result, which is determined up to the added constant of each
/// connected component.
///
/// The solution of `K(u) = k_target` need not exist for the triangulation at hand
/// (a deformed triangle may degenerate, which for large changes of curvature calls
/// for edge flips); the call then fails with [`Error::ConvergenceFailure`] and `u`
/// holds the last accepted iterate.
///
/// # Errors
/// * [`Error::OutOfRange`] if `tol` is not positive and finite, or if `k_target`
///   violates Gauss-Bonnet on a connected component (see [`ricci_flow_step`]);
/// * [`Error::ConvergenceFailure`] if the tolerance is not reached within `max_iter`
///   iterations or the line search cannot reduce the error;
/// * the errors of [`ricci_flow_step`] for malformed input.
pub fn ricci_flow(
    mesh: &DiscreteMesh,
    u: &mut [f64],
    k_target: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<RicciFlowSummary> {
    if !(tol.is_finite() && tol > 0.0) {
        return Err(Error::out_of_range("tol", tol, f64::MIN_POSITIVE, f64::MAX));
    }
    let surface = RicciSurface::new(mesh)?;
    surface.check_factors(u)?;
    surface.check_target(k_target)?;
    let mut error = surface.error(k_target, &surface.curvature(u)?);
    let (mut largest, mut euclid) = norms(&error);
    for iteration in 0..max_iter {
        if largest <= tol {
            return Ok(RicciFlowSummary {
                iterations: iteration,
                residual: largest,
            });
        }
        let direction = surface.newton_direction(u, &error)?;
        let mut alpha = 1.0_f64;
        let mut accepted = None;
        for _ in 0..MAX_HALVINGS {
            let trial: Vec<f64> = u
                .iter()
                .zip(direction.iter())
                .map(|(x, d)| x + alpha * d)
                .collect();
            if let Ok(k) = surface.curvature(&trial) {
                let trial_error = surface.error(k_target, &k);
                let (trial_largest, trial_euclid) = norms(&trial_error);
                if trial_euclid <= (1.0 - SUFFICIENT_DECREASE * alpha) * euclid {
                    accepted = Some((trial, trial_error, trial_largest, trial_euclid));
                    break;
                }
            }
            alpha *= 0.5;
        }
        let Some((trial, trial_error, trial_largest, trial_euclid)) = accepted else {
            return Err(Error::convergence_failure(
                "discrete Ricci flow",
                iteration,
                largest,
            ));
        };
        u.copy_from_slice(&trial);
        error = trial_error;
        largest = trial_largest;
        euclid = trial_euclid;
    }
    if largest <= tol {
        return Ok(RicciFlowSummary {
            iterations: max_iter,
            residual: largest,
        });
    }
    Err(Error::convergence_failure(
        "discrete Ricci flow",
        max_iter,
        largest,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discrete_geometry::test_support::{
        XorShift, discrete, open_cylinder, spherical_cap, torus,
    };
    use crate::discrete_geometry::{unit_icosahedron, unit_tetrahedron};
    use crate::test_meshes::{flat_grid, icosphere};

    fn perturbed_icosahedron(seed: u64, amount: f64) -> DiscreteMesh {
        let mut m = unit_icosahedron();
        let mut rng = XorShift::new(seed);
        for p in &mut m.vertices {
            let s = 1.0 + amount * rng.signed();
            *p = [p[0] * s, p[1] * s, p[2] * s];
        }
        m
    }

    fn random_factors(n: usize, amplitude: f64, seed: u64) -> Vec<f64> {
        let mut rng = XorShift::new(seed);
        (0..n).map(|_| amplitude * rng.signed()).collect()
    }

    fn sum(values: &[f64]) -> f64 {
        values.iter().sum()
    }

    fn largest_error(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
    }

    #[test]
    fn corner_angles_follow_the_law_of_cosines() {
        let equilateral = triangle_corners([1.0, 1.0, 1.0]);
        for angle in equilateral.angles {
            assert!((angle - PI / 3.0).abs() < 1e-14);
        }
        // Sides 3, 4, 5 opposite corners 0, 1, 2: the right angle is at corner 2.
        let right = triangle_corners([3.0, 4.0, 5.0]);
        assert!((right.angles[2] - PI / 2.0).abs() < 1e-14);
        assert!((right.angles[0] - (3.0_f64 / 4.0).atan()).abs() < 1e-14);
        assert!((right.angles[1] - (4.0_f64 / 3.0).atan()).abs() < 1e-14);
        // 4 A = 24, in units of the squared longest side 25.
        assert!((right.four_area - 24.0 / 25.0).abs() < 1e-14);
        let mut rng = XorShift::new(2);
        for _ in 0..200 {
            let l = [0.5 + rng.unit(), 0.5 + rng.unit(), 0.5 + rng.unit()];
            let corners = triangle_corners(l);
            assert!((sum(&corners.angles) - PI).abs() < 1e-13, "{l:?}");
            // The longest side is opposite the largest angle.
            let longest = (0..3).max_by(|&i, &j| l[i].total_cmp(&l[j])).unwrap_or(0);
            let widest = (0..3)
                .max_by(|&i, &j| corners.angles[i].total_cmp(&corners.angles[j]))
                .unwrap_or(0);
            assert_eq!(longest, widest);
        }
    }

    #[test]
    fn violated_triangle_inequality_gives_the_degenerate_angles() {
        // Side 3 is longer than 1 + 1: the angle opposite it is pi, the others vanish.
        let too_long = triangle_corners([1.0, 1.0, 3.0]);
        assert_eq!(too_long.four_area, 0.0);
        assert_eq!(too_long.angles[0], 0.0);
        assert_eq!(too_long.angles[1], 0.0);
        assert!((too_long.angles[2] - PI).abs() < 1e-15);
        // Exactly flat.
        let flat = triangle_corners([1.0, 1.0, 2.0]);
        assert_eq!(flat.angles[0], 0.0);
        assert!((flat.angles[2] - PI).abs() < 1e-15);
        // Almost flat but valid: tiny angles are still accurate.
        let sliver = triangle_corners([1.0, 1.0, 2.0 - 1e-12]);
        assert!(sliver.angles[0] > 0.0 && sliver.angles[0] < 1e-5);
        assert!((sum(&sliver.angles) - PI).abs() < 1e-14);
    }

    #[test]
    fn curvature_of_the_mesh_metric_is_its_angle_defect() {
        let m = discrete(icosphere(2));
        let k = ricci_curvature(&m, &vec![0.0; m.num_vertices()]).expect("valid mesh");
        let integrated: Vec<f64> = m
            .gaussian_curvature()
            .iter()
            .zip(m.mixed_voronoi_areas().iter())
            .map(|(kv, av)| kv * av)
            .collect();
        assert!(largest_error(&k, &integrated) < 1e-13);
        assert!((sum(&k) - 4.0 * PI).abs() < 1e-11);
        // The regular icosahedron has the angle defect pi / 3 at every vertex.
        let ico = unit_icosahedron();
        let k = ricci_curvature(&ico, &[0.0; 12]).expect("valid mesh");
        assert!(k.iter().all(|x| (x - PI / 3.0).abs() < 1e-13));
    }

    #[test]
    fn gauss_bonnet_holds_for_every_conformal_factor() {
        // Closed sphere and torus, a disc (cap with one boundary loop) and an annulus.
        let cases = [
            (discrete(icosphere(2)), 4.0 * PI),
            (unit_tetrahedron(), 4.0 * PI),
            (discrete(torus(16, 8, 3.0, 1.0)), 0.0),
            (discrete(spherical_cap(5, 20, 1.2)), 2.0 * PI),
            (discrete(flat_grid(6, 5, 1.0)), 2.0 * PI),
            (discrete(open_cylinder(12, 4, 1.0, 2.0)), 0.0),
        ];
        for (k, (m, expected)) in cases.iter().enumerate() {
            for seed in 0..3 {
                let u = random_factors(m.num_vertices(), 0.3, 10 * k as u64 + seed + 1);
                let curvature = ricci_curvature(m, &u).expect("valid mesh");
                assert!(
                    (sum(&curvature) - expected).abs() < 1e-11,
                    "case {k}: total curvature {} vs {expected}",
                    sum(&curvature)
                );
            }
        }
    }

    #[test]
    fn boundary_vertices_use_the_straight_angle() {
        // A flat grid has zero curvature everywhere except at the corners of its
        // boundary: interior vertices 2 pi - 2 pi, edge vertices pi - pi.
        let m = discrete(flat_grid(4, 4, 1.0));
        let k = ricci_curvature(&m, &vec![0.0; m.num_vertices()]).expect("valid mesh");
        for j in 0..=4 {
            for i in 0..=4 {
                let v = j * 5 + i;
                let corner = (i == 0 || i == 4) && (j == 0 || j == 4);
                if corner {
                    assert!((k[v] - PI / 2.0).abs() < 1e-13, "corner {v}: {}", k[v]);
                } else {
                    assert!(k[v].abs() < 1e-13, "vertex {v}: {}", k[v]);
                }
            }
        }
    }

    #[test]
    fn uniform_shift_rescales_the_metric_without_changing_the_curvature() {
        let m = discrete(icosphere(1));
        let u = random_factors(m.num_vertices(), 0.2, 3);
        let shifted: Vec<f64> = u.iter().map(|x| x + 0.7).collect();
        let a = ricci_curvature(&m, &u).expect("valid mesh");
        let b = ricci_curvature(&m, &shifted).expect("valid mesh");
        assert!(largest_error(&a, &b) < 1e-12);
        let la = conformal_edge_lengths(&m, &u).expect("valid mesh");
        let lb = conformal_edge_lengths(&m, &shifted).expect("valid mesh");
        for (x, y) in la.iter().zip(lb.iter()) {
            assert!((y.2 / x.2 - 0.7_f64.exp()).abs() < 1e-12);
        }
    }

    #[test]
    fn edge_lengths_follow_the_conformal_factors() {
        let m = discrete(icosphere(1));
        let u = random_factors(m.num_vertices(), 0.4, 8);
        let lengths = conformal_edge_lengths(&m, &u).expect("valid mesh");
        assert_eq!(lengths.len(), 120);
        for &(a, b, l) in &lengths {
            assert!(a < b);
            let expected =
                super::super::dist3(m.vertices[a], m.vertices[b]) * ((u[a] + u[b]) / 2.0).exp();
            assert!((l - expected).abs() < 1e-13 * expected, "edge ({a}, {b})");
        }
        assert!(
            lengths
                .windows(2)
                .all(|w| (w[0].0, w[0].1) < (w[1].0, w[1].1))
        );
    }

    #[test]
    fn jacobian_is_the_cotangent_laplacian_and_matches_finite_differences() {
        let m = discrete(icosphere(1));
        let n = m.num_vertices();
        let u = random_factors(n, 0.2, 5);
        let surface = RicciSurface::new(&m).expect("valid mesh");
        let mut dense = vec![vec![0.0_f64; n]; n];
        for &(i, j, value) in &surface.jacobian_entries(&u).expect("finite factors") {
            dense[i][j] += value;
        }
        let h = 1e-6;
        for j in 0..n {
            let mut up = u.clone();
            let mut down = u.clone();
            up[j] += h;
            down[j] -= h;
            let kp = surface.curvature(&up).expect("finite factors");
            let km = surface.curvature(&down).expect("finite factors");
            for i in 0..n {
                let derivative = (kp[i] - km[i]) / (2.0 * h);
                assert!(
                    (dense[i][j] - derivative).abs() < 1e-6,
                    "dK[{i}]/du[{j}]: analytic {} vs finite differences {derivative}",
                    dense[i][j]
                );
            }
        }
        for (i, row) in dense.iter().enumerate() {
            assert!(row.iter().sum::<f64>().abs() < 1e-12, "row {i}");
            for (j, value) in row.iter().enumerate() {
                assert!(
                    (value - dense[j][i]).abs() < 1e-14,
                    "asymmetric at ({i}, {j})"
                );
            }
        }
    }

    #[test]
    fn a_step_changes_the_curvature_and_moves_it_towards_the_target() {
        let m = perturbed_icosahedron(4, 0.2);
        let target = vec![PI / 3.0; 12];
        let mut u = vec![0.0; 12];
        let before = ricci_curvature(&m, &u).expect("valid mesh");
        let error_before = largest_error(&before, &target);
        assert!(
            error_before > 0.1,
            "the perturbation is too small: {error_before}"
        );
        let reported = ricci_flow_step(&m, &mut u, &target, 0.1).expect("valid input");
        assert!((reported - error_before).abs() < 1e-14);
        assert!(
            u.iter().any(|x| x.abs() > 1e-3),
            "the conformal factors did not move"
        );
        let after = ricci_curvature(&m, &u).expect("valid mesh");
        // The old code never changed the curvature.
        assert!(largest_error(&before, &after) > 1e-3);
        let norm = |k: &[f64]| -> f64 {
            k.iter()
                .zip(target.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        assert!(norm(&after) < norm(&before));
    }

    #[test]
    fn flow_lowers_the_factor_where_the_curvature_is_too_high() {
        let m = perturbed_icosahedron(9, 0.25);
        let target = vec![PI / 3.0; 12];
        let k = ricci_curvature(&m, &[0.0; 12]).expect("valid mesh");
        let mut u = vec![0.0; 12];
        ricci_flow_step(&m, &mut u, &target, 0.05).expect("valid input");
        for v in 0..12 {
            // du = step * (Kbar - K): opposite in sign to the excess curvature.
            assert!((u[v] - 0.05 * (target[v] - k[v])).abs() < 1e-14);
            assert!(u[v] * (k[v] - target[v]) <= 0.0);
        }
    }

    #[test]
    fn explicit_flow_converges_to_the_uniform_curvature() {
        let m = perturbed_icosahedron(4, 0.2);
        let target = vec![PI / 3.0; 12];
        let euclid = |u: &[f64]| -> f64 {
            let k = ricci_curvature(&m, u).expect("valid mesh");
            k.iter()
                .zip(target.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        let mut u = vec![0.0; 12];
        let initial = euclid(&u);
        let mut last_error = f64::INFINITY;
        for _ in 0..400 {
            last_error = ricci_flow_step(&m, &mut u, &target, 0.1).expect("valid input");
        }
        assert!(last_error < 1e-8, "explicit flow stalled at {last_error}");
        assert!(euclid(&u) < 1e-7 * initial);
    }

    #[test]
    fn newton_flow_reaches_the_uniform_curvature_on_a_perturbed_icosahedron() {
        let m = perturbed_icosahedron(4, 0.2);
        let target = vec![4.0 * PI / 12.0; 12];
        let mut u = vec![0.0; 12];
        let initial = largest_error(&ricci_curvature(&m, &u).expect("valid mesh"), &target);
        assert!(initial > 0.1);
        let summary = ricci_flow(&m, &mut u, &target, 1e-12, 30).expect("converges");
        assert!(summary.iterations <= 10, "{summary:?}");
        assert!(summary.residual <= 1e-12);
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(largest_error(&k, &target) < 1e-10);
        let spread = u.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - u.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            spread > 1e-2,
            "the conformal factors should differ: {spread}"
        );
        // Every deformed triangle is a genuine triangle.
        let surface = RicciSurface::new(&m).expect("valid mesh");
        let scales = surface.edge_scales(&u).expect("finite factors");
        for f in 0..m.num_faces() {
            assert!(triangle_corners(surface.sides(f, &scales)).four_area > 0.0);
        }
    }

    #[test]
    fn newton_flow_recovers_a_known_conformal_factor() {
        let m = discrete(icosphere(1));
        let n = m.num_vertices();
        let exact = random_factors(n, 0.15, 17);
        let target = ricci_curvature(&m, &exact).expect("valid mesh");
        let mut u = vec![0.0; n];
        ricci_flow(&m, &mut u, &target, 1e-12, 30).expect("converges");
        // The factors are determined up to a constant (the rescaling of the surface).
        let shift = u[0] - exact[0];
        for v in 0..n {
            assert!((u[v] - shift - exact[v]).abs() < 1e-8, "vertex {v}");
        }
    }

    #[test]
    fn newton_flow_flattens_a_spherical_cap_to_a_round_disc() {
        let m = discrete(spherical_cap(6, 24, 1.2));
        let n = m.num_vertices();
        let table = EdgeTable::build(n, &m.triangles).expect("manifold");
        let boundary = table.boundary_vertices(n);
        let boundary_count = boundary.iter().filter(|&&b| b).count();
        assert_eq!(boundary_count, 24);
        // Flat inside; the whole turning of 2 pi is spread evenly over the boundary.
        let target: Vec<f64> = boundary
            .iter()
            .map(|&b| {
                if b {
                    2.0 * PI / boundary_count as f64
                } else {
                    0.0
                }
            })
            .collect();
        let mut u = vec![0.0; n];
        let initial = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(largest_error(&initial, &target) > 0.05);
        let summary = ricci_flow(&m, &mut u, &target, 1e-11, 40).expect("converges");
        assert!(summary.residual <= 1e-11);
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        for v in 0..n {
            assert!((k[v] - target[v]).abs() < 1e-10, "vertex {v}: {}", k[v]);
        }
        // The rotational symmetry is kept: all vertices of a ring get the same factor.
        for r in 1..=6 {
            let ring: Vec<f64> = (0..24).map(|i| u[1 + (r - 1) * 24 + i]).collect();
            let spread = ring.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - ring.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(spread < 1e-7, "ring {r} spread {spread}");
        }
    }

    #[test]
    fn newton_flow_flattens_a_torus() {
        let m = discrete(torus(16, 8, 3.0, 1.0));
        let n = m.num_vertices();
        let mut u = vec![0.0; n];
        let target = vec![0.0; n];
        let initial = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(initial.iter().any(|k| k.abs() > 0.1));
        let summary = ricci_flow(&m, &mut u, &target, 1e-11, 40).expect("converges");
        assert!(summary.iterations <= 12, "{summary:?}");
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(k.iter().all(|x| x.abs() < 1e-10));
    }

    #[test]
    fn target_must_satisfy_gauss_bonnet_on_every_component() {
        let m = discrete(icosphere(1));
        let n = m.num_vertices();
        let mut u = vec![0.0; n];
        let wrong = vec![PI / 3.0; n];
        assert!(matches!(
            ricci_flow(&m, &mut u, &wrong, 1e-8, 10),
            Err(Error::OutOfRange { .. })
        ));
        assert!(matches!(
            ricci_flow_step(&m, &mut u, &wrong, 0.1),
            Err(Error::OutOfRange { .. })
        ));
        assert!(
            u.iter().all(|x| *x == 0.0),
            "a rejected step must not change u"
        );
        // Two spheres: the total is right, but one sphere has 4 pi + 1 and the other 4 pi - 1.
        let (mut vertices, mut faces) = icosphere(0);
        let (more_vertices, more_faces) = icosphere(0);
        vertices.extend(more_vertices.iter().map(|p| [p[0] + 5.0, p[1], p[2]]));
        faces.extend(more_faces.iter().map(|t| [t[0] + 12, t[1] + 12, t[2] + 12]));
        let two = DiscreteMesh::new(vertices, faces);
        let mut target = vec![PI / 3.0; 24];
        target[0] += 1.0;
        target[12] -= 1.0;
        assert!(matches!(
            ricci_flow(&two, &mut [0.0; 24], &target, 1e-8, 10),
            Err(Error::OutOfRange { .. })
        ));
        // Without the imbalance each sphere is solved on its own.
        let mut u = vec![0.0; 24];
        let balanced = vec![PI / 3.0; 24];
        ricci_flow(&two, &mut u, &balanced, 1e-12, 10).expect("converges");
    }

    #[test]
    fn invalid_input_is_reported() {
        let m = discrete(icosphere(0));
        let target = vec![PI / 3.0; 12];
        let mut u = vec![0.0; 12];
        assert!(matches!(
            ricci_flow_step(&m, &mut u[..5], &target, 0.1),
            Err(Error::LengthMismatch { .. })
        ));
        assert!(matches!(
            ricci_flow_step(&m, &mut u, &target[..5], 0.1),
            Err(Error::LengthMismatch { .. })
        ));
        for bad_step in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                ricci_flow_step(&m, &mut u, &target, bad_step),
                Err(Error::OutOfRange { .. })
            ));
        }
        for bad_tol in [0.0, -1e-3, f64::NAN] {
            assert!(matches!(
                ricci_flow(&m, &mut u, &target, bad_tol, 5),
                Err(Error::OutOfRange { .. })
            ));
        }
        let mut nan_u = u.clone();
        nan_u[3] = f64::NAN;
        assert!(matches!(
            ricci_flow_step(&m, &mut nan_u, &target, 0.1),
            Err(Error::General(_))
        ));
        let mut nan_target = target.clone();
        nan_target[2] = f64::NAN;
        assert!(matches!(
            ricci_flow_step(&m, &mut u, &nan_target, 0.1),
            Err(Error::General(_))
        ));
        assert!(matches!(
            ricci_curvature(&m, &[1000.0; 12]),
            Err(Error::General(_))
        ));
        // A triangle that references a missing vertex, a repeated vertex (zero side),
        // a non-finite vertex and a non-manifold edge.
        let mut bad = m.clone();
        bad.triangles[1][2] = 40;
        assert!(matches!(
            ricci_curvature(&bad, &u),
            Err(Error::InvalidMesh { .. })
        ));
        let mut bad = m.clone();
        bad.vertices[2][0] = f64::NAN;
        assert!(matches!(
            ricci_curvature(&bad, &u),
            Err(Error::InvalidMesh { .. })
        ));
        let mut bad = m.clone();
        bad.triangles[0] = [0, 0, 1];
        assert!(matches!(
            ricci_curvature(&bad, &u),
            Err(Error::DegenerateGeometry { .. } | Error::InvalidMesh { .. })
        ));
        let mut bad = m.clone();
        bad.vertices.push([3.0, 3.0, 3.0]);
        let (a, b) = (bad.triangles[0][0], bad.triangles[0][1]);
        bad.triangles.push([a, b, 12]);
        assert!(matches!(
            ricci_curvature(&bad, &[0.0; 13]),
            Err(Error::InvalidMesh { .. })
        ));
    }

    #[test]
    fn meshes_without_triangles_have_nothing_to_flow() {
        let none = DiscreteMesh::new(Vec::new(), Vec::new());
        let summary = ricci_flow(&none, &mut [], &[], 1e-9, 5).expect("nothing to do");
        assert_eq!(summary.iterations, 0);
        assert_eq!(
            ricci_flow_step(&none, &mut [], &[], 0.1).expect("nothing to do"),
            0.0
        );
        assert!(
            ricci_curvature(&none, &[])
                .expect("nothing to do")
                .is_empty()
        );
        assert!(
            conformal_edge_lengths(&none, &[])
                .expect("nothing to do")
                .is_empty()
        );
        let points = DiscreteMesh::new(vec![[0.0; 3], [1.0, 0.0, 0.0]], Vec::new());
        assert_eq!(
            ricci_curvature(&points, &[0.3, -0.2]).expect("valid mesh"),
            vec![0.0, 0.0]
        );
    }

    #[test]
    fn too_few_iterations_are_reported_as_a_convergence_failure() {
        let m = perturbed_icosahedron(4, 0.2);
        let target = vec![PI / 3.0; 12];
        let mut u = vec![0.0; 12];
        match ricci_flow(&m, &mut u, &target, 1e-12, 1) {
            Err(Error::ConvergenceFailure {
                iterations,
                residual,
                ..
            }) => {
                assert_eq!(iterations, 1);
                assert!(residual > 1e-12);
            }
            other => panic!("expected a convergence failure, got {other:?}"),
        }
        // The first iterate was kept and is closer to the target than the start.
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(largest_error(&k, &target) < 0.3);
        // Zero iterations succeed when the input already satisfies the tolerance.
        let ico = unit_icosahedron();
        let summary =
            ricci_flow(&ico, &mut [0.0; 12], &target, 1e-9, 0).expect("already converged");
        assert_eq!(summary.iterations, 0);
    }

    #[test]
    fn vertices_without_triangles_are_left_alone() {
        let mut m = perturbed_icosahedron(4, 0.2);
        m.vertices.push([9.0, 9.0, 9.0]);
        let mut target = vec![PI / 3.0; 12];
        target.push(123.0);
        let mut u = vec![0.0; 13];
        u[12] = 0.25;
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        assert_eq!(k[12], 0.0);
        let summary = ricci_flow(&m, &mut u, &target, 1e-12, 20).expect("converges");
        assert!(summary.iterations >= 1, "the mesh is not uniform yet");
        assert_eq!(u[12], 0.25);
        assert!(u[..12].iter().any(|x| x.abs() > 1e-3));
        let mut u = vec![0.0; 13];
        ricci_flow_step(&m, &mut u, &target, 0.1).expect("valid input");
        assert_eq!(u[12], 0.0);
    }

    #[test]
    fn newton_flow_is_fast_on_a_larger_mesh() {
        let m = discrete(icosphere(4));
        let n = m.num_vertices();
        let exact = random_factors(n, 0.1, 31);
        let target = ricci_curvature(&m, &exact).expect("valid mesh");
        let mut u = vec![0.0; n];
        let start = std::time::Instant::now();
        let summary = ricci_flow(&m, &mut u, &target, 1e-10, 30).expect("converges");
        let elapsed = start.elapsed().as_secs_f64();
        assert!(elapsed < 60.0, "{n} vertices took {elapsed:.2} s");
        assert!(summary.iterations <= 10, "{summary:?}");
        let k = ricci_curvature(&m, &u).expect("valid mesh");
        assert!(largest_error(&k, &target) < 1e-9);
    }
}
