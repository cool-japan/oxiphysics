// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Detailed cloth simulation using Position-Based Dynamics (PBD).
//!
//! Provides a standalone cloth mesh with stretch, shear, and bend constraints
//! using plain `[f64; 3]` arrays (no external math dependencies).

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// ClothVertex
// ---------------------------------------------------------------------------

/// A single vertex (particle) in the cloth mesh.
#[derive(Debug, Clone)]
pub struct ClothVertex {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass of the vertex (kg).
    pub mass: f64,
    /// When `true` the vertex is kinematically fixed.
    pub pinned: bool,
    /// Texture coordinate.
    pub uv: [f64; 2],
}

impl ClothVertex {
    /// Create a new cloth vertex.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            pinned: false,
            uv: [0.0; 2],
        }
    }
}

// ---------------------------------------------------------------------------
// ClothConstraintType
// ---------------------------------------------------------------------------

/// Semantic type of a cloth constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClothConstraintType {
    /// Structural constraint that resists stretching along edges.
    Stretch,
    /// Diagonal constraint that resists shear deformation.
    Shear,
    /// Bend constraint across a shared edge connecting two triangles.
    Bend,
}

// ---------------------------------------------------------------------------
// ClothConstraint
// ---------------------------------------------------------------------------

/// A distance-type constraint between two cloth vertices.
#[derive(Debug, Clone)]
pub struct ClothConstraint {
    /// Index of the first vertex.
    pub i: usize,
    /// Index of the second vertex.
    pub j: usize,
    /// Rest (natural) length of the constraint (m).
    pub rest_length: f64,
    /// Stiffness coefficient in \[0, 1\].  1 = fully rigid.
    pub stiffness: f64,
    /// Semantic type of this constraint.
    pub constraint_type: ClothConstraintType,
}

// ---------------------------------------------------------------------------
// ClothMesh
// ---------------------------------------------------------------------------

/// A rectangular cloth mesh driven by PBD.
#[derive(Debug, Clone)]
pub struct ClothMesh {
    /// All vertices of the cloth.
    pub vertices: Vec<ClothVertex>,
    /// All constraints (stretch, shear, bend).
    pub constraints: Vec<ClothConstraint>,
    /// Number of vertices along the width direction.
    pub width: usize,
    /// Number of vertices along the height direction.
    pub height: usize,
}

// ---------------------------------------------------------------------------
// create_cloth_grid
// ---------------------------------------------------------------------------

/// Create a flat rectangular cloth grid lying in the XZ plane.
///
/// * `w`, `h` – vertex count in each dimension (>= 2).
/// * `spacing` – uniform spacing between adjacent vertices (m).
/// * `mass` – mass per vertex (kg).
/// * `stiffness` – structural stretch constraint stiffness (0–1).
pub fn create_cloth_grid(w: usize, h: usize, spacing: f64, mass: f64, stiffness: f64) -> ClothMesh {
    let mut vertices = Vec::with_capacity(w * h);
    for row in 0..h {
        for col in 0..w {
            let mut v = ClothVertex::new([col as f64 * spacing, 0.0, row as f64 * spacing], mass);
            v.uv = [
                col as f64 / (w - 1).max(1) as f64,
                row as f64 / (h - 1).max(1) as f64,
            ];
            vertices.push(v);
        }
    }

    let mut constraints = Vec::new();

    // Horizontal stretch
    for row in 0..h {
        for col in 0..w - 1 {
            let i = row * w + col;
            let j = row * w + col + 1;
            let rest = spacing;
            constraints.push(ClothConstraint {
                i,
                j,
                rest_length: rest,
                stiffness,
                constraint_type: ClothConstraintType::Stretch,
            });
        }
    }
    // Vertical stretch
    for row in 0..h - 1 {
        for col in 0..w {
            let i = row * w + col;
            let j = (row + 1) * w + col;
            let rest = spacing;
            constraints.push(ClothConstraint {
                i,
                j,
                rest_length: rest,
                stiffness,
                constraint_type: ClothConstraintType::Stretch,
            });
        }
    }

    ClothMesh {
        vertices,
        constraints,
        width: w,
        height: h,
    }
}

// ---------------------------------------------------------------------------
// add_shear_constraints
// ---------------------------------------------------------------------------

/// Add diagonal shear constraints to an existing cloth mesh.
pub fn add_shear_constraints(cloth: &mut ClothMesh, stiffness: f64) {
    let w = cloth.width;
    let h = cloth.height;
    let diag = (2.0f64).sqrt() * {
        // Estimate spacing from first horizontal constraint if available, else 1.
        cloth
            .constraints
            .iter()
            .find(|c| c.constraint_type == ClothConstraintType::Stretch)
            .map(|c| c.rest_length)
            .unwrap_or(1.0)
    };
    for row in 0..h - 1 {
        for col in 0..w - 1 {
            let i = row * w + col;
            let j = (row + 1) * w + col + 1;
            cloth.constraints.push(ClothConstraint {
                i,
                j,
                rest_length: diag,
                stiffness,
                constraint_type: ClothConstraintType::Shear,
            });
            let i2 = row * w + col + 1;
            let j2 = (row + 1) * w + col;
            cloth.constraints.push(ClothConstraint {
                i: i2,
                j: j2,
                rest_length: diag,
                stiffness,
                constraint_type: ClothConstraintType::Shear,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// add_bend_constraints
// ---------------------------------------------------------------------------

/// Add long-range bend (flexion) constraints to the cloth mesh.
///
/// Connects every vertex to its vertex two steps away in the horizontal and
/// vertical directions.
pub fn add_bend_constraints(cloth: &mut ClothMesh, stiffness: f64) {
    let w = cloth.width;
    let h = cloth.height;
    // Horizontal bend (skip one vertex)
    for row in 0..h {
        for col in 0..w.saturating_sub(2) {
            let i = row * w + col;
            let j = row * w + col + 2;
            let rest = len3(sub3(cloth.vertices[j].position, cloth.vertices[i].position));
            cloth.constraints.push(ClothConstraint {
                i,
                j,
                rest_length: rest,
                stiffness,
                constraint_type: ClothConstraintType::Bend,
            });
        }
    }
    // Vertical bend
    for row in 0..h.saturating_sub(2) {
        for col in 0..w {
            let i = row * w + col;
            let j = (row + 2) * w + col;
            let rest = len3(sub3(cloth.vertices[j].position, cloth.vertices[i].position));
            cloth.constraints.push(ClothConstraint {
                i,
                j,
                rest_length: rest,
                stiffness,
                constraint_type: ClothConstraintType::Bend,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// apply_wind_force
// ---------------------------------------------------------------------------

/// Apply a uniform wind force to all un-pinned cloth vertices for one time step.
///
/// Wind is applied proportional to the projection of the wind velocity onto
/// the local surface normal, scaled by `dt` and vertex mass.
pub fn apply_wind_force(cloth: &mut ClothMesh, wind: [f64; 3], dt: f64) {
    let w = cloth.width;
    let h = cloth.height;
    // Accumulate per-vertex forces from triangle normals.
    let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; cloth.vertices.len()];
    for row in 0..h - 1 {
        for col in 0..w - 1 {
            let i00 = row * w + col;
            let i10 = row * w + col + 1;
            let i01 = (row + 1) * w + col;
            // Triangle (i00, i10, i01)
            let p0 = cloth.vertices[i00].position;
            let p1 = cloth.vertices[i10].position;
            let p2 = cloth.vertices[i01].position;
            let e1 = sub3(p1, p0);
            let e2 = sub3(p2, p0);
            let normal = cross3(e1, e2);
            let area = len3(normal) * 0.5;
            let n_dot_w = dot3(normal, wind);
            if area < 1e-30 || n_dot_w.abs() < 1e-30 {
                continue;
            }
            let f_tri = scale3(wind, n_dot_w * area / 3.0);
            forces[i00] = add3(forces[i00], f_tri);
            forces[i10] = add3(forces[i10], f_tri);
            forces[i01] = add3(forces[i01], f_tri);
        }
    }
    for (idx, v) in cloth.vertices.iter_mut().enumerate() {
        if !v.pinned && v.mass > 1e-30 {
            let impulse = scale3(forces[idx], dt / v.mass);
            v.velocity = add3(v.velocity, impulse);
        }
    }
}

// ---------------------------------------------------------------------------
// resolve_stretch_constraint
// ---------------------------------------------------------------------------

/// Resolve a single stretch constraint in the cloth mesh using PBD projection.
///
/// Runs `n_iters` Gauss–Seidel iterations on the constraint at index `c_idx`.
pub fn resolve_stretch_constraint(cloth: &mut ClothMesh, c_idx: usize, n_iters: usize) {
    for _ in 0..n_iters {
        let c = &cloth.constraints[c_idx];
        let i = c.i;
        let j = c.j;
        let rest = c.rest_length;
        let k = c.stiffness;

        let pi = cloth.vertices[i].position;
        let pj = cloth.vertices[j].position;
        let delta = sub3(pj, pi);
        let dist = len3(delta);
        if dist < 1e-15 {
            continue;
        }
        let correction = (dist - rest) / dist;
        let wi = if cloth.vertices[i].pinned {
            0.0
        } else {
            1.0 / cloth.vertices[i].mass.max(1e-30)
        };
        let wj = if cloth.vertices[j].pinned {
            0.0
        } else {
            1.0 / cloth.vertices[j].mass.max(1e-30)
        };
        let sum_w = wi + wj;
        if sum_w < 1e-30 {
            continue;
        }
        let corr = scale3(delta, k * correction / sum_w);
        if !cloth.vertices[i].pinned {
            cloth.vertices[i].position = add3(pi, scale3(corr, wi));
        }
        if !cloth.vertices[j].pinned {
            cloth.vertices[j].position = sub3(pj, scale3(corr, wj));
        }
    }
}

// ---------------------------------------------------------------------------
// cloth_step_pbd
// ---------------------------------------------------------------------------

/// Advance the cloth by one PBD time step.
///
/// * Applies gravity to velocities.
/// * Integrates positions.
/// * Resolves all constraints for `n_iters` iterations.
/// * Updates velocities from position change.
pub fn cloth_step_pbd(cloth: &mut ClothMesh, gravity: [f64; 3], dt: f64, n_iters: usize) {
    // Store old positions and apply external forces.
    let old_positions: Vec<[f64; 3]> = cloth.vertices.iter().map(|v| v.position).collect();

    for v in cloth.vertices.iter_mut() {
        if v.pinned {
            continue;
        }
        // Velocity update with gravity.
        v.velocity = add3(v.velocity, scale3(gravity, dt));
        // Position prediction.
        v.position = add3(v.position, scale3(v.velocity, dt));
    }

    // Constraint resolution.
    for _ in 0..n_iters {
        for c_idx in 0..cloth.constraints.len() {
            let c = cloth.constraints[c_idx].clone();
            let pi = cloth.vertices[c.i].position;
            let pj = cloth.vertices[c.j].position;
            let delta = sub3(pj, pi);
            let dist = len3(delta);
            if dist < 1e-15 {
                continue;
            }
            let correction = (dist - c.rest_length) / dist;
            let wi = if cloth.vertices[c.i].pinned {
                0.0
            } else {
                1.0 / cloth.vertices[c.i].mass.max(1e-30)
            };
            let wj = if cloth.vertices[c.j].pinned {
                0.0
            } else {
                1.0 / cloth.vertices[c.j].mass.max(1e-30)
            };
            let sum_w = wi + wj;
            if sum_w < 1e-30 {
                continue;
            }
            let corr = scale3(delta, c.stiffness * correction / sum_w);
            if !cloth.vertices[c.i].pinned {
                cloth.vertices[c.i].position = add3(pi, scale3(corr, wi));
            }
            if !cloth.vertices[c.j].pinned {
                cloth.vertices[c.j].position = sub3(pj, scale3(corr, wj));
            }
        }
    }

    // Velocity update from position change (PBD: v = (x_new - x_old) / dt).
    let inv_dt = if dt > 1e-30 { 1.0 / dt } else { 0.0 };
    for (idx, v) in cloth.vertices.iter_mut().enumerate() {
        if !v.pinned {
            v.velocity = scale3(sub3(v.position, old_positions[idx]), inv_dt);
        }
    }
}

// ---------------------------------------------------------------------------
// cloth_sphere_collision
// ---------------------------------------------------------------------------

/// Resolve cloth–sphere collision by pushing penetrating vertices outside the sphere.
pub fn cloth_sphere_collision(cloth: &mut ClothMesh, sphere_center: [f64; 3], sphere_radius: f64) {
    for v in cloth.vertices.iter_mut() {
        if v.pinned {
            continue;
        }
        let diff = sub3(v.position, sphere_center);
        let dist = len3(diff);
        if dist < sphere_radius && dist > 1e-15 {
            // Push vertex to surface.
            let normal = scale3(diff, 1.0 / dist);
            v.position = add3(sphere_center, scale3(normal, sphere_radius));
            // Remove inward velocity component.
            let vn = dot3(v.velocity, normal);
            if vn < 0.0 {
                v.velocity = sub3(v.velocity, scale3(normal, vn));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. create_cloth_grid: correct vertex count.
    #[test]
    fn test_create_cloth_grid_vertex_count() {
        let cloth = create_cloth_grid(5, 4, 0.1, 1.0, 1.0);
        assert_eq!(cloth.vertices.len(), 20);
    }

    // 2. create_cloth_grid: width and height stored correctly.
    #[test]
    fn test_create_cloth_grid_dimensions() {
        let cloth = create_cloth_grid(6, 3, 0.1, 1.0, 1.0);
        assert_eq!(cloth.width, 6);
        assert_eq!(cloth.height, 3);
    }

    // 3. create_cloth_grid: no pinned vertices by default.
    #[test]
    fn test_create_cloth_grid_no_pins() {
        let cloth = create_cloth_grid(3, 3, 0.1, 1.0, 1.0);
        assert!(cloth.vertices.iter().all(|v| !v.pinned));
    }

    // 4. create_cloth_grid: first vertex at origin.
    #[test]
    fn test_create_cloth_grid_first_vertex() {
        let cloth = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        let p = cloth.vertices[0].position;
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }

    // 5. create_cloth_grid: stretch constraint count = h*(w-1) + w*(h-1).
    #[test]
    fn test_create_cloth_grid_constraint_count() {
        let w = 4usize;
        let h = 3usize;
        let cloth = create_cloth_grid(w, h, 0.1, 1.0, 1.0);
        let expected = h * (w - 1) + w * (h - 1);
        assert_eq!(cloth.constraints.len(), expected);
    }

    // 6. create_cloth_grid: all stretch constraints have positive rest_length.
    #[test]
    fn test_create_cloth_grid_rest_lengths_positive() {
        let cloth = create_cloth_grid(4, 4, 0.5, 1.0, 1.0);
        assert!(cloth.constraints.iter().all(|c| c.rest_length > 0.0));
    }

    // 7. add_shear_constraints: adds 2*(w-1)*(h-1) constraints.
    #[test]
    fn test_add_shear_constraint_count() {
        let mut cloth = create_cloth_grid(4, 3, 0.1, 1.0, 1.0);
        let before = cloth.constraints.len();
        add_shear_constraints(&mut cloth, 0.8);
        let added = cloth.constraints.len() - before;
        let expected = 2 * (cloth.width - 1) * (cloth.height - 1);
        assert_eq!(added, expected);
    }

    // 8. add_shear_constraints: all added constraints are type Shear.
    #[test]
    fn test_add_shear_types() {
        let mut cloth = create_cloth_grid(3, 3, 0.1, 1.0, 1.0);
        let before = cloth.constraints.len();
        add_shear_constraints(&mut cloth, 0.5);
        for c in &cloth.constraints[before..] {
            assert_eq!(c.constraint_type, ClothConstraintType::Shear);
        }
    }

    // 9. add_bend_constraints: adds constraints of type Bend.
    #[test]
    fn test_add_bend_types() {
        let mut cloth = create_cloth_grid(4, 4, 0.1, 1.0, 1.0);
        let before = cloth.constraints.len();
        add_bend_constraints(&mut cloth, 0.3);
        for c in &cloth.constraints[before..] {
            assert_eq!(c.constraint_type, ClothConstraintType::Bend);
        }
    }

    // 10. add_bend_constraints: rest lengths are positive.
    #[test]
    fn test_add_bend_rest_lengths() {
        let mut cloth = create_cloth_grid(4, 4, 0.2, 1.0, 1.0);
        add_bend_constraints(&mut cloth, 0.3);
        let bend: Vec<_> = cloth
            .constraints
            .iter()
            .filter(|c| c.constraint_type == ClothConstraintType::Bend)
            .collect();
        assert!(!bend.is_empty());
        assert!(bend.iter().all(|c| c.rest_length > 0.0));
    }

    // 11. apply_wind_force: zero wind → no velocity change.
    #[test]
    fn test_wind_zero() {
        let mut cloth = create_cloth_grid(3, 3, 0.1, 1.0, 1.0);
        apply_wind_force(&mut cloth, [0.0, 0.0, 0.0], 0.01);
        assert!(cloth.vertices.iter().all(|v| v.velocity == [0.0, 0.0, 0.0]));
    }

    // 12. apply_wind_force: non-zero wind → at least one vertex gains velocity.
    #[test]
    fn test_wind_nonzero() {
        let mut cloth = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        apply_wind_force(&mut cloth, [0.0, 10.0, 0.0], 0.01);
        let any_moved = cloth.vertices.iter().any(|v| v.velocity[1].abs() > 1e-15);
        assert!(
            any_moved,
            "Wind should impart velocity to at least one vertex"
        );
    }

    // 13. apply_wind_force: pinned vertices do not move.
    #[test]
    fn test_wind_pinned_stays() {
        let mut cloth = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        cloth.vertices[0].pinned = true;
        apply_wind_force(&mut cloth, [0.0, 50.0, 0.0], 0.1);
        assert_eq!(cloth.vertices[0].velocity, [0.0, 0.0, 0.0]);
    }

    // 14. resolve_stretch_constraint: constraint enforced after many iterations.
    #[test]
    fn test_resolve_stretch_convergence() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        // Stretch first two vertices apart.
        cloth.vertices[1].position = [3.0, 0.0, 0.0];
        resolve_stretch_constraint(&mut cloth, 0, 200);
        let i = cloth.constraints[0].i;
        let j = cloth.constraints[0].j;
        let pi = cloth.vertices[i].position;
        let pj = cloth.vertices[j].position;
        let dist = len3(sub3(pj, pi));
        let rest = cloth.constraints[0].rest_length;
        assert!((dist - rest).abs() < 0.05, "dist={dist} should be ≈ {rest}");
    }

    // 15. cloth_step_pbd: gravity causes unpinned vertices to fall.
    #[test]
    fn test_pbd_step_gravity() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            cloth_step_pbd(&mut cloth, [0.0, -9.81, 0.0], dt, 3);
        }
        let y0 = cloth.vertices[0].position[1];
        assert!(y0 < -0.1, "Cloth vertex should fall under gravity, y={y0}");
    }

    // 16. cloth_step_pbd: pinned vertices stay in place.
    #[test]
    fn test_pbd_step_pinned_fixed() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        cloth.vertices[0].pinned = true;
        let orig = cloth.vertices[0].position;
        cloth_step_pbd(&mut cloth, [0.0, -9.81, 0.0], 1.0 / 60.0, 5);
        assert_eq!(cloth.vertices[0].position, orig);
    }

    // 17. cloth_step_pbd: velocity is updated from position change.
    #[test]
    fn test_pbd_step_velocity_updated() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        cloth_step_pbd(&mut cloth, [0.0, -9.81, 0.0], 1.0 / 60.0, 1);
        // At least one unpinned vertex should have non-zero velocity.
        assert!(cloth.vertices.iter().any(|v| v.velocity[1].abs() > 1e-10));
    }

    // 18. cloth_sphere_collision: vertices inside sphere are pushed out.
    #[test]
    fn test_sphere_collision_pushes_out() {
        let mut cloth = create_cloth_grid(2, 2, 0.1, 1.0, 1.0);
        // Place all vertices at the sphere centre.
        for v in cloth.vertices.iter_mut() {
            v.position = [0.0, 0.0, 0.0];
            v.position[0] += 1e-6; // tiny offset to avoid zero-length direction
        }
        cloth_sphere_collision(&mut cloth, [0.0, 0.0, 0.0], 1.0);
        for v in &cloth.vertices {
            let dist = len3(sub3(v.position, [0.0, 0.0, 0.0]));
            assert!(
                dist >= 1.0 - 1e-9,
                "Vertex should be on or outside sphere: dist={dist}"
            );
        }
    }

    // 19. cloth_sphere_collision: vertices outside sphere are untouched.
    #[test]
    fn test_sphere_collision_no_change_outside() {
        let mut cloth = create_cloth_grid(3, 3, 2.0, 1.0, 1.0);
        // All vertices at z > 0, x > 0, at least some far from sphere at (−10, 0, 0).
        let before: Vec<_> = cloth.vertices.iter().map(|v| v.position).collect();
        cloth_sphere_collision(&mut cloth, [-10.0, 0.0, 0.0], 0.1);
        for (b, v) in before.iter().zip(cloth.vertices.iter()) {
            assert_eq!(*b, v.position, "Vertex outside sphere should not move");
        }
    }

    // 20. ClothConstraintType: variants are distinct.
    #[test]
    fn test_constraint_type_distinct() {
        assert_ne!(ClothConstraintType::Stretch, ClothConstraintType::Shear);
        assert_ne!(ClothConstraintType::Stretch, ClothConstraintType::Bend);
        assert_ne!(ClothConstraintType::Shear, ClothConstraintType::Bend);
    }

    // 21. ClothVertex: default velocity is zero.
    #[test]
    fn test_cloth_vertex_default_velocity() {
        let v = ClothVertex::new([1.0, 2.0, 3.0], 0.5);
        assert_eq!(v.velocity, [0.0, 0.0, 0.0]);
    }

    // 22. ClothVertex: not pinned by default.
    #[test]
    fn test_cloth_vertex_not_pinned() {
        let v = ClothVertex::new([0.0, 0.0, 0.0], 1.0);
        assert!(!v.pinned);
    }

    // 23. create_cloth_grid: spacing is reflected in vertex positions.
    #[test]
    fn test_create_cloth_grid_spacing() {
        let cloth = create_cloth_grid(3, 2, 0.5, 1.0, 1.0);
        let p0 = cloth.vertices[0].position;
        let p1 = cloth.vertices[1].position;
        let dx = p1[0] - p0[0];
        assert!((dx - 0.5).abs() < 1e-12);
    }

    // 24. resolve_stretch_constraint: pinned-pinned pair → no change.
    #[test]
    fn test_resolve_stretch_both_pinned() {
        let mut cloth = create_cloth_grid(2, 1, 1.0, 1.0, 1.0);
        cloth.vertices[0].pinned = true;
        cloth.vertices[1].pinned = true;
        cloth.vertices[1].position = [5.0, 0.0, 0.0]; // stretch far
        let before = cloth.vertices[1].position;
        resolve_stretch_constraint(&mut cloth, 0, 50);
        assert_eq!(
            cloth.vertices[1].position, before,
            "Pinned vertices must not move"
        );
    }

    // 25. cloth_step_pbd: n_iters=0 still integrates velocity.
    #[test]
    fn test_pbd_zero_iters() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        cloth_step_pbd(&mut cloth, [0.0, -9.81, 0.0], 0.01, 0);
        // Vertices should have moved (gravity applied).
        let y = cloth.vertices[0].position[1];
        assert!(
            y < 0.0,
            "Gravity should move vertices even with 0 iters, y={y}"
        );
    }

    // 26. create_cloth_grid: UV coordinates assigned for corner vertices.
    #[test]
    fn test_create_cloth_grid_uv() {
        let cloth = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        // Bottom-left (row=0,col=0) → uv=(0,0).
        assert_eq!(cloth.vertices[0].uv, [0.0, 0.0]);
        // Bottom-right (row=0,col=2) → uv=(1,0).
        assert!((cloth.vertices[2].uv[0] - 1.0).abs() < 1e-9);
    }

    // 27. add_bend_constraints: 2x2 grid adds no bend constraints (too small).
    #[test]
    fn test_bend_2x2_no_bends() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 1.0);
        let before = cloth.constraints.len();
        add_bend_constraints(&mut cloth, 0.5);
        // 2x2: w-2=0 horizontal bends, h-2=0 vertical bends.
        assert_eq!(
            cloth.constraints.len(),
            before,
            "2x2 grid should add no bend constraints"
        );
    }

    // 28. cloth_sphere_collision: pinned vertices are not moved by collision.
    #[test]
    fn test_sphere_collision_pinned() {
        let mut cloth = create_cloth_grid(2, 2, 0.01, 1.0, 1.0);
        for v in cloth.vertices.iter_mut() {
            v.position = [1e-6, 0.0, 0.0];
            v.pinned = true;
        }
        let before: Vec<_> = cloth.vertices.iter().map(|v| v.position).collect();
        cloth_sphere_collision(&mut cloth, [0.0, 0.0, 0.0], 5.0);
        for (b, v) in before.iter().zip(cloth.vertices.iter()) {
            assert_eq!(*b, v.position);
        }
    }

    // 29. cloth_step_pbd: multiple steps accumulate displacement.
    #[test]
    fn test_pbd_accumulates_displacement() {
        let mut cloth = create_cloth_grid(2, 2, 1.0, 1.0, 0.0); // stiffness=0
        let dt = 0.01;
        for _ in 0..10 {
            cloth_step_pbd(&mut cloth, [0.0, -9.81, 0.0], dt, 1);
        }
        let y = cloth.vertices[0].position[1];
        let expected = -0.5 * 9.81 * (10.0 * dt) * (10.0 * dt);
        // Should be in the right ballpark (constraints may shift it slightly).
        assert!(
            y < -0.01,
            "Expected significant sag, got y={y} (expected ≈ {expected:.4})"
        );
    }

    // 30. apply_wind_force: velocity scales with dt.
    #[test]
    fn test_wind_velocity_scales_dt() {
        let mut c1 = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        let mut c2 = create_cloth_grid(3, 3, 1.0, 1.0, 1.0);
        apply_wind_force(&mut c1, [0.0, 10.0, 0.0], 0.01);
        apply_wind_force(&mut c2, [0.0, 10.0, 0.0], 0.02);
        let v1: f64 = c1.vertices.iter().map(|v| v.velocity[1]).sum();
        let v2: f64 = c2.vertices.iter().map(|v| v.velocity[1]).sum();
        // v2 should be roughly twice v1 (linear in dt).
        assert!((v2 - 2.0 * v1).abs() < v1.abs() * 0.01 + 1e-15);
    }
}
