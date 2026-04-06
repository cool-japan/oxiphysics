//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CollisionPlane, CollisionSphere};

/// Compute the cross product of two 3-vectors.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Compute the dot product of two 3-vectors.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Return the Euclidean length of a 3-vector.
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Normalize a 3-vector.  Returns the zero vector when the length is too small.
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}
/// Component-wise addition of two 3-vectors.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Component-wise subtraction `a - b`.
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scalar multiplication.
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Signed volume of a tetrahedron: `(1/6) * det([p1-p0, p2-p0, p3-p0])`.
///
/// The sign encodes orientation; use `.abs()` for physical volume.
pub fn tet_volume_signed(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    let e3 = sub3(p3, p0);
    dot3(e1, cross3(e2, e3)) / 6.0
}
/// Centroid of a tetrahedron (average of four vertices).
pub fn tet_centroid(positions: &[[f64; 3]; 4]) -> [f64; 3] {
    let sum = positions.iter().fold([0.0; 3], |acc, &p| add3(acc, p));
    scale3(sum, 0.25)
}
/// Compute the work done by internal pressure during a volume change.
///
/// W = p * (V_current - V_rest)
#[allow(dead_code)]
pub fn pressure_volume_work(pressure: f64, v_current: f64, v_rest: f64) -> f64 {
    pressure * (v_current - v_rest)
}
/// Compute the internal pressure using the ideal gas law:
///
/// p = n R T / V  →  p = p0 * V0 / V  (isothermal)
///
/// where `p0` and `v0` are the reference pressure and volume.
#[allow(dead_code)]
pub fn ideal_gas_pressure(p0: f64, v0: f64, v_current: f64) -> f64 {
    if v_current.abs() < 1e-15 {
        return f64::MAX * 0.5;
    }
    p0 * v0 / v_current
}
/// Compute the internal pressure using the adiabatic gas law:
///
/// p = p0 * (V0/V)^γ
///
/// where γ is the heat capacity ratio (1.4 for diatomic gas).
#[allow(dead_code)]
pub fn adiabatic_gas_pressure(p0: f64, v0: f64, v_current: f64, gamma: f64) -> f64 {
    if v_current.abs() < 1e-15 {
        return f64::MAX * 0.5;
    }
    p0 * (v0 / v_current).powf(gamma)
}
/// Compute the volume constraint gradient (Jacobian) for a tetrahedron.
///
/// Returns the four gradient vectors `dV/dp_i` for each vertex.
/// These are the cross-product normals of opposite faces, scaled by 1/6.
#[allow(dead_code)]
pub fn tet_volume_gradient(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
) -> [[f64; 3]; 4] {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    let e3 = sub3(p3, p0);
    let g0 = scale3(cross3(sub3(p2, p1), sub3(p3, p1)), 1.0 / 6.0);
    let g1 = scale3(cross3(e2, e3), 1.0 / 6.0);
    let g2 = scale3(cross3(e3, e1), 1.0 / 6.0);
    let g3 = scale3(cross3(e1, e2), 1.0 / 6.0);
    [g0, g1, g2, g3]
}
/// Compute the Lagrange multiplier lambda for the volume constraint.
///
/// lambda = -stiffness * C / (Σ w_i |g_i|²)
///
/// where C = V_current - V_rest.
#[allow(dead_code)]
pub fn volume_constraint_lambda(
    gradients: &[[f64; 3]; 4],
    inv_masses: &[f64; 4],
    constraint: f64,
    stiffness: f64,
) -> f64 {
    let denom: f64 = gradients
        .iter()
        .zip(inv_masses.iter())
        .map(|(g, &w)| w * dot3(*g, *g))
        .sum();
    if denom < 1e-15 {
        return 0.0;
    }
    -stiffness * constraint / denom
}
/// Compute the surface area of a tetrahedron (sum of 4 triangle face areas).
#[allow(dead_code)]
pub fn tet_surface_area(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let triangle_area = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| -> f64 {
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        len3(cross3(ab, ac)) * 0.5
    };
    triangle_area(p0, p1, p2)
        + triangle_area(p0, p1, p3)
        + triangle_area(p0, p2, p3)
        + triangle_area(p1, p2, p3)
}
/// Compute the inradius of a tetrahedron: r = 3V / A.
#[allow(dead_code)]
pub fn tet_inradius(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let vol = tet_volume_signed(p0, p1, p2, p3).abs();
    let area = tet_surface_area(p0, p1, p2, p3);
    if area < 1e-15 {
        return 0.0;
    }
    3.0 * vol / area
}
/// Compute the total signed volume of a closed triangulated surface using the
/// divergence theorem. Each triangle contributes `(v0 · (v1 × v2)) / 6`.
#[allow(dead_code)]
pub fn signed_surface_volume(triangles: &[[usize; 3]], positions: &[[f64; 3]]) -> f64 {
    let mut vol = 0.0_f64;
    for tri in triangles {
        let v0 = positions[tri[0]];
        let v1 = positions[tri[1]];
        let v2 = positions[tri[2]];
        let cross = cross3(sub3(v1, v0), sub3(v2, v0));
        vol += dot3(v0, cross);
    }
    vol / 6.0
}
/// Gradient of the signed surface volume with respect to vertex `idx`.
/// Used in pressurized constraint corrections.
#[allow(dead_code)]
pub fn signed_surface_volume_gradient(
    idx: usize,
    triangles: &[[usize; 3]],
    positions: &[[f64; 3]],
) -> [f64; 3] {
    let mut grad = [0.0_f64; 3];
    for tri in triangles {
        let (i0, i1, i2) = (tri[0], tri[1], tri[2]);
        if i0 == idx {
            let c = cross3(positions[i1], positions[i2]);
            grad[0] += c[0];
            grad[1] += c[1];
            grad[2] += c[2];
        }
        if i1 == idx {
            let c = cross3(positions[i2], positions[i0]);
            grad[0] += c[0];
            grad[1] += c[1];
            grad[2] += c[2];
        }
        if i2 == idx {
            let c = cross3(positions[i0], positions[i1]);
            grad[0] += c[0];
            grad[1] += c[1];
            grad[2] += c[2];
        }
    }
    [grad[0] / 6.0, grad[1] / 6.0, grad[2] / 6.0]
}
/// Apply a pressurized volume constraint to a closed triangulated surface.
///
/// `rest_volume`  – target volume (e.g. initial inflated volume)
/// `pressure`     – over-pressure multiplier (1.0 = maintain rest volume exactly)
/// `stiffness`    – constraint stiffness in \[0, 1\]
/// `inv_masses`   – inverse masses of each vertex
/// `positions`    – mutable vertex positions
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn solve_pressurized_volume_constraint(
    triangles: &[[usize; 3]],
    rest_volume: f64,
    pressure: f64,
    stiffness: f64,
    inv_masses: &[f64],
    positions: &mut [[f64; 3]],
) {
    let target = rest_volume * pressure;
    let current = signed_surface_volume(triangles, positions);
    let c = current - target;
    if c.abs() < 1e-15 {
        return;
    }
    let n = positions.len();
    let grads: Vec<[f64; 3]> = (0..n)
        .map(|i| signed_surface_volume_gradient(i, triangles, positions))
        .collect();
    let denom: f64 = (0..n)
        .map(|i| {
            let g = grads[i];
            inv_masses[i] * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2])
        })
        .sum();
    if denom.abs() < 1e-30 {
        return;
    }
    let lambda = -stiffness * c / denom;
    for i in 0..n {
        let g = grads[i];
        positions[i][0] += lambda * inv_masses[i] * g[0];
        positions[i][1] += lambda * inv_masses[i] * g[1];
        positions[i][2] += lambda * inv_masses[i] * g[2];
    }
    let _ = grads;
}
/// Energy of the volume preservation penalty:
/// `E = 0.5 * k * (V - V0)^2`
#[allow(dead_code)]
pub fn volume_preservation_energy(
    tetrahedra: &[[usize; 4]],
    rest_volumes: &[f64],
    positions: &[[f64; 3]],
    stiffness: f64,
) -> f64 {
    let mut energy = 0.0_f64;
    for (t, tet) in tetrahedra.iter().enumerate() {
        let tet_pos = [
            positions[tet[0]],
            positions[tet[1]],
            positions[tet[2]],
            positions[tet[3]],
        ];
        let v = tet_volume_signed(tet_pos[0], tet_pos[1], tet_pos[2], tet_pos[3]).abs();
        let diff = v - rest_volumes[t];
        energy += 0.5 * stiffness * diff * diff;
    }
    energy
}
/// Gradient of the volume preservation penalty energy with respect to each vertex.
/// Returns a vector of per-vertex force vectors (negative gradient).
#[allow(dead_code)]
pub fn volume_preservation_gradient(
    tetrahedra: &[[usize; 4]],
    rest_volumes: &[f64],
    positions: &[[f64; 3]],
    stiffness: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    for (t, tet) in tetrahedra.iter().enumerate() {
        let p = [
            positions[tet[0]],
            positions[tet[1]],
            positions[tet[2]],
            positions[tet[3]],
        ];
        let v = tet_volume_signed(p[0], p[1], p[2], p[3]).abs();
        let sign = if tet_volume_signed(p[0], p[1], p[2], p[3]) >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let diff = v - rest_volumes[t];
        let coeff = stiffness * diff * sign;
        let dv_dp0 = scale3(cross3(sub3(p[1], p[3]), sub3(p[2], p[3])), 1.0 / 6.0);
        let dv_dp1 = scale3(cross3(sub3(p[2], p[3]), sub3(p[0], p[3])), 1.0 / 6.0);
        let dv_dp2 = scale3(cross3(sub3(p[0], p[3]), sub3(p[1], p[3])), 1.0 / 6.0);
        let dv_dp3 = scale3(cross3(sub3(p[1], p[0]), sub3(p[2], p[0])), 1.0 / 6.0);
        for k in 0..3 {
            grad[tet[0]][k] += coeff * dv_dp0[k];
            grad[tet[1]][k] += coeff * dv_dp1[k];
            grad[tet[2]][k] += coeff * dv_dp2[k];
            grad[tet[3]][k] += coeff * dv_dp3[k];
        }
    }
    grad
}
/// Boyle's law spring force: gas inside a deformable volume.
/// Returns the pressure given current and rest volume, and rest pressure.
/// `P * V = P0 * V0`  →  `P = P0 * V0 / V`
#[allow(dead_code)]
pub fn boyle_pressure(rest_pressure: f64, rest_volume: f64, current_volume: f64) -> f64 {
    if current_volume.abs() < 1e-30 {
        return 0.0;
    }
    rest_pressure * rest_volume / current_volume
}
/// Adiabatic gas law spring:
/// `P * V^γ = const`  →  `P = P0 * (V0/V)^γ`
#[allow(dead_code)]
pub fn adiabatic_pressure(
    rest_pressure: f64,
    rest_volume: f64,
    current_volume: f64,
    gamma: f64,
) -> f64 {
    if current_volume.abs() < 1e-30 {
        return 0.0;
    }
    rest_pressure * (rest_volume / current_volume).powf(gamma)
}
/// Apply Boyle's law gas pressure as a surface force to a closed triangle mesh.
/// Each triangle face receives an outward normal force proportional to pressure.
/// `positions` are updated in-place; `inv_masses` are per-vertex.
#[allow(dead_code)]
pub fn apply_gas_pressure_force(
    triangles: &[[usize; 3]],
    rest_pressure: f64,
    rest_volume: f64,
    inv_masses: &[f64],
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    dt: f64,
) {
    let current_volume = signed_surface_volume(triangles, positions).abs();
    let pressure = boyle_pressure(rest_pressure, rest_volume, current_volume);
    for tri in triangles {
        let v0 = positions[tri[0]];
        let v1 = positions[tri[1]];
        let v2 = positions[tri[2]];
        let e01 = sub3(v1, v0);
        let e02 = sub3(v2, v0);
        let area_normal = cross3(e01, e02);
        let area = len3(area_normal) * 0.5;
        let force_magnitude = pressure * area;
        let normal = if len3(area_normal) > 1e-30 {
            normalize3(area_normal)
        } else {
            [0.0, 0.0, 0.0]
        };
        for &vi in tri {
            let f = scale3(normal, force_magnitude / 3.0);
            velocities[vi][0] += f[0] * inv_masses[vi] * dt;
            velocities[vi][1] += f[1] * inv_masses[vi] * dt;
            velocities[vi][2] += f[2] * inv_masses[vi] * dt;
        }
    }
}
/// Resolve soft (penalty-based) collision between particles and a half-space.
/// Particles on the wrong side of the plane are pushed back with a force
/// proportional to penetration depth.
///
/// `stiffness` – penalty spring stiffness (force per unit penetration per unit mass)
#[allow(dead_code)]
pub fn resolve_plane_collision(
    plane: &CollisionPlane,
    stiffness: f64,
    restitution: f64,
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    inv_masses: &[f64],
    dt: f64,
) {
    let n = plane.normal;
    let p0 = plane.point;
    for i in 0..positions.len() {
        if inv_masses[i] == 0.0 {
            continue;
        }
        let diff = sub3(positions[i], p0);
        let dist = dot3(diff, n);
        if dist < 0.0 {
            let depth = -dist;
            positions[i][0] += n[0] * depth;
            positions[i][1] += n[1] * depth;
            positions[i][2] += n[2] * depth;
            let vn = dot3(velocities[i], n);
            if vn < 0.0 {
                velocities[i][0] -= (1.0 + restitution) * vn * n[0];
                velocities[i][1] -= (1.0 + restitution) * vn * n[1];
                velocities[i][2] -= (1.0 + restitution) * vn * n[2];
            }
            let f = scale3(n, stiffness * depth);
            velocities[i][0] += f[0] * inv_masses[i] * dt;
            velocities[i][1] += f[1] * inv_masses[i] * dt;
            velocities[i][2] += f[2] * inv_masses[i] * dt;
        }
    }
}
/// Resolve soft collision between particles and a solid sphere.
#[allow(dead_code)]
pub fn resolve_sphere_collision(
    sphere: &CollisionSphere,
    stiffness: f64,
    restitution: f64,
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    inv_masses: &[f64],
    dt: f64,
) {
    for i in 0..positions.len() {
        if inv_masses[i] == 0.0 {
            continue;
        }
        let diff = sub3(positions[i], sphere.center);
        let dist = len3(diff);
        if dist < sphere.radius && dist > 1e-30 {
            let depth = sphere.radius - dist;
            let normal = normalize3(diff);
            positions[i][0] += normal[0] * depth;
            positions[i][1] += normal[1] * depth;
            positions[i][2] += normal[2] * depth;
            let vn = dot3(velocities[i], normal);
            if vn < 0.0 {
                velocities[i][0] -= (1.0 + restitution) * vn * normal[0];
                velocities[i][1] -= (1.0 + restitution) * vn * normal[1];
                velocities[i][2] -= (1.0 + restitution) * vn * normal[2];
            }
            let f = scale3(normal, stiffness * depth);
            velocities[i][0] += f[0] * inv_masses[i] * dt;
            velocities[i][1] += f[1] * inv_masses[i] * dt;
            velocities[i][2] += f[2] * inv_masses[i] * dt;
        }
    }
}
/// Resolve soft collision between all pairs of particles with a given radius.
#[allow(dead_code)]
pub fn resolve_particle_collisions(
    particle_radius: f64,
    stiffness: f64,
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    inv_masses: &[f64],
    dt: f64,
) {
    let n = positions.len();
    let diameter = 2.0 * particle_radius;
    for i in 0..n {
        for j in (i + 1)..n {
            if inv_masses[i] == 0.0 && inv_masses[j] == 0.0 {
                continue;
            }
            let diff = sub3(positions[i], positions[j]);
            let dist = len3(diff);
            if dist < diameter && dist > 1e-30 {
                let depth = diameter - dist;
                let normal = normalize3(diff);
                let w_total = inv_masses[i] + inv_masses[j];
                if w_total < 1e-30 {
                    continue;
                }
                let wi = inv_masses[i] / w_total;
                let wj = inv_masses[j] / w_total;
                for k in 0..3 {
                    positions[i][k] += normal[k] * depth * wi;
                    positions[j][k] -= normal[k] * depth * wj;
                }
                let f = scale3(normal, stiffness * depth);
                for k in 0..3 {
                    velocities[i][k] += f[k] * inv_masses[i] * dt;
                    velocities[j][k] -= f[k] * inv_masses[j] * dt;
                }
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::volume::TetrahedralSoftBody;
    use crate::volume::VolumeParticle;
    fn regular_tet() -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
        let p0 = [1.0, 1.0, 1.0];
        let p1 = [1.0, -1.0, -1.0];
        let p2 = [-1.0, 1.0, -1.0];
        let p3 = [-1.0, -1.0, 1.0];
        (p0, p1, p2, p3)
    }
    fn make_particle(pos: [f64; 3], fixed: bool) -> VolumeParticle {
        VolumeParticle {
            position: pos,
            prev_position: pos,
            mass: 1.0,
            fixed,
        }
    }
    #[test]
    fn test_tet_volume_signed_regular() {
        let (p0, p1, p2, p3) = regular_tet();
        let expected = 8.0 / 3.0;
        let vol = tet_volume_signed(p0, p1, p2, p3).abs();
        assert!(
            (vol - expected).abs() < 1e-10,
            "Regular tet volume should be {expected}, got {vol}"
        );
    }
    #[test]
    fn test_solve_edge_constraints_rest_length() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        body.particles[3].position = [0.0, 0.0, 2.0];
        for _ in 0..200 {
            body.solve_edge_constraints();
        }
        let rest_lengths = body.rest_lengths.clone();
        for (i, j, rest) in rest_lengths {
            let cur = len3(sub3(body.particles[i].position, body.particles[j].position));
            assert!(
                (cur - rest).abs() < 0.1,
                "Edge ({i},{j}): current length {cur:.4} should be close to rest {rest:.4}"
            );
        }
    }
    #[test]
    fn test_solve_volume_constraints_corrects_perturbation() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let v_rest = body.rest_volumes[0];
        body.particles[3].position = [0.0, 0.0, 0.2];
        for _ in 0..500 {
            body.solve_volume_constraints();
        }
        let v_cur = TetrahedralSoftBody::tet_volume(&body.particles, &[0, 1, 2, 3]);
        let err = ((v_cur - v_rest) / v_rest).abs();
        assert!(
            err < 0.10,
            "Volume after correction should be close to rest ({v_rest:.4}), got {v_cur:.4} (err={err:.3})"
        );
    }
    #[test]
    fn test_fixed_particle_does_not_move() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut particles: Vec<VolumeParticle> =
            pts.iter().map(|&p| make_particle(p, false)).collect();
        particles[0].fixed = true;
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let initial_pos = body.particles[0].position;
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..60 {
            body.step(1.0 / 60.0, gravity);
        }
        assert_eq!(
            body.particles[0].position, initial_pos,
            "Fixed particle must not move under gravity/constraints"
        );
    }
    #[test]
    fn test_total_volume_equals_sum_of_tets() {
        let particles = vec![
            make_particle([0.0, 0.0, 0.0], false),
            make_particle([1.0, 0.0, 0.0], false),
            make_particle([0.0, 1.0, 0.0], false),
            make_particle([0.0, 0.0, 1.0], false),
            make_particle([1.0, 1.0, 1.0], false),
        ];
        let tetrahedra = vec![[0, 1, 2, 3], [1, 2, 3, 4]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        let total = body.total_volume();
        let sum: f64 = body
            .tetrahedra
            .iter()
            .map(|tet| TetrahedralSoftBody::tet_volume(&body.particles, tet))
            .sum();
        assert!(
            (total - sum).abs() < 1e-14,
            "total_volume() {total} should equal sum of individual volumes {sum}"
        );
    }
    #[test]
    fn test_pressure_volume_work() {
        let w = pressure_volume_work(100.0, 1.5, 1.0);
        assert!(
            (w - 50.0).abs() < 1e-10,
            "pV work = 100 * 0.5 = 50, got {w}"
        );
        let w2 = pressure_volume_work(100.0, 0.8, 1.0);
        assert!(w2 < 0.0, "Compression should give negative work");
    }
    #[test]
    fn test_ideal_gas_pressure() {
        let p = ideal_gas_pressure(100.0, 1.0, 0.5);
        assert!(
            (p - 200.0).abs() < 1e-10,
            "p should be 200 when volume halved, got {p}"
        );
        let p2 = ideal_gas_pressure(100.0, 1.0, 1.0);
        assert!(
            (p2 - 100.0).abs() < 1e-10,
            "p should stay 100 at same volume, got {p2}"
        );
    }
    #[test]
    fn test_adiabatic_gas_pressure() {
        let gamma = 1.4;
        let p = adiabatic_gas_pressure(100.0, 1.0, 0.5, gamma);
        let expected = 100.0 * 2.0_f64.powf(gamma);
        assert!((p - expected).abs() < 1e-8, "Expected {expected}, got {p}");
    }
    #[test]
    fn test_tet_volume_gradient_finite_difference() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let grads = tet_volume_gradient(p0, p1, p2, p3);
        let v0 = tet_volume_signed(p0, p1, p2, p3).abs();
        let eps = 1e-6;
        let p3_shifted = [p3[0], p3[1], p3[2] + eps];
        let v_shifted = tet_volume_signed(p0, p1, p2, p3_shifted).abs();
        let fd = (v_shifted - v0) / eps;
        assert!(
            (fd - grads[3][2]).abs() < 1e-4,
            "FD check: dV/dp3z = {fd}, gradient = {}",
            grads[3][2]
        );
    }
    #[test]
    fn test_volume_constraint_lambda() {
        let grads = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ];
        let inv_masses = [1.0, 1.0, 1.0, 0.0];
        let constraint = 0.5;
        let stiffness = 1.0;
        let lambda = volume_constraint_lambda(&grads, &inv_masses, constraint, stiffness);
        let expected = -0.5 / 3.0;
        assert!(
            (lambda - expected).abs() < 1e-12,
            "Expected {expected}, got {lambda}"
        );
    }
    #[test]
    fn test_solve_volume_with_pressure() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        body.particles[3].position = [0.0, 0.0, 0.2];
        for _ in 0..200 {
            body.solve_volume_with_pressure(100.0);
        }
        let v = TetrahedralSoftBody::tet_volume(&body.particles, &[0, 1, 2, 3]);
        let v_rest = body.rest_volumes[0];
        let err = ((v - v_rest) / v_rest).abs();
        assert!(
            err < 0.2,
            "Pressure should help restore volume: rest={v_rest:.4}, cur={v:.4}"
        );
    }
    #[test]
    fn test_volume_ratios() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        let ratios = body.volume_ratios();
        assert_eq!(ratios.len(), 1);
        assert!(
            (ratios[0] - 1.0).abs() < 1e-10,
            "Initial volume ratio should be 1.0, got {}",
            ratios[0]
        );
    }
    #[test]
    fn test_tet_surface_area() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let area = tet_surface_area(p0, p1, p2, p3);
        let expected = 1.5 + 3.0_f64.sqrt() / 2.0;
        assert!(
            (area - expected).abs() < 1e-10,
            "Expected area {expected:.4}, got {area:.4}"
        );
    }
    #[test]
    fn test_tet_inradius() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let p3 = [0.0, 0.0, 1.0];
        let r = tet_inradius(p0, p1, p2, p3);
        assert!(r > 0.0, "Inradius should be positive");
        assert!(r < 1.0, "Inradius should be less than edge length");
    }
    #[test]
    fn test_step_with_gas_law_stable() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 120.0;
        for _ in 0..60 {
            body.step_with_gas_law(dt, gravity, 100.0);
        }
        for p in &body.particles {
            for &x in &p.position {
                assert!(x.is_finite(), "positions should stay finite");
            }
        }
    }
    #[test]
    fn test_tet_centroid() {
        let positions: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 4.0],
        ];
        let c = tet_centroid(&positions);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_total_volume_centroid_agrees() {
        let pts = [
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_particle(p, false)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        let v1 = body.total_volume();
        let v2 = body.total_volume_centroid();
        assert!(
            (v1 - v2).abs() < v1 * 0.5 + 1e-10,
            "total_volume={v1}, centroid_volume={v2}"
        );
    }
}
