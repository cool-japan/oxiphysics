//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{IbMarker, IbMarker3D};

/// 4-point Peskin (Roma) delta function.
///
/// Returns `(1/4)(1 + cos(π r/2))` for `|r| ≤ 2`, else `0`.
pub fn delta_function(r: f64) -> f64 {
    if r.abs() <= 2.0 {
        0.25 * (1.0 + (PI * r / 2.0).cos())
    } else {
        0.0
    }
}
/// Peskin 4-point regularized delta function.
///
/// The standard Peskin 4-point kernel with compact support \[-2, 2\]:
///
/// ```text
/// phi(r) = (1/8)(3 - 2|r| + sqrt(1 + 4|r| - 4r^2))   if 0 <= |r| <= 1
/// phi(r) = (1/8)(5 - 2|r| - sqrt(-7 + 12|r| - 4r^2))  if 1 <= |r| <= 2
/// phi(r) = 0                                             otherwise
/// ```
pub fn peskin_delta_4pt(r: f64) -> f64 {
    let abs_r = r.abs();
    if abs_r >= 2.0 {
        0.0
    } else if abs_r <= 1.0 {
        0.125 * (3.0 - 2.0 * abs_r + (1.0 + 4.0 * abs_r - 4.0 * abs_r * abs_r).max(0.0).sqrt())
    } else {
        0.125 * (5.0 - 2.0 * abs_r - (-7.0 + 12.0 * abs_r - 4.0 * abs_r * abs_r).max(0.0).sqrt())
    }
}
/// 3-point delta function with compact support \[-1.5, 1.5\].
///
/// ```text
/// phi(r) = (1/3)(1 + sqrt(-3r^2 + 1))   if |r| <= 0.5
/// phi(r) = (1/6)(5 - 3|r| - sqrt(-3(1-|r|)^2 + 1))  if 0.5 <= |r| <= 1.5
/// phi(r) = 0                              otherwise
/// ```
pub fn delta_3pt(r: f64) -> f64 {
    let abs_r = r.abs();
    if abs_r >= 1.5 {
        0.0
    } else if abs_r <= 0.5 {
        (1.0 / 3.0) * (1.0 + (-3.0 * abs_r * abs_r + 1.0).max(0.0).sqrt())
    } else {
        let t = 1.0 - abs_r;
        (1.0 / 6.0) * (5.0 - 3.0 * abs_r - (-3.0 * t * t + 1.0).max(0.0).sqrt())
    }
}
/// Interpolate fluid velocity at each marker position using the delta function.
///
/// # Arguments
/// * `markers` – Lagrangian markers
/// * `fluid_u` – fluid velocity field, indexed as `fluid_u[y * nx + x]`
/// * `nx`, `ny` – grid dimensions
/// * `dx` – grid spacing
pub fn interpolate_fluid_velocity(
    markers: &[IbMarker],
    fluid_u: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<[f64; 2]> {
    let mut result = vec![[0.0_f64; 2]; markers.len()];
    for (m, marker) in markers.iter().enumerate() {
        let xm = marker.position[0];
        let ym = marker.position[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        let mut vel = [0.0_f64; 2];
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let dy_val = delta_function(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let dx_val = delta_function(rx);
                let w = dx_val * dy_val * dx * dx;
                let k = iy as usize * nx + ix as usize;
                vel[0] += w * fluid_u[k][0];
                vel[1] += w * fluid_u[k][1];
            }
        }
        result[m] = vel;
    }
    result
}
/// Interpolate fluid velocity using the Peskin 4-point delta function.
pub fn interpolate_fluid_velocity_peskin(
    markers: &[IbMarker],
    fluid_u: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<[f64; 2]> {
    let mut result = vec![[0.0_f64; 2]; markers.len()];
    for (m, marker) in markers.iter().enumerate() {
        let xm = marker.position[0];
        let ym = marker.position[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        let mut vel = [0.0_f64; 2];
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let dy_val = peskin_delta_4pt(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let dx_val = peskin_delta_4pt(rx);
                let w = dx_val * dy_val * dx * dx;
                let k = iy as usize * nx + ix as usize;
                vel[0] += w * fluid_u[k][0];
                vel[1] += w * fluid_u[k][1];
            }
        }
        result[m] = vel;
    }
    result
}
/// Spread marker forces to the Eulerian fluid grid using the delta function.
///
/// # Arguments
/// * `markers` – Lagrangian markers carrying body forces
/// * `nx`, `ny` – grid dimensions
/// * `dx` – grid spacing
///
/// Returns force-per-unit-volume `[fx, fy]` at each grid node.
pub fn spread_force(markers: &[IbMarker], nx: usize, ny: usize, dx: f64) -> Vec<[f64; 2]> {
    let mut grid_force = vec![[0.0_f64; 2]; nx * ny];
    for marker in markers {
        let xm = marker.position[0];
        let ym = marker.position[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let dy_val = delta_function(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let dx_val = delta_function(rx);
                let w = dx_val * dy_val;
                let k = iy as usize * nx + ix as usize;
                grid_force[k][0] += w * marker.force[0];
                grid_force[k][1] += w * marker.force[1];
            }
        }
    }
    grid_force
}
/// Spread marker forces using the Peskin 4-point delta function.
pub fn spread_force_peskin(markers: &[IbMarker], nx: usize, ny: usize, dx: f64) -> Vec<[f64; 2]> {
    let mut grid_force = vec![[0.0_f64; 2]; nx * ny];
    for marker in markers {
        let xm = marker.position[0];
        let ym = marker.position[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let dy_val = peskin_delta_4pt(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let dx_val = peskin_delta_4pt(rx);
                let w = dx_val * dy_val;
                let k = iy as usize * nx + ix as usize;
                grid_force[k][0] += w * marker.force[0];
                grid_force[k][1] += w * marker.force[1];
            }
        }
    }
    grid_force
}
/// Move markers by their current velocity over time step `dt`.
pub fn update_marker_positions(markers: &mut [IbMarker], dt: f64) {
    for marker in markers.iter_mut() {
        marker.position[0] += marker.velocity[0] * dt;
        marker.position[1] += marker.velocity[1] * dt;
    }
}
/// Compute spring forces on markers toward their rest positions.
///
/// `F = -stiffness * (x - x_rest)` for each component.
pub fn compute_elastic_forces(
    markers: &mut [IbMarker],
    rest_positions: &[[f64; 2]],
    stiffness: f64,
) {
    for (m, marker) in markers.iter_mut().enumerate() {
        if m < rest_positions.len() {
            marker.force[0] = -stiffness * (marker.position[0] - rest_positions[m][0]);
            marker.force[1] = -stiffness * (marker.position[1] - rest_positions[m][1]);
        }
    }
}
/// Compute membrane forces including tension and bending.
///
/// The membrane tension force uses a spring model between consecutive markers.
/// The bending force penalizes curvature changes between consecutive segments.
///
/// # Arguments
/// * `markers` – Lagrangian markers (ordered along the membrane)
/// * `rest_lengths` – rest length between consecutive marker pairs
/// * `tension_stiffness` – spring constant for tension
/// * `bending_stiffness` – spring constant for bending resistance
pub fn compute_membrane_forces(
    markers: &mut [IbMarker],
    rest_lengths: &[f64],
    tension_stiffness: f64,
    bending_stiffness: f64,
) {
    let n = markers.len();
    if n < 2 {
        return;
    }
    for m in markers.iter_mut() {
        m.force = [0.0, 0.0];
    }
    for i in 0..n - 1 {
        let dx = markers[i + 1].position[0] - markers[i].position[0];
        let dy = markers[i + 1].position[1] - markers[i].position[1];
        let length = (dx * dx + dy * dy).sqrt();
        if length < 1e-14 {
            continue;
        }
        let rest = if i < rest_lengths.len() {
            rest_lengths[i]
        } else {
            length
        };
        let strain = length - rest;
        let fx = tension_stiffness * strain * dx / length;
        let fy = tension_stiffness * strain * dy / length;
        markers[i].force[0] += fx;
        markers[i].force[1] += fy;
        markers[i + 1].force[0] -= fx;
        markers[i + 1].force[1] -= fy;
    }
    if n >= 3 {
        for i in 1..n - 1 {
            let bx = markers[i - 1].position[0] - 2.0 * markers[i].position[0]
                + markers[i + 1].position[0];
            let by = markers[i - 1].position[1] - 2.0 * markers[i].position[1]
                + markers[i + 1].position[1];
            markers[i].force[0] += bending_stiffness * bx;
            markers[i].force[1] += bending_stiffness * by;
        }
    }
}
/// Compute direct forcing (penalty method) forces on markers.
///
/// The force is computed to enforce a desired velocity at the IB:
///
/// `F = (u_desired - u_interpolated) * rho / dt`
///
/// This is used for rigid-body IB where the boundary velocity is known.
pub fn compute_penalty_forces(
    markers: &mut [IbMarker],
    desired_velocity: &[[f64; 2]],
    interpolated_velocity: &[[f64; 2]],
    rho: f64,
    dt: f64,
) {
    for (i, marker) in markers.iter_mut().enumerate() {
        if i < desired_velocity.len() && i < interpolated_velocity.len() {
            marker.force[0] = rho * (desired_velocity[i][0] - interpolated_velocity[i][0]) / dt;
            marker.force[1] = rho * (desired_velocity[i][1] - interpolated_velocity[i][1]) / dt;
        }
    }
}
/// Add body force to a single LBM distribution function using Guo's scheme.
///
/// `fi += (1 - 1/(2τ)) * wi * (ci · F) / cs² * dt`
///
/// # Arguments
/// * `_grid_force` – force vector at this node (not used here; caller passes `ci · F` externally)
/// * `f_eq` – equilibrium component `wi * (ci · F) / cs²`
/// * `dt` – time step
pub fn ib_forcing_term(_grid_force: &[[f64; 2]], f_eq: f64, dt: f64) -> f64 {
    let factor = 0.5;
    factor * f_eq * dt
}
/// 3D delta function (tensor product of 1D functions).
///
/// δ³(r) = δ(rx) * δ(ry) * δ(rz) where δ is the 1D Roma kernel.
pub fn delta_function_3d(r: [f64; 3]) -> f64 {
    delta_function(r[0]) * delta_function(r[1]) * delta_function(r[2])
}
/// Spread 3D IBM force to Eulerian grid using the Roma delta function.
///
/// Returns a flat `Vec<[f64; 3]>` of length `nx * ny * nz`, indexed as
/// `k = z * ny * nx + y * nx + x`.
pub fn spread_force_3d(
    markers: &[IbMarker3D],
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
) -> Vec<[f64; 3]> {
    let mut grid_force = vec![[0.0_f64; 3]; nx * ny * nz];
    for m in markers {
        let xi = (m.position[0] / dx).floor() as isize;
        let yi = (m.position[1] / dx).floor() as isize;
        let zi = (m.position[2] / dx).floor() as isize;
        for iz in (zi - 2)..=(zi + 2) {
            for iy in (yi - 2)..=(yi + 2) {
                for ix in (xi - 2)..=(xi + 2) {
                    if ix < 0
                        || iy < 0
                        || iz < 0
                        || ix >= nx as isize
                        || iy >= ny as isize
                        || iz >= nz as isize
                    {
                        continue;
                    }
                    let rx = (m.position[0] - ix as f64 * dx) / dx;
                    let ry = (m.position[1] - iy as f64 * dx) / dx;
                    let rz = (m.position[2] - iz as f64 * dx) / dx;
                    let d = delta_function_3d([rx, ry, rz]);
                    let idx = iz as usize * ny * nx + iy as usize * nx + ix as usize;
                    for (gf, mf) in grid_force[idx].iter_mut().zip(m.force.iter()) {
                        *gf += mf * d / (dx * dx * dx);
                    }
                }
            }
        }
    }
    grid_force
}
/// Interpolate 3D fluid velocity to marker positions.
pub fn interpolate_fluid_velocity_3d(
    markers: &[IbMarker3D],
    fluid_u: &[[f64; 3]],
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
) -> Vec<[f64; 3]> {
    let mut marker_u = vec![[0.0_f64; 3]; markers.len()];
    for (k, m) in markers.iter().enumerate() {
        let xi = (m.position[0] / dx).floor() as isize;
        let yi = (m.position[1] / dx).floor() as isize;
        let zi = (m.position[2] / dx).floor() as isize;
        for iz in (zi - 2)..=(zi + 2) {
            for iy in (yi - 2)..=(yi + 2) {
                for ix in (xi - 2)..=(xi + 2) {
                    if ix < 0
                        || iy < 0
                        || iz < 0
                        || ix >= nx as isize
                        || iy >= ny as isize
                        || iz >= nz as isize
                    {
                        continue;
                    }
                    let rx = (m.position[0] - ix as f64 * dx) / dx;
                    let ry = (m.position[1] - iy as f64 * dx) / dx;
                    let rz = (m.position[2] - iz as f64 * dx) / dx;
                    let d = delta_function_3d([rx, ry, rz]);
                    let idx = iz as usize * ny * nx + iy as usize * nx + ix as usize;
                    for dim in 0..3 {
                        marker_u[k][dim] += fluid_u[idx][dim] * d * dx * dx * dx;
                    }
                }
            }
        }
    }
    marker_u
}
/// Generate Lagrangian markers for a 3D sphere surface.
///
/// Uses a Fibonacci lattice for approximately uniform coverage.
/// Returns `n` markers on a sphere of radius `r` centred at `(cx, cy, cz)`.
pub fn sphere_markers_3d(cx: f64, cy: f64, cz: f64, r: f64, n: usize) -> Vec<IbMarker3D> {
    let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
    (0..n)
        .map(|i| {
            let theta = std::f64::consts::PI * 2.0 * i as f64 / golden_ratio;
            let phi = (1.0 - 2.0 * (i as f64 + 0.5) / n as f64).acos();
            let x = cx + r * phi.sin() * theta.cos();
            let y = cy + r * phi.sin() * theta.sin();
            let z = cz + r * phi.cos();
            IbMarker3D::new(x, y, z)
        })
        .collect()
}
/// Compute the area element ds for a set of markers on a sphere of radius r.
///
/// ds = 4π r² / n (uniform approximation).
pub fn sphere_marker_area(r: f64, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    4.0 * std::f64::consts::PI * r * r / n as f64
}
/// Generate markers for a cylinder aligned with z-axis.
///
/// Returns `n_theta * n_z` markers uniformly distributed on the cylinder
/// surface from z=z0 to z=z1, radius r, centred at (cx, cy).
pub fn cylinder_markers_3d(
    cx: f64,
    cy: f64,
    z0: f64,
    z1: f64,
    r: f64,
    n_theta: usize,
    n_z: usize,
) -> Vec<IbMarker3D> {
    let mut markers = Vec::with_capacity(n_theta * n_z);
    for iz in 0..n_z {
        let z = z0 + (z1 - z0) * iz as f64 / (n_z.max(2) - 1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            markers.push(IbMarker3D::new(
                cx + r * theta.cos(),
                cy + r * theta.sin(),
                z,
            ));
        }
    }
    markers
}
/// Compute lift and drag from IBM grid forces.
///
/// Given the Eulerian body force field `f_grid` (length nx*ny, units force/volume)
/// and a reference velocity direction `(cos_alpha, sin_alpha)`, returns
/// `(drag, lift)` in physical units (force density summed over all cells * dx²).
pub fn compute_lift_drag(
    f_grid: &[[f64; 2]],
    cos_alpha: f64,
    sin_alpha: f64,
    dx: f64,
) -> (f64, f64) {
    let area = dx * dx;
    let mut drag = 0.0;
    let mut lift = 0.0;
    for f in f_grid {
        drag += (f[0] * cos_alpha + f[1] * sin_alpha) * area;
        lift += (-f[0] * sin_alpha + f[1] * cos_alpha) * area;
    }
    (drag, lift)
}
/// Compute lift coefficient: CL = 2 * L / (ρ * U² * A).
pub fn lift_coefficient(lift: f64, rho: f64, u_inf: f64, area: f64) -> f64 {
    if rho.abs() < 1e-30 || u_inf.abs() < 1e-30 || area.abs() < 1e-30 {
        return 0.0;
    }
    2.0 * lift / (rho * u_inf * u_inf * area)
}
/// Compute drag coefficient: CD = 2 * D / (ρ * U² * A).
pub fn drag_coefficient(drag: f64, rho: f64, u_inf: f64, area: f64) -> f64 {
    if rho.abs() < 1e-30 || u_inf.abs() < 1e-30 || area.abs() < 1e-30 {
        return 0.0;
    }
    2.0 * drag / (rho * u_inf * u_inf * area)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::immersed_boundary::types::*;
    #[test]
    fn test_delta_at_zero() {
        let d = delta_function(0.0);
        assert!((d - 0.5).abs() < 1e-14, "delta(0) = {d}");
    }
    #[test]
    fn test_delta_at_boundary() {
        assert!(delta_function(2.0).abs() < 1e-14, "delta(2) should be 0");
        assert!(delta_function(-2.0).abs() < 1e-14, "delta(-2) should be 0");
    }
    #[test]
    fn test_delta_outside_support() {
        assert_eq!(delta_function(3.0), 0.0);
        assert_eq!(delta_function(-3.0), 0.0);
        assert_eq!(delta_function(100.0), 0.0);
    }
    #[test]
    fn test_delta_symmetric() {
        for &r in &[0.5, 1.0, 1.5, 1.9] {
            let pos = delta_function(r);
            let neg = delta_function(-r);
            assert!((pos - neg).abs() < 1e-14, "delta not symmetric at r={r}");
        }
    }
    #[test]
    fn test_delta_at_one() {
        let d = delta_function(1.0);
        assert!((d - 0.25).abs() < 1e-14, "delta(1) = {d}");
    }
    #[test]
    fn test_delta_partition_of_unity_approx() {
        let offset = 0.3_f64;
        let sum: f64 = (-3..=3).map(|i| delta_function(i as f64 - offset)).sum();
        assert!((sum - 1.0).abs() < 0.05, "partition of unity: sum={sum}");
    }
    #[test]
    fn test_ib_marker_defaults() {
        let m = IbMarker::new(1.0, 2.0);
        assert_eq!(m.velocity, [0.0, 0.0]);
        assert_eq!(m.force, [0.0, 0.0]);
        assert_eq!(m.position, [1.0, 2.0]);
    }
    #[test]
    fn test_ibm_config() {
        let cfg = IbmConfig::new(0.1, 0.01, 1.0);
        assert!((cfg.delta_x - 0.1).abs() < 1e-14);
        assert!((cfg.delta_t - 0.01).abs() < 1e-14);
        assert!((cfg.rho_f - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_ib_circle_markers_on_circle() {
        let circ = IbCircle::new(5.0, 5.0, 2.0, 8, 1.0);
        let markers = circ.to_markers();
        assert_eq!(markers.len(), 8);
        for m in &markers {
            let dx = m.position[0] - 5.0;
            let dy = m.position[1] - 5.0;
            let r = (dx * dx + dy * dy).sqrt();
            assert!((r - 2.0).abs() < 1e-12, "marker not on circle: r={r}");
        }
    }
    #[test]
    fn test_ib_circle_first_marker() {
        let circ = IbCircle::new(0.0, 0.0, 3.0, 4, 1.0);
        let markers = circ.to_markers();
        assert!((markers[0].position[0] - 3.0).abs() < 1e-12);
        assert!((markers[0].position[1]).abs() < 1e-12);
    }
    #[test]
    fn test_update_marker_positions() {
        let mut markers = vec![IbMarker {
            position: [1.0, 2.0],
            velocity: [0.5, -0.3],
            force: [0.0, 0.0],
        }];
        update_marker_positions(&mut markers, 0.1);
        assert!((markers[0].position[0] - 1.05).abs() < 1e-14);
        assert!((markers[0].position[1] - 1.97).abs() < 1e-14);
    }
    #[test]
    fn test_elastic_force_at_rest() {
        let mut markers = vec![IbMarker::new(1.0, 2.0)];
        let rest = vec![[1.0, 2.0]];
        compute_elastic_forces(&mut markers, &rest, 100.0);
        assert!(markers[0].force[0].abs() < 1e-14);
        assert!(markers[0].force[1].abs() < 1e-14);
    }
    #[test]
    fn test_elastic_force_displacement() {
        let mut markers = vec![IbMarker::new(2.0, 3.0)];
        let rest = vec![[1.0, 2.0]];
        compute_elastic_forces(&mut markers, &rest, 10.0);
        assert!((markers[0].force[0] - (-10.0)).abs() < 1e-14);
        assert!((markers[0].force[1] - (-10.0)).abs() < 1e-14);
    }
    #[test]
    fn test_spread_force_single_marker() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let mut marker = IbMarker::new(4.0, 4.0);
        marker.force = [1.0, 0.0];
        let grid_force = spread_force(&[marker], nx, ny, dx);
        let total_fx: f64 = grid_force.iter().map(|f| f[0]).sum();
        assert!(
            total_fx > 0.0,
            "spread force should be positive: {total_fx}"
        );
    }
    #[test]
    fn test_spread_force_zero_force() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let marker = IbMarker::new(4.0, 4.0);
        let grid_force = spread_force(&[marker], nx, ny, dx);
        let total: f64 = grid_force.iter().map(|f| f[0] + f[1]).sum();
        assert!(total.abs() < 1e-15);
    }
    #[test]
    fn test_interpolate_uniform_field() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let fluid_u = vec![[0.3_f64, 0.1_f64]; nx * ny];
        let markers = vec![IbMarker::new(3.5, 3.5)];
        let interp = interpolate_fluid_velocity(&markers, &fluid_u, nx, ny, dx);
        assert!(
            (interp[0][0] - 0.3).abs() < 0.01,
            "interp ux = {}",
            interp[0][0]
        );
        assert!(
            (interp[0][1] - 0.1).abs() < 0.01,
            "interp uy = {}",
            interp[0][1]
        );
    }
    #[test]
    fn test_interpolate_zero_field() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let fluid_u = vec![[0.0_f64, 0.0_f64]; nx * ny];
        let markers = vec![IbMarker::new(4.0, 4.0)];
        let interp = interpolate_fluid_velocity(&markers, &fluid_u, nx, ny, dx);
        assert!(interp[0][0].abs() < 1e-15);
        assert!(interp[0][1].abs() < 1e-15);
    }
    #[test]
    fn test_ib_forcing_linear() {
        let dummy_force: Vec<[f64; 2]> = vec![];
        let dt = 0.1;
        let a = ib_forcing_term(&dummy_force, 1.0, dt);
        let b = ib_forcing_term(&dummy_force, 2.0, dt);
        assert!((b - 2.0 * a).abs() < 1e-14, "forcing not linear");
    }
    #[test]
    fn test_ib_forcing_dt_linear() {
        let dummy_force: Vec<[f64; 2]> = vec![];
        let a = ib_forcing_term(&dummy_force, 1.0, 0.1);
        let b = ib_forcing_term(&dummy_force, 1.0, 0.2);
        assert!((b - 2.0 * a).abs() < 1e-14, "forcing not linear in dt");
    }
    #[test]
    fn test_ibm_stats_zero_force() {
        let markers = vec![IbMarker::new(0.0, 0.0), IbMarker::new(1.0, 0.0)];
        let rest = vec![[0.0, 0.0], [1.0, 0.0]];
        let stats = IbmStats::compute(&markers, &rest);
        assert!(stats.total_force[0].abs() < 1e-14);
        assert!(stats.total_force[1].abs() < 1e-14);
        assert!(stats.max_marker_displacement < 1e-14);
    }
    #[test]
    fn test_ibm_stats_total_force() {
        let mut m1 = IbMarker::new(0.0, 0.0);
        m1.force = [1.0, 2.0];
        let mut m2 = IbMarker::new(1.0, 0.0);
        m2.force = [3.0, -1.0];
        let markers = vec![m1, m2];
        let rest = vec![[0.0, 0.0], [1.0, 0.0]];
        let stats = IbmStats::compute(&markers, &rest);
        assert!((stats.total_force[0] - 4.0).abs() < 1e-14);
        assert!((stats.total_force[1] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_ibm_stats_max_displacement() {
        let markers = vec![IbMarker::new(1.0, 0.0), IbMarker::new(3.0, 4.0)];
        let rest = vec![[0.0, 0.0], [0.0, 0.0]];
        let stats = IbmStats::compute(&markers, &rest);
        assert!((stats.max_marker_displacement - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_ib_circle_single_marker() {
        let circ = IbCircle::new(2.0, 3.0, 1.5, 1, 1.0);
        let markers = circ.to_markers();
        assert_eq!(markers.len(), 1);
        assert!((markers[0].position[0] - 3.5).abs() < 1e-12);
        assert!((markers[0].position[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_elastic_force_direction() {
        let mut markers = vec![IbMarker::new(5.0, 0.0)];
        let rest = vec![[0.0, 0.0]];
        compute_elastic_forces(&mut markers, &rest, 1.0);
        assert!(markers[0].force[0] < 0.0, "force should be leftward");
    }
    #[test]
    fn test_spread_force_y_component() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let mut marker = IbMarker::new(4.0, 4.0);
        marker.force = [0.0, 2.0];
        let grid_force = spread_force(&[marker], nx, ny, dx);
        let total_fy: f64 = grid_force.iter().map(|f| f[1]).sum();
        assert!(
            total_fy > 0.0,
            "spread y-force should be positive: {total_fy}"
        );
        let total_fx: f64 = grid_force.iter().map(|f| f[0]).sum();
        assert!(
            total_fx.abs() < 1e-15,
            "x-component should be zero: {total_fx}"
        );
    }
    #[test]
    fn test_ib_circle_markers_zero_vel_force() {
        let circ = IbCircle::new(0.0, 0.0, 1.0, 4, 10.0);
        for m in circ.to_markers() {
            assert_eq!(m.velocity, [0.0, 0.0]);
            assert_eq!(m.force, [0.0, 0.0]);
        }
    }
    #[test]
    fn test_delta_at_half() {
        let r = 0.5;
        let expected = 0.25 * (1.0 + (PI * r / 2.0).cos());
        let got = delta_function(r);
        assert!((got - expected).abs() < 1e-14, "delta(0.5) = {got}");
    }
    #[test]
    fn test_spread_force_linearity() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let mut m1 = IbMarker::new(4.0, 4.0);
        m1.force = [1.0, 0.0];
        let mut m2 = IbMarker::new(4.0, 4.0);
        m2.force = [1.0, 0.0];
        let g1 = spread_force(&[m1.clone()], nx, ny, dx);
        let g2 = spread_force(&[m1, m2], nx, ny, dx);
        let sum1: f64 = g1.iter().map(|f| f[0]).sum();
        let sum2: f64 = g2.iter().map(|f| f[0]).sum();
        assert!(
            (sum2 - 2.0 * sum1).abs() < 1e-12,
            "linearity failed: {sum1} {sum2}"
        );
    }
    #[test]
    fn test_peskin_delta_at_zero() {
        let d = peskin_delta_4pt(0.0);
        assert!(d > 0.0, "peskin delta at 0 should be positive");
        assert!(
            (d - 0.5).abs() < 1e-12,
            "peskin delta(0) = {d}, expected 0.5"
        );
    }
    #[test]
    fn test_peskin_delta_symmetric() {
        for &r in &[0.3, 0.8, 1.2, 1.7] {
            let pos = peskin_delta_4pt(r);
            let neg = peskin_delta_4pt(-r);
            assert!(
                (pos - neg).abs() < 1e-14,
                "peskin delta not symmetric at r={r}"
            );
        }
    }
    #[test]
    fn test_peskin_delta_outside() {
        assert_eq!(peskin_delta_4pt(2.5), 0.0);
        assert_eq!(peskin_delta_4pt(-3.0), 0.0);
    }
    #[test]
    fn test_delta_3pt_at_zero() {
        let d = delta_3pt(0.0);
        assert!((d - 2.0 / 3.0).abs() < 1e-12, "3pt delta(0) = {d}");
    }
    #[test]
    fn test_delta_3pt_symmetric() {
        for &r in &[0.2, 0.7, 1.2] {
            let pos = delta_3pt(r);
            let neg = delta_3pt(-r);
            assert!(
                (pos - neg).abs() < 1e-14,
                "3pt delta not symmetric at r={r}"
            );
        }
    }
    #[test]
    fn test_delta_3pt_outside() {
        assert_eq!(delta_3pt(2.0), 0.0);
        assert_eq!(delta_3pt(-2.0), 0.0);
    }
    #[test]
    fn test_membrane_at_rest() {
        let mut markers = vec![IbMarker::new(0.0, 0.0), IbMarker::new(1.0, 0.0)];
        let rest_lengths = vec![1.0];
        compute_membrane_forces(&mut markers, &rest_lengths, 10.0, 1.0);
        for m in &markers {
            assert!(m.force[0].abs() < 1e-12, "force should be ~0 at rest");
            assert!(m.force[1].abs() < 1e-12, "force should be ~0 at rest");
        }
    }
    #[test]
    fn test_membrane_tension_stretched() {
        let mut markers = vec![IbMarker::new(0.0, 0.0), IbMarker::new(2.0, 0.0)];
        let rest_lengths = vec![1.0];
        compute_membrane_forces(&mut markers, &rest_lengths, 10.0, 0.0);
        assert!(markers[0].force[0] > 0.0, "marker 0 should be pulled right");
        assert!(markers[1].force[0] < 0.0, "marker 1 should be pulled left");
        assert!(
            (markers[0].force[0] + markers[1].force[0]).abs() < 1e-12,
            "Newton III violation"
        );
    }
    #[test]
    fn test_penalty_forces() {
        let mut markers = vec![IbMarker::new(1.0, 0.0)];
        let desired = vec![[0.5, 0.0]];
        let interpolated = vec![[0.3, 0.0]];
        let rho = 1.0;
        let dt = 0.01;
        compute_penalty_forces(&mut markers, &desired, &interpolated, rho, dt);
        assert!(
            (markers[0].force[0] - 20.0).abs() < 1e-10,
            "penalty force x = {}",
            markers[0].force[0]
        );
    }
    #[test]
    fn test_multi_ib_total_markers() {
        let mut sys = MultiIbSystem::new();
        sys.add_circle(&IbCircle::new(5.0, 5.0, 1.0, 10, 1.0));
        sys.add_circle(&IbCircle::new(10.0, 5.0, 1.0, 8, 1.0));
        assert_eq!(sys.total_markers(), 18);
        assert_eq!(sys.bodies.len(), 2);
    }
    #[test]
    fn test_multi_ib_spread_forces() {
        let mut sys = MultiIbSystem::new();
        sys.add_circle(&IbCircle::new(4.0, 4.0, 1.0, 4, 1.0));
        for m in sys.bodies[0].iter_mut() {
            m.position[0] += 0.1;
        }
        sys.compute_all_forces();
        let grid_force = sys.spread_all_forces(8, 8, 1.0);
        let total_fx: f64 = grid_force.iter().map(|f| f[0]).sum();
        assert!(
            total_fx.abs() > 0.0,
            "multi-IB should produce nonzero grid force"
        );
    }
    #[test]
    fn test_peskin_interpolate_uniform() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let fluid_u = vec![[0.2_f64, 0.05_f64]; nx * ny];
        let markers = vec![IbMarker::new(3.5, 3.5)];
        let interp = interpolate_fluid_velocity_peskin(&markers, &fluid_u, nx, ny, dx);
        assert!(
            (interp[0][0] - 0.2).abs() < 0.02,
            "peskin interp ux = {}",
            interp[0][0]
        );
        assert!(
            (interp[0][1] - 0.05).abs() < 0.02,
            "peskin interp uy = {}",
            interp[0][1]
        );
    }
    #[test]
    fn test_peskin_spread_force() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let mut marker = IbMarker::new(4.0, 4.0);
        marker.force = [1.0, 0.0];
        let grid_force = spread_force_peskin(&[marker], nx, ny, dx);
        let total_fx: f64 = grid_force.iter().map(|f| f[0]).sum();
        assert!(total_fx > 0.0, "peskin spread force should be positive");
    }
}
#[cfg(test)]
mod tests_extended_ibm {
    use super::*;
    use crate::immersed_boundary::types::*;
    #[test]
    fn test_ib_marker_3d_new() {
        let m = IbMarker3D::new(1.0, 2.0, 3.0);
        assert_eq!(m.position, [1.0, 2.0, 3.0]);
        assert_eq!(m.velocity, [0.0, 0.0, 0.0]);
        assert_eq!(m.force, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_ib_marker_3d_distance() {
        let a = IbMarker3D::new(0.0, 0.0, 0.0);
        let b = IbMarker3D::new(1.0, 0.0, 0.0);
        let d = a.distance_to(&b);
        assert!((d - 1.0).abs() < 1e-14, "distance = {d}");
    }
    #[test]
    fn test_delta_function_3d_origin() {
        let d = delta_function_3d([0.0, 0.0, 0.0]);
        assert!((d - 0.125).abs() < 1e-14, "delta_3d(0,0,0) = {d}");
    }
    #[test]
    fn test_delta_function_3d_outside() {
        let d = delta_function_3d([3.0, 0.0, 0.0]);
        assert_eq!(d, 0.0, "outside support → 0");
    }
    #[test]
    fn test_sphere_markers_3d_count() {
        let markers = sphere_markers_3d(0.0, 0.0, 0.0, 1.0, 20);
        assert_eq!(markers.len(), 20, "should have 20 markers");
    }
    #[test]
    fn test_sphere_markers_3d_radius() {
        let r = 0.5;
        let markers = sphere_markers_3d(0.0, 0.0, 0.0, r, 50);
        for m in &markers {
            let dist = (m.position[0] * m.position[0]
                + m.position[1] * m.position[1]
                + m.position[2] * m.position[2])
                .sqrt();
            assert!((dist - r).abs() < 1e-10, "radius = {dist}, expected {r}");
        }
    }
    #[test]
    fn test_sphere_marker_area() {
        let area = sphere_marker_area(1.0, 100);
        let expected = 4.0 * std::f64::consts::PI / 100.0;
        assert!((area - expected).abs() < 1e-13, "area = {area}");
    }
    #[test]
    fn test_cylinder_markers_count() {
        let markers = cylinder_markers_3d(0.0, 0.0, 0.0, 1.0, 0.5, 8, 4);
        assert_eq!(markers.len(), 32, "8 × 4 = 32 markers");
    }
    #[test]
    fn test_feedback_forcing_proportional() {
        let mut ff = FeedbackForcing::new(100.0, 0.0, 1);
        let u_interp = [[0.0, 0.0]];
        let u_target = [[1.0, 0.0]];
        let forces = ff.compute(&u_interp, &u_target, 0.01);
        assert!(
            (forces[0][0] - 100.0).abs() < 1e-10,
            "F_x = {}",
            forces[0][0]
        );
    }
    #[test]
    fn test_feedback_forcing_zero_error() {
        let mut ff = FeedbackForcing::new(100.0, 50.0, 1);
        let u = [[0.3_f64, 0.1_f64]];
        let forces = ff.compute(&u, &u, 0.01);
        assert!(forces[0][0].abs() < 1e-14, "zero error → zero force");
        assert!(forces[0][1].abs() < 1e-14, "zero error → zero force y");
    }
    #[test]
    fn test_feedback_forcing_reset() {
        let mut ff = FeedbackForcing::new(100.0, 10.0, 2);
        ff.integral[0] = [5.0, 3.0];
        ff.integral[1] = [-1.0, 2.0];
        ff.reset();
        for s in &ff.integral {
            assert_eq!(*s, [0.0, 0.0], "integral should be zero after reset");
        }
    }
    #[test]
    fn test_lift_drag_uniform_force_along_x() {
        let f_grid = vec![[1.0_f64, 0.0_f64]; 4];
        let (drag, lift) = compute_lift_drag(&f_grid, 1.0, 0.0, 1.0);
        assert!((drag - 4.0).abs() < 1e-12, "drag = {drag}");
        assert!(lift.abs() < 1e-12, "lift = {lift}");
    }
    #[test]
    fn test_lift_coefficient_formula() {
        let cl = lift_coefficient(1.0, 1.0, 1.0, 1.0);
        assert!((cl - 2.0).abs() < 1e-14, "CL = {cl}");
    }
    #[test]
    fn test_drag_coefficient_zero_velocity() {
        let cd = drag_coefficient(1.0, 1.2, 0.0, 1.0);
        assert_eq!(cd, 0.0, "CD with zero velocity");
    }
    #[test]
    fn test_spread_force_3d_nonzero() {
        let mut m = IbMarker3D::new(4.0, 4.0, 4.0);
        m.force = [1.0, 0.0, 0.0];
        let f = spread_force_3d(&[m], 8, 8, 8, 1.0);
        let total_fx: f64 = f.iter().map(|v| v[0]).sum();
        assert!(
            total_fx > 0.0,
            "3D spread force should be nonzero: {total_fx}"
        );
    }
}
/// 2D Peskin regularized delta function (product of 1D Peskin kernels).
///
/// δ₂(r) = δ_1D(rx) * δ_1D(ry) using the 4-point Roma kernel.
pub fn peskin_delta_2d(rx: f64, ry: f64) -> f64 {
    peskin_delta_4pt(rx) * peskin_delta_4pt(ry)
}
/// 3D Peskin regularized delta function (product of 1D Peskin kernels).
pub fn peskin_delta_3d(rx: f64, ry: f64, rz: f64) -> f64 {
    peskin_delta_4pt(rx) * peskin_delta_4pt(ry) * peskin_delta_4pt(rz)
}
/// 1D "3-point" Yang delta function (compact support \[-1.5, 1.5\]).
///
/// Slightly different from delta_3pt; uses a cos-based kernel:
/// φ(r) = (1 + cos(π|r|/1.5)) / 3  for |r| ≤ 1.5, else 0.
pub fn yang_delta_1d(r: f64) -> f64 {
    let abs_r = r.abs();
    if abs_r >= 1.5 {
        0.0
    } else {
        (1.0 + (PI * abs_r / 1.5).cos()) / 3.0
    }
}
