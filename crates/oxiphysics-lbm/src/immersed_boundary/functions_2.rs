//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::peskin_delta_4pt;
use std::f64::consts::PI;

use super::types::IbmBody;

/// Uhlmann direct-forcing spread.
///
/// Distributes marker forces `F_k` to fluid nodes using a given 1D delta
/// kernel `phi`. Returns the Eulerian body force `f(x)` as a flat Vec.
///
/// Uses a 4-cell stencil (same as `spread_force`).
pub fn spread_force_uhlmann(
    marker_positions: &[[f64; 2]],
    marker_forces: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<[f64; 2]> {
    assert_eq!(marker_positions.len(), marker_forces.len());
    let mut grid_force = vec![[0.0_f64; 2]; nx * ny];
    for (pos, force) in marker_positions.iter().zip(marker_forces.iter()) {
        let xm = pos[0];
        let ym = pos[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let phy = peskin_delta_4pt(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let phx = peskin_delta_4pt(rx);
                let w = phx * phy;
                let k = iy as usize * nx + ix as usize;
                grid_force[k][0] += w * force[0];
                grid_force[k][1] += w * force[1];
            }
        }
    }
    grid_force
}
/// Interpolate fluid velocity at arbitrary positions using the Peskin kernel.
///
/// Separate from `interpolate_fluid_velocity` — takes raw position arrays
/// rather than `IbMarker` slices, for use with `IbmBody`.
pub fn interpolate_velocity_at_positions(
    positions: &[[f64; 2]],
    fluid_u: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
) -> Vec<[f64; 2]> {
    let mut result = vec![[0.0_f64; 2]; positions.len()];
    for (k, pos) in positions.iter().enumerate() {
        let xm = pos[0];
        let ym = pos[1];
        let ix0 = ((xm / dx).floor() as i64) - 1;
        let iy0 = ((ym / dx).floor() as i64) - 1;
        let mut vel = [0.0_f64; 2];
        for jy in 0..4_i64 {
            let iy = iy0 + jy;
            if iy < 0 || iy >= ny as i64 {
                continue;
            }
            let ry = (ym - iy as f64 * dx) / dx;
            let phy = peskin_delta_4pt(ry);
            for jx in 0..4_i64 {
                let ix = ix0 + jx;
                if ix < 0 || ix >= nx as i64 {
                    continue;
                }
                let rx = (xm - ix as f64 * dx) / dx;
                let phx = peskin_delta_4pt(rx);
                let w = phx * phy * dx * dx;
                let idx = iy as usize * nx + ix as usize;
                vel[0] += w * fluid_u[idx][0];
                vel[1] += w * fluid_u[idx][1];
            }
        }
        result[k] = vel;
    }
    result
}
/// Uhlmann direct forcing step.
///
/// Given the desired boundary velocity `u_b` and the interpolated fluid
/// velocity `u_f` at marker positions, computes the IBM forcing:
///
/// F_k = (u_b_k - u_f_k) / dt  (per unit volume after spreading)
///
/// Returns a vector of forces, one per marker.
pub fn uhlmann_direct_forcing(
    u_boundary: &[[f64; 2]],
    u_fluid_at_boundary: &[[f64; 2]],
    dt: f64,
) -> Vec<[f64; 2]> {
    assert_eq!(u_boundary.len(), u_fluid_at_boundary.len());
    u_boundary
        .iter()
        .zip(u_fluid_at_boundary.iter())
        .map(|(&ub, &uf)| [(ub[0] - uf[0]) / dt, (ub[1] - uf[1]) / dt])
        .collect()
}
/// Apply Uhlmann IBM: full step for one body.
///
/// 1. Interpolate fluid velocity to boundary points.
/// 2. Compute direct forcing forces.
/// 3. Spread forces back to fluid grid.
///
/// Returns the Eulerian forcing field `f(x)`.
pub fn uhlmann_ibm_step(
    body_positions: &[[f64; 2]],
    body_velocities: &[[f64; 2]],
    fluid_u: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
    dt: f64,
) -> Vec<[f64; 2]> {
    let u_interp = interpolate_velocity_at_positions(body_positions, fluid_u, nx, ny, dx);
    let forces = uhlmann_direct_forcing(body_velocities, &u_interp, dt);
    spread_force_uhlmann(body_positions, &forces, nx, ny, dx)
}
/// Advance an IbmBody position using its current velocity (explicit Euler).
///
/// Separate from `IbmBody::advance` — does not apply forces, just kinematics.
pub fn advance_body_kinematics(body: &mut IbmBody, dt: f64) {
    body.centroid[0] += body.velocity[0] * dt;
    body.centroid[1] += body.velocity[1] * dt;
    body.angle += body.angular_velocity * dt;
}
/// Oscillating body position for a sinusoidal motion.
///
/// x(t) = x0 + A * sin(2π f t),  y fixed.
pub fn oscillating_body_position(x0: f64, amplitude: f64, freq: f64, t: f64) -> f64 {
    x0 + amplitude * (2.0 * PI * freq * t).sin()
}
/// Oscillating body velocity (derivative of `oscillating_body_position`).
///
/// dx/dt = A * 2π f * cos(2π f t)
pub fn oscillating_body_velocity(amplitude: f64, freq: f64, t: f64) -> f64 {
    amplitude * 2.0 * PI * freq * (2.0 * PI * freq * t).cos()
}
/// Compute force on immersed body by integrating surface forces.
///
/// F_body = -Σ_k F_k * dA_k  (reaction to IBM forcing)
pub fn body_force_from_surface(surface_forces: &[[f64; 2]], ds: f64) -> [f64; 2] {
    let mut fx = 0.0_f64;
    let mut fy = 0.0_f64;
    for &[sfx, sfy] in surface_forces {
        fx += sfx * ds;
        fy += sfy * ds;
    }
    [-fx, -fy]
}
/// Compute torque on immersed body from surface forces.
///
/// τ = Σ_k (r_k × F_k) * dA_k  (2D cross product z-component)
pub fn body_torque_from_surface(
    surface_positions: &[[f64; 2]],
    surface_forces: &[[f64; 2]],
    centroid: [f64; 2],
    ds: f64,
) -> f64 {
    -surface_positions
        .iter()
        .zip(surface_forces.iter())
        .map(|(&[px, py], &[fx, fy])| {
            let rx = px - centroid[0];
            let ry = py - centroid[1];
            (rx * fy - ry * fx) * ds
        })
        .sum::<f64>()
}
/// Compute lift and drag on body from surface pressure.
///
/// Uses the momentum exchange method.
/// `cos_alpha`, `sin_alpha` define the flow direction.
pub fn body_lift_drag_from_pressure(
    surface_pressures: &[f64],
    surface_normals: &[[f64; 2]],
    surface_areas: &[f64],
    cos_alpha: f64,
    sin_alpha: f64,
) -> (f64, f64) {
    let n = surface_pressures.len();
    assert_eq!(surface_normals.len(), n);
    assert_eq!(surface_areas.len(), n);
    let mut drag = 0.0_f64;
    let mut lift = 0.0_f64;
    for i in 0..n {
        let p = surface_pressures[i];
        let nx = surface_normals[i][0];
        let ny = surface_normals[i][1];
        let da = surface_areas[i];
        let dfx = -p * nx * da;
        let dfy = -p * ny * da;
        drag += dfx * cos_alpha + dfy * sin_alpha;
        lift += -dfx * sin_alpha + dfy * cos_alpha;
    }
    (drag, lift)
}
#[cfg(test)]
mod tests_ibm_body {
    use super::super::functions::{peskin_delta_2d, peskin_delta_3d, yang_delta_1d};
    use super::*;
    use crate::immersed_boundary::types::*;
    #[test]
    fn test_ibm_body_new() {
        let b = IbmBody::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(b.centroid, [1.0, 2.0]);
        assert_eq!(b.velocity, [0.0, 0.0]);
        assert_eq!(b.mass, 3.0);
        assert_eq!(b.moment_of_inertia, 4.0);
    }
    #[test]
    fn test_ibm_body_add_surface_point() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.add_surface_point(1.0, 0.0);
        b.add_surface_point(0.0, 1.0);
        assert_eq!(b.surface_offsets.len(), 2);
        assert_eq!(b.surface_forces.len(), 2);
    }
    #[test]
    fn test_ibm_body_circle_surface_on_circle() {
        let r = 1.5;
        let b = IbmBody::circle(0.0, 0.0, r, 12, 1.0);
        assert_eq!(b.surface_offsets.len(), 12);
        for &[ox, oy] in &b.surface_offsets {
            let dist = (ox * ox + oy * oy).sqrt();
            assert!((dist - r).abs() < 1e-12, "surface point dist = {dist}");
        }
    }
    #[test]
    fn test_ibm_body_surface_positions_no_rotation() {
        let mut b = IbmBody::new(2.0, 3.0, 1.0, 1.0);
        b.add_surface_point(1.0, 0.0);
        let pos = b.surface_positions();
        assert!((pos[0][0] - 3.0).abs() < 1e-12);
        assert!((pos[0][1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_ibm_body_surface_velocities_no_rotation() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.velocity = [0.5, 0.3];
        b.add_surface_point(1.0, 0.0);
        let vels = b.surface_velocities();
        assert!((vels[0][0] - 0.5).abs() < 1e-12);
        assert!((vels[0][1] - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_ibm_body_net_force_zero() {
        let b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        let f = b.net_force();
        assert_eq!(f, [0.0, 0.0]);
    }
    #[test]
    fn test_ibm_body_net_force_sum() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.add_surface_point(1.0, 0.0);
        b.add_surface_point(-1.0, 0.0);
        b.surface_forces[0] = [2.0, 1.0];
        b.surface_forces[1] = [3.0, -1.0];
        let f = b.net_force();
        assert!((f[0] - 5.0).abs() < 1e-14);
        assert!(f[1].abs() < 1e-14);
    }
    #[test]
    fn test_ibm_body_net_torque_symmetric() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.add_surface_point(1.0, 0.0);
        b.add_surface_point(-1.0, 0.0);
        b.surface_forces[0] = [1.0, 0.0];
        b.surface_forces[1] = [-1.0, 0.0];
        let tau = b.net_torque();
        assert!(tau.abs() < 1e-12, "torque = {tau}");
    }
    #[test]
    fn test_ibm_body_advance_linear() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.add_surface_point(0.0, 1.0);
        b.surface_forces[0] = [1.0, 0.0];
        b.advance(0.1);
        assert!(
            b.velocity[0] > 0.0,
            "velocity should increase: {}",
            b.velocity[0]
        );
        assert!(
            b.centroid[0] > 0.0,
            "centroid should move: {}",
            b.centroid[0]
        );
    }
    #[test]
    fn test_peskin_delta_2d_at_origin() {
        let d = peskin_delta_2d(0.0, 0.0);
        assert!((d - 0.25).abs() < 1e-12, "peskin_delta_2d(0,0) = {d}");
    }
    #[test]
    fn test_peskin_delta_2d_outside_support() {
        let d = peskin_delta_2d(3.0, 0.0);
        assert_eq!(d, 0.0, "outside support");
    }
    #[test]
    fn test_peskin_delta_3d_at_origin() {
        let d = peskin_delta_3d(0.0, 0.0, 0.0);
        assert!((d - 0.125).abs() < 1e-12, "peskin_delta_3d(0,0,0) = {d}");
    }
    #[test]
    fn test_yang_delta_at_zero() {
        let d = yang_delta_1d(0.0);
        assert!((d - 2.0 / 3.0).abs() < 1e-12, "yang_delta(0) = {d}");
    }
    #[test]
    fn test_yang_delta_symmetric() {
        for &r in &[0.3, 0.8, 1.2] {
            let pos = yang_delta_1d(r);
            let neg = yang_delta_1d(-r);
            assert!(
                (pos - neg).abs() < 1e-14,
                "yang delta not symmetric at r={r}"
            );
        }
    }
    #[test]
    fn test_yang_delta_outside_support() {
        assert_eq!(yang_delta_1d(2.0), 0.0);
        assert_eq!(yang_delta_1d(-2.0), 0.0);
    }
    #[test]
    fn test_uhlmann_direct_forcing_zero_error() {
        let ub = vec![[0.5, 0.3]];
        let uf = vec![[0.5, 0.3]];
        let f = uhlmann_direct_forcing(&ub, &uf, 0.01);
        assert!(f[0][0].abs() < 1e-14, "zero error → zero force");
        assert!(f[0][1].abs() < 1e-14);
    }
    #[test]
    fn test_uhlmann_direct_forcing_value() {
        let ub = vec![[1.0, 0.0]];
        let uf = vec![[0.8, 0.0]];
        let f = uhlmann_direct_forcing(&ub, &uf, 0.1);
        assert!((f[0][0] - 2.0).abs() < 1e-12, "force = {}", f[0][0]);
    }
    #[test]
    fn test_uhlmann_ibm_step_nonzero() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let dt = 0.01;
        let fluid_u = vec![[0.0_f64; 2]; nx * ny];
        let body_pos = vec![[4.0, 4.0]];
        let body_vel = vec![[1.0, 0.0]];
        let f_grid = uhlmann_ibm_step(&body_pos, &body_vel, &fluid_u, nx, ny, dx, dt);
        let total_fx: f64 = f_grid.iter().map(|f| f[0]).sum();
        assert!(
            total_fx.abs() > 0.0,
            "Uhlmann IBM should produce nonzero force"
        );
    }
    #[test]
    fn test_interpolate_velocity_at_positions_uniform() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0;
        let fluid_u = vec![[0.3_f64, 0.1_f64]; nx * ny];
        let positions = vec![[3.5, 3.5]];
        let interp = interpolate_velocity_at_positions(&positions, &fluid_u, nx, ny, dx);
        assert!((interp[0][0] - 0.3).abs() < 0.02, "ux = {}", interp[0][0]);
        assert!((interp[0][1] - 0.1).abs() < 0.02, "uy = {}", interp[0][1]);
    }
    #[test]
    fn test_goldstein_forcing_proportional() {
        let mut gf = GoldsteinForcing::new(1000.0, 0.0, 1);
        let u_des = [[1.0, 0.0]];
        let u_int = [[0.9, 0.0]];
        let f = gf.compute_forces(&u_des, &u_int, 0.01);
        assert!((f[0][0] - 100.0).abs() < 1e-10, "F_x = {}", f[0][0]);
    }
    #[test]
    fn test_goldstein_forcing_zero_error() {
        let mut gf = GoldsteinForcing::new(1000.0, 100.0, 1);
        let u = [[0.5_f64, 0.2_f64]];
        let f = gf.compute_forces(&u, &u, 0.01);
        assert!(f[0][0].abs() < 1e-14, "zero error → zero force");
    }
    #[test]
    fn test_goldstein_forcing_reset() {
        let mut gf = GoldsteinForcing::new(1000.0, 100.0, 2);
        gf.error_integral[0] = [3.0, -2.0];
        gf.error_integral[1] = [1.0, 5.0];
        gf.reset();
        for s in &gf.error_integral {
            assert_eq!(*s, [0.0, 0.0], "integral should be zero after reset");
        }
    }
    #[test]
    fn test_advance_body_kinematics() {
        let mut b = IbmBody::new(0.0, 0.0, 1.0, 1.0);
        b.velocity = [1.0, 0.5];
        b.angular_velocity = 0.1;
        advance_body_kinematics(&mut b, 0.1);
        assert!((b.centroid[0] - 0.1).abs() < 1e-12);
        assert!((b.centroid[1] - 0.05).abs() < 1e-12);
        assert!((b.angle - 0.01).abs() < 1e-12);
    }
    #[test]
    fn test_oscillating_body_position_at_zero() {
        let x = oscillating_body_position(1.0, 0.5, 1.0, 0.0);
        assert!((x - 1.0).abs() < 1e-12, "at t=0: x = {x}");
    }
    #[test]
    fn test_oscillating_body_position_amplitude() {
        let x_max = oscillating_body_position(0.0, 0.5, 1.0, 0.25);
        assert!((x_max - 0.5).abs() < 1e-12, "max position = {x_max}");
    }
    #[test]
    fn test_oscillating_body_velocity_at_zero() {
        let v = oscillating_body_velocity(0.5, 1.0, 0.0);
        let expected = 0.5 * 2.0 * PI * 1.0;
        assert!((v - expected).abs() < 1e-10, "v at t=0 = {v}");
    }
    #[test]
    fn test_body_force_from_surface_reaction() {
        let surface_forces = vec![[1.0_f64, 0.0_f64], [1.0, 0.0]];
        let f = body_force_from_surface(&surface_forces, 1.0);
        assert!((f[0] - (-2.0)).abs() < 1e-14, "body fx = {}", f[0]);
        assert!(f[1].abs() < 1e-14);
    }
    #[test]
    fn test_body_torque_from_surface_zero() {
        let pos = vec![[1.0_f64, 0.0_f64], [-1.0, 0.0]];
        let forces = vec![[1.0_f64, 0.0_f64], [-1.0, 0.0]];
        let tau = body_torque_from_surface(&pos, &forces, [0.0, 0.0], 1.0);
        assert!(tau.abs() < 1e-12, "torque = {tau}");
    }
    #[test]
    fn test_body_torque_nonzero() {
        let pos = vec![[1.0_f64, 0.0_f64]];
        let forces = vec![[0.0_f64, 1.0_f64]];
        let tau = body_torque_from_surface(&pos, &forces, [0.0, 0.0], 1.0);
        assert!((tau - (-1.0)).abs() < 1e-12, "torque = {tau}");
    }
    #[test]
    fn test_body_lift_drag_from_pressure_basic() {
        let pressures = [1.0_f64];
        let normals = [[1.0_f64, 0.0_f64]];
        let areas = [1.0_f64];
        let (drag, lift) = body_lift_drag_from_pressure(&pressures, &normals, &areas, 1.0, 0.0);
        assert!(lift.abs() < 1e-12, "lift = {lift}");
        assert!(drag.abs() > 0.0, "drag should be nonzero: {drag}");
    }
    #[test]
    fn test_spread_force_uhlmann_nonzero() {
        let positions = vec![[4.0_f64, 4.0_f64]];
        let forces = vec![[1.0_f64, 0.0_f64]];
        let f = spread_force_uhlmann(&positions, &forces, 8, 8, 1.0);
        let total_fx: f64 = f.iter().map(|v| v[0]).sum();
        assert!(
            total_fx > 0.0,
            "Uhlmann spread force should be nonzero: {total_fx}"
        );
    }
    #[test]
    fn test_spread_force_uhlmann_zero_force() {
        let positions = vec![[4.0_f64, 4.0_f64]];
        let forces = vec![[0.0_f64, 0.0_f64]];
        let f = spread_force_uhlmann(&positions, &forces, 8, 8, 1.0);
        let total: f64 = f.iter().map(|v| v[0] + v[1]).sum();
        assert!(total.abs() < 1e-15, "zero force → zero grid force: {total}");
    }
}
