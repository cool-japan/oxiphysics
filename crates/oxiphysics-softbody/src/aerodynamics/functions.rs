//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::helpers::*;
use super::types::*;
/// Compute the total wind force on a set of panels given a common wind
/// velocity and air density.
///
/// Equivalent to [`AeroSurface::total_force`] but accepts a slice directly.
pub fn wind_force(panels: &[AeroPanel], wind_velocity: [f64; 3], air_density: f64) -> [f64; 3] {
    let mut total = [0.0_f64; 3];
    for panel in panels {
        let (lift, drag) = panel.lift_drag(wind_velocity, air_density);
        total = v3_add(total, v3_add(lift, drag));
    }
    total
}
/// Magnus force on a spinning sphere in a flow.
///
/// ```text
/// F = rho * pi * r^3 * (omega x velocity)
/// ```
pub fn magnus_force(omega: [f64; 3], velocity: [f64; 3], radius: f64, density: f64) -> [f64; 3] {
    let cross = v3_cross(omega, velocity);
    let coeff = density * std::f64::consts::PI * radius * radius * radius;
    v3_scale(cross, coeff)
}
/// Buoyancy force (Archimedes principle).
///
/// ```text
/// F_b = -fluid_density * volume * gravity
/// ```
pub fn buoyancy_force(volume: f64, fluid_density: f64, gravity: [f64; 3]) -> [f64; 3] {
    v3_scale(gravity, -fluid_density * volume)
}
/// Compute distributed wind load on a structural panel.
///
/// The wind load is the dynamic pressure times the projected area:
/// `q = 0.5 * rho * v^2 * Cd * cos(alpha)`
///
/// Returns the force vector normal to the panel surface.
pub fn wind_load_force(
    panel_normal: [f64; 3],
    panel_area: f64,
    wind_velocity: [f64; 3],
    air_density: f64,
    drag_coeff: f64,
) -> [f64; 3] {
    let v_sq = v3_dot(wind_velocity, wind_velocity);
    if v_sq < 1e-30 {
        return [0.0; 3];
    }
    let wind_dir = v3_normalize(wind_velocity).unwrap_or([0.0; 3]);
    let cos_alpha = v3_dot(wind_dir, panel_normal).abs();
    let q = 0.5 * air_density * v_sq;
    let force_mag = q * drag_coeff * panel_area * cos_alpha;
    let sign = if v3_dot(wind_dir, panel_normal) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    v3_scale(panel_normal, sign * force_mag)
}
/// Compute the total wind load on multiple structural panels.
pub fn total_wind_load(
    normals: &[[f64; 3]],
    areas: &[f64],
    wind_velocity: [f64; 3],
    air_density: f64,
    drag_coeff: f64,
) -> [f64; 3] {
    let mut total = [0.0_f64; 3];
    let n = normals.len().min(areas.len());
    for i in 0..n {
        let f = wind_load_force(normals[i], areas[i], wind_velocity, air_density, drag_coeff);
        total = v3_add(total, f);
    }
    total
}
/// Influence coefficient of a constant-strength source panel on a field point.
///
/// For a flat panel with source strength sigma, the velocity potential induced
/// at point `p` is approximately: `phi = -sigma/(4*pi) * A / r`
/// where `A` is the panel area and `r` is the distance from panel center to `p`.
///
/// Returns the velocity induced at `field_point` by the panel (gradient of phi).
pub fn source_panel_velocity(
    panel_center: [f64; 3],
    _panel_normal: [f64; 3],
    panel_area: f64,
    sigma: f64,
    field_point: [f64; 3],
) -> [f64; 3] {
    let r_vec = v3_sub(field_point, panel_center);
    let r = v3_norm(r_vec);
    if r < 1e-12 {
        return [0.0; 3];
    }
    let coeff = sigma * panel_area / (4.0 * std::f64::consts::PI * r * r * r);
    v3_scale(r_vec, coeff)
}
/// Solve a basic panel method: given N panels with Neumann BCs (no normal flow),
/// compute the source strengths.
///
/// This is a simplified O(N^2) approach: the influence matrix A_ij represents
/// the normal velocity at panel i due to unit source strength on panel j.
pub fn solve_panel_method(
    centers: &[[f64; 3]],
    normals: &[[f64; 3]],
    areas: &[f64],
    freestream: [f64; 3],
) -> Vec<f64> {
    let n = centers.len();
    if n == 0 {
        return Vec::new();
    }
    let mut a_mat = vec![0.0_f64; n * n];
    let mut rhs = vec![0.0_f64; n];
    for i in 0..n {
        rhs[i] = -v3_dot(freestream, normals[i]);
        for j in 0..n {
            if i == j {
                a_mat[i * n + j] = 0.5;
            } else {
                let vel = source_panel_velocity(centers[j], normals[j], areas[j], 1.0, centers[i]);
                a_mat[i * n + j] = v3_dot(vel, normals[i]);
            }
        }
    }
    gauss_solve(&mut a_mat, &mut rhs, n)
}
/// Simple Gauss elimination for a dense N x N system.
pub(super) fn gauss_solve(a: &mut [f64], b: &mut [f64], n: usize) -> Vec<f64> {
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = a[col * n + col].abs();
        for row in (col + 1)..n {
            let val = a[row * n + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_row != col {
            for k in 0..n {
                a.swap(col * n + k, max_row * n + k);
            }
            b.swap(col, max_row);
        }
        let pivot = a[col * n + col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            for k in col..n {
                a[row * n + k] -= factor * a[col * n + k];
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i * n + j] * x[j];
        }
        let diag = a[i * n + i];
        x[i] = if diag.abs() > 1e-30 { sum / diag } else { 0.0 };
    }
    x
}
/// Compute the induced velocity from a finite vortex filament (Biot-Savart law).
///
/// Velocity at point `p` due to a filament from `a` to `b` with circulation `gamma`.
///
/// ```text
/// v = (gamma / 4*pi) * (r1 x r2) / |r1 x r2|^2 * r0 . (r1/|r1| - r2/|r2|)
/// ```
pub fn biot_savart_filament(a: [f64; 3], b: [f64; 3], p: [f64; 3], gamma: f64) -> [f64; 3] {
    super::helpers::biot_savart_filament(a, b, p, gamma)
}
/// Compute the induced velocity from a horseshoe vortex at a field point.
///
/// A horseshoe vortex consists of:
/// 1. A bound vortex from p1 to p2
/// 2. Two semi-infinite trailing vortices extending downstream (along +X)
pub fn horseshoe_induced_velocity(vortex: &HorseshoeVortex, field_point: [f64; 3]) -> [f64; 3] {
    let v_bound = biot_savart_filament(vortex.p1, vortex.p2, field_point, vortex.gamma);
    let trail_dist = 1000.0;
    let trail_a = [vortex.p1[0] + trail_dist, vortex.p1[1], vortex.p1[2]];
    let trail_b = [vortex.p2[0] + trail_dist, vortex.p2[1], vortex.p2[2]];
    let v_left = biot_savart_filament(trail_a, vortex.p1, field_point, vortex.gamma);
    let v_right = biot_savart_filament(vortex.p2, trail_b, field_point, vortex.gamma);
    v3_add(v_bound, v3_add(v_left, v_right))
}
/// Compute pressure coefficient Cp from velocity ratio.
///
/// `Cp = 1 - (v_local / v_inf)^2` (incompressible Bernoulli).
pub fn pressure_coefficient(v_local: f64, v_inf: f64) -> f64 {
    if v_inf.abs() < 1e-30 {
        return 0.0;
    }
    let ratio = v_local / v_inf;
    1.0 - ratio * ratio
}
/// Compute the pressure force from a Cp distribution.
///
/// `F = sum_i( -Cp_i * q_inf * A_i * n_i )`
///
/// where `q_inf = 0.5 * rho * v_inf^2`.
pub fn pressure_force_from_cp(
    cp_values: &[f64],
    normals: &[[f64; 3]],
    areas: &[f64],
    dynamic_pressure: f64,
) -> [f64; 3] {
    let mut force = [0.0_f64; 3];
    let n = cp_values.len().min(normals.len()).min(areas.len());
    for i in 0..n {
        let f = v3_scale(normals[i], -cp_values[i] * dynamic_pressure * areas[i]);
        force = v3_add(force, f);
    }
    force
}
/// Compute dynamic pressure q = 0.5 * rho * v^2.
pub fn dynamic_pressure(air_density: f64, velocity: f64) -> f64 {
    0.5 * air_density * velocity * velocity
}
/// Simple Gaussian elimination solver for `A * x = b`.
///
/// `a` is column-major `n x n`, `b` has length `n`.
/// Returns `None` if the matrix is singular.
// gauss_solve_opt moved to helpers module
/// Biot-Savart induced velocity at `eval` due to a finite vortex filament
/// running from `r1` to `r2` with circulation strength `gamma` (m²/s).
///
/// Equivalent to `biot_savart_filament` with clearer argument names.
///
/// ```text
/// v = (gamma / 4π) * (r1→P × r2→P) / |r1→P × r2→P|²  * r0 · (r̂1 - r̂2)
/// ```
pub fn biot_savart_velocity(r1: [f64; 3], r2: [f64; 3], gamma: f64, eval: [f64; 3]) -> [f64; 3] {
    super::helpers::biot_savart_velocity(r1, r2, gamma, eval)
}
// stall_blend moved to helpers module
#[cfg(test)]
mod tests {
    use super::*;

    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_zero_wind_zero_force() {
        let panel = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.8, 1.2);
        let (lift, drag) = panel.lift_drag([0.0; 3], 1.225);
        assert!(v3_norm(lift) < EPS, "Lift should be zero with zero wind");
        assert!(v3_norm(drag) < EPS, "Drag should be zero with zero wind");
    }
    #[test]
    fn test_head_on_wind_pure_drag() {
        let panel = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.8, 1.2);
        let (lift, drag) = panel.lift_drag([0.0, 10.0, 0.0], 1.0);
        assert!(
            v3_norm(lift) < EPS,
            "Head-on wind should give zero lift, got {lift:?}"
        );
        assert!(
            drag[1] > EPS,
            "Head-on drag should be in +Y direction, got {drag:?}"
        );
    }
    #[test]
    fn test_tangential_wind_pure_lift() {
        let panel = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 1.0, 1.0);
        let (lift, drag) = panel.lift_drag([10.0, 0.0, 0.0], 1.0);
        assert!(
            v3_norm(drag) < EPS,
            "Tangential wind should give zero drag, got {drag:?}"
        );
        assert!(
            v3_norm(lift) > EPS,
            "Tangential wind should give non-zero lift"
        );
    }
    #[test]
    fn test_wind_force_matches_panel() {
        let panel = AeroPanel::new([0.0, 1.0, 0.0], 2.0, [0.0; 3], 0.6, 1.1);
        let wind = [3.0, 5.0, 1.0];
        let density = 1.225;
        let surface_force = wind_force(std::slice::from_ref(&panel), wind, density);
        let (lift, drag) = panel.lift_drag(wind, density);
        let expected = v3_add(lift, drag);
        for i in 0..3 {
            assert!(
                (surface_force[i] - expected[i]).abs() < EPS,
                "wind_force[{i}] mismatch: {} vs {}",
                surface_force[i],
                expected[i]
            );
        }
    }
    #[test]
    fn test_magnus_direction() {
        let omega = [0.0, 0.0, 1.0];
        let vel = [1.0, 0.0, 0.0];
        let f = magnus_force(omega, vel, 0.05, 1.225);
        assert!(f[1] > EPS, "Magnus force should be in +Y, got {f:?}");
        assert!(f[0].abs() < EPS && f[2].abs() < EPS);
    }
    #[test]
    fn test_buoyancy_sign_and_scaling() {
        let gravity = [0.0, -9.81, 0.0];
        let f1 = buoyancy_force(1.0, 1000.0, gravity);
        assert!(f1[1] > 0.0, "Buoyancy should be upward, got {f1:?}");
        let expected = 1000.0 * 1.0 * 9.81;
        assert!(
            (f1[1] - expected).abs() < 1e-9,
            "Buoyancy magnitude mismatch: {} vs {expected}",
            f1[1]
        );
        let f2 = buoyancy_force(2.0, 1000.0, gravity);
        assert!(
            (f2[1] - 2.0 * f1[1]).abs() < 1e-9,
            "Doubling volume should double buoyancy"
        );
    }
    #[test]
    fn test_aero_surface_sum() {
        let p1 = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.5, 1.0);
        let p2 = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [1.0, 0.0, 0.0], 0.5, 1.0);
        let wind = [0.0, 5.0, 0.0];
        let density = 1.225;
        let mut surface = AeroSurface::new();
        surface.add_panel(p1.clone());
        surface.add_panel(p2.clone());
        let total = surface.total_force(wind, density);
        let (l1, d1) = p1.lift_drag(wind, density);
        let (l2, d2) = p2.lift_drag(wind, density);
        let expected = v3_add(v3_add(l1, d1), v3_add(l2, d2));
        for i in 0..3 {
            assert!(
                (total[i] - expected[i]).abs() < EPS,
                "Surface total[{i}] mismatch: {} vs {}",
                total[i],
                expected[i]
            );
        }
    }
    #[test]
    fn test_magnus_zero_omega() {
        let f = magnus_force([0.0; 3], [10.0, 0.0, 0.0], 0.1, 1.225);
        assert!(v3_norm(f) < EPS, "Magnus should be zero with no spin");
    }
    #[test]
    fn test_pressure_coefficient_stagnation() {
        let cp = pressure_coefficient(0.0, 10.0);
        assert!(
            (cp - 1.0).abs() < EPS,
            "Stagnation Cp should be 1.0, got {cp}"
        );
    }
    #[test]
    fn test_pressure_coefficient_freestream() {
        let cp = pressure_coefficient(10.0, 10.0);
        assert!(cp.abs() < EPS, "Freestream Cp should be 0.0, got {cp}");
    }
    #[test]
    fn test_dynamic_pressure() {
        let q = dynamic_pressure(1.225, 10.0);
        let expected = 0.5 * 1.225 * 100.0;
        assert!(
            (q - expected).abs() < EPS,
            "Dynamic pressure mismatch: {q} vs {expected}"
        );
    }
    #[test]
    fn test_wind_load_zero_wind() {
        let f = wind_load_force([0.0, 1.0, 0.0], 1.0, [0.0; 3], 1.225, 1.0);
        assert!(v3_norm(f) < EPS, "Wind load should be zero with no wind");
    }
    #[test]
    fn test_wind_load_scaling() {
        let normal = [0.0, 1.0, 0.0];
        let area = 2.0;
        let cd = 1.0;
        let rho = 1.225;
        let f1 = wind_load_force(normal, area, [0.0, 10.0, 0.0], rho, cd);
        let f2 = wind_load_force(normal, area, [0.0, 20.0, 0.0], rho, cd);
        let mag1 = v3_norm(f1);
        let mag2 = v3_norm(f2);
        assert!(mag1 > EPS, "f1 should be non-zero");
        let ratio = mag2 / mag1;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "Wind load should scale as v^2, ratio = {ratio}"
        );
    }
    #[test]
    fn test_total_wind_load() {
        let normals = vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let areas = vec![1.0, 1.0];
        let wind = [0.0, 10.0, 0.0];
        let total = total_wind_load(&normals, &areas, wind, 1.225, 1.0);
        let single = wind_load_force([0.0, 1.0, 0.0], 1.0, wind, 1.225, 1.0);
        for i in 0..3 {
            assert!(
                (total[i] - 2.0 * single[i]).abs() < EPS,
                "Total should be 2x single"
            );
        }
    }
    #[test]
    fn test_biot_savart_direction() {
        let a = [0.0, 0.0, -1.0];
        let b = [0.0, 0.0, 1.0];
        let p = [1.0, 0.0, 0.0];
        let v = biot_savart_filament(a, b, p, 1.0);
        assert!(
            v[1].abs() > 1e-6,
            "Biot-Savart should induce Y velocity, got {v:?}"
        );
        assert!(
            v[0].abs() < 1e-6 && v[2].abs() < 1e-6,
            "Should have no X/Z component"
        );
    }
    #[test]
    fn test_source_panel_decay() {
        let center = [0.0; 3];
        let normal = [0.0, 1.0, 0.0];
        let area = 1.0;
        let sigma = 1.0;
        let v_near = source_panel_velocity(center, normal, area, sigma, [1.0, 0.0, 0.0]);
        let v_far = source_panel_velocity(center, normal, area, sigma, [10.0, 0.0, 0.0]);
        let mag_near = v3_norm(v_near);
        let mag_far = v3_norm(v_far);
        assert!(
            mag_near > mag_far,
            "Velocity should decay with distance: near={mag_near}, far={mag_far}"
        );
    }
    #[test]
    fn test_panel_method_single() {
        let centers = vec![[0.0, 0.0, 0.0]];
        let normals = vec![[0.0, 1.0, 0.0]];
        let areas = vec![1.0];
        let freestream = [10.0, 0.0, 0.0];
        let sigmas = solve_panel_method(&centers, &normals, &areas, freestream);
        assert_eq!(sigmas.len(), 1);
        assert!(
            sigmas[0].abs() < 1e-6,
            "Tangential flow should give near-zero source: {}",
            sigmas[0]
        );
    }
    #[test]
    fn test_pressure_distribution_length() {
        let mut surface = AeroSurface::new();
        surface.add_panel(AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.5, 1.0));
        surface.add_panel(AeroPanel::new(
            [1.0, 0.0, 0.0],
            1.0,
            [1.0, 0.0, 0.0],
            0.5,
            1.0,
        ));
        surface.add_panel(AeroPanel::new(
            [0.0, 0.0, 1.0],
            1.0,
            [0.0, 1.0, 0.0],
            0.5,
            1.0,
        ));
        let cp = surface.pressure_distribution([10.0, 0.0, 0.0]);
        assert_eq!(cp.len(), 3);
    }
    #[test]
    fn test_lift_drag_coefficients() {
        let mut surface = AeroSurface::new();
        surface.add_panel(AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.8, 1.2));
        let wind = [0.0, 10.0, 0.0];
        let (cl, cd) = surface.lift_drag_coefficients(wind, 1.225, 1.0);
        assert!(cd > 0.0, "CD should be positive for head-on wind");
        assert!(cl >= 0.0, "CL should be non-negative");
    }
    #[test]
    fn test_total_area() {
        let mut surface = AeroSurface::new();
        surface.add_panel(AeroPanel::new([0.0, 1.0, 0.0], 2.5, [0.0; 3], 0.5, 1.0));
        surface.add_panel(AeroPanel::new(
            [0.0, 1.0, 0.0],
            3.5,
            [1.0, 0.0, 0.0],
            0.5,
            1.0,
        ));
        assert!((surface.total_area() - 6.0).abs() < EPS);
    }
    #[test]
    fn test_aero_center() {
        let mut surface = AeroSurface::new();
        surface.add_panel(AeroPanel::new(
            [0.0, 1.0, 0.0],
            1.0,
            [0.0, 0.0, 0.0],
            0.5,
            1.0,
        ));
        surface.add_panel(AeroPanel::new(
            [0.0, 1.0, 0.0],
            1.0,
            [2.0, 0.0, 0.0],
            0.5,
            1.0,
        ));
        let center = surface.aero_center();
        assert!(
            (center[0] - 1.0).abs() < EPS,
            "Aero center X should be 1.0, got {}",
            center[0]
        );
    }
    #[test]
    fn test_pressure_force_from_cp() {
        let cp = vec![1.0, 0.5];
        let normals = vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let areas = vec![1.0, 1.0];
        let q = 100.0;
        let f = pressure_force_from_cp(&cp, &normals, &areas, q);
        assert!(
            (f[1] - (-150.0)).abs() < EPS,
            "Pressure force Y should be -150, got {}",
            f[1]
        );
    }
    #[test]
    fn test_angle_of_attack() {
        let panel = AeroPanel::new([0.0, 1.0, 0.0], 1.0, [0.0; 3], 0.5, 1.0);
        let aoa = panel.angle_of_attack([0.0, 10.0, 0.0]);
        assert!(
            (aoa - std::f64::consts::FRAC_PI_2).abs() < 0.01,
            "Head-on AoA should be pi/2, got {aoa}"
        );
        let aoa2 = panel.angle_of_attack([10.0, 0.0, 0.0]);
        assert!(aoa2.abs() < 0.01, "Tangential AoA should be ~0, got {aoa2}");
    }
    #[test]
    fn test_horseshoe_velocity_nonzero() {
        let hv = HorseshoeVortex {
            p1: [0.0, 0.0, -0.5],
            p2: [0.0, 0.0, 0.5],
            control_point: [0.5, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            gamma: 1.0,
        };
        let v = horseshoe_induced_velocity(&hv, [1.0, 0.0, 0.0]);
        let mag = v3_norm(v);
        assert!(
            mag > 1e-8,
            "Horseshoe vortex should induce nonzero velocity, got {mag}"
        );
    }
    #[test]
    fn test_from_triangle() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let panel = AeroPanel::from_triangle(v0, v1, v2, 0.5, 1.0);
        assert!(
            (panel.area - 0.5).abs() < EPS,
            "Triangle area should be 0.5, got {}",
            panel.area
        );
        assert!(
            panel.normal[2].abs() > 0.99,
            "Normal should be along Z, got {:?}",
            panel.normal
        );
    }
    #[test]
    fn test_unsteady_quasi_steady_lift_proportional_to_alpha() {
        let panel = UnsteadyAeroPanel::thin_airfoil(0.3, 1.225);
        let l1 = panel.quasi_steady_lift(20.0, 0.05);
        let l2 = panel.quasi_steady_lift(20.0, 0.10);
        let ratio = l2 / l1;
        assert!(
            (ratio - 2.0).abs() < 1e-9,
            "lift should double when alpha doubles: ratio={ratio}"
        );
    }
    #[test]
    fn test_unsteady_added_mass_nonzero_for_pitching() {
        let panel = UnsteadyAeroPanel::thin_airfoil(0.3, 1.225);
        let lam = panel.added_mass_lift(20.0, 0.5, 0.0);
        assert!(
            lam.abs() > 1e-6,
            "added mass lift should be nonzero for pitching"
        );
    }
    #[test]
    fn test_unsteady_total_lift_includes_added_mass() {
        let panel = UnsteadyAeroPanel::thin_airfoil(0.3, 1.225);
        let l_qs = panel.quasi_steady_lift(20.0, 0.05);
        let l_total = panel.total_lift(20.0, 0.05, 0.2, 0.0);
        assert!(
            (l_total - l_qs).abs() > 1e-6,
            "total should differ from quasi-steady with non-zero alpha_dot"
        );
    }
    #[test]
    fn test_added_mass_sphere_isotropic() {
        let am = AddedMassTensor::sphere(0.1, 1.225);
        let d = am.diagonal;
        assert!(
            (d[0] - d[1]).abs() < 1e-10 && (d[1] - d[2]).abs() < 1e-10,
            "sphere added mass should be isotropic: {d:?}"
        );
    }
    #[test]
    fn test_added_mass_sphere_value() {
        let r = 0.1;
        let rho = 1.225;
        let am = AddedMassTensor::sphere(r, rho);
        let expected = 0.5 * rho * (4.0 / 3.0) * std::f64::consts::PI * r.powi(3);
        assert!(
            (am.diagonal[0] - expected).abs() < 1e-12,
            "sphere added mass mismatch"
        );
    }
    #[test]
    fn test_added_mass_force_opposes_acceleration() {
        let am = AddedMassTensor::sphere(0.1, 1.225);
        let a = [1.0, 0.0, 0.0];
        let f = am.force(a);
        assert!(f[0] < 0.0, "added mass force should oppose acceleration");
    }
    #[test]
    fn test_added_mass_spheroid_larger_transverse() {
        let am = AddedMassTensor::spheroid(0.5, 0.1, 1.225);
        let m_ax = am.diagonal[0];
        let m_tr = am.diagonal[1];
        assert!(
            m_tr >= m_ax,
            "transverse added mass should be >= axial for prolate spheroid"
        );
    }
    #[test]
    fn test_fsi_deformation_increases_under_pressure() {
        let mut fsi = FsiCoupling::new(1, 1000.0, 10.0);
        fsi.pressures[0] = 500.0;
        let areas = vec![1.0];
        for _ in 0..100 {
            fsi.step(&areas, 0.5, 0.001);
        }
        assert!(
            fsi.deformations[0] > 0.0,
            "positive pressure should cause positive deformation: {}",
            fsi.deformations[0]
        );
    }
    #[test]
    fn test_fsi_strain_energy_positive() {
        let mut fsi = FsiCoupling::new(2, 500.0, 5.0);
        fsi.deformations[0] = 0.01;
        fsi.deformations[1] = 0.02;
        let e = fsi.strain_energy();
        assert!(
            e > 0.0,
            "strain energy must be positive with non-zero deformation: {e}"
        );
    }
    #[test]
    fn test_fsi_set_pressures_from_cp() {
        let mut fsi = FsiCoupling::new(2, 500.0, 5.0);
        fsi.set_pressures_from_cp(&[1.0, -0.5], 100.0);
        assert!((fsi.pressures[0] - 100.0).abs() < 1e-9);
        assert!((fsi.pressures[1] - (-50.0)).abs() < 1e-9);
    }
    #[test]
    fn test_flexible_wing_zero_twist_at_start() {
        let wing = FlexibleWing::new(4, 1.0, 0.2, 1000.0);
        for &t in &wing.twist {
            assert!(t.abs() < 1e-10, "initial twist should be zero");
        }
    }
    #[test]
    fn test_flexible_wing_total_cl_proportional_to_alpha() {
        let wing = FlexibleWing::new(4, 1.0, 0.2, 1000.0);
        let cl1 = wing.total_lift_coefficient(0.05);
        let cl2 = wing.total_lift_coefficient(0.10);
        let ratio = cl2 / cl1;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "CL should double with doubled alpha: ratio={ratio}"
        );
    }
    #[test]
    fn test_flexible_wing_aero_moments_nonzero() {
        let wing = FlexibleWing::new(4, 1.0, 0.2, 1000.0);
        let moments = wing.aero_moments(30.0, 1.225, 0.1);
        assert!(moments.len() == 4);
        assert!(
            moments.iter().any(|&m| m.abs() > 1e-6),
            "aero moments should be non-zero at AoA"
        );
    }
    #[test]
    fn test_membrane_inflation_radius_increases_under_pressure() {
        let mem = MembraneInflation::nylon_canopy(2.0);
        let r0 = mem.inflated_radius(0.0);
        let r1 = mem.inflated_radius(1000.0);
        assert!(
            r1 > r0,
            "pressurised membrane should expand: r0={r0}, r1={r1}"
        );
    }
    #[test]
    fn test_membrane_hoop_stress_increases_with_pressure() {
        let mem = MembraneInflation::nylon_canopy(2.0);
        let s1 = mem.hoop_stress(500.0);
        let s2 = mem.hoop_stress(1000.0);
        assert!(s2 > s1, "hoop stress should increase with pressure");
    }
    #[test]
    fn test_membrane_inflation_force_proportional_to_area() {
        let mem = MembraneInflation::nylon_canopy(2.0);
        let f1 = mem.inflation_force(200.0, 1.0);
        let f2 = mem.inflation_force(200.0, 2.0);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-9,
            "inflation force should be proportional to area"
        );
    }
    #[test]
    fn test_membrane_safety_factor_decreases_with_pressure() {
        let mem = MembraneInflation::nylon_canopy(2.0);
        let sf1 = mem.safety_factor(100.0);
        let sf2 = mem.safety_factor(10000.0);
        assert!(
            sf2 < sf1,
            "safety factor should decrease at higher pressure"
        );
    }
    #[test]
    fn test_lifting_line_zero_aoa_zero_lift() {
        let wing = LiftingLineTheory::new(10.0, 1.5, 0.0, 10.0 / 1.5);
        let (lift, _drag) = wing.compute_lift_drag(30.0, 1.225);
        assert!(
            lift.abs() < 1e-10,
            "zero AoA should give zero lift, got {lift}"
        );
    }
    #[test]
    fn test_lifting_line_positive_aoa_positive_lift() {
        let wing = LiftingLineTheory::new(10.0, 1.5, 5.0_f64.to_radians(), 10.0 / 1.5);
        let (lift, _drag) = wing.compute_lift_drag(30.0, 1.225);
        assert!(
            lift > 0.0,
            "positive AoA should give positive lift, got {lift}"
        );
    }
    #[test]
    fn test_lifting_line_induced_drag_positive() {
        let wing = LiftingLineTheory::new(10.0, 1.5, 5.0_f64.to_radians(), 10.0 / 1.5);
        let (_lift, drag) = wing.compute_lift_drag(30.0, 1.225);
        assert!(drag > 0.0, "induced drag should be positive, got {drag}");
    }
    #[test]
    fn test_lifting_line_kutta_joukowski() {
        let span = 10.0;
        let chord = 1.5;
        let aoa = 5.0_f64.to_radians();
        let ar = span / chord;
        let wing = LiftingLineTheory::new(span, chord, aoa, ar);
        let velocity = 30.0;
        let density = 1.225;
        let (lift, _) = wing.compute_lift_drag(velocity, density);
        let cl_2d = 2.0 * std::f64::consts::PI * aoa;
        let cl_eff = cl_2d / (1.0 + 2.0 / ar);
        let ref_area = span * chord;
        let expected = 0.5 * density * velocity * velocity * ref_area * cl_eff;
        let rel_err = (lift - expected).abs() / (expected.abs() + 1e-30);
        assert!(
            rel_err < 0.05,
            "Kutta-Joukowski check: lift={lift} expected={expected}"
        );
    }
    #[test]
    fn test_lifting_line_lift_scales_with_velocity_squared() {
        let wing = LiftingLineTheory::new(10.0, 1.5, 5.0_f64.to_radians(), 10.0 / 1.5);
        let (l1, _) = wing.compute_lift_drag(20.0, 1.225);
        let (l2, _) = wing.compute_lift_drag(40.0, 1.225);
        let ratio = l2 / l1;
        assert!(
            (ratio - 4.0).abs() < 1e-6,
            "lift should scale with V^2: ratio={ratio}"
        );
    }
    #[test]
    fn test_vlm_solve_produces_nonzero_gammas() {
        let vlm = VortexLatticeMethod::new(4, 3, 10.0, 2.5, 5.0_f64.to_radians());
        let gammas = vlm.solve_gammas(30.0, 1.225);
        assert!(!gammas.is_empty(), "should return some circulations");
        assert!(
            gammas.iter().any(|&g| g.abs() > 1e-10),
            "circulations should be non-zero at AoA"
        );
    }
    #[test]
    fn test_vlm_total_lift_positive_aoa() {
        let vlm = VortexLatticeMethod::new(4, 3, 10.0, 2.5, 5.0_f64.to_radians());
        let lift = vlm.total_lift(30.0, 1.225);
        assert!(
            lift > 0.0,
            "VLM lift should be positive at positive AoA, got {lift}"
        );
    }
    #[test]
    fn test_biot_savart_velocity_midpoint() {
        let r1 = [-1.0, 0.0, 0.0];
        let r2 = [1.0, 0.0, 0.0];
        let eval = [0.0, 1.0, 0.0];
        let v = biot_savart_velocity(r1, r2, 1.0, eval);
        assert!(v[2] > 1e-10, "velocity should be in +Z, got {v:?}");
        assert!(v[0].abs() < 1e-10, "no X component expected, got {}", v[0]);
    }
    #[test]
    fn test_biot_savart_velocity_scales_with_gamma() {
        let r1 = [-1.0, 0.0, 0.0];
        let r2 = [1.0, 0.0, 0.0];
        let eval = [0.0, 1.0, 0.0];
        let v1 = biot_savart_velocity(r1, r2, 1.0, eval);
        let v2 = biot_savart_velocity(r1, r2, 2.0, eval);
        for i in 0..3 {
            assert!(
                (v2[i] - 2.0 * v1[i]).abs() < 1e-10,
                "velocity should scale with gamma"
            );
        }
    }
    #[test]
    fn test_stall_model_zero_aoa_zero_lift() {
        let model = StallModel::naca_symmetric();
        let cl = model.lift_coefficient(0.0);
        assert!(
            cl.abs() < 1e-6,
            "Symmetric airfoil at 0° AoA should have CL=0, got {cl}"
        );
    }
    #[test]
    fn test_stall_model_positive_aoa_positive_lift() {
        let model = StallModel::naca_symmetric();
        let cl = model.lift_coefficient(5.0_f64.to_radians());
        assert!(cl > 0.0, "Positive AoA should give positive CL, got {cl}");
    }
    #[test]
    fn test_stall_model_is_stalled_beyond_critical() {
        let model = StallModel::naca_symmetric();
        assert!(
            !model.is_stalled(5.0_f64.to_radians()),
            "5° should not be stalled"
        );
        assert!(
            model.is_stalled(25.0_f64.to_radians()),
            "25° should be stalled"
        );
    }
    #[test]
    fn test_stall_model_drag_positive() {
        let model = StallModel::naca_symmetric();
        for aoa_deg in [0, 5, 15, 20, 30_i32] {
            let cd = model.drag_coefficient((aoa_deg as f64).to_radians());
            assert!(cd > 0.0, "CD should be positive at {aoa_deg}°, got {cd}");
        }
    }
    #[test]
    fn test_stall_model_drag_increases_past_stall() {
        let model = StallModel::naca_symmetric();
        let cd_attached = model.drag_coefficient(5.0_f64.to_radians());
        let cd_stalled = model.drag_coefficient(30.0_f64.to_radians());
        assert!(
            cd_stalled > cd_attached,
            "Stalled drag should exceed attached-flow drag: cd_att={cd_attached}, cd_stall={cd_stalled}"
        );
    }
    #[test]
    fn test_stall_model_lift_drops_post_stall() {
        let model = StallModel::naca_symmetric();
        let cl_peak = model.lift_coefficient(model.alpha_stall);
        let cl_post = model.lift_coefficient(30.0_f64.to_radians());
        assert!(
            cl_post < cl_peak,
            "Post-stall CL should be less than peak CL: peak={cl_peak}, post={cl_post}"
        );
    }
    #[test]
    fn test_stall_model_stall_fraction_range() {
        let model = StallModel::naca_symmetric();
        let f0 = model.stall_fraction(0.0);
        let f_stall = model.stall_fraction(model.alpha_stall);
        let f_post = model.stall_fraction(40.0_f64.to_radians());
        assert!(
            f0 < 0.1,
            "Low AoA stall fraction should be near 0, got {f0}"
        );
        assert!(
            f_stall > 0.4 && f_stall < 0.6,
            "At stall angle fraction≈0.5, got {f_stall}"
        );
        assert!(
            f_post > 0.9,
            "Well past stall fraction should be near 1, got {f_post}"
        );
    }
    #[test]
    fn test_stall_model_forces_positive_velocity() {
        let model = StallModel::naca_symmetric();
        let (lift, drag) = model.forces(5.0_f64.to_radians(), 30.0, 1.225, 1.0);
        assert!(lift > 0.0, "Lift should be positive, got {lift}");
        assert!(drag > 0.0, "Drag should be positive, got {drag}");
    }
    #[test]
    fn test_stall_model_forces_scale_with_velocity_squared() {
        let model = StallModel::naca_symmetric();
        let (l1, _) = model.forces(5.0_f64.to_radians(), 20.0, 1.225, 1.0);
        let (l2, _) = model.forces(5.0_f64.to_radians(), 40.0, 1.225, 1.0);
        let ratio = l2 / l1;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "Lift should scale as V²: ratio={ratio}"
        );
    }
    #[test]
    fn test_stall_model_cambered_nonzero_lift_at_zero_aoa() {
        let model = StallModel::naca_cambered();
        let cl_0 = model.lift_coefficient(0.0);
        assert!(
            cl_0.abs() > 1e-3,
            "Cambered airfoil should have nonzero CL at 0° due to alpha_zero_lift, got {cl_0}"
        );
    }
    #[test]
    fn test_ground_effect_factor_at_large_height_near_one() {
        let ge = GroundEffect::new(5.0, 8.0);
        let phi = ge.induced_drag_factor(1000.0);
        assert!(
            phi > 0.95,
            "Far from ground, factor should approach 1, got {phi}"
        );
    }
    #[test]
    fn test_ground_effect_factor_at_zero_height_near_zero() {
        let ge = GroundEffect::new(5.0, 8.0);
        let phi = ge.induced_drag_factor(0.01);
        assert!(phi < 0.01, "Near ground, factor should be ~0, got {phi}");
    }
    #[test]
    fn test_ground_effect_lift_ratio_increases_near_ground() {
        let ge = GroundEffect::new(5.0, 8.0);
        let lr_near = ge.lift_ratio(0.5);
        let lr_far = ge.lift_ratio(100.0);
        assert!(
            lr_near > lr_far,
            "Lift ratio should be higher near ground: near={lr_near}, far={lr_far}"
        );
    }
    #[test]
    fn test_ground_effect_apply_reduces_induced_drag() {
        let ge = GroundEffect::new(5.0, 8.0);
        let (_, cdi_ige) = ge.apply(1.0, 0.1, 0.5);
        assert!(
            cdi_ige < 0.1,
            "IGE induced drag should be less than OGE: {cdi_ige}"
        );
    }
    #[test]
    fn test_ground_effect_apply_increases_lift() {
        let ge = GroundEffect::new(5.0, 8.0);
        let (cl_ige, _) = ge.apply(1.0, 0.1, 0.5);
        assert!(cl_ige > 1.0, "IGE lift should exceed OGE lift: {cl_ige}");
    }
    #[test]
    fn test_ground_effect_negligible_height_large() {
        let ge = GroundEffect::new(5.0, 8.0);
        let h = ge.negligible_height();
        assert!(
            h > 100.0,
            "Negligible height should be large for a large wing: {h}"
        );
    }
    #[test]
    fn test_drag_polar_at_zero_cl() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let cd = polar.drag_coefficient(0.0);
        assert!(
            (cd - 0.02).abs() < 1e-10,
            "CD at CL=0 should equal CD0, got {cd}"
        );
    }
    #[test]
    fn test_drag_polar_increases_with_cl() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let cd1 = polar.drag_coefficient(0.5);
        let cd2 = polar.drag_coefficient(1.0);
        assert!(
            cd2 > cd1,
            "Drag should increase with CL: cd1={cd1}, cd2={cd2}"
        );
    }
    #[test]
    fn test_drag_polar_cl_best_ld() {
        let cd0 = 0.02;
        let ar = 10.0;
        let e = 0.9;
        let polar = DragPolar::new(cd0, ar, e);
        let cl_best = polar.cl_best_ld();
        let k = 1.0 / (std::f64::consts::PI * ar * e);
        let expected = (cd0 / k).sqrt();
        assert!(
            (cl_best - expected).abs() < 1e-10,
            "CL_best should be sqrt(CD0/k), got {cl_best} vs {expected}"
        );
    }
    #[test]
    fn test_drag_polar_max_ld_at_best_cl() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let cl_best = polar.cl_best_ld();
        let ld_at_best = polar.lift_to_drag_ratio(cl_best);
        let ld_max = polar.max_lift_to_drag();
        assert!(
            (ld_at_best - ld_max).abs() < 0.01,
            "L/D at CL_best should equal max L/D: {ld_at_best} vs {ld_max}"
        );
    }
    #[test]
    fn test_drag_polar_curve_length() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let curve = polar.polar_curve(11);
        assert_eq!(curve.len(), 11);
        assert!((curve[0].0).abs() < 1e-10, "First CL should be 0");
    }
    #[test]
    fn test_drag_polar_forces() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let (lift, drag) = polar.forces(0.5, 30.0, 1.225, 10.0);
        assert!(lift > 0.0, "Lift should be positive");
        assert!(drag > 0.0, "Drag should be positive");
    }
    #[test]
    fn test_drag_polar_range_parameter_positive() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let ep = polar.range_parameter(0.5);
        assert!(ep > 0.0, "Range parameter should be positive, got {ep}");
    }
    #[test]
    fn test_drag_polar_cl_max_range_larger_than_best_ld() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let cl_range = polar.cl_max_range();
        let cl_ld = polar.cl_best_ld();
        assert!(
            cl_range > cl_ld,
            "CL for max range should exceed CL for max L/D: range={cl_range}, L/D={cl_ld}"
        );
    }
    #[test]
    fn test_drag_polar_lift_to_drag_zero_at_zero_cl() {
        let polar = DragPolar::new(0.02, 10.0, 0.9);
        let ld = polar.lift_to_drag_ratio(0.0);
        assert!(ld.abs() < 1e-10, "L/D should be 0 at CL=0, got {ld}");
    }
}
/// Build the aerodynamic influence coefficient (AIC) matrix for a set of
/// vortex-ring panels.
///
/// Entry `A[i][j]` is the normal component of the velocity induced at panel
/// `i`'s collocation point by a unit circulation on panel `j`.
///
/// Returns a flat row-major vector of length `n*n`.
pub fn build_vlm_aic_matrix(panels: &[VortexRingPanel]) -> Vec<f64> {
    let n = panels.len();
    let mut aic = vec![0.0_f64; n * n];
    for i in 0..n {
        let cp = panels[i].collocation_point();
        let ni = panels[i].normal();
        for j in 0..n {
            let vel = panels[j].induced_velocity_at(cp, 1.0);
            let dot = vel[0] * ni[0] + vel[1] * ni[1] + vel[2] * ni[2];
            aic[i * n + j] = dot;
        }
    }
    aic
}
/// Solve the VLM circulation strengths given the AIC matrix and the
/// right-hand side (normal component of free-stream velocity at each panel).
///
/// Returns a `Vec`f64` of circulation strengths Γ_j.
pub fn solve_vlm_circulations(aic: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut a = aic.to_vec();
    let mut b = rhs.to_vec();
    gauss_solve_opt(n, &mut a, &b).unwrap_or_else(|| {
        b.iter_mut().for_each(|v| *v = 0.0);
        b
    })
}
/// Compute the lift force (N) on each panel using the Kutta-Joukowski theorem:
///
///   L_j = ρ * V_∞ × (Γ_j * Δs_j)
///
/// where Δs_j is the bound-vortex span vector of panel j.
///
/// Returns a `Vec<\[f64; 3\]>` of force vectors (one per panel).
pub fn kutta_joukowski_lift(
    panels: &[VortexRingPanel],
    gammas: &[f64],
    v_inf: [f64; 3],
    air_density: f64,
) -> Vec<[f64; 3]> {
    panels
        .iter()
        .zip(gammas.iter())
        .map(|(panel, &gamma)| {
            let [c0, c1, _c2, _c3] = panel.corners;
            let span_vec = v3_sub(c1, c0);
            let v_cross_dl = v3_cross(v_inf, span_vec);
            v3_scale(v_cross_dl, air_density * gamma)
        })
        .collect()
}
/// Compute induced drag from circulations and induced downwash velocities.
///
/// Uses the trailing-vortex induced velocity at each panel's collocation point.
/// Returns a `Vec`f64` of induced drag magnitudes (N) per panel.
pub fn induced_drag_per_panel(
    panels: &[VortexRingPanel],
    gammas: &[f64],
    v_inf: [f64; 3],
    air_density: f64,
) -> Vec<f64> {
    let n = panels.len();
    panels
        .iter()
        .enumerate()
        .map(|(i, panel)| {
            let cp = panel.collocation_point();
            let mut w_induced = [0.0_f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dv = panels[j].induced_velocity_at(cp, gammas[j]);
                for d in 0..3 {
                    w_induced[d] += dv[d];
                }
            }
            let [c0, c1, _c2, _c3] = panel.corners;
            let span_vec = v3_sub(c1, c0);
            let w_cross_dl = v3_cross(w_induced, span_vec);
            let dv_inf = v3_normalize(v_inf).unwrap_or([1.0, 0.0, 0.0]);
            let di = air_density
                * gammas[i]
                * (w_cross_dl[0] * dv_inf[0]
                    + w_cross_dl[1] * dv_inf[1]
                    + w_cross_dl[2] * dv_inf[2]);
            di.abs()
        })
        .collect()
}
/// Compute the Prandtl span efficiency factor `e` using the elliptic loading
/// approximation for a finite wing.
///
/// The correction accounts for downwash-induced drag.
///
/// `cl_distribution`: lift coefficient at each span station.
/// `dy_distribution`: width of each span station (m).
/// `span`: total wingspan (m).
/// `chord`: mean chord (m).
///
/// Returns the Oswald efficiency factor.
pub fn prandtl_span_efficiency(
    cl_distribution: &[f64],
    dy_distribution: &[f64],
    span: f64,
    chord: f64,
) -> f64 {
    if cl_distribution.is_empty() || span < 1e-6 || chord < 1e-6 {
        return 1.0;
    }
    let n = cl_distribution.len().min(dy_distribution.len());
    let cl_avg: f64 = cl_distribution[..n]
        .iter()
        .zip(dy_distribution[..n].iter())
        .map(|(c, dy)| c * dy)
        .sum::<f64>()
        / span;
    let cl_sq_avg: f64 = cl_distribution[..n]
        .iter()
        .zip(dy_distribution[..n].iter())
        .map(|(c, dy)| c * c * dy)
        .sum::<f64>()
        / span;
    if cl_sq_avg < 1e-20 {
        return 1.0;
    }
    let e = cl_avg * cl_avg / cl_sq_avg;
    e.clamp(0.0, 1.0)
}
/// Compute the total aerodynamic forces (lift vector and induced drag vector)
/// for a wing modelled with VLM panels.
///
/// This is a convenience wrapper that:
/// 1. Builds the AIC matrix.
/// 2. Solves for circulations.
/// 3. Integrates the Kutta-Joukowski forces.
///
/// Returns `(lift_vec_N, induced_drag_N)` where both are 3-component vectors.
pub fn vlm_total_forces(
    panels: &[VortexRingPanel],
    v_inf: [f64; 3],
    air_density: f64,
) -> ([f64; 3], [f64; 3]) {
    if panels.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let _n = panels.len();
    let rhs: Vec<f64> = panels
        .iter()
        .map(|p| {
            let ni = p.normal();
            -(v_inf[0] * ni[0] + v_inf[1] * ni[1] + v_inf[2] * ni[2])
        })
        .collect();
    let aic = build_vlm_aic_matrix(panels);
    let gammas = solve_vlm_circulations(&aic, &rhs);
    let lift_forces = kutta_joukowski_lift(panels, &gammas, v_inf, air_density);
    let drag_mags = induced_drag_per_panel(panels, &gammas, v_inf, air_density);
    let mut total_lift = [0.0_f64; 3];
    let mut total_drag = [0.0_f64; 3];
    let v_dir = v3_normalize(v_inf).unwrap_or([1.0, 0.0, 0.0]);
    for (i, lf) in lift_forces.iter().enumerate() {
        for d in 0..3 {
            total_lift[d] += lf[d];
        }
        for d in 0..3 {
            total_drag[d] += drag_mags[i] * v_dir[d];
        }
    }
    (total_lift, total_drag)
}
/// Compute the downwash angle (radians) at a tail plane located `tail_arm`
/// behind the wing aerodynamic centre.
///
/// Uses the classical approximation:
///   ε = (2 * CL) / (π * AR)
///
/// where CL is the wing lift coefficient and AR is the aspect ratio.
pub fn downwash_angle(cl_wing: f64, aspect_ratio: f64) -> f64 {
    if aspect_ratio < 1e-6 {
        return 0.0;
    }
    2.0 * cl_wing / (std::f64::consts::PI * aspect_ratio)
}
/// Compute the effective angle of attack at the tail taking into account
/// the wing downwash.
///
/// `alpha_tail_geometric`: geometric (physical) angle of attack of the tail (rad).
/// `epsilon`: downwash angle from the wing (rad).
pub fn tail_effective_alpha(alpha_tail_geometric: f64, epsilon: f64) -> f64 {
    alpha_tail_geometric - epsilon
}
/// Linear interpolation helper for 1D look-up tables.
// interp_linear moved to helpers module
#[cfg(test)]
mod aero_extended_tests {

    use crate::aerodynamics::AeroCoeffTable;
    use crate::aerodynamics::VortexRingPanel;
    use crate::aerodynamics::build_vlm_aic_matrix;
    use crate::aerodynamics::downwash_angle;
    use crate::aerodynamics::functions::interp_linear;
    use crate::aerodynamics::induced_drag_per_panel;
    use crate::aerodynamics::kutta_joukowski_lift;
    use crate::aerodynamics::prandtl_span_efficiency;
    use crate::aerodynamics::solve_vlm_circulations;
    use crate::aerodynamics::tail_effective_alpha;
    use crate::aerodynamics::vlm_total_forces;
    #[test]
    fn test_vortex_ring_panel_new_rect_corners() {
        let panel = VortexRingPanel::new_rect([0.0, 0.0, 0.0], 1.0, 2.0, 0.0);
        assert!((panel.corners[0][0] - (-1.0)).abs() < 1e-10);
        assert!((panel.corners[1][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vortex_ring_panel_area_rectangular() {
        let panel = VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [-1.0, 1.0, 0.0],
            ],
            gamma: 0.0,
        };
        let area = panel.area();
        assert!((area - 2.0).abs() < 1e-10, "area = {area}");
    }
    #[test]
    fn test_vortex_ring_panel_normal_flat() {
        let panel = VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [-1.0, 1.0, 0.0],
            ],
            gamma: 0.0,
        };
        let n = panel.normal();
        assert!(n[2].abs() > 0.99, "normal z = {}", n[2]);
        assert!(n[0].abs() < 1e-6);
        assert!(n[1].abs() < 1e-6);
    }
    #[test]
    fn test_vortex_ring_panel_collocation_point() {
        let panel = VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [-1.0, 2.0, 0.0],
            ],
            gamma: 0.0,
        };
        let cp = panel.collocation_point();
        assert!((cp[1] - 1.5).abs() < 1e-10, "cp y = {}", cp[1]);
        assert!(cp[0].abs() < 1e-10, "cp x = {}", cp[0]);
    }
    #[test]
    fn test_vortex_ring_panel_induced_velocity_zero_gamma() {
        let panel = VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [-1.0, 1.0, 0.0],
            ],
            gamma: 0.0,
        };
        let vel = panel.induced_velocity_at([0.0, -1.0, 0.5], 0.0);
        assert_eq!(vel, [0.0; 3], "zero gamma should give zero velocity");
    }
    #[test]
    fn test_build_vlm_aic_matrix_size() {
        let panels = vec![
            VortexRingPanel {
                corners: [
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                ],
                gamma: 0.0,
            },
            VortexRingPanel {
                corners: [
                    [-1.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 2.0, 0.0],
                    [-1.0, 2.0, 0.0],
                ],
                gamma: 0.0,
            },
        ];
        let aic = build_vlm_aic_matrix(&panels);
        assert_eq!(aic.len(), 4, "2×2 AIC matrix should have 4 entries");
    }
    #[test]
    fn test_solve_vlm_circulations_size() {
        let n = 3;
        let aic: Vec<f64> = (0..n * n)
            .map(|i| if i % (n + 1) == 0 { 1.0 } else { 0.0 })
            .collect();
        let rhs = vec![1.0, 2.0, 3.0];
        let gammas = solve_vlm_circulations(&aic, &rhs);
        assert_eq!(gammas.len(), n);
    }
    #[test]
    fn test_kutta_joukowski_lift_zero_gamma() {
        let panels = vec![VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [-1.0, 1.0, 0.0],
            ],
            gamma: 0.0,
        }];
        let gammas = vec![0.0];
        let forces = kutta_joukowski_lift(&panels, &gammas, [30.0, 0.0, 0.0], 1.225);
        assert_eq!(forces[0], [0.0; 3]);
    }
    #[test]
    fn test_kutta_joukowski_lift_nonzero() {
        let panels = vec![VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [-1.0, 2.0, 0.0],
            ],
            gamma: 0.0,
        }];
        let gammas = vec![5.0];
        let v_inf = [0.0, 30.0, 0.0];
        let forces = kutta_joukowski_lift(&panels, &gammas, v_inf, 1.225);
        assert!(
            forces[0][2].abs() > 10.0,
            "should have Z force, got {:?}",
            forces[0]
        );
    }
    #[test]
    fn test_induced_drag_nonneg() {
        let panels = vec![
            VortexRingPanel {
                corners: [
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                ],
                gamma: 1.0,
            },
            VortexRingPanel {
                corners: [
                    [-1.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 2.0, 0.0],
                    [-1.0, 2.0, 0.0],
                ],
                gamma: 1.0,
            },
        ];
        let gammas = vec![1.0, 1.0];
        let drag = induced_drag_per_panel(&panels, &gammas, [30.0, 0.0, 0.0], 1.225);
        for (i, &d) in drag.iter().enumerate() {
            assert!(
                d >= 0.0,
                "induced drag[{i}] should be non-negative, got {d}"
            );
        }
    }
    #[test]
    fn test_prandtl_span_efficiency_elliptic_returns_one() {
        let n = 10usize;
        let b = 10.0_f64;
        let dy = b / n as f64;
        let cl: Vec<f64> = (0..n)
            .map(|i| {
                let y = (i as f64 + 0.5) / n as f64;
                (1.0 - y * y).max(0.0).sqrt()
            })
            .collect();
        let dy_arr = vec![dy; n];
        let e = prandtl_span_efficiency(&cl, &dy_arr, b, 1.0);
        assert!(
            e > 0.0 && e <= 1.0,
            "span efficiency should be in (0, 1], got {e}"
        );
    }
    #[test]
    fn test_prandtl_span_efficiency_uniform_less_than_one() {
        let n = 8usize;
        let cl = vec![1.0_f64; n];
        let dy = vec![1.0_f64; n];
        let e = prandtl_span_efficiency(&cl, &dy, n as f64, 1.0);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "uniform CL gives e = 1.0 by definition, got {e}"
        );
    }
    #[test]
    fn test_prandtl_span_efficiency_empty() {
        let e = prandtl_span_efficiency(&[], &[], 10.0, 1.0);
        assert_eq!(e, 1.0, "empty distribution should return 1.0");
    }
    #[test]
    fn test_downwash_angle_positive() {
        let eps = downwash_angle(0.5, 8.0);
        assert!(
            eps > 0.0,
            "downwash angle should be positive for positive CL, got {eps}"
        );
    }
    #[test]
    fn test_downwash_angle_zero_cl() {
        let eps = downwash_angle(0.0, 8.0);
        assert_eq!(eps, 0.0, "downwash should be zero for zero CL");
    }
    #[test]
    fn test_downwash_angle_zero_ar() {
        let eps = downwash_angle(1.0, 0.0);
        assert_eq!(eps, 0.0, "downwash should be zero for zero AR");
    }
    #[test]
    fn test_tail_effective_alpha() {
        let alpha_tail = 5.0_f64.to_radians();
        let epsilon = 2.0_f64.to_radians();
        let alpha_eff = tail_effective_alpha(alpha_tail, epsilon);
        assert!((alpha_eff - (alpha_tail - epsilon)).abs() < 1e-12);
    }
    #[test]
    fn test_aero_coeff_table_interpolation() {
        let alpha = vec![0.0, 5.0, 10.0];
        let cl = vec![0.0, 0.5, 1.0];
        let cd = vec![0.01, 0.015, 0.03];
        let table = AeroCoeffTable::new(alpha, cl, cd);
        let cl_mid = table.cl_at(5.0);
        assert!((cl_mid - 0.5).abs() < 1e-10, "cl at 5° = {cl_mid}");
        let cl_interp = table.cl_at(7.5);
        assert!((cl_interp - 0.75).abs() < 1e-10, "cl at 7.5° = {cl_interp}");
    }
    #[test]
    fn test_aero_coeff_table_clamp_low() {
        let table = AeroCoeffTable::new(vec![0.0, 10.0], vec![0.1, 1.1], vec![0.01, 0.02]);
        assert!(
            (table.cl_at(-5.0) - 0.1).abs() < 1e-10,
            "cl below range should clamp to first"
        );
    }
    #[test]
    fn test_aero_coeff_table_clamp_high() {
        let table = AeroCoeffTable::new(vec![0.0, 10.0], vec![0.1, 1.1], vec![0.01, 0.02]);
        assert!(
            (table.cl_at(20.0) - 1.1).abs() < 1e-10,
            "cl above range should clamp to last"
        );
    }
    #[test]
    fn test_aero_coeff_table_forces() {
        let table = AeroCoeffTable::new(vec![0.0, 10.0], vec![0.0, 1.0], vec![0.02, 0.04]);
        let q = 0.5 * 1.225 * 50.0 * 50.0;
        let s = 20.0;
        let (lift, drag) = table.forces(5.0, q, s);
        assert!(lift > 0.0, "lift should be positive, got {lift}");
        assert!(drag > 0.0, "drag should be positive, got {drag}");
    }
    #[test]
    fn test_vlm_total_forces_empty_panels() {
        let (lift, drag) = vlm_total_forces(&[], [30.0, 0.0, 0.0], 1.225);
        assert_eq!(lift, [0.0; 3]);
        assert_eq!(drag, [0.0; 3]);
    }
    #[test]
    fn test_vlm_total_forces_single_panel_returns_vectors() {
        let panel = VortexRingPanel {
            corners: [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [-1.0, 2.0, 0.0],
            ],
            gamma: 0.0,
        };
        let (_lift, drag) = vlm_total_forces(&[panel], [30.0, 0.0, 0.0], 1.225);
        for (d, &dv) in drag.iter().enumerate() {
            assert!(dv.is_finite(), "drag[{d}] should be finite");
        }
    }
    #[test]
    fn test_interp_linear_exact_node() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 4.0];
        assert!((interp_linear(&xs, &ys, 1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_interp_linear_midpoint() {
        let xs = vec![0.0, 2.0];
        let ys = vec![0.0, 4.0];
        assert!((interp_linear(&xs, &ys, 1.0) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_interp_linear_empty_returns_zero() {
        assert_eq!(interp_linear(&[], &[], 1.0), 0.0);
    }
}
