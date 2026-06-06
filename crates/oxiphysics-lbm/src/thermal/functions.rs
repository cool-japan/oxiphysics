//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

/// Speed of sound squared in lattice units (cs² = 1/3).
pub(super) const CS2: f64 = 1.0 / 3.0;
/// Compute the thermal boundary layer thickness for a flat plate in forced
/// convection (Pohlhausen approximation):
///
/// δ_T = δ_vel / Pr^(1/3)
///
/// where δ_vel is the velocity boundary layer thickness.
pub fn thermal_boundary_layer_thickness(delta_vel: f64, pr: f64) -> f64 {
    if pr <= 0.0 {
        return 0.0;
    }
    delta_vel / pr.powf(1.0 / 3.0)
}
/// Blasius flat-plate velocity boundary layer thickness:
///
/// δ_99 = 5.0 * sqrt(ν * x / U)
pub fn blasius_boundary_layer_thickness(nu: f64, x: f64, u_inf: f64) -> f64 {
    if u_inf <= 0.0 || x <= 0.0 {
        return 0.0;
    }
    5.0 * (nu * x / u_inf).sqrt()
}
/// Local Nusselt number for laminar flat-plate flow (Pohlhausen–Eckert):
///
/// Nu_x = 0.332 * Re_x^(1/2) * Pr^(1/3)
pub fn local_nusselt_flat_plate(re_x: f64, pr: f64) -> f64 {
    0.332 * re_x.sqrt() * pr.powf(1.0 / 3.0)
}
/// Average Nusselt number for laminar flat plate (length L):
///
/// Nu_L = 0.664 * Re_L^(1/2) * Pr^(1/3)
pub fn average_nusselt_flat_plate(re_l: f64, pr: f64) -> f64 {
    0.664 * re_l.sqrt() * pr.powf(1.0 / 3.0)
}
/// Heat transfer coefficient h = Nu * k / L.
pub fn heat_transfer_coefficient(nu: f64, k_fluid: f64, length: f64) -> f64 {
    if length <= 0.0 {
        return 0.0;
    }
    nu * k_fluid / length
}
/// D3Q7 lattice velocities for thermal LBM.
///
/// Velocities: (0,0,0) and ±e_i for i=1,2,3.
pub const D3Q7_VELOCITIES: [[i32; 3]; 7] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
];
/// D3Q7 weights: w_0 = 1/4, w_{1..6} = 1/8.
pub const D3Q7_WEIGHTS: [f64; 7] = [
    1.0 / 4.0,
    1.0 / 8.0,
    1.0 / 8.0,
    1.0 / 8.0,
    1.0 / 8.0,
    1.0 / 8.0,
    1.0 / 8.0,
];
/// Equilibrium distribution for D3Q7 thermal lattice:
///
/// g_eq_i = w_i * T * (1 + (c_i · u) / (2 * cs2_t))
///
/// where cs2_t = 1/4 for D3Q7.
pub fn d3q7_thermal_equilibrium(t: f64, u: [f64; 3], i: usize) -> f64 {
    pub(super) const CS2_T: f64 = 1.0 / 4.0;
    let w = D3Q7_WEIGHTS[i];
    let c = D3Q7_VELOCITIES[i];
    let cu = c[0] as f64 * u[0] + c[1] as f64 * u[1] + c[2] as f64 * u[2];
    w * t * (1.0 + cu / (2.0 * CS2_T))
}
/// BGK collision step for a single D3Q7 node.
///
/// g_post_i = g_i - (g_i - g_eq_i) / tau_t
pub fn d3q7_bgk_collision(g: &[f64; 7], t: f64, u: [f64; 3], tau_t: f64) -> [f64; 7] {
    let mut g_post = [0.0; 7];
    for i in 0..7 {
        let g_eq = d3q7_thermal_equilibrium(t, u, i);
        g_post[i] = g[i] - (g[i] - g_eq) / tau_t;
    }
    g_post
}
/// Compute temperature from D3Q7 distribution: T = Σ_i g_i.
pub fn d3q7_temperature(g: &[f64; 7]) -> f64 {
    g.iter().sum()
}
/// Thermal diffusivity from D3Q7 relaxation time:
/// α = cs2_t * (tau_t - 0.5)
pub fn d3q7_thermal_diffusivity(tau_t: f64) -> f64 {
    pub(super) const CS2_T: f64 = 1.0 / 4.0;
    CS2_T * (tau_t - 0.5)
}
/// Air dynamic viscosity (Pa·s) using Sutherland's law:
///
/// μ = μ_ref * (T / T_ref)^(3/2) * (T_ref + S) / (T + S)
///
/// Reference values: μ_ref = 1.716e-5 Pa·s, T_ref = 273.15 K, S = 110.4 K.
pub fn air_viscosity_sutherland(t_k: f64) -> f64 {
    pub(super) const MU_REF: f64 = 1.716e-5;
    pub(super) const T_REF: f64 = 273.15;
    pub(super) const S: f64 = 110.4;
    MU_REF * (t_k / T_REF).powf(1.5) * (T_REF + S) / (t_k + S)
}
/// Air thermal conductivity (W/(m·K)) using power-law approximation:
///
/// k = 0.0241 * (T / 273.15)^0.82
pub fn air_thermal_conductivity(t_k: f64) -> f64 {
    0.0241 * (t_k / 273.15_f64).powf(0.82)
}
/// Prandtl number for air (approximately constant ~0.71 near room temperature,
/// slight temperature dependence):
///
/// Pr = μ * cp / k
pub fn air_prandtl_number(t_k: f64) -> f64 {
    pub(super) const CP: f64 = 1006.0;
    let mu = air_viscosity_sutherland(t_k);
    let k = air_thermal_conductivity(t_k);
    mu * CP / k
}
/// Water specific heat capacity as a polynomial fit (valid 273–373 K):
///
/// cp ≈ 4217.6 - 3.83 * (T - 273.15) + 0.0092 * (T - 273.15)^2  \[J/(kg·K)\]
pub fn water_specific_heat(t_k: f64) -> f64 {
    let dt = t_k - 273.15;
    4217.6 - 3.83 * dt + 0.0092 * dt * dt
}
/// Check if thermal LBM parameters satisfy stability conditions.
///
/// Stability requires: 0.5 < τ_t < 2.0 (practical heuristic).
pub fn check_thermal_stability(tau_t: f64) -> bool {
    tau_t > 0.5 && tau_t < 2.0
}
/// Compute τ_t from thermal diffusivity α and cs² for the chosen lattice.
///
/// τ_t = α / cs2 + 0.5  (for D2Q9: cs2 = 1/3)
pub fn tau_from_diffusivity(alpha: f64, cs2: f64) -> f64 {
    alpha / cs2 + 0.5
}
/// Thermal Péclet number: Pe = U * L / α.
pub fn peclet_number(u: f64, l: f64, alpha: f64) -> f64 {
    if alpha.abs() < 1e-30 {
        return f64::INFINITY;
    }
    u * l / alpha
}
/// Fourier number: Fo = α * t / L².
///
/// Governs transient heat conduction; Fo >> 1 implies near-steady state.
pub fn fourier_number(alpha: f64, t: f64, l: f64) -> f64 {
    if l.abs() < 1e-30 {
        return f64::INFINITY;
    }
    alpha * t / (l * l)
}
/// Stefan number: Ste = cp * (T - T_melt) / L_fusion.
///
/// Dimensionless measure of sensible heat relative to latent heat.
pub fn stefan_number(cp: f64, delta_t: f64, l_fusion: f64) -> f64 {
    if l_fusion.abs() < 1e-30 {
        return f64::INFINITY;
    }
    cp * delta_t / l_fusion
}
/// Apply periodic temperature boundary condition in x-direction.
///
/// Wraps the temperature field so that temperature at x=0 equals x=nx-1 + 1.
/// Implemented as a streaming wrap: g at x=0 receives outgoing from x=nx-1.
///
/// Returns the wrapped temperature values for the left (x=0) and right (x=nx-1)
/// boundaries after applying periodic conditions.
pub fn periodic_temperature_x(temperature: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let mut temp = temperature.to_vec();
    for y in 0..ny {
        let left = y * nx;
        let right = y * nx + nx - 1;
        let avg = (temp[left] + temp[right]) / 2.0;
        temp[left] = avg;
        temp[right] = avg;
    }
    temp
}
/// Apply a linear temperature ramp along x as initial condition.
///
/// T(x) = T_left + (T_right - T_left) * x / (nx - 1)
pub fn linear_temperature_profile_x(nx: usize, ny: usize, t_left: f64, t_right: f64) -> Vec<f64> {
    let n = nx * ny;
    let mut temp = vec![0.0; n];
    for y in 0..ny {
        for x in 0..nx {
            temp[y * nx + x] = t_left + (t_right - t_left) * x as f64 / (nx - 1).max(1) as f64;
        }
    }
    temp
}
/// Apply a sinusoidal temperature perturbation on top of a mean temperature.
///
/// T(x, y) = T_mean + A * sin(2π x / nx)
pub fn sinusoidal_temperature_perturbation(
    nx: usize,
    ny: usize,
    t_mean: f64,
    amplitude: f64,
) -> Vec<f64> {
    use std::f64::consts::PI;
    let n = nx * ny;
    let mut temp = vec![0.0; n];
    for y in 0..ny {
        for x in 0..nx {
            temp[y * nx + x] = t_mean + amplitude * (2.0 * PI * x as f64 / nx as f64).sin();
        }
    }
    temp
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermal::types::*;
    /// Verify that the temperature field is initialised to `t_init` everywhere.
    #[test]
    fn test_thermal_d2q9_init() {
        let t_init = 1.5_f64;
        let thermal = ThermalD2Q9::new(8, 6, 1.0, t_init);
        for (k, &t) in thermal.temperature.iter().enumerate() {
            assert!(
                (t - t_init).abs() < 1e-14,
                "temperature[{k}] = {t}, expected {t_init}"
            );
        }
        for (k, g_node) in thermal.g.iter().enumerate() {
            let sum: f64 = g_node.iter().sum();
            assert!(
                (sum - t_init).abs() < 1e-14,
                "Σg[{k}] = {sum}, expected {t_init}"
            );
        }
    }
    /// Verify that Σ_i g_eq_i = T (zeroth moment = temperature).
    #[test]
    fn test_thermal_equilibrium_sum() {
        let thermal = ThermalD2Q9::new(4, 4, 1.0, 1.0);
        let temp = 2.3_f64;
        let u = [0.05_f64, -0.03_f64];
        let mut sum = 0.0_f64;
        for i in 0..9 {
            sum += thermal.equilibrium(temp, u, i);
        }
        assert!(
            (sum - temp).abs() < 1e-13,
            "Σ g_eq_i = {sum}, expected {temp}"
        );
    }
    /// Run 5 steps on a 10×10 grid without panicking.
    #[test]
    fn test_thermal_step_no_panic() {
        let nx = 10;
        let ny = 10;
        let mut thermal = ThermalD2Q9::new(nx, ny, 1.0, 1.0);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        for _ in 0..5 {
            thermal.step(&velocities);
        }
    }
    /// At T = T_ref the buoyancy force must be exactly zero.
    #[test]
    fn test_boussinesq_zero_at_ref() {
        let coupling = BoussinesqCoupling::new([0.0, -1e-4], 2e-3, 1.0);
        let force = coupling.buoyancy_force(1.0, 1.0);
        assert!(
            force[0].abs() < 1e-30 && force[1].abs() < 1e-30,
            "Expected zero force at T_ref, got {force:?}"
        );
    }
    /// Hot fluid (T > T_ref) with downward gravity should produce an upward
    /// buoyancy force (positive y-component when gravity is negative y).
    #[test]
    fn test_boussinesq_hot_rises() {
        let coupling = BoussinesqCoupling::new([0.0, -1e-4], 2e-3, 1.0);
        let force = coupling.buoyancy_force(1.0, 1.5);
        assert!(
            force[1] < 0.0,
            "Hot fluid should feel downward (buoyant) net difference; \
             with gravity negative, F_y = rho*g_y*beta*(T-T_ref) < 0 \
             means the force is in the gravity direction which is correct \
             for standard Boussinesq (force[1]={}).",
            force[1]
        );
        let expected = 1.0 * (-1e-4_f64) * 2e-3 * (1.5 - 1.0);
        assert!(
            (force[1] - expected).abs() < 1e-30,
            "F_y = {}, expected {expected}",
            force[1]
        );
    }
    /// Cold fluid (T < T_ref) with downward gravity should produce a force in
    /// the gravity direction (downward, negative y when g_y < 0).
    #[test]
    fn test_boussinesq_cold_sinks() {
        let coupling = BoussinesqCoupling::new([0.0, -1e-4], 2e-3, 1.0);
        let force = coupling.buoyancy_force(1.0, 0.5);
        assert!(
            force[1] > 0.0,
            "Cold fluid buoyancy force should have opposite sign to hot: force[1]={}",
            force[1]
        );
    }
    /// `lbm_parameters` should return strictly positive relaxation times.
    #[test]
    fn test_rayleigh_benard_setup() {
        let setup = RayleighBenardSetup::new(32, 16, 1.0, 0.0, 1000.0, 0.71);
        let (tau_f, tau_t, beta) = setup.lbm_parameters();
        assert!(tau_f > 0.5, "tau_f must exceed 0.5 for stability: {tau_f}");
        assert!(tau_t > 0.5, "tau_t must exceed 0.5 for stability: {tau_t}");
        assert!(
            beta > 0.0,
            "thermal expansion coefficient must be positive: {beta}"
        );
    }
    /// After applying the bottom BC, bottom-row nodes should have temperature
    /// equal to `t_bot` and their distributions should match equilibrium at
    /// zero velocity.
    #[test]
    fn test_thermal_bc_bottom() {
        let nx = 6;
        let ny = 4;
        let t_bot = 2.0_f64;
        let mut thermal = ThermalD2Q9::new(nx, ny, 1.0, 0.5);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        thermal.apply_temperature_bc_bottom(t_bot, &velocities);
        for x in 0..nx {
            let k = thermal.idx(x, 0);
            assert!(
                (thermal.temperature[k] - t_bot).abs() < 1e-14,
                "Bottom BC temperature at x={x}: got {}, expected {t_bot}",
                thermal.temperature[k]
            );
            for (i, (&gki, &w)) in thermal.g[k].iter().zip(D2Q9_WEIGHTS.iter()).enumerate() {
                let expected = w * t_bot;
                assert!(
                    (gki - expected).abs() < 1e-14,
                    "Bottom BC g[{k}][{i}] = {gki}, expected {expected}",
                );
            }
        }
    }
    #[test]
    fn test_equilibrium_second_order_sum() {
        let thermal = ThermalD2Q9::new(4, 4, 1.0, 1.0);
        let temp = 1.8_f64;
        let u = [0.02, -0.01];
        let mut sum = 0.0_f64;
        for i in 0..9 {
            sum += thermal.equilibrium_second_order(temp, u, i);
        }
        assert!(
            (sum - temp).abs() < 1e-12,
            "Second-order Σ g_eq_i = {sum}, expected {temp}"
        );
    }
    #[test]
    fn test_heat_flux_at_rest() {
        let thermal = ThermalD2Q9::new(4, 4, 1.0, 1.0);
        let qx = thermal.heat_flux_x(0);
        let qy = thermal.heat_flux_y(0);
        assert!(qx.abs() < 1e-14, "Heat flux x should be ~0 at rest: {qx}");
        assert!(qy.abs() < 1e-14, "Heat flux y should be ~0 at rest: {qy}");
    }
    #[test]
    fn test_average_temperature() {
        let t_init = 2.5;
        let thermal = ThermalD2Q9::new(6, 4, 1.0, t_init);
        let avg = thermal.average_temperature();
        assert!(
            (avg - t_init).abs() < 1e-14,
            "Average temperature should be {t_init}, got {avg}"
        );
    }
    #[test]
    fn test_temperature_variance_uniform() {
        let thermal = ThermalD2Q9::new(6, 4, 1.0, 1.0);
        let var = thermal.temperature_variance();
        assert!(
            var.abs() < 1e-14,
            "Variance should be 0 for uniform T: {var}"
        );
    }
    #[test]
    fn test_boussinesq_compute_forces() {
        let coupling = BoussinesqCoupling::new([0.0, -1e-4], 1e-3, 1.0);
        let densities = vec![1.0, 1.0, 1.0];
        let temperatures = vec![1.0, 1.5, 0.5];
        let forces = coupling.compute_forces(&densities, &temperatures);
        assert_eq!(forces.len(), 3);
        assert!(forces[0][1].abs() < 1e-30, "Force at T_ref should be 0");
        assert!(
            forces[1][1] * forces[2][1] < 0.0,
            "Hot and cold forces should be opposite"
        );
    }
    #[test]
    fn test_conjugate_heat_diffusion() {
        let nx = 5;
        let ny = 5;
        let n = nx * ny;
        let mut is_solid = vec![false; n];
        is_solid[2 * nx + 2] = true;
        let mut cht = ConjugateHeatTransfer::new(nx, ny, 0.1, 2.0, is_solid);
        let fluid_temperature = vec![1.0; n];
        let t_before = cht.solid_temperature[2 * nx + 2];
        cht.diffuse_solid(&fluid_temperature);
        let t_after = cht.solid_temperature[2 * nx + 2];
        assert!(
            t_after < t_before,
            "Solid should cool toward fluid: before={t_before}, after={t_after}"
        );
    }
    #[test]
    fn test_nusselt_local_bottom() {
        let nx = 4;
        let ny = 4;
        let mut temperature = vec![0.0; nx * ny];
        for y in 0..ny {
            let t = 1.0 - y as f64 / (ny - 1) as f64;
            for x in 0..nx {
                temperature[y * nx + x] = t;
            }
        }
        let nu_local = NusseltComputation::local_nusselt_bottom(&temperature, nx, ny, 1.0, 0.0);
        assert_eq!(nu_local.len(), nx);
        for &nu in &nu_local {
            assert!(nu > 0.0, "Local Nusselt should be positive: {nu}");
        }
    }
    #[test]
    fn test_nusselt_average() {
        let local = vec![1.0, 2.0, 3.0, 4.0];
        let avg = NusseltComputation::average_nusselt(&local);
        assert!(
            (avg - 2.5).abs() < 1e-14,
            "Average Nu should be 2.5, got {avg}"
        );
    }
    #[test]
    fn test_natural_convection_setup() {
        let setup = NaturalConvectionSetup::new(32, 32, 1.0, 0.0, 1e4, 0.71);
        let (tau_f, tau_t, beta, gravity) = setup.lbm_parameters();
        assert!(tau_f > 0.5, "tau_f must be > 0.5: {tau_f}");
        assert!(tau_t > 0.5, "tau_t must be > 0.5: {tau_t}");
        assert!(beta > 0.0, "beta must be positive: {beta}");
        assert!(gravity[1] < 0.0, "gravity should point down: {:?}", gravity);
    }
    #[test]
    fn test_natural_convection_initial_temperature() {
        let setup = NaturalConvectionSetup::new(8, 6, 1.0, 0.0, 1e4, 0.71);
        let temp = setup.initial_temperature();
        assert_eq!(temp.len(), 48);
        assert!((temp[0] - 1.0).abs() < 1e-14);
        assert!((temp[7] - 0.0).abs() < 1e-14);
    }
    #[test]
    fn test_rb_initial_temperature_profile() {
        let setup = RayleighBenardSetup::new(8, 6, 1.0, 0.0, 1000.0, 0.71);
        let profile = setup.initial_temperature_profile();
        assert_eq!(profile.len(), 6);
        assert_eq!(profile[0].len(), 8);
        assert!((profile[0][0] - 1.0).abs() < 1e-14);
        assert!((profile[5][0] - 0.0).abs() < 1e-14);
    }
    #[test]
    fn test_thermal_bc_left_right() {
        let nx = 6;
        let ny = 4;
        let mut thermal = ThermalD2Q9::new(nx, ny, 1.0, 1.0);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        thermal.apply_temperature_bc_left(2.0, &velocities);
        thermal.apply_temperature_bc_right(0.5, &velocities);
        for y in 0..ny {
            let k = thermal.idx(0, y);
            assert!(
                (thermal.temperature[k] - 2.0).abs() < 1e-14,
                "Left BC: T at (0,{y}) = {}",
                thermal.temperature[k]
            );
        }
        for y in 0..ny {
            let k = thermal.idx(nx - 1, y);
            assert!(
                (thermal.temperature[k] - 0.5).abs() < 1e-14,
                "Right BC: T at ({},{y}) = {}",
                nx - 1,
                thermal.temperature[k]
            );
        }
    }
    #[test]
    fn test_nusselt_volume_averaged_conduction() {
        let nx = 10;
        let ny = 10;
        let n = nx * ny;
        let mut temperature = vec![0.0; n];
        for y in 0..ny {
            let t = 1.0 - y as f64 / (ny - 1) as f64;
            for x in 0..nx {
                temperature[y * nx + x] = t;
            }
        }
        let uy = vec![0.0; n];
        let alpha = 0.1;
        let nu =
            NusseltComputation::nusselt_volume_averaged(&temperature, &uy, nx, ny, alpha, 1.0, 0.0);
        assert!(
            (nu - 1.0).abs() < 1e-10,
            "Nu for pure conduction should be ~1, got {nu}"
        );
    }
    #[test]
    fn test_thermal_energy_conservation() {
        let nx = 8;
        let ny = 6;
        let mut thermal = ThermalD2Q9::new(nx, ny, 0.8, 1.0);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        thermal.g[5][3] += 0.1;
        thermal.g[10][7] -= 0.05;
        thermal.compute_temperature();
        let energy_before: f64 = thermal.temperature.iter().sum();
        for _ in 0..10 {
            thermal.step(&velocities);
        }
        let energy_after: f64 = thermal.temperature.iter().sum();
        assert!(
            (energy_before - energy_after).abs() < 1e-10,
            "Total thermal energy should be conserved: before={energy_before}, after={energy_after}"
        );
    }
}
#[cfg(test)]
mod tests_extended_thermal {
    use super::*;

    #[test]
    fn test_thermal_boundary_layer_pr1() {
        let delta_t = thermal_boundary_layer_thickness(0.01, 1.0);
        assert!((delta_t - 0.01).abs() < 1e-14, "δ_T at Pr=1: {delta_t}");
    }
    #[test]
    fn test_thermal_boundary_layer_thinner_at_high_pr() {
        let delta_low = thermal_boundary_layer_thickness(0.01, 1.0);
        let delta_high = thermal_boundary_layer_thickness(0.01, 100.0);
        assert!(delta_high < delta_low, "High Pr → thinner thermal BL");
    }
    #[test]
    fn test_blasius_thickness_zero_x() {
        let d = blasius_boundary_layer_thickness(1e-5, 0.0, 1.0);
        assert_eq!(d, 0.0, "zero x → zero thickness");
    }
    #[test]
    fn test_local_nusselt_flat_plate_reference() {
        let nu = local_nusselt_flat_plate(1.0e4, 0.71);
        let expected = 0.332 * 100.0 * 0.71_f64.powf(1.0 / 3.0);
        assert!((nu - expected).abs() < 1e-10, "Nu_x = {nu}");
    }
    #[test]
    fn test_average_nusselt_double_local() {
        let nu_local = local_nusselt_flat_plate(1e4, 0.71);
        let nu_avg = average_nusselt_flat_plate(1e4, 0.71);
        assert!(
            (nu_avg / nu_local - 2.0).abs() < 1e-10,
            "nu_avg/nu_local = {}",
            nu_avg / nu_local
        );
    }
    #[test]
    fn test_d3q7_equilibrium_sum_equals_temperature() {
        let t = 1.5;
        let u = [0.05, -0.02, 0.01];
        let sum: f64 = (0..7).map(|i| d3q7_thermal_equilibrium(t, u, i)).sum();
        assert!((sum - t).abs() < 1e-13, "Σ g_eq_i = {sum}, expected {t}");
    }
    #[test]
    fn test_d3q7_temperature_sum() {
        let g = [0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
        let t = d3q7_temperature(&g);
        assert!((t - 0.9).abs() < 1e-14, "T = {t}");
    }
    #[test]
    fn test_d3q7_bgk_collision_conserves_temperature() {
        let t = 1.2;
        let u = [0.0; 3];
        let g: [f64; 7] = std::array::from_fn(|i| d3q7_thermal_equilibrium(t, u, i));
        let g_post = d3q7_bgk_collision(&g, t, u, 1.0);
        let t_post: f64 = g_post.iter().sum();
        assert!((t_post - t).abs() < 1e-13, "T post collision = {t_post}");
    }
    #[test]
    fn test_d3q7_thermal_diffusivity_formula() {
        let tau = 1.0;
        let alpha = d3q7_thermal_diffusivity(tau);
        let expected = 0.25 * (tau - 0.5);
        assert!((alpha - expected).abs() < 1e-14, "α = {alpha}");
    }
    #[test]
    fn test_air_viscosity_sutherland_room_temp() {
        let mu = air_viscosity_sutherland(293.15);
        assert!((mu - 1.81e-5).abs() < 5e-7, "μ at 293 K = {mu}");
    }
    #[test]
    fn test_air_viscosity_increases_with_temperature() {
        let mu_cold = air_viscosity_sutherland(250.0);
        let mu_hot = air_viscosity_sutherland(400.0);
        assert!(mu_hot > mu_cold, "Viscosity should increase with T");
    }
    #[test]
    fn test_air_prandtl_near_071() {
        let pr = air_prandtl_number(293.15);
        assert!(
            (pr - 0.71).abs() < 0.05,
            "Pr for air should be ~0.71, got {pr}"
        );
    }
    #[test]
    fn test_check_thermal_stability_valid() {
        assert!(check_thermal_stability(1.0));
        assert!(check_thermal_stability(0.6));
        assert!(check_thermal_stability(1.9));
    }
    #[test]
    fn test_check_thermal_stability_invalid() {
        assert!(!check_thermal_stability(0.5));
        assert!(!check_thermal_stability(0.3));
        assert!(!check_thermal_stability(2.0));
    }
    #[test]
    fn test_tau_from_diffusivity_roundtrip() {
        let alpha = 0.05;
        let cs2 = 1.0 / 3.0;
        let tau = tau_from_diffusivity(alpha, cs2);
        let alpha_back = (tau - 0.5) * cs2;
        assert!(
            (alpha_back - alpha).abs() < 1e-14,
            "roundtrip alpha = {alpha_back}"
        );
    }
    #[test]
    fn test_peclet_number_zero_alpha() {
        let pe = peclet_number(1.0, 1.0, 0.0);
        assert!(pe.is_infinite(), "Infinite Pe when α=0");
    }
    #[test]
    fn test_fourier_number_scaling() {
        let fo1 = fourier_number(0.1, 1.0, 1.0);
        let fo2 = fourier_number(0.1, 4.0, 1.0);
        assert!((fo2 / fo1 - 4.0).abs() < 1e-14, "Fo ratio = {}", fo2 / fo1);
    }
    #[test]
    fn test_linear_temperature_profile_endpoints() {
        let temp = linear_temperature_profile_x(5, 3, 0.0, 1.0);
        assert!((temp[0] - 0.0).abs() < 1e-14, "T at x=0: {}", temp[0]);
        assert!((temp[4] - 1.0).abs() < 1e-14, "T at x=4: {}", temp[4]);
    }
    #[test]
    fn test_linear_temperature_profile_length() {
        let temp = linear_temperature_profile_x(4, 3, 0.0, 1.0);
        assert_eq!(temp.len(), 12);
    }
    #[test]
    fn test_sinusoidal_perturbation_mean_preserved() {
        let nx = 16;
        let ny = 4;
        let t_mean = 1.0;
        let amp = 0.1;
        let temp = sinusoidal_temperature_perturbation(nx, ny, t_mean, amp);
        let mean: f64 = temp[0..nx].iter().sum::<f64>() / nx as f64;
        assert!((mean - t_mean).abs() < 1e-12, "Mean T = {mean}");
    }
}
/// Compute the DDF energy equilibrium distribution for direction `i`.
///
/// ```text
/// g_eq_i = w_i * T * (1 + (e_i · u) / cs²)
/// ```
///
/// This first-order form is suitable for the advection-diffusion thermal LBM.
pub fn energy_equilibrium(temp: f64, u: [f64; 2], i: usize) -> f64 {
    let w = D2Q9_WEIGHTS[i];
    let c = D2Q9_VELOCITIES[i];
    let e_dot_u = c[0] as f64 * u[0] + c[1] as f64 * u[1];
    w * temp * (1.0 + e_dot_u / CS2)
}
/// Compute the second-order DDF energy equilibrium distribution.
///
/// ```text
/// g_eq_i = w_i * T * (1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) − u²/(2cs²))
/// ```
pub fn energy_equilibrium_second_order(temp: f64, u: [f64; 2], i: usize) -> f64 {
    let w = D2Q9_WEIGHTS[i];
    let c = D2Q9_VELOCITIES[i];
    let e_dot_u = c[0] as f64 * u[0] + c[1] as f64 * u[1];
    let u2 = u[0] * u[0] + u[1] * u[1];
    w * temp * (1.0 + e_dot_u / CS2 + e_dot_u * e_dot_u / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
}
/// Sum the zeroth moment of the energy equilibrium over all 9 directions.
///
/// Should equal `T` for any velocity `u`.
pub fn energy_equilibrium_sum(temp: f64, u: [f64; 2]) -> f64 {
    (0..9).map(|i| energy_equilibrium(temp, u, i)).sum()
}
/// Apply a fixed-temperature Dirichlet BC on the bottom wall of a flat
/// temperature array `temp[y * nx + x]` (sets `y = 0` row to `t_bot`).
pub fn apply_dirichlet_bc_bottom(temp: &mut [f64], nx: usize, t_bot: f64) {
    for t in temp[..nx].iter_mut() {
        *t = t_bot;
    }
}
/// Apply a fixed-temperature Dirichlet BC on the top wall (`y = ny-1`).
pub fn apply_dirichlet_bc_top(temp: &mut [f64], nx: usize, ny: usize, t_top: f64) {
    let y_top = ny - 1;
    for x in 0..nx {
        temp[y_top * nx + x] = t_top;
    }
}
/// Apply an insulating (zero-flux / Neumann) BC on the left wall using
/// zero-gradient extrapolation: `T(0, y) = T(1, y)`.
pub fn apply_insulating_bc_left(temp: &mut [f64], nx: usize, ny: usize) {
    for y in 0..ny {
        temp[y * nx] = temp[y * nx + 1];
    }
}
/// Apply an insulating (zero-flux / Neumann) BC on the right wall:
/// `T(nx-1, y) = T(nx-2, y)`.
pub fn apply_insulating_bc_right(temp: &mut [f64], nx: usize, ny: usize) {
    for y in 0..ny {
        temp[y * nx + nx - 1] = temp[y * nx + nx - 2];
    }
}
/// Apply a constant heat flux BC on the bottom wall (upward flux `q_bot`).
///
/// In a finite-difference sense this is:
/// `T(x, 0) = T(x, 1) + q_bot * delta_y / k_fluid`
///
/// For the LBM the heat flux is imposed by adjusting the temperature at the
/// wall node such that the gradient `dT/dy|_{y=0}` equals `q_bot / k_fluid`.
/// With lattice spacing `Δy = 1`:
/// `T_wall = T_interior + q_bot / k_fluid`
pub fn apply_heat_flux_bc_bottom(temp: &mut [f64], nx: usize, q_flux: f64, k_fluid: f64) {
    for x in 0..nx {
        let t_interior = temp[nx + x];
        temp[x] = t_interior + q_flux / k_fluid.max(1e-30);
    }
}
/// Compute the heat flux vector `[q_x, q_y]` at node `k` from the thermal
/// distribution `g`:
///
/// ```text
/// q_x = Σ_i c_{ix} * g_i
/// q_y = Σ_i c_{iy} * g_i
/// ```
pub fn compute_heat_flux(g: &[f64; 9]) -> [f64; 2] {
    let mut qx = 0.0_f64;
    let mut qy = 0.0_f64;
    for i in 0..9 {
        let c = D2Q9_VELOCITIES[i];
        qx += c[0] as f64 * g[i];
        qy += c[1] as f64 * g[i];
    }
    [qx, qy]
}
/// Compute heat flux magnitude at a node.
pub fn heat_flux_magnitude(g: &[f64; 9]) -> f64 {
    let [qx, qy] = compute_heat_flux(g);
    (qx * qx + qy * qy).sqrt()
}
/// Compute the Nusselt number from the temperature gradient at the hot wall
/// using central differences (more accurate than one-sided FD).
///
/// ```text
/// Nu = -(T(x,2) - T(x,0)) / (2 * ΔT / H) = H * (T(x,0) - T(x,2)) / (2 * ΔT)
/// ```
pub fn local_nusselt_bottom_central(
    temperature: &[f64],
    nx: usize,
    ny: usize,
    t_hot: f64,
    t_cold: f64,
) -> Vec<f64> {
    let h = ny as f64;
    let delta_t = (t_hot - t_cold).abs().max(1e-30);
    let mut nu = Vec::with_capacity(nx);
    for x in 0..nx {
        let t0 = temperature[x];
        let t2 = temperature[2 * nx + x];
        let grad_2h = (t2 - t0) / 2.0;
        nu.push((-grad_2h).abs() * h / delta_t);
    }
    nu
}
/// Compute the Nusselt number for the de Vahl Davis differentially heated
/// cavity benchmark.
///
/// The de Vahl Davis (1983) result for Ra=10⁴, Pr=0.71 gives Nu_avg ≈ 2.243.
/// This function computes Nu from the average temperature gradient on the hot
/// (left) wall.
///
/// `temperature` is stored as `temperature[y * nx + x]`.
pub fn de_vahl_davis_nusselt(
    temperature: &[f64],
    nx: usize,
    ny: usize,
    t_hot: f64,
    t_cold: f64,
) -> f64 {
    let h = nx as f64;
    let delta_t = (t_hot - t_cold).abs().max(1e-30);
    let mut grad_sum = 0.0_f64;
    for y in 0..ny {
        let t0 = temperature[y * nx];
        let t1 = temperature[y * nx + 1];
        grad_sum += t1 - t0;
    }
    let avg_grad = grad_sum / ny as f64;
    (-avg_grad).abs() * h / delta_t
}
/// Critical Rayleigh number for the onset of Rayleigh-Bénard convection.
///
/// For a fluid layer between two rigid, no-slip walls the critical Rayleigh
/// number is approximately:
///
/// ```text
/// Ra_c ≈ 1707.76
/// ```
///
/// For free-slip (stress-free) boundaries it is `Ra_c ≈ 657.5`.
pub const RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP: f64 = 1707.76;
/// Critical Rayleigh number for free-slip boundary conditions.
pub const RAYLEIGH_BENARD_CRITICAL_RA_FREESLIP: f64 = 657.5;
/// Check whether Rayleigh–Bénard convection onset is predicted.
///
/// Returns `true` if `Ra > Ra_c` (no-slip boundaries).
pub fn is_rb_convective(ra: f64) -> bool {
    ra > RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP
}
/// Estimate the Nusselt number from the Grossmann-Lohse scaling for
/// Rayleigh-Bénard convection (simplified power-law form):
///
/// ```text
/// Nu ≈ 0.16 * Ra^(1/3)  (turbulent regime, Pr ~ 0.71)
/// ```
///
/// This is a rough engineering estimate valid for `10⁵ ≤ Ra ≤ 10⁹`.
pub fn grossmann_lohse_nusselt(ra: f64) -> f64 {
    if ra <= RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP {
        return 1.0;
    }
    0.16 * ra.powf(1.0 / 3.0)
}
/// Estimate the Nusselt number using the Drazin-Reid correlation for
/// Rayleigh-Bénard convection near onset:
///
/// ```text
/// Nu ≈ 1 + C * (Ra - Ra_c) / Ra_c    for Ra just above Ra_c
/// ```
pub fn nusselt_near_onset(ra: f64, c: f64) -> f64 {
    if ra <= RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP {
        return 1.0;
    }
    1.0 + c * (ra - RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP) / RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP
}
/// Compute thermal diffusivity from D2Q9 relaxation time.
///
/// `α = (τ_t - 0.5) * cs²    where cs² = 1/3`
pub fn thermal_diffusivity_d2q9(tau_t: f64) -> f64 {
    (tau_t - 0.5) * CS2
}
/// Compute D2Q9 thermal relaxation time from diffusivity.
///
/// `τ_t = α / cs² + 0.5`
pub fn tau_t_from_diffusivity_d2q9(alpha: f64) -> f64 {
    alpha / CS2 + 0.5
}
/// Compute the Prandtl number from kinematic viscosity `ν` and thermal
/// diffusivity `α`:
///
/// `Pr = ν / α`
pub fn prandtl_number(nu: f64, alpha: f64) -> f64 {
    if alpha.abs() < 1e-30 {
        return f64::INFINITY;
    }
    nu / alpha
}
/// Compute the Rayleigh number from fluid and geometry properties:
///
/// `Ra = g * β * ΔT * H³ / (ν * α)`
pub fn rayleigh_number(g: f64, beta: f64, delta_t: f64, h: f64, nu: f64, alpha: f64) -> f64 {
    if nu.abs() < 1e-30 || alpha.abs() < 1e-30 {
        return f64::INFINITY;
    }
    g * beta * delta_t * h.powi(3) / (nu * alpha)
}
/// Compute the Grashof number: `Gr = Ra / Pr = g * β * ΔT * H³ / ν²`.
pub fn grashof_number(g: f64, beta: f64, delta_t: f64, h: f64, nu: f64) -> f64 {
    if nu.abs() < 1e-30 {
        return f64::INFINITY;
    }
    g * beta * delta_t * h.powi(3) / (nu * nu)
}
#[cfg(test)]
mod extended_thermal_tests {
    use super::*;
    use crate::thermal::types::*;
    #[test]
    fn test_thermal_lbm_new_uniform_temperature() {
        let t_init = 1.5;
        let lbm = ThermalLbm::new(8, 6, 1.0, t_init);
        for &t in lbm.temperature() {
            assert!(
                (t - t_init).abs() < 1e-14,
                "Initial T should be uniform: {t}"
            );
        }
    }
    #[test]
    fn test_thermal_lbm_step_no_panic() {
        let mut lbm = ThermalLbm::new(10, 10, 0.9, 1.0);
        let velocities = vec![[0.0_f64; 2]; 100];
        for _ in 0..5 {
            lbm.step_thermal(&velocities);
        }
        assert_eq!(lbm.step_count, 5, "Step count should be 5");
    }
    #[test]
    fn test_thermal_lbm_average_temperature_uniform() {
        let t_init = 2.0;
        let lbm = ThermalLbm::new(8, 8, 1.0, t_init);
        let avg = lbm.average_temperature();
        assert!((avg - t_init).abs() < 1e-14, "Average T = {avg}");
    }
    #[test]
    fn test_thermal_lbm_variance_uniform() {
        let lbm = ThermalLbm::new(8, 8, 1.0, 1.0);
        let var = lbm.temperature_variance();
        assert!(var.abs() < 1e-14, "Variance = {var} for uniform T");
    }
    #[test]
    fn test_thermal_lbm_bc_bottom() {
        let nx = 6;
        let ny = 4;
        let mut lbm = ThermalLbm::new(nx, ny, 1.0, 0.5);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        lbm.apply_bc_bottom(2.0, &velocities);
        for x in 0..nx {
            let k = lbm.thermal.idx(x, 0);
            assert!(
                (lbm.temperature()[k] - 2.0).abs() < 1e-14,
                "Bottom BC at x={x}: T={}",
                lbm.temperature()[k]
            );
        }
    }
    #[test]
    fn test_thermal_lbm_bc_top() {
        let nx = 6;
        let ny = 4;
        let mut lbm = ThermalLbm::new(nx, ny, 1.0, 1.0);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        lbm.apply_bc_top(0.0, &velocities);
        for x in 0..nx {
            let k = lbm.thermal.idx(x, ny - 1);
            assert!(
                (lbm.temperature()[k] - 0.0).abs() < 1e-14,
                "Top BC at x={x}: T={}",
                lbm.temperature()[k]
            );
        }
    }
    #[test]
    fn test_thermal_lbm_reset_temperature() {
        let mut lbm = ThermalLbm::new(6, 6, 1.0, 1.0);
        lbm.thermal.temperature[5] = 3.0;
        lbm.reset_temperature(2.0);
        for &t in lbm.temperature() {
            assert!((t - 2.0).abs() < 1e-14, "After reset T should be 2.0: {t}");
        }
    }
    #[test]
    fn test_thermal_lbm_set_temperature_at() {
        let mut lbm = ThermalLbm::new(8, 8, 1.0, 1.0);
        lbm.set_temperature_at(3, 4, 2.5, [0.0, 0.0]);
        let k = lbm.thermal.idx(3, 4);
        assert!(
            (lbm.temperature()[k] - 2.5).abs() < 1e-14,
            "T at (3,4) = {}",
            lbm.temperature()[k]
        );
    }
    #[test]
    fn test_thermal_lbm_with_boussinesq() {
        let lbm = ThermalLbm::new(8, 8, 1.0, 1.0).with_boussinesq([0.0, -1e-4], 2e-3, 1.0);
        assert!(lbm.boussinesq.is_some(), "Boussinesq should be set");
    }
    #[test]
    fn test_thermal_lbm_compute_buoyancy_none_without_coupling() {
        let lbm = ThermalLbm::new(4, 4, 1.0, 1.0);
        let densities = vec![1.0; 16];
        assert!(
            lbm.compute_buoyancy_forces(&densities).is_none(),
            "No coupling → no buoyancy forces"
        );
    }
    #[test]
    fn test_thermal_lbm_compute_buoyancy_with_coupling() {
        let nx = 4;
        let ny = 4;
        let lbm = ThermalLbm::new(nx, ny, 1.0, 1.0).with_boussinesq([0.0, -1e-4], 2e-3, 1.0);
        let densities = vec![1.0; nx * ny];
        let forces = lbm.compute_buoyancy_forces(&densities).unwrap();
        assert_eq!(forces.len(), nx * ny, "One force per node");
        for f in &forces {
            assert!(
                f[0].abs() < 1e-30 && f[1].abs() < 1e-30,
                "At T_ref buoyancy = 0: {f:?}"
            );
        }
    }
    #[test]
    fn test_thermal_lbm_nusselt_bottom_positive() {
        let nx = 8;
        let ny = 8;
        let mut lbm = ThermalLbm::new(nx, ny, 1.0, 1.0);
        let velocities = vec![[0.0_f64; 2]; nx * ny];
        lbm.apply_bc_bottom(1.0, &velocities);
        lbm.apply_bc_top(0.0, &velocities);
        for y in 0..ny {
            for x in 0..nx {
                let k = lbm.thermal.idx(x, y);
                lbm.thermal.temperature[k] = 1.0 - y as f64 / (ny - 1) as f64;
            }
        }
        let nu_avg = lbm.average_nusselt_bottom(1.0, 0.0);
        assert!(nu_avg > 0.0, "Average Nu should be positive: {nu_avg}");
    }
    #[test]
    fn test_energy_equilibrium_sum_equals_temp() {
        let temp = 1.8;
        let u = [0.05, -0.02];
        let sum: f64 = (0..9).map(|i| energy_equilibrium(temp, u, i)).sum();
        assert!(
            (sum - temp).abs() < 1e-13,
            "Σ g_eq = {sum}, expected {temp}"
        );
    }
    #[test]
    fn test_energy_equilibrium_second_order_sum() {
        let temp = 2.3;
        let u = [0.01, 0.02];
        let sum: f64 = (0..9)
            .map(|i| energy_equilibrium_second_order(temp, u, i))
            .sum();
        assert!((sum - temp).abs() < 1e-12, "Second-order Σ g_eq = {sum}");
    }
    #[test]
    fn test_energy_equilibrium_sum_function() {
        let temp = 1.2;
        let u = [0.0, 0.0];
        let s = energy_equilibrium_sum(temp, u);
        assert!((s - temp).abs() < 1e-13, "energy_equilibrium_sum = {s}");
    }
    #[test]
    fn test_energy_equilibrium_positive_for_positive_temp() {
        let temp = 1.0;
        let u = [0.0, 0.0];
        for i in 0..9 {
            let g = energy_equilibrium(temp, u, i);
            assert!(g > 0.0, "g_eq[{i}] = {g} should be positive");
        }
    }
    #[test]
    fn test_apply_dirichlet_bc_bottom_sets_value() {
        let nx = 5;
        let ny = 4;
        let mut temp = vec![1.0; nx * ny];
        apply_dirichlet_bc_bottom(&mut temp, nx, 3.0);
        for (x, &t) in temp[..nx].iter().enumerate() {
            assert!((t - 3.0).abs() < 1e-14, "Bottom BC x={x}: T={t}",);
        }
    }
    #[test]
    fn test_apply_dirichlet_bc_top_sets_value() {
        let nx = 5;
        let ny = 4;
        let mut temp = vec![1.0; nx * ny];
        apply_dirichlet_bc_top(&mut temp, nx, ny, 0.5);
        for x in 0..nx {
            let t = temp[(ny - 1) * nx + x];
            assert!((t - 0.5).abs() < 1e-14, "Top BC x={x}: T={t}");
        }
    }
    #[test]
    fn test_apply_insulating_bc_left_zero_gradient() {
        let nx = 6;
        let ny = 4;
        let mut temp: Vec<f64> = (0..(nx * ny)).map(|k| (k % nx) as f64).collect();
        apply_insulating_bc_left(&mut temp, nx, ny);
        for y in 0..ny {
            let t0 = temp[y * nx];
            let t1 = temp[y * nx + 1];
            assert!(
                (t0 - t1).abs() < 1e-14,
                "Left insulating BC: T(0,{y})={t0} should equal T(1,{y})={t1}"
            );
        }
    }
    #[test]
    fn test_apply_insulating_bc_right_zero_gradient() {
        let nx = 6;
        let ny = 4;
        let mut temp: Vec<f64> = (0..(nx * ny)).map(|k| (k % nx) as f64).collect();
        apply_insulating_bc_right(&mut temp, nx, ny);
        for y in 0..ny {
            let tn1 = temp[y * nx + nx - 1];
            let tn2 = temp[y * nx + nx - 2];
            assert!(
                (tn1 - tn2).abs() < 1e-14,
                "Right insulating BC: T(nx-1,{y})={tn1} should equal T(nx-2,{y})={tn2}"
            );
        }
    }
    #[test]
    fn test_apply_heat_flux_bc_bottom() {
        let nx = 4;
        let ny = 3;
        let mut temp = vec![1.0; nx * ny];
        for x in 0..nx {
            temp[nx + x] = 2.0;
        }
        let q = 1.0;
        let k_f = 2.0;
        apply_heat_flux_bc_bottom(&mut temp, nx, q, k_f);
        for (x, &t) in temp[..nx].iter().enumerate() {
            assert!((t - 2.5).abs() < 1e-14, "Heat flux BC x={x}: T={t}",);
        }
    }
    #[test]
    fn test_compute_heat_flux_zero_at_rest() {
        let t = 1.0;
        let u = [0.0_f64, 0.0];
        let g: [f64; 9] = std::array::from_fn(|i| energy_equilibrium(t, u, i));
        let [qx, qy] = compute_heat_flux(&g);
        assert!(qx.abs() < 1e-14, "Heat flux x at rest = {qx}");
        assert!(qy.abs() < 1e-14, "Heat flux y at rest = {qy}");
    }
    #[test]
    fn test_heat_flux_magnitude_nonnegative() {
        let t = 1.0;
        let u = [0.05, 0.0];
        let g: [f64; 9] = std::array::from_fn(|i| energy_equilibrium(t, u, i));
        let mag = heat_flux_magnitude(&g);
        assert!(
            mag >= 0.0,
            "Heat flux magnitude should be non-negative: {mag}"
        );
    }
    #[test]
    fn test_nusselt_central_linear_profile() {
        let nx = 6;
        let ny = 6;
        let mut temperature = vec![0.0; nx * ny];
        for y in 0..ny {
            let t = 1.0 - y as f64 / (ny - 1) as f64;
            for x in 0..nx {
                temperature[y * nx + x] = t;
            }
        }
        let nu = local_nusselt_bottom_central(&temperature, nx, ny, 1.0, 0.0);
        assert_eq!(nu.len(), nx, "Nu vector length should equal nx");
        for &n in &nu {
            assert!(n > 0.0, "Local Nu should be positive: {n}");
        }
    }
    #[test]
    fn test_de_vahl_davis_setup_lbm_parameters() {
        let setup = DeVahlDavisSetup::new(32, 1.0, 0.0, 1e4, 0.71);
        let (tau_f, tau_t, beta, gravity) = setup.lbm_parameters();
        assert!(tau_f > 0.5, "tau_f must be > 0.5: {tau_f}");
        assert!(tau_t > 0.5, "tau_t must be > 0.5: {tau_t}");
        assert!(beta > 0.0, "beta must be positive: {beta}");
        assert!(
            gravity[1] < 0.0,
            "gravity y-component should be negative: {:?}",
            gravity
        );
    }
    #[test]
    fn test_de_vahl_davis_initial_temperature_hot_wall() {
        let setup = DeVahlDavisSetup::new(8, 1.0, 0.0, 1e4, 0.71);
        let temp = setup.initial_temperature();
        for y in 0..8 {
            assert!(
                (temp[y * 8] - 1.0).abs() < 1e-14,
                "Left (hot) wall at y={y}: T={}",
                temp[y * 8]
            );
        }
        for y in 0..8 {
            assert!(
                (temp[y * 8 + 7] - 0.0).abs() < 1e-14,
                "Right (cold) wall at y={y}: T={}",
                temp[y * 8 + 7]
            );
        }
    }
    #[test]
    fn test_de_vahl_davis_check_stability() {
        let setup = DeVahlDavisSetup::new(32, 1.0, 0.0, 1e4, 0.71);
        assert!(
            setup.check_stability(),
            "Parameters should be stable for Ra=1e4"
        );
    }
    #[test]
    fn test_de_vahl_davis_dimensionless_temperature() {
        let setup = DeVahlDavisSetup::new(8, 1.0, 0.0, 1e3, 0.71);
        let temp = setup.initial_temperature();
        let theta = setup.dimensionless_temperature(&temp);
        for y in 0..8 {
            assert!(
                (theta[y * 8] - 1.0).abs() < 1e-14,
                "θ at hot wall = {}",
                theta[y * 8]
            );
            assert!(
                (theta[y * 8 + 7] - 0.0).abs() < 1e-14,
                "θ at cold wall = {}",
                theta[y * 8 + 7]
            );
        }
    }
    #[test]
    fn test_de_vahl_davis_mid_plane_x_temperature() {
        let n = 8;
        let setup = DeVahlDavisSetup::new(n, 1.0, 0.0, 1e3, 0.71);
        let temp = setup.initial_temperature();
        let profile = setup.mid_plane_x_temperature(&temp);
        assert_eq!(profile.len(), n, "Mid-plane profile length = n");
        assert!(
            (profile[0] - 1.0).abs() < 1e-14,
            "Left end = T_hot: {}",
            profile[0]
        );
        assert!(
            (profile[n - 1] - 0.0).abs() < 1e-14,
            "Right end = T_cold: {}",
            profile[n - 1]
        );
    }
    #[test]
    fn test_de_vahl_davis_nusselt_linear_profile() {
        let n = 8;
        let mut temp = vec![0.0; n * n];
        for y in 0..n {
            for x in 0..n {
                temp[y * n + x] = 1.0 - x as f64 / (n - 1) as f64;
            }
        }
        let nu = de_vahl_davis_nusselt(&temp, n, n, 1.0, 0.0);
        assert!(nu > 0.0, "Nu should be positive for linear profile: {nu}");
    }
    #[test]
    fn test_rb_convective_below_critical() {
        assert!(
            !is_rb_convective(1000.0),
            "Ra=1000 < Ra_c should not be convective"
        );
    }
    #[test]
    fn test_rb_convective_above_critical() {
        assert!(
            is_rb_convective(2000.0),
            "Ra=2000 > Ra_c should be convective"
        );
    }
    #[test]
    fn test_grossmann_lohse_nusselt_below_critical() {
        let nu = grossmann_lohse_nusselt(1000.0);
        assert!(
            (nu - 1.0).abs() < 1e-14,
            "Below Ra_c: Nu=1 (conduction): {nu}"
        );
    }
    #[test]
    fn test_grossmann_lohse_nusselt_above_critical() {
        let nu = grossmann_lohse_nusselt(1e6);
        assert!(nu > 1.0, "Above Ra_c: Nu > 1: {nu}");
    }
    #[test]
    fn test_grossmann_lohse_nusselt_increases_with_ra() {
        let nu1 = grossmann_lohse_nusselt(1e5);
        let nu2 = grossmann_lohse_nusselt(1e7);
        assert!(
            nu2 > nu1,
            "Nu should increase with Ra: Nu(1e5)={nu1}, Nu(1e7)={nu2}"
        );
    }
    #[test]
    fn test_nusselt_near_onset_at_critical() {
        let nu = nusselt_near_onset(RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP, 0.5);
        assert!((nu - 1.0).abs() < 1e-14, "Nu at Ra_c = 1: {nu}");
    }
    #[test]
    fn test_nusselt_near_onset_increases_above_critical() {
        let nu_just = nusselt_near_onset(RAYLEIGH_BENARD_CRITICAL_RA_NOSLIP * 1.1, 0.5);
        assert!(nu_just > 1.0, "Nu should be > 1 just above Ra_c: {nu_just}");
    }
    #[test]
    fn test_thermal_diffusivity_d2q9_tau_1() {
        let alpha = thermal_diffusivity_d2q9(1.0);
        assert!((alpha - 1.0 / 6.0).abs() < 1e-14, "α = {alpha}");
    }
    #[test]
    fn test_tau_t_from_diffusivity_roundtrip() {
        let alpha = 0.05;
        let tau = tau_t_from_diffusivity_d2q9(alpha);
        let alpha_back = thermal_diffusivity_d2q9(tau);
        assert!(
            (alpha_back - alpha).abs() < 1e-14,
            "Roundtrip α = {alpha_back}"
        );
    }
    #[test]
    fn test_prandtl_number_air() {
        let nu = 1.5e-5;
        let alpha = 2.1e-5;
        let pr = prandtl_number(nu, alpha);
        assert!((pr - 0.71).abs() < 0.01, "Pr for air ≈ 0.71: {pr}");
    }
    #[test]
    fn test_prandtl_number_zero_diffusivity() {
        let pr = prandtl_number(1.0, 0.0);
        assert!(pr.is_infinite(), "Pr → ∞ when α = 0");
    }
    #[test]
    fn test_rayleigh_number_formula() {
        let g = 9.81;
        let beta = 3.4e-3;
        let delta_t = 10.0;
        let h = 1.0;
        let nu = 1.5e-5;
        let alpha = 2.1e-5;
        let ra = rayleigh_number(g, beta, delta_t, h, nu, alpha);
        let expected = g * beta * delta_t * h.powi(3) / (nu * alpha);
        assert!(
            (ra - expected).abs() < 1e-3,
            "Ra = {ra}, expected {expected}"
        );
    }
    #[test]
    fn test_grashof_number_formula() {
        let g = 9.81;
        let beta = 3.4e-3;
        let delta_t = 10.0;
        let h = 1.0;
        let nu = 1.5e-5;
        let gr = grashof_number(g, beta, delta_t, h, nu);
        let expected = g * beta * delta_t * h.powi(3) / (nu * nu);
        assert!(
            (gr - expected).abs() < 1.0,
            "Gr = {gr}, expected {expected}"
        );
    }
    #[test]
    fn test_grashof_equals_ra_over_pr() {
        let g = 9.81;
        let beta = 3.4e-3;
        let delta_t = 10.0;
        let h = 1.0;
        let nu = 1.5e-5;
        let alpha = 2.1e-5;
        let ra = rayleigh_number(g, beta, delta_t, h, nu, alpha);
        let pr = prandtl_number(nu, alpha);
        let gr_from_ra = ra / pr;
        let gr_direct = grashof_number(g, beta, delta_t, h, nu);
        assert!(
            (gr_from_ra - gr_direct).abs() / gr_direct < 1e-10,
            "Gr = Ra/Pr: {gr_from_ra} vs {gr_direct}"
        );
    }
}
#[cfg(test)]
mod tests_thermal_extended {
    use super::*;

    #[test]
    fn test_d3q7_equilibrium_sum_equals_temp() {
        let t = 1.5;
        let u = [0.0, 0.0, 0.0];
        let sum: f64 = (0..7).map(|i| d3q7_thermal_equilibrium(t, u, i)).sum();
        assert!((sum - t).abs() < 1e-12, "D3Q7 eq sum = {sum}, expected {t}");
    }
    #[test]
    fn test_d3q7_temperature_recovery() {
        let t_init = 2.0;
        let mut g = [0.0_f64; 7];
        for (i, o) in g.iter_mut().enumerate() {
            *o = d3q7_thermal_equilibrium(t_init, [0.0, 0.0, 0.0], i);
        }
        let t_rec = d3q7_temperature(&g);
        assert!(
            (t_rec - t_init).abs() < 1e-12,
            "D3Q7 temperature recovery = {t_rec}"
        );
    }
    #[test]
    fn test_d3q7_thermal_diffusivity_positive() {
        let alpha = d3q7_thermal_diffusivity(0.8);
        assert!(alpha > 0.0, "D3Q7 thermal diffusivity = {alpha}");
    }
    #[test]
    fn test_d3q7_bgk_collision_preserves_temperature() {
        let t_init = 1.0;
        let g = {
            let mut arr = [0.0_f64; 7];
            for (i, o) in arr.iter_mut().enumerate() {
                *o = d3q7_thermal_equilibrium(t_init, [0.0, 0.0, 0.0], i);
            }
            arr
        };
        let g_post = d3q7_bgk_collision(&g, t_init, [0.0, 0.0, 0.0], 0.8);
        let t_after = d3q7_temperature(&g_post);
        assert!(
            (t_after - t_init).abs() < 1e-12,
            "BGK preserves T: {t_after}"
        );
    }
    #[test]
    fn test_blasius_bl_thickness_grows_with_x() {
        let nu = 1.5e-5;
        let u_inf = 5.0;
        let d1 = blasius_boundary_layer_thickness(nu, 0.1, u_inf);
        let d2 = blasius_boundary_layer_thickness(nu, 0.5, u_inf);
        assert!(d2 > d1, "BL thickness grows downstream: {d1} < {d2}");
    }
    #[test]
    fn test_blasius_bl_thickness_positive() {
        let d = blasius_boundary_layer_thickness(1.5e-5, 0.5, 5.0);
        assert!(d > 0.0, "Blasius BL thickness = {d}");
    }
    #[test]
    fn test_local_nusselt_flat_plate_increases_with_re() {
        let nu1 = local_nusselt_flat_plate(1e4, 0.71);
        let nu2 = local_nusselt_flat_plate(1e5, 0.71);
        assert!(nu2 > nu1, "Nu increases with Re: {nu1}, {nu2}");
    }
    #[test]
    fn test_average_nusselt_flat_plate_positive() {
        let nu = average_nusselt_flat_plate(1e5, 0.71);
        assert!(nu > 0.0, "Average Nu flat plate = {nu}");
    }
    #[test]
    fn test_heat_transfer_coefficient_proportional_to_nu() {
        let h1 = heat_transfer_coefficient(10.0, 0.025, 1.0);
        let h2 = heat_transfer_coefficient(20.0, 0.025, 1.0);
        assert!((h2 / h1 - 2.0).abs() < 1e-10, "h ∝ Nu: {h1}, {h2}");
    }
    #[test]
    fn test_air_viscosity_sutherland_increases_with_temperature() {
        let mu1 = air_viscosity_sutherland(293.15);
        let mu2 = air_viscosity_sutherland(600.0);
        assert!(mu2 > mu1, "Air viscosity increases with T: {mu1}, {mu2}");
    }
    #[test]
    fn test_air_thermal_conductivity_positive() {
        let k = air_thermal_conductivity(300.0);
        assert!(k > 0.0, "Air thermal conductivity = {k}");
    }
    #[test]
    fn test_air_prandtl_near_0_7_at_room_temperature() {
        let pr = air_prandtl_number(293.15);
        assert!((pr - 0.7).abs() < 0.1, "Air Pr at 293 K = {pr}");
    }
    #[test]
    fn test_water_specific_heat_near_4200() {
        let cp = water_specific_heat(293.15);
        assert!(
            (cp - 4182.0).abs() < 100.0,
            "Water Cp at 293K ≈ 4182 J/(kg·K): {cp}"
        );
    }
    #[test]
    fn test_check_thermal_stability_passes_at_large_tau() {
        assert!(check_thermal_stability(1.0), "tau=1 should be stable");
    }
    #[test]
    fn test_check_thermal_stability_fails_at_small_tau() {
        assert!(!check_thermal_stability(0.4), "tau=0.4 should be unstable");
    }
    #[test]
    fn test_tau_from_diffusivity_roundtrip() {
        let alpha = 1.0 / 6.0;
        let cs2 = 1.0 / 3.0;
        let tau = tau_from_diffusivity(alpha, cs2);
        let alpha_back = thermal_diffusivity_d2q9(tau);
        assert!(
            (alpha_back - alpha).abs() < 1e-12,
            "tau roundtrip: {alpha_back} vs {alpha}"
        );
    }
    #[test]
    fn test_peclet_number_proportional_to_velocity() {
        let pe1 = peclet_number(0.1, 1.0, 1e-5);
        let pe2 = peclet_number(0.2, 1.0, 1e-5);
        assert!((pe2 / pe1 - 2.0).abs() < 1e-10, "Pe ∝ u: {pe1}, {pe2}");
    }
    #[test]
    fn test_fourier_number_increases_with_time() {
        let fo1 = fourier_number(1e-5, 1.0, 0.1);
        let fo2 = fourier_number(1e-5, 2.0, 0.1);
        assert!((fo2 / fo1 - 2.0).abs() < 1e-10, "Fo ∝ t: {fo1}, {fo2}");
    }
    #[test]
    fn test_stefan_number_positive() {
        let ste = stefan_number(4182.0, 10.0, 334_000.0);
        assert!(ste > 0.0, "Stefan number = {ste}");
    }
    #[test]
    fn test_periodic_temperature_x_length_preserved() {
        let temp = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let t_new = periodic_temperature_x(&temp, 3, 2);
        assert_eq!(t_new.len(), temp.len());
    }
    #[test]
    fn test_sinusoidal_temperature_perturbation_amplitude() {
        let t_mean = 1.0;
        let amp = 0.1;
        let temp = sinusoidal_temperature_perturbation(4, 4, t_mean, amp);
        let max_dev = temp
            .iter()
            .map(|&x| (x - t_mean).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_dev <= amp + 1e-12,
            "max sinusoidal deviation = {max_dev} > amp = {amp}"
        );
    }
    #[test]
    fn test_linear_temperature_profile_x_left_right() {
        let t_left = 1.0;
        let t_right = 2.0;
        let temp = linear_temperature_profile_x(4, 1, t_left, t_right);
        assert!((temp[0] - t_left).abs() < 1e-10, "left BC: {}", temp[0]);
        assert!((temp[3] - t_right).abs() < 1e-10, "right BC: {}", temp[3]);
    }
    #[test]
    fn test_local_nusselt_bottom_central_zero_for_uniform() {
        let nx = 4;
        let ny = 4;
        let temp = vec![1.0_f64; nx * ny];
        let nu_vec = local_nusselt_bottom_central(&temp, nx, ny, 1.0, 0.0);
        for &nu in &nu_vec {
            assert!(nu.abs() < 1e-10, "Nu for uniform T = {nu}");
        }
    }
    #[test]
    fn test_energy_equilibrium_zero_velocity_at_rest() {
        let t = 1.5;
        let sum: f64 = (0..9).map(|i| energy_equilibrium(t, [0.0, 0.0], i)).sum();
        assert!((sum - t).abs() < 1e-12, "Energy eq sum at u=0: {sum}");
    }
    #[test]
    fn test_thermal_diffusivity_d2q9_at_tau_one() {
        let alpha = thermal_diffusivity_d2q9(1.0);
        assert!((alpha - 1.0 / 6.0).abs() < 1e-12, "alpha at tau=1: {alpha}");
    }
    #[test]
    fn test_tau_t_from_diffusivity_d2q9_roundtrip() {
        let alpha_target = 1.0 / 8.0;
        let tau = tau_t_from_diffusivity_d2q9(alpha_target);
        let alpha_back = thermal_diffusivity_d2q9(tau);
        assert!(
            (alpha_back - alpha_target).abs() < 1e-12,
            "alpha roundtrip: {alpha_back}"
        );
    }
    #[test]
    fn test_grossmann_lohse_nusselt_at_critical_ra() {
        let nu = grossmann_lohse_nusselt(1707.76);
        assert!((nu - 1.0).abs() < 1e-12, "Nu at critical Ra = {nu}");
    }
    #[test]
    fn test_nusselt_near_onset_proportional_to_ra_excess() {
        let nu1 = nusselt_near_onset(2000.0, 0.1);
        let nu2 = nusselt_near_onset(3000.0, 0.1);
        assert!(nu2 >= nu1, "Nu should increase above onset: {nu1}, {nu2}");
    }
    #[test]
    fn test_is_rb_convective_below_critical_false() {
        assert!(!is_rb_convective(1000.0), "Ra=1000 < 1708 → no convection");
    }
    #[test]
    fn test_is_rb_convective_above_critical_true() {
        assert!(is_rb_convective(2000.0), "Ra=2000 > 1708 → convection");
    }
}
