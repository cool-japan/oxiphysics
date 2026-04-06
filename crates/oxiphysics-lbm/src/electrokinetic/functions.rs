//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains, clippy::needless_range_loop)]
/// Boltzmann constant (J/K)
pub const K_B: f64 = 1.380649e-23;
/// Elementary charge (C)
pub const E_CHARGE: f64 = 1.602176634e-19;
/// Permittivity of free space (F/m)
pub const EPSILON_0: f64 = 8.854187817e-12;
/// Avogadro's number (mol^-1)
pub(super) const N_A: f64 = 6.022e23;
#[cfg(test)]
mod tests {

    use crate::electrokinetic::types::*;
    fn water_params() -> ElectrolyteParams {
        ElectrolyteParams::new(298.15, 80.0, 1e-3)
    }
    #[test]
    fn test_ion_mobility_positive_valence() {
        let na = IonSpecies::new("Na+", 1, 1.33e-9, 100.0);
        let mu = na.mobility(298.15);
        assert!(mu > 0.0, "Na+ mobility = {mu}");
    }
    #[test]
    fn test_debye_length_decreases_with_ionic_strength() {
        let params = water_params();
        let lambda_low = params.debye_length(1.0);
        let lambda_high = params.debye_length(100.0);
        assert!(
            lambda_low > lambda_high,
            "low={lambda_low}, high={lambda_high}"
        );
    }
    #[test]
    fn test_helmholtz_smoluchowski_sign() {
        let params = water_params();
        let zeta = -0.05;
        let eof = ElectroosmosticFlow::new(zeta, params);
        let u = eof.helmholtz_smoluchowski_velocity(1000.0);
        assert!(u > 0.0, "u = {u}");
        assert!(u * (zeta * 1000.0) < 0.0);
    }
    #[test]
    fn test_poisson_solver_converges() {
        let nx = 10;
        let ny = 10;
        let rho_charge = vec![1e-3_f64; nx * ny];
        let mut phi = vec![0.0_f64; nx * ny];
        let iters =
            PoissonSolver::solve_poisson_2d(&mut phi, &rho_charge, nx, ny, 1e-6, 80.0, 5000, 1e-10);
        assert!(iters < 5000, "Did not converge (used {iters})");
        assert!(phi[5 * nx + 5].abs() > 0.0);
    }
    #[test]
    fn test_charge_density_cancels() {
        let conc = 10.0_f64;
        let na = IonSpecies::new("Na+", 1, 1.33e-9, conc);
        let cl = IonSpecies::new("Cl-", -1, 2.03e-9, conc);
        let concentrations = vec![vec![conc; 4], vec![conc; 4]];
        let rho = PoissonSolver::charge_density(&concentrations, &[na, cl]);
        for &r in &rho {
            assert!(r.abs() < 1e-10, "rho = {r}");
        }
    }
    #[test]
    fn test_electrophoretic_velocity_direction() {
        let params = water_params();
        let ep = ElectrophoreticMobility::new(500e-9, -0.05, params);
        let v = ep.velocity(1000.0, 1e6);
        assert!(v < 0.0, "v = {v}");
    }
    #[test]
    fn test_pb_flat_surface_at_wall() {
        let phi = PoissonBoltzmann::flat_surface_potential(-0.05, 0.0, 10e-9);
        assert!((phi - (-0.05)).abs() < 1e-12, "phi(0) = {phi}");
    }
    #[test]
    fn test_pb_flat_surface_decays() {
        let phi_near = PoissonBoltzmann::flat_surface_potential(-0.05, 5e-9, 10e-9);
        let phi_far = PoissonBoltzmann::flat_surface_potential(-0.05, 50e-9, 10e-9);
        assert!(
            phi_near.abs() > phi_far.abs(),
            "near={phi_near}, far={phi_far}"
        );
    }
    #[test]
    fn test_pb_ionic_strength() {
        let na = IonSpecies::new("Na+", 1, 1.33e-9, 100.0);
        let cl = IonSpecies::new("Cl-", -1, 2.03e-9, 100.0);
        let i = PoissonBoltzmann::ionic_strength(&[na, cl]);
        assert!((i - 100.0).abs() < 1e-10, "I = {i}");
    }
    #[test]
    fn test_pb_solve_1d_boundary() {
        let zeta = -0.025;
        let phi = PoissonBoltzmann::solve_1d(zeta, 10e-9, 100, 1e-9);
        assert!((phi[0] - zeta).abs() < 1e-12, "phi[0] = {}", phi[0]);
        assert!(
            phi[99].abs() < 1e-6,
            "phi at far end should be ~0: {}",
            phi[99]
        );
    }
    #[test]
    fn test_pb_solve_1d_monotone() {
        let phi = PoissonBoltzmann::solve_1d(-0.05, 10e-9, 50, 1e-9);
        for i in 0..49 {
            assert!(
                phi[i].abs() >= phi[i + 1].abs() - 1e-14,
                "Not monotone at {i}"
            );
        }
    }
    #[test]
    fn test_debye_layer_thin_edl() {
        assert!(DebyeLayerDiagnostics::is_thin_edl(1e-9, 100e-9));
        assert!(!DebyeLayerDiagnostics::is_thin_edl(10e-9, 20e-9));
    }
    #[test]
    fn test_overlap_parameter() {
        let op = DebyeLayerDiagnostics::overlap_parameter(10e-9, 100e-9);
        assert!((op - 0.1).abs() < 1e-12, "op = {op}");
    }
    #[test]
    fn test_surface_charge_from_zeta() {
        let sigma = DebyeLayerDiagnostics::surface_charge_from_zeta(-0.025, 10e-9, 80.0);
        assert!(
            sigma < 0.0,
            "sigma should be negative for negative zeta: {sigma}"
        );
    }
    #[test]
    fn test_streaming_potential_coefficient_sign() {
        let coeff = StreamingPotential::streaming_potential_coefficient(-0.05, 80.0, 1e-3, 0.1);
        assert!(coeff > 0.0, "coeff = {coeff}");
    }
    #[test]
    fn test_streaming_current_sign() {
        let i_s = StreamingPotential::streaming_current(-0.05, 80.0, 1e-3, 1000.0, 50e-6, 0.01);
        assert!(i_s > 0.0, "I_stream = {i_s}");
    }
    #[test]
    fn test_coupling_coefficient() {
        let l12 = StreamingPotential::coupling_coefficient(-0.05, 80.0, 1e-3);
        assert!(l12 < 0.0, "L12 should be negative for negative zeta: {l12}");
    }
    #[test]
    fn test_slit_channel_velocity_at_center() {
        let params = water_params();
        let eof = ElectroosmosticFlow::new(-0.025, params);
        let u_center = eof.slit_channel_velocity(1000.0, 0.0, 50e-6, 10e-9);
        let u_eo = eof.helmholtz_smoluchowski_velocity(1000.0);
        assert!(
            (u_center / u_eo - 1.0).abs() < 0.01,
            "u_center = {u_center}, u_eo = {u_eo}"
        );
    }
    #[test]
    fn test_slit_channel_velocity_at_wall_zero() {
        let params = water_params();
        let eof = ElectroosmosticFlow::new(-0.025, params);
        let h = 100e-9;
        let lambda_d = 50e-9;
        let u_wall = eof.slit_channel_velocity(1000.0, h, h, lambda_d);
        assert!(u_wall.abs() < 1e-10, "u at wall = {u_wall}");
    }
    #[test]
    fn test_slit_channel_average_less_than_hs() {
        let params = water_params();
        let eof = ElectroosmosticFlow::new(-0.025, params);
        let u_avg = eof.slit_channel_average_velocity(1000.0, 50e-6, 10e-9);
        let u_eo = eof.helmholtz_smoluchowski_velocity(1000.0);
        assert!(u_avg < u_eo, "Average velocity should be <= HS velocity");
        assert!(u_avg > 0.0, "Average should be positive");
    }
    #[test]
    fn test_stokes_drag_proportional() {
        let params = water_params();
        let ep = ElectrophoreticMobility::new(500e-9, -0.025, params);
        let f1 = ep.stokes_drag(1.0);
        let f2 = ep.stokes_drag(2.0);
        assert!(
            (f2 - 2.0 * f1).abs() < 1e-20,
            "Drag should be linear in velocity"
        );
    }
    #[test]
    fn test_sedimentation_velocity_positive() {
        let params = water_params();
        let ep = ElectrophoreticMobility::new(1e-6, -0.025, params);
        let v = ep.sedimentation_velocity(2000.0, 1000.0, 9.81);
        assert!(v > 0.0, "Heavier particle should sediment: v = {v}");
    }
    #[test]
    fn test_peclet_number() {
        let params = water_params();
        let ep = ElectrophoreticMobility::new(500e-9, -0.025, params);
        let pe = ep.peclet_number(1e-4, 1e-9);
        assert!((pe - 0.05).abs() < 1e-8, "Pe = {pe}");
    }
}
/// Helmholtz-Smoluchowski streaming potential.
///
/// E_stream = (epsilon * zeta / (mu * sigma_bulk)) * dP
///
/// # Arguments
/// * `d_p`        - pressure difference (Pa)
/// * `epsilon`    - absolute permittivity (F/m)
/// * `zeta`       - zeta potential (V)
/// * `mu`         - dynamic viscosity (Pa·s)
/// * `sigma_bulk` - bulk electrical conductivity (S/m)
pub fn streaming_potential(d_p: f64, epsilon: f64, zeta: f64, mu: f64, sigma_bulk: f64) -> f64 {
    epsilon * zeta * d_p / (mu * sigma_bulk)
}
/// Linearised Poisson-Boltzmann profile.
///
/// psi(x) = psi_0 * exp(-kappa * x)
///
/// # Arguments
/// * `psi_0`  - surface potential (V)
/// * `kappa`  - inverse Debye length (m^-1)
/// * `x`      - positions (m)
pub fn poisson_boltzmann_1d(psi_0: f64, kappa: f64, x: &[f64]) -> Vec<f64> {
    x.iter().map(|&xi| psi_0 * (-kappa * xi).exp()).collect()
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::electrokinetic::types::*;
    fn water_params() -> ElectrolyteParams {
        ElectrolyteParams::new(298.15, 80.0, 1e-3)
    }
    #[test]
    fn test_edl_energy_decays_with_separation() {
        let e1 = ElectricDoubleLayerEnergy::derjaguin_energy_per_area(
            100.0e3, 298.15, 10e-9, -0.025, 1, 5e-9,
        );
        let e2 = ElectricDoubleLayerEnergy::derjaguin_energy_per_area(
            100.0e3, 298.15, 10e-9, -0.025, 1, 20e-9,
        );
        assert!(e2 < e1, "EDL energy should decay with separation");
    }
    #[test]
    fn test_disjoining_pressure_positive() {
        let p =
            ElectricDoubleLayerEnergy::disjoining_pressure(100.0e3, 298.15, 10e-9, -0.025, 1, 5e-9);
        assert!(p > 0.0, "Disjoining pressure should be positive: {p}");
    }
    #[test]
    fn test_eof_body_force_direction() {
        let eof_force = ElectroosmoticBodyForce::new([1000.0, 0.0, 0.0]);
        let f = eof_force.force_at_cell(1e-3);
        assert!(f[0] > 0.0, "Force should be in E-field direction");
        assert!(f[1].abs() < 1e-15);
    }
    #[test]
    fn test_eof_body_force_apply_2d() {
        let eof_force = ElectroosmoticBodyForce::new([1000.0, 0.0, 0.0]);
        let rho_e = vec![1e-3, -1e-3, 0.0];
        let mut forces = vec![[0.0_f64; 3]; 3];
        eof_force.apply_2d(&rho_e, &mut forces);
        assert!(forces[0][0] > 0.0);
        assert!(forces[1][0] < 0.0);
        assert!(forces[2][0].abs() < 1e-30);
    }
    #[test]
    fn test_zeta_from_streaming_potential() {
        let _params = water_params();
        let zeta_input = -0.025;
        let coeff =
            StreamingPotential::streaming_potential_coefficient(zeta_input, 80.0, 1e-3, 0.1);
        let zeta_est =
            ZetaPotentialEstimator::from_streaming_potential_coefficient(coeff, 1e-3, 0.1, 80.0);
        assert!(
            (zeta_est - zeta_input).abs() < 1e-12,
            "zeta_est = {zeta_est}"
        );
    }
    #[test]
    fn test_zeta_from_mobility() {
        let zeta_input = -0.025_f64;
        let ep = ElectrophoreticMobility::new(500e-9, zeta_input, water_params());
        let kappa = 1.0 / 10e-9;
        let mu_e = ep.mobility(kappa);
        let zeta_est = ZetaPotentialEstimator::from_electrophoretic_mobility(mu_e, 1e-3, 80.0);
        assert!(
            zeta_est < 0.0,
            "zeta should be negative for negative mobility"
        );
    }
    #[test]
    fn test_zeta_from_surface_charge() {
        let sigma = -0.01_f64;
        let zeta = ZetaPotentialEstimator::from_surface_charge(sigma, 10e-9, 80.0);
        assert!(zeta < 0.0, "negative sigma should give negative zeta");
    }
    #[test]
    fn test_transient_eof_startup() {
        let nu = 1e-6_f64;
        let h = 50e-6_f64;
        let u_eo = 1e-4_f64;
        let t = 0.001;
        let u = TransientEof::startup_velocity_center(u_eo, t, nu, h);
        assert!(
            u > 0.0 && u < u_eo,
            "startup velocity should be between 0 and u_eo: {u}"
        );
    }
    #[test]
    fn test_transient_eof_dimensionless_time() {
        let t_star = TransientEof::dimensionless_time(0.01, 1e-6, 50e-6);
        assert!(t_star > 0.0);
    }
    #[test]
    fn test_transient_eof_time_to_steady_state() {
        let t99 = TransientEof::time_to_steady_state(1e-6, 50e-6);
        assert!(t99 > 0.0, "t99 = {t99}");
        let expected = 0.466 * (50e-6 * 50e-6) / 1e-6;
        assert!(
            (t99 - expected).abs() / expected < 0.05,
            "t99 = {t99}, expected ~ {expected}"
        );
    }
    #[test]
    fn test_diffusio_osmosis_velocity_sign() {
        let u = DiffusioOsmosis::velocity(-0.025, 298.15, 1e-3, 80.0, 1.0);
        assert!(
            u < 0.0,
            "Negative zeta + positive grad_ln_c → negative velocity: {u}"
        );
    }
    #[test]
    fn test_lewis_number() {
        let le = DiffusioOsmosis::lewis_number(1.4e-7, 1e-9);
        assert!(le > 100.0, "Le = {le}");
    }
    #[test]
    fn test_activity_coefficient_less_than_one() {
        let gamma = DebyeHuckel::activity_coefficient(1, 100.0, 3e-10);
        assert!(gamma < 1.0 && gamma > 0.0, "gamma = {gamma}");
    }
    #[test]
    fn test_activity_coefficient_approaches_one_at_zero_strength() {
        let gamma = DebyeHuckel::activity_coefficient(1, 1e-6, 3e-10);
        assert!(
            (gamma - 1.0).abs() < 0.01,
            "gamma at low I should be ~1, got {gamma}"
        );
    }
    #[test]
    fn test_activity_coefficient_higher_valence_lower() {
        let g1 = DebyeHuckel::activity_coefficient(1, 10.0, 3e-10);
        let g2 = DebyeHuckel::activity_coefficient(2, 10.0, 3e-10);
        assert!(
            g2 < g1,
            "divalent ion should have smaller activity coefficient"
        );
    }
    #[test]
    fn test_diffusive_flux_uniform_concentration() {
        let c = vec![1.0; 16];
        let flux = NernstPlanckSolver::diffusive_flux(&c, 4, 4, 0.1, 1e-9);
        for &f in &flux {
            assert!(f.abs() < 1e-20, "uniform concentration flux = {f}");
        }
    }
    #[test]
    fn test_migrative_flux_uniform_potential() {
        let c = vec![1.0; 9];
        let phi = vec![0.5; 9];
        let flux = NernstPlanckSolver::migrative_flux(&c, &phi, 3, 3, 0.1, 1e-8);
        for &f in &flux {
            assert!(f.abs() < 1e-25, "uniform potential flux = {f}");
        }
    }
    #[test]
    fn test_update_concentrations_decreasing_gradient() {
        let nx = 5;
        let ny = 5;
        let n = nx * ny;
        let mut c = vec![1.0; n];
        c[12] = 2.0;
        let flux = NernstPlanckSolver::diffusive_flux(&c, nx, ny, 0.1, 1e-9);
        let c_before = c[12];
        NernstPlanckSolver::update_concentrations(&mut c, &flux, nx, ny, 0.001, 0.1);
        assert!(
            c[12] < c_before,
            "Peak should decrease after diffusion step"
        );
    }
    #[test]
    fn test_charge_density_nonzero() {
        let mg2 = IonSpecies::new("Mg2+", 2, 0.7e-9, 10.0);
        let concentrations = vec![vec![10.0; 4]];
        let rho = PoissonSolver::charge_density(&concentrations, &[mg2]);
        for &r in &rho {
            assert!(
                r > 0.0,
                "divalent cation should give positive charge density: {r}"
            );
        }
    }
    #[test]
    fn test_debye_length_from_ionic_decreases_with_concentration() {
        let eps = EPSILON_0 * 80.0;
        let ld1 = DoubleLayer::debye_length_from_ionic(1.0, 1.0, 298.15, eps);
        let ld2 = DoubleLayer::debye_length_from_ionic(100.0, 1.0, 298.15, eps);
        assert!(
            ld1 > ld2,
            "Debye length should decrease with ionic concentration: {ld1} > {ld2}"
        );
    }
    #[test]
    fn test_double_layer_surface_charge_sign() {
        let eps = EPSILON_0 * 80.0;
        let dl = DoubleLayer::new(10e-9, -0.025, eps);
        let sigma = dl.surface_charge_density();
        assert!(
            sigma > 0.0,
            "sigma should be positive for negative zeta: {sigma}"
        );
    }
    #[test]
    fn test_eo_velocity_direction_matches_field() {
        let eo = ElectroosmoticVelocity::new(1e-8);
        let e = [1000.0, 0.0, 0.0];
        let v = eo.velocity(e);
        assert!(
            v[0] > 0.0,
            "EO velocity should align with positive field: {}",
            v[0]
        );
        assert!(v[1].abs() < 1e-30);
        assert!(v[2].abs() < 1e-30);
    }
    #[test]
    fn test_eo_velocity_opposes_negative_field() {
        let eo = ElectroosmoticVelocity::new(1e-8);
        let e = [-500.0, 0.0, 0.0];
        let v = eo.velocity(e);
        assert!(
            v[0] < 0.0,
            "EO velocity should follow field direction: {}",
            v[0]
        );
    }
    #[test]
    fn test_streaming_potential_sign() {
        let eps = EPSILON_0 * 80.0;
        let sp = streaming_potential(1000.0, eps, 0.025, 1e-3, 0.1);
        assert!(sp > 0.0, "streaming potential should be positive: {sp}");
    }
    #[test]
    fn test_streaming_potential_negative_zeta() {
        let eps = EPSILON_0 * 80.0;
        let sp = streaming_potential(1000.0, eps, -0.025, 1e-3, 0.1);
        assert!(
            sp < 0.0,
            "negative zeta → negative streaming potential: {sp}"
        );
    }
    #[test]
    fn test_linearized_pb_at_surface() {
        let psi_0 = -0.025;
        let kappa = 1.0 / 10e-9;
        let x = [0.0];
        let psi = poisson_boltzmann_1d(psi_0, kappa, &x);
        assert!(
            (psi[0] - psi_0).abs() < 1e-15,
            "psi(0) should equal psi_0: {}",
            psi[0]
        );
    }
    #[test]
    fn test_linearized_pb_exponential_decay() {
        let psi_0 = -0.025;
        let kappa = 1.0 / 10e-9;
        let x = [0.0, 10e-9, 20e-9, 50e-9];
        let psi = poisson_boltzmann_1d(psi_0, kappa, &x);
        for i in 0..x.len() - 1 {
            assert!(
                psi[i].abs() > psi[i + 1].abs(),
                "PB profile should decay: psi[{}]={}, psi[{}]={}",
                i,
                psi[i],
                i + 1,
                psi[i + 1]
            );
        }
    }
}
/// D2Q9 lattice constants for electrokinetic LBM.
pub(super) const EK_CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
pub(super) const EK_CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
pub(super) const EK_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
pub(super) const EK_CS2: f64 = 1.0 / 3.0;
#[cfg(test)]
mod tests_electrokinetic_extended {
    use super::*;
    use crate::electrokinetic::types::*;
    #[test]
    fn test_nonlinear_pb_boundary_condition() {
        let zeta = -0.025;
        let n_bulk = 100.0;
        let z = 1.0;
        let eps = EPSILON_0 * 80.0;
        let t = 298.15;
        let phi = NonlinearPoissonBoltzmann::solve_1d(zeta, n_bulk, z, eps, t, 50, 1e-9);
        assert!(
            (phi[0] - zeta).abs() < 1e-12,
            "Nonlinear PB: phi[0] should equal zeta: {}",
            phi[0]
        );
        assert!(
            phi[49].abs() < 1e-4,
            "Far-field phi should approach 0: {}",
            phi[49]
        );
    }
    #[test]
    fn test_nonlinear_pb_monotone_decay() {
        let zeta = -0.025;
        let phi = NonlinearPoissonBoltzmann::solve_1d(
            zeta,
            100.0,
            1.0,
            EPSILON_0 * 80.0,
            298.15,
            40,
            1e-9,
        );
        for i in 0..39 {
            assert!(
                phi[i].abs() >= phi[i + 1].abs() - 1e-12,
                "Nonlinear PB not monotone at i={i}"
            );
        }
    }
    #[test]
    fn test_nonlinear_pb_charge_density_sign() {
        let zeta = -0.05;
        let phi = NonlinearPoissonBoltzmann::solve_1d(
            zeta,
            100.0,
            1.0,
            EPSILON_0 * 80.0,
            298.15,
            30,
            1e-9,
        );
        let rho = NonlinearPoissonBoltzmann::charge_density_field(&phi, 100.0, 1.0, 298.15);
        assert!(
            rho[0] > 0.0,
            "Positive charge near negative surface: {}",
            rho[0]
        );
    }
    #[test]
    fn test_onsager_symmetry() {
        let mat =
            OnsagerCouplingMatrix::from_physical(1e-12, 1e-3, -0.025, EPSILON_0 * 80.0, 0.1, 10e-9);
        assert!(mat.is_symmetric(), "Onsager matrix should be symmetric");
    }
    #[test]
    fn test_onsager_figure_of_merit_positive() {
        let mat =
            OnsagerCouplingMatrix::from_physical(1e-12, 1e-3, -0.025, EPSILON_0 * 80.0, 0.1, 10e-9);
        let z = mat.figure_of_merit();
        assert!(z > 0.0, "Figure of merit should be positive: {z}");
    }
    #[test]
    fn test_onsager_max_efficiency_between_0_and_1() {
        let mat = OnsagerCouplingMatrix {
            l11: 1e-10,
            l12: 1e-8,
            l22: 0.1,
        };
        let eta = mat.max_efficiency();
        assert!(
            eta >= 0.0 && eta <= 1.0,
            "Efficiency should be in [0,1]: {eta}"
        );
    }
    #[test]
    fn test_onsager_volume_flux_zero_field() {
        let mat = OnsagerCouplingMatrix {
            l11: 1e-10,
            l12: 1e-8,
            l22: 0.1,
        };
        let j = mat.volume_flux(0.0, 0.0);
        assert!(j.abs() < 1e-30, "Zero flux at zero driving: {j}");
    }
    #[test]
    fn test_eo_pump_max_flow_direction() {
        let pump = ElectroosmosticPump::new(-0.025, EPSILON_0 * 80.0, 1e-3, 1e-12, 1e-6);
        let q = pump.max_flow_rate(1000.0);
        assert!(
            q < 0.0,
            "Negative zeta + positive E should give negative flow: {q}"
        );
    }
    #[test]
    fn test_eo_pump_flow_at_max_backpressure_zero() {
        let pump = ElectroosmosticPump::new(-0.025, EPSILON_0 * 80.0, 1e-3, 1e-12, 1e-6);
        let e = 1000.0;
        let dp_max = pump.max_pressure(e);
        let q = pump.flow_rate(e, dp_max);
        assert!(
            q.abs() < 1e-20,
            "Flow at max back-pressure should be ~0: {q}"
        );
    }
    #[test]
    fn test_eo_pump_max_flow_proportional_to_field() {
        let pump = ElectroosmosticPump::new(-0.025, EPSILON_0 * 80.0, 1e-3, 1e-12, 1e-6);
        let q1 = pump.max_flow_rate(1000.0);
        let q2 = pump.max_flow_rate(2000.0);
        assert!(
            (q2 / q1 - 2.0).abs() < 1e-10,
            "Flow should be proportional to field: {}",
            q2 / q1
        );
    }
    #[test]
    fn test_edl_capacitance_gouy_chapman_positive() {
        let eps = EPSILON_0 * 80.0;
        let kappa = 1.0 / 10e-9;
        let c = ElectricDoubleLayerCapacitance::gouy_chapman(eps, kappa, -0.025, 1.0, 298.15);
        assert!(c > 0.0, "EDL capacitance should be positive: {c}");
    }
    #[test]
    fn test_edl_capacitance_helmholtz() {
        let eps = EPSILON_0 * 80.0;
        let d_h = 0.4e-9;
        let c_h = ElectricDoubleLayerCapacitance::helmholtz_layer(eps, d_h);
        let expected = eps / d_h;
        assert!(
            (c_h - expected).abs() < 1e-10 * expected.abs(),
            "Helmholtz C mismatch"
        );
    }
    #[test]
    fn test_edl_capacitance_series_less_than_components() {
        let c_h: f64 = 0.2;
        let c_d: f64 = 0.5;
        let c_t = ElectricDoubleLayerCapacitance::total_series(c_h, c_d);
        assert!(
            c_t < c_h.min(c_d),
            "Series capacitance should be less than either component"
        );
    }
    #[test]
    fn test_electroviscous_ratio_gte_one() {
        let eps = EPSILON_0 * 80.0;
        let ratio = ElectroviscousEffect::apparent_viscosity_ratio(eps, -0.025, 1e-3, 1e-9);
        assert!(
            ratio >= 1.0,
            "Apparent viscosity ratio should be >= 1: {ratio}"
        );
    }
    #[test]
    fn test_electroviscous_larger_zeta_larger_ratio() {
        let eps = EPSILON_0 * 80.0;
        let r1 = ElectroviscousEffect::apparent_viscosity_ratio(eps, -0.025, 1e-3, 1e-9);
        let r2 = ElectroviscousEffect::apparent_viscosity_ratio(eps, -0.050, 1e-3, 1e-9);
        assert!(
            r2 > r1,
            "Larger |zeta| should give higher apparent viscosity: {r2} vs {r1}"
        );
    }
    #[test]
    fn test_eklbm_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let mut lbm = ElectrokineticLbm::new(nx, ny, 1.0 / 6.0);
        lbm.set_uniform(1.0, [0.0, 0.0]);
        lbm.rho_e = vec![1e-4; nx * ny];
        lbm.set_electric_field(0.01, 0.0);
        let mass_before = lbm.total_mass();
        for _ in 0..10 {
            lbm.step();
        }
        let mass_after = lbm.total_mass();
        assert!(
            (mass_before - mass_after).abs() / mass_before < 1e-6,
            "EK-LBM mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_eklbm_electroosmotic_flow_direction() {
        let nx = 10;
        let ny = 10;
        let mut lbm = ElectrokineticLbm::new(nx, ny, 1.0 / 6.0);
        lbm.set_uniform(1.0, [0.0, 0.0]);
        lbm.rho_e = vec![1e-3; nx * ny];
        lbm.set_electric_field(0.05, 0.0);
        for _ in 0..50 {
            lbm.step();
        }
        let u_avg = lbm.average_velocity();
        assert!(u_avg > 0.0, "EK-LBM should develop positive flow: {u_avg}");
    }
    #[test]
    fn test_eklbm_potential_update() {
        let nx = 5;
        let ny = 5;
        let mut lbm = ElectrokineticLbm::new(nx, ny, 1.0 / 6.0);
        lbm.rho_e = vec![1e-6; nx * ny];
        let eps = EPSILON_0 * 80.0;
        lbm.update_potential_gauss_seidel(eps, 1.0);
        for &p in &lbm.phi {
            assert!(p.is_finite(), "Potential should be finite");
        }
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_Re_proportional_to_velocity() {
        let re1 = ElectrokineticNumbers::Re(1000.0, 1e-4, 1e-3, 1e-3);
        let re2 = ElectrokineticNumbers::Re(1000.0, 2e-4, 1e-3, 1e-3);
        assert!(
            (re2 / re1 - 2.0).abs() < 1e-10,
            "Re should double with velocity"
        );
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_Du_number_positive() {
        let du = ElectrokineticNumbers::Du(1e-8, 0.1, 1e-6);
        assert!(du > 0.0, "Dukhin number should be positive: {du}");
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_kappa_L_scaling() {
        let kl1 = ElectrokineticNumbers::kappa_L(10e-9, 1e-6);
        let kl2 = ElectrokineticNumbers::kappa_L(10e-9, 2e-6);
        assert!(
            (kl2 / kl1 - 2.0).abs() < 1e-10,
            "kappa*L should scale with L"
        );
    }
    #[test]
    fn test_stern_layer_potential_drop_sign() {
        let stern = SternLayer::new(0.4e-9, EPSILON_0 * 5.0, 0.01);
        let drop = stern.potential_drop();
        assert!(
            drop > 0.0,
            "Positive sigma should give positive potential drop: {drop}"
        );
    }
    #[test]
    fn test_stern_layer_capacitance() {
        let stern = SternLayer::new(0.4e-9, EPSILON_0 * 5.0, 0.0);
        let c = stern.capacitance();
        let expected = EPSILON_0 * 5.0 / 0.4e-9;
        assert!(
            (c - expected).abs() < 1e-3 * expected.abs(),
            "Stern capacitance mismatch: {c} vs {expected}"
        );
    }
    #[test]
    fn test_stern_ohp_potential() {
        let stern = SternLayer::new(0.4e-9, EPSILON_0 * 5.0, 0.01);
        let psi_0: f64 = -0.05;
        let psi_ohp = stern.ohp_potential(psi_0);
        assert!(
            psi_ohp.abs() <= psi_0.abs() + stern.potential_drop().abs() + 1e-15,
            "OHP potential out of expected range: {psi_ohp}"
        );
    }
    #[test]
    fn test_boltzmann_concentration_at_zero() {
        let c_bulk = 100.0;
        let c = BoltzmannIonDistribution::concentration(c_bulk, 1.0, 0.0, 298.15);
        assert!(
            (c - c_bulk).abs() < 1e-10,
            "Concentration at zero potential should equal bulk: {c}"
        );
    }
    #[test]
    fn test_boltzmann_cation_enrichment_near_negative_surface() {
        let c_bulk = 100.0;
        let psi_negative = -0.025;
        let c = BoltzmannIonDistribution::concentration(c_bulk, 1.0, psi_negative, 298.15);
        assert!(
            c > c_bulk,
            "Cation concentration should increase near negative surface: {c}"
        );
    }
    #[test]
    fn test_boltzmann_anion_depletion_near_negative_surface() {
        let c_bulk = 100.0;
        let psi_negative = -0.025;
        let c = BoltzmannIonDistribution::concentration(c_bulk, -1.0, psi_negative, 298.15);
        assert!(
            c < c_bulk,
            "Anion concentration should decrease near negative surface: {c}"
        );
    }
    #[test]
    fn test_boltzmann_net_charge_density_near_negative_wall() {
        let psi = -0.025;
        let rho_e = BoltzmannIonDistribution::net_charge_density_symmetric(100.0, psi, 298.15);
        assert!(
            rho_e > 0.0,
            "Net charge near negative wall should be positive: {rho_e}"
        );
    }
    #[test]
    fn test_electric_field_from_potential_linear() {
        let a = 100.0;
        let dx = 0.01;
        let phi: Vec<f64> = (0..5).map(|i| a * i as f64 * dx).collect();
        let e = BoltzmannIonDistribution::electric_field_from_potential(&phi, dx);
        for i in 1..4 {
            assert!(
                (e[i] + a).abs() < 1e-6,
                "E field at {i} should be ~{}, got {}",
                -a,
                e[i]
            );
        }
    }
    #[test]
    fn test_streaming_potential_sign_with_negative_zeta() {
        let sp = streaming_potential(1e3, 7.08e-10, -0.05, 1e-3, 1e-3);
        assert!(
            sp < 0.0,
            "streaming potential should be negative for zeta<0: {sp}"
        );
    }
    #[test]
    fn test_streaming_potential_proportional_to_zeta() {
        let sp1 = streaming_potential(1e3, 7.08e-10, -0.025, 1e-3, 1e-3);
        let sp2 = streaming_potential(1e3, 7.08e-10, -0.050, 1e-3, 1e-3);
        assert!(
            (sp2 / sp1 - 2.0).abs() < 1e-10,
            "streaming potential ∝ zeta: {sp1}, {sp2}"
        );
    }
    #[test]
    fn test_poisson_boltzmann_1d_at_surface() {
        let psi = poisson_boltzmann_1d(-0.025, 1e8, &[0.0]);
        assert!((psi[0] + 0.025).abs() < 1e-14, "PB at surface: {}", psi[0]);
    }
    #[test]
    fn test_poisson_boltzmann_1d_decays_monotonically() {
        let kappa = 1e8;
        let xs: Vec<f64> = (0..5).map(|i| i as f64 * 1e-9).collect();
        let psi = poisson_boltzmann_1d(-0.05, kappa, &xs);
        for i in 1..psi.len() {
            assert!(
                psi[i].abs() <= psi[i - 1].abs(),
                "PB potential should decay: psi[{i}]={} > psi[{}]={}",
                psi[i].abs(),
                i - 1,
                psi[i - 1].abs()
            );
        }
    }
    #[test]
    fn test_electrophoresis_particle_velocity_direction() {
        let p = ElectrophoresisParticle::new(100e-9, 0.05, 1e-3);
        let v = p.electrophoretic_velocity(1000.0);
        assert!(v > 0.0, "positive zeta + positive E → positive v: {v}");
    }
    #[test]
    fn test_electrophoresis_particle_velocity_zero_at_zero_e() {
        let p = ElectrophoresisParticle::new(100e-9, 0.05, 1e-3);
        let v = p.electrophoretic_velocity(0.0);
        assert!(v.abs() < 1e-30, "velocity zero at E=0: {v}");
    }
    #[test]
    fn test_derjaguin_energy_decays_with_separation() {
        let energy1 = ElectricDoubleLayerEnergy::derjaguin_energy_per_area(
            6.022e23, 298.0, 10e-9, -0.025, 1, 1e-9,
        );
        let energy2 = ElectricDoubleLayerEnergy::derjaguin_energy_per_area(
            6.022e23, 298.0, 10e-9, -0.025, 1, 5e-9,
        );
        assert!(
            energy1 > energy2,
            "DLVO energy decays with D: {energy1}, {energy2}"
        );
    }
    #[test]
    fn test_zeta_from_surface_charge_scales_with_lambda_d() {
        let zeta1 = ZetaPotentialEstimator::from_surface_charge(0.01, 10e-9, 80.0);
        let zeta2 = ZetaPotentialEstimator::from_surface_charge(0.01, 20e-9, 80.0);
        assert!(
            (zeta2 / zeta1 - 2.0).abs() < 1e-10,
            "zeta ∝ lambda_d: {zeta1}, {zeta2}"
        );
    }
    #[test]
    fn test_zeta_from_electrophoretic_mobility_positive() {
        let zeta = ZetaPotentialEstimator::from_electrophoretic_mobility(1e-8, 1e-3, 80.0);
        assert!(zeta > 0.0, "zeta from positive mobility: {zeta}");
    }
    #[test]
    fn test_transient_eof_steady_state_approached() {
        let u_eo = 0.01;
        let u = TransientEof::startup_velocity_center(u_eo, 1e10, 1e-6, 1e-4);
        assert!((u - u_eo).abs() < 1e-6, "Transient EOF at large t = {u}");
    }
    #[test]
    fn test_transient_eof_zero_at_t_zero() {
        let u = TransientEof::startup_velocity_center(0.01, 0.0, 1e-6, 1e-4);
        assert!(u.abs() < 1e-14, "Transient EOF at t=0 = {u}");
    }
    #[test]
    fn test_transient_eof_dimensionless_time_positive() {
        let t_star = TransientEof::dimensionless_time(1.0, 1e-6, 1e-4);
        assert!(t_star > 0.0, "Dimensionless time = {t_star}");
    }
    #[test]
    fn test_transient_eof_time_to_steady_state_positive() {
        let t99 = TransientEof::time_to_steady_state(1e-6, 1e-4);
        assert!(t99 > 0.0, "Time to steady state = {t99}");
    }
    #[test]
    fn test_double_layer_debye_length_field() {
        let dl = DoubleLayer::new(10e-9, -0.025, 7.08e-10);
        assert!(
            (dl.debye_length - 10e-9).abs() < 1e-20,
            "Debye length = {}",
            dl.debye_length
        );
    }
    #[test]
    fn test_double_layer_surface_charge_negative_zeta() {
        let dl = DoubleLayer::new(10e-9, -0.025, 7.08e-10);
        let sigma = dl.surface_charge_density();
        assert!(
            sigma > 0.0,
            "Positive surface charge for negative zeta: {sigma}"
        );
    }
    #[test]
    fn test_double_layer_from_ionic_positive() {
        let lambda_d = DoubleLayer::debye_length_from_ionic(0.1, 1.0, 298.0, 7.08e-10);
        assert!(
            lambda_d > 0.0 && lambda_d.is_finite(),
            "Debye length from ionic = {lambda_d}"
        );
    }
    #[test]
    fn test_electroosmotic_velocity_proportional_to_e_field() {
        let eof = ElectroosmoticVelocity::new(-3.54e-8);
        let u1 = eof.velocity([100.0, 0.0, 0.0]);
        let u2 = eof.velocity([200.0, 0.0, 0.0]);
        assert!(
            (u2[0] / u1[0] - 2.0).abs() < 1e-10,
            "EOF ∝ E: u1={}, u2={}",
            u1[0],
            u2[0]
        );
    }
    #[test]
    fn test_electroosmotic_velocity_direction_follows_mu_eo() {
        let eof = ElectroosmoticVelocity::new(-3.54e-8);
        let u = eof.velocity([100.0, 0.0, 0.0]);
        assert!(u[0] < 0.0, "Negative mu_eo → negative velocity: {}", u[0]);
    }
    #[test]
    fn test_eo_body_force_zero_at_zero_charge_density() {
        let bf = ElectroosmoticBodyForce::new([100.0, 0.0, 0.0]);
        let f = bf.force_at_cell(0.0);
        assert!(
            f[0].abs() < 1e-30,
            "Zero charge → zero body force: {}",
            f[0]
        );
    }
    #[test]
    fn test_eo_body_force_direction_follows_e_field() {
        let bf = ElectroosmoticBodyForce::new([100.0, 0.0, 0.0]);
        let f = bf.force_at_cell(1000.0);
        assert!(f[0] > 0.0, "Positive charge+E → positive force: {}", f[0]);
    }
    #[test]
    fn test_nernst_planck_diffusive_flux_zero_uniform() {
        let c = vec![1.0_f64; 16];
        let flux = NernstPlanckSolver::diffusive_flux(&c, 4, 4, 0.01, 1e-9);
        let max_flux = flux.iter().map(|&x| x.abs()).fold(0.0_f64, f64::max);
        assert!(
            max_flux < 1e-20,
            "Uniform concentration → zero diffusive flux: {max_flux}"
        );
    }
    #[test]
    fn test_ion_species_new_preserves_fields() {
        let ion = IonSpecies::new("Na+", 1, 1.33e-9, 0.15);
        assert_eq!(ion.valence, 1);
        assert!((ion.concentration_bulk - 0.15).abs() < 1e-15);
    }
    #[test]
    fn test_ion_species_mobility_scales_with_diffusivity() {
        let ion1 = IonSpecies::new("A", 1, 1e-9, 0.1);
        let ion2 = IonSpecies::new("B", 1, 2e-9, 0.1);
        let m1 = ion1.mobility(298.0);
        let m2 = ion2.mobility(298.0);
        assert!(
            (m2 / m1 - 2.0).abs() < 1e-10,
            "Mobility ∝ diffusivity: {m1}, {m2}"
        );
    }
    #[test]
    fn test_poisson_solver_charge_density_symmetric() {
        let na = IonSpecies::new("Na+", 1, 1.33e-9, 0.1);
        let cl = IonSpecies::new("Cl-", -1, 2.03e-9, 0.1);
        let c_na = vec![0.1_f64; 4];
        let c_cl = vec![0.1_f64; 4];
        let rho = PoissonSolver::charge_density(&[c_na, c_cl], &[na, cl]);
        let max_rho = rho.iter().map(|&x| x.abs()).fold(0.0_f64, f64::max);
        assert!(
            max_rho < 1e-20,
            "Equal Na/Cl → zero charge density: {max_rho}"
        );
    }
    #[test]
    fn test_debye_layer_overlap_at_half_width_is_one() {
        let overlap = DebyeLayerDiagnostics::overlap_parameter(5e-9, 5e-9);
        assert!(
            (overlap - 1.0).abs() < 1e-14,
            "overlap at lambda_d = h: {overlap}"
        );
    }
    #[test]
    fn test_is_thin_edl_for_small_lambda_d() {
        assert!(
            DebyeLayerDiagnostics::is_thin_edl(1e-9, 1e-6),
            "1 nm << 1 um → thin EDL"
        );
    }
    #[test]
    fn test_is_not_thin_edl_for_large_lambda_d() {
        assert!(
            !DebyeLayerDiagnostics::is_thin_edl(1e-6, 1e-6),
            "lambda_d = h → not thin EDL"
        );
    }
    #[test]
    fn test_electrolyte_debye_length_finite() {
        let params = ElectrolyteParams::new(298.0, 80.0, 1e-3);
        let lambda_d = params.debye_length(0.1);
        assert!(
            lambda_d > 0.0 && lambda_d.is_finite(),
            "Debye length = {lambda_d}"
        );
    }
    #[test]
    fn test_ionic_concentration_field_new() {
        let field = IonicConcentrationField::new(4, 4, 2);
        assert_eq!(field.nx, 4);
        assert_eq!(field.ny, 4);
        assert_eq!(field.n_species, 2);
    }
    #[test]
    fn test_ionic_concentration_field_idx_first() {
        let field = IonicConcentrationField::new(5, 5, 2);
        assert_eq!(field.idx(0, 0), 0);
    }
    #[test]
    fn test_ionic_concentration_field_idx_last() {
        let field = IonicConcentrationField::new(5, 5, 2);
        assert_eq!(field.idx(4, 4), 24);
    }
}
