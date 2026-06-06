//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::phase_field::types::*;
    fn make_params(mobility: f64, interface_width: f64, surface_tension: f64) -> PhaseFieldParams {
        PhaseFieldParams {
            mobility,
            interface_width,
            surface_tension,
            kappa: surface_tension * interface_width,
            alpha: 1.0,
            n_components: 2,
        }
    }
    #[test]
    fn test_phase_field_volume_conservation() {
        let nx = 32;
        let ny = 32;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(16.0, 16.0, 6.0);
        let vf_initial = pf.volume_fraction();
        let dt = 0.001;
        for _ in 0..10 {
            pf.step(dt);
        }
        let vf_final = pf.volume_fraction();
        assert!(
            (vf_final - vf_initial).abs() < 0.1 * vf_initial.max(1e-10),
            "Volume fraction not conserved: initial={vf_initial}, final={vf_final}"
        );
    }
    #[test]
    fn test_phase_field_circle_initialization() {
        let nx = 64;
        let ny = 64;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(32.0, 32.0, 10.0);
        let k_center = pf.idx(32, 32);
        assert!(
            pf.phi[k_center] > 0.9,
            "phi at center should be ~1, got {}",
            pf.phi[k_center]
        );
        let k_corner = pf.idx(0, 0);
        assert!(
            pf.phi[k_corner] < 0.1,
            "phi at far corner should be ~0, got {}",
            pf.phi[k_corner]
        );
    }
    #[test]
    fn test_allen_cahn_tanh_interface() {
        let nx = 4;
        let ny = 64;
        let epsilon = 2.0;
        let y_interface = 32.0;
        let dx = 1.0;
        let mut ac = AllenCahn::new(nx, ny, epsilon, 1.0);
        ac.initialize_tanh_interface(y_interface, dx);
        let k = ac.idx(0, 32);
        assert!(
            (ac.phi[k] - 0.5).abs() < 0.05,
            "phi at interface row should be ~0.5, got {}",
            ac.phi[k]
        );
        let k_bottom = ac.idx(0, 0);
        assert!(
            ac.phi[k_bottom] < 0.1,
            "phi far below interface should be ~0, got {}",
            ac.phi[k_bottom]
        );
        let k_top = ac.idx(0, ny - 1);
        assert!(
            ac.phi[k_top] > 0.9,
            "phi far above interface should be ~1, got {}",
            ac.phi[k_top]
        );
    }
    #[test]
    fn test_allen_cahn_step_stability() {
        let nx = 16;
        let ny = 32;
        let epsilon = 2.0;
        let dx = 1.0;
        let mut ac = AllenCahn::new(nx, ny, epsilon, 0.1);
        ac.initialize_tanh_interface(16.0, dx);
        let dt = 0.01;
        for _ in 0..50 {
            ac.step(dt, dx);
        }
        for (k, &phi) in ac.phi.iter().enumerate() {
            assert!(
                (-0.1..=1.1).contains(&phi),
                "phi[{k}] = {phi} out of bounds [-0.1, 1.1]"
            );
        }
    }
    #[test]
    fn test_interface_energy_positive() {
        let nx = 32;
        let ny = 32;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(16.0, 16.0, 8.0);
        let energy = pf.interface_energy();
        assert!(
            energy > 0.0,
            "Interface energy should be positive for non-trivial phi, got {energy}"
        );
    }
    #[test]
    fn test_total_free_energy_positive() {
        let nx = 32;
        let ny = 32;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(16.0, 16.0, 8.0);
        let fe = pf.total_free_energy();
        assert!(fe > 0.0, "Total free energy should be positive: {fe}");
    }
    #[test]
    fn test_gradient_uniform_field() {
        let nx = 16;
        let ny = 16;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        for p in pf.phi.iter_mut() {
            *p = 0.5;
        }
        let grad = pf.compute_gradient();
        for (k, &gk) in grad.iter().enumerate() {
            assert!(
                gk[0].abs() < 1e-14 && gk[1].abs() < 1e-14,
                "Uniform field should have zero gradient at k={k}"
            );
        }
    }
    #[test]
    fn test_laplacian_uniform_field() {
        let nx = 16;
        let ny = 16;
        let params = make_params(0.01, 2.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        for p in pf.phi.iter_mut() {
            *p = 0.7;
        }
        let lap = pf.compute_laplacian();
        for (k, &l) in lap.iter().enumerate() {
            assert!(
                l.abs() < 1e-14,
                "Uniform field: Laplacian[{k}] = {l} should be 0"
            );
        }
    }
    #[test]
    fn test_interface_normal_unit() {
        let nx = 64;
        let ny = 64;
        let params = make_params(0.01, 3.0, 1.0);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(32.0, 32.0, 10.0);
        let normals = pf.compute_interface_normal();
        let grad = pf.compute_gradient();
        for k in 0..normals.len() {
            let mag_grad = (grad[k][0] * grad[k][0] + grad[k][1] * grad[k][1]).sqrt();
            if mag_grad > 1e-6 {
                let mag_n = (normals[k][0] * normals[k][0] + normals[k][1] * normals[k][1]).sqrt();
                assert!(
                    (mag_n - 1.0).abs() < 1e-10,
                    "Interface normal should be unit: |n|={mag_n} at k={k}"
                );
            }
        }
    }
    #[test]
    fn test_constant_mobility() {
        assert!((constant_mobility(0.5, 0.3) - 0.5).abs() < 1e-14);
        assert!((constant_mobility(0.5, 0.9) - 0.5).abs() < 1e-14);
    }
    #[test]
    fn test_degenerate_mobility_bulk() {
        assert!(degenerate_mobility(1.0, 0.0).abs() < 1e-14);
        assert!(degenerate_mobility(1.0, 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_degenerate_mobility_interface() {
        let m = degenerate_mobility(1.0, 0.5);
        assert!((m - 0.25).abs() < 1e-14, "M(0.5) = 0.25, got {m}");
    }
    #[test]
    fn test_body_force_uniform_mu() {
        let nx = 8;
        let ny = 8;
        let n = nx * ny;
        let phi = vec![0.5; n];
        let mu = vec![1.0; n];
        let force = phase_field_body_force(&phi, &mu, nx, ny);
        for (k, f) in force.iter().enumerate() {
            assert!(
                f[0].abs() < 1e-14 && f[1].abs() < 1e-14,
                "Uniform mu should give zero force at k={k}"
            );
        }
    }
    #[test]
    fn test_advection_uniform_phi() {
        let nx = 8;
        let ny = 8;
        let n = nx * ny;
        let phi = vec![0.5; n];
        let ux = vec![0.1; n];
        let uy = vec![0.05; n];
        let adv = advect_phase_field(&phi, &ux, &uy, nx, ny);
        for (k, &a) in adv.iter().enumerate() {
            assert!(
                a.abs() < 1e-14,
                "Advection of uniform phi should be zero: adv[{k}]={a}"
            );
        }
    }
    #[test]
    fn test_density_from_phase_field() {
        let phi = vec![0.0, 0.5, 1.0];
        let rho = density_from_phase_field(&phi, 1.0, 1000.0);
        assert!((rho[0] - 1.0).abs() < 1e-14);
        assert!((rho[1] - 500.5).abs() < 1e-14);
        assert!((rho[2] - 1000.0).abs() < 1e-14);
    }
    #[test]
    fn test_find_interface_cells() {
        let nx = 8;
        let ny = 8;
        let mut phi = vec![0.0; nx * ny];
        for y in ny / 2..ny {
            for x in 0..nx {
                phi[y * nx + x] = 1.0;
            }
        }
        let interface = find_interface_cells(&phi, nx, ny, 0.5);
        assert!(!interface.is_empty(), "Should find interface cells");
        for &k in &interface {
            let y = k / nx;
            assert!(
                y == ny / 2 - 1 || y == ny / 2,
                "Interface cell at unexpected y={y}"
            );
        }
    }
    #[test]
    fn test_allen_cahn_circle() {
        let nx = 64;
        let ny = 64;
        let mut ac = AllenCahn::new(nx, ny, 2.0, 1.0);
        ac.initialize_circle(32.0, 32.0, 10.0);
        let k_center = ac.idx(32, 32);
        assert!(
            ac.phi[k_center] > 0.9,
            "Center should be ~1, got {}",
            ac.phi[k_center]
        );
        let k_corner = ac.idx(0, 0);
        assert!(
            ac.phi[k_corner] < 0.1,
            "Far corner should be ~0, got {}",
            ac.phi[k_corner]
        );
    }
    #[test]
    fn test_allen_cahn_advection_stability() {
        let nx = 16;
        let ny = 32;
        let mut ac = AllenCahn::new(nx, ny, 2.0, 0.1);
        ac.initialize_tanh_interface(16.0, 1.0);
        let n = nx * ny;
        let ux = vec![0.01; n];
        let uy = vec![0.0; n];
        let dt = 0.01;
        for _ in 0..20 {
            ac.step_with_advection(dt, &ux, &uy);
        }
        for (k, &phi) in ac.phi.iter().enumerate() {
            assert!(
                (-0.2..=1.2).contains(&phi),
                "phi[{k}] = {phi} out of reasonable bounds"
            );
        }
    }
    #[test]
    fn test_allen_cahn_interface_length() {
        let nx = 64;
        let ny = 64;
        let mut ac = AllenCahn::new(nx, ny, 2.0, 1.0);
        ac.initialize_circle(32.0, 32.0, 10.0);
        let len = ac.interface_length();
        assert!(len > 0.0, "Interface length should be positive: {len}");
    }
    #[test]
    fn test_surface_tension_force_nonzero() {
        let nx = 64;
        let ny = 64;
        let params = make_params(0.01, 3.0, 0.072);
        let mut pf = PhaseField::new(nx, ny, params);
        pf.initialize_circle(32.0, 32.0, 10.0);
        let forces = pf.surface_tension_force();
        let total_mag: f64 = forces
            .iter()
            .map(|f| (f[0] * f[0] + f[1] * f[1]).sqrt())
            .sum();
        assert!(
            total_mag > 0.0,
            "Surface tension forces should be non-zero for a circle"
        );
    }
    #[test]
    fn test_measure_interface_width() {
        let nx = 4;
        let ny = 128;
        let mut ac = AllenCahn::new(nx, ny, 4.0, 1.0);
        ac.initialize_tanh_interface(64.0, 1.0);
        let width = measure_interface_width(&ac.phi, nx, ny, 0);
        assert!(width > 0.0, "Interface width should be positive: {width}");
        assert!(width < 40.0, "Interface width should be finite: {width}");
    }
    #[test]
    fn test_soy_mu_bulk_at_equilibrium() {
        let soy = SoyModel::new(-0.2, 0.1, 0.01);
        let phi_eq = soy.coexistence_values().1;
        let mu = soy.mu_bulk(phi_eq);
        assert!(
            mu.abs() < 1e-12,
            "mu_bulk at equilibrium should be ~0: {mu}"
        );
    }
    #[test]
    fn test_soy_coexistence_symmetric() {
        let soy = SoyModel::new(-0.2, 0.1, 0.01);
        let (phi_m, phi_p) = soy.coexistence_values();
        assert!(phi_m < 0.0 && phi_p > 0.0, "phi_m < 0, phi_p > 0");
        assert!((phi_m + phi_p).abs() < 1e-14, "Symmetric: phi_m = -phi_p");
    }
    #[test]
    fn test_soy_interface_width_positive() {
        let soy = SoyModel::new(-0.2, 0.1, 0.01);
        let xi = soy.interface_width_xi();
        assert!(xi > 0.0, "Interface width xi should be positive: {xi}");
    }
    #[test]
    fn test_soy_surface_tension_positive() {
        let soy = SoyModel::new(-0.2, 0.1, 0.01);
        let sigma = soy.surface_tension();
        assert!(sigma > 0.0, "Surface tension should be positive: {sigma}");
    }
    #[test]
    fn test_soy_bulk_pressure_ideal_limit() {
        let soy = SoyModel::new(0.0, 0.0, 0.01);
        let p = soy.bulk_pressure(1.0, 0.5);
        assert!((p - 1.0 / 3.0).abs() < 1e-12, "Ideal gas limit: p = {p}");
    }
    #[test]
    fn test_soy_equilibrium_radius_finite() {
        let soy = SoyModel::new(-0.2, 0.1, 0.01);
        let r = soy.equilibrium_radius(0.01);
        assert!(r > 0.0 && r.is_finite(), "Equilibrium radius = {r}");
    }
    #[test]
    fn test_convective_cahn_hilliard_conserves_mass() {
        let nx = 16;
        let ny = 16;
        let n = nx * ny;
        let phi: Vec<f64> = (0..n).map(|k| 0.05 * ((k as f64 * 0.4).sin())).collect();
        let ux = vec![0.01_f64; n];
        let uy = vec![0.0_f64; n];
        let mass_before: f64 = phi.iter().sum();
        let phi_new =
            convective_cahn_hilliard_step(&phi, &ux, &uy, nx, ny, 1e-4, -0.1, 0.1, 0.01, 0.1);
        let mass_after: f64 = phi_new.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "Mass conservation: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_convective_cahn_hilliard_no_nan() {
        let nx = 8;
        let ny = 8;
        let n = nx * ny;
        let phi = vec![0.3_f64; n];
        let ux = vec![0.0_f64; n];
        let uy = vec![0.0_f64; n];
        let phi_new =
            convective_cahn_hilliard_step(&phi, &ux, &uy, nx, ny, 0.01, -0.1, 0.1, 0.01, 0.1);
        for (k, &p) in phi_new.iter().enumerate() {
            assert!(!p.is_nan(), "phi[{k}] is NaN");
            assert!(!p.is_infinite(), "phi[{k}] is infinite");
        }
    }
    #[test]
    fn test_cahn_number_typical() {
        let cn = cahn_number(2.0, 100.0);
        assert!(
            (cn - 0.02).abs() < 1e-12,
            "Cahn number = {cn}, expected 0.02"
        );
    }
    #[test]
    fn test_phase_field_peclet_finite() {
        let pe = phase_field_peclet(0.1, 100.0, 0.01, -0.1);
        assert!(pe > 0.0 && pe.is_finite(), "Peclet = {pe}");
    }
    #[test]
    fn test_max_grid_spacing_resolves_interface() {
        let dx = max_grid_spacing(0.01, -0.1, 4);
        let xi = (0.01_f64 / 0.1_f64).sqrt();
        assert!(
            (dx - xi / 4.0).abs() < 1e-10,
            "dx_max = {dx}, expected {}",
            xi / 4.0
        );
    }
    #[test]
    fn test_max_grid_spacing_stable_phase() {
        let dx = max_grid_spacing(0.01, 0.1, 4);
        assert_eq!(dx, f64::INFINITY, "No interface for A > 0");
    }
    #[test]
    fn test_double_well_free_energy_minima() {
        let a = -0.2_f64;
        let b = 0.1_f64;
        let phi_eq = (-a / b).sqrt();
        let f_eq = double_well_free_energy(phi_eq, a, b);
        let f_0 = double_well_free_energy(0.0, a, b);
        assert!(f_eq < f_0, "f at equilibrium ({f_eq}) < f at phi=0 ({f_0})");
    }
    #[test]
    fn test_critical_nucleus_radius_positive() {
        let r = critical_nucleus_radius(0.01, 0.001);
        assert!(r > 0.0 && r.is_finite(), "Critical radius = {r}");
    }
    #[test]
    fn test_nucleation_energy_barrier_at_critical() {
        let sigma = 0.01;
        let dfv = 0.001;
        let r_crit = critical_nucleus_radius(sigma, dfv);
        let dg_crit = nucleation_energy_barrier(r_crit, dfv, sigma);
        let expected = critical_nucleation_barrier(sigma, dfv);
        assert!(
            (dg_crit - expected).abs() / expected < 1e-10,
            "dG_crit = {dg_crit}, expected {expected}"
        );
    }
    #[test]
    fn test_tolman_correction_reduces_surface_tension() {
        let sigma0 = 0.072;
        let r = 1e-9;
        let delta = 1e-10;
        let sigma_curved = tolman_surface_tension(sigma0, r, delta);
        assert!(
            sigma_curved < sigma0,
            "Tolman correction should reduce sigma: {sigma_curved}"
        );
    }
    #[test]
    fn test_mixture_density_pure_phases() {
        let rho = mixture_density(&[1.0, 0.0], &[1000.0, 1.0]);
        assert!((rho - 1000.0).abs() < 1e-12, "Pure liquid: rho = {rho}");
        let rho_v = mixture_density(&[0.0, 1.0], &[1000.0, 1.0]);
        assert!((rho_v - 1.0).abs() < 1e-12, "Pure vapor: rho = {rho_v}");
    }
    #[test]
    fn test_mixture_density_50_50() {
        let rho = mixture_density(&[0.5, 0.5], &[1000.0, 1.0]);
        assert!((rho - 500.5).abs() < 1e-12, "50/50 mix: rho = {rho}");
    }
    #[test]
    fn test_mixture_viscosity() {
        let mu = mixture_viscosity(&[0.6, 0.4], &[1e-3, 1e-5]);
        let expected = 0.6 * 1e-3 + 0.4 * 1e-5;
        assert!(
            (mu - expected).abs() < 1e-20,
            "mu = {mu}, expected {expected}"
        );
    }
    #[test]
    fn test_partition_of_unity_after_enforce() {
        let n = 4;
        let mut phi_a = vec![0.7_f64; n];
        let mut phi_b = vec![0.4_f64; n];
        let mut fields = vec![phi_a.clone(), phi_b.clone()];
        enforce_partition_of_unity(&mut fields, n);
        for (k, (&f0, &f1)) in fields[0].iter().zip(fields[1].iter()).enumerate() {
            let s = f0 + f1;
            assert!(
                (s - 1.0).abs() < 1e-12,
                "Partition of unity at k={k}: sum = {s}"
            );
        }
        let _ = phi_a.iter_mut();
        let _ = phi_b.iter_mut();
    }
    #[test]
    fn test_partition_of_unity_error_exact() {
        let phi_a = vec![0.5_f64; 4];
        let phi_b = vec![0.5_f64; 4];
        let err = partition_of_unity_error(&[&phi_a, &phi_b], 4);
        assert!(err.abs() < 1e-14, "Exact partition: error = {err}");
    }
    #[test]
    fn test_cahn_hilliard_free_energy_uniform() {
        let nx = 8;
        let ny = 8;
        let phi = vec![0.5_f64; nx * ny];
        let a = -0.1;
        let b = 0.1;
        let kappa = 0.01;
        let f = cahn_hilliard_free_energy(&phi, nx, ny, a, b, kappa);
        let f_bulk_expected = double_well_free_energy(0.5, a, b) * (nx * ny) as f64;
        assert!(
            (f - f_bulk_expected).abs() < 1e-12,
            "Uniform phi: F = {f}, expected {f_bulk_expected}"
        );
    }
    #[test]
    fn test_cahn_hilliard_chemical_potential_uniform() {
        let nx = 8;
        let ny = 8;
        let phi_val = 0.5;
        let phi = vec![phi_val; nx * ny];
        let a = -0.1;
        let b = 0.1;
        let kappa = 0.01;
        let mu = cahn_hilliard_chemical_potential(&phi, nx, ny, a, b, kappa);
        let expected = a * phi_val + b * phi_val.powi(3);
        for (k, &m) in mu.iter().enumerate() {
            assert!(
                (m - expected).abs() < 1e-12,
                "mu[{k}] = {m}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_cahn_hilliard_chemical_potential_at_equilibrium() {
        let a = -0.2_f64;
        let b = 0.1_f64;
        let phi_eq = (-a / b).sqrt();
        let nx = 4;
        let ny = 4;
        let phi = vec![phi_eq; nx * ny];
        let mu = cahn_hilliard_chemical_potential(&phi, nx, ny, a, b, 0.01);
        for (k, &m) in mu.iter().enumerate() {
            assert!(
                m.abs() < 1e-12,
                "mu at equilibrium should be ~0: mu[{k}] = {m}"
            );
        }
    }
    #[test]
    fn test_phase_field_lbm_interface_width_tanh_profile() {
        let nx = 1;
        let ny = 64;
        let epsilon = 4.0;
        let pf = PhaseFieldLbm::new(nx, ny, epsilon, 0.01);
        let mut phi = vec![0.0_f64; nx * ny];
        let y_c = 32.0;
        for (y, o) in phi.iter_mut().enumerate() {
            *o = 0.5 * (1.0 + ((y as f64 - y_c) / (2.0 * epsilon)).tanh());
        }
        let w = pf.compute_interface_width(&phi, 0, 0.1, 0.9);
        assert!(w > 4.0, "Interface width should be > 4 lattice units: {w}");
        assert!(
            w < 40.0,
            "Interface width should be < 40 lattice units: {w}"
        );
    }
    #[test]
    fn test_phase_field_lbm_interface_width_uniform_zero() {
        let pf = PhaseFieldLbm::new(8, 8, 2.0, 0.01);
        let phi = vec![0.0_f64; 64];
        let w = pf.compute_interface_width(&phi, 0, 0.1, 0.9);
        assert_eq!(w, 0.0, "Uniform phi=0 should give width 0: {w}");
    }
    #[test]
    fn test_phase_field_lbm_interface_width_uniform_one() {
        let pf = PhaseFieldLbm::new(8, 8, 2.0, 0.01);
        let phi = vec![1.0_f64; 64];
        let w = pf.compute_interface_width(&phi, 0, 0.1, 0.9);
        assert_eq!(w, 0.0, "Uniform phi=1 should give width 0: {w}");
    }
    #[test]
    fn test_phase_field_lbm_contact_angle_90_deg_symmetry() {
        let nx = 8;
        let ny = 8;
        let pf = PhaseFieldLbm::new(nx, ny, 2.0, 0.01);
        let mut phi = vec![0.5_f64; nx * ny];
        let phi_before: Vec<f64> = phi.clone();
        pf.apply_contact_angle(&mut phi, 90.0, 0);
        for x in 0..nx {
            let k = x;
            let k_above = nx + x;
            assert!(
                (phi[k] - phi_before[k_above]).abs() < 1e-12,
                "90° BC: phi[{x},0] should equal phi[{x},1]: {} vs {}",
                phi[k],
                phi_before[k_above]
            );
        }
    }
    #[test]
    fn test_phase_field_lbm_contact_angle_0_deg_wetting() {
        let nx = 4;
        let ny = 8;
        let pf = PhaseFieldLbm::new(nx, ny, 2.0, 0.01);
        let mut phi = vec![0.5_f64; nx * ny];
        pf.apply_contact_angle(&mut phi, 0.0, 0);
        for &p in &phi[..nx] {
            assert!(
                (0.0..=1.0).contains(&p),
                "phi at wall must be in [0,1]: {p}"
            );
        }
    }
    #[test]
    fn test_phase_field_lbm_capillary_pressure_uniform_phi_zero() {
        let pf = PhaseFieldLbm::new(8, 8, 2.0, 0.01);
        let phi = vec![0.5_f64; 64];
        let dp = pf.compute_capillary_pressure(&phi, 1e-6);
        for (k, &p) in dp.iter().enumerate() {
            assert!(
                p.abs() < 1e-10,
                "Uniform phi → zero capillary pressure at cell {k}: {p}"
            );
        }
    }
    #[test]
    fn test_phase_field_lbm_capillary_pressure_circular_droplet_sign() {
        let nx = 32;
        let ny = 32;
        let epsilon = 2.0;
        let sigma = 0.01;
        let pf = PhaseFieldLbm::new(nx, ny, epsilon, sigma);
        let cx = 16.0;
        let cy = 16.0;
        let r = 8.0;
        let mut phi = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let dist = (dx * dx + dy * dy).sqrt();
                phi[y * nx + x] = 0.5 * (1.0 - ((dist - r) / (2.0 * epsilon)).tanh());
            }
        }
        let dp = pf.compute_capillary_pressure(&phi, 1e-6);
        let nonzero_count = dp.iter().filter(|&&p| p.abs() > 1e-10).count();
        assert!(
            nonzero_count > 0,
            "Circular droplet should produce nonzero capillary pressure"
        );
    }
    #[test]
    fn test_phase_field_lbm_capillary_pressure_scales_with_sigma() {
        let nx = 16;
        let ny = 16;
        let epsilon = 2.0;
        let pf1 = PhaseFieldLbm::new(nx, ny, epsilon, 0.01);
        let pf2 = PhaseFieldLbm::new(nx, ny, epsilon, 0.02);
        let mut phi = vec![0.0_f64; nx * ny];
        for y in 0..ny {
            for x in 0..nx {
                let r = ((x as f64 - 8.0).powi(2) + (y as f64 - 8.0).powi(2)).sqrt();
                phi[y * nx + x] = 0.5 * (1.0 - ((r - 4.0) / (2.0 * epsilon)).tanh());
            }
        }
        let dp1 = pf1.compute_capillary_pressure(&phi, 1e-6);
        let dp2 = pf2.compute_capillary_pressure(&phi, 1e-6);
        let sum1: f64 = dp1.iter().map(|p| p.abs()).sum();
        let sum2: f64 = dp2.iter().map(|p| p.abs()).sum();
        assert!(
            (sum2 - 2.0 * sum1).abs() / (sum1 + 1e-30) < 1e-10,
            "Doubling sigma should double capillary pressure: sum1={sum1}, sum2={sum2}"
        );
    }
    #[test]
    fn test_phase_field_lbm_new_stores_params() {
        let pf = PhaseFieldLbm::new(10, 20, 3.0, 0.05);
        assert_eq!(pf.nx, 10);
        assert_eq!(pf.ny, 20);
        assert!((pf.epsilon - 3.0).abs() < 1e-12);
        assert!((pf.sigma - 0.05).abs() < 1e-12);
    }
    #[test]
    fn test_phase_field_lbm_interface_width_monotone_with_epsilon() {
        let nx = 1;
        let ny = 128;
        let y_c = 64.0;
        let build_phi = |eps: f64| -> Vec<f64> {
            (0..ny)
                .map(|y| 0.5 * (1.0 + ((y as f64 - y_c) / (2.0 * eps)).tanh()))
                .collect()
        };
        let pf1 = PhaseFieldLbm::new(nx, ny, 2.0, 0.01);
        let pf2 = PhaseFieldLbm::new(nx, ny, 6.0, 0.01);
        let w1 = pf1.compute_interface_width(&build_phi(2.0), 0, 0.1, 0.9);
        let w2 = pf2.compute_interface_width(&build_phi(6.0), 0, 0.1, 0.9);
        assert!(
            w2 > w1,
            "Wider epsilon should give wider interface: w1={w1}, w2={w2}"
        );
    }
}
