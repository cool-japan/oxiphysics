//! Unit tests for [`super`] (SPH sediment transport).

use super::*;

const TOL: f64 = 1.0e-10;

#[test]
fn test_settling_velocity_stokes() {
    let ws = settling_velocity_stokes(0.001, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    // Should be positive and physically reasonable
    assert!(ws > 0.0);
    // Stokes: ws = (s-1)*g*d^2 / (18*nu)
    let expected = (2650.0 / 1000.0 - 1.0) * 9.81 * 1.0e-6 / (18.0 * 1.0e-6);
    assert!((ws - expected).abs() < TOL);
}

#[test]
fn test_settling_velocity_ferguson_church() {
    // 0.5 mm quartz in water: R = 1.65,
    // w = R g d^2 / (18 nu + sqrt(0.75 R g d^3)) = 0.0710 m/s.
    let ws = settling_velocity_ferguson_church(0.5e-3, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    assert!((ws - 0.0710).abs() < 1.0e-3, "w = {ws}");
    // The drag correction makes it slower than Stokes' law for this grain.
    let ws_stokes = settling_velocity_stokes(0.5e-3, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    assert!(
        ws < ws_stokes,
        "F-C {ws} must stay below Stokes {ws_stokes}"
    );
}

#[test]
fn test_settling_velocity_ferguson_church_closed_form() {
    for &d in &[2.0e-5_f64, 1.0e-4, 5.0e-4, 2.0e-3, 1.0e-2] {
        let r = RHO_S_DEFAULT / RHO_F_DEFAULT - 1.0;
        let expected = r * G * d * d / (18.0 * NU_WATER + (0.75 * r * G * d * d * d).sqrt());
        let ws = settling_velocity_ferguson_church(d, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
        assert!(
            (ws - expected).abs() <= 1.0e-12 * expected,
            "d = {d}: {ws} vs {expected}"
        );
    }
}

#[test]
fn test_settling_velocity_ferguson_church_limits() {
    // Fine grains: Stokes regime.
    let fine = settling_velocity_ferguson_church(1.0e-5, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    let fine_stokes = settling_velocity_stokes(1.0e-5, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    assert!(
        (fine / fine_stokes - 1.0).abs() < 0.02,
        "{fine} vs {fine_stokes}"
    );
    // Coarse gravel: Newton regime w = sqrt(R g d / 0.75).
    let d = 0.05;
    let coarse = settling_velocity_ferguson_church(d, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    let newton = (1.65 * G * d / 0.75).sqrt();
    assert!((coarse / newton - 1.0).abs() < 0.01, "{coarse} vs {newton}");
}

#[test]
fn test_settling_velocity_ferguson_church_degenerate_inputs() {
    assert_eq!(
        settling_velocity_ferguson_church(0.0, 2650.0, 1000.0, 1.0e-6),
        0.0
    );
    assert_eq!(
        settling_velocity_ferguson_church(1.0e-3, 2650.0, 1000.0, 0.0),
        0.0
    );
    assert_eq!(
        settling_velocity_ferguson_church(1.0e-3, 2650.0, 0.0, 1.0e-6),
        0.0
    );
    assert_eq!(
        settling_velocity_ferguson_church(1.0e-3, 1000.0, 1000.0, 1.0e-6),
        0.0
    );
    // A grain lighter than the fluid rises with the mirrored velocity.
    let sink = settling_velocity_ferguson_church(1.0e-3, 1100.0, 1000.0, 1.0e-6);
    let rise = settling_velocity_ferguson_church(1.0e-3, 900.0, 1000.0, 1.0e-6);
    assert!(sink > 0.0 && rise < 0.0);
    assert!((rise + sink).abs() < 1.0e-12 * sink, "{rise} vs {sink}");
}

#[test]
fn test_settling_velocity_hindered() {
    let ws0 = 0.05;
    // Zero concentration => full velocity
    assert!((settling_velocity_hindered(ws0, 0.0, 4.65) - ws0).abs() < TOL);
    // High concentration => much lower
    let ws_high = settling_velocity_hindered(ws0, 0.5, 4.65);
    assert!(ws_high < ws0);
    assert!(ws_high > 0.0);
    // Full concentration => zero
    assert!((settling_velocity_hindered(ws0, 1.0, 4.65)).abs() < TOL);
}

#[test]
fn test_dimensionless_diameter() {
    let d_star = dimensionless_diameter(0.5e-3, RHO_S_DEFAULT, RHO_F_DEFAULT, NU_WATER);
    assert!(d_star > 1.0);
    assert!(d_star < 100.0); // D* for 0.5mm sand is ~ 12.6
}

#[test]
fn test_shields_parameter() {
    let theta = shields_parameter(1.0, RHO_S_DEFAULT, RHO_F_DEFAULT, 0.001);
    let expected = 1.0 / ((2650.0 - 1000.0) * 9.81 * 0.001);
    assert!((theta - expected).abs() < 1.0e-8);
}

#[test]
fn test_critical_shields_soulsby_whitehouse() {
    // For very small D* => approaches 0.3
    let theta_cr_small = critical_shields_soulsby_whitehouse(0.1);
    assert!(theta_cr_small > 0.2);
    // For large D* => approaches 0.055
    let theta_cr_large = critical_shields_soulsby_whitehouse(100.0);
    assert!(theta_cr_large > 0.04 && theta_cr_large < 0.07);
}

#[test]
fn test_is_motion_initiated() {
    assert!(is_motion_initiated(0.06, 0.047));
    assert!(!is_motion_initiated(0.03, 0.047));
}

#[test]
fn test_bed_shear_stress() {
    let tau = bed_shear_stress(1000.0, 0.003, 1.0);
    assert!((tau - 3.0).abs() < TOL);
}

#[test]
fn test_friction_coeff_manning() {
    let cf = friction_coeff_manning(0.025, 1.0);
    // Cf = g * n^2 / h^(1/3) = 9.81 * 0.000625 / 1.0 = 0.0061...
    let expected = 9.81 * 0.025 * 0.025;
    assert!((cf - expected).abs() < 1.0e-8);
}

#[test]
fn test_friction_velocity() {
    let u_star = friction_velocity(4.0, 1000.0);
    assert!((u_star - 0.0632455532).abs() < 1.0e-6);
}

#[test]
fn test_bedload_mpm() {
    // Below threshold => 0
    assert!((bedload_mpm(0.04, 0.047)).abs() < TOL);
    // Above threshold => positive
    let phi = bedload_mpm(0.1, 0.047);
    let expected = 8.0 * (0.1 - 0.047_f64).powf(1.5);
    assert!((phi - expected).abs() < 1.0e-10);
}

#[test]
fn test_bedload_van_rijn() {
    let phi = bedload_van_rijn(0.1, 0.047, 12.0);
    assert!(phi > 0.0);
    // Below threshold => 0
    assert!((bedload_van_rijn(0.03, 0.047, 12.0)).abs() < TOL);
}

#[test]
fn test_bedload_einstein() {
    let phi = bedload_einstein(0.1);
    assert!(phi > 0.0);
    assert!(phi.is_finite());
    assert!((bedload_einstein(0.0)).abs() < TOL);
}

#[test]
fn test_slope_correction() {
    let theta_cr = slope_correction(0.047, 0.0, std::f64::consts::FRAC_PI_6); // flat bed
    assert!((theta_cr - 0.047).abs() < 1.0e-6);
    // Steep slope => lower critical
    let theta_steep = slope_correction(0.047, 0.3, std::f64::consts::FRAC_PI_6);
    assert!(theta_steep < 0.047);
}

#[test]
fn test_rouse_number() {
    let z = rouse_number(0.01, 0.05);
    // Z = 0.01 / (0.41 * 0.05)
    let expected = 0.01 / (0.41 * 0.05);
    assert!((z - expected).abs() < 1.0e-10);
}

#[test]
fn test_rouse_profile() {
    let c = rouse_profile(0.5, 1.0, 0.05, 0.01, 1.0);
    assert!(c > 0.0);
    assert!(c <= 0.01);
    // At reference height => c_a
    let c_ref = rouse_profile(0.05, 1.0, 0.05, 0.01, 1.0);
    assert!((c_ref - 0.01).abs() < 1.0e-6);
}

#[test]
fn test_rouse_depth_averaged() {
    let c_avg = rouse_depth_averaged(1.0, 0.05, 0.01, 1.0, 100);
    assert!(c_avg > 0.0);
    assert!(c_avg < 0.01); // averaged should be less than reference
}

#[test]
fn test_erosion_deposition() {
    // Erosion above critical
    let e = erosion_partheniades(1.0, 0.5, 0.001);
    assert!((e - 0.001).abs() < 1.0e-10);
    // No erosion below critical
    assert!((erosion_partheniades(0.3, 0.5, 0.001)).abs() < TOL);
    // Deposition below critical
    let d = deposition_krone(0.01, 0.005, 0.1, 0.5);
    assert!(d > 0.0);
    // No deposition above critical
    assert!((deposition_krone(0.01, 0.005, 0.6, 0.5)).abs() < TOL);
}

#[test]
fn test_exner_bed_change() {
    let dz = exner_bed_change(0.0, 0.001, 0.0005, 0.4);
    // dz = 1/(1-0.4) * (0.001 - 0.0005) = 1.6667 * 0.0005
    let expected = (1.0 / 0.6) * 0.0005;
    assert!((dz - expected).abs() < 1.0e-10);
}

#[test]
fn test_sediment_particle_creation() {
    let p = SedimentSphParticle::new_fluid([1.0, 2.0, 3.0], [0.1, 0.0, 0.0], 0.01, 1000.0, 0.05);
    assert!(!p.is_bed);
    assert!((p.concentration).abs() < TOL);
    let b = SedimentSphParticle::new_bed([1.0, 2.0, 0.0], 0.01, 2650.0, 0.05);
    assert!(b.is_bed);
    assert!((b.concentration - 0.6).abs() < TOL);
}

#[test]
fn test_mixture_density() {
    let p = SedimentSphParticle {
        pos: [0.0; 3],
        vel: [0.0; 3],
        mass: 1.0,
        density: 1000.0,
        concentration: 0.1,
        diffusivity: 0.0,
        h: 0.05,
        is_bed: false,
    };
    let rho_mix = p.mixture_density(2650.0, 1000.0);
    let expected = 0.1 * 2650.0 + 0.9 * 1000.0;
    assert!((rho_mix - expected).abs() < TOL);
}

#[test]
fn test_scour_config() {
    let sc = ScourConfig::cylinder([0.0, 0.0, 0.0], 1.0);
    // At the structure surface
    let amp_surface = sc.local_amplification([0.5, 0.0, 0.0]);
    assert!((amp_surface - 2.5).abs() < 1.0e-6);
    // Far away
    let amp_far = sc.local_amplification([10.0, 0.0, 0.0]);
    assert!((amp_far - 1.0).abs() < 1.0e-6);
    // Equilibrium depth
    let d_se = sc.equilibrium_depth(1.0);
    assert!((d_se - 2.4).abs() < 1.0e-6);
}

#[test]
fn test_hiding_exposure() {
    // Same size as mean => xi = 1
    let xi = hiding_exposure_egiazaroff(0.001, 0.001);
    assert!((xi - 1.0).abs() < 1.0e-6);
    // Finer grain => xi > 1 (harder to move)
    let xi_fine = hiding_exposure_egiazaroff(0.0005, 0.001);
    assert!(xi_fine > 1.0);
}

#[test]
fn test_multi_fraction_bedload() {
    let fracs = vec![
        GrainFraction {
            diameter: 0.5e-3,
            rho_s: 2650.0,
            fraction: 0.5,
        },
        GrainFraction {
            diameter: 2.0e-3,
            rho_s: 2650.0,
            fraction: 0.5,
        },
    ];
    let rates = multi_fraction_bedload(&fracs, 5.0, 1000.0, 0.047);
    assert_eq!(rates.len(), 2);
    // Finer fraction should have higher transport rate
    assert!(rates[0] > rates[1]);
}

#[test]
fn test_flocculation() {
    let params = FlocculationParams::estuarine_mud();
    let d_eq = params.equilibrium_floc_diameter(5.0);
    assert!(d_eq > params.d_primary);
    assert!(d_eq <= params.d_max);
    // Low shear => larger flocs
    let d_eq_low = params.equilibrium_floc_diameter(1.0);
    let d_eq_high = params.equilibrium_floc_diameter(100.0);
    assert!(d_eq_low > d_eq_high);
}

#[test]
fn test_floc_density() {
    let params = FlocculationParams::estuarine_mud();
    let rho_floc = params.floc_density(0.5e-3, 2650.0, 1000.0);
    // Floc should be less dense than primary grain
    assert!(rho_floc < 2650.0);
    assert!(rho_floc > 1000.0);
}

#[test]
fn test_solver_creation() {
    let mut solver = SedimentSphSolver::new(SedimentSphConfig::default());
    solver.add_fluid_particle([0.0, 0.0, 0.1], [1.0, 0.0, 0.0], 0.01, 0.05);
    solver.add_bed_particle([0.0, 0.0, 0.0], 0.01, 0.05);
    assert_eq!(solver.fluid_count(), 1);
    assert_eq!(solver.bed_count(), 1);
}

#[test]
fn test_avalanche_correction() {
    let dz = avalanche_correction(1.0, 0.0, 1.0, std::f64::consts::FRAC_PI_6, 0.5);
    // Slope = 1.0, tan(30 deg) ~ 0.577
    // Excess = (1.0 - 0.577) * 1.0 ~ 0.423
    // dz = -0.5 * 0.423 * 0.5 ~ -0.106
    assert!(dz < 0.0);
}

#[test]
fn test_classify_transport() {
    assert_eq!(classify_transport(3.0), TransportMode::BedLoad);
    assert_eq!(classify_transport(2.0), TransportMode::Mixed);
    assert_eq!(classify_transport(1.0), TransportMode::SuspendedLoad);
    assert_eq!(classify_transport(0.5), TransportMode::WashLoad);
}

#[test]
fn test_cerc_longshore() {
    let q = cerc_longshore_transport(0.39, 1.0, 0.2, 2650.0, 1000.0, 0.4, 0.78);
    assert!(q > 0.0);
    // Zero wave height => zero transport
    let q0 = cerc_longshore_transport(0.39, 0.0, 0.2, 2650.0, 1000.0, 0.4, 0.78);
    assert!(q0.abs() < TOL);
}

#[test]
fn test_bed_interaction_force() {
    // Above bed => no force
    let f = bed_interaction_force(
        [0.0, 0.0, 1.0],
        [0.0; 3],
        0.0,
        [0.0, 0.0, 1.0],
        1000.0,
        10.0,
    );
    assert!((len3(f)).abs() < TOL);
    // Below bed => repulsion
    let f2 = bed_interaction_force(
        [0.0, 0.0, -0.01],
        [0.0; 3],
        0.0,
        [0.0, 0.0, 1.0],
        1000.0,
        10.0,
    );
    assert!(f2[2] > 0.0);
}

#[test]
fn test_submerged_weight() {
    let w = submerged_weight(0.001, 2650.0, 1000.0);
    let vol = PI / 6.0 * 1.0e-9;
    let expected = 1650.0 * vol * 9.81;
    assert!((w - expected).abs() < 1.0e-12);
}

#[test]
fn test_turbidity_front() {
    let g_prime = reduced_gravity(0.01, 2650.0, 1000.0);
    assert!(g_prime > 0.0);
    let u_f = turbidity_front_velocity(g_prime, 1.0, 0.7);
    assert!(u_f > 0.0);
}

#[test]
fn test_entrainment_coefficient() {
    let e0 = entrainment_coefficient(0.0);
    assert!((e0 - 0.075).abs() < 1.0e-6);
    let e_high = entrainment_coefficient(1.0);
    assert!(e_high < e0);
}

#[test]
fn test_armour_ratio() {
    assert!((armour_ratio(2.0e-3, 1.0e-3) - 2.0).abs() < TOL);
    assert!(is_armoured(2.0e-3, 1.0e-3, 1.5));
    assert!(!is_armoured(1.0e-3, 1.0e-3, 1.5));
}

#[test]
fn test_bed_column_update() {
    let mut col = BedColumn::new(0.0, 0.0, 0.0);
    col.update(0.01);
    assert!((col.z_bed - 0.01).abs() < TOL);
    assert!((col.cumulative_deposition - 0.01).abs() < TOL);
    col.update(-0.005);
    assert!((col.z_bed - 0.005).abs() < TOL);
    assert!((col.cumulative_erosion - 0.005).abs() < TOL);
}

#[test]
fn test_morphological_step() {
    let mut solver = SedimentSphSolver::new(SedimentSphConfig::default());
    // Five clear-water particles within 2h of the column, all 0.05 m above the bed.
    for i in 0..5 {
        let x = i as f64 * 0.02;
        solver.add_fluid_particle([x, 0.0, 0.05], [0.5, 0.0, 0.0], 0.01, 0.05);
    }
    solver.add_bed_column(0.05, 0.0, 0.0);
    let dt = 0.01;
    // The 0.5 m/s current is far above the threshold of motion.
    let tau = solver.compute_bed_shear(0.5, 0.05);
    let theta = shields_parameter(tau, RHO_S_DEFAULT, RHO_F_DEFAULT, solver.config.d50);
    assert!(theta > solver.config.theta_cr, "theta = {theta}");
    solver.morphological_step(dt);
    // A lone column has no flux divergence and the water carries no sediment, so
    // the bed erodes at the pick-up rate E = alpha w_s T^1.5 (Exner, 1 / (1 - p)).
    let cfg = &solver.config;
    let w_s = settling_velocity_ferguson_church(cfg.d50, cfg.rho_s, cfg.rho_f, NU_WATER);
    let stage = (theta - cfg.theta_cr) / cfg.theta_cr;
    let expected = -cfg.erosion_coeff * w_s * stage.powf(1.5) * dt / (1.0 - cfg.porosity);
    let z = solver.bed_columns[0].z_bed;
    assert!(z < 0.0, "the bed must erode, got {z}");
    assert!(
        (z - expected).abs() <= 1.0e-12 * expected.abs(),
        "{z} vs {expected}"
    );
    assert!((solver.time - dt).abs() < 1.0e-15);
}

#[test]
fn test_ripple_and_dune() {
    let rh = ripple_height_van_rijn(0.5e-3, 12.0, 3.0);
    assert!(rh > 0.0);
    let rl = ripple_wavelength(rh);
    assert!((rl - 7.0 * rh).abs() < TOL);
    let dh = dune_height_yalin(2.0, 0.5);
    assert!(dh > 0.0);
    // Supercritical => no dunes
    assert!((dune_height_yalin(2.0, 1.5)).abs() < TOL);
}

#[test]
fn test_sph_diffusion_pair() {
    use crate::kernel::{CubicSplineKernel, SphKernel};
    let xi = [0.0, 0.0, 0.0];
    let xj = [0.05, 0.0, 0.0];
    let rate = sph_diffusion_pair(0.1, 0.2, 1.0e-3, 1.0e-3, xi, xj, 0.1, 1.0e-4);
    // c_j > c_i: sediment diffuses from the richer neighbour into particle i.
    assert!(rate > 0.0, "rate = {rate}");
    // (eps_i + eps_j) (c_i - c_j) r W'(r) / (r^2 + eta^2)
    let r = 0.05;
    let expected = 2.0e-3 * (0.1 - 0.2) * r * CubicSplineKernel.grad_w(r, 0.1) / (r * r + 1.0e-4);
    assert!((rate - expected).abs() <= 1.0e-12 * expected.abs());
    // The reverse pair loses sediment and a uniform field does not diffuse.
    let reverse = sph_diffusion_pair(0.2, 0.1, 1.0e-3, 1.0e-3, xi, xj, 0.1, 1.0e-4);
    assert!(reverse < 0.0);
    assert!((reverse + rate).abs() <= 1.0e-12 * rate);
    assert_eq!(
        sph_diffusion_pair(0.3, 0.3, 1.0e-3, 1.0e-3, xi, xj, 0.1, 1.0e-4),
        0.0
    );
}

#[test]
fn test_sph_diffusion_rate_flattens_a_peak_and_conserves_sediment() {
    let (h, spacing, n) = (0.1, 0.05, 9_usize);
    let positions: Vec<[f64; 3]> = (0..n).map(|k| [k as f64 * spacing, 0.0, 0.0]).collect();
    let concentrations: Vec<f64> = positions
        .iter()
        .map(|p| (-((p[0] - 0.2) / 0.08).powi(2)).exp())
        .collect();
    let diffusivities = vec![1.0e-3; n];
    let masses = vec![1.0; n];
    let densities = vec![1000.0; n];
    let all: Vec<usize> = (0..n).collect();
    let rates: Vec<f64> = (0..n)
        .map(|i| {
            sph_diffusion_rate(
                i,
                &positions,
                &concentrations,
                &diffusivities,
                &masses,
                &densities,
                &all,
                h,
                i,
            )
        })
        .collect();
    let peak = n / 2;
    assert!(
        rates[peak] < 0.0,
        "peak must lose sediment: {}",
        rates[peak]
    );
    assert!(
        rates[0] > 0.0 && rates[n - 1] > 0.0,
        "tails gain: {rates:?}"
    );
    // Every pair appears twice with opposite sign, so the volume-weighted total vanishes.
    let total: f64 = rates.iter().map(|r| r * (masses[0] / densities[0])).sum();
    let scale: f64 = rates
        .iter()
        .map(|r| (r * (masses[0] / densities[0])).abs())
        .sum();
    assert!(
        total.abs() <= 1.0e-12 * scale,
        "total = {total}, scale = {scale}"
    );
}

#[test]
fn test_total_suspended_load() {
    let particles = vec![
        SedimentSphParticle {
            pos: [0.0; 3],
            vel: [0.0; 3],
            mass: 1.0,
            density: 1000.0,
            concentration: 0.05,
            diffusivity: 0.0,
            h: 0.05,
            is_bed: false,
        },
        SedimentSphParticle {
            pos: [0.0; 3],
            vel: [0.0; 3],
            mass: 1.0,
            density: 1000.0,
            concentration: 0.03,
            diffusivity: 0.0,
            h: 0.05,
            is_bed: false,
        },
        SedimentSphParticle {
            pos: [0.0; 3],
            vel: [0.0; 3],
            mass: 1.0,
            density: 2650.0,
            concentration: 0.6,
            diffusivity: 0.0,
            h: 0.05,
            is_bed: true,
        },
    ];
    let tsl = total_suspended_load(&particles);
    // Only fluid particles: (0.05 * 1/1000) + (0.03 * 1/1000) = 8e-5
    let expected = 0.05 * 0.001 + 0.03 * 0.001;
    assert!((tsl - expected).abs() < 1.0e-10);

    let avg = average_concentration(&particles);
    assert!((avg - 0.04).abs() < 1.0e-10);
}

#[test]
fn test_transport_capacity_eh() {
    let qt = transport_capacity_eh(1.0, 0.5e-3, 2650.0, 1000.0, 0.003);
    assert!(qt > 0.0);
    assert!(qt.is_finite());
    // Zero velocity => zero transport
    assert!((transport_capacity_eh(0.0, 0.5e-3, 2650.0, 1000.0, 0.003)).abs() < TOL);
}

// ── Exner grid mass conservation (E5) ────────────────────────────────────

#[test]
fn test_exner_grid_mass_conservation() {
    // 4×4 grid with uniform outward flux: sum of bed elevation changes
    // should match -sum(outflux * dt) / (1 - porosity) within 1%.
    let nx = 4_usize;
    let ny = 4_usize;
    let dx = 1.0_f64;
    let dy = 1.0_f64;
    let porosity = 0.4_f64;
    let dt = 0.1_f64;

    let mut grid = ExnerGrid::new(nx, ny, dx, dy, porosity, 0.0);

    // Set a constant x-flux of 0.01 m²/s everywhere and zero y-flux
    let q_val = 0.01_f64;
    for v in &mut grid.flux_x {
        *v = q_val;
    }

    let eta_before: f64 = grid.bed_elevation.iter().sum();
    grid.compute_bed_change(dt);
    let eta_after: f64 = grid.bed_elevation.iter().sum();

    let total_change = eta_after - eta_before; // Σ dη

    // For a uniform flux q_val in x, the interior divergence ∂q_x/∂x = 0
    // (forward difference of equal values).  At x-boundaries (i = nx-1),
    // the backward difference is also 0.  So the bed should not change.
    // This validates the finite-difference scheme is consistent.
    assert!(
        total_change.abs() < 1e-12,
        "uniform flux should give zero divergence: total_change={total_change}"
    );
}

#[test]
fn test_exner_grid_nonzero_flux_divergence() {
    // A linearly increasing flux qx[j*nx+i] = (i+1) * 0.001 gives
    // ∂qx/∂x = 0.001/dx everywhere.  The bed should erode uniformly.
    let nx = 5_usize;
    let ny = 3_usize;
    let dx = 1.0_f64;
    let dy = 1.0_f64;
    let porosity = 0.4_f64;
    let dt = 1.0_f64;
    let lambda_p = 1.0 - porosity;

    let mut grid = ExnerGrid::new(nx, ny, dx, dy, porosity, 0.0);

    let q_slope = 0.001_f64;
    for j in 0..ny {
        for i in 0..nx {
            grid.flux_x[j * nx + i] = (i + 1) as f64 * q_slope;
        }
    }

    let eta_before: Vec<f64> = grid.bed_elevation.clone();
    grid.compute_bed_change(dt);

    // Interior cells (i < nx-1): forward diff = q_slope/dx
    // Exner: dη = -dt * q_slope / (dx * (1-λ_p))
    let expected_deta_interior = -dt * q_slope / (dx * lambda_p);

    for j in 0..ny {
        for i in 0..nx - 1 {
            let deta = grid.bed_elevation[j * nx + i] - eta_before[j * nx + i];
            assert!(
                (deta - expected_deta_interior).abs() < 1e-12,
                "interior cell ({i},{j}): deta={deta}, expected={expected_deta_interior}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Exner / morphological_step
// ---------------------------------------------------------------------------

/// Height of the single fluid particle above each flat-bed column (m).
const PARTICLE_HEIGHT: f64 = 0.05;
/// Smoothing length of the fluid particles: 2h is below the column spacing,
/// so every column only sees the particle straight above it.
const PARTICLE_H: f64 = 0.04;
/// Spacing of the bed columns (m).
const COLUMN_SPACING: f64 = 0.1;

fn set_concentration(solver: &mut SedimentSphSolver, concentration: f64) {
    if let Some(p) = solver.particles.last_mut() {
        p.concentration = concentration;
    }
}

/// One bed column per `xs` entry at elevation 0 and one fluid particle above
/// it with the given horizontal velocity and concentration.
fn flat_bed_solver(
    config: SedimentSphConfig,
    xs: &[f64],
    velocities: &[[f64; 3]],
    concentrations: &[f64],
) -> SedimentSphSolver {
    let mut solver = SedimentSphSolver::new(config);
    for ((&x, &vel), &c) in xs.iter().zip(velocities).zip(concentrations) {
        solver.add_bed_column(x, 0.0, 0.0);
        solver.add_fluid_particle([x, 0.0, PARTICLE_HEIGHT], vel, 0.01, PARTICLE_H);
        set_concentration(&mut solver, c);
    }
    solver
}

fn bedload_only() -> SedimentSphConfig {
    SedimentSphConfig {
        enable_suspended: false,
        ..SedimentSphConfig::default()
    }
}

fn suspended_only() -> SedimentSphConfig {
    SedimentSphConfig {
        enable_bedload: false,
        ..SedimentSphConfig::default()
    }
}

fn row_positions(n: usize) -> Vec<f64> {
    (0..n).map(|k| k as f64 * COLUMN_SPACING).collect()
}

/// Bed-load flux magnitude (m^2/s) for a near-bed speed, via the public helpers.
fn bedload_flux_magnitude(solver: &SedimentSphSolver, speed: f64) -> f64 {
    let tau = solver.compute_bed_shear(speed, PARTICLE_HEIGHT);
    let phi = solver.compute_bedload_rate(tau);
    let cfg = &solver.config;
    bedload_volumetric(phi, cfg.d50, cfg.rho_s, cfg.rho_f)
}

fn bed_levels(solver: &SedimentSphSolver) -> Vec<f64> {
    solver.bed_columns.iter().map(|c| c.z_bed).collect()
}

#[test]
fn test_uniform_flow_over_a_flat_bed_changes_nothing_by_bedload() {
    let xs = row_positions(5);
    let velocities = vec![[0.6, 0.0, 0.0]; 5];
    let mut solver = flat_bed_solver(bedload_only(), &xs, &velocities, &[0.0; 5]);
    // The flow is well above the threshold of motion, so there is real transport.
    let tau = solver.compute_bed_shear(0.6, PARTICLE_HEIGHT);
    assert!(
        solver.compute_bedload_rate(tau) > 1.0,
        "flow must move sediment"
    );
    solver.morphological_step(0.01);
    for (k, z) in bed_levels(&solver).iter().enumerate() {
        assert!(z.abs() < 1.0e-15, "column {k} moved to {z}");
    }
    assert!((solver.time - 0.01).abs() < 1.0e-15);
}

#[test]
fn test_uniform_flow_over_a_2d_bed_changes_nothing_by_bedload() {
    let mut solver = SedimentSphSolver::new(bedload_only());
    for i in 0..4 {
        for j in 0..4 {
            let (x, y) = (i as f64 * COLUMN_SPACING, j as f64 * COLUMN_SPACING);
            solver.add_bed_column(x, y, 0.0);
            solver.add_fluid_particle([x, y, PARTICLE_HEIGHT], [0.4, 0.3, 0.0], 0.01, PARTICLE_H);
        }
    }
    solver.morphological_step(0.01);
    assert!(bed_levels(&solver).iter().all(|z| z.abs() < 1.0e-15));
}

#[test]
fn test_subcritical_flow_with_clear_water_changes_nothing() {
    let xs = row_positions(3);
    let velocities = vec![[0.05, 0.0, 0.0]; 3];
    let mut solver = flat_bed_solver(SedimentSphConfig::default(), &xs, &velocities, &[0.0; 3]);
    let tau = solver.compute_bed_shear(0.05, PARTICLE_HEIGHT);
    assert!(shields_parameter(tau, 2650.0, 1000.0, 0.5e-3) < THETA_CR_DEFAULT);
    solver.morphological_step(0.01);
    assert!(bed_levels(&solver).iter().all(|z| z.abs() < 1.0e-18));
}

#[test]
fn test_accelerating_flow_erodes_and_decelerating_flow_deposits() {
    let xs = row_positions(5);
    let speeds: Vec<f64> = (0..5).map(|k| 0.50 + 0.05 * k as f64).collect();
    let velocities: Vec<[f64; 3]> = speeds.iter().map(|&u| [u, 0.0, 0.0]).collect();
    let dt = 0.01;
    let mut solver = flat_bed_solver(bedload_only(), &xs, &velocities, &[0.0; 5]);
    let q: Vec<f64> = speeds
        .iter()
        .map(|&u| bedload_flux_magnitude(&solver, u))
        .collect();
    solver.morphological_step(dt);
    let porosity_factor = 1.0 / (1.0 - solver.config.porosity);
    let z = bed_levels(&solver);
    // Interior columns: central difference of the flux over the 2 * spacing stencil.
    for k in 1..4 {
        let div = (q[k + 1] - q[k - 1]) / (2.0 * COLUMN_SPACING);
        let expected = -porosity_factor * div * dt;
        assert!(
            div > 0.0 && z[k] < 0.0,
            "column {k}: accelerating flow must erode"
        );
        assert!(
            (z[k] - expected).abs() <= 1.0e-12 * expected.abs(),
            "column {k}: dz = {}, expected {expected}",
            z[k]
        );
    }
    // End columns use the one-sided difference towards their only neighbour.
    let first = -porosity_factor * (q[1] - q[0]) / COLUMN_SPACING * dt;
    let last = -porosity_factor * (q[4] - q[3]) / COLUMN_SPACING * dt;
    assert!((z[0] - first).abs() <= 1.0e-12 * first.abs());
    assert!((z[4] - last).abs() <= 1.0e-12 * last.abs());

    // Reversing the velocity ramp reverses the sign: the bed accretes.
    let reversed: Vec<[f64; 3]> = velocities.iter().rev().copied().collect();
    let mut decelerating = flat_bed_solver(bedload_only(), &xs, &reversed, &[0.0; 5]);
    decelerating.morphological_step(dt);
    for (k, z_dec) in bed_levels(&decelerating).iter().enumerate().take(4).skip(1) {
        assert!(*z_dec > 0.0, "column {k}: decelerating flow must deposit");
        assert!((z_dec + z[4 - k]).abs() <= 1.0e-12 * z_dec.abs());
    }
    assert!(decelerating.total_deposition() > 0.0);
    assert!(solver.total_erosion() > 0.0);
}

#[test]
fn test_flux_direction_follows_the_flow_in_the_horizontal_plane() {
    // Flow along +y with a speed that grows along y erodes the middle column.
    let mut solver = SedimentSphSolver::new(bedload_only());
    let speeds = [0.5, 0.6, 0.7];
    for (k, &u) in speeds.iter().enumerate() {
        let y = k as f64 * COLUMN_SPACING;
        solver.add_bed_column(0.0, y, 0.0);
        solver.add_fluid_particle([0.0, y, PARTICLE_HEIGHT], [0.0, u, 0.0], 0.01, PARTICLE_H);
    }
    let q: Vec<f64> = speeds
        .iter()
        .map(|&u| bedload_flux_magnitude(&solver, u))
        .collect();
    solver.morphological_step(0.01);
    let div = (q[2] - q[0]) / (2.0 * COLUMN_SPACING);
    let expected = -div * 0.01 / (1.0 - solver.config.porosity);
    let z_mid = solver.bed_columns[1].z_bed;
    assert!(z_mid < 0.0);
    assert!(
        (z_mid - expected).abs() <= 1.0e-12 * expected.abs(),
        "{z_mid} vs {expected}"
    );
}

#[test]
fn test_vertical_velocity_does_not_scour_the_bed() {
    let mut solver = flat_bed_solver(
        SedimentSphConfig::default(),
        &[0.0],
        &[[0.0, 0.0, -2.0]],
        &[0.0],
    );
    solver.morphological_step(0.01);
    assert_eq!(solver.bed_columns[0].z_bed, 0.0);
}

#[test]
fn test_clear_water_deposition_follows_the_settling_flux() {
    let dt = 0.01;
    let c = 0.1;
    let mut solver = flat_bed_solver(suspended_only(), &[0.0], &[[0.0; 3]], &[c]);
    solver.morphological_step(dt);
    let w_s = settling_velocity_ferguson_church(0.5e-3, 2650.0, 1000.0, NU_WATER);
    // Still water: tau_b = 0, so D = w_s c; the bed rises by D dt / (1 - p).
    let expected = w_s * c * dt / (1.0 - 0.4);
    let z = solver.bed_columns[0].z_bed;
    assert!(z > 0.0);
    assert!(
        (z - expected).abs() <= 1.0e-12 * expected,
        "{z} vs {expected}"
    );
    assert!((solver.total_deposition() - z).abs() < 1.0e-18);
}

#[test]
fn test_strong_current_erodes_at_the_pickup_rate() {
    let dt = 0.01;
    let speed = 0.6;
    let mut solver = flat_bed_solver(suspended_only(), &[0.0], &[[speed, 0.0, 0.0]], &[0.0]);
    let tau = solver.compute_bed_shear(speed, PARTICLE_HEIGHT);
    let theta = shields_parameter(tau, 2650.0, 1000.0, 0.5e-3);
    assert!(theta > THETA_CR_DEFAULT);
    solver.morphological_step(dt);
    let w_s = settling_velocity_ferguson_church(0.5e-3, 2650.0, 1000.0, NU_WATER);
    let stage = (theta - THETA_CR_DEFAULT) / THETA_CR_DEFAULT;
    let expected = -0.01 * w_s * stage.powf(1.5) * dt / (1.0 - 0.4);
    let z = solver.bed_columns[0].z_bed;
    assert!(z < 0.0);
    assert!(
        (z - expected).abs() <= 1.0e-12 * expected.abs(),
        "{z} vs {expected}"
    );
    assert!((solver.total_erosion() + z).abs() < 1.0e-18);
}

#[test]
fn test_erosion_coefficient_and_deposition_threshold_come_from_the_config() {
    let speed = 0.6;
    let erode = |alpha: f64| {
        let cfg = SedimentSphConfig {
            erosion_coeff: alpha,
            ..suspended_only()
        };
        let mut solver = flat_bed_solver(cfg, &[0.0], &[[speed, 0.0, 0.0]], &[0.0]);
        solver.morphological_step(0.01);
        solver.bed_columns[0].z_bed
    };
    let (z1, z2) = (erode(0.01), erode(0.02));
    assert!(z1 < 0.0 && (z2 / z1 - 2.0).abs() < 1.0e-12, "{z1}, {z2}");

    // Deposition at a shear stress between 0 and tau_cd: D = w_s c (1 - tau_b / tau_cd).
    let weak = 0.1;
    let deposit = |tau_cd: f64| {
        let cfg = SedimentSphConfig {
            tau_cd,
            ..suspended_only()
        };
        let mut solver = flat_bed_solver(cfg, &[0.0], &[[weak, 0.0, 0.0]], &[0.2]);
        let tau_b = solver.compute_bed_shear(weak, PARTICLE_HEIGHT);
        solver.morphological_step(0.01);
        (solver.bed_columns[0].z_bed, tau_b)
    };
    let (z_low, tau_b) = deposit(0.5);
    let (z_high, _) = deposit(1.0);
    let w_s = settling_velocity_ferguson_church(0.5e-3, 2650.0, 1000.0, NU_WATER);
    for (z, tau_cd) in [(z_low, 0.5), (z_high, 1.0)] {
        let expected = w_s * 0.2 * (1.0 - tau_b / tau_cd) * 0.01 / (1.0 - 0.4);
        assert!(tau_b < tau_cd);
        assert!(
            (z - expected).abs() <= 1.0e-12 * expected,
            "{z} vs {expected}"
        );
    }
    assert!(
        z_high > z_low,
        "a higher tau_cd keeps depositing at larger stress"
    );
}

#[test]
fn test_concentration_average_uses_the_same_particles_as_the_velocity_average() {
    let mut solver = SedimentSphSolver::new(suspended_only());
    solver.add_bed_column(0.0, 0.0, 0.0);
    // Inside the support of the column (|dx| = 0.01 < 2h = 0.08).
    solver.add_fluid_particle([0.01, 0.0, PARTICLE_HEIGHT], [0.0; 3], 0.01, PARTICLE_H);
    set_concentration(&mut solver, 0.1);
    // Outside the support (0.09 > 0.08) though within 0.1 m of the column.
    solver.add_fluid_particle([0.09, 0.0, PARTICLE_HEIGHT], [0.0; 3], 0.01, PARTICLE_H);
    set_concentration(&mut solver, 0.5);
    solver.morphological_step(0.01);
    let w_s = settling_velocity_ferguson_church(0.5e-3, 2650.0, 1000.0, NU_WATER);
    let expected = w_s * 0.1 * 0.01 / (1.0 - 0.4);
    let z = solver.bed_columns[0].z_bed;
    assert!(
        (z - expected).abs() <= 1.0e-12 * expected,
        "{z} vs {expected}"
    );
}

#[test]
fn test_water_depth_comes_from_the_highest_particle_above_the_column() {
    let (u, c) = (0.1, 0.2);
    let build = |top: Option<f64>| {
        let mut solver = flat_bed_solver(suspended_only(), &[0.0], &[[u, 0.0, 0.0]], &[c]);
        if let Some(z) = top {
            // Above the near-bed sampling height: only sets the water depth.
            solver.add_fluid_particle([0.0, 0.0, z], [3.0, 0.0, 0.0], 0.01, PARTICLE_H);
            set_concentration(&mut solver, 0.9);
        }
        solver.morphological_step(0.01);
        solver.bed_columns[0].z_bed
    };
    let shallow = build(None);
    let deep = build(Some(0.4));
    let w_s = settling_velocity_ferguson_church(0.5e-3, 2650.0, 1000.0, NU_WATER);
    let cf = friction_coeff_manning(0.025, 0.4);
    let tau_b = bed_shear_stress(1000.0, cf, u);
    assert!(
        tau_b < 0.1,
        "depth 0.4 m keeps the stress below tau_cd: {tau_b}"
    );
    let expected = w_s * c * (1.0 - tau_b / 0.1) * 0.01 / (1.0 - 0.4);
    assert!(
        (deep - expected).abs() <= 1.0e-12 * expected,
        "{deep} vs {expected}"
    );
    // In 5 cm of water the same current exerts more than tau_cd: no deposition.
    assert_eq!(shallow, 0.0);
}

#[test]
fn test_columns_without_nearby_fluid_are_not_touched() {
    let mut solver = SedimentSphSolver::new(SedimentSphConfig::default());
    solver.add_bed_column(0.0, 0.0, 0.2);
    // Fluid above the near-bed sampling height, and fluid far away horizontally.
    solver.add_fluid_particle([0.0, 0.0, 0.9], [1.0, 0.0, 0.0], 0.01, PARTICLE_H);
    solver.add_fluid_particle([5.0, 5.0, 0.25], [1.0, 0.0, 0.0], 0.01, PARTICLE_H);
    solver.add_bed_particle([0.0, 0.0, 0.21], 0.01, PARTICLE_H);
    solver.morphological_step(0.02);
    assert_eq!(solver.bed_columns[0].z_bed, 0.2);
    assert_eq!(solver.total_erosion(), 0.0);
    assert_eq!(solver.total_deposition(), 0.0);
    assert!((solver.time - 0.02).abs() < 1.0e-15);
}

#[test]
fn test_morfac_and_porosity_scale_the_bed_change() {
    let run = |morfac: f64, porosity: f64| {
        let cfg = SedimentSphConfig {
            morfac,
            porosity,
            ..suspended_only()
        };
        let mut solver = flat_bed_solver(cfg, &[0.0], &[[0.0; 3]], &[0.1]);
        solver.morphological_step(0.01);
        solver.bed_columns[0].z_bed
    };
    let base = run(1.0, 0.0);
    assert!(base > 0.0);
    assert!((run(5.0, 0.0) / base - 5.0).abs() < 1.0e-12);
    assert!((run(1.0, 0.4) / base - 1.0 / 0.6).abs() < 1.0e-12);
}

#[test]
fn test_disabled_processes_contribute_nothing() {
    let cfg = SedimentSphConfig {
        enable_bedload: false,
        enable_suspended: false,
        ..SedimentSphConfig::default()
    };
    let velocities = [[0.4, 0.0, 0.0], [0.8, 0.0, 0.0]];
    let mut solver = flat_bed_solver(cfg, &[0.0, COLUMN_SPACING], &velocities, &[0.3, 0.3]);
    solver.morphological_step(0.01);
    assert!(bed_levels(&solver).iter().all(|z| *z == 0.0));
}

#[test]
fn test_column_radius_controls_the_flux_stencil() {
    let xs = row_positions(3);
    let velocities = [[0.5, 0.0, 0.0], [0.6, 0.0, 0.0], [0.9, 0.0, 0.0]];
    let narrow = SedimentSphConfig {
        column_radius: 0.05,
        ..bedload_only()
    };
    let mut isolated = flat_bed_solver(narrow, &xs, &velocities, &[0.0; 3]);
    isolated.morphological_step(0.01);
    assert!(
        bed_levels(&isolated).iter().all(|z| *z == 0.0),
        "no column inside the radius means no flux divergence"
    );
    let mut coupled = flat_bed_solver(bedload_only(), &xs, &velocities, &[0.0; 3]);
    coupled.morphological_step(0.01);
    assert!(coupled.bed_columns[1].z_bed < 0.0);
}

#[test]
fn test_bucketed_fluid_search_matches_a_brute_force_scan() {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1_u64 << 53) as f64
    };
    let mut particles = Vec::new();
    for _ in 0..400 {
        let pos = [next() * 2.0 - 0.5, next() * 2.0 - 0.5, next() * 0.2];
        let mut p =
            SedimentSphParticle::new_fluid(pos, [0.0; 3], 0.01, 1000.0, 0.01 + 0.05 * next());
        p.is_bed = next() < 0.1;
        particles.push(p);
    }
    let buckets = FluidBuckets::build(&particles);
    for _ in 0..60 {
        let (cx, cy) = (next() * 2.0 - 0.5, next() * 2.0 - 0.5);
        let inside = |p: &SedimentSphParticle| {
            !p.is_bed && (p.pos[0] - cx).powi(2) + (p.pos[1] - cy).powi(2) < (2.0 * p.h).powi(2)
        };
        let mut expected: Vec<usize> = (0..particles.len())
            .filter(|&i| inside(&particles[i]))
            .collect();
        let mut found: Vec<usize> = buckets
            .candidates(cx, cy)
            .filter(|&i| inside(&particles[i]))
            .collect();
        expected.sort_unstable();
        found.sort_unstable();
        assert_eq!(found, expected, "column at ({cx}, {cy})");
    }
}

#[test]
fn test_least_squares_divergence_is_exact_for_linear_fields() {
    // q = q0 + A x with div q = A_xx + A_yy, sampled on a 3 x 3 grid.
    let (q0, a) = ([0.3, -0.2], [[2.0, 0.5], [-1.5, 4.0]]);
    let flux = |x: [f64; 2]| {
        [
            q0[0] + a[0][0] * x[0] + a[0][1] * x[1],
            q0[1] + a[1][0] * x[0] + a[1][1] * x[1],
        ]
    };
    let mut points = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let x = [0.1 * i as f64, 0.1 * j as f64];
            points.push((x, flux(x)));
        }
    }
    for (c, &centre) in points.iter().enumerate() {
        let others: Vec<_> = points
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != c)
            .map(|(_, s)| *s)
            .collect();
        let div = least_squares_divergence(centre, &others);
        assert!((div - 6.0).abs() < 1.0e-12, "centre {c}: div = {div}");
    }
}

#[test]
fn test_least_squares_divergence_on_degenerate_stencils() {
    let centre = ([0.0, 0.0], [1.0, 2.0]);
    // No neighbours.
    assert_eq!(least_squares_divergence(centre, &[]), 0.0);
    // One neighbour along x resolves dq_x/dx only.
    let div = least_squares_divergence(centre, &[([0.1, 0.0], [1.5, 7.0])]);
    assert!((div - 5.0).abs() < 1.0e-12, "{div}");
    // Neighbours on a diagonal line: q varies as (s, 2s) along s, direction (1,1)/sqrt2.
    let diag = ([0.0, 0.0], [0.0, 0.0]);
    let line = [([0.1, 0.1], [0.1, 0.2]), ([-0.1, -0.1], [-0.1, -0.2])];
    // dq/ds . s_hat where dq/dx_along = (1, 2) per unit of x = y: div = (1 + 2) / 2.
    let div = least_squares_divergence(diag, &line);
    assert!((div - 1.5).abs() < 1.0e-12, "{div}");
    // Coincident neighbours carry no geometric information.
    assert_eq!(
        least_squares_divergence(centre, &[([0.0, 0.0], [3.0, 3.0])]),
        0.0
    );
}
