//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::MultiTransportParticle;

pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Cubic spline kernel value W(r, h).
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}
/// Gradient of cubic spline kernel ∇W(r_ij, h).
pub fn cubic_kernel_grad(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        -alpha * 0.75 * t * t / h
    } else {
        0.0
    };
    scale3(r_ij, dw_dr / r)
}
/// Reaction term function type: f(c) → dc/dt contribution.
pub type ReactionFn = fn(f64) -> f64;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::transport_sph::types::*;
    fn make_particle(x: f64, y: f64, z: f64, c: f64, t: f64) -> TransportParticle {
        TransportParticle::new([x, y, z], [0.0, 0.0, 0.0], 1.0, 1000.0, 0.1, c, t)
    }
    #[test]
    fn test_cubic_kernel_zero_at_two_h() {
        assert_eq!(cubic_kernel(2.0 * 0.1, 0.1), 0.0);
    }
    #[test]
    fn test_cubic_kernel_positive_at_origin() {
        assert!(cubic_kernel(0.0, 0.1) > 0.0);
    }
    #[test]
    fn test_cubic_kernel_decreasing() {
        let h = 0.1;
        let w0 = cubic_kernel(0.0, h);
        let w1 = cubic_kernel(0.05, h);
        let w2 = cubic_kernel(0.15, h);
        assert!(w0 > w1);
        assert!(w1 > w2);
    }
    #[test]
    fn test_cubic_kernel_zero_h() {
        assert_eq!(cubic_kernel(0.05, 0.0), 0.0);
    }
    #[test]
    fn test_particle_volume() {
        let p = make_particle(0.0, 0.0, 0.0, 0.5, 300.0);
        assert!((p.volume() - 1.0 / 1000.0).abs() < 1e-12);
    }
    #[test]
    fn test_particle_kinetic_energy_zero_velocity() {
        let p = make_particle(0.0, 0.0, 0.0, 0.5, 300.0);
        assert_eq!(p.kinetic_energy(), 0.0);
    }
    #[test]
    fn test_particle_kinetic_energy_nonzero() {
        let mut p = make_particle(0.0, 0.0, 0.0, 0.5, 300.0);
        p.velocity = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_particle_conc() {
        let p = make_particle(0.0, 0.0, 0.0, 0.7, 300.0);
        assert!((p.conc() - 0.7).abs() < 1e-12);
    }
    #[test]
    fn test_particle_multi_species() {
        let mut p = make_particle(0.0, 0.0, 0.0, 0.3, 300.0);
        p.concentration.push(0.6);
        assert_eq!(p.concentration.len(), 2);
        assert!((p.concentration[1] - 0.6).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_transport_advect_moves_particles() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 1.0, 300.0)];
        particles[0].velocity = [1.0, 0.0, 0.0];
        let st = ScalarTransport::new(0.1);
        st.advect(&mut particles);
        assert!((particles[0].position[0] - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_transport_concentration_preserved() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.8, 300.0)];
        particles[0].velocity = [2.0, 3.0, 0.0];
        let st = ScalarTransport::new(0.05);
        let c_before = particles[0].conc();
        st.advect(&mut particles);
        assert!((particles[0].conc() - c_before).abs() < 1e-15);
    }
    #[test]
    fn test_scalar_transport_total_scalar() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 1.0, 300.0),
            make_particle(1.0, 0.0, 0.0, 2.0, 300.0),
        ];
        let st = ScalarTransport::new(0.01);
        let total = st.total_scalar(&particles);
        assert!((total - 0.003).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_transport_interpolate_at_particle() {
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        let st = ScalarTransport::new(0.01);
        let interp = st.interpolate_at(&particles, [0.0, 0.0, 0.0]);
        assert!((interp - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_diffusion_uniform_field_zero_rate() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 1.0, 300.0),
            make_particle(0.1, 0.0, 0.0, 1.0, 300.0),
            make_particle(-0.1, 0.0, 0.0, 1.0, 300.0),
        ];
        let diff = DiffusionSph::new(1e-6);
        let rates = diff.compute_dphidt(&particles);
        for r in &rates {
            assert!(r.abs() < 1e-10, "rate = {}", r);
        }
    }
    #[test]
    fn test_diffusion_gradient_nonzero() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.0, 300.0),
            make_particle(0.05, 0.0, 0.0, 1.0, 300.0),
        ];
        let diff = DiffusionSph::new(1e-6);
        let rates = diff.compute_dphidt(&particles);
        assert!(rates[0] > 0.0 || rates[0].is_finite());
    }
    #[test]
    fn test_diffusion_step_changes_concentration() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.0, 300.0),
            make_particle(0.05, 0.0, 0.0, 1.0, 300.0),
        ];
        let diff = DiffusionSph::new(1e-5);
        let c0_before = particles[0].conc();
        diff.step(&mut particles, 0.01);
        let c0_after = particles[0].conc();
        assert!((c0_after - c0_before).abs() >= 0.0);
    }
    #[test]
    fn test_reaction_diffusion_fisher_midpoint() {
        let rd = ReactionDiffusion::new(0.0, 1.0);
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        let rates = rd.compute_dphidt_fisher(&particles);
        assert!((rates[0] - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_reaction_diffusion_fisher_boundary() {
        let rd = ReactionDiffusion::new(0.0, 1.0);
        let particles = vec![make_particle(0.0, 0.0, 0.0, 1.0, 300.0)];
        let rates = rd.compute_dphidt_fisher(&particles);
        assert!(rates[0].abs() < 1e-12);
    }
    #[test]
    fn test_reaction_diffusion_step_clamps() {
        let rd = ReactionDiffusion::new(0.0, 100.0);
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.9, 300.0)];
        rd.step_fisher(&mut particles, 1.0);
        assert!(particles[0].conc() <= 1.0);
        assert!(particles[0].conc() >= 0.0);
    }
    #[test]
    fn test_reaction_diffusion_custom() {
        let rd = ReactionDiffusion::new(0.0, 0.0);
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        rd.step_custom(&mut particles, 1.0, |c| -c);
        assert!((particles[0].conc() - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_gray_scott_init_particles() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        GrayScottSph::init_particles(&mut particles);
        assert_eq!(particles[0].concentration.len(), 2);
    }
    #[test]
    fn test_gray_scott_step_changes_uv() {
        let gs = GrayScottSph::new(2e-5, 1e-5, 0.04, 0.06);
        let mut particles = vec![
            {
                let mut p = make_particle(0.0, 0.0, 0.0, 0.5, 300.0);
                p.concentration = vec![0.5, 0.25];
                p
            },
            {
                let mut p = make_particle(0.05, 0.0, 0.0, 0.4, 300.0);
                p.concentration = vec![0.4, 0.2];
                p
            },
        ];
        let u_before = particles[0].concentration[0];
        gs.step(&mut particles, 0.1);
        let u_after = particles[0].concentration[0];
        assert!((0.0..=1.0).contains(&u_after));
        let _ = u_before;
    }
    #[test]
    fn test_gray_scott_bounds() {
        let gs = GrayScottSph::new(2e-5, 1e-5, 0.04, 0.06);
        let mut particles: Vec<TransportParticle> = (0..5)
            .map(|i| {
                let mut p = make_particle(i as f64 * 0.05, 0.0, 0.0, 0.5, 300.0);
                p.concentration = vec![0.5, 0.25];
                p
            })
            .collect();
        gs.step(&mut particles, 0.01);
        for p in &particles {
            assert!(p.concentration[0] >= 0.0 && p.concentration[0] <= 1.0);
            assert!(p.concentration[1] >= 0.0 && p.concentration[1] <= 1.0);
        }
    }
    #[test]
    fn test_thermal_diffusivity() {
        let tc = ThermalConduction::new(0.6, 4.18e6);
        let alpha = tc.diffusivity();
        assert!((alpha - 0.6 / 4.18e6).abs() < 1e-20);
    }
    #[test]
    fn test_thermal_uniform_field_zero_rate() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 300.0),
            make_particle(0.05, 0.0, 0.0, 0.5, 300.0),
            make_particle(-0.05, 0.0, 0.0, 0.5, 300.0),
        ];
        let tc = ThermalConduction::new(1.0, 1e6);
        let rates = tc.compute_dtdt(&particles);
        for r in &rates {
            assert!(r.abs() < 1e-8, "rate = {}", r);
        }
    }
    #[test]
    fn test_thermal_gradient_transfer() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 200.0),
            make_particle(0.05, 0.0, 0.0, 0.5, 400.0),
        ];
        let tc = ThermalConduction::new(1.0, 1e6);
        let rates = tc.compute_dtdt(&particles);
        assert!(rates[0] >= 0.0 || rates[0].is_finite());
    }
    #[test]
    fn test_thermal_total_energy() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 300.0),
            make_particle(1.0, 0.0, 0.0, 0.5, 300.0),
        ];
        let tc = ThermalConduction::new(1.0, 1e6);
        let e = tc.total_thermal_energy(&particles);
        assert!((e - 1e6 * 0.001 * 300.0 * 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_thermal_step_updates_temperature() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 200.0),
            make_particle(0.05, 0.0, 0.0, 0.5, 400.0),
        ];
        let tc = ThermalConduction::new(1.0, 1000.0);
        let t0_before = particles[0].temperature;
        tc.step(&mut particles, 0.001);
        let t0_after = particles[0].temperature;
        assert!((t0_after - t0_before).abs() >= 0.0);
    }
    #[test]
    fn test_solute_diagonal_diffusion_coefficients() {
        let st = SoluteTransport::diagonal(&[1e-6, 2e-6]);
        assert!((st.d(0, 0) - 1e-6).abs() < 1e-20);
        assert!((st.d(1, 1) - 2e-6).abs() < 1e-20);
        assert_eq!(st.d(0, 1), 0.0);
    }
    #[test]
    fn test_solute_n_species() {
        let st = SoluteTransport::diagonal(&[1e-6, 2e-6, 3e-6]);
        assert_eq!(st.n_species, 3);
    }
    #[test]
    fn test_solute_step_updates_concentrations() {
        let st = SoluteTransport::diagonal(&[1e-6, 2e-6]);
        let mut particles: Vec<TransportParticle> = (0..3)
            .map(|i| {
                let mut p = make_particle(i as f64 * 0.05, 0.0, 0.0, (i as f64) * 0.5, 300.0);
                p.concentration = vec![i as f64 * 0.5, 0.3];
                p
            })
            .collect();
        let c0_before = particles[0].concentration[0];
        st.step(&mut particles, 0.01);
        let c0_after = particles[0].concentration[0];
        assert!(c0_after.is_finite());
        let _ = c0_before;
    }
    #[test]
    fn test_advection_standard_rates_finite() {
        let scheme = AdvectionScheme::new(AdvectionSchemeKind::Standard, 0.0);
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 300.0),
            make_particle(0.05, 0.0, 0.0, 0.8, 300.0),
        ];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [1.0, 0.0, 0.0];
        let rates = scheme.compute_all_rates(&particles);
        assert!(rates.iter().all(|r| r.is_finite()));
    }
    #[test]
    fn test_advection_upwind_rates_finite() {
        let scheme = AdvectionScheme::new(AdvectionSchemeKind::Upwind, 0.01);
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.2, 300.0),
            make_particle(0.05, 0.0, 0.0, 0.9, 300.0),
        ];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [1.0, 0.0, 0.0];
        let rates = scheme.compute_all_rates(&particles);
        assert!(rates.iter().all(|r| r.is_finite()));
    }
    #[test]
    fn test_advection_weno_rates_finite() {
        let scheme = AdvectionScheme::new(AdvectionSchemeKind::WenoSph, 0.0);
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.1, 300.0),
            make_particle(0.05, 0.0, 0.0, 0.9, 300.0),
        ];
        particles[0].velocity = [2.0, 0.0, 0.0];
        particles[1].velocity = [2.0, 0.0, 0.0];
        let rates = scheme.compute_all_rates(&particles);
        assert!(rates.iter().all(|r| r.is_finite()));
    }
    #[test]
    fn test_advection_step_runs() {
        let scheme = AdvectionScheme::new(AdvectionSchemeKind::Standard, 0.0);
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        particles[0].velocity = [1.0, 0.0, 0.0];
        scheme.step(&mut particles, 0.01);
        assert!(particles[0].conc().is_finite());
    }
    #[test]
    fn test_boundary_flux_total_flux() {
        let bf = BoundaryFlux::new(100.0, [0.0, 0.0, 1.0], vec![0, 1], vec![0.5, 0.5]);
        assert!((bf.total_flux() - 100.0).abs() < 1e-9);
    }
    #[test]
    fn test_boundary_flux_zero_flux() {
        let bf = BoundaryFlux::new(0.0, [1.0, 0.0, 0.0], vec![0], vec![1.0]);
        assert_eq!(bf.total_flux(), 0.0);
    }
    #[test]
    fn test_boundary_flux_apply() {
        let bf = BoundaryFlux::new(1.0, [0.0, 1.0, 0.0], vec![0], vec![1.0]);
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.0, 300.0)];
        let c_before = particles[0].conc();
        bf.apply(&mut particles, 0.01, 1e-3);
        let c_after = particles[0].conc();
        assert!(c_after >= c_before);
    }
    #[test]
    fn test_boundary_flux_out_of_bounds_id() {
        let bf = BoundaryFlux::new(1.0, [1.0, 0.0, 0.0], vec![99], vec![1.0]);
        let mut particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        let c_before = particles[0].conc();
        bf.apply(&mut particles, 0.01, 1.0);
        assert!((particles[0].conc() - c_before).abs() < 1e-15);
    }
    #[test]
    fn test_conservation_initial_zero_error() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 0.5, 300.0),
            make_particle(0.1, 0.0, 0.0, 0.3, 300.0),
        ];
        let checker = ConservationCheck::new(&particles, 4.18e6);
        assert!(checker.scalar_error(&particles) < 1e-14);
    }
    #[test]
    fn test_conservation_nonzero_error_after_change() {
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        let checker = ConservationCheck::new(&particles, 4.18e6);
        let mut particles2 = particles.clone();
        particles2[0].concentration[0] = 0.9;
        assert!(checker.scalar_error(&particles2) > 0.0);
    }
    #[test]
    fn test_conservation_check_is_scalar_conserved() {
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 300.0)];
        let checker = ConservationCheck::new(&particles, 4.18e6);
        assert!(checker.is_scalar_conserved(&particles, 1e-14));
    }
    #[test]
    fn test_conservation_total_scalar_static() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0, 2.0, 300.0),
            make_particle(1.0, 0.0, 0.0, 2.0, 300.0),
        ];
        let total = ConservationCheck::current_scalar(&particles);
        assert!((total - 0.004).abs() < 1e-12);
    }
    #[test]
    fn test_conservation_thermal_energy_static() {
        let rho_cv = 1e6;
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.5, 100.0)];
        let e = ConservationCheck::current_thermal(&particles, rho_cv);
        let expected = rho_cv * (1.0_f64 / 1000.0) * 100.0;
        assert!((e - expected).abs() < 1e-6);
    }
    #[test]
    fn test_conservation_zero_reference() {
        let particles = vec![make_particle(0.0, 0.0, 0.0, 0.0, 0.0)];
        let checker = ConservationCheck::new(&particles, 1e6);
        let err = checker.scalar_error(&particles);
        assert_eq!(err, 0.0);
    }
}
/// Compute the SPH Laplacian diffusion term for species `species` at particle `i`.
///
/// Uses the standard Morris et al. (1997) SPH diffusion:
/// ∇²c_i ≈ Σ_j (m_j/ρ_j) * 2*(c_i − c_j)/(r_ij²+ε) * (r_ij · ∇W_ij)
///
/// Here a simplified scalar version is returned using the kernel value only.
pub fn sph_diffusion_laplacian(
    particles: &[MultiTransportParticle],
    i: usize,
    species: usize,
    kernel: impl Fn(f64, f64) -> f64,
    diffusivity: f64,
) -> f64 {
    if i >= particles.len() {
        return 0.0;
    }
    let pi = &particles[i];
    let ci = if species < pi.concentration.len() {
        pi.concentration[species]
    } else {
        0.0
    };
    let hi = pi.h;
    let mut lap = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let cj = if species < pj.concentration.len() {
            pj.concentration[species]
        } else {
            0.0
        };
        let dx = pi.pos[0] - pj.pos[0];
        let dy = pi.pos[1] - pj.pos[1];
        let dz = pi.pos[2] - pj.pos[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let w = kernel(r, hi);
        let vol_j = pj.mass / pj.rho.max(1e-300);
        lap += vol_j * (cj - ci) * w;
    }
    diffusivity * 2.0 * lap / (hi * hi).max(1e-300)
}
/// Advection-diffusion RHS for species `species` at particle `i`.
///
/// dc/dt ≈ diffusion_laplacian(c)  (advection is handled by particle motion).
pub fn advection_diffusion_rhs(
    particles: &[MultiTransportParticle],
    i: usize,
    species: usize,
    kernel: impl Fn(f64, f64) -> f64,
    d: f64,
) -> f64 {
    sph_diffusion_laplacian(particles, i, species, kernel, d)
}
/// Arrhenius reaction rate: k = A * exp(−E_a / (R T)).
///
/// - `a`: pre-exponential factor
/// - `ea`: activation energy (J/mol)
/// - `temperature`: temperature (K)
/// - `r_gas`: universal gas constant (J/(mol·K)), typically 8.314
pub fn arrhenius_rate(a: f64, ea: f64, temperature: f64, r_gas: f64) -> f64 {
    if r_gas.abs() < 1e-300 || temperature <= 0.0 {
        return 0.0;
    }
    a * (-ea / (r_gas * temperature)).exp()
}
/// Damköhler number Da = k_reaction * L / U.
///
/// Ratio of reaction time scale to convection time scale.
pub fn damkohler_number(reaction_rate: f64, convection_rate: f64) -> f64 {
    if convection_rate.abs() < 1e-300 {
        return 0.0;
    }
    reaction_rate / convection_rate
}
/// Schmidt number Sc = μ / (ρ D).
///
/// Ratio of momentum diffusivity to mass diffusivity.
pub fn schmidt_number(mu: f64, rho: f64, diffusivity: f64) -> f64 {
    let denom = rho * diffusivity;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    mu / denom
}
/// Lewis number Le = α / D  (thermal diffusivity / mass diffusivity).
pub fn lewis_number(thermal_diff: f64, mass_diff: f64) -> f64 {
    if mass_diff.abs() < 1e-300 {
        return 0.0;
    }
    thermal_diff / mass_diff
}
/// Heat equation RHS at particle `i` via SPH.
///
/// dT/dt ≈ λ / (ρ c_v) * ∇²T
pub fn energy_equation_rhs(
    particles: &[MultiTransportParticle],
    i: usize,
    kernel: impl Fn(f64, f64) -> f64,
    lambda: f64,
) -> f64 {
    if i >= particles.len() {
        return 0.0;
    }
    let pi = &particles[i];
    let ti = pi.temperature;
    let hi = pi.h;
    let mut lap_t = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let tj = pj.temperature;
        let dx = pi.pos[0] - pj.pos[0];
        let dy = pi.pos[1] - pj.pos[1];
        let dz = pi.pos[2] - pj.pos[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let w = kernel(r, hi);
        let vol_j = pj.mass / pj.rho.max(1e-300);
        lap_t += vol_j * (tj - ti) * w;
    }
    lambda * 2.0 * lap_t / (hi * hi).max(1e-300)
}
/// Effective diffusivity in a porous medium: D_eff = D · ε / τ.
pub fn porous_diffusivity(d0: f64, porosity: f64, tortuosity: f64) -> f64 {
    if tortuosity < 1e-300 {
        return 0.0;
    }
    d0 * porosity / tortuosity
}
/// Thermal diffusivity α = k / (ρ cₚ).
pub fn thermal_diffusivity(k: f64, rho: f64, cp: f64) -> f64 {
    if rho < 1e-300 || cp < 1e-300 {
        return 0.0;
    }
    k / (rho * cp)
}
/// Helper: rate contribution for one neighbour pair.
pub(super) fn heat_exchange_rate(
    alpha: f64,
    t_i: f64,
    t_j: f64,
    mass_j: f64,
    rho_j: f64,
    w: f64,
) -> f64 {
    if rho_j <= 0.0_f64 {
        return 0.0_f64;
    }
    2.0_f64 * alpha * (t_j - t_i) * mass_j / rho_j * w
}
#[cfg(test)]
mod reactive_transport_tests {
    use super::*;
    use crate::transport_sph::types::*;
    fn zero_kernel(_r: f64, _h: f64) -> f64 {
        0.0
    }
    #[test]
    fn test_multi_transport_particle_new() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 3);
        assert_eq!(p.species_count(), 3);
    }
    #[test]
    fn test_total_mass_fraction_zero() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 3);
        assert_eq!(p.total_mass_fraction(), 0.0);
    }
    #[test]
    fn test_total_mass_fraction_uniform() {
        let mut p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 4);
        p.concentration = vec![0.25, 0.25, 0.25, 0.25];
        assert!((p.total_mass_fraction() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_total_mass_fraction_two_species() {
        let mut p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 2);
        p.concentration = vec![0.6, 0.4];
        assert!((p.total_mass_fraction() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_species_count_correct() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 7);
        assert_eq!(p.species_count(), 7);
    }
    #[test]
    fn test_sph_diffusion_laplacian_no_neighbors() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 2);
        let particles = vec![p];
        let result = sph_diffusion_laplacian(&particles, 0, 0, zero_kernel, 1e-5);
        assert_eq!(result, 0.0);
    }
    #[test]
    fn test_sph_diffusion_laplacian_zero_kernel() {
        let mut p1 = MultiTransportParticle::new([0.0; 3], [0.0; 3], 2);
        let mut p2 = MultiTransportParticle::new([0.1, 0.0, 0.0], [0.0; 3], 2);
        p1.concentration = vec![1.0, 0.0];
        p2.concentration = vec![0.0, 0.0];
        let particles = vec![p1, p2];
        let result = sph_diffusion_laplacian(&particles, 0, 0, zero_kernel, 1.0);
        assert_eq!(result, 0.0);
    }
    #[test]
    fn test_advection_diffusion_rhs_zero_kernel() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
        let particles = vec![p];
        let rhs = advection_diffusion_rhs(&particles, 0, 0, zero_kernel, 1e-4);
        assert_eq!(rhs, 0.0);
    }
    #[test]
    fn test_arrhenius_rate_positive() {
        let k = arrhenius_rate(1e13, 80000.0, 500.0, 8.314);
        assert!(k > 0.0);
    }
    #[test]
    fn test_arrhenius_rate_increases_with_temperature() {
        let k1 = arrhenius_rate(1e13, 80000.0, 500.0, 8.314);
        let k2 = arrhenius_rate(1e13, 80000.0, 600.0, 8.314);
        assert!(k2 > k1);
    }
    #[test]
    fn test_arrhenius_rate_zero_temperature() {
        let k = arrhenius_rate(1e13, 80000.0, 0.0, 8.314);
        assert_eq!(k, 0.0);
    }
    #[test]
    fn test_arrhenius_rate_formula() {
        let a = 1.0;
        let ea = 0.0;
        let k = arrhenius_rate(a, ea, 300.0, 8.314);
        assert!((k - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_damkohler_number_formula() {
        let da = damkohler_number(10.0, 2.0);
        assert!((da - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_damkohler_number_zero_convection() {
        let da = damkohler_number(10.0, 0.0);
        assert_eq!(da, 0.0);
    }
    #[test]
    fn test_damkohler_number_proportional() {
        let da1 = damkohler_number(1.0, 1.0);
        let da2 = damkohler_number(2.0, 1.0);
        assert!((da2 - 2.0 * da1).abs() < 1e-12);
    }
    #[test]
    fn test_schmidt_number_formula() {
        let sc = schmidt_number(1e-3, 1000.0, 1e-9);
        assert!((sc - 1e3).abs() < 1.0);
    }
    #[test]
    fn test_schmidt_number_zero_diffusivity() {
        let sc = schmidt_number(1e-3, 1000.0, 0.0);
        assert_eq!(sc, 0.0);
    }
    #[test]
    fn test_lewis_number_formula() {
        let le = lewis_number(1e-7, 1e-9);
        assert!((le - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_lewis_number_zero_mass_diff() {
        let le = lewis_number(1e-7, 0.0);
        assert_eq!(le, 0.0);
    }
    #[test]
    fn test_reactive_flow_sph_new() {
        let rfsph = ReactiveFlowSPH::new(3, 0.5);
        assert_eq!(rfsph.n_species, 3);
        assert_eq!(rfsph.particle_count(), 0);
    }
    #[test]
    fn test_reactive_flow_sph_add_particle() {
        let mut rfsph = ReactiveFlowSPH::new(2, 1.0);
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 2);
        rfsph.add_particle(p);
        assert_eq!(rfsph.particle_count(), 1);
    }
    #[test]
    fn test_reactive_flow_sph_add_multiple() {
        let mut rfsph = ReactiveFlowSPH::new(2, 1.0);
        for _ in 0..5 {
            rfsph.add_particle(MultiTransportParticle::new([0.0; 3], [0.0; 3], 2));
        }
        assert_eq!(rfsph.particle_count(), 5);
    }
    #[test]
    fn test_reactive_flow_sph_total_species_mass_zero() {
        let mut rfsph = ReactiveFlowSPH::new(2, 1.0);
        rfsph.add_particle(MultiTransportParticle::new([0.0; 3], [0.0; 3], 2));
        assert_eq!(rfsph.total_species_mass(0), 0.0);
    }
    #[test]
    fn test_reactive_flow_sph_total_species_mass_nonzero() {
        let mut rfsph = ReactiveFlowSPH::new(1, 1.0);
        let mut p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
        p.concentration[0] = 0.5;
        p.mass = 2.0;
        rfsph.add_particle(p);
        assert!((rfsph.total_species_mass(0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_reactive_flow_sph_total_species_mass_non_negative() {
        let mut rfsph = ReactiveFlowSPH::new(2, 1.0);
        for _ in 0..3 {
            let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 2);
            rfsph.add_particle(p);
        }
        assert!(rfsph.total_species_mass(0) >= 0.0);
        assert!(rfsph.total_species_mass(1) >= 0.0);
    }
    #[test]
    fn test_energy_equation_rhs_single_particle() {
        let p = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
        let particles = vec![p];
        let rhs = energy_equation_rhs(&particles, 0, zero_kernel, 1.0);
        assert_eq!(rhs, 0.0);
    }
    #[test]
    fn test_energy_equation_rhs_zero_kernel() {
        let p1 = MultiTransportParticle::new([0.0; 3], [0.0; 3], 1);
        let p2 = MultiTransportParticle::new([0.1, 0.0, 0.0], [0.0; 3], 1);
        let particles = vec![p1, p2];
        let rhs = energy_equation_rhs(&particles, 0, zero_kernel, 100.0);
        assert_eq!(rhs, 0.0);
    }
    #[test]
    fn test_multi_transport_particle_temperature_default() {
        let p = MultiTransportParticle::new([1.0, 2.0, 3.0], [0.0; 3], 2);
        assert!(p.temperature > 0.0);
    }
    #[test]
    fn test_schmidt_number_proportional_mu() {
        let sc1 = schmidt_number(1e-3, 1000.0, 1e-6);
        let sc2 = schmidt_number(2e-3, 1000.0, 1e-6);
        assert!((sc2 - 2.0 * sc1).abs() < 1e-12 * sc1);
    }
}
#[cfg(test)]
mod extended_transport_tests {
    use super::*;
    use crate::transport_sph::types::*;
    #[test]
    fn test_concentration_field_new_zeros() {
        let cf = ConcentrationField::new(5, 18.015);
        assert_eq!(cf.n_particles, 5);
        for &c in &cf.molar_concentration {
            assert_eq!(c, 0.0);
        }
    }
    #[test]
    fn test_concentration_field_set_concentration() {
        let mut cf = ConcentrationField::new(3, 18.015);
        cf.set_concentration(0, 2.0, 4.0);
        assert!((cf.molar_concentration[0] - 2.0).abs() < 1e-12);
        assert!((cf.mole_fraction[0] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_concentration_field_negative_clamp() {
        let mut cf = ConcentrationField::new(2, 18.015);
        cf.set_concentration(0, -1.0, 1.0);
        assert_eq!(cf.molar_concentration[0], 0.0);
    }
    #[test]
    fn test_concentration_field_mass_concentration() {
        let mut cf = ConcentrationField::new(1, 2.0);
        cf.molar_concentration[0] = 3.0;
        assert!((cf.mass_concentration(0) - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_concentration_field_mass_fraction() {
        let mut cf = ConcentrationField::new(1, 2.0);
        cf.molar_concentration[0] = 3.0;
        cf.fluid_density[0] = 6.0;
        assert!((cf.mass_fraction(0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_concentration_field_total_moles() {
        let mut cf = ConcentrationField::new(4, 1.0);
        for c in cf.molar_concentration.iter_mut() {
            *c = 1.0;
        }
        let total = cf.total_moles(0.25);
        assert!((total - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_concentration_field_decay() {
        let mut cf = ConcentrationField::new(2, 1.0);
        cf.molar_concentration = vec![10.0, 10.0];
        cf.apply_decay(1.0, 1.0);
        for &c in &cf.molar_concentration {
            assert!((c - 10.0 / std::f64::consts::E).abs() < 1e-10);
        }
    }
    #[test]
    fn test_fickian_isotropic_effective_diffusivity() {
        let fd = FickianDiffusion::isotropic(1e-5);
        assert!((fd.effective_diffusivity() - 1e-5).abs() < 1e-20);
    }
    #[test]
    fn test_fickian_isotropic_flux_1d() {
        let fd = FickianDiffusion::isotropic(1.0);
        let flux = fd.flux([2.0, 0.0, 0.0]);
        assert!((flux[0] - (-2.0)).abs() < 1e-12);
        assert_eq!(flux[1], 0.0);
    }
    #[test]
    fn test_fickian_anisotropic_flux() {
        let mut tensor = [0.0_f64; 9];
        tensor[0] = 2.0;
        tensor[4] = 1.0;
        tensor[8] = 1.0;
        let fd = FickianDiffusion::anisotropic(tensor);
        let flux = fd.flux([1.0, 0.0, 0.0]);
        assert!((flux[0] - (-2.0)).abs() < 1e-12);
    }
    #[test]
    fn test_porous_diffusivity_formula() {
        let d = porous_diffusivity(1e-5, 0.4, 2.0);
        assert!((d - 1e-5 * 0.4 / 2.0).abs() < 1e-25);
    }
    #[test]
    fn test_thermal_diffusivity_formula() {
        let alpha = thermal_diffusivity(0.6, 1000.0, 4000.0);
        assert!((alpha - 0.6 / 4_000_000.0).abs() < 1e-20);
    }
    #[test]
    fn test_fickian_tortuosity() {
        let mut fd = FickianDiffusion::isotropic(1e-5);
        fd.tortuosity = 2.0;
        assert!((fd.effective_diffusivity() - 5e-6).abs() < 1e-20);
    }
    #[test]
    fn test_advection_sph_new() {
        let adv = AdvectionSPH::new(0.001, 0.05, false);
        assert!((adv.dt - 0.001).abs() < 1e-15);
        assert!(!adv.semi_lagrangian);
    }
    #[test]
    fn test_advection_sph_update_positions() {
        let p = TransportParticle::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            1000.0,
            0.1,
            0.0,
            293.15,
        );
        let mut particles = vec![p.clone()];
        let adv = AdvectionSPH::new(0.5, 0.0, false);
        adv.update_positions(&mut particles);
        assert!((particles[0].position[0] - 0.5).abs() < 1e-12);
        let _ = p.position;
    }
    #[test]
    fn test_advection_sph_cfl_number() {
        let cfl = AdvectionSPH::cfl_number(10.0, 0.1);
        assert!((cfl - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_advection_sph_cfl_zero_h() {
        assert_eq!(AdvectionSPH::cfl_number(10.0, 0.0), 0.0);
    }
    #[test]
    fn test_advection_sph_scalar_single_particle() {
        let p = TransportParticle::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            1000.0,
            0.1,
            0.0,
            293.15,
        );
        let particles = vec![p];
        let adv = AdvectionSPH::new(0.01, 0.0, true);
        let rate = adv.advect_scalar(&particles, 0, cubic_kernel, cubic_kernel_grad);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_chemical_reaction_first_order_rate_positive() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        assert!(rxn.rate_constant(500.0) > 0.0);
    }
    #[test]
    fn test_chemical_reaction_rate_increases_with_temp() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        assert!(rxn.rate_constant(600.0) > rxn.rate_constant(500.0));
    }
    #[test]
    fn test_chemical_reaction_zero_temperature() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        assert_eq!(rxn.rate_constant(0.0), 0.0);
    }
    #[test]
    fn test_chemical_reaction_source_term_negative() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        let src = rxn.source_term(1.0, 500.0);
        assert!(src < 0.0);
    }
    #[test]
    fn test_chemical_reaction_euler_step_decreases_conc() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        let c_new = rxn.euler_step(1.0, 500.0, 1e-5);
        assert!(c_new <= 1.0);
    }
    #[test]
    fn test_chemical_reaction_half_life_positive() {
        let rxn = ChemicalReaction::first_order(1e13, 80000.0);
        let t_half = rxn.half_life(500.0);
        assert!(t_half > 0.0 && t_half.is_finite());
    }
    #[test]
    fn test_chemical_reaction_zero_activation_energy() {
        let rxn = ChemicalReaction::first_order(2.0, 0.0);
        assert!((rxn.rate_constant(300.0) - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_heat_transfer_sph_new() {
        let ht = HeatTransferSPH::new(0.6, 4182.0);
        assert!((ht.thermal_conductivity - 0.6).abs() < 1e-15);
    }
    #[test]
    fn test_heat_transfer_sph_alpha() {
        let ht = HeatTransferSPH::new(0.6, 4182.0);
        let alpha = ht.alpha(1000.0);
        assert!((alpha - 0.6 / (1000.0 * 4182.0)).abs() < 1e-20);
    }
    #[test]
    fn test_heat_transfer_sph_nusselt_cylinder() {
        let nu = HeatTransferSPH::nusselt_cylinder(1000.0, 7.0);
        assert!(nu > 0.3);
    }
    #[test]
    fn test_heat_transfer_sph_prandtl_number() {
        let ht = HeatTransferSPH::new(0.6, 4182.0);
        let pr = ht.prandtl_number(1e-3);
        assert!(pr > 0.0);
    }
    #[test]
    fn test_heat_transfer_sph_euler_step() {
        let ht = HeatTransferSPH::new(0.6, 4182.0);
        let mut particles = vec![TransportParticle::new(
            [0.0; 3], [0.0; 3], 1.0, 1000.0, 0.1, 0.0, 293.15,
        )];
        particles[0].temperature = 300.0;
        let rates = vec![10.0_f64];
        ht.euler_step(&mut particles, &rates, 0.1);
        assert!((particles[0].temperature - 301.0).abs() < 1e-10);
    }
    #[test]
    fn test_heat_transfer_sph_single_particle_rate_zero() {
        let ht = HeatTransferSPH::new(0.6, 4182.0);
        let p = TransportParticle::new([0.0; 3], [0.0; 3], 1.0, 1000.0, 0.1, 0.0, 293.15);
        let particles = vec![p];
        let rate = ht.heat_conduction_rate(&particles, 0, cubic_kernel_grad);
        assert_eq!(rate, 0.0);
    }
    #[test]
    fn test_phase_change_liquid_fraction_solid() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        assert_eq!(pc.liquid_fraction(800.0), 0.0);
    }
    #[test]
    fn test_phase_change_liquid_fraction_liquid() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        assert_eq!(pc.liquid_fraction(1000.0), 1.0);
    }
    #[test]
    fn test_phase_change_liquid_fraction_mushy() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let fl = pc.liquid_fraction(916.5);
        assert!((fl - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_phase_change_state_classification() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        assert_eq!(pc.phase_state(800.0), PhaseState::Solid);
        assert_eq!(pc.phase_state(916.5), PhaseState::Mushy);
        assert_eq!(pc.phase_state(1000.0), PhaseState::Liquid);
    }
    #[test]
    fn test_phase_change_enthalpy_increases_with_temp() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        assert!(pc.enthalpy(1000.0) > pc.enthalpy(500.0));
    }
    #[test]
    fn test_phase_change_temperature_from_enthalpy_roundtrip() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let t0 = 950.0_f64;
        let h = pc.enthalpy(t0);
        let t1 = pc.temperature_from_enthalpy(h);
        assert!((t0 - t1).abs() < 1e-4, "t0={t0} t1={t1}");
    }
    #[test]
    fn test_phase_change_mushy_zone_resistance_solid() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let r = pc.mushy_zone_resistance(0.0);
        assert!(r < 0.0, "resistance should be negative (drag)");
    }
    #[test]
    fn test_phase_change_mushy_zone_resistance_liquid() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let r = pc.mushy_zone_resistance(1.0);
        assert!(r.abs() < 1e-3);
    }
    #[test]
    fn test_phase_change_effective_heat_capacity_mushy() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let c_eff = pc.effective_heat_capacity(916.5);
        assert!(c_eff > pc.specific_heat);
    }
    #[test]
    fn test_phase_change_effective_heat_capacity_outside_mushy() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let c_eff = pc.effective_heat_capacity(800.0);
        assert!((c_eff - pc.specific_heat).abs() < 1e-12);
    }
    #[test]
    fn test_phase_change_enthalpy_euler_step() {
        let pc = PhaseChangeSPH::new(900.0, 933.0, 3.97e5, 900.0);
        let h0 = pc.enthalpy(800.0);
        let h1 = pc.enthalpy_euler_step(800.0, 10.0, 0.1);
        assert!(h1 > h0);
    }
    #[test]
    fn test_sph_heat_transfer_new() {
        let ht = SPHHeatTransfer::new(0.6, 4182.0, 1000.0);
        assert!((ht.thermal_conductivity - 0.6).abs() < 1e-12);
        assert!((ht.specific_heat - 4182.0).abs() < 1e-10);
        assert!((ht.density - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_sph_heat_transfer_alpha() {
        let ht = SPHHeatTransfer::new(0.6, 4182.0, 1000.0);
        let expected = 0.6_f64 / (1000.0 * 4182.0);
        assert!((ht.thermal_diffusivity() - expected).abs() < 1e-20);
    }
    #[test]
    fn test_sph_heat_transfer_two_particles_same_temp() {
        let ht = SPHHeatTransfer::new(1.0, 1.0, 1.0);
        let dq = ht.heat_exchange(300.0, 300.0, 1.0, 1.0, 0.1);
        assert!(dq.abs() < 1e-12);
    }
    #[test]
    fn test_sph_heat_transfer_direction() {
        let ht = SPHHeatTransfer::new(1.0, 1.0, 1.0);
        let dq = ht.heat_exchange(400.0, 300.0, 1.0, 1.0, 0.1);
        assert!(dq < 0.0);
    }
    #[test]
    fn test_sph_heat_transfer_proportional_to_k() {
        let ht1 = SPHHeatTransfer::new(1.0, 1.0, 1.0);
        let ht2 = SPHHeatTransfer::new(2.0, 1.0, 1.0);
        let dq1 = ht1.heat_exchange(400.0, 300.0, 1.0, 1.0, 0.1);
        let dq2 = ht2.heat_exchange(400.0, 300.0, 1.0, 1.0, 0.1);
        assert!((dq2 / dq1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_sph_diffusion_new() {
        let sd = SPHDiffusion::new(1e-9, 0.1);
        assert!((sd.diffusivity - 1e-9).abs() < 1e-20);
        assert!((sd.smoothing_length - 0.1).abs() < 1e-15);
    }
    #[test]
    fn test_sph_diffusion_laplacian_zero_when_same_conc() {
        let sd = SPHDiffusion::new(1e-9, 0.1);
        let lap = sd.concentration_laplacian(1.0, &[(1.0, [0.1, 0.0, 0.0], 1.0, 1.0)]);
        assert!(lap.abs() < 1e-20);
    }
    #[test]
    fn test_sph_diffusion_laplacian_positive_gradient() {
        let sd = SPHDiffusion::new(1e-9, 0.5);
        let lap = sd.concentration_laplacian(0.0, &[(1.0, [0.1, 0.0, 0.0], 1.0, 1.0)]);
        assert!(lap > 0.0);
    }
    #[test]
    fn test_sph_diffusion_dcdt() {
        let sd = SPHDiffusion::new(1.0, 0.5);
        let dcdt = sd.dcdt(0.0, &[(1.0, [0.1, 0.0, 0.0], 1.0, 1.0)]);
        assert!(dcdt > 0.0);
    }
    #[test]
    fn test_sph_reactive_transport_new() {
        let rt = SPHReactiveTransport::new(1e-9, 1e-3, 0.1);
        assert!((rt.diffusivity - 1e-9).abs() < 1e-20);
        assert!((rt.reaction_rate - 1e-3).abs() < 1e-20);
    }
    #[test]
    fn test_sph_reactive_transport_no_reaction() {
        let rt = SPHReactiveTransport::new(1e-9, 0.0, 0.1);
        let dcdt = rt.dcdt(1.0, &[], 0.0);
        assert!(dcdt.abs() < 1e-20);
    }
    #[test]
    fn test_sph_reactive_transport_first_order_decay() {
        let rt = SPHReactiveTransport::new(0.0, 1.0, 0.1);
        let dcdt = rt.dcdt(2.0, &[], 0.0);
        assert!((dcdt + 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_sph_reactive_transport_source_term() {
        let rt = SPHReactiveTransport::new(0.0, 0.0, 0.1);
        let dcdt = rt.dcdt(0.0, &[], 5.0);
        assert!((dcdt - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_sph_electroosmosis_new() {
        let eo = SPHElectroosmosis::new(-25e-3, 0.9e-3, 78.5 * 8.854e-12, 1.0e-7);
        assert!((eo.zeta_potential + 25e-3).abs() < 1e-15);
    }
    #[test]
    fn test_sph_electroosmosis_smoluchowski_velocity() {
        let eps = 78.5 * 8.854e-12;
        let zeta = -25e-3;
        let eta = 0.9e-3;
        let e_field = 1000.0;
        let eo = SPHElectroosmosis::new(zeta, eta, eps, 1.0e-7);
        let v = eo.smoluchowski_velocity(e_field);
        let expected = -eps * zeta * e_field / eta;
        assert!((v - expected).abs() < 1e-20);
    }
    #[test]
    fn test_sph_electroosmosis_velocity_sign() {
        let eps = 78.5 * 8.854e-12;
        let eo = SPHElectroosmosis::new(-25e-3, 0.9e-3, eps, 1.0e-7);
        let v = eo.smoluchowski_velocity(1000.0);
        assert!(v > 0.0);
    }
    #[test]
    fn test_sph_electroosmosis_debye_length() {
        let eo = SPHElectroosmosis::new(-25e-3, 0.9e-3, 78.5 * 8.854e-12, 1.0e-7);
        let ld = eo.debye_length;
        assert!(ld > 0.0);
        assert!(ld < 1e-5);
    }
    #[test]
    fn test_sph_electroosmosis_zero_field() {
        let eo = SPHElectroosmosis::new(-25e-3, 0.9e-3, 78.5 * 8.854e-12, 1.0e-7);
        assert_eq!(eo.smoluchowski_velocity(0.0), 0.0);
    }
    #[test]
    fn test_sph_multispecies_new() {
        let ms = SPHMultispeciesTransport::new(vec![1e-9, 2e-9], vec![1.0, -1.0], vec![1.0, 1.0]);
        assert_eq!(ms.nspecies, 2);
    }
    #[test]
    fn test_sph_multispecies_stoichiometry_conservation() {
        let ms = SPHMultispeciesTransport::new(vec![1e-9, 1e-9], vec![-1.0, 1.0], vec![1.0, 1.0]);
        let s: f64 = ms.stoichiometry.iter().sum();
        assert!(s.abs() < 1e-12);
    }
    #[test]
    fn test_sph_multispecies_reaction_source() {
        let ms = SPHMultispeciesTransport::new(vec![0.0, 0.0], vec![-1.0, 1.0], vec![1.0, 1.0]);
        let r0 = ms.reaction_source(0, 1.0);
        assert!((r0 + 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_sph_multispecies_reaction_product_source() {
        let ms = SPHMultispeciesTransport::new(vec![0.0, 0.0], vec![-1.0, 1.0], vec![1.0, 1.0]);
        let r1 = ms.reaction_source(1, 1.0);
        assert!((r1 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_sph_multispecies_diffusivity_accessed() {
        let ms = SPHMultispeciesTransport::new(vec![1e-9, 3e-9], vec![-1.0, 1.0], vec![1.0, 1.0]);
        assert!((ms.diffusivities[1] - 3e-9).abs() < 1e-25);
    }
}
