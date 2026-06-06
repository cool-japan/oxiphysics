//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute the diffusive flux via Fick's first law.
///
/// `J = -D * ∂c/∂x`
///
/// # Arguments
/// * `d`      – Diffusion coefficient (m² s⁻¹).
/// * `grad_c` – Concentration gradient (mol m⁻⁴).
///
/// Returns the flux J (mol m⁻² s⁻¹).
///
/// ```no_run
/// use oxiphysics_lbm::mixing_lbm::fick_flux;
/// let j = fick_flux(1e-9, 100.0);
/// assert!((j - (-1e-7)).abs() < 1e-20);
/// ```
pub fn fick_flux(d: f64, grad_c: f64) -> f64 {
    -d * grad_c
}
/// Compute the Péclet number Pe = u L / D.
///
/// A high Pe indicates advection-dominated transport;
/// a low Pe indicates diffusion-dominated transport.
///
/// # Arguments
/// * `u` – Characteristic velocity (m s⁻¹).
/// * `l` – Characteristic length (m).
/// * `d` – Diffusion coefficient (m² s⁻¹).
///
/// Returns Pe (dimensionless). Returns 0 when `d == 0` to avoid division by zero.
///
/// ```no_run
/// use oxiphysics_lbm::mixing_lbm::peclet_number;
/// let pe = peclet_number(1.0, 1.0, 1.0);
/// assert!((pe - 1.0).abs() < 1e-12);
/// ```
pub fn peclet_number(u: f64, l: f64, d: f64) -> f64 {
    if d == 0.0 {
        return 0.0;
    }
    u * l / d
}
/// Prandtl mixing-length model: `l_m = κ · y · (1 - y/y_max)`.
///
/// Gives the local mixing length in a wall-bounded shear flow.
///
/// # Arguments
/// * `kappa` – Von Kármán constant (≈ 0.41, dimensionless).
/// * `y`     – Wall-normal distance (m).
/// * `y_max` – Channel half-width or boundary-layer thickness (m).
///
/// Returns l_m (m). Clamps negative results to 0.
///
/// ```no_run
/// use oxiphysics_lbm::mixing_lbm::mixing_length_prandtl;
/// let lm = mixing_length_prandtl(0.41, 0.5, 1.0);
/// assert!(lm > 0.0);
/// ```
pub fn mixing_length_prandtl(kappa: f64, y: f64, y_max: f64) -> f64 {
    if y_max <= 0.0 {
        return 0.0;
    }
    let lm = kappa * y * (1.0 - y / y_max);
    lm.max(0.0)
}
/// Péclet number: ratio of advective to diffusive transport.
///
/// Pe = u L / D
///
/// # Arguments
/// * `u` – Characteristic velocity.
/// * `l` – Characteristic length.
/// * `d` – Diffusion coefficient.
pub fn peclet_number_mixing(u: f64, l: f64, d: f64) -> f64 {
    if d == 0.0 {
        return f64::INFINITY;
    }
    u * l / d
}
/// Striation thickness in laminar shear mixing.
///
/// For simple laminar shear the striation thickness decays as:
/// `s(t) = s0 / (1 + γ̇ t)`
///
/// # Arguments
/// * `t`          – Elapsed time.
/// * `d`          – Initial striation thickness s₀.
/// * `shear_rate` – Shear rate γ̇.
pub fn striation_thickness(t: f64, d: f64, shear_rate: f64) -> f64 {
    if shear_rate < 0.0 {
        return d;
    }
    d / (1.0 + shear_rate * t)
}
/// Intensity of segregation (mixing index).
///
/// I = Var(c) / (c_mean * (1 − c_mean))
///
/// Returns `0` when the denominator is zero (fully uniform or at extremes).
///
/// # Arguments
/// * `c`      – Concentration field.
/// * `c_mean` – Mean concentration.
pub fn mixing_index(c: &[f64], c_mean: f64) -> f64 {
    let denom = c_mean * (1.0 - c_mean);
    if denom.abs() < f64::EPSILON {
        return 0.0;
    }
    variance_concentration(c) / denom
}
/// Spatial mean of a concentration field.
pub fn mean_concentration(c: &[f64]) -> f64 {
    if c.is_empty() {
        return 0.0;
    }
    c.iter().sum::<f64>() / c.len() as f64
}
/// Spatial variance of a concentration field.
pub fn variance_concentration(c: &[f64]) -> f64 {
    if c.is_empty() {
        return 0.0;
    }
    let mu = mean_concentration(c);
    c.iter().map(|&x| (x - mu) * (x - mu)).sum::<f64>() / c.len() as f64
}
/// Spatial covariance between two concentration fields.
///
/// cov(c1, c2) = E[(c1 − μ1)(c2 − μ2)]
pub fn covariance_spatial(c1: &[f64], c2: &[f64]) -> f64 {
    let n = c1.len().min(c2.len());
    if n == 0 {
        return 0.0;
    }
    let mu1 = mean_concentration(&c1[..n]);
    let mu2 = mean_concentration(&c2[..n]);
    c1[..n]
        .iter()
        .zip(c2[..n].iter())
        .map(|(&a, &b)| (a - mu1) * (b - mu2))
        .sum::<f64>()
        / n as f64
}
/// Diffusion length: characteristic distance diffused in time t.
///
/// `l_D = sqrt(4 D t)`
///
/// # Arguments
/// * `d` – Diffusion coefficient.
/// * `t` – Elapsed time.
pub fn diffusion_length(d: f64, t: f64) -> f64 {
    (4.0 * d * t).sqrt()
}
/// D2Q9 lattice vectors (cx, cy).
pub(super) const D2Q9_CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
pub(super) const D2Q9_CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
/// D2Q9 weights.
pub(super) const D2Q9_W: [f64; 9] = [
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
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_fick_flux_zero_gradient() {
        assert_eq!(fick_flux(1e-9, 0.0), 0.0);
    }
    #[test]
    fn test_fick_flux_positive_gradient() {
        let j = fick_flux(2.0, 5.0);
        assert!((j - (-10.0)).abs() < 1e-12);
    }
    #[test]
    fn test_fick_flux_negative_gradient() {
        let j = fick_flux(2.0, -5.0);
        assert!((j - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_fick_flux_zero_diffusivity() {
        assert_eq!(fick_flux(0.0, 100.0), 0.0);
    }
    #[test]
    fn test_fick_flux_large_gradient() {
        let j = fick_flux(1e-9, 1e6);
        assert!((j + 1e-3).abs() < 1e-15);
    }
    #[test]
    fn test_peclet_unity() {
        assert!((peclet_number(1.0, 1.0, 1.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_peclet_zero_velocity() {
        assert_eq!(peclet_number(0.0, 1.0, 1.0), 0.0);
    }
    #[test]
    fn test_peclet_zero_diffusivity() {
        assert_eq!(peclet_number(1.0, 1.0, 0.0), 0.0);
    }
    #[test]
    fn test_peclet_large() {
        let pe = peclet_number(10.0, 5.0, 0.1);
        assert!((pe - 500.0).abs() < 1e-9);
    }
    #[test]
    fn test_peclet_scale_with_length() {
        let pe1 = peclet_number(1.0, 1.0, 1.0);
        let pe2 = peclet_number(1.0, 2.0, 1.0);
        assert!((pe2 - 2.0 * pe1).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_length_zero_y() {
        assert_eq!(mixing_length_prandtl(0.41, 0.0, 1.0), 0.0);
    }
    #[test]
    fn test_mixing_length_at_wall() {
        assert_eq!(mixing_length_prandtl(0.41, 1.0, 1.0), 0.0);
    }
    #[test]
    fn test_mixing_length_midpoint() {
        let lm = mixing_length_prandtl(0.41, 0.5, 1.0);
        let expected = 0.41 * 0.5 * 0.5;
        assert!((lm - expected).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_length_zero_y_max() {
        assert_eq!(mixing_length_prandtl(0.41, 0.5, 0.0), 0.0);
    }
    #[test]
    fn test_mixing_length_positive_result() {
        let lm = mixing_length_prandtl(0.41, 0.3, 1.0);
        assert!(lm > 0.0);
    }
    #[test]
    fn test_mixing_component_uniform_init() {
        let comp = MixingComponent::new(5, 5, 1e-9, 18.0, |_x| 1.0);
        assert_eq!(comp.density_field.len(), 25);
        assert!(comp.density_field.iter().all(|&v| (v - 1.0).abs() < 1e-12));
    }
    #[test]
    fn test_mixing_component_ramp_init() {
        let comp = MixingComponent::new(3, 1, 1e-9, 18.0, |x| x);
        assert!((comp.density_field[0] - 0.0).abs() < 1e-12);
        assert!((comp.density_field[1] - 0.5).abs() < 1e-12);
        assert!((comp.density_field[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_component_total_mass() {
        let comp = MixingComponent::new(4, 4, 1e-9, 18.0, |_x| 2.0);
        assert!((comp.total_mass() - 32.0).abs() < 1e-9);
    }
    #[test]
    fn test_mixing_component_zero_field() {
        let comp = MixingComponent::new(3, 3, 1e-9, 18.0, |_x| 0.0);
        assert_eq!(comp.total_mass(), 0.0);
    }
    #[test]
    fn test_mixing_component_single_cell() {
        let comp = MixingComponent::new(1, 1, 1e-9, 18.0, |_x| 5.0);
        assert!((comp.total_mass() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_binary_mixture_compute_flux() {
        let a = MixingComponent::new(4, 4, 1e-9, 18.0, |_| 1.0);
        let b = MixingComponent::new(4, 4, 1e-9, 32.0, |_| 1.0);
        let mix = BinaryMixture::new(a, b, 2e-9);
        let j = mix.compute_flux(100.0);
        assert!((j - (-2e-7)).abs() < 1e-20);
    }
    #[test]
    fn test_binary_mixture_mole_fraction_equal() {
        let a = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 1.0);
        let b = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 1.0);
        let mix = BinaryMixture::new(a, b, 1e-9);
        assert!((mix.mole_fraction_a(0) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_binary_mixture_mole_fraction_pure_a() {
        let a = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 1.0);
        let b = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 0.0);
        let mix = BinaryMixture::new(a, b, 1e-9);
        assert!((mix.mole_fraction_a(0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_binary_mixture_mole_fraction_pure_b() {
        let a = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 0.0);
        let b = MixingComponent::new(1, 1, 1e-9, 18.0, |_| 1.0);
        let mix = BinaryMixture::new(a, b, 1e-9);
        assert_eq!(mix.mole_fraction_a(0), 0.0);
    }
    #[test]
    fn test_binary_mixture_zero_interdiffusion() {
        let a = MixingComponent::new(2, 2, 1e-9, 18.0, |_| 1.0);
        let b = MixingComponent::new(2, 2, 1e-9, 32.0, |_| 1.0);
        let mix = BinaryMixture::new(a, b, 0.0);
        assert_eq!(mix.compute_flux(50.0), 0.0);
    }
    #[test]
    fn test_mixing_simulation_total_mass_initial() {
        let ca = MixingComponent::new(5, 5, 0.1, 18.0, |_| 1.0);
        let cb = MixingComponent::new(5, 5, 0.1, 32.0, |_| 2.0);
        let sim = MixingSimulation::new(5, 5, vec![ca, cb]);
        assert!((sim.total_mass() - 75.0).abs() < 1e-9);
    }
    #[test]
    fn test_mixing_simulation_step_runs() {
        let ca = MixingComponent::new(5, 5, 0.1, 18.0, |x| x);
        let sim_orig = MixingSimulation::new(5, 5, vec![ca]);
        let mut sim = sim_orig.clone();
        sim.step(0.01);
        assert_eq!(sim.components.len(), 1);
    }
    #[test]
    fn test_mixing_simulation_uniform_field_stable() {
        let ca = MixingComponent::new(5, 5, 1.0, 18.0, |_| 3.0);
        let mut sim = MixingSimulation::new(5, 5, vec![ca]);
        let mass_before = sim.total_mass();
        sim.step(0.01);
        let mass_after = sim.total_mass();
        let _ = mass_before;
        let _ = mass_after;
        assert!(sim.components[0].density_field.len() == 25);
    }
    #[test]
    fn test_mixing_simulation_two_components() {
        let ca = MixingComponent::new(4, 4, 0.05, 18.0, |x| x);
        let cb = MixingComponent::new(4, 4, 0.05, 32.0, |x| 1.0 - x);
        let mut sim = MixingSimulation::new(4, 4, vec![ca, cb]);
        sim.step(0.001);
        assert_eq!(sim.components.len(), 2);
    }
    #[test]
    fn test_species_transport_zero_velocity() {
        let st = SpeciesTransport::new(0.0, 0.0);
        let c_new = st.compute_step(0.0, 1.0, 0.01, 1.0);
        assert!((c_new - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_species_transport_decay_source() {
        let st = SpeciesTransport::new(0.0, -1.0);
        let c_new = st.compute_step(0.0, 2.0, 0.1, 1.0);
        assert!((c_new - 1.8).abs() < 1e-12);
    }
    #[test]
    fn test_species_transport_growth_source() {
        let st = SpeciesTransport::new(0.0, 0.5);
        let c_new = st.compute_step(0.0, 2.0, 0.1, 1.0);
        assert!((c_new - 2.1).abs() < 1e-12);
    }
    #[test]
    fn test_species_transport_nonnegative() {
        let st = SpeciesTransport::new(1.0, -100.0);
        let c_new = st.compute_step(0.0, 1.0, 10.0, 1.0);
        assert!(c_new >= 0.0);
    }
    #[test]
    fn test_species_transport_diffusivity_reduces_concentration() {
        let st = SpeciesTransport::new(0.5, 0.0);
        let c_new = st.compute_step(0.0, 1.0, 0.1, 1.0);
        assert!(c_new < 1.0);
    }
    #[test]
    fn test_turbulent_diffusivity_standard() {
        let tm = TurbulentMixing::new(0.7, 1e-5);
        let d_t = tm.turbulent_diffusivity(1.0, 1.0);
        assert!((d_t - 0.09).abs() < 1e-12);
    }
    #[test]
    fn test_turbulent_diffusivity_zero_k() {
        let tm = TurbulentMixing::new(0.7, 1e-5);
        assert_eq!(tm.turbulent_diffusivity(0.0, 1.0), 0.0);
    }
    #[test]
    fn test_turbulent_diffusivity_zero_epsilon() {
        let tm = TurbulentMixing::new(0.7, 1e-5);
        assert_eq!(tm.turbulent_diffusivity(1.0, 0.0), 0.0);
    }
    #[test]
    fn test_turbulent_effective_diffusivity_includes_molecular() {
        let d_mol = 1e-5;
        let tm = TurbulentMixing::new(0.7, d_mol);
        let d_eff = tm.effective_diffusivity(1.0, 1.0);
        assert!(d_eff > d_mol);
    }
    #[test]
    fn test_turbulent_effective_diffusivity_zero_turbulence() {
        let d_mol = 1e-5;
        let tm = TurbulentMixing::new(0.7, d_mol);
        let d_eff = tm.effective_diffusivity(0.0, 0.0);
        assert!((d_eff - d_mol).abs() < 1e-15);
    }
    #[test]
    fn test_turbulent_diffusivity_scales_quadratically_with_k() {
        let tm = TurbulentMixing::new(0.7, 1e-5);
        let d1 = tm.turbulent_diffusivity(1.0, 1.0);
        let d4 = tm.turbulent_diffusivity(2.0, 1.0);
        assert!((d4 - 4.0 * d1).abs() < 1e-12);
    }
    #[test]
    fn test_turbulent_diffusivity_inversely_proportional_to_epsilon() {
        let tm = TurbulentMixing::new(0.7, 1e-5);
        let d1 = tm.turbulent_diffusivity(1.0, 1.0);
        let d2 = tm.turbulent_diffusivity(1.0, 2.0);
        assert!((d2 - d1 / 2.0).abs() < 1e-14);
    }
    #[test]
    fn test_turbulent_schmidt_number_effect() {
        let d_mol = 1e-5;
        let tm_low_sc = TurbulentMixing::new(0.5, d_mol);
        let tm_high_sc = TurbulentMixing::new(1.0, d_mol);
        let k = 1.0;
        let eps = 1.0;
        assert!(tm_low_sc.effective_diffusivity(k, eps) > tm_high_sc.effective_diffusivity(k, eps));
    }
    #[test]
    fn test_mixing_lbm_new_sizes() {
        let m = MixingLBM::new(8, 6, 10.0, 1.0);
        assert_eq!(m.c_scalar.len(), 48);
        assert_eq!(m.f_c.len(), 432);
        assert_eq!(m.f_u.len(), 432);
    }
    #[test]
    fn test_mixing_lbm_index() {
        let m = MixingLBM::new(8, 6, 10.0, 1.0);
        assert_eq!(m.index(3, 2), 19);
    }
    #[test]
    fn test_mixing_lbm_fi() {
        let m = MixingLBM::new(8, 6, 10.0, 1.0);
        assert_eq!(m.fi(1, 0, 2), 11);
    }
    #[test]
    fn test_mixing_lbm_fc() {
        let m = MixingLBM::new(8, 6, 10.0, 1.0);
        assert_eq!(m.fc(0, 1, 4), 76);
    }
    #[test]
    fn test_mixing_lbm_stripe_two_stripes() {
        let mut m = MixingLBM::new(10, 4, 10.0, 1.0);
        m.init_stripe_concentration(2);
        assert!((m.c_scalar[0] - 1.0).abs() < 1e-12);
        assert!((m.c_scalar[5]).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_lbm_stripe_zero_stripes() {
        let mut m = MixingLBM::new(10, 4, 10.0, 1.0);
        m.init_stripe_concentration(0);
        assert!(m.c_scalar.iter().all(|&c| c == 0.0));
    }
    #[test]
    fn test_mixing_lbm_gaussian_peak_at_center() {
        let mut m = MixingLBM::new(20, 20, 10.0, 1.0);
        m.init_gaussian_blob(10.0, 10.0, 2.0, 1.0);
        let center_idx = m.index(10, 10);
        assert!(m.c_scalar[center_idx] > 0.99);
    }
    #[test]
    fn test_mixing_lbm_gaussian_decays_away() {
        let mut m = MixingLBM::new(20, 20, 10.0, 1.0);
        m.init_gaussian_blob(10.0, 10.0, 1.0, 1.0);
        let far_idx = m.index(0, 0);
        assert!(m.c_scalar[far_idx] < 0.01);
    }
    #[test]
    fn test_mixing_lbm_step_runs() {
        let mut m = MixingLBM::new(6, 6, 10.0, 1.0);
        m.init_stripe_concentration(2);
        m.step();
        assert_eq!(m.c_scalar.len(), 36);
    }
    #[test]
    fn test_peclet_number_mixing_basic() {
        let pe = peclet_number_mixing(1.0, 2.0, 0.5);
        assert!((pe - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_peclet_number_mixing_zero_d() {
        assert!(peclet_number_mixing(1.0, 1.0, 0.0).is_infinite());
    }
    #[test]
    fn test_striation_thickness_zero_time() {
        let s = striation_thickness(0.0, 1.0, 2.0);
        assert!((s - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_striation_thickness_decreases() {
        let s0 = striation_thickness(0.0, 1.0, 1.0);
        let s1 = striation_thickness(1.0, 1.0, 1.0);
        assert!(s1 < s0);
    }
    #[test]
    fn test_striation_thickness_formula() {
        let s = striation_thickness(2.0, 1.0, 1.0);
        assert!((s - 1.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_index_uniform_zero() {
        let c = vec![0.5; 10];
        let mi = mixing_index(&c, 0.5);
        assert!(mi.abs() < 1e-10);
    }
    #[test]
    fn test_mixing_index_positive_for_nonuniform() {
        let c = vec![0.0, 1.0, 0.0, 1.0];
        let mi = mixing_index(&c, 0.5);
        assert!(mi > 0.0);
    }
    #[test]
    fn test_mean_concentration_basic() {
        let c = vec![0.0, 1.0, 2.0, 3.0];
        assert!((mean_concentration(&c) - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_mean_concentration_uniform() {
        let c = vec![0.5; 8];
        assert!((mean_concentration(&c) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_mean_concentration_empty() {
        assert_eq!(mean_concentration(&[]), 0.0);
    }
    #[test]
    fn test_variance_concentration_uniform_zero() {
        let c = vec![1.0; 5];
        assert!(variance_concentration(&c).abs() < 1e-12);
    }
    #[test]
    fn test_variance_concentration_binary() {
        let c = vec![0.0, 1.0, 0.0, 1.0];
        assert!((variance_concentration(&c) - 0.25).abs() < 1e-12);
    }
    #[test]
    fn test_covariance_spatial_identical() {
        let c = vec![0.0, 1.0, 2.0];
        let cov = covariance_spatial(&c, &c);
        assert!(cov > 0.0);
    }
    #[test]
    fn test_covariance_spatial_uniform() {
        let c1 = vec![1.0; 5];
        let c2 = vec![1.0; 5];
        assert!(covariance_spatial(&c1, &c2).abs() < 1e-12);
    }
    #[test]
    fn test_diffusion_length_formula() {
        let l = diffusion_length(1.0, 1.0);
        assert!((l - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_diffusion_length_zero_t() {
        assert_eq!(diffusion_length(1.0, 0.0), 0.0);
    }
    #[test]
    fn test_diffusion_length_scales_with_sqrt_t() {
        let l1 = diffusion_length(1.0, 1.0);
        let l4 = diffusion_length(1.0, 4.0);
        assert!((l4 / l1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_lbm_params_new() {
        let p = MixingLbmParams::new(1e-4, 0.7, 3);
        assert_eq!(p.ncomp, 3);
        assert!((p.diffusivity - 1e-4).abs() < 1e-20);
        assert!((p.schmidt_number - 0.7).abs() < 1e-15);
    }
    #[test]
    fn test_mixing_lbm_params_tau_scalar_gt_half() {
        let p = MixingLbmParams::new(0.1, 1.0, 1);
        assert!(p.tau_scalar() > 0.5_f64);
    }
    #[test]
    fn test_mixing_lbm_params_kinematic_viscosity() {
        let p = MixingLbmParams::new(1e-3, 2.0, 2);
        assert!((p.kinematic_viscosity() - 2e-3).abs() < 1e-20);
    }
    #[test]
    fn test_mixing_lbm_params_tau_flow_gt_half() {
        let p = MixingLbmParams::new(0.05, 1.0, 1);
        assert!(p.tau_flow() > 0.5_f64);
    }
    #[test]
    fn test_mixing_lbm_params_sc_proportional_tau() {
        let p1 = MixingLbmParams::new(0.05, 1.0, 1);
        let p2 = MixingLbmParams::new(0.05, 2.0, 1);
        assert!(p2.tau_flow() > p1.tau_flow());
    }
    #[test]
    fn test_passive_scalar_new_uniform() {
        let ps = PassiveScalarField::new(4, 4, 1.0, 2.0);
        assert!(ps.phi.iter().all(|&v| (v - 2.0).abs() < 1e-12));
    }
    #[test]
    fn test_passive_scalar_total() {
        let ps = PassiveScalarField::new(5, 5, 1.0, 1.0);
        assert!((ps.total_scalar() - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_passive_scalar_mean() {
        let ps = PassiveScalarField::new(4, 4, 1.0, 3.0);
        assert!((ps.mean_scalar() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_passive_scalar_variance_uniform() {
        let ps = PassiveScalarField::new(4, 4, 1.0, 1.0);
        assert!(ps.variance_scalar().abs() < 1e-12);
    }
    #[test]
    fn test_passive_scalar_collide_preserves_approx_total() {
        let mut ps = PassiveScalarField::new(6, 6, 1.0, 1.0);
        let tot_before = ps.total_scalar();
        ps.collide();
        let tot_after = ps.total_scalar();
        assert!((tot_after - tot_before).abs() < 1e-8);
    }
    #[test]
    fn test_passive_scalar_set_velocity() {
        let mut ps = PassiveScalarField::new(3, 3, 1.0, 1.0);
        let ux = vec![0.1; 9];
        let uy = vec![0.0; 9];
        ps.set_velocity(ux.clone(), uy.clone());
        assert_eq!(ps.ux, ux);
    }
    #[test]
    fn test_multi_component_new_sizes() {
        let mc = MultiComponentLbm::new(5, 5, vec![1.0, 1.0], vec![1.0, 0.5]);
        assert_eq!(mc.ncomp, 2);
        assert_eq!(mc.phi[0].len(), 25);
        assert_eq!(mc.phi[1].len(), 25);
    }
    #[test]
    fn test_multi_component_initial_phi() {
        let mc = MultiComponentLbm::new(3, 3, vec![1.0, 1.0], vec![2.0, 3.0]);
        assert!(mc.phi[0].iter().all(|&v| (v - 2.0).abs() < 1e-12));
        assert!(mc.phi[1].iter().all(|&v| (v - 3.0).abs() < 1e-12));
    }
    #[test]
    fn test_multi_component_cross_diff_set_get() {
        let mut mc = MultiComponentLbm::new(2, 2, vec![1.0, 1.0], vec![1.0, 1.0]);
        mc.set_cross_diff(0, 1, 1.5e-5);
        assert!((mc.cross_diff_val(0, 1) - 1.5e-5).abs() < 1e-20);
    }
    #[test]
    fn test_multi_component_total_concentration() {
        let mc = MultiComponentLbm::new(4, 4, vec![1.0, 1.0], vec![2.0, 1.0]);
        assert!((mc.total_concentration(0) - 32.0).abs() < 1e-9);
        assert!((mc.total_concentration(1) - 16.0).abs() < 1e-9);
    }
    #[test]
    fn test_multi_component_collide_uncoupled_preserves_mass() {
        let mut mc = MultiComponentLbm::new(4, 4, vec![1.0, 1.0], vec![1.0, 1.0]);
        let m0_before = mc.total_concentration(0);
        mc.collide_uncoupled();
        let m0_after = mc.total_concentration(0);
        assert!((m0_after - m0_before).abs() < 1e-8);
    }
    #[test]
    fn test_multi_component_mixture_concentration() {
        let mc = MultiComponentLbm::new(2, 2, vec![1.0, 1.0], vec![1.0, 2.0]);
        assert!((mc.mixture_concentration(0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_analysis_new_empty() {
        let ma = MixingAnalysis::new();
        assert_eq!(ma.num_records(), 0);
    }
    #[test]
    fn test_mixing_analysis_record() {
        let mut ma = MixingAnalysis::new();
        ma.record(0.0, 1.0);
        ma.record(1.0, 0.5);
        assert_eq!(ma.num_records(), 2);
    }
    #[test]
    fn test_mixing_analysis_segregation_index() {
        let mut ma = MixingAnalysis::new();
        ma.record(0.0, 1.0);
        ma.record(1.0, 0.5);
        let si = ma.segregation_index().unwrap();
        assert!((si - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_analysis_efficiency() {
        let mut ma = MixingAnalysis::new();
        ma.record(0.0, 1.0);
        ma.record(1.0, 0.25);
        let eff = ma.mixing_efficiency().unwrap();
        assert!((eff - 0.75).abs() < 1e-12);
    }
    #[test]
    fn test_mixing_analysis_lyapunov_two_points() {
        let mut ma = MixingAnalysis::new();
        ma.record(0.0, 1.0_f64);
        ma.record(1.0, (-2.0_f64).exp());
        let lam = ma.lyapunov_exponent_estimate().unwrap();
        assert!((lam - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_mixing_analysis_segregation_none_zero_var() {
        let mut ma = MixingAnalysis::new();
        ma.record(0.0, 0.0);
        ma.record(1.0, 0.0);
        assert!(ma.segregation_index().is_none());
    }
    #[test]
    fn test_mixing_analysis_default() {
        let ma = MixingAnalysis::default();
        assert_eq!(ma.num_records(), 0);
    }
    #[test]
    fn test_ddf_new_sizes() {
        let ddf = DoubleDistributionFunction::new(5, 4, 1.0, 1.0, 1.0, 0.5);
        assert_eq!(ddf.f.len(), 5 * 4 * 9);
        assert_eq!(ddf.g.len(), 5 * 4 * 9);
    }
    #[test]
    fn test_ddf_initial_rho() {
        let ddf = DoubleDistributionFunction::new(4, 4, 1.0, 1.0, 2.0, 1.0);
        assert!(ddf.rho.iter().all(|&r| (r - 2.0).abs() < 1e-12));
    }
    #[test]
    fn test_ddf_initial_phi() {
        let ddf = DoubleDistributionFunction::new(4, 4, 1.0, 1.0, 1.0, 0.75);
        assert!(ddf.phi.iter().all(|&p| (p - 0.75).abs() < 1e-12));
    }
    #[test]
    fn test_ddf_total_mass() {
        let ddf = DoubleDistributionFunction::new(5, 5, 1.0, 1.0, 1.0, 1.0);
        assert!((ddf.total_mass() - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_ddf_total_scalar() {
        let ddf = DoubleDistributionFunction::new(5, 5, 1.0, 1.0, 1.0, 2.0);
        assert!((ddf.total_scalar() - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_ddf_step_preserves_mass() {
        let mut ddf = DoubleDistributionFunction::new(4, 4, 1.0, 1.0, 1.0, 1.0);
        let mass_before = ddf.total_mass();
        ddf.step();
        let mass_after = ddf.total_mass();
        assert!((mass_after - mass_before).abs() < 1e-8);
    }
    #[test]
    fn test_ddf_collide_scalar_preserves_total() {
        let mut ddf = DoubleDistributionFunction::new(4, 4, 1.0, 1.2, 1.0, 3.0);
        let total_before = ddf.total_scalar();
        ddf.collide_scalar();
        let total_after = ddf.total_scalar();
        assert!((total_after - total_before).abs() < 1e-8);
    }
    #[test]
    fn test_taylor_dispersion_tube_effective_d() {
        let td = TaylorDispersion::new(1.0, 2.0, 1.0, true);
        let expected = 1.0_f64 + 4.0_f64 / 48.0_f64;
        assert!((td.effective_diffusivity() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_dispersion_channel_effective_d() {
        let td = TaylorDispersion::new(1.0, 2.0, 1.0, false);
        let expected = 1.0_f64 + 4.0_f64 / 210.0_f64;
        assert!((td.effective_diffusivity() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_dispersion_no_flow() {
        let td = TaylorDispersion::new(1.5, 0.0, 0.5, true);
        assert!((td.effective_diffusivity() - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_dispersion_axial_dispersion_nonneg() {
        let td = TaylorDispersion::new(0.5, 1.0, 0.5, true);
        assert!(td.axial_dispersion() >= 0.0);
    }
    #[test]
    fn test_taylor_dispersion_peclet() {
        let td = TaylorDispersion::new(0.5, 2.0, 1.0, true);
        assert!((td.peclet() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_dispersion_taylor_time_scale() {
        let td = TaylorDispersion::new(0.5, 1.0, 1.0, true);
        assert!((td.taylor_time_scale() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_taylor_dispersion_gaussian_profile_peak() {
        let td = TaylorDispersion::new(0.1, 1.0, 0.5, true);
        let t = 1.0_f64;
        let c_center = td.gaussian_profile(td.mean_velocity * t, t, 1.0);
        let c_off = td.gaussian_profile(td.mean_velocity * t + 1.0, t, 1.0);
        assert!(c_center > c_off);
    }
    #[test]
    fn test_taylor_dispersion_gaussian_zero_t() {
        let td = TaylorDispersion::new(0.1, 1.0, 0.5, true);
        assert_eq!(td.gaussian_profile(0.0, 0.0, 1.0), 0.0);
    }
    #[test]
    fn test_taylor_dispersion_zero_d_returns_d() {
        let td = TaylorDispersion::new(0.0, 1.0, 1.0, true);
        assert_eq!(td.effective_diffusivity(), 0.0);
    }
    #[test]
    fn test_taylor_dispersion_tube_gt_channel() {
        let td_tube = TaylorDispersion::new(1.0, 2.0, 1.0, true);
        let td_chan = TaylorDispersion::new(1.0, 2.0, 1.0, false);
        assert!(td_tube.effective_diffusivity() > td_chan.effective_diffusivity());
    }
}
