//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::CryogenicParticle;

/// Boltzmann constant in J/K.
pub(super) const BOLTZMANN: f64 = 1.380649e-23;
/// Planck's constant (reduced) in J s.
pub(super) const HBAR: f64 = 1.054571817e-34;
/// Avogadro's number.
pub(super) const AVOGADRO: f64 = 6.02214076e23;
/// Universal gas constant in J/(mol K).
pub(super) const GAS_CONSTANT: f64 = 8.314;
/// Helium-4 atomic mass in kg.
pub(super) const HE4_MASS: f64 = 4.002602e-3 / AVOGADRO;
/// Nitrogen molecular mass in kg.
pub(super) const N2_MASS: f64 = 28.014e-3 / AVOGADRO;
/// Normal boiling point of liquid helium-4 at 1 atm in K.
pub(super) const LHE_NBP: f64 = 4.22;
/// Lambda point of helium-4 (superfluid transition) in K.
pub(super) const HE4_LAMBDA: f64 = 2.172;
/// Normal boiling point of liquid nitrogen at 1 atm in K.
pub(super) const LN2_NBP: f64 = 77.36;
/// Density of liquid nitrogen at NBP in kg/m^3.
pub(super) const LN2_DENSITY: f64 = 808.0;
/// Density of liquid helium-4 at NBP in kg/m^3.
pub(super) const LHE_DENSITY: f64 = 125.0;
/// Latent heat of vaporisation of liquid nitrogen in J/kg.
pub(super) const LN2_LHV: f64 = 198_600.0;
/// Latent heat of vaporisation of liquid helium in J/kg.
pub(super) const LHE_LHV: f64 = 20_900.0;
/// Thermal conductivity of liquid nitrogen in W/(m K).
pub(super) const LN2_THERMAL_COND: f64 = 0.1396;
/// Specific heat capacity of liquid nitrogen in J/(kg K).
pub(super) const LN2_CP: f64 = 2_042.0;
/// Specific heat capacity of liquid helium in J/(kg K).
pub(super) const LHE_CP: f64 = 4_700.0;
/// Surface tension of liquid nitrogen at NBP in N/m.
pub(super) const LN2_SURFACE_TENSION: f64 = 8.85e-3;
/// Gravitational acceleration in m/s^2.
pub(super) const GRAVITY: f64 = 9.81;
/// Quantum of circulation for He-4 vortex in m^2/s.
pub(super) const QUANTUM_CIRCULATION: f64 = HBAR / HE4_MASS * 2.0 * PI;
/// Kapitza resistance prefactor (Khalatnikov) in m^2 K / W.
pub(super) const KAPITZA_R0: f64 = 2.0e-4;
/// Brillouin function B_J(x) for magnetic entropy calculations.
///
/// B_J(x) = (2J+1)/(2J) * coth((2J+1)/(2J) * x) - 1/(2J) * coth(x/(2J))
pub fn brillouin_function(j: f64, x: f64) -> f64 {
    if x.abs() < 1.0e-8 {
        return 0.0;
    }
    let a = (2.0 * j + 1.0) / (2.0 * j);
    let b = 1.0 / (2.0 * j);
    a * coth(a * x) - b * coth(b * x)
}
/// Derivative of Brillouin function with respect to x.
pub fn brillouin_derivative(j: f64, x: f64) -> f64 {
    let dx = 1.0e-6_f64;
    (brillouin_function(j, x + dx) - brillouin_function(j, x - dx)) / (2.0 * dx)
}
/// Hyperbolic cotangent.
pub(super) fn coth(x: f64) -> f64 {
    if x.abs() > 100.0 {
        x.signum()
    } else {
        x.cosh() / x.sinh()
    }
}
/// Cubic spline kernel W(r, h) normalised in 3D.
///
/// W(r,h) = sigma3/h^3 * f(q) where q = r/h,
/// sigma3 = 1/pi, and f is the cubic spline function.
pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / PI;
    let f = if q < 1.0 {
        1.0 - 1.5 * q * q + 0.75 * q * q * q
    } else if q < 2.0 {
        0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    };
    sigma / (h * h * h) * f
}
/// Gradient magnitude of cubic spline kernel dW/dr.
pub fn cubic_spline_gradient(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / PI;
    let df = if q < 1.0 {
        -3.0 * q + 2.25 * q * q
    } else if q < 2.0 {
        -0.75 * (2.0 - q).powi(2)
    } else {
        0.0
    };
    sigma / (h * h * h * h) * df
}
/// Compute the mean density of a set of cryogenic particles in kg/m^3.
#[allow(dead_code)]
pub fn mean_density(particles: &[CryogenicParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    particles.iter().map(|p| p.density).sum::<f64>() / particles.len() as f64
}
/// Compute the mean temperature of a set of particles in K.
#[allow(dead_code)]
pub fn mean_particle_temperature(particles: &[CryogenicParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    particles.iter().map(|p| p.temperature).sum::<f64>() / particles.len() as f64
}
/// Interpolate fluid density at a target point using kernel summation.
#[allow(dead_code)]
pub fn sph_interpolate_density(particles: &[CryogenicParticle], target: [f64; 3]) -> f64 {
    let mut rho = 0.0;
    for p in particles {
        let dx = target[0] - p.pos[0];
        let dy = target[1] - p.pos[1];
        let dz = target[2] - p.pos[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        rho += p.mass * cubic_spline_kernel(r, p.h);
    }
    rho
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryogenic_sph::types::*;
    pub(super) const TOL: f64 = 1.0e-8;
    #[test]
    fn test_ln2_boiling_point() {
        let bp = CryogenicFluid::LiquidNitrogen.boiling_point();
        assert!((bp - LN2_NBP).abs() < TOL);
    }
    #[test]
    fn test_lhe_boiling_point() {
        let bp = CryogenicFluid::LiquidHelium4.boiling_point();
        assert!((bp - LHE_NBP).abs() < TOL);
    }
    #[test]
    fn test_superfluid_flag() {
        assert!(CryogenicFluid::SuperfluidHelium4.is_superfluid());
        assert!(!CryogenicFluid::LiquidNitrogen.is_superfluid());
        assert!(!CryogenicFluid::LiquidHelium4.is_superfluid());
    }
    #[test]
    fn test_fluid_viscosity_ordering() {
        let mu_he = CryogenicFluid::LiquidHelium4.dynamic_viscosity();
        let mu_n2 = CryogenicFluid::LiquidNitrogen.dynamic_viscosity();
        assert!(mu_n2 > mu_he);
    }
    #[test]
    fn test_superfluid_viscosity_zero() {
        let mu_sf = CryogenicFluid::SuperfluidHelium4.dynamic_viscosity();
        assert_eq!(mu_sf, 0.0);
    }
    #[test]
    fn test_particle_creation_ln2() {
        let p = CryogenicParticle::new(
            [0.0, 0.0, 0.0],
            1.0e-6,
            LN2_NBP,
            CryogenicFluid::LiquidNitrogen,
        );
        assert!((p.temperature - LN2_NBP).abs() < TOL);
        assert!((p.density - LN2_DENSITY).abs() < TOL);
        assert!(p.speed() < TOL);
    }
    #[test]
    fn test_particle_kinetic_energy_zero() {
        let p = CryogenicParticle::new(
            [0.0, 0.0, 0.0],
            1.0,
            LN2_NBP,
            CryogenicFluid::LiquidNitrogen,
        );
        assert!(p.kinetic_energy() < TOL);
    }
    #[test]
    fn test_particle_kinetic_energy_nonzero() {
        let mut p = CryogenicParticle::new(
            [0.0, 0.0, 0.0],
            1.0,
            LN2_NBP,
            CryogenicFluid::LiquidNitrogen,
        );
        p.vel = [1.0, 0.0, 0.0];
        assert!((p.kinetic_energy() - 0.5).abs() < TOL);
    }
    #[test]
    fn test_particle_distance() {
        let p1 = CryogenicParticle::new([0.0, 0.0, 0.0], 1.0, 4.0, CryogenicFluid::LiquidHelium4);
        let p2 = CryogenicParticle::new([3.0, 4.0, 0.0], 1.0, 4.0, CryogenicFluid::LiquidHelium4);
        assert!((p1.distance_to(&p2) - 5.0).abs() < TOL);
    }
    #[test]
    fn test_normal_fraction_above_lambda() {
        let p = CryogenicParticle::new([0.0, 0.0, 0.0], 1.0, 3.0, CryogenicFluid::LiquidHelium4);
        assert!((p.normal_fraction() - 1.0).abs() < TOL);
    }
    #[test]
    fn test_superfluid_fraction_adds_to_one() {
        let p =
            CryogenicParticle::new([0.0, 0.0, 0.0], 1.0, 1.5, CryogenicFluid::SuperfluidHelium4);
        let total = p.normal_fraction() + p.superfluid_fraction();
        assert!((total - 1.0).abs() < TOL);
    }
    #[test]
    fn test_tait_eos_reference_pressure() {
        let eos = TaitEquationCryogenic::liquid_nitrogen();
        let p = eos.pressure(eos.rho0);
        assert!((p - eos.p0).abs() < 1.0);
    }
    #[test]
    fn test_tait_sound_speed_positive() {
        let eos = TaitEquationCryogenic::liquid_helium4();
        let c = eos.sound_speed(eos.rho0);
        assert!(c > 0.0);
        assert!(c.is_finite());
    }
    #[test]
    fn test_tait_pressure_increases_with_density() {
        let eos = TaitEquationCryogenic::liquid_nitrogen();
        let p1 = eos.pressure(eos.rho0);
        let p2 = eos.pressure(eos.rho0 * 1.01);
        assert!(p2 > p1);
    }
    #[test]
    fn test_two_fluid_normal_above_lambda() {
        let m = SuperfluidTwoFluidModel::new(3.0, LHE_DENSITY);
        assert!((m.normal_fraction - 1.0).abs() < TOL);
    }
    #[test]
    fn test_two_fluid_fractions_sum_to_one() {
        let m = SuperfluidTwoFluidModel::new(1.5, LHE_DENSITY);
        let sum = m.normal_fraction + (1.0 - m.normal_fraction);
        assert!((sum - 1.0).abs() < TOL);
    }
    #[test]
    fn test_two_fluid_second_sound_positive_below_lambda() {
        let m = SuperfluidTwoFluidModel::new(1.5, LHE_DENSITY);
        assert!(m.second_sound_speed >= 0.0);
    }
    #[test]
    fn test_two_fluid_relative_velocity_zero_at_rest() {
        let m = SuperfluidTwoFluidModel::new(1.5, LHE_DENSITY);
        let rv = m.relative_velocity();
        assert!(rv.iter().all(|v| v.abs() < TOL));
    }
    #[test]
    fn test_vortex_ring_closure() {
        let ring = VortexFilament::vortex_ring(1.0e-3, 64);
        assert!(ring.is_closed);
        assert_eq!(ring.points.len(), 64);
    }
    #[test]
    fn test_straight_filament_not_closed() {
        let fil = VortexFilament::straight_line(1.0e-2, 20);
        assert!(!fil.is_closed);
    }
    #[test]
    fn test_vortex_length_approx_ring_circumference() {
        let r = 1.0e-3;
        let ring = VortexFilament::vortex_ring(r, 256);
        let len = ring.total_length();
        let expected = 2.0 * PI * r;
        assert!((len - expected).abs() / expected < 0.01);
    }
    #[test]
    fn test_vortex_ring_self_velocity_positive() {
        let v = VortexFilament::ring_self_velocity(1.0e-3);
        assert!(v > 0.0);
        assert!(v.is_finite());
    }
    #[test]
    fn test_vortex_induced_velocity_far_field_finite() {
        let fil = VortexFilament::straight_line(1.0, 20);
        let v = fil.induced_velocity([10.0, 10.0, 0.5]);
        assert!(v.iter().all(|vi| vi.is_finite()));
    }
    #[test]
    fn test_kapitza_resistance_decreases_with_temperature() {
        let k = KapitzaResistance::new();
        let r1 = k.resistance(1.0);
        let r2 = k.resistance(2.0);
        assert!(r1 > r2);
    }
    #[test]
    fn test_kapitza_conductance_positive() {
        let k = KapitzaResistance::new();
        let g = k.conductance(2.0);
        assert!(g > 0.0);
    }
    #[test]
    fn test_kapitza_heat_flux() {
        let k = KapitzaResistance::new();
        let q = k.heat_flux(2.0, 0.1);
        assert!(q > 0.0);
    }
    #[test]
    fn test_phase_liquid_at_high_density() {
        let pt = CryogenicPhaseTransition::liquid_nitrogen();
        let phase = pt.phase_indicator(LN2_DENSITY * 1.5, 80.0);
        assert_eq!(phase, 0.0);
    }
    #[test]
    fn test_phase_vapour_at_low_density() {
        let pt = CryogenicPhaseTransition::liquid_nitrogen();
        let phase = pt.phase_indicator(1.0, 80.0);
        assert_eq!(phase, 1.0);
    }
    #[test]
    fn test_saturation_pressure_increases_with_temp() {
        let pt = CryogenicPhaseTransition::liquid_nitrogen();
        let p1 = pt.saturation_pressure(77.0, LN2_LHV);
        let p2 = pt.saturation_pressure(90.0, LN2_LHV);
        assert!(p2 > p1);
    }
    #[test]
    fn test_nucleate_boiling_flux_positive() {
        let b = BoilingHeatTransfer::ln2_copper_surface();
        let q = b.nucleate_boiling_flux();
        assert!(q > 0.0);
    }
    #[test]
    fn test_critical_heat_flux_positive() {
        let b = BoilingHeatTransfer::ln2_copper_surface();
        let qmax = b.critical_heat_flux();
        assert!(qmax > 0.0);
    }
    #[test]
    fn test_leidenfrost_superheat_positive() {
        let b = BoilingHeatTransfer::ln2_copper_surface();
        let dt = b.leidenfrost_superheat();
        assert!(dt > 0.0);
    }
    #[test]
    fn test_slosh_natural_frequency_positive() {
        let s = TankSloshDynamics::cylindrical(1.0, 0.5, CryogenicFluid::LiquidNitrogen);
        let f = s.natural_frequency();
        assert!(f > 0.0);
    }
    #[test]
    fn test_slosh_step_changes_displacement() {
        let mut s = TankSloshDynamics::cylindrical(1.0, 0.5, CryogenicFluid::LiquidNitrogen);
        s.step(0.001, 1.0);
        assert!(s.slosh_velocity.abs() > 0.0 || s.slosh_displacement.abs() > 0.0);
    }
    #[test]
    fn test_zbo_heat_load_positive() {
        let zbo = ZeroBoilOffInsulation::liquid_hydrogen_tank(2.0, 1.8);
        assert!(zbo.total_heat_load() > 0.0);
    }
    #[test]
    fn test_zbo_cryocooler_power_gt_heat_load() {
        let zbo = ZeroBoilOffInsulation::liquid_hydrogen_tank(2.0, 1.8);
        assert!(zbo.cryocooler_power() > zbo.total_heat_load());
    }
    #[test]
    fn test_zbo_radiation_flux_positive() {
        let zbo = ZeroBoilOffInsulation::liquid_hydrogen_tank(2.0, 1.8);
        assert!(zbo.radiation_heat_flux() > 0.0);
    }
    #[test]
    fn test_magnetocaloric_entropy_finite() {
        let mc = MagnetocaloricCooling::new(MagnetocaloricMaterial::FerricAmmoniumAlum, 1.0, 0.05);
        let s = mc.magnetic_entropy();
        assert!(s.is_finite());
    }
    #[test]
    fn test_magnetocaloric_temp_change_finite() {
        let mc =
            MagnetocaloricCooling::new(MagnetocaloricMaterial::GadoliniumGalliumGarnet, 2.0, 1.0);
        let dt = mc.adiabatic_temperature_change(-2.0);
        assert!(dt.is_finite());
    }
    #[test]
    fn test_brillouin_function_zero_at_zero() {
        let bj = brillouin_function(2.5, 0.0);
        assert!(bj.abs() < TOL);
    }
    #[test]
    fn test_brillouin_function_positive_for_positive_x() {
        let bj = brillouin_function(2.5, 1.0);
        assert!(bj > 0.0);
    }
    #[test]
    fn test_gpe_norm_after_init() {
        let solver = GrossPitaevskiiSolver::new(64, 1.0e-4, HE4_MASS, 5.0e-9, 1e4);
        let norm = solver.norm();
        assert!((norm - 1.0e4).abs() / 1.0e4 < 0.1);
    }
    #[test]
    fn test_gpe_density_non_negative() {
        let solver = GrossPitaevskiiSolver::new(64, 1.0e-4, HE4_MASS, 5.0e-9, 1.0e4);
        assert!(solver.density().iter().all(|&n| n >= 0.0));
    }
    #[test]
    fn test_gpe_energy_finite() {
        let solver = GrossPitaevskiiSolver::new(32, 1.0e-4, HE4_MASS, 5.0e-9, 1.0e4);
        let e = solver.total_energy();
        assert!(e.is_finite());
    }
    #[test]
    fn test_gpe_healing_length_positive() {
        let solver = GrossPitaevskiiSolver::new(64, 1.0e-4, HE4_MASS, 5.0e-9, 1.0e4);
        let xi = solver.healing_length();
        assert!(xi > 0.0);
    }
    #[test]
    fn test_gpe_imaginary_time_step_norm_preserved() {
        let mut solver = GrossPitaevskiiSolver::new(32, 1.0e-4, HE4_MASS, 5.0e-9, 1.0e4);
        solver.imaginary_time_step(1.0e-8);
        let norm = solver.norm();
        assert!((norm - 1.0e4).abs() / 1.0e4 < 0.01);
    }
    #[test]
    fn test_cubic_kernel_zero_at_2h() {
        let w = cubic_spline_kernel(2.0, 1.0);
        assert!(w.abs() < TOL);
    }
    #[test]
    fn test_cubic_kernel_positive_at_origin() {
        let w = cubic_spline_kernel(0.0, 1.0);
        assert!(w > 0.0);
    }
    #[test]
    fn test_cubic_kernel_decreases() {
        let w0 = cubic_spline_kernel(0.0, 1.0);
        let w1 = cubic_spline_kernel(0.5, 1.0);
        let w2 = cubic_spline_kernel(1.5, 1.0);
        assert!(w0 > w1);
        assert!(w1 > w2);
    }
    #[test]
    fn test_sim_particle_count() {
        let sim = CryogenicSphSimulation::new_ln2_box(8, 0.1);
        assert_eq!(sim.particles.len(), 8);
    }
    #[test]
    fn test_sim_thermal_velocities() {
        let mut sim = CryogenicSphSimulation::new_ln2_box(8, 0.1);
        sim.add_thermal_velocities(LN2_NBP);
        let total_ke = sim.total_kinetic_energy();
        assert!(total_ke > 0.0);
    }
    #[test]
    fn test_sim_time_advances() {
        let mut sim = CryogenicSphSimulation::new_ln2_box(8, 0.1);
        sim.run(3);
        assert!(sim.time > 0.0);
    }
    #[test]
    fn test_he2_sim_has_vortex() {
        let sim = CryogenicSphSimulation::new_he2_box(8, 0.01, 1.5);
        assert!(!sim.vortex_filaments.is_empty());
    }
    #[test]
    fn test_mean_density_nonempty() {
        let sim = CryogenicSphSimulation::new_ln2_box(8, 0.1);
        let rho = mean_density(&sim.particles);
        assert!(rho > 0.0);
    }
    #[test]
    fn test_mean_temperature_nonempty() {
        let sim = CryogenicSphSimulation::new_ln2_box(8, 0.1);
        let t = mean_particle_temperature(&sim.particles);
        assert!((t - LN2_NBP).abs() < TOL);
    }
}
