//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
use std::f64::consts::PI;

use super::types::{ColorParticle, MultiphaseFluidParticle, MultiphaseParticle};

/// Compute the smoothed color function gradient for all particles using SPH.
///
/// ∇C_i = Σⱼ (mⱼ/ρⱼ) · (Cⱼ − Cᵢ) · ∇W_ij
pub fn compute_color_gradients(particles: &mut Vec<ColorParticle>, h: f64) {
    let n = particles.len();
    let positions: Vec<[f64; 2]> = particles.iter().map(|p| p.pos).collect();
    let colors: Vec<f64> = particles.iter().map(|p| p.color).collect();
    let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let densities: Vec<f64> = particles.iter().map(|p| p.density).collect();
    let mut grads = vec![[0.0f64; 2]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let r = (dx * dx + dy * dy).sqrt();
            if r >= h {
                continue;
            }
            let (gwx, gwy) = cubic_spline_gradient(dx, dy, r, h);
            let w = masses[j] / densities[j].max(1e-30) * (colors[j] - colors[i]);
            grads[i][0] += w * gwx;
            grads[i][1] += w * gwy;
        }
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.color_grad = grads[i];
        p.update_normal();
    }
}
/// Cubic spline kernel gradient in 2D.
pub(super) fn cubic_spline_gradient(dx: f64, dy: f64, r: f64, h: f64) -> (f64, f64) {
    if r < 1e-14 {
        return (0.0, 0.0);
    }
    let q = r / h;
    let factor = 10.0 / (7.0 * PI * h * h);
    let dw_dr = if q < 1.0 {
        factor * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        factor * (-0.75 * (2.0 - q).powi(2)) / h
    } else {
        0.0
    };
    (dw_dr * dx / r, dw_dr * dy / r)
}
/// Cubic spline kernel value in 2D.
pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let factor = 10.0 / (7.0 * PI * h * h);
    if q < 1.0 {
        factor * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        factor * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}
/// Compute the colour function `C_i` for each particle using an SPH kernel sum.
///
/// `C_i = Σ_j (m_j / ρ_j) W(|r_ij|, h) × c_j`
/// where `c_j = 1` if particle j belongs to phase 0, 0 otherwise.
///
/// A simplified version is used here: each particle is assigned `color = 1`
/// if it belongs to phase 0, then the value is smoothed over the neighbourhood.
///
/// # Arguments
/// * `particles` – slice of [`MultiphaseParticle`]
/// * `h`         – smoothing length \[m\]
///
/// # Returns
/// Vector of colour values (one per particle, values in \[0, 1\]).
pub fn compute_color_function(particles: &[MultiphaseParticle], h: f64) -> Vec<f64> {
    let n = particles.len();
    let mut colors = vec![0.0_f64; n];
    for i in 0..n {
        let mut w_sum = 0.0_f64;
        let mut wc_sum = 0.0_f64;
        for j in 0..n {
            let r = particles[i].dist(&particles[j]);
            if r < h {
                let w = cubic_spline_kernel(r, h);
                let c_j = if particles[j].phase_id == 0 { 1.0 } else { 0.0 };
                w_sum += w;
                wc_sum += w * c_j;
            }
        }
        colors[i] = if w_sum > 1e-14 {
            wc_sum / w_sum
        } else if particles[i].phase_id == 0 {
            1.0
        } else {
            0.0
        };
    }
    colors
}
/// Compute interface curvature `κ_i` from colour-function values.
///
/// Uses a finite-difference Laplacian approximation on the ordered colour
/// array:
///
/// ```text
/// κ_i ≈ (C_{i+1} - 2 C_i + C_{i-1}) / Δs²
/// ```
///
/// This is a 1-D approximation; for 3-D SPH the divergence of the normalised
/// colour gradient would be used.
///
/// # Arguments
/// * `color` – colour-function values ordered by particle index
///
/// # Returns
/// Vector of curvature values \[1/m\] (same length as `color`).
pub fn compute_curvature(color: &[f64]) -> Vec<f64> {
    let n = color.len();
    if n == 0 {
        return Vec::new();
    }
    let mut kappa = vec![0.0_f64; n];
    for i in 1..(n.saturating_sub(1)) {
        kappa[i] = color[i + 1] - 2.0 * color[i] + color[i - 1];
    }
    if n >= 2 {
        kappa[0] = color[1] - color[0];
        kappa[n - 1] = color[n - 1] - color[n - 2];
    }
    kappa
}
/// SPH kernel gradient magnitude |∇W(r, h)|.
///
/// Derivative of the cubic spline used in [`compute_color_function`].
pub(super) fn sph_kernel_gradient_magnitude(r: f64, h: f64) -> f64 {
    let q = r / h;
    let norm = 10.0 / (7.0 * std::f64::consts::PI * h * h);
    if q < 1.0 {
        norm * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        norm * (-0.75 * (2.0 - q).powi(2)) / h
    } else {
        0.0
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiphase_sph::types::*;
    #[test]
    fn phase_properties_density_ratio() {
        let water = PhaseProperties::water();
        let air = PhaseProperties::air();
        let ratio = water.density_ratio(&air);
        assert!((ratio - 1000.0 / 1.225).abs() < 1.0);
    }
    #[test]
    fn atwood_number_water_air() {
        let water = PhaseProperties::water();
        let air = PhaseProperties::air();
        let at = PhaseProperties::atwood_number(&water, &air);
        assert!(
            at > 0.99,
            "water/air Atwood number should be close to 1: {at}"
        );
    }
    #[test]
    fn color_particle_liquid_phase() {
        let p = ColorParticle::new_liquid(0.0, 0.0, 1e-3);
        assert_eq!(p.phase(), Phase::Liquid);
        assert!(!p.is_interface());
    }
    #[test]
    fn color_particle_gas_phase() {
        let p = ColorParticle::new_gas(0.0, 0.0, 1e-3);
        assert_eq!(p.phase(), Phase::Gas);
    }
    #[test]
    fn color_particle_interface_detection() {
        let mut p = ColorParticle::new_liquid(0.0, 0.0, 1e-3);
        p.color = 0.5;
        assert!(p.is_interface());
    }
    #[test]
    fn color_particle_interpolated_density() {
        let p = ColorParticle::new_liquid(0.0, 0.0, 1e-3);
        let rho = p.interpolated_density(1000.0, 1.225);
        assert!((rho - 1000.0).abs() < 1.0, "liquid color=1 → rho_l: {rho}");
    }
    #[test]
    fn color_particle_update_normal() {
        let mut p = ColorParticle::new_liquid(0.0, 0.0, 1e-3);
        p.color_grad = [3.0, 4.0];
        p.update_normal();
        let norm = (p.normal[0].powi(2) + p.normal[1].powi(2)).sqrt();
        assert!(
            (norm - 1.0).abs() < 1e-10,
            "normal should be unit vector: {norm}"
        );
    }
    #[test]
    fn cubic_spline_kernel_zero_at_cutoff() {
        let w = cubic_spline_kernel(2.0, 1.0);
        assert!(w.abs() < 1e-10, "kernel should be zero at r=2h: {w}");
    }
    #[test]
    fn cubic_spline_kernel_positive_inside() {
        let w = cubic_spline_kernel(0.5, 1.0);
        assert!(w > 0.0, "kernel should be positive inside support: {w}");
    }
    #[test]
    fn cubic_spline_kernel_max_at_zero() {
        let w0 = cubic_spline_kernel(0.0, 1.0);
        let w1 = cubic_spline_kernel(0.5, 1.0);
        assert!(w0 > w1, "kernel should peak at origin");
    }
    #[test]
    fn cahn_hilliard_tanh_profile_interface() {
        let field = CahnHilliardField::new_1d(101, 1.0, 0.1, 1e-6);
        let pos = field.interface_position();
        assert!(pos.is_some(), "tanh profile should have a zero crossing");
        let x0 = pos.unwrap();
        assert!((x0 - 0.5).abs() < 0.1, "interface near center: {x0}");
    }
    #[test]
    fn cahn_hilliard_bulk_energy_minima_at_pm1() {
        let field = CahnHilliardField::new_1d(10, 1.0, 0.1, 1e-6);
        let f_m1 = field.bulk_free_energy(-1.0);
        let f_p1 = field.bulk_free_energy(1.0);
        let f_0 = field.bulk_free_energy(0.0);
        assert!(f_m1 < f_0, "φ=±1 should be energy minima");
        assert!(f_p1 < f_0, "φ=±1 should be energy minima");
    }
    #[test]
    fn cahn_hilliard_step_no_panic() {
        let mut field = CahnHilliardField::new_1d(51, 1.0, 0.05, 1e-8);
        field.step(1e-5);
    }
    #[test]
    fn cahn_hilliard_volume_fraction_half() {
        let field = CahnHilliardField::new_1d(101, 1.0, 0.1, 1e-6);
        let vf = field.volume_fraction();
        assert!(
            (vf - 0.5).abs() < 0.1,
            "symmetric tanh should give ~50% volume: {vf}"
        );
    }
    #[test]
    fn csf_delta_zero_outside_interface() {
        let csf = CsfSurfaceTension::new(0.072, 0.01);
        let d = csf.delta(0.1);
        assert!(d.abs() < 1e-10, "delta should be zero outside: {d}");
    }
    #[test]
    fn csf_delta_positive_inside() {
        let csf = CsfSurfaceTension::new(0.072, 0.1);
        let d = csf.delta(0.05);
        assert!(d > 0.0, "delta should be positive inside interface: {d}");
    }
    #[test]
    fn csf_laplace_pressure_sphere() {
        let csf = CsfSurfaceTension::new(0.072, 0.01);
        let dp = csf.laplace_pressure_3d(1e-3);
        let expected = 2.0 * 0.072 / 1e-3;
        assert!(
            (dp - expected).abs() < 1.0,
            "Laplace pressure: {dp} vs {expected}"
        );
    }
    #[test]
    fn csf_weber_capillary_numbers() {
        let csf = CsfSurfaceTension::new(0.072, 0.01);
        let we = csf.weber_number(1000.0, 0.1, 1e-3);
        let ca = csf.capillary_number(1e-3, 0.1);
        assert!(we > 0.0);
        assert!(ca > 0.0);
        assert!(ca < 1.0, "Ca should be small for slow flow: {ca}");
    }
    #[test]
    fn boussinesq_zero_temp_diff_zero_force() {
        let model = BoussinesqModel::water(9.81);
        let f = model.buoyancy_force(model.t_ref, 0.0);
        assert!(f.abs() < 1e-8, "no temp diff → no buoyancy: {f}");
    }
    #[test]
    fn boussinesq_hot_lighter() {
        let model = BoussinesqModel::water(9.81);
        let rho_hot = model.effective_density(model.t_ref + 10.0);
        let rho_cold = model.effective_density(model.t_ref);
        assert!(rho_hot < rho_cold, "hot fluid should be lighter");
    }
    #[test]
    fn boussinesq_rayleigh_number_positive() {
        let model = BoussinesqModel::water(9.81);
        let ra = model.rayleigh_number(20.0, 0.1, 1e-6, 1.4e-7);
        assert!(ra > 0.0);
    }
    #[test]
    fn rt_atwood_number_water_air() {
        let rt = RayleighTaylorInstability {
            rho_heavy: 1000.0,
            rho_light: 1.225,
            g: 9.81,
            sigma: 0.072,
            nu_heavy: 1e-6,
            nu_light: 1.5e-5,
        };
        let at = rt.atwood_number();
        assert!(at > 0.99, "water/air At ≈ 1: {at}");
    }
    #[test]
    fn rt_growth_rate_positive_unstable() {
        let rt = RayleighTaylorInstability {
            rho_heavy: 2000.0,
            rho_light: 1000.0,
            g: 9.81,
            sigma: 0.0,
            nu_heavy: 1e-6,
            nu_light: 1e-6,
        };
        let gamma = rt.growth_rate_inviscid(100.0);
        assert!(
            gamma > 0.0,
            "unstable interface should have positive growth rate: {gamma}"
        );
    }
    #[test]
    fn rt_surface_tension_stabilizes_short_waves() {
        let rt = RayleighTaylorInstability {
            rho_heavy: 2000.0,
            rho_light: 1000.0,
            g: 9.81,
            sigma: 1.0,
            nu_heavy: 1e-6,
            nu_light: 1e-6,
        };
        let kc = rt.critical_wavenumber();
        let gamma_above = rt.growth_rate_with_tension(kc * 2.0);
        assert_eq!(gamma_above, 0.0, "waves above k_c should be stabilized");
    }
    #[test]
    fn rayleigh_plesset_initial_equilibrium() {
        let rp = RayleighPlesset::new(1e-4, 1e5, 1000.0, 1e-3, 0.072);
        let rdd = rp.radius_ddot();
        assert!(rdd.abs() < 1e12, "initial R̈ should be bounded: {rdd}");
    }
    #[test]
    fn rayleigh_plesset_minnaert_frequency_positive() {
        let rp = RayleighPlesset::new(1e-3, 1e5, 1000.0, 1e-3, 0.072);
        let omega = rp.minnaert_frequency();
        assert!(
            omega > 0.0,
            "Minnaert frequency should be positive: {omega}"
        );
    }
    #[test]
    fn rayleigh_plesset_collapse_time_positive() {
        let rp = RayleighPlesset::new(1e-3, 1e5, 1000.0, 1e-3, 0.072);
        let tc = rp.collapse_time();
        assert!(tc > 0.0 && tc < 1.0, "collapse time: {tc}");
    }
    #[test]
    fn rayleigh_plesset_step_no_nan() {
        let mut rp = RayleighPlesset::new(1e-3, 1e5, 1000.0, 1e-3, 0.072);
        rp.step(1e-8);
        assert!(!rp.radius.is_nan() && !rp.radius_dot.is_nan());
    }
    #[test]
    fn droplet_overlap_detection() {
        let d1 = Droplet::new(0, 0.0, 0.0, 1.0, 1000.0, 0.072);
        let d2 = Droplet::new(1, 1.5, 0.0, 1.0, 1000.0, 0.072);
        assert!(
            d1.overlaps(&d2),
            "droplets at distance 1.5 < 2.0 should overlap"
        );
        let d3 = Droplet::new(2, 3.0, 0.0, 1.0, 1000.0, 0.072);
        assert!(
            !d1.overlaps(&d3),
            "droplets at distance 3.0 > 2.0 should not overlap"
        );
    }
    #[test]
    fn droplet_coalescence_conserves_volume() {
        let d1 = Droplet::new(0, 0.0, 0.0, 1.0, 1000.0, 0.072);
        let d2 = Droplet::new(1, 2.5, 0.0, 1.0, 1000.0, 0.072);
        let merged = Droplet::coalesce(&d1, &d2);
        let vol_before = d1.volume + d2.volume;
        assert!(
            (merged.volume - vol_before).abs() / vol_before < 1e-8,
            "volume should be conserved: {} vs {}",
            merged.volume,
            vol_before
        );
    }
    #[test]
    fn droplet_breakup_conserves_volume() {
        let d = Droplet::new(0, 0.0, 0.0, 1.0, 1000.0, 0.072);
        let (d1, d2) = d.breakup([1.0, 0.0]);
        let vol_after = d1.volume + d2.volume;
        assert!(
            (vol_after - d.volume).abs() / d.volume < 1e-8,
            "breakup should conserve volume: {vol_after} vs {}",
            d.volume
        );
    }
    #[test]
    fn droplet_kinetic_energy_zero_at_rest() {
        let d = Droplet::new(0, 0.0, 0.0, 0.5, 1000.0, 0.072);
        assert!(d.kinetic_energy().abs() < 1e-20);
    }
    #[test]
    fn contact_angle_young_equation() {
        let model = ContactAngleModel::new(0.072, 0.040, 0.030);
        let cos_theta = (0.040 - 0.030) / 0.072;
        assert!((model.theta_y.cos() - cos_theta).abs() < 1e-8);
    }
    #[test]
    fn contact_angle_hydrophilic_water_glass() {
        let model = ContactAngleModel::new(0.072, 0.070, 0.030);
        assert!(model.is_hydrophilic(), "glass should be hydrophilic");
    }
    #[test]
    fn contact_angle_capillary_length() {
        let model = ContactAngleModel::new(0.072, 0.040, 0.030);
        let lc = model.capillary_length(1000.0, 9.81);
        let expected = (0.072_f64 / (1000.0_f64 * 9.81_f64)).sqrt();
        assert!((lc - expected).abs() < 1e-8, "capillary length: {lc}");
    }
    #[test]
    fn contact_angle_work_of_adhesion_positive() {
        let model = ContactAngleModel::new(0.072, 0.040, 0.030);
        let w = model.work_of_adhesion();
        assert!(w > 0.0, "work of adhesion should be positive: {w}");
    }
    #[test]
    fn vof_step_init_conserves_volume() {
        let mut vof = VolumeFractionField::new_1d(101, 1.0, 0.0);
        vof.step_init(0.5);
        let v1 = vof.total_volume();
        vof.set_uniform_velocity(0.0);
        vof.advect_step(1e-3);
        let v2 = vof.total_volume();
        assert!(
            (v1 - v2).abs() / v1 < 1e-8,
            "zero velocity should conserve volume: {v1} vs {v2}"
        );
    }
    #[test]
    fn vof_interface_location_step() {
        let mut vof = VolumeFractionField::new_1d(101, 1.0, 0.0);
        vof.step_init(0.5);
        let pos = vof.interface_location();
        assert!(pos.is_some(), "step function should have interface");
        let x0 = pos.unwrap();
        assert!((x0 - 0.5).abs() < 0.02, "interface at 0.5: {x0}");
    }
    #[test]
    fn multiphase_sim_counts() {
        let mut sim = MultiphaseSphSimulation::new(0.1, 0.072, 1000.0, 1.225);
        sim.add_liquid_particle(0.0, 0.0, 1e-3);
        sim.add_liquid_particle(0.1, 0.0, 1e-3);
        sim.add_gas_particle(0.0, 0.1, 1e-4);
        assert_eq!(sim.n_liquid(), 2);
        assert_eq!(sim.n_gas(), 1);
    }
    #[test]
    fn multiphase_sim_gravity_accelerates_liquid_down() {
        let mut sim = MultiphaseSphSimulation::new(0.1, 0.072, 1000.0, 1.225);
        sim.add_liquid_particle(0.5, 1.0, 1e-3);
        sim.dt = 0.001;
        sim.step();
        let vy = sim.particles[0].vel[1];
        assert!(vy < 0.0, "gravity should push particle downward: vy = {vy}");
    }
    #[test]
    fn multiphase_sim_kinetic_energy_increases_under_gravity() {
        let mut sim = MultiphaseSphSimulation::new(0.1, 0.072, 1000.0, 1.225);
        sim.add_liquid_particle(0.5, 1.0, 1e-3);
        sim.dt = 0.001;
        let ke0 = sim.kinetic_energy();
        for _ in 0..10 {
            sim.step();
        }
        let ke1 = sim.kinetic_energy();
        assert!(ke1 > ke0, "KE should increase under gravity");
    }
    #[test]
    fn phase_desc_id_stored() {
        let p = PhaseDesc::new(2, 1000.0, 1e-3, 0.072);
        assert_eq!(p.id, 2);
    }
    #[test]
    fn phase_desc_density_stored() {
        let p = PhaseDesc::new(0, 800.0, 0.01, 0.03);
        assert_eq!(p.density, 800.0);
    }
    #[test]
    fn phase_desc_capillary_length_positive() {
        let p = PhaseDesc::new(0, 1000.0, 1e-3, 0.072);
        assert!(p.capillary_length(9.81) > 0.0);
    }
    #[test]
    fn phase_desc_capillary_length_water_approx_2_7mm() {
        let p = PhaseDesc::new(0, 1000.0, 1e-3, 0.072);
        let lc = p.capillary_length(9.81);
        assert!(
            (lc - 0.00271).abs() < 0.0001,
            "capillary length ≈ 2.71 mm, got {lc}"
        );
    }
    #[test]
    fn multiphase_particle_phase0_color_one() {
        let p = MultiphaseParticle::new([0.0, 0.0, 0.0], 0, 1e-3);
        assert_eq!(p.color_function, 1.0);
    }
    #[test]
    fn multiphase_particle_phase1_color_zero() {
        let p = MultiphaseParticle::new([0.0, 0.0, 0.0], 1, 1e-3);
        assert_eq!(p.color_function, 0.0);
    }
    #[test]
    fn multiphase_particle_at_rest_zero_ke() {
        let p = MultiphaseParticle::new([1.0, 2.0, 3.0], 0, 0.5);
        assert_eq!(p.kinetic_energy(), 0.0);
    }
    #[test]
    fn multiphase_particle_dist_self_zero() {
        let p = MultiphaseParticle::new([1.0, 2.0, 3.0], 0, 1e-3);
        assert_eq!(p.dist(&p.clone()), 0.0);
    }
    #[test]
    fn multiphase_particle_dist_pythagorean() {
        let p1 = MultiphaseParticle::new([0.0, 0.0, 0.0], 0, 1e-3);
        let p2 = MultiphaseParticle::new([3.0, 4.0, 0.0], 1, 1e-3);
        assert!((p1.dist(&p2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn interface_model_equality() {
        assert_eq!(InterfaceModel::ColorFunction, InterfaceModel::ColorFunction);
        assert_ne!(InterfaceModel::ColorFunction, InterfaceModel::PhaseField);
    }
    #[test]
    fn surface_tension_force_equality() {
        assert_eq!(SurfaceTensionForce::Csf, SurfaceTensionForce::Csf);
        assert_ne!(SurfaceTensionForce::Csf, SurfaceTensionForce::Pairwise);
    }
    #[test]
    fn multiphase_sph_add_phases_and_particles() {
        let mut sim = MultiphaseSPH::new(0.1, [0.0, -9.81, 0.0], 0.001);
        let id0 = sim.add_phase(1000.0, 1e-3, 0.072);
        let id1 = sim.add_phase(1.225, 1.8e-5, 0.072);
        sim.add_particle([0.0, 0.0, 0.0], id0, 1e-3);
        sim.add_particle([0.2, 0.2, 0.0], id1, 1e-5);
        assert_eq!(sim.n_particles(), 2);
    }
    #[test]
    fn multiphase_sph_phase_ids_sequential() {
        let mut sim = MultiphaseSPH::new(0.1, [0.0, -9.81, 0.0], 0.001);
        let id0 = sim.add_phase(1000.0, 1e-3, 0.072);
        let id1 = sim.add_phase(800.0, 0.01, 0.03);
        assert_eq!(id0, 0);
        assert_eq!(id1, 1);
    }
    #[test]
    fn multiphase_sph_step_no_panic() {
        let mut sim = MultiphaseSPH::new(0.15, [0.0, -9.81, 0.0], 0.001);
        let id0 = sim.add_phase(1000.0, 1e-3, 0.072);
        let id1 = sim.add_phase(1.225, 1.8e-5, 0.0);
        sim.add_particle([0.0, 0.0, 0.0], id0, 1e-3);
        sim.add_particle([0.05, 0.05, 0.0], id1, 1e-5);
        sim.step();
        assert!(sim.particles[0].position[1] != 0.0 || sim.particles[0].velocity[1] != 0.0);
    }
    #[test]
    fn multiphase_sph_gravity_accelerates_particle() {
        let mut sim = MultiphaseSPH::new(0.1, [0.0, -9.81, 0.0], 0.001);
        let id0 = sim.add_phase(1000.0, 1e-3, 0.0);
        sim.add_particle([0.5, 1.0, 0.0], id0, 1e-3);
        sim.step();
        assert!(
            sim.particles[0].velocity[1] < 0.0,
            "gravity should accelerate downward"
        );
    }
    #[test]
    fn multiphase_sph_kinetic_energy_increases() {
        let mut sim = MultiphaseSPH::new(0.1, [0.0, -9.81, 0.0], 0.001);
        let id0 = sim.add_phase(1000.0, 1e-3, 0.0);
        sim.add_particle([0.5, 1.0, 0.0], id0, 1e-3);
        let ke0 = sim.kinetic_energy();
        for _ in 0..5 {
            sim.step();
        }
        let ke1 = sim.kinetic_energy();
        assert!(ke1 > ke0);
    }
    #[test]
    fn color_function_single_phase0_particle() {
        let particles = vec![MultiphaseParticle::new([0.0, 0.0, 0.0], 0, 1e-3)];
        let c = compute_color_function(&particles, 0.5);
        assert!((c[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn color_function_single_phase1_particle() {
        let particles = vec![MultiphaseParticle::new([0.0, 0.0, 0.0], 1, 1e-3)];
        let c = compute_color_function(&particles, 0.5);
        assert!((c[0]).abs() < 1e-10);
    }
    #[test]
    fn color_function_length_matches_particles() {
        let particles = vec![
            MultiphaseParticle::new([0.0, 0.0, 0.0], 0, 1e-3),
            MultiphaseParticle::new([1.0, 0.0, 0.0], 1, 1e-3),
        ];
        let c = compute_color_function(&particles, 0.2);
        assert_eq!(c.len(), 2);
    }
    #[test]
    fn color_function_empty_returns_empty() {
        let c = compute_color_function(&[], 0.1);
        assert!(c.is_empty());
    }
    #[test]
    fn color_function_values_in_range() {
        let particles = vec![
            MultiphaseParticle::new([0.0, 0.0, 0.0], 0, 1e-3),
            MultiphaseParticle::new([0.05, 0.0, 0.0], 0, 1e-3),
            MultiphaseParticle::new([0.1, 0.0, 0.0], 1, 1e-3),
            MultiphaseParticle::new([0.15, 0.0, 0.0], 1, 1e-3),
        ];
        let c = compute_color_function(&particles, 0.12);
        for &ci in &c {
            assert!(
                (0.0..=1.0).contains(&ci) || (ci - 1.0).abs() < 1e-10 || ci.abs() < 1e-10,
                "color value {ci} out of range"
            );
        }
    }
    #[test]
    fn curvature_empty_returns_empty() {
        let k = compute_curvature(&[]);
        assert!(k.is_empty());
    }
    #[test]
    fn curvature_single_element_returns_zero() {
        let k = compute_curvature(&[0.5]);
        assert_eq!(k.len(), 1);
        assert_eq!(k[0], 0.0);
    }
    #[test]
    fn curvature_length_matches_input() {
        let colors = vec![1.0, 0.8, 0.5, 0.2, 0.0];
        let k = compute_curvature(&colors);
        assert_eq!(k.len(), 5);
    }
    #[test]
    fn curvature_flat_interface_zero() {
        let colors: Vec<f64> = (0..5).map(|i| i as f64 * 0.25).collect();
        let k = compute_curvature(&colors);
        for &ki in &k[1..k.len() - 1] {
            assert!(
                ki.abs() < 1e-12,
                "linear profile should give zero curvature, got {ki}"
            );
        }
    }
    #[test]
    fn curvature_parabola_constant() {
        let colors: Vec<f64> = (0..6).map(|i| (i * i) as f64).collect();
        let k = compute_curvature(&colors);
        for &ki in &k[1..k.len() - 1] {
            assert!(
                (ki - 2.0).abs() < 1e-10,
                "parabola should give curvature 2, got {ki}"
            );
        }
    }
}
#[allow(dead_code)]
pub(super) fn cubic_kernel_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        alpha * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}
#[allow(dead_code)]
pub(super) fn cubic_kernel_grad_3d(rij: [f64; 3], r: f64, h: f64) -> [f64; 3] {
    if r < 1e-30 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dq = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        alpha * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    };
    let dw_dr = dw_dq / h;
    [dw_dr * rij[0] / r, dw_dr * rij[1] / r, dw_dr * rij[2] / r]
}
/// Compute the SPH colour-function estimate at particle `i`.
///
/// `C_i = Σ_j (m_j / ρ_j) C_j W(|x_i - x_j|, h)`
///
/// # Arguments
/// * `particles` – all particles
/// * `i`         – index of the query particle
/// * `h`         – kernel smoothing length \[m\]
pub fn color_function(particles: &[MultiphaseFluidParticle], i: usize, h: f64) -> f64 {
    let xi = particles[i].position;
    let mut c_sum = 0.0_f64;
    for p in particles {
        let dx = xi[0] - p.position[0];
        let dy = xi[1] - p.position[1];
        let dz = xi[2] - p.position[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let w = cubic_kernel_3d(r, h);
        c_sum += (p.mass / p.density.max(1e-30)) * p.color_field * w;
    }
    c_sum
}
/// Estimate the interface normal at particle `i` as ∇C_i (not normalized).
///
/// # Arguments
/// * `particles` – all particles
/// * `i`         – query particle index
/// * `h`         – smoothing length \[m\]
pub fn interface_normal(particles: &[MultiphaseFluidParticle], i: usize, h: f64) -> [f64; 3] {
    let xi = particles[i].position;
    let mut grad = [0.0_f64; 3];
    for p in particles {
        let rij = [
            xi[0] - p.position[0],
            xi[1] - p.position[1],
            xi[2] - p.position[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        let dw = cubic_kernel_grad_3d(rij, r, h);
        let vol_j = p.mass / p.density.max(1e-30);
        grad[0] += vol_j * p.color_field * dw[0];
        grad[1] += vol_j * p.color_field * dw[1];
        grad[2] += vol_j * p.color_field * dw[2];
    }
    grad
}
/// Compute the interface curvature κ = ∇ · n̂ from the gradient of the normal.
///
/// `κ = tr(∂n̂/∂x) = Σ_d (∂n̂_d/∂x_d)` approximated from `grad_normal`.
///
/// # Arguments
/// * `normal`      – unit normal vector n̂
/// * `grad_normal` – Jacobian of n̂: `grad_normal[d][k] = ∂n̂_d/∂x_k`
pub fn interface_curvature(normal: [f64; 3], grad_normal: [[f64; 3]; 3]) -> f64 {
    let n_mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if n_mag < 1e-30 {
        return 0.0;
    }
    -(grad_normal[0][0] + grad_normal[1][1] + grad_normal[2][2])
}
/// Pressure jump across a curved interface: `Δp = σ κ` (Young-Laplace).
///
/// # Arguments
/// * `curvature` – interface curvature κ \[1/m\]
/// * `sigma`     – surface-tension coefficient σ \[N/m\]
pub fn pressure_jump_young_laplace(curvature: f64, sigma: f64) -> f64 {
    sigma * curvature
}
/// Compute the Continuum Surface Force (CSF) surface-tension force on particle `i`.
///
/// `f_st = σ κ_i n_i / ρ_i`  (force per unit mass → acceleration contribution)
///
/// # Arguments
/// * `particles` – all particles
/// * `i`         – query particle index
/// * `h`         – smoothing length \[m\]
/// * `sigma`     – surface-tension coefficient \[N/m\]
pub fn surface_tension_csf(
    particles: &[MultiphaseFluidParticle],
    i: usize,
    h: f64,
    sigma: f64,
) -> [f64; 3] {
    let n_i = interface_normal(particles, i, h);
    let n_mag = (n_i[0] * n_i[0] + n_i[1] * n_i[1] + n_i[2] * n_i[2]).sqrt();
    if n_mag < 1e-30 {
        return [0.0; 3];
    }
    let xi = particles[i].position;
    let mut div_n = 0.0_f64;
    for p in particles {
        let rij = [
            xi[0] - p.position[0],
            xi[1] - p.position[1],
            xi[2] - p.position[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        let dw = cubic_kernel_grad_3d(rij, r, h);
        let vol_j = p.mass / p.density.max(1e-30);
        let n_j = [
            p.color_field * rij[0].signum(),
            p.color_field * rij[1].signum(),
            0.0,
        ];
        div_n += vol_j * (n_j[0] * dw[0] + n_j[1] * dw[1] + n_j[2] * dw[2]);
    }
    let kappa = -div_n;
    let rho_i = particles[i].density.max(1e-30);
    let coeff = sigma * kappa / rho_i;
    let n_hat = [n_i[0] / n_mag, n_i[1] / n_mag, n_i[2] / n_mag];
    [coeff * n_hat[0], coeff * n_hat[1], coeff * n_hat[2]]
}
/// Compute the volume fraction of particles with the given phase id.
///
/// Returns N_phase / N_total, or 0.0 if there are no particles.
///
/// # Arguments
/// * `particles` – all particles
/// * `phase_id`  – target phase identifier
pub fn phase_fraction(particles: &[MultiphaseFluidParticle], phase_id: u32) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let count = particles.iter().filter(|p| p.phase == phase_id).count();
    count as f64 / particles.len() as f64
}
/// Classify a particle's phase based on its density relative to a threshold.
///
/// Returns `0` if `density ≥ density_threshold`, `1` otherwise.
///
/// # Arguments
/// * `density`           – particle density \[kg/m³\]
/// * `density_threshold` – threshold separating heavy (phase 0) from light (phase 1)
pub fn classify_particle_phase(density: f64, density_threshold: f64) -> u32 {
    if density >= density_threshold { 0 } else { 1 }
}
/// Multiphase equation of state: Tait-type EOS per phase.
///
/// Phase 0: `p = k1 * ((ρ/ρ₀)^γ - 1)`
/// Phase 1: `p = k2 * ((ρ/ρ₀)^γ - 1)`
/// where `ρ₀ = 1000.0` \[kg/m³\] (reference).
///
/// # Arguments
/// * `density`  – particle density ρ \[kg/m³\]
/// * `phase`    – phase identifier
/// * `k1`       – stiffness for phase 0
/// * `k2`       – stiffness for phase 1
/// * `gamma`    – polytropic exponent γ
pub fn multiphase_eos(density: f64, phase: u32, k1: f64, k2: f64, gamma: f64) -> f64 {
    let rho0 = 1000.0_f64;
    let k = if phase == 0 { k1 } else { k2 };
    k * ((density / rho0).powf(gamma) - 1.0)
}
#[cfg(test)]
mod mp_ext_tests {
    use super::*;
    use crate::multiphase_sph::types::*;
    fn make_two_phase_particles() -> Vec<MultiphaseFluidParticle> {
        let mut v = Vec::new();
        for i in 0..5 {
            let mut p = MultiphaseFluidParticle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 1000.0, 0);
            p.color_field = 1.0;
            p.density = 1000.0;
            v.push(p);
        }
        for i in 5..10 {
            let mut p =
                MultiphaseFluidParticle::new([i as f64 * 0.1, 0.0, 0.0], 0.001225, 1.225, 1);
            p.color_field = 0.0;
            p.density = 1.225;
            v.push(p);
        }
        v
    }
    #[test]
    fn test_multiphase_fluid_particle_fields() {
        let p = MultiphaseFluidParticle::new([1.0, 2.0, 3.0], 0.5, 1000.0, 0);
        assert!((p.position[0] - 1.0).abs() < 1e-12);
        assert!((p.mass - 0.5).abs() < 1e-12);
        assert_eq!(p.phase, 0);
        assert!((p.color_field - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_multiphase_fluid_particle_phase1_color() {
        let p = MultiphaseFluidParticle::new([0.0, 0.0, 0.0], 1.0, 1.0, 1);
        assert!((p.color_field - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_multiphase_params_fields() {
        let mp = MultiphaseParams::new(0.07, 0.01, 800.0);
        assert!((mp.surface_tension_coeff - 0.07).abs() < 1e-12);
        assert!((mp.interface_width - 0.01).abs() < 1e-12);
        assert!((mp.density_ratio - 800.0).abs() < 1e-12);
    }
    #[test]
    fn test_color_function_single_particle() {
        let p = MultiphaseFluidParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0, 0);
        let particles = vec![p];
        let c = color_function(&particles, 0, 0.1);
        assert!(c >= 0.0);
    }
    #[test]
    fn test_color_function_phase0_higher() {
        let particles = make_two_phase_particles();
        let c0 = color_function(&particles, 0, 0.5);
        let c9 = color_function(&particles, 9, 0.5);
        assert!(c0 >= c9, "phase-0 color {c0} should >= phase-1 color {c9}");
    }
    #[test]
    fn test_interface_normal_zero_same_phase() {
        let particles: Vec<MultiphaseFluidParticle> = (0..4)
            .map(|i| {
                let mut p =
                    MultiphaseFluidParticle::new([i as f64 * 0.01, 0.0, 0.0], 1.0, 1000.0, 0);
                p.color_field = 1.0;
                p
            })
            .collect();
        let n = interface_normal(&particles, 1, 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(mag.is_finite());
    }
    #[test]
    fn test_interface_normal_returns_3d() {
        let particles = make_two_phase_particles();
        let n = interface_normal(&particles, 4, 0.3);
        assert_eq!(n.len(), 3);
        assert!(n[0].is_finite() && n[1].is_finite() && n[2].is_finite());
    }
    #[test]
    fn test_interface_curvature_zero_normal() {
        let kappa = interface_curvature([0.0, 0.0, 0.0], [[0.0; 3]; 3]);
        assert_eq!(kappa, 0.0);
    }
    #[test]
    fn test_interface_curvature_from_divergence() {
        let n = [1.0, 0.0, 0.0];
        let gn = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let kappa = interface_curvature(n, gn);
        assert!((kappa + 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_pressure_jump_young_laplace_positive() {
        let dp = pressure_jump_young_laplace(10.0, 0.07);
        assert!((dp - 0.7).abs() < 1e-12);
    }
    #[test]
    fn test_pressure_jump_young_laplace_zero_curvature() {
        let dp = pressure_jump_young_laplace(0.0, 0.07);
        assert_eq!(dp, 0.0);
    }
    #[test]
    fn test_pressure_jump_scales_with_sigma() {
        let dp1 = pressure_jump_young_laplace(5.0, 1.0);
        let dp2 = pressure_jump_young_laplace(5.0, 2.0);
        assert!((dp2 - 2.0 * dp1).abs() < 1e-12);
    }
    #[test]
    fn test_phase_fraction_half() {
        let particles = make_two_phase_particles();
        let f0 = phase_fraction(&particles, 0);
        let f1 = phase_fraction(&particles, 1);
        assert!((f0 - 0.5).abs() < 1e-12);
        assert!((f1 - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_phase_fraction_all_same() {
        let particles: Vec<MultiphaseFluidParticle> = (0..5)
            .map(|_| MultiphaseFluidParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0, 0))
            .collect();
        assert!((phase_fraction(&particles, 0) - 1.0).abs() < 1e-12);
        assert!((phase_fraction(&particles, 1) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_phase_fraction_empty() {
        let particles: Vec<MultiphaseFluidParticle> = vec![];
        assert_eq!(phase_fraction(&particles, 0), 0.0);
    }
    #[test]
    fn test_classify_particle_phase_above() {
        assert_eq!(classify_particle_phase(1000.0, 500.0), 0);
    }
    #[test]
    fn test_classify_particle_phase_below() {
        assert_eq!(classify_particle_phase(100.0, 500.0), 1);
    }
    #[test]
    fn test_classify_particle_phase_at_threshold() {
        assert_eq!(classify_particle_phase(500.0, 500.0), 0);
    }
    #[test]
    fn test_multiphase_eos_reference_density() {
        let p0 = multiphase_eos(1000.0, 0, 1e5, 1e4, 7.0);
        assert!(p0.abs() < 1e-8);
    }
    #[test]
    fn test_multiphase_eos_phase0_vs_phase1() {
        let p0 = multiphase_eos(1100.0, 0, 1e5, 1e4, 7.0);
        let p1 = multiphase_eos(1100.0, 1, 1e5, 1e4, 7.0);
        assert!((p0 / p1 - 10.0).abs() < 1e-5);
    }
    #[test]
    fn test_multiphase_eos_increases_with_density() {
        let p1 = multiphase_eos(1000.0, 0, 1e5, 1e4, 7.0);
        let p2 = multiphase_eos(1100.0, 0, 1e5, 1e4, 7.0);
        assert!(p2 > p1);
    }
    #[test]
    fn test_surface_tension_csf_finite() {
        let particles = make_two_phase_particles();
        let f = surface_tension_csf(&particles, 4, 0.3, 0.07);
        assert!(f[0].is_finite() && f[1].is_finite() && f[2].is_finite());
    }
    #[test]
    fn test_surface_tension_csf_single_phase_small() {
        let particles: Vec<MultiphaseFluidParticle> = (0..5)
            .map(|i| {
                let mut p =
                    MultiphaseFluidParticle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 1000.0, 0);
                p.color_field = 1.0;
                p
            })
            .collect();
        let f = surface_tension_csf(&particles, 2, 0.5, 0.07);
        assert!(f[0].is_finite());
    }
    #[test]
    fn test_color_function_empty_neighborhood() {
        let mut particles = vec![
            MultiphaseFluidParticle::new([0.0, 0.0, 0.0], 1.0, 1000.0, 0),
            MultiphaseFluidParticle::new([100.0, 0.0, 0.0], 1.0, 1000.0, 1),
        ];
        particles[0].color_field = 1.0;
        particles[1].color_field = 0.0;
        let c = color_function(&particles, 0, 0.1);
        assert!(c >= 0.0 && c.is_finite());
    }
    #[test]
    fn test_multiphase_eos_negative_density_gives_negative_pressure() {
        let p = multiphase_eos(100.0, 0, 1e5, 1e4, 7.0);
        assert!(p < 0.0);
    }
    #[test]
    fn test_classify_many_particles() {
        let densities = [200.0, 600.0, 1000.0, 1200.0, 400.0];
        let threshold = 500.0;
        let classes: Vec<u32> = densities
            .iter()
            .map(|&d| classify_particle_phase(d, threshold))
            .collect();
        assert_eq!(classes[0], 1);
        assert_eq!(classes[1], 0);
        assert_eq!(classes[2], 0);
        assert_eq!(classes[3], 0);
        assert_eq!(classes[4], 1);
    }
    #[test]
    fn test_phase_fraction_unknown_phase_zero() {
        let particles = make_two_phase_particles();
        assert_eq!(phase_fraction(&particles, 99), 0.0);
    }
    #[test]
    fn test_multiphase_fluid_particle_velocity_default_zero() {
        let p = MultiphaseFluidParticle::new([1.0, 2.0, 3.0], 1.0, 1000.0, 0);
        assert_eq!(p.velocity, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_multiphase_params_density_ratio() {
        let mp = MultiphaseParams::new(0.07, 0.005, 800.0);
        assert!((mp.density_ratio - 800.0).abs() < 1e-12);
    }
    #[test]
    fn test_interface_curvature_one_axis() {
        let n = [1.0, 0.0, 0.0];
        let gn = [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let kappa = interface_curvature(n, gn);
        assert!((kappa + 2.0).abs() < 1e-12);
    }
}
/// Cubic spline kernel gradient magnitude |dW/dr| for 3-D.
pub(super) fn cubic_spline_kernel_grad(r: f64, h: f64) -> f64 {
    let q = r / h;
    let norm = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        norm * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        norm * (-0.75 * (2.0 - q).powi(2)) / h
    } else {
        0.0
    }
}
#[cfg(test)]
mod new_struct_tests {

    use crate::multiphase_sph::types::*;
    #[test]
    fn msp_phase0_color_one() {
        let p = MultiphaseSphParticle::new([0.0; 3], 0, 1e-3, 1000.0);
        assert!((p.color_function - 1.0).abs() < 1e-12);
    }
    #[test]
    fn msp_phase1_color_zero() {
        let p = MultiphaseSphParticle::new([0.0; 3], 1, 1e-3, 1.225);
        assert!(p.color_function.abs() < 1e-12);
    }
    #[test]
    fn msp_at_rest_zero_ke() {
        let p = MultiphaseSphParticle::new([1.0, 2.0, 3.0], 0, 0.5, 1000.0);
        assert_eq!(p.kinetic_energy(), 0.0);
    }
    #[test]
    fn msp_is_interface() {
        let mut p = MultiphaseSphParticle::new([0.0; 3], 0, 1e-3, 1000.0);
        p.color_function = 0.5;
        assert!(p.is_interface());
    }
    #[test]
    fn msp_not_interface_bulk_liquid() {
        let p = MultiphaseSphParticle::new([0.0; 3], 0, 1e-3, 1000.0);
        assert!(!p.is_interface());
    }
    #[test]
    fn msp_interpolate_density() {
        let p = MultiphaseSphParticle::new([0.0; 3], 0, 1e-3, 1000.0);
        let rho = p.interpolate(1000.0, 1.225);
        assert!((rho - 1000.0).abs() < 1e-9);
    }
    #[test]
    fn msp_dist_pythagorean() {
        let p1 = MultiphaseSphParticle::new([0.0; 3], 0, 1e-3, 1000.0);
        let p2 = MultiphaseSphParticle::new([3.0, 4.0, 0.0], 1, 1e-3, 1.225);
        assert!((p1.dist_to(&p2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn csf_smooth_color_uniform_phase0() {
        let particles: Vec<MultiphaseSphParticle> = (0..4)
            .map(|i| MultiphaseSphParticle::new([i as f64 * 0.05, 0.0, 0.0], 0, 1e-3, 1000.0))
            .collect();
        let csf = ColorFunctionSph::new(0.072, 0.3);
        let c = csf.smooth_color(&particles, 1);
        assert!(
            (c - 1.0).abs() < 0.1,
            "uniform phase-0 should give C≈1: {c}"
        );
    }
    #[test]
    fn csf_laplace_pressure_sphere() {
        let csf = ColorFunctionSph::new(0.072, 0.01);
        let dp = csf.laplace_pressure(1e-3);
        let expected = 2.0 * 0.072 / 1e-3;
        assert!((dp - expected).abs() < 1.0, "Laplace: {dp} vs {expected}");
    }
    #[test]
    fn csf_color_gradient_finite() {
        let mut particles: Vec<MultiphaseSphParticle> = (0..6)
            .map(|i| {
                let phase = if i < 3 { 0 } else { 1 };
                MultiphaseSphParticle::new(
                    [i as f64 * 0.05, 0.0, 0.0],
                    phase,
                    1e-3,
                    if phase == 0 { 1000.0 } else { 1.225 },
                )
            })
            .collect();
        for p in particles.iter_mut() {
            p.color_function = if p.phase_label == 0 { 1.0 } else { 0.0 };
        }
        let csf = ColorFunctionSph::new(0.072, 0.2);
        let grad = csf.color_gradient(&particles, 2);
        assert!(grad[0].is_finite() && grad[1].is_finite() && grad[2].is_finite());
    }
    #[test]
    fn mp_density_ratio_water_air() {
        let mpd = MultiphaseDensity::new(0.01, 1000.0, 1.225);
        let ratio = mpd.density_ratio();
        assert!((ratio - 1000.0 / 1.225).abs() < 1.0, "ratio: {ratio}");
    }
    #[test]
    fn mp_interpolated_density_pure_phases() {
        let mpd = MultiphaseDensity::new(0.01, 1000.0, 1.225);
        assert!((mpd.interpolated_density(1.0) - 1000.0).abs() < 1e-10);
        assert!((mpd.interpolated_density(0.0) - 1.225).abs() < 1e-10);
    }
    #[test]
    fn mp_interpolated_density_interface() {
        let mpd = MultiphaseDensity::new(0.01, 1000.0, 1000.0);
        assert!((mpd.interpolated_density(0.5) - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn mp_compute_density_positive() {
        let particles: Vec<MultiphaseSphParticle> = (0..5)
            .map(|i| MultiphaseSphParticle::new([i as f64 * 0.01, 0.0, 0.0], 0, 1e-3, 1000.0))
            .collect();
        let mpd = MultiphaseDensity::new(0.1, 1000.0, 1.225);
        let rho = mpd.compute_density(&particles, 2);
        assert!(rho > 0.0, "SPH density should be positive: {rho}");
    }
    #[test]
    fn mp_update_all_no_nan() {
        let mut particles: Vec<MultiphaseSphParticle> = (0..5)
            .map(|i| MultiphaseSphParticle::new([i as f64 * 0.02, 0.0, 0.0], 0, 1e-3, 1000.0))
            .collect();
        let mpd = MultiphaseDensity::new(0.2, 1000.0, 1.225);
        mpd.update_all(&mut particles);
        for p in &particles {
            assert!(p.density.is_finite() && p.density > 0.0, "density NaN/zero");
        }
    }
    #[test]
    fn meos_at_reference_density_zero_pressure() {
        let eos = MultiphaseEquationOfState::water_air();
        let p = eos.pressure(1000.0, 0);
        assert!(p.abs() < 1e-6, "p at ρ₀ should be ~0: {p}");
    }
    #[test]
    fn meos_pressure_increases_with_density() {
        let eos = MultiphaseEquationOfState::water_air();
        let p1 = eos.pressure(1000.0, 0);
        let p2 = eos.pressure(1100.0, 0);
        assert!(p2 > p1);
    }
    #[test]
    fn meos_sound_speed_positive() {
        let eos = MultiphaseEquationOfState::water_air();
        assert!(eos.sound_speed(0) > 0.0);
        assert!(eos.sound_speed(1) > 0.0);
    }
    #[test]
    fn meos_density_ratio_water_air() {
        let eos = MultiphaseEquationOfState::water_air();
        let ratio = eos.density_ratio();
        assert!((ratio - 1000.0 / 1.225).abs() < 1.0);
    }
    #[test]
    fn meos_cfl_timestep_positive() {
        let eos = MultiphaseEquationOfState::water_air();
        let dt = eos.cfl_timestep(0.01, 0.3);
        assert!(dt > 0.0 && dt.is_finite());
    }
    #[test]
    fn meos_update_pressures_no_nan() {
        let eos = MultiphaseEquationOfState::water_air();
        let mut particles: Vec<MultiphaseSphParticle> = (0..4)
            .map(|i| {
                let phase = if i < 2 { 0 } else { 1 };
                let rho = if phase == 0 { 1000.0 } else { 1.225 };
                MultiphaseSphParticle::new([i as f64 * 0.1, 0.0, 0.0], phase, 1e-3, rho)
            })
            .collect();
        eos.update_pressures(&mut particles);
        for p in &particles {
            assert!(p.pressure.is_finite(), "pressure NaN detected");
        }
    }
    #[test]
    fn bdm_initial_radius_stored() {
        let b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        assert!((b.radius - 1e-3).abs() < 1e-14);
    }
    #[test]
    fn bdm_gas_pressure_at_r0() {
        let b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        let p_b = b.gas_pressure();
        assert!(
            (p_b - b.p_b0).abs() < 1.0,
            "gas pressure at R₀: {p_b} vs {}",
            b.p_b0
        );
    }
    #[test]
    fn bdm_minnaert_frequency_positive() {
        let b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        let omega = b.minnaert_frequency();
        assert!(omega > 0.0, "Minnaert frequency: {omega}");
    }
    #[test]
    fn bdm_collapse_time_positive() {
        let b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        let tc = b.collapse_time();
        assert!(tc > 0.0 && tc < 1.0, "collapse time: {tc}");
    }
    #[test]
    fn bdm_step_no_nan() {
        let mut b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        b.step(1e-9, 1e5);
        assert!(!b.radius.is_nan() && !b.radius_dot.is_nan());
    }
    #[test]
    fn bdm_radiated_pressure_finite() {
        let b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        let p = b.radiated_pressure_at([0.01, 0.0, 0.0]);
        assert!(p.is_finite());
    }
    #[test]
    fn bdm_bubble_expands_under_negative_pressure() {
        let mut b = BubbleDynamicsModel::new(1e-3, 1e5, 1000.0, 1e-3, 0.072, [0.0; 3]);
        let p_ambient = 0.5e5;
        b.step(1e-8, p_ambient);
        assert!(b.radius.is_finite() && b.radius > 0.0);
    }
    #[test]
    fn phase_trans_evap_flux_positive_below_psat() {
        let pt = PhaseTransitionSph::water_steam();
        let flux = pt.evaporation_flux(373.15, 50000.0);
        assert!(flux > 0.0, "evaporation flux should be positive: {flux}");
    }
    #[test]
    fn phase_trans_evap_flux_zero_above_psat() {
        let pt = PhaseTransitionSph::water_steam();
        let flux = pt.evaporation_flux(373.15, 200000.0);
        assert_eq!(flux, 0.0);
    }
    #[test]
    fn phase_trans_cond_flux_positive_above_psat() {
        let pt = PhaseTransitionSph::water_steam();
        let flux = pt.condensation_flux(373.15, 200000.0);
        assert!(flux > 0.0, "condensation flux should be positive: {flux}");
    }
    #[test]
    fn phase_trans_net_flux_sign() {
        let pt = PhaseTransitionSph::water_steam();
        let evap = pt.net_mass_flux(373.15, 50000.0);
        let cond = pt.net_mass_flux(373.15, 200000.0);
        assert!(evap > 0.0, "net flux below psat should be positive (evap)");
        assert!(cond <= 0.0, "net flux above psat should be ≤0 (cond)");
    }
    #[test]
    fn phase_trans_clausius_clapeyron_at_tsat() {
        let pt = PhaseTransitionSph::water_steam();
        let p = pt.clausius_clapeyron_pressure(pt.t_sat);
        assert!(
            (p - pt.p_sat).abs() / pt.p_sat < 1e-10,
            "p at T_sat should equal p_sat: {p}"
        );
    }
    #[test]
    fn phase_trans_clausius_clapeyron_increases_with_temp() {
        let pt = PhaseTransitionSph::water_steam();
        let p1 = pt.clausius_clapeyron_pressure(360.0);
        let p2 = pt.clausius_clapeyron_pressure(400.0);
        assert!(p2 > p1, "p_sat should increase with temperature");
    }
    #[test]
    fn phase_trans_energy_sink_negative() {
        let pt = PhaseTransitionSph::water_steam();
        let q = pt.energy_sink(1.0);
        assert!(q < 0.0, "evaporation should be an energy sink: {q}");
    }
    #[test]
    fn phase_trans_jakob_number_positive() {
        let pt = PhaseTransitionSph::water_steam();
        let ja = pt.jakob_number(1000.0, 4186.0, 10.0, 0.6);
        assert!(ja > 0.0, "Jakob number should be positive: {ja}");
    }
    #[test]
    fn phase_trans_custom_params_stored() {
        let pt = PhaseTransitionSph::new(1e6, 0.2, 0.3, 0.04, 300.0, 5000.0);
        assert!((pt.latent_heat - 1e6).abs() < 1.0);
        assert!((pt.t_sat - 300.0).abs() < 1e-10);
    }
}
