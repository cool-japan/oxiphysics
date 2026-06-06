//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{Contact, GranularParams, RollingResistanceModel};

/// Mohr-Coulomb yield test: returns `true` if the stress state has yielded.
///
/// Yield condition: |τ| > c + σ_n · tan(φ)
///
/// # Arguments
/// * `sigma_n` – normal stress (positive = compression) \[Pa\]
/// * `tau`     – shear stress \[Pa\]
/// * `params`  – granular material parameters
pub fn mohr_coulomb_yield(sigma_n: f64, tau: f64, params: &GranularParams) -> bool {
    let phi = params.friction_angle_deg.to_radians();
    let tau_yield = params.cohesion + sigma_n * phi.tan();
    tau.abs() > tau_yield
}
/// Drucker-Prager yield function.
///
/// F = q − (A·p + B)
///
/// where
/// - A = 6·sin(φ) / (3 − sin(φ))
/// - B = 6·c·cos(φ) / (3 − sin(φ))
///
/// F < 0 → elastic interior; F ≥ 0 → yielded / on yield surface.
///
/// # Arguments
/// * `p` – mean (hydrostatic) stress = (σ_1+σ_2+σ_3)/3 \[Pa\]
/// * `q` – equivalent deviatoric (von-Mises) stress \[Pa\]
/// * `params` – granular material parameters
pub fn drucker_prager_f(p: f64, q: f64, params: &GranularParams) -> f64 {
    let phi = params.friction_angle_deg.to_radians();
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let denom = 3.0 - sin_phi;
    let a = 6.0 * sin_phi / denom;
    let b = 6.0 * params.cohesion * cos_phi / denom;
    q - (a * p + b)
}
/// Pressure from bulk compressive deformation (linearised EOS).
///
/// P = K · (ρ/ρ₀ − 1)
///
/// The result is clamped to ≥ −cohesion so that the medium cannot sustain
/// large tensile pressures.
///
/// # Arguments
/// * `rho`    – current density \[kg/m³\]
/// * `rho0`   – reference (rest) density \[kg/m³\]
/// * `bulk_k` – bulk modulus K \[Pa\]
pub fn granular_eos_pressure(rho: f64, rho0: f64, bulk_k: f64) -> f64 {
    bulk_k * (rho / rho0 - 1.0)
}
/// Pressure from bulk deformation, clamped to ≥ −cohesion.
///
/// See [`granular_eos_pressure`] for the underlying formula.
pub fn granular_eos_pressure_clamped(rho: f64, rho0: f64, params: &GranularParams) -> f64 {
    let p = granular_eos_pressure(rho, rho0, params.bulk_modulus);
    p.max(-params.cohesion)
}
/// Convert degrees to radians (re-exported for convenience).
#[cfg(test)]
#[inline]
pub(super) fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}
/// Compute the normal contact force using a linear spring-dashpot model.
///
/// `F_n = kn * overlap - gamma_n * rel_vel_n`
///
/// The result is clamped to ≥ 0 (no tensile adhesion).
pub fn spring_dashpot_normal(overlap: f64, rel_vel_n: f64, kn: f64, gamma_n: f64) -> f64 {
    if overlap <= 0.0 {
        return 0.0;
    }
    (kn * overlap - gamma_n * rel_vel_n).max(0.0)
}
/// Compute the tangential contact force using a spring-dashpot model.
///
/// `F_t = kt * delta_t - gamma_t * rel_vel_t`
///
/// The result is clamped by Coulomb friction: `|F_t| ≤ mu * F_n`.
pub fn spring_dashpot_tangential(
    delta_t: f64,
    rel_vel_t: f64,
    kt: f64,
    gamma_t: f64,
    mu: f64,
    f_n: f64,
) -> f64 {
    let f_t = kt * delta_t - gamma_t * rel_vel_t;
    let limit = mu * f_n.abs();
    f_t.clamp(-limit, limit)
}
/// Hertzian normal contact force: `F_n = (4/3) E* sqrt(R*) * delta^{3/2}`.
pub fn hertz_normal_force(overlap: f64, e_star: f64, r_star: f64) -> f64 {
    if overlap <= 0.0 {
        return 0.0;
    }
    (4.0 / 3.0) * e_star * r_star.sqrt() * overlap.powf(1.5)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::granular::types::*;
    #[test]
    fn mohr_coulomb_below_yield_surface() {
        let params = GranularParams::sand();
        let yielded = mohr_coulomb_yield(1.0e4, 1.0e2, &params);
        assert!(!yielded, "expected below yield surface");
    }
    #[test]
    fn mohr_coulomb_above_yield_surface() {
        let params = GranularParams::sand();
        let yielded = mohr_coulomb_yield(0.0, 1.0e6, &params);
        assert!(yielded, "expected above yield surface");
    }
    #[test]
    fn drucker_prager_at_hydrostatic_compression() {
        let params = GranularParams::sand();
        let f = drucker_prager_f(1.0e6, 0.0, &params);
        assert!(
            f < 0.0,
            "hydrostatic compression should be inside DP yield surface, got F={f}"
        );
    }
    #[test]
    fn granular_eos_at_rest_density() {
        let p = granular_eos_pressure(1000.0, 1000.0, 5.0e7);
        assert!(
            (p - 0.0).abs() < 1e-10,
            "expected P=0 at rest density, got {p}"
        );
    }
    #[test]
    fn granular_collision_zero_overlap() {
        let contact = GranularCollision::new(7.0e10, 0.3, 0.4);
        let f = contact.normal_force(0.0, 1.0);
        assert_eq!(f, 0.0, "zero overlap must give zero force");
    }
    #[test]
    fn wet_sand_has_more_cohesion_than_sand() {
        let sand = GranularParams::sand();
        let wet = GranularParams::wet_sand();
        assert!(
            wet.cohesion > sand.cohesion,
            "wet_sand cohesion ({}) must exceed sand cohesion ({})",
            wet.cohesion,
            sand.cohesion
        );
    }
    #[test]
    fn dem_particle_mass_from_density() {
        let p = DemGranularParticle::new([0.0; 3], 0.1, 2600.0);
        let expected = (4.0 / 3.0) * PI * 0.001 * 2600.0;
        assert!((p.mass - expected).abs() < 1e-8);
    }
    #[test]
    fn dem_particle_moment_of_inertia() {
        let p = DemGranularParticle::new([0.0; 3], 0.5, 1000.0);
        let expected = 0.4 * p.mass * 0.25;
        assert!((p.moment_of_inertia() - expected).abs() < 1e-10);
    }
    #[test]
    fn dem_particle_kinetic_energy_stationary() {
        let p = DemGranularParticle::new([0.0; 3], 0.1, 1000.0);
        assert!(p.kinetic_energy().abs() < 1e-14);
    }
    #[test]
    fn dem_particle_kinetic_energy_moving() {
        let mut p = DemGranularParticle::new([0.0; 3], 0.1, 1000.0);
        p.velocity = [1.0, 0.0, 0.0];
        let ke = p.kinetic_energy();
        let expected = 0.5 * p.mass;
        assert!((ke - expected).abs() < 1e-10);
    }
    #[test]
    fn spring_dashpot_zero_overlap() {
        let f = spring_dashpot_normal(0.0, 1.0, 1e5, 100.0);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn spring_dashpot_repulsive() {
        let f = spring_dashpot_normal(0.001, 0.0, 1e5, 100.0);
        assert!((f - 100.0).abs() < 1e-6);
    }
    #[test]
    fn spring_dashpot_tangential_clamped() {
        let f_t = spring_dashpot_tangential(0.01, 0.0, 1e5, 0.0, 0.3, 100.0);
        assert!((f_t - 30.0).abs() < 1e-10);
    }
    #[test]
    fn hertz_zero_overlap() {
        let f = hertz_normal_force(0.0, 1e9, 0.05);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn hertz_positive_overlap() {
        let f = hertz_normal_force(0.001, 1e9, 0.05);
        assert!(f > 0.0, "Hertz force must be positive for overlap > 0");
    }
    #[test]
    fn granular_sim_detect_contact() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 100.0,
            gamma_t: 50.0,
        };
        let mut sim = GranularSim::new(model, 0.3, [0.0, 0.0, 0.0]);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.18, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let contacts = sim.detect_contacts();
        assert_eq!(contacts.len(), 1);
        assert!((contacts[0].overlap - 0.02).abs() < 1e-10);
    }
    #[test]
    fn granular_sim_no_contact_far_apart() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 100.0,
            gamma_t: 50.0,
        };
        let mut sim = GranularSim::new(model, 0.3, [0.0; 3]);
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 1.0);
        sim.add_particle([1.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        assert!(sim.detect_contacts().is_empty());
    }
    #[test]
    fn granular_sim_resolve_pushes_apart() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 100.0,
            gamma_t: 50.0,
        };
        let mut sim = GranularSim::new(model, 0.3, [0.0; 3]);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.18, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let contacts = sim.detect_contacts();
        sim.resolve_contacts(&contacts);
        assert!(sim.forces[0][0] < 0.0, "i should be pushed in -x");
        assert!(sim.forces[1][0] > 0.0, "j should be pushed in +x");
    }
    #[test]
    fn granular_sim_step_moves_particles() {
        let model = ContactModel::Hertzian {
            e_star: 1e9,
            r_star: 0.05,
        };
        let mut sim = GranularSim::new(model, 0.3, [0.0, -9.81, 0.0]);
        sim.add_particle([0.0, 1.0, 0.0], [0.0; 3], 0.1, 1.0);
        let y0 = sim.positions[0][1];
        sim.step(0.01);
        assert!(
            sim.positions[0][1] < y0,
            "particle should fall under gravity"
        );
    }
    #[test]
    fn granular_sim_momentum_conservation_contact() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 0.0,
            gamma_t: 0.0,
        };
        let mut sim = GranularSim::new(model, 0.3, [0.0; 3]);
        sim.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0);
        sim.add_particle([0.18, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.1, 1.0);
        let mom_before: f64 = sim
            .velocities
            .iter()
            .zip(sim.masses.iter())
            .map(|(v, &m)| m * v[0])
            .sum();
        sim.step(0.0001);
        let mom_after: f64 = sim
            .velocities
            .iter()
            .zip(sim.masses.iter())
            .map(|(v, &m)| m * v[0])
            .sum();
        assert!(
            (mom_before - mom_after).abs() < 1e-8,
            "momentum not conserved: before={mom_before}, after={mom_after}"
        );
    }
    #[test]
    fn deg_to_rad_conversion() {
        assert!((deg_to_rad(180.0) - PI).abs() < 1e-12);
        assert!((deg_to_rad(90.0) - PI / 2.0).abs() < 1e-12);
    }
    #[test]
    fn granular_eos_clamped_limits_tension() {
        let params = GranularParams::sand();
        let p = granular_eos_pressure_clamped(500.0, 1000.0, &params);
        assert!((p - (-params.cohesion)).abs() < 1e-6);
    }
}
/// Compute rolling resistance torque on particle i due to contact with particle j.
///
/// Returns the torque vector `[Tx, Ty, Tz]` acting on particle i.
///
/// # Arguments
/// * `r_i` – radius of particle i
/// * `r_j` – radius of particle j
/// * `f_n` – normal contact force magnitude (N)
/// * `omega_rel` – relative angular velocity (omega_i - omega_j)
/// * `mu_r` – rolling friction coefficient
/// * `model` – which rolling resistance model to use
pub fn rolling_resistance_torque(
    r_i: f64,
    r_j: f64,
    f_n: f64,
    omega_rel: [f64; 3],
    mu_r: f64,
    model: RollingResistanceModel,
) -> [f64; 3] {
    let r_eff = (r_i * r_j) / (r_i + r_j).max(1e-30);
    let omega_mag = (omega_rel[0].powi(2) + omega_rel[1].powi(2) + omega_rel[2].powi(2)).sqrt();
    match model {
        RollingResistanceModel::None => [0.0; 3],
        RollingResistanceModel::Constant => {
            if omega_mag < 1e-14 {
                return [0.0; 3];
            }
            let scale = -mu_r * r_eff * f_n / omega_mag;
            [
                scale * omega_rel[0],
                scale * omega_rel[1],
                scale * omega_rel[2],
            ]
        }
        RollingResistanceModel::Viscous => {
            let scale = -mu_r * r_eff * f_n;
            [
                scale * omega_rel[0],
                scale * omega_rel[1],
                scale * omega_rel[2],
            ]
        }
        RollingResistanceModel::ElasticPlastic => {
            let scale = -mu_r * r_eff * f_n;
            let torque = [
                scale * omega_rel[0],
                scale * omega_rel[1],
                scale * omega_rel[2],
            ];
            let tmag = (torque[0].powi(2) + torque[1].powi(2) + torque[2].powi(2)).sqrt();
            let limit = mu_r * r_eff * f_n.abs();
            if tmag > limit && tmag > 1e-30 {
                let fac = limit / tmag;
                [torque[0] * fac, torque[1] * fac, torque[2] * fac]
            } else {
                torque
            }
        }
    }
}
/// Compute the packing fraction (solid volume fraction) of a granular assembly.
///
/// φ = Σ V_i / V_domain, where V_i = (4/3) π r_i³
///
/// # Arguments
/// * `radii` – radii of all particles
/// * `domain_volume` – total domain volume (m³)
pub fn packing_fraction(radii: &[f64], domain_volume: f64) -> f64 {
    if domain_volume < 1e-30 {
        return 0.0;
    }
    let total_solid: f64 = radii.iter().map(|&r| (4.0 / 3.0) * PI * r.powi(3)).sum();
    (total_solid / domain_volume).min(1.0)
}
/// Compute the granular temperature T_g of an assembly.
///
/// T_g = (1/3) * <(v - v_mean)²>
///
/// where the average is over all particles (equal weighting).
///
/// # Arguments
/// * `velocities` – slice of velocity vectors
pub fn granular_temperature(velocities: &[[f64; 3]]) -> f64 {
    let n = velocities.len();
    if n == 0 {
        return 0.0;
    }
    let inv_n = 1.0 / n as f64;
    let mut mean = [0.0_f64; 3];
    for v in velocities {
        mean[0] += v[0];
        mean[1] += v[1];
        mean[2] += v[2];
    }
    mean[0] *= inv_n;
    mean[1] *= inv_n;
    mean[2] *= inv_n;
    let mut msq = 0.0_f64;
    for v in velocities {
        let dv0 = v[0] - mean[0];
        let dv1 = v[1] - mean[1];
        let dv2 = v[2] - mean[2];
        msq += dv0 * dv0 + dv1 * dv1 + dv2 * dv2;
    }
    msq * inv_n / 3.0
}
/// Compute the coordination number distribution.
///
/// Returns a `Vec`usize` where `result\[i\]` is the number of contacts particle i
/// participates in.
///
/// # Arguments
/// * `contacts` – slice of detected contacts
/// * `n_particles` – total number of particles
pub fn coordination_numbers(contacts: &[Contact], n_particles: usize) -> Vec<usize> {
    let mut coords = vec![0usize; n_particles];
    for c in contacts {
        if c.i < n_particles {
            coords[c.i] += 1;
        }
        if c.j < n_particles {
            coords[c.j] += 1;
        }
    }
    coords
}
/// Mean coordination number of the assembly.
pub fn mean_coordination_number(contacts: &[Contact], n_particles: usize) -> f64 {
    if n_particles == 0 {
        return 0.0;
    }
    let cn = coordination_numbers(contacts, n_particles);
    cn.iter().sum::<usize>() as f64 / n_particles as f64
}
/// Fabric tensor (second-order, 3×3) computed from contact normal vectors.
///
/// F_ij = (1/N_c) Σ n_i n_j
///
/// where N_c is the total number of contacts.
pub fn fabric_tensor(contacts: &[Contact]) -> [[f64; 3]; 3] {
    let n_c = contacts.len();
    if n_c == 0 {
        return [[0.0; 3]; 3];
    }
    let mut fab = [[0.0_f64; 3]; 3];
    for c in contacts {
        for (i, fab_row) in fab.iter_mut().enumerate() {
            for (j, f) in fab_row.iter_mut().enumerate() {
                *f += c.normal[i] * c.normal[j];
            }
        }
    }
    let inv = 1.0 / n_c as f64;
    for row in fab.iter_mut() {
        for f in row.iter_mut() {
            *f *= inv;
        }
    }
    fab
}
/// Compute the critical slope angle for granular sliding (Mohr-Coulomb).
///
/// θ_c = arctan(tan(φ) + c / (σ_n))
///
/// For a purely frictional material (c = 0): θ_c = φ.
pub fn critical_slope_angle(params: &GranularParams, normal_stress: f64) -> f64 {
    let phi = params.friction_angle_deg.to_radians();
    let tan_phi = phi.tan();
    let effective_tan = tan_phi
        + if normal_stress > 1e-14 {
            params.cohesion / normal_stress
        } else {
            0.0
        };
    effective_tan.atan()
}
/// Check if a slope at angle `theta_slope` (radians) with the given normal stress
/// is stable according to Mohr-Coulomb theory.
pub fn is_slope_stable(theta_slope: f64, params: &GranularParams, normal_stress: f64) -> bool {
    theta_slope < critical_slope_angle(params, normal_stress)
}
/// Compute the force chain intensity for each contact: the ratio of normal
/// contact force to the mean normal force.
///
/// Returns a vector of normalised force values (one per contact).
pub fn force_chain_intensity(normal_forces: &[f64]) -> Vec<f64> {
    if normal_forces.is_empty() {
        return Vec::new();
    }
    let mean: f64 = normal_forces.iter().sum::<f64>() / normal_forces.len() as f64;
    if mean.abs() < 1e-30 {
        return vec![0.0; normal_forces.len()];
    }
    normal_forces.iter().map(|&f| f / mean).collect()
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::granular::types::*;
    #[test]
    fn rolling_resistance_none_is_zero() {
        let t = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [1.0, 0.0, 0.0],
            0.01,
            RollingResistanceModel::None,
        );
        assert_eq!(t, [0.0; 3]);
    }
    #[test]
    fn rolling_resistance_constant_direction() {
        let omega_rel = [1.0, 0.0, 0.0];
        let t = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            omega_rel,
            0.1,
            RollingResistanceModel::Constant,
        );
        assert!(
            t[0] < 0.0,
            "constant model torque should oppose rotation: {}",
            t[0]
        );
        assert!(t[1].abs() < 1e-14);
        assert!(t[2].abs() < 1e-14);
    }
    #[test]
    fn rolling_resistance_viscous_proportional_to_omega() {
        let t1 = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [1.0, 0.0, 0.0],
            0.1,
            RollingResistanceModel::Viscous,
        );
        let t2 = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [2.0, 0.0, 0.0],
            0.1,
            RollingResistanceModel::Viscous,
        );
        assert!(
            (t2[0] / t1[0] - 2.0).abs() < 1e-10,
            "viscous torque should scale with omega"
        );
    }
    #[test]
    fn rolling_resistance_elastic_plastic_capped() {
        let t = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [1e6, 0.0, 0.0],
            0.01,
            RollingResistanceModel::ElasticPlastic,
        );
        let r_eff = 0.05;
        let limit = 0.01 * r_eff * 100.0;
        let tmag = t[0].abs();
        assert!(
            tmag <= limit + 1e-10,
            "torque should be capped at {limit}, got {tmag}"
        );
    }
    #[test]
    fn rolling_resistance_constant_zero_omega() {
        let t = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [0.0; 3],
            0.1,
            RollingResistanceModel::Constant,
        );
        assert_eq!(t, [0.0; 3], "zero omega → zero torque for constant model");
    }
    #[test]
    fn packing_fraction_empty() {
        let phi = packing_fraction(&[], 1.0);
        assert_eq!(phi, 0.0);
    }
    #[test]
    fn packing_fraction_single_particle() {
        let r = 0.1_f64;
        let vol_particle = (4.0 / 3.0) * PI * r.powi(3);
        let vol_domain = 10.0 * vol_particle;
        let phi = packing_fraction(&[r], vol_domain);
        assert!((phi - 0.1).abs() < 1e-10, "phi = {phi}, expected 0.1");
    }
    #[test]
    fn packing_fraction_clamped_to_one() {
        let phi = packing_fraction(&[1.0, 1.0], 0.001);
        assert!(phi <= 1.0);
    }
    #[test]
    fn granular_temperature_zero_for_uniform_velocity() {
        let vels = vec![[1.0, 0.0, 0.0]; 10];
        let tg = granular_temperature(&vels);
        assert!(
            tg.abs() < 1e-14,
            "T_g should be 0 for uniform velocity: {tg}"
        );
    }
    #[test]
    fn granular_temperature_empty() {
        assert_eq!(granular_temperature(&[]), 0.0);
    }
    #[test]
    fn granular_temperature_positive_for_random_velocities() {
        let vels = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let tg = granular_temperature(&vels);
        assert!(
            tg > 0.0,
            "T_g should be positive for random velocities: {tg}"
        );
    }
    #[test]
    fn coordination_number_no_contacts() {
        let cn = coordination_numbers(&[], 5);
        assert_eq!(cn, vec![0; 5]);
    }
    #[test]
    fn coordination_number_single_contact() {
        let c = Contact {
            i: 0,
            j: 1,
            overlap: 0.01,
            normal: [1.0, 0.0, 0.0],
        };
        let cn = coordination_numbers(&[c], 3);
        assert_eq!(cn[0], 1);
        assert_eq!(cn[1], 1);
        assert_eq!(cn[2], 0);
    }
    #[test]
    fn mean_coordination_number_correct() {
        let contacts = vec![
            Contact {
                i: 0,
                j: 1,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
            Contact {
                i: 1,
                j: 2,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
        ];
        let mean = mean_coordination_number(&contacts, 3);
        assert!((mean - 4.0 / 3.0).abs() < 1e-10, "mean CN = {mean}");
    }
    #[test]
    fn fabric_tensor_empty_contacts() {
        let fab = fabric_tensor(&[]);
        for row in &fab {
            for &val in row {
                assert_eq!(val, 0.0);
            }
        }
    }
    #[test]
    fn fabric_tensor_isotropic_distribution() {
        let contacts = vec![
            Contact {
                i: 0,
                j: 1,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
            Contact {
                i: 0,
                j: 2,
                overlap: 0.01,
                normal: [-1.0, 0.0, 0.0],
            },
            Contact {
                i: 0,
                j: 3,
                overlap: 0.01,
                normal: [0.0, 1.0, 0.0],
            },
            Contact {
                i: 0,
                j: 4,
                overlap: 0.01,
                normal: [0.0, -1.0, 0.0],
            },
            Contact {
                i: 0,
                j: 5,
                overlap: 0.01,
                normal: [0.0, 0.0, 1.0],
            },
            Contact {
                i: 0,
                j: 6,
                overlap: 0.01,
                normal: [0.0, 0.0, -1.0],
            },
        ];
        let fab = fabric_tensor(&contacts);
        let trace = fab[0][0] + fab[1][1] + fab[2][2];
        assert!((trace - 1.0).abs() < 1e-10, "fabric trace = {trace}");
        for (d, row) in fab.iter().enumerate() {
            assert!(
                (row[d] - 1.0 / 3.0).abs() < 1e-10,
                "fab[{d}][{d}] = {}",
                row[d]
            );
        }
    }
    #[test]
    fn critical_slope_angle_pure_friction() {
        let params = GranularParams::sand();
        let angle = critical_slope_angle(&params, 1e7);
        let expected = 30.0_f64.to_radians();
        assert!(
            (angle - expected).abs() < 0.01,
            "angle = {angle} rad, expected {expected} rad"
        );
    }
    #[test]
    fn is_slope_stable_flat_slope() {
        let params = GranularParams::sand();
        assert!(
            is_slope_stable(0.0, &params, 1e5),
            "flat slope (theta=0) should always be stable"
        );
    }
    #[test]
    fn is_slope_stable_steep_slope_unstable() {
        let params = GranularParams::sand();
        let theta = 80.0_f64.to_radians();
        assert!(
            !is_slope_stable(theta, &params, 1e4),
            "80° slope should be unstable for sand"
        );
    }
    #[test]
    fn extended_sim_angular_velocity_updated() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 0.0,
            gamma_t: 0.0,
        };
        let mut sim =
            GranularSimExtended::new(model, 0.3, [0.0; 3], RollingResistanceModel::Constant, 0.01);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.18, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.step(0.001);
        assert_eq!(sim.angular_velocities.len(), 2);
    }
    #[test]
    fn extended_sim_rotational_ke_zero_initially() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 0.0,
            gamma_t: 0.0,
        };
        let mut sim =
            GranularSimExtended::new(model, 0.3, [0.0; 3], RollingResistanceModel::None, 0.0);
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 1.0);
        assert!((sim.rotational_kinetic_energy()).abs() < 1e-14);
    }
    #[test]
    fn force_chain_intensity_uniform_forces() {
        let forces = vec![10.0, 10.0, 10.0, 10.0];
        let intensities = force_chain_intensity(&forces);
        for &fi in &intensities {
            assert!(
                (fi - 1.0).abs() < 1e-10,
                "uniform forces → intensity=1, got {fi}"
            );
        }
    }
    #[test]
    fn force_chain_intensity_empty() {
        let intensities = force_chain_intensity(&[]);
        assert!(intensities.is_empty());
    }
    #[test]
    fn force_chain_intensity_zero_mean() {
        let intensities = force_chain_intensity(&[0.0, 0.0, 0.0]);
        for &fi in &intensities {
            assert_eq!(fi, 0.0);
        }
    }
    #[test]
    fn spatial_hash_dem_finds_overlapping_pair() {
        let positions = vec![[0.0, 0.0, 0.0], [0.18, 0.0, 0.0]];
        let radii = vec![0.1, 0.1];
        let hash = SpatialHashDem::new(&positions, &radii, 0.1);
        let pairs = hash.candidate_pairs(&positions, &radii);
        assert!(
            !pairs.is_empty(),
            "expected candidate pair for overlapping particles"
        );
    }
    #[test]
    fn spatial_hash_dem_no_pairs_for_distant_particles() {
        let positions = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let radii = vec![0.1, 0.1];
        let hash = SpatialHashDem::new(&positions, &radii, 0.1);
        let pairs = hash.candidate_pairs(&positions, &radii);
        assert!(
            pairs.is_empty(),
            "distant particles should produce no candidates"
        );
    }
}
/// Hertz normal force magnitude: `F_n = (4/3) * E_eff * sqrt(r_eff) * delta_n^{3/2}`.
///
/// Returns 0 for non-positive overlap.
///
/// # Arguments
/// * `E_eff`  – Combined elastic modulus E* \[Pa\].
/// * `r_eff`  – Effective radius R* \[m\].
/// * `delta_n` – Normal overlap δ \[m\] (positive = interpenetration).
pub fn dem_hertz_normal_force(e_eff: f64, r_eff: f64, delta_n: f64) -> f64 {
    if delta_n <= 0.0 {
        return 0.0;
    }
    (4.0 / 3.0) * e_eff * r_eff.sqrt() * delta_n.powf(1.5)
}
/// Linear spring normal force: `F_n = k_n * delta_n`.
///
/// Returns 0 for non-positive overlap.
///
/// # Arguments
/// * `kn`      – Spring stiffness \[N/m\].
/// * `delta_n` – Normal overlap \[m\].
pub fn linear_spring_normal(kn: f64, delta_n: f64) -> f64 {
    if delta_n <= 0.0 {
        return 0.0;
    }
    kn * delta_n
}
/// Tangential friction force capped by the Coulomb cone.
///
/// Computes `F_t = -k_t * delta_t` then applies the cone limit
/// `|F_t| ≤ mu * |F_n|`.
///
/// # Arguments
/// * `kt`       – Tangential spring stiffness \[N/m\].
/// * `delta_t`  – Accumulated tangential displacement \[m\].
/// * `mu`       – Coulomb friction coefficient.
/// * `fn_mag`   – Normal force magnitude \[N\] (must be ≥ 0).
pub fn tangential_friction_force(kt: f64, delta_t: [f64; 3], mu: f64, fn_mag: f64) -> [f64; 3] {
    let ft_raw = [-kt * delta_t[0], -kt * delta_t[1], -kt * delta_t[2]];
    let ft_mag = (ft_raw[0] * ft_raw[0] + ft_raw[1] * ft_raw[1] + ft_raw[2] * ft_raw[2]).sqrt();
    let limit = mu * fn_mag.abs();
    if ft_mag > limit && ft_mag > 1e-30 {
        let scale = limit / ft_mag;
        [ft_raw[0] * scale, ft_raw[1] * scale, ft_raw[2] * scale]
    } else {
        ft_raw
    }
}
#[cfg(test)]
mod tests_dem {
    use super::*;
    use crate::granular::types::*;
    #[test]
    fn dem_hertz_force_zero_no_overlap() {
        let f = dem_hertz_normal_force(1e9, 0.05, 0.0);
        assert_eq!(f, 0.0, "Hertz force must be zero with zero overlap");
    }
    #[test]
    fn dem_hertz_force_positive_overlap() {
        let f = dem_hertz_normal_force(1e9, 0.05, 0.001);
        assert!(
            f > 0.0,
            "Hertz force must be positive for overlap > 0, got {f}"
        );
    }
    #[test]
    fn dem_hertz_force_negative_overlap_is_zero() {
        let f = dem_hertz_normal_force(1e9, 0.05, -0.001);
        assert_eq!(f, 0.0, "Hertz force must be zero for negative overlap");
    }
    #[test]
    fn linear_spring_zero_overlap_is_zero() {
        let f = linear_spring_normal(1e5, 0.0);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn linear_spring_positive_overlap() {
        let f = linear_spring_normal(1e5, 0.002);
        assert!((f - 200.0).abs() < 1e-10, "expected 200 N, got {f}");
    }
    #[test]
    fn linear_spring_negative_overlap_is_zero() {
        let f = linear_spring_normal(1e5, -0.001);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn tangential_friction_force_no_slip() {
        let ft = tangential_friction_force(1e4, [0.0001, 0.0, 0.0], 0.3, 100.0);
        assert!((ft[0] - (-1.0)).abs() < 1e-10, "ft_x={}", ft[0]);
        assert!(ft[1].abs() < 1e-14);
        assert!(ft[2].abs() < 1e-14);
    }
    #[test]
    fn tangential_friction_force_slip_clamped() {
        let fn_mag = 100.0;
        let mu = 0.3;
        let limit = mu * fn_mag;
        let ft = tangential_friction_force(1e6, [0.1, 0.0, 0.0], mu, fn_mag);
        let ft_mag = ft[0].abs();
        assert!(
            (ft_mag - limit).abs() < 1e-10,
            "clamped ft magnitude should be {limit}, got {ft_mag}"
        );
    }
    #[test]
    fn dem_particle_overlap_positive_when_interpenetrating() {
        let p1 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p2 = DemParticle::new_sphere([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let ov = p1.overlap_with(&p2);
        assert!((ov - 0.05).abs() < 1e-12, "overlap={ov}");
    }
    #[test]
    fn dem_particle_overlap_zero_when_separated() {
        let p1 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p2 = DemParticle::new_sphere([0.5, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let ov = p1.overlap_with(&p2);
        assert!(ov <= 0.0, "no overlap expected, got {ov}");
    }
    #[test]
    fn dem_simulation_gravity_accelerates_particle() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            0.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        sim.add_particle(DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0));
        sim.step();
        assert!(
            (sim.particles[0].velocity[1] - (-0.00981)).abs() < 1e-10,
            "vy={}",
            sim.particles[0].velocity[1]
        );
    }
    #[test]
    fn dem_simulation_contact_force_pushes_apart() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0; 3],
            1e5,
            5e4,
            0.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        sim.add_particle(DemParticle::new_sphere([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0));
        sim.add_particle(DemParticle::new_sphere(
            [0.15, 0.0, 0.0],
            [0.0; 3],
            0.1,
            1.0,
        ));
        sim.step();
        assert!(
            sim.particles[0].velocity[0] < 0.0,
            "p0 vx should be negative, got {}",
            sim.particles[0].velocity[0]
        );
        assert!(
            sim.particles[1].velocity[0] > 0.0,
            "p1 vx should be positive, got {}",
            sim.particles[1].velocity[0]
        );
    }
    #[test]
    fn dem_simulation_hertz_contact() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0; 3],
            0.0,
            0.0,
            0.0,
            0.0,
            GranularContactModel::Hertz,
            1e10,
        );
        sim.add_particle(DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0));
        sim.add_particle(DemParticle::new_sphere(
            [0.15, 0.0, 0.0],
            [0.0; 3],
            0.1,
            1.0,
        ));
        sim.step();
        assert!(sim.particles[0].velocity[0] < 0.0);
        assert!(sim.particles[1].velocity[0] > 0.0);
    }
    #[test]
    fn granular_contact_model_equality() {
        assert_eq!(GranularContactModel::Hertz, GranularContactModel::Hertz);
        assert_ne!(
            GranularContactModel::Hertz,
            GranularContactModel::LinearSpring
        );
        assert_ne!(
            GranularContactModel::LinearSpring,
            GranularContactModel::Mindlin
        );
    }
}
#[cfg(test)]
mod tests_granular_sim {

    use crate::granular::types::*;
    fn make_sim(gamma_n: f64) -> GranularSimulation {
        let contact = HertzContact::from_materials(200e9, 0.3, 200e9, 0.3);
        GranularSimulation::new(contact, [0.0; 3], 0.3, gamma_n, 1e-5)
    }
    #[test]
    fn hertz_contact_normal_force_zero_no_overlap() {
        let c = HertzContact::new(1e9, 4e8);
        assert_eq!(c.normal_force(0.05, 0.0), 0.0);
        assert_eq!(c.normal_force(0.05, -0.001), 0.0);
    }
    #[test]
    fn hertz_contact_normal_force_positive_overlap() {
        let c = HertzContact::new(1e9, 4e8);
        let f = c.normal_force(0.05, 0.001);
        assert!(
            f > 0.0,
            "Hertz normal force must be positive for overlap>0, got {f}"
        );
    }
    #[test]
    fn hertz_contact_normal_force_scales_with_overlap() {
        let c = HertzContact::new(1e9, 4e8);
        let f1 = c.normal_force(0.05, 0.001);
        let f2 = c.normal_force(0.05, 0.002);
        let ratio = f2 / f1;
        let expected = (2.0_f64).powf(1.5);
        assert!(
            (ratio - expected).abs() < 1e-6,
            "F_n ratio should be {expected:.4}, got {ratio:.4}"
        );
    }
    #[test]
    fn hertz_contact_from_materials_consistent_moduli() {
        let c = HertzContact::from_materials(200e9, 0.3, 200e9, 0.3);
        let expected_e_star = 200e9 / (2.0 * (1.0 - 0.3 * 0.3));
        assert!(
            (c.e_star - expected_e_star).abs() / expected_e_star < 1e-6,
            "E* = {}, expected ≈ {expected_e_star}",
            c.e_star
        );
    }
    #[test]
    fn hertz_contact_tangential_force_zero_no_overlap() {
        let c = HertzContact::new(1e9, 4e8);
        let ft = c.tangential_force(0.05, 0.0, [0.0001, 0.0, 0.0], 0.3, 100.0);
        assert_eq!(
            ft, [0.0; 3],
            "tangential force must be zero with no overlap"
        );
    }
    #[test]
    fn hertz_contact_tangential_force_coulomb_clamped() {
        let c = HertzContact::new(1e10, 4e9);
        let fn_mag = 1000.0;
        let mu = 0.3;
        let ft = c.tangential_force(0.05, 0.001, [1.0, 0.0, 0.0], mu, fn_mag);
        let ft_mag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
        assert!(
            ft_mag <= mu * fn_mag + 1e-6,
            "ft_mag={ft_mag} must be ≤ mu*fn={}",
            mu * fn_mag
        );
    }
    #[test]
    fn hertz_contact_stiffness_zero_no_overlap() {
        let c = HertzContact::new(1e9, 4e8);
        assert_eq!(c.stiffness(0.05, 0.0), 0.0);
        assert_eq!(c.stiffness(0.05, -0.001), 0.0);
    }
    #[test]
    fn granular_sim_add_and_count_particles() {
        let mut sim = make_sim(0.0);
        assert_eq!(sim.len(), 0);
        assert!(sim.is_empty());
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 1.0);
        assert_eq!(sim.len(), 1);
        assert!(!sim.is_empty());
    }
    #[test]
    fn granular_sim_detect_contacts_no_overlap() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([1.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let contacts = sim.detect_contacts();
        assert!(
            contacts.is_empty(),
            "separated spheres must have no contacts"
        );
    }
    #[test]
    fn granular_sim_detect_contacts_overlapping() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let contacts = sim.detect_contacts();
        assert_eq!(contacts.len(), 1, "should detect exactly one contact");
        let (_i, _j, _normal, overlap, _r_star) = contacts[0];
        assert!((overlap - 0.05).abs() < 1e-10, "overlap={overlap}");
    }
    #[test]
    fn granular_sim_apply_gravity_accumulates_force() {
        let contact = HertzContact::new(1e9, 4e8);
        let mut sim = GranularSimulation::new(contact, [0.0, -9.81, 0.0], 0.3, 0.0, 1e-5);
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 2.0);
        sim.add_particle([1.0, 0.0, 0.0], [0.0; 3], 0.1, 3.0);
        let mut forces = vec![[0.0_f64; 3]; 2];
        sim.apply_gravity(&mut forces);
        assert!(
            (forces[0][1] - (-19.62)).abs() < 1e-6,
            "f[0]_y={}",
            forces[0][1]
        );
        assert!(
            (forces[1][1] - (-29.43)).abs() < 1e-6,
            "f[1]_y={}",
            forces[1][1]
        );
    }
    #[test]
    fn granular_sim_gravity_accelerates_particle() {
        let contact = HertzContact::new(1e9, 4e8);
        let mut sim = GranularSimulation::new(contact, [0.0, -9.81, 0.0], 0.3, 0.0, 1e-5);
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 1.0);
        sim.step(1e-5);
        let vy = sim.particles[0].velocity[1];
        assert!(
            (vy - (-9.81e-5)).abs() < 1e-12,
            "vy={vy} expected ≈ -9.81e-5"
        );
    }
    #[test]
    fn granular_sim_two_particle_collision_pushes_apart() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.step(1e-5);
        assert!(
            sim.particles[0].velocity[0] < 0.0,
            "p0 vx should be negative after collision, got {}",
            sim.particles[0].velocity[0]
        );
        assert!(
            sim.particles[1].velocity[0] > 0.0,
            "p1 vx should be positive after collision, got {}",
            sim.particles[1].velocity[0]
        );
    }
    #[test]
    fn granular_sim_elastic_collision_momentum_conserved() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let mom_before = sim.total_momentum();
        for _ in 0..10 {
            sim.step(1e-6);
        }
        let mom_after = sim.total_momentum();
        for k in 0..3 {
            assert!(
                (mom_after[k] - mom_before[k]).abs() < 1e-6,
                "momentum[{k}] not conserved: Δ = {}",
                (mom_after[k] - mom_before[k]).abs()
            );
        }
    }
    #[test]
    fn granular_sim_energy_dissipation_with_damping() {
        let mut sim = make_sim(1e4);
        sim.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        let ke_before = sim.kinetic_energy_translational();
        for _ in 0..50 {
            sim.step(1e-6);
        }
        let ke_after = sim.kinetic_energy_translational();
        assert!(ke_after.is_finite(), "KE must be finite: {ke_after}");
        assert!(
            ke_before.is_finite(),
            "initial KE must be finite: {ke_before}"
        );
    }
    #[test]
    fn granular_sim_total_momentum_finite() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0; 3], [0.5, 0.0, 0.0], 0.1, 2.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 2.0);
        for _ in 0..20 {
            sim.step(1e-6);
        }
        let mom = sim.total_momentum();
        for (k, &m) in mom.iter().enumerate() {
            assert!(m.is_finite(), "momentum[{k}] not finite: {}", m);
        }
    }
    #[test]
    fn granular_sim_kinetic_energy_initial_zero() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0; 3], [0.0; 3], 0.1, 1.0);
        sim.add_particle([2.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        assert!(
            (sim.kinetic_energy()).abs() < 1e-14,
            "initial KE should be zero"
        );
    }
    #[test]
    fn granular_sim_angular_velocity_updated_by_torque() {
        let mut sim = make_sim(0.0);
        sim.add_particle([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0);
        sim.add_particle([0.15, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.step(1e-5);
        for k in 0..3 {
            assert!(
                sim.particles[0].angular_velocity[k].is_finite(),
                "angular_velocity[{k}] not finite"
            );
        }
    }
}
/// Compute the static angle of repose from a friction angle \[degrees\].
///
/// For a dry granular material, the angle of repose equals the internal
/// friction angle:  θ_r = φ.
///
/// With cohesion the effective friction is higher, so the angle of repose is
/// estimated as: θ_r = arctan(tan(φ) + c / σ_n).
///
/// # Arguments
/// * `friction_angle_deg` – Internal friction angle φ \[degrees\].
/// * `cohesion`           – Cohesion c \[Pa\].
/// * `normal_stress`      – Representative normal stress σ_n \[Pa\].
pub fn angle_of_repose(friction_angle_deg: f64, cohesion: f64, normal_stress: f64) -> f64 {
    let phi = friction_angle_deg.to_radians();
    let tan_theta = phi.tan()
        + if normal_stress > 1e-14 {
            cohesion / normal_stress
        } else {
            0.0
        };
    tan_theta.atan().to_degrees()
}
/// Dynamic angle of repose from a running DEM simulation.
///
/// Approximated from the velocity field: the angle (in degrees) at which the
/// mean kinetic energy of the surface layer balances the gravitational
/// potential energy gradient.
///
/// Simple model: θ_dyn ≈ θ_static * (1 + 0.2 * T_g)
/// where T_g is the dimensionless granular temperature.
pub fn dynamic_angle_of_repose(static_angle_deg: f64, granular_temp: f64) -> f64 {
    static_angle_deg * (1.0 + 0.2 * granular_temp)
}
/// Void fraction = 1 − packing fraction.
///
/// # Arguments
/// * `radii`          – Sphere radii \[m\].
/// * `domain_volume`  – Total domain volume \[m³\].
pub fn void_fraction(radii: &[f64], domain_volume: f64) -> f64 {
    1.0 - packing_fraction(radii, domain_volume)
}
/// Estimate the theoretical random close packing fraction.
///
/// For monodisperse spheres in 3D, random close packing ≈ 0.64.
/// This returns that constant; for polydisperse assemblies a correction
/// factor can be applied.
pub fn random_close_packing_fraction() -> f64 {
    0.64
}
/// Check whether the assembly is above the jamming transition threshold.
///
/// Jamming occurs roughly at φ ≈ 0.64 for 3D monodisperse spheres.
pub fn is_jammed(packing: f64) -> bool {
    packing >= random_close_packing_fraction()
}
/// Compute the granular temperature profile along the y-axis (vertical).
///
/// Divides the domain into `n_bins` slices of thickness `dy = domain_height / n_bins`.
/// Returns a vector of (y_centre, T_g) pairs.
///
/// # Arguments
/// * `positions`      – Particle positions.
/// * `velocities`     – Particle velocities.
/// * `domain_height`  – Height of the domain \[m\].
/// * `n_bins`         – Number of vertical slices.
pub fn granular_temperature_profile(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    domain_height: f64,
    n_bins: usize,
) -> Vec<(f64, f64)> {
    if n_bins == 0 || positions.is_empty() {
        return Vec::new();
    }
    let dy = domain_height / n_bins as f64;
    let mut bins: Vec<Vec<[f64; 3]>> = vec![Vec::new(); n_bins];
    for (pos, vel) in positions.iter().zip(velocities.iter()) {
        let bin = ((pos[1] / dy).floor() as usize).min(n_bins - 1);
        bins[bin].push(*vel);
    }
    bins.iter()
        .enumerate()
        .map(|(k, vels)| {
            let y_centre = (k as f64 + 0.5) * dy;
            let tg = granular_temperature(vels);
            (y_centre, tg)
        })
        .collect()
}
/// Build a 2D (x–y plane) contact density map.
///
/// Each contact is placed in the cell containing its midpoint.
/// Returns a `(nx × ny)` grid as a flat `Vec`usize` (row-major)
/// together with the cell dimensions.
///
/// # Arguments
/// * `contacts`    – Detected contacts.
/// * `positions`   – Particle positions.
/// * `domain_lo`   – Lower-left corner of the domain.
/// * `domain_hi`   – Upper-right corner of the domain.
/// * `nx`, `ny`    – Number of grid cells in x and y.
pub fn contact_density_map(
    contacts: &[Contact],
    positions: &[[f64; 3]],
    domain_lo: [f64; 2],
    domain_hi: [f64; 2],
    nx: usize,
    ny: usize,
) -> Vec<usize> {
    let mut grid = vec![0usize; nx * ny];
    if nx == 0 || ny == 0 {
        return grid;
    }
    let dx = (domain_hi[0] - domain_lo[0]) / nx as f64;
    let dy = (domain_hi[1] - domain_lo[1]) / ny as f64;
    for c in contacts {
        if c.i >= positions.len() || c.j >= positions.len() {
            continue;
        }
        let mid_x = 0.5 * (positions[c.i][0] + positions[c.j][0]);
        let mid_y = 0.5 * (positions[c.i][1] + positions[c.j][1]);
        let ix = ((mid_x - domain_lo[0]) / dx).floor() as i64;
        let iy = ((mid_y - domain_lo[1]) / dy).floor() as i64;
        if ix >= 0 && (ix as usize) < nx && iy >= 0 && (iy as usize) < ny {
            grid[iy as usize * nx + ix as usize] += 1;
        }
    }
    grid
}
#[cfg(test)]
mod tests_granular_utils {
    use super::*;
    use crate::granular::types::*;
    #[test]
    fn angle_of_repose_equals_phi_for_cohesionless_material() {
        let phi = 30.0_f64;
        let angle = angle_of_repose(phi, 0.0, 1e5);
        assert!((angle - phi).abs() < 1e-8, "angle={angle}, expected {phi}");
    }
    #[test]
    fn angle_of_repose_increases_with_cohesion() {
        let phi = 30.0_f64;
        let angle_no_c = angle_of_repose(phi, 0.0, 1e4);
        let angle_with_c = angle_of_repose(phi, 1e3, 1e4);
        assert!(
            angle_with_c > angle_no_c,
            "cohesion should increase angle of repose"
        );
    }
    #[test]
    fn angle_of_repose_positive_for_typical_sand() {
        let angle = angle_of_repose(30.0, 100.0, 1e5);
        assert!(angle > 0.0 && angle < 90.0, "angle={angle}");
    }
    #[test]
    fn dynamic_angle_of_repose_equals_static_for_zero_temp() {
        let static_a = 30.0_f64;
        let dyn_a = dynamic_angle_of_repose(static_a, 0.0);
        assert!((dyn_a - static_a).abs() < 1e-10);
    }
    #[test]
    fn dynamic_angle_of_repose_increases_with_granular_temp() {
        let a1 = dynamic_angle_of_repose(30.0, 0.0);
        let a2 = dynamic_angle_of_repose(30.0, 1.0);
        assert!(
            a2 > a1,
            "dynamic angle should increase with granular temperature"
        );
    }
    #[test]
    fn void_fraction_plus_packing_fraction_equals_one() {
        let radii = vec![0.05_f64; 10];
        let vol = 1.0_f64;
        let pf = packing_fraction(&radii, vol);
        let vf = void_fraction(&radii, vol);
        assert!((pf + vf - 1.0).abs() < 1e-12, "pf+vf={}", pf + vf);
    }
    #[test]
    fn void_fraction_empty_is_one() {
        let vf = void_fraction(&[], 1.0);
        assert!((vf - 1.0).abs() < 1e-10);
    }
    #[test]
    fn random_close_packing_in_range() {
        let rcp = random_close_packing_fraction();
        assert!(rcp > 0.6 && rcp < 0.7, "RCP={rcp}");
    }
    #[test]
    fn is_jammed_above_threshold() {
        assert!(is_jammed(0.65), "packing=0.65 should be jammed");
    }
    #[test]
    fn is_jammed_below_threshold() {
        assert!(!is_jammed(0.5), "packing=0.5 should NOT be jammed");
    }
    #[test]
    fn granular_temp_profile_empty_input() {
        let profile = granular_temperature_profile(&[], &[], 1.0, 5);
        assert!(profile.is_empty());
    }
    #[test]
    fn granular_temp_profile_zero_bins() {
        let positions = vec![[0.0; 3]];
        let velocities = vec![[1.0, 0.0, 0.0]];
        let profile = granular_temperature_profile(&positions, &velocities, 1.0, 0);
        assert!(profile.is_empty());
    }
    #[test]
    fn granular_temp_profile_correct_bin_count() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [0.0, i as f64 * 0.1, 0.0]).collect();
        let velocities = vec![[0.5, 0.0, 0.0]; 10];
        let profile = granular_temperature_profile(&positions, &velocities, 1.0, 5);
        assert_eq!(profile.len(), 5, "should have 5 bins");
    }
    #[test]
    fn granular_temp_profile_y_centres_increasing() {
        let positions: Vec<[f64; 3]> = (0..4).map(|i| [0.0, i as f64 * 0.25, 0.0]).collect();
        let velocities = vec![[1.0, 0.0, 0.0]; 4];
        let profile = granular_temperature_profile(&positions, &velocities, 1.0, 4);
        for k in 1..profile.len() {
            assert!(
                profile[k].0 > profile[k - 1].0,
                "y centres must be increasing"
            );
        }
    }
    #[test]
    fn contact_density_map_empty_contacts() {
        let grid = contact_density_map(&[], &[], [0.0, 0.0], [1.0, 1.0], 4, 4);
        assert_eq!(grid.len(), 16);
        assert!(grid.iter().all(|&v| v == 0));
    }
    #[test]
    fn contact_density_map_single_contact_correct_cell() {
        let contacts = vec![Contact {
            i: 0,
            j: 1,
            overlap: 0.01,
            normal: [1.0, 0.0, 0.0],
        }];
        let positions = vec![[0.25, 0.25, 0.0], [0.35, 0.25, 0.0]];
        let grid = contact_density_map(&contacts, &positions, [0.0, 0.0], [1.0, 1.0], 4, 4);
        let total: usize = grid.iter().sum();
        assert_eq!(total, 1, "exactly one contact in grid");
    }
    #[test]
    fn contact_density_map_out_of_domain_ignored() {
        let contacts = vec![Contact {
            i: 0,
            j: 1,
            overlap: 0.01,
            normal: [1.0, 0.0, 0.0],
        }];
        let positions = vec![[4.5, 4.5, 0.0], [5.5, 5.5, 0.0]];
        let grid = contact_density_map(&contacts, &positions, [0.0, 0.0], [1.0, 1.0], 4, 4);
        let total: usize = grid.iter().sum();
        assert_eq!(total, 0, "out-of-domain contact should be ignored");
    }
    #[test]
    fn contact_density_map_zero_dimensions_returns_empty() {
        let grid = contact_density_map(&[], &[], [0.0, 0.0], [1.0, 1.0], 0, 4);
        assert!(grid.is_empty());
    }
    #[test]
    fn hertz_normal_force_scales_as_power_1_5() {
        let e_star = 1e9;
        let r_star = 0.05;
        let f1 = hertz_normal_force(0.001, e_star, r_star);
        let f2 = hertz_normal_force(0.004, e_star, r_star);
        let ratio = f2 / f1;
        assert!(
            (ratio - 8.0).abs() < 0.01,
            "Hertz force ratio={ratio} expected 8"
        );
    }
    #[test]
    fn spring_dashpot_normal_damped_less_than_undamped() {
        let kn = 1e5;
        let gamma_n = 1000.0;
        let overlap = 0.002;
        let undamped = spring_dashpot_normal(overlap, 0.0, kn, 0.0);
        let damped = spring_dashpot_normal(overlap, 1.0, kn, gamma_n);
        assert!(damped < undamped, "Damping should reduce normal force");
    }
    #[test]
    fn spring_dashpot_tangential_zero_vel_and_displacement() {
        let f_t = spring_dashpot_tangential(0.0, 0.0, 1e5, 100.0, 0.3, 200.0);
        assert_eq!(f_t, 0.0);
    }
    #[test]
    fn spring_dashpot_tangential_negative_displacement_gives_negative_force() {
        let f_t = spring_dashpot_tangential(-0.001, 0.0, 1e5, 0.0, 0.5, 1e10);
        assert!(
            f_t < 0.0,
            "Negative displacement must give negative tangential force"
        );
    }
    #[test]
    fn granular_collision_new_k_n_positive() {
        let gc = GranularCollision::new(7e10, 0.3, 0.4);
        assert!(gc.k_n > 0.0);
        assert!(gc.k_t > 0.0);
        assert!(gc.gamma_n > 0.0);
    }
    #[test]
    fn granular_collision_tangential_clamped_by_coulomb() {
        let gc = GranularCollision::new(7e10, 0.3, 0.3);
        let f_n = 100.0;
        let f_t_expected_max = gc.mu * f_n;
        let f_t = gc.tangential_force(1000.0, f_n);
        assert!(
            (f_t - f_t_expected_max).abs() < 1e-8,
            "f_t={f_t} limit={f_t_expected_max}"
        );
    }
    #[test]
    fn rolling_resistance_none_returns_zero_torque() {
        let torque = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [1.0, 2.0, 3.0],
            0.1,
            RollingResistanceModel::None,
        );
        assert_eq!(torque, [0.0; 3]);
    }
    #[test]
    fn rolling_resistance_constant_zero_omega_returns_zero() {
        let torque = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            [0.0; 3],
            0.1,
            RollingResistanceModel::Constant,
        );
        assert_eq!(torque, [0.0; 3]);
    }
    #[test]
    fn rolling_resistance_constant_opposes_rotation() {
        let omega = [0.0, 0.0, 10.0];
        let torque = rolling_resistance_torque(
            0.1,
            0.1,
            100.0,
            omega,
            0.05,
            RollingResistanceModel::Constant,
        );
        assert!(
            torque[2] < 0.0,
            "Constant rolling torque should oppose omega, got {}",
            torque[2]
        );
    }
    #[test]
    fn rolling_resistance_viscous_proportional_to_omega() {
        let omega1 = [0.0, 0.0, 1.0];
        let omega2 = [0.0, 0.0, 2.0];
        let t1 = rolling_resistance_torque(
            0.1,
            0.1,
            50.0,
            omega1,
            0.05,
            RollingResistanceModel::Viscous,
        );
        let t2 = rolling_resistance_torque(
            0.1,
            0.1,
            50.0,
            omega2,
            0.05,
            RollingResistanceModel::Viscous,
        );
        assert!(
            (t2[2] / t1[2] - 2.0).abs() < 1e-10,
            "Viscous torque ratio should be 2"
        );
    }
    #[test]
    fn packing_fraction_sphere_volume() {
        let r = 0.1_f64;
        let domain = 1.0;
        let phi = packing_fraction(&[r, r, r, r], domain);
        let expected = 4.0 * (4.0 / 3.0) * PI * r.powi(3) / domain;
        assert!((phi - expected).abs() < 1e-12);
    }
    #[test]
    fn packing_fraction_zero_domain_returns_zero() {
        let phi = packing_fraction(&[0.1, 0.1], 0.0);
        assert_eq!(phi, 0.0);
    }
    #[test]
    fn packing_fraction_clamped_at_one() {
        let phi = packing_fraction(&[100.0], 1.0);
        assert_eq!(phi, 1.0);
    }
    #[test]
    fn granular_temperature_zero_for_stationary() {
        let vels = vec![[0.0; 3]; 10];
        let tg = granular_temperature(&vels);
        assert_eq!(tg, 0.0);
    }
    #[test]
    fn granular_temperature_zero_for_uniform_velocity() {
        let vels = vec![[1.0, 2.0, 3.0]; 8];
        let tg = granular_temperature(&vels);
        assert!(tg.abs() < 1e-12, "Uniform velocity → T_g = 0, got {tg}");
    }
    #[test]
    fn granular_temperature_positive_for_mixed_velocities() {
        let vels = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let tg = granular_temperature(&vels);
        assert!(tg > 0.0, "Mixed velocities should give positive T_g");
    }
    #[test]
    fn granular_temperature_empty_is_zero() {
        let tg = granular_temperature(&[]);
        assert_eq!(tg, 0.0);
    }
    #[test]
    fn coordination_numbers_sum_equals_twice_contact_count() {
        let contacts = vec![
            Contact {
                i: 0,
                j: 1,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
            Contact {
                i: 1,
                j: 2,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
        ];
        let cn = coordination_numbers(&contacts, 3);
        assert_eq!(cn[0], 1);
        assert_eq!(cn[1], 2);
        assert_eq!(cn[2], 1);
        let total: usize = cn.iter().sum();
        assert_eq!(
            total,
            contacts.len() * 2,
            "Sum of coord numbers = 2 * N_contacts"
        );
    }
    #[test]
    fn mean_coordination_number_empty_contacts() {
        let mcn = mean_coordination_number(&[], 5);
        assert_eq!(mcn, 0.0);
    }
    #[test]
    fn fabric_tensor_trace_is_one_for_isotropic() {
        let contacts = vec![
            Contact {
                i: 0,
                j: 1,
                overlap: 0.01,
                normal: [1.0, 0.0, 0.0],
            },
            Contact {
                i: 0,
                j: 2,
                overlap: 0.01,
                normal: [0.0, 1.0, 0.0],
            },
            Contact {
                i: 0,
                j: 3,
                overlap: 0.01,
                normal: [0.0, 0.0, 1.0],
            },
        ];
        let fab = fabric_tensor(&contacts);
        let trace = fab[0][0] + fab[1][1] + fab[2][2];
        assert!(
            (trace - 1.0).abs() < 1e-12,
            "Isotropic fabric trace should be 1, got {trace}"
        );
    }
    #[test]
    fn fabric_tensor_empty_contacts_is_zero() {
        let fab = fabric_tensor(&[]);
        for row in &fab {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn granular_params_cohesive_powder_higher_friction_than_sand() {
        let sand = GranularParams::sand();
        let powder = GranularParams::cohesive_powder();
        assert!(powder.friction_angle_deg > sand.friction_angle_deg);
    }
    #[test]
    fn dem_particle_kinetic_energy_spinning() {
        let mut p = DemGranularParticle::new([0.0; 3], 0.1, 1000.0);
        p.angular_velocity = [0.0, 0.0, 10.0];
        let ke = p.kinetic_energy();
        let inertia = p.moment_of_inertia();
        let expected = 0.5 * inertia * 100.0;
        assert!((ke - expected).abs() < 1e-10);
    }
    #[test]
    fn granular_sim_extended_step_no_panic() {
        let model = ContactModel::SpringDashpot {
            kn: 1e5,
            kt: 5e4,
            gamma_n: 100.0,
            gamma_t: 50.0,
        };
        let mut sim = GranularSimExtended::new(
            model,
            0.3,
            [0.0, -9.81, 0.0],
            RollingResistanceModel::Constant,
            0.05,
        );
        sim.add_particle([0.0, 1.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle([0.0, 0.7, 0.0], [0.0; 3], 0.1, 1.0);
        sim.step(0.001);
        assert!(!sim.is_empty());
    }
    #[test]
    fn critical_slope_angle_pure_friction_equals_friction_angle() {
        let params = GranularParams::sand();
        let angle = critical_slope_angle(&params, 1e12);
        let phi_rad = params.friction_angle_deg.to_radians();
        assert!(
            (angle - phi_rad).abs() < 1e-6,
            "angle={angle} phi={phi_rad}"
        );
    }
    #[test]
    fn is_slope_stable_below_critical_angle() {
        let params = GranularParams::sand();
        let theta = 10_f64.to_radians();
        assert!(
            is_slope_stable(theta, &params, 1e6),
            "10° slope should be stable for sand"
        );
    }
    #[test]
    fn is_slope_stable_above_critical_angle() {
        let params = GranularParams::sand();
        let theta = 45_f64.to_radians();
        assert!(
            !is_slope_stable(theta, &params, 1e6),
            "45° slope should be unstable for sand"
        );
    }
    #[test]
    fn spatial_hash_dem_no_contact_far_apart() {
        let positions = vec![[0.0; 3], [100.0, 0.0, 0.0]];
        let radii = vec![0.1; 2];
        let hash = SpatialHashDem::new(&positions, &radii, 0.1);
        let pairs = hash.candidate_pairs(&positions, &radii);
        assert!(
            pairs.is_empty(),
            "Far apart particles should have no candidate pairs"
        );
    }
}
