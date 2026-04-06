#![allow(clippy::needless_range_loop)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{
    AdiabaticEos, ArtificialViscosity, AstroParticle, BBox3, CoolingFunction, GravitationalSph,
    OctreeNode,
};

/// Gravitational constant G in CGS: 6.674e-8 cm^3 g^-1 s^-2.
pub const G_CGS: f64 = 6.674e-8;
/// Boltzmann constant k_B in CGS: 1.381e-16 erg K^-1.
pub const K_BOLTZ_CGS: f64 = 1.381e-16;
/// Proton mass in CGS: 1.673e-24 g.
pub const M_PROTON_CGS: f64 = 1.673e-24;
/// Solar mass in CGS: 1.989e33 g.
pub const M_SUN_CGS: f64 = 1.989e33;
/// Parsec in CGS: 3.086e18 cm.
pub const PARSEC_CGS: f64 = 3.086e18;
/// Year in seconds.
pub const YEAR_S: f64 = 3.156e7;
/// Dot product of two 3-vectors.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean norm of a 3-vector.
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Distance between two 3D points.
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    norm3(d)
}
/// Vector subtraction.
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Cross product of two 3-vectors.
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Cubic spline kernel W(r, h) in 3D.
///
/// Normalized for 3D: sigma = 1 / (pi * h^3).
pub fn cubic_spline_w(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        0.0
    }
}
/// Gradient magnitude of the cubic spline kernel dW/dr.
pub fn cubic_spline_dw(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h * h);
    if q < 1e-14 {
        0.0
    } else if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t)
    } else {
        0.0
    }
}
/// Compute SPH density for all particles using direct summation.
pub fn compute_density(particles: &mut [AstroParticle]) {
    let n = particles.len();
    let positions: Vec<[f64; 3]> = particles.iter().map(|p| p.pos).collect();
    let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let hs: Vec<f64> = particles.iter().map(|p| p.h).collect();
    for i in 0..n {
        let mut rho = 0.0;
        for j in 0..n {
            let r = dist3(positions[i], positions[j]);
            let h_avg = 0.5 * (hs[i] + hs[j]);
            rho += masses[j] * cubic_spline_w(r, h_avg);
        }
        particles[i].density = rho;
    }
}
/// Compute velocity divergence for a particle.
pub fn compute_div_v(idx: usize, particles: &[AstroParticle]) -> f64 {
    let pi = &particles[idx];
    let mut div_v = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == idx {
            continue;
        }
        let rij = sub3(pi.pos, pj.pos);
        let r = norm3(rij);
        if r < 1e-30 {
            continue;
        }
        let h_avg = 0.5 * (pi.h + pj.h);
        let dw = cubic_spline_dw(r, h_avg);
        let vij = sub3(pi.vel, pj.vel);
        let grad_w_factor = dw / r;
        let v_dot_r = dot3(vij, rij);
        if pj.density > 1e-30 {
            div_v -= pj.mass / pj.density * v_dot_r * grad_w_factor;
        }
    }
    div_v
}
/// Compute curl of velocity magnitude for a particle.
pub fn compute_curl_v_mag(idx: usize, particles: &[AstroParticle]) -> f64 {
    let pi = &particles[idx];
    let mut curl = [0.0; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == idx {
            continue;
        }
        let rij = sub3(pi.pos, pj.pos);
        let r = norm3(rij);
        if r < 1e-30 || pj.density < 1e-30 {
            continue;
        }
        let h_avg = 0.5 * (pi.h + pj.h);
        let dw = cubic_spline_dw(r, h_avg);
        let vij = sub3(pi.vel, pj.vel);
        let grad_factor = pj.mass / pj.density * dw / r;
        let cross = cross3(vij, rij);
        curl[0] += grad_factor * cross[0];
        curl[1] += grad_factor * cross[1];
        curl[2] += grad_factor * cross[2];
    }
    norm3(curl)
}
/// Leapfrog (kick-drift-kick) integrator for astrophysical SPH.
///
/// Advances particles by one time step dt.
pub fn leapfrog_step(
    particles: &mut [AstroParticle],
    dt: f64,
    gravity: &GravitationalSph,
    eos: &AdiabaticEos,
    viscosity: &ArtificialViscosity,
) {
    let _ = viscosity;
    let n = particles.len();
    for i in 0..n {
        if particles[i].is_sink {
            continue;
        }
        for d in 0..3 {
            particles[i].vel[d] += 0.5 * dt * particles[i].acc[d];
        }
        particles[i].internal_energy += 0.5 * dt * particles[i].du_dt;
        if particles[i].internal_energy < 0.0 {
            particles[i].internal_energy = 0.0;
        }
    }
    for i in 0..n {
        for d in 0..3 {
            particles[i].pos[d] += dt * particles[i].vel[d];
        }
    }
    compute_density(particles);
    eos.apply(particles);
    for p in particles.iter_mut() {
        p.acc = [0.0; 3];
    }
    gravity.compute_gravity(particles);
    for i in 0..n {
        if particles[i].is_sink {
            continue;
        }
        for d in 0..3 {
            particles[i].vel[d] += 0.5 * dt * particles[i].acc[d];
        }
        particles[i].internal_energy += 0.5 * dt * particles[i].du_dt;
        if particles[i].internal_energy < 0.0 {
            particles[i].internal_energy = 0.0;
        }
    }
}
/// Compute an adaptive time step based on CFL and acceleration criteria.
pub fn adaptive_timestep(particles: &[AstroParticle], gamma: f64, cfl: f64) -> f64 {
    let mut dt_min = f64::MAX;
    for p in particles {
        if p.is_sink {
            continue;
        }
        let cs = p.sound_speed(gamma);
        let v = p.speed();
        let denom = cs + v;
        if denom > 1e-30 {
            let dt_cfl = cfl * p.h / denom;
            dt_min = dt_min.min(dt_cfl);
        }
        let a_mag = norm3(p.acc);
        if a_mag > 1e-30 {
            let dt_acc = (p.h / a_mag).sqrt();
            dt_min = dt_min.min(dt_acc);
        }
    }
    dt_min
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::astrophysical_sph::types::*;
    #[test]
    fn test_adiabatic_eos_pressure() {
        let eos = AdiabaticEos::monatomic();
        let p = eos.pressure(1.0, 1.0);
        assert!((p - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_eos_sound_speed() {
        let eos = AdiabaticEos::monatomic();
        let cs = eos.sound_speed(1.0);
        let expected = (5.0_f64 / 3.0 * 2.0 / 3.0).sqrt();
        assert!((cs - expected).abs() < 1e-10);
    }
    #[test]
    fn test_eos_temperature_roundtrip() {
        let eos = AdiabaticEos::monatomic();
        let t = 1000.0;
        let mu = 1.0;
        let u = eos.internal_energy_from_temperature(t, mu);
        let t_back = eos.temperature(u, mu);
        assert!((t_back - t).abs() / t < 1e-10);
    }
    #[test]
    fn test_eos_polytropic() {
        let eos = AdiabaticEos::monatomic();
        let k = 1.0;
        let rho = 2.0;
        let p = eos.polytropic_pressure(rho, k);
        assert!((p - rho.powf(5.0 / 3.0)).abs() < 1e-10);
    }
    #[test]
    fn test_jeans_mass_formula() {
        let cs = 1.0;
        let rho = 1.0;
        let g = 1.0;
        let jeans = JeansInstability::new(cs, rho, g);
        let lj = jeans.jeans_length();
        let mj = jeans.jeans_mass();
        let expected_lj = (PI).sqrt();
        assert!((lj - expected_lj).abs() < 1e-10);
        let expected_mj = (PI / 6.0) * expected_lj.powi(3);
        assert!((mj - expected_mj).abs() < 1e-10);
    }
    #[test]
    fn test_free_fall_time() {
        let rho = 1.0;
        let g = 1.0;
        let jeans = JeansInstability::new(1.0, rho, g);
        let tff = jeans.free_fall_time();
        let expected = (3.0 * PI / 32.0).sqrt();
        assert!((tff - expected).abs() < 1e-10);
    }
    #[test]
    fn test_jeans_wavenumber() {
        let cs = 2.0;
        let rho = 1.0;
        let g = 1.0;
        let jeans = JeansInstability::new(cs, rho, g);
        let kj = jeans.jeans_wavenumber();
        let expected = (4.0 * PI * g * rho).sqrt() / cs;
        assert!((kj - expected).abs() < 1e-10);
    }
    #[test]
    fn test_jeans_stability() {
        let jeans = JeansInstability::new(1.0, 1.0, 1.0);
        let kj = jeans.jeans_wavenumber();
        assert!(jeans.is_unstable(kj * 0.5));
        assert!(!jeans.is_unstable(kj * 2.0));
    }
    #[test]
    fn test_jeans_growth_rate() {
        let jeans = JeansInstability::new(1.0, 1.0, 1.0);
        let kj = jeans.jeans_wavenumber();
        let omega = jeans.growth_rate(kj);
        assert!(omega.abs() < 1e-6, "growth rate at kJ: {omega}");
    }
    #[test]
    fn test_gravitational_potential() {
        let grav = GravitationalSph::new(1.0, 0.0);
        let phi = grav.potential_pair(1.0, 1.0);
        assert!((phi - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_gravitational_acceleration() {
        let grav = GravitationalSph::new(1.0, 0.0);
        let acc = grav.acceleration_pair([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        assert!((acc[0] - (-1.0)).abs() < 1e-10);
        assert!(acc[1].abs() < 1e-10);
    }
    #[test]
    fn test_gravitational_softening() {
        let grav = GravitationalSph::new(1.0, 0.1);
        let acc = grav.acceleration_pair([0.0, 0.0, 0.0], [0.001, 0.0, 0.0], 1.0);
        assert!(acc[0].is_finite());
    }
    #[test]
    fn test_artificial_viscosity_approaching() {
        let visc = ArtificialViscosity::standard();
        let pi_ij = visc.compute_pi(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            1.0,
            1.0,
            1.0,
            1.0,
            0.5,
            0.5,
        );
        assert!(pi_ij > 0.0);
    }
    #[test]
    fn test_artificial_viscosity_receding() {
        let visc = ArtificialViscosity::standard();
        let pi_ij = visc.compute_pi(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1.0,
            1.0,
            1.0,
            1.0,
            0.5,
            0.5,
        );
        assert!(pi_ij.abs() < 1e-14);
    }
    #[test]
    fn test_balsara_switch() {
        let visc = ArtificialViscosity::standard();
        let f = visc.balsara_switch(-10.0, 0.01, 1.0, 0.1);
        assert!(f > 0.9);
        let f = visc.balsara_switch(-0.01, 10.0, 1.0, 0.1);
        assert!(f < 0.1);
    }
    #[test]
    fn test_von_neumann_richtmyer() {
        let visc = ArtificialViscosity::standard();
        let q = visc.von_neumann_richtmyer(1.0, -2.0, 0.1);
        assert!((q - 0.08).abs() < 1e-10);
    }
    #[test]
    fn test_keplerian_velocity() {
        let disk = DiskDynamics::new(1.0, 1.0);
        let v = disk.keplerian_velocity(1.0);
        assert!((v - 1.0).abs() < 1e-10);
        let v4 = disk.keplerian_velocity(4.0);
        assert!((v4 - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_toomre_q() {
        let disk = DiskDynamics::new(1.0, 1.0);
        let q = disk.toomre_q(1.0, 1.0, 1.0);
        let expected = 1.0 / PI;
        assert!((q - expected).abs() < 1e-10);
    }
    #[test]
    fn test_toomre_stability() {
        let disk = DiskDynamics::new(1.0, 1.0);
        assert!(disk.is_toomre_stable(100.0, 1.0, 1.0));
        assert!(!disk.is_toomre_stable(0.01, 10.0, 1.0));
    }
    #[test]
    fn test_orbital_period() {
        let disk = DiskDynamics::new(1.0, 1.0);
        let t = disk.orbital_period(1.0);
        assert!((t - 2.0 * PI).abs() < 1e-10);
    }
    #[test]
    fn test_scale_height() {
        let disk = DiskDynamics::new(1.0, 1.0);
        let h = disk.scale_height(1.0, 1.0);
        assert!((h - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_virial_ratio() {
        let analysis = AstroAnalysis::new(1.0, 0.01);
        let particles = vec![
            AstroParticle::new(1.0, [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], 0.0, 0.1),
            AstroParticle::new(1.0, [-0.5, 0.0, 0.0], [0.0, -0.5, 0.0], 0.0, 0.1),
        ];
        let ratio = analysis.virial_ratio(&particles);
        assert!(ratio > 0.0 && ratio.is_finite());
    }
    #[test]
    fn test_energy_conservation_check() {
        let analysis = AstroAnalysis::new(1.0, 0.01);
        let err = analysis.energy_conservation(-10.0, -10.01);
        assert!(err < 0.01);
    }
    #[test]
    fn test_center_of_mass() {
        let analysis = AstroAnalysis::new(1.0, 0.01);
        let particles = vec![
            AstroParticle::new(1.0, [0.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
            AstroParticle::new(1.0, [2.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
        ];
        let com = analysis.center_of_mass(&particles);
        assert!((com[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_total_angular_momentum() {
        let analysis = AstroAnalysis::new(1.0, 0.01);
        let particles = vec![AstroParticle::new(
            1.0,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.0,
            0.1,
        )];
        let l = analysis.total_angular_momentum(&particles);
        assert!((l[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_particle_kinetic_energy() {
        let p = AstroParticle::new(2.0, [0.0; 3], [3.0, 4.0, 0.0], 0.0, 0.1);
        let ke = p.kinetic_energy();
        assert!((ke - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_particle_thermal_energy() {
        let p = AstroParticle::new(3.0, [0.0; 3], [0.0; 3], 5.0, 0.1);
        assert!((p.thermal_energy() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_cubic_spline_normalization() {
        let h = 1.0;
        let n = 100;
        let dr = 2.0 * h / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) * dr;
            integral += cubic_spline_w(r, h) * 4.0 * PI * r * r * dr;
        }
        assert!(
            (integral - 1.0).abs() < 0.05,
            "Kernel normalization: {integral}"
        );
    }
    #[test]
    fn test_density_computation() {
        let mut particles = vec![
            AstroParticle::new(1.0, [0.0, 0.0, 0.0], [0.0; 3], 1.0, 1.0),
            AstroParticle::new(1.0, [0.5, 0.0, 0.0], [0.0; 3], 1.0, 1.0),
        ];
        compute_density(&mut particles);
        assert!(particles[0].density > 0.0);
        assert!(particles[1].density > 0.0);
    }
    #[test]
    fn test_sink_particle_creation() {
        let sf = StarFormation::new(10.0, 0.1, 4.0, 0.5);
        let gas = AstroParticle::new(1.0, [0.0; 3], [0.0; 3], 0.01, 0.1);
        let sink = sf.create_sink(&gas);
        assert!(sink.is_sink);
        assert!((sink.mass - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_star_formation_rate() {
        let sf = StarFormation::new(10.0, 0.1, 4.0, 0.1);
        let rate = sf.star_formation_rate(1.0, 1.0);
        let tff = (3.0 * PI / 32.0).sqrt();
        let expected = 0.1 / tff;
        assert!((rate - expected).abs() < 1e-10);
    }
    #[test]
    fn test_alpha_viscosity() {
        let disk = DiskDynamics::new(1.0, 1.0);
        let nu = disk.alpha_viscosity(0.1, 1.0, 1.0);
        assert!((nu - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_half_mass_radius() {
        let analysis = AstroAnalysis::new(1.0, 0.01);
        let particles = vec![
            AstroParticle::new(1.0, [0.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
            AstroParticle::new(1.0, [1.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
            AstroParticle::new(1.0, [2.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
            AstroParticle::new(1.0, [3.0, 0.0, 0.0], [0.0; 3], 0.0, 0.1),
        ];
        let r_half = analysis.half_mass_radius(&particles);
        assert!(r_half > 0.0);
    }
    #[test]
    fn test_bonnor_ebert_mass() {
        let jeans = JeansInstability::new(1.0, 1.0, 1.0);
        let m_be = jeans.bonnor_ebert_mass(1.0);
        assert!(m_be > 0.0);
        assert!(m_be.is_finite());
    }
    #[test]
    fn test_adaptive_timestep() {
        let mut p = AstroParticle::new(1.0, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.1);
        p.density = 1.0;
        p.pressure = 2.0 / 3.0;
        p.acc = [0.0, 0.0, -1.0];
        let dt = adaptive_timestep(&[p], 5.0 / 3.0, 0.3);
        assert!(dt > 0.0 && dt.is_finite());
    }
}
/// Build a Barnes-Hut octree from particle positions and masses.
#[allow(dead_code)]
pub fn build_barnes_hut_tree(
    positions: &[[f64; 3]],
    masses: &[f64],
    max_depth: usize,
) -> OctreeNode {
    assert_eq!(positions.len(), masses.len());
    let mut lo = [f64::MAX; 3];
    let mut hi = [f64::MIN; 3];
    for p in positions {
        for d in 0..3 {
            if p[d] < lo[d] {
                lo[d] = p[d];
            }
            if p[d] > hi[d] {
                hi[d] = p[d];
            }
        }
    }
    for d in 0..3 {
        let margin = 0.01 * (hi[d] - lo[d]).max(1e-10);
        lo[d] -= margin;
        hi[d] += margin;
    }
    let bbox = BBox3::new(lo, hi);
    let mut root = OctreeNode::new(bbox);
    for (i, (p, &m)) in positions.iter().zip(masses.iter()).enumerate() {
        root.insert(i, *p, m, max_depth);
    }
    root
}
/// Compute gravitational accelerations on all particles using Barnes-Hut.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn barnes_hut_gravity(
    positions: &[[f64; 3]],
    masses: &[f64],
    theta: f64,
    softening: f64,
    g: f64,
    max_depth: usize,
) -> Vec<[f64; 3]> {
    let tree = build_barnes_hut_tree(positions, masses, max_depth);
    positions
        .iter()
        .map(|p| tree.gravity_at(*p, theta, softening, g))
        .collect()
}
/// Compute the cooling rate per unit volume \[erg cm^-3 s^-1\] given
/// temperature T \[K\] and number density n \[cm^-3\].
#[allow(dead_code)]
pub fn cooling_rate(func: CoolingFunction, temperature: f64, number_density: f64) -> f64 {
    let n2 = number_density * number_density;
    match func {
        CoolingFunction::None => 0.0,
        CoolingFunction::Bremsstrahlung => 1.43e-27 * temperature.sqrt() * n2,
        CoolingFunction::OpticallyThin => {
            if temperature < 1e4 {
                0.0
            } else if temperature < 1e5 {
                1e-22 * (temperature / 1e5).powf(6.0) * n2
            } else if temperature < 1e7 {
                1e-22 * n2
            } else {
                1.43e-27 * temperature.sqrt() * n2
            }
        }
        CoolingFunction::AtomicLine => {
            let peak_t = 1e5;
            let sigma = 0.5;
            let log_ratio = (temperature / peak_t).ln();
            let amplitude = 2e-22 * n2;
            amplitude * (-0.5 * (log_ratio / sigma).powi(2)).exp()
        }
        CoolingFunction::Tabulated => {
            let log_t = temperature.log10();
            let lambda = if log_t < 4.0 {
                0.0
            } else if log_t < 4.5 {
                1e-24 * 10.0_f64.powf(2.0 * (log_t - 4.0))
            } else if log_t < 5.5 {
                1e-22
            } else if log_t < 7.0 {
                1e-22 * 10.0_f64.powf(-0.5 * (log_t - 5.5))
            } else {
                1.43e-27 * temperature.sqrt()
            };
            lambda * n2
        }
    }
}
/// Compute the cooling time t_cool = (3/2) n k_B T / Λ.
#[allow(dead_code)]
pub fn cooling_time(func: CoolingFunction, temperature: f64, number_density: f64) -> f64 {
    let lambda = cooling_rate(func, temperature, number_density);
    if lambda <= 0.0 {
        f64::INFINITY
    } else {
        1.5 * number_density * K_BOLTZ_CGS * temperature / lambda
    }
}
/// Apply radiative cooling to internal energy over a timestep dt.
///
/// Uses implicit (backward-Euler style) subcycling if dt > t_cool.
#[allow(dead_code)]
pub fn apply_cooling(
    func: CoolingFunction,
    internal_energy: f64,
    density: f64,
    mu: f64,
    gamma: f64,
    dt: f64,
) -> f64 {
    if matches!(func, CoolingFunction::None) {
        return internal_energy;
    }
    let n = density / (mu * M_PROTON_CGS);
    let temperature = (gamma - 1.0) * mu * M_PROTON_CGS * internal_energy / K_BOLTZ_CGS;
    if temperature < 10.0 {
        return internal_energy;
    }
    let lambda = cooling_rate(func, temperature, n);
    if lambda <= 0.0 {
        return internal_energy;
    }
    let du = lambda * dt / density;

    (internal_energy - du).max(K_BOLTZ_CGS * 10.0 / ((gamma - 1.0) * mu * M_PROTON_CGS))
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::astrophysical_sph::types::*;
    #[test]
    fn test_bbox_centre() {
        let bb = BBox3::new([0.0; 3], [2.0; 3]);
        let c = bb.centre();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_bbox_contains() {
        let bb = BBox3::new([0.0; 3], [1.0; 3]);
        assert!(bb.contains([0.5, 0.5, 0.5]));
        assert!(!bb.contains([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_bbox_octant() {
        let bb = BBox3::new([0.0; 3], [2.0; 3]);
        assert_eq!(bb.octant([0.5, 0.5, 0.5]), 0);
        assert_eq!(bb.octant([1.5, 0.5, 0.5]), 1);
        assert_eq!(bb.octant([0.5, 1.5, 0.5]), 2);
        assert_eq!(bb.octant([1.5, 1.5, 1.5]), 7);
    }
    #[test]
    fn test_build_barnes_hut_tree() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0, 1.0, 1.0];
        let tree = build_barnes_hut_tree(&positions, &masses, 20);
        assert!((tree.total_mass - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_barnes_hut_gravity_two_particles() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let accs = barnes_hut_gravity(&positions, &masses, 0.5, 0.01, 1.0, 20);
        assert_eq!(accs.len(), 2);
        assert!(accs[0][0] > 0.0);
        assert!(accs[1][0] < 0.0);
    }
    #[test]
    fn test_barnes_hut_gravity_symmetry() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let accs = barnes_hut_gravity(&positions, &masses, 0.0, 0.01, 1.0, 20);
        let mag0 =
            (accs[0][0] * accs[0][0] + accs[0][1] * accs[0][1] + accs[0][2] * accs[0][2]).sqrt();
        let mag1 =
            (accs[1][0] * accs[1][0] + accs[1][1] * accs[1][1] + accs[1][2] * accs[1][2]).sqrt();
        assert!((mag0 - mag1).abs() / mag0.max(1e-30) < 0.1);
    }
    #[test]
    fn test_cooling_none() {
        let rate = cooling_rate(CoolingFunction::None, 1e6, 1.0);
        assert!((rate - 0.0).abs() < 1e-30);
    }
    #[test]
    fn test_cooling_bremsstrahlung_positive() {
        let rate = cooling_rate(CoolingFunction::Bremsstrahlung, 1e7, 1.0);
        assert!(rate > 0.0);
    }
    #[test]
    fn test_cooling_bremsstrahlung_scales_with_sqrt_t() {
        let r1 = cooling_rate(CoolingFunction::Bremsstrahlung, 1e6, 1.0);
        let r2 = cooling_rate(CoolingFunction::Bremsstrahlung, 4e6, 1.0);
        let ratio = r2 / r1;
        assert!((ratio - 2.0).abs() < 0.01);
    }
    #[test]
    fn test_cooling_optically_thin_below_threshold() {
        let rate = cooling_rate(CoolingFunction::OpticallyThin, 1e3, 1.0);
        assert!((rate - 0.0).abs() < 1e-30);
    }
    #[test]
    fn test_cooling_time_finite() {
        let t_cool = cooling_time(CoolingFunction::Bremsstrahlung, 1e7, 1.0);
        assert!(t_cool > 0.0 && t_cool.is_finite());
    }
    #[test]
    fn test_apply_cooling_reduces_energy() {
        let u0 = 1e12;
        let u1 = apply_cooling(
            CoolingFunction::Bremsstrahlung,
            u0,
            1e-24,
            0.6,
            5.0 / 3.0,
            1e10,
        );
        assert!(u1 < u0);
    }
    #[test]
    fn test_apply_cooling_none_preserves() {
        let u0 = 1e12;
        let u1 = apply_cooling(CoolingFunction::None, u0, 1e-24, 0.6, 5.0 / 3.0, 1e10);
        assert!((u1 - u0).abs() < 1e-10);
    }
    #[test]
    fn test_sn_standard_params() {
        let sn = SupernovaFeedback::standard();
        assert!((sn.energy_per_sn - 1e51).abs() < 1e40);
        assert!(sn.coupling_efficiency > 0.0);
    }
    #[test]
    fn test_sn_energy_injected() {
        let sn = SupernovaFeedback::standard();
        let e = sn.energy_injected(M_SUN_CGS);
        assert!(e > 0.0);
    }
    #[test]
    fn test_sn_mass_returned() {
        let sn = SupernovaFeedback::standard();
        let m = sn.mass_returned(M_SUN_CGS);
        assert!((m - 0.1 * M_SUN_CGS).abs() < 1e20);
    }
    #[test]
    fn test_sn_thermal_per_neighbour() {
        let sn = SupernovaFeedback::standard();
        let e_per = sn.thermal_energy_per_neighbour(M_SUN_CGS, 32);
        assert!(e_per > 0.0);
        let e_total = sn.energy_injected(M_SUN_CGS);
        assert!((e_per - e_total / 32.0).abs() < 1e30);
    }
    #[test]
    fn test_sn_sedov_taylor_radius() {
        let sn = SupernovaFeedback::standard();
        let r = sn.sedov_taylor_radius(1e51, 1e-24, 1e4 * YEAR_S);
        assert!(r > 0.0 && r.is_finite());
    }
    #[test]
    fn test_cosmology_planck() {
        let cosmo = Cosmology::planck2018();
        assert!((cosmo.omega_m + cosmo.omega_lambda + cosmo.omega_r - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_redshift_round_trip() {
        let cosmo = Cosmology::planck2018();
        let z = 2.0;
        let a = cosmo.scale_factor_from_z(z);
        let z_back = cosmo.redshift(a);
        assert!((z - z_back).abs() < 1e-12);
    }
    #[test]
    fn test_comoving_physical_round_trip() {
        let cosmo = Cosmology::planck2018();
        let phys = [1.0, 2.0, 3.0];
        let a = 0.5;
        let comov = cosmo.to_comoving(phys, a);
        let phys_back = cosmo.to_physical(comov, a);
        for d in 0..3 {
            assert!((phys[d] - phys_back[d]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_hubble_at_present() {
        let cosmo = Cosmology::planck2018();
        let h1 = cosmo.hubble(1.0);
        assert!((h1 - cosmo.h0).abs() / cosmo.h0 < 0.01);
    }
    #[test]
    fn test_critical_density_positive() {
        let cosmo = Cosmology::planck2018();
        let rho_c = cosmo.critical_density(1.0);
        assert!(rho_c > 0.0 && rho_c.is_finite());
    }
    #[test]
    fn test_comoving_distance_positive() {
        let cosmo = Cosmology::planck2018();
        let d = cosmo.comoving_distance(1.0, 1000);
        assert!(d > 0.0 && d.is_finite());
    }
    #[test]
    fn test_polytrope_pressure() {
        let eos = PolytropicEos::new(1.0, 5.0 / 3.0);
        let p = eos.pressure(2.0);
        let expected = 2.0_f64.powf(5.0 / 3.0);
        assert!((p - expected).abs() < 1e-12);
    }
    #[test]
    fn test_polytrope_sound_speed() {
        let eos = PolytropicEos::new(1.0, 5.0 / 3.0);
        let cs = eos.sound_speed(1.0);
        assert!((cs - (5.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_polytrope_density_from_pressure() {
        let eos = PolytropicEos::new(2.0, 1.5);
        let rho = 3.0;
        let p = eos.pressure(rho);
        let rho_back = eos.density_from_pressure(p);
        assert!((rho - rho_back).abs() < 1e-10);
    }
    #[test]
    fn test_polytrope_from_index() {
        let eos = PolytropicEos::from_index(1.0, 1.5);
        assert!((eos.gamma_poly - 5.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_polytrope_isothermal() {
        let eos = PolytropicEos::isothermal(100.0);
        assert!((eos.gamma_poly - 1.0).abs() < 1e-12);
        let p = eos.pressure(2.0);
        assert!((p - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_polytrope_lane_emden_scale() {
        let eos = PolytropicEos::new(1.0, 5.0 / 3.0);
        let r = eos.lane_emden_scale(1.0, 1.0);
        assert!(r > 0.0 && r.is_finite());
    }
    #[test]
    fn test_polytrope_internal_energy() {
        let eos = PolytropicEos::new(1.0, 5.0 / 3.0);
        let u = eos.internal_energy(1.0);
        assert!((u - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_cooling_tabulated_hot() {
        let rate = cooling_rate(CoolingFunction::Tabulated, 1e8, 1.0);
        assert!(rate > 0.0);
    }
    #[test]
    fn test_cooling_atomic_line_peak() {
        let r_peak = cooling_rate(CoolingFunction::AtomicLine, 1e5, 1.0);
        let r_off = cooling_rate(CoolingFunction::AtomicLine, 1e3, 1.0);
        assert!(r_peak > r_off);
    }
    #[test]
    fn test_sn_momentum_injection() {
        let sn = SupernovaFeedback::standard();
        let p = sn.momentum_injection(1e51, 1e-24);
        assert!(p > 0.0 && p.is_finite());
    }
    #[test]
    fn test_peculiar_velocity() {
        let cosmo = Cosmology::new(70.0, 0.3, 0.7);
        let pv = cosmo.peculiar_velocity([100.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0);
        assert!((pv[0] - 30.0).abs() < 1.0);
    }
}
