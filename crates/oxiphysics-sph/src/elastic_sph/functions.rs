//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::{ConstitutiveModel, ElasticParticle};

#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Matrix-matrix product C = A·B (3×3).
pub(super) fn mat_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Transpose of a 3×3 matrix.
pub(super) fn mat_transpose(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}
/// Trace of a 3×3 matrix.
pub(super) fn mat_trace(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}
/// Determinant of a 3×3 matrix.
pub(super) fn mat_det(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}
/// Inverse of a 3×3 matrix.  Returns identity if the determinant is too small.
pub(super) fn mat_inv(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = mat_det(a);
    if det.abs() < 1e-300 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let inv_det = 1.0 / det;
    [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
        ],
        [
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
        ],
    ]
}
/// Identity 3×3 matrix.
pub(super) const IDENTITY3: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// Add two 3×3 matrices.
pub(super) fn mat_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}
/// Scale a 3×3 matrix by a scalar.
pub(super) fn mat_scale(a: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}
/// Matrix-vector product y = A·x (3×3 · 3).
pub(super) fn mat_vec(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    [
        a[0][0] * x[0] + a[0][1] * x[1] + a[0][2] * x[2],
        a[1][0] * x[0] + a[1][1] * x[1] + a[1][2] * x[2],
        a[2][0] * x[0] + a[2][1] * x[1] + a[2][2] * x[2],
    ]
}
/// Cubic-spline SPH kernel W(r, h).
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
/// Gradient of W: ∇W = (dW/dr) · (rij / r).
///
/// Returns the vector gradient of the kernel evaluated at `rij = xi - xj`.
pub fn kernel_gradient(rij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(rij);
    if r < 1e-300 || h < 1e-300 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        -alpha / h * 0.75 * t * t
    } else {
        return [0.0; 3];
    };
    scale3(rij, dw_dr / r)
}
/// Compute the first Piola-Kirchhoff stress P and Cauchy stress σ for a particle.
///
/// Updates `particle.pk1_stress` and `particle.stress` in-place.
///
/// # Arguments
/// * `particle` – elastic particle whose `deformation_gradient` is current
/// * `model`    – constitutive model to use
pub fn compute_stress(particle: &mut ElasticParticle, model: &ConstitutiveModel) {
    let f = particle.deformation_gradient;
    let mu = particle.shear_modulus();
    let lambda = particle.lame_lambda();
    let j = mat_det(f).max(1e-15);
    match model {
        ConstitutiveModel::StVenantKirchhoff => {
            let ft = mat_transpose(f);
            let ftf = mat_mul(ft, f);
            let mut green = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    green[i][k] = 0.5 * (ftf[i][k] - IDENTITY3[i][k]);
                }
            }
            let tr_e = mat_trace(green);
            let mut s = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    s[i][k] = lambda * tr_e * IDENTITY3[i][k] + 2.0 * mu * green[i][k];
                }
            }
            particle.pk1_stress = mat_mul(f, s);
            let fs = mat_mul(f, s);
            let fsft = mat_mul(fs, ft);
            particle.stress = mat_scale(fsft, 1.0 / j);
        }
        ConstitutiveModel::NeoHookean => {
            let f_inv = mat_inv(f);
            let f_inv_t = mat_transpose(f_inv);
            let ln_j = j.ln();
            let mut p = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    p[i][k] = mu * (f[i][k] - f_inv_t[i][k]) + lambda * ln_j * f_inv_t[i][k];
                }
            }
            particle.pk1_stress = p;
            let pft = mat_mul(p, mat_transpose(f));
            particle.stress = mat_scale(pft, 1.0 / j);
        }
    }
}
/// Compute the SPH acceleration from the stress divergence for all particles.
///
/// Uses the symmetric SPH approximation:
///
/// ```text
/// aᵢ = Σⱼ mⱼ (Pᵢ/ρᵢ² + Pⱼ/ρⱼ²) · ∇Wᵢⱼ
/// ```
///
/// where `Pᵢ` is the first Piola-Kirchhoff stress tensor.
/// Results are written into `particle.acceleration`.
///
/// # Arguments
/// * `particles` – mutable slice of all elastic particles
pub fn stress_divergence(particles: &mut [ElasticParticle]) {
    let n = particles.len();
    let pos: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
    let mass: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let rho: Vec<f64> = particles.iter().map(|p| p.density).collect();
    let h_sph: Vec<f64> = particles.iter().map(|p| p.h).collect();
    let pk1: Vec<[[f64; 3]; 3]> = particles.iter().map(|p| p.pk1_stress).collect();
    let mut accel = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let pi_rho2 = mat_scale(pk1[i], 1.0 / (rho[i] * rho[i]));
        for j in 0..n {
            if i == j {
                continue;
            }
            let rij = sub3(pos[i], pos[j]);
            let h_avg = 0.5 * (h_sph[i] + h_sph[j]);
            let grad_w = kernel_gradient(rij, h_avg);
            if grad_w == [0.0; 3] {
                continue;
            }
            let pj_rho2 = mat_scale(pk1[j], 1.0 / (rho[j] * rho[j]));
            let sum_p = mat_add(pi_rho2, pj_rho2);
            let force = mat_vec(sum_p, grad_w);
            let mj = mass[j];
            for d in 0..3 {
                accel[i][d] -= mj * force[d];
            }
        }
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.acceleration = accel[i];
    }
}
/// Hourglass (zero-energy mode) control correction.
///
/// Applies an artificial hourglass viscosity following Randles & Libersky (1996):
///
/// ```text
/// aᵢ_hg = α_hg · cₛ · Σⱼ mⱼ/ρⱼ · (vᵢ−vⱼ) · W(rᵢⱼ, h)
/// ```
///
/// Adds the hourglass correction to `particle.acceleration`.
///
/// # Arguments
/// * `particles`   – mutable slice of elastic particles
/// * `alpha_hg`    – hourglass coefficient (typically 0.01–0.1)
/// * `sound_speed` – characteristic wave speed (m/s)
pub fn hourglass_control(particles: &mut [ElasticParticle], alpha_hg: f64, sound_speed: f64) {
    let n = particles.len();
    let pos: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
    let vel: Vec<[f64; 3]> = particles.iter().map(|p| p.velocity).collect();
    let mass: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let rho: Vec<f64> = particles.iter().map(|p| p.density).collect();
    let h_sph: Vec<f64> = particles.iter().map(|p| p.h).collect();
    let mut corrections = vec![[0.0f64; 3]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let rij = sub3(pos[i], pos[j]);
            let r = len3(rij);
            let h_avg = 0.5 * (h_sph[i] + h_sph[j]);
            let w = cubic_kernel(r, h_avg);
            if w == 0.0 {
                continue;
            }
            let dv = sub3(vel[i], vel[j]);
            let factor = alpha_hg * sound_speed * mass[j] / rho[j] * w;
            for d in 0..3 {
                corrections[i][d] += factor * dv[d];
            }
        }
    }
    for (i, p) in particles.iter_mut().enumerate() {
        for d in 0..3 {
            p.acceleration[d] += corrections[i][d];
        }
    }
}
/// P-wave (compressional) speed in a linear elastic medium.
///
/// ```text
/// c_P = √((λ + 2μ) / ρ)
/// ```
///
/// # Arguments
/// * `youngs_modulus` – E (Pa)
/// * `poissons_ratio` – ν (dimensionless)
/// * `density`        – ρ (kg/m³)
///
/// Returns c_P in m/s.
pub fn p_wave_speed(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
    if density <= 0.0 {
        return 0.0;
    }
    let nu = poissons_ratio;
    let e = youngs_modulus;
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    ((lambda + 2.0 * mu) / density).sqrt()
}
/// S-wave (shear) speed in a linear elastic medium.
///
/// ```text
/// c_S = √(μ / ρ)
/// ```
///
/// # Arguments
/// * `youngs_modulus` – E (Pa)
/// * `poissons_ratio` – ν (dimensionless)
/// * `density`        – ρ (kg/m³)
///
/// Returns c_S in m/s.
pub fn s_wave_speed(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
    if density <= 0.0 {
        return 0.0;
    }
    let mu = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
    (mu / density).sqrt()
}
/// Returns `(c_P, c_S)` — both wave speeds.
///
/// See [`p_wave_speed`] and [`s_wave_speed`] for details.
pub fn elastic_wave_speed(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> (f64, f64) {
    (
        p_wave_speed(youngs_modulus, poissons_ratio, density),
        s_wave_speed(youngs_modulus, poissons_ratio, density),
    )
}
/// Maximum principal stress fracture criterion.
///
/// Returns `true` if the maximum principal stress of `sigma` exceeds
/// `tensile_strength`.
///
/// The eigenvalues of the symmetric 3×3 Cauchy stress are estimated using
/// Gershgorin's circle theorem upper bound.  A full eigendecomposition is
/// avoided to keep the implementation allocation-free.
///
/// # Arguments
/// * `sigma`            – Cauchy stress tensor (3×3, Pa)
/// * `tensile_strength` – σ_c (Pa)
pub fn fracture_criterion(sigma: [[f64; 3]; 3], tensile_strength: f64) -> bool {
    for i in 0..3 {
        let mut off_sum = 0.0f64;
        for j in 0..3 {
            if j != i {
                off_sum += sigma[i][j].abs();
            }
        }
        if sigma[i][i] + off_sum > tensile_strength {
            return true;
        }
    }
    false
}
/// Apply fracture criterion to all particles and mark failures.
///
/// Sets `particle.fractured = true` for any particle whose maximum principal
/// stress exceeds `tensile_strength`.
///
/// # Arguments
/// * `particles`        – mutable slice of elastic particles
/// * `tensile_strength` – σ_c (Pa)
pub fn apply_fracture(particles: &mut [ElasticParticle], tensile_strength: f64) {
    for p in particles.iter_mut() {
        if fracture_criterion(p.stress, tensile_strength) {
            p.fractured = true;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_identity_deformation_gradient() {
        let p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        assert_eq!(p.deformation_gradient, IDENTITY3);
    }
    #[test]
    fn test_new_zero_stress() {
        let p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        for row in &p.stress {
            for v in row {
                assert_eq!(*v, 0.0);
            }
        }
    }
    #[test]
    fn test_new_not_fractured() {
        let p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        assert!(!p.fractured);
    }
    #[test]
    fn test_shear_modulus_steel() {
        let p = ElasticParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        let expected = 200e9 / (2.0 * 1.3);
        assert!((p.shear_modulus() - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_lame_lambda_positive() {
        let p = ElasticParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert!(p.lame_lambda() > 0.0);
    }
    #[test]
    fn test_mat_mul_identity() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let c = mat_mul(a, IDENTITY3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_mat_det_identity_is_one() {
        assert!((mat_det(IDENTITY3) - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_mat_inv_identity() {
        let inv = mat_inv(IDENTITY3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((inv[i][j] - IDENTITY3[i][j]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_mat_inv_roundtrip() {
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]];
        let ainv = mat_inv(a);
        let should_be_id = mat_mul(a, ainv);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (should_be_id[i][j] - IDENTITY3[i][j]).abs() < 1e-12,
                    "roundtrip failed at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_mat_transpose() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = mat_transpose(a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - a[j][i]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_stvk_identity_f_zero_stress() {
        let mut p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        compute_stress(&mut p, &ConstitutiveModel::StVenantKirchhoff);
        for row in &p.pk1_stress {
            for v in row {
                assert!(
                    v.abs() < 1e-6,
                    "StVK: expected zero PK1 stress with F=I: {v}"
                );
            }
        }
    }
    #[test]
    fn test_neohookean_identity_f_zero_stress() {
        let mut p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        compute_stress(&mut p, &ConstitutiveModel::NeoHookean);
        for row in &p.pk1_stress {
            for v in row {
                assert!(
                    v.abs() < 1e-6,
                    "NeoHookean: expected zero PK1 stress with F=I: {v}"
                );
            }
        }
    }
    #[test]
    fn test_stvk_uniaxial_stretch_nonzero_stress() {
        let mut p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        p.deformation_gradient[0][0] = 1.01;
        compute_stress(&mut p, &ConstitutiveModel::StVenantKirchhoff);
        assert!(
            p.stress[0][0].abs() > 1e3,
            "StVK: non-zero stress expected under stretch"
        );
    }
    #[test]
    fn test_neohookean_compression_nonzero_stress() {
        let mut p = ElasticParticle::new([0.0; 3], 1e-3, 2700.0, 0.01, 70e9, 0.33);
        for i in 0..3 {
            for k in 0..3 {
                p.deformation_gradient[i][k] = if i == k { 0.99 } else { 0.0 };
            }
        }
        compute_stress(&mut p, &ConstitutiveModel::NeoHookean);
        assert!(
            p.stress[0][0] < 0.0,
            "NeoHookean: should have compressive stress"
        );
    }
    #[test]
    fn test_p_wave_faster_than_s_wave() {
        let (cp, cs) = elastic_wave_speed(200e9, 0.3, 7800.0);
        assert!(cp > cs, "P-wave must be faster than S-wave");
    }
    #[test]
    fn test_wave_speed_steel_approximate() {
        let (cp, _cs) = elastic_wave_speed(200e9, 0.3, 7800.0);
        assert!(
            cp > 5000.0 && cp < 7000.0,
            "steel P-wave speed out of range: {cp}"
        );
    }
    #[test]
    fn test_wave_speed_zero_density() {
        let (cp, cs) = elastic_wave_speed(200e9, 0.3, 0.0);
        assert_eq!(cp, 0.0);
        assert_eq!(cs, 0.0);
    }
    #[test]
    fn test_s_wave_speed_formula() {
        let e = 200e9f64;
        let nu = 0.3f64;
        let rho = 7800.0f64;
        let mu = e / (2.0 * (1.0 + nu));
        let cs_expected = (mu / rho).sqrt();
        let cs = s_wave_speed(e, nu, rho);
        assert!((cs - cs_expected).abs() / cs_expected < 1e-12);
    }
    #[test]
    fn test_fracture_not_triggered_below_strength() {
        let sigma = [[1e6, 0.0, 0.0], [0.0, 1e6, 0.0], [0.0, 0.0, 1e6]];
        assert!(!fracture_criterion(sigma, 500e6));
    }
    #[test]
    fn test_fracture_triggered_above_strength() {
        let sigma = [[600e6, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        assert!(fracture_criterion(sigma, 500e6));
    }
    #[test]
    fn test_fracture_zero_stress() {
        let sigma = [[0.0; 3]; 3];
        assert!(!fracture_criterion(sigma, 500e6));
    }
    #[test]
    fn test_fracture_shear_contribution() {
        let sigma = [[300e6, 300e6, 0.0], [300e6, 0.0, 0.0], [0.0, 0.0, 0.0]];
        assert!(fracture_criterion(sigma, 500e6));
    }
    #[test]
    fn test_apply_fracture_marks_failed_particle() {
        let mut particles = vec![ElasticParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3,
        )];
        particles[0].stress[0][0] = 600e6;
        apply_fracture(&mut particles, 500e6);
        assert!(particles[0].fractured);
    }
    #[test]
    fn test_apply_fracture_no_mark_below_strength() {
        let mut particles = vec![ElasticParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3,
        )];
        particles[0].stress[0][0] = 100e6;
        apply_fracture(&mut particles, 500e6);
        assert!(!particles[0].fractured);
    }
    #[test]
    fn test_stress_divergence_two_particles_finite() {
        let mut particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
            ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ];
        stress_divergence(&mut particles);
        for p in &particles {
            assert!(p.acceleration.iter().all(|v| v.is_finite()));
        }
    }
    #[test]
    fn test_stress_divergence_antisymmetry() {
        let mut particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
            ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ];
        for p in &mut particles {
            p.pk1_stress[0][0] = 1e6;
        }
        stress_divergence(&mut particles);
        for d in 0..3 {
            assert!(
                (particles[0].acceleration[d] + particles[1].acceleration[d]).abs() < 1e-6,
                "momentum should be conserved in direction {d}"
            );
        }
    }
    #[test]
    fn test_hourglass_control_runs() {
        let mut particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
            ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ];
        particles[0].velocity = [1.0, 0.0, 0.0];
        hourglass_control(&mut particles, 0.05, 5000.0);
        assert!(particles[0].acceleration.iter().all(|v| v.is_finite()));
    }
    #[test]
    fn test_hourglass_equal_velocities_no_correction() {
        let mut particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
            ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ];
        particles[0].velocity = [2.0, 1.0, 0.5];
        particles[1].velocity = [2.0, 1.0, 0.5];
        hourglass_control(&mut particles, 0.05, 5000.0);
        for p in &particles {
            for v in &p.acceleration {
                assert!(v.abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_kernel_gradient_direction() {
        let rij = [0.03, 0.0, 0.0];
        let g = kernel_gradient(rij, 0.05);
        assert!(
            g[0] < 0.0,
            "gradient x should be negative (kernel falls with r): {}",
            g[0]
        );
        assert!(g[1].abs() < 1e-14);
        assert!(g[2].abs() < 1e-14);
    }
    #[test]
    fn test_kernel_gradient_zero_r() {
        let g = kernel_gradient([0.0; 3], 0.05);
        assert_eq!(g, [0.0; 3]);
    }
    #[test]
    fn test_kernel_gradient_beyond_support() {
        let g = kernel_gradient([1.0, 0.0, 0.0], 0.05);
        assert_eq!(g, [0.0; 3]);
    }
}
/// 3×3 matrix multiply: C = A · B (row-major).
pub fn mat_mul_3x3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Transpose a 3×3 matrix.
pub fn transpose_3x3(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}
/// Trace of a 3×3 tensor (sum of diagonal elements).
pub fn trace_tensor(t: &[[f64; 3]; 3]) -> f64 {
    t[0][0] + t[1][1] + t[2][2]
}
/// Frobenius norm of a 3×3 tensor: sqrt(Σ tᵢⱼ²).
pub fn frobenius_norm(t: &[[f64; 3]; 3]) -> f64 {
    let mut s = 0.0_f64;
    for row in t {
        for &v in row {
            s += v * v;
        }
    }
    s.sqrt()
}
/// Compute the Green-Lagrange strain tensor E = (Fᵀ F − I) / 2.
///
/// # Arguments
/// * `f` – deformation gradient (3×3)
pub fn green_lagrange_strain(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ft = transpose_3x3(f);
    let ftf = mat_mul_3x3(&ft, f);
    let mut e = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            e[i][j] = 0.5 * (ftf[i][j] - delta);
        }
    }
    e
}
/// Linear elastic (Cauchy) stress via Hooke's law:
/// σ = λ tr(ε) I + 2μ ε
///
/// where λ and μ are Lamé parameters derived from Young's modulus `e_mod`
/// and Poisson's ratio `nu`.
///
/// # Arguments
/// * `strain` – infinitesimal strain tensor ε (3×3)
/// * `e_mod`  – Young's modulus (Pa)
/// * `nu`     – Poisson's ratio
pub fn linear_elastic_stress(strain: &[[f64; 3]; 3], e_mod: f64, nu: f64) -> [[f64; 3]; 3] {
    let lambda = e_mod * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e_mod / (2.0 * (1.0 + nu));
    let tr = trace_tensor(strain);
    let mut sigma = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            sigma[i][j] = lambda * tr * delta + 2.0 * mu * strain[i][j];
        }
    }
    sigma
}
/// Von Mises stress: σ_vm = sqrt(3/2 · s:s) where s = σ − (tr σ / 3) I.
///
/// # Arguments
/// * `stress` – Cauchy stress tensor σ (3×3, Pa)
pub fn von_mises_stress(stress: &[[f64; 3]; 3]) -> f64 {
    let tr = trace_tensor(stress) / 3.0;
    let mut s = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            s[i][j] = stress[i][j] - tr * delta;
        }
    }
    let s_inner: f64 = s.iter().flatten().map(|v| v * v).sum();
    (1.5 * s_inner).sqrt()
}
/// Compute the SPH approximation of the deformation gradient for particle `i`.
///
/// F ≈ Σⱼ (xⱼ − xᵢ) ⊗ ∇W · (Xⱼ − Xᵢ)  (simplified version assuming unit kernel correction)
///
/// Returns the identity matrix when there are no neighbours within range,
/// consistent with an undeformed reference configuration.
///
/// # Arguments
/// * `particles` – slice of elastic particles
/// * `i`         – index of the particle of interest
/// * `kernel`    – SPH kernel function W(r, h) → f64
pub fn deformation_gradient_sph(
    particles: &[ElasticParticle],
    i: usize,
    kernel: impl Fn(f64, f64) -> f64,
) -> [[f64; 3]; 3] {
    if i >= particles.len() {
        return IDENTITY3;
    }
    let pi = &particles[i];
    let xi = pi.position;
    let xi_ref = pi.ref_position;
    let h = pi.h;
    let mut f = [[0.0f64; 3]; 3];
    let mut has_neighbour = false;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let xj = pj.position;
        let xj_ref = pj.ref_position;
        let rij_ref: [f64; 3] = [
            xj_ref[0] - xi_ref[0],
            xj_ref[1] - xi_ref[1],
            xj_ref[2] - xi_ref[2],
        ];
        let r_ref = (rij_ref[0].powi(2) + rij_ref[1].powi(2) + rij_ref[2].powi(2)).sqrt();
        let w = kernel(r_ref, h);
        if w == 0.0 {
            continue;
        }
        has_neighbour = true;
        let dx: [f64; 3] = [xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]];
        let vol_j = pj.mass / pj.density.max(1e-30);
        for row in 0..3 {
            for col in 0..3 {
                f[row][col] += dx[row] * rij_ref[col] * w * vol_j;
            }
        }
    }
    if !has_neighbour {
        return IDENTITY3;
    }
    f
}
/// Hourglass stabilisation correction for particle `i`.
///
/// Adds a penalty force proportional to the difference between the
/// SPH-interpolated displacement and the displacement predicted by the
/// local deformation gradient (Randles & Libersky 2000).
///
/// # Arguments
/// * `particles`  – slice of elastic particles
/// * `i`          – particle index
/// * `alpha_hg`   – hourglass coefficient (dimensionless, typically 0.01–0.1)
pub fn hourglass_correction(particles: &[ElasticParticle], i: usize, alpha_hg: f64) -> [f64; 3] {
    if i >= particles.len() {
        return [0.0; 3];
    }
    let pi = &particles[i];
    let mut corr = [0.0f64; 3];
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij: [f64; 3] = [
            pj.ref_position[0] - pi.ref_position[0],
            pj.ref_position[1] - pi.ref_position[1],
            pj.ref_position[2] - pi.ref_position[2],
        ];
        let r = (rij[0].powi(2) + rij[1].powi(2) + rij[2].powi(2)).sqrt();
        let w = cubic_kernel(r, pi.h);
        if w == 0.0 {
            continue;
        }
        let f = pi.deformation_gradient;
        let du_pred: [f64; 3] = [
            f[0][0] * rij[0] + f[0][1] * rij[1] + f[0][2] * rij[2] - rij[0],
            f[1][0] * rij[0] + f[1][1] * rij[1] + f[1][2] * rij[2] - rij[1],
            f[2][0] * rij[0] + f[2][1] * rij[1] + f[2][2] * rij[2] - rij[2],
        ];
        let dv: [f64; 3] = [
            pj.velocity[0] - pi.velocity[0],
            pj.velocity[1] - pi.velocity[1],
            pj.velocity[2] - pi.velocity[2],
        ];
        let vol_j = pj.mass / pj.density.max(1e-30);
        for d in 0..3 {
            corr[d] += alpha_hg * (du_pred[d] - dv[d]) * w * vol_j;
        }
    }
    corr
}
#[cfg(test)]
mod tests_elastic_solver {
    use super::*;
    use crate::elastic_sph::types::*;
    #[test]
    fn test_mat_mul_3x3_identity() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let c = mat_mul_3x3(&a, &IDENTITY3);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - a[i][j]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_mat_mul_3x3_zero() {
        let z = [[0.0f64; 3]; 3];
        let a = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let c = mat_mul_3x3(&a, &z);
        for row in c {
            for v in row {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn test_mat_mul_3x3_associative() {
        let a = [[1.0, 2.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]];
        let b = [[2.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 2.0]];
        let c = [[1.0, 1.0, 1.0], [0.0, 2.0, 0.0], [1.0, 0.0, 1.0]];
        let ab_c = mat_mul_3x3(&mat_mul_3x3(&a, &b), &c);
        let a_bc = mat_mul_3x3(&a, &mat_mul_3x3(&b, &c));
        for i in 0..3 {
            for j in 0..3 {
                assert!((ab_c[i][j] - a_bc[i][j]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_transpose_3x3_correct() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = transpose_3x3(&a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - a[j][i]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_transpose_3x3_symmetric_unchanged() {
        let s = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];
        let t = transpose_3x3(&s);
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - s[i][j]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_transpose_twice_is_identity() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let tt = transpose_3x3(&transpose_3x3(&a));
        for i in 0..3 {
            for j in 0..3 {
                assert!((tt[i][j] - a[i][j]).abs() < 1e-14);
            }
        }
    }
    #[test]
    fn test_trace_identity_is_3() {
        assert!((trace_tensor(&IDENTITY3) - 3.0).abs() < 1e-14);
    }
    #[test]
    fn test_trace_diagonal() {
        let d = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        assert!((trace_tensor(&d) - 9.0).abs() < 1e-14);
    }
    #[test]
    fn test_trace_zero_matrix() {
        let z = [[0.0f64; 3]; 3];
        assert_eq!(trace_tensor(&z), 0.0);
    }
    #[test]
    fn test_frobenius_identity_is_sqrt3() {
        let f = frobenius_norm(&IDENTITY3);
        assert!((f - 3.0_f64.sqrt()).abs() < 1e-14);
    }
    #[test]
    fn test_frobenius_zero_is_zero() {
        let z = [[0.0f64; 3]; 3];
        assert_eq!(frobenius_norm(&z), 0.0);
    }
    #[test]
    fn test_frobenius_nonneg() {
        let a = [[1.0, -2.0, 3.0], [-4.0, 5.0, -6.0], [7.0, -8.0, 9.0]];
        assert!(frobenius_norm(&a) > 0.0);
    }
    #[test]
    fn test_gl_strain_identity_f_is_zero() {
        let e = green_lagrange_strain(&IDENTITY3);
        for row in e {
            for v in row {
                assert!(v.abs() < 1e-14, "GL strain from F=I should be 0, got {v}");
            }
        }
    }
    #[test]
    fn test_gl_strain_symmetric() {
        let f = [[1.1, 0.05, 0.0], [0.02, 0.98, 0.0], [0.0, 0.0, 1.0]];
        let e = green_lagrange_strain(&f);
        for i in 0..3 {
            for j in 0..3 {
                assert!((e[i][j] - e[j][i]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_elastic_stress_zero_strain_zero_stress() {
        let eps = [[0.0f64; 3]; 3];
        let sigma = linear_elastic_stress(&eps, 200e9, 0.3);
        for row in sigma {
            for v in row {
                assert!(v.abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_elastic_stress_symmetric() {
        let eps = [
            [0.01, 0.005, 0.0],
            [0.005, 0.02, 0.001],
            [0.0, 0.001, 0.015],
        ];
        let sigma = linear_elastic_stress(&eps, 70e9, 0.33);
        for i in 0..3 {
            for j in 0..3 {
                assert!((sigma[i][j] - sigma[j][i]).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_elastic_stress_uniaxial() {
        let mut eps = [[0.0f64; 3]; 3];
        eps[0][0] = 0.001;
        let sigma = linear_elastic_stress(&eps, 200e9, 0.3);
        assert!(sigma[0][0] > 0.0);
    }
    #[test]
    fn test_von_mises_zero_stress() {
        let s = [[0.0f64; 3]; 3];
        assert!((von_mises_stress(&s)).abs() < 1e-10);
    }
    #[test]
    fn test_von_mises_isotropic_hydrostatic_zero() {
        let p = 1e6;
        let s = [[p, 0.0, 0.0], [0.0, p, 0.0], [0.0, 0.0, p]];
        assert!(von_mises_stress(&s) < 1e-6, "hydrostatic → von Mises = 0");
    }
    #[test]
    fn test_von_mises_pure_shear() {
        let tau = 1e6;
        let s = [[0.0, tau, 0.0], [tau, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let expected = (3.0f64).sqrt() * tau;
        let got = von_mises_stress(&s);
        assert!((got - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_von_mises_nonneg() {
        let s = [[1e6, 2e5, 0.0], [2e5, -5e5, 1e5], [0.0, 1e5, 3e5]];
        assert!(von_mises_stress(&s) >= 0.0);
    }
    #[test]
    fn test_elastic_solver_empty() {
        let solver = ElasticSolver::new(200e9, 0.3);
        assert_eq!(solver.particle_count(), 0);
    }
    #[test]
    fn test_elastic_solver_add_particle() {
        let mut solver = ElasticSolver::new(200e9, 0.3);
        let p = ElasticParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        solver.add_particle(p);
        assert_eq!(solver.particle_count(), 1);
    }
    #[test]
    fn test_elastic_solver_strain_energy_undeformed_zero() {
        let mut solver = ElasticSolver::new(200e9, 0.3);
        solver.add_particle(ElasticParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3,
        ));
        let se = solver.strain_energy();
        assert!(
            se.abs() < 1e-6,
            "undeformed strain energy should be 0: {se}"
        );
    }
    #[test]
    fn test_hourglass_correction_empty() {
        let particles: Vec<ElasticParticle> = Vec::new();
        let corr = hourglass_correction(&particles, 0, 0.05);
        assert_eq!(corr, [0.0; 3]);
    }
    #[test]
    fn test_hourglass_correction_single_particle() {
        let particles = vec![ElasticParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.05, 200e9, 0.3,
        )];
        let corr = hourglass_correction(&particles, 0, 0.05);
        for v in corr {
            assert_eq!(v, 0.0);
        }
    }
    #[test]
    fn test_hourglass_correction_finite() {
        let particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 7800.0, 0.05, 200e9, 0.3),
            ElasticParticle::new([0.02, 0.0, 0.0], 1e-3, 7800.0, 0.05, 200e9, 0.3),
        ];
        let corr = hourglass_correction(&particles, 0, 0.05);
        for v in corr {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn test_deformation_gradient_single_particle_is_zero() {
        let particles = vec![ElasticParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.05, 200e9, 0.3,
        )];
        let f = deformation_gradient_sph(&particles, 0, cubic_kernel);
        assert_eq!(f, IDENTITY3);
    }
}
#[cfg(test)]
mod tests_new {

    use crate::elastic_sph::types::*;
    #[test]
    fn test_elastic_sph_particle_new_zero_stress() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        for row in &p.stress {
            for v in row {
                assert_eq!(*v, 0.0);
            }
        }
    }
    #[test]
    fn test_elastic_sph_particle_shear_modulus() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        let expected = 200e9 / (2.0 * 1.3);
        assert!((p.shear_modulus() - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_elastic_sph_particle_bulk_modulus_positive() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert!(p.bulk_modulus() > 0.0);
    }
    #[test]
    fn test_elastic_sph_particle_pressure_zero_stress() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert!(p.pressure().abs() < 1e-10);
    }
    #[test]
    fn test_elastic_sph_particle_von_mises_zero() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert!(p.von_mises_eq().abs() < 1e-10);
    }
    #[test]
    fn test_elastic_sph_particle_not_fractured() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert!(!p.fractured);
    }
    #[test]
    fn test_elastic_sph_particle_plastic_strain_zero() {
        let p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        assert_eq!(p.plastic_strain, 0.0);
    }
    #[test]
    fn test_strain_rate_tensor_symmetric() {
        let raw = [[1.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 0.5]];
        let srt = StrainRateTensor::from_tensor(raw);
        let t = srt.tensor;
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (t[i][j] - t[j][i]).abs() < 1e-14,
                    "not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_strain_rate_tensor_volumetric_isotropic() {
        let a = 2.0;
        let raw = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]];
        let srt = StrainRateTensor::from_tensor(raw);
        assert!((srt.volumetric - 3.0 * a).abs() < 1e-14);
    }
    #[test]
    fn test_strain_rate_tensor_zero_input() {
        let srt = StrainRateTensor::from_tensor([[0.0; 3]; 3]);
        assert_eq!(srt.volumetric, 0.0);
        assert_eq!(srt.j2, 0.0);
    }
    #[test]
    fn test_strain_rate_tensor_j2_pure_shear() {
        let s = 0.1;
        let raw = [[0.0, s, 0.0], [s, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let srt = StrainRateTensor::from_tensor(raw);
        assert!((srt.j2 - s * s).abs() < 1e-14);
    }
    #[test]
    fn test_strain_rate_tensor_equivalent_positive() {
        let raw = [[0.01, 0.0, 0.0], [0.0, -0.005, 0.0], [0.0, 0.0, -0.005]];
        let srt = StrainRateTensor::from_tensor(raw);
        assert!(srt.equivalent() > 0.0);
    }
    #[test]
    fn test_strain_rate_tensor_compute_sph_empty() {
        let particles: Vec<ElasticSphParticle> = vec![];
        let srt = StrainRateTensor::compute_sph(&particles, 0);
        assert_eq!(srt.volumetric, 0.0);
    }
    #[test]
    fn test_strain_rate_tensor_compute_sph_single_particle() {
        let particles = vec![ElasticSphParticle::new(
            [0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3,
        )];
        let srt = StrainRateTensor::compute_sph(&particles, 0);
        assert_eq!(srt.volumetric, 0.0);
    }
    #[test]
    fn test_sph_elastic_stress_shear_modulus_from_params() {
        let ses = SphElasticStress::from_elastic_params(200e9, 0.3);
        let expected_g = 200e9 / 2.6;
        assert!((ses.shear_modulus - expected_g).abs() / expected_g < 1e-12);
    }
    #[test]
    fn test_sph_elastic_stress_bulk_modulus_from_params() {
        let ses = SphElasticStress::from_elastic_params(200e9, 0.3);
        let expected_k = 200e9 / (3.0 * 0.4);
        assert!((ses.bulk_modulus - expected_k).abs() / expected_k < 1e-12);
    }
    #[test]
    fn test_jaumann_increment_zero_rate_gives_zero() {
        let ses = SphElasticStress::from_elastic_params(200e9, 0.3);
        let sigma = [[1e6, 0.0, 0.0], [0.0, 1e6, 0.0], [0.0, 0.0, 1e6]];
        let ds = ses.jaumann_stress_increment(sigma, [[0.0; 3]; 3], [[0.0; 3]; 3], 1e-5);
        for row in &ds {
            for v in row {
                assert!(v.abs() < 1e-6, "zero strain rate → zero increment: {v}");
            }
        }
    }
    #[test]
    fn test_jaumann_increment_uniaxial_stretch_positive() {
        let ses = SphElasticStress::from_elastic_params(200e9, 0.3);
        let mut eps_dot = [[0.0f64; 3]; 3];
        eps_dot[0][0] = 0.001;
        let ds = ses.jaumann_stress_increment([[0.0; 3]; 3], eps_dot, [[0.0; 3]; 3], 1.0);
        assert!(
            ds[0][0] > 0.0,
            "uniaxial stretch → positive σ₁₁ increment: {}",
            ds[0][0]
        );
    }
    #[test]
    fn test_update_particle_stress_no_motion_unchanged() {
        let ses = SphElasticStress::from_elastic_params(200e9, 0.3);
        let mut p = ElasticSphParticle::new([0.0; 3], 1e-3, 7800.0, 0.01, 200e9, 0.3);
        p.stress[0][0] = 1e6;
        p.stress[1][1] = 1e6;
        p.stress[2][2] = 1e6;
        ses.update_particle_stress(&mut p, [[0.0; 3]; 3], 1e-5);
        assert!(
            (p.stress[0][0] - 1e6).abs() < 1.0,
            "stress unchanged: {}",
            p.stress[0][0]
        );
    }
    #[test]
    fn test_jc_yield_stress_at_zero_plastic_strain() {
        let jc = JohnsonCookModel::new_steel();
        let ys = jc.yield_stress(0.0, 1.0, 293.0);
        assert!((ys - jc.a).abs() / jc.a < 1e-12);
    }
    #[test]
    fn test_jc_yield_stress_increases_with_plastic_strain() {
        let jc = JohnsonCookModel::new_steel();
        let ys0 = jc.yield_stress(0.0, 1.0, 293.0);
        let ys1 = jc.yield_stress(0.1, 1.0, 293.0);
        assert!(ys1 > ys0, "yield stress must increase with plastic strain");
    }
    #[test]
    fn test_jc_yield_stress_decreases_with_temperature() {
        let jc = JohnsonCookModel::new_steel();
        let ys_cold = jc.yield_stress(0.1, 1.0, 293.0);
        let ys_hot = jc.yield_stress(0.1, 1.0, 1200.0);
        assert!(
            ys_hot < ys_cold,
            "yield stress must decrease with temperature"
        );
    }
    #[test]
    fn test_jc_yield_stress_at_melt_is_zero() {
        let jc = JohnsonCookModel::new_steel();
        let ys = jc.yield_stress(0.0, 1.0, jc.t_melt);
        assert!(ys.abs() < 1e-6, "yield stress at T_melt should be 0: {ys}");
    }
    #[test]
    fn test_jc_radial_return_no_yield_unchanged() {
        let jc = JohnsonCookModel::new_steel();
        let s_trial = [[1e5, 0.0, 0.0], [0.0, 1e5, 0.0], [0.0, 0.0, 0.0]];
        let g = 80e9_f64;
        let (s_ret, d_eps) = jc.radial_return(s_trial, 0.0, 1.0, 293.0, g);
        assert!(d_eps.abs() < 1e-20, "no plastic flow expected: {d_eps}");
        for i in 0..3 {
            for j in 0..3 {
                assert!((s_ret[i][j] - s_trial[i][j]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_jc_radial_return_plastic_flow_limits_stress() {
        let jc = JohnsonCookModel::new_steel();
        let big = 2e9;
        let s_trial = [
            [big, 0.0, 0.0],
            [0.0, -big / 2.0, 0.0],
            [0.0, 0.0, -big / 2.0],
        ];
        let g = 80e9_f64;
        let (s_ret, d_eps) = jc.radial_return(s_trial, 0.0, 1.0, 293.0, g);
        let s_ret_norm: f64 = s_ret.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
        let s_trial_norm: f64 = s_trial.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            s_ret_norm < s_trial_norm,
            "returned stress must be smaller than trial"
        );
        assert!(d_eps > 0.0, "plastic strain increment must be positive");
    }
    #[test]
    fn test_elastic_wave_speed_p_faster_than_s() {
        let (cp, cs) = ElasticWaveSpeed::both(200e9, 0.3, 7800.0);
        assert!(cp > cs, "P-wave must be faster than S-wave");
    }
    #[test]
    fn test_elastic_wave_speed_steel_p_wave_range() {
        let cp = ElasticWaveSpeed::p_wave(200e9, 0.3, 7800.0);
        assert!(
            cp > 5000.0 && cp < 7000.0,
            "steel P-wave out of range: {cp}"
        );
    }
    #[test]
    fn test_elastic_wave_speed_zero_density_gives_zero() {
        let cp = ElasticWaveSpeed::p_wave(200e9, 0.3, 0.0);
        assert_eq!(cp, 0.0);
    }
    #[test]
    fn test_elastic_wave_speed_courant_dt_positive() {
        let dt = ElasticWaveSpeed::courant_dt(0.01, 200e9, 0.3, 7800.0, 0.3);
        assert!(dt > 0.0 && dt.is_finite());
    }
    #[test]
    fn test_elastic_wave_speed_rayleigh_less_than_s() {
        let cr = ElasticWaveSpeed::rayleigh_wave(200e9, 0.3, 7800.0);
        let cs = ElasticWaveSpeed::s_wave(200e9, 0.3, 7800.0);
        assert!(cr < cs, "Rayleigh speed must be less than S-wave speed");
    }
    #[test]
    fn test_elastic_wave_speed_bulk_wave_positive() {
        let cb = ElasticWaveSpeed::bulk_wave(200e9, 0.3, 7800.0);
        assert!(cb > 0.0);
    }
    #[test]
    fn test_sph_fracture_damage_rate_below_threshold_zero() {
        let frac = SphFracture::new_granite();
        assert_eq!(frac.damage_rate(frac.sigma_threshold * 0.5), 0.0);
    }
    #[test]
    fn test_sph_fracture_damage_rate_above_threshold_positive() {
        let frac = SphFracture::new_granite();
        let rate = frac.damage_rate(frac.sigma_threshold * 2.0);
        assert!(
            rate > 0.0,
            "damage rate above threshold must be positive: {rate}"
        );
    }
    #[test]
    fn test_sph_fracture_advance_damage_clamped_to_one() {
        let frac = SphFracture::new_granite();
        let d = frac.advance_damage(0.99, frac.sigma_threshold * 10.0, 1e10);
        assert!((d - 1.0).abs() < 1e-14, "damage must be clamped to 1: {d}");
    }
    #[test]
    fn test_sph_fracture_advance_damage_increases() {
        let frac = SphFracture::new_granite();
        let d0 = 0.1;
        let d1 = frac.advance_damage(d0, frac.sigma_threshold * 2.0, 1e-6);
        assert!(d1 >= d0, "damage must not decrease: {d0} → {d1}");
    }
    #[test]
    fn test_sph_fracture_effective_stress_zero_damage() {
        let sigma = [[1e6, 0.0, 0.0], [0.0, 1e6, 0.0], [0.0, 0.0, 1e6]];
        let frac = SphFracture::new_granite();
        let eff = frac.effective_stress(sigma, 0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!((eff[i][j] - sigma[i][j]).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_sph_fracture_effective_stress_full_damage_zero() {
        let sigma = [[1e6, 0.0, 0.0], [0.0, 1e6, 0.0], [0.0, 0.0, 1e6]];
        let frac = SphFracture::new_granite();
        let eff = frac.effective_stress(sigma, 1.0);
        for row in &eff {
            for v in row {
                assert!(v.abs() < 1e-6, "fully damaged → zero effective stress: {v}");
            }
        }
    }
    #[test]
    fn test_sph_fracture_fragment_size_finite() {
        let frac = SphFracture::new_granite();
        let s = frac.fragment_size(1e6, 2700.0, 1e4);
        assert!(
            s.is_finite() && s > 0.0,
            "fragment size must be finite positive: {s}"
        );
    }
    #[test]
    fn test_sph_fracture_fragment_velocity_positive() {
        let frac = SphFracture::new_granite();
        let v = frac.fragment_velocity(2700.0, 1e4);
        assert!(v > 0.0, "fragment velocity must be positive: {v}");
    }
    #[test]
    fn test_sph_fracture_fragment_size_zero_strain_rate_infinite() {
        let frac = SphFracture::new_granite();
        let s = frac.fragment_size(1e6, 2700.0, 0.0);
        assert!(
            s.is_infinite(),
            "zero strain rate → infinite fragment size: {s}"
        );
    }
}
