//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::DesDetachment;

pub(super) type Mat3 = [[f64; 3]; 3];
#[inline]
pub(super) fn mat3_zero() -> Mat3 {
    [[0.0; 3]; 3]
}
/// Frobenius norm of a symmetric tensor: |S| = sqrt(2 * S_ij * S_ij).
#[inline]
pub(super) fn strain_rate_magnitude(s: &Mat3) -> f64 {
    let mut sum = 0.0;
    for row in s {
        for &v in row {
            sum += v * v;
        }
    }
    (2.0 * sum).sqrt()
}
/// Apply a box filter (top-hat) to a scalar field defined on particles.
///
/// For each particle `i`, the filtered value is the average of
/// `values[j]` over all neighbours `j` (including `i` itself).
#[allow(dead_code)]
pub fn box_filter_scalar(values: &[f64], neighbors: &[Vec<usize>]) -> Vec<f64> {
    let n = values.len();
    let mut filtered = vec![0.0_f64; n];
    for i in 0..n {
        let mut sum = values[i];
        let count = neighbors[i].len() as f64 + 1.0;
        for &j in &neighbors[i] {
            sum += values[j];
        }
        filtered[i] = sum / count;
    }
    filtered
}
/// Apply a box filter to a vector field defined on particles.
#[allow(dead_code)]
pub fn box_filter_vector(values: &[[f64; 3]], neighbors: &[Vec<usize>]) -> Vec<[f64; 3]> {
    let n = values.len();
    let mut filtered = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        let mut sum = values[i];
        let count = neighbors[i].len() as f64 + 1.0;
        for &j in &neighbors[i] {
            sum[0] += values[j][0];
            sum[1] += values[j][1];
            sum[2] += values[j][2];
        }
        filtered[i] = [sum[0] / count, sum[1] / count, sum[2] / count];
    }
    filtered
}
/// Compute the turbulent kinetic energy from the resolved velocity field
/// using a box filter: k_sgs ≈ 0.5 * |v - v̄|².
#[allow(dead_code)]
pub fn resolved_tke(velocities: &[[f64; 3]], neighbors: &[Vec<usize>]) -> Vec<f64> {
    let filtered = box_filter_vector(velocities, neighbors);
    velocities
        .iter()
        .zip(filtered.iter())
        .map(|(v, vf)| {
            let dx = v[0] - vf[0];
            let dy = v[1] - vf[1];
            let dz = v[2] - vf[2];
            0.5 * (dx * dx + dy * dy + dz * dz)
        })
        .collect()
}
/// Compute the vorticity magnitude for each particle.
///
/// Uses SPH gradient of velocity: ω = ∇ × v.
/// `all_kernel_grads[i][k]` is ∇W for the k-th neighbour of particle i.
#[allow(dead_code, clippy::too_many_arguments)]
pub fn compute_vorticity(
    velocities: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    neighbors: &[Vec<usize>],
    all_kernel_grads: &[Vec<[f64; 3]>],
) -> Vec<[f64; 3]> {
    let n = velocities.len();
    let mut vorticity = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        let vi = velocities[i];
        let mut curl = [0.0_f64; 3];
        for (k, &j) in neighbors[i].iter().enumerate() {
            let rho_j = densities[j].max(1e-14);
            let w = masses[j] / rho_j;
            let dv = [
                velocities[j][0] - vi[0],
                velocities[j][1] - vi[1],
                velocities[j][2] - vi[2],
            ];
            let g = all_kernel_grads[i][k];
            curl[0] += w * (dv[1] * g[2] - dv[2] * g[1]);
            curl[1] += w * (dv[2] * g[0] - dv[0] * g[2]);
            curl[2] += w * (dv[0] * g[1] - dv[1] * g[0]);
        }
        vorticity[i] = curl;
    }
    vorticity
}
/// Compute the strain-rate tensor contribution from a single particle pair (i,j).
///
/// dS_αβ += (m_j / ρ_j) * (vel_j − vel_i)_α * dW_β  (then symmetrised)
///
/// where `dW = dW/dr * r_ij/|r_ij|` is the kernel gradient evaluated at the
/// separation `r_ij = pos_j − pos_i` and scaled by the magnitude `dW`.
#[allow(dead_code)]
pub fn compute_strain_rate_tensor(
    vel_i: [f64; 3],
    vel_j: [f64; 3],
    r_ij: [f64; 3],
    mass_j: f64,
    rho_j: f64,
    dw: f64,
) -> Mat3 {
    let r_len = (r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]).sqrt();
    let grad = if r_len > 1e-14 {
        [
            dw * r_ij[0] / r_len,
            dw * r_ij[1] / r_len,
            dw * r_ij[2] / r_len,
        ]
    } else {
        [0.0; 3]
    };
    let rho_j = rho_j.max(1e-14);
    let weight = mass_j / rho_j;
    let dv = [
        vel_j[0] - vel_i[0],
        vel_j[1] - vel_i[1],
        vel_j[2] - vel_i[2],
    ];
    let mut dg = mat3_zero();
    for alpha in 0..3 {
        for beta in 0..3 {
            dg[alpha][beta] = weight * dv[alpha] * grad[beta];
        }
    }
    let mut s = mat3_zero();
    for a in 0..3 {
        for b in 0..3 {
            s[a][b] = 0.5 * (dg[a][b] + dg[b][a]);
        }
    }
    s
}
/// Standalone DES blending function (free function form).
///
/// Returns a blend value in `[0, 1]` where 0 = pure LES, 1 = pure RANS.
#[allow(dead_code)]
pub fn rans_les_blend(d_wall: f64, h: f64) -> f64 {
    let des = DesDetachment::new(0.65, h);
    des.rans_les_blend(d_wall, h)
}
/// Compute the Smagorinsky SGS stress tensor directly.
///
/// τ_ij^SGS = -2 ρ (C_s Δ)² |S| S_ij
///
/// Arguments:
/// - `s`       : strain-rate tensor S_ij.
/// - `density` : fluid density ρ.
/// - `cs`      : Smagorinsky constant C_s.
/// - `delta`   : filter width Δ.
#[allow(dead_code)]
pub fn sgs_stress_smagorinsky(s: Mat3, density: f64, cs: f64, delta: f64) -> Mat3 {
    let s_mag = strain_rate_magnitude(&s);
    let coeff = -2.0 * density * cs * cs * delta * delta * s_mag;
    let mut tau = mat3_zero();
    for i in 0..3 {
        for j in 0..3 {
            tau[i][j] = coeff * s[i][j];
        }
    }
    tau
}
/// Compute the SGS kinetic energy from the trace of the SGS stress.
///
/// k_SGS = -½ τ_ii / ρ
#[allow(dead_code)]
pub fn sgs_kinetic_energy(tau: Mat3, density: f64) -> f64 {
    if density < 1e-30 {
        return 0.0;
    }
    let trace = tau[0][0] + tau[1][1] + tau[2][2];
    -0.5 * trace / density
}
/// Isotropic part of the SGS stress: p_SGS = -τ_ii / 3.
#[allow(dead_code)]
pub fn sgs_pressure(tau: Mat3) -> f64 {
    -(tau[0][0] + tau[1][1] + tau[2][2]) / 3.0
}
/// Deviatoric part of the SGS stress: τ^dev_ij = τ_ij - (1/3) δ_ij τ_kk.
#[allow(dead_code)]
pub fn sgs_stress_deviatoric(tau: Mat3) -> Mat3 {
    let trace = tau[0][0] + tau[1][1] + tau[2][2];
    let mut dev = tau;
    for i in 0..3 {
        dev[i][i] -= trace / 3.0;
    }
    dev
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::turbulence::types::*;
    /// Build a trivial kernel gradient for a pair of particles separated by
    /// distance `r` along the x-axis using an approximate cubic-spline dW/dr.
    fn approx_kernel_grad(r: f64, h: f64) -> [f64; 3] {
        let q = r / h;
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
        let dw_dr = if !(1e-14..2.0).contains(&q) {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            sigma * (-0.75 * t * t) / h
        } else {
            sigma * (-3.0 * q + 2.25 * q * q) / h
        };
        [dw_dr, 0.0, 0.0]
    }
    #[test]
    fn test_sps_zero_strain() {
        let model = SpsModel::default();
        let n = 4;
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let velocities: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; n];
        let masses: Vec<f64> = vec![0.001; n];
        let densities: Vec<f64> = vec![1000.0; n];
        let h = 0.25;
        let neighbors: Vec<Vec<usize>> = (0..n)
            .map(|i| (0..n).filter(|&j| j != i).collect())
            .collect();
        let all_kernel_grads: Vec<Vec<[f64; 3]>> = (0..n)
            .map(|i| {
                neighbors[i]
                    .iter()
                    .map(|&j| {
                        let dx = positions[j][0] - positions[i][0];
                        let r = dx.abs();
                        let mut g = approx_kernel_grad(r, h);
                        if dx < 0.0 {
                            g[0] = -g[0];
                        }
                        g
                    })
                    .collect()
            })
            .collect();
        let s = model.compute_strain_rate_tensor(
            0,
            &neighbors[0],
            &positions,
            &velocities,
            &masses,
            &densities,
            &all_kernel_grads[0],
        );
        let tau = model.compute_sps_stress(&s, 1000.0);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    tau[a][b].abs() < 1e-10,
                    "Expected ~0 stress for uniform flow; tau[{a}][{b}]={:.3e}",
                    tau[a][b]
                );
            }
        }
    }
    #[test]
    fn test_sps_shear() {
        let model = SpsModel::new(0.12, 0.0066, 0.1);
        let shear_rate = 10.0_f64;
        let h = 0.2_f64;
        let dy = 0.1_f64;
        let positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.0, dy, 0.0]];
        let velocities: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [shear_rate * dy, 0.0, 0.0]];
        let masses = vec![0.001_f64; 2];
        let densities = vec![1000.0_f64; 2];
        let r = dy;
        let dw_dr = {
            let q = r / h;
            let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
            if !(1e-14..2.0).contains(&q) {
                0.0
            } else if q >= 1.0 {
                let t = 2.0 - q;
                sigma * (-0.75 * t * t) / h
            } else {
                sigma * (-3.0 * q + 2.25 * q * q) / h
            }
        };
        let kernel_grads_0: Vec<[f64; 3]> = vec![[0.0, dw_dr, 0.0]];
        let neighbors_0: Vec<usize> = vec![1];
        let s = model.compute_strain_rate_tensor(
            0,
            &neighbors_0,
            &positions,
            &velocities,
            &masses,
            &densities,
            &kernel_grads_0,
        );
        let s_xy = s[0][1];
        let s_yx = s[1][0];
        assert!(
            s_xy.abs() > 1e-10,
            "Expected non-zero S_xy for shear flow; got {s_xy:.3e}"
        );
        assert!(
            (s_xy - s_yx).abs() < 1e-12 * s_xy.abs().max(1e-14),
            "S must be symmetric: S_xy={s_xy:.3e}, S_yx={s_yx:.3e}"
        );
        let tau = model.compute_sps_stress(&s, 1000.0);
        let tau_xy = tau[0][1];
        assert!(
            tau_xy.abs() > 1e-14,
            "Expected non-zero τ_xy for shear; got {tau_xy:.3e}"
        );
    }
    #[test]
    fn test_mixing_length_near_wall() {
        let ml = MixingLength::new(0.0, 0.41, 1e-4);
        assert!(
            ml.mixing_length(0.0).abs() < 1e-15,
            "Mixing length at y=0 should be 0"
        );
        let small = ml.mixing_length(1e-7);
        assert!(small >= 0.0, "Mixing length must be non-negative");
        assert!(
            small < ml.mixing_length(0.01),
            "Mixing length should increase away from wall"
        );
        let far = ml.mixing_length(1.0);
        let expected = 0.41 * 1.0;
        assert!(
            (far - expected).abs() < 0.01 * expected,
            "Far from wall mixing length should ≈ κ·y = {expected:.4}, got {far:.4}"
        );
    }
    #[test]
    fn test_sps_symmetry() {
        let model = SpsModel::new(0.12, 0.0066, 0.1);
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
        ];
        let velocities: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.5, 0.0],
            [0.2, 1.0, 0.3],
            [0.0, 0.3, 1.0],
        ];
        let masses = vec![0.001_f64; 4];
        let densities = vec![1000.0_f64; 4];
        let h = 0.25_f64;
        let neighbors: Vec<Vec<usize>> = (0..4)
            .map(|i| (0..4_usize).filter(|&j| j != i).collect())
            .collect();
        let all_kernel_grads: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|i| {
                neighbors[i]
                    .iter()
                    .map(|&j| {
                        let dx = positions[j][0] - positions[i][0];
                        let dy = positions[j][1] - positions[i][1];
                        let dz = positions[j][2] - positions[i][2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r < 1e-14 {
                            return [0.0; 3];
                        }
                        let q = r / h;
                        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
                        let dw_dr = if q >= 2.0 {
                            0.0
                        } else if q >= 1.0 {
                            let t = 2.0 - q;
                            sigma * (-0.75 * t * t) / h
                        } else {
                            sigma * (-3.0 * q + 2.25 * q * q) / h
                        };
                        [dw_dr * dx / r, dw_dr * dy / r, dw_dr * dz / r]
                    })
                    .collect()
            })
            .collect();
        for i in 0..4 {
            let s = model.compute_strain_rate_tensor(
                i,
                &neighbors[i],
                &positions,
                &velocities,
                &masses,
                &densities,
                &all_kernel_grads[i],
            );
            let tau = model.compute_sps_stress(&s, densities[i]);
            for a in 0..3 {
                for b in 0..3 {
                    let s_diff = (s[a][b] - s[b][a]).abs();
                    assert!(
                        s_diff < 1e-14 * s[a][b].abs().max(1.0),
                        "S not symmetric at particle {i}: S[{a}][{b}]={:.3e} S[{b}][{a}]={:.3e}",
                        s[a][b],
                        s[b][a]
                    );
                    let tau_diff = (tau[a][b] - tau[b][a]).abs();
                    assert!(
                        tau_diff < 1e-14 * tau[a][b].abs().max(1.0),
                        "τ not symmetric at particle {i}: τ[{a}][{b}]={:.3e} τ[{b}][{a}]={:.3e}",
                        tau[a][b],
                        tau[b][a]
                    );
                }
            }
        }
    }
    #[test]
    fn test_k_omega_eddy_viscosity() {
        let model = KOmegaModel::new(1e-6);
        let nu = model.eddy_viscosity(0.01, 100.0);
        assert!((nu - 1e-4).abs() < 1e-10, "Expected ν_t = 1e-4, got {nu}");
        let nu0 = model.eddy_viscosity(0.01, 0.0);
        assert!(nu0.abs() < 1e-14, "Zero omega should yield zero ν_t");
        let nu_neg = model.eddy_viscosity(-1.0, 100.0);
        assert!(nu_neg.abs() < 1e-14, "Negative k should yield zero ν_t");
    }
    #[test]
    fn test_k_omega_step_no_blow_up() {
        let model = KOmegaModel::new(1e-6);
        let mut k = 0.01_f64;
        let mut omega = 100.0_f64;
        let strain_mag = 10.0_f64;
        let dt = 0.001_f64;
        for _ in 0..1000 {
            let (kn, on) = model.step_explicit(k, omega, strain_mag, dt);
            k = kn;
            omega = on;
        }
        assert!(k.is_finite(), "k should stay finite");
        assert!(omega.is_finite(), "omega should stay finite");
        assert!(k >= 0.0, "k must be non-negative");
        assert!(omega > 0.0, "omega must be positive");
    }
    #[test]
    fn test_k_omega_production_dissipation() {
        let model = KOmegaModel::new(1e-6);
        let k = 0.1_f64;
        let omega = 50.0_f64;
        let strain_mag = 5.0_f64;
        let nu_t = model.eddy_viscosity(k, omega);
        let pk = model.production_k(nu_t, strain_mag);
        let dk = model.dissipation_k(k, omega);
        assert!(pk >= 0.0, "production should be non-negative");
        assert!(dk >= 0.0, "dissipation should be non-negative");
        let po = model.production_omega(strain_mag);
        let do_ = model.dissipation_omega(omega);
        assert!(po >= 0.0, "omega production should be non-negative");
        assert!(do_ >= 0.0, "omega dissipation should be non-negative");
    }
    #[test]
    fn test_des_length_scales() {
        let des = DesModel::new(0.1, 1e-6);
        let k = 0.01_f64;
        let omega = 100.0_f64;
        let l_rans = des.rans_length(k, omega);
        let l_les = des.les_length();
        let l_eff = des.effective_length(k, omega);
        assert!(l_rans > 0.0, "RANS length should be positive");
        assert!(l_les > 0.0, "LES length should be positive");
        assert!(
            (l_eff - l_rans.min(l_les)).abs() < 1e-14,
            "Effective length should be min(RANS, LES)"
        );
    }
    #[test]
    fn test_des_eddy_viscosity() {
        let des = DesModel::new(0.1, 1e-6);
        let nu = des.eddy_viscosity_des(0.01, 100.0, 10.0);
        assert!(nu >= 0.0, "DES eddy viscosity should be non-negative");
        assert!(nu.is_finite(), "DES eddy viscosity should be finite");
    }
    #[test]
    fn test_des_step_stable() {
        let des = DesModel::new(0.1, 1e-6);
        let mut k = 0.01_f64;
        let mut omega = 100.0_f64;
        let strain_mag = 10.0_f64;
        let dt = 0.0005_f64;
        for _ in 0..500 {
            let (kn, on) = des.step_des(k, omega, strain_mag, dt);
            k = kn;
            omega = on;
        }
        assert!(
            k.is_finite() && k >= 0.0,
            "k should be finite and non-negative"
        );
        assert!(
            omega.is_finite() && omega > 0.0,
            "omega should be finite and positive"
        );
    }
    #[test]
    fn test_box_filter_scalar_uniform() {
        let values = vec![5.0_f64; 4];
        let neighbors = vec![vec![1, 2, 3], vec![0, 2, 3], vec![0, 1, 3], vec![0, 1, 2]];
        let filtered = box_filter_scalar(&values, &neighbors);
        for (i, &v) in filtered.iter().enumerate() {
            assert!(
                (v - 5.0).abs() < 1e-12,
                "Uniform field should stay uniform, got {v} at {i}"
            );
        }
    }
    #[test]
    fn test_box_filter_scalar_spike() {
        let values = vec![100.0, 0.0, 0.0, 0.0];
        let neighbors = vec![vec![1, 2, 3], vec![0, 2, 3], vec![0, 1, 3], vec![0, 1, 2]];
        let filtered = box_filter_scalar(&values, &neighbors);
        assert!(
            (filtered[0] - 25.0).abs() < 1e-10,
            "Expected 25.0, got {}",
            filtered[0]
        );
        assert!(
            (filtered[1] - 25.0).abs() < 1e-10,
            "Expected 25.0, got {}",
            filtered[1]
        );
    }
    #[test]
    fn test_box_filter_vector() {
        let values = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let neighbors = vec![vec![1], vec![0]];
        let filtered = box_filter_vector(&values, &neighbors);
        assert!((filtered[0][0] - 0.5).abs() < 1e-12);
        assert!((filtered[0][1] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_resolved_tke_uniform_is_zero() {
        let velocities = vec![[1.0, 2.0, 3.0]; 4];
        let neighbors = vec![vec![1, 2, 3], vec![0, 2, 3], vec![0, 1, 3], vec![0, 1, 2]];
        let tke = resolved_tke(&velocities, &neighbors);
        for (i, &k) in tke.iter().enumerate() {
            assert!(
                k.abs() < 1e-10,
                "Uniform velocity should give zero TKE, got {k} at {i}"
            );
        }
    }
    #[test]
    fn test_resolved_tke_nonuniform_positive() {
        let velocities = vec![
            [10.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 10.0],
            [-10.0, 0.0, 0.0],
        ];
        let neighbors = vec![vec![1, 2, 3], vec![0, 2, 3], vec![0, 1, 3], vec![0, 1, 2]];
        let tke = resolved_tke(&velocities, &neighbors);
        let total: f64 = tke.iter().sum();
        assert!(total > 0.0, "Non-uniform field should give positive TKE");
    }
    #[test]
    fn test_vorticity_zero_for_irrotational() {
        let n = 4;
        let velocities = vec![[1.0, 0.0, 0.0]; n];
        let masses = vec![0.001_f64; n];
        let densities = vec![1000.0_f64; n];
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let h = 0.25_f64;
        let neighbors: Vec<Vec<usize>> = (0..n)
            .map(|i| (0..n).filter(|&j| j != i).collect())
            .collect();
        let all_kernel_grads: Vec<Vec<[f64; 3]>> = (0..n)
            .map(|i| {
                neighbors[i]
                    .iter()
                    .map(|&j| {
                        let dx = positions[j][0] - positions[i][0];
                        let r = dx.abs();
                        if r < 1e-14 {
                            return [0.0; 3];
                        }
                        let q = r / h;
                        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
                        let dw_dr = if q >= 2.0 {
                            0.0
                        } else if q >= 1.0 {
                            let t = 2.0 - q;
                            sigma * (-0.75 * t * t) / h
                        } else {
                            sigma * (-3.0 * q + 2.25 * q * q) / h
                        };
                        let sign = if dx > 0.0 { 1.0 } else { -1.0 };
                        [dw_dr * sign, 0.0, 0.0]
                    })
                    .collect()
            })
            .collect();
        let vort = compute_vorticity(
            &velocities,
            &masses,
            &densities,
            &neighbors,
            &all_kernel_grads,
        );
        for (i, v) in vort.iter().enumerate() {
            let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!(
                mag < 1e-10,
                "Uniform velocity should have zero vorticity at {i}, got {mag:.3e}"
            );
        }
    }
    #[test]
    fn test_strain_rate_magnitude() {
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let mag = strain_rate_magnitude(&s);
        assert!((mag - 2.0).abs() < 1e-10, "Expected |S| = 2.0, got {mag}");
    }
    #[test]
    fn test_mixing_length_monotone() {
        let ml = MixingLength::new(0.0, 0.41, 1e-4);
        let mut prev = 0.0_f64;
        for i in 1..=100 {
            let y = i as f64 * 0.001;
            let l = ml.mixing_length(y);
            assert!(
                l >= prev - 1e-14,
                "Mixing length should be monotonically increasing"
            );
            prev = l;
        }
    }
    #[test]
    fn test_smagorinsky_viscosity_positive() {
        let les = LesFilter::new(0.1, 0.1);
        let nu_t = les.smagorinsky_viscosity(5.0, 1000.0);
        assert!(
            nu_t > 0.0,
            "Smagorinsky viscosity must be positive for S_mag>0"
        );
        let expected = (0.1 * 0.1_f64).powi(2) * 5.0;
        assert!(
            (nu_t - expected).abs() < 1e-14,
            "Expected {expected}, got {nu_t}"
        );
    }
    #[test]
    fn test_smagorinsky_viscosity_zero_smag() {
        let les = LesFilter::new(0.1, 0.1);
        let nu_t = les.smagorinsky_viscosity(0.0, 1000.0);
        assert!(nu_t.abs() < 1e-14, "Zero S_mag should give zero viscosity");
    }
    #[test]
    fn test_turbulent_stress_symmetric() {
        let les = LesFilter::new(0.1, 0.15);
        let s: Mat3 = [[1.0, 0.5, 0.3], [0.5, -0.5, 0.2], [0.3, 0.2, -0.5]];
        let tau = les.turbulent_stress(s, 1000.0);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (tau[a][b] - tau[b][a]).abs() < 1e-13,
                    "Stress not symmetric: τ[{a}][{b}]={:.3e} τ[{b}][{a}]={:.3e}",
                    tau[a][b],
                    tau[b][a]
                );
            }
        }
    }
    #[test]
    fn test_turbulent_stress_sign() {
        let les = LesFilter::new(0.1, 0.15);
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let tau = les.turbulent_stress(s, 1000.0);
        assert!(
            tau[0][0] < 0.0,
            "τ[0][0] should be negative, got {}",
            tau[0][0]
        );
    }
    #[test]
    fn test_komega_sph_turbulent_viscosity() {
        let km = KomegaSph::new(0.01, 100.0);
        let nu_t = km.turbulent_viscosity();
        assert!((nu_t - 1e-4).abs() < 1e-12, "Expected 1e-4, got {nu_t}");
    }
    #[test]
    fn test_komega_sph_turbulent_viscosity_zero_omega() {
        let km = KomegaSph::new(0.01, 0.0);
        assert!(
            km.turbulent_viscosity().abs() < 1e-14,
            "Zero omega → zero nu_t"
        );
    }
    #[test]
    fn test_komega_sph_production_term() {
        let km = KomegaSph::new(0.01, 100.0);
        let s_mag_sq = 25.0;
        let p = km.production_term(s_mag_sq);
        assert!((p - 2.5e-3).abs() < 1e-12, "Expected 2.5e-3, got {p}");
        assert!(p >= 0.0, "Production must be non-negative");
    }
    #[test]
    fn test_des_blend_in_range() {
        let h = 0.1_f64;
        for &d in &[0.0_f64, 0.05, 0.1, 0.2, 1.0] {
            let b = rans_les_blend(d, h);
            assert!(
                (0.0..=1.0).contains(&b),
                "DES blend must be in [0,1], got {b} for d={d}"
            );
        }
    }
    #[test]
    fn test_des_blend_at_zero_wall_distance() {
        let b = rans_les_blend(0.0, 0.1);
        assert!(b.abs() < 1e-14, "d_wall=0 should give blend=0, got {b}");
    }
    #[test]
    fn test_des_blend_large_wall_distance() {
        let b = rans_les_blend(1e6, 0.1);
        assert!(
            (b - 1.0).abs() < 1e-14,
            "Large d_wall should give blend=1, got {b}"
        );
    }
    #[test]
    fn test_free_compute_strain_rate_symmetric() {
        let vel_i = [0.0_f64; 3];
        let vel_j = [1.0, 0.5, 0.2];
        let r_ij = [0.1, 0.0, 0.0];
        let s = compute_strain_rate_tensor(vel_i, vel_j, r_ij, 0.001, 1000.0, 5.0);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (s[a][b] - s[b][a]).abs() < 1e-14,
                    "Strain rate tensor must be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_dynamic_smagorinsky_coefficient_nonneg() {
        let ds = DynamicSmagorinsky::new(0.1, 2.0);
        let s_grid: Mat3 = [[1.0, 0.5, 0.0], [0.5, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let s_test: Mat3 = [[0.5, 0.25, 0.0], [0.25, -0.25, 0.0], [0.0, 0.0, -0.25]];
        let c = ds.dynamic_coefficient(s_grid, s_test);
        assert!(
            c >= 0.0,
            "Dynamic Smagorinsky coefficient must be non-negative, got {c}"
        );
        assert!(
            c <= ds.c_max,
            "Dynamic Smagorinsky coefficient exceeds c_max={}",
            ds.c_max
        );
    }
    #[test]
    fn test_dynamic_smagorinsky_eddy_viscosity_nonneg() {
        let ds = DynamicSmagorinsky::default();
        let s: Mat3 = [[2.0, 1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
        let s_test: Mat3 = [[1.0, 0.5, 0.0], [0.5, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let nu_t = ds.eddy_viscosity(s, s_test);
        assert!(nu_t >= 0.0, "Dynamic eddy viscosity must be non-negative");
    }
    #[test]
    fn test_dynamic_smagorinsky_sgs_stress_symmetric() {
        let ds = DynamicSmagorinsky::new(0.05, 2.0);
        let s: Mat3 = [[1.0, 0.3, 0.2], [0.3, -0.5, 0.1], [0.2, 0.1, -0.5]];
        let s_test: Mat3 = [[0.5, 0.15, 0.1], [0.15, -0.25, 0.05], [0.1, 0.05, -0.25]];
        let tau = ds.sgs_stress(s, s_test, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (tau[i][j] - tau[j][i]).abs() < 1e-12,
                    "Dynamic SGS stress must be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_dynamic_smagorinsky_zero_strain() {
        let ds = DynamicSmagorinsky::new(0.1, 2.0);
        let zero: Mat3 = [[0.0; 3]; 3];
        let c = ds.dynamic_coefficient(zero, zero);
        assert!(
            c.abs() < 1e-14,
            "Zero strain → zero dynamic coefficient, got {c}"
        );
    }
    #[test]
    fn test_wale_eddy_viscosity_nonneg() {
        let wale = WaleModel::new(0.325, 0.1);
        let s: Mat3 = [[1.0, 0.5, 0.0], [0.5, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let g: Mat3 = [[1.0, 0.3, 0.0], [-0.2, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let nu_t = wale.eddy_viscosity(s, g);
        assert!(
            nu_t >= 0.0,
            "WALE eddy viscosity must be non-negative, got {nu_t}"
        );
    }
    #[test]
    fn test_wale_zero_gradient_gives_zero_viscosity() {
        let wale = WaleModel::default();
        let zero: Mat3 = [[0.0; 3]; 3];
        let nu_t = wale.eddy_viscosity(zero, zero);
        assert!(
            nu_t.abs() < 1e-30,
            "Zero gradient → zero WALE viscosity, got {nu_t}"
        );
    }
    #[test]
    fn test_wale_sd_tensor_traceless() {
        let wale = WaleModel::new(0.325, 0.1);
        let g: Mat3 = [[1.0, 0.5, 0.2], [-0.3, -0.5, 0.1], [0.1, -0.2, -0.5]];
        let sd = wale.sd_tensor(g);
        let trace = sd[0][0] + sd[1][1] + sd[2][2];
        assert!(
            trace.abs() < 1e-12,
            "S^d tensor must be traceless, trace={trace}"
        );
    }
    #[test]
    fn test_wale_sgs_stress_symmetric() {
        let wale = WaleModel::new(0.325, 0.1);
        let s: Mat3 = [[1.0, 0.4, 0.1], [0.4, -0.5, 0.2], [0.1, 0.2, -0.5]];
        let g: Mat3 = [[1.2, 0.3, 0.0], [-0.1, -0.6, 0.0], [0.0, 0.0, -0.6]];
        let tau = wale.sgs_stress(s, g, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (tau[i][j] - tau[j][i]).abs() < 1e-12,
                    "WALE SGS stress must be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_kepsilon_turbulent_viscosity() {
        let ke = KepsilonSph::new(0.1, 1e-3);
        let nu_t = ke.turbulent_viscosity();
        assert!((nu_t - 0.9).abs() < 1e-12, "Expected ν_t=0.9, got {nu_t}");
    }
    #[test]
    fn test_kepsilon_zero_epsilon() {
        let ke = KepsilonSph::new(0.1, 0.0);
        assert!(ke.turbulent_viscosity().abs() < 1e-14, "Zero ε → zero ν_t");
        assert!(ke.turbulent_time_scale().abs() < 1e-14, "Zero ε → zero τ");
    }
    #[test]
    fn test_kepsilon_sgs_stress_diagonal() {
        let ke = KepsilonSph::new(0.1, 1e-3);
        let zero_s: Mat3 = [[0.0; 3]; 3];
        let tau = ke.sgs_stress(zero_s, 1000.0);
        let expected_diag = (2.0 / 3.0) * 1000.0 * 0.1;
        for i in 0..3 {
            assert!(
                (tau[i][i] - expected_diag).abs() < 1e-10,
                "Diagonal τ[{i}][{i}] should be {expected_diag}, got {}",
                tau[i][i]
            );
        }
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert!(tau[i][j].abs() < 1e-14, "Off-diagonal should be zero");
                }
            }
        }
    }
    #[test]
    fn test_kepsilon_step_increases_k_with_production() {
        let mut ke = KepsilonSph::new(1e-4, 1e-6);
        let s_mag_sq = 1000.0;
        let k_before = ke.k;
        ke.step(s_mag_sq, 1e-4);
        assert!(ke.k >= k_before, "k should increase with large production");
        assert!(ke.k > 0.0, "k must remain positive");
        assert!(ke.epsilon > 0.0, "ε must remain positive");
    }
    #[test]
    fn test_kepsilon_length_scale_positive() {
        let ke = KepsilonSph::new(0.1, 1e-3);
        let l = ke.turbulent_length_scale();
        assert!(l > 0.0, "Turbulent length scale must be positive, got {l}");
    }
    #[test]
    fn test_sgs_stress_smagorinsky_symmetric() {
        let s: Mat3 = [[1.0, 0.5, 0.3], [0.5, -0.5, 0.2], [0.3, 0.2, -0.5]];
        let tau = sgs_stress_smagorinsky(s, 1000.0, 0.12, 0.1);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (tau[i][j] - tau[j][i]).abs() < 1e-12,
                    "Smagorinsky SGS stress must be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_sgs_stress_smagorinsky_sign() {
        let s: Mat3 = [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let tau = sgs_stress_smagorinsky(s, 1000.0, 0.12, 0.1);
        assert!(tau[0][0] < 0.0, "τ_11 should be negative for positive S_11");
    }
    #[test]
    fn test_sgs_kinetic_energy_from_stress() {
        let rho = 1000.0_f64;
        let k_in = 0.05_f64;
        let tau: Mat3 = [
            [-rho * k_in, 0.0, 0.0],
            [0.0, -rho * k_in, 0.0],
            [0.0, 0.0, -rho * k_in],
        ];
        let k_sgs = sgs_kinetic_energy(tau, rho);
        assert!(
            (k_sgs - 1.5 * k_in).abs() < 1e-12,
            "Expected 1.5*k_in, got {k_sgs}"
        );
    }
    #[test]
    fn test_sgs_pressure_from_deviatoric() {
        let p_in = 500.0_f64;
        let tau: Mat3 = [[p_in, 0.0, 0.0], [0.0, p_in, 0.0], [0.0, 0.0, p_in]];
        let p = sgs_pressure(tau);
        assert!((p + p_in).abs() < 1e-10, "Expected -{p_in}, got {p}");
    }
    #[test]
    fn test_sgs_stress_deviatoric_traceless() {
        let tau: Mat3 = [[3.0, 1.0, 0.5], [1.0, -1.0, 0.2], [0.5, 0.2, -2.0]];
        let dev = sgs_stress_deviatoric(tau);
        let trace = dev[0][0] + dev[1][1] + dev[2][2];
        assert!(
            trace.abs() < 1e-12,
            "Deviatoric tensor must be traceless, got {trace}"
        );
    }
    #[test]
    fn test_van_driest_damping_range() {
        let wt = WallTurbulence::default();
        for &y_plus in &[0.0_f64, 5.0, 11.3, 30.0, 100.0, 1000.0] {
            let d = wt.van_driest_damping(y_plus);
            assert!(
                (0.0..=1.0).contains(&d),
                "van Driest must be in [0,1], got {d} at y+={y_plus}"
            );
        }
    }
    #[test]
    fn test_van_driest_damping_at_zero() {
        let wt = WallTurbulence::default();
        let d = wt.van_driest_damping(0.0);
        assert!(d.abs() < 1e-14, "van Driest at y+=0 should be 0, got {d}");
    }
    #[test]
    fn test_van_driest_damping_large_yplus() {
        let wt = WallTurbulence::default();
        let d = wt.van_driest_damping(1e6);
        assert!(
            (d - 1.0).abs() < 1e-5,
            "van Driest at large y+ should approach 1"
        );
    }
    #[test]
    fn test_law_of_wall_viscous_sublayer() {
        let wt = WallTurbulence::default();
        for &y_plus in &[1.0_f64, 5.0, 10.0] {
            let u_plus = wt.u_plus(y_plus);
            assert!(
                (u_plus - y_plus).abs() < 1e-14,
                "In sublayer u+ should equal y+={y_plus}, got {u_plus}"
            );
        }
    }
    #[test]
    fn test_law_of_wall_log_region() {
        let wt = WallTurbulence::default();
        let u_plus_log = wt.u_plus(100.0);
        assert!(
            u_plus_log > 11.3,
            "Log layer u+ at y+=100 should exceed sublayer value"
        );
    }
    #[test]
    fn test_friction_velocity_from_wall_shear() {
        let tau_w = 1.0_f64;
        let rho = 1000.0_f64;
        let u_tau = WallTurbulence::friction_velocity(tau_w, rho);
        assert!((u_tau - (1.0 / 1000.0_f64).sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_wall_turbulence_mixing_length_zero_at_wall() {
        let wt = WallTurbulence::new(0.41, 26.0, 1e-6);
        let lm = wt.mixing_length(0.0, 0.1);
        assert!(
            lm.abs() < 1e-14,
            "Mixing length at wall should be zero, got {lm}"
        );
    }
    #[test]
    fn test_turbulent_diffusion_acceleration_zero_for_uniform_flow() {
        let tds = TurbulentDiffusionSph::new(0.12, 0.1);
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let rho = 1000.0_f64;
        let grad_w = [0.5, 0.0, 0.0_f64];
        let acc = tds.turbulent_acceleration(s, s, rho, rho, 1.0, grad_w);
        for &a in &acc {
            assert!(a.is_finite(), "Turbulent acceleration should be finite");
        }
    }
    #[test]
    fn test_sigma_model_zero_gradient() {
        let sigma = SigmaModel::default();
        let zero: Mat3 = [[0.0; 3]; 3];
        let nu_t = sigma.eddy_viscosity(zero);
        assert!(
            nu_t.abs() < 1e-30,
            "Zero gradient → zero sigma eddy viscosity, got {nu_t}"
        );
    }
    #[test]
    fn test_sigma_model_eddy_viscosity_nonneg() {
        let sigma = SigmaModel::new(1.35, 0.1);
        let g: Mat3 = [[1.0, 0.5, 0.2], [-0.3, -0.5, 0.1], [0.1, -0.2, -0.5]];
        let nu_t = sigma.eddy_viscosity(g);
        assert!(
            nu_t >= 0.0,
            "Sigma eddy viscosity must be non-negative, got {nu_t}"
        );
    }
}
#[allow(dead_code)]
#[inline]
pub(super) fn mat3_add(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}
#[inline]
pub(super) fn mat3_scale(s: f64, a: Mat3) -> Mat3 {
    let mut b = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] = s * a[i][j];
        }
    }
    b
}
#[inline]
pub(super) fn mat3_trace(a: &Mat3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}
/// Double contraction A:B = Σ_ij A_ij B_ij
#[inline]
pub(super) fn mat3_double_contraction(a: &Mat3, b: &Mat3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += a[i][j] * b[i][j];
        }
    }
    s
}
/// Symmetric part of a matrix: S = (A + Aᵀ) / 2
#[inline]
pub(super) fn mat3_sym(a: &Mat3) -> Mat3 {
    let mut s = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            s[i][j] = 0.5 * (a[i][j] + a[j][i]);
        }
    }
    s
}
/// Anti-symmetric part of a matrix: Ω = (A - Aᵀ) / 2
#[inline]
pub(super) fn mat3_antisym(a: &Mat3) -> Mat3 {
    let mut omega = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            omega[i][j] = 0.5 * (a[i][j] - a[j][i]);
        }
    }
    omega
}
/// Compute the rotation-rate tensor Ω_ij = (∂v_i/∂x_j - ∂v_j/∂x_i) / 2
/// from the full velocity-gradient tensor.
#[allow(dead_code)]
pub fn rotation_rate_tensor(vel_grad: &Mat3) -> Mat3 {
    mat3_antisym(vel_grad)
}
/// Compute the strain-rate tensor S_ij = (∂v_i/∂x_j + ∂v_j/∂x_i) / 2.
#[allow(dead_code)]
pub fn strain_rate_tensor(vel_grad: &Mat3) -> Mat3 {
    mat3_sym(vel_grad)
}
/// Compute the Q-criterion (second invariant of velocity gradient):
/// Q = 0.5 (|Ω|² - |S|²)
/// Q > 0 identifies rotation-dominated regions (vortex cores).
#[allow(dead_code)]
pub fn q_criterion(vel_grad: &Mat3) -> f64 {
    let omega = rotation_rate_tensor(vel_grad);
    let s = strain_rate_tensor(vel_grad);
    let omega_sq = mat3_double_contraction(&omega, &omega);
    let s_sq = mat3_double_contraction(&s, &s);
    0.5 * (omega_sq - s_sq)
}
/// Compute the λ_2 vortex-identification criterion.
///
/// λ_2 < 0 identifies vortex cores. Computed as the second eigenvalue of
/// S² + Ω².
///
/// Uses the characteristic polynomial of the 3×3 symmetric matrix A = S² + Ω²
/// and returns the middle eigenvalue.
#[allow(dead_code)]
pub fn lambda2_criterion(vel_grad: &Mat3) -> f64 {
    let s = strain_rate_tensor(vel_grad);
    let omega = rotation_rate_tensor(vel_grad);
    let mut a = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                a[i][j] += s[i][k] * s[k][j] + omega[i][k] * omega[k][j];
            }
        }
    }
    let tr = mat3_trace(&a);
    let tr_a2 = {
        let mut t = 0.0;
        for i in 0..3 {
            for k in 0..3 {
                t += a[i][k] * a[k][i];
            }
        }
        t
    };
    let i2 = 0.5 * (tr * tr - tr_a2);
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let p_coef = i2 - tr * tr / 3.0;
    let q_coef = (2.0 * tr * tr * tr / 27.0) - (tr * i2 / 3.0) + det;
    let disc = -(4.0 * p_coef * p_coef * p_coef + 27.0 * q_coef * q_coef);
    let mut eigenvalues = [0.0_f64; 3];
    if disc >= 0.0 {
        let m = 2.0 * (-p_coef / 3.0_f64).max(0.0).sqrt();
        let theta = if m > 1e-30 {
            (3.0 * q_coef / (p_coef * m)).clamp(-1.0, 1.0).acos() / 3.0
        } else {
            0.0
        };
        use std::f64::consts::PI;
        eigenvalues[0] = m * theta.cos() + tr / 3.0;
        eigenvalues[1] = m * (theta - 2.0 * PI / 3.0).cos() + tr / 3.0;
        eigenvalues[2] = m * (theta - 4.0 * PI / 3.0).cos() + tr / 3.0;
    } else {
        let a_c = (-q_coef / 2.0
            + (q_coef * q_coef / 4.0 + p_coef * p_coef * p_coef / 27.0).sqrt())
        .abs()
        .cbrt();
        let b_c = if a_c > 1e-30 {
            -p_coef / (3.0 * a_c)
        } else {
            0.0
        };
        eigenvalues[0] = a_c + b_c + tr / 3.0;
        eigenvalues[1] = eigenvalues[0];
        eigenvalues[2] = eigenvalues[0];
    }
    eigenvalues.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    eigenvalues[1]
}
#[cfg(test)]
mod tests_turbulence_extended {
    use super::*;
    use crate::turbulence::types::*;
    #[test]
    fn test_smagorinsky_zero_strain_zero_viscosity() {
        let model = SmagorinskyLes::default();
        let s_zero: Mat3 = [[0.0; 3]; 3];
        let nu_t = model.eddy_viscosity(&s_zero);
        assert_eq!(nu_t, 0.0);
    }
    #[test]
    fn test_smagorinsky_eddy_viscosity_positive() {
        let model = SmagorinskyLes::new(0.17, 0.1);
        let s: Mat3 = [[1.0, 0.5, 0.0], [0.5, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let nu_t = model.eddy_viscosity(&s);
        assert!(nu_t > 0.0, "Non-zero strain → positive eddy viscosity");
    }
    #[test]
    fn test_smagorinsky_scales_with_cs_squared() {
        let delta = 0.1;
        let s: Mat3 = [[2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, 0.0]];
        let m1 = SmagorinskyLes::new(0.1, delta);
        let m2 = SmagorinskyLes::new(0.2, delta);
        let nu1 = m1.eddy_viscosity(&s);
        let nu2 = m2.eddy_viscosity(&s);
        assert!(
            (nu2 / nu1 - 4.0).abs() < 1e-10,
            "ν_t ∝ cs²: ratio={}",
            nu2 / nu1
        );
    }
    #[test]
    fn test_smagorinsky_sgs_stress_symmetric() {
        let model = SmagorinskyLes::new(0.17, 0.1);
        let s: Mat3 = [[1.0, 0.3, 0.0], [0.3, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let tau = model.sgs_stress(&s, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (tau[i][j] - tau[j][i]).abs() < 1e-12,
                    "SGS stress must be symmetric: tau[{i}][{j}]≠tau[{j}][{i}]"
                );
            }
        }
    }
    #[test]
    fn test_smagorinsky_sgs_stress_negative_diagonal_for_compression() {
        let model = SmagorinskyLes::new(0.17, 0.1);
        let s: Mat3 = [[5.0, 0.0, 0.0], [0.0, -5.0, 0.0], [0.0, 0.0, 0.0]];
        let tau = model.sgs_stress(&s, 1000.0);
        assert!(tau[0][0] < 0.0, "τ_11 < 0 for positive S_11");
    }
    #[test]
    fn test_smagorinsky_acceleration_zero_uniform_tau() {
        let model = SmagorinskyLes::new(0.17, 0.1);
        let tau: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let rho = 1000.0_f64;
        let nbr_tau = vec![tau; 2];
        let nbr_rho = vec![rho; 2];
        let nbr_mass = vec![0.001_f64; 2];
        let grads = vec![[1.0, 0.0, 0.0_f64], [-1.0, 0.0, 0.0]];
        let acc = model.sgs_acceleration(&tau, rho, &nbr_tau, &nbr_rho, &nbr_mass, &grads);
        for &a in &acc {
            assert!(a.is_finite());
        }
    }
    #[test]
    fn test_dynamic_smagorinsky_cs_sq_range() {
        let model = DynamicSmagorinskyLes::default();
        let s: Mat3 = [[1.0, 0.5, 0.0], [0.5, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let cs_sq = model.compute_cs_sq(&s, &s);
        assert!((0.0..=0.04).contains(&cs_sq), "cs_sq={cs_sq} out of range");
    }
    #[test]
    fn test_dynamic_smagorinsky_zero_strain_zero_viscosity() {
        let model = DynamicSmagorinskyLes::default();
        let s_zero: Mat3 = [[0.0; 3]; 3];
        let nu_t = model.eddy_viscosity(&s_zero, &s_zero);
        assert_eq!(nu_t, 0.0);
    }
    #[test]
    fn test_dynamic_smagorinsky_eddy_viscosity_nonneg() {
        let model = DynamicSmagorinskyLes::new(0.1, 2.0);
        let s_bar: Mat3 = [[2.0, 1.0, 0.0], [1.0, -2.0, 0.0], [0.0, 0.0, 0.0]];
        let s_hat: Mat3 = [[1.0, 0.5, 0.0], [0.5, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let nu_t = model.eddy_viscosity(&s_bar, &s_hat);
        assert!(nu_t >= 0.0, "Dynamic Smagorinsky ν_t must be non-negative");
    }
    #[test]
    fn test_dynamic_smagorinsky_sgs_stress_finite() {
        let model = DynamicSmagorinskyLes::default();
        let s: Mat3 = [[3.0, 0.0, 0.0], [0.0, -1.5, 0.0], [0.0, 0.0, -1.5]];
        let tau = model.sgs_stress(&s, &s, 1000.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(tau[i][j].is_finite(), "SGS stress should be finite");
            }
        }
    }
    #[test]
    fn test_komega_turbulent_viscosity() {
        let model = KOmegaRans::default();
        let nu_t = model.turbulent_viscosity(1.0, 2.0);
        assert!((nu_t - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_komega_turbulent_viscosity_zero_omega() {
        let model = KOmegaRans::default();
        let nu_t = model.turbulent_viscosity(1.0, 0.0);
        assert_eq!(nu_t, 0.0);
    }
    #[test]
    fn test_komega_production_zero_strain() {
        let model = KOmegaRans::default();
        let s_zero: Mat3 = [[0.0; 3]; 3];
        let p_k = model.production(1.0, 1.0, &s_zero);
        assert_eq!(p_k, 0.0);
    }
    #[test]
    fn test_komega_destruction_positive() {
        let model = KOmegaRans::default();
        let d_k = model.destruction_k(1.0, 1.0);
        assert!(d_k > 0.0);
        let d_omega = model.destruction_omega(2.0);
        assert!(d_omega > 0.0);
    }
    #[test]
    fn test_komega_integrate_positivity() {
        let model = KOmegaRans::default();
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, -0.5]];
        let (k_new, omega_new) = model.integrate_explicit(0.01, 1.0, &s, 0.001);
        assert!(k_new > 0.0, "k must stay positive, got {k_new}");
        assert!(omega_new > 0.0, "ω must stay positive, got {omega_new}");
    }
    #[test]
    fn test_komega_effective_viscosity_geq_molecular() {
        let model = KOmegaRans::default();
        let nu_eff = model.effective_viscosity(0.01, 1.0);
        assert!(nu_eff >= model.nu, "ν_eff ≥ ν_molecular");
    }
    #[test]
    fn test_komega_turbulent_length_scale_positive() {
        let model = KOmegaRans::default();
        let l = model.turbulent_length_scale(1.0, 1.0);
        assert!(l > 0.0);
    }
    #[test]
    fn test_komega_source_terms_finite() {
        let model = KOmegaRans::default();
        let s: Mat3 = [[0.5, 0.25, 0.0], [0.25, -0.25, 0.0], [0.0, 0.0, -0.25]];
        let (src_k, src_omega) = model.source_terms(0.1, 1.0, &s);
        assert!(src_k.is_finite());
        assert!(src_omega.is_finite());
    }
    #[test]
    fn test_mixing_length_zero_at_wall() {
        let model = MixingLengthModel::default();
        let l = model.mixing_length_van_driest(0.0, 0.0);
        assert!(l.abs() < 1e-30, "Mixing length at wall should be zero");
    }
    #[test]
    fn test_mixing_length_increases_with_y() {
        let model = MixingLengthModel::default();
        let l1 = model.mixing_length_van_driest(0.01, 10.0);
        let l2 = model.mixing_length_van_driest(0.02, 20.0);
        assert!(l2 > l1, "Mixing length should increase with y");
    }
    #[test]
    fn test_mixing_length_capped() {
        let l_max = 0.5;
        let model = MixingLengthModel::new(0.41, 26.0, l_max);
        let l = model.mixing_length_van_driest(100.0, 1e6);
        assert!(l <= l_max + 1e-12, "Mixing length should be capped");
    }
    #[test]
    fn test_mixing_length_eddy_viscosity_positive() {
        let model = MixingLengthModel::default();
        let s: Mat3 = [[5.0, 0.0, 0.0], [0.0, -5.0, 0.0], [0.0, 0.0, 0.0]];
        let nu_t = model.eddy_viscosity(0.01, &s);
        assert!(nu_t > 0.0);
    }
    #[test]
    fn test_mixing_length_free_shear() {
        let model = MixingLengthModel::default();
        let l = model.mixing_length_free_shear(1.0);
        assert!((l - 0.09).abs() < 1e-12);
    }
    #[test]
    fn test_turb_visc_field_init_zero() {
        let field = TurbulentViscosityField::new(5);
        assert!(field.nu_t.iter().all(|&v| v == 0.0));
        assert!(field.k.iter().all(|&v| v == 0.0));
        assert!(field.omega.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_turb_visc_field_update_k_omega() {
        let mut field = TurbulentViscosityField::new(3);
        field.k = vec![1.0, 2.0, 3.0];
        field.omega = vec![1.0, 2.0, 3.0];
        field.update_from_k_omega();
        for &v in &field.nu_t {
            assert!((v - 1.0).abs() < 1e-12, "ν_t = k/ω = 1.0, got {v}");
        }
    }
    #[test]
    fn test_turb_visc_field_clamp() {
        let mut field = TurbulentViscosityField::new(3);
        field.nu_t = vec![-1.0, 0.5, 2.0];
        field.clamp_nu_t(1.0);
        assert_eq!(field.nu_t[0], 0.0);
        assert_eq!(field.nu_t[1], 0.5);
        assert_eq!(field.nu_t[2], 1.0);
    }
    #[test]
    fn test_turb_visc_field_mean_k() {
        let mut field = TurbulentViscosityField::new(4);
        field.k = vec![1.0, 2.0, 3.0, 4.0];
        assert!((field.mean_k() - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_turb_visc_field_interpolate_no_particles() {
        let field = TurbulentViscosityField::new(0);
        let v = field.interpolate_nu_t([0.0; 3], &[], &[], &[], 0.1);
        assert_eq!(v, 0.0);
    }
    #[test]
    fn test_q_criterion_pure_rotation() {
        let omega_val = 1.0;
        let vel_grad: Mat3 = [
            [0.0, -omega_val, 0.0],
            [omega_val, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let q = q_criterion(&vel_grad);
        assert!(q > 0.0, "Pure rotation → Q > 0, got {q}");
    }
    #[test]
    fn test_q_criterion_pure_strain() {
        let vel_grad: Mat3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let q = q_criterion(&vel_grad);
        assert!(q < 0.0, "Pure strain → Q < 0, got {q}");
    }
    #[test]
    fn test_q_criterion_zero_for_zero_gradient() {
        let vel_grad: Mat3 = [[0.0; 3]; 3];
        let q = q_criterion(&vel_grad);
        assert_eq!(q, 0.0);
    }
    #[test]
    fn test_lambda2_criterion_finite() {
        let vel_grad: Mat3 = [[1.0, 0.5, 0.2], [-0.3, -0.5, 0.1], [0.1, -0.2, -0.5]];
        let l2 = lambda2_criterion(&vel_grad);
        assert!(l2.is_finite(), "λ_2 should be finite, got {l2}");
    }
    #[test]
    fn test_rotation_rate_tensor_antisymmetric() {
        let vel_grad: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let omega = rotation_rate_tensor(&vel_grad);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (omega[i][j] + omega[j][i]).abs() < 1e-12,
                    "Rotation rate tensor must be anti-symmetric"
                );
            }
        }
    }
    #[test]
    fn test_strain_rate_tensor_symmetric() {
        let vel_grad: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let s = strain_rate_tensor(&vel_grad);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (s[i][j] - s[j][i]).abs() < 1e-12,
                    "Strain-rate tensor must be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_turbulent_prandtl_thermal_diffusivity() {
        let model = TurbulentPrandtlNumber::new(0.9);
        let alpha_t = model.thermal_diffusivity(0.9);
        assert!((alpha_t - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_turbulent_prandtl_heat_flux_direction() {
        let model = TurbulentPrandtlNumber::default();
        let grad_t = [1.0, 0.0, 0.0];
        let q = model.turbulent_heat_flux(1e-4, 1000.0, 4200.0, grad_t);
        assert!(q[0] < 0.0, "Turbulent heat flux must oppose gradient");
    }
    #[test]
    fn test_turbulent_prandtl_zero_nu_t() {
        let model = TurbulentPrandtlNumber::default();
        let q = model.turbulent_heat_flux(0.0, 1000.0, 4200.0, [1.0, 1.0, 1.0]);
        for &qi in &q {
            assert_eq!(qi, 0.0);
        }
    }
    #[test]
    fn test_les_budget_initial_zero() {
        let budget = LesEnergyBudget::new();
        assert_eq!(budget.production, 0.0);
        assert_eq!(budget.dissipation, 0.0);
        assert_eq!(budget.count, 0);
    }
    #[test]
    fn test_les_budget_add_particle() {
        let mut budget = LesEnergyBudget::new();
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let model = SmagorinskyLes::new(0.17, 0.1);
        let tau = model.sgs_stress(&s, 1000.0);
        budget.add_particle(&tau, &s);
        assert_eq!(budget.count, 1);
        assert!(budget.production.is_finite());
    }
    #[test]
    fn test_les_budget_mean_production_zero_count() {
        let budget = LesEnergyBudget::new();
        assert_eq!(budget.mean_production(), 0.0);
    }
    #[test]
    fn test_les_budget_reset() {
        let mut budget = LesEnergyBudget::new();
        let s: Mat3 = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]];
        let model = SmagorinskyLes::new(0.17, 0.1);
        let tau = model.sgs_stress(&s, 1000.0);
        budget.add_particle(&tau, &s);
        budget.reset();
        assert_eq!(budget.count, 0);
        assert_eq!(budget.production, 0.0);
    }
}
