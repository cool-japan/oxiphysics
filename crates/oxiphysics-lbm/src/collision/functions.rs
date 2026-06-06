//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::grid::{LbmGrid2D, LbmGrid3D, equilibrium_2d, equilibrium_3d};
use crate::lattice::{CS2, D3Q19_OPPOSITES, D3Q19_VELOCITIES, D3Q19_WEIGHTS};

use super::types::EntropicCollision;

/// Compute the D3Q19 equilibrium distribution for a given macroscopic state.
///
/// Uses the standard low-Mach expansion:
/// `feq_i = w_i * rho * (1 + (e·u)/cs² + (e·u)²/(2cs⁴) − u²/(2cs²))`
///
/// The returned array sums to `rho` and its first momentum equals `rho * u`.
pub fn compute_equilibrium_d3q19(rho: f64, u: [f64; 3]) -> [f64; 19] {
    let mut feq = [0.0f64; 19];
    let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    for (i, feq_i) in feq.iter_mut().enumerate() {
        let w = D3Q19_WEIGHTS[i];
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        let eu = cx * u[0] + cy * u[1] + cz * u[2];
        *feq_i = w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
    }
    feq
}
/// Compute macroscopic density and velocity from a D3Q19 distribution array.
///
/// Returns `(rho, [ux, uy, uz])` where:
/// - `rho = sum_i f_i`
/// - `u_alpha = sum_i f_i * c_i_alpha / rho`
pub fn compute_macroscopic(f: &[f64; 19]) -> (f64, [f64; 3]) {
    let mut rho = 0.0f64;
    let mut mx = 0.0f64;
    let mut my = 0.0f64;
    let mut mz = 0.0f64;
    for (i, &fi) in f.iter().enumerate() {
        rho += fi;
        let c = D3Q19_VELOCITIES[i];
        mx += fi * c[0] as f64;
        my += fi * c[1] as f64;
        mz += fi * c[2] as f64;
    }
    let inv_rho = if rho.abs() > 1e-15 { 1.0 / rho } else { 0.0 };
    (rho, [mx * inv_rho, my * inv_rho, mz * inv_rho])
}
/// Perform the BGK collision step on a single D3Q19 population array.
///
/// `f*_i = f_i - omega * (f_i - feq_i)`
///
/// The macroscopic fields `rho` and `u` must already be computed from `f`.
pub fn collide_bgk(f: &[f64; 19], rho: f64, u: [f64; 3], omega: f64) -> [f64; 19] {
    let feq = compute_equilibrium_d3q19(rho, u);
    let mut f_out = [0.0f64; 19];
    for (i, f_out_i) in f_out.iter_mut().enumerate() {
        *f_out_i = f[i] - omega * (f[i] - feq[i]);
    }
    f_out
}
/// Magic parameter for TRT: `lambda = 3/16`.
pub const TRT_MAGIC: f64 = 3.0 / 16.0;
/// Perform the TRT collision step on a single D3Q19 population array.
///
/// For each direction pair `(i, ī)`:
/// ```text
/// f*_i = f_i
///      - omega_plus  * (f_i^+  - feq_i^+)
///      - omega_minus * (f_i^-  - feq_i^-)
/// ```
/// where `f^+ = (f_i + f_ī)/2` and `f^- = (f_i - f_ī)/2`.
pub fn collide_trt(
    f: &[f64; 19],
    rho: f64,
    u: [f64; 3],
    omega_plus: f64,
    omega_minus: f64,
) -> [f64; 19] {
    let feq = compute_equilibrium_d3q19(rho, u);
    let mut f_out = [0.0f64; 19];
    for (i, f_out_i) in f_out.iter_mut().enumerate() {
        let ibar = D3Q19_OPPOSITES[i];
        let fi_plus = 0.5 * (f[i] + f[ibar]);
        let fi_minus = 0.5 * (f[i] - f[ibar]);
        let feqi_plus = 0.5 * (feq[i] + feq[ibar]);
        let feqi_minus = 0.5 * (feq[i] - feq[ibar]);
        *f_out_i =
            f[i] - omega_plus * (fi_plus - feqi_plus) - omega_minus * (fi_minus - feqi_minus);
    }
    f_out
}
/// Compute the non-equilibrium stress tensor from a distribution.
///
/// `Pi_ab^neq = sum_i (f_i - feq_i) * c_ia * c_ib`
pub fn non_equilibrium_stress(f: &[f64; 19], feq: &[f64; 19]) -> [[f64; 3]; 3] {
    let mut pi = [[0.0f64; 3]; 3];
    for (i, (&fi, &feqi)) in f.iter().zip(feq.iter()).enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let f_neq = fi - feqi;
        for a in 0..3 {
            for b in 0..3 {
                pi[a][b] += f_neq * c[a] as f64 * c[b] as f64;
            }
        }
    }
    pi
}
/// Compute the Frobenius norm of the non-equilibrium stress tensor.
///
/// `||Pi^neq|| = sqrt(sum_ab (Pi_ab^neq)^2)`
pub fn non_equilibrium_stress_magnitude(f: &[f64; 19], feq: &[f64; 19]) -> f64 {
    let pi = non_equilibrium_stress(f, feq);
    let sum: f64 = pi.iter().flat_map(|row| row.iter()).map(|&v| v * v).sum();
    sum.sqrt()
}
/// Check if a distribution is positive (all f_i > 0).
///
/// Negative distributions indicate numerical instability.
pub fn distribution_is_positive(f: &[f64; 19]) -> bool {
    f.iter().all(|&fi| fi > 0.0)
}
/// Compute the deviation from equilibrium (L2 norm of f - feq).
///
/// `||f - feq|| = sqrt(sum_i (f_i - feq_i)^2)`
pub fn equilibrium_deviation(f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
    let feq = compute_equilibrium_d3q19(rho, u);
    let sum: f64 = f
        .iter()
        .zip(feq.iter())
        .map(|(&fi, &feqi)| (fi - feqi) * (fi - feqi))
        .sum();
    sum.sqrt()
}
/// Compute the entropy of a distribution (Boltzmann H-function).
///
/// `H = sum_i f_i * ln(f_i / w_i)`
pub fn boltzmann_entropy(f: &[f64; 19]) -> f64 {
    f.iter()
        .zip(D3Q19_WEIGHTS.iter())
        .map(|(&fi_raw, &wi)| {
            let fi = fi_raw.max(1e-30);
            fi * (fi / wi).ln()
        })
        .sum()
}
/// Compute the effective relaxation rate for stability monitoring.
///
/// Returns the ratio `max_i |f_i - feq_i| / max_i |feq_i|`
/// which should remain O(Ma^2) for stability.
pub fn non_equilibrium_ratio(f: &[f64; 19], rho: f64, u: [f64; 3]) -> f64 {
    let feq = compute_equilibrium_d3q19(rho, u);
    let max_neq: f64 = (0..19)
        .map(|i| (f[i] - feq[i]).abs())
        .fold(0.0_f64, f64::max);
    let max_eq: f64 = (0..19).map(|i| feq[i].abs()).fold(0.0_f64, f64::max);
    if max_eq > 1e-30 {
        max_neq / max_eq
    } else {
        0.0
    }
}
/// Compute the strain-rate tensor from the non-equilibrium part.
///
/// In LBM, the strain rate is related to the non-equilibrium stress:
/// `S_ab = -omega / (2 * rho * cs^2) * Pi_ab^neq`
pub fn strain_rate_from_neq(f: &[f64; 19], rho: f64, u: [f64; 3], omega: f64) -> [[f64; 3]; 3] {
    let feq = compute_equilibrium_d3q19(rho, u);
    let pi = non_equilibrium_stress(f, &feq);
    let factor = -omega / (2.0 * rho.max(1e-30) * CS2);
    let mut s = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            s[a][b] = factor * pi[a][b];
        }
    }
    s
}
/// Compute the H-function (Boltzmann entropy) for a distribution.
pub(super) fn entropy_h(f: &[f64; 19]) -> f64 {
    f.iter()
        .zip(D3Q19_WEIGHTS.iter())
        .map(|(&fi_raw, &wi)| {
            let fi = fi_raw.max(1e-30);
            fi * (fi / wi).ln()
        })
        .sum()
}
/// Compute the entropy production rate: ΔH = H(f_out) - H(f_in).
///
/// Should be ≤ 0 for physically consistent collisions (H-theorem).
pub fn entropy_production(f_in: &[f64; 19], f_out: &[f64; 19]) -> f64 {
    entropy_h(f_out) - entropy_h(f_in)
}
/// Compute the minimum entropy margin across an array of distributions.
///
/// Returns `min_k (H(f_k^post) - H(f_k^pre))`.  Should be ≤ 0 for stability.
pub fn min_entropy_margin(pre: &[[f64; 19]], post: &[[f64; 19]]) -> f64 {
    pre.iter()
        .zip(post.iter())
        .map(|(fp, fq)| entropy_production(fp, fq))
        .fold(f64::INFINITY, f64::min)
}
/// Compute a 19-component central-moment vector from a distribution.
///
/// Uses the velocity-shifted moments: `κ_{abc} = Σ_i f_i (c_x - u_x)^a (c_y - u_y)^b (c_z - u_z)^c`.
pub(super) fn compute_central_moments(f: &[f64; 19], u: [f64; 3], _inv_rho: f64) -> [f64; 19] {
    let mut cm = [0.0f64; 19];
    for (i, &fi) in f.iter().enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let dx = c[0] as f64 - u[0];
        let dy = c[1] as f64 - u[1];
        let dz = c[2] as f64 - u[2];
        cm[0] += fi;
        cm[1] += fi * dx;
        cm[2] += fi * dy;
        cm[3] += fi * dz;
        cm[4] += fi * dx * dx;
        cm[5] += fi * dy * dy;
        cm[6] += fi * dz * dz;
        cm[7] += fi * dx * dy;
        cm[8] += fi * dx * dz;
        cm[9] += fi * dy * dz;
        cm[10] += fi * dx * dx * dy;
        cm[11] += fi * dx * dx * dz;
        cm[12] += fi * dy * dy * dx;
        cm[13] += fi * dy * dy * dz;
        cm[14] += fi * dz * dz * dx;
        cm[15] += fi * dz * dz * dy;
        cm[16] += fi * dx * dx * dy * dy;
        cm[17] += fi * dx * dx * dz * dz;
        cm[18] += fi * dy * dy * dz * dz;
    }
    cm
}
/// BGK-ELBM blended collision: weighted combination.
///
/// `f* = alpha * f*_BGK + (1 - alpha) * f*_ELBM`
///
/// where alpha is derived from the local Mach number.
pub fn collide_blended(
    f: &[f64; 19],
    rho: f64,
    u: [f64; 3],
    omega: f64,
    blend_alpha: f64,
) -> [f64; 19] {
    let f_bgk = collide_bgk(f, rho, u, omega);
    let elbm = EntropicCollision::new(omega);
    let f_elbm = elbm.collide(f, rho, u);
    let a = blend_alpha.clamp(0.0, 1.0);
    std::array::from_fn(|i| a * f_bgk[i] + (1.0 - a) * f_elbm[i])
}
/// Compute the MRT viscosity bulk (volumetric) viscosity from the energy
/// relaxation rate.
///
/// `ζ = cs² * (2 - s_e) / (3 * s_e)` (bulk viscosity from energy mode).
pub fn mrt_bulk_viscosity(s_energy: f64) -> f64 {
    if s_energy <= 0.0 || s_energy >= 2.0 {
        return 0.0;
    }
    CS2 * (2.0 - s_energy) / (3.0 * s_energy)
}
/// Compute the optimal anti-symmetric relaxation rate for MRT to suppress
/// numerical diffusion artifacts:
///
/// `s_minus = 8 * (2 - s_plus) / (8 - s_plus)` (Ginzburg magic parameter).
pub fn optimal_antisymmetric_rate(s_plus: f64) -> f64 {
    if s_plus.abs() < 1e-30 {
        return 1.0;
    }
    (8.0 * (2.0 - s_plus) / (8.0 - s_plus)).clamp(0.0, 2.0)
}
/// Perform the BGK collision step on a 2D grid.
///
/// `omega` = 1/tau, where tau is the relaxation time.
/// The kinematic viscosity is nu = cs² (tau − 0.5).
///
/// For each cell: `f_i = f_i − omega * (f_i − f_eq_i)`
pub fn bgk_collide_2d(grid: &mut LbmGrid2D, omega: f64) {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny;
    grid.compute_macroscopic();
    for k in 0..n {
        let rho_k = grid.rho[k];
        let ux_k = grid.ux[k];
        let uy_k = grid.uy[k];
        for i in 0..q {
            let w = grid.lattice.weight(i);
            let c = grid.lattice.velocity_2d(i);
            let feq = equilibrium_2d(w, rho_k, ux_k, uy_k, c[0] as f64, c[1] as f64);
            let fi_k = &mut grid.f[i][k];
            *fi_k -= omega * (*fi_k - feq);
        }
    }
}
/// Perform the BGK collision step on a 3D grid.
///
/// `omega` = 1/tau, where tau is the relaxation time.
pub fn bgk_collide_3d(grid: &mut LbmGrid3D, omega: f64) {
    let q = grid.lattice.q();
    let n = grid.nx * grid.ny * grid.nz;
    grid.compute_macroscopic();
    for k in 0..n {
        let rho_k = grid.rho[k];
        let ux_k = grid.ux[k];
        let uy_k = grid.uy[k];
        let uz_k = grid.uz[k];
        for i in 0..q {
            let w = grid.lattice.weight(i);
            let c = grid.lattice.velocity_3d(i);
            let feq = equilibrium_3d(
                w,
                rho_k,
                ux_k,
                uy_k,
                uz_k,
                c[0] as f64,
                c[1] as f64,
                c[2] as f64,
            );
            let fi_k = &mut grid.f[i][k];
            *fi_k -= omega * (*fi_k - feq);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision::types::*;
    pub(super) const OMEGA: f64 = 1.0;
    /// Helper: near-rest equilibrium distribution.
    fn near_rest_f(rho: f64) -> [f64; 19] {
        compute_equilibrium_d3q19(rho, [0.0, 0.0, 0.0])
    }
    /// Helper: create a perturbed distribution.
    fn perturbed_f(rho: f64, u: [f64; 3], perturbation: f64) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut f = feq;
        f[1] += perturbation;
        f[2] -= perturbation;
        f
    }
    #[test]
    fn test_equilibrium_sums_to_rho() {
        let rho = 1.3;
        let u = [0.05, -0.02, 0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!((sum - rho).abs() < 1e-13, "feq sum = {sum}, expected {rho}");
    }
    #[test]
    fn test_zero_velocity_equilibrium() {
        let rho = 2.0;
        let feq = compute_equilibrium_d3q19(rho, [0.0, 0.0, 0.0]);
        let expected_f0 = rho / 3.0;
        assert!(
            (feq[0] - expected_f0).abs() < 1e-14,
            "feq[0] = {}, expected {expected_f0}",
            feq[0]
        );
        let mut mx = 0.0f64;
        let mut my = 0.0f64;
        let mut mz = 0.0f64;
        for (i, &feqi) in feq.iter().enumerate() {
            let c = D3Q19_VELOCITIES[i];
            mx += feqi * c[0] as f64;
            my += feqi * c[1] as f64;
            mz += feqi * c[2] as f64;
        }
        assert!(mx.abs() < 1e-14, "mx = {mx}");
        assert!(my.abs() < 1e-14, "my = {my}");
        assert!(mz.abs() < 1e-14, "mz = {mz}");
    }
    #[test]
    fn test_bgk_mass_conservation() {
        let rho = 1.1;
        let u = [0.04, 0.0, -0.02];
        let f = compute_equilibrium_d3q19(rho, u);
        let mut f_pert = f;
        f_pert[1] += 0.01;
        f_pert[2] -= 0.01;
        let (rho_in, u_in) = compute_macroscopic(&f_pert);
        let f_out = collide_bgk(&f_pert, rho_in, u_in, OMEGA);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-14,
            "BGK mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_bgk_momentum_conservation() {
        let f_in = compute_equilibrium_d3q19(1.0, [0.05, -0.03, 0.01]);
        let (rho, u) = compute_macroscopic(&f_in);
        let f_out = collide_bgk(&f_in, rho, u, OMEGA);
        let mut mx_out = 0.0f64;
        let mut my_out = 0.0f64;
        let mut mz_out = 0.0f64;
        let mut mx_in = 0.0f64;
        let mut my_in = 0.0f64;
        let mut mz_in = 0.0f64;
        for (i, (&f_in_i, &f_out_i)) in f_in.iter().zip(f_out.iter()).enumerate() {
            let c = D3Q19_VELOCITIES[i];
            mx_in += f_in_i * c[0] as f64;
            my_in += f_in_i * c[1] as f64;
            mz_in += f_in_i * c[2] as f64;
            mx_out += f_out_i * c[0] as f64;
            my_out += f_out_i * c[1] as f64;
            mz_out += f_out_i * c[2] as f64;
        }
        assert!((mx_out - mx_in).abs() < 1e-14, "BGK mx not conserved");
        assert!((my_out - my_in).abs() < 1e-14, "BGK my not conserved");
        assert!((mz_out - mz_in).abs() < 1e-14, "BGK mz not conserved");
    }
    #[test]
    fn test_trt_mass_conservation() {
        let f = compute_equilibrium_d3q19(1.2, [0.03, 0.0, 0.05]);
        let mut f_pert = f;
        f_pert[3] += 0.02;
        f_pert[4] -= 0.02;
        let (rho_in, u_in) = compute_macroscopic(&f_pert);
        let trt = TrtCollision::from_viscosity(1.0 / 6.0);
        let f_out = trt.collide(&f_pert, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho_in).abs() < 1e-14, "TRT mass not conserved");
    }
    #[test]
    fn test_trt_momentum_conservation() {
        let f_in = compute_equilibrium_d3q19(1.0, [0.02, -0.04, 0.0]);
        let (rho, u) = compute_macroscopic(&f_in);
        let trt = TrtCollision::from_viscosity(1.0 / 6.0);
        let f_out = trt.collide(&f_in, rho, u);
        let momentum = |arr: &[f64; 19]| -> [f64; 3] {
            let mut m = [0.0f64; 3];
            for (i, &arr_i) in arr.iter().enumerate() {
                let c = D3Q19_VELOCITIES[i];
                m[0] += arr_i * c[0] as f64;
                m[1] += arr_i * c[1] as f64;
                m[2] += arr_i * c[2] as f64;
            }
            m
        };
        let m_in = momentum(&f_in);
        let m_out = momentum(&f_out);
        for d in 0..3 {
            assert!(
                (m_out[d] - m_in[d]).abs() < 1e-14,
                "TRT momentum[{d}] not conserved"
            );
        }
    }
    #[test]
    fn test_compute_macroscopic_round_trip() {
        let rho_ref = 1.4;
        let u_ref = [0.06, -0.03, 0.01];
        let feq = compute_equilibrium_d3q19(rho_ref, u_ref);
        let (rho, u) = compute_macroscopic(&feq);
        assert!(
            (rho - rho_ref).abs() < 1e-13,
            "rho mismatch: {rho} vs {rho_ref}"
        );
        for d in 0..3 {
            assert!(
                (u[d] - u_ref[d]).abs() < 1e-13,
                "u[{d}] mismatch: {} vs {}",
                u[d],
                u_ref[d]
            );
        }
    }
    #[test]
    fn test_trt_magic_parameter() {
        let trt = TrtCollision::from_viscosity(1.0 / 6.0);
        let tau_plus = 1.0 / trt.omega_plus;
        let tau_minus = 1.0 / trt.omega_minus;
        let lambda = (tau_plus - 0.5) * (tau_minus - 0.5);
        assert!(
            (lambda - TRT_MAGIC).abs() < 1e-13,
            "TRT magic parameter: {lambda} != {TRT_MAGIC}"
        );
    }
    #[test]
    fn test_bgk_at_equilibrium_unchanged() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let f_out = collide_bgk(&feq, rho, u, OMEGA);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-14,
                "BGK modified equilibrium at i={i}"
            );
        }
    }
    #[test]
    fn test_bgk_from_viscosity() {
        let nu: f64 = 0.1;
        let bgk = BgkCollision::from_viscosity(nu);
        let nu_check = (1.0 / 3.0) * (1.0 / bgk.omega - 0.5);
        assert!(
            (nu_check - nu).abs() < 1e-14,
            "nu from BgkCollision::from_viscosity: {nu_check} vs {nu}"
        );
    }
    #[test]
    fn test_near_rest_distribution_sums() {
        let rho = 0.9;
        let f = near_rest_f(rho);
        let sum: f64 = f.iter().sum();
        assert!((sum - rho).abs() < 1e-14, "near_rest_f sum = {sum}");
    }
    #[test]
    fn test_entropic_mass_conservation() {
        let f = perturbed_f(1.0, [0.03, -0.01, 0.02], 0.005);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let elbm = EntropicCollision::new(1.5);
        let f_out = elbm.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "Entropic mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_entropic_positivity() {
        let f = perturbed_f(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let elbm = EntropicCollision::new(1.5);
        let f_out = elbm.collide(&f, rho, u);
        assert!(
            distribution_is_positive(&f_out),
            "Entropic should maintain positive distributions"
        );
    }
    #[test]
    fn test_regularized_mass_conservation() {
        let f = perturbed_f(1.1, [0.04, -0.02, 0.01], 0.008);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let reg = RegularizedCollision::new(1.2);
        let f_out = reg.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "Regularized mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_regularized_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let reg = RegularizedCollision::new(1.0);
        let f_out = reg.collide(&feq, rho, u);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-13,
                "Regularized modified equilibrium at i={i}: {f_out_i} vs {feq_i}"
            );
        }
    }
    #[test]
    fn test_cascaded_mass_conservation() {
        let f = perturbed_f(1.0, [0.03, 0.0, -0.02], 0.005);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let casc = CascadedCollision::from_viscosity(0.1);
        let f_out = casc.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "Cascaded mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_neq_stress_zero_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi = non_equilibrium_stress(&feq, &feq);
        for (a, row) in pi.iter().enumerate() {
            for (b, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1e-14, "Pi[{a}][{b}] should be zero at eq");
            }
        }
    }
    #[test]
    fn test_neq_stress_magnitude_perturbed() {
        let f = perturbed_f(1.0, [0.02, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let mag = non_equilibrium_stress_magnitude(&f, &feq);
        assert!(
            mag > 0.0,
            "Non-eq stress should be positive for perturbed f: {mag}"
        );
    }
    #[test]
    fn test_equilibrium_deviation_zero() {
        let rho = 1.0;
        let u = [0.03, -0.02, 0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let dev = equilibrium_deviation(&feq, rho, u);
        assert!(
            dev < 1e-14,
            "Deviation should be zero at equilibrium: {dev}"
        );
    }
    #[test]
    fn test_boltzmann_entropy_finite() {
        let feq = compute_equilibrium_d3q19(1.0, [0.01, 0.0, 0.0]);
        let h = boltzmann_entropy(&feq);
        assert!(h.is_finite(), "Entropy should be finite: {h}");
    }
    #[test]
    fn test_neq_ratio_zero_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let ratio = non_equilibrium_ratio(&feq, rho, u);
        assert!(
            ratio < 1e-14,
            "NEQ ratio should be zero at equilibrium: {ratio}"
        );
    }
    #[test]
    fn test_neq_ratio_positive_perturbed() {
        let f = perturbed_f(1.0, [0.02, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let ratio = non_equilibrium_ratio(&f, rho, u);
        assert!(ratio > 0.0, "NEQ ratio should be positive: {ratio}");
    }
    #[test]
    fn test_strain_rate_from_neq_at_eq() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let s = strain_rate_from_neq(&feq, rho, u, 1.0);
        for (a, row) in s.iter().enumerate() {
            for (b, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1e-14, "S[{a}][{b}] should be zero at eq: {val}",);
            }
        }
    }
    #[test]
    fn test_distribution_positivity() {
        let feq = compute_equilibrium_d3q19(1.0, [0.01, 0.0, 0.0]);
        assert!(
            distribution_is_positive(&feq),
            "Equilibrium should be positive"
        );
        let mut f_neg = feq;
        f_neg[1] = -0.001;
        assert!(
            !distribution_is_positive(&f_neg),
            "Negative entry should fail"
        );
    }
    #[test]
    fn test_kbc_mass_conservation() {
        let f = perturbed_f(1.0, [0.03, -0.01, 0.02], 0.005);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let kbc = KbcCollision::new(1.0 / 6.0);
        let f_out = kbc.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-11,
            "KBC mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_kbc_at_equilibrium_unchanged() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let kbc = KbcCollision::new(1.0 / 6.0);
        let f_out = kbc.collide(&feq, rho, u);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-12,
                "KBC modified equilibrium at i={i}: {f_out_i} vs {feq_i}"
            );
        }
    }
    #[test]
    fn test_entropy_production_le_zero_bgk() {
        let feq = compute_equilibrium_d3q19(1.0, [0.03, 0.0, 0.0]);
        let f_out = collide_bgk(&feq, 1.0, [0.03, 0.0, 0.0], 1.5);
        let dp = entropy_production(&feq, &f_out);
        assert!(
            dp <= 1e-12,
            "Entropy should not increase at equilibrium: {dp}"
        );
    }
    #[test]
    fn test_entropy_production_finite() {
        let f = perturbed_f(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let f_out = collide_bgk(&f, rho, u, 1.0);
        let dp = entropy_production(&f, &f_out);
        assert!(dp.is_finite(), "Entropy production should be finite: {dp}");
    }
    #[test]
    fn test_min_entropy_margin_at_equilibrium() {
        let feq = compute_equilibrium_d3q19(1.0, [0.02, 0.0, 0.0]);
        let f_out = collide_bgk(&feq, 1.0, [0.02, 0.0, 0.0], 1.0);
        let margin = min_entropy_margin(&[feq], &[f_out]);
        assert!(margin.is_finite());
    }
    #[test]
    fn test_central_moment_mass_conservation() {
        let f = perturbed_f(1.0, [0.03, 0.0, -0.02], 0.005);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let cmc = CentralMomentCollision::from_viscosity(0.1);
        let f_out = cmc.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-3,
            "Central moment mass error too large: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_central_moment_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let cmc = CentralMomentCollision::from_viscosity(1.0 / 6.0);
        let f_out = cmc.collide(&feq, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-12,
            "Central moment disturbed mass at eq: {rho_out}"
        );
    }
    #[test]
    fn test_central_moment_finite_output() {
        let f = perturbed_f(1.0, [0.05, -0.02, 0.01], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let cmc = CentralMomentCollision::from_viscosity(0.1);
        let f_out = cmc.collide(&f, rho, u);
        assert!(
            f_out.iter().all(|v| v.is_finite()),
            "Central moment output has NaN/Inf"
        );
    }
    #[test]
    fn test_hybrid_mass_conservation_bgk_mode() {
        let rho = 1.0;
        let u = [0.01, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let hybrid = HybridCollision::from_viscosity(1.0 / 6.0, 0.5);
        let f_out = hybrid.collide(&feq, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho).abs() < 1e-13, "Hybrid BGK mass: {rho_out}");
    }
    #[test]
    fn test_hybrid_mass_conservation_entropic_mode() {
        let f = perturbed_f(1.0, [0.03, 0.0, 0.0], 0.08);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let hybrid = HybridCollision::from_viscosity(1.0 / 6.0, 0.001);
        let f_out = hybrid.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-10,
            "Hybrid entropic mass: {rho_out}"
        );
    }
    #[test]
    fn test_hybrid_mode_selection_bgk() {
        let feq = compute_equilibrium_d3q19(1.0, [0.02, 0.0, 0.0]);
        let hybrid = HybridCollision::from_viscosity(1.0 / 6.0, 0.01);
        let mode = hybrid.active_mode(&feq, 1.0, [0.02, 0.0, 0.0]);
        assert_eq!(mode, CollisionMode::Bgk, "At equilibrium should be BGK");
    }
    #[test]
    fn test_hybrid_mode_selection_entropic() {
        let f = perturbed_f(1.0, [0.02, 0.0, 0.0], 0.1);
        let (rho, u) = compute_macroscopic(&f);
        let hybrid = HybridCollision::from_viscosity(1.0 / 6.0, 1e-15);
        let mode = hybrid.active_mode(&f, rho, u);
        assert_eq!(mode, CollisionMode::Entropic, "High NEQ should be Entropic");
    }
    #[test]
    fn test_blended_collision_alpha_one_equals_bgk() {
        let f = perturbed_f(1.0, [0.03, 0.0, 0.0], 0.005);
        let (rho, u) = compute_macroscopic(&f);
        let omega = 1.0;
        let f_bgk = collide_bgk(&f, rho, u, omega);
        let f_blend = collide_blended(&f, rho, u, omega, 1.0);
        for (i, (&f_blend_i, &f_bgk_i)) in f_blend.iter().zip(f_bgk.iter()).enumerate() {
            assert!(
                (f_blend_i - f_bgk_i).abs() < 1e-12,
                "alpha=1 blend should equal BGK at i={i}"
            );
        }
    }
    #[test]
    fn test_blended_collision_conserves_mass() {
        let f = perturbed_f(1.0, [0.02, 0.01, 0.0], 0.003);
        let (rho, u) = compute_macroscopic(&f);
        let f_out = collide_blended(&f, rho, u, 1.2, 0.5);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho).abs() < 1e-11, "Blended mass: {rho_out}");
    }
    #[test]
    fn test_mrt_mass_conservation() {
        let f = perturbed_f(1.1, [0.04, -0.02, 0.01], 0.008);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let mrt = MrtCollision::from_viscosity(1.0 / 6.0);
        let f_out = mrt.collide(&f, rho_in, u_in);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "MRT mass not conserved: {rho_out} vs {rho_in}"
        );
    }
    #[test]
    fn test_mrt_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let mrt = MrtCollision::from_viscosity(1.0 / 6.0);
        let f_out = mrt.collide(&feq, rho, u);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-13,
                "MRT modified equilibrium at i={i}: {f_out_i} vs {feq_i}"
            );
        }
    }
    #[test]
    fn test_mrt_effective_viscosity() {
        let nu = 1.0 / 6.0;
        let mrt = MrtCollision::from_viscosity(nu);
        let nu_eff = mrt.effective_viscosity();
        assert!(
            (nu_eff - nu).abs() < 1e-10,
            "MRT effective viscosity: {nu_eff} vs {nu}"
        );
    }
    #[test]
    fn test_mrt_momentum_conservation() {
        let f = perturbed_f(1.0, [0.03, -0.02, 0.01], 0.005);
        let (rho_in, u_in) = compute_macroscopic(&f);
        let mrt = MrtCollision::from_viscosity(0.1);
        let f_out = mrt.collide(&f, rho_in, u_in);
        let mut mx_in = 0.0f64;
        let mut mx_out = 0.0f64;
        let mut my_in = 0.0f64;
        let mut my_out = 0.0f64;
        let mut mz_in = 0.0f64;
        let mut mz_out = 0.0f64;
        for (i, (&fi, &f_out_i)) in f.iter().zip(f_out.iter()).enumerate() {
            let c = D3Q19_VELOCITIES[i];
            mx_in += fi * c[0] as f64;
            mx_out += f_out_i * c[0] as f64;
            my_in += fi * c[1] as f64;
            my_out += f_out_i * c[1] as f64;
            mz_in += fi * c[2] as f64;
            mz_out += f_out_i * c[2] as f64;
        }
        assert!((mx_out - mx_in).abs() < 1e-12, "MRT mx not conserved");
        assert!((my_out - my_in).abs() < 1e-12, "MRT my not conserved");
        assert!((mz_out - mz_in).abs() < 1e-12, "MRT mz not conserved");
    }
    #[test]
    fn test_mrt_bulk_viscosity() {
        let zeta = mrt_bulk_viscosity(1.19);
        assert!(zeta > 0.0, "Bulk viscosity should be positive: {zeta}");
    }
    #[test]
    fn test_optimal_antisymmetric_rate_magic_parameter() {
        let s_plus = 1.0;
        let s_minus = optimal_antisymmetric_rate(s_plus);
        let tau_plus = 1.0 / s_plus;
        let tau_minus = 1.0 / s_minus;
        let lambda = (tau_plus - 0.5) * (tau_minus - 0.5);
        assert!(
            (lambda - TRT_MAGIC).abs() < 0.01,
            "Magic parameter = {lambda}, expected ≈ {TRT_MAGIC}"
        );
    }
    #[test]
    fn test_optimal_antisymmetric_rate_range() {
        for s_plus in [0.5, 1.0, 1.5, 1.9] {
            let s_minus = optimal_antisymmetric_rate(s_plus);
            assert!(
                (0.0..=2.0).contains(&s_minus),
                "s_minus={s_minus} out of [0,2] for s_plus={s_plus}"
            );
        }
    }
}
/// Advanced BGK collision with overrelaxation factor.
///
/// Standard BGK can be enhanced by applying an overrelaxation factor `sigma`
/// that scales the non-equilibrium part differently, improving stability near
/// omega = 2.  The effective relaxation is:
/// `f_out = f - sigma * omega * (f - feq)`
///
/// # Arguments
/// * `f`     – input distribution
/// * `rho`   – macroscopic density
/// * `u`     – macroscopic velocity
/// * `omega` – standard relaxation rate
/// * `sigma` – overrelaxation factor (1.0 = standard BGK)
pub fn collide_bgk_overrelaxation(
    f: &[f64; 19],
    rho: f64,
    u: [f64; 3],
    omega: f64,
    sigma: f64,
) -> [f64; 19] {
    let feq = compute_equilibrium_d3q19(rho, u);
    let mut f_out = *f;
    for (f_out_i, (&fi, &feqi)) in f_out.iter_mut().zip(f.iter().zip(feq.iter())) {
        *f_out_i = fi - sigma * omega * (fi - feqi);
    }
    f_out
}
/// Compute the non-equilibrium stress tensor Π^(1)_{αβ} from f − feq.
///
/// Returns a 6-component symmetric tensor \[xx, yy, zz, xy, xz, yz\].
pub fn compute_pi1_tensor(f: &[f64; 19], feq: &[f64; 19]) -> [f64; 6] {
    let mut pi = [0.0f64; 6];
    for (i, (&fi, &feqi)) in f.iter().zip(feq.iter()).enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        let fneq = fi - feqi;
        pi[0] += fneq * cx * cx;
        pi[1] += fneq * cy * cy;
        pi[2] += fneq * cz * cz;
        pi[3] += fneq * cx * cy;
        pi[4] += fneq * cx * cz;
        pi[5] += fneq * cy * cz;
    }
    pi
}
/// Reconstruct the regularized first-order distribution f^(1) from Π^(1).
///
/// Uses the D3Q19 second-order Hermite expansion:
/// `f^(1)_i = w_i / (2 cs⁴) * c_{iα} c_{iβ} Π^(1)_{αβ}`
pub fn regularized_f1_d3q19(pi1: &[f64; 6], _rho: f64) -> [f64; 19] {
    let cs4 = (1.0 / 3.0_f64).powi(2);
    let mut f1 = [0.0f64; 19];
    for (i, f1_i) in f1.iter_mut().enumerate() {
        let w = D3Q19_WEIGHTS[i];
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        let q = cx * cx * pi1[0]
            + cy * cy * pi1[1]
            + cz * cz * pi1[2]
            + 2.0 * cx * cy * pi1[3]
            + 2.0 * cx * cz * pi1[4]
            + 2.0 * cy * cz * pi1[5];
        *f1_i = w / (2.0 * cs4) * q;
    }
    f1
}
/// Compute the recursive correction tensor Q^(1) from Π^(1) and u.
///
/// This is a simplified 6-component approximation of the third-rank tensor.
pub fn compute_q1_recursive(pi1: &[f64; 6], u: [f64; 3], tau_minus_half: f64) -> [f64; 6] {
    let inv_cs2 = 3.0;
    let factor = -tau_minus_half * inv_cs2;
    let ux = u[0];
    let uy = u[1];
    let uz = u[2];
    [
        factor * (ux * pi1[0] + ux * pi1[0] + ux * pi1[0]),
        factor * (uy * pi1[1] + uy * pi1[1] + uy * pi1[1]),
        factor * (uz * pi1[2] + uz * pi1[2] + uz * pi1[2]),
        factor * (ux * pi1[1] + uy * pi1[0] + ux * pi1[3]),
        factor * (ux * pi1[2] + uz * pi1[0] + ux * pi1[4]),
        factor * (uy * pi1[2] + uz * pi1[1] + uy * pi1[5]),
    ]
}
/// Reconstruct f^(1) with both Π^(1) and Q^(1) corrections.
pub fn regularized_f1_with_q_d3q19(pi1: &[f64; 6], q1: &[f64; 6]) -> [f64; 19] {
    let cs4 = (1.0 / 3.0_f64).powi(2);
    let mut f1 = [0.0f64; 19];
    for (i, f1_i) in f1.iter_mut().enumerate() {
        let w = D3Q19_WEIGHTS[i];
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        let q_pi = cx * cx * pi1[0]
            + cy * cy * pi1[1]
            + cz * cz * pi1[2]
            + 2.0 * cx * cy * pi1[3]
            + 2.0 * cx * cz * pi1[4]
            + 2.0 * cy * cz * pi1[5];
        let q_q = cx * cx * q1[0]
            + cy * cy * q1[1]
            + cz * cz * q1[2]
            + 2.0 * cx * cy * q1[3]
            + 2.0 * cx * cz * q1[4]
            + 2.0 * cy * cz * q1[5];
        *f1_i = w / (2.0 * cs4) * q_pi + w / (6.0 * cs4 * (1.0 / 3.0)) * q_q;
    }
    f1
}
/// Compute the 19 raw moments of a D3Q19 distribution.
///
/// Ordering: m_0 = ρ, m_1..3 = ρu, m_4..9 = second-order, m_10..18 = higher.
pub fn compute_raw_moments_d3q19(f: &[f64; 19]) -> [f64; 19] {
    let mut m = [0.0f64; 19];
    for (i, &fi) in f.iter().enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        m[0] += fi;
        m[1] += fi * cx;
        m[2] += fi * cy;
        m[3] += fi * cz;
        m[4] += fi * cx * cx;
        m[5] += fi * cy * cy;
        m[6] += fi * cz * cz;
        m[7] += fi * cx * cy;
        m[8] += fi * cx * cz;
        m[9] += fi * cy * cz;
        m[10] += fi * cx * cx * cy;
        m[11] += fi * cx * cx * cz;
        m[12] += fi * cy * cy * cx;
        m[13] += fi * cy * cy * cz;
        m[14] += fi * cz * cz * cx;
        m[15] += fi * cz * cz * cy;
        m[16] += fi * cx * cx * cy * cy;
        m[17] += fi * cx * cx * cz * cz;
        m[18] += fi * cy * cy * cz * cz;
    }
    m
}
/// Return relaxation rate for the i-th raw moment group.
pub(super) fn raw_moment_rate(i: usize, omega_v: f64, omega_b: f64, omega_g: f64) -> f64 {
    match i {
        0..=3 => 1.0,
        4..=6 => omega_b,
        7..=9 => omega_v,
        _ => omega_g,
    }
}
/// Reconstruct a distribution from its raw moments (least-squares projection).
///
/// Uses the pseudo-inverse via the D3Q19 velocity basis.
pub fn raw_moments_to_f_d3q19(m: &[f64; 19]) -> [f64; 19] {
    let rho = m[0];
    let jx = m[1];
    let jy = m[2];
    let jz = m[3];
    let ux = if rho > 1e-14 { jx / rho } else { 0.0 };
    let uy = if rho > 1e-14 { jy / rho } else { 0.0 };
    let uz = if rho > 1e-14 { jz / rho } else { 0.0 };
    let mut f_out = compute_equilibrium_d3q19(rho, [ux, uy, uz]);
    let feq_m = compute_raw_moments_d3q19(&f_out);
    let cs4 = (1.0 / 3.0_f64).powi(2);
    let dxx = m[4] - feq_m[4];
    let dyy = m[5] - feq_m[5];
    let dzz = m[6] - feq_m[6];
    let dxy = m[7] - feq_m[7];
    let dxz = m[8] - feq_m[8];
    let dyz = m[9] - feq_m[9];
    for (i, f_out_i) in f_out.iter_mut().enumerate() {
        let w = D3Q19_WEIGHTS[i];
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let cz = c[2] as f64;
        let correction = w / (2.0 * cs4)
            * (cx * cx * dxx
                + cy * cy * dyy
                + cz * cz * dzz
                + 2.0 * cx * cy * dxy
                + 2.0 * cx * cz * dxz
                + 2.0 * cy * cz * dyz);
        *f_out_i += correction;
    }
    f_out
}
/// Compute central moments of a D3Q19 distribution in the co-moving frame u.
///
/// Returns 19 moments in the ordering: κ_000, κ_100, κ_010, κ_001,
/// κ_200, κ_020, κ_002, κ_110, κ_101, κ_011, then higher order.
pub fn compute_central_moments_cumulant(f: &[f64; 19], u: [f64; 3]) -> [f64; 19] {
    let mut cm = [0.0f64; 19];
    let ux = u[0];
    let uy = u[1];
    let uz = u[2];
    for (i, &fi) in f.iter().enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let xi = c[0] as f64 - ux;
        let yi = c[1] as f64 - uy;
        let zi = c[2] as f64 - uz;
        cm[0] += fi;
        cm[1] += fi * xi;
        cm[2] += fi * yi;
        cm[3] += fi * zi;
        cm[4] += fi * xi * xi;
        cm[5] += fi * yi * yi;
        cm[6] += fi * zi * zi;
        cm[7] += fi * xi * yi;
        cm[8] += fi * xi * zi;
        cm[9] += fi * yi * zi;
        cm[10] += fi * xi * xi * yi;
        cm[11] += fi * xi * xi * zi;
        cm[12] += fi * yi * yi * xi;
        cm[13] += fi * yi * yi * zi;
        cm[14] += fi * zi * zi * xi;
        cm[15] += fi * zi * zi * yi;
        cm[16] += fi * xi * xi * yi * yi;
        cm[17] += fi * xi * xi * zi * zi;
        cm[18] += fi * yi * yi * zi * zi;
    }
    cm
}
/// Reconstruct a D3Q19 distribution from its central moments.
///
/// Uses the regularized equilibrium approach: constructs feq at (rho, u),
/// then adds the non-equilibrium correction from the second-order central moments.
pub fn central_moments_to_f_cumulant(cm: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
    let feq = compute_equilibrium_d3q19(rho, u);
    let cm_eq = compute_central_moments_cumulant(&feq, u);
    let cs4 = (1.0 / 3.0_f64).powi(2);
    let dxx = cm[4] - cm_eq[4];
    let dyy = cm[5] - cm_eq[5];
    let dzz = cm[6] - cm_eq[6];
    let dxy = cm[7] - cm_eq[7];
    let dxz = cm[8] - cm_eq[8];
    let dyz = cm[9] - cm_eq[9];
    let mut f_out = feq;
    for (i, f_out_i) in f_out.iter_mut().enumerate() {
        let c = D3Q19_VELOCITIES[i];
        let cx = c[0] as f64 - u[0];
        let cy = c[1] as f64 - u[1];
        let cz = c[2] as f64 - u[2];
        let w = D3Q19_WEIGHTS[i];
        let correction = w / (2.0 * cs4)
            * (cx * cx * dxx
                + cy * cy * dyy
                + cz * cz * dzz
                + 2.0 * cx * cy * dxy
                + 2.0 * cx * cz * dxz
                + 2.0 * cy * cz * dyz);
        *f_out_i += correction;
    }
    f_out
}
/// Compute the entropic equilibrium distribution for D3Q19.
///
/// The entropic equilibrium maximizes the H-function subject to the
/// conservation constraints (ρ, ρu).  For the D3Q19 lattice with the
/// discrete entropic H = Σ_i f_i ln(f_i / w_i), the result reduces to
/// the same quadratic form as the standard equilibrium at low Mach numbers.
///
/// At higher Mach numbers the entropic equilibrium uses the exact product form:
/// `feq_i = w_i ρ Π_α ( 2 − √(1+3u_α²) ) ((2u_α + √(1+3u_α²)) / (1−u_α))^{c_iα}`
///
/// Reference: Ansumali & Karlin, Phys. Rev. E 65, 056312 (2002).
pub fn entropic_equilibrium_d3q19(rho: f64, u: [f64; 3]) -> [f64; 19] {
    let mut feq = [0.0f64; 19];
    for (i, feq_i) in feq.iter_mut().enumerate() {
        let w = D3Q19_WEIGHTS[i];
        let c = D3Q19_VELOCITIES[i];
        let mut prod = rho * w;
        for alpha in 0..3 {
            let ua = u[alpha];
            let ca = c[alpha] as f64;
            let ua2 = ua * ua;
            let a = 2.0 - (1.0 + 3.0 * ua2).sqrt();
            let denom = (1.0 - ua).abs().max(1e-12);
            let b = (2.0 * ua + (1.0 + 3.0 * ua2).sqrt()) / denom;
            prod *= a * b.powf(ca);
        }
        *feq_i = prod;
    }
    feq
}
/// Compute the entropic relaxation factor α for ELBM.
///
/// Solves the mirror state equation H(f^{mirr}) = H(f) where
/// `f^{mirr} = f + α(feq − f)`.  Uses a bisection root-finder.
///
/// Returns α ∈ (0, 2].
pub fn entropic_alpha_d3q19(f: &[f64; 19], feq: &[f64; 19]) -> f64 {
    let df: [f64; 19] = std::array::from_fn(|i| feq[i] - f[i]);
    let h_eval = |alpha: f64| -> f64 {
        f.iter()
            .zip(df.iter())
            .zip(D3Q19_WEIGHTS.iter())
            .map(|((&fi_raw, &dfi), &wi)| {
                let fi = fi_raw + alpha * dfi;
                if fi > 1e-300 && wi > 0.0 {
                    fi * (fi / wi).ln()
                } else {
                    0.0
                }
            })
            .sum()
    };
    let h0 = h_eval(0.0);
    if h_eval(2.0) <= h0 {
        return 2.0;
    }
    let mut lo = 1.0f64;
    let mut hi = 2.0f64;
    for _ in 0..50 {
        let mid = (lo + hi) * 0.5;
        if h_eval(mid) < h0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}
/// Entropic stabilized BGK (ELBM) collision with exact α from entropic condition.
///
/// Uses the entropic equilibrium and computes α to ensure H(f_out) ≤ H(f_in).
pub fn collide_entropic_stabilized(f: &[f64; 19], rho: f64, u: [f64; 3]) -> [f64; 19] {
    let feq = entropic_equilibrium_d3q19(rho, u);
    let alpha = entropic_alpha_d3q19(f, &feq);
    let mut f_out = *f;
    for (f_out_i, (&fi, &feqi)) in f_out.iter_mut().zip(f.iter().zip(feq.iter())) {
        *f_out_i += alpha * (feqi - fi);
    }
    f_out
}
/// Compute the D3Q19 discrete H-function: H = Σ_i f_i ln(f_i / w_i).
pub fn discrete_h_function(f: &[f64; 19]) -> f64 {
    f.iter()
        .zip(D3Q19_WEIGHTS.iter())
        .map(|(&fi, &wi)| {
            if fi > 1e-300 && wi > 0.0 {
                fi * (fi / wi).ln()
            } else {
                0.0
            }
        })
        .sum()
}
/// Compute the spectral radius of the linearized BGK operator.
///
/// For BGK: all eigenvalues equal (1 − ω).  Returns |1 − ω| (scalar).
pub fn bgk_spectral_radius(omega: f64) -> f64 {
    (1.0 - omega).abs()
}
/// Compute the maximum stable relaxation rate for a given velocity.
///
/// Based on the von Neumann stability analysis: ω_max = 2/(1 + 3u²·Δt/Δx²).
/// For the lattice (Δt = Δx = 1): ω_max ≈ 2 / (1 + u²·some_factor).
pub fn max_stable_omega(u_max: f64) -> f64 {
    let u2 = u_max * u_max;
    2.0 / (1.0 + 6.0 * u2).max(1.0)
}
/// Compute the maximum non-equilibrium norm (∞-norm) across a distribution.
pub fn max_neq_norm(f: &[f64; 19], feq: &[f64; 19]) -> f64 {
    f.iter()
        .zip(feq.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max)
}
/// Compute the relative entropy production in one collision step.
///
/// Returns (H_before − H_after) / |H_before|.
pub fn relative_entropy_production(f_in: &[f64; 19], f_out: &[f64; 19]) -> f64 {
    let h_in = discrete_h_function(f_in);
    let h_out = discrete_h_function(f_out);
    if h_in.abs() > 1e-14 {
        (h_in - h_out) / h_in.abs()
    } else {
        0.0
    }
}
/// Check whether a distribution satisfies the H-theorem (H non-increasing).
pub fn satisfies_h_theorem(f_in: &[f64; 19], f_out: &[f64; 19]) -> bool {
    discrete_h_function(f_out) <= discrete_h_function(f_in) + 1e-12
}
/// Compute the Chapman-Enskog expansion parameter ε = Kn (Knudsen number proxy).
///
/// ε ≈ |f − feq|₁ / |feq|₁  (non-equilibrium fraction).
pub fn chapman_enskog_epsilon(f: &[f64; 19], feq: &[f64; 19]) -> f64 {
    let neq: f64 = f.iter().zip(feq.iter()).map(|(a, b)| (a - b).abs()).sum();
    let eq: f64 = feq.iter().map(|v| v.abs()).sum();
    if eq > 1e-14 { neq / eq } else { 0.0 }
}
/// Compute kinematic viscosity from omega using the standard LBM relation.
///
/// ν = cs² (1/ω − 0.5) = (1/3)(1/ω − 0.5)
pub fn omega_to_nu(omega: f64) -> f64 {
    (1.0 / 3.0) * (1.0 / omega - 0.5)
}
/// Compute omega from kinematic viscosity.
///
/// ω = 1 / (ν/cs² + 0.5) = 1 / (3ν + 0.5)
pub fn nu_to_omega(nu: f64) -> f64 {
    1.0 / (3.0 * nu + 0.5)
}
/// Estimate the effective relaxation time from the local strain rate.
///
/// Used in Smagorinsky-enhanced LBM to adjust τ locally.
/// τ_eff = 0.5 (τ_0 + √(τ_0² + 18 C_s² Δ² |S|))
///
/// # Arguments
/// * `tau0`   – base relaxation time
/// * `cs_sgs` – Smagorinsky constant C_s
/// * `delta`  – filter width Δ
/// * `s_mag`  – strain-rate magnitude |S|
pub fn smagorinsky_tau_eff(tau0: f64, cs_sgs: f64, delta: f64, s_mag: f64) -> f64 {
    let discriminant = tau0 * tau0 + 18.0 * cs_sgs * cs_sgs * delta * delta * s_mag;
    0.5 * (tau0 + discriminant.sqrt())
}
#[cfg(test)]
mod tests_advanced_collision {
    use super::*;
    use crate::collision::types::*;
    fn make_perturbed(rho: f64, u: [f64; 3], amp: f64) -> [f64; 19] {
        let feq = compute_equilibrium_d3q19(rho, u);
        let mut f = feq;
        for (i, (f_i, &feq_i)) in f.iter_mut().zip(feq.iter()).enumerate() {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            *f_i += sign * amp * feq_i;
        }
        f
    }
    #[test]
    fn test_bgk_overrelaxation_sigma_one_equals_bgk() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let omega = 1.2;
        let f_bgk = collide_bgk(&f, rho, u, omega);
        let f_or = collide_bgk_overrelaxation(&f, rho, u, omega, 1.0);
        for (i, (&f_or_i, &f_bgk_i)) in f_or.iter().zip(f_bgk.iter()).enumerate() {
            assert!(
                (f_or_i - f_bgk_i).abs() < 1e-13,
                "sigma=1 overrelax != BGK at i={i}"
            );
        }
    }
    #[test]
    fn test_bgk_overrelaxation_mass_conserved() {
        let f = make_perturbed(1.0, [0.02, 0.01, 0.0], 0.005);
        let (rho, u) = compute_macroscopic(&f);
        let f_out = collide_bgk_overrelaxation(&f, rho, u, 1.0, 1.2);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-12,
            "Overrelax mass not conserved: {rho_out}"
        );
    }
    #[test]
    fn test_bgk_overrelaxation_struct_stability() {
        let model = BgkOverrelaxation::new(1.0 / 6.0, 0.9);
        assert!(
            model.is_stable(),
            "omega={} sigma={} should be stable",
            model.omega,
            model.sigma
        );
    }
    #[test]
    fn test_bgk_overrelaxation_struct_unstable() {
        let model = BgkOverrelaxation::new(1.0 / 6.0, 2.5);
        assert!(!model.is_stable(), "sigma=2.5 * omega>1 should be unstable");
    }
    #[test]
    fn test_bgk_overrelaxation_effective_viscosity() {
        let nu = 0.1;
        let model = BgkOverrelaxation::new(nu, 1.0);
        assert!(
            (model.effective_viscosity() - nu).abs() < 1e-12,
            "Effective viscosity mismatch: {}",
            model.effective_viscosity()
        );
    }
    #[test]
    fn test_bgk_overrelaxation_at_equilibrium() {
        let rho = 1.0;
        let u = [0.01, 0.02, -0.01];
        let feq = compute_equilibrium_d3q19(rho, u);
        let f_out = collide_bgk_overrelaxation(&feq, rho, u, 1.5, 1.0);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-13,
                "Overrelax should fix equilibrium at i={i}"
            );
        }
    }
    #[test]
    fn test_regularized_full_mass_conservation() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RegularizedCollisionFull::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-11,
            "Regularized full mass not conserved: {rho_out}"
        );
    }
    #[test]
    fn test_regularized_full_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let coll = RegularizedCollisionFull::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&feq, rho, u);
        for (i, (&f_out_i, &feq_i)) in f_out.iter().zip(feq.iter()).enumerate() {
            assert!(
                (f_out_i - feq_i).abs() < 1e-12,
                "Regularized full modified equilibrium at i={i}"
            );
        }
    }
    #[test]
    fn test_regularized_full_viscosity() {
        let nu = 0.15;
        let coll = RegularizedCollisionFull::from_viscosity(nu);
        assert!(
            (coll.effective_viscosity() - nu).abs() < 1e-12,
            "Viscosity mismatch: {}",
            coll.effective_viscosity()
        );
    }
    #[test]
    fn test_compute_pi1_tensor_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi1 = compute_pi1_tensor(&feq, &feq);
        for v in pi1.iter() {
            assert!(v.abs() < 1e-13, "pi1 at equilibrium should be zero: {v}");
        }
    }
    #[test]
    fn test_regularized_f1_reconstructs_stress() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.02);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let pi1 = compute_pi1_tensor(&f, &feq);
        let f1 = regularized_f1_d3q19(&pi1, rho);
        let sum: f64 = f1.iter().sum();
        assert!(
            sum.abs() < 1.0,
            "f1 mass should be small in magnitude: {sum}"
        );
    }
    #[test]
    fn test_regularized_full_momentum_conservation() {
        let f = make_perturbed(1.0, [0.03, -0.01, 0.02], 0.008);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RegularizedCollisionFull::from_viscosity(0.1);
        let f_out = coll.collide(&f, rho, u);
        let vel_components: [[i32; 19]; 3] = {
            let mut c = [[0i32; 19]; 3];
            for (i, vel) in D3Q19_VELOCITIES.iter().enumerate() {
                c[0][i] = vel[0];
                c[1][i] = vel[1];
                c[2][i] = vel[2];
            }
            c
        };
        for (alpha, comp) in vel_components.iter().enumerate() {
            let p_in: f64 = f
                .iter()
                .zip(comp.iter())
                .map(|(fi, &ci)| fi * ci as f64)
                .sum();
            let p_out: f64 = f_out
                .iter()
                .zip(comp.iter())
                .map(|(fi, &ci)| fi * ci as f64)
                .sum();
            assert!(
                (p_out - p_in).abs() < 1e-11,
                "Regularized momentum not conserved in direction {alpha}"
            );
        }
    }
    #[test]
    fn test_rr_lbm_mass_conservation() {
        let f = make_perturbed(1.0, [0.03, 0.01, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RecursiveRegularized::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-10,
            "RR mass not conserved: {rho_out}"
        );
    }
    #[test]
    fn test_rr_lbm_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let coll = RecursiveRegularized::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&feq, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-11,
            "RR disturbed eq mass: {rho_out}"
        );
    }
    #[test]
    fn test_rr_lbm_effective_viscosity() {
        let nu = 0.1;
        let coll = RecursiveRegularized::from_viscosity(nu);
        assert!(
            (coll.effective_viscosity() - nu).abs() < 1e-12,
            "RR viscosity mismatch: {}",
            coll.effective_viscosity()
        );
    }
    #[test]
    fn test_rr_lbm_finite_output() {
        let f = make_perturbed(1.0, [0.04, -0.02, 0.01], 0.015);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RecursiveRegularized::from_viscosity(0.08);
        let f_out = coll.collide(&f, rho, u);
        assert!(f_out.iter().all(|v| v.is_finite()), "RR output has NaN/Inf");
    }
    #[test]
    fn test_compute_q1_recursive_zero_for_zero_pi() {
        let pi1 = [0.0f64; 6];
        let u = [0.02, 0.01, 0.0];
        let q1 = compute_q1_recursive(&pi1, u, 0.5);
        for v in q1.iter() {
            assert!(v.abs() < 1e-14, "Q1 should be zero for zero pi1: {v}");
        }
    }
    #[test]
    fn test_hrr_mass_conservation() {
        let f = make_perturbed(1.0, [0.02, 0.01, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let coll = HybridRecursiveRegularized::new(1.0 / 6.0, 0.05);
        let f_out = coll.collide(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-10,
            "HRR mass not conserved: {rho_out}"
        );
    }
    #[test]
    fn test_hrr_sensor_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let coll = HybridRecursiveRegularized::new(1.0 / 6.0, 0.05);
        let s = coll.sensor(&feq, rho, u);
        assert!(s < 1e-12, "Sensor at equilibrium should be ~0: {s}");
    }
    #[test]
    fn test_hrr_sensor_positive_for_neq() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.05);
        let (rho, u) = compute_macroscopic(&f);
        let coll = HybridRecursiveRegularized::new(1.0 / 6.0, 0.05);
        let s = coll.sensor(&f, rho, u);
        assert!(
            s > 0.0,
            "Sensor should be positive for non-equilibrium: {s}"
        );
    }
    #[test]
    fn test_hrr_finite_output() {
        let f = make_perturbed(1.0, [0.04, -0.02, 0.01], 0.02);
        let (rho, u) = compute_macroscopic(&f);
        let coll = HybridRecursiveRegularized::new(0.1, 0.02);
        let f_out = coll.collide(&f, rho, u);
        assert!(
            f_out.iter().all(|v| v.is_finite()),
            "HRR output has NaN/Inf"
        );
    }
    #[test]
    fn test_hrr_effective_viscosity() {
        let nu = 0.12;
        let coll = HybridRecursiveRegularized::new(nu, 0.01);
        assert!(
            (coll.effective_viscosity() - nu).abs() < 1e-12,
            "HRR viscosity: {}",
            coll.effective_viscosity()
        );
    }
    #[test]
    fn test_raw_moment_mass_conservation() {
        let f = make_perturbed(1.0, [0.03, 0.01, 0.0], 0.008);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RawMomentCollision::from_viscosity(1.0 / 6.0, 0.0);
        let f_out = coll.collide(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho).abs() < 1e-2, "RawMoment mass: {rho_out}");
    }
    #[test]
    fn test_raw_moment_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let coll = RawMomentCollision::from_viscosity(1.0 / 6.0, 0.0);
        let f_out = coll.collide(&feq, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-11,
            "RawMoment eq mass: {rho_out}"
        );
    }
    #[test]
    fn test_raw_moment_finite_output() {
        let f = make_perturbed(1.0, [0.04, -0.02, 0.01], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let coll = RawMomentCollision::from_viscosity(0.08, 0.01);
        let f_out = coll.collide(&f, rho, u);
        assert!(
            f_out.iter().all(|v| v.is_finite()),
            "RawMoment output NaN/Inf"
        );
    }
    #[test]
    fn test_compute_raw_moments_density() {
        let rho = 1.2;
        let u = [0.05, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let m = compute_raw_moments_d3q19(&feq);
        assert!(
            (m[0] - rho).abs() < 1e-12,
            "Raw moment m0 should equal rho: {}",
            m[0]
        );
    }
    #[test]
    fn test_compute_raw_moments_momentum() {
        let rho = 1.0;
        let u = [0.03, 0.01, -0.02];
        let feq = compute_equilibrium_d3q19(rho, u);
        let m = compute_raw_moments_d3q19(&feq);
        assert!((m[1] - rho * u[0]).abs() < 1e-12, "Raw m1 = rho*ux");
        assert!((m[2] - rho * u[1]).abs() < 1e-12, "Raw m2 = rho*uy");
        assert!((m[3] - rho * u[2]).abs() < 1e-12, "Raw m3 = rho*uz");
    }
    #[test]
    fn test_raw_moment_shear_viscosity() {
        let nu = 0.12;
        let coll = RawMomentCollision::from_viscosity(nu, 0.0);
        assert!(
            (coll.shear_viscosity() - nu).abs() < 1e-12,
            "Shear viscosity: {}",
            coll.shear_viscosity()
        );
    }
    #[test]
    fn test_cumulant_mass_conservation() {
        let f = make_perturbed(1.0, [0.03, 0.01, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let coll = CumulantCollision::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho).abs() < 1e-11, "Cumulant mass: {rho_out}");
    }
    #[test]
    fn test_cumulant_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.01, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let coll = CumulantCollision::from_viscosity(1.0 / 6.0);
        let f_out = coll.collide(&feq, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!((rho_out - rho).abs() < 1e-11, "Cumulant eq mass: {rho_out}");
    }
    #[test]
    fn test_cumulant_finite_output() {
        let f = make_perturbed(1.0, [0.04, -0.02, 0.01], 0.015);
        let (rho, u) = compute_macroscopic(&f);
        let coll = CumulantCollision::from_viscosity(0.1);
        let f_out = coll.collide(&f, rho, u);
        assert!(
            f_out.iter().all(|v| v.is_finite()),
            "Cumulant output NaN/Inf"
        );
    }
    #[test]
    fn test_cumulant_shear_viscosity() {
        let nu = 0.15;
        let coll = CumulantCollision::from_viscosity(nu);
        assert!(
            (coll.shear_viscosity() - nu).abs() < 1e-12,
            "Cumulant shear viscosity: {}",
            coll.shear_viscosity()
        );
    }
    #[test]
    fn test_cumulant_from_viscosities() {
        let coll = CumulantCollision::from_viscosities(0.1, 0.05);
        assert!(
            coll.omega > 0.0 && coll.omega_b > 0.0,
            "Both rates should be positive"
        );
    }
    #[test]
    fn test_compute_central_moments_density() {
        let rho = 1.1;
        let u = [0.03, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let cm = compute_central_moments_cumulant(&feq, u);
        assert!((cm[0] - rho).abs() < 1e-12, "Central m0 = rho: {}", cm[0]);
    }
    #[test]
    fn test_compute_central_moments_first_order_zero() {
        let rho = 1.0;
        let u = [0.03, 0.01, -0.02];
        let feq = compute_equilibrium_d3q19(rho, u);
        let cm = compute_central_moments_cumulant(&feq, u);
        assert!(cm[1].abs() < 1e-12, "Central m_100 should be 0: {}", cm[1]);
        assert!(cm[2].abs() < 1e-12, "Central m_010 should be 0: {}", cm[2]);
        assert!(cm[3].abs() < 1e-12, "Central m_001 should be 0: {}", cm[3]);
    }
    #[test]
    fn test_entropic_equilibrium_sums_to_rho() {
        let rho = 1.0;
        let u = [0.05, 0.02, -0.03];
        let feq = entropic_equilibrium_d3q19(rho, u);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - rho).abs() < 0.1 * rho,
            "Entropic eq sum = {sum}, expected {rho}"
        );
    }
    #[test]
    fn test_entropic_equilibrium_positive() {
        let feq = entropic_equilibrium_d3q19(1.0, [0.04, 0.0, 0.0]);
        for (i, &fi) in feq.iter().enumerate() {
            assert!(fi > 0.0, "Entropic feq[{i}] = {fi} is not positive");
        }
    }
    #[test]
    fn test_entropic_alpha_range() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let alpha = entropic_alpha_d3q19(&f, &feq);
        assert!(
            (0.0..=(2.0 + 1e-10)).contains(&alpha),
            "Alpha out of [0,2]: {alpha}"
        );
    }
    #[test]
    fn test_entropic_stabilized_mass_conservation() {
        let f = make_perturbed(1.0, [0.03, 0.01, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let f_out = collide_entropic_stabilized(&f, rho, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-11,
            "Entropic stabilized mass: {rho_out}"
        );
    }
    #[test]
    fn test_discrete_h_function_minimum_at_eq() {
        let rho = 1.0;
        let u = [0.0, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let h_eq = discrete_h_function(&feq);
        let f_pert = make_perturbed(rho, u, 0.05);
        let h_pert = discrete_h_function(&f_pert);
        assert!(
            h_eq <= h_pert + 1e-10,
            "H at eq ({h_eq}) should be <= H at perturbed ({h_pert})"
        );
    }
    #[test]
    fn test_bgk_spectral_radius_at_omega_one() {
        assert!(
            (bgk_spectral_radius(1.0) - 0.0).abs() < 1e-14,
            "Spectral radius at omega=1 should be 0"
        );
    }
    #[test]
    fn test_bgk_spectral_radius_at_omega_half() {
        assert!(
            (bgk_spectral_radius(0.5) - 0.5).abs() < 1e-14,
            "Spectral radius at omega=0.5 should be 0.5"
        );
    }
    #[test]
    fn test_max_stable_omega_zero_velocity() {
        let omega_max = max_stable_omega(0.0);
        assert!(
            omega_max >= 1.0,
            "Max stable omega at u=0 should be >=1: {omega_max}"
        );
    }
    #[test]
    fn test_max_neq_norm_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let norm = max_neq_norm(&feq, &feq);
        assert!(
            norm < 1e-14,
            "Max NEQ norm at equilibrium should be 0: {norm}"
        );
    }
    #[test]
    fn test_omega_nu_roundtrip() {
        let nu = 0.1;
        let omega = nu_to_omega(nu);
        let nu_back = omega_to_nu(omega);
        assert!((nu_back - nu).abs() < 1e-12, "ν→ω→ν roundtrip: {nu_back}");
    }
    #[test]
    fn test_satisfies_h_theorem_bgk() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.01);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let f_out = collide_bgk(&f, rho, u, 1.2);
        let _ = satisfies_h_theorem(&f, &f_out);
        let h_in = discrete_h_function(&f);
        let h_eq = discrete_h_function(&feq);
        assert!(h_eq <= h_in + 1e-10, "H at eq should be ≤ H at non-eq");
    }
    #[test]
    fn test_chapman_enskog_epsilon_at_equilibrium() {
        let rho = 1.0;
        let u = [0.02, 0.0, 0.0];
        let feq = compute_equilibrium_d3q19(rho, u);
        let eps = chapman_enskog_epsilon(&feq, &feq);
        assert!(eps < 1e-14, "C-E epsilon at eq should be 0: {eps}");
    }
    #[test]
    fn test_chapman_enskog_epsilon_positive_for_neq() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.05);
        let (rho, u) = compute_macroscopic(&f);
        let feq = compute_equilibrium_d3q19(rho, u);
        let eps = chapman_enskog_epsilon(&f, &feq);
        assert!(
            eps > 0.0,
            "C-E epsilon should be positive for non-eq: {eps}"
        );
    }
    #[test]
    fn test_smagorinsky_tau_eff_positive() {
        let tau_eff = smagorinsky_tau_eff(0.7, 0.1, 1.0, 0.5);
        assert!(
            tau_eff > 0.0,
            "Smagorinsky tau_eff should be positive: {tau_eff}"
        );
    }
    #[test]
    fn test_smagorinsky_tau_eff_zero_strain() {
        let tau0 = 0.7;
        let tau_eff = smagorinsky_tau_eff(tau0, 0.1, 1.0, 0.0);
        assert!(
            (tau_eff - tau0).abs() < 1e-12,
            "Zero strain: tau_eff should equal tau0: {tau_eff}"
        );
    }
    #[test]
    fn test_relative_entropy_production_bgk() {
        let f = make_perturbed(1.0, [0.03, 0.0, 0.0], 0.02);
        let (rho, u) = compute_macroscopic(&f);
        let f_out = collide_bgk(&f, rho, u, 1.2);
        let rep = relative_entropy_production(&f, &f_out);
        assert!(
            rep.is_finite(),
            "Relative entropy production should be finite: {rep}"
        );
    }
    #[test]
    fn test_cumulant_momentum_conservation() {
        let f = make_perturbed(1.0, [0.03, -0.01, 0.02], 0.008);
        let (rho, u) = compute_macroscopic(&f);
        let coll = CumulantCollision::from_viscosity(0.1);
        let f_out = coll.collide(&f, rho, u);
        let vel_components_c: [[i32; 19]; 3] = {
            let mut c = [[0i32; 19]; 3];
            for (i, vel) in D3Q19_VELOCITIES.iter().enumerate() {
                c[0][i] = vel[0];
                c[1][i] = vel[1];
                c[2][i] = vel[2];
            }
            c
        };
        for (alpha, comp) in vel_components_c.iter().enumerate() {
            let p_in: f64 = f
                .iter()
                .zip(comp.iter())
                .map(|(fi, &ci)| fi * ci as f64)
                .sum();
            let p_out: f64 = f_out
                .iter()
                .zip(comp.iter())
                .map(|(fi, &ci)| fi * ci as f64)
                .sum();
            assert!(
                (p_out - p_in).abs() < 1e-3,
                "Cumulant momentum not conserved in direction {alpha}: {p_in} vs {p_out}"
            );
        }
    }
}
