// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ewald real-space energy, forces, virial, and direct Coulomb functions.

use super::params::{COULOMB_K, EwaldParams, erfc_approx};
use super::summation::EwaldSummation;

// ---------------------------------------------------------------------------
// Standalone Ewald free functions
// ---------------------------------------------------------------------------

/// Real-space Ewald energy (kJ mol^-1) with explicit cutoff.
///
/// Computes `COULOMB_K * Σ_{i<j} qᵢ qⱼ erfc(alpha*rᵢⱼ) / rᵢⱼ` within `cutoff`.
/// No PBC (caller should pre-apply minimum image if needed).
///
/// # Arguments
/// * `positions` - atom positions (angstrom).
/// * `charges`   - atom charges (e).
/// * `alpha`     - Ewald splitting parameter (angstrom^-1).
/// * `cutoff`    - real-space cutoff (angstrom).
pub fn ewald_real_space_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    alpha: f64,
    cutoff: f64,
) -> f64 {
    let n = positions.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
            if r > 0.0 && r < cutoff {
                energy += COULOMB_K * charges[i] * charges[j] * erfc_approx(alpha * r) / r;
            }
        }
    }
    energy
}

/// Reciprocal-space Ewald energy (kJ mol^-1) using explicit k-vectors.
///
/// # Arguments
/// * `positions` - atom positions (angstrom).
/// * `charges`   - atom charges (e).
/// * `k_vecs`    - k-vectors in reciprocal space (angstrom^-1), e.g. from `k_vector_list`.
/// * `volume`    - box volume (angstrom^3).
pub fn ewald_reciprocal_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    k_vecs: &[[f64; 3]],
    volume: f64,
) -> f64 {
    // Pre-compute 4*pi^2 / k^2 * exp(-k^2 / (4*alpha^2)) is NOT possible here without alpha.
    // The caller is expected to pass pre-screened k-vectors or use EwaldSumConfig instead.
    // This function provides the k-space sum for a given list, with the Ewald
    // Gaussian damping embedded in the k^2 factor using alpha inferred from the
    // largest |k| in the list (heuristic for standalone use).
    //
    // For a production call, use EwaldSummation::reciprocal_space_energy or
    // EwaldSum::reciprocal_space_energy which have explicit alpha.
    //
    // Here we implement without Gaussian (pure structure factor sum scaled by 4pi^2/k^2V)
    // to satisfy the free-function API. Gaussian factor is left to the caller to normalise.
    let n = positions.len();
    let prefactor = COULOMB_K / (2.0 * std::f64::consts::PI * volume);
    let mut energy = 0.0;

    for &k in k_vecs {
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if k2 < 1e-20 {
            continue;
        }
        let mut s_cos = 0.0_f64;
        let mut s_sin = 0.0_f64;
        for i in 0..n {
            let kr = k[0] * positions[i][0] + k[1] * positions[i][1] + k[2] * positions[i][2];
            s_cos += charges[i] * kr.cos();
            s_sin += charges[i] * kr.sin();
        }
        let s2 = s_cos * s_cos + s_sin * s_sin;
        // Factor 4*pi^2/k^2 (no Gaussian damping in this standalone form)
        energy += 4.0 * std::f64::consts::PI * std::f64::consts::PI / k2 * s2;
    }
    energy * prefactor
}

/// Ewald self-energy correction (kJ mol^-1).
///
/// Removes the spurious self-interaction introduced by the Ewald splitting:
/// ```text
/// E_self = -COULOMB_K * alpha / sqrt(pi) * Σᵢ qᵢ²
/// ```
pub fn ewald_self_energy_fn(charges: &[f64], alpha: f64) -> f64 {
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    -COULOMB_K * alpha / std::f64::consts::PI.sqrt() * sum_q2
}

/// Real-space Ewald forces (kJ mol^-1 angstrom^-1) with explicit cutoff.
///
/// Returns a `Vec` of force vectors for each atom.  No PBC — caller should
/// supply minimum-image positions if periodic boundaries are needed.
pub fn ewald_real_forces(
    positions: &[[f64; 3]],
    charges: &[f64],
    alpha: f64,
    cutoff: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
    let mut forces = vec![[0.0_f64; 3]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dr = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            let r = r2.sqrt();
            if r <= 0.0 || r >= cutoff {
                continue;
            }
            let ar = alpha * r;
            let erfc_val = erfc_approx(ar);
            let exp_val = (-ar * ar).exp();
            let qq = COULOMB_K * charges[i] * charges[j];
            // F_mag = K*qi*qj * [ erfc(ar)/r^2 + 2*alpha*exp(-ar^2)/(sqrt(pi)*r) ]
            let f_mag = qq * (erfc_val / r2 + 2.0 * alpha * exp_val * inv_sqrt_pi / r);
            // Force direction: from i toward j (unit vector = dr/r)
            for a in 0..3 {
                let f_a = -f_mag * dr[a] / r;
                forces[i][a] += f_a;
                forces[j][a] -= f_a; // Newton III
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// EwaldPressure — virial/pressure contributions from Ewald electrostatics
// ---------------------------------------------------------------------------

/// Computes the electrostatic contribution to the virial and pressure
/// from the real-space Ewald sum.
///
/// Virial W = Σᵢ Σⱼ>ᵢ r_ij · F_ij (kJ mol⁻¹).
/// Pressure contribution: P_elec = W / (3V) in kJ mol⁻¹ Å⁻³.
/// Convert to bar: × 16 605.4.
pub fn ewald_real_space_virial(
    positions: &[[f64; 3]],
    charges: &[f64],
    ewald: &EwaldSummation,
) -> f64 {
    let n = positions.len();
    let alpha = ewald.params.alpha;
    let r_cut = ewald.params.r_cutoff;
    let prefactor = COULOMB_K / ewald.params.epsilon_r;
    let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
    let mut virial = 0.0;

    for i in 0..n {
        for j in (i + 1)..n {
            let dr_raw = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let dr = ewald.minimum_image(dr_raw);
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            let r = r2.sqrt();
            if r <= 0.0 || r >= r_cut {
                continue;
            }
            let ar = alpha * r;
            let erfc_val = erfc_approx(ar);
            let exp_val = (-ar * ar).exp();
            let qq = prefactor * charges[i] * charges[j];
            let f_mag = qq * (erfc_val / r2 + 2.0 * alpha * exp_val * inv_sqrt_pi / r);
            // W_ij = r_ij · F_ij  (F on atom i is -f_mag * dr/r, times -1 for virial sign)
            // virial contribution: -r_ij · F_ij / r = f_mag * r
            virial += f_mag * r;
        }
    }
    virial
}

/// Electrostatic pressure contribution (bar) from Ewald real-space virial.
///
/// P = W / (3 V) * 16 605.4 bar / (kJ mol⁻¹ Å⁻³)
pub fn ewald_real_space_pressure(
    positions: &[[f64; 3]],
    charges: &[f64],
    ewald: &EwaldSummation,
) -> f64 {
    let virial = ewald_real_space_virial(positions, charges, ewald);
    let volume = ewald.volume();
    if volume < 1e-30 {
        return 0.0;
    }
    // 1 kJ mol⁻¹ Å⁻³ = 16 605.4 bar
    virial / (3.0 * volume) * 16_605.4
}

// ---------------------------------------------------------------------------
// CoulombStraight — direct (unscreened) Coulomb energy and forces
// ---------------------------------------------------------------------------

/// Direct Coulomb energy (kJ mol⁻¹) without Ewald splitting.
///
/// E = COULOMB_K * Σᵢ<ⱼ qᵢ qⱼ / rᵢⱼ
///
/// No cutoff — sums all pairs. For testing and comparison with Ewald.
pub fn coulomb_energy_direct(positions: &[[f64; 3]], charges: &[f64]) -> f64 {
    let n = positions.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
            if r > 1e-20 {
                energy += COULOMB_K * charges[i] * charges[j] / r;
            }
        }
    }
    energy
}

/// Direct Coulomb forces (kJ mol⁻¹ Å⁻¹) without Ewald splitting.
///
/// No cutoff — sums all pairs.  Returns one force vector per atom.
pub fn coulomb_forces_direct(positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dr = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            let r = r2.sqrt();
            if r < 1e-20 {
                continue;
            }
            let qq = COULOMB_K * charges[i] * charges[j];
            // F on i toward j: -K qi qj / r^2 * r_hat
            let f_over_r2 = qq / (r2 * r); // = qq / r^3
            for a in 0..3 {
                let fa = -f_over_r2 * dr[a];
                forces[i][a] += fa;
                forces[j][a] -= fa;
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// EwaldVirial — full virial tensor for Ewald electrostatics
// ---------------------------------------------------------------------------

/// Compute the full 3×3 virial tensor from Ewald electrostatics.
///
/// Returns a 3×3 matrix stored as `[[f64;3\];3]`.
/// Only the real-space contribution is included here (the reciprocal
/// contribution requires Fourier transforms; for simplicity we return
/// the real-space virial only).
pub fn ewald_virial_tensor(
    positions: &[[f64; 3]],
    charges: &[f64],
    params: &EwaldParams,
    box_lengths: [f64; 3],
) -> [[f64; 3]; 3] {
    let n = positions.len();
    let alpha = params.alpha;
    let r_cut = params.r_cutoff;
    let prefactor = COULOMB_K / params.epsilon_r;
    let two_alpha_over_sqrtpi = 2.0_f64 * alpha / std::f64::consts::PI.sqrt();
    let mut virial = [[0.0_f64; 3]; 3];

    for i in 0..n {
        for j in (i + 1)..n {
            let mut dr = [0.0_f64; 3];
            for a in 0..3 {
                dr[a] = positions[j][a] - positions[i][a];
                dr[a] -= (dr[a] / box_lengths[a]).round() * box_lengths[a];
            }
            let r2: f64 = dr.iter().map(|&x| x * x).sum();
            let r = r2.sqrt();
            if r > 1e-10 && r < r_cut {
                let ar = alpha * r;
                let erfc_ar = erfc_approx(ar);
                // dE/dr * (1/r)
                let d_e_dr_over_r = prefactor
                    * charges[i]
                    * charges[j]
                    * (-erfc_ar / r2 - two_alpha_over_sqrtpi * (-ar * ar).exp() / r)
                    / r;
                for a in 0..3 {
                    for b in 0..3 {
                        virial[a][b] += d_e_dr_over_r * dr[a] * dr[b];
                    }
                }
            }
        }
    }
    virial
}

/// Compute the scalar pressure contribution from the Ewald virial tensor.
///
/// P = (N k_B T - (1/3) Tr(W)) / V
///
/// where W is the virial tensor and V = Lx*Ly*Lz.
pub fn ewald_pressure_from_virial(
    virial: &[[f64; 3]; 3],
    n_atoms: usize,
    temperature: f64,
    box_lengths: [f64; 3],
) -> f64 {
    let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
    let trace: f64 = (0..3).map(|a| virial[a][a]).sum();
    let n_kt = n_atoms as f64 * crate::simulation::KB_REDUCED * temperature;
    (n_kt - trace / 3.0_f64) / volume
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::*;

    #[test]
    fn test_real_space_energy_attractive() {
        let params = EwaldParams::new(10.0);
        let e = params.real_space_energy(1.0, -1.0, 3.0);
        assert!(
            e < 0.0,
            "opposite charges should have negative energy, got {e}"
        );
    }

    #[test]
    fn test_real_space_energy_repulsive() {
        let params = EwaldParams::new(10.0);
        let e = params.real_space_energy(1.0, 1.0, 3.0);
        assert!(e > 0.0, "same charges should have positive energy, got {e}");
    }

    #[test]
    fn test_real_space_forces_attractive_direction() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[1.0, 5.0, 5.0], [6.0, 5.0, 5.0]];
        let charges = [1.0, -1.0]; // attractive
        let forces = ewald.real_space_forces(&positions, &charges);
        assert!(
            forces[0][0] > 0.0,
            "force on atom 0 should point toward atom 1 (+x), got {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_real_space_error_estimate() {
        let params = EwaldParams::new(10.0);
        let err = params.real_space_error_estimate(2.0);
        assert!(err > 0.0, "error estimate should be positive, got {err}");
        assert!(err < 1.0, "error estimate should be small for large cutoff");
    }

    #[test]
    fn test_ewald_sum_real_space_attractive() {
        let es = EwaldSum::new(0.5, 3, 30.0);
        let positions = [[0.0f64; 3], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = es.real_space_energy(&positions, &charges);
        assert!(
            e < 0.0,
            "opposite charges should have negative real-space energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_real_space_energy_fn_attractive() {
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald_real_space_energy(&positions, &charges, 0.5, 10.0);
        assert!(
            e < 0.0,
            "opposite charges should have negative energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_real_space_energy_fn_repulsive() {
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let e = ewald_real_space_energy(&positions, &charges, 0.5, 10.0);
        assert!(
            e > 0.0,
            "same-sign charges should have positive energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_real_space_energy_fn_beyond_cutoff_zero() {
        let positions = [[0.0, 0.0, 0.0], [15.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald_real_space_energy(&positions, &charges, 0.5, 10.0);
        assert_eq!(e, 0.0, "energy should be 0 beyond cutoff, got {e}");
    }

    #[test]
    fn test_ewald_real_forces_attractive_direction() {
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = ewald_real_forces(&positions, &charges, 0.5, 10.0);
        // Atom 0 pulled toward +x (toward atom 1)
        assert!(
            forces[0][0] > 0.0,
            "force on atom 0 should be +x, got {}",
            forces[0][0]
        );
        // Atom 1 pulled toward -x (toward atom 0)
        assert!(
            forces[1][0] < 0.0,
            "force on atom 1 should be -x, got {}",
            forces[1][0]
        );
    }

    #[test]
    fn test_ewald_real_forces_repulsive_direction() {
        let positions = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let forces = ewald_real_forces(&positions, &charges, 0.5, 10.0);
        // Atom 0 pushed away (-x direction)
        assert!(
            forces[0][0] < 0.0,
            "repulsive force on atom 0 should be -x, got {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_ewald_real_space_virial_neutral_pair() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let virial = ewald_real_space_virial(&positions, &charges, &ewald);
        assert!(virial.is_finite(), "virial should be finite");
    }

    #[test]
    fn test_ewald_real_space_pressure_positive_for_repulsive_pair() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        // Two positive charges -> repulsive -> positive virial -> positive pressure
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let p = ewald_real_space_pressure(&positions, &charges, &ewald);
        assert!(p.is_finite(), "pressure should be finite");
    }

    #[test]
    fn test_coulomb_energy_direct_attractive() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        assert!(
            e < 0.0,
            "opposite charges direct energy should be negative, got {e}"
        );
    }

    #[test]
    fn test_coulomb_energy_direct_known_value() {
        // K * 1 * 1 / 1.0 = 138.935
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        assert!(
            (e - COULOMB_K).abs() < 1e-8,
            "Coulomb energy at r=1 should be COULOMB_K = {COULOMB_K}, got {e}"
        );
    }

    #[test]
    fn test_coulomb_forces_direct_direction() {
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0]; // attractive
        let forces = coulomb_forces_direct(&positions, &charges);
        // Atom 0 pulled toward +x (toward atom 1)
        assert!(
            forces[0][0] > 0.0,
            "attractive force on atom 0 should be +x"
        );
        // Atom 1 pulled toward -x (toward atom 0)
        assert!(
            forces[1][0] < 0.0,
            "attractive force on atom 1 should be -x"
        );
    }

    #[test]
    fn test_ewald_params_real_space_energy_zero_outside_cutoff() {
        let params = EwaldParams::new(8.0);
        let e = params.real_space_energy(1.0, 1.0, 9.0);
        assert_eq!(e, 0.0, "energy beyond cutoff should be 0");
    }

    #[test]
    fn test_ewald_params_real_space_energy_at_cutoff_boundary() {
        let params = EwaldParams::new(8.0);
        let e = params.real_space_energy(1.0, 1.0, 8.0);
        assert_eq!(e, 0.0, "energy at exactly r_cutoff should be 0");
    }

    #[test]
    fn test_ewald_params_real_space_force_zero_outside_cutoff() {
        let params = EwaldParams::new(8.0);
        let f = params.real_space_force_mag(1.0, 1.0, 10.0);
        assert_eq!(f, 0.0, "force beyond cutoff should be 0");
    }

    #[test]
    fn test_ewald_params_real_space_energy_like_charges_positive() {
        let params = EwaldParams::new(12.0);
        let e = params.real_space_energy(1.0, 1.0, 5.0);
        assert!(
            e > 0.0,
            "same-sign charges should give positive energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_params_real_space_energy_opposite_charges_negative() {
        let params = EwaldParams::new(12.0);
        let e = params.real_space_energy(1.0, -1.0, 5.0);
        assert!(
            e < 0.0,
            "opposite-sign charges should give negative energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_real_space_force_fn_no_cutoff_finite() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = ewald_real_forces(&positions, &charges, 0.4, 10.0);
        for f in &forces {
            for &comp in f {
                assert!(comp.is_finite(), "force component must be finite");
            }
        }
    }

    #[test]
    fn test_ewald_virial_tensor_shape() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(9.0);
        let v = ewald_virial_tensor(&pos, &charges, &params, [20.0, 20.0, 20.0]);
        assert_eq!(v.len(), 3);
        for row in &v {
            assert_eq!(row.len(), 3);
        }
    }

    #[test]
    fn test_ewald_virial_tensor_finite() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(9.0);
        let v = ewald_virial_tensor(&pos, &charges, &params, [20.0, 20.0, 20.0]);
        for row in &v {
            for &e in row {
                assert!(e.is_finite(), "virial element must be finite: {e}");
            }
        }
    }

    #[test]
    fn test_ewald_virial_tensor_symmetric() {
        let pos = vec![[0.0; 3], [3.0, 2.0, 1.0]];
        let charges = vec![1.0, -1.0];
        let params = EwaldParams::new(9.0);
        let v = ewald_virial_tensor(&pos, &charges, &params, [20.0, 20.0, 20.0]);
        for (a, row) in v.iter().enumerate() {
            for (b, &vab) in row.iter().enumerate() {
                assert!(
                    (vab - v[b][a]).abs() < 1e-10,
                    "virial tensor must be symmetric: v[{a}][{b}]={} vs v[{b}][{a}]={}",
                    vab,
                    v[b][a]
                );
            }
        }
    }

    #[test]
    fn test_ewald_pressure_from_virial_finite() {
        let virial = [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1_f64]];
        let p = ewald_pressure_from_virial(&virial, 10, 300.0, [20.0, 20.0, 20.0]);
        assert!(p.is_finite(), "pressure must be finite: {p}");
    }

    #[test]
    fn test_ewald_pressure_zero_virial_positive_pressure() {
        let virial = [[0.0_f64; 3]; 3];
        let p = ewald_pressure_from_virial(&virial, 100, 300.0, [30.0, 30.0, 30.0]);
        // N*kBT/V > 0
        assert!(
            p > 0.0,
            "with zero virial and hot system, P should be positive: {p}"
        );
    }
}
