//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::types::PressureCouplingMode;

/// Compute pressure from kinetic energy and virial.
///
/// Uses the virial theorem:  P = (2·KE + W) / (3·V)
///
/// # Arguments
/// * `volume` – simulation-box volume.
/// * `kinetic_energy` – total kinetic energy of the system.
/// * `virial` – scalar virial W = Σ rᵢ · fᵢ.
pub fn virial_pressure(volume: f64, kinetic_energy: f64, virial: f64) -> f64 {
    (2.0 * kinetic_energy + virial) / (3.0 * volume)
}
/// Compute the scalar virial W = Σ rᵢ · fᵢ from atomic positions and forces.
///
/// # Arguments
/// * `positions` – slice of `[x, y, z]` position arrays.
/// * `forces`    – slice of `[fx, fy, fz]` force arrays (same length).
pub fn compute_virial(positions: &[[f64; 3]], forces: &[[f64; 3]]) -> f64 {
    positions
        .iter()
        .zip(forces.iter())
        .map(|(r, f)| r[0] * f[0] + r[1] * f[1] + r[2] * f[2])
        .sum()
}
/// Trait for barostats that control system pressure.
pub trait Barostat: Send + Sync {
    /// Apply pressure control, potentially rescaling positions and box dimensions.
    fn apply(&mut self, atoms: &mut AtomSet, pbox: &mut PeriodicBox, target_pressure: f64, dt: f64);
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::BerendsenBarostat;
    use crate::ParrinelloRahmanBarostat;
    use crate::barostat::AnisotropicBerendsenBarostat;
    use crate::barostat::MonteCarloBarostat;
    use crate::barostat::MtkBarostat;
    use crate::barostat::VolumeTracker;
    use oxiphysics_core::Vec3;
    #[test]
    fn test_no_barostat_does_nothing() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0, 0.0, 0);
        let mut pbox = PeriodicBox::cubic(10.0);
        let orig_dims = pbox.dims;
        let orig_pos = atoms.positions[0];
        let mut baro = NoBarostat::new();
        baro.apply(&mut atoms, &mut pbox, 1.0, 0.001);
        assert_eq!(pbox.dims, orig_dims);
        assert_eq!(atoms.positions[0], orig_pos);
    }
    #[test]
    fn test_berendsen_barostat_compresses() {
        let box_l = 10.0_f64;
        let mut pbox = PeriodicBox::cubic(box_l);
        let mut atoms = AtomSet::new();
        let n = 20_usize;
        for i in 0..n {
            let x = (i as f64) * (box_l / n as f64);
            let sign = if i % 2 == 0 { 1.0_f64 } else { -1.0_f64 };
            atoms.add_atom(
                Vec3::new(x, x % box_l, 0.0),
                Vec3::new(sign * 5.0, 0.0, 0.0),
                1.0,
                0.0,
                0,
            );
        }
        let ke = atoms.kinetic_energy();
        let vol_init = pbox.volume();
        let p_kinetic = 2.0 * ke / (3.0 * vol_init);
        let target_pressure = p_kinetic * 10.0;
        let mut baro = BerendsenBarostat::new(target_pressure, 1.0, 0.01);
        for _ in 0..100 {
            Barostat::apply(&mut baro, &mut atoms, &mut pbox, target_pressure, 0.01);
        }
        let vol_final = pbox.volume();
        assert!(
            vol_final < vol_init,
            "Berendsen barostat should compress: V_init={vol_init:.3}, V_final={vol_final:.3}"
        );
    }
    #[test]
    fn test_virial_pressure_ideal_gas() {
        let ke = 150.0;
        let vol = 1000.0;
        let p = virial_pressure(vol, ke, 0.0);
        let expected = 2.0 * ke / (3.0 * vol);
        assert!(
            (p - expected).abs() < 1e-12,
            "virial_pressure ideal gas: got {p}, expected {expected}"
        );
    }
    #[test]
    fn test_compute_virial_equal_opposite() {
        let positions = [[1.0_f64, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let forces = [[2.0_f64, 0.0, 0.0], [-2.0, 0.0, 0.0]];
        let w = compute_virial(&positions, &forces);
        assert!(
            (w - (-4.0)).abs() < 1e-12,
            "compute_virial: got {w}, expected -4"
        );
    }
    #[test]
    fn test_berendsen_scale_factor_over_pressure() {
        let baro = BerendsenBarostat::new(1.0, 1.0, 0.01);
        let mu_under = baro.scale_factor(0.5, 0.01);
        assert!(
            mu_under < 1.0,
            "under-pressure should shrink box: mu={mu_under}"
        );
        let mu_exact = baro.scale_factor(1.0, 0.01);
        assert!(
            (mu_exact - 1.0).abs() < 1e-12,
            "at target pressure scale_factor should be 1.0, got {mu_exact}"
        );
    }
    #[test]
    fn test_mc_barostat_zero_delta() {
        let baro = MonteCarloBarostat::new(1.0);
        let log_p = baro.acceptance_log_prob(0.0, 0.0, 100, 300.0, 1.380649e-23);
        assert!(
            log_p.abs() < 1e-12,
            "ΔE=0, ΔV=0 should give log_prob=0, got {log_p}"
        );
    }
    #[test]
    fn test_pressure_tensor_scalar_ideal_gas() {
        let velocities = [[1.0, 0.0, 0.0]];
        let masses = [1.0];
        let positions = [[0.0, 0.0, 0.0]];
        let forces = [[0.0, 0.0, 0.0]];
        let volume = 10.0;
        let pt = PressureTensor::compute(volume, &velocities, &masses, &positions, &forces);
        let scalar = pt.scalar_pressure();
        assert!(
            (scalar - 1.0 / 30.0).abs() < 1e-12,
            "scalar pressure: got {scalar}, expected {}",
            1.0 / 30.0
        );
    }
    #[test]
    fn test_pressure_tensor_diagonal() {
        let velocities = [[1.0, 2.0, 3.0]];
        let masses = [2.0];
        let positions = [[0.0, 0.0, 0.0]];
        let forces = [[0.0, 0.0, 0.0]];
        let volume = 1.0;
        let pt = PressureTensor::compute(volume, &velocities, &masses, &positions, &forces);
        let d = pt.diagonal();
        assert!((d[0] - 2.0).abs() < 1e-12);
        assert!((d[1] - 8.0).abs() < 1e-12);
        assert!((d[2] - 18.0).abs() < 1e-12);
    }
    #[test]
    fn test_pressure_tensor_anisotropy() {
        let velocities = [[1.0, 2.0, 3.0]];
        let masses = [1.0];
        let positions = [[0.0, 0.0, 0.0]];
        let forces = [[0.0, 0.0, 0.0]];
        let volume = 1.0;
        let pt = PressureTensor::compute(volume, &velocities, &masses, &positions, &forces);
        let a = pt.anisotropy();
        assert!((a - 8.0).abs() < 1e-12, "anisotropy: got {a}, expected 8");
    }
    #[test]
    fn test_pressure_tensor_off_diagonal() {
        let velocities = [[1.0, 0.0, 0.0]];
        let masses = [1.0];
        let positions = [[1.0, 0.0, 0.0]];
        let forces = [[0.0, 1.0, 0.0]];
        let volume = 1.0;
        let pt = PressureTensor::compute(volume, &velocities, &masses, &positions, &forces);
        let od = pt.off_diagonal();
        assert!((od[0] - 1.0).abs() < 1e-12, "P_xy: got {}", od[0]);
    }
    #[test]
    fn test_virial_tensor_two_atoms() {
        let positions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let forces = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let w = PressureTensor::virial_tensor(&positions, &forces);
        assert!((w[0][0] - 2.0).abs() < 1e-12);
        assert!((w[1][1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_volume_tracker_empty() {
        let vt = VolumeTracker::new(100);
        assert_eq!(vt.count, 0);
        assert!((vt.mean() - 0.0).abs() < 1e-15);
        assert!((vt.variance() - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_volume_tracker_mean() {
        let mut vt = VolumeTracker::new(100);
        vt.record(10.0);
        vt.record(20.0);
        vt.record(30.0);
        assert!((vt.mean() - 20.0).abs() < 1e-12, "mean: got {}", vt.mean());
    }
    #[test]
    fn test_volume_tracker_variance() {
        let mut vt = VolumeTracker::new(100);
        for _ in 0..10 {
            vt.record(5.0);
        }
        assert!(
            vt.variance().abs() < 1e-12,
            "var for constant: {}",
            vt.variance()
        );
    }
    #[test]
    fn test_volume_tracker_min_max() {
        let mut vt = VolumeTracker::new(100);
        vt.record(100.0);
        vt.record(50.0);
        vt.record(200.0);
        assert!((vt.min_volume - 50.0).abs() < 1e-12);
        assert!((vt.max_volume - 200.0).abs() < 1e-12);
    }
    #[test]
    fn test_volume_tracker_relative_fluctuation() {
        let mut vt = VolumeTracker::new(100);
        for _ in 0..5 {
            vt.record(100.0);
        }
        assert!(vt.relative_fluctuation().abs() < 1e-12);
    }
    #[test]
    fn test_volume_tracker_reset() {
        let mut vt = VolumeTracker::new(100);
        vt.record(10.0);
        vt.record(20.0);
        vt.reset();
        assert_eq!(vt.count, 0);
        assert!(vt.history().is_empty());
    }
    #[test]
    fn test_volume_tracker_compressibility_estimate() {
        let mut vt = VolumeTracker::new(100);
        vt.record(100.0);
        vt.record(110.0);
        let kappa = vt.compressibility_estimate(1.0);
        assert!(kappa > 0.0, "compressibility should be positive");
    }
    #[test]
    fn test_volume_tracker_history_bounded() {
        let mut vt = VolumeTracker::new(3);
        for i in 0..10 {
            vt.record(i as f64);
        }
        assert_eq!(vt.history().len(), 3);
        assert_eq!(vt.count, 10);
    }
    #[test]
    fn test_anisotropic_berendsen_uniform() {
        let baro = AnisotropicBerendsenBarostat::new(1.0, 1.0, 0.01);
        let mu = baro.scale_factors([1.0, 1.0, 1.0], 0.01);
        for k in 0..3 {
            assert!(
                (mu[k] - 1.0).abs() < 1e-12,
                "axis {k}: mu={}, expected 1.0",
                mu[k]
            );
        }
    }
    #[test]
    fn test_anisotropic_berendsen_different_axes() {
        let baro =
            AnisotropicBerendsenBarostat::anisotropic([1.0, 2.0, 3.0], 1.0, [0.01, 0.01, 0.01]);
        let mu = baro.scale_factors([1.0, 2.0, 3.0], 0.01);
        for k in 0..3 {
            assert!((mu[k] - 1.0).abs() < 1e-12, "axis {k}: mu={}", mu[k]);
        }
    }
    #[test]
    fn test_anisotropic_berendsen_apply() {
        let baro =
            AnisotropicBerendsenBarostat::anisotropic([2.0, 2.0, 2.0], 1.0, [0.01, 0.01, 0.01]);
        let mut positions = [[1.0, 1.0, 1.0]];
        let mut box_len = [10.0, 10.0, 10.0];
        baro.apply_anisotropic(&mut positions, &mut box_len, [1.0, 1.0, 1.0], 0.01);
        for k in 0..3 {
            assert!(
                box_len[k] < 10.0,
                "axis {k}: box should shrink, got {}",
                box_len[k]
            );
        }
    }
    #[test]
    fn test_mtk_barostat_creation() {
        let mtk = MtkBarostat::new(1.0, 100.0, 10.0, 300.0, 100);
        assert!((mtk.epsilon_dot - 0.0).abs() < 1e-15);
        assert!((mtk.xi_baro - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_mtk_barostat_force_at_target() {
        let mtk = MtkBarostat::new(1.0, 100.0, 10.0, 300.0, 100);
        let g = mtk.barostat_force(1000.0, 1.0, 0.0);
        assert!(g.abs() < 1e-12, "force at target: {g}");
    }
    #[test]
    fn test_mtk_barostat_kinetic_energy() {
        let mut mtk = MtkBarostat::new(1.0, 100.0, 10.0, 300.0, 100);
        mtk.epsilon_dot = 0.5;
        let ke = mtk.kinetic_energy();
        assert!((ke - 12.5).abs() < 1e-12, "KE: {ke}");
    }
    #[test]
    fn test_mtk_scale_factor_zero_velocity() {
        let mtk = MtkBarostat::new(1.0, 100.0, 10.0, 300.0, 100);
        let s = mtk.scale_factor(0.001);
        assert!((s - 1.0).abs() < 1e-12, "scale at zero vel: {s}");
    }
    #[test]
    fn test_mtk_barostat_step() {
        let mut mtk = MtkBarostat::new(1.0, 100.0, 10.0, 300.0, 10);
        let mut positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut box_lengths = [10.0, 10.0, 10.0];
        mtk.apply_step(&mut positions, &mut box_lengths, 1000.0, 2.0, 100.0, 0.001);
        let changed = box_lengths[0] != 10.0 || box_lengths[1] != 10.0 || box_lengths[2] != 10.0;
        assert!(changed, "MTK step should change box dimensions");
    }
    #[test]
    fn test_mc_barostat_propose_at_half() {
        let baro = MonteCarloBarostat::new(1.0);
        let v_new = baro.propose_volume_change(100.0, 0.5);
        assert!(
            (v_new - 100.0).abs() < 1e-12,
            "at rng=0.5 volume should be unchanged, got {v_new}"
        );
    }
    #[test]
    fn test_mc_barostat_propose_expand() {
        let baro = MonteCarloBarostat::new(1.0);
        let v_new = baro.propose_volume_change(100.0, 0.9);
        assert!(v_new > 100.0, "rng>0.5 should expand, got {v_new}");
    }
    #[test]
    fn test_mc_barostat_propose_compress() {
        let baro = MonteCarloBarostat::new(1.0);
        let v_new = baro.propose_volume_change(100.0, 0.1);
        assert!(v_new < 100.0, "rng<0.5 should compress, got {v_new}");
    }
    #[test]
    fn test_pr_box_velocity_at_target() {
        let pr = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
        let v_new = pr.box_velocity_update(0.0, 1.0, 1000.0, 0.001);
        assert!(v_new.abs() < 1e-12, "at target: v={v_new}");
    }
    #[test]
    fn test_pr_box_velocity_above_target() {
        let pr = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
        let v_new = pr.box_velocity_update(0.0, 2.0, 1000.0, 0.001);
        assert!(v_new > 0.0, "above target: v={v_new}");
    }
}
/// Compute the 3×3 virial tensor from pair interactions.
///
/// W_ab = Σ_{i<j} r_{ij,a} * f_{ij,b}
///
/// where r_ij is the vector from j to i, and f_ij is the force on i from j.
///
/// # Arguments
/// * `pair_rij`   – `[dx, dy, dz]` for each pair (i to j).
/// * `pair_forces`– force on atom i from atom j for each pair.
pub fn virial_tensor_from_pairs(pair_rij: &[[f64; 3]], pair_forces: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let mut w = [[0.0f64; 3]; 3];
    for (r, f) in pair_rij.iter().zip(pair_forces.iter()) {
        for a in 0..3 {
            for b in 0..3 {
                w[a][b] += r[a] * f[b];
            }
        }
    }
    w
}
/// Compute pressure tensor diagonal components (P_xx, P_yy, P_zz) from
/// a known pressure tensor.
///
/// Returns `[P_xx, P_yy, P_zz]`.
pub fn pressure_tensor_diagonal(p_tensor: &[[f64; 3]; 3]) -> [f64; 3] {
    [p_tensor[0][0], p_tensor[1][1], p_tensor[2][2]]
}
/// Compute instantaneous pressure along each axis from the kinetic energy
/// and virial tensor.
///
/// P_alpha = (2 * KE_alpha + W_alpha_alpha) / V
///
/// where KE_alpha = Σ_i m_i v_{i,alpha}^2 / 2.
pub fn component_pressures(
    velocities: &[[f64; 3]],
    masses: &[f64],
    virial_diag: [f64; 3],
    volume: f64,
) -> [f64; 3] {
    let n = velocities.len();
    let mut ke = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            ke[k] += 0.5 * masses[i] * velocities[i][k] * velocities[i][k];
        }
    }
    let inv_v = if volume > 1e-30 { 1.0 / volume } else { 0.0 };
    [
        (2.0 * ke[0] + virial_diag[0]) * inv_v,
        (2.0 * ke[1] + virial_diag[1]) * inv_v,
        (2.0 * ke[2] + virial_diag[2]) * inv_v,
    ]
}
/// Description string for a pressure coupling mode.
#[allow(dead_code)]
pub fn pressure_coupling_mode_name(mode: PressureCouplingMode) -> &'static str {
    match mode {
        PressureCouplingMode::None => "none",
        PressureCouplingMode::Isotropic => "isotropic",
        PressureCouplingMode::SemiIsotropic => "semi-isotropic",
        PressureCouplingMode::Anisotropic => "anisotropic",
    }
}
#[cfg(test)]
mod extended_barostat_tests {
    use super::super::types::*;
    use super::*;
    use oxiphysics_core::Vec3;
    #[test]
    fn test_semi_iso_at_target_pressure() {
        let baro = SemiIsotropicBerendsenBarostat::new(1.0, 1.0, 1.0, 0.01, 0.01);
        let mu_xy = baro.scale_xy(1.0, 0.01);
        let mu_z = baro.scale_z(1.0, 0.01);
        assert!((mu_xy - 1.0).abs() < 1e-12, "mu_xy at target: {mu_xy}");
        assert!((mu_z - 1.0).abs() < 1e-12, "mu_z at target: {mu_z}");
    }
    #[test]
    fn test_semi_iso_over_pressure_xy_shrinks() {
        let baro = SemiIsotropicBerendsenBarostat::new(2.0, 1.0, 1.0, 0.01, 0.01);
        let mu_xy = baro.scale_xy(1.0, 0.01);
        assert!(
            mu_xy < 1.0,
            "lower P_xy than target should shrink box: {mu_xy}"
        );
    }
    #[test]
    fn test_semi_iso_over_pressure_z_shrinks() {
        let baro = SemiIsotropicBerendsenBarostat::new(1.0, 2.0, 1.0, 0.01, 0.01);
        let mu_z = baro.scale_z(1.0, 0.01);
        assert!(mu_z < 1.0, "lower P_z than target should shrink z: {mu_z}");
    }
    #[test]
    fn test_semi_iso_apply_scales_z_only() {
        let baro = SemiIsotropicBerendsenBarostat::new(1.0, 2.0, 1.0, 0.0, 0.01);
        let mut positions = [[1.0, 1.0, 1.0]];
        let mut box_lengths = [10.0, 10.0, 10.0];
        baro.apply_scaling(&mut positions, &mut box_lengths, 1.0, 1.0, 0.01);
        assert!((box_lengths[0] - 10.0).abs() < 1e-12);
        assert!((box_lengths[1] - 10.0).abs() < 1e-12);
        assert!(
            box_lengths[2] != 10.0,
            "z should have changed: {}",
            box_lengths[2]
        );
    }
    #[test]
    fn test_semi_iso_apply_xy_independent_of_z() {
        let baro = SemiIsotropicBerendsenBarostat::new(2.0, 0.5, 1.0, 0.01, 0.01);
        let mut positions = [[1.0, 1.0, 1.0]];
        let mut box_lengths = [10.0, 10.0, 10.0];
        baro.apply_scaling(&mut positions, &mut box_lengths, 1.0, 1.0, 0.01);
        assert!(
            box_lengths[2] > 10.0,
            "z should expand when P_z > target: {}",
            box_lengths[2]
        );
    }
    #[test]
    fn test_box_tensor_cubic_volume() {
        let bt = BoxTensor::cubic(5.0);
        let vol = bt.volume();
        assert!((vol - 125.0).abs() < 1e-10, "cubic 5^3 volume: {vol}");
    }
    #[test]
    fn test_box_tensor_cell_lengths_cubic() {
        let bt = BoxTensor::cubic(3.0);
        let lengths = bt.cell_lengths();
        for l in &lengths {
            assert!((l - 3.0).abs() < 1e-10, "length: {l}");
        }
    }
    #[test]
    fn test_box_tensor_isotropic_pressure() {
        let p = [[2.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 6.0]];
        let p_iso = BoxTensor::isotropic_pressure(&p);
        assert!((p_iso - 4.0).abs() < 1e-12, "isotropic pressure: {p_iso}");
    }
    #[test]
    fn test_box_tensor_advance_at_target() {
        let mut bt = BoxTensor::cubic(5.0);
        bt.p_ref = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p_current = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let h_before = bt.h;
        bt.advance(&p_current, 0.001);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (bt.h[i][j] - h_before[i][j]).abs() < 1e-12,
                    "H[{i}][{j}] should not change at target pressure"
                );
            }
        }
    }
    #[test]
    fn test_box_tensor_advance_off_target() {
        let mut bt = BoxTensor::cubic(5.0);
        bt.mass = 1.0;
        bt.p_ref = [[0.0; 3]; 3];
        let p_current = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        bt.advance(&p_current, 0.001);
        assert!(bt.hdot[0][0].abs() > 0.0, "hdot should be nonzero");
    }
    #[test]
    fn test_virial_tensor_from_pairs_single() {
        let pairs_r = [[1.0, 0.0, 0.0]];
        let pairs_f = [[2.0, 3.0, 0.0]];
        let w = virial_tensor_from_pairs(&pairs_r, &pairs_f);
        assert!((w[0][0] - 2.0).abs() < 1e-12, "W_xx: {}", w[0][0]);
        assert!((w[0][1] - 3.0).abs() < 1e-12, "W_xy: {}", w[0][1]);
        assert!(w[2][2].abs() < 1e-12, "W_zz should be 0: {}", w[2][2]);
    }
    #[test]
    fn test_virial_tensor_from_pairs_two() {
        let pairs_r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let pairs_f = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let w = virial_tensor_from_pairs(&pairs_r, &pairs_f);
        assert!((w[0][0] - 1.0).abs() < 1e-12);
        assert!((w[1][1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_component_pressures_x_only() {
        let velocities = [[2.0, 0.0, 0.0]];
        let masses = [1.0];
        let virial_diag = [0.0, 0.0, 0.0];
        let p = component_pressures(&velocities, &masses, virial_diag, 1.0);
        assert!((p[0] - 4.0).abs() < 1e-12, "P_x: {}", p[0]);
        assert!(p[1].abs() < 1e-12, "P_y should be 0: {}", p[1]);
        assert!(p[2].abs() < 1e-12, "P_z should be 0: {}", p[2]);
    }
    #[test]
    fn test_component_pressures_virial_contribution() {
        let velocities = [[0.0, 0.0, 0.0]];
        let masses = [1.0];
        let virial_diag = [3.0, 6.0, 9.0];
        let p = component_pressures(&velocities, &masses, virial_diag, 1.0);
        assert!((p[0] - 3.0).abs() < 1e-12);
        assert!((p[1] - 6.0).abs() < 1e-12);
        assert!((p[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_mc_rng_creation() {
        let mc = McBarostatRng::new(1.0);
        assert_eq!(mc.n_attempted, 0);
        assert_eq!(mc.n_accepted, 0);
    }
    #[test]
    fn test_mc_rng_acceptance_rate_zero_attempts() {
        let mc = McBarostatRng::new(1.0);
        assert!((mc.acceptance_rate() - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_mc_rng_many_attempts() {
        let mut mc = McBarostatRng::with_seed(1.0, 0.01, 42);
        let v0 = 1000.0;
        for _ in 0..200 {
            mc.try_volume_move(v0, 0.0, 100, 300.0, 1.0);
        }
        assert_eq!(mc.n_attempted, 200);
        assert!(mc.n_accepted > 0, "some moves should be accepted");
        let rate = mc.acceptance_rate();
        assert!(rate > 0.0 && rate <= 1.0, "rate must be in (0,1]: {rate}");
    }
    #[test]
    fn test_mc_rng_reset_statistics() {
        let mut mc = McBarostatRng::with_seed(1.0, 0.01, 99);
        for _ in 0..50 {
            mc.try_volume_move(100.0, 0.0, 10, 300.0, 1.0);
        }
        mc.reset_statistics();
        assert_eq!(mc.n_attempted, 0);
        assert_eq!(mc.n_accepted, 0);
    }
    #[test]
    fn test_mc_rng_negative_delta_energy_accepted() {
        let mut mc = McBarostatRng::with_seed(1.0, 0.001, 111);
        let result = mc.try_volume_move(1000.0, -1e10, 10, 300.0, 1.0);
        assert!(result.is_some(), "hugely favourable move must be accepted");
    }
    #[test]
    fn test_mc_rng_large_positive_energy_rejected() {
        let mut mc = McBarostatRng::with_seed(1.0, 0.001, 222);
        let mut n_accepted = 0u32;
        for _ in 0..20 {
            if mc.try_volume_move(1000.0, 1e30, 10, 300.0, 1.0).is_some() {
                n_accepted += 1;
            }
        }
        assert_eq!(
            n_accepted, 0,
            "hugely unfavourable moves must always be rejected"
        );
    }
    #[test]
    fn test_pressure_coupling_mode_names() {
        assert_eq!(
            pressure_coupling_mode_name(PressureCouplingMode::None),
            "none"
        );
        assert_eq!(
            pressure_coupling_mode_name(PressureCouplingMode::Isotropic),
            "isotropic"
        );
        assert_eq!(
            pressure_coupling_mode_name(PressureCouplingMode::SemiIsotropic),
            "semi-isotropic"
        );
        assert_eq!(
            pressure_coupling_mode_name(PressureCouplingMode::Anisotropic),
            "anisotropic"
        );
    }
    #[test]
    fn test_pressure_coupling_mode_equality() {
        assert_eq!(
            PressureCouplingMode::Isotropic,
            PressureCouplingMode::Isotropic
        );
        assert_ne!(
            PressureCouplingMode::Isotropic,
            PressureCouplingMode::Anisotropic
        );
    }
    #[test]
    fn test_pressure_tensor_diagonal_fn() {
        let p = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let d = pressure_tensor_diagonal(&p);
        assert!((d[0] - 1.0).abs() < 1e-12);
        assert!((d[1] - 5.0).abs() < 1e-12);
        assert!((d[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_no_barostat_is_noop_with_atoms() {
        let mut atoms = AtomSet::new();
        atoms.add_atom(Vec3::new(5.0, 5.0, 5.0), Vec3::zeros(), 1.0, 0.0, 0);
        let mut pbox = PeriodicBox::cubic(10.0);
        let dim_before = pbox.dims;
        let mut baro = NoBarostat::new();
        Barostat::apply(&mut baro, &mut atoms, &mut pbox, 2.0, 0.01);
        assert_eq!(pbox.dims, dim_before, "NoBarostat must not change box");
    }
    /// At constant volume (V_new = V_old) the Jacobian vanishes and the result
    /// is just delta_pot + P₀·ΔV = delta_pot.
    #[test]
    fn test_mc_volume_change_energy_no_volume_change() {
        let baro = MonteCarloBarostat::new(1.0);
        let v = 1000.0;
        let du = baro.compute_volume_change_energy(5.0, v, v, 100, 300.0, 1.0);
        assert!((du - 5.0).abs() < 1e-12, "no volume change: du={du}");
    }
    /// Energy must be finite for any reasonable inputs.
    #[test]
    fn test_mc_volume_change_energy_is_finite() {
        let baro = MonteCarloBarostat::new(1.5);
        let du = baro.compute_volume_change_energy(-10.0, 800.0, 820.0, 50, 300.0, 1.380649e-23);
        assert!(
            du.is_finite(),
            "volume change energy must be finite, got {du}"
        );
    }
    /// Expansion (V_new > V_old) should increase the energy via the P₀·ΔV term
    /// when delta_pot = 0.
    #[test]
    fn test_mc_volume_change_energy_expansion_positive() {
        let baro = MonteCarloBarostat::new(1.0);
        let v_old = 1000.0;
        let v_new = 1010.0;
        let du = baro.compute_volume_change_energy(0.0, v_old, v_new, 1, 1e10, 1.0);
        assert!(du < 0.0, "large T Jacobian should dominate: du={du}");
    }
    /// Isotropic box: all three elements of the scaling matrix must be equal.
    #[test]
    fn test_berendsen_scaling_matrix_isotropic() {
        let baro = BerendsenBarostat::new(1.0, 1.0, 0.01);
        let mu_vec = baro.compute_scaling_matrix(0.8, 0.01);
        let mu = baro.scale_factor(0.8, 0.01);
        for k in 0..3 {
            assert!(
                (mu_vec[k] - mu).abs() < 1e-14,
                "scaling_matrix[{k}] = {}, expected {mu}",
                mu_vec[k]
            );
        }
    }
    /// At target pressure, scaling matrix must be the identity \[1, 1, 1\].
    #[test]
    fn test_berendsen_scaling_matrix_at_target_pressure() {
        let baro = BerendsenBarostat::new(2.0, 1.0, 0.05);
        let mu_vec = baro.compute_scaling_matrix(2.0, 0.01);
        for k in 0..3 {
            assert!(
                (mu_vec[k] - 1.0).abs() < 1e-12,
                "at target P, scaling_matrix[{k}] must be 1; got {}",
                mu_vec[k]
            );
        }
    }
    /// Scaling matrix elements must be positive for reasonable inputs.
    #[test]
    fn test_berendsen_scaling_matrix_positive() {
        let baro = BerendsenBarostat::new(1.0, 0.5, 0.01);
        let mu_vec = baro.compute_scaling_matrix(10.0, 0.001);
        for k in 0..3 {
            assert!(
                mu_vec[k] > 0.0,
                "scaling_matrix[{k}] must be positive; got {}",
                mu_vec[k]
            );
        }
    }
    /// At target pressure the cell derivative must be zero.
    #[test]
    fn test_pr_cell_derivative_at_target() {
        let pr = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
        let dh_dl = pr.compute_cell_derivative(1.0, 1000.0);
        assert!(
            dh_dl.abs() < 1e-12,
            "dH/dL must be 0 at target pressure; got {dh_dl}"
        );
    }
    /// Over-pressure (P > P₀) should yield negative cell derivative
    /// (box wants to expand to reduce pressure).
    #[test]
    fn test_pr_cell_derivative_sign_over_pressure() {
        let pr = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
        let dh_dl = pr.compute_cell_derivative(2.0, 1000.0);
        assert!(
            dh_dl < 0.0,
            "over-pressure: dH/dL should be negative; got {dh_dl}"
        );
    }
    /// Under-pressure (P < P₀) should yield positive cell derivative
    /// (box wants to shrink to raise pressure).
    #[test]
    fn test_pr_cell_derivative_sign_under_pressure() {
        let pr = ParrinelloRahmanBarostat::new(2.0, 1.0, 100.0);
        let dh_dl = pr.compute_cell_derivative(1.0, 1000.0);
        assert!(
            dh_dl > 0.0,
            "under-pressure: dH/dL should be positive; got {dh_dl}"
        );
    }
    /// Cell derivative must scale linearly with volume.
    #[test]
    fn test_pr_cell_derivative_scales_with_volume() {
        let pr = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
        let d1 = pr.compute_cell_derivative(2.0, 1000.0);
        let d2 = pr.compute_cell_derivative(2.0, 2000.0);
        assert!(
            (d2 / d1 - 2.0).abs() < 1e-12,
            "dH/dL should scale linearly with V; ratio={}",
            d2 / d1
        );
    }
}
