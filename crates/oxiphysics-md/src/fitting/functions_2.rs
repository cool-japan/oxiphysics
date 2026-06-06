//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Orientational order parameter P₂ for a set of bond vectors.
///
/// Same as nematic_order_parameter.
pub fn p2_order_parameter(vectors: &[[f64; 3]], director: [f64; 3]) -> f64 {
    if vectors.is_empty() {
        return 0.0;
    }
    let d_norm =
        (director[0] * director[0] + director[1] * director[1] + director[2] * director[2]).sqrt();
    if d_norm < 1e-30 {
        return 0.0;
    }
    let d = [
        director[0] / d_norm,
        director[1] / d_norm,
        director[2] / d_norm,
    ];
    let cos_vals: Vec<f64> = vectors
        .iter()
        .map(|v| v[0] * d[0] + v[1] * d[1] + v[2] * d[2])
        .collect();
    s2_order_parameter(&cos_vals)
}
/// TM-score-like structural alignment score (simplified).
///
/// TM-score ∈ (0, 1]: 1 = perfect match.
/// d₀ = 1.24 * (N - 15)^(1/3) - 1.8 (empirical).
///
/// Score = 1/N * Σ_i 1/(1 + (dᵢ/d₀)²)
pub fn tm_score_simplified(
    coords_a: &[[f64; 3]],
    coords_b: &[[f64; 3]],
    chain_length: usize,
) -> f64 {
    let n = coords_a.len().min(coords_b.len());
    if n == 0 || chain_length == 0 {
        return 0.0;
    }
    let d0 = if chain_length > 15 {
        1.24 * ((chain_length as f64 - 15.0).cbrt()) - 1.8
    } else {
        0.5
    };
    let d0 = d0.max(0.5);
    let score: f64 = (0..n)
        .map(|i| {
            let dx = coords_a[i][0] - coords_b[i][0];
            let dy = coords_a[i][1] - coords_b[i][1];
            let dz = coords_a[i][2] - coords_b[i][2];
            let d2 = dx * dx + dy * dy + dz * dz;
            1.0 / (1.0 + d2 / (d0 * d0))
        })
        .sum::<f64>();
    score / n as f64
}
/// GDT_TS (Global Distance Test, Total Score) — simplified version.
///
/// Fraction of Cα atoms within cutoff distances 1, 2, 4, 8 Å.
pub fn gdt_ts(coords_a: &[[f64; 3]], coords_b: &[[f64; 3]]) -> f64 {
    let n = coords_a.len().min(coords_b.len());
    if n == 0 {
        return 0.0;
    }
    let cutoffs = [1.0, 2.0, 4.0, 8.0_f64];
    let scores: Vec<f64> = cutoffs
        .iter()
        .map(|&cut| {
            let within = (0..n)
                .filter(|&i| {
                    let dx = coords_a[i][0] - coords_b[i][0];
                    let dy = coords_a[i][1] - coords_b[i][1];
                    let dz = coords_a[i][2] - coords_b[i][2];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    d <= cut
                })
                .count();
            within as f64 / n as f64
        })
        .collect();
    scores.iter().sum::<f64>() / 4.0
}
/// MaxSub score: fraction of atoms within 3.5 Å.
pub fn maxsub_score(coords_a: &[[f64; 3]], coords_b: &[[f64; 3]]) -> f64 {
    let n = coords_a.len().min(coords_b.len());
    if n == 0 {
        return 0.0;
    }
    let within = (0..n)
        .filter(|&i| {
            let dx = coords_a[i][0] - coords_b[i][0];
            let dy = coords_a[i][1] - coords_b[i][1];
            let dz = coords_a[i][2] - coords_b[i][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            d <= 3.5
        })
        .count();
    within as f64 / n as f64
}
/// Compute RMSD of each frame in a trajectory relative to the first frame.
pub fn trajectory_rmsd(trajectory: &[Vec<[f64; 3]>]) -> Vec<f64> {
    if trajectory.is_empty() {
        return vec![];
    }
    let reference = &trajectory[0];
    trajectory
        .iter()
        .map(|frame| atom_rmsd(reference.as_slice(), frame.as_slice()))
        .collect()
}
/// Running average of RMSD trajectory (smoothing).
pub fn rmsd_running_average(rmsds: &[f64], window: usize) -> Vec<f64> {
    if rmsds.is_empty() || window == 0 {
        return vec![];
    }
    let w = window.min(rmsds.len());
    rmsds
        .windows(w)
        .map(|win| win.iter().sum::<f64>() / w as f64)
        .collect()
}
/// Plateau RMSD: mean of last N frames.
pub fn plateau_rmsd(rmsds: &[f64], n_tail: usize) -> f64 {
    if rmsds.is_empty() || n_tail == 0 {
        return 0.0;
    }
    let start = rmsds.len().saturating_sub(n_tail);
    let tail = &rmsds[start..];
    if tail.is_empty() {
        return 0.0;
    }
    tail.iter().sum::<f64>() / tail.len() as f64
}
/// Compute pair-wise energy covariance matrix.
///
/// Given energy values for N components across M frames,
/// returns N×N covariance matrix.
pub fn energy_covariance_matrix(energy_matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = energy_matrix.len();
    if m == 0 {
        return vec![];
    }
    let n = energy_matrix[0].len();
    if n == 0 {
        return vec![];
    }
    let means: Vec<f64> = (0..n)
        .map(|j| {
            energy_matrix
                .iter()
                .map(|row| if j < row.len() { row[j] } else { 0.0 })
                .sum::<f64>()
                / m as f64
        })
        .collect();
    let mut cov = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let c = energy_matrix
                .iter()
                .map(|row| {
                    let ri = if i < row.len() { row[i] } else { 0.0 };
                    let rj = if j < row.len() { row[j] } else { 0.0 };
                    (ri - means[i]) * (rj - means[j])
                })
                .sum::<f64>();
            cov[i][j] = c / m as f64;
        }
    }
    cov
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use std::f64::consts::PI;
    fn make_trajectory(n_frames: usize, n_atoms: usize, oscillation: f64) -> Vec<Vec<[f64; 3]>> {
        (0..n_frames)
            .map(|t| {
                (0..n_atoms)
                    .map(|i| [i as f64 + oscillation * (t as f64 * 0.1).sin(), 0.0, 0.0])
                    .collect()
            })
            .collect()
    }
    #[test]
    fn test_atom_rmsd_identical() {
        let coords: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        assert!(atom_rmsd(&coords, &coords).abs() < 1e-12);
    }
    #[test]
    fn test_atom_rmsd_shifted() {
        let a: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let b: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let r = atom_rmsd(&a, &b);
        assert!((r - 1.0).abs() < 1e-10, "RMSD of uniform shift = 1: {r}");
    }
    #[test]
    fn test_atom_rmsd_empty() {
        assert_eq!(atom_rmsd(&[], &[]), 0.0);
    }
    #[test]
    fn test_atom_rmsd_non_negative() {
        let a: Vec<[f64; 3]> = vec![[0.0, 1.0, 2.0]];
        let b: Vec<[f64; 3]> = vec![[3.0, 4.0, 5.0]];
        assert!(atom_rmsd(&a, &b) >= 0.0);
    }
    #[test]
    fn test_rmsf_static_trajectory_zero() {
        let traj: Vec<Vec<[f64; 3]>> = vec![vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]; 5];
        let rmsf = rmsf_per_atom(&traj);
        for &v in &rmsf {
            assert!(v.abs() < 1e-10, "Static trajectory → zero RMSF: {v}");
        }
    }
    #[test]
    fn test_rmsf_length_equals_n_atoms() {
        let traj = make_trajectory(10, 4, 0.5);
        let rmsf = rmsf_per_atom(&traj);
        assert_eq!(rmsf.len(), 4, "RMSF length should equal n_atoms");
    }
    #[test]
    fn test_rmsf_positive_for_oscillation() {
        let traj = make_trajectory(20, 3, 1.0);
        let rmsf = rmsf_per_atom(&traj);
        for &v in &rmsf {
            assert!(v >= 0.0, "RMSF must be non-negative: {v}");
        }
        let total: f64 = rmsf.iter().sum();
        assert!(total > 0.0, "Some atoms should have positive RMSF");
    }
    #[test]
    fn test_mean_rmsf_empty() {
        assert_eq!(mean_rmsf(&[]), 0.0);
    }
    #[test]
    fn test_mean_rmsf_positive() {
        let traj = make_trajectory(10, 3, 0.5);
        let m = mean_rmsf(&traj);
        assert!(m >= 0.0, "Mean RMSF must be non-negative: {m}");
    }
    #[test]
    fn test_autocorr_zero_lag_equals_one() {
        let traj = make_trajectory(20, 2, 1.0);
        let acf = displacement_autocorrelation(&traj, 0, 5);
        assert!(
            (acf[0] - 1.0).abs() < 1e-10,
            "ACF at lag 0 must be 1: {}",
            acf[0]
        );
    }
    #[test]
    fn test_autocorr_length_equals_max_lag() {
        let traj = make_trajectory(15, 3, 0.5);
        let acf = displacement_autocorrelation(&traj, 0, 5);
        assert_eq!(acf.len(), 5, "ACF length should equal max_lag");
    }
    #[test]
    fn test_autocorr_empty_trajectory() {
        let acf = displacement_autocorrelation(&[], 0, 5);
        assert!(acf.is_empty());
    }
    #[test]
    fn test_autocorr_static_trajectory() {
        let traj: Vec<Vec<[f64; 3]>> = vec![vec![[1.0, 0.0, 0.0]]; 5];
        let acf = displacement_autocorrelation(&traj, 0, 3);
        for &v in &acf {
            assert!(
                v.abs() < 1e-10,
                "Static ACF should be zero (normalized): {v}"
            );
        }
    }
    #[test]
    fn test_s2_all_aligned() {
        let cos_vals = vec![1.0, 1.0, 1.0, 1.0];
        let s2 = s2_order_parameter(&cos_vals);
        assert!((s2 - 1.0).abs() < 1e-10, "All aligned: S₂ = 1, got {s2}");
    }
    #[test]
    fn test_s2_all_perpendicular() {
        let cos_vals = vec![0.0, 0.0, 0.0, 0.0];
        let s2 = s2_order_parameter(&cos_vals);
        assert!(
            (s2 - (-0.5)).abs() < 1e-10,
            "All perpendicular: S₂ = -0.5, got {s2}"
        );
    }
    #[test]
    fn test_s2_random_range() {
        let cos_vals: Vec<f64> = (0..100).map(|i| i as f64 * 0.02 - 1.0).collect();
        let s2 = s2_order_parameter(&cos_vals);
        assert!(
            (-0.5 - 1e-10..=1.0 + 1e-10).contains(&s2),
            "S₂ out of range: {s2}"
        );
    }
    #[test]
    fn test_s2_empty() {
        assert_eq!(s2_order_parameter(&[]), 0.0);
    }
    #[test]
    fn test_p2_all_along_director() {
        let director = [0.0, 0.0, 1.0];
        let vecs: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0]; 5];
        let p2 = p2_order_parameter(&vecs, director);
        assert!(
            (p2 - 1.0).abs() < 1e-10,
            "All along director: P₂ = 1, got {p2}"
        );
    }
    #[test]
    fn test_p2_range() {
        let director = [1.0, 0.0, 0.0];
        let vecs: Vec<[f64; 3]> = (0..10)
            .map(|i| {
                let a = i as f64 * PI / 10.0;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let p2 = p2_order_parameter(&vecs, director);
        assert!(
            (-0.5 - 1e-10..=1.0 + 1e-10).contains(&p2),
            "P₂ out of range: {p2}"
        );
    }
    #[test]
    fn test_tm_score_identical_structures() {
        let coords: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let score = tm_score_simplified(&coords, &coords, 10);
        assert!(
            (score - 1.0).abs() < 1e-10,
            "Identical → TM-score = 1: {score}"
        );
    }
    #[test]
    fn test_tm_score_range() {
        let a: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let b: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 + 5.0, 0.0, 0.0]).collect();
        let score = tm_score_simplified(&a, &b, 10);
        assert!(
            (0.0..=1.0).contains(&score),
            "TM-score out of [0,1]: {score}"
        );
    }
    #[test]
    fn test_tm_score_empty() {
        let score = tm_score_simplified(&[], &[], 10);
        assert_eq!(score, 0.0);
    }
    #[test]
    fn test_tm_score_decreases_with_distance() {
        let a: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]];
        let b_near: Vec<[f64; 3]> = vec![[0.1, 0.0, 0.0]];
        let b_far: Vec<[f64; 3]> = vec![[10.0, 0.0, 0.0]];
        let s_near = tm_score_simplified(&a, &b_near, 20);
        let s_far = tm_score_simplified(&a, &b_far, 20);
        assert!(
            s_near > s_far,
            "Nearer structure should have higher TM-score: {s_near} vs {s_far}"
        );
    }
    #[test]
    fn test_gdt_ts_identical() {
        let coords: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let score = gdt_ts(&coords, &coords);
        assert!(
            (score - 1.0).abs() < 1e-10,
            "Identical → GDT_TS = 1: {score}"
        );
    }
    #[test]
    fn test_gdt_ts_empty() {
        assert_eq!(gdt_ts(&[], &[]), 0.0);
    }
    #[test]
    fn test_gdt_ts_range() {
        let a: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let b: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 + 2.5, 0.0, 0.0]).collect();
        let score = gdt_ts(&a, &b);
        assert!((0.0..=1.0).contains(&score), "GDT_TS out of range: {score}");
    }
    #[test]
    fn test_maxsub_identical() {
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let score = maxsub_score(&coords, &coords);
        assert!(
            (score - 1.0).abs() < 1e-10,
            "Identical → maxsub = 1: {score}"
        );
    }
    #[test]
    fn test_maxsub_empty() {
        assert_eq!(maxsub_score(&[], &[]), 0.0);
    }
    #[test]
    fn test_maxsub_large_deviation() {
        let a: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]];
        let b: Vec<[f64; 3]> = vec![[10.0, 0.0, 0.0]];
        let score = maxsub_score(&a, &b);
        assert_eq!(score, 0.0, "Large deviation → 0 maxsub score");
    }
    #[test]
    fn test_trajectory_rmsd_first_frame_zero() {
        let traj = make_trajectory(5, 3, 0.5);
        let rmsds = trajectory_rmsd(&traj);
        assert!(
            rmsds[0].abs() < 1e-10,
            "RMSD of first frame to itself = 0: {}",
            rmsds[0]
        );
    }
    #[test]
    fn test_trajectory_rmsd_length() {
        let traj = make_trajectory(8, 3, 0.5);
        let rmsds = trajectory_rmsd(&traj);
        assert_eq!(
            rmsds.len(),
            8,
            "RMSD trajectory length should equal n_frames"
        );
    }
    #[test]
    fn test_trajectory_rmsd_non_negative() {
        let traj = make_trajectory(10, 4, 1.0);
        let rmsds = trajectory_rmsd(&traj);
        for &r in &rmsds {
            assert!(r >= 0.0, "RMSD must be non-negative: {r}");
        }
    }
    #[test]
    fn test_running_average_constant_series() {
        let rmsds = vec![2.0; 10];
        let avg = rmsd_running_average(&rmsds, 3);
        for &v in &avg {
            assert!((v - 2.0).abs() < 1e-10, "Running avg of const = const: {v}");
        }
    }
    #[test]
    fn test_running_average_length() {
        let rmsds = vec![1.0; 10];
        let avg = rmsd_running_average(&rmsds, 3);
        assert_eq!(
            avg.len(),
            10 - 3 + 1,
            "Running average length: {}",
            avg.len()
        );
    }
    #[test]
    fn test_plateau_rmsd() {
        let rmsds: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let p = plateau_rmsd(&rmsds, 5);
        assert!((p - 17.0).abs() < 1e-10, "Plateau RMSD = 17: {p}");
    }
    #[test]
    fn test_plateau_rmsd_empty() {
        assert_eq!(plateau_rmsd(&[], 5), 0.0);
    }
    #[test]
    fn test_energy_cov_diagonal_non_negative() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, (10 - i) as f64]).collect();
        let cov = energy_covariance_matrix(&data);
        assert_eq!(cov.len(), 2);
        assert!(
            cov[0][0] >= 0.0,
            "Variance must be non-negative: {}",
            cov[0][0]
        );
        assert!(
            cov[1][1] >= 0.0,
            "Variance must be non-negative: {}",
            cov[1][1]
        );
    }
    #[test]
    fn test_energy_cov_symmetry() {
        let data: Vec<Vec<f64>> = (0..10)
            .map(|i| vec![i as f64, (i as f64).sin(), (i as f64).cos()])
            .collect();
        let cov = energy_covariance_matrix(&data);
        for (i, row) in cov.iter().enumerate() {
            for (j, &cij) in row.iter().enumerate() {
                assert!(
                    (cij - cov[j][i]).abs() < 1e-10,
                    "Covariance matrix must be symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_energy_cov_constant_series_zero_variance() {
        let data: Vec<Vec<f64>> = vec![vec![1.0, 2.0]; 5];
        let cov = energy_covariance_matrix(&data);
        assert!(
            cov[0][0].abs() < 1e-10,
            "Constant series → zero variance: {}",
            cov[0][0]
        );
        assert!(
            cov[1][1].abs() < 1e-10,
            "Constant series → zero variance: {}",
            cov[1][1]
        );
    }
    #[test]
    fn test_energy_cov_empty() {
        let cov = energy_covariance_matrix(&[]);
        assert!(cov.is_empty());
    }
    #[test]
    fn test_nematic_order_empty() {
        assert_eq!(nematic_order_parameter(&[]), 0.0);
    }
    #[test]
    fn test_nematic_order_all_parallel() {
        let vecs: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; 5];
        let q = nematic_order_parameter(&vecs);
        assert!((q - 1.0).abs() < 1e-10, "All parallel → Q = 1: {q}");
    }
    #[test]
    fn test_nematic_order_range() {
        let vecs: Vec<[f64; 3]> = (0..8)
            .map(|i| {
                let a = i as f64 * PI / 4.0;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let q = nematic_order_parameter(&vecs);
        assert!(
            (-0.5 - 1e-10..=1.0 + 1e-10).contains(&q),
            "Nematic Q out of range: {q}"
        );
    }
    #[test]
    fn test_atom_rmsd_3d_displacement() {
        let a = vec![[0.0, 0.0, 0.0_f64]];
        let b = vec![[3.0, 4.0, 0.0_f64]];
        let r = atom_rmsd(&a, &b);
        assert!(
            (r - 5.0).abs() < 1e-10,
            "RMSD = 5 for 3-4-5 displacement: {r}"
        );
    }
    #[test]
    fn test_rmsf_oscillating_atom_positive() {
        let traj: Vec<Vec<[f64; 3]>> = (0..20)
            .map(|t| vec![[(t as f64 * 0.1).sin(), 0.0, 0.0], [1.0, 0.0, 0.0]])
            .collect();
        let rmsf = rmsf_per_atom(&traj);
        assert!(
            rmsf[0] > rmsf[1],
            "Oscillating atom should have larger RMSF: {} vs {}",
            rmsf[0],
            rmsf[1]
        );
    }
}
