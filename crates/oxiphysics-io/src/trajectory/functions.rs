//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TrajectoryFrame;

/// Extract timestep and time from an XYZ comment line.
/// Expected format: `Timestep N time=T`  (falls back to 0 / 0.0 if not matched).
pub(super) fn parse_xyz_comment(comment: &str) -> (u64, f64) {
    let mut timestep = 0u64;
    let mut time = 0.0f64;
    for token in comment.split_whitespace() {
        if let Some(rest) = token.strip_prefix("time=") {
            if let Ok(t) = rest.parse::<f64>() {
                time = t;
            }
        } else if let Ok(ts) = token.parse::<u64>() {
            timestep = ts;
        }
    }
    (timestep, time)
}
/// Count the byte length of the first `n` lines in `s` (including newlines).
pub(super) fn count_line_bytes(s: &str, n: usize) -> usize {
    let mut count = 0;
    let mut bytes = 0;
    for ch in s.chars() {
        if count >= n {
            break;
        }
        bytes += ch.len_utf8();
        if ch == '\n' {
            count += 1;
        }
    }
    bytes
}
/// Compute the root-mean-square deviation (RMSD) between two sets of positions.
///
/// Both frames must have the same number of atoms.  Returns 0.0 if the frames
/// are identical.
pub fn compute_rmsd(frame_a: &TrajectoryFrame, frame_b: &TrajectoryFrame) -> f64 {
    assert_eq!(
        frame_a.n_atoms(),
        frame_b.n_atoms(),
        "RMSD requires equal atom counts"
    );
    if frame_a.n_atoms() == 0 {
        return 0.0;
    }
    let sum_sq: f64 = frame_a
        .positions
        .iter()
        .zip(frame_b.positions.iter())
        .map(|(a, b)| {
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            let dz = a[2] - b[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum_sq / frame_a.n_atoms() as f64).sqrt()
}
/// Compute the center of mass for a set of positions and corresponding masses.
///
/// `masses` must have the same length as `positions`.
pub fn center_of_mass(positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
    assert_eq!(
        positions.len(),
        masses.len(),
        "positions and masses must match"
    );
    let total_mass: f64 = masses.iter().sum();
    if total_mass == 0.0 {
        return [0.0; 3];
    }
    let mut com = [0.0f64; 3];
    for (pos, &m) in positions.iter().zip(masses.iter()) {
        com[0] += m * pos[0];
        com[1] += m * pos[1];
        com[2] += m * pos[2];
    }
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;
    com
}
/// Compute the nearest rotation matrix to `h` using iterative polar
/// decomposition: R_{k+1} = (R_k + (R_k^{-T})) / 2.
///
/// Returns a 3×3 rotation matrix as `[[f64; 3\]; 3]`.
pub(super) fn polar_rotation_3x3(h: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = h;
    for _ in 0..50 {
        let inv_t = mat3_inv_transpose(r);
        let mut next = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                next[i][j] = 0.5 * (r[i][j] + inv_t[i][j]);
            }
        }
        let mut diff = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let d = next[i][j] - r[i][j];
                diff += d * d;
            }
        }
        r = next;
        if diff < 1e-20 {
            break;
        }
    }
    r
}
/// Compute the inverse-transpose of a 3×3 matrix.
pub(super) fn mat3_inv_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-30 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let inv_det = 1.0 / det;
    let c00 = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    let c01 = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * inv_det;
    let c02 = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    let c10 = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * inv_det;
    let c11 = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    let c12 = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * inv_det;
    let c20 = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    let c21 = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * inv_det;
    let c22 = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;
    [[c00, c10, c20], [c01, c11, c21], [c02, c12, c22]]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::trajectory::types::*;
    use oxiphysics_core::math::Vec3;
    fn make_positions(n: usize, offset: f64) -> Vec<Vec3> {
        (0..n)
            .map(|i| Vec3::new(i as f64 + offset, 0.0, 0.0))
            .collect()
    }
    #[test]
    fn test_trajectory_frame_count() {
        let mut traj = TrajectoryWriter::new();
        assert_eq!(traj.frame_count(), 0);
        traj.add_frame(0.0, &make_positions(4, 0.0));
        traj.add_frame(0.5, &make_positions(4, 1.0));
        traj.add_frame(1.0, &make_positions(4, 2.0));
        assert_eq!(traj.frame_count(), 3);
    }
    #[test]
    fn test_trajectory_xyz_format() {
        let n = 3;
        let mut traj = TrajectoryWriter::new();
        traj.add_frame(0.0, &make_positions(n, 0.0));
        traj.add_frame(1.0, &make_positions(n, 10.0));
        let mut buf = Vec::new();
        traj.write_xyz(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines.len(),
            2 * (n + 2),
            "unexpected line count in XYZ output"
        );
        assert_eq!(lines[0].trim(), n.to_string());
        assert_eq!(lines[n + 2].trim(), n.to_string());
    }
    #[test]
    fn test_trajectory_xyz_coord_lines() {
        let positions = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0)];
        let mut traj = TrajectoryWriter::new();
        traj.add_frame(0.0, &positions);
        let mut buf = Vec::new();
        traj.write_xyz(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let coord_lines: Vec<&str> = s.lines().filter(|l| l.starts_with('X')).collect();
        assert_eq!(coord_lines.len(), 2);
        assert!(
            coord_lines[0].contains('1')
                && coord_lines[0].contains('2')
                && coord_lines[0].contains('3')
        );
    }
    #[test]
    fn test_trajectory_xdmf_output() {
        let mut traj = TrajectoryWriter::new();
        traj.add_frame(0.0, &make_positions(2, 0.0));
        traj.add_frame(1.0, &make_positions(2, 5.0));
        let mut buf = Vec::new();
        traj.write_xdmf(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("<?xml"), "missing XML declaration");
        assert!(s.contains("Temporal"), "missing Temporal collection");
    }
    #[test]
    fn test_trajectory_default() {
        let traj = TrajectoryWriter::default();
        assert_eq!(traj.frame_count(), 0);
    }
    fn make_frame(timestep: u64, time: f64, n: usize) -> TrajectoryFrame {
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, i as f64 * 0.5, 0.0]).collect();
        let atom_types: Vec<String> = (0..n)
            .map(|i| {
                if i % 2 == 0 {
                    "C".to_string()
                } else {
                    "H".to_string()
                }
            })
            .collect();
        TrajectoryFrame::new(timestep, time, positions, atom_types)
    }
    #[test]
    fn test_xyz_write_frame_format() {
        let frame = make_frame(0, 0.0, 3);
        let s = XyzWriter::write_frame(&frame);
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "3");
        assert!(lines[1].contains("Timestep"));
        assert_eq!(lines.len(), 5);
        assert!(lines[2].starts_with("C ") || lines[2].starts_with("C\t"));
    }
    #[test]
    fn test_xyz_roundtrip_single_frame() {
        let frame = make_frame(42, 1.5, 4);
        let s = XyzWriter::write_frame(&frame);
        let parsed = XyzReader::read_frame(&s).unwrap();
        assert_eq!(parsed.timestep, 42);
        assert!((parsed.time - 1.5).abs() < 1e-9);
        assert_eq!(parsed.n_atoms(), 4);
        for i in 0..4 {
            assert!((parsed.positions[i][0] - frame.positions[i][0]).abs() < 1e-10);
            assert!((parsed.positions[i][1] - frame.positions[i][1]).abs() < 1e-10);
            assert!((parsed.positions[i][2] - frame.positions[i][2]).abs() < 1e-10);
            assert_eq!(parsed.atom_types[i], frame.atom_types[i]);
        }
    }
    #[test]
    fn test_xyz_multi_frame_write_read() {
        let frames = vec![
            make_frame(0, 0.0, 3),
            make_frame(1, 0.5, 3),
            make_frame(2, 1.0, 3),
        ];
        let s = XyzWriter::write_frames(&frames);
        let parsed = XyzReader::read_all_frames(&s).unwrap();
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed[0].timestep, 0);
        assert_eq!(parsed[1].timestep, 1);
        assert_eq!(parsed[2].timestep, 2);
        assert!((parsed[2].time - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_xyz_read_frame_bad_input() {
        assert!(XyzReader::read_frame("").is_err());
        assert!(XyzReader::read_frame("not_a_number\ncomment\n").is_err());
    }
    #[test]
    fn test_rmsd_identical_frames() {
        let frame = make_frame(0, 0.0, 5);
        let rmsd = compute_rmsd(&frame, &frame);
        assert!(rmsd.abs() < 1e-15);
    }
    #[test]
    fn test_rmsd_known_displacement() {
        let frame_a = TrajectoryFrame::new(
            0,
            0.0,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec!["C".to_string(), "C".to_string()],
        );
        let frame_b = TrajectoryFrame::new(
            1,
            0.0,
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec!["C".to_string(), "C".to_string()],
        );
        let rmsd = compute_rmsd(&frame_a, &frame_b);
        assert!((rmsd - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_center_of_mass_equal_masses() {
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let masses = vec![1.0, 1.0, 1.0];
        let com = center_of_mass(&positions, &masses);
        assert!((com[0] - 1.0).abs() < 1e-10);
        assert!((com[1] - 2.0 / 3.0).abs() < 1e-10);
        assert!(com[2].abs() < 1e-10);
    }
    #[test]
    fn test_center_of_mass_unequal_masses() {
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let masses = vec![3.0, 1.0];
        let com = center_of_mass(&positions, &masses);
        assert!((com[0] - 0.75).abs() < 1e-10);
    }
    #[test]
    fn test_lammps_dump_frame_format() {
        let frame = make_frame(10, 0.1, 2);
        let dump = TrajLammpsWriter::write_dump_frame(&frame, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
        assert!(dump.contains("ITEM: TIMESTEP"));
        assert!(dump.contains("10\n"));
        assert!(dump.contains("ITEM: NUMBER OF ATOMS"));
        assert!(dump.contains("2\n"));
        assert!(dump.contains("ITEM: BOX BOUNDS"));
        assert!(dump.contains("ITEM: ATOMS id type x y z"));
        assert!(dump.contains("1 C ") || dump.contains("1 C\t"));
    }
    #[test]
    fn test_resample_empty() {
        let result = TrajectoryResampler::resample_uniform(&[], 1.0);
        assert!(result.is_empty());
    }
    #[test]
    fn test_resample_single_frame() {
        let frames = vec![make_frame(0, 0.0, 3)];
        let result = TrajectoryResampler::resample_uniform(&frames, 1.0);
        assert_eq!(result.len(), 1);
    }
    #[test]
    fn test_resample_uniform() {
        let frames = vec![make_frame(0, 0.0, 3), make_frame(10, 1.0, 3)];
        let result = TrajectoryResampler::resample_uniform(&frames, 0.5);
        assert_eq!(result.len(), 3);
        assert!((result[0].time - 0.0).abs() < 1e-10);
        assert!((result[1].time - 0.5).abs() < 1e-10);
        assert!((result[2].time - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_subsample() {
        let frames = vec![
            make_frame(0, 0.0, 2),
            make_frame(1, 0.1, 2),
            make_frame(2, 0.2, 2),
            make_frame(3, 0.3, 2),
            make_frame(4, 0.4, 2),
        ];
        let result = TrajectoryResampler::subsample(&frames, 2);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].timestep, 0);
        assert_eq!(result[1].timestep, 2);
        assert_eq!(result[2].timestep, 4);
    }
    #[test]
    fn test_subsample_every_one() {
        let frames = vec![make_frame(0, 0.0, 2), make_frame(1, 0.1, 2)];
        let result = TrajectoryResampler::subsample(&frames, 1);
        assert_eq!(result.len(), 2);
    }
    #[test]
    fn test_concatenate_empty() {
        let result = TrajectoryConcatenator::concatenate(&[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_concatenate_two_trajectories() {
        let traj1 = vec![make_frame(0, 0.0, 2), make_frame(1, 1.0, 2)];
        let traj2 = vec![make_frame(0, 0.0, 2), make_frame(1, 1.0, 2)];
        let result = TrajectoryConcatenator::concatenate(&[traj1, traj2]);
        assert_eq!(result.len(), 4);
        assert!((result[2].time - 1.0).abs() < 1e-10);
        assert!((result[3].time - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_merge_sorted() {
        let traj1 = vec![make_frame(0, 0.0, 2), make_frame(2, 2.0, 2)];
        let traj2 = vec![make_frame(1, 1.0, 2), make_frame(3, 3.0, 2)];
        let result = TrajectoryConcatenator::merge_sorted(&[traj1, traj2]);
        assert_eq!(result.len(), 4);
        assert!((result[0].time - 0.0).abs() < 1e-10);
        assert!((result[1].time - 1.0).abs() < 1e-10);
        assert!((result[2].time - 2.0).abs() < 1e-10);
        assert!((result[3].time - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_filter_time_range() {
        let frames = vec![
            make_frame(0, 0.0, 2),
            make_frame(1, 0.5, 2),
            make_frame(2, 1.0, 2),
            make_frame(3, 1.5, 2),
        ];
        let result = TrajectoryFilter::time_range(&frames, 0.5, 1.0);
        assert_eq!(result.len(), 2);
        assert!((result[0].time - 0.5).abs() < 1e-10);
        assert!((result[1].time - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_filter_atom_types() {
        let frame = make_frame(0, 0.0, 4);
        let result = TrajectoryFilter::filter_atom_types(&[frame], &["C"]);
        assert_eq!(result[0].n_atoms(), 2);
        assert!(result[0].atom_types.iter().all(|t| t == "C"));
    }
    #[test]
    fn test_filter_remove_empty() {
        let mut frames = vec![make_frame(0, 0.0, 2), make_frame(1, 0.5, 0)];
        frames[1].positions.clear();
        frames[1].atom_types.clear();
        let result = TrajectoryFilter::remove_empty(&frames);
        assert_eq!(result.len(), 1);
    }
    #[test]
    fn test_wrap_positions() {
        let mut positions = [[1.5, 2.5, -0.5]];
        PeriodicImageHandler::wrap_positions(&mut positions, [1.0, 1.0, 1.0]);
        assert!((positions[0][0] - 0.5).abs() < 1e-10);
        assert!((positions[0][1] - 0.5).abs() < 1e-10);
        assert!((positions[0][2] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_minimum_image_distance() {
        let d = PeriodicImageHandler::minimum_image_distance(
            [0.1, 0.0, 0.0],
            [0.9, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        );
        assert!(
            (d - 0.2).abs() < 1e-10,
            "Minimum image distance should be 0.2, got {d}"
        );
    }
    #[test]
    fn test_minimum_image_same_point() {
        let d = PeriodicImageHandler::minimum_image_distance(
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [1.0, 1.0, 1.0],
        );
        assert!(d.abs() < 1e-14);
    }
    #[test]
    fn test_unwrap_trajectory() {
        let mut frames = vec![
            TrajectoryFrame::new(0, 0.0, vec![[0.1, 0.0, 0.0]], vec!["C".to_string()]),
            TrajectoryFrame::new(1, 1.0, vec![[0.9, 0.0, 0.0]], vec!["C".to_string()]),
        ];
        PeriodicImageHandler::unwrap_trajectory(&mut frames, [1.0, 1.0, 1.0]);
        let dx = frames[1].positions[0][0] - frames[0].positions[0][0];
        assert!(
            dx.abs() < 0.5,
            "Unwrapped displacement should be small, got {dx}"
        );
    }
    #[test]
    fn test_mean_positions() {
        let frames = vec![
            TrajectoryFrame::new(0, 0.0, vec![[0.0, 0.0, 0.0]], vec!["C".to_string()]),
            TrajectoryFrame::new(1, 1.0, vec![[2.0, 4.0, 6.0]], vec!["C".to_string()]),
        ];
        let mean = TrajectoryStatistics::mean_positions(&frames);
        assert_eq!(mean.len(), 1);
        assert!((mean[0][0] - 1.0).abs() < 1e-10);
        assert!((mean[0][1] - 2.0).abs() < 1e-10);
        assert!((mean[0][2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_rmsf_constant_positions() {
        let frame = TrajectoryFrame::new(0, 0.0, vec![[1.0, 2.0, 3.0]], vec!["C".to_string()]);
        let frames = vec![frame.clone(), frame.clone(), frame];
        let rmsf = TrajectoryStatistics::rmsf(&frames);
        assert_eq!(rmsf.len(), 1);
        assert!(
            rmsf[0].abs() < 1e-14,
            "RMSF should be 0 for constant positions"
        );
    }
    #[test]
    fn test_radius_of_gyration() {
        let frame = TrajectoryFrame::new(
            0,
            0.0,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec!["C".to_string(), "C".to_string()],
        );
        let masses = [1.0, 1.0];
        let rg = TrajectoryStatistics::radius_of_gyration(&frame, &masses);
        assert!((rg - 0.5).abs() < 1e-10, "Rg should be 0.5, got {rg}");
    }
    #[test]
    fn test_end_to_end_distance() {
        let frame = TrajectoryFrame::new(
            0,
            0.0,
            vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
            vec!["C".to_string(), "C".to_string()],
        );
        let d = TrajectoryStatistics::end_to_end_distance(&frame);
        assert!(
            (d - 5.0).abs() < 1e-10,
            "End-to-end distance should be 5.0, got {d}"
        );
    }
    #[test]
    fn test_end_to_end_single_atom() {
        let frame = TrajectoryFrame::new(0, 0.0, vec![[1.0, 2.0, 3.0]], vec!["C".to_string()]);
        assert_eq!(TrajectoryStatistics::end_to_end_distance(&frame), 0.0);
    }
    #[test]
    fn test_mean_positions_empty() {
        let mean = TrajectoryStatistics::mean_positions(&[]);
        assert!(mean.is_empty());
    }
    #[test]
    fn test_rmsf_empty() {
        let rmsf = TrajectoryStatistics::rmsf(&[]);
        assert!(rmsf.is_empty());
    }
    fn make_diffusing_frames(n_atoms: usize, n_frames: usize, dt: f64) -> Vec<TrajectoryFrame> {
        (0..n_frames)
            .map(|fi| {
                let t = fi as f64 * dt;
                let positions: Vec<[f64; 3]> = (0..n_atoms).map(|_| [t, 0.0, 0.0]).collect();
                let types: Vec<String> = (0..n_atoms).map(|_| "C".to_string()).collect();
                TrajectoryFrame::new(fi as u64, t, positions, types)
            })
            .collect()
    }
    #[test]
    fn test_msd_zero_at_lag0() {
        let frames = make_diffusing_frames(5, 10, 1.0);
        let msd = MsdCalculator::compute(&frames, 5);
        assert_eq!(msd[0].0, 0);
        assert!(msd[0].1.abs() < 1e-10, "MSD at lag=0 should be 0");
    }
    #[test]
    fn test_msd_ballistic() {
        let frames = make_diffusing_frames(1, 6, 1.0);
        let msd = MsdCalculator::compute(&frames, 3);
        assert!((msd[1].1 - 1.0).abs() < 1e-8, "MSD at lag=1: {}", msd[1].1);
        assert!((msd[2].1 - 4.0).abs() < 1e-8, "MSD at lag=2: {}", msd[2].1);
    }
    #[test]
    fn test_msd_empty_frames() {
        let msd = MsdCalculator::compute(&[], 5);
        assert_eq!(msd.len(), 1);
        assert_eq!(msd[0], (0, 0.0));
    }
    #[test]
    fn test_diffusion_coefficient_zero_slope() {
        let msd: Vec<(usize, f64)> = (0..5).map(|i| (i, 1.0)).collect();
        let d = MsdCalculator::diffusion_coefficient(&msd, 1.0);
        assert!(d.abs() < 1.0, "D should be near 0 for flat MSD, got {d}");
    }
    fn make_velocity_frames(n_atoms: usize, n_frames: usize) -> Vec<VelocityFrame> {
        (0..n_frames)
            .map(|fi| {
                let velocities: Vec<[f64; 3]> = (0..n_atoms).map(|_| [1.0, 0.0, 0.0]).collect();
                VelocityFrame::new(fi as u64, fi as f64, velocities)
            })
            .collect()
    }
    #[test]
    fn test_vacf_constant_velocity() {
        let frames = make_velocity_frames(3, 5);
        let vacf = VacfCalculator::compute_normalized(&frames, 3);
        for (_, cv) in &vacf {
            assert!(
                (cv - 1.0).abs() < 1e-10,
                "VACF should be 1 for constant velocity, got {cv}"
            );
        }
    }
    #[test]
    fn test_vacf_empty() {
        let vacf = VacfCalculator::compute_normalized(&[], 3);
        assert!(vacf.is_empty());
    }
    #[test]
    fn test_vacf_integrate_constant() {
        let vacf: Vec<(usize, f64)> = (0..5).map(|i| (i, 1.0)).collect();
        let d = VacfCalculator::integrate_vacf(&vacf, 0.5);
        assert!((d - 2.0 / 3.0).abs() < 1e-10, "D={d}");
    }
    fn make_two_atom_frame(r: f64) -> TrajectoryFrame {
        TrajectoryFrame::new(
            0,
            0.0,
            vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]],
            vec!["C".to_string(), "C".to_string()],
        )
    }
    #[test]
    fn test_bond_length_known() {
        let frame = make_two_atom_frame(2.5);
        let bl = BondLengthAnalyser::bond_length(&frame, 0, 1);
        assert!((bl - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_bond_length_time_series() {
        let frames = vec![
            make_two_atom_frame(1.5),
            make_two_atom_frame(2.0),
            make_two_atom_frame(2.5),
        ];
        let ts = BondLengthAnalyser::time_series(&frames, 0, 1);
        assert_eq!(ts.len(), 3);
        assert!((ts[0] - 1.5).abs() < 1e-10);
        assert!((ts[2] - 2.5).abs() < 1e-10);
    }
    #[test]
    fn test_bond_length_mean_std() {
        let frames = vec![
            make_two_atom_frame(1.0),
            make_two_atom_frame(2.0),
            make_two_atom_frame(3.0),
        ];
        let (mean, _std) = BondLengthAnalyser::mean_and_std(&frames, 0, 1);
        assert!((mean - 2.0).abs() < 1e-10, "mean={mean}");
    }
    #[test]
    fn test_bond_angle_90() {
        let frame = TrajectoryFrame::new(
            0,
            0.0,
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
        );
        let angle = BondLengthAnalyser::bond_angle_deg(&frame, 0, 1, 2);
        assert!((angle - 90.0).abs() < 1e-8, "angle={angle}");
    }
    #[test]
    fn test_bond_angle_180() {
        let frame = TrajectoryFrame::new(
            0,
            0.0,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
        );
        let angle = BondLengthAnalyser::bond_angle_deg(&frame, 0, 1, 2);
        assert!((angle - 180.0).abs() < 1e-8, "angle={angle}");
    }
    #[test]
    fn test_rdf_empty_frames() {
        let (r, g) = RdfCalculator::compute(&[], 5.0, 10, [10.0; 3]);
        assert!(r.is_empty());
        assert!(g.is_empty());
    }
    #[test]
    fn test_rdf_two_atoms() {
        let frame = make_two_atom_frame(1.0);
        let (r, g) = RdfCalculator::compute(&[frame], 5.0, 50, [10.0; 3]);
        assert_eq!(r.len(), 50);
        assert!(g[0] < 1e-6);
        let bin_1 = (1.0 / (5.0 / 50.0)) as usize;
        assert!(
            g[bin_1] > 0.0,
            "g(r) at r~1 should be > 0, got {}",
            g[bin_1]
        );
    }
    #[test]
    fn test_to_flat_xyz() {
        let frame = make_two_atom_frame(3.0);
        let flat = TrajectoryConverter::to_flat_xyz(&frame);
        assert_eq!(flat.len(), 6);
        assert_eq!(flat[0], 0.0);
        assert_eq!(flat[3], 3.0);
    }
    #[test]
    fn test_from_flat_xyz_roundtrip() {
        let frame = make_two_atom_frame(5.0);
        let flat = TrajectoryConverter::to_flat_xyz(&frame);
        let types = frame.atom_types.clone();
        let recovered = TrajectoryConverter::from_flat_xyz(&flat, types, 0, 0.0).unwrap();
        assert_eq!(recovered.n_atoms(), 2);
        assert!((recovered.positions[1][0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_from_flat_xyz_error_non_multiple_of_3() {
        let flat = vec![1.0, 2.0];
        let result = TrajectoryConverter::from_flat_xyz(&flat, vec!["C".to_string()], 0, 0.0);
        assert!(result.is_err());
    }
    #[test]
    fn test_frame_displacement() {
        let fa = make_two_atom_frame(0.0);
        let fb = make_two_atom_frame(3.0);
        let disp = TrajectoryConverter::frame_displacement(&fa, &fb);
        assert_eq!(disp.len(), 2);
        assert!((disp[0][0]).abs() < 1e-10);
        assert!((disp[1][0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_translate_frame() {
        let frame = make_two_atom_frame(1.0);
        let shifted = TrajectoryConverter::translate(&frame, [10.0, 0.0, 0.0]);
        assert!((shifted.positions[0][0] - 10.0).abs() < 1e-10);
        assert!((shifted.positions[1][0] - 11.0).abs() < 1e-10);
    }
    #[test]
    fn test_scale_frame() {
        let frame = make_two_atom_frame(2.0);
        let scaled = TrajectoryConverter::scale(&frame, 3.0);
        assert!((scaled.positions[1][0] - 6.0).abs() < 1e-10);
    }
}
#[cfg(test)]
mod tests_trajectory_analysis {

    use crate::trajectory::types::*;

    fn atom3_frame(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> TrajectoryFrame {
        TrajectoryFrame::new(
            0,
            0.0,
            vec![a, b, c],
            vec!["C".to_string(), "C".to_string(), "C".to_string()],
        )
    }
    #[test]
    fn test_rmsd_trajectory_identical_is_zero() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let frames = vec![ref_frame.clone()];
        let rmsd = Trajectory::compute_rmsd_trajectory(&frames, &ref_frame);
        assert_eq!(rmsd.len(), 1);
        assert!(
            rmsd[0].abs() < 1e-12,
            "RMSD of identical frames should be 0"
        );
    }
    #[test]
    fn test_rmsd_trajectory_known_displacement() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let shifted = atom3_frame([3.0, 0.0, 0.0], [4.0, 0.0, 0.0], [5.0, 0.0, 0.0]);
        let rmsd = Trajectory::compute_rmsd_trajectory(&[shifted], &ref_frame);
        assert!((rmsd[0] - 3.0).abs() < 1e-10, "RMSD={}", rmsd[0]);
    }
    #[test]
    fn test_rmsd_trajectory_empty_frames() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let rmsd = Trajectory::compute_rmsd_trajectory(&[], &ref_frame);
        assert!(rmsd.is_empty());
    }
    #[test]
    fn test_rmsd_trajectory_atom_count_mismatch_is_nan() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let other = TrajectoryFrame::new(0, 0.0, vec![[0.0, 0.0, 0.0]], vec!["C".to_string()]);
        let rmsd = Trajectory::compute_rmsd_trajectory(&[other], &ref_frame);
        assert!(rmsd[0].is_nan(), "Mismatched counts should give NaN");
    }
    #[test]
    fn test_rg_symmetric_triangle() {
        let sqrt3 = 3.0_f64.sqrt();
        let frame = atom3_frame([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, sqrt3, 0.0]);
        let rg = Trajectory::compute_radius_of_gyration(&[frame]);
        let expected = 2.0 / sqrt3;
        assert!(
            (rg[0] - expected).abs() < 1e-8,
            "Rg={} expected={}",
            rg[0],
            expected
        );
    }
    #[test]
    fn test_rg_empty_frames() {
        let rg = Trajectory::compute_radius_of_gyration(&[]);
        assert!(rg.is_empty());
    }
    #[test]
    fn test_rg_single_atom_is_zero() {
        let frame = TrajectoryFrame::new(0, 0.0, vec![[1.0, 2.0, 3.0]], vec!["H".to_string()]);
        let rg = Trajectory::compute_radius_of_gyration(&[frame]);
        assert!(rg[0].abs() < 1e-12, "Single atom Rg should be 0");
    }
    #[test]
    fn test_align_preserves_atom_count_and_types() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let moved = atom3_frame([5.0, 5.0, 0.0], [6.0, 5.0, 0.0], [5.0, 6.0, 0.0]);
        let aligned = Trajectory::align_to_reference(&[moved], &ref_frame);
        assert_eq!(aligned.len(), 1);
        assert_eq!(aligned[0].n_atoms(), 3);
        assert_eq!(aligned[0].atom_types, ref_frame.atom_types);
    }
    #[test]
    fn test_align_pure_translation_reduces_rmsd() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let moved = atom3_frame([10.0, 10.0, 0.0], [11.0, 10.0, 0.0], [10.0, 11.0, 0.0]);
        let aligned = Trajectory::align_to_reference(std::slice::from_ref(&moved), &ref_frame);
        let rmsd_before = Trajectory::compute_rmsd_trajectory(&[moved], &ref_frame);
        let rmsd_after = Trajectory::compute_rmsd_trajectory(&aligned, &ref_frame);
        assert!(
            rmsd_after[0] < rmsd_before[0],
            "Alignment should reduce RMSD: before={} after={}",
            rmsd_before[0],
            rmsd_after[0]
        );
    }
    #[test]
    fn test_align_mismatch_returns_unchanged() {
        let ref_frame = atom3_frame([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let single = TrajectoryFrame::new(0, 0.0, vec![[99.0, 0.0, 0.0]], vec!["X".to_string()]);
        let aligned = Trajectory::align_to_reference(std::slice::from_ref(&single), &ref_frame);
        assert!((aligned[0].positions[0][0] - 99.0).abs() < 1e-10);
    }
}
