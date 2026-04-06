//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use crate::netcdf::types::*;

#[cfg(test)]
mod tests_netcdf_ext {
    use super::*;
    #[test]
    fn frame_n_atoms_correct() {
        let f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]; 5],
            velocities: None,
            box_lengths: None,
        };
        assert_eq!(f.n_atoms(), 5);
    }
    #[test]
    fn frame_centre_of_mass_symmetric() {
        let f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: None,
        };
        let com = f.centre_of_mass();
        assert!((com[0]).abs() < 1e-12);
    }
    #[test]
    fn frame_rmsd_self_is_zero() {
        let f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: None,
            box_lengths: None,
        };
        assert!(f.rmsd_from(&f) < 1e-12);
    }
    #[test]
    fn frame_rmsd_translated_exact() {
        let f0 = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]],
            velocities: None,
            box_lengths: None,
        };
        let f1 = TrajectoryFrame {
            time_ps: 1.0,
            positions: vec![[1.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: None,
        };
        assert!((f1.rmsd_from(&f0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn frame_kinetic_energy_no_velocity_is_zero() {
        let f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]],
            velocities: None,
            box_lengths: None,
        };
        assert!((f.kinetic_energy(1.0)).abs() < 1e-12);
    }
    #[test]
    fn frame_kinetic_energy_unit_mass_unit_velocity() {
        let f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]],
            velocities: Some(vec![[1.0, 0.0, 0.0]]),
            box_lengths: None,
        };
        assert!((f.kinetic_energy(1.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn frame_translate_moves_positions() {
        let mut f = TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: None,
        };
        f.translate([1.0, 2.0, 3.0]);
        assert_eq!(f.positions[0], [1.0, 2.0, 3.0]);
    }
    #[test]
    fn builder_frame_count_after_add() {
        let mut b = NetcdfTrajectoryBuilder::new();
        b.add_frame(TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]],
            velocities: None,
            box_lengths: None,
        });
        b.add_frame(TrajectoryFrame {
            time_ps: 1.0,
            positions: vec![[1.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: None,
        });
        assert_eq!(b.frame_count(), 2);
    }
    #[test]
    fn builder_rmsd_series_first_zero() {
        let mut b = NetcdfTrajectoryBuilder::new();
        b.add_frame(TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]],
            velocities: None,
            box_lengths: None,
        });
        b.add_frame(TrajectoryFrame {
            time_ps: 1.0,
            positions: vec![[1.0, 0.0, 0.0]],
            velocities: None,
            box_lengths: None,
        });
        let rmsd = b.rmsd_series();
        assert!((rmsd[0]).abs() < 1e-12);
        assert!((rmsd[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn builder_time_series_correct() {
        let mut b = NetcdfTrajectoryBuilder::new();
        for t in [0.0_f64, 0.5, 1.0] {
            b.add_frame(TrajectoryFrame {
                time_ps: t,
                positions: vec![],
                velocities: None,
                box_lengths: None,
            });
        }
        let ts = b.time_series();
        assert_eq!(ts, vec![0.0, 0.5, 1.0]);
    }
    #[test]
    fn builder_write_cdl_contains_conventions() {
        let b = NetcdfTrajectoryBuilder::new()
            .with_title("Test")
            .in_angstroms();
        let mut buf = Vec::new();
        b.write_cdl(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("AMBER"));
        assert!(s.contains("angstrom"));
        assert!(s.contains("Test"));
    }
    #[test]
    fn builder_n_atoms_from_first_frame() {
        let mut b = NetcdfTrajectoryBuilder::new();
        b.add_frame(TrajectoryFrame {
            time_ps: 0.0,
            positions: vec![[0.0; 3]; 7],
            velocities: None,
            box_lengths: None,
        });
        assert_eq!(b.n_atoms(), 7);
    }
    #[test]
    fn dimspec_fixed_cdl() {
        let d = NetcdfDimSpec::fixed("atom", 100);
        let s = d.to_cdl();
        assert!(s.contains("atom"));
        assert!(s.contains("100"));
        assert!(!s.contains("UNLIMITED"));
    }
    #[test]
    fn dimspec_unlimited_cdl() {
        let d = NetcdfDimSpec::unlimited("frame", 50);
        let s = d.to_cdl();
        assert!(s.contains("UNLIMITED"));
        assert!(s.contains("50"));
    }
    #[test]
    fn write_cdl_dimensions_has_header() {
        let dims = vec![
            NetcdfDimSpec::fixed("x", 10),
            NetcdfDimSpec::unlimited("t", 5),
        ];
        let mut buf = Vec::new();
        write_cdl_dimensions(&mut buf, &dims).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("dimensions:"));
        assert!(s.contains("UNLIMITED"));
    }
    #[test]
    fn angstroms_to_nm_single() {
        let p = angstroms_to_nm(&[[10.0, 0.0, 0.0]]);
        assert!((p[0][0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn nm_to_angstroms_roundtrip() {
        let orig = vec![[1.5_f64, 2.0, 0.5]];
        let rt = angstroms_to_nm(&nm_to_angstroms(&orig));
        assert!((rt[0][0] - orig[0][0]).abs() < 1e-12);
    }
    #[test]
    fn time_vector_length_and_values() {
        let tv = time_vector(1.0, 0.5, 4);
        assert_eq!(tv.len(), 4);
        assert!((tv[3] - 2.5).abs() < 1e-12);
    }
    #[test]
    fn msd_1d_first_zero() {
        let msd = msd_1d(&[1.0, 2.0, 3.0]);
        assert!((msd[0]).abs() < 1e-12);
        assert!((msd[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn slice_mean_basic() {
        let m = slice_mean(&[1.0, 2.0, 3.0, 4.0]);
        assert!((m - 2.5).abs() < 1e-12);
    }
    #[test]
    fn slice_variance_uniform_zero() {
        let v = slice_variance(&[5.0, 5.0, 5.0]);
        assert!(v.abs() < 1e-12);
    }
    #[test]
    fn autocorrelation_lag0_is_one() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ac = autocorrelation(&data, 0);
        assert!((ac - 1.0).abs() < 1e-10);
    }
    #[test]
    fn autocorrelation_large_lag_zero() {
        let data = vec![1.0, 2.0, 3.0];
        let ac = autocorrelation(&data, 10);
        assert_eq!(ac, 0.0);
    }
}
