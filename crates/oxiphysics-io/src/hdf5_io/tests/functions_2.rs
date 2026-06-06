//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod new_io_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_read_velocities() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        write_velocities(&mut f, "frame0", &v).unwrap();
        let b = read_velocities(&f, "frame0").unwrap();
        assert_eq!(b, v);
    }
    #[test]
    fn test_write_read_masses() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![1.008, 15.999];
        write_masses(&mut f, "atoms", &m).unwrap();
        assert_eq!(read_masses(&f, "atoms").unwrap(), m);
    }
    #[test]
    fn test_write_read_box_vectors() {
        let mut f = Hdf5File::new("t.h5");
        let bv: [f64; 9] = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0];
        write_box_vectors(&mut f, "cell", &bv).unwrap();
        assert_eq!(read_box_vectors(&f, "cell").unwrap(), bv);
    }
    #[test]
    fn test_write_read_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.4, -0.8, 0.4];
        write_charges(&mut f, "atoms", &q).unwrap();
        assert_eq!(read_charges(&f, "atoms").unwrap(), q);
    }
    #[test]
    fn test_write_read_rdf() {
        let mut f = Hdf5File::new("t.h5");
        let r: Vec<f64> = (0..5).map(|i| i as f64 * 0.1).collect();
        let gr: Vec<f64> = r.iter().map(|&x| x * x).collect();
        write_rdf(&mut f, "obs", &r, &gr).unwrap();
        let (r2, gr2) = read_rdf(&f, "obs").unwrap();
        assert_eq!(r2, r);
        assert_eq!(gr2, gr);
    }
    #[test]
    fn test_write_read_msd() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0];
        let m = vec![0.0, 0.5, 2.0];
        write_msd(&mut f, "msd", &t, &m).unwrap();
        let (t2, m2) = read_msd(&f, "msd").unwrap();
        assert_eq!(t2, t);
        assert_eq!(m2, m);
    }
    #[test]
    fn test_neighbour_list_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let mut nl = NeighbourList::new(3);
        nl.row_ptr = vec![0, 2, 3, 3];
        nl.col_idx = vec![1, 2, 0];
        nl.distances = vec![2.5, 3.1, 2.5];
        write_neighbour_list(&mut f, "nl", &nl).unwrap();
        let b = read_neighbour_list(&f, "nl").unwrap();
        assert_eq!(b.col_idx, nl.col_idx);
        assert_eq!(b.distances, nl.distances);
    }
    #[test]
    fn test_neighbour_list_lookup() {
        let mut nl = NeighbourList::new(2);
        nl.row_ptr = vec![0, 1, 2];
        nl.col_idx = vec![1, 0];
        nl.distances = vec![3.0, 3.0];
        assert_eq!(nl.neighbours(0), &[1usize]);
        assert_eq!(nl.neighbour_distances(1), &[3.0_f64]);
    }
    #[test]
    fn test_sparse_coo_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let mut m = SparseCoo::new(3, 3);
        m.push(0, 0, 1.0);
        m.push(1, 2, 5.0);
        write_sparse_coo(&mut f, "sp", &m).unwrap();
        let b = read_sparse_coo(&f, "sp").unwrap();
        assert_eq!(b.nnz(), 2);
        assert_eq!(b.data[1], 5.0);
    }
    #[test]
    fn test_md_run_metadata_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let meta = MdRunMetadata {
            title: "Water Box".into(),
            n_atoms: 3000,
            n_steps: 1_000_000,
            dt_fs: 2.0,
            temperature_k: 300.0,
            pressure_bar: 1.0,
            ensemble: "NPT".into(),
            force_field: "TIP4P".into(),
        };
        meta.write_to(&mut f).unwrap();
        let b = MdRunMetadata::read_from(&f).unwrap();
        assert_eq!(b.title, "Water Box");
        assert_eq!(b.n_atoms, 3000);
        assert!((b.temperature_k - 300.0).abs() < 1e-12);
    }
    #[test]
    fn test_append_trajectory_frames() {
        let mut f = Hdf5File::new("t.h5");
        let p0 = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let p1 = vec![[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]];
        append_trajectory_frame(&mut f, "traj", 0, &p0).unwrap();
        append_trajectory_frame(&mut f, "traj", 1, &p1).unwrap();
        assert_eq!(count_trajectory_frames(&f, "traj"), 2);
        assert_eq!(read_trajectory_frame(&f, "traj", 0).unwrap(), p0);
    }
    #[test]
    fn test_write_pes_samples() {
        let mut f = Hdf5File::new("t.h5");
        let c: Vec<f64> = (0..6).map(|i| i as f64).collect();
        write_pes_samples(&mut f, "pes", 1, 2, &c, &[-10.5]).unwrap();
        let e = f
            .open_dataset("pes", "energies")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(e, vec![-10.5]);
    }
    #[test]
    fn test_write_aimd_step() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0]];
        let frc = vec![[0.1, -0.2, 0.3]];
        write_aimd_step(
            &mut f,
            "aimd",
            0,
            &pos,
            &frc,
            -42.0,
            &[0.1, 0.2, 0.3, 0.0, 0.0, 0.0],
        )
        .unwrap();
        let g = f.open_group("aimd").unwrap();
        assert!(g.groups.contains_key("step_00000000"));
    }
    #[test]
    fn test_write_nnp_training_batch() {
        let mut f = Hdf5File::new("t.h5");
        write_nnp_training_batch(
            &mut f,
            "nnp",
            2,
            2,
            &[0.1, 0.2, 0.3, 0.4],
            &[-5.0],
            &[0.1, -0.1, 0.0, -0.1, 0.1, 0.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("nnp", "descriptors")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }
    #[test]
    fn test_write_histogram() {
        let mut f = Hdf5File::new("t.h5");
        write_histogram(
            &mut f,
            "h",
            "angle",
            &[0.0, 1.0, 2.0, 3.0],
            &[5.0, 10.0, 3.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("h", "angle_counts")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![5.0, 10.0, 3.0]
        );
    }
    #[test]
    fn test_write_parameter_sweep() {
        let mut f = Hdf5File::new("t.h5");
        write_parameter_sweep(
            &mut f,
            "sw",
            "temperature",
            &[200.0, 300.0],
            "diffusivity",
            &[1e-9, 2e-9],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("sw", "temperature")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![200.0, 300.0]
        );
    }
    #[test]
    fn test_write_remd_metadata() {
        let mut f = Hdf5File::new("t.h5");
        write_remd_metadata(&mut f, "remd", &[300.0, 350.0], &[0.2, 0.3]).unwrap();
        assert_eq!(
            f.open_dataset("remd", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![300.0, 350.0]
        );
    }
    #[test]
    fn test_write_free_energy_profile() {
        let mut f = Hdf5File::new("t.h5");
        write_free_energy_profile(&mut f, "pmf", &[0.0, 0.5, 1.0], &[0.0, 2.0, 0.5]).unwrap();
        assert_eq!(
            f.open_dataset("pmf", "pmf").unwrap().read_f64().unwrap(),
            vec![0.0, 2.0, 0.5]
        );
    }
    #[test]
    fn test_write_mc_run() {
        let mut f = Hdf5File::new("t.h5");
        write_mc_run(
            &mut f,
            "mc",
            &[0, 1, 2],
            &[-10.0, -10.5, -10.5],
            &[true, true, false],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("mc", "accepted")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![1.0, 1.0, 0.0]
        );
    }
    #[test]
    fn test_write_structure_factor() {
        let mut f = Hdf5File::new("t.h5");
        write_structure_factor(&mut f, "sq", &[0.5, 1.0], &[1.0, 2.5]).unwrap();
        assert_eq!(
            f.open_dataset("sq", "sq").unwrap().read_f64().unwrap(),
            vec![1.0, 2.5]
        );
    }
    #[test]
    fn test_write_dielectric_tensor() {
        let mut f = Hdf5File::new("t.h5");
        let eps: [f64; 9] = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
        write_dielectric_tensor(&mut f, "dft", &eps).unwrap();
        assert!(
            (f.open_dataset("dft", "dielectric_tensor")
                .unwrap()
                .read_f64()
                .unwrap()[0]
                - 2.0)
                .abs()
                < 1e-12
        );
    }
    #[test]
    fn test_write_born_charges() {
        let mut f = Hdf5File::new("t.h5");
        let born: Vec<f64> = (0..18).map(|i| i as f64).collect();
        write_born_charges(&mut f, "born", 2, &born).unwrap();
        assert_eq!(
            f.open_dataset("born", "born_charges")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            18
        );
    }
    #[test]
    fn test_write_covariance_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let cov = vec![1.0, 0.5, 0.5, 1.0];
        write_covariance_matrix(&mut f, "cov", 2, &cov).unwrap();
        assert_eq!(
            f.open_dataset("cov", "covariance")
                .unwrap()
                .read_f64()
                .unwrap(),
            cov
        );
    }
    #[test]
    fn test_write_format_version() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 1, 2, 3, "OxiPhysics").unwrap();
        let (maj, min, pat, c) = read_format_version(&f).unwrap();
        assert_eq!((maj, min, pat), (1, 2, 3));
        assert_eq!(c, "OxiPhysics");
    }
    #[test]
    fn test_write_elastic_constants() {
        let mut f = Hdf5File::new("t.h5");
        let mut cij = [0.0_f64; 36];
        cij[0] = 200.0;
        write_elastic_constants(&mut f, "el", &cij).unwrap();
        let b = read_elastic_constants(&f, "el").unwrap();
        assert!((b[0] - 200.0).abs() < 1e-12);
    }
    #[test]
    fn test_write_molecular_orbitals() {
        let mut f = Hdf5File::new("t.h5");
        write_molecular_orbitals(&mut f, "mo", &[-10.0, -5.0], &[2.0, 0.0]).unwrap();
        assert_eq!(
            f.open_dataset("mo", "mo_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-10.0, -5.0]
        );
    }
    #[test]
    fn test_write_homo_lumo_gap() {
        let mut f = Hdf5File::new("t.h5");
        write_homo_lumo_gap(&mut f, "elec", 1.4).unwrap();
        let g = f.open_group("elec").unwrap();
        assert_eq!(g.attributes["homo_lumo_gap_eV"], AttrValue::Float64(1.4));
    }
    #[test]
    fn test_write_pdos() {
        let mut f = Hdf5File::new("t.h5");
        write_pdos(&mut f, "dos", "Fe", &[-1.0, 0.0, 1.0], &[0.2, 0.6, 0.2]).unwrap();
        assert!(f.open_group("dos").unwrap().groups.contains_key("pdos_Fe"));
    }
    #[test]
    fn test_write_gcmc_run() {
        let mut f = Hdf5File::new("t.h5");
        write_gcmc_run(&mut f, "gcmc", &[0, 100], &[10, 11], &[-50.0, -52.0], -4.5).unwrap();
        let g = f.open_group("gcmc").unwrap();
        assert_eq!(g.attributes["chemical_potential"], AttrValue::Float64(-4.5));
    }
    #[test]
    fn test_write_thermal_conductivity() {
        let mut f = Hdf5File::new("t.h5");
        write_thermal_conductivity(
            &mut f,
            "kappa",
            &[0.0, 0.1],
            &[1.0, 0.8],
            1.5,
            &[1.4, 1.5, 1.6],
        )
        .unwrap();
        let g = f.open_group("kappa").unwrap();
        assert_eq!(g.attributes["kappa"], AttrValue::Float64(1.5));
    }
    #[test]
    fn test_total_dataset_count() {
        let mut f = Hdf5File::new("t.h5");
        write_masses(&mut f, "atoms", &[1.0, 2.0]).unwrap();
        write_charges(&mut f, "atoms", &[0.5, -0.5]).unwrap();
        assert!(total_dataset_count(&f) >= 2);
    }
    #[test]
    fn test_list_top_level_groups_sorted() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 0, 1, 0, "t").unwrap();
        write_masses(&mut f, "atoms", &[1.0]).unwrap();
        let groups = list_top_level_groups(&f);
        assert!(groups.windows(2).all(|w| w[0] <= w[1]));
    }
    #[test]
    fn test_write_surface_energies() {
        let mut f = Hdf5File::new("t.h5");
        write_surface_energies(&mut f, "surf", &["100", "111"], &[1.5, 1.2]).unwrap();
        let g = f.open_group("surf").unwrap();
        assert_eq!(g.attributes["111"], AttrValue::Float64(1.2));
    }
    #[test]
    fn test_write_velocity_field_2d() {
        let mut f = Hdf5File::new("t.h5");
        let u = vec![1.0_f64; 4];
        let v = u.clone();
        write_velocity_field_2d(&mut f, "flow", 2, 2, &u, &v).unwrap();
        assert_eq!(f.open_dataset("flow", "u").unwrap().read_f64().unwrap(), u);
    }
    #[test]
    fn test_write_fem_mesh() {
        let mut f = Hdf5File::new("t.h5");
        let nodes = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let elements = vec![0usize, 1, 2];
        write_fem_mesh(&mut f, "mesh", 3, &nodes, 1, 3, &elements).unwrap();
        let g = f.open_group("mesh").unwrap();
        assert_eq!(g.attributes["n_nodes"], AttrValue::Int32(3));
    }
    #[test]
    fn test_write_lbm_populations() {
        let mut f = Hdf5File::new("t.h5");
        let pop = vec![1.0 / 9.0_f64; 36];
        write_lbm_populations(&mut f, "lbm", 2, 2, &pop).unwrap();
        assert_eq!(
            f.open_dataset("lbm", "populations")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            36
        );
    }
    #[test]
    fn test_write_mlp_weights_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let w0 = vec![0.1, 0.2, 0.3, 0.4];
        let b0 = vec![0.0, 0.0];
        let w1 = vec![1.0, -1.0];
        let b1 = vec![0.5];
        write_mlp_weights(&mut f, "net", &[(&w0, &b0, 2, 2), (&w1, &b1, 1, 2)]).unwrap();
        let (rw, rb) = read_mlp_layer(&f, "net", 1).unwrap();
        assert_eq!(rw, w1);
        assert_eq!(rb, b1);
    }
    #[test]
    fn test_write_hparam_trial() {
        let mut f = Hdf5File::new("t.h5");
        write_hparam_trial(&mut f, "search", 0, &[("lr", 1e-3), ("dropout", 0.1)], 0.95).unwrap();
        let g = f.open_group("search").unwrap();
        let sg = g.groups.get("trial_000000").unwrap();
        assert_eq!(sg.attributes["val_metric"], AttrValue::Float64(0.95));
    }
    #[test]
    fn test_write_band_structure() {
        let mut f = Hdf5File::new("t.h5");
        let k: Vec<f64> = (0..6).map(|i| i as f64 * 0.1).collect();
        let e: Vec<f64> = (0..6).map(|i| i as f64 - 3.0).collect();
        write_band_structure(&mut f, "band", 2, 3, &k, &e).unwrap();
        assert_eq!(
            f.open_dataset("band", "eigenvalues")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_sparse_coo_empty() {
        let m = SparseCoo::new(5, 5);
        assert_eq!(m.nnz(), 0);
    }
    #[test]
    fn test_write_polymer_ete() {
        let mut f = Hdf5File::new("t.h5");
        let ete = vec![[10.0, 0.0, 0.0]];
        write_polymer_ete(&mut f, "poly", &ete).unwrap();
        assert_eq!(
            f.open_dataset("poly", "end_to_end")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            3
        );
    }
    #[test]
    fn test_write_ir_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let wn: Vec<f64> = (400..404).map(|i| i as f64).collect();
        let ab = vec![0.1; 4];
        write_ir_spectrum(&mut f, "ir", &wn, &ab).unwrap();
        assert_eq!(
            f.open_dataset("ir", "wavenumbers")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }
    #[test]
    fn test_write_point_cloud_labels() {
        let mut f = Hdf5File::new("t.h5");
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        write_point_cloud(&mut f, "cloud", &pts, Some(&[0i64, 1])).unwrap();
        assert_eq!(
            f.open_dataset("cloud", "labels")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![0.0, 1.0]
        );
    }
    #[test]
    fn test_write_arrhenius_data() {
        let mut f = Hdf5File::new("t.h5");
        write_arrhenius_data(&mut f, "kin", &[500.0, 600.0], &[1e3, 5e3], 0.8).unwrap();
        let g = f.open_group("kin").unwrap();
        assert_eq!(g.attributes["Ea_eV"], AttrValue::Float64(0.8));
    }
    #[test]
    fn test_write_feature_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let feat: Vec<f64> = (0..6).map(|i| i as f64).collect();
        write_feature_matrix(&mut f, "ds", 3, 2, &feat, Some(&[0.0, 1.0, 0.0])).unwrap();
        assert_eq!(
            f.open_dataset("ds", "labels").unwrap().read_f64().unwrap(),
            vec![0.0, 1.0, 0.0]
        );
    }
    #[test]
    fn test_trajectory_frame_count_empty() {
        let f = Hdf5File::new("t.h5");
        assert_eq!(count_trajectory_frames(&f, "no_group"), 0);
    }
    #[test]
    fn test_write_virial_coefficients() {
        let mut f = Hdf5File::new("t.h5");
        write_virial_coefficients(&mut f, "virial", &[-15.0, 200.0]).unwrap();
        assert_eq!(
            f.open_dataset("virial", "virial_coefficients")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-15.0, 200.0]
        );
    }
    #[test]
    fn test_write_order_parameter() {
        let mut f = Hdf5File::new("t.h5");
        let q6 = vec![0.1, 0.2, 0.15];
        write_order_parameter(&mut f, "ops", "q6", &q6).unwrap();
        assert_eq!(f.open_dataset("ops", "q6").unwrap().read_f64().unwrap(), q6);
    }
    #[test]
    fn test_write_stress_tensor() {
        let mut f = Hdf5File::new("t.h5");
        let v = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3];
        write_stress_tensor(&mut f, "sl", 42, &v).unwrap();
        assert!(
            f.open_group("sl")
                .unwrap()
                .datasets
                .contains_key("stress_00000042")
        );
    }
    #[test]
    fn test_write_temperature_series() {
        let mut f = Hdf5File::new("t.h5");
        let t: Vec<f64> = (0..5).map(|i| 300.0 + i as f64).collect();
        write_temperature_series(&mut f, "thermo", &t).unwrap();
        assert_eq!(read_scalar_series(&f, "thermo", "temperature").unwrap(), t);
    }
    #[test]
    fn test_write_pressure_series() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![1.0, 1.01, 0.99];
        write_pressure_series(&mut f, "thermo", &p).unwrap();
        assert_eq!(read_scalar_series(&f, "thermo", "pressure").unwrap(), p);
    }
}
#[cfg(test)]
mod quantum_chem_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_phonon_band_structure() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![0.0; 6];
        let freq = vec![0.0, 5.0, 10.0, 15.0, 20.0, 25.0];
        write_phonon_band_structure(&mut f, "phonon", 2, 3, &q, &freq).unwrap();
        assert_eq!(
            f.open_dataset("phonon", "frequencies_meV")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_write_phonon_group_velocities() {
        let mut f = Hdf5File::new("t.h5");
        let vg = vec![1.0; 18];
        write_phonon_group_velocities(&mut f, "phonon", 2, 3, &vg).unwrap();
        assert_eq!(
            f.open_dataset("phonon", "group_velocities")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            18
        );
    }
    #[test]
    fn test_write_lattice_parameters() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![100.0, 200.0, 300.0];
        let a = vec![5.0, 5.01, 5.02];
        let b = a.clone();
        let c = a.clone();
        write_lattice_parameters(&mut f, "thermal", &t, &a, &b, &c).unwrap();
        assert_eq!(
            f.open_dataset("thermal", "temperatures")
                .unwrap()
                .read_f64()
                .unwrap(),
            t
        );
    }
    #[test]
    fn test_write_euler_angles() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        write_euler_angles(&mut f, "orient", &e).unwrap();
        assert_eq!(
            f.open_dataset("orient", "euler_angles")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_write_electron_density() {
        let mut f = Hdf5File::new("t.h5");
        let rho = vec![0.1_f64; 8];
        write_electron_density(&mut f, "dft", 2, 2, 2, &rho).unwrap();
        assert_eq!(
            f.open_dataset("dft", "rho_e")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            8
        );
    }
    #[test]
    fn test_write_magnetic_moments() {
        let mut f = Hdf5File::new("t.h5");
        let m = vec![[0.0, 0.0, 2.0], [0.0, 0.0, -2.0]];
        write_magnetic_moments(&mut f, "mag", &m).unwrap();
        assert_eq!(
            f.open_dataset("mag", "magnetic_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_write_tight_binding_hamiltonian() {
        let mut f = Hdf5File::new("t.h5");
        let h = vec![0.0; 4];
        write_tight_binding_hamiltonian(&mut f, "tb", 2, &h, &h).unwrap();
        assert_eq!(
            f.open_dataset("tb", "H_real")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            4
        );
    }
    #[test]
    fn test_write_neb_path() {
        let mut f = Hdf5File::new("t.h5");
        let pos0 = vec![[0.0, 0.0, 0.0]];
        let pos1 = vec![[0.5, 0.0, 0.0]];
        let pos2 = vec![[1.0, 0.0, 0.0]];
        let images = vec![pos0, pos1, pos2];
        let e = vec![-10.0, -8.0, -10.0];
        write_neb_path(&mut f, "neb", &images, &e).unwrap();
        let e_back = read_neb_energies(&f, "neb").unwrap();
        assert_eq!(e_back, e);
    }
    #[test]
    fn test_write_transition_state() {
        let mut f = Hdf5File::new("t.h5");
        write_transition_state(&mut f, "ts", -10.0, -5.0, -9.0).unwrap();
        let g = f.open_group("ts").unwrap();
        assert_eq!(g.attributes["barrier"], AttrValue::Float64(5.0));
    }
    #[test]
    fn test_write_gw_corrections() {
        let mut f = Hdf5File::new("t.h5");
        let dft = vec![-3.0, -1.0, 2.0];
        let gw = vec![-3.5, -1.2, 2.3];
        write_gw_corrections(&mut f, "gw", &dft, &gw).unwrap();
        assert_eq!(
            f.open_dataset("gw", "gw_energies")
                .unwrap()
                .read_f64()
                .unwrap(),
            gw
        );
    }
    #[test]
    fn test_write_optical_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 3.0];
        let eps2 = vec![0.1, 0.5, 0.2];
        write_optical_spectrum(&mut f, "opt", &e, &eps2).unwrap();
        assert_eq!(
            f.open_dataset("opt", "epsilon2")
                .unwrap()
                .read_f64()
                .unwrap(),
            eps2
        );
    }
    #[test]
    fn test_write_thermochemistry() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![300.0, 400.0];
        let h = vec![-200.0, -195.0];
        let s = vec![0.1, 0.15];
        let g = vec![-230.0, -255.0];
        write_thermochemistry(&mut f, "thermo", &t, &h, &s, &g).unwrap();
        assert_eq!(
            f.open_dataset("thermo", "enthalpy")
                .unwrap()
                .read_f64()
                .unwrap(),
            h
        );
    }
    #[test]
    fn test_write_force_constants() {
        let mut f = Hdf5File::new("t.h5");
        let ifc = vec![0.1_f64; 36];
        write_force_constants(&mut f, "fc", 2, &ifc).unwrap();
        assert_eq!(
            f.open_dataset("fc", "force_constants")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            36
        );
    }
    #[test]
    fn test_write_dipole_trajectory() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        write_dipole_trajectory(&mut f, "dipole", &d).unwrap();
        assert_eq!(
            f.open_dataset("dipole", "dipole_moments")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            6
        );
    }
    #[test]
    fn test_write_polarisability() {
        let mut f = Hdf5File::new("t.h5");
        let alpha = vec![1.0_f64; 9];
        write_polarisability(&mut f, "pol", &alpha).unwrap();
        assert_eq!(
            f.open_dataset("pol", "polarisability")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }
    #[test]
    fn test_write_bader_charges() {
        let mut f = Hdf5File::new("t.h5");
        let q = vec![7.5, 0.5, 0.5];
        write_bader_charges(&mut f, "bader", &q).unwrap();
        assert_eq!(
            f.open_dataset("bader", "bader_charges")
                .unwrap()
                .read_f64()
                .unwrap(),
            q
        );
    }
    #[test]
    fn test_write_adsorption_energies() {
        let mut f = Hdf5File::new("t.h5");
        write_adsorption_energies(&mut f, "ads", &["fcc", "hcp", "top"], &[-1.5, -1.3, -0.8])
            .unwrap();
        let g = f.open_group("ads").unwrap();
        assert_eq!(g.attributes["fcc"], AttrValue::Float64(-1.5));
    }
    #[test]
    fn test_write_coverage_profile() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 1.0, 2.0];
        let cov = vec![0.0, 0.5, 0.9];
        write_coverage_profile(&mut f, "kinetics", &t, &cov, "CO").unwrap();
        assert_eq!(
            f.open_dataset("kinetics/CO", "coverage")
                .unwrap()
                .read_f64()
                .unwrap(),
            cov
        );
    }
    #[test]
    fn test_write_wannier_centres() {
        let mut f = Hdf5File::new("t.h5");
        let c = vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]];
        let sp = vec![1.2, 1.3];
        write_wannier_centres(&mut f, "wannier", &c, &sp).unwrap();
        assert_eq!(
            f.open_dataset("wannier", "wannier_spreads")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }
    #[test]
    fn test_write_fermi_surface() {
        let mut f = Hdf5File::new("t.h5");
        let k = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let w = vec![1.0, 2.0];
        write_fermi_surface(&mut f, "fs", &k, &w).unwrap();
        assert_eq!(
            f.open_dataset("fs", "fs_weights")
                .unwrap()
                .read_f64()
                .unwrap(),
            w
        );
    }
}
#[cfg(test)]
mod misc_final_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_dose_rate_profile() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![0.0, 1.0, 2.0];
        let r = vec![100.0, 50.0, 25.0];
        write_dose_rate_profile(&mut f, "shielding", &d, &r).unwrap();
        assert_eq!(
            f.open_dataset("shielding", "dose_Gy_hr")
                .unwrap()
                .read_f64()
                .unwrap(),
            r
        );
    }
    #[test]
    fn test_write_stopping_power() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![1.0, 2.0, 5.0, 10.0];
        let sp = vec![200.0, 150.0, 80.0, 50.0];
        write_stopping_power(&mut f, "sp", &e, &sp).unwrap();
        assert_eq!(
            f.open_dataset("sp", "stopping_power")
                .unwrap()
                .read_f64()
                .unwrap(),
            sp
        );
    }
    #[test]
    fn test_write_vacf() {
        let mut f = Hdf5File::new("t.h5");
        let t = vec![0.0, 0.1, 0.2];
        let v = vec![1.0, 0.6, 0.2];
        write_vacf(&mut f, "vd", &t, &v).unwrap();
        assert_eq!(f.open_dataset("vd", "vacf").unwrap().read_f64().unwrap(), v);
    }
    #[test]
    fn test_write_power_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let fr: Vec<f64> = (0..4).map(|i| i as f64).collect();
        let d: Vec<f64> = fr.iter().map(|&x| x * 0.1).collect();
        write_power_spectrum(&mut f, "dos2", &fr, &d).unwrap();
        assert_eq!(
            f.open_dataset("dos2", "frequencies")
                .unwrap()
                .read_f64()
                .unwrap(),
            fr
        );
    }
    #[test]
    fn test_write_bonds_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        let bonds = vec![(0usize, 1), (1, 2)];
        write_bonds(&mut f, "topo", &bonds).unwrap();
        assert_eq!(read_bonds(&f, "topo").unwrap(), bonds);
    }
    #[test]
    fn test_write_charge_density() {
        let mut f = Hdf5File::new("t.h5");
        let data: Vec<f64> = (0..8).map(|i| i as f64).collect();
        write_charge_density(&mut f, "grid", [2, 2, 2], &data).unwrap();
        let (dims, back) = read_charge_density(&f, "grid").unwrap();
        assert_eq!(dims, vec![2, 2, 2]);
        assert_eq!(back, data);
    }
    #[test]
    fn test_write_sph_particles() {
        let mut f = Hdf5File::new("t.h5");
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        write_sph_particles(&mut f, "sph", &pos, &[0.5, 0.5], &[1000.0, 1000.0]).unwrap();
        assert_eq!(
            f.open_dataset("sph", "h").unwrap().read_f64().unwrap(),
            vec![0.5, 0.5]
        );
    }
    #[test]
    fn test_write_nodal_displacements() {
        let mut f = Hdf5File::new("t.h5");
        let d = vec![0.0, 0.01, 0.0, 0.0, 0.02, 0.0];
        write_nodal_displacements(&mut f, "fem", 2, &d).unwrap();
        assert_eq!(
            f.open_dataset("fem", "displacements")
                .unwrap()
                .read_f64()
                .unwrap(),
            d
        );
    }
    #[test]
    fn test_write_sensor_recording() {
        let mut f = Hdf5File::new("t.h5");
        let data: Vec<f64> = (0..20).map(|i| i as f64 * 0.01).collect();
        write_sensor_recording(&mut f, "acc", 4, 5, 1000.0, &data).unwrap();
        let g = f.open_group("acc").unwrap();
        assert_eq!(g.attributes["n_channels"], AttrValue::Int32(4));
    }
    #[test]
    fn test_write_image_stack() {
        let mut f = Hdf5File::new("t.h5");
        let px: Vec<f64> = (0..16).map(|i| i as f64).collect();
        write_image_stack(&mut f, "imgs", 1, 4, 4, &px).unwrap();
        assert_eq!(
            f.open_dataset("imgs", "images")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            16
        );
    }
    #[test]
    fn test_write_covariance_matrix_3x3() {
        let mut f = Hdf5File::new("t.h5");
        let cov: Vec<f64> = (0..9).map(|i| i as f64).collect();
        write_covariance_matrix(&mut f, "cov3", 3, &cov).unwrap();
        assert_eq!(
            f.open_dataset("cov3", "covariance")
                .unwrap()
                .read_f64()
                .unwrap()
                .len(),
            9
        );
    }
    #[test]
    fn test_write_gcmc_n_particles() {
        let mut f = Hdf5File::new("t.h5");
        write_gcmc_run(&mut f, "g2", &[0], &[5], &[-10.0], -3.0).unwrap();
        assert_eq!(
            f.open_dataset("g2", "n_particles")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![5.0]
        );
    }
    #[test]
    fn test_write_read_format_version_default() {
        let mut f = Hdf5File::new("t.h5");
        write_format_version(&mut f, 0, 1, 0, "test").unwrap();
        let (maj, min, pat, c) = read_format_version(&f).unwrap();
        assert_eq!((maj, min, pat), (0, 1, 0));
        assert_eq!(c, "test");
    }
    #[test]
    fn test_write_velocity_field_3d() {
        let mut f = Hdf5File::new("t.h5");
        let n = 8;
        let u = vec![1.0_f64; n];
        let v = vec![2.0_f64; n];
        let w = vec![3.0_f64; n];
        write_velocity_field_3d(&mut f, "fl3d", 2, 2, 2, &u, &v, &w).unwrap();
        assert_eq!(f.open_dataset("fl3d", "w").unwrap().read_f64().unwrap(), w);
    }
    #[test]
    fn test_write_vorticity_field() {
        let mut f = Hdf5File::new("t.h5");
        let omega = vec![0.1, -0.2, 0.3, 0.0];
        write_vorticity_field(&mut f, "vort", 2, 2, &omega).unwrap();
        assert_eq!(
            f.open_dataset("vort", "vorticity")
                .unwrap()
                .read_f64()
                .unwrap(),
            omega
        );
    }
    #[test]
    fn test_write_pressure_field() {
        let mut f = Hdf5File::new("t.h5");
        let p = vec![1.0, 1.5, 1.5, 1.0];
        write_pressure_field(&mut f, "ns", 2, 2, &p).unwrap();
        assert_eq!(
            f.open_dataset("ns", "pressure")
                .unwrap()
                .read_f64()
                .unwrap(),
            p
        );
    }
    #[test]
    fn test_write_diffusion_coefficient() {
        let mut f = Hdf5File::new("t.h5");
        write_diffusion_coefficient(&mut f, "tr", "Na", 1.33e-9, "m^2/s").unwrap();
        let g = f.open_group("tr").unwrap();
        let ds = g.datasets.get("D_Na").unwrap();
        if let DataStorage::Float64(ref v) = ds.data {
            assert!((v[0] - 1.33e-9).abs() < 1e-20);
        }
    }
}
#[cfg(test)]
mod comp_chem_tests {
    use crate::hdf5_io::*;
    #[test]
    fn test_write_cluster_expansion_eci() {
        let mut f = Hdf5File::new("t.h5");
        write_cluster_expansion_eci(&mut f, "ce", &[1, 2, 3], &[-0.5, 0.1, 0.05]).unwrap();
        assert_eq!(
            f.open_dataset("ce", "eci_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-0.5, 0.1, 0.05]
        );
    }
    #[test]
    fn test_write_enumerated_structures() {
        let mut f = Hdf5File::new("t.h5");
        write_enumerated_structures(
            &mut f,
            "enum",
            &[0, 1, 2],
            &[0.0, 0.5, 1.0],
            &[-1.0, -1.5, -1.0],
        )
        .unwrap();
        assert_eq!(
            f.open_dataset("enum", "concentrations")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![0.0, 0.5, 1.0]
        );
    }
    #[test]
    fn test_write_eos_fit() {
        let mut f = Hdf5File::new("t.h5");
        let v = vec![60.0, 65.0, 70.0, 75.0];
        let e = vec![-4.0, -4.5, -4.3, -4.0];
        write_eos_fit(&mut f, "eos", &v, &e, 68.0, -4.5, 160.0, 4.0).unwrap();
        let g = f.open_group("eos").unwrap();
        assert_eq!(g.attributes["B0_GPa"], AttrValue::Float64(160.0));
    }
    #[test]
    fn test_write_convex_hull() {
        let mut f = Hdf5File::new("t.h5");
        let x = vec![0.0, 0.5, 1.0];
        let e = vec![0.0, -0.1, 0.0];
        let s = vec![true, true, true];
        write_convex_hull(&mut f, "hull", &x, &e, &s).unwrap();
        let sm = f
            .open_dataset("hull", "stable_mask")
            .unwrap()
            .read_f64()
            .unwrap();
        assert_eq!(sm, vec![1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_write_solvation_energies() {
        let mut f = Hdf5File::new("t.h5");
        let solvents = vec!["water".to_string(), "methanol".to_string()];
        let dg = vec![-5.0, -3.5];
        write_solvation_energies(&mut f, "solv", &solvents, &dg).unwrap();
        assert_eq!(
            f.open_dataset("solv", "dG_solv_kcal_mol")
                .unwrap()
                .read_f64()
                .unwrap(),
            dg
        );
    }
    #[test]
    fn test_write_excitations() {
        let mut f = Hdf5File::new("t.h5");
        let e = vec![3.5, 4.2, 5.1];
        let osc = vec![0.5, 0.1, 0.8];
        write_excitations(&mut f, "tddft", &e, &osc).unwrap();
        assert_eq!(
            f.open_dataset("tddft", "oscillator_strengths")
                .unwrap()
                .read_f64()
                .unwrap(),
            osc
        );
    }
    #[test]
    fn test_write_nmr_shielding() {
        let mut f = Hdf5File::new("t.h5");
        let iso = vec![25.0, 120.0];
        let aniso = vec![5.0, 30.0];
        write_nmr_shielding(&mut f, "nmr", &iso, &aniso).unwrap();
        assert_eq!(
            f.open_dataset("nmr", "sigma_iso")
                .unwrap()
                .read_f64()
                .unwrap(),
            iso
        );
    }
    #[test]
    fn test_write_hyperfine_coupling() {
        let mut f = Hdf5File::new("t.h5");
        let a = vec![50.0, 10.0, 5.0];
        write_hyperfine_coupling(&mut f, "hfc", &a).unwrap();
        assert_eq!(
            f.open_dataset("hfc", "A_iso_MHz")
                .unwrap()
                .read_f64()
                .unwrap(),
            a
        );
    }
    #[test]
    fn test_write_nics() {
        let mut f = Hdf5File::new("t.h5");
        let h = vec![0.0, 1.0];
        let n = vec![-10.5, -5.2];
        write_nics(&mut f, "arom", &h, &n).unwrap();
        assert_eq!(
            f.open_dataset("arom", "NICS_values")
                .unwrap()
                .read_f64()
                .unwrap(),
            n
        );
    }
    #[test]
    fn test_write_charge_flow_matrix() {
        let mut f = Hdf5File::new("t.h5");
        let dq = vec![0.0, 0.1, -0.1, 0.0];
        write_charge_flow_matrix(&mut f, "cflow", 2, &dq).unwrap();
        assert_eq!(
            f.open_dataset("cflow", "charge_flow")
                .unwrap()
                .read_f64()
                .unwrap(),
            dq
        );
    }
    #[test]
    fn test_write_fukui_functions() {
        let mut f = Hdf5File::new("t.h5");
        let fp = vec![0.3, 0.7];
        let fm = vec![0.4, 0.6];
        let fz = vec![0.35, 0.65];
        write_fukui_functions(&mut f, "fukui", &fp, &fm, &fz).unwrap();
        assert_eq!(
            f.open_dataset("fukui", "f_plus")
                .unwrap()
                .read_f64()
                .unwrap(),
            fp
        );
    }
    #[test]
    fn test_write_potential_energy_series() {
        let mut f = Hdf5File::new("t.h5");
        let e: Vec<f64> = (0..5).map(|i| -(i as f64) * 10.0).collect();
        write_potential_energy_series(&mut f, "thermo", &e).unwrap();
        assert_eq!(
            read_scalar_series(&f, "thermo", "potential_energy").unwrap(),
            e
        );
    }
    #[test]
    fn test_write_von_mises_stress() {
        let mut f = Hdf5File::new("t.h5");
        let s = vec![150e6, 200e6, 180e6];
        write_von_mises_stress(&mut f, "result", &s).unwrap();
        assert_eq!(
            f.open_dataset("result", "von_mises")
                .unwrap()
                .read_f64()
                .unwrap(),
            s
        );
    }
    #[test]
    fn test_write_radius_of_gyration() {
        let mut f = Hdf5File::new("t.h5");
        let rg = vec![3.5, 3.6, 3.4];
        write_radius_of_gyration(&mut f, "poly", &rg).unwrap();
        assert_eq!(
            f.open_dataset("poly", "rg").unwrap().read_f64().unwrap(),
            rg
        );
    }
    #[test]
    fn test_write_raman_spectrum() {
        let mut f = Hdf5File::new("t.h5");
        let wn = vec![200.0, 400.0, 800.0];
        let i = vec![100.0, 500.0, 200.0];
        write_raman_spectrum(&mut f, "raman", &wn, &i).unwrap();
        assert_eq!(
            f.open_dataset("raman", "intensities")
                .unwrap()
                .read_f64()
                .unwrap(),
            i
        );
    }
    #[test]
    fn test_write_virial_roundtrip() {
        let mut f = Hdf5File::new("t.h5");
        write_virial_coefficients(&mut f, "v2", &[-15.0, 200.0]).unwrap();
        assert_eq!(
            f.open_dataset("v2", "virial_coefficients")
                .unwrap()
                .read_f64()
                .unwrap(),
            vec![-15.0, 200.0]
        );
    }
}
